/**
 * Coupled Boussinesq convection, GPU-resident, with live controls.
 *
 * The whole pipeline (buoyancy assembly → biharmonic Stokes → u = ∇×ψ →
 * semi-Lagrangian/BFECC advection → implicit diffusion → temperature map, with
 * ψ-isocontour streamlines and the mesh as optional overlays) runs in WebGPU
 * compute and fragment shaders.
 * Nothing crosses back to the host in the frame loop; the diagnostics readout is
 * an asynchronous poll of a GPU-side reduction, deliberately off the frame's
 * dependency chain.
 *
 * Resolution and dt are fixed per run: dt is an accuracy knob, not a
 * stability limit, and sizing it from the advective CFL would need a max-|u|
 * readback every frame. Both are adjustable from the pane — at their true cost,
 * announced rather than hidden.
 */

import { GpuSimulation } from "./gpu/sim";
import { gammaFor } from "./solver/rheology";
import { buildPane, defaultState, MESH, PRESETS, VISCOSITY, type State } from "./ui/controls";

const RI = 0.55, RO = 1.0;

const el = (id: string) => document.getElementById(id)!;

/** Step rate for the readout, in whichever direction reads naturally. */
const rate = (v: number) =>
  v >= 1 ? `${v} step${v > 1 ? "s" : ""}/frame` : `1 step / ${Math.round(1 / v)} frames`;
const notice = (msg: string): void => {
  el("msg").textContent = msg;
  el("msg").setAttribute("data-show", "");
};

async function main(): Promise<void> {
  if (!navigator.gpu) return notice("WebGPU is unavailable in this browser.");
  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) return notice("No suitable GPU adapter was found.");
  const device = await adapter.requestDevice();

  const canvas = el("view") as HTMLCanvasElement;
  const ctx = canvas.getContext("webgpu")!;
  const format = navigator.gpu.getPreferredCanvasFormat();
  ctx.configure({ device, format, alphaMode: "opaque" });

  const resize = (): void => {
    const s = Math.min(devicePixelRatio, 2);
    const side = Math.max(1, (Math.min(canvas.clientWidth, canvas.clientHeight) * s) | 0);
    canvas.width = canvas.height = side;
  };
  new ResizeObserver(resize).observe(canvas);
  resize();

  const state = defaultState();
  const log = el("log");
  let sim: GpuSimulation | null = null;

  // Fractional step rates need an accumulator: at 1/16 the loop steps on one
  // frame in sixteen. `speed` is bounded by 16, so at most that many steps fall
  // due in a frame. Declared before `build`, which clears it — a rebuilt solver
  // must not inherit a partial credit earned by the one it replaced.
  let carry = 0;

  /**
   * Build (or rebuild) the solver. Factorising the radial operators in f64 and
   * compiling every pipeline blocks for a second or two, so the notice is
   * painted first and the old solver's buffers are released before the new
   * ones are claimed.
   */
  const build = async (s: State): Promise<void> => {
    const p = s.resolution;
    notice(`building ${p} — factorising radial operators, compiling pipelines…`);
    await new Promise(requestAnimationFrame);
    sim?.destroy();
    sim = null;
    const { nr, na, gnr, gna } = PRESETS[p];
    const { variable, strainRate } = VISCOSITY[s.viscosity];
    const next = GpuSimulation.create(device, format, {
      nr, na, gnr, gna, ri: RI, ro: RO,
      Ra: 10 ** s.logRa, dt: s.dt,
      levels: s.contours, lineW: s.lineWidth, mesh: MESH[s.mesh],
      variable, gamma: gammaFor(10 ** s.logContrast), iters: s.iters,
      n: strainRate ? s.n : 1, picard: s.picard,
    });
    next.reseed(0.05, s.wavenumber);
    sim = next;
    carry = 0;
    el("msg").removeAttribute("data-show");
  };

  await build(state);

  /** Announce, yield a frame so the notice paints, then run the f64 job. */
  const announce = async (msg: string, job: () => void): Promise<void> => {
    notice(msg);
    await new Promise(requestAnimationFrame);
    job();
    el("msg").removeAttribute("data-show");
  };

  buildPane(state, {
    onRa: (v) => { if (sim) sim.Ra = v; },
    onDt: (v) => sim?.setDt(v),
    onStreamlines: (levels, lineW) => sim?.setStreamlines(levels, lineW),
    onMesh: (m) => { if (sim) sim.mesh = MESH[m]; },
    onReseed: () => sim?.reseed(0.05, state.wavenumber),
    onResolution: (p) => {
      // dt is resolution-dependent (finer grids want smaller steps), so a
      // preset carries its own; adopt it and reflect that back into the pane.
      state.dt = PRESETS[p].dt;
      void build(state);
    },
    // Entering or leaving the Krylov tier changes which buffers exist, so that
    // is a rebuild. Moving between the two variable laws is not: n = 1 collapses
    // the power law exactly, so it is one uniform write — see `VISCOSITY` in
    // presets.ts.
    onViscosity: (v) => {
      const { variable, strainRate } = VISCOSITY[v];
      if (!sim || variable !== sim.o.variable) return void build(state);
      sim.n = strainRate ? state.n : 1;
    },
    onContrast: (log10) => {
      // Captured, not re-read after the await: `announce` yields a frame so the
      // notice paints, and a rebuild in that gap would leave `sim` pointing at a
      // solver whose buffers this call did not intend to touch.
      const s = sim;
      if (s?.o.variable)
        void announce("re-inverting the μ̄(r) radial blocks…",
          () => { if (sim === s) s.setGamma(gammaFor(10 ** log10)); });
    },
    onIters: (n) => { if (sim) sim.iters = n; },
    onExponent: (n) => { if (sim) sim.n = n; },
    onPicard: (n) => { if (sim) sim.picard = n; },
  });

  let frames = 0, fps = 0, last = performance.now();
  const frame = (): void => {
    if (sim) {
      if (state.paused) {
        carry = 0;
      } else {
        carry += state.speed;
        const due = Math.floor(carry);
        carry -= due;
        for (let n = 0; n < due; n++) sim.step();
      }
      sim.draw(ctx.getCurrentTexture().createView());

      if (++frames % 15 === 0) {
        sim.pollStats();
        const now = performance.now();
        fps = (15 * 1000) / (now - last);
        last = now;
      }
      // The reductions are one poll behind at worst, and NaN until the first
      // lands — say so rather than printing a number that is not there yet.
      const n = (v: number, d = 4) => (Number.isNaN(v) ? "—" : v.toFixed(d));
      const { nuInner, nuOuter, psiMax } = sim.stats;
      const power = sim.o.variable && sim.n !== 1;
      const law = sim.o.variable
        ? `μ(T${power ? ", ε̇" : ""}), contrast 10^${state.logContrast.toFixed(2)}` +
          (power ? `, n = ${sim.n}` : "")
        : "constant viscosity";
      log.textContent =
        `2-D spherical annulus · Boussinesq convection · WebGPU\n` +
        `Ra = ${(10 ** state.logRa).toExponential(2)}   ${law}   free-slip\n` +
        `ψ ${sim.nr}×${sim.na} splines   T ${sim.gnr}×${sim.gna} grid   dt = ${sim.dt.toExponential(1)}\n\n` +
        `step ${String(sim.steps).padStart(6)}   t = ${sim.time.toFixed(4)}   ` +
        `${fps.toFixed(0)} fps   ${rate(state.speed)}${state.paused ? "   paused" : ""}\n` +
        `Nu   inner ${n(nuInner)}   outer ${n(nuOuter)}\n` +
        `max |ψ| ${n(psiMax, 3)}` +
        // The budget, not a residual: see `pollStats` on why a residual is not a
        // convergence diagnostic for this operator once ψ is stored in f32.
        (sim.o.variable ? `   ${sim.iters} CG iterations/solve` : "") +
        (power && sim.picard > 1 ? ` × ${sim.picard} Picard sweeps` : "");
    }
    requestAnimationFrame(frame);
  };
  requestAnimationFrame(frame);
}

main();
