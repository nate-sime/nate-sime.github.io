/**
 * Coupled Boussinesq convection, GPU-resident, with live controls.
 *
 * The whole pipeline (buoyancy assembly → biharmonic Stokes → u = ∇×ψ →
 * semi-Lagrangian/BFECC advection → implicit diffusion → temperature map, with
 * ψ-isocontour streamlines and the mesh as optional overlays) runs in WebGPU
 * compute and fragment shaders.
 * Nothing crosses back to the host in the frame loop; the diagnostics readout is
 * an asynchronous poll of a GPU-side reduction, deliberately off the frame's
 * dependency chain. Those polls also accumulate into the two corner traces —
 * Nusselt number (`ui/nuplot.ts`) and RMS velocity (`ui/rmsplot.ts`) — which
 * are the same reductions read as time series rather than as instants.
 *
 * Resolution and dt are fixed per run: dt is an accuracy knob, not a
 * stability limit, and sizing it from the advective CFL would need a max-|u|
 * readback every frame. Both are adjustable from the pane — at their true cost,
 * announced rather than hidden.
 */

import { GpuSimulation } from "./gpu/sim";
import { boundaryNames } from "./geometry";
import { gammaFor } from "./solver/rheology";
import {
  buildPane, defaultState, geometryFor, MESH, PRESETS, VISCOSITY, type State,
} from "./ui/controls";
import { dimensionalTime, referenceNote } from "./ui/dimensional";
import { NusseltPlot } from "./ui/nuplot";
import { RmsPlot } from "./ui/rmsplot";

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
  const nu = new NusseltPlot(el("nu"), state.nuWindow, geometryFor(state).kind);
  // Shares the Nu plot's window control (`state.nuWindow`) rather than getting
  // a slider of its own: both are the same frame loop's poll, at the same
  // cadence, and a second "how much of the run" control next to the first
  // would offer the reader two knobs with nothing to distinguish them.
  const rms = new RmsPlot(el("rms"), state.nuWindow, geometryFor(state).kind);
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
    notice(`building ${s.geometry}, ${p} — factorising radial operators, `
      + `compiling pipelines…`);
    await new Promise(requestAnimationFrame);
    sim?.destroy();
    sim = null;
    const { nr, na, gnr, gna } = PRESETS[p];
    const { variable, strainRate } = VISCOSITY[s.viscosity];
    const geom = geometryFor(s);
    const next = GpuSimulation.create(device, format, {
      nr, na, gnr, gna, geom,
      Ra: 10 ** s.logRa, dt: s.dt,
      levels: s.contours, lineW: s.lineWidth, mesh: MESH[s.mesh],
      variable, gamma: gammaFor(10 ** s.logContrast),
      cz: gammaFor(10 ** s.logDepthContrast), iters: s.iters,
      n: strainRate ? s.n : 1, picard: s.picard,
    });
    next.reseed(0.05, s.wavenumber);
    sim = next;
    carry = 0;
    // The trace belonged to the solver just destroyed. `NuTrace.push` would drop
    // it anyway once the clock came back from zero, but that is a floor, not the
    // behaviour to rely on: clearing here empties the panel with the notice
    // rather than a second later, when the first poll of the new run lands.
    // The geometry goes with it: the panel's clock and its two series names are
    // both properties of the domain, not of the samples.
    nu.setGeometry(geom.kind);
    nu.clear();
    rms.setGeometry(geom.kind);
    rms.clear();
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
    onNuWindow: (steps) => { nu.setWindow(steps); rms.setWindow(steps); },
    onReseed: () => {
      sim?.reseed(0.05, state.wavenumber);
      nu.clear();
      rms.clear();
    },
    onResolution: (p) => {
      // dt is resolution-dependent (finer grids want smaller steps), so a
      // preset carries its own; adopt it and reflect that back into the pane.
      state.dt = PRESETS[p].dt;
      void build(state);
    },
    // The metric is compiled into the shaders and the box length reaches the
    // knot vector, so neither half of the domain is anything but a full rebuild.
    onGeometry: () => void build(state),
    // Entering or leaving the Krylov tier changes which buffers exist, so that
    // is a rebuild. Moving between the two variable laws is not: n = 1 collapses
    // the power law exactly, so it is one uniform write — see `VISCOSITY` in
    // presets.ts.
    onViscosity: (v) => {
      const { variable, strainRate } = VISCOSITY[v];
      if (!sim || variable !== sim.o.variable) return void build(state);
      sim.n = strainRate ? state.n : 1;
    },
    onContrast: (log10, log10Depth) => {
      // Captured, not re-read after the await: `announce` yields a frame so the
      // notice paints, and a rebuild in that gap would leave `sim` pointing at a
      // solver whose buffers this call did not intend to touch.
      const s = sim;
      if (s?.o.variable)
        void announce("re-inverting the μ̄(r) radial blocks…", () => {
          if (sim === s)
            s.setContrast(gammaFor(10 ** log10), gammaFor(10 ** log10Depth));
        });
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
      const { nuInner, nuOuter, psiMax, vrms, at, atStep } = sim.stats;
      // Offered every frame rather than on the poll cadence: `stats` is replaced
      // asynchronously, so there is no frame the host can name as the one a
      // reading arrived on. The trace drops the NaNs before the first readback
      // and the repeats of a sample already held, so this is the same series
      // either way, at most one frame sooner.
      nu.push({ t: at, step: atStep, inner: nuInner, outer: nuOuter });
      rms.push({ t: at, step: atStep, v: vrms });
      const power = sim.o.variable && sim.n !== 1;
      // Named only when it is doing something: at c = 0 the depth term is
      // exactly absent, and a "d" in the law with "× 10^0.00" after it would
      // describe a dependence the solver does not have.
      const deep = sim.o.variable && sim.cz !== 0;
      const law = sim.o.variable
        ? `μ(T${deep ? ", d" : ""}${power ? ", ε̇" : ""}), ` +
          `contrast 10^${state.logContrast.toFixed(2)}` +
          (deep ? ` × 10^${state.logDepthContrast.toFixed(2)} with depth` : "") +
          (power ? `, n = ${sim.n}` : "")
        : "constant viscosity";
      const g = sim.o.geom;
      const bn = boundaryNames(g.kind);
      // The domain, said as its own dimensions: a radius ratio for the annulus,
      // a width × depth for the box. Both are what the reader would have to
      // measure off the canvas otherwise, and the box's is the one control on
      // the pane whose effect is a number rather than a picture.
      //
      // A walled box names the period it is *solved* on as well as the width it
      // is drawn at. The two differ by the mirror, and a reader comparing the
      // resolution line below against the picture would otherwise be off by two.
      const domain = g.kind === "annulus"
        ? `2-D spherical annulus r ${g.lo} … ${g.hi}`
        : `2-D Cartesian box ${g.width} × ${g.hi - g.lo}   ` +
          (g.walls === "free-slip"
            ? `free-slip walls at x = 0, ${g.width} (mirrored, period ${g.span})`
            : `periodic in x`);
      log.textContent =
        `${domain} · Boussinesq convection · WebGPU\n` +
        `Ra = ${(10 ** state.logRa).toExponential(2)}   ${law}   free-slip\n` +
        `ψ ${sim.nr}×${sim.na} splines   T ${sim.gnr}×${sim.gna} grid   dt = ${sim.dt.toExponential(1)}\n` +
        // The scaling behind every dimensional figure below, stated where the
        // rest of the configuration is. Without it "133 Gyr" is a number and not
        // a quantity — and the reference is a display assumption, not something
        // the solver knows (see ui/dimensional.ts).
        `${referenceNote(g.kind)}\n\n` +
        `step ${String(sim.steps).padStart(6)}   t = ${sim.time.toFixed(4)} = ` +
        `${dimensionalTime(g.kind, sim.time)}   ` +
        `${fps.toFixed(0)} fps   ${rate(state.speed)}${state.paused ? "   paused" : ""}\n` +
        `Nu   ${bn.inner} ${n(nuInner)}   ${bn.outer} ${n(nuOuter)}   v_rms ${n(vrms, 3)}\n` +
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
