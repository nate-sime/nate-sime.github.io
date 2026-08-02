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
 * Resolution is fixed per run — a rebuild, announced. dt is not: it is an
 * accuracy knob, not a stability limit, and is sized every poll from the
 * advective CFL via a GPU-side max-|u| reduction (`cflSource` in `gpu/wgsl.ts`)
 * read back asynchronously and turned into a step by `adaptiveDt`, which
 * hysteresis-gates the f64 refactorisation `setDt` costs. The pane exposes
 * what actually drives it — a target Courant number — plus the cap it is held
 * under.
 */

import { adaptiveDt } from "./adaptiveDt";
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

  // ---- zoom / pan --------------------------------------------------------
  //
  // Tracked on the host, not read back from the GPU: `sim.setView` is a plain
  // uniform write, so the JS side is the only place that knows the current
  // zoom/pan between frames. A rebuild replaces the solver but not the view
  // the reader had chosen — `applyView` below re-sends it to whatever `sim`
  // currently points at, which is also what makes wheel/drag handlers safe to
  // register once, before the first `build`.
  const ZOOM_MIN = 1, ZOOM_MAX = 40;
  let view = { zoom: 1, panX: 0, panY: 0 };
  const applyView = (): void => sim?.setView(view.zoom, view.panX, view.panY);
  const resetView = (): void => {
    view = { zoom: 1, panX: 0, panY: 0 };
    applyView();
  };
  // however far zoomed in, keep at least the reciprocal fraction of the fixed
  // [-halfExtent, halfExtent] window in view — panning past that would show
  // nothing but background on that side no matter which way the domain sits.
  const clampPan = (he: number): void => {
    const max = he * (1 - 1 / view.zoom);
    view.panX = Math.min(max, Math.max(-max, view.panX));
    view.panY = Math.min(max, Math.max(-max, view.panY));
  };
  // Screen → world, the same map the vertex shader applies (see wgsl.ts `vs`):
  // p = ndc * halfExtent/zoom + pan. NDC flips y — WebGPU clip space is y-up,
  // the DOM is y-down.
  const ndcOf = (e: PointerEvent | WheelEvent): [number, number] => {
    const r = canvas.getBoundingClientRect();
    return [
      ((e.clientX - r.left) / r.width) * 2 - 1,
      1 - ((e.clientY - r.top) / r.height) * 2,
    ];
  };

  canvas.addEventListener("wheel", (e) => {
    if (!sim) return;
    e.preventDefault();
    const [nx, ny] = ndcOf(e);
    const he = sim.halfExtent;
    const z0 = view.zoom;
    const z1 = Math.min(ZOOM_MAX, Math.max(ZOOM_MIN, z0 * Math.exp(-e.deltaY * 0.0015)));
    // Shift the pan by what the zoom change alone would have moved the point
    // under the cursor, so that point is what stays under it.
    view.panX += nx * he * (1 / z0 - 1 / z1);
    view.panY += ny * he * (1 / z0 - 1 / z1);
    view.zoom = z1;
    clampPan(he);
    applyView();
    canvas.style.cursor = view.zoom > 1 ? "grab" : "default";
  }, { passive: false });

  let dragging = false, lastX = 0, lastY = 0;
  canvas.addEventListener("pointerdown", (e) => {
    if (!sim || view.zoom <= 1) return;
    dragging = true;
    lastX = e.clientX; lastY = e.clientY;
    canvas.setPointerCapture(e.pointerId);
    canvas.style.cursor = "grabbing";
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!dragging || !sim) return;
    const r = canvas.getBoundingClientRect();
    const he = sim.halfExtent;
    const dx = e.clientX - lastX, dy = e.clientY - lastY;
    lastX = e.clientX; lastY = e.clientY;
    view.panX -= (dx * 2 / r.width) * he / view.zoom;
    view.panY += (dy * 2 / r.height) * he / view.zoom;
    clampPan(he);
    applyView();
  });
  const endDrag = (e: PointerEvent): void => {
    dragging = false;
    canvas.releasePointerCapture(e.pointerId);
    canvas.style.cursor = view.zoom > 1 ? "grab" : "default";
  };
  canvas.addEventListener("pointerup", endDrag);
  canvas.addEventListener("pointercancel", endDrag);
  canvas.addEventListener("dblclick", () => {
    resetView();
    canvas.style.cursor = "default";
  });

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
      // Seeded at the cap: `maxSpeed` is not known until the constructor's own
      // Stokes solve has run, so there is no CFL-implied value yet to start
      // from, and the cap is the safe upper bound for whatever that turns out
      // to be. The first poll shrinks it to the real adaptive step.
      Ra: 10 ** s.logRa, dt: s.dtMax,
      levels: s.contours, lineW: s.lineWidth, mesh: MESH[s.mesh],
      colormap: s.colormap,
      variable, gamma: gammaFor(10 ** s.logContrast),
      cz: gammaFor(10 ** s.logDepthContrast), iters: s.iters,
      n: strainRate ? s.n : 1, picard: s.picard,
      tackley: VISCOSITY[s.viscosity].tackley,
      brandenburg: VISCOSITY[s.viscosity].brandenburg,
      sigmaY: s.sigmaY, sigmaB: s.sigmaB, etaStar: s.etaStar,
      aUpper: s.aUpper, aLower: s.aLower, d0: s.d0,
    });
    next.reseed(0.05, s.wavenumber);
    sim = next;
    carry = 0;
    // A rebuilt solver starts at the identity view either way (see the
    // defaults in `GpuSimulation`'s constructor); resetting the tracked state
    // to match is what keeps a stale zoom from being re-applied to a domain
    // that may just have changed shape.
    resetView();
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
    onStreamlines: (levels, lineW) => sim?.setStreamlines(levels, lineW),
    onMesh: (m) => { if (sim) sim.mesh = MESH[m]; },
    onColormap: (v) => sim?.setColormap(v),
    onNuWindow: (steps) => { nu.setWindow(steps); rms.setWindow(steps); },
    onReseed: () => {
      sim?.reseed(0.05, state.wavenumber);
      nu.clear();
      rms.clear();
    },
    onResolution: (p) => {
      // The dt cap is resolution-dependent (finer grids want a lower ceiling),
      // so a preset carries its own; adopt it and reflect that back into the
      // pane. The Courant number is a property of the accuracy target, not the
      // grid, so it is untouched.
      state.dtMax = PRESETS[p].dtMax;
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
      const { variable, strainRate, tackley, brandenburg } = VISCOSITY[v];
      // Tackley and Brandenburg's pointwise laws are each a different GPU
      // kernel (no `sref` pass, a different Params layout use), so entering or
      // leaving either is a rebuild even though it stays inside the Krylov
      // tier — see `presets.ts`.
      if (!sim || variable !== sim.o.variable || tackley !== sim.o.tackley
          || brandenburg !== sim.o.brandenburg)
        return void build(state);
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
    onSigmaY: (v) => { if (sim) sim.sigmaY = v; },
    onSigmaB: (v) => { if (sim) sim.sigmaB = v; },
    onEtaStar: (v) => { if (sim) sim.etaStar = v; },
    onBrandenburgProfile: (aUpper, aLower, d0) => {
      const s = sim;
      if (s?.o.brandenburg)
        void announce("re-inverting the μ̄(r) radial blocks…", () => {
          if (sim === s) s.setBrandenburgProfile(aUpper, aLower, d0);
        });
    },
    onResetView: () => { resetView(); canvas.style.cursor = "default"; },
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

      // Every frame, not on the FPS cadence below: `pollStats` is guarded by
      // its own `statPending` flag and off the dependency chain either way, and
      // `speed` can run several steps per frame — polling less often would lag
      // the adaptive step by that many frames on top of the poll's own latency.
      sim.pollStats();
      const next = adaptiveDt(sim.dt, state.dtMax, state.courant, sim.stats.maxSpeed);
      if (next !== null) sim.setDt(next);

      if (++frames % 15 === 0) {
        const now = performance.now();
        fps = (15 * 1000) / (now - last);
        last = now;
      }
      // The reductions are one poll behind at worst, and NaN until the first
      // lands — say so rather than printing a number that is not there yet.
      const n = (v: number, d = 4) => (Number.isNaN(v) ? "—" : v.toFixed(d));
      const { nuInner, nuOuter, psiMax, vrms, maxSpeed, at, atStep } = sim.stats;
      // Offered every frame rather than on the poll cadence: `stats` is replaced
      // asynchronously, so there is no frame the host can name as the one a
      // reading arrived on. The trace drops the NaNs before the first readback
      // and the repeats of a sample already held, so this is the same series
      // either way, at most one frame sooner.
      nu.push({ t: at, step: atStep, inner: nuInner, outer: nuOuter });
      rms.push({ t: at, step: atStep, v: vrms });
      const power = sim.o.variable && (sim.n !== 1 || sim.o.tackley);
      // Named only when it is doing something: at c = 0 the depth term is
      // exactly absent, and a "d" in the law with "× 10^0.00" after it would
      // describe a dependence the solver does not have.
      const deep = sim.o.variable && sim.cz !== 0;
      const law = sim.o.tackley
        ? `Tackley η(T, d, ε̇), σ_Y = ${sim.sigmaY.toFixed(2)}, ` +
          `σ_b = ${sim.sigmaB.toFixed(2)}, η* = ${sim.etaStar.toExponential(1)}`
        : sim.o.brandenburg
          ? `Brandenburg η(T, d, ε̇), b = ${sim.gamma.toFixed(2)}, ` +
            `c = ${sim.cz.toFixed(2)}, A₀ = ${sim.aUpper.toFixed(2)}/` +
            `${sim.aLower.toFixed(2)} at d₀ = ${sim.d0.toFixed(3)}, ` +
            `σ_Y = ${sim.sigmaY.toFixed(2)}, σ_b = ${sim.sigmaB.toFixed(2)}, ` +
            `η* = ${sim.etaStar.toExponential(1)}`
          : sim.o.variable
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
        `ψ ${sim.nr}×${sim.na} splines   T ${sim.gnr}×${sim.gna} grid   ` +
        `dt = ${sim.dt.toExponential(1)}   Co = ${n(sim.dt * maxSpeed, 2)}\n` +
        // The scaling behind every dimensional figure below, stated where the
        // rest of the configuration is. Without it "133 Gyr" is a number and not
        // a quantity — and the reference is a display assumption, not something
        // the solver knows (see ui/dimensional.ts).
        `${referenceNote()}\n\n` +
        `step ${String(sim.steps).padStart(6)}   t = ${sim.time.toFixed(4)} = ` +
        `${dimensionalTime(sim.time)}   ` +
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
