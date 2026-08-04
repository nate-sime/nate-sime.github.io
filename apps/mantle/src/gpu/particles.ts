/**
 * The tracer cloud, resident on the GPU: a pathline overlay when it is only
 * drawn, and a marker-in-cell discretisation of a chemically distinct layer
 * when its composition is also projected into the buoyancy load. Everything
 * here is the validated CPU twin (`solver/particles.ts`) run on the device —
 * same seeding, same RK2 push, same cloud-in-cell projection — so that a
 * cloud started from the same named seed traces the same pathlines and
 * carries the same dense layer either way.
 *
 * A `GpuParticles` is built and torn down against a *running* simulation
 * (`sim.particles = new GpuParticles(sim, opts)`, `sim.particles?.destroy()`)
 * rather than being a construction-time option of `GpuSimulation` itself.
 * Turning the tracer overlay on or off is then a few tens of milliseconds —
 * three small shader modules, not the multi-second f64 factorisation a
 * resolution or viscosity-tier change pays for — and every buffer this class
 * owns is discarded and rebuilt on `destroy`, so a stale cloud can never
 * survive past the run that seeded it.
 *
 * Two things are *borrowed* from the host simulation rather than owned here:
 * the velocity field (`knots`, `psi`) a tracer is pushed through, and the
 * step's own dt (`passA` — the same buffer the temperature advection reads,
 * so a tracer's path and the field it is drawn over always advance through
 * the same instant of the flow). The one thing that flows the other way is
 * the composition grid: when this class is attached, the host's buoyancy
 * load is repointed at *this* class's `part`/`C` buffers
 * (`GpuParticleHost.bindBuoyancy`) so that `T − buoy·C` reads a live
 * projection; detaching repoints it at the host's own zeroed placeholders,
 * so a thermal-only run never depends on this file having run at all.
 */

import type { Geometry } from "../geometry";
import {
  seedParticles, depthOf, DEFAULT_LAYER_DEPTH, TARGET_PARTICLES_PER_CELL,
  SPECIES_CONDITIONS, PARTICLE_TINT,
  type SpeciesConditionName, type TintMode, type TintEntry,
} from "../particles";
import * as w from "./wgsl";

const S = Float32Array.BYTES_PER_ELEMENT;

export interface GpuParticleOptions {
  /** Tracer count. Defaults to `TARGET_PARTICLES_PER_CELL` per composition-grid cell. */
  count?: number;
  /** Names the draw — the same seed reproduces the same cloud, on the CPU or the GPU. */
  seed?: number;
  /** Initial composition profile; see `SPECIES_CONDITIONS` in `particles.ts`. */
  species?: SpeciesConditionName;
  /** Thickness of the dense layer, as a fraction of the mantle's depth. */
  layerDepth?: number;
  /** Composition-grid resolution. Default: half the temperature grid, each direction. */
  cnr?: number;
  cna?: number;
  /** How a tracer is coloured on screen; see `PARTICLE_TINT`. Defaults to *initial depth*. */
  tint?: TintMode;
  /** Dot radius, in screen pixels. */
  radius?: number;
  /** Dot opacity, 0–1. */
  opacity?: number;
  /** Time scale the *age* colour mode saturates against — nondimensional simulation time. */
  tau?: number;
}

/**
 * The handful of resources `GpuParticles` needs from a live simulation:
 * enough to push a tracer through the current flow at the current step's dt,
 * to read the colour modes that track a tracer's surroundings, and to hand
 * the buoyancy load a place to read its composition back from.
 * `GpuSimulation` satisfies this without saying so — structurally, the way
 * every bind group in this codebase is matched by position rather than by name.
 */
export interface GpuParticleHost {
  readonly device: GPUDevice;
  readonly o: { geom: Geometry; gnr: number; gna: number };
  /** The canvas/texture format particles must draw into — the same one the host's own render pipeline targets. */
  readonly format: GPUTextureFormat;
  /** The current step's dt (`passA`'s own value) — what the *age* colour mode's clock advances by. */
  readonly dt: number;
  buffer(name: string): GPUBuffer;
  bindBuoyancy(part: GPUBuffer, C: GPUBuffer): void;
  unbindBuoyancy(): void;
}

/** The one-time colour a static tint mode paints on at seeding. Dynamic modes (temperature, speed, age) are left at 0 for the first push to fill in. */
const seededTint = (mode: TintMode, geom: Geometry, r: number, phi: number, c: number): number => {
  switch (mode) {
    case "uniform": return 0;
    case "initial depth": return depthOf(geom, r);
    case "initial φ": return phi / geom.span;
    case "species": return c;   // the coupled mode: its colour *is* the composition
    default: return 0;
  }
};

export class GpuParticles {
  readonly count: number;
  readonly cnr: number;
  readonly cna: number;
  private readonly device: GPUDevice;
  private readonly geom: Geometry;
  private readonly tint: TintMode;
  private readonly species: SpeciesConditionName;
  private readonly layerDepth: number;
  private readonly buf: Record<string, GPUBuffer> = {};
  private readonly pipe: Record<string, GPUComputePipeline> = {};
  private readonly bind: Record<string, GPUBindGroup> = {};
  private render!: GPURenderPipeline;
  private renderBind!: GPUBindGroup;
  /** Elapsed simulation time since the cloud was last (re)seeded — what the *age* colour mode reads, normalised by `tau`. */
  private age = 0;
  private readonly tau: number;
  private seedValue: number;
  private readonly partData: ArrayBuffer;
  private readonly partInt: Int32Array;
  private readonly partFloat: Float32Array;

  constructor(private readonly host: GpuParticleHost, o: GpuParticleOptions = {}) {
    const { device } = host;
    this.device = device;
    this.geom = host.o.geom;
    this.cnr = o.cnr ?? Math.floor((host.o.gnr + 1) / 2);
    this.cna = o.cna ?? host.o.gna / 2;
    this.count = o.count ?? TARGET_PARTICLES_PER_CELL * this.cnr * this.cna;
    this.tint = o.tint ?? "initial depth";
    this.species = o.species ?? "dense basal layer";
    this.layerDepth = o.layerDepth ?? DEFAULT_LAYER_DEPTH;
    this.tau = o.tau ?? 1;
    this.seedValue = o.seed ?? 1;

    const ST = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST;
    this.buf.Prt = device.createBuffer({ size: this.count * 4 * S, usage: ST });
    this.buf.num = device.createBuffer({ size: this.cnr * this.cna * S, usage: ST });
    this.buf.den = device.createBuffer({ size: this.cnr * this.cna * S, usage: ST });
    this.buf.C = device.createBuffer({ size: this.cnr * this.cna * S, usage: ST });

    // `Part` — see wgsl.ts's own header on why this is a small uniform of its
    // own rather than a few more fields on the shared `PARAMS` block: none of
    // it (grid shape, dot radius/opacity, canvas side, the age clock) has
    // anything to do with the temperature/Stokes pipelines it sits beside.
    this.partData = new ArrayBuffer(48);
    this.partInt = new Int32Array(this.partData, 0, 4);
    this.partFloat = new Float32Array(this.partData, 16, 8);
    this.partInt.set([this.count, this.cnr, this.cna, 0]);
    this.partFloat.set([o.radius ?? 2.5, o.opacity ?? 0.85, 1, 0, this.tau, 0, 0, 0]);
    this.buf.part = device.createBuffer({
      size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST,
    });
    this.syncPart();

    this.buildComputePipelines();
    this.buildRenderPipeline();
    host.bindBuoyancy(this.buf.part, this.buf.C);
    this.seed(this.seedValue);
  }

  private syncPart(): void {
    this.device.queue.writeBuffer(this.buf.part, 0, this.partData);
  }

  private buildComputePipelines(): void {
    const { device, geom } = this;
    const tintEntry: TintEntry = PARTICLE_TINT[this.tint];
    const wantsT = (tintEntry.wgsl ?? "").includes("sample_T");
    const wantsStat = (tintEntry.wgsl ?? "").includes("stat[");

    const compute = (code: string) => device.createComputePipeline({
      layout: "auto",
      compute: { module: device.createShaderModule({ code }), entryPoint: "main" },
    });
    const group = (pipeline: GPUComputePipeline, res: GPUBuffer[]) =>
      device.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: res.map((buffer, binding) => ({ binding, resource: { buffer } })),
      });

    this.pipe.push = compute(w.pushSource(geom, tintEntry));
    const pushRes = [
      this.host.buffer("params"), this.buf.part, this.host.buffer("passA"),
      this.host.buffer("knots"), this.host.buffer("psi"), this.buf.Prt,
    ];
    if (wantsT) pushRes.push(this.host.buffer("T"));
    else if (wantsStat) pushRes.push(this.host.buffer("stat"));
    this.bind.push = group(this.pipe.push, pushRes);

    this.pipe.clear = compute(w.clearSource());
    this.bind.clear = group(this.pipe.clear, [this.buf.part, this.buf.num, this.buf.den]);

    this.pipe.scatter = compute(w.scatterSource());
    this.bind.scatter = group(this.pipe.scatter,
      [this.host.buffer("params"), this.buf.part, this.buf.Prt, this.buf.num, this.buf.den]);

    this.pipe.normalise = compute(w.normaliseSource(geom));
    this.bind.normalise = group(this.pipe.normalise,
      [this.buf.part, this.buf.num, this.buf.den, this.buf.C]);
  }

  private buildRenderPipeline(): void {
    const { device } = this;
    const entry: TintEntry = PARTICLE_TINT[this.tint];
    const module = device.createShaderModule({
      code: w.particleRenderSource(this.geom, entry.colormap, entry.discrete ?? false),
    });
    this.render = device.createRenderPipeline({
      layout: "auto",
      vertex: { module, entryPoint: "vs" },
      fragment: {
        module, entryPoint: "fs",
        targets: [{
          format: this.host.format,
          // Order-dependent, and the order is buffer order — arbitrary, and
          // invisible at a 2–3 px dot (see `particles.ts` §3 in the plan).
          blend: {
            color: { srcFactor: "src-alpha", dstFactor: "one-minus-src-alpha", operation: "add" },
            alpha: { srcFactor: "one", dstFactor: "one-minus-src-alpha", operation: "add" },
          },
        }],
      },
      primitive: { topology: "triangle-strip" },
    });
    this.renderBind = device.createBindGroup({
      layout: this.render.getBindGroupLayout(0),
      entries: [this.host.buffer("params"), this.buf.part, this.buf.Prt]
        .map((buffer, binding) => ({ binding, resource: { buffer } })),
    });
  }

  private dispatch(pass: GPUComputePassEncoder, name: string, n: number): void {
    pass.setPipeline(this.pipe[name]);
    pass.setBindGroup(0, this.bind[name]);
    pass.dispatchWorkgroups(Math.ceil(n / w.WG));
  }

  /**
   * Advance every tracer one step along the flow that produced the current
   * ψ, and refresh the colour modes that track a tracer's present
   * surroundings. Also advances the *age* colour mode's clock by the host's
   * current dt — cheap enough (one 48-byte uniform write) to do unconditionally
   * rather than asking the caller whether that mode happens to be selected.
   * Pass a live compute pass to dispatch inline with a step already in
   * progress (`GpuSimulation.step` — see its own ordering note on why this
   * belongs *before* the Stokes re-solve); call with no argument to push once
   * on its own, self-submitted encoder.
   */
  push(pass?: GPUComputePassEncoder): void {
    this.age += this.host.dt;
    this.partFloat[3] = this.age;   // pt.age
    this.syncPart();
    if (pass) { this.dispatch(pass, "push", this.count); return; }
    const enc = this.device.createCommandEncoder();
    const p = enc.beginComputePass();
    this.dispatch(p, "push", this.count);
    p.end();
    this.device.queue.submit([enc.finish()]);
  }

  /**
   * Project the cloud's composition onto the (coarser) composition grid —
   * clear the fixed-point accumulators, scatter every tracer's `c` by
   * cloud-in-cell weight, then divide back out. Three dispatches, always in
   * that order, since the second depends on the first having zeroed the
   * accumulators and the third on the second having filled them. Same
   * calling convention as `push`: a pass to embed in, or none to run alone.
   */
  project(pass?: GPUComputePassEncoder): void {
    if (pass) {
      this.dispatch(pass, "clear", this.cnr * this.cna);
      this.dispatch(pass, "scatter", this.count);
      this.dispatch(pass, "normalise", this.cnr * this.cna);
      return;
    }
    const enc = this.device.createCommandEncoder();
    const p = enc.beginComputePass();
    this.dispatch(p, "clear", this.cnr * this.cna);
    this.dispatch(p, "scatter", this.count);
    this.dispatch(p, "normalise", this.cnr * this.cna);
    p.end();
    this.device.queue.submit([enc.finish()]);
  }

  /** Draw every tracer into an already-open render pass, over whatever it already holds — see `GpuSimulation.draw`. */
  draw(pass: GPURenderPassEncoder): void {
    pass.setPipeline(this.render);
    pass.setBindGroup(0, this.renderBind);
    pass.draw(4, this.count);
  }

  /**
   * Draw a fresh cloud, area-uniform, painted with the species profile and
   * the seeded colour modes — what pressing "reseed" means for the tracer
   * overlay, and what the constructor itself does once. Projects immediately
   * (self-submitted; see `project`), so the buoyancy load never reads a
   * composition grid belonging to the previous cloud.
   */
  seed(seedValue = this.seedValue): void {
    this.seedValue = seedValue;
    this.age = 0;
    this.partFloat[3] = 0;   // pt.age
    this.syncPart();

    const { geom, count } = this;
    const drawn = seedParticles(geom, count, seedValue);
    const cond = SPECIES_CONDITIONS[this.species];
    const packed = new Float32Array(count * 4);
    for (let i = 0; i < count; i++) {
      const r = drawn.r[i], phi = drawn.phi[i];
      const c = cond.composition(geom, r, this.layerDepth);
      packed[4 * i] = r;
      packed[4 * i + 1] = phi;
      packed[4 * i + 2] = seededTint(this.tint, geom, r, phi, c);
      packed[4 * i + 3] = c;
    }
    this.device.queue.writeBuffer(this.buf.Prt, 0, packed);
    this.project();
  }

  /** The canvas side length in pixels, so the vertex shader can turn a pixel radius into NDC — see `GpuSimulation`'s resize hook. */
  setViewport(side: number): void {
    this.partFloat[2] = side;
    this.syncPart();
  }

  /** Dot radius (px) and opacity — a pure uniform write, like `GpuSimulation.setView`. */
  setStyle(radius: number, opacity: number): void {
    this.partFloat[0] = radius;
    this.partFloat[1] = opacity;
    this.syncPart();
  }

  /** Copy one of this cloud's own buffers to the host — verification only, as `GpuSimulation.read` is for the main pipeline. */
  async read(name: string): Promise<Float32Array> {
    const src = this.buf[name];
    const dst = this.device.createBuffer({
      size: src.size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });
    const enc = this.device.createCommandEncoder();
    enc.copyBufferToBuffer(src, 0, dst, 0, src.size);
    this.device.queue.submit([enc.finish()]);
    await dst.mapAsync(GPUMapMode.READ);
    const out = new Float32Array(dst.getMappedRange().slice(0));
    dst.unmap();
    dst.destroy();
    return out;
  }

  /** Release every buffer this class owns and hand the buoyancy load back to the host's zeroed placeholders. */
  destroy(): void {
    for (const b of Object.values(this.buf)) b.destroy();
    this.host.unbindBuoyancy();
  }
}
