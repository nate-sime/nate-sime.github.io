/**
 * The 3D cutaway-globe view — a purely cosmetic alternate render of the same
 * (r, φ) field `GpuSimulation.draw` shows flat, reached by a button and an
 * animated transition rather than a second simulation. See PLAN-3d-view.md.
 *
 * Built and torn down against a *live* simulation, the way `GpuParticles` is
 * (`new Globe3D(sim, colormap)` / `globe.destroy()`), rather than being a
 * construction-time option of `GpuSimulation` itself — it borrows `params`
 * and `T` from the host and owns nothing the solver depends on, so it can be
 * discarded and rebuilt on every `main.ts` rebuild for free.
 *
 * **The transition is a camera move, not a cross-fade.** `progress` runs 0→1
 * (2D pose → hero pose) and drives everything linearly: `az`/`el` orbit in,
 * `persp` blends the ray/projection formulas in `wgsl.ts` from orthographic
 * to perspective, and `reveal` fades in the exterior shell and core. The cut
 * faces are not faded at all — they read the same `T` buffer the flat view
 * does, at `progress = 0` framed to match it exactly (see `toggle`), so the
 * one thing that must look continuous across the whole animation always does.
 */

import type { ColormapName } from "../colormaps";
import type { EarthTexture } from "./earthTexture";
import type { GpuParticles } from "./particles";
import * as w from "./wgsl";

/** The slice of `GpuSimulation` this class needs — satisfied structurally, the way `GpuParticleHost` is. */
export interface GlobeHost {
  readonly device: GPUDevice;
  readonly format: GPUTextureFormat;
  buffer(name: string): GPUBuffer;
}

/** The 2D view's own screen↔world map, at the moment `toggle` is called — what `progress = 0` reproduces. */
export interface View2D {
  halfExtent: number;
  zoom: number;
  panX: number;
  panY: number;
}

const DURATION_MS = 1300;
const WEDGE_W = (100 * Math.PI) / 180;
const HERO = { az: 0.55, el: 0.32, dist: 6.5, fov: (46 * Math.PI) / 180 };
const DIST_MIN = 3.2, DIST_MAX = 14;
const EL_LIMIT = 1.35;

/** Smoothstep, doing double duty as an easing curve — the same shape `meshLine`'s fades use. */
const ease = (t: number): number => t * t * (3 - 2 * t);

export class Globe3D {
  private readonly device: GPUDevice;
  private readonly buf: Record<string, GPUBuffer> = {};
  private scenePipe!: GPURenderPipeline;
  private sceneBind!: GPUBindGroup;
  private particlePipe!: GPURenderPipeline;
  private particleBind: GPUBindGroup | null = null;
  private lastParticles: GpuParticles | null = null;
  private depthTex: GPUTexture | null = null;
  private depthView: GPUTextureView | null = null;
  private side = 1;

  private mode: "2d" | "3d" = "2d";
  private animating = false;
  private animStart = 0;
  private fromProgress = 0;
  private target = 0;
  private progress = 0;   // 0 = 2D-matching pose, 1 = hero pose

  private az = 0;
  private el = 0;
  private userAz = 0;
  private userEl = 0;
  private dist = HERO.dist;
  private readonly fov = HERO.fov;
  private readonly wedgeW = WEDGE_W;
  private orthoHalf = 1;
  private panX = 0;
  private panY = 0;

  private readonly params = new ArrayBuffer(48);
  private readonly gf = new Float32Array(this.params);

  /**
   * `earth` is owned by `main.ts`, not this class — it is fetched and
   * decoded once for the app's whole lifetime (`earthTexture.ts`) and handed
   * to every `Globe3D` a rebuild constructs, the same borrowed-not-duplicated
   * relationship this class already has with `host`'s `T` buffer. `destroy`
   * below never touches it.
   */
  constructor(private readonly host: GlobeHost, colormap: ColormapName, private readonly earth: EarthTexture) {
    this.device = host.device;
    this.buf.globe = this.device.createBuffer({
      size: this.params.byteLength,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.buildScenePipeline(colormap);
    this.buildParticlePipeline(colormap, false);
    this.syncGlobe();
  }

  get viewMode(): "2d" | "3d" { return this.mode; }
  get inTransition(): boolean { return this.animating; }

  private buildScenePipeline(colormap: ColormapName): void {
    const module = this.device.createShaderModule({ code: w.globeSource(colormap) });
    this.scenePipe = this.device.createRenderPipeline({
      layout: "auto",
      vertex: { module, entryPoint: "vs" },
      fragment: { module, entryPoint: "fs", targets: [{ format: this.host.format }] },
      depthStencil: { format: "depth32float", depthWriteEnabled: true, depthCompare: "less" },
    });
    this.sceneBind = this.device.createBindGroup({
      layout: this.scenePipe.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.host.buffer("params") } },
        { binding: 1, resource: { buffer: this.buf.globe } },
        { binding: 2, resource: { buffer: this.host.buffer("T") } },
        { binding: 3, resource: this.earth.view },
        { binding: 4, resource: this.earth.sampler },
      ],
    });
  }

  private buildParticlePipeline(colormap: ColormapName, discrete: boolean): void {
    const module = this.device.createShaderModule({ code: w.globeParticleSource(colormap, discrete) });
    this.particlePipe = this.device.createRenderPipeline({
      layout: "auto",
      vertex: { module, entryPoint: "vs" },
      fragment: {
        module, entryPoint: "fs",
        targets: [{
          format: this.host.format,
          blend: {
            color: { srcFactor: "src-alpha", dstFactor: "one-minus-src-alpha", operation: "add" },
            alpha: { srcFactor: "one", dstFactor: "one-minus-src-alpha", operation: "add" },
          },
        }],
      },
      primitive: { topology: "triangle-strip" },
      // No depth write: tracers never need to occlude one another, only the
      // scene pass's shell/core/cut-face surfaces — same call `particles.ts`
      // §3 makes for the flat view's own draw-order blending.
      depthStencil: { format: "depth32float", depthWriteEnabled: false, depthCompare: "less" },
    });
    // The bind group is rebuilt lazily against whatever `GpuParticles` is
    // attached when `draw` next runs — invalidate it so a colormap change
    // does not draw one frame through a stale pipeline's old layout.
    this.particleBind = null;
    this.lastParticles = null;
  }

  /** Swap the temperature colour map — a pipeline rebuild, exactly `GpuSimulation.setColormap`'s own reason: the control points are compiled into the shader, not read from a buffer. */
  setColormap(colormap: ColormapName, discrete = false): void {
    this.buildScenePipeline(colormap);
    this.buildParticlePipeline(colormap, discrete);
  }

  /** Depth buffer sized to the canvas — call from the same `ResizeObserver` that sizes the canvas itself. */
  setViewport(side: number): void {
    if (side === this.side && this.depthTex) return;
    this.side = side;
    this.depthTex?.destroy();
    this.depthTex = this.device.createTexture({
      size: [side, side],
      format: "depth32float",
      usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
    this.depthView = this.depthTex.createView();
  }

  /**
   * Start (or reverse) the transition. Going toward 3D captures the 2D
   * view's own `halfExtent`/zoom/pan, so `progress = 0`'s orthographic ray
   * origins (`rayOrigin` in `wgsl.ts`) reproduce exactly the screen↔world map
   * the flat canvas was just showing — see this class's own header.
   *
   * `instant` skips the tween and lands straight on the target pose — for
   * `main.ts`'s rebuild, which opens straight in 3D and has no prior flat
   * frame on screen for an animated pan-out to read as continuous with.
   */
  toggle(from2D: View2D, instant = false): void {
    this.mode = this.mode === "2d" ? "3d" : "2d";
    this.target = this.mode === "3d" ? 1 : 0;
    this.fromProgress = this.progress;
    this.animStart = performance.now();
    this.animating = !instant;
    if (instant) this.progress = this.target;
    if (this.mode === "3d") {
      this.orthoHalf = from2D.halfExtent / from2D.zoom;
      this.panX = from2D.panX;
      this.panY = from2D.panY;
    }
  }

  /** Orbit by a screen-space drag delta (NDC units) — only meaningful once settled in 3D; see `main.ts`'s pointer wiring. */
  orbit(dxNdc: number, dyNdc: number): void {
    this.userAz -= dxNdc * 2.4;
    this.userEl = Math.min(EL_LIMIT, Math.max(-EL_LIMIT, this.userEl + dyNdc * 2.4));
  }

  /** Dolly by a wheel delta — same sign convention as `main.ts`'s 2D zoom handler. */
  zoomBy(deltaY: number): void {
    this.dist = Math.min(DIST_MAX, Math.max(DIST_MIN, this.dist * Math.exp(deltaY * 0.0015)));
  }

  /** Advance the transition tween and refresh the uniform — call once per frame, whether or not `mode` is settled. */
  tick(now: number): void {
    if (this.animating) {
      const frac = Math.min(1, (now - this.animStart) / DURATION_MS);
      this.progress = this.fromProgress + (this.target - this.fromProgress) * ease(frac);
      if (frac >= 1) { this.animating = false; this.progress = this.target; }
    }
    this.az = (HERO.az + this.userAz) * this.progress;
    this.el = (HERO.el + this.userEl) * this.progress;
    this.syncGlobe();
  }

  private syncGlobe(): void {
    this.gf.set([
      this.az, this.el, this.dist, this.fov,
      this.wedgeW, this.progress, this.orthoHalf, this.panX,
      this.panY, this.progress, this.earth.available ? 1 : 0, 0,
    ]);
    this.device.queue.writeBuffer(this.buf.globe, 0, this.gf);
  }

  private bindParticles(particles: GpuParticles): void {
    if (particles === this.lastParticles && this.particleBind) return;
    this.lastParticles = particles;
    this.particleBind = this.device.createBindGroup({
      layout: this.particlePipe.getBindGroupLayout(0),
      entries: [this.buf.globe, particles.buffer("part"), particles.buffer("Prt")]
        .map((buffer, binding) => ({ binding, resource: { buffer } })),
    });
  }

  /** Draw the scene, and the tracer overlay over it if one is attached — mirrors `GpuSimulation.draw`. */
  draw(view: GPUTextureView, particles: GpuParticles | null): void {
    if (!this.depthView) return;
    const enc = this.device.createCommandEncoder();
    const pass = enc.beginRenderPass({
      colorAttachments: [{
        view, clearValue: { r: 0.02, g: 0.02, b: 0.047, a: 1 },
        loadOp: "clear", storeOp: "store",
      }],
      depthStencilAttachment: {
        view: this.depthView, depthClearValue: 1, depthLoadOp: "clear", depthStoreOp: "store",
      },
    });
    pass.setPipeline(this.scenePipe);
    pass.setBindGroup(0, this.sceneBind);
    pass.draw(3);
    if (particles) {
      this.bindParticles(particles);
      pass.setPipeline(this.particlePipe);
      pass.setBindGroup(0, this.particleBind!);
      pass.draw(4, particles.count);
    }
    pass.end();
    this.device.queue.submit([enc.finish()]);
  }

  /** Release every GPU allocation this class owns. Required before building a replacement (`main.ts`'s `build`, on every rebuild). */
  destroy(): void {
    for (const b of Object.values(this.buf)) b.destroy();
    this.depthTex?.destroy();
  }
}
