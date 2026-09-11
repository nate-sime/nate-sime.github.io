/** Lightweight, deliberately not-to-scale exterior-only travel stage. */
import type { PlanetDefinition } from "../planets";
import { ease } from "./globe";

export type SolarShot = "departure" | "overview" | "arrival";

const source = /* wgsl */`
struct U { progress: f32, shot: f32, departure: vec4f, arrival: vec4f, }
@group(0) @binding(0) var<uniform> u: U;
struct Out { @builtin(position) pos: vec4f, @location(0) uv: vec2f, }
@vertex fn vs(@builtin(vertex_index) i: u32) -> Out {
  let p = array<vec2f, 3>(vec2f(-1., -1.), vec2f(3., -1.), vec2f(-1., 3.));
  var o: Out; o.pos = vec4f(p[i], 0., 1.); o.uv = p[i] * .5 + .5; return o;
}
fn disc(uv: vec2f, c: vec2f, r: f32, colour: vec3f) -> vec3f {
  let d = length(uv - c); let edge = smoothstep(r, r - .008, d);
  let light = clamp(1.15 - d / max(r, .001) * .55 + (uv.y - c.y) * .35 / max(r, .001), .25, 1.2);
  return colour * edge * light;
}
@fragment fn fs(in: Out) -> @location(0) vec4f {
  var col = vec3f(.008, .009, .025) + .018 * vec3f(in.uv.y, in.uv.y, in.uv.y);
  let k = u.progress;
  if (u.shot == 1.) {
    col += disc(in.uv, vec2f(.18, .57), .12, vec3f(1.0, .55, .12));
    col += disc(in.uv, vec2f(.51, .52), .035, u.departure.rgb);
    col += disc(in.uv, vec2f(.78, .41), .031, u.arrival.rgb);
  } else {
    let from = mix(vec2f(.5), vec2f(.20, .52), k);
    let to = mix(vec2f(.80, .42), vec2f(.5), k);
    let c = select(from, to, u.shot == 2.);
    let colour = select(u.departure.rgb, u.arrival.rgb, u.shot == 2.);
    col += disc(in.uv, c, mix(.44, .08, k), colour);
  }
  return vec4f(col, 1.);
}`;

export class SolarSystemScene {
  private readonly uniform: GPUBuffer;
  private readonly pipeline: GPURenderPipeline;
  private readonly bind: GPUBindGroup;
  private readonly data = new Float32Array(12);
  private shot: SolarShot = "departure";
  private progress = 0;
  private from: PlanetDefinition | null = null;
  private to: PlanetDefinition | null = null;
  constructor(private readonly device: GPUDevice, format: GPUTextureFormat) {
    this.uniform = device.createBuffer({ size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    const module = device.createShaderModule({ code: source });
    this.pipeline = device.createRenderPipeline({ layout: "auto", vertex: { module, entryPoint: "vs" }, fragment: { module, entryPoint: "fs", targets: [{ format }] } });
    this.bind = device.createBindGroup({ layout: this.pipeline.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: this.uniform } }] });
  }
  show(shot: SolarShot, progress: number, from: PlanetDefinition, to: PlanetDefinition): void {
    this.shot = shot; this.progress = ease(Math.min(1, Math.max(0, progress))); this.from = from; this.to = to;
  }
  draw(view: GPUTextureView): void {
    if (!this.from || !this.to) return;
    const tint = (p: PlanetDefinition) => p.visual.surface === "venus-procedural" ? [0.94, .58, .16] : [.20, .48, .82];
    this.data.set([this.progress, this.shot === "overview" ? 1 : this.shot === "arrival" ? 2 : 0, ...tint(this.from), 1, ...tint(this.to), 1]);
    this.device.queue.writeBuffer(this.uniform, 0, this.data);
    const e = this.device.createCommandEncoder(); const pass = e.beginRenderPass({ colorAttachments: [{ view, clearValue: { r: .008, g: .009, b: .025, a: 1 }, loadOp: "clear", storeOp: "store" }] });
    pass.setPipeline(this.pipeline); pass.setBindGroup(0, this.bind); pass.draw(3); pass.end(); this.device.queue.submit([e.finish()]);
  }
  destroy(): void { this.uniform.destroy(); }
}
