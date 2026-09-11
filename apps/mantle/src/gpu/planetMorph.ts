/**
 * Exterior-only bridge between two solver-backed cutaways. It owns no solver
 * buffer: the destination simulation is deliberately not built until this
 * centred whole-planet morph has completed.
 */
import { radiiFor, type PlanetDefinition } from "../planets";
import type { SurfaceTexture } from "./surfaceAssets";
import { ease } from "./globe";

const source = /* wgsl */`
struct U { blend: f32, srcRadius: f32, dstRadius: f32, srcKind: f32,
  dstKind: f32, eyeAz: f32, eyeEl: f32, eyeDist: f32,
  srcTilt: f32, dstTilt: f32, _pad0: f32, _pad1: f32, }
@group(0) @binding(0) var<uniform> u: U;
@group(0) @binding(1) var srcTex: texture_2d<f32>;
@group(0) @binding(2) var srcSamp: sampler;
@group(0) @binding(3) var dstTex: texture_2d<f32>;
@group(0) @binding(4) var dstSamp: sampler;
struct Out { @builtin(position) pos: vec4f, @location(0) uv: vec2f, }
@vertex fn vs(@builtin(vertex_index) i: u32) -> Out {
  let p = array<vec2f, 3>(vec2f(-1., -1.), vec2f(3., -1.), vec2f(-1., 3.));
  var o: Out; o.pos = vec4f(p[i], 0., 1.); o.uv = p[i] * .5 + .5; return o;
}
fn venusSurface(n: vec3f) -> vec3f {
  let bands = .08 * sin(n.y * 18. + n.x * 9.);
  return vec3f(.72 + bands, .31 + bands * .45, .07);
}
fn sampled(tex: texture_2d<f32>, samp: sampler, n: vec3f, kind: f32, tilt: f32) -> vec3f {
  // Identical to Globe3D surface UV: north remains at the texture's top
  // throughout the exterior morph instead of flipping between renderers.
  let tilted = vec3f(n.x, cos(tilt) * n.y - sin(tilt) * n.z,
    sin(tilt) * n.y + cos(tilt) * n.z);
  let longitude = atan2(-tilted.x, -tilted.z);
  let latitude = asin(clamp(tilted.y, -1., 1.));
  let mapped = textureSampleLevel(tex, samp,
    vec2f(longitude / 6.2831853 + .5, .5 - latitude / 3.1415926), 0.).rgb;
  return select(mapped, venusSurface(n), kind > .5);
}
@fragment fn fs(in: Out) -> @location(0) vec4f {
  let radius = mix(u.srcRadius, u.dstRadius, u.blend);
  var bg = vec3f(.008, .009, .025) + .018 * vec3f(in.uv.y);
  // Reconstruct the same world-space surface normal as Globe3D's camera.
  // Keeping this pose (including reader orbit) makes the image stay put at
  // the exact frame the live globe hands over to the exterior bridge.
  let ca = cos(u.eyeAz); let sa = sin(u.eyeAz);
  let ce = cos(u.eyeEl); let se = sin(u.eyeEl);
  let eyeDir = vec3f(ce * sa, se, ce * ca);
  let forward = -eyeDir;
  let right = normalize(cross(forward, vec3f(0., 1., 0.)));
  let up = cross(right, forward);
  // Same perspective ray and sphere intersection as Globe3D. This avoids a
  // texture swim at handoff even when the reader has orbited or zoomed.
  let ndc = in.uv * 2. - 1.;
  let rayOrigin = eyeDir * u.eyeDist;
  let fov = tan(.5 * .80285146); // Globe3D HERO.fov = 46 degrees
  let rayDir = normalize(forward + ndc.x * fov * right + ndc.y * fov * up);
  let b = dot(rayOrigin, rayDir);
  let disc = b * b - (dot(rayOrigin, rayOrigin) - radius * radius);
  if (disc <= 0.) { return vec4f(bg, 1.); }
  let n = normalize(rayOrigin + (-b - sqrt(disc)) * rayDir);
  let light = .30 + .70 * max(dot(n, normalize(vec3f(.45, .55, .70))), 0.);
  let src = sampled(srcTex, srcSamp, n, u.srcKind, u.srcTilt);
  let dst = sampled(dstTex, dstSamp, n, u.dstKind, u.dstTilt);
  return vec4f(mix(src, dst, u.blend) * light, 1.);
}`;

export class PlanetMorphScene {
  private readonly uniform: GPUBuffer;
  private readonly pipeline: GPURenderPipeline;
  private readonly bind: GPUBindGroup;
  private readonly data = new Float32Array(12);
  private blend = 0;
  constructor(private readonly device: GPUDevice, format: GPUTextureFormat,
    sourceTexture: SurfaceTexture, destinationTexture: SurfaceTexture,
    sourcePlanet: PlanetDefinition, destinationPlanet: PlanetDefinition,
    orientation: readonly [number, number, number] = [.55, .32, 6.5]) {
    this.uniform = device.createBuffer({ size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    const module = device.createShaderModule({ code: source });
    this.pipeline = device.createRenderPipeline({ layout: "auto", vertex: { module, entryPoint: "vs" }, fragment: { module, entryPoint: "fs", targets: [{ format }] } });
    this.bind = device.createBindGroup({ layout: this.pipeline.getBindGroupLayout(0), entries: [
      { binding: 0, resource: { buffer: this.uniform } },
      { binding: 1, resource: sourceTexture.view }, { binding: 2, resource: sourceTexture.sampler },
      { binding: 3, resource: destinationTexture.view }, { binding: 4, resource: destinationTexture.sampler },
    ] });
    const radius = (p: PlanetDefinition) => radiiFor(p).ro;
    this.data.set([0, radius(sourcePlanet), radius(destinationPlanet),
      sourcePlanet.visual.surface === "venus-procedural" ? 1 : 0,
      destinationPlanet.visual.surface === "venus-procedural" ? 1 : 0,
      orientation[0], orientation[1], orientation[2],
      sourcePlanet.visual.axialTiltDeg * Math.PI / 180,
      destinationPlanet.visual.axialTiltDeg * Math.PI / 180]);
  }
  setBlend(value: number): void { this.blend = ease(Math.min(1, Math.max(0, value))); }
  draw(view: GPUTextureView): void {
    this.data[0] = this.blend;
    this.device.queue.writeBuffer(this.uniform, 0, this.data);
    const enc = this.device.createCommandEncoder();
    const pass = enc.beginRenderPass({ colorAttachments: [{ view, clearValue: { r: .008, g: .009, b: .025, a: 1 }, loadOp: "clear", storeOp: "store" }] });
    pass.setPipeline(this.pipeline); pass.setBindGroup(0, this.bind); pass.draw(3); pass.end();
    this.device.queue.submit([enc.finish()]);
  }
  destroy(): void { this.uniform.destroy(); }
}
