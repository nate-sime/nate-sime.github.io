// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Exterior-only bridge between two solver-backed cutaways. It owns no solver
 * buffer: the destination simulation is deliberately not built until this
 * centred whole-planet morph has completed.
 */
import { radiiFor, type PlanetDefinition } from "../planets";
import { SURFACE_MATERIALS, type SurfaceTexture } from "./surfaceAssets";
import { ease, sunDirection } from "./globe";
import { milkyWayCube } from "./sky";
import { FALLBACK_SURFACE, MILKY_WAY, NOISE3, SPACE } from "./wgsl";

/**
 * A planet's atmosphere as the shader's three vec4s — (colour, strength),
 * (twilight, opacity), (thickness, reach) — zero for an airless body, which
 * is what Globe3D's own uniform carries for one.
 */
const atmosphereOf = (planet: PlanetDefinition): number[][] => {
  const a = planet.visual.atmosphere;
  return a
    ? [[...a.color, a.strength], [...a.twilight, a.opacity], [a.thickness, a.twilightReach, 0, 0]]
    : [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
};

const fallbackKind = (planet: PlanetDefinition): number =>
  planet.visual.surface.startsWith("venus") ? 1 :
  planet.visual.surface.startsWith("mars") ? 2 : 3;

const source = NOISE3 + MILKY_WAY + SPACE + FALLBACK_SURFACE + /* wgsl */`
struct U { blend: f32, srcRadius: f32, dstRadius: f32, srcKind: f32,
  dstKind: f32, eyeAz: f32, eyeEl: f32, eyeDist: f32,
  srcTilt: f32, dstTilt: f32, srcOffset: f32, dstOffset: f32,
  srcFeather: f32, dstFeather: f32, pixelUv: f32, _pad1: f32,
  srcAtmosphere: vec4f, dstAtmosphere: vec4f,
  srcTwilight: vec4f, dstTwilight: vec4f,
  srcShape: vec4f, dstShape: vec4f,
  srcSun: vec4f, dstSun: vec4f, }
@group(0) @binding(0) var<uniform> u: U;
@group(0) @binding(1) var srcTex: texture_2d<f32>;
@group(0) @binding(2) var srcSamp: sampler;
@group(0) @binding(3) var dstTex: texture_2d<f32>;
@group(0) @binding(4) var dstSamp: sampler;
@group(0) @binding(5) var skyTex: texture_cube<f32>;
@group(0) @binding(6) var skySamp: sampler;
struct Out { @builtin(position) pos: vec4f, @location(0) uv: vec2f, }
@vertex fn vs(@builtin(vertex_index) i: u32) -> Out {
  let p = array<vec2f, 3>(vec2f(-1., -1.), vec2f(3., -1.), vec2f(-1., 3.));
  var o: Out; o.pos = vec4f(p[i], 0., 1.); o.uv = p[i] * .5 + .5; return o;
}
fn sampled(tex: texture_2d<f32>, samp: sampler, n: vec3f, kind: f32, tilt: f32, offset: f32, feather: f32) -> vec3f {
  // Identical to Globe3D surface UV: north remains at the texture's top
  // throughout the exterior morph instead of flipping between renderers.
  let tilted = vec3f(n.x, cos(tilt) * n.y - sin(tilt) * n.z,
    sin(tilt) * n.y + cos(tilt) * n.z);
  let longitude = atan2(-tilted.x, -tilted.z);
  let latitude = asin(clamp(tilted.y, -1., 1.));
  let uv = vec2f(fract(longitude / 6.2831853 + .5 + offset), .5 - latitude / 3.1415926);
  let side = min(uv.x, 1.0 - uv.x);
  let inset = max(side, 0.002);
  let innerU = select(inset, 1.0 - inset, uv.x > 0.5);
  let oppositeU = 1.0 - innerU;
  let inner = textureSampleLevel(tex, samp, vec2f(innerU, uv.y), 0.).rgb;
  let opposite = textureSampleLevel(tex, samp, vec2f(oppositeU, uv.y), 0.).rgb;
  let mapped = select(textureSampleLevel(tex, samp, uv, 0.).rgb,
    mix((inner + opposite) * .5, inner, smoothstep(0.002, max(feather, 0.0021), side)),
    feather > 0.0 && side < feather);
  // Linear albedo, as Globe3D's outerAlbedo — including its fallbacks, so a
  // missing image looks the same on both sides of the handoff.
  if (kind < .5) { return mapped; }
  return fallbackAlbedo(n, tilted, kind);
}
struct Ray { O: vec3f, D: vec3f, }
fn ray(uv: vec2f) -> Ray {
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
  let ndc = uv * 2. - 1.;
  let fov = tan(.5 * .80285146); // Globe3D HERO.fov = 46 degrees
  return Ray(eyeDir * u.eyeDist, normalize(forward + ndc.x * fov * right + ndc.y * fov * up));
}
/** One ray's display colour: the blended planet, or the sky past it. \`pix\` is the pixel's own, borrowed by an edge pixel's second ray. */
fn shade(r: Ray, radius: f32, atmo: Atmo, pix: f32) -> vec3f {
  let b = dot(r.O, r.D);
  let disc = b * b - (dot(r.O, r.O) - radius * radius);
  if (disc <= 0.) {
    let sky = atmosphereHalo(skyRadiance(r.D, pix), r.O, r.D, radius, atmo);
    return toSrgb(shoulder(sky));
  }
  let n = normalize(r.O + (-b - sqrt(disc)) * r.D);
  let src = sampled(srcTex, srcSamp, n, u.srcKind, u.srcTilt, u.srcOffset, u.srcFeather);
  let dst = sampled(dstTex, dstSamp, n, u.dstKind, u.dstTilt, u.dstOffset, u.dstFeather);
  return toSrgb(shoulder(throughAtmosphere(sunlit(mix(src, dst, u.blend), n), n, r.D, atmo)));
}
@fragment fn fs(in: Out) -> @location(0) vec4f {
  // Each end's sun is that planet's own at the reader's time of day (its
  // tilt sets the spin axis), so the light lands where each globe has it.
  sunDir = normalize(mix(u.srcSun.xyz, u.dstSun.xyz, u.blend));
  let radius = mix(u.srcRadius, u.dstRadius, u.blend);
  let r = ray(in.uv);
  // Same sky and sun as Globe3D's background, so the handoff has no seam.
  let pix = length(fwidth(r.D));
  // Atmosphere morphs with the body, so Venus's haze thins into Mars's
  // rather than switching at either end of the bridge.
  let a0 = mix(u.srcAtmosphere, u.dstAtmosphere, u.blend);
  let a1 = mix(u.srcTwilight, u.dstTwilight, u.blend);
  let a2 = mix(u.srcShape, u.dstShape, u.blend);
  let atmo = Atmo(a0.rgb, a0.w, a1.rgb, a1.w, a2.x, a2.y);
  // The silhouette, anti-aliased as Globe3D's edges are (see its \`fs\`):
  // the limb is where the ray's closest approach equals the radius.
  let limb = radius - length(r.O - dot(r.O, r.D) * r.D);
  let grad = vec2f(dpdx(limb), dpdy(limb));
  var col = shade(r, radius, atmo, pix);
  let len = length(grad);
  let s = limb / max(len, 1e-12);
  if (len > 1e-12 && abs(s) < 0.5) {
    // uv y runs up, framebuffer y down.
    let across = in.uv - (s + sign(s) * 0.75) * grad / len * vec2f(u.pixelUv, -u.pixelUv);
    col = mix(shade(ray(across), radius, atmo, pix), col, 0.5 + abs(s));
  }
  return vec4f(col, 1.);
}`;

export class PlanetMorphScene {
  private readonly uniform: GPUBuffer;
  private readonly pipeline: GPURenderPipeline;
  private readonly bind: GPUBindGroup;
  private readonly data = new Float32Array(48);
  private blend = 0;
  private readonly tilts: [number, number];
  constructor(private readonly device: GPUDevice, format: GPUTextureFormat,
    sourceTexture: SurfaceTexture, destinationTexture: SurfaceTexture,
    sourcePlanet: PlanetDefinition, destinationPlanet: PlanetDefinition,
    orientation: readonly [number, number, number] = [.55, .32, 6.5],
    sunHours = 0) {
    this.uniform = device.createBuffer({ size: this.data.byteLength, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    const module = device.createShaderModule({ code: source });
    this.pipeline = device.createRenderPipeline({ layout: "auto", vertex: { module, entryPoint: "vs" }, fragment: { module, entryPoint: "fs", targets: [{ format }] } });
    const sky = milkyWayCube(device);
    this.bind = device.createBindGroup({ layout: this.pipeline.getBindGroupLayout(0), entries: [
      { binding: 0, resource: { buffer: this.uniform } },
      { binding: 1, resource: sourceTexture.view }, { binding: 2, resource: sourceTexture.sampler },
      { binding: 3, resource: destinationTexture.view }, { binding: 4, resource: destinationTexture.sampler },
      { binding: 5, resource: sky.view }, { binding: 6, resource: sky.sampler },
    ] });
    const radius = (p: PlanetDefinition) => radiiFor(p).ro;
    this.data.set([0, radius(sourcePlanet), radius(destinationPlanet),
      // `kind` selects a missing-texture fallback, never a loaded planet.
      // The photographic path is always sampled and blended directly.
      sourceTexture.available ? 0 : fallbackKind(sourcePlanet),
      destinationTexture.available ? 0 : fallbackKind(destinationPlanet),
      orientation[0], orientation[1], orientation[2],
      sourcePlanet.visual.axialTiltDeg * Math.PI / 180,
      destinationPlanet.visual.axialTiltDeg * Math.PI / 180,
      // Keep a material's selected map orientation fixed across the handoff,
      // including Mars's source-join placement on the far side of the globe.
      SURFACE_MATERIALS[sourcePlanet.visual.surface].longitudeOffsetTurns,
      SURFACE_MATERIALS[destinationPlanet.visual.surface].longitudeOffsetTurns,
      SURFACE_MATERIALS[sourcePlanet.visual.surface].seamFeatherTurns,
      SURFACE_MATERIALS[destinationPlanet.visual.surface].seamFeatherTurns, 0, 0]);
    const [srcAir, dstAir] = [atmosphereOf(sourcePlanet), atmosphereOf(destinationPlanet)];
    this.data.set(srcAir.flatMap((v, i) => [...v, ...dstAir[i]]), 16);
    this.tilts = [sourcePlanet.visual.axialTiltDeg, destinationPlanet.visual.axialTiltDeg];
    this.setSunHours(sunHours);
  }
  /** Follow the time-of-day control while the bridge is on screen — see `Globe3D.setLighting`. */
  setSunHours(hours: number): void {
    this.data.set([...sunDirection(this.tilts[0], hours), 0, ...sunDirection(this.tilts[1], hours), 0], 40);
  }
  setBlend(value: number): void { this.blend = ease(Math.min(1, Math.max(0, value))); }
  /** Draw into `target`, whose size sets the silhouette's anti-aliasing step. */
  draw(target: GPUTexture): void {
    this.data[0] = this.blend;
    this.data[14] = 1 / target.width;
    this.device.queue.writeBuffer(this.uniform, 0, this.data);
    const enc = this.device.createCommandEncoder();
    const pass = enc.beginRenderPass({ colorAttachments: [{ view: target.createView(), clearValue: { r: .008, g: .009, b: .025, a: 1 }, loadOp: "clear", storeOp: "store" }] });
    pass.setPipeline(this.pipeline); pass.setBindGroup(0, this.bind); pass.draw(3); pass.end();
    this.device.queue.submit([enc.finish()]);
  }
  destroy(): void { this.uniform.destroy(); }
}
