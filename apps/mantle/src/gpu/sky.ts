// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The Milky Way, rendered once into a cube map rather than once per pixel per
 * frame (`MILKY_WAY` in `wgsl.ts` has the why). It is fixed to the world, so
 * nothing ever invalidates it: one texture per device, shared by every
 * `Globe3D` a rebuild constructs and by the planet bridge, and never freed —
 * like the surface images, it lives as long as the app does.
 *
 * Only the smooth band is baked. Stars stay live (`starLayer`): baked, they
 * would blur or swell wherever a cube texel is larger than a screen pixel.
 * The sun moves with the time of day, so it stays live too.
 */

import { MILKY_WAY, NOISE3 } from "./wgsl";

/**
 * Texels along a face edge — about 0.006 rad each, several to the period of
 * the band's finest octave that still has any amplitude; linear filtering
 * covers the rest. 6 × 256² × 8 bytes is 3 MB.
 */
export const SKY_SIZE = 256;

export interface SkyTexture {
  readonly view: GPUTextureView;
  readonly sampler: GPUSampler;
}

/**
 * The world direction through the centre of texel (x, y) on cube face
 * `face` — WebGPU's face order (+X, −X, +Y, −Y, +Z, −Z) and orientation, row
 * 0 at the top. Exported for the test that checks a sample against
 * `milkyWay` evaluated directly.
 */
export const CUBE_DIRECTION = /* wgsl */ `
fn cubeDirection(face: u32, texel: vec2f, size: f32) -> vec3f {
  let s = texel.x / size * 2.0 - 1.0;
  let t = texel.y / size * 2.0 - 1.0;
  switch (face) {
    case 0u: { return normalize(vec3f(1.0, -t, -s)); }
    case 1u: { return normalize(vec3f(-1.0, -t, s)); }
    case 2u: { return normalize(vec3f(s, 1.0, t)); }
    case 3u: { return normalize(vec3f(s, -1.0, -t)); }
    case 4u: { return normalize(vec3f(s, -t, 1.0)); }
    default: { return normalize(vec3f(-s, -t, -1.0)); }
  }
}
`;

// One fullscreen triangle per face; the face index rides in as the instance
// index, so the six passes need no uniform buffer between them.
const bakeSource = NOISE3 + MILKY_WAY + CUBE_DIRECTION + /* wgsl */ `
struct Out { @builtin(position) pos: vec4f, @location(0) @interpolate(flat) face: u32, }
@vertex fn vs(@builtin(vertex_index) i: u32, @builtin(instance_index) face: u32) -> Out {
  let p = array<vec2f, 3>(vec2f(-1., -1.), vec2f(3., -1.), vec2f(-1., 3.));
  var o: Out; o.pos = vec4f(p[i], 0., 1.); o.face = face; return o;
}
@fragment fn fs(in: Out) -> @location(0) vec4f {
  return vec4f(milkyWay(cubeDirection(in.face, in.pos.xy, ${SKY_SIZE.toFixed(1)})), 1.0);
}`;

const baked = new WeakMap<GPUDevice, SkyTexture>();

/** The baked Milky Way for `device`, rendering it on first use. */
export function milkyWayCube(device: GPUDevice): SkyTexture {
  const cached = baked.get(device);
  if (cached) return cached;
  // Half-float: the band peaks near 0.01 in linear light, too faint for
  // eight bits to hold its gradations without banding.
  const format: GPUTextureFormat = "rgba16float";
  const texture = device.createTexture({
    size: [SKY_SIZE, SKY_SIZE, 6], format,
    usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING,
  });
  const module = device.createShaderModule({ code: bakeSource });
  const pipeline = device.createRenderPipeline({
    layout: "auto",
    vertex: { module, entryPoint: "vs" },
    fragment: { module, entryPoint: "fs", targets: [{ format }] },
  });
  const enc = device.createCommandEncoder();
  for (let face = 0; face < 6; face++) {
    const pass = enc.beginRenderPass({ colorAttachments: [{
      view: texture.createView({ dimension: "2d", baseArrayLayer: face, arrayLayerCount: 1 }),
      loadOp: "clear", storeOp: "store", clearValue: { r: 0, g: 0, b: 0, a: 1 },
    }] });
    pass.setPipeline(pipeline);
    pass.draw(3, 1, 0, face);
    pass.end();
  }
  device.queue.submit([enc.finish()]);
  const sky: SkyTexture = {
    view: texture.createView({ dimension: "cube" }),
    sampler: device.createSampler({ magFilter: "linear", minFilter: "linear" }),
  };
  baked.set(device, sky);
  return sky;
}
