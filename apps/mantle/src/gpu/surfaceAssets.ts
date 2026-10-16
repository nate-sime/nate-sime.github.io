// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Exterior surface assets for the cosmetic planet renderers.
 *
 * A surface image is deliberately separate from a solver: it can fail to load,
 * be replaced, or be absent altogether without changing a single numerical
 * parameter. Every material therefore has a bindable 1×1 fallback and a tint
 * for the procedural shader path.
 */

export type SurfaceMaterialId = "earth-daymap" | "venus-procedural" | "venus-magellan" | "mars-viking";

export type ProceduralSurface = "earth" | "venus" | "mars";

export interface SurfaceMaterial {
  readonly id: SurfaceMaterialId;
  /** Optional because a deliberately procedural material is still complete. */
  readonly imageUrl?: string;
  /** Multiplies both the image and procedural fallback in the cutaway shader. */
  readonly tint: readonly [number, number, number];
  /** Rotates an equirectangular image so its source join stays off the hero view. */
  readonly longitudeOffsetTurns: number;
  /** Width of a symmetric blend across an imperfect left/right image join. */
  readonly seamFeatherTurns: number;
  /**
   * Grayscale artificial-light map on `imageUrl`'s own equirectangular grid,
   * glowing on the night side (`cityLights` in `wgsl.ts`). Like the day image
   * it may fail to load, and then the night side is simply dark.
   */
  readonly nightLightsUrl?: string;
  readonly procedural: ProceduralSurface;
  readonly attribution: string;
}

export const SURFACE_MATERIALS: Record<SurfaceMaterialId, SurfaceMaterial> = {
  "earth-daymap": {
    id: "earth-daymap",
    imageUrl: `${import.meta.env.BASE_URL}earth-daymap.jpg`,
    tint: [1, 1, 1],
    longitudeOffsetTurns: 0,
    seamFeatherTurns: 0,
    // NASA Black Marble 2016 (0.1°, 3600×1800), lights only: its moonlit
    // land and ocean are blue-cast, city light warm or white, so each pixel
    // keeps how far mean(r, g) exceeds b / 1.65, rescaled to reach 1 at a
    // white core, less a 3% floor; then resampled to the day map's 2048×1024.
    nightLightsUrl: `${import.meta.env.BASE_URL}earth-night.webp`,
    procedural: "earth",
    attribution: "NASA Visible Earth, Blue Marble (2002), and NASA Earth Observatory, Black Marble (2016), public domain",
  },
  "venus-procedural": {
    id: "venus-procedural",
    tint: [1, 1, 1],
    longitudeOffsetTurns: 0,
    seamFeatherTurns: 0,
    procedural: "venus",
    attribution: "Procedural Venus surface; no surface imagery is presented as data",
  },
  "venus-magellan": {
    id: "venus-magellan",
    imageUrl: `${import.meta.env.BASE_URL}venus-magellan.webp`,
    tint: [1, 1, 1],
    longitudeOffsetTurns: 0,
    seamFeatherTurns: 0,
    procedural: "venus",
    attribution: "NASA/JPL-Caltech Magellan radar-derived Venus texture",
  },
  "mars-viking": {
    id: "mars-viking",
    imageUrl: `${import.meta.env.BASE_URL}mars-viking.webp`,
    tint: [1, 1, 1],
    // This map's left/right source join has a visible tonal mismatch. Place
    // it on the far side of the default 3-D view and feather it continuously.
    longitudeOffsetTurns: 0.5,
    seamFeatherTurns: 0.035,
    procedural: "mars",
    attribution: "NASA/JPL-Caltech Viking imagery processed at USGS",
  },
};

export interface SurfaceTexture {
  readonly view: GPUTextureView;
  readonly sampler: GPUSampler;
  readonly available: boolean;
  /** The material's night lights, sampled with `sampler`; a black 1×1 when it has none or they failed to load. */
  readonly nightView: GPUTextureView;
  readonly nightAvailable: boolean;
}

/** A material's decoded images — either may be missing. */
export interface SurfaceImages {
  readonly day: ImageBitmap | null;
  readonly night: ImageBitmap | null;
}

const fetchImage = async (url: string | undefined): Promise<ImageBitmap | null> => {
  if (!url) return null;
  try {
    const res = await fetch(url);
    if (!res.ok) return null;
    return await createImageBitmap(await res.blob());
  } catch {
    return null;
  }
};

/** Start a fetch/decode early; each image is `null` on any failure. No GPU device is needed. */
export async function fetchSurfaceImages(material: SurfaceMaterial): Promise<SurfaceImages> {
  const [day, night] = await Promise.all([fetchImage(material.imageUrl), fetchImage(material.nightLightsUrl)]);
  return { day, night };
}

const MIP_BLIT = /* wgsl */ `
@group(0) @binding(0) var src: texture_2d<f32>;
@group(0) @binding(1) var samp: sampler;
struct Out { @builtin(position) pos: vec4f, @location(0) uv: vec2f, }
@vertex fn vs(@builtin(vertex_index) i: u32) -> Out {
  let p = array<vec2f, 3>(vec2f(-1., -1.), vec2f(3., -1.), vec2f(-1., 3.));
  var o: Out; o.pos = vec4f(p[i], 0., 1.); o.uv = vec2f(p[i].x * .5 + .5, .5 - p[i].y * .5); return o;
}
@fragment fn fs(in: Out) -> @location(0) vec4f { return textureSampleLevel(src, samp, in.uv, 0.); }`;

const mipBlits = new WeakMap<GPUDevice, { pipeline: GPURenderPipeline; sampler: GPUSampler }>();

/**
 * Fill every mip level of an `rgba8unorm-srgb` texture from level 0. Each
 * texel is one bilinear tap at the shared corner of the four above it — their
 * exact box average, for these power-of-two maps — and through sRGB views, so
 * it averages light rather than its encoding: a gamma-space average would dim
 * sparse bright points, like city lights, as they shrink.
 */
export function generateMips(device: GPUDevice, texture: GPUTexture): void {
  let blit = mipBlits.get(device);
  if (!blit) {
    const module = device.createShaderModule({ code: MIP_BLIT });
    blit = {
      pipeline: device.createRenderPipeline({
        layout: "auto",
        vertex: { module, entryPoint: "vs" },
        fragment: { module, entryPoint: "fs", targets: [{ format: "rgba8unorm-srgb" }] },
      }),
      sampler: device.createSampler({ magFilter: "linear", minFilter: "linear" }),
    };
    mipBlits.set(device, blit);
  }
  const enc = device.createCommandEncoder();
  for (let level = 1; level < texture.mipLevelCount; level++) {
    const pass = enc.beginRenderPass({ colorAttachments: [{
      view: texture.createView({ baseMipLevel: level, mipLevelCount: 1 }),
      loadOp: "clear", storeOp: "store", clearValue: { r: 0, g: 0, b: 0, a: 1 },
    }] });
    pass.setPipeline(blit.pipeline);
    pass.setBindGroup(0, device.createBindGroup({
      layout: blit.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: texture.createView({ baseMipLevel: level - 1, mipLevelCount: 1 }) },
        { binding: 1, resource: blit.sampler },
      ],
    }));
    pass.draw(3);
    pass.end();
  }
  device.queue.submit([enc.finish()]);
}

/**
 * sRGB-tagged, so sampling returns linear colour: the exterior shaders light
 * in linear space and re-encode (`SPACE` in `wgsl.ts`). The night map is
 * encoded the same way, so its lights decode to linear radiance too, and is
 * the one map given mips: they pick its level for the pixel (`equirectLod`)
 * and supply its glow (`cityLights`). The day images are still sampled at
 * full resolution only.
 */
const upload = (device: GPUDevice, bitmap: ImageBitmap, mips = false): GPUTextureView => {
  const mipLevelCount = mips ? Math.floor(Math.log2(Math.max(bitmap.width, bitmap.height))) + 1 : 1;
  const texture = device.createTexture({
    size: [bitmap.width, bitmap.height], format: "rgba8unorm-srgb", mipLevelCount,
    usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
  });
  device.queue.copyExternalImageToTexture({ source: bitmap }, { texture }, [bitmap.width, bitmap.height]);
  bitmap.close();
  if (mips) generateMips(device, texture);
  return texture.createView();
};

const blank = (device: GPUDevice): GPUTextureView => {
  const dummy = device.createTexture({
    size: [1, 1], format: "rgba8unorm",
    usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
  });
  device.queue.writeTexture({ texture: dummy }, new Uint8Array([0, 0, 0, 255]), { bytesPerRow: 4 }, [1, 1]);
  return dummy.createView();
};

/** Turn decoded images (or their absence) into the resource every cutaway binds. */
export function toSurfaceTexture(device: GPUDevice, images: SurfaceImages | null): SurfaceTexture {
  const sampler = device.createSampler({
    addressModeU: "repeat", addressModeV: "clamp-to-edge",
    magFilter: "linear", minFilter: "linear", mipmapFilter: "linear",
  });
  const day = images?.day ?? null, night = images?.night ?? null;
  return {
    view: day ? upload(device, day) : blank(device), sampler, available: day !== null,
    nightView: night ? upload(device, night, true) : blank(device), nightAvailable: night !== null,
  };
}
