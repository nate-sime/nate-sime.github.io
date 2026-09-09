/**
 * Exterior surface assets for the cosmetic planet renderers.
 *
 * A surface image is deliberately separate from a solver: it can fail to load,
 * be replaced, or be absent altogether without changing a single numerical
 * parameter. Every material therefore has a bindable 1×1 fallback and a tint
 * for the procedural shader path.
 */

export type SurfaceMaterialId = "earth-daymap" | "venus-procedural";

export type ProceduralSurface = "earth" | "venus";

export interface SurfaceMaterial {
  readonly id: SurfaceMaterialId;
  /** Optional because a deliberately procedural material is still complete. */
  readonly imageUrl?: string;
  /** Multiplies both the image and procedural fallback in the cutaway shader. */
  readonly tint: readonly [number, number, number];
  readonly procedural: ProceduralSurface;
  readonly attribution: string;
}

export const SURFACE_MATERIALS: Record<SurfaceMaterialId, SurfaceMaterial> = {
  "earth-daymap": {
    id: "earth-daymap",
    imageUrl: `${import.meta.env.BASE_URL}earth-daymap.jpg`,
    tint: [1, 1, 1],
    procedural: "earth",
    attribution: "NASA Visible Earth, Blue Marble (2002), public domain",
  },
  "venus-procedural": {
    id: "venus-procedural",
    tint: [1, 1, 1],
    procedural: "venus",
    attribution: "Procedural Venus surface; no surface imagery is presented as data",
  },
};

export interface SurfaceTexture {
  readonly view: GPUTextureView;
  readonly sampler: GPUSampler;
  readonly available: boolean;
}

/** Start a fetch/decode early; `null` on any failure. No GPU device is needed. */
export async function fetchSurfaceImage(material: SurfaceMaterial): Promise<ImageBitmap | null> {
  if (!material.imageUrl) return null;
  try {
    const res = await fetch(material.imageUrl);
    if (!res.ok) return null;
    return await createImageBitmap(await res.blob());
  } catch {
    return null;
  }
}

/** Turn a decoded bitmap (or its absence) into the resource every cutaway binds. */
export function toSurfaceTexture(device: GPUDevice, bitmap: ImageBitmap | null): SurfaceTexture {
  const sampler = device.createSampler({
    addressModeU: "repeat", addressModeV: "clamp-to-edge",
    magFilter: "linear", minFilter: "linear",
  });
  if (bitmap) {
    const texture = device.createTexture({
      size: [bitmap.width, bitmap.height], format: "rgba8unorm",
      usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
    });
    device.queue.copyExternalImageToTexture({ source: bitmap }, { texture }, [bitmap.width, bitmap.height]);
    bitmap.close();
    return { view: texture.createView(), sampler, available: true };
  }
  const dummy = device.createTexture({
    size: [1, 1], format: "rgba8unorm",
    usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
  });
  device.queue.writeTexture({ texture: dummy }, new Uint8Array([0, 0, 0, 255]), { bytesPerRow: 4 }, [1, 1]);
  return { view: dummy.createView(), sampler, available: false };
}
