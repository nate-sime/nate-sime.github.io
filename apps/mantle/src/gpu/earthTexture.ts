/**
 * The globe view's exterior shell, textured with a real photographic map of
 * Earth rather than the procedural fbm "continents" `globe.ts` falls back to
 * (`shadeOuter` in `wgsl.ts`, gated on `Globe.hasEarth`).
 *
 * `public/earth-daymap.jpg` is NASA's "Blue Marble: Land Surface, Shallow
 * Water, and Shaded Topography" (Visible Earth, 11 Feb 2002; MODIS/Terra
 * observations, Jun–Sep 2001), a US-government work and so public domain —
 * see the on-screen credit `main.ts` shows whenever the globe is textured
 * with it, which is the citation this asset needs, not a separate licence
 * file. https://visibleearth.nasa.gov/images/57752
 *
 * Fetching and decoding is split from GPU upload (`fetchEarthImage` takes no
 * device) so `main.ts` can start the network request before the adapter/
 * device negotiation that has to happen anyway, rather than paying for both
 * in sequence. Either step failing — offline, blocked, a decode error — is
 * caught here and turned into `available: false` plus a harmless 1×1 stand-in
 * texture, never a thrown error: the caller always gets something bindable,
 * and the shader reads the flag to fall back to the procedural shell.
 */

export interface EarthTexture {
  readonly view: GPUTextureView;
  readonly sampler: GPUSampler;
  readonly available: boolean;
}

const EARTH_IMAGE_URL = `${import.meta.env.BASE_URL}earth-daymap.jpg`;

/** Start the fetch/decode early; `null` on any failure. No GPU device needed yet. */
export async function fetchEarthImage(): Promise<ImageBitmap | null> {
  try {
    const res = await fetch(EARTH_IMAGE_URL);
    if (!res.ok) return null;
    return await createImageBitmap(await res.blob());
  } catch {
    return null;
  }
}

/** Turn a decoded bitmap (or its absence) into the bindable resource `Globe3D` samples. */
export function toEarthTexture(device: GPUDevice, bitmap: ImageBitmap | null): EarthTexture {
  // Longitude wraps (`repeat`); latitude does not (`clamp-to-edge` short of
  // the poles reads as a solid cap rather than smearing texel colour across
  // the seam a `repeat` wrap would draw there instead).
  const sampler = device.createSampler({
    addressModeU: "repeat", addressModeV: "clamp-to-edge",
    magFilter: "linear", minFilter: "linear",
  });
  if (bitmap) {
    const texture = device.createTexture({
      size: [bitmap.width, bitmap.height],
      format: "rgba8unorm",
      // RENDER_ATTACHMENT: copyExternalImageToTexture requires it even though
      // nothing here ever renders to this texture directly.
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
