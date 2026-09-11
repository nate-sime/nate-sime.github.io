import { afterEach, describe, expect, it } from "vitest";
import { gpuDevice, gpuErrors } from "./gpu";
import { PlanetMorphScene } from "../src/gpu/planetMorph";
import { toSurfaceTexture } from "../src/gpu/surfaceAssets";
import { EARTH, MARS, VENUS } from "../src/planets";

const device = await gpuDevice();
afterEach(() => expect(gpuErrors.splice(0).join(" | ")).toBe(""));

(device ? describe : describe.skip)("planet exterior morph", () => {
  it("draws both endpoints and an interpolated geometry/material frame", async () => {
    const fallback = toSurfaceTexture(device!, null);
    const scene = new PlanetMorphScene(device!, "rgba8unorm", fallback, fallback, EARTH, VENUS);
    const tex = device!.createTexture({ size: [64, 64], format: "rgba8unorm", usage: GPUTextureUsage.RENDER_ATTACHMENT });
    for (const blend of [0, .5, 1]) { scene.setBlend(blend); scene.draw(tex.createView()); }
    const mars = new PlanetMorphScene(device!, "rgba8unorm", fallback, fallback, VENUS, MARS);
    for (const blend of [0, .5, 1]) { mars.setBlend(blend); mars.draw(tex.createView()); }
    await device!.queue.onSubmittedWorkDone();
    tex.destroy(); scene.destroy(); mars.destroy();
  });
});
