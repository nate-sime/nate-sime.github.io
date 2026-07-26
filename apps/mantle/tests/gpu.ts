/**
 * Headless WebGPU for the test suite (Dawn, via the `webgpu` package).
 *
 * `npm test` must stay runnable on a machine — or a CI runner — with no usable
 * adapter, so this resolves to `null` rather than throwing and the GPU suites
 * skip themselves. Everything the GPU tests assert is a *parity* claim against
 * the CPU reference, which is exercised unconditionally elsewhere.
 */

let cached: Promise<GPUDevice | null> | undefined;

/**
 * Keeps the `GPU` and `GPUAdapter` wrappers reachable for the process lifetime.
 * Dawn's node bindings do not root them from the device, so once they become
 * unreachable the finaliser frees native state the live device still points at
 * and the next GPU call segfaults — reliably, but only after an idle period long
 * enough for a GC, which makes it look like a timing bug in the code under test.
 */
const roots: unknown[] = [];

/** Validation faults seen since the last check. */
export const gpuErrors: string[] = [];

export function gpuDevice(): Promise<GPUDevice | null> {
  cached ??= (async () => {
    try {
      const { create, globals } = await import("webgpu");
      Object.assign(globalThis, globals);
      const gpu = create([]);
      const adapter = await gpu.requestAdapter();
      const device = (await adapter?.requestDevice()) ?? null;
      roots.push(gpu, adapter, device);
      // Collect rather than throw: this fires inside a native callback, where any
      // escaping exception is a napi fatal error that takes the worker process
      // down and hides the message that caused it.
      device?.addEventListener("uncapturederror", (e) => {
        try {
          gpuErrors.push(String((e as GPUUncapturedErrorEvent)?.error?.message ?? e));
        } catch {
          gpuErrors.push("uncaptured GPU error");
        }
      });
      return device;
    } catch {
      return null;
    }
  })();
  return cached;
}

/** Max |a − b| over two flat fields. */
export function maxDiff(a: ArrayLike<number>, b: ArrayLike<number>): number {
  let m = 0;
  for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
}
