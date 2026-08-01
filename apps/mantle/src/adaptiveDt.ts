/**
 * Host-lagged adaptive dt: the CFL-implied step, gated by a hysteresis band
 * against the cost of `GpuSimulation.setDt` — a ~60k f64 flop refactorisation
 * of the diffusion operator, not a uniform write.
 *
 * The lag is deliberate, not a rounding compromise. `maxSpeed` is read back
 * asynchronously through `pollStats`, so it is always a frame or more stale by
 * the time it reaches here; Stokes is quasi-static and max|u| evolves on the
 * convective timescale, so that staleness is far inside the 10–20% margin a
 * Courant number already carries. Being briefly generous costs accuracy, not
 * stability — dt here is an accuracy knob, never a stability limit (see
 * `gpu/sim.ts`).
 */
export function adaptiveDt(
  current: number, dtMax: number, courant: number, maxSpeed: number, band = 0.15,
): number | null {
  if (!(maxSpeed > 0)) return null;   // no reading yet, or a genuinely still flow
  const target = Math.min(dtMax, courant / maxSpeed);
  // A target pinned at the cap is not a CFL estimate to damp against jitter —
  // it is the user's own ceiling, set explicitly on the pane. Banding it the
  // same way as a moving target left it unreachable: once `current` drifted
  // within the band below `dtMax`, raising the Courant number past the point
  // where the cap binds had no further effect, silently.
  if (target === dtMax) return current === dtMax ? null : target;
  return Math.abs(target / current - 1) > band ? target : null;
}
