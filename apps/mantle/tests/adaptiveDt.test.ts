/**
 * Pure-function tests for the host-lagged hysteresis band, decoupled from the
 * GPU: the invariant that `setDt` fires only outside the band belongs here,
 * not in a spy on a buffer write (see PLAN-adaptive-dt.md).
 */

import { describe, it, expect } from "vitest";
import { adaptiveDt } from "../src/adaptiveDt";

describe("adaptiveDt", () => {
  it("returns null before any reading has arrived, or when the flow is still", () => {
    expect(adaptiveDt(1e-4, 5e-4, 1, NaN)).toBeNull();
    expect(adaptiveDt(1e-4, 5e-4, 1, 0)).toBeNull();
  });

  it("caps the CFL-implied step at dtMax", () => {
    expect(adaptiveDt(1e-4, 5e-4, 1, 1)).toBe(5e-4); // implied 1.0 far exceeds the cap
  });

  it("stays put inside the hysteresis band", () => {
    // implied dt = 1 / 9000 ≈ 1.111e-4 — 11% above current, inside the 15% band
    expect(adaptiveDt(1e-4, 5e-4, 1, 9000)).toBeNull();
  });

  it("moves once the implied step drifts past the band", () => {
    // implied dt = 1 / 20000 = 5e-5 — 50% below current, well outside the band
    expect(adaptiveDt(1e-4, 5e-4, 1, 20000)).toBeCloseTo(5e-5);
  });

  it("honours a narrower band when asked", () => {
    // implied dt ≈ 1.111e-4 is 11% away — outside a 5% band, inside the default 15%
    expect(adaptiveDt(1e-4, 5e-4, 1, 9000, 0.05)).toBeCloseTo(1 / 9000);
    expect(adaptiveDt(1e-4, 5e-4, 1, 9000, 0.15)).toBeNull();
  });

  it("always reaches dtMax once the CFL target exceeds it, even inside the band", () => {
    // current = 9.5e-5, dtMax = 1e-4 is only 5.3% away — inside the default 15%
    // band, so a plain hysteresis check would strand dt there forever, and no
    // amount of raising `courant` past this point could ever close the gap.
    const maxSpeed = 20000;   // courant / maxSpeed = 1.5e-4 > dtMax, so the cap binds
    expect(adaptiveDt(9.5e-5, 1e-4, 3, maxSpeed)).toBe(1e-4);
  });

  it("does not re-signal a change once dt has actually reached dtMax", () => {
    const maxSpeed = 20000;
    expect(adaptiveDt(1e-4, 1e-4, 3, maxSpeed)).toBeNull();
  });
});
