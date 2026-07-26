/** Spline space: velocity recovery order and exact divergence-free property. */

import { describe, it, expect } from "vitest";
import { verify } from "../src/verify";

const rows = verify();
const last = rows[rows.length - 1];

describe("stream-function spline space", () => {
  it("recovers u = ∇×ψ at order p = 3", () => {
    // ψ interpolates at O(h^{p+1}); the curl costs one derivative.
    expect(Number(last.order)).toBeGreaterThan(2.8);
    expect(Number(last.order)).toBeLessThan(3.3);
  });

  it("is divergence-free to round-off at every resolution", () => {
    // Structural (div curl ≡ 0), so mesh-independent — not an asymptotic claim.
    for (const r of rows) expect(Number(r.maxDiv)).toBeLessThan(1e-12);
  });
});
