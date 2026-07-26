/** Constant-μ Stokes: convergence, free-slip BC, symbols, k=0 kernel. */

import { describe, it, expect } from "vitest";
import { clampedAxis, periodicAxis, Field } from "../src/spline";
import { StokesSolver, loadVector } from "../src/solver/stokes";
import { azimuthalSymbols } from "../src/solver/operators";
import { stokesConvergence, source, Ri, Ro, K } from "../src/mms";

describe("biharmonic Stokes solve", () => {
  it("converges at the optimal rate h^(p+1) = h⁴", () => {
    const rows = stokesConvergence([12, 24, 48]);
    expect(Number(rows[rows.length - 1].order)).toBeGreaterThan(3.8);
    expect(Number(rows[rows.length - 1].err)).toBeLessThan(1e-7);
  });

  // The decisive check that free-slip — not ω = 0 — is what the
  // discrete problem actually imposes. Being *natural* to the dissipation form it
  // holds only in the limit, so the discriminating assertion is the **rate**:
  // ε_rφ involves ψ'', which converges at h^(p−1) = h². An ω = 0 implementation
  // would instead plateau at the nonzero limit 2u_φ/r, failing this outright.
  it("satisfies the curved-boundary free-slip condition ε_rφ → 0 at O(h²)", () => {
    const err = [16, 32, 64].map((nr) => {
      const rAx = clampedAxis(nr, Ri, Ro), aAx = periodicAxis(64);
      const psi = new StokesSolver(rAx, aAx).solve(loadVector(rAx, aAx, source));
      const f = new Field(rAx, aAx);
      for (let i = 0; i < nr; i++) f.c[i].set(psi[i]);

      let e = 0;
      for (let j = 0; j < 64; j++) {
        const phi = (2 * Math.PI * j) / 64;
        e = Math.max(e, Math.abs(f.shearRate(Ri, phi)), Math.abs(f.shearRate(Ro, phi)));
      }
      return e;
    });
    for (let i = 1; i < err.length; i++)
      expect(Math.log2(err[i - 1] / err[i])).toBeGreaterThan(1.8);
  });

  it("keeps the solution exactly divergence-free", () => {
    const rAx = clampedAxis(32, Ri, Ro), aAx = periodicAxis(64);
    const psi = new StokesSolver(rAx, aAx).solve(loadVector(rAx, aAx, source));
    const f = new Field(rAx, aAx);
    for (let i = 0; i < 32; i++) f.c[i].set(psi[i]);

    let d = 0;
    for (let i = 1; i < 40; i++)
      for (let j = 0; j < 60; j++)
        d = Math.max(d, Math.abs(f.divergence(
          Ri + ((Ro - Ri) * i) / 40, (2 * Math.PI * j) / 60)));
    expect(d).toBeLessThan(1e-12);
  });

  // The per-mode multiplier is the discrete B-spline symbol. It must
  // be *consistent* with k² for resolved modes — using k² itself would not be.
  it("has azimuthal symbols consistent with k² for resolved modes", () => {
    const [S0, S1] = azimuthalSymbols(periodicAxis(128));
    for (const k of [1, 2, 4, 8]) {
      expect(S1[k] / S0[k]).toBeCloseTo(k * k, 3);
      expect(S1[128 - k] / S0[128 - k]).toBeCloseTo(k * k, 3); // Hermitian pair
    }
    expect(S1[0]).toBeCloseTo(0, 12); // constants are in the kernel of ∫N'N'
  });

  it("leaves the k = 0 mode at rest (singular block, zero forcing)", () => {
    const rAx = clampedAxis(24, Ri, Ro), aAx = periodicAxis(64);
    const psi = new StokesSolver(rAx, aAx).solve(loadVector(rAx, aAx, source));
    // Source ∝ cos(Kφ) with K ≠ 0 ⇒ no azimuthal mean flow may be generated.
    for (let i = 0; i < 24; i++) {
      const mean = psi[i].reduce((a, b) => a + b, 0) / 64;
      expect(Math.abs(mean)).toBeLessThan(1e-12);
    }
    expect(K).not.toBe(0);
  });
});
