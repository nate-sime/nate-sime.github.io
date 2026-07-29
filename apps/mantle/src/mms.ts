/**
 * Manufactured solution for constant-μ Stokes.
 *
 * Take ψ* = F(s) cos(Kφ), s = (r−r_i)/L, with F = s³(1−s)³. The triple zeros at
 * s = 0, 1 give F = F' = F'' = 0 at both radii, so ψ* satisfies **both physical
 * boundary conditions homogeneously**: ψ* = 0 and ψ*_rr − ψ*_r/r = 0. This
 * matters — an MMS carrying inhomogeneous boundary data would not exercise the
 * free-slip condition at all.
 *
 * The load is the strong-form source S = μΔ²ψ*, using
 *
 *   Δ²[F(r)cos kφ] = [F'''' + 2F'''/r − (1+2k²)F''/r²
 *                     + (1+2k²)F'/r³ + (k⁴−4k²)F/r⁴] cos kφ,
 *
 * valid since ⟨S, v⟩ = ∫(∇×f)_z v dx = ∫ f·u[v] dx = a(ψ*, v) for v vanishing on
 * ∂Ω. Driving the dissipation form with the biharmonic source therefore tests
 * the derivation of a(·,·) against the PDE, not merely its own assembly.
 */

import { clampedAxis, periodicAxis, Field } from "./spline";
import { StokesSolver, loadVector } from "./solver/stokes";

export const Ri = 0.55, Ro = 1.0, K = 4, MU = 1;
const L = Ro - Ri;
const s = (r: number) => (r - Ri) / L;

// F = s³ − 3s⁴ + 3s⁵ − s⁶ and its s-derivatives.
const F0 = (t: number) => t ** 3 - 3 * t ** 4 + 3 * t ** 5 - t ** 6;
const F1 = (t: number) => 3 * t ** 2 - 12 * t ** 3 + 15 * t ** 4 - 6 * t ** 5;
const F2 = (t: number) => 6 * t - 36 * t ** 2 + 60 * t ** 3 - 30 * t ** 4;
const F3 = (t: number) => 6 - 72 * t + 180 * t ** 2 - 120 * t ** 3;
const F4 = (t: number) => -72 + 360 * t - 360 * t ** 2;

export const psiEx = (r: number, phi: number) => F0(s(r)) * Math.cos(K * phi);

export function source(r: number, phi: number): number {
  const t = s(r);
  const f0 = F0(t), f1 = F1(t) / L, f2 = F2(t) / L ** 2,
    f3 = F3(t) / L ** 3, f4 = F4(t) / L ** 4;
  const k2 = K * K, k4 = k2 * k2;
  const b = f4 + (2 * f3) / r - ((1 + 2 * k2) * f2) / r ** 2
    + ((1 + 2 * k2) * f1) / r ** 3 + ((k4 - 4 * k2) * f0) / r ** 4;
  return MU * b * Math.cos(K * phi);
}

export interface Row { nr: number; na: number; h: string; err: string; order: string; }

export function stokesConvergence(levels = [12, 24, 48, 96]): Row[] {
  const rows: Row[] = [];
  let prev = 0;
  for (const nr of levels) {
    const na = 32 * (nr / levels[0]);
    const rAx = clampedAxis(nr, Ri, Ro), aAx = periodicAxis(na);
    const psi = new StokesSolver(rAx, aAx, () => MU).solve(loadVector(rAx, aAx, source));

    const f = new Field(rAx, aAx);
    for (let i = 0; i < nr; i++) f.c[i].set(psi[i]);

    let err = 0;
    for (let i = 1; i < 60; i++) {
      const r = Ri + (L * i) / 60;
      for (let j = 0; j < 120; j++) {
        const phi = (2 * Math.PI * j) / 120;
        err = Math.max(err, Math.abs(f.eval(r, phi).psi - psiEx(r, phi)));
      }
    }
    const order = prev ? Math.log2(prev / err) : NaN;
    prev = err;
    rows.push({
      nr, na,
      h: (L / (nr - 3)).toExponential(2),
      err: err.toExponential(2),
      order: isNaN(order) ? "  — " : order.toFixed(2),
    });
  }
  return rows;
}
