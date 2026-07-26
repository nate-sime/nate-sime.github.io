/**
 * Spline-space verification. A manufactured stream function
 *
 *   ψ(r, φ) = sin(π s) cos(m φ),   s = (r − r_i)/(r_o − r_i),
 *
 * is interpolated into the spline space; the recovered velocity u = ∇×ψ is
 * compared to its closed form under mesh refinement. Expected behaviour:
 *   • ‖u_h − u‖_∞ = O(h^p)  (ψ interpolates at O(h^{p+1}), one derivative lost),
 *   • ‖∇·u_h‖_∞ ≈ machine ε (exact divergence-free, structural).
 */

import { clampedAxis, periodicAxis, Field } from "./spline";

export const Ri = 0.55, Ro = 1.0; // Earth-like radius ratio
const M = 4, L = Ro - Ri;
const s = (r: number) => (r - Ri) / L;

export const psiEx = (r: number, phi: number) => Math.sin(Math.PI * s(r)) * Math.cos(M * phi);
const urEx = (r: number, phi: number) => -(M / r) * Math.sin(Math.PI * s(r)) * Math.sin(M * phi);
const upEx = (r: number, phi: number) => -(Math.PI / L) * Math.cos(Math.PI * s(r)) * Math.cos(M * phi);

export interface Row { nr: number; na: number; h: string; errU: string; order: string; maxDiv: string; }

export function verify(): Row[] {
  const rows: Row[] = [];
  let prev = 0;
  for (const K of [1, 2, 4, 8]) {
    const nr = 8 * K, na = 12 * K;
    const f = new Field(clampedAxis(nr, Ri, Ro), periodicAxis(na));
    f.interpolate(psiEx);

    let errU = 0, maxDiv = 0;
    const NR = 121, NA = 240;
    for (let i = 1; i < NR; i++) {
      const r = Ri + (L * i) / NR;
      for (let j = 0; j < NA; j++) {
        const phi = (2 * Math.PI * j) / NA;
        const v = f.velocity(r, phi);
        errU = Math.max(errU, Math.abs(v.ur - urEx(r, phi)), Math.abs(v.up - upEx(r, phi)));
        maxDiv = Math.max(maxDiv, Math.abs(f.divergence(r, phi)));
      }
    }
    const h = L / (nr - 3);
    const order = prev ? Math.log2(prev / errU) : NaN;
    prev = errU;
    rows.push({
      nr, na,
      h: h.toExponential(2),
      errU: errU.toExponential(2),
      order: isNaN(order) ? "  — " : order.toFixed(2),
      maxDiv: maxDiv.toExponential(1),
    });
  }
  return rows;
}
