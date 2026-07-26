/** 5-point Gauss–Legendre: exact to degree 9. Integrands here are degree-6
 *  spline products times smooth rational weights in r, so this is ample. */

const X = [0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
const W = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];

export interface QPoint { x: number; w: number; }

export function gauss(a: number, b: number): QPoint[] {
  const c = (a + b) / 2, h = (b - a) / 2;
  return X.map((xi, i) => ({ x: c + h * xi, w: h * W[i] }));
}
