// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Just enough 3-D for the introduction (`views/intro.ts`): a camera that turns
 * about the vertical, mild perspective, and flat-shaded quads drawn back to
 * front into a 2-D canvas — the painter's algorithm. A beam is a few hundred
 * quads and a plate about a thousand, so sorting them every frame costs
 * nothing worth a WebGL context.
 *
 * World axes: x along the span, y across it, z up.
 */

export type Vec3 = readonly [number, number, number];
type RGB = readonly [number, number, number];

export interface Camera {
  /** Turn about z, then tilt toward the viewer about the screen's horizontal. */
  yaw: number;
  pitch: number;
}

export interface Quad {
  readonly p: readonly [Vec3, Vec3, Vec3, Vec3];
  readonly rgb: RGB;
  /** Drawn straight after the quad's fill, in screen space: what lies on it, so it is hidden with it. */
  readonly then?: (ctx: CanvasRenderingContext2D, to: (v: Vec3) => [number, number]) => void;
}

/** Light from over the viewer's left shoulder, in world axes. */
const LIGHT: Vec3 = (() => {
  const l = [-0.35, -0.55, 0.76], n = Math.hypot(...l);
  return l.map((v) => v / n) as unknown as Vec3;
})();

/** World to camera: x right, y up the screen, z away from the viewer. */
function view(c: Camera, [x, y, z]: Vec3): Vec3 {
  const cy = Math.cos(c.yaw), sy = Math.sin(c.yaw), cp = Math.cos(c.pitch), sp = Math.sin(c.pitch);
  const X = cy * x - sy * y, Y = sy * x + cy * y;
  // Looking down by `pitch`: what is farther back shows higher, what is higher is nearer.
  return [X, cp * z + sp * Y, cp * Y - sp * z];
}

export interface Projector {
  readonly to: (v: Vec3) => [number, number];
  readonly depth: (v: Vec3) => number;
}

/**
 * A projection that fits `bounds` (world-space corners of what will be drawn)
 * into a w × h box with a margin, centred, for this camera.
 */
export function projector(c: Camera, bounds: readonly Vec3[], w: number, h: number, margin = 14): Projector {
  const D = 5; // eye distance, in world units: mild perspective
  const persp = (v: Vec3): [number, number] => {
    const [x, y, z] = view(c, v), f = D / (D + z);
    return [x * f, y * f];
  };
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const b of bounds) {
    const [x, y] = persp(b);
    x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y);
  }
  const s = Math.min((w - 2 * margin) / (x1 - x0), (h - 2 * margin) / (y1 - y0));
  const ox = w / 2 - (s * (x0 + x1)) / 2, oy = h / 2 + (s * (y0 + y1)) / 2;
  return {
    to: (v) => {
      const [x, y] = persp(v);
      return [ox + s * x, oy - s * y];
    },
    depth: (v) => view(c, v)[2],
  };
}

/** Fill every quad, farthest first, each shaded by how squarely it faces the light. */
export function paintQuads(ctx: CanvasRenderingContext2D, pr: Projector, quads: readonly Quad[]): void {
  const order = quads.map((q, i) => ({ i, d: q.p.reduce((s, v) => s + pr.depth(v), 0) }));
  order.sort((a, b) => b.d - a.d);
  ctx.lineJoin = "round";
  for (const { i } of order) {
    const q = quads[i], [a, b, c, d] = q.p;
    // Normal from the diagonals; either side may face the light.
    const u = [c[0] - a[0], c[1] - a[1], c[2] - a[2]], v = [d[0] - b[0], d[1] - b[1], d[2] - b[2]];
    const n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    const len = Math.hypot(n[0], n[1], n[2]) || 1;
    const lit = 0.5 + 0.5 * Math.abs((n[0] * LIGHT[0] + n[1] * LIGHT[1] + n[2] * LIGHT[2]) / len);
    const fill = `rgb(${(q.rgb[0] * lit) | 0}, ${(q.rgb[1] * lit) | 0}, ${(q.rgb[2] * lit) | 0})`;
    ctx.beginPath();
    for (const [k, p] of q.p.entries()) {
      const [x, y] = pr.to(p);
      if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y);
    }
    ctx.closePath();
    ctx.fillStyle = fill;
    ctx.fill();
    // A hairline of the same colour closes the seams antialiasing leaves between neighbours.
    ctx.strokeStyle = fill;
    ctx.lineWidth = 0.6;
    ctx.stroke();
    q.then?.(ctx, pr.to);
  }
}

/**
 * The zero crossings of a bilinear cell with corner values f (counter-
 * clockwise from (0, 0)), as segments in the cell's unit square: marching
 * squares, with a saddle split by the cell's mean as `heatmap.contour` does.
 */
export function cellZeros(f: readonly [number, number, number, number]): [number, number, number, number][] {
  const c = [[0, 0], [1, 0], [1, 1], [0, 1]];
  const pts: [number, number][] = [];
  for (let e = 0; e < 4; e++) {
    const f0 = f[e], f1 = f[(e + 1) % 4];
    if ((f0 < 0) !== (f1 < 0)) {
      const t = f0 / (f0 - f1), a = c[e], b = c[(e + 1) % 4];
      pts.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
    }
  }
  if (pts.length === 2) return [[...pts[0], ...pts[1]]];
  if (pts.length === 4) {
    const centre = (f[0] + f[1] + f[2] + f[3]) / 4, isolate0 = (centre < 0) !== (f[0] < 0);
    const [p, q, r, s] = isolate0 ? [pts[0], pts[3], pts[1], pts[2]] : [pts[0], pts[1], pts[2], pts[3]];
    return [[...p, ...q], [...r, ...s]];
  }
  return [];
}

/** Drag on a canvas to turn its camera; `onTurn` asks for a frame. */
export function orbit(canvas: HTMLCanvasElement, cam: Camera, onTurn: () => void): void {
  let last: [number, number] | null = null;
  canvas.addEventListener("pointerdown", (e) => { last = [e.clientX, e.clientY]; canvas.setPointerCapture(e.pointerId); });
  canvas.addEventListener("pointerup", () => { last = null; });
  canvas.addEventListener("pointercancel", () => { last = null; });
  canvas.addEventListener("pointermove", (e) => {
    if (!last) return;
    cam.yaw += (e.clientX - last[0]) * 0.008;
    cam.pitch = Math.min(1.35, Math.max(0.05, cam.pitch + (e.clientY - last[1]) * 0.006));
    last = [e.clientX, e.clientY];
    onTurn();
  });
}
