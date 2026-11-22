// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// The top of the slab, on the exporter's fixed grid of radial lines. Every
// keyframe shares the grid, so a frame between two keyframes is a per-node
// blend: radius, slab-top temperature, and the signed distance to the slab's
// edge, whose zero contour the fragment shader cuts along — the slab grows
// smoothly rather than a keyframe at a time. Depth contours (every 10 km, and
// heavier every 50 km, as in the published figure) are drawn in the same
// shader from the interpolated depth, so they cost no geometry.

import {
  BufferAttribute,
  BufferGeometry,
  Color,
  DoubleSide,
  Mesh,
  MeshStandardMaterial,
  type Texture,
} from "three";
import type { Dataset } from "./data";
import { type Centre, toScene } from "./frame";

export type SlabColour = "titanium" | "temperature" | "depth";

/** Radius of the domain's top surface, km: depth is measured from here. */
export const TOP_RADIUS = 6370;
/** Depth the depth colour map spans, km (the slab top reaches ~224 km). */
export const DEPTH_SPAN = 225;

export class Slab {
  readonly mesh: Mesh<BufferGeometry, MeshStandardMaterial>;
  private readonly dirs: Float64Array;
  private readonly pos: Float32Array;
  private readonly nrm: Float32Array;
  private readonly edge: Float32Array;
  private readonly temp: Float32Array;
  private readonly depth: Float32Array;
  private readonly uniforms = {
    uMode: { value: 0 },
    uContours: { value: 1 },
    uDepthMax: { value: DEPTH_SPAN },
    uTempLut: { value: null as Texture | null },
    uDepthLut: { value: null as Texture | null },
  };

  constructor(
    private readonly data: Dataset,
    private readonly centre: Centre,
    tempLut: Texture,
    depthLut: Texture,
  ) {
    const s = data.manifest.surface;
    const { nx, ny } = s;
    const N = nx * ny;
    this.dirs = new Float64Array(3 * N);
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const xi = s.origin[0] + s.spacing * i;
        const eta = s.origin[1] + s.spacing * j;
        const len = Math.hypot(xi, eta, s.R0);
        const o = 3 * (j * nx + i);
        this.dirs[o] = xi / len;
        this.dirs[o + 1] = eta / len;
        this.dirs[o + 2] = s.R0 / len;
      }
    }
    this.pos = new Float32Array(3 * N);
    this.nrm = new Float32Array(3 * N);
    this.edge = new Float32Array(N);
    this.temp = new Float32Array(N);
    this.depth = new Float32Array(N);

    // Two triangles per cell, wound so the front face looks up (cross of the
    // grid's +i and +j directions is the outward radial).
    const index = new Uint32Array(6 * (nx - 1) * (ny - 1));
    let q = 0;
    for (let j = 0; j < ny - 1; j++) {
      for (let i = 0; i < nx - 1; i++) {
        const a = j * nx + i;
        const b = a + 1;
        const c = a + nx;
        const d = c + 1;
        index.set([a, b, d, a, d, c], q);
        q += 6;
      }
    }
    const g = new BufferGeometry();
    g.setIndex(new BufferAttribute(index, 1));
    g.setAttribute("position", new BufferAttribute(this.pos, 3));
    g.setAttribute("normal", new BufferAttribute(this.nrm, 3));
    g.setAttribute("aEdge", new BufferAttribute(this.edge, 1));
    g.setAttribute("aTemp", new BufferAttribute(this.temp, 1));
    g.setAttribute("aDepth", new BufferAttribute(this.depth, 1));

    this.uniforms.uTempLut.value = tempLut;
    this.uniforms.uDepthLut.value = depthLut;
    const mat = new MeshStandardMaterial({ side: DoubleSide, transparent: true });
    mat.onBeforeCompile = (shader) => {
      Object.assign(shader.uniforms, this.uniforms);
      shader.vertexShader = shader.vertexShader
        .replace(
          "#include <common>",
          `#include <common>
attribute float aEdge;
attribute float aTemp;
attribute float aDepth;
varying float vEdge;
varying float vTemp;
varying float vDepth;`,
        )
        .replace(
          "#include <begin_vertex>",
          `#include <begin_vertex>
vEdge = aEdge;
vTemp = aTemp;
vDepth = aDepth;`,
        );
      shader.fragmentShader = shader.fragmentShader
        .replace(
          "#include <common>",
          `#include <common>
uniform sampler2D uTempLut;
uniform sampler2D uDepthLut;
uniform int uMode;
uniform float uContours;
uniform float uDepthMax;
varying float vEdge;
varying float vTemp;
varying float vDepth;
// Coverage of a contour line of the given width (px) through v = k * spacing.
float isoline(float v, float spacing, float px) {
  float f = v / spacing;
  float d = abs(fract(f + 0.5) - 0.5) / max(fwidth(f), 1e-6);
  return 1.0 - smoothstep(0.5 * px - 0.5, 0.5 * px + 0.5, d);
}`,
        )
        .replace(
          "#include <color_fragment>",
          `#include <color_fragment>
if (vEdge < 0.0) discard;
vec3 base = diffuseColor.rgb;
if (uMode == 1) base = texture2D(uTempLut, vec2(clamp(vTemp, 0.0, 1.0), 0.5)).rgb;
else if (uMode == 2) base = texture2D(uDepthLut, vec2(clamp(vDepth / uDepthMax, 0.0, 1.0), 0.5)).rgb;
float lines = max(0.55 * isoline(vDepth, 10.0, 1.0), isoline(vDepth, 50.0, 2.2)) * uContours;
diffuseColor.rgb = mix(base, uMode == 0 ? vec3(1.0) : vec3(0.12), lines);`,
        );
    };
    mat.customProgramCacheKey = () => "slab-top";
    this.mesh = new Mesh(g, mat);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = 10;
    this.setColour("titanium");
  }

  /** Blend keyframes k and k + 1 with weight s on the later. */
  update(k: number, s: number): void {
    const d = this.data;
    const m = d.manifest.surface;
    const { nx, ny } = m;
    const N = nx * ny;
    const a = k * N;
    const b = a + N;
    const [lo, hi] = m.radius.range;
    const rScale = (hi - lo) / 65535;
    const eScale = m.edgeDistance.kmPerUnit;
    const R = d.radius;
    const T = d.slabTemperature;
    const E = d.edge;
    const w = 1 - s;
    const dirs = this.dirs;
    for (let p = 0; p < N; p++) {
      const r = lo + (R[a + p] * w + R[b + p] * s) * rScale;
      const o = 3 * p;
      toScene(this.centre, r * dirs[o], r * dirs[o + 1], r * dirs[o + 2], this.pos, o);
      this.edge[p] = (E[a + p] * w + E[b + p] * s) * eScale;
      this.temp[p] = (T[a + p] * w + T[b + p] * s) / 255;
      this.depth[p] = TOP_RADIUS - r;
    }
    // Normals by central differences along the grid (one-sided at its border).
    const P = this.pos;
    for (let j = 0; j < ny; j++) {
      const jl = Math.max(j - 1, 0) * nx;
      const jr = Math.min(j + 1, ny - 1) * nx;
      for (let i = 0; i < nx; i++) {
        const il = 3 * (j * nx + Math.max(i - 1, 0));
        const ir = 3 * (j * nx + Math.min(i + 1, nx - 1));
        const dl = 3 * (jl + i);
        const dr = 3 * (jr + i);
        const ux = P[ir] - P[il], uy = P[ir + 1] - P[il + 1], uz = P[ir + 2] - P[il + 2];
        const vx = P[dr] - P[dl], vy = P[dr + 1] - P[dl + 1], vz = P[dr + 2] - P[dl + 2];
        let nx_ = uy * vz - uz * vy;
        let ny_ = uz * vx - ux * vz;
        let nz_ = ux * vy - uy * vx;
        const len = Math.hypot(nx_, ny_, nz_) || 1;
        nx_ /= len;
        ny_ /= len;
        nz_ /= len;
        const o = 3 * (j * nx + i);
        this.nrm[o] = nx_;
        this.nrm[o + 1] = ny_;
        this.nrm[o + 2] = nz_;
      }
    }
    const g = this.mesh.geometry;
    for (const name of ["position", "normal", "aEdge", "aTemp", "aDepth"]) g.getAttribute(name).needsUpdate = true;
  }

  setColour(c: SlabColour): void {
    const mat = this.mesh.material;
    this.uniforms.uMode.value = c === "titanium" ? 0 : c === "temperature" ? 1 : 2;
    if (c === "titanium") {
      // ParaView's "Metal_Titanium_brushed", as far as a rasteriser goes.
      mat.color = new Color(0xa4aab3);
      mat.metalness = 0.7;
      mat.roughness = 0.36;
    } else {
      mat.color = new Color(0xffffff);
      mat.metalness = 0.0;
      mat.roughness = 0.65;
    }
  }

  setOpacity(o: number): void {
    const mat = this.mesh.material;
    mat.opacity = o;
    mat.transparent = o < 1;
    mat.depthWrite = true;
    mat.needsUpdate = true;
  }

  setContours(on: boolean): void {
    this.uniforms.uContours.value = on ? 1 : 0;
  }
}
