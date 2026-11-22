// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// Temperature on the planes y = -200, 0, 200 km. Each is one quad over the
// sampled rectangle, textured from a layer of a 2D array texture holding every
// keyframe of every slice; the shader blends the two keyframes either side of
// the current time and cuts along the domain's edge (see prepareSlices).

import {
  BufferAttribute,
  BufferGeometry,
  DataArrayTexture,
  DoubleSide,
  Group,
  LinearFilter,
  Mesh,
  RGFormat,
  ShaderMaterial,
  type Texture,
  UnsignedByteType,
} from "three";
import { type Dataset, prepareSlices } from "./data";
import { type Centre, toScene } from "./frame";

const vertexShader = /* glsl */ `
varying vec2 vUv;
void main() {
  vUv = uv;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const fragmentShader = /* glsl */ `
precision highp sampler2DArray;
uniform sampler2DArray uTex;
uniform sampler2D uLut;
uniform float uA;
uniform float uB;
uniform float uS;
varying vec2 vUv;
void main() {
  vec2 v = mix(texture(uTex, vec3(vUv, uA)).rg, texture(uTex, vec3(vUv, uB)).rg, uS) * 255.0;
  if (v.g < 128.0) discard;
  float f = clamp((v.r - 1.0) / 254.0, 0.0, 1.0);
  gl_FragColor = vec4(texture2D(uLut, vec2(f, 0.5)).rgb, 1.0);
  #include <colorspace_fragment>
}`;

export class Slices {
  readonly group = new Group();
  readonly meshes: Mesh<BufferGeometry, ShaderMaterial>[] = [];
  private readonly count: number;

  constructor(data: Dataset, centre: Centre, lut: Texture) {
    const m = data.manifest.slices;
    const K = data.manifest.keyframes.length;
    this.count = m.y.length;
    const tex = new DataArrayTexture(prepareSlices(data.slices, m.nx, m.nz, K * this.count), m.nx, m.nz, K * this.count);
    tex.format = RGFormat;
    tex.type = UnsignedByteType;
    tex.minFilter = tex.magFilter = LinearFilter;
    tex.generateMipmaps = false;
    tex.unpackAlignment = 1;
    tex.needsUpdate = true;

    // Texel centres sit half a texel in from the texture's edge; the quad's
    // corners are the first and last samples, so its uv range is inset to match.
    const u0 = 0.5 / m.nx, u1 = 1 - u0, v0 = 0.5 / m.nz, v1 = 1 - v0;
    for (const y of m.y) {
      const pos = new Float32Array(12);
      toScene(centre, m.x[0], y, m.z[0], pos, 0);
      toScene(centre, m.x[1], y, m.z[0], pos, 3);
      toScene(centre, m.x[0], y, m.z[1], pos, 6);
      toScene(centre, m.x[1], y, m.z[1], pos, 9);
      const g = new BufferGeometry();
      g.setAttribute("position", new BufferAttribute(pos, 3));
      g.setAttribute("uv", new BufferAttribute(new Float32Array([u0, v0, u1, v0, u0, v1, u1, v1]), 2));
      g.setIndex([0, 1, 3, 0, 3, 2]);
      const mat = new ShaderMaterial({
        uniforms: { uTex: { value: tex }, uLut: { value: lut }, uA: { value: 0 }, uB: { value: 0 }, uS: { value: 0 } },
        vertexShader,
        fragmentShader,
        side: DoubleSide,
      });
      const mesh = new Mesh(g, mat);
      mesh.frustumCulled = false;
      this.meshes.push(mesh);
      this.group.add(mesh);
    }
  }

  update(k: number, s: number): void {
    this.meshes.forEach((mesh, i) => {
      const u = mesh.material.uniforms;
      u.uA.value = k * this.count + i;
      u.uB.value = (k + 1) * this.count + i;
      u.uS.value = s;
    });
  }
}
