// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { execSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { defineConfig } from "vite";

// `version` is the one manually-bumped number; the hash and date are stamped
// on at build time, as in apps/mantle, so every bundle that leaves the repo is
// traceable to the commit and day it was built from.
const { version } = JSON.parse(
  readFileSync(new URL("./package.json", import.meta.url), "utf-8"),
) as { version: string };

const commit = (() => {
  try {
    return execSync("git rev-parse --short HEAD").toString().trim();
  } catch {
    return "unknown";
  }
})();
const date = new Date().toISOString().slice(0, 10);

// Served under /assets/subduction-3d/ on the Jekyll site; `build` emits directly
// into the assets folder that GitHub Pages serves. public/data/ (written by
// tools/export_data.py) is copied across as it is.
export default defineConfig({
  base: "/assets/subduction-3d/",
  define: { __APP_VERSION__: JSON.stringify(`${version}+${date}.${commit}`) },
  // three.js is most of the one chunk (~185 kB gzipped); splitting it buys nothing.
  build: { outDir: "../../assets/subduction-3d", emptyOutDir: true, chunkSizeWarningLimit: 900 },
});
