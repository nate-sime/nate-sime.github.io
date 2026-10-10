// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
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

// Served under /assets/mc-sim/ on the Jekyll site; `build` emits directly into
// the assets folder that GitHub Pages serves.
export default defineConfig({
  base: "/assets/mc-sim/",
  define: { __APP_VERSION__: JSON.stringify(`${version}+${date}.${commit}`) },
  build: { outDir: "../../assets/mc-sim", emptyOutDir: true },
  // The Monte Carlo workers are ES modules: they import the solver like the page does.
  worker: { format: "es" },
});
