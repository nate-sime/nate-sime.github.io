// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    // A stray `jekyll build` run from this directory copies the suite into
    // `_site/`, where Vitest would find and run it a second time.
    exclude: ["**/node_modules/**", "**/dist/**", "**/_site/**"],
  },
});
