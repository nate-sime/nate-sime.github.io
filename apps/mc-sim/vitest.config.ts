// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { defineConfig } from "vitest/config";

// Verification runs will be real solves (convergence studies, sampling runs),
// not unit-test stubs — they need far more than the 5 s default.
export default defineConfig({
  test: {
    testTimeout: 180_000,
    hookTimeout: 180_000,
    // A stray `jekyll build` run from this directory copies the suite into
    // `_site/`, where Vitest would find and run it a second time.
    exclude: ["**/node_modules/**", "**/dist/**", "**/_site/**"],
  },
});
