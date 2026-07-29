import { defineConfig } from "vitest/config";

// Verification runs are real solves (convergence studies, runs to steady state),
// not unit-test stubs — they need far more than the 5 s default.
export default defineConfig({
  test: {
    testTimeout: 180_000,
    hookTimeout: 180_000,
    // A stray `jekyll build` run from this directory copies the suite into
    // `_site/`, where Vitest happily finds and runs it a second time — silently
    // doubling the run and inflating the count. It is git-ignored, so nothing
    // else notices.
    exclude: ["**/node_modules/**", "**/dist/**", "**/_site/**"],
  },
});
