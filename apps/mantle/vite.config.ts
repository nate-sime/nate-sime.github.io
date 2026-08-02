import { readFileSync } from "node:fs";
import { defineConfig } from "vite";

// The single source of the version shown in the readout (main.ts) — bump it
// here, not at the display site, so `npm version` stays the one command that
// changes it.
const { version } = JSON.parse(
  readFileSync(new URL("./package.json", import.meta.url), "utf-8"),
) as { version: string };

// Served under /assets/mantle/ on the Jekyll site; `build` emits directly into
// the committed assets folder that GitHub Pages serves.
export default defineConfig({
  base: "/assets/mantle/",
  define: { __APP_VERSION__: JSON.stringify(version) },
  build: { outDir: "../../assets/mantle", emptyOutDir: true },
});
