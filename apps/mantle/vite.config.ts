import { execSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { defineConfig } from "vite";

// `version` is the one manually-bumped number (`npm version patch`, etc.).
// The hash and date are stamped on at build time instead, so every artefact
// that actually left this repo — a local `npm run build` or the one
// deploy.yml runs on push to master — is traceable to the commit and day it
// was built from, without anyone having to remember to bump anything.
const { version } = JSON.parse(
  readFileSync(new URL("./package.json", import.meta.url), "utf-8"),
) as { version: string };

const commit = (() => {
  try {
    return execSync("git rev-parse --short HEAD").toString().trim();
  } catch {
    // No .git (e.g. a source tarball) — the readout still shows a version,
    // just without build provenance.
    return "unknown";
  }
})();
const date = new Date().toISOString().slice(0, 10);

// Served under /assets/mantle/ on the Jekyll site; `build` emits directly into
// the committed assets folder that GitHub Pages serves.
export default defineConfig({
  base: "/assets/mantle/",
  define: { __APP_VERSION__: JSON.stringify(`${version}+${date}.${commit}`) },
  build: { outDir: "../../assets/mantle", emptyOutDir: true },
});
