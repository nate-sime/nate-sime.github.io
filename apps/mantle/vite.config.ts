import { defineConfig } from "vite";

// Served under /assets/mantle/ on the Jekyll site; `build` emits directly into
// the committed assets folder that GitHub Pages serves.
export default defineConfig({
  base: "/assets/mantle/",
  build: { outDir: "../../assets/mantle", emptyOutDir: true },
});
