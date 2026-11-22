# Subduction-3D

An interactive viewer for the 3D subduction-zone run of
[Sime et al. (2023)](http://arxiv.org/abs/2309.12895), Fig. 16 — the slab top
beneath the Chilean–Argentinean border over 11 Myr, with temperature on three
vertical slices and tracers with ~1 Myr pathlines. It replaces the rendered
video that was embedded in the research page, and is built with three.js
(WebGL 2) and Tweakpane.

Nothing is simulated here: everything is precomputed from the run's output by
`tools/export_data.py` and played back, blending between saved frames.

## Data

The source is the run `2a_pushforward_n100_fine_xc75.0`: 100 output steps
(0–99) over 11 Myr, each with the slab-top surface (`surface_data_*`) and the
temperature and velocity of the whole domain (`temperature_*`, `velocity_*`,
~23 GB in all). The `n25` and `n50` runs of the same model kept only the slab
top, so they cannot give the slices or the tracers. The domain follows the
slab and is remeshed at every step.

Regenerate `public/data/` with ParaView's Python (5.11 or later; it bundles
VTK, numpy, h5py and matplotlib):

    pvpython tools/export_data.py <path>/2a_pushforward_n100_fine_xc75.0 public/data

It takes a couple of minutes and writes ~2.6 MB of gzip-compressed arrays and
a `manifest.json` describing them. Its stages are `surface`, `volume`,
`tracers` and `pack`; each extraction is cached under `tools/.export-cache/`
(ignored by git), so `--only pack` repacks without reading the run again.

| file | what |
|---|---|
| `surface-*.bin.gz` | the slab top every 3rd step, resampled onto one fixed grid of radial lines (gnomonic coordinates), so frames can be blended: radius, slab-top temperature, and the signed distance to the slab's edge, whose zero contour is where the slab is cut off |
| `slices.bin.gz` | temperature on the planes y = −200, 0, 200 km every 3rd step, 2.5 km texels |
| `outline.bin.gz` | the domain's feature edges at every step, simplified |
| `tracers.bin.gz` | ParaView's `ParticleTracer` through the velocity at every step, seeded as in the published animation (`statev2temp.pvsm`): 77 points on the plane x = 0, reseeded every 25 steps |

The colour map, seeds, slice planes and initial camera are those of the
ParaView state that rendered the original animation.

## Layout

    src/
      main.ts       entry point: load, scene, controls, the clock
      data.ts       manifest, decoding, keyframe brackets, tracer trails, slice textures
      frame.ts      model (km, Earth-centred, up = +z) to scene coordinates (up = +y)
      slab.ts       the slab top: per-node blending, edge cut and depth contours in the shader
      slices.ts     the temperature slices: one 2D-array texture of every keyframe
      tracers.ts    tracer spheres and fading pathlines
      outline.ts    the domain's outline, per step
      colormap.ts   ParaView's "Reds" and the depth ramp
    tests/          decoding, bookkeeping, and checks on the exported data itself
    tools/          export_data.py

## Develop

    npm install
    npm run dev        # http://localhost:5173/assets/subduction-3d/
    npm test
    npm run build      # into ../../assets/subduction-3d/ (untracked; deploy.yml builds it)
