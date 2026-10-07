# Monte Carlo for vibrating structures — implementation plan

A teaching app, built the way `apps/mantle` is, that shows Monte Carlo error
estimation for functions and functionals of vibrating systems (beams, then
plates) discretised with arbitrary-order spline finite elements. It is meant to
teach Monte Carlo, hierarchies of discretisation error, and multilevel Monte
Carlo (MLMC).

## What carries over from mantle

- **Scaffold:** Vite + TypeScript + Tweakpane + vitest. The build goes to
  `../../assets/<app>/`, with a Jekyll page that embeds it and the existing
  deploy workflow.
- **Spline code:** the Piegl & Tiller basis routines in `src/spline.ts`. That
  file fixes the degree with `export const P = 3` and is written around a
  periodic axis, so it needs **rewriting for arbitrary degree p and adjustable
  continuity**. Copy it into the new app rather than share it, so mantle can't
  break.
- **Small helpers:** `quad.ts` (Gauss), `linalg.ts` (LU), the tour and overlay
  system (`ui/tour.ts`, `ui/tours.ts`), and the time-series plot
  (`ui/nuplot.ts`), which can be reused for the convergence and MLMC plots.
- **Testing habit:** CPU reference first, checked against analytic answers,
  with real measured numbers written into the README.

## Stages

Each stage ends with something you can run and check.

**Status (2026-10-07):** stages 0–3 done — see the README for what was
measured. Decided: CPU f64 reference, structured for a later WebGPU port;
solves nondimensional, with a dimensional/nondimensional display toggle.

### 0. Scaffold
Create `apps/<name>/` mirroring mantle's config, plus a placeholder site page.
Done when `npm run dev`, `npm test` and `npm run build` all work.

### 1. 1D splines of any degree
- Clamped knot vectors with chosen degree p and continuity (h-, p- and
  k-refinement).
- Basis derivatives up to order 2+; (p+1)-point Gauss quadrature.
- *Tests:* partition of unity, derivatives against finite differences, exact
  reproduction of polynomials up to degree p.

### 2. Euler–Bernoulli beam (deterministic)
- Build the stiffness matrix K and mass matrix M. C¹ splines are the natural fit
  for 4th-order problems; this is the first selling point of the app.
- Ends: clamped, simply supported, cantilever.
- Static solve (banded Cholesky); free vibration `Kφ = λMφ`.
- *Tests:* analytic natural frequencies (the βL roots) and a manufactured static
  solution. Rates: about h^{p+1} in L², and h^{2(p−1)} for eigenvalues.

### 3. Hierarchy of discretisation error
- Level sequence h_ℓ = h₀·2^{−ℓ}. Measure how the quantity of interest Q_ℓ
  changes, |Q_ℓ − Q_{ℓ−1}| ~ h^α, and how cost grows, ~ h^{−γ}.
- First teaching view: convergence plots across h, p and k refinement.
- Show the known spline-versus-C⁰ high-frequency spectrum result
  (Cottrell/Hughes "outlier" modes).

### 4. Random inputs
- Seeded PRNG (PCG or xoshiro) so every sample can be reproduced.
- Stiffness EI(x, ω) as a lognormal field from a truncated Karhunen–Loève
  expansion; optionally a random load too.
- **Key requirement:** each sample ω must be usable on every level, because MLMC
  needs the same ω on the fine and coarse solves.

### 5. Plain Monte Carlo
- Quantities of interest:
  - functions: mean and variance of the deflection field, with bands;
  - functionals: tip deflection, compliance, first natural frequency, response
    amplitude at a forcing frequency.
- Histogram, running mean with CLT confidence intervals, the σ/√N rate.
- Teaching point: MSE = bias² (from discretisation) + variance (from sampling).
- Samples are independent, so run them in Web Workers.

### 6. Multilevel Monte Carlo
- Telescoping estimator with coupled corrections Y_ℓ = Q_ℓ − Q_{ℓ−1}; optimal
  N_ℓ allocation; Giles' adaptive algorithm that estimates α, β and γ while it
  runs.
- Standard plots: log of mean and variance against level, N_ℓ against level,
  cost against ε for MC versus MLMC.
- *Tests:* the telescoping sum is unbiased; with a fixed seed the measured β
  matches stage 3.

### 7. Kirchhoff plate
- Rectangular tensor-product splines (C¹ is again the natural fit).
- Verify a simply supported plate against the Navier series solution; add
  clamped edges.
- Mode shapes give Chladni-like patterns.
- 2D KL random field.
- Cost now matters: banded or sparse Cholesky, shift-invert Lanczos for a few
  eigenpairs. Consider WebGPU only if Workers turn out too slow — the
  parallelism is across samples, not within one solve.

### 8. Visualisation and UI
Animate the vibrating mode or response, overlay the mesh for each level, and
show the histogram and MLMC diagnostic panels next to the solver view.

### 9. Tours and write-up
Four tours:
- Monte Carlo and the √N rate
- discretisation hierarchies (h, p and k refinement)
- the bias–variance split
- MLMC and its complexity gain

Plus the formulation page on the site.

## Pitfalls to design around early

- **Eigenvalue quantities:** past the first mode, eigenvalues can cross between
  samples, so "the k-th frequency" isn't a smooth function of ω, and coarse and
  fine levels may pick different modes. Start with λ₁, or track modes by
  matching shapes.
- **Response near resonance:** variance blows up unless there is damping.
  Include Rayleigh damping from the start.
- **Coarse levels and the random field:** if the coarsest mesh can't resolve the
  KL modes, the early level variances V_ℓ won't decay. Pick h₀ to match the
  correlation length, and make that a visible lesson.

## Open decisions

- What to call the app (working name: mc-sim).
- GPU precision for batched solves: f32 cannot hold a fourth-order stiffness
  (λ_max/λ₁ ~ 10⁶–10⁹ at 64 elements); mixed precision with refinement, or
  double-f32 throughout — decide at stage 5.
