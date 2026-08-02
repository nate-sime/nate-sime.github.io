/**
 * The selected viscosity law, written out under the list that selects it.
 *
 * The sliders name quantities that do not appear in the equation they act on.
 * `log₁₀ contrast` and `log₁₀ depth contrast` set γ and c = ln(contrast) —
 * different numbers, an order of magnitude smaller — and `power-law n` enters
 * only through the exponent (1−n)/n, whose *sign* is the whole behaviour:
 * negative for n > 1, which is why the fluid thins where it shears. None of
 * those mappings is guessable from a slider, so the folder states the law it is
 * currently solving and names, symbol by symbol, which control sets what and
 * what that symbol's value is right now.
 *
 * DOM-free on purpose. The formulas and the control names are the part worth
 * regression-testing — a legend that names a slider the pane no longer has is
 * exactly the failure this file exists to prevent — and that test should not
 * need a browser. `controls.ts` renders it.
 */

import { EPS_MIN, gammaFor, TACKLEY_TRANSITION_DEPTH } from "../solver/rheology";
import { LABELS, type State, type ViscosityName } from "./presets";

/** A symbol in the formula, and the control that sets it. */
export interface EquationParam {
  /** The symbol, exactly as it is written above. */
  sym: string;
  /** Label of the slider that sets it, verbatim from the pane. */
  control: string;
  /** Its current value *in the equation's terms*, not the slider's. */
  value: (s: State) => string;
}

export interface Equation {
  /** The law, one formula per line. `^{…}` marks a superscript. */
  lines: string[];
  /** The symbols the pane controls, in the order they are written. */
  params: EquationParam[];
  /** What the remaining symbols are; they are set by the flow, not the user. */
  note?: string;
}

const GAMMA: EquationParam = {
  sym: "γ",
  control: LABELS.contrast,
  value: (s) => gammaFor(10 ** s.logContrast).toFixed(2),
};

/**
 * The depth term's coefficient. Written `c` because that is the letter the
 * Blankenbach cases state it with, and displayed as γ is — the slider carries a
 * log₁₀ and the equation an ln.
 */
const C: EquationParam = {
  sym: "c",
  control: LABELS.depth,
  value: (s) => gammaFor(10 ** s.logDepthContrast).toFixed(2),
};

const N: EquationParam = {
  sym: "n",
  control: LABELS.n,
  value: (s) => String(s.n),
};

/** What `d` is, said wherever the depth term appears. */
const DEPTH = "d is depth, 0 at the cold boundary and 1 at the hot one, so c > 0 "
  + "stiffens the deep interior.";

const SIGMA_Y: EquationParam = {
  sym: "σ_Y",
  control: LABELS.sigmaY,
  value: (s) => s.sigmaY.toFixed(2),
};

const SIGMA_B: EquationParam = {
  sym: "σ_b",
  control: LABELS.sigmaB,
  value: (s) => s.sigmaB.toFixed(2),
};

const ETA_STAR: EquationParam = {
  sym: "η*",
  control: LABELS.etaStar,
  value: (s) => s.etaStar.toExponential(1),
};

export const EQUATION: Record<ViscosityName, Equation> = {
  "constant": {
    lines: ["μ = 1"],
    params: [],
    note: "Nothing below applies; the Stokes operator has no coefficient to vary.",
  },
  "μ(T, d)": {
    lines: ["μ(T, d) = exp(−γ (T − ½) + c (d − ½))"],
    params: [GAMMA, C],
    // Centred on T = ½ and d = ½, so each contrast is what its slider says it
    // is: a ratio across the layer, not a rescaling of the whole viscosity.
    note: `γ = ln(contrast) = ln(μ|T=0 ÷ μ|T=1), c = ln(depth contrast). ${DEPTH}`,
  },
  "μ(T, d, ε̇)": {
    lines: [
      "μ(T, d, ε̇) = exp(−γ (T − ½) + c (d − ½)) · ŝ^{(1−n)/n}",
      "ŝ = (ε̇ + δ) / G",
      "clamped to [e^{−(γ+c)/2}, e^{+(γ+c)/2}]",
    ],
    params: [GAMMA, C, N],
    note: `ε̇ is the second invariant of the strain rate, G its geometric mean `
      + `over the domain and δ = ${EPS_MIN} · rms ε̇; both come from the flow. `
      + `The exponent is negative for n > 1, so shear thins. ${DEPTH}`,
  },
  "Tackley": {
    lines: [
      "μ_lin(T, d) = A₀(d) exp(27.631/(T + 0.88)) · 5.86052e-13",
      "μ_plast(d, ε̇) = η* + (σ_Y + σ_b d)/ε̇",
      "μ = (μ_lin⁻¹ + μ_plast⁻¹)⁻¹",
    ],
    params: [SIGMA_Y, SIGMA_B, ETA_STAR],
    note: `Diffusion creep in parallel with Bingham yielding — the weaker branch `
      + `sets μ. A₀ = 1 above 670 km, 30 below (d = ${TACKLEY_TRANSITION_DEPTH.toFixed(3)}). `
      + `ε̇ is the strain rate itself, unnormalised: the yield stress is an `
      + `absolute threshold. ${DEPTH}`,
  },
  "Tosi": {
    lines: [
      "μ_lin(T, d) = exp(−γ (T − ½) + c (d − ½))",
      "μ_plast(d, ε̇) = η* + (σ_Y + σ_b d)/ε̇",
      "μ = 2 (μ_lin⁻¹ + μ_plast⁻¹)⁻¹",
    ],
    params: [GAMMA, C, SIGMA_Y, SIGMA_B, ETA_STAR],
    note: `Tosi et al. (2015): the μ(T, d) exponential (unclamped) harmonically `
      + `averaged with the same Bingham yielding Tackley states — the weaker `
      + `branch dominates, and the factor of 2 is what makes it a genuine `
      + `average of the two branches' conductances rather than half of one. `
      + `ε̇ is the strain rate itself, unnormalised: the yield `
      + `stress is an absolute threshold. ${DEPTH}`,
  },
};

/**
 * Split a formula into runs of ordinary and superscript text. `^{…}` is the
 * only markup, and an unclosed one is left as literal text rather than being
 * treated as an error: the strings above are the only input this ever sees.
 */
export const parseFormula = (line: string): { text: string; sup: boolean }[] => {
  const re = /\^\{([^}]*)\}/g;
  const out: { text: string; sup: boolean }[] = [];
  let i = 0;
  for (let m = re.exec(line); m; m = re.exec(line)) {
    if (m.index > i) out.push({ text: line.slice(i, m.index), sup: false });
    out.push({ text: m[1], sup: true });
    i = m.index + m[0].length;
  }
  if (i < line.length) out.push({ text: line.slice(i), sup: false });
  return out;
};
