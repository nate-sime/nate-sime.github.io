/**
 * The selected viscosity law, written out under the list that selects it.
 *
 * The two sliders name quantities that do not appear in the equation they act
 * on. `log₁₀ contrast` sets γ = ln(contrast) — a different number, an order of
 * magnitude smaller — and `power-law n` enters only through the exponent
 * (1−n)/n, whose *sign* is the whole behaviour: negative for n > 1, which is
 * why the fluid thins where it shears. Neither mapping is guessable from a
 * slider, so the folder states the law it is currently solving and names, symbol
 * by symbol, which control sets what and what that symbol's value is right now.
 *
 * DOM-free on purpose. The formulas and the control names are the part worth
 * regression-testing — a legend that names a slider the pane no longer has is
 * exactly the failure this file exists to prevent — and that test should not
 * need a browser. `controls.ts` renders it.
 */

import { EPS_MIN, gammaFor } from "../solver/rheology";
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

const N: EquationParam = {
  sym: "n",
  control: LABELS.n,
  value: (s) => String(s.n),
};

export const EQUATION: Record<ViscosityName, Equation> = {
  "constant": {
    lines: ["μ = 1"],
    params: [],
    note: "Nothing below applies; the Stokes operator has no coefficient to vary.",
  },
  "μ(T)": {
    lines: ["μ(T) = exp(−γ (T − ½))"],
    params: [GAMMA],
    // Centred on T = ½, so the contrast is what the slider says it is: the ratio
    // across the layer, not a rescaling of the whole viscosity.
    note: "γ = ln(contrast) = ln(μ|T=0 ÷ μ|T=1).",
  },
  "μ(T, ε̇)": {
    lines: [
      "μ(T, ε̇) = exp(−γ (T − ½)) · ŝ^{(1−n)/n}",
      "ŝ = (ε̇ + δ) / G",
      "clamped to [e^{−γ/2}, e^{+γ/2}]",
    ],
    params: [GAMMA, N],
    note: `ε̇ is the second invariant of the strain rate, G its geometric mean `
      + `over the domain and δ = ${EPS_MIN} · rms ε̇; both come from the flow. `
      + `The exponent is negative for n > 1, so shear thins.`,
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
