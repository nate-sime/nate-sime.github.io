/**
 * The equation shown under the viscosity list. No DOM: what is checked here is
 * that the legend and the pane agree, and that is a property of the tables, not
 * of the rendering.
 *
 * The failure this guards against is silent and permanent — a law gains a knob,
 * or a slider is renamed, and the legend goes on naming a control that is not
 * there or omitting one that is. Nothing throws, the equation simply starts
 * lying, and only a reader who already knew the rheology would notice.
 */

import { describe, it, expect } from "vitest";
import { EQUATION, parseFormula } from "../src/ui/equation";
import { LABELS, VISCOSITY, defaultState, type ViscosityName } from "../src/ui/presets";
import { gammaFor } from "../src/solver/rheology";

const laws = Object.keys(VISCOSITY) as ViscosityName[];
const text = (line: string) => parseFormula(line).map((r) => r.text).join("");

describe("viscosity equation", () => {
  it.each(laws)("%s is written out", (law) => {
    const eq = EQUATION[law];
    expect(eq.lines.length).toBeGreaterThan(0);
    expect(eq.lines[0]).toContain("μ");
  });

  it.each(laws)("%s names every symbol it parameterises", (law) => {
    const eq = EQUATION[law];
    const formula = eq.lines.map(text).join(" ");
    for (const p of eq.params) expect(formula).toContain(p.sym);
  });

  // The legend's whole job. A label that is not the slider's is worse than no
  // legend: it sends the reader looking for a control the pane does not have.
  it.each(laws)("%s points at controls the pane actually shows", (law) => {
    for (const p of EQUATION[law].params)
      expect(Object.values(LABELS) as string[]).toContain(p.control);
  });

  /**
   * Tied to `VISCOSITY`, not written out twice: the same table decides which
   * sliders `controls.ts` enables, so a law that enables a knob without the
   * equation containing its symbol would leave the user adjusting a quantity
   * the displayed law does not mention.
   */
  it.each(laws)("%s parameterises exactly what its tier enables", (law) => {
    const { variable, strainRate, tackley, tosi } = VISCOSITY[law];
    const by = Object.fromEntries(EQUATION[law].params.map((p) => [p.control, p]));
    // Tackley is parameterised by σ_Y, σ_b, η* instead of γ, c, n — a different
    // law on the same Krylov tier, not a further case of the power law's
    // contrast/depth/n (see `presets.ts`). Tosi keeps contrast/depth (they are
    // the paper's own γ_T, γ_z) but likewise has no n, and shares σ_Y, σ_b, η*
    // with Tackley.
    expect(LABELS.contrast in by).toBe(variable && !tackley);
    // The depth contrast is enabled by the *tier*, not by the power law: it is a
    // slider on both non-Tackley variable laws, which is the distinction
    // `VISCOSITY` draws and the reason the list has five entries rather than
    // nine.
    expect(LABELS.depth in by).toBe(variable && !tackley);
    expect(LABELS.n in by).toBe(strainRate && !tackley && !tosi);
    expect(LABELS.sigmaY in by).toBe(tackley || tosi);
    expect(LABELS.sigmaB in by).toBe(tackley || tosi);
    expect(LABELS.etaStar in by).toBe(tackley || tosi);
  });

  it("shows γ, which is not the number on the contrast slider", () => {
    const s = { ...defaultState(), viscosity: "μ(T, d)" as const, logContrast: 3 };
    const [gamma] = EQUATION["μ(T, d)"].params;
    expect(gamma.sym).toBe("γ");
    expect(gamma.value(s)).toBe(gammaFor(1e3).toFixed(2));   // 6.91, not 3
    expect(gamma.value(s)).not.toBe(String(s.logContrast));
  });

  // Same mapping, second slider — and the two must not be confused for each
  // other, which is exactly what a legend copied from γ would do.
  it("shows c from the depth slider, and 0 when there is no depth term", () => {
    const c = EQUATION["μ(T, d)"].params.find((p) => p.control === LABELS.depth)!;
    expect(c.sym).toBe("c");
    expect(c.value({ ...defaultState(), logDepthContrast: 2 }))
      .toBe(gammaFor(100).toFixed(2));                       // 4.61, not 2
    // The default: no depth dependence at all, said as a number rather than by
    // the term quietly vanishing from the equation.
    expect(c.value(defaultState())).toBe("0.00");
    expect(defaultState().logDepthContrast).toBe(0);
  });

  it("shows n as the slider sets it", () => {
    const s = { ...defaultState(), n: 3.25 };
    const n = EQUATION["μ(T, d, ε̇)"].params.find((p) => p.control === LABELS.n)!;
    expect(n.value(s)).toBe("3.25");
  });
});

describe("formula markup", () => {
  it("strips to the formula without its markup", () => {
    expect(text("μ = ŝ^{(1−n)/n}")).toBe("μ = ŝ(1−n)/n");
    expect(parseFormula("μ = ŝ^{(1−n)/n}")).toEqual([
      { text: "μ = ŝ", sup: false },
      { text: "(1−n)/n", sup: true },
    ]);
  });

  it("handles several superscripts, and none", () => {
    expect(parseFormula("[e^{−γ/2}, e^{+γ/2}]").filter((r) => r.sup))
      .toEqual([{ text: "−γ/2", sup: true }, { text: "+γ/2", sup: true }]);
    expect(parseFormula("μ = 1")).toEqual([{ text: "μ = 1", sup: false }]);
    expect(parseFormula("")).toEqual([]);
  });

  // Every run is rendered, so a `^` or `{` surviving into one would be printed
  // literally in the pane.
  it.each(laws)("%s leaves no markup in the rendered runs", (law) => {
    for (const line of EQUATION[law].lines)
      for (const run of parseFormula(line)) expect(run.text).not.toMatch(/[\^{}]/);
  });
});
