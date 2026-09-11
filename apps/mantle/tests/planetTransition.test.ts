import { describe, expect, it } from "vitest";
import { PlanetTransitionController, nextPhase, phaseDuration } from "../src/ui/planetTransition";

describe("planet transition controller", () => {
  it("has one explicit legal route back to focused", () => {
    expect(nextPhase("closing")).toBe("departing");
    expect(nextPhase("departing")).toBe("overview");
    expect(nextPhase("overview")).toBe("arriving");
    expect(nextPhase("arriving")).toBe("rebuilding");
    expect(nextPhase("rebuilding")).toBe("revealing");
    expect(nextPhase("revealing")).toBe("focused");
  });

  it("makes the final request win over stale callbacks", () => {
    const t = new PlanetTransitionController();
    const earthToVenus = t.request();
    const venusToEarth = t.request();
    expect(t.advance(earthToVenus)).toBeNull();
    expect(t.phase).toBe("closing");
    expect(t.finish(earthToVenus)).toBe(false);
    expect(t.finish(venusToEarth)).toBe(true);
    expect(t.phase).toBe("focused");
  });

  it("keeps the rebuild hold but reduces every camera phase to a frame", () => {
    expect(phaseDuration("departing", true)).toBe(1);
    expect(phaseDuration("rebuilding", true)).toBe(0);
    expect(phaseDuration("departing", false)).toBeGreaterThan(1);
  });
});
