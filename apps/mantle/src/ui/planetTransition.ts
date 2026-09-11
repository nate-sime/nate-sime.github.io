/**
 * Pure timing/state model for the illustrative planetary trip.  Keeping this
 * independent of WebGPU makes cancellation and reduced motion testable, and
 * prevents camera coordinates from becoming implicit lifecycle state.
 */
export type PlanetTransitionPhase =
  | "focused" | "closing" | "departing" | "overview" | "arriving"
  | "rebuilding" | "revealing";

export const TRANSITION_PHASES: readonly PlanetTransitionPhase[] = [
  "closing", "departing", "overview", "arriving", "rebuilding", "revealing",
];

export const phaseDuration = (phase: PlanetTransitionPhase, reduced: boolean): number => {
  if (reduced) return phase === "rebuilding" ? 0 : 1;
  switch (phase) {
    case "closing": return 420;
    case "departing": return 720;
    case "overview": return 520;
    case "arriving": return 720;
    case "revealing": return 480;
    default: return 0;
  }
};

export const nextPhase = (phase: PlanetTransitionPhase): PlanetTransitionPhase => {
  const at = TRANSITION_PHASES.indexOf(phase);
  return at < 0 || at === TRANSITION_PHASES.length - 1 ? "focused" : TRANSITION_PHASES[at + 1];
};

/** A monotonically increasing request token makes the most recent selection win. */
export class PlanetTransitionController {
  private token = 0;
  private phase_: PlanetTransitionPhase = "focused";
  get phase(): PlanetTransitionPhase { return this.phase_; }
  get traveling(): boolean { return this.phase_ !== "focused"; }
  request(): number { this.phase_ = "closing"; return ++this.token; }
  isCurrent(token: number): boolean { return token === this.token; }
  advance(token: number): PlanetTransitionPhase | null {
    if (!this.isCurrent(token)) return null;
    this.phase_ = nextPhase(this.phase_);
    return this.phase_;
  }
  finish(token: number): boolean {
    if (!this.isCurrent(token)) return false;
    this.phase_ = "focused";
    return true;
  }
}
