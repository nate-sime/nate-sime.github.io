// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * View 2: the two structures side by side — the beam (`beam.ts`) on the left,
 * the Kirchhoff plate (`plate.ts`) on the right — each solved for the same
 * load, degree and continuity. The top row is their static deflections, the
 * bottom row what they do in motion: the same mode number, or the steady
 * response to the same forcing. The readout is the beam's, then the plate's.
 */

import type { Figure } from "../figure";
import type { State } from "../state";
import { renderBeam } from "./beam";
import { renderPlate } from "./plate";
import type { ViewResult } from "./view";

export function renderStructures(fig: Figure, st: State, t: number): ViewResult {
  const [beamStatic, plateStatic, beamMotion, plateMotion] = fig.panels(4, { cols: 2 });
  const beam = renderBeam([beamStatic, beamMotion], st, t);
  const plate = renderPlate([plateStatic, plateMotion], st, t);
  return {
    readout: ["BEAM", beam.readout, "", "PLATE", plate.readout].join("\n"),
    animate: beam.animate || plate.animate,
  };
}
