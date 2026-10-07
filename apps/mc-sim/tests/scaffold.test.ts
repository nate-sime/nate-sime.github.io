// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The bundle's address is written in three places that nothing else ties
 * together: Vite's `base` (every asset URL in the bundle), its `outDir` (where
 * the files land on disk) and the site page's iframe. Any one drifting gives a
 * blank frame on the deployed page and a working dev server, so the mismatch
 * would only show after a push.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { describe, expect, it } from "vitest";
import config from "../vite.config";

const app = resolve(fileURLToPath(new URL("..", import.meta.url)));
const site = resolve(app, "../..");

describe("scaffold", () => {
  const base = config.base!;

  it("builds into the folder its base URL names", () => {
    expect(base).toBe("/assets/mc-sim/");
    expect(resolve(app, config.build!.outDir!)).toBe(resolve(site, base.slice(1)));
  });

  it("is embedded by the site page at that URL", () => {
    const page = readFileSync(resolve(site, "mc-sim.md"), "utf-8");
    expect(page).toContain(`<iframe src="${base}"`);
  });

  it("is kept out of git, as a build artefact", () => {
    const ignore = readFileSync(resolve(site, ".gitignore"), "utf-8").split(/\r?\n/);
    expect(ignore).toContain(base.slice(1));
  });
});
