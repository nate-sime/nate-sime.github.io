// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// Entry point: load the exported run, build the scene, and drive it from one
// clock. Time is a fractional output step t in [0, steps - 1]; the slab and the
// slices blend between the keyframes around t, the tracers between the steps
// around it, and the outline shows the nearest step.

import {
  Color,
  DirectionalLight,
  NeutralToneMapping,
  Object3D,
  PerspectiveCamera,
  PMREMGenerator,
  Scene,
  SRGBColorSpace,
  WebGLRenderer,
} from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { RoomEnvironment } from "three/addons/environments/RoomEnvironment.js";
import { Pane } from "tweakpane";
import { cssGradient, DEPTH, lutTexture, REDS } from "./colormap";
import { bracket, load } from "./data";
import { centreOf, toSceneInPlace, vec } from "./frame";
import { Outline } from "./outline";
import { DEPTH_SPAN, Slab, type SlabColour } from "./slab";
import { Slices } from "./slices";
import { Tracers } from "./tracers";

const $ = (id: string) => document.getElementById(id) as HTMLElement;
const msg = $("msg");

/** The published animation's frame, 844 × 539: narrower views back off to keep its width. */
const FIT_ASPECT = 1.5;
const MAX_TAIL_MYR = 3;

async function main(): Promise<void> {
  const canvas = $("view") as HTMLCanvasElement;
  let renderer: WebGLRenderer;
  try {
    renderer = new WebGLRenderer({ canvas, antialias: true });
  } catch {
    msg.textContent = "This browser could not start WebGL 2, which the viewer needs.";
    return;
  }
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.outputColorSpace = SRGBColorSpace;
  renderer.toneMapping = NeutralToneMapping;

  const data = await load(`${import.meta.env.BASE_URL}data/`, (f) => {
    msg.textContent = `Loading the model… ${Math.round(100 * f)}%`;
  });
  const m = data.manifest;
  const myrPerStep = m.myrPerStep;
  const totalMyr = (m.steps - 1) * myrPerStep;

  // Centre the scene on the domain as it is over the whole run.
  const centre = centreOf(data.outline);
  toSceneInPlace(centre, data.outline);
  toSceneInPlace(centre, data.tracers);

  // Redraw only when something changed; `lastT` is the step last applied.
  let dirty = true;
  let lastT = NaN;

  const scene = new Scene();
  scene.background = new Color(0xffffff);
  const pmrem = new PMREMGenerator(renderer);
  scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
  scene.environmentIntensity = 0.9;

  const camera = new PerspectiveCamera(m.camera.fov, 1, 5, 20000);
  // A key light that turns with the view, up and to the left of the eye.
  const key = new DirectionalLight(0xffffff, 1.6);
  key.position.set(-0.6, 1.0, 0.7);
  const keyTarget = new Object3D();
  keyTarget.position.set(0, 0, -1);
  key.target = keyTarget;
  camera.add(key, keyTarget);
  scene.add(camera);

  const tempLut = lutTexture(REDS);
  const depthLut = lutTexture(DEPTH);
  const slab = new Slab(data, centre, tempLut, depthLut);
  const slices = new Slices(data, centre, tempLut);
  const tracers = new Tracers(data, 4.5, Math.ceil(MAX_TAIL_MYR / myrPerStep));
  const outline = new Outline(data);
  scene.add(slices.group, tracers.spheres, tracers.tails, outline.lines, slab.mesh);

  // --- view --------------------------------------------------------------
  const controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.minDistance = 250;
  controls.maxDistance = 8000;
  const pvPos = vec(centre, m.camera.position);
  const pvDir = vec(centre, m.camera.target).sub(pvPos).normalize();
  const pvDist = -pvPos.dot(pvDir); // to the point on the line of sight nearest the domain's centre

  function resetView(): void {
    // The published view's direction and distance, but aimed at (and orbiting)
    // the middle of the domain rather than ParaView's focal point, which sat
    // well outside it; up is kept radial.
    const aspect = Math.max(camera.aspect, 0.2);
    // A little further back than ParaView's 2.6k km: aimed at the domain's
    // middle, its near corner would otherwise touch the frame's bottom edge.
    const dist = 1.12 * pvDist * Math.max(1, FIT_ASPECT / aspect);
    camera.up.set(0, 1, 0);
    camera.position.copy(pvDir).multiplyScalar(-dist);
    controls.target.set(0, 0, 0);
    camera.lookAt(controls.target);
    controls.update();
    dirty = true;
  }

  // Wheel-zoom only once the reader has clicked into the frame, so scrolling
  // the page past an embed doesn't get captured by it.
  const embedded = window.self !== window.top;
  const touch = window.matchMedia("(pointer: coarse)").matches;
  const hint = $("hint");
  controls.enableZoom = !embedded;
  canvas.addEventListener("pointerdown", () => {
    controls.enableZoom = true;
    hint.style.opacity = "0";
  });
  document.documentElement.addEventListener("pointerleave", () => {
    if (embedded) controls.enableZoom = false;
  });
  hint.textContent = touch
    ? "drag to rotate · pinch to zoom"
    : `drag to rotate · ${embedded ? "click, then scroll" : "scroll"} to zoom`;
  hint.title = touch ? "two fingers to pan" : "right-drag to pan";

  // --- state and controls ------------------------------------------------
  const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const P = {
    myr: reduceMotion ? totalMyr : 0,
    playing: !reduceMotion,
    speed: 1.5,
    colour: "titanium" as SlabColour,
    opacity: 0.85,
    contours: true,
    slices: m.slices.y.map(() => true),
    tracers: true,
    tail: 1.0,
    outline: true,
  };

  // Folded away in the page's frame, where it would cover the figure; open
  // when the viewer has a wide window to itself.
  const pane = new Pane({ title: "controls", container: $("pane"), expanded: !embedded && window.innerWidth >= 700 });
  // Two of each: the pane's play button and time slider, and the timeline's
  // along the bottom, which is always visible.
  const playButton = pane.addButton({ title: "play" });
  const playEl = $("play") as HTMLButtonElement;
  const scrub = $("scrub") as HTMLInputElement;
  scrub.max = String(totalMyr);
  const setPlaying = (on: boolean) => {
    P.playing = on;
    const label = on ? "pause" : "play";
    playButton.title = label;
    playEl.classList.toggle("playing", on);
    playEl.setAttribute("aria-label", label);
    playEl.title = label;
  };
  setPlaying(P.playing);
  const togglePlay = () => {
    if (!P.playing && P.myr >= totalMyr - 1e-6) P.myr = 0;
    setPlaying(!P.playing);
    timeBinding.refresh(); // during playback it only catches up every 100 ms
    wake(); // buttons emit "click", not the pane's "change"
  };
  playButton.on("click", togglePlay);
  playEl.addEventListener("click", togglePlay);

  const timeBinding = pane.addBinding(P, "myr", { label: "time (Myr)", min: 0, max: totalMyr, step: 0.01 });
  timeBinding.on("change", () => (dirty = true));
  // Grabbing either slider pauses, so playback doesn't snatch it back.
  timeBinding.element.addEventListener("pointerdown", () => setPlaying(false), true);
  const seek = (myr: number) => {
    setPlaying(false);
    P.myr = Math.min(totalMyr, Math.max(0, myr));
    scrub.value = String(P.myr);
    timeBinding.refresh();
    dirty = true;
    wake();
  };
  scrub.addEventListener("pointerdown", () => setPlaying(false));
  scrub.addEventListener("input", () => seek(Number(scrub.value)));
  // Keys step by a tenth of a Myr (a million years with Page Up/Down), which
  // a range with step="any" — kept so the thumb glides during playback —
  // would otherwise leave to the browser.
  scrub.addEventListener("keydown", (e) => {
    const steps: Record<string, number> = { ArrowLeft: -0.1, ArrowDown: -0.1, ArrowRight: 0.1, ArrowUp: 0.1, PageDown: -1, PageUp: 1 };
    let to: number;
    if (e.key in steps) to = Math.round((P.myr + steps[e.key]) * 10) / 10;
    else if (e.key === "Home") to = 0;
    else if (e.key === "End") to = totalMyr;
    else return;
    e.preventDefault();
    seek(to);
  });
  pane.addBinding(P, "speed", { label: "speed (Myr/s)", min: 0.25, max: 4, step: 0.05 });

  const fSlab = pane.addFolder({ title: "slab top" });
  fSlab
    .addBinding(P, "colour", {
      label: "colour",
      options: { titanium: "titanium", "slab-top temperature": "temperature", depth: "depth" },
    })
    .on("change", (e) => {
      slab.setColour(e.value);
      $("legend-depth").hidden = e.value !== "depth";
      dirty = true;
    });
  fSlab.addBinding(P, "opacity", { min: 0.2, max: 1, step: 0.05 }).on("change", (e) => {
    slab.setOpacity(e.value);
    dirty = true;
  });
  fSlab.addBinding(P, "contours", { label: "depth contours" }).on("change", (e) => {
    slab.setContours(e.value);
    dirty = true;
  });

  const fSlices = pane.addFolder({ title: "temperature slices" });
  m.slices.y.forEach((y, i) => {
    const holder = { on: true };
    fSlices.addBinding(holder, "on", { label: `y = ${y > 0 ? "+" : ""}${y} km` }).on("change", (e) => {
      slices.meshes[i].visible = e.value;
      dirty = true;
    });
  });

  const fTracers = pane.addFolder({ title: "tracers" });
  fTracers.addBinding(P, "tracers", { label: "show" }).on("change", (e) => {
    tracers.setVisible(e.value);
    dirty = true;
  });
  fTracers.addBinding(P, "tail", { label: "tail (Myr)", min: 0, max: MAX_TAIL_MYR, step: 0.1 }).on("change", () => {
    lastT = NaN;
    dirty = true;
  });

  const fView = pane.addFolder({ title: "view", expanded: false });
  fView.addBinding(P, "outline", { label: "domain outline" }).on("change", (e) => {
    outline.lines.visible = e.value;
    dirty = true;
  });
  fView.addButton({ title: "reset view" }).on("click", () => {
    resetView();
    wake();
  });

  slab.setOpacity(P.opacity);

  // --- HUD ---------------------------------------------------------------
  const timeEl = $("time");
  $("about").textContent = "Sime et al. (2023), fig. 16";
  $("hud").title = `subduction-3d v${__APP_VERSION__}`;
  const [tMin, tMax] = m.temperatureRange;
  $("bar-temp").style.background = cssGradient(REDS);
  $("ticks-temp").innerHTML = [300, 700, 1100, 1500]
    .map((T) => `<span style="left:${(100 * (T - tMin)) / (tMax - tMin)}%">${T}</span>`)
    .join("");
  $("bar-depth").style.background = cssGradient(DEPTH);
  $("ticks-depth").innerHTML = [0, 50, 100, 150, 200]
    .map((d) => `<span style="left:${(100 * d) / DEPTH_SPAN}%">${d}</span>`)
    .join("");

  // --- clock -------------------------------------------------------------
  let holdUntil = 0;
  let lastFrame = NaN;
  let lastPaneRefresh = 0;

  function apply(): void {
    const t = Math.min(m.steps - 1, Math.max(0, P.myr / myrPerStep));
    if (t === lastT) return;
    lastT = t;
    const { k, s } = bracket(m.keyframes, t);
    slab.update(k, s);
    slices.update(k, s);
    tracers.update(t, P.tail / myrPerStep);
    outline.update(t);
    const label = `${P.myr.toFixed(1)} Myr`;
    timeEl.textContent = label;
    scrub.value = String(P.myr);
    scrub.setAttribute("aria-valuetext", label);
  }

  function resize(): void {
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    tracers.setResolution(w, h);
    outline.setResolution(w, h);
    dirty = true;
  }

  // At most one loop. `raf` stays non-zero from scheduling a frame until it
  // has run, including while it runs: a wake() from inside a frame (the time
  // slider's refresh, the controls' damping) must not start a second loop.
  const RUNNING = -1;
  let raf = 0;
  let onScreen = true;
  function frame(now: number): void {
    raf = RUNNING;
    // Time only from frame timestamps: in an iframe they need not share an
    // origin with performance.now(). The first frame after a wake steps 0.
    const dt = Number.isNaN(lastFrame) ? 0 : Math.min(0.1, Math.max(0, (now - lastFrame) / 1000));
    lastFrame = now;
    if (P.playing && now >= holdUntil) {
      if (P.myr >= totalMyr) P.myr = 0;
      P.myr = Math.min(totalMyr, P.myr + dt * P.speed);
      if (P.myr >= totalMyr) holdUntil = now + 1500; // rest on the last frame before looping
      dirty = true;
      if (now - lastPaneRefresh > 100) {
        timeBinding.refresh();
        lastPaneRefresh = now;
      }
    }
    if (dirty) apply();
    const moved = controls.update();
    if (moved || dirty) {
      renderer.render(scene, camera);
      dirty = false;
    }
    // Run while playing or while the view is still moving (damping); otherwise
    // sleep until an input calls wake().
    raf = (P.playing || moved) && onScreen && !document.hidden ? requestAnimationFrame(frame) : 0;
  }
  const wake = () => {
    if (!raf && onScreen && !document.hidden) {
      lastFrame = NaN;
      raf = requestAnimationFrame(frame);
    }
  };
  controls.addEventListener("change", () => {
    dirty = true;
    wake();
  });
  pane.on("change", wake);
  // Stop drawing while the frame is scrolled out of view or the tab is hidden.
  new IntersectionObserver((entries) => {
    onScreen = entries.some((e) => e.isIntersecting);
    wake();
  }).observe(canvas);
  document.addEventListener("visibilitychange", wake);

  new ResizeObserver(() => {
    resize();
    wake();
  }).observe(canvas);
  resize();
  resetView();
  apply();
  msg.hidden = true;
  $("timeline").hidden = false;
  wake();
}

main().catch((e: unknown) => {
  msg.hidden = false;
  msg.textContent = `Could not load the model: ${e instanceof Error ? e.message : String(e)}`;
  console.error(e);
});
