/**
 * The "Acknowledgements" chip and its popover, in the same `#hud` column as
 * the caption and "Report a bug or request a feature" link (index.html) — a
 * small, always-present toggle rather than a folder buried behind "advanced
 * controls" in the Tweakpane pane, so it costs nothing on screen when closed
 * and one click when it isn't.
 */

interface Credit {
  name: string;
  url: string;
}

const CREDITS: readonly Credit[] = [
  { name: "Peter van Keken", url: "https://carnegiescience.edu/bio/dr-peter-van-keken" },
  { name: "Cian Wilson", url: "https://carnegiescience.edu/bio/dr-cian-wilson" },
];

/**
 * One dot colour per row, cycling through the same validated categorical
 * palette the corner trace legends use (`nusselt.ts`'s inner/outer,
 * `rms.ts`'s two series) — reused rather than picked fresh so a panel with no
 * scientific content still reads as part of the same system as the ones that
 * do.
 */
const DOT_COLOURS = ["#d95926", "#3987e5", "#199e70", "#8b5cf6"];

/** Wires the static `#ack-toggle`/`#ack-panel` pair (index.html) together. */
export function buildAcknowledgements(toggle: HTMLElement, panel: HTMLElement): void {
  // No separate title row inside the panel — the toggle button's own label
  // ("Acknowledgements") already says what this is, and `#ack`'s squared-off
  // shared corner (index.html) is what ties the two together visually, so a
  // second copy of the word right below it would only repeat that. This
  // strapline instead says *why* the names below are here.
  const intro = document.createElement("p");
  intro.className = "ack-intro";
  intro.textContent = "Thanks to the following for their feedback and advice";
  panel.append(intro);

  const list = document.createElement("ul");
  CREDITS.forEach(({ name, url }, i) => {
    const li = document.createElement("li");
    const dot = document.createElement("span");
    dot.className = "ack-dot";
    dot.style.background = DOT_COLOURS[i % DOT_COLOURS.length];
    const a = document.createElement("a");
    a.href = url;
    a.textContent = name;
    a.target = "_blank";
    a.rel = "noopener";
    const arrow = document.createElement("span");
    arrow.className = "ack-arrow";
    arrow.textContent = "↗";
    arrow.setAttribute("aria-hidden", "true");
    li.append(dot, a, arrow);
    list.append(li);
  });

  panel.append(list);

  // Fades/grows in and out rather than snapping, like a real drop-down — but
  // still `display: none` at rest so it never sits in #hud's flex column
  // taking up a gap it isn't using. The two are sequenced across the
  // transition: opening sets `display: block` immediately, so there's
  // something painted for the very next frame's `.ack-open` to transition
  // from; closing removes `.ack-open` first and only reaches back for
  // `display: none` once `transitionend` confirms the fade actually finished.
  // The toggle button gets the same class at the same time — no transition
  // of its own, just squaring off its bottom corners so it and the panel
  // read as one continuous shape while open (index.html).
  let open = false;
  const show = (): void => {
    panel.style.display = "block";
    requestAnimationFrame(() => panel.classList.add("ack-open"));
    toggle.classList.add("ack-open");
    open = true;
  };
  const hide = (): void => {
    if (!open) return;
    open = false;
    panel.classList.remove("ack-open");
    toggle.classList.remove("ack-open");
    const onEnd = (e: TransitionEvent): void => {
      if (e.propertyName !== "opacity") return;
      panel.style.display = "none";
      panel.removeEventListener("transitionend", onEnd);
    };
    panel.addEventListener("transitionend", onEnd);
  };

  toggle.addEventListener("click", (e) => {
    e.stopPropagation();
    if (open) hide(); else show();
  });
  // A popover that only closes by re-clicking its own toggle reads as stuck
  // open the moment a reader clicks anywhere else to get back to the sim.
  document.addEventListener("click", (e) => {
    if (open && !panel.contains(e.target as Node)) hide();
  });
  document.addEventListener("keydown", (e) => {
    if (open && e.key === "Escape") hide();
  });
}
