/**
 * The "Acknowledgements" chip and its popover, in the same `#hud` column as
 * the caption and "Report a bug" link (index.html) — a small, always-present
 * toggle rather than a folder buried behind "advanced controls" in the
 * Tweakpane pane, so it costs nothing on screen when closed and one click
 * when it isn't.
 */

interface Credit {
  name: string;
  url: string;
}

const CREDITS: readonly Credit[] = [
  { name: "Peter van Keken", url: "https://carnegiescience.edu/bio/dr-peter-van-keken" },
  { name: "Cian Wilson", url: "https://carnegiescience.edu/bio/dr-cian-wilson" },
  { name: "Tim Jones", url: "https://www.fleetspace.com/team/timothy-jones" },
];

/** Wires the static `#ack-toggle`/`#ack-panel` pair (index.html) together. */
export function buildAcknowledgements(toggle: HTMLElement, panel: HTMLElement): void {
  const list = document.createElement("ul");
  for (const { name, url } of CREDITS) {
    const li = document.createElement("li");
    const a = document.createElement("a");
    a.href = url;
    a.textContent = name;
    a.target = "_blank";
    a.rel = "noopener";
    li.append(a);
    list.append(li);
  }
  panel.append(list);

  toggle.addEventListener("click", () => {
    if (panel.hasAttribute("data-show")) panel.removeAttribute("data-show");
    else panel.setAttribute("data-show", "");
  });
}
