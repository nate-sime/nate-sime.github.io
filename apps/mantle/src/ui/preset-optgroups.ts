/**
 * Cosmetic-only: groups the flat "try an example" dropdown into visual
 * submenus via native `<optgroup>`, since Tweakpane's list view always
 * renders a flat `<select>`/`<option>` list with no grouping of its own —
 * checked directly against the installed version: `ListParamsOptions`'s own
 * type is a flat array or a flat `{[text]: value}` object (`@tweakpane/core`'s
 * `common/params.d.ts`), and `ListView`'s DOM builder (`common/view/list.js`)
 * appends every option in one flat loop, nothing nested anywhere.
 *
 * Deliberately isolated in its own file rather than inlined into
 * `controls.ts`: this reaches into Tweakpane's *rendered* DOM the same way
 * the `dtMax`/`dtInitial`/`courant` log-sliders elsewhere in that file
 * already do, but unlike those (which add a capability Tweakpane's own
 * control cannot express at all), this is purely a visual regrouping of a
 * control that already works correctly without it — the kind of thing
 * Tweakpane could grow native support for later, or that could simply be
 * dropped if it turns out not to be worth the DOM reach-in. One file, one
 * function, called from one place: removing this feature is deleting the
 * `applyOptgroups(...)` call site and this file, nothing else.
 *
 * Applied once, right after the binding is built. The preset dropdown's own
 * options are a static table (`QUICK_STARTS`/`BENCHMARKS` merged in
 * `controls.ts`), set once at pane construction and never reassigned after —
 * so there is no later rebuild this would need to be re-applied against.
 */

/** Option value → the `<optgroup>` label it belongs under; absent stays top-level. */
export type GroupMap = Record<string, string>;

/**
 * Wrap runs of consecutive `<option>`s whose value is in `groups` into
 * `<optgroup>` elements, preserving the select's original order — an
 * ungrouped option (or a grouped one not adjacent to its own kind) is left
 * exactly where it was. Safe to call on an already-grouped select: reading
 * `select.options` flattens across existing `<optgroup>`s, so a second call
 * re-derives the same grouping rather than nesting `<optgroup>` in
 * `<optgroup>`.
 */
export function applyOptgroups(select: HTMLSelectElement, groups: GroupMap): void {
  const options = Array.from(select.options);
  let current: HTMLOptGroupElement | null = null;
  let currentLabel: string | null = null;
  for (const opt of options) {
    const label = groups[opt.value] ?? null;
    if (label === null) {
      select.appendChild(opt);
      current = null;
      currentLabel = null;
      continue;
    }
    if (label !== currentLabel || !current) {
      current = select.ownerDocument.createElement("optgroup");
      current.label = label;
      select.appendChild(current);
      currentLabel = label;
    }
    current.appendChild(opt);
  }
}

/**
 * Derive a `GroupMap` from a "Family case" naming convention — "Blankenbach
 * 1a" groups under "Blankenbach", "Tosi 4" under "Tosi", "van Keken 1a" under
 * "van Keken" — rather than hand-listing every current (and future) case.
 * Adding "van Keken 1b" to `BENCHMARKS` groups itself with no change here.
 *
 * A name with nothing to split on (no trailing whitespace-separated token) is
 * left out of the map, hence ungrouped — there is no case to guess a family
 * from. Intended for `Object.keys(BENCHMARKS)` specifically, not the quick
 * starts: their friendly names ("Sluggish mantle") would false-positive
 * against this same pattern.
 */
export function deriveGroups(names: readonly string[]): GroupMap {
  const groups: GroupMap = {};
  for (const name of names) {
    const m = /^(.+)\s+\S+$/.exec(name);
    if (m) groups[name] = m[1];
  }
  return groups;
}
