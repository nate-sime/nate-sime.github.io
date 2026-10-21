// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The checkpoint library: the pane's "checkpoints" folder, and the dialog
 * behind it where tags are named, kept in a file and started from.
 *
 * The library lives in memory; the file is what it is opened from and saved
 * to. Where the browser has the File System Access API (Chromium, which is
 * also where WebGPU is most at home) "open…" keeps the file's handle and
 * "save" writes back to it in place. Elsewhere "open…" is a file input and
 * "save" a download, and "save as…" — which would be the same download — is
 * not offered. Either way an unsaved library says so, in the pane and in the
 * dialog, and leaving the page asks first.
 *
 * Tagging captures at the click, then opens the dialog with the new name
 * selected — so the tag is the field that was on screen when the button was
 * pressed, not whatever the run had become by the time a name was typed.
 *
 * Every number shown, every rule about what may be reseeded or restored, and
 * the file format itself are in `checkpoint.ts`; this file is DOM only. The
 * solver is reached through `CheckpointHooks`, which `main.ts` provides.
 */

import type { FolderApi } from "tweakpane";
import { COLORMAPS, type Stop } from "../colormaps";
import {
  CHECKPOINT_EXTENSION, CheckpointError, decodeCheckpoints, defaultName, describe,
  domainMismatch, encodeCheckpoints, restorePatch, thumbnail, uniqueName,
  type Checkpoint, type Domain,
} from "../checkpoint";
import type { State } from "./presets";

export interface CheckpointHooks {
  /** Read the running state — everything a tag holds but its name. Rejects when there is no run. */
  capture(): Promise<Omit<Checkpoint, "name">>;
  /** Start the run over from a tag's field, under the pane's current settings. */
  reseed(c: Checkpoint): Promise<void>;
  /** Put a tag's settings back, then continue its run from its own clock. */
  restore(c: Checkpoint): Promise<void>;
  /** The live run's domain — null while there is none, mid-rebuild. */
  domain(): Domain | null;
  /** Why the run cannot be touched just now, or null when it can. */
  busy(): string | null;
  state(): State;
}

export interface CheckpointsHandle {
  /**
   * Re-read what each tag can do and how it is drawn — after a rebuild has
   * changed the domain, or the colour map has changed the thumbnails.
   */
  refresh(): void;
}

// The File System Access API, declared here rather than relied on from the
// DOM library, which does not carry the pickers — and only as much of it as
// is used.
interface WritableFile { write(data: Blob): Promise<void>; close(): Promise<void> }
interface FileHandle { readonly name: string; getFile(): Promise<File>; createWritable(): Promise<WritableFile> }
interface PickerOptions {
  id?: string;
  suggestedName?: string;
  types: { description: string; accept: Record<string, string[]> }[];
}
interface FilePickers {
  showOpenFilePicker?(o: PickerOptions): Promise<FileHandle[]>;
  showSaveFilePicker?(o: PickerOptions): Promise<FileHandle>;
}

const PICKER: PickerOptions = {
  // Remembers the last folder used for these files, apart from other sites' pickers.
  id: "mantle-checkpoints",
  types: [{
    description: "Mantle checkpoints",
    accept: { "application/octet-stream": [CHECKPOINT_EXTENSION] },
  }],
};
const DEFAULT_FILE = `checkpoints${CHECKPOINT_EXTENSION}`;
/** Thumbnail raster size; drawn at half that in CSS pixels, so it stays sharp at 2× DPR. */
const THUMB = 96;
/** How long a "delete?" stays armed for its second click. */
const ARM_MS = 3000;

const isAbort = (e: unknown): boolean => e instanceof DOMException && e.name === "AbortError";
const messageOf = (e: unknown): string => (e instanceof Error ? e.message : String(e));
const plural = (n: number, word: string): string => `${n} ${word}${n === 1 ? "" : "s"}`;

function el<K extends keyof HTMLElementTagNameMap>(
  tag: K, cls?: string, text?: string,
): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function button(text: string, onClick: () => void): HTMLButtonElement {
  const b = el("button", undefined, text);
  b.type = "button";
  b.addEventListener("click", onClick);
  return b;
}

interface Row {
  name: HTMLInputElement;
  reseed: HTMLButtonElement;
  restore: HTMLButtonElement;
  del: HTMLButtonElement;
}

export function buildCheckpoints(folder: FolderApi, hooks: CheckpointHooks): CheckpointsHandle {
  const pickers = window as unknown as FilePickers;
  const inPlace = typeof pickers.showOpenFilePicker === "function"
    && typeof pickers.showSaveFilePicker === "function";

  let entries: Checkpoint[] = [];
  let fileName = DEFAULT_FILE;
  /** The open file, when the browser can write back to it. */
  let handle: FileHandle | null = null;
  let dirty = false;
  /** An action is in flight; every button that starts another stands down. */
  let working = false;
  const names = (except?: Checkpoint): string[] =>
    entries.filter((c) => c !== except).map((c) => c.name);

  // ---- the folder ----------------------------------------------------------

  const status = { library: "" };
  folder.addButton({ title: "tag state" }).on("click", () => void tag());
  folder.addButton({ title: "manage checkpoints…" }).on("click", () => openDialog());
  // A monitor refreshed by hand (`interval: 0`): it only changes when this
  // file changes it, so there is nothing for a timer to poll.
  const libraryLine = folder.addBinding(status, "library", { readonly: true, interval: 0 });

  // ---- the dialog ------------------------------------------------------------

  const dialog = el("dialog", "pane-dialog checkpoint-dialog");
  dialog.setAttribute("aria-labelledby", "checkpoint-title");
  const title = el("h2", undefined, "Checkpoints");
  title.id = "checkpoint-title";
  const intro = el("p", undefined,
    "Tag the running temperature field to start from it later. “reseed” starts the run over " +
    "from a tag under the current settings; “restore” puts the tag's settings back too, and " +
    "carries on from where it was tagged. Tracers are not saved.");
  const fileInfo = el("p", "checkpoint-file");
  const fileActions = el("div", "checkpoint-file-actions");
  const newButton = button("new", () => void startNew());
  const openButton = button("open…", () => void open(false));
  const importButton = button("import…", () => void open(true));
  const saveButton = button("save", () => void save(false));
  const saveAsButton = button("save as…", () => void save(true));
  openButton.title = "Replace this library with one from a file.";
  importButton.title = "Add the tags from another file to this library.";
  saveButton.title = inPlace ? "Write this library to its file." : "Download this library as a file.";
  saveAsButton.hidden = !inPlace;
  fileActions.append(newButton, openButton, importButton, saveButton, saveAsButton);

  const confirmBar = el("div", "checkpoint-confirm");
  confirmBar.hidden = true;
  const confirmText = el("span");
  const discardButton = button("discard", () => settle?.(true));
  const keepButton = button("keep", () => settle?.(false));
  confirmBar.append(confirmText, discardButton, keepButton);
  /** Resolves the pending confirmation, if there is one. */
  let settle: ((discard: boolean) => void) | null = null;

  const list = el("ol", "checkpoint-list");
  const empty = el("p", "checkpoint-empty", "No tags yet — “tag state” adds the field on screen.");
  const statusLine = el("p", "checkpoint-status");
  statusLine.setAttribute("role", "status");
  const actions = el("div", "dialog-actions");
  const tagButton = button("tag state", () => void tag());
  const closeButton = button("close", () => dialog.close());
  actions.append(tagButton, closeButton);

  // Only used where there is no picker API to open a file with.
  const fileInput = el("input");
  fileInput.type = "file";
  fileInput.accept = CHECKPOINT_EXTENSION;
  fileInput.hidden = true;

  dialog.append(title, intro, fileInfo, fileActions, confirmBar, list, empty, statusLine,
    actions, fileInput);
  document.body.append(dialog);
  // Non-modal, as "create a planet" is — the run stays visible, and tagging
  // it again stays a click away — so Escape is not wired up by the browser.
  dialog.addEventListener("keydown", (e) => { if (e.key === "Escape") dialog.close(); });
  dialog.addEventListener("close", () => settle?.(false));

  const say = (msg: string, error = false): void => {
    statusLine.textContent = msg;
    statusLine.classList.toggle("error", error);
  };

  // ---- rows ------------------------------------------------------------------

  const rows = new Map<Checkpoint, Row>();
  const thumbs = new WeakMap<Checkpoint, { stops: readonly Stop[]; px: Uint8ClampedArray<ArrayBuffer> }>();
  const thumbOf = (c: Checkpoint, stops: readonly Stop[]): Uint8ClampedArray<ArrayBuffer> => {
    const held = thumbs.get(c);
    if (held?.stops === stops) return held.px;
    const px = thumbnail(c, THUMB, stops);
    thumbs.set(c, { stops, px });
    return px;
  };

  const rename = (c: Checkpoint, input: HTMLInputElement): void => {
    const wanted = input.value.trim();
    if (wanted && wanted !== c.name) {
      c.name = uniqueName(wanted, names(c));
      dirty = true;
      sync();
    }
    input.value = c.name;
  };

  const remove = (c: Checkpoint): void => {
    entries = entries.filter((e) => e !== c);
    dirty = true;
    render();
    say(`deleted “${c.name}”`);
  };

  /** Rebuild the list. `sync` is the cheap half, for when only states change. */
  const render = (): void => {
    rows.clear();
    list.replaceChildren();
    const stops = COLORMAPS[hooks.state().colormap];
    for (const c of entries) {
      const li = el("li", "checkpoint-row");
      const canvas = el("canvas", "checkpoint-thumb");
      canvas.width = canvas.height = THUMB;
      canvas.getContext("2d")?.putImageData(new ImageData(thumbOf(c, stops), THUMB, THUMB), 0, 0);
      const body = el("div", "checkpoint-body");
      const name = el("input", "checkpoint-name");
      name.type = "text";
      name.value = c.name;
      name.maxLength = 120;
      name.spellcheck = false;
      name.setAttribute("aria-label", "checkpoint name");
      // "change" is the commit — Enter or leaving the field — so a name being
      // typed is never half-applied.
      name.addEventListener("change", () => rename(c, name));
      name.addEventListener("keydown", (e) => {
        if (e.key === "Enter") name.blur();
        if (e.key === "Escape") {
          name.value = c.name;
          name.blur();
          e.stopPropagation();   // revert the name, not close the dialog
        }
      });
      const summary = el("p", "checkpoint-summary", describe(c));
      const when = c.created ? new Date(c.created).toLocaleString() : "an unknown time";
      summary.title = `Tagged ${when}${c.appVersion ? ` by v${c.appVersion}` : ""}.`;
      const reseed = button("reseed", () => void act(`reseeded from “${c.name}”`, () => hooks.reseed(c)));
      const restore = button("restore", () => void act(`restored “${c.name}”`, () => hooks.restore(c)));
      // Two clicks: a tag can be hours of run, and this one is not undone by
      // anything short of reopening a saved file.
      const del = button("delete", () => {
        if (del.dataset.armed) return remove(c);
        del.dataset.armed = "yes";
        del.textContent = "delete?";
        window.setTimeout(() => {
          delete del.dataset.armed;
          del.textContent = "delete";
        }, ARM_MS);
      });
      const buttons = el("div", "checkpoint-actions");
      buttons.append(reseed, restore, del);
      body.append(name, summary, buttons);
      li.append(canvas, body);
      list.append(li);
      rows.set(c, { name, reseed, restore, del });
    }
    sync();
  };

  /** Every enabled/disabled state, reason and count — no DOM rebuilt, so focus stays put. */
  const sync = (): void => {
    const domain = hooks.domain();
    const state = hooks.state();
    const blocked = working || hooks.busy() !== null;
    for (const [c, row] of rows) {
      const why = domain ? domainMismatch(c.domain, domain) : null;
      row.reseed.disabled = blocked || !domain || why !== null;
      row.reseed.title = !domain ? "No run is live just now."
        : why ? `This tag is ${why}. “restore” brings its settings with it.`
          : "Start the run over from this field, under the current settings.";
      const r = restorePatch(c, state);
      row.restore.disabled = blocked || "problem" in r;
      row.restore.title = "problem" in r ? `Can't restore this tag: ${r.problem}.`
        : "Put this tag's settings back and carry on from where it was tagged.";
      row.del.disabled = working;
    }
    for (const b of [tagButton, newButton, openButton, importButton, saveButton, saveAsButton])
      b.disabled = working;
    const summary = `${plural(entries.length, "tag")}${dirty ? ", unsaved" : ""}`;
    fileInfo.textContent = `${fileName} — ${summary}`;
    // The count and the unsaved marker first: the pane's field is narrow, and
    // they are what it is there to say.
    status.library = `${summary} · ${fileName}`;
    libraryLine.refresh();
    empty.hidden = entries.length > 0;
  };

  const openDialog = (): void => {
    if (!dialog.open) dialog.show();
    render();
  };

  // ---- actions ---------------------------------------------------------------

  /** Run one solver-touching action at a time, and say how it went. */
  const act = async (done: string, job: () => Promise<void>): Promise<void> => {
    if (working) return;
    const why = hooks.busy();
    if (why) return say(why, true);
    working = true;
    sync();
    say("working…");
    try {
      await job();
      say(done);
    } catch (e) {
      say(messageOf(e), true);
    } finally {
      working = false;
      sync();
    }
  };

  const tag = async (): Promise<void> => {
    if (working) return;
    const why = hooks.busy();
    if (why) {
      openDialog();
      return say(why, true);
    }
    working = true;
    try {
      const c = await hooks.capture();
      const entry: Checkpoint = { ...c, name: defaultName(c.state, c.domain, c.clock.time, names()) };
      entries.push(entry);
      dirty = true;
      openDialog();
      const input = rows.get(entry)?.name;
      input?.scrollIntoView({ block: "nearest" });
      input?.focus();
      input?.select();
      say(`tagged “${entry.name}” — type to rename it, or keep the default`);
    } catch (e) {
      openDialog();
      say(`Couldn't tag the state: ${messageOf(e)}`, true);
    } finally {
      working = false;
      sync();
    }
  };

  /** Ask before an action throws away unsaved changes; true to go ahead. */
  const confirmDiscard = (doing: string): Promise<boolean> => {
    if (!dirty) return Promise.resolve(true);
    settle?.(false);
    return new Promise((resolve) => {
      confirmText.textContent = `This library has unsaved changes. Discard them and ${doing}?`;
      confirmBar.hidden = false;
      settle = (discard) => {
        settle = null;
        confirmBar.hidden = true;
        resolve(discard);
      };
    });
  };

  const startNew = async (): Promise<void> => {
    if (!(await confirmDiscard("start a new library"))) return;
    entries = [];
    fileName = DEFAULT_FILE;
    handle = null;
    dirty = false;
    render();
    say("new library");
  };

  /** The fallback picker. Resolves null if the browser reports a cancel; may never settle if it cannot. */
  const pickFile = (): Promise<File | null> => new Promise((resolve) => {
    fileInput.value = "";
    fileInput.onchange = () => resolve(fileInput.files?.[0] ?? null);
    fileInput.oncancel = () => resolve(null);
    fileInput.click();
  });

  /** "open…" replaces the library; "import…" (`merge`) adds a file's tags to it. */
  const open = async (merge: boolean): Promise<void> => {
    if (!merge && !(await confirmDiscard("open another"))) return;
    let file: File | null;
    let picked: FileHandle | null = null;
    try {
      if (pickers.showOpenFilePicker) {
        [picked] = await pickers.showOpenFilePicker(PICKER);
        file = await picked.getFile();
      } else file = await pickFile();
    } catch (e) {
      if (!isAbort(e)) say(`Couldn't open the file: ${messageOf(e)}`, true);
      return;
    }
    if (!file) return;
    let read: Checkpoint[];
    try {
      read = decodeCheckpoints(await file.arrayBuffer());
    } catch (e) {
      say(`Couldn't open ${file.name}: ${e instanceof CheckpointError ? e.message : messageOf(e)}.`, true);
      return;
    }
    if (!merge) {
      entries = [];
      fileName = file.name;
      handle = inPlace ? picked : null;
    }
    // A hand-edited file can repeat a name; the library cannot.
    for (const c of read) {
      c.name = uniqueName(c.name, names());
      entries.push(c);
    }
    dirty = merge ? dirty || read.length > 0 : false;
    render();
    say(merge ? `imported ${plural(read.length, "tag")} from ${file.name}`
      : `opened ${file.name}: ${plural(read.length, "tag")}`);
  };

  const download = (blob: Blob): void => {
    const a = el("a");
    a.href = URL.createObjectURL(blob);
    a.download = fileName;
    a.click();
    // Revoked late rather than at once: the download starts asynchronously.
    window.setTimeout(() => URL.revokeObjectURL(a.href), 60_000);
  };

  /** "save" writes back to the open file where it can; "save as…" (`as`) always asks where. */
  const save = async (as: boolean): Promise<void> => {
    if (working) return;
    // Encoded before any await, so the picker still runs inside the click's
    // user activation, which a browser requires of it.
    const blob = new Blob([encodeCheckpoints(entries)], { type: "application/octet-stream" });
    try {
      if (inPlace) {
        if (as || !handle) {
          handle = await pickers.showSaveFilePicker!({ ...PICKER, suggestedName: fileName });
          fileName = handle.name;
        }
        const out = await handle.createWritable();
        await out.write(blob);
        await out.close();
      } else download(blob);
    } catch (e) {
      if (!isAbort(e)) say(`Couldn't save: ${messageOf(e)}`, true);
      return;
    }
    dirty = false;
    sync();
    say(`${inPlace ? "saved" : "downloaded"} ${plural(entries.length, "tag")} to ${fileName}`);
  };

  // Hours of run in an unsaved tag are worth one question before the page goes.
  window.addEventListener("beforeunload", (e) => {
    if (!dirty) return;
    e.preventDefault();
    e.returnValue = "";
  });

  sync();
  return { refresh: () => { if (dialog.open) render(); else sync(); } };
}
