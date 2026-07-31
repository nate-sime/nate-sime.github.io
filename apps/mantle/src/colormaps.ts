/**
 * Colour maps for the temperature field: five RGB control points evenly
 * spaced over [0, 1], shared between the WGSL fragment shader (`gpu/wgsl.ts`,
 * which embeds them as shader constants) and the pane's colour-bar legend
 * (`ui/colorbar.ts`, which reads the same numbers into a CSS gradient) — one
 * table, so the bar drawn beside the field can never disagree with it.
 *
 * `inferno` is the map the app has always drawn in: perceptually uniform and
 * monotone in lightness, so the field reads correctly in greyscale and stays
 * legible with colour-vision deficiency. The rest trade that guarantee for
 * familiarity (`coolwarm`) or hue range (`turbo`) — both still legible, just
 * not by lightness alone.
 */

export type Stop = readonly [number, number, number];

export const COLORMAPS = {
  inferno: [
    [0.0, 0.0, 0.016], [0.341, 0.063, 0.431], [0.737, 0.216, 0.329],
    [0.976, 0.557, 0.035], [0.988, 1.0, 0.643],
  ],
  viridis: [
    [0.267, 0.005, 0.329], [0.229, 0.322, 0.545], [0.128, 0.567, 0.551],
    [0.369, 0.789, 0.383], [0.993, 0.906, 0.144],
  ],
  plasma: [
    [0.050, 0.030, 0.528], [0.494, 0.012, 0.658], [0.798, 0.280, 0.469],
    [0.973, 0.585, 0.254], [0.940, 0.975, 0.131],
  ],
  magma: [
    [0.001, 0.000, 0.014], [0.317, 0.072, 0.485], [0.716, 0.215, 0.475],
    [0.972, 0.469, 0.361], [0.987, 0.991, 0.749],
  ],
  coolwarm: [
    [0.230, 0.299, 0.754], [0.552, 0.690, 0.996], [0.865, 0.865, 0.865],
    [0.955, 0.639, 0.505], [0.706, 0.016, 0.150],
  ],
  turbo: [
    [0.190, 0.072, 0.232], [0.136, 0.592, 0.921], [0.610, 0.908, 0.294],
    [0.990, 0.564, 0.133], [0.480, 0.016, 0.011],
  ],
} as const satisfies Record<string, readonly Stop[]>;

export type ColormapName = keyof typeof COLORMAPS;
