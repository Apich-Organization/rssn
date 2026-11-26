//! Scientific colormaps for scalar field visualization.

/// Perceptual and diverging colormaps for heatmaps, PNGs, and animated GIFs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Colormap {
    /// Viridis (perceptually uniform, dark purple -> teal -> yellow).
    #[default]
    Viridis,
    /// Plasma (perceptually uniform, dark blue -> magenta -> yellow).
    Plasma,
    /// Inferno (perceptually uniform, black -> red -> yellow).
    Inferno,
    /// Magma (perceptually uniform, black -> purple -> cream).
    Magma,
    /// Turbo (high-contrast rainbow).
    Turbo,
    /// Coolwarm (diverging blue -> white -> red).
    Coolwarm,
    /// Linear grayscale (black -> white).
    Grayscale,
}

impl Colormap {
    /// Samples the colormap at `t` in `[0.0, 1.0]` and returns `(r, g, b)`.
    #[must_use]
    pub fn sample_rgb(
        self,
        t: f64,
    ) -> (u8, u8, u8) {
        let t = if t.is_nan() {
            0.0
        } else {
            t.clamp(0.0, 1.0)
        };
        match self {
            | Self::Grayscale => {
                let v = (t * 255.0).round() as u8;
                (v, v, v)
            },
            | Self::Viridis => {
                interp_stops(
                    t,
                    &[
                        (68, 1, 84),
                        (59, 82, 139),
                        (33, 145, 140),
                        (94, 201, 98),
                        (253, 231, 37),
                    ],
                )
            },
            | Self::Plasma => {
                interp_stops(
                    t,
                    &[
                        (13, 8, 135),
                        (126, 3, 168),
                        (204, 71, 120),
                        (248, 149, 64),
                        (240, 249, 33),
                    ],
                )
            },
            | Self::Inferno => {
                interp_stops(
                    t,
                    &[
                        (0, 0, 4),
                        (87, 16, 110),
                        (188, 55, 84),
                        (249, 142, 9),
                        (252, 255, 164),
                    ],
                )
            },
            | Self::Magma => {
                interp_stops(
                    t,
                    &[
                        (0, 0, 4),
                        (81, 18, 124),
                        (183, 55, 121),
                        (252, 137, 97),
                        (252, 253, 191),
                    ],
                )
            },
            | Self::Turbo => {
                interp_stops(
                    t,
                    &[
                        (48, 18, 59),
                        (40, 188, 235),
                        (164, 252, 60),
                        (251, 126, 33),
                        (122, 4, 3),
                    ],
                )
            },
            | Self::Coolwarm => {
                interp_stops(
                    t,
                    &[
                        (59, 76, 192),
                        (141, 176, 254),
                        (221, 221, 221),
                        (244, 154, 123),
                        (180, 4, 38),
                    ],
                )
            },
        }
    }

    /// Samples the colormap as a `#rrggbb` hex string.
    #[must_use]
    pub fn sample_hex(
        self,
        t: f64,
    ) -> String {
        let (r, g, b) = self.sample_rgb(t);
        format!("#{r:02x}{g:02x}{b:02x}")
    }
}

fn interp_stops(
    t: f64,
    stops: &[(u8, u8, u8)],
) -> (u8, u8, u8) {
    let n = stops.len() - 1;
    let scaled = t * (n as f64);
    let idx = (scaled.floor() as usize).min(n - 1);
    let frac = scaled - (idx as f64);
    let (r0, g0, b0) = stops[idx];
    let (r1, g1, b1) = stops[idx + 1];
    let lerp = |a: u8, b: u8| -> u8 {
        (f64::from(a) + (f64::from(b) - f64::from(a)) * frac)
            .round()
            .clamp(0.0, 255.0) as u8
    };
    (lerp(r0, r1), lerp(g0, g1), lerp(b0, b1))
}
