//! Standalone SVG plots.
//!
//! Every function writes one self-contained `.svg` file and needs no
//! plotting dependency. Functions are given as plain Rust closures or data;
//! [`plot_term_2d`] plots a term of the expression graph.
//!
//! 3D plots use one fixed oblique projection (yaw 0.6, pitch 0.5 radians):
//! surfaces are painted back to front as height-coloured quads, curves and
//! vector fields as projected lines inside a wireframe bounding box.
//!
//! # Limitations
//!
//! - **Output formats.** Most functions write self-contained SVG; the
//!   [`plot_scalar_field_pgm`] / [`plot_scalar_field_ppm`] pair writes 8-bit
//!   PGM/PPM rasters through the pure-Rust [`super::raster`] writer. There is
//!   no PNG, JPEG or other compressed raster export: that needs an image
//!   encoder, i.e. a new dependency.
//! - **Plot types.** 2-D: function, series, expression term, vector field,
//!   scaled quiver ([`plot_quiver_2d`]), contour lines ([`plot_contour_2d`]),
//!   filled contour bands ([`plot_contour_bands_2d`]) and streamlines
//!   ([`plot_streamlines_2d`]). 3-D: surface, parametric curve, polyline and
//!   vector field, all with one fixed oblique projection (no interactive
//!   rotation or real depth cueing). Missing entirely: log axes, subplots,
//!   histograms/box plots, error bars, chart-like band/area plots, contour
//!   labels and any interactive output.
//! - **Contours** come from marching squares on a regular sample grid: cells
//!   with a non-finite corner are skipped, plateaus exactly at a level are not
//!   outlined, and the caller chooses the levels and grid. Saddle cells use
//!   the bilinear saddle value (the mean of the four corners) to choose a
//!   connection, so output is deterministic but may differ from a finer grid.
//!   Filled bands are rasterised as one clipped triangle per cell half: the
//!   result is an exact partition, but its SVG size grows as
//!   `O(grid^2 * (levels + 1))` and there is no smoothing, polygon merging or
//!   contour labelling.
//! - **Streamlines** are one-sided (forward) RK4 paths from explicit seeds;
//!   there is no adaptive step, no automatic seeding/density rule, no
//!   backward/two-sided integration and no arrowheads.
//! - **Style** is fixed: one palette, a small colour bar for surfaces/heat
//!   maps, no configuration of axes, fonts, line widths or texts beyond the
//!   legend names.

use std::collections::HashMap;
use std::path::Path;

use ndarray::Array2;

use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;

const WIDTH: f64 = 800.0;
const HEIGHT: f64 = 600.0;
const MARGIN_LEFT: f64 = 70.0;
const MARGIN_RIGHT: f64 = 30.0;
const MARGIN_TOP: f64 = 30.0;
const MARGIN_BOTTOM: f64 = 50.0;
const YAW: f64 = 0.6;
const PITCH: f64 = 0.5;
const PALETTE: [&str; 8] = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf",
];

/// A growing SVG document.
struct Svg {
    body: String,
}

impl Svg {
    const fn new() -> Self {
        Self { body: String::new() }
    }

    fn push(
        &mut self,
        element: &str,
    ) {
        self.body.push_str(element);
        self.body.push('\n');
    }

    fn line(
        &mut self,
        from: (f64, f64),
        to: (f64, f64),
        stroke: &str,
        width: f64,
    ) {
        self.push(&format!(
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{stroke}" stroke-width="{width}"/>"#,
            from.0, from.1, to.0, to.1
        ));
    }

    fn polyline(
        &mut self,
        points: &[(f64, f64)],
        stroke: &str,
    ) {
        let list: Vec<String> = points
            .iter()
            .map(|p| format!("{:.2},{:.2}", p.0, p.1))
            .collect();
        self.push(&format!(
            r#"<polyline fill="none" stroke="{stroke}" stroke-width="1.5" points="{}"/>"#,
            list.join(" ")
        ));
    }

    fn polygon(
        &mut self,
        points: &[(f64, f64)],
        fill: &str,
        stroke: &str,
    ) {
        let list: Vec<String> = points
            .iter()
            .map(|p| format!("{:.2},{:.2}", p.0, p.1))
            .collect();
        self.push(&format!(
            r#"<polygon points="{}" fill="{fill}" stroke="{stroke}" stroke-width="0.4"/>"#,
            list.join(" ")
        ));
    }

    fn rect(
        &mut self,
        origin: (f64, f64),
        size: (f64, f64),
        fill: &str,
    ) {
        self.push(&format!(
            r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{fill}" stroke="{fill}" stroke-width="0.3"/>"#,
            origin.0, origin.1, size.0, size.1
        ));
    }

    fn text(
        &mut self,
        at: (f64, f64),
        anchor: &str,
        content: &str,
    ) {
        self.push(&format!(
            r##"<text x="{:.2}" y="{:.2}" text-anchor="{anchor}" font-family="sans-serif" font-size="12" fill="#222">{}</text>"##,
            at.0,
            at.1,
            escape(content)
        ));
    }

    fn finish(self) -> String {
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{WIDTH}\" height=\"{HEIGHT}\" viewBox=\"0 0 {WIDTH} {HEIGHT}\">\n<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n{}</svg>\n",
            self.body
        )
    }
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn write_svg(
    svg: Svg,
    path: &Path,
) -> Result<(), String> {
    std::fs::write(path, svg.finish()).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

fn format_tick(v: f64) -> String {
    if v == 0.0 {
        return "0".to_owned();
    }
    if v.abs() < 1e-3 || v.abs() >= 1e5 {
        return format!("{v:.2e}");
    }
    let text = format!("{v:.4}");
    text.trim_end_matches('0').trim_end_matches('.').to_owned()
}

/// Round tick positions covering `[lo, hi]`.
fn ticks(
    lo: f64,
    hi: f64,
) -> Vec<f64> {
    let span = hi - lo;
    if !(span.is_finite() && span > 0.0) {
        return vec![lo];
    }
    let raw = span / 6.0;
    let magnitude = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * magnitude)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * magnitude);
    let mut out = Vec::new();
    let mut k = (lo / step).ceil();
    while k * step <= hi + step * 1e-9 && out.len() < 64 {
        out.push(k * step);
        k += 1.0;
    }
    out
}

fn check_range(
    range: (f64, f64),
    what: &str,
) -> Result<(), String> {
    if range.0.is_finite() && range.1.is_finite() && range.0 < range.1 {
        Ok(())
    } else {
        Err(format!("invalid {what} range ({}, {})", range.0, range.1))
    }
}

/// Min and max of the finite values, widened when they coincide.
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
fn extent(values: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
    let (lo, hi) = values
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
            (a.min(v), b.max(v))
        });
    if lo > hi {
        None
    } else if lo == hi {
        Some((lo - 1.0, hi + 1.0))
    } else {
        Some((lo, hi))
    }
}

fn padded(range: (f64, f64)) -> (f64, f64) {
    let pad = (range.1 - range.0) * 0.05;
    (range.0 - pad, range.1 + pad)
}

/// A viridis-like colour map as RGB; `t` is clamped to `[0, 1]`.
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn colormap_rgb(t: f64) -> (u8, u8, u8) {
    const STOPS: [(f64, f64, f64); 5] = [
        (68.0, 1.0, 84.0),
        (59.0, 82.0, 139.0),
        (33.0, 145.0, 140.0),
        (94.0, 201.0, 98.0),
        (253.0, 231.0, 37.0),
    ];
    let t = if t.is_finite() {
        t.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let scaled = t * 4.0;
    let i = (scaled.floor() as usize).min(3);
    let f = scaled - i as f64;
    let (a, b) = (STOPS[i], STOPS[i + 1]);
    let mix = |x: f64, y: f64| (x + (y - x) * f).round() as u8;
    (mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

/// A viridis-like colour map as an SVG `#rrggbb` string; `t` is clamped to
/// `[0, 1]`.
fn colormap(t: f64) -> String {
    let (r, g, b) = colormap_rgb(t);
    format!("#{r:02x}{g:02x}{b:02x}")
}

fn unit(
    v: f64,
    range: (f64, f64),
) -> f64 {
    if range.1 > range.0 {
        (v - range.0) / (range.1 - range.0)
    } else {
        0.5
    }
}

/// Axes, ticks and grid of a 2D plot.
struct Frame {
    x: (f64, f64),
    y: (f64, f64),
}

impl Frame {
    fn px(
        &self,
        x: f64,
    ) -> f64 {
        MARGIN_LEFT + unit(x, self.x) * (WIDTH - MARGIN_LEFT - MARGIN_RIGHT)
    }

    fn py(
        &self,
        y: f64,
    ) -> f64 {
        HEIGHT - MARGIN_BOTTOM - unit(y, self.y) * (HEIGHT - MARGIN_TOP - MARGIN_BOTTOM)
    }

    fn at(
        &self,
        p: (f64, f64),
    ) -> (f64, f64) {
        (self.px(p.0), self.py(p.1))
    }

    fn draw(
        &self,
        svg: &mut Svg,
    ) {
        let (left, right) = (MARGIN_LEFT, WIDTH - MARGIN_RIGHT);
        let (top, bottom) = (MARGIN_TOP, HEIGHT - MARGIN_BOTTOM);
        for t in ticks(self.x.0, self.x.1) {
            let px = self.px(t);
            svg.line((px, top), (px, bottom), "#e6e6e6", 1.0);
            svg.line((px, bottom), (px, bottom + 5.0), "#444", 1.0);
            svg.text((px, bottom + 19.0), "middle", &format_tick(t));
        }
        for t in ticks(self.y.0, self.y.1) {
            let py = self.py(t);
            svg.line((left, py), (right, py), "#e6e6e6", 1.0);
            svg.line((left - 5.0, py), (left, py), "#444", 1.0);
            svg.text((left - 8.0, py + 4.0), "end", &format_tick(t));
        }
        svg.line((left, bottom), (right, bottom), "#222", 1.5);
        svg.line((left, top), (left, bottom), "#222", 1.5);
    }
}

/// Splits a point list at non-finite values and draws each run.
fn draw_series(
    svg: &mut Svg,
    frame: &Frame,
    points: &[(f64, f64)],
    stroke: &str,
) {
    let mut run: Vec<(f64, f64)> = Vec::new();
    for &p in points {
        if p.0.is_finite() && p.1.is_finite() {
            run.push(frame.at(p));
        } else {
            if run.len() > 1 {
                svg.polyline(&run, stroke);
            }
            run.clear();
        }
    }
    if run.len() > 1 {
        svg.polyline(&run, stroke);
    }
}

fn linspace(
    range: (f64, f64),
    n: usize,
) -> Vec<f64> {
    let last = (n.max(2) - 1) as f64;
    (0..n.max(2))
        .map(|i| range.0 + (range.1 - range.0) * i as f64 / last)
        .collect()
}

/// Plots `y = f(x)` over `range` with `samples` points. Non-finite values
/// interrupt the curve.
///
/// # Errors
/// Fails for an invalid range, fewer than two samples, a function that is
/// nowhere finite, or when the file cannot be written.
pub fn plot_function_2d(
    f: impl Fn(f64) -> f64,
    range: (f64, f64),
    samples: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(range, "x")?;
    if samples < 2 {
        return Err("at least two samples are required".to_owned());
    }
    let points: Vec<(f64, f64)> = linspace(range, samples)
        .into_iter()
        .map(|x| (x, f(x)))
        .collect();
    plot_series_2d(&[("f(x)".to_owned(), points)], path)
}

/// Samples from which [`plot_function_2d_par`] / [`plot_term_2d`] evaluate concurrently.
const PLOT_PAR_SAMPLES: usize = 256;

/// [`plot_function_2d`] for a `Sync` function: the samples are evaluated concurrently (collected in
/// order, so the picture is identical) from `PLOT_PAR_SAMPLES` points on.
///
/// # Errors
/// As [`plot_function_2d`].
pub fn plot_function_2d_par(
    f: impl Fn(f64) -> f64 + Sync + Send,
    range: (f64, f64),
    samples: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(range, "x")?;
    if samples < 2 {
        return Err("at least two samples are required".to_owned());
    }
    let xs = linspace(range, samples);
    let ys = crate::kernels::par::map_slice(&xs, PLOT_PAR_SAMPLES, |&x| f(x));
    let points: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
    plot_series_2d(&[("f(x)".to_owned(), points)], path)
}

/// [`plot_surface_3d`] for a `Sync` function: the rows of the grid are evaluated concurrently.
///
/// # Errors
/// As [`plot_surface_3d`].
pub fn plot_surface_3d_par(
    f: impl Fn(f64, f64) -> f64 + Sync + Send,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    let (xs, ys) = (linspace(x_range, grid), linspace(y_range, grid));
    let z: Vec<Vec<f64>> =
        crate::kernels::par::map_slice(&ys, 16, |&y| xs.iter().map(|&x| f(x, y)).collect());
    let mut svg = Svg::new();
    draw_surface(&mut svg, &xs, &ys, &z)?;
    write_svg(svg, path)
}

/// Plots the term `node` as a function of the symbol `var`, evaluating it
/// with [`Graph::eval`]. Points where the term does not evaluate to a
/// number are left out.
///
/// # Errors
/// Same as [`plot_function_2d`]; additionally fails when the term has free
/// symbols other than `var`, since then no point evaluates.
pub fn plot_term_2d(
    graph: &mut Graph,
    node: NodeId,
    var: &str,
    range: (f64, f64),
    samples: usize,
    path: &Path,
) -> Result<(), String> {
    let symbol = graph.interner_mut().symbol(var);
    let value_at = |graph: &Graph, x: f64| {
        let mut env = Env::numeric(0.0);
        env.bind(symbol, x);
        graph.eval(node, &env).unwrap_or(f64::NAN)
    };
    let graph: &Graph = graph;
    plot_function_2d_par(|x| value_at(graph, x), range, samples, path)
}

/// Plots several named point series on common axes, with a legend.
///
/// # Errors
/// Fails when there is no series, no finite point, or the file cannot be
/// written.
pub fn plot_series_2d(
    series: &[(String, Vec<(f64, f64)>)],
    path: &Path,
) -> Result<(), String> {
    if series.is_empty() {
        return Err("no data series provided".to_owned());
    }
    let finite = |p: &&(f64, f64)| p.0.is_finite() && p.1.is_finite();
    let all = || series.iter().flat_map(|(_, pts)| pts.iter().filter(finite));
    let x = extent(all().map(|p| p.0)).ok_or("no finite data points")?;
    let y = extent(all().map(|p| p.1)).ok_or("no finite data points")?;
    let frame = Frame { x, y: padded(y) };
    let mut svg = Svg::new();
    frame.draw(&mut svg);
    for (i, (name, points)) in series.iter().enumerate() {
        let color = PALETTE[i % PALETTE.len()];
        draw_series(&mut svg, &frame, points, color);
        let ly = MARGIN_TOP + 14.0 + 16.0 * i as f64;
        let lx = WIDTH - MARGIN_RIGHT - 120.0;
        svg.line((lx, ly - 4.0), (lx + 22.0, ly - 4.0), color, 2.0);
        svg.text((lx + 28.0, ly), "start", name);
    }
    write_svg(svg, path)
}

fn arrow(
    svg: &mut Svg,
    from: (f64, f64),
    to: (f64, f64),
    color: &str,
) {
    svg.line(from, to, color, 1.2);
    let (dx, dy) = (to.0 - from.0, to.1 - from.1);
    let len = dx.hypot(dy);
    if len < 1e-9 {
        return;
    }
    let (ux, uy) = (dx / len, dy / len);
    let head = (len * 0.3).min(7.0);
    let base = (to.0 - ux * head, to.1 - uy * head);
    svg.polygon(
        &[
            to,
            (base.0 - uy * head * 0.4, base.1 + ux * head * 0.4),
            (base.0 + uy * head * 0.4, base.1 - ux * head * 0.4),
        ],
        color,
        color,
    );
}

/// Plots the 2D vector field `f(x, y) = (u, v)` as arrows on a
/// `grid` x `grid` lattice, coloured by magnitude.
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, or when the file cannot be
/// written.
pub fn plot_vector_field_2d(
    f: impl Fn(f64, f64) -> (f64, f64),
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    let (xs, ys) = (linspace(x_range, grid), linspace(y_range, grid));
    let field: Vec<(f64, f64, f64, f64)> = ys
        .iter()
        .flat_map(|&y| xs.iter().map(move |&x| (x, y)))
        .map(|(x, y)| {
            let (u, v) = f(x, y);
            (x, y, u, v)
        })
        .filter(|&(_, _, u, v)| u.is_finite() && v.is_finite())
        .collect();
    let max = field
        .iter()
        .map(|&(_, _, u, v)| u.hypot(v))
        .fold(0.0, f64::max);
    let frame = Frame {
        x: x_range,
        y: y_range,
    };
    let mut svg = Svg::new();
    frame.draw(&mut svg);
    let cell_x = (WIDTH - MARGIN_LEFT - MARGIN_RIGHT) / (grid - 1) as f64;
    let cell_y = (HEIGHT - MARGIN_TOP - MARGIN_BOTTOM) / (grid - 1) as f64;
    let reach = 0.85 * cell_x.min(cell_y);
    for (x, y, u, v) in field {
        let m = u.hypot(v);
        if m == 0.0 || max == 0.0 {
            continue;
        }
        let from = frame.at((x, y));
        let to = (from.0 + u / max * reach, from.1 - v / max * reach);
        arrow(&mut svg, from, to, &colormap(m / max));
    }
    write_svg(svg, path)
}

/// Plots the 2-D vector field `f(x, y) = (u, v)` as a quiver of scaled arrows.
///
/// There is one arrow per grid sample and its length is proportional to the
/// field magnitude relative to the largest sampled magnitude, so zero vectors
/// are omitted and long arrows mark the strong regions.
///
/// `scale` is the fraction of the smaller grid cell used by the longest arrow;
/// it must be finite and positive.
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, a non-finite or non-positive
/// `scale`, or when the file cannot be written.
#[allow(clippy::float_cmp)] // exact zero is the "field vanishes" sentinel
pub fn plot_quiver_2d(
    f: impl Fn(f64, f64) -> (f64, f64),
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    scale: f64,
    path: &Path,
) -> Result<(), String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    if !(scale.is_finite() && scale > 0.0) {
        return Err(format!("invalid quiver scale ({scale})"));
    }
    let (xs, ys) = (linspace(x_range, grid), linspace(y_range, grid));
    let field: Vec<(f64, f64, f64, f64)> = ys
        .iter()
        .flat_map(|&y| xs.iter().map(move |&x| (x, y)))
        .map(|(x, y)| {
            let (u, v) = f(x, y);
            (x, y, u, v)
        })
        .filter(|&(_, _, u, v)| u.is_finite() && v.is_finite())
        .collect();
    let max = field
        .iter()
        .map(|&(_, _, u, v)| u.hypot(v))
        .fold(0.0, f64::max);
    let frame = Frame {
        x: x_range,
        y: y_range,
    };
    let mut svg = Svg::new();
    frame.draw(&mut svg);
    if max > 0.0 {
        let cell_x = (WIDTH - MARGIN_LEFT - MARGIN_RIGHT) / (grid - 1) as f64;
        let cell_y = (HEIGHT - MARGIN_TOP - MARGIN_BOTTOM) / (grid - 1) as f64;
        let reach = scale * cell_x.min(cell_y);
        for (x, y, u, v) in field {
            let m = u.hypot(v);
            if m == 0.0 {
                continue;
            }
            let part = reach * (m / max);
            let from = frame.at((x, y));
            let to = (from.0 + u / m * part, from.1 - v / m * part);
            arrow(&mut svg, from, to, PALETTE[0]);
        }
    }
    write_svg(svg, path)
}

/// One contour line of a scalar field at a fixed level, in data coordinates.
///
/// A closed contour repeats its first point as its last point.
#[derive(Clone, Debug, PartialEq)]
pub struct ContourLine {
    /// The level this line belongs to.
    pub level: f64,
    /// The polyline, in plot coordinates.
    pub points: Vec<(f64, f64)>,
}

/// Interpolates the crossing of `level` on the edge `p -> q`.
///
/// Returns the endpoint itself when the interpolation parameter is exactly
/// `0` or `1`, so a sample that lies exactly on the level is reproduced
/// bit-exactly by every cell that touches it.
fn crossing(
    p: (f64, f64),
    vp: f64,
    q: (f64, f64),
    vq: f64,
    level: f64,
) -> (f64, f64) {
    let denom = vq - vp;
    let t = if denom == 0.0 {
        0.5
    } else {
        ((level - vp) / denom).clamp(0.0, 1.0)
    };
    if t <= 0.0 {
        return p;
    }
    if t >= 1.0 {
        return q;
    }
    (p.0 + (q.0 - p.0) * t, p.1 + (q.1 - p.1) * t)
}

/// One marching-squares segment: the two cell-edge crossings it connects.
type ContourSegment = ((f64, f64), (f64, f64));

/// Marching-squares segments of one cell.
///
/// `edges` are the crossings on the bottom, right, top and left edge; `code`
/// has one bit per corner above the level (bottom-left is bit 0, counter-
/// clockwise). In the two saddle cases (5 and 10) the bilinear saddle value,
/// the mean of the four corners, decides how the cell is connected.
fn cell_segments(
    edges: [(f64, f64); 4],
    code: u8,
    center_above: bool,
) -> Vec<ContourSegment> {
    let pair = |a: usize, b: usize| vec![(edges[a], edges[b])];
    match code {
        | 1 | 14 => pair(3, 0),
        | 2 | 13 => pair(0, 1),
        | 3 | 12 => pair(3, 1),
        | 4 | 11 => pair(1, 2),
        | 6 | 9 => pair(0, 2),
        | 7 | 8 => pair(3, 2),
        | 5 | 10 => {
            // The bilinear saddle value decides the connection: for case 5 the
            // connected branch is chosen when the centre is above the level,
            // for case 10 when it is below.
            if (code == 5) == center_above {
                let mut out = pair(0, 1);
                out.extend(pair(2, 3));
                out
            } else {
                let mut out = pair(3, 0);
                out.extend(pair(1, 2));
                out
            }
        },
        | _ => Vec::new(),
    }
}

/// Canonical key of a point: exact bit pattern with `-0.0` normalised, so
/// points computed in different cells still compare equal.
fn point_key(p: (f64, f64)) -> String {
    format!("{:?},{:?}", p.0 + 0.0, p.1 + 0.0)
}

/// Joins segments that share an endpoint into polylines.
///
/// Endpoints are compared by [`point_key`], and shared grid edges are computed
/// by one canonical formula, so this is deterministic. Closed contours keep
/// their first point repeated at the end; segments are otherwise chained in
/// scanning order.
fn chain_segments(segments: &[ContourSegment]) -> Vec<Vec<(f64, f64)>> {
    let key = point_key;
    let mut adjacency: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, &(a, b)) in segments.iter().enumerate() {
        adjacency.entry(key(a)).or_default().push(i);
        adjacency.entry(key(b)).or_default().push(i);
    }
    let mut used = vec![false; segments.len()];
    let mut lines: Vec<Vec<(f64, f64)>> = Vec::new();
    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let (a, b) = segments[start];
        let mut points = vec![a, b];
        loop {
            let tail = key(*points.last().unwrap());
            let next = adjacency
                .get(&tail)
                .and_then(|ids| ids.iter().copied().find(|&i| !used[i]));
            let Some(i) = next else {
                break;
            };
            used[i] = true;
            let (p, q) = segments[i];
            if key(p) == tail {
                points.push(q);
            } else {
                points.push(p);
            }
        }
        loop {
            let head = key(points[0]);
            let next = adjacency
                .get(&head)
                .and_then(|ids| ids.iter().copied().find(|&i| !used[i]));
            let Some(i) = next else {
                break;
            };
            used[i] = true;
            let (p, q) = segments[i];
            if key(q) == head {
                points.insert(0, p);
            } else {
                points.insert(0, q);
            }
        }
        lines.push(points);
    }
    lines
}

/// Extracts iso-lines of `f` over a rectangle with marching squares on a
/// regular `grid` x `grid` sample lattice.
///
/// Returns one polyline per connected contour component. Components are
/// ordered by level (in the given order) and then by the scan order of their
/// first segment, so the result is deterministic. Non-finite levels are
/// ignored; cells with a non-finite corner are skipped, and cells whose
/// corners are all exactly on a level are not outlined. A sample exactly on a
/// level is treated as below it, so the contour passes through that sample
/// rather than through a degenerate cell edge (zero-length segments produced by
/// cells the contour merely touches are dropped). Saddle cells are resolved
/// with the bilinear saddle value (the mean of the four corners).
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, or an empty `levels` slice.
pub fn contour_lines_2d(
    f: impl Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    levels: &[f64],
) -> Result<Vec<ContourLine>, String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    if levels.is_empty() {
        return Err("at least one contour level is required".to_owned());
    }
    let xs = linspace(x_range, grid);
    let ys = linspace(y_range, grid);
    let z: Vec<Vec<f64>> = ys
        .iter()
        .map(|&y| xs.iter().map(|&x| f(x, y)).collect())
        .collect();
    let mut lines: Vec<ContourLine> = Vec::new();
    for &level in levels {
        if !level.is_finite() {
            continue;
        }
        let horizontal = |row: usize, col: usize| {
            crossing(
                (xs[col], ys[row]),
                z[row][col],
                (xs[col + 1], ys[row]),
                z[row][col + 1],
                level,
            )
        };
        let vertical = |row: usize, col: usize| {
            crossing(
                (xs[col], ys[row]),
                z[row][col],
                (xs[col], ys[row + 1]),
                z[row + 1][col],
                level,
            )
        };
        let mut segments = Vec::new();
        for i in 0..grid - 1 {
            for j in 0..grid - 1 {
                let corners = [z[i][j], z[i][j + 1], z[i + 1][j + 1], z[i + 1][j]];
                if !corners.iter().all(|v| v.is_finite()) {
                    continue;
                }
                let mut code = 0u8;
                for (bit, &v) in corners.iter().enumerate() {
                    // Strictly above: a sample exactly on the level is treated
                    // as below, so crossings land on the sample itself.
                    if v > level {
                        code |= 1 << bit;
                    }
                }
                // Shared edges are computed by the same call on both sides, so
                // neighbouring cells produce bit-identical crossing points.
                let edges = [
                    horizontal(i, j),
                    vertical(i, j + 1),
                    horizontal(i + 1, j),
                    vertical(i, j),
                ];
                let center_above = corners.iter().sum::<f64>() / 4.0 > level;
                segments.extend(cell_segments(edges, code, center_above));
            }
        }
        // A contour that only touches a sample degenerates to a zero-length
        // segment in the neighbouring cell; dropping it keeps tangent points
        // connected through the two real segments that meet at the sample.
        segments.retain(|&(a, b)| point_key(a) != point_key(b));
        for points in chain_segments(&segments) {
            lines.push(ContourLine { level, points });
        }
    }
    Ok(lines)
}

/// Plots the contour lines of `f` over the rectangle as SVG: one polyline per
/// connected contour component, coloured per level.
///
/// # Errors
/// As [`contour_lines_2d`], plus file-writing errors.
pub fn plot_contour_2d(
    f: impl Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    levels: &[f64],
    path: &Path,
) -> Result<(), String> {
    let lines = contour_lines_2d(f, x_range, y_range, grid, levels)?;
    let frame = Frame {
        x: x_range,
        y: y_range,
    };
    let mut svg = Svg::new();
    frame.draw(&mut svg);
    let mut seen: Vec<u64> = Vec::new();
    for line in &lines {
        let bits = line.level.to_bits();
        let color = match seen.iter().position(|&b| b == bits) {
            | Some(i) => PALETTE[i % PALETTE.len()],
            | None => {
                seen.push(bits);
                PALETTE[(seen.len() - 1) % PALETTE.len()]
            },
        };
        draw_series(&mut svg, &frame, &line.points, color);
    }
    write_svg(svg, path)
}

/// One filled contour band: a scalar range and the filled polygons covering
/// it, in data coordinates.
///
/// Polygons are closed rings (the last vertex is implicit) with at least three
/// vertices. Adjacent polygons share edges exactly, so drawing them with the
/// same fill leaves no seams or double-fill artefacts.
#[derive(Clone, Debug, PartialEq)]
pub struct ContourBand {
    /// Lower bound of the band (`f64::NEG_INFINITY` below the first level).
    pub low: f64,
    /// Upper bound of the band (`f64::INFINITY` above the last level).
    pub high: f64,
    /// The filled polygons of the band.
    pub polygons: Vec<Vec<(f64, f64)>>,
}

/// A polygon vertex carrying the linear interpolant of the scalar field, so
/// clipping can intersect on the field value instead of the coordinates.
type BandVertex = (f64, f64, f64);

/// Sutherland-Hodgman clip of a convex polygon to one side of `value == level`.
fn clip_values(
    poly: &[BandVertex],
    level: f64,
    keep_above: bool,
) -> Vec<BandVertex> {
    let inside = |v: f64| {
        if keep_above {
            v >= level
        } else {
            v <= level
        }
    };
    let mut out = Vec::with_capacity(poly.len() + 2);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let (a_in, b_in) = (inside(a.2), inside(b.2));
        if a_in {
            out.push(a);
        }
        if a_in != b_in {
            let denom = b.2 - a.2;
            let t = if denom == 0.0 {
                0.5
            } else {
                ((level - a.2) / denom).clamp(0.0, 1.0)
            };
            out.push((a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t, level));
        }
    }
    out
}

/// Clips a triangle to the strip `low <= value <= high`; an infinite bound
/// leaves that side untouched.
fn clip_band(
    poly: &[BandVertex],
    low: f64,
    high: f64,
) -> Vec<BandVertex> {
    let mut out = poly.to_vec();
    if low.is_finite() {
        out = clip_values(&out, low, true);
    }
    if high.is_finite() {
        out = clip_values(&out, high, false);
    }
    out
}

/// Twice the signed area of a polygon (shoelace formula).
fn band_area2(points: &[BandVertex]) -> f64 {
    let mut acc = 0.0;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        acc += a.0 * b.1 - b.0 * a.1;
    }
    acc
}

/// Validates contour boundary levels: non-empty, finite and strictly
/// increasing.
fn check_levels(levels: &[f64]) -> Result<(), String> {
    if levels.is_empty() {
        return Err("at least one contour level is required".to_owned());
    }
    if levels.iter().any(|v| !v.is_finite()) {
        return Err("contour levels must be finite".to_owned());
    }
    if levels.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err("contour levels must be strictly increasing".to_owned());
    }
    Ok(())
}

/// Splits the sampled rectangle into filled level bands.
///
/// There is one band per gap between consecutive `levels`, plus the two open
/// outer bands `(-inf, levels[0])` and `(levels[last], +inf)`. Each grid cell
/// is split on its diagonal into two triangles; the scalar field is linearly
/// interpolated on each triangle and the triangle is clipped to the band with
/// Sutherland-Hodgman. That makes every polygon convex and the bands an exact
/// partition of the rectangle (no overlap and no double fill), with band
/// boundaries piecewise linear on the grid diagonals. Polygons are returned in
/// band order and then row-major cell order, so the result is deterministic.
///
/// Cells with a non-finite corner are skipped. Degenerate zero-area clipped
/// polygons (a level touching a sample) are dropped.
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, an empty `levels` slice, or levels
/// that are not finite and strictly increasing.
#[allow(clippy::float_cmp)] // exact zero area identifies a degenerate clipped polygon
pub fn contour_bands_2d(
    f: impl Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    levels: &[f64],
) -> Result<Vec<ContourBand>, String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    check_levels(levels)?;
    let xs = linspace(x_range, grid);
    let ys = linspace(y_range, grid);
    let z: Vec<Vec<f64>> = ys
        .iter()
        .map(|&y| xs.iter().map(|&x| f(x, y)).collect())
        .collect();
    let edges: Vec<f64> = std::iter::once(f64::NEG_INFINITY)
        .chain(levels.iter().copied())
        .chain(std::iter::once(f64::INFINITY))
        .collect();
    let mut bands: Vec<ContourBand> = edges
        .windows(2)
        .map(|pair| {
            ContourBand {
                low: pair[0],
                high: pair[1],
                polygons: Vec::new(),
            }
        })
        .collect();
    for i in 0..grid - 1 {
        for j in 0..grid - 1 {
            let corners: [BandVertex; 4] = [
                (xs[j], ys[i], z[i][j]),
                (xs[j + 1], ys[i], z[i][j + 1]),
                (xs[j + 1], ys[i + 1], z[i + 1][j + 1]),
                (xs[j], ys[i + 1], z[i + 1][j]),
            ];
            if !corners.iter().all(|corner| corner.2.is_finite()) {
                continue;
            }
            for triangle in [[0usize, 1, 2], [0, 2, 3]] {
                let tri = [
                    corners[triangle[0]],
                    corners[triangle[1]],
                    corners[triangle[2]],
                ];
                for band in &mut bands {
                    let clipped = clip_band(&tri, band.low, band.high);
                    if clipped.len() >= 3 && band_area2(&clipped) != 0.0 {
                        band.polygons
                            .push(clipped.iter().map(|vertex| (vertex.0, vertex.1)).collect());
                    }
                }
            }
        }
    }
    Ok(bands)
}

/// Plots filled contour bands of `f` over the rectangle as SVG: one filled
/// polygon per clipped triangle, coloured by band index (low to high).
///
/// The axes and grid are drawn on top of the fills.
///
/// # Errors
/// As [`contour_bands_2d`], plus file-writing errors.
pub fn plot_contour_bands_2d(
    f: impl Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    levels: &[f64],
    path: &Path,
) -> Result<(), String> {
    let bands = contour_bands_2d(f, x_range, y_range, grid, levels)?;
    let frame = Frame {
        x: x_range,
        y: y_range,
    };
    let mut svg = Svg::new();
    let last = bands.len().saturating_sub(1).max(1) as f64;
    for (index, band) in bands.iter().enumerate() {
        let color = colormap(index as f64 / last);
        for polygon in &band.polygons {
            let screen: Vec<(f64, f64)> = polygon.iter().map(|&p| frame.at(p)).collect();
            svg.polygon(&screen, &color, &color);
        }
    }
    frame.draw(&mut svg);
    write_svg(svg, path)
}

/// One classical fourth-order Runge-Kutta step of the planar field.
fn rk4_step(
    field: &impl Fn(f64, f64) -> (f64, f64),
    x: f64,
    y: f64,
    h: f64,
) -> Option<(f64, f64)> {
    let (k1x, k1y) = field(x, y);
    let (k2x, k2y) = field(x + 0.5 * h * k1x, y + 0.5 * h * k1y);
    let (k3x, k3y) = field(x + 0.5 * h * k2x, y + 0.5 * h * k2y);
    let (k4x, k4y) = field(x + h * k3x, y + h * k3y);
    let next = (
        x + h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
        y + h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y),
    );
    (next.0.is_finite() && next.1.is_finite()).then_some(next)
}

/// Integrates the planar vector field `f(x, y) = (u, v)` forward from each
/// seed with fourth-order Runge-Kutta and returns one polyline per seed.
///
/// A path starts at its seed and stops as soon as a step would leave the
/// rectangle, the field becomes non-finite or exactly zero, or `max_steps`
/// steps were taken. Integration is forward only (no backward or two-sided
/// lines).
///
/// # Errors
/// Fails for an invalid range, a non-positive or non-finite `step`,
/// `max_steps == 0`, an empty seed list, or a non-finite seed outside the
/// rectangle.
#[allow(clippy::float_cmp)] // exact zero is the "field vanishes" sentinel
pub fn streamline_paths_2d(
    field: impl Fn(f64, f64) -> (f64, f64),
    x_range: (f64, f64),
    y_range: (f64, f64),
    seeds: &[(f64, f64)],
    step: f64,
    max_steps: usize,
) -> Result<Vec<Vec<(f64, f64)>>, String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if !(step.is_finite() && step > 0.0) {
        return Err(format!("invalid streamline step ({step})"));
    }
    if max_steps == 0 {
        return Err("streamlines need at least one step".to_owned());
    }
    if seeds.is_empty() {
        return Err("no streamline seeds provided".to_owned());
    }
    let inside = |p: (f64, f64)| {
        // Bind the bounds first: comparing `p.0`/`p.1` against tuple fields
        // directly trips `clippy::suspicious_operation_groupings` in the
        // reverse reading, even though the chain below is the intended one.
        let (x0, x1) = x_range;
        let (y0, y1) = y_range;
        p.0 >= x0 && p.0 <= x1 && p.1 >= y0 && p.1 <= y1
    };
    let mut paths = Vec::with_capacity(seeds.len());
    for &seed in seeds {
        if !(seed.0.is_finite() && seed.1.is_finite() && inside(seed)) {
            return Err(format!(
                "streamline seed ({}, {}) is outside the ranges",
                seed.0, seed.1
            ));
        }
        let mut points = vec![seed];
        let (mut x, mut y) = seed;
        for _ in 0..max_steps {
            let (u, v) = field(x, y);
            if !u.is_finite() || !v.is_finite() || (u == 0.0 && v == 0.0) {
                break;
            }
            let Some(next) = rk4_step(&field, x, y, step) else {
                break;
            };
            if !inside(next) {
                break;
            }
            points.push(next);
            (x, y) = next;
        }
        paths.push(points);
    }
    Ok(paths)
}

/// Plots streamlines of the planar vector field as SVG, one polyline per seed.
///
/// # Errors
/// As [`streamline_paths_2d`], plus file-writing errors.
pub fn plot_streamlines_2d(
    field: impl Fn(f64, f64) -> (f64, f64),
    x_range: (f64, f64),
    y_range: (f64, f64),
    seeds: &[(f64, f64)],
    step: f64,
    max_steps: usize,
    path: &Path,
) -> Result<(), String> {
    let paths = streamline_paths_2d(field, x_range, y_range, seeds, step, max_steps)?;
    let frame = Frame {
        x: x_range,
        y: y_range,
    };
    let mut svg = Svg::new();
    frame.draw(&mut svg);
    for (i, points) in paths.iter().enumerate() {
        draw_series(&mut svg, &frame, points, PALETTE[i % PALETTE.len()]);
    }
    write_svg(svg, path)
}

/// Fixed oblique projection of the unit cube `[-1, 1]^3` to the plane;
/// returns `(screen x, screen y, depth)`.
fn project(p: [f64; 3]) -> (f64, f64, f64) {
    let (sy, cy) = YAW.sin_cos();
    let (sp, cp) = PITCH.sin_cos();
    let x1 = p[0] * cy - p[1] * sy;
    let y1 = p[0] * sy + p[1] * cy;
    (x1, p[2] * cp + y1 * sp, y1)
}

/// Maps data coordinates into the unit cube and onto the page.
struct Scene {
    ranges: [(f64, f64); 3],
    scale: f64,
    offset: (f64, f64),
}

impl Scene {
    fn new(ranges: [(f64, f64); 3]) -> Self {
        let corners: Vec<(f64, f64, f64)> = (0..8)
            .map(|i| {
                project([
                    if i & 1 == 0 { -1.0 } else { 1.0 },
                    if i & 2 == 0 { -1.0 } else { 1.0 },
                    if i & 4 == 0 { -1.0 } else { 1.0 },
                ])
            })
            .collect();
        let (x0, x1) = extent(corners.iter().map(|c| c.0)).unwrap_or((-1.0, 1.0));
        let (y0, y1) = extent(corners.iter().map(|c| c.1)).unwrap_or((-1.0, 1.0));
        let avail = (WIDTH - 2.0 * 90.0, HEIGHT - 2.0 * 70.0);
        let scale = (avail.0 / (x1 - x0)).min(avail.1 / (y1 - y0));
        Self {
            ranges,
            scale,
            offset: (
                WIDTH / 2.0 - scale * (x0 + x1) / 2.0,
                HEIGHT / 2.0 + scale * (y0 + y1) / 2.0,
            ),
        }
    }

    fn cube(
        &self,
        p: [f64; 3],
    ) -> [f64; 3] {
        [
            2.0 * unit(p[0], self.ranges[0]) - 1.0,
            2.0 * unit(p[1], self.ranges[1]) - 1.0,
            2.0 * unit(p[2], self.ranges[2]) - 1.0,
        ]
    }

    /// Page position and depth of a data point.
    fn at(
        &self,
        p: [f64; 3],
    ) -> (f64, f64, f64) {
        let (x, y, d) = project(self.cube(p));
        (
            self.offset.0 + self.scale * x,
            self.offset.1 - self.scale * y,
            d,
        )
    }

    fn draw_box(
        &self,
        svg: &mut Svg,
    ) {
        let corner = |i: usize| {
            let pick = |axis: usize| {
                if i >> axis & 1 == 0 {
                    self.ranges[axis].0
                } else {
                    self.ranges[axis].1
                }
            };
            let (x, y, _) = self.at([pick(0), pick(1), pick(2)]);
            (x, y)
        };
        for i in 0..8usize {
            for axis in 0..3 {
                let j = i | (1 << axis);
                if j != i {
                    svg.line(
                        corner(i),
                        corner(j),
                        if i == 0 { "#222" } else { "#bbb" },
                        if i == 0 { 1.5 } else { 0.8 },
                    );
                }
            }
        }
        for (axis, name) in ["x", "y", "z"].iter().enumerate() {
            let (lo, hi) = (corner(0), corner(1 << axis));
            svg.text(
                (lo.0, lo.1 + 16.0),
                "middle",
                &format_tick(self.ranges[axis].0),
            );
            svg.text(
                (hi.0, hi.1 + 16.0),
                "middle",
                &format!("{name} = {}", format_tick(self.ranges[axis].1)),
            );
        }
    }
}

/// A projected quad: depth, outline, height fraction.
type Quad = (f64, Vec<(f64, f64)>, f64);

/// Paints the grid `z[row][col]` over the given axes as coloured quads.
fn draw_surface(
    svg: &mut Svg,
    xs: &[f64],
    ys: &[f64],
    z: &[Vec<f64>],
) -> Result<(), String> {
    let zr = extent(z.iter().flatten().copied()).ok_or("no finite surface values")?;
    let xr = (
        xs.first().copied().unwrap_or(0.0),
        xs.last().copied().unwrap_or(1.0),
    );
    let yr = (
        ys.first().copied().unwrap_or(0.0),
        ys.last().copied().unwrap_or(1.0),
    );
    let scene = Scene::new([xr, yr, zr]);
    scene.draw_box(svg);
    let mut quads: Vec<Quad> = Vec::new();
    for (r, pair) in z.windows(2).enumerate() {
        for c in 0..xs.len().saturating_sub(1) {
            let corners = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)];
            let pts: Vec<(f64, f64, f64, f64)> = corners
                .iter()
                .filter_map(|&(i, j)| {
                    let v = if i == r {
                        pair[0][j]
                    } else {
                        pair[1][j]
                    };
                    v.is_finite().then(|| {
                        let (px, py, d) = scene.at([xs[j], ys[i], v]);
                        (px, py, d, v)
                    })
                })
                .collect();
            if pts.len() == 4 {
                let depth = pts.iter().map(|p| p.2).sum::<f64>() / 4.0;
                let mean = pts.iter().map(|p| p.3).sum::<f64>() / 4.0;
                quads.push((
                    depth,
                    pts.iter().map(|p| (p.0, p.1)).collect(),
                    unit(mean, zr),
                ));
            }
        }
    }
    quads.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (_, poly, t) in quads {
        svg.polygon(&poly, &colormap(t), "#00000033");
    }
    colorbar(svg, zr);
    Ok(())
}

fn colorbar(
    svg: &mut Svg,
    range: (f64, f64),
) {
    let (x, top, h) = (WIDTH - 24.0, 60.0, 200.0);
    let steps = 40;
    for i in 0..steps {
        let t = 1.0 - (f64::from(i) + 0.5) / f64::from(steps);
        svg.rect(
            (x, top + h * f64::from(i) / f64::from(steps)),
            (10.0, h / f64::from(steps) + 0.5),
            &colormap(t),
        );
    }
    svg.text((x + 5.0, top - 6.0), "middle", &format_tick(range.1));
    svg.text((x + 5.0, top + h + 14.0), "middle", &format_tick(range.0));
}

/// Plots the surface `z = f(x, y)` on a `grid` x `grid` lattice as height
/// coloured cells in a fixed oblique projection. Non-finite cells are left
/// out.
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, no finite value, or when the
/// file cannot be written.
pub fn plot_surface_3d(
    f: impl Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(x_range, "x")?;
    check_range(y_range, "y")?;
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    let (xs, ys) = (linspace(x_range, grid), linspace(y_range, grid));
    let z: Vec<Vec<f64>> = ys
        .iter()
        .map(|&y| xs.iter().map(|&x| f(x, y)).collect())
        .collect();
    let mut svg = Svg::new();
    draw_surface(&mut svg, &xs, &ys, &z)?;
    write_svg(svg, path)
}

/// Plots a 2D array as a surface; the column index is `x`, the row index
/// is `y` and the entry is the height.
///
/// # Errors
/// Fails when the array has fewer than 2 rows or columns, holds no finite
/// value, or when the file cannot be written.
pub fn plot_surface_2d(
    data: &Array2<f64>,
    path: &Path,
) -> Result<(), String> {
    let (rows, cols) = data.dim();
    if rows < 2 || cols < 2 {
        return Err("a surface needs at least a 2 x 2 array".to_owned());
    }
    let xs: Vec<f64> = (0..cols).map(|i| i as f64).collect();
    let ys: Vec<f64> = (0..rows).map(|i| i as f64).collect();
    let z: Vec<Vec<f64>> = data.rows().into_iter().map(|r| r.to_vec()).collect();
    let mut svg = Svg::new();
    draw_surface(&mut svg, &xs, &ys, &z)?;
    write_svg(svg, path)
}

fn draw_path_3d(
    svg: &mut Svg,
    scene: &Scene,
    points: &[[f64; 3]],
    color: &str,
) {
    let mut run: Vec<(f64, f64)> = Vec::new();
    for p in points {
        if p.iter().all(|v| v.is_finite()) {
            let (x, y, _) = scene.at(*p);
            run.push((x, y));
        } else {
            if run.len() > 1 {
                svg.polyline(&run, color);
            }
            run.clear();
        }
    }
    if run.len() > 1 {
        svg.polyline(&run, color);
    }
}

fn bounds_3d(points: &[[f64; 3]]) -> Result<[(f64, f64); 3], String> {
    let axis = |k: usize| {
        extent(points.iter().map(move |p| p[k])).ok_or_else(|| "no finite points".to_owned())
    };
    Ok([axis(0)?, axis(1)?, axis(2)?])
}

/// Plots the space curve `t -> (x, y, z)` over `range` with `samples`
/// points.
///
/// # Errors
/// Fails for an invalid range, fewer than two samples, no finite point, or
/// when the file cannot be written.
#[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
pub fn plot_parametric_curve_3d(
    f: impl Fn(f64) -> (f64, f64, f64),
    range: (f64, f64),
    samples: usize,
    path: &Path,
) -> Result<(), String> {
    check_range(range, "t")?;
    if samples < 2 {
        return Err("at least two samples are required".to_owned());
    }
    let points: Vec<[f64; 3]> = linspace(range, samples)
        .into_iter()
        .map(|t| {
            let (x, y, z) = f(t);
            [x, y, z]
        })
        .collect();
    plot_3d_path_from_points(&points, path)
}

/// Plots a polyline through the given 3D points.
///
/// # Errors
/// Fails when fewer than two points are given, none is finite, or the file
/// cannot be written.
pub fn plot_3d_path_from_points(
    points: &[[f64; 3]],
    path: &Path,
) -> Result<(), String> {
    if points.len() < 2 {
        return Err("a path needs at least two points".to_owned());
    }
    let scene = Scene::new(bounds_3d(points)?);
    let mut svg = Svg::new();
    scene.draw_box(&mut svg);
    draw_path_3d(&mut svg, &scene, points, PALETTE[0]);
    write_svg(svg, path)
}

/// A projected arrow: depth, tail, tip, relative magnitude.
type Arrow = (f64, (f64, f64), (f64, f64), f64);

/// Plots the 3D vector field `f(x, y, z) = (u, v, w)` as arrows on a
/// `grid`^3 lattice over the three ranges, coloured by magnitude.
///
/// # Errors
/// Fails for an invalid range, `grid < 2`, or when the file cannot be
/// written.
pub fn plot_vector_field_3d(
    f: impl Fn(f64, f64, f64) -> (f64, f64, f64),
    ranges: [(f64, f64); 3],
    grid: usize,
    path: &Path,
) -> Result<(), String> {
    for (range, name) in ranges.iter().zip(["x", "y", "z"]) {
        check_range(*range, name)?;
    }
    if grid < 2 {
        return Err("the grid needs at least two points per axis".to_owned());
    }
    let axes: Vec<Vec<f64>> = ranges.iter().map(|&r| linspace(r, grid)).collect();
    let mut field: Vec<([f64; 3], [f64; 3])> = Vec::new();
    for &z in &axes[2] {
        for &y in &axes[1] {
            for &x in &axes[0] {
                let (u, v, w) = f(x, y, z);
                if u.is_finite() && v.is_finite() && w.is_finite() {
                    field.push(([x, y, z], [u, v, w]));
                }
            }
        }
    }
    let max = field
        .iter()
        .map(|(_, d)| d[0].hypot(d[1]).hypot(d[2]))
        .fold(0.0, f64::max);
    let scene = Scene::new(ranges);
    let mut svg = Svg::new();
    scene.draw_box(&mut svg);
    let reach = 0.8
        * ranges
            .iter()
            .map(|r| r.1 - r.0)
            .fold(f64::INFINITY, f64::min)
        / (grid - 1) as f64;
    let mut arrows: Vec<Arrow> = field
        .iter()
        .filter_map(|(p, d)| {
            let m = d[0].hypot(d[1]).hypot(d[2]);
            if m == 0.0 || max == 0.0 {
                return None;
            }
            let tip = [
                p[0] + d[0] / max * reach,
                p[1] + d[1] / max * reach,
                p[2] + d[2] / max * reach,
            ];
            let (x0, y0, depth) = scene.at(*p);
            let (x1, y1, _) = scene.at(tip);
            Some((depth, (x0, y0), (x1, y1), m / max))
        })
        .collect();
    arrows.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (_, from, to, t) in arrows {
        arrow(&mut svg, from, to, &colormap(t));
    }
    colorbar(&mut svg, (0.0, max));
    write_svg(svg, path)
}

/// Plots a 2D array as a heat map: row 0 on top, one coloured cell per
/// entry, with a colour bar.
///
/// # Errors
/// Fails for an empty array, no finite value, or when the file cannot be
/// written.
pub fn plot_heatmap_2d(
    data: &Array2<f64>,
    path: &Path,
) -> Result<(), String> {
    let (rows, cols) = data.dim();
    if rows == 0 || cols == 0 {
        return Err("the array is empty".to_owned());
    }
    let range = extent(data.iter().copied()).ok_or("no finite values")?;
    let (w, h) = (
        WIDTH - MARGIN_LEFT - MARGIN_RIGHT - 40.0,
        HEIGHT - MARGIN_TOP - MARGIN_BOTTOM,
    );
    let (cw, ch) = (w / cols as f64, h / rows as f64);
    let mut svg = Svg::new();
    for ((r, c), &v) in data.indexed_iter() {
        let fill = if v.is_finite() {
            colormap(unit(v, range))
        } else {
            "#dddddd".to_owned()
        };
        svg.rect(
            (MARGIN_LEFT + c as f64 * cw, MARGIN_TOP + r as f64 * ch),
            (cw, ch),
            &fill,
        );
    }
    svg.line(
        (MARGIN_LEFT, MARGIN_TOP),
        (MARGIN_LEFT + w, MARGIN_TOP),
        "#222",
        1.0,
    );
    svg.line(
        (MARGIN_LEFT, MARGIN_TOP + h),
        (MARGIN_LEFT + w, MARGIN_TOP + h),
        "#222",
        1.0,
    );
    svg.line(
        (MARGIN_LEFT, MARGIN_TOP),
        (MARGIN_LEFT, MARGIN_TOP + h),
        "#222",
        1.0,
    );
    svg.line(
        (MARGIN_LEFT + w, MARGIN_TOP),
        (MARGIN_LEFT + w, MARGIN_TOP + h),
        "#222",
        1.0,
    );
    svg.text((MARGIN_LEFT, MARGIN_TOP + h + 16.0), "middle", "0");
    svg.text(
        (MARGIN_LEFT + w, MARGIN_TOP + h + 16.0),
        "middle",
        &cols.to_string(),
    );
    svg.text((MARGIN_LEFT - 8.0, MARGIN_TOP + 4.0), "end", "0");
    svg.text(
        (MARGIN_LEFT - 8.0, MARGIN_TOP + h + 4.0),
        "end",
        &rows.to_string(),
    );
    colorbar(&mut svg, range);
    write_svg(svg, path)
}

/// Writes a 2-D array as an 8-bit grayscale PGM (P5) raster, row 0 on top and
/// column 0 on the left.
///
/// Finite entries are linearly normalised to `0..=255` (rounded); a non-finite
/// entry becomes `0`, and a constant finite array becomes mid-grey (`128`).
/// This is the raster counterpart of [`plot_heatmap_2d`]: no axes, colour bar
/// or annotations are added.
///
/// # Errors
/// Fails for an empty array, an array without a finite entry, or when the file
/// cannot be written.
pub fn plot_scalar_field_pgm(
    data: &Array2<f64>,
    path: &Path,
) -> Result<(), String> {
    let (rows, cols) = data.dim();
    if rows == 0 || cols == 0 {
        return Err("the array is empty".to_owned());
    }
    let range = extent(data.iter().copied()).ok_or("no finite values")?;
    let samples: Vec<u8> = data
        .iter()
        .map(|&value| {
            if value.is_finite() {
                (unit(value, range) * 255.0).round() as u8
            } else {
                0
            }
        })
        .collect();
    super::raster::write_pgm_file(path, cols, rows, &samples)
}

/// Writes a 2-D array as an 8-bit RGB PPM (P6) raster coloured with the same
/// viridis-like palette as [`plot_heatmap_2d`], row 0 on top and column 0 on
/// the left.
///
/// Finite entries are normalised over their range and mapped through the
/// palette; a non-finite entry becomes light grey (`#dddddd`), matching the
/// SVG heat map. No axes, colour bar or annotations are added.
///
/// # Errors
/// Fails for an empty array, an array without a finite entry, or when the file
/// cannot be written.
pub fn plot_scalar_field_ppm(
    data: &Array2<f64>,
    path: &Path,
) -> Result<(), String> {
    let (rows, cols) = data.dim();
    if rows == 0 || cols == 0 {
        return Err("the array is empty".to_owned());
    }
    let range = extent(data.iter().copied()).ok_or("no finite values")?;
    let samples: Vec<u8> = data
        .iter()
        .flat_map(|&value| {
            let (r, g, b) = if value.is_finite() {
                colormap_rgb(unit(value, range))
            } else {
                (0xdd, 0xdd, 0xdd)
            };
            [r, g, b]
        })
        .collect();
    super::raster::write_ppm_file(path, cols, rows, &samples)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn out(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("rssn_plot_test_{}_{name}.svg", std::process::id()))
    }

    fn read_and_remove(path: &std::path::Path) -> String {
        let s = std::fs::read_to_string(path).unwrap();
        let _r = std::fs::remove_file(path);
        s
    }

    fn polylines(svg: &str) -> usize {
        svg.matches("<polyline").count()
    }

    #[test]
    fn single_function_has_one_polyline() {
        let p = out("fn");
        plot_function_2d(|x| x * x, (-1.0, 1.0), 50, &p).unwrap();
        let svg = read_and_remove(&p);
        assert!(svg.starts_with("<svg"));
        assert_eq!(polylines(&svg), 1);
    }

    #[test]
    fn series_polyline_count_matches() {
        let p = out("series");
        let series: Vec<(String, Vec<(f64, f64)>)> = (0..3)
            .map(|k| {
                (
                    format!("s{k}"),
                    (0..10).map(|i| (f64::from(i), f64::from(i * k))).collect(),
                )
            })
            .collect();
        plot_series_2d(&series, &p).unwrap();
        let svg = read_and_remove(&p);
        assert!(svg.starts_with("<svg"));
        assert_eq!(polylines(&svg), 3);
    }

    #[test]
    fn term_plot_evaluates_graph() {
        let mut g = Graph::new();
        let n = g.parse("x^2 + 1").unwrap();
        let p = out("term");
        plot_term_2d(&mut g, n, "x", (0.0, 3.0), 30, &p).unwrap();
        let svg = read_and_remove(&p);
        assert!(svg.starts_with("<svg"));
        assert_eq!(polylines(&svg), 1);
    }

    #[test]
    fn invalid_inputs_are_errors() {
        let p = out("bad");
        assert!(plot_function_2d(|x| x, (1.0, 1.0), 10, &p).is_err());
        assert!(plot_function_2d(|x| x, (0.0, 1.0), 1, &p).is_err());
        assert!(plot_series_2d(&[], &p).is_err());
        assert!(plot_function_2d(|_| f64::NAN, (0.0, 1.0), 10, &p).is_err());
        assert!(!p.exists());
    }

    #[test]
    fn other_plots_write_svg() {
        let m = Array2::from_shape_fn((4, 5), |(r, c)| (r * c) as f64);
        let cases: Vec<(&str, Result<(), String>)> = vec![
            (
                "vf2",
                plot_vector_field_2d(|x, y| (-y, x), (-1.0, 1.0), (-1.0, 1.0), 5, &out("vf2")),
            ),
            (
                "s3",
                plot_surface_3d(|x, y| x * y, (-1.0, 1.0), (-1.0, 1.0), 6, &out("s3")),
            ),
            ("s2", plot_surface_2d(&m, &out("s2"))),
            (
                "curve",
                plot_parametric_curve_3d(|t| (t.cos(), t.sin(), t), (0.0, 6.0), 40, &out("curve")),
            ),
            (
                "path",
                plot_3d_path_from_points(&[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]], &out("path")),
            ),
            (
                "vf3",
                plot_vector_field_3d(|x, y, z| (y, z, x), [(-1.0, 1.0); 3], 3, &out("vf3")),
            ),
            ("heat", plot_heatmap_2d(&m, &out("heat"))),
        ];
        for (name, r) in cases {
            r.unwrap();
            let svg = read_and_remove(&out(name));
            assert!(svg.starts_with("<svg"), "{name}");
            assert!(svg.trim_end().ends_with("</svg>"), "{name}");
        }
    }
}
