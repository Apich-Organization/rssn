//! Standalone SVG plots.
//!
//! Every function writes one self-contained `.svg` file and needs no
//! plotting dependency. Functions are given as plain Rust closures or data;
//! [`plot_term_2d`] plots a term of the expression graph.
//!
//! 3D plots use one fixed oblique projection (yaw 0.6, pitch 0.5 radians):
//! surfaces are painted back to front as height-coloured quads, curves and
//! vector fields as projected lines inside a wireframe bounding box.

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
        let list: Vec<String> = points.iter().map(|p| format!("{:.2},{:.2}", p.0, p.1)).collect();
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
        let list: Vec<String> = points.iter().map(|p| format!("{:.2},{:.2}", p.0, p.1)).collect();
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
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
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
fn extent(values: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
    let (lo, hi) = values
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
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

/// A viridis-like colour map; `t` is clamped to `[0, 1]`.
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn colormap(t: f64) -> String {
    const STOPS: [(f64, f64, f64); 5] = [
        (68.0, 1.0, 84.0),
        (59.0, 82.0, 139.0),
        (33.0, 145.0, 140.0),
        (94.0, 201.0, 98.0),
        (253.0, 231.0, 37.0),
    ];
    let t = if t.is_finite() { t.clamp(0.0, 1.0) } else { 0.0 };
    let scaled = t * 4.0;
    let i = (scaled.floor() as usize).min(3);
    let f = scaled - i as f64;
    let (a, b) = (STOPS[i], STOPS[i + 1]);
    let mix = |x: f64, y: f64| (x + (y - x) * f).round() as u8;
    format!("#{:02x}{:02x}{:02x}", mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

fn unit(
    v: f64,
    range: (f64, f64),
) -> f64 {
    if range.1 > range.0 { (v - range.0) / (range.1 - range.0) } else { 0.5 }
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
    let points: Vec<(f64, f64)> = linspace(range, samples).into_iter().map(|x| (x, f(x))).collect();
    plot_series_2d(&[("f(x)".to_owned(), points)], path)
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
    plot_function_2d(|x| value_at(graph, x), range, samples, path)
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
        &[to, (base.0 - uy * head * 0.4, base.1 + ux * head * 0.4), (base.0 + uy * head * 0.4, base.1 - ux * head * 0.4)],
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
    let max = field.iter().map(|&(_, _, u, v)| u.hypot(v)).fold(0.0, f64::max);
    let frame = Frame { x: x_range, y: y_range };
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
            offset: (WIDTH / 2.0 - scale * (x0 + x1) / 2.0, HEIGHT / 2.0 + scale * (y0 + y1) / 2.0),
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
        (self.offset.0 + self.scale * x, self.offset.1 - self.scale * y, d)
    }

    fn draw_box(
        &self,
        svg: &mut Svg,
    ) {
        let corner = |i: usize| {
            let pick = |axis: usize| if i >> axis & 1 == 0 { self.ranges[axis].0 } else { self.ranges[axis].1 };
            let (x, y, _) = self.at([pick(0), pick(1), pick(2)]);
            (x, y)
        };
        for i in 0..8usize {
            for axis in 0..3 {
                let j = i | (1 << axis);
                if j != i {
                    svg.line(corner(i), corner(j), if i == 0 { "#222" } else { "#bbb" }, if i == 0 { 1.5 } else { 0.8 });
                }
            }
        }
        for (axis, name) in ["x", "y", "z"].iter().enumerate() {
            let (lo, hi) = (corner(0), corner(1 << axis));
            svg.text((lo.0, lo.1 + 16.0), "middle", &format_tick(self.ranges[axis].0));
            svg.text((hi.0, hi.1 + 16.0), "middle", &format!("{name} = {}", format_tick(self.ranges[axis].1)));
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
    let xr = (xs.first().copied().unwrap_or(0.0), xs.last().copied().unwrap_or(1.0));
    let yr = (ys.first().copied().unwrap_or(0.0), ys.last().copied().unwrap_or(1.0));
    let scene = Scene::new([xr, yr, zr]);
    scene.draw_box(svg);
    let mut quads: Vec<Quad> = Vec::new();
    for (r, pair) in z.windows(2).enumerate() {
        for c in 0..xs.len().saturating_sub(1) {
            let corners = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)];
            let pts: Vec<(f64, f64, f64, f64)> = corners
                .iter()
                .filter_map(|&(i, j)| {
                    let v = if i == r { pair[0][j] } else { pair[1][j] };
                    v.is_finite().then(|| {
                        let (px, py, d) = scene.at([xs[j], ys[i], v]);
                        (px, py, d, v)
                    })
                })
                .collect();
            if pts.len() == 4 {
                let depth = pts.iter().map(|p| p.2).sum::<f64>() / 4.0;
                let mean = pts.iter().map(|p| p.3).sum::<f64>() / 4.0;
                quads.push((depth, pts.iter().map(|p| (p.0, p.1)).collect(), unit(mean, zr)));
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
        svg.rect((x, top + h * f64::from(i) / f64::from(steps)), (10.0, h / f64::from(steps) + 0.5), &colormap(t));
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
    let z: Vec<Vec<f64>> = ys.iter().map(|&y| xs.iter().map(|&x| f(x, y)).collect()).collect();
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
    let axis = |k: usize| extent(points.iter().map(move |p| p[k])).ok_or_else(|| "no finite points".to_owned());
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
    let max = field.iter().map(|(_, d)| d[0].hypot(d[1]).hypot(d[2])).fold(0.0, f64::max);
    let scene = Scene::new(ranges);
    let mut svg = Svg::new();
    scene.draw_box(&mut svg);
    let reach = 0.8 * ranges.iter().map(|r| r.1 - r.0).fold(f64::INFINITY, f64::min) / (grid - 1) as f64;
    let mut arrows: Vec<Arrow> = field
        .iter()
        .filter_map(|(p, d)| {
            let m = d[0].hypot(d[1]).hypot(d[2]);
            if m == 0.0 || max == 0.0 {
                return None;
            }
            let tip = [p[0] + d[0] / max * reach, p[1] + d[1] / max * reach, p[2] + d[2] / max * reach];
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
    let (w, h) = (WIDTH - MARGIN_LEFT - MARGIN_RIGHT - 40.0, HEIGHT - MARGIN_TOP - MARGIN_BOTTOM);
    let (cw, ch) = (w / cols as f64, h / rows as f64);
    let mut svg = Svg::new();
    for ((r, c), &v) in data.indexed_iter() {
        let fill = if v.is_finite() { colormap(unit(v, range)) } else { "#dddddd".to_owned() };
        svg.rect((MARGIN_LEFT + c as f64 * cw, MARGIN_TOP + r as f64 * ch), (cw, ch), &fill);
    }
    svg.line((MARGIN_LEFT, MARGIN_TOP), (MARGIN_LEFT + w, MARGIN_TOP), "#222", 1.0);
    svg.line((MARGIN_LEFT, MARGIN_TOP + h), (MARGIN_LEFT + w, MARGIN_TOP + h), "#222", 1.0);
    svg.line((MARGIN_LEFT, MARGIN_TOP), (MARGIN_LEFT, MARGIN_TOP + h), "#222", 1.0);
    svg.line((MARGIN_LEFT + w, MARGIN_TOP), (MARGIN_LEFT + w, MARGIN_TOP + h), "#222", 1.0);
    svg.text((MARGIN_LEFT, MARGIN_TOP + h + 16.0), "middle", "0");
    svg.text((MARGIN_LEFT + w, MARGIN_TOP + h + 16.0), "middle", &cols.to_string());
    svg.text((MARGIN_LEFT - 8.0, MARGIN_TOP + 4.0), "end", "0");
    svg.text((MARGIN_LEFT - 8.0, MARGIN_TOP + h + 4.0), "end", &rows.to_string());
    colorbar(&mut svg, range);
    write_svg(svg, path)
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
            .map(|k| (format!("s{k}"), (0..10).map(|i| (f64::from(i), f64::from(i * k))).collect()))
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
            ("vf2", plot_vector_field_2d(|x, y| (-y, x), (-1.0, 1.0), (-1.0, 1.0), 5, &out("vf2"))),
            ("s3", plot_surface_3d(|x, y| x * y, (-1.0, 1.0), (-1.0, 1.0), 6, &out("s3"))),
            ("s2", plot_surface_2d(&m, &out("s2"))),
            ("curve", plot_parametric_curve_3d(|t| (t.cos(), t.sin(), t), (0.0, 6.0), 40, &out("curve"))),
            ("path", plot_3d_path_from_points(&[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]], &out("path"))),
            ("vf3", plot_vector_field_3d(|x, y, z| (y, z, x), [(-1.0, 1.0); 3], 3, &out("vf3"))),
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
