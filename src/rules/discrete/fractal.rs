//! Fractals and chaos: escape-time sets, one- and two-dimensional maps,
//! strange attractors, iterated function systems and fractal dimensions.
//!
//! Symbolic operators take formulas in a variable (`f` in `x`) and use the
//! calculus and solve rule sets (`diff`, `solve`); the numeric ones wrap
//! [`crate::kernels::fractal_geometry_and_chaos`] and need literal numbers.
//! Everything is deterministic: where the legacy code had a random start or
//! chaos game, the seed is an argument.
//!
//! | operator | value |
//! |---|---|
//! | `mandelbrot_escape(cr, ci, max_iter)`, `julia_escape(zr, zi, cr, ci, max_iter)` | escape time of `z^2 + c` (`max_iter` inside the set) |
//! | `burning_ship_escape(cr, ci, max_iter)`, `multibrot_escape(d, cr, ci, max_iter)` | escape times of the burning ship and of `z^d + c` |
//! | `newton_fractal_root(coeffs, zr, zi, max_iter)` | Newton's method for the polynomial with coefficients (highest degree first) from `zr + i zi`: `list(re, im, iterations)` of the root reached, `false` without convergence |
//! | `mandelbrot_iterate(z, c)`, `mandelbrot_orbit(c, n)` | `z^2 + c`; the symbolic orbit `0, c, c^2 + c, ...` of length `n + 1` (`c` may be symbolic) |
//! | `mandelbrot_fixed_points(c)`, `mandelbrot_stability(z)` | the two (complex) fixed points `(1 -+ sqrt(1 - 4 c)) / 2`; `|2 z|`, the multiplier of `z^2 + c` |
//! | `complex_map_fixed_points(f, z)`, `complex_map_stability(f, z, z0)` | real solutions of `f(z) = z`; `abs(f'(z0))` |
//! | `map_fixed_points(f, x)`, `map_stability(f, x, x0)` | solutions of `f(x) = x`; `f'(x0)` |
//! | `lyapunov_exponent(f, x, x0, n)` | `(1/n) sum ln abs(f'(x_k))` along the orbit of `x0`: a float when `x0` and the map are numeric, else the symbolic average |
//! | `logistic_iterate(r, x0, n)`, `logistic_bifurcation(r_min, r_max, steps, transient, keep)`, `logistic_lyapunov(r, x0, n)` | logistic map `r x (1 - x)`: orbit of length `n + 1`; list of `list(r, x)`; Lyapunov exponent after a transient of 100 steps |
//! | `lorenz()`, `lorenz(sigma, rho, beta)` | the vector field `list(sigma (y - x), x (rho - z) - y, x y - beta z)` in `x`, `y`, `z` (parameters optional) |
//! | `lorenz_orbit(p0, dt, n[, sigma, rho, beta])`, `rossler_orbit(p0, dt, n, a, b, c)` | integrated trajectories, lists of `list(x, y, z)` |
//! | `henon_orbit(p0, n, a, b)`, `tinkerbell_orbit(p0, n, a, b, c, d)` | planar map orbits |
//! | `lorenz_lyapunov(p0, dt, n[, sigma, rho, beta])` | largest Lyapunov exponent of the Lorenz system |
//! | `ifs_apply(maps, vars, point)` | the images of `point` under every map (a list of coordinate formulas in `vars`) |
//! | `ifs_generate(name_or_maps, n, seed)` | chaos game: `sierpinski`, `barnsley` or a list of `list(a, b, c, d, e, f)` affine maps (`x' = a x + b y + e`, `y' = c x + d y + f`, equal probabilities) |
//! | `similarity_dimension(list(r1, ...))`, `moran_dimension(list(r1, ...))` | `-ln(N) / ln(r)` for equal ratios, else the Moran equation `sum r_i^D = 1` in the symbol `D`; its root by bisection |
//! | `box_counting(points, k)`, `box_counting(points, list(sizes))` | box-counting dimension with `k` automatic scales, or the log-log slope over the given box sizes |
//! | `correlation_dimension(points, k)`, `correlation_dimension(points, list(radii))` | correlation dimension, likewise |
//! | `orbit_density(points, nx, ny, xmin, xmax, ymin, ymax)`, `orbit_entropy(density)` | occupancy histogram and its Shannon entropy |

use num_complex::Complex;

use super::apply;
use super::def;
use super::def_request;
use super::div;
use super::float;
use super::idx;
use super::items;
use super::neg;
use super::pow;
use super::prod;
use super::rows;
use super::sum;
use super::V;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::RuleError;
use crate::kernels::fractal_geometry_and_chaos as fg;

fn arg<const N: usize>(a: &[NodeId]) -> Option<[NodeId; N]> {
    <[NodeId; N]>::try_from(a).ok()
}

fn iter_count(
    g: &Graph,
    n: NodeId,
) -> Option<u32> {
    u32::try_from(idx(g, n)?).ok()
}

fn f64s<const N: usize>(
    g: &Graph,
    a: [NodeId; N],
) -> Option<[f64; N]> {
    let mut out = [0.0; N];
    for (slot, n) in out.iter_mut().zip(a) {
        *slot = float(g, n)?;
    }
    Some(out)
}

fn pairs(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<(f64, f64)>> {
    rows(g, n)?
        .into_iter()
        .map(|r| match r.as_slice() {
            | [x, y] => Some((float(g, *x)?, float(g, *y)?)),
            | _ => None,
        })
        .collect()
}

fn triple(
    g: &Graph,
    n: NodeId,
) -> Option<(f64, f64, f64)> {
    match items(g, n)?.as_slice() {
        | [x, y, z] => Some((float(g, *x)?, float(g, *y)?, float(g, *z)?)),
        | _ => None,
    }
}

fn pair(
    g: &Graph,
    n: NodeId,
) -> Option<(f64, f64)> {
    match items(g, n)?.as_slice() {
        | [x, y] => Some((float(g, *x)?, float(g, *y)?)),
        | _ => None,
    }
}

fn point_list(points: &[(f64, f64)]) -> V {
    V::List(points.iter().map(|&(x, y)| V::List(vec![V::Float(x), V::Float(y)])).collect())
}

fn triple_list(points: &[(f64, f64, f64)]) -> V {
    V::List(points.iter().map(|&(x, y, z)| V::List(vec![V::Float(x), V::Float(y), V::Float(z)])).collect())
}

// ----------------------------------------------------------------------
// Escape-time sets
// ----------------------------------------------------------------------

fn mandelbrot_escape(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cr, ci, n] = arg(a)?;
    let g = &*cx.graph;
    Some(V::Int(fg::mandelbrot_escape_time(float(g, cr)?, float(g, ci)?, iter_count(g, n)?).into()))
}

fn julia_escape(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [zr, zi, cr, ci, n] = arg(a)?;
    let g = &*cx.graph;
    let [zr, zi, cr, ci] = f64s(g, [zr, zi, cr, ci])?;
    Some(V::Int(fg::julia_escape_time(zr, zi, cr, ci, iter_count(g, n)?).into()))
}

fn burning_ship_escape(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cr, ci, n] = arg(a)?;
    let g = &*cx.graph;
    let [cx0, cy] = f64s(g, [cr, ci])?;
    let max = iter_count(g, n)?;
    let (mut zx, mut zy, mut k) = (0.0_f64, 0.0_f64, 0);
    while zx * zx + zy * zy <= 4.0 && k < max {
        let next = zx * zx - zy * zy + cx0;
        zy = 2.0 * zx.abs() * zy.abs() + cy;
        zx = next;
        k += 1;
    }
    Some(V::Int(k.into()))
}

fn multibrot_escape(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [d, cr, ci, n] = arg(a)?;
    let g = &*cx.graph;
    let [d, cr, ci] = f64s(g, [d, cr, ci])?;
    if d <= 1.0 {
        return None;
    }
    let max = iter_count(g, n)?;
    let c = Complex::new(cr, ci);
    let radius = 2.0_f64.max(c.norm().powf(1.0 / (d - 1.0)));
    let (mut z, mut k) = (Complex::new(0.0, 0.0), 0);
    while z.norm() <= radius && k < max {
        z = z.powf(d) + c;
        k += 1;
    }
    Some(V::Int(k.into()))
}

fn newton_fractal_root(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [coeffs, zr, zi, n] = arg(a)?;
    let g = &*cx.graph;
    let p: Vec<f64> = items(g, coeffs)?.into_iter().map(|c| float(g, c)).collect::<Option<_>>()?;
    if p.len() < 2 || p[0] == 0.0 {
        return None;
    }
    let [zr, zi] = f64s(g, [zr, zi])?;
    let max = iter_count(g, n)?;
    let degree = p.len() - 1;
    let mut z = Complex::new(zr, zi);
    for k in 0..max {
        // Horner for p and p'.
        let (mut value, mut slope) = (Complex::new(p[0], 0.0), Complex::new(0.0, 0.0));
        for &c in &p[1..] {
            slope = slope * z + value;
            value = value * z + c;
        }
        if slope.norm() < 1e-14 {
            return Some(V::Bool(false));
        }
        let step = value / slope;
        z -= step;
        if step.norm() < 1e-12 * (1.0 + z.norm()) {
            let _ = degree;
            return Some(V::List(vec![V::Float(z.re), V::Float(z.im), V::Int((k + 1).into())]));
        }
    }
    Some(V::Bool(false))
}

// ----------------------------------------------------------------------
// Symbolic dynamics
// ----------------------------------------------------------------------

fn add(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let s = sum(cx.graph, &[a, b]);
    cx.simplify(s)
}

/// `term` with `var` replaced by `value`, simplified.
fn subst(
    cx: &mut Cx<'_>,
    term: NodeId,
    var: NodeId,
    value: NodeId,
) -> NodeId {
    let replaced = cx.graph.substitute(term, var, value);
    cx.simplify(replaced)
}

/// `f'`, simplified (the derivative is taken by the calculus rule set).
fn derivative(
    cx: &mut Cx<'_>,
    f: NodeId,
    var: NodeId,
) -> Option<NodeId> {
    let d = apply(cx.graph, "diff", &[f, var])?;
    Some(cx.simplify(d))
}

fn mandelbrot_step(
    cx: &mut Cx<'_>,
    z: NodeId,
    c: NodeId,
) -> NodeId {
    let two = cx.graph.int(2);
    let sq = pow(cx.graph, z, two);
    add(cx, sq, c)
}

fn mandelbrot_iterate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [z, c] = arg(a)?;
    Some(V::Node(mandelbrot_step(cx, z, c)))
}

fn mandelbrot_orbit(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [c, n] = arg(a)?;
    let n = idx(cx.graph, n)?;
    let mut z = cx.graph.int(0);
    let mut out = vec![z];
    for _ in 0..n {
        z = mandelbrot_step(cx, z, c);
        out.push(z);
    }
    Some(V::nodes(&out))
}

fn mandelbrot_fixed_points(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [c] = arg(a)?;
    let one = cx.graph.int(1);
    let four_c = {
        let four = cx.graph.int(4);
        prod(cx.graph, &[four, c])
    };
    let nfc = neg(cx.graph, four_c);
    let disc = add(cx, one, nfc);
    let root = apply(cx.graph, "sqrt", &[disc])?;
    let root = cx.simplify(root);
    let half = {
        let two = cx.graph.int(2);
        div(cx.graph, one, two)
    };
    let minus_root = neg(cx.graph, root);
    let low = sum(cx.graph, &[one, minus_root]);
    let high = sum(cx.graph, &[one, root]);
    let low = prod(cx.graph, &[half, low]);
    let high = prod(cx.graph, &[half, high]);
    let (low, high) = (cx.simplify(low), cx.simplify(high));
    Some(V::nodes(&[low, high]))
}

fn mandelbrot_stability(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [z] = arg(a)?;
    let two = cx.graph.int(2);
    let m = prod(cx.graph, &[two, z]);
    let abs = apply(cx.graph, "abs", &[m])?;
    Some(V::Node(cx.simplify(abs)))
}

/// The request `solve(f = x, x)`.
fn fixed_points(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, x] = arg(a)?;
    cx.graph.as_symbol(x)?;
    let eq = cx.graph.node(core::EQ, &[f, x]);
    Some(V::Node(apply(cx.graph, "solve", &[eq, x])?))
}

/// `f'(point)`, optionally under `abs`.
fn stability(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    absolute: bool,
) -> Option<V> {
    let [f, x, p] = arg(a)?;
    cx.graph.as_symbol(x)?;
    let d = derivative(cx, f, x)?;
    let v = subst(cx, d, x, p);
    if absolute {
        let abs = apply(cx.graph, "abs", &[v])?;
        return Some(V::Node(cx.simplify(abs)));
    }
    Some(V::Node(v))
}

/// `(1/n) sum ln abs(f'(x_k))`, evaluated in floating point.
fn numeric_lyapunov(
    cx: &mut Cx<'_>,
    f: NodeId,
    df: NodeId,
    x: NodeId,
    x0: f64,
    n: usize,
) -> Option<f64> {
    let symbol = cx.graph.as_symbol(x)?;
    let mut env = Env::numeric(0.0);
    let mut cur = x0;
    let mut total = 0.0;
    for _ in 0..n {
        env.bind(symbol, cur);
        total += cx.graph.eval(df, &env)?.abs().ln();
        cur = cx.graph.eval(f, &env)?;
    }
    let value = total / n as f64;
    value.is_finite().then_some(value)
}

fn lyapunov_exponent(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, x, x0, n] = arg(a)?;
    cx.graph.as_symbol(x)?;
    let n = idx(cx.graph, n)?;
    if n == 0 {
        return None;
    }
    let df = derivative(cx, f, x)?;
    if let Some(start) = float(cx.graph, x0) {
        if let Some(v) = numeric_lyapunov(cx, f, df, x, start, n) {
            return Some(V::Float(v));
        }
    }
    // Symbolic average of ln|f'| along the orbit.
    let mut cur = x0;
    let mut terms = Vec::new();
    for _ in 0..n {
        let slope = subst(cx, df, x, cur);
        let abs = apply(cx.graph, "abs", &[slope])?;
        let ln = apply(cx.graph, "ln", &[abs])?;
        terms.push(cx.simplify(ln));
        cur = subst(cx, f, x, cur);
    }
    let total = sum(cx.graph, &terms);
    let count = cx.graph.int(i64::try_from(n).ok()?);
    let avg = div(cx.graph, total, count);
    Some(V::Node(cx.simplify(avg)))
}

fn lorenz(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (sigma, rho, beta) = match a {
        | [] => (cx.graph.sym("sigma"), cx.graph.sym("rho"), cx.graph.sym("beta")),
        | [s, r, b] => (*s, *r, *b),
        | _ => return None,
    };
    let (x, y, z) = (cx.graph.sym("x"), cx.graph.sym("y"), cx.graph.sym("z"));
    let g = &mut *cx.graph;
    let neg_x = neg(g, x);
    let neg_y = neg(g, y);
    let neg_z = neg(g, z);
    let y_minus_x = sum(g, &[y, neg_x]);
    let dx = prod(g, &[sigma, y_minus_x]);
    let rho_minus_z = sum(g, &[rho, neg_z]);
    let xr = prod(g, &[x, rho_minus_z]);
    let dy = sum(g, &[xr, neg_y]);
    let xy = prod(g, &[x, y]);
    let bz = prod(g, &[beta, z]);
    let nbz = neg(g, bz);
    let dz = sum(g, &[xy, nbz]);
    let field: Vec<NodeId> = [dx, dy, dz].iter().map(|&t| cx.simplify(t)).collect();
    Some(V::nodes(&field))
}

// ----------------------------------------------------------------------
// Numeric maps and attractors
// ----------------------------------------------------------------------

fn logistic_iterate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [r, x0, n] = arg(a)?;
    let g = &*cx.graph;
    let orbit = fg::logistic_map_iterate(float(g, x0)?, float(g, r)?, idx(g, n)?);
    Some(V::List(orbit.into_iter().map(V::Float).collect()))
}

fn logistic_bifurcation(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [lo, hi, steps, transient, keep] = arg(a)?;
    let g = &*cx.graph;
    let steps = idx(g, steps)?;
    if steps < 2 {
        return None;
    }
    let data = fg::logistic_bifurcation((float(g, lo)?, float(g, hi)?), steps, idx(g, transient)?, idx(g, keep)?, 0.5);
    Some(point_list(&data))
}

fn logistic_lyapunov(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [r, x0, n] = arg(a)?;
    let g = &*cx.graph;
    let n = idx(g, n)?;
    (n > 0).then(|| V::Float(fg::lyapunov_exponent_logistic(float(g, r).unwrap_or(f64::NAN), float(g, x0).unwrap_or(f64::NAN), 100, n)))
        .filter(|v| matches!(v, V::Float(x) if x.is_finite()))
}

fn lorenz_orbit(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let g = &*cx.graph;
    let (p0, dt, n) = (triple(g, *a.first()?)?, float(g, *a.get(1)?)?, idx(g, *a.get(2)?)?);
    let orbit = match a.len() {
        | 3 => fg::generate_lorenz_attractor(p0, dt, n),
        | 6 => {
            let [s, r, b] = f64s(g, [a[3], a[4], a[5]])?;
            fg::generate_lorenz_attractor_custom(p0, dt, n, s, r, b)
        },
        | _ => return None,
    };
    Some(triple_list(&orbit))
}

fn lorenz_lyapunov(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let g = &*cx.graph;
    let (p0, dt, n) = (triple(g, *a.first()?)?, float(g, *a.get(1)?)?, idx(g, *a.get(2)?)?);
    let [s, r, b] = match a.len() {
        | 3 => [10.0, 28.0, 8.0 / 3.0],
        | 6 => f64s(g, [a[3], a[4], a[5]])?,
        | _ => return None,
    };
    let v = fg::lyapunov_exponent_lorenz(p0, dt, n, s, r, b);
    v.is_finite().then_some(V::Float(v))
}

fn rossler_orbit(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p0, dt, n, pa, pb, pc] = arg(a)?;
    let g = &*cx.graph;
    let [pa, pb, pc] = f64s(g, [pa, pb, pc])?;
    Some(triple_list(&fg::generate_rossler_attractor(triple(g, p0)?, float(g, dt)?, idx(g, n)?, pa, pb, pc)))
}

fn henon_orbit(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p0, n, pa, pb] = arg(a)?;
    let g = &*cx.graph;
    let [pa, pb] = f64s(g, [pa, pb])?;
    Some(point_list(&fg::generate_henon_map(pair(g, p0)?, idx(g, n)?, pa, pb)))
}

fn tinkerbell_orbit(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p0, n, pa, pb, pc, pd] = arg(a)?;
    let g = &*cx.graph;
    let [pa, pb, pc, pd] = f64s(g, [pa, pb, pc, pd])?;
    Some(point_list(&fg::generate_tinkerbell_map(pair(g, p0)?, idx(g, n)?, pa, pb, pc, pd)))
}

// ----------------------------------------------------------------------
// Iterated function systems
// ----------------------------------------------------------------------

fn ifs_apply(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [maps, vars, point] = arg(a)?;
    let vars = items(cx.graph, vars)?;
    let point = items(cx.graph, point)?;
    if vars.len() != point.len() || vars.iter().any(|&v| cx.graph.as_symbol(v).is_none()) {
        return None;
    }
    // Substitute through fresh placeholders so that a point that mentions
    // the variables is not substituted twice.
    let holders: Vec<NodeId> = (0..vars.len()).map(|k| cx.graph.sym(&format!("$ifs{k}"))).collect();
    let mut images = Vec::new();
    for map in items(cx.graph, maps)? {
        let coords = items(cx.graph, map).unwrap_or_else(|| vec![map]);
        let mut image = Vec::new();
        for coord in coords {
            let mut t = coord;
            for (&v, &h) in vars.iter().zip(&holders) {
                t = cx.graph.substitute(t, v, h);
            }
            for (&h, &p) in holders.iter().zip(&point) {
                t = cx.graph.substitute(t, h, p);
            }
            image.push(cx.simplify(t));
        }
        images.push(V::nodes(&image));
    }
    Some(V::List(images))
}

type Maps = Vec<fg::AffineTransform2D>;

/// A named preset or a list of affine maps with equal probabilities.
fn read_ifs(
    g: &Graph,
    n: NodeId,
) -> Option<(Maps, Vec<f64>)> {
    if let Some(l) = items(g, n) {
        let maps: Maps = l
            .into_iter()
            .map(|m| {
                let v: Vec<f64> = items(g, m)?.into_iter().map(|c| float(g, c)).collect::<Option<_>>()?;
                let [a, b, c, d, e, f] = <[f64; 6]>::try_from(v).ok()?;
                Some(fg::AffineTransform2D::new(a, b, c, d, e, f))
            })
            .collect::<Option<_>>()?;
        if maps.is_empty() {
            return None;
        }
        let p = vec![1.0 / maps.len() as f64; maps.len()];
        return Some((maps, p));
    }
    match g.display(n).trim_matches('"') {
        | "sierpinski" => Some(fg::sierpinski_triangle_ifs()),
        | "barnsley" => Some(fg::barnsley_fern_ifs()),
        | _ => None,
    }
}

fn ifs_generate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [spec, n, seed] = arg(a)?;
    let g = &*cx.graph;
    let (maps, probs) = read_ifs(g, spec)?;
    let n = idx(g, n)?;
    let mut state = super::small(g, seed)?.unsigned_abs() ^ 0x9e37_79b9_7f4a_7c15;
    let mut next = move || {
        state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1_u64 << 53) as f64
    };
    let total: f64 = probs.iter().sum();
    let mut p = (next(), next());
    let mut out = Vec::with_capacity(n);
    for k in 0..n + 20 {
        let r = next() * total;
        let mut acc = 0.0;
        let pick = probs
            .iter()
            .position(|&q| {
                acc += q;
                r <= acc
            })
            .unwrap_or(maps.len() - 1);
        p = maps[pick].apply(p);
        if k >= 20 {
            out.push(p);
        }
    }
    Some(point_list(&out))
}

// ----------------------------------------------------------------------
// Dimensions
// ----------------------------------------------------------------------

fn similarity_dimension(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [ratios] = arg(a)?;
    let r = items(cx.graph, ratios)?;
    let first = *r.first()?;
    let g = &mut *cx.graph;
    if r.iter().all(|&x| x == first) {
        // D = ln(N) / ln(1 / r).
        let count = g.int(i64::try_from(r.len()).ok()?);
        let top = apply(g, "ln", &[count])?;
        let inv = super::inv(g, first);
        let inv = cx.simplify(inv);
        let bottom = apply(cx.graph, "ln", &[inv])?;
        let q = div(cx.graph, top, bottom);
        return Some(V::Node(cx.simplify(q)));
    }
    let d = g.sym("D");
    let terms: Vec<NodeId> = r.iter().map(|&x| pow(g, x, d)).collect();
    let total = sum(g, &terms);
    let one = g.int(1);
    let eq = g.node(core::EQ, &[total, one]);
    Some(V::Node(eq))
}

fn moran_dimension(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [ratios] = arg(a)?;
    let g = &*cx.graph;
    let r: Vec<f64> = items(g, ratios)?.into_iter().map(|x| float(g, x)).collect::<Option<_>>()?;
    if r.is_empty() || r.iter().any(|&x| x <= 0.0 || x >= 1.0 || !x.is_finite()) {
        return None;
    }
    let f = |d: f64| r.iter().map(|x| x.powf(d)).sum::<f64>() - 1.0;
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    while f(hi) > 0.0 {
        hi *= 2.0;
        if hi > 1e6 {
            return None;
        }
    }
    for _ in 0..200 {
        let mid = f64::midpoint(lo, hi);
        if f(mid) > 0.0 { lo = mid } else { hi = mid }
    }
    Some(V::Float(f64::midpoint(lo, hi)))
}

/// The least-squares slope of `ys` against `xs`.
fn slope(
    xs: &[f64],
    ys: &[f64],
) -> Option<f64> {
    let n = xs.len() as f64;
    if xs.len() < 2 {
        return None;
    }
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    (sxx > 0.0).then(|| sxy / sxx).filter(|s| s.is_finite())
}

fn box_counting(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [points, scales] = arg(a)?;
    let g = &*cx.graph;
    let pts = pairs(g, points)?;
    if pts.len() < 2 {
        return None;
    }
    if let Some(k) = idx(g, scales) {
        let d = fg::box_counting_dimension(&pts, k);
        return d.is_finite().then_some(V::Float(d));
    }
    let sizes: Vec<f64> = items(g, scales)?.into_iter().map(|s| float(g, s)).collect::<Option<_>>()?;
    if sizes.iter().any(|&s| s <= 0.0) {
        return None;
    }
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for s in sizes {
        let mut boxes: Vec<(i64, i64)> = pts.iter().map(|&(x, y)| ((x / s).floor() as i64, (y / s).floor() as i64)).collect();
        boxes.sort_unstable();
        boxes.dedup();
        xs.push((1.0 / s).ln());
        ys.push((boxes.len() as f64).ln());
    }
    slope(&xs, &ys).map(V::Float)
}

fn correlation_dimension(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [points, radii] = arg(a)?;
    let g = &*cx.graph;
    let pts = pairs(g, points)?;
    if pts.len() < 2 {
        return None;
    }
    if let Some(k) = idx(g, radii) {
        let d = fg::correlation_dimension(&pts, k);
        return d.is_finite().then_some(V::Float(d));
    }
    let radii: Vec<f64> = items(g, radii)?.into_iter().map(|s| float(g, s)).collect::<Option<_>>()?;
    let total = (pts.len() * (pts.len() - 1) / 2) as f64;
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for r in radii {
        if r <= 0.0 {
            return None;
        }
        let mut close = 0_usize;
        for (i, p) in pts.iter().enumerate() {
            close += pts[i + 1..].iter().filter(|q| (p.0 - q.0).hypot(p.1 - q.1) < r).count();
        }
        if close > 0 {
            xs.push(r.ln());
            ys.push((close as f64 / total).ln());
        }
    }
    slope(&xs, &ys).map(V::Float)
}

fn orbit_density(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [points, nx, ny, x0, x1, y0, y1] = arg(a)?;
    let g = &*cx.graph;
    let [x0, x1, y0, y1] = f64s(g, [x0, x1, y0, y1])?;
    let (nx, ny) = (idx(g, nx)?, idx(g, ny)?);
    if nx == 0 || ny == 0 || x1 <= x0 || y1 <= y0 {
        return None;
    }
    let d = fg::orbit_density(&pairs(g, points)?, nx, ny, (x0, x1), (y0, y1));
    Some(V::List(d.into_iter().map(|row| V::List(row.into_iter().map(V::uint).collect())).collect()))
}

fn orbit_entropy(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [density] = arg(a)?;
    let g = &*cx.graph;
    let d: Vec<Vec<usize>> = rows(g, density)?
        .into_iter()
        .map(|r| r.into_iter().map(|c| idx(g, c)).collect::<Option<_>>())
        .collect::<Option<_>>()?;
    Some(V::Float(fg::orbit_entropy(&d)))
}

/// Registers the fractal and chaos operators.
pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "mandelbrot_escape", Arity::Fixed(3), mandelbrot_escape)?;
    def(i, "julia_escape", Arity::Fixed(5), julia_escape)?;
    def(i, "burning_ship_escape", Arity::Fixed(3), burning_ship_escape)?;
    def(i, "multibrot_escape", Arity::Fixed(4), multibrot_escape)?;
    def(i, "newton_fractal_root", Arity::Fixed(4), newton_fractal_root)?;
    def_request(i, "mandelbrot_iterate", Arity::Fixed(2), mandelbrot_iterate)?;
    def(i, "mandelbrot_orbit", Arity::Fixed(2), mandelbrot_orbit)?;
    def(i, "mandelbrot_fixed_points", Arity::Fixed(1), mandelbrot_fixed_points)?;
    def_request(i, "mandelbrot_stability", Arity::Fixed(1), mandelbrot_stability)?;
    def_request(i, "complex_map_fixed_points", Arity::Fixed(2), fixed_points)?;
    def_request(i, "complex_map_stability", Arity::Fixed(3), |cx, a| stability(cx, a, true))?;
    def_request(i, "map_fixed_points", Arity::Fixed(2), fixed_points)?;
    def_request(i, "map_stability", Arity::Fixed(3), |cx, a| stability(cx, a, false))?;
    def_request(i, "lyapunov_exponent", Arity::Fixed(4), lyapunov_exponent)?;
    def(i, "logistic_iterate", Arity::Fixed(3), logistic_iterate)?;
    def(i, "logistic_bifurcation", Arity::Fixed(5), logistic_bifurcation)?;
    def(i, "logistic_lyapunov", Arity::Fixed(3), logistic_lyapunov)?;
    def(i, "lorenz", Arity::Variadic, lorenz)?;
    def(i, "lorenz_orbit", Arity::Variadic, lorenz_orbit)?;
    def(i, "lorenz_lyapunov", Arity::Variadic, lorenz_lyapunov)?;
    def(i, "rossler_orbit", Arity::Fixed(6), rossler_orbit)?;
    def(i, "henon_orbit", Arity::Fixed(4), henon_orbit)?;
    def(i, "tinkerbell_orbit", Arity::Fixed(6), tinkerbell_orbit)?;
    def(i, "ifs_apply", Arity::Fixed(3), ifs_apply)?;
    def(i, "ifs_generate", Arity::Fixed(3), ifs_generate)?;
    def_request(i, "similarity_dimension", Arity::Fixed(1), similarity_dimension)?;
    def(i, "moran_dimension", Arity::Fixed(1), moran_dimension)?;
    def(i, "box_counting", Arity::Fixed(2), box_counting)?;
    def(i, "correlation_dimension", Arity::Fixed(2), correlation_dimension)?;
    def(i, "orbit_density", Arity::Fixed(7), orbit_density)?;
    def(i, "orbit_entropy", Arity::Fixed(1), orbit_entropy)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;

    /// Every number in `text`, in order.
    fn floats(text: &str) -> Vec<f64> {
        text.split(|c: char| !(c.is_ascii_digit() || matches!(c, '.' | '-' | 'e')))
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    fn value(src: &str) -> f64 {
        let v = floats(&s(src));
        assert_eq!(v.len(), 1, "{src}");
        v[0]
    }

    fn close(
        src: &str,
        expected: &[f64],
    ) {
        let got = floats(&s(src));
        assert_eq!(got.len(), expected.len(), "{src}: {got:?}");
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() < 1e-9, "{src}: {g} vs {e}");
        }
    }

    #[test]
    fn escape_times() {
        assert_eq!(s("mandelbrot_escape(0, 0, 50)"), "50");
        assert_eq!(s("mandelbrot_escape(-1, 0, 30)"), "30", "period-2 cycle stays bounded");
        assert_eq!(s("mandelbrot_escape(1, 1, 50)"), "2");
        assert_eq!(s("mandelbrot_escape(2, 0, 50)"), "2");
        assert_eq!(s("mandelbrot_escape(1/4, 0, 100)"), "100");
        assert_eq!(s("mandelbrot_escape(0, 0, -1)"), "mandelbrot_escape(0, 0, -1)");
        assert_eq!(s("mandelbrot_escape(c, 0, 5)"), "mandelbrot_escape(c, 0, 5)");
        assert_eq!(s("julia_escape(0, 0, -0.4, 0.6, 50)"), "26");
        assert_eq!(s("julia_escape(2, 2, 0, 0, 50)"), "0");
        assert_eq!(s("julia_escape(0, 0, 0, 0, 20)"), "20");
        assert_eq!(s("burning_ship_escape(0, 0, 40)"), "40");
        assert_eq!(s("burning_ship_escape(2, 2, 40)"), "1");
        assert_eq!(s("multibrot_escape(3, 0, 0, 40)"), "40");
        assert_eq!(s("multibrot_escape(3, 2, 2, 40)"), "1");
        assert_eq!(s("multibrot_escape(2, 1, 1, 50)"), s("mandelbrot_escape(1, 1, 50)"));
        assert_eq!(s("multibrot_escape(1, 0, 0, 5)"), "multibrot_escape(1, 0, 0, 5)");
    }

    #[test]
    fn newton_fractal_finds_roots() {
        close("newton_fractal_root(list(1, 0, 0, -1), 2, 0, 50)", &[1.0, 0.0, 7.0]);
        close("newton_fractal_root(list(1, 0, 1), 1, 1, 50)", &[0.0, 1.0, 7.0]);
        close("newton_fractal_root(list(1, 0, 1), 1, -1, 50)", &[0.0, -1.0, 7.0]);
        assert_eq!(s("newton_fractal_root(list(1, 0, 1), 1, 0, 50)"), "false", "p'(1) is fine but the real line never reaches a root");
        assert_eq!(s("newton_fractal_root(list(1, 0, -1), 0, 0, 50)"), "false", "p'(0) = 0");
        assert_eq!(s("newton_fractal_root(list(0, 1), 1, 0, 50)"), "newton_fractal_root(list(0, 1), 1, 0, 50)");
    }

    #[test]
    fn mandelbrot_orbits_are_symbolic() {
        assert_eq!(s("mandelbrot_orbit(1, 4)"), "list(0, 1, 2, 5, 26)");
        assert_eq!(s("mandelbrot_orbit(-1, 5)"), "list(0, -1, 0, -1, 0, -1)");
        assert_eq!(s("mandelbrot_orbit(c, 0)"), "list(0)");
        assert_eq!(s("mandelbrot_orbit(c, 3)"), "list(0, c, c^2 + c, (c^2 + c)^2 + c)");
        assert_eq!(s("mandelbrot_iterate(z, c)"), "z^2 + c");
        assert_eq!(s("mandelbrot_iterate(3, 1)"), "10");
    }

    #[test]
    fn mandelbrot_fixed_points_and_stability() {
        assert_eq!(s("mandelbrot_fixed_points(0)"), "list(0, 1)");
        assert_eq!(s("mandelbrot_fixed_points(-2)"), "list(-1, 2)");
        assert_eq!(s("mandelbrot_fixed_points(c)"), "list(1/2 - 1/2*sqrt(1 - 4*c), 1/2*(sqrt(1 - 4*c) + 1))");
        assert_eq!(s("mandelbrot_stability(1/2)"), "1");
        assert_eq!(s("mandelbrot_stability(z)"), "abs(2*z)");
    }

    #[test]
    fn map_fixed_points_and_stability() {
        assert_eq!(s("map_fixed_points(3*x*(1-x), x)"), "list(0, 2/3)");
        assert_eq!(s("map_fixed_points(r*x*(1-x), x)"), "list(0, (r - 1)/r)");
        assert_eq!(s("complex_map_fixed_points(x^2 - 2, x)"), "list(-1, 2)");
        assert_eq!(s("map_stability(r*x*(1-x), x, 0)"), "r");
        assert_eq!(s("map_stability(r*x*(1-x), x, 1 - 1/r)"), "2 - r");
        assert_eq!(s("map_stability(4*x*(1-x), x, 1/2)"), "0");
        assert_eq!(s("complex_map_stability(z^2, z, 1)"), "2");
        assert_eq!(s("complex_map_stability(z^2, z, -3)"), "6");
        assert_eq!(s("complex_map_stability(z^2, z, w)"), "abs(2*w)");
    }

    #[test]
    fn lyapunov_exponents_of_maps() {
        // Numeric orbit.
        assert!((value("lyapunov_exponent(4*x*(1-x), x, 0.3, 200)") - 2.0_f64.ln()).abs() < 0.05);
        assert!((value("lyapunov_exponent(2*x, x, 1/3, 5)") - 2.0_f64.ln()).abs() < 1e-12);
        assert!(value("lyapunov_exponent(x/2, x, 1, 7)") < 0.0);
        // Symbolic average.
        assert_eq!(s("lyapunov_exponent(2*x, x, y, 5)"), "ln(2)");
        assert_eq!(s("lyapunov_exponent(a*x, x, y, 3)"), "ln(abs(a))");
    }

    #[test]
    fn logistic_map() {
        assert_eq!(s("logistic_iterate(2, 0.5, 3)"), "list(0.5, 0.5, 0.5, 0.5)");
        close("logistic_iterate(2.5, 0.25, 2)", &[0.25, 0.46875, 0.62255859375]);
        assert_eq!(floats(&s("logistic_iterate(3.7, 0.2, 10)")).len(), 11);
        // Below r = 3 every start settles on the fixed point 1 - 1/r.
        let b = floats(&s("logistic_bifurcation(2.5, 2.5001, 2, 200, 3)"));
        assert_eq!(b.len(), 12);
        assert!(b.chunks(2).all(|p| (p[1] - 0.6).abs() < 1e-3), "{b:?}");
        assert!(s("logistic_bifurcation(1, 2, 1, 0, 1)").starts_with("logistic_bifurcation("));
        assert!((value("logistic_lyapunov(4, 0.3, 500)") - 2.0_f64.ln()).abs() < 0.1);
        assert!((value("logistic_lyapunov(2.5, 0.3, 500)") - 0.5_f64.ln()).abs() < 1e-6);
    }

    #[test]
    fn lorenz_system_and_attractors() {
        assert_eq!(s("lorenz()"), "list(sigma*(y - x), x*(rho - z) - y, x*y - beta*z)");
        assert_eq!(s("lorenz(10, 28, 8/3)"), "list(10*y - 10*x, x*(28 - z) - y, x*y - 8/3*z)");
        let std = s("lorenz_orbit(list(1, 1, 1), 0.01, 2)");
        assert_eq!(std, s("lorenz_orbit(list(1, 1, 1), 0.01, 2, 10, 28, 8/3)"));
        close("lorenz_orbit(list(1, 1, 1), 0.01, 2)", &[1.0, 1.26, 0.9833333333333333, 1.026, 1.5175666666666667, 0.9697111111111111]);
        close("rossler_orbit(list(1, 1, 1), 0.01, 1, 0.2, 0.2, 5.7)", &[0.98, 1.012, 0.955]);
        close("henon_orbit(list(0, 0), 2, 1.4, 0.3)", &[1.0, 0.0, -0.4, 0.3]);
        assert_eq!(floats(&s("tinkerbell_orbit(list(-0.72, -0.64), 5, 0.9, -0.6013, 2, 0.5)")).len(), 10);
        let l = value("lorenz_lyapunov(list(1, 1, 1), 0.01, 5000)");
        assert!(l > 0.3 && l < 2.0, "{l}");
        assert_eq!(s("lorenz_lyapunov(list(1, 1, 1), 0.01, 100)"), s("lorenz_lyapunov(list(1, 1, 1), 0.01, 100, 10, 28, 8/3)"));
        assert!(s("lorenz_orbit(list(1, 1), 0.01, 2)").starts_with("lorenz_orbit("));
    }

    #[test]
    fn ifs_application_and_chaos_game() {
        assert_eq!(s("ifs_apply(list(list(x/2, y/2), list(x/2 + 1/2, y/2)), list(x, y), list(1, 1))"), "list(list(1/2, 1/2), list(1, 1/2))");
        assert_eq!(s("ifs_apply(list(list(x/2, y/2), list(x/2 + 1/2, y/2)), list(x, y), list(y, x))"), "list(list(1/2*y, 1/2*x), list(1/2*y + 1/2, 1/2*x))");
        assert_eq!(s("ifs_apply(list(2*x, x + 1), list(x), list(3))"), "list(list(6), list(4))");
        assert!(s("ifs_apply(list(list(x, y)), list(x, y), list(1))").starts_with("ifs_apply("));
        // Chaos game points stay inside the attractor's bounding box.
        let tri = floats(&s("ifs_generate(sierpinski, 200, 7)"));
        assert_eq!(tri.len(), 400);
        assert!(tri.chunks(2).all(|p| (0.0..=1.0).contains(&p[0]) && (0.0..=1.0).contains(&p[1]) && p[1] <= 2.0 * p[0].min(1.0 - p[0]) + 1e-9));
        assert_eq!(s("ifs_generate(sierpinski, 10, 7)"), s("ifs_generate(sierpinski, 10, 7)"));
        assert_ne!(s("ifs_generate(sierpinski, 10, 7)"), s("ifs_generate(sierpinski, 10, 8)"));
        let fern = floats(&s("ifs_generate(barnsley, 100, 1)"));
        assert!(fern.chunks(2).all(|p| p[0].abs() < 3.0 && (0.0..10.5).contains(&p[1])));
        let custom = floats(&s("ifs_generate(list(list(0.5, 0, 0, 0.5, 0, 0), list(0.5, 0, 0, 0.5, 0.5, 0.5)), 50, 3)"));
        assert_eq!(custom.len(), 100);
        assert!(s("ifs_generate(unknown, 3, 1)").starts_with("ifs_generate("));
    }

    #[test]
    fn similarity_dimensions() {
        assert_eq!(s("similarity_dimension(list(1/3, 1/3))"), "ln(2)/ln(3)");
        assert_eq!(s("similarity_dimension(list(1/2, 1/2, 1/2))"), "ln(3)/ln(2)");
        assert_eq!(s("similarity_dimension(list(1/2, 1/2, 1/2, 1/2))"), "ln(4)/ln(2)");
        assert_eq!(s("similarity_dimension(list(1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3))"), "ln(8)/ln(3)");
        let eq = s("similarity_dimension(list(1/2, 1/3))");
        assert!(eq.contains('D') && eq.contains('='), "{eq}");
        assert!((value("moran_dimension(list(1/3, 1/3))") - 2.0_f64.ln() / 3.0_f64.ln()).abs() < 1e-12);
        assert!((value("moran_dimension(list(1/2, 1/2, 1/2))") - 3.0_f64.ln() / 2.0_f64.ln()).abs() < 1e-12);
        let d = value("moran_dimension(list(1/2, 1/4))");
        assert!((0.5_f64.powf(d) + 0.25_f64.powf(d) - 1.0).abs() < 1e-12);
        assert_eq!(s("moran_dimension(list(2, 1/2))"), "moran_dimension(list(2, 1/2))");
    }

    fn line(n: usize) -> String {
        let pts: Vec<String> = (0..n).map(|k| format!("list({k}/{n}, {k}/{n})")).collect();
        format!("list({})", pts.join(", "))
    }

    fn grid(n: usize) -> String {
        let pts: Vec<String> = (0..n).flat_map(|i| (0..n).map(move |j| format!("list({i}/{n}, {j}/{n})"))).collect();
        format!("list({})", pts.join(", "))
    }

    #[test]
    fn box_counting_and_correlation_dimensions() {
        let l = line(64);
        assert!((value(&format!("box_counting({l}, list(0.5, 0.25, 0.125, 0.0625))")) - 1.0).abs() < 0.1);
        let g = grid(16);
        assert!((value(&format!("box_counting({g}, list(0.5, 0.25, 0.125))")) - 2.0).abs() < 0.2);
        assert!((value(&format!("correlation_dimension({l}, list(0.1, 0.2, 0.4))")) - 1.0).abs() < 0.3);
        assert!((value(&format!("correlation_dimension({g}, list(0.1, 0.2, 0.3))")) - 2.0).abs() < 0.5);
        // Automatic scales (kernel).
        assert!((value(&format!("box_counting({l}, 4)")) - 1.0).abs() < 0.5);
        assert!(value(&format!("correlation_dimension({g}, 5)")) > 1.0);
        assert!(s("box_counting(list(list(0, 0)), list(1))").starts_with("box_counting("));
        assert!(s("box_counting(list(list(0, 0), list(1, 1)), list(0))").starts_with("box_counting("));
    }

    #[test]
    fn orbit_density_and_entropy() {
        assert_eq!(s("orbit_density(list(list(0.1, 0.1), list(0.9, 0.9), list(0.2, 0.1)), 2, 2, 0, 1, 0, 1)"), "list(list(2, 0), list(0, 1))");
        assert!((value("orbit_entropy(list(list(1, 1), list(1, 1)))") - 4.0_f64.ln()).abs() < 1e-12);
        assert_eq!(value("orbit_entropy(list(list(5, 0), list(0, 0)))"), 0.0);
        assert!(s("orbit_density(list(list(0, 0)), 0, 2, 0, 1, 0, 1)").starts_with("orbit_density("));
    }
}
