//! Operators with polynomial eigenfunctions (Bochner's classification) and
//! singular endpoints: Hermite, Laguerre and Jacobi type.
//!
//! `κ (a(x) u_xx + b(x) u_x) + c u` with `a` a polynomial of degree at most
//! two, `b` of degree at most one, and no boundary conditions (regularity)
//! on the natural interval:
//!
//! | `a` | domain | weight `w = p/a`, `p = exp ∫ b/a` |
//! |---|---|---|
//! | constant, `b` linear | the line | `exp((b₁x²/2 + b₀x)/a)` (Hermite) |
//! | linear `a₁(x - x₀)` | `(x₀, ∞)` or `(-∞, x₀)` | `|x - x₀|^α e^{(b₁/a₁) x}` (Laguerre) |
//! | quadratic with real roots | between them | `(x - r₁)^α (r₂ - x)^β` (Jacobi) |
//!
//! The endpoint classification is the usual one: the endpoint is
//! limit-circle for the weights above exactly when the exponent satisfies
//! `α > -1` (resp. the Gaussian/exponential decay at infinity makes the
//! polynomials square integrable, the endpoint being limit-point), and the
//! polynomial eigenfunctions are then the solutions bounded (respectively
//! of polynomial growth) there, so that no boundary condition is imposed.
//! The monic eigenpolynomial of degree `n` has the eigenvalue
//! `λ_n = a₂ n (n-1) + b₁ n`; its coefficients follow from the downward
//! recurrence
//! `[a₂k(k-1) + b₁k - λ_n] c_k = -[a₁(k+1)k + b₀(k+1)] c_{k+1} - a₀(k+2)(k+1) c_{k+2}`.
//! The modes up to degree ten are written out, so the data must be a
//! polynomial of at most that degree (as for Legendre expansions).

use super::Conditions;
use super::Problem;
use super::modes::exact_value;
use super::spectral::Axis;
use super::spectral::Engine;
use super::spectral::Item;
use super::spectral::Role;
use super::spectral::Shape;
use super::spectral::Symbol;
use super::spectral::Time;
use super::spectral::evolution_time;
use super::spectral::finish;
use super::spectral::time_items;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::infinity;
use super::util::is_zero_number;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;

/// Highest degree written out.
pub(super) const MAX_DEGREE: i64 = 10;

/// The coefficients `[c₀, c₁, …]` of a polynomial in `x`, if it is one of
/// degree at most `max`.
fn polynomial(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    max: usize,
) -> Option<Vec<NodeId>> {
    let symbol = cx.graph.symbol_of(x)?;
    let zero = cx.graph.int(0);
    let mut out = if cx.graph.depends_on(cx.graph.find(f), symbol) {
        let mut gens = Gens::default();
        let poly = from_term(cx.graph, &mut gens, f, Limits { terms: 256, exponent: 8 })?;
        let g = gens.find(cx.graph, x)?;
        let mut parts = Vec::new();
        for part in &poly.coefficients_in(g) {
            let term = to_term(cx.graph, &gens, part);
            if cx.graph.depends_on(cx.graph.find(term), symbol) {
                return None;
            }
            parts.push(term);
        }
        parts
    } else {
        vec![f]
    };
    while out.len() > max + 1 {
        let last = out.pop()?;
        if !cx.is_zero(last) {
            return None;
        }
    }
    while out.len() < max + 1 {
        out.push(zero);
    }
    Some(out)
}

/// `k · node`.
fn times(
    cx: &mut Cx<'_>,
    node: NodeId,
    k: i64,
) -> NodeId {
    let factor = cx.graph.int(k);
    mul(cx.graph, &[factor, node])
}

/// An exact small fraction for `v`, or the float.
fn snapped(
    cx: &mut Cx<'_>,
    v: f64,
) -> NodeId {
    for q in 1..=24_i64 {
        #[allow(clippy::cast_precision_loss)]
        let scaled = (v * q as f64).round();
        #[allow(clippy::cast_precision_loss)]
        if (v * q as f64 - scaled).abs() < 1e-10 {
            #[allow(clippy::cast_possible_truncation)]
            if let Some(n) = Number::fraction(scaled as i64, q) {
                return cx.graph.num(n);
            }
        }
    }
    cx.graph.num(Number::from(v))
}

/// The natural domain and weight (as nodes) of the normalised operator
/// `â u'' + b̂ u'`, or `None` if it is not of a classical type.
#[allow(clippy::many_single_char_names, clippy::too_many_lines)]
fn classify(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: [f64; 3],
    b: [f64; 2],
) -> Option<(NodeId, NodeId, NodeId)> {
    let eps = 1e-12;
    let oo = infinity(cx)?;
    let minus_oo = neg(cx.graph, oo);
    let exp = |cx: &mut Cx<'_>, e: NodeId| call(cx, "exp", &[e]);
    if a[2].abs() < eps && a[1].abs() < eps {
        // Hermite: a₀ > 0 and b₁ < 0.
        if a[0] <= 0.0 || b[1] >= 0.0 {
            return None;
        }
        let (quad, lin) = (snapped(cx, b[1] / (2.0 * a[0])), snapped(cx, b[0] / a[0]));
        let x2 = powi(cx.graph, x, 2);
        let t1 = mul(cx.graph, &[quad, x2]);
        let t2 = mul(cx.graph, &[lin, x]);
        let e = add(cx.graph, &[t1, t2]);
        let weight = exp(cx, e)?;
        return Some((minus_oo, oo, weight));
    }
    if a[2].abs() < eps {
        // Laguerre: root x₀ of the linear coefficient.
        let x0 = -a[0] / a[1];
        let k = b[1] / a[1];
        let alpha = (b[1] * x0 + b[0]) / a[1] - 1.0;
        if k >= 0.0 || alpha <= -1.0 {
            return None;
        }
        let x0_node = snapped(cx, x0);
        let (gap, lo, hi) = if a[1] > 0.0 { (sub(cx.graph, x, x0_node), x0_node, oo) } else { (sub(cx.graph, x0_node, x), minus_oo, x0_node) };
        let alpha_node = snapped(cx, alpha);
        let power = pow(cx.graph, gap, alpha_node);
        let k_node = snapped(cx, k);
        let e = mul(cx.graph, &[k_node, x]);
        let decay = exp(cx, e)?;
        let weight = mul(cx.graph, &[power, decay]);
        return Some((lo, hi, weight));
    }
    // Jacobi: two real roots, positive between them.
    let disc = a[1] * a[1] - 4.0 * a[2] * a[0];
    if disc <= 0.0 || a[2] >= 0.0 {
        return None;
    }
    let root = disc.sqrt();
    let (r1, r2) = {
        let (p, q) = ((-a[1] - root) / (2.0 * a[2]), (-a[1] + root) / (2.0 * a[2]));
        (p.min(q), p.max(q))
    };
    let slope = |r: f64| 2.0 * a[2] * r + a[1];
    let c1 = (b[1] * r1 + b[0]) / slope(r1);
    let c2 = (b[1] * r2 + b[0]) / slope(r2);
    if c1 <= 0.0 || c2 <= 0.0 {
        return None;
    }
    let (n1, n2) = (snapped(cx, r1), snapped(cx, r2));
    let g1 = sub(cx.graph, x, n1);
    let g2 = sub(cx.graph, n2, x);
    let e1 = snapped(cx, c1 - 1.0);
    let p1 = pow(cx.graph, g1, e1);
    let e2 = snapped(cx, c2 - 1.0);
    let p2 = pow(cx.graph, g2, e2);
    let weight = mul(cx.graph, &[p1, p2]);
    Some((n1, n2, weight))
}

/// Separation of variables for the operators above.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear || p.dimension() != 2 {
        return None;
    }
    let time_index = evolution_time(cx, p)?;
    let j = 1 - time_index;
    let x = p.vars[j];
    let n = 2;
    if conditions.0.iter().any(|c| c.on == j) {
        return None;
    }
    let (ca, cb) = (p.coefficient(cx.graph, &p.unit(j, 2)), p.coefficient(cx.graph, &p.unit(j, 1)));
    let c0 = p.coefficient(cx.graph, &vec![0; n]);
    let (a1, a2) = (p.coefficient(cx.graph, &p.unit(time_index, 1)), p.coefficient(cx.graph, &p.unit(time_index, 2)));
    if cx.is_zero(ca) || !p.constant(cx.graph, c0) || !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
        return None;
    }
    for (index, c) in &p.linear {
        let known = *index == p.unit(j, 1) || *index == p.unit(j, 2) || *index == vec![0; n] || *index == p.unit(time_index, 1) || *index == p.unit(time_index, 2);
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let a_poly = polynomial(cx, ca, x, 2)?;
    let b_poly = polynomial(cx, cb, x, 1)?;
    let value = |cx: &Cx<'_>, node: NodeId| exact_value(cx, node);
    let af: Vec<f64> = a_poly.iter().map(|&c| value(cx, c)).collect::<Option<_>>()?;
    let bf: Vec<f64> = b_poly.iter().map(|&c| value(cx, c)).collect::<Option<_>>()?;
    // At least one coefficient must depend on x (else the box solvers apply).
    if af[1].abs() < 1e-12 && af[2].abs() < 1e-12 && bf[1].abs() < 1e-12 {
        return None;
    }
    // The evolution is well posed for `κ c₂ / d > 0`-normalised operators:
    // `d u_t = -c₂ u_xx - …` (residual form), so `κ = -sign(d)`.
    let time_coefficient = if is_zero_number(cx.graph, a2) { a1 } else { a2 };
    let sign = super::util::sample(cx.graph, time_coefficient, 0)?;
    let kappa = if sign > 0.0 { -1.0_f64 } else { 1.0 };
    let ah = [kappa * af[0], kappa * af[1], kappa * af[2]];
    let bh = [kappa * bf[0], kappa * bf[1]];
    let chosen = classify(cx, x, ah, bh).map(|found| (kappa, found));
    let (kappa, (lo, hi, weight)) = chosen?;
    let kappa_node = cx.graph.int(if kappa > 0.0 { 1 } else { -1 });
    // Normalised coefficients as nodes (exact).
    let scaled = |cx: &mut Cx<'_>, c: NodeId| {
        let m = mul(cx.graph, &[kappa_node, c]);
        cx.simplify(m)
    };
    let (ah0, ah1, ah2) = (scaled(cx, a_poly[0]), scaled(cx, a_poly[1]), scaled(cx, a_poly[2]));
    let (bh0, bh1) = (scaled(cx, b_poly[0]), scaled(cx, b_poly[1]));
    let mut modes = Vec::new();
    for degree in 0..=MAX_DEGREE {
        let lambda = {
            let t1 = times(cx, ah2, degree * (degree - 1));
            let t2 = times(cx, bh1, degree);
            let l = add(cx.graph, &[t1, t2]);
            cx.simplify(l)
        };
        let mut coefficients = vec![cx.graph.int(0); usize::try_from(degree).ok()? + 3];
        coefficients[usize::try_from(degree).ok()?] = cx.graph.int(1);
        for k in (0..degree).rev() {
            let ku = usize::try_from(k).ok()?;
            let denominator = {
                let t1 = times(cx, ah2, k * (k - 1));
                let t2 = times(cx, bh1, k);
                let m = neg(cx.graph, lambda);
                add(cx.graph, &[t1, t2, m])
            };
            let denominator = cx.simplify(denominator);
            if cx.is_zero(denominator) {
                return None;
            }
            let numerator = {
                let u1 = times(cx, ah1, (k + 1) * k);
                let u2 = times(cx, bh0, k + 1);
                let l1 = add(cx.graph, &[u1, u2]);
                let t1 = mul(cx.graph, &[l1, coefficients[ku + 1]]);
                let u3 = times(cx, ah0, (k + 2) * (k + 1));
                let t2 = mul(cx.graph, &[u3, coefficients[ku + 2]]);
                add(cx.graph, &[t1, t2])
            };
            let c = div(cx, numerator, denominator);
            let c = neg(cx.graph, c);
            coefficients[ku] = cx.simplify(c);
        }
        let mut terms = Vec::new();
        for (power, &c) in coefficients.iter().enumerate().take(usize::try_from(degree).ok()? + 1) {
            if cx.is_zero(c) {
                continue;
            }
            let xp = powi(cx.graph, x, i64::try_from(power).ok()?);
            terms.push(mul(cx.graph, &[c, xp]));
        }
        let mode = add(cx.graph, &terms);
        let mode = cx.simplify(mode);
        let square = powi(cx.graph, mode, 2);
        let body = mul(cx.graph, &[weight, square]);
        let norm = defint(cx, body, x, lo, hi)?;
        let norm = cx.simplify(norm);
        let eigen = {
            let m = neg(cx.graph, lambda);
            cx.simplify(m)
        };
        modes.push((mode, eigen, norm));
    }
    let order = if is_zero_number(cx.graph, a2) { u32::from(!is_zero_number(cx.graph, a1)) } else { 2 };
    let time = Time { var: p.vars[time_index], order, a2, a1 };
    let spatial = vec![(vec![1], kappa_node), (vec![0], c0)];
    let axes = vec![Axis { var: x, shape: Shape::Orthogonal { start: lo, stop: hi, weight, modes } }];
    let mut items = time_items(cx, p, conditions, Some(time_index), order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item { role: Role::Source, h: p.source, special: None });
    }
    let engine = Engine { axes, time: Some(time), symbol: Symbol::Cartesian(spatial), face: None, exterior: None };
    finish(cx, p, conditions, &engine, &items, None, None)
}
