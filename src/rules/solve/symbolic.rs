//! Real roots of polynomials whose coefficients contain parameters.
//!
//! The polynomial is first factored over `Q` in the unknown *and* the
//! parameters (so `(x - a)(x - b)` and `x^3 + a x^2 + a x + 1` split into
//! factors that are linear or quadratic in `x`). Each irreducible factor is
//! then solved by the first applicable method:
//!
//! * a polynomial in `x^k` (`k > 1`) through `u = x^k`: the roots of the
//!   polynomial in `u`, then real `k`-th roots (this covers `x^n = c`,
//!   biquadratics and `x^6 + a x^3 + b`);
//! * degree one and two by the linear and quadratic formulas;
//! * degree three by Cardano's formula (real cube roots, valid where the
//!   discriminant is non-negative) together with the trigonometric form
//!   (valid where it is not);
//! * a palindromic quartic through `u = x + 1/x`;
//! * any other quartic by Ferrari's method with a resolvent-cubic root that
//!   is positive (so the square roots in the formula are defined whenever
//!   the quartic has the real roots in question).
//!
//! The formulas are real only on part of the parameter space; the caller
//! keeps every formula that is a root wherever it is defined and discards
//! those that are wrong at a spot check.

use super::product;
use super::reciprocal;
use crate::graph::op::core;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::multifactor;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;

const MAX_DEPTH: usize = 4;

fn cap() -> usize {
    Limits::default().terms
}

/// The real `k`-th root of `u` (`k` odd) as a term: `u^(1/k)` where `u` is
/// known not to be negative, and `u |u|^(1/k - 1)` otherwise.
pub(super) fn real_root(
    graph: &mut Graph,
    u: NodeId,
    k: i64,
) -> Option<NodeId> {
    let exponent = graph.num(Number::fraction(1, k)?);
    let known_sign = graph.number_of(u).map(Number::to_f64);
    let nonnegative = match known_sign {
        | Some(v) => v >= 0.0,
        | None => graph.facts(u).has(Facts::NONNEGATIVE),
    };
    if nonnegative {
        return Some(graph.node(core::POW, &[u, exponent]));
    }
    if known_sign.is_some() {
        let minus_one = graph.int(-1);
        let negated = graph.node(core::MUL, &[minus_one, u]);
        let root = graph.node(core::POW, &[negated, exponent]);
        return Some(graph.node(core::MUL, &[minus_one, root]));
    }
    let abs = graph.ops().lookup("abs")?;
    let magnitude = graph.node(abs, &[u]);
    let shifted = graph.num(Number::fraction(1 - k, k)?);
    let scale = graph.node(core::POW, &[magnitude, shifted]);
    Some(product(graph, &[u, scale]))
}

/// Roots of `poly` as a polynomial in the generator `g` (all other
/// generators are parameters).
pub(super) fn roots(
    graph: &mut Graph,
    gens: &Gens,
    poly: &Poly,
    g: u32,
    depth: usize,
) -> Option<Vec<NodeId>> {
    if depth > MAX_DEPTH {
        return None;
    }
    let vars = poly.support();
    let factored = if super::factorable(poly) { multifactor::factor(poly, &vars) } else { None };
    let pieces: Vec<Poly> = match factored {
        | Some((_, factors)) => factors.into_iter().map(|(f, _)| f).collect(),
        | None => vec![poly.clone()],
    };
    let mut out = Vec::new();
    for piece in pieces {
        if piece.degree_in(g) == 0 {
            continue;
        }
        out.extend(irreducible(graph, gens, &piece, g, depth)?);
    }
    Some(out)
}

fn irreducible(
    graph: &mut Graph,
    gens: &Gens,
    poly: &Poly,
    g: u32,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let coefficients = poly.coefficients_in(g);
    let degree = coefficients.len().checked_sub(1)?;
    // Powers of the generator common to all terms: x = 0 is a root.
    let low = coefficients.iter().position(|c| !c.is_zero())?;
    let mut out = Vec::new();
    if low > 0 {
        out.push(graph.int(0));
    }
    let coefficients = coefficients.get(low..)?;
    let degree = degree - low;
    if degree == 0 {
        return Some(out);
    }
    // A polynomial in g^k.
    let mut k = 0_usize;
    for (i, c) in coefficients.iter().enumerate() {
        if i > 0 && !c.is_zero() {
            k = num_integer::Integer::gcd(&k, &i);
        }
    }
    if k > 1 {
        let reduced: Vec<Poly> = coefficients.iter().step_by(k).cloned().collect();
        let mut u_poly = Poly::zero();
        for (j, c) in reduced.iter().enumerate() {
            let power = Poly::generator(g).pow(u32::try_from(j).ok()?, cap())?;
            u_poly = u_poly.add(&c.mul(&power, cap())?);
        }
        let k_i = i64::try_from(k).ok()?;
        for r in roots(graph, gens, &u_poly, g, depth + 1)? {
            if k % 2 == 1 {
                out.push(real_root(graph, r, k_i)?);
            } else {
                let exponent = graph.num(Number::fraction(1, k_i)?);
                let magnitude = graph.node(core::POW, &[r, exponent]);
                let minus_one = graph.int(-1);
                out.push(graph.node(core::MUL, &[minus_one, magnitude]));
                out.push(magnitude);
            }
        }
        return Some(out);
    }
    let terms: Vec<NodeId> = coefficients.iter().map(|c| to_term(graph, gens, c)).collect();
    match degree {
        | 1 => {
            let minus_one = graph.int(-1);
            let inverse = reciprocal(graph, *terms.get(1)?);
            out.push(product(graph, &[minus_one, *terms.first()?, inverse]));
        },
        | 2 => out.extend(quadratic(graph, gens, coefficients)?),
        | 3 => out.extend(cubic(graph, gens, coefficients)?),
        | 4 => {
            let found = match palindromic_quartic(graph, gens, coefficients, depth) {
                | Some(found) => found,
                | None => ferrari(graph, gens, coefficients)?,
            };
            out.extend(found);
        },
        | _ => return None,
    }
    Some(out)
}

fn quadratic(
    graph: &mut Graph,
    gens: &Gens,
    c: &[Poly],
) -> Option<Vec<NodeId>> {
    let (c0, b, a) = (c.first()?, c.get(1)?, c.get(2)?);
    let four_ac = a.mul(c0, cap())?.scale(&Number::from(4));
    let discriminant_poly = b.mul(b, cap())?.sub(&four_ac);
    // A negative constant discriminant: no real roots.
    if discriminant_poly.as_constant().is_some_and(|d| d.to_f64() < 0.0) {
        return Some(Vec::new());
    }
    let discriminant = to_term(graph, gens, &discriminant_poly);
    let half = graph.num(Number::fraction(1, 2)?);
    let root = graph.node(core::POW, &[discriminant, half]);
    let neg_b = to_term(graph, gens, &b.neg());
    let two_a = to_term(graph, gens, &a.scale(&Number::from(2)));
    let inverse = reciprocal(graph, two_a);
    let mut out = Vec::new();
    for sign in [-1, 1] {
        let s = graph.int(sign);
        let signed_root = product(graph, &[s, root]);
        let numerator = graph.node(core::ADD, &[neg_b, signed_root]);
        out.push(product(graph, &[numerator, inverse]));
    }
    Some(out)
}

/// `a x^3 + b x^2 + c x + d`: with `x = t - b/(3a)` and `t^3 + p t + q`,
/// `p = (3ac - b^2)/(3a^2)`, `q = (2b^3 - 9abc + 27a^2 d)/(27a^3)`.
fn cubic(
    graph: &mut Graph,
    gens: &Gens,
    c: &[Poly],
) -> Option<Vec<NodeId>> {
    let (d, cc, b, a) = (c.first()?, c.get(1)?, c.get(2)?, c.get(3)?);
    let n = |v: i64| Number::from(v);
    let mul = |x: &Poly, y: &Poly| x.mul(y, cap());
    let a2 = mul(a, a)?;
    let b2 = mul(b, b)?;
    let p_num = mul(a, cc)?.scale(&n(3)).sub(&b2);
    let q_num = mul(&b2, b)?
        .scale(&n(2))
        .sub(&mul(&mul(a, b)?, cc)?.scale(&n(9)))
        .add(&mul(&a2, d)?.scale(&n(27)));
    let a3 = mul(&a2, a)?;
    let p_term = {
        let top = to_term(graph, gens, &p_num);
        let bottom = to_term(graph, gens, &a2.scale(&n(3)));
        let inverse = reciprocal(graph, bottom);
        product(graph, &[top, inverse])
    };
    let q_term = {
        let top = to_term(graph, gens, &q_num);
        let bottom = to_term(graph, gens, &a3.scale(&n(27)));
        let inverse = reciprocal(graph, bottom);
        product(graph, &[top, inverse])
    };
    let shift = {
        let top = to_term(graph, gens, &b.neg());
        let bottom = to_term(graph, gens, &a.scale(&n(3)));
        let inverse = reciprocal(graph, bottom);
        product(graph, &[top, inverse])
    };
    let (two, three, minus_one) = (graph.int(2), graph.int(3), graph.int(-1));
    let half = graph.num(Number::fraction(1, 2)?);
    let third = graph.num(Number::fraction(1, 3)?);
    // Cardano: t = cbrt(-q/2 + s) + cbrt(-q/2 - s), s = sqrt(q^2/4 + p^3/27).
    let q_half = product(graph, &[half, q_term]);
    let neg_q_half = product(graph, &[minus_one, q_half]);
    let q_sq = graph.node(core::POW, &[q_term, two]);
    let p_cube = graph.node(core::POW, &[p_term, three]);
    let quarter = graph.num(Number::fraction(1, 4)?);
    let twenty_seventh = graph.num(Number::fraction(1, 27)?);
    let q_part = product(graph, &[quarter, q_sq]);
    let p_part = product(graph, &[twenty_seventh, p_cube]);
    let discriminant = graph.node(core::ADD, &[q_part, p_part]);
    let s = graph.node(core::POW, &[discriminant, half]);
    let neg_s = product(graph, &[minus_one, s]);
    let u1 = graph.node(core::ADD, &[neg_q_half, s]);
    let u2 = graph.node(core::ADD, &[neg_q_half, neg_s]);
    if p_num.is_zero() {
        // t^3 = -q: one real root, and the trigonometric form would be a
        // spurious 0 * undefined.
        let minus_q = product(graph, &[minus_one, q_term]);
        let t = real_root(graph, minus_q, 3)?;
        return Some(vec![graph.node(core::ADD, &[t, shift])]);
    }
    let (r1, r2) = (real_root(graph, u1, 3)?, real_root(graph, u2, 3)?);
    let t = graph.node(core::ADD, &[r1, r2]);
    let mut out = vec![graph.node(core::ADD, &[t, shift])];
    // Trigonometric form: t_k = 2 sqrt(-p/3) cos(acos((3q/2p) sqrt(-3/p))/3 - 2 pi k/3).
    let (acos, cos, pi) = (graph.ops().lookup("acos")?, graph.ops().lookup("cos")?, graph.ops().lookup("pi")?);
    let minus_p_third = product(graph, &[minus_one, third, p_term]);
    let amplitude_root = graph.node(core::POW, &[minus_p_third, half]);
    let amplitude = product(graph, &[two, amplitude_root]);
    let inverse_p = reciprocal(graph, p_term);
    let minus_three_over_p = product(graph, &[minus_one, three, inverse_p]);
    let sqrt_part = graph.node(core::POW, &[minus_three_over_p, half]);
    let three_halves = graph.num(Number::fraction(3, 2)?);
    let argument = product(graph, &[three_halves, q_term, inverse_p, sqrt_part]);
    let angle = graph.node(acos, &[argument]);
    let base_angle = product(graph, &[third, angle]);
    let pi_node = graph.node(pi, &[]);
    for k in 0..3_i64 {
        let offset = graph.num(Number::fraction(-2 * k, 3)?);
        let shift_angle = product(graph, &[offset, pi_node]);
        let total = graph.node(core::ADD, &[base_angle, shift_angle]);
        let wave = graph.node(cos, &[total]);
        let t = product(graph, &[amplitude, wave]);
        out.push(graph.node(core::ADD, &[t, shift]));
    }
    Some(out)
}

/// `a x^4 + b x^3 + c x^2 + b x + a` through `u = x + 1/x`:
/// `a u^2 + b u + (c - 2a) = 0`, then `x^2 - u x + 1 = 0`.
fn palindromic_quartic(
    graph: &mut Graph,
    gens: &Gens,
    c: &[Poly],
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (e, d, cc, b, a) = (c.first()?, c.get(1)?, c.get(2)?, c.get(3)?, c.get(4)?);
    if e != a || d != b {
        return None;
    }
    let _ = depth;
    let constant = cc.sub(&a.scale(&Number::from(2)));
    let quadratic_u = quadratic(graph, gens, &[constant, b.clone(), a.clone()])?;
    let (two, minus_one, four) = (graph.int(2), graph.int(-1), graph.int(-4));
    let half = graph.num(Number::fraction(1, 2)?);
    let mut out = Vec::new();
    for u in quadratic_u {
        let u_sq = graph.node(core::POW, &[u, two]);
        let disc = graph.node(core::ADD, &[u_sq, four]);
        let root = graph.node(core::POW, &[disc, half]);
        for sign in [minus_one, graph.int(1)] {
            let signed = product(graph, &[sign, root]);
            let top = graph.node(core::ADD, &[u, signed]);
            out.push(product(graph, &[half, top]));
        }
    }
    Some(out)
}

/// Ferrari/Descartes for `a x^4 + b x^3 + c x^2 + d x + e` (ascending
/// coefficients in `c`). With `x = y - b/(4a)` the depressed quartic
/// `y^4 + p y^2 + q y + r` factors as
/// `(y^2 - s y + p/2 + m + q/(2s)) (y^2 + s y + p/2 + m - q/(2s))`
/// where `s = sqrt(2m)` and `m` is a root of the resolvent cubic
/// `8 m^3 + 8 p m^2 + (2 p^2 - 8 r) m - q^2`. That cubic is negative at 0
/// and positive at infinity, so its largest real root is positive and
/// `s` is real; the largest root is Cardano's expression when the cubic
/// has one real root and the trigonometric `k = 0` branch when it has
/// three, and both are offered (each is defined only in its own regime).
fn ferrari(
    graph: &mut Graph,
    gens: &Gens,
    c: &[Poly],
) -> Option<Vec<NodeId>> {
    let (e, d, cc, b, a) = (c.first()?, c.get(1)?, c.get(2)?, c.get(3)?, c.get(4)?);
    let n = |v: i64| Number::from(v);
    let mul = |x: &Poly, y: &Poly| x.mul(y, cap());
    let a2 = mul(a, a)?;
    let a3 = mul(&a2, a)?;
    let a4 = mul(&a3, a)?;
    let b2 = mul(b, b)?;
    let b3 = mul(&b2, b)?;
    let pn = mul(a, cc)?.scale(&n(8)).sub(&b2.scale(&n(3)));
    let qn = b3
        .sub(&mul(&mul(a, b)?, cc)?.scale(&n(4)))
        .add(&mul(&a2, d)?.scale(&n(8)));
    let rn = mul(&b2, &b2)?
        .scale(&n(-3))
        .add(&mul(&a3, e)?.scale(&n(256)))
        .sub(&mul(&mul(&a2, b)?, d)?.scale(&n(64)))
        .add(&mul(&mul(a, &b2)?, cc)?.scale(&n(16)));
    // m = nu / (8 a^2) with nu^3 + P nu^2 + (P^2 - R)/4 nu - Q^2 = 0.
    let quarter = Number::fraction(1, 4)?;
    let linear = mul(&pn, &pn)?.sub(&rn).scale(&quarter);
    let minus_q2 = mul(&qn, &qn)?.neg();
    let resolvent = cubic(graph, gens, &[minus_q2, linear, pn.clone(), Poly::constant(Number::from(1))])?;
    let nus: Vec<NodeId> = resolvent.iter().take(2).copied().collect();
    let (half, minus_one, two) = (graph.num(Number::fraction(1, 2)?), graph.int(-1), graph.int(2));
    let ratio_to = |graph: &mut Graph, top: &Poly, bottom: &Poly| {
        let top = to_term(graph, gens, top);
        let bottom = to_term(graph, gens, bottom);
        let inverse = reciprocal(graph, bottom);
        product(graph, &[top, inverse])
    };
    let p_term = ratio_to(graph, &pn, &a2.scale(&n(8)));
    let q_term = ratio_to(graph, &qn, &a3.scale(&n(8)));
    let shift = ratio_to(graph, &b.neg(), &a.scale(&n(4)));
    let inv_8a2 = {
        let bottom = to_term(graph, gens, &a2.scale(&n(8)));
        reciprocal(graph, bottom)
    };
    let _ = a4;
    let mut out = Vec::new();
    for nu in nus {
        let m = product(graph, &[nu, inv_8a2]);
        let two_m = product(graph, &[two, m]);
        let s = graph.node(core::POW, &[two_m, half]);
        let inv_s = reciprocal(graph, s);
        for sign in [-1_i64, 1] {
            let sg = graph.int(sign);
            let lead = product(graph, &[half, sg, s]);
            // radicand = -2m - 2p - 2 sign q / s
            let minus_two = graph.int(-2);
            let t1 = product(graph, &[minus_two, m]);
            let t2 = product(graph, &[minus_two, p_term]);
            let t3 = product(graph, &[minus_two, sg, q_term, inv_s]);
            let radicand = graph.node(core::ADD, &[t1, t2, t3]);
            let root = graph.node(core::POW, &[radicand, half]);
            for eps in [-1_i64, 1] {
                let ep = graph.int(eps);
                let tail = product(graph, &[half, ep, root]);
                let y = graph.node(core::ADD, &[lead, tail]);
                out.push(graph.node(core::ADD, &[y, shift]));
            }
        }
    }
    let _ = minus_one;
    Some(out)
}
