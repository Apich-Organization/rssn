//! Sub-solvers tried when the equation is not a polynomial in a single
//! generator of the unknown, each reducing it to a problem the main
//! solver handles; the caller verifies every candidate against the
//! original equation, so these may over-generate.
//!
//! * **Higher-degree polynomials** over `Q`: irreducible cubics by
//!   Cardano's formula (the trigonometric form for three real roots),
//!   biquadratic quartics by `u = x²`, quartics whose Ferrari resolvent
//!   has a rational root by Ferrari's method, and otherwise the real roots
//!   numerically (isolated exactly by Sturm sequences, then refined).
//! * **Absolute values**: both signs of the argument.
//! * **Trigonometric equations** in one argument (after expanding sums
//!   and integer multiples): the Weierstrass substitution
//!   `t = tan(a/2)` makes them rational in `t`; `a = 2 atan(t)` and the
//!   point `a = π` give the solutions in `(-π, π]`.
//! * **Exponential–polynomial equations** `α + β x + γ e^{dx} = 0` and
//!   `α + β x e^{dx} = 0` by the Lambert W function (both real branches),
//!   `x` with `ln x` through `x = e^y`.
//! * **Radical equations**: powers `x^(p/q)` through `x = u^q`, and a
//!   single square root by isolating it and squaring.
//! * **Systems** that are not polynomial: elimination of an unknown that
//!   one equation determines, and recursion on the rest.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::difference;
use super::product;
use super::rational;
use super::reciprocal;
use super::solve_for;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::univariate;

fn value(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.eval(node, &Env::numeric(0.0)).filter(|v| v.is_finite())
}

fn power(
    graph: &mut Graph,
    base: NodeId,
    num: i64,
    den: i64,
) -> Option<NodeId> {
    let e = graph.num(Number::fraction(num, den)?);
    Some(graph.node(core::POW, &[base, e]))
}

/// The real cube root of a constant term, as `sign · |a|^(1/3)`.
fn real_cbrt(
    graph: &mut Graph,
    a: NodeId,
) -> Option<NodeId> {
    let v = value(graph, a)?;
    if v >= 0.0 {
        return power(graph, a, 1, 3);
    }
    let minus_one = graph.int(-1);
    let negated = graph.node(core::MUL, &[minus_one, a]);
    let root = power(graph, negated, 1, 3)?;
    Some(graph.node(core::MUL, &[minus_one, root]))
}

/// Real roots of an irreducible polynomial over `Q` of degree ≥ 3 (as
/// `(approximation, term)` pairs).
pub(super) fn higher_degree(
    graph: &mut Graph,
    coefficients: &[BigRational],
) -> Option<Vec<(f64, NodeId)>> {
    match coefficients.len() {
        | 4 => cubic(graph, coefficients),
        | 5 => quartic(graph, coefficients).or_else(|| Some(numeric_roots(graph, coefficients))),
        | _ => Some(numeric_roots(graph, coefficients)),
    }
}

/// Real roots by exact isolation and bisection to double precision.
fn numeric_roots(
    graph: &mut Graph,
    coefficients: &[BigRational],
) -> Vec<(f64, NodeId)> {
    let width = BigRational::new(BigInt::one(), BigInt::from(10_u64.pow(16)));
    let intervals = crate::rules::poly::algebra::isolate(coefficients);
    let mut out = Vec::with_capacity(intervals.len());
    for interval in intervals {
        let (a, b) = crate::rules::poly::algebra::refine_interval(coefficients, interval, &width);
        let mid = Number::rat((a + b) / BigRational::from_integer(BigInt::from(2))).to_f64();
        out.push((mid, graph.float(mid)));
    }
    out
}

fn cubic(
    graph: &mut Graph,
    c: &[BigRational],
) -> Option<Vec<(f64, NodeId)>> {
    let (d, cc, b, a) = (&c[0], &c[1], &c[2], &c[3]);
    let three = BigRational::from_integer(BigInt::from(3));
    // x = t - b/(3a): t³ + p t + q = 0.
    let shift = -(b / (&three * a));
    let p = (&three * a * cc - b * b) / (&three * a * a);
    let twenty_seven = BigRational::from_integer(BigInt::from(27));
    let nine = BigRational::from_integer(BigInt::from(9));
    let two = BigRational::from_integer(BigInt::from(2));
    let q = (&two * b * b * b - &nine * a * b * cc + &twenty_seven * a * a * d) / (&twenty_seven * a * a * a);
    let four = BigRational::from_integer(BigInt::from(4));
    let delta = -(&four * &p * &p * &p + &twenty_seven * &q * &q);
    let shift_node = rational(graph, shift.clone());
    let shift_f = Number::rat(shift).to_f64();
    let mut out = Vec::new();
    if delta.is_positive() {
        // t_k = 2 √(-p/3) cos(acos((3q/2p) √(-3/p))/3 - 2πk/3)
        let (acos, cos, pi) = (graph.ops().lookup("acos")?, graph.ops().lookup("cos")?, graph.ops().lookup("pi")?);
        let minus_p_third = rational(graph, -(&p / &three));
        let m = power(graph, minus_p_third, 1, 2)?;
        let two_node = graph.int(2);
        let amplitude = product(graph, &[two_node, m]);
        let minus_three_over_p = rational(graph, -(&three / &p));
        let root = power(graph, minus_three_over_p, 1, 2)?;
        let factor = rational(graph, &three * &q / (&two * &p));
        let argument = product(graph, &[factor, root]);
        let angle = graph.node(acos, &[argument]);
        let third = graph.num(Number::fraction(1, 3)?);
        let base_angle = product(graph, &[third, angle]);
        let pi = graph.node(pi, &[]);
        for k in 0..3_i64 {
            let offset = graph.num(Number::fraction(-2 * k, 3)?);
            let shift_angle = product(graph, &[offset, pi]);
            let total = graph.node(core::ADD, &[base_angle, shift_angle]);
            let wave = graph.node(cos, &[total]);
            let t = product(graph, &[amplitude, wave]);
            let x = graph.node(core::ADD, &[t, shift_node]);
            out.push((value(graph, x).unwrap_or(shift_f), x));
        }
    } else {
        // One real root: Cardano.
        let half_q = rational(graph, -(&q / &two));
        let radicand = rational(graph, &q * &q / &four + &p * &p * &p / &twenty_seven);
        let s = power(graph, radicand, 1, 2)?;
        let a1 = graph.node(core::ADD, &[half_q, s]);
        let minus_one = graph.int(-1);
        let neg_s = graph.node(core::MUL, &[minus_one, s]);
        let a2 = graph.node(core::ADD, &[half_q, neg_s]);
        let (r1, r2) = (real_cbrt(graph, a1)?, real_cbrt(graph, a2)?);
        let x = graph.node(core::ADD, &[r1, r2, shift_node]);
        out.push((value(graph, x)?, x));
    }
    Some(out)
}

fn quartic(
    graph: &mut Graph,
    c: &[BigRational],
) -> Option<Vec<(f64, NodeId)>> {
    let (e, d, cc, b, a) = (&c[0], &c[1], &c[2], &c[3], &c[4]);
    let mut out = Vec::new();
    if b.is_zero() && d.is_zero() {
        // a u² + c u + e with u = x².
        let disc = cc * cc - BigRational::from_integer(BigInt::from(4)) * a * e;
        if disc.is_negative() {
            return Some(out);
        }
        let disc_node = rational(graph, disc);
        let root = power(graph, disc_node, 1, 2)?;
        let two_a = rational(graph, BigRational::from_integer(BigInt::from(2)) * a);
        let inverse = reciprocal(graph, two_a);
        let minus_c = rational(graph, -cc.clone());
        for sign in [-1, 1] {
            let s = graph.int(sign);
            let signed = product(graph, &[s, root]);
            let top = graph.node(core::ADD, &[minus_c, signed]);
            let u = product(graph, &[top, inverse]);
            if value(graph, u).is_some_and(|v| v >= 0.0) {
                let x = power(graph, u, 1, 2)?;
                let minus_one = graph.int(-1);
                let neg = graph.node(core::MUL, &[minus_one, x]);
                let v = value(graph, x)?;
                out.push((-v, neg));
                out.push((v, x));
            }
        }
        return Some(out);
    }
    // Ferrari with a rational root of the resolvent: depressed
    // y⁴ + p y² + q y + r, x = y - b/(4a).
    let four = BigRational::from_integer(BigInt::from(4));
    let (b1, c1, d1, e1) = (b / a, cc / a, d / a, e / a);
    let three = BigRational::from_integer(BigInt::from(3));
    let eight = BigRational::from_integer(BigInt::from(8));
    let p = &c1 - &three * &b1 * &b1 / &eight;
    let q = &d1 - &b1 * &c1 / BigRational::from_integer(BigInt::from(2)) + &b1 * &b1 * &b1 / &eight;
    let r = &e1 - &b1 * &d1 / &four + &b1 * &b1 * &c1 / BigRational::from_integer(BigInt::from(16))
        - &three * &b1 * &b1 * &b1 * &b1 / BigRational::from_integer(BigInt::from(256));
    // 8m³ + 8p m² + (2p² - 8r) m - q² = 0
    let two = BigRational::from_integer(BigInt::from(2));
    let resolvent = vec![-(&q * &q), &two * &p * &p - &eight * &r, &eight * &p, eight];
    let m = univariate::rational_roots(&resolvent).into_iter().find(Signed::is_positive)?;
    let shift = -(&b1 / &four);
    let two_m = rational(graph, &two * &m);
    let s = power(graph, two_m, 1, 2)?;
    // y² ∓ s y + (p/2 + m ± q/(2s)) = 0
    let base = rational(graph, &p / &two + &m);
    let q_half = rational(graph, &q / &two);
    let inverse_s = reciprocal(graph, s);
    let q_term = product(graph, &[q_half, inverse_s]);
    let shift_node = rational(graph, shift);
    for sign in [-1_i64, 1] {
        // y² + sign·s·y + (base - sign·q_term) = 0
        let sgn = graph.int(sign);
        let linear = product(graph, &[sgn, s]);
        let minus_sgn = graph.int(-sign);
        let signed_q = product(graph, &[minus_sgn, q_term]);
        let constant = graph.node(core::ADD, &[base, signed_q]);
        // y = (-linear ± √(linear² - 4 constant))/2
        let two_node = graph.int(2);
        let sq = graph.node(core::POW, &[linear, two_node]);
        let minus_four = graph.int(-4);
        let fc = product(graph, &[minus_four, constant]);
        let disc = graph.node(core::ADD, &[sq, fc]);
        if value(graph, disc).is_none_or(|v| v < -1e-12) {
            continue;
        }
        let root = power(graph, disc, 1, 2)?;
        let half = graph.num(Number::fraction(1, 2)?);
        let minus_one = graph.int(-1);
        let neg_linear = product(graph, &[minus_one, linear]);
        for pm in [-1, 1] {
            let pm_node = graph.int(pm);
            let signed_root = product(graph, &[pm_node, root]);
            let top = graph.node(core::ADD, &[neg_linear, signed_root]);
            let y = product(graph, &[half, top]);
            let x = graph.node(core::ADD, &[y, shift_node]);
            out.push((value(graph, x)?, x));
        }
    }
    Some(out)
}

/// Every node of `term` with operator `op`.
pub(super) fn nodes_with(
    graph: &Graph,
    term: NodeId,
    pred: impl Fn(&Graph, NodeId) -> bool,
) -> Vec<NodeId> {
    let mut out = Vec::new();
    let mut stack = vec![term];
    while let Some(n) = stack.pop() {
        if pred(graph, n) && !out.contains(&n) {
            out.push(n);
        }
        stack.extend_from_slice(graph.children(n));
    }
    out
}

/// Entry point for equations with several generators of the unknown.
pub(super) fn multi_generator(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    let depends = |graph: &Graph, n: NodeId| graph.depends_on(graph.find(n), symbol);
    for method in [
        absolute_values,
        trigonometric,
        hyperbolic,
        lambert,
        self_power,
        logarithm_sum,
        logarithmic,
        radicals,
        radicals_general,
        inverse_trig,
        log_both_sides,
        substitution,
    ] {
        if let Some(found) = method(graph, term, x, depth) {
            return Some(found);
        }
    }
    let _ = depends;
    None
}

/// Hyperbolic functions written through `exp`.
fn hyperbolic(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    let rewritten = super::normalize::hyperbolic_to_exponentials(graph, term, symbol)?;
    solve_for(graph, rewritten, x, depth + 1)
}

/// `|g|` replaced by `g` and by `-g`.
fn absolute_values(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let abs = graph.ops().lookup("abs")?;
    let symbol = graph.symbol_of(x)?;
    let found = nodes_with(graph, term, |g, n| g.op(n) == abs && g.depends_on(g.find(n), symbol));
    let &first = found.first()?;
    let inner = *graph.children(first).first()?;
    let mut out = Vec::new();
    for sign in [1, -1] {
        let s = graph.int(sign);
        let replacement = graph.node(core::MUL, &[s, inner]);
        let case = graph.replace_subterm(term, first, replacement);
        out.extend(solve_for(graph, case, x, depth + 1).unwrap_or_default());
    }
    Some(out)
}

/// Linear trigonometric equations (see the trig module), else the
/// Weierstrass substitution for equations in sin, cos, tan of one
/// argument of any shape.
fn trigonometric(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    if let Some(found) = super::trig::solve(graph, term, x, super::general_mode()) {
        return Some(found);
    }
    weierstrass(graph, term, x, depth)
}

/// Weierstrass substitution for equations in sin, cos, tan of one
/// argument (expanded first).
fn weierstrass(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (sin, cos, tan, atan, pi) = (
        graph.ops().lookup("sin")?,
        graph.ops().lookup("cos")?,
        graph.ops().lookup("tan")?,
        graph.ops().lookup("atan")?,
        graph.ops().lookup("pi")?,
    );
    let symbol = graph.symbol_of(x)?;
    let expanded = crate::rules::poly::algebra::expand_trig_term(graph, term).unwrap_or(term);
    let trig = nodes_with(graph, expanded, |g, n| {
        (g.op(n) == sin || g.op(n) == cos || g.op(n) == tan) && g.depends_on(g.find(n), symbol)
    });
    let mut arguments: Vec<NodeId> = trig.iter().filter_map(|&n| graph.children(n).first().copied()).collect();
    arguments.dedup();
    let &[a] = arguments.as_slice() else {
        return None;
    };
    let t_symbol = graph.interner_mut().fresh_symbol("t");
    let t = graph.symbol_node(t_symbol);
    let one = graph.int(1);
    let two = graph.int(2);
    let minus_one = graph.int(-1);
    let t2 = graph.node(core::POW, &[t, two]);
    let one_plus = graph.node(core::ADD, &[one, t2]);
    let neg_t2 = graph.node(core::MUL, &[minus_one, t2]);
    let one_minus = graph.node(core::ADD, &[one, neg_t2]);
    let inv_plus = reciprocal(graph, one_plus);
    let inv_minus = reciprocal(graph, one_minus);
    let s_val = product(graph, &[two, t, inv_plus]);
    let c_val = product(graph, &[one_minus, inv_plus]);
    let t_val = product(graph, &[two, t, inv_minus]);
    let mut substituted = expanded;
    for (op, replacement) in [(sin, s_val), (cos, c_val), (tan, t_val)] {
        let node = graph.node(op, &[a]);
        substituted = graph.replace_subterm(substituted, node, replacement);
    }
    if graph.depends_on(graph.find(substituted), symbol) {
        return None;
    }
    let mut out = Vec::new();
    for t_value in solve_for(graph, substituted, t, depth + 1)? {
        let angle = graph.node(atan, &[t_value]);
        let target = product(graph, &[two, angle]);
        let equation = difference(graph, a, target);
        out.extend(solve_for(graph, equation, x, depth + 1).unwrap_or_default());
    }
    // a = π, where t = tan(a/2) is infinite.
    let pi = graph.node(pi, &[]);
    let at_pi = difference(graph, a, pi);
    out.extend(solve_for(graph, at_pi, x, depth + 1).unwrap_or_default());
    Some(out)
}

/// `α + β x + γ e^{u}` and `α + β x e^{u}` with `u = d x + e`.
fn lambert(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    _depth: usize,
) -> Option<Vec<NodeId>> {
    let (exp, w0, wm1) = (graph.ops().lookup("exp")?, graph.ops().lookup("lambertw")?, graph.ops().lookup("lambertw_m1")?);
    let symbol = graph.symbol_of(x)?;
    let term = super::normalize::powers_to_exp(graph, term, symbol).unwrap_or(term);
    let exps = nodes_with(graph, term, |g, n| g.op(n) == exp && g.depends_on(g.find(n), symbol));
    let &[e_node] = exps.as_slice() else {
        return None;
    };
    let u = *graph.children(e_node).first()?;
    // u = d x + e
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let u_poly = crate::rules::poly::repr::from_term(graph, &mut gens, u, Limits::default())?;
    if u_poly.degree_in(gx) != 1 || u_poly.support().iter().any(|&g| g != gx && gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol))) {
        return None;
    }
    let coefficients = u_poly.coefficients_in(gx);
    let d = crate::rules::poly::repr::to_term(graph, &gens, coefficients.get(1)?);
    let e_shift = crate::rules::poly::repr::to_term(graph, &gens, coefficients.first()?);
    // The equation as a polynomial in x and E = exp(u).
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let ge = gens.index(graph, e_node);
    let fraction = crate::rules::poly::ratio(graph, &mut gens, term, Limits::default())?;
    let numer = fraction.numer;
    for g in numer.support() {
        if g != gx && g != ge && gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)) {
            return None;
        }
    }
    let coefficient = |graph: &mut Graph, i: u32, j: u32| -> NodeId {
        let mut acc = crate::rules::poly::repr::Poly::zero();
        for (mono, c) in numer.terms() {
            let ex = mono.iter().find(|&&(g, _)| g == gx).map_or(0, |&(_, e)| e);
            let ee = mono.iter().find(|&&(g, _)| g == ge).map_or(0, |&(_, e)| e);
            if ex == i && ee == j {
                let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(g, _)| g != gx && g != ge).collect();
                acc = acc.add(&crate::rules::poly::repr::Poly::monomial(rest, c.clone()));
            }
        }
        crate::rules::poly::repr::to_term(graph, &gens, &acc)
    };
    let pattern: Vec<(u32, u32)> = numer
        .terms()
        .map(|(mono, _)| {
            (
                mono.iter().find(|&&(g, _)| g == gx).map_or(0, |&(_, e)| e),
                mono.iter().find(|&&(g, _)| g == ge).map_or(0, |&(_, e)| e),
            )
        })
        .collect();
    let allowed = |set: &[(u32, u32)]| pattern.iter().all(|p| set.contains(p));
    let minus_one = graph.int(-1);
    let e_e = graph.node(exp, &[e_shift]);
    let candidates = if allowed(&[(0, 0), (1, 0), (0, 1)]) && pattern.contains(&(1, 0)) && pattern.contains(&(0, 1)) {
        // α + β x + γ e^e e^{dx}: x = -α/β - W(z)/d, z = d γ e^e / β · e^{-dα/β}
        let (alpha, beta, gamma) = (coefficient(graph, 0, 0), coefficient(graph, 1, 0), coefficient(graph, 0, 1));
        let inv_beta = reciprocal(graph, beta);
        let ratio = product(graph, &[alpha, inv_beta]);
        let neg_d_ratio = product(graph, &[minus_one, d, ratio]);
        let decay = graph.node(exp, &[neg_d_ratio]);
        let z = product(graph, &[d, gamma, e_e, inv_beta, decay]);
        let neg_ratio = product(graph, &[minus_one, ratio]);
        let inv_d = reciprocal(graph, d);
        branches(graph, z, w0, wm1)
            .into_iter()
            .map(|w| {
                let part = product(graph, &[minus_one, w, inv_d]);
                graph.node(core::ADD, &[neg_ratio, part])
            })
            .collect::<Vec<_>>()
    } else if allowed(&[(0, 0), (1, 1)]) && pattern.contains(&(1, 1)) {
        // α + β x e^e e^{dx} = 0: x = W(z)/d, z = -d α e^{-e} / β
        let (alpha, beta) = (coefficient(graph, 0, 0), coefficient(graph, 1, 1));
        let inv_beta = reciprocal(graph, beta);
        let neg_e = product(graph, &[minus_one, e_shift]);
        let e_inv = graph.node(exp, &[neg_e]);
        let z = product(graph, &[minus_one, d, alpha, e_inv, inv_beta]);
        let inv_d = reciprocal(graph, d);
        branches(graph, z, w0, wm1).into_iter().map(|w| product(graph, &[w, inv_d])).collect()
    } else {
        power_exponential(graph, &pattern, &coefficient, d, e_shift, w0, wm1, exp)?
    };
    Some(candidates)
}

/// `α + β x^m e^{dx+e} = 0` and `α x^m + γ e^{dx+e} = 0` for an integer
/// `m >= 2`: taking `m`-th roots leaves `x e^{s x} = c`, so `x = W(s c)/s`.
#[allow(clippy::too_many_arguments)]
fn power_exponential(
    graph: &mut Graph,
    pattern: &[(u32, u32)],
    coefficient: &dyn Fn(&mut Graph, u32, u32) -> NodeId,
    d: NodeId,
    e_shift: NodeId,
    w0: crate::graph::OpId,
    wm1: crate::graph::OpId,
    exp: crate::graph::OpId,
) -> Option<Vec<NodeId>> {
    let m = pattern.iter().map(|p| p.0).max()?;
    if m < 2 {
        return None;
    }
    let only = |set: &[(u32, u32)]| pattern.iter().all(|p| set.contains(p));
    let a_form = only(&[(0, 0), (m, 1)]) && pattern.contains(&(m, 1));
    let b_form = only(&[(m, 0), (0, 1)]) && pattern.contains(&(m, 0)) && pattern.contains(&(0, 1));
    if !a_form && !b_form {
        return None;
    }
    let m_i = i64::from(m);
    let minus_one = graph.int(-1);
    let k = if a_form {
        let (alpha, beta) = (coefficient(graph, 0, 0), coefficient(graph, m, 1));
        let inv = reciprocal(graph, beta);
        product(graph, &[minus_one, alpha, inv])
    } else {
        let (alpha, gamma) = (coefficient(graph, m, 0), coefficient(graph, 0, 1));
        let inv = reciprocal(graph, alpha);
        product(graph, &[minus_one, gamma, inv])
    };
    let inv_m = graph.num(Number::fraction(1, m_i)?);
    let d_over_m = product(graph, &[d, inv_m]);
    let e_over_m = product(graph, &[e_shift, inv_m]);
    let (s, shift) = if a_form {
        let neg = product(graph, &[minus_one, e_over_m]);
        (d_over_m, graph.node(exp, &[neg]))
    } else {
        (product(graph, &[minus_one, d_over_m]), graph.node(exp, &[e_over_m]))
    };
    let root = if m % 2 == 1 {
        vec![super::symbolic::real_root(graph, k, m_i)?]
    } else {
        let magnitude = graph.node(core::POW, &[k, inv_m]);
        vec![magnitude, product(graph, &[minus_one, magnitude])]
    };
    let mut out = Vec::new();
    let inverse_s = reciprocal(graph, s);
    for r in root {
        let c = product(graph, &[r, shift]);
        let z = product(graph, &[s, c]);
        for w in branches(graph, z, w0, wm1) {
            out.push(product(graph, &[w, inverse_s]));
        }
    }
    Some(out)
}

/// `W₀(z)`, and `W₋₁(z)` when `z` is known to lie in `(-1/e, 0)`.
fn branches(
    graph: &mut Graph,
    z: NodeId,
    w0: crate::graph::OpId,
    wm1: crate::graph::OpId,
) -> Vec<NodeId> {
    let mut out = vec![graph.node(w0, &[z])];
    // Symbolic arguments: the second branch is real only on part of the
    // parameter space, where the final check keeps it.
    if value(graph, z).is_none_or(|v| v < 0.0 && v > -(-1.0_f64).exp()) {
        out.push(graph.node(wm1, &[z]));
    }
    out
}

/// `x^x = c`: `x ln x = ln c`, so `ln x = W(ln c)` and `x = e^(W(ln c))`
/// (both real branches when `-1/e < ln c < 0`).
fn self_power(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (exp, ln, w0, wm1) =
        (graph.ops().lookup("exp")?, graph.ops().lookup("ln")?, graph.ops().lookup("lambertw")?, graph.ops().lookup("lambertw_m1")?);
    let x_to_x = graph.node(core::POW, &[x, x]);
    if !super::super::ode::occurs_in(graph, term, x_to_x) {
        return None;
    }
    let w_symbol = graph.interner_mut().fresh_symbol("w");
    let w = graph.symbol_node(w_symbol);
    let replaced = graph.replace_subterm(term, x_to_x, w);
    if super::super::ode::occurs_in(graph, replaced, x) {
        return None;
    }
    let mut out = Vec::new();
    for c in solve_for(graph, replaced, w, depth + 1)? {
        if value(graph, c).is_some_and(|v| v <= 0.0) {
            continue;
        }
        let log_c = graph.node(ln, &[c]);
        for branch in branches(graph, log_c, w0, wm1) {
            out.push(graph.node(exp, &[branch]));
        }
    }
    Some(out)
}

/// `Σ k_i ln(g_i(x)) + C = 0` with integers `k_i` and `C` free of `x`:
/// `Π g_i^(k_i) = e^(-C)`, solved as a rational equation; candidates where
/// some `g_i` is not positive are discarded (real logarithms).
fn logarithm_sum(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (exp, ln) = (graph.ops().lookup("exp")?, graph.ops().lookup("ln")?);
    let symbol = graph.symbol_of(x)?;
    let summands = if graph.op(term) == core::ADD { graph.children(term).to_vec() } else { vec![term] };
    let mut logs: Vec<(NodeId, i64)> = Vec::new();
    let mut constants = Vec::new();
    for t in summands {
        if !graph.depends_on(graph.find(t), symbol) {
            constants.push(t);
            continue;
        }
        let (k, inner) = match (graph.op(t), graph.children(t)) {
            | (op, &[g]) if op == ln => (1, g),
            | (op, &[c, l]) if op == core::MUL && graph.op(l) == ln => (graph.number_of(c)?.to_i64()?, *graph.children(l).first()?),
            | (op, &[l, c]) if op == core::MUL && graph.op(l) == ln => (graph.number_of(c)?.to_i64()?, *graph.children(l).first()?),
            | _ => return None,
        };
        logs.push((inner, k));
    }
    if logs.len() < 2 && constants.is_empty() {
        return None;
    }
    let (mut positive, mut negative) = (Vec::new(), Vec::new());
    for &(g, k) in &logs {
        let e = graph.int(k.abs());
        let p = graph.node(core::POW, &[g, e]);
        if k > 0 { positive.push(p) } else { negative.push(p) }
    }
    let constant = if constants.is_empty() { graph.int(0) } else { graph.node(core::ADD, &constants) };
    let minus_one = graph.int(-1);
    let minus_c = graph.node(core::MUL, &[minus_one, constant]);
    let rhs_scale = graph.node(exp, &[minus_c]);
    let one = graph.int(1);
    let left = if positive.is_empty() { one } else { product(graph, &positive) };
    let mut right_factors = negative;
    right_factors.push(rhs_scale);
    let right = product(graph, &right_factors);
    let equation = difference(graph, left, right);
    let mut out = Vec::new();
    for candidate in solve_for(graph, equation, x, depth + 1)? {
        let all_positive = logs.iter().all(|&(g, _)| {
            let at = graph.substitute(g, x, candidate);
            graph.eval(at, &Env::numeric(0.0)).is_none_or(|v| v > 0.0)
        });
        if all_positive {
            out.push(candidate);
        }
    }
    Some(out)
}

/// `x` together with `ln(x)`: `x = e^y`.
fn logarithmic(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (ln, exp) = (graph.ops().lookup("ln")?, graph.ops().lookup("exp")?);
    let log_x = graph.node(ln, &[x]);
    if !super::super::ode::occurs_in(graph, term, log_x) {
        return None;
    }
    let y_symbol = graph.interner_mut().fresh_symbol("y");
    let y = graph.symbol_node(y_symbol);
    let e_y = graph.node(exp, &[y]);
    let replaced = graph.replace_subterm(term, log_x, y);
    let replaced = graph.replace_subterm(replaced, x, e_y);
    let mut out = Vec::new();
    for value in solve_for(graph, replaced, y, depth + 1)? {
        out.push(graph.node(exp, &[value]));
    }
    Some(out)
}

/// Fractional powers of `x` through `x = u^q`, or one square root
/// isolated and squared.
fn radicals(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    let roots = nodes_with(graph, term, |g, n| {
        g.op(n) == core::POW
            && g.children(n).get(1).and_then(|&e| g.number_of(e)).is_some_and(|e| !e.is_integer())
            && g.depends_on(g.find(n), symbol)
    });
    if roots.is_empty() {
        return None;
    }
    // Powers of x itself: x = u^q with q the lcm of the denominators.
    if roots.iter().all(|&r| graph.children(r).first() == Some(&x)) {
        let mut q = BigInt::one();
        for &r in &roots {
            let e = graph.number_of(*graph.children(r).get(1)?)?.to_rational()?;
            q = num_integer::Integer::lcm(&q, e.denom());
        }
        let q_i = i64::try_from(q).ok()?;
        let u_symbol = graph.interner_mut().fresh_symbol("u");
        let u = graph.symbol_node(u_symbol);
        let mut replaced = term;
        for &r in &roots {
            let e = graph.number_of(*graph.children(r).get(1)?)?.to_rational()?;
            let k = e * BigRational::from_integer(BigInt::from(q_i));
            let k = graph.num(Number::rat(k));
            let pow = graph.node(core::POW, &[u, k]);
            replaced = graph.replace_subterm(replaced, r, pow);
        }
        let q_node = graph.int(q_i);
        let x_as = graph.node(core::POW, &[u, q_node]);
        replaced = graph.replace_subterm(replaced, x, x_as);
        let mut out = Vec::new();
        for value in solve_for(graph, replaced, u, depth + 1)? {
            if value_of_nonnegative(graph, value) {
                out.push(graph.node(core::POW, &[value, q_node]));
            }
        }
        return Some(out);
    }
    // Square roots isolated one at a time: with every root a generator and
    // s_i² reduced to f_i, the numerator is A + B s for the first root s,
    // and A² - B² f is free of s; the recursive call removes the rest.
    // Candidates introduced by squaring are dropped by `solve_for`.
    let half = BigRational::new(BigInt::one(), BigInt::from(2));
    let mut gens = Gens::default();
    gens.index(graph, x);
    let mut square_roots: Vec<(u32, NodeId)> = Vec::new();
    for &r in &roots {
        if graph.number_of(*graph.children(r).get(1)?)?.to_rational()? != half {
            return None;
        }
        let g = gens.index(graph, r);
        square_roots.push((g, *graph.children(r).first()?));
    }
    let fraction = crate::rules::poly::ratio(graph, &mut gens, term, Limits::default())?;
    let mut radicands = Vec::with_capacity(square_roots.len());
    for &(g, f) in &square_roots {
        let f = crate::rules::poly::best(graph, f)?;
        let p = crate::rules::poly::repr::from_term(graph, &mut gens, f, Limits::default())?;
        if square_roots.iter().any(|&(h, _)| p.degree_in(h) > 0) {
            return None;
        }
        radicands.push((g, p));
    }
    let reduce = |p: &Poly| -> Option<Poly> {
        let mut p = p.clone();
        for (g, f) in &radicands {
            let mut out = Poly::zero();
            for (k, c) in p.coefficients_in(*g).iter().enumerate() {
                let k = u32::try_from(k).ok()?;
                let mut term = c.mul(&f.pow(k / 2, 4096)?, 4096)?;
                if k % 2 == 1 {
                    term = term.mul(&Poly::generator(*g), 4096)?;
                }
                out = out.add(&term);
            }
            p = out;
        }
        Some(p)
    };
    let numer = reduce(&fraction.numer)?;
    let (g, f) = radicands.iter().find(|(g, _)| numer.degree_in(*g) == 1)?;
    let parts = numer.coefficients_in(*g);
    let [a, b] = parts.as_slice() else {
        return None;
    };
    let squared = a.mul(a, 4096)?.sub(&b.mul(b, 4096)?.mul(f, 4096)?);
    let squared = reduce(&squared)?;
    let squared = crate::rules::poly::repr::to_term(graph, &gens, &squared);
    solve_for(graph, squared, x, depth + 1)
}

fn value_of_nonnegative(
    graph: &Graph,
    node: NodeId,
) -> bool {
    value(graph, node).is_none_or(|v| v >= 0.0)
}

/// Radicals of any index: each radical `r = b^(p/q)` is a generator with the
/// relation `r^q = b^p`, and the resultant with respect to `r` removes it
/// from the numerator (one radical at a time). Square roots alone are left
/// to [`radicals`].
fn radicals_general(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    let roots = nodes_with(graph, term, |g, n| {
        g.op(n) == core::POW
            && g.children(n).get(1).and_then(|&e| g.number_of(e)).is_some_and(|e| !e.is_integer())
            && g.depends_on(g.find(n), symbol)
    });
    if roots.is_empty() || roots.len() > 3 {
        return None;
    }
    let limits = Limits::default();
    let mut gens = Gens::default();
    gens.index(graph, x);
    let mut radicals: Vec<(u32, i64, i64, NodeId)> = Vec::new();
    let mut max_index = 0;
    for &r in &roots {
        let &[base, e] = graph.children(r) else {
            return None;
        };
        let e = graph.number_of(e)?.to_rational()?;
        let (p, q) = (i64::try_from(e.numer().clone()).ok()?, i64::try_from(e.denom().clone()).ok()?);
        max_index = max_index.max(q);
        radicals.push((gens.index(graph, r), p, q, base));
    }
    if max_index < 3 {
        return None;
    }
    let fraction = crate::rules::poly::ratio(graph, &mut gens, term, limits)?;
    let mut numer = fraction.numer;
    for &(g, p, q, base) in &radicals {
        let base = crate::rules::poly::best(graph, base)?;
        let base_poly = crate::rules::poly::repr::from_term(graph, &mut gens, base, limits)?;
        if radicals.iter().any(|&(h, ..)| base_poly.degree_in(h) > 0) {
            return None;
        }
        let power = base_poly.pow(u32::try_from(p.unsigned_abs()).ok()?, limits.terms)?;
        let g_q = Poly::generator(g).pow(u32::try_from(q).ok()?, limits.terms)?;
        let relation = if p > 0 {
            g_q.sub(&power)
        } else {
            g_q.mul(&power, limits.terms)?.sub(&Poly::constant(Number::from(1)))
        };
        if numer.degree_in(g) == 0 {
            continue;
        }
        numer = super::elim::resultant(&numer, &relation, g)?;
    }
    if numer.is_zero() {
        return None;
    }
    let reduced = crate::rules::poly::repr::to_term(graph, &gens, &numer);
    solve_for(graph, reduced, x, depth + 1)
}

/// Equations `±asin/acos/atan(g1) ± asin/acos/atan(g2) + c = 0`: the
/// tangent of both sides gives an algebraic equation (the final check
/// removes what the tangent introduced).
fn inverse_trig(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (asin, acos, atan, tan) =
        (graph.ops().lookup("asin")?, graph.ops().lookup("acos")?, graph.ops().lookup("atan")?, graph.ops().lookup("tan")?);
    let symbol = graph.symbol_of(x)?;
    let summands = if graph.op(term) == core::ADD { graph.children(term).to_vec() } else { vec![term] };
    let mut inverse: Vec<(i64, NodeId)> = Vec::new();
    let mut constants = Vec::new();
    for t in summands {
        if !graph.depends_on(graph.find(t), symbol) {
            constants.push(t);
            continue;
        }
        let (sign, h) = match (graph.op(t), graph.children(t).to_vec().as_slice()) {
            | (op, _) if op == asin || op == acos || op == atan => (1, t),
            // `±1 * h` with the factors in either order.
            | (op, &[a, b]) if op == core::MUL => {
                let unit = |n: NodeId| graph.number_of(n).and_then(Number::to_i64).filter(|v| v.abs() == 1);
                match (unit(a), unit(b)) {
                    | (Some(sign), None) => (sign, b),
                    | (None, Some(sign)) => (sign, a),
                    | _ => return None,
                }
            },
            | _ => return None,
        };
        let hop = graph.op(h);
        if hop != asin && hop != acos && hop != atan {
            return None;
        }
        inverse.push((sign, h));
    }
    let [(s1, h1), (s2, h2)] = inverse.as_slice() else {
        return None;
    };
    let tangent = |graph: &mut Graph, h: NodeId| -> Option<NodeId> {
        let g = *graph.children(h).first()?;
        let (one, minus_one, two) = (graph.int(1), graph.int(-1), graph.int(2));
        let half = graph.num(Number::fraction(1, 2)?);
        let g_sq = graph.node(core::POW, &[g, two]);
        let neg = graph.node(core::MUL, &[minus_one, g_sq]);
        let complement = graph.node(core::ADD, &[one, neg]);
        Some(if graph.op(h) == asin {
            let root = graph.node(core::POW, &[complement, half]);
            let inv_root = reciprocal(graph, root);
            product(graph, &[g, inv_root])
        } else if graph.op(h) == acos {
            let root = graph.node(core::POW, &[complement, half]);
            let inv = reciprocal(graph, g);
            product(graph, &[root, inv])
        } else {
            g
        })
    };
    let (t1, t2) = (tangent(graph, *h1)?, tangent(graph, *h2)?);
    let s = s1 * s2;
    // s1 (h1 + s h2) + c = 0  =>  h1 + s h2 + s1 c = 0
    let constant = if constants.is_empty() {
        None
    } else {
        let sum = if constants.len() == 1 { constants[0] } else { graph.node(core::ADD, &constants) };
        let signed = graph.int(*s1);
        Some(product(graph, &[signed, sum]))
    };
    let tc = match constant {
        | Some(c) => {
            let t = graph.node(tan, &[c]);
            if value(graph, t).is_none_or(|v| v.abs() > 1e9) {
                return None;
            }
            Some(t)
        },
        | None => None,
    };
    let s_node = graph.int(s);
    // t1 - s t1 t2 tc + s t2 + tc
    let mut parts = vec![t1, product(graph, &[s_node, t2])];
    if let Some(tc) = tc {
        let minus_s = graph.int(-s);
        parts.push(product(graph, &[minus_s, t1, t2, tc]));
        parts.push(tc);
    }
    let equation = graph.node(core::ADD, &parts);
    solve_for(graph, equation, x, depth + 1)
}

/// `ln` of every factor and power taken apart (positive arguments
/// assumed): `ln(a^b c) = b ln a + ln c`.
fn expand_log(
    graph: &mut Graph,
    node: NodeId,
    ln: crate::graph::OpId,
    exp: crate::graph::OpId,
) -> NodeId {
    let children = graph.children(node).to_vec();
    if graph.op(node) == core::MUL {
        let parts: Vec<NodeId> = children.iter().map(|&c| expand_log(graph, c, ln, exp)).collect();
        return graph.node(core::ADD, &parts);
    }
    if graph.op(node) == core::POW
        && let &[base, power] = children.as_slice() {
            let log = expand_log(graph, base, ln, exp);
            return product(graph, &[power, log]);
        }
    if graph.op(node) == exp
        && let Some(&a) = children.first() {
            return a;
        }
    graph.node(ln, &[node])
}

/// `t1 + t2 = 0` where a term has an exponent that involves the unknown:
/// `ln|t1| = ln|t2|`, with the logarithms expanded.
fn log_both_sides(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let (ln, exp) = (graph.ops().lookup("ln")?, graph.ops().lookup("exp")?);
    let symbol = graph.symbol_of(x)?;
    if graph.op(term) != core::ADD {
        return None;
    }
    let &[t1, t2] = graph.children(term) else {
        return None;
    };
    let variable_exponent = nodes_with(graph, term, |g, n| {
        g.op(n) == core::POW && g.children(n).get(1).is_some_and(|&e| g.depends_on(g.find(e), symbol))
    });
    if variable_exponent.is_empty() {
        return None;
    }
    // Numeric sign of a term: the product of its numeric factors.
    let sign_of = |graph: &Graph, t: NodeId| -> f64 {
        match graph.number_of(t) {
            | Some(n) => n.to_f64(),
            | None if graph.op(t) == core::MUL => {
                graph.children(t).iter().filter_map(|&c| graph.number_of(c)).map(Number::to_f64).product()
            },
            | None => 1.0,
        }
    };
    let strip = |graph: &mut Graph, t: NodeId| -> NodeId {
        // Remove numeric factors from a product, keeping their magnitude.
        if graph.op(t) != core::MUL {
            return t;
        }
        let kids = graph.children(t).to_vec();
        let magnitude: Vec<NodeId> = kids
            .iter()
            .map(|&c| match graph.number_of(c) {
                | Some(n) if n.to_f64() < 0.0 => graph.num(n.neg()),
                | _ => c,
            })
            .collect();
        graph.node(core::MUL, &magnitude)
    };
    // t1 = -t2: the signs of t1 and -t2 must agree (numeric factors).
    let (s1, s2) = (sign_of(graph, t1), -sign_of(graph, t2));
    if s1 * s2 < 0.0 {
        return None;
    }
    let (m1, m2) = (strip(graph, t1), strip(graph, t2));
    let left = expand_log(graph, m1, ln, exp);
    let right = expand_log(graph, m2, ln, exp);
    let equation = difference(graph, left, right);
    solve_for(graph, equation, x, depth + 1)
}

fn size(
    graph: &Graph,
    node: NodeId,
) -> usize {
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if seen.insert(n) {
            stack.extend_from_slice(graph.children(n));
        }
    }
    seen.len()
}

/// A common subterm `f(x)` through which the unknown enters the equation:
/// with `u = f(x)` the equation is solved for `u` and each root inverted.
fn substitution(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    let mut candidates: Vec<NodeId> = Vec::new();
    for n in nodes_with(graph, term, |g, n| !g.children(n).is_empty()) {
        let op = graph.op(n);
        if op == core::ADD || op == core::MUL {
            continue;
        }
        for &c in graph.children(n) {
            if c != x && !graph.children(c).is_empty() && graph.depends_on(graph.find(c), symbol) && !candidates.contains(&c) {
                candidates.push(c);
            }
        }
    }
    candidates.sort_by_key(|&c| std::cmp::Reverse(size(graph, c)));
    for f in candidates {
        let u_symbol = graph.interner_mut().fresh_symbol("u");
        let u = graph.symbol_node(u_symbol);
        let mut replaced = term;
        // x^(km) = u^k when f = x^m.
        if let (core::POW, &[base, m]) = (graph.op(f), graph.children(f))
            && base == x
                && let Some(m) = graph.number_of(m).and_then(Number::to_i64).filter(|&m| m >= 2) {
                    let mut rewrite = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
                        if graph.op(node) != core::POW {
                            return None;
                        }
                        let &[b, e] = children else {
                            return None;
                        };
                        let k = graph.number_of(e).and_then(Number::to_i64)?;
                        (b == x && k % m == 0).then(|| {
                            let power = graph.int(k / m);
                            graph.node(core::POW, &[u, power])
                        })
                    };
                    replaced = super::normalize::map_term(graph, replaced, &mut rewrite);
                }
        replaced = graph.replace_subterm(replaced, f, u);
        if graph.depends_on(graph.find(replaced), symbol) {
            continue;
        }
        let Some(values) = solve_for(graph, replaced, u, depth + 1) else {
            continue;
        };
        let mut out = Vec::new();
        let mut complete = true;
        for value in values {
            let equation = difference(graph, f, value);
            match solve_for(graph, equation, x, depth + 1) {
                | Some(found) => out.extend(found),
                | None => {
                    complete = false;
                    break;
                },
            }
        }
        if complete {
            return Some(out);
        }
    }
    None
}
