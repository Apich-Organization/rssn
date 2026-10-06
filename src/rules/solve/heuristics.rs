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
fn nodes_with(
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
    for method in [absolute_values, trigonometric, lambert, logarithmic, radicals] {
        if let Some(found) = method(graph, term, x, depth) {
            return Some(found);
        }
    }
    let _ = depends;
    None
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

/// Weierstrass substitution for equations in sin, cos, tan of one
/// argument (expanded first).
fn trigonometric(
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
        return None;
    };
    Some(candidates)
}

/// `W₀(z)`, and `W₋₁(z)` when `z` is known to lie in `(-1/e, 0)`.
fn branches(
    graph: &mut Graph,
    z: NodeId,
    w0: crate::graph::OpId,
    wm1: crate::graph::OpId,
) -> Vec<NodeId> {
    let mut out = vec![graph.node(w0, &[z])];
    if value(graph, z).is_some_and(|v| v < 0.0 && v > -(-1.0_f64).exp()) {
        out.push(graph.node(wm1, &[z]));
    }
    out
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
    // A single square root: A + B √f = 0  ⇒  A² - B² f = 0.
    let &[root] = roots.as_slice() else {
        return None;
    };
    let exponent = graph.number_of(*graph.children(root).get(1)?)?.to_rational()?;
    if exponent != BigRational::new(BigInt::one(), BigInt::from(2)) {
        return None;
    }
    let f = *graph.children(root).first()?;
    let r_symbol = graph.interner_mut().fresh_symbol("r");
    let r = graph.symbol_node(r_symbol);
    let replaced = graph.replace_subterm(term, root, r);
    let mut gens = Gens::default();
    let gr = gens.index(graph, r);
    let fraction = crate::rules::poly::ratio(graph, &mut gens, replaced, Limits::default())?;
    let parts = fraction.numer.coefficients_in(gr);
    let [a, b] = parts.as_slice() else {
        return None;
    };
    let (a, b) = (crate::rules::poly::repr::to_term(graph, &gens, a), crate::rules::poly::repr::to_term(graph, &gens, b));
    let two = graph.int(2);
    let a2 = graph.node(core::POW, &[a, two]);
    let b2 = graph.node(core::POW, &[b, two]);
    let b2f = product(graph, &[b2, f]);
    let squared = difference(graph, a2, b2f);
    solve_for(graph, squared, x, depth + 1)
}

fn value_of_nonnegative(
    graph: &Graph,
    node: NodeId,
) -> bool {
    value(graph, node).is_none_or(|v| v >= 0.0)
}

/// Systems that are not polynomial: eliminate an unknown that one
/// equation determines, recurse, and back-substitute.
pub(super) fn eliminate(
    graph: &mut Graph,
    equations: &[NodeId],
    unknowns: &[NodeId],
    depth: usize,
) -> Option<Vec<Vec<NodeId>>> {
    if depth > 4 {
        return None;
    }
    if unknowns.is_empty() {
        return Some(vec![Vec::new()]);
    }
    let exprs: Vec<NodeId> = equations.iter().map(|&e| super::as_expression(graph, e)).collect();
    for (i, &expr) in exprs.iter().enumerate() {
        for (j, &u) in unknowns.iter().enumerate() {
            let Some(values) = solve_for(graph, expr, u, 0) else {
                continue;
            };
            if values.is_empty() {
                continue;
            }
            let rest_eqs: Vec<NodeId> = exprs.iter().enumerate().filter(|&(k, _)| k != i).map(|(_, &e)| e).collect();
            let rest_unknowns: Vec<NodeId> = unknowns.iter().enumerate().filter(|&(k, _)| k != j).map(|(_, &v)| v).collect();
            let mut out = Vec::new();
            for value in values {
                let substituted: Vec<NodeId> = rest_eqs.iter().map(|&e| graph.substitute(e, u, value)).collect();
                let tails = if rest_unknowns.is_empty() {
                    // Every remaining equation must hold.
                    let ok = substituted.iter().all(|&e| value_of_zero(graph, e));
                    if ok { vec![Vec::new()] } else { Vec::new() }
                } else {
                    // Underdetermined or unsolvable: give up rather than
                    // report a partial solution set.
                    eliminate(graph, &substituted, &rest_unknowns, depth + 1)?
                };
                for tail in tails {
                    // Insert u's value (with the others substituted) at j.
                    let mut full = tail.clone();
                    let mut v = value;
                    for (k, &other) in rest_unknowns.iter().enumerate() {
                        v = graph.substitute(v, other, tail[k]);
                    }
                    full.insert(j, v);
                    out.push(full);
                }
            }
            return (!out.is_empty()).then_some(out);
        }
    }
    None
}

fn value_of_zero(
    graph: &Graph,
    e: NodeId,
) -> bool {
    value(graph, e).is_none_or(|v| v.abs() < 1e-9)
}
