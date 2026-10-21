//! Linear equations of order three and more with polynomial coefficients.
//!
//! * **Generalised hypergeometric matching.** Written in `θ = x d/dx`, an
//!   operator with only two monomials, `x^{j0} P_0(θ) + x^{j1} P_1(θ)`
//!   (`d = j1 - j0 > 0`), has the series solutions
//!   `x^ρ pFq(a; b; c x^d)` for the rational roots `ρ` of `P_0`, with
//!   `a_i = (ρ - σ_i)/d` (`σ_i` the roots of `P_1`), `b_j = 1 + (ρ - ρ_j)/d`
//!   (the other roots of `P_0`) and `c = -lead(P_1) d^{deg P_1} / (lead(P_0)
//!   d^{deg P_0})`. All roots must be rational and none of the `b_j` may be
//!   a non-positive integer (that would call for logarithms).
//! * **Hyperexponential solutions and reduction of order.** A solution
//!   `y = P(x) Π f^ρ exp(E)` (`P` a polynomial, `f` the irreducible factors
//!   of the leading coefficient, `E` rational with poles at their roots and
//!   a polynomial part) is sought from the local formal solutions
//!   ([`super::formal`]) at every singular point and at infinity: they give
//!   the candidates for each `ρ` and the exponential part, the behaviour at
//!   infinity fixes the degree of `P`, and `P` solves a linear system.
//!   Then `L = M (D - y'/y)` is an exact right division, the quotient `M`
//!   of order one less is solved by the full solver (so Kovacic, the
//!   hypergeometric families and this method again apply), and each
//!   solution `z` of `M` gives `y1 ∫ z / y1`.
//!
//! Not done: the eigenring (Beke/van Hoeij) factorisation into factors of
//! order greater than one, left factors, and solutions whose exponents are
//! irrational.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::Problem;
use super::add;
use super::coefficients_in;
use super::explicit;
use super::formal::formal_solutions;
use super::formal::Point as FormalPoint;
use super::integrate;
use super::kovacic::polynomial_node;
use super::kovacic::Rat;
use super::mul;
use super::occurs;
use super::reduce::chain;
use super::reduce::fresh_function;
use super::solve_equation;
use crate::graph::op::core;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::rational_function_in;
use crate::rules::poly::univariate;

type Q = BigRational;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn trim(mut p: Vec<Q>) -> Vec<Q> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

fn rat(
    num: &[Q],
    den: &[Q],
) -> Option<Rat> {
    Rat::new(num, den)
}

fn zero_rat() -> Option<Rat> {
    rat(&[], &[Q::one()])
}

// ---------------------------------------------------------------------
// Extraction.

/// The polynomial coefficients `a_0, …, a_n` of a homogeneous linear
/// equation of order `n >= 3` (denominators cleared).
fn operator(
    cx: &mut Cx<'_>,
    problem: &Problem,
) -> Option<Vec<Vec<Q>>> {
    let order = problem.order();
    if order < 3 {
        return None;
    }
    let mut rest = problem.expr;
    let mut coefficients = vec![NodeId::NONE; order + 1];
    for k in (0..=order).rev() {
        let parts = coefficients_in(cx.graph, rest, *problem.stand.get(k)?)?;
        match *parts.as_slice() {
            | [constant, linear] => {
                *coefficients.get_mut(k)? = linear;
                rest = constant;
            },
            | [constant] => {
                *coefficients.get_mut(k)? = cx.graph.int(0);
                rest = constant;
            },
            | _ => return None,
        }
    }
    if !cx.is_zero(rest) || problem.stand.iter().any(|&s| occurs(cx.graph, rest, s) || coefficients.iter().any(|&c| occurs(cx.graph, c, s))) {
        return None;
    }
    let mut fractions = Vec::new();
    for &c in &coefficients {
        fractions.push(rational_function_in(cx.graph, c, problem.x)?);
    }
    let mut out = Vec::new();
    for (k, (numerator, _)) in fractions.iter().enumerate() {
        let mut value = numerator.clone();
        for (j, (_, denominator)) in fractions.iter().enumerate() {
            if j != k {
                value = univariate::mul(&value, denominator);
            }
        }
        out.push(trim(value));
    }
    (out.last().is_some_and(|l| !l.is_empty())).then_some(out)
}

// ---------------------------------------------------------------------
// Hypergeometric matching.

fn falling_theta(i: usize) -> Vec<Q> {
    let mut poly = vec![Q::one()];
    for m in 0..i {
        poly = univariate::mul(&poly, &[-q(i64::try_from(m).unwrap_or(0)), Q::one()]);
    }
    poly
}

/// `x^j P_j(θ)` pieces of the operator, `j -> P_j`.
fn euler_pieces(op: &[Vec<Q>]) -> std::collections::BTreeMap<i64, Vec<Q>> {
    let mut out: std::collections::BTreeMap<i64, Vec<Q>> = std::collections::BTreeMap::new();
    for (i, a) in op.iter().enumerate() {
        let ff = falling_theta(i);
        for (k, c) in a.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let j = i64::try_from(k).unwrap_or(0) - i64::try_from(i).unwrap_or(0);
            let scaled: Vec<Q> = ff.iter().map(|v| v * c).collect();
            let entry = out.entry(j).or_default();
            *entry = trim(univariate::add(entry, &scaled));
        }
    }
    out.retain(|_, p| !p.is_empty());
    out
}

fn rational_roots_of(p: &[Q]) -> Option<Vec<Q>> {
    let (_, factors) = univariate::factor(p);
    let mut roots = Vec::new();
    for (f, multiplicity) in factors {
        let [c0, c1] = f.as_slice() else { return None };
        let root = Q::new(-c0.clone(), c1.clone());
        for _ in 0..multiplicity {
            roots.push(root.clone());
        }
    }
    Some(roots)
}

/// One generalised hypergeometric series solution `x^ρ pFq(a; b; c x^d)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pfq {
    /// The exponent `ρ`.
    pub rho: Q,
    /// The numerator parameters.
    pub a: Vec<Q>,
    /// The denominator parameters.
    pub b: Vec<Q>,
    /// The argument `c x^d`: coefficient.
    pub c: Q,
    /// The argument `c x^d`: exponent.
    pub d: i64,
}

/// The `n` series solutions of an operator of order `n` that is
/// hypergeometric in `θ`, if it is.
#[must_use]
pub fn hypergeometric_solutions(op: &[Vec<Q>]) -> Option<Vec<Pfq>> {
    let n = op.len().checked_sub(1)?;
    let pieces = euler_pieces(op);
    let mut items = pieces.iter();
    let ((&j0, p0), (&j1, p1)) = (items.next()?, items.next()?);
    if items.next().is_some() || p0.len() != n + 1 || p1.len() > n + 1 {
        return None;
    }
    let d = j1 - j0;
    let dq = q(d);
    let roots0 = rational_roots_of(p0)?;
    let roots1 = rational_roots_of(p1)?;
    let lead0 = p0.last()?.clone();
    let lead1 = p1.last()?.clone();
    let power = |base: &Q, e: usize| (0..e).fold(Q::one(), |acc, _| acc * base);
    let c = -(&lead1 * power(&dq, roots1.len())) / (&lead0 * power(&dq, roots0.len()));
    let mut out = Vec::new();
    for (index, rho) in roots0.iter().enumerate() {
        // The same root twice would need logarithms.
        if roots0.iter().enumerate().any(|(i, r)| i != index && r == rho) {
            return None;
        }
        let a: Vec<Q> = roots1.iter().map(|s| (rho - s) / &dq).collect();
        let mut b = Vec::new();
        for (i, other) in roots0.iter().enumerate() {
            if i == index {
                continue;
            }
            let value = Q::one() + (rho - other) / &dq;
            if value.is_integer() && !value.is_positive() {
                return None;
            }
            b.push(value);
        }
        out.push(Pfq { rho: rho.clone(), a, b, c: c.clone(), d });
    }
    Some(out)
}

fn hypergeometric_node(
    cx: &mut Cx<'_>,
    x: NodeId,
    series: &Pfq,
) -> Option<NodeId> {
    let op = cx.graph.ops().lookup("hyperpfq")?;
    let node_of = |cx: &mut Cx<'_>, v: &Q| cx.graph.num(Number::rat(v.clone()));
    let mut arguments = vec![cx.graph.int(i64::try_from(series.a.len()).ok()?)];
    for v in series.a.iter().chain(&series.b) {
        arguments.push(node_of(cx, v));
    }
    let power = if series.d == 1 {
        x
    } else {
        let e = cx.graph.int(series.d);
        cx.graph.node(core::POW, &[x, e])
    };
    let c = node_of(cx, &series.c);
    arguments.push(mul(cx.graph, &[c, power]));
    let function = cx.graph.node(op, &arguments);
    if series.rho.is_zero() {
        return Some(function);
    }
    let rho = node_of(cx, &series.rho);
    let prefactor = cx.graph.node(core::POW, &[x, rho]);
    Some(mul(cx.graph, &[prefactor, function]))
}

/// `Σ_k (a)_k / (b)_k z^k / k!` (the generalised hypergeometric series),
/// for `ops` registration. Arguments: `p, a_1 … a_p, b_1 … b_q, z`.
#[must_use]
pub fn hyper_pfq(args: &[f64]) -> f64 {
    let Some(&p_f) = args.first() else { return f64::NAN };
    let p = p_f as usize;
    let Some(&z) = args.last() else { return f64::NAN };
    if args.len() < p + 2 {
        return f64::NAN;
    }
    let a = &args[1..=p];
    let b = &args[p + 1..args.len() - 1];
    let (mut term, mut sum) = (1.0_f64, 1.0_f64);
    for k in 0..5000_u32 {
        let kf = f64::from(k);
        let numerator: f64 = a.iter().map(|ai| ai + kf).product();
        let denominator: f64 = b.iter().map(|bj| bj + kf).product::<f64>() * (kf + 1.0);
        if denominator == 0.0 {
            return f64::NAN;
        }
        term *= numerator / denominator * z;
        sum += term;
        if term.abs() < 1e-17 * sum.abs() || term == 0.0 {
            break;
        }
    }
    sum
}

// ---------------------------------------------------------------------
// Hyperexponential solutions.

/// A hyperexponential solution `P(x) Π f^ρ exp(E)`.
#[derive(Clone)]
pub struct Hyper {
    /// The polynomial `P`.
    pub poly: Vec<Q>,
    /// The factors `f` with exponents `ρ`.
    pub factors: Vec<(Vec<Q>, Q)>,
    /// The exponent `E`.
    pub exponent: Rat,
}

impl Hyper {
    /// `y'/y`.
    fn log_derivative(&self) -> Option<Rat> {
        let mut total = self.exponent.derivative()?;
        for (f, rho) in &self.factors {
            total = total.add(&rat(&univariate::derivative(f), f)?.scale(rho)?)?;
        }
        total.add(&rat(&univariate::derivative(&self.poly), &self.poly)?)
    }
}

fn nullspace(
    mut rows: Vec<Vec<Q>>,
    columns: usize,
) -> Vec<Vec<Q>> {
    let mut pivots: Vec<usize> = Vec::new();
    let mut r = 0;
    for c in 0..columns {
        let Some(found) = (r..rows.len()).find(|&i| rows.get(i).and_then(|row| row.get(c)).is_some_and(|v| !v.is_zero())) else { continue };
        rows.swap(r, found);
        let pivot = rows.get(r).and_then(|row| row.get(c)).cloned().unwrap_or_else(Q::one);
        if let Some(row) = rows.get_mut(r) {
            for v in row.iter_mut() {
                *v /= &pivot;
            }
        }
        let pivot_row = rows.get(r).cloned().unwrap_or_default();
        for (i, row) in rows.iter_mut().enumerate() {
            if i == r {
                continue;
            }
            let factor = row.get(c).cloned().unwrap_or_else(Q::zero);
            if factor.is_zero() {
                continue;
            }
            for (v, p) in row.iter_mut().zip(&pivot_row) {
                *v -= &factor * p;
            }
        }
        pivots.push(c);
        r += 1;
    }
    let mut basis = Vec::new();
    for free in (0..columns).filter(|c| !pivots.contains(c)) {
        let mut vector = vec![Q::zero(); columns];
        if let Some(slot) = vector.get_mut(free) {
            *slot = Q::one();
        }
        for (row, &pivot) in pivots.iter().enumerate() {
            let value = rows.get(row).and_then(|r| r.get(free)).cloned().unwrap_or_else(Q::zero);
            if let Some(slot) = vector.get_mut(pivot) {
                *slot = -value;
            }
        }
        basis.push(vector);
    }
    basis
}

/// `L(P E)/E` for `P = x^j`: `Σ_i a_i (D + ω)^i x^j`.
fn apply_conjugated(
    op: &[Vec<Q>],
    omega: &Rat,
    j: usize,
) -> Option<Rat> {
    let mut monomial = vec![Q::zero(); j];
    monomial.push(Q::one());
    let mut current = rat(&monomial, &[Q::one()])?;
    let mut total = zero_rat()?;
    for (i, a) in op.iter().enumerate() {
        if i > 0 {
            current = current.derivative()?.add(&omega.mul(&current)?)?;
        }
        total = total.add(&rat(a, &[Q::one()])?.mul(&current)?)?;
    }
    Some(total)
}

/// The polynomials `P` of degree at most `d` with `L(P E) = 0`, for
/// `ω = E'/E`.
fn polynomial_solutions(
    op: &[Vec<Q>],
    omega: &Rat,
    d: usize,
) -> Option<Vec<Vec<Q>>> {
    let mut images = Vec::new();
    for j in 0..=d {
        images.push(apply_conjugated(op, omega, j)?);
    }
    // Common denominator.
    let mut common = vec![Q::one()];
    for image in &images {
        let g = univariate::gcd(&common, &image.den);
        let (quotient, _) = univariate::divrem(&univariate::mul(&common, &image.den), &g)?;
        common = quotient;
    }
    let mut numerators = Vec::new();
    for image in &images {
        let (factor, _) = univariate::divrem(&common, &image.den)?;
        numerators.push(univariate::mul(&image.num, &factor));
    }
    let height = numerators.iter().map(Vec::len).max().unwrap_or(0);
    let rows: Vec<Vec<Q>> = (0..height)
        .map(|m| numerators.iter().map(|n| n.get(m).cloned().unwrap_or_else(Q::zero)).collect())
        .collect();
    Some(nullspace(rows, d + 1))
}

/// Local candidates `(E, ρ)` at the roots of the irreducible factor `f`.
fn finite_candidates(
    op: &[Vec<Q>],
    f: &[Q],
) -> Vec<(Rat, Q)> {
    let mut out: Vec<(Rat, Q)> = Vec::new();
    let zero = zero_rat();
    let push = |out: &mut Vec<(Rat, Q)>, e: Rat, rho: Q| {
        if !out.iter().any(|(e2, r2)| *r2 == rho && e2.num == e.num && e2.den == e.den) {
            out.push((e, rho));
        }
    };
    if let Some(zero) = &zero {
        push(&mut out, zero.clone(), Q::zero());
    }
    for solution in formal_solutions(op, &FormalPoint::Finite(f.to_vec()), 2) {
        if solution.ramification != 1 {
            continue;
        }
        let Some(rho) = solution.field.as_rational(&solution.rho) else { continue };
        let mut exponent = zero_rat();
        let linear = f.len() == 2;
        let mut valid = true;
        for (s, kappa) in &solution.exponent {
            let Some(k) = solution.field.as_rational(kappa) else {
                valid = false;
                break;
            };
            if !linear {
                valid = false;
                break;
            }
            // κ (x - c)^{-s}, with c = -f_0.
            let mut den = vec![Q::one()];
            for _ in 0..*s {
                den = univariate::mul(&den, f);
            }
            exponent = exponent.and_then(|e| e.add(&rat(&[k], &den)?));
        }
        if let (true, Some(e)) = (valid, exponent) {
            push(&mut out, e, rho);
        }
    }
    out
}

/// Candidates `(E, ρ)` at infinity (`E` a polynomial in `x`, `ρ` the
/// exponent of `1/x`).
fn infinity_candidates(op: &[Vec<Q>]) -> Vec<(Rat, Q)> {
    let mut out: Vec<(Rat, Q)> = Vec::new();
    for solution in formal_solutions(op, &FormalPoint::Infinity, 2) {
        if solution.ramification != 1 {
            continue;
        }
        let Some(rho) = solution.field.as_rational(&solution.rho) else { continue };
        let mut coefficients: Vec<Q> = Vec::new();
        let mut valid = true;
        for (s, kappa) in &solution.exponent {
            let (Some(k), Ok(s)) = (solution.field.as_rational(kappa), usize::try_from(*s)) else {
                valid = false;
                break;
            };
            while coefficients.len() <= s {
                coefficients.push(Q::zero());
            }
            if let Some(slot) = coefficients.get_mut(s) {
                *slot = k;
            }
        }
        if !valid {
            continue;
        }
        let Some(e) = rat(&coefficients, &[Q::one()]) else { continue };
        if !out.iter().any(|(e2, r2)| *r2 == rho && e2.num == e.num) {
            out.push((e, rho));
        }
    }
    out
}

/// Hyperexponential solutions of `Σ a_i D^i`.
#[must_use]
pub fn hyperexponential_solutions(op: &[Vec<Q>]) -> Vec<Hyper> {
    let mut out = Vec::new();
    let Some(lead) = op.last() else { return out };
    let (_, factors) = univariate::factor(lead);
    let monic = |f: &[BigInt]| -> Vec<Q> {
        let l = Q::from_integer(f.last().cloned().unwrap_or_else(BigInt::one));
        f.iter().map(|c| Q::from_integer(c.clone()) / &l).collect()
    };
    let singular: Vec<Vec<Q>> = factors.iter().filter(|(f, _)| f.len() > 1).map(|(f, _)| monic(f)).collect();
    let local: Vec<Vec<(Rat, Q)>> = singular.iter().map(|f| finite_candidates(op, f)).collect();
    let infinity = infinity_candidates(op);
    if infinity.is_empty() {
        return out;
    }
    let total: usize = local.iter().map(Vec::len).product::<usize>() * infinity.len();
    if total > 5000 {
        return out;
    }
    let mut chosen = vec![0_usize; local.len()];
    for (e_inf, rho_inf) in &infinity {
        chosen.fill(0);
        loop {
            let mut shift = Q::zero();
            let mut exponent = e_inf.clone();
            let mut selected: Vec<(Vec<Q>, Q)> = Vec::new();
            let mut ok = true;
            for ((f, options), &i) in singular.iter().zip(local.iter()).zip(&chosen) {
                let Some((e, rho)) = options.get(i) else {
                    ok = false;
                    break;
                };
                shift += rho * q(i64::try_from(f.len() - 1).unwrap_or(0));
                match exponent.add(e) {
                    | Some(sum) => exponent = sum,
                    | None => ok = false,
                }
                selected.push((f.clone(), rho.clone()));
            }
            let degree = -rho_inf - shift;
            if ok && degree.is_integer() && !degree.is_negative() && degree.to_integer() <= BigInt::from(40) {
                let d = usize::try_from(degree.to_integer().to_i64().unwrap_or(0)).unwrap_or(0);
                if let Some(omega) = hyper_omega(&selected, &exponent)
                    && let Some(basis) = polynomial_solutions(op, &omega, d)
                {
                    for p in basis {
                        let poly = trim(p);
                        if !poly.is_empty() {
                            out.push(Hyper { poly, factors: selected.clone(), exponent: exponent.clone() });
                        }
                    }
                }
            }
            let mut position = 0;
            loop {
                if position == local.len() {
                    break;
                }
                let count = local.get(position).map_or(0, Vec::len);
                let Some(slot) = chosen.get_mut(position) else { break };
                *slot += 1;
                if *slot < count {
                    break;
                }
                *slot = 0;
                position += 1;
            }
            if position == local.len() {
                break;
            }
        }
        if !out.is_empty() {
            break;
        }
    }
    out
}

/// `ω = E' + Σ ρ f'/f`.
fn hyper_omega(
    factors: &[(Vec<Q>, Q)],
    exponent: &Rat,
) -> Option<Rat> {
    let mut total = exponent.derivative()?;
    for (f, rho) in factors {
        total = total.add(&rat(&univariate::derivative(f), f)?.scale(rho)?)?;
    }
    Some(total)
}

// ---------------------------------------------------------------------
// Right division and reduction of order.

/// The quotient `M` of `L = M (D - ω)` (coefficients of `M` as rational
/// functions, ascending order), when the division is exact.
fn right_quotient(
    op: &[Vec<Q>],
    omega: &Rat,
) -> Option<Vec<Rat>> {
    let n = op.len().checked_sub(1)?;
    let mut remainder: Vec<Rat> = op.iter().map(|a| rat(a, &[Q::one()])).collect::<Option<_>>()?;
    // ω, ω', … for the binomial expansion of D^{i-1} ∘ ω.
    let mut derivatives = vec![omega.clone()];
    for _ in 1..n {
        let next = derivatives.last()?.derivative()?;
        derivatives.push(next);
    }
    let mut quotient: Vec<Rat> = vec![zero_rat()?; n];
    for i in (1..=n).rev() {
        let lead = remainder.get(i)?.clone();
        *quotient.get_mut(i - 1)? = lead.clone();
        *remainder.get_mut(i)? = zero_rat()?;
        // r_{i-1-k} += lead * C(i-1, k) ω^(k).
        let mut binomial = BigInt::one();
        for k in 0..i {
            let term = lead.mul(derivatives.get(k)?)?.scale(&Q::from_integer(binomial.clone()))?;
            let slot = remainder.get_mut(i - 1 - k)?;
            *slot = slot.add(&term)?;
            binomial = binomial * BigInt::from(i - 1 - k) / BigInt::from(k + 1);
        }
    }
    remainder.first()?.num.is_empty().then_some(quotient)
}

fn hyper_node(
    cx: &mut Cx<'_>,
    x: NodeId,
    hyper: &Hyper,
) -> Option<NodeId> {
    let mut factors = Vec::new();
    for (f, rho) in &hyper.factors {
        if rho.is_zero() {
            continue;
        }
        let base = polynomial_node(cx, x, f);
        let exponent = cx.graph.num(Number::rat(rho.clone()));
        factors.push(cx.graph.node(core::POW, &[base, exponent]));
    }
    if !hyper.exponent.num.is_empty() {
        let (n, d) = (polynomial_node(cx, x, &hyper.exponent.num), polynomial_node(cx, x, &hyper.exponent.den));
        let ratio = super::div(cx.graph, n, d);
        factors.push(super::exp_of(cx, ratio)?);
    }
    factors.push(polynomial_node(cx, x, &hyper.poly));
    let product = mul(cx.graph, &factors);
    Some(cx.simplify(product))
}

/// The solver for linear equations of order three and more.
pub(super) fn high_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    let op = operator(cx, problem)?;
    let x = problem.x;
    let order = op.len() - 1;
    // Generalised hypergeometric functions.
    if let Some(series) = hypergeometric_solutions(&op) {
        let mut terms = Vec::new();
        for s in &series {
            let node = hypergeometric_node(cx, x, s)?;
            let c = problem.constant(cx.graph);
            terms.push(mul(cx.graph, &[c, node]));
        }
        if terms.len() == order {
            let solution = add(cx.graph, &terms);
            return Some(explicit(cx, problem, solution));
        }
    }
    // Hyperexponential solution and reduction of order.
    let found = hyperexponential_solutions(&op);
    let hyper = found.first()?;
    let y1 = hyper_node(cx, x, hyper)?;
    let omega = hyper.log_derivative()?;
    let quotient = right_quotient(&op, &omega)?;
    // Clear denominators of M.
    let mut common = vec![Q::one()];
    for m in &quotient {
        let g = univariate::gcd(&common, &m.den);
        common = univariate::divrem(&univariate::mul(&common, &m.den), &g)?.0;
    }
    let polynomials: Vec<Vec<Q>> = quotient
        .iter()
        .map(|m| {
            let (factor, _) = univariate::divrem(&common, &m.den)?;
            Some(trim(univariate::mul(&m.num, &factor)))
        })
        .collect::<Option<_>>()?;
    let z = fresh_function(cx.graph, "z", x);
    let derivatives = chain(cx.graph, z, x, order - 1)?;
    let mut terms = Vec::new();
    for (coefficient, derivative) in polynomials.iter().zip(&derivatives) {
        if coefficient.is_empty() {
            continue;
        }
        let c = polynomial_node(cx, x, coefficient);
        terms.push(mul(cx.graph, &[c, *derivative]));
    }
    let lhs = add(cx.graph, &terms);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[lhs, zero]);
    let (sub_problem, answer) = solve_equation(cx, equation, z, depth + 1)?;
    let rhs = super::reduce::explicit_rhs(cx.graph, answer, z)?;
    let count = super::reduce::constants_in(cx.graph, rhs).max(sub_problem.order());
    // The solutions z_j: set one constant to 1 and the others to 0.
    let mut basis = Vec::new();
    for j in 1..=count {
        let mut value = rhs;
        for k in 1..=count {
            let symbol = cx.graph.sym(&format!("C{k}"));
            let replacement = cx.graph.int(i64::from(k == j));
            value = cx.graph.substitute(value, symbol, replacement);
        }
        let value = cx.simplify(value);
        if !cx.is_zero(value) {
            basis.push(value);
        }
    }
    if basis.len() + 1 != order {
        return None;
    }
    let c1 = problem.constant(cx.graph);
    let mut total = vec![mul(cx.graph, &[c1, y1])];
    for z_j in basis {
        let integrand = super::div(cx.graph, z_j, y1);
        let integrand = cx.simplify(integrand);
        let inner = match integrate(cx, integrand, x) {
            | Some(found) => found,
            | None => {
                let op = cx.graph.ops().lookup("integral")?;
                cx.graph.node(op, &[integrand, x])
            },
        };
        let c = problem.constant(cx.graph);
        total.push(mul(cx.graph, &[c, y1, inner]));
    }
    let solution = add(cx.graph, &total);
    Some(explicit(cx, problem, solution))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qs(v: &[i64]) -> Vec<Q> {
        v.iter().map(|&c| q(c)).collect()
    }

    fn frac(
        a: i64,
        b: i64,
    ) -> Q {
        Q::new(BigInt::from(a), BigInt::from(b))
    }

    /// The D-form of a polynomial in θ = x d/dx with `x^k` shifts:
    /// θ^i = Σ S(i, m) x^m D^m (Stirling numbers of the second kind).
    fn theta_polynomial(
        coefficients: &[Q],
        shift: usize,
        op: &mut Vec<Vec<Q>>,
    ) {
        let n = coefficients.len().saturating_sub(1);
        // stirling[i][m]
        let mut stirling = vec![vec![Q::zero(); n + 2]; n + 2];
        stirling[0][0] = Q::one();
        for i in 1..=n {
            for m in 1..=i {
                let value = &stirling[i - 1][m - 1] + q(m as i64) * &stirling[i - 1][m];
                stirling[i][m] = value;
            }
        }
        while op.len() < n + 1 {
            op.push(Vec::new());
        }
        for (i, c) in coefficients.iter().enumerate() {
            for m in 0..=i {
                let s = stirling[i][m].clone();
                if s.is_zero() {
                    continue;
                }
                // c s x^m D^m times x^shift.
                let degree = m + shift;
                let slot = &mut op[m];
                while slot.len() <= degree {
                    slot.push(Q::zero());
                }
                slot[degree] += c * s;
            }
        }
    }

    #[test]
    fn hypergeometric_operator_is_matched() {
        // θ (θ - 1/2)(θ - 2/3) - x (θ + 1)(θ + 1/4)(θ + 1/5): 3F2 with
        // numerator parameters 1, 1/4, 1/5 and (b) = 1/2, 1/3 (for ρ = 0).
        let times = |a: &[Q], b: &[Q]| univariate::mul(a, b);
        let p0 = times(&times(&[Q::zero(), Q::one()], &[-frac(1, 2), Q::one()]), &[-frac(2, 3), Q::one()]);
        let p1 = times(&times(&[q(1), Q::one()], &[frac(1, 4), Q::one()]), &[frac(1, 5), Q::one()]);
        let p1 = p1.iter().map(|c| -c).collect::<Vec<_>>();
        let mut op = Vec::new();
        theta_polynomial(&p0, 0, &mut op);
        theta_polynomial(&p1, 1, &mut op);
        let series = hypergeometric_solutions(&op).expect("matched");
        assert_eq!(series.len(), 3);
        let zero = series.iter().find(|s| s.rho.is_zero()).expect("rho = 0");
        let mut a = zero.a.clone();
        a.sort();
        assert_eq!(a, vec![frac(1, 5), frac(1, 4), q(1)]);
        assert_eq!(zero.c, Q::one());
        assert_eq!(zero.d, 1);
        // b = 1 + (0 - 1/2), 1 + (0 - 2/3) = 1/2, 1/3.
        let mut b = zero.b.clone();
        b.sort();
        assert_eq!(b, vec![frac(1, 3), frac(1, 2)]);
        // The series satisfies the equation: x^3 y''' coefficients via D-form at x = 0.3.
        let numeric = |x: f64| -> f64 {
            let args: Vec<f64> = std::iter::once(zero.a.len() as f64)
                .chain(zero.a.iter().map(|v| v.to_f64().unwrap_or(f64::NAN)))
                .chain(zero.b.iter().map(|v| v.to_f64().unwrap_or(f64::NAN)))
                .chain(std::iter::once(x))
                .collect();
            hyper_pfq(&args)
        };
        let x = 0.3;
        let h = 1e-2;
        let samples: Vec<f64> = (-2..=2).map(|k| numeric(x + f64::from(k) * h)).collect();
        let d1 = (samples[0] - 8.0 * samples[1] + 8.0 * samples[3] - samples[4]) / (12.0 * h);
        let d2 = (-samples[0] + 16.0 * samples[1] - 30.0 * samples[2] + 16.0 * samples[3] - samples[4]) / (12.0 * h * h);
        let d3 = (-samples[0] + 2.0 * samples[1] - 2.0 * samples[3] + samples[4]) / (2.0 * h * h * h);
        let value = |c: &[Q]| c.iter().rev().fold(0.0, |acc, v| acc * x + v.to_f64().unwrap_or(f64::NAN));
        let residual = value(&op[0]) * samples[2] + value(&op[1]) * d1 + value(&op[2]) * d2 + value(&op[3]) * d3;
        assert!(residual.abs() < 2e-4, "residual {residual}");
    }

    #[test]
    fn hyperexponential_solution_of_a_third_order_operator() {
        // (D^2 - x)(D - 1) = D^3 - D^2 - x D + x: the solution e^x.
        let op = vec![qs(&[0, 1]), qs(&[0, -1]), qs(&[-1]), qs(&[1])];
        let found = hyperexponential_solutions(&op);
        let exp = found.iter().find(|h| h.exponent.num == qs(&[0, 1]) && h.poly == qs(&[1]));
        assert!(exp.is_some(), "e^x not found");
        // Right division by D - 1 leaves D^2 - x.
        let omega = rat(&qs(&[1]), &qs(&[1])).expect("omega");
        let quotient = right_quotient(&op, &omega).expect("exact");
        assert_eq!(quotient.len(), 3);
        assert_eq!(quotient[0].num, qs(&[0, -1]));
        assert!(quotient[1].num.is_empty());
        assert_eq!(quotient[2].num, qs(&[1]));
    }

    #[test]
    fn polynomial_and_algebraic_hyperexponential_solutions() {
        // x^3 y''' - 3 x^2 y'' + 6 x y' - 6 y = 0 is Euler: solutions x, x^2, x^3.
        let op = vec![qs(&[-6]), qs(&[0, 6]), qs(&[0, 0, -3]), qs(&[0, 0, 0, 1])];
        let found = hyperexponential_solutions(&op);
        assert!(!found.is_empty());
        for h in &found {
            // Each is c x^m.
            assert_eq!(h.exponent.num.len(), 0);
        }
    }
}
