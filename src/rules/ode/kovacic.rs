//! Kovacic's algorithm, case 1: hyperexponential solutions of linear
//! second-order equations with rational coefficients.
//!
//! `a2 y'' + a1 y' + a0 y = 0` is brought to the normal form `u'' = r u`
//! with `y = u exp(-½ ∫ a1/a2)`, `r = p²/4 + p'/2 - q` (`p = a1/a2`,
//! `q = a0/a2`). Kovacic's first case looks for a solution
//! `u = P exp(∫ω)` with `P` a polynomial and `ω` a rational function built
//! from the poles of `r`: at a pole of order `1` the contribution is
//! `1/(x - c)`; of order `2` it is `α/(x - c)` with `α = (1 ± sqrt(1 + 4b))/2`
//! (`b` the coefficient of `(x - c)^-2`); of even order `2ν` it is the
//! truncated Laurent series `±[√r]_c` plus `α/(x - c)`, `α = (±b/a + ν)/2`;
//! at infinity the same with the polynomial part of `√r`. For each choice of
//! signs the degree `d = α_∞ - Σ α_c` must be a non-negative integer, and
//! the monic `P` of degree `d` with `P'' + 2ωP' + (ω' + ω² - r)P = 0` is
//! found by exact linear algebra.
//!
//! The poles must be rational and the square roots that occur rational.
//! Cases 2 and 3 (solutions whose logarithmic derivative is algebraic of
//! degree 2, 4, 6, 12) are not implemented, so an equation with only such
//! solutions is not solved here. The second solution comes from reduction
//! of order in the caller.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::add;
use super::integrate;
use super::mul;
use super::neg;
use crate::graph::op::core;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::rational_function_in;
use crate::rules::poly::univariate;

type Q = BigRational;

/// A rational function `num / den` over `Q`.
#[derive(Clone)]
struct Rat {
    num: Vec<Q>,
    den: Vec<Q>,
}

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn trim(mut p: Vec<Q>) -> Vec<Q> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

impl Rat {
    fn new(
        num: &[Q],
        den: &[Q],
    ) -> Option<Self> {
        let (num, den) = (trim(num.to_vec()), trim(den.to_vec()));
        if den.is_empty() {
            return None;
        }
        if num.is_empty() {
            return Some(Self { num, den: vec![Q::one()] });
        }
        let g = univariate::gcd(&num, &den);
        let (num, _) = univariate::divrem(&num, &g)?;
        let (den, _) = univariate::divrem(&den, &g)?;
        let lead = den.last()?.clone();
        Some(Self { num: num.iter().map(|c| c / &lead).collect(), den: den.iter().map(|c| c / &lead).collect() })
    }

    fn mul(
        &self,
        o: &Self,
    ) -> Option<Self> {
        Self::new(&univariate::mul(&self.num, &o.num), &univariate::mul(&self.den, &o.den))
    }

    fn div(
        &self,
        o: &Self,
    ) -> Option<Self> {
        Self::new(&univariate::mul(&self.num, &o.den), &univariate::mul(&self.den, &o.num))
    }

    fn add(
        &self,
        o: &Self,
    ) -> Option<Self> {
        let n = univariate::add(&univariate::mul(&self.num, &o.den), &univariate::mul(&o.num, &self.den));
        Self::new(&n, &univariate::mul(&self.den, &o.den))
    }

    fn scale(
        &self,
        c: &Q,
    ) -> Option<Self> {
        Self::new(&self.num.iter().map(|x| x * c).collect::<Vec<_>>(), &self.den)
    }

    fn derivative(&self) -> Option<Self> {
        let n = univariate::sub(
            &univariate::mul(&univariate::derivative(&self.num), &self.den),
            &univariate::mul(&self.num, &univariate::derivative(&self.den)),
        );
        Self::new(&n, &univariate::mul(&self.den, &self.den))
    }

    fn eval(
        &self,
        x: &Q,
    ) -> Option<Q> {
        let d = univariate::eval(&self.den, x);
        (!d.is_zero()).then(|| univariate::eval(&self.num, x) / d)
    }
}

fn rational_sqrt(value: &Q) -> Option<Q> {
    if value.is_negative() {
        return None;
    }
    let (n, d) = (num_integer::Roots::sqrt(value.numer()), num_integer::Roots::sqrt(value.denom()));
    (&n * &n == *value.numer() && &d * &d == *value.denom()).then(|| Q::new(n, d))
}

/// Power series `numerator / denominator` to `n` terms.
fn series_quotient(
    numerator: &[Q],
    denominator: &[Q],
    n: usize,
) -> Option<Vec<Q>> {
    let lead = denominator.first()?.clone();
    if lead.is_zero() {
        return None;
    }
    let mut out: Vec<Q> = Vec::with_capacity(n);
    for k in 0..n {
        let mut value = numerator.get(k).cloned().unwrap_or_else(Q::zero);
        for j in 1..=k {
            value -= denominator.get(j).cloned().unwrap_or_else(Q::zero) * out.get(k - j)?;
        }
        out.push(value / &lead);
    }
    Some(out)
}

/// `p(t + h)`.
fn shift(
    p: &[Q],
    h: &Q,
) -> Vec<Q> {
    let mut work = p.to_vec();
    let n = work.len();
    for i in 0..n {
        for j in (i..n.saturating_sub(1)).rev() {
            let next = work.get(j + 1).cloned().unwrap_or_else(Q::zero);
            if let Some(slot) = work.get_mut(j) {
                *slot = &*slot + h * next;
            }
        }
    }
    work
}

/// One choice at a pole: coefficients of `(x - c)^-j`, `j = 1, 2, ...`
/// (index `j - 1`).
type PoleOption = Vec<Q>;

struct Pole {
    c: Q,
    options: Vec<PoleOption>,
}

/// One choice at infinity: the polynomial part (ascending) and `α_∞`.
struct InfinityOption {
    polynomial: Vec<Q>,
    alpha: Q,
}

fn square_root_pair(b: &Q) -> Option<(Q, Q)> {
    let s = rational_sqrt(&(Q::one() + q(4) * b))?;
    Some(((Q::one() + &s) / q(2), (Q::one() - &s) / q(2)))
}

/// The options at the pole `c` of order `order`, from the series of `r`.
fn pole_options(
    r: &Rat,
    c: &Q,
    order: usize,
) -> Option<Vec<PoleOption>> {
    if order == 1 {
        return Some(vec![vec![Q::one()]]);
    }
    if order % 2 == 1 {
        return None;
    }
    let nu = order / 2;
    // r(c + t) t^order = N(c + t) / D~(c + t).
    let numerator = shift(&r.num, c);
    let denominator = shift(&r.den, c);
    let rest: Vec<Q> = denominator.iter().skip(order).cloned().collect();
    if denominator.iter().take(order).any(|d| !d.is_zero()) {
        return None;
    }
    let s = series_quotient(&numerator, &rest, order + 1)?;
    if nu == 1 {
        let (plus, minus) = square_root_pair(s.first()?)?;
        return Some(vec![vec![plus], vec![minus]]);
    }
    let a = rational_sqrt(s.first()?)?;
    if a.is_zero() {
        return None;
    }
    // σ_i of the square root of the series.
    let mut sigma: Vec<Q> = vec![a.clone()];
    for i in 1..nu {
        let mut value = s.get(i)?.clone();
        for j in 1..i {
            value -= sigma.get(j)? * sigma.get(i - j)?;
        }
        sigma.push(value / (q(2) * &a));
    }
    // b = s_{nu-1} - Σ_{j=1}^{nu-2} σ_j σ_{nu-1-j}.
    let mut b = s.get(nu - 1)?.clone();
    for j in 1..nu - 1 {
        b -= sigma.get(j)? * sigma.get(nu - 1 - j)?;
    }
    let mut out = Vec::new();
    for sign in [1_i64, -1] {
        let sg = q(sign);
        let alpha = (&sg * &b / &a + Q::from_integer(BigInt::from(nu))) / q(2);
        // Coefficient of t^-j: j = nu - i for the σ_i, i = 0..nu-2.
        let mut coefficients = vec![Q::zero(); nu];
        for (i, sig) in sigma.iter().enumerate().take(nu - 1) {
            *coefficients.get_mut(nu - i - 1)? = &sg * sig;
        }
        *coefficients.first_mut()? = alpha;
        out.push(coefficients);
    }
    Some(out)
}

fn infinity_options(r: &Rat) -> Option<Vec<InfinityOption>> {
    let (dn, dd) = (r.num.len().checked_sub(1)?, r.den.len().checked_sub(1)?);
    let order = i64::try_from(dd).ok()? - i64::try_from(dn).ok()?;
    let lead = r.num.last()?.clone() / r.den.last()?;
    if order > 2 {
        return Some(vec![
            InfinityOption { polynomial: Vec::new(), alpha: Q::zero() },
            InfinityOption { polynomial: Vec::new(), alpha: Q::one() },
        ]);
    }
    if order == 2 {
        let (plus, minus) = square_root_pair(&lead)?;
        return Some(vec![
            InfinityOption { polynomial: Vec::new(), alpha: plus },
            InfinityOption { polynomial: Vec::new(), alpha: minus },
        ]);
    }
    if order % 2 != 0 {
        return None;
    }
    let nu = usize::try_from(-order / 2).ok()?;
    let a = rational_sqrt(&lead)?;
    if a.is_zero() {
        return None;
    }
    // r = x^(2 nu) Σ ρ_i x^-i: ratio of the reversed polynomials.
    let reversed = |p: &[Q]| -> Vec<Q> { p.iter().rev().cloned().collect() };
    let rho = series_quotient(&reversed(&r.num), &reversed(&r.den), nu + 2)?;
    let mut sigma: Vec<Q> = vec![a.clone()];
    for i in 1..=nu {
        let mut value = rho.get(i)?.clone();
        for j in 1..i {
            value -= sigma.get(j)? * sigma.get(i - j)?;
        }
        sigma.push(value / (q(2) * &a));
    }
    let mut b = rho.get(nu + 1)?.clone();
    for j in 1..=nu {
        b -= sigma.get(j)? * sigma.get(nu + 1 - j)?;
    }
    let mut out = Vec::new();
    for sign in [1_i64, -1] {
        let sg = q(sign);
        let alpha = (&sg * &b / &a - Q::from_integer(BigInt::from(nu))) / q(2);
        // x^(nu - i) with coefficient σ_i: ascending index nu - i.
        let mut polynomial = vec![Q::zero(); nu + 1];
        for (i, sig) in sigma.iter().enumerate() {
            *polynomial.get_mut(nu - i)? = &sg * sig;
        }
        out.push(InfinityOption { polynomial, alpha });
    }
    Some(out)
}

/// `ω(x)` and `ω'(x)` for the chosen options.
fn omega(
    polynomial: &[Q],
    poles: &[(&Q, &PoleOption)],
    x: &Q,
) -> Option<(Q, Q)> {
    let mut value = univariate::eval(polynomial, x);
    let mut slope = univariate::eval(&univariate::derivative(polynomial), x);
    for (c, coefficients) in poles {
        let t = x - *c;
        if t.is_zero() {
            return None;
        }
        let inverse = t.recip();
        let mut power = inverse.clone();
        for (j, coefficient) in coefficients.iter().enumerate() {
            // coefficient * t^-(j+1); derivative -(j+1) * coefficient * t^-(j+2).
            value += coefficient * &power;
            slope -= Q::from_integer(BigInt::from(j + 1)) * coefficient * &power * &inverse;
            power *= &inverse;
        }
    }
    Some((value, slope))
}

/// Gaussian elimination for `rows · p = rhs` (consistent, possibly
/// overdetermined); `None` when inconsistent or underdetermined.
fn solve_linear(
    mut rows: Vec<Vec<Q>>,
    mut rhs: Vec<Q>,
    unknowns: usize,
) -> Option<Vec<Q>> {
    let m = rows.len();
    let mut pivot_row = 0;
    let mut pivots = Vec::new();
    for col in 0..unknowns {
        let found = (pivot_row..m).find(|&r| rows.get(r).is_some_and(|row| row.get(col).is_some_and(|v| !v.is_zero())))?;
        rows.swap(pivot_row, found);
        rhs.swap(pivot_row, found);
        let pivot = rows.get(pivot_row)?.get(col)?.clone();
        for v in rows.get_mut(pivot_row)?.iter_mut() {
            *v /= &pivot;
        }
        *rhs.get_mut(pivot_row)? /= &pivot;
        let prow = rows.get(pivot_row)?.clone();
        let prhs = rhs.get(pivot_row)?.clone();
        for r in 0..m {
            if r == pivot_row {
                continue;
            }
            let factor = rows.get(r)?.get(col)?.clone();
            if factor.is_zero() {
                continue;
            }
            for (v, p) in rows.get_mut(r)?.iter_mut().zip(&prow) {
                *v -= &factor * p;
            }
            *rhs.get_mut(r)? -= &factor * &prhs;
        }
        pivots.push(col);
        pivot_row += 1;
    }
    // Remaining rows must be consistent.
    for r in pivot_row..m {
        if !rhs.get(r)?.is_zero() {
            return None;
        }
    }
    Some(pivots.iter().enumerate().map(|(i, _)| rhs.get(i).cloned().unwrap_or_else(Q::zero)).collect())
}

/// Solutions of `a2 y'' + a1 y' + a0 y = 0` by Kovacic's cases 1 and 2: one
/// hyperexponential solution, or the two solutions of case 2.
pub(super) fn solutions(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
) -> Vec<NodeId> {
    let Some((p, r)) = normal_form(cx, a, x) else { return Vec::new() };
    if let Some(y) = case_one(cx, x, &p, &r) {
        return vec![y];
    }
    case_two(cx, x, &p, &r).unwrap_or_default()
}

/// `p = a1/a2` and `r = p²/4 + p'/2 - a0/a2`.
fn normal_form(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
) -> Option<(Rat, Rat)> {
    let mut fractions = Vec::new();
    for &c in a {
        let (n, d) = rational_function_in(cx.graph, c, x)?;
        fractions.push(Rat::new(&n, &d)?);
    }
    let (a0, a1, a2) = (fractions.first()?, fractions.get(1)?, fractions.get(2)?);
    let p = a1.div(a2)?;
    let qq = a0.div(a2)?;
    let r = p.mul(&p)?.scale(&Q::new(BigInt::one(), BigInt::from(4)))?
        .add(&p.derivative()?.scale(&Q::new(BigInt::one(), BigInt::from(2)))?)?
        .add(&qq.scale(&q(-1))?)?;
    if r.num.is_empty() {
        return None;
    }
    Some((p, r))
}

/// Case 1: `u = P exp(∫ω)` with `ω` rational.
fn case_one(
    cx: &mut Cx<'_>,
    x: NodeId,
    p: &Rat,
    r: &Rat,
) -> Option<NodeId> {
    // Poles.
    let (_, factors) = univariate::factor(&r.den);
    let mut poles: Vec<Pole> = Vec::new();
    for (factor, multiplicity) in factors {
        let [c0, c1] = factor.as_slice() else { return None };
        let c = Q::new(-c0.clone(), c1.clone());
        let order = usize::try_from(multiplicity).ok()?;
        poles.push(Pole { options: pole_options(r, &c, order)?, c });
    }
    let infinity = infinity_options(r)?;
    // Enumerate the families.
    let total: usize = poles.iter().map(|p| p.options.len()).product::<usize>() * infinity.len();
    if total > 4096 {
        return None;
    }
    let mut chosen = vec![0_usize; poles.len()];
    for inf in &infinity {
        chosen.fill(0);
        loop {
            let alpha_sum: Q = poles
                .iter()
                .zip(&chosen)
                .map(|(p, &i)| p.options.get(i).and_then(|o| o.first()).cloned().unwrap_or_else(Q::zero))
                .sum();
            let d = &inf.alpha - alpha_sum;
            if d.denom().is_one() && !d.is_negative() && d.to_integer() <= BigInt::from(60) {
                let degree = usize::try_from(d.to_integer().to_i64()?).ok()?;
                let selection: Vec<(&Q, &PoleOption)> =
                    poles.iter().zip(&chosen).filter_map(|(p, &i)| Some((&p.c, p.options.get(i)?))).collect();
                if let Some(polynomial) = find_polynomial(r, &inf.polynomial, &selection, degree)
                    && let Some(node) = build(cx, x, p, &inf.polynomial, &selection, &polynomial)
                {
                    return Some(node);
                }
            }
            // Next combination.
            let mut position = 0;
            loop {
                if position == poles.len() {
                    break;
                }
                let count = poles.get(position)?.options.len();
                let slot = chosen.get_mut(position)?;
                *slot += 1;
                if *slot < count {
                    break;
                }
                *slot = 0;
                position += 1;
            }
            if position == poles.len() {
                break;
            }
        }
    }
    None
}

/// The monic polynomial `P` of degree `d` with
/// `P'' + 2ωP' + (ω' + ω² - r)P = 0`.
fn find_polynomial(
    r: &Rat,
    omega_polynomial: &[Q],
    poles: &[(&Q, &PoleOption)],
    d: usize,
) -> Option<Vec<Q>> {
    let size: usize = poles.iter().map(|(_, o)| o.len()).sum::<usize>() + omega_polynomial.len();
    let samples = d + 12 + 2 * size;
    let mut points = Vec::new();
    let mut k = 0_i64;
    while points.len() < samples {
        let x0 = q(2 + k) + Q::new(BigInt::one(), BigInt::from(3));
        if poles.iter().all(|(c, _)| x0 != **c) && r.eval(&x0).is_some() {
            points.push(x0);
        }
        k += 1;
    }
    let mut rows: Vec<Vec<Q>> = Vec::new();
    let mut rhs: Vec<Q> = Vec::new();
    for x0 in &points {
        let (w, dw) = omega(omega_polynomial, poles, x0)?;
        let residual_factor = &dw + &w * &w - r.eval(x0)?;
        // L[x^i](x0).
        let apply = |i: usize| -> Q {
            let ii = Q::from_integer(BigInt::from(i));
            let mut value = &residual_factor * pow(x0, i);
            if i >= 1 {
                value += q(2) * &w * &ii * pow(x0, i - 1);
            }
            if i >= 2 {
                value += &ii * (&ii - Q::one()) * pow(x0, i - 2);
            }
            value
        };
        rows.push((0..d).map(&apply).collect());
        rhs.push(-apply(d));
    }
    let mut coefficients = if d == 0 {
        // Only the consistency of L[1] = 0.
        if rhs.iter().any(|v| !v.is_zero()) {
            return None;
        }
        Vec::new()
    } else {
        solve_linear(rows, rhs, d)?
    };
    coefficients.push(Q::one());
    Some(coefficients)
}

fn pow(
    x: &Q,
    n: usize,
) -> Q {
    (0..n).fold(Q::one(), |acc, _| acc * x)
}

/// `y1 = P exp(∫ω) exp(-½ ∫p)` as a term.
fn build(
    cx: &mut Cx<'_>,
    x: NodeId,
    p: &Rat,
    omega_polynomial: &[Q],
    poles: &[(&Q, &PoleOption)],
    polynomial: &[Q],
) -> Option<NodeId> {
    let ln = cx.graph.ops().lookup("ln")?;
    let node_of = |cx: &mut Cx<'_>, v: &Q| cx.graph.num(Number::rat(v.clone()));
    let mut exponent_terms: Vec<NodeId> = Vec::new();
    for (i, c) in omega_polynomial.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let coefficient = node_of(cx, &(c / Q::from_integer(BigInt::from(i + 1))));
        let e = cx.graph.int(i64::try_from(i + 1).ok()?);
        let power = cx.graph.node(core::POW, &[x, e]);
        exponent_terms.push(mul(cx.graph, &[coefficient, power]));
    }
    for (c, coefficients) in poles {
        let shifted = {
            let minus_c = node_of(cx, &-(*c).clone());
            cx.graph.node(core::ADD, &[x, minus_c])
        };
        for (j, coefficient) in coefficients.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            if j == 0 {
                let log = cx.graph.node(ln, &[shifted]);
                let cn = node_of(cx, coefficient);
                exponent_terms.push(mul(cx.graph, &[cn, log]));
            } else {
                // coefficient t^-(j+1) integrates to coefficient t^-j / (-j).
                let scale = coefficient / Q::from_integer(BigInt::from(-i64::try_from(j).ok()?));
                let cn = node_of(cx, &scale);
                let e = cx.graph.int(-i64::try_from(j).ok()?);
                let power = cx.graph.node(core::POW, &[shifted, e]);
                exponent_terms.push(mul(cx.graph, &[cn, power]));
            }
        }
    }
    // -1/2 ∫ p.
    let p_node = {
        let mut gens = crate::rules::poly::repr::Gens::default();
        let gx = gens.index(cx.graph, x);
        let to = |cx: &mut Cx<'_>, gens: &crate::rules::poly::repr::Gens, v: &[Q]| {
            let numbers: Vec<Number> = v.iter().cloned().map(Number::rat).collect();
            crate::rules::poly::repr::to_term(cx.graph, gens, &crate::rules::poly::repr::Poly::from_univariate(gx, &numbers))
        };
        let (n, d) = (to(cx, &gens, &p.num), to(cx, &gens, &p.den));
        super::div(cx.graph, n, d)
    };
    let integral = integrate(cx, p_node, x)?;
    let half = cx.graph.num(Number::fraction(-1, 2)?);
    exponent_terms.push(mul(cx.graph, &[half, integral]));
    let exponent = add(cx.graph, &exponent_terms);
    let exponential = super::exp_of(cx, exponent)?;
    // P(x).
    let mut terms = Vec::new();
    for (i, c) in polynomial.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let cn = node_of(cx, c);
        let e = cx.graph.int(i64::try_from(i).ok()?);
        let power = if i == 0 { cx.graph.int(1) } else { cx.graph.node(core::POW, &[x, e]) };
        terms.push(mul(cx.graph, &[cn, power]));
    }
    let p_poly = add(cx.graph, &terms);
    let product = mul(cx.graph, &[p_poly, exponential]);
    let _ = neg;
    Some(cx.simplify(product))
}

/// The integers among `2 + k sqrt(1 + 4b)`, `k = 0, ±2` (the exponent set
/// at a pole of order two).
fn order_two_set(b: &Q) -> Vec<i64> {
    let mut out = vec![2_i64];
    if let Some(s) = rational_sqrt(&(Q::one() + q(4) * b)) {
        for k in [2_i64, -2] {
            let value = q(2) + q(k) * &s;
            if value.denom().is_one()
                && let Some(v) = value.to_integer().to_i64()
                && !out.contains(&v)
            {
                out.push(v);
            }
        }
    }
    out
}

/// Case 2: `u = exp(∫ω)` with `ω` a root of `ω² - φ ω + (φ'/2 + φ²/2 - r) = 0`,
/// `φ = θ + P'/P` rational, `θ = ½ Σ e_c/(x - c)`. The two roots give two
/// independent solutions.
fn case_two(
    cx: &mut Cx<'_>,
    x: NodeId,
    p: &Rat,
    r: &Rat,
) -> Option<Vec<NodeId>> {
    // Poles with their orders and the coefficient of (x - c)^-2.
    let (_, factors) = univariate::factor(&r.den);
    let mut sets: Vec<(Q, Vec<i64>)> = Vec::new();
    let mut interesting = false;
    for (factor, multiplicity) in factors {
        let [c0, c1] = factor.as_slice() else { return None };
        let c = Q::new(-c0.clone(), c1.clone());
        let order = usize::try_from(multiplicity).ok()?;
        let set = match order {
            | 1 => vec![4],
            | 2 => {
                interesting = true;
                let numerator = shift(&r.num, &c);
                let denominator = shift(&r.den, &c);
                let rest: Vec<Q> = denominator.iter().skip(2).cloned().collect();
                let s = series_quotient(&numerator, &rest, 1)?;
                order_two_set(s.first()?)
            },
            | n => {
                interesting = true;
                vec![i64::try_from(n).ok()?]
            },
        };
        sets.push((c, set));
    }
    if !interesting {
        return None;
    }
    let (dn, dd) = (r.num.len().checked_sub(1)?, r.den.len().checked_sub(1)?);
    let order_inf = i64::try_from(dd).ok()? - i64::try_from(dn).ok()?;
    let infinity_set: Vec<i64> = match order_inf.cmp(&2) {
        | std::cmp::Ordering::Greater => vec![0, 2, 4],
        | std::cmp::Ordering::Equal => order_two_set(&(r.num.last()?.clone() / r.den.last()?)),
        | std::cmp::Ordering::Less => vec![order_inf],
    };
    let total: usize = sets.iter().map(|(_, s)| s.len()).product::<usize>() * infinity_set.len();
    if total > 4096 {
        return None;
    }
    let mut chosen = vec![0_usize; sets.len()];
    for &e_inf in &infinity_set {
        chosen.fill(0);
        loop {
            let sum: i64 = sets.iter().zip(&chosen).map(|((_, s), &i)| s.get(i).copied().unwrap_or(0)).sum();
            let twice = e_inf - sum;
            if twice >= 0 && twice % 2 == 0 && twice / 2 <= 40 {
                let degree = usize::try_from(twice / 2).ok()?;
                let exponents: Vec<(&Q, i64)> = sets.iter().zip(&chosen).map(|((c, s), &i)| (c, s.get(i).copied().unwrap_or(0))).collect();
                if let Some(polynomial) = find_third_order(r, &exponents, degree)
                    && let Some(found) = build_case_two(cx, x, p, r, &exponents, &polynomial)
                {
                    return Some(found);
                }
            }
            let mut position = 0;
            loop {
                if position == sets.len() {
                    break;
                }
                let count = sets.get(position)?.1.len();
                let slot = chosen.get_mut(position)?;
                *slot += 1;
                if *slot < count {
                    break;
                }
                *slot = 0;
                position += 1;
            }
            if position == sets.len() {
                break;
            }
        }
    }
    None
}

/// `θ`, `θ'`, `θ''` at `x0`.
fn theta_values(
    exponents: &[(&Q, i64)],
    x0: &Q,
) -> Option<(Q, Q, Q)> {
    let (mut t, mut t1, mut t2) = (Q::zero(), Q::zero(), Q::zero());
    for (c, e) in exponents {
        let d = x0 - *c;
        if d.is_zero() {
            return None;
        }
        let inv = d.recip();
        let e = q(*e);
        t += &e * &inv / q(2);
        t1 -= &e * &inv * &inv / q(2);
        t2 += &e * &inv * &inv * &inv;
    }
    Some((t, t1, t2))
}

/// The monic `P` of degree `d` with
/// the third-order equation of Kovacic's case 2 (`P` with three derivatives).
fn find_third_order(
    r: &Rat,
    exponents: &[(&Q, i64)],
    d: usize,
) -> Option<Vec<Q>> {
    let dr = r.derivative()?;
    let samples = d + 14 + 3 * exponents.len();
    let mut points = Vec::new();
    let mut k = 0_i64;
    while points.len() < samples {
        let x0 = q(2 + k) + Q::new(BigInt::one(), BigInt::from(3));
        if exponents.iter().all(|(c, _)| x0 != **c) && r.eval(&x0).is_some() {
            points.push(x0);
        }
        k += 1;
    }
    let mut rows: Vec<Vec<Q>> = Vec::new();
    let mut rhs: Vec<Q> = Vec::new();
    for x0 in &points {
        let (t, t1, t2) = theta_values(exponents, x0)?;
        let (rv, drv) = (r.eval(x0)?, dr.eval(x0)?);
        let c2 = q(3) * &t * &t + q(3) * &t1 - q(4) * &rv;
        let c3 = &t2 + q(3) * &t * &t1 + &t * &t * &t - q(4) * &rv * &t - q(2) * &drv;
        let apply = |i: usize| -> Q {
            let ii = Q::from_integer(BigInt::from(i));
            let mut value = &c3 * pow(x0, i);
            if i >= 1 {
                value += &c2 * &ii * pow(x0, i - 1);
            }
            if i >= 2 {
                value += q(3) * &t * &ii * (&ii - Q::one()) * pow(x0, i - 2);
            }
            if i >= 3 {
                value += &ii * (&ii - Q::one()) * (&ii - q(2)) * pow(x0, i - 3);
            }
            value
        };
        rows.push((0..d).map(&apply).collect());
        rhs.push(-apply(d));
    }
    let mut coefficients = if d == 0 {
        if rhs.iter().any(|v| !v.is_zero()) {
            return None;
        }
        Vec::new()
    } else {
        solve_linear(rows, rhs, d)?
    };
    coefficients.push(Q::one());
    Some(coefficients)
}

/// Both solutions of case 2 as terms.
fn build_case_two(
    cx: &mut Cx<'_>,
    x: NodeId,
    p: &Rat,
    r: &Rat,
    exponents: &[(&Q, i64)],
    polynomial: &[Q],
) -> Option<Vec<NodeId>> {
    let ln = cx.graph.ops().lookup("ln")?;
    // φ = θ + P'/P and R = 4 r - φ² - 2 φ'.
    let mut theta = Rat::new(&[], &[Q::one()])?;
    for (c, e) in exponents {
        let term = Rat::new(&[q(*e) / q(2)], &[-(*c).clone(), Q::one()])?;
        theta = theta.add(&term)?;
    }
    let ratio = Rat::new(&univariate::derivative(polynomial), polynomial)?;
    let phi = theta.add(&ratio)?;
    let discriminant = r
        .scale(&q(4))?
        .add(&phi.mul(&phi)?.scale(&q(-1))?)?
        .add(&phi.derivative()?.scale(&q(-2))?)?;
    let to_node = |cx: &mut Cx<'_>, v: &[Q]| {
        let mut gens = crate::rules::poly::repr::Gens::default();
        let gx = gens.index(cx.graph, x);
        let numbers: Vec<Number> = v.iter().cloned().map(Number::rat).collect();
        crate::rules::poly::repr::to_term(cx.graph, &gens, &crate::rules::poly::repr::Poly::from_univariate(gx, &numbers))
    };
    let node_of = |cx: &mut Cx<'_>, v: &Q| cx.graph.num(Number::rat(v.clone()));
    // Common part of the exponent: Σ e_c/4 ln(x - c) + ½ ln P - ½ ∫ p.
    let mut common: Vec<NodeId> = Vec::new();
    for (c, e) in exponents {
        let minus_c = node_of(cx, &-(*c).clone());
        let shifted = cx.graph.node(core::ADD, &[x, minus_c]);
        let log = cx.graph.node(ln, &[shifted]);
        let coefficient = node_of(cx, &(q(*e) / q(4)));
        common.push(mul(cx.graph, &[coefficient, log]));
    }
    if polynomial.len() > 1 {
        let poly_node = to_node(cx, polynomial);
        let log = cx.graph.node(ln, &[poly_node]);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        common.push(mul(cx.graph, &[half, log]));
    }
    let p_node = {
        let (n, d) = (to_node(cx, &p.num), to_node(cx, &p.den));
        super::div(cx.graph, n, d)
    };
    let integral = integrate(cx, p_node, x)?;
    let minus_half = cx.graph.num(Number::fraction(-1, 2)?);
    common.push(mul(cx.graph, &[minus_half, integral]));
    // ± ½ ∫ sqrt(R), with sqrt written factor by factor.
    let root = {
        let (nc, nf) = univariate::factor(&discriminant.num);
        let (dc, df) = univariate::factor(&discriminant.den);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let lead = node_of(cx, &(nc / dc));
        let mut factors = vec![cx.graph.node(core::POW, &[lead, half])];
        for (list, sign) in [(nf, 1_i64), (df, -1)] {
            for (f, multiplicity) in list {
                let polynomial: Vec<Q> = f.iter().cloned().map(Q::from_integer).collect();
                let base = to_node(cx, &polynomial);
                let exponent = node_of(cx, &Q::new(BigInt::from(sign * i64::from(multiplicity)), BigInt::from(2)));
                factors.push(cx.graph.node(core::POW, &[base, exponent]));
            }
        }
        let product = mul(cx.graph, &factors);
        cx.simplify(product)
    };
    let root_integral = match integrate(cx, root, x) {
        | Some(found) => found,
        | None => {
            let op = cx.graph.ops().lookup("integral")?;
            cx.graph.node(op, &[root, x])
        },
    };
    let mut out = Vec::new();
    for sign in [1_i64, -1] {
        let coefficient = cx.graph.num(Number::rat(Q::new(BigInt::from(sign), BigInt::from(2))));
        let mut terms = common.clone();
        terms.push(mul(cx.graph, &[coefficient, root_integral]));
        let exponent = add(cx.graph, &terms);
        let y = super::exp_of(cx, exponent)?;
        out.push(cx.simplify(y));
    }
    let _ = neg;
    Some(out)
}
