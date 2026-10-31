//! Combinatorics: factorials, binomial coefficients, counting sequences and
//! the identities among them.
//!
//! Operators are reduced by exact kernels on literal integers. Results grow
//! super-exponentially, so each kernel declines inputs whose answer would
//! run to tens of thousands of bits; the term then stays symbolic.
//!
//! Core operators: `factorial`, `double_factorial`, `binomial`,
//! `multinomial`, `fibonacci`, `lucas`, `catalan`, `bell`, `stirling1`,
//! `stirling2`, `partitions`, `derangements`, `harmonic`, `permutations`,
//! `falling`, `rising`, `inclusion_exclusion`, `rsolve` and `gf_coeffs`.
//! The module also provides `lah`, `eulerian`, `narayana`, `motzkin`,
//! `schroeder`, `fubini`, `central_binomial`, `ballot`, `involutions`,
//! `rencontres`, `partitions_k`, `partitions_max`, `partitions_distinct`,
//! `partition_list`, `compositions`, `weak_compositions`, `necklaces`,
//! `bracelets`, `bernoulli_number`, `euler_number`, `zigzag`,
//! `permutation_rank`, `permutation_unrank`, `inversions`,
//! `permutation_sign`, `cycle_type`, `permutation_order`, `hook_lengths`
//! and `syt_count`.

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use statrs::function::gamma::ln_gamma;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimRule;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::graph::rule::Installer;


use super::arith::arith;
use super::number_theory::Value;
use super::number_theory::attach_dimensions;
use super::number_theory::exact;
use super::number_theory::exact_lists;
use super::number_theory::procedure;
use super::number_theory::residual;
use super::poly::best;
use super::poly::ratio;
use super::poly::repr::Gens;
use super::poly::repr::Limits;
use super::poly::repr::from_term;
use super::poly::univariate;

/// Largest supported argument of `factorial` (a 54,000-bit result).
const MAX_FACTORIAL: u64 = 5000;
/// Bound on the bits a multiplicative kernel may accumulate.
const MAX_BITS: u64 = 60_000;

/// The combinatorics rule set.
#[must_use]
pub fn combinatorics() -> RuleSet {
    RuleSet::new("combinatorics", install)
        .needs(arith())
        .needs(super::calculus::calculus())
}

/// The value of `v` if it is a natural number not above `cap`.
fn nat(
    v: &BigInt,
    cap: u64,
) -> Option<u64> {
    v.to_u64().filter(|&x| x <= cap)
}

fn factorial(n: u64) -> BigInt {
    (2..=n).fold(BigInt::one(), |acc, k| acc * k)
}

fn double_factorial(n: u64) -> BigInt {
    (1..=n)
        .rev()
        .step_by(2)
        .fold(BigInt::one(), |acc, k| acc * k)
}

/// `binomial(n, k)` for `0 <= k <= n`, by the multiplicative formula.
pub(crate) fn choose(
    n: &BigInt,
    k: u64,
) -> Option<BigInt> {
    let k = n.to_u64().map_or(k, |n| k.min(n.saturating_sub(k)));
    if k.saturating_mul(n.bits()) > MAX_BITS {
        return None;
    }
    let base = n - k;
    let mut acc = BigInt::one();
    for i in 1..=k {
        // Each partial product is itself a binomial coefficient, so the
        // division is exact.
        acc = acc * (&base + i) / i;
    }
    Some(acc)
}

/// Generalised `binomial(n, k)` for every integer `n`; zero for `k < 0`.
fn binomial(
    n: &BigInt,
    k: &BigInt,
) -> Option<BigInt> {
    if k.is_negative() {
        return Some(BigInt::zero());
    }
    if n.is_negative() {
        // (-m choose k) = (-1)^k (m + k - 1 choose k)
        let k = nat(k, MAX_FACTORIAL)?;
        let c = choose(&(BigInt::from(k) - n - 1), k)?;
        return Some(if k % 2 == 0 { c } else { -c });
    }
    if k > n {
        return Some(BigInt::zero());
    }
    choose(n, k.to_u64()?)
}

/// `x (x - step) ... ` over `k` factors.
fn product_run(
    x: &BigInt,
    k: &BigInt,
    step: i8,
) -> Option<BigInt> {
    let k = nat(k, MAX_FACTORIAL)?;
    if k.saturating_mul(
        x.bits()
            .max(1)
            .saturating_add(u64::from(64 - k.leading_zeros())),
    ) > MAX_BITS
    {
        return None;
    }
    let mut acc = BigInt::one();
    let mut term = x.clone();
    for _ in 0..k {
        acc *= &term;
        term += i32::from(step);
    }
    Some(acc)
}

/// `(F(n), F(n + 1))` by fast doubling.
fn fibonacci_pair(n: u64) -> (BigInt, BigInt) {
    if n == 0 {
        return (BigInt::zero(), BigInt::one());
    }
    let (a, b) = fibonacci_pair(n / 2);
    let c = &a * (&b * 2 - &a);
    let d = &a * &a + &b * &b;
    if n.is_multiple_of(2) {
        (c, d)
    } else {
        (d.clone(), c + d)
    }
}

fn bell(n: u64) -> BigInt {
    // Bell's triangle: each row starts with the end of the previous one.
    let mut row = vec![BigInt::one()];
    for _ in 0..n {
        let mut next = vec![row.last().cloned().unwrap_or_default()];
        for above in &row {
            let value = next.last().cloned().unwrap_or_default() + above;
            next.push(value);
        }
        row = next;
    }
    row.into_iter().next().unwrap_or_default()
}

/// Row `n` of the triangle `t(n, j) = t(n-1, j-1) + weight(n, j) * t(n-1, j)`
/// started from `t(0, 0) = 1`.
fn stirling_row(
    n: u64,
    weight: fn(u64, u64) -> u64,
) -> Vec<BigInt> {
    let mut row = vec![BigInt::one()];
    for m in 1..=n {
        let mut next = vec![BigInt::zero(); row.len() + 1];
        for (j, value) in row.iter().enumerate() {
            let w = weight(m, u64::try_from(j).unwrap_or(0));
            // Contribution to t(m, j + 1) and to t(m, j).
            if let Some(slot) = next.get_mut(j + 1) {
                *slot += value;
            }
            if let Some(slot) = next.get_mut(j) {
                *slot += value * w;
            }
        }
        row = next;
    }
    row
}

fn partitions(n: usize) -> BigInt {
    // Euler's pentagonal number recurrence.
    let mut p = vec![BigInt::one()];
    for m in 1..=n {
        let mut total = BigInt::zero();
        let mut j = 1_usize;
        loop {
            let (g1, g2) = (j * (3 * j - 1) / 2, j * (3 * j + 1) / 2);
            if g1 > m {
                break;
            }
            let mut part = p[m - g1].clone();
            if g2 <= m {
                part += &p[m - g2];
            }
            if j % 2 == 1 {
                total += part;
            } else {
                total -= part;
            }
            j += 1;
        }
        p.push(total);
    }
    p.pop().unwrap_or_default()
}

fn derangements(n: u64) -> BigInt {
    let (mut prev, mut cur) = (BigInt::one(), BigInt::zero());
    if n == 0 {
        return prev;
    }
    for k in 2..=n {
        let next = (&prev + &cur) * (k - 1);
        prev = cur;
        cur = next;
    }
    cur
}

fn harmonic(n: u64) -> BigRational {
    (1..=n).fold(BigRational::zero(), |acc, k| {
        acc + BigRational::new(BigInt::one(), BigInt::from(k))
    })
}

fn multinomial(args: &[Vec<BigInt>]) -> Option<Value> {
    let [parts] = args else {
        return None;
    };
    let mut total = 0_u64;
    let mut acc = BigInt::one();
    for part in parts {
        let k = nat(part, MAX_FACTORIAL)?;
        total = total.checked_add(k).filter(|&t| t <= MAX_FACTORIAL)?;
        acc *= choose(&BigInt::from(total), k)?;
    }
    Some(Value::Int(acc))
}

/// `Gamma(n + 1)` on reals: the value `n!` extends to.
fn factorial_eval(args: &[f64]) -> f64 {
    match args.first() {
        | Some(&x) if x > -1.0 => ln_gamma(x + 1.0).exp(),
        | _ => f64::NAN,
    }
}

/// `Gamma(n + 1) / (Gamma(k + 1) Gamma(n - k + 1))`, zero where the last
/// factor has a pole, and the usual extension to negative integers `n`.
fn binomial_eval(args: &[f64]) -> f64 {
    let (Some(&n), Some(&k)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if n.is_nan() || k.is_nan() {
        return f64::NAN;
    }
    if k < 0.0 {
        return 0.0;
    }
    if n < 0.0 {
        if n.fract() != 0.0 || k.fract() != 0.0 {
            return f64::NAN;
        }
        let sign = if k.rem_euclid(2.0) == 0.0 {
            1.0
        } else {
            -1.0
        };
        return sign * binomial_eval(&[k - n - 1.0, k]);
    }
    let m = n - k + 1.0;
    if m <= 0.0 && m.fract() == 0.0 {
        return 0.0;
    }
    if m < 0.0 {
        return f64::NAN;
    }
    (ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(m)).exp()
}

/// `Gamma(n + 1) / Gamma(n - k + 1)`, zero where the denominator has a pole.
fn permutations_eval(args: &[f64]) -> f64 {
    let (Some(&n), Some(&k)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    let m = n - k + 1.0;
    if n < 0.0 || m < 0.0 && m.fract() != 0.0 {
        return f64::NAN;
    }
    if m <= 0.0 {
        return 0.0;
    }
    (ln_gamma(n + 1.0) - ln_gamma(m)).exp()
}

/// Most terms of a generating function that are expanded.
const MAX_SERIES: u64 = 2000;
/// Most derivatives taken when a series has to be found by differentiation.
const MAX_TAYLOR: u64 = 16;
/// Highest order of a recurrence that is solved.
const MAX_ORDER: usize = 24;
/// Largest exponent of a rational power computed exactly.
const MAX_POWER: i64 = 256;

fn int_literal(
    graph: &Graph,
    node: NodeId,
) -> Option<BigInt> {
    match graph.number_of(node)? {
        | Number::Int(v) => Some(v.clone()),
        | _ => None,
    }
}

fn rat_pow(
    base: &BigRational,
    exponent: i64,
) -> Option<BigRational> {
    if exponent.abs() > MAX_POWER || (exponent < 0 && base.is_zero()) {
        return None;
    }
    let magnitude = (0..exponent.unsigned_abs()).fold(BigRational::one(), |acc, _| acc * base);
    Some(if exponent < 0 {
        magnitude.recip()
    } else {
        magnitude
    })
}

/// `p(n + shift)` for the polynomial with ascending coefficients `p`.
fn poly_shift(
    p: &[BigRational],
    shift: &BigRational,
) -> Option<Vec<BigRational>> {
    let mut out = vec![BigRational::zero(); p.len()];
    for (k, pk) in p.iter().enumerate() {
        for l in 0..=k {
            let c = choose(&BigInt::from(k), u64::try_from(l).ok()?)?;
            let power = rat_pow(shift, i64::try_from(k - l).ok()?)?;
            if let Some(slot) = out.get_mut(l) {
                *slot += pk * BigRational::from_integer(c) * power;
            }
        }
    }
    Some(out)
}

/// `a + b sqrt(d)`: an element of the rationals extended by one square
/// root, with `d` kept alongside by the caller.
#[derive(Clone, Debug, PartialEq)]
struct Surd {
    a: BigRational,
    b: BigRational,
}

impl Surd {
    fn rational(a: BigRational) -> Self {
        Self {
            a,
            b: BigRational::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    fn sub(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            a: &self.a - &other.a,
            b: &self.b - &other.b,
        }
    }

    fn mul(
        &self,
        other: &Self,
        d: &BigRational,
    ) -> Self {
        Self {
            a: &self.a * &other.a + &self.b * &other.b * d,
            b: &self.a * &other.b + &self.b * &other.a,
        }
    }

    fn inverse(
        &self,
        d: &BigRational,
    ) -> Option<Self> {
        let norm = &self.a * &self.a - &self.b * &self.b * d;
        if norm.is_zero() {
            return None;
        }
        Some(Self {
            a: &self.a / &norm,
            b: -&self.b / &norm,
        })
    }

    fn pow(
        &self,
        exponent: u64,
        d: &BigRational,
    ) -> Self {
        (0..exponent).fold(Self::rational(BigRational::one()), |acc, _| {
            acc.mul(self, d)
        })
    }

    fn approx(
        &self,
        d: &BigRational,
    ) -> f64 {
        self.b.to_f64().unwrap_or(0.0).mul_add(
            d.to_f64().unwrap_or(0.0).sqrt(),
            self.a.to_f64().unwrap_or(0.0),
        )
    }
}

/// Solves the square system `rows * x = rhs` by Gauss-Jordan elimination.
fn solve_system(
    mut rows: Vec<Vec<Surd>>,
    mut rhs: Vec<Surd>,
    d: &BigRational,
) -> Option<Vec<Surd>> {
    let n = rows.len();
    for col in 0..n {
        let pivot = (col..n).find(|&r| {
            rows.get(r)
                .and_then(|row| row.get(col))
                .is_some_and(|v| !v.is_zero())
        })?;
        rows.swap(col, pivot);
        rhs.swap(col, pivot);
        let inverse = rows.get(col)?.get(col)?.inverse(d)?;
        let pivot_row = rows.get(col)?.clone();
        let pivot_rhs = rhs.get(col)?.clone();
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = rows.get(r)?.get(col)?.mul(&inverse, d);
            if factor.is_zero() {
                continue;
            }
            for (slot, above) in rows.get_mut(r)?.iter_mut().zip(&pivot_row) {
                *slot = slot.sub(&above.mul(&factor, d));
            }
            let slot = rhs.get_mut(r)?;
            *slot = slot.sub(&pivot_rhs.mul(&factor, d));
        }
    }
    (0..n)
        .map(|r| Some(rhs.get(r)?.mul(&rows.get(r)?.get(r)?.inverse(d)?, d)))
        .collect()
}

/// A root of the characteristic polynomial with its multiplicity.
struct Root {
    value: Surd,
    multiplicity: u32,
}

/// The roots of the characteristic polynomial with ascending coefficients
/// `c`, and the radicand of the one square root they need (zero if they
/// are all rational). `None` if some root is not in a quadratic extension
/// that is real.
fn characteristic_roots(c: &[BigRational]) -> Option<(Vec<Root>, BigRational)> {
    let mut roots = Vec::new();
    let mut radicand = BigRational::zero();
    for (factor, multiplicity) in univariate::factor(c).1 {
        match factor.as_slice() {
            | [c0, c1] => {
                roots.push(Root {
                    value: Surd::rational(BigRational::new(-c0.clone(), c1.clone())),
                    multiplicity,
                });
            },
            | [c0, c1, c2] => {
                let discriminant = c1 * c1 - BigInt::from(4) * c2 * c0;
                if !discriminant.is_positive() {
                    return None;
                }
                // discriminant = square * free, with `free` the radicand.
                let (square, free) = match super::number_theory::factor(&discriminant) {
                    | Some(primes) => {
                        primes
                            .iter()
                            .fold((BigInt::one(), BigInt::one()), |(s, f), &(p, e)| {
                                let p = BigInt::from(p);
                                (s * p.pow(e / 2), if e % 2 == 1 { f * p } else { f })
                            })
                    },
                    | None => (BigInt::one(), discriminant),
                };
                let d = BigRational::from_integer(free);
                if !radicand.is_zero() && radicand != d {
                    return None;
                }
                radicand = d;
                let two_a = BigRational::from_integer(BigInt::from(2) * c2);
                let centre = -BigRational::from_integer(c1.clone()) / &two_a;
                let spread = BigRational::from_integer(square) / two_a;
                for sign in [-1, 1] {
                    roots.push(Root {
                        value: Surd {
                            a: centre.clone(),
                            b: spread.clone() * BigRational::from_integer(BigInt::from(sign)),
                        },
                        multiplicity,
                    });
                }
            },
            | _ => return None,
        }
    }
    roots.sort_by(|x, y| {
        x.value
            .approx(&radicand)
            .total_cmp(&y.value.approx(&radicand))
    });
    Some((roots, radicand))
}

/// What a generator of a recurrence stands for.
enum Role {
    /// `a(n + offset)`.
    Unknown(i64),
    /// The index `n` itself.
    Index,
    /// `factor * base^n`.
    Power {
        factor: BigRational,
        base: BigRational,
    },
}

/// `arg = slope * n + shift` with rational constants.
fn affine(
    graph: &mut Graph,
    arg: NodeId,
    n: NodeId,
) -> Option<(BigRational, BigRational)> {
    let mut gens = Gens::default();
    let gn = gens.index(graph, n);
    let poly = from_term(graph, &mut gens, arg, Limits::default())?;
    let coefficients = poly.univariate_in(gn)?;
    match coefficients.as_slice() {
        | [shift] => Some((BigRational::zero(), shift.to_rational()?)),
        | [shift, slope] => Some((slope.to_rational()?, shift.to_rational()?)),
        | _ => None,
    }
}

fn role(
    graph: &mut Graph,
    generator: NodeId,
    function: NodeId,
    n: NodeId,
) -> Option<Role> {
    if graph.same(generator, n) {
        return Some(Role::Index);
    }
    let children = graph.children(generator).to_vec();
    match (graph.op(generator), children.as_slice()) {
        | (core::APPLY, &[f, arg]) if graph.same(f, function) => {
            let (slope, shift) = affine(graph, arg, n)?;
            (slope.is_one() && shift.is_integer())
                .then(|| shift.to_integer().to_i64())
                .flatten()
                .map(Role::Unknown)
        },
        | (core::POW, &[base, exponent]) => {
            let base = graph
                .number_of(base)?
                .to_rational()
                .filter(|b| !b.is_zero())?;
            let (slope, shift) = affine(graph, exponent, n)?;
            if !slope.is_integer() || !shift.is_integer() || slope.is_zero() {
                return None;
            }
            Some(Role::Power {
                factor: rat_pow(&base, shift.to_integer().to_i64()?)?,
                base: rat_pow(&base, slope.to_integer().to_i64()?)?,
            })
        },
        | _ => None,
    }
}

fn product_node(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn sum_node(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

/// The term `a + b * d^(1/2)`.
fn surd_node(
    graph: &mut Graph,
    value: &Surd,
    d: &BigRational,
) -> NodeId {
    let rational = graph.num(Number::rat(value.a.clone()));
    if value.b.is_zero() {
        return rational;
    }
    let radicand = graph.num(Number::rat(d.clone()));
    let half = graph.num(Number::rat(BigRational::new(
        BigInt::one(),
        BigInt::from(2),
    )));
    let root = graph.node(core::POW, &[radicand, half]);
    let scale = graph.num(Number::rat(value.b.clone()));
    let irrational = graph.node(core::MUL, &[scale, root]);
    if value.a.is_zero() {
        irrational
    } else {
        graph.node(core::ADD, &[rational, irrational])
    }
}

/// The coefficients of `sum_j c_j B^j (n + j)^e`, ascending.
fn shifted_operator(
    c: &[BigRational],
    base: &BigRational,
    e: usize,
) -> Option<Vec<BigRational>> {
    let mut unit = vec![BigRational::zero(); e + 1];
    *unit.last_mut()? = BigRational::one();
    let mut out = vec![BigRational::zero(); e + 1];
    for (j, cj) in c.iter().enumerate() {
        let j = i64::try_from(j).ok()?;
        let weight = cj * rat_pow(base, j)?;
        let shifted = poly_shift(&unit, &BigRational::from_integer(BigInt::from(j)))?;
        for (slot, value) in out.iter_mut().zip(shifted) {
            *slot += &weight * value;
        }
    }
    Some(out)
}

/// Solves the linear recurrence `equation` for the function `target`,
/// with the values `initial` of the first terms.
fn solve_recurrence(
    graph: &mut Graph,
    equation: NodeId,
    target: NodeId,
    initial: NodeId,
) -> Option<NodeId> {
    let target = best(graph, target)?;
    if graph.op(target) != core::APPLY {
        return None;
    }
    let &[function, n] = graph.children(target) else {
        return None;
    };
    graph.as_symbol(function)?;
    graph.as_symbol(n)?;
    if graph.op(initial) != core::LIST {
        return None;
    }
    let initial: Vec<BigRational> = graph
        .children(initial)
        .iter()
        .map(|&v| graph.number_of(v).and_then(Number::to_rational))
        .collect::<Option<_>>()?;

    let expr = residual(graph, equation);
    let term = best(graph, expr)?;
    let mut gens = Gens::default();
    gens.index(graph, n);
    let poly = from_term(graph, &mut gens, term, Limits::default())?;
    let mut roles = Vec::new();
    for g in 0..u32::try_from(gens.len()).ok()? {
        roles.push(role(graph, gens.node(g)?, function, n)?);
    }

    // Split into the recurrence proper (offset -> coefficient) and the
    // forcing term (base -> polynomial in n), moved to the right.
    let mut recurrence: BTreeMap<i64, BigRational> = BTreeMap::new();
    let mut forcing: BTreeMap<BigRational, Vec<BigRational>> = BTreeMap::new();
    for (mono, coeff) in poly.terms() {
        let mut c = coeff.to_rational()?;
        let (mut unknown, mut degree, mut base) = (None, 0_usize, BigRational::one());
        for &(g, e) in mono {
            match roles.get(usize::try_from(g).ok()?)? {
                | Role::Unknown(offset) => {
                    if e != 1 || unknown.is_some() {
                        return None;
                    }
                    unknown = Some(*offset);
                },
                | Role::Index => degree += usize::try_from(e).ok()?,
                | Role::Power { factor, base: b } => {
                    let e = i64::from(e);
                    c *= rat_pow(factor, e)?;
                    base *= rat_pow(b, e)?;
                },
            }
        }
        if let Some(offset) = unknown {
            if degree != 0 || !base.is_one() {
                return None;
            }
            *recurrence.entry(offset).or_insert_with(BigRational::zero) += c;
        } else {
            if degree > 64 {
                return None;
            }
            let slot = forcing.entry(base).or_default();
            if slot.len() <= degree {
                slot.resize(degree + 1, BigRational::zero());
            }
            *slot.get_mut(degree)? -= c;
        }
    }
    recurrence.retain(|_, c| !c.is_zero());
    let low = *recurrence.keys().next()?;
    let order = usize::try_from(recurrence.keys().next_back()? - low)
        .ok()
        .filter(|&o| (1..=MAX_ORDER).contains(&o))?;
    let c: Vec<BigRational> = (0..=order)
        .map(|j| {
            recurrence
                .get(&(low + i64::try_from(j).unwrap_or(0)))
                .cloned()
                .unwrap_or_else(BigRational::zero)
        })
        .collect();

    let (roots, radicand) = characteristic_roots(&c)?;
    let basis: Vec<(u32, &Root)> = roots
        .iter()
        .flat_map(|r| (0..r.multiplicity).map(move |i| (i, r)))
        .collect();
    if basis.len() != order {
        return None;
    }

    // Particular solutions, term by term, in the index n' = n + low.
    let shift = BigRational::from_integer(BigInt::from(-low));
    let mut particular: Vec<(BigRational, Vec<BigRational>)> = Vec::new();
    for (base, p) in forcing {
        let scale = rat_pow(&base, -low)?;
        let mut p: Vec<BigRational> = poly_shift(&p, &shift)?
            .into_iter()
            .map(|v| v * &scale)
            .collect();
        while p.last().is_some_and(Zero::is_zero) {
            p.pop();
        }
        if p.is_empty() {
            continue;
        }
        let s = roots
            .iter()
            .find(|r| r.value.b.is_zero() && r.value.a == base)
            .map_or(0, |r| usize::try_from(r.multiplicity).unwrap_or(0));
        let unknowns = p.len();
        let mut columns = Vec::with_capacity(unknowns);
        for u in 0..unknowns {
            columns.push(shifted_operator(&c, &base, s + u)?);
        }
        let rows: Vec<Vec<Surd>> = (0..unknowns)
            .map(|l| {
                columns
                    .iter()
                    .map(|col| {
                        Surd::rational(col.get(l).cloned().unwrap_or_else(BigRational::zero))
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<Surd> = p.iter().cloned().map(Surd::rational).collect();
        let q = solve_system(rows, rhs, &BigRational::zero())?;
        let mut coefficients = vec![BigRational::zero(); s + unknowns];
        for (u, value) in q.into_iter().enumerate() {
            *coefficients.get_mut(s + u)? = value.a;
        }
        particular.push((base, coefficients));
    }
    let particular_at = |k: u64| -> Option<BigRational> {
        let mut total = BigRational::zero();
        for (base, coefficients) in &particular {
            let mut value = BigRational::zero();
            for (e, q) in coefficients.iter().enumerate() {
                value += q * rat_pow(
                    &BigRational::from_integer(BigInt::from(k)),
                    i64::try_from(e).ok()?,
                )?;
            }
            total += value * rat_pow(base, i64::try_from(k).ok()?)?;
        }
        Some(total)
    };

    // Constants: from the initial values, or free symbols.
    let values: Option<Vec<Surd>> = if initial.is_empty() {
        None
    } else if initial.len() == order {
        let mut rows = Vec::with_capacity(order);
        let mut rhs = Vec::with_capacity(order);
        for (k, start) in initial.iter().enumerate() {
            let k = u64::try_from(k).ok()?;
            rows.push(
                basis
                    .iter()
                    .map(|&(i, root)| {
                        let index = BigRational::from_integer(BigInt::from(k));
                        Surd::rational(
                            rat_pow(&index, i64::from(i)).unwrap_or_else(BigRational::zero),
                        )
                        .mul(&root.value.pow(k, &radicand), &radicand)
                    })
                    .collect::<Vec<_>>(),
            );
            rhs.push(Surd::rational(start - particular_at(k)?));
        }
        Some(solve_system(rows, rhs, &radicand)?)
    } else {
        return None;
    };
    let constants: Vec<NodeId> = match &values {
        | Some(v) => v.iter().map(|c| surd_node(graph, c, &radicand)).collect(),
        | None => (1..=order).map(|k| graph.sym(&format!("C{k}"))).collect(),
    };

    let mut parts = Vec::new();
    for (&(i, root), constant) in basis.iter().zip(&constants) {
        if graph.number_of(*constant).is_some_and(Number::is_zero) {
            continue;
        }
        let mut factors = vec![*constant];
        if i > 0 {
            let exponent = graph.int(i64::from(i));
            factors.push(graph.node(core::POW, &[n, exponent]));
        }
        if !(root.value.b.is_zero() && root.value.a.is_one()) {
            let base = surd_node(graph, &root.value, &radicand);
            factors.push(graph.node(core::POW, &[base, n]));
        }
        parts.push(product_node(graph, &factors));
    }
    for (base, coefficients) in &particular {
        for (e, q) in coefficients.iter().enumerate() {
            if q.is_zero() {
                continue;
            }
            let mut factors = vec![graph.num(Number::rat(q.clone()))];
            if e > 0 {
                let exponent = graph.int(i64::try_from(e).ok()?);
                factors.push(graph.node(core::POW, &[n, exponent]));
            }
            if !base.is_one() {
                let base = graph.num(Number::rat(base.clone()));
                factors.push(graph.node(core::POW, &[base, n]));
            }
            parts.push(product_node(graph, &factors));
        }
    }
    Some(sum_node(graph, &parts))
}

fn rsolve(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[equation, target, initial] = args else {
        return Outcome::Pass;
    };
    solve_recurrence(cx.graph, equation, target, initial).map_or(Outcome::Pass, Outcome::Equal)
}

/// The first `count` Taylor coefficients at zero of a rational function.
fn rational_series(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    count: usize,
) -> Option<Vec<BigRational>> {
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let fraction = ratio(graph, &mut gens, term, Limits::default())?;
    let exact = |poly: &super::poly::repr::Poly| -> Option<Vec<BigRational>> {
        poly.univariate_in(gx)?
            .iter()
            .map(Number::to_rational)
            .collect()
    };
    let (a, b) = (exact(&fraction.numer)?, exact(&fraction.denom)?);
    let b0 = b.first().filter(|v| !v.is_zero())?;
    let mut series: Vec<BigRational> = Vec::with_capacity(count);
    for k in 0..count {
        let mut value = a.get(k).cloned().unwrap_or_else(BigRational::zero);
        for (j, bj) in b.iter().enumerate().skip(1).take(k) {
            value -= bj * series.get(k - j)?;
        }
        series.push(value / b0);
    }
    Some(series)
}

/// Taylor coefficients at zero by repeated differentiation, exact values
/// only; needs the `diff` request of the calculus rules.
fn derivative_series(
    cx: &mut Cx<'_>,
    term: NodeId,
    x: NodeId,
    count: usize,
) -> Option<Vec<BigRational>> {
    let diff = cx.graph.ops().lookup("diff")?;
    let zero = cx.graph.int(0);
    let (mut current, mut factorial) = (term, BigInt::one());
    let mut out = Vec::with_capacity(count);
    for k in 0..count {
        if k > 0 {
            let request = cx.graph.node(diff, &[current, x]);
            current = cx.simplify(request);
            factorial *= k;
        }
        let at_zero = cx.graph.substitute(current, x, zero);
        let value = cx.simplify(at_zero);
        let r = cx.graph.number_of(value)?.to_rational()?;
        out.push(r / BigRational::from_integer(factorial.clone()));
    }
    Some(out)
}

fn gf_coefficients(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[f, x, count] = args else {
        return Outcome::Pass;
    };
    let Some(count) = int_literal(cx.graph, count)
        .and_then(|c| c.to_u64())
        .filter(|&c| c <= MAX_SERIES)
    else {
        return Outcome::Pass;
    };
    if cx.graph.as_symbol(x).is_none() {
        return Outcome::Pass;
    }
    let Some(term) = best(cx.graph, f) else {
        return Outcome::Pass;
    };
    let len = usize::try_from(count).unwrap_or(0);
    let series = rational_series(cx.graph, term, x, len).or_else(|| {
        if count > MAX_TAYLOR {
            return None;
        }
        derivative_series(cx, term, x, len)
    });
    let Some(series) = series else {
        return Outcome::Pass;
    };
    let nodes: Vec<NodeId> = series
        .into_iter()
        .map(|c| cx.graph.num(Number::rat(c)))
        .collect();
    Outcome::Equal(cx.graph.node(core::LIST, &nodes))
}

/// `|A_1 ∪ ... ∪ A_n|` from the sizes of the intersections, level by
/// level: `sum_1 - sum_2 + sum_3 - ...`.
fn inclusion_exclusion(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[levels] = args else {
        return Outcome::Pass;
    };
    let graph = &mut *cx.graph;
    if graph.op(levels) != core::LIST {
        return Outcome::Pass;
    }
    let mut terms = Vec::new();
    for (level, &sizes) in graph.children(levels).to_vec().iter().enumerate() {
        if graph.op(sizes) != core::LIST {
            return Outcome::Pass;
        }
        let sign = graph.int(if level % 2 == 0 { 1 } else { -1 });
        for &size in &graph.children(sizes).to_vec() {
            terms.push(graph.node(core::MUL, &[sign, size]));
        }
    }
    Outcome::Equal(sum_node(graph, &terms))
}

/// The least `p` dividing the length such that the list repeats every `p`
/// items; the length itself when there is no shorter period.
fn period(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[sequence] = args else {
        return Outcome::Pass;
    };
    let graph = &mut *cx.graph;
    if graph.op(sequence) != core::LIST {
        return Outcome::Pass;
    }
    let items = graph.children(sequence).to_vec();
    let n = items.len();
    let found = (1..=n).find(|&p| {
        n.is_multiple_of(p)
            && items
                .iter()
                .zip(items.iter().skip(p))
                .all(|(&a, &b)| graph.same(a, b))
    });
    match found.and_then(|p| i64::try_from(p).ok()) {
        | Some(p) => Outcome::Equal(graph.int(p)),
        | None => Outcome::Pass,
    }
}

/// `expand_binomial(a, b, n)`: for a literal `n`, every term of the
/// expansion has dimension `a^n`; for a symbolic `n` the expansion is only
/// dimensionally coherent when `a` and `b` are dimensionless, but the
/// identity itself has the dimension of `a`.
#[allow(clippy::unnecessary_wraps)]
fn binomial_expansion_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let base = args.first().cloned().unwrap_or_else(DimTerm::dimensionless);
    let exponent = ctx
        .graph
        .children(ctx.node)
        .get(2)
        .and_then(|&node| ctx.graph.number_of(node))
        .and_then(Number::to_rational)
        .filter(|exponent| !exponent.is_negative() && exponent.denom().is_one());
    Ok(exponent.map_or_else(|| base.clone(), |exponent| base.pow(&exponent)))
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));
    let set = "combinatorics";

    exact(i, set, unary("factorial").eval(factorial_eval), |a| {
        match a {
            | [n] => Some(Value::Int(factorial(nat(n, MAX_FACTORIAL)?))),
            | _ => None,
        }
    })?;
    exact(i, set, unary("double_factorial"), |a| {
        match a {
            | [n] => Some(Value::Int(double_factorial(nat(n, 2 * MAX_FACTORIAL)?))),
            | _ => None,
        }
    })?;
    exact(i, set, binary("binomial").eval(binomial_eval), |a| {
        match a {
            | [n, k] => binomial(n, k).map(Value::Int),
            | _ => None,
        }
    })?;
    exact_lists(i, set, unary("multinomial"), multinomial)?;
    install_ext(i)?;
    exact(
        i,
        set,
        binary("permutations").eval(permutations_eval),
        |a| {
            match a {
                | [n, k] if !n.is_negative() && !k.is_negative() => {
                    if k > n {
                        Some(Value::Int(BigInt::zero()))
                    } else {
                        product_run(n, k, -1).map(Value::Int)
                    }
                },
                | _ => None,
            }
        },
    )?;
    exact(i, set, binary("rising"), |a| {
        match a {
            | [x, k] => product_run(x, k, 1).map(Value::Int),
            | _ => None,
        }
    })?;
    exact(i, set, binary("falling"), |a| {
        match a {
            | [x, k] => product_run(x, k, -1).map(Value::Int),
            | _ => None,
        }
    })?;
    exact(i, set, unary("catalan"), |a| {
        match a {
            | [n] => {
                let n = nat(n, MAX_FACTORIAL / 2)?;
                let c = choose(&BigInt::from(2 * n), n)?;
                Some(Value::Int(c / (n + 1)))
            },
            | _ => None,
        }
    })?;
    exact(i, set, unary("fibonacci"), |a| {
        match a {
            | [n] => Some(Value::Int(fibonacci_pair(nat(n, 100_000)?).0)),
            | _ => None,
        }
    })?;
    exact(i, set, unary("lucas"), |a| {
        match a {
            | [n] => {
                let (f, next) = fibonacci_pair(nat(n, 100_000)?);
                Some(Value::Int(next * 2 - f))
            },
            | _ => None,
        }
    })?;
    exact(i, set, unary("bell"), |a| {
        match a {
            | [n] => Some(Value::Int(bell(nat(n, 1000)?))),
            | _ => None,
        }
    })?;
    exact(i, set, binary("stirling1"), |a| {
        match a {
            | [n, k] => {
                let (n, k) = (nat(n, 1000)?, k.to_u64()?);
                let row = stirling_row(n, |m, _| m - 1);
                Some(Value::Int(
                    row.get(usize::try_from(k).ok()?)
                        .cloned()
                        .unwrap_or_default(),
                ))
            },
            | _ => None,
        }
    })?;
    exact(i, set, binary("stirling2"), |a| {
        match a {
            | [n, k] => {
                let (n, k) = (nat(n, 1000)?, k.to_u64()?);
                let row = stirling_row(n, |_, j| j);
                Some(Value::Int(
                    row.get(usize::try_from(k).ok()?)
                        .cloned()
                        .unwrap_or_default(),
                ))
            },
            | _ => None,
        }
    })?;
    exact(i, set, unary("partitions"), |a| {
        match a {
            | [n] => Some(Value::Int(partitions(usize::try_from(nat(n, 5000)?).ok()?))),
            | _ => None,
        }
    })?;
    exact(i, set, unary("derangements"), |a| {
        match a {
            | [n] => Some(Value::Int(derangements(nat(n, MAX_FACTORIAL)?))),
            | _ => None,
        }
    })?;
    exact(
        i,
        set,
        unary("harmonic").eval(|a| {
            match a {
                // H(x) = ψ(x + 1) + γ, the analytic continuation.
                | [x] => {
                    crate::kernels::special::digamma_numerical(x + 1.0)
                        + std::f64::consts::EULER_GAMMA
                },
                | _ => f64::NAN,
            }
        }),
        |a| {
            match a {
                | [n] => Some(Value::Rat(harmonic(nat(n, 1000)?))),
                | _ => None,
            }
        },
    )?;

    procedure(
        i,
        set,
        OpDescriptor::new("rsolve", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(100),
        Tier::Reduce,
        rsolve,
    )?;
    procedure(
        i,
        set,
        OpDescriptor::new("gf_coeffs", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(100),
        Tier::Reduce,
        gf_coefficients,
    )?;
    procedure(
        i,
        set,
        unary("inclusion_exclusion"),
        Tier::Normalize,
        inclusion_exclusion,
    )?;
    procedure(i, set, unary("period"), Tier::Normalize, period)?;

    // Combinatorial numbers count configurations: they are pure numbers,
    // their multi-part results are lists of integers, and the expansion
    // identity is homogeneous in its base.
    attach_dimensions(
        i,
        &[
            ("factorial", DimRule::Dimensionless),
            ("double_factorial", DimRule::Dimensionless),
            ("binomial", DimRule::Dimensionless),
            ("multinomial", DimRule::Dimensionless),
            ("lah", DimRule::Dimensionless),
            ("eulerian", DimRule::Dimensionless),
            ("narayana", DimRule::Dimensionless),
            ("motzkin", DimRule::Dimensionless),
            ("schroeder", DimRule::Dimensionless),
            ("fubini", DimRule::Dimensionless),
            ("central_binomial", DimRule::Dimensionless),
            ("ballot", DimRule::Dimensionless),
            ("involutions", DimRule::Dimensionless),
            ("rencontres", DimRule::Dimensionless),
            ("partitions_k", DimRule::Dimensionless),
            ("partitions_max", DimRule::Dimensionless),
            ("partitions_distinct", DimRule::Dimensionless),
            ("partition_list", DimRule::List),
            ("compositions", DimRule::Dimensionless),
            ("weak_compositions", DimRule::Dimensionless),
            ("necklaces", DimRule::Dimensionless),
            ("bracelets", DimRule::Dimensionless),
            ("bernoulli_number", DimRule::Dimensionless),
            ("zigzag", DimRule::Dimensionless),
            ("euler_number", DimRule::Dimensionless),
            ("permutation_rank", DimRule::Dimensionless),
            ("permutation_unrank", DimRule::Dimensionless),
            ("inversions", DimRule::Dimensionless),
            ("permutation_sign", DimRule::Dimensionless),
            ("cycle_type", DimRule::List),
            ("permutation_order", DimRule::Dimensionless),
            ("hook_lengths", DimRule::List),
            ("syt_count", DimRule::Dimensionless),
            ("permutations", DimRule::Dimensionless),
            ("rising", DimRule::Dimensionless),
            ("falling", DimRule::Dimensionless),
            ("catalan", DimRule::Dimensionless),
            ("fibonacci", DimRule::Dimensionless),
            ("lucas", DimRule::Dimensionless),
            ("bell", DimRule::Dimensionless),
            ("stirling1", DimRule::Dimensionless),
            ("stirling2", DimRule::Dimensionless),
            ("partitions", DimRule::Dimensionless),
            ("derangements", DimRule::Dimensionless),
            ("harmonic", DimRule::Dimensionless),
            ("rsolve", DimRule::NonMetric),
            ("gf_coeffs", DimRule::List),
            ("inclusion_exclusion", DimRule::Dimensionless),
            ("period", DimRule::Dimensionless),
        ],
    )?;

    // All rules hold for the continuous extensions used as numeric
    // semantics, hence for every real argument satisfying the guards.
    i.rewrites(
        Tier::Normalize,
        &[
            "combinatorics/binomial-0: binomial(?n, 0) => 1",
            "combinatorics/binomial-1: binomial(?n, 1) => ?n if nonnegative(?n)",
            "combinatorics/binomial-self: binomial(?n, ?n) => 1 if nonnegative(?n)",
            "combinatorics/permutations-0: permutations(?n, 0) => 1 if nonnegative(?n)",
            "combinatorics/permutations-1: permutations(?n, 1) => ?n if nonnegative(?n)",
            "combinatorics/permutations-self: permutations(?n, ?n) => factorial(?n) \
             if nonnegative(?n)",
            "combinatorics/binomial-permutations: binomial(?n, ?k) * factorial(?k) \
             => permutations(?n, ?k) if nonnegative(?n), nonnegative(?k)",
            "combinatorics/factorial-ratio: factorial(?n + 1) / factorial(?n) => ?n + 1 \
             if nonnegative(?n)",
            "combinatorics/factorial-ratio-inverse: factorial(?n) / factorial(?n + 1) => 1/(?n + 1) \
             if nonnegative(?n)",
            "combinatorics/binomial-ratio-k: binomial(?n, ?k + 1) / binomial(?n, ?k) => (?n - ?k)/(?k + 1) \
             if integer(?k), nonnegative(?k), nonnegative(?n)",
            "combinatorics/binomial-ratio-n: binomial(?n + 1, ?k) / binomial(?n, ?k) => (?n + 1)/(?n + 1 - ?k) \
             if integer(?k), nonnegative(?k), nonnegative(?n)",
        ],
    )?;
    i.define(&["expand_binomial(a, b, n) := sum(binomial(n, j)*a^(n - j)*b^j, j, 0, n)"])?;
    let expand = i.graph().ops().lookup("expand_binomial").ok_or_else(|| {
        RuleError::Invalid {
            rule: "dimensions/expand_binomial".to_owned(),
            reason: "operator is not registered",
        }
    })?;
    i.dim(expand, DimRule::Custom(binomial_expansion_dimension));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[combinatorics()], src)
    }

    fn value(src: &str) -> BigInt {
        s(src)
            .parse()
            .unwrap_or_else(|_| panic!("`{src}` did not reduce to an integer"))
    }

    #[test]
    fn factorials() {
        assert_eq!(s("factorial(0)"), "1");
        assert_eq!(s("factorial(1)"), "1");
        assert_eq!(s("factorial(5)"), "120");
        assert_eq!(s("factorial(20)"), "2432902008176640000");
        assert_eq!(s("factorial(30)"), "265252859812191058636308480000000");
        assert_eq!(s("factorial(-1)"), "factorial(-1)");
        assert_eq!(s("factorial(1000000)"), "factorial(1000000)");
        assert_eq!(s("double_factorial(0)"), "1");
        assert_eq!(s("double_factorial(7)"), "105");
        assert_eq!(s("double_factorial(8)"), "384");
        assert_eq!(s("double_factorial(10)"), "3840");
        assert_eq!(s("double_factorial(-3)"), "double_factorial(-3)");
    }

    #[test]
    fn binomials() {
        assert_eq!(s("binomial(5, 2)"), "10");
        assert_eq!(s("binomial(5, 0)"), "1");
        assert_eq!(s("binomial(5, 5)"), "1");
        assert_eq!(s("binomial(0, 0)"), "1");
        assert_eq!(s("binomial(5, 7)"), "0");
        assert_eq!(s("binomial(5, -1)"), "0");
        assert_eq!(s("binomial(-3, 2)"), "6");
        assert_eq!(s("binomial(-2, 3)"), "-4");
        assert_eq!(s("binomial(100, 50)"), "100891344545564193334812497256");
        assert_eq!(s("binomial(100, 3)"), "161700");
        assert_eq!(
            s("binomial(10^20, 2)"),
            "4999999999999999999950000000000000000000"
        );
        assert_eq!(
            s("binomial(10^30, 10^9)"),
            "binomial(1000000000000000000000000000000, 1000000000)"
        );
        assert_eq!(s("multinomial(list(2, 3, 4))"), "1260");
        assert_eq!(s("multinomial(list())"), "1");
        assert_eq!(s("multinomial(list(5))"), "1");
        assert_eq!(s("multinomial(list(1, 1, 1, 1))"), "24");
        assert_eq!(s("multinomial(list(2, -1))"), "multinomial(list(2, -1))");
    }

    #[test]
    fn falling_and_rising() {
        assert_eq!(s("permutations(5, 2)"), "20");
        assert_eq!(s("permutations(5, 0)"), "1");
        assert_eq!(s("permutations(5, 5)"), "120");
        assert_eq!(s("permutations(3, 5)"), "0");
        assert_eq!(s("permutations(0, 0)"), "1");
        assert_eq!(s("rising(3, 4)"), "360");
        assert_eq!(s("rising(-2, 4)"), "0");
        assert_eq!(s("rising(7, 0)"), "1");
        assert_eq!(s("falling(5, 3)"), "60");
        assert_eq!(s("falling(-1, 3)"), "-6");
        assert_eq!(s("falling(3, 5)"), "0");
        assert_eq!(s("falling(5, -1)"), "falling(5, -1)");
    }

    #[test]
    fn counting_sequences() {
        let expect = |op: &str, values: &[&str]| {
            for (n, v) in values.iter().enumerate() {
                assert_eq!(s(&format!("{op}({n})")), *v, "{op}({n})");
            }
        };
        expect("catalan", &["1", "1", "2", "5", "14", "42", "132", "429"]);
        expect("fibonacci", &["0", "1", "1", "2", "3", "5", "8", "13"]);
        expect("lucas", &["2", "1", "3", "4", "7", "11", "18", "29"]);
        expect(
            "bell",
            &[
                "1", "1", "2", "5", "15", "52", "203", "877", "4140", "21147", "115975",
            ],
        );
        expect(
            "partitions",
            &["1", "1", "2", "3", "5", "7", "11", "15", "22", "30", "42"],
        );
        expect(
            "derangements",
            &["1", "0", "1", "2", "9", "44", "265", "1854"],
        );
        assert_eq!(s("catalan(30)"), "3814986502092304");
        assert_eq!(s("fibonacci(100)"), "354224848179261915075");
        assert_eq!(
            s("fibonacci(200)"),
            "280571172992510140037611932413038677189525"
        );
        assert_eq!(s("lucas(100)"), "792070839848372253127");
        assert_eq!(s("bell(20)"), "51724158235372");
        assert_eq!(s("partitions(50)"), "204226");
        assert_eq!(s("partitions(100)"), "190569292");
        assert_eq!(s("derangements(10)"), "1334961");
        assert_eq!(s("catalan(-1)"), "catalan(-1)");
        assert_eq!(s("fibonacci(-1)"), "fibonacci(-1)");
        assert_eq!(s("bell(-1)"), "bell(-1)");
        assert_eq!(s("partitions(-1)"), "partitions(-1)");
    }

    #[test]
    fn stirling_numbers() {
        assert_eq!(s("stirling1(0, 0)"), "1");
        assert_eq!(s("stirling1(3, 0)"), "0");
        assert_eq!(s("stirling1(4, 2)"), "11");
        assert_eq!(s("stirling1(5, 2)"), "50");
        assert_eq!(s("stirling1(5, 5)"), "1");
        assert_eq!(s("stirling1(10, 3)"), "1172700");
        assert_eq!(s("stirling1(3, 5)"), "0");
        assert_eq!(s("stirling2(0, 0)"), "1");
        assert_eq!(s("stirling2(4, 2)"), "7");
        assert_eq!(s("stirling2(5, 2)"), "15");
        assert_eq!(s("stirling2(6, 3)"), "90");
        assert_eq!(s("stirling2(10, 3)"), "9330");
        assert_eq!(s("stirling2(3, 5)"), "0");
        assert_eq!(s("stirling2(5, -1)"), "stirling2(5, -1)");
    }

    #[test]
    fn harmonic_numbers() {
        assert_eq!(s("harmonic(0)"), "0");
        assert_eq!(s("harmonic(1)"), "1");
        assert_eq!(s("harmonic(2)"), "3/2");
        assert_eq!(s("harmonic(5)"), "137/60");
        assert_eq!(s("harmonic(10)"), "7381/2520");
    }

    #[test]
    fn symbolic_arguments_stay() {
        let sets = [combinatorics()];
        for (src, expected) in [
            ("factorial(n)", "factorial(n)"),
            ("binomial(n, 2)", "binomial(n, 2)"),
            ("fibonacci(n)", "fibonacci(n)"),
            ("stirling2(3, k)", "stirling2(3, k)"),
            ("multinomial(list(a, 2))", "multinomial(list(a, 2))"),
            ("harmonic(n)", "harmonic(n)"),
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(reduced);
            assert_eq!(text, expected);
        }
    }

    #[test]
    fn rewrites_fire() {
        let sets = [combinatorics()];
        let run =
            |src: &str, assume: &[(&str, crate::graph::Facts)]| reduce_with(&sets, src, assume).0;
        let nonneg = [("n", crate::graph::Facts::NONNEGATIVE)];
        assert_eq!(run("binomial(n, 0)", &[]), "1");
        assert_eq!(run("binomial(n, 1)", &nonneg), "n");
        assert_eq!(run("binomial(n, 1)", &[]), "binomial(n, 1)");
        assert_eq!(run("binomial(n, n)", &nonneg), "1");
        assert_eq!(run("binomial(n, n)", &[]), "binomial(n, n)");
        assert_eq!(run("permutations(n, 0)", &nonneg), "1");
        assert_eq!(run("permutations(n, 1)", &nonneg), "n");
        assert_eq!(run("permutations(n, n)", &nonneg), "factorial(n)");
        assert_eq!(run("factorial(n + 1) / factorial(n)", &nonneg), "n + 1");
        assert_eq!(
            run("factorial(n + 1) / factorial(n)", &[]),
            "factorial(n + 1)/factorial(n)"
        );
        let both = [
            ("n", crate::graph::Facts::NONNEGATIVE),
            ("k", crate::graph::Facts::NONNEGATIVE),
        ];
        assert_eq!(
            run("binomial(n, k) * factorial(k)", &both),
            "permutations(n, k)"
        );
    }

    #[test]
    fn float_semantics() {
        let sets = [combinatorics()];
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs().max(1.0);
        assert!(close(
            eval(&sets, "factorial(x)", &[("x", 10.0)]),
            3_628_800.0
        ));
        assert!(close(
            eval(&sets, "factorial(x)", &[("x", 0.5)]),
            0.886_226_925_452_758
        ));
        assert!(eval(&sets, "factorial(x)", &[("x", -2.0)]).is_nan());
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 10.0), ("y", 3.0)]),
            120.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 3.0), ("y", 5.0)]),
            0.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", -3.0), ("y", 2.0)]),
            6.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 4.0), ("y", -1.0)]),
            0.0
        ));
        assert!(close(
            eval(&sets, "permutations(x, y)", &[("x", 5.0), ("y", 2.0)]),
            20.0
        ));
    }

    #[test]
    fn fibonacci_recurrence() {
        for n in 0..60 {
            let sum = value(&format!("fibonacci({n})")) + value(&format!("fibonacci({})", n + 1));
            assert_eq!(sum, value(&format!("fibonacci({})", n + 2)), "n = {n}");
        }
        // d'Ocagne / Cassini: F(n-1) F(n+1) - F(n)^2 = (-1)^n.
        for n in 1..40_i32 {
            let lhs = value(&format!("fibonacci({})", n - 1))
                * value(&format!("fibonacci({})", n + 1))
                - value(&format!("fibonacci({n})")).pow(2);
            assert_eq!(
                lhs,
                BigInt::from(if n % 2 == 0 { 1 } else { -1 }),
                "n = {n}"
            );
        }
        // Lucas numbers satisfy the same recurrence and L(n) = F(n-1) + F(n+1).
        for n in 1..40 {
            assert_eq!(
                value(&format!("lucas({n})")),
                value(&format!("fibonacci({})", n - 1)) + value(&format!("fibonacci({})", n + 1)),
                "n = {n}"
            );
        }
    }

    #[test]
    fn binomial_recurrence() {
        // Pascal's rule.
        for n in 1..14 {
            for k in 0..=n {
                let lhs = value(&format!("binomial({n}, {k})"));
                let rhs = value(&format!("binomial({}, {})", n - 1, k - 1))
                    + value(&format!("binomial({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "C({n}, {k})");
                assert_eq!(lhs, value(&format!("binomial({n}, {})", n - k)), "symmetry");
            }
        }
        // C(n, k) = n! / (k! (n-k)!) and row sums 2^n.
        for n in 0..12 {
            let mut total = BigInt::zero();
            for k in 0..=n {
                let c = value(&format!("binomial({n}, {k})"));
                assert_eq!(
                    c * value(&format!("factorial({k})")) * value(&format!("factorial({})", n - k)),
                    value(&format!("factorial({n})"))
                );
                total += value(&format!("binomial({n}, {k})"));
            }
            assert_eq!(total, BigInt::from(2).pow(n), "row {n}");
        }
    }

    #[test]
    fn stirling_and_bell_recurrences() {
        for n in 1..10 {
            for k in 1..=n {
                let lhs = value(&format!("stirling2({n}, {k})"));
                let rhs = value(&format!("stirling2({}, {})", n - 1, k - 1))
                    + BigInt::from(k) * value(&format!("stirling2({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "S({n}, {k})");
                let lhs = value(&format!("stirling1({n}, {k})"));
                let rhs = value(&format!("stirling1({}, {})", n - 1, k - 1))
                    + BigInt::from(n - 1) * value(&format!("stirling1({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "c({n}, {k})");
            }
            // Bell numbers are row sums of the Stirling numbers of the second kind.
            let total: BigInt = (0..=n)
                .map(|k| value(&format!("stirling2({n}, {k})")))
                .sum();
            assert_eq!(total, value(&format!("bell({n})")));
        }
    }

    #[test]
    fn partition_recurrence() {
        // n p(n) = sum_{k=1..n} sigma(k) p(n-k), with sigma from number theory.
        let sets = [combinatorics(), crate::rules::number_theory()];
        let v = |src: &str| -> BigInt {
            simplify(&sets, src)
                .parse()
                .unwrap_or_else(|_| panic!("`{src}` did not reduce to an integer"))
        };
        for n in 1..40 {
            let rhs: BigInt = (1..=n)
                .map(|k| v(&format!("divisor_sum({k})")) * v(&format!("partitions({})", n - k)))
                .sum();
            assert_eq!(
                BigInt::from(n) * v(&format!("partitions({n})")),
                rhs,
                "n = {n}"
            );
        }
        // Derangements: D(n) = n D(n-1) + (-1)^n.
        for n in 1..30_i32 {
            let sign = if n % 2 == 0 { 1 } else { -1 };
            assert_eq!(
                value(&format!("derangements({n})")),
                BigInt::from(n) * value(&format!("derangements({})", n - 1)) + sign
            );
        }
    }

    /// The closed form of a recurrence, as a function of `n`.
    fn closed_form(src: &str) -> impl Fn(f64) -> f64 + use<> {
        let text = s(src);
        move |n| eval(&[combinatorics()], &text, &[("n", n)])
    }

    fn close(
        a: f64,
        b: f64,
    ) {
        assert!((a - b).abs() <= 1e-9 * (1.0 + b.abs()), "{a} != {b}");
    }

    #[test]
    fn linear_recurrences_with_distinct_rational_roots() {
        assert_eq!(
            s("rsolve(a(n) = 5*a(n-1) - 6*a(n-2), a(n), list(1, 4))"),
            "2*3^n - 2^n"
        );
        assert_eq!(s("rsolve(a(n+1) = 2*a(n), a(n), list(3))"), "3*2^n");
        assert_eq!(
            s("rsolve(a(n+2) = 3*a(n+1) - 2*a(n), a(n), list(0, 1))"),
            "2^n - 1"
        );
        // A negative root and a fractional root.
        let a = closed_form("rsolve(a(n+2) = -a(n+1) + 6*a(n), a(n), list(1, 0))");
        let b = closed_form("rsolve(a(n+2) = a(n+1)/2 + a(n)/2, a(n), list(1, 2))");
        let (mut p, mut q) = ((1.0, 0.0), (1.0, 2.0));
        for n in 0..12 {
            close(a(f64::from(n)), p.0);
            close(b(f64::from(n)), q.0);
            p = (p.1, -p.1 + 6.0 * p.0);
            q = (q.1, q.1 / 2.0 + q.0 / 2.0);
        }
        // Without initial values the constants stay free.
        assert_eq!(
            s("rsolve(a(n+2) = 3*a(n+1) - 2*a(n), a(n), list())"),
            "C2*2^n + C1"
        );
    }

    #[test]
    fn linear_recurrences_with_repeated_roots() {
        assert_eq!(
            s("rsolve(a(n+2) = 4*a(n+1) - 4*a(n), a(n), list(1, 4))"),
            "n*2^n + 2^n"
        );
        // (r - 1)^3: quadratic growth.
        let a = closed_form("rsolve(a(n+3) = 3*a(n+2) - 3*a(n+1) + a(n), a(n), list(1, 2, 5))");
        for n in 0..10 {
            let n = f64::from(n);
            close(a(n), n * n + 1.0);
        }
        // (r - 2)^2 (r + 1)
        let b = closed_form("rsolve(a(n+3) = 3*a(n+2) - 4*a(n), a(n), list(1, 0, 2))");
        let mut window = [1.0, 0.0, 2.0];
        for n in 0..12 {
            close(b(f64::from(n)), window[0]);
            let next = 3.0 * window[2] - 4.0 * window[0];
            window = [window[1], window[2], next];
        }
    }

    #[test]
    fn linear_recurrences_with_irrational_roots() {
        // Binet's formula.
        let f = closed_form("rsolve(a(n+2) = a(n+1) + a(n), a(n), list(0, 1))");
        for n in 0..30 {
            close(f(f64::from(n)), fib(n));
            close(
                f(f64::from(n)),
                s(&format!("fibonacci({n})")).parse().unwrap_or(f64::NAN),
            );
        }
        // Lucas numbers satisfy the same recurrence.
        let l = closed_form("rsolve(a(n+2) = a(n+1) + a(n), a(n), list(2, 1))");
        for n in 0..25 {
            close(
                l(f64::from(n)),
                s(&format!("lucas({n})")).parse().unwrap_or(f64::NAN),
            );
        }
        // sqrt(2) roots, and a radicand with a square factor: r^2 - 2r - 7.
        let g = closed_form("rsolve(a(n+2) = 2*a(n+1) + 7*a(n), a(n), list(1, 3))");
        let mut w = (1.0, 3.0);
        for n in 0..10 {
            close(g(f64::from(n)), w.0);
            w = (w.1, 2.0 * w.1 + 7.0 * w.0);
        }
        // Mixed with a rational root: (r - 2)(r^2 - r - 1).
        let h = closed_form("rsolve(a(n+3) = 3*a(n+2) - a(n+1) - 2*a(n), a(n), list(1, 1, 4))");
        let mut w = [1.0, 1.0, 4.0];
        for n in 0..12 {
            close(h(f64::from(n)), w[0]);
            let next = 3.0 * w[2] - w[1] - 2.0 * w[0];
            w = [w[1], w[2], next];
        }
    }

    fn fib(n: i32) -> f64 {
        (0..n).fold((0.0, 1.0), |(a, b), _| (b, a + b)).0
    }

    #[test]
    fn inhomogeneous_recurrences() {
        assert_eq!(
            s("rsolve(a(n+1) = a(n) + n, a(n), list(0))"),
            "1/2*n^2 - 1/2*n"
        );
        assert_eq!(s("rsolve(a(n+1) = 2*a(n) + 3^n, a(n), list(1))"), "3^n");
        // Resonance: the forcing base is a characteristic root.
        assert_eq!(
            s("rsolve(a(n+1) = 2*a(n) + 2^n, a(n), list(1))"),
            "1/2*n*2^n + 2^n"
        );
        // Polynomial times exponential, with shifted arguments.
        let a = closed_form("rsolve(a(n+2) = 3*a(n+1) - 2*a(n) + 2^(n+1) + n^2, a(n), list(1, 1))");
        let mut w = (1.0, 1.0);
        for n in 0..14 {
            close(a(f64::from(n)), w.0);
            let k = f64::from(n);
            w = (w.1, 3.0 * w.1 - 2.0 * w.0 + 2.0_f64.powf(k + 1.0) + k * k);
        }
        // Constant forcing on a root at one (double resonance) and a
        // combination of two bases.
        let b = closed_form("rsolve(a(n+2) = 2*a(n+1) - a(n) + 6, a(n), list(0, 0))");
        for n in 0..10 {
            let n = f64::from(n);
            close(b(n), 3.0 * n * n - 3.0 * n);
        }
        let c = closed_form("rsolve(a(n+1) = 2*a(n) + 3^n + 5*7^n, a(n), list(2))");
        let mut x = 2.0;
        for n in 0..10 {
            close(c(f64::from(n)), x);
            x = 2.0 * x + 3.0_f64.powi(n) + 5.0 * 7.0_f64.powi(n);
        }
    }

    #[test]
    fn recurrence_solutions_satisfy_the_recurrence() {
        // Substituting the closed form back into the equation, at real
        // arguments (not just at the integers the initial values pin).
        for (equation, coefficients) in [
            ("a(n+2) = a(n+1) + a(n)", (1.0, 1.0)),
            ("a(n+2) = 2*a(n+1) + 7*a(n)", (2.0, 7.0)),
            ("a(n+2) = 4*a(n+1) - 4*a(n)", (4.0, -4.0)),
            ("a(n+2) = 3*a(n+1) - 2*a(n)", (3.0, -2.0)),
        ] {
            let a = closed_form(&format!("rsolve({equation}, a(n), list(2, 3))"));
            for n in [0.0, 1.0, 2.0, 5.0] {
                close(
                    a(n + 2.0),
                    coefficients.0 * a(n + 1.0) + coefficients.1 * a(n),
                );
            }
        }
        // With forcing: a(n+1) - 2 a(n) = 3^n + n.
        let a = closed_form("rsolve(a(n+1) = 2*a(n) + 3^n + n, a(n), list(5))");
        for n in [0.0, 1.0, 2.0, 4.0, 7.0] {
            close(a(n + 1.0) - 2.0 * a(n), 3.0_f64.powf(n) + n);
        }
        // The numeric kernel unrolls the same recurrence.
        let a = closed_form("rsolve(a(n) = 3*a(n-1) - a(n-2), a(n), list(1, 5))");
        for n in 0..10_u32 {
            let want = crate::kernels::combinatorics::solve_recurrence_numerical(
                &[3.0, -1.0],
                &[1.0, 5.0],
                n as usize,
            )
            .unwrap_or(f64::NAN);
            close(a(f64::from(n)), want);
        }
    }

    #[test]
    fn unsolvable_recurrences_stay_requests() {
        let sets = [combinatorics()];
        for src in [
            // nonlinear
            "rsolve(a(n+1) = a(n) + a(n)^2, a(n), list(1))",
            // non-constant coefficients
            "rsolve(a(n+1) = n*a(n), a(n), list(1))",
            // wrong number of initial values
            "rsolve(a(n+2) = a(n+1) + a(n), a(n), list(1))",
            // symbolic initial values
            "rsolve(a(n+1) = 2*a(n), a(n), list(c))",
            // complex characteristic roots
            "rsolve(a(n+2) = -a(n), a(n), list(1, 0))",
            // an equation that is not a recurrence
            "rsolve(a(n) = 3, a(n), list())",
            // foreign function or symbol in the forcing term
            "rsolve(a(n+1) = a(n) + b(n), a(n), list(1))",
            "rsolve(a(n+1) = a(n) + x, a(n), list(1))",
            // not a function application
            "rsolve(a(n+1) = a(n), n, list(1))",
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(!reduced, "{src} => {text}");
        }
    }

    #[test]
    fn series_of_generating_functions() {
        let sets = [combinatorics(), crate::rules::calculus()];
        let g = |src: &str| simplify(&sets, src);
        assert_eq!(
            g("gf_coeffs(1/(1 - x - x^2), x, 8)"),
            "list(1, 1, 2, 3, 5, 8, 13, 21)"
        );
        assert_eq!(g("gf_coeffs(1/(1 - x), x, 4)"), "list(1, 1, 1, 1)");
        assert_eq!(
            g("gf_coeffs((1 + x)^5, x, 8)"),
            "list(1, 5, 10, 10, 5, 1, 0, 0)"
        );
        assert_eq!(g("gf_coeffs(x/(1 - x)^2, x, 5)"), "list(0, 1, 2, 3, 4)");
        assert_eq!(
            g("gf_coeffs(1/((1 - x)*(1 - 2*x)), x, 5)"),
            "list(1, 3, 7, 15, 31)"
        );
        assert_eq!(g("gf_coeffs(1/(1 - x - x^2), x, 0)"), "list()");
        // Not rational: found by differentiation.
        assert_eq!(g("gf_coeffs(exp(x), x, 5)"), "list(1, 1, 1/2, 1/6, 1/24)");
        assert_eq!(
            g("gf_coeffs(sin(x), x, 6)"),
            "list(0, 1, 0, -1/6, 0, 1/120)"
        );
        assert_eq!(
            g("gf_coeffs((1 - 4*x)^(-1/2), x, 5)"),
            "list(1, 2, 6, 20, 70)"
        );
        // sqrt(1 - 4x) = 1 - 2 * sum C(n-1) x^n with the Catalan numbers.
        assert_eq!(
            g("gf_coeffs((1 - 4*x)^(1/2), x, 5)"),
            "list(1, -2, -2, -4, -10)"
        );
        let expected: Vec<String> = (0..4).map(|n| s(&format!("2*catalan({n})"))).collect();
        assert_eq!(expected, ["2", "2", "4", "10"]);
        // Symbolic coefficients or unsuitable arguments stay requests.
        for src in [
            "gf_coeffs(1/(1 - a*x), x, 3)",
            "gf_coeffs(ln(x), x, 3)",
            "gf_coeffs(1/x, x, 3)",
            "gf_coeffs(1/(1 - x), 2, 3)",
            "gf_coeffs(1/(1 - x), x, n)",
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(!reduced, "{src} => {text}");
        }
        // Calculus is a dependency, so transcendental functions expand too.
        let (text, reduced) = reduce_with(&[combinatorics()], "gf_coeffs(exp(x), x, 3)", &[]);
        assert!(reduced && text == "list(1, 1, 1/2)", "{text}");
    }

    #[test]
    fn inclusion_exclusion_counts_unions() {
        assert_eq!(
            s("inclusion_exclusion(list(list(10, 20, 30), list(5, 5, 5), list(1)))"),
            "46"
        );
        assert_eq!(s("inclusion_exclusion(list(list(7)))"), "7");
        assert_eq!(s("inclusion_exclusion(list())"), "0");
        assert_eq!(
            s("inclusion_exclusion(list(list(a, b), list(c)))"),
            "a + b - c"
        );
        // Numbers not divisible by 2, 3 or 5 up to 30 (Euler's totient of 30
        // times one): |2 or 3 or 5| = 15 + 10 + 6 - 5 - 3 - 2 + 1 = 22.
        assert_eq!(
            s("30 - inclusion_exclusion(list(list(15, 10, 6), list(5, 3, 2), list(1)))"),
            "8"
        );
        assert_eq!(s("inclusion_exclusion(x)"), "inclusion_exclusion(x)");
        assert_eq!(
            s("inclusion_exclusion(list(1, 2))"),
            "inclusion_exclusion(list(1, 2))"
        );
    }

    #[test]
    fn periods_of_sequences() {
        assert_eq!(s("period(list(1, 2, 1, 2))"), "2");
        assert_eq!(s("period(list(1, 2, 3, 1, 2, 3, 1, 2, 3))"), "3");
        assert_eq!(s("period(list(4, 4, 4))"), "1");
        assert_eq!(s("period(list(1, 2, 3))"), "3");
        assert_eq!(s("period(list(1, 2, 1))"), "3");
        assert_eq!(s("period(list(a, b, a, b, a, b))"), "2");
        assert_eq!(s("period(list(a))"), "1");
        assert_eq!(s("period(list())"), "period(list())");
        assert_eq!(s("period(x)"), "period(x)");
    }

    #[test]
    fn symbolic_binomial_expansion() {
        let rules = crate::rules::standard();
        assert_eq!(simplify(&rules, "expand_binomial(1, 1, n)"), "2^n");
        assert_eq!(
            simplify(&rules, "expand_binomial(x, 1, 3)"),
            "x^3 + 3*x^2 + 3*x + 1"
        );
        assert_eq!(
            simplify(&rules, "expand_binomial(x, y, 2)"),
            "x^2 + 2*x*y + y^2"
        );
        assert_eq!(simplify(&rules, "expand_binomial(2, 1, n)"), "3^n");
        // Inside a form request the definition is unfolded first.
        assert_eq!(
            simplify(&rules, "expand(expand_binomial(x, 1, 3))"),
            "x^3 + 3*x^2 + 3*x + 1"
        );
    }

    #[test]
    fn binomial_expansion_is_the_expand_request() {
        let sets = [combinatorics(), crate::rules::poly()];
        let e = |src: &str| simplify(&sets, src);
        assert_eq!(
            e("expand((a + b)^4)"),
            "a^4 + 4*a^3*b + 6*a^2*b^2 + 4*a*b^3 + b^4"
        );
        assert_eq!(e("expand((x - 1)^3)"), "x^3 - 3*x^2 + 3*x - 1");
        // The coefficients are the binomial coefficients of the kernel.
        for n in 0..8_u32 {
            let row: Vec<String> = (0..=n).map(|k| s(&format!("binomial({n}, {k})"))).collect();
            let expanded = e(&format!("expand((x + 1)^{n})"));
            for (k, c) in row.iter().enumerate() {
                let power = n as usize - k;
                let term = match (c.as_str(), power) {
                    | (c, 0) => c.to_owned(),
                    | ("1", 1) => "x".to_owned(),
                    | (c, 1) => format!("{c}*x"),
                    | ("1", p) => format!("x^{p}"),
                    | (c, p) => format!("{c}*x^{p}"),
                };
                assert!(expanded.contains(&term), "{expanded} lacks {term}");
            }
        }
    }

    #[test]
    fn counting_functions_match_the_legacy_ones() {
        // Old symbolic and numerical entry points, by value.
        assert_eq!(s("binomial(10, 3)"), "120");
        assert_eq!(s("permutations(10, 3)"), "720");
        assert_eq!(s("catalan(5)"), "42");
        assert_eq!(s("stirling2(5, 2)"), "15");
        assert_eq!(s("bell(6)"), "203");
        assert_eq!(s("rising(3, 4)"), "360");
        assert_eq!(s("falling(6, 3)"), "120");
        for n in 0..10_u64 {
            close(
                s(&format!("catalan({n})")).parse().unwrap_or(f64::NAN),
                crate::kernels::combinatorics::catalan(n),
            );
        }
    }
}

// More counting sequences and combinatorial structures: Lah, Eulerian,
// Narayana, Motzkin and Schröder numbers, restricted partitions and
// compositions, necklaces, Bernoulli and Euler numbers, and permutation
// statistics (rank, unrank, inversions, cycle type) and Young tableaux.
//
// | operator | value |
// |---|---|
// | `lah(n, k)` | Lah number `C(n-1, k-1) n!/k!`, the ways to split `n` items into `k` ordered lists |
// | `eulerian(n, k)` | Eulerian number: permutations of `n` with `k` ascents |
// | `narayana(n, k)` | `C(n, k) C(n, k-1) / n` |
// | `motzkin(n)`, `schroeder(n)` | Motzkin and large Schröder numbers |
// | `fubini(n)` | the ordered Bell number `Σ k! S(n, k)` |
// | `central_binomial(n)` | `C(2n, n)` |
// | `ballot(a, b)` | lattice paths with `a` up and `b` down steps staying non-negative |
// | `involutions(n)` | permutations with `p² = id` |
// | `rencontres(n, k)` | permutations of `n` with exactly `k` fixed points |
// | `partitions_k(n, k)` | partitions of `n` into exactly `k` parts |
// | `partitions_max(n, k)` | partitions of `n` with every part at most `k` |
// | `partitions_distinct(n)` | partitions of `n` into distinct parts |
// | `partition_list(n)` | the partitions themselves, in reverse lexicographic order (`n <= 40`) |
// | `compositions(n, k)`, `weak_compositions(n, k)` | ordered sums of `k` positive (non-negative) integers equal to `n` |
// | `necklaces(n, k)`, `bracelets(n, k)` | `k`-colourings of `n` beads up to rotation (and reflection), by Burnside's lemma |
// | `bernoulli_number(n)` | the Bernoulli number `B_n` (`B_1 = -1/2`) |
// | `euler_number(n)`, `zigzag(n)` | the secant numbers `E_n` (signed, zero for odd `n`) and the up-down numbers |
// | `permutation_rank(p)`, `permutation_unrank(n, r)` | lexicographic rank (from 0) of a permutation of `1..n` and its inverse |
// | `inversions(p)`, `permutation_sign(p)`, `cycle_type(p)`, `permutation_order(p)` | inversion count, sign, sorted cycle lengths and order |
// | `hook_lengths(lambda)`, `syt_count(lambda)` | hook lengths of a Young diagram and the number of standard tableaux (hook length formula) |


const SET: &str = "combinatorics";
const MAX_N: u64 = 2000;

fn big(v: u64) -> BigInt {
    BigInt::from(v)
}

fn binom(
    n: u64,
    k: u64,
) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    choose(&big(n), k).unwrap_or_default()
}

fn eulerian(
    n: u64,
    k: u64,
) -> BigInt {
    let mut row = vec![BigInt::one()];
    for m in 2..=n {
        let mut next = vec![BigInt::zero(); row.len() + 1];
        for (j, slot) in next.iter_mut().enumerate() {
            let j64 = u64::try_from(j).unwrap_or(0);
            let above = row.get(j).cloned().unwrap_or_default();
            let left = if j == 0 {
                BigInt::zero()
            } else {
                row.get(j - 1).cloned().unwrap_or_default()
            };
            *slot = above * (j64 + 1) + left * (m - j64);
        }
        row = next;
    }
    usize::try_from(k)
        .ok()
        .and_then(|k| row.get(k).cloned())
        .unwrap_or_default()
}

fn motzkin(n: u64) -> BigInt {
    let (mut previous, mut current) = (BigInt::one(), BigInt::one());
    if n == 0 {
        return previous;
    }
    for m in 2..=n {
        let next = (&current * (2 * m + 1) + &previous * (3 * m - 3)) / (m + 2);
        previous = current;
        current = next;
    }
    current
}

fn schroeder(n: u64) -> BigInt {
    let (mut previous, mut current) = (BigInt::one(), BigInt::from(2));
    if n == 0 {
        return previous;
    }
    for m in 2..=n {
        let next = (&current * (6 * m - 3) - &previous * (m - 2)) / (m + 1);
        previous = current;
        current = next;
    }
    current
}

fn fubini(n: u64) -> BigInt {
    let row = stirling_row(n, |_, j| j);
    row.iter()
        .enumerate()
        .map(|(k, s)| s * factorial(u64::try_from(k).unwrap_or(0)))
        .sum()
}

/// Partitions of `n` into exactly `k` parts, by `p(n, k) = p(n-1, k-1) + p(n-k, k)`.
fn partitions_exact(
    n: usize,
    k: usize,
) -> BigInt {
    let mut table = vec![vec![BigInt::zero(); k + 1]; n + 1];
    table[0][0] = BigInt::one();
    for m in 1..=n {
        for j in 1..=k.min(m) {
            let a = table[m - 1][j - 1].clone();
            let b = if m >= j {
                table[m - j][j].clone()
            } else {
                BigInt::zero()
            };
            table[m][j] = a + b;
        }
    }
    table[n][k].clone()
}

fn partitions_distinct(n: usize) -> BigInt {
    let mut ways = vec![BigInt::zero(); n + 1];
    ways[0] = BigInt::one();
    for part in 1..=n {
        for total in (part..=n).rev() {
            let add = ways[total - part].clone();
            ways[total] += add;
        }
    }
    ways[n].clone()
}

fn partition_list(n: u64) -> Vec<Vec<BigInt>> {
    fn go(
        remaining: u64,
        bound: u64,
        prefix: &mut Vec<u64>,
        out: &mut Vec<Vec<BigInt>>,
    ) {
        if remaining == 0 {
            out.push(prefix.iter().map(|&p| BigInt::from(p)).collect());
            return;
        }
        for part in (1..=remaining.min(bound)).rev() {
            prefix.push(part);
            go(remaining - part, part, prefix, out);
            prefix.pop();
        }
    }
    let mut out = Vec::new();
    go(n, n, &mut Vec::new(), &mut out);
    out
}

fn necklaces(
    n: u64,
    k: &BigInt,
) -> Option<BigInt> {
    if n == 0 {
        return Some(BigInt::one());
    }
    let mut total = BigInt::zero();
    for d in 1..=n {
        if n.is_multiple_of(d) {
            let phi = (1..=d).filter(|&x| gcd_u64(x, d) == 1).count();
            total +=
                BigInt::from(phi) * num_traits::pow::pow(k.clone(), usize::try_from(n / d).ok()?);
        }
    }
    Some(total / n)
}

const fn gcd_u64(
    mut a: u64,
    mut b: u64,
) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn bracelets(
    n: u64,
    k: &BigInt,
) -> Option<BigInt> {
    if n == 0 {
        return Some(BigInt::one());
    }
    let rotations = necklaces(n, k)?;
    let exponent = |e: u64| num_traits::pow::pow(k.clone(), usize::try_from(e).ok().unwrap_or(0));
    let reflections = if n % 2 == 1 {
        BigInt::from(n) * exponent(n.div_ceil(2))
    } else {
        BigInt::from(n / 2) * (exponent(n / 2 + 1) + exponent(n / 2))
    };
    // (rotations * n + reflections) / 2n
    Some((rotations * n + reflections) / (2 * n))
}

/// Bernoulli numbers by the Akiyama–Tanigawa transform (`B_1 = -1/2`).
fn bernoulli(n: usize) -> BigRational {
    let mut a: Vec<BigRational> = Vec::with_capacity(n + 1);
    for m in 0..=n {
        a.push(BigRational::new(BigInt::one(), BigInt::from(m + 1)));
        for j in (1..=m).rev() {
            let value =
                (a[j - 1].clone() - a[j].clone()) * BigRational::from_integer(BigInt::from(j));
            a[j - 1] = value;
        }
    }
    let value = a.first().cloned().unwrap_or_default();
    if n == 1 { -value } else { value }
}

/// Up-down (Euler zigzag) numbers by the Seidel–Entringer triangle.
fn zigzag(n: usize) -> BigInt {
    let mut row = vec![BigInt::one()];
    for m in 1..=n {
        let mut next = vec![BigInt::zero()];
        for k in 1..=m {
            let value = next[k - 1].clone() + row[m - k].clone();
            next.push(value);
        }
        row = next;
    }
    row.last().cloned().unwrap_or_default()
}

/// A permutation of `1..n` read from a list of integers.
fn permutation(items: &[BigInt]) -> Option<Vec<usize>> {
    let n = items.len();
    let mut seen = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for v in items {
        let v = usize::try_from(v.to_u64()?).ok()?;
        if v == 0 || v > n || seen[v - 1] {
            return None;
        }
        seen[v - 1] = true;
        out.push(v - 1);
    }
    Some(out)
}

fn cycle_lengths(perm: &[usize]) -> Vec<u64> {
    let mut seen = vec![false; perm.len()];
    let mut lengths = Vec::new();
    for start in 0..perm.len() {
        if seen[start] {
            continue;
        }
        let (mut at, mut length) = (start, 0_u64);
        while !seen[at] {
            seen[at] = true;
            at = perm[at];
            length += 1;
        }
        lengths.push(length);
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    lengths
}

fn inversions(perm: &[usize]) -> u64 {
    let mut total = 0_u64;
    for (i, &a) in perm.iter().enumerate() {
        total += u64::try_from(perm[i + 1..].iter().filter(|&&b| b < a).count()).unwrap_or(0);
    }
    total
}

fn partition_shape(items: &[BigInt]) -> Option<Vec<usize>> {
    let shape: Vec<usize> = items
        .iter()
        .map(|v| v.to_u64().and_then(|v| usize::try_from(v).ok()))
        .collect::<Option<_>>()?;
    let valid = shape.windows(2).all(|w| w[0] >= w[1]) && shape.iter().all(|&p| p > 0);
    (valid && shape.len() <= 60 && shape.iter().sum::<usize>() <= 400).then_some(shape)
}

fn hooks(shape: &[usize]) -> Vec<Vec<u64>> {
    shape
        .iter()
        .enumerate()
        .map(|(row, &length)| {
            (0..length)
                .map(|col| {
                    let arm = length - col - 1;
                    let leg = shape[row + 1..]
                        .iter()
                        .filter(|&&below| below > col)
                        .count();
                    u64::try_from(arm + leg + 1).unwrap_or(0)
                })
                .collect()
        })
        .collect()
}

#[allow(clippy::too_many_lines)]
fn install_ext(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));

    exact(i, SET, binary("lah"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if n == 0 || k == 0 || k > n {
            return Some(Value::Int(if n == 0 && k == 0 {
                BigInt::one()
            } else {
                BigInt::zero()
            }));
        }
        Some(Value::Int(
            binom(n - 1, k - 1) * factorial(n) / factorial(k),
        ))
    })?;
    exact(i, SET, binary("eulerian"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 500)?, nat(k, 500)?);
        if n == 0 {
            return Some(Value::Int(if k == 0 {
                BigInt::one()
            } else {
                BigInt::zero()
            }));
        }
        Some(Value::Int(eulerian(n, k)))
    })?;
    exact(i, SET, binary("narayana"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if n == 0 || k == 0 || k > n {
            return Some(Value::Int(BigInt::zero()));
        }
        Some(Value::Int(binom(n, k) * binom(n, k - 1) / n))
    })?;
    exact(i, SET, unary("motzkin"), |a| {
        match a {
            | [n] => Some(Value::Int(motzkin(nat(n, MAX_N)?))),
            | _ => None,
        }
    })?;
    exact(i, SET, unary("schroeder"), |a| {
        match a {
            | [n] => Some(Value::Int(schroeder(nat(n, MAX_N)?))),
            | _ => None,
        }
    })?;
    exact(i, SET, unary("fubini"), |a| {
        match a {
            | [n] => Some(Value::Int(fubini(nat(n, 400)?))),
            | _ => None,
        }
    })?;
    exact(i, SET, unary("central_binomial"), |a| {
        match a {
            | [n] => Some(Value::Int(binom(2 * nat(n, MAX_N)?, nat(n, MAX_N)?))),
            | _ => None,
        }
    })?;
    exact(i, SET, binary("ballot"), |a| {
        let [up, down] = a else { return None };
        let (up, down) = (nat(up, MAX_N)?, nat(down, MAX_N)?);
        if down > up {
            return Some(Value::Int(BigInt::zero()));
        }
        Some(Value::Int(
            binom(up + down, down) * (up - down + 1) / (up + 1),
        ))
    })?;
    exact(i, SET, unary("involutions"), |a| {
        let [n] = a else { return None };
        let n = nat(n, MAX_N)?;
        let (mut previous, mut current) = (BigInt::one(), BigInt::one());
        for m in 2..=n {
            let next = &current + &previous * (m - 1);
            previous = current;
            current = next;
        }
        Some(Value::Int(if n == 0 { previous } else { current }))
    })?;
    exact(i, SET, binary("rencontres"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if k > n {
            return Some(Value::Int(BigInt::zero()));
        }
        let m = n - k;
        let (mut previous, mut current) = (BigInt::one(), BigInt::zero());
        let derangements = match m {
            | 0 => previous,
            | _ => {
                for j in 2..=m {
                    let next = (&previous + &current) * (j - 1);
                    previous = current;
                    current = next;
                }
                current
            },
        };
        Some(Value::Int(binom(n, k) * derangements))
    })?;
    exact(i, SET, binary("partitions_k"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 600)?, nat(k, 600)?);
        Some(Value::Int(partitions_exact(
            usize::try_from(n).ok()?,
            usize::try_from(k).ok()?,
        )))
    })?;
    exact(i, SET, binary("partitions_max"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 600)?, nat(k, 600)?);
        // Conjugation: parts at most k <-> at most k parts.
        let (n, k) = (usize::try_from(n).ok()?, usize::try_from(k).ok()?);
        let total: BigInt = (0..=k).map(|j| partitions_exact(n, j)).sum();
        Some(Value::Int(total))
    })?;
    exact(i, SET, unary("partitions_distinct"), |a| {
        match a {
            | [n] => {
                Some(Value::Int(partitions_distinct(
                    usize::try_from(nat(n, 2000)?).ok()?,
                )))
            },
            | _ => None,
        }
    })?;
    exact(i, SET, unary("partition_list"), |a| {
        match a {
            | [n] => Some(Value::Nested(partition_list(nat(n, 40)?))),
            | _ => None,
        }
    })?;
    exact(i, SET, binary("compositions"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        Some(Value::Int(if n == 0 || k == 0 {
            if n == 0 && k == 0 {
                BigInt::one()
            } else {
                BigInt::zero()
            }
        } else {
            binom(n - 1, k - 1)
        }))
    })?;
    exact(i, SET, binary("weak_compositions"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        Some(Value::Int(if k == 0 {
            if n == 0 {
                BigInt::one()
            } else {
                BigInt::zero()
            }
        } else {
            binom(n + k - 1, k - 1)
        }))
    })?;
    exact(i, SET, binary("necklaces"), |a| {
        let [n, k] = a else { return None };
        necklaces(nat(n, 400)?, k).map(Value::Int)
    })?;
    exact(i, SET, binary("bracelets"), |a| {
        let [n, k] = a else { return None };
        bracelets(nat(n, 400)?, k).map(Value::Int)
    })?;
    exact(i, SET, unary("bernoulli_number"), |a| {
        let [n] = a else { return None };
        Some(Value::Rat(bernoulli(usize::try_from(nat(n, 300)?).ok()?)))
    })?;
    exact(i, SET, unary("zigzag"), |a| {
        let [n] = a else { return None };
        Some(Value::Int(zigzag(usize::try_from(nat(n, 1000)?).ok()?)))
    })?;
    exact(i, SET, unary("euler_number"), |a| {
        let [n] = a else { return None };
        let n = nat(n, 1000)?;
        if n % 2 == 1 {
            return Some(Value::Int(BigInt::zero()));
        }
        let z = zigzag(usize::try_from(n).ok()?);
        Some(Value::Int(if n % 4 == 2 { -z } else { z }))
    })?;
    exact_lists(i, SET, unary("permutation_rank"), |a| {
        let [items] = a else { return None };
        let perm = permutation(items)?;
        let n = perm.len();
        let mut rank = BigInt::zero();
        for (idx, &value) in perm.iter().enumerate() {
            let smaller = perm[idx + 1..]
                .iter()
                .filter(|&&later| later < value)
                .count();
            rank += BigInt::from(smaller) * factorial(u64::try_from(n - idx - 1).ok()?);
        }
        Some(Value::Int(rank))
    })?;
    exact(i, SET, binary("permutation_unrank"), |a| {
        let [n, rank] = a else { return None };
        let n = nat(n, 500)?;
        if rank.is_negative() || *rank >= factorial(n) {
            return None;
        }
        let mut pool: Vec<u64> = (1..=n).collect();
        let mut rank = rank.clone();
        let mut out = Vec::new();
        for place in (0..n).rev() {
            let base = factorial(place);
            let index = (&rank / &base).to_usize()?;
            rank %= base;
            out.push(big(pool.remove(index)));
        }
        Some(Value::Flat(out))
    })?;
    exact_lists(i, SET, unary("inversions"), |a| {
        let [items] = a else { return None };
        Some(Value::Int(big(inversions(&permutation(items)?))))
    })?;
    exact_lists(i, SET, unary("permutation_sign"), |a| {
        let [items] = a else { return None };
        let even = inversions(&permutation(items)?).is_multiple_of(2);
        Some(Value::Int(BigInt::from(if even { 1 } else { -1 })))
    })?;
    exact_lists(i, SET, unary("cycle_type"), |a| {
        let [items] = a else { return None };
        Some(Value::Flat(
            cycle_lengths(&permutation(items)?)
                .into_iter()
                .map(BigInt::from)
                .collect(),
        ))
    })?;
    exact_lists(i, SET, unary("permutation_order"), |a| {
        let [items] = a else { return None };
        let lcm = cycle_lengths(&permutation(items)?)
            .into_iter()
            .fold(BigInt::one(), |acc, l| {
                &acc / num_integer::gcd(acc.clone(), BigInt::from(l)) * l
            });
        Some(Value::Int(lcm))
    })?;
    exact_lists(i, SET, unary("hook_lengths"), |a| {
        let [items] = a else { return None };
        let shape = partition_shape(items)?;
        Some(Value::Nested(
            hooks(&shape)
                .into_iter()
                .map(|row| row.into_iter().map(BigInt::from).collect())
                .collect(),
        ))
    })?;
    exact_lists(i, SET, unary("syt_count"), |a| {
        let [items] = a else { return None };
        let shape = partition_shape(items)?;
        let cells: usize = shape.iter().sum();
        let denominator: BigInt = hooks(&shape)
            .into_iter()
            .flatten()
            .map(BigInt::from)
            .product();
        Some(Value::Int(
            factorial(u64::try_from(cells).ok()?) / denominator,
        ))
    })?;
    Ok(())
}

#[cfg(test)]
mod ext_tests {
    use crate::rules::combinatorics::combinatorics;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[combinatorics()], src)
    }

    #[test]
    fn classical_triangles_and_sequences() {
        assert_eq!(s("lah(4, 2)"), "36");
        assert_eq!(s("lah(5, 1)"), "120");
        assert_eq!(s("eulerian(4, 1)"), "11");
        assert_eq!(s("eulerian(5, 2)"), "66");
        assert_eq!(s("narayana(4, 2)"), "6");
        assert_eq!(s("motzkin(6)"), "51");
        assert_eq!(s("schroeder(4)"), "90");
        assert_eq!(s("fubini(4)"), "75");
        assert_eq!(s("central_binomial(5)"), "252");
        assert_eq!(s("ballot(4, 2)"), "9");
        assert_eq!(s("involutions(5)"), "26");
        assert_eq!(s("rencontres(5, 2)"), "20");
    }

    #[test]
    fn restricted_partitions_and_compositions() {
        assert_eq!(s("partitions_k(7, 3)"), "4");
        assert_eq!(s("partitions_max(7, 3)"), "8");
        assert_eq!(s("partitions_distinct(10)"), "10");
        assert_eq!(
            s("partition_list(4)"),
            "list(list(4), list(3, 1), list(2, 2), list(2, 1, 1), list(1, 1, 1, 1))"
        );
        assert_eq!(s("compositions(5, 3)"), "6");
        assert_eq!(s("weak_compositions(5, 3)"), "21");
        assert_eq!(s("necklaces(6, 2)"), "14");
        assert_eq!(s("bracelets(6, 2)"), "13");
    }

    #[test]
    fn bernoulli_and_euler_numbers() {
        assert_eq!(s("bernoulli_number(0)"), "1");
        assert_eq!(s("bernoulli_number(1)"), "-1/2");
        assert_eq!(s("bernoulli_number(2)"), "1/6");
        assert_eq!(s("bernoulli_number(12)"), "-691/2730");
        assert_eq!(s("bernoulli_number(13)"), "0");
        assert_eq!(s("zigzag(6)"), "61");
        assert_eq!(s("euler_number(4)"), "5");
        assert_eq!(s("euler_number(6)"), "-61");
    }

    #[test]
    fn permutation_statistics_and_tableaux() {
        assert_eq!(s("permutation_rank(list(2, 3, 1))"), "3");
        assert_eq!(s("permutation_unrank(3, 3)"), "list(2, 3, 1)");
        assert_eq!(s("inversions(list(3, 1, 2))"), "2");
        assert_eq!(s("permutation_sign(list(2, 1, 3))"), "-1");
        assert_eq!(s("cycle_type(list(2, 1, 4, 5, 3))"), "list(3, 2)");
        assert_eq!(s("permutation_order(list(2, 1, 4, 5, 3))"), "6");
        assert_eq!(
            s("hook_lengths(list(3, 2))"),
            "list(list(4, 3, 1), list(2, 1))"
        );
        assert_eq!(s("syt_count(list(3, 2))"), "5");
        assert_eq!(s("syt_count(list(3, 2, 1))"), "16");
    }
}
