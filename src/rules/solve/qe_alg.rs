//! Exact univariate polynomials over `Q` and real algebraic numbers.
//!
//! A [`RealAlg`] is a real algebraic number represented by
//!
//! * a square-free polynomial `p` over `Q` with positive leading
//!   coefficient, and
//! * a rational interval `[lo, hi]` that contains exactly one real root of
//!   `p` (a point interval when that root is rational).
//!
//! `p` need not be irreducible: the root is selected by the interval, and
//! every question is decided without floating point. The sign of an
//! arbitrary `f` over `Q` at the root `a` is found by reducing modulo `p`:
//! `f(a) = 0` exactly when `gcd(f, p)` has a root in the isolating interval
//! (a Sturm count), and otherwise interval Horner evaluation is refined
//! until its enclosure of `f(a)` excludes zero. Comparison refines both
//! intervals; when they overlap, the two numbers are equal exactly when the
//! gcd of the defining polynomials has a common root in the intersection,
//! again a Sturm count. This is the number type used by the cylindrical
//! algebraic decomposition and by the one-variable decisions in the tests.
//!
//! The univariate helpers (division, gcd, Sturm chains, isolation,
//! resultants) are exact over [`BigRational`]; the Sturm chains are kept
//! primitive so that coefficient growth stays manageable.

// The algebraic-number API is documented and exercised by the unit tests;
// comparison, isolation and refinement are part of the representation's
// contract even where the QE hot path does not use them.
#![allow(
    dead_code,
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;

/// Bound on the bisection steps of a refinement.
const MAX_REFINE: usize = 4096;
/// Bound on the steps of the isolation loop.
const MAX_ISOLATE: usize = 200_000;
/// Bound on the number of terms of an intermediate univariate product.
const MAX_TERMS: usize = 100_000;

/// Removes trailing zero coefficients.
pub(crate) fn trim(p: &mut Vec<BigRational>) {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
}

/// Degree of a polynomial, `None` for zero.
#[must_use]
pub(crate) fn degree(p: &[BigRational]) -> Option<usize> {
    let mut n = p.len();
    while n > 0 && p.get(n.saturating_sub(1)).is_some_and(Zero::is_zero) {
        n = n.saturating_sub(1);
    }
    (n > 0).then(|| n.saturating_sub(1))
}

/// Sum.
#[must_use]
pub(crate) fn add(
    a: &[BigRational],
    b: &[BigRational],
) -> Vec<BigRational> {
    let n = a.len().max(b.len());
    let mut out = vec![BigRational::zero(); n];
    for (i, slot) in out.iter_mut().enumerate() {
        let x = a.get(i).cloned().unwrap_or_else(BigRational::zero);
        let y = b.get(i).cloned().unwrap_or_else(BigRational::zero);
        *slot = x + y;
    }
    trim(&mut out);
    out
}

/// Difference.
#[must_use]
pub(crate) fn sub(
    a: &[BigRational],
    b: &[BigRational],
) -> Vec<BigRational> {
    let n = a.len().max(b.len());
    let mut out = vec![BigRational::zero(); n];
    for (i, slot) in out.iter_mut().enumerate() {
        let x = a.get(i).cloned().unwrap_or_else(BigRational::zero);
        let y = b.get(i).cloned().unwrap_or_else(BigRational::zero);
        *slot = x - y;
    }
    trim(&mut out);
    out
}

/// Product, or `None` above the term cap.
#[must_use]
pub(crate) fn mul(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<Vec<BigRational>> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let n = a.len().saturating_add(b.len()).saturating_sub(1);
    if n > MAX_TERMS {
        return None;
    }
    let mut out = vec![BigRational::zero(); n];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i.saturating_add(j)) {
                *slot += x * y;
            }
        }
    }
    trim(&mut out);
    Some(out)
}

/// Scale by a rational.
#[must_use]
pub(crate) fn scale(
    a: &[BigRational],
    factor: &BigRational,
) -> Vec<BigRational> {
    if factor.is_zero() {
        return Vec::new();
    }
    let mut out: Vec<BigRational> = a.iter().map(|c| c * factor).collect();
    trim(&mut out);
    out
}

/// Quotient and remainder of `a` by `b` over `Q`.
#[must_use]
pub(crate) fn divrem(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(Vec<BigRational>, Vec<BigRational>)> {
    let bn = degree(b)?;
    let b_lead = b.get(bn)?.clone();
    if b_lead.is_zero() {
        return None;
    }
    let mut rest = a.to_vec();
    trim(&mut rest);
    if degree(&rest).is_none() {
        return Some((Vec::new(), Vec::new()));
    }
    let an = degree(&rest)?;
    if an < bn {
        return Some((Vec::new(), rest));
    }
    let mut quotient = vec![BigRational::zero(); an.saturating_sub(bn).saturating_add(1)];
    for step in 0..=an.saturating_sub(bn) {
        let power = an.saturating_sub(bn).saturating_sub(step);
        let rest_degree = match degree(&rest) {
            | Some(d) if d >= bn => d,
            | _ => break,
        };
        let shift = rest_degree.saturating_sub(bn);
        let lead = rest.get(rest_degree)?.clone();
        if b_lead.is_zero() {
            return None;
        }
        let factor = lead / &b_lead;
        if let Some(slot) = quotient.get_mut(shift) {
            *slot = factor.clone();
        }
        let scaled = scale(b, &factor);
        let shifted: Vec<BigRational> = std::iter::repeat_n(BigRational::zero(), shift)
            .chain(scaled)
            .collect();
        rest = sub(&rest, &shifted);
        trim(&mut rest);
        if step > MAX_TERMS {
            return None;
        }
        if power == 0 {
            break;
        }
    }
    Some((quotient, rest))
}

/// Monic form; the zero polynomial stays zero.
#[must_use]
pub(crate) fn monic(p: &[BigRational]) -> Vec<BigRational> {
    match degree(p) {
        | None => Vec::new(),
        | Some(n) => {
            match p.get(n) {
                | Some(lead) if !lead.is_zero() => scale(p, &lead.recip()),
                | _ => p.to_vec(),
            }
        },
    }
}

/// Monic greatest common divisor over `Q` (Euclid).
#[must_use]
pub(crate) fn gcd(
    a: &[BigRational],
    b: &[BigRational],
) -> Vec<BigRational> {
    let mut x = monic(a);
    let mut y = monic(b);
    let mut steps = 0_usize;
    while degree(&y).is_some() {
        let Some((_, r)) = divrem(&x, &y) else {
            break;
        };
        x = y;
        y = monic(&r);
        steps = steps.saturating_add(1);
        if steps > MAX_TERMS {
            break;
        }
    }
    x
}

/// Derivative.
#[must_use]
pub(crate) fn derivative(p: &[BigRational]) -> Vec<BigRational> {
    let mut out = Vec::new();
    for (i, c) in p.iter().enumerate().skip(1) {
        let factor = BigRational::from_integer(BigInt::from(i));
        out.push(c * factor);
    }
    trim(&mut out);
    out
}

/// The square-free part `p / gcd(p, p')`, monic.
#[must_use]
pub(crate) fn square_free_part(p: &[BigRational]) -> Vec<BigRational> {
    let d = derivative(p);
    let g = gcd(p, &d);
    if degree(&g).is_some() {
        divrem(p, &g).map_or_else(|| monic(p), |(q, _)| monic(&q))
    } else {
        monic(p)
    }
}

/// Value at a rational point (Horner).
#[must_use]
pub(crate) fn eval(
    p: &[BigRational],
    x: &BigRational,
) -> BigRational {
    let mut acc = BigRational::zero();
    for c in p.iter().rev() {
        acc = acc * x + c;
    }
    acc
}

fn interval_mul(
    a: (&BigRational, &BigRational),
    b: (&BigRational, &BigRational),
) -> (BigRational, BigRational) {
    let products = [a.0 * b.0, a.0 * b.1, a.1 * b.0, a.1 * b.1];
    let mut lo = products.first().cloned().unwrap_or_else(BigRational::zero);
    let mut hi = lo.clone();
    for p in products.iter().skip(1) {
        if *p < lo {
            lo = p.clone();
        }
        if *p > hi {
            hi = p.clone();
        }
    }
    (lo, hi)
}

/// An enclosure of the values of `p` on `[lo, hi]` by interval Horner.
#[must_use]
pub(crate) fn eval_interval(
    p: &[BigRational],
    lo: &BigRational,
    hi: &BigRational,
) -> (BigRational, BigRational) {
    let mut acc = (BigRational::zero(), BigRational::zero());
    for c in p.iter().rev() {
        let (a, b) = interval_mul((&acc.0, &acc.1), (lo, hi));
        acc = (a + c, b + c);
    }
    if acc.0 > acc.1 {
        (acc.1, acc.0)
    } else {
        acc
    }
}

/// Keeps a Sturm chain entry primitive with positive leading coefficient.
/// The sign variations of a chain are unchanged by a positive factor.
#[must_use]
pub(crate) fn primitive(p: &[BigRational]) -> Vec<BigRational> {
    let Some(n) = degree(p) else {
        return Vec::new();
    };
    let mut denominator = BigInt::one();
    for c in p {
        denominator = denominator.lcm(c.denom());
    }
    let mut numerator = BigInt::one();
    let mut integers: Vec<BigInt> = Vec::with_capacity(n.saturating_add(1));
    for c in p {
        let scaled = c.numer() * (&denominator / c.denom());
        numerator = numerator.gcd(&scaled);
        integers.push(scaled);
    }
    if numerator.is_zero() {
        return Vec::new();
    }
    // Scaling by a negative factor would change the signs of the Sturm
    // variations, so the sign is deliberately not normalised here.
    integers
        .into_iter()
        .map(|c| BigRational::from_integer(c / &numerator))
        .collect()
}

/// The Sturm chain `p, p', -rem(p_{i-1}, p_i), ...` with primitive entries.
#[must_use]
pub(crate) fn sturm_chain(p: &[BigRational]) -> Vec<Vec<BigRational>> {
    let mut chain = vec![primitive(p)];
    let d = primitive(&derivative(p));
    if degree(&d).is_some() {
        chain.push(d);
    }
    let mut steps = 0_usize;
    loop {
        let n = chain.len();
        let (Some(a), Some(b)) = (
            chain.get(n.saturating_sub(2)),
            chain.get(n.saturating_sub(1)),
        ) else {
            break;
        };
        let Some((_, rem)) = divrem(a, b) else {
            break;
        };
        if degree(&rem).is_none() {
            break;
        }
        chain.push(primitive(&scale(
            &rem,
            &BigRational::from_integer(BigInt::from(-1)),
        )));
        steps = steps.saturating_add(1);
        if steps > MAX_TERMS {
            break;
        }
    }
    chain.retain(|q| degree(q).is_some());
    chain
}

fn sign_of(value: &BigRational) -> i32 {
    if value.is_zero() {
        0
    } else if value.is_positive() {
        1
    } else {
        -1
    }
}

/// Sign of `q` at `x`, or at `±∞` when `x` is `None`.
fn sign_at(
    q: &[BigRational],
    at: Option<&BigRational>,
    positive: bool,
) -> i32 {
    match at {
        | Some(x) => sign_of(&eval(q, x)),
        | None => {
            let lead = q.last().map_or(0, sign_of);
            let odd = q.len().is_multiple_of(2);
            if positive || !odd {
                lead
            } else {
                -lead
            }
        },
    }
}

fn variations(
    chain: &[Vec<BigRational>],
    at: Option<&BigRational>,
    positive: bool,
) -> usize {
    let signs: Vec<i32> = chain
        .iter()
        .map(|q| sign_at(q, at, positive))
        .filter(|&s| s != 0)
        .collect();
    signs.windows(2).filter(|w| w.first() != w.get(1)).count()
}

/// Distinct real roots of the square-free `p` in `(a, b]` (`None` ends are
/// infinite).
#[must_use]
pub(crate) fn sturm_count(
    chain: &[Vec<BigRational>],
    a: Option<&BigRational>,
    b: Option<&BigRational>,
) -> usize {
    variations(chain, a, false).saturating_sub(variations(chain, b, true))
}

/// Distinct real roots in `[a, b]`, counting the endpoints.
#[must_use]
pub(crate) fn count_closed(
    p: &[BigRational],
    a: &BigRational,
    b: &BigRational,
) -> usize {
    if a > b {
        return 0;
    }
    if a == b {
        return usize::from(eval(p, a).is_zero());
    }
    let chain = sturm_chain(p);
    let mut count = sturm_count(&chain, Some(a), Some(b));
    if eval(p, a).is_zero() {
        count = count.saturating_add(1);
    }
    count
}

/// Distinct real roots in `(a, b)`.
#[must_use]
pub(crate) fn count_open(
    p: &[BigRational],
    a: &BigRational,
    b: &BigRational,
) -> usize {
    if a >= b {
        return 0;
    }
    let chain = sturm_chain(p);
    let upper = sturm_count(&chain, Some(a), Some(b));
    if eval(p, b).is_zero() {
        upper.saturating_sub(1)
    } else {
        upper
    }
}

/// Cauchy bound: every real root has absolute value below this.
#[must_use]
pub(crate) fn cauchy_bound(p: &[BigRational]) -> BigRational {
    let Some(n) = degree(p) else {
        return BigRational::one();
    };
    let Some(lead) = p.get(n) else {
        return BigRational::one();
    };
    let mut max = BigRational::zero();
    for c in p.iter().take(n) {
        let ratio = (c / lead).abs();
        if ratio > max {
            max = ratio;
        }
    }
    max + BigRational::one()
}

/// Bisection step that keeps the unique root of the square-free `p` in
/// `[lo, hi]`.
#[must_use]
pub(crate) fn refine_once(
    p: &[BigRational],
    lo: &BigRational,
    hi: &BigRational,
) -> (BigRational, BigRational) {
    if lo >= hi {
        return (lo.clone(), hi.clone());
    }
    let two = BigRational::from_integer(BigInt::from(2));
    let mid = (lo + hi) / &two;
    let at_mid = eval(p, &mid);
    if at_mid.is_zero() {
        return (mid.clone(), mid);
    }
    let at_lo = eval(p, lo);
    if at_lo.is_negative() ^ at_mid.is_negative() {
        (lo.clone(), mid)
    } else {
        (mid, hi.clone())
    }
}

/// Isolating intervals of the distinct real roots of `p`, ascending. A
/// square-free polynomial is assumed; its square-free part is taken anyway.
#[must_use]
pub(crate) fn isolate_roots(p: &[BigRational]) -> Vec<(BigRational, BigRational)> {
    let p = square_free_part(p);
    if degree(&p).is_none() {
        return Vec::new();
    }
    let bound = cauchy_bound(&p);
    let minus = -bound.clone();
    let chain = sturm_chain(&p);
    let mut out: Vec<(BigRational, BigRational)> = Vec::new();
    let mut stack = vec![(minus, bound)];
    let mut work = 0_usize;
    while let Some((a, b)) = stack.pop() {
        work = work.saturating_add(1);
        if work > MAX_ISOLATE {
            break;
        }
        let count = sturm_count(&chain, Some(&a), Some(&b));
        match count {
            | 0 => {},
            | 1 => out.push((a, b)),
            | _ => {
                let two = BigRational::from_integer(BigInt::from(2));
                let mid = (&a + &b) / &two;
                if eval(&p, &mid).is_zero() {
                    out.push((mid.clone(), mid.clone()));
                    let left = count_open(&p, &a, &mid);
                    let right = count_open(&p, &mid, &b);
                    if left > 0 {
                        stack.push((a, mid.clone()));
                    }
                    if right > 0 {
                        stack.push((mid, b));
                    }
                } else {
                    let left = sturm_count(&chain, Some(&a), Some(&mid));
                    let right = sturm_count(&chain, Some(&mid), Some(&b));
                    if left > 0 {
                        stack.push((a, mid.clone()));
                    }
                    if right > 0 {
                        stack.push((mid, b));
                    }
                }
            },
        }
    }
    for interval in &mut out {
        // Sharpen the isolating intervals; a few bisections are cheap and
        // make the later algebraic lifting much better conditioned.
        for _ in 0..8 {
            let (lo, hi) = refine_once(&p, &interval.0, &interval.1);
            *interval = (lo, hi);
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    out
}

/// The Sylvester resultant of `p` and `q`.
#[must_use]
pub(crate) fn resultant(
    p: &[BigRational],
    q: &[BigRational],
) -> Option<BigRational> {
    let m = degree(p)?;
    let n = degree(q)?;
    if m == 0 && n == 0 {
        return Some(BigRational::one());
    }
    let size = m.saturating_add(n);
    let mut rows: Vec<Vec<BigRational>> = Vec::with_capacity(size);
    let descending_p: Vec<BigRational> =
        p.iter().take(m.saturating_add(1)).rev().cloned().collect();
    let descending_q: Vec<BigRational> =
        q.iter().take(n.saturating_add(1)).rev().cloned().collect();
    for shift in 0..n {
        let mut row = vec![BigRational::zero(); size];
        for (k, c) in descending_p.iter().enumerate() {
            if let Some(slot) = row.get_mut(shift.saturating_add(k)) {
                *slot = c.clone();
            }
        }
        rows.push(row);
    }
    for shift in 0..m {
        let mut row = vec![BigRational::zero(); size];
        for (k, c) in descending_q.iter().enumerate() {
            if let Some(slot) = row.get_mut(shift.saturating_add(k)) {
                *slot = c.clone();
            }
        }
        rows.push(row);
    }
    determinant(&mut rows)
}

/// Determinant by Gaussian elimination over `Q`.
fn determinant(rows: &mut [Vec<BigRational>]) -> Option<BigRational> {
    let n = rows.len();
    let mut det = BigRational::one();
    let mut flip = false;
    for k in 0..n {
        let pivot = (k..n).find(|&r| {
            rows.get(r)
                .and_then(|row| row.get(k))
                .is_some_and(|v| !v.is_zero())
        });
        let Some(pivot) = pivot else {
            // A singular matrix: the determinant is zero.
            return Some(BigRational::zero());
        };
        if pivot != k {
            rows.swap(k, pivot);
            flip = !flip;
        }
        let pivot_value = rows.get(k)?.get(k)?.clone();
        if pivot_value.is_zero() {
            return None;
        }
        det *= &pivot_value;
        for r in k.saturating_add(1)..n {
            let factor = rows.get(r)?.get(k)?.clone() / &pivot_value;
            if factor.is_zero() {
                continue;
            }
            let pivot_row = rows.get(k)?.clone();
            let row = rows.get_mut(r)?;
            for c in k..n {
                let value = pivot_row.get(c)?.clone();
                if let Some(slot) = row.get_mut(c) {
                    *slot -= &factor * &value;
                }
            }
        }
    }
    if flip {
        det = -det;
    }
    Some(det)
}

/// A real algebraic number: a square-free defining polynomial over `Q`
/// together with an isolating interval for one of its real roots.
///
/// See the module documentation for the exactness arguments.
#[derive(Clone, Debug)]
pub(crate) struct RealAlg {
    /// Square-free, positive leading coefficient, ascending order.
    poly: Vec<BigRational>,
    /// Isolating interval; `lo == hi` for a rational root.
    lo: BigRational,
    hi: BigRational,
}

impl RealAlg {
    /// A number given by a defining polynomial and an isolating interval.
    /// Returns `None` when the interval does not isolate exactly one root.
    #[must_use]
    pub(crate) fn from_root(
        poly: &[BigRational],
        lo: BigRational,
        hi: BigRational,
    ) -> Option<Self> {
        let sf = square_free_part(poly);
        let d = degree(&sf)?;
        if d == 0 || lo > hi {
            return None;
        }
        let count = count_closed(&sf, &lo, &hi);
        if count != 1 {
            return None;
        }
        let positive = sf.get(d).is_some_and(num_traits::Signed::is_positive);
        let sf = if positive {
            sf
        } else {
            scale(&sf, &BigRational::from_integer(BigInt::from(-1)))
        };
        // A degree-one polynomial has an exact rational root; normalise the
        // interval to a point so that the number is recognised as rational.
        if d == 1 {
            let c0 = sf.first()?.clone();
            let c1 = sf.get(1)?.clone();
            let value = -c0 / c1;
            if value < lo || value > hi {
                return None;
            }
            return Some(Self {
                poly: sf,
                lo: value.clone(),
                hi: value,
            });
        }
        Some(Self { poly: sf, lo, hi })
    }

    /// The defining polynomial.
    #[must_use]
    pub(crate) fn polynomial(&self) -> &[BigRational] {
        &self.poly
    }

    /// The isolating interval.
    #[must_use]
    pub(crate) fn interval(&self) -> (&BigRational, &BigRational) {
        (&self.lo, &self.hi)
    }

    /// Whether the root is rational and known exactly.
    #[must_use]
    pub(crate) fn rational_value(&self) -> Option<BigRational> {
        (self.lo == self.hi).then(|| self.lo.clone())
    }

    /// One bisection step towards the root.
    pub(crate) fn refine_once(&mut self) {
        let (lo, hi) = refine_once(&self.poly, &self.lo, &self.hi);
        self.lo = lo;
        self.hi = hi;
    }

    /// Shrinks the interval below `width` (a no-op for a point interval).
    pub(crate) fn refine_to(
        &mut self,
        width: &BigRational,
    ) {
        for _ in 0..MAX_REFINE {
            if self.lo == self.hi || &self.hi - &self.lo <= *width {
                return;
            }
            self.refine_once();
        }
    }

    /// The sign of `f` at this number, exactly. `None` only when the
    /// refinement budget is exhausted (never for a correctly constructed
    /// number).
    #[must_use]
    pub(crate) fn sign_of(
        &self,
        f: &[BigRational],
    ) -> Option<i32> {
        let mut copy = self.clone();
        copy.sign_of_mut(f)
    }

    /// The sign of `f` at this number, refining the stored interval in
    /// place so that repeated questions get cheaper.
    #[must_use]
    pub(crate) fn sign_of_mut(
        &mut self,
        f: &[BigRational],
    ) -> Option<i32> {
        if degree(f).is_none() {
            return Some(0);
        }
        let g = gcd(f, &self.poly);
        if degree(&g).is_some() && count_closed(&g, &self.lo, &self.hi) > 0 {
            return Some(0);
        }
        for _ in 0..MAX_REFINE {
            if self.lo == self.hi {
                return Some(sign_of(&eval(f, &self.lo)));
            }
            let (min, max) = eval_interval(f, &self.lo, &self.hi);
            if min.is_positive() {
                return Some(1);
            }
            if max.is_negative() {
                return Some(-1);
            }
            let (next_lo, next_hi) = refine_once(&self.poly, &self.lo, &self.hi);
            if next_lo == self.lo && next_hi == self.hi {
                return None;
            }
            self.lo = next_lo;
            self.hi = next_hi;
        }
        None
    }

    /// Exact comparison with another real algebraic number.
    #[must_use]
    pub(crate) fn compare(
        &self,
        other: &Self,
    ) -> Option<std::cmp::Ordering> {
        let mut a = self.clone();
        let mut b = other.clone();
        for _ in 0..MAX_REFINE {
            if a.hi < b.lo {
                return Some(std::cmp::Ordering::Less);
            }
            if b.hi < a.lo {
                return Some(std::cmp::Ordering::Greater);
            }
            let lo = a.lo.clone().max(b.lo.clone());
            let hi = a.hi.clone().min(b.hi.clone());
            let g = gcd(&a.poly, &b.poly);
            if degree(&g).is_some() && count_closed(&g, &lo, &hi) > 0 {
                return Some(std::cmp::Ordering::Equal);
            }
            if a.lo == a.hi && b.lo == b.hi {
                return Some(a.lo.cmp(&b.lo));
            }
            a.refine_once();
            b.refine_once();
        }
        None
    }

    /// The index of this root among the real roots of its defining
    /// polynomial, ascending.
    #[must_use]
    pub(crate) fn root_index(&self) -> usize {
        let chain = sturm_chain(&self.poly);
        if self.lo == self.hi {
            sturm_count(&chain, None, Some(&self.lo)).saturating_sub(1)
        } else {
            sturm_count(&chain, None, Some(&self.lo))
        }
    }

    /// The `root_of(k, ...)` term (or an exact rational / radical for degree
    /// one and two).
    #[must_use]
    pub(crate) fn to_term(
        &self,
        graph: &mut Graph,
    ) -> Option<NodeId> {
        match degree(&self.poly) {
            | None | Some(0) => None,
            | Some(1) => {
                let c0 = self.poly.first()?.clone();
                let c1 = self.poly.get(1)?.clone();
                Some(graph.num(Number::rat(-c0 / c1)))
            },
            | Some(2) => self.quadratic_term(graph),
            | Some(_) => {
                let index = self.root_index();
                let (integers, _) = integer_coefficients(&self.poly)?;
                let op = graph.ops().lookup("root_of")?;
                let mut args = vec![graph.int(i64::try_from(index).ok()?)];
                for c in integers {
                    args.push(graph.num(Number::Int(c)));
                }
                Some(graph.node(op, &args))
            },
        }
    }

    /// `(-b ± sqrt(b² - 4ac)) / (2a)` as a term.
    fn quadratic_term(
        &self,
        graph: &mut Graph,
    ) -> Option<NodeId> {
        let c0 = self.poly.first()?.clone();
        let c1 = self.poly.get(1)?.clone();
        let c2 = self.poly.get(2)?.clone();
        let discriminant = &c1 * &c1 - BigRational::from_integer(BigInt::from(4)) * &c2 * &c0;
        if discriminant.is_negative() {
            return None;
        }
        let chain = sturm_chain(&self.poly);
        let below = sturm_count(&chain, None, Some(&self.lo));
        let sign = if below == 0 { -1_i64 } else { 1 };
        let denominator = BigRational::from_integer(BigInt::from(2)) * &c2;
        let centre = graph.num(Number::rat(-c1 / denominator.clone()));
        let sqrt = rational_sqrt(&discriminant);
        if let Some(root) = sqrt {
            let candidate = (-self.poly.get(1)?.clone()
                + BigRational::from_integer(BigInt::from(sign)) * root)
                / denominator;
            return Some(graph.num(Number::rat(candidate)));
        }
        let half = Number::fraction(1, 2)?;
        let radicand = graph.num(Number::rat(discriminant));
        let exponential = graph.num(half);
        let radical = graph.node(core::POW, &[radicand, exponential]);
        let coefficient = graph.num(Number::rat(
            BigRational::from_integer(BigInt::from(sign)) / denominator,
        ));
        let offset = graph.node(core::MUL, &[coefficient, radical]);
        Some(graph.node(core::ADD, &[centre, offset]))
    }
}

/// `Some(r)` when the non-negative `value` is `r²` for a rational `r`.
fn rational_sqrt(value: &BigRational) -> Option<BigRational> {
    let num = value.numer().sqrt();
    let den = value.denom().sqrt();
    if &num * &num == *value.numer() && &den * &den == *value.denom() {
        Some(BigRational::new(num, den))
    } else {
        None
    }
}

/// Primitive integer coefficients with a positive leading coefficient.
fn integer_coefficients(p: &[BigRational]) -> Option<(Vec<BigInt>, BigInt)> {
    let n = degree(p)?;
    let mut denominator = BigInt::one();
    for c in p.iter().take(n.saturating_add(1)) {
        denominator = denominator.lcm(c.denom());
    }
    let mut integers: Vec<BigInt> = Vec::with_capacity(n.saturating_add(1));
    let mut content = BigInt::zero();
    for c in p.iter().take(n.saturating_add(1)) {
        let value = c.numer() * (&denominator / c.denom());
        content = content.gcd(&value);
        integers.push(value);
    }
    if content.is_zero() {
        return None;
    }
    let mut integers: Vec<BigInt> = integers.into_iter().map(|c| c / &content).collect();
    if integers.last().is_some_and(Signed::is_negative) {
        for c in &mut integers {
            *c = -c.clone();
        }
    }
    Some((integers, content))
}

#[cfg(test)]
mod tests {
    use crate::graph::Env;
    use crate::graph::Graph;

    use super::*;

    fn q(
        n: i64,
        d: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(n), BigInt::from(d))
    }

    fn poly(coeffs: &[i64]) -> Vec<BigRational> {
        coeffs
            .iter()
            .map(|&c| BigRational::from_integer(BigInt::from(c)))
            .collect()
    }

    #[test]
    fn division_and_gcd() {
        let a = poly(&[-1, 0, 1]);
        let b = poly(&[-1, 1]);
        let (q1, r) = divrem(&a, &b).unwrap_or_default();
        assert_eq!(q1, poly(&[1, 1]));
        assert!(degree(&r).is_none());
        assert_eq!(gcd(&a, &b), poly(&[-1, 1]));
        // The square-free part of (x - 1)^2 is x - 1.
        assert_eq!(square_free_part(&poly(&[1, -2, 1])), poly(&[-1, 1]));
        assert_eq!(square_free_part(&poly(&[0, 0, 1])), poly(&[0, 1]));
    }

    #[test]
    fn isolation_and_sturm() {
        let p = poly(&[-2, 0, 1]);
        let roots = isolate_roots(&p);
        assert_eq!(roots.len(), 2);
        let first = roots.first().cloned().unwrap_or_default();
        let second = roots.get(1).cloned().unwrap_or_default();
        assert!(first.0 < q(0, 1) && first.1 <= q(0, 1), "first {first:?}");
        assert!(second.0 >= q(0, 1), "second {second:?}");
        for interval in [&first, &second] {
            let at_lo = eval(&p, &interval.0);
            let at_hi = eval(&p, &interval.1);
            assert!(at_lo.is_negative() != at_hi.is_negative(), "{interval:?}");
        }
        let chain = sturm_chain(&p);
        assert_eq!(sturm_count(&chain, Some(&q(-2, 1)), Some(&q(2, 1))), 2);
        let r = poly(&[-1, 0, 1]);
        assert_eq!(count_closed(&r, &q(-1, 1), &q(1, 1)), 2);
        assert_eq!(count_open(&r, &q(-1, 1), &q(1, 1)), 0);
    }

    #[test]
    fn signed_square_roots() {
        let p = poly(&[-2, 0, 1]);
        let lo = RealAlg::from_root(&p, q(-2, 1), q(-1, 1)).unwrap_or_else(|| panic!());
        assert_eq!(lo.sign_of(&poly(&[-2, 0, 1])), Some(0));
        assert_eq!(lo.sign_of(&poly(&[1, 0])), Some(1));
        assert_eq!(lo.sign_of(&poly(&[0, 1])), Some(-1));
        let value = {
            let mut g = Graph::new();
            let node = lo.to_term(&mut g).unwrap_or_else(|| panic!());
            g.eval(node, &Env::numeric(0.0)).unwrap_or_else(|| panic!())
        };
        assert!((value + 2.0_f64.sqrt()).abs() < 1e-12);
        let hi = RealAlg::from_root(&p, q(1, 1), q(2, 1)).unwrap_or_else(|| panic!());
        assert_eq!(hi.sign_of(&poly(&[0, 1])), Some(1));
        assert_eq!(lo.compare(&hi), Some(std::cmp::Ordering::Less));
        let c = RealAlg::from_root(&p, q(1, 1), q(2, 1)).unwrap_or_else(|| panic!());
        assert_eq!(c.compare(&hi), Some(std::cmp::Ordering::Equal));
        assert_eq!(hi.root_index(), 1);
        assert_eq!(lo.root_index(), 0);
        let value = {
            let mut g = Graph::new();
            let node = hi.to_term(&mut g).unwrap_or_else(|| panic!());
            g.eval(node, &Env::numeric(0.0)).unwrap_or_else(|| panic!())
        };
        assert!((value - 2.0_f64.sqrt()).abs() < 1e-12);
    }
}
