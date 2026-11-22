//! Dense univariate polynomials over the rationals and their factorisation.
//!
//! Coefficient vectors are in ascending degree with no trailing zeros; the
//! empty vector is the zero polynomial. Factorisation over `Q` follows the
//! classical route: make the polynomial primitive over `Z`, split off
//! repeated factors (Yun), and factor each square-free part with
//! Berlekamp–Zassenhaus — factor modulo a small prime (distinct-degree and
//! Cantor–Zassenhaus equal-degree splitting), Hensel-lift the modular
//! factors past the Landau–Mignotte bound, and recombine them into true
//! factors by trial division.
//!
//! Parallelism: for polynomials of degree 24 and more the modular images (the first five admissible
//! primes) are factored concurrently, and with 12 or more modular factors the Zassenhaus
//! recombination tests consecutive subsets in batches of 128, taking the first hit in subset order,
//! so the factorisation is identical to the sequential one.

use std::sync::Arc;
use std::sync::OnceLock;

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_bigint::Sign;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::kernels::finite_field::BarrettMod;
use crate::kernels::finite_field::MontgomeryPrimeField;
use crate::kernels::finite_field::is_prime_u64;
use crate::kernels::par;
use crate::kernels::reuse::global_pool;

/// Coefficients over `Q`, ascending degree.
pub type QPoly = Vec<BigRational>;
/// Coefficients over `Z`, ascending degree.
pub type ZPoly = Vec<BigInt>;

fn trim<T: Zero>(mut p: Vec<T>) -> Vec<T> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

fn trim_slice<T: Zero>(mut p: &[T]) -> &[T] {
    while let Some(last) = p.last() {
        if last.is_zero() {
            p = &p[..p.len() - 1];
        } else {
            break;
        }
    }
    p
}

/// Degree, or `None` for the zero polynomial.
#[must_use]
pub const fn degree<T>(p: &[T]) -> Option<usize> {
    p.len().checked_sub(1)
}

/// Sum of two polynomials over `Q`.
#[must_use]
pub fn add(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let mut out: QPoly = a.to_vec();
    out.resize(a.len().max(b.len()), BigRational::zero());
    for (slot, x) in out.iter_mut().zip(b) {
        *slot += x;
    }
    trim(out)
}

/// `a - b` over `Q`.
#[must_use]
pub fn sub(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let negated: QPoly = b.iter().map(|x| -x).collect();
    add(a, &negated)
}

fn clear_denominators(p: &[BigRational]) -> (BigInt, ZPoly) {
    let mut den = BigInt::one();
    for c in p {
        let d = c.denom();
        if !d.is_one() {
            den = den.lcm(d);
        }
    }
    if den.is_one() {
        let z = p.iter().map(|c| c.numer().clone()).collect();
        return (den, z);
    }
    let z = p
        .iter()
        .map(|c| {
            if c.is_zero() {
                BigInt::zero()
            } else if c.denom().is_one() {
                c.numer() * &den
            } else {
                c.numer() * (&den / c.denom())
            }
        })
        .collect();
    (den, z)
}

/// Product over `Q`.
///
/// Clears common denominators when both operands have at least 4 coefficients,
/// reducing to fast integer polynomial multiplication [`z_mul`] and $O(n)$
/// final normalization.
#[must_use]
pub fn mul(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let a = trim_slice(a);
    let b = trim_slice(b);
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() == 1 {
        let s = &a[0];
        return trim(b.iter().map(|y| s * y).collect());
    }
    if b.len() == 1 {
        let s = &b[0];
        return trim(a.iter().map(|x| x * s).collect());
    }
    if a.len() < 4 || b.len() < 4 {
        let mut out = vec![BigRational::zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            if x.is_zero() {
                continue;
            }
            for (j, y) in b.iter().enumerate() {
                if !y.is_zero() {
                    out[i + j] += x * y;
                }
            }
        }
        return trim(out);
    }
    let (den_a, za) = clear_denominators(a);
    let (den_b, zb) = clear_denominators(b);
    let prod = z_mul(&za, &zb);
    let den_prod = den_a * den_b;
    if den_prod.is_one() {
        trim(prod.into_iter().map(BigRational::from_integer).collect())
    } else {
        trim(
            prod.into_iter()
                .map(|c| {
                    if c.is_zero() {
                        BigRational::zero()
                    } else {
                        BigRational::new(c, den_prod.clone())
                    }
                })
                .collect(),
        )
    }
}

fn divrem_schoolbook(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let a = trim_slice(a);
    let b = trim_slice(b);
    let lead = b.last()?;
    if a.len() < b.len() {
        return Some((Vec::new(), a.to_vec()));
    }
    let mut rem: QPoly = a.to_vec();
    let q_len = a.len() - b.len() + 1;
    let mut quo = vec![BigRational::zero(); q_len];
    for k in (0..q_len).rev() {
        let factor = &rem[k + b.len() - 1] / lead;
        if !factor.is_zero() {
            for (j, y) in b[..b.len() - 1].iter().enumerate() {
                if !y.is_zero() {
                    rem[k + j] -= &factor * y;
                }
            }
            rem[k + b.len() - 1] = BigRational::zero();
        }
        quo[k] = factor;
    }
    rem.truncate(b.len() - 1);
    Some((trim(quo), trim(rem)))
}

fn newton_series_inv(
    rev_b: &[BigRational],
    k: usize,
) -> QPoly {
    if k == 0 {
        return Vec::new();
    }
    let mut g = vec![BigRational::one()];
    let mut cur = 1_usize;
    while cur < k {
        let next = (cur * 2).min(k);
        let b_trunc = &rev_b[..rev_b.len().min(next)];
        let mut bg = mul(b_trunc, &g);
        bg.truncate(next);
        if bg.len() > cur {
            let err = &bg[cur..];
            let g_trunc = &g[..g.len().min(next - cur)];
            let mut delta = mul(g_trunc, err);
            delta.truncate(next - cur);
            g.resize(cur, BigRational::zero());
            for c in delta {
                g.push(-c);
            }
            g = trim(g);
        }
        cur = next;
    }
    g
}

fn divrem_newton(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let a = trim_slice(a);
    let b = trim_slice(b);
    let lead_b = b.last()?;
    if a.len() < b.len() {
        return Some((Vec::new(), a.to_vec()));
    }
    if b.len() == 1 {
        let quo = a.iter().map(|c| c / lead_b).collect();
        return Some((trim(quo), Vec::new()));
    }
    let k = a.len() - b.len() + 1;
    let b_monic = monic(b);
    let rev_b: QPoly = b_monic.into_iter().rev().collect();
    let inv_rev_b = newton_series_inv(&rev_b, k);
    let rev_a_trunc = trim(a.iter().rev().take(k).cloned().collect());
    let mut rev_q = mul(&rev_a_trunc, &inv_rev_b);
    rev_q.resize(k, BigRational::zero());
    let quo_monic: QPoly = rev_q.into_iter().rev().collect();
    let quo = if lead_b.is_one() {
        trim(quo_monic)
    } else {
        trim(quo_monic.into_iter().map(|c| c / lead_b).collect())
    };
    let prod = mul(&quo, b);
    let rem = sub(a, &prod);
    Some((quo, rem))
}

/// Fast quotient and remainder of `a` by `b` over `Q`, using Newton-Raphson
/// power-series inversion for high degrees and schoolbook division for small
/// degrees. Returns `None` when `b` is zero.
#[must_use]
pub fn fast_divrem(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let a = trim_slice(a);
    let b = trim_slice(b);
    if b.is_empty() {
        return None;
    }
    if a.len() < b.len() {
        return Some((Vec::new(), a.to_vec()));
    }
    let q_len = a.len() - b.len() + 1;
    if b.len() >= 16 && q_len >= 16 {
        divrem_newton(a, b)
    } else {
        divrem_schoolbook(a, b)
    }
}

/// Quotient and remainder of `a` by `b` over `Q`. `None` when `b` is zero.
#[must_use]
pub fn divrem(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    fast_divrem(a, b)
}

/// Divides by the leading coefficient.
#[must_use]
pub fn monic(p: &[BigRational]) -> QPoly {
    let p = trim_slice(p);
    match p.last() {
        | Some(lead) if lead.is_one() => p.to_vec(),
        | Some(lead) => p.iter().map(|c| c / lead).collect(),
        | None => Vec::new(),
    }
}

fn z_to_monic_q(p: &[BigInt]) -> QPoly {
    let p = trim_slice(p);
    let Some(lead) = p.last() else {
        return Vec::new();
    };
    if lead.is_one() {
        return p.iter().cloned().map(BigRational::from_integer).collect();
    }
    p.iter()
        .map(|c| {
            if c.is_zero() {
                BigRational::zero()
            } else {
                BigRational::new(c.clone(), lead.clone())
            }
        })
        .collect()
}

/// Monic greatest common divisor over `Q`.
///
/// Routes polynomials of degree $\ge 2$ through [`primitive`] and [`z_gcd`]
/// to eliminate rational Euclidean coefficient swell.
#[must_use]
pub fn gcd(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let a = trim_slice(a);
    let b = trim_slice(b);
    if a.is_empty() {
        return monic(b);
    }
    if b.is_empty() {
        return monic(a);
    }
    if a.len() == 1 || b.len() == 1 {
        return vec![BigRational::one()];
    }
    if a.len() >= 3 && b.len() >= 3 {
        let (_, za) = primitive(a);
        let (_, zb) = primitive(b);
        let gz = z_gcd(&za, &zb);
        return z_to_monic_q(&gz);
    }
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    while !y.is_empty() {
        let Some((_, r)) = divrem(&x, &y) else {
            break;
        };
        x = y;
        y = r;
    }
    monic(&x)
}

/// The Sylvester resultant over `Q`. Coefficients are in ascending degree.
/// Returns `None` for a zero polynomial, including an all-zero coefficient
/// list. Two nonzero constants have resultant one.
#[must_use]
pub(crate) fn resultant(
    p: &[BigRational],
    q: &[BigRational],
) -> Option<BigRational> {
    let (mut left, mut right) = (trim(p.to_vec()), trim(q.to_vec()));
    if left.is_empty() || right.is_empty() {
        return None;
    }

    let mut value = BigRational::one();
    while right.len() > 1 {
        let (m, n) = (degree(&left)?, degree(&right)?);
        if m % 2 == 1 && n % 2 == 1 {
            value = -value;
        }
        if m < n {
            std::mem::swap(&mut left, &mut right);
            continue;
        }
        let (_, remainder) = divrem(&left, &right)?;
        let Some(r) = degree(&remainder) else {
            return Some(BigRational::zero());
        };
        // At roots of right, left and remainder agree. Reducing left's
        // degree changes the power of right's leading coefficient by m - r.
        let lead = right.last()?.clone();
        value *= num_traits::pow(lead, m - r);
        (left, right) = (right, remainder);
    }
    Some(value * num_traits::pow(right.first()?.clone(), degree(&left)?))
}

/// Derivative.
#[must_use]
pub fn derivative(p: &[BigRational]) -> QPoly {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * BigRational::from_integer(BigInt::from(i)))
        .collect()
}

/// Value at `x` (Horner).
#[must_use]
pub fn eval(
    p: &[BigRational],
    x: &BigRational,
) -> BigRational {
    p.iter()
        .rev()
        .fold(BigRational::zero(), |acc, c| acc * x + c)
}

enum SubproductNode {
    Leaf {
        poly: QPoly,
    },
    Branch {
        poly: QPoly,
        left: Box<SubproductNode>,
        right: Box<SubproductNode>,
    },
}

impl SubproductNode {
    fn poly(&self) -> &[BigRational] {
        match self {
            | Self::Leaf { poly } | Self::Branch { poly, .. } => poly,
        }
    }
}

fn build_subproduct_tree(xs: &[BigRational]) -> SubproductNode {
    if xs.len() == 1 {
        return SubproductNode::Leaf {
            poly: vec![-xs[0].clone(), BigRational::one()],
        };
    }
    let mid = xs.len() / 2;
    let left = build_subproduct_tree(&xs[..mid]);
    let right = build_subproduct_tree(&xs[mid..]);
    let poly = mul(left.poly(), right.poly());
    SubproductNode::Branch {
        poly,
        left: Box::new(left),
        right: Box::new(right),
    }
}

fn eval_subproduct_tree(
    node: &SubproductNode,
    rem: &[BigRational],
    out: &mut Vec<BigRational>,
) {
    match node {
        | SubproductNode::Leaf { poly } => {
            if rem.is_empty() {
                out.push(BigRational::zero());
            } else if rem.len() == 1 {
                out.push(rem[0].clone());
            } else {
                // `poly = [-x_i, 1]`, so `rem mod (x - x_i) = rem(x_i)`.
                out.push(eval(rem, &(-&poly[0])));
            }
        },
        | SubproductNode::Branch { left, right, .. } => {
            let left_poly = left.poly();
            if rem.len() >= left_poly.len() {
                let left_rem = fast_divrem(rem, left_poly).map_or_else(Vec::new, |(_, r)| r);
                eval_subproduct_tree(left, &left_rem, out);
            } else {
                eval_subproduct_tree(left, rem, out);
            }
            let right_poly = right.poly();
            if rem.len() >= right_poly.len() {
                let right_rem = fast_divrem(rem, right_poly).map_or_else(Vec::new, |(_, r)| r);
                eval_subproduct_tree(right, &right_rem, out);
            } else {
                eval_subproduct_tree(right, rem, out);
            }
        },
    }
}

fn multipoint_eval_tree(
    p: &[BigRational],
    xs: &[BigRational],
) -> Vec<BigRational> {
    if xs.is_empty() {
        return Vec::new();
    }
    let tree = build_subproduct_tree(xs);
    let root_rem = if p.len() >= tree.poly().len() {
        fast_divrem(p, tree.poly()).map_or_else(Vec::new, |(_, r)| r)
    } else {
        p.to_vec()
    };
    let mut out = Vec::with_capacity(xs.len());
    eval_subproduct_tree(&tree, &root_rem, &mut out);
    out
}

/// Evaluates `p` at every point in `xs` using a balanced subproduct tree and
/// fast division for $O(M(n) \log n)$ complexity on large point sets and
/// degrees (Horner evaluation for small inputs).
#[must_use]
pub fn multipoint_eval(
    p: &[BigRational],
    xs: &[BigRational],
) -> Vec<BigRational> {
    if xs.is_empty() {
        return Vec::new();
    }
    let p = trim_slice(p);
    if p.is_empty() {
        return vec![BigRational::zero(); xs.len()];
    }
    if p.len() == 1 {
        return vec![p[0].clone(); xs.len()];
    }
    if xs.len() < 8 || p.len() < 8 {
        return xs.iter().map(|x| eval(p, x)).collect();
    }
    multipoint_eval_tree(p, xs)
}

/// Square-free decomposition (Yun): `p = c * prod f_i^i` with the `f_i`
/// monic, square-free and pairwise coprime. Returns `(f_i, i)` for the
/// non-trivial `f_i`.
#[must_use]
pub fn square_free(p: &[BigRational]) -> Vec<(QPoly, u32)> {
    let mut out = Vec::new();
    let f = monic(p);
    if degree(&f).is_none_or(|d| d == 0) {
        return out;
    }
    let df = derivative(&f);
    let g = gcd(&f, &df);
    let exact = |a: &[BigRational], b: &[BigRational]| divrem(a, b).map_or_default(|(q, _)| q);
    let mut c = exact(&f, &g);
    let mut d = sub(&exact(&df, &g), &derivative(&c));
    let mut multiplicity = 1_u32;
    while degree(&c).is_some_and(|deg| deg > 0) {
        let part = gcd(&c, &d);
        if degree(&part).is_some_and(|deg| deg > 0) {
            out.push((part.clone(), multiplicity));
        }
        c = exact(&c, &part);
        d = sub(&exact(&d, &part), &derivative(&c));
        multiplicity += 1;
    }
    out
}

/// Splits a rational polynomial into its content and primitive integer
/// part with positive leading coefficient: `p = content * primitive`.
#[must_use]
pub fn primitive(p: &[BigRational]) -> (BigRational, ZPoly) {
    let p = trim_slice(p);
    if p.is_empty() {
        return (BigRational::zero(), Vec::new());
    }
    let (denominators, integers) = clear_denominators(p);
    let mut content = integers.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    if content.is_zero() {
        return (BigRational::zero(), Vec::new());
    }
    if integers.last().is_some_and(Signed::is_negative) {
        content = -content;
    }
    let prim = integers.iter().map(|c| c / &content).collect();
    (BigRational::new(content, denominators), prim)
}

fn exact_div_schoolbook(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<ZPoly> {
    let lead = b.last()?;
    let mut rem: ZPoly = a.to_vec();
    let q_len = a.len() - b.len() + 1;
    let mut quo = vec![BigInt::zero(); q_len];
    for k in (0..q_len).rev() {
        let (factor, leftover) = rem[k + b.len() - 1].div_rem(lead);
        if !leftover.is_zero() {
            return None;
        }
        if !factor.is_zero() {
            for (j, y) in b[..b.len() - 1].iter().enumerate() {
                if !y.is_zero() {
                    rem[k + j] -= &factor * y;
                }
            }
            rem[k + b.len() - 1] = BigInt::zero();
        }
        quo[k] = factor;
    }
    rem[..b.len() - 1]
        .iter()
        .all(Zero::is_zero)
        .then(|| trim(quo))
}

fn newton_series_inv_monic_z(
    rev_b: &[BigInt],
    k: usize,
) -> ZPoly {
    if k == 0 {
        return Vec::new();
    }
    let mut g = vec![BigInt::one()];
    let mut cur = 1_usize;
    while cur < k {
        let next = (cur * 2).min(k);
        let b_trunc = &rev_b[..rev_b.len().min(next)];
        let mut bg = z_mul(b_trunc, &g);
        bg.truncate(next);
        if bg.len() > cur {
            let err = &bg[cur..];
            let g_trunc = &g[..g.len().min(next - cur)];
            let mut delta = z_mul(g_trunc, err);
            delta.truncate(next - cur);
            g.resize(cur, BigInt::zero());
            for c in delta {
                g.push(-c);
            }
            g = trim(g);
        }
        cur = next;
    }
    g
}

fn exact_div_monic_newton(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<ZPoly> {
    let lead_b = b.last()?;
    let flip = lead_b.is_negative();
    let k = a.len() - b.len() + 1;
    let rev_b: ZPoly = if flip {
        b.iter().rev().map(|c| -c).collect()
    } else {
        b.iter().rev().cloned().collect()
    };
    let inv_rev_b = newton_series_inv_monic_z(&rev_b, k);
    let rev_a_trunc = trim(if flip {
        a.iter().rev().take(k).map(|c| -c).collect()
    } else {
        a.iter().rev().take(k).cloned().collect()
    });
    let mut rev_q = z_mul(&rev_a_trunc, &inv_rev_b);
    rev_q.resize(k, BigInt::zero());
    let quo = trim(rev_q.into_iter().rev().collect());
    (z_mul(&quo, b) == a).then_some(quo)
}

/// Exact quotient `a / b` over `Z`, or `None` if `b` does not divide `a`.
///
/// Applies $O(1)/O(n)$ early rejection filters (degree, trailing zero
/// valuation, leading/trailing coefficient divisibility, evaluations at
/// $\pm 1$, and modular trial division modulo a 31-bit prime) before
/// computing the quotient.
#[must_use]
pub fn exact_div(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<ZPoly> {
    let a = trim_slice(a);
    let b = trim_slice(b);
    let lead_b = b.last()?;
    if a.is_empty() {
        return Some(Vec::new());
    }
    if a.len() < b.len() {
        return None;
    }
    if b.len() == 1 {
        let mut quo = Vec::with_capacity(a.len());
        for c in a {
            let (q, r) = c.div_rem(lead_b);
            if !r.is_zero() {
                return None;
            }
            quo.push(q);
        }
        return Some(trim(quo));
    }

    // 1. Trailing-zero (x^k) valuation filter.
    let tz_a = a.iter().take_while(|c| c.is_zero()).count();
    let tz_b = b.iter().take_while(|c| c.is_zero()).count();
    if tz_a < tz_b {
        return None;
    }
    let shift = tz_a - tz_b;
    let a = &a[tz_a..];
    let b = &b[tz_b..];
    if a.len() < b.len() {
        return None;
    }
    if b.len() == 1 {
        let mut quo = vec![BigInt::zero(); shift];
        for c in a {
            let (q, r) = c.div_rem(lead_b);
            if !r.is_zero() {
                return None;
            }
            quo.push(q);
        }
        return Some(trim(quo));
    }

    // 2. O(1) leading and trailing nonzero coefficient divisibility checks.
    let lead_a = &a[a.len() - 1];
    if lead_a.bits() < lead_b.bits() || !(lead_a % lead_b).is_zero() {
        return None;
    }
    let trail_a = &a[0];
    let trail_b = &b[0];
    if trail_a.bits() < trail_b.bits() || !(trail_a % trail_b).is_zero() {
        return None;
    }

    // 3. O(n) integer evaluation divisibility checks at x = 1 and x = -1.
    let mut a1 = BigInt::zero();
    let mut am1 = BigInt::zero();
    for (i, c) in a.iter().enumerate() {
        a1 += c;
        if i % 2 == 0 {
            am1 += c;
        } else {
            am1 -= c;
        }
    }
    let mut b1 = BigInt::zero();
    let mut bm1 = BigInt::zero();
    for (i, c) in b.iter().enumerate() {
        b1 += c;
        if i % 2 == 0 {
            bm1 += c;
        } else {
            bm1 -= c;
        }
    }
    if b1.is_zero() {
        if !a1.is_zero() {
            return None;
        }
    } else if !a1.is_zero() && (a1.bits() < b1.bits() || !(a1 % &b1).is_zero()) {
        return None;
    }
    if bm1.is_zero() {
        if !am1.is_zero() {
            return None;
        }
    } else if !am1.is_zero() && (am1.bits() < bm1.bits() || !(am1 % &bm1).is_zero()) {
        return None;
    }

    // 4. Modular F_p trial division check for moderate/large degrees.
    if b.len() >= 6 && a.len() >= 12 {
        const CHECK_PRIMES: [u64; 5] = [
            2_147_483_647,
            2_147_483_629,
            2_147_483_587,
            2_147_483_579,
            2_147_483_563,
        ];
        if let Some(&p) = CHECK_PRIMES.iter().find(|&&p| !(lead_b % p).is_zero()) {
            let fp = Fp::new(p);
            let ap = fp.reduce(a);
            let bp = fp.reduce(b);
            if !fp.rem(&ap, &bp).is_empty() {
                return None;
            }
        }
    }

    let q_len = a.len() - b.len() + 1;
    let quo = if b.len() >= 48 && q_len >= 48 && lead_b.magnitude().is_one() {
        exact_div_monic_newton(a, b)?
    } else {
        exact_div_schoolbook(a, b)?
    };

    if shift == 0 || quo.is_empty() {
        Some(quo)
    } else {
        let mut shifted = Vec::with_capacity(shift + quo.len());
        shifted.resize(shift, BigInt::zero());
        shifted.extend(quo);
        Some(shifted)
    }
}

// ----------------------------------------------------------------------
// High-speed integer polynomial multiplication:
// Schoolbook -> Karatsuba -> 3-Prime 62-bit Montgomery NTT + Garner CRT / Kronecker
// ----------------------------------------------------------------------

const KARATSUBA_Z_MIN_LEN: usize = 24;
const NTT_Z_MIN_LEN: usize = 96;
const NTT_PAR_MIN_LEN: usize = 512;
const KRONECKER_MAX_PACKED_BITS: u64 = 256_000_000;

// Three verified 62-bit NTT-friendly primes (< 2^62):
// P1 = 29 * 2^57 + 1 (g = 3), P2 = 27 * 2^56 + 1 (g = 5), P3 = 163 * 2^54 + 1 (g = 3).
const NTT_P1: u64 = 4_179_340_454_199_820_289;
const NTT_G1: u64 = 3;
const NTT_P2: u64 = 1_945_555_039_024_054_273;
const NTT_G2: u64 = 5;
const NTT_P3: u64 = 2_936_346_957_045_563_393;
const NTT_G3: u64 = 3;

const INV_P1_MOD_P2: u64 = 313_799_199_842_589_398;
const INV_P1P2_MOD_P3: u64 = 2_350_625_047_563_747_447;
const P1_P2_U128: u128 = (NTT_P1 as u128) * (NTT_P2 as u128);
const HALF_P1_P2_U128: u128 = P1_P2_U128 / 2;
const HALF_P3: u64 = NTT_P3 / 2;

fn half_p123() -> &'static BigUint {
    static HALF: OnceLock<BigUint> = OnceLock::new();
    HALF.get_or_init(|| {
        let prod = BigUint::from(NTT_P1) * BigUint::from(NTT_P2) * BigUint::from(NTT_P3);
        prod >> 1
    })
}

fn p123_biguint() -> &'static BigUint {
    static PROD: OnceLock<BigUint> = OnceLock::new();
    PROD.get_or_init(|| BigUint::from(NTT_P1) * BigUint::from(NTT_P2) * BigUint::from(NTT_P3))
}

fn p1p2_bigint() -> &'static BigInt {
    static P1P2: OnceLock<BigInt> = OnceLock::new();
    P1P2.get_or_init(|| BigInt::from(P1_P2_U128))
}

#[derive(Clone, Debug)]
struct NttTwiddles {
    mont: MontgomeryPrimeField,
    len: usize,
    roots_mont: Vec<u64>,
    inv_roots_mont: Vec<u64>,
    n_inv_mont: u64,
}

fn build_ntt_twiddles(
    prime: u64,
    prim_root: u64,
    len: usize,
) -> Option<NttTwiddles> {
    if len == 0 || !len.is_power_of_two() || prime <= 2 {
        return None;
    }
    let len_u64 = u64::try_from(len).ok()?;
    if !(prime - 1).is_multiple_of(len_u64) {
        return None;
    }
    let mont = MontgomeryPrimeField::new(prime)?;
    let g_mont = mont.to_mont(prim_root);
    let omega_mont = mont.pow_mont(g_mont, (prime - 1) / len_u64);
    let omega_inv = mont.inv(mont.from_mont(omega_mont))?;
    let omega_inv_mont = mont.to_mont(omega_inv);
    let half = (len / 2).max(1);
    let mut roots_mont = Vec::with_capacity(half);
    let mut inv_roots_mont = Vec::with_capacity(half);
    let mut cur = mont.one_mont();
    let mut cur_inv = mont.one_mont();
    for _ in 0..half {
        roots_mont.push(cur);
        inv_roots_mont.push(cur_inv);
        cur = mont.mul_mont(cur, omega_mont);
        cur_inv = mont.mul_mont(cur_inv, omega_inv_mont);
    }
    let n_inv = mont.inv(len_u64)?;
    let n_inv_mont = mont.to_mont(n_inv);
    Some(NttTwiddles {
        mont,
        len,
        roots_mont,
        inv_roots_mont,
        n_inv_mont,
    })
}

fn cached_ntt_twiddles(
    prime: u64,
    prim_root: u64,
    len: usize,
) -> Option<Arc<NttTwiddles>> {
    let cache = global_pool()
        .named_plan_cache::<(u64, usize), NttTwiddles>("univariate_ntt_62bit_twiddles", 64);
    cache
        .try_get_or_insert_with((prime, len), || {
            build_ntt_twiddles(prime, prim_root, len)
                .map(Arc::new)
                .ok_or(())
        })
        .ok()
}

#[cfg(target_pointer_width = "64")]
fn ntt_dif_mont_usize(
    values: &mut [usize],
    tw: &NttTwiddles,
) {
    let n = tw.len;
    if values.len() < n || n <= 1 {
        return;
    }
    let mont = tw.mont;
    let roots = &tw.roots_mont;
    let mut half = n / 2;
    let mut stride = 1_usize;
    while half >= 1 {
        let step = half * 2;
        let mut base = 0_usize;
        while base < n {
            for j in 0..half {
                let u = values[base + j] as u64;
                let v = values[base + j + half] as u64;
                let w = roots[j * stride];
                values[base + j] = mont.add_mont(u, v) as usize;
                values[base + j + half] = mont.mul_mont(mont.sub_mont(u, v), w) as usize;
            }
            base += step;
        }
        half /= 2;
        stride *= 2;
    }
}

#[cfg(target_pointer_width = "64")]
fn ntt_dit_inv_mont_usize(
    values: &mut [usize],
    tw: &NttTwiddles,
) {
    let n = tw.len;
    if values.len() < n || n <= 1 {
        return;
    }
    let mont = tw.mont;
    let inv_roots = &tw.inv_roots_mont;
    let mut half = 1_usize;
    let mut stride = n / 2;
    while half < n {
        let step = half * 2;
        let mut base = 0_usize;
        while base < n {
            for j in 0..half {
                let u = values[base + j] as u64;
                let w = inv_roots[j * stride];
                let v = mont.mul_mont(values[base + j + half] as u64, w);
                values[base + j] = mont.add_mont(u, v) as usize;
                values[base + j + half] = mont.sub_mont(u, v) as usize;
            }
            base += step;
        }
        half = step;
        stride /= 2;
    }
    let n_inv_mont = tw.n_inv_mont;
    for slot in &mut values[..n] {
        *slot = mont.mul_mont(*slot as u64, n_inv_mont) as usize;
    }
}

#[cfg(not(target_pointer_width = "64"))]
fn ntt_dif_mont_u64(
    values: &mut [u64],
    tw: &NttTwiddles,
) {
    let n = tw.len;
    if values.len() < n || n <= 1 {
        return;
    }
    let mont = tw.mont;
    let roots = &tw.roots_mont;
    let mut half = n / 2;
    let mut stride = 1_usize;
    while half >= 1 {
        let step = half * 2;
        let mut base = 0_usize;
        while base < n {
            for j in 0..half {
                let u = values[base + j];
                let v = values[base + j + half];
                let w = roots[j * stride];
                values[base + j] = mont.add_mont(u, v);
                values[base + j + half] = mont.mul_mont(mont.sub_mont(u, v), w);
            }
            base += step;
        }
        half /= 2;
        stride *= 2;
    }
}

#[cfg(not(target_pointer_width = "64"))]
fn ntt_dit_inv_mont_u64(
    values: &mut [u64],
    tw: &NttTwiddles,
) {
    let n = tw.len;
    if values.len() < n || n <= 1 {
        return;
    }
    let mont = tw.mont;
    let inv_roots = &tw.inv_roots_mont;
    let mut half = 1_usize;
    let mut stride = n / 2;
    while half < n {
        let step = half * 2;
        let mut base = 0_usize;
        while base < n {
            for j in 0..half {
                let u = values[base + j];
                let w = inv_roots[j * stride];
                let v = mont.mul_mont(values[base + j + half], w);
                values[base + j] = mont.add_mont(u, v);
                values[base + j + half] = mont.sub_mont(u, v);
            }
            base += step;
        }
        half = step;
        stride /= 2;
    }
    let n_inv_mont = tw.n_inv_mont;
    for slot in &mut values[..n] {
        *slot = mont.mul_mont(*slot, n_inv_mont);
    }
}

#[inline]
fn bigint_mod_u64(
    c: &BigInt,
    p: u64,
    barrett: BarrettMod,
    mont: MontgomeryPrimeField,
) -> u64 {
    if c.is_zero() {
        return 0;
    }
    if let Some(v) = c.to_i128() {
        let r = barrett.reduce_u128(v.unsigned_abs());
        return if v < 0 && r != 0 { p - r } else { r };
    }
    let mut rem = 0_u64;
    let r2_mod = mont.r2_mod();
    let mut pow_r_mont = mont.one_mont();
    for limb in c.magnitude().iter_u64_digits() {
        let limb_mod = barrett.reduce_u64(limb);
        let term = mont.redc(u128::from(limb_mod) * u128::from(pow_r_mont));
        rem = barrett.add(rem, term);
        pow_r_mont = mont.mul_mont(pow_r_mont, r2_mod);
    }
    if c.is_negative() && rem != 0 {
        p - rem
    } else {
        rem
    }
}

fn ntt_convolve_with<FA, FB>(
    a_len: usize,
    b_len: usize,
    prime: u64,
    prim_root: u64,
    n: usize,
    out_len: usize,
    is_square: bool,
    fill_a: FA,
    fill_b: FB,
) -> Option<Vec<u64>>
where
    FA: FnOnce(&mut [u64], MontgomeryPrimeField),
    FB: FnOnce(&mut [u64], MontgomeryPrimeField),
{
    let tw = cached_ntt_twiddles(prime, prim_root, n)?;
    let mont = tw.mont;

    #[cfg(target_pointer_width = "64")]
    {
        let pool = global_pool();
        let mut fa = pool.scratch_usize(n);
        fa.resize(n, 0);
        {
            // Fill into temporary u64 loop or directly via closure over &mut [u64]?
            // Convert in-place using a small adapter or write directly into usize slice.
            let mut tmp_a = vec![0_u64; a_len];
            fill_a(&mut tmp_a, mont);
            for (dst, &v) in fa[..a_len].iter_mut().zip(&tmp_a) {
                *dst = v as usize;
            }
        }
        ntt_dif_mont_usize(&mut fa, &tw);
        if is_square {
            for slot in &mut fa[..n] {
                let x = *slot as u64;
                *slot = mont.mul_mont(x, x) as usize;
            }
        } else {
            let mut fb = pool.scratch_usize(n);
            fb.resize(n, 0);
            {
                let mut tmp_b = vec![0_u64; b_len];
                fill_b(&mut tmp_b, mont);
                for (dst, &v) in fb[..b_len].iter_mut().zip(&tmp_b) {
                    *dst = v as usize;
                }
            }
            ntt_dif_mont_usize(&mut fb, &tw);
            for (da, &db) in fa[..n].iter_mut().zip(&fb[..n]) {
                *da = mont.mul_mont(*da as u64, db as u64) as usize;
            }
        }
        ntt_dit_inv_mont_usize(&mut fa, &tw);
        let out = fa[..out_len]
            .iter()
            .map(|&v| mont.from_mont(v as u64))
            .collect();
        Some(out)
    }

    #[cfg(not(target_pointer_width = "64"))]
    {
        let mut fa = vec![0_u64; n];
        fill_a(&mut fa[..a_len], mont);
        ntt_dif_mont_u64(&mut fa, &tw);
        if is_square {
            for slot in &mut fa[..n] {
                *slot = mont.mul_mont(*slot, *slot);
            }
        } else {
            let mut fb = vec![0_u64; n];
            fill_b(&mut fb[..b_len], mont);
            ntt_dif_mont_u64(&mut fb, &tw);
            for (da, &db) in fa[..n].iter_mut().zip(&fb[..n]) {
                *da = mont.mul_mont(*da, db);
            }
        }
        ntt_dit_inv_mont_u64(&mut fa, &tw);
        let out = fa[..out_len].iter().map(|&v| mont.from_mont(v)).collect();
        Some(out)
    }
}

fn z_mul_schoolbook(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() == 1 {
        let x = &a[0];
        return trim(b.iter().map(|y| x * y).collect());
    }
    if b.len() == 1 {
        let y = &b[0];
        return trim(a.iter().map(|x| x * y).collect());
    }
    let (short, long) = if a.len() <= b.len() {
        (a, b)
    } else {
        (b, a)
    };
    let mut out = vec![BigInt::zero(); short.len() + long.len() - 1];
    for (i, x) in short.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in long.iter().enumerate() {
            if !y.is_zero() {
                out[i + j] += x * y;
            }
        }
    }
    trim(out)
}

fn z_add_slices(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let mut out = a.to_vec();
    out.resize(a.len().max(b.len()), BigInt::zero());
    for (slot, y) in out.iter_mut().zip(b) {
        if !y.is_zero() {
            *slot += y;
        }
    }
    trim(out)
}

fn z_mul_karatsuba_rec(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len().min(b.len()) < KARATSUBA_Z_MIN_LEN {
        return z_mul_schoolbook(a, b);
    }
    if a.len() > b.len() {
        return z_mul_karatsuba_rec(b, a);
    }
    let na = a.len();
    let nb = b.len();
    let m = nb / 2;
    if na <= m {
        let (b0, b1) = b.split_at(m);
        let p0 = z_mul_karatsuba_rec(a, trim_slice(b0));
        let p1 = z_mul_karatsuba_rec(a, trim_slice(b1));
        let total_len = (na + nb - 1).max(p0.len()).max(m + p1.len());
        let mut out = vec![BigInt::zero(); total_len];
        for (slot, val) in out.iter_mut().zip(p0) {
            *slot += val;
        }
        for (slot, val) in out[m..].iter_mut().zip(p1) {
            *slot += val;
        }
        return trim(out);
    }

    let (a0, a1) = a.split_at(m);
    let (b0, b1) = b.split_at(m);
    let z0 = z_mul_karatsuba_rec(trim_slice(a0), trim_slice(b0));
    let z2 = z_mul_karatsuba_rec(trim_slice(a1), trim_slice(b1));
    let sum_a = z_add_slices(a0, a1);
    let sum_b = z_add_slices(b0, b1);
    let mut mid = z_mul_karatsuba_rec(&sum_a, &sum_b);

    let sub_len = z0.len().max(z2.len());
    if mid.len() < sub_len {
        mid.resize(sub_len, BigInt::zero());
    }
    for (slot, val) in mid.iter_mut().zip(&z0) {
        if !val.is_zero() {
            *slot -= val;
        }
    }
    for (slot, val) in mid.iter_mut().zip(&z2) {
        if !val.is_zero() {
            *slot -= val;
        }
    }
    let mid = trim(mid);

    let total_len = (na + nb - 1)
        .max(z0.len())
        .max(m + mid.len())
        .max(2 * m + z2.len());
    let mut out = vec![BigInt::zero(); total_len];
    for (slot, val) in out.iter_mut().zip(z0) {
        *slot += val;
    }
    for (slot, val) in out[m..].iter_mut().zip(mid) {
        *slot += val;
    }
    for (slot, val) in out[2 * m..].iter_mut().zip(z2) {
        *slot += val;
    }
    trim(out)
}

fn z_mul_ntt_crt(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<ZPoly> {
    let out_len = a.len().checked_add(b.len())?.checked_sub(1)?;
    let n = out_len.checked_next_power_of_two()?;
    let is_square = std::ptr::eq(a.as_ptr(), b.as_ptr()) && a.len() == b.len();
    let primes = [(NTT_P1, NTT_G1), (NTT_P2, NTT_G2), (NTT_P3, NTT_G3)];

    let run_one = |idx: usize| -> Option<Vec<u64>> {
        let (prime, prim_root) = primes[idx];
        let barrett = BarrettMod::new(prime)?;
        ntt_convolve_with(
            a.len(),
            b.len(),
            prime,
            prim_root,
            n,
            out_len,
            is_square,
            |dst, mont| {
                for (slot, c) in dst.iter_mut().zip(a) {
                    let r = bigint_mod_u64(c, prime, barrett, mont);
                    *slot = mont.to_mont(r);
                }
            },
            |dst, mont| {
                for (slot, c) in dst.iter_mut().zip(b) {
                    let r = bigint_mod_u64(c, prime, barrett, mont);
                    *slot = mont.to_mont(r);
                }
            },
        )
    };

    let residues: Vec<Option<Vec<u64>>> = if n >= NTT_PAR_MIN_LEN {
        par::map_indexed(3, 2, run_one)
    } else {
        (0..3).map(run_one).collect()
    };

    let mut iter = residues.into_iter();
    let r1s = iter.next()??;
    let r2s = iter.next()??;
    let r3s = iter.next()??;

    let barrett_p2 = BarrettMod::new(NTT_P2)?;
    let barrett_p3 = BarrettMod::new(NTT_P3)?;
    let mont_p2 = MontgomeryPrimeField::new(NTT_P2)?;
    let mont_p3 = MontgomeryPrimeField::new(NTT_P3)?;
    let inv12_r = mont_p2.to_mont(INV_P1_MOD_P2);
    let inv123_r = mont_p3.to_mont(INV_P1P2_MOD_P3);
    let p1p2_big = p1p2_bigint();

    let mut out = Vec::with_capacity(out_len);
    for i in 0..out_len {
        let r1 = r1s[i];
        let r2 = r2s[i];
        let r3 = r3s[i];

        let r1_mod_p2 = barrett_p2.reduce_u64(r1);
        let diff12 = barrett_p2.sub(r2, r1_mod_p2);
        let k2 = mont_p2.redc(u128::from(diff12) * u128::from(inv12_r));
        let x12 = u128::from(r1) + u128::from(k2) * u128::from(NTT_P1);

        let x12_mod_p3 = barrett_p3.reduce_u128(x12);
        let diff3 = barrett_p3.sub(r3, x12_mod_p3);
        let k3 = mont_p3.redc(u128::from(diff3) * u128::from(inv123_r));

        let val = if k3 == 0 {
            BigInt::from(x12)
        } else if k3 == NTT_P3 - 1 {
            -BigInt::from(P1_P2_U128 - x12)
        } else if k3 < HALF_P3 || (k3 == HALF_P3 && x12 <= HALF_P1_P2_U128) {
            BigInt::from(x12) + BigInt::from(k3) * p1p2_big
        } else {
            -(BigInt::from(NTT_P3 - k3) * p1p2_big - BigInt::from(x12))
        };
        out.push(val);
    }
    Some(trim(out))
}

fn pack_kronecker_limbs(
    slice: &[BigInt],
    words_per_coeff: usize,
    flip: bool,
) -> Vec<u32> {
    let mut limbs = vec![0_u32; slice.len() * words_per_coeff];
    let mut borrow = false;
    for (idx, c) in slice.iter().enumerate() {
        let dst = &mut limbs[idx * words_per_coeff..(idx + 1) * words_per_coeff];
        if c.is_zero() {
            if borrow {
                dst.fill(u32::MAX);
            }
            continue;
        }
        let is_neg = c.is_negative() ^ flip;
        if !is_neg {
            let mut b_val = u32::from(borrow);
            let mut digits = c.magnitude().iter_u32_digits();
            for slot in dst.iter_mut() {
                let d = digits.next().unwrap_or(0);
                let (res, underflow) = d.overflowing_sub(b_val);
                *slot = res;
                b_val = u32::from(underflow);
            }
            borrow = false;
        } else {
            let mut carry = u64::from(!borrow);
            let mut digits = c.magnitude().iter_u32_digits();
            for slot in dst.iter_mut() {
                let d = digits.next().unwrap_or(0);
                let sum = u64::from(!d) + carry;
                *slot = sum as u32;
                carry = sum >> 32;
            }
            borrow = true;
        }
    }
    limbs
}

fn unpack_kronecker_limbs(
    prod_digits: &[u32],
    out_len: usize,
    words_per_coeff: usize,
    flip_out: bool,
) -> ZPoly {
    let mut out = Vec::with_capacity(out_len);
    let mut borrow = false;
    let mut work = vec![0_u32; words_per_coeff];
    for k in 0..out_len {
        let start = k * words_per_coeff;
        let end = (start + words_per_coeff).min(prod_digits.len());
        work.fill(0);
        if start < end {
            work[..end - start].copy_from_slice(&prod_digits[start..end]);
        }
        let mut overflow = false;
        if borrow {
            let mut carry = 1_u32;
            for slot in &mut work {
                let (res, c) = slot.overflowing_add(carry);
                *slot = res;
                carry = u32::from(c);
                if carry == 0 {
                    break;
                }
            }
            if carry != 0 {
                overflow = true;
            }
        }
        if overflow {
            out.push(BigInt::zero());
            borrow = true;
        } else if (work[words_per_coeff - 1] >> 31) != 0 {
            let mut carry = 1_u64;
            for slot in &mut work {
                let sum = u64::from(!*slot) + carry;
                *slot = sum as u32;
                carry = sum >> 32;
            }
            let mag = BigUint::from_slice(&work);
            let sign = if mag.is_zero() {
                Sign::NoSign
            } else if flip_out {
                Sign::Plus
            } else {
                Sign::Minus
            };
            out.push(BigInt::from_biguint(sign, mag));
            borrow = true;
        } else {
            let mag = BigUint::from_slice(&work);
            let sign = if mag.is_zero() {
                Sign::NoSign
            } else if flip_out {
                Sign::Minus
            } else {
                Sign::Plus
            };
            out.push(BigInt::from_biguint(sign, mag));
            borrow = false;
        }
    }
    trim(out)
}

fn z_mul_kronecker(
    a: &[BigInt],
    b: &[BigInt],
    words_per_coeff: usize,
) -> ZPoly {
    let flip_a = a.last().is_some_and(Signed::is_negative);
    let flip_b = b.last().is_some_and(Signed::is_negative);
    let flip_out = flip_a ^ flip_b;

    let is_square = std::ptr::eq(a.as_ptr(), b.as_ptr()) && a.len() == b.len();
    let limbs_a = pack_kronecker_limbs(a, words_per_coeff, flip_a);
    let big_a = BigUint::new(limbs_a);
    let prod = if is_square {
        &big_a * &big_a
    } else {
        let limbs_b = pack_kronecker_limbs(b, words_per_coeff, flip_b);
        let big_b = BigUint::new(limbs_b);
        big_a * big_b
    };
    let prod_digits = prod.to_u32_digits();
    let out_len = a.len() + b.len() - 1;
    unpack_kronecker_limbs(&prod_digits, out_len, words_per_coeff, flip_out)
}

/// Product of two dense integer polynomials over `Z` (ascending degree).
///
/// Selects across three tiers by operand length and coefficient magnitude:
/// 1. Zero-skipping schoolbook (`min_len < 24`).
/// 2. Recursive Karatsuba (`24 <= min_len < 96`, or extreme coefficient bit-lengths).
/// 3. 3-prime 62-bit Montgomery NTT + signed Garner CRT (parallelized across the
///    3 primes when transform size `n >= 512`), or Kronecker bit-packing
///    substitution into two `BigInt`s when coefficients exceed the 183-bit CRT bound.
#[must_use]
pub fn z_mul(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let a = trim_slice(a);
    let b = trim_slice(b);
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let shift_a = a.iter().take_while(|c| c.is_zero()).count();
    let shift_b = b.iter().take_while(|c| c.is_zero()).count();
    let a = &a[shift_a..];
    let b = &b[shift_b..];
    let shift = shift_a + shift_b;
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }

    let min_len = a.len().min(b.len());
    let prod = if min_len < KARATSUBA_Z_MIN_LEN {
        z_mul_schoolbook(a, b)
    } else if min_len < NTT_Z_MIN_LEN {
        z_mul_karatsuba_rec(a, b)
    } else {
        let max_bits_a = a.iter().map(BigInt::bits).max().unwrap_or(0);
        let max_bits_b = b.iter().map(BigInt::bits).max().unwrap_or(0);
        let min_len_bits = u64::from(usize::BITS - min_len.leading_zeros());
        let approx_bound_bits = max_bits_a
            .saturating_add(max_bits_b)
            .saturating_add(min_len_bits);

        let fits_3_primes = if approx_bound_bits <= 182 {
            true
        } else if approx_bound_bits <= 186 {
            let max_a = a
                .iter()
                .map(BigInt::magnitude)
                .max()
                .cloned()
                .unwrap_or_default();
            let max_b = b
                .iter()
                .map(BigInt::magnitude)
                .max()
                .cloned()
                .unwrap_or_default();
            max_a * max_b * BigUint::from(min_len) < *half_p123()
        } else {
            false
        };

        if fits_3_primes {
            z_mul_ntt_crt(a, b).unwrap_or_else(|| z_mul_karatsuba_rec(a, b))
        } else {
            let slot_bits = approx_bound_bits.saturating_add(2);
            let words_per_coeff = usize::try_from(slot_bits.div_ceil(32))
                .unwrap_or(usize::MAX)
                .max(1);
            let out_len = a.len() + b.len() - 1;
            let total_packed_bits = slot_bits.saturating_mul((a.len() + b.len()) as u64);
            if out_len.checked_mul(words_per_coeff).is_some()
                && total_packed_bits <= KRONECKER_MAX_PACKED_BITS
            {
                z_mul_kronecker(a, b, words_per_coeff)
            } else {
                z_mul_karatsuba_rec(a, b)
            }
        }
    };

    if shift == 0 || prod.is_empty() {
        prod
    } else {
        let mut out = Vec::with_capacity(shift + prod.len());
        out.resize(shift, BigInt::zero());
        out.extend(prod);
        out
    }
}

// ----------------------------------------------------------------------
// Brown's modular GCD + CRT + Collins Subresultant PRS over Z[x]
// ----------------------------------------------------------------------

const GCD_PRIMES_62: [u64; 32] = [
    4_179_340_454_199_820_289_u64,
    1_945_555_039_024_054_273_u64,
    2_936_346_957_045_563_393_u64,
    4_611_686_018_427_387_847_u64,
    4_611_686_018_427_387_817_u64,
    4_611_686_018_427_387_787_u64,
    4_611_686_018_427_387_761_u64,
    4_611_686_018_427_387_751_u64,
    4_611_686_018_427_387_737_u64,
    4_611_686_018_427_387_733_u64,
    4_611_686_018_427_387_709_u64,
    4_611_686_018_427_387_701_u64,
    4_611_686_018_427_387_631_u64,
    4_611_686_018_427_387_617_u64,
    4_611_686_018_427_387_587_u64,
    4_611_686_018_427_387_461_u64,
    4_611_686_018_427_387_421_u64,
    4_611_686_018_427_387_409_u64,
    4_611_686_018_427_387_329_u64,
    4_611_686_018_427_387_323_u64,
    4_611_686_018_427_387_301_u64,
    4_611_686_018_427_387_271_u64,
    4_611_686_018_427_387_241_u64,
    4_611_686_018_427_387_139_u64,
    4_611_686_018_427_387_131_u64,
    4_611_686_018_427_387_127_u64,
    4_611_686_018_427_387_113_u64,
    4_611_686_018_427_387_091_u64,
    4_611_686_018_427_387_073_u64,
    4_611_686_018_427_386_981_u64,
    4_611_686_018_427_386_923_u64,
    4_611_686_018_427_386_911_u64,
];

fn next_62bit_prime(
    index: &mut usize,
    fallback_cursor: &mut u64,
) -> u64 {
    if let Some(&p) = GCD_PRIMES_62.get(*index) {
        *index += 1;
        return p;
    }
    *index += 1;
    loop {
        let cand = *fallback_cursor;
        *fallback_cursor = cand.saturating_sub(2);
        if is_prime_u64(cand) {
            return cand;
        }
    }
}

fn z_content_and_primitive(p: &[BigInt]) -> (BigInt, ZPoly) {
    let p = trim_slice(p);
    if p.is_empty() {
        return (BigInt::zero(), Vec::new());
    }
    let content = p.iter().fold(BigInt::zero(), |acc, x| acc.gcd(x));
    if content.is_zero() {
        return (BigInt::zero(), Vec::new());
    }
    let signed_content = if p.last().is_some_and(Signed::is_negative) {
        -content.clone()
    } else {
        content.clone()
    };
    let prim = if signed_content.is_one() {
        p.to_vec()
    } else {
        p.iter().map(|x| x / &signed_content).collect()
    };
    (content, prim)
}

fn z_pos_lead(p: &[BigInt]) -> ZPoly {
    let p = trim_slice(p);
    if p.last().is_some_and(Signed::is_negative) {
        p.iter().map(|c| -c).collect()
    } else {
        p.to_vec()
    }
}

fn z_pseudo_rem(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let a = trim_slice(a);
    let b = trim_slice(b);
    let Some(lead) = b.last() else {
        return a.to_vec();
    };
    if a.len() < b.len() {
        return a.to_vec();
    }
    let mut rem = a.to_vec();
    let q_len = a.len() - b.len() + 1;
    for k in (0..q_len).rev() {
        let top = rem[k + b.len() - 1].clone();
        for c in &mut rem {
            *c *= lead;
        }
        if !top.is_zero() {
            for (j, y) in b.iter().enumerate() {
                if !y.is_zero() {
                    rem[k + j] -= &top * y;
                }
            }
        }
    }
    trim(rem)
}

fn z_gcd_subresultant(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let (mut u, mut v) = if a.len() >= b.len() {
        (a.to_vec(), b.to_vec())
    } else {
        (b.to_vec(), a.to_vec())
    };
    let mut g = BigInt::one();
    let mut h = BigInt::one();
    while !v.is_empty() {
        let delta = u.len() - v.len();
        let r = z_pseudo_rem(&u, &v);
        if r.is_empty() {
            return z_content_and_primitive(&v).1;
        }
        if r.len() == 1 {
            return vec![BigInt::one()];
        }
        u = v;
        let denom = &g * num_traits::pow(h.clone(), delta);
        v = r.into_iter().map(|c| c / &denom).collect();
        g = u.last().cloned().unwrap_or_else(BigInt::one);
        if delta == 1 {
            h = g.clone();
        } else if delta > 1 {
            let num = num_traits::pow(g.clone(), delta);
            let den = num_traits::pow(h, delta - 1);
            h = num / den;
        }
    }
    vec![BigInt::one()]
}

fn crt_step_z(
    acc_residues: &mut ZPoly,
    modulus: &mut BigInt,
    new_residues: &[u64],
    p: u64,
) -> bool {
    let p_big = BigInt::from(p);
    if modulus.is_zero() {
        *acc_residues = new_residues.iter().map(|&c| BigInt::from(c)).collect();
        *modulus = p_big;
        return true;
    }
    let fp = Fp::new(p);
    let m_mod_p = fp.mont.map_or_else(
        || modulus.mod_floor(&p_big).to_u64().unwrap_or(0),
        |mont| bigint_mod_u64(modulus, p, fp.barrett, mont),
    );
    let m_inv = fp.inv(m_mod_p);
    let new_modulus = &*modulus * &p_big;
    let half_old = &*modulus >> 1;
    let half_new = &new_modulus >> 1;

    if acc_residues.len() < new_residues.len() {
        acc_residues.resize(new_residues.len(), BigInt::zero());
    }
    let mut changed = false;
    for i in 0..acc_residues.len() {
        let r1 = &acc_residues[i];
        let old_sym = {
            let pos = r1.mod_floor(modulus);
            if pos > half_old {
                pos - &*modulus
            } else {
                pos
            }
        };
        let r2 = new_residues.get(i).copied().unwrap_or(0);
        let r1_mod_p = fp.mont.map_or_else(
            || r1.mod_floor(&p_big).to_u64().unwrap_or(0),
            |mont| bigint_mod_u64(r1, p, fp.barrett, mont),
        );
        let diff = fp.barrett.sub(r2, r1_mod_p);
        let k = fp.barrett.mul(diff, m_inv);
        let updated = if k == 0 {
            r1.clone()
        } else {
            r1 + &*modulus * BigInt::from(k)
        };
        let new_sym = {
            let pos = updated.mod_floor(&new_modulus);
            if pos > half_new {
                pos - &new_modulus
            } else {
                pos
            }
        };
        if new_sym != old_sym {
            changed = true;
        }
        acc_residues[i] = updated;
    }
    *modulus = new_modulus;
    changed
}

fn centered_z_poly(
    residues: &[BigInt],
    modulus: &BigInt,
) -> ZPoly {
    let half = modulus >> 1;
    let out = residues
        .iter()
        .map(|r| {
            let pos = r.mod_floor(modulus);
            if pos > half {
                pos - modulus
            } else {
                pos
            }
        })
        .collect();
    trim(out)
}

fn mignotte_bound_bits_z(
    a: &[BigInt],
    b: &[BigInt],
    gamma: &BigInt,
) -> u64 {
    let norm_bits = |p: &[BigInt]| -> u64 {
        let max_b = p.iter().map(BigInt::bits).max().unwrap_or(1);
        let len_b = u64::from(usize::BITS - p.len().leading_zeros());
        max_b + len_b + 2
    };
    let deg = a.len().min(b.len()).saturating_sub(1) as u64;
    let min_norm = norm_bits(a).min(norm_bits(b));
    gamma.bits() + deg + min_norm + 4
}

/// Greatest common divisor over `Z[x]` with positive leading coefficient.
///
/// Extracts integer contents and common $x^k$ factors, computes the primitive
/// GCD via Brown's 62-bit modular GCD + Garner CRT + [`exact_div`] verification,
/// and falls back to Collins' subresultant PRS if modular bounds are exhausted.
#[must_use]
pub fn z_gcd(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let a = trim_slice(a);
    let b = trim_slice(b);
    if a.is_empty() {
        return z_pos_lead(b);
    }
    if b.is_empty() {
        return z_pos_lead(a);
    }

    let (ca, pa_full) = z_content_and_primitive(a);
    let (cb, pb_full) = z_content_and_primitive(b);
    let c_gcd = ca.gcd(&cb);

    let shift_a = pa_full.iter().take_while(|c| c.is_zero()).count();
    let shift_b = pb_full.iter().take_while(|c| c.is_zero()).count();
    let common_shift = shift_a.min(shift_b);
    let pa = &pa_full[shift_a..];
    let pb = &pb_full[shift_b..];

    let finish = |mut prim: ZPoly| -> ZPoly {
        if prim.is_empty() {
            prim.push(BigInt::one());
        }
        if common_shift > 0 {
            let mut shifted = vec![BigInt::zero(); common_shift];
            shifted.extend(prim);
            prim = shifted;
        }
        if !c_gcd.is_one() {
            for c in &mut prim {
                *c *= &c_gcd;
            }
        }
        trim(prim)
    };

    if pa.len() <= 1 || pb.len() <= 1 {
        return finish(vec![BigInt::one()]);
    }
    if pa == pb {
        return finish(pa.to_vec());
    }
    if pa.len() >= pb.len() && exact_div(pa, pb).is_some() {
        return finish(pb.to_vec());
    }
    if pb.len() > pa.len() && exact_div(pb, pa).is_some() {
        return finish(pa.to_vec());
    }

    let lc_a = &pa[pa.len() - 1];
    let lc_b = &pb[pb.len() - 1];
    let gamma = lc_a.gcd(lc_b);
    let bound_bits = mignotte_bound_bits_z(pa, pb, &gamma);

    let mut best_deg = pa.len().min(pb.len()) - 1;
    let mut crt_residues: ZPoly = Vec::new();
    let mut crt_modulus = BigInt::zero();
    let mut prime_idx = 0_usize;
    let mut fallback_cursor = 4_611_686_018_427_386_909_u64;
    let mut primes_used = 0_usize;
    let max_deg = pa.len().max(pb.len());
    let par_batch = if max_deg >= 32 && crate::execution::should_split_work(4, max_deg * max_deg) {
        4
    } else {
        1
    };

    for round in 0..32_usize {
        // Test the first prime alone so small-coefficient GCDs finish in one step.
        let batch_size = if round == 0 { 1 } else { par_batch };
        let mut primes = Vec::with_capacity(batch_size);
        while primes.len() < batch_size {
            let p = next_62bit_prime(&mut prime_idx, &mut fallback_cursor);
            let fp = Fp::new(p);
            let la_mod = fp
                .mont
                .map_or(0, |m| bigint_mod_u64(lc_a, p, fp.barrett, m));
            let lb_mod = fp
                .mont
                .map_or(0, |m| bigint_mod_u64(lc_b, p, fp.barrett, m));
            if la_mod != 0 && lb_mod != 0 {
                primes.push(p);
            }
        }

        let eval_prime = |&p: &u64| -> (u64, PPoly) {
            let fp = Fp::new(p);
            let ap = fp.reduce(pa);
            let bp = fp.reduce(pb);
            let gp = fp.gcd(&ap, &bp);
            let gamma_p = fp.mont.map_or_else(
                || gamma.mod_floor(&BigInt::from(p)).to_u64().unwrap_or(0),
                |m| bigint_mod_u64(&gamma, p, fp.barrett, m),
            );
            let scaled = gp.into_iter().map(|c| fp.barrett.mul(c, gamma_p)).collect();
            (p, scaled)
        };

        let images = if primes.len() > 1 {
            par::map_slice(&primes, 1, eval_prime)
        } else {
            primes.iter().map(eval_prime).collect()
        };

        for (p, gp) in images {
            let deg = gp.len().saturating_sub(1);
            if deg == 0 {
                return finish(vec![BigInt::one()]);
            }
            if deg > best_deg {
                continue;
            }
            if deg < best_deg {
                best_deg = deg;
                crt_residues.clear();
                crt_modulus = BigInt::zero();
                primes_used = 0;
            }
            let changed = crt_step_z(&mut crt_residues, &mut crt_modulus, &gp, p);
            primes_used += 1;

            if primes_used == 1 || !changed || crt_modulus.bits() > bound_bits {
                let candidate = centered_z_poly(&crt_residues, &crt_modulus);
                let (_, prim_cand) = z_content_and_primitive(&candidate);
                if prim_cand.len() == best_deg + 1
                    && exact_div(pa, &prim_cand).is_some()
                    && exact_div(pb, &prim_cand).is_some()
                {
                    return finish(prim_cand);
                }
                if crt_modulus.bits() > bound_bits + 128 {
                    crt_residues.clear();
                    crt_modulus = BigInt::zero();
                    primes_used = 0;
                }
            }
        }
    }

    finish(z_gcd_subresultant(pa, pb))
}

// ----------------------------------------------------------------------
// Arithmetic modulo a prime p (< 2^62) with Barrett/Montgomery + Half-GCD + NTT
// ----------------------------------------------------------------------

/// Polynomials over `F_p` (`2 <= p < 2^62`), ascending degree.
type PPoly = Vec<u64>;

const FP_KARATSUBA_MIN_LEN: usize = 32;
const FP_NTT_MIN_LEN: usize = 128;
const FP_HGCD_THRESHOLD: usize = 24;

#[derive(Copy, Clone)]
struct Fp {
    p: u64,
    barrett: BarrettMod,
    mont: Option<MontgomeryPrimeField>,
}

#[derive(Clone)]
struct PolyMatrixFp {
    m00: PPoly,
    m01: PPoly,
    m10: PPoly,
    m11: PPoly,
    neg_det: bool,
}

impl PolyMatrixFp {
    fn identity() -> Self {
        Self {
            m00: vec![1],
            m01: Vec::new(),
            m10: Vec::new(),
            m11: vec![1],
            neg_det: false,
        }
    }

    fn mul(
        &self,
        rhs: &Self,
        fp: Fp,
    ) -> Self {
        let m00 = fp.add(&fp.mul(&self.m00, &rhs.m00), &fp.mul(&self.m01, &rhs.m10));
        let m01 = fp.add(&fp.mul(&self.m00, &rhs.m01), &fp.mul(&self.m01, &rhs.m11));
        let m10 = fp.add(&fp.mul(&self.m10, &rhs.m00), &fp.mul(&self.m11, &rhs.m10));
        let m11 = fp.add(&fp.mul(&self.m10, &rhs.m01), &fp.mul(&self.m11, &rhs.m11));
        Self {
            m00,
            m01,
            m10,
            m11,
            neg_det: self.neg_det ^ rhs.neg_det,
        }
    }

    fn step_quotient(
        &mut self,
        q: &[u64],
        fp: Fp,
    ) {
        let next_m00 = fp.add(&fp.mul(&self.m00, q), &self.m01);
        let next_m10 = fp.add(&fp.mul(&self.m10, q), &self.m11);
        self.m01 = std::mem::replace(&mut self.m00, next_m00);
        self.m11 = std::mem::replace(&mut self.m10, next_m10);
        self.neg_det = !self.neg_det;
    }

    fn apply_inverse(
        &self,
        a: &[u64],
        b: &[u64],
        fp: Fp,
    ) -> (PPoly, PPoly) {
        let t0 = fp.sub(&fp.mul(&self.m11, a), &fp.mul(&self.m01, b));
        let t1 = fp.sub(&fp.mul(&self.m00, b), &fp.mul(&self.m10, a));
        if self.neg_det {
            (fp.sub(&[], &t0), fp.sub(&[], &t1))
        } else {
            (t0, t1)
        }
    }
}

impl Fp {
    fn new(p: u64) -> Self {
        let valid_p = p.max(2);
        let Some(barrett) = BarrettMod::new(valid_p) else {
            unreachable!()
        };
        Self {
            p: valid_p,
            barrett,
            mont: MontgomeryPrimeField::new(valid_p),
        }
    }

    #[inline]
    const fn inv(
        self,
        x: u64,
    ) -> u64 {
        match self.barrett.inv(x) {
            | Some(v) => v,
            | None => 0,
        }
    }

    fn reduce(
        self,
        p: &[BigInt],
    ) -> PPoly {
        if let Some(mont) = self.mont {
            return trim(
                p.iter()
                    .map(|c| bigint_mod_u64(c, self.p, self.barrett, mont))
                    .collect(),
            );
        }
        let modulus = BigInt::from(self.p);
        trim(
            p.iter()
                .map(|c| c.mod_floor(&modulus).to_u64().unwrap_or(0))
                .collect(),
        )
    }

    fn add(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let mut out = a.to_vec();
        out.resize(a.len().max(b.len()), 0);
        for (slot, &y) in out.iter_mut().zip(b) {
            *slot = self.barrett.add(*slot, y);
        }
        trim(out)
    }

    fn sub(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let mut out = a.to_vec();
        out.resize(a.len().max(b.len()), 0);
        for (slot, &y) in out.iter_mut().zip(b) {
            *slot = self.barrett.sub(*slot, y);
        }
        trim(out)
    }

    fn mul_schoolbook(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let out_len = a.len() + b.len() - 1;
        if self.p < (1_u64 << 31) {
            let mut acc = vec![0_u128; out_len];
            let (short, long) = if a.len() <= b.len() {
                (a, b)
            } else {
                (b, a)
            };
            for (i, &x) in short.iter().enumerate() {
                if x == 0 {
                    continue;
                }
                let x128 = u128::from(x);
                for (j, &y) in long.iter().enumerate() {
                    if y != 0 {
                        acc[i + j] += x128 * u128::from(y);
                    }
                }
            }
            let out = acc
                .into_iter()
                .map(|v| self.barrett.reduce_u128(v))
                .collect();
            return trim(out);
        }
        let mut out = vec![0_u64; out_len];
        let (short, long) = if a.len() <= b.len() {
            (a, b)
        } else {
            (b, a)
        };
        for (i, &x) in short.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in long.iter().enumerate() {
                if y != 0 {
                    let prod = self.barrett.mul(x, y);
                    out[i + j] = self.barrett.add(out[i + j], prod);
                }
            }
        }
        trim(out)
    }

    fn mul_karatsuba(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let a = trim_slice(a);
        let b = trim_slice(b);
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        if a.len().min(b.len()) < FP_KARATSUBA_MIN_LEN {
            return self.mul_schoolbook(a, b);
        }
        if a.len() > b.len() {
            return self.mul_karatsuba(b, a);
        }
        let m = b.len() / 2;
        if a.len() <= m {
            let (b0, b1) = b.split_at(m);
            let p0 = self.mul_karatsuba(a, trim_slice(b0));
            let p1 = self.mul_karatsuba(a, trim_slice(b1));
            let total = (a.len() + b.len() - 1).max(p0.len()).max(m + p1.len());
            let mut out = vec![0_u64; total];
            for (slot, &v) in out.iter_mut().zip(&p0) {
                *slot = self.barrett.add(*slot, v);
            }
            for (slot, &v) in out[m..].iter_mut().zip(&p1) {
                *slot = self.barrett.add(*slot, v);
            }
            return trim(out);
        }

        let (a0, a1) = a.split_at(m);
        let (b0, b1) = b.split_at(m);
        let z0 = self.mul_karatsuba(trim_slice(a0), trim_slice(b0));
        let z2 = self.mul_karatsuba(trim_slice(a1), trim_slice(b1));
        let sum_a = self.add(a0, a1);
        let sum_b = self.add(b0, b1);
        let z1_full = self.mul_karatsuba(&sum_a, &sum_b);
        let mid = self.sub(&self.sub(&z1_full, &z0), &z2);

        let total = (a.len() + b.len() - 1)
            .max(z0.len())
            .max(m + mid.len())
            .max(2 * m + z2.len());
        let mut out = vec![0_u64; total];
        for (slot, &v) in out.iter_mut().zip(&z0) {
            *slot = self.barrett.add(*slot, v);
        }
        for (slot, &v) in out[m..].iter_mut().zip(&mid) {
            *slot = self.barrett.add(*slot, v);
        }
        for (slot, &v) in out[2 * m..].iter_mut().zip(&z2) {
            *slot = self.barrett.add(*slot, v);
        }
        trim(out)
    }

    fn mul_ntt(
        self,
        a: &[u64],
        b: &[u64],
    ) -> Option<PPoly> {
        let out_len = a.len().checked_add(b.len())?.checked_sub(1)?;
        let n = out_len.checked_next_power_of_two()?;
        let min_len = a.len().min(b.len());

        if self.mont.is_some()
            && let Ok(n_u64) = u64::try_from(n)
            && (self.p - 1).is_multiple_of(n_u64)
            && let Some(mut prod) = crate::kernels::finite_field::poly_mul_ntt_mod(a, b, self.p)
        {
            return Some(trim(std::mem::take(&mut prod)));
        }

        let is_square = std::ptr::eq(a.as_ptr(), b.as_ptr()) && a.len() == b.len();

        // 2-prime NTT for p < 2^31 (max coeff < min_len * 2^62 < P1 * P2).
        if self.p < (1_u64 << 31) && min_len < (1_usize << 58) {
            let run_one = |prime: u64, prim_root: u64| -> Option<Vec<u64>> {
                let barrett = BarrettMod::new(prime)?;
                ntt_convolve_with(
                    a.len(),
                    b.len(),
                    prime,
                    prim_root,
                    n,
                    out_len,
                    is_square,
                    |dst, mont| {
                        for (slot, &c) in dst.iter_mut().zip(a) {
                            *slot = mont.to_mont(barrett.reduce_u64(c));
                        }
                    },
                    |dst, mont| {
                        for (slot, &c) in dst.iter_mut().zip(b) {
                            *slot = mont.to_mont(barrett.reduce_u64(c));
                        }
                    },
                )
            };
            let r1s = run_one(NTT_P1, NTT_G1)?;
            let r2s = run_one(NTT_P2, NTT_G2)?;
            let barrett_p2 = BarrettMod::new(NTT_P2)?;
            let mont_p2 = MontgomeryPrimeField::new(NTT_P2)?;
            let inv12_r = mont_p2.to_mont(INV_P1_MOD_P2);
            let mut out = Vec::with_capacity(out_len);
            for (&r1, &r2) in r1s.iter().zip(&r2s) {
                let diff = barrett_p2.sub(r2, barrett_p2.reduce_u64(r1));
                let k2 = mont_p2.redc(u128::from(diff) * u128::from(inv12_r));
                let x12 = u128::from(r1) + u128::from(k2) * u128::from(NTT_P1);
                out.push(self.barrett.reduce_u128(x12));
            }
            return Some(trim(out));
        }

        // 3-prime NTT for larger p < 2^62 when unsigned sum fits in P1*P2*P3.
        if self.p < (1_u64 << 62) {
            let max_coeff = BigUint::from(self.p - 1);
            let max_prod = &max_coeff * &max_coeff * BigUint::from(min_len);
            if max_prod >= *p123_biguint() {
                return None;
            }
            let primes = [(NTT_P1, NTT_G1), (NTT_P2, NTT_G2), (NTT_P3, NTT_G3)];
            let run_one = |idx: usize| -> Option<Vec<u64>> {
                let (prime, prim_root) = primes[idx];
                let barrett = BarrettMod::new(prime)?;
                ntt_convolve_with(
                    a.len(),
                    b.len(),
                    prime,
                    prim_root,
                    n,
                    out_len,
                    is_square,
                    |dst, mont| {
                        for (slot, &c) in dst.iter_mut().zip(a) {
                            *slot = mont.to_mont(barrett.reduce_u64(c));
                        }
                    },
                    |dst, mont| {
                        for (slot, &c) in dst.iter_mut().zip(b) {
                            *slot = mont.to_mont(barrett.reduce_u64(c));
                        }
                    },
                )
            };
            let residues: Vec<Option<Vec<u64>>> = if n >= NTT_PAR_MIN_LEN {
                par::map_indexed(3, 2, run_one)
            } else {
                (0..3).map(run_one).collect()
            };
            let mut iter = residues.into_iter();
            let r1s = iter.next()??;
            let r2s = iter.next()??;
            let r3s = iter.next()??;

            let barrett_p2 = BarrettMod::new(NTT_P2)?;
            let barrett_p3 = BarrettMod::new(NTT_P3)?;
            let mont_p2 = MontgomeryPrimeField::new(NTT_P2)?;
            let mont_p3 = MontgomeryPrimeField::new(NTT_P3)?;
            let inv12_r = mont_p2.to_mont(INV_P1_MOD_P2);
            let inv123_r = mont_p3.to_mont(INV_P1P2_MOD_P3);
            let p1p2_mod_p = self.barrett.reduce_u128(P1_P2_U128);
            let mut out = Vec::with_capacity(out_len);
            for i in 0..out_len {
                let r1 = r1s[i];
                let r2 = r2s[i];
                let r3 = r3s[i];
                let diff12 = barrett_p2.sub(r2, barrett_p2.reduce_u64(r1));
                let k2 = mont_p2.redc(u128::from(diff12) * u128::from(inv12_r));
                let x12 = u128::from(r1) + u128::from(k2) * u128::from(NTT_P1);
                let diff3 = barrett_p3.sub(r3, barrett_p3.reduce_u128(x12));
                let k3 = mont_p3.redc(u128::from(diff3) * u128::from(inv123_r));
                let x12_mod_p = self.barrett.reduce_u128(x12);
                let hi_mod_p = self.barrett.mul(k3, p1p2_mod_p);
                out.push(self.barrett.add(x12_mod_p, hi_mod_p));
            }
            return Some(trim(out));
        }

        None
    }

    fn mul(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let a = trim_slice(a);
        let b = trim_slice(b);
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let min_len = a.len().min(b.len());
        if min_len < FP_KARATSUBA_MIN_LEN {
            self.mul_schoolbook(a, b)
        } else if min_len < FP_NTT_MIN_LEN {
            self.mul_karatsuba(a, b)
        } else {
            self.mul_ntt(a, b)
                .unwrap_or_else(|| self.mul_karatsuba(a, b))
        }
    }

    fn divrem_schoolbook(
        self,
        a: &[u64],
        b: &[u64],
    ) -> (PPoly, PPoly) {
        let Some(&lead) = b.last() else {
            return (Vec::new(), a.to_vec());
        };
        if a.len() < b.len() {
            return (Vec::new(), a.to_vec());
        }
        let inv = self.inv(lead);
        let mut rem = a.to_vec();
        let q_len = a.len() - b.len() + 1;
        let mut quo = vec![0_u64; q_len];
        for k in (0..q_len).rev() {
            let top = rem[k + b.len() - 1];
            if top == 0 {
                continue;
            }
            let factor = self.barrett.mul(top, inv);
            quo[k] = factor;
            for (j, &y) in b[..b.len() - 1].iter().enumerate() {
                if y != 0 {
                    let prod = self.barrett.mul(factor, y);
                    rem[k + j] = self.barrett.sub(rem[k + j], prod);
                }
            }
            rem[k + b.len() - 1] = 0;
        }
        rem.truncate(b.len() - 1);
        (trim(quo), trim(rem))
    }

    fn divrem_newton(
        self,
        a: &[u64],
        b: &[u64],
    ) -> (PPoly, PPoly) {
        let lead = b[b.len() - 1];
        let inv_lead = self.inv(lead);
        let k = a.len() - b.len() + 1;
        let b_monic: PPoly = b.iter().map(|&c| self.barrett.mul(c, inv_lead)).collect();
        let rev_b: PPoly = b_monic.into_iter().rev().collect();

        let mut g = vec![1_u64];
        let mut cur = 1_usize;
        while cur < k {
            let next = (cur * 2).min(k);
            let b_trunc = &rev_b[..rev_b.len().min(next)];
            let mut bg = self.mul(b_trunc, &g);
            bg.truncate(next);
            if bg.len() > cur {
                let err = &bg[cur..];
                let g_trunc = &g[..g.len().min(next - cur)];
                let mut delta = self.mul(g_trunc, err);
                delta.truncate(next - cur);
                g.resize(cur, 0);
                for &c in &delta {
                    g.push(self.barrett.neg(c));
                }
                g = trim(g);
            }
            cur = next;
        }

        let rev_a_trunc = trim(a.iter().rev().take(k).copied().collect());
        let mut rev_q = self.mul(&rev_a_trunc, &g);
        rev_q.resize(k, 0);
        let quo = trim(
            rev_q
                .into_iter()
                .rev()
                .map(|c| self.barrett.mul(c, inv_lead))
                .collect(),
        );
        let rem = self.sub(a, &self.mul(&quo, b));
        (quo, rem)
    }

    /// Quotient and remainder; `b` must be non-zero.
    fn divrem(
        self,
        a: &[u64],
        b: &[u64],
    ) -> (PPoly, PPoly) {
        let a = trim_slice(a);
        let b = trim_slice(b);
        let Some(&lead) = b.last() else {
            return (Vec::new(), a.to_vec());
        };
        if a.len() < b.len() {
            return (Vec::new(), a.to_vec());
        }
        if b.len() == 1 {
            let inv = self.inv(lead);
            let quo = trim(a.iter().map(|&c| self.barrett.mul(c, inv)).collect());
            return (quo, Vec::new());
        }
        let q_len = a.len() - b.len() + 1;
        if b.len() >= 32 && q_len >= 16 {
            self.divrem_newton(a, b)
        } else {
            self.divrem_schoolbook(a, b)
        }
    }

    fn rem(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        self.divrem(a, b).1
    }

    fn monic(
        self,
        p: &[u64],
    ) -> PPoly {
        let p = trim_slice(p);
        match p.last() {
            | Some(&lead) if lead == 1 => p.to_vec(),
            | Some(&lead) => {
                let inv = self.inv(lead);
                p.iter().map(|&c| self.barrett.mul(c, inv)).collect()
            },
            | None => Vec::new(),
        }
    }

    fn hgcd(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PolyMatrixFp {
        let deg_a = a.len().saturating_sub(1);
        let m = deg_a / 2;
        if b.is_empty() || b.len() - 1 <= m {
            return PolyMatrixFp::identity();
        }
        if deg_a < FP_HGCD_THRESHOLD {
            let mut mat = PolyMatrixFp::identity();
            let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
            while !r1.is_empty() && r1.len() - 1 > m {
                let (q, r2) = self.divrem(&r0, &r1);
                mat.step_quotient(&q, self);
                r0 = r1;
                r1 = r2;
            }
            return mat;
        }

        let a_hi = &a[m..];
        let b_hi = if b.len() > m { &b[m..] } else { &[] };
        let mut r_mat = self.hgcd(a_hi, b_hi);
        let (mut a1, mut b1) = r_mat.apply_inverse(a, b, self);
        if b1.is_empty() || b1.len() - 1 <= m {
            return r_mat;
        }

        let (q, r) = self.divrem(&a1, &b1);
        r_mat.step_quotient(&q, self);
        a1 = b1;
        b1 = r;
        if b1.is_empty() || b1.len() - 1 <= m {
            return r_mat;
        }

        let deg_a1 = a1.len() - 1;
        let k = (2 * m).saturating_sub(deg_a1);
        let a1_hi = if a1.len() > k {
            &a1[k..]
        } else {
            &[]
        };
        let b1_hi = if b1.len() > k {
            &b1[k..]
        } else {
            &[]
        };
        let s_mat = self.hgcd(a1_hi, b1_hi);
        r_mat.mul(&s_mat, self)
    }

    fn gcd(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let (mut x, mut y) = (trim(a.to_vec()), trim(b.to_vec()));
        if x.len() < y.len() {
            std::mem::swap(&mut x, &mut y);
        }
        while !y.is_empty() {
            if x.len() <= FP_HGCD_THRESHOLD || x.len() - y.len() > x.len() / 2 {
                let r = self.rem(&x, &y);
                x = y;
                y = r;
                continue;
            }
            if x.len() == y.len() {
                let r = self.rem(&x, &y);
                x = y;
                y = r;
                continue;
            }
            let mat = self.hgcd(&x, &y);
            let (nx, ny) = mat.apply_inverse(&x, &y, self);
            x = nx;
            y = ny;
            if !y.is_empty() {
                let r = self.rem(&x, &y);
                x = y;
                y = r;
            }
        }
        self.monic(&x)
    }

    /// `(s, t)` with `s*a + t*b = 1`, for coprime `a`, `b`.
    fn bezout(
        self,
        a: &[u64],
        b: &[u64],
    ) -> Option<(PPoly, PPoly)> {
        let (mut r0, mut r1) = (trim(a.to_vec()), trim(b.to_vec()));
        let mut total = PolyMatrixFp::identity();
        while !r1.is_empty() {
            if r0.len() > r1.len()
                && r0.len() > FP_HGCD_THRESHOLD
                && r0.len() - r1.len() <= r0.len() / 2
            {
                let m = self.hgcd(&r0, &r1);
                let (nr0, nr1) = m.apply_inverse(&r0, &r1, self);
                total = total.mul(&m, self);
                r0 = nr0;
                r1 = nr1;
                if r1.is_empty() {
                    break;
                }
            }
            let (q, r2) = self.divrem(&r0, &r1);
            total.step_quotient(&q, self);
            r0 = r1;
            r1 = r2;
        }
        let &[unit] = r0.as_slice() else {
            return None;
        };
        let inv = self.inv(unit);
        let (s_raw, t_raw) = if total.neg_det {
            (self.sub(&[], &total.m11), total.m01)
        } else {
            (total.m11, self.sub(&[], &total.m01))
        };
        let scale = |poly: PPoly| -> PPoly {
            trim(poly.into_iter().map(|c| self.barrett.mul(c, inv)).collect())
        };
        let mut s = scale(s_raw);
        let mut t = scale(t_raw);
        // Ensure canonical minimal-degree Bezout cofactors: deg(t) < deg(a).
        if !a.is_empty() && t.len() >= a.len() {
            let (q, r) = self.divrem(&t, a);
            t = r;
            s = self.add(&s, &self.mul(&q, b));
        }
        Some((s, t))
    }

    fn derivative(
        self,
        p: &[u64],
    ) -> PPoly {
        trim(
            p.iter()
                .enumerate()
                .skip(1)
                .map(|(i, &c)| self.barrett.mul(c, self.barrett.reduce_u64(i as u64)))
                .collect(),
        )
    }

    /// `base^exp mod modulus`.
    fn pow_mod(
        self,
        base: &[u64],
        exp: &BigUint,
        modulus: &[u64],
    ) -> PPoly {
        let mut acc: PPoly = vec![1];
        let base = self.rem(base, modulus);
        for bit in (0..exp.bits()).rev() {
            acc = self.rem(&self.mul(&acc, &acc), modulus);
            if exp.bit(bit) {
                acc = self.rem(&self.mul(&acc, &base), modulus);
            }
        }
        acc
    }

    /// Splits a monic square-free polynomial into products of irreducibles
    /// of equal degree: returns `(product, degree)` pairs.
    fn distinct_degree(
        self,
        f: &[u64],
    ) -> Vec<(PPoly, usize)> {
        let mut out = Vec::new();
        let mut f = f.to_vec();
        let x: PPoly = vec![0, 1];
        let mut h = x.clone();
        let mut d = 0_usize;
        let p = BigUint::from(self.p);
        while degree(&f).is_some_and(|deg| deg >= 2 * (d + 1)) {
            d += 1;
            h = self.pow_mod(&h, &p, &f);
            let g = self.gcd(&self.sub(&h, &x), &f);
            if degree(&g).is_some_and(|deg| deg > 0) {
                f = self.divrem(&f, &g).0;
                h = self.rem(&h, &f);
                out.push((g, d));
            }
        }
        if let Some(deg) = degree(&f).filter(|&deg| deg > 0) {
            out.push((f, deg));
        }
        out
    }

    /// Splits a product of irreducibles of degree `d` (Cantor–Zassenhaus,
    /// odd `p`).
    fn equal_degree(
        self,
        f: &[u64],
        d: usize,
        seed: &mut u64,
        out: &mut Vec<PPoly>,
    ) {
        let n = degree(f).unwrap_or(0);
        if n <= d {
            if n > 0 {
                out.push(self.monic(f));
            }
            return;
        }
        let exponent =
            (BigUint::from(self.p).pow(u32::try_from(d).unwrap_or(u32::MAX)) - 1_u32) >> 1;
        loop {
            let trial: PPoly = trim(
                (0..n)
                    .map(|_| {
                        *seed = seed
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        self.barrett.reduce_u64(*seed >> 33)
                    })
                    .collect(),
            );
            if degree(&trial).is_none_or(|deg| deg == 0) {
                continue;
            }
            let shared = self.gcd(&trial, f);
            let split = if degree(&shared).is_some_and(|deg| deg > 0) {
                shared
            } else {
                let power = self.pow_mod(&trial, &exponent, f);
                self.gcd(&self.sub(&power, &[1]), f)
            };
            if degree(&split).is_some_and(|deg| deg > 0 && deg < n) {
                let other = self.divrem(f, &split).0;
                self.equal_degree(&split, d, seed, out);
                self.equal_degree(&other, d, seed, out);
                return;
            }
        }
    }

    /// Monic irreducible factors of a monic square-free polynomial.
    fn factor(
        self,
        f: &[u64],
    ) -> Vec<PPoly> {
        let mut out = Vec::new();
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64 ^ self.p;
        for (product, d) in self.distinct_degree(f) {
            self.equal_degree(&product, d, &mut seed, &mut out);
        }
        out
    }
}

// ----------------------------------------------------------------------
// Quadratic Tree Hensel Lifting and Pruned Recombination
// ----------------------------------------------------------------------

fn symmetric(
    c: &BigInt,
    modulus: &BigInt,
) -> BigInt {
    let r = c.mod_floor(modulus);
    if &r * 2 > *modulus {
        r - modulus
    } else {
        r
    }
}

fn lift_coeffs(p: &[u64]) -> ZPoly {
    p.iter().map(|&c| BigInt::from(c)).collect()
}

fn z_add_mod(
    a: &[BigInt],
    b: &[BigInt],
    m: &BigInt,
) -> ZPoly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let av = a.get(i).map_or_else(BigInt::zero, Clone::clone);
        let bv = b.get(i).map_or_else(BigInt::zero, Clone::clone);
        out.push((av + bv).mod_floor(m));
    }
    trim(out)
}

fn z_sub_mod(
    a: &[BigInt],
    b: &[BigInt],
    m: &BigInt,
) -> ZPoly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let av = a.get(i).map_or_else(BigInt::zero, Clone::clone);
        let bv = b.get(i).map_or_else(BigInt::zero, Clone::clone);
        out.push((av - bv).mod_floor(m));
    }
    trim(out)
}

/// Divides `a` by monic `g` in `(Z / m Z)[x]`, performing only $O(n)$ `mod_floor` reductions.
fn z_divrem_monic_mod(
    a: &[BigInt],
    g: &[BigInt],
    m: &BigInt,
) -> (ZPoly, ZPoly) {
    let a = trim_slice(a);
    if a.len() < g.len() {
        let rem = trim(a.iter().map(|c| c.mod_floor(m)).collect());
        return (Vec::new(), rem);
    }
    let mut rem = a.to_vec();
    let q_len = a.len() - g.len() + 1;
    let mut quo = vec![BigInt::zero(); q_len];
    for k in (0..q_len).rev() {
        let factor = rem[k + g.len() - 1].mod_floor(m);
        if !factor.is_zero() {
            for (j, y) in g[..g.len() - 1].iter().enumerate() {
                if !y.is_zero() {
                    rem[k + j] -= &factor * y;
                }
            }
            rem[k + g.len() - 1] = BigInt::zero();
        }
        quo[k] = factor;
    }
    rem.truncate(g.len() - 1);
    for slot in &mut rem {
        *slot = slot.mod_floor(m);
    }
    (trim(quo), trim(rem))
}

/// Computes the quadratic lifting modulus ladder `p^{e_0}, ..., p^K >= target`
/// where `e_{i+1} <= 2 * e_i`.
fn hensel_modulus_ladder(
    p: u64,
    target: &BigInt,
) -> Vec<BigInt> {
    let p_big = BigInt::from(p);
    let mut k = 1_u32;
    let mut m = p_big.clone();
    while m < *target {
        m *= &p_big;
        k += 1;
    }
    let mut exps = Vec::new();
    let mut cur = k;
    while cur > 1 {
        exps.push(cur);
        cur = cur.div_ceil(2);
    }
    exps.push(1);
    exps.reverse();
    exps.into_iter().map(|e| p_big.pow(e)).collect()
}

/// Lifts `f = g*h (mod p)` quadratically along `ladder = [p, ..., p^K]`,
/// keeping `g` monic and `h` carrying the leading coefficient of `f`.
fn hensel_pair_ladder(
    fp: Fp,
    f: &[BigInt],
    g: &[u64],
    h: &[u64],
    ladder: &[BigInt],
) -> Option<(ZPoly, ZPoly)> {
    let (s_p, t_p) = fp.bezout(g, h)?;
    let mut big_g = lift_coeffs(g);
    let mut big_h = lift_coeffs(h);
    let mut s = lift_coeffs(&s_p);
    let mut t = lift_coeffs(&t_p);
    let f_lead = f.last()?.clone();
    if let Some(slot) = big_h.last_mut() {
        *slot = f_lead.mod_floor(&ladder[0]);
    }
    let steps = ladder.len();
    for step_idx in 1..steps {
        let next_m = &ladder[step_idx];
        if let Some(slot) = big_h.last_mut() {
            *slot = f_lead.mod_floor(next_m);
        }
        let gh = z_mul(&big_g, &big_h);
        let e = z_sub_mod(f, &gh, next_m);
        if !e.is_empty() {
            let te = z_mul(&t, &e);
            let (_, dg) = z_divrem_monic_mod(&te, &big_g, next_m);
            let rhs = z_sub_mod(&e, &z_mul(&big_h, &dg), next_m);
            let (dh, rem) = z_divrem_monic_mod(&rhs, &big_g, next_m);
            if !rem.is_empty() {
                return None;
            }
            big_g = z_add_mod(&big_g, &dg, next_m);
            big_h = z_add_mod(&big_h, &dh, next_m);
        }
        if step_idx + 1 < steps {
            let sg_th = z_add_mod(&z_mul(&s, &big_g), &z_mul(&t, &big_h), next_m);
            let e_bez = z_sub_mod(&[BigInt::one()], &sg_th, next_m);
            if !e_bez.is_empty() {
                let (_, dt) = z_divrem_monic_mod(&z_mul(&t, &e_bez), &big_g, next_m);
                let rhs_b = z_sub_mod(&e_bez, &z_mul(&big_h, &dt), next_m);
                let (ds, rem_b) = z_divrem_monic_mod(&rhs_b, &big_g, next_m);
                if !rem_b.is_empty() {
                    return None;
                }
                s = z_add_mod(&s, &ds, next_m);
                t = z_add_mod(&t, &dt, next_m);
            }
        }
    }
    Some((big_g, big_h))
}

fn fp_product_balanced(
    fp: Fp,
    factors: &[PPoly],
) -> PPoly {
    match factors {
        | [] => vec![1],
        | [only] => only.clone(),
        | [a, b] => fp.mul(a, b),
        | _ => {
            let mid = factors.len() / 2;
            let left = fp_product_balanced(fp, &factors[..mid]);
            let right = fp_product_balanced(fp, &factors[mid..]);
            fp.mul(&left, &right)
        },
    }
}

fn hensel_tree(
    fp: Fp,
    f_cur: &[BigInt],
    factors: &[PPoly],
    ladder: &[BigInt],
    final_modulus: &BigInt,
) -> Option<Vec<ZPoly>> {
    if factors.is_empty() {
        return Some(Vec::new());
    }
    if factors.len() == 1 {
        let lead = f_cur.last()?.clone();
        let lifted = if lead.is_one() {
            f_cur.iter().map(|c| symmetric(c, final_modulus)).collect()
        } else {
            let inverse = lead.extended_gcd(final_modulus).x.mod_floor(final_modulus);
            f_cur
                .iter()
                .map(|c| symmetric(&(c * &inverse), final_modulus))
                .collect()
        };
        return Some(vec![lifted]);
    }
    let mid = factors.len() / 2;
    let (left_f, right_f) = factors.split_at(mid);
    let g_p = fp_product_balanced(fp, left_f);
    let right_prod = fp_product_balanced(fp, right_f);
    let lead_p = fp.reduce(&[f_cur.last()?.clone()]);
    let h_p = fp.mul(&lead_p, &right_prod);
    let (big_g, big_h) = hensel_pair_ladder(fp, f_cur, &g_p, &h_p, ladder)?;
    let mut out = hensel_tree(fp, &big_g, left_f, ladder, final_modulus)?;
    let right_lifted = hensel_tree(fp, &big_h, right_f, ladder, final_modulus)?;
    out.extend(right_lifted);
    Some(out)
}

/// Lifts all modular factors of `f` to monic factors modulo a power of `p`
/// that is at least `target`. Returns the lifted factors and that modulus.
fn hensel_all(
    fp: Fp,
    f: &[BigInt],
    factors: &[PPoly],
    target: &BigInt,
) -> Option<(Vec<ZPoly>, BigInt)> {
    let ladder = hensel_modulus_ladder(fp.p, target);
    let modulus = ladder.last()?.clone();
    let lifted = hensel_tree(fp, f, factors, &ladder, &modulus)?;
    Some((lifted, modulus))
}

fn z_primitive(p: &[BigInt]) -> ZPoly {
    let p = trim_slice(p);
    let mut content = p.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    if content.is_zero() {
        return Vec::new();
    }
    if p.last().is_some_and(Signed::is_negative) {
        content = -content;
    }
    p.iter().map(|c| c / &content).collect()
}

fn z_product_mod_balanced(
    pool: &[ZPoly],
    subset: &[usize],
    modulus: &BigInt,
) -> ZPoly {
    match subset {
        | [] => vec![BigInt::one()],
        | &[idx] => pool[idx].clone(),
        | &[a, b] => {
            z_mul(&pool[a], &pool[b])
                .into_iter()
                .map(|c| symmetric(&c, modulus))
                .collect()
        },
        | _ => {
            let mid = subset.len() / 2;
            let left = z_product_mod_balanced(pool, &subset[..mid], modulus);
            let right = z_product_mod_balanced(pool, &subset[mid..], modulus);
            z_mul(&left, &right)
                .into_iter()
                .map(|c| symmetric(&c, modulus))
                .collect()
        },
    }
}

/// Degree from which the modular images are factored in parallel.
const MODULAR_PAR_DEGREE: usize = 24;
/// Smallest number of modular factors for which recombination tests subsets in parallel.
const RECOMBINE_PAR_POOL: usize = 12;
/// Consecutive subsets tested together by the parallel recombination.
const RECOMBINE_BATCH: usize = 128;

const PRIMES: [u64; 24] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

/// Irreducible factors over `Z` of a primitive square-free polynomial with
/// positive leading coefficient.
fn factor_square_free(f: &[BigInt]) -> Vec<ZPoly> {
    let n = degree(f).unwrap_or(0);
    if n <= 1 {
        return if n == 1 {
            vec![f.to_vec()]
        } else {
            Vec::new()
        };
    }
    let Some(lead) = f.last() else {
        return Vec::new();
    };
    // Choose the prime giving the fewest modular factors among a few
    // admissible ones: fewer factors, fewer recombinations.
    let mut best: Option<(Fp, Vec<PPoly>)> = None;
    let admissible = |p: u64| -> bool {
        let fp = Fp::new(p);
        if (lead % p).is_zero() {
            return false;
        }
        let reduced = fp.reduce(f);
        let shared = fp.gcd(&reduced, &fp.derivative(&reduced));
        degree(&shared) == Some(0)
    };
    if n >= MODULAR_PAR_DEGREE && crate::execution::effective_threads() > 1 {
        // factor modulo the first five admissible primes in parallel; the
        // earliest minimum wins, exactly as in the sequential scan below
        let chosen: Vec<u64> = PRIMES
            .iter()
            .copied()
            .filter(|&p| admissible(p))
            .take(5)
            .collect();
        let images = crate::kernels::par::map_slice(&chosen, 2, |&p| {
            let fp = Fp::new(p);
            (fp, fp.factor(&fp.monic(&fp.reduce(f))))
        });
        for (fp, factors) in images {
            if best
                .as_ref()
                .is_none_or(|(_, old)| factors.len() < old.len())
            {
                best = Some((fp, factors));
            }
        }
    } else {
        let mut tried = 0;
        for &p in &PRIMES {
            if !admissible(p) {
                continue;
            }
            let fp = Fp::new(p);
            let factors = fp.factor(&fp.monic(&fp.reduce(f)));
            if best
                .as_ref()
                .is_none_or(|(_, old)| factors.len() < old.len())
            {
                best = Some((fp, factors));
            }
            tried += 1;
            if tried == 5 || best.as_ref().is_some_and(|(_, fs)| fs.len() == 1) {
                break;
            }
        }
    }
    let Some((fp, modular)) = best else {
        return vec![f.to_vec()];
    };
    if modular.len() == 1 {
        return vec![f.to_vec()];
    }
    // Landau–Mignotte: a factor's coefficients are bounded by 2^n * |f|_2;
    // candidates are scaled by the leading coefficient.
    let norm_squared: BigInt = f.iter().map(|c| c * c).sum();
    let bound = (norm_squared.sqrt() + 1) * (BigInt::one() << n) * lead.abs() * 2;
    let Some((lifted, modulus)) = hensel_all(fp, f, &modular, &bound) else {
        return vec![f.to_vec()];
    };

    let half_mod = &modulus >> 1;
    let mut out = Vec::new();
    let mut work: ZPoly = f.to_vec();
    let mut pool: Vec<ZPoly> = lifted;
    let mut size = 1;
    while 2 * size <= pool.len() {
        let mut found = None;
        let mut choice: Vec<usize> = (0..size).collect();
        let work_lead = work.last().cloned().unwrap_or_else(BigInt::one);
        let scaled: ZPoly = work.iter().map(|c| c * &work_lead).collect();
        let scaled_const = scaled.first().cloned().unwrap_or_else(BigInt::zero);
        // The first successful subset in lexicographic order wins. Large
        // pools test a batch of consecutive subsets in parallel and take the
        // first hit of the batch in order, so the factor found is the same.
        let batch = if pool.len() >= RECOMBINE_PAR_POOL && crate::execution::effective_threads() > 1
        {
            RECOMBINE_BATCH
        } else {
            1
        };
        let mut exhausted = false;
        while !exhausted && found.is_none() {
            let mut subsets: Vec<Vec<usize>> = Vec::with_capacity(batch);
            while subsets.len() < batch {
                let mut viable = true;
                if !scaled_const.is_zero() {
                    let mut c0 = work_lead.clone();
                    for &idx in &choice {
                        if let Some(p0) = pool[idx].first() {
                            c0 = (&c0 * p0).mod_floor(&modulus);
                        } else {
                            c0 = BigInt::zero();
                            break;
                        }
                    }
                    if c0 > half_mod {
                        c0 -= &modulus;
                    }
                    if c0.is_zero() || !(&scaled_const % &c0).is_zero() {
                        viable = false;
                    }
                }
                if viable {
                    subsets.push(choice.clone());
                }
                // Next subset of `size` indices in lexicographic order.
                let mut i = size;
                let more = loop {
                    if i == 0 {
                        break false;
                    }
                    i -= 1;
                    if choice[i] != i + pool.len() - size {
                        break true;
                    }
                };
                if !more {
                    exhausted = true;
                    break;
                }
                choice[i] += 1;
                for j in i + 1..size {
                    choice[j] = choice[j - 1] + 1;
                }
            }
            if subsets.is_empty() {
                continue;
            }
            let test = |subset: &Vec<usize>| -> Option<ZPoly> {
                let prod_mod = z_product_mod_balanced(&pool, subset, &modulus);
                let candidate_scaled: ZPoly = prod_mod
                    .iter()
                    .map(|c| symmetric(&(c * &work_lead), &modulus))
                    .collect();
                let candidate = z_primitive(&candidate_scaled);
                (!candidate.is_empty() && exact_div(&scaled, &candidate).is_some())
                    .then_some(candidate)
            };
            let results = if subsets.len() >= 2 {
                crate::kernels::par::map_slice(&subsets, 2, test)
            } else {
                subsets.iter().map(test).collect()
            };
            if let Some((k, candidate)) = results
                .into_iter()
                .enumerate()
                .find_map(|(k, r)| r.map(|c| (k, c)))
            {
                found = Some((subsets[k].clone(), candidate));
            }
        }
        match found {
            | Some((choice, factor)) => {
                if let Some(quotient) = exact_div(&work, &factor) {
                    work = quotient;
                } else {
                    // The candidate divides lead*work; the primitive part of
                    // the quotient is what remains.
                    let lead_now = work.last().cloned().unwrap_or_else(BigInt::one);
                    let scaled: ZPoly = work.iter().map(|c| c * &lead_now).collect();
                    work = z_primitive(&exact_div(&scaled, &factor).unwrap_or_default());
                }
                out.push(factor);
                for &i in choice.iter().rev() {
                    pool.remove(i);
                }
            },
            | None => size += 1,
        }
    }
    if degree(&work).is_some_and(|d| d > 0) {
        out.push(z_primitive(&work));
    }
    out
}

/// Complete factorisation over `Q`.
///
/// `p = content * prod factor_i^e_i` with
/// each factor primitive over `Z`, irreducible, and with positive leading
/// coefficient. Factors are ordered by degree, then by coefficients.
#[must_use]
pub fn factor(p: &[BigRational]) -> (BigRational, Vec<(ZPoly, u32)>) {
    let (content, _) = primitive(p);
    let mut out = Vec::new();
    for (part, multiplicity) in square_free(p) {
        let (_, integer_part) = primitive(&part);
        for factor in factor_square_free(&integer_part) {
            out.push((factor, multiplicity));
        }
    }
    out.sort_by(|(a, ea), (b, eb)| {
        a.len()
            .cmp(&b.len())
            .then_with(|| a.cmp(b))
            .then_with(|| ea.cmp(eb))
    });
    // The content of the input relative to the product of primitive
    // factors: compare leading coefficients.
    let lead_of_product = out.iter().fold(BigInt::one(), |acc, (f, e)| {
        acc * num_traits::pow(f.last().cloned().unwrap_or_else(BigInt::one), *e as usize)
    });
    let unit = match p.last() {
        | Some(lead) => lead / BigRational::from_integer(lead_of_product),
        | None => content,
    };
    (unit, out)
}

/// The rational roots of `p`, each once, ascending.
#[must_use]
pub fn rational_roots(p: &[BigRational]) -> Vec<BigRational> {
    let mut roots: Vec<BigRational> = factor(p)
        .1
        .into_iter()
        .filter_map(|(f, _)| {
            match f.as_slice() {
                | [c0, c1] => Some(BigRational::new(-c0.clone(), c1.clone())),
                | _ => None,
            }
        })
        .collect();
    roots.sort();
    roots.dedup();
    roots
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(coeffs: &[i64]) -> QPoly {
        coeffs
            .iter()
            .map(|&c| BigRational::from_integer(BigInt::from(c)))
            .collect()
    }

    fn z(coeffs: &[i64]) -> ZPoly {
        coeffs.iter().map(|&c| BigInt::from(c)).collect()
    }

    #[test]
    fn resultant_preserves_orientation_and_non_monic_scaling() {
        let third = BigRational::new(BigInt::one(), BigInt::from(3));
        let five_thirds = BigRational::new(BigInt::from(5), BigInt::from(3));
        let p = [third, BigRational::one()];
        let r = q(&[-1, 2]);
        assert_eq!(resultant(&p, &r), Some(-&five_thirds));
        assert_eq!(resultant(&r, &p), Some(five_thirds));
        assert_eq!(
            resultant(&q(&[-5, 13, 6]), &q(&[3, 4])),
            Some(BigRational::from_integer(BigInt::from(-182)))
        );
        assert_eq!(
            resultant(&q(&[-5, 13, 6]), &q(&[-1, 3])),
            Some(BigRational::zero())
        );
    }

    #[test]
    fn resultant_distinguishes_constants_zero_and_trailing_padding() {
        assert_eq!(resultant(&q(&[3, 0]), &q(&[5])), Some(BigRational::one()));
        assert_eq!(
            resultant(&q(&[3]), &q(&[-2, 0, 0, 1, 0])),
            Some(BigRational::from_integer(BigInt::from(27)))
        );
        assert_eq!(
            resultant(&q(&[-2, 0, 0, 1]), &q(&[-3])),
            Some(BigRational::from_integer(BigInt::from(-27)))
        );
        assert_eq!(resultant(&q(&[]), &q(&[1, 1])), None);
        assert_eq!(resultant(&q(&[1, 1]), &q(&[0, 0])), None);
    }

    fn product(
        unit: &BigRational,
        factors: &[(ZPoly, u32)],
    ) -> QPoly {
        let mut out: QPoly = vec![unit.clone()];
        for (f, e) in factors {
            let fq: QPoly = f
                .iter()
                .map(|c| BigRational::from_integer(c.clone()))
                .collect();
            for _ in 0..*e {
                out = mul(&out, &fq);
            }
        }
        out
    }

    #[test]
    fn division_and_gcd() {
        // (x^2 - 1) / (x - 1) = x + 1
        assert_eq!(
            divrem(&q(&[-1, 0, 1]), &q(&[-1, 1])),
            Some((q(&[1, 1]), Vec::new()))
        );
        assert_eq!(
            divrem(&q(&[1, 0, 1]), &q(&[-1, 1])),
            Some((q(&[1, 1]), q(&[2])))
        );
        assert_eq!(divrem(&q(&[1]), &[]), None);
        // gcd((x-1)(x+2), (x-1)(x+3)) = x - 1
        assert_eq!(
            gcd(
                &mul(&q(&[-1, 1]), &q(&[2, 1])),
                &mul(&q(&[-1, 1]), &q(&[3, 1]))
            ),
            q(&[-1, 1])
        );
        assert_eq!(gcd(&q(&[1, 1]), &q(&[2, 1])), q(&[1]));
    }

    #[test]
    fn square_free_decomposition() {
        // (x - 1)^2 * (x + 2)^3 * (x + 5)
        let a = q(&[-1, 1]);
        let b = q(&[2, 1]);
        let c = q(&[5, 1]);
        let p = mul(&mul(&mul(&a, &a), &mul(&mul(&b, &b), &b)), &c);
        assert_eq!(square_free(&p), vec![(c, 1), (a, 2), (b, 3)]);
        assert_eq!(square_free(&q(&[7])), []);
    }

    #[test]
    fn factorisation_reproduces_the_input() {
        let cases: [&[i64]; 9] = [
            &[-1, 0, 1],             // x^2 - 1
            &[1, 0, 1],              // x^2 + 1 (irreducible)
            &[-1, 0, 0, 0, 1],       // x^4 - 1
            &[1, 0, 0, 0, 1],        // x^4 + 1 (irreducible, splits mod every p)
            &[-2, 0, 1],             // x^2 - 2
            &[6, -5, -2, 1],         // (x-1)(x+2)(x-3)
            &[1, 1, 1, 1, 1, 1, 1],  // cyclotomic 7
            &[-1, 0, 0, 0, 0, 0, 1], // x^6 - 1
            &[4, 8, 5, 1],           // (x+1)(x+2)^2
        ];
        for coeffs in cases {
            let p = q(coeffs);
            let (unit, factors) = factor(&p);
            assert_eq!(
                product(&unit, &factors),
                p,
                "factors of {coeffs:?}: {factors:?}"
            );
            for (f, _) in &factors {
                assert!(f.last().is_some_and(Signed::is_positive));
            }
        }
    }

    #[test]
    fn known_factorisations() {
        let degrees = |coeffs: &[i64]| -> Vec<(usize, u32)> {
            factor(&q(coeffs))
                .1
                .iter()
                .map(|(f, e)| (f.len() - 1, *e))
                .collect()
        };
        assert_eq!(degrees(&[-1, 0, 1]), vec![(1, 1), (1, 1)]);
        assert_eq!(degrees(&[1, 0, 1]), vec![(2, 1)]);
        assert_eq!(degrees(&[1, 0, 0, 0, 1]), vec![(4, 1)]);
        assert_eq!(degrees(&[-1, 0, 0, 0, 1]), vec![(1, 1), (1, 1), (2, 1)]);
        assert_eq!(
            degrees(&[-1, 0, 0, 0, 0, 0, 1]),
            vec![(1, 1), (1, 1), (2, 1), (2, 1)]
        );
        assert_eq!(degrees(&[4, 8, 5, 1]), vec![(1, 1), (1, 2)]);
        // x^8 + x^4 + 1 = (x^2 + x + 1)(x^2 - x + 1)(x^4 - x^2 + 1)
        assert_eq!(
            degrees(&[1, 0, 0, 0, 1, 0, 0, 0, 1]),
            vec![(2, 1), (2, 1), (4, 1)]
        );
    }

    #[test]
    fn non_monic_and_rational_inputs() {
        // 6x^2 + 5x + 1 = (2x + 1)(3x + 1)
        let (unit, factors) = factor(&q(&[1, 5, 6]));
        assert_eq!(unit, BigRational::one());
        assert_eq!(factors, vec![(z(&[1, 2]), 1), (z(&[1, 3]), 1)]);
        // x^2/2 - 1/2 = 1/2 (x - 1)(x + 1)
        let half = BigRational::new(BigInt::from(1), BigInt::from(2));
        let p: QPoly = vec![-half.clone(), BigRational::zero(), half.clone()];
        let (unit, factors) = factor(&p);
        assert_eq!(unit, half);
        assert_eq!(product(&unit, &factors), p);
        // -x^2 + 1
        let (unit, factors) = factor(&q(&[1, 0, -1]));
        assert_eq!(unit, -BigRational::one());
        assert_eq!(product(&unit, &factors), q(&[1, 0, -1]));
    }

    #[test]
    fn products_of_random_factors_are_recovered() {
        // Multiply irreducibles together and check the factor multiset.
        let pieces: [&[i64]; 5] = [&[3, 1], &[-7, 2], &[1, 1, 1], &[2, 0, 1], &[-3, 0, 0, 1]];
        let mut seed = 12_345_u64;
        for _ in 0..20 {
            let mut p = q(&[1]);
            let mut expected: Vec<ZPoly> = Vec::new();
            for piece in pieces {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                if (seed >> 40).is_multiple_of(2) {
                    p = mul(&p, &q(piece));
                    expected.push(z(piece));
                }
            }
            let (_, factors) = factor(&p);
            let mut got: Vec<ZPoly> = factors.into_iter().map(|(f, _)| f).collect();
            got.sort();
            expected.sort();
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn roots() {
        // (2x - 1)(x + 3)(x^2 + 1)
        let p = mul(&mul(&q(&[-1, 2]), &q(&[3, 1])), &q(&[1, 0, 1]));
        let roots = rational_roots(&p);
        assert_eq!(
            roots,
            vec![
                BigRational::from_integer(BigInt::from(-3)),
                BigRational::new(BigInt::from(1), BigInt::from(2))
            ]
        );
        assert!(roots.iter().all(|r| eval(&p, r).is_zero()));
    }
}
