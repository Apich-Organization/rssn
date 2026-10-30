//! Hermite reduction in `Q(x)(y)`, `y^2 = D`, `D` square-free.
//!
//! The integrand is written `(A + B y)/C`.  For an irreducible square-free
//! factor `w` of `C` with multiplicity `k >= 2` and `gcd(w, D) = 1` the
//! derivative
//!
//! ```text
//! d/dx [ (U + V y) / w^(k-1) ]
//!     = [ (U + V y)' w - (k-1) (U + V y) w' ] / w^k
//! ```
//!
//! has the same principal part at `w` as `(A + B y)/w^k` exactly when
//!
//! ```text
//! A / (C/w^k) + (k-1) w' U = 0   (mod w)
//! B / (C/w^k) + (k-1) w' V = 0   (mod w)
//! ```
//!
//! which determines `U` and `V` modulo `w` uniquely over `Q`, because `w`
//! is square-free and `k-1 != 0`.  Subtracting the derivative lowers the
//! multiplicity of `w` by one, so the recursion terminates.  The condition
//! `gcd(w, D) = 1` keeps `y' = D'/(2D) y` regular modulo `w`; a repeated
//! factor that divides `D` is a ramified place and is left for the caller's
//! boundary case.

use num_rational::BigRational;

use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

use super::field;
use super::field::Element;
use super::field::Frac;

/// The most Hermite steps taken before giving up on a pathologically large
/// input.
const MAX_STEPS: usize = 64;
/// Degree beyond which the reduction gives up rather than grow without
/// bound.
const MAX_DEGREE: usize = 400;

/// One subtracted exact derivative `(U + V y) / w^power`.
#[derive(Clone, Debug)]
pub(super) struct Piece {
    /// Coefficient of `1`.
    pub u: QPoly,
    /// Coefficient of `y`.
    pub v: QPoly,
    /// The irreducible factor.
    pub w: QPoly,
    /// The exponent of `w` in the denominator (`k - 1`).
    pub power: u32,
}

/// The outcome of the reduction.
#[derive(Clone, Debug)]
pub(super) struct Reduced {
    /// Exact derivatives subtracted from the integrand.
    pub pieces: Vec<Piece>,
    /// Numerator of the remaining integrand.
    pub a: QPoly,
    /// Coefficient of `y` of the remaining integrand.
    pub b: QPoly,
    /// Its denominator.
    pub c: QPoly,
    /// Number of reduction steps taken.
    pub steps: usize,
}

/// Converts the factor returned by `univariate::factor` to a monic
/// polynomial over `Q`.
fn monic_factor(f: &[num_bigint::BigInt]) -> QPoly {
    let q: QPoly = f
        .iter()
        .map(|c| BigRational::from_integer(c.clone()))
        .collect::<QPoly>();
    field::monic(&q)
}

/// The derivative of a rational function, exactly.
fn frac_derivative(f: &Frac) -> Option<Frac> {
    let numerator = field::sub(
        &field::mul(&field::derivative(&f.num), &f.den),
        &field::mul(&f.num, &field::derivative(&f.den)),
    );
    let denominator = field::mul(&f.den, &f.den);
    Frac::new(&numerator, &denominator)
}

/// The derivative of `a + b y` in `Q(x)(y)`, exactly.
fn element_derivative(
    element: &Element,
    d: &[BigRational],
) -> Option<Element> {
    let a = frac_derivative(&element.a)?;
    let bprime = frac_derivative(&element.b)?;
    let two_d = field::scale(d, &BigRational::from_integer(num_bigint::BigInt::from(2)));
    let dp = Frac::polynomial(&field::derivative(d));
    let tail = element.b.mul(&dp)?.div(&Frac::polynomial(&two_d))?;
    Some(Element {
        a,
        b: bprime.add(&tail)?,
    })
}

/// Hermite reduction of `(a + b y)/c` with respect to `d`.
pub(super) fn reduce(
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> Reduced {
    let Some(mut element) = Element::quotient(a, b, c) else {
        return Reduced {
            pieces: Vec::new(),
            a: a.to_vec(),
            b: b.to_vec(),
            c: c.to_vec(),
            steps: 0,
        };
    };
    let mut pieces = Vec::new();
    let mut steps = 0_usize;
    while steps < MAX_STEPS {
        let Some((cur_a, cur_b, cur_c)) = element.to_quotient() else {
            break;
        };
        if field::degree(&cur_c).is_none_or(|k| k > MAX_DEGREE) {
            break;
        }
        let (_, factors) = univariate::factor(&cur_c);
        let mut chosen: Option<(QPoly, u32)> = None;
        for (f, e) in &factors {
            if *e < 2 {
                continue;
            }
            let w = monic_factor(f);
            let shared = field::gcd(&w, d);
            if field::degree(&shared).is_none_or(|k| k == 0) {
                chosen = Some((w, *e));
                break;
            }
        }
        let Some((w, k)) = chosen else {
            break;
        };
        let wk = field::pow(&w, k);
        let Some(c1) = field::div_exact(&cur_c, &wk) else {
            break;
        };
        let wp = field::derivative(&w);
        let s = field::scale(
            &wp,
            &BigRational::from_integer(num_bigint::BigInt::from(-i64::from(k) + 1)),
        );
        let (Some(inv_s), Some(inv_c1)) = (field::inverse_mod(&s, &w), field::inverse_mod(&c1, &w))
        else {
            break;
        };
        let Some(t) = field::rem(&field::mul(&inv_s, &inv_c1), &w) else {
            break;
        };
        let (Some(u), Some(v)) = (
            field::rem(&field::mul(&cur_a, &t), &w),
            field::rem(&field::mul(&cur_b, &t), &w),
        ) else {
            break;
        };
        // Derivative of (u + v y)/w^(k-1), over w^k.
        let up = field::derivative(&u);
        let vp = field::derivative(&v);
        let k_minus_1 = BigRational::from_integer(num_bigint::BigInt::from(i64::from(k) - 1));
        let p = field::sub(
            &field::mul(&up, &w),
            &field::scale(&field::mul(&u, &wp), &k_minus_1),
        );
        let q_poly = field::sub(
            &field::mul(&vp, &w),
            &field::scale(&field::mul(&v, &wp), &k_minus_1),
        );
        let two_d = field::scale(d, &BigRational::from_integer(num_bigint::BigInt::from(2)));
        let q_num = field::add(
            &field::mul(&q_poly, &two_d),
            &field::mul(&field::mul(&v, &w), &field::derivative(d)),
        );
        let Some(deriv_a) = Frac::new(&p, &wk) else {
            break;
        };
        let den_b = field::mul(&two_d, &wk);
        let Some(deriv_b) = Frac::new(&q_num, &den_b) else {
            break;
        };
        // Exact check that the piece really differentiates to the quantity
        // subtracted from the integrand.
        let denominator = field::pow(&w, k.saturating_sub(1));
        let (Some(piece_a), Some(piece_b)) =
            (Frac::new(&u, &denominator), Frac::new(&v, &denominator))
        else {
            break;
        };
        let piece = Element {
            a: piece_a,
            b: piece_b,
        };
        let Some(exact) = element_derivative(&piece, d) else {
            break;
        };
        if exact.a != deriv_a || exact.b != deriv_b {
            break;
        }
        let (Some(next_a), Some(next_b)) = (element.a.sub(&deriv_a), element.b.sub(&deriv_b))
        else {
            break;
        };
        let next = Element { a: next_a, b: next_b };
        // Confirm the multiplicity really dropped before accepting the step:
        // the new denominator must no longer be divisible by `w^k`.
        let Some((_, _, next_c)) = next.to_quotient() else {
            break;
        };
        if field::div_exact(&next_c, &wk).is_some() {
            break;
        }
        element = next;
        pieces.push(Piece {
            u,
            v,
            w,
            power: k.saturating_sub(1),
        });
        steps = steps.saturating_add(1);
    }
    match element.to_quotient() {
        | Some((a, b, c)) => Reduced {
            pieces,
            a,
            b,
            c,
            steps,
        },
        | None => Reduced {
            pieces: Vec::new(),
            a: a.to_vec(),
            b: b.to_vec(),
            c: c.to_vec(),
            steps: 0,
        },
    }
}
