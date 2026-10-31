//! Exact arithmetic in `Q(x)(y)`, `y^2 = D`: dense polynomials over `Q`,
//! rational functions and the two-dimensional algebra over them.
//!
//! Everything here is exact (`BigRational`); no floating point is used
//! until a candidate antiderivative is checked numerically.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::rules::poly::univariate as uni;
use crate::rules::poly::univariate::QPoly;

/// `1`.
pub(super) fn one() -> BigRational {
    BigRational::from_integer(BigInt::from(1))
}

/// Degree, or `None` for the zero polynomial.
pub(super) const fn degree(p: &[BigRational]) -> Option<usize> {
    p.len().checked_sub(1)
}

/// Sum.
pub(super) fn add(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    uni::add(a, b)
}

/// Difference.
pub(super) fn sub(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    uni::sub(a, b)
}

/// Product.
pub(super) fn mul(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    uni::mul(a, b)
}

/// Additive inverse.
pub(super) fn neg(a: &[BigRational]) -> QPoly {
    a.iter().map(|c| -c).collect()
}

/// Product with a number.
pub(super) fn scale(
    a: &[BigRational],
    c: &BigRational,
) -> QPoly {
    if c.is_zero() {
        return Vec::new();
    }
    let mut out: QPoly = a.iter().map(|x| x * c).collect();
    while out.last().is_some_and(Zero::is_zero) {
        out.pop();
    }
    out
}

/// Derivative.
pub(super) fn derivative(p: &[BigRational]) -> QPoly {
    uni::derivative(p)
}

/// Monic form.
pub(super) fn monic(p: &[BigRational]) -> QPoly {
    uni::monic(p)
}

/// The leading coefficient.
pub(super) const fn leading(p: &[BigRational]) -> Option<&BigRational> {
    p.last()
}

/// `a mod m`.
pub(super) fn rem(
    a: &[BigRational],
    m: &[BigRational],
) -> Option<QPoly> {
    let (_, r) = uni::divrem(a, m)?;
    Some(r)
}

/// The exact quotient `a / b`, or `None` when `b` does not divide `a`.
pub(super) fn div_exact(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<QPoly> {
    if b.is_empty() {
        return None;
    }
    let (q, r) = uni::divrem(a, b)?;
    r.is_empty().then_some(q)
}

/// Quotient and remainder.
pub(super) fn divrem(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    uni::divrem(a, b)
}

/// Monic greatest common divisor.
pub(super) fn gcd(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    uni::gcd(a, b)
}

/// `a^e` for a small exponent.
pub(super) fn pow(
    a: &[BigRational],
    mut e: u32,
) -> QPoly {
    let mut out: QPoly = vec![one()];
    let mut base: QPoly = a.to_vec();
    while e > 0 {
        if e & 1 == 1 {
            out = mul(&out, &base);
        }
        e >>= 1;
        if e > 0 {
            base = mul(&base, &base);
        }
    }
    out
}

/// The inverse of `a` modulo the polynomial `m`, or `None` when they are
/// not coprime.
pub(super) fn inverse_mod(
    a: &[BigRational],
    m: &[BigRational],
) -> Option<QPoly> {
    let a = rem(a, m)?;
    if a.is_empty() {
        return None;
    }
    let (mut r0, mut r1) = (a, m.to_vec());
    let (mut s0, mut s1): (QPoly, QPoly) = (vec![one()], Vec::new());
    let mut guard = 0_u32;
    while !r1.is_empty() {
        guard = guard.checked_add(1)?;
        if guard > 512 {
            return None;
        }
        let (q, r) = divrem(&r0, &r1)?;
        let s2 = sub(&s0, &mul(&q, &s1));
        (r0, r1) = (r1, r);
        (s0, s1) = (s1, s2);
    }
    match r0.as_slice() {
        | [c] if !c.is_zero() => {
            let inv = s0.iter().map(|x| x / c).collect::<QPoly>();
            rem(&inv, m)
        },
        | _ => None,
    }
}

/// Whether the rational number is a square of a rational number.
pub(super) fn rational_sqrt(q: &BigRational) -> Option<BigRational> {
    if q.is_negative() {
        return None;
    }
    let (n, d) = (q.numer().abs(), q.denom().abs());
    let (rn, rd) = (n.sqrt(), d.sqrt());
    if &rn * &rn == n && &rd * &rd == d {
        Some(BigRational::new(rn, rd))
    } else {
        None
    }
}

/// A rational function `num / den` in lowest terms, with `den` monic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Frac {
    /// Numerator.
    pub num: QPoly,
    /// Denominator, monic and coprime to the numerator.
    pub den: QPoly,
}

impl Frac {
    /// The zero rational function.
    pub(super) fn zero() -> Self {
        Self {
            num: Vec::new(),
            den: vec![one()],
        }
    }

    /// The constant `c`.
    pub(super) fn constant(c: BigRational) -> Self {
        Self {
            num: if c.is_zero() {
                Vec::new()
            } else {
                vec![c]
            },
            den: vec![one()],
        }
    }

    /// The polynomial `p`.
    pub(super) fn polynomial(p: &[BigRational]) -> Self {
        Self {
            num: p.to_vec(),
            den: vec![one()],
        }
    }

    /// `num / den` in lowest terms, or `None` for a zero denominator.
    pub(super) fn new(
        num: &[BigRational],
        den: &[BigRational],
    ) -> Option<Self> {
        let lead = leading(den)?;
        if lead.is_zero() {
            return None;
        }
        let inverse = lead.recip();
        let mut n = scale(num, &inverse);
        let mut d = monic(den);
        let g = gcd(&n, &d);
        if degree(&g).is_some_and(|k| k > 0) {
            n = div_exact(&n, &g)?;
            d = div_exact(&d, &g)?;
        }
        Some(Self { num: n, den: d })
    }

    /// Whether this is the zero function.
    pub(super) const fn is_zero(&self) -> bool {
        self.num.is_empty()
    }

    /// Sum.
    pub(super) fn add(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let n = add(&mul(&self.num, &other.den), &mul(&other.num, &self.den));
        let d = mul(&self.den, &other.den);
        Self::new(&n, &d)
    }

    /// Product.
    pub(super) fn mul(
        &self,
        other: &Self,
    ) -> Option<Self> {
        Self::new(&mul(&self.num, &other.num), &mul(&self.den, &other.den))
    }

    /// Difference.
    pub(super) fn sub(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let n = sub(&mul(&self.num, &other.den), &mul(&other.num, &self.den));
        let d = mul(&self.den, &other.den);
        Self::new(&n, &d)
    }

    /// Quotient.
    pub(super) fn div(
        &self,
        other: &Self,
    ) -> Option<Self> {
        if other.is_zero() {
            return None;
        }
        Self::new(&mul(&self.num, &other.den), &mul(&self.den, &other.num))
    }

    /// Additive inverse.
    pub(super) fn neg(&self) -> Self {
        Self {
            num: neg(&self.num),
            den: self.den.clone(),
        }
    }

    /// Product with a polynomial.
    pub(super) fn scale_poly(
        &self,
        p: &[BigRational],
    ) -> Option<Self> {
        Self::new(&mul(&self.num, p), &self.den)
    }
}

/// An element `a + b*y` of `Q(x)(y)`, `y^2 = D`.
#[derive(Clone, Debug)]
pub(super) struct Element {
    /// The coefficient of `1`.
    pub a: Frac,
    /// The coefficient of `y`.
    pub b: Frac,
}

impl Element {
    /// Zero.
    pub(super) fn zero() -> Self {
        Self {
            a: Frac::zero(),
            b: Frac::zero(),
        }
    }

    /// The constant `c`.
    pub(super) fn constant(c: BigRational) -> Self {
        Self {
            a: Frac::constant(c),
            b: Frac::zero(),
        }
    }

    /// `num / den` with both parts polynomials, over the common denominator
    /// `den`.
    pub(super) fn quotient(
        a: &[BigRational],
        b: &[BigRational],
        den: &[BigRational],
    ) -> Option<Self> {
        Some(Self {
            a: Frac::new(a, den)?,
            b: Frac::new(b, den)?,
        })
    }

    /// Sum.
    pub(super) fn add(
        &self,
        other: &Self,
    ) -> Option<Self> {
        Some(Self {
            a: self.a.add(&other.a)?,
            b: self.b.add(&other.b)?,
        })
    }

    /// Product, using `y^2 = D`.
    pub(super) fn mul(
        &self,
        other: &Self,
        d: &[BigRational],
    ) -> Option<Self> {
        let aa = self.a.mul(&other.a)?;
        let bb = self.b.mul(&other.b)?;
        let cross = self.a.mul(&other.b)?.add(&self.b.mul(&other.a)?)?;
        let d_frac = Frac::polynomial(d);
        Some(Self {
            a: aa.add(&bb.mul(&d_frac)?)?,
            b: cross,
        })
    }

    /// Difference (used by the exactness checks in the tests).
    #[cfg(test)]
    pub(super) fn sub(
        &self,
        other: &Self,
    ) -> Option<Self> {
        Some(Self {
            a: self.a.sub(&other.a)?,
            b: self.b.sub(&other.b)?,
        })
    }

    /// The inverse `(a - b y) / (a^2 - b^2 D)`, or `None` at zero.
    pub(super) fn inverse(
        &self,
        d: &[BigRational],
    ) -> Option<Self> {
        let aa = self.a.mul(&self.a)?;
        let bb = self.b.mul(&self.b)?;
        let norm = aa.sub(&bb.mul(&Frac::polynomial(d))?)?;
        if norm.is_zero() {
            return None;
        }
        Some(Self {
            a: self.a.div(&norm)?,
            b: self.b.neg().div(&norm)?,
        })
    }

    /// `self^n` for an integer `n`.
    pub(super) fn pow(
        &self,
        n: i64,
        d: &[BigRational],
    ) -> Option<Self> {
        let mut e = n.unsigned_abs();
        let mut base = if n < 0 {
            self.inverse(d)?
        } else {
            self.clone()
        };
        let mut out = Self::constant(one());
        while e > 0 {
            if e & 1 == 1 {
                out = out.mul(&base, d)?;
            }
            e >>= 1;
            if e > 0 {
                base = base.mul(&base, d)?;
            }
        }
        Some(out)
    }

    /// The element in the form `(A + B y)/C` with `C` monic and
    /// `gcd(gcd(A, B), C) = 1`.
    pub(super) fn to_quotient(&self) -> Option<(QPoly, QPoly, QPoly)> {
        let g = gcd(&self.a.den, &self.b.den);
        let l = div_exact(&mul(&self.a.den, &self.b.den), &g)?;
        let ca = div_exact(&l, &self.a.den)?;
        let cb = div_exact(&l, &self.b.den)?;
        let mut a = mul(&self.a.num, &ca);
        let mut b = mul(&self.b.num, &cb);
        let mut c = l;
        let g1 = gcd(&a, &b);
        let g2 = gcd(&g1, &c);
        if degree(&g2).is_some_and(|k| k > 0) {
            a = div_exact(&a, &g2)?;
            b = div_exact(&b, &g2)?;
            c = div_exact(&c, &g2)?;
        }
        let lead = leading(&c)?.clone();
        a = scale(&a, &lead.recip());
        b = scale(&b, &lead.recip());
        c = monic(&c);
        Some((a, b, c))
    }
}
