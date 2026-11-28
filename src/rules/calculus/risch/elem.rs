//! Elements of a tower of monomial extensions over `Q`.
//!
//! Depth 0 is `Q`. Depth `d >= 1` is `K_{d-1}(theta_d)` where `theta_d` is a
//! monomial with `theta_d' = a*theta_d + b` for `a`, `b` in `K_{d-1}` (so a
//! logarithm has `a = 0`, an exponential has `b = 0`, and the variable `x`
//! itself is the primitive monomial `a = 0`, `b = 1`).
//!
//! An element is kept as a reduced fraction of polynomials in `theta_d` with
//! coefficients of depth `d - 1`.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use num_bigint::BigInt;

use super::Q;
use super::flat::Flat;
use super::poly;

/// A coefficient field the polynomial algorithms work over.
pub trait Coeff: Clone + PartialEq + std::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn add(
        &self,
        other: &Self,
    ) -> Self;
    fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }
    fn neg(&self) -> Self;
    fn mul(
        &self,
        other: &Self,
    ) -> Self;
    fn inv(&self) -> Option<Self>;
    fn pow_i(
        &self,
        e: i64,
    ) -> Option<Self>;
    fn from_q(value: &Q) -> Self;
    /// The value as a rational number, when the element is a constant.
    fn as_q(&self) -> Option<Q>;
    /// `d/dx`.
    fn deriv(
        &self,
        tower: &super::tower::Tower,
    ) -> Self;
}

impl Coeff for Q {
    fn zero() -> Self {
        Q::from_integer(BigInt::from(0))
    }

    fn one() -> Self {
        Q::from_integer(BigInt::from(1))
    }

    fn is_zero(&self) -> bool {
        num_traits::Zero::is_zero(self)
    }

    fn is_one(&self) -> bool {
        num_traits::One::is_one(self)
    }

    fn add(
        &self,
        other: &Self,
    ) -> Self {
        self + other
    }

    fn neg(&self) -> Self {
        -self
    }

    fn mul(
        &self,
        other: &Self,
    ) -> Self {
        self * other
    }

    fn inv(&self) -> Option<Self> {
        if num_traits::Zero::is_zero(self) {
            None
        } else {
            Some(self.recip())
        }
    }

    fn pow_i(
        &self,
        e: i64,
    ) -> Option<Self> {
        if e >= 0 {
            Some(num_traits::pow::pow(self.clone(), usize::try_from(e).ok()?))
        } else {
            self.inv()
                .map(|v| num_traits::pow::pow(v, usize::try_from(-e).unwrap_or(0)))
        }
    }

    fn from_q(value: &Q) -> Self {
        value.clone()
    }

    fn as_q(&self) -> Option<Q> {
        Some(self.clone())
    }

    fn deriv(
        &self,
        _tower: &super::tower::Tower,
    ) -> Self {
        Self::zero()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Elem {
    /// A rational number.
    Const(Q),
    /// A reduced fraction of polynomials in `theta_depth`.
    Ext(Ext),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ext {
    pub depth: usize,
    pub num: Vec<Elem>,
    pub den: Vec<Elem>,
}

impl Elem {
    pub fn depth(&self) -> usize {
        match self {
            | Self::Const(_) => 0,
            | Self::Ext(e) => e.depth,
        }
    }

    pub fn zero_at(depth: usize) -> Self {
        if depth == 0 {
            Self::Const(Q::from_integer(BigInt::from(0)))
        } else {
            Self::Ext(Ext {
                depth,
                num: Vec::new(),
                den: vec![Self::one_at(depth - 1)],
            })
        }
    }

    pub fn one_at(depth: usize) -> Self {
        if depth == 0 {
            Self::Const(Q::from_integer(BigInt::from(1)))
        } else {
            let inner = Self::one_at(depth - 1);
            Self::Ext(Ext {
                depth,
                num: vec![inner.clone()],
                den: vec![inner],
            })
        }
    }

    /// The generator `theta_depth`.
    pub fn gen_at(depth: usize) -> Self {
        if depth == 0 {
            return Self::zero_at(0);
        }
        let zero = Self::zero_at(depth - 1);
        let one = Self::one_at(depth - 1);
        Self::Ext(Ext {
            depth,
            num: vec![zero, one.clone()],
            den: vec![one],
        })
    }

    /// Raises an element of a lower level into `depth`.
    pub fn embed(
        &self,
        depth: usize,
    ) -> Self {
        if self.depth() >= depth {
            return self.clone();
        }
        if self.is_zero() {
            return Self::zero_at(depth);
        }
        let mut out = self.clone();
        for level in out.depth() + 1..=depth {
            let one = Self::one_at(level - 1);
            out = Self::Ext(Ext {
                depth: level,
                num: vec![out],
                den: vec![one],
            });
        }
        out
    }

    /// Numerator and denominator as polynomials in the top generator.
    pub fn parts(&self) -> Option<(&[Self], &[Self])> {
        match self {
            | Self::Ext(e) => Some((&e.num, &e.den)),
            | Self::Const(_) => None,
        }
    }

    pub fn is_polynomial(&self) -> bool {
        match self {
            | Self::Ext(e) => e.den.len() == 1 && e.den[0].is_one(),
            | Self::Const(_) => true,
        }
    }

    /// The numerator of the top-level fraction when the denominator is one.
    pub fn polynomial_parts(&self) -> Option<&[Self]> {
        match self {
            | Self::Ext(e) if e.den.len() == 1 && e.den[0].is_one() => Some(&e.num),
            | _ => None,
        }
    }

    /// `(numerator, denominator)` as sparse polynomials in `x, theta_1, ...`.
    pub fn flatten(&self) -> (Flat, Flat) {
        match self {
            | Self::Const(q) => (Flat::constant(q.clone()), Flat::one()),
            | Self::Ext(e) => {
                let var = u16::try_from(e.depth - 1).unwrap_or(0);
                let mut num = (Flat::zero(), Flat::one());
                for (i, c) in e.num.iter().enumerate() {
                    let (cn, cd) = c.flatten();
                    let power = Flat::monomial(var, u32::try_from(i).unwrap_or(0));
                    // num += (cn/cd) * theta^i
                    num.0 = Flat::add(
                        &Flat::mul(&num.0, &cd),
                        &Flat::mul(&Flat::mul(&cn, &power), &num.1),
                    );
                    num.1 = Flat::mul(&num.1, &cd);
                }
                let mut den = (Flat::zero(), Flat::one());
                for (i, c) in e.den.iter().enumerate() {
                    let (cn, cd) = c.flatten();
                    let power = Flat::monomial(var, u32::try_from(i).unwrap_or(0));
                    den.0 = Flat::add(
                        &Flat::mul(&den.0, &cd),
                        &Flat::mul(&Flat::mul(&cn, &power), &den.1),
                    );
                    den.1 = Flat::mul(&den.1, &cd);
                }
                // (num.0/num.1) / (den.0/den.1)
                let n = Flat::mul(&num.0, &den.1);
                let d = Flat::mul(&num.1, &den.0);
                super::flat::reduce_pair(&n, &d)
            },
        }
    }

    /// Reduces `num/den` into canonical form.
    pub fn ratio(
        depth: usize,
        num: Vec<Self>,
        den: Vec<Self>,
    ) -> Option<Self> {
        if depth == 0 {
            let a = num.first()?.as_q()?;
            let b = den.first()?.as_q()?;
            if b.is_zero() {
                return None;
            }
            return Some(Self::Const(a / b));
        }
        let num = poly::trim(num);
        let den = poly::trim(den);
        if den.is_empty() {
            return None;
        }
        if den.len() == 1 {
            if den[0].is_one() {
                return Some(Self::Ext(Ext {
                    depth,
                    num,
                    den: vec![Self::one_at(depth - 1)],
                }));
            }
            let inv = den[0].inv()?;
            return Some(Self::Ext(Ext {
                depth,
                num: poly::scale(&num, &inv),
                den: vec![Self::one_at(depth - 1)],
            }));
        }
        let g = poly::gcd(&num, &den);
        let (num, den) = if poly::degree(&g).is_some_and(|d| d > 0) {
            (
                poly::exact_div(&num, &g).unwrap_or(num),
                poly::exact_div(&den, &g).unwrap_or(den),
            )
        } else {
            (num, den)
        };
        if den.len() == 1 {
            let inv = den[0].inv()?;
            return Some(Self::Ext(Ext {
                depth,
                num: poly::scale(&num, &inv),
                den: vec![Self::one_at(depth - 1)],
            }));
        }
        match den.last().and_then(Coeff::inv) {
            | Some(inv) => {
                Some(Self::Ext(Ext {
                    depth,
                    num: poly::scale(&num, &inv),
                    den: poly::scale(&den, &inv),
                }))
            },
            | None => None,
        }
    }
}

impl Coeff for Elem {
    fn zero() -> Self {
        Self::Const(Q::from_integer(BigInt::from(0)))
    }

    fn one() -> Self {
        Self::Const(Q::from_integer(BigInt::from(1)))
    }

    fn is_zero(&self) -> bool {
        match self {
            | Self::Const(q) => num_traits::Zero::is_zero(q),
            | Self::Ext(e) => e.num.iter().all(Coeff::is_zero),
        }
    }

    fn is_one(&self) -> bool {
        self.as_q().is_some_and(|q| num_traits::One::is_one(&q))
    }

    fn add(
        &self,
        other: &Self,
    ) -> Self {
        if other.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return other.clone();
        }
        let depth = self.depth().max(other.depth());
        let (a, b) = (self.embed(depth), other.embed(depth));
        if depth == 0 {
            let (Self::Const(x), Self::Const(y)) = (&a, &b) else {
                unreachable!("depth zero elements are constants")
            };
            return Self::Const(x + y);
        }
        let (Some((an, ad)), Some((bn, bd))) = (a.parts(), b.parts()) else {
            unreachable!("depth is positive")
        };
        if an.is_empty() {
            return b;
        }
        if bn.is_empty() {
            return a;
        }
        let num = poly::add(&poly::mul(an, bd), &poly::mul(bn, ad));
        let den = poly::mul(ad, bd);
        Self::ratio(depth, num, den).unwrap_or_else(|| a.clone())
    }

    fn neg(&self) -> Self {
        match self {
            | Self::Const(q) => Self::Const(-q),
            | Self::Ext(e) => {
                Self::Ext(Ext {
                    depth: e.depth,
                    num: poly::neg(&e.num),
                    den: e.den.clone(),
                })
            },
        }
    }

    fn mul(
        &self,
        other: &Self,
    ) -> Self {
        if self.is_one() {
            return other.clone();
        }
        if other.is_one() {
            return self.clone();
        }
        if self.is_zero() || other.is_zero() {
            return Self::zero_at(self.depth().max(other.depth()));
        }
        let depth = self.depth().max(other.depth());
        let (a, b) = (self.embed(depth), other.embed(depth));
        if depth == 0 {
            let (Self::Const(x), Self::Const(y)) = (&a, &b) else {
                unreachable!("depth zero elements are constants")
            };
            return Self::Const(x * y);
        }
        let (Some((an, ad)), Some((bn, bd))) = (a.parts(), b.parts()) else {
            unreachable!("depth is positive")
        };
        if an.is_empty() || bn.is_empty() {
            return Self::zero_at(depth);
        }
        let num = poly::mul(an, bn);
        let den = poly::mul(ad, bd);
        Self::ratio(depth, num, den).unwrap_or_else(|| Self::zero_at(depth))
    }

    fn inv(&self) -> Option<Self> {
        match self {
            | Self::Const(q) => {
                if num_traits::Zero::is_zero(q) {
                    None
                } else {
                    Some(Self::Const(q.recip()))
                }
            },
            | Self::Ext(e) => {
                if self.is_zero() {
                    return None;
                }
                Self::ratio(e.depth, e.den.clone(), e.num.clone())
            },
        }
    }

    fn pow_i(
        &self,
        e: i64,
    ) -> Option<Self> {
        let base = if e < 0 {
            self.inv()?
        } else {
            self.clone()
        };
        let mut acc = Self::one_at(base.depth());
        for _ in 0..e.unsigned_abs() {
            acc = acc.mul(&base);
        }
        Some(acc)
    }

    fn from_q(value: &Q) -> Self {
        Self::Const(value.clone())
    }

    fn as_q(&self) -> Option<Q> {
        match self {
            | Self::Const(q) => Some(q.clone()),
            | Self::Ext(e) => {
                let num = poly::trim_ref(&e.num);
                let den = poly::trim_ref(&e.den);
                if num.is_empty() {
                    return Some(Q::from_integer(BigInt::from(0)));
                }
                if num.len() == 1 && den.len() == 1 {
                    let a = num[0].as_q()?;
                    let b = den[0].as_q()?;
                    if b.is_zero() {
                        None
                    } else {
                        Some(a / b)
                    }
                } else {
                    None
                }
            },
        }
    }

    fn deriv(
        &self,
        tower: &super::tower::Tower,
    ) -> Self {
        match self {
            | Self::Const(_) => Self::zero(),
            | Self::Ext(e) => {
                let d = e.depth;
                let (a, b) = tower.monomial(d);
                let num_d = poly_derivative(&e.num, &a, &b, tower);
                let den_d = poly_derivative(&e.den, &a, &b, tower);
                let num = poly::sub(&poly::mul(&num_d, &e.den), &poly::mul(&e.num, &den_d));
                let den = poly::mul(&e.den, &e.den);
                Self::ratio(d, num, den).unwrap_or_else(Self::zero)
            },
        }
    }
}

/// `d/dx` of a polynomial in `theta` with `theta' = a*theta + b`.
fn poly_derivative(
    p: &[Elem],
    a: &Elem,
    b: &Elem,
    tower: &super::tower::Tower,
) -> Vec<Elem> {
    let mut out = Vec::with_capacity(p.len());
    for (i, c) in p.iter().enumerate() {
        let mut term = c.deriv(tower);
        if i > 0 {
            let k = Elem::Const(Q::from_integer(BigInt::from(i)));
            term = term.add(&k.mul(&a.mul(c)));
        }
        if let Some(next) = p.get(i + 1) {
            let k = Elem::Const(Q::from_integer(BigInt::from(i + 1)));
            term = term.add(&k.mul(&b.mul(next)));
        }
        out.push(term);
    }
    poly::trim(out)
}
