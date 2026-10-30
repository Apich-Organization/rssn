//! The numeric tower carried by literal leaves.
//!
//! Exact values (`Int`, `Rat`) stay exact under every operation that can be
//! performed exactly; anything else degrades to `Float`. This is the single
//! place where "exact phase" and "numeric phase" arithmetic meet.

use std::cmp::Ordering;
use std::fmt;
use std::hash::Hash;
use std::hash::Hasher;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

/// A literal number: exact integer, exact rational, or binary float.
#[derive(Clone, Debug)]
pub enum Number {
    /// Arbitrary precision integer.
    Int(BigInt),
    /// Arbitrary precision rational, always with denominator `> 1`.
    Rat(BigRational),
    /// IEEE-754 double. Never produced from exact inputs by exact operations.
    Float(f64),
}

impl PartialEq for Number {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        match (self, other) {
            | (Self::Int(a), Self::Int(b)) => a == b,
            | (Self::Rat(a), Self::Rat(b)) => a == b,
            | (Self::Float(a), Self::Float(b)) => a.to_bits() == b.to_bits(),
            | _ => false,
        }
    }
}

impl Eq for Number {}

impl Hash for Number {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {
        match self {
            | Self::Int(a) => {
                0_u8.hash(state);
                a.hash(state);
            },
            | Self::Rat(a) => {
                1_u8.hash(state);
                a.hash(state);
            },
            | Self::Float(a) => {
                2_u8.hash(state);
                a.to_bits().hash(state);
            },
        }
    }
}

impl fmt::Display for Number {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Int(a) => write!(f, "{a}"),
            | Self::Rat(a) => write!(f, "{}/{}", a.numer(), a.denom()),
            | Self::Float(a) => write!(f, "{a}"),
        }
    }
}

impl From<i64> for Number {
    fn from(v: i64) -> Self {
        Self::Int(BigInt::from(v))
    }
}

impl From<f64> for Number {
    fn from(v: f64) -> Self {
        Self::Float(v)
    }
}

impl From<BigInt> for Number {
    fn from(v: BigInt) -> Self {
        Self::Int(v)
    }
}

impl From<BigRational> for Number {
    fn from(v: BigRational) -> Self {
        Self::rat(v)
    }
}

impl Number {
    /// Builds a number from a rational, collapsing integral values to `Int`.
    #[must_use]
    pub fn rat(r: BigRational) -> Self {
        if r.is_integer() {
            Self::Int(r.to_integer())
        } else {
            Self::Rat(r)
        }
    }

    /// Builds the exact fraction `n / d`. Returns `None` when `d == 0`.
    #[must_use]
    pub fn fraction(
        n: i64,
        d: i64,
    ) -> Option<Self> {
        if d == 0 {
            return None;
        }
        Some(Self::rat(BigRational::new(
            BigInt::from(n),
            BigInt::from(d),
        )))
    }

    /// Returns `true` for `Int` and `Rat`.
    #[must_use]
    pub const fn is_exact(&self) -> bool {
        !matches!(self, Self::Float(_))
    }

    /// Returns `true` for integers.
    #[must_use]
    pub const fn is_integer(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    /// Returns `true` when the value equals zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match self {
            | Self::Int(a) => a.is_zero(),
            | Self::Rat(a) => a.is_zero(),
            | Self::Float(a) => *a == 0.0,
        }
    }

    /// Returns `true` when the value equals one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        match self {
            | Self::Int(a) => a.is_one(),
            | Self::Rat(_) => false,
            | Self::Float(a) => (*a - 1.0).abs() == 0.0,
        }
    }

    /// Returns `true` when the value is strictly negative.
    #[must_use]
    pub fn is_negative(&self) -> bool {
        match self {
            | Self::Int(a) => a.is_negative(),
            | Self::Rat(a) => a.is_negative(),
            | Self::Float(a) => *a < 0.0,
        }
    }

    /// Nearest `f64`. Huge exact values saturate to infinity.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        match self {
            | Self::Int(a) => a.to_f64().unwrap_or(f64::NAN),
            | Self::Rat(a) => a.to_f64().unwrap_or(f64::NAN),
            | Self::Float(a) => *a,
        }
    }

    /// The value as an `i64` if it is an integer in range.
    #[must_use]
    pub fn to_i64(&self) -> Option<i64> {
        match self {
            | Self::Int(a) => a.to_i64(),
            | _ => None,
        }
    }

    /// The exact value as a rational; `None` for floats.
    #[must_use]
    pub fn to_rational(&self) -> Option<BigRational> {
        match self {
            | Self::Int(a) => Some(BigRational::from_integer(a.clone())),
            | Self::Rat(a) => Some(a.clone()),
            | Self::Float(_) => None,
        }
    }

    fn binary(
        &self,
        other: &Self,
        int: impl FnOnce(&BigInt, &BigInt) -> BigInt,
        rat: impl FnOnce(BigRational, BigRational) -> BigRational,
        float: impl FnOnce(f64, f64) -> f64,
    ) -> Self {
        match (self, other) {
            | (Self::Int(a), Self::Int(b)) => Self::Int(int(a, b)),
            | _ => match (self.to_rational(), other.to_rational()) {
                | (Some(a), Some(b)) => Self::rat(rat(a, b)),
                | _ => Self::Float(float(self.to_f64(), other.to_f64())),
            },
        }
    }

    /// Sum.
    #[must_use]
    pub fn add(
        &self,
        other: &Self,
    ) -> Self {
        self.binary(other, |a, b| a + b, |a, b| a + b, |a, b| a + b)
    }

    /// Product.
    #[must_use]
    pub fn mul(
        &self,
        other: &Self,
    ) -> Self {
        self.binary(other, |a, b| a * b, |a, b| a * b, |a, b| a * b)
    }

    /// Additive inverse.
    #[must_use]
    pub fn neg(&self) -> Self {
        match self {
            | Self::Int(a) => Self::Int(-a),
            | Self::Rat(a) => Self::Rat(-a),
            | Self::Float(a) => Self::Float(-a),
        }
    }

    /// Multiplicative inverse. `None` for an exact zero.
    #[must_use]
    pub fn recip(&self) -> Option<Self> {
        match self {
            | Self::Float(a) => Some(Self::Float(1.0 / a)),
            | _ => {
                let r = self.to_rational()?;
                if r.is_zero() {
                    None
                } else {
                    Some(Self::rat(r.recip()))
                }
            },
        }
    }

    /// Exponentiation.
    ///
    /// Exact when the base is exact and the exponent is a machine-sized
    /// integer, or a fraction `p/q` and the base is a positive perfect
    /// `q`-th power; a float result when either operand is a float. Returns
    /// `None` for `0^negative` and when the exact result would be an
    /// irrational algebraic number, which must stay symbolic.
    #[must_use]
    pub fn pow(
        &self,
        exp: &Self,
    ) -> Option<Self> {
        if let (Some(base), Some(e)) = (self.to_rational(), exp.to_i64()) {
            let mag = u32::try_from(e.unsigned_abs()).ok()?;
            // Refuse to build astronomically large exact powers.
            if mag > 4096 {
                return None;
            }
            if e < 0 && base.is_zero() {
                return None;
            }
            let p = num_traits::pow(base, mag as usize);
            return Some(Self::rat(if e < 0 { p.recip() } else { p }));
        }
        if let (Some(base), Self::Rat(e)) = (self.to_rational(), exp) {
            return exact_root(&base, e);
        }
        if self.is_exact() && exp.is_exact() {
            return None;
        }
        Some(Self::Float(self.to_f64().powf(exp.to_f64())))
    }

    /// Total order used for canonical sorting: by value, exact before float
    /// on ties.
    #[must_use]
    pub fn total_cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        match (self.to_rational(), other.to_rational()) {
            | (Some(a), Some(b)) => a.cmp(&b),
            | _ => self
                .to_f64()
                .total_cmp(&other.to_f64())
                .then_with(|| other.is_exact().cmp(&self.is_exact())),
        }
    }
}

/// `base^(p/q)` when `base` is a positive rational whose numerator and
/// denominator are perfect `q`-th powers.
///
/// Negative bases are refused: their principal root is not real.
fn exact_root(
    base: &BigRational,
    exp: &BigRational,
) -> Option<Number> {
    if !base.is_positive() {
        return None;
    }
    let q = exp.denom().to_u32().filter(|&q| q <= 64)?;
    let p = exp.numer().to_i64()?;
    let mag = u32::try_from(p.unsigned_abs())
        .ok()
        .filter(|&m| m <= 4096)?;
    let root = |n: &BigInt| {
        let r = n.nth_root(q);
        (num_traits::pow(r.clone(), q as usize) == *n).then_some(r)
    };
    let rooted = BigRational::new(root(base.numer())?, root(base.denom())?);
    let powered = num_traits::pow(rooted, mag as usize);
    Some(Number::rat(if p < 0 {
        powered.recip()
    } else {
        powered
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_roots_are_exact() {
        let frac = |n, d| Number::fraction(n, d).unwrap_or_else(|| Number::from(0));
        assert_eq!(Number::from(4).pow(&frac(1, 2)), Some(Number::from(2)));
        assert_eq!(Number::from(8).pow(&frac(2, 3)), Some(Number::from(4)));
        assert_eq!(frac(9, 4).pow(&frac(-1, 2)), Some(frac(2, 3)));
        assert_eq!(Number::from(2).pow(&frac(1, 2)), None);
        assert_eq!(
            Number::from(-8).pow(&frac(1, 3)),
            None,
            "principal cube root of -8 is complex"
        );
    }

    #[test]
    fn exact_arithmetic_stays_exact() {
        let third = Number::fraction(1, 3).unwrap_or_else(|| Number::from(0));
        let sum = third.add(&third).add(&third);
        assert_eq!(sum, Number::from(1));
        assert!(sum.is_integer());
    }

    #[test]
    fn float_contaminates() {
        let r = Number::from(1).add(&Number::from(0.5));
        assert_eq!(r, Number::from(1.5));
        assert!(!r.is_exact());
    }

    #[test]
    fn pow_rules() {
        let two = Number::from(2);
        assert_eq!(two.pow(&Number::from(10)), Some(Number::from(1024)));
        assert_eq!(two.pow(&Number::from(-1)), Number::fraction(1, 2));
        assert_eq!(Number::from(0).pow(&Number::from(-1)), None);
        // sqrt(2) must stay symbolic.
        assert_eq!(
            two.pow(&Number::fraction(1, 2).unwrap_or_else(|| Number::from(0))),
            None
        );
        assert_eq!(
            Number::from(4.0).pow(&Number::from(0.5)),
            Some(Number::from(2.0))
        );
    }

    #[test]
    fn recip_zero() {
        assert_eq!(Number::from(0).recip(), None);
        assert_eq!(Number::from(4).recip(), Number::fraction(1, 4));
    }
}
