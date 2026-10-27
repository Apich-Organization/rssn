//! Linear differential operators `Σ c_i(x) D^i` with rational-function
//! coefficients as elements of the Ore ring `Q(x)[D]`: composition, right
//! division with remainder, and the greatest common right divisor.
//!
//! These are the tools of the eigenring factorisation of
//! [`super::linear`]: an element `R` of the eigenring of `L` satisfies
//! `L R ≡ 0 (mod L)` (the remainder of `L ∘ R` on right division by `L`
//! vanishes), and `gcrd(L, g(R))` for a factor `g` of the minimal polynomial
//! of `R` modulo `L` is a right factor of `L`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

use super::kovacic::Rat;

type Q = BigRational;

/// An operator `Σ c[i] D^i`.
#[derive(Clone)]
pub(super) struct Ore {
    /// The coefficients of `D^0, D^1, …` (no trailing zeros).
    pub(super) c: Vec<Rat>,
}

fn zero() -> Option<Rat> {
    Rat::new(&[], &[Q::one()])
}

impl Ore {
    /// An operator from coefficients, with trailing zeros removed.
    pub(super) fn new(mut c: Vec<Rat>) -> Self {
        while c.last().is_some_and(|r| r.num.is_empty()) {
            c.pop();
        }
        Self { c }
    }

    /// From polynomial coefficients.
    pub(super) fn from_polynomials(p: &[Vec<Q>]) -> Option<Self> {
        let one = [Q::one()];
        Some(Self::new(p.iter().map(|a| Rat::new(a, &one)).collect::<Option<_>>()?))
    }

    /// The order, `None` for the zero operator.
    pub(super) const fn order(&self) -> Option<usize> {
        self.c.len().checked_sub(1)
    }

    /// `self - other`.
    pub(super) fn sub(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let n = self.c.len().max(other.c.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.c.get(i).cloned().map_or_else(zero, Some)?;
            let b = other.c.get(i).cloned().map_or_else(zero, Some)?;
            out.push(a.add(&b.scale(&-Q::one())?)?);
        }
        Some(Self::new(out))
    }

    /// `self + other`.
    pub(super) fn add(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let n = self.c.len().max(other.c.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.c.get(i).cloned().map_or_else(zero, Some)?;
            let b = other.c.get(i).cloned().map_or_else(zero, Some)?;
            out.push(a.add(&b)?);
        }
        Some(Self::new(out))
    }

    /// The operator times a rational number.
    pub(super) fn scale(
        &self,
        k: &Q,
    ) -> Option<Self> {
        Some(Self::new(self.c.iter().map(|a| a.scale(k)).collect::<Option<_>>()?))
    }

    /// The composition `self ∘ other`.
    pub(super) fn mul(
        &self,
        other: &Self,
    ) -> Option<Self> {
        if self.c.is_empty() || other.c.is_empty() {
            return Some(Self::new(Vec::new()));
        }
        let mut out: Vec<Rat> = vec![zero()?; self.c.len() + other.c.len() - 1];
        for (i, a) in self.c.iter().enumerate() {
            for (j, b) in other.c.iter().enumerate() {
                // a D^i ∘ b D^j = a Σ_k C(i, k) b^(k) D^(i - k + j).
                let mut derivative = b.clone();
                let mut binomial = BigInt::one();
                for k in 0..=i {
                    let term = a.mul(&derivative)?.scale(&Q::from_integer(binomial.clone()))?;
                    let slot = out.get_mut(i - k + j)?;
                    *slot = slot.add(&term)?;
                    derivative = derivative.derivative()?;
                    binomial = binomial * BigInt::from(i - k) / BigInt::from(k + 1);
                }
            }
        }
        Some(Self::new(out))
    }

    /// Right division: `self = quotient ∘ divisor + remainder`, the order of
    /// the remainder below that of the divisor.
    pub(super) fn right_divide(
        &self,
        divisor: &Self,
    ) -> Option<(Self, Self)> {
        let ob = divisor.order()?;
        let lead = divisor.c.last()?.clone();
        let mut remainder = self.clone();
        let mut quotient = Self::new(Vec::new());
        while let Some(or) = remainder.order() {
            if or < ob {
                break;
            }
            let e = or - ob;
            let q = remainder.c.last()?.div(&lead)?;
            let mut term = vec![zero()?; e];
            term.push(q);
            let term = Self::new(term);
            remainder = remainder.sub(&term.mul(divisor)?)?;
            quotient = quotient.add(&term)?;
            // Guard against a non-terminating cancellation failure.
            if remainder.order().is_some_and(|o| o >= or) {
                return None;
            }
        }
        Some((quotient, remainder))
    }

    /// The monic greatest common right divisor.
    pub(super) fn right_gcd(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let (mut a, mut b) = (self.clone(), other.clone());
        while b.order().is_some() {
            let (_, r) = a.right_divide(&b)?;
            a = b;
            b = r;
        }
        let lead = a.c.last()?.clone();
        let inverse = Rat::new(&[Q::one()], &[Q::one()])?.div(&lead)?;
        Some(Self::new(a.c.iter().map(|x| x.mul(&inverse)).collect::<Option<_>>()?))
    }

    /// Whether the operator is the zero operator.
    pub(super) const fn is_zero(&self) -> bool {
        self.c.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(BigInt::from(n))
    }

    fn op(coefficients: &[&[i64]]) -> Ore {
        let rows: Vec<Vec<Q>> = coefficients.iter().map(|c| c.iter().map(|&v| q(v)).collect()).collect();
        Ore::from_polynomials(&rows).expect("operator")
    }

    #[test]
    fn composition_division_and_gcrd() {
        // (D + x) ∘ (D - 1) = D^2 + (x - 1) D - x.
        let a = op(&[&[0, 1], &[1]]);
        let b = op(&[&[-1], &[1]]);
        let product = a.mul(&b).expect("product");
        let expected = op(&[&[0, -1], &[-1, 1], &[1]]);
        assert!(product.sub(&expected).expect("difference").is_zero());
        // Right division by the right factor leaves no remainder and returns the left one.
        let (quotient, remainder) = product.right_divide(&b).expect("division");
        assert!(remainder.is_zero());
        assert!(quotient.sub(&a).expect("difference").is_zero());
        // gcrd((D + x)(D - 1), (D - x)(D - 1)) = D - 1.
        let other = op(&[&[0, -1], &[1]]).mul(&b).expect("product");
        let gcd = product.right_gcd(&other).expect("gcrd");
        assert!(gcd.sub(&b).expect("difference").is_zero());
    }
}
