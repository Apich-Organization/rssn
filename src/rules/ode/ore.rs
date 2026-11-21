//! Linear differential operators `Σ c_i(x) D^i` with rational-function
//! coefficients as elements of the Ore ring `Q(x)[D]`: composition, right
//! division with remainder, and the greatest common right divisor.
//!
//! These are the tools of the eigenring factorisation of
//! [`super::linear`]: an element `R` of the eigenring of `L` satisfies
//! `L R ≡ 0 (mod L)` (the remainder of `L ∘ R` on right division by `L`
//! vanishes), and `gcrd(L, g(R))` for a factor `g` of the minimal polynomial
//! of `R` modulo `L` is a right factor of `L`.
//!
//! Composition and right division reuse derivatives of the right operand's
//! coefficients. Large sets of derivative ladders can use the active CPU pool.

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

/// Right coefficient derivatives up to `order`, stopping at the first zero.
fn derivative_ladders(
    coefficients: &[Rat],
    order: usize,
) -> Option<Vec<Vec<Rat>>> {
    let ladder = |coefficient: &Rat| {
        if coefficient.num.is_empty() {
            return Some(Vec::new());
        }
        let mut derivatives = vec![coefficient.clone()];
        while derivatives.len() <= order {
            let next = derivatives.last()?.derivative()?;
            if next.num.is_empty() {
                break;
            }
            derivatives.push(next);
        }
        Some(derivatives)
    };
    if order >= 4
        && coefficients
            .iter()
            .filter(|coefficient| {
                !coefficient.num.is_empty()
                    && (coefficient.num.len() > 1 || coefficient.den.len() > 1)
            })
            .take(4)
            .count()
            == 4
    {
        crate::kernels::par::map_slice(coefficients, 4, ladder)
            .into_iter()
            .collect()
    } else {
        coefficients.iter().map(ladder).collect()
    }
}

/// Adds `coefficient * D^degree` composed on the right with the prepared operator.
fn add_monomial_product(
    out: &mut [Rat],
    coefficient: &Rat,
    degree: usize,
    right: &[Vec<Rat>],
) -> Option<()> {
    for (j, derivatives) in right.iter().enumerate() {
        // a D^i ∘ b D^j = a Σ_k C(i, k) b^(k) D^(i - k + j).
        let mut binomial = BigInt::one();
        for (k, derivative) in derivatives.iter().take(degree + 1).enumerate() {
            let mut term = coefficient.mul(derivative)?;
            if !binomial.is_one() {
                term = term.scale(&Q::from_integer(binomial.clone()))?;
            }
            let slot = out.get_mut(degree - k + j)?;
            *slot = if slot.num.is_empty() {
                term
            } else {
                slot.add(&term)?
            };
            if k == degree {
                break;
            }
            binomial = binomial * BigInt::from(degree - k) / BigInt::from(k + 1);
        }
    }
    Some(())
}

impl Ore {
    /// An operator from coefficients, with trailing zeros removed.
    pub(super) fn new(c: Vec<Rat>) -> Self {
        let mut operator = Self { c };
        operator.trim();
        operator
    }

    fn trim(&mut self) {
        while self.c.last().is_some_and(|r| r.num.is_empty()) {
            self.c.pop();
        }
    }

    /// From polynomial coefficients.
    pub(super) fn from_polynomials(p: &[Vec<Q>]) -> Option<Self> {
        let one = [Q::one()];
        Some(Self::new(
            p.iter().map(|a| Rat::new(a, &one)).collect::<Option<_>>()?,
        ))
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
        Some(Self::new(
            self.c.iter().map(|a| a.scale(k)).collect::<Option<_>>()?,
        ))
    }

    /// The composition `self ∘ other`.
    pub(super) fn mul(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let Some(order) = self.order().filter(|_| !other.is_zero()) else {
            return Some(Self::new(Vec::new()));
        };
        let right = derivative_ladders(&other.c, order)?;
        let mut out: Vec<Rat> = vec![zero()?; self.c.len() + other.c.len() - 1];
        for (degree, coefficient) in self.c.iter().enumerate().filter(|(_, a)| !a.num.is_empty()) {
            add_monomial_product(&mut out, coefficient, degree, &right)?;
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
        let Some(max_degree) = self.order().and_then(|order| order.checked_sub(ob)) else {
            return Some((Self::new(Vec::new()), self.clone()));
        };
        let lead = divisor.c.last()?;
        let right = derivative_ladders(&divisor.c, max_degree)?;
        let negative_one = -Q::one();
        let mut remainder = self.clone();
        let mut quotient = vec![zero()?; max_degree + 1];
        while let Some(or) = remainder.order() {
            if or < ob {
                break;
            }
            let e = or - ob;
            let q = remainder.c.last()?.div(lead)?;
            add_monomial_product(&mut remainder.c, &q.scale(&negative_one)?, e, &right)?;
            remainder.trim();
            // The remainder order strictly decreases, so each quotient coefficient is assigned once.
            *quotient.get_mut(e)? = q;
            // Guard against a non-terminating cancellation failure.
            if remainder.order().is_some_and(|o| o >= or) {
                return None;
            }
        }
        Some((Self::new(quotient), remainder))
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
        Some(Self::new(
            a.c.iter().map(|x| x.mul(&inverse)).collect::<Option<_>>()?,
        ))
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
        let rows: Vec<Vec<Q>> = coefficients
            .iter()
            .map(|c| c.iter().map(|&v| q(v)).collect())
            .collect();
        Ore::from_polynomials(&rows).expect("operator")
    }

    fn rat(
        numerator: &[i64],
        denominator: &[i64],
    ) -> Rat {
        let rationals =
            |coefficients: &[i64]| coefficients.iter().copied().map(q).collect::<Vec<_>>();
        Rat::new(&rationals(numerator), &rationals(denominator)).expect("rational function")
    }

    fn equal(
        left: &Ore,
        right: &Ore,
    ) {
        assert_eq!(left.c.len(), right.c.len());
        for (a, b) in left.c.iter().zip(&right.c) {
            assert_eq!(a.num, b.num);
            assert_eq!(a.den, b.den);
        }
    }

    fn apply(
        operator: &Ore,
        input: &Rat,
    ) -> Rat {
        let mut value = zero().unwrap();
        let mut differentiated = input.clone();
        for coefficient in &operator.c {
            value = value
                .add(&coefficient.mul(&differentiated).unwrap())
                .unwrap();
            differentiated = differentiated.derivative().unwrap();
        }
        value
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

    #[test]
    fn leibniz_composition_is_noncommutative_and_truncates_zero_derivatives() {
        let d4 = op(&[&[], &[], &[], &[], &[1]]);
        let cube = op(&[&[0, 0, 0, 1]]);
        let expected = op(&[&[], &[24], &[0, 36], &[0, 0, 12], &[0, 0, 0, 1]]);
        equal(&d4.mul(&cube).unwrap(), &expected);
        equal(
            &cube.mul(&d4).unwrap(),
            &op(&[&[], &[], &[], &[], &[0, 0, 0, 1]]),
        );
    }

    #[test]
    fn rational_composition_agrees_with_separate_operator_actions() {
        let left = Ore::new(vec![
            rat(&[1, 2], &[1, 1]),
            zero().unwrap(),
            rat(&[-1, 1], &[1]),
        ]);
        let right = Ore::new(vec![rat(&[1], &[2, 1]), rat(&[2, 0, 1], &[1, 0, 1])]);
        let product = left.mul(&right).unwrap();
        for input in [rat(&[1, 0, 2], &[1]), rat(&[-2, 1], &[1, 1])] {
            let direct = apply(&product, &input);
            let sequential = apply(&left, &apply(&right, &input));
            assert_eq!(direct.num, sequential.num);
            assert_eq!(direct.den, sequential.den);
        }
    }

    #[test]
    fn right_division_recovers_a_sparse_quotient_and_small_remainder() {
        let divisor = Ore::new(vec![
            rat(&[1], &[2, 1]),
            zero().unwrap(),
            rat(&[1, 1], &[1]),
        ]);
        let mut terms = vec![zero().unwrap(); 7];
        terms[1] = rat(&[0, 1], &[1]);
        terms[6] = rat(&[3], &[1, 1]);
        let quotient = Ore::new(terms);
        let remainder = Ore::new(vec![rat(&[1], &[0, 1]), rat(&[2, -1], &[1])]);
        let dividend = quotient.mul(&divisor).unwrap().add(&remainder).unwrap();
        let (computed, rest) = dividend.right_divide(&divisor).unwrap();
        equal(&computed, &quotient);
        equal(&rest, &remainder);
        let (computed, rest) = remainder.right_divide(&divisor).unwrap();
        assert!(computed.is_zero());
        equal(&rest, &remainder);
        let zero = op(&[]);
        assert!(dividend.right_divide(&zero).is_none());
        let (computed, rest) = zero.right_divide(&divisor).unwrap();
        assert!(computed.is_zero() && rest.is_zero());
    }

    #[test]
    fn prepared_derivatives_match_across_worker_counts() {
        let left = op(&[&[1], &[], &[0, 1], &[], &[], &[1]]);
        let right = op(&[&[0, 1], &[1, 0, 1], &[], &[-1, 1], &[1, 0, 0, 1]]);
        let run = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let product = left.mul(&right).unwrap();
                    let (quotient, remainder) = product.right_divide(&right).unwrap();
                    equal(&quotient, &left);
                    assert!(remainder.is_zero());
                    product
                })
        };
        equal(&run(1), &run(4));
    }

    #[test]
    fn right_gcd_retains_a_shared_rational_right_factor() {
        let factor = Ore::new(vec![rat(&[-1], &[0, 1]), rat(&[1], &[1])]);
        let first = op(&[&[0, 1], &[1]]).mul(&factor).unwrap();
        let second = op(&[&[0, -1], &[1]]).mul(&factor).unwrap();
        equal(&first.right_gcd(&second).unwrap(), &factor);
        equal(&second.right_gcd(&first).unwrap(), &factor);
        equal(&factor.right_gcd(&op(&[])).unwrap(), &factor);
        assert!(op(&[]).right_gcd(&op(&[])).is_none());
    }
}
