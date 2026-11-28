//! Sparse multivariate polynomials over the rationals.
//!
//! This is the exact counterpart of [`crate::rules::poly::repr::Poly`] used by
//! the cylindrical algebraic decomposition: coefficients are
//! [`BigRational`]s (never [`crate::graph::Number::Float`]), there are no
//! opaque generators, and the variable indices are a fixed dense range
//! `0..n`. Every operation that can grow the number of terms takes a cap and
//! returns `None` instead of exhausting memory.
//!
//! The representation is a map from monomials to non-zero coefficients. A
//! monomial is a vector of `(variable, exponent)` pairs sorted by variable,
//! with positive exponents, so the empty vector is the constant monomial
//! `1`.

// The polynomial API is exercised by the decomposition and by the unit
// tests; the handful of views not on the QE hot path are kept because they
// document the representation and back the tests.
#![allow(
    dead_code,
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

/// A monomial: `(variable, exponent)` pairs sorted by variable.
pub(crate) type Mono = Vec<(u32, u32)>;

/// `a * b` for monomials.
fn mono_mul(
    a: &Mono,
    b: &Mono,
) -> Mono {
    let mut out = Vec::with_capacity(a.len().saturating_add(b.len()));
    let (mut i, mut j) = (0, 0);
    loop {
        match (a.get(i), b.get(j)) {
            | (Some(&(ga, ea)), Some(&(gb, eb))) => {
                match ga.cmp(&gb) {
                    | std::cmp::Ordering::Equal => {
                        out.push((ga, ea.saturating_add(eb)));
                        i = i.saturating_add(1);
                        j = j.saturating_add(1);
                    },
                    | std::cmp::Ordering::Less => {
                        out.push((ga, ea));
                        i = i.saturating_add(1);
                    },
                    | std::cmp::Ordering::Greater => {
                        out.push((gb, eb));
                        j = j.saturating_add(1);
                    },
                }
            },
            | (Some(&x), None) => {
                out.push(x);
                i = i.saturating_add(1);
            },
            | (None, Some(&y)) => {
                out.push(y);
                j = j.saturating_add(1);
            },
            | (None, None) => return out,
        }
    }
}

/// Lexicographic order: the exponent of the smallest variable decides first.
fn lex_cmp(
    a: &Mono,
    b: &Mono,
) -> std::cmp::Ordering {
    let (mut i, mut j) = (0, 0);
    loop {
        match (a.get(i), b.get(j)) {
            | (Some(&(ga, ea)), Some(&(gb, eb))) => {
                match ga.cmp(&gb) {
                    | std::cmp::Ordering::Equal => {
                        if ea != eb {
                            return ea.cmp(&eb);
                        }
                        i = i.saturating_add(1);
                        j = j.saturating_add(1);
                    },
                    | std::cmp::Ordering::Less => return std::cmp::Ordering::Greater,
                    | std::cmp::Ordering::Greater => return std::cmp::Ordering::Less,
                }
            },
            | (Some(_), None) => return std::cmp::Ordering::Greater,
            | (None, Some(_)) => return std::cmp::Ordering::Less,
            | (None, None) => return std::cmp::Ordering::Equal,
        }
    }
}

/// `a / b` if the monomial `b` divides `a`.
fn mono_div(
    a: &Mono,
    b: &Mono,
) -> Option<Mono> {
    let mut out = a.clone();
    for &(g, e) in b {
        let slot = out.iter_mut().find(|(h, _)| *h == g)?;
        slot.1 = slot.1.checked_sub(e)?;
    }
    out.retain(|&(_, e)| e > 0);
    Some(out)
}

/// A sparse polynomial with rational coefficients.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct MPoly {
    terms: BTreeMap<Mono, BigRational>,
}

impl MPoly {
    /// The zero polynomial.
    #[must_use]
    pub(crate) fn zero() -> Self {
        Self::default()
    }

    /// A constant polynomial.
    #[must_use]
    pub(crate) fn constant(value: BigRational) -> Self {
        let mut p = Self::zero();
        if !value.is_zero() {
            p.terms.insert(Vec::new(), value);
        }
        p
    }

    /// The polynomial `x_index`.
    #[must_use]
    pub(crate) fn variable(index: u32) -> Self {
        Self::monomial(vec![(index, 1)], BigRational::one())
    }

    /// A single term.
    #[must_use]
    pub(crate) fn monomial(
        mono: Mono,
        coeff: BigRational,
    ) -> Self {
        let mut p = Self::zero();
        if !coeff.is_zero() {
            p.terms.insert(mono, coeff);
        }
        p
    }

    /// Whether this is the zero polynomial.
    #[must_use]
    pub(crate) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Number of terms.
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.terms.len()
    }

    /// The terms.
    pub(crate) fn terms(&self) -> impl Iterator<Item = (&Mono, &BigRational)> {
        self.terms.iter()
    }

    /// The value when the polynomial is constant.
    #[must_use]
    pub(crate) fn as_constant(&self) -> Option<BigRational> {
        match self.terms.len() {
            | 0 => Some(BigRational::zero()),
            | 1 => self.terms.get(&Vec::new()).cloned(),
            | _ => None,
        }
    }

    fn add_term(
        &mut self,
        mono: Mono,
        coeff: &BigRational,
    ) {
        match self.terms.get_mut(&mono) {
            | Some(slot) => {
                let sum = &*slot + coeff;
                if sum.is_zero() {
                    self.terms.remove(&mono);
                } else {
                    *slot = sum;
                }
            },
            | None => {
                if !coeff.is_zero() {
                    self.terms.insert(mono, coeff.clone());
                }
            },
        }
    }

    /// Sum.
    #[must_use]
    pub(crate) fn add(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = self.clone();
        for (mono, coeff) in &other.terms {
            out.add_term(mono.clone(), coeff);
        }
        out
    }

    /// Additive inverse.
    #[must_use]
    pub(crate) fn neg(&self) -> Self {
        Self {
            terms: self.terms.iter().map(|(m, c)| (m.clone(), -c)).collect(),
        }
    }

    /// Difference.
    #[must_use]
    pub(crate) fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }

    /// Product with a rational.
    #[must_use]
    pub(crate) fn scale(
        &self,
        factor: &BigRational,
    ) -> Self {
        if factor.is_zero() {
            return Self::zero();
        }
        Self {
            terms: self
                .terms
                .iter()
                .map(|(m, c)| (m.clone(), c * factor))
                .collect(),
        }
    }

    /// Product. `None` when the result would exceed `cap` terms.
    #[must_use]
    pub(crate) fn mul(
        &self,
        other: &Self,
        cap: usize,
    ) -> Option<Self> {
        let mut out = Self::zero();
        for (ma, ca) in &self.terms {
            for (mb, cb) in &other.terms {
                out.add_term(mono_mul(ma, mb), &(ca * cb));
            }
            if out.terms.len() > cap {
                return None;
            }
        }
        Some(out)
    }

    /// Power by repeated squaring. `None` when a product exceeds `cap`.
    #[must_use]
    pub(crate) fn pow(
        &self,
        mut exponent: u32,
        cap: usize,
    ) -> Option<Self> {
        let mut result = Self::constant(BigRational::one());
        let mut base = self.clone();
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.mul(&base, cap)?;
            }
            exponent >>= 1;
            if exponent > 0 {
                base = base.mul(&base, cap)?;
            }
        }
        Some(result)
    }

    /// The leading term in the lexicographic order of the variables.
    #[must_use]
    pub(crate) fn leading(&self) -> Option<(&Mono, &BigRational)> {
        self.terms.iter().max_by(|a, b| lex_cmp(a.0, b.0))
    }

    /// The coefficient of the numerically largest monomial.
    #[must_use]
    pub(crate) fn leading_coefficient(&self) -> Option<&BigRational> {
        self.leading().map(|(_, c)| c)
    }

    /// The polynomial with a positive leading coefficient (the sign of a
    /// projection polynomial is irrelevant).
    #[must_use]
    pub(crate) fn positive_leading(&self) -> Self {
        match self.leading_coefficient() {
            | Some(c) if c.is_negative() => self.neg(),
            | _ => self.clone(),
        }
    }

    /// Degree in `variable`; zero for the zero polynomial.
    #[must_use]
    pub(crate) fn degree_in(
        &self,
        variable: u32,
    ) -> u32 {
        self.terms
            .keys()
            .filter_map(|m| m.iter().find(|&&(g, _)| g == variable).map(|&(_, e)| e))
            .max()
            .unwrap_or(0)
    }

    /// Total degree; zero for the zero polynomial.
    #[must_use]
    pub(crate) fn total_degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|m| m.iter().fold(0_u32, |acc, &(_, e)| acc.saturating_add(e)))
            .max()
            .unwrap_or(0)
    }

    /// The variables that actually occur.
    #[must_use]
    pub(crate) fn support(&self) -> Vec<u32> {
        let mut vars: Vec<u32> = self
            .terms
            .keys()
            .flat_map(|m| m.iter().map(|&(g, _)| g))
            .collect();
        vars.sort_unstable();
        vars.dedup();
        vars
    }

    /// Views the polynomial as univariate in `variable`: entry `k` is the
    /// coefficient of `variable^k`.
    #[must_use]
    pub(crate) fn coefficients_in(
        &self,
        variable: u32,
    ) -> Vec<Self> {
        let degree = usize::try_from(self.degree_in(variable)).unwrap_or(usize::MAX);
        let mut out = vec![Self::zero(); degree.saturating_add(1)];
        for (mono, coeff) in &self.terms {
            let power = mono
                .iter()
                .find(|&&(g, _)| g == variable)
                .map_or(0, |&(_, e)| e);
            let rest: Mono = mono
                .iter()
                .copied()
                .filter(|&(g, _)| g != variable)
                .collect();
            if let Some(slot) = out.get_mut(usize::try_from(power).unwrap_or(usize::MAX)) {
                slot.add_term(rest, coeff);
            }
        }
        out
    }

    /// The univariate coefficients when only `variable` occurs.
    #[must_use]
    pub(crate) fn univariate_in(
        &self,
        variable: u32,
    ) -> Option<Vec<BigRational>> {
        self.coefficients_in(variable)
            .iter()
            .map(Self::as_constant)
            .collect()
    }

    /// Builds a polynomial in `variable` from coefficients by degree.
    #[must_use]
    pub(crate) fn from_univariate(
        variable: u32,
        coeffs: &[BigRational],
    ) -> Self {
        let mut out = Self::zero();
        for (power, coeff) in coeffs.iter().enumerate() {
            let exponent = u32::try_from(power).unwrap_or(u32::MAX);
            let mono = if exponent == 0 {
                Vec::new()
            } else {
                vec![(variable, exponent)]
            };
            out.add_term(mono, coeff);
        }
        out
    }

    /// Partial derivative with respect to `variable`.
    #[must_use]
    pub(crate) fn derivative(
        &self,
        variable: u32,
    ) -> Self {
        let mut out = Self::zero();
        for (mono, coeff) in &self.terms {
            let Some(&(_, power)) = mono.iter().find(|&&(g, _)| g == variable) else {
                continue;
            };
            if power == 0 {
                continue;
            }
            let lowered: Mono = mono
                .iter()
                .filter_map(|&(g, e)| {
                    if g == variable {
                        (e > 1).then_some((g, e.saturating_sub(1)))
                    } else {
                        Some((g, e))
                    }
                })
                .collect();
            out.add_term(
                lowered,
                &(coeff * BigRational::from_integer(BigInt::from(power))),
            );
        }
        out
    }

    /// The exact quotient `self / divisor`, or `None` if `divisor` does not
    /// divide `self` or the division needs more than `cap` steps.
    #[must_use]
    pub(crate) fn div_exact(
        &self,
        divisor: &Self,
        cap: usize,
    ) -> Option<Self> {
        let (lead_mono, lead_coeff) = divisor.leading()?;
        if lead_coeff.is_zero() {
            return None;
        }
        let (lead_mono, lead_coeff) = (lead_mono.clone(), lead_coeff.clone());
        let mut quotient = Self::zero();
        let mut rest = self.clone();
        for _ in 0..cap {
            let Some((mono, coeff)) = rest.leading() else {
                return Some(quotient);
            };
            let (mono, coeff) = (mono.clone(), coeff.clone());
            let factor = mono_div(&mono, &lead_mono)?;
            let term = Self::monomial(factor, &coeff / &lead_coeff);
            rest = rest.sub(&term.mul(divisor, cap.saturating_mul(4))?);
            quotient = quotient.add(&term);
        }
        rest.is_zero().then_some(quotient)
    }

    /// Evaluates the polynomial at rational values, in variable order.
    #[must_use]
    pub(crate) fn eval_rational(
        &self,
        values: &[BigRational],
    ) -> Option<BigRational> {
        let mut total = BigRational::zero();
        for (mono, coeff) in &self.terms {
            let mut term = coeff.clone();
            for &(variable, exponent) in mono {
                let value = values.get(usize::try_from(variable).ok()?)?;
                term *= num_traits::pow(value.clone(), usize::try_from(exponent).ok()?);
            }
            total += term;
        }
        Some(total)
    }

    /// Substitutes the rational values in `values` for their variables and
    /// leaves the others symbolic.
    #[must_use]
    pub(crate) fn substitute_rational(
        &self,
        values: &BTreeMap<u32, BigRational>,
    ) -> Self {
        let mut out = Self::zero();
        for (mono, coeff) in &self.terms {
            let mut c = coeff.clone();
            let mut rest: Mono = Vec::with_capacity(mono.len());
            for &(variable, exponent) in mono {
                match values.get(&variable) {
                    | Some(value) => {
                        c *= num_traits::pow(
                            value.clone(),
                            usize::try_from(exponent).unwrap_or(usize::MAX),
                        );
                    },
                    | None => rest.push((variable, exponent)),
                }
            }
            out.add_term(rest, &c);
        }
        out
    }

    /// Applies a variable renaming: `map[i]` is the new index of variable
    /// `i`. Variables not listed are kept.
    #[must_use]
    pub(crate) fn renamed(
        &self,
        map: &[u32],
    ) -> Self {
        let mut out = Self::zero();
        for (mono, coeff) in &self.terms {
            let mapped: Mono = mono
                .iter()
                .map(|&(g, e)| {
                    (
                        map.get(usize::try_from(g).unwrap_or(usize::MAX))
                            .copied()
                            .unwrap_or(g),
                        e,
                    )
                })
                .collect();
            let mut sorted = mapped;
            sorted.sort_unstable_by_key(|&(g, _)| g);
            out.add_term(sorted, coeff);
        }
        out
    }

    /// The polynomial with every variable shifted by `offset`.
    #[must_use]
    pub(crate) fn shifted_variables(
        &self,
        offset: u32,
    ) -> Self {
        let mut out = Self::zero();
        for (mono, coeff) in &self.terms {
            let mapped: Mono = mono
                .iter()
                .map(|&(g, e)| (g.saturating_add(offset), e))
                .collect();
            out.add_term(mapped, coeff);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(
        n: i64,
        d: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(n), BigInt::from(d))
    }

    #[test]
    fn arithmetic() {
        let x = MPoly::variable(0);
        let y = MPoly::variable(1);
        let one = MPoly::constant(BigRational::one());
        let two = MPoly::constant(BigRational::from_integer(BigInt::from(2)));
        let a = x.add(&one);
        let b = y.sub(&x);
        let cap = 1_000;
        let left = a.mul(&b, cap).and_then(|p| p.mul(&a, cap));
        let right = a.pow(2, cap).and_then(|p| p.mul(&b, cap));
        assert_eq!(left, right);
        assert!(a.sub(&a).is_zero());
        assert_eq!(
            x.add(&one).pow(2, cap),
            Some(
                x.pow(2, cap)
                    .unwrap_or_default()
                    .add(&two.mul(&x, cap).unwrap_or_default())
                    .add(&one)
            )
        );
        assert_eq!(a.scale(&q(1, 2)).eval_rational(&[q(1, 2)]), Some(q(3, 4)));
    }

    #[test]
    fn univariate_views() {
        let x = MPoly::variable(0);
        let y = MPoly::variable(1);
        let p = x
            .pow(2, 100)
            .unwrap_or_default()
            .mul(&y, 100)
            .unwrap_or_default()
            .add(
                &x.mul(
                    &MPoly::constant(BigRational::from_integer(BigInt::from(5))),
                    100,
                )
                .unwrap_or_default(),
            );
        assert_eq!(p.degree_in(0), 2);
        assert_eq!(p.degree_in(1), 1);
        assert_eq!(p.total_degree(), 3);
        let coeffs = p.coefficients_in(0);
        assert_eq!(coeffs.len(), 3);
        assert_eq!(
            coeffs.first().and_then(MPoly::as_constant),
            Some(BigRational::zero())
        );
        assert_eq!(
            coeffs.get(1).and_then(MPoly::as_constant),
            Some(BigRational::from_integer(BigInt::from(5)))
        );
        assert_eq!(
            coeffs.get(2).and_then(|c| c.univariate_in(1)),
            Some(vec![BigRational::zero(), BigRational::one()])
        );
        assert_eq!(p.support(), vec![0, 1]);
    }

    #[test]
    fn exact_division() {
        let x = MPoly::variable(0);
        let y = MPoly::variable(1);
        let one = MPoly::constant(BigRational::one());
        let a = x.add(&y);
        let b = x.sub(&y);
        let product = a.mul(&b, 1_000).unwrap_or_default();
        assert_eq!(product.div_exact(&a, 1_000), Some(b));
        assert_eq!(product.div_exact(&x, 1_000), None);
        assert_eq!(one.div_exact(&one, 10), Some(one));
    }

    #[test]
    fn derivatives() {
        let x = MPoly::variable(0);
        let y = MPoly::variable(1);
        let p = x
            .pow(2, 100)
            .unwrap_or_default()
            .mul(&y, 100)
            .unwrap_or_default()
            .add(
                &x.mul(&y.pow(2, 100).unwrap_or_default(), 100)
                    .unwrap_or_default(),
            );
        let dx = p.derivative(0);
        assert_eq!(dx.degree_in(0), 1);
        assert_eq!(dx.degree_in(1), 2);
        let dy = p.derivative(1);
        assert_eq!(dy.degree_in(0), 2);
        assert_eq!(dy.degree_in(1), 1);
        assert!(p.derivative(2).is_zero());
    }

    #[test]
    fn limits_are_enforced() {
        let x = MPoly::variable(0);
        let y = MPoly::variable(1);
        let z = MPoly::variable(2);
        let sum = x.add(&y).add(&z);
        assert!(sum.pow(40, 100).is_none());
        assert!(sum.pow(2, 100).is_some());
    }
}
