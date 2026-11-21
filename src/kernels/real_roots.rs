//! Numerical Sturm isolation and refinement for real roots of a polynomial.
//!
//! [`find_roots`] uses a Cauchy bound for the complete search. Use
//! [`find_roots_in`] when a finite interval is already known. Both use a
//! floating Sturm chain for the interior of each interval.

use crate::kernels::polynomial::Polynomial;
use num_rational::BigRational;
use num_traits::Zero;

/// Generates the Sturm sequence for a numerical polynomial.
///
/// The Sturm sequence is a sequence of polynomials derived from the original polynomial
/// and its derivative. It is used to determine the number of distinct real roots
/// of a polynomial in a given interval.
///
/// # Arguments
/// * `poly` - The input polynomial as a `Polynomial`.
///
/// # Returns
/// A `Vec<Polynomial>` representing the Sturm sequence.
#[must_use]
pub fn sturm_sequence(poly: &Polynomial) -> Vec<Polynomial> {
    let mut seq = Vec::new();

    if poly.coeffs.is_empty() || (poly.coeffs.len() == 1 && poly.coeffs[0] == 0.0) {
        return seq;
    }

    seq.push(poly.clone());

    let p1 = poly.derivative();

    if p1.coeffs.is_empty() || (p1.coeffs.len() == 1 && p1.coeffs[0] == 0.0) {
        return seq;
    }

    seq.push(p1);

    let mut i = 1;

    while seq[i].coeffs.len() > 1 || seq[i].coeffs[0] != 0.0 {
        let p_prev = &seq[i - 1];

        let p_curr = &seq[i];

        let (_, mut remainder) = p_prev.clone().long_division(p_curr);

        if remainder.coeffs.is_empty()
            || (remainder.coeffs.len() == 1 && remainder.coeffs[0] == 0.0)
        {
            break;
        }

        for c in &mut remainder.coeffs {
            *c = -*c;
        }

        seq.push(remainder);

        i += 1;
    }

    seq
}

// Count distinct roots using a square-free, positively rescaled chain. The
// public Sturm sequence retains its original coefficient conventions; removing
// its common polynomial factor prevents cancellation in repeated-root powers
// from moving a sign change by sqrt(machine epsilon).
fn counting_sequence(poly: &Polynomial) -> Result<Vec<Polynomial>, String> {
    // Binary rescaling preserves coefficient ratios unless a coefficient
    // underflows. Arbitrary positive scaling can move an exact split root.
    let power_of_two = |magnitude: f64| {
        let exponent_bits = magnitude.to_bits() & 0x7ff0_0000_0000_0000;
        if exponent_bits == 0 {
            f64::MIN_POSITIVE
        } else {
            f64::from_bits(exponent_bits)
        }
    };
    let scale = poly
        .coeffs
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if !scale.is_finite() || super::exact_eq(scale, 0.0) {
        return Ok(Vec::new());
    }
    let scale = power_of_two(scale);
    let normalized = Polynomial::new(poly.coeffs.iter().map(|value| value / scale).collect());
    if super::exact_eq(normalized.coeffs[0], 0.0) {
        return Err("polynomial leading coefficient underflows during Sturm scaling".into());
    }
    let mut sequence = sturm_sequence(&normalized);
    if let Some(gcd) = sequence.last().filter(|p| p.coeffs.len() > 1).cloned() {
        sequence = sequence
            .into_iter()
            .map(|p| p.long_division(&gcd).0)
            .collect();
    }
    for polynomial in &mut sequence {
        let scale = polynomial
            .coeffs
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0, f64::max);
        if scale > 0.0 && scale.is_finite() {
            let scale = power_of_two(scale);
            for coefficient in &mut polynomial.coeffs {
                *coefficient /= scale;
            }
        }
    }
    if sequence.iter().any(|p| {
        p.coeffs.first().is_none_or(|&c| super::exact_eq(c, 0.0))
            || p.coeffs.iter().any(|c| !c.is_finite())
    }) {
        return Err("polynomial Sturm chain is not representable".into());
    }
    Ok(sequence)
}

struct SturmTerm {
    polynomial: Polynomial,
    zero_order: usize,
}

impl SturmTerm {
    fn new(mut polynomial: Polynomial) -> Self {
        let retained = polynomial
            .coeffs
            .iter()
            .rposition(|&c| !super::exact_eq(c, 0.0))
            .map_or(0, |last| last + 1);
        let zero_order = polynomial.coeffs.len() - retained;
        polynomial.coeffs.truncate(retained);
        Self {
            polynomial,
            zero_order,
        }
    }

    fn sign_at(
        &self,
        x: f64,
    ) -> i8 {
        if self.zero_order != 0 && super::exact_eq(x, 0.0) {
            return 0;
        }
        // Evaluate the other factor without multiplying by a power of x.
        // This preserves the sign next to zero even at subnormal x.
        let value = self.polynomial.eval(x);
        let sign = i8::from(value > 0.0) - i8::from(value < 0.0);
        if x < 0.0 && self.zero_order % 2 == 1 {
            -sign
        } else {
            sign
        }
    }
}

#[derive(Clone, Copy)]
struct Endpoint {
    x: f64,
    variations: usize,
    is_root: bool,
}

#[derive(Clone, Copy)]
struct Bracket {
    left: Endpoint,
    right: Endpoint,
}

impl Bracket {
    fn count(self) -> Result<usize, String> {
        self.left
            .variations
            .checked_sub(self.right.variations)
            .ok_or_else(|| "polynomial Sturm variations increased over an interval".into())
    }

    const fn midpoint(self) -> f64 {
        f64::midpoint(self.left.x, self.right.x)
    }

    fn width(self) -> f64 {
        self.right.x - self.left.x
    }
}

/// Certifies zero at a caller-supplied floating endpoint using the original
/// binary coefficients. Well-separated endpoints need only floating Horner.
fn exact_endpoint_root(
    poly: &Polynomial,
    x: f64,
) -> bool {
    let Some(&leading) = poly.coeffs.first() else {
        return false;
    };
    if super::exact_eq(x, 0.0) {
        return poly.coeffs.last().is_some_and(|&c| super::exact_eq(c, 0.0));
    }
    // The positive Horner recurrence bounds the magnitude of the monomials.
    // Leave extra room for accumulation and subnormal rounding. A subnormal
    // leading coefficient needs the exact path because its initial rounding
    // error could be magnified by a large endpoint.
    let (value, magnitude) = poly.coeffs.iter().fold((0.0_f64, 0.0_f64), |(v, m), &c| {
        (v.mul_add(x, c), m.mul_add(x.abs(), c.abs()))
    });
    let terms = poly.coeffs.len() as f64;
    let uncertainty = terms * (8.0 * f64::EPSILON) * magnitude + terms * f64::from_bits(1);
    if leading.is_normal() && value.is_finite() && value.abs() > uncertainty {
        return false;
    }
    if !magnitude.is_finite()
        && root_bound(poly).is_ok_and(|bound| x.abs() > bound.next_up().next_up())
    {
        return false;
    }
    let Some(argument) = BigRational::from_float(x) else {
        return false;
    };
    poly.coeffs
        .iter()
        .try_fold(BigRational::zero(), |value, &coefficient| {
            BigRational::from_float(coefficient).map(|coefficient| value * &argument + coefficient)
        })
        .is_some_and(|value| value.is_zero())
}

struct SturmCounter<'a> {
    original: &'a Polynomial,
    terms: Vec<SturmTerm>,
}

fn validate_coefficients(poly: &Polynomial) -> Result<(), String> {
    if poly
        .coeffs
        .iter()
        .any(|coefficient| !coefficient.is_finite())
    {
        return Err("polynomial coefficients must be finite".into());
    }
    if poly
        .coeffs
        .first()
        .is_some_and(|&c| super::exact_eq(c, 0.0))
    {
        return Err("Leading coefficient cannot be zero.".into());
    }
    Ok(())
}

impl<'a> SturmCounter<'a> {
    fn new(poly: &'a Polynomial) -> Result<Self, String> {
        validate_coefficients(poly)?;
        let terms = counting_sequence(poly)?
            .into_iter()
            .map(SturmTerm::new)
            .collect();
        Ok(Self {
            original: poly,
            terms,
        })
    }

    fn endpoint(
        &self,
        x: f64,
        exact_root: Option<bool>,
    ) -> Endpoint {
        let Some((first, rest)) = self.terms.split_first() else {
            return Endpoint {
                x,
                variations: 0,
                is_root: false,
            };
        };
        let mut previous = if exact_root == Some(true) {
            0
        } else {
            first.sign_at(x)
        };
        let is_root = exact_root.unwrap_or(previous == 0);
        let mut variations = 0;
        for term in rest {
            let sign = term.sign_at(x);
            if sign != 0 {
                variations += usize::from(previous != 0 && sign != previous);
                previous = sign;
            }
        }
        Endpoint {
            x,
            variations,
            is_root,
        }
    }

    fn certified_endpoint(
        &self,
        x: f64,
    ) -> Endpoint {
        self.endpoint(x, Some(exact_endpoint_root(self.original, x)))
    }

    fn bracket(
        &self,
        (left, right): (f64, f64),
    ) -> Bracket {
        let left = self.certified_endpoint(left);
        let right = if super::exact_eq(left.x, right) {
            left
        } else {
            self.certified_endpoint(right)
        };
        Bracket { left, right }
    }

    // Sturm counts exclude the lower endpoint and include the upper one.
    fn visit_isolated(
        &self,
        bracket: Bracket,
        precision: f64,
        mut accept: impl FnMut(Bracket) -> Result<(), String>,
    ) -> Result<(), String> {
        let mut pending = vec![bracket];
        let mut subdivisions = 0;
        while let Some(bracket) = pending.pop() {
            let count = bracket.count()?;
            if count == 0 {
                continue;
            }
            if count == 1 || bracket.width() < precision {
                accept(bracket)?;
                continue;
            }
            let midpoint = bracket.midpoint();
            if super::exact_eq(midpoint, bracket.left.x)
                || super::exact_eq(midpoint, bracket.right.x)
            {
                return Err("roots cannot be separated at requested floating precision".into());
            }
            subdivisions += 1;
            if subdivisions > 131072 {
                return Err("root isolation work budget exhausted".into());
            }
            let midpoint = self.certified_endpoint(midpoint);
            pending.push(Bracket {
                left: midpoint,
                right: bracket.right,
            });
            pending.push(Bracket {
                left: bracket.left,
                right: midpoint,
            });
        }
        Ok(())
    }

    fn closed_interval(
        &self,
        bracket: Bracket,
    ) -> Result<(f64, f64), String> {
        let mut left = bracket.left;
        if left.is_root {
            let successor = left.x.next_up();
            if successor > bracket.right.x {
                return Err("roots cannot be enclosed by separate floating intervals".into());
            }
            let successor = self.certified_endpoint(successor);
            let skipped = Bracket {
                left,
                right: successor,
            }
            .count()?;
            if skipped > usize::from(successor.is_root) {
                return Err("roots cannot be enclosed by separate floating intervals".into());
            }
            left = successor;
        }
        Ok((left.x, bracket.right.x))
    }

    fn refine(
        &self,
        mut bracket: Bracket,
        tolerance: f64,
    ) -> Result<f64, String> {
        if bracket.count()? == 0 {
            return Err("the root interval contains no root".into());
        }
        if bracket.right.is_root {
            return Ok(bracket.right.x);
        }
        // Counts, unlike residual signs, also locate roots of even multiplicity.
        for _ in 0..4096 {
            let midpoint = bracket.midpoint();
            if super::exact_eq(midpoint, bracket.left.x) {
                return Ok(bracket.right.x);
            }
            if bracket.width() <= tolerance || super::exact_eq(midpoint, bracket.right.x) {
                return Ok(midpoint);
            }
            let endpoint = self.endpoint(midpoint, None);
            let lower = Bracket {
                left: bracket.left,
                right: endpoint,
            };
            let upper = Bracket {
                left: endpoint,
                right: bracket.right,
            };
            let lower_count = lower.count()?;
            let upper_count = upper.count()?;
            if lower_count != 0 {
                bracket = lower;
            } else if upper_count != 0 {
                bracket = upper;
            } else {
                return Err("polynomial Sturm count was lost during refinement".into());
            }
        }
        Err("root refinement work budget exhausted".into())
    }

    fn find(
        &self,
        bracket: Bracket,
        tolerance: f64,
    ) -> Result<Vec<f64>, String> {
        let mut roots = Vec::new();
        if bracket.left.is_root {
            roots.push(bracket.left.x);
        }
        self.visit_isolated(bracket, tolerance, |interval| {
            roots.push(self.refine(interval, tolerance)?);
            Ok(())
        })?;
        Ok(roots)
    }
}

fn valid_precision(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn valid_interval((left, right): (f64, f64)) -> bool {
    left.is_finite() && right.is_finite() && left <= right
}

fn whole_interval(poly: &Polynomial) -> Result<(f64, f64), String> {
    let bound = root_bound(poly)?;
    if !bound.is_finite() {
        return Err("polynomial root bound is not representable".into());
    }
    Ok((-bound, bound))
}

/// Isolates the numerical real roots of a polynomial in ascending order.
///
/// Intervals are closed. A bracket may be wider than `precision` if its
/// Sturm count is one. Roots that stay together in a bracket narrower than
/// `precision` are returned as one cluster interval.
///
/// # Errors
/// Invalid coefficients or precision, an unrepresentable root bound or
/// Sturm chain, or failure to separate adjacent brackets in floating point.
pub fn isolate_real_roots(
    poly: &Polynomial,
    precision: f64,
) -> Result<Vec<(f64, f64)>, String> {
    if !valid_precision(precision) {
        return Err("root isolation precision must be finite and positive".into());
    }
    let interval = whole_interval(poly)?;
    let counter = SturmCounter::new(poly)?;
    let bracket = counter.bracket(interval);
    let mut intervals = Vec::new();
    if bracket.left.is_root {
        intervals.push((bracket.left.x, bracket.left.x));
    }
    counter.visit_isolated(bracket, precision, |interval| {
        intervals.push(counter.closed_interval(interval)?);
        Ok(())
    })?;
    Ok(intervals)
}

/// Computes an upper bound for the absolute value of the real roots (Cauchy's bound).
pub(crate) fn root_bound(poly: &Polynomial) -> Result<f64, String> {
    if poly.coeffs.is_empty() {
        return Ok(1.0);
    }

    let lc = poly.coeffs[0];

    if lc == 0.0 {
        return Err("Leading coefficient \
             cannot be zero."
            .to_string());
    }

    let max_coeff = poly
        .coeffs
        .iter()
        .skip(1)
        .map(|c| c.abs())
        .fold(0.0, f64::max);

    Ok(1.0 + max_coeff / lc.abs())
}

/// Refines one isolated root in a closed interval by Sturm bisection.
///
/// Returns NaN for invalid endpoints, tolerance, polynomial or root counts.
/// If the width cannot be reduced further in floating point, returns an
/// endpoint of the last bracket.
#[must_use]
pub fn refine_root_bisection(
    poly: &Polynomial,
    interval: (f64, f64),
    tolerance: f64,
) -> f64 {
    if !valid_interval(interval) || !valid_precision(tolerance) {
        return f64::NAN;
    }
    if super::exact_eq(interval.0, interval.1) {
        return if validate_coefficients(poly).is_ok() && exact_endpoint_root(poly, interval.0) {
            interval.0
        } else {
            f64::NAN
        };
    }
    let Ok(counter) = SturmCounter::new(poly) else {
        return f64::NAN;
    };
    let bracket = counter.bracket(interval);
    if bracket.left.is_root {
        return bracket.left.x;
    }
    counter.refine(bracket, tolerance).unwrap_or(f64::NAN)
}

/// Finds numerical real roots of the polynomial in ascending order.
///
/// Root isolation and refinement reuse the same Sturm chain. When distinct
/// roots stay together in a bracket narrower than `tolerance`, the method
/// returns one representative for that cluster.
///
/// # Errors
/// Invalid coefficients or tolerance, an unrepresentable root bound or
/// Sturm chain, or failure to isolate or refine a bracket.
pub fn find_roots(
    poly: &Polynomial,
    tolerance: f64,
) -> Result<Vec<f64>, String> {
    if !valid_precision(tolerance) {
        return Err("root tolerance must be finite and positive".into());
    }
    let interval = whole_interval(poly)?;
    let counter = SturmCounter::new(poly)?;
    counter.find(counter.bracket(interval), tolerance)
}

/// Finds numerical real roots in the finite closed interval `[left, right]`.
///
/// Includes roots at either endpoint and returns results in ascending order.
/// Endpoint checks use the original binary coefficients. Only the requested
/// interval is searched, with the tolerance and cluster behavior of [`find_roots`].
///
/// # Errors
/// Non-finite or reversed endpoints, invalid coefficients or tolerance,
/// or a Sturm chain that cannot resolve the requested interval.
pub fn find_roots_in(
    poly: &Polynomial,
    interval: (f64, f64),
    tolerance: f64,
) -> Result<Vec<f64>, String> {
    if !valid_interval(interval) {
        return Err("root interval endpoints must be finite and nondecreasing".into());
    }
    if !valid_precision(tolerance) {
        return Err("root tolerance must be finite and positive".into());
    }
    if super::exact_eq(interval.0, interval.1) {
        validate_coefficients(poly)?;
        return Ok(exact_endpoint_root(poly, interval.0)
            .then_some(interval.0)
            .into_iter()
            .collect());
    }
    let counter = SturmCounter::new(poly)?;
    counter.find(counter.bracket(interval), tolerance)
}
