//! Complex root localization, bounds, polynomial solvers, argument principle,
//! Delves-Lyness, Rouché theorem checks, and 1D complex homotopy continuation.

use crate::kernels::complex_ball::ComplexBall;
use crate::kernels::complex_ball::ComplexBallJet;
use crate::kernels::complex_ball::certify_root_1d;
use num_complex::Complex;
use std::f64::consts::PI;

use super::reuse::Scratch;
use super::reuse::global_pool;

fn eval_poly_and_deriv<const DERIVATIVES: u8>(
    coeffs: &[Complex<f64>],
    z: Complex<f64>,
) -> (Complex<f64>, Complex<f64>, Complex<f64>) {
    let mut p = Complex::new(0.0, 0.0);
    let mut dp = Complex::new(0.0, 0.0);
    let mut ddp = Complex::new(0.0, 0.0);
    for &c in coeffs.iter().rev() {
        if DERIVATIVES > 1 {
            ddp = ddp * z + dp * 2.0;
        }
        if DERIVATIVES > 0 {
            dp = dp * z + p;
        }
        p = p * z + c;
    }
    (p, dp, ddp)
}

fn trim_leading_zeros(coeffs: &[Complex<f64>]) -> &[Complex<f64>] {
    coeffs
        .iter()
        .rposition(|c| c.re != 0.0 || c.im != 0.0)
        .map_or(&[], |degree| &coeffs[..=degree])
}

fn finite_polynomial(coeffs: &[Complex<f64>]) -> Option<(Complex<f64>, &[Complex<f64>])> {
    let (&leading, lower) = trim_leading_zeros(coeffs).split_last()?;
    if lower.is_empty()
        || lower
            .iter()
            .chain(std::iter::once(&leading))
            .any(|c| !c.re.is_finite() || !c.im.is_finite())
    {
        return None;
    }
    Some((leading, lower))
}

/// Computes a log-magnitude even when the complex norm would overflow or
/// round to a subnormal with fewer significant bits.
fn log2_magnitude(coefficient: Complex<f64>) -> f64 {
    let scale = coefficient.re.abs().max(coefficient.im.abs());
    if scale == 0.0 {
        return f64::NEG_INFINITY;
    }
    scale.log2()
        + (coefficient.re / scale)
            .hypot(coefficient.im / scale)
            .log2()
}

/// Classic Cauchy upper bound on the modulus of all roots of
/// $p(z) = \sum_{k=0}^n `c_k` z^k$ (ascending coefficients):
/// $R = 1 + \max_{0 \le k < n} |`c_k` / `c_n`|$.
#[must_use]
pub fn cauchy_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let (leading, lower) = finite_polynomial(coeffs)?;
    let leading_norm = leading.norm();
    let largest = lower.iter().map(|c| c.norm()).fold(0.0, f64::max);
    if largest == 0.0 {
        return Some(1.0);
    }
    let ratio = if leading_norm.is_normal() && largest.is_normal() {
        largest / leading_norm
    } else {
        let max_log = lower
            .iter()
            .copied()
            .map(log2_magnitude)
            .fold(f64::NEG_INFINITY, f64::max);
        (max_log - log2_magnitude(leading)).exp2()
    };
    Some(1.0 + ratio)
}

/// Computes the sharp Cauchy radius, the positive solution of
/// $|`c_n`| R^n = \sum_{k=0}^{n-1} |`c_k`| R^k$ for a non-monomial polynomial.
/// Coefficients are ascending. A monomial has radius zero.
#[must_use]
pub fn cauchy_sharp_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let (leading, lower) = finite_polynomial(coeffs)?;
    let degree = lower.len();
    let mut magnitudes = super::reuse::f64_scratch(degree);
    magnitudes.clear();
    magnitudes.extend(lower.iter().map(|c| c.norm()));
    let lower_max = super::simd::linf_norm_f64(&magnitudes);
    if lower_max == 0.0 {
        return Some(0.0);
    }
    let leading_norm = leading.norm();
    let scale = lower_max.max(leading_norm);
    let normalized_leading = leading_norm / scale;
    // In the ordinary case, compute norms once and evaluate only a
    // reciprocal Horner recurrence during the search (no powers of R).
    if leading_norm.is_normal()
        && normalized_leading.is_normal()
        && magnitudes.iter_mut().all(|m| {
            let normalized = *m / scale;
            let representable = *m == 0.0 || (m.is_normal() && normalized.is_normal());
            *m = normalized;
            representable
        })
    {
        return Some(bracket_positive_radius(|radius| {
            let reciprocal = radius.recip();
            let sum = if reciprocal.is_normal() {
                magnitudes.iter().fold(0.0, |acc, m| (acc + m) * reciprocal)
            } else {
                magnitudes.iter().fold(0.0, |acc, m| (acc + m) / radius)
            };
            sum <= normalized_leading
        }));
    }
    // A common linear scale would lose nonzero coefficients. Keep
    // log-magnitude ratios and sum only terms small enough to fit.
    let leading_log = log2_magnitude(leading);
    magnitudes.clear();
    magnitudes.extend(lower.iter().map(|&c| log2_magnitude(c) - leading_log));
    Some(bracket_positive_radius(|radius| {
        let log_radius = radius.log2();
        let mut sum = 0.0;
        for (power, &log_ratio) in (1..=degree).zip(magnitudes.iter().rev()) {
            let exponent = (-(power as f64)).mul_add(log_radius, log_ratio);
            if exponent > 0.0 {
                return false;
            }
            sum += exponent.exp2();
            if sum > 1.0 {
                return false;
            }
        }
        true
    }))
}

fn bracket_positive_radius(contains: impl Fn(f64) -> bool) -> f64 {
    if !contains(f64::MAX) {
        return f64::INFINITY;
    }
    // Bit patterns order the nonnegative binary64 values. Searching them
    // brackets the root to adjacent floats even across extreme exponents.
    let mut lower = 0;
    let mut upper = f64::MAX.to_bits();
    while upper - lower > 1 {
        let mid = u64::midpoint(lower, upper);
        if contains(f64::from_bits(mid)) {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    f64::from_bits(upper)
}

/// Fujiwara upper bound on the modulus of all polynomial roots:
/// $R = 2 \max_{1 \le k \le n} \left|\frac{c_{n-k}}{c_n}\right|^{1/k}$ with factor $1/2$ on $k=n$.
#[must_use]
pub fn fujiwara_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let (leading, lower) = finite_polynomial(coeffs)?;
    let degree = lower.len();
    let leading_norm = leading.norm();
    let mut leading_log = None;
    let mut bound = 0.0_f64;
    for (offset, &coefficient) in lower.iter().rev().enumerate() {
        let norm = coefficient.norm();
        if norm == 0.0 {
            continue;
        }
        let k = offset + 1;
        let ratio = norm / leading_norm;
        let normal = leading_norm.is_normal() && norm.is_normal() && ratio.is_normal();
        let discount = if k == degree { 0.5 } else { 1.0 };
        let discounted = ratio * discount;
        let term = if normal && k == 1 {
            if degree == 1 {
                ratio
            } else {
                2.0 * ratio
            }
        } else if normal && discounted.is_normal() {
            2.0 * discounted.powf((k as f64).recip())
        } else {
            let base_log = if normal {
                ratio.log2()
            } else {
                log2_magnitude(coefficient)
                    - *leading_log.get_or_insert_with(|| log2_magnitude(leading))
            };
            let log_discount = if k == degree { 1.0 } else { 0.0 };
            (1.0 + (base_log - log_discount) / k as f64)
                .exp2()
                .max(f64::from_bits(1))
        };
        bound = bound.max(term);
    }
    Some(bound)
}

/// Eneström–Kakeya annulus bounds $(r_{\min}, r_{\max})$ for a polynomial
/// with strictly positive real ascending coefficients $`a_0`, \dots, `a_n` > 0$.
#[must_use]
pub fn enestrom_kakeya_bounds(real_positive_coeffs: &[f64]) -> Option<(f64, f64)> {
    let (&first, remaining) = real_positive_coeffs.split_first()?;
    if remaining.is_empty() || !first.is_finite() || first <= 0.0 {
        return None;
    }
    let mut previous = first;
    let mut lower = f64::INFINITY;
    let mut upper = f64::from_bits(1);
    for &coefficient in remaining {
        if !coefficient.is_finite() || coefficient <= 0.0 {
            return None;
        }
        let ratio = previous / coefficient;
        lower = lower.min(ratio.min(f64::MAX));
        upper = upper.max(ratio);
        previous = coefficient;
    }
    Some((lower, upper))
}

fn initial_circle_guesses(c: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = c.len().saturating_sub(1);
    let r = cauchy_sharp_root_bound(c).unwrap_or(1.0).max(0.5);
    let center = if n >= 1 {
        let c_n = c.get(n).copied().unwrap_or(Complex::new(1.0, 0.0));
        let c_n1 = c.get(n.saturating_sub(1)).copied().unwrap_or_default();
        -divide(c_n1, c_n) / (n as f64)
    } else {
        Complex::new(0.0, 0.0)
    };
    (0..n)
        .map(|k| {
            let angle = 2.0 * PI * (k as f64) / (n as f64) + 0.137;
            center + Complex::from_polar(r * 0.6, angle)
        })
        .collect()
}

/// Scales by the largest coefficient component without changing the roots.
/// Returns an error when the polynomial is zero or scaling loses a coefficient.
pub(crate) fn normalized_coefficients(
    coeffs: &[Complex<f64>]
) -> Result<Scratch<'static, Complex<f64>>, String> {
    let c = trim_leading_zeros(coeffs);
    if c.is_empty() || c.iter().any(|coefficient| !coefficient.is_finite()) {
        return Err("expected a nonzero polynomial with finite coefficients".into());
    }
    let scale = c
        .iter()
        .map(|c| c.re.abs().max(c.im.abs()))
        .fold(0.0_f64, f64::max);
    let mut normalized = global_pool().scratch_complex64(c.len());
    for &coefficient in c {
        let scaled = coefficient / scale;
        if lost_component(coefficient, scaled) {
            return Err("coefficient range exceeds normalized binary64 precision".into());
        }
        normalized.push(scaled);
    }
    Ok(normalized)
}

fn lost_component(
    original: Complex<f64>,
    scaled: Complex<f64>,
) -> bool {
    (original.re != 0.0 && scaled.re == 0.0) || (original.im != 0.0 && scaled.im == 0.0)
}

fn divide(
    numerator: Complex<f64>,
    denominator: Complex<f64>,
) -> Complex<f64> {
    let scale = denominator.re.abs().max(denominator.im.abs());
    (numerator / scale) / (denominator / scale)
}

/// Computes a coefficientwise relative backward error by Horner recurrence.
/// For |z| > 1 it evaluates p(z)/z^n via the reversed polynomial.
pub(crate) fn backward_error(
    coeffs: &[Complex<f64>],
    z: Complex<f64>,
) -> f64 {
    if !z.is_finite() {
        return f64::NAN;
    }
    let zero = Complex::new(0.0, 0.0);
    let Some(start) = coeffs.iter().position(|&coefficient| coefficient != zero) else {
        return f64::NAN;
    };
    if z == zero && start > 0 {
        return 0.0;
    }
    let coeffs = &coeffs[start..];
    let reversed = z.norm() > 1.0;
    let argument = if reversed {
        divide(Complex::new(1.0, 0.0), z)
    } else {
        z
    };
    let decay = argument.norm();
    let mut value = zero;
    let mut magnitude = 0.0_f64;
    let evaluate = |&coefficient: &Complex<f64>| {
        value = value * argument + coefficient;
        magnitude = magnitude * decay + coefficient.norm();
    };
    if reversed {
        coeffs.iter().for_each(evaluate);
    } else {
        coeffs.iter().rev().for_each(evaluate);
    }
    value.norm() / magnitude
}

fn validate_iteration(
    max_iter: usize,
    tol: f64,
) -> Result<(), String> {
    if max_iter == 0 || !tol.is_finite() || tol <= 0.0 {
        return Err("expected a positive iteration count and finite positive tolerance".into());
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum SimultaneousMethod {
    DurandKerner,
    AberthEhrlich,
}

impl SimultaneousMethod {
    fn correction(
        self,
        monic: &[Complex<f64>],
        roots: &[Complex<f64>],
        i: usize,
    ) -> Option<Complex<f64>> {
        let zi = roots[i];
        let (value, derivative, _) = match self {
            | Self::DurandKerner => eval_poly_and_deriv::<0>(monic, zi),
            | Self::AberthEhrlich => eval_poly_and_deriv::<1>(monic, zi),
        };
        if value == Complex::new(0.0, 0.0) {
            return Some(value);
        }
        let neighbors = roots
            .iter()
            .enumerate()
            .filter_map(|(j, &root)| (i != j).then_some(zi - root));
        match self {
            | Self::DurandKerner => {
                let denominator: Complex<f64> = neighbors.product();
                (denominator != Complex::new(0.0, 0.0)).then(|| divide(value, denominator))
            },
            | Self::AberthEhrlich => {
                if derivative == Complex::new(0.0, 0.0) {
                    return None;
                }
                let w = divide(value, derivative);
                let sum = neighbors
                    .map(|diff| divide(Complex::new(1.0, 0.0), diff))
                    .sum::<Complex<f64>>();
                let denominator = Complex::new(1.0, 0.0) - w * sum;
                (denominator != Complex::new(0.0, 0.0)).then(|| divide(w, denominator))
            },
        }
    }
}

fn quadratic_roots(
    coefficients: &[Complex<f64>],
    tol: f64,
) -> Option<[Complex<f64>; 2]> {
    let [constant, linear, quadratic] = *coefficients else {
        return None;
    };
    let discriminant = (linear * linear - (4.0 * quadratic) * constant).sqrt();
    let plus = linear + discriminant;
    let minus = linear - discriminant;
    let q = -0.5
        * if plus.norm() >= minus.norm() {
            plus
        } else {
            minus
        };
    let roots = [divide(q, quadratic), divide(constant, q)];
    roots
        .iter()
        .all(|&root| backward_error(coefficients, root) <= tol)
        .then_some(roots)
}

fn simultaneous_roots(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
    method: SimultaneousMethod,
) -> Result<Vec<Complex<f64>>, String> {
    validate_iteration(max_iter, tol)?;
    let coefficients = normalized_coefficients(coeffs)?;
    let n = coefficients.len() - 1;
    let zero = Complex::new(0.0, 0.0);
    let zero_count = coefficients.iter().take_while(|&&c| c == zero).count();
    let (leading, lower) = coefficients[zero_count..]
        .split_last()
        .ok_or("missing leading coefficient")?;
    if lower.is_empty() {
        return Ok(vec![zero; zero_count]);
    }
    if lower.len() == 1 {
        let root = -divide(lower[0], *leading);
        if !root.is_finite() {
            return Err("linear root is outside binary64 range".into());
        }
        let mut roots = vec![zero; n];
        roots[zero_count] = root;
        return Ok(roots);
    }
    if let Some(quadratic) = quadratic_roots(&coefficients[zero_count..], tol) {
        let mut roots = vec![zero; zero_count];
        roots.extend(quadratic);
        return Ok(roots);
    }
    let mut monic = global_pool().scratch_complex64(lower.len() + 1);
    monic.extend(lower.iter().map(|&c| divide(c, *leading)));
    monic.push(Complex::new(1.0, 0.0));
    if monic.iter().any(|c| !c.is_finite())
        || lower
            .iter()
            .zip(monic.iter())
            .any(|(&c, &m)| c != zero && m == zero)
    {
        return Err("monic polynomial is outside binary64 range".into());
    }
    let mut roots = initial_circle_guesses(&monic);
    for sweep in 0..max_iter {
        if sweep.is_multiple_of(16) && crate::cancellation::is_cancelled() {
            return Err("polynomial root iteration cancelled".into());
        }
        let mut max_step = 0.0_f64;
        let mut stalled = false;
        for i in 0..roots.len() {
            let Some(delta) = method.correction(&monic, &roots, i) else {
                stalled = true;
                continue;
            };
            let next = roots[i] - delta;
            let step = delta.norm();
            if !next.is_finite() || !step.is_finite() {
                return Err("polynomial root iteration left binary64 range".into());
            }
            max_step = max_step.max(step);
            roots[i] = next;
        }
        if max_step <= tol
            && roots
                .iter()
                .all(|&root| backward_error(&monic, root) <= tol)
        {
            roots.resize(n, zero);
            return Ok(roots);
        }
        if stalled && max_step == 0.0 {
            return Err("polynomial root iteration stalled above the residual tolerance".into());
        }
    }
    Err(format!(
        "polynomial roots did not converge within {max_iter} iterations"
    ))
}

/// Finds approximate roots with multiplicities from ascending coefficients using Durand–Kerner.
///
/// # Errors
///
/// Returns an error for an empty or zero polynomial, invalid controls, nonfinite iterates,
/// cancellation, or failure to converge. A nonzero constant has no roots.
pub fn durand_kerner(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    simultaneous_roots(coeffs, max_iter, tol, SimultaneousMethod::DurandKerner)
}

/// Finds approximate roots with multiplicities from ascending coefficients using Aberth–Ehrlich.
///
/// # Errors
///
/// Returns an error for an empty or zero polynomial, invalid controls, nonfinite iterates,
/// cancellation, or failure to converge. A nonzero constant has no roots.
pub fn aberth_ehrlich(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    simultaneous_roots(coeffs, max_iter, tol, SimultaneousMethod::AberthEhrlich)
}

/// Finds an approximate root of an ascending polynomial starting from `z0` via Laguerre's method.
/// `tol` controls the relative backward error of the returned candidate.
///
/// # Errors
///
/// Returns an error for a constant or zero polynomial, invalid controls, nonfinite iterates,
/// cancellation, a vanishing Laguerre denominator, or failure to converge.
pub fn laguerre_root(
    coeffs: &[Complex<f64>],
    z0: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String> {
    validate_iteration(max_iter, tol)?;
    let c = normalized_coefficients(coeffs)?;
    let n = c.len() - 1;
    if n == 0 || !z0.is_finite() {
        return Err("Laguerre requires a nonconstant polynomial and a finite seed".into());
    }
    if n == 1 {
        let root = -divide(c[0], c[1]);
        if root.is_finite() && backward_error(&c, root) <= tol {
            return Ok(root);
        }
        return Err("linear root exceeds binary64 range or requested tolerance".into());
    }
    let n_f = n as f64;
    let mut z = z0;
    for step in 0..max_iter {
        if step.is_multiple_of(16) && crate::cancellation::is_cancelled() {
            return Err("Laguerre iteration cancelled".into());
        }
        if backward_error(&c, z) <= tol {
            return Ok(z);
        }
        let (p, dp, ddp) = eval_poly_and_deriv::<2>(&c, z);
        let g = divide(dp, p);
        let g_squared = g * g;
        let h = g_squared - divide(ddp, p);
        let disc = ((n_f - 1.0) * (n_f * h - g_squared)).sqrt();
        let plus = g + disc;
        let minus = g - disc;
        let denominator = if plus.norm() >= minus.norm() {
            plus
        } else {
            minus
        };
        if denominator == Complex::new(0.0, 0.0) || !denominator.is_finite() {
            return Err("Laguerre denominator vanished or exceeded binary64 range".into());
        }
        z -= divide(Complex::new(n_f, 0.0), denominator);
    }
    if backward_error(&c, z) <= tol {
        return Ok(z);
    }
    Err(format!(
        "Laguerre did not converge within {max_iter} iterations"
    ))
}

fn checked_sample(z: Complex<f64>) -> Result<Complex<f64>, String> {
    z.is_finite()
        .then_some(z)
        .ok_or_else(|| "root callback returned a nonfinite value".into())
}

fn root_iteration_options(
    max_iter: usize,
    tol: f64,
) -> Result<(), String> {
    if max_iter == 0 || !tol.is_finite() || tol <= 0.0 {
        Err("expected a positive iteration bound and finite positive residual tolerance".into())
    } else {
        Ok(())
    }
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "common scaling uses finite floats and a checked nonzero scale"
)]
fn rescale_root_samples(values: [Complex<f64>; 3]) -> [Complex<f64>; 3] {
    let scale = values
        .iter()
        .fold(0.0_f64, |bound, z| bound.max(z.re.abs()).max(z.im.abs()));
    if scale > 0.0 {
        values.map(|z| z / scale)
    } else {
        values
    }
}

/// Finds a complex root of an analytic function using Muller's quadratic interpolation.
///
/// # Errors
///
/// Rejects invalid iteration parameters or seeds, singular interpolation,
/// nonfinite callback values, cancellation or failure to meet the residual tolerance.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "Muller interpolation operates on finite complex samples and checks singular updates"
)]
pub fn muller_root<F>(
    f: F,
    mut z0: Complex<f64>,
    mut z1: Complex<f64>,
    mut z2: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    root_iteration_options(max_iter, tol)?;
    if [z0, z1, z2].iter().any(|point| !point.is_finite()) {
        return Err("Muller seeds must be finite".into());
    }
    if crate::cancellation::is_cancelled() {
        return Err("Muller solve cancelled".into());
    }
    let mut f0 = checked_sample(f(z0))?;
    let mut f1 = checked_sample(f(z1))?;
    let mut f2 = checked_sample(f(z2))?;
    if f2.norm() <= tol {
        return Ok(z2);
    }
    let zero = Complex::new(0.0, 0.0);
    for _ in 0..max_iter {
        if crate::cancellation::is_cancelled() {
            return Err("Muller solve cancelled".into());
        }
        let h1 = z1 - z0;
        let h2 = z2 - z1;
        let span = z2 - z0;
        if h1 == zero || h2 == zero || span == zero {
            return Err("Muller interpolation needs distinct points".into());
        }
        // A common value scale preserves the interpolation step while keeping
        // the discriminant away from avoidable overflow and underflow.
        let [a, b, c] = rescale_root_samples([f0, f1, f2]);
        let d1 = divide(b - a, h1);
        let d2 = divide(c - b, h2);
        let quadratic = divide(d2 - d1, span);
        let linear = quadratic * h2 + d2;
        let rad = (linear * linear - quadratic * c * 4.0).sqrt();
        let den_plus = linear + rad;
        let den_minus = linear - rad;
        let den = if den_plus.norm() >= den_minus.norm() {
            den_plus
        } else {
            den_minus
        };
        if den == zero || !den.is_finite() {
            return Err("Muller denominator is singular".into());
        }
        let next = z2 - divide(2.0 * c, den);
        if !next.is_finite() || next == z2 {
            return Err("Muller update stagnated above the residual tolerance".into());
        }
        z0 = z1;
        z1 = z2;
        z2 = next;
        f0 = f1;
        f1 = f2;
        // The two previous samples are already available after shifting.
        f2 = checked_sample(f(next))?;
        if f2.norm() <= tol {
            return Ok(next);
        }
    }
    Err(format!(
        "Muller failed to converge within {max_iter} iterations"
    ))
}

/// Finds a complex root using Halley's cubically convergent method.
///
/// # Errors
///
/// Rejects invalid iteration parameters, singular updates, nonfinite callback
/// values, cancellation or failure to meet the residual tolerance.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "Halley iteration operates on finite complex samples and checks singular updates"
)]
pub fn halley_root<F, DF, DDF>(
    f: F,
    df: DF,
    ddf: DDF,
    mut z: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
    DDF: Fn(Complex<f64>) -> Complex<f64>,
{
    root_iteration_options(max_iter, tol)?;
    if !z.is_finite() {
        return Err("Halley seed must be finite".into());
    }
    if crate::cancellation::is_cancelled() {
        return Err("Halley solve cancelled".into());
    }
    let zero = Complex::new(0.0, 0.0);
    let mut fz = checked_sample(f(z))?;
    for _ in 0..max_iter {
        if fz.norm() <= tol {
            return Ok(z);
        }
        if crate::cancellation::is_cancelled() {
            return Err("Halley solve cancelled".into());
        }
        let [f0, f1, f2] =
            rescale_root_samples([fz, checked_sample(df(z))?, checked_sample(ddf(z))?]);
        let den = f1 * f1 * 2.0 - f0 * f2;
        if den == zero || !den.is_finite() {
            return Err("Halley denominator is singular".into());
        }
        let next = z - divide(f0 * f1 * 2.0, den);
        if !next.is_finite() || next == z {
            return Err("Halley update stagnated above the residual tolerance".into());
        }
        z = next;
        fz = checked_sample(f(z))?;
    }
    if fz.norm() <= tol {
        Ok(z)
    } else {
        Err(format!(
            "Halley failed to converge within {max_iter} iterations"
        ))
    }
}

pub(crate) fn validate_circle(
    center: Complex<f64>,
    radius: f64,
    samples: usize,
) -> Result<(), String> {
    if !center.is_finite() || !radius.is_finite() || radius <= 0.0 || samples < 8 {
        return Err(
            "expected a finite circle of positive radius and at least eight samples".into(),
        );
    }
    if [center.re, center.im].into_iter().any(|coordinate| {
        let lower = coordinate - radius;
        let upper = coordinate + radius;
        !(lower.is_finite() && upper.is_finite() && lower < coordinate && coordinate < upper)
    }) {
        return Err("circle radius does not resolve around this center in binary64".into());
    }
    Ok(())
}

/// Computes power-sum moments $`s_k` = \frac{1}{2\pi i} \oint_{|z-c|=r} z^k \frac{f'(z)}{f(z)} dz$
/// for $k = 0, \dots, \text{max\_moment}$ using the periodic trapezoidal rule.
///
/// # Errors
///
/// Returns an error for invalid circle geometry, insufficient samples, cancellation,
/// a near-zero function value on the circle, or a nonfinite callback or moment.
pub fn argument_principle_moments<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    max_moment: usize,
    n_samples: usize,
) -> Result<Vec<Complex<f64>>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    validate_circle(center, radius, n_samples)?;
    let Some(count) = max_moment.checked_add(1) else {
        return Err("moment order exceeds addressable capacity".into());
    };
    if count.saturating_mul(n_samples) > 33_554_432 {
        return Err("moment grid exceeds the work limit".into());
    }
    let mut moments = vec![Complex::new(0.0, 0.0); count];
    visit_circle_weights(f, df, center, radius, n_samples, |_, z, weight| {
        add_moment_sample(&mut moments, z, weight);
    })?;
    if moments.iter().any(|moment| !moment.is_finite()) {
        return Err("circle moments are outside binary64 range".into());
    }
    Ok(moments)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the moment recurrence uses finite IEEE complex values and checks the resulting moments"
)]
fn add_moment_sample(
    moments: &mut [Complex<f64>],
    z: Complex<f64>,
    weight: Complex<f64>,
) {
    let mut z_pow = Complex::new(1.0, 0.0);
    for slot in moments {
        *slot += weight * z_pow;
        z_pow *= z;
    }
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "fixed quadrature indices fit in usize and finite complex results are checked by callers"
)]
fn visit_circle_weights(
    f: impl Fn(Complex<f64>) -> Complex<f64>,
    df: impl Fn(Complex<f64>) -> Complex<f64>,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
    mut sample: impl FnMut(Complex<f64>, Complex<f64>, Complex<f64>),
) -> Result<(), String> {
    let n_f = n_samples as f64;
    for j in 0..n_samples {
        if j.is_multiple_of(64) && crate::cancellation::is_cancelled() {
            return Err("circle sampling cancelled".into());
        }
        let theta = 2.0 * PI * (j as f64) / n_f;
        let offset = Complex::from_polar(radius, theta);
        let z = center + offset;
        if !z.is_finite() {
            return Err("circle sample is outside binary64 range".into());
        }
        let fz = f(z);
        if !fz.is_finite() {
            return Err("circle callback returned a nonfinite function value".into());
        }
        if fz.norm() < 1e-14 {
            return Err("Function has a zero/near-zero on the integration circle".to_string());
        }
        let log_deriv = divide(df(z), fz);
        if !log_deriv.is_finite() {
            return Err("circle callback returned a nonfinite logarithmic derivative".into());
        }
        // dz = i * offset * dtheta, divided by (2 * pi * i) -> offset / n_samples
        let weight = log_deriv * offset / n_f;
        sample(offset, z, weight);
    }
    Ok(())
}

/// Computes the argument principle winding number `(Z - P) = \text{round}(\text{Re}(s_0))`.
///
/// # Errors
///
/// Returns an error when moment evaluation fails or the winding estimate is too far from an integer.
pub fn argument_principle_winding<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<isize, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    let m = argument_principle_moments(f, df, center, radius, 0, n_samples)?;
    let s0 = m
        .first()
        .copied()
        .ok_or_else(|| "Missing moment".to_string())?;
    let integral = s0.re.round();
    if !s0.is_finite() || (s0 - integral).norm() > 1e-2 {
        return Err(format!(
            "argument principle estimate is not near an integer: {s0}"
        ));
    }
    if integral < isize::MIN as f64 || integral >= -(isize::MIN as f64) {
        return Err("winding estimate exceeds signed integer range".into());
    }
    Ok(integral as isize)
}

/// Finds zeros of an analytic function in a disk using Delves–Lyness localization.
///
/// Uses argument-principle moments, Newton’s identities to form the monic polynomial, and
/// Aberth–Ehrlich iteration to extract its zeros. The function must be pole-free in the disk.
///
/// # Errors
///
/// Returns an error when moment evaluation or the polynomial solve fails, or the winding number is invalid.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "moment indexing uses checked and saturating sizes and complex arithmetic propagates numeric failures"
)]
pub fn delves_lyness_roots<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<Vec<Complex<f64>>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    // Keep one sample of each callback per point, and form moments only
    // after the sampled winding establishes the degree of the polynomial.
    validate_circle(center, radius, n_samples)?;
    if n_samples > 1_048_576 {
        return Err("circle grid exceeds the sampling budget".into());
    }
    let pool = global_pool();
    let mut coordinates = pool.scratch_complex64(n_samples);
    let mut weights = pool.scratch_complex64(n_samples);
    let mut s0 = Complex::new(0.0, 0.0);
    visit_circle_weights(&f, &df, center, radius, n_samples, |offset, _, weight| {
        coordinates.push(offset / radius);
        weights.push(weight);
        s0 += weight;
    })?;
    let rounded = s0.re.round();
    if !s0.is_finite() || (s0 - rounded).norm() > 1e-2 {
        return Err(format!("Non-integer argument principle integral s0 = {s0}"));
    }
    if rounded < 0.0 {
        return Err("Negative winding number indicates poles inside contour".into());
    }
    let n = rounded as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    if n > 4096 || n.saturating_mul(n_samples) > 33_554_432 {
        return Err("root moment grid exceeds the work limit".into());
    }
    let mut shifted_moments = vec![Complex::new(0.0, 0.0); n.saturating_add(1)];
    for (&point, &weight) in coordinates.iter().zip(weights.iter()) {
        add_moment_sample(&mut shifted_moments, point, weight);
    }
    // Newton identities: k * e_k = sum_{j=1}^k (-1)^{j-1} * e_{k-j} * s_j
    let mut e = vec![Complex::new(0.0, 0.0); n.saturating_add(1)];
    if let Some(e0) = e.first_mut() {
        *e0 = Complex::new(1.0, 0.0);
    }
    for k in 1..=n {
        let mut acc = Complex::new(0.0, 0.0);
        for j in 1..=k {
            let sign = if j % 2 == 1 { 1.0 } else { -1.0 };
            let e_prev = e.get(k.saturating_sub(j)).copied().unwrap_or_default();
            let s_j = shifted_moments.get(j).copied().unwrap_or_default();
            acc += e_prev * s_j * sign;
        }
        if let Some(ek) = e.get_mut(k) {
            *ek = acc / (k as f64);
        }
    }
    // Monic polynomial in ascending powers: sum_{m=0}^n (-1)^{n-m} e_{n-m} w^m
    let mut poly = Vec::with_capacity(n.saturating_add(1));
    for m in 0..=n {
        let idx = n.saturating_sub(m);
        let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
        let ek = e.get(idx).copied().unwrap_or_default();
        poly.push(ek * sign);
    }
    let shifted_roots = aberth_ehrlich(&poly, 200, 1e-13)?;
    // Polish each root with Newton steps on f
    let mut polished: Vec<Complex<f64>> = shifted_roots
        .into_iter()
        .map(|w| center + radius * w)
        .collect();
    for z in &mut polished {
        for _ in 0..8 {
            if crate::cancellation::is_cancelled() {
                return Err("Delves–Lyness refinement cancelled".into());
            }
            let fz = checked_sample(f(*z))?;
            if fz.norm() < 1e-14 {
                break;
            }
            let dfz = checked_sample(df(*z))?;
            if dfz == Complex::new(0.0, 0.0) {
                break;
            }
            let next = *z - divide(fz, dfz);
            if !next.is_finite() || next == *z || (next - center).norm() >= radius {
                break;
            }
            *z = next;
        }
    }
    Ok(polished)
}

/// Samples Rouché's condition $|f(z)| > |g(z)|$ on $|z - \text{center}| = \text{radius}$.
///
/// # Errors
///
/// The result describes the sampled points.
/// Returns an error for invalid geometry, insufficient samples, cancellation, or a nonfinite sample.
pub fn rouche_check_circle<F, G>(
    f: F,
    g: G,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<bool, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    G: Fn(Complex<f64>) -> Complex<f64>,
{
    validate_circle(center, radius, n_samples)?;
    let mut dominates = true;
    for j in 0..n_samples {
        if j.is_multiple_of(64) && crate::cancellation::is_cancelled() {
            return Err("circle sampling cancelled".into());
        }
        let theta = 2.0 * PI * (j as f64) / (n_samples as f64);
        let z = center + Complex::from_polar(radius, theta);
        if !z.is_finite() {
            return Err("circle sample is outside binary64 range".into());
        }
        let left = f(z);
        let right = g(z);
        if !left.is_finite() || !right.is_finite() {
            return Err("Rouché callback returned a nonfinite sample".into());
        }
        let scale = left
            .re
            .abs()
            .max(left.im.abs())
            .max(right.re.abs())
            .max(right.im.abs());
        dominates &= scale > 0.0 && (left / scale).norm() > (right / scale).norm();
    }
    Ok(dominates)
}

/// Tracks a complex homotopy root from `t = 0` to `t = 1`.
///
/// Starts at `z(0) = z_0` for `H(z(t), t) = 0`, using an adaptive RK4 predictor and Newton
/// corrector with damping near a singular target.
///
/// # Errors
///
/// Returns an error for a zero step count, a singular Jacobian, or failure to reach and correct the endpoint within the step budget.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "the bounded predictor and corrector use checked floating-point values"
)]
pub fn homotopy_track_1d<H, DHZ, DHT>(
    h: H,
    dh_dz: DHZ,
    dh_dt: DHT,
    z0: Complex<f64>,
    n_steps: usize,
    newton_iters: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    H: Fn(Complex<f64>, f64) -> Complex<f64>,
    DHZ: Fn(Complex<f64>, f64) -> Complex<f64>,
    DHT: Fn(Complex<f64>, f64) -> Complex<f64>,
{
    root_iteration_options(newton_iters, tol)?;
    if n_steps == 0 || !z0.is_finite() {
        return Err("homotopy needs a finite seed and a positive step count".into());
    }
    if crate::cancellation::is_cancelled() {
        return Err("homotopy cancelled".into());
    }
    if checked_sample(h(z0, 0.0))?.norm() > tol {
        return Err("initial homotopy point fails the residual tolerance".into());
    }
    let mut dt = 1.0 / (n_steps as f64);
    let min_dt = 1e-6;
    let mut t = 0.0_f64;
    let mut z = z0;
    for _ in 0..50_000 {
        if crate::cancellation::is_cancelled() {
            return Err("homotopy cancelled".into());
        }
        if t >= 1.0 {
            return Ok(z);
        }
        let step = dt.min(1.0 - t);
        let t_next = (t + step).min(1.0);

        // RK4 predictor for dz/dt = -(dH/dt)/(dH/dz)
        let vel = |zz: Complex<f64>, tt: f64| -> Option<Complex<f64>> {
            let hz = dh_dz(zz, tt);
            if !hz.is_finite() || hz == Complex::new(0.0, 0.0) {
                None
            } else {
                let direction = divide(-dh_dt(zz, tt), hz);
                direction.is_finite().then_some(direction)
            }
        };
        let z_pred = if let Some(k1) = vel(z, t)
            && let Some(k2) = vel(z + 0.5 * step * k1, t + 0.5 * step)
            && let Some(k3) = vel(z + 0.5 * step * k2, t + 0.5 * step)
            && let Some(k4) = vel(z + step * k3, t_next)
        {
            z + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        } else if t_next >= 0.99 {
            z
        } else if step > min_dt {
            dt *= 0.5;
            continue;
        } else {
            return Err(format!("Singular Jacobian along homotopy at t = {t}"));
        };

        // Newton corrector at t_next
        let mut z_corr = z_pred;
        let mut converged = false;
        for _ in 0..newton_iters {
            let val = checked_sample(h(z_corr, t_next))?;
            if val.norm() < tol {
                converged = true;
                break;
            }
            let dval = checked_sample(dh_dz(z_corr, t_next))?;
            if dval == Complex::new(0.0, 0.0) {
                break;
            }
            z_corr -= divide(val, dval);
        }
        if !converged {
            converged = checked_sample(h(z_corr, t_next))?.norm() <= tol;
            if !converged {
                if step > min_dt {
                    dt *= 0.5;
                    continue;
                }
                return Err(format!("homotopy correction failed at t = {t_next}"));
            }
        }
        z = z_corr;
        t = t_next;
        if converged && dt < 0.1 {
            dt = (dt * 1.25).min(0.1);
        }
    }
    if t >= 1.0 {
        Ok(z)
    } else {
        Err("Homotopy step budget exceeded".into())
    }
}

// Scaling only by a power of two preserves the encoded polynomial. A
// round-trip check detects lost subnormal components before certification.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "finitely bounded scaling uses a nonzero power of two and verifies round-trip exactness"
)]
fn exact_enclosure_coefficients(coeffs: &[Complex<f64>]) -> Option<Scratch<'static, Complex<f64>>> {
    let max_component = coeffs.iter().fold(0.0_f64, |scale, value| {
        scale.max(value.re.abs()).max(value.im.abs())
    });
    let scale_bits = max_component.to_bits() & 0x7ff0_0000_0000_0000;
    let scale = f64::from_bits(scale_bits).max(f64::MIN_POSITIVE);
    let mut verification = global_pool().scratch_complex64(coeffs.len());
    for &coefficient in coeffs {
        let value = coefficient / scale;
        let recovered = value * scale;
        if recovered.re.to_bits() != coefficient.re.to_bits()
            || recovered.im.to_bits() != coefficient.im.to_bits()
        {
            return None;
        }
        verification.push(value);
    }
    Some(verification)
}

fn enclosure_jet(
    point: ComplexBallJet,
    coeffs: &[Complex<f64>],
) -> ComplexBallJet {
    coeffs.iter().rev().fold(
        ComplexBallJet::constant(ComplexBall::zero()),
        |value, &coefficient| {
            value
                .mul(&point)
                .add(&ComplexBallJet::constant(ComplexBall::of(coefficient)))
        },
    )
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "binary64 guesses only propose disks; Krawczyk verifies every accepted certificate"
)]
fn seeded_certificates(
    coeffs: &[Complex<f64>],
    tol: f64,
) -> Vec<ComplexBall> {
    let Ok(guesses) = aberth_ehrlich(coeffs, 256, 1e-13) else {
        return Vec::new();
    };
    let degree = guesses.len();
    let mut certified = Vec::with_capacity(degree);
    for guess in guesses {
        if crate::cancellation::is_cancelled() {
            break;
        }
        let mut radius = tol.max(32.0 * f64::EPSILON * guess.norm().max(1.0));
        for _ in 0..8 {
            let Some(disk) = ComplexBall::new(guess, radius) else {
                break;
            };
            if let Some(proven) = certify_root_1d(|jet| Some(enclosure_jet(jet, coeffs)), disk) {
                if !certified
                    .iter()
                    .any(|existing: &ComplexBall| existing.overlaps(&proven))
                {
                    certified.push(proven);
                }
                break;
            }
            radius *= 2.0;
        }
    }
    certified
}

/// Certifies disjoint disks, each isolating a simple root of the supplied complex polynomial.
///
/// Coefficients are in ascending order. The result can be partial if some roots
/// are multiple or too near each other, the precision or search limit is reached,
/// or evaluation is cancelled.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "numeric candidate disks are accepted only after verified interval checks"
)]
pub fn isolate_roots_verified(
    coeffs: &[Complex<f64>],
    tol: f64,
) -> Vec<ComplexBall> {
    let coeffs = trim_leading_zeros(coeffs);
    if coeffs.len() < 2
        || coeffs.iter().any(|coefficient| !coefficient.is_finite())
        || !tol.is_finite()
        || tol <= 0.0
        || crate::cancellation::is_cancelled()
    {
        return Vec::new();
    }
    let degree = coeffs.len() - 1;
    let exact = exact_enclosure_coefficients(coeffs);
    let coeffs = exact.as_deref().map_or(coeffs, Vec::as_slice);
    let mut certified = if degree <= 256 {
        seeded_certificates(coeffs, tol)
    } else {
        Vec::new()
    };
    if certified.len() == degree {
        return certified;
    }
    let Some(bound) = cauchy_root_bound(coeffs) else {
        return certified;
    };
    let radius = bound * 1.05 * 2.0_f64.sqrt();
    let mut stack = vec![(Complex::new(0.0, 0.0), radius)];
    for _ in 0..20_000 {
        if crate::cancellation::is_cancelled() {
            break;
        }
        let Some((center, radius)) = stack.pop() else {
            break;
        };
        let Some(ball) = ComplexBall::new(center, radius) else {
            continue;
        };
        let mut p_ball = ComplexBall::zero();
        for &coefficient in coeffs.iter().rev() {
            p_ball = p_ball * ball + ComplexBall::of(coefficient);
        }
        if !p_ball.contains_zero() {
            continue;
        }
        if let Some(proven) = certify_root_1d(|jet| Some(enclosure_jet(jet, coeffs)), ball) {
            if !certified.iter().any(|existing| existing.overlaps(&proven)) {
                certified.push(proven);
                if certified.len() == degree {
                    break;
                }
            }
            continue;
        }
        if radius <= tol.max(1e-9) {
            continue;
        }
        let half = 0.5 * radius;
        let offset = half / 2.0_f64.sqrt();
        for (sx, sy) in [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0)] {
            stack.push((
                Complex::new(center.re + sx * offset, center.im + sy * offset),
                half * 1.05,
            ));
        }
    }
    certified
}
