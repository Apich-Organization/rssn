//! Complex root localization, bounds, polynomial solvers, argument principle,
//! Delves-Lyness, Rouché theorem checks, and 1D complex homotopy continuation.

use num_complex::Complex;
use std::f64::consts::PI;

fn eval_poly_and_deriv(
    coeffs: &[Complex<f64>],
    z: Complex<f64>,
) -> (Complex<f64>, Complex<f64>, Complex<f64>) {
    let mut p = Complex::new(0.0, 0.0);
    let mut dp = Complex::new(0.0, 0.0);
    let mut ddp = Complex::new(0.0, 0.0);
    for &c in coeffs.iter().rev() {
        ddp = ddp * z + dp * 2.0;
        dp = dp * z + p;
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
        -c_n1 / (c_n * (n as f64))
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

/// Finds all roots of a polynomial with ascending coefficients using the Durand–Kerner method.
///
/// # Errors
///
/// Returns an error for an empty coefficient slice or when the iterations fail to converge.
pub fn durand_kerner(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    let cn = c
        .last()
        .copied()
        .ok_or_else(|| "Empty coefficients".to_string())?;
    let monic: Vec<Complex<f64>> = c.iter().map(|&v| v / cn).collect();
    let mut roots = initial_circle_guesses(&monic);
    for _ in 0..max_iter {
        let mut max_step = 0.0_f64;
        for i in 0..n {
            let zi = *roots.get(i).ok_or_else(|| "Index error".to_string())?;
            let (p_val, _, _) = eval_poly_and_deriv(&monic, zi);
            let mut denom = Complex::new(1.0, 0.0);
            for (j, &zj) in roots.iter().enumerate() {
                if i != j {
                    denom *= zi - zj;
                }
            }
            if denom.norm() < 1e-30 {
                continue;
            }
            let delta = p_val / denom;
            max_step = max_step.max(delta.norm());
            if let Some(slot) = roots.get_mut(i) {
                *slot = zi - delta;
            }
        }
        if max_step < tol {
            return Ok(roots);
        }
    }
    Err(format!(
        "Durand-Kerner failed to converge within {max_iter} iterations"
    ))
}

/// Finds all roots of a polynomial with ascending coefficients using the Aberth–Ehrlich method.
///
/// # Errors
///
/// Returns an error for an empty coefficient slice or when the iterations fail to converge.
pub fn aberth_ehrlich(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    let cn = c
        .last()
        .copied()
        .ok_or_else(|| "Empty coefficients".to_string())?;
    let monic: Vec<Complex<f64>> = c.iter().map(|&v| v / cn).collect();
    let mut roots = initial_circle_guesses(&monic);
    for _ in 0..max_iter {
        let mut max_step = 0.0_f64;
        for i in 0..n {
            let zi = *roots.get(i).ok_or_else(|| "Index error".to_string())?;
            let (p_val, dp_val, _) = eval_poly_and_deriv(&monic, zi);
            if p_val.norm() < tol * 1e-3 {
                continue;
            }
            if dp_val.norm() < 1e-30 {
                if let Some(slot) = roots.get_mut(i) {
                    *slot = zi + Complex::new(1e-5, 1.3e-5);
                }
                continue;
            }
            let w = p_val / dp_val;
            let mut sum_inv = Complex::new(0.0, 0.0);
            for (j, &zj) in roots.iter().enumerate() {
                if i != j {
                    let diff = zi - zj;
                    if diff.norm() > 1e-30 {
                        sum_inv += Complex::new(1.0, 0.0) / diff;
                    }
                }
            }
            let denom = Complex::new(1.0, 0.0) - w * sum_inv;
            let delta = if denom.norm() > 1e-30 {
                w / denom
            } else {
                w
            };
            max_step = max_step.max(delta.norm());
            if let Some(slot) = roots.get_mut(i) {
                *slot = zi - delta;
            }
        }
        if max_step < tol {
            return Ok(roots);
        }
    }
    Err(format!(
        "Durand-Kerner failed to converge within {max_iter} iterations"
    ))
}

/// Finds a single root of a polynomial starting from `z0` via Laguerre's method.
///
/// # Errors
///
/// Returns an error for a constant polynomial, a vanishing Laguerre denominator, or nonconvergence.
pub fn laguerre_root(
    coeffs: &[Complex<f64>],
    z0: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Err("Constant polynomial has no roots".to_string());
    }
    let n_f = n as f64;
    let mut z = z0;
    for _ in 0..max_iter {
        let (p, dp, ddp) = eval_poly_and_deriv(c, z);
        if p.norm() < tol {
            return Ok(z);
        }
        let g = dp / p;
        let g_squared = g * g;
        let h = g_squared - ddp / p;
        let disc = ((n_f - 1.0) * (n_f * h - g_squared)).sqrt();
        let d_plus = g + disc;
        let d_minus = g - disc;
        let denom = if d_plus.norm() >= d_minus.norm() {
            d_plus
        } else {
            d_minus
        };
        if denom.norm() < 1e-30 {
            return Err("Laguerre denominator vanished".to_string());
        }
        let a = Complex::new(n_f, 0.0) / denom;
        z -= a;
        if a.norm() < tol {
            return Ok(z);
        }
    }
    Err(format!(
        "Laguerre failed to converge within {max_iter} iterations"
    ))
}

/// Finds a complex root of an arbitrary analytic function `f` using Muller's quadratic interpolation.
///
/// # Errors
///
/// Returns an error if the interpolation denominator vanishes or the iterations fail to converge.
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
    for _ in 0..max_iter {
        let f0 = f(z0);
        let f1 = f(z1);
        let f2 = f(z2);
        if f2.norm() < tol {
            return Ok(z2);
        }
        let h1 = z1 - z0;
        let h2 = z2 - z1;
        if h1.norm() < 1e-30 || h2.norm() < 1e-30 {
            return Ok(z2);
        }
        let d1 = (f1 - f0) / h1;
        let d2 = (f2 - f1) / h2;
        let a = (d2 - d1) / (h2 + h1);
        let b = a * h2 + d2;
        let c = f2;
        let rad = (b * b - a * c * 4.0).sqrt();
        let den_plus = b + rad;
        let den_minus = b - rad;
        let den = if den_plus.norm() >= den_minus.norm() {
            den_plus
        } else {
            den_minus
        };
        if den.norm() < 1e-30 {
            return Err("Muller denominator vanished".to_string());
        }
        let dz = -c * 2.0 / den;
        let z3 = z2 + dz;
        if dz.norm() < tol {
            return Ok(z3);
        }
        z0 = z1;
        z1 = z2;
        z2 = z3;
    }
    Err(format!(
        "Muller failed to converge within {max_iter} iterations"
    ))
}

/// Finds a complex root using Halley's cubically convergent method.
///
/// # Errors
///
/// Returns an error if the iterations fail to converge.
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
    for _ in 0..max_iter {
        let fz = f(z);
        if fz.norm() < tol {
            return Ok(z);
        }
        let dfz = df(z);
        let ddfz = ddf(z);
        let denom = dfz * dfz * 2.0 - fz * ddfz;
        if denom.norm() < 1e-30 {
            return Err("Halley denominator vanished".to_string());
        }
        let dz = (fz * dfz * 2.0) / denom;
        z -= dz;
        if dz.norm() < tol {
            return Ok(z);
        }
    }
    Err(format!(
        "Halley failed to converge within {max_iter} iterations"
    ))
}

/// Computes power-sum moments $`s_k` = \frac{1}{2\pi i} \oint_{|z-c|=r} z^k \frac{f'(z)}{f(z)} dz$
/// for $k = 0, \dots, \text{max\_moment}$ using the periodic trapezoidal rule.
///
/// # Errors
///
/// Returns an error for an invalid radius or sample count, or if the function is nearly zero on the circle.
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
    if radius <= 0.0 || n_samples < 8 {
        return Err("Invalid circle radius or sample count".to_string());
    }
    let mut moments = vec![Complex::new(0.0, 0.0); max_moment.saturating_add(1)];
    let n_f = n_samples as f64;
    for j in 0..n_samples {
        let theta = 2.0 * PI * (j as f64) / n_f;
        let offset = Complex::from_polar(radius, theta);
        let z = center + offset;
        let fz = f(z);
        if fz.norm() < 1e-14 {
            return Err("Function has a zero/near-zero on the integration circle".to_string());
        }
        let log_deriv = df(z) / fz;
        // dz = i * offset * dtheta, divided by (2 * pi * i) -> offset / n_samples
        let weight = log_deriv * offset / n_f;
        let mut z_pow = Complex::new(1.0, 0.0);
        for slot in &mut moments {
            *slot += weight * z_pow;
            z_pow *= z;
        }
    }
    Ok(moments)
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
    if s0.im.abs() > 1e-2 {
        return Err(format!(
            "Non-integer winding number imaginary part: {}",
            s0.im
        ));
    }
    Ok(s0.re.round() as isize)
}

/// Finds zeros of an analytic function in a disk using Delves–Lyness localization.
///
/// Uses argument-principle moments, Newton’s identities to form the monic polynomial, and
/// Aberth–Ehrlich iteration to extract its zeros. The function must be pole-free in the disk.
///
/// # Errors
///
/// Returns an error when moment evaluation or the polynomial solve fails, or the winding number is invalid.
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
    let s0_vec = argument_principle_moments(&f, &df, center, radius, 0, n_samples)?;
    let s0 = s0_vec
        .first()
        .copied()
        .ok_or_else(|| "Missing s0".to_string())?;
    if s0.im.abs() > 1e-2 || (s0.re - s0.re.round()).abs() > 1e-2 {
        return Err(format!("Non-integer argument principle integral s0 = {s0}"));
    }
    let count_signed = s0.re.round() as isize;
    if count_signed < 0 {
        return Err("Negative winding number indicates poles inside contour".to_string());
    }
    let n = count_signed as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    // Shift origin to center for numerical stability of moments
    let shifted_moments = argument_principle_moments(
        |w| f(center + w),
        |w| df(center + w),
        Complex::new(0.0, 0.0),
        radius,
        n,
        n_samples,
    )?;
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
    let mut polished: Vec<Complex<f64>> = shifted_roots.into_iter().map(|w| center + w).collect();
    for z in &mut polished {
        for _ in 0..8 {
            let fz = f(*z);
            if fz.norm() < 1e-14 {
                break;
            }
            let dfz = df(*z);
            if dfz.norm() < 1e-14 {
                break;
            }
            *z -= fz / dfz;
        }
    }
    Ok(polished)
}

/// Checks Rouché's condition $|f(z)| > |g(z)|$ on $|z - \text{center}| = \text{radius}$.
///
/// # Errors
///
/// Returns an error for an invalid circle radius or sample count.
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
    if radius <= 0.0 || n_samples < 8 {
        return Err("Invalid circle radius or sample count".to_string());
    }
    let mut min_margin = f64::INFINITY;
    for j in 0..n_samples {
        let theta = 2.0 * PI * (j as f64) / (n_samples as f64);
        let z = center + Complex::from_polar(radius, theta);
        let margin = f(z).norm() - g(z).norm();
        min_margin = min_margin.min(margin);
    }
    Ok(min_margin > 0.0)
}

/// Tracks a complex homotopy root from `t = 0` to `t = 1`.
///
/// Starts at `z(0) = z_0` for `H(z(t), t) = 0`, using an adaptive RK4 predictor and Newton
/// corrector with damping near a singular target.
///
/// # Errors
///
/// Returns an error for a zero step count, a singular Jacobian, or failure to reach and correct the endpoint within the step budget.
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
    if n_steps == 0 {
        return Err("n_steps must be positive".to_string());
    }
    let mut dt = 1.0 / (n_steps as f64);
    let min_dt = 1e-6;
    let mut t = 0.0_f64;
    let mut z = z0;
    let mut guard = 0_usize;
    while t < 1.0 - 1e-14 {
        guard += 1;
        if guard > 50_000 {
            return Err("Homotopy step budget exceeded".to_string());
        }
        let step = dt.min(1.0 - t);
        let t_next = t + step;

        // RK4 predictor for dz/dt = -(dH/dt)/(dH/dz)
        let vel = |zz: Complex<f64>, tt: f64| -> Option<Complex<f64>> {
            let hz = dh_dz(zz, tt);
            if hz.norm() < 1e-14 {
                None
            } else {
                Some(-dh_dt(zz, tt) / hz)
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
        for _ in 0..newton_iters.max(6) {
            let val = h(z_corr, t_next);
            if val.norm() < tol {
                converged = true;
                break;
            }
            let dval = dh_dz(z_corr, t_next);
            if dval.norm() < 1e-14 {
                break;
            }
            z_corr -= val / dval;
        }
        if !converged && t_next < 1.0 - 1e-6 && step > min_dt {
            dt *= 0.5;
            continue;
        }
        z = z_corr;
        t = t_next;
        if converged && dt < 0.1 {
            dt = (dt * 1.25).min(0.1);
        }
    }
    if h(z, 1.0).norm() > tol.sqrt().max(1e-6) {
        return Err("Homotopy corrector failed to converge at t = 1".to_string());
    }
    Ok(z)
}

/// Verified root isolation for a complex polynomial inside `|Re(z)|, |Im(z)| <= R`
/// using [`ComplexBall`](crate::kernels::complex_ball::ComplexBall) interval exclusion,
/// quadtree subdivision, and 1D Krawczyk certification.
#[must_use]
pub fn isolate_roots_verified(
    coeffs: &[Complex<f64>],
    tol: f64,
) -> Vec<crate::kernels::complex_ball::ComplexBall> {
    use crate::kernels::complex_ball::ComplexBall;
    use crate::kernels::complex_ball::ComplexBallJet;
    use crate::kernels::complex_ball::certify_root_1d;

    let Some(bound) = cauchy_root_bound(coeffs) else {
        return Vec::new();
    };
    let r = bound * 1.05;
    let mut stack = vec![(Complex::new(0.0, 0.0), r * 2.0_f64.sqrt())];
    let mut certified: Vec<ComplexBall> = Vec::new();
    let mut iters = 0_usize;
    while let Some((center, rad)) = stack.pop() {
        iters += 1;
        if iters > 20_000 {
            break;
        }
        let Some(ball) = ComplexBall::new(center, rad) else {
            continue;
        };
        // Evaluate polynomial over ball via Horner
        let mut p_ball = ComplexBall::zero();
        for &c in coeffs.iter().rev() {
            p_ball = p_ball * ball + ComplexBall::of(c);
        }
        if !p_ball.contains_zero() {
            continue;
        }
        if let Some(cert) = certify_root_1d(
            |z_jet| {
                let mut acc = ComplexBallJet::constant(ComplexBall::zero());
                for &c in coeffs.iter().rev() {
                    acc = acc
                        .mul(&z_jet)
                        .add(&ComplexBallJet::constant(ComplexBall::of(c)));
                }
                Some(acc)
            },
            ball,
        ) {
            if !certified.iter().any(|existing| existing.overlaps(&cert)) {
                certified.push(cert);
            }
            continue;
        }
        if rad <= tol.max(1e-9) {
            continue;
        }
        let half = 0.5 * rad;
        let offset = 0.5 * rad / 2.0_f64.sqrt();
        for (sx, sy) in [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0)] {
            stack.push((
                Complex::new(center.re + sx * offset, center.im + sy * offset),
                half * 1.05,
            ));
        }
    }
    certified
}
