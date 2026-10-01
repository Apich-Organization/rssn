//! # Indefinite sums and products
//!
//! Numerical anti-differences of a function `f` given as a closure. The
//! indefinite sum `F(x) = Σ_{k=h}^{x-1} f(k)` (Nörlund principal solution)
//! is the solution of `F(x + S) - F(x) = f(x)` with `F(h) = 0`, defined for
//! real, non-integer `x` as well. Two strategies are used:
//!
//! 1. **Euler–Maclaurin** after shifting by a large integer `N`, which
//!    turns the problem into a smooth one far from the start:
//!    `F(x) = Σ_{j<N} [f(j) - f(x+j)] + F(x+N) - F(N)` with the difference
//!    `F(b) - F(a) = ∫_a^b f - ½ (f(b) - f(a)) + Σ B_{2j}/(2j)! (f^{(2j-1)}(b)
//!    - f^{(2j-1)}(a))` evaluated at `a = N`, `b = x + N`. It is accepted
//!    when two different shifts agree.
//! 2. **Taylor series**: `f` is expanded about the middle of the unit
//!    interval after `h`, and the antidifference of each power is a Bernoulli
//!    polynomial, `Δ⁻¹ u^m = (B_{m+1}(u) - B_{m+1}(0)) / (m + 1)`. It is exact
//!    for polynomials and entire functions such as exponentials, and is
//!    reduced to the unit interval by exact integer sums.
//!
//! Indefinite products follow from `Π f = exp Σ ln f`.

use crate::kernels::integrate::gauss_legendre_quadrature;
use crate::kernels::special::bernoulli_poly;

/// Options of the indefinite sum.
#[derive(Debug, Clone)]
pub struct IndefiniteSumConfig {
    /// Normalisation point: `F(h) = 0`.
    pub h: f64,
    /// Step `S` of the recurrence `F(x + S) - F(x) = f(x)`; the sum runs over
    /// `h, h + S, ...`.
    pub step: f64,
    /// Relative agreement required between the two Euler–Maclaurin shifts.
    pub tolerance: f64,
    /// Smallest shift `N` (in steps) of the Euler–Maclaurin strategy; a
    /// second evaluation uses `N + 10`.
    pub shift: usize,
    /// Number of Taylor terms of the series strategy.
    pub taylor_terms: usize,
}

impl Default for IndefiniteSumConfig {
    fn default() -> Self {
        Self {
            h: 0.0,
            step: 1.0,
            tolerance: 1e-9,
            shift: 30,
            taylor_terms: 14,
        }
    }
}

/// Taylor coefficients `c_m` of `f` about `center`.
///
/// They satisfy `f(t) ≈ Σ c_m (t -
/// center)^m` and come from the interpolating polynomial of degree `terms - 1`
/// through Chebyshev nodes on `[center - radius, center + radius]`. This is
/// accurate for moderate `terms` (up to about 16), unlike repeated finite
/// differences, and converges to the Taylor series for analytic `f`.
#[must_use]
pub fn numeric_taylor_coefficients(
    f: impl Fn(f64) -> f64,
    center: f64,
    radius: f64,
    terms: usize,
) -> Vec<f64> {
    #[allow(clippy::cast_precision_loss)]
    let nodes: Vec<f64> = (0..terms)
        .map(|j| (std::f64::consts::PI * (2 * j + 1) as f64 / (2 * terms) as f64).cos())
        .collect();
    // Solve V a = y with V[j][m] = s_j^m, by Gauss elimination.
    let mut a: Vec<Vec<f64>> = nodes
        .iter()
        .map(|&s| {
            let mut row: Vec<f64> = (0..terms)
                .scan(1.0, |power, _| {
                    let current = *power;
                    *power *= s;
                    Some(current)
                })
                .collect();
            row.push(f(center + radius * s));
            row
        })
        .collect();
    for col in 0..terms {
        let pivot = (col..terms)
            .max_by(|&p, &q| a[p][col].abs().total_cmp(&a[q][col].abs()))
            .unwrap_or(col);
        a.swap(col, pivot);
        let diagonal = a[col][col];
        if diagonal == 0.0 {
            return vec![f64::NAN; terms];
        }
        let pivot_row = a[col].clone();
        for row in a.iter_mut().skip(col + 1) {
            let factor = row[col] / diagonal;
            for (value, p) in row.iter_mut().zip(&pivot_row).skip(col) {
                *value -= factor * p;
            }
        }
    }
    let mut scaled = vec![0.0; terms];
    for m in (0..terms).rev() {
        let tail: f64 = (m + 1..terms).map(|k| a[m][k] * scaled[k]).sum();
        scaled[m] = (a[m][terms] - tail) / a[m][m];
    }
    let mut scale = 1.0;
    scaled
        .into_iter()
        .map(|c| {
            let coefficient = c / scale;
            scale *= radius;
            coefficient
        })
        .collect()
}

/// Antidifference of the power series `Σ coefficients[m] (t - center)^m`,
/// normalised to zero at `t = center`: `Σ c_m (B_{m+1}(u) - B_{m+1}(0)) /
/// (m + 1)` with `u = x - center`. Exact for polynomials.
#[must_use]
pub fn series_antidiff(
    coefficients: &[f64],
    center: f64,
    x: f64,
) -> f64 {
    let u = x - center;
    coefficients
        .iter()
        .enumerate()
        .map(|(m, &c)| {
            #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
            let n = m as u32 + 1;
            c * (bernoulli_poly(n, u) - bernoulli_poly(n, 0.0)) / f64::from(n)
        })
        .sum()
}

/// Euler–Maclaurin difference `F(b) - F(a)` for unit step, with the first
/// and third derivatives of `g` at the endpoints taken by finite differences.
fn euler_maclaurin_difference(
    g: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    panels: usize,
) -> f64 {
    const H: f64 = 0.05;
    let d1 = |x: f64| (-g(x + 2.0 * H) + 8.0 * g(x + H) - 8.0 * g(x - H) + g(x - 2.0 * H)) / (12.0 * H);
    let d3 = |x: f64| {
        (g(x + 2.0 * H) - 2.0 * g(x + H) + 2.0 * g(x - H) - g(x - 2.0 * H)) / (2.0 * H.powi(3))
    };
    // Composite 5-point Gauss-Legendre on panels of width `(b - a) / panels`.
    #[allow(clippy::cast_precision_loss)]
    let width = (b - a) / panels as f64;
    let integral: f64 = (0..panels)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let left = a + width * i as f64;
            gauss_legendre_quadrature(g, (left, left + width))
        })
        .sum();
    // B_2 / 2! = 1/12, B_4 / 4! = -1/720.
    integral - 0.5 * (g(b) - g(a)) + (d1(b) - d1(a)) / 12.0 - (d3(b) - d3(a)) / 720.0
}

/// Unit-step indefinite sum from 0 by shifted Euler–Maclaurin.
fn shifted_euler_maclaurin(
    g: &impl Fn(f64) -> f64,
    u: f64,
    shift: usize,
) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let n = shift as f64;
    let head: f64 = (0..shift)
        .map(|j| {
            #[allow(clippy::cast_precision_loss)]
            let j = j as f64;
            g(j) - g(u + j)
        })
        .sum();
    // Panels of width at most a quarter of a step.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let panels = ((u.abs() * 4.0).ceil() as usize).max(1);
    head + euler_maclaurin_difference(g, n, u + n, panels)
}

/// Unit-step indefinite sum from 0 for `u` reduced by exact integer sums
/// to `[0, 1)`, then the Taylor/Bernoulli antidifference.
fn taylor_strategy(
    g: &impl Fn(f64) -> f64,
    u: f64,
    terms: usize,
) -> f64 {
    let whole = u.floor();
    let y = u - whole;
    let coefficients = numeric_taylor_coefficients(g, 0.5, 0.5, terms);
    let base = series_antidiff(&coefficients, 0.5, y) - series_antidiff(&coefficients, 0.5, 0.0);
    #[allow(clippy::cast_possible_truncation)]
    let n = whole as i64;
    if n >= 0 {
        base + (0..n).map(|j| {
            #[allow(clippy::cast_precision_loss)]
            let j = j as f64;
            g(y + j)
        }).sum::<f64>()
    } else {
        base - (n..0).map(|j| {
            #[allow(clippy::cast_precision_loss)]
            let j = j as f64;
            g(y + j)
        }).sum::<f64>()
    }
}

/// Evaluates the indefinite sum `F(x) = Σ_{k=h}^{x-1} f(k)` (steps of
/// `config.step`, `F(h) = 0`) at a real `x`. For `x` an integer number of
/// steps from `h` this is the plain sum.
///
/// # Errors
///
/// Returns an error if the step is zero or non-finite, if `f` is not finite
/// where it is needed, or if neither strategy gives a result.
pub fn indefinite_sum(
    f: impl Fn(f64) -> f64,
    x: f64,
    config: &IndefiniteSumConfig,
) -> Result<f64, String> {
    if config.step == 0.0 || !config.step.is_finite() {
        return Err("step must be finite and non-zero".to_string());
    }
    let g = |u: f64| f(config.step.mul_add(u, config.h));
    let u = (x - config.h) / config.step;
    if (u - u.round()).abs() < 1e-12 {
        #[allow(clippy::cast_possible_truncation)]
        let n = u.round() as i64;
        let sum = |range: std::ops::Range<i64>| -> f64 {
            #[allow(clippy::cast_precision_loss)]
            range.map(|k| g(k as f64)).sum()
        };
        let value = if n >= 0 { sum(0..n) } else { -sum(n..0) };
        return if value.is_finite() { Ok(value) } else { Err("f is not finite on the sum".to_string()) };
    }
    let tolerance = config.tolerance.max(1e-14);
    let first = shifted_euler_maclaurin(&g, u, config.shift);
    let second = shifted_euler_maclaurin(&g, u, config.shift + 10);
    if first.is_finite() && second.is_finite() && (first - second).abs() <= tolerance * first.abs().max(1.0) {
        return Ok(second);
    }
    let taylor = taylor_strategy(&g, u, config.taylor_terms);
    if taylor.is_finite() {
        Ok(taylor)
    } else {
        Err("indefinite sum: no strategy produced a finite value".to_string())
    }
}

/// Evaluates the indefinite product `F(x) = Π_{k=h}^{x-1} f(k)` (`F(h) = 1`)
/// as `exp Σ ln f`, for a real `x`.
///
/// # Errors
///
/// Returns an error if `f` is not positive and finite where needed, or if
/// the indefinite sum of `ln f` fails.
pub fn indefinite_product(
    f: impl Fn(f64) -> f64,
    x: f64,
    config: &IndefiniteSumConfig,
) -> Result<f64, String> {
    let log_f = |t: f64| {
        let value = f(t);
        if value > 0.0 { value.ln() } else { f64::NAN }
    };
    let log_sum = indefinite_sum(log_f, x, config)
        .map_err(|e| format!("indefinite product needs f > 0: {e}"))?;
    Ok(log_sum.exp())
}
