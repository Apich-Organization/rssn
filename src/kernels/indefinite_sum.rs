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
//!    `F(b) - F(a) = ∫_a^b f - ½ (f(b) - f(a)) +
//!    Σ B_{2j}/(2j)! (f^{(2j-1)}(b) - f^{(2j-1)}(a))` evaluated at `a = N`,
//!    `b = x + N`. It is accepted when two different shifts agree.
//! 2. **Taylor series**: `f` is expanded about the middle of the unit
//!    interval after `h`, and the antidifference of each power is a Bernoulli
//!    polynomial, `Δ⁻¹ u^m = (B_{m+1}(u) - B_{m+1}(0)) / (m + 1)`. It is exact
//!    for polynomials and entire functions such as exponentials, and is
//!    reduced to the unit interval by exact integer sums.
//!
//! The [`PolynomialAntidifference`] compiles the polynomial formula for repeated
//! evaluations and nearby intervals. Indefinite products follow from
//! `Π f = exp Σ ln f`.

use super::exact_eq;
use super::numeric_common::scaled_difference;
use crate::kernels::integrate::gauss_legendre_quadrature;
use crate::kernels::special::bernoulli_number;

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
    /// Number of Taylor terms of the series strategy, at most 64 when used.
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

/// A polynomial antidifference compiled for repeated evaluation.
///
/// For the input polynomial `P(x) = Σ coefficients[m] (x - center)^m`,
/// this stores the coefficients of the unique polynomial `F` satisfying
/// `F(x + 1) - F(x) = P(x)` and `F(center) = 0`. Construction takes
/// quadratic work in the degree, plus the cost of computing the required
/// Bernoulli numbers. Each subsequent evaluation takes linear work.
///
/// The power basis can be ill-conditioned at high degree or far from the
/// center. Use finite arguments and coefficients at representable scales.
#[derive(Debug, Clone)]
pub struct PolynomialAntidifference {
    center: f64,
    // The input constant contributes a separate linear term of F.
    constant: f64,
    // Ascending coefficients of the remaining part of F(center + u) / u.
    polynomial: Vec<f64>,
}

impl PolynomialAntidifference {
    /// Compile the antidifference of the given ascending power coefficients.
    #[must_use]
    pub fn new(
        coefficients: &[f64],
        center: f64,
    ) -> Self {
        let mut polynomial = vec![0.0; coefficients.len()];
        compile_antidifference(coefficients, &mut polynomial);
        Self {
            center,
            constant: coefficients.first().copied().unwrap_or(0.0),
            polynomial,
        }
    }

    /// Evaluate the antidifference at one point.
    #[must_use]
    pub fn evaluate(
        &self,
        x: f64,
    ) -> f64 {
        evaluate_antidifference(self.constant, &self.polynomial, x - self.center)
    }

    /// Evaluate `F(to) - F(from)` without subtracting two rounded values.
    #[must_use]
    pub fn difference(
        &self,
        from: f64,
        to: f64,
    ) -> f64 {
        self.difference_of_offsets(from - self.center, to - self.center, to - from)
    }

    pub(crate) fn difference_of_offsets(
        &self,
        start: f64,
        end: f64,
        width: f64,
    ) -> f64 {
        let Some((&leading, tail)) = self.polynomial.split_last() else {
            return 0.0;
        };
        if exact_eq(width.abs(), 1.0)
            && normalized_unit_step(start, end)
            && self.polynomial.iter().all(|c| c.is_finite())
        {
            return width * self.constant;
        }
        let mut value = leading;
        let mut slope = 0.0;
        // Divided-difference Horner evaluation, including F's zero constant.
        for &c in tail.iter().rev() {
            slope = value + start * slope;
            value = c + end * value;
        }
        width * self.constant + width * (value + start * slope)
    }

    pub(crate) fn coefficients_are_finite(&self) -> bool {
        self.constant.is_finite()
            && self
                .polynomial
                .iter()
                .all(|coefficient| coefficient.is_finite())
    }
}

/// Antidifference of the power series `Σ coefficients[m] (t - center)^m`,
/// normalised to zero at `t = center`. For repeated evaluation or differences
/// between nearby points, use [`PolynomialAntidifference`].
#[must_use]
pub fn series_antidiff(
    coefficients: &[f64],
    center: f64,
    x: f64,
) -> f64 {
    // Keep the usual Taylor expansion on the stack in the one-shot path.
    let mut scratch = [0.0; 16];
    let Some(polynomial) = scratch.get_mut(..coefficients.len()) else {
        return PolynomialAntidifference::new(coefficients, center).evaluate(x);
    };
    compile_antidifference(coefficients, polynomial);
    let constant = coefficients.first().copied().unwrap_or(0.0);
    evaluate_antidifference(constant, polynomial, x - center)
}

fn compile_antidifference(
    coefficients: &[f64],
    polynomial: &mut [f64],
) {
    // The constant input term has its own linear antidifference.
    for (index, (input, output)) in coefficients
        .iter()
        .zip(&mut *polynomial)
        .enumerate()
        .skip(1)
    {
        *output = input / (index as f64 + 1.0);
    }
    for (order, _) in coefficients.iter().enumerate().skip(1) {
        let Ok(k) = u32::try_from(order) else {
            polynomial.fill(f64::NAN);
            return;
        };
        let bernoulli = bernoulli_number(k);
        if bernoulli == 0.0 {
            continue;
        }
        // In (B_{m+1}(u) - B_{m+1}(0))/(m+1), the weight of u^j
        // for m=j+k-1 is B_k * binomial(j+k-1,k)/j. At j=1 that
        // is B_k; the ratio of consecutive weights is (j+k)/(j+1).
        let mut weight = bernoulli;
        for (index, (input, output)) in coefficients
            .iter()
            .skip(order)
            .zip(&mut *polynomial)
            .enumerate()
        {
            *output += input * weight;
            let j = index as f64 + 1.0;
            weight *= (j + f64::from(k)) / (j + 1.0);
        }
    }
}

fn normalized_unit_step(
    start: f64,
    end: f64,
) -> bool {
    (start == 0.0 && exact_eq(end, 1.0)) || (end == 0.0 && exact_eq(start, 1.0))
}

fn evaluate_antidifference(
    constant: f64,
    polynomial: &[f64],
    u: f64,
) -> f64 {
    let Some((&leading, tail)) = polynomial.split_last() else {
        return 0.0;
    };
    // F(center + 1) = P(center) exactly. Avoid letting a large higher-order
    // term that vanishes there obscure the constant input contribution.
    if exact_eq(u, 1.0) && polynomial.iter().all(|c| c.is_finite()) {
        return constant;
    }
    u * constant + u * tail.iter().rev().fold(leading, |value, &c| c + u * value)
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
    let d1 =
        |x: f64| (-g(x + 2.0 * H) + 8.0 * g(x + H) - 8.0 * g(x - H) + g(x - 2.0 * H)) / (12.0 * H);
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
    panels: usize,
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
    head + euler_maclaurin_difference(g, n, u + n, panels)
}

const DEFAULT_EVALUATION_BUDGET: usize = 5_000_000;
const MAX_TAYLOR_TERMS: usize = 64;
const MAX_EXACT_LATTICE_INDEX: u64 = 1 << f64::MANTISSA_DIGITS;

/// A conservative evaluation count for both shifts. The quadrature has
/// five nodes per panel; each shift samples two endpoint values and two
/// four-point derivative formulas at both endpoints.
fn euler_work(
    u: f64,
    shift: usize,
) -> Option<(usize, usize, usize)> {
    #[allow(clippy::cast_sign_loss)]
    let panels = usize::try_from((u.abs() * 4.0).ceil() as u64).ok()?.max(1);
    let second_shift = shift.checked_add(10)?;
    if u64::try_from(second_shift).ok()? > MAX_EXACT_LATTICE_INDEX {
        return None;
    }
    let head = shift.checked_add(second_shift)?.checked_mul(2)?;
    let cost = head.checked_add(panels.checked_mul(10)?)?.checked_add(36)?;
    Some((second_shift, panels, cost))
}

fn lattice_sum(
    g: &impl Fn(f64) -> f64,
    start: f64,
    n: i64,
) -> f64 {
    let sum = |range: std::ops::Range<i64>| -> f64 {
        #[allow(clippy::cast_precision_loss)]
        range.map(|k| g(start + k as f64)).sum()
    };
    if n >= 0 {
        sum(0..n)
    } else {
        -sum(n..0)
    }
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
    let base = PolynomialAntidifference::new(&coefficients, 0.5).difference(0.0, y);
    #[allow(clippy::cast_possible_truncation)]
    let n = whole as i64;
    base + lattice_sum(g, y, n)
}

/// Evaluates the indefinite sum `F(x) = Σ_{k=h}^{x-1} f(k)` (steps of
/// `config.step`, `F(h) = 0`) at a real `x`. For `x` an integer number of
/// steps from `h` this is the plain sum.
///
/// Up to 5,000,000 calls to `f` are permitted. Use
/// [`indefinite_sum_with_budget`] to set a different limit.
///
/// # Errors
///
/// See [`indefinite_sum_with_budget`] for errors.
pub fn indefinite_sum(
    f: impl Fn(f64) -> f64,
    x: f64,
    config: &IndefiniteSumConfig,
) -> Result<f64, String> {
    indefinite_sum_with_budget(f, x, config, DEFAULT_EVALUATION_BUDGET)
}

/// Evaluate the indefinite sum with an explicit callback budget.
///
/// The exact lattice path takes one callback per step. For a fractional
/// offset, both Euler–Maclaurin shifts are reserved together and Taylor
/// fallback may use the remaining budget. If the shifts exceed the budget,
/// Taylor can use it in full. A finite normalization point itself needs no
/// calls. The callback receives only finite coordinates.
///
/// # Errors
///
/// Returns an error for non-finite endpoints, zero or non-finite step, or an
/// offset outside the exact floating-point lattice range. It also errors if
/// no strategy fits the budget or can produce a finite sum. A Taylor fallback
/// that is needed requires at most 64 terms.
pub fn indefinite_sum_with_budget(
    f: impl Fn(f64) -> f64,
    x: f64,
    config: &IndefiniteSumConfig,
    max_evaluations: usize,
) -> Result<f64, String> {
    if config.step == 0.0 || !config.step.is_finite() {
        return Err("step must be finite and non-zero".to_string());
    }
    if !x.is_finite() || !config.h.is_finite() {
        return Err("indefinite sum: endpoints must be finite".to_string());
    }
    let u = scaled_difference(x, config.h, config.step);
    if !u.is_finite() || u.abs() > MAX_EXACT_LATTICE_INDEX as f64 {
        return Err("indefinite sum: offset exceeds the exact Float lattice range".to_string());
    }
    let coordinate = |u: f64| config.step.mul_add(u, config.h);
    let nearest = u.round();
    if (u - nearest).abs() < 1e-12 {
        #[allow(clippy::cast_possible_truncation)]
        let n = nearest as i64;
        if usize::try_from(n.unsigned_abs())
            .ok()
            .is_none_or(|cost| cost > max_evaluations)
        {
            return Err(
                "indefinite sum: evaluation budget is too small for the interval".to_string(),
            );
        }
        // The physical lattice is monotone, so a finite outermost sample
        // guarantees the intermediate samples are finite as well.
        let last = if n > 0 { n - 1 } else { n };
        if last != 0 && !coordinate(last as f64).is_finite() {
            return Err("indefinite sum: lattice sample is outside the Float range".to_string());
        }
        let value = lattice_sum(&|u| f(coordinate(u)), 0.0, n);
        return if value.is_finite() {
            Ok(value)
        } else {
            Err("f is not finite on the sum".to_string())
        };
    }
    // Every sample is at most the larger shift plus the span and a margin
    // for endpoint derivatives. Double the broad bound to cover rounding
    // in composite-panel coordinates. A finite absolute affine bound lets
    // all callbacks on the common path run without a per-sample guard.
    let extent = 2.0 * (u.abs() + config.shift as f64 + 12.0);
    if config
        .step
        .abs()
        .mul_add(extent, config.h.abs())
        .is_finite()
    {
        fractional_sum(&|u| f(coordinate(u)), u, config, max_evaluations)
    } else {
        fractional_sum(
            &|u| {
                let x = coordinate(u);
                if x.is_finite() {
                    f(x)
                } else {
                    f64::NAN
                }
            },
            u,
            config,
            max_evaluations,
        )
    }
}

fn fractional_sum(
    g: &impl Fn(f64) -> f64,
    u: f64,
    config: &IndefiniteSumConfig,
    max_evaluations: usize,
) -> Result<f64, String> {
    let mut remaining = max_evaluations;
    if let Some((second_shift, panels, cost)) = euler_work(u, config.shift)
        && let Some(rest) = remaining.checked_sub(cost)
    {
        remaining = rest;
        let tolerance = config.tolerance.max(1e-14);
        let first = shifted_euler_maclaurin(g, u, config.shift, panels);
        let second = shifted_euler_maclaurin(g, u, second_shift, panels);
        if first.is_finite()
            && second.is_finite()
            && (first - second).abs() <= tolerance * first.abs().max(1.0)
        {
            return Ok(second);
        }
    }
    if config.taylor_terms > MAX_TAYLOR_TERMS {
        return Err("indefinite sum: Taylor fallback requires at most 64 terms".to_string());
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let steps = u.floor().abs() as u64;
    let taylor_cost = usize::try_from(steps)
        .ok()
        .and_then(|steps| steps.checked_add(config.taylor_terms));
    if taylor_cost.is_none_or(|cost| cost > remaining) {
        return Err("indefinite sum: evaluation budget is too small for the interval".to_string());
    }
    let taylor = taylor_strategy(g, u, config.taylor_terms);
    if taylor.is_finite() {
        Ok(taylor)
    } else {
        Err("indefinite sum: no strategy produced a finite value".to_string())
    }
}

/// Evaluates the indefinite product `F(x) = Π_{k=h}^{x-1} f(k)` (`F(h) = 1`)
/// as `exp Σ ln f`, for a real `x`. Up to 5,000,000 calls to `f` are permitted.
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
    indefinite_product_with_budget(f, x, config, DEFAULT_EVALUATION_BUDGET)
}

/// Evaluate the indefinite product with an explicit callback budget.
///
/// # Errors
///
/// Returns the errors of [`indefinite_sum_with_budget`] or an error if `f`
/// is not positive and finite where needed.
pub fn indefinite_product_with_budget(
    f: impl Fn(f64) -> f64,
    x: f64,
    config: &IndefiniteSumConfig,
    max_evaluations: usize,
) -> Result<f64, String> {
    let log_f = |t: f64| {
        let value = f(t);
        if value > 0.0 {
            value.ln()
        } else {
            f64::NAN
        }
    };
    let log_sum = indefinite_sum_with_budget(log_f, x, config, max_evaluations)
        .map_err(|e| format!("indefinite product needs f > 0: {e}"))?;
    Ok(log_sum.exp())
}
