//! Numerical methods for the theory of univalent functions.
//!
//! Includes Bieberbach–de Branges checks, the Fekete–Szegő functional, Hankel determinants,
//! Koebe growth, distortion and rotation bounds, Schwarzian derivatives, the Nehari criterion,
//! Grunsky matrices, and the Loewner slit equation.

use crate::kernels::complex::circle_spectrum;
use crate::kernels::complex_roots::validate_circle;
use crate::kernels::complex_scalar::Complex as ScaledComplex;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::fft_scratch_for;
use crate::kernels::reuse::global_pool;

use num_complex::Complex;
use std::f64::consts::PI;

fn stable_divide(
    value: Complex<f64>,
    divisor: Complex<f64>,
) -> Complex<f64> {
    let quotient =
        ScaledComplex::new(value.re, value.im).div(ScaledComplex::new(divisor.re, divisor.im));
    Complex::new(quotient.re, quotient.im)
}

/// Sharp radius of starlikeness for the class $\mathcal{S}$: $\tanh(\pi / 4) \approx 0.65579$.
pub const S_RADIUS_OF_STARLIKENESS: f64 = 0.655_794_202_632_672_4;

/// Sharp radius of convexity for the class $\mathcal{S}$: $2 - \sqrt{3} \approx 0.26795$.
pub const S_RADIUS_OF_CONVEXITY: f64 = 0.267_949_192_431_122_7;

/// Checks the Bieberbach–de Branges necessary condition $|`a_n` / `a_1`| \le n$ ($n \ge 2$)
/// for a normalized univalent function $f(z) = `a_1` z + `a_2` z^2 + \dots$ (with `coeffs[k]` = $`a_k`$).
#[must_use]
pub fn check_bieberbach_de_branges(
    taylor_coeffs: &[Complex<f64>],
    tol: f64,
) -> bool {
    if !tol.is_finite() || tol < 0.0 || taylor_coeffs.iter().any(|z| !z.is_finite()) {
        return false;
    }
    let Some(&a1) = taylor_coeffs.get(1).filter(|z| **z != Complex::default()) else {
        return false;
    };
    taylor_coeffs
        .iter()
        .enumerate()
        .skip(2)
        .all(|(n, &coefficient)| stable_divide(coefficient, a1).norm() <= n as f64 + tol)
}

/// Computes the Fekete–Szegő functional $|`a_3` - \mu `a_2^2`|$.
#[must_use]
pub fn fekete_szego_functional(
    a2: Complex<f64>,
    a3: Complex<f64>,
    mu: Complex<f64>,
) -> f64 {
    (a3 - mu * a2 * a2).norm()
}

/// Computes the $q \times q$ Hankel determinant $`H_q(n)` = \det [c_{n + i + j}]_{i,j=0}^{q-1}$.
///
/// # Errors
///
/// Returns an error if the supplied series has too few coefficients for the requested determinant.
pub fn hankel_determinant(
    coeffs: &[Complex<f64>],
    q: usize,
    n: usize,
) -> Result<Complex<f64>, String> {
    if q == 0 {
        return Ok(Complex::new(1.0, 0.0));
    }
    let last = q - 1;
    let needed = last.checked_mul(2).and_then(|degree| n.checked_add(degree));
    let Some(needed) = needed.filter(|&degree| degree < coeffs.len()) else {
        return Err("Insufficient coefficients for requested Hankel determinant".into());
    };
    if coeffs[n..=needed].iter().any(|value| !value.is_finite()) {
        return Err("Hankel coefficients must be finite".into());
    }
    let elements = q.checked_mul(q).ok_or("Hankel matrix size overflow")?;
    let mut matrix = fft_scratch_for(elements);
    for offset in 0..q {
        matrix.extend_from_slice(&coeffs[(n + offset)..(n + offset + q)]);
    }
    let mut determinant = Complex::from(1.0);
    for column in 0..q {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let pivot = (column..q)
            .max_by(|&i, &j| {
                matrix[i * q + column]
                    .norm()
                    .total_cmp(&matrix[j * q + column].norm())
            })
            .ok_or("Missing Hankel pivot")?;
        let diagonal = matrix[pivot * q + column];
        if diagonal == Complex::default() {
            return Ok(Complex::default());
        }
        if column != pivot {
            for j in column..q {
                matrix.swap(column * q + j, pivot * q + j);
            }
            determinant = -determinant;
        }
        determinant *= diagonal;
        let (above, below) = matrix.split_at_mut((column + 1) * q);
        let pivot_row = &above[(column * q)..];
        for row in below.chunks_exact_mut(q) {
            let factor = stable_divide(row[column], diagonal);
            for (cell, &pivot_value) in row[(column + 1)..]
                .iter_mut()
                .zip(&pivot_row[(column + 1)..])
            {
                *cell -= factor * pivot_value;
            }
        }
    }
    if determinant.is_finite() {
        Ok(determinant)
    } else {
        Err("Hankel determinant exceeds Float range".into())
    }
}

/// Koebe growth theorem bounds for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $\frac{r}{(1+r)^2} \le |f(z)| \le \frac{r}{(1-r)^2}$.
#[must_use]
pub fn koebe_growth_bounds(r: f64) -> Option<(f64, f64)> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let lo = r / ((1.0 + r) * (1.0 + r));
    let hi = r / ((1.0 - r) * (1.0 - r));
    Some((lo, hi))
}

/// Koebe distortion theorem bounds for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $\frac{1-r}{(1+r)^3} \le |f'(z)| \le \frac{1+r}{(1-r)^3}$.
#[must_use]
pub fn koebe_distortion_bounds(r: f64) -> Option<(f64, f64)> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let lo = (1.0 - r) / (1.0 + r).powi(3);
    let hi = (1.0 + r) / (1.0 - r).powi(3);
    Some((lo, hi))
}

/// Returns Goluzin’s sharp rotation bound on the circle of radius `r`.
///
/// Goluzin's sharp rotation theorem bound on $|\arg f'(z)|$ for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $4 \arcsin(r)$ for $r \le 1/\sqrt{2}$, and $\pi + \ln\frac{r^2}{1 - r^2}$ for $1/\sqrt{2} < r < 1$.
#[must_use]
pub fn koebe_rotation_bound(r: f64) -> Option<f64> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let thresh = std::f64::consts::FRAC_1_SQRT_2;
    if r <= thresh {
        Some(4.0 * r.asin())
    } else {
        let complement = (-r).mul_add(r, 1.0);
        let twice_square_minus_one = (r + r).mul_add(r, -1.0);
        Some(PI + (twice_square_minus_one / complement).ln_1p())
    }
}

pub(crate) const SCHWARZIAN_SAMPLES: usize = 64;
pub(crate) const RADIUS_ITERATIONS: usize = 40;

/// Computes the Schwarzian derivative $\{f, z\} = \frac{f'''(z)}{f'(z)} - \frac{3}{2} \left(\frac{f''(z)}{f'(z)}\right)^2$
/// using Cauchy–FFT derivatives on a circle of radius `h`.
///
/// # Errors
///
/// Returns an error when derivative evaluation fails or the first derivative vanishes.
pub fn schwarzian_derivative_numeric<F>(
    f: &F,
    z: Complex<f64>,
    h: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    validate_circle(z, h, SCHWARZIAN_SAMPLES)?;
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let spectrum = circle_spectrum(f, z, h, SCHWARZIAN_SAMPLES);
    let [f0, f1, f2, f3, ..] = spectrum.as_slice() else {
        return Err("Schwarzian spectrum lost four low coefficients".into());
    };
    if ![f0, f1, f2, f3].into_iter().all(|value| value.is_finite()) {
        return Err("Schwarzian callback samples and Fourier coefficients must be finite".into());
    }
    let scale = |v: &Complex<f64>| v.re.abs().max(v.im.abs());
    if scale(f1) <= (64.0 * f64::EPSILON) * scale(f0).max(scale(f2)).max(scale(f3)) {
        return Err(
            "Critical or numerically unresolved first derivative in Schwarzian stencil".into(),
        );
    }
    let second = stable_divide(*f2, *f1);
    let third = stable_divide(*f3, *f1);
    // Both derivative ratios carry two inverse powers of the circle radius here.
    let value = (((third - second * second) * 6.0) / h) / h;
    if value.is_finite() {
        Ok(value)
    } else {
        Err("Schwarzian derivative exceeds Float range".into())
    }
}

/// Checks Nehari's sufficient univalence criterion $(1 - |z|^2)^2 |\{f, z\}| \le 2$
/// across sample radii in $(0, 1)$ and `n_angles` angles per radius.
///
/// # Errors
///
/// Returns an error for invalid sample radii, a zero angular sample count, or failed Schwarzian evaluation.
pub fn check_nehari_univalence<F>(
    f: &F,
    sample_radii: &[f64],
    n_angles: usize,
    h: f64,
) -> Result<bool, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    if n_angles == 0 || sample_radii.is_empty() || !h.is_finite() || h <= 0.0 {
        return Err(
            "Nehari grid needs nonempty angle/radius samples and a positive finite step".into(),
        );
    }
    for &r in sample_radii {
        if !(0.0..1.0).contains(&r) || r + h >= 1.0 {
            return Err("Every Nehari sampling circle must fit in the open unit disk".into());
        }
    }
    let directions = AngularGrid::new(n_angles);
    for &r in sample_radii {
        let weight = (1.0 - r) * (1.0 + r);
        for direction in directions.directions() {
            let z = r * direction;
            let value = schwarzian_derivative_numeric(f, z, h)?;
            if weight * weight * value.norm() > 2.0 + 1e-9 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

enum AngularGrid {
    Cached(std::sync::Arc<[Complex<f64>]>),
    Local(Scratch<'static, Complex<f64>>),
}

impl AngularGrid {
    fn new(count: usize) -> Self {
        if count.is_power_of_two() && count <= 65_536 {
            Self::Cached(global_pool().unit_roots(count))
        } else {
            let mut storage = fft_scratch_for(count);
            storage.extend((0..count).map(|index| {
                Complex::from_polar(1.0, std::f64::consts::TAU * (index as f64) / (count as f64))
                    .conj()
            }));
            Self::Local(storage)
        }
    }

    fn directions(&self) -> impl Iterator<Item = Complex<f64>> + '_ {
        let roots = match self {
            | Self::Cached(roots) => roots.as_ref(),
            | Self::Local(roots) => roots.as_slice(),
        };
        roots.iter().map(Complex::conj)
    }
}

fn radius_by_bisection(
    maximum: f64,
    count: usize,
    positive: impl Fn(Complex<f64>) -> Option<bool>,
) -> Option<f64> {
    if count == 0 || !(0.0..=1.0).contains(&maximum) {
        return None;
    }
    let directions = AngularGrid::new(count);
    let mut low = 0.0;
    let mut high = maximum;
    for _ in 0..RADIUS_ITERATIONS {
        crate::cancellation::checkpoint().ok()?;
        let middle = f64::midpoint(low, high);
        let mut holds = true;
        for direction in directions.directions() {
            if !positive(middle * direction)? {
                holds = false;
                break;
            }
        }
        if holds {
            low = middle;
        } else {
            high = middle;
        }
    }
    Some(low)
}

/// Numerically estimates the radius of starlikeness $\sup \{ r \in (0, r_{\max}] : \text{Re}(z f'(z)/f(z)) > 0 \text{ on } |z|=r \}$
/// via bisection.
pub fn radius_of_starlikeness_numeric<F, DF>(
    f: &F,
    df: &DF,
    r_max: f64,
    n_angles: usize,
) -> f64
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    radius_by_bisection(r_max, n_angles, |z| {
        let denominator = f(z);
        let numerator = df(z);
        if !denominator.is_finite() || !numerator.is_finite() {
            return None;
        }
        if denominator == Complex::default() {
            return Some(false);
        }
        let criterion = z * stable_divide(numerator, denominator);
        criterion.is_finite().then_some(criterion.re > 0.0)
    })
    .unwrap_or(f64::NAN)
}

/// Numerically estimates the radius of convexity $\sup \{ r \in (0, r_{\max}] : 1 + \text{Re}(z f''(z)/f'(z)) > 0 \text{ on } |z|=r \}$
/// via bisection.
pub fn radius_of_convexity_numeric<DF, DDF>(
    df: &DF,
    ddf: &DDF,
    r_max: f64,
    n_angles: usize,
) -> f64
where
    DF: Fn(Complex<f64>) -> Complex<f64>,
    DDF: Fn(Complex<f64>) -> Complex<f64>,
{
    radius_by_bisection(r_max, n_angles, |z| {
        let denominator = df(z);
        let numerator = ddf(z);
        if !denominator.is_finite() || !numerator.is_finite() {
            return None;
        }
        if denominator == Complex::default() {
            return Some(false);
        }
        let criterion = z * stable_divide(numerator, denominator);
        criterion.is_finite().then_some(1.0 + criterion.re > 0.0)
    })
    .unwrap_or(f64::NAN)
}

/// Computes the Grunsky matrix of a Taylor series.
///
/// Computes the $M \times M$ Grunsky matrix $G_{m,n} = \sqrt{mn} \, c_{m,n}$ for $1 \le m, n \le M$
/// from the Taylor coefficients $f(z) = z + `a_2` z^2 + \dots$ (`taylor_coeffs[k]` = $`a_k`$, $`a_1` = 1$),
/// defined via $\ln \frac{f(z) - f(\zeta)}{z - \zeta} = -\sum_{m,n \ge 1} c_{m,n} z^m \zeta^n$.
///
/// # Errors
///
/// Returns an error when the Taylor coefficients are missing, have a zero linear term, or do not reach the requested order.
pub fn grunsky_matrix(
    taylor_coeffs: &[Complex<f64>],
    order: usize,
) -> Result<Vec<Vec<Complex<f64>>>, String> {
    if order == 0 {
        return Ok(Vec::new());
    }
    let required = order
        .checked_mul(2)
        .and_then(|degree| degree.checked_add(2));
    let Some(required) = required.filter(|&needed| needed <= taylor_coeffs.len()) else {
        return Err("Grunsky series needs at least 2 * order + 2 Taylor coefficients".into());
    };
    let a1 = taylor_coeffs[1];
    if a1 == Complex::default() || taylor_coeffs[..required].iter().any(|v| !v.is_finite()) {
        return Err("Grunsky coefficients must be finite with a nonzero linear term".into());
    }
    let mut normalized = fft_scratch_for(required - 1);
    normalized.push(Complex::from(1.0));
    normalized.extend(
        taylor_coeffs[2..required]
            .iter()
            .map(|&v| stable_divide(v, a1)),
    );
    if normalized.iter().any(|v| !v.is_finite()) {
        return Err("Normalized Grunsky coefficient exceeds Float range".into());
    }

    // Q=(f(z)-f(w))/(z-w), with Q[p,q]=normalized[p+q] and Q[0,0]=1.
    // Solve Q*F=z*Q_z for F=z*(log Q)_z by increasing row and column.
    let width = order + 1;
    let elements = width
        .checked_mul(width)
        .ok_or("Grunsky matrix size overflow")?;
    let mut derivative = fft_scratch_for(elements);
    derivative.resize(elements, Complex::default());
    for m in 1..=order {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let (above, rest) = derivative.split_at_mut(m * width);
        let output = &mut rest[..width];
        let m_float = m as f64;
        for n in 0..width {
            let left: Complex<f64> = output[..n]
                .iter()
                .rev()
                .zip(&normalized[1..=n])
                .map(|(&coefficient, &weight)| coefficient * weight)
                .sum();
            let mut value = m_float * normalized[m + n] - left;
            for (p, row) in above.chunks_exact(width).rev().take(m - 1).enumerate() {
                let first = p + 1;
                let product: Complex<f64> = row[..=n]
                    .iter()
                    .rev()
                    .zip(&normalized[first..=(first + n)])
                    .map(|(&coefficient, &weight)| coefficient * weight)
                    .sum();
                value -= product;
            }
            if !value.is_finite() {
                return Err("Grunsky recurrence exceeds Float range".into());
            }
            output[n] = value;
        }
    }
    // F[m,n] = m*[log Q][m,n]. The mathematical result is symmetric.
    let mut result = vec![vec![Complex::default(); order]; order];
    for m in 1..=order {
        for n in m..=order {
            let value = -(n as f64 / m as f64).sqrt() * derivative[m * width + n];
            if !value.is_finite() {
                return Err("Grunsky matrix exceeds Float range".into());
            }
            result[m - 1][n - 1] = value;
            result[n - 1][m - 1] = value;
        }
    }
    Ok(result)
}

/// Verifies the Grunsky inequality $\left|\sum_{m,n=1}^M G_{m,n} `x_m` `x_n\right`| \le \sum_{m=1}^M |`x_m|^2` + \text{tol}$.
///
/// # Errors
///
/// Returns an error if the supplied vector length does not match the matrix dimension.
pub fn check_grunsky_inequality(
    grunsky: &[Vec<Complex<f64>>],
    test_vec: &[Complex<f64>],
    tol: f64,
) -> Result<bool, String> {
    let m = grunsky.len();
    if test_vec.len() != m {
        return Err("Vector length must match Grunsky matrix dimension".to_string());
    }
    if !tol.is_finite()
        || tol < 0.0
        || test_vec.iter().any(|v| !v.is_finite())
        || grunsky
            .iter()
            .any(|row| row.len() != m || row.iter().any(|v| !v.is_finite()))
    {
        return Err("Grunsky check needs a finite square matrix, finite test vector and nonnegative tolerance".into());
    }
    let vector_scale = test_vec.iter().fold(1.0_f64, |scale, value| {
        scale.max(value.re.abs()).max(value.im.abs())
    });
    let mut normalized = fft_scratch_for(m);
    normalized.extend(test_vec.iter().map(|&value| value / vector_scale));
    let norm2: f64 = normalized.iter().map(Complex::norm_sqr).sum();
    let quad: Complex<f64> = grunsky
        .iter()
        .zip(normalized.iter())
        .map(|(row, &xi)| {
            let product: Complex<f64> = row
                .iter()
                .zip(normalized.iter())
                .map(|(&value, &xj)| value * xj)
                .sum();
            xi * product
        })
        .sum();
    if !quad.is_finite() {
        return Ok(false);
    }
    if quad.norm() > norm2 + (tol / vector_scale) / vector_scale {
        return Ok(false);
    }
    // Also estimate spectral norm ||G||_2 via power iteration on G^* G
    if m > 0 {
        let mut v = fft_scratch_for(m);
        v.resize(m, Complex::new(1.0 / (m as f64).sqrt(), 0.0));
        let mut gv = fft_scratch_for(m);
        let mut ghgv = fft_scratch_for(m);
        ghgv.resize(m, Complex::default());
        let mut sigma_sq = 0.0_f64;
        for _ in 0..24 {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            gv.clear();
            gv.extend(grunsky.iter().map(|row| {
                row.iter()
                    .zip(v.iter())
                    .map(|(&entry, &vj)| entry * vj)
                    .sum::<Complex<f64>>()
            }));
            ghgv.fill(Complex::default());
            for (row, &gvi) in grunsky.iter().zip(gv.iter()) {
                for (entry, &gij) in ghgv.iter_mut().zip(row) {
                    *entry += gij.conj() * gvi;
                }
            }
            let norm = ghgv
                .iter()
                .map(num_complex::Complex::norm_sqr)
                .sum::<f64>()
                .sqrt();
            if !norm.is_finite() {
                return Ok(false);
            }
            if norm < 1e-30 {
                sigma_sq = 0.0;
                break;
            }
            sigma_sq = norm;
            for (slot, &value) in v.iter_mut().zip(ghgv.iter()) {
                *slot = value / norm;
            }
        }
        if sigma_sq.sqrt() > 1.0 + tol {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Advances one RK4 step of the radial Loewner slit ODE
/// $\frac{dw}{dt} = w \frac{e^{i\kappa} + w}{e^{i\kappa} - w}$ with constant driving `driving_kappa` over `dt`.
#[must_use]
pub fn loewner_slit_step(
    w: Complex<f64>,
    driving_kappa: f64,
    dt: f64,
) -> Complex<f64> {
    let u = Complex::from_polar(1.0, driving_kappa);
    if (u - w).norm() < 1e-12 {
        return w;
    }
    loewner_step_checked(w, driving_kappa, dt).unwrap_or(Complex::new(f64::NAN, f64::NAN))
}

pub(crate) fn loewner_step_checked(
    w: Complex<f64>,
    driving_kappa: f64,
    dt: f64,
) -> Option<Complex<f64>> {
    let u = Complex::from_polar(1.0, driving_kappa);
    let rhs = |z: Complex<f64>| {
        let divisor = u - z;
        if divisor.norm() < 1e-12 {
            return None;
        }
        let result = z * stable_divide(u + z, divisor);
        result.is_finite().then_some(result)
    };
    let half_step = 0.5 * dt;
    let k1 = rhs(w)?;
    let k2 = rhs(w + k1 * half_step)?;
    let k3 = rhs(w + k2 * half_step)?;
    let k4 = rhs(w + k3 * dt)?;
    let result = w + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
    result.is_finite().then_some(result)
}

/// Checks Becker's sufficient univalence criterion `(1 - |z|^2) * |z * f''(z)/f'(z)| <= 1`.
#[must_use]
pub fn check_becker_univalence(
    z: Complex<f64>,
    df: Complex<f64>,
    ddf: Complex<f64>,
) -> bool {
    let r2 = z.norm_sqr();
    if r2 >= 1.0
        || !z.is_finite()
        || !df.is_finite()
        || !ddf.is_finite()
        || df == Complex::default()
    {
        return false;
    }
    (1.0 - r2) * (z * stable_divide(ddf, df)).norm() <= 1.0
}
