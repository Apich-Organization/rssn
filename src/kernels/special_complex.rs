//! # Special functions of a complex argument
//!
//! The real-argument routines of [`special`](super::special) cover the
//! gamma family, Bessel functions and Lambert W on the real line. This module
//! continues the most used of them into the complex plane:
//!
//! * [`ln_gamma_complex`], [`gamma_complex`], [`digamma_complex`]: Stirling's
//!   series for large arguments, recurrence shifts, and scaled reflection
//!   near the negative real axis.
//!   `ln_gamma_complex` is the analytic branch with
//!   `ln Γ(z + 1) = ln Γ(z) + ln z` (the one `SciPy` calls `loggamma`),
//!   continuous off the non-positive real axis.
//! * [`lambert_w_complex`]: every branch `k`, by Halley iteration from the
//!   starting values of Corless, Gonnet, Hare, Jeffrey and Knuth, with
//!   logarithmic iteration at extreme magnitudes.
//! * [`bessel_j_complex`], [`bessel_i_complex`]: real order `ν`, complex `z`;
//!   the ascending series for moderate `|z|` and Hankel's asymptotic
//!   expansion beyond it. Cancellation in the ascending series can reduce
//!   accuracy when `ν²` is close to `2 |z|`.
//!
//! No function panics; invalid input gives NaN components, poles give
//! infinity.

use std::f64::consts::PI;

use num_complex::Complex64;

const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const INF: Complex64 = Complex64::new(f64::INFINITY, 0.0);
const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

const fn is_nan(z: Complex64) -> bool {
    z.re.is_nan() || z.im.is_nan()
}

const fn is_inf(z: Complex64) -> bool {
    z.re.is_infinite() || z.im.is_infinite()
}

#[allow(clippy::float_cmp)] // exact integer test: poles sit exactly on the non-positive integers
fn is_pole(z: Complex64) -> bool {
    z.im == 0.0 && z.re <= 0.0 && z.re == z.re.floor()
}

/// `B_{2k} / (2k (2k - 1))` for `k = 1..=8`, the Stirling coefficients.
const STIRLING: [f64; 8] = [
    1.0 / 12.0,
    -1.0 / 360.0,
    1.0 / 1260.0,
    -1.0 / 1680.0,
    1.0 / 1188.0,
    -691.0 / 360_360.0,
    1.0 / 156.0,
    -3617.0 / 122_400.0,
];

/// `ln Γ(z)` after shifting right or for a large imaginary part.
fn stirling(z: Complex64) -> Complex64 {
    let zi = scaled_reciprocal(z);
    let z2 = zi * zi;
    let mut term = zi;
    let mut sum = ZERO;
    for c in STIRLING {
        sum += term * c;
        term *= z2;
    }
    let coefficient = scaled_log(z) - 1.0;
    let scale = z.re.abs().max(z.im.abs());
    let leading = if scale > f64::MAX / 1024.0 {
        let mut leading = (((z - 0.5) / scale) * coefficient) * scale;
        if z.re > 0.0 && coefficient.im.abs() < f64::MIN_POSITIVE {
            // Here arg(z) ≈ Im(z)/Re(z), a ratio that can underflow.
            leading.im = z.im * (coefficient.re + 1.0);
        }
        leading
    } else {
        (z - 0.5) * coefficient
    };
    leading + (0.5 * (2.0 * PI).ln() - 0.5) + sum
}

/// Sine and cosine of `πx`, with `x` already reduced to `[-1/2, 1/2]`.
fn sin_cos_pi_fraction(x: f64) -> (f64, f64) {
    if x.abs() > 0.25 {
        let (cos, sin) = (PI * (0.5 - x.abs())).sin_cos();
        (sin.copysign(x), cos)
    } else {
        (PI * x).sin_cos()
    }
}

/// Returns the scaled sine of `π(z - n)`, its logarithmic scale, and `n`,
/// the nearest integer to `Re z`.
fn scaled_sin_pi(z: Complex64) -> (Complex64, f64, f64) {
    let n = z.re.round();
    let (sin, cos) = sin_cos_pi_fraction(z.re - n);
    let b = PI * z.im;
    if b.abs() < 1.0 {
        (Complex64::new(sin * b.cosh(), cos * b.sinh()), 0.0, n)
    } else {
        let q = (-2.0 * b.abs()).exp();
        let direction = Complex64::new(sin * (1.0 + q), (cos * (1.0 - q)).copysign(b));
        (direction, b.abs() - std::f64::consts::LN_2, n)
    }
}

/// The analytic branch of `ln Γ(z)`, defined and continuous on the plane cut
/// along the non-positive real axis, with `ln Γ(z + 1) = ln Γ(z) + ln z`.
///
/// Non-positive integers (the poles) give `+∞`.
#[must_use]
pub fn ln_gamma_complex(z: Complex64) -> Complex64 {
    if is_nan(z) {
        return NAN;
    }
    if is_inf(z) {
        return if z.re == f64::INFINITY && z.im == 0.0 {
            INF
        } else {
            NAN
        };
    }
    if is_pole(z) {
        return INF;
    }
    // Far from the real cut, the exponentially small reflection correction
    // is negligible and Stirling avoids subtracting large reflected terms.
    if z.im.abs() >= 64.0 {
        return stirling(z);
    }
    if z.re < 0.5 {
        let (sine, scale, n) = scaled_sin_pi(z);
        let log_sine = scaled_log(sine);
        // Moving one interval to the right decreases the continuous sine
        // phase by π above the cut, and increases it by π below the cut.
        let log_sine = Complex64::new(log_sine.re + scale, log_sine.im - PI.copysign(z.im) * n);
        return PI.ln() - ln_gamma_complex(1.0 - z) - log_sine;
    }
    // Shift to Re >= 12, accumulating principal logarithms.
    let mut acc = ZERO;
    let mut w = z;
    while w.re < 12.0 {
        acc += scaled_log(w);
        w += 1.0;
    }
    stirling(w) - acc
}

/// `Γ(z)`; poles (non-positive integers) give `+∞`.
#[must_use]
pub fn gamma_complex(z: Complex64) -> Complex64 {
    if is_nan(z) {
        return NAN;
    }
    if is_pole(z) {
        return INF;
    }
    if z.im == 0.0 && z.re.abs() <= 170.0 {
        // The dedicated real routine is exact on the real axis.
        return Complex64::new(crate::kernels::special::gamma_numerical(z.re), 0.0);
    }
    if is_inf(z) {
        return if z.re == f64::INFINITY && z.im == 0.0 {
            INF
        } else {
            NAN
        };
    }
    if z.re < 0.5 && z.im.abs() < 64.0 {
        // Reflect the magnitude in log space and retain the sine direction
        // directly, including tiny imaginary parts near the real axis.
        let (sine, scale, n) = scaled_sin_pi(z);
        let reflected = ln_gamma_complex(1.0 - z);
        let magnitude_log = PI.ln() - scale - scaled_log(sine).re - reflected.re;
        if magnitude_log < -746.0 {
            return ZERO;
        }
        let direction = sine / sine.re.abs().max(sine.im.abs());
        let sign = if n % 2.0 == 0.0 { 1.0 } else { -1.0 };
        let phase = sign * unit_complex(-reflected.im) * direction.conj() / direction.norm();
        return scale_real_exponential(phase, magnitude_log);
    }
    let log = ln_gamma_complex(z);
    if log.re < -746.0 {
        return ZERO;
    }
    scale_real_exponential(unit_complex(log.im), log.re)
}

/// Principal logarithm with a scaled magnitude near the range limits.
fn scaled_log(z: Complex64) -> Complex64 {
    let scale = z.re.abs().max(z.im.abs());
    let log_abs = if (f64::MIN_POSITIVE..=f64::MAX / 2.0).contains(&scale) {
        z.norm().ln()
    } else {
        scale.ln() + (z / scale).norm().ln()
    };
    Complex64::new(log_abs, z.arg())
}

/// Reciprocal for a finite nonzero argument, without squaring its magnitude.
fn scaled_reciprocal(z: Complex64) -> Complex64 {
    let scale = z.re.abs().max(z.im.abs());
    if (1e-150..=1e150).contains(&scale) {
        return 1.0 / z;
    }
    let normalized = z / scale;
    normalized.conj() / normalized.norm_sqr() / scale
}

fn unit_complex(phase: f64) -> Complex64 {
    let (sin, cos) = phase.sin_cos();
    Complex64::new(cos, sin)
}

/// Apply a real exponential to a complex prefactor, scaling each component
/// before it reaches the range limits of `f64`.
fn scale_real_exponential(
    factor: Complex64,
    exponent: f64,
) -> Complex64 {
    if exponent.abs() < 700.0 {
        return factor * exponent.exp();
    }
    let half_scale = (0.5 * exponent).exp();
    let scale = |component: f64| {
        if component == 0.0 {
            component
        } else if half_scale.is_infinite() {
            (exponent + component.abs().ln()).exp().copysign(component)
        } else {
            component * half_scale * half_scale
        }
    };
    Complex64::new(scale(factor.re), scale(factor.im))
}

/// `cot(π z)` for a finite argument away from its poles.
fn cot_pi(z: Complex64) -> Complex64 {
    let (sin, cos) = sin_cos_pi_fraction(z.re - z.re.round());
    let b = PI * z.im;
    if b.abs() >= 1.0 {
        // Divide both numerator and denominator by exp(2 |b|).
        let q = (-2.0 * b.abs()).exp();
        let denominator = 1.0 + q * q - 2.0 * q * (cos * cos - sin * sin);
        return Complex64::new(4.0 * q * sin * cos, -z.im.signum() * (1.0 - q * q)) / denominator;
    }
    let sinh = b.sinh();
    let scale = sin.abs().max(sinh.abs());
    let (s, h) = (sin / scale, sinh / scale);
    Complex64::new(s * cos, -h * b.cosh()) / (s * s + h * h) / scale
}

/// The digamma function `ψ(z) = Γ'(z) / Γ(z)`; poles give NaN.
#[must_use]
pub fn digamma_complex(z: Complex64) -> Complex64 {
    // ψ(w) ~ ln w - 1/(2w) - Σ B_{2k} / (2k w^{2k})
    const B: [f64; 7] = [
        1.0 / 12.0,
        -1.0 / 120.0,
        1.0 / 252.0,
        -1.0 / 240.0,
        1.0 / 132.0,
        -691.0 / 32_760.0,
        1.0 / 12.0,
    ];
    if is_nan(z) || is_inf(z) || is_pole(z) {
        return NAN;
    }
    if z.re < 0.0 {
        // ψ(z) = ψ(1 - z) - π cot(π z)
        return digamma_complex(1.0 - z) - PI * cot_pi(z);
    }
    let mut acc = ZERO;
    let mut w = z;
    while w.re < 12.0 {
        acc += scaled_reciprocal(w);
        w += 1.0;
    }
    let wi = scaled_reciprocal(w);
    let w2 = wi * wi;
    let mut term = w2;
    let mut sum = ZERO;
    for c in B {
        sum += term * c;
        term *= w2;
    }
    scaled_log(w) - 0.5 * wi - sum - acc
}

/// Solve `w + log(w) = log(z) + 2πik` without forming either side of `w exp(w) = z`.
fn lambert_w_logarithmic(target: Complex64) -> Complex64 {
    let mut w = target - target.ln();
    for _ in 0..12 {
        let mut log_w = w.ln();
        if w.re < 0.0 && log_w.im.is_sign_negative() != target.im.is_sign_negative() {
            // Keep the same lip of log's cut when an iterate crosses the real axis.
            log_w.im += std::f64::consts::TAU.copysign(target.im);
        }
        let step = (w + log_w - target) / (1.0 + 1.0 / w);
        w -= step;
        if step.norm() <= f64::EPSILON * w.norm() {
            break;
        }
    }
    w
}

/// The branch `k` of the Lambert W function: the solution of `w e^w = z`
/// in the strip of branch `k`.
///
/// `W_0(0) = 0` and `W_k(0) = -∞` for `k != 0`; the branch point `-1/e`
/// gives `-1`.
#[must_use]
pub fn lambert_w_complex(
    z: Complex64,
    k: i32,
) -> Complex64 {
    const INV_E: f64 = 0.367_879_441_171_442_33;
    const INV_E_RESIDUAL: f64 = -1.242_875_367_278_836_3e-17;
    if is_nan(z) {
        return NAN;
    }
    if z.re == 0.0 && z.im == 0.0 {
        return if k == 0 {
            ZERO
        } else {
            Complex64::new(f64::NEG_INFINITY, 0.0)
        };
    }
    if is_inf(z) {
        return INF;
    }
    let branch_offset = Complex64::new(z.re + INV_E, z.im);
    let coalescing_branch = k == 0 || (k == -1 && z.im >= 0.0) || (k == 1 && z.im < 0.0);
    if coalescing_branch && branch_offset == ZERO {
        return Complex64::new(-1.0, 0.0);
    }
    let phase_shift = Complex64::new(0.0, std::f64::consts::TAU * f64::from(k));
    let scale = z.re.abs().max(z.im.abs());
    if scale >= 1e100 || (k != 0 && scale < 1e-100) {
        return lambert_w_logarithmic(scaled_log(z) + phase_shift);
    }
    let e = std::f64::consts::E;
    let mut w = if coalescing_branch && branch_offset.norm() < 0.3 {
        // Keep the low bits of exp(-1) when the argument is one ulp from the branch point.
        let mut p = (2.0 * e * (branch_offset + INV_E_RESIDUAL)).sqrt();
        if k != 0 {
            p = -p;
        }
        let correction = -43.0 / 540.0 + p * (769.0 / 17_280.0);
        let correction = -1.0 / 3.0 + p * (11.0 / 72.0 + p * correction);
        let expansion = -1.0 + p * (1.0 + p * correction);
        if p.norm() < 0.004 {
            // The omitted term is at most about one ulp in W. Iterating
            // would amplify rounding in the residual near the double root.
            return expansion;
        }
        expansion
    } else if k == 0 && z.norm() < 0.3 {
        z - z * z + 1.5 * z * z * z
    } else if k == 0 && z.re >= -INV_E && z.norm() < 2.5 {
        // ln(1 + z) is a safe start away from the branch point.
        (1.0 + z).ln()
    } else {
        let l = z.ln() + phase_shift;
        l - l.ln()
    };
    for _ in 0..100 {
        let ew = w.exp();
        let f = w * ew - z;
        let wp1 = w + 1.0;
        if wp1.norm() == 0.0 {
            break;
        }
        let step = f / (ew * wp1 - (w + 2.0) * f / (2.0 * wp1));
        if !step.re.is_finite() || !step.im.is_finite() {
            break;
        }
        w -= step;
        if step.norm() <= 1e-15 * (1.0 + w.norm()) {
            break;
        }
    }
    w
}

/// `exp(iπx/2)` with exact values at whole quarter-turns.
fn unit_pi_half(x: f64) -> Complex64 {
    let x = x % 4.0;
    let k = x.round();
    let unit = unit_complex(std::f64::consts::FRAC_PI_2 * (x - k));
    match (k as i32).rem_euclid(4) {
        | 0 => unit,
        | 1 => Complex64::new(-unit.im, unit.re),
        | 2 => -unit,
        | _ => Complex64::new(unit.im, -unit.re),
    }
}

/// Half the sum of the two Hankel expansions (DLMF 10.17.5), with the
/// additional phase `order_rotation * πν/2` applied before exponential scaling.
fn bessel_j_asymptotic(
    nu: f64,
    mut z: Complex64,
    mut order_rotation: f64,
) -> Complex64 {
    if z.re < 0.0 {
        // Continue from the corresponding lip of the cut into the right half plane.
        let side = if z.im == 0.0 {
            1.0
        } else {
            z.im.signum()
        };
        order_rotation += 2.0 * side;
        z = -z;
    }
    let inverse = scaled_reciprocal(z);
    let order_square = 0.5 * nu * nu;
    let one = Complex64::new(1.0, 0.0);
    let (mut p, mut q) = (one, one);
    let (mut term_p, mut term_q) = (one, one);
    let mut last = f64::INFINITY;
    for k in 1..80 {
        let kf = f64::from(k);
        let factor = ((order_square - 0.5 * (kf - 0.5).powi(2)) / kf) * inverse;
        term_p *= factor * I;
        term_q *= factor * -I;
        let mag = term_p.norm().max(term_q.norm());
        if mag > last {
            break;
        }
        last = mag;
        p += term_p;
        q += term_q;
        if mag < 1e-17 {
            break;
        }
    }
    let prefactor = (inverse * (0.5 / PI)).sqrt();
    let order = nu % 4.0;
    if z.re == 0.0 {
        // Jν(±iy) = exp(±iπν/2) Iν(y). Keep the terms on this ray even
        // when one component of the Hankel sum overflows.
        let (growing, decaying) = if z.im > 0.0 {
            (q.re, p.re)
        } else {
            (p.re, q.re)
        };
        let common = unit_pi_half((order_rotation + z.im.signum()) * order) * prefactor.norm();
        return scale_real_exponential(common * growing, z.im.abs())
            + scale_real_exponential(
                common * (-unit_pi_half(2.0 * order).im * decaying),
                -z.im.abs(),
            );
    }
    let phase = unit_complex(z.re) * unit_complex(-std::f64::consts::FRAC_PI_4);
    let h1 = prefactor * phase * unit_pi_half((order_rotation - 1.0) * order) * p;
    let h2 = prefactor * phase.conj() * unit_pi_half((order_rotation + 1.0) * order) * q;
    scale_real_exponential(h1, -z.im) + scale_real_exponential(h2, z.im)
}

/// Ascending series `(z/2)^ν Σ (±(z/2)²)^k / (k! Γ(k + ν + 1))` for orders
/// `ν` that are not negative integers.
fn bessel_series(
    nu: f64,
    z: Complex64,
    alternating: bool,
) -> Complex64 {
    let half = z * 0.5;
    let q = if alternating {
        -(half * half)
    } else {
        half * half
    };
    let g = crate::kernels::special::gamma_numerical(nu + 1.0);
    let mut term = Complex64::new(1.0 / g, 0.0);
    let mut sum = term;
    let mut k = 1.0_f64;
    for _ in 0..500 {
        term *= q / (k * (k + nu));
        sum += term;
        if term.norm() <= 1e-17 * sum.norm().max(f64::MIN_POSITIVE) {
            break;
        }
        k += 1.0;
    }
    sum * half.powf(nu)
}

fn integer_order(nu: f64) -> Option<i64> {
    let r = nu.round();
    ((nu - r).abs() < 1e-14 && r.abs() < 1e9).then_some(r as i64)
}

const fn bessel_invalid(
    nu: f64,
    z: Complex64,
) -> bool {
    nu.is_nan() || is_nan(z) || is_inf(z)
}

fn bessel_at_origin(nu: f64) -> Complex64 {
    if nu == 0.0 {
        one()
    } else if nu > 0.0 {
        ZERO
    } else {
        INF
    }
}

const fn one() -> Complex64 {
    Complex64::new(1.0, 0.0)
}

fn use_asymptotic(
    nu: f64,
    z: Complex64,
) -> bool {
    let magnitude = z.norm();
    magnitude >= 14.0 && magnitude > 0.5 * nu * nu
}

/// The Bessel function of the first kind `J_ν(z)` for real order `ν` and
/// complex `z` (principal branch of `z^ν`, cut along the negative real axis).
#[must_use]
pub fn bessel_j_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if bessel_invalid(nu, z) {
        return NAN;
    }
    if let Some(n) = integer_order(nu)
        && n < 0
    {
        // J_{-n} = (-1)^n J_n
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_j_complex(nu.abs(), z);
    }
    if z.re == 0.0 && z.im == 0.0 {
        return bessel_at_origin(nu);
    }
    if use_asymptotic(nu, z) {
        return bessel_j_asymptotic(nu, z, 0.0);
    }
    bessel_series(nu, z, true)
}

/// The modified Bessel function of the first kind `I_ν(z)` for real order
/// `ν` and complex `z` (principal branch).
#[must_use]
pub fn bessel_i_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if bessel_invalid(nu, z) {
        return NAN;
    }
    if let Some(n) = integer_order(nu)
        && n < 0
    {
        // I_{-n} = I_n
        return bessel_i_complex(nu.abs(), z);
    }
    if z.re == 0.0 && z.im == 0.0 {
        return bessel_at_origin(nu);
    }
    if use_asymptotic(nu, z) {
        // DLMF 10.27.6: I_ν(z) = e^{∓νπi/2} J_ν(z e^{±πi/2}),
        // upper sign for -π <= arg z <= π/2.
        let direction = if z.re < 0.0 && z.im.is_sign_positive() {
            -1.0
        } else {
            1.0
        };
        return bessel_j_asymptotic(nu, (direction * I) * z, -direction);
    }
    bessel_series(nu, z, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[track_caller]
    fn assert_close(
        value: Complex64,
        expected: Complex64,
    ) {
        for (actual, target) in [(value.re, expected.re), (value.im, expected.im)] {
            assert!(
                (actual - target).abs() <= 8e-15 * target.abs(),
                "got {value}, expected {expected}"
            );
        }
    }

    #[test]
    fn log_gamma_reflection_preserves_the_unwrapped_phase() {
        for (z, want) in [
            (
                Complex64::new(-1e8, 1.2),
                Complex64::new(-1_742_068_086.456_018, -314_159_244.824_958_75),
            ),
            (
                Complex64::new(-4.25, 0.0),
                Complex64::new(-2.070_072_434_257_324, -5.0 * PI),
            ),
        ] {
            assert_close(ln_gamma_complex(z), want);
            assert_close(ln_gamma_complex(z.conj()), want.conj());
        }
        let large = ln_gamma_complex(Complex64::new(-1.7e308, 7.5e305));
        assert_eq!(large.re, f64::NEG_INFINITY);
        assert!((large.im / -1.775_621_007_399_342_4e306 - 1.0).abs() < 3e-14);
    }

    #[test]
    fn gamma_reflection_preserves_subnormal_and_mixed_overflow_components() {
        let small = gamma_complex(Complex64::new(-172.2, 0.25));
        assert!(
            small
                .re
                .to_bits()
                .abs_diff(2.185_107_940_575_647e-312_f64.to_bits())
                <= 4
        );
        assert!(
            small
                .im
                .to_bits()
                .abs_diff((-4.507_904_050_014_684e-312_f64).to_bits())
                <= 4
        );
        let near_cut = gamma_complex(Complex64::new(-0.3, 1e-200));
        assert!((near_cut.re / -4.326_851_108_825_193 - 1.0).abs() < 3e-14);
        assert!((near_cut.im / -9.143_976_763_306_549e-200 - 1.0).abs() < 3e-14);
        let mixed = gamma_complex(Complex64::new(172.0, 1.2));
        assert_eq!(mixed.re, f64::INFINITY);
        assert!((mixed.im / -1.352_655_263_960_584_5e308 - 1.0).abs() < 1e-13);
        let huge = Complex64::new(1.7e308, 1e-200);
        assert_eq!(
            gamma_complex(huge),
            Complex64::new(f64::INFINITY, f64::INFINITY)
        );
        assert_eq!(
            gamma_complex(huge.conj()),
            Complex64::new(f64::INFINITY, f64::NEG_INFINITY)
        );
    }

    #[test]
    fn bessel_hankel_scaling_preserves_finite_components() {
        let value = Complex64::new(1.335_146_345_194_718_2e307, -2.076_157_839_179_489_5e307);
        assert_close(bessel_j_complex(0.0, Complex64::new(1.0, 712.0)), value);
        assert_close(
            bessel_i_complex(0.0, Complex64::new(712.0, 1.0)),
            value.conj(),
        );

        // The norm overflows here, while both components fit after rotation.
        let z = Complex64::new(std::f64::consts::FRAC_PI_4, 714.25);
        let value = Complex64::new(1.654_315_352_658_447_2e308, -1.652_496_610_211_949_7e308);
        assert_close(bessel_j_complex(0.0, z), value);
        assert_close(
            bessel_i_complex(0.0, Complex64::new(z.im, z.re)),
            value.conj(),
        );
        assert_eq!(bessel_i_complex(0.0, Complex64::new(z.im, 0.0)), INF);

        assert_close(
            bessel_j_complex(0.7, Complex64::new(1e308, 1.0)),
            Complex64::new(8.700_095_277_546_95e-155, -6.634_780_775_020_227e-155),
        );
    }

    #[test]
    fn bessel_i_keeps_the_two_lips_of_the_negative_real_cut() {
        let z = Complex64::new(-712.0, 0.0);
        let value = Complex64::new(0.0, 2.467_977_432_400_639_6e307);
        assert_close(bessel_i_complex(0.5, z), value);
        assert_close(bessel_i_complex(0.5, z.conj()), value.conj());
    }

    #[test]
    fn digamma_near_poles_and_at_large_complex_arguments() {
        // High-precision references for the exact binary64 inputs.
        for (z, want) in [
            (
                Complex64::new(-1.0, 300.0),
                Complex64::new(5.703_794_511_570_811, 1.575_796_289_758_271_1),
            ),
            (
                Complex64::new(-1.6e308, 1.7e308),
                Complex64::new(710.044_016_446_503_6, 2.325_900_730_273_57),
            ),
            (
                Complex64::new(-128.0, 1e-200),
                Complex64::new(4.855_931_427_687_639_5, 1e200),
            ),
            (
                Complex64::new(0.0, 1e-300),
                Complex64::new(-std::f64::consts::EULER_GAMMA, 1e300),
            ),
            (
                Complex64::new(-1e12 - 0.125, 1e-4),
                Complex64::new(35.215_491_596_208_65, 0.006_739_383_341_793_43),
            ),
        ] {
            assert_close(digamma_complex(z), want);
            assert_close(digamma_complex(z.conj()), want.conj());
        }
    }

    #[test]
    fn lambert_w_logarithmic_iteration_at_extreme_magnitudes() {
        for (z, branch, want) in [
            (
                Complex64::new(1e308, 1e307),
                0,
                Complex64::new(702.646_330_118_942_5, 0.099_527_006_542_680_08),
            ),
            (
                Complex64::new(1.7e308, 1.7e308),
                2,
                Complex64::new(703.517_138_694_697_5, 13.332_819_383_684_028),
            ),
            (
                Complex64::new(-1.7e308, 0.0),
                -1,
                Complex64::new(703.171_226_513_655_5, -3.137_131_278_683_244),
            ),
            (
                Complex64::new(1e-200, 1e-200),
                1,
                Complex64::new(-466.315_342_725_542_4, 3.935_430_034_499_216_7),
            ),
        ] {
            assert_close(lambert_w_complex(z, branch), want);
        }
    }

    #[test]
    fn lambert_w_keeps_the_coalescing_branches_and_cut_lips() {
        let branch_point = -1.0 / std::f64::consts::E;
        let z = Complex64::new(branch_point, 1e-20);
        assert_close(
            lambert_w_complex(z, 0),
            Complex64::new(-0.999_999_999_996_693_1, 8.220_080_379_985_646e-9),
        );
        let secondary = Complex64::new(-1.000_000_000_003_306_9, -8.220_080_380_021_89e-9);
        assert_close(lambert_w_complex(z, -1), secondary);
        assert_close(lambert_w_complex(z.conj(), 1), secondary.conj());
        let cut = Complex64::new(-0.7, 1e-20);
        let cut_value = Complex64::new(-0.564_874_052_922_317_8, 1.094_260_913_312_513_8);
        assert_close(lambert_w_complex(cut, 0), cut_value);
        assert_close(lambert_w_complex(cut.conj(), 0), cut_value.conj());
        assert_close(
            lambert_w_complex(Complex64::new(branch_point.next_up(), 0.0), -1),
            Complex64::new(-1.000_000_015_304_254_3, 0.0),
        );
        assert_eq!(
            lambert_w_complex(Complex64::new(branch_point, 0.0), -1),
            Complex64::new(-1.0, 0.0)
        );
    }
}
