//! Complex gamma, reciprocal gamma, polygamma, beta, incomplete gamma/beta,
//! Pochhammer, binomial, and Barnes G functions.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

use crate::kernels::hypergeometric::hyp2f1_cplx;
use crate::kernels::special_complex::digamma_complex;
use crate::kernels::special_complex::gamma_complex;
use crate::kernels::special_complex::ln_gamma_complex;
use crate::kernels::special_ext::polylog_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

fn is_nonpos_int(z: Complex64) -> bool {
    z.im == 0.0 && z.re <= 0.0 && (z.re - z.re.round()).abs() == 0.0
}

/// Exp-safe `sin(pi * z) / pi` in log-space when `|Im(z)|` is large,
/// returning `ln(sin(pi*z) / pi)`.
fn ln_sin_pi_over_pi(z: Complex64) -> Complex64 {
    let y = z.im;
    if y.abs() <= 20.0 {
        return ((PI * z).sin() / PI).ln();
    }
    // sin(pi*z) = (e^{i*pi*z} - e^{-i*pi*z}) / (2i)
    if y > 0.0 {
        let two_i_pz = Complex64::new(-2.0 * PI * y, 2.0 * PI * z.re);
        let base = Complex64::new(PI * y - (2.0 * PI).ln(), -PI * z.re + 0.5 * PI);
        base + (ONE - two_i_pz.exp()).ln()
    } else {
        let two_i_pz = Complex64::new(2.0 * PI * y, -2.0 * PI * z.re);
        let base = Complex64::new(-PI * y - (2.0 * PI).ln(), PI * z.re - 0.5 * PI);
        base + (ONE - two_i_pz.exp()).ln()
    }
}

/// Exp-safe `cot(pi * z)`.
fn cot_pi(z: Complex64) -> Complex64 {
    let py = PI * z.im;
    if py > 20.0 {
        return Complex64::new(0.0, -1.0);
    }
    if py < -20.0 {
        return Complex64::new(0.0, 1.0);
    }
    let pz = PI * z;
    pz.cos() / pz.sin()
}

/// Reciprocal gamma function `1 / Gamma(z)`, which is entire on `C` with
/// zeros at the non-positive integers.
#[must_use]
pub fn rgamma_complex(z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if is_nonpos_int(z) {
        return ZERO;
    }
    if z.re < 0.5 {
        if z.im.abs() > 20.0 || z.re < -150.0 {
            return (ln_sin_pi_over_pi(z) + ln_gamma_complex(ONE - z)).exp();
        }
        let s = (Complex64::new(PI, 0.0) * z).sin();
        return s * gamma_complex(ONE - z) / PI;
    }
    let g = gamma_complex(z);
    if g.re.is_infinite() || g.im.is_infinite() {
        return (-ln_gamma_complex(z)).exp();
    }
    ONE / g
}

/// Computes `d^n/dz^n (pi * cot(pi * z))` via polynomial recurrence in `t = cot(pi * z)`:
/// `P_0(t) = t`, `P_{k+1}(t) = -pi * (1 + t^2) * P_k'(t)`.
fn d_n_pi_cot_pi(
    n: u32,
    z: Complex64,
) -> Complex64 {
    let t = cot_pi(z);
    // Represent P_k(t) as coefficients c[0] + c[1]*t + ... + c[k+1]*t^{k+1}.
    let mut coeffs = vec![ZERO, ONE]; // t
    for _ in 0..n {
        let mut next = vec![ZERO; coeffs.len() + 1];
        for (j, &cj) in coeffs.iter().enumerate().skip(1) {
            let deriv_coeff = cj * (j as f64);
            next[j - 1] += deriv_coeff;
            next[j + 1] += deriv_coeff;
        }
        for c in &mut next {
            *c *= -PI;
        }
        coeffs = next;
    }
    let mut val = ZERO;
    for &c in coeffs.iter().rev() {
        val = val * t + c;
    }
    PI * val
}

/// Complex polygamma function `psi^(n)(z) = d^n/dz^n digamma(z)`.
#[must_use]
pub fn polygamma_complex(
    n: u32,
    mut z: Complex64,
) -> Complex64 {
    // Asymptotic expansion:
    // psi^(n)(z) ~ (-1)^(n+1) * [ (n-1)!/z^n + n!/(2 z^(n+1)) + sum_{k>=1} B_{2k}/(2k)! * (2k+n-1)! / z^(2k+n) ]
    // B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, B_8 = -1/30, B_10 = 5/66, B_12 = -691/2730.
    const BERNOULLI_EVEN: [f64; 6] = [
        1.0 / 6.0,
        -1.0 / 30.0,
        1.0 / 42.0,
        -1.0 / 30.0,
        5.0 / 66.0,
        -691.0 / 2730.0,
    ];

    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if n == 0 {
        return digamma_complex(z);
    }
    if is_nonpos_int(z) {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    // General reflection for Re(z) < 0:
    // (-1)^n psi^(n)(1 - z) - psi^(n)(z) = d^n/dz^n (pi * cot(pi * z))
    if z.re < 0.0 {
        let sign_n = if n.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        let refl = d_n_pi_cot_pi(n, z);
        return sign_n * polygamma_complex(n, ONE - z) - refl;
    }

    let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
    let mut fact = 1.0_f64;
    for k in 1..=n {
        fact *= f64::from(k);
    }

    // Shift z right until Re(z) >= 10.
    let mut sum = ZERO;
    while z.re < 10.0 {
        sum += ONE / z.powi(i32::try_from(n + 1).unwrap_or(i32::MAX));
        z += ONE;
    }

    let inv_z = ONE / z;
    let inv_z2 = inv_z * inv_z;
    let nf = f64::from(n);
    // (n-1)! / z^n + n! / (2 * z^(n+1)) + sum_{k>=1} B_{2k}/(2k)! * (2k+n-1)! / z^(2k+n)
    let mut base_term = fact * inv_z.powi(i32::try_from(n).unwrap_or(i32::MAX));
    let mut asym = base_term / nf + 0.5 * base_term * inv_z;

    // For k = 1..=6: term_k = B_{2k} * [(2k+n-1)! / ((2k)! * n!)] * (n! / z^(n + 2k))
    let mut ratio = 1.0_f64;
    for (idx, &b2k) in BERNOULLI_EVEN.iter().enumerate() {
        let k = (idx + 1) as f64;
        let two_k = 2.0 * k;
        if idx == 0 {
            ratio = f64::midpoint(nf, 1.0);
        } else {
            ratio *= (two_k + nf - 2.0) * (two_k + nf - 1.0) / ((two_k - 1.0) * two_k);
        }
        base_term *= inv_z2;
        asym += b2k * ratio * base_term;
    }

    sign * (fact * sum + asym)
}

/// Complex natural log-beta function `ln B(a, b) = ln Gamma(a) + ln Gamma(b) - ln Gamma(a + b)`.
#[must_use]
pub fn ln_beta_complex(
    a: Complex64,
    b: Complex64,
) -> Complex64 {
    ln_gamma_complex(a) + ln_gamma_complex(b) - ln_gamma_complex(a + b)
}

/// Complex beta function `B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b)`.
#[must_use]
pub fn beta_complex(
    a: Complex64,
    b: Complex64,
) -> Complex64 {
    if is_nonpos_int(a) || is_nonpos_int(b) {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    if is_nonpos_int(a + b) {
        return ZERO;
    }
    ln_beta_complex(a, b).exp()
}

/// Complex Pochhammer symbol `(z)_w = Gamma(z + w) / Gamma(z)`.
#[must_use]
pub fn pochhammer_complex(
    z: Complex64,
    w: Complex64,
) -> Complex64 {
    if w == ZERO {
        return ONE;
    }
    if w.im == 0.0 && w.re > 0.0 && w.re.fract() == 0.0 && w.re <= 20.0 {
        let n = w.re as usize;
        let mut prod = ONE;
        for k in 0..n {
            prod *= z + (k as f64);
        }
        return prod;
    }
    if is_nonpos_int(z) {
        return if is_nonpos_int(z + w) {
            // Ratio of residues (-1)^(n - m) * n! / m!
            let n = (-z.re).round();
            let m = (-(z + w).re).round();
            let sign = if ((n - m).abs() as i64) % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            Complex64::new(
                sign * (ln_gamma_complex(Complex64::new(n + 1.0, 0.0))
                    - ln_gamma_complex(Complex64::new(m + 1.0, 0.0)))
                .re
                .exp(),
                0.0,
            )
        } else {
            ZERO
        };
    }
    if is_nonpos_int(z + w) {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    (ln_gamma_complex(z + w) - ln_gamma_complex(z)).exp()
}

/// Complex binomial coefficient `C(n, k) = Gamma(n + 1) / (Gamma(k + 1) * Gamma(n - k + 1))`.
#[must_use]
pub fn binomial_complex(
    n: Complex64,
    k: Complex64,
) -> Complex64 {
    if k.im == 0.0 && k.re.fract() == 0.0 {
        if k.re < 0.0 {
            return ZERO;
        }
        if k.re <= 30.0 {
            let m = k.re as usize;
            let mut prod = ONE;
            for j in 0..m {
                prod *= (n - (j as f64)) / ((j + 1) as f64);
            }
            return prod;
        }
    }
    if is_nonpos_int(k + ONE) || is_nonpos_int(n - k + ONE) {
        return ZERO;
    }
    if is_nonpos_int(n + ONE) {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    (ln_gamma_complex(n + ONE) - ln_gamma_complex(k + ONE) - ln_gamma_complex(n - k + ONE)).exp()
}

/// Lower incomplete gamma function `gamma(a, z) = int_0^z t^(a-1) e^(-t) dt`.
#[must_use]
pub fn lower_incomplete_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return ZERO;
    }
    if z.norm() <= a.norm() + 2.0 || (z.re < 0.0 && z.norm() < 22.0) {
        // Power series: z^a * exp(-z) * sum_{k>=0} z^k / (a * (a+1) * ... * (a+k))
        let mut term = ONE / a;
        let mut sum = term;
        for k in 1..400 {
            term *= z / (a + f64::from(k));
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        (a * z.ln() - z + sum.ln()).exp()
    } else {
        gamma_complex(a) - upper_incomplete_gamma_complex(a, z)
    }
}

/// Upper incomplete gamma function `Gamma(a, z) = int_z^oo t^(a-1) e^(-t) dt`.
#[must_use]
pub fn upper_incomplete_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return gamma_complex(a);
    }
    if z.norm() <= a.norm() + 2.0 || (z.re < 0.0 && z.norm() < 22.0) {
        return gamma_complex(a) - lower_incomplete_gamma_complex(a, z);
    }
    // Asymptotic expansion when |z| >= 22 near (-inf, 0]: avoids divergent Legendre CF on cut.
    if z.re < 0.0 && z.im.abs() < 4.0 {
        let mut sum = ONE;
        let mut term = ONE;
        let mut best_sum = sum;
        let mut min_term = f64::INFINITY;
        for k in 1..=60 {
            term *= (a - f64::from(k)) / z;
            let tn = term.norm();
            if tn > min_term {
                break;
            }
            min_term = tn;
            sum += term;
            best_sum = sum;
            if tn <= 1e-16 * sum.norm() {
                break;
            }
        }
        return ((a - ONE) * z.ln() - z + best_sum.ln()).exp();
    }
    // Legendre continued fraction via modified Lentz's method:
    let tiny = 1e-30_f64;
    let mut b = z + ONE - a;
    let mut c = Complex64::new(1.0 / tiny, 0.0);
    let mut d = if b.norm() < tiny {
        Complex64::new(1.0 / tiny, 0.0)
    } else {
        ONE / b
    };
    let mut h = d;
    for k in 1..400 {
        let kf = f64::from(k);
        let an = -kf * (kf - a);
        b += 2.0;
        d = an * d + b;
        if d.norm() < tiny {
            d = Complex64::new(tiny, 0.0);
        }
        c = b + an / c;
        if c.norm() < tiny {
            c = Complex64::new(tiny, 0.0);
        }
        d = ONE / d;
        let del = d * c;
        h *= del;
        if (del - ONE).norm() < 1e-15 {
            break;
        }
    }
    (a * z.ln() - z + h.ln()).exp()
}

/// Regularized lower incomplete gamma `P(a, z) = gamma(a, z) / Gamma(a)`.
#[must_use]
pub fn reg_lower_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
    if is_nonpos_int(a) {
        return ONE;
    }
    lower_incomplete_gamma_complex(a, z) * rgamma_complex(a)
}

/// Regularized upper incomplete gamma `Q(a, z) = Gamma(a, z) / Gamma(a)`.
#[must_use]
pub fn reg_upper_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
    if is_nonpos_int(a) {
        return ZERO;
    }
    upper_incomplete_gamma_complex(a, z) * rgamma_complex(a)
}

/// Incomplete beta function `B_z(a, b) = (z^a / a) * 2F1(a, 1 - b; a + 1; z)`.
#[must_use]
pub fn incomplete_beta_complex(
    a: Complex64,
    b: Complex64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return ZERO;
    }
    if z == ONE {
        return beta_complex(a, b);
    }
    (z.powc(a) / a) * hyp2f1_cplx(a, ONE - b, a + ONE, z)
}

/// Regularized incomplete beta `I_z(a, b) = B_z(a, b) / B(a, b)`.
#[must_use]
pub fn reg_beta_complex(
    a: Complex64,
    b: Complex64,
    z: Complex64,
) -> Complex64 {
    incomplete_beta_complex(a, b, z) / beta_complex(a, b)
}

/// Analytic logarithm of the Barnes G-function, real on the positive axis.
/// Signed imaginary zero selects the side of the non-positive real cut.
/// Non-positive integers give negative infinity.
#[must_use]
pub fn ln_barnes_g_complex(mut z: Complex64) -> Complex64 {
    const ZETA_PRIME_MINUS_1: f64 = -0.165_421_143_700_450_92;
    const LN_2PI: f64 = 1.837_877_066_409_345_3;
    const SHIFT_TARGET: f64 = 8.0;
    const MAX_SHIFTS: usize = 16;
    const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);

    // DLMF 5.17.5 reduced by Stirling's expansion of w*log(Gamma(w + 1)).
    // The coefficient of w^(-2k) is B_{2k+2} / (4k(k+1)).
    const COEFFS: [f64; 8] = [
        -1.0 / 240.0,
        1.0 / 1008.0,
        -1.0 / 1440.0,
        1.0 / 1056.0,
        -691.0 / 327_600.0,
        1.0 / 144.0,
        -3617.0 / 114_240.0,
        43_867.0 / 229_824.0,
    ];

    if !z.re.is_finite() || !z.im.is_finite() {
        return NAN;
    }
    if is_nonpos_int(z) {
        return Complex64::new(f64::NEG_INFINITY, 0.0);
    }
    if z == ONE || z == Complex64::new(2.0, 0.0) {
        return ZERO;
    }
    let ln_factor = if z.re <= -SHIFT_TARGET {
        // Kinkelin reflection for N = 1 - z in the upper half-plane.
        // Subtracting the left asymptotic series leaves only
        // N*log(1 - q) - i*Li_2(q)/(2*pi), where q = exp(2*pi*i*N).
        let height = z.im.abs();
        let log_radius = -TAU * height;
        let radius = log_radius.exp();
        // Reduce the real part before multiplying by 2*pi, also when
        // 1 - Re(z) would round away a fractional part.
        let phase = -TAU * (z.re - z.re.round());
        let (sin, cos) = phase.sin_cos();
        let q = radius * Complex64::new(cos, sin);
        let sin_half = (0.5 * phase).sin();
        let one_minus_q = Complex64::new(
            -log_radius.exp_m1() + 2.0 * radius * sin_half * sin_half,
            -q.im,
        );
        let dilog = polylog_complex(Complex64::new(2.0, 0.0), q);
        let correction = Complex64::new(1.0 - z.re, height) * one_minus_q.ln()
            + Complex64::new(dilog.im / TAU, -dilog.re / TAU);
        if z.im.is_sign_negative() {
            correction
        } else {
            correction.conj()
        }
    } else {
        let mut recurrence = ZERO;
        for _ in 0..MAX_SHIFTS {
            if z.re >= SHIFT_TARGET {
                break;
            }
            let next = z.re + 1.0;
            if next <= z.re {
                return NAN;
            }
            recurrence -= ln_gamma_complex(z);
            // Updating the real component preserves the signed cut side.
            z.re = next;
        }
        recurrence
    };
    let w = z - 1.0;
    let ln_w = w.ln();
    let w2 = w * w;
    let mut ln_g =
        0.5 * w2 * ln_w - 0.75 * w2 + 0.5 * w * LN_2PI - (1.0 / 12.0) * ln_w + ZETA_PRIME_MINUS_1;
    let inv_w2 = ONE / w2;
    let mut term_pow = inv_w2;
    for &c in &COEFFS {
        ln_g += c * term_pow;
        term_pow *= inv_w2;
    }
    ln_g + ln_factor
}

/// Complex Barnes G-function `G(z)` satisfying `G(z + 1) = Gamma(z) * G(z)`, `G(1) = 1`.
#[must_use]
pub fn barnes_g_complex(z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if is_nonpos_int(z) {
        return ZERO;
    }
    if z == ONE || z == Complex64::new(2.0, 0.0) {
        return ONE;
    }
    ln_barnes_g_complex(z).exp()
}
