//! # Numerical Special Functions
//!
//! This module provides numerical implementations of mathematical special functions.
//! It includes gamma functions, beta functions, error functions, Bessel functions,
//! orthogonal polynomials, and other commonly used special functions.

use std::f64::consts::FRAC_2_PI;

use statrs::function::beta::beta;
use statrs::function::beta::ln_beta;
use statrs::function::erf::erf;
use statrs::function::erf::erfc;
use statrs::function::gamma::digamma;
use statrs::function::gamma::gamma;
use statrs::function::gamma::ln_gamma;

// ============================================================================
// Gamma Functions
// ============================================================================

/// Computes the gamma function, Γ(x).
/// Γ(n) = (n-1)! for positive integers.
#[must_use]
pub fn gamma_numerical(x: f64) -> f64 {
    gamma(x)
}

/// Computes the natural logarithm of the gamma function, ln(Γ(x)).
#[must_use]
pub fn ln_gamma_numerical(x: f64) -> f64 {
    ln_gamma(x)
}

/// Computes the digamma function, ψ(x) = d/dx ln(Γ(x)) = Γ'(x)/Γ(x).
#[must_use]
pub fn digamma_numerical(x: f64) -> f64 {
    digamma(x)
}

/// Computes the lower incomplete gamma function, γ(s, x) = ∫₀ˣ t^(s-1) e^(-t) dt.
/// Uses series expansion for computation.
#[must_use]
pub fn lower_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {
    if x < 0.0 || s <= 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return 0.0;
    }

    // Use series expansion: γ(s,x) = x^s * e^(-x) * Σ(x^n / Γ(s+n+1))
    // Or regularized: P(s,x) = γ(s,x) / Γ(s)
    // Using statrs regularized gamma if available, otherwise series
    let mut sum = 0.0;

    let mut term = 1.0 / s;

    sum += term;

    for n in 1..200 {
        term *= x / (s + f64::from(n));

        sum += term;

        if term.abs() < 1e-15 * sum.abs() {
            break;
        }
    }

    x.powf(s) * (-x).exp() * sum
}

/// Computes the upper incomplete gamma function, Γ(s, x) = ∫ₓ^∞ t^(s-1) e^(-t) dt.
#[must_use]
pub fn upper_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {
    gamma(s) - lower_incomplete_gamma(s, x)
}

/// Computes the regularized lower incomplete gamma function, P(s, x) = γ(s, x) / Γ(s).
#[must_use]
pub fn regularized_lower_gamma(
    s: f64,
    x: f64,
) -> f64 {
    lower_incomplete_gamma(s, x) / gamma(s)
}

/// Computes the regularized upper incomplete gamma function, Q(s, x) = Γ(s, x) / Γ(s).
#[must_use]
pub fn regularized_upper_gamma(
    s: f64,
    x: f64,
) -> f64 {
    1.0 - regularized_lower_gamma(s, x)
}

// ============================================================================
// Beta Functions
// ============================================================================

/// Computes the beta function, B(a, b) = Γ(a)Γ(b)/Γ(a+b).
#[must_use]
pub fn beta_numerical(
    a: f64,
    b: f64,
) -> f64 {
    beta(a, b)
}

/// Computes the natural logarithm of the beta function, ln(B(a, b)).
#[must_use]
pub fn ln_beta_numerical(
    a: f64,
    b: f64,
) -> f64 {
    ln_beta(a, b)
}

/// Computes the incomplete beta function, B(x; a, b) = ∫₀ˣ t^(a-1) (1-t)^(b-1) dt.
#[must_use]
pub fn incomplete_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {
    if !(0.0..=1.0).contains(&x) || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return 0.0;
    }

    if (x - 1.0).abs() < 1e-15 {
        return beta(a, b);
    }

    // Use series expansion
    regularized_beta(x, a, b) * beta(a, b)
}

/// Computes the regularized incomplete beta function, `I_x(a`, b) = B(x; a, b) / B(a, b).
/// Uses continued fraction expansion for better convergence.
#[must_use]
pub fn regularized_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {
    if !(0.0..=1.0).contains(&x) || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return 0.0;
    }

    if (x - 1.0).abs() < 1e-15 {
        return 1.0;
    }

    // Use symmetry for better convergence
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularized_beta(1.0 - x, b, a);
    }

    // Continued fraction expansion
    let bt = if x == 0.0 || (x - 1.0).abs() < 1e-15 {
        0.0
    } else {
        b.mul_add(
            (1.0 - x).ln(),
            a.mul_add(x.ln(), ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)),
        )
        .exp()
    };

    let mut am = 1.0;

    let mut bm = 1.0;

    let mut az = 1.0;

    let qab = a + b;

    let qap = a + 1.0;

    let qam = a - 1.0;

    let mut bz = 1.0 - qab * x / qap;

    for m in 1..200 {
        let em = f64::from(m);

        let tem = em + em;

        let d = em * (b - em) * x / ((qam + tem) * (a + tem));

        let ap = az + d * am;

        let bp = bz + d * bm;

        let d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));

        let app = ap + d * az;

        let bpp = bp + d * bz;

        let aold = az;

        am = ap / bpp;

        bm = bp / bpp;

        az = app / bpp;

        bz = 1.0;

        if (az - aold).abs() < 1e-14 * az.abs() {
            return bt * az / a;
        }
    }

    bt * az / a
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the error function, erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt.
#[must_use]
pub fn erf_numerical(x: f64) -> f64 {
    erf(x)
}

/// Computes the complementary error function, erfc(x) = 1 - erf(x).
#[must_use]
pub fn erfc_numerical(x: f64) -> f64 {
    erfc(x)
}

/// Computes the inverse error function, erf⁻¹(x).
///
/// Starts from Winitzki's closed-form approximation and polishes with Halley
/// iterations on `erfc(y) = 1 - x` (for x ≥ 0.5, where `1 - x` is exact and
/// the tail keeps full relative precision) or on `erf(y) = x` otherwise.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn inverse_erf_numerical(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }

    if x <= -1.0 {
        return if x == -1.0 {
            f64::NEG_INFINITY
        } else {
            f64::NAN
        };
    }

    if x >= 1.0 {
        return if x == 1.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };
    }

    if x.abs() < 1e-15 {
        // erf(y) = 2y/sqrt(pi) to first order.
        return x * std::f64::consts::PI.sqrt() / 2.0;
    }

    let sign = if x < 0.0 { -1.0 } else { 1.0 };

    let x = x.abs();

    // Winitzki (2008): a = 0.147 gives a relative error below 2e-3.
    let a = 0.147;

    let ln1 = ((1.0 - x) * (1.0 + x)).ln();

    let t = 2.0 / (std::f64::consts::PI * a) + ln1 / 2.0;

    #[allow(clippy::suspicious_operation_groupings)] // t^2 - ln1/a is the intended discriminant
    let mut y = ((t * t - ln1 / a).sqrt() - t).sqrt();

    let two_over_sqrt_pi = 2.0 / std::f64::consts::PI.sqrt();

    // Halley step for f(y) with f''/f' = -2y (true for both erf - x and
    // erfc - c because f' = ±(2/√π) e^(-y²)): y -= h / (1 + h y), h = f / f'.
    for _ in 0..8 {
        let deriv = two_over_sqrt_pi * (-y * y).exp();

        let h = if x >= 0.5 {
            (erfc(y) - (1.0 - x)) / -deriv
        } else {
            (erf(y) - x) / deriv
        };

        let step = h / (1.0 + h * y);

        y -= step;

        if step.abs() <= 1e-16 * y.abs() {
            break;
        }
    }

    sign * y
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Euler-Mascheroni constant γ.
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Below this |x| the alternating power series are used for `J` and `Y`
/// (largest term ~ I₀(5) ≈ 27, so cancellation costs under one digit).
const BESSEL_SERIES_MAX: f64 = 5.0;

/// Above this |x| the Hankel asymptotic expansion is used for `J` and `Y`; its
/// smallest term is about e^(-2x) ≈ 2e-22 there, below double precision.
const BESSEL_ASYMPTOTIC_MIN: f64 = 25.0;

/// Below this |x| the (all-positive, hence cancellation-free) power series is
/// used for `I`; above it the asymptotic expansion (smallest term ~ e^(-2x)).
const BESSEL_I_SERIES_MAX: f64 = 30.0;

/// Power series `Σ (-1)^k (x²/4)^k / (k! (k+ν)!)` for ν = 0 or 1 times
/// `(x/2)^ν` (Abramowitz & Stegun 9.1.10): Jν(x) for small |x|.
fn bessel_j_series(
    nu: u32,
    x: f64,
) -> f64 {
    let q = x * x / 4.0;

    let mut term = if nu == 0 { 1.0 } else { x / 2.0 };

    let mut sum = term;

    for k in 1..200 {
        let kf = f64::from(k);

        term *= -q / (kf * (kf + f64::from(nu)));

        sum += term;

        if term.abs() < 1e-18 * sum.abs().max(1e-300) {
            break;
        }
    }

    sum
}

/// Hankel asymptotic expansion (A&S 9.2.5-9.2.10) returning `(P, Q)` for
/// order `nu` at large positive `x`, summed until the terms stop decreasing.
fn hankel_pq(
    nu: u32,
    x: f64,
) -> (f64, f64) {
    let mu = 4.0 * f64::from(nu) * f64::from(nu);

    let mut p = 1.0;

    let mut q = 0.0;

    // a_k(ν) / x^k with a_k = a_{k-1} (μ - (2k-1)²) / (8k)
    let mut term = 1.0;

    let mut last = f64::INFINITY;

    for k in 1..100 {
        let kf = f64::from(k);

        term *= (mu - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf * x);

        if term.abs() > last {
            break;
        }

        last = term.abs();

        // Signs follow (-1)^floor(k/2) for P (even k) and Q (odd k).
        match k % 4 {
            | 0 => p += term,
            | 1 => q += term,
            | 2 => p -= term,
            | _ => q -= term,
        }

        if term.abs() < 1e-17 {
            break;
        }
    }

    (p, q)
}

/// Bessel `Jν` (ν = 0, 1) for large positive `x` from the Hankel expansion.
fn bessel_j_asymptotic(
    nu: u32,
    x: f64,
) -> f64 {
    let (p, q) = hankel_pq(nu, x);

    let chi = x - (f64::from(nu) / 2.0 + 0.25) * std::f64::consts::PI;

    (FRAC_2_PI / x).sqrt() * p.mul_add(chi.cos(), -q * chi.sin())
}

/// Bessel `Yν` (ν = 0, 1) for large positive `x` from the Hankel expansion.
fn bessel_y_asymptotic(
    nu: u32,
    x: f64,
) -> f64 {
    let (p, q) = hankel_pq(nu, x);

    let chi = x - (f64::from(nu) / 2.0 + 0.25) * std::f64::consts::PI;

    (FRAC_2_PI / x).sqrt() * p.mul_add(chi.sin(), q * chi.cos())
}

/// `J₀(x), J₁(x), …, J_m(x)` for 5 ≤ x ≤ 25 by Miller's backward recurrence
/// `J_{n-1} = (2n/x) J_n - J_{n+1}`, normalised with `J₀ + 2 Σ J_{2k} = 1`
/// (A&S 9.1.46). Backward recurrence is stable for n > x, and the start index
/// is far beyond x so the (arbitrary) starting values are damped out.
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn bessel_j_miller(
    x: f64,
    m: usize,
) -> Vec<f64> {
    debug_assert!(m.is_multiple_of(2));

    let start = m + 2 * (x as usize + 30);

    let start = start + start % 2;

    let mut j = vec![0.0; start + 2];

    j[start] = 1e-30;

    for n in (1..=start).rev() {
        j[n - 1] = (2.0 * n as f64 / x).mul_add(j[n], -j[n + 1]);

        if j[n - 1].abs() > 1e250 {
            // Rescale to avoid overflow; only ratios matter before normalising.
            for v in &mut j[(n - 1)..] {
                *v *= 1e-250;
            }
        }
    }

    let norm = j[0] + 2.0 * (1..=start / 2).map(|k| j[2 * k]).sum::<f64>();

    j.truncate(m + 1);

    for v in &mut j {
        *v /= norm;
    }

    j
}

/// Number of Bessel orders needed for the Neumann series for `Y` at `x`.
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
const fn neumann_order(x: f64) -> usize {
    let m = x as usize + 40;

    m + m % 2
}

/// Computes the Bessel function of the first kind, J₀(x).
///
/// Accurate to roughly 1e-15 (absolute): power series for |x| < 5, Miller
/// backward recurrence for 5 ≤ |x| ≤ 25 and the Hankel asymptotic expansion
/// beyond.
#[must_use]
pub fn bessel_j0(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }

    let ax = x.abs();

    if ax < BESSEL_SERIES_MAX {
        bessel_j_series(0, ax)
    } else if ax <= BESSEL_ASYMPTOTIC_MIN {
        bessel_j_miller(ax, 2)[0]
    } else if ax.is_finite() {
        bessel_j_asymptotic(0, ax)
    } else {
        0.0
    }
}

/// Computes the Bessel function of the first kind, J₁(x) (odd in `x`).
///
/// Same method and accuracy as [`bessel_j0`].
#[must_use]
pub fn bessel_j1(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }

    let ax = x.abs();

    let v = if ax < BESSEL_SERIES_MAX {
        bessel_j_series(1, ax)
    } else if ax <= BESSEL_ASYMPTOTIC_MIN {
        bessel_j_miller(ax, 2)[1]
    } else if ax.is_finite() {
        bessel_j_asymptotic(1, ax)
    } else {
        0.0
    };

    if x < 0.0 { -v } else { v }
}

/// Computes the Bessel function of the second kind, Y₀(x), for x ≥ 0
/// (`NaN` for x < 0, `-∞` at 0).
///
/// Series (A&S 9.1.13) for x < 5, the Neumann expansion
/// `Y₀ = (2/π)(ln(x/2) + γ) J₀ - (4/π) Σ (-1)^k J_{2k}/k` (A&S 9.1.88) for
/// 5 ≤ x ≤ 25 and the Hankel expansion beyond; accurate to about 1e-14.
#[must_use]
pub fn bessel_y0(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return f64::NEG_INFINITY;
    }

    if x < BESSEL_SERIES_MAX {
        let q = x * x / 4.0;

        let mut term = 1.0;

        let mut harmonic = 0.0;

        let mut sum = 0.0;

        for k in 1..200 {
            let kf = f64::from(k);

            term *= -q / (kf * kf);

            harmonic += 1.0 / kf;

            // (-1)^(k+1) H_k (x²/4)^k / (k!)²
            let t = -term * harmonic;

            sum += t;

            if t.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
        }

        FRAC_2_PI * (((x / 2.0).ln() + EULER_GAMMA) * bessel_j_series(0, x) + sum)
    } else if x <= BESSEL_ASYMPTOTIC_MIN {
        let m = neumann_order(x);

        let j = bessel_j_miller(x, m);

        let mut sum = 0.0;

        for k in 1..=m / 2 {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };

            sum += sign * j[2 * k] / k as f64;
        }

        FRAC_2_PI * ((x / 2.0).ln() + EULER_GAMMA).mul_add(j[0], -2.0 * sum)
    } else if x.is_finite() {
        bessel_y_asymptotic(0, x)
    } else {
        0.0
    }
}

/// Computes the Bessel function of the second kind, Y₁(x), for x ≥ 0
/// (`NaN` for x < 0, `-∞` at 0).
///
/// Series (A&S 9.1.11) for x < 5; for 5 ≤ x ≤ 25 the derivative of the
/// Neumann expansion of Y₀ (`Y₁ = -Y₀'`, `J₀' = -J₁`,
/// `J_{2k}' = (J_{2k-1} - J_{2k+1})/2`); Hankel expansion beyond.
/// Accurate to about 1e-14.
#[must_use]
pub fn bessel_y1(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return f64::NEG_INFINITY;
    }

    if x < BESSEL_SERIES_MAX {
        let q = x * x / 4.0;

        // (x/2)^(2k+1) / (k! (k+1)!) with alternating sign, and
        // ψ(k+1) + ψ(k+2) = -2γ + H_k + H_{k+1}.
        let mut term = x / 2.0;

        let mut h_k = 0.0;

        let mut h_k1 = 1.0;

        let mut sum = term * (-2.0 * EULER_GAMMA + h_k + h_k1);

        for k in 1..200 {
            let kf = f64::from(k);

            term *= -q / (kf * (kf + 1.0));

            h_k += 1.0 / kf;

            h_k1 += 1.0 / (kf + 1.0);

            let t = term * (-2.0 * EULER_GAMMA + h_k + h_k1);

            sum += t;

            if t.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
        }

        FRAC_2_PI * (((x / 2.0).ln() * bessel_j_series(1, x)) - 1.0 / x)
            - sum / std::f64::consts::PI
    } else if x <= BESSEL_ASYMPTOTIC_MIN {
        let m = neumann_order(x);

        let j = bessel_j_miller(x, m + 2);

        let mut sum = 0.0;

        for k in 1..=m / 2 {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };

            sum += sign * (j[2 * k - 1] - j[2 * k + 1]) / (2.0 * k as f64);
        }

        // Y₀' = (2/π)[J₀/x - (ln(x/2) + γ) J₁ - 2 Σ (-1)^k (J_{2k-1} - J_{2k+1})/(2k)]
        let y0_prime = FRAC_2_PI * (j[0] / x - ((x / 2.0).ln() + EULER_GAMMA) * j[1] - 2.0 * sum);

        -y0_prime
    } else if x.is_finite() {
        bessel_y_asymptotic(1, x)
    } else {
        0.0
    }
}

/// `Σ (x²/4)^k / (k! (k+ν)!) · (x/2)^ν` for ν = 0, 1: `Iν(x)` by its power
/// series (A&S 9.6.10); all terms are positive so there is no cancellation.
fn bessel_i_series(
    nu: u32,
    x: f64,
) -> f64 {
    let q = x * x / 4.0;

    let mut term = if nu == 0 { 1.0 } else { x / 2.0 };

    let mut sum = term;

    for k in 1..500 {
        let kf = f64::from(k);

        term *= q / (kf * (kf + f64::from(nu)));

        sum += term;

        if term < 1e-17 * sum {
            break;
        }
    }

    sum
}

/// Asymptotic expansion `Iν(x) ~ e^x / sqrt(2πx) Σ (-1)^k a_k(ν) / x^k`
/// (A&S 9.7.1) for large positive `x`.
fn bessel_i_asymptotic(
    nu: u32,
    x: f64,
) -> f64 {
    let mu = 4.0 * f64::from(nu) * f64::from(nu);

    let mut term = 1.0;

    let mut sum = 1.0;

    let mut last = f64::INFINITY;

    for k in 1..100 {
        let kf = f64::from(k);

        term *= -(mu - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf * x);

        if term.abs() > last {
            break;
        }

        last = term.abs();

        sum += term;

        if term.abs() < 1e-17 {
            break;
        }
    }

    // exp(x/2) twice avoids overflow of the intermediate for x near 709.
    let e = (x / 2.0).exp();

    (e / (2.0 * std::f64::consts::PI * x).sqrt()) * e * sum
}

/// Computes the modified Bessel function of the first kind, I₀(x).
///
/// Power series for |x| ≤ 30 and the asymptotic expansion beyond; relative
/// accuracy about 1e-15.
#[must_use]
pub fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();

    if ax <= BESSEL_I_SERIES_MAX {
        bessel_i_series(0, ax)
    } else {
        bessel_i_asymptotic(0, ax)
    }
}

/// Computes the modified Bessel function of the first kind, I₁(x) (odd in `x`).
///
/// Same method and accuracy as [`bessel_i0`].
#[must_use]
pub fn bessel_i1(x: f64) -> f64 {
    let ax = x.abs();

    let v = if ax <= BESSEL_I_SERIES_MAX {
        bessel_i_series(1, ax)
    } else {
        bessel_i_asymptotic(1, ax)
    };

    if x < 0.0 { -v } else { v }
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Computes the Legendre polynomial `P_n(x)` using recurrence relation.
#[must_use]
pub fn legendre_p(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    if n == 1 {
        return x;
    }

    let mut p_prev = 1.0;

    let mut p_curr = x;

    for k in 2..=n {
        let p_next =
            (f64::from(2 * k - 1) * x).mul_add(p_curr, -(f64::from(k - 1) * p_prev)) / f64::from(k);

        p_prev = p_curr;

        p_curr = p_next;
    }

    p_curr
}

/// Computes the Chebyshev polynomial of the first kind `T_n(x)`.
#[must_use]
pub fn chebyshev_t(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    if n == 1 {
        return x;
    }

    let mut t_prev = 1.0;

    let mut t_curr = x;

    for _ in 2..=n {
        let t_next = (2.0 * x).mul_add(t_curr, -t_prev);

        t_prev = t_curr;

        t_curr = t_next;
    }

    t_curr
}

/// Computes the Chebyshev polynomial of the second kind `U_n(x)`.
#[must_use]
pub fn chebyshev_u(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    if n == 1 {
        return 2.0 * x;
    }

    let mut u_prev = 1.0;

    let mut u_curr = 2.0 * x;

    for _ in 2..=n {
        let u_next = (2.0 * x).mul_add(u_curr, -u_prev);

        u_prev = u_curr;

        u_curr = u_next;
    }

    u_curr
}

/// Computes the (physicists') Hermite polynomial `H_n(x)`.
#[must_use]
pub fn hermite_h(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    if n == 1 {
        return 2.0 * x;
    }

    let mut h_prev = 1.0;

    let mut h_curr = 2.0 * x;

    for k in 2..=n {
        let h_next = (2.0 * x).mul_add(h_curr, -(2.0 * f64::from(k - 1) * h_prev));

        h_prev = h_curr;

        h_curr = h_next;
    }

    h_curr
}

/// Computes the Laguerre polynomial `L_n(x)`.
#[must_use]
pub fn laguerre_l(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    if n == 1 {
        return 1.0 - x;
    }

    let mut l_prev = 1.0;

    let mut l_curr = 1.0 - x;

    for k in 2..=n {
        let l_next = (f64::from(2 * k - 1) - x) * l_curr / f64::from(k)
            - f64::from(k - 1) * l_prev / f64::from(k);

        l_prev = l_curr;

        l_curr = l_next;
    }

    l_curr
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes the factorial n!
#[must_use]
pub fn factorial(n: u64) -> f64 {
    if n <= 1 {
        return 1.0;
    }

    gamma((n + 1) as f64)
}

/// Computes the double factorial n!!
/// n!! = n * (n-2) * (n-4) * ... * (1 or 2)
#[must_use]
pub fn double_factorial(n: u64) -> f64 {
    if n <= 1 {
        return 1.0;
    }

    let mut result = 1.0;

    let mut k = n;

    while k > 1 {
        result *= k as f64;

        k -= 2;
    }

    result
}

/// Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!)
#[must_use]
pub fn binomial(
    n: u64,
    k: u64,
) -> f64 {
    if k > n {
        return 0.0;
    }

    if k == 0 || k == n {
        return 1.0;
    }

    // Use gamma for large values to avoid overflow
    gamma((n + 1) as f64) / (gamma((k + 1) as f64) * gamma((n - k + 1) as f64))
}

/// Computes the Riemann zeta function ζ(s) for real s > 1
/// (`+∞` at s = 1, `NaN` for s < 1).
///
/// Uses Euler-Maclaurin summation (via [`hurwitz_zeta`] with q = 1): the first
/// terms directly, the tail by its integral plus Bernoulli-number corrections;
/// accurate to about 1e-15.
#[must_use]
pub fn riemann_zeta(s: f64) -> f64 {
    if s.is_nan() {
        return f64::NAN;
    }

    if (s - 1.0).abs() < 1e-15 {
        return f64::INFINITY;
    }

    if s < 1.0 {
        return f64::NAN; // Analytic continuation is not implemented.
    }

    hurwitz_zeta(s, 1.0)
}

/// Computes the sinc function sinc(x) = sin(πx) / (πx).
#[must_use]
pub fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-10 {
        return 1.0;
    }

    let px = std::f64::consts::PI * x;

    px.sin() / px
}

/// Computes the logit function logit(p) = ln(p / (1-p)).
#[must_use]
pub fn logit(p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return f64::NAN;
    }

    (p / (1.0 - p)).ln()
}

/// Computes the logistic (sigmoid) function σ(x) = 1 / (1 + e^(-x)).
#[must_use]
pub fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Computes the softplus function softplus(x) = ln(1 + e^x).
#[must_use]
pub fn softplus(x: f64) -> f64 {
    if x > 20.0 {
        x // Avoid overflow
    } else if x < -20.0 {
        x.exp()
    } else {
        x.exp().ln_1p()
    }
}

/// Computes the n-th Bernoulli number B_n.
#[must_use]
pub fn bernoulli_number(n: u32) -> f64 {
    // Precomputed first few Bernoulli numbers
    let precomputed = [
        1.0,               // B_0
        -0.5,              // B_1
        1.0 / 6.0,         // B_2
        0.0,               // B_3
        -1.0 / 30.0,       // B_4
        0.0,               // B_5
        1.0 / 42.0,        // B_6
        0.0,               // B_7
        -1.0 / 30.0,       // B_8
        0.0,               // B_9
        5.0 / 66.0,        // B_10
        0.0,               // B_11
        -691.0 / 2730.0,   // B_12
        0.0,               // B_13
        7.0 / 6.0,         // B_14
        0.0,               // B_15
        -3617.0 / 510.0,   // B_16
        0.0,               // B_17
        43867.0 / 798.0,   // B_18
        0.0,               // B_19
        -174_611.0 / 330.0, // B_20
    ];

    if (n as usize) < precomputed.len() {
        return precomputed[n as usize];
    }

    if n % 2 == 1 {
        return 0.0;
    }

    // Dynamic computation using recurrence:
    // B_m = -1/(m+1) * Σ_{k=0}^{m-1} (m+1 choose k) B_k
    let mut b = vec![0.0; (n + 1) as usize];
    b[..precomputed.len()].copy_from_slice(&precomputed[..]);

    for m in precomputed.len()..=(n as usize) {
        if m % 2 == 1 {
            b[m] = 0.0;
            continue;
        }
        let mut sum = 0.0;
        let m_plus_1 = (m + 1) as u64;
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for k in 0..m {
            sum += binomial(m_plus_1, k as u64) * b[k];
        }
        b[m] = -sum / (m_plus_1 as f64);
    }

    b[n as usize]
}

/// Computes the n-th Bernoulli polynomial B_n(x).
#[must_use]
pub fn bernoulli_poly(
    n: u32,
    x: f64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let mut sum = 0.0;
    for k in 0..=n {
        let coeff = binomial(u64::from(n), u64::from(k)) * bernoulli_number(k);
        sum += coeff * x.powi((n - k) as i32);
    }
    sum
}

/// Computes the Hurwitz zeta function ζ(s, q) for real s and q.
/// Defined as Σ_{n=0}^∞ 1/(n+q)^s.
#[must_use]
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn hurwitz_zeta(
    s: f64,
    q: f64,
) -> f64 {
    if q <= 0.0 {
        // If q is a non-positive integer, the function is undefined/infinite.
        if q == q.round() {
            return f64::NAN;
        }
        // Use the recurrence relation: ζ(s, q) = ζ(s, q+1) + q^(-s)
        return hurwitz_zeta(s, q + 1.0) + q.powf(-s);
    }

    if (s - 1.0).abs() < 1e-15 {
        return f64::INFINITY;
    }

    // If s is a non-positive integer, we can use Bernoulli polynomials:
    // ζ(-n, q) = -B_{n+1}(q) / (n+1)
    if s <= 0.0 && s == s.round() {
        let n = (-s) as u32;
        return -bernoulli_poly(n + 1, q) / f64::from(n + 1);
    }

    // For other cases, use Euler-Maclaurin summation.
    // Choose M such that M + q >= 15.0
    let m = if q >= 15.0 {
        0
    } else {
        (15.0 - q).ceil() as usize
    };

    let mut sum = 0.0;
    for j in 0..m {
        sum += (j as f64).mul_add(1.0, q).powf(-s);
    }

    let s_val = (m as f64) + q;

    // Euler-Maclaurin terms:
    // term_1 = S^(1-s) / (s - 1)
    let term1 = s_val.powf(1.0 - s) / (s - 1.0);
    // term_2 = S^(-s) / 2
    let term2 = 0.5 * s_val.powf(-s);

    // Sum the Bernoulli number correction terms
    let bernoulli_coeffs = [
        1.0 / 12.0,                    // B_2 / 2!
        -1.0 / 720.0,                  // B_4 / 4!
        1.0 / 30240.0,                 // B_6 / 6!
        -1.0 / 1_209_600.0,              // B_8 / 8!
        1.0 / 47_900_160.0,              // B_10 / 10!
        -691.0 / 1_307_674_368_000.0,      // B_12 / 12!
        1.0 / 74_724_249_600.0,           // B_14 / 14!
        -3617.0 / 10_670_622_842_880_000.0, // B_16 / 16!
    ];

    let mut correction = 0.0;
    let mut rising_factorial = s;
    let mut rising_next = 1.0;
    let mut s_pow = s_val.powf(-s - 1.0); // S^(-s-1)
    let s_val_sq_inv = 1.0 / (s_val * s_val);

    for &coeff in &bernoulli_coeffs {
        let term = coeff * rising_factorial * s_pow;
        if !term.is_finite() || term.abs() < 1e-16 * (sum + term1 + term2).abs() {
            if term.is_finite() {
                correction += term;
            }
            break;
        }
        correction += term;
        // Next Euler-Maclaurin term needs s (s+1) ... (s + 2j): the two new
        // factors advance with j.
        rising_factorial *= (s + rising_next) * (s + rising_next + 1.0);
        rising_next += 2.0;
        s_pow *= s_val_sq_inv;
    }

    sum + term1 + term2 + correction
}

/// Computes the polygamma function ψ^(n)(z) for real z > 0 and integer n >= 0.
/// The 0-th polygamma is the digamma function.
#[must_use]
pub fn polygamma_numerical(
    n: u32,
    z: f64,
) -> f64 {
    if n == 0 {
        return digamma_numerical(z);
    }

    if z <= 0.0 {
        return f64::NAN;
    }

    let factor = if n.is_multiple_of(2) { -1.0 } else { 1.0 };
    let n_fact = factorial(u64::from(n));
    factor * n_fact * hurwitz_zeta(f64::from(n + 1), z)
}

// ============================================================================
// Bessel Functions of Arbitrary Real Order
// ============================================================================

/// `∫_a^b f` to near machine precision with adaptive Gauss–Kronrod.
fn quad(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
) -> f64 {
    crate::kernels::integrate::gauss_kronrod(f, a, b, 1e-15, 4_000).value
}

/// Upper end for `∫_0^∞ g(t) e^{-x sinh t}` (or `cosh t`): where the
/// exponential has decayed below 1e-30 even against a growth `e^{|ν| t}`.
fn decay_end(
    x: f64,
    nu: f64,
) -> f64 {
    let mut t = 1.0_f64;
    while x * t.sinh() - nu.abs() * t < 80.0 && t < 60.0 {
        t *= 1.5;
    }
    t
}

/// Whether `nu` is an integer.
fn is_integer_order(nu: f64) -> bool {
    nu.fract() == 0.0 && nu.is_finite()
}

/// Bessel function of the first kind `J_ν(x)` for real order `ν` and
/// real `x` (`x < 0` only for integer orders, through `J_n(-x) = (-1)^n
/// J_n(x)`).
///
/// Integer orders 0 and 1 use [`bessel_j0`]/[`bessel_j1`]; other orders
/// Bessel's integral `(1/π)∫_0^π cos(νθ - x sin θ) dθ - (sin νπ/π)
/// ∫_0^∞ e^{-x sinh t - νt} dt`, evaluated by adaptive quadrature —
/// accurate to about 1e-13 absolute for moderate arguments.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn bessel_j(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        if !is_integer_order(nu) {
            return f64::NAN;
        }
        let sign = if nu.rem_euclid(2.0) == 0.0 { 1.0 } else { -1.0 };
        return sign * bessel_j(nu, -x);
    }
    if x == 0.0 {
        return if nu == 0.0 {
            1.0
        } else if nu > 0.0 || is_integer_order(nu) {
            0.0
        } else {
            f64::INFINITY
        };
    }
    if nu == 0.0 {
        return bessel_j0(x);
    }
    if nu == 1.0 {
        return bessel_j1(x);
    }
    if is_integer_order(nu) && nu < 0.0 {
        let sign = if nu.rem_euclid(2.0) == 0.0 { 1.0 } else { -1.0 };
        return sign * bessel_j(-nu, x);
    }
    let pi = std::f64::consts::PI;
    let main = quad(|theta| (nu * theta - x * theta.sin()).cos(), 0.0, pi) / pi;
    if is_integer_order(nu) {
        return main;
    }
    let tail = quad(|t| (-x * t.sinh() - nu * t).exp(), 0.0, decay_end(x, nu));
    (nu * pi).sin().mul_add(-tail / pi, main)
}

/// Bessel function of the second kind `Y_ν(x)` for real order and
/// `x > 0`, by `(1/π)∫_0^π sin(x sin θ - νθ) dθ - (1/π)∫_0^∞ (e^{νt} +
/// e^{-νt} cos νπ) e^{-x sinh t} dt`.
#[must_use]
pub fn bessel_y(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if nu == 0.0 {
        return bessel_y0(x);
    }
    let pi = std::f64::consts::PI;
    let main = quad(|theta| (x * theta.sin() - nu * theta).sin(), 0.0, pi) / pi;
    let c = (nu * pi).cos();
    let tail = quad(|t| ((nu * t).exp() + (-nu * t).exp() * c) * (-x * t.sinh()).exp(), 0.0, decay_end(x, nu));
    main - tail / pi
}

/// Modified Bessel function of the first kind `I_ν(x)`, by
/// `(1/π)∫_0^π e^{x cos θ} cos νθ dθ - (sin νπ/π)∫_0^∞ e^{-x cosh t - νt}
/// dt` (`x < 0` for integer orders by parity).
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn bessel_i(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        if !is_integer_order(nu) {
            return f64::NAN;
        }
        let sign = if nu.rem_euclid(2.0) == 0.0 { 1.0 } else { -1.0 };
        return sign * bessel_i(nu, -x);
    }
    if x == 0.0 {
        return if nu == 0.0 { 1.0 } else if nu > 0.0 || is_integer_order(nu) { 0.0 } else { f64::INFINITY };
    }
    if nu == 0.0 {
        return bessel_i0(x);
    }
    if nu == 1.0 {
        return bessel_i1(x);
    }
    let order = if is_integer_order(nu) { nu.abs() } else { nu };
    let pi = std::f64::consts::PI;
    let main = quad(|theta| (x * theta.cos()).exp() * (order * theta).cos(), 0.0, pi) / pi;
    if is_integer_order(order) {
        return main;
    }
    let tail = quad(|t| (-x * t.cosh() - order * t).exp(), 0.0, decay_end(x, order));
    (order * pi).sin().mul_add(-tail / pi, main)
}

/// Modified Bessel function of the second kind `K_ν(x)` for `x > 0`, by
/// `∫_0^∞ e^{-x cosh t} cosh νt dt`.
#[must_use]
pub fn bessel_k(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::INFINITY;
    }
    quad(|t| (-x * t.cosh()).exp() * (nu * t).cosh(), 0.0, decay_end(x, nu))
}

/// `K₀(x)`.
#[must_use]
pub fn bessel_k0(x: f64) -> f64 {
    bessel_k(0.0, x)
}

/// `K₁(x)`.
#[must_use]
pub fn bessel_k1(x: f64) -> f64 {
    bessel_k(1.0, x)
}

/// The imaginary error function `erfi(x) = -i erf(ix) = (2/√π) ∫_0^x e^{t²}
/// dt`, from its everywhere-convergent power series (all terms of one
/// sign, so no cancellation).
#[must_use]
pub fn erfi(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.abs() > 26.7 {
        return x.signum() * f64::INFINITY;
    }
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    let mut k = 0.0_f64;
    loop {
        k += 1.0;
        term *= x2 / k;
        let contribution = term / 2.0_f64.mul_add(k, 1.0);
        sum += contribution;
        if contribution.abs() <= 1e-17 * sum.abs() {
            break;
        }
    }
    sum * 2.0 / std::f64::consts::PI.sqrt()
}

/// The inverse of `erfc`: `y` with `erfc(y) = p`, `0 < p < 2`. Newton
/// steps on `erfc` refine the inverse of `erf(1 - p)`, which keeps full
/// relative accuracy for small `p`.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn inverse_erfc(p: f64) -> f64 {
    if !(0.0..=2.0).contains(&p) || p.is_nan() {
        return f64::NAN;
    }
    if p == 0.0 {
        return f64::INFINITY;
    }
    if p == 2.0 {
        return f64::NEG_INFINITY;
    }
    let mut y = inverse_erf_numerical(1.0 - p);
    if !y.is_finite() {
        // 1 - p rounded to ±1: start from the tail asymptotics.
        y = (-(p / 2.0).ln()).sqrt() * if p < 1.0 { 1.0 } else { -1.0 };
    }
    for _ in 0..8 {
        let f = erfc_numerical(y) - p;
        let derivative = -2.0 / std::f64::consts::PI.sqrt() * (-y * y).exp();
        if derivative == 0.0 {
            break;
        }
        let step = f / derivative;
        y -= step;
        if step.abs() <= 1e-16 * y.abs().max(1.0) {
            break;
        }
    }
    y
}

/// `ln(n!)` for real `n > -1`, through `ln Γ(n + 1)`.
#[must_use]
pub fn ln_factorial(n: f64) -> f64 {
    ln_gamma_numerical(n + 1.0)
}

/// The generalised Laguerre polynomial `L_n^{(α)}(x)` by the three-term
/// recurrence.
#[must_use]
pub fn generalized_laguerre(
    n: u32,
    alpha: f64,
    x: f64,
) -> f64 {
    let mut previous = 1.0;
    if n == 0 {
        return previous;
    }
    let mut current = 1.0 + alpha - x;
    for k in 1..n {
        let k = f64::from(k);
        let next = ((2.0f64.mul_add(k, 1.0) + alpha - x) * current - (k + alpha) * previous) / (k + 1.0);
        previous = current;
        current = next;
    }
    current
}

#[cfg(test)]
mod arbitrary_order_bessel_tests {
    use super::*;

    #[test]
    fn integer_orders_against_reference_values() {
        // Abramowitz & Stegun tables (ten digits).
        let close = |a: f64, b: f64| (a - b).abs() < 1e-9;
        assert!(close(bessel_j(2.0, 1.0), 0.114_903_484_9));
        assert!(close(bessel_j(5.0, 10.0), -0.234_061_528_2));
        assert!(close(bessel_j(1.0, 3.0), 0.339_058_958_5));
        assert!(close(bessel_y(1.0, 1.0), -0.781_212_821_3));
        assert!(close(bessel_i(2.0, 1.0), 0.135_747_669_8));
        assert!(close(bessel_k0(1.0), 0.421_024_438_2));
        assert!(close(bessel_k1(2.0), 0.139_865_881_8));
        // Recurrences: C_{ν+1} = (2ν/x) C_ν - C_{ν-1} for J and Y,
        // K_{ν+1} = K_{ν-1} + (2ν/x) K_ν, I_{ν+1} = I_{ν-1} - (2ν/x) I_ν.
        for &(nu, x) in &[(1.0, 5.0), (2.5, 3.0), (3.0, 0.7)] {
            let r = 2.0 * nu / x;
            assert!((bessel_y(nu + 1.0, x) - (r * bessel_y(nu, x) - bessel_y(nu - 1.0, x))).abs() < 1e-10 * bessel_y(nu + 1.0, x).abs().max(1.0));
            assert!((bessel_j(nu + 1.0, x) - (r * bessel_j(nu, x) - bessel_j(nu - 1.0, x))).abs() < 1e-12);
            assert!((bessel_k(nu + 1.0, x) - (bessel_k(nu - 1.0, x) + r * bessel_k(nu, x))).abs() < 1e-10 * bessel_k(nu + 1.0, x));
            assert!((bessel_i(nu + 1.0, x) - (bessel_i(nu - 1.0, x) - r * bessel_i(nu, x))).abs() < 1e-12);
        }
        // Parity for negative arguments and orders.
        assert!((bessel_j(3.0, -2.0) + bessel_j(3.0, 2.0)).abs() < 1e-15);
        assert!((bessel_j(-3.0, 2.0) + bessel_j(3.0, 2.0)).abs() < 1e-14);
    }

    #[test]
    fn half_integer_orders_are_elementary() {
        for &x in &[0.3, 1.0, 4.5, 12.0] {
            let s = (2.0 / (std::f64::consts::PI * x)).sqrt();
            assert!((bessel_j(0.5, x) - s * x.sin()).abs() < 1e-12, "J(1/2, {x})");
            assert!((bessel_j(-0.5, x) - s * x.cos()).abs() < 1e-12, "J(-1/2, {x})");
            assert!((bessel_y(0.5, x) + s * x.cos()).abs() < 1e-12, "Y(1/2, {x})");
            assert!((bessel_i(0.5, x) - s * x.sinh()).abs() < 1e-11 * x.sinh().max(1.0), "I(1/2, {x})");
            let k = (std::f64::consts::PI / (2.0 * x)).sqrt() * (-x).exp();
            assert!((bessel_k(0.5, x) - k).abs() < 1e-13, "K(1/2, {x})");
        }
    }

    #[test]
    fn wronskian() {
        // J_ν Y_ν' - J_ν' Y_ν = 2/(πx) through J_{ν+1} Y_ν - J_ν Y_{ν+1}.
        for &(nu, x) in &[(0.3, 2.0), (2.5, 7.0), (4.0, 1.5)] {
            let w = bessel_j(nu + 1.0, x) * bessel_y(nu, x) - bessel_j(nu, x) * bessel_y(nu + 1.0, x);
            assert!((w - 2.0 / (std::f64::consts::PI * x)).abs() < 1e-11, "ν = {nu}, x = {x}: {w}");
        }
    }

    #[test]
    fn error_function_relatives() {
        assert!((erfi(1.0) - 1.650_425_758_797_542_8).abs() < 1e-14);
        assert!((erfi(-0.5) + 0.614_952_094_696_510_9).abs() < 1e-14);
        for &p in &[1e-10, 0.01, 0.5, 1.0, 1.7] {
            assert!((erfc_numerical(inverse_erfc(p)) - p).abs() < 1e-12 * p.max(1e-3), "{p}");
        }
        assert!((ln_factorial(10.0) - 3_628_800.0_f64.ln()).abs() < 1e-12);
        // L_2^{(1)}(x) = (x² - 6x + 6)/2.
        assert!((generalized_laguerre(2, 1.0, 0.7) - f64::midpoint(0.49 - 4.2, 6.0)).abs() < 1e-14);
    }
}

/// The Lambert W function: `W(x) e^{W(x)} = x`.
///
/// `upper = true` is the principal branch `W₀` (defined for `x ≥ -1/e`); `upper = false` the
/// branch `W₋₁` (for `-1/e ≤ x < 0`). NaN outside the domain. Halley's
/// iteration from a branch-appropriate start.
#[must_use]
pub fn lambert_w(
    x: f64,
    upper: bool,
) -> f64 {
    let branch_point = -(-1.0_f64).exp();
    if !x.is_finite() || x < branch_point || (!upper && x >= 0.0) {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    if (x - branch_point).abs() < 1e-300 {
        return -1.0;
    }
    // Starting points: the series about the branch point near -1/e, a
    // logarithmic estimate elsewhere.
    let p = (2.0 * (1.0 + std::f64::consts::E * x)).max(0.0).sqrt();
    let mut w = if upper {
        if x < 0.25 {
            -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p
        } else {
            let l = x.ln_1p();
            l - l.max(1e-300).ln_1p() * 0.5
        }
    } else if x > -0.25 {
        let l = (-x).ln();
        l - (-l).ln()
    } else {
        -1.0 - p - p * p / 3.0 - 11.0 / 72.0 * p * p * p
    };
    for _ in 0..64 {
        let e = w.exp();
        let f = w * e - x;
        let wp1 = w + 1.0;
        if wp1.abs() < 1e-300 {
            break;
        }
        let step = f / (e * wp1 - (w + 2.0) * f / (2.0 * wp1));
        w -= step;
        if step.abs() <= 1e-15 * w.abs().max(1.0) {
            break;
        }
    }
    w
}

#[cfg(test)]
mod lambert_tests {
    use super::lambert_w;

    #[test]
    fn lambert_w_branches() {
        for &x in &[-0.3, -0.1, 0.5, 1.0, 10.0, 1e6] {
            let w = lambert_w(x, true);
            assert!((w * w.exp() - x).abs() < 1e-10 * x.abs().max(1.0), "W0({x}) = {w}");
        }
        assert!((lambert_w(1.0, true) - 0.567_143_290_409_784).abs() < 1e-14);
        for &x in &[-0.35, -0.2, -0.01] {
            let w = lambert_w(x, false);
            assert!(w <= -1.0 && (w * w.exp() - x).abs() < 1e-10, "W-1({x}) = {w}");
        }
        assert!(lambert_w(-1.0, true).is_nan() && lambert_w(0.5, false).is_nan());
    }
}
