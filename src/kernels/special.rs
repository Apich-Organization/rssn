//! # Numerical Special Functions
//!
//! This module provides numerical implementations of mathematical special functions.
//! It includes gamma functions, beta functions, error functions, Bessel functions,
//! orthogonal polynomials, and other commonly used special functions.
//!
//! ## Limitations
//!
//! [`bessel_j`] and [`bessel_y`] use Bessel's integral, evaluated by adaptive
//! quadrature, outside the uniform region `ν >= 50` with
//! `|x/ν − 1| <= 1/4`. There the integrand oscillates `O(ν)` times and the
//! quadrature cost grows with `ν`; inside that region both functions switch to
//! the uniform Airy-type (Debye) expansion of DLMF 10.20.4-10.20.5, which is
//! accurate for large order and through the turning point `x = ν`. What remains
//! limited: (i) the quadrature branch is unchanged and still costs `O(ν)`
//! integrand evaluations outside the band; (ii) [`bessel_j`]/[`bessel_y`] accept
//! real `ν` only — complex order is not implemented; (iii) the Airy expansion
//! is evaluated with `A_0..A_2` and `B_0..B_2`, so its relative error grows
//! like `ν^{-6}` in the worst case (about `1e-14` at `ν = 50`, `1e-16` at
//! `ν = 100`); (iv) for `x` far below `ν` the exponentially small `J_ν` and
//! exponentially large `Y_ν` underflow/overflow in `f64` regardless of method.

use std::f64::consts::FRAC_2_PI;

use num_complex::Complex64;

use statrs::function::beta::beta;
use statrs::function::beta::ln_beta;
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
    if x < 0.0 && x.is_finite() {
        // ln |Gamma(x)| for negative x by the reflection formula (+inf at the poles)
        let s = (std::f64::consts::PI * x).sin().abs();
        if crate::kernels::exact_eq(x, x.floor()) || crate::kernels::exact_eq(s, 0.0) {
            return f64::INFINITY;
        }
        return (std::f64::consts::PI / s).ln() - ln_gamma(1.0 - x);
    }
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
    if x.is_nan() {
        return f64::NAN;
    }
    if x.abs() < 1.5 {
        return erf_series(x);
    }
    let c = erfc_tail(x.abs());
    if x > 0.0 { 1.0 - c } else { c - 1.0 }
}

/// Maclaurin-type series `erf(x) = 2/sqrt(pi) e^(-x^2) sum 2^n x^(2n+1)/(2n+1)!!`
/// with only positive terms (no cancellation), accurate to a few ulps for `|x| < 1.5`.
fn erf_series(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for n in 1..200 {
        term *= 2.0 * x2 / f64::from(2 * n + 1);
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * (-x2).exp() * sum
}

/// `erfc(x)` for `x >= 1.5` by the continued fraction
/// `e^(-x^2)/sqrt(pi) / (x + (1/2)/(x + 1/(x + (3/2)/(x + ...))))`.
fn erfc_tail(x: f64) -> f64 {
    let mut f = x;
    for k in (1..=300).rev() {
        f = x + f64::from(k) * 0.5 / f;
    }
    (-x * x).exp() / (std::f64::consts::PI.sqrt() * f)
}

/// Computes the complementary error function, erfc(x) = 1 - erf(x).
#[must_use]
pub fn erfc_numerical(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 2.0 - erfc_numerical(-x);
    }
    if x < 1.5 {
        1.0 - erf_series(x)
    } else {
        erfc_tail(x)
    }
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
            (erfc_numerical(y) - (1.0 - x)) / -deriv
        } else {
            (erf_numerical(y) - x) / deriv
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
const EULER_GAMMA: f64 = std::f64::consts::EULER_GAMMA;

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
    if x.is_infinite() {
        return f64::INFINITY;
    }
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
    if x.is_infinite() {
        return x;
    }
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

/// Computes the Riemann zeta function ζ(s) for real `s`.
///
/// * `s > 1`: Euler-Maclaurin summation (via [`hurwitz_zeta`] with q = 1):
///   the first terms directly, the tail by its integral plus
///   Bernoulli-number corrections; accurate to about 1e-15.
/// * `s = 1`: `+∞` (the pole).
/// * `-1/2 < s < 1`: the alternating eta series
///   `η(s) = (1 - 2^(1-s)) ζ(s)` accelerated with Borwein's algorithm
///   (relative error about `1e-15`); `s = 0` is exactly `-1/2`.
/// * `s <= -1/2`: the analytic continuation by the functional equation
///   `ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)`. The trivial zeros at the
///   negative even integers are returned as exactly `0`, and
///   `ζ(-1) = -1/12`, `ζ(-3) = 1/120`, ... follow from the same formula.
///
/// Returns `NaN` for `NaN` or `-∞`, and `±∞` where the value overflows
/// (roughly `s < -170`).
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against sentinel / integer-valued input is intended
pub fn riemann_zeta(s: f64) -> f64 {
    if s.is_nan() || s == f64::NEG_INFINITY {
        return f64::NAN;
    }

    if (s - 1.0).abs() < 1e-15 {
        return f64::INFINITY;
    }

    if s > 1.0 {
        return hurwitz_zeta(s, 1.0);
    }

    if s == 0.0 {
        return -0.5;
    }

    if s > -0.5 {
        return zeta_via_eta(s);
    }

    // s <= -0.5: trivial zeros at negative even integers.
    if s == s.round() && (s / 2.0) == (s / 2.0).round() {
        return 0.0;
    }

    // sin(pi s / 2) with exact argument reduction modulo 4.
    let r = s - 4.0 * (s / 4.0).floor();
    let sin_part = (std::f64::consts::FRAC_PI_2 * r).sin();
    let one_minus_s = 1.0 - s;
    let zeta_reflected = hurwitz_zeta(one_minus_s, 1.0);
    let g = gamma(one_minus_s);
    if g.is_infinite() {
        // Overflow: only the sign is meaningful.
        return (sin_part.signum() * zeta_reflected.signum()) * f64::INFINITY;
    }
    (2.0_f64).powf(s) * std::f64::consts::PI.powf(s - 1.0) * sin_part * g * zeta_reflected
}

/// ζ(s) for `-1/2 < s < 1` through Borwein's accelerated alternating series for
/// the Dirichlet eta function.
fn zeta_via_eta(s: f64) -> f64 {
    const N: usize = 32;
    // d_k = N * sum_{i<=k} (N+i-1)! 4^i / ((N-i)! (2i)!), built from term ratios.
    let n = N as f64;
    let mut d = [0.0_f64; N + 1];
    let mut term = 1.0_f64;
    let mut acc = 1.0_f64;
    d[0] = acc;
    for i in 0..N {
        let fi = i as f64;
        term *= (n + fi) * 4.0 * (n - fi) / ((2.0 * fi + 1.0) * (2.0 * fi + 2.0));
        acc += term;
        d[i + 1] = acc;
    }
    let dn = d[N];
    let mut sum = 0.0;
    for (k, dk) in d.iter().enumerate().take(N) {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        sum += sign * (dk - dn) / ((k + 1) as f64).powf(s);
    }
    let eta = -sum / dn;
    // 1 - 2^(1-s) computed without cancellation near s = 1
    let denom = -((1.0 - s) * std::f64::consts::LN_2).exp_m1();
    eta / denom
}

/// Computes the sinc function sinc(x) = sin(πx) / (πx).
#[must_use]
pub fn sinc(x: f64) -> f64 {
    if x.is_infinite() {
        return 0.0;
    }
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

/// Computes the n-th Bernoulli number `B_n`.
#[must_use]
pub fn bernoulli_number(n: u32) -> f64 {
    // Precomputed first few Bernoulli numbers
    let precomputed = [
        1.0,                // B_0
        -0.5,               // B_1
        1.0 / 6.0,          // B_2
        0.0,                // B_3
        -1.0 / 30.0,        // B_4
        0.0,                // B_5
        1.0 / 42.0,         // B_6
        0.0,                // B_7
        -1.0 / 30.0,        // B_8
        0.0,                // B_9
        5.0 / 66.0,         // B_10
        0.0,                // B_11
        -691.0 / 2730.0,    // B_12
        0.0,                // B_13
        7.0 / 6.0,          // B_14
        0.0,                // B_15
        -3617.0 / 510.0,    // B_16
        0.0,                // B_17
        43867.0 / 798.0,    // B_18
        0.0,                // B_19
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
        #[allow(clippy::needless_range_loop)]
        // index is used for more than one array / arithmetic; iterator form would not be clearer
        for k in 0..m {
            sum += binomial(m_plus_1, k as u64) * b[k];
        }
        b[m] = -sum / (m_plus_1 as f64);
    }

    b[n as usize]
}

/// Computes the n-th Bernoulli polynomial `B_n(x)`.
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
        1.0 / 12.0,                         // B_2 / 2!
        -1.0 / 720.0,                       // B_4 / 4!
        1.0 / 30240.0,                      // B_6 / 6!
        -1.0 / 1_209_600.0,                 // B_8 / 8!
        1.0 / 47_900_160.0,                 // B_10 / 10!
        -691.0 / 1_307_674_368_000.0,       // B_12 / 12!
        1.0 / 74_724_249_600.0,             // B_14 / 14!
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

    let factor = if n.is_multiple_of(2) {
        -1.0
    } else {
        1.0
    };
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
        let sign = if nu.rem_euclid(2.0) == 0.0 {
            1.0
        } else {
            -1.0
        };
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
        let sign = if nu.rem_euclid(2.0) == 0.0 {
            1.0
        } else {
            -1.0
        };
        return sign * bessel_j(-nu, x);
    }
    if let Some((j, _)) = uniform_bessel_jy(nu, x) {
        return j;
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
    if let Some((_, y)) = uniform_bessel_jy(nu, x) {
        return y;
    }
    let pi = std::f64::consts::PI;
    let main = quad(|theta| (x * theta.sin() - nu * theta).sin(), 0.0, pi) / pi;
    let c = (nu * pi).cos();
    let tail = quad(
        |t| ((nu * t).exp() + (-nu * t).exp() * c) * (-x * t.sinh()).exp(),
        0.0,
        decay_end(x, nu),
    );
    main - tail / pi
}

// ----------------------------------------------------------------------------
// Uniform Airy-type (Debye) expansion for large order
// ----------------------------------------------------------------------------

/// Debye polynomial `U_0(p) = 1` (DLMF 10.41.10), coefficients ascending in `p`.
const DEBYE_U0: &[f64] = &[1.0];

/// Debye polynomial `U_1(p) = (3p − 5p³)/24`.
const DEBYE_U1: &[f64] = &[0.0, 0.125, 0.0, -0.208_333_333_333_333_34];

/// Debye polynomial `U_2(p) = (81p² − 462p⁴ + 385p⁶)/1152`.
const DEBYE_U2: &[f64] = &[
    0.0,
    0.0,
    0.070_312_5,
    0.0,
    -0.401_041_666_666_666_7,
    0.0,
    0.334_201_388_888_888_9,
];

/// Debye polynomial `U_3(p)` (DLMF 10.41.10).
const DEBYE_U3: &[f64] = &[
    0.0,
    0.0,
    0.0,
    0.073_242_187_5,
    0.0,
    -0.891_210_937_5,
    0.0,
    1.846_462_673_611_111_2,
    0.0,
    -1.025_812_596_450_617_3,
];

/// Debye polynomial `U_4(p)` (from the DLMF 10.41.9 recurrence).
const DEBYE_U4: &[f64] = &[
    0.0,
    0.0,
    0.0,
    0.0,
    0.112_152_099_609_375,
    0.0,
    -2.364_086_914_062_5,
    0.0,
    8.789_123_535_156_25,
    0.0,
    -11.207_002_616_222_994,
    0.0,
    4.669_584_423_426_247,
];

/// Debye polynomial `U_5(p)` (from the DLMF 10.41.9 recurrence).
const DEBYE_U5: &[f64] = &[
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.227_108_001_708_984_4,
    0.0,
    -7.368_794_359_479_632,
    0.0,
    42.534_998_745_388_46,
    0.0,
    -91.818_241_543_24,
    0.0,
    84.636_217_674_600_73,
    0.0,
    -28.212_072_558_200_245,
];

/// Airy asymptotic coefficient `u_1 = 5/72` (DLMF 9.7.2).
const AIRY_U1: f64 = 0.069_444_444_444_444_45;

/// Airy asymptotic coefficient `u_2 = 385/10368`.
const AIRY_U2: f64 = 0.037_133_487_654_320_99;

/// Airy asymptotic coefficient `u_3 = 85085/2239488`.
const AIRY_U3: f64 = 0.037_993_059_127_800_64;

/// Airy asymptotic coefficient `u_4 = 37182145/644972544`.
const AIRY_U4: f64 = 0.057_649_190_412_669_72;

/// Airy asymptotic coefficient `u_5 = 5391411025/46438023168`.
const AIRY_U5: f64 = 0.116_099_064_025_515_41;

/// Airy asymptotic coefficient `v_1 = −7/72`.
const AIRY_V1: f64 = -0.097_222_222_222_222_22;

/// Airy asymptotic coefficient `v_2 = −455/10368`.
const AIRY_V2: f64 = -0.043_885_030_864_197_53;

/// Airy asymptotic coefficient `v_3 = −95095/2239488`.
const AIRY_V3: f64 = -0.042_462_830_789_894_83;

/// Airy asymptotic coefficient `v_4 = −40415375/644972544`.
const AIRY_V4: f64 = -0.062_662_163_492_032_31;

/// Local expansion coefficients of `A_1(ζ(w))` about `w = 1 − z = 0`
/// (DLMF 10.20.10), used only for `|w| < 0.05` where the defining sum
/// cancels catastrophically.
const A1_TAYLOR: &[f64] = &[
    -0.004_444_444_444_444_444,
    -0.001_844_155_844_155_844_2,
    0.000_568_120_768_120_768_1,
    0.001_681_378_656_616_743_6,
    0.001_867_440_421_392_210_7,
    0.001_613_301_057_994_372_2,
    0.001_231_773_155_394_427_6,
];

/// Local expansion coefficients of `A_2(ζ(w))` about `w = 0`.
const A2_TAYLOR: &[f64] = &[
    0.000_693_735_541_354_589,
    0.000_464_483_490_365_843_3,
    -0.000_289_036_254_605_598_1,
    -0.000_874_764_943_953_730_5,
    -0.001_029_716_376_135_103_6,
    -0.000_836_857_330_452_832_8,
    -0.000_488_910_822_108_250_4,
];

/// Local expansion coefficients of `B_0(ζ(w))` about `w = 0`
/// (`B_0(0) = 2^{1/3}/70`).
const B0_TAYLOR: &[f64] = &[
    0.017_998_872_141_355_33,
    0.011_199_298_221_287_762,
    0.005_940_406_978_601_43,
    0.002_867_672_451_639_014,
    0.001_233_918_905_254_176,
    0.000_416_925_067_849_407_2,
    0.000_033_017_300_252_851_02,
];

/// Local expansion coefficients of `B_1(ζ(w))` about `w = 0`.
const B1_TAYLOR: &[f64] = &[
    -0.001_492_829_532_134_291_7,
    -0.001_756_409_419_092_778_6,
    -0.001_133_461_488_741_749_1,
    -0.000_346_910_909_813_837_3,
    0.000_227_525_161_050_399_68,
    0.000_517_641_456_930_732_4,
    0.000_589_061_778_735_897_8,
];

/// Local expansion coefficients of `B_2(ζ(w))` about `w = 0`.
const B2_TAYLOR: &[f64] = &[
    0.000_552_213_076_721_292_8,
    0.000_895_865_163_104_769_3,
    0.000_670_150_034_415_699_5,
    0.000_101_662_633_619_006_3,
    -0.000_440_863_451_261_136_6,
    -0.000_739_630_822_719_961,
    -0.000_767_454_454_768_502_5,
];

/// Evaluates a real polynomial given in ascending coefficient order.
fn horner_real(
    coeffs: &[f64],
    x: f64,
) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc.mul_add(x, c))
}

/// Evaluates a complex polynomial given in ascending coefficient order.
fn horner_complex(
    coeffs: &[f64],
    x: Complex64,
) -> Complex64 {
    coeffs
        .iter()
        .rev()
        .fold(Complex64::new(0.0, 0.0), |acc, &c| {
            acc * x + Complex64::new(c, 0.0)
        })
}

/// Debye polynomial `U_k(p)` for `k <= 5`.
fn debye_u(
    k: usize,
    p: Complex64,
) -> Complex64 {
    let coeffs = match k {
        | 0 => DEBYE_U0,
        | 1 => DEBYE_U1,
        | 2 => DEBYE_U2,
        | 3 => DEBYE_U3,
        | 4 => DEBYE_U4,
        | _ => DEBYE_U5,
    };
    horner_complex(coeffs, p)
}

/// `ζ(z)` of DLMF 10.20.2-10.20.3 for `z = x/ν`.
///
/// The two integrals lose all digits by cancellation as `z → 1`; the first
/// three terms of the local series `ζ = 2^{1/3} w (1 + 3w/10 + 32w²/175)`
/// are used there instead.
fn debye_zeta(
    z: f64,
    w: f64,
) -> f64 {
    if w.abs() < 3.0e-3 {
        let c = 2.0_f64.cbrt();
        let w2 = w * w;
        return c
            * w
            * (1.0
                + 0.3 * w
                + 0.182_857_142_857_142_86 * w2
                + 0.131_682_539_682_539_68 * w2 * w
                + 0.102_636_487_322_201_61 * w2 * w2);
    }
    if z < 1.0 {
        let s = (1.0 / z).acosh() - (1.0 - z * z).sqrt();
        return (1.5 * s).powf(2.0 / 3.0);
    }
    let s = (z * z - 1.0).sqrt() - (1.0 / z).acos();
    -(1.5 * s).powf(2.0 / 3.0)
}

/// The five nontrivial coefficients `A_1, A_2, B_0, B_1, B_2` of the uniform
/// expansion at scaled argument `z = x/ν` (`w = 1 − z`).
///
/// For `|w| < 0.05` the local polynomials above are used. Outside, the
/// defining sums of DLMF 10.20.10-10.20.11 are evaluated in `f64` complex
/// arithmetic: `p = (1 − z²)^{-1/2}` is imaginary for `z > 1`, and the
/// imaginary parts of the two sums cancel.
fn debye_coefficients(
    w: f64,
    zeta: f64,
) -> (f64, f64, f64, f64, f64) {
    if w.abs() < 0.05 {
        return (
            horner_real(A1_TAYLOR, w),
            horner_real(A2_TAYLOR, w),
            horner_real(B0_TAYLOR, w),
            horner_real(B1_TAYLOR, w),
            horner_real(B2_TAYLOR, w),
        );
    }
    let zc = Complex64::new(zeta, 0.0);
    let sqrt_zeta = zc.sqrt();
    let zm = Complex64::new(1.0, 0.0) / (zc * sqrt_zeta); // ζ^{-3/2}
    let zm2 = zm * zm;
    let zm3 = zm2 * zm;
    let zm4 = zm3 * zm;
    let zm5 = zm4 * zm;
    let zis = Complex64::new(1.0, 0.0) / sqrt_zeta; // ζ^{-1/2}
    let p = Complex64::new(1.0, 0.0) / Complex64::new(w * (2.0 - w), 0.0).sqrt();
    let (u0, u1, u2, u3, u4, u5) = (
        debye_u(0, p),
        debye_u(1, p),
        debye_u(2, p),
        debye_u(3, p),
        debye_u(4, p),
        debye_u(5, p),
    );
    let a1 = u2 + zm * u1 * AIRY_V1 * 1.5 + zm2 * u0 * AIRY_V2 * 2.25;
    let a2 = u4
        + zm * u3 * AIRY_V1 * 1.5
        + zm2 * u2 * AIRY_V2 * 2.25
        + zm3 * u1 * AIRY_V3 * 3.375
        + zm4 * u0 * AIRY_V4 * 5.062_5;
    let b0 = -zis * (u1 + zm * u0 * AIRY_U1 * 1.5);
    let b1 = -zis
        * (u3 + zm * u2 * AIRY_U1 * 1.5 + zm2 * u1 * AIRY_U2 * 2.25 + zm3 * u0 * AIRY_U3 * 3.375);
    let b2 = -zis
        * (u5
            + zm * u4 * AIRY_U1 * 1.5
            + zm2 * u3 * AIRY_U2 * 2.25
            + zm3 * u2 * AIRY_U3 * 3.375
            + zm4 * u1 * AIRY_U4 * 5.062_5
            + zm5 * u0 * AIRY_U5 * 7.593_75);
    (a1.re, a2.re, b0.re, b1.re, b2.re)
}

/// Uniform Airy-type (Debye) expansion of `J_ν(x)` and `Y_ν(x)` (DLMF
/// 10.20.4-10.20.5), valid for large `ν` and uniformly through `x = ν`.
///
/// Returns `None` outside `ν >= 50`, `x > 0`, `|x/ν − 1| <= 1/4`, where the
/// quadrature branches of [`bessel_j`]/[`bessel_y`] are used instead. Inside
/// the band the expansion is accurate to about `ν^{-6}` (worst case) and is
/// the only branch that remains accurate as `ν` grows: the quadrature needs
/// `O(ν)` oscillations and eventually loses all digits.
#[allow(clippy::float_cmp)] // w == 0 is exactly the turning point x = nu
fn uniform_bessel_jy(
    nu: f64,
    x: f64,
) -> Option<(f64, f64)> {
    if nu.is_nan() || nu < 50.0 || !nu.is_finite() || x <= 0.0 {
        return None;
    }
    // Compute w = 1 - x/ν from (ν - x)/ν: the direct form loses every digit
    // of w as x -> ν (Sterbenz's lemma makes the subtraction exact).
    let w = (nu - x) / nu;
    let z = 1.0 - w;
    if w.abs() > 0.25 {
        return None;
    }
    let (zeta, pref) = if w == 0.0 {
        (0.0, 2.0_f64.cbrt())
    } else {
        let zeta = debye_zeta(z, w);
        let pref = (4.0 * zeta / (w * (2.0 - w))).powf(0.25);
        (zeta, pref)
    };
    let (a1, a2, b0, b1, b2) = debye_coefficients(w, zeta);
    let nu2 = nu * nu;
    let sum_a = 1.0 + a1 / nu2 + a2 / (nu2 * nu2);
    let sum_b = b0 + b1 / nu2 + b2 / (nu2 * nu2);
    let tau = nu.powf(2.0 / 3.0) * zeta;
    // The K-function representation used by `airy_ai` for x > 0 loses its
    // meaning as x -> 0+ (the integrand stops decaying); the Taylor series
    // from the origin is exact there and is used for |tau| < 1/2.
    let (ai, aip, bi, bip) = if tau.abs() < 0.5 {
        let (ai, aip) = airy_taylor(tau, AIRY_AI0, -AIRY_AIP0);
        let (bi, bip) = airy_taylor(tau, 3.0_f64.sqrt() * AIRY_AI0, 3.0_f64.sqrt() * AIRY_AIP0);
        (ai, aip, bi, bip)
    } else {
        (
            airy_ai(tau),
            airy_ai_prime(tau),
            airy_bi(tau),
            airy_bi_prime(tau),
        )
    };
    let inv_cbrt = 1.0 / nu.cbrt();
    let inv_five_thirds = 1.0 / nu.powf(5.0 / 3.0);
    let j = pref * (ai * inv_cbrt * sum_a + aip * inv_five_thirds * sum_b);
    let y = -pref * (bi * inv_cbrt * sum_a + bip * inv_five_thirds * sum_b);
    if j.is_nan() || y.is_nan() {
        return None;
    }
    Some((j, y))
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
        let sign = if nu.rem_euclid(2.0) == 0.0 {
            1.0
        } else {
            -1.0
        };
        return sign * bessel_i(nu, -x);
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
        return bessel_i0(x);
    }
    if nu == 1.0 {
        return bessel_i1(x);
    }
    let order = if is_integer_order(nu) {
        nu.abs()
    } else {
        nu
    };
    let pi = std::f64::consts::PI;
    let main = quad(
        |theta| (x * theta.cos()).exp() * (order * theta).cos(),
        0.0,
        pi,
    ) / pi;
    if is_integer_order(order) {
        return main;
    }
    let tail = quad(
        |t| (-x * t.cosh() - order * t).exp(),
        0.0,
        decay_end(x, order),
    );
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
    quad(
        |t| (-x * t.cosh()).exp() * (nu * t).cosh(),
        0.0,
        decay_end(x, nu),
    )
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
        let next =
            ((2.0f64.mul_add(k, 1.0) + alpha - x) * current - (k + alpha) * previous) / (k + 1.0);
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
            assert!(
                (bessel_y(nu + 1.0, x) - (r * bessel_y(nu, x) - bessel_y(nu - 1.0, x))).abs()
                    < 1e-10 * bessel_y(nu + 1.0, x).abs().max(1.0)
            );
            assert!(
                (bessel_j(nu + 1.0, x) - (r * bessel_j(nu, x) - bessel_j(nu - 1.0, x))).abs()
                    < 1e-12
            );
            assert!(
                (bessel_k(nu + 1.0, x) - (bessel_k(nu - 1.0, x) + r * bessel_k(nu, x))).abs()
                    < 1e-10 * bessel_k(nu + 1.0, x)
            );
            assert!(
                (bessel_i(nu + 1.0, x) - (bessel_i(nu - 1.0, x) - r * bessel_i(nu, x))).abs()
                    < 1e-12
            );
        }
        // Parity for negative arguments and orders.
        assert!((bessel_j(3.0, -2.0) + bessel_j(3.0, 2.0)).abs() < 1e-15);
        assert!((bessel_j(-3.0, 2.0) + bessel_j(3.0, 2.0)).abs() < 1e-14);
    }

    #[test]
    fn half_integer_orders_are_elementary() {
        for &x in &[0.3, 1.0, 4.5, 12.0] {
            let s = (2.0 / (std::f64::consts::PI * x)).sqrt();
            assert!(
                (bessel_j(0.5, x) - s * x.sin()).abs() < 1e-12,
                "J(1/2, {x})"
            );
            assert!(
                (bessel_j(-0.5, x) - s * x.cos()).abs() < 1e-12,
                "J(-1/2, {x})"
            );
            assert!(
                (bessel_y(0.5, x) + s * x.cos()).abs() < 1e-12,
                "Y(1/2, {x})"
            );
            assert!(
                (bessel_i(0.5, x) - s * x.sinh()).abs() < 1e-11 * x.sinh().max(1.0),
                "I(1/2, {x})"
            );
            let k = (std::f64::consts::PI / (2.0 * x)).sqrt() * (-x).exp();
            assert!((bessel_k(0.5, x) - k).abs() < 1e-13, "K(1/2, {x})");
        }
    }

    #[test]
    fn wronskian() {
        // J_ν Y_ν' - J_ν' Y_ν = 2/(πx) through J_{ν+1} Y_ν - J_ν Y_{ν+1}.
        for &(nu, x) in &[(0.3, 2.0), (2.5, 7.0), (4.0, 1.5)] {
            let w =
                bessel_j(nu + 1.0, x) * bessel_y(nu, x) - bessel_j(nu, x) * bessel_y(nu + 1.0, x);
            assert!(
                (w - 2.0 / (std::f64::consts::PI * x)).abs() < 1e-11,
                "ν = {nu}, x = {x}: {w}"
            );
        }
    }

    #[test]
    fn error_function_relatives() {
        assert!((erfi(1.0) - 1.650_425_758_797_542_8).abs() < 1e-14);
        assert!((erfi(-0.5) + 0.614_952_094_696_510_9).abs() < 1e-14);
        for &p in &[1e-10, 0.01, 0.5, 1.0, 1.7] {
            assert!(
                (erfc_numerical(inverse_erfc(p)) - p).abs() < 1e-12 * p.max(1e-3),
                "{p}"
            );
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
            assert!(
                (w * w.exp() - x).abs() < 1e-10 * x.abs().max(1.0),
                "W0({x}) = {w}"
            );
        }
        assert!((lambert_w(1.0, true) - 0.567_143_290_409_784).abs() < 1e-14);
        for &x in &[-0.35, -0.2, -0.01] {
            let w = lambert_w(x, false);
            assert!(
                w <= -1.0 && (w * w.exp() - x).abs() < 1e-10,
                "W-1({x}) = {w}"
            );
        }
        assert!(lambert_w(-1.0, true).is_nan() && lambert_w(0.5, false).is_nan());
    }
}

#[cfg(test)]
mod hurwitz_tests {
    use super::hurwitz_zeta;

    #[test]
    fn hurwitz_values() {
        let pi = std::f64::consts::PI;
        assert!((hurwitz_zeta(2.0, 1.0) - pi * pi / 6.0).abs() < 1e-13);
        assert!((hurwitz_zeta(2.0, 0.5) - pi * pi / 2.0).abs() < 1e-12);
        // ζ(-1, a) = -B₂(a)/2 = -(a² - a + 1/6)/2
        let a = 0.3_f64;
        assert!((hurwitz_zeta(-1.0, a) + f64::midpoint(a * a - a, 1.0 / 6.0)).abs() < 1e-10);
        // Σ_{k=1}^{4} k^0.5 = ζ(-0.5, 1) - ζ(-0.5, 5)
        let direct: f64 = (1..=4).map(|k| f64::from(k).sqrt()).sum();
        assert!((hurwitz_zeta(-0.5, 1.0) - hurwitz_zeta(-0.5, 5.0) - direct).abs() < 1e-10);
    }
}

/// `E1(x) = ∫_x^∞ e^(-t)/t dt` for `x > 0`: the power series for
/// `x <= 1`, Lentz's continued fraction beyond.
#[must_use]
pub fn exp_integral_e1(x: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    }
    if x == f64::INFINITY {
        return 0.0;
    }
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x <= 1.0 {
        let (mut sum, mut term) = (0.0, 1.0);
        for k in 1..200 {
            let k = f64::from(k);
            term *= -x / k;
            let add = -term / k;
            sum += add;
            if add.abs() < 1e-17 * sum.abs() {
                break;
            }
        }
        return -EULER_GAMMA - x.ln() + sum;
    }
    // E1(x) = e^(-x) / (x + 1 - 1/(x + 3 - 4/(x + 5 - ...)))
    let tiny = 1e-300;
    let mut b = x + 1.0;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..500 {
        let i = f64::from(i);
        let a = -i * i;
        b += 2.0;
        d = 1.0 / a.mul_add(d, b);
        c = b + a / c;
        let delta = c * d;
        h *= delta;
        if (delta - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h * (-x).exp()
}

/// The exponential integral `Ei(x) = -PV ∫_(-x)^∞ e^(-t)/t dt`.
#[must_use]
pub fn exp_integral_ei(x: f64) -> f64 {
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == f64::NEG_INFINITY {
        return 0.0;
    }
    if x.is_nan() || x == 0.0 {
        return if x == 0.0 {
            f64::NEG_INFINITY
        } else {
            f64::NAN
        };
    }
    if x < 0.0 {
        return -exp_integral_e1(-x);
    }
    if x < 40.0 {
        let (mut sum, mut term) = (0.0, 1.0);
        for k in 1..400 {
            let k = f64::from(k);
            term *= x / k;
            let add = term / k;
            sum += add;
            if add < 1e-17 * sum {
                break;
            }
        }
        return EULER_GAMMA + x.ln() + sum;
    }
    // Asymptotic: e^x/x Σ k!/x^k.
    let (mut sum, mut term) = (1.0, 1.0);
    for k in 1..60 {
        let next = term * f64::from(k) / x;
        if next > term {
            break;
        }
        term = next;
        sum += term;
        if term < 1e-17 {
            break;
        }
    }
    x.exp() / x * sum
}

/// The logarithmic integral `li(x) = Ei(ln x)` for `x > 0`.
#[must_use]
pub fn log_integral(x: f64) -> f64 {
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    exp_integral_ei(x.ln())
}

/// `(Si(x), Ci(x))` for `x > 0`: power series for `x <= 2`, the continued
/// fraction of `E1(i x)` beyond (Numerical Recipes' `cisi`).
fn si_ci(x: f64) -> (f64, f64) {
    use num_complex::Complex64;
    if x <= 2.0 {
        let (mut si, mut ci) = (0.0, 0.0);
        let (mut odd, mut even) = (x, 1.0);
        for k in 0..60 {
            let n = 2 * k + 1;
            if k > 0 {
                odd *= -x * x / (f64::from(n - 1) * f64::from(n));
            }
            si += odd / f64::from(n);
            let m = 2 * k + 2;
            even *= -x * x / (f64::from(m - 1) * f64::from(m));
            ci += even / f64::from(m);
            if odd.abs() < 1e-18 && even.abs() < 1e-18 {
                break;
            }
        }
        return (si, EULER_GAMMA + x.ln() + ci);
    }
    let tiny = 1e-300;
    let mut b = Complex64::new(1.0, x);
    let mut c = Complex64::new(1.0 / tiny, 0.0);
    let mut d = Complex64::new(1.0, 0.0) / b;
    let mut h = d;
    for i in 1..500 {
        let i = f64::from(i);
        let a = -i * i;
        b += 2.0;
        d = Complex64::new(1.0, 0.0) / (d * a + b);
        c = b + Complex64::new(a, 0.0) / c;
        let delta = c * d;
        h *= delta;
        if (delta.re - 1.0).abs() + delta.im.abs() < 1e-16 {
            break;
        }
    }
    h *= Complex64::new(x.cos(), -x.sin());
    (std::f64::consts::FRAC_PI_2 + h.im, -h.re)
}

/// The sine integral `Si(x) = ∫_0^x sin(t)/t dt` (odd).
#[must_use]
pub fn sine_integral(x: f64) -> f64 {
    if x.is_infinite() {
        return x.signum() * std::f64::consts::FRAC_PI_2;
    }
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    let s = si_ci(x.abs()).0;
    if x < 0.0 { -s } else { s }
}

/// The cosine integral `Ci(x) = γ + ln x + ∫_0^x (cos t - 1)/t dt` for
/// `x > 0`.
#[must_use]
pub fn cosine_integral(x: f64) -> f64 {
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return 0.0;
    }
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    si_ci(x).1
}

/// The hyperbolic sine integral `Shi(x) = ∫_0^x sinh(t)/t dt` (odd).
#[must_use]
pub fn sinh_integral(x: f64) -> f64 {
    if x.is_infinite() {
        return x;
    }
    if x.is_nan() {
        return f64::NAN;
    }
    let a = x.abs();
    let s = if a < 1.0 {
        let (mut sum, mut term) = (0.0, a);
        for k in 0..40 {
            let n = 2 * k + 1;
            if k > 0 {
                term *= a * a / (f64::from(n - 1) * f64::from(n));
            }
            sum += term / f64::from(n);
        }
        sum
    } else {
        f64::midpoint(exp_integral_ei(a), exp_integral_e1(a))
    };
    if x < 0.0 { -s } else { s }
}

/// The hyperbolic cosine integral `Chi(x) = γ + ln x + ∫_0^x (cosh t - 1)/t dt`
/// for `x > 0`.
#[must_use]
pub fn cosh_integral(x: f64) -> f64 {
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    (exp_integral_ei(x) - exp_integral_e1(x)) / 2.0
}

#[cfg(test)]
mod exponential_integral_tests {
    use super::*;

    #[test]
    fn reference_values() {
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12 * (1.0 + b.abs());
        assert!(close(exp_integral_ei(1.0), 1.895_117_816_355_936_8));
        assert!(close(exp_integral_ei(-1.0), -0.219_383_934_395_520_27));
        assert!(close(exp_integral_ei(50.0), 1.058_563_689_713_169e20));
        assert!(close(exp_integral_e1(2.5), 0.024_914_917_870_269_73));
        assert!(close(sine_integral(1.0), 0.946_083_070_367_183));
        assert!(close(sine_integral(10.0), 1.658_347_594_218_874));
        assert!(close(cosine_integral(1.0), 0.337_403_922_900_968_1));
        assert!(close(cosine_integral(10.0), -0.045_456_433_004_455_37));
        assert!(close(log_integral(10.0), 6.165_599_504_787_297));
        assert!(close(sinh_integral(1.0), 1.057_250_875_375_728_5));
        assert!(close(cosh_integral(1.0), 0.837_866_940_980_208_2));
        assert!(close(sinh_integral(3.0), 4.973_440_475_859_807));
    }
}

/// `Ai(0)`.
const AIRY_AI0: f64 = 0.355_028_053_887_817_2;
/// `-Ai'(0)`.
const AIRY_AIP0: f64 = 0.258_819_403_792_806_8;

/// `K_ν(z) = ∫_0^∞ e^(-z cosh t) cosh(ν t) dt` for `z > 0`, by the
/// trapezoidal rule (exponentially convergent for this integrand).
fn bessel_k_fractional(
    nu: f64,
    z: f64,
) -> f64 {
    let h = 0.05;
    let mut sum = 0.5 * (-z).exp();
    let mut k = 1.0_f64;
    loop {
        let t = k * h;
        let term = (-z * t.cosh()).exp() * (nu * t).cosh();
        sum += term;
        if term < 1e-18 * sum || k > 4000.0 {
            break;
        }
        k += 1.0;
    }
    sum * h
}

/// `(y, y')` at `x` of the solution of `y'' = x y` with `y(0) = y0`,
/// `y'(0) = dy0`, by Taylor steps (the Airy equation's coefficients follow
/// `a_(n+2) = (x0 a_n + a_(n-1)) / ((n+2)(n+1))` at any centre `x0`).
fn airy_taylor(
    x: f64,
    y0: f64,
    dy0: f64,
) -> (f64, f64) {
    let steps = (x.abs() / 0.5).ceil().max(1.0);
    let h = x / steps;
    let (mut y, mut dy, mut x0) = (y0, dy0, 0.0_f64);
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let steps_count = steps as usize;
    for _ in 0..steps_count {
        let mut a = [0.0_f64; 64];
        a[0] = y;
        a[1] = dy;
        a[2] = x0 * y / 2.0;
        for n in 1..62 {
            a[n + 2] = x0.mul_add(a[n], a[n - 1]) / (((n + 2) * (n + 1)) as f64);
        }
        let (mut value, mut slope, mut power) = (0.0, 0.0, 1.0);
        for (n, &c) in a.iter().enumerate() {
            value += c * power;
            if n + 1 < a.len() {
                slope += (n + 1) as f64 * a[n + 1] * power;
            }
            power *= h;
        }
        y = value;
        dy = slope;
        x0 += h;
    }
    (y, dy)
}

/// Smallest `-x` from which the oscillatory asymptotic expansions of the Airy functions are used
/// (`zeta = 2/3 |x|^(3/2) >= 21`; the optimal truncation error is about `e^(-2 zeta)`).
const AIRY_ASYMPTOTIC_FROM: f64 = 10.0;

/// `x` above which `Bi(x)` overflows (`exp(2/3 x^(3/2)) > 1.8e308`).
const AIRY_BI_OVERFLOW: f64 = 105.0;

/// `(Ai, Ai', Bi, Bi')` at `x <= -10` from the expansions of DLMF 9.7.9-9.7.12.
fn airy_negative_asymptotic(x: f64) -> (f64, f64, f64, f64) {
    let z = -x;
    let zeta = 2.0 / 3.0 * z.powf(1.5);
    let phase = zeta - std::f64::consts::FRAC_PI_4;
    let (sn, cs) = phase.sin_cos();
    // u_k and v_k = -(6k+1)/(6k-1) u_k; the series are asymptotic: stop at the smallest term
    let (mut p, mut q, mut pd, mut qd) = (1.0_f64, 0.0_f64, 1.0_f64, 0.0_f64);
    let mut u = 1.0_f64;
    let mut last = f64::INFINITY;
    for k in 1..=60_u32 {
        let kf = f64::from(k);
        u *= (6.0 * kf - 5.0) * (6.0 * kf - 3.0) * (6.0 * kf - 1.0)
            / ((2.0 * kf - 1.0) * 216.0 * kf);
        let v = -(6.0 * kf + 1.0) / (6.0 * kf - 1.0) * u;
        let t = u / zeta.powi(k as i32);
        if t.abs() >= last {
            break;
        }
        last = t.abs();
        let sign = if (k / 2) % 2 == 0 {
            1.0
        } else {
            -1.0
        };
        if k % 2 == 0 {
            p += sign * t;
            pd += sign * v / zeta.powi(k as i32);
        } else {
            q += sign * t;
            qd += sign * v / zeta.powi(k as i32);
        }
        if t.abs() < 1e-17 {
            break;
        }
    }
    let amp = 1.0 / (std::f64::consts::PI.sqrt() * z.powf(0.25));
    let damp = z.powf(0.25) / std::f64::consts::PI.sqrt();
    let ai = amp * (cs * p + sn * q);
    let bi = amp * (-sn * p + cs * q);
    let aip = damp * (sn * pd - cs * qd);
    let bip = damp * (cs * pd + sn * qd);
    (ai, aip, bi, bip)
}

/// The Airy function `Ai(x)`.
///
/// `Ai(x) -> 0` as `x -> +-inf`; for `x <= -10` the oscillatory asymptotic expansion is used
/// (the Taylor continuation from the origin is only used on `(-10, 0]`).
#[must_use]
pub fn airy_ai(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return 0.0;
    }
    if x > 0.0 {
        let zeta = 2.0 / 3.0 * x.powf(1.5);
        return (x / 3.0).sqrt() * bessel_k_fractional(1.0 / 3.0, zeta) / std::f64::consts::PI;
    }
    if x <= -AIRY_ASYMPTOTIC_FROM {
        return airy_negative_asymptotic(x).0;
    }
    airy_taylor(x, AIRY_AI0, -AIRY_AIP0).0
}

/// `Ai'(x)` (no limit at `-inf`: NaN there).
#[must_use]
pub fn airy_ai_prime(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return 0.0;
    }
    if x > 0.0 {
        let zeta = 2.0 / 3.0 * x.powf(1.5);
        return -x / (std::f64::consts::PI * 3.0_f64.sqrt()) * bessel_k_fractional(2.0 / 3.0, zeta);
    }
    if x <= -AIRY_ASYMPTOTIC_FROM {
        return airy_negative_asymptotic(x).1;
    }
    airy_taylor(x, AIRY_AI0, -AIRY_AIP0).1
}

/// The Airy function `Bi(x)`: `+inf` from `x > 105` (overflow) and at `+inf`, `0` at `-inf`.
#[must_use]
pub fn airy_bi(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x > AIRY_BI_OVERFLOW {
        return f64::INFINITY;
    }
    if x == f64::NEG_INFINITY {
        return 0.0;
    }
    if x <= -AIRY_ASYMPTOTIC_FROM {
        return airy_negative_asymptotic(x).2;
    }
    let r3 = 3.0_f64.sqrt();
    airy_taylor(x, r3 * AIRY_AI0, r3 * AIRY_AIP0).0
}

/// `Bi'(x)` (`+inf` from `x > 105`; no limit at `-inf`: NaN there).
#[must_use]
pub fn airy_bi_prime(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return f64::NAN;
    }
    if x > AIRY_BI_OVERFLOW {
        return f64::INFINITY;
    }
    if x <= -AIRY_ASYMPTOTIC_FROM {
        return airy_negative_asymptotic(x).3;
    }
    let r3 = 3.0_f64.sqrt();
    airy_taylor(x, r3 * AIRY_AI0, r3 * AIRY_AIP0).1
}

#[cfg(test)]
mod airy_tests {
    use super::*;

    #[test]
    fn reference_values() {
        let close =
            |a: f64, b: f64| (a - b).abs() < 1e-11 * (1.0 + b.abs()) || (a / b - 1.0).abs() < 1e-10;
        assert!(close(airy_ai(1.0), 0.135_292_416_312_881_4));
        assert!(close(airy_ai(5.0), 1.083_444_281_360_744e-4));
        assert!(close(airy_ai(-5.0), 0.350_761_009_024_114));
        assert!(close(airy_bi(1.0), 1.207_423_594_952_871));
        assert!(close(airy_bi(-5.0), -0.138_369_134_901_600_5));
        assert!(close(airy_ai_prime(1.0), -0.159_147_441_296_793_2));
        assert!(close(airy_bi_prime(0.0), 0.448_288_357_353_826_4));
        // The Wronskian Ai Bi' - Ai' Bi = 1/π.
        for x in [-7.3, -1.2, 0.4, 2.9, 6.1] {
            let w = airy_ai(x) * airy_bi_prime(x) - airy_ai_prime(x) * airy_bi(x);
            assert!((w - 1.0 / std::f64::consts::PI).abs() < 1e-10, "{x}: {w}");
        }
    }
}

#[cfg(test)]
mod debye_tests {
    use super::*;

    fn rel_err(
        a: f64,
        b: f64,
    ) -> f64 {
        (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
    }

    #[test]
    fn uniform_expansion_matches_reference_values() {
        // mpmath 50-digit values at ν = 100 (the old quadrature agrees to 1e-14).
        let cases = [
            (90.0, 2.602_130_581_996_329e-3, -2.830_777_138_718_563),
            (95.0, 2.315_076_800_942_797e-2, -4.576_220_049_590_498e-1),
            (99.0, 7.768_716_170_045_94e-2, -2.010_721_995_738_357e-1),
            (100.0, 9.636_667_329_586_116e-2, -1.669_214_114_175_765e-1),
            (101.0, 1.148_013_214_278_991e-1, -1.332_273_838_156_164e-1),
            (105.0, 1.358_350_278_036_409e-1, 2.633_762_899_822_202e-2),
            (110.0, -5.385_144_819_503_075e-2, 1.043_543_190_047_037e-1),
        ];
        for (x, j, y) in cases {
            assert!(rel_err(bessel_j(100.0, x), j) < 5e-13, "J(100, {x})");
            assert!(rel_err(bessel_y(100.0, x), y) < 5e-13, "Y(100, {x})");
        }
    }

    #[test]
    fn recurrence_holds_across_turning_point() {
        let nu = 100.0;
        for x in [76.0, 85.0, 95.0, 99.0, 100.0, 101.0, 105.0, 115.0, 124.0] {
            let jm = bessel_j(nu - 1.0, x);
            let j0 = bessel_j(nu, x);
            let jp = bessel_j(nu + 1.0, x);
            let lhs = jm + jp;
            let rhs = 2.0 * nu / x * j0;
            assert!(
                rel_err(lhs, rhs) < 5e-13,
                "J recurrence at {x}: {lhs} vs {rhs}"
            );

            let ym = bessel_y(nu - 1.0, x);
            let y0 = bessel_y(nu, x);
            let yp = bessel_y(nu + 1.0, x);
            let lhs = ym + yp;
            let rhs = 2.0 * nu / x * y0;
            assert!(
                rel_err(lhs, rhs) < 5e-13,
                "Y recurrence at {x}: {lhs} vs {rhs}"
            );
        }
    }

    #[test]
    fn wronskian_holds_across_turning_point() {
        let nu = 100.0;
        for x in [76.0, 85.0, 95.0, 99.0, 100.0, 101.0, 105.0, 115.0, 124.0] {
            // J' = (J_{ν−1} − J_{ν+1})/2, likewise for Y (DLMF 10.6.1).
            let j = bessel_j(nu, x);
            let y = bessel_y(nu, x);
            let jp = (bessel_j(nu - 1.0, x) - bessel_j(nu + 1.0, x)) / 2.0;
            let yp = (bessel_y(nu - 1.0, x) - bessel_y(nu + 1.0, x)) / 2.0;
            let w = j * yp - jp * y;
            let expected = 2.0 / (std::f64::consts::PI * x);
            assert!(
                rel_err(w, expected) < 1e-12,
                "Wronskian at {x}: {w} vs {expected}"
            );
        }
    }

    #[test]
    fn sweep_through_turning_point_is_smooth() {
        let nu = 100.0;
        let mut worst = 0.0_f64;
        let mut x = 90.0;
        while x <= 110.0 {
            let jm = bessel_j(nu - 1.0, x);
            let j0 = bessel_j(nu, x);
            let jp = bessel_j(nu + 1.0, x);
            assert!(j0.is_finite(), "J(100, {x}) is not finite");
            let residual = rel_err(jm + jp, 2.0 * nu / x * j0);
            worst = worst.max(residual);
            assert!(residual < 2e-12, "recurrence residual {residual} at {x}");
            x += 0.1;
        }
        assert!(worst < 2e-12, "worst recurrence residual {worst}");
    }

    #[test]
    fn seam_between_uniform_and_quadrature_is_continuous() {
        // The uniform band is |x/ν − 1| <= 1/4, so 75 and 125 are the seams.
        // Both sides are checked against mpmath, which is what rules out a jump.
        let cases = [
            // Just outside the band: quadrature, accurate to ~2e-9 relative on
            // this exponentially small J (absolute error ~1e-16).
            (
                74.999,
                7.441_273_517_697_171e-8,
                -6.469_552_632_950_964e4,
                1e-8,
            ),
            // On the seam the uniform path is used; it is accurate to 1e-12.
            (
                75.0,
                7.447_900_590_590_388e-8,
                -6.463_907_226_123_209e4,
                1e-12,
            ),
            (
                125.0,
                -8.359_856_456_811_825e-2,
                -3.865_789_123_650_626e-2,
                1e-12,
            ),
            (
                125.001,
                -8.357_441_484_005_472e-2,
                -3.870_764_520_013_322e-2,
                1e-8,
            ),
        ];
        for (x, j, y, tol) in cases {
            assert!(rel_err(bessel_j(100.0, x), j) < tol, "J seam at {x}");
            assert!(rel_err(bessel_y(100.0, x), y) < tol, "Y seam at {x}");
        }
    }

    #[test]
    fn nicholson_integral_matches_jy() {
        // J_ν(x)² + Y_ν(x)² = (8/π²) ∫_0^∞ K_0(2x sinh t) cosh(2νt) dt
        // (Watson, §13.72).  The integrand is logarithmically singular at t = 0
        // and decays like exp(2νt − 2x sinh t), so 2 is a safe upper limit here.
        let nu = 100.0_f64;
        let x = 100.0_f64;
        let integral = quad(
            |t| bessel_k_fractional(0.0, 2.0 * x * t.sinh()) * (2.0 * nu * t).cosh(),
            0.0,
            2.0,
        );
        let rhs = 8.0 / (std::f64::consts::PI * std::f64::consts::PI) * integral;
        let lhs = bessel_j(nu, x).powi(2) + bessel_y(nu, x).powi(2);
        assert!(rel_err(lhs, rhs) < 1e-12, "{lhs} vs {rhs}");
    }

    #[test]
    fn large_order_limit_is_the_airy_value() {
        // DLMF 10.20.4 at z = 1: ζ = 0, prefactor (4ζ/(1−z²))^{1/4} -> 2^{1/3},
        // A_0(0) = 1, so J_ν(ν) ~ 2^{1/3} Ai(0) ν^{-1/3}.  The B_0 term is
        // O(ν^{-5/3}) with B_0(0) = 2^{1/3}/70.
        let nu = 1.0e6_f64;
        let constant = 2.0_f64.cbrt() * airy_ai(0.0);
        let scaled = bessel_j(nu, nu) * nu.cbrt();
        assert!(rel_err(scaled, constant) < 1e-6, "{scaled} vs {constant}");
    }

    #[test]
    fn outside_band_and_small_orders_are_bit_identical() {
        // Values produced by the quadrature branch before the uniform path was
        // added; the guard must leave them untouched.
        assert_eq!(
            bessel_j(100.0, 74.0).to_bits(),
            3.017_594_391_878_188_9e-8_f64.to_bits()
        );
        assert_eq!(
            bessel_y(100.0, 74.0).to_bits(),
            (-1.568_827_309_478_636_9e5_f64).to_bits()
        );
        assert_eq!(
            bessel_j(100.0, 126.0).to_bits(),
            (-4.630_129_654_029_885_3e-2_f64).to_bits()
        );
        assert_eq!(
            bessel_y(100.0, 126.0).to_bits(),
            (-7.846_463_206_197_208_7e-2_f64).to_bits()
        );
        assert_eq!(
            bessel_j(2.5, 5.0).to_bits(),
            2.403_772_011_113_173_3e-1_f64.to_bits()
        );
        assert_eq!(
            bessel_y(2.5, 5.0).to_bits(),
            2.943_723_749_617_924_2e-1_f64.to_bits()
        );
        assert_eq!(
            bessel_j(10.0, 20.0).to_bits(),
            1.864_825_580_239_450_3e-1_f64.to_bits()
        );
    }
}
