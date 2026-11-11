//! # Special functions of a complex argument
//!
//! The real-argument routines of [`special`](super::special) cover the
//! gamma family, Bessel functions and Lambert W on the real line. This module
//! continues the most used of them into the complex plane:
//!
//! * [`ln_gamma_complex`], [`gamma_complex`], [`digamma_complex`]: Stirling's
//!   series after shifting the argument to `Re z >= 12` by the recurrence,
//!   reflection for `Re z < 0.5` (about `1e-14` relative, away from poles).
//!   `ln_gamma_complex` is the analytic branch with
//!   `ln Γ(z + 1) = ln Γ(z) + ln z` (the one SciPy calls `loggamma`),
//!   continuous off the non-positive real axis.
//! * [`lambert_w_complex`]: every branch `k`, by Halley iteration from the
//!   starting values of Corless, Gonnet, Hare, Jeffrey and Knuth.
//! * [`bessel_j_complex`], [`bessel_i_complex`]: real order `ν`, complex `z`;
//!   the ascending series for moderate `|z|` and Hankel's asymptotic
//!   expansion beyond it. Accuracy degrades (about `1e-9`) for orders with
//!   `ν² > 2 |z|` when `|z|` is near the series/asymptotic transition.
//!
//! No function panics; invalid input gives NaN components, poles give
//! infinity.

use std::f64::consts::PI;

use num_complex::Complex64;

const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const INF: Complex64 = Complex64::new(f64::INFINITY, 0.0);
const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

fn is_nan(z: Complex64) -> bool {
    z.re.is_nan() || z.im.is_nan()
}

fn is_inf(z: Complex64) -> bool {
    z.re.is_infinite() || z.im.is_infinite()
}

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

/// `ln Γ(z)` for `Re z >= 12`.
fn stirling(z: Complex64) -> Complex64 {
    let zi = 1.0 / z;
    let z2 = zi * zi;
    let mut term = zi;
    let mut sum = ZERO;
    for c in STIRLING {
        sum += term * c;
        term *= z2;
    }
    (z - 0.5) * z.ln() - z + 0.5 * (2.0 * PI).ln() + sum
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
    if z.re.abs() > 1e7 {
        // Too far to shift by the recurrence; Stirling alone is accurate for
        // a large positive real part, the negative side is unsupported.
        return if z.re > 0.0 {
            stirling(z)
        } else {
            NAN
        };
    }
    // Shift to Re >= 12, accumulating principal logarithms.
    let mut acc = ZERO;
    let mut w = z;
    while w.re < 12.0 {
        acc += w.ln();
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
    if z.re < 0.5 {
        // Γ(z) = π / (sin(πz) Γ(1 - z))
        let s = (PI * z).sin();
        if s.re == 0.0 && s.im == 0.0 {
            return INF;
        }
        return PI / (s * gamma_complex(1.0 - z));
    }
    ln_gamma_complex(z).exp()
}

/// The digamma function `ψ(z) = Γ'(z) / Γ(z)`; poles give NaN.
#[must_use]
pub fn digamma_complex(z: Complex64) -> Complex64 {
    if is_nan(z) || is_inf(z) || is_pole(z) {
        return NAN;
    }
    if z.re < 0.0 {
        // ψ(z) = ψ(1 - z) - π cot(π z)
        let pz = PI * z;
        return digamma_complex(1.0 - z) - PI * pz.cos() / pz.sin();
    }
    let mut acc = ZERO;
    let mut w = z;
    while w.re < 12.0 {
        acc += 1.0 / w;
        w += 1.0;
    }
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
    let wi = 1.0 / w;
    let w2 = wi * wi;
    let mut term = w2;
    let mut sum = ZERO;
    for c in B {
        sum += term * c;
        term *= w2;
    }
    w.ln() - 0.5 * wi - sum - acc
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
    let e = std::f64::consts::E;
    let near_branch = (z + 1.0 / e).norm() < 0.3;
    let mut w = if near_branch && (k == 0 || (k == -1 && z.im >= 0.0) || (k == 1 && z.im < 0.0)) {
        let mut p = (2.0 * (e * z + 1.0)).sqrt();
        if k != 0 {
            p = -p;
        }
        -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p
    } else if k == 0 && z.norm() < 0.3 {
        z - z * z + 1.5 * z * z * z
    } else if k == 0 && z.norm() < 2.5 {
        // ln(1 + z) is a safe start away from the branch point.
        (1.0 + z).ln()
    } else {
        let l = z.ln() + Complex64::new(0.0, 2.0 * PI * f64::from(k));
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

/// Hankel's asymptotic expansion of `H1_ν(z)` and `H2_ν(z)` (DLMF 10.17.5),
/// for large `|z|` with `|arg z| < π`.
fn hankel(
    nu: f64,
    z: Complex64,
) -> (Complex64, Complex64) {
    let mu = 4.0 * nu * nu;
    let one = Complex64::new(1.0, 0.0);
    let (mut p, mut q) = (one, one);
    let (mut term_p, mut term_q) = (one, one);
    let mut last = f64::INFINITY;
    for k in 1..80 {
        let kf = f64::from(k);
        let factor = (mu - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf * z);
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
    let chi = z - (0.5 * nu + 0.25) * PI;
    let pref = (2.0 / (PI * z)).sqrt();
    (pref * (I * chi).exp() * p, pref * (-I * chi).exp() * q)
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

fn bessel_invalid(
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
    z.norm() >= 14.0 && z.norm() > nu * nu / 2.0
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
    if let Some(n) = integer_order(nu) {
        if n < 0 {
            // J_{-n} = (-1)^n J_n
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            return sign * bessel_j_complex(nu.abs(), z);
        }
    }
    if z.re == 0.0 && z.im == 0.0 {
        return bessel_at_origin(nu);
    }
    if use_asymptotic(nu, z) {
        if z.im == 0.0 && z.re < 0.0 {
            // J_ν(-x) = e^{iπν} J_ν(x) on the upper lip of the cut.
            return (I * PI * nu).exp() * bessel_j_complex(nu, -z);
        }
        let (h1, h2) = hankel(nu, z);
        return 0.5 * (h1 + h2);
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
    if let Some(n) = integer_order(nu) {
        if n < 0 {
            // I_{-n} = I_n
            return bessel_i_complex(nu.abs(), z);
        }
    }
    if z.re == 0.0 && z.im == 0.0 {
        return bessel_at_origin(nu);
    }
    if use_asymptotic(nu, z) {
        // DLMF 10.27.6: I_ν(z) = e^{∓νπi/2} J_ν(z e^{±πi/2}),
        // upper sign for -π <= arg z <= π/2.
        let rot = if z.im > 0.0 && z.re < 0.0 {
            -I
        } else {
            I
        };
        return (-rot * 0.5 * nu * PI).exp() * bessel_j_complex(nu, rot * z);
    }
    bessel_series(nu, z, false)
}
