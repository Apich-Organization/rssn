//! Complex gamma, reciprocal gamma, polygamma, beta, incomplete gamma/beta,
//! Pochhammer, binomial, and Barnes G functions.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::hypergeometric::hyp2f1_cplx;
use crate::kernels::special_complex::digamma_complex;
use crate::kernels::special_complex::gamma_complex;
use crate::kernels::special_complex::ln_gamma_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

fn is_nonpos_int(z: Complex64) -> bool {
    z.im == 0.0 && z.re <= 0.0 && (z.re - z.re.round()).abs() == 0.0
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
        let s = (Complex64::new(PI, 0.0) * z).sin();
        return s * gamma_complex(ONE - z) / PI;
    }
    let g = gamma_complex(z);
    if g.re.is_infinite() || g.im.is_infinite() {
        return ZERO;
    }
    ONE / g
}

/// Complex polygamma function `psi^(n)(z) = d^n/dz^n digamma(z)`.
#[must_use]
pub fn polygamma_complex(
    n: u32,
    mut z: Complex64,
) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if n == 0 {
        return digamma_complex(z);
    }
    if is_nonpos_int(z) {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    // Reflection for Re(z) < 0.5: psi^(n)(1 - z) + (-1)^(n+1) psi^(n)(z) = (-1)^n d^n/dz^n (pi cot(pi z))
    // Use finite-step recurrence on z + m for small Re(z), and for Re(z) < 0 use reflection for n=1,2
    // or shift z -> 1-z for n=1: psi'(z) + psi'(1-z) = pi^2 / sin^2(pi z).
    if z.re < 0.0 && n == 1 {
        let sin_pz = (PI * z).sin();
        let refl = Complex64::new(PI * PI, 0.0) / (sin_pz * sin_pz);
        return refl - polygamma_complex(1, ONE - z);
    }
    if z.re < 0.0 && n == 2 {
        let pz = PI * z;
        let sin_pz = pz.sin();
        let cos_pz = pz.cos();
        let refl = Complex64::new(-2.0 * PI * PI * PI, 0.0) * cos_pz / (sin_pz * sin_pz * sin_pz);
        return polygamma_complex(2, ONE - z) + refl;
    }

    let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
    let mut fact = 1.0_f64;
    for k in 1..=n {
        fact *= f64::from(k);
    }

    // Shift z right until |z| >= 10 and Re(z) >= 10.
    let mut sum = ZERO;
    while z.re < 10.0 && z.norm() < 12.0 {
        sum += ONE / z.powi(i32::try_from(n + 1).unwrap_or(i32::MAX));
        z += ONE;
    }

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
            ratio = (nf + 1.0) / 2.0;
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
    gamma_complex(a) * gamma_complex(b) * rgamma_complex(a + b)
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
    gamma_complex(z + w) * rgamma_complex(z)
}

/// Complex binomial coefficient `C(n, k) = Gamma(n + 1) / (Gamma(k + 1) * Gamma(n - k + 1))`.
#[must_use]
pub fn binomial_complex(
    n: Complex64,
    k: Complex64,
) -> Complex64 {
    gamma_complex(n + ONE) * rgamma_complex(k + ONE) * rgamma_complex(n - k + ONE)
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
    if z.norm() <= a.norm() + 2.0 || z.re < 0.0 && z.norm() < 15.0 {
        // Power series: z^a * exp(-z) * sum_{k>=0} z^k / (a * (a+1) * ... * (a+k))
        let mut term = ONE / a;
        let mut sum = term;
        for k in 1..400 {
            term *= z / (a + (k as f64));
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        z.powc(a) * (-z).exp() * sum
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
    if z.norm() <= a.norm() + 2.0 || (z.re < 0.0 && z.im.abs() < 1.0 && z.norm() < 15.0) {
        return gamma_complex(a) - lower_incomplete_gamma_complex(a, z);
    }
    // Legendre continued fraction via modified Lentz's method:
    // Gamma(a, z) = z^a * exp(-z) * 1 / (z + (1 - a) / (1 + 1 / (z + (2 - a) / (1 + 2 / ...))))
    // Equivalent form: b_0 = z + 1 - a, a_k = k * (a - k), b_k = z + 2k + 1 - a
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
        let kf = k as f64;
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
    z.powc(a) * (-z).exp() * h
}

/// Regularized lower incomplete gamma `P(a, z) = gamma(a, z) / Gamma(a)`.
#[must_use]
pub fn reg_lower_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
    lower_incomplete_gamma_complex(a, z) * rgamma_complex(a)
}

/// Regularized upper incomplete gamma `Q(a, z) = Gamma(a, z) / Gamma(a)`.
#[must_use]
pub fn reg_upper_gamma_complex(
    a: Complex64,
    z: Complex64,
) -> Complex64 {
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

/// Complex Barnes G-function `G(z)` satisfying `G(z + 1) = Gamma(z) * G(z)`, `G(1) = 1`.
#[must_use]
pub fn barnes_g_complex(mut z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if is_nonpos_int(z) {
        return ZERO;
    }
    if z == ONE || z == Complex64::new(2.0, 0.0) {
        return ONE;
    }
    // Shift z right until Re(z) >= 8 using G(z) = G(z + 1) / Gamma(z).
    let mut factor = ONE;
    while z.re < 8.0 {
        factor /= gamma_complex(z);
        z += ONE;
    }
    // Asymptotic for ln G(z + 1):
    // ln G(z + 1) ~ z^2/2 * ln(z) - 3/4 * z^2 + z/2 * ln(2*pi) - 1/12 * ln(z) + zeta'(-1)
    //             + sum_{k>=1} B_{2k+2} / (2k * (2k+1) * (2k+2) * z^(2k))
    // where zeta'(-1) = 1/12 - ln(A), Glaisher-Kinkelin A = 1.28242712910062263687534256886979.
    const ZETA_PRIME_MINUS_1: f64 = -0.165_421_143_700_450_92;
    const LN_2PI: f64 = 1.837_877_066_409_345_3;
    let w = z - ONE; // G(z) = G(w + 1)
    let ln_w = w.ln();
    let w2 = w * w;
    let mut ln_g =
        0.5 * w2 * ln_w - 0.75 * w2 + 0.5 * w * LN_2PI - (1.0 / 12.0) * ln_w + ZETA_PRIME_MINUS_1;
    let inv_w2 = ONE / w2;
    let mut term_pow = inv_w2;
    // k=1: B_4 / (2*3*4) = (-1/30)/24 = -1/720
    ln_g += (-1.0 / 720.0) * term_pow;
    term_pow *= inv_w2;
    // k=2: B_6 / (4*5*6) = (1/42)/120 = 1/5040
    ln_g += (1.0 / 5040.0) * term_pow;
    term_pow *= inv_w2;
    // k=3: B_8 / (6*7*8) = (-1/30)/336 = -1/10080
    ln_g += (-1.0 / 10080.0) * term_pow;
    ln_g.exp() * factor
}
