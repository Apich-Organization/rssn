//! Complex Bessel Y/K, Hankel H1/H2, spherical Bessel, Airy Ai/Bi (and
//! derivatives), Kelvin, Struve, Anger, and Weber functions.

use std::f64::consts::FRAC_1_SQRT_2;
use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_segment_checked;
use crate::kernels::integrate::gauss_kronrod_any;
use crate::kernels::special_complex::bessel_i_complex;
use crate::kernels::special_complex::bessel_j_asymptotic_complex_order;
use crate::kernels::special_complex::bessel_j_complex;
use crate::kernels::special_complex::gamma_complex;
use crate::kernels::special_complex::ln_gamma_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);
const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
// The integer series initializes its finite sum with (n - 1)! in f64.
const MAX_FACTORIAL_ORDER: f64 = 171.0;
use std::f64::consts::EULER_GAMMA;

/// Complex Bessel function of the second kind `Y_nu(z)` for real order `nu`.
/// Returns NaN if the integer series would start with an overflowing factorial.
#[must_use]
pub fn bessel_y_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO || !nu.is_finite() {
        return NAN;
    }
    let nu_round = nu.round();
    if (nu - nu_round).abs() > 1e-10 {
        let s = (nu * PI).sin();
        let c = (nu * PI).cos();
        return (bessel_j_complex(nu, z) * c - bessel_j_complex(-nu, z)) / s;
    }
    let order = nu_round.abs();
    let sign = if nu_round < 0.0 && nu_round % 2.0 != 0.0 {
        -1.0
    } else {
        1.0
    };

    // For large |z|, use Hankel asymptotic expansion for Y_n(z).
    if z.norm() > 25.0 + 1.5 * order {
        let chi = z - (0.5 * order + 0.25) * PI;
        let mu = 4.0 * order * order;
        let eight_z = 8.0 * z;
        let mut p = ONE;
        let mut q = ZERO;
        let mut term = ONE;
        for k in 1..=12 {
            let kf = f64::from(k);
            let num = mu - (2.0 * kf - 1.0) * (2.0 * kf - 1.0);
            term *= num / (kf * eight_z);
            let s = if (k / 2) % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            if k % 2 == 1 {
                q += s * term;
            } else {
                p += s * term;
            }
            if term.norm() < 1e-16 {
                break;
            }
        }
        let scale = (2.0 / (PI * z)).sqrt();
        return sign * scale * (p * chi.sin() + q * chi.cos());
    }

    // DLMF 10.8.1 ascending series for integer n >= 0:
    if order > MAX_FACTORIAL_ORDER {
        return NAN;
    }
    let n_abs = order as u32;
    let half_z = 0.5 * z;
    let half_z2 = half_z * half_z;
    let j_n = bessel_j_complex(f64::from(n_abs), z);
    let mut sum1 = ZERO;
    if n_abs > 0 {
        let mut fact_n_m_1 = 1.0_f64;
        for k in 1..n_abs {
            fact_n_m_1 *= f64::from(k);
        }
        let mut term = Complex64::new(fact_n_m_1, 0.0);
        for k in 0..n_abs {
            if k > 0 {
                let kf = f64::from(k);
                let nm_k = f64::from(n_abs - k);
                term = term * half_z2 / (kf * nm_k);
            }
            sum1 += term;
        }
        sum1 *= -half_z.powi(-(n_abs as i32));
    }

    // Second sum: (-1)^k * (psi(k+1) + psi(n+k+1)) * (z/2)^(2k) / (k! * (n+k)!)
    let mut h_k = -EULER_GAMMA;
    let mut h_nk = -EULER_GAMMA;
    let mut fact_nk = 1.0_f64;
    for j in 1..=n_abs {
        h_nk += 1.0 / f64::from(j);
        fact_nk *= f64::from(j);
    }
    let mut term2 = ONE / fact_nk;
    let mut sum2 = term2 * (h_k + h_nk);
    for k in 1..80 {
        let kf = f64::from(k);
        h_k += 1.0 / kf;
        h_nk += 1.0 / (f64::from(n_abs) + kf);
        term2 *= -half_z2 / (kf * (f64::from(n_abs) + kf));
        let add = term2 * (h_k + h_nk);
        sum2 += add;
        if add.norm() <= 1e-16 * (1.0 + sum2.norm()) {
            break;
        }
    }
    sum2 *= -half_z.powi(n_abs as i32);

    let y_n = (2.0 / PI) * (half_z.ln() * j_n) + (sum1 + sum2) / PI;
    sign * y_n
}

/// Complex modified Bessel function of the second kind `K_nu(z)` for real order `nu`.
/// Returns NaN if the integer series would start with an overflowing factorial.
#[must_use]
pub fn bessel_k_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO || !nu.is_finite() {
        return NAN;
    }
    let nu_abs = nu.abs();
    // For large |z|, use Hankel asymptotic expansion for K_nu(z) = sqrt(pi/(2z)) e^(-z) sum ...
    if z.norm() > 12.0 + 1.5 * nu_abs {
        if z.arg().abs() >= PI - 0.15 {
            // Analytic continuation across (-inf, 0]: K_nu(z e^{+-i pi}) = e^{-+i pi nu} K_nu(z) -+ i pi I_nu(z)
            let (z_refl, sgn) = if z.im.is_sign_positive() {
                (-z, 1.0)
            } else {
                (-z, -1.0)
            };
            let phase = Complex64::from_polar(1.0, -sgn * nu_abs * PI);
            return phase * bessel_k_complex(nu_abs, z_refl)
                - sgn * I * PI * bessel_i_complex(nu_abs, z_refl);
        }
        let mu = 4.0 * nu_abs * nu_abs;
        let eight_z = 8.0 * z;
        let mut sum = ONE;
        let mut term = ONE;
        for k in 1..=25 {
            let kf = f64::from(k);
            let next_term = term * (mu - (2.0 * kf - 1.0) * (2.0 * kf - 1.0)) / (kf * eight_z);
            if next_term.norm() > term.norm() && k > 2 {
                break;
            }
            term = next_term;
            sum += term;
            if term.norm() < 1e-16 * sum.norm() {
                break;
            }
        }
        return (PI / (2.0 * z)).sqrt() * (-z).exp() * sum;
    }
    // Connection with Hankel H_nu^(1)(i * z):
    // K_nu(z) = (pi/2) * i^(nu+1) * H_nu^(1)(i*z) when Re(z) >= 0 or Im(z) < 0;
    // or for non-integer nu: K_nu(z) = (pi/2) * (I_{-nu}(z) - I_nu(z)) / sin(nu*pi).
    let nu_round = nu_abs.round();
    if (nu_abs - nu_round).abs() > 1e-8 {
        let s = (nu_abs * PI).sin();
        return (0.5 * PI / s) * (bessel_i_complex(-nu_abs, z) - bessel_i_complex(nu_abs, z));
    }
    // Integer n >= 0: use small-order series (DLMF 10.31.1) for |z| <= 25.
    if nu_round > MAX_FACTORIAL_ORDER {
        return NAN;
    }
    let n = nu_round as u32;
    let half_z = 0.5 * z;
    let half_z2 = half_z * half_z;
    let i_n = bessel_i_complex(f64::from(n), z);
    let mut sum1 = ZERO;
    if n > 0 {
        let mut fact = 1.0_f64;
        for k in 1..n {
            fact *= f64::from(k);
        }
        let mut term = Complex64::new(fact, 0.0);
        for k in 0..n {
            if k > 0 {
                term = -term * half_z2 / (f64::from(k) * f64::from(n - k));
            }
            sum1 += term;
        }
        sum1 *= 0.5 * half_z.powi(-(n as i32));
    }
    let mut h_k = -EULER_GAMMA;
    let mut h_nk = -EULER_GAMMA;
    let mut fact_nk = 1.0_f64;
    for j in 1..=n {
        h_nk += 1.0 / f64::from(j);
        fact_nk *= f64::from(j);
    }
    let mut term2 = ONE / fact_nk;
    let mut sum2 = term2 * (h_k + h_nk);
    for k in 1..80 {
        let kf = f64::from(k);
        h_k += 1.0 / kf;
        h_nk += 1.0 / (f64::from(n) + kf);
        term2 *= half_z2 / (kf * (f64::from(n) + kf));
        let add = term2 * (h_k + h_nk);
        sum2 += add;
        if add.norm() <= 1e-16 * (1.0 + sum2.norm()) {
            break;
        }
    }
    let sign_n = if n.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    };
    sum2 *= sign_n * 0.5 * half_z.powi(n as i32);
    let ln_term = -sign_n * half_z.ln() * i_n;
    sum1 + ln_term + sum2
}

/// Exponentially scaled `exp(z) * K_nu(z)` computed without `exp(z)` overflow for large `Re(z)`.
/// Returns NaN if the integer series would start with an overflowing factorial.
#[must_use]
pub fn bessel_ke_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO || !nu.is_finite() {
        return NAN;
    }
    let nu_abs = nu.abs();
    if z.norm() > 12.0 + 1.5 * nu_abs && z.arg().abs() < PI - 0.15 {
        let mu = 4.0 * nu_abs * nu_abs;
        let eight_z = 8.0 * z;
        let mut sum = ONE;
        let mut term = ONE;
        for k in 1..=25 {
            let kf = f64::from(k);
            let next_term = term * (mu - (2.0 * kf - 1.0) * (2.0 * kf - 1.0)) / (kf * eight_z);
            if next_term.norm() > term.norm() && k > 2 {
                break;
            }
            term = next_term;
            sum += term;
            if term.norm() < 1e-16 * sum.norm() {
                break;
            }
        }
        return (PI / (2.0 * z)).sqrt() * sum;
    }
    z.exp() * bessel_k_complex(nu, z)
}

/// Complex Hankel function of the first kind `H_nu^(1)(z) = J_nu(z) + i * Y_nu(z)`.
/// Uses `K_nu(-i*z)` connection in the upper half-plane (`Im(z) > 1.0`) to avoid cancellation.
#[must_use]
pub fn hankel_h1_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z.im > 1.0 && z.norm() > 2.0 {
        let phase = Complex64::from_polar(1.0, -0.5 * PI * nu);
        return (2.0 / (PI * I)) * phase * bessel_k_complex(nu, -I * z);
    }
    bessel_j_complex(nu, z) + I * bessel_y_complex(nu, z)
}

/// Complex Hankel function of the second kind `H_nu^(2)(z) = J_nu(z) - i * Y_nu(z)`.
/// Uses `K_nu(i*z)` connection in the lower half-plane (`Im(z) < -1.0`) to avoid cancellation.
#[must_use]
pub fn hankel_h2_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z.im < -1.0 && z.norm() > 2.0 {
        let phase = Complex64::from_polar(1.0, 0.5 * PI * nu);
        return -(2.0 / (PI * I)) * phase * bessel_k_complex(nu, I * z);
    }
    bessel_j_complex(nu, z) - I * bessel_y_complex(nu, z)
}

/// Spherical Bessel function `j_n(z) = sqrt(pi / (2z)) * J_{n + 1/2}(z)`.
#[must_use]
pub fn spherical_jn_complex(
    n: i32,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return if n == 0 { ONE } else { ZERO };
    }
    (0.5 * PI / z).sqrt() * bessel_j_complex(f64::from(n) + 0.5, z)
}

/// Spherical Bessel function `y_n(z) = sqrt(pi / (2z)) * Y_{n + 1/2}(z)`.
#[must_use]
pub fn spherical_yn_complex(
    n: i32,
    z: Complex64,
) -> Complex64 {
    (0.5 * PI / z).sqrt() * bessel_y_complex(f64::from(n) + 0.5, z)
}

/// Modified spherical Bessel function `i_n(z) = sqrt(pi / (2z)) * I_{n + 1/2}(z)`.
#[must_use]
pub fn spherical_in_complex(
    n: i32,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return if n == 0 { ONE } else { ZERO };
    }
    (0.5 * PI / z).sqrt() * bessel_i_complex(f64::from(n) + 0.5, z)
}

/// Modified spherical Bessel function `k_n(z) = sqrt(2 / (pi*z)) * K_{n + 1/2}(z)` (DLMF 10.47.9: `sqrt(pi/(2z)) * K_{n+1/2}(z)`).
#[must_use]
pub fn spherical_kn_complex(
    n: i32,
    z: Complex64,
) -> Complex64 {
    (0.5 * PI / z).sqrt() * bessel_k_complex(f64::from(n) + 0.5, z)
}

/// Complex Airy function `(Ai(z), Ai'(z))` via power series near the origin and
/// Bessel connection `(1/pi)*sqrt(z/3)*K_{1/3}(zeta)` / reflection for large `|z|`.
#[must_use]
pub fn airy_ai_complex(z: Complex64) -> (Complex64, Complex64) {
    // Power series converges rapidly for |z| <= 7.2 when Re(z) <= 2, or |z| <= 5.5 in RHP.
    if z.norm() <= 5.5 || (z.re <= 1.0 && z.norm() <= 7.2) {
        return airy_series(z, true);
    }
    // Connection with K_{1/3}, K_{2/3}:
    // For |arg(z)| <= 2*pi/3: Ai(z) = (1/pi) * sqrt(z/3) * K_{1/3}(zeta), zeta = (2/3) * z^(3/2)
    // Ai'(z) = -(z / (pi * sqrt(3))) * K_{2/3}(zeta)
    // For |arg(z)| > 2*pi/3 (near negative real axis): Ai(-w) = e^{i*pi/3} Ai(-w e^{i*pi/3}) + e^{-i*pi/3} Ai(-w e^{-i*pi/3})
    // Wait: DLMF 9.2.12: Ai(z) + e^{2pi i/3} Ai(z e^{2pi i/3}) + e^{-2pi i/3} Ai(z e^{-2pi i/3}) = 0!
    // So when |arg(z)| > 2*pi/3, rotate z by +- 2*pi/3 so both rotated args lie in |arg| <= 2*pi/3!
    let arg = z.arg();
    if arg.abs() > 2.0 * PI / 3.0 - 0.05 {
        let omega = Complex64::from_polar(1.0, 2.0 * PI / 3.0);
        let omega_conj = omega.conj();
        let (a1, ap1) = airy_ai_right_sector(z * omega);
        let (a2, ap2) = airy_ai_right_sector(z * omega_conj);
        let ai = -(omega * a1 + omega_conj * a2);
        let aip = -(omega * omega * ap1 + omega_conj * omega_conj * ap2);
        return (ai, aip);
    }
    airy_ai_right_sector(z)
}

fn airy_ai_right_sector(z: Complex64) -> (Complex64, Complex64) {
    let zeta = (2.0 / 3.0) * z * z.sqrt();
    let ai = (1.0 / PI) * (z / 3.0).sqrt() * bessel_k_complex(1.0 / 3.0, zeta);
    let aip = -(z / (PI * 3.0_f64.sqrt())) * bessel_k_complex(2.0 / 3.0, zeta);
    (ai, aip)
}

/// Complex Airy function `(Bi(z), Bi'(z))`.
#[must_use]
pub fn airy_bi_complex(z: Complex64) -> (Complex64, Complex64) {
    if z.norm() <= 6.5 {
        return airy_series(z, false);
    }
    // DLMF 9.2.11: Bi(z) = i * (e^{-2pi i/3} Ai(z e^{-2pi i/3}) - e^{2pi i/3} Ai(z e^{2pi i/3}))
    // Wait: let's check DLMF 9.2.11:
    // Bi(z) = e^{-pi i/6} Ai(z e^{-2pi i/3}) + e^{pi i/6} Ai(z e^{2pi i/3})
    // Wait: at z = 0: Ai(0) * (e^{-pi i/6} + e^{pi i/6}) = Ai(0) * sqrt(3) = Bi(0)! Exact!
    // And Bi'(z) = e^{-pi i/6 - 2pi i/3} Ai'(z e^{-2pi i/3}) + e^{pi i/6 + 2pi i/3} Ai'(z e^{2pi i/3})
    //            = -i * Ai'(z e^{-2pi i/3}) + i * Ai'(z e^{2pi i/3})!
    // Check at z = 0: Bi'(0) = 2 * Im(Ai'(0) * (-i)?) Wait: -i + i = 0?!
    // Let's check: e^{-i*pi/6} * e^{-2i*pi/3} = e^{-5i*pi/6}, e^{i*pi/6} * e^{2i*pi/3} = e^{5i*pi/6}.
    // At z = 0: (e^{-5i*pi/6} + e^{5i*pi/6}) * Ai'(0) = -sqrt(3) * Ai'(0) = +Bi'(0)! (since Ai'(0) < 0 and Bi'(0) > 0).
    let w_plus = Complex64::from_polar(1.0, 2.0 * PI / 3.0);
    let w_minus = w_plus.conj();
    let e_plus = Complex64::from_polar(1.0, PI / 6.0);
    let e_minus = e_plus.conj();
    let (a_m, ap_m) = airy_ai_complex(z * w_minus);
    let (a_p, ap_p) = airy_ai_complex(z * w_plus);
    let bi = e_minus * a_m + e_plus * a_p;
    let bip = (e_minus * w_minus) * ap_m + (e_plus * w_plus) * ap_p;
    (bi, bip)
}

fn airy_series(
    z: Complex64,
    is_ai: bool,
) -> (Complex64, Complex64) {
    // f(z) = sum_{k>=0} c_k z^{3k} / (3k)! , g(z) = sum_{k>=0} d_k z^{3k+1} / (3k+1)!
    // Ai(z) = c1 * f(z) - c2 * g(z), Bi(z) = sqrt(3) * (c1 * f(z) + c2 * g(z))
    const C1: f64 = 0.355_028_053_887_817_2; // 1 / (3^(2/3) * Gamma(2/3))
    const C2: f64 = 0.258_819_403_792_806_8; // 1 / (3^(1/3) * Gamma(1/3))
    let z3 = z * z * z;
    let mut f = ONE;
    let mut fp = ZERO;
    let mut g = z;
    let mut gp = ONE;
    let mut term_f = ONE;
    let mut term_g = z;
    for k in 1..40 {
        let k3 = f64::from(3 * k);
        let term_fp = term_f * z * z / (k3 - 1.0);
        fp += term_fp;
        term_f = term_fp * z / k3;
        f += term_f;

        let term_gp = term_g * z * z / k3;
        gp += term_gp;
        term_g = term_gp * z / (k3 + 1.0);
        g += term_g;

        if term_f.norm() + term_g.norm() < 1e-17 {
            break;
        }
    }
    let _ = z3;
    if is_ai {
        (C1 * f - C2 * g, C1 * fp - C2 * gp)
    } else {
        let s3 = 3.0_f64.sqrt();
        (s3 * (C1 * f + C2 * g), s3 * (C1 * fp + C2 * gp))
    }
}

/// Complex Kelvin functions `(ber_nu(z), bei_nu(z), ker_nu(z), kei_nu(z))`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "phases and continuations use IEEE complex arithmetic"
)]
pub fn kelvin_complex(
    nu: f64,
    z: Complex64,
) -> (Complex64, Complex64, Complex64, Complex64) {
    let e34 = Complex64::new(-FRAC_1_SQRT_2, FRAC_1_SQRT_2);
    let e14 = Complex64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2);
    let phase = Complex64::from_polar(1.0, -0.5 * nu * PI);
    let angle = z.arg();
    let j_rotated = |rotation: Complex64, turn: f64| {
        let point = z * rotation;
        let sheet = angle + turn - point.arg();
        let value = bessel_j_complex(nu, point);
        if sheet.abs() > PI {
            Complex64::from_polar(1.0, sheet.signum() * 2.0 * PI * nu) * value
        } else {
            value
        }
    };
    let k_rotated = |rotation: Complex64, turn: f64| {
        let point = z * rotation;
        let expected = angle + turn;
        if (expected - point.arg()).abs() > PI {
            // DLMF 10.34.2 at m=±1 continues continuously along the rotation.
            let side = expected.signum();
            let reflected = -point;
            Complex64::from_polar(1.0, -side * PI * nu) * bessel_k_complex(nu, reflected)
                - I * PI * side * bessel_i_complex(nu, reflected)
        } else {
            bessel_k_complex(nu, point)
        }
    };
    let j_plus = j_rotated(e34, 0.75 * PI);
    let k_plus = phase * k_rotated(e14, 0.25 * PI);
    if z.im == 0.0 && z.re >= 0.0 {
        // The negative rotations are complex conjugates on the positive real axis.
        return (
            j_plus.re.into(),
            j_plus.im.into(),
            k_plus.re.into(),
            k_plus.im.into(),
        );
    }
    let j_minus = j_rotated(e34.conj(), -0.75 * PI);
    let k_minus = phase.conj() * k_rotated(e14.conj(), -0.25 * PI);
    (
        0.5 * (j_plus + j_minus),
        -0.5 * I * (j_plus - j_minus),
        0.5 * (k_plus + k_minus),
        -0.5 * I * (k_plus - k_minus),
    )
}

/// Complex Struve function `H_nu(z)` via ascending power series (using `rgamma`
/// to handle negative half-integer `nu`) and large-`|z|` asymptotic expansion.
#[must_use]
pub fn struve_h_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    if z.norm() > 25.0 + 1.5 * nu.abs() {
        // DLMF 11.6.1: H_nu(z) ~ Y_nu(z) + (1/pi) * sum_{k>=0} Gamma(k + 1/2) / Gamma(nu + 1/2 - k) * (z/2)^(nu - 2k - 1)
        let half_z = 0.5 * z;
        let inv_hz2 = ONE / (half_z * half_z);
        let mut sum = ZERO;
        let mut term = crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(
            Complex64::new(nu + 0.5, 0.0),
        ) * PI.sqrt();
        for k in 0..16 {
            if k > 0 {
                let kf = f64::from(k);
                let next_term = term * (kf - 0.5) * (nu + 0.5 - kf) * inv_hz2;
                if next_term.norm() > term.norm() && k > 2 {
                    break;
                }
                term = next_term;
            }
            sum += term;
            if term.norm() <= 1e-16 * (1.0 + sum.norm()) {
                break;
            }
        }
        return bessel_y_complex(nu, z)
            + (1.0 / PI) * half_z.powc(Complex64::new(nu - 1.0, 0.0)) * sum;
    }
    let half_z = 0.5 * z;
    let half_z2 = half_z * half_z;
    let mut sum = ZERO;
    let mut g1 = PI.sqrt() * 0.5; // Gamma(3/2)
    let mut rg2 = crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(Complex64::new(
        nu + 1.5,
        0.0,
    ));
    let mut z_pow = ONE;
    for k in 0..150 {
        if k > 0 {
            let kf = f64::from(k);
            g1 *= kf + 0.5;
            rg2 /= nu + kf + 0.5;
            z_pow *= -half_z2;
        }
        let term = z_pow * rg2 / g1;
        sum += term;
        if k > 4 && term.norm() <= 1e-16 * (1.0 + sum.norm()) {
            break;
        }
    }
    half_z.powc(Complex64::new(nu + 1.0, 0.0)) * sum
}

/// Scorer function `Hi(z) = (1/pi) * int_0^oo exp(-t^3/3 + z*t) dt`.
#[must_use]
pub fn airy_hi_complex(z: Complex64) -> Complex64 {
    if z.norm() > 8.0 && z.arg().abs() > PI / 3.0 {
        return airy_bi_complex(z).0 - airy_gi_complex(z);
    }
    let mut sum = ZERO;
    let mut fact = 1.0_f64;
    let mut z_pow = ONE;
    let three_third = 3.0_f64.cbrt();
    for k in 0..80 {
        if k > 0 {
            fact *= f64::from(k);
            z_pow *= z;
        }
        let g = gamma_complex(Complex64::new((f64::from(k) + 1.0) / 3.0, 0.0));
        let scale = 3.0_f64.powf((f64::from(k) - 2.0) / 3.0);
        let term = g * scale * z_pow / fact;
        sum += term;
        if k > 6 && term.norm() <= 1e-16 * (1.0 + sum.norm()) {
            break;
        }
    }
    let _ = three_third;
    sum / PI
}

/// Scorer function `Gi(z) = (1/pi) * int_0^oo sin(t^3/3 + z*t) dt = Bi(z) - Hi(z)`.
#[must_use]
pub fn airy_gi_complex(z: Complex64) -> Complex64 {
    if z.norm() > 8.0 && z.arg().abs() > PI / 3.0 {
        // DLMF 9.12.26: Gi(z) ~ 1/(pi*z) * sum_{k>=0} (3k)! / (3^k * k! * z^(3k))
        let z3 = z * z * z;
        let mut sum = ONE;
        let mut term = ONE;
        for k in 1..15 {
            let kf = f64::from(k);
            let next_term = term * (3.0 * kf - 2.0) * (3.0 * kf - 1.0) / z3;
            if next_term.norm() > term.norm() {
                break;
            }
            term = next_term;
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        return sum / (PI * z);
    }
    airy_bi_complex(z).0 - airy_hi_complex(z)
}

/// Complex-order Bessel function of the first kind `J_nu(z)` on the principal branch.
#[must_use]
pub fn bessel_j_complex_order(
    nu: Complex64,
    z: Complex64,
) -> Complex64 {
    if nu.im == 0.0 {
        return bessel_j_complex(nu.re, z);
    }
    bessel_order(nu, z, false)
}

/// Complex-order modified Bessel function `I_nu(z)` on the principal branch.
#[must_use]
pub fn bessel_i_complex_order(
    nu: Complex64,
    z: Complex64,
) -> Complex64 {
    if nu.im == 0.0 {
        return bessel_i_complex(nu.re, z);
    }
    bessel_order(nu, z, true)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the floating-point series preserves error status and uses a finite term budget"
)]
fn bessel_order(
    nu: Complex64,
    z: Complex64,
    modified: bool,
) -> Complex64 {
    if !nu.is_finite() || !z.is_finite() || crate::cancellation::is_cancelled() {
        return NAN;
    }
    if z == ZERO {
        return if nu.re > 0.0 { ZERO } else { NAN };
    }
    let magnitude = z.norm();
    if magnitude >= 17.0 && magnitude > 0.5 * nu.norm_sqr() {
        if modified {
            // Rotate I into J in the right half-plane while preserving the
            // principal argument and selecting the corresponding cut lip.
            let side = z.im.signum();
            return bessel_j_asymptotic_complex_order(nu, (-side * I) * z, side);
        }
        return bessel_j_asymptotic_complex_order(nu, z, 0.0);
    }
    // DLMF 10.2.2 and 10.25.2 differ only in the sign of the series ratio.
    // Compute each seed with the power and reciprocal gamma in one exponential.
    let half = 0.5 * z;
    let logarithm = z.ln() - std::f64::consts::LN_2;
    let square = if modified {
        half * half
    } else {
        -half * half
    };
    let order_power = nu * logarithm;
    let mut term = (order_power - ln_gamma_complex(nu + ONE)).exp();
    let mut sum = term;
    let mut residual = ZERO;
    // A nearly zero denominator occurs close to a negative integer order. Restart
    // that term with its analytic coefficient, bypassing division by a tiny number.
    let near_integer = (-nu.re).round();
    let restart = if (1.0..=512.0).contains(&near_integer) && (nu + near_integer).norm_sqr() < 1e-8
    {
        near_integer as u32
    } else {
        0
    };
    for k in 1..=512 {
        let kf = f64::from(k);
        if k == restart {
            let mut log_power = (nu + 2.0 * kf) * logarithm;
            if !modified && k % 2 == 1 {
                log_power += I * PI;
            }
            term = (log_power
                - ln_gamma_complex(nu + (kf + 1.0))
                - ln_gamma_complex(Complex64::new(kf + 1.0, 0.0)).re)
                .exp();
        } else {
            term *= square / (kf * (nu + kf));
        }
        let corrected = term - residual;
        let next = sum + corrected;
        residual = (next - sum) - corrected;
        sum = next;
        let next_index = kf + 1.0;
        if k >= restart
            && term.norm() <= f64::EPSILON * sum.norm()
            && square.norm() <= 0.5 * next_index * (nu + next_index).norm()
            && next_index + nu.re > 0.0
        {
            return sum;
        }
        if !sum.is_finite() {
            return NAN;
        }
    }
    NAN
}

/// Complex modified Struve function `L_nu(z) = -i * exp(-i*nu*pi/2) * H_nu(i*z)`.
#[must_use]
pub fn struve_l_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    -I * Complex64::from_polar(1.0, -0.5 * nu * PI) * struve_h_complex(nu, I * z)
}

/// Complex Anger function `J_nu^A(z) = (1/pi) * int_0^pi cos(nu*t - z*sin(t)) dt`.
#[must_use]
pub fn anger_j_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    anger_weber(nu, z, f64::cos, Complex64::cos)
}

/// Complex Weber function `E_nu(z) = (1/pi) * int_0^pi sin(nu*t - z*sin(t)) dt`.
#[must_use]
pub fn weber_e_complex(
    nu: f64,
    z: Complex64,
) -> Complex64 {
    anger_weber(nu, z, f64::sin, Complex64::sin)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the integrand uses finite math and the shared quadrature bounds its work"
)]
fn anger_weber(
    nu: f64,
    z: Complex64,
    real_wave: fn(f64) -> f64,
    wave: fn(Complex64) -> Complex64,
) -> Complex64 {
    if !nu.is_finite() || !z.is_finite() || crate::cancellation::is_cancelled() {
        return NAN;
    }
    let phase = |theta: f64| wave(Complex64::new(nu * theta, 0.0) - z * theta.sin());
    let panels = ((nu.abs() + z.re.abs()).ceil() as usize).clamp(24, 4096);
    if z.im == 0.0 {
        // For real arguments the imaginary component vanishes identically.
        let result = gauss_kronrod_any(
            |t| real_wave(nu * t - z.re * t.sin()),
            0.0,
            PI,
            1e-12,
            panels,
        );
        return Complex64::new(result.value / PI, 0.0);
    }
    let tolerance = 1e-12 * z.im.abs().exp();
    let segment = ContourSegment::Line {
        start: ZERO,
        end: PI.into(),
    };
    integrate_segment_checked(|t| phase(t.re), segment, tolerance, panels).map_or(NAN, |result| {
        if result.converged {
            result.value / PI
        } else {
            NAN
        }
    })
}
