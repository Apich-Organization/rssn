//! Complex orthogonal polynomials, spherical harmonics, Wright omega,
//! Whittaker M/W, parabolic cylinder D, and 0F1 confluent hypergeometric limit.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::hypergeometric::hyp_u_cplx;
use crate::kernels::hypergeometric::hyp1f1_cplx;

use crate::kernels::special_complex::scaled_reciprocal;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

/// Complex Legendre polynomial `P_n(z)`.
#[must_use]
pub fn legendre_p_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return z;
    }
    let mut p0 = ONE;
    let mut p1 = z;
    for k in 2..=n {
        let kf = f64::from(k);
        let next = ((2.0 * kf - 1.0) * z * p1 - (kf - 1.0) * p0) / kf;
        p0 = p1;
        p1 = next;
    }
    p1
}

/// Complex associated Legendre function `P_l^m(z)` for `0 <= m <= l`.
#[must_use]
pub fn assoc_legendre_p_complex(
    l: u32,
    m: u32,
    z: Complex64,
) -> Complex64 {
    if m > l {
        return ZERO;
    }
    let mut pmm = ONE;
    if m > 0 {
        let somx2 = (ONE - z * z).sqrt();
        let mut fact = 1.0_f64;
        for _ in 1..=m {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if l == m {
        return pmm;
    }
    let mut pmmp1 = z * (2.0 * f64::from(m) + 1.0) * pmm;
    if l == m + 1 {
        return pmmp1;
    }
    let mut pll = ZERO;
    for ll in (m + 2)..=l {
        let lf = f64::from(ll);
        let mf = f64::from(m);
        pll = ((2.0 * lf - 1.0) * z * pmmp1 - (lf + mf - 1.0) * pmm) / (lf - mf);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    pll
}

/// Complex Chebyshev polynomial of the first kind `T_n(z)`.
#[must_use]
pub fn chebyshev_t_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    chebyshev_complex(n, z, z)
}

/// Complex Chebyshev polynomial of the second kind `U_n(z)`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "initial value for the Chebyshev recurrence"
)]
pub fn chebyshev_u_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    chebyshev_complex(n, z, 2.0 * z)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "complex polynomial recurrence"
)]
fn chebyshev_complex(
    n: u32,
    z: Complex64,
    first: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    let double_z = 2.0 * z;
    let mut prev = ONE;
    let mut current = first;
    for _ in 2..=n {
        let next = double_z * current - prev;
        prev = current;
        current = next;
    }
    current
}

/// Complex physicists' Hermite polynomial `H_n(z)`.
#[must_use]
pub fn hermite_h_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return 2.0 * z;
    }
    let mut h0 = ONE;
    let mut h1 = 2.0 * z;
    for k in 2..=n {
        let next = 2.0 * z * h1 - 2.0 * f64::from(k - 1) * h0;
        h0 = h1;
        h1 = next;
    }
    h1
}

/// Complex generalized Laguerre polynomial `L_n^(alpha)(z)`.
#[must_use]
pub fn laguerre_l_complex(
    n: u32,
    alpha: Complex64,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return ONE + alpha - z;
    }
    let mut l0 = ONE;
    let mut l1 = ONE + alpha - z;
    for k in 2..=n {
        let kf = f64::from(k);
        let next = ((2.0 * kf - 1.0 + alpha - z) * l1 - (kf - 1.0 + alpha) * l0) / kf;
        l0 = l1;
        l1 = next;
    }
    l1
}

/// Complex Jacobi polynomial `P_n^(alpha, beta)(z)`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "complex polynomial recurrence"
)]
pub fn jacobi_p_complex(
    n: u32,
    alpha: Complex64,
    beta: Complex64,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    let mut p0 = ONE;
    let mut p1 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * z);
    if n == 1 {
        return p1;
    }
    let ab = alpha + beta;
    for k in 2..=n {
        let kf = f64::from(k);
        let two_k_ab = 2.0 * kf + ab;
        let degree_sum = kf + ab;
        let leading_factor = two_k_ab - 2.0;
        let tolerance = 16.0 * f64::EPSILON * (2.0 * kf + ab.norm());
        if degree_sum.norm() <= tolerance || leading_factor.norm() <= tolerance {
            return jacobi_binomial_sum(n, alpha, beta, z);
        }
        let c1 = 2.0 * kf * degree_sum * leading_factor;
        let c2 = (two_k_ab - 1.0) * (alpha * alpha - beta * beta + two_k_ab * (two_k_ab - 2.0) * z);
        let c3 = 2.0 * (kf - 1.0 + alpha) * (kf - 1.0 + beta) * two_k_ab;
        let next = (c2 * p1 - c3 * p0) / c1;
        p0 = p1;
        p1 = next;
    }
    p1
}

// DLMF 18.5.8 remains defined when the three-term recurrence has a
// vanishing leading coefficient. Use homogeneous Horner evaluation with
// the smaller of (z - 1)/2 and (z + 1)/2 as the explicit power.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "keeps the singular recurrence fallback separate"
)]
fn jacobi_binomial_sum(
    n: u32,
    mut alpha: Complex64,
    mut beta: Complex64,
    mut z: Complex64,
) -> Complex64 {
    let reflected = z.re < 0.0;
    if reflected {
        std::mem::swap(&mut alpha, &mut beta);
        z = -z;
    }
    let minus = 0.5 * (z - ONE);
    let plus = 0.5 * (z + ONE);
    let mut alpha_coefficients = vec![ONE];
    let mut alpha_coefficient = ONE;
    for k in 1..=n {
        alpha_coefficient *= (f64::from(n - k + 1) + alpha) / f64::from(k);
        alpha_coefficients.push(alpha_coefficient);
    }
    let mut beta_coefficient = ONE;
    let mut minus_power = ONE;
    let mut result = ZERO;
    for (k, alpha_coefficient) in (0..=n).zip(alpha_coefficients.into_iter().rev()) {
        result = result * plus + alpha_coefficient * beta_coefficient * minus_power;
        beta_coefficient *= (f64::from(n - k) + beta) / (f64::from(k) + 1.0);
        minus_power *= minus;
    }
    if reflected && n % 2 == 1 {
        -result
    } else {
        result
    }
}

/// Complex Gegenbauer (ultraspherical) polynomial `C_n^(lambda)(z)`.
#[must_use]
pub fn gegenbauer_c_complex(
    n: u32,
    lambda: Complex64,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return 2.0 * lambda * z;
    }
    let mut c0 = ONE;
    let mut c1 = 2.0 * lambda * z;
    for k in 2..=n {
        let kf = f64::from(k);
        let next = (2.0 * (kf + lambda - 1.0) * z * c1 - (kf + 2.0 * lambda - 2.0) * c0) / kf;
        c0 = c1;
        c1 = next;
    }
    c1
}

/// Complex spherical harmonic `Y_l^m(theta, phi)`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "normalized complex harmonic recurrence"
)]
pub fn spherical_harmonic_complex(
    l: u32,
    m: i32,
    theta: Complex64,
    phi: Complex64,
) -> Complex64 {
    let m_abs = m.unsigned_abs();
    if m_abs > l {
        return ZERO;
    }
    let z = theta.cos();
    let sin_theta = (ONE - z * z).sqrt();
    // Keep the associated polynomial normalized at every degree. Computing
    // P_l^m and the inverse factorial ratio separately can overflow and
    // underflow even when their product is near one.
    let mut diagonal = ONE / (4.0 * PI).sqrt();
    for order in 1..=m_abs {
        let twice_order = 2.0 * f64::from(order);
        diagonal *= -((twice_order + 1.0) / twice_order).sqrt() * sin_theta;
    }
    let polar = if l == m_abs {
        diagonal
    } else {
        let mut prev = diagonal;
        let mut current = (2.0 * f64::from(m_abs) + 3.0).sqrt() * z * diagonal;
        let mf = f64::from(m_abs);
        for degree in (m_abs + 1)..l {
            let k = f64::from(degree) + 1.0;
            let denominator = (k - mf) * (k + mf);
            let a = (((2.0 * k + 1.0) * (2.0 * k - 1.0)) / denominator).sqrt();
            let b = (((2.0 * k + 1.0) / (2.0 * k - 3.0)) * ((k - 1.0 - mf) * (k - 1.0 + mf))
                / denominator)
                .sqrt();
            let next = a * z * current - b * prev;
            prev = current;
            current = next;
        }
        current
    };
    let val = polar * (I * f64::from(m) * phi).exp();
    if m < 0 && m_abs % 2 == 1 {
        -val
    } else {
        val
    }
}

/// Complex Wright omega function `omega(z)` solving `omega + ln(omega) = z`
/// on the principal branch (Corless-Jeffrey Halley iteration).
#[must_use]
pub fn wright_omega_complex(z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if z.re < -100.0 {
        return z.exp();
    }
    if z.re > 100.0 {
        let l = z.ln();
        let mut w = z - l + l / z;
        for _ in 0..8 {
            let r = z - w - w.ln();
            let wp1 = ONE + w;
            let t = w * r;
            let denom = wp1 * wp1 + 0.5 * t;
            if denom.norm() < 1e-30 {
                break;
            }
            let delta = (t * wp1) / denom;
            w += delta;
            if delta.norm() <= 1e-15 * (1.0 + w.norm()) {
                break;
            }
        }
        return w;
    }
    let mut w = if z.re < -2.0 {
        z.exp()
    } else if z.norm() < 2.0 {
        Complex64::new(0.567_143_290_409_783_9, 0.0) + 0.36 * z
    } else {
        let l = z.ln();
        z - l
    };
    if w == ZERO {
        w = Complex64::new(1e-12, 0.0);
    }
    for _ in 0..32 {
        let r = z - w - w.ln();
        let wp1 = ONE + w;
        let t = w * r;
        let denom = wp1 * wp1 + 0.5 * t;
        let delta = (t * wp1) / denom;
        w += delta;
        if delta.norm() <= 1e-15 * (1.0 + w.norm()) {
            break;
        }
    }
    w
}

/// Confluent hypergeometric limit function `0F1(; b; z)`.
#[must_use]
pub fn hyp0f1_complex(
    b: Complex64,
    z: Complex64,
) -> Complex64 {
    if b.im == 0.0 && b.re <= 0.0 && (b.re - b.re.round()).abs() == 0.0 {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if z.norm() > 25.0 && b.im == 0.0 {
        // 0F1(; b; z) = Gamma(b) * (sqrt(z))^(1-b) * I_{b-1}(2*sqrt(z))
        let two_sq = 2.0 * z.sqrt();
        let nu = b.re - 1.0;
        let i_val = crate::kernels::special_complex::bessel_i_complex(nu, two_sq);
        let pref = crate::kernels::special_complex::gamma_complex(b) * z.sqrt().powc(ONE - b);
        return pref * i_val;
    }
    let mut sum = ONE;
    let mut term = ONE;
    for k in 1..300 {
        let kf = f64::from(k);
        term *= z / (kf * (b + (kf - 1.0)));
        sum += term;
        if term.norm() <= 1e-16 * (1.0 + sum.norm()) {
            break;
        }
    }
    sum
}

/// Whittaker function `M_{kappa, mu}(z) = exp(-z/2) * z^(mu + 1/2) * 1F1(mu - kappa + 1/2; 2*mu + 1; z)`.
#[must_use]
pub fn whittaker_m_complex(
    kappa: Complex64,
    mu: Complex64,
    z: Complex64,
) -> Complex64 {
    let half = Complex64::new(0.5, 0.0);
    (-0.5 * z).exp() * z.powc(mu + half) * hyp1f1_cplx(mu - kappa + half, 2.0 * mu + ONE, z)
}

/// Whittaker function `W_{kappa, mu}(z) = exp(-z/2) * z^(mu + 1/2) * U(mu - kappa + 1/2, 2*mu + 1, z)`.
#[must_use]
pub fn whittaker_w_complex(
    kappa: Complex64,
    mu: Complex64,
    z: Complex64,
) -> Complex64 {
    let half = Complex64::new(0.5, 0.0);
    (-0.5 * z).exp() * z.powc(mu + half) * hyp_u_cplx(mu - kappa + half, 2.0 * mu + ONE, z)
}

/// Parabolic cylinder function `D_nu(z)` via asymptotic expansion for large `|z|`
/// in `|arg(z)| < 3*pi/4` and `1F1` connection formula with `rgamma` elsewhere.
#[must_use]
pub fn parabolic_cylinder_d_complex(
    nu: Complex64,
    z: Complex64,
) -> Complex64 {
    if z.norm() > 8.0 && z.arg().abs() < 0.72 * PI {
        // DLMF 12.9.1: D_nu(z) ~ exp(-z^2/4) * z^nu * sum_{k>=0} (-1)^k (-nu)_{2k} / (k! (2 z^2)^k)
        let z2 = z * z;
        let inv_2z2 = 0.5 / z2;
        let mut sum = ONE;
        let mut term = ONE;
        for k in 1..=25 {
            let kf = f64::from(k);
            let a1 = 2.0 * kf - 2.0 - nu;
            let a2 = 2.0 * kf - 1.0 - nu;
            let next_term = -term * a1 * a2 * inv_2z2 / kf;
            if next_term.norm() > term.norm() && k > 2 {
                break;
            }
            term = next_term;
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        return (-0.25 * z2).exp() * z.powc(nu) * sum;
    }
    let z2 = z * z;
    let pref = Complex64::new(2.0, 0.0).powc(0.5 * nu) * (-0.25 * z2).exp();
    let sqrt_pi = PI.sqrt();
    let rg1 = crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(0.5 * (ONE - nu));
    let rg2 = crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(-0.5 * nu);
    let t1 = sqrt_pi * rg1 * hyp1f1_cplx(-0.5 * nu, Complex64::new(0.5, 0.0), 0.5 * z2);
    let t2 = 2.0_f64.sqrt()
        * sqrt_pi
        * rg2
        * z
        * hyp1f1_cplx(0.5 * (ONE - nu), Complex64::new(1.5, 0.0), 0.5 * z2);
    pref * (t1 - t2)
}

/// Complex Mittag-Leffler function `E_{alpha, beta}(z) = sum_{k>=0} z^k / Gamma(alpha * k + beta)`
/// for `Re(alpha) > 0`.
#[must_use]
pub fn mittag_leffler_complex(
    alpha: Complex64,
    beta: Complex64,
    z: Complex64,
) -> Complex64 {
    if z == ZERO {
        return crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(beta);
    }
    if (alpha - ONE).norm() < 1e-14 && (beta - ONE).norm() < 1e-14 {
        return z.exp();
    }
    let mut sum = ZERO;
    for k in 0..250 {
        let arg = alpha * f64::from(k) + beta;
        let rg = crate::kernels::special_complex_ext::gamma_beta::rgamma_complex(arg);
        let term = z.powi(k) * rg;
        sum += term;
        if k > 5 && term.norm() <= 1e-15 * (1.0 + sum.norm()) {
            break;
        }
    }
    sum
}

/// Principal complex Legendre function of the second kind `Q_n(z)` for integer degree.
/// The signed imaginary zero selects the upper or lower boundary value of the cut `[-1, 1]`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "the bounded recurrence handles floating-point branch selection and normalizes the minimal solution"
)]
pub fn legendre_q_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    if !z.is_finite() || (z.im == 0.0 && (z.re.abs() - 1.0).abs() == 0.0) {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if z.norm_sqr() >= 16.0 {
        return legendre_q_inverse_series(n, z);
    }
    let plus = Complex64::new(z.re + 1.0, z.im);
    let minus = Complex64::new(z.re - 1.0, z.im);
    let q0 = 0.5 * (plus.ln() - minus.ln());
    if n == 0 {
        return q0;
    }
    let root = plus.sqrt() * minus.sqrt();
    let first = z + root;
    let second = z - root;
    let dominant = if first.norm_sqr() > second.norm_sqr() {
        first
    } else {
        second
    };
    let decay = dominant.norm().ln();
    if f64::from(n) * decay > 4.0 {
        // Backward continued fractions select the minimal Q solution away
        // from the cut, where forward recurrence amplifies the dominant P.
        let padding = (24.0 / decay).ceil().clamp(8.0, 32_768.0) as u32;
        let end = n.saturating_add(padding);
        let mut ratio = ONE / dominant;
        let mut product = ONE;
        for k in (1..=end).rev() {
            let kf = f64::from(k);
            ratio = kf / ((2.0 * kf + 1.0) * z - (kf + 1.0) * ratio);
            if k <= n {
                product *= ratio;
            }
        }
        return q0 * product;
    }
    let mut previous = q0;
    let mut current = z * q0 - ONE;
    for k in 1..n {
        let kf = f64::from(k);
        let next = ((2.0 * kf + 1.0) * z * current - kf * previous) / (kf + 1.0);
        previous = current;
        current = next;
    }
    current
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "a separate inverse-power regime avoids overflow and cancellation in the recurrence"
)]
fn legendre_q_inverse_series(
    n: u32,
    z: Complex64,
) -> Complex64 {
    // DLMF 14.3.7 for integral degree and zero order, analytically
    // continued off the cut. The explicit power preserves the parity at -z.
    let inverse = scaled_reciprocal(z);
    let mut leading = inverse;
    for k in 1..=n {
        let kf = f64::from(k);
        leading *= (kf / (2.0 * kf + 1.0)) * inverse;
        if leading == ZERO {
            return ZERO;
        }
    }
    let square = inverse * inverse;
    let degree = f64::from(n);
    let mut term = leading;
    let mut sum = term;
    for k in 1..=2048 {
        let kf = f64::from(k);
        let ratio =
            ((0.5 * (degree - 1.0) + kf) * (0.5 * degree + kf)) / ((degree + 0.5 + kf) * kf);
        term *= square * ratio;
        sum += term;
        if ratio * square.norm() < 0.5 && term.norm() <= f64::EPSILON * sum.norm() {
            return sum;
        }
    }
    Complex64::new(f64::NAN, f64::NAN)
}
