//! Complex orthogonal polynomials, spherical harmonics, Wright omega,
//! Whittaker M/W, parabolic cylinder D, and 0F1 confluent hypergeometric limit.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::hypergeometric::hyp_u_cplx;
use crate::kernels::hypergeometric::hyp1f1_cplx;
use crate::kernels::special_complex::gamma_complex;

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
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return z;
    }
    let mut t0 = ONE;
    let mut t1 = z;
    for _ in 2..=n {
        let next = 2.0 * z * t1 - t0;
        t0 = t1;
        t1 = next;
    }
    t1
}

/// Complex Chebyshev polynomial of the second kind `U_n(z)`.
#[must_use]
pub fn chebyshev_u_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return ONE;
    }
    if n == 1 {
        return 2.0 * z;
    }
    let mut u0 = ONE;
    let mut u1 = 2.0 * z;
    for _ in 2..=n {
        let next = 2.0 * z * u1 - u0;
        u0 = u1;
        u1 = next;
    }
    u1
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
    for k in 2..=n {
        let kf = f64::from(k);
        let ab = alpha + beta;
        let two_k_ab = 2.0 * kf + ab;
        let c1 = 2.0 * kf * (kf + ab) * (two_k_ab - 2.0);
        let c2 = (two_k_ab - 1.0) * (alpha * alpha - beta * beta + two_k_ab * (two_k_ab - 2.0) * z);
        let c3 = 2.0 * (kf - 1.0 + alpha) * (kf - 1.0 + beta) * two_k_ab;
        let next = (c2 * p1 - c3 * p0) / c1;
        p0 = p1;
        p1 = next;
    }
    p1
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
    let mut fact_ratio = 1.0_f64;
    for k in (l - m_abs + 1)..=(l + m_abs) {
        fact_ratio /= f64::from(k);
    }
    let norm = ((2.0 * f64::from(l) + 1.0) / (4.0 * PI) * fact_ratio).sqrt();
    let plm = assoc_legendre_p_complex(l, m_abs, theta.cos());
    let val = norm * plm * (I * f64::from(m) * phi).exp();
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
        return z - l + l / z;
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
    let mut sum = ONE;
    let mut term = ONE;
    for k in 1..300 {
        let kf = k as f64;
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

/// Parabolic cylinder function `D_nu(z)` via `1F1` connection formula (DLMF 12.7.12).
#[must_use]
pub fn parabolic_cylinder_d_complex(
    nu: Complex64,
    z: Complex64,
) -> Complex64 {
    let z2 = z * z;
    let pref = Complex64::new(2.0, 0.0).powc(0.5 * nu) * (-0.25 * z2).exp();
    let sqrt_pi = PI.sqrt();
    let t1 = (sqrt_pi / gamma_complex(0.5 * (ONE - nu)))
        * hyp1f1_cplx(-0.5 * nu, Complex64::new(0.5, 0.0), 0.5 * z2);
    let t2 = (2.0_f64.sqrt() * sqrt_pi / gamma_complex(-0.5 * nu))
        * z
        * hyp1f1_cplx(0.5 * (ONE - nu), Complex64::new(1.5, 0.0), 0.5 * z2);
    pref * (t1 - t2)
}
