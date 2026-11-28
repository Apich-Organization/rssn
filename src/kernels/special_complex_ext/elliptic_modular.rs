//! Complex elliptic integrals, functions, and modular forms.
//!
//! Includes Carlson symmetric integrals, complete and incomplete elliptic integrals, Jacobi
//! elliptic and theta functions, Weierstrass functions, the Dedekind eta function, and the Klein
//! j-invariant.

use std::f64::consts::PI;

use num_complex::Complex64;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

/// Complex Carlson symmetric elliptic integral `R_F(x, y, z)`.
#[must_use]
pub fn carlson_rf_complex(
    mut x: Complex64,
    mut y: Complex64,
    mut z: Complex64,
) -> Complex64 {
    for _ in 0..80 {
        let sx = x.sqrt();
        let sy = y.sqrt();
        let sz = z.sqrt();
        let lam = sx * sy + sy * sz + sz * sx;
        x = 0.25 * (x + lam);
        y = 0.25 * (y + lam);
        z = 0.25 * (z + lam);
        let a = (x + y + z) / 3.0;
        let dx = (a - x) / a;
        let dy = (a - y) / a;
        let dz = (a - z) / a;
        if dx.norm().max(dy.norm()).max(dz.norm()) < 1e-4 {
            let dx_dy = dx * dy;
            let e2 = dx_dy - dz * dz;
            let e3 = dx_dy * dz;
            let poly = ONE - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0;
            return poly / a.sqrt();
        }
    }
    let a = (x + y + z) / 3.0;
    ONE / a.sqrt()
}

/// Complex Carlson degenerate integral `R_C(x, y) = R_F(x, y, y)`.
#[must_use]
pub fn carlson_rc_complex(
    x: Complex64,
    y: Complex64,
) -> Complex64 {
    carlson_rf_complex(x, y, y)
}

/// Complex Carlson symmetric elliptic integral `R_D(x, y, z)`.
#[must_use]
pub fn carlson_rd_complex(
    mut x: Complex64,
    mut y: Complex64,
    mut z: Complex64,
) -> Complex64 {
    let mut sum = ZERO;
    let mut fac = ONE;
    for _ in 0..80 {
        let sx = x.sqrt();
        let sy = y.sqrt();
        let sz = z.sqrt();
        let lam = sx * sy + sy * sz + sz * sx;
        sum += fac / (sz * (z + lam));
        fac *= 0.25;
        x = 0.25 * (x + lam);
        y = 0.25 * (y + lam);
        z = 0.25 * (z + lam);
        let a = (x + y + 3.0 * z) / 5.0;
        let dx = (a - x) / a;
        let dy = (a - y) / a;
        let dz = (a - z) / a;
        if dx.norm().max(dy.norm()).max(dz.norm()) < 1e-4 {
            let ea = dx * dy;
            let eb = dz * dz;
            let ec = ea - eb;
            let ed = ea - 6.0 * eb;
            let ee = ed + ec + ec;
            let poly = ONE + (-3.0 / 14.0) * ed + (1.0 / 6.0) * ee * dz + (9.0 / 88.0) * ed * ed
                - (9.0 / 22.0) * dz * dz * ec;
            return 3.0 * sum + fac * poly / (a * a.sqrt());
        }
    }
    3.0 * sum
}

/// Complex Carlson symmetric elliptic integral `R_G(x, y, z)`.
#[must_use]
pub fn carlson_rg_complex(
    x: Complex64,
    y: Complex64,
    z: Complex64,
) -> Complex64 {
    let principal_sqrt = |value: Complex64| {
        let root = crate::kernels::complex_scalar::Complex::new(value.re, value.im).sqrt();
        Complex64::new(root.re, root.im)
    };
    let (x, mut y, mut z) = if y == ZERO {
        (y, x, z)
    } else if z == ZERO {
        (z, x, y)
    } else {
        (x, y, z)
    };
    if x == ZERO {
        // DLMF 19.20.4. Two zeros are allowed for R_G but not R_F or R_D.
        if y == ZERO {
            return 0.5 * principal_sqrt(z);
        }
        if z == ZERO {
            return 0.5 * principal_sqrt(y);
        }
        // The equal-argument limit also retains the selected signed side of
        // the negative real axis through the principal square root.
        if y == z && y.im.is_sign_negative() == z.im.is_sign_negative() {
            return (PI / 4.0) * principal_sqrt(z);
        }
        // For one zero, keep the larger nonzero magnitude in the R_D
        // denominator and preserve each argument's original branch side.
        if y.norm() > z.norm() {
            std::mem::swap(&mut y, &mut z);
        }
    }
    // DLMF 19.21.10 uses the product of the principal square roots.
    let roots = if x == ZERO {
        ZERO
    } else {
        principal_sqrt(x) * principal_sqrt(y) / principal_sqrt(z)
    };
    0.5 * (z * carlson_rf_complex(x, y, z) - (x - z) * (y - z) * carlson_rd_complex(x, y, z) / 3.0
        + roots)
}

/// Complex Carlson symmetric elliptic integral `R_J(x, y, z, p)`.
#[must_use]
pub fn carlson_rj_complex(
    mut x: Complex64,
    mut y: Complex64,
    mut z: Complex64,
    mut p: Complex64,
) -> Complex64 {
    let mut sum = ZERO;
    let mut fac = ONE;
    let delta = (p - x) * (p - y) * (p - z);
    for _ in 0..80 {
        let sx = x.sqrt();
        let sy = y.sqrt();
        let sz = z.sqrt();
        let sp = p.sqrt();
        let lam = sx * sy + sy * sz + sz * sx;
        let d = (sp + sx) * (sp + sy) * (sp + sz);
        let e = fac * fac * fac * delta / (d * d);
        sum += fac * carlson_rc_complex(ONE, ONE + e) / d;
        fac *= 0.25;
        x = 0.25 * (x + lam);
        y = 0.25 * (y + lam);
        z = 0.25 * (z + lam);
        p = 0.25 * (p + lam);
        let a = (x + y + z + 2.0 * p) / 5.0;
        let dx = (a - x) / a;
        let dy = (a - y) / a;
        let dz = (a - z) / a;
        let dp = (a - p) / a;
        if dx.norm().max(dy.norm()).max(dz.norm()).max(dp.norm()) < 1e-4 {
            let ea = dx * (dy + dz) + dy * dz;
            let eb = dx * dy * dz;
            let ec = dp * dp;
            let ed = ea - 3.0 * ec;
            let ee = eb + 2.0 * dp * (ea - ec);
            let poly =
                ONE - 3.0 * ed / 14.0 + ee / 6.0 + 9.0 * ed * ed / 88.0 - 3.0 * eb * dp / 11.0;
            return 3.0 * sum + fac * poly / (a * a.sqrt());
        }
    }
    3.0 * sum
}

/// Complete elliptic integral of the first kind `K(m)` for complex parameter `m`.
#[must_use]
pub fn elliptic_k_complex(m: Complex64) -> Complex64 {
    carlson_rf_complex(ZERO, ONE - m, ONE)
}

/// Complete elliptic integral of the second kind `E(m)` for complex parameter `m`.
#[must_use]
pub fn elliptic_e_complex(m: Complex64) -> Complex64 {
    2.0 * carlson_rg_complex(ZERO, ONE - m, ONE)
}

/// Incomplete elliptic integral of the first kind `F(phi | m)`.
#[must_use]
pub fn elliptic_f_complex(
    phi: Complex64,
    m: Complex64,
) -> Complex64 {
    let k_turn = (phi.re / PI).round();
    let phi0 = Complex64::new(phi.re - k_turn * PI, phi.im);
    let s = phi0.sin();
    let c = phi0.cos();
    let s2 = s * s;
    let base = s * carlson_rf_complex(c * c, ONE - m * s2, ONE);
    if k_turn == 0.0 {
        base
    } else {
        base + (2.0 * k_turn) * elliptic_k_complex(m)
    }
}

/// Incomplete elliptic integral of the second kind `E(phi | m)`.
#[must_use]
pub fn elliptic_e_inc_complex(
    phi: Complex64,
    m: Complex64,
) -> Complex64 {
    let k_turn = (phi.re / PI).round();
    let phi0 = Complex64::new(phi.re - k_turn * PI, phi.im);
    let s = phi0.sin();
    let c = phi0.cos();
    let s2 = s * s;
    let c2 = c * c;
    let one_m_ms2 = ONE - m * s2;
    let base = s
        * (carlson_rf_complex(c2, one_m_ms2, ONE)
            - (m * s2 / 3.0) * carlson_rd_complex(c2, one_m_ms2, ONE));
    if k_turn == 0.0 {
        base
    } else {
        base + (2.0 * k_turn) * elliptic_e_complex(m)
    }
}

/// Incomplete elliptic integral of the third kind `Pi(n; phi | m)`.
#[must_use]
pub fn elliptic_pi_complex(
    n: Complex64,
    phi: Complex64,
    m: Complex64,
) -> Complex64 {
    let k_turn = (phi.re / PI).round();
    let phi0 = Complex64::new(phi.re - k_turn * PI, phi.im);
    let s = phi0.sin();
    let c = phi0.cos();
    let s2 = s * s;
    let c2 = c * c;
    let one_m_ms2 = ONE - m * s2;
    let base = s
        * (carlson_rf_complex(c2, one_m_ms2, ONE)
            + (n * s2 / 3.0) * carlson_rj_complex(c2, one_m_ms2, ONE, ONE - n * s2));
    if k_turn == 0.0 {
        base
    } else {
        let complete = carlson_rf_complex(ZERO, ONE - m, ONE)
            + (n / 3.0) * carlson_rj_complex(ZERO, ONE - m, ONE, ONE - n);
        base + (2.0 * k_turn) * complete
    }
}

/// Complex Jacobi elliptic functions `(sn(u | m), cn(u | m), dn(u | m))`
/// via the arithmetic-geometric mean (Landen descending transform).
#[must_use]
pub fn jacobi_elliptic_complex(
    u: Complex64,
    m: Complex64,
) -> (Complex64, Complex64, Complex64) {
    const MAX_N: usize = 16;

    if m == ZERO {
        return (u.sin(), u.cos(), ONE);
    }
    if m == ONE {
        let th = u.tanh();
        let sech = ONE / u.cosh();
        return (th, sech, sech);
    }
    let mut a = [ZERO; MAX_N];
    let mut c = [ZERO; MAX_N];
    a[0] = ONE;
    let mut b = (ONE - m).sqrt();
    c[0] = m.sqrt();
    let mut n = 0_usize;
    while n + 1 < MAX_N && c.get(n).copied().unwrap_or(ZERO).norm() > 1e-15 {
        let an = a.get(n).copied().unwrap_or(ONE);
        let next_a = 0.5 * (an + b);
        let next_b = (an * b).sqrt();
        let next_c = 0.5 * (an - b);
        n += 1;
        if let Some(slot) = a.get_mut(n) {
            *slot = next_a;
        }
        b = next_b;
        if let Some(slot) = c.get_mut(n) {
            *slot = next_c;
        }
    }
    let scale = 2.0_f64.powi(i32::try_from(n).unwrap_or(0));
    let mut phi = scale * a.get(n).copied().unwrap_or(ONE) * u;
    for k in (1..=n).rev() {
        let ck = c.get(k).copied().unwrap_or(ZERO);
        let ak = a.get(k).copied().unwrap_or(ONE);
        phi = 0.5 * (phi + ((ck / ak) * phi.sin()).asin());
    }
    let sn = phi.sin();
    let cn = phi.cos();
    let dn = (ONE - m * sn * sn).sqrt();
    (sn, cn, dn)
}

/// Complex Jacobi theta functions `theta_k(z, q)` for `k in 1..=4` and `|q| < 1`.
#[must_use]
pub fn jacobi_theta_complex(
    k: u8,
    z: Complex64,
    q: Complex64,
) -> Complex64 {
    if k == 0 || k > 4 || q.norm() >= 1.0 {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    // Jacobi imaginary transformation when |q| > 0.7 (Im(tau) small):
    // tau = ln(q) / (i * pi), tau' = -1 / tau, q' = exp(i * pi * tau').
    if q.norm() > 0.7 && q != ZERO {
        let tau = q.ln() / (I * PI);
        let tau_p = -ONE / tau;
        let q_p = (I * PI * tau_p).exp();
        if q_p.norm() < q.norm() {
            let z_p = z * tau_p;
            let pref = (-I * tau_p).sqrt() * (I * tau_p * z * z / PI).exp();
            return match k {
                | 1 => -I * pref * jacobi_theta_complex(1, z_p, q_p),
                | 2 => pref * jacobi_theta_complex(4, z_p, q_p),
                | 3 => pref * jacobi_theta_complex(3, z_p, q_p),
                | 4 => pref * jacobi_theta_complex(2, z_p, q_p),
                | _ => Complex64::new(f64::NAN, f64::NAN),
            };
        }
    }
    // Reduce z.re modulo pi (or 2*pi for theta_1, theta_2) to keep sin/cos accurate.
    let z = if k <= 2 {
        let n_turn = (z.re / (2.0 * PI)).round();
        Complex64::new(z.re - n_turn * 2.0 * PI, z.im)
    } else {
        let n_turn = (z.re / PI).round();
        Complex64::new(z.re - n_turn * PI, z.im)
    };
    match k {
        | 1 => {
            // 2 q^(1/4) sum_{n>=0} (-1)^n q^{n(n+1)} sin((2n+1) z)
            let mut sum = ZERO;
            for n in 0..80 {
                let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                let q_pow = q.powi(n * (n + 1));
                let term = sign * q_pow * (f64::from(2 * n + 1) * z).sin();
                sum += term;
                if term.norm() < 1e-16 * (1.0 + sum.norm()) {
                    break;
                }
            }
            2.0 * q.powc(Complex64::new(0.25, 0.0)) * sum
        },
        | 2 => {
            // 2 q^(1/4) sum_{n>=0} q^{n(n+1)} cos((2n+1) z)
            let mut sum = ZERO;
            for n in 0..80 {
                let q_pow = q.powi(n * (n + 1));
                let term = q_pow * (f64::from(2 * n + 1) * z).cos();
                sum += term;
                if term.norm() < 1e-16 * (1.0 + sum.norm()) {
                    break;
                }
            }
            2.0 * q.powc(Complex64::new(0.25, 0.0)) * sum
        },
        | 3 => {
            // 1 + 2 sum_{n>=1} q^{n^2} cos(2 n z)
            let mut sum = ONE;
            for n in 1..80 {
                let q_pow = q.powi(n * n);
                let term = 2.0 * q_pow * (f64::from(2 * n) * z).cos();
                sum += term;
                if term.norm() < 1e-16 * (1.0 + sum.norm()) {
                    break;
                }
            }
            sum
        },
        | 4 => {
            // 1 + 2 sum_{n>=1} (-1)^n q^{n^2} cos(2 n z)
            let mut sum = ONE;
            for n in 1..80 {
                let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                let q_pow = q.powi(n * n);
                let term = 2.0 * sign * q_pow * (f64::from(2 * n) * z).cos();
                sum += term;
                if term.norm() < 1e-16 * (1.0 + sum.norm()) {
                    break;
                }
            }
            sum
        },
        | _ => Complex64::new(f64::NAN, f64::NAN),
    }
}

/// Derivative `d/dz theta_1(z, q)` at `z` (used for Weierstrass elliptic functions).
fn jacobi_theta1_prime(
    z: Complex64,
    q: Complex64,
) -> Complex64 {
    let mut sum = ZERO;
    for n in 0..80 {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        let m = f64::from(2 * n + 1);
        let q_pow = q.powi(n * (n + 1));
        let term = sign * m * q_pow * (m * z).cos();
        sum += term;
        if term.norm() < 1e-16 * (1.0 + sum.norm()) {
            break;
        }
    }
    2.0 * q.powc(Complex64::new(0.25, 0.0)) * sum
}

fn jacobi_theta1_triple_prime_zero(q: Complex64) -> Complex64 {
    let mut sum = ZERO;
    for n in 0..80 {
        let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
        let m = f64::from(2 * n + 1);
        let q_pow = q.powi(n * (n + 1));
        let term = sign * m * m * m * q_pow;
        sum += term;
        if term.norm() < 1e-16 * (1.0 + sum.norm()) {
            break;
        }
    }
    2.0 * q.powc(Complex64::new(0.25, 0.0)) * sum
}

/// Normalizes half-periods `(omega1, omega2)` so `Im(omega2 / omega1) > 0` and
/// returns `(tau, q, pi / (2 * omega1))`.
fn weierstrass_params(
    omega1: Complex64,
    omega2: Complex64,
) -> (Complex64, Complex64, Complex64) {
    let mut tau = omega2 / omega1;
    if tau.im < 0.0 {
        tau = -tau;
    }
    let q = (I * PI * tau).exp();
    let scale = Complex64::new(PI, 0.0) / (2.0 * omega1);
    (tau, q, scale)
}

/// Weierstrass elliptic function `wp(z; omega1, omega2)` via Jacobi theta functions.
#[must_use]
pub fn weierstrass_p_complex(
    z: Complex64,
    omega1: Complex64,
    omega2: Complex64,
) -> Complex64 {
    let (_, q, scale) = weierstrass_params(omega1, omega2);
    let w = scale * z;
    let t1 = jacobi_theta_complex(1, w, q);
    let t2 = jacobi_theta_complex(2, w, q);
    let _t2_0 = jacobi_theta_complex(2, ZERO, q);
    let t3_0 = jacobi_theta_complex(3, ZERO, q);
    let t4_0 = jacobi_theta_complex(4, ZERO, q);
    let ratio = t3_0 * t4_0 * t2 / t1;
    let e2 = (t3_0.powi(4) + t4_0.powi(4)) / 3.0;
    scale * scale * (ratio * ratio + e2)
}

/// Weierstrass elliptic derivative `wp'(z; omega1, omega2)` via central complex step or theta.
#[must_use]
pub fn weierstrass_p_prime_complex(
    z: Complex64,
    omega1: Complex64,
    omega2: Complex64,
) -> Complex64 {
    let (_, q, scale) = weierstrass_params(omega1, omega2);
    let w = scale * z;
    let t1 = jacobi_theta_complex(1, w, q);
    let t2 = jacobi_theta_complex(2, w, q);
    let t3 = jacobi_theta_complex(3, w, q);
    let t4 = jacobi_theta_complex(4, w, q);
    let t2_0 = jacobi_theta_complex(2, ZERO, q);
    let t3_0 = jacobi_theta_complex(3, ZERO, q);
    let t4_0 = jacobi_theta_complex(4, ZERO, q);
    let t1p_0 = t2_0 * t3_0 * t4_0; // theta1'(0) = theta2(0)*theta3(0)*theta4(0)
    let _ = jacobi_theta1_prime(ZERO, q);
    -2.0 * scale * scale * scale * (t1p_0 * t1p_0 * t1p_0) * (t2 * t3 * t4)
        / (t2_0 * t3_0 * t4_0 * t1 * t1 * t1)
}

/// Weierstrass invariants `(g2, g3)` for half-periods `(omega1, omega2)`.
#[must_use]
pub fn weierstrass_invariants_complex(
    omega1: Complex64,
    omega2: Complex64,
) -> (Complex64, Complex64) {
    let (_, q, scale) = weierstrass_params(omega1, omega2);
    let t2_0 = jacobi_theta_complex(2, ZERO, q);
    let t3_0 = jacobi_theta_complex(3, ZERO, q);
    let t4_0 = jacobi_theta_complex(4, ZERO, q);
    let a = t2_0.powi(4);
    let b = t3_0.powi(4);
    let c = t4_0.powi(4);
    // e1 = (b + c)/3 * scale^2, e2 = (a - c)/3 * scale^2, e3 = -(a + b)/3 * scale^2
    let s2 = scale * scale;
    let e1 = s2 * (b + c) / 3.0;
    let e2 = s2 * (a - c) / 3.0;
    let e3 = -s2 * (a + b) / 3.0;
    let g2 = -4.0 * (e1 * e2 + e2 * e3 + e3 * e1);
    let g3 = 4.0 * e1 * e2 * e3;
    (g2, g3)
}

/// Weierstrass zeta function `zeta_W(z; omega1, omega2)`.
#[must_use]
pub fn weierstrass_zeta_complex(
    z: Complex64,
    omega1: Complex64,
    omega2: Complex64,
) -> Complex64 {
    let (_, q, scale) = weierstrass_params(omega1, omega2);
    let w = scale * z;
    let t1 = jacobi_theta_complex(1, w, q);
    let t1p = jacobi_theta1_prime(w, q);
    let t1p0 = jacobi_theta1_prime(ZERO, q);
    let t1ppp0 = jacobi_theta1_triple_prime_zero(q);
    let eta1 = -(scale / 6.0) * (t1ppp0 / t1p0);
    (eta1 / scale) * w + scale * (t1p / t1)
}

/// Weierstrass sigma function `sigma(z; omega1, omega2)`.
#[must_use]
pub fn weierstrass_sigma_complex(
    z: Complex64,
    omega1: Complex64,
    omega2: Complex64,
) -> Complex64 {
    let (_, q, scale) = weierstrass_params(omega1, omega2);
    let w = scale * z;
    let t1 = jacobi_theta_complex(1, w, q);
    let t1p0 = jacobi_theta1_prime(ZERO, q);
    let t1ppp0 = jacobi_theta1_triple_prime_zero(q);
    let eta1 = -(scale / 6.0) * (t1ppp0 / t1p0);
    let quad = (0.5 * eta1 / (scale * scale)) * w * w;
    (quad.exp() / scale) * (t1 / t1p0)
}

/// Dedekind eta function `eta(tau)` on the upper half-plane `Im(tau) > 0`
/// using modular reduction to the fundamental domain and Euler's pentagonal theorem.
#[must_use]
pub fn dedekind_eta_complex(mut tau: Complex64) -> Complex64 {
    if tau.im <= 0.0 || !tau.re.is_finite() || !tau.im.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let mut factor = ONE;
    for _ in 0..64 {
        let shift = tau.re.round();
        if shift != 0.0 {
            tau -= shift;
            factor *= Complex64::from_polar(1.0, PI * shift / 12.0);
        }
        let n2 = tau.norm_sqr();
        if n2 >= 1.0 - 1e-12 {
            break;
        }
        // eta(-1/tau) = sqrt(-i*tau) * eta(tau) => eta(tau) = eta(-1/tau) / sqrt(-i*tau)
        factor /= (-I * tau).sqrt();
        tau = -ONE / tau;
    }
    let q = (2.0 * PI * I * tau).exp();
    let q_24 = (PI * I * tau / 12.0).exp();
    let mut sum = ONE;
    for n in 1..60 {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        let e1 = n * (3 * n - 1) / 2;
        let e2 = n * (3 * n + 1) / 2;
        let term = sign * (q.powi(e1) + q.powi(e2));
        sum += term;
        if term.norm() < 1e-16 * sum.norm() {
            break;
        }
    }
    factor * q_24 * sum
}

/// Klein modular j-invariant `j(tau)` on the upper half-plane `Im(tau) > 0`.
#[must_use]
pub fn klein_j_complex(mut tau: Complex64) -> Complex64 {
    if tau.im <= 0.0 || !tau.re.is_finite() || !tau.im.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    // Modularly reduce tau to fundamental domain |tau| >= 1, |Re(tau)| <= 0.5 since j is SL_2(Z)-invariant!
    for _ in 0..64 {
        tau.re -= tau.re.round();
        if tau.norm_sqr() >= 1.0 - 1e-12 {
            break;
        }
        tau = -ONE / tau;
    }
    let q = (I * PI * tau).exp();
    let t2 = jacobi_theta_complex(2, ZERO, q);
    let t3 = jacobi_theta_complex(3, ZERO, q);
    let lam = (t2 / t3).powi(4);
    let lam_squared = lam * lam;
    let complementary = ONE - lam;
    let num = (complementary + lam_squared).powi(3);
    let den = lam_squared * complementary * complementary;
    256.0 * num / den
}

/// Normalized Eisenstein series `E_4(tau) = 1 + 240 * sum_{n>=1} n^3 q^n / (1 - q^n)`
/// with `q = exp(2*pi*i*tau)` for `Im(tau) > 0`.
#[must_use]
pub fn eisenstein_e4_complex(mut tau: Complex64) -> Complex64 {
    if tau.im <= 0.0 || !tau.re.is_finite() || !tau.im.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let mut factor = ONE;
    for _ in 0..64 {
        tau.re -= tau.re.round();
        if tau.norm_sqr() >= 1.0 - 1e-12 {
            break;
        }
        factor *= tau.powi(4);
        tau = -ONE / tau;
    }
    let q = (2.0 * PI * I * tau).exp();
    let mut sum = ZERO;
    let mut qn = q;
    for n in 1..60 {
        let nf = f64::from(n);
        let term = (nf * nf * nf) * qn / (ONE - qn);
        sum += term;
        if term.norm() < 1e-16 * (1.0 + sum.norm()) {
            break;
        }
        qn *= q;
    }
    (ONE + 240.0 * sum) / factor
}

/// Normalized Eisenstein series `E_6(tau) = 1 - 504 * sum_{n>=1} n^5 q^n / (1 - q^n)`
/// with `q = exp(2*pi*i*tau)` for `Im(tau) > 0`.
#[must_use]
pub fn eisenstein_e6_complex(mut tau: Complex64) -> Complex64 {
    if tau.im <= 0.0 || !tau.re.is_finite() || !tau.im.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let mut factor = ONE;
    for _ in 0..64 {
        tau.re -= tau.re.round();
        if tau.norm_sqr() >= 1.0 - 1e-12 {
            break;
        }
        factor *= tau.powi(6);
        tau = -ONE / tau;
    }
    let q = (2.0 * PI * I * tau).exp();
    let mut sum = ZERO;
    let mut qn = q;
    for n in 1..60 {
        let nf = f64::from(n);
        let term = nf.powi(5) * qn / (ONE - qn);
        sum += term;
        if term.norm() < 1e-16 * (1.0 + sum.norm()) {
            break;
        }
        qn *= q;
    }
    (ONE - 504.0 * sum) / factor
}

/// Modular discriminant `Delta(tau) = (2*pi)^12 * eta(tau)^24` for `Im(tau) > 0`.
#[must_use]
pub fn modular_discriminant_complex(tau: Complex64) -> Complex64 {
    let eta = dedekind_eta_complex(tau);
    (2.0 * PI).powi(12) * eta.powi(24)
}

/// Modular lambda function `lambda(tau) = (theta_2(0, q) / theta_3(0, q))^4`, `q = exp(i*pi*tau)`.
#[must_use]
pub fn modular_lambda_complex(tau: Complex64) -> Complex64 {
    if tau.im <= 0.0 || !tau.re.is_finite() || !tau.im.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let q = (I * PI * tau).exp();
    let t2 = jacobi_theta_complex(2, ZERO, q);
    let t3 = jacobi_theta_complex(3, ZERO, q);
    (t2 / t3).powi(4)
}
