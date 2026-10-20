//! # Elliptic integrals of the third kind and Jacobi-function helpers
//!
//! The closed-form solutions of the free asymmetric top ([`crate::sim::symplectic`]) and of
//! bound Kerr geodesics ([`crate::sim::relativity`]) need, besides the Jacobi functions and the
//! integrals of the first two kinds of [`crate::kernels::special_ext`], the elliptic integral of
//! the **third kind** and the integral of `1/(1 - n sn^2)^2`. They are provided here through
//! Carlson's symmetric forms.
//!
//! | Function | Meaning |
//! |---|---|
//! | [`carlson_rc`], [`carlson_rj`] | Carlson `R_C(x, y)` and `R_J(x, y, z, p)` (`p > 0`) |
//! | [`elliptic_pi_complete`] | `Pi(n, m) = int_0^{pi/2} dt / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))`, `n < 1` |
//! | [`elliptic_pi`] | incomplete `Pi(n; phi, m)` for any real `phi` (quasi-periodic continuation) |
//! | [`jacobi_int_inv_sq`] | `int_0^u du / (1 - n sn^2(u, m))^2` in closed form |
//! | [`jacobi_int_sn2`] | `int_0^u sn^2` in terms of `E` and `F` |
//!
//! Conventions: the *parameter* `m = k^2` is used everywhere. The characteristic must satisfy
//! `n sin^2 phi < 1` along the integration path (no singularity on the path); this holds for
//! all uses in this crate (`n < 1`, or `n < 0`).

use std::f64::consts::PI;

use crate::kernels::special_ext::{carlson_rf, elliptic_e_inc, elliptic_f, elliptic_k};

/// Carlson's degenerate integral `R_C(x, y)` (`x >= 0`, `y != 0`; `y < 0` gives the Cauchy
/// principal value).
#[must_use]
pub fn carlson_rc(
    x: f64,
    y: f64,
) -> f64 {
    let (mut xt, mut yt, w) = if y > 0.0 {
        (x, y, 1.0)
    } else {
        let xt = x - y;
        (xt, -y, x.sqrt() / xt.sqrt())
    };
    let mut ave;
    let mut s;
    let mut iter = 0;
    loop {
        let lam = 2.0 * xt.sqrt() * yt.sqrt() + yt;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        ave = (xt + yt + yt) / 3.0;
        s = (yt - ave) / ave;
        iter += 1;
        if s.abs() < 5e-4 || iter > 200 {
            break;
        }
    }
    w * (1.0 + s * s * (0.3 + s * (1.0 / 7.0 + s * (0.375 + s * 9.0 / 22.0)))) / ave.sqrt()
}

/// Carlson's `R_J(x, y, z, p)` for `x, y, z >= 0` (at most one zero) and `p > 0`.
#[must_use]
pub fn carlson_rj(
    x: f64,
    y: f64,
    z: f64,
    p: f64,
) -> f64 {
    let (c1, c2, c3, c4) = (3.0 / 14.0, 1.0 / 3.0, 3.0 / 22.0, 3.0 / 26.0);
    let (c5, c6, c7, c8) = (0.75 * c3, 1.5 * c4, 0.5 * c2, c3 + c3);
    let (mut xt, mut yt, mut zt, mut pt) = (x, y, z, p);
    let (mut sum, mut fac) = (0.0, 1.0);
    let (mut delx, mut dely, mut delz, mut delp, mut ave);
    let mut iter = 0;
    loop {
        let (sx, sy, sz) = (xt.sqrt(), yt.sqrt(), zt.sqrt());
        let lam = sx * (sy + sz) + sy * sz;
        let alpha = (pt * (sx + sy + sz) + sx * sy * sz).powi(2);
        let beta = pt * (pt + lam).powi(2);
        sum += fac * carlson_rc(alpha, beta);
        fac *= 0.25;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        zt = 0.25 * (zt + lam);
        pt = 0.25 * (pt + lam);
        ave = 0.2 * (xt + yt + zt + pt + pt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        delp = (ave - pt) / ave;
        iter += 1;
        if delx.abs().max(dely.abs()).max(delz.abs()).max(delp.abs()) < 5e-4 || iter > 200 {
            break;
        }
    }
    let ea = delx * (dely + delz) + dely * delz;
    let eb = delx * dely * delz;
    let ec = delp * delp;
    let ed = ea - 3.0 * ec;
    let ee = eb + 2.0 * delp * (ea - ec);
    3.0 * sum
        + fac
            * (1.0
                + ed * (-c1 + c5 * ed - c6 * ee)
                + eb * (c7 + delp * (-c8 + delp * c4))
                + delp * ea * (c2 - delp * c3)
                - c2 * delp * ec)
            / (ave * ave.sqrt())
}

/// Complete elliptic integral of the third kind `Pi(n, m)` with `n < 1`, `m < 1`.
#[must_use]
pub fn elliptic_pi_complete(
    n: f64,
    m: f64,
) -> f64 {
    carlson_rf(0.0, 1.0 - m, 1.0) + n / 3.0 * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n)
}

/// Incomplete elliptic integral of the third kind
/// `Pi(n; phi, m) = int_0^phi dt / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))` for any real `phi`
/// (`Pi(phi + j pi) = Pi(phi) + 2 j Pi(n, m)`); the integrand must not be singular, i.e.
/// `n < 1` (any `phi`) or `n sin^2 t < 1` on the path.
#[must_use]
pub fn elliptic_pi(
    n: f64,
    phi: f64,
    m: f64,
) -> f64 {
    let j = (phi / PI).round();
    let r = phi - j * PI;
    let (s, c) = r.sin_cos();
    let delta2 = 1.0 - m * s * s;
    let base = s * carlson_rf(c * c, delta2, 1.0)
        + n / 3.0 * s * s * s * carlson_rj(c * c, delta2, 1.0, 1.0 - n * s * s);
    if j == 0.0 {
        base
    } else {
        base + 2.0 * j * elliptic_pi_complete(n, m)
    }
}

/// `int_0^u sn^2(v, m) dv` for the amplitude `phi = am(u)` (any real `phi`):
/// `(u - E(phi, m)) / m` with `u = F(phi, m)`.
#[must_use]
pub fn jacobi_int_sn2(
    phi: f64,
    m: f64,
) -> f64 {
    if m < 1e-12 {
        // series: int sin^2 = phi/2 - sin(2 phi)/4 + O(m)
        return 0.5 * phi - 0.25 * (2.0 * phi).sin();
    }
    (elliptic_f(phi, m) - elliptic_e_inc(phi, m)) / m
}

/// `int_0^u dv / (1 - n sn^2(v, m))^2` with `u = F(phi, m)`, in closed form
/// (Byrd-Friedman 336.02) through `E`, `F` and `Pi`; requires `n != 1` and `n != m`.
#[must_use]
pub fn jacobi_int_inv_sq(
    n: f64,
    phi: f64,
    m: f64,
) -> f64 {
    let u = elliptic_f(phi, m);
    let e = elliptic_e_inc(phi, m);
    let p = elliptic_pi(n, phi, m);
    let (s, c) = phi.sin_cos();
    let dn = (1.0 - m * s * s).sqrt();
    // sn cn dn at u, with the sign of the branch of phi: cn = cos(phi)
    let sncndn = s * c * dn;
    let pref = 1.0 / (2.0 * (n - 1.0) * (m - n));
    pref * (n * e + (m - n) * u + (2.0 * n * m + 2.0 * n - n * n - 3.0 * m) * p
        - n * n * sncndn / (1.0 - n * s * s))
}
