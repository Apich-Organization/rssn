//! Elliptic integrals of the third kind and Jacobi-function integrals.
//!
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
//! | [`EllipticPiPlan`] | reuses the complete integral for fixed parameters and varied amplitudes |
//! | [`jacobi_int_inv_sq`] | `int_0^u du / (1 - n sn^2(u, m))^2` in closed form |
//! | [`jacobi_int_sn2`] | `int_0^u sn^2` in terms of `E` and `F` |
//!
//! The parameter convention is `m = k^2`. Real paths for the third-kind integral
//! satisfy both `m sin^2 phi < 1` and `n sin^2 phi < 1` along the integration interval.

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::LN_2;
use std::f64::consts::PI;

use crate::kernels::numeric_common::scale_product;

use crate::kernels::special_ext::carlson_rf;
use crate::kernels::special_ext::elliptic_e_inc;
use crate::kernels::special_ext::elliptic_f;

/// Carlson's degenerate integral `R_C(x, y)` (`x >= 0`, `y != 0`; `y < 0` gives the Cauchy
/// principal value).
///
/// Finite inputs are handled by the real elementary formulas. Invalid
/// arguments return `NaN`.
#[must_use]
pub fn carlson_rc(
    x: f64,
    y: f64,
) -> f64 {
    if !x.is_finite() || !y.is_finite() || x < 0.0 || y == 0.0 {
        return f64::NAN;
    }
    if y < 0.0 {
        let root_x = x.sqrt();
        let root_y = (-y).sqrt();
        asinh_ratio(root_x, root_y) / root_x.hypot(root_y)
    } else if x > y {
        let d = (x - y).sqrt();
        asinh_ratio(d, y.sqrt()) / d
    } else if x < y {
        let d = (y - x).sqrt();
        d.atan2(x.sqrt()) / d
    } else {
        1.0 / y.sqrt()
    }
}

fn asinh_ratio(
    numerator: f64,
    denominator: f64,
) -> f64 {
    let quotient = numerator / denominator;
    if quotient.is_finite() {
        quotient.asinh()
    } else {
        numerator.ln() - denominator.ln() + LN_2
    }
}

const fn binary_scale(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

/// Carlson's `R_J(x, y, z, p)` for `x, y, z >= 0` (at most one zero) and `p > 0`.
///
/// Finite arguments are evaluated using homogeneity and duplication.
/// Invalid inputs return `NaN`.
#[must_use]
pub fn carlson_rj(
    x: f64,
    y: f64,
    z: f64,
    p: f64,
) -> f64 {
    // The direct duplication contains cubic products of its arguments.
    const DIRECT_MIN: f64 = binary_scale(-250);
    const DIRECT_MAX: f64 = binary_scale(250);

    let values = [x, y, z, p];
    if p <= 0.0
        || values
            .iter()
            .any(|&value| !value.is_finite() || value < 0.0)
        || values.iter().filter(|&&value| value == 0.0).count() > 1
    {
        return f64::NAN;
    }
    let maximum = values.into_iter().fold(0.0, f64::max);
    let minimum = values
        .into_iter()
        .filter(|&value| value > 0.0)
        .fold(f64::INFINITY, f64::min);
    if minimum >= DIRECT_MIN && maximum <= DIRECT_MAX {
        return rj_direct(values);
    }
    if minimum / maximum >= DIRECT_MIN {
        let unit = rj_direct(values.map(|value| value / maximum));
        return scale_product(unit, &[], &[maximum, maximum.sqrt()]).unwrap_or(f64::INFINITY);
    }

    let roots = values.map(f64::sqrt);
    // Accumulate with headroom before rounding at the subnormal boundary.
    let scaled = rj_from_roots(roots, 512);
    if scaled.is_finite() {
        scaled * binary_scale(-512)
    } else {
        rj_from_roots(roots, 0)
    }
}

fn rj_direct([mut xt, mut yt, mut zt, mut pt]: [f64; 4]) -> f64 {
    let (mut sum, mut fac) = (0.0, 1.0);
    let (mut delx, mut dely, mut delz, mut delp, mut ave);
    let mut iter = 0;
    loop {
        let (sx, sy, sz) = (xt.sqrt(), yt.sqrt(), zt.sqrt());
        let lam = sx * (sy + sz) + sy * sz;
        let alpha = (pt * (sx + sy + sz) + sx * sy * sz).powi(2);
        let beta = pt * (pt + lam).powi(2);
        sum += fac * rc_duplication(alpha, beta);
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
        if delx.abs().max(dely.abs()).max(delz.abs()).max(delp.abs()) < 5e-4 || iter >= 320 {
            break;
        }
    }
    3.0 * sum + fac * rj_correction([delx, dely, delz, delp]) / (ave * ave.sqrt())
}

fn rc_duplication(
    mut x: f64,
    mut y: f64,
) -> f64 {
    // The ordinary RJ iteration keeps both positive arguments in range.
    // Its late auxiliary RC calls are nearly equal before duplication.
    loop {
        let average = (x + y + y) / 3.0;
        let delta = (y - average) / average;
        if delta.abs() < 1.5e-3 {
            return (1.0
                + delta
                    * delta
                    * (0.3 + delta * (1.0 / 7.0 + delta * (0.375 + delta * 9.0 / 22.0))))
                / average.sqrt();
        }
        let lambda = 2.0 * x.sqrt() * y.sqrt() + y;
        x = 0.25 * (x + lambda);
        y = 0.25 * (y + lambda);
    }
}

fn rj_from_roots(
    mut roots: [f64; 4],
    output_shift: i32,
) -> f64 {
    let mut sum = 0.0;
    // Each duplication quarters differences across the binary64 exponent range.
    for step in 0..1100 {
        let [x, y, z, p] = roots;
        let divisors = [p + x, p + y, p + z];
        let [px, py, pz] = divisors.map(|d| p / d);
        let [qx, qy, qz] = [x / divisors[0], y / divisors[1], z / divisors[2]];
        // After one RC duplication its arguments are proportional to
        // 1 and 2b/(a+b), with a+b = (p+x)(p+y)(p+z) in root coordinates.
        let ratio = 2.0 * (px * py * pz + px * qy * qz + qx * py * qz + qx * qy * pz);
        let exponent = output_shift - 2 * step;
        sum += rj_weighted(6.0 * carlson_rc(1.0, ratio), exponent, &divisors);
        if !sum.is_finite() {
            return sum;
        }

        let [sx, sy, sz] = [x.sqrt(), y.sqrt(), z.sqrt()];
        let lambda = (sx * sy).hypot(sx * sz).hypot(sy * sz);
        roots = roots.map(|root| 0.5 * root.hypot(lambda));
        let scale = roots.into_iter().fold(0.0, f64::max);
        let [x, y, z, p] = roots.map(|root| (root / scale).powi(2));
        let average = 0.2 * (x + y + z + p + p);
        let deltas = [x, y, z, p].map(|value| (average - value) / average);
        if deltas.iter().all(|delta| delta.abs() < 5e-4) {
            let root_average = scale * average.sqrt();
            return sum + rj_weighted(rj_correction(deltas), exponent - 2, &[root_average; 3]);
        }
    }
    f64::NAN
}

fn rj_weighted(
    value: f64,
    exponent: i32,
    divisors: &[f64; 3],
) -> f64 {
    let first = exponent.clamp(-1022, 1023);
    let second = (exponent - first).clamp(-1022, 1023);
    scale_product(
        value,
        &[
            binary_scale(first),
            binary_scale(second),
            binary_scale(exponent - first - second),
        ],
        divisors,
    )
    .unwrap_or(f64::INFINITY)
}

fn rj_correction([delx, dely, delz, delp]: [f64; 4]) -> f64 {
    let (c1, c2, c3, c4) = (3.0 / 14.0, 1.0 / 3.0, 3.0 / 22.0, 3.0 / 26.0);
    let (c5, c6, c7, c8) = (0.75 * c3, 1.5 * c4, 0.5 * c2, c3 + c3);
    let ea = delx * (dely + delz) + dely * delz;
    let eb = delx * dely * delz;
    let ec = delp * delp;
    let ed = ea - 3.0 * ec;
    let ee = eb + 2.0 * delp * (ea - ec);
    1.0 + ed * (-c1 + c5 * ed - c6 * ee)
        + eb * (c7 + delp * (-c8 + delp * c4))
        + delp * ea * (c2 - delp * c3)
        - c2 * delp * ec
}

/// Complete elliptic integral of the third kind `Pi(n, m)` with `n < 1`, `m < 1`.
#[must_use]
pub fn elliptic_pi_complete(
    n: f64,
    m: f64,
) -> f64 {
    if m == 0.0 && n < 1.0 {
        return FRAC_PI_2 / (1.0 - n).sqrt();
    }
    if n == 0.0 {
        return carlson_rf(0.0, 1.0 - m, 1.0);
    }
    if let Some((complement, scale)) = pi_negative_reduction(n, m) {
        return FRAC_PI_2 / scale
            - complement / 3.0 * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - complement);
    }
    carlson_rf(0.0, 1.0 - m, 1.0) + n / 3.0 * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n)
}

/// Incomplete elliptic integral of the third kind.
///
/// `Pi(n; phi, m) = int_0^phi dt / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))` for any real `phi`
/// (`Pi(phi + j pi) = Pi(phi) + 2 j Pi(n, m)`). Continuation across whole
/// periods requires `n < 1` and `m < 1`. A principal-interval path is real
/// when both `n sin^2 t < 1` and `m sin^2 t < 1` throughout the interval.
#[must_use]
pub fn elliptic_pi(
    n: f64,
    phi: f64,
    m: f64,
) -> f64 {
    elliptic_pi_with_period(n, phi, m, || elliptic_pi_complete(n, m))
}

fn elliptic_pi_with_period(
    n: f64,
    phi: f64,
    m: f64,
    complete: impl FnOnce() -> f64,
) -> f64 {
    let j = (phi / PI).round();
    let r = phi - j * PI;
    if r.abs().to_bits() == FRAC_PI_2.to_bits() && n < 1.0 && m < 1.0 {
        return (2.0 * j + r.signum()) * complete();
    }
    let (s, c) = r.sin_cos();
    let delta2 = 1.0 - m * s * s;
    let base = if m == 0.0 && n < 1.0 {
        let k = (1.0 - n).sqrt();
        (k * s).atan2(c) / k
    } else if let Some((complement, scale)) = pi_negative_reduction(n, m) {
        let elementary = (scale * s).atan2(c * delta2.sqrt()) / scale;
        elementary
            - complement / 3.0
                * s
                * s
                * s
                * carlson_rj(c * c, delta2, 1.0, 1.0 - complement * s * s)
    } else {
        let first_kind = s * carlson_rf(c * c, delta2, 1.0);
        if n == 0.0 {
            first_kind
        } else {
            first_kind + n / 3.0 * s * s * s * carlson_rj(c * c, delta2, 1.0, 1.0 - n * s * s)
        }
    };
    if j == 0.0 {
        base
    } else {
        base + 2.0 * j * complete()
    }
}

fn pi_negative_reduction(
    n: f64,
    m: f64,
) -> Option<(f64, f64)> {
    if n >= -1.0 {
        return None;
    }
    let complement = m / n;
    // For the two characteristics n and m/n, their integrands minus the
    // first-kind integrand combine into the derivative of an arctangent.
    // Evaluate the small F - Pi(m/n) term directly through RJ.
    (complement.abs() <= 0.5).then(|| (complement, (1.0 - n).sqrt() * (1.0 - complement).sqrt()))
}

/// Prepared third-kind elliptic integral for varying amplitudes.
///
/// The characteristic and parameter are fixed. The complete integral is
/// evaluated once for the quasi-periodic continuation outside the principal
/// interval.
#[derive(Clone, Copy, Debug)]
pub struct EllipticPiPlan {
    n: f64,
    m: f64,
    complete: f64,
}

impl EllipticPiPlan {
    /// Prepares a complete integral and its incomplete continuation.
    ///
    /// # Errors
    ///
    /// The characteristic or parameter is non-finite or greater than or
    /// equal to one, or the computed complete integral is unrepresentable.
    pub fn new(
        n: f64,
        m: f64,
    ) -> Result<Self, &'static str> {
        if !n.is_finite() || !m.is_finite() || n >= 1.0 || m >= 1.0 {
            return Err("elliptic Pi requires finite n and m below one");
        }
        let complete = elliptic_pi_complete(n, m);
        if !complete.is_finite() || complete <= 0.0 {
            return Err("complete elliptic Pi is not representable");
        }
        Ok(Self { n, m, complete })
    }

    /// Returns the complete integral for the prepared parameters.
    #[must_use]
    pub const fn complete(&self) -> f64 {
        self.complete
    }

    /// Evaluates the incomplete integral at amplitude `phi` in radians.
    #[must_use]
    pub fn evaluate(
        &self,
        phi: f64,
    ) -> f64 {
        elliptic_pi_with_period(self.n, phi, self.m, || self.complete)
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
