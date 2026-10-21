//! # Exact bound Kerr geodesics in Mino time and black-hole shadows
//!
//! Timelike geodesics of the Kerr metric (Boyer-Lindquist coordinates, `G = c = 1`, mass `M`,
//! spin `a`) with the constants of motion energy `E`, axial angular momentum `L` and Carter
//! constant `Q`. In **Mino time** `d lambda = d tau / Sigma` the radial and polar motions
//! decouple,
//!
//! ```text
//! (dr/dlambda)^2 = R(r) = [E (r^2 + a^2) - a L]^2 - Delta [r^2 + (L - a E)^2 + Q]
//! (dz/dlambda)^2 = Q - (Q + L^2 + a^2 (1 - E^2)) z^2 + a^2 (1 - E^2) z^4,   z = cos(theta)
//! ```
//!
//! and for bound orbits they are solved in closed form (Fujita and Hikida 2009):
//! `r(lambda)` is a Moebius function of `sn^2` with the parameter
//! `k_r^2 = (r1 - r2)(r3 - r4)/((r1 - r3)(r2 - r4))`, `z(lambda) = z_- sn(...)`. The coordinate
//! time and azimuth `t(lambda)`, `phi(lambda)` split into radial and polar parts; every part is a
//! combination of `lambda`, elliptic integrals of the first, second and third kind
//! ([`crate::sim::elliptic`]), evaluated for any `lambda` (no accumulation of integration error).
//!
//! | Capability | Notes |
//! |---|---|
//! | [`KerrOrbit::new`] | bound (`E < 1`) orbits with four real radial roots `r1 > r2 > r3 >= r4`, `r2` above the horizon; Schwarzschild (`a = 0`), equatorial (`Q = 0`) and inclined orbits |
//! | [`KerrOrbit::from_apsides`] | constants of motion from periapsis, apoapsis and the turning latitude |
//! | [`KerrOrbit::state`] | `(t, r, theta, phi)` at any Mino time |
//! | [`KerrOrbit::frequencies`] | Mino frequencies `Upsilon_r, Upsilon_theta, Upsilon_phi`, `Gamma` and the coordinate-time frequencies |
//! | [`shadow_boundary_analytic`] | Bardeen's critical curve of the shadow seen by a distant observer |
//! | [`trace_shadow`] | backward ray tracing of null geodesics (Mino-time integration) to locate the shadow boundary numerically |
//!
//! Limitations: the closed form covers bound timelike orbits (`mu = 1`, `E < 1`) without
//! radiation reaction and excludes the separatrix and circular orbits (`r1 = r2`, where the
//! radial motion degenerates to a point) and extremal spin; ray tracing is numerical.

use std::f64::consts::PI;

use crate::kernels::special_ext::{elliptic_f, elliptic_k, jacobi_elliptic};
use crate::sim::elliptic::{elliptic_pi, jacobi_int_inv_sq, jacobi_int_sn2};

/// Bound Kerr geodesic in closed form.
#[derive(Clone, Debug)]
pub struct KerrOrbit {
    /// Black-hole mass.
    pub m: f64,
    /// Spin parameter `a = J/M`.
    pub a: f64,
    /// Energy per unit rest mass.
    pub energy: f64,
    /// Axial angular momentum.
    pub lz: f64,
    /// Carter constant.
    pub carter: f64,
    /// Radial roots `r1 > r2 > r3 >= r4` of `R(r)`.
    pub roots: [f64; 4],
    /// Parameter `k_r^2` of the radial Jacobi functions.
    pub m_r: f64,
    /// Parameter `k_theta^2` of the polar Jacobi functions.
    pub m_th: f64,
    /// Turning value `z_- = cos(theta_min)` of the polar motion.
    pub z_minus: f64,
    g_r: f64,
    rate_r: f64,
    rate_th: f64,
    u0: f64,
    w0: f64,
    phi0: f64,
}

/// Coefficients `[c4, c3, c2, c1, c0]` of `R(r)` for a unit-mass particle.
fn radial_poly(
    m: f64,
    a: f64,
    e: f64,
    l: f64,
    q: f64,
) -> [f64; 5] {
    let k = (l - a * e).powi(2) + q;
    [
        e * e - 1.0,
        2.0 * m,
        a * a * (e * e - 1.0) - l * l - q,
        2.0 * m * k,
        (e * a * a - a * l).powi(2) - a * a * k,
    ]
}

fn poly_eval(
    c: &[f64; 5],
    x: f64,
) -> f64 {
    c.iter().fold(0.0, |acc, &ck| acc * x + ck)
}

/// Real roots of the quartic by Durand-Kerner iteration; `None` unless all four are real.
fn quartic_real_roots(c: &[f64; 5]) -> Option<[f64; 4]> {
    if c[0] == 0.0 {
        return None;
    }
    let mon: Vec<f64> = c.iter().map(|x| x / c[0]).collect();
    let radius = 1.0 + mon[1..].iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    let mut z: Vec<(f64, f64)> = (0..4)
        .map(|k| {
            let ang = 0.4 + 2.0 * PI * f64::from(k) / 4.0;
            (0.6 * radius * ang.cos(), 0.6 * radius * ang.sin())
        })
        .collect();
    let cmul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let cdiv = |a: (f64, f64), b: (f64, f64)| {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    };
    for _ in 0..500 {
        let mut delta = 0.0_f64;
        for i in 0..4 {
            let mut p = (0.0, 0.0);
            for &ck in &mon {
                p = cmul(p, z[i]);
                p.0 += ck;
            }
            let mut den = (1.0, 0.0);
            for j in 0..4 {
                if j != i {
                    den = cmul(den, (z[i].0 - z[j].0, z[i].1 - z[j].1));
                }
            }
            let step = cdiv(p, den);
            z[i].0 -= step.0;
            z[i].1 -= step.1;
            delta = delta.max(step.0.abs() + step.1.abs());
        }
        if delta < 1e-15 * radius {
            break;
        }
    }
    let mut out = [0.0; 4];
    for (i, zi) in z.iter().enumerate() {
        if zi.1.abs() > 1e-7 * radius.max(1.0) {
            return None;
        }
        out[i] = zi.0;
    }
    // polish with Newton on the real polynomial
    for r in &mut out {
        for _ in 0..4 {
            let f = poly_eval(c, *r);
            let df = ((4.0 * c[0] * *r + 3.0 * c[1]) * *r + 2.0 * c[2]) * *r + c[3];
            if df != 0.0 {
                *r -= f / df;
            }
        }
    }
    out.sort_by(|a, b| b.total_cmp(a));
    Some(out)
}

/// Partial-fraction data of `(c1 r + c0)/Delta`: `[(gamma_+, r_+), (gamma_-, r_-)]`.
fn lin_over_delta(
    c1: f64,
    c0: f64,
    rp: f64,
    rm: f64,
) -> [(f64, f64); 2] {
    [((c1 * rp + c0) / (rp - rm), rp), ((c1 * rm + c0) / (rm - rp), rm)]
}

impl KerrOrbit {
    /// Closed-form bound geodesic with constants `(E, L, Q)` starting at `(r0, theta0, phi0)`;
    /// `radial_out` / `polar_down` give the signs of `dr/dlambda` and of `dtheta/dlambda`
    /// (`polar_down`: theta increasing, i.e. moving toward the south pole).
    ///
    /// # Errors
    /// Fails unless `E < 1`, the radial polynomial has four real roots with `r1 > r2 > r3`,
    /// `r2` is outside the horizon, `r0` lies in `[r2, r1]` and the polar motion is allowed.
    pub fn new(
        m: f64,
        a: f64,
        (energy, lz, carter): (f64, f64, f64),
        (r0, theta0, phi0): (f64, f64, f64),
        (radial_out, polar_down): (bool, bool),
    ) -> Result<Self, String> {
        if m <= 0.0 || a.abs() >= m || energy >= 1.0 || energy <= 0.0 || carter < 0.0 {
            return Err("need 0 < E < 1, Q >= 0, M > 0 and |a| < M".into());
        }
        let c = radial_poly(m, a, energy, lz, carter);
        let roots = quartic_real_roots(&c).ok_or("radial polynomial has complex roots: not a bound orbit")?;
        let [r1, r2, r3, r4] = roots;
        let rh = m + (m * m - a * a).sqrt();
        if !(r1 > r2 && r2 > r3 && r3 >= r4 - 1e-12) || r2 <= rh {
            return Err("bound orbit needs r1 > r2 > r3 and r2 above the horizon (plunging or circular)".into());
        }
        if r0 < r2 * (1.0 - 1e-12) || r0 > r1 * (1.0 + 1e-12) {
            return Err("initial radius outside the radial turning points".into());
        }
        let m_r = (r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4));
        let g_r = 2.0 / ((r1 - r3) * (r2 - r4)).sqrt();
        let rate_r = (1.0 - energy * energy).sqrt() / g_r;
        // radial phase
        let s0 = ((r1 - r3) * (r0 - r2) / ((r1 - r2) * (r0 - r3))).clamp(0.0, 1.0);
        let psi0 = s0.sqrt().asin();
        let amp = if radial_out { psi0 } else { PI - psi0 };
        let u0 = elliptic_f(amp, m_r);
        // polar motion
        let beta = a * a * (1.0 - energy * energy);
        let b = carter + lz * lz + beta;
        let disc = (b * b - 4.0 * beta * carter).max(0.0).sqrt();
        let (z_minus, m_th, rate_th) = if carter == 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let zm2 = 2.0 * carter / (b + disc);
            (zm2.sqrt(), beta * zm2 * zm2 / carter, (carter / zm2).sqrt())
        };
        let z0 = theta0.cos();
        if z0.abs() > z_minus * (1.0 + 1e-12) + 1e-14 {
            return Err("initial latitude outside the polar turning points".into());
        }
        let w0 = if z_minus == 0.0 {
            0.0
        } else {
            let sn0 = (z0 / z_minus).clamp(-1.0, 1.0);
            let cn0 = (1.0 - sn0 * sn0).max(0.0).sqrt();
            // dz/dlambda has the sign of cn; theta increasing means z decreasing
            let cn0 = if polar_down { -cn0 } else { cn0 };
            elliptic_f(sn0.atan2(cn0), m_th)
        };
        Ok(Self { m, a, energy, lz, carter, roots, m_r, m_th, z_minus, g_r, rate_r, rate_th, u0, w0, phi0 })
    }

    /// Orbit specified by periapsis `r_p`, apoapsis `r_a` and the turning latitude `z_minus =
    /// cos(theta_min)` (`0` = equatorial); `prograde` selects the sign of `L`. The constants
    /// `(E, L, Q)` are found by Newton iteration. The orbit starts at periapsis on the
    /// equator-crossing side of the polar motion.
    ///
    /// # Errors
    /// Fails if the iteration does not converge or the resulting orbit is not bound.
    pub fn from_apsides(
        m: f64,
        a: f64,
        (r_p, r_a): (f64, f64),
        z_minus: f64,
        prograde: bool,
    ) -> Result<Self, String> {
        if r_p <= 0.0 || r_a <= r_p || !(0.0..1.0).contains(&z_minus) {
            return Err("need 0 < r_p < r_a and 0 <= z_minus < 1".into());
        }
        let p = 2.0 * r_p * r_a / (r_p + r_a);
        let e_ecc = (r_a - r_p) / (r_a + r_p);
        let sgn = if prograde { 1.0 } else { -1.0 };
        // Newton on (E, L, Q): R(r_p) = R(r_a) = 0 and Theta(z_minus) = 0
        let resid = |x: &[f64; 3]| -> [f64; 3] {
            let c = radial_poly(m, a, x[0], x[1], x[2]);
            let beta = a * a * (1.0 - x[0] * x[0]);
            let z2 = z_minus * z_minus;
            [
                poly_eval(&c, r_p),
                poly_eval(&c, r_a),
                x[2] - (x[2] + x[1] * x[1] + beta) * z2 + beta * z2 * z2,
            ]
        };
        let l_guess = sgn * (m * p).sqrt() * (1.0 + 3.0 * m / p).sqrt();
        let q_guess = l_guess * l_guess * z_minus * z_minus / (1.0 - z_minus * z_minus);
        let mut x = [(1.0 - (1.0 - e_ecc * e_ecc) * m / p).sqrt(), l_guess / (1.0 + q_guess / (l_guess * l_guess + 1e-300)).sqrt(), q_guess];
        let scale = [1.0, 1.0 + x[1].abs(), 1.0 + q_guess];
        for _ in 0..80 {
            let f = resid(&x);
            let fn_ = f.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if fn_ < 1e-13 * (1.0 + r_a.powi(4)) {
                break;
            }
            let mut jac = [[0.0; 3]; 3];
            for k in 0..3 {
                let eps = 1e-7 * scale[k];
                let mut xp = x;
                xp[k] += eps;
                let mut xm = x;
                xm[k] -= eps;
                let (fp, fm) = (resid(&xp), resid(&xm));
                for i in 0..3 {
                    jac[i][k] = (fp[i] - fm[i]) / (2.0 * eps);
                }
            }
            let dx = solve3_sys(&jac, &f).ok_or("singular Jacobian while solving for the constants of motion")?;
            // damped update
            let mut lam = 1.0;
            for _ in 0..20 {
                let trial = [x[0] - lam * dx[0], x[1] - lam * dx[1], x[2] - lam * dx[2]];
                let ft = resid(&trial);
                if trial[0] > 0.0 && trial[0] < 1.0 && trial[2] >= 0.0 && ft.iter().map(|v| v.abs()).fold(0.0, f64::max) < fn_ {
                    x = trial;
                    break;
                }
                lam *= 0.5;
            }
        }
        if resid(&x).iter().map(|v| v.abs()).fold(0.0, f64::max) > 1e-9 * (1.0 + r_a.powi(4)) {
            return Err("constants of motion did not converge".into());
        }
        let roots = quartic_real_roots(&radial_poly(m, a, x[0], x[1], x[2])).ok_or("complex radial roots: not a bound orbit")?;
        Self::new(m, a, (x[0], x[1], x[2]), (roots[1], std::f64::consts::FRAC_PI_2, 0.0), (true, true))
    }

    /// Mino period of the radial motion.
    #[must_use]
    pub fn radial_period(&self) -> f64 {
        2.0 * elliptic_k(self.m_r) / self.rate_r
    }

    /// Mino period of the polar motion (`None` for equatorial orbits).
    #[must_use]
    pub fn polar_period(&self) -> Option<f64> {
        (self.z_minus > 0.0).then(|| 4.0 * elliptic_k(self.m_th) / self.rate_th)
    }

    /// Continuous Jacobi amplitude at argument `u` (parameter `m`).
    fn amplitude(
        u: f64,
        m: f64,
    ) -> f64 {
        let (sn, cn, _) = jacobi_elliptic(u, m);
        let phi_mod = sn.atan2(cn);
        let kk = elliptic_k(m);
        let wrap = ((u - elliptic_f(phi_mod, m)) / (4.0 * kk)).round();
        phi_mod + 2.0 * PI * wrap
    }

    /// Radius at Mino time `lambda` (measured from the start of the orbit).
    #[must_use]
    pub fn radius(
        &self,
        lambda: f64,
    ) -> f64 {
        let [r1, r2, r3, _] = self.roots;
        let u = self.u0 + self.rate_r * lambda;
        let (sn, _, _) = jacobi_elliptic(u, self.m_r);
        let s = sn * sn;
        ((r1 - r3) * r2 - s * (r1 - r2) * r3) / ((r1 - r3) - s * (r1 - r2))
    }

    /// `cos(theta)` at Mino time `lambda`.
    #[must_use]
    pub fn cos_theta(
        &self,
        lambda: f64,
    ) -> f64 {
        if self.z_minus == 0.0 {
            return 0.0;
        }
        let w = self.w0 + self.rate_th * lambda;
        let (sn, _, _) = jacobi_elliptic(w, self.m_th);
        self.z_minus * sn
    }

    /// Radial parts `(T_r, Phi_r)` of the antiderivatives of `dt/dlambda` and `dphi/dlambda`
    /// from the orbit start to Mino time `lambda`.
    fn radial_parts(
        &self,
        lambda: f64,
    ) -> (f64, f64) {
        let [r1, r2, r3, _] = self.roots;
        let (m, a, e, l) = (self.m, self.a, self.energy, self.lz);
        let disc = (m * m - a * a).sqrt();
        let (rp, rm) = (m + disc, m - disc);
        let (aa, bb, cc, dd) = ((r1 - r3) * r2, -(r1 - r2) * r3, r1 - r3, -(r1 - r2));
        let u1 = self.u0 + self.rate_r * lambda;
        let psi1 = Self::amplitude(u1, self.m_r);
        let psi0 = Self::amplitude(self.u0, self.m_r);
        let mm = self.m_r;
        let du = u1 - self.u0;
        // basic integrals over u from u0 to u1
        let pi_n = |n: f64| elliptic_pi(n, psi1, mm) - elliptic_pi(n, psi0, mm);
        let pi2_n = |n: f64| jacobi_int_inv_sq(n, psi1, mm) - jacobi_int_inv_sq(n, psi0, mm);
        let i1 = du;
        let n_r = -dd / cc;
        let kap = aa - bb * cc / dd;
        let i_r = (bb / dd) * du + (kap / cc) * pi_n(n_r);
        let i_r2 = (bb / dd).powi(2) * du + 2.0 * (bb / dd) * kap / cc * pi_n(n_r) + (kap / cc).powi(2) * pi2_n(n_r);
        let i_inv = |rho: f64| {
            let p0 = aa - rho * cc;
            let p1 = bb - rho * dd;
            (dd / p1) * du + ((cc - dd * p0 / p1) / p0) * pi_n(-p1 / p0)
        };
        let scale = self.g_r / (1.0 - e * e).sqrt();
        // T_r = E r^2 + 2 M E r + (E a^2 - a L + 4 M^2 E) + 2 M c(r)/Delta
        let c1 = 4.0 * m * m * e - a * l;
        let c0 = -2.0 * m * e * a * a;
        let mut t_r = e * i_r2 + 2.0 * m * e * i_r + (e * a * a - a * l + 4.0 * m * m * e) * i1;
        for (gam, rho) in lin_over_delta(c1, c0, rp, rm) {
            t_r += 2.0 * m * gam * i_inv(rho);
        }
        // Phi_r = a E + a (2 M E r - a L)/Delta
        let mut ph_r = a * e * i1;
        for (gam, rho) in lin_over_delta(2.0 * m * e, -a * l, rp, rm) {
            ph_r += a * gam * i_inv(rho);
        }
        (scale * t_r, scale * ph_r)
    }

    /// Polar parts `(T_theta, Phi_theta)` of the antiderivatives.
    fn polar_parts(
        &self,
        lambda: f64,
    ) -> (f64, f64) {
        let (a, e, l) = (self.a, self.energy, self.lz);
        if self.z_minus == 0.0 {
            // equatorial: dt/dlambda = a(L - aE) ; dphi/dlambda = L - a E
            return (a * (l - a * e) * lambda, (l - a * e) * lambda);
        }
        let w1 = self.w0 + self.rate_th * lambda;
        let psi1 = Self::amplitude(w1, self.m_th);
        let psi0 = Self::amplitude(self.w0, self.m_th);
        let dw = w1 - self.w0;
        let zm2 = self.z_minus * self.z_minus;
        let pi_part = elliptic_pi(zm2, psi1, self.m_th) - elliptic_pi(zm2, psi0, self.m_th);
        let sn2_part = jacobi_int_sn2(psi1, self.m_th) - jacobi_int_sn2(psi0, self.m_th);
        let inv = 1.0 / self.rate_th;
        let phi_th = inv * (l * pi_part - a * e * dw);
        let t_th = inv * ((a * l - a * a * e) * dw + a * a * e * zm2 * sn2_part);
        (t_th, phi_th)
    }

    /// State `(t, r, theta, phi)` at Mino time `lambda` (`t = 0` and `phi = phi0` at the start).
    #[must_use]
    pub fn state(
        &self,
        lambda: f64,
    ) -> (f64, f64, f64, f64) {
        let (tr, pr) = self.radial_parts(lambda);
        let (tt, pt) = self.polar_parts(lambda);
        let theta = self.cos_theta(lambda).clamp(-1.0, 1.0).acos();
        (tr + tt, self.radius(lambda), theta, self.phi0 + pr + pt)
    }

    /// Fundamental frequencies: `(Upsilon_r, Upsilon_theta, Upsilon_phi, Gamma)` in Mino time
    /// (`Upsilon_theta = 0` for equatorial orbits, where `theta` is constant) and the
    /// coordinate-time frequencies `(Omega_r, Omega_theta, Omega_phi)`.
    #[must_use]
    pub fn frequencies(&self) -> ([f64; 4], [f64; 3]) {
        let lam_r = self.radial_period();
        let (tr, pr) = self.radial_parts(lam_r);
        let (gam_r, phi_r) = (tr / lam_r, pr / lam_r);
        let (gam_th, phi_th, ups_th) = self.polar_period().map_or_else(
            || {
                let (tt, pt) = self.polar_parts(1.0);
                (tt, pt, 0.0)
            },
            |lam_th| {
                let (tt, pt) = self.polar_parts(lam_th);
                (tt / lam_th, pt / lam_th, 2.0 * PI / lam_th)
            },
        );
        let ups_r = 2.0 * PI / lam_r;
        let gamma = gam_r + gam_th;
        let ups_phi = phi_r + phi_th;
        ([ups_r, ups_th, ups_phi, gamma], [ups_r / gamma, ups_th / gamma, ups_phi / gamma])
    }
}

fn solve3_sys(
    a: &[[f64; 3]; 3],
    b: &[f64; 3],
) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-300 {
        return None;
    }
    let mut x = [0.0; 3];
    for c in 0..3 {
        let mut mm = *a;
        for r in 0..3 {
            mm[r][c] = b[r];
        }
        let d = mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1]) - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
            + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]);
        x[c] = d / det;
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// Black-hole shadow
// ---------------------------------------------------------------------------

/// Radius of the circular equatorial photon orbit of Kerr (`prograde` or retrograde).
#[must_use]
pub fn kerr_photon_orbit_radius(
    m: f64,
    a: f64,
    prograde: bool,
) -> f64 {
    let s = if prograde { -1.0 } else { 1.0 };
    2.0 * m * (1.0 + (2.0 / 3.0 * (s * a / m).clamp(-1.0, 1.0).acos()).cos())
}

/// Impact parameters `(xi, eta) = (L/E, Q/E^2)` of the spherical photon orbit of radius `r`.
fn photon_orbit_constants(
    m: f64,
    a: f64,
    r: f64,
) -> (f64, f64) {
    let aa = a * a;
    let xi = (r * r * (3.0 * m - r) - aa * (r + m)) / (a * (r - m));
    let eta = r.powi(3) * (4.0 * a * a * m - r * (r - 3.0 * m).powi(2)) / (a * a * (r - m).powi(2));
    (xi, eta)
}

/// Bardeen's critical curve of the Kerr shadow.
///
/// The boundary of the shadow of a Kerr black hole on the observer's
/// image plane (celestial coordinates `(alpha, beta)` of an observer at infinity and inclination
/// `theta_o`), as `2 n` points (upper half then lower half). For `a = 0` it is the circle of
/// radius `3 sqrt(3) M`.
///
/// # Errors
/// Fails for `|a| >= M`, an observer on the axis or `n < 2`.
pub fn shadow_boundary_analytic(
    m: f64,
    a: f64,
    theta_o: f64,
    n: usize,
) -> Result<Vec<(f64, f64)>, String> {
    if a.abs() >= m || theta_o.sin().abs() < 1e-6 || n < 2 {
        return Err("shadow needs |a| < M, an off-axis observer and n >= 2".into());
    }
    let mut out = Vec::with_capacity(2 * n);
    if a.abs() < 1e-9 {
        let rad = 27.0_f64.sqrt() * m;
        for k in 0..2 * n {
            let ang = 2.0 * PI * k as f64 / (2 * n) as f64;
            out.push((rad * ang.cos(), rad * ang.sin()));
        }
        return Ok(out);
    }
    let (so, co) = theta_o.sin_cos();
    let beta2 = |r: f64| {
        let (xi, eta) = photon_orbit_constants(m, a, r);
        eta + a * a * co * co - xi * xi * co * co / (so * so)
    };
    // beta^2 >= 0 between the two zeros inside [r_pro, r_retro]
    let (r_lo, r_hi) = (kerr_photon_orbit_radius(m, a.abs(), true), kerr_photon_orbit_radius(m, a.abs(), false));
    let samples = 20_000;
    let rs: Vec<f64> = (0..=samples).map(|k| r_lo + (r_hi - r_lo) * f64::from(k) / f64::from(samples)).collect();
    let first = rs.iter().position(|&r| beta2(r) >= 0.0).ok_or("no photon-orbit window for this observer")?;
    let last = rs.iter().rposition(|&r| beta2(r) >= 0.0).ok_or("no photon-orbit window for this observer")?;
    let refine = |mut lo: f64, mut hi: f64| {
        // beta2(lo) < 0 <= beta2(hi) (or the reverse): bisection on the sign change
        let flo = beta2(lo) >= 0.0;
        for _ in 0..200 {
            let mid = f64::midpoint(lo, hi);
            if (beta2(mid) >= 0.0) == flo {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        f64::midpoint(lo, hi)
    };
    let ra = if first > 0 { refine(rs[first - 1], rs[first]) } else { rs[first] };
    let rb = if last + 1 < rs.len() { refine(rs[last], rs[last + 1]) } else { rs[last] };
    let mut upper = Vec::with_capacity(n);
    for k in 0..n {
        let r = ra + (rb - ra) * k as f64 / (n - 1) as f64;
        let (xi, _) = photon_orbit_constants(m, a, r);
        let b = beta2(r).max(0.0).sqrt();
        upper.push((-xi / so, b));
    }
    let mut lower: Vec<(f64, f64)> = upper.iter().rev().map(|&(x, y)| (x, -y)).collect();
    out.append(&mut upper);
    out.append(&mut lower);
    Ok(out)
}

/// Outcome of a backward-traced light ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RayOutcome {
    /// The ray fell through the horizon.
    Captured,
    /// The ray escaped to the observer radius at the given celestial position.
    Escaped {
        /// Polar angle of the arrival direction on the celestial sphere.
        theta: f64,
        /// Azimuth of the arrival direction.
        phi: f64,
    },
}

/// Traces a null geodesic backwards from the observer.
///
/// The impact parameters are `xi = L/E`, `eta = Q/E^2` from the observer at
/// `(r_o, theta_o)` inwards (initial direction: decreasing `r`, polar motion toward larger
/// `theta` if `polar_down`) by integrating the Mino-time equations with RK4.
///
/// # Errors
/// Fails if the ray is not allowed at the observer (`R(r_o) < 0` or `Theta(theta_o) < 0`) or if
/// the step budget is exhausted.
pub fn trace_ray(
    m: f64,
    a: f64,
    theta_o: f64,
    r_o: f64,
    (xi, eta): (f64, f64),
    polar_down: bool,
) -> Result<RayOutcome, String> {
    let rh = m + (m * m - a * a).sqrt();
    let rpot = |r: f64| {
        let rr = r * r + a * a;
        (rr - a * xi).powi(2) - (rr - 2.0 * m * r) * ((xi - a).powi(2) + eta)
    };
    let z0 = theta_o.cos();
    let zp2 = eta - (eta + xi * xi - a * a) * z0 * z0 - a * a * z0.powi(4);
    if rpot(r_o) < 0.0 || zp2 < -1e-12 {
        return Err("ray is not allowed at the observer".into());
    }
    let rhs = |y: &[f64; 5]| -> [f64; 5] {
        let (r, z) = (y[0], y[2]);
        let k = (xi - a).powi(2) + eta;
        let rr = r * r + a * a;
        let rpp = 0.5 * (4.0 * r * (rr - a * xi) - (2.0 * r - 2.0 * m) * k);
        let zpp = -(eta + xi * xi - a * a) * z - 2.0 * a * a * z.powi(3);
        let delta = r * r - 2.0 * m * r + a * a;
        let p = rr - a * xi;
        let phid = a * p / delta - a + xi / (1.0 - z * z).max(1e-12);
        [y[1], rpp, y[3], zpp, phid]
    };
    let mut y = [r_o, -rpot(r_o).sqrt(), z0, if polar_down { -1.0 } else { 1.0 } * zp2.max(0.0).sqrt(), 0.0];
    let mut turned = false;
    for _ in 0..4_000_000 {
        let h = 0.01 * (1.0_f64).min(1.0 / y[0]);
        let k1 = rhs(&y);
        let mut t = y;
        for i in 0..5 {
            t[i] = y[i] + 0.5 * h * k1[i];
        }
        let k2 = rhs(&t);
        for i in 0..5 {
            t[i] = y[i] + 0.5 * h * k2[i];
        }
        let k3 = rhs(&t);
        for i in 0..5 {
            t[i] = y[i] + h * k3[i];
        }
        let k4 = rhs(&t);
        for i in 0..5 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        // project the velocities back onto the first integrals (keeps the large-radius branch stable)
        y[1] = y[1].signum() * rpot(y[0]).max(0.0).sqrt();
        let zz = y[2] * y[2];
        y[3] = y[3].signum() * (eta - (eta + xi * xi - a * a) * zz - a * a * zz * zz).max(0.0).sqrt();
        if y[0] <= rh * (1.0 + 1e-9) {
            return Ok(RayOutcome::Captured);
        }
        if y[1] > 0.0 {
            turned = true;
        }
        if turned && y[0] >= r_o {
            return Ok(RayOutcome::Escaped { theta: y[2].clamp(-1.0, 1.0).acos(), phi: y[4] });
        }
    }
    Err("ray tracing step budget exhausted".into())
}

/// Shadow boundary found by ray tracing.
///
/// For each of `n` directions in the image plane the
/// critical impact radius separating captured from escaping rays is located by bisection.
/// Returns `(alpha, beta)` points of the observer's image plane (observer at inclination
/// `theta_o`, far away at `r_o`).
///
/// # Errors
/// Propagates [`trace_ray`] failures.
pub fn trace_shadow(
    m: f64,
    a: f64,
    theta_o: f64,
    r_o: f64,
    n: usize,
) -> Result<Vec<(f64, f64)>, String> {
    let (so, co) = theta_o.sin_cos();
    let captured = |rho: f64, psi: f64| -> Result<bool, String> {
        let (al, be) = (rho * psi.cos(), rho * psi.sin());
        let xi = -al * so;
        let eta = be * be + (al * al - a * a) * co * co;
        Ok(trace_ray(m, a, theta_o, r_o, (xi, eta), be < 0.0)? == RayOutcome::Captured)
    };
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let psi = 2.0 * PI * (k as f64 + 0.5) / n as f64;
        let (mut lo, mut hi) = (1e-3 * m, 12.0 * m);
        for _ in 0..60 {
            let mid = f64::midpoint(lo, hi);
            if captured(mid, psi)? {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let rho = f64::midpoint(lo, hi);
        out.push((rho * psi.cos(), rho * psi.sin()));
    }
    Ok(out)
}
