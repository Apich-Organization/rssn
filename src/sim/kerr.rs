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
//! | [`KerrSchildGeodesic`] | geodesics in ingoing Kerr-Schild Cartesian coordinates (Hamiltonian form, adaptive Dormand-Prince): regular across the horizon and through the poles/axis; checked against the exact Mino orbits, the Schwarzschild free-fall proper time and conserved `E`, `L_z`, `Q` across the horizon |
//! | [`trace_shadow`] | backward ray tracing of null geodesics (Mino-time integration) to locate the shadow boundary numerically |
//!
//! Limitations: the closed form covers bound timelike orbits (`mu = 1`, `E < 1`) without
//! radiation reaction and excludes the separatrix and circular orbits (`r1 = r2`, where the
//! radial motion degenerates to a point) and extremal spin; ray tracing is numerical.

use std::f64::consts::PI;

use crate::kernels::special_ext::elliptic_f;
use crate::kernels::special_ext::elliptic_k;
use crate::kernels::special_ext::jacobi_elliptic;
use crate::sim::elliptic::elliptic_pi;
use crate::sim::elliptic::jacobi_int_inv_sq;
use crate::sim::elliptic::jacobi_int_sn2;

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
    [
        ((c1 * rp + c0) / (rp - rm), rp),
        ((c1 * rm + c0) / (rm - rp), rm),
    ]
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
        let roots = quartic_real_roots(&c)
            .ok_or("radial polynomial has complex roots: not a bound orbit")?;
        let [r1, r2, r3, r4] = roots;
        let rh = m + (m * m - a * a).sqrt();
        if !(r1 > r2 && r2 > r3 && r3 >= r4 - 1e-12) || r2 <= rh {
            return Err(
                "bound orbit needs r1 > r2 > r3 and r2 above the horizon (plunging or circular)"
                    .into(),
            );
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
        let amp = if radial_out {
            psi0
        } else {
            PI - psi0
        };
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
        Ok(Self {
            m,
            a,
            energy,
            lz,
            carter,
            roots,
            m_r,
            m_th,
            z_minus,
            g_r,
            rate_r,
            rate_th,
            u0,
            w0,
            phi0,
        })
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
        let mut x = [
            (1.0 - (1.0 - e_ecc * e_ecc) * m / p).sqrt(),
            l_guess / (1.0 + q_guess / (l_guess * l_guess + 1e-300)).sqrt(),
            q_guess,
        ];
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
            let dx = solve3_sys(&jac, &f)
                .ok_or("singular Jacobian while solving for the constants of motion")?;
            // damped update
            let mut lam = 1.0;
            for _ in 0..20 {
                let trial = [x[0] - lam * dx[0], x[1] - lam * dx[1], x[2] - lam * dx[2]];
                let ft = resid(&trial);
                if trial[0] > 0.0
                    && trial[0] < 1.0
                    && trial[2] >= 0.0
                    && ft.iter().map(|v| v.abs()).fold(0.0, f64::max) < fn_
                {
                    x = trial;
                    break;
                }
                lam *= 0.5;
            }
        }
        if resid(&x).iter().map(|v| v.abs()).fold(0.0, f64::max) > 1e-9 * (1.0 + r_a.powi(4)) {
            return Err("constants of motion did not converge".into());
        }
        let roots = quartic_real_roots(&radial_poly(m, a, x[0], x[1], x[2]))
            .ok_or("complex radial roots: not a bound orbit")?;
        Self::new(
            m,
            a,
            (x[0], x[1], x[2]),
            (roots[1], std::f64::consts::FRAC_PI_2, 0.0),
            (true, true),
        )
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
        let i_r2 = (bb / dd).powi(2) * du
            + 2.0 * (bb / dd) * kap / cc * pi_n(n_r)
            + (kap / cc).powi(2) * pi2_n(n_r);
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
        (
            [ups_r, ups_th, ups_phi, gamma],
            [ups_r / gamma, ups_th / gamma, ups_phi / gamma],
        )
    }
}

fn solve3_sys(
    a: &[[f64; 3]; 3],
    b: &[f64; 3],
) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
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
        let d = mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
            - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
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
    let (r_lo, r_hi) = (
        kerr_photon_orbit_radius(m, a.abs(), true),
        kerr_photon_orbit_radius(m, a.abs(), false),
    );
    let samples = 20_000;
    let rs: Vec<f64> = (0..=samples)
        .map(|k| r_lo + (r_hi - r_lo) * f64::from(k) / f64::from(samples))
        .collect();
    let first = rs
        .iter()
        .position(|&r| beta2(r) >= 0.0)
        .ok_or("no photon-orbit window for this observer")?;
    let last = rs
        .iter()
        .rposition(|&r| beta2(r) >= 0.0)
        .ok_or("no photon-orbit window for this observer")?;
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
    let ra = if first > 0 {
        refine(rs[first - 1], rs[first])
    } else {
        rs[first]
    };
    let rb = if last + 1 < rs.len() {
        refine(rs[last], rs[last + 1])
    } else {
        rs[last]
    };
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
    let mut y = [
        r_o,
        -rpot(r_o).sqrt(),
        z0,
        if polar_down { -1.0 } else { 1.0 } * zp2.max(0.0).sqrt(),
        0.0,
    ];
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
        y[3] = y[3].signum()
            * (eta - (eta + xi * xi - a * a) * zz - a * a * zz * zz)
                .max(0.0)
                .sqrt();
        if y[0] <= rh * (1.0 + 1e-9) {
            return Ok(RayOutcome::Captured);
        }
        if y[1] > 0.0 {
            turned = true;
        }
        if turned && y[0] >= r_o {
            return Ok(RayOutcome::Escaped {
                theta: y[2].clamp(-1.0, 1.0).acos(),
                phi: y[4],
            });
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

// ---------------------------------------------------------------------------
// Horizon-penetrating geodesics in ingoing Kerr-Schild coordinates
// ---------------------------------------------------------------------------

use crate::sim::hermite::Dual;
use crate::sim::hermite::Scalar;

/// Kerr radius `r(x, y, z)` of the oblate-spheroidal coordinates of the Kerr-Schild form.
fn ks_radius<T: Scalar>(
    a2: f64,
    x: &[T; 3],
) -> T {
    let rr = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - T::from_f64(a2);
    let disc = rr * rr + T::from_f64(4.0 * a2) * x[2] * x[2];
    ((rr + disc.sqrt()) * T::from_f64(0.5)).sqrt()
}

/// Hamiltonian `1/2 g^{mu nu} p_mu p_nu = 1/2 [ -p_t^2 + p^2 - 2 H (-p_t + l . p)^2 ]` of the
/// ingoing Kerr-Schild metric `g = eta + 2 H l l`, `H = M r^3/(r^4 + a^2 z^2)`.
fn ks_hamiltonian<T: Scalar>(
    m: f64,
    a: f64,
    x: &[T; 3],
    p: &[T; 4],
) -> T {
    let a2 = a * a;
    let r = ks_radius(a2, x);
    let r2a2 = r * r + T::from_f64(a2);
    let z = x[2];
    let hh = T::from_f64(m) * r * r * r / (r * r * r * r + T::from_f64(a2) * z * z);
    let lx = (r * x[0] + T::from_f64(a) * x[1]) / r2a2;
    let ly = (r * x[1] - T::from_f64(a) * x[0]) / r2a2;
    let lz = z / r;
    let w = -p[0] + lx * p[1] + ly * p[2] + lz * p[3];
    let p2 = p[1] * p[1] + p[2] * p[2] + p[3] * p[3];
    (-p[0] * p[0] + p2 - T::from_f64(2.0) * hh * w * w) * T::from_f64(0.5)
}

/// Azimuthal shift `psi(r) = int a/Delta dr` between Boyer-Lindquist and ingoing Kerr azimuth
/// (`phi~ = phi + psi`), valid outside the horizon.
fn ks_azimuth_shift(
    m: f64,
    a: f64,
    r: f64,
) -> f64 {
    let d = (m * m - a * a).sqrt();
    if a == 0.0 || d <= 0.0 {
        return 0.0;
    }
    let (rp, rm) = (m + d, m - d);
    a / (rp - rm) * ((r - rp) / (r - rm)).ln()
}

/// Phase-space state of a geodesic in ingoing Kerr-Schild Cartesian coordinates `(t, x, y, z)`
/// with covariant momenta `p_mu`; regular at the horizon, at the poles and on the axis.
#[derive(Clone, Copy, Debug)]
pub struct KerrSchildGeodesic {
    /// Mass.
    pub m: f64,
    /// Spin parameter.
    pub a: f64,
    /// Coordinates `(t, x, y, z)` (`t` is the ingoing Kerr time).
    pub x: [f64; 4],
    /// Covariant momenta `(p_t, p_x, p_y, p_z)`.
    pub p: [f64; 4],
    /// Affine parameter.
    pub tau: f64,
}

impl KerrSchildGeodesic {
    /// Kerr radius of the current position.
    #[must_use]
    pub fn radius(&self) -> f64 {
        ks_radius(self.a * self.a, &[self.x[1], self.x[2], self.x[3]])
    }

    /// Mass-shell value `1/2 g^{mu nu} p_mu p_nu` (equals `-mu^2/2`).
    #[must_use]
    pub fn hamiltonian(&self) -> f64 {
        ks_hamiltonian(self.m, self.a, &[self.x[1], self.x[2], self.x[3]], &self.p)
    }

    /// Boyer-Lindquist `(r, theta, phi)` of the position (outside the horizon).
    #[must_use]
    pub fn boyer_lindquist(&self) -> (f64, f64, f64) {
        let r = self.radius();
        let phit = self.x[2].atan2(self.x[1]) - self.a.atan2(r);
        let phi = phit - ks_azimuth_shift(self.m, self.a, r);
        (r, (self.x[3] / r).clamp(-1.0, 1.0).acos(), phi)
    }

    /// Conserved energy `E = -p_t`.
    #[must_use]
    pub fn energy(&self) -> f64 {
        -self.p[0]
    }

    /// Conserved axial angular momentum `L_z = x p_y - y p_x`.
    #[must_use]
    pub fn angular_momentum(&self) -> f64 {
        self.x[1] * self.p[2] - self.x[2] * self.p[1]
    }

    fn derivative(
        &self,
        y: &[f64; 8],
    ) -> [f64; 8] {
        let (m, a) = (self.m, self.a);
        let x = [y[1], y[2], y[3]];
        let p = [y[4], y[5], y[6], y[7]];
        let a2 = a * a;
        let r = ks_radius(a2, &x);
        let r2a2 = r * r + a2;
        let hh = m * r * r * r / (r * r * r * r + a2 * x[2] * x[2]);
        let l = [
            (r * x[0] + a * x[1]) / r2a2,
            (r * x[1] - a * x[0]) / r2a2,
            x[2] / r,
        ];
        let w = -p[0] + l[0] * p[1] + l[1] * p[2] + l[2] * p[3];
        let mut d = [0.0; 8];
        d[0] = -p[0] + 2.0 * hh * w;
        for i in 0..3 {
            d[1 + i] = p[1 + i] - 2.0 * hh * l[i] * w;
        }
        for c in 0..3 {
            let xd: [Dual; 3] = std::array::from_fn(|k| Dual {
                v: x[k],
                d: if k == c { 1.0 } else { 0.0 },
            });
            let pd: [Dual; 4] = std::array::from_fn(|k| Dual { v: p[k], d: 0.0 });
            d[5 + c] = -ks_hamiltonian(m, a, &xd, &pd).d;
        }
        d
    }

    /// Advances the geodesic by the affine parameter `dtau_total` with an adaptive
    /// Dormand-Prince 5(4) integrator (relative tolerance `tol`). Returns `false` if the
    /// radius fell below `r_stop` (the state is left at the last accepted step).
    pub fn advance(
        &mut self,
        dtau_total: f64,
        tol: f64,
        r_stop: f64,
    ) -> bool {
        let a: [[f64; 6]; 6] = [
            [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
            [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
            [
                19372.0 / 6561.0,
                -25360.0 / 2187.0,
                64448.0 / 6561.0,
                -212.0 / 729.0,
                0.0,
                0.0,
            ],
            [
                9017.0 / 3168.0,
                -355.0 / 33.0,
                46732.0 / 5247.0,
                49.0 / 176.0,
                -5103.0 / 18656.0,
                0.0,
            ],
            [
                35.0 / 384.0,
                0.0,
                500.0 / 1113.0,
                125.0 / 192.0,
                -2187.0 / 6784.0,
                11.0 / 84.0,
            ],
        ];
        let b5 = [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
            0.0,
        ];
        let b4 = [
            5179.0 / 57600.0,
            0.0,
            7571.0 / 16695.0,
            393.0 / 640.0,
            -92097.0 / 339200.0,
            187.0 / 2100.0,
            1.0 / 40.0,
        ];
        let mut y: [f64; 8] = [
            self.x[0], self.x[1], self.x[2], self.x[3], self.p[0], self.p[1], self.p[2], self.p[3],
        ];
        let mut t = 0.0;
        let mut h = dtau_total.signum() * 0.01_f64.min(dtau_total.abs());
        let mut steps = 0usize;
        while (dtau_total - t) * dtau_total.signum() > 1e-14 * dtau_total.abs() {
            steps += 1;
            if steps > 5_000_000 {
                break;
            }
            if (t + h - dtau_total) * dtau_total.signum() > 0.0 {
                h = dtau_total - t;
            }
            let mut k = [[0.0; 8]; 7];
            k[0] = self.derivative(&y);
            for s in 1..7 {
                let mut ys = y;
                for j in 0..s {
                    let coef = if s < 6 { a[s - 1][j] } else { b5[j] };
                    for c in 0..8 {
                        ys[c] += h * coef * k[j][c];
                    }
                }
                k[s] = self.derivative(&ys);
            }
            let mut y5 = y;
            let mut err = 0.0_f64;
            for c in 0..8 {
                let (mut d5, mut d4) = (0.0, 0.0);
                for j in 0..7 {
                    d5 += b5[j] * k[j][c];
                    d4 += b4[j] * k[j][c];
                }
                y5[c] += h * d5;
                let sc = tol * (1.0 + y[c].abs().max(y5[c].abs()));
                err = err.max((h * (d5 - d4)).abs() / sc);
            }
            if err <= 1.0 || h.abs() < 1e-12 {
                y = y5;
                t += h;
                self.tau += h;
                let rr = ks_radius(self.a * self.a, &[y[1], y[2], y[3]]);
                if rr < r_stop {
                    self.x = [y[0], y[1], y[2], y[3]];
                    self.p = [y[4], y[5], y[6], y[7]];
                    return false;
                }
            }
            let fac = if err == 0.0 {
                5.0
            } else {
                (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
            };
            h *= fac;
        }
        self.x = [y[0], y[1], y[2], y[3]];
        self.p = [y[4], y[5], y[6], y[7]];
        true
    }

    /// Starts a geodesic from Boyer-Lindquist data: position `(r, theta, phi)`, constants of
    /// motion `(E, L, Q)`, mass-shell value `mu2` (1 timelike, 0 null) and the signs of `dr/d tau`
    /// and `d theta / d tau`.
    ///
    /// # Errors
    /// Fails if the position is inside the horizon, on the axis, or the motion is forbidden there.
    pub fn from_constants(
        m: f64,
        a: f64,
        (energy, lz, carter): (f64, f64, f64),
        (r, theta, phi): (f64, f64, f64),
        mu2: f64,
        (radial_out, polar_down): (bool, bool),
    ) -> Result<Self, String> {
        let delta = r * r - 2.0 * m * r + a * a;
        let (st, ct) = theta.sin_cos();
        if delta <= 0.0 || st.abs() < 1e-9 {
            return Err(
                "start outside the horizon and away from the axis (use Cartesian data otherwise)"
                    .into(),
            );
        }
        let rpot = (energy * (r * r + a * a) - a * lz).powi(2)
            - delta * (mu2 * r * r + (lz - a * energy).powi(2) + carter);
        let tpot = carter - ct * ct * (a * a * (mu2 - energy * energy) + lz * lz / (st * st));
        let scale = (energy * (r * r + a * a)).powi(2).max(1.0);
        if rpot < -1e-9 * scale || tpot < -1e-9 {
            return Err("motion is forbidden at the requested position".into());
        }
        let pr = (if radial_out { 1.0 } else { -1.0 })
            * (if rpot < 1e-13 * scale {
                0.0
            } else {
                rpot
            })
            .sqrt()
            / delta;
        let pth = (if polar_down { 1.0 } else { -1.0 }) * tpot.max(0.0).sqrt();
        let p_bl = [-energy, pr, pth, lz];
        // BL -> ingoing KS Cartesian
        let phit = phi + ks_azimuth_shift(m, a, r);
        let (sp, cp) = phit.sin_cos();
        let x = st * (r * cp - a * sp);
        let y = st * (r * sp + a * cp);
        let z = r * ct;
        let jac = [
            [1.0, 0.0, 0.0, 0.0],
            [
                2.0 * m * r / delta,
                st * cp - y * a / delta,
                st * sp + x * a / delta,
                ct,
            ],
            [0.0, ct * (r * cp - a * sp), ct * (r * sp + a * cp), -r * st],
            [0.0, -y, x, 0.0],
        ];
        // p_bl[alpha] = sum_mu jac[alpha][mu] p_ks[mu]
        let p_ks = solve4_lin(&jac, &p_bl).ok_or("singular coordinate transformation")?;
        Ok(Self {
            m,
            a,
            x: [0.0, x, y, z],
            p: p_ks,
            tau: 0.0,
        })
    }

    /// Carter constant `Q = p_theta^2 + cos^2 theta [a^2 (mu^2 - E^2) + L^2/sin^2 theta]` of the
    /// current state (off the axis).
    #[must_use]
    pub fn carter_constant(&self) -> f64 {
        let (m, a) = (self.m, self.a);
        let r = self.radius();
        let (x, y, z) = (self.x[1], self.x[2], self.x[3]);
        let ct = z / r;
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        let delta = r * r - 2.0 * m * r + a * a;
        // phi~ from x + i y = (r + i a) sin(theta) exp(i phi~)
        let phit = y.atan2(x) - a.atan2(r);
        let (sp, cp) = phit.sin_cos();
        // d x^mu / d theta (BL) is the third row of the Jacobian
        let jth = [0.0, ct * (r * cp - a * sp), ct * (r * sp + a * cp), -r * st];
        let pth: f64 = (0..4).map(|k| jth[k] * self.p[k]).sum();
        let e = self.energy();
        let l = self.angular_momentum();
        let mu2 = -2.0 * self.hamiltonian();
        let _ = delta;
        pth * pth + ct * ct * (a * a * (mu2 - e * e) + l * l / (st * st))
    }
}

/// Solves `j^T x = b`-type system `sum_mu jac[alpha][mu] x[mu] = b[alpha]` by Gaussian elimination.
fn solve4_lin(
    jac: &[[f64; 4]; 4],
    b: &[f64; 4],
) -> Option<[f64; 4]> {
    let mut a = *jac;
    let mut x = *b;
    for c in 0..4 {
        let mut piv = c;
        for r in c + 1..4 {
            if a[r][c].abs() > a[piv][c].abs() {
                piv = r;
            }
        }
        if a[piv][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, piv);
        x.swap(c, piv);
        for r in c + 1..4 {
            let f = a[r][c] / a[c][c];
            for k in c..4 {
                a[r][k] -= f * a[c][k];
            }
            x[r] -= f * x[c];
        }
    }
    for c in (0..4).rev() {
        for k in c + 1..4 {
            x[c] -= a[c][k] * x[k];
        }
        x[c] /= a[c][c];
    }
    Some(x)
}
