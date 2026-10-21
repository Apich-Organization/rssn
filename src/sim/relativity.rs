//! # Kerr-Newman geodesics (Schwarzschild, Kerr, Reissner-Nordstrom)
//!
//! Timelike and null geodesics of the Kerr-Newman family in Boyer-Lindquist
//! coordinates `(t, r, theta, phi)` from the Hamiltonian
//! `H = 1/2 g^{mu nu} p_mu p_nu = -mu^2/2`. `t` and `phi` are cyclic, so
//! `E = -p_t` and `L = p_phi` are conserved exactly; the reduced system
//! `(r, theta, p_r, p_theta)` is advanced by the implicit midpoint rule
//! (symplectic, symmetric, second order) with Newton iteration, followed by
//! an optional projection of `(p_r, p_theta)` onto the mass-shell constraint.
//!
//! | Capability | Notes |
//! |---|---|
//! | Spacetimes | Schwarzschild (`a = Q = 0`), Reissner-Nordstrom (`a = 0`), Kerr (`Q = 0`), Kerr-Newman |
//! | Orbits | arbitrary timelike/null geodesics, circular equatorial orbits ([`KerrNewman::circular_orbit`]) |
//! | Characteristic radii | [`KerrNewman::isco`], [`KerrNewman::photon_orbit`], [`KerrNewman::horizon`] |
//! | Events | horizon-crossing (`r <= r_stop`, located by bisection on the step), escape (`r >= r_max`) |
//! | Dimensions | 3+1 (a single physical spacetime; this module is not dimension generic) |
//!
//! Closed-form bound orbits (Mino time, Jacobi and elliptic integrals), shadow ray tracing and
//! horizon- and pole-penetrating geodesics in Kerr-Schild coordinates are in [`crate::sim::kerr`].
//!
//! Limitations: Boyer-Lindquist coordinates are singular at the horizon, so
//! integration stops at a user radius slightly outside it; the poles
//! `theta = 0, pi` are excluded.

/// Kerr-Newman spacetime with mass `m`, spin parameter `a` (`J / M`) and charge `q` (geometric units `G = c = 1`).
#[derive(Clone, Copy, Debug)]
pub struct KerrNewman {
    /// Mass.
    pub m: f64,
    /// Spin parameter.
    pub a: f64,
    /// Charge.
    pub q: f64,
}

/// Why an integration stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Termination {
    /// Requested number of steps completed.
    Steps,
    /// `r` fell to the stop radius (horizon crossing).
    Horizon,
    /// `r` exceeded the escape radius.
    Escaped,
}

/// Phase-space state of a geodesic.
#[derive(Clone, Copy, Debug)]
pub struct Geodesic {
    /// Spacetime.
    pub metric: KerrNewman,
    /// Coordinate time.
    pub t: f64,
    /// Radius.
    pub r: f64,
    /// Polar angle.
    pub theta: f64,
    /// Azimuth.
    pub phi: f64,
    /// Radial momentum `p_r`.
    pub pr: f64,
    /// Polar momentum `p_theta`.
    pub ptheta: f64,
    /// Energy `E = -p_t`.
    pub energy: f64,
    /// Angular momentum `L = p_phi`.
    pub lz: f64,
    /// Affine parameter.
    pub lambda: f64,
    /// Mass-shell value `mu^2` (1 timelike, 0 null).
    pub mu2: f64,
}

impl KerrNewman {
    /// Schwarzschild black hole.
    #[must_use]
    pub const fn schwarzschild(m: f64) -> Self {
        Self { m, a: 0.0, q: 0.0 }
    }

    /// Outer horizon radius `M + sqrt(M^2 - a^2 - Q^2)` (`None` for a naked singularity).
    #[must_use]
    pub fn horizon(&self) -> Option<f64> {
        let d = self.m * self.m - self.a * self.a - self.q * self.q;
        (d >= 0.0).then(|| self.m + d.sqrt())
    }

    fn delta(
        &self,
        r: f64,
    ) -> f64 {
        r * r - 2.0 * self.m * r + self.a * self.a + self.q * self.q
    }

    /// Hamiltonian `1/2 g^{mu nu} p p` at `(r, theta)` with momenta `(p_r, p_theta)`, energy and angular momentum.
    #[must_use]
    pub fn hamiltonian(
        &self,
        r: f64,
        th: f64,
        pr: f64,
        pth: f64,
        e: f64,
        l: f64,
    ) -> f64 {
        let a = self.a;
        let (s, c) = th.sin_cos();
        let sig = r * r + a * a * c * c;
        let del = self.delta(r);
        let w = 2.0 * self.m * r - self.q * self.q;
        let gtt = -((r * r + a * a).powi(2) - del * a * a * s * s) / (sig * del);
        let gtp = -a * w / (sig * del);
        let gpp = (del - a * a * s * s) / (sig * del * s * s);
        0.5 * (gtt * e * e - 2.0 * gtp * e * l + gpp * l * l + del / sig * pr * pr + pth * pth / sig)
    }

    /// Hamilton's equations for `(r, theta, p_r, p_theta)` plus `(dt, dphi)`: returns `([4], [2])`.
    fn rhs(
        &self,
        y: &[f64; 4],
        e: f64,
        l: f64,
    ) -> ([f64; 4], [f64; 2]) {
        let (r, th, pr, pth) = (y[0], y[1], y[2], y[3]);
        let a = self.a;
        let (s, c) = th.sin_cos();
        let sig = r * r + a * a * c * c;
        let del = self.delta(r);
        let w = 2.0 * self.m * r - self.q * self.q;
        let gtt = -((r * r + a * a).powi(2) - del * a * a * s * s) / (sig * del);
        let gtp = -a * w / (sig * del);
        let gpp = (del - a * a * s * s) / (sig * del * s * s);
        let ham = |rr: f64, tt: f64| self.hamiltonian(rr, tt, pr, pth, e, l);
        let dr = 3e-4 * r.abs().max(1.0);
        let dt = 1e-3;
        let d5 = |f: &dyn Fn(f64) -> f64, x: f64, h: f64| (-f(x + 2.0 * h) + 8.0 * f(x + h) - 8.0 * f(x - h) + f(x - 2.0 * h)) / (12.0 * h);
        let dh_dr = d5(&|x| ham(x, th), r, dr);
        let dh_dth = d5(&|x| ham(r, x), th, dt);
        (
            [del / sig * pr, pth / sig, -dh_dr, -dh_dth],
            [-gtt * e + gtp * l, -gtp * e + gpp * l],
        )
    }

    /// Circular equatorial orbit `(E, L, Omega)` at radius `r` (`None` if none exists).
    #[must_use]
    pub fn circular_orbit(
        &self,
        r: f64,
        prograde: bool,
    ) -> Option<(f64, f64, f64)> {
        let (m, a, q) = (self.m, self.a, self.q);
        let w = 2.0 * m / r - q * q / (r * r);
        let wp = -2.0 * m / (r * r) + 2.0 * q * q / (r * r * r);
        let gtt = -(1.0 - w);
        let gtp = -a * w;
        let gpp = r * r + a * a + a * a * w;
        let (gttp, gtpp, gppp) = (wp, -a * wp, 2.0 * r + a * a * wp);
        let (b2, ac) = (gtpp * gtpp, gttp * gppp);
        let disc = b2 - ac;
        if disc < 0.0 {
            return None;
        }
        let sgn = if prograde { 1.0 } else { -1.0 };
        let omega = (-gtpp + sgn * disc.sqrt()) / gppp;
        let norm = -(gtt + 2.0 * gtp * omega + gpp * omega * omega);
        if norm <= 0.0 {
            return None;
        }
        let ut = 1.0 / norm.sqrt();
        Some((-(gtt + gtp * omega) * ut, (gtp + gpp * omega) * ut, omega))
    }

    /// Radius of the circular photon orbit (prograde or retrograde), by bisection.
    #[must_use]
    pub fn photon_orbit(
        &self,
        prograde: bool,
    ) -> Option<f64> {
        let (m, a, q) = (self.m, self.a, self.q);
        // g_tt + 2 g_tphi Omega + g_phiphi Omega^2 = 0 at the circular-orbit Omega
        let f = |r: f64| -> f64 {
            let w = 2.0 * m / r - q * q / (r * r);
            let wp = -2.0 * m / (r * r) + 2.0 * q * q / (r * r * r);
            let (gtt, gtp, gpp) = (-(1.0 - w), -a * w, r * r + a * a + a * a * w);
            let (gttp, gtpp, gppp) = (wp, -a * wp, 2.0 * r + a * a * wp);
            let (b2, ac) = (gtpp * gtpp, gttp * gppp);
            let disc = b2 - ac;
            let sgn = if prograde { 1.0 } else { -1.0 };
            let om = (-gtpp + sgn * disc.max(0.0).sqrt()) / gppp;
            gtt + 2.0 * gtp * om + gpp * om * om
        };
        let rh = self.horizon()?;
        let (mut lo, mut hi) = (rh * 1.0000001 + 1e-9, 10.0 * m + 10.0 * (a.abs() + q.abs()));
        // timelike orbits (f<0) outside, spacelike inside: find sign change from + to -
        let n = 4000;
        let mut prev = lo;
        let mut found = false;
        for k in 1..=n {
            let x = lo + (hi - lo) * f64::from(k) / f64::from(n);
            if f(prev) > 0.0 && f(x) <= 0.0 {
                lo = prev;
                hi = x;
                found = true;
                break;
            }
            prev = x;
        }
        if !found {
            return None;
        }
        for _ in 0..200 {
            let mid = f64::midpoint(lo, hi);
            if f(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(f64::midpoint(lo, hi))
    }

    /// Innermost stable circular orbit: the minimum of `E(r)` over timelike circular orbits.
    #[must_use]
    pub fn isco(
        &self,
        prograde: bool,
    ) -> Option<f64> {
        let rp = self.photon_orbit(prograde)?;
        let e = |r: f64| self.circular_orbit(r, prograde).map_or(f64::INFINITY, |o| o.0);
        let mut lo = rp * 1.0001;
        let mut hi = 30.0 * self.m;
        // golden section on a unimodal E(r)
        let g = 0.618_033_988_749_894_9;
        for _ in 0..200 {
            let c = hi - g * (hi - lo);
            let d = lo + g * (hi - lo);
            if e(c) < e(d) {
                hi = d;
            } else {
                lo = c;
            }
        }
        Some(f64::midpoint(lo, hi))
    }
}

/// Initial data of a geodesic: position, conserved `E`, `L`, polar momentum, radial direction and `mu^2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodesicInit {
    /// Boyer-Lindquist radius.
    pub r: f64,
    /// Polar angle.
    pub theta: f64,
    /// Azimuth.
    pub phi: f64,
    /// Conserved energy.
    pub energy: f64,
    /// Conserved axial angular momentum.
    pub lz: f64,
    /// Polar momentum.
    pub ptheta: f64,
    /// Sign (`+-1`) of the radial momentum.
    pub radial_sign: f64,
    /// Mass-shell value `mu^2` (1 for massive, 0 for null).
    pub mu2: f64,
}

impl Geodesic {
    /// Starts a geodesic at `(r, theta, phi)` with conserved `E`, `L`, polar momentum and radial direction
    /// (`radial_sign = +-1`), solving the mass-shell constraint for `p_r`.
    ///
    /// # Errors
    /// Fails when the constraint has no real solution (forbidden region).
    #[allow(clippy::too_many_arguments)] // kept for API compatibility; prefer `from_init`
    pub fn new(
        metric: KerrNewman,
        r: f64,
        theta: f64,
        phi: f64,
        energy: f64,
        lz: f64,
        ptheta: f64,
        radial_sign: f64,
        mu2: f64,
    ) -> Result<Self, String> {
        Self::from_init(metric, &GeodesicInit { r, theta, phi, energy, lz, ptheta, radial_sign, mu2 })
    }

    /// As [`Geodesic::new`] with the initial data bundled in a [`GeodesicInit`].
    ///
    /// # Errors
    /// Fails when the constraint has no real solution (forbidden region).
    pub fn from_init(
        metric: KerrNewman,
        init: &GeodesicInit,
    ) -> Result<Self, String> {
        let GeodesicInit { r, theta, phi, energy, lz, ptheta, radial_sign, mu2 } = *init;
        let mut g = Self { metric, t: 0.0, r, theta, phi, pr: 0.0, ptheta, energy, lz, lambda: 0.0, mu2 };
        let h0 = metric.hamiltonian(r, theta, 0.0, ptheta, energy, lz);
        let kcoef = metric.hamiltonian(r, theta, 1.0, ptheta, energy, lz) - h0;
        let need = -0.5 * mu2 - h0;
        if need < 0.0 || kcoef <= 0.0 {
            return Err("mass-shell constraint has no real radial momentum here".into());
        }
        g.pr = radial_sign.signum() * (need / kcoef).sqrt();
        Ok(g)
    }

    /// Circular equatorial timelike orbit at radius `r`.
    ///
    /// # Errors
    /// Fails if no circular orbit exists at `r`.
    pub fn circular(
        metric: KerrNewman,
        r: f64,
        prograde: bool,
    ) -> Result<Self, String> {
        let (e, l, _) = metric.circular_orbit(r, prograde).ok_or("no circular orbit at this radius")?;
        Ok(Self { metric, t: 0.0, r, theta: std::f64::consts::FRAC_PI_2, phi: 0.0, pr: 0.0, ptheta: 0.0, energy: e, lz: l, lambda: 0.0, mu2: 1.0 })
    }

    /// Value of the constraint residual `H + mu^2/2` (zero on the mass shell).
    #[must_use]
    pub fn constraint(&self) -> f64 {
        self.metric.hamiltonian(self.r, self.theta, self.pr, self.ptheta, self.energy, self.lz) + 0.5 * self.mu2
    }

    fn project(&mut self) {
        let m = &self.metric;
        let v = m.hamiltonian(self.r, self.theta, 0.0, 0.0, self.energy, self.lz);
        let kin = m.hamiltonian(self.r, self.theta, self.pr, self.ptheta, self.energy, self.lz) - v;
        let target = -0.5 * self.mu2 - v;
        if kin > 1e-14 && target > 0.0 {
            let s = (target / kin).sqrt();
            self.pr *= s;
            self.ptheta *= s;
        }
    }

    fn midpoint(
        &self,
        h: f64,
    ) -> Result<([f64; 4], [f64; 2]), String> {
        let y0 = [self.r, self.theta, self.pr, self.ptheta];
        let m = &self.metric;
        let (e, l) = (self.energy, self.lz);
        let resid = |y1: &[f64; 4]| -> [f64; 4] {
            let mid = [f64::midpoint(y0[0], y1[0]), f64::midpoint(y0[1], y1[1]), f64::midpoint(y0[2], y1[2]), f64::midpoint(y0[3], y1[3])];
            let (f, _) = m.rhs(&mid, e, l);
            [y1[0] - y0[0] - h * f[0], y1[1] - y0[1] - h * f[1], y1[2] - y0[2] - h * f[2], y1[3] - y0[3] - h * f[3]]
        };
        let (f0, _) = m.rhs(&y0, e, l);
        let mut y1 = [y0[0] + h * f0[0], y0[1] + h * f0[1], y0[2] + h * f0[2], y0[3] + h * f0[3]];
        let mut prev = f64::INFINITY;
        for it in 0..60 {
            let f = resid(&y1);
            let scale = y1.iter().map(|v| v.abs()).fold(1.0, f64::max);
            let res = f.iter().map(|v| v.abs()).fold(0.0, f64::max) / scale;
            // the finite-difference gradient limits the attainable residual to ~1e-11
            if res < 1e-10 || (it >= 2 && res < 1e-7 && res >= 0.5 * prev) {
                let mid = [f64::midpoint(y0[0], y1[0]), f64::midpoint(y0[1], y1[1]), f64::midpoint(y0[2], y1[2]), f64::midpoint(y0[3], y1[3])];
                let (_, tp) = m.rhs(&mid, e, l);
                return Ok((y1, tp));
            }
            prev = res;
            let mut jac = [[0.0; 4]; 4];
            for k in 0..4 {
                let eps = 1e-7 * y1[k].abs().max(1.0);
                let mut yp = y1;
                yp[k] += eps;
                let fp = resid(&yp);
                for i in 0..4 {
                    jac[i][k] = (fp[i] - f[i]) / eps;
                }
            }
            let dy = solve4(jac, f).ok_or("singular Jacobian in implicit midpoint step")?;
            for k in 0..4 {
                y1[k] -= dy[k];
            }
        }
        Err("implicit midpoint iteration did not converge".into())
    }

    /// One implicit-midpoint step of affine length `h` (optionally projecting onto the mass shell).
    ///
    /// # Errors
    /// Fails if the Newton iteration does not converge.
    pub fn step(
        &mut self,
        h: f64,
        project: bool,
    ) -> Result<(), String> {
        let (y1, tp) = self.midpoint(h)?;
        self.r = y1[0];
        self.theta = y1[1];
        self.pr = y1[2];
        self.ptheta = y1[3];
        self.t += h * tp[0];
        self.phi += h * tp[1];
        self.lambda += h;
        if project {
            self.project();
        }
        Ok(())
    }

    /// Integrates until `n` steps are done, `r <= r_stop` (horizon) or `r >= r_max`; the horizon event is
    /// located to `1e-12` by bisecting the length of the last step.
    ///
    /// # Errors
    /// Propagates step failures.
    pub fn integrate(
        &mut self,
        h: f64,
        n: usize,
        project: bool,
        r_stop: f64,
        r_max: f64,
    ) -> Result<Termination, String> {
        for _ in 0..n {
            let before = *self;
            self.step(h, project)?;
            if self.r <= r_stop {
                let (mut lo, mut hi) = (0.0, h);
                let mut best = *self;
                for _ in 0..50 {
                    let mid = f64::midpoint(lo, hi);
                    let mut trial = before;
                    trial.step(mid, project)?;
                    if trial.r <= r_stop {
                        hi = mid;
                        best = trial;
                    } else {
                        lo = mid;
                    }
                }
                *self = best;
                return Ok(Termination::Horizon);
            }
            if self.r >= r_max {
                return Ok(Termination::Escaped);
            }
        }
        Ok(Termination::Steps)
    }
}

fn solve4(
    mut a: [[f64; 4]; 4],
    b: [f64; 4],
) -> Option<[f64; 4]> {
    let mut x = b;
    for c in 0..4 {
        let mut p = c;
        for i in c + 1..4 {
            if a[i][c].abs() > a[p][c].abs() {
                p = i;
            }
        }
        if a[p][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, p);
        x.swap(c, p);
        for i in c + 1..4 {
            let f = a[i][c] / a[c][c];
            for k in c..4 {
                a[i][k] -= f * a[c][k];
            }
            x[i] -= f * x[c];
        }
    }
    for i in (0..4).rev() {
        for k in i + 1..4 {
            x[i] -= a[i][k] * x[k];
        }
        x[i] /= a[i][i];
    }
    Some(x)
}
