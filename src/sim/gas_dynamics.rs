//! Compressible gas dynamics (Euler equations, HLLC, MUSCL).
//!
//! # Compressible gas dynamics
//!
//! The one-dimensional Euler equations of an ideal gas,
//! `∂_t (ρ, ρu, E) + ∂_x (ρu, ρu² + p, u(E + p)) = 0` with
//! `p = (γ - 1)(E - ρu²/2)`:
//!
//! * a finite-volume solver with the **HLLC** approximate Riemann solver
//!   (Toro, Spruce & Speares), **MUSCL** reconstruction of the primitive
//!   variables with a minmod or van Leer limiter, and the third-order
//!   strong-stability-preserving Runge–Kutta method (Shu–Osher), the time
//!   step from a CFL condition; transmissive or reflective boundaries;
//! * the **exact Riemann solver** (Newton iteration for the star pressure,
//!   sampling of the self-similar solution), used for validation and for
//!   shock-tube reference solutions.
//!
//! # Coverage audit
//!
//! | Equation | 1D | 2D | 3D | Boundaries |
//! |---|---|---|---|---|
//! | Euler, HLLC + MUSCL (minmod, van Leer) | yes | no | no | transmissive, reflective |
//!
//! [`crate::sim::euler_nd`] extends this to 2D/3D, WENO5, Roe/Rusanov fluxes,
//! periodic and fixed-state faces and shallow water.

/// Primitive state `(ρ, u, p)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Primitive {
    /// Density.
    pub rho: f64,
    /// Velocity.
    pub u: f64,
    /// Pressure.
    pub p: f64,
}

/// Slope limiter for the MUSCL reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Limiter {
    /// First order (no reconstruction).
    None,
    /// Minmod (most dissipative TVD limiter).
    Minmod,
    /// Van Leer (smooth harmonic-mean limiter).
    VanLeer,
}

/// Boundary treatment at both ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    /// Zero-gradient (outflow).
    Transmissive,
    /// Solid wall (velocity mirrored).
    Reflective,
}

type Conserved = [f64; 3];

fn to_conserved(
    w: Primitive,
    gamma: f64,
) -> Conserved {
    [
        w.rho,
        w.rho * w.u,
        w.p / (gamma - 1.0) + 0.5 * w.rho * w.u * w.u,
    ]
}

fn to_primitive(
    q: Conserved,
    gamma: f64,
) -> Primitive {
    let rho = q[0];
    let u = q[1] / rho;
    Primitive {
        rho,
        u,
        p: (gamma - 1.0) * (q[2] - 0.5 * rho * u * u),
    }
}

fn physical_flux(
    w: Primitive,
    gamma: f64,
) -> Conserved {
    let e = w.p / (gamma - 1.0) + 0.5 * w.rho * w.u * w.u;
    [w.rho * w.u, w.rho * w.u * w.u + w.p, w.u * (e + w.p)]
}

fn sound_speed(
    w: Primitive,
    gamma: f64,
) -> f64 {
    (gamma * w.p / w.rho).sqrt()
}

/// The HLLC numerical flux between the states `l` and `r`.
#[must_use]
pub fn hllc_flux(
    l: Primitive,
    r: Primitive,
    gamma: f64,
) -> [f64; 3] {
    let (cl, cr) = (sound_speed(l, gamma), sound_speed(r, gamma));
    // Davis wave-speed estimates.
    let sl = (l.u - cl).min(r.u - cr);
    let sr = (l.u + cl).max(r.u + cr);
    let (ql, qr) = (to_conserved(l, gamma), to_conserved(r, gamma));
    let (fl, fr) = (physical_flux(l, gamma), physical_flux(r, gamma));
    if sl >= 0.0 {
        return fl;
    }
    if sr <= 0.0 {
        return fr;
    }
    let s_star = (r.p - l.p + l.rho * l.u * (sl - l.u) - r.rho * r.u * (sr - r.u))
        / (l.rho * (sl - l.u) - r.rho * (sr - r.u));
    let star = |w: Primitive, q: Conserved, s: f64| -> Conserved {
        let factor = w.rho * (s - w.u) / (s - s_star);
        let energy = q[2] / w.rho + (s_star - w.u) * (s_star + w.p / (w.rho * (s - w.u)));
        [factor, factor * s_star, factor * energy]
    };
    if s_star >= 0.0 {
        let qs = star(l, ql, sl);
        [
            fl[0] + sl * (qs[0] - ql[0]),
            fl[1] + sl * (qs[1] - ql[1]),
            fl[2] + sl * (qs[2] - ql[2]),
        ]
    } else {
        let qs = star(r, qr, sr);
        [
            fr[0] + sr * (qs[0] - qr[0]),
            fr[1] + sr * (qs[1] - qr[1]),
            fr[2] + sr * (qs[2] - qr[2]),
        ]
    }
}

fn limited_slope(
    limiter: Limiter,
    a: f64,
    b: f64,
) -> f64 {
    match limiter {
        | Limiter::None => 0.0,
        | Limiter::Minmod => {
            if a * b <= 0.0 {
                0.0
            } else if a.abs() < b.abs() {
                a
            } else {
                b
            }
        },
        | Limiter::VanLeer => {
            if a * b <= 0.0 {
                0.0
            } else {
                2.0 * a * b / (a + b)
            }
        },
    }
}

/// The right-hand side `-∂_x F` of the semi-discrete scheme.
fn residual(
    q: &[Conserved],
    dx: f64,
    gamma: f64,
    limiter: Limiter,
    boundary: Boundary,
) -> Vec<Conserved> {
    let n = q.len();
    // Two ghost cells on each side.
    let ghost = |w: Primitive| {
        match boundary {
            | Boundary::Transmissive => w,
            | Boundary::Reflective => Primitive { u: -w.u, ..w },
        }
    };
    let mut w: Vec<Primitive> = Vec::with_capacity(n + 4);
    let prim: Vec<Primitive> = q.iter().map(|&c| to_primitive(c, gamma)).collect();
    w.push(ghost(prim[1.min(n - 1)]));
    w.push(ghost(prim[0]));
    w.extend_from_slice(&prim);
    w.push(ghost(prim[n - 1]));
    w.push(ghost(prim[n.saturating_sub(2)]));
    // Limited slopes of (ρ, u, p) in every cell, ghosts included.
    let slope = |i: usize| -> [f64; 3] {
        let (a, b, c) = (w[i - 1], w[i], w[i + 1]);
        [
            limited_slope(limiter, b.rho - a.rho, c.rho - b.rho),
            limited_slope(limiter, b.u - a.u, c.u - b.u),
            limited_slope(limiter, b.p - a.p, c.p - b.p),
        ]
    };
    let face = |i: usize, side: f64| -> Primitive {
        let s = slope(i);
        let w = w[i];
        let candidate = Primitive {
            rho: w.rho + side * 0.5 * s[0],
            u: w.u + side * 0.5 * s[1],
            p: w.p + side * 0.5 * s[2],
        };
        // Fall back to first order where reconstruction leaves the
        // physical region.
        if candidate.rho > 0.0 && candidate.p > 0.0 {
            candidate
        } else {
            w
        }
    };
    // Fluxes through the n + 1 interfaces between cells i+1 and i+2.
    let fluxes: Vec<Conserved> = (0..=n)
        .map(|k| hllc_flux(face(k + 1, 1.0), face(k + 2, -1.0), gamma))
        .collect();
    (0..n)
        .map(|i| [0, 1, 2].map(|c| -(fluxes[i + 1][c] - fluxes[i][c]) / dx))
        .collect()
}

/// Solves the 1D Euler equations (HLLC, MUSCL, SSP-RK3).
///
/// The domain is `[x0, x1]` with `cells` cells, the initial data `initial(x)`,
/// the final time `t_end`, the CFL number `cfl` and the reconstruction
/// `limiter`. Returns the cell centres and the primitive states.
///
/// This entry point delegates to the dimension-generic engine
/// [`crate::sim::euler_nd`] (1D case); grids with fewer than six cells or a
/// CFL number outside `(0, 1)` use the compact reference implementation.
#[must_use]
#[allow(clippy::too_many_arguments, clippy::option_if_let_else)]
pub fn solve_euler_1d<I>(
    initial: I,
    x0: f64,
    x1: f64,
    cells: usize,
    t_end: f64,
    gamma: f64,
    cfl: f64,
    limiter: Limiter,
    boundary: Boundary,
) -> (Vec<f64>, Vec<Primitive>)
where
    I: Fn(f64) -> Primitive,
{
    use crate::sim::euler_nd::FvConfig;
    use crate::sim::euler_nd::FvGrid;
    use crate::sim::euler_nd::GasBc;
    use crate::sim::euler_nd::Reconstruction;
    use crate::sim::euler_nd::SlopeLimiter;
    use crate::sim::euler_nd::State;
    use crate::sim::euler_nd::solve;
    if cells < 6 || !(cfl > 0.0 && cfl < 1.0) {
        return solve_euler_1d_reference(
            initial, x0, x1, cells, t_end, gamma, cfl, limiter, boundary,
        );
    }
    let Ok(grid) = FvGrid::new(&[cells], &[x0], &[x1]) else {
        return solve_euler_1d_reference(
            initial, x0, x1, cells, t_end, gamma, cfl, limiter, boundary,
        );
    };
    let mut cfg = FvConfig::euler(gamma);
    cfg.cfl = cfl;
    cfg.recon = match limiter {
        | Limiter::None => Reconstruction::FirstOrder,
        | Limiter::Minmod => Reconstruction::Muscl(SlopeLimiter::Minmod),
        | Limiter::VanLeer => Reconstruction::Muscl(SlopeLimiter::VanLeer),
    };
    cfg = cfg.with_axis(
        0,
        match boundary {
            | Boundary::Transmissive => GasBc::Transmissive,
            | Boundary::Reflective => GasBc::Reflective,
        },
    );
    let init = |x: [f64; 3]| -> State {
        let w = initial(x[0]);
        [w.rho, w.u, 0.0, 0.0, w.p]
    };
    match solve(&grid, &cfg, &init, t_end) {
        | Ok(sol) => {
            let centres: Vec<f64> = (0..cells).map(|i| grid.centre(i)[0]).collect();
            let prim = (0..cells)
                .map(|i| {
                    let w = sol.primitive(i);
                    Primitive {
                        rho: w[0],
                        u: w[1],
                        p: w[4],
                    }
                })
                .collect();
            (centres, prim)
        },
        | Err(_) => {
            solve_euler_1d_reference(initial, x0, x1, cells, t_end, gamma, cfl, limiter, boundary)
        },
    }
}

/// Compact reference implementation (small grids and unusual CFL numbers).
/// Solves the Euler equations on `[x0, x1]` with `cells` cells from the
/// initial data `initial(x)` up to time `t_end`:
///
/// SSP-RK3, HLLC fluxes, MUSCL reconstruction with `limiter`, CFL number
/// `cfl`. Returns the cell centres and the primitive states.
#[must_use]
#[allow(clippy::too_many_arguments)]
fn solve_euler_1d_reference<I>(
    initial: I,
    x0: f64,
    x1: f64,
    cells: usize,
    t_end: f64,
    gamma: f64,
    cfl: f64,
    limiter: Limiter,
    boundary: Boundary,
) -> (Vec<f64>, Vec<Primitive>)
where
    I: Fn(f64) -> Primitive,
{
    let dx = (x1 - x0) / cells as f64;
    let centres: Vec<f64> = (0..cells).map(|i| x0 + (i as f64 + 0.5) * dx).collect();
    let mut q: Vec<Conserved> = centres
        .iter()
        .map(|&x| to_conserved(initial(x), gamma))
        .collect();
    let mut t = 0.0;
    loop {
        crate::sim::models_geodesic_relativity::native_charge();
        if t >= t_end - 1e-14 {
            break;
        }
        let max_speed = q
            .iter()
            .map(|&c| {
                let w = to_primitive(c, gamma);
                w.u.abs() + sound_speed(w, gamma)
            })
            .fold(0.0, f64::max);
        let dt = (cfl * dx / max_speed).min(t_end - t);
        let add = |a: &[Conserved], b: &[Conserved], wa: f64, wb: f64, wc: f64| -> Vec<Conserved> {
            a.iter()
                .zip(b)
                .map(|(x, y)| [0, 1, 2].map(|c| wa * x[c] + wb * y[c] * wc))
                .collect()
        };
        let l0 = residual(&q, dx, gamma, limiter, boundary);
        let q1 = add(&q, &l0, 1.0, 1.0, dt);
        let l1 = residual(&q1, dx, gamma, limiter, boundary);
        let q2: Vec<Conserved> = add(&q, &add(&q1, &l1, 1.0, 1.0, dt), 0.75, 0.25, 1.0);
        let l2 = residual(&q2, dx, gamma, limiter, boundary);
        q = add(&q, &add(&q2, &l2, 1.0, 1.0, dt), 1.0 / 3.0, 2.0 / 3.0, 1.0);
        t += dt;
    }
    (
        centres,
        q.into_iter().map(|c| to_primitive(c, gamma)).collect(),
    )
}

/// Pressure function of one side and its derivative (Toro, ch. 4).
fn side_function(
    p: f64,
    w: Primitive,
    gamma: f64,
    c: f64,
) -> (f64, f64) {
    if p > w.p {
        // Shock.
        let a = 2.0 / ((gamma + 1.0) * w.rho);
        let b = (gamma - 1.0) / (gamma + 1.0) * w.p;
        let root = (a / (p + b)).sqrt();
        ((p - w.p) * root, root * (1.0 - 0.5 * (p - w.p) / (p + b)))
    } else {
        // Rarefaction.
        let ratio = p / w.p;
        let exponent = (gamma - 1.0) / (2.0 * gamma);
        (
            2.0 * c / (gamma - 1.0) * (ratio.powf(exponent) - 1.0),
            ratio.powf(-(gamma + 1.0) / (2.0 * gamma)) / (w.rho * c),
        )
    }
}

// A common rescaling of density and pressure leaves all velocities unchanged.
fn riemann_reference(
    left: Primitive,
    right: Primitive,
) -> (Primitive, Primitive, f64) {
    let reference = left.p.max(right.p);
    // At low scales the legacy Newton pressure floor dominates. At high
    // scales the product of the shock density and pressure factors can exceed
    // the range of f64 even though their square root is representable.
    if !(1e-4..=1e150).contains(&reference) && reference > 0.0 && reference.is_finite() {
        let scale = |state: Primitive| {
            Primitive {
                rho: state.rho / reference,
                p: state.p / reference,
                ..state
            }
        };
        let l = scale(left);
        let r = scale(right);
        if [l.rho, l.p, r.rho, r.p]
            .into_iter()
            .all(|value| value > 0.0 && value.is_finite())
        {
            return (l, r, reference);
        }
    }
    (left, right, 1.0)
}

fn star_pressure(
    l: Primitive,
    r: Primitive,
    gamma: f64,
    cl: f64,
    cr: f64,
) -> f64 {
    let z = (gamma - 1.0) / (2.0 * gamma);
    let estimate = ((cl + cr - 0.5 * (gamma - 1.0) * (r.u - l.u))
        / (cl / l.p.powf(z) + cr / r.p.powf(z)))
    .powf(1.0 / z);

    // This estimate is exact when both waves are rarefactions. Keep the
    // positive pressure arbitrarily close to vacuum, below the Newton floor.
    if estimate < 1e-10 && estimate <= l.p.min(r.p) {
        return estimate;
    }

    let mut p = estimate.max(1e-10);
    for _ in 0..100 {
        let (fl, dl) = side_function(p, l, gamma, cl);
        let (fr, dr) = side_function(p, r, gamma, cr);
        let next = (p - (fl + fr + r.u - l.u) / (dl + dr)).max(1e-12);
        let done = (next - p).abs() < 1e-14 * (next + p);
        p = next;
        if done {
            break;
        }
    }
    p
}

#[derive(Clone, Copy, Debug)]
struct RiemannSide {
    initial: Primitive,
    star: Primitive,
    sound: f64,
    head: f64,
    tail: f64,
}

impl RiemannSide {
    fn with_star<const LEFT: bool>(
        initial: Primitive,
        pressure: f64,
        sound: f64,
        gamma: f64,
        p_star: f64,
        star: Primitive,
    ) -> Self {
        let ratio = p_star / pressure;
        let (head, tail, rho) = if p_star > pressure {
            let k = (gamma - 1.0) / (gamma + 1.0);
            let jump = sound
                * ((gamma + 1.0) / (2.0 * gamma) * p_star / pressure
                    + (gamma - 1.0) / (2.0 * gamma))
                    .sqrt();
            let speed = if LEFT {
                initial.u - jump
            } else {
                initial.u + jump
            };
            // A shock occupies one similarity velocity and has no fan.
            (speed, speed, initial.rho * (ratio + k) / (k * ratio + 1.0))
        } else {
            let (c_star, rho) = if ratio.to_bits() == 1.0f64.to_bits() {
                (sound, initial.rho)
            } else {
                (
                    sound * ratio.powf((gamma - 1.0) / (2.0 * gamma)),
                    initial.rho * ratio.powf(1.0 / gamma),
                )
            };
            (
                if LEFT {
                    initial.u - sound
                } else {
                    initial.u + sound
                },
                if LEFT {
                    star.u - c_star
                } else {
                    star.u + c_star
                },
                rho,
            )
        };
        Self {
            initial,
            sound,
            head,
            tail,
            star: Primitive { rho, ..star },
        }
    }

    fn toward_vacuum<const LEFT: bool>(
        initial: Primitive,
        sound: f64,
        tail: f64,
    ) -> Self {
        Self {
            initial,
            sound,
            head: if LEFT {
                initial.u - sound
            } else {
                initial.u + sound
            },
            tail,
            star: vacuum(tail),
        }
    }

    fn sample<const LEFT: bool>(
        &self,
        gamma: f64,
        xi: f64,
    ) -> Primitive {
        if (LEFT && xi <= self.head) || (!LEFT && xi >= self.head) {
            return self.initial;
        }
        if (LEFT && xi >= self.tail) || (!LEFT && xi <= self.tail) {
            return self.star;
        }

        let c = self.sound;
        let u0 = self.initial.u;
        let (u, cf) = if LEFT {
            (
                2.0 / (gamma + 1.0) * (c + 0.5 * (gamma - 1.0) * u0 + xi),
                2.0 / (gamma + 1.0) * (c + 0.5 * (gamma - 1.0) * (u0 - xi)),
            )
        } else {
            (
                2.0 / (gamma + 1.0) * (-c + 0.5 * (gamma - 1.0) * u0 + xi),
                2.0 / (gamma + 1.0) * (c - 0.5 * (gamma - 1.0) * (u0 - xi)),
            )
        };
        // Rounding immediately inside a vacuum tail can make cf negative.
        let ratio = (cf / c).max(0.0);
        Primitive {
            rho: self.initial.rho * ratio.powf(2.0 / (gamma - 1.0)),
            u,
            p: self.initial.p * ratio.powf(2.0 * gamma / (gamma - 1.0)),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum RiemannInterface {
    Contact(f64),
    Vacuum { left: f64, right: f64 },
}

enum WaveLocation {
    Left,
    Vacuum,
    Right,
}

impl RiemannInterface {
    fn locate(
        self,
        xi: f64,
    ) -> WaveLocation {
        match self {
            | Self::Contact(speed) => {
                if xi <= speed {
                    WaveLocation::Left
                } else {
                    WaveLocation::Right
                }
            },
            | Self::Vacuum { left, right } => {
                if xi < left {
                    WaveLocation::Left
                } else if xi > right {
                    WaveLocation::Right
                } else {
                    WaveLocation::Vacuum
                }
            },
        }
    }
}

struct StarSolution {
    initial: [Primitive; 2],
    pressure: [f64; 2],
    sound: [f64; 2],
    star_pressure: f64,
    output_pressure: f64,
    interface: RiemannInterface,
    gamma: f64,
}

impl StarSolution {
    #[inline]
    fn solve(
        left: Primitive,
        right: Primitive,
        gamma: f64,
    ) -> Self {
        let (l, r, scale) = riemann_reference(left, right);
        let (cl, cr) = (sound_speed(l, gamma), sound_speed(r, gamma));
        if r.u > l.u {
            let left_vacuum = l.u + 2.0 * cl / (gamma - 1.0);
            let right_vacuum = r.u - 2.0 * cr / (gamma - 1.0);
            if left_vacuum <= right_vacuum {
                return Self {
                    initial: [left, right],
                    pressure: [l.p, r.p],
                    sound: [cl, cr],
                    star_pressure: 0.0,
                    output_pressure: 0.0,
                    interface: RiemannInterface::Vacuum {
                        left: left_vacuum,
                        right: right_vacuum,
                    },
                    gamma,
                };
            }
        }

        let star_pressure = star_pressure(l, r, gamma, cl, cr);
        let (fl, _) = side_function(star_pressure, l, gamma, cl);
        let (fr, _) = side_function(star_pressure, r, gamma, cr);
        let u_star = f64::midpoint(l.u, r.u) + 0.5 * (fr - fl);
        Self {
            initial: [left, right],
            pressure: [l.p, r.p],
            sound: [cl, cr],
            star_pressure,
            output_pressure: star_pressure * scale,
            interface: RiemannInterface::Contact(u_star),
            gamma,
        }
    }

    fn side<const LEFT: bool>(&self) -> RiemannSide {
        let i = usize::from(!LEFT);
        match self.interface {
            | RiemannInterface::Contact(speed) => {
                RiemannSide::with_star::<LEFT>(
                    self.initial[i],
                    self.pressure[i],
                    self.sound[i],
                    self.gamma,
                    self.star_pressure,
                    Primitive {
                        rho: 0.0,
                        u: speed,
                        p: self.output_pressure,
                    },
                )
            },
            | RiemannInterface::Vacuum { left, right } => {
                RiemannSide::toward_vacuum::<LEFT>(
                    self.initial[i],
                    self.sound[i],
                    if LEFT { left } else { right },
                )
            },
        }
    }
}

const fn vacuum(xi: f64) -> Primitive {
    Primitive {
        rho: 0.0,
        u: xi,
        p: 0.0,
    }
}

/// Prepared exact Riemann solution for two one-dimensional ideal-gas states.
///
/// Preparation determines the waves and star states once. The input states
/// require finite velocities, positive finite densities and pressures, and a
/// finite adiabatic index greater than one. A common change of density and
/// pressure scale preserves the wave velocities.
#[derive(Clone, Copy, Debug)]
pub struct ExactRiemann {
    left: RiemannSide,
    right: RiemannSide,
    interface: RiemannInterface,
    gamma: f64,
}

impl ExactRiemann {
    /// Solves the wave pattern for initial states separated by a discontinuity.
    #[must_use]
    pub fn new(
        left: Primitive,
        right: Primitive,
        gamma: f64,
    ) -> Self {
        let star = StarSolution::solve(left, right, gamma);
        Self {
            left: star.side::<true>(),
            right: star.side::<false>(),
            interface: star.interface,
            gamma,
        }
    }

    /// Samples at the self-similar velocity `xi = x/t`.
    ///
    /// At a contact this returns the left star state, and at a shock it
    /// returns the undisturbed state. In a vacuum, including its fan tails,
    /// density and pressure are zero and velocity is continued as `xi`.
    #[must_use]
    pub fn sample(
        &self,
        xi: f64,
    ) -> Primitive {
        match self.interface.locate(xi) {
            | WaveLocation::Left => self.left.sample::<true>(self.gamma, xi),
            | WaveLocation::Vacuum => vacuum(xi),
            | WaveLocation::Right => self.right.sample::<false>(self.gamma, xi),
        }
    }

    /// Samples a series of self-similar velocities in their input order.
    /// Large queries share the active CPU pool without altering any state value.
    #[must_use]
    pub fn sample_batch(
        &self,
        similarities: &[f64],
    ) -> Vec<Primitive> {
        crate::kernels::par::map_slice(similarities, 8192, |&xi| self.sample(xi))
    }
}

/// The exact solution of the Riemann problem `l | r` at `ξ = x/t`.
///
/// Use [`ExactRiemann`] to sample the same initial states at several points.
#[must_use]
pub fn exact_riemann(
    l: Primitive,
    r: Primitive,
    gamma: f64,
    xi: f64,
) -> Primitive {
    let star = StarSolution::solve(l, r, gamma);
    match star.interface.locate(xi) {
        | WaveLocation::Left => star.side::<true>().sample::<true>(gamma, xi),
        | WaveLocation::Vacuum => vacuum(xi),
        | WaveLocation::Right => star.side::<false>().sample::<false>(gamma, xi),
    }
}

/// Sod's shock tube on `[0, 1]` at `t`: numerical and exact densities at
/// the cell centres.
#[must_use]
pub fn sod_shock_tube(
    cells: usize,
    t: f64,
    limiter: Limiter,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let gamma = 1.4;
    let left = Primitive {
        rho: 1.0,
        u: 0.0,
        p: 1.0,
    };
    let right = Primitive {
        rho: 0.125,
        u: 0.0,
        p: 0.1,
    };
    let initial = |x| if x < 0.5 { left } else { right };
    let (x, w) = solve_euler_1d(
        initial,
        0.0,
        1.0,
        cells,
        t,
        gamma,
        0.5,
        limiter,
        Boundary::Transmissive,
    );
    let exact = if t == 0.0 {
        x.iter().map(|&xc| initial(xc).rho).collect()
    } else {
        let solution = ExactRiemann::new(left, right, gamma);
        x.iter()
            .map(|&xc| solution.sample((xc - 0.5) / t).rho)
            .collect()
    };
    (x, w.into_iter().map(|s| s.rho).collect(), exact)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_riemann_star_state() {
        // Toro's test 1 (Sod): p* = 0.30313, u* = 0.92745.
        let l = Primitive {
            rho: 1.0,
            u: 0.0,
            p: 1.0,
        };
        let r = Primitive {
            rho: 0.125,
            u: 0.0,
            p: 0.1,
        };
        let star = exact_riemann(l, r, 1.4, 0.5);
        assert!(
            (star.p - 0.303_13).abs() < 1e-4 && (star.u - 0.927_45).abs() < 1e-4,
            "{star:?}"
        );
        // Left star state (between rarefaction tail and contact) and right
        // star state (between contact and shock).
        assert!((star.rho - 0.426_32).abs() < 1e-4, "{star:?}");
        assert!((exact_riemann(l, r, 1.4, 1.2).rho - 0.265_57).abs() < 1e-4);
    }

    #[test]
    fn sod_tube_converges_to_the_exact_solution() {
        let error = |cells: usize, limiter: Limiter| {
            let (_, numeric, exact) = sod_shock_tube(cells, 0.2, limiter);
            numeric
                .iter()
                .zip(&exact)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                / cells as f64
        };
        let (coarse, fine) = (error(100, Limiter::VanLeer), error(400, Limiter::VanLeer));
        assert!(fine < 0.006 && fine < coarse * 0.6, "{coarse} {fine}");
        assert!(
            error(200, Limiter::Minmod) < error(200, Limiter::None),
            "MUSCL beats first order"
        );
    }

    #[test]
    fn conservation_with_walls() {
        // Total mass and energy are conserved in a closed tube.
        let gamma = 1.4;
        let init = |x: f64| {
            Primitive {
                rho: 1.0 + 0.5 * (-(x - 0.5).powi(2) / 0.01).exp(),
                u: 0.0,
                p: 1.0,
            }
        };
        let (_, w0) = solve_euler_1d(
            init,
            0.0,
            1.0,
            200,
            0.0,
            gamma,
            0.5,
            Limiter::Minmod,
            Boundary::Reflective,
        );
        let (_, w1) = solve_euler_1d(
            init,
            0.0,
            1.0,
            200,
            0.3,
            gamma,
            0.5,
            Limiter::Minmod,
            Boundary::Reflective,
        );
        let totals = |w: &[Primitive]| {
            let q: Vec<[f64; 3]> = w.iter().map(|&s| to_conserved(s, gamma)).collect();
            (
                q.iter().map(|c| c[0]).sum::<f64>(),
                q.iter().map(|c| c[2]).sum::<f64>(),
            )
        };
        let ((m0, e0), (m1, e1)) = (totals(&w0), totals(&w1));
        assert!(
            (m0 - m1).abs() < 1e-10 * m0 && (e0 - e1).abs() < 1e-10 * e0,
            "{m0} {m1} {e0} {e1}"
        );
    }
}
