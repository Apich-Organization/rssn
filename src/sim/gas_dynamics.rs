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
    [w.rho, w.rho * w.u, w.p / (gamma - 1.0) + 0.5 * w.rho * w.u * w.u]
}

fn to_primitive(
    q: Conserved,
    gamma: f64,
) -> Primitive {
    let rho = q[0];
    let u = q[1] / rho;
    Primitive { rho, u, p: (gamma - 1.0) * (q[2] - 0.5 * rho * u * u) }
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
    let s_star = (r.p - l.p + l.rho * l.u * (sl - l.u) - r.rho * r.u * (sr - r.u)) / (l.rho * (sl - l.u) - r.rho * (sr - r.u));
    let star = |w: Primitive, q: Conserved, s: f64| -> Conserved {
        let factor = w.rho * (s - w.u) / (s - s_star);
        let energy = q[2] / w.rho + (s_star - w.u) * (s_star + w.p / (w.rho * (s - w.u)));
        [factor, factor * s_star, factor * energy]
    };
    if s_star >= 0.0 {
        let qs = star(l, ql, sl);
        [fl[0] + sl * (qs[0] - ql[0]), fl[1] + sl * (qs[1] - ql[1]), fl[2] + sl * (qs[2] - ql[2])]
    } else {
        let qs = star(r, qr, sr);
        [fr[0] + sr * (qs[0] - qr[0]), fr[1] + sr * (qs[1] - qr[1]), fr[2] + sr * (qs[2] - qr[2])]
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
    let ghost = |w: Primitive| match boundary {
        | Boundary::Transmissive => w,
        | Boundary::Reflective => Primitive { u: -w.u, ..w },
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
        let candidate = Primitive { rho: w.rho + side * 0.5 * s[0], u: w.u + side * 0.5 * s[1], p: w.p + side * 0.5 * s[2] };
        // Fall back to first order where reconstruction leaves the
        // physical region.
        if candidate.rho > 0.0 && candidate.p > 0.0 { candidate } else { w }
    };
    // Fluxes through the n + 1 interfaces between cells i+1 and i+2.
    let fluxes: Vec<Conserved> = (0..=n).map(|k| hllc_flux(face(k + 1, 1.0), face(k + 2, -1.0), gamma)).collect();
    (0..n).map(|i| [0, 1, 2].map(|c| -(fluxes[i + 1][c] - fluxes[i][c]) / dx)).collect()
}

/// Solves the Euler equations on `[x0, x1]` with `cells` cells from the
/// initial data `initial(x)` up to time `t_end`:
///
/// SSP-RK3, HLLC fluxes, MUSCL reconstruction with `limiter`, CFL number
/// `cfl`. Returns the cell centres and the primitive states.
#[must_use]
#[allow(clippy::too_many_arguments)]
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
    let dx = (x1 - x0) / cells as f64;
    let centres: Vec<f64> = (0..cells).map(|i| x0 + (i as f64 + 0.5) * dx).collect();
    let mut q: Vec<Conserved> = centres.iter().map(|&x| to_conserved(initial(x), gamma)).collect();
    let mut t = 0.0;
    loop {
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
            a.iter().zip(b).map(|(x, y)| [0, 1, 2].map(|c| wa * x[c] + wb * y[c] * wc)).collect()
        };
        let l0 = residual(&q, dx, gamma, limiter, boundary);
        let q1 = add(&q, &l0, 1.0, 1.0, dt);
        let l1 = residual(&q1, dx, gamma, limiter, boundary);
        let q2: Vec<Conserved> = add(&q, &add(&q1, &l1, 1.0, 1.0, dt), 0.75, 0.25, 1.0);
        let l2 = residual(&q2, dx, gamma, limiter, boundary);
        q = add(&q, &add(&q2, &l2, 1.0, 1.0, dt), 1.0 / 3.0, 2.0 / 3.0, 1.0);
        t += dt;
    }
    (centres, q.into_iter().map(|c| to_primitive(c, gamma)).collect())
}

/// Pressure function of one side and its derivative (Toro, ch. 4).
fn side_function(
    p: f64,
    w: Primitive,
    gamma: f64,
) -> (f64, f64) {
    let c = sound_speed(w, gamma);
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
        (2.0 * c / (gamma - 1.0) * (ratio.powf(exponent) - 1.0), ratio.powf(-(gamma + 1.0) / (2.0 * gamma)) / (w.rho * c))
    }
}

/// The exact solution of the Riemann problem `l | r` at `ξ = x/t`.
#[must_use]
pub fn exact_riemann(
    l: Primitive,
    r: Primitive,
    gamma: f64,
    xi: f64,
) -> Primitive {
    // Star pressure by Newton from the two-rarefaction guess.
    let (cl, cr) = (sound_speed(l, gamma), sound_speed(r, gamma));
    let z = (gamma - 1.0) / (2.0 * gamma);
    let mut p = ((cl + cr - 0.5 * (gamma - 1.0) * (r.u - l.u)) / (cl / l.p.powf(z) + cr / r.p.powf(z))).powf(1.0 / z).max(1e-10);
    for _ in 0..100 {
        let (fl, dl) = side_function(p, l, gamma);
        let (fr, dr) = side_function(p, r, gamma);
        let next = (p - (fl + fr + r.u - l.u) / (dl + dr)).max(1e-12);
        let done = (next - p).abs() < 1e-14 * (next + p);
        p = next;
        if done {
            break;
        }
    }
    let (fl, _) = side_function(p, l, gamma);
    let (fr, _) = side_function(p, r, gamma);
    let u_star = f64::midpoint(l.u, r.u) + 0.5 * (fr - fl);
    let g = gamma;
    let star_density = |w: Primitive| -> f64 {
        if p > w.p {
            let ratio = p / w.p;
            let k = (g - 1.0) / (g + 1.0);
            w.rho * (ratio + k) / (k * ratio + 1.0)
        } else {
            w.rho * (p / w.p).powf(1.0 / g)
        }
    };
    if xi <= u_star {
        // Left of the contact.
        let c = cl;
        if p > l.p {
            let s = l.u - c * ((g + 1.0) / (2.0 * g) * p / l.p + (g - 1.0) / (2.0 * g)).sqrt();
            if xi <= s { l } else { Primitive { rho: star_density(l), u: u_star, p } }
        } else {
            let head = l.u - c;
            let c_star = c * (p / l.p).powf((g - 1.0) / (2.0 * g));
            let tail = u_star - c_star;
            if xi <= head {
                l
            } else if xi >= tail {
                Primitive { rho: star_density(l), u: u_star, p }
            } else {
                let u = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * l.u + xi);
                let cf = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * (l.u - xi));
                let rho = l.rho * (cf / c).powf(2.0 / (g - 1.0));
                Primitive { rho, u, p: l.p * (cf / c).powf(2.0 * g / (g - 1.0)) }
            }
        }
    } else {
        let c = cr;
        if p > r.p {
            let s = r.u + c * ((g + 1.0) / (2.0 * g) * p / r.p + (g - 1.0) / (2.0 * g)).sqrt();
            if xi >= s { r } else { Primitive { rho: star_density(r), u: u_star, p } }
        } else {
            let head = r.u + c;
            let c_star = c * (p / r.p).powf((g - 1.0) / (2.0 * g));
            let tail = u_star + c_star;
            if xi >= head {
                r
            } else if xi <= tail {
                Primitive { rho: star_density(r), u: u_star, p }
            } else {
                let u = 2.0 / (g + 1.0) * (-c + 0.5 * (g - 1.0) * r.u + xi);
                let cf = 2.0 / (g + 1.0) * (c - 0.5 * (g - 1.0) * (r.u - xi));
                let rho = r.rho * (cf / c).powf(2.0 / (g - 1.0));
                Primitive { rho, u, p: r.p * (cf / c).powf(2.0 * g / (g - 1.0)) }
            }
        }
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
    let left = Primitive { rho: 1.0, u: 0.0, p: 1.0 };
    let right = Primitive { rho: 0.125, u: 0.0, p: 0.1 };
    let (x, w) = solve_euler_1d(|x| if x < 0.5 { left } else { right }, 0.0, 1.0, cells, t, gamma, 0.5, limiter, Boundary::Transmissive);
    let exact: Vec<f64> = x.iter().map(|&xc| exact_riemann(left, right, gamma, (xc - 0.5) / t).rho).collect();
    (x, w.into_iter().map(|s| s.rho).collect(), exact)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_riemann_star_state() {
        // Toro's test 1 (Sod): p* = 0.30313, u* = 0.92745.
        let l = Primitive { rho: 1.0, u: 0.0, p: 1.0 };
        let r = Primitive { rho: 0.125, u: 0.0, p: 0.1 };
        let star = exact_riemann(l, r, 1.4, 0.5);
        assert!((star.p - 0.303_13).abs() < 1e-4 && (star.u - 0.927_45).abs() < 1e-4, "{star:?}");
        // Left star state (between rarefaction tail and contact) and right
        // star state (between contact and shock).
        assert!((star.rho - 0.426_32).abs() < 1e-4, "{star:?}");
        assert!((exact_riemann(l, r, 1.4, 1.2).rho - 0.265_57).abs() < 1e-4);
    }

    #[test]
    fn sod_tube_converges_to_the_exact_solution() {
        let error = |cells: usize, limiter: Limiter| {
            let (_, numeric, exact) = sod_shock_tube(cells, 0.2, limiter);
            numeric.iter().zip(&exact).map(|(a, b)| (a - b).abs()).sum::<f64>() / cells as f64
        };
        let (coarse, fine) = (error(100, Limiter::VanLeer), error(400, Limiter::VanLeer));
        assert!(fine < 0.006 && fine < coarse * 0.6, "{coarse} {fine}");
        assert!(error(200, Limiter::Minmod) < error(200, Limiter::None), "MUSCL beats first order");
    }

    #[test]
    fn conservation_with_walls() {
        // Total mass and energy are conserved in a closed tube.
        let gamma = 1.4;
        let init = |x: f64| Primitive { rho: 1.0 + 0.5 * (-(x - 0.5).powi(2) / 0.01).exp(), u: 0.0, p: 1.0 };
        let (_, w0) = solve_euler_1d(init, 0.0, 1.0, 200, 0.0, gamma, 0.5, Limiter::Minmod, Boundary::Reflective);
        let (_, w1) = solve_euler_1d(init, 0.0, 1.0, 200, 0.3, gamma, 0.5, Limiter::Minmod, Boundary::Reflective);
        let totals = |w: &[Primitive]| {
            let q: Vec<[f64; 3]> = w.iter().map(|&s| to_conserved(s, gamma)).collect();
            (q.iter().map(|c| c[0]).sum::<f64>(), q.iter().map(|c| c[2]).sum::<f64>())
        };
        let ((m0, e0), (m1, e1)) = (totals(&w0), totals(&w1));
        assert!((m0 - m1).abs() < 1e-10 * m0 && (e0 - e1).abs() < 1e-10 * e0, "{m0} {m1} {e0} {e1}");
    }
}
