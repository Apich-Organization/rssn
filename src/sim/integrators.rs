//! Symplectic and stiff (implicit, Rosenbrock) integrators.
//!
//! # Structure-preserving and stiff integrators
//!
//! * **Symplectic** integrators for separable Hamiltonians
//!   `H(q, p) = T(p) + V(q)` (here `T = |p|²/2m`): Störmer–Verlet (order
//!   2) and the Forest–Ruth / Yoshida composition (order 4). They conserve
//!   a modified energy, so the energy error stays bounded over very long
//!   runs instead of drifting.
//! * **Stiff** integrators for `y' = f(t, y)`: backward Euler and the
//!   implicit midpoint rule (Newton iterations with a finite-difference
//!   Jacobian), and an adaptive L-stable Rosenbrock method (ROS2, with the
//!   linearly implicit Euler step as error estimate), which needs one
//!   Jacobian and two linear solves per step and no Newton iteration.
//!
//! Linear systems are solved by Gaussian elimination with partial pivoting
//! (the systems are small: one row per state variable).
//!
//! ## Coverage
//!
//! | Method | Order | Dimension | Notes |
//! |---|---|---|---|
//! | Verlet, Yoshida-4 (`verlet_step`, `yoshida4_step`) | 2, 4 | any (state slices) | separable `H` |
//! | Backward Euler, implicit midpoint | 1, 2 | any | Newton with FD Jacobian |
//! | Rosenbrock ROS2 | 2 | any | adaptive, L-stable |
//!
//! Higher-order (6, 8) and Blanes-Moan schemes, time-transformed and
//! regularised integrators and rigid-body splittings are in
//! [`crate::sim::symplectic`].
//! [`symplectic_step_with_workspace`] reuses a caller-provided acceleration buffer.

use crate::kernels::exact_eq;
use crate::kernels::numeric_common::scale_product;

/// Solves `a x = b` by Gaussian elimination with partial pivoting;
/// `None` for a (numerically) singular matrix.
/// The Jacobian `∂f/∂y` at `(t, y)` by central differences.
fn jacobian<F>(
    f: &F,
    t: f64,
    y: &[f64],
) -> Vec<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let mut jac = vec![vec![0.0; n]; n];
    let mut shifted = y.to_vec();
    for j in 0..n {
        let h = 1e-7 * (1.0 + y[j].abs());
        shifted[j] = y[j] + h;
        let plus = f(t, &shifted);
        shifted[j] = y[j] - h;
        let minus = f(t, &shifted);
        shifted[j] = y[j];
        for i in 0..n {
            jac[i][j] = (plus[i] - minus[i]) / (2.0 * h);
        }
    }
    jac
}

fn verlet_kick(
    p: &mut [f64],
    a: &[f64],
    mass: f64,
    dt: f64,
) {
    let half_dt = 0.5 * dt;
    let coefficient = half_dt * mass;
    if (half_dt.is_normal() && coefficient.is_normal()) || !valid_step_scales(mass, dt) {
        for (pi, ai) in p.iter_mut().zip(a) {
            *pi += coefficient * ai;
        }
    } else {
        let factors = [0.5, dt.abs(), mass];
        for (pi, &ai) in p.iter_mut().zip(a) {
            let kick = if ai.is_finite() {
                step_product(ai, dt, &factors, &[])
            } else {
                coefficient * ai
            };
            *pi += kick;
        }
    }
}

fn verlet_drift(
    q: &mut [f64],
    p: &[f64],
    mass: f64,
    dt: f64,
) {
    let (p, _) = p.split_at(q.len().min(p.len()));
    if exact_eq(mass, 1.0) {
        for (qi, pi) in q.iter_mut().zip(p) {
            *qi += dt * pi;
        }
    } else if valid_step_scales(mass, dt) && !all_normal_drift_products(p, dt) {
        let factors = [dt.abs()];
        let divisors = [mass];
        for (qi, &pi) in q.iter_mut().zip(p) {
            let product = dt * pi;
            let drift = if product.is_normal() || exact_eq(pi, 0.0) || !pi.is_finite() {
                product / mass
            } else {
                step_product(pi, dt, &factors, &divisors)
            };
            *qi += drift;
        }
    } else {
        for (qi, pi) in q.iter_mut().zip(p) {
            *qi += dt * pi / mass;
        }
    }
}

fn all_normal_drift_products(
    p: &[f64],
    dt: f64,
) -> bool {
    let exponent = dt.abs().to_bits() >> 52;
    if exponent == 0 {
        return false;
    }
    // With biased exponents E and F the exact product is in
    // [2^(E+F-2046), 2^(E+F-2044)). Bound both significands by that interval.
    let lower = 1024_u64.saturating_sub(exponent).max(1) << 52;
    let upper = (3068 - exponent).min(2047) << 52;
    let width = upper - lower;
    // Use a reduction without early exit so the integer range check can vectorize.
    p.iter().fold(true, |safe, p| {
        let magnitude = p.to_bits() & 0x7fff_ffff_ffff_ffff;
        safe & (magnitude.wrapping_sub(lower) < width || magnitude == 0)
    })
}

fn valid_step_scales(
    mass: f64,
    dt: f64,
) -> bool {
    mass.is_finite() && mass > 0.0 && dt.is_finite() && !exact_eq(dt, 0.0)
}

#[cold]
fn step_product(
    value: f64,
    dt: f64,
    factors: &[f64],
    divisors: &[f64],
) -> f64 {
    let signed = if dt.is_sign_negative() {
        -value
    } else {
        value
    };
    scale_product(signed, factors, divisors).unwrap_or_else(|_| f64::INFINITY.copysign(signed))
}

fn symplectic_step_by(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    intervals: &[f64],
    mut kick: impl FnMut(&[f64], &mut [f64], f64),
) {
    for &h in intervals {
        kick(q, p, h);
        verlet_drift(q, p, mass, h);
        kick(q, p, h);
    }
}

fn yoshida_intervals(dt: f64) -> [f64; 3] {
    let cbrt2 = 2.0_f64.cbrt();
    let w1 = 1.0 / (2.0 - cbrt2);
    let w0 = -cbrt2 * w1;
    [w1 * dt, w0 * dt, w1 * dt]
}

/// One step of the Störmer–Verlet (leapfrog) method for `q'' = a(q)`
/// (`a` the acceleration, `-∇V/m`): kick, drift, kick.
pub fn verlet_step<A>(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    dt: f64,
    acceleration: &A,
) where
    A: Fn(&[f64]) -> Vec<f64>,
{
    symplectic_step_by(q, p, mass, &[dt], |q, p, h| {
        verlet_kick(p, &acceleration(q), mass, h);
    });
}

/// One step of the fourth-order Forest–Ruth (Yoshida) composition of
/// Verlet steps with weights `w1, w0, w1`, `w1 = 1/(2 - 2^(1/3))`.
pub fn yoshida4_step<A>(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    dt: f64,
    acceleration: &A,
) where
    A: Fn(&[f64]) -> Vec<f64>,
{
    symplectic_step_by(q, p, mass, &yoshida_intervals(dt), |q, p, h| {
        verlet_kick(p, &acceleration(q), mass, h);
    });
}

/// Integrates a separable Hamiltonian system for `steps` steps with the
/// symplectic method of the given `order` (2: Verlet, 4: Yoshida),
/// returning the states `(q, p)` after every step.
#[must_use]
pub fn symplectic_trajectory<A>(
    q0: &[f64],
    p0: &[f64],
    mass: f64,
    dt: f64,
    steps: usize,
    order: u32,
    acceleration: A,
) -> Vec<(Vec<f64>, Vec<f64>)>
where
    A: Fn(&[f64]) -> Vec<f64>,
{
    let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
    let mut out = Vec::with_capacity(steps + 1);
    out.push((q.clone(), p.clone()));
    for _ in 0..steps {
        if order >= 4 {
            yoshida4_step(&mut q, &mut p, mass, dt, &acceleration);
        } else {
            verlet_step(&mut q, &mut p, mass, dt, &acceleration);
        }
        out.push((q.clone(), p.clone()));
    }
    out
}

/// Symplectic integrator scheme for separable Hamiltonians `H(q, p) = T(p) + V(q)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymplecticIntegratorMethod {
    /// Second-order Störmer–Verlet (leapfrog).
    Verlet,
    /// Fourth-order Forest–Ruth / Yoshida composition.
    Yoshida4,
}

/// Advances a separable Hamiltonian by one signed step using an acceleration buffer.
///
/// `q` is position, `p` is momentum, and `mass` is a shared positive scalar mass.
/// The callback must fill every entry of `work` with the acceleration at its position argument.
/// The integrator makes two callback calls for Verlet and six for Yoshida4, in kick order.
/// The workspace is reused by all of these calls and may be retained between steps.
///
/// # Errors
/// Invalid dimensions, mass or unrepresentable step intervals are rejected before
/// touching `q`, `p` or `work`.
/// Once started, the step mutates `q` and `p` in place and propagates the callback's values.
pub fn symplectic_step_with_workspace(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    dt: f64,
    method: SymplecticIntegratorMethod,
    work: &mut [f64],
    mut acceleration: impl FnMut(&[f64], &mut [f64]),
) -> Result<(), String> {
    if q.is_empty() || q.len() != p.len() || q.len() != work.len() {
        return Err(
            "positions, momenta and acceleration workspace must have the same nonzero length"
                .into(),
        );
    }
    if !mass.is_finite() || mass <= 0.0 || !dt.is_finite() {
        return Err("symplectic step requires positive finite mass and a finite time step".into());
    }
    let single = [dt];
    let yoshida;
    let intervals: &[f64] = match method {
        | SymplecticIntegratorMethod::Verlet => &single,
        | SymplecticIntegratorMethod::Yoshida4 => {
            yoshida = yoshida_intervals(dt);
            if yoshida.iter().any(|h| !h.is_finite()) {
                return Err("Yoshida substep is not representable".into());
            }
            &yoshida
        },
    };
    symplectic_step_by(q, p, mass, intervals, |q, p, h| {
        acceleration(q, work);
        verlet_kick(p, work, mass, h);
    });
    Ok(())
}

/// Trajectory returned by [`solve_symplectic`].
#[derive(Debug, Clone, PartialEq)]
pub struct SymplecticTrajectory {
    /// Sample times.
    pub t: Vec<f64>,
    /// Generalised positions at each sample time.
    pub q: Vec<Vec<f64>>,
    /// Generalised momenta at each sample time.
    pub p: Vec<Vec<f64>>,
    /// Number of rejected adaptive steps (zero in fixed-step mode).
    pub rejected: usize,
}

/// Unified symplectic integrator supporting both fixed and adaptive step sizes
/// (adaptive via Richardson step-doubling error control).
///
/// # Errors
/// Returns an error when inputs are invalid, `mass <= 0`, states become
/// non-finite, or the adaptive step budget or step floor is exceeded.
#[allow(clippy::float_cmp)]
pub fn solve_symplectic<A>(
    q0: &[f64],
    p0: &[f64],
    mass: f64,
    t_span: (f64, f64),
    method: SymplecticIntegratorMethod,
    step_mode: crate::kernels::ode::OdeStepMode,
    acceleration: A,
) -> Result<SymplecticTrajectory, String>
where
    A: Fn(&[f64]) -> Vec<f64>,
{
    let (t0, t1) = t_span;
    if q0.is_empty()
        || q0.len() != p0.len()
        || !mass.is_finite()
        || mass <= 0.0
        || !t0.is_finite()
        || !t1.is_finite()
        || q0.iter().chain(p0).any(|v| !v.is_finite())
    {
        return Err("Invalid symplectic integration inputs.".to_string());
    }
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(SymplecticTrajectory {
            t: vec![t0],
            q: vec![q0.to_vec()],
            p: vec![p0.to_vec()],
            rejected: 0,
        });
    }
    let advance_one = |q: &mut Vec<f64>, p: &mut Vec<f64>, dt: f64| {
        match method {
            | SymplecticIntegratorMethod::Verlet => verlet_step(q, p, mass, dt, &acceleration),
            | SymplecticIntegratorMethod::Yoshida4 => yoshida4_step(q, p, mass, dt, &acceleration),
        }
    };
    match step_mode {
        | crate::kernels::ode::OdeStepMode::Fixed { dt } => {
            if !dt.is_finite() || dt <= 0.0 {
                return Err("Fixed step size dt must be positive and finite.".to_string());
            }
            let steps = ((t1 - t0).abs() / dt).ceil().max(1.0) as usize;
            let h = (t1 - t0) / steps as f64;
            let mut q = q0.to_vec();
            let mut p = p0.to_vec();
            let mut ts = Vec::with_capacity(steps.saturating_add(1));
            let mut qs = Vec::with_capacity(steps.saturating_add(1));
            let mut ps = Vec::with_capacity(steps.saturating_add(1));
            ts.push(t0);
            qs.push(q.clone());
            ps.push(p.clone());
            for step_idx in 1..=steps {
                crate::cancellation::checkpoint().map_err(|e| e.to_string())?;
                advance_one(&mut q, &mut p, h);
                if q.iter().chain(&p).any(|v| !v.is_finite()) {
                    return Err(
                        "Non-finite state encountered in symplectic integration.".to_string()
                    );
                }
                let tn = if step_idx == steps {
                    t1
                } else {
                    t0 + step_idx as f64 * h
                };
                ts.push(tn);
                qs.push(q.clone());
                ps.push(p.clone());
            }
            Ok(SymplecticTrajectory {
                t: ts,
                q: qs,
                p: ps,
                rejected: 0,
            })
        },
        | crate::kernels::ode::OdeStepMode::Adaptive { .. } => {
            let opts = step_mode.to_ode_options().map_err(|e| e.to_string())?;
            let order = match method {
                | SymplecticIntegratorMethod::Verlet => 2.0_f64,
                | SymplecticIntegratorMethod::Yoshida4 => 4.0_f64,
            };
            let denom = 2.0_f64.powf(order) - 1.0;
            let dir = if t1 > t0 { 1.0 } else { -1.0 };
            let mut h = opts
                .h0
                .unwrap_or_else(|| ((t1 - t0).abs() * 0.01).max(1e-6))
                .abs()
                .clamp(opts.hmin.max(1e-14), opts.hmax)
                * dir;
            let mut t = t0;
            let mut q = q0.to_vec();
            let mut p = p0.to_vec();
            let mut ts = vec![t0];
            let mut qs = vec![q.clone()];
            let mut ps = vec![p.clone()];
            let mut rejected = 0usize;
            let mut attempts = 0usize;
            while (t1 - t) * dir > 0.0 {
                crate::cancellation::checkpoint().map_err(|e| e.to_string())?;
                if attempts >= opts.max_steps {
                    return Err("Symplectic adaptive step limit exceeded.".to_string());
                }
                attempts = attempts.saturating_add(1);
                if (t + h - t1) * dir > 0.0 {
                    h = t1 - t;
                }
                if !h.is_finite()
                    || crate::kernels::exact_eq(t + h, t)
                    || (h.abs() < opts.hmin && (t1 - t).abs() > opts.hmin)
                {
                    return Err("Symplectic adaptive step size underflow.".to_string());
                }
                let (mut qf, mut pf) = (q.clone(), p.clone());
                advance_one(&mut qf, &mut pf, h);
                let (mut qh, mut ph) = (q.clone(), p.clone());
                advance_one(&mut qh, &mut ph, 0.5 * h);
                advance_one(&mut qh, &mut ph, 0.5 * h);
                if qh
                    .iter()
                    .chain(&ph)
                    .chain(&qf)
                    .chain(&pf)
                    .any(|v| !v.is_finite())
                {
                    h *= 0.5;
                    rejected = rejected.saturating_add(1);
                    continue;
                }
                let n_tot = q.len().saturating_add(p.len()) as f64;
                let sq_sum: f64 = q
                    .iter()
                    .zip(&qh)
                    .zip(&qf)
                    .chain(p.iter().zip(&ph).zip(&pf))
                    .map(|((&y0, &yh), &yf)| {
                        let err_i = (yh - yf).abs() / denom;
                        let scale = opts.atol + opts.rtol * y0.abs().max(yh.abs());
                        (err_i / scale).powi(2)
                    })
                    .sum();
                let en = (sq_sum / n_tot).sqrt();
                if en.is_nan() {
                    return Err("Non-finite error estimate in symplectic integration.".to_string());
                }
                if en <= 1.0 {
                    t += h;
                    q = qh;
                    p = ph;
                    ts.push(t);
                    qs.push(q.clone());
                    ps.push(p.clone());
                } else {
                    rejected = rejected.saturating_add(1);
                }
                let factor = if en == 0.0 {
                    5.0
                } else {
                    (0.9 * en.powf(-1.0 / (order + 1.0))).clamp(0.2, 5.0)
                };
                h = (h.abs() * factor).min(opts.hmax) * dir;
            }
            Ok(SymplecticTrajectory {
                t: ts,
                q: qs,
                p: ps,
                rejected,
            })
        },
    }
}

/// Solves `g(z) = 0` by Newton's method from `z0` with the
/// finite-difference Jacobian of `g`.
fn newton<G>(
    g: &G,
    z0: Vec<f64>,
) -> Option<Vec<f64>>
where
    G: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut z = z0;
    for _ in 0..50 {
        let r = g(0.0, &z);
        let norm: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-12 * (1.0 + z.iter().map(|v| v.abs()).fold(0.0, f64::max)) {
            return Some(z);
        }
        let jac = jacobian(g, 0.0, &z);
        let residual: Vec<_> = r.iter().map(|v| -v).collect();
        let delta = crate::kernels::dense::solve_linear_rows(&jac, &residual).ok()?;
        for (zi, di) in z.iter_mut().zip(&delta) {
            *zi += di;
        }
        if delta.iter().map(|v| v.abs()).fold(0.0, f64::max) < 1e-14 {
            return Some(z);
        }
    }
    Some(z)
}

/// One backward Euler step: `y1 = y0 + h f(t + h, y1)`.
pub fn backward_euler_step<F>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
) -> Option<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let residual = |_: f64, z: &[f64]| -> Vec<f64> {
        let fz = f(t + h, z);
        z.iter()
            .zip(y)
            .zip(&fz)
            .map(|((zi, yi), fi)| zi - yi - h * fi)
            .collect()
    };
    newton(&residual, y.to_vec())
}

/// One implicit midpoint step: `y1 = y0 + h f(t + h/2, (y0 + y1)/2)`
/// (symplectic and A-stable, order 2).
pub fn implicit_midpoint_step<F>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
) -> Option<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let residual = |_: f64, z: &[f64]| -> Vec<f64> {
        let mid: Vec<f64> = z.iter().zip(y).map(|(zi, yi)| 0.5 * (zi + yi)).collect();
        let fm = f(t + 0.5 * h, &mid);
        z.iter()
            .zip(y)
            .zip(&fm)
            .map(|((zi, yi), fi)| zi - yi - h * fi)
            .collect()
    };
    newton(&residual, y.to_vec())
}

/// Integrates `y' = f(t, y)` from `t0` to `t1` with the adaptive Rosenbrock
/// method ROS2 (L-stable, order 2; `γ = 1 + 1/√2`), the error estimated
/// against the linearly implicit Euler step.
///
/// Returns the accepted `(t, y)` pairs, or `None` when a linear system is
/// singular or the step size underflows.
pub fn rosenbrock_solve<F>(
    f: F,
    t0: f64,
    y0: &[f64],
    t1: f64,
    tolerance: f64,
) -> Option<Vec<(f64, Vec<f64>)>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let gamma = 1.0 + 1.0 / 2.0_f64.sqrt();
    let n = y0.len();
    let (mut t, mut y) = (t0, y0.to_vec());
    let mut h = ((t1 - t0) / 100.0).abs().max(1e-12) * (t1 - t0).signum();
    let mut out = vec![(t, y.clone())];
    let mut guard = 0_usize;
    loop {
        if (t1 - t) * h.signum() <= 1e-14 * (1.0 + t1.abs()) {
            break;
        }
        guard += 1;
        if guard > 1_000_000 || h.abs() < 1e-14 * (1.0 + t.abs()) {
            return None;
        }
        if (t + h - t1) * h.signum() > 0.0 {
            h = t1 - t;
        }
        let jac = jacobian(&f, t, &y);
        // W = I - γ h J
        let w: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| f64::from(u8::from(i == j)) - gamma * h * jac[i][j])
                    .collect()
            })
            .collect();
        let f0 = f(t, &y);
        if f0.len() != n || f0.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let k1 = crate::kernels::dense::solve_linear_rows(&w, &f0).ok()?;
        let y_mid: Vec<f64> = y.iter().zip(&k1).map(|(yi, ki)| yi + h * ki).collect();
        let f1 = f(t + h, &y_mid);
        if f1.len() != n || f1.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let rhs: Vec<f64> = f1.iter().zip(&k1).map(|(fi, ki)| fi - 2.0 * ki).collect();
        let k2 = crate::kernels::dense::solve_linear_rows(&w, &rhs).ok()?;
        let next: Vec<f64> = y
            .iter()
            .zip(&k1)
            .zip(&k2)
            .map(|((yi, a), b)| yi + 1.5 * h * a + 0.5 * h * b)
            .collect();
        // Error against the first-order (linearly implicit Euler) solution.
        let error = next
            .iter()
            .zip(&y_mid)
            .zip(&y)
            .map(|((a, b), y0)| ((a - b) / (tolerance * (1.0 + y0.abs().max(a.abs())))).powi(2))
            .sum::<f64>()
            .sqrt()
            / (n as f64).sqrt();
        if error <= 1.0 || !error.is_finite() && h.abs() < 1e-10 {
            t += h;
            y = next;
            out.push((t, y.clone()));
        }
        let factor = if error.is_finite() {
            (0.9 / error.max(1e-10).sqrt()).clamp(0.2, 5.0)
        } else {
            0.2
        };
        h *= factor;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symplectic_energy_stays_bounded() {
        // Kepler orbit with eccentricity 0.5: energy error bounded over
        // 100 periods for both methods, and far smaller at order 4.
        let accel = |q: &[f64]| {
            let r3 = (q[0] * q[0] + q[1] * q[1]).powf(1.5);
            vec![-q[0] / r3, -q[1] / r3]
        };
        let energy = |q: &[f64], p: &[f64]| 0.5 * p[0].hypot(p[1]).powi(2) - 1.0 / q[0].hypot(q[1]);
        let (q0, p0) = ([0.5, 0.0], [0.0, 3.0_f64.sqrt()]);
        let e0 = energy(&q0, &p0);
        let period = 2.0 * std::f64::consts::PI; // a = 1
        let dt = period / 400.0;
        let steps = 400 * 100;
        let mut worst = [0.0_f64; 2];
        for (k, order) in [2_u32, 4].into_iter().enumerate() {
            let path = symplectic_trajectory(&q0, &p0, 1.0, dt, steps, order, accel);
            for (q, p) in path.iter().step_by(97) {
                worst[k] = worst[k].max((energy(q, p) - e0).abs());
            }
            // At order 4 the orbit still closes after 100 periods (Verlet's
            // orbit precesses).
            if order == 4 {
                let (q, _) = path.last().unwrap();
                assert!((q[0] - 0.5).abs() < 1e-3 && q[1].abs() < 1e-2, "{q:?}");
            }
        }
        assert!(worst[0] < 2e-2 && worst[1] < 1e-5, "{worst:?}");
    }

    #[test]
    fn stiff_problems() {
        // Robertson's chemical kinetics: very stiff; mass is conserved and
        // the known values at t = 40 are matched.
        let robertson = |_: f64, y: &[f64]| {
            vec![
                -0.04 * y[0] + 1e4 * y[1] * y[2],
                0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1],
                3e7 * y[1] * y[1],
            ]
        };
        let path = rosenbrock_solve(robertson, 0.0, &[1.0, 0.0, 0.0], 40.0, 1e-6).unwrap();
        let (_, y) = path.last().unwrap();
        assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-6, "{y:?}");
        assert!(
            (y[0] - 0.715_8).abs() < 2e-3 && (y[2] - 0.284_1).abs() < 2e-3,
            "{y:?}"
        );
        assert!(path.len() < 2000, "{} steps", path.len());
        // y' = -1000 (y - cos t): backward Euler and implicit midpoint stay
        // stable with h = 0.1 (explicit Euler would explode).
        let f = |t: f64, y: &[f64]| vec![-1000.0 * (y[0] - t.cos())];
        let (mut a, mut b) = (vec![0.0], vec![0.0]);
        for k in 0..50 {
            let t = 0.1 * f64::from(k);
            a = backward_euler_step(&f, t, &a, 0.1).unwrap();
            b = implicit_midpoint_step(&f, t, &b, 0.1).unwrap();
        }
        assert!(
            (a[0] - 5.0_f64.cos()).abs() < 1e-2 && b[0].abs() < 2.0,
            "{a:?} {b:?}"
        );
    }
}
