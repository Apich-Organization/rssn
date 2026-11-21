//! # Runge-Kutta ODE Solvers
//!
//! This module provides a collection of Runge-Kutta methods for solving systems of first-order
//! Ordinary Differential Equations (ODEs). It includes both fixed-step and adaptive step-size
//! algorithms to cater to different accuracy and performance requirements.
//!
//! # Overview
//!
//! Runge-Kutta methods are a family of implicit and explicit iterative methods used in temporal discretization
//! for the approximate solutions of ordinary differential equations. This module focuses on explicit methods,
//! which are suitable for non-stiff problems.
//!
//! Key components include:
//! - **ODE System Interface**: The `OdeSystem` trait allows users to define their own differential equations.
//! - **Fixed-Step Solvers**: `solve_rk4` implements the classic 4th-order Runge-Kutta method.
//! - **Adaptive Solvers**:
//!     - `DormandPrince54`: A high-order method (5th order with 4th order error estimate), often the default choice for general purposes.
//!     - `CashKarp45`: A 5th order method with a 4th order embedded error estimate.
//!     - `BogackiShampine23`: A 3rd order method with a 2nd order error estimate, efficient for lower accuracy needs.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/rkm_solver_comparison.png)
//!
//! ## Coverage
//!
//! | Solver | Order | Dimension | Boundary / initial conditions |
//! |---|---|---|---|
//! | `solve_rk4` | 4 | any (`OdeSystem::dim`) | initial value problems |
//! | `DormandPrince54`, `CashKarp45`, `BogackiShampine23` | 5, 5, 3 | any | initial value problems, adaptive step |
//!
//! Stiff problems: [`crate::sim::integrators`]. Hamiltonian problems with
//! geometric structure: [`crate::sim::symplectic`]. Boundary value problems
//! are not handled by this module.

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::par::for_chunks_mut;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::f64_scratch;

/// Defines the interface for a system of first-order ODEs: dy/dt = f(t, y).
pub trait OdeSystem: Sync + Send {
    /// The dimension of the system (number of equations).
    fn dim(&self) -> usize;

    /// Evaluates the function f(t, y) and stores the result in `dy`.
    fn eval(
        &self,
        t: f64,
        y: &[f64],
        dy: &mut [f64],
    );
}

/// Solves an ODE system using the classic 4th-order Runge-Kutta method with a fixed step size.
///
/// This method is a widely used, robust, and relatively accurate explicit method for
/// approximating the solutions of ordinary differential equations.
///
/// # Arguments
/// * `system` - The ODE system to solve, implementing the `OdeSystem` trait.
/// * `y0` - The initial state vector.
/// * `t_span` - A tuple `(t_start, t_end)` specifying the time interval.
/// * `dt` - The fixed time step.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution at each time step.
pub fn solve_rk4<S: OdeSystem + Sync>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Vec<(f64, Vec<f64>)> {
    let (t_start, t_end) = t_span;

    let steps: usize = (((t_end - t_start) / dt).ceil() as i64)
        .try_into()
        .unwrap_or(0);

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let dim = system.dim();

    let mut k1 = vec![0.0; dim];

    let mut k2 = vec![0.0; dim];

    let mut k3 = vec![0.0; dim];

    let mut k4 = vec![0.0; dim];

    let mut y_temp = vec![0.0; dim];

    for _ in 0..steps {
        let mut current_dt = dt;

        if t + current_dt > t_end {
            current_dt = t_end - t;
        }

        if current_dt <= 0.0 {
            break;
        }

        system.eval(t, &y, &mut k1);

        // y_temp = y + 0.5 * current_dt * k1
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k1)
            .for_each(|((yt, &yi), &k1i)| {
                *yt = (0.5 * current_dt).mul_add(k1i, yi);
            });

        system.eval(0.5f64.mul_add(current_dt, t), &y_temp, &mut k2);

        // y_temp = y + 0.5 * current_dt * k2
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k2)
            .for_each(|((yt, &yi), &k2i)| {
                *yt = (0.5 * current_dt).mul_add(k2i, yi);
            });

        system.eval(0.5f64.mul_add(current_dt, t), &y_temp, &mut k3);

        // y_temp = y + current_dt * k3
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k3)
            .for_each(|((yt, &yi), &k3i)| {
                *yt = current_dt.mul_add(k3i, yi);
            });

        system.eval(t + current_dt, &y_temp, &mut k4);

        // y = y + (current_dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        y.par_iter_mut()
            .zip(&k1)
            .zip(&k2)
            .zip(&k3)
            .zip(&k4)
            .for_each(|((((yi, &k1i), &k2i), &k3i), &k4i)| {
                *yi += (current_dt / 6.0) * (2.0f64.mul_add(k3i, 2.0f64.mul_add(k2i, k1i)) + k4i);
            });

        t += current_dt;

        history.push((t, y.clone()));
    }

    history
}

/// Adaptive Runge-Kutta solver using Dormand-Prince 5(4).
pub struct DormandPrince54 {
    c: [f64; 7],
    a: [[f64; 6]; 6],
    b5: [f64; 7],
    b4: [f64; 7],
}

impl Default for DormandPrince54 {
    fn default() -> Self {
        Self::new()
    }
}

impl DormandPrince54 {
    /// Creates a new Dormand-Prince 5(4) solver.
    #[must_use]
    pub fn new() -> Self {
        Self {
            c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0],
            a: [
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
            ],
            b5: [
                35.0 / 384.0,
                0.0,
                500.0 / 1113.0,
                125.0 / 192.0,
                -2187.0 / 6784.0,
                11.0 / 84.0,
                0.0,
            ],
            b4: [
                5179.0 / 57600.0,
                0.0,
                7571.0 / 16695.0,
                393.0 / 640.0,
                -92097.0 / 339_200.0,
                187.0 / 2100.0,
                1.0 / 40.0,
            ],
        }
    }

    /// Solves an ODE system forward in time using the Dormand-Prince 5(4) method.
    ///
    /// The output includes the initial state and every accepted step. The
    /// last time is exactly the requested endpoint when it is reached.
    /// Invalid controls, non-finite rates, or a step that cannot advance
    /// return the accepted prefix. Tolerances are nonnegative and at
    /// least one must be positive.
    pub fn solve<S: OdeSystem + Sync>(
        &self,
        system: &S,
        y0: &[f64],
        t_span: (f64, f64),
        mut dt: f64,
        tol: (f64, f64),
    ) -> Vec<(f64, Vec<f64>)> {
        let (mut t, t_end) = t_span;
        let (rtol, atol) = tol;
        let mut history = vec![(t, y0.to_vec())];
        if !t.is_finite()
            || !t_end.is_finite()
            || t >= t_end
            || !dt.is_finite()
            || dt <= 0.0
            || !rtol.is_finite()
            || !atol.is_finite()
            || rtol < 0.0
            || atol < 0.0
            || (rtol == 0.0 && atol == 0.0)
            || y0.is_empty()
            || y0.len() != system.dim()
            || y0.iter().any(|v| !v.is_finite())
        {
            return history;
        }

        // Keep the ordinary minimum step and scale it down for short spans.
        let minimum_dt = 1e-12 * (t_end - t).min(1.0);
        let mut y = y0.to_vec();
        let mut work = DormandWorkspace::new(y.len());
        while t < t_end {
            let mut next = t + dt;
            if next > t_end {
                next = t_end;
                dt = t_end - t;
            }
            if !dt.is_finite() || dt <= 0.0 || next <= t {
                break;
            }

            let Some(error) = work.attempt(self, system, &y, (t, next, dt), tol) else {
                break;
            };
            if error <= 1.0 {
                t = next;
                std::mem::swap(&mut y, &mut work.high_order);
                history.push((t, y.clone()));
            }

            let factor = (0.9 * (1.0 / error).powf(0.2)).clamp(0.1, 4.0);
            dt *= factor;
            if dt == 0.0 || dt < minimum_dt {
                break;
            }
        }

        history
    }
}

const DORMAND_BLOCK: usize = 512;
const DORMAND_PARALLEL_MIN: usize = 16_384;

struct DormandWorkspace {
    derivatives: [Vec<f64>; 7],
    stage: Scratch<'static, f64>,
    high_order: Vec<f64>,
}

impl DormandWorkspace {
    fn new(dim: usize) -> Self {
        let mut stage = f64_scratch(dim);
        stage.resize(dim, 0.0);
        Self {
            derivatives: std::array::from_fn(|_| vec![0.0; dim]),
            stage,
            high_order: vec![0.0; dim],
        }
    }

    // None indicates an invalid callback. An infinite error is a rejected
    // trial based on finite rates and can be retried with a smaller step.
    fn attempt<S: OdeSystem>(
        &mut self,
        solver: &DormandPrince54,
        system: &S,
        y: &[f64],
        (t, next, dt): (f64, f64, f64),
        (rtol, atol): (f64, f64),
    ) -> Option<f64> {
        system.eval(t, y, &mut self.derivatives[0]);
        if self.derivatives[0].iter().any(|d| !d.is_finite()) {
            return None;
        }
        let parallel_grain = if y.len() >= DORMAND_PARALLEL_MIN {
            2
        } else {
            usize::MAX
        };
        for i in 1..7 {
            let (earlier, current) = self.derivatives.split_at_mut(i);
            let weights = &solver.a[i - 1];
            let scaled = weights.map(|a| dt * a);
            for_chunks_mut(
                &mut self.stage,
                DORMAND_BLOCK,
                parallel_grain,
                |block, out| {
                    let start = block * DORMAND_BLOCK;
                    let end = start + out.len();
                    out.copy_from_slice(&y[start..end]);
                    for ((&a, &w), rate) in weights.iter().zip(&scaled).zip(earlier.iter()) {
                        if a == 0.0 {
                            continue;
                        }
                        let rate = &rate[start..end];
                        if w.is_normal() {
                            for (v, &k) in out.iter_mut().zip(rate) {
                                *v += w * k;
                            }
                        } else {
                            for (v, &k) in out.iter_mut().zip(rate) {
                                *v += (dt * k) * a;
                            }
                        }
                    }
                },
            );
            if self.stage.iter().any(|value| !value.is_finite()) {
                return Some(f64::INFINITY);
            }
            let c = solver.c[i];
            let stage_time = if c >= 1.0 {
                next
            } else {
                c.mul_add(dt, t).min(next)
            };
            system.eval(stage_time, &self.stage, &mut current[0]);
            if current[0].iter().any(|d| !d.is_finite()) {
                return None;
            }
        }

        let mut squared_error = 0.0;
        // Fill both embedded estimates in cache-sized blocks, then sum the
        // error ratios in the original global component order.
        for (block, (high, low)) in self
            .high_order
            .chunks_mut(DORMAND_BLOCK)
            .zip(self.stage.chunks_mut(DORMAND_BLOCK))
            .enumerate()
        {
            let start = block * DORMAND_BLOCK;
            let end = start + high.len();
            let state = &y[start..end];
            high.copy_from_slice(state);
            low.copy_from_slice(state);
            for (j, rate) in self.derivatives.iter().enumerate() {
                let (b5, b4) = (solver.b5[j], solver.b4[j]);
                for ((high, low), &k) in high.iter_mut().zip(low.iter_mut()).zip(&rate[start..end])
                {
                    let increment = dt * k;
                    *high += increment * b5;
                    *low += increment * b4;
                }
            }
            for ((&old, &high), &low) in state.iter().zip(high.iter()).zip(low.iter()) {
                if !high.is_finite() || !low.is_finite() {
                    return Some(f64::INFINITY);
                }
                squared_error += dormand_error_ratio(old, high, low, rtol, atol).powi(2);
            }
        }

        Some((squared_error / y.len() as f64).sqrt())
    }
}

fn dormand_error_ratio(
    old: f64,
    high: f64,
    low: f64,
    rtol: f64,
    atol: f64,
) -> f64 {
    let amplitude = old.abs().max(high.abs());
    let scale = amplitude.mul_add(rtol, atol);
    let difference = high - low;
    if difference == 0.0 {
        return 0.0;
    }
    if difference.is_finite() && scale.is_finite() && scale > 0.0 {
        return difference / scale;
    }

    // Rescale by finite state magnitudes when either subtraction or the
    // tolerance scale overflows, or the tolerance product underflows.
    let reference = amplitude.max(low.abs());
    let normalized_scale = (amplitude / reference).mul_add(rtol, atol / reference);
    let normalized_difference = if difference.is_finite() {
        difference / reference
    } else {
        high / reference - low / reference
    };
    normalized_difference / normalized_scale
}

/// Adaptive Runge-Kutta solver using Cash-Karp 4(5).
pub struct CashKarp45 {
    c: [f64; 6],
    a: [[f64; 5]; 5],
    b5: [f64; 6],
    b4: [f64; 6],
}

impl Default for CashKarp45 {
    fn default() -> Self {
        Self {
            c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0],
            a: [
                [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0],
                [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0],
                [3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0.0, 0.0],
                [-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0.0],
                [
                    1631.0 / 55296.0,
                    175.0 / 512.0,
                    575.0 / 13824.0,
                    44275.0 / 110_592.0,
                    253.0 / 4096.0,
                ],
            ],
            b5: [
                37.0 / 378.0,
                0.0,
                250.0 / 621.0,
                125.0 / 594.0,
                0.0,
                512.0 / 1771.0,
            ],
            b4: [
                2825.0 / 27648.0,
                0.0,
                18575.0 / 48384.0,
                13525.0 / 55296.0,
                277.0 / 14336.0,
                1.0 / 4.0,
            ],
        }
    }
}

impl CashKarp45 {
    /// Solves an ODE system using the Cash-Karp 4(5) method.
    pub fn solve<S: OdeSystem + Sync>(
        &self,
        system: &S,
        y0: &[f64],
        t_span: (f64, f64),
        mut dt: f64,
        tol: (f64, f64),
    ) -> Vec<(f64, Vec<f64>)> {
        let (t_start, t_end) = t_span;

        let (rtol, atol) = tol;

        let mut t = t_start;

        let mut y = y0.to_vec();

        let mut history = vec![(t, y.clone())];

        let dim = system.dim();

        let mut k = vec![vec![0.0; dim]; 6];

        loop {
            if (t_end - t).abs() <= 1e-15 {
                break;
            }

            if t + dt > t_end {
                dt = t_end - t;
            }

            system.eval(t, &y, &mut k[0]);

            for i in 1..6 {
                let mut y_temp = y.clone();

                for (j, kj) in k.iter().enumerate().take(i) {
                    let a_val = self.a[i - 1][j];

                    if a_val != 0.0 {
                        y_temp.par_iter_mut().zip(kj).for_each(|(yt, &kj_i)| {
                            *yt += dt * a_val * kj_i;
                        });
                    }
                }

                system.eval(self.c[i].mul_add(dt, t), &y_temp, &mut k[i]);
            }

            let mut error = 0.0;

            for i in 0..dim {
                let mut y5_i = y[i];

                let mut y4_i = y[i];

                for (j, kj) in k.iter().enumerate().take(6) {
                    y5_i += dt * kj[i] * self.b5[j];

                    y4_i += dt * kj[i] * self.b4[j];
                }

                let scale = y[i].abs().max(y5_i.abs()).mul_add(rtol, atol);

                error += ((y5_i - y4_i) / scale).powi(2);
            }

            error = (error / dim as f64).sqrt();
            // Invalid stages cannot produce an accepted path. In particular a
            // bounded host oracle returns NaN after its first failure; allowing
            // NaN step sizes here would otherwise leave the loop running forever.
            if !error.is_finite() || !dt.is_finite() || dt <= 0.0 || t + dt <= t {
                break;
            }

            let factor = (0.9 * (1.0 / error).powf(0.2)).clamp(0.1, 4.0);

            if error <= 1.0 {
                t += dt;

                y.par_iter_mut().enumerate().for_each(|(i, yi)| {
                    for (j, kj) in k.iter().enumerate().take(6) {
                        *yi += dt * kj[i] * self.b5[j];
                    }
                });

                history.push((t, y.clone()));
            }

            dt *= factor;

            if dt < 1e-12 {
                break;
            }
        }

        history
    }
}

/// Adaptive Runge-Kutta solver using Bogacki-Shampine 3(2).
/// Efficient for low-accuracy requirements.
pub struct BogackiShampine23 {
    c: [f64; 4],
    a: [[f64; 3]; 3],
    b3: [f64; 4],
    b2: [f64; 4],
}

impl Default for BogackiShampine23 {
    fn default() -> Self {
        Self {
            c: [0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0],
            a: [
                [1.0 / 2.0, 0.0, 0.0],
                [0.0, 3.0 / 4.0, 0.0],
                [2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0],
            ],
            b3: [2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0],
            b2: [7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0],
        }
    }
}

impl BogackiShampine23 {
    /// Solves an ODE system using the Bogacki-Shampine 2(3) method.
    pub fn solve<S: OdeSystem + Sync>(
        &self,
        system: &S,
        y0: &[f64],
        t_span: (f64, f64),
        mut dt: f64,
        tol: (f64, f64),
    ) -> Vec<(f64, Vec<f64>)> {
        let (t_start, t_end) = t_span;

        let (rtol, atol) = tol;

        let mut t = t_start;

        let mut y = y0.to_vec();

        let mut history = vec![(t, y.clone())];

        let dim = system.dim();

        let mut k = vec![vec![0.0; dim]; 4];

        loop {
            if (t_end - t).abs() <= 1e-15 {
                break;
            }

            if t + dt > t_end {
                dt = t_end - t;
            }

            system.eval(t, &y, &mut k[0]);

            for i in 1..4 {
                let mut y_temp = y.clone();

                for (j, kj) in k.iter().enumerate().take(i) {
                    let a_val = self.a[i - 1][j];

                    if a_val != 0.0 {
                        y_temp.par_iter_mut().zip(kj).for_each(|(yt, &kj_i)| {
                            *yt += dt * a_val * kj_i;
                        });
                    }
                }

                system.eval(self.c[i].mul_add(dt, t), &y_temp, &mut k[i]);
            }

            let mut error = 0.0;

            for i in 0..dim {
                let mut y3_i = y[i];

                let mut y2_i = y[i];

                for (j, kj) in k.iter().enumerate().take(4) {
                    y3_i += dt * kj[i] * self.b3[j];

                    y2_i += dt * kj[i] * self.b2[j];
                }

                let scale = y[i].abs().max(y3_i.abs()).mul_add(rtol, atol);

                error += ((y3_i - y2_i) / scale).powi(2);
            }

            error = (error / dim as f64).sqrt();
            // Invalid stages cannot produce an accepted path. In particular a
            // bounded host oracle returns NaN after its first failure; allowing
            // NaN step sizes here would otherwise leave the loop running forever.
            if !error.is_finite() || !dt.is_finite() || dt <= 0.0 || t + dt <= t {
                break;
            }

            let factor = (0.9 * (1.0 / error).powf(0.33)).clamp(0.1, 4.0);

            if error <= 1.0 {
                t += dt;

                y.par_iter_mut().enumerate().for_each(|(i, yi)| {
                    for (j, kj) in k.iter().enumerate().take(4) {
                        *yi += dt * kj[i] * self.b3[j];
                    }
                });

                history.push((t, y.clone()));
            }

            dt *= factor;

            if dt < 1e-12 {
                break;
            }
        }

        history
    }
}

// ============================================================================
// ODE Systems
// ============================================================================

/// The Lorenz attractor system.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LorenzSystem {
    /// The sigma parameter.
    pub sigma: f64,
    /// The rho parameter.
    pub rho: f64,
    /// The beta parameter.
    pub beta: f64,
}

impl OdeSystem for LorenzSystem {
    fn dim(&self) -> usize {
        3
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = self.sigma * (y[1] - y[0]);

        dy[1] = y[0].mul_add(self.rho - y[2], -y[1]);

        dy[2] = y[0].mul_add(y[1], -(self.beta * y[2]));
    }
}

/// A damped harmonic oscillator (y'' + 2*zeta*omega*y' + omega^2*y = 0).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DampedOscillatorSystem {
    /// The angular frequency.
    pub omega: f64,
    /// The damping ratio.
    pub zeta: f64,
}

impl OdeSystem for DampedOscillatorSystem {
    fn dim(&self) -> usize {
        2
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = y[1];

        dy[1] = (-2.0 * self.zeta * self.omega).mul_add(y[1], -(self.omega.powi(2) * y[0]));
    }
}

/// Van der Pol oscillator system.
/// y'' - mu(1 - y^2)y' + y = 0
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VanDerPolSystem {
    /// The mu parameter.
    pub mu: f64,
}

impl OdeSystem for VanDerPolSystem {
    fn dim(&self) -> usize {
        2
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = y[1];

        dy[1] = (self.mu * y[0].mul_add(-y[0], 1.0)).mul_add(y[1], -y[0]);
    }
}

/// Lotka-Volterra predator-prey system.
/// dx/dt = alpha*x - beta*x*y
/// dy/dt = delta*x*y - gamma*y
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LotkaVolterraSystem {
    /// The alpha parameter.
    pub alpha: f64,
    /// The beta parameter.
    pub beta: f64,
    /// The delta parameter.
    pub delta: f64,
    /// The gamma parameter.
    pub gamma: f64,
}

impl OdeSystem for LotkaVolterraSystem {
    fn dim(&self) -> usize {
        2
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = self.alpha.mul_add(y[0], -(self.beta * y[0] * y[1]));

        dy[1] = (self.delta * y[0]).mul_add(y[1], -(self.gamma * y[1]));
    }
}

/// Simple pendulum system.
/// theta'' + (g/L)sin(theta) = 0
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PendulumSystem {
    /// The acceleration due to gravity.
    pub g: f64,
    /// The length of the pendulum.
    pub l: f64,
}

impl OdeSystem for PendulumSystem {
    fn dim(&self) -> usize {
        2
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = y[1];

        dy[1] = -(self.g / self.l) * y[0].sin();
    }
}

// ============================================================================
// Scenarios
// ============================================================================

/// Simulates the Lorenz attractor system.
#[must_use]
pub fn simulate_lorenz_attractor_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = LorenzSystem {
        sigma: 10.0,
        rho: 28.0,
        beta: 8.0 / 3.0,
    };

    let y0 = &[1.0, 1.0, 1.0];

    let t_span = (0.0, 50.0);

    let dt_initial = 0.01;

    let tolerance = (1e-6, 1e-6);

    let solver = DormandPrince54::new();

    solver.solve(&system, y0, t_span, dt_initial, tolerance)
}

/// Simulates a damped harmonic oscillator.
#[must_use]
pub fn simulate_damped_oscillator_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = DampedOscillatorSystem {
        omega: 1.0,
        zeta: 0.15,
    };

    let y0 = &[1.0, 0.0];

    let t_span = (0.0, 40.0);

    let dt = 0.1;

    solve_rk4(&system, y0, t_span, dt)
}

/// Simulates the Van der Pol oscillator.
#[must_use]
pub fn simulate_vanderpol_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = VanDerPolSystem { mu: 1.0 };

    let y0 = &[2.0, 0.0];

    let t_span = (0.0, 20.0);

    let dt_initial = 0.1;

    let tolerance = (1e-6, 1e-6);

    let solver = CashKarp45::default();

    solver.solve(&system, y0, t_span, dt_initial, tolerance)
}

/// Simulates the Lotka-Volterra predator-prey system.
#[must_use]
pub fn simulate_lotka_volterra_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = LotkaVolterraSystem {
        alpha: 1.5,
        beta: 1.0,
        delta: 1.0,
        gamma: 3.0,
    };

    let y0 = &[10.0, 5.0];

    let t_span = (0.0, 15.0);

    let dt_initial = 0.01;

    let tolerance = (1e-6, 1e-6);

    let solver = BogackiShampine23::default();

    solver.solve(&system, y0, t_span, dt_initial, tolerance)
}
