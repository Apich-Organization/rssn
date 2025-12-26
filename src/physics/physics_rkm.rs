use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Defines the interface for a system of first-order ODEs: dy/dt = f(t, y).
pub trait OdeSystem: Sync + Send {
    /// The dimension of the system (number of equations).
    fn dim(&self) -> usize;
    /// Evaluates the function f(t, y) and stores the result in `dy`.
    fn eval(&self, t: f64, y: &[f64], dy: &mut [f64]);
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
    let steps = ((t_end - t_start) / dt).ceil() as usize;
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

    while t < t_end {
        let mut current_dt = dt;
        if t + current_dt > t_end {
            current_dt = t_end - t;
        }

        system.eval(t, &y, &mut k1);

        // y_temp = y + 0.5 * current_dt * k1
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k1)
            .for_each(|((yt, &yi), &k1i)| {
                *yt = yi + 0.5 * current_dt * k1i;
            });

        system.eval(t + 0.5 * current_dt, &y_temp, &mut k2);

        // y_temp = y + 0.5 * current_dt * k2
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k2)
            .for_each(|((yt, &yi), &k2i)| {
                *yt = yi + 0.5 * current_dt * k2i;
            });

        system.eval(t + 0.5 * current_dt, &y_temp, &mut k3);

        // y_temp = y + current_dt * k3
        y_temp
            .par_iter_mut()
            .zip(&y)
            .zip(&k3)
            .for_each(|((yt, &yi), &k3i)| {
                *yt = yi + current_dt * k3i;
            });

        system.eval(t + current_dt, &y_temp, &mut k4);

        // y = y + (current_dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        y.par_iter_mut()
            .zip(&k1)
            .zip(&k2)
            .zip(&k3)
            .zip(&k4)
            .for_each(|((((yi, &k1i), &k2i), &k3i), &k4i)| {
                *yi += (current_dt / 6.0) * (k1i + 2.0 * k2i + 2.0 * k3i + k4i);
            });

        t += current_dt;
        history.push((t, y.clone()));
    }
    history
}

/// Adaptive Runge-Kutta solver using Dormand-Prince 5(4).
#[derive(Default)]
pub struct DormandPrince54 {
    c: [f64; 7],
    a: [[f64; 6]; 6],
    b5: [f64; 7],
    b4: [f64; 7],
}

impl DormandPrince54 {
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
                -92097.0 / 339200.0,
                187.0 / 2100.0,
                1.0 / 40.0,
            ],
        }
    }

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
        let mut k = vec![vec![0.0; dim]; 7];

        while t < t_end {
            if t + dt > t_end {
                dt = t_end - t;
            }

            system.eval(t, &y, &mut k[0]);

            for i in 1..7 {
                let mut y_temp = y.clone();
                for (j, _vars) in k.iter().enumerate().take(i) {
                    let a_val = self.a[i - 1][j];
                    if a_val != 0.0 {
                        y_temp.par_iter_mut().zip(&k[j]).for_each(|(yt, &kj)| {
                            *yt += dt * a_val * kj;
                        });
                    }
                }
                system.eval(t + self.c[i] * dt, &y_temp, &mut k[i]);
            }

            let mut error = 0.0;
            for i in 0..dim {
                let mut y5_i = y[i];
                let mut y4_i = y[i];
                for j in 0..7 {
                    y5_i += dt * k[j][i] * self.b5[j];
                    y4_i += dt * k[j][i] * self.b4[j];
                }
                let scale = atol + y[i].abs().max(y5_i.abs()) * rtol;
                error += ((y5_i - y4_i) / scale).powi(2);
            }
            error = (error / dim as f64).sqrt();

            let factor = (0.9 * (1.0 / error).powf(0.2)).clamp(0.1, 4.0);
            if error <= 1.0 {
                t += dt;
                y.par_iter_mut().enumerate().for_each(|(i, yi)| {
                    for j in 0..7 {
                        *yi += dt * k[j][i] * self.b5[j];
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
                    44275.0 / 110592.0,
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

        while t < t_end {
            if t + dt > t_end {
                dt = t_end - t;
            }

            system.eval(t, &y, &mut k[0]);

            for i in 1..6 {
                let mut y_temp = y.clone();
                for j in 0..i {
                    let a_val = self.a[i - 1][j];
                    if a_val != 0.0 {
                        y_temp.par_iter_mut().zip(&k[j]).for_each(|(yt, &kj)| {
                            *yt += dt * a_val * kj;
                        });
                    }
                }
                system.eval(t + self.c[i] * dt, &y_temp, &mut k[i]);
            }

            let mut error = 0.0;
            for i in 0..dim {
                let mut y5_i = y[i];
                let mut y4_i = y[i];
                for j in 0..6 {
                    y5_i += dt * k[j][i] * self.b5[j];
                    y4_i += dt * k[j][i] * self.b4[j];
                }
                let scale = atol + y[i].abs().max(y5_i.abs()) * rtol;
                error += ((y5_i - y4_i) / scale).powi(2);
            }
            error = (error / dim as f64).sqrt();

            let factor = (0.9 * (1.0 / error).powf(0.2)).clamp(0.1, 4.0);
            if error <= 1.0 {
                t += dt;
                y.par_iter_mut().enumerate().for_each(|(i, yi)| {
                    for j in 0..6 {
                        *yi += dt * k[j][i] * self.b5[j];
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

        while t < t_end {
            if t + dt > t_end {
                dt = t_end - t;
            }

            system.eval(t, &y, &mut k[0]);

            for i in 1..4 {
                let mut y_temp = y.clone();
                for j in 0..i {
                    let a_val = self.a[i - 1][j];
                    if a_val != 0.0 {
                        y_temp.par_iter_mut().zip(&k[j]).for_each(|(yt, &kj)| {
                            *yt += dt * a_val * kj;
                        });
                    }
                }
                system.eval(t + self.c[i] * dt, &y_temp, &mut k[i]);
            }

            let mut error = 0.0;
            for i in 0..dim {
                let mut y3_i = y[i];
                let mut y2_i = y[i];
                for j in 0..4 {
                    y3_i += dt * k[j][i] * self.b3[j];
                    y2_i += dt * k[j][i] * self.b2[j];
                }
                let scale = atol + y[i].abs().max(y3_i.abs()) * rtol;
                error += ((y3_i - y2_i) / scale).powi(2);
            }
            error = (error / dim as f64).sqrt();

            let factor = (0.9 * (1.0 / error).powf(0.33)).clamp(0.1, 4.0);
            if error <= 1.0 {
                t += dt;
                y.par_iter_mut().enumerate().for_each(|(i, yi)| {
                    for j in 0..4 {
                        *yi += dt * k[j][i] * self.b3[j];
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
    pub sigma: f64,
    pub rho: f64,
    pub beta: f64,
}

impl OdeSystem for LorenzSystem {
    fn dim(&self) -> usize {
        3
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.sigma * (y[1] - y[0]);
        dy[1] = y[0] * (self.rho - y[2]) - y[1];
        dy[2] = y[0] * y[1] - self.beta * y[2];
    }
}

/// A damped harmonic oscillator (y'' + 2*zeta*omega*y' + omega^2*y = 0).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DampedOscillatorSystem {
    pub omega: f64,
    pub zeta: f64,
}

impl OdeSystem for DampedOscillatorSystem {
    fn dim(&self) -> usize {
        2
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = -2.0 * self.zeta * self.omega * y[1] - self.omega.powi(2) * y[0];
    }
}

/// Van der Pol oscillator system.
/// y'' - mu(1 - y^2)y' + y = 0
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VanDerPolSystem {
    pub mu: f64,
}

impl OdeSystem for VanDerPolSystem {
    fn dim(&self) -> usize {
        2
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = self.mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
    }
}

/// Lotka-Volterra predator-prey system.
/// dx/dt = alpha*x - beta*x*y
/// dy/dt = delta*x*y - gamma*y
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LotkaVolterraSystem {
    pub alpha: f64,
    pub beta: f64,
    pub delta: f64,
    pub gamma: f64,
}

impl OdeSystem for LotkaVolterraSystem {
    fn dim(&self) -> usize {
        2
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.alpha * y[0] - self.beta * y[0] * y[1];
        dy[1] = self.delta * y[0] * y[1] - self.gamma * y[1];
    }
}

/// Simple pendulum system.
/// theta'' + (g/L)sin(theta) = 0
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PendulumSystem {
    pub g: f64,
    pub l: f64,
}

impl OdeSystem for PendulumSystem {
    fn dim(&self) -> usize {
        2
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = -(self.g / self.l) * y[0].sin();
    }
}

// ============================================================================
// Scenarios
// ============================================================================

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

pub fn simulate_vanderpol_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = VanDerPolSystem { mu: 1.0 };
    let y0 = &[2.0, 0.0];
    let t_span = (0.0, 20.0);
    let dt_initial = 0.1;
    let tolerance = (1e-6, 1e-6);
    let solver = CashKarp45::default();
    solver.solve(&system, y0, t_span, dt_initial, tolerance)
}

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
