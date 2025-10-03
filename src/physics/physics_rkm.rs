// src/physics/physics_rkm.rs
// A module for solving ordinary differential equations (ODEs) using Runge-Kutta methods.

//use rayon::prelude::*;

/// Defines the interface for a system of first-order ODEs: dy/dt = f(t, y).
pub trait OdeSystem {
    /// The dimension of the system (number of equations).
    fn dim(&self) -> usize;
    /// Evaluates the function f(t, y) and stores the result in `dy`.
    fn eval(&self, t: f64, y: &[f64], dy: &mut [f64]);
}

// --- Fixed-Step 4th-Order Runge-Kutta (RK4) Solver ---

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

    for _ in 0..steps {
        // k1 = f(t, y)
        system.eval(t, &y, &mut k1);

        // k2 = f(t + dt/2, y + dt/2 * k1)
        y_temp
            .iter_mut()
            .zip(&y)
            .zip(&k1)
            .for_each(|((yt, &yi), &k1i)| *yt = yi + 0.5 * dt * k1i);
        system.eval(t + 0.5 * dt, &y_temp, &mut k2);

        // k3 = f(t + dt/2, y + dt/2 * k2)
        y_temp
            .iter_mut()
            .zip(&y)
            .zip(&k2)
            .for_each(|((yt, &yi), &k2i)| *yt = yi + 0.5 * dt * k2i);
        system.eval(t + 0.5 * dt, &y_temp, &mut k3);

        // k4 = f(t + dt, y + dt * k3)
        y_temp
            .iter_mut()
            .zip(&y)
            .zip(&k3)
            .for_each(|((yt, &yi), &k3i)| *yt = yi + dt * k3i);
        system.eval(t + dt, &y_temp, &mut k4);

        // y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..dim {
            y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        t += dt;
        history.push((t, y.clone()));
    }

    history
}

// --- Adaptive-Step Dormand-Prince (RKDP) Solver ---

/// An adaptive step-size solver using the Dormand-Prince 5(4) pair (also known as `ode45`).
pub struct DormandPrince54 {
    // Butcher Tableau coefficients
    c: [f64; 7],
    a: [[f64; 6]; 6],
    b5: [f64; 7], // 5th order solution
    b4: [f64; 7], // 4th order solution (for error estimation)
}

impl Default for DormandPrince54 {
    fn default() -> Self {
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
}

impl DormandPrince54 {
    /// Solves an ODE system using the adaptive Dormand-Prince 5(4) Runge-Kutta method.
    ///
    /// This method (often referred to as `ode45`) uses a pair of Runge-Kutta methods
    /// of different orders (5th and 4th) to estimate the local truncation error and
    /// adaptively adjust the step size `dt` to maintain a desired tolerance.
    ///
    /// # Arguments
    /// * `system` - The ODE system to solve, implementing the `OdeSystem` trait.
    /// * `y0` - The initial state vector.
    /// * `t_span` - A tuple `(t_start, t_end)` specifying the time interval.
    /// * `dt` - The initial time step.
    /// * `tol` - A tuple `(rtol, atol)` specifying the relative and absolute tolerances for error control.
    ///
    /// # Returns
    /// A `Vec` of tuples `(time, state_vector)` representing the solution path.
    pub fn solve<S: OdeSystem + Sync>(
        &self,
        system: &S,
        y0: &[f64],
        t_span: (f64, f64),
        mut dt: f64,
        tol: (f64, f64), // (rtol, atol)
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

            // Calculate k stages
            system.eval(t, &y, &mut k[0]);
            for i in 1..7 {
                let mut y_temp = y.clone();
                for j in 0..i {
                    let a_val = self.a[i - 1][j];
                    if a_val != 0.0 {
                        y_temp
                            .iter_mut()
                            .zip(&k[j])
                            .for_each(|(yt, &kj)| *yt += dt * a_val * kj);
                    }
                }
                system.eval(t + self.c[i] * dt, &y_temp, &mut k[i]);
            }

            // Calculate error and new step size
            let mut error = 0.0;
            for i in 0..dim {
                let y5_i = y[i]
                    + dt * k
                        .iter()
                        .map(|ki| ki[i] * self.b5[k.iter().position(|x| x == ki).unwrap()])
                        .sum::<f64>();
                let y4_i = y[i]
                    + dt * k
                        .iter()
                        .map(|ki| ki[i] * self.b4[k.iter().position(|x| x == ki).unwrap()])
                        .sum::<f64>();
                let scale = atol + y[i].abs().max(y5_i.abs()) * rtol;
                error += ((y5_i - y4_i) / scale).powi(2);
            }
            error = (error / dim as f64).sqrt();

            let factor = (0.9 * (1.0 / error).powf(0.2)).min(4.0).max(0.1);

            if error <= 1.0 {
                // Step accepted
                t += dt;
                for i in 0..dim {
                    y[i] += dt
                        * k.iter()
                            .map(|ki| ki[i] * self.b5[k.iter().position(|x| x == ki).unwrap()])
                            .sum::<f64>();
                }
                history.push((t, y.clone()));
            }

            dt *= factor;
        }
        history
    }
}

// --- Example Scenarios ---

/// The Lorenz attractor system.
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
pub struct DampedOscillatorSystem {
    pub omega: f64, // Natural frequency
    pub zeta: f64,  // Damping ratio
}

impl OdeSystem for DampedOscillatorSystem {
    fn dim(&self) -> usize {
        2
    }
    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        // y[0] = position, y[1] = velocity
        dy[0] = y[1];
        dy[1] = -2.0 * self.zeta * self.omega * y[1] - self.omega.powi(2) * y[0];
    }
}

/// Solves the Lorenz attractor system using the adaptive Dormand-Prince solver.
///
/// The Lorenz system is a set of three ordinary differential equations known for its
/// chaotic solutions for certain parameter values. This scenario demonstrates the
/// adaptive solver's ability to handle complex, non-linear dynamics.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.
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

    let solver = DormandPrince54::default();
    solver.solve(&system, y0, t_span, dt_initial, tolerance)
}

/// Solves the damped harmonic oscillator system using the fixed-step RK4 solver.
///
/// This scenario models a mass-spring-damper system, demonstrating the RK4 solver's
/// ability to accurately capture oscillatory behavior with damping.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.
pub fn simulate_damped_oscillator_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = DampedOscillatorSystem {
        omega: 1.0,
        zeta: 0.15,
    };
    let y0 = &[1.0, 0.0]; // Initial position = 1, initial velocity = 0
    let t_span = (0.0, 40.0);
    let dt = 0.1;

    solve_rk4(&system, y0, t_span, dt)
}
