use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::physics::physics_rkm::OdeSystem;

#[derive(Clone, Debug, Serialize, Deserialize)]

pub struct EulerSolverConfig {
    pub dt : f64,
}

/// Solves an ODE system using the forward (explicit) Euler method.
///
/// This method is simple to implement but can be unstable and inaccurate for large step sizes,
/// especially for stiff differential equations. It is a first-order method.
///
/// # Arguments
/// * `system` - The ODE system to solve, implementing the `OdeSystem` trait.
/// * `y0` - The initial state vector.
/// * `t_span` - A tuple `(t_start, t_end)` specifying the time interval.
/// * `dt` - The fixed time step.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.

pub fn solve_forward_euler<S : OdeSystem>(
    system : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
) -> Vec<(f64, Vec<f64>)> {

    let (t_start, t_end) = t_span;

    let steps = ((t_end - t_start) / dt).ceil() as usize;

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let dim = system.dim();

    let mut dy = vec![0.0; dim];

    for _ in 0 .. steps {

        system.eval(t, &y, &mut dy);

        y.par_iter_mut()
            .zip(&dy)
            .for_each(|(yi, &dyi)| {

                *yi += dt * dyi;
            });

        t += dt;

        history.push((t, y.clone()));
    }

    history
}

/// Solves an ODE system using the explicit midpoint method (Modified Euler).

pub fn solve_midpoint_euler<S : OdeSystem>(
    system : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
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

    let mut y_mid = vec![0.0; dim];

    for _ in 0 .. steps {

        system.eval(t, &y, &mut k1);

        y_mid
            .par_iter_mut()
            .zip(&y)
            .zip(&k1)
            .for_each(
                |((ym, &yi), &k1i)| {

                    *ym = yi + 0.5 * dt * k1i;
                },
            );

        system.eval(
            t + 0.5 * dt,
            &y_mid,
            &mut k2,
        );

        y.par_iter_mut()
            .zip(&k2)
            .for_each(|(yi, &k2i)| {

                *yi += dt * k2i;
            });

        t += dt;

        history.push((t, y.clone()));
    }

    history
}

/// Solves an ODE system using Heun's method (Improved Euler).

pub fn solve_heun_euler<S : OdeSystem>(
    system : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
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

    let mut y_predict = vec![0.0; dim];

    for _ in 0 .. steps {

        system.eval(t, &y, &mut k1);

        y_predict
            .par_iter_mut()
            .zip(&y)
            .zip(&k1)
            .for_each(
                |((yp, &yi), &k1i)| {

                    *yp = yi + dt * k1i;
                },
            );

        system.eval(
            t + dt,
            &y_predict,
            &mut k2,
        );

        y.par_iter_mut()
            .zip(&k1)
            .zip(&k2)
            .for_each(
                |((yi, &k1i), &k2i)| {

                    *yi += 0.5 * dt * (k1i + k2i);
                },
            );

        t += dt;

        history.push((t, y.clone()));
    }

    history
}

/// Defines a special kind of ODE system for mechanics: dv/dt = a(x), dx/dt = v.
/// The state vector is partitioned into `[positions..., velocities...]`.

pub trait MechanicalSystem {
    /// The number of spatial dimensions (e.g., 1, 2, 3).

    fn spatial_dim(&self) -> usize;

    /// Evaluates the acceleration `a(x)` given the positions `x`.

    fn eval_acceleration(
        &self,
        x : &[f64],
        a : &mut [f64],
    );
}

/// Solves a second-order ODE system (typically mechanical systems) using the semi-implicit Euler method.
///
/// This method is often used in mechanics for its energy conservation properties (it is symplectic).
/// It is particularly well-suited for systems where the state vector can be partitioned into
/// positions and velocities, and accelerations depend only on positions.
/// It assumes the state vector `y` is ordered as `[x0, x1, ..., v0, v1, ...]`.
///
/// # Arguments
/// * `system` - The mechanical system to solve, implementing the `MechanicalSystem` trait.
/// * `y0` - The initial state `[initial_positions..., initial_velocities...]`.
/// * `t_span` - The time interval `(t_start, t_end)`.
/// * `dt` - The fixed time step.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.

pub fn solve_semi_implicit_euler<S : MechanicalSystem>(
    system : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
) -> Result<Vec<(f64, Vec<f64>)>, String> {

    let s_dim = system.spatial_dim();

    if y0.len() != 2 * s_dim {

        return Err("State vector length must be twice the spatial dimension.".to_string());
    }

    let (t_start, t_end) = t_span;

    let steps = ((t_end - t_start) / dt).ceil() as usize;

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let mut a = vec![0.0; s_dim];

    for _ in 0 .. steps {

        let (x, v) = y.split_at_mut(s_dim);

        system.eval_acceleration(x, &mut a);

        v.par_iter_mut()
            .zip(&a)
            .for_each(|(vi, &ai)| {

                *vi += dt * ai;
            });

        x.par_iter_mut()
            .zip(&*v)
            .for_each(|(xi, &vi)| {

                *xi += dt * vi;
            });

        t += dt;

        history.push((t, y.clone()));
    }

    Ok(history)
}

use crate::physics::physics_rkm::DampedOscillatorSystem;

/// Solves the damped harmonic oscillator with the less stable forward Euler method.
///
/// This scenario demonstrates the numerical energy gain over time that can occur
/// with the forward Euler method when solving oscillatory systems, highlighting
/// its limitations for such problems.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.

pub fn simulate_oscillator_forward_euler_scenario() -> Vec<(f64, Vec<f64>)> {

    let system = DampedOscillatorSystem {
        omega : 2.0 * std::f64::consts::PI,
        zeta : 0.0,
    };

    let y0 = &[1.0, 0.0];

    let t_span = (0.0, 10.0);

    let dt = 0.01;

    solve_forward_euler(
        &system,
        y0,
        t_span,
        dt,
    )
}

/// A simple 2D orbital system (e.g., planet around a star).
#[derive(Clone, Debug, Serialize, Deserialize)]

pub struct OrbitalSystem {
    pub gravitational_constant : f64,
    pub star_mass : f64,
}

impl MechanicalSystem for OrbitalSystem {
    fn spatial_dim(&self) -> usize {

        2
    }

    fn eval_acceleration(
        &self,
        x : &[f64],
        a : &mut [f64],
    ) {

        let (px, py) = (x[0], x[1]);

        let dist_sq = px.powi(2) + py.powi(2);

        let dist_cubed = dist_sq
            .sqrt()
            .powi(3);

        let force_magnitude = -self.gravitational_constant * self.star_mass / dist_cubed;

        a[0] = force_magnitude * px;

        a[1] = force_magnitude * py;
    }
}

/// Solves a 2D orbital mechanics problem using the energy-conserving semi-implicit Euler method.
///
/// This scenario simulates the orbit of a celestial body around a much more massive star,
/// demonstrating the semi-implicit Euler method's ability to preserve energy and produce
/// stable, realistic orbits over long simulation times.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.

pub fn simulate_gravity_semi_implicit_euler_scenario() -> Result<Vec<(f64, Vec<f64>)>, String> {

    let system = OrbitalSystem {
        gravitational_constant : 1.0,
        star_mass : 1000.0,
    };

    let y0 = &[10.0, 0.0, 0.0, 30.0];

    let t_span = (0.0, 2.5);

    let dt = 0.001;

    solve_semi_implicit_euler(
        &system,
        y0,
        t_span,
        dt,
    )
}

use crate::numerical::matrix::Matrix;

/// Defines a linear ODE system: dy/dt = A * y.

pub trait LinearOdeSystem {
    /// The dimension of the system.

    fn dim(&self) -> usize;

    /// Provides the system matrix A.

    fn get_matrix(&self) -> Matrix<f64>;
}

/// Solves a linear ODE system y' = Ay using the backward (implicit) Euler method.
/// This method is A-stable, making it excellent for stiff equations.
/// The update rule is y_{n+1} = (I - dt*A)^-1 * y_n.
///
/// # Returns
/// A `Result` containing the solution path, or an error string if matrix inversion fails.
/// Solves a linear ODE system `y' = Ay` using the backward (implicit) Euler method.
///
/// This method is A-stable, making it excellent for stiff equations where explicit methods
/// would require prohibitively small time steps. The update rule is `y_{n+1} = (I - dt*A)^-1 * y_n`.
///
/// # Arguments
/// * `system` - The linear ODE system to solve, implementing the `LinearOdeSystem` trait.
/// * `y0` - The initial state vector.
/// * `t_span` - The time interval `(t_start, t_end)`.
/// * `dt` - The fixed time step.
///
/// # Returns
/// A `Result` containing the solution path as `Vec<(f64, Vec<f64>)>`, or an error string if matrix inversion fails.

pub fn solve_backward_euler_linear<S : LinearOdeSystem>(
    system : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
) -> Result<Vec<(f64, Vec<f64>)>, String> {

    let (t_start, t_end) = t_span;

    let steps = ((t_end - t_start) / dt).ceil() as usize;

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let dim = system.dim();

    let a = system.get_matrix();

    let identity = Matrix::identity(dim);

    let m = identity - (a * dt);

    let m_inv = m
        .inverse()
        .ok_or("Matrix (I - dt*A) is not invertible.")?;

    for _ in 0 .. steps {

        let y_matrix = Matrix::new(dim, 1, y.clone());

        let y_new_matrix = m_inv.clone() * y_matrix;

        y.clone_from(&y_new_matrix.get_cols()[0]);

        t += dt;

        history.push((t, y.clone()));
    }

    Ok(history)
}

/// A stiff system where one component decays much faster than the other.

pub struct StiffDecaySystem;

impl LinearOdeSystem for StiffDecaySystem {
    fn dim(&self) -> usize {

        2
    }

    fn get_matrix(&self) -> Matrix<f64> {

        Matrix::new(
            2,
            2,
            vec![
                -20.0, 0.0, 0.0, -0.5,
            ],
        )
    }
}

/// Solves a stiff ODE system, demonstrating the stability of the backward Euler method.
///
/// This scenario highlights how the backward Euler method can handle stiff equations
/// with relatively large time steps, whereas a forward method would become unstable.
///
/// # Returns
/// A `Result` containing the solution path as `Vec<(f64, Vec<f64>)>`, or an error string if matrix inversion fails.

pub fn simulate_stiff_decay_scenario() -> Result<Vec<(f64, Vec<f64>)>, String> {

    let system = StiffDecaySystem;

    let y0 = &[1.0, 1.0];

    let t_span = (0.0, 5.0);

    let dt = 0.2;

    solve_backward_euler_linear(
        &system,
        y0,
        t_span,
        dt,
    )
}
