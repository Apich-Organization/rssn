//! # Classical Mechanics Module
//!
//! This module provides functions and structures related to classical mechanics,
//! covering concepts from both Newtonian and Lagrangian/Hamiltonian mechanics.
//! It includes tools for kinematics, dynamics (Newton's laws), energy, momentum,
//! and the formulation of equations of motion using variational principles.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;

/// # Kinematics State
///
/// Represents the kinematic state of a particle, including its position,
/// velocity, and acceleration vectors.
#[derive(Debug, Clone)]
pub struct Kinematics {
    /// The position vector of the particle (e.g., `r(t)`).
    pub position: Expr,
    /// The velocity vector, `v = dr/dt`.
    pub velocity: Expr,
    /// The acceleration vector, `a = dv/dt = d²r/dt²`.
    pub acceleration: Expr,
}

impl Kinematics {
    /// Creates a new `Kinematics` state from a given position expression.
    ///
    /// Velocity and acceleration are automatically derived by taking the first and
    /// second time derivatives of the position, respectively. The differentiation
    /// is performed with respect to a variable named "t".
    ///
    /// # Arguments
    /// * `position` - An `Expr` representing the position vector of the particle.
    pub fn new(position: Expr) -> Self {
        // Velocity is the first time derivative of position.
        let velocity = differentiate(&position, "t");
        // Acceleration is the first time derivative of velocity.
        let acceleration = differentiate(&velocity, "t");
        Self {
            position,
            velocity,
            acceleration,
        }
    }
}

/// Calculates the force `F` using Newton's second law, `F = m * a`.
///
/// # Arguments
/// * `mass` - An `Expr` representing the mass `m` of the object.
/// * `acceleration` - An `Expr` representing the acceleration `a` of the object.
///
/// # Returns
/// An `Expr` for the force `F`.
pub fn newtons_second_law(mass: Expr, acceleration: Expr) -> Expr {
    Expr::Mul(Box::new(mass), Box::new(acceleration))
}

/// Calculates the momentum `p` of an object, `p = m * v`.
///
/// # Arguments
/// * `mass` - An `Expr` representing the mass `m` of the object.
/// * `velocity` - An `Expr` representing the velocity `v` of the object.
///
/// # Returns
/// An `Expr` for the momentum `p`.
pub fn momentum(mass: Expr, velocity: Expr) -> Expr {
    Expr::Mul(Box::new(mass), Box::new(velocity))
}

/// Calculates the kinetic energy `T` of an object, `T = 0.5 * m * v^2`.
///
/// # Arguments
/// * `mass` - An `Expr` representing the mass `m` of the object.
/// * `velocity` - An `Expr` representing the velocity `v` of the object.
///
/// # Returns
/// An `Expr` for the kinetic energy `T`.
pub fn kinetic_energy(mass: Expr, velocity: Expr) -> Expr {
    Expr::Mul(
        Box::new(Expr::Constant(0.5)),
        Box::new(Expr::Mul(
            Box::new(mass),
            Box::new(Expr::Power(
                Box::new(velocity),
                Box::new(Expr::Constant(2.0)),
            )),
        )),
    )
}

/// Calculates the Lagrangian `L` of a system, defined as `L = T - V`,
/// where `T` is the kinetic energy and `V` is the potential energy.
///
/// # Arguments
/// * `kinetic_energy` - An `Expr` for the kinetic energy `T`.
/// * `potential_energy` - An `Expr` for the potential energy `V`.
///
/// # Returns
/// An `Expr` for the Lagrangian `L`.
pub fn lagrangian(kinetic_energy: Expr, potential_energy: Expr) -> Expr {
    Expr::Sub(Box::new(kinetic_energy), Box::new(potential_energy))
}

/// Calculates the Hamiltonian `H` of a system, defined as `H = T + V`,
/// where `T` is the kinetic energy and `V` is the potential energy.
/// Note: This is only true for a specific class of systems (scleronomic and holonomic).
/// In general, H is derived from the Lagrangian via a Legendre transform.
///
/// # Arguments
/// * `kinetic_energy` - An `Expr` for the kinetic energy `T`.
/// * `potential_energy` - An `Expr` for the potential energy `V`.
///
/// # Returns
/// An `Expr` for the Hamiltonian `H`.
pub fn hamiltonian(kinetic_energy: Expr, potential_energy: Expr) -> Expr {
    Expr::Add(Box::new(kinetic_energy), Box::new(potential_energy))
}

/// Computes the left-hand side of the Euler-Lagrange equation.
///
/// The Euler-Lagrange equation is a fundamental equation in classical mechanics
/// derived from the principle of least action. It describes the path of a system
/// in terms of its generalized coordinates `q` and `q_dot`.
/// Formula: `d/dt (∂L/∂(q_dot)) - dL/dq = 0`.
///
/// # Arguments
/// * `lagrangian` - The Lagrangian `L` of the system.
/// * `q` - The generalized coordinate `q`.
/// * `q_dot` - The generalized velocity `dq/dt`.
///
/// # Returns
/// An `Expr` representing the left-hand side of the Euler-Lagrange equation.
pub fn euler_lagrange_equation(lagrangian: &Expr, q: &Expr, q_dot: &Expr) -> Expr {
    // Partial derivative of L with respect to q_dot.
    let dl_dq_dot = differentiate(lagrangian, &q_dot.to_string());
    // Total time derivative of the above result.
    let d_dt_dl_dq_dot = differentiate(&dl_dq_dot, "t");
    // Partial derivative of L with respect to q.
    let dl_dq = differentiate(lagrangian, &q.to_string());
    // Combine to form the equation.
    Expr::Sub(Box::new(d_dt_dl_dq_dot), Box::new(dl_dq))
}

/// Calculates the Poisson bracket `{f, g}` of two functions `f(q, p, t)` and `g(q, p, t)`.
///
/// The Poisson bracket is a fundamental concept in Hamiltonian mechanics, describing the
/// time evolution of a function on the phase space.
/// It is defined as: `{f, g} = (∂f/∂q)(∂g/∂p) - (∂f/∂p)(∂g/∂q)`.
///
/// # Arguments
/// * `f` - An `Expr` for the first function.
/// * `g` - An `Expr` for the second function.
/// * `q` - The name of the canonical position coordinate.
/// * `p` - The name of the canonical momentum coordinate.
///
/// # Returns
/// An `Expr` representing the Poisson bracket `{f, g}`.
pub fn poisson_bracket(f: &Expr, g: &Expr, q: &str, p: &str) -> Expr {
    let df_dq = differentiate(f, q);
    let dg_dp = differentiate(g, p);
    let df_dp = differentiate(f, p);
    let dg_dq = differentiate(g, q);

    let term1 = Expr::Mul(Box::new(df_dq), Box::new(dg_dp));
    let term2 = Expr::Mul(Box::new(df_dp), Box::new(dg_dq));

    Expr::Sub(Box::new(term1), Box::new(term2))
}
