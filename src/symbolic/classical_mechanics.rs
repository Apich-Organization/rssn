//! # Classical Mechanics Module
//!
//! This module provides functions and structures related to classical mechanics,
//! covering concepts from Newtonian, Lagrangian, and Hamiltonian mechanics.
//! It includes tools for kinematics, dynamics, energy, momentum, and rotational motion.
//!
//! ## Key Formulas
//! - **Newton's Second Law**: $F = m a$
//! - **Kinetic Energy**: $T = \frac{1}{2} m v^2$
//! - **Lagrangian**: $L = T - V$
//! - **Hamiltonian**: $H = T + V$

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::vector::Vector;
use crate::symbolic::vector_calculus::{line_integral_vector, ParametricCurve};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// # Kinematics State
///
/// Represents the kinematic state of a particle, including its position,
/// velocity, and acceleration vectors.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// second time derivatives of the position with respect to time `t`.
    pub fn new(position: Expr, t_var: &str) -> Self {
        let velocity = differentiate(&position, t_var);
        let acceleration = differentiate(&velocity, t_var);
        Self {
            position,
            velocity,
            acceleration,
        }
    }
}

/// Calculates the force `F` using Newton's second law, `F = m * a`.
///
/// ## Example
/// ```rust
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::classical_mechanics::newtons_second_law;
///
/// let m = Expr::new_variable("m");
/// let a = Expr::new_variable("a");
/// let f = newtons_second_law(&m, &a); // Result: m * a
/// ```
pub fn newtons_second_law(mass: &Expr, acceleration: &Expr) -> Expr {
    simplify(&Expr::new_mul(mass.clone(), acceleration.clone()))
}

/// Calculates the momentum `p` of an object, `p = m * v`.
pub fn momentum(mass: &Expr, velocity: &Expr) -> Expr {
    simplify(&Expr::new_mul(mass.clone(), velocity.clone()))
}

/// Calculates the kinetic energy `T` of an object, `T = 1/2 * m * v^2`.
pub fn kinetic_energy(mass: &Expr, velocity: &Expr) -> Expr {
    simplify(&Expr::new_mul(
        Expr::Constant(0.5),
        Expr::new_mul(mass.clone(), Expr::new_pow(velocity.clone(), Expr::Constant(2.0))),
    ))
}

/// Calculates the gravitational potential energy near Earth's surface, `V = m * g * h`.
pub fn potential_energy_gravity_uniform(mass: &Expr, height: &Expr, g: &Expr) -> Expr {
    simplify(&Expr::new_mul(
        mass.clone(),
        Expr::new_mul(g.clone(), height.clone()),
    ))
}

/// Calculates the universal gravitational potential energy, `V = -G * m1 * m2 / r`.
pub fn potential_energy_gravity_universal(m1: &Expr, m2: &Expr, r: &Expr, g_constant: &Expr) -> Expr {
    simplify(&Expr::new_neg(Arc::new(Expr::new_div(
        Expr::new_mul(g_constant.clone(), Expr::new_mul(m1.clone(), m2.clone())),
        r.clone(),
    ))))
}

/// Calculates the potential energy of a spring, `V = 1/2 * k * x^2`.
pub fn potential_energy_spring(k: &Expr, x: &Expr) -> Expr {
    simplify(&Expr::new_mul(
        Expr::Constant(0.5),
        Expr::new_mul(k.clone(), Expr::new_pow(x.clone(), Expr::Constant(2.0))),
    ))
}

/// Calculates the mechanical work done by a constant force, `W = F · d`.
pub fn work_constant_force(force: &Vector, displacement: &Vector) -> Expr {
    force.dot(displacement)
}

/// Calculates the work done by a variable force field along a path.
pub fn work_line_integral(force_field: &Vector, path: &ParametricCurve) -> Expr {
    line_integral_vector(force_field, path)
}

/// Calculates the power delivered by a force, `P = F · v`.
pub fn power(force: &Vector, velocity: &Vector) -> Expr {
    force.dot(velocity)
}

/// Calculates the torque, `τ = r × F`.
pub fn torque(r: &Vector, force: &Vector) -> Vector {
    r.cross(force)
}

/// Calculates the angular momentum, `L = r × p`.
pub fn angular_momentum(r: &Vector, p: &Vector) -> Vector {
    r.cross(p)
}

/// Calculates the centripetal acceleration, `a_c = v^2 / r`.
pub fn centripetal_acceleration(velocity: &Expr, radius: &Expr) -> Expr {
    simplify(&Expr::new_div(
        Expr::new_pow(velocity.clone(), Expr::Constant(2.0)),
        radius.clone(),
    ))
}

/// Calculates the moment of inertia for a point mass, `I = m * r^2`.
pub fn moment_of_inertia_point_mass(mass: &Expr, radius: &Expr) -> Expr {
    simplify(&Expr::new_mul(
        mass.clone(),
        Expr::new_pow(radius.clone(), Expr::Constant(2.0)),
    ))
}

/// Calculates the rotational kinetic energy, `T_rot = 1/2 * I * ω^2`.
pub fn rotational_kinetic_energy(moment_of_inertia: &Expr, angular_velocity: &Expr) -> Expr {
    simplify(&Expr::new_mul(
        Expr::Constant(0.5),
        Expr::new_mul(
            moment_of_inertia.clone(),
            Expr::new_pow(angular_velocity.clone(), Expr::Constant(2.0)),
        ),
    ))
}

/// Calculates the Lagrangian `L = T - V`.
pub fn lagrangian(kinetic_energy: &Expr, potential_energy: &Expr) -> Expr {
    simplify(&Expr::new_sub(kinetic_energy.clone(), potential_energy.clone()))
}

/// Calculates the Hamiltonian `H = T + V`.
pub fn hamiltonian(kinetic_energy: &Expr, potential_energy: &Expr) -> Expr {
    simplify(&Expr::new_add(kinetic_energy.clone(), potential_energy.clone()))
}

/// Computes the left-hand side of the Euler-Lagrange equation.
/// Formula: `d/dt (∂L/∂(q_dot)) - dL/dq = 0`.
pub fn euler_lagrange_equation(lagrangian: &Expr, q: &str, q_dot: &str, t_var: &str) -> Expr {
    // 1. Partial derivative wrt q
    let dl_dq = differentiate(lagrangian, q);
    
    // 2. Partial derivative wrt q_dot
    let dl_dq_dot = differentiate(lagrangian, q_dot);
    
    // 3. Before taking d/dt, we must ensure q and q_dot are substituted with 
    // their symbolic time-dependent forms if we want d/dt to be non-zero.
    let q_time = Expr::Variable(q.to_string()); // In a more advanced version, this could be q(t)
    let q_dot_time = Expr::Derivative(Arc::new(q_time), t_var.to_string());
    
    let dl_dq_dot_time = crate::symbolic::calculus::substitute_expr(
        &dl_dq_dot,
        &Expr::Variable(q_dot.to_string()),
        &q_dot_time
    );
    
    let d_dt_dl_dq_dot = differentiate(&dl_dq_dot_time, t_var);
    
    simplify(&Expr::new_sub(d_dt_dl_dq_dot, dl_dq))
}

/// Calculates the Poisson bracket `{f, g} = (∂f/∂q)(∂g/∂p) - (∂f/∂p)(∂g/∂q)`.
pub fn poisson_bracket(f: &Expr, g: &Expr, q: &str, p: &str) -> Expr {
    let df_dq = differentiate(f, q);
    let dg_dp = differentiate(g, p);
    let df_dp = differentiate(f, p);
    let dg_dq = differentiate(g, q);
    simplify(&Expr::new_sub(
        Expr::new_mul(df_dq, dg_dp),
        Expr::new_mul(df_dp, dg_dq),
    ))
}
