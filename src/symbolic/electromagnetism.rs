//! # Electromagnetism
//!
//! This module provides symbolic representations of fundamental concepts and equations
//! in classical electromagnetism, including Maxwell's equations and field potentials.
//! It allows for symbolic manipulation and derivation of relationships between electric
//! and magnetic fields, charge densities, and current densities.

use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::{
    calculus::differentiate,
    vector::{curl, divergence, gradient},
};

/// Represents Maxwell's equations in their differential form.
///
/// This struct holds four fundamental equations of classical electromagnetism,
/// which describe how electric and magnetic fields are generated and altered
/// by each other and by charges and currents.
#[derive(Debug, Clone)]
pub struct MaxwellEquations {
    /// Gauss's Law for Electricity: ∇ · E = ρ / ε₀
    /// Relates the divergence of the electric field to the charge density.
    pub gauss_law_electric: Expr,

    /// Gauss's Law for Magnetism: ∇ · B = 0
    /// States that there are no magnetic monopoles.
    pub gauss_law_magnetic: Expr,

    /// Faraday's Law of Induction: ∇ × E = -∂B/∂t
    /// Describes how a time-varying magnetic field creates a circulating electric field.
    pub faradays_law: Expr,

    /// Ampère-Maxwell Law: ∇ × B = μ₀J + μ₀ε₀(∂E/∂t)
    /// Relates the curl of the magnetic field to the current density and the rate of change
    /// of the electric field.
    pub amperes_law: Expr,
}

impl MaxwellEquations {
    /// Creates a new set of Maxwell's equations from the given fields and sources.
    ///
    /// The equations are constructed symbolically based on the provided electric field `E`,
    /// magnetic field `B`, charge density `ρ`, and current density `J`.
    ///
    /// # Arguments
    /// * `e_field` - The electric field vector `E(x, y, z, t)`.
    /// * `b_field` - The magnetic field vector `B(x, y, z, t)`.
    /// * `rho` - The charge density `ρ(x, y, z, t)`.
    /// * `j_field` - The current density vector `J(x, y, z, t)`.
    ///
    /// Note: `epsilon_0` (permittivity of free space) and `mu_0` (permeability of free space)
    /// are represented as symbolic variables.
    pub fn new(e_field: &Vector, b_field: &Vector, rho: &Expr, j_field: &Vector) -> Self {
        // Gauss's Law for Electricity: ∇ · E - ρ / ε₀ = 0
        let gauss_law_electric = Expr::Sub(
            Box::new(divergence(e_field, ("x", "y", "z"))),
            Box::new(Expr::Div(
                Box::new(rho.clone()),
                Box::new(Expr::Variable("epsilon_0".to_string())),
            )),
        );

        // Gauss's Law for Magnetism: ∇ · B = 0
        let gauss_law_magnetic = divergence(b_field, ("x", "y", "z"));

        // Faraday's Law: ∇ × E + ∂B/∂t = 0
        let faradays_law = Expr::Add(
            Box::new(curl(e_field, ("x", "y", "z")).to_expr()),
            Box::new(differentiate(&b_field.to_expr(), "t")),
        );

        // Ampère-Maxwell Law: ∇ × B - (μ₀J + μ₀ε₀(∂E/∂t)) = 0
        let term1 = Expr::Mul(
            Box::new(Expr::Variable("mu_0".to_string())),
            Box::new(j_field.to_expr()),
        );
        let term2 = Expr::Mul(
            Box::new(Expr::Variable("mu_0".to_string())),
            Box::new(Expr::Mul(
                Box::new(Expr::Variable("epsilon_0".to_string())),
                Box::new(differentiate(&e_field.to_expr(), "t")),
            )),
        );
        let amperes_law = Expr::Sub(
            Box::new(curl(b_field, ("x", "y", "z")).to_expr()),
            Box::new(Expr::Add(Box::new(term1), Box::new(term2))),
        );

        Self {
            gauss_law_electric,
            gauss_law_magnetic,
            faradays_law,
            amperes_law,
        }
    }
}

/// Calculates the electric field `E` from the scalar electric potential `V`.
///
/// The relationship is given by `E = -∇V`, where `∇` is the gradient operator.
/// This is valid for static electric fields (electrostatics).
///
/// # Arguments
/// * `potential` - The scalar potential `V(x, y, z)`.
///
/// # Returns
/// A `Vector` representing the electric field `E`.
pub fn electric_field_from_potential(potential: &Expr) -> Vector {
    gradient(potential, ("x", "y", "z")).scalar_mul(&Expr::Constant(-1.0))
}

/// Calculates the magnetic field `B` from the magnetic vector potential `A`.
///
/// The relationship is given by `B = ∇ × A`, where `∇ ×` is the curl operator.
///
/// # Arguments
/// * `vector_potential` - The vector potential `A(x, y, z)`.
///
/// # Returns
/// A `Vector` representing the magnetic field `B`.
pub fn magnetic_field_from_vector_potential(vector_potential: &Vector) -> Vector {
    curl(vector_potential, ("x", "y", "z"))
}
