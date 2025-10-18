//! # Relativity Module
//!
//! This module provides symbolic tools for special and general relativity.
//! It includes functions for calculating the Lorentz factor, performing Lorentz
//! transformations, mass-energy equivalence, and representing Einstein's field
//! equations and the geodesic equation.
use crate::symbolic::core::Expr;
use crate::symbolic::tensor::MetricTensor;
/// Calculates the Lorentz factor, `γ = 1 / sqrt(1 - v^2/c^2)`.
///
/// The Lorentz factor quantifies the relativistic effects (time dilation, length contraction)
/// that occur when an object moves at a significant fraction of the speed of light.
///
/// # Arguments
/// * `velocity` - The velocity `v` of the object.
///
/// # Returns
/// An `Expr` representing the Lorentz factor `γ`.
pub fn lorentz_factor(velocity: Expr) -> Expr {
    let c = Expr::Variable("c".to_string());
    let v_squared = Expr::new_pow(velocity, Expr::Constant(2.0));
    let c_squared = Expr::new_pow(c.clone(), Expr::Constant(2.0));
    let ratio = Expr::new_div(v_squared, c_squared);
    let one_minus_ratio = Expr::new_sub(Expr::Constant(1.0), ratio);
    let sqrt_expr = Expr::new_pow(one_minus_ratio, Expr::Constant(0.5));
    Expr::new_div(Expr::Constant(1.0), sqrt_expr)
}
/// Performs a Lorentz transformation for a single coordinate.
///
/// The Lorentz transformation describes how measurements of space and time by two observers
/// are related, especially when one observer is moving at a constant velocity relative to the other.
/// Formulas: `x' = γ * (x - v*t)` and `t' = γ * (t - v*x/c^2)`.
///
/// # Arguments
/// * `x` - The spatial coordinate.
/// * `t` - The time coordinate.
/// * `velocity` - The relative velocity `v` between the frames.
///
/// # Returns
/// A tuple `(x_prime, t_prime)` representing the transformed coordinates.
pub fn lorentz_transformation(x: Expr, t: Expr, velocity: Expr) -> (Expr, Expr) {
    let gamma = lorentz_factor(velocity.clone());
    let c = Expr::Variable("c".to_string());
    let term1_x = Expr::new_mul(velocity.clone(), t.clone());
    let inner_x = Expr::new_sub(x.clone(), term1_x);
    let x_prime = Expr::new_mul(gamma.clone(), inner_x);
    let term1_t = Expr::new_div(
        Expr::new_mul(velocity, x),
        Expr::new_pow(c, Expr::Constant(2.0)),
    );
    let inner_t = Expr::new_sub(t, term1_t);
    let t_prime = Expr::new_mul(gamma, inner_t);
    (x_prime, t_prime)
}
/// Calculates the mass-energy equivalence, `E = m * c^2`.
///
/// This famous equation from special relativity states that mass and energy are
/// interchangeable and are two forms of the same thing.
///
/// # Arguments
/// * `mass` - The mass `m` of the object.
///
/// # Returns
/// An `Expr` representing the energy `E`.
pub fn mass_energy_equivalence(mass: Expr) -> Expr {
    let c = Expr::Variable("c".to_string());
    let c_squared = Expr::new_pow(c, Expr::Constant(2.0));
    Expr::new_mul(mass, c_squared)
}
/// Represents Einstein's field equations, `G_μν = (8πG/c^4) * T_μν`.
///
/// These equations form the core of general relativity, describing how spacetime
/// is curved by matter and energy. `G_μν` is the Einstein tensor, and `T_μν` is the
/// stress-energy tensor. `G` is Newton's gravitational constant, and `c` is the speed of light.
/// (Note: `c^4` is often absorbed into `G` for simplicity in some contexts).
///
/// # Arguments
/// * `ricci_tensor` - The Ricci tensor `R_μν`.
/// * `scalar_curvature` - The scalar curvature `R`.
/// * `metric_tensor` - The `MetricTensor` `g_μν`.
/// * `stress_energy_tensor` - The stress-energy tensor `T_μν`.
///
/// # Returns
/// An `Expr` representing the symbolic Einstein field equations.
pub fn einstein_field_equations(
    ricci_tensor: Expr,
    scalar_curvature: Expr,
    metric_tensor: &MetricTensor,
    stress_energy_tensor: Expr,
) -> Result<Expr, String> {
    let g_const = Expr::Variable("G".to_string());
    let pi = Expr::Variable("pi".to_string());
    let term1 = ricci_tensor;
    let term2 = Expr::new_mul(
        Expr::Constant(0.5),
        Expr::new_mul(scalar_curvature, metric_tensor.g.to_matrix_expr()?),
    );
    let einstein_tensor = Expr::new_sub(term1, term2);
    let rhs = Expr::new_mul(
        Expr::Constant(8.0),
        Expr::new_mul(pi, Expr::new_mul(g_const, stress_energy_tensor)),
    );
    Ok(Expr::new_sub(einstein_tensor, rhs))
}
/// Represents the geodesic equation.
///
/// The geodesic equation describes the path of a particle in curved spacetime.
/// It is given by: `d²x^μ/dτ² + Γ^μ_αβ * (dx^α/dτ) * (dx^β/dτ) = 0`,
/// where `Γ^μ_αβ` are the Christoffel symbols.
///
/// # Arguments
/// * `christoffel_symbols` - The Christoffel symbols `Γ^μ_αβ`.
/// * `position_vec` - The four-position vector `x^μ`.
/// * `tau` - The proper time parameter `τ`.
///
/// # Returns
/// An `Expr` representing the symbolic geodesic equation.
pub fn geodesic_equation(christoffel_symbols: Expr, position_vec: Expr, tau: &str) -> Expr {
    let d2x_dtau2 = crate::symbolic::calculus::differentiate(
        &crate::symbolic::calculus::differentiate(&position_vec, tau),
        tau,
    );
    let christoffel_term = Expr::new_apply(christoffel_symbols, position_vec);
    Expr::new_add(d2x_dtau2, christoffel_term)
}
