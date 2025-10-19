//! # Quantum Field Theory
//!
//! This module provides symbolic representations of fundamental concepts and equations
//! in quantum field theory (QFT). It includes symbolic Lagrangians for QED (Quantum
//! Electrodynamics) and QCD (Quantum Chromodynamics), propagators for various particles,
//! and a high-level representation of scattering cross-section calculations.
use crate::symbolic::core::Expr;
/// Represents the QED (Quantum Electrodynamics) Lagrangian for a fermion field `psi`
/// interacting with a photon field `A_mu`.
///
/// The QED Lagrangian describes the dynamics of electrons and photons and their interactions.
/// It is a high-level symbolic representation.
///
/// # Arguments
/// * `psi_bar` - The adjoint fermion field `ψ̄`.
/// * `psi` - The fermion field `ψ`.
/// * `a_mu` - The photon field `A_μ`.
/// * `mass` - The mass `m` of the fermion.
/// * `_charge` - The electric charge `e` (currently unused in this high-level representation).
///
/// # Returns
/// An `Expr` representing the QED Lagrangian.
pub fn qed_lagrangian(psi_bar: &Expr, psi: &Expr, a_mu: Expr, mass: Expr, _charge: Expr) -> Expr {
    let dirac_term = Expr::new_apply(
        Expr::Variable("DiracTerm".to_string()),
        Expr::Tuple(vec![psi_bar.clone(), psi.clone(), a_mu.clone(), mass]),
    );
    let field_strength_term =
        Expr::new_apply(Expr::Variable("FieldStrengthTerm".to_string()), a_mu);
    Expr::new_sub(dirac_term, field_strength_term)
}
/// Represents the QCD (Quantum Chromodynamics) Lagrangian for a quark field `psi`
/// interacting with a gluon field `A_mu^a`.
///
/// The QCD Lagrangian describes the strong interaction between quarks and gluons.
/// It is a high-level symbolic representation.
///
/// # Arguments
/// * `psi_bar` - The adjoint quark field `ψ̄`.
/// * `psi` - The quark field `ψ`.
/// * `a_mu_a` - The gluon field `A_μ^a`.
/// * `mass` - The mass `m` of the quark.
///
/// # Returns
/// An `Expr` representing the QCD Lagrangian.
pub fn qcd_lagrangian(psi_bar: Expr, psi: Expr, a_mu_a: Expr, mass: Expr) -> Expr {
    let quark_term = Expr::new_apply(
        Expr::Variable("QuarkTerm".to_string()),
        Expr::Tuple(vec![psi_bar, psi, a_mu_a.clone(), mass]),
    );
    let gluon_term = Expr::new_apply(Expr::Variable("GluonFieldStrengthTerm".to_string()), a_mu_a);
    Expr::new_sub(quark_term, gluon_term)
}
/// Represents a propagator for a particle in quantum field theory.
///
/// A propagator describes the amplitude for a particle to travel between two points
/// or to transition between two states. It is typically a function of momentum and mass.
///
/// # Arguments
/// * `momentum` - The momentum `p` of the particle.
/// * `mass` - The mass `m` of the particle.
/// * `is_fermion` - A boolean indicating if the particle is a fermion (true) or a boson (false).
///
/// # Returns
/// An `Expr` representing the symbolic propagator.
pub fn propagator(momentum: Expr, mass: Expr, is_fermion: bool) -> Expr {
    let p_squared = Expr::new_pow(momentum.clone(), Expr::Constant(2.0));
    let m_squared = Expr::new_pow(mass.clone(), Expr::Constant(2.0));
    let denominator = Expr::new_sub(p_squared, m_squared);
    if is_fermion {
        let gamma_dot_p = Expr::new_apply(Expr::Variable("GammaDotP".to_string()), momentum);
        let numerator = Expr::new_add(gamma_dot_p, mass);
        Expr::new_div(numerator, denominator)
    } else {
        let numerator = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
        Expr::new_div(numerator, denominator)
    }
}
/// Symbolic representation of a scattering cross-section calculation.
///
/// The scattering cross-section `σ` is a measure of the probability that two particles
/// will scatter off each other. It is typically proportional to the square of the
/// scattering matrix element `|M|^2`, divided by a flux factor and multiplied by a phase space factor.
///
/// # Arguments
/// * `matrix_element` - The scattering matrix element `M`.
/// * `flux_factor` - The flux factor.
/// * `phase_space_factor` - The phase space factor.
///
/// # Returns
/// An `Expr` representing the symbolic scattering cross-section.
pub fn scattering_cross_section(
    matrix_element: Expr,
    flux_factor: Expr,
    phase_space_factor: Expr,
) -> Expr {
    let m_squared = Expr::new_pow(matrix_element, Expr::Constant(2.0));
    let term1 = Expr::new_div(m_squared, flux_factor);
    Expr::new_mul(term1, phase_space_factor)
}
