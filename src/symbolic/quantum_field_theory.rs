//! # Quantum Field Theory
//!
//! This module provides symbolic representations of fundamental concepts and equations
//! in quantum field theory (QFT). It includes symbolic Lagrangians for QED, QCD,
//! and scalar fields, as well as propagators and scattering formulas.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use std::sync::Arc;

/// Computes the Dirac adjoint of a fermion field: `ψ̄ = ψ†γ⁰`.
pub fn dirac_adjoint(psi: &Expr) -> Expr {
    let gamma_0 = Expr::new_variable("gamma_0");
    simplify(&Expr::new_mul(psi.clone(), gamma_0))
}

/// Computes the Feynman slash notation: `A̸ = γμ Aμ`.
pub fn feynman_slash(v_mu: &Expr) -> Expr {
    let gamma_mu = Expr::new_variable("gamma_mu");
    simplify(&Expr::new_mul(gamma_mu, v_mu.clone()))
}

/// Lagrangian density for a free real scalar field (Klein-Gordon):
/// `L = 1/2 (∂μϕ ∂μϕ - m²ϕ²)`.
pub fn scalar_field_lagrangian(phi: &Expr, m: &Expr) -> Expr {
    let half = Expr::Constant(0.5);
    let d_mu_phi = Expr::new_variable("partial_mu_phi");
    let d_mu_phi_sq = Expr::new_pow(d_mu_phi, Expr::Constant(2.0));
    let m_sq = Expr::new_pow(m.clone(), Expr::Constant(2.0));
    let phi_sq = Expr::new_pow(phi.clone(), Expr::Constant(2.0));
    let mass_term = Expr::new_mul(m_sq, phi_sq);
    let diff = Expr::new_sub(d_mu_phi_sq, mass_term);
    simplify(&Expr::new_mul(half, diff))
}

/// Lagrangian density for Quantum Electrodynamics (QED):
/// `L = ψ̄(iD̸ - m)ψ - 1/4 Fμν Fμν`
/// where `Dμ = ∂μ + ieAμ` and `Fμν = ∂μAν - ∂νAμ`.
pub fn qed_lagrangian(psi_bar: &Expr, psi: &Expr, a_mu: &Expr, m: &Expr, e: &Expr) -> Expr {
    let i = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
    let partial_slash = feynman_slash(&Expr::new_variable("partial_mu"));
    let a_slash = feynman_slash(a_mu);
    
    // iD_slash = i*gamma_mu*(partial_mu + i*e*A_mu) = i*partial_slash - e*A_slash
    let id_slash = Expr::new_sub(
        Expr::new_mul(i, partial_slash),
        Expr::new_mul(e.clone(), a_slash)
    );
    
    let dirac_part = Expr::new_mul(psi_bar.clone(), Expr::new_mul(Expr::new_sub(id_slash, m.clone()), psi.clone()));
    
    let f_mu_nu = Expr::new_variable("F_mu_nu");
    let f_sq = Expr::new_pow(f_mu_nu, Expr::Constant(2.0));
    let gauge_part = Expr::new_mul(Expr::Constant(-0.25), f_sq);
    
    simplify(&Expr::new_add(dirac_part, gauge_part))
}

/// Lagrangian density for Quantum Chromodynamics (QCD):
/// `L = Σ ψ̄_i (iD̸ - m)_ij ψ_j - 1/4 G^a_μν G^a_μν`.
pub fn qcd_lagrangian(psi_bar: &Expr, psi: &Expr, g_mu: &Expr, m: &Expr, gs: &Expr) -> Expr {
    let i = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
    let partial_slash = feynman_slash(&Expr::new_variable("partial_mu"));
    let g_slash = feynman_slash(g_mu);
    
    // D_mu = ∂_mu - i*gs*A_mu^a T^a
    // iD_slash = i*partial_slash + gs*A_slash
    let id_slash = Expr::new_add(
        Expr::new_mul(i, partial_slash),
        Expr::new_mul(gs.clone(), g_slash)
    );
    
    let quark_part = Expr::new_mul(psi_bar.clone(), Expr::new_mul(Expr::new_sub(id_slash, m.clone()), psi.clone()));
    
    let g_strength = Expr::new_variable("G_mu_nu_a");
    let g_sq = Expr::new_pow(g_strength, Expr::Constant(2.0));
    let gluon_part = Expr::new_mul(Expr::Constant(-0.25), g_sq);
    
    simplify(&Expr::new_add(quark_part, gluon_part))
}

/// Represents a propagator for a particle in QFT.
///
/// For a scalar: `i / (p² - m² + iε)`
/// For a fermion: `i(p̸ + m) / (p² - m² + iε)`
pub fn propagator(p: &Expr, m: &Expr, is_fermion: bool) -> Expr {
    let i = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
    let p_sq = Expr::new_pow(p.clone(), Expr::Constant(2.0));
    let m_sq = Expr::new_pow(m.clone(), Expr::Constant(2.0));
    let eps = Expr::new_mul(i.clone(), Expr::new_variable("epsilon"));
    let denominator = Expr::new_add(Expr::new_sub(p_sq, m_sq), eps);
    
    if is_fermion {
        let p_slash = feynman_slash(p);
        let numerator = Expr::new_mul(i, Expr::new_add(p_slash, m.clone()));
        simplify(&Expr::new_div(numerator, denominator))
    } else {
        simplify(&Expr::new_div(i, denominator))
    }
}

/// Scattering cross-section: `dσ ∝ |M|² / (flux) * dΦ`.
pub fn scattering_cross_section(matrix_element: &Expr, flux: &Expr, phase_space: &Expr) -> Expr {
    let m_sq = Expr::new_pow(Expr::new_abs(matrix_element.clone()), Expr::Constant(2.0));
    simplify(&Expr::new_mul(Expr::new_div(m_sq, flux.clone()), phase_space.clone()))
}

/// Feynman propagator in position space (symbolic integral representation).
pub fn feynman_propagator_position_space(x: &Expr, y: &Expr, m: &Expr) -> Expr {
    let p = Expr::new_variable("p");
    let prop_p = propagator(&p, m, false);
    let diff = Expr::new_sub(x.clone(), y.clone());
    let exponent = Expr::new_mul(Expr::new_complex(Expr::Constant(0.0), Expr::Constant(-1.0)), Expr::new_mul(p.clone(), diff));
    let integrand = Expr::new_mul(prop_p, Expr::new_exp(exponent));
    
    simplify(&Expr::Integral {
        integrand: Arc::new(integrand),
        var: Arc::new(p),
        lower_bound: Arc::new(Expr::NegativeInfinity),
        upper_bound: Arc::new(Expr::Infinity),
    })
}

