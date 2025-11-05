//! # Solid-State Physics
//!
//! This module provides symbolic tools for solid-state physics, including representations
//! of crystal lattices, Bloch's theorem for electron wave functions in periodic potentials,
//! and simplified energy band models.
use crate::symbolic::core::Expr;
/// Represents a crystal lattice with basis vectors.
#[derive(Clone, Debug)]
pub struct CrystalLattice {
    pub basis_vectors: Vec<Expr>,
}
impl CrystalLattice {
    /// Creates a new crystal lattice.
    ///
    /// # Arguments
    /// * `basis_vectors` - A vector of `Expr` representing the basis vectors of the lattice.
    ///
    /// # Returns
    /// A new `CrystalLattice` instance.
    pub const fn new(basis_vectors: Vec<Expr>) -> Self {
        Self { basis_vectors }
    }
}
/// Represents Bloch's Theorem: `ψ_k(r) = exp(i*k*r) * u_k(r)`.
///
/// Bloch's theorem states that the wave function of an electron in a periodic potential
/// (like in a crystal lattice) can be expressed as a product of a plane wave `exp(i*k*r)`
/// and a periodic function `u_k(r)` that has the same periodicity as the lattice.
///
/// # Arguments
/// * `k_vector` - The wave vector `k`.
/// * `r_vector` - The position vector `r`.
/// * `periodic_function` - The periodic function `u_k(r)`.
///
/// # Returns
/// An `Expr` representing the Bloch wave function `ψ_k(r)`.
pub fn bloch_theorem(k_vector: Expr, r_vector: Expr, periodic_function: Expr) -> Expr {
    let i = Expr::new_complex(Expr::Constant(0.0), Expr::Constant(1.0));
    let ikr = Expr::new_mul(i, Expr::new_mul(k_vector, r_vector));
    let exp_term = Expr::new_exp(ikr);
    Expr::new_mul(exp_term, periodic_function)
}
/// Represents a simple energy band model, typically using the parabolic band approximation.
///
/// The energy band `E(k)` describes the allowed energy levels for electrons in a crystal
/// as a function of their wave vector `k`. A common approximation is `E(k) = E_c + (hbar^2 * k^2) / (2 * m*)`,
/// where `E_c` is the conduction band minimum, `hbar` is the reduced Planck constant,
/// and `m*` is the effective mass.
///
/// # Arguments
/// * `k_vector` - The wave vector `k`.
/// * `effective_mass` - The effective mass `m*`.
/// * `band_gap` - The band gap energy `E_c` (or other reference energy).
///
/// # Returns
/// An `Expr` representing the symbolic energy band.
pub fn energy_band(k_vector: Expr, effective_mass: Expr, band_gap: Expr) -> Expr {
    let hbar = Expr::Variable("hbar".to_string());
    let hbar_sq = Expr::new_pow(hbar, Expr::Constant(2.0));
    let k_sq = Expr::new_pow(k_vector, Expr::Constant(2.0));
    let numerator = Expr::new_mul(hbar_sq, k_sq);
    let denominator = Expr::new_mul(Expr::Constant(2.0), effective_mass);
    let kinetic_term = Expr::new_div(numerator, denominator);
    Expr::new_add(band_gap, kinetic_term)
}
