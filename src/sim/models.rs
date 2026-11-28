//! # Physics Simulations
//!
//! This module contains various physics simulations.
//!
//! ## Coverage of the particle/ODE models
//!
//! `ising_statistical` (2-D lattice) and `geodesic_relativity`
//! (Schwarzschild, equatorial) are fixed-geometry models; the general
//! engines live in [`crate::sim::spin_models`] and
//! [`crate::sim::relativity`].

/// FDTD electrodynamics simulation.
pub use super::models_fdtd_electrodynamics as fdtd_electrodynamics;
/// Geodesic relativity simulation.
pub use super::models_geodesic_relativity as geodesic_relativity;
/// GPE superfluidity simulation.
pub use super::models_gpe_superfluidity as gpe_superfluidity;
/// Ising statistical simulation.
pub use super::models_ising_statistical as ising_statistical;
/// Linear elasticity simulation.
pub use super::models_linear_elasticity as linear_elasticity;
/// Navier-Stokes fluid simulation.
pub use super::models_navier_stokes_fluid as navier_stokes_fluid;
/// Schrodinger quantum simulation.
pub use super::models_schrodinger_quantum as schrodinger_quantum;
