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
pub mod fdtd_electrodynamics;
/// Geodesic relativity simulation.
pub mod geodesic_relativity;
/// GPE superfluidity simulation.
pub mod gpe_superfluidity;
/// Ising statistical simulation.
pub mod ising_statistical;
/// Linear elasticity simulation.
pub mod linear_elasticity;
/// Navier-Stokes fluid simulation.
pub mod navier_stokes_fluid;
/// Schrodinger quantum simulation.
pub mod schrodinger_quantum;
