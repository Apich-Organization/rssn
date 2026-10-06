//! # Simulation
//!
//! Time stepping, discretised fields and particle systems.
//!
//! This is the part of
//! scientific computing that is not an identity transformation. A
//! simulation consumes closed-form terms (compiled through a
//! [`Backend`](crate::backend::Backend)) and produces data, not terms.

// The grid kernels index raw buffers inside parallel loops.
#![allow(unsafe_code)]

pub mod classical;
/// Compressible gas dynamics (Euler equations, HLLC, MUSCL).
pub mod gas_dynamics;
/// Symplectic and stiff (implicit, Rosenbrock) integrators.
pub mod integrators;
pub mod models;
/// Boundary element method.
pub mod physics_bem;
pub mod physics_cfd;
/// Crank–Nicolson time stepping.
pub mod physics_cnm;
/// Euler-type integrators for mechanical systems.
pub mod physics_em;
pub mod physics_fdm;
pub mod physics_fea;
pub mod physics_fem;
pub mod physics_fvm;
pub mod physics_md;
/// Meshless methods (smoothed-particle hydrodynamics).
pub mod physics_mm;
/// Multigrid methods.
pub mod physics_mtm;
pub mod physics_rkm;
/// Spectral methods.
pub mod physics_sm;
pub mod scenario;
