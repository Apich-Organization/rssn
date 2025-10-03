//! # Physics Simulation & Numerical Methods Module
//!
//! This module is dedicated to the implementation of various numerical methods for solving problems
//! in physics and engineering. It serves as a collection of powerful techniques for simulating
//! physical phenomena governed by differential equations.
//!
//! ## Core Simulation Techniques
//!
//! The module is organized by numerical method, with each submodule providing a specific approach
//! to discretization and solving. This allows users to choose the most appropriate tool for their
//! specific problem, whether it involves fluid dynamics, electromagnetism, quantum mechanics, or
//! structural analysis.
//!
//! - **`physics_sim`**: Contains high-level, ready-to-run simulation examples that demonstrate
//!   how to use the numerical methods in this module to model real-world physical systems.
//! - **Finite Difference/Volume/Element Methods**: `physics_fdm`, `physics_fvm`, and `physics_fem`
//!   are foundational for solving PDEs by discretizing a domain into a grid or mesh.
//! - **Time-Evolution Methods**: `physics_cnm`, `physics_em`, and `physics_rkm` provide various
//!   schemes for solving time-dependent problems, offering different trade-offs between
//!   accuracy, stability, and computational cost.
//! - **Advanced & Specialized Methods**: Includes `physics_bem` (for boundary-value problems),
//!   `physics_sm` (for high-accuracy solutions on simple geometries), `physics_mm` (for mesh-free
//!   simulations), and `physics_mtm` (for accelerating the convergence of large linear systems).

/// The Boundary Element Method (BEM) for solving potential problems.
pub mod physics_bem;
/// The Crank-Nicolson method for solving time-dependent partial differential equations.
pub mod physics_cnm;
/// The Euler methods (Forward, Backward, Semi-Implicit) for solving ordinary differential equations.
pub mod physics_em;
/// The Finite Difference Method (FDM) for solving differential equations on a grid.
pub mod physics_fdm;
/// The Finite Element Method (FEM) for solving PDEs, especially on complex domains.
pub mod physics_fem;
/// The Finite Volume Method (FVM), ideal for problems involving conservation laws.
pub mod physics_fvm;
/// Meshfree Methods, such as Smoothed-Particle Hydrodynamics (SPH).
pub mod physics_mm;
/// Multigrid Methods (MTM) for efficiently solving large linear systems.
pub mod physics_mtm;
/// Runge-Kutta methods for solving ordinary differential equations with high accuracy.
pub mod physics_rkm;
/// A collection of high-level physics simulation scenarios.
pub mod physics_sim;
/// Spectral Methods (SM) for solving PDEs using global basis functions like Fourier series.
pub mod physics_sm;
