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

/// Common boundary-condition description shared by the continuum solvers.
pub mod boundary;
pub mod classical;
/// Finite elements (P1/P2) on unstructured triangle meshes, adaptive.
pub mod fem_triangles;
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
/// Iterative 3D BEM with ACA-compressed operators and GMRES.
pub mod bem_aca;
/// Linear elasticity on triangle/tetrahedron meshes (P1/P2).
pub mod elasticity_fem;
/// P1/P2 finite elements on interval, triangle and tetrahedron meshes.
pub mod fem_simplex;
/// Incompressible Navier-Stokes projection solver on a MAC grid (2D/3D).
pub mod ns_mac;
/// Yee FDTD for Maxwell in 1D-3D with PML and TF/SF.
pub mod yee_fdtd;
/// Pseudo-spectral Fourier and Chebyshev solvers (1D-3D).
pub mod spectral;
/// Split-step Fourier Schrodinger/GPE solver (1D-3D, imaginary time).
pub mod split_step;
/// Nodal discontinuous Galerkin (1D/2D) for conservation laws.
pub mod dg;
/// High-order finite volumes for Euler and shallow water in 1D-3D.
pub mod euler_nd;
/// Geometric multigrid (V/W/FMG) for Poisson/Helmholtz in 1D/2D/3D.
pub mod multigrid;
pub mod physics_rkm;
/// Spectral methods.
pub mod physics_sm;
pub mod scenario;
/// Dimension-generic structured-grid finite differences (1D/2D/3D).
pub mod structured;
