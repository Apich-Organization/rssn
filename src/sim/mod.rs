//! # Simulation
//!
//! Time stepping, discretised fields and particle systems.
//!
//! This is the part of
//! scientific computing that is not an identity transformation. A
//! simulation consumes closed-form terms (compiled through a
//! [`Backend`](crate::backend::Backend)) and produces data, not terms.
//!
//! ## Coverage of the particle / ODE modules
//!
//! | Module | Dimensions | Boundaries | Highlights |
//! |---|---|---|---|
//! | [`md_engine`] | any (const generic) | open, periodic, reflecting, thermal, Lees-Edwards | cell+Verlet lists, LJ/Morse/Coulomb, bonded, NHC/Langevin/Berendsen, Berendsen/MTK barostat, SHAKE/RATTLE, RDF/MSD/VACF |
//! | [`ewald`] | 3 | periodic | Ewald, smooth PME |
//! | [`nbody`] | 1-3 (FMM: 2) | open | Barnes-Hut, 2-D FMM |
//! | [`symplectic`] | any | - | order 2-8 schemes, TTL, KS, rigid body |
//! | [`spin_models`] | any | periodic / open | Ising, Potts, O(N); Metropolis, heat bath, Wolff, SW, PT, WL |
//! | [`relativity`] | 3+1 | horizon / escape events | Kerr-Newman geodesics |
//! | [`sph`] | any | open, periodic, reflecting, fixed particles | cubic / Wendland kernels |
//! | `classical`, `physics_md`, `physics_rkm`, `integrators` | fixed or runtime | see module docs | legacy / ODE tools |

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
/// Dissipative particle dynamics (2D/3D).
pub mod dpd;
/// Linear acoustics / scalar wave equation with PML (1D-3D).
pub mod acoustics;
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
/// Named JSON scenarios of the continuum and particle engines.
pub mod scenario_engines;
/// Ewald summation and particle-mesh Ewald.
pub mod ewald;
/// Dimension-generic molecular dynamics engine.
pub mod md_engine;
/// Barnes-Hut and fast multipole N-body solvers.
pub mod nbody;
/// Kerr-Newman geodesics.
pub mod relativity;
/// Smoothed-particle hydrodynamics in 1-3 dimensions.
pub mod sph;
/// Lattice spin models (Ising, Potts, XY, Heisenberg) with cluster updates.
pub mod spin_models;
/// Composition symplectic integrators, TTL, KS regularisation, rigid bodies.
pub mod symplectic;
/// Dimension-generic structured-grid finite differences (1D/2D/3D).
pub mod structured;
