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
//! | [`md_engine`] | any (const generic) | open, periodic, reflecting, thermal, Lees-Edwards | cell+Verlet lists, LJ/Morse/Coulomb, bonded, NHC/Langevin/Berendsen, Berendsen/MTK and flexible-cell Parrinello-Rahman/MTTK barostats on triclinic cells, SHAKE/RATTLE (also with barostat and Langevin), RDF/MSD/VACF |
//! | [`md_rigid`] | 3 | periodic | rigid molecules: COM + quaternion, symmetric-top / NO_SQUISH rotor splitting, molecular virial |
//! | [`ewald`] | 3 | periodic (orthorhombic and triclinic) | Ewald, smooth PME |
//! | [`pic`] | 1 | periodic | electrostatic `PIC`: `CIC` and four-point cubic shapes, spectral Poisson |
//! | [`nbody`], [`fmm`] | 1-3 (FMM: 2, 3) | open (3-D FMM also periodic) | Barnes-Hut, adaptive 2-D/3-D FMM with softened far field, periodic triclinic FMM |
//! | [`sse`] | any bipartite graph | periodic / open | spin-1/2 XXZ quantum Monte Carlo: stochastic series expansion with directed loops |
//! | [`hermite`] | 3 | open | 4th/6th-order Hermite with block time steps, 1PN EIH terms, KS-regularised binaries |
//! | [`symplectic`] | any | - | order 2-8 schemes, TTL (composed to order 8), KS, rigid body (exact symmetric-top flows) |
//! | [`chain`] | 3 | open | Mikkola-Aarseth chain regularisation, algorithmic regularisation + GBS extrapolation |
//! | [`spin_models`] | any | periodic / open | Ising, Potts, O(N); Metropolis, heat bath (also XY/Heisenberg), Wolff, SW, PT (O(N) and Potts), WL (Ising and continuous spins) |
//! | [`relativity`], [`kerr`], [`elliptic`] | 3+1 | horizon / escape events | Kerr-Newman geodesics; exact bound and numerical unbound Kerr geodesics in Mino time, shadows; third-kind elliptic integrals |
//! | [`sph`] | any | open, periodic, reflecting, fixed particles | cubic / Wendland kernels |
//! | `classical`, `physics_md`, `physics_rkm`, `integrators` | fixed or runtime | see module docs | legacy / ODE tools |

// The grid kernels index raw buffers inside parallel loops.

pub mod acoustics;
pub mod bem_aca;
pub mod bem_helmholtz;
pub mod boundary;
pub mod chain;
pub mod classical;
pub mod cutcell;
pub mod dg;
pub mod dg3d;
pub mod dimensional;
pub mod dpd;
pub mod elasticity_fem;
pub mod elliptic;
pub mod euler_amr;
pub mod euler_nd;
pub mod ewald;
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::trivially_copy_pass_by_ref
)]
pub mod fdm_nd;
pub mod fem_simplex;
pub mod fem_triangles;
pub mod fmm;
pub mod gas_dynamics;
pub mod hermite;
pub mod hyperelastic;
pub mod integrators;
pub mod kerr;
pub mod md_engine;
pub mod md_rigid;
pub mod models;
pub mod models_fdtd_electrodynamics;
pub mod models_geodesic_relativity;
pub mod models_gpe_superfluidity;
pub mod models_ising_statistical;
pub mod models_linear_elasticity;
pub mod models_navier_stokes_fluid;
pub mod models_schrodinger_quantum;
pub mod multigrid;
pub mod nbody;
pub mod ns_mac;
pub mod physics_bem;
pub mod physics_cfd;
pub mod physics_cnm;
pub mod physics_em;
pub mod physics_fdm;
pub mod physics_fea;
pub mod physics_fem;
pub mod physics_fvm;
pub mod physics_md;
/// Meshless methods (smoothed-particle hydrodynamics).
pub mod physics_mm;
pub mod physics_mtm;
pub mod physics_rkm;
pub mod physics_sm;
pub mod pic;
pub mod relativity;
pub mod scenario;
pub mod scenario_engines;
pub mod schema;
pub mod spectral;
pub mod sph;
pub mod spin_models;
pub mod split_step;
/// Stochastic series expansion (directed loops) for spin-1/2 XXZ models.
pub mod sse;
pub mod structured;
pub mod symplectic;
pub mod yee_fdtd;

/// Stochastic integration, jump diffusions, Levy processes, and multilevel Monte Carlo.
pub mod stochastic;

/// Structured FEM, FVM, discontinuous Galerkin, and spectral entry points.
pub mod pde_api;
