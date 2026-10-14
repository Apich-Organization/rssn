//! Integration tests for `rssn::sim`.

// Numeric tests index several arrays in step; index loops read closest to
// the formulas they check.
#![allow(clippy::needless_range_loop)]

#[path = "sim/classical.rs"]
mod classical;
#[path = "sim/models_fdtd_electrodynamics.rs"]
mod models_fdtd_electrodynamics;
#[path = "sim/models_geodesic_relativity.rs"]
mod models_geodesic_relativity;
#[path = "sim/models_gpe_superfluidity.rs"]
mod models_gpe_superfluidity;
#[path = "sim/models_ising_statistical.rs"]
mod models_ising_statistical;
#[path = "sim/models_linear_elasticity.rs"]
mod models_linear_elasticity;
#[path = "sim/models_navier_stokes_fluid.rs"]
mod models_navier_stokes_fluid;
#[path = "sim/models_schrodinger_quantum.rs"]
mod models_schrodinger_quantum;
#[path = "sim/physics_bem.rs"]
mod physics_bem;
#[path = "sim/physics_cfd.rs"]
mod physics_cfd;
#[path = "sim/physics_cnm.rs"]
mod physics_cnm;
#[path = "sim/physics_em.rs"]
mod physics_em;
#[path = "sim/physics_fdm.rs"]
mod physics_fdm;
#[path = "sim/physics_fea.rs"]
mod physics_fea;
#[path = "sim/physics_fem.rs"]
mod physics_fem;
#[path = "sim/physics_fvm.rs"]
mod physics_fvm;
#[path = "sim/physics_md.rs"]
mod physics_md;
#[path = "sim/physics_mm.rs"]
mod physics_mm;
#[path = "sim/physics_mtm.rs"]
mod physics_mtm;
#[path = "sim/physics_rkm.rs"]
mod physics_rkm;
#[path = "sim/multigrid.rs"]
mod multigrid;
#[path = "sim/euler_nd.rs"]
mod euler_nd;
#[path = "sim/dg.rs"]
mod dg;
#[path = "sim/spectral.rs"]
mod spectral;
#[path = "sim/split_step.rs"]
mod split_step;
#[path = "sim/yee_fdtd.rs"]
mod yee_fdtd;
#[path = "sim/ns_mac.rs"]
mod ns_mac;
#[path = "sim/fem_simplex.rs"]
mod fem_simplex;
#[path = "sim/elasticity_fem.rs"]
mod elasticity_fem;
#[path = "sim/bem_aca.rs"]
mod bem_aca;
#[path = "sim/continuum_ext.rs"]
mod continuum_ext;
#[path = "sim/continuum_ext2.rs"]
mod continuum_ext2;
#[path = "sim/continuum_ext3.rs"]
mod continuum_ext3;
#[path = "sim/scenario_engines.rs"]
mod scenario_engines;
#[path = "sim/structured.rs"]
mod structured;
#[path = "sim/physics_sm.rs"]
mod physics_sm;
#[path = "sim/symplectic.rs"]
mod symplectic;
#[path = "sim/spin_models.rs"]
mod spin_models;
#[path = "sim/relativity.rs"]
mod relativity;
#[path = "sim/ewald.rs"]
mod ewald;
#[path = "sim/md_engine.rs"]
mod md_engine;
#[path = "sim/nbody.rs"]
mod nbody;
#[path = "sim/sph.rs"]
mod sph;
