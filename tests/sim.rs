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
#[path = "sim/physics_sm.rs"]
mod physics_sm;
