//! Integration tests for `rssn::sim`.

#[path = "sim/classical.rs"]
mod classical;
#[path = "sim/physics_cfd.rs"]
mod physics_cfd;
#[path = "sim/physics_fea.rs"]
mod physics_fea;
#[path = "sim/physics_md.rs"]
mod physics_md;
#[path = "sim/physics_bem.rs"]
mod physics_bem;
#[path = "sim/physics_cnm.rs"]
mod physics_cnm;
#[path = "sim/physics_em.rs"]
mod physics_em;
#[path = "sim/physics_fdm.rs"]
mod physics_fdm;
