//! Markov chain Monte Carlo samplers and convergence diagnostics.
//!
//! Every sampler takes an explicit seed and uses the small internal
//! [`Pcg32`] generator, so a run can be repeated exactly. The samplers return
//! a [`Chain`] of post-burn-in, thinned draws together with the target log
//! density at each draw; [`diagnostics`] holds the usual mixing and convergence
//! measures.

pub mod chain;
pub mod diagnostics;
pub mod gibbs;
pub mod hmc;
pub mod mh;
pub mod rng;

pub use chain::Chain;
pub use chain::MhOptions;
pub use diagnostics::autocorrelation;
pub use diagnostics::effective_sample_size;
pub use diagnostics::geweke;
pub use diagnostics::rhat;
pub use gibbs::gibbs;
pub use hmc::HmcOptions;
pub use hmc::hmc;
pub use mh::metropolis_hastings;
pub use mh::random_walk_metropolis;
pub use rng::Pcg32;
