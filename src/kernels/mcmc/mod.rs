//! Markov chain Monte Carlo samplers and convergence diagnostics.
//!
//! Every sampler takes an explicit seed and uses the small internal
//! [`Pcg32`] generator, so a run can be repeated exactly. The samplers return
//! a [`Chain`] of post-burn-in, thinned draws together with the target log
//! density at each draw; [`diagnostics`] holds the usual mixing and convergence
//! measures.
//!
//! [`Pcg32`]: crate::kernels::mcmc::Pcg32
//! [`Chain`]: crate::kernels::mcmc::Chain
//! [`diagnostics`]: crate::kernels::mcmc::diagnostics

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
pub use gibbs::GibbsOptions;
pub use gibbs::gibbs;
pub use gibbs::gibbs_with_log_posterior;
pub use hmc::HmcOptions;
pub use hmc::hmc;
pub use mh::metropolis_hastings;
pub use mh::random_walk_metropolis;
pub use rng::Pcg32;
