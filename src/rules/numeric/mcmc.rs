//! Markov chain Monte Carlo diagnostics.
//!
//! | operator | value |
//! |---|---|
//! | `rhat(chains)` | the split Gelman–Rubin potential scale reduction of several chains of the same parameter |
//! | `effective_sample_size(chain)` | the effective sample size of one chain, by the initial positive sequence of its autocorrelations |
//! | `chain_autocorrelation(chain, max_lag)` | the sample autocorrelation of a chain |
//! | `geweke(chain)` | the Geweke z-score of the first tenth against the last half |
//!
//! The samplers themselves — random-walk Metropolis–Hastings, Hamiltonian Monte
//! Carlo with leapfrog integration, and a component-wise Gibbs sweep — are
//! `crate::kernels::mcmc`. They take a Rust closure for the log density, which
//! a list of numbers cannot carry, so they are a library API rather than
//! operators, exactly as the sparse SQP subproblem is. What a term *can*
//! carry is the draws, and the diagnostics are what a host checks them with:
//! `rhat` near one, an effective size near the draw count, a Geweke score
//! inside two.
//!
//! `chains` is a list of chains (each a list of numbers) for `rhat`, and a
//! single list of numbers for the others.

use super::bridge::Data;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::require;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;

/// Every operator of this family.
reviewed_entries! {
    (
        "rhat", "t", rhat, mcmc::rhat,
        "`rhat(chains)`: the split Gelman–Rubin potential scale reduction of a list of chains of the same parameter; values near 1 indicate the chains have mixed."
    ),
    (
        "effective_sample_size", "v", effective_sample_size, mcmc::effective_sample_size,
        "`effective_sample_size(chain)`: the effective number of independent draws in `chain`, from the initial positive sequence of its autocorrelations."
    ),
    (
        "chain_autocorrelation", "vd", chain_autocorrelation, mcmc::autocorrelation,
        "`chain_autocorrelation(chain, max_lag)`: the sample autocorrelation of `chain` for lags `0..=max_lag`."
    ),
    (
        "geweke", "v", geweke, mcmc::geweke,
        "`geweke(chain)`: the Geweke convergence z-score of the first tenth against the last half of `chain`."
    ),
}

fn rhat(args: &[Data]) -> Option<Data> {
    let Some(Data::Tree(chains)) = args.first() else {
        return bad(
            0,
            DiagKind::Shape,
            format_args!("expected a list of chains"),
        );
    };
    require(
        chains.len() >= 2,
        Some(0),
        DiagKind::Domain,
        format_args!("at least two chains are needed, got {}", chains.len()),
    )?;
    let mut samples = Vec::with_capacity(chains.len());
    for (index, chain) in chains.iter().enumerate() {
        let Some(values) = numbers_of(chain) else {
            return bad(
                0,
                DiagKind::Shape,
                format_args!("chain {index} must be a list of numbers"),
            );
        };
        samples.push(values);
    }
    match crate::kernels::mcmc::rhat(&samples) {
        | Ok(value) => super::bridge::num(value),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

/// The numbers of a list literal.
fn numbers_of(data: &Data) -> Option<Vec<f64>> {
    match data {
        | Data::Array { shape, data } if shape.len() == 1 => Some(data.clone()),
        | _ => None,
    }
}

fn effective_sample_size(args: &[Data]) -> Option<Data> {
    let chain = vector(args, 0)?;
    match crate::kernels::mcmc::effective_sample_size(chain) {
        | Ok(value) => super::bridge::num(value),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn chain_autocorrelation(args: &[Data]) -> Option<Data> {
    let chain = vector(args, 0)?;
    let lags = count(args, 1, 4096)?;
    super::bridge::vec_out(crate::kernels::mcmc::autocorrelation(chain, lags))
}

fn geweke(args: &[Data]) -> Option<Data> {
    let chain = vector(args, 0)?;
    match crate::kernels::mcmc::geweke(chain) {
        | Ok(value) => super::bridge::num(value),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
