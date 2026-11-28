//! Bounded seed-only MCMC specializations with typed host callbacks.

use super::codecs::Arg;
use super::codecs::Output;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeOutput;
use super::native_object_api::error;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::TensorValue;
use crate::kernels::mcmc::GibbsOptions;
use crate::kernels::mcmc::HmcOptions;
use crate::kernels::mcmc::MhOptions;
use crate::kernels::mcmc::{
    self,
};
use std::cell::RefCell;

const MAX_COORDINATES: usize = 32;
const MAX_DRAWS: usize = 8_192;
const MAX_ITERATIONS: usize = 65_536;
const MAX_WORK: usize = 1_048_576;

struct Config {
    draws: usize,
    burn_in: usize,
    thin: usize,
    seed: u64,
    iterations: usize,
}

impl NativeArg for Config {
    fn ty() -> TypeExpr {
        TypeExpr::Record(
            ["draws", "burn_in", "thin", "seed"]
                .into_iter()
                .map(|key| (key.into(), TypeExpr::Int))
                .collect(),
        )
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let draws = <usize as Arg>::decode(api::field(value, "draws")?, index)?;
        let burn_in = <usize as Arg>::decode(api::field(value, "burn_in")?, index)?;
        let thin = <usize as Arg>::decode(api::field(value, "thin")?, index)?;
        let seed = <u64 as Arg>::decode(api::field(value, "seed")?, index)?;
        let iterations = draws.saturating_mul(thin).saturating_add(burn_in);
        if !(1..=MAX_DRAWS).contains(&draws) || thin == 0 || iterations > MAX_ITERATIONS {
            return Err(OpError::BudgetExceeded);
        }
        Ok(Self {
            draws,
            burn_in,
            thin,
            seed,
            iterations,
        })
    }
}

impl Config {
    fn work(
        &self,
        init: &[f64],
        effort: usize,
    ) -> Result<(), OpError> {
        if !(1..=MAX_COORDINATES).contains(&init.len()) {
            return Err(error("MCMC supports one through 32 normalized coordinates"));
        }
        if self
            .iterations
            .saturating_mul(init.len())
            .saturating_mul(effort)
            > MAX_WORK
            || self.draws.saturating_mul(init.len()) > api::VALUE_LIMIT.saturating_sub(MAX_DRAWS)
        {
            return Err(OpError::BudgetExceeded);
        }
        Ok(())
    }

    const fn mh(
        &self,
        proposal_sd: Vec<f64>,
    ) -> MhOptions {
        MhOptions {
            proposal_sd,
            burn_in: self.burn_in,
            thin: self.thin,
            seed: self.seed,
        }
    }
}

const LOG: u8 = 0;
const GRADIENT: u8 = 1;
const CONDITIONAL: u8 = 2;

struct Callback<const ROLE: u8>(CallbackRef);

impl<const ROLE: u8> NativeArg for Callback<ROLE> {
    fn ty() -> TypeExpr {
        let state = api::tensor_type(1, TypeExpr::Float);
        let params = if ROLE == CONDITIONAL {
            vec![TypeExpr::Int, state.clone(), TypeExpr::Float]
        } else {
            vec![state.clone()]
        };
        let result = if ROLE == GRADIENT {
            state
        } else {
            TypeExpr::Float
        };
        TypeExpr::Function(params, Box::new(result))
    }

    fn decode(
        value: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Function(handle) = value else {
            return Err(error("expected a typed MCMC callback"));
        };
        let signature = ctx.callback_signature(*handle)?;
        let TypeExpr::Function(params, result) = Self::ty() else {
            unreachable!("MCMC callback uses a function type")
        };
        if signature.params != params || signature.ret != *result {
            return Err(error(
                "MCMC callback has an incompatible normalized signature",
            ));
        }
        Ok(Self(*handle))
    }
}

struct Eval<'a, 'b> {
    ctx: &'a CallCtx<'b>,
    first_failure: RefCell<Option<OpError>>,
}
impl<'a, 'b> Eval<'a, 'b> {
    const fn new(ctx: &'a CallCtx<'b>) -> Self {
        Self {
            ctx,
            first_failure: RefCell::new(None),
        }
    }

    fn run<T>(
        &self,
        evaluate: impl FnOnce() -> Result<T, OpError>,
    ) -> Option<T> {
        if self.first_failure.borrow().is_some() {
            return None;
        }
        match crate::cancellation::checkpoint()
            .map_err(OpError::from)
            .and_then(|()| evaluate())
        {
            | Ok(value) => Some(value),
            | Err(failure) => {
                *self.first_failure.borrow_mut() = Some(failure);
                None
            },
        }
    }

    fn scalar(
        &self,
        handle: CallbackRef,
        args: &[Value],
    ) -> Result<f64, OpError> {
        let output = self.ctx.invoke_callback(handle, args)?;
        api::same_dimension(&api::dimension(&output)?, &DimExpr::Dimless)?;
        <f64 as Arg>::decode(&output, 0)
    }

    fn log(
        &self,
        function: &Callback<LOG>,
        state: &[f64],
    ) -> f64 {
        self.run(|| {
            let output = self
                .ctx
                .invoke_callback(function.0, &[state_value(state)?])?;
            api::same_dimension(&api::dimension(&output)?, &DimExpr::Dimless)?;
            if output
                .as_f64()
                .is_some_and(|value| value.is_infinite() && value.is_sign_negative())
            {
                return Ok(f64::NEG_INFINITY);
            }
            <f64 as Arg>::decode(&output, 0)
        })
        .unwrap_or(f64::NAN)
    }

    fn gradient(
        &self,
        function: &Callback<GRADIENT>,
        state: &[f64],
    ) -> Vec<f64> {
        self.run(|| {
            let output = self
                .ctx
                .invoke_callback(function.0, &[state_value(state)?])?;
            api::same_dimension(&api::dimension(&output)?, &DimExpr::Dimless)?;
            <Vec<f64> as Arg>::decode(&output, 0)
        })
        .unwrap_or_default()
    }

    fn conditional(
        &self,
        function: &Callback<CONDITIONAL>,
        coordinate: usize,
        state: &[f64],
        u: f64,
    ) -> f64 {
        self.run(|| {
            self.scalar(
                function.0,
                &[Value::int(coordinate), state_value(state)?, Value::Float(u)],
            )
        })
        .unwrap_or(f64::NAN)
    }

    fn finish(
        &self,
        chain: Result<mcmc::Chain, String>,
        keep_log: bool,
    ) -> Result<ChainResult, OpError> {
        if let Some(failure) = self.first_failure.borrow_mut().take() {
            return Err(failure);
        }
        chain
            .map(|chain| ChainResult { chain, keep_log })
            .map_err(error)
    }
}

fn state_value(state: &[f64]) -> Result<Value, OpError> {
    <Vec<f64> as Output>::encode(state.to_vec())
}

struct ChainResult {
    chain: mcmc::Chain,
    keep_log: bool,
}
impl NativeOutput for ChainResult {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("draws".into(), api::tensor_type(2, TypeExpr::Float)),
            ("log_posterior".into(), api::tensor_type(1, TypeExpr::Float)),
            ("accepted".into(), TypeExpr::Int),
            ("proposed".into(), TypeExpr::Int),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let chain = self.chain;
        let dimensions = chain.draws.first().map_or(0, Vec::len);
        let shape = vec![chain.draws.len(), dimensions];
        let data = chain
            .draws
            .into_iter()
            .flatten()
            .map(<f64 as Output>::encode)
            .collect::<Result<_, _>>()?;
        let log = if self.keep_log {
            chain.log_posterior
        } else {
            Vec::new()
        };
        Ok(api::record(vec![
            (
                "draws",
                Value::Tensor(TensorValue {
                    shape,
                    elem: TypeExpr::Float,
                    data,
                }),
            ),
            ("log_posterior", <Vec<f64> as Output>::encode(log)?),
            ("accepted", Value::int(chain.accepted)),
            ("proposed", Value::int(chain.proposed)),
        ]))
    }
}

type Metropolis =
    fn(&dyn Fn(&[f64]) -> f64, &[f64], usize, &MhOptions) -> Result<mcmc::Chain, String>;

fn run_metropolis(
    context: &CallCtx<'_>,
    method: Metropolis,
    log: &Callback<LOG>,
    init: &[f64],
    proposal_sd: Vec<f64>,
    config: &Config,
) -> Result<ChainResult, OpError> {
    config.work(init, 1)?;
    let eval = Eval::new(context);
    let result = method(
        &|state| eval.log(log, state),
        init,
        config.draws,
        &config.mh(proposal_sd),
    );
    eval.finish(result, true)
}

object_binding!(
    metropolis, "num_mcmc_metropolis_hastings", "crate::kernels::mcmc::mh::metropolis_hastings",
    "Seeded Metropolis–Hastings with normalized coordinates, positive Gaussian proposal scales and a typed log-target function.",
    true, (log: Callback<LOG>, initial: Vec<f64>, proposal_sd: Vec<f64>, config: Config) -> ChainResult,
    |ctx, log, initial, proposal_sd, config| run_metropolis(ctx, mcmc::metropolis_hastings, &log, &initial, proposal_sd, &config)
);
object_binding!(
    random_walk, "num_mcmc_random_walk_metropolis", "crate::kernels::mcmc::mh::random_walk_metropolis",
    "Seeded random-walk Metropolis with normalized coordinates, fixed Gaussian proposal scales and a typed log-target function.",
    true, (log: Callback<LOG>, initial: Vec<f64>, proposal_sd: Vec<f64>, config: Config) -> ChainResult,
    |ctx, log, initial, proposal_sd, config| run_metropolis(ctx, mcmc::random_walk_metropolis, &log, &initial, proposal_sd, &config)
);
object_binding!(
    hmc, "num_mcmc_hmc", "crate::kernels::mcmc::hmc::hmc",
    "Seeded Hamiltonian Monte Carlo on normalized coordinates with typed log-target and gradient functions and a fixed leapfrog schedule.",
    true, (log: Callback<LOG>, gradient: Callback<GRADIENT>, initial: Vec<f64>, step_size: f64, leapfrog_steps: usize, config: Config) -> ChainResult,
    |ctx, log, gradient, initial, step_size, leapfrog_steps, config| {
        if !(1..=128).contains(&leapfrog_steps) {
            return Err(error("HMC requires one through 128 leapfrog steps per proposal"));
        }
        config.work(&initial, leapfrog_steps.saturating_add(2))?;
        let eval = Eval::new(ctx);
        let result = mcmc::hmc(
            &|state| eval.log(&log, state), &|state| eval.gradient(&gradient, state), &initial,
            config.draws, &HmcOptions {
                step_size, leapfrog_steps, burn_in: config.burn_in,
                thin: config.thin, seed: config.seed,
            },
        );
        eval.finish(result, true)
    }
);

fn run_gibbs(
    ctx: &CallCtx<'_>,
    conditional: &Callback<CONDITIONAL>,
    log: Option<&Callback<LOG>>,
    initial: &[f64],
    config: &Config,
) -> Result<ChainResult, OpError> {
    config.work(initial, initial.len().saturating_add(1))?;
    let eval = Eval::new(ctx);
    let callbacks = (0..initial.len())
        .map(|index| {
            let eval = &eval;
            move |state: &[f64], u: f64| eval.conditional(conditional, index, state, u)
        })
        .collect::<Vec<_>>();
    let conditional_refs = callbacks
        .iter()
        .map(|f| f as &dyn Fn(&[f64], f64) -> f64)
        .collect::<Vec<_>>();
    let result = log.map_or_else(
        || {
            mcmc::gibbs(
                &conditional_refs,
                initial,
                config.draws,
                &config.mh(Vec::new()),
            )
        },
        |log| {
            mcmc::gibbs_with_log_posterior(
                &conditional_refs,
                |state| eval.log(log, state),
                initial,
                config.draws,
                &GibbsOptions {
                    burn_in: config.burn_in,
                    thin: config.thin,
                    seed: config.seed,
                },
            )
        },
    );
    eval.finish(result, log.is_some())
}

object_binding!(
    gibbs, "num_mcmc_gibbs", "crate::kernels::mcmc::gibbs::gibbs",
    "Seeded Gibbs sweeps with a typed conditional sampler taking a zero-based coordinate, normalized state and a uniform draw. No joint log values are recorded.",
    true, (conditional: Callback<CONDITIONAL>, initial: Vec<f64>, config: Config) -> ChainResult,
    |ctx, conditional, initial, config| run_gibbs(ctx, &conditional, None, &initial, &config)
);
object_binding!(
    gibbs_with_log, "num_mcmc_gibbs_with_log_posterior", "crate::kernels::mcmc::gibbs::gibbs_with_log_posterior",
    "Seeded Gibbs sweeps with a typed coordinate conditional and joint log target. The log posterior is evaluated at every retained draw.",
    true, (conditional: Callback<CONDITIONAL>, log: Callback<LOG>, initial: Vec<f64>, config: Config) -> ChainResult,
    |ctx, conditional, log, initial, config| run_gibbs(ctx, &conditional, Some(&log), &initial, &config)
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps the sampler entries beside their checked adapters"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &metropolis::BINDING,
            &random_walk::BINDING,
            &hmc::BINDING,
            &gibbs::BINDING,
            &gibbs_with_log::BINDING,
        ],
        "numerical_stats",
    )
}
