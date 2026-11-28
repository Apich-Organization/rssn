//! Deterministic MCMC streams and checked convergence diagnostics.

use super::codecs::Arg;
use super::native_object_api::Handle;
use super::native_object_api::MeasuredVector;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::error;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernels::dense::Mat;
use crate::kernels::mcmc;
use std::mem::size_of;

struct PcgStream(mcmc::Pcg32);

impl NativeType for PcgStream {
    const KIND: &'static str = "rssn/native/mcmc/pcg32/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

object_binding!(
    pcg_new, "num_mcmc_pcg32_new", "crate::kernels::mcmc::rng::Pcg32::new",
    "Creates a deterministic PCG32 stream from a u64 seed and a 63-bit stream selector.",
    true, (seed: u64, stream: u64) -> NewObject<PcgStream>,
    |_, seed, stream| {
        if stream > i64::MAX.unsigned_abs() {
            return Err(error("PCG32 stream selector must fit in 63 bits"));
        }
        Ok(NewObject::normalized(PcgStream(mcmc::Pcg32::new(seed, stream))))
    }
);

object_binding!(
    pcg_seeded, "num_mcmc_pcg32_seeded", "crate::kernels::mcmc::rng::Pcg32::seeded",
    "Creates a deterministic PCG32 stream with the default stream selector.",
    true, (seed: u64) -> NewObject<PcgStream>,
    |_, seed| Ok(NewObject::normalized(PcgStream(mcmc::Pcg32::seeded(seed))))
);

object_binding!(
    pcg_word, "num_mcmc_pcg32_next_u32", "crate::kernels::mcmc::rng::Pcg32::next_u32",
    "Advances the specified PCG32 stream and returns its next unsigned 32-bit output.",
    true, (stream: Handle<PcgStream>) -> u32,
    |ctx, stream| stream.with_mut(ctx, |stream| Ok(stream.0.next_u32()))
);

object_binding!(
    pcg_uniform, "num_mcmc_pcg32_next_f64", "crate::kernels::mcmc::rng::Pcg32::next_f64",
    "Advances the specified PCG32 stream and returns a 53-bit uniform draw in [0, 1).",
    true, (stream: Handle<PcgStream>) -> f64,
    |ctx, stream| stream.with_mut(ctx, |stream| Ok(stream.0.next_f64()))
);

object_binding!(
    pcg_normal, "num_mcmc_pcg32_normal", "crate::kernels::mcmc::rng::Pcg32::normal",
    "Draws a standard normal from the specified PCG32 stream, preserving its cached Box–Muller variate.",
    true, (stream: Handle<PcgStream>) -> f64,
    |ctx, stream| stream.with_mut(ctx, |stream| Ok(stream.0.normal()))
);

const MAX_CHAIN: usize = 16_384;
const MAX_ESS_CHAIN: usize = 4_096;
const MAX_TOTAL_SAMPLES: usize = 131_072;
const MAX_LAG_WORK: usize = 16_777_216;

type Chain = MeasuredVector<MAX_CHAIN>;
type EssChain = MeasuredVector<MAX_ESS_CHAIN>;

struct Chains(Mat);
impl NativeArg for Chains {
    fn ty() -> TypeExpr {
        let measured = TypeExpr::Float.with_dim(DimExpr::Var("Samples".into()));
        api::tensor_type(2, measured)
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(error("split Rhat requires a rank-two tensor of chains"));
        };
        let [rows, draws] = tensor.shape.as_slice() else {
            return Err(error("split Rhat requires a rank-two tensor of chains"));
        };
        if !(2..=128).contains(rows)
            || !(4..=MAX_CHAIN).contains(draws)
            || tensor.data.len() > MAX_TOTAL_SAMPLES
        {
            return Err(error("split Rhat exceeds the chain or draw limits"));
        }
        api::dimension(value)?;
        <Mat as Arg>::decode(value, index).map(Self)
    }
}

fn normalize(samples: &mut [f64]) -> Result<(), OpError> {
    let scale = samples.iter().copied().map(f64::abs).fold(0.0, f64::max);
    if !(scale > 0.0 && scale.is_finite()) {
        return Err(error("chain has no finite, nonzero scale"));
    }
    for value in samples {
        *value /= scale;
    }
    Ok(())
}

object_binding!(
    split_rhat, "num_mcmc_split_rhat", "",
    "Split Rhat for rank-two measured draws, with one equally long chain per row. All chains measure the same quantity. The result is dimensionless.",
    false, (chains: Chains) -> f64,
    |_, mut chains| {
        normalize(&mut chains.0.data)?;
        let rows = chains.0.data.chunks_exact(chains.0.cols).map(<[f64]>::to_vec).collect::<Vec<_>>();
        mcmc::rhat(&rows).map_err(error)
    }
);

object_binding!(
    ess, "num_mcmc_effective_sample_size", "",
    "Initial positive sequence estimate of effective sample size from a measured chain of at most 4096 draws.",
    false, (chain: EssChain) -> f64,
    |_, mut chain| {
        normalize(&mut chain.values)?;
        mcmc::effective_sample_size(&chain.values).map_err(error)
    }
);

object_binding!(
    autocorrelation, "num_mcmc_autocorrelation", "",
    "Dimensionless sample autocorrelations at lags zero through max_lag. The chain must vary and the lag and attempted work must fit the limits.",
    false, (chain: Chain, max_lag: usize) -> Vec<f64>,
    |_, mut chain, max_lag| {
        let n = chain.values.len();
        if n < 2 || max_lag >= n || max_lag.saturating_mul(n) > MAX_LAG_WORK {
            return Err(error("autocorrelation lag or attempted work exceeds the limits"));
        }
        normalize(&mut chain.values)?;
        let result = mcmc::autocorrelation(&chain.values, max_lag);
        if result.iter().any(|value| !value.is_finite()) {
            return Err(error("autocorrelation requires a varying chain"));
        }
        Ok(result)
    }
);

object_binding!(
    geweke, "num_mcmc_geweke", "",
    "Dimensionless Geweke z-score comparing the first tenth against the last half of one measured chain.",
    false, (chain: Chain) -> f64,
    |_, mut chain| {
        normalize(&mut chain.values)?;
        mcmc::geweke(&chain.values).map_err(error)
    }
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps the typed diagnostics and stream entries together"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    const NAMESPACE: &str = "numerical_stats";
    api::install_with_dimension(
        i,
        &[
            &pcg_new::BINDING,
            &pcg_seeded::BINDING,
            &autocorrelation::BINDING,
        ],
        NAMESPACE,
        DimRule::Container,
    )?;
    api::install_with_dimension(
        i,
        &[
            &pcg_word::BINDING,
            &pcg_uniform::BINDING,
            &pcg_normal::BINDING,
            &split_rhat::BINDING,
            &ess::BINDING,
            &geweke::BINDING,
        ],
        NAMESPACE,
        DimRule::Dimensionless,
    )
}
