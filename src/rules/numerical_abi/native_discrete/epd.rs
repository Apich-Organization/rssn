//! Dimensioned sampled Euler–Poisson–Darboux solutions.

use super::Arg;
use super::CallCtx;
use super::DimExpr;
use super::Installer;
use super::NativeArg;
use super::NativeOutput;
use super::OpError;
use super::RuleError;
use super::TypeExpr;
use super::Value;
use super::api;
use super::error;
use super::object_binding;
use super::quantities;
use crate::rules::pde::riemann::epd_goursat_samples;

fn measured_tensor(
    rank: u64,
    domain: char,
) -> TypeExpr {
    api::tensor_type(
        rank,
        TypeExpr::Float.with_dim(DimExpr::Var(domain.to_string())),
    )
}

struct MeasuredSamples<const DOMAIN: char>(api::MeasuredVector<128>);

impl<const DOMAIN: char> NativeArg for MeasuredSamples<DOMAIN> {
    fn ty() -> TypeExpr {
        measured_tensor(1, DOMAIN)
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        <api::MeasuredVector<128> as NativeArg>::decode(value, index, ctx).map(Self)
    }
}

impl<const DOMAIN: char> NativeOutput for MeasuredSamples<DOMAIN> {
    fn ty() -> TypeExpr {
        measured_tensor(1, DOMAIN)
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(quantities(self.0.values, &self.0.dimension)?.0)
    }
}

struct Targets {
    data: Vec<[f64; 2]>,
    dimension: DimExpr,
}

impl NativeArg for Targets {
    fn ty() -> TypeExpr {
        measured_tensor(2, 'X')
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(error("Goursat targets must be a matrix of [r, s] pairs"));
        };
        if tensor.shape.first().is_some_and(|&count| count > 256) {
            return Err(error("Goursat exceeds the 256-target budget"));
        }
        Ok(Self {
            data: <Vec<[f64; 2]> as Arg>::decode(value, index)?,
            dimension: api::dimension(value)?,
        })
    }
}

object_binding!(
    sampled_goursat,
    "num_native_epd_goursat_samples",
    "crate::rules::pde::riemann::epd_goursat_samples",
    "Evaluates Euler–Poisson–Darboux Goursat data from two piecewise-linear sampled characteristics at [r, s] targets. Both axes and target coordinates share a physical dimension. The samples share the output dimension. The first samples define the common corner. Supports at most 128 knots per characteristic and 256 targets.",
    false,
    (order: f64, rho: MeasuredSamples<'X'>, phi: MeasuredSamples<'Y'>, sigma: MeasuredSamples<'X'>, psi: MeasuredSamples<'Y'>, targets: Targets) -> MeasuredSamples<'Y'>,
    |_, order, MeasuredSamples(rho), MeasuredSamples(phi), MeasuredSamples(sigma), MeasuredSamples(psi), targets| {
        api::same_dimension(&rho.dimension, &sigma.dimension)?;
        api::same_dimension(&rho.dimension, &targets.dimension)?;
        api::same_dimension(&phi.dimension, &psi.dimension)?;
        let values = epd_goursat_samples(order, &rho.values, &phi.values, &sigma.values, &psi.values, &targets.data)
            .map_err(|e| error(e.to_string()))?;
        Ok(MeasuredSamples(api::MeasuredVector { values, dimension: phi.dimension }))
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, &[&sampled_goursat::BINDING], "numerical_vector")
}
