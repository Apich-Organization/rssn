//! Session-owned statrs distributions exposed through dimensionless queries.
use std::mem::size_of;

use super::codecs::Output;
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::Value;
use crate::kernels::distributions_ext::StatrsDistribution;

use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::object_binding;

const MAX_PARAMETERS: usize = 4096;

impl NativeType for StatrsDistribution {
    const KIND: &'static str = "rssn/native/statrs/distribution/v1";
    const MAX_BYTES: usize = 64 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if let Self::Categorical { weights } = self
            && weights.len() > MAX_PARAMETERS
        {
            return Err(error(
                "a statrs distribution accepts at most 4096 parameters",
            ));
        }
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        match self {
            | Self::Categorical { weights } => {
                size_of::<Self>()
                    .saturating_add(weights.capacity().saturating_mul(size_of::<f64>()))
            },
            | _ => size_of::<Self>(),
        }
    }
}

struct NonnegativeDensity(f64);

impl Output for NonnegativeDensity {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn encode(self) -> Result<Value, OpError> {
        if self.0 >= 0.0 {
            Ok(Value::Float(self.0))
        } else {
            Err(error(
                "density or probability mass is undefined at this point",
            ))
        }
    }
}

object_binding!(
    statrs_new,
    "num_statrs_distribution_new",
    "crate::kernels::distributions_ext::StatrsDistribution::from_name",
    "Creates a checked statrs distribution from its kebab-case Rust name and ordered dimensionless parameters. Use exact parameter counts. The categorical family takes a vector of weights.",
    true,
    (name: String, parameters: Vec<f64>) -> NewObject<StatrsDistribution>,
    |_, name, parameters| {
        if parameters.len() > super::MAX_PARAMETERS {
            return Err(error("a statrs distribution accepts at most 4096 parameters"));
        }
        let distribution = StatrsDistribution::from_name(&name, &parameters).map_err(error)?;
        // This call constructs and validates both continuous and discrete families.
        distribution.pdf(0.0).map_err(error)?;
        Ok(NewObject::normalized(distribution))
    }
);

object_binding!(
    statrs_pdf,
    "num_statrs_distribution_pdf",
    "crate::kernels::distributions_ext::StatrsDistribution::pdf",
    "Evaluates the dimensionless density or discrete mass. A density singularity may return positive infinity.",
    true,
    (distribution: Handle<StatrsDistribution>, x: f64) -> NonnegativeDensity,
    |ctx, distribution, x| distribution.with(ctx, |d| d.pdf(x).map(NonnegativeDensity).map_err(error))
);

macro_rules! statrs_query {
    ($fun:ident, $name:literal, $method:ident, $source:literal, $doc:literal, ($($arg:ident: $ty:ty),*) -> $out:ty) => {
        object_binding!(
            $fun, $name, $source, $doc, true,
            (distribution: Handle<StatrsDistribution> $(, $arg: $ty)*) -> $out,
            |ctx, distribution $(, $arg)*| distribution.with(ctx, |d| d.$method($($arg),*).map_err(error))
        );
    }
}

statrs_query!(statrs_cdf, "num_statrs_distribution_cdf", cdf, "crate::kernels::distributions_ext::StatrsDistribution::cdf", "Cumulative probability at the dimensionless query point.", (x: f64) -> f64);
statrs_query!(statrs_sf, "num_statrs_distribution_sf", sf, "crate::kernels::distributions_ext::StatrsDistribution::sf", "Upper-tail probability evaluated directly to avoid subtracting a cdf rounded to one.", (x: f64) -> f64);
statrs_query!(statrs_mean, "num_statrs_distribution_mean", mean, "crate::kernels::distributions_ext::StatrsDistribution::mean", "Returns a finite distribution mean when defined.", () -> f64);
statrs_query!(statrs_variance, "num_statrs_distribution_variance", variance, "crate::kernels::distributions_ext::StatrsDistribution::variance", "Returns a finite distribution variance when defined.", () -> f64);
statrs_query!(statrs_entropy, "num_statrs_distribution_entropy", entropy, "crate::kernels::distributions_ext::StatrsDistribution::entropy", "Returns the discrete or differential entropy in nats when defined.", () -> f64);

object_binding!(
    statrs_discrete,
    "num_statrs_distribution_is_discrete",
    "crate::kernels::distributions_ext::StatrsDistribution::is_discrete",
    "Whether the distribution is discrete.",
    true,
    (distribution: Handle<StatrsDistribution>) -> bool,
    |ctx, distribution| distribution.with(ctx, |d| Ok(d.is_discrete()))
);

object_binding!(
    statrs_signed,
    "num_statrs_distribution_is_signed_discrete",
    "crate::kernels::distributions_ext::StatrsDistribution::is_signed_discrete",
    "Whether the distribution uses signed integer support.",
    true,
    (distribution: Handle<StatrsDistribution>) -> bool,
    |ctx, distribution| distribution.with(ctx, |d| Ok(d.is_signed_discrete()))
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[
            &statrs_new::BINDING,
            &statrs_pdf::BINDING,
            &statrs_cdf::BINDING,
            &statrs_sf::BINDING,
            &statrs_mean::BINDING,
            &statrs_variance::BINDING,
            &statrs_entropy::BINDING,
            &statrs_discrete::BINDING,
            &statrs_signed::BINDING,
        ],
        "numerical_stats",
        DimRule::Dimensionless,
    )
}
