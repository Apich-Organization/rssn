//! Distribution enum specializations with explicit sample reference measures.
#![allow(clippy::wildcard_imports)] // Shared typed method helpers.
use super::*;
use crate::kernels::distributions::Distribution;
use crate::kernels::distributions::MultivariateDistribution;
use crate::kernels::random::Rng;
struct Univariate {
    value: Distribution,
    scale: f64,
}
struct Multivariate {
    value: MultivariateDistribution,
    reference: f64,
}
impl NativeType for Univariate {
    const KIND: &'static str = "rssn/native/discrete/distribution/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        self.value.validate().map_err(error)?;
        if !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(error("invalid distribution reference scale"));
        }
        Ok(())
    }
}
impl NativeType for Multivariate {
    const KIND: &'static str = "rssn/native/discrete/multivariate/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        if self.value.dim() == 0
            || self.value.dim() > 32
            || !self.reference.is_finite()
            || self.reference <= 0.0
        {
            return Err(error("invalid multivariate dimension/reference scale"));
        }
        let mean = self.value.mean().map_err(error)?;
        let log = self.value.log_pdf(&mean).map_err(error)?;
        if log.is_nan() {
            return Err(error("invalid multivariate parameters"));
        }
        Ok(())
    }
}
fn parameter(
    value: &Input,
    name: &str,
) -> Result<f64, OpError> {
    api::normalized(&value.0)?;
    numeric(api::field(&value.0, name)?, 0)
}
fn distribution(
    kind: &str,
    parameters: &Input,
) -> Result<Distribution, OpError> {
    let get = |name| parameter(parameters, name);
    let result = match kind {
        | "Normal" => {
            Distribution::Normal {
                mu: get("mu")?,
                sigma: get("sigma")?,
            }
        },
        | "Uniform" => {
            Distribution::Uniform {
                lo: get("lo")?,
                hi: get("hi")?,
            }
        },
        | "Exponential" => Distribution::Exponential { rate: get("rate")? },
        | "Gamma" => {
            Distribution::Gamma {
                shape: get("shape")?,
                rate: get("rate")?,
            }
        },
        | "Beta" => {
            Distribution::Beta {
                alpha: get("alpha")?,
                beta: get("beta")?,
            }
        },
        | "ChiSquared" => Distribution::ChiSquared { df: get("df")? },
        | "StudentT" => Distribution::StudentT { df: get("df")? },
        | "Poisson" => {
            Distribution::Poisson {
                lambda: get("lambda")?,
            }
        },
        | "Binomial" => {
            Distribution::Binomial {
                n: <u64 as Arg>::decode(api::field(&parameters.0, "n")?, 0)?,
                p: get("p")?,
            }
        },
        | "LogNormal" => {
            Distribution::LogNormal {
                mu: get("mu")?,
                sigma: get("sigma")?,
            }
        },
        | "Weibull" => {
            Distribution::Weibull {
                shape: get("shape")?,
                scale: get("scale")?,
            }
        },
        | "Cauchy" => {
            Distribution::Cauchy {
                loc: get("loc")?,
                scale: get("scale")?,
            }
        },
        | "F" => {
            Distribution::F {
                df1: get("df1")?,
                df2: get("df2")?,
            }
        },
        | "NegativeBinomial" => {
            Distribution::NegativeBinomial {
                r: get("r")?,
                p: get("p")?,
            }
        },
        | "Geometric" => Distribution::Geometric { p: get("p")? },
        | _ => return Err(error("unknown native distribution variant")),
    };
    result.validate().map_err(error)?;
    Ok(result)
}
const fn discrete(value: &Distribution) -> bool {
    matches!(
        value,
        Distribution::Poisson { .. }
            | Distribution::Binomial { .. }
            | Distribution::NegativeBinomial { .. }
            | Distribution::Geometric { .. }
    )
}
fn query(
    a: &Handle<Univariate>,
    ctx: &CallCtx<'_>,
    x: &Input,
) -> Result<(f64, DimExpr), OpError> {
    let d = object_dimension(a, ctx)?;
    api::same_dimension(&d, &api::dimension(&x.0)?)?;
    Ok((api::extended_numeric(&x.0)?, d))
}
fn multi_query(
    a: &Handle<Multivariate>,
    ctx: &CallCtx<'_>,
    x: &Input,
) -> Result<(Vec<f64>, DimExpr), OpError> {
    let d = object_dimension(a, ctx)?;
    api::same_dimension(&d, &api::dimension(&x.0)?)?;
    let x = <Vec<f64> as Arg>::decode(&x.0, 0)?;
    a.with(ctx, |o| {
        if x.len() == o.value.dim() {
            Ok(())
        } else {
            Err(error("multivariate observation dimension mismatch"))
        }
    })?;
    Ok((x, d))
}
object_binding!(distribution_new,"num_native_distribution_new","","Construct a native distribution enum in normalized sample coordinates, then restore an explicit positive sample reference scale. All15 named variants are supported. Continuous density has inverse sample units; discrete probability mass remains dimensionless.",false,(kind:String,parameters:Input,sample_scale:Input)->NewObject<Univariate>,|_ctx:&CallCtx<'_>,kind:String,parameters:Input,sample_scale:Input|{let(scale,d)=api::reference(&sample_scale.0)?;let value=distribution(&kind,&parameters)?;Ok(NewObject::new(Univariate{value,scale},vec![d]))});
object_binding!(distribution_cdf,"num_native_distribution_cdf","crate::kernels::distributions::Distribution::cdf","Native CDF of a sample quantity sharing the constructor reference domain.",false,(a:Handle<Univariate>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<Univariate>,x:Input|{let(x,_)=query(&a,ctx,&x)?;a.with(ctx,|o|Ok(o.value.cdf(x/o.scale)))});
object_binding!(distribution_pdf,"num_native_distribution_pdf","crate::kernels::distributions::Distribution::pdf","Native normalized density restored by the explicit reference scale. For discrete variants this is probability mass and has no physical units.",false,(a:Handle<Univariate>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Univariate>,x:Input|{let(x,d)=query(&a,ctx,&x)?;a.with(ctx,|o|{let density=o.value.pdf(x/o.scale);if discrete(&o.value){extended_quantity(density,DimExpr::Dimless)}else{extended_quantity(density/o.scale,canonical(&(DimExpr::Inverse(Box::new(d))))?)}})});
object_binding!(distribution_mean,"num_native_distribution_mean","crate::kernels::distributions::Distribution::mean","Native mean with sample units; mathematically undefined NaN moments report a domain error.",false,(a:Handle<Univariate>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Univariate>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.value.mean()*o.scale,d))});
object_binding!(distribution_variance,"num_native_distribution_variance","crate::kernels::distributions::Distribution::variance","Native variance with squared sample units; undefined moments report an explicit domain error.",false,(a:Handle<Univariate>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Univariate>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|extended_quantity(o.value.variance()*o.scale*o.scale,d))});
object_binding!(distribution_quantile,"num_native_distribution_quantile","crate::kernels::distributions::Distribution::quantile","Native inverse CDF restored to sample units. Valid infinite endpoint quantiles are preserved.",false,(a:Handle<Univariate>,probability:f64)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Univariate>,probability:f64|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.value.quantile(probability).map_err(error)?*o.scale,d))});
object_binding!(distribution_sample,"num_native_distribution_sample","crate::kernels::distributions::Distribution::sample","Native FnMut uniform-source specialization using the shared session-owned Rng. Both distribution and random stream kinds are checked; result retains sample units.",false,(a:Handle<Univariate>,rng:Handle<Rng>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Univariate>,rng:Handle<Rng>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|rng.with_mut(ctx,|rng|extended_quantity(o.value.sample(&mut||rng.uniform()).map_err(error)?*o.scale,d)))});
object_binding!(multivariate_normal_new,"num_native_multivariate_normal_new","","Construct a native multivariate normal from homogeneous mean units and squared covariance units. A positive reference length/scale with mean units defines the dimensionless log-density reference measure. Dimension is bounded to32.",false,(mu:Input,cov:Input,reference:Input)->NewObject<Multivariate>,|_ctx:&CallCtx<'_>,mu:Input,cov:Input,reference:Input|{let(r,d)=api::reference(&reference.0)?;api::same_dimension(&d,&api::dimension(&mu.0)?)?;api::same_dimension(&canonical(&(d.clone().pow(Exp::Rat(Rational::int(2)))))?,&api::dimension(&cov.0)?)?;let mu=<Vec<f64> as Arg>::decode(&mu.0,0)?;limit(mu.len(),32)?;let matrix=<Vec<Vec<f64>> as Arg>::decode(&cov.0,1)?;if matrix.len()!=mu.len()||matrix.iter().any(|row|row.len()!=mu.len()){return Err(error("multivariate covariance shape disagrees with mean"));}let value=MultivariateDistribution::MultivariateNormal{mu,cov:matrix.into_iter().flatten().collect()};Ok(NewObject::new(Multivariate{value,reference:r},vec![d]))});
object_binding!(dirichlet_new,"num_native_dirichlet_new","","Construct the native dimensionless Dirichlet simplex distribution with positive concentration parameters.",true,(alpha:Vec<f64>)->NewObject<Multivariate>,|_ctx:&CallCtx<'_>,alpha:Vec<f64>|{limit(alpha.len(),32)?;Ok(NewObject::new(Multivariate{value:MultivariateDistribution::Dirichlet{alpha},reference:1.0},vec![DimExpr::Dimless]))});
object_binding!(multivariate_dim,"num_native_multivariate_dim","crate::kernels::distributions::MultivariateDistribution::dim","Native dimension query.",false,(a:Handle<Multivariate>)->usize,|ctx:&CallCtx<'_>,a:Handle<Multivariate>|a.with(ctx,|o|Ok(o.value.dim())));
object_binding!(multivariate_mean,"num_native_multivariate_mean","crate::kernels::distributions::MultivariateDistribution::mean","Native mean vector retains coordinate units.",false,(a:Handle<Multivariate>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivariate>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|quantities(o.value.mean().map_err(error)?,&(d)))});
object_binding!(multivariate_cov,"num_native_multivariate_cov","crate::kernels::distributions::MultivariateDistribution::cov","Native covariance matrix retains squared coordinate units.",false,(a:Handle<Multivariate>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivariate>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|{let n=o.value.dim();quantities_matrix(o.value.cov().map_err(error)?.chunks_exact(n).map(<[f64]>::to_vec).collect(),&(d))})});
object_binding!(multivariate_pdf,"num_native_multivariate_pdf","crate::kernels::distributions::MultivariateDistribution::pdf","Native density with inverse coordinate units raised to the checked distribution dimension; Dirichlet remains dimensionless.",false,(a:Handle<Multivariate>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivariate>,x:Input|{let(x,d)=multi_query(&a,ctx,&x)?;a.with(ctx,|o|extended_quantity(o.value.pdf(&x).map_err(error)?,canonical(&(d.pow(Exp::Rat(Rational::int(-(o.value.dim()as i64))))))?))});
object_binding!(multivariate_log_pdf,"num_native_multivariate_log_pdf","crate::kernels::distributions::MultivariateDistribution::log_pdf","Native log density relative to the constructor explicit reference measure, avoiding a logarithm of dimensioned density.",false,(a:Handle<Multivariate>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<Multivariate>,x:Input|{let(x,_)=multi_query(&a,ctx,&x)?;a.with(ctx,|o|Ok(o.value.log_pdf(&x).map_err(error)?+(o.value.dim()as f64)*o.reference.ln()))});
object_binding!(multivariate_sample,"num_native_multivariate_sample","crate::kernels::distributions::MultivariateDistribution::sample","Native FnMut uniform-source specialization using the shared authenticated random stream; output retains coordinate units.",false,(a:Handle<Multivariate>,rng:Handle<Rng>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivariate>,rng:Handle<Rng>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|rng.with_mut(ctx,|rng|quantities(o.value.sample(&mut||rng.uniform()).map_err(error)?,&(d))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &distribution_new::BINDING,
            &distribution_cdf::BINDING,
            &distribution_pdf::BINDING,
            &distribution_mean::BINDING,
            &distribution_variance::BINDING,
            &distribution_quantile::BINDING,
            &distribution_sample::BINDING,
            &multivariate_normal_new::BINDING,
            &dirichlet_new::BINDING,
            &multivariate_dim::BINDING,
            &multivariate_mean::BINDING,
            &multivariate_cov::BINDING,
            &multivariate_pdf::BINDING,
            &multivariate_log_pdf::BINDING,
            &multivariate_sample::BINDING,
        ],
        "numerical_vector",
    )
}
