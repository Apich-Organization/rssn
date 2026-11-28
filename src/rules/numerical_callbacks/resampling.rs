//! Seeded resampling with independent sample/statistic dimensions and budgets.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RuleError;
use super::Signature;
use super::ToPrimitive;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::value::TensorValue;
use crate::kernels::random::Rng;
use crate::kernels::stats_ext;

#[derive(Clone, Copy)]
enum Resampling {
    Bootstrap,
    Jackknife,
}
impl EvalDyn for Resampling {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let bootstrap = matches!(self, Self::Bootstrap);
        let mut params = vec![tensor_type(1, d.clone())];
        if bootstrap {
            params.extend([TypeExpr::Int, TypeExpr::Int]);
        }
        params.push(TypeExpr::Int);
        params.push(TypeExpr::Function(
            vec![tensor_type(1, d)],
            Box::new(dimensional_type(e.clone())),
        ));
        let fields = if bootstrap {
            vec![
                ("estimate".into(), dimensional_type(e.clone())),
                ("mean".into(), dimensional_type(e.clone())),
                ("std_error".into(), dimensional_type(e.clone())),
                ("ci_low".into(), dimensional_type(e.clone())),
                ("ci_high".into(), dimensional_type(e.clone())),
                ("replicates".into(), tensor_type(1, e)),
            ]
        } else {
            vec![
                ("estimate".into(), dimensional_type(e.clone())),
                ("bias".into(), dimensional_type(e.clone())),
                ("std_error".into(), dimensional_type(e.clone())),
                ("jackknife_estimates".into(), tensor_type(1, e)),
            ]
        };
        let mut signature = Signature::new(params, TypeExpr::Record(fields));
        signature.type_params = ["D", "E"]
            .into_iter()
            .map(|name| {
                TypeParam {
                    name: name.into(),
                    kind: TypeParamKind::Dim,
                    bounds: vec![],
                }
            })
            .collect();
        effectful(signature)
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let bootstrap = matches!(self, Self::Bootstrap);
        let arity = if bootstrap { 5 } else { 3 };
        let at = arity - 1;
        if args.len() != arity {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(arity as u8),
                got: args.len(),
            });
        }
        let (callback, sig) =
            CallbackRun::new(ctx, &args[at], at, 1, self.signature().params[at].clone())?;
        let sample_dimension = vector_dimension(&sig.params[0], None)?;
        let statistic_dimension = scalar_dimension(&sig.ret)?;
        let sample = vector(&args[0], &sample_dimension, None, 0, ctx)?;
        if sample.is_empty() || (!bootstrap && sample.len() < 2) {
            return Err(OpError::Domain(
                "resampling needs a nonempty sample; jackknife needs at least two values".into(),
            ));
        }
        if !bootstrap && matches!(&sig.params[0], TypeExpr::Tensor { shape: Some(_), .. }) {
            return Err(OpError::Domain("jackknife statistic must accept vectors of both original and leave-one-out lengths; use a rank-one callback domain without a fixed extent".into()));
        }
        if bootstrap {
            vector_dimension(&sig.params[0], Some(sample.len()))?;
        }
        let repetitions = if bootstrap {
            positive_budget(&args[1], 1, "bootstrap replicate count")?
        } else {
            sample.len()
        };
        let budget = positive_budget(&args[at - 1], at - 1, "statistic evaluation budget")?;
        if repetitions.checked_add(1).is_none_or(|n| n > budget) {
            return Err(OpError::BudgetExceeded);
        }
        if sample
            .len()
            .checked_mul(repetitions)
            .is_none_or(|n| n > 100_000_000)
            || sample.len().saturating_add(repetitions) > 8_388_608
        {
            return Err(OpError::Domain(
                "resampling exceeds the finite work/storage budget".into(),
            ));
        }
        let statistic = |values: &[f64]| {
            let input = Value::Tensor(TensorValue {
                shape: vec![values.len()],
                data: values.iter().copied().map(Value::Float).collect(),
                elem: dimensional_type(sample_dimension.clone()),
            });
            callback
                .evaluate(&[input], |value| {
                    check_dimension(ctx, &value, &statistic_dimension, at)?;
                    numeric(&value, at)
                })
                .ok_or_else(|| stats_ext::StatsExtError::Domain("statistic callback failed".into()))
        };
        let fields = if bootstrap {
            let Value::Int(seed) = &args[2] else {
                return Err(OpError::Domain(
                    "bootstrap seed must be an exact nonnegative u64 integer".into(),
                ));
            };
            let seed = seed
                .to_u64()
                .ok_or_else(|| OpError::Domain("bootstrap seed must fit u64".into()))?;
            let mut rng = Rng::new(seed);
            let result =
                stats_ext::bootstrap(&sample, repetitions, &statistic, &mut || rng.uniform());
            callback.finish()?;
            let result = result.map_err(|error| OpError::Domain(error.to_string()))?;
            vec![
                (
                    "estimate".into(),
                    dimensional_float(result.estimate, statistic_dimension.clone()),
                ),
                (
                    "mean".into(),
                    dimensional_float(result.mean, statistic_dimension.clone()),
                ),
                (
                    "std_error".into(),
                    dimensional_float(result.std_error, statistic_dimension.clone()),
                ),
                (
                    "ci_low".into(),
                    dimensional_float(result.ci_low, statistic_dimension.clone()),
                ),
                (
                    "ci_high".into(),
                    dimensional_float(result.ci_high, statistic_dimension.clone()),
                ),
                (
                    "replicates".into(),
                    output_tensor(
                        vec![result.replicates.len()],
                        result.replicates,
                        statistic_dimension.clone(),
                    )?,
                ),
            ]
        } else {
            let result = stats_ext::jackknife(&sample, &statistic);
            callback.finish()?;
            let result = result.map_err(|error| OpError::Domain(error.to_string()))?;
            vec![
                (
                    "estimate".into(),
                    dimensional_float(result.estimate, statistic_dimension.clone()),
                ),
                (
                    "bias".into(),
                    dimensional_float(result.bias, statistic_dimension.clone()),
                ),
                (
                    "std_error".into(),
                    dimensional_float(result.std_error, statistic_dimension.clone()),
                ),
                (
                    "jackknife_estimates".into(),
                    output_tensor(
                        vec![result.jackknife_estimates.len()],
                        result.jackknife_estimates,
                        statistic_dimension,
                    )?,
                ),
            ]
        };
        Ok(Value::Record(RecordValue { fields }))
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        installer,
        "num_bootstrap",
        "crate::kernels::stats_ext::bootstrap",
        "Seeded nonparametric bootstrap: sample vector, positive replicate count, nonnegative u64 seed, statistic-evaluation budget, statistic callback. The callback accepts a sample-unit vector and returns a finite scalar with independently declared statistic units; every summary/replicate retains those units. Budget must cover the original estimate and every replicate. Work is capped at 100000000 resampled entries and storage at 8388608 values. Percentile confidence bounds are 95 percent and are statistical estimates, not coverage guarantees.",
        5,
        Resampling::Bootstrap,
    )?;
    register(
        installer,
        "num_jackknife",
        "crate::kernels::stats_ext::jackknife",
        "Leave-one-out jackknife: sample vector, statistic-evaluation budget, statistic callback. The callback domain must accept variable vector lengths; its output units independently determine estimate/bias/standard-error/replicate units. Budget covers the original estimate plus every leave-one-out sample. Work/storage caps match bootstrap. Callback errors propagate unchanged and stop further host calls.",
        3,
        Resampling::Jackknife,
    )
}
