//! Prepared single-level wavelet transformations with measured samples.
use super::Installer;
use super::RuleError;
use super::native_object_api as api;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::wavelets::Boundary;
use crate::kernels::wavelets::DwtPlan;
use crate::kernels::wavelets::Wavelet;
use api::Dynamic;
use api::Handle;
use api::MeasuredVector;
use api::NativeType;
use api::NewObject;
use api::object_binding;

const MAX_SAMPLES: usize = 262_144;

impl NativeType for DwtPlan {
    const KIND: &'static str = "rssn/native/dwt-plan/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        if self.signal_len() > MAX_SAMPLES {
            return Err(api::error("wavelet signal exceeds the sample budget"));
        }
        Ok(())
    }
}

type Signal = MeasuredVector<MAX_SAMPLES>;

fn measured(
    samples: &[f64],
    dimension: &DimExpr,
) -> Result<Value, OpError> {
    Ok(Value::Tensor(TensorValue {
        shape: vec![samples.len()],
        elem: TypeExpr::Float.with_dim(dimension.clone()),
        data: samples
            .iter()
            .map(|&x| api::quantity(x, dimension.clone()))
            .collect::<Result<_, _>>()?,
    }))
}

object_binding!(
    dwt_plan_new,
    "num_analysis_wavelets_dwt_plan_new",
    "crate::kernels::wavelets::DwtPlan::new",
    "Prepares one-level analysis and synthesis for a fixed nonzero signal length. Use a wavelet such as Haar or Db10 and one of Periodic, Zero, Symmetric or Reflect. Shares filter taps and boundary layout across independent signals.",
    true,
    (wavelet: Wavelet, boundary: Boundary, samples: usize) -> NewObject<DwtPlan>,
    |_, wavelet, boundary, samples| {
        if samples > MAX_SAMPLES {
            return Err(api::error("wavelet signal exceeds the sample budget"));
        }
        DwtPlan::new(wavelet, boundary, samples).map(NewObject::normalized)
            .map_err(|error| api::error(format!("{error:?}")))
    }
);

object_binding!(
    dwt_signal_len,
    "num_analysis_wavelets_dwt_plan_signal_len",
    "crate::kernels::wavelets::DwtPlan::signal_len",
    "Returns the signal length fixed by the plan and used for synthesis output.",
    true,
    (plan: Handle<DwtPlan>) -> usize,
    |ctx, plan| plan.with(ctx, |plan| Ok(plan.signal_len()))
);

object_binding!(
    dwt_coeff_len,
    "num_analysis_wavelets_dwt_plan_coeff_len",
    "crate::kernels::wavelets::DwtPlan::coeff_len",
    "Returns the number of entries in each approximation and detail band for this signal length and boundary.",
    true,
    (plan: Handle<DwtPlan>) -> usize,
    |ctx, plan| plan.with(ctx, |plan| Ok(plan.coeff_len()))
);

object_binding!(
    dwt_analyze,
    "num_analysis_wavelets_dwt_plan_analyze",
    "crate::kernels::wavelets::DwtPlan::analyze_into",
    "Computes one level of approximation and detail for the plan's input length. The two vectors both retain the physical dimension of the original signal.",
    false,
    (plan: Handle<DwtPlan>, signal: Signal) -> Dynamic,
    |ctx, plan, signal| plan.with(ctx, |plan| {
        let mut approximation = f64_scratch(plan.coeff_len());
        let mut detail = f64_scratch(plan.coeff_len());
        approximation.resize(plan.coeff_len(), 0.0);
        detail.resize(plan.coeff_len(), 0.0);
        plan.analyze_into(&signal.values, &mut approximation, &mut detail)
            .map_err(|error| api::error(format!("{error:?}")))?;
        Ok(Dynamic(api::record(vec![
            ("approximation", measured(&approximation, &signal.dimension)?),
            ("detail", measured(&detail, &signal.dimension)?),
        ])))
    })
);

object_binding!(
    dwt_synthesize,
    "num_analysis_wavelets_dwt_plan_synthesize",
    "crate::kernels::wavelets::DwtPlan::synthesize_into",
    "Reconstructs the original signal length from two equal-dimension coefficient bands of the plan's required length. The output retains the shared dimension.",
    false,
    (plan: Handle<DwtPlan>, approximation: Signal, detail: Signal) -> Dynamic,
    |ctx, plan, approximation, detail| {
        api::same_dimension(&approximation.dimension, &detail.dimension)?;
        plan.with(ctx, |plan| {
            let mut signal = f64_scratch(plan.signal_len());
            signal.resize(plan.signal_len(), 0.0);
            plan.synthesize_into(&approximation.values, &detail.values, &mut signal)
                .map_err(|error| api::error(format!("{error:?}")))?;
            Ok(Dynamic(measured(&signal, &approximation.dimension)?))
        })
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &dwt_plan_new::BINDING,
            &dwt_signal_len::BINDING,
            &dwt_coeff_len::BINDING,
            &dwt_analyze::BINDING,
            &dwt_synthesize::BINDING,
        ],
        "numerical_vector",
    )
}
