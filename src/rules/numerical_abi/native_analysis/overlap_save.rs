//! Session-owned overlap-save plans for measured real signals.

use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::fftconv::OverlapSavePlan;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::MeasuredVector;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};

const TAPS_LIMIT: usize = 65_536;
const ROWS_LIMIT: usize = 4096;
const VALUES_LIMIT: usize = object::VALUE_LIMIT;

impl NativeType for OverlapSavePlan {
    const KIND: &'static str = "rssn/native/fft/overlap-save/v1";
    const MAX_BYTES: usize = 24 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        self.retained_budget()
    }
}

fn result_dimension(
    handle: &Handle<OverlapSavePlan>,
    ctx: &CallCtx<'_>,
    samples: DimExpr,
) -> Result<DimExpr, OpError> {
    let taps = handle
        .dimensions(ctx)?
        .into_iter()
        .next()
        .ok_or_else(|| object::error("overlap-save plan lacks tap units"))?;
    Ok(DimExpr::Mul(vec![samples, taps]))
}

fn result_len(
    plan: &OverlapSavePlan,
    input_len: usize,
) -> Result<usize, OpError> {
    plan.output_len(input_len).map_err(object::error)
}

fn value_tensor(
    values: Vec<f64>,
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.into_iter().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(dimension),
    })
}

object_binding!(
    new,
    "num_native_overlap_save_plan_new",
    "crate::kernels::fftconv::OverlapSavePlan::new",
    "Prepare the spectrum of a measured real filter with at most 65536 taps. The immutable handle retains the tap units and can be reused for several recordings in this session.",
    false,
    (taps: MeasuredVector<{TAPS_LIMIT}>) -> NewObject<OverlapSavePlan>,
    |_, taps| {
        let plan = OverlapSavePlan::new(&taps.values).map_err(object::error)?;
        Ok(NewObject::new(plan, vec![taps.dimension]))
    }
);

object_binding!(
    apply,
    "num_native_overlap_save_plan_apply",
    "crate::kernels::fftconv::OverlapSavePlan::apply",
    "Compute the full linear convolution of a real recording with the prepared filter. The output has sample units times tap units. Empty recordings give an empty measured tensor. The result has at most 262144 values.",
    false,
    (plan: Handle<OverlapSavePlan>, samples: MeasuredVector<{VALUES_LIMIT}>) -> Dynamic,
    |ctx, plan, samples| {
        let dimension = result_dimension(&plan, ctx, samples.dimension)?;
        let result = plan.with(ctx, |plan| {
            if result_len(plan, samples.values.len())? > VALUES_LIMIT {
                return Err(object::error("overlap-save result exceeds the value budget"));
            }
            plan.apply(&samples.values).map_err(object::error)
        })?;
        Ok(Dynamic(value_tensor(result, dimension)))
    }
);

object_binding!(
    apply_batch,
    "num_native_overlap_save_plan_apply_batch",
    "crate::kernels::fftconv::OverlapSavePlan::apply_batch",
    "Convolve up to 4096 real recordings sharing one sample dimension. Supply flat samples and zero-based row offsets; repeated offsets represent empty recordings. Returns output offsets and flat values with sample units times tap units. The whole batch has at most 262144 output values.",
    false,
    (plan: Handle<OverlapSavePlan>, samples: MeasuredVector<{VALUES_LIMIT}>, offsets: Vec<usize>) -> Dynamic,
    |ctx, plan, samples, offsets| {
        if offsets.first() != Some(&0)
            || offsets.last() != Some(&samples.values.len())
            || offsets.len() > ROWS_LIMIT + 1
            || offsets.windows(2).any(|bounds| bounds[0] > bounds[1])
        {
            return Err(object::error("invalid overlap-save batch row offsets"));
        }
        let dimension = result_dimension(&plan, ctx, samples.dimension)?;
        let mut result_offsets = Vec::with_capacity(offsets.len());
        result_offsets.push(Value::int(0));
        let results = plan.with(ctx, |plan| {
            let mut total = 0usize;
            let mut rows = Vec::with_capacity(offsets.len().saturating_sub(1));
            for bounds in offsets.windows(2) {
                let row = &samples.values[bounds[0]..bounds[1]];
                total = total.saturating_add(result_len(plan, row.len())?);
                if total > VALUES_LIMIT {
                    return Err(object::error("overlap-save batch result exceeds the value budget"));
                }
                rows.push(row);
                result_offsets.push(Value::int(total));
            }
            plan.apply_batch(&rows).map_err(object::error)
        })?;
        let output_offsets = Value::Tensor(TensorValue {
            shape: vec![result_offsets.len()],
            data: result_offsets,
            elem: TypeExpr::Int,
        });
        let values = value_tensor(results.into_iter().flatten().collect(), dimension);
        Ok(Dynamic(object::record(vec![("offsets", output_offsets), ("values", values)])))
    }
);

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(
        installer,
        &[&new::BINDING, &apply::BINDING, &apply_batch::BINDING],
        "numerical_vector",
    )
}
