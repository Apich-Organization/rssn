//! Physical scaling of the legacy TM FDTD model and explicit native NPY output.
use super::native_object_api::Dynamic;
use super::native_object_api::Input;
use super::native_object_api::NativeOutput;
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
use crate::kernel::value::TensorValue;
use crate::sim::models_fdtd_electrodynamics as fdtd;
fn electric() -> DimExpr {
    DimExpr::Mul(vec![
        DimExpr::base("Mass"),
        DimExpr::base("Length"),
        DimExpr::base("Current").pow(-1),
        DimExpr::base("Time").pow(-3),
    ])
}
fn reference(
    v: &Value,
    d: &DimExpr,
) -> Result<f64, OpError> {
    let (x, actual) = api::reference(v)?;
    api::same_dimension(&actual, d)?;
    Ok(x)
}
fn field(
    values: impl IntoIterator<Item = f64>,
    shape: Vec<usize>,
    scale: f64,
    dimension: DimExpr,
) -> Result<Value, OpError> {
    let mut data = Vec::with_capacity(shape.iter().product());
    for x in values {
        data.push(api::quantity(
            api::scale_product(x, &[scale], &[])?,
            dimension.clone(),
        )?);
    }
    Ok(Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float.with_dim(dimension),
        data,
    }))
}
struct PhysicalModel {
    parameters: fdtd::FdtdParameters,
    dt: f64,
    electric_scale: f64,
}

struct FinalField {
    time: f64,
    electric_field: Value,
}

impl NativeOutput for FinalField {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            (
                "time".into(),
                TypeExpr::Float.with_dim(DimExpr::base("Time")),
            ),
            (
                "electric_field".into(),
                api::tensor_type(2, TypeExpr::Float.with_dim(electric())),
            ),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("time", api::quantity(self.time, DimExpr::base("Time"))?),
            ("electric_field", self.electric_field),
        ]))
    }
}

fn prepare(
    [width, height]: [usize; 2],
    steps: usize,
    source: &[usize],
    frequency: &Input,
    [length_reference, speed_reference, electric_reference]: [&Input; 3],
    sampled: bool,
) -> Result<PhysicalModel, OpError> {
    let minimum = if sampled { 2 } else { 1 };
    let [row, column] = *source else {
        return Err(api::error("FDTD source requires two grid indices"));
    };
    if width < minimum || height < minimum || row >= width || column >= height || steps > 65_536 {
        return Err(api::error("FDTD dimensions, source or step bounds invalid"));
    }
    let cells = width
        .checked_mul(height)
        .filter(|n| *n <= 16_384)
        .ok_or_else(|| api::error("FDTD cell budget exceeded"))?;
    if cells.checked_mul(steps).is_none_or(|n| n > 10_000_000)
        || (sampled
            && cells
                .checked_mul(steps.div_ceil(5))
                .is_none_or(|n| n > 65_536))
    {
        return Err(api::error(
            "FDTD work or snapshot allocation budget exceeded",
        ));
    }
    let length = reference(&length_reference.0, &DimExpr::base("Length"))?;
    let speed = reference(
        &speed_reference.0,
        &DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
    )?;
    let electric_scale = reference(&electric_reference.0, &electric())?;
    let dt = api::scale_product(0.5, &[length], &[speed])?;
    if dt <= 0.0 {
        return Err(api::error("physical FDTD timestep is not representable"));
    }
    api::same_dimension(
        &api::dimension(&frequency.0)?,
        &DimExpr::base("Time").pow(-1),
    )?;
    let frequency = crate::rules::numerical_abi::numeric(&frequency.0, 0)?;
    if frequency < 0.0 {
        return Err(api::error("source frequency must be nonnegative"));
    }
    let cycles = frequency * dt;
    // A product of finite binary floats that overflows is an exact integer in step units.
    let source_freq = if cycles.is_infinite() {
        0.0
    } else {
        cycles
    };
    Ok(PhysicalModel {
        parameters: fdtd::FdtdParameters {
            width,
            height,
            time_steps: steps,
            source_pos: (row, column),
            source_freq,
        },
        dt,
        electric_scale,
    })
}

object_binding!(
    run, "num_native_fdtd_model_run", "crate::sim::models_fdtd_electrodynamics::run_fdtd_simulation",
    "Advance the legacy two-dimensional TM model from zero fields. The length and wave-speed references define the physical time step L/(2c). The physical source frequency is sampled at that step. Returns the measured electric field after steps 1, 6, 11 and so on, and the corresponding times. The soft source precedes the ten-cell boundary damping.",
    false,
    (width: usize, height: usize, steps: usize, source: Vec<usize>, frequency: Input, length_reference: Input, speed_reference: Input, electric_reference: Input) -> Dynamic,
    |_ctx, width, height, steps, source, frequency, length_reference, speed_reference, electric_reference| {
        let m = prepare([width, height], steps, &source, &frequency,
                        [&length_reference, &speed_reference, &electric_reference], true)?;
        let snapshots = steps.div_ceil(5);
        let values = fdtd::run_fdtd_simulation(&m.parameters);
        if values.len() != snapshots {
            return Err(api::error("FDTD snapshot count mismatch"));
        }
        let times = (0..snapshots).map(|i| (5 * i + 1) as f64);
        let values = values.into_iter().flat_map(IntoIterator::into_iter);
        Ok(Dynamic(api::record(vec![
            ("times", field(times, vec![snapshots], m.dt, DimExpr::base("Time"))?),
            ("electric_field", field(values, vec![snapshots, width, height], m.electric_scale, electric())?),
        ])))
    }
);
object_binding!(
    run_final, "num_native_fdtd_model_run_final",
    "crate::sim::models_fdtd_electrodynamics::run_fdtd_final_field",
    "Advance the legacy two-dimensional TM model and return its electric field at the last requested step. The field has shape (width, height), and time is steps multiplied by L/(2c). Zero steps return the initial zero field. The final-state mode admits one-cell axes and uses the same source, wave-speed and damping rules as the snapshot mode.",
    false,
    (width: usize, height: usize, steps: usize, source: Vec<usize>, frequency: Input, length_reference: Input, speed_reference: Input, electric_reference: Input) -> FinalField,
    |_ctx, width, height, steps, source, frequency, length_reference, speed_reference, electric_reference| {
        let m = prepare([width, height], steps, &source, &frequency,
                        [&length_reference, &speed_reference, &electric_reference], false)?;
        let end_time = api::scale_product(steps as f64, &[m.dt], &[])?;
        let values = fdtd::run_fdtd_final_field(&m.parameters).map_err(api::error)?;
        Ok(FinalField {
            time: end_time,
            electric_field: field(values, vec![width, height], m.electric_scale, electric())?,
        })
    }
);
object_binding!(save,"num_native_fdtd_model_save","crate::sim::models_fdtd_electrodynamics::simulate_and_save_final_state","IO: call the actual legacy1D zero-initial-field FDTD NPY writer. The explicit native-normalized contract states that the NPY numeric array has the native dimensionless representation; NPY carries no unit metadata. Grid/work/path bounds precede native allocation and file IO.",true,(cells:usize,steps:usize,filename:String,representation:String)->bool,|_ctx,cells,steps,filename,representation|{
 if representation!="native-normalized"||!(2..=65536).contains(&cells)||steps>65536||cells.checked_mul(steps).is_none_or(|n|n>10000000)||filename.is_empty()||filename.len()>4096{return Err(api::error("native FDTD file contract or grid/work/path bounds invalid"));}
 fdtd::simulate_and_save_final_state(cells,steps,&filename).map_err(api::error)?;Ok(true)
});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[&run::BINDING, &run_final::BINDING, &save::BINDING],
        "numerical_vector",
    )
}
