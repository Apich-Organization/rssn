//! Physical scaling of the legacy TM FDTD model and explicit native NPY output.
use super::native_object_api::Dynamic;
use super::native_object_api::Input;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
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
    data: Vec<f64>,
    shape: Vec<usize>,
    scale: f64,
    dimension: DimExpr,
) -> Result<Value, OpError> {
    let data = data
        .into_iter()
        .map(|x| api::quantity(api::scale_product(x, &[scale], &[])?, dimension.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float.with_dim(dimension),
        data,
    }))
}
object_binding!(run,"num_native_fdtd_model_run","crate::sim::models_fdtd_electrodynamics::run_fdtd_simulation","Run the actual legacy2D TM model. Explicit positive Length, wave-speed and electric-field references define grid pitch L and physical timestep L/(2c), matching native0.5updates. Physical source frequency is converted to cycles per native step. Returns measured field snapshots and their post-update times. Native10cell boundary damping is retained; use Yee for general boundary/material models.",false,(width:usize,height:usize,steps:usize,source:Vec<usize>,frequency:Input,length_reference:Input,speed_reference:Input,electric_reference:Input)->Dynamic,|_ctx,width,height,steps,source,frequency,length_reference,speed_reference,electric_reference|{
 if width<2||height<2||steps>65536||source.len()!=2||source[0]>=width||source[1]>=height{return Err(api::error("legacy FDTD dimensions/source/step bounds invalid"));}
 let cells=width.checked_mul(height).filter(|n|*n<=16384).ok_or_else(||api::error("legacy FDTD cell budget exceeded"))?;
 let snapshots=steps.div_ceil(5);if cells.checked_mul(steps).is_none_or(|n|n>10000000)||cells.checked_mul(snapshots).is_none_or(|n|n>65536){return Err(api::error("legacy FDTD work/snapshot allocation budget exceeded"));}
 let l=reference(&length_reference.0,&DimExpr::base("Length"))?;
 let c=reference(&speed_reference.0,&DimExpr::Mul(vec![DimExpr::base("Length"),DimExpr::base("Time").pow(-1)]))?;
 let e=reference(&electric_reference.0,&electric())?;let dt=api::scale_product(0.5,&[l],&[c])?;if dt<=0.0{return Err(api::error("physical FDTD timestep is not representable"));}
 api::scale_product(steps as f64,&[dt],&[])?;
 api::same_dimension(&api::dimension(&frequency.0)?,&DimExpr::base("Time").pow(-1))?;
 let f=crate::rules::numerical_abi::numeric(&frequency.0,0)?;if f<0.0{return Err(api::error("source frequency must be nonnegative"));}
 let parameters=fdtd::FdtdParameters{width,height,time_steps:steps,source_pos:(source[0],source[1]),source_freq:api::scale_product(f,&[dt],&[])?};
 let values=fdtd::run_fdtd_simulation(&parameters);if values.len()!=snapshots{return Err(api::error("legacy FDTD snapshot count mismatch"));}
 let times=(0..snapshots).map(|i|(5*i+1)as f64).collect();
 let values=values.into_iter().flat_map(IntoIterator::into_iter).collect();
 Ok(Dynamic(api::record(vec![("times",field(times,vec![snapshots],dt,DimExpr::base("Time"))?),("electric_field",field(values,vec![snapshots,width,height],e,electric())?)])))
});
object_binding!(save,"num_native_fdtd_model_save","crate::sim::models_fdtd_electrodynamics::simulate_and_save_final_state","IO: call the actual legacy1D zero-initial-field FDTD NPY writer. The explicit native-normalized contract states that the NPY numeric array has the native dimensionless representation; NPY carries no unit metadata. Grid/work/path bounds precede native allocation and file IO.",true,(cells:usize,steps:usize,filename:String,representation:String)->bool,|_ctx,cells,steps,filename,representation|{
 if representation!="native-normalized"||!(2..=65536).contains(&cells)||steps>65536||cells.checked_mul(steps).is_none_or(|n|n>10000000)||filename.is_empty()||filename.len()>4096{return Err(api::error("native FDTD file contract or grid/work/path bounds invalid"));}
 fdtd::simulate_and_save_final_state(cells,steps,&filename).map_err(api::error)?;Ok(true)
});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, &[&run::BINDING, &save::BINDING], "numerical_vector")
}
