//! Remaining stateless native physics functions with measured SI carriers.
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
fn scalar(
    v: &Value,
    d: &DimExpr,
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(v)?, d)?;
    super::numeric(v, 0)
}
object_binding!(oscillator,"num_native_damped_harmonic_oscillator","crate::sim::classical::damped_harmonic_oscillator","Actual native damped oscillator. Amplitude/result retain an arbitrary physical domain; angular frequency and damping have inverse-Time units, phase is dimensionless, elapsed time is Time.",false,(amplitude:Input,omega0:Input,gamma:Input,phase:f64,time:Input)->Dynamic,|_ctx,amplitude,omega0,gamma,phase,time|{let d=api::dimension(&amplitude.0)?;let a=scalar(&amplitude.0,&d)?;let frequency=DimExpr::base("Time").pow(-1);let omega=scalar(&omega0.0,&frequency)?;let gamma=scalar(&gamma.0,&frequency)?;let time=scalar(&time.0,&DimExpr::base("Time"))?;if omega<0.0||gamma<0.0{return Err(api::error("frequency and damping must be nonnegative"));}Ok(Dynamic(api::quantity(crate::sim::classical::damped_harmonic_oscillator(a,omega,gamma,phase,time),d)?))});
object_binding!(diag_cell,"num_native_ewald_diag_cell","crate::sim::ewald::diag_cell","Construct the actual diagonal Ewald cell from three measured lengths; returns a measured3x3tensor.",false,(lengths:Input)->Dynamic,|_ctx,lengths|{let Value::Tensor(t)=&lengths.0 else{return Err(api::error("expected three cell lengths"));};
if t.shape!=[3]||t.data.len()!=3{return Err(api::error("cell length shape mismatch"));}let d=DimExpr::base("Length");api::same_dimension(&api::dimension(&lengths.0)?,&d)?;let mut v=[0.0;3];for(i,x)in t.data.iter().enumerate(){v[i]=crate::rules::numerical_abi::numeric(x,0)?;if v[i]<=0.0{return Err(api::error("cell lengths must be positive"));}}let cell=crate::sim::ewald::diag_cell(v);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3,3],elem:TypeExpr::Float.with_dim(d.clone()),data:cell.into_iter().flatten().map(|x|api::quantity(x,d.clone())).collect::<Result<_,_>>()?})))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[&oscillator::BINDING, &diag_cell::BINDING],
        "numerical_vector",
    )
}
