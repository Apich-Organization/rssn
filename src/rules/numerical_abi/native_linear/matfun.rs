//! Sparse and matrix-free matrix exponential actions with dimensionless A*t.
use super::super::native_object_api::Maybe;
use super::Array;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::RuleError;
use super::Value;
use super::WORK;
use super::api;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::sparse::Csr;
use super::sparse::result;
use super::sparse::vector;
use crate::kernels::matfun as mf;
fn exposure(
    a: &Csr,
    t: &Value,
    d: &DimExpr,
) -> Result<f64, OpError> {
    api::same_dimension(&product(d.clone(), api::dimension(t)?)?, &DimExpr::Dimless)?;
    let tau = crate::rules::numerical_abi::numeric(t, 2)? * a.scale;
    if !tau.is_finite() || tau.abs() * mf::csr_norm1(&a.data) > 128.0 {
        return Err(err("matrix exponential exposure/work budget exceeded"));
    }
    Ok(tau)
}
object_binding!(csr_norm1,"num_linear_csr_norm1","crate::kernels::matfun::csr_norm1","Native sparse matrix 1-norm preserving the matrix entry reference.",true,(a:Handle<Csr>)->Dynamic,|ctx,a|{let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(mf::csr_norm1(&a.data)*a.scale,d)?)))});
object_binding!(expm_multiply,"num_linear_expm_multiply","crate::kernels::matfun::expm_multiply","Native sparse exponential action exp(tA)v. Matrix and exposure dimensions must multiply to one; vector units are preserved. Exposure/operator size are bounded before native stepping.",false,(a:Handle<Csr>,v:Handle<Array>,t:Input)->NewObject<Array>,|ctx,a,v,t|{let da=dim(ctx,&a.value)?;let dv=dim(ctx,&v.value)?;a.with(ctx,|a|v.with(ctx,|v|{if a.data.rows!=a.data.cols||a.data.rows>512{return Err(err("exponential action needs bounded square matrix"));}let input=vector(v,a.data.cols)?;result(mf::expm_multiply(&a.data,input,exposure(a,&t.0,&da)?),v.scale,dv)}))});
object_binding!(expm_krylov,"num_linear_expm_krylov","crate::kernels::matfun::expm_krylov","Native matrix-free Krylov exponential action with explicit dimensionless A*t and preserved vector units. None reports native failure to meet its approximation controls; no algorithm substitution occurs.",false,(a:Handle<Csr>,v:Handle<Array>,t:Input,space:usize,tolerance:Normalized<f64>)->Maybe<NewObject<Array>>,|ctx,a,v,t,space,tolerance|{let da=dim(ctx,&a.value)?;let dv=dim(ctx,&v.value)?;a.with(ctx,|a|v.with(ctx,|v|{if a.data.rows!=a.data.cols||a.data.rows>512||space==0||space>64||tolerance.0<=0.0||a.data.values.len().checked_mul(space).is_none_or(|w|w>WORK){return Err(err("Krylov exponential shape/control/work budget exceeded"));}let input=vector(v,a.data.cols)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);Ok(Maybe(mf::expm_krylov(apply,input,exposure(a,&t.0,&da)?,space,tolerance.0).map(|data|result(data,v.scale,dv)).transpose()?))}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &csr_norm1::BINDING,
            &expm_multiply::BINDING,
            &expm_krylov::BINDING,
        ],
        "numerical_matrix",
    )
}
