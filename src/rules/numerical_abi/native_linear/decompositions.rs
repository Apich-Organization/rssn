//! Native CP/Tucker/HOSVD construction and reconstruction.
use super::Array;
use super::ArrayD;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::IxDyn;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::array;
use super::check_work;
use super::count;
use super::dim;
use super::err;
use super::object_binding;
use super::tensor;
use crate::kernels::dense::Mat;
impl NativeType for tensor::CpOptions {
    const KIND: &'static str = "rssn/native/cp-options/v1";
    const MAX_BYTES: usize = 256;

    fn validate(&self) -> Result<(), OpError> {
        if self.rank == 0
            || self.rank > 64
            || self.max_sweeps == 0
            || self.max_sweeps > 2048
            || !self.tol.is_finite()
            || self.tol < 0.0
            || !self.divergence_limit.is_finite()
            || self.divergence_limit <= 0.0
        {
            return Err(err("invalid CP options or work budget"));
        }
        Ok(())
    }
}
struct Cp {
    data: tensor::CpDecomposition,
    scale: f64,
}
struct Hosvd {
    data: tensor::Hosvd,
    scale: f64,
}
struct Tucker {
    data: tensor::Tucker,
    scale: f64,
}
fn factor_sizes(factors: &[Mat]) -> Result<usize, OpError> {
    if factors.len() > 8 {
        return Err(err("tensor decomposition order budget exceeded"));
    }
    let n = factors
        .iter()
        .try_fold(0usize, |n, m| n.checked_add(m.data.len()))
        .ok_or_else(|| err("factor size overflow"))?;
    if n > LIMIT {
        return Err(err("tensor decomposition factor storage budget exceeded"));
    }
    Ok(n)
}
impl NativeType for Cp {
    const KIND: &'static str = "rssn/native/cp/v1";
    const MAX_BYTES: usize = LIMIT * 24 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        factor_sizes(&self.data.factors)?;
        count(self.data.reconstruction.shape())?;
        Ok(())
    }
}
impl NativeType for Hosvd {
    const KIND: &'static str = "rssn/native/hosvd/v1";
    const MAX_BYTES: usize = LIMIT * 24 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        factor_sizes(&self.data.factors)?;
        count(self.data.core.shape())?;
        Ok(())
    }
}
impl NativeType for Tucker {
    const KIND: &'static str = "rssn/native/tucker/v1";
    const MAX_BYTES: usize = LIMIT * 24 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        factor_sizes(&self.data.factors)?;
        count(self.data.core.shape())?;
        Ok(())
    }
}
fn matrices(
    ctx: &CallCtx<'_>,
    input: &Input,
) -> Result<Vec<Mat>, OpError> {
    let Value::Tensor(t) = &input.0 else {
        return Err(err("expected factor handle vector"));
    };
    if t.shape != [t.data.len()] || t.data.is_empty() || t.data.len() > 8 {
        return Err(err("factor vector requires 1..=8 handles"));
    }
    let mut out = Vec::new();
    for h in &t.data {
        let d = dim(ctx, h)?;
        api::same_dimension(&d, &DimExpr::Dimless)?;
        out.push(api::get::<Array, _>(ctx, h, Array::KIND, |a| {
            if a.data.ndim() != 2 || a.scale.to_bits() != 1.0f64.to_bits() {
                return Err(err(
                    "factor must be dimensionless matrix with unit numeric reference",
                ));
            }
            Ok(Mat {
                rows: a.data.shape()[0],
                cols: a.data.shape()[1],
                data: a.data.iter().copied().collect(),
            })
        })?);
    }
    factor_sizes(&out)?;
    Ok(out)
}
fn factors(
    ctx: &CallCtx<'_>,
    data: &[Mat],
) -> Result<Dynamic, OpError> {
    factor_sizes(data)?;
    let out = data
        .iter()
        .map(|f| {
            api::store(
                ctx,
                Array::KIND,
                vec![DimExpr::Dimless],
                Array {
                    data: ArrayD::from_shape_vec(IxDyn(&[f.rows, f.cols]), f.data.clone())
                        .map_err(err)?,
                    scale: 1.0,
                },
                Array::MAX_BYTES,
            )
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![out.len()],
        elem: TypeExpr::Opaque(Array::KIND.into()),
        data: out,
    })))
}
object_binding!(decomp_cp_options_new,"num_linear_cp_options_new","crate::kernels::tensor::CpOptions::new","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(rank:usize)->NewObject<tensor::CpOptions>,|ctx,rank|{let _=ctx;if rank==0||rank>64{return Err(err("CP rank budget exceeded"));}Ok(NewObject::normalized(tensor::CpOptions::new(rank)))});
object_binding!(decomp_cp_options,"num_linear_cp_options","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(options:Handle<tensor::CpOptions>,sweeps:usize,tolerance:f64,seed:u64,divergence_limit:f64)->NewObject<tensor::CpOptions>,|ctx,options,sweeps,tolerance,seed,divergence_limit|{let _=ctx;options.with(ctx,|o|{let mut o=o.clone();o.max_sweeps=sweeps;o.tol=tolerance;o.seed=seed;o.divergence_limit=divergence_limit;o.validate()?;Ok(NewObject::normalized(o))})});
object_binding!(decomp_cp_als,"num_linear_cp_als","crate::kernels::tensor::cp_als","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Array>,options:Handle<tensor::CpOptions>)->NewObject<Cp>,|ctx,a,options|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|options.with(ctx,|o|{check_work(a.data.shape(),o.rank,o.max_sweeps)?;Ok(NewObject::new(Cp{data:tensor::cp_als(&a.data,o).map_err(err)?,scale:a.scale},vec![d]))}))});
object_binding!(decomp_hosvd,"num_linear_hosvd","crate::kernels::tensor::hosvd","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Array>,ranks:Vec<usize>)->NewObject<Hosvd>,|ctx,a,ranks|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{check_work(a.data.shape(),ranks.iter().copied().max().unwrap_or(1),1)?;Ok(NewObject::new(Hosvd{data:tensor::hosvd(&a.data,if ranks.is_empty(){None}else{Some(&ranks)}).map_err(err)?,scale:a.scale},vec![d]))})});
object_binding!(decomp_hooi,"num_linear_hooi","crate::kernels::tensor::hooi","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Array>,ranks:Vec<usize>,sweeps:usize,tolerance:f64)->NewObject<Tucker>,|ctx,a,ranks,sweeps,tolerance|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{check_work(a.data.shape(),ranks.iter().copied().max().unwrap_or(1),sweeps)?;Ok(NewObject::new(Tucker{data:tensor::hooi(&a.data,&ranks,sweeps,tolerance).map_err(err)?,scale:a.scale},vec![d]))})});
object_binding!(decomp_hosvd_reconstruct,"num_linear_hosvd_reconstruct","crate::kernels::tensor::Hosvd::reconstruct","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Hosvd>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let shape=a.data.factors.iter().map(|f|f.rows).collect::<Vec<_>>();count(&shape)?;Ok(array(a.data.reconstruct().map_err(err)?,a.scale,d))})});
object_binding!(decomp_tucker_reconstruct_method,"num_linear_tucker_reconstruct_method","crate::kernels::tensor::Tucker::reconstruct","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Tucker>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let shape=a.data.factors.iter().map(|f|f.rows).collect::<Vec<_>>();count(&shape)?;Ok(array(a.data.reconstruct().map_err(err)?,a.scale,d))})});
object_binding!(decomp_cp_factors,"num_linear_cp_factors","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Cp>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|factors(ctx,&a.data.factors))});
object_binding!(decomp_hosvd_factors,"num_linear_hosvd_factors","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Hosvd>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|factors(ctx,&a.data.factors))});
object_binding!(decomp_tucker_factors,"num_linear_tucker_factors","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Tucker>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|factors(ctx,&a.data.factors))});
object_binding!(decomp_cp_reconstruction,"num_linear_cp_reconstruction","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Cp>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(a.data.reconstruction.clone(),a.scale,d)))});
object_binding!(decomp_cp_weights,"num_linear_cp_weights","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Cp>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(ArrayD::from_shape_vec(IxDyn(&[a.data.weights.len()]),a.data.weights.clone()).map_err(err)?,a.scale,d)))});
object_binding!(decomp_hosvd_core,"num_linear_hosvd_core","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Hosvd>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(a.data.core.clone(),a.scale,d)))});
object_binding!(decomp_tucker_core,"num_linear_tucker_core","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Tucker>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(a.data.core.clone(),a.scale,d)))});
object_binding!(decomp_cp_diagnostics,"num_linear_cp_diagnostics","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Cp>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::record(vec![("residual",api::quantity(a.data.residual*a.scale,d)?),("relative_residual",a.data.relative_residual.encode()?),("sweeps",a.data.sweeps.encode()?),("converged",a.data.converged.encode()?),("stagnated",a.data.stagnated.encode()?),("diverged",a.data.diverged.encode()?),("rank_deficient",a.data.rank_deficient.encode()?)]))))});
object_binding!(decomp_tucker_diagnostics,"num_linear_tucker_diagnostics","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Tucker>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(api::record(vec![("ranks",a.data.ranks.clone().encode()?),("sweeps",a.data.sweeps.encode()?),("relative_error",a.data.relative_error.encode()?),("converged",a.data.converged.encode()?),("error_history",a.data.error_history.clone().encode()?)]))))});
object_binding!(decomp_hosvd_ranks,"num_linear_hosvd_ranks","","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(a:Handle<Hosvd>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.ranks.clone()))});
object_binding!(decomp_tucker_reconstruct,"num_linear_tucker_reconstruct","crate::kernels::tensor::tucker_reconstruct","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(core:Handle<Array>,factors:Input)->NewObject<Array>,|ctx,core,factors|{let _=ctx;let d=dim(ctx,&core.value)?;let factors=matrices(ctx,&factors)?;count(&factors.iter().map(|f|f.rows).collect::<Vec<_>>())?;core.with(ctx,|a|Ok(array(tensor::tucker_reconstruct(&a.data,&factors).map_err(err)?,a.scale,d)))});
object_binding!(decomp_cp_reconstruct,"num_linear_cp_reconstruct","crate::kernels::tensor::cp_reconstruct","Native tensor decomposition with explicit reference scale and bounded factor/order/sweep storage. Factors are normalized; reconstructed tensors retain the input physical domain.",true,(weights:Handle<Array>,factors:Input)->NewObject<Array>,|ctx,weights,factors|{let _=ctx;let d=dim(ctx,&weights.value)?;let factors=matrices(ctx,&factors)?;count(&factors.iter().map(|f|f.rows).collect::<Vec<_>>())?;weights.with(ctx,|a|{if a.data.ndim()!=1{return Err(err("CP weights must have rank one"));}Ok(array(tensor::cp_reconstruct(a.data.as_slice().ok_or_else(||err("CP weights not contiguous"))?,&factors).map_err(err)?,a.scale,d))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &decomp_cp_options_new::BINDING,
    &decomp_cp_options::BINDING,
    &decomp_cp_als::BINDING,
    &decomp_hosvd::BINDING,
    &decomp_hooi::BINDING,
    &decomp_hosvd_reconstruct::BINDING,
    &decomp_tucker_reconstruct_method::BINDING,
    &decomp_cp_factors::BINDING,
    &decomp_hosvd_factors::BINDING,
    &decomp_tucker_factors::BINDING,
    &decomp_cp_reconstruction::BINDING,
    &decomp_cp_weights::BINDING,
    &decomp_hosvd_core::BINDING,
    &decomp_tucker_core::BINDING,
    &decomp_cp_diagnostics::BINDING,
    &decomp_tucker_diagnostics::BINDING,
    &decomp_hosvd_ranks::BINDING,
    &decomp_tucker_reconstruct::BINDING,
    &decomp_cp_reconstruct::BINDING,
];
