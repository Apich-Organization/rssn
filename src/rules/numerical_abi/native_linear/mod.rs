//! Checked native linear and tensor objects with explicit physical reference scales.
//! Every stored tensor is normalized; its positive scale and unit domain are retained
//! with the session capability. Factors and compressed cores are dimensionless.
use super::codecs::Arg;
use super::codecs::Output;
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernels::tensor;
use crate::kernels::tensor_train;
use ndarray::ArrayD;
use ndarray::IxDyn;
const LIMIT: usize = 262_144;
const WORK: usize = 8_388_608;
fn err(e: impl std::fmt::Display) -> OpError {
    api::error(e)
}
fn count(shape: &[usize]) -> Result<usize, OpError> {
    if shape.len() > 8 || shape.iter().any(|n| *n > LIMIT) {
        return Err(err("native tensor rank/extent budget exceeded"));
    }
    let n = shape
        .iter()
        .try_fold(1usize, |n, m| n.checked_mul(*m))
        .ok_or_else(|| err("native tensor extent overflow"))?;
    if n > LIMIT {
        return Err(err("native tensor allocation budget exceeded"));
    }
    Ok(n)
}
fn dim(
    ctx: &CallCtx<'_>,
    h: &Value,
) -> Result<DimExpr, OpError> {
    api::dimensions(ctx, h)?
        .into_iter()
        .next()
        .ok_or_else(|| err("native tensor dimension missing"))
}
fn product(
    a: DimExpr,
    b: DimExpr,
) -> Result<DimExpr, OpError> {
    crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Mul(vec![a, b]))?,
    )
}
fn scale(
    a: f64,
    b: f64,
) -> Result<f64, OpError> {
    let n = a * b;
    if !n.is_finite() || n <= 0.0 {
        return Err(err("native tensor reference scale overflow/underflow"));
    }
    Ok(n)
}
struct Input(Value);
impl NativeArg for Input {
    fn ty() -> TypeExpr {
        api::dynamic_type()
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::bounded(v)?;
        Ok(Self(v.clone()))
    }
}
pub(super) struct Normalized<T>(pub T);
impl<T: Arg> NativeArg for Normalized<T> {
    fn ty() -> TypeExpr {
        T::ty()
    }

    fn decode(
        v: &Value,
        i: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(v)?;
        Ok(Self(T::decode(v, i)?))
    }
}
#[derive(Clone)]
struct Array {
    data: ArrayD<f64>,
    scale: f64,
}
impl NativeType for Array {
    const KIND: &'static str = "rssn/native/array/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 512;

    fn validate(&self) -> Result<(), OpError> {
        count(self.data.shape())?;
        if !self.scale.is_finite() || self.scale <= 0.0 || self.data.iter().any(|v| !v.is_finite())
        {
            return Err(err(
                "invalid finite normalized tensor or positive reference scale",
            ));
        }
        Ok(())
    }
}
fn array(
    data: ArrayD<f64>,
    scale: f64,
    d: DimExpr,
) -> NewObject<Array> {
    NewObject::new(Array { data, scale }, vec![d])
}
fn packed(
    shape: Vec<usize>,
    data: Vec<f64>,
    reference: Value,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("shape", shape.encode()?),
        ("data", data.encode()?),
        ("reference", reference),
    ]))
}
fn check_work(
    shape: &[usize],
    rank: usize,
    steps: usize,
) -> Result<(), OpError> {
    let n = count(shape)?;
    if rank == 0
        || rank > 64
        || steps == 0
        || steps > 2048
        || n.checked_mul(rank)
            .and_then(|x| x.checked_mul(steps))
            .is_none_or(|n| n > WORK)
    {
        return Err(err("native tensor decomposition work budget exceeded"));
    }
    Ok(())
}
struct Train {
    data: tensor_train::TensorTrain,
    scale: f64,
}
impl NativeType for Train {
    const KIND: &'static str = "rssn/native/tensor-train/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.shape().len() > 8
            || self.data.storage_len() > LIMIT
            || self.data.shape().iter().any(|n| *n > LIMIT)
            || !self.scale.is_finite()
            || self.scale <= 0.0
        {
            return Err(err("tensor train exceeds retained size/rank/scale budget"));
        }
        Ok(())
    }
}
impl NativeType for tensor_train::TrainCore {
    const KIND: &'static str = "rssn/native/train-core/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 256;

    fn validate(&self) -> Result<(), OpError> {
        count(&self.shape())?;
        Ok(())
    }
}
fn core_budget(c: &tensor_train::TrainCore) -> Result<(), OpError> {
    count(&c.shape())?;
    if c.shape()[0] > 64 || c.shape()[2] > 64 {
        return Err(err("tensor train bond rank budget exceeded"));
    }
    Ok(())
}
fn train_result(
    data: tensor_train::TensorTrain,
    scale: f64,
    d: DimExpr,
) -> NewObject<Train> {
    NewObject::new(Train { data, scale }, vec![d])
}
object_binding!(linear_array_new,"num_linear_array_new","crate::kernels::tensor::TensorData::to_arrayd","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",false,(shape:Vec<usize>,data:Normalized<Vec<f64>>,reference:Input)->NewObject<Array>,|ctx,shape,data,reference|{let _=ctx;if count(&shape)?!=data.0.len(){return Err(err("tensor shape/storage mismatch"));}let (s,d)=api::reference(&reference.0)?;let data=tensor::TensorData{shape,data:data.0}.to_arrayd().map_err(err)?;Ok(array(data,s,d))});
object_binding!(linear_array_data,"num_linear_array_data","","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(packed(a.data.shape().to_vec(),a.data.iter().copied().collect(),api::quantity(a.scale,d)?)?)))});
object_binding!(linear_array_norm,"num_linear_array_norm","crate::kernels::tensor::norm","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(tensor::norm(&a.data)*a.scale,d)?)))});
object_binding!(linear_outer,"num_linear_outer","crate::kernels::tensor::outer_product","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,b:Handle<Array>)->NewObject<Array>,|ctx,a,b|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|{if a.data.len().checked_mul(b.data.len()).is_none_or(|n|n>LIMIT){return Err(err("outer product allocation budget exceeded"));}Ok(array(tensor::outer_product(&a.data,&b.data).map_err(err)?,scale(a.scale,b.scale)?,d))}))});
object_binding!(linear_inner,"num_linear_inner","crate::kernels::tensor::inner_product","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,b:Handle<Array>)->Dynamic,|ctx,a,b|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|{Ok(Dynamic(api::quantity(tensor::inner_product(&a.data,&b.data).map_err(err)?*scale(a.scale,b.scale)?,d)?))}))});
object_binding!(linear_tensordot,"num_linear_tensordot","crate::kernels::tensor::tensordot","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,b:Handle<Array>,axes_a:Vec<usize>,axes_b:Vec<usize>)->NewObject<Array>,|ctx,a,b,axes_a,axes_b|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|{let aa=(0..a.data.ndim()).filter(|i|!axes_a.contains(i)).map(|i|a.data.shape()[i]);let bb=(0..b.data.ndim()).filter(|i|!axes_b.contains(i)).map(|i|b.data.shape()[i]);count(&aa.chain(bb).collect::<Vec<_>>())?;Ok(array(tensor::tensordot(&a.data,&b.data,&axes_a,&axes_b).map_err(err)?,scale(a.scale,b.scale)?,d))}))});
object_binding!(linear_contract,"num_linear_contract","crate::kernels::tensor::contract","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,axis1:usize,axis2:usize)->NewObject<Array>,|ctx,a,axis1,axis2|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(tensor::contract(&a.data,axis1,axis2).map_err(err)?,a.scale,d)))});
object_binding!(linear_tensor_vec_mul,"num_linear_tensor_vec_mul","crate::kernels::tensor::tensor_vec_mul","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,b:Handle<Array>)->NewObject<Array>,|ctx,a,b|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|{if b.data.ndim()!=1{return Err(err("tensor-vector multiplication requires rank-one vector"));}Ok(array(tensor::tensor_vec_mul(&a.data,b.data.as_slice().ok_or_else(||err("vector is not contiguous"))?).map_err(err)?,scale(a.scale,b.scale)?,d))}))});
object_binding!(linear_tt_svd,"num_linear_tt_svd","crate::kernels::tensor_train::tt_svd","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Array>,tolerance:f64,max_rank:usize)->NewObject<Train>,|ctx,a,tolerance,max_rank|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{check_work(a.data.shape(),max_rank,1)?;let data=a.data.as_slice().ok_or_else(||err("tensor not contiguous"))?;Ok(train_result(tensor_train::tt_svd(data,a.data.shape(),tolerance,Some(max_rank)).map_err(err)?,a.scale,d))})});
object_binding!(linear_train_core_new,"num_linear_train_core_new","crate::kernels::tensor_train::TrainCore::new","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(shape:Vec<usize>,data:Vec<f64>)->NewObject<tensor_train::TrainCore>,|ctx,shape,data|{let _=ctx;if shape.len()!=3||count(&shape)?!=data.len(){return Err(err("train core shape/storage mismatch"));}let shape:[usize;3]=shape.try_into().map_err(|_|err("train core needs three axes"))?;let core=tensor_train::TrainCore::new(shape,data).map_err(err)?;core_budget(&core)?;Ok(NewObject::normalized(core))});
object_binding!(linear_train_core_shape,"num_linear_train_core_shape","crate::kernels::tensor_train::TrainCore::shape","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(c:Handle<tensor_train::TrainCore>)->Vec<usize>,|ctx,c|{let _=ctx;c.with(ctx,|c|Ok(c.shape().to_vec()))});
object_binding!(linear_train_core_data,"num_linear_train_core_data","crate::kernels::tensor_train::TrainCore::data","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(c:Handle<tensor_train::TrainCore>)->Vec<f64>,|ctx,c|{let _=ctx;c.with(ctx,|c|Ok(c.data().to_vec()))});
object_binding!(linear_train_from_cores,"num_linear_train_from_cores","crate::kernels::tensor_train::TensorTrain::from_cores","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",false,(cores:Input,reference:Input)->NewObject<Train>,|ctx,cores,reference|{let _=ctx;let Value::Tensor(t)=&cores.0 else{return Err(err("expected train core capability vector"));};
if t.shape!=[t.data.len()]||t.data.is_empty()||t.data.len()>8{return Err(err("tensor train order must be one through eight"));}let cores=t.data.iter().map(|v|api::get::<tensor_train::TrainCore,_>(ctx,v,tensor_train::TrainCore::KIND,|c|{core_budget(c)?;Ok(c.clone())})).collect::<Result<Vec<_>,_>>()?;let(s,d)=api::reference(&reference.0)?;Ok(train_result(tensor_train::TensorTrain::from_cores(cores).map_err(err)?,s,d))});
object_binding!(linear_train_cores,"num_linear_train_cores","crate::kernels::tensor_train::TensorTrain::cores","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let data=a.data.cores().iter().map(|c|api::store(ctx,tensor_train::TrainCore::KIND,vec![],c.clone(),tensor_train::TrainCore::MAX_BYTES)).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![data.len()],elem:TypeExpr::Opaque(tensor_train::TrainCore::KIND.into()),data})))})});
object_binding!(linear_train_shape,"num_linear_train_shape","crate::kernels::tensor_train::TensorTrain::shape","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.shape().to_vec()))});
object_binding!(linear_train_storage_len,"num_linear_train_storage_len","crate::kernels::tensor_train::TensorTrain::storage_len","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.storage_len()))});
object_binding!(linear_train_tolerance_met,"num_linear_train_tolerance_met","crate::kernels::tensor_train::TensorTrain::tolerance_met","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.tolerance_met()))});
object_binding!(linear_train_norm,"num_linear_train_norm","crate::kernels::tensor_train::TensorTrain::frobenius_norm","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(a.data.frobenius_norm()*a.scale,d)?)))});
object_binding!(linear_train_truncation_error,"num_linear_train_truncation_error","crate::kernels::tensor_train::TensorTrain::truncation_error","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(a.data.truncation_error()*a.scale,d)?)))});
object_binding!(linear_train_entry,"num_linear_train_entry","crate::kernels::tensor_train::TensorTrain::entry","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,indices:Vec<usize>)->Dynamic,|ctx,a,indices|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(a.data.entry(&indices).map_err(err)?*a.scale,d)?)))});
object_binding!(linear_train_to_dense,"num_linear_train_to_dense","crate::kernels::tensor_train::TensorTrain::to_dense","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{count(a.data.shape())?;let data=a.data.to_dense().map_err(err)?;Ok(array(ArrayD::from_shape_vec(IxDyn(a.data.shape()),data).map_err(err)?,a.scale,d))})});
object_binding!(linear_train_to_dense_with_limit,"num_linear_train_to_dense_with_limit","crate::kernels::tensor_train::TensorTrain::to_dense_with_limit","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,max_elements:usize)->NewObject<Array>,|ctx,a,max_elements|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{count(a.data.shape())?;if max_elements>LIMIT{return Err(err("dense limit exceeds allocation budget"));}let data=a.data.to_dense_with_limit(max_elements).map_err(err)?;Ok(array(ArrayD::from_shape_vec(IxDyn(a.data.shape()),data).map_err(err)?,a.scale,d))})});
object_binding!(linear_train_add,"num_linear_train_add","crate::kernels::tensor_train::TensorTrain::add","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,b:Handle<Train>)->NewObject<Train>,|ctx,a,b|{let _=ctx;let da=dim(ctx,&a.value)?;let db=dim(ctx,&b.value)?;api::same_dimension(&da,&db)?;a.with(ctx,|a|b.with(ctx,|b|{if a.data.shape()!=b.data.shape(){return Err(err("tensor train shapes differ"));}
if a.data.storage_len().checked_mul(b.data.storage_len()).is_none_or(|n|n>WORK){return Err(err("tensor train operation work budget exceeded"));}
if a.scale.to_bits()!=b.scale.to_bits(){return Err(err("addition needs equal numeric reference scales; explicitly renormalize first"));}Ok(train_result(a.data.add(&b.data).map_err(err)?,a.scale,da))}))});
object_binding!(linear_train_hadamard,"num_linear_train_hadamard","crate::kernels::tensor_train::TensorTrain::hadamard","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,b:Handle<Train>)->NewObject<Train>,|ctx,a,b|{let _=ctx;let da=dim(ctx,&a.value)?;let db=dim(ctx,&b.value)?;a.with(ctx,|a|b.with(ctx,|b|{if a.data.shape()!=b.data.shape(){return Err(err("tensor train shapes differ"));}
if a.data.storage_len().checked_mul(b.data.storage_len()).is_none_or(|n|n>WORK){return Err(err("tensor train operation work budget exceeded"));}Ok(train_result(a.data.hadamard(&b.data).map_err(err)?,scale(a.scale,b.scale)?,product(da,db)?))}))});
object_binding!(linear_train_inner,"num_linear_train_inner","crate::kernels::tensor_train::TensorTrain::inner","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,b:Handle<Train>)->Dynamic,|ctx,a,b|{let _=ctx;let da=dim(ctx,&a.value)?;let db=dim(ctx,&b.value)?;a.with(ctx,|a|b.with(ctx,|b|{if a.data.shape()!=b.data.shape(){return Err(err("tensor train shapes differ"));}
if a.data.storage_len().checked_mul(b.data.storage_len()).is_none_or(|n|n>WORK){return Err(err("tensor train operation work budget exceeded"));}let value=a.data.inner(&b.data).map_err(err)?;Ok(Dynamic(api::quantity(api::scale_product(value,&[a.scale,b.scale],&[])?,product(da,db)?)?))}))});
object_binding!(
    linear_train_inner_workspace_len,
    "num_linear_train_inner_workspace_len",
    "crate::kernels::tensor_train::TensorTrain::inner_workspace_len",
    "Reports the number of Float work slots for contracting two equal-shape tensor trains. Counts depend on bond ranks. A single core or two rank-one trains require no work slots.",
    true,
    (a: Handle<Train>, b: Handle<Train>) -> usize,
    |ctx, a, b| a.with(ctx, |a| b.with(ctx, |b| a.data.inner_workspace_len(&b.data).map_err(err)))
);
object_binding!(linear_train_round,"num_linear_train_round","crate::kernels::tensor_train::TensorTrain::round","Native numerical method. Normalized coordinates and positive physical reference scales are checked; native object kind/session/lifetime and allocation/work budgets are enforced.",true,(a:Handle<Train>,tolerance:f64,max_rank:usize)->NewObject<Train>,|ctx,a,tolerance,max_rank|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if max_rank==0||max_rank>64{return Err(err("tensor train rank budget exceeded"));}Ok(train_result(a.data.round(tolerance,Some(max_rank)).map_err(err)?,a.scale,d))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")?;
    decompositions::install(i)?;
    sparse::install(i)?;
    exact::install(i)?;
    lowrank::install(i)?;
    interpolation::install(i)?;
    contractions::install(i)?;
    schur::install(i)?;
    multigrid::install(i)?;
    eigen::install(i)?;
    csc::install(i)?;
    matrix_field::install(i)?;
    matfun::install(i)
}
static BINDINGS: &[&api::Binding] = &[
    &linear_solve_linear_system::BINDING,
    &linear_array_new::BINDING,
    &linear_array_data::BINDING,
    &linear_array_norm::BINDING,
    &linear_outer::BINDING,
    &linear_inner::BINDING,
    &linear_tensordot::BINDING,
    &linear_contract::BINDING,
    &linear_tensor_vec_mul::BINDING,
    &linear_tt_svd::BINDING,
    &linear_train_core_new::BINDING,
    &linear_train_core_shape::BINDING,
    &linear_train_core_data::BINDING,
    &linear_train_from_cores::BINDING,
    &linear_train_cores::BINDING,
    &linear_train_shape::BINDING,
    &linear_train_storage_len::BINDING,
    &linear_train_tolerance_met::BINDING,
    &linear_train_norm::BINDING,
    &linear_train_truncation_error::BINDING,
    &linear_train_entry::BINDING,
    &linear_train_to_dense::BINDING,
    &linear_train_to_dense_with_limit::BINDING,
    &linear_train_add::BINDING,
    &linear_train_hadamard::BINDING,
    &linear_train_inner::BINDING,
    &linear_train_inner_workspace_len::BINDING,
    &linear_train_round::BINDING,
];

mod decompositions;

mod sparse;

mod exact;

object_binding!(linear_solve_linear_system,"num_linear_solve_linear_system","crate::kernels::solve::solve_linear_system","Native RREF linear solve for normalized matrix/RHS arrays with independent physical references. Particular solutions carry RHS/matrix units; nullspace coefficients are dimensionless and scale the free physical parameters. Reports Unique, Parametric or NoSolution.",true,(a:Handle<Array>,rhs:Handle<Array>)->Dynamic,|ctx,a,rhs|{
 let d=product(dim(ctx,&rhs.value)?,crate::rules::numerical_abi::runtime_dimensions::expression(&crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Inverse(Box::new(dim(ctx,&a.value)?)))?)?)?;
 a.with(ctx,|a|rhs.with(ctx,|b|{if a.data.ndim()!=2||a.data.shape().iter().any(|n|*n>128)||b.data.ndim()!=1||b.data.len()!=a.data.shape()[0]{return Err(err("linear system shape/work budget exceeded"));}let m=crate::kernels::matrix::Matrix::new(a.data.shape()[0],a.data.shape()[1],a.data.iter().copied().collect());let solution=crate::kernels::solve::solve_linear_system(&m,b.data.as_slice().ok_or_else(||err("RHS not contiguous"))?).map_err(err)?;let (kind,particular,basis_shape,basis)=match solution{crate::kernels::solve::LinearSolution::Unique(v)=>("Unique",v,vec![a.data.shape()[1],0],vec![]),crate::kernels::solve::LinearSolution::Parametric{particular,null_space_basis}=>("Parametric",particular,vec![null_space_basis.rows(),null_space_basis.cols()],null_space_basis.into_data()),crate::kernels::solve::LinearSolution::NoSolution=>("NoSolution",vec![],vec![a.data.shape()[1],0],vec![])};Ok(Dynamic(api::record(vec![("kind",Value::String(kind.into())),("particular",particular.encode()?),("solution_reference",api::quantity(scale(b.scale,1.0/a.scale)?,d)?),("basis_shape",basis_shape.encode()?),("basis",basis.encode()?)])))}))
});

mod lowrank;

mod interpolation;

mod contractions;

mod schur;

mod multigrid;

mod eigen;

mod csc;

mod matfun;
mod matrix_field;
