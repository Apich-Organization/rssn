//! Explicit f64 monomorphizations of native Field backend hooks.
use super::super::native_object_api::Maybe;
use super::Array;
use super::ArrayD;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::IxDyn;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Value;
use super::api;
use super::array;
use super::dim;
use super::err;
use super::object_binding;
use super::packed;
use super::product;
use super::scale;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernels::matrix::FaerDecompositionResult as ResultKind;
use crate::kernels::matrix::FaerDecompositionType as Kind;
use crate::kernels::matrix::Field;
use crate::kernels::matrix::Matrix;
fn matrix(a: &Array) -> Result<Matrix<f64>, OpError> {
    let m = super::lowrank::matrix(a)?;
    Ok(Matrix::new(m.rows, m.cols, m.data))
}
fn restore(
    m: Matrix<f64>,
    reference: f64,
    d: DimExpr,
) -> Result<NewObject<Array>, OpError> {
    Ok(array(
        ArrayD::from_shape_vec(IxDyn(&[m.rows(), m.cols()]), m.into_data()).map_err(err)?,
        reference,
        d,
    ))
}
fn packed_matrix(
    m: Matrix<f64>,
    reference: f64,
    d: DimExpr,
) -> Result<Value, OpError> {
    packed(
        vec![m.rows(), m.cols()],
        m.into_data(),
        api::quantity(reference, d)?,
    )
}
fn inverted(d: DimExpr) -> Result<DimExpr, OpError> {
    crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Inverse(Box::new(
            d,
        )))?,
    )
}
object_binding!(scalar_invertible,"num_linear_field_is_invertible","crate::kernels::matrix::Field::is_invertible","Native f64 Field invertibility of a finite measured scalar.",false,(value:Input)->bool,|ctx,value|{let _=ctx;Ok(<f64 as Field>::is_invertible(&crate::rules::numerical_abi::numeric(&value.0,0)?))});
object_binding!(scalar_inverse,"num_linear_field_inverse","crate::kernels::matrix::Field::inverse","Native f64 Field inverse; output dimension is the inverse of the measured scalar dimension.",false,(value:Input)->Dynamic,|ctx,value|{let _=ctx;let x=crate::rules::numerical_abi::numeric(&value.0,0)?;let d=inverted(api::dimension(&value.0)?)?;Ok(Dynamic(api::quantity(<f64 as Field>::inverse(&x).map_err(err)?,d)?))});
object_binding!(scalar_magnitude,"num_linear_field_magnitude","crate::kernels::matrix::Field::magnitude","Native f64 Field magnitude preserves the scalar's physical dimension.",false,(value:Input)->Dynamic,|ctx,value|{let _=ctx;let x=crate::rules::numerical_abi::numeric(&value.0,0)?;Ok(Dynamic(api::quantity(<f64 as Field>::magnitude(&x),api::dimension(&value.0)?)?))});
object_binding!(faer_mul,"num_linear_field_faer_mul","crate::kernels::matrix::Field::faer_mul","Actual native f64 Field faer multiplication hook; product retains the independent matrix references.",true,(a:Handle<Array>,b:Handle<Array>)->Maybe<NewObject<Array>>,|ctx,a,b|{let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|Ok(Maybe(<f64 as Field>::faer_mul(&matrix(a)?,&matrix(b)?).map(|m|restore(m,scale(a.scale,b.scale)?,d)).transpose()?))))});
object_binding!(faer_inverse,"num_linear_field_faer_inverse","crate::kernels::matrix::Field::faer_inverse","Actual native f64 Field faer inverse hook; absent inverse is None, with physical inverse entry reference.",true,(a:Handle<Array>)->Maybe<NewObject<Array>>,|ctx,a|{let d=inverted(dim(ctx,&a.value)?)?;a.with(ctx,|a|Ok(Maybe(<f64 as Field>::faer_inverse(&matrix(a)?).map(|m|restore(m,1.0/a.scale,d)).transpose()?))) });
object_binding!(faer_solve,"num_linear_field_faer_solve","crate::kernels::matrix::Field::faer_solve","Actual native f64 Field faer multi-RHS solve hook; solution entries retain RHS/matrix reference dimensions.",true,(a:Handle<Array>,b:Handle<Array>)->Maybe<NewObject<Array>>,|ctx,a,b|{let d=product(dim(ctx,&b.value)?,inverted(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|b.with(ctx,|b|Ok(Maybe(<f64 as Field>::faer_solve(&matrix(a)?,&matrix(b)?).map(|m|restore(m,scale(b.scale,1.0/a.scale)?,d)).transpose()?))))});
object_binding!(faer_decompose,"num_linear_field_faer_decompose","crate::kernels::matrix::Field::faer_decompose","Actual f64 Field backend decomposition hook. Every returned factor carries its own normalized values and physical reference; orthogonal factors are dimensionless, Cholesky factors use sqrt(matrix reference), singular/eigenvalues use matrix reference.",true,(a:Handle<Array>,kind:String)->Maybe<Dynamic>,|ctx,a,kind|{
 let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let kind=match kind.as_str(){"Lu"=>Kind::Lu,"Qr"=>Kind::Qr,"Cholesky"=>Kind::Cholesky,"Svd"=>Kind::Svd,"EigenSymmetric"=>Kind::EigenSymmetric,_=>return Err(err("unknown native Field decomposition"))};let Some(out)=<f64 as Field>::faer_decompose(&0.0,&matrix(a)?,kind)else{return Ok(Maybe(None));};let value=match out{ResultKind::Lu{l,u,p}=>api::record(vec![("kind",Value::String("Lu".into())),("l",packed_matrix(l,1.0,DimExpr::Dimless)?),("u",packed_matrix(u,a.scale,d)?),("p",p.encode()?)]),ResultKind::Qr{q,r}=>api::record(vec![("kind",Value::String("Qr".into())),("q",packed_matrix(q,1.0,DimExpr::Dimless)?),("r",packed_matrix(r,a.scale,d)?)]),ResultKind::Cholesky{l}=>{let d=crate::rules::numerical_abi::runtime_dimensions::expression(&crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Pow(Box::new(d),Exp::Rat(Rational::new(1,2))))?)?;api::record(vec![("kind",Value::String("Cholesky".into())),("l",packed_matrix(l,a.scale.sqrt(),d)?)])},ResultKind::Svd{u,s,v}=>api::record(vec![("kind",Value::String("Svd".into())),("u",packed_matrix(u,1.0,DimExpr::Dimless)?),("s",packed(vec![s.len()],s,api::quantity(a.scale,d)?)?),("v",packed_matrix(v,1.0,DimExpr::Dimless)?)]),ResultKind::EigenSymmetric{values,vectors}=>api::record(vec![("kind",Value::String("EigenSymmetric".into())),("values",packed(vec![values.len()],values,api::quantity(a.scale,d)?)?),("vectors",packed_matrix(vectors,1.0,DimExpr::Dimless)?)])};Ok(Maybe(Some(Dynamic(value))))})
});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &scalar_invertible::BINDING,
            &scalar_inverse::BINDING,
            &scalar_magnitude::BINDING,
            &faer_mul::BINDING,
            &faer_inverse::BINDING,
            &faer_solve::BINDING,
            &faer_decompose::BINDING,
        ],
        "numerical_matrix",
    )
}
