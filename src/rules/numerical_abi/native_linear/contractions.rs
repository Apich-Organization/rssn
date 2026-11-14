//! Genuine scalar callback monomorphizations for general contraction and simplex integration.
use super::super::native_backend::compiled_snapshot;
use super::Array;
use super::ArrayD;
use super::DimExpr;
use super::Dynamic;
use super::Input;
use super::Installer;
use super::IxDyn;
use super::LIMIT;
use super::NewObject;
use super::Normalized;
use super::Output;
use super::RuleError;
use super::api;
use super::array;
use super::err;
use super::object_binding;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
object_binding!(rectangular,"num_linear_contract_rectangular","crate::kernels::contractions::contract_rectangular","Native generic contraction specialized to f64. The compiled scalar callback receives row, inner index, column, has_previous, previous in that order; index/presence inputs are dimensionless and previous/output use the explicit reference dimension. Indices are zero based and coherent physical values are normalized only at the callback boundary.",false,(rows:usize,inner:usize,cols:usize,zero:Normalized<f64>,program:Input,reference:Input)->NewObject<Array>,|ctx,rows,inner,cols,zero,program,reference|{
 if rows>512||inner>512||cols>512||rows.checked_mul(cols).is_none_or(|n|n>LIMIT)||rows.checked_mul(inner).and_then(|n|n.checked_mul(cols)).is_none_or(|n|n>1_000_000){return Err(err("rectangular contraction work/allocation budget exceeded"));}
 let (scale,d)=api::reference(&reference.0)?;let budget=rows*inner*cols;let callback=compiled_snapshot(ctx,&program.0,budget.max(1))?;if callback.input_dimensions.len()!=5{return Err(err("contraction callback needs five inputs"));}for unit in &callback.input_dimensions[..4]{api::same_dimension(unit,&DimExpr::Dimless)?;}api::same_dimension(&callback.input_dimensions[4],&d)?;api::same_dimension(&callback.output_dimension,&d)?;
 let values=crate::kernels::contractions::contract_rectangular(rows,inner,cols,||zero.0,|row,index,column,previous:Option<f64>|{callback.evaluate_coherent(&[row as f64,index as f64,column as f64,if previous.is_some(){1.0}else{0.0},previous.unwrap_or(0.0)*scale]).map(|y|y/scale)}).map_err(|e|err(format!("{e:?}")))?;
 Ok(array(ArrayD::from_shape_vec(IxDyn(&[rows,cols]),values).map_err(err)?,scale,d))
});
object_binding!(simplex,"num_linear_integrate_simplex","crate::kernels::simplex_quadrature::integrate_simplex","Native adaptive simplex integration of a retained compiled scalar field, with explicit independent coordinate and ordinate references. Each vertex is normalized by the coordinate reference; value/error dimensions are ordinate times coordinate^dimension. Native evaluation and dimension limits are preserved with a checked finite callback/work budget.",false,(vertices:Normalized<Vec<Vec<f64>>>,coordinate_reference:Input,value_reference:Input,program:Input,absolute_tolerance:Normalized<f64>,relative_tolerance:Normalized<f64>,max_evaluations:usize)->Dynamic,|ctx,vertices,coordinate_reference,value_reference,program,absolute_tolerance,relative_tolerance,max_evaluations|{
 let n=vertices.0.len().checked_sub(1).ok_or_else(||err("simplex needs vertices"))?;if !(1..=20).contains(&n)||vertices.0.iter().any(|v|v.len()!=n)||max_evaluations==0||max_evaluations>10_000_000{return Err(err("simplex shape/evaluation budget invalid"));}
 let (x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;let callback=compiled_snapshot(ctx,&program.0,max_evaluations)?;if callback.input_dimensions.len()!=n{return Err(err("simplex callback arity differs from coordinate dimension"));}for input in &callback.input_dimensions{api::same_dimension(input,&dx)?;}api::same_dimension(&callback.output_dimension,&dy)?;
 let function=|p:&[f64]|{let coherent=p.iter().map(|v|v*x).collect::<Vec<_>>();callback.sample(&coherent)/y};
 let out=crate::kernels::simplex_quadrature::integrate_simplex(function,&vertices.0,absolute_tolerance.0,relative_tolerance.0,max_evaluations);callback.finish()?;let out=out.map_err(|e|err(format!("{e:?}")))?;let dim=crate::rules::numerical_abi::runtime_dimensions::expression(&crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Mul(vec![dy,DimExpr::Pow(Box::new(dx),Exp::Rat(Rational::new(n as i64,1)))]))?)?;let scale=y*x.powi(n as i32);if !scale.is_finite()||scale<=0.0{return Err(err("simplex reference volume overflow/underflow"));}
 Ok(Dynamic(api::record(vec![("value",api::quantity(out.value*scale,dim.clone())?),("error",api::quantity(out.error*scale,dim)?),("evaluations",out.evaluations.encode()?),("converged",out.converged.encode()?)])))
});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[&rectangular::BINDING, &simplex::BINDING],
        "numerical_matrix",
    )
}
