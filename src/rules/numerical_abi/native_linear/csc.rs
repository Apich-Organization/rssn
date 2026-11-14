//! Native CSC sparse kernels; public historical CSR names retain their actual CSC implementation.
use super::Array;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::RuleError;
use super::WORK;
use super::api;
use super::array;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::scale;
use super::sparse::Csr;
use super::sparse::inverse;
use super::sparse::result;
use super::sparse::triplets;
use super::sparse::vector;
use crate::kernels::sparse as sp;
struct Csc {
    data: sp::SparseMat,
    scale: f64,
}
impl NativeType for Csc {
    const KIND: &'static str = "rssn/native/csc/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.nrows() > 512
            || self.data.ncols() > 512
            || self.data.val().len() > LIMIT
            || self.data.val().iter().any(|v| !v.is_finite())
            || !self.scale.is_finite()
            || self.scale <= 0.0
        {
            return Err(err(
                "invalid CSC finite storage/reference or allocation budget",
            ));
        }
        Ok(())
    }
}
object_binding!(csc_csc_from_triplets,"num_linear_csc_from_triplets","crate::kernels::sparse::csr_from_triplets","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",false,(rows:usize,cols:usize,ri:Vec<usize>,ci:Vec<usize>,values:Normalized<Vec<f64>>,reference:Input)->NewObject<Csc>,|ctx,rows,cols,ri,ci,values,reference|{let _=ctx;if rows>512||cols>512{return Err(err("CSC extent budget exceeded"));}let trip=triplets(rows,cols,ri,ci,values.0)?;let(s,d)=api::reference(&reference.0)?;Ok(NewObject::new(Csc{data:sp::csr_from_triplets(rows,cols,&trip),scale:s},vec![d]))});
object_binding!(csc_csc_from_csr_data,"num_linear_csc_from_csr_data","crate::kernels::sparse::SparseMatrixData::to_csmat","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csr>)->NewObject<Csc>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if a.data.rows>512||a.data.cols>512{return Err(err("CSC extent budget exceeded"));}let data=sp::SparseMatrixData{rows:a.data.rows,cols:a.data.cols,indptr:a.data.indptr.clone(),indices:a.data.indices.clone(),data:a.data.values.clone()};Ok(NewObject::new(Csc{data:data.to_csmat(),scale:a.scale},vec![d]))})});
object_binding!(csc_csc_from_dense,"num_linear_csc_from_dense","crate::kernels::sparse::to_csr","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Array>)->NewObject<Csc>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if a.data.ndim()!=2||a.data.shape().iter().any(|n|*n>512){return Err(err("CSC source must be bounded rank-two tensor"));}Ok(NewObject::new(Csc{data:sp::to_csr(&a.data),scale:a.scale},vec![d]))})});
object_binding!(csc_csc_to_dense,"num_linear_csc_to_dense","crate::kernels::sparse::to_dense","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(array(sp::to_dense(&a.data).into_dyn(),a.scale,d)))});
object_binding!(csc_csc_transpose,"num_linear_csc_transpose","crate::kernels::sparse::transpose","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->NewObject<Csc>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(NewObject::new(Csc{data:sp::transpose(&a.data),scale:a.scale},vec![d])))});
object_binding!(csc_csc_frobenius_norm,"num_linear_csc_frobenius_norm","crate::kernels::sparse::frobenius_norm","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(sp::frobenius_norm(&a.data)*a.scale,d)?)))});
object_binding!(csc_csc_l1_norm,"num_linear_csc_l1_norm","crate::kernels::sparse::l1_norm","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(sp::l1_norm(&a.data)*a.scale,d)?)))});
object_binding!(csc_csc_linf_norm,"num_linear_csc_linf_norm","crate::kernels::sparse::linf_norm","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(sp::linf_norm(&a.data)*a.scale,d)?)))});
object_binding!(csc_csc_trace,"num_linear_csc_trace","crate::kernels::sparse::trace","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::quantity(sp::trace(&a.data).map_err(err)?*a.scale,d)?)))});
object_binding!(csc_csc_rank,"num_linear_csc_rank","crate::kernels::sparse::rank","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|{if a.data.nrows()>128||a.data.ncols()>128{return Err(err("CSC rank work budget exceeded"));}Ok(sp::rank(&a.data))})});
object_binding!(csc_csc_is_diagonal,"num_linear_csc_is_diagonal","crate::kernels::sparse::is_diagonal","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(sp::is_diagonal(&a.data)))});
object_binding!(csc_csc_is_symmetric,"num_linear_csc_is_symmetric","crate::kernels::sparse::is_symmetric","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>,epsilon:f64)->bool,|ctx,a,epsilon|{let _=ctx;if epsilon<0.0{return Err(err("symmetry tolerance must be nonnegative"));}a.with(ctx,|a|Ok(sp::is_symmetric(&a.data,epsilon)))});
object_binding!(csc_csc_matvec,"num_linear_csc_matvec","crate::kernels::sparse::sp_mat_vec_mul","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>,x:Handle<Array>)->NewObject<Array>,|ctx,a,x|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&x.value)?)?;a.with(ctx,|a|x.with(ctx,|x|{let data=sp::sp_mat_vec_mul(&a.data,vector(x,a.data.ncols())?).map_err(err)?;result(data,scale(a.scale,x.scale)?,d)}))});
object_binding!(csc_csc_cg,"num_linear_csc_cg","crate::kernels::sparse::solve_conjugate_gradient","Actual native CSC sparse kernel with explicit entry reference, preserving historical CSR function names without transposing logical rows/columns. Physical matrix/vector domains and normalized tolerances are checked; work and allocations are bounded.",true,(a:Handle<Csc>,rhs:Handle<Array>,initial:Option<Vec<f64>>,iterations:usize,tolerance:f64)->NewObject<Array>,|ctx,a,rhs,initial,iterations,tolerance|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{if a.data.nrows()!=a.data.ncols()||iterations==0||iterations>4096||tolerance<=0.0||a.data.val().len().checked_mul(iterations).is_none_or(|w|w>WORK)||initial.as_ref().is_some_and(|v|v.len()!=a.data.ncols()){return Err(err("CSC CG shape/options/work budget exceeded"));}let rhs=ndarray::Array1::from_vec(vector(r,a.data.nrows())?.to_vec());let initial=initial.map(ndarray::Array1::from_vec);let out=sp::solve_conjugate_gradient(&a.data,&rhs,initial.as_ref(),iterations,tolerance).map_err(err)?;result(out.to_vec(),scale(r.scale,1.0/a.scale)?,d)}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &csc_csc_from_triplets::BINDING,
    &csc_csc_from_csr_data::BINDING,
    &csc_csc_from_dense::BINDING,
    &csc_csc_to_dense::BINDING,
    &csc_csc_transpose::BINDING,
    &csc_csc_frobenius_norm::BINDING,
    &csc_csc_l1_norm::BINDING,
    &csc_csc_linf_norm::BINDING,
    &csc_csc_trace::BINDING,
    &csc_csc_rank::BINDING,
    &csc_csc_is_diagonal::BINDING,
    &csc_csc_is_symmetric::BINDING,
    &csc_csc_matvec::BINDING,
    &csc_csc_cg::BINDING,
];
