//! Native sparse matrix methods retain an explicit entry reference scale.
use super::Array;
use super::ArrayD;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::IxDyn;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Value;
use super::WORK;
use super::api;
use super::array;
use super::count;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::scale;
use crate::kernels::bsr;
use crate::kernels::krylov;
#[derive(Clone)]
pub(super) struct Csr {
    pub data: krylov::Csr,
    pub scale: f64,
}
impl NativeType for Csr {
    const KIND: &'static str = "rssn/native/csr/v1";
    const MAX_BYTES: usize = LIMIT * 24 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        let a = &self.data;
        if a.rows > 4096
            || a.cols > 4096
            || a.values.len() > LIMIT
            || a.indices.len() != a.values.len()
            || a.indptr.len() != a.rows + 1
            || a.indptr.first() != Some(&0)
            || a.indptr.last() != Some(&a.values.len())
            || a.indptr.windows(2).any(|w| w[0] > w[1])
            || a.indices.iter().any(|i| *i >= a.cols)
            || a.values.iter().any(|x| !x.is_finite())
            || !self.scale.is_finite()
            || self.scale <= 0.0
        {
            return Err(err(
                "invalid native CSR representation or allocation budget",
            ));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Bsr {
    data: bsr::Bsr,
    scale: f64,
}
impl NativeType for Bsr {
    const KIND: &'static str = "rssn/native/bsr/v1";
    const MAX_BYTES: usize = LIMIT * 24 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.data.validate().map_err(err)?;
        if self.data.nrows > 4096
            || self.data.ncols > 4096
            || self.data.stored_len() > LIMIT
            || self.data.br > 32
            || self.data.bc > 32
            || self.data.values.iter().any(|x| !x.is_finite())
        {
            return Err(err("block sparse allocation/rank budget exceeded"));
        }
        Ok(())
    }
}
pub(super) fn triplets(
    rows: usize,
    cols: usize,
    ri: Vec<usize>,
    ci: Vec<usize>,
    values: Vec<f64>,
) -> Result<Vec<(usize, usize, f64)>, OpError> {
    if rows > 4096
        || cols > 4096
        || ri.len() != ci.len()
        || ci.len() != values.len()
        || values.len() > LIMIT
        || ri.iter().any(|i| *i >= rows)
        || ci.iter().any(|i| *i >= cols)
    {
        return Err(err("sparse triplet shape/index/allocation budget exceeded"));
    }
    Ok(ri
        .into_iter()
        .zip(ci)
        .zip(values)
        .map(|((r, c), v)| (r, c, v))
        .collect())
}
pub(super) fn vector(
    a: &Array,
    n: usize,
) -> Result<&[f64], OpError> {
    if a.data.ndim() != 1 || a.data.len() != n {
        return Err(err("sparse/vector shape mismatch"));
    }
    a.data
        .as_slice()
        .ok_or_else(|| err("vector not contiguous"))
}
pub(super) fn result(
    data: Vec<f64>,
    scale: f64,
    d: DimExpr,
) -> Result<NewObject<Array>, OpError> {
    count(&[data.len()])?;
    Ok(array(
        ArrayD::from_shape_vec(IxDyn(&[data.len()]), data).map_err(err)?,
        scale,
        d,
    ))
}
pub(super) fn inverse(a: DimExpr) -> Result<DimExpr, OpError> {
    crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Inverse(Box::new(
            a,
        )))?,
    )
}
fn square(a: &krylov::Csr) -> Result<(), OpError> {
    if a.rows != a.cols || a.rows == 0 || a.rows > 512 {
        return Err(err(
            "factorization needs a nonempty square matrix of side at most 512",
        ));
    }
    Ok(())
}
object_binding!(sparse_csr_new,"num_linear_csr_new","crate::kernels::krylov::Csr::from_triplets","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",false,(rows:usize,cols:usize,ri:Vec<usize>,ci:Vec<usize>,values:Normalized<Vec<f64>>,reference:Input)->NewObject<Csr>,|ctx,rows,cols,ri,ci,values,reference|{let _=ctx;let trip=triplets(rows,cols,ri,ci,values.0)?;let(s,d)=api::reference(&reference.0)?;Ok(NewObject::new(Csr{data:krylov::Csr::from_triplets(rows,cols,&trip),scale:s},vec![d]))});
object_binding!(sparse_csr_matvec,"num_linear_csr_matvec","crate::kernels::krylov::Csr::matvec","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>,x:Handle<Array>)->NewObject<Array>,|ctx,a,x|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&x.value)?)?;a.with(ctx,|a|x.with(ctx,|x|{let input=vector(x,a.data.cols)?;let mut y=vec![0.0;a.data.rows];a.data.matvec(input,&mut y);result(y,scale(a.scale,x.scale)?,d)}))});
object_binding!(sparse_csr_bandwidth,"num_linear_csr_bandwidth","crate::kernels::krylov::Csr::bandwidth","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(a.data.bandwidth())})});
object_binding!(sparse_csr_rcm_ordering,"num_linear_csr_rcm_ordering","crate::kernels::krylov::Csr::rcm_ordering","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|{square(&a.data)?;Ok(a.data.rcm_ordering())})});
object_binding!(sparse_csr_min_degree_ordering,"num_linear_csr_min_degree_ordering","crate::kernels::krylov::Csr::min_degree_ordering","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|{square(&a.data)?;Ok(a.data.min_degree_ordering())})});
object_binding!(sparse_csr_amd_ordering,"num_linear_csr_amd_ordering","crate::kernels::krylov::Csr::amd_ordering","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|{square(&a.data)?;Ok(a.data.amd_ordering())})});
object_binding!(sparse_csr_nested_dissection_ordering,"num_linear_csr_nested_dissection_ordering","crate::kernels::krylov::Csr::nested_dissection_ordering","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|{square(&a.data)?;Ok(a.data.nested_dissection_ordering())})});
object_binding!(sparse_csr_diagonal,"num_linear_csr_diagonal","crate::kernels::krylov::Csr::diagonal","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let data=a.data.diagonal();result(data,a.scale,d)})});
object_binding!(sparse_csr_jacobi_inverse,"num_linear_csr_jacobi_inverse","crate::kernels::krylov::Csr::jacobi_inverse","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>)->NewObject<Array>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let data=a.data.jacobi_inverse();result(data,1.0/a.scale,inverse(d)?)})});
object_binding!(sparse_csr_permute_symmetric,"num_linear_csr_permute_symmetric","crate::kernels::krylov::Csr::permute_symmetric","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>,permutation:Vec<usize>)->NewObject<Csr>,|ctx,a,permutation|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;let mut p=permutation.clone();p.sort_unstable();if p!=(0..a.data.rows).collect::<Vec<_>>(){return Err(err("expected a complete permutation"));}Ok(NewObject::new(Csr{data:a.data.permute_symmetric(&permutation),scale:a.scale},vec![d]))})});
object_binding!(sparse_bsr_from_csr,"num_linear_bsr_from_csr","crate::kernels::bsr::Bsr::from_csr","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>,br:usize,bc:usize)->NewObject<Bsr>,|ctx,a,br,bc|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if br==0||bc==0||br>32||bc>32||a.data.values.len().checked_mul(br).and_then(|n|n.checked_mul(bc)).is_none_or(|n|n>LIMIT){return Err(err("block expansion budget exceeded"));}Ok(NewObject::new(Bsr{data:bsr::Bsr::from_csr(&a.data,br,bc).map_err(err)?,scale:a.scale},vec![d]))})});
object_binding!(sparse_bsr_from_triplets,"num_linear_bsr_from_triplets","crate::kernels::bsr::Bsr::from_triplets","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",false,(rows:usize,cols:usize,br:usize,bc:usize,ri:Vec<usize>,ci:Vec<usize>,values:Normalized<Vec<f64>>,reference:Input)->NewObject<Bsr>,|ctx,rows,cols,br,bc,ri,ci,values,reference|{let _=ctx;let trip=triplets(rows,cols,ri,ci,values.0)?;if br==0||bc==0||br>32||bc>32||trip.len().checked_mul(br).and_then(|n|n.checked_mul(bc)).is_none_or(|n|n>LIMIT){return Err(err("block expansion budget exceeded"));}let(s,d)=api::reference(&reference.0)?;Ok(NewObject::new(Bsr{data:bsr::Bsr::from_triplets(rows,cols,br,bc,&trip).map_err(err)?,scale:s},vec![d]))});
object_binding!(sparse_bsr_zeros,"num_linear_bsr_zeros","crate::kernels::bsr::Bsr::zeros","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",false,(rows:usize,cols:usize,br:usize,bc:usize,reference:Input)->NewObject<Bsr>,|ctx,rows,cols,br,bc,reference|{let _=ctx;if rows>4096||cols>4096||br>32||bc>32{return Err(err("block sparse extent budget exceeded"));}let(s,d)=api::reference(&reference.0)?;Ok(NewObject::new(Bsr{data:bsr::Bsr::zeros(rows,cols,br,bc).map_err(err)?,scale:s},vec![d]))});
object_binding!(sparse_bsr_block_cols,"num_linear_bsr_block_cols","crate::kernels::bsr::Bsr::block_cols","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.block_cols()))});
object_binding!(sparse_bsr_block_rows,"num_linear_bsr_block_rows","crate::kernels::bsr::Bsr::block_rows","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.block_rows()))});
object_binding!(sparse_bsr_block_count,"num_linear_bsr_block_count","crate::kernels::bsr::Bsr::block_count","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.block_count()))});
object_binding!(sparse_bsr_stored_len,"num_linear_bsr_stored_len","crate::kernels::bsr::Bsr::stored_len","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.stored_len()))});
object_binding!(sparse_bsr_nnz,"num_linear_bsr_nnz","crate::kernels::bsr::Bsr::nnz","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.nnz()))});
object_binding!(sparse_bsr_fill_ratio,"num_linear_bsr_fill_ratio","crate::kernels::bsr::Bsr::fill_ratio","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.fill_ratio()))});
object_binding!(sparse_bsr_validate,"num_linear_bsr_validate","crate::kernels::bsr::Bsr::validate","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{a.data.validate().map_err(err)?;Ok(true)})});
object_binding!(sparse_bsr_transpose,"num_linear_bsr_transpose","crate::kernels::bsr::Bsr::transpose","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->NewObject<Bsr>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(NewObject::new(Bsr{data:a.data.transpose(),scale:a.scale},vec![d])))});
object_binding!(sparse_bsr_to_csr,"num_linear_bsr_to_csr","crate::kernels::bsr::Bsr::to_csr","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->NewObject<Csr>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(NewObject::new(Csr{data:a.data.to_csr(),scale:a.scale},vec![d])))});
object_binding!(sparse_bsr_diagonal_blocks,"num_linear_bsr_diagonal_blocks","crate::kernels::bsr::Bsr::diagonal_blocks","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::record(vec![("normalized",a.data.diagonal_blocks().map_err(err)?.encode()?),("reference",api::quantity(a.scale,d)?)]))))});
object_binding!(sparse_bsr_matvec,"num_linear_bsr_matvec","crate::kernels::bsr::Bsr::matvec","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>,x:Handle<Array>)->NewObject<Array>,|ctx,a,x|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&x.value)?)?;a.with(ctx,|a|x.with(ctx,|x|{let input=vector(x,a.data.ncols)?;let mut y=vec![0.0;a.data.nrows];a.data.matvec(input,&mut y);result(y,scale(a.scale,x.scale)?,d)}))});
object_binding!(sparse_bsr_matvec_transpose,"num_linear_bsr_matvec_transpose","crate::kernels::bsr::Bsr::matvec_transpose","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>,x:Handle<Array>)->NewObject<Array>,|ctx,a,x|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&x.value)?)?;a.with(ctx,|a|x.with(ctx,|x|{let input=vector(x,a.data.nrows)?;let mut y=vec![0.0;a.data.ncols];a.data.matvec_transpose(input,&mut y);result(y,scale(a.scale,x.scale)?,d)}))});
object_binding!(sparse_bsr_matmul,"num_linear_bsr_matmul","crate::kernels::bsr::Bsr::matmul","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Bsr>,b:Handle<Bsr>)->NewObject<Bsr>,|ctx,a,b|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&b.value)?)?;a.with(ctx,|a|b.with(ctx,|b|{if a.data.nrows.checked_mul(b.data.ncols).is_none_or(|n|n>LIMIT)||a.data.stored_len().checked_mul(b.data.stored_len()).is_none_or(|n|n>WORK){return Err(err("sparse product allocation/work budget exceeded"));}Ok(NewObject::new(Bsr{data:a.data.matmul(&b.data).map_err(err)?,scale:scale(a.scale,b.scale)?},vec![d]))}))});
object_binding!(sparse_choose_block_size,"num_linear_choose_block_size","crate::kernels::bsr::choose_block_size","Native sparse matrix operation with retained entry reference and physical unit domain. Indices are zero based; normalized array storage, capability kind/session and size/work budgets are checked.",true,(a:Handle<Csr>,candidates:Vec<[usize;2]>)->Vec<usize>,|ctx,a,candidates|{let _=ctx;if candidates.len()>64||candidates.iter().flatten().any(|n|*n==0||*n>32){return Err(err("block candidate budget exceeded"));}a.with(ctx,|a|{let candidates=candidates.iter().map(|x|(x[0],x[1])).collect::<Vec<_>>();let(r,c)=bsr::choose_block_size(&a.data,&candidates);Ok(vec![r,c])})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &sparse_ilut_stats::BINDING,
    &sparse_sparse_native_matvec::BINDING,
    &sparse_sparse_native_lu_solve::BINDING,
    &sparse_sparse_native_cholesky_solve::BINDING,
    &sparse_sparse_native_solve_auto::BINDING,
    &sparse_sparse_native_cg::BINDING,
    &sparse_sparse_native_cg_precond::BINDING,
    &sparse_sparse_native_krylov_precond::BINDING,
    &sparse_cg_with::BINDING,
    &sparse_bicgstab_with::BINDING,
    &sparse_gmres_with::BINDING,
    &sparse_fgmres::BINDING,
    &sparse_bicgstab_l::BINDING,
    &sparse_idr_s::BINDING,
    &sparse_idr_s_with::BINDING,
    &sparse_csr_ilu0::BINDING,
    &sparse_csr_ic0::BINDING,
    &sparse_auto_new::BINDING,
    &sparse_auto_preconditioner::BINDING,
    &sparse_csr_sparse_lu::BINDING,
    &sparse_csr_sparse_lu_with::BINDING,
    &sparse_csr_ilut::BINDING,
    &sparse_bsr_block_jacobi::BINDING,
    &sparse_bsr_bilu0::BINDING,
    &sparse_ilu0_apply::BINDING,
    &sparse_ic0_apply::BINDING,
    &sparse_auto_apply::BINDING,
    &sparse_ilut_apply::BINDING,
    &sparse_block_jacobi_apply::BINDING,
    &sparse_bilu0_apply::BINDING,
    &sparse_as_preconditioner_apply::BINDING,
    &sparse_preconditioner_apply::BINDING,
    &sparse_block_jacobi_block_size::BINDING,
    &sparse_auto_kind::BINDING,
    &sparse_sparse_lu_dim::BINDING,
    &sparse_sparse_lu_nnz::BINDING,
    &sparse_sparse_lu_permutation::BINDING,
    &sparse_sparse_lu_solve::BINDING,
    &sparse_identity_precond::BINDING,
    &sparse_csr_new::BINDING,
    &sparse_csr_matvec::BINDING,
    &sparse_csr_bandwidth::BINDING,
    &sparse_csr_rcm_ordering::BINDING,
    &sparse_csr_min_degree_ordering::BINDING,
    &sparse_csr_amd_ordering::BINDING,
    &sparse_csr_nested_dissection_ordering::BINDING,
    &sparse_csr_diagonal::BINDING,
    &sparse_csr_jacobi_inverse::BINDING,
    &sparse_csr_permute_symmetric::BINDING,
    &sparse_bsr_from_csr::BINDING,
    &sparse_bsr_from_triplets::BINDING,
    &sparse_bsr_zeros::BINDING,
    &sparse_bsr_block_cols::BINDING,
    &sparse_bsr_block_rows::BINDING,
    &sparse_bsr_block_count::BINDING,
    &sparse_bsr_stored_len::BINDING,
    &sparse_bsr_nnz::BINDING,
    &sparse_bsr_fill_ratio::BINDING,
    &sparse_bsr_validate::BINDING,
    &sparse_bsr_transpose::BINDING,
    &sparse_bsr_to_csr::BINDING,
    &sparse_bsr_diagonal_blocks::BINDING,
    &sparse_bsr_matvec::BINDING,
    &sparse_bsr_matvec_transpose::BINDING,
    &sparse_bsr_matmul::BINDING,
    &sparse_choose_block_size::BINDING,
];

struct Factor<T> {
    data: T,
    n: usize,
    scale: f64,
}
macro_rules! native_factor {
    ($t:path,$kind:literal) => {
        impl NativeType for Factor<$t> {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = 8 * 1024 * 1024;

            fn validate(&self) -> Result<(), OpError> {
                if self.n == 0 || self.n > 512 || !self.scale.is_finite() || self.scale <= 0.0 {
                    return Err(err(
                        "invalid factor dimension/reference or storage budget",
                    ));
                }
                Ok(())
            }
        }
    };
}
native_factor!(krylov::Ilu0, "rssn/native/ilu0/v1");
native_factor!(krylov::Ic0, "rssn/native/ic0/v1");
native_factor!(
    krylov::AutoPreconditioner,
    "rssn/native/auto-preconditioner/v1"
);
native_factor!(krylov::SparseLu, "rssn/native/sparse-lu/v1");
native_factor!(krylov::Ilut, "rssn/native/ilut/v1");
native_factor!(bsr::BlockJacobi, "rssn/native/block-jacobi/v1");
native_factor!(bsr::Bilu0, "rssn/native/bilu0/v1");
fn ordering(name: &str) -> Result<krylov::SparseOrdering, OpError> {
    Ok(match name {
        | "Natural" => krylov::SparseOrdering::Natural,
        | "Rcm" => krylov::SparseOrdering::Rcm,
        | "MinDegree" => krylov::SparseOrdering::MinDegree,
        | "Amd" => krylov::SparseOrdering::Amd,
        | "NestedDissection" => krylov::SparseOrdering::NestedDissection,
        | _ => return Err(err("unknown sparse ordering")),
    })
}
object_binding!(sparse_csr_ilu0,"num_linear_csr_ilu0","crate::kernels::krylov::Csr::ilu0","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>)->api::Maybe<NewObject<Factor<krylov::Ilu0>>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;Ok(api::Maybe(a.data.ilu0().map(|data|NewObject::new(Factor{data,n:a.data.rows,scale:a.scale},vec![d]))))})});
object_binding!(sparse_csr_ic0,"num_linear_csr_ic0","crate::kernels::krylov::Csr::ic0","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>)->api::Maybe<NewObject<Factor<krylov::Ic0>>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;Ok(api::Maybe(a.data.ic0().map(|data|NewObject::new(Factor{data,n:a.data.rows,scale:a.scale},vec![d]))))})});
object_binding!(sparse_auto_new,"num_linear_auto_new","crate::kernels::krylov::AutoPreconditioner::new","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>)->NewObject<Factor<krylov::AutoPreconditioner>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;Ok(NewObject::new(Factor{data:krylov::AutoPreconditioner::new(&a.data),n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(sparse_auto_preconditioner,"num_linear_auto_preconditioner","crate::kernels::krylov::auto_preconditioner","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>)->NewObject<Factor<krylov::AutoPreconditioner>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;Ok(NewObject::new(Factor{data:krylov::auto_preconditioner(&a.data),n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(sparse_csr_sparse_lu,"num_linear_csr_sparse_lu","crate::kernels::krylov::Csr::sparse_lu","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>,ordering_name:String)->NewObject<Factor<krylov::SparseLu>>,|ctx,a,ordering_name|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;Ok(NewObject::new(Factor{data:a.data.sparse_lu(ordering(&ordering_name)?).map_err(|e|err(format!("{e:?}")))?,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(sparse_csr_sparse_lu_with,"num_linear_csr_sparse_lu_with","crate::kernels::krylov::Csr::sparse_lu_with","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>,ordering_name:String,threshold:f64)->NewObject<Factor<krylov::SparseLu>>,|ctx,a,ordering_name,threshold|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;if !(0.0..=1.0).contains(&threshold){return Err(err("pivot threshold must be in [0,1]"));}Ok(NewObject::new(Factor{data:a.data.sparse_lu_with(ordering(&ordering_name)?,threshold).map_err(|e|err(format!("{e:?}")))?,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(sparse_csr_ilut,"num_linear_csr_ilut","crate::kernels::krylov::Csr::ilut","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Csr>,drop_tolerance:f64,max_entries:usize,pivot_tolerance:f64)->NewObject<Factor<krylov::Ilut>>,|ctx,a,drop_tolerance,max_entries,pivot_tolerance|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{square(&a.data)?;if max_entries>128{return Err(err("ILUT fill budget exceeded"));}let options=krylov::IlutOptions{drop_tolerance,max_entries_per_side:max_entries,pivot_tolerance};let data=a.data.ilut(options).map_err(|e|err(format!("{e:?}")))?;Ok(NewObject::new(Factor{data,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(sparse_bsr_block_jacobi,"num_linear_bsr_block_jacobi","crate::kernels::bsr::Bsr::block_jacobi","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Bsr>)->NewObject<Factor<bsr::BlockJacobi>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if a.data.nrows!=a.data.ncols||a.data.nrows==0||a.data.nrows>512{return Err(err("block factorization needs bounded square matrix"));}Ok(NewObject::new(Factor{data:a.data.block_jacobi().map_err(err)?,n:a.data.nrows,scale:a.scale},vec![d]))})});
object_binding!(sparse_bsr_bilu0,"num_linear_bsr_bilu0","crate::kernels::bsr::Bsr::bilu0","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Bsr>)->NewObject<Factor<bsr::Bilu0>>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{if a.data.nrows!=a.data.ncols||a.data.nrows==0||a.data.nrows>512{return Err(err("block factorization needs bounded square matrix"));}Ok(NewObject::new(Factor{data:a.data.bilu0().map_err(err)?,n:a.data.nrows,scale:a.scale},vec![d]))})});
object_binding!(sparse_ilu0_apply,"num_linear_ilu0_apply","crate::kernels::krylov::Ilu0::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::Ilu0>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_ic0_apply,"num_linear_ic0_apply","crate::kernels::krylov::Ic0::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::Ic0>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_auto_apply,"num_linear_auto_apply","crate::kernels::krylov::AutoPreconditioner::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_ilut_apply,"num_linear_ilut_apply","crate::kernels::krylov::Ilut::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::Ilut>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out).map_err(|e|err(format!("{e:?}")))?;result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_block_jacobi_apply,"num_linear_block_jacobi_apply","crate::kernels::bsr::BlockJacobi::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<bsr::BlockJacobi>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_bilu0_apply,"num_linear_bilu0_apply","crate::kernels::bsr::Bilu0::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<bsr::Bilu0>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_as_preconditioner_apply,"num_linear_as_preconditioner_apply","crate::kernels::krylov::AutoPreconditioner::as_preconditioner","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];let p=a.data.as_preconditioner();p.apply(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_preconditioner_apply,"num_linear_preconditioner_apply","crate::kernels::krylov::Preconditioner<'_>::apply","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];let p=a.data.as_preconditioner();krylov::Preconditioner::apply(&p,input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_block_jacobi_block_size,"num_linear_block_jacobi_block_size","crate::kernels::bsr::BlockJacobi::block_size","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<bsr::BlockJacobi>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.block_size()))});
object_binding!(sparse_auto_kind,"num_linear_auto_kind","crate::kernels::krylov::AutoPreconditioner::kind","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::AutoPreconditioner>>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(format!("{:?}",a.data.kind())))});
object_binding!(sparse_sparse_lu_dim,"num_linear_sparse_lu_dim","crate::kernels::krylov::SparseLu::dim","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::SparseLu>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.dim()))});
object_binding!(sparse_sparse_lu_nnz,"num_linear_sparse_lu_nnz","crate::kernels::krylov::SparseLu::nnz","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::SparseLu>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.nnz()))});
object_binding!(sparse_sparse_lu_permutation,"num_linear_sparse_lu_permutation","crate::kernels::krylov::SparseLu::permutation","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::SparseLu>>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.permutation().to_vec()))});
object_binding!(sparse_sparse_lu_solve,"num_linear_sparse_lu_solve","crate::kernels::krylov::SparseLu::solve","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(a:Handle<Factor<krylov::SparseLu>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;result(a.data.solve(input),scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_identity_precond,"num_linear_identity_precond","crate::kernels::krylov::identity_precond","Native sparse factor/preconditioner with its matrix entry reference. Application returns a normalized result with physical reference RHS/matrix; shapes, unit domains and factor work are checked.",true,(rhs:Handle<Array>)->NewObject<Array>,|ctx,rhs|{let _=ctx;let d=dim(ctx,&rhs.value)?;rhs.with(ctx,|r|{let input=vector(r,r.data.len())?;let mut out=vec![0.0;input.len()];krylov::identity_precond(input,&mut out);result(out,r.scale,d)})});

fn borrowed(a: &krylov::Csr) -> crate::kernels::sparse_native::Csr<'_> {
    crate::kernels::sparse_native::Csr {
        rows: a.rows,
        cols: a.cols,
        row_ptr: &a.indptr,
        col_idx: &a.indices,
        values: &a.values,
    }
}
fn iterative_budget(
    a: &krylov::Csr,
    tolerance: f64,
    iterations: usize,
) -> Result<(), OpError> {
    square(a)?;
    if tolerance <= 0.0
        || iterations == 0
        || iterations > 4096
        || a.values
            .len()
            .checked_mul(iterations)
            .is_none_or(|n| n > WORK)
    {
        return Err(err("iterative solver tolerance/work budget exceeded"));
    }
    Ok(())
}
fn preconditioner(name: &str) -> Result<crate::kernels::sparse_native::PrecondKind, OpError> {
    use crate::kernels::sparse_native::PrecondKind as K;
    match name {
        | "Identity" => Ok(K::Identity),
        | "Jacobi" => Ok(K::Jacobi),
        | "BlockJacobi" => Ok(K::BlockJacobi),
        | "Ilu0" => Ok(K::Ilu0),
        | "Ic0" => Ok(K::Ic0),
        | "Auto" => Ok(K::Auto),
        | _ => Err(err("unknown native preconditioner")),
    }
}
fn iter_output(
    out: krylov::IterResult,
    reference: f64,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("solution", out.x.encode()?),
        ("reference", api::quantity(reference, d)?),
        ("iterations", out.iterations.encode()?),
        ("residual", out.residual.encode()?),
        ("converged", out.converged.encode()?),
    ])))
}
object_binding!(sparse_sparse_native_matvec,"num_linear_sparse_native_matvec","crate::kernels::sparse_native::Csr<'_>::matvec","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,x:Handle<Array>)->NewObject<Array>,|ctx,a,x|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&x.value)?)?;a.with(ctx,|a|x.with(ctx,|x|{let out=borrowed(&a.data).matvec(vector(x,a.data.cols)?).map_err(err)?;result(out,scale(a.scale,x.scale)?,d)}))});
object_binding!(sparse_sparse_native_lu_solve,"num_linear_sparse_native_lu_solve","crate::kernels::sparse_native::lu_solve","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{square(&a.data)?;let input=vector(b,a.data.rows)?;let out=crate::kernels::sparse_native::lu_solve(&borrowed(&a.data),input).map_err(err)?;result(out,scale(b.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_sparse_native_cholesky_solve,"num_linear_sparse_native_cholesky_solve","crate::kernels::sparse_native::cholesky_solve","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{square(&a.data)?;let input=vector(b,a.data.rows)?;let out=crate::kernels::sparse_native::cholesky_solve(&borrowed(&a.data),input).map_err(err)?;result(out,scale(b.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_sparse_native_solve_auto,"num_linear_sparse_native_solve_auto","crate::kernels::sparse_native::solve_auto","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>)->Dynamic,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{square(&a.data)?;let input=vector(b,a.data.rows)?;let out=crate::kernels::sparse_native::solve_auto(&borrowed(&a.data),input).map_err(err)?;Ok(Dynamic(api::record(vec![("solution",out.0.encode()?),("reference",api::quantity(scale(b.scale,1.0/a.scale)?,d)?),("algorithm",Value::String(out.1.into()))])))}))});
object_binding!(sparse_sparse_native_cg,"num_linear_sparse_native_cg","crate::kernels::sparse_native::cg","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,tolerance:f64,max_iterations:usize)->Dynamic,|ctx,a,rhs,tolerance,max_iterations|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;let input=vector(b,a.data.rows)?;let out=crate::kernels::sparse_native::cg(&borrowed(&a.data),input,tolerance,max_iterations).map_err(err)?;Ok(Dynamic(api::record(vec![("solution",out.x.encode()?),("reference",api::quantity(scale(b.scale,1.0/a.scale)?,d)?),("iterations",out.iterations.encode()?),("residual",out.residual.encode()?),("converged",out.converged.encode()?)])))}))});
object_binding!(sparse_sparse_native_cg_precond,"num_linear_sparse_native_cg_precond","crate::kernels::sparse_native::cg_precond","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,tolerance:f64,max_iterations:usize,kind:String,block_rows:usize,block_cols:usize)->Dynamic,|ctx,a,rhs,tolerance,max_iterations,kind,block_rows,block_cols|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;let input=vector(b,a.data.rows)?;let block=if block_rows==0&&block_cols==0{None}else{if block_rows==0||block_cols==0||block_rows>32||block_cols>32{return Err(err("block shape budget invalid"));}Some((block_rows,block_cols))};let out=crate::kernels::sparse_native::cg_precond(&borrowed(&a.data),input,tolerance,max_iterations,preconditioner(&kind)?,block).map_err(err)?;Ok(Dynamic(api::record(vec![("solution",out.x.encode()?),("reference",api::quantity(scale(b.scale,1.0/a.scale)?,d)?),("iterations",out.iterations.encode()?),("residual",out.residual.encode()?),("converged",out.converged.encode()?),("preconditioner",Value::String(out.precond.into()))])))}))});
object_binding!(sparse_sparse_native_krylov_precond,"num_linear_sparse_native_krylov_precond","crate::kernels::sparse_native::krylov_precond","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,tolerance:f64,max_iterations:usize,kind:String,block_rows:usize,block_cols:usize,method_name:String,restart:usize)->Dynamic,|ctx,a,rhs,tolerance,max_iterations,kind,block_rows,block_cols,method_name,restart|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;let input=vector(b,a.data.rows)?;let block=if block_rows==0&&block_cols==0{None}else{if block_rows==0||block_cols==0||block_rows>32||block_cols>32{return Err(err("block shape budget invalid"));}Some((block_rows,block_cols))};let method=match method_name.as_str(){"Bicgstab"=>crate::kernels::sparse_native::Krylov::Bicgstab,"Gmres"=>{if restart==0||restart>64{return Err(err("GMRES restart budget exceeded"));}crate::kernels::sparse_native::Krylov::Gmres{restart}},_=>return Err(err("Krylov method must be Gmres or Bicgstab"))};let out=crate::kernels::sparse_native::krylov_precond(&borrowed(&a.data),input,tolerance,max_iterations,method,preconditioner(&kind)?,block).map_err(err)?;Ok(Dynamic(api::record(vec![("solution",out.x.encode()?),("reference",api::quantity(scale(b.scale,1.0/a.scale)?,d)?),("iterations",out.iterations.encode()?),("residual",out.residual.encode()?),("converged",out.converged.encode()?),("preconditioner",Value::String(out.precond.into()))])))}))});
object_binding!(sparse_cg_with,"num_linear_cg_with","crate::kernels::krylov::cg_with","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,p:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize)->Dynamic,|ctx,a,p,rhs,initial,tolerance,max_iterations|{let _=ctx;api::same_dimension(&dim(ctx,&a.value)?,&dim(ctx,&p.value)?)?;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|p.with(ctx,|p|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if p.n!=a.data.rows||p.scale.to_bits()!=a.scale.to_bits()||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("preconditioner/initial frame mismatch"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=p.data.as_preconditioner();let out=krylov::cg_with(apply,&pre,input,initial.as_deref(),tolerance,max_iterations);iter_output(out,scale(b.scale,1.0/a.scale)?,d)})))});
object_binding!(sparse_bicgstab_with,"num_linear_bicgstab_with","crate::kernels::krylov::bicgstab_with","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,p:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize)->Dynamic,|ctx,a,p,rhs,initial,tolerance,max_iterations|{let _=ctx;api::same_dimension(&dim(ctx,&a.value)?,&dim(ctx,&p.value)?)?;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|p.with(ctx,|p|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if p.n!=a.data.rows||p.scale.to_bits()!=a.scale.to_bits()||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("preconditioner/initial frame mismatch"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=p.data.as_preconditioner();let out=krylov::bicgstab_with(apply,&pre,input,initial.as_deref(),tolerance,max_iterations);iter_output(out,scale(b.scale,1.0/a.scale)?,d)})))});
object_binding!(sparse_gmres_with,"num_linear_gmres_with","crate::kernels::krylov::gmres_with","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,p:Handle<Factor<krylov::AutoPreconditioner>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize,restart:usize)->Dynamic,|ctx,a,p,rhs,initial,tolerance,max_iterations,restart|{let _=ctx;api::same_dimension(&dim(ctx,&a.value)?,&dim(ctx,&p.value)?)?;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|p.with(ctx,|p|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if p.n!=a.data.rows||p.scale.to_bits()!=a.scale.to_bits()||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("preconditioner/initial frame mismatch"));}
if restart==0||restart>64{return Err(err("GMRES restart budget exceeded"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=p.data.as_preconditioner();let out=krylov::gmres_with(apply,&pre,input,initial.as_deref(),restart,tolerance,max_iterations);iter_output(out,scale(b.scale,1.0/a.scale)?,d)})))});
object_binding!(sparse_fgmres,"num_linear_fgmres","crate::kernels::krylov_ext::fgmres","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize,space:usize)->Dynamic,|ctx,a,rhs,initial,tolerance,max_iterations,space|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if space==0||space>64||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("Krylov space/initial shape budget invalid"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=krylov::identity_precond;let out=crate::kernels::krylov_ext::fgmres(apply,pre,input,initial.as_deref(),space,tolerance,max_iterations);iter_output(out,scale(b.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_bicgstab_l,"num_linear_bicgstab_l","crate::kernels::krylov_ext::bicgstab_l","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize,space:usize)->Dynamic,|ctx,a,rhs,initial,tolerance,max_iterations,space|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if space==0||space>64||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("Krylov space/initial shape budget invalid"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=krylov::identity_precond;let out=crate::kernels::krylov_ext::bicgstab_l(apply,pre,space,input,initial.as_deref(),tolerance,max_iterations);iter_output(out,scale(b.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_idr_s,"num_linear_idr_s","crate::kernels::krylov_ext::idr_s","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize,space:usize,seed:u64)->Dynamic,|ctx,a,rhs,initial,tolerance,max_iterations,space,seed|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if space==0||space>64||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("Krylov space/initial shape budget invalid"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=krylov::identity_precond;#[allow(deprecated,reason="Explicit compatibility coverage must call the native deprecated entry point.")]let out=crate::kernels::krylov_ext::idr_s(apply,pre,space,input,initial.as_deref(),tolerance,max_iterations,seed);iter_output(out,scale(b.scale,1.0/a.scale)?,d)}))});
object_binding!(sparse_idr_s_with,"num_linear_idr_s_with","crate::kernels::krylov_ext::idr_s_with","Native sparse solver with explicit matrix/RHS reference domains. Solutions retain RHS/matrix units; initial vectors use that normalized solution reference. Native method/preconditioner choices, residuals and convergence reports are preserved; shape and work budgets are checked.",true,(a:Handle<Csr>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_iterations:usize,space:usize,seed:u64)->Dynamic,|ctx,a,rhs,initial,tolerance,max_iterations,space,seed|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|b|{iterative_budget(&a.data,tolerance,max_iterations)?;if space==0||space>64||initial.as_ref().is_some_and(|v|v.len()!=a.data.rows){return Err(err("Krylov space/initial shape budget invalid"));}let input=vector(b,a.data.rows)?;let apply=|x:&[f64],y:&mut[f64]|a.data.matvec(x,y);let pre=krylov::identity_precond;let options=crate::kernels::krylov_ext::IdrOptions{s:space,tol:tolerance,max_iter:max_iterations,seed};let out=crate::kernels::krylov_ext::idr_s_with(apply,pre,input,initial.as_deref(),&options);iter_output(out,scale(b.scale,1.0/a.scale)?,d)}))});

object_binding!(sparse_ilut_stats,"num_linear_ilut_stats","crate::kernels::krylov::Ilut::stats","Actual retained ILUT storage/drop diagnostics, with dimensionless counts.",true,(a:Handle<Factor<krylov::Ilut>>)->Dynamic,|ctx,a|a.with(ctx,|a|{let s=a.data.stats();Ok(Dynamic(api::record(vec![("lower_entries",s.lower_entries.encode()?),("upper_entries",s.upper_entries.encode()?),("dropped_entries",s.dropped_entries.encode()?)])))}));
