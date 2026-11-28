//! Bounded exact sparse arithmetic, specialized to the rationals and checked prime fields.
use super::super::native_object_api::Maybe;
use super::Dynamic;
use super::Handle;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Value;
use super::api;
use super::err;
use super::object_binding;
use crate::kernels::sparse_exact::Domain;
use crate::kernels::sparse_exact::{
    self as se,
};
use num_bigint::BigInt;
use num_rational::BigRational;
#[path = "exact_prepared.rs"]
mod prepared;
fn bounded_q(q: &BigRational) -> Result<(), OpError> {
    if q.numer().bits() > 16384 || q.denom().bits() > 16384 {
        return Err(err("exact coefficient bit budget exceeded"));
    }
    Ok(())
}
fn exact_shape(
    r: usize,
    c: usize,
) -> Result<(), OpError> {
    if r > 32 || c > 32 {
        return Err(err("exact sparse side budget is 32"));
    }
    Ok(())
}
impl NativeType for se::SparseMat<BigRational> {
    const KIND: &'static str = "rssn/native/sparse-rational/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        exact_shape(self.rows, self.cols)?;
        for (_, _, q) in self.triplets() {
            bounded_q(q)?;
        }
        Ok(())
    }
}
impl NativeType for se::Echelon<BigRational> {
    const KIND: &'static str = "rssn/native/sparse-echelon-rational/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.cols > 32 || self.rows.len() > 32 || self.ops.len() > 1024 {
            return Err(err("exact echelon size budget exceeded"));
        }
        for row in &self.rows {
            for (_, q) in row {
                bounded_q(q)?;
            }
        }
        for (_, _, q) in &self.ops {
            bounded_q(q)?;
        }
        Ok(())
    }
}
#[derive(Clone)]
struct PrimeMatrix {
    data: se::SparseMat<u64>,
    field: se::Gfp,
}
impl NativeType for PrimeMatrix {
    const KIND: &'static str = "rssn/native/sparse-prime/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.rows > 256 || self.data.cols > 256 || self.data.nnz() > 65536 {
            return Err(err("prime sparse size budget exceeded"));
        }
        prime(self.field.p)?;
        Ok(())
    }
}
fn prime(p: u64) -> Result<se::Gfp, OpError> {
    if !se::is_prime_u64(p) {
        return Err(err("finite-field modulus must be prime"));
    }
    Ok(se::Gfp::new(p))
}
fn exact_flat(a: &se::SparseMat<BigRational>) -> Result<Dynamic, OpError> {
    let dense = a.to_dense(&mut se::Rationals);
    Ok(Dynamic(api::record(vec![
        ("rows", a.rows.encode()?),
        ("cols", a.cols.encode()?),
        (
            "values",
            dense.into_iter().flatten().collect::<Vec<_>>().encode()?,
        ),
    ])))
}
object_binding!(exact_rational_domain_zero,"num_linear_rational_domain_zero","crate::kernels::sparse_exact::Domain::zero","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,()->BigRational,|ctx,|{let _=ctx;Ok(se::Rationals.zero())});
object_binding!(exact_rational_domain_one,"num_linear_rational_domain_one","crate::kernels::sparse_exact::Domain::one","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,()->BigRational,|ctx,|{let _=ctx;Ok(se::Rationals.one())});
object_binding!(exact_rational_domain_is_zero,"num_linear_rational_domain_is_zero","crate::kernels::sparse_exact::Domain::is_zero","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational)->bool,|ctx,a|{let _=ctx;bounded_q(&a)?;Ok(se::Rationals.is_zero(&a))});
object_binding!(exact_rational_domain_neg,"num_linear_rational_domain_neg","crate::kernels::sparse_exact::Domain::neg","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational)->BigRational,|ctx,a|{let _=ctx;bounded_q(&a)?;Ok(se::Rationals.neg(&a))});
object_binding!(exact_rational_domain_add,"num_linear_rational_domain_add","crate::kernels::sparse_exact::Domain::add","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational,b:BigRational)->BigRational,|ctx,a,b|{let _=ctx;bounded_q(&a)?;bounded_q(&b)?;Ok(se::Rationals.add(&a,&b))});
object_binding!(exact_rational_domain_sub,"num_linear_rational_domain_sub","crate::kernels::sparse_exact::Domain::sub","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational,b:BigRational)->BigRational,|ctx,a,b|{let _=ctx;bounded_q(&a)?;bounded_q(&b)?;Ok(se::Rationals.sub(&a,&b))});
object_binding!(exact_rational_domain_mul,"num_linear_rational_domain_mul","crate::kernels::sparse_exact::Domain::mul","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational,b:BigRational)->BigRational,|ctx,a,b|{let _=ctx;bounded_q(&a)?;bounded_q(&b)?;Ok(se::Rationals.mul(&a,&b))});
object_binding!(exact_rational_domain_div,"num_linear_rational_domain_div","crate::kernels::sparse_exact::Domain::div","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:BigRational,b:BigRational)->Option<BigRational>,|ctx,a,b|{let _=ctx;bounded_q(&a)?;bounded_q(&b)?;Ok(se::Rationals.div(&a,&b))});
object_binding!(exact_gfp_mul_mod,"num_linear_gfp_mul_mod","crate::kernels::sparse_exact::Gfp::mul_mod","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(p:u64,a:u64,b:u64)->u64,|ctx,p,a,b|{let _=ctx;let f=prime(p)?;Ok(f.mul_mod(a,b))});
object_binding!(exact_gfp_pow_mod,"num_linear_gfp_pow_mod","crate::kernels::sparse_exact::Gfp::pow_mod","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(p:u64,a:u64,exponent:u64)->u64,|ctx,p,a,exponent|{let _=ctx;let f=prime(p)?;Ok(f.pow_mod(a,exponent))});
object_binding!(exact_gfp_inv,"num_linear_gfp_inv","crate::kernels::sparse_exact::Gfp::inv","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(p:u64,a:u64)->Option<u64>,|ctx,p,a|{let _=ctx;let f=prime(p)?;Ok(f.inv(a))});
object_binding!(exact_gfp_reduce,"num_linear_gfp_reduce","crate::kernels::sparse_exact::Gfp::reduce","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(p:u64,a:BigInt)->u64,|ctx,p,a|{let _=ctx;let f=prime(p)?;if a.bits()>16384{return Err(err("integer bit budget exceeded"));}Ok(f.reduce(&a))});
object_binding!(exact_word_primes,"num_linear_word_primes","crate::kernels::sparse_exact::word_primes","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(count:usize)->Vec<u64>,|ctx,count|{let _=ctx;if count>1024{return Err(err("prime iterator count budget exceeded"));}Ok(se::word_primes().take(count).collect())});
object_binding!(exact_rational_sparse_zeros,"num_linear_rational_sparse_zeros","crate::kernels::sparse_exact::SparseMat<E>::zeros","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(rows:usize,cols:usize)->NewObject<se::SparseMat<BigRational>>,|ctx,rows,cols|{let _=ctx;exact_shape(rows,cols)?;Ok(NewObject::normalized(se::SparseMat::zeros(rows,cols)))});
object_binding!(exact_rational_sparse_from_triplets,"num_linear_rational_sparse_from_triplets","crate::kernels::sparse_exact::SparseMat<E>::from_triplets","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(rows:usize,cols:usize,ri:Vec<usize>,ci:Vec<usize>,values:Vec<BigRational>)->NewObject<se::SparseMat<BigRational>>,|ctx,rows,cols,ri,ci,values|{let _=ctx;exact_shape(rows,cols)?;if ri.len()!=ci.len()||ri.len()!=values.len()||values.len()>1024{return Err(err("exact triplet shape/allocation budget exceeded"));}for q in &values{bounded_q(q)?;}let trip=ri.into_iter().zip(ci).zip(values).map(|((r,c),v)|(r,c,v));Ok(NewObject::normalized(se::SparseMat::from_triplets(&mut se::Rationals,rows,cols,trip).map_err(err)?))});
object_binding!(exact_rational_sparse_from_dense,"num_linear_rational_sparse_from_dense","crate::kernels::sparse_exact::SparseMat<E>::from_dense","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(rows:usize,cols:usize,values:Vec<BigRational>)->NewObject<se::SparseMat<BigRational>>,|ctx,rows,cols,values|{let _=ctx;exact_shape(rows,cols)?;if rows.checked_mul(cols)!=Some(values.len()){return Err(err("exact dense shape/storage mismatch"));}for q in &values{bounded_q(q)?;}if rows==0||cols==0{return Ok(NewObject::normalized(se::SparseMat::zeros(rows,cols)));}let dense=values.chunks(cols).map(<[num_rational::BigRational]>::to_vec).collect::<Vec<_>>();Ok(NewObject::normalized(se::SparseMat::from_dense(&mut se::Rationals,&dense)))});
object_binding!(exact_rational_sparse_to_dense,"num_linear_rational_sparse_to_dense","crate::kernels::sparse_exact::SparseMat<E>::to_dense","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,exact_flat)});
object_binding!(exact_rational_sparse_triplets,"num_linear_rational_sparse_triplets","crate::kernels::sparse_exact::SparseMat<E>::triplets","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let trip=a.triplets().collect::<Vec<_>>();Ok(Dynamic(api::record(vec![("rows",trip.iter().map(|(r,_,_)|*r).collect::<Vec<_>>().encode()?),("cols",trip.iter().map(|(_,c,_)|*c).collect::<Vec<_>>().encode()?),("values",trip.iter().map(|(_,_,v)|(*v).clone()).collect::<Vec<_>>().encode()?)])))})});
object_binding!(exact_rational_sparse_nnz,"num_linear_rational_sparse_nnz","crate::kernels::sparse_exact::SparseMat<E>::nnz","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.nnz()))});
object_binding!(exact_rational_sparse_transpose,"num_linear_rational_sparse_transpose","crate::kernels::sparse_exact::SparseMat<E>::transpose","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->NewObject<se::SparseMat<BigRational>>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(NewObject::normalized(a.transpose())))});
object_binding!(exact_rational_sparse_scale,"num_linear_rational_sparse_scale","crate::kernels::sparse_exact::SparseMat<E>::scale","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,factor:BigRational)->NewObject<se::SparseMat<BigRational>>,|ctx,a,factor|{let _=ctx;bounded_q(&factor)?;a.with(ctx,|a|Ok(NewObject::normalized(a.scale(&mut se::Rationals,&factor))))});
object_binding!(exact_rational_sparse_add,"num_linear_rational_sparse_add","crate::kernels::sparse_exact::SparseMat<E>::add","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,b:Handle<se::SparseMat<BigRational>>)->Maybe<NewObject<se::SparseMat<BigRational>>>,|ctx,a,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|Ok(Maybe(a.add(&mut se::Rationals,b).map(NewObject::normalized)))))});
object_binding!(exact_rational_sparse_matmul,"num_linear_rational_sparse_matmul","crate::kernels::sparse_exact::SparseMat<E>::matmul","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,b:Handle<se::SparseMat<BigRational>>)->Maybe<NewObject<se::SparseMat<BigRational>>>,|ctx,a,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|Ok(Maybe(a.matmul(&mut se::Rationals,b).map(NewObject::normalized)))))});
object_binding!(exact_rational_sparse_matvec,"num_linear_rational_sparse_matvec","crate::kernels::sparse_exact::SparseMat<E>::matvec","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,x:Vec<BigRational>)->Option<Vec<BigRational>>,|ctx,a,x|{let _=ctx;for q in &x{bounded_q(q)?;}a.with(ctx,|a|Ok(a.matvec(&mut se::Rationals,&x)))});
object_binding!(exact_rational_sparse_det,"num_linear_rational_sparse_det","crate::kernels::sparse_exact::det","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->Option<BigRational>,|ctx,a|{let _=ctx;a.with(ctx,|a|{let result=se::det(&mut se::Rationals,a);Ok(result)})});
object_binding!(exact_rational_sparse_rank,"num_linear_rational_sparse_rank","crate::kernels::sparse_exact::rank","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|{let result=se::rank(&mut se::Rationals,a);Ok(result)})});
object_binding!(exact_rational_sparse_nullspace,"num_linear_rational_sparse_nullspace","crate::kernels::sparse_exact::nullspace","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let result=se::nullspace(&mut se::Rationals,a);let rows=result.len();let cols=result.first().map_or(a.cols,Vec::len);Ok(Dynamic(api::record(vec![("rows",rows.encode()?),("cols",cols.encode()?),("values",result.into_iter().flatten().collect::<Vec<_>>().encode()?)])))})});
object_binding!(exact_rational_sparse_det_modular,"num_linear_rational_sparse_det_modular","crate::kernels::sparse_exact::det_rational","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>)->Option<BigRational>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(se::det_rational(a)))});
object_binding!(exact_rational_sparse_solve_dixon,"num_linear_rational_sparse_solve_dixon","crate::kernels::sparse_exact::solve_dixon","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,b:Vec<BigRational>)->Option<Vec<BigRational>>,|ctx,a,b|{let _=ctx;for q in &b{bounded_q(q)?;}a.with(ctx,|a|Ok(se::solve_dixon(a,&b)))});
object_binding!(exact_rational_sparse_solve,"num_linear_rational_sparse_solve","crate::kernels::sparse_exact::solve","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,b:Vec<BigRational>)->Dynamic,|ctx,a,b|{let _=ctx;for q in &b{bounded_q(q)?;}a.with(ctx,|a|{let result=se::solve(&mut se::Rationals,a,&b).ok_or_else(||err("exact RHS length differs from matrix rows"))?;let(tag,particular,basis)=match result{se::Solution::Inconsistent=>("Inconsistent",vec![],vec![]),se::Solution::Family{particular,nullspace}=>("Family",particular,nullspace)};let rows=basis.len();Ok(Dynamic(api::record(vec![("kind",Value::String(tag.into())),("particular",particular.encode()?),("basis_rows",rows.encode()?),("basis_values",basis.into_iter().flatten().collect::<Vec<_>>().encode()?)])))})});
object_binding!(exact_rational_sparse_echelon,"num_linear_rational_sparse_echelon","crate::kernels::sparse_exact::echelon","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::SparseMat<BigRational>>,pivot_cols:usize,record:bool)->NewObject<se::Echelon<BigRational>>,|ctx,a,pivot_cols,record|{let _=ctx;a.with(ctx,|a|{if pivot_cols>a.cols{return Err(err("pivot column count exceeds matrix shape"));}Ok(NewObject::normalized(se::echelon(&mut se::Rationals,a,pivot_cols,record)))})});
object_binding!(exact_echelon_rank,"num_linear_echelon_rank","crate::kernels::sparse_exact::Echelon<E>::rank","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::Echelon<BigRational>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.rank()))});
object_binding!(exact_echelon_odd_permutation,"num_linear_echelon_odd_permutation","crate::kernels::sparse_exact::Echelon<E>::odd_permutation","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<se::Echelon<BigRational>>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.odd_permutation()))});
object_binding!(exact_rational_reconstruct,"num_linear_rational_reconstruct","crate::kernels::sparse_exact::rational_reconstruct","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(residue:BigInt,modulus:BigInt)->Option<BigRational>,|ctx,residue,modulus|{let _=ctx;if residue.bits()>16384||modulus.bits()>16384||modulus<=BigInt::from(0){return Err(err("positive modulus and bounded integers required"));}Ok(se::rational_reconstruct(&residue,&modulus))});
object_binding!(exact_prime_sparse_new,"num_linear_prime_sparse_new","","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(p:u64,rows:usize,cols:usize,ri:Vec<usize>,ci:Vec<usize>,values:Vec<u64>)->NewObject<PrimeMatrix>,|ctx,p,rows,cols,ri,ci,values|{let _=ctx;let mut f=prime(p)?;if rows>256||cols>256||ri.len()!=ci.len()||ri.len()!=values.len()||values.len()>65536{return Err(err("prime matrix shape/allocation budget exceeded"));}let data=se::SparseMat::from_triplets(&mut f,rows,cols,ri.into_iter().zip(ci).zip(values).map(|((r,c),v)|(r,c,v%p))).map_err(err)?;Ok(NewObject::normalized(PrimeMatrix{data,field:f}))});
object_binding!(exact_berlekamp_massey,"num_linear_berlekamp_massey","crate::kernels::sparse_exact::berlekamp_massey","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(sequence:Vec<u64>,p:u64)->Vec<u64>,|ctx,sequence,p|{let _=ctx;let f=prime(p)?;if sequence.len()>4096||sequence.iter().any(|v|*v>=p){return Err(err("finite field sequence must have at most 4096 reduced residues"));}Ok(se::berlekamp_massey(&sequence,f))});
object_binding!(exact_minpoly_wiedemann,"num_linear_minpoly_wiedemann","crate::kernels::sparse_exact::minpoly_wiedemann","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<PrimeMatrix>,trials:usize,seed:u64)->Option<Vec<u64>>,|ctx,a,trials,seed|{let _=ctx;if trials==0||trials>32{return Err(err("Wiedemann trial budget exceeded"));}a.with(ctx,|a|Ok(se::minpoly_wiedemann(&a.data,a.field,trials,seed)))});
object_binding!(exact_solve_wiedemann,"num_linear_solve_wiedemann","crate::kernels::sparse_exact::solve_wiedemann","Native exact sparse/domain operation specialized to BigRational or a validated prime field. Coefficients are normalized dimensionless exact values; finite shape, coefficient-bit and work budgets are enforced.",true,(a:Handle<PrimeMatrix>,b:Vec<u64>,trials:usize,seed:u64)->Option<Vec<u64>>,|ctx,a,b,trials,seed|{let _=ctx;if trials==0||trials>32{return Err(err("Wiedemann trial budget exceeded"));}a.with(ctx,|a|{if b.iter().any(|v|*v>=a.field.p){return Err(err("finite field RHS residues must be reduced"));}Ok(se::solve_wiedemann(&a.data,&b,a.field,trials,seed))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")?;
    prepared::install(i)
}
static BINDINGS: &[&api::Binding] = &[
    &exact_rational_domain_zero::BINDING,
    &exact_rational_domain_one::BINDING,
    &exact_rational_domain_is_zero::BINDING,
    &exact_rational_domain_neg::BINDING,
    &exact_rational_domain_add::BINDING,
    &exact_rational_domain_sub::BINDING,
    &exact_rational_domain_mul::BINDING,
    &exact_rational_domain_div::BINDING,
    &exact_gfp_mul_mod::BINDING,
    &exact_gfp_pow_mod::BINDING,
    &exact_gfp_inv::BINDING,
    &exact_gfp_reduce::BINDING,
    &exact_word_primes::BINDING,
    &exact_rational_sparse_zeros::BINDING,
    &exact_rational_sparse_from_triplets::BINDING,
    &exact_rational_sparse_from_dense::BINDING,
    &exact_rational_sparse_to_dense::BINDING,
    &exact_rational_sparse_triplets::BINDING,
    &exact_rational_sparse_nnz::BINDING,
    &exact_rational_sparse_transpose::BINDING,
    &exact_rational_sparse_scale::BINDING,
    &exact_rational_sparse_add::BINDING,
    &exact_rational_sparse_matmul::BINDING,
    &exact_rational_sparse_matvec::BINDING,
    &exact_rational_sparse_det::BINDING,
    &exact_rational_sparse_rank::BINDING,
    &exact_rational_sparse_nullspace::BINDING,
    &exact_rational_sparse_det_modular::BINDING,
    &exact_rational_sparse_solve_dixon::BINDING,
    &exact_rational_sparse_solve::BINDING,
    &exact_rational_sparse_echelon::BINDING,
    &exact_echelon_rank::BINDING,
    &exact_echelon_odd_permutation::BINDING,
    &exact_rational_reconstruct::BINDING,
    &exact_prime_sparse_new::BINDING,
    &exact_berlekamp_massey::BINDING,
    &exact_minpoly_wiedemann::BINDING,
    &exact_solve_wiedemann::BINDING,
];
