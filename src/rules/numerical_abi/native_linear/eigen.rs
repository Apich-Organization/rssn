//! Native matrix-free eigensolver monomorphizations over checked dense operators.
use super::Array;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Normalized;
use super::OpError;
use super::Output;
use super::RuleError;
use super::WORK;
use super::api;
use super::dim;
use super::err;
use super::lowrank::matrix;
use super::object_binding;
use super::product;
use super::scale;
use crate::kernels::dense;
use crate::kernels::eigs;
use crate::kernels::eigs_ext as ex;
fn opts(
    n: usize,
    k: usize,
    tolerance: f64,
    ncv: usize,
    restarts: usize,
    seed: u64,
) -> Result<eigs::EigsOptions, OpError> {
    if n == 0
        || n > 128
        || k == 0
        || k > n
        || tolerance <= 0.0
        || ncv > 128
        || restarts == 0
        || restarts > 512
        || n.checked_mul(n)
            .and_then(|w| w.checked_mul(restarts))
            .is_none_or(|w| w > WORK)
    {
        return Err(err("eigensolver dimensions/tolerance/work budget exceeded"));
    }
    Ok(eigs::EigsOptions {
        ncv,
        tol: tolerance,
        max_restarts: restarts,
        seed,
    })
}
fn which(name: &str) -> Result<eigs::Which, OpError> {
    Ok(match name {
        | "LargestMagnitude" => eigs::Which::LargestMagnitude,
        | "SmallestMagnitude" => eigs::Which::SmallestMagnitude,
        | "LargestReal" => eigs::Which::LargestReal,
        | "SmallestReal" => eigs::Which::SmallestReal,
        | "LargestImag" => eigs::Which::LargestImag,
        | _ => return Err(err("unknown eigenspectrum selection")),
    })
}
fn output(
    out: eigs::EigsResult,
    reference: f64,
    dim: DimExpr,
    residual_ref: f64,
    residual_dim: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("values", out.values.encode()?),
        ("vectors_re", out.vectors_re.encode()?),
        ("vectors_im", out.vectors_im.encode()?),
        ("residuals", out.residuals.encode()?),
        ("converged", out.converged.encode()?),
        ("restarts", out.restarts.encode()?),
        ("matvecs", out.matvecs.encode()?),
        ("eigenvalue_reference", api::quantity(reference, dim)?),
        (
            "residual_reference",
            api::quantity(residual_ref, residual_dim)?,
        ),
    ])))
}
fn square(a: &Array) -> Result<dense::Mat, OpError> {
    let m = matrix(a)?;
    if m.rows != m.cols || m.rows == 0 {
        return Err(err("eigensolver operator must be nonempty square"));
    }
    Ok(m)
}
fn inverse_dim(d: DimExpr) -> Result<DimExpr, OpError> {
    crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Inverse(Box::new(
            d,
        )))?,
    )
}
#[allow(
    clippy::too_many_arguments,
    reason = "checks the independent native Jacobi-Davidson controls together"
)]
fn jd(
    n: usize,
    tolerance: f64,
    iterations: usize,
    inner: usize,
    min_dim: usize,
    max_dim: usize,
    harmonic: bool,
    seed: u64,
) -> Result<ex::JdOptions, OpError> {
    if tolerance <= 0.0
        || iterations == 0
        || iterations > 512
        || inner == 0
        || inner > 128
        || min_dim == 0
        || min_dim >= max_dim
        || max_dim > n
        || n.checked_mul(iterations)
            .and_then(|w| w.checked_mul(inner))
            .is_none_or(|w| w > WORK)
    {
        return Err(err("Jacobi-Davidson control/work budget exceeded"));
    }
    Ok(ex::JdOptions {
        tol: tolerance,
        max_iter: iterations,
        min_dim,
        max_dim,
        inner_iters: inner,
        inner_tol: 1e-3,
        harmonic,
        seed,
    })
}
object_binding!(eigen_arnoldi,"num_linear_arnoldi","crate::kernels::eigs::arnoldi","Builds a dimensionless Krylov basis and the Hessenberg matrix for the normalized linear operator. The result contains basis, hessenberg, beta and the measured operator_reference. Multiply the Hessenberg coefficients and residual beta by that reference to restore physical units.",true,(a:Handle<Array>,start:Vec<f64>,space:usize)->Dynamic,|ctx,a,start,space|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;if start.len()!=m.rows||space==0||space>m.rows{return Err(err("Arnoldi start/space dimension invalid"));}let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let(v,h,beta)=eigs::arnoldi(apply,&start,space);Ok(Dynamic(api::record(vec![("basis",v.encode()?),("hessenberg",h.encode()?),("beta",beta.encode()?),("operator_reference",api::quantity(a.scale,d)?)])))})});
object_binding!(eigen_eigs,"num_linear_eigs","crate::kernels::eigs::eigs","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,selection:String,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,k,selection,tolerance,ncv,restarts,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let out=eigs::eigs(apply,m.rows,k,which(&selection)?,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_sym,"num_linear_eigs_sym","crate::kernels::eigs::eigs_sym","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,selection:String,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,k,selection,tolerance,ncv,restarts,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let out=eigs::eigs_sym(apply,m.rows,k,which(&selection)?,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_krylov_schur,"num_linear_eigs_krylov_schur","crate::kernels::eigs_ext::eigs_krylov_schur","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,selection:String,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,k,selection,tolerance,ncv,restarts,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let out=ex::eigs_krylov_schur(apply,m.rows,k,which(&selection)?,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_sym_block,"num_linear_eigs_sym_block","crate::kernels::eigs_ext::eigs_sym_block","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,selection:String,tolerance:f64,ncv:usize,restarts:usize,seed:u64,block:usize)->Dynamic,|ctx,a,k,selection,tolerance,ncv,restarts,seed,block|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;if block==0||block>k{return Err(err("Lanczos block size invalid"));}let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let out=ex::eigs_sym_block(apply,m.rows,k,which(&selection)?,block,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_shift_invert,"num_linear_eigs_shift_invert","crate::kernels::eigs::eigs_shift_invert","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,sigma:f64,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,k,sigma,tolerance,ncv,restarts,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let mut m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;let n=m.rows;for i in 0..n{m.set(i,i,m.at(i,i)-sigma);}let lu=dense::lu_factor(&m).map_err(err)?;let solve=|x:&[f64],y:&mut[f64]|y.copy_from_slice(&lu.solve(x));let out=eigs::eigs_shift_invert(solve,n,k,sigma,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_sym_shift_invert,"num_linear_eigs_sym_shift_invert","crate::kernels::eigs::eigs_sym_shift_invert","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,sigma:f64,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,k,sigma,tolerance,ncv,restarts,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let mut m=square(a)?;let options=opts(m.rows,k,tolerance,ncv,restarts,seed)?;let n=m.rows;for i in 0..n{m.set(i,i,m.at(i,i)-sigma);}let lu=dense::lu_factor(&m).map_err(err)?;let solve=|x:&[f64],y:&mut[f64]|y.copy_from_slice(&lu.solve(x));let out=eigs::eigs_sym_shift_invert(solve,n,k,sigma,&options).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_b_is_positive_definite,"num_linear_b_is_positive_definite","crate::kernels::eigs_ext::b_is_positive_definite","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(b:Handle<Array>,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->bool,|ctx,b,tolerance,ncv,restarts,seed|{let _=ctx;b.with(ctx,|b|{let m=square(b)?;let options=opts(m.rows,1,tolerance,ncv,restarts,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);Ok(ex::b_is_positive_definite(&apply,m.rows,&options))})});
object_binding!(eigen_eigs_gen_shift_invert,"num_linear_eigs_gen_shift_invert","crate::kernels::eigs_ext::eigs_gen_shift_invert","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,b:Handle<Array>,k:usize,sigma:f64,block:usize,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,b,k,sigma,block,tolerance,ncv,restarts,seed|{let _=ctx;let da=dim(ctx,&a.value)?;let d=product(da.clone(),inverse_dim(dim(ctx,&b.value)?)?)?;a.with(ctx,|a|b.with(ctx,|b|{let am=square(a)?;let bm=square(b)?;if am.rows!=bm.rows||block==0||block>k{return Err(err("generalized eigensystem shapes/block invalid"));}let options=ex::GenEigsOptions{sigma,block,eigs:opts(am.rows,k,tolerance,ncv,restarts,seed)?};let mut shifted=am.clone();for i in 0..shifted.data.len(){shifted.data[i]-=sigma*bm.data[i];}let lu=dense::lu_factor(&shifted).map_err(err)?;let apply_a=|x:&[f64],y:&mut[f64]|am.matvec_into(x, y);let apply_b=|x:&[f64],y:&mut[f64]|bm.matvec_into(x, y);let solve=|x:&[f64],y:&mut[f64]|y.copy_from_slice(&lu.solve(x));let out=ex::eigs_gen_shift_invert(apply_a,apply_b,solve,am.rows,k,&options).map_err(|e|err(format!("{e:?}")))?;output(out,scale(a.scale,1.0/b.scale)?,d,a.scale,da)}))});
object_binding!(eigen_eigs_gen_sym_shift_invert,"num_linear_eigs_gen_sym_shift_invert","crate::kernels::eigs_ext::eigs_gen_sym_shift_invert","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,b:Handle<Array>,k:usize,sigma:f64,block:usize,tolerance:f64,ncv:usize,restarts:usize,seed:u64)->Dynamic,|ctx,a,b,k,sigma,block,tolerance,ncv,restarts,seed|{let _=ctx;let da=dim(ctx,&a.value)?;let d=product(da.clone(),inverse_dim(dim(ctx,&b.value)?)?)?;a.with(ctx,|a|b.with(ctx,|b|{let am=square(a)?;let bm=square(b)?;if am.rows!=bm.rows||block==0||block>k{return Err(err("generalized eigensystem shapes/block invalid"));}let options=ex::GenEigsOptions{sigma,block,eigs:opts(am.rows,k,tolerance,ncv,restarts,seed)?};let mut shifted=am.clone();for i in 0..shifted.data.len(){shifted.data[i]-=sigma*bm.data[i];}let lu=dense::lu_factor(&shifted).map_err(err)?;let apply_a=|x:&[f64],y:&mut[f64]|am.matvec_into(x, y);let apply_b=|x:&[f64],y:&mut[f64]|bm.matvec_into(x, y);let solve=|x:&[f64],y:&mut[f64]|y.copy_from_slice(&lu.solve(x));let out=ex::eigs_gen_sym_shift_invert(apply_a,apply_b,solve,am.rows,k,&options).map_err(|e|err(format!("{e:?}")))?;output(out,scale(a.scale,1.0/b.scale)?,d,a.scale,da)}))});
object_binding!(eigen_eigs_jd,"num_linear_eigs_jd","crate::kernels::eigs_ext::eigs_jd","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,target_real:f64,target_imag:f64,tolerance:f64,iterations:usize,inner:usize,min_dim:usize,max_dim:usize,harmonic:bool,seed:u64)->Dynamic,|ctx,a,k,target_real,target_imag,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;if k==0||k>m.rows{return Err(err("invalid requested eigenpair count"));}let options=jd(m.rows,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let out=ex::eigs_jd(apply,m.rows,k,(target_real,target_imag),&options,None).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_jd_sym,"num_linear_eigs_jd_sym","crate::kernels::eigs_ext::eigs_jd_sym","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,target_real:f64,target_imag:f64,tolerance:f64,iterations:usize,inner:usize,min_dim:usize,max_dim:usize,harmonic:bool,seed:u64)->Dynamic,|ctx,a,k,target_real,target_imag,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;if k==0||k>m.rows{return Err(err("invalid requested eigenpair count"));}let options=jd(m.rows,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);if !crate::kernels::exact_eq(target_imag,0.0){return Err(err("symmetric Jacobi-Davidson target must be real"));}let out=ex::eigs_jd_sym(apply,m.rows,k,target_real,&options,None).map_err(|e|err(format!("{e:?}")))?;output(out,a.scale,d.clone(),a.scale,d)})});
object_binding!(eigen_eigs_jd_two_sided,"num_linear_eigs_jd_two_sided","crate::kernels::eigs_ext::eigs_jd_two_sided","Native eigensolver over a checked operator matrix, retaining explicit entry reference scales. Eigenvectors are normalized native rays; eigenvalue and residual references are returned separately (generalized eigenvalues use A/B). Native algorithm identity, convergence reports and work limits are preserved.",true,(a:Handle<Array>,k:usize,target_real:f64,target_imag:f64,tolerance:f64,iterations:usize,inner:usize,min_dim:usize,max_dim:usize,harmonic:bool,seed:u64)->Dynamic,|ctx,a,k,target_real,target_imag,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=square(a)?;if k==0||k>m.rows{return Err(err("invalid requested eigenpair count"));}let options=jd(m.rows,tolerance,iterations,inner,min_dim,max_dim,harmonic,seed)?;let apply=|x:&[f64],y:&mut[f64]|m.matvec_into(x, y);let mt=m.transpose();let transpose=|x:&[f64],y:&mut[f64]|mt.matvec_into(x, y);let out=ex::eigs_jd_two_sided(apply,transpose,m.rows,k,(target_real,target_imag),&options).map_err(|e|err(format!("{e:?}")))?;Ok(Dynamic(api::record(vec![("values",out.values.encode()?),("right_re",out.right_re.encode()?),("right_im",out.right_im.encode()?),("left_re",out.left_re.encode()?),("left_im",out.left_im.encode()?),("right_residuals",out.right_residuals.encode()?),("left_residuals",out.left_residuals.encode()?),("condition",out.condition.encode()?),("converged",out.converged.encode()?),("iterations",out.iterations.encode()?),("matvecs",out.matvecs.encode()?),("operator_reference",api::quantity(a.scale,d)?)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &complex_schur::BINDING,
    &reorder_schur::BINDING,
    &eigen_arnoldi::BINDING,
    &eigen_eigs::BINDING,
    &eigen_eigs_sym::BINDING,
    &eigen_eigs_krylov_schur::BINDING,
    &eigen_eigs_sym_block::BINDING,
    &eigen_eigs_shift_invert::BINDING,
    &eigen_eigs_sym_shift_invert::BINDING,
    &eigen_b_is_positive_definite::BINDING,
    &eigen_eigs_gen_shift_invert::BINDING,
    &eigen_eigs_gen_sym_shift_invert::BINDING,
    &eigen_eigs_jd::BINDING,
    &eigen_eigs_jd_sym::BINDING,
    &eigen_eigs_jd_two_sided::BINDING,
];
object_binding!(complex_schur,"num_linear_complex_schur","crate::kernels::eigs_ext::complex_schur","Native complex Schur decomposition of normalized coefficients with a positive physical entry reference; q is dimensionless and t uses the retained entry reference.",false,(a:Normalized<Vec<Vec<num_complex::Complex64>>>,reference:Input)->Dynamic,|ctx,a,reference|{let _=ctx;let n=a.0.len();if n>64||a.0.iter().any(|r|r.len()!=n){return Err(err("complex Schur needs bounded square matrix"));}let(s,d)=api::reference(&reference.0)?;let out=ex::complex_schur(&a.0).ok_or_else(||err("native complex Schur iteration did not converge"))?;Ok(Dynamic(api::record(vec![("t",out.0.encode()?),("q",out.1.encode()?),("reference",api::quantity(s,d)?)])))});
object_binding!(reorder_schur,"num_linear_reorder_schur","crate::kernels::eigs_ext::reorder_schur","Native reordered complex Schur form; a retained compiled real key accepts coherent real/imaginary eigenvalue parts and returns a dimensionless ordering key. Both returned matrices retain native layouts and the t entry reference.",false,(t:Normalized<Vec<Vec<num_complex::Complex64>>>,q:Normalized<Vec<Vec<num_complex::Complex64>>>,count:usize,program:Input,reference:Input,budget:usize)->Dynamic,|ctx,t,q,count,program,reference,budget|{let n=t.0.len();if n>64||count>n||q.0.len()!=n||t.0.iter().chain(&q.0).any(|r|r.len()!=n){return Err(err("Schur reorder shape/count budget invalid"));}let(s,d)=api::reference(&reference.0)?;let key=crate::rules::numerical_abi::native_backend::compiled_snapshot(ctx,&program.0,budget)?;if key.input_dimensions.len()!=2{return Err(err("complex eigenvalue key needs two inputs"));}for input in &key.input_dimensions{api::same_dimension(input,&d)?;}api::same_dimension(&key.output_dimension,&DimExpr::Dimless)?;let keyfn=|z:num_complex::Complex64|key.sample(&[z.re*s,z.im*s]);let(mut t,mut q)=(t.0,q.0);ex::reorder_schur(&mut t,&mut q,count,&keyfn);key.finish()?;Ok(Dynamic(api::record(vec![("t",t.encode()?),("q",q.encode()?),("reference",api::quantity(s,d)?)])))});
