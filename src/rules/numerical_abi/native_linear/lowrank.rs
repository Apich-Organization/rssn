//! Native low-rank factor and matrix-free operator specializations.
use super::super::native_backend::compiled_snapshot;
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
use super::TypeExpr;
use super::Value;
use super::api;
use super::array;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::scale;
use crate::kernels::dense::Mat;
use crate::kernels::lowrank;
impl NativeType for lowrank::Gaussian {
    const KIND: &'static str = "rssn/native/gaussian/v1";
    const MAX_BYTES: usize = 256;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
struct Pivot {
    data: lowrank::PivotedQr,
    scale: f64,
}
impl NativeType for Pivot {
    const KIND: &'static str = "rssn/native/pivoted-qr/v1";
    const MAX_BYTES: usize = 2 * LIMIT * 8 + 2048;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.q.rows > 128 || self.data.r.cols > 128 {
            return Err(err("pivoted QR work budget exceeded"));
        }
        Ok(())
    }
}
pub(super) fn matrix(a: &Array) -> Result<Mat, OpError> {
    if a.data.ndim() != 2 || a.data.shape().iter().any(|n| *n > 128) {
        return Err(err("dense operator shape/work budget exceeded"));
    }
    Ok(Mat {
        rows: a.data.shape()[0],
        cols: a.data.shape()[1],
        data: a.data.iter().copied().collect(),
    })
}
fn sresult(
    out: crate::kernels::dense::Svd,
    reference: f64,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("u", out.u.encode()?),
        ("s", out.s.encode()?),
        ("v", out.v.encode()?),
        ("singular_value_reference", api::quantity(reference, d)?),
    ])))
}

impl super::super::native_object_api::NativeOutput for lowrank::ColumnInterpolation {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("columns".into(), <Vec<usize> as Output>::ty()),
            ("coefficients".into(), <Mat as Output>::ty()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("columns", Output::encode(self.columns)?),
            ("coefficients", Output::encode(self.coefficients)?),
        ]))
    }
}

object_binding!(
    lowrank_interpolative_decomposition,
    "num_linear_interpolative_decomposition",
    "crate::kernels::lowrank::interpolative_decomposition",
    "Approximates a matrix using a subset of its original columns and dimensionless interpolation coefficients. Returns zero-based column indices and a coefficient matrix in the original column order.",
    true,
    (a: Handle<Array>, max_rank: usize) -> lowrank::ColumnInterpolation,
    |ctx: &CallCtx<'_>, a: Handle<Array>, max_rank: usize| {
        a.with(ctx, |a| lowrank::interpolative_decomposition(&matrix(a)?, max_rank).map_err(err))
    }
);
object_binding!(lowrank_gaussian_new,"num_linear_gaussian_new","crate::kernels::lowrank::Gaussian::new","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(seed:u64)->NewObject<lowrank::Gaussian>,|ctx,seed|{let _=ctx;Ok(NewObject::normalized(lowrank::Gaussian::new(seed)))});
object_binding!(lowrank_gaussian_sample,"num_linear_gaussian_sample","crate::kernels::lowrank::Gaussian::sample","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(generator:Handle<lowrank::Gaussian>)->f64,|ctx,generator|{let _=ctx;generator.with_mut(ctx,|g|Ok(g.sample()))});
object_binding!(lowrank_pivoted_qr,"num_linear_pivoted_qr","crate::kernels::lowrank::qr_pivoted","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(a:Handle<Array>)->NewObject<Pivot>,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(NewObject::new(Pivot{data:lowrank::qr_pivoted(&matrix(a)?).map_err(err)?,scale:a.scale},vec![d])))});
object_binding!(lowrank_pivoted_qr_rank,"num_linear_pivoted_qr_rank","crate::kernels::lowrank::PivotedQr::rank","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(a:Handle<Pivot>,tolerance:Option<f64>)->usize,|ctx,a,tolerance|{let _=ctx;if tolerance.is_some_and(|t|t<0.0){return Err(err("rank tolerance must be nonnegative"));}a.with(ctx,|a|Ok(a.data.rank(tolerance)))});
object_binding!(lowrank_pivoted_qr_solve,"num_linear_pivoted_qr_solve","crate::kernels::lowrank::PivotedQr::solve","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(a:Handle<Pivot>,rhs:Handle<Array>,tolerance:Option<f64>)->NewObject<Array>,|ctx,a,rhs,tolerance|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,crate::rules::numerical_abi::runtime_dimensions::expression(&crate::rules::numerical_abi::runtime_dimensions::canonical(&DimExpr::Inverse(Box::new(dim(ctx,&a.value)?)))?)?)?;if tolerance.is_some_and(|t|t<0.0){return Err(err("rank tolerance must be nonnegative"));}a.with(ctx,|a|rhs.with(ctx,|b|{if b.data.ndim()!=1{return Err(err("RHS must be a vector"));}let out=a.data.solve(b.data.as_slice().ok_or_else(||err("RHS not contiguous"))?,tolerance).map_err(err)?;Ok(array(ArrayD::from_shape_vec(IxDyn(&[out.len()]),out).map_err(err)?,scale(b.scale,1.0/a.scale)?,d))}))});
object_binding!(
    lowrank_randomized_svd,
    "num_linear_randomized_svd",
    "crate::kernels::lowrank::randomized_svd",
    "Approximates the leading singular triplets of a finite matrix with a seeded Gaussian sketch. Oversampling widens the sketch and power iterations improve separation of nearby singular values. The columns of u and v are dimensionless. Multiply each entry of s by singular_value_reference to obtain a singular value in the matrix's physical units.",
    true,
    (a: Handle<Array>, rank: usize, oversample: usize, power_iterations: usize, seed: u64) -> Dynamic,
    |ctx, a, rank, oversample, power_iterations, seed| {
        let d = dim(ctx, &a.value)?;
        a.with(ctx, |a| {
            let m = matrix(a)?;
            if rank == 0 || rank > 64 || oversample > 64 || power_iterations > 16 {
                return Err(err("randomized SVD work budget exceeded"));
            }
            let mt = m.transpose();
            let apply = |x: &[f64], y: &mut [f64]| m.matvec_into(x, y);
            let transpose = |x: &[f64], y: &mut [f64]| mt.matvec_into(x, y);
            let options = lowrank::RandomizedOptions { oversample, power_iterations, seed };
            let out = lowrank::randomized_svd(&apply, &transpose, m.rows, m.cols, rank, &options)
                .map_err(err)?;
            sresult(out, a.scale, d)
        })
    }
);
object_binding!(lowrank_lobpcg,"num_linear_lobpcg","crate::kernels::lowrank::lobpcg","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(a:Handle<Array>,rank:usize,tolerance:f64,iterations:usize,seed:u64)->Dynamic,|ctx,a,rank,tolerance,iterations,seed|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{let m=matrix(a)?;if m.rows!=m.cols||rank==0||iterations==0||iterations>1024||rank>32||tolerance<=0.0{return Err(err("LOBPCG shape/iteration budget exceeded"));}let apply=|x:&[f64],y:&mut[f64]|y.copy_from_slice(&m.matvec(x));let out=lowrank::lobpcg(&apply,None,m.rows,rank,tolerance,iterations,seed).map_err(err)?;Ok(Dynamic(api::record(vec![("values",out.values.encode()?),("vectors",out.vectors.encode()?),("residuals",out.residuals.encode()?),("iterations",out.iterations.encode()?),("eigenvalue_reference",api::quantity(a.scale,d)?)])))})});
object_binding!(lowrank_anderson,"num_linear_anderson","crate::kernels::lowrank::anderson","Actual native low-rank/operator method. Matrices and state arrays retain explicit physical references; orthonormal factors are dimensionless. Shapes, iterations, callbacks and allocations are bounded.",true,(start:Handle<Array>,programs:Input,window:usize,damping:f64,tolerance:f64,max_evaluations:usize)->Dynamic,|ctx,start,programs,window,damping,tolerance,max_evaluations|{let _=ctx;let d=dim(ctx,&start.value)?;start.with(ctx,|a|{if a.data.ndim()!=1||a.data.len()>32||window>32||max_evaluations==0||max_evaluations>4096||tolerance<=0.0||damping<=0.0||damping>1.0{return Err(err("Anderson shape/options/work budget exceeded"));}let Value::Tensor(t)=&programs.0 else{return Err(err("expected vector of compiled coordinate maps"));};
if t.shape!=[a.data.len()]{return Err(err("one compiled scalar output is required per state coordinate"));}let mut callbacks=Vec::new();for p in &t.data{let s=compiled_snapshot(ctx,p,max_evaluations)?;if s.input_dimensions.len()!=a.data.len(){return Err(err("Anderson map arity differs from state size"));}for input in &s.input_dimensions{api::same_dimension(input,&d)?;}api::same_dimension(&s.output_dimension,&d)?;callbacks.push(s);}let map=|x:&[f64]|{let coherent=x.iter().map(|x|x*a.scale).collect::<Vec<_>>();callbacks.iter().map(|c|c.sample(&coherent)/a.scale).collect::<Vec<_>>()};let out=lowrank::anderson(&map,a.data.as_slice().ok_or_else(||err("state not contiguous"))?,window,damping,tolerance,max_evaluations);for c in callbacks{c.finish()?;}let out=out.map_err(err)?;Ok(Dynamic(api::record(vec![("point",out.x.encode()?),("reference",api::quantity(a.scale,d.clone())?),("residual",api::quantity(out.residual*a.scale,d)?),("evaluations",out.evaluations.encode()?)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &lowrank_gaussian_new::BINDING,
    &lowrank_gaussian_sample::BINDING,
    &lowrank_pivoted_qr::BINDING,
    &lowrank_pivoted_qr_rank::BINDING,
    &lowrank_pivoted_qr_solve::BINDING,
    &lowrank_interpolative_decomposition::BINDING,
    &lowrank_randomized_svd::BINDING,
    &lowrank_lobpcg::BINDING,
    &lowrank_anderson::BINDING,
];
