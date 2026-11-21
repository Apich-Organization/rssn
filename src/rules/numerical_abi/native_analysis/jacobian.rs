//! Native normalized Jacobian/mass/factor objects; RHS objects retain physical frames.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::RecordValue;
use crate::kernels::dense::Mat;
use crate::kernels::jacobian as j;
use crate::kernels::krylov::Csr;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_analysis_ode_common::Budget;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeFrame;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeTape;
use crate::rules::numerical_abi::native_analysis_ode_common::RateFn;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::Maybe;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use std::sync::Arc;
#[derive(Clone)]
pub(in crate::rules) enum MatrixData {
    Dense(Mat),
    Sparse(Csr),
}
#[derive(Clone)]
pub(in crate::rules) struct MatrixObject {
    pub data: MatrixData,
    pub frame: Option<OdeFrame>,
    pub evaluations: usize,
}
impl NativeType for MatrixObject {
    const KIND: &'static str = "rssn/native/analysis/jacobian/matrix/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        let (rows, cols) = self.shape();
        if rows > 256 || cols > 256 {
            return Err(api::error(
                "Jacobian/mass matrix dimension budget exceeds256",
            ));
        }
        match &self.data {
            | MatrixData::Dense(m) => {
                if m.rows.checked_mul(m.cols) != Some(m.data.len())
                    || m.data.iter().any(|v| !v.is_finite())
                {
                    return Err(api::error("invalid normalized dense matrix"));
                }
            },
            | MatrixData::Sparse(c) => {
                if c.indptr.len() != c.rows + 1
                    || c.indptr.first() != Some(&0)
                    || c.indptr.last() != Some(&c.values.len())
                    || c.indices.len() != c.values.len()
                    || c.values.len() > 65536
                    || c.indptr.windows(2).any(|w| w[0] > w[1])
                    || c.indices.iter().any(|i| *i >= c.cols)
                    || c.values.iter().any(|v| !v.is_finite())
                {
                    return Err(api::error("invalid normalized sparse matrix"));
                }
            },
        }
        Ok(())
    }
}
impl MatrixObject {
    pub const fn shape(&self) -> (usize, usize) {
        match &self.data {
            | MatrixData::Dense(m) => (m.rows, m.cols),
            | MatrixData::Sparse(c) => (c.rows, c.cols),
        }
    }

    pub const fn mass(&self) -> j::Mass<'_> {
        match &self.data {
            | MatrixData::Dense(m) => j::Mass::Dense(m),
            | MatrixData::Sparse(c) => j::Mass::Sparse(c),
        }
    }

    pub fn jac(&self) -> j::Jac {
        match &self.data {
            | MatrixData::Dense(m) => j::Jac::Dense(m.clone()),
            | MatrixData::Sparse(c) => j::Jac::Sparse(c.clone()),
        }
    }
}
pub(in crate::rules) struct Triplets(pub Vec<(usize, usize, f64)>);
impl Arg for Triplets {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("rows".into(), <Vec<usize> as Arg>::ty()),
            ("cols".into(), <Vec<usize> as Arg>::ty()),
            ("values".into(), <Vec<f64> as Arg>::ty()),
        ])
    }

    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        let rows = <Vec<usize> as Arg>::decode(api::field(v, "rows")?, i)?;
        let cols = <Vec<usize> as Arg>::decode(api::field(v, "cols")?, i)?;
        let values = <Vec<f64> as Arg>::decode(api::field(v, "values")?, i)?;
        if rows.len() != cols.len() || rows.len() != values.len() || rows.len() > 65536 {
            return Err(api::error(
                "triplet arrays differ in length or exceed65536entries",
            ));
        }
        Ok(Self(
            rows.into_iter()
                .zip(cols)
                .zip(values)
                .map(|((i, j), v)| (i, j, v))
                .collect(),
        ))
    }
}
impl Output for Triplets {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let (mut rows, mut cols, mut values) = (Vec::new(), Vec::new(), Vec::new());
        for (i, j, v) in self.0 {
            rows.push(i);
            cols.push(j);
            values.push(v);
        }
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("rows".into(), Output::encode(rows)?),
                ("cols".into(), Output::encode(cols)?),
                ("values".into(), Output::encode(values)?),
            ],
        }))
    }
}
const fn matrix(data: MatrixData) -> NewObject<MatrixObject> {
    NewObject::normalized(MatrixObject {
        data,
        frame: None,
        evaluations: 0,
    })
}
fn triplet_bounds(
    trip: &Triplets,
    n: usize,
) -> Result<(), OpError> {
    if n == 0 || n > 512 || trip.0.iter().any(|(i, j, _)| *i >= n || *j >= n) {
        return Err(api::error(
            "factor dimensions/triplet indices exceed bounds",
        ));
    }
    Ok(())
}
object_binding!(dense_new,"num_analysis_jacobian_dense","","Retain a native normalized dense Jacobian/mass matrix. Measured quantities require explicit prior normalization; RHS-derived Jacobians retain their physical reference frame separately.",true,(data:Mat)->NewObject<MatrixObject>,|_ctx:&CallCtx<'_>,data:Mat|Ok(matrix(MatrixData::Dense(data))));
object_binding!(sparse_new,"num_analysis_jacobian_sparse","","Retain checked native normalized CSR data from exact row/column triplets.",true,(rows:usize,cols:usize,triplets:Triplets)->NewObject<MatrixObject>,|_ctx:&CallCtx<'_>,rows:usize,cols:usize,triplets:Triplets|{if rows>256||cols>256||triplets.0.iter().any(|(i,j,_)|*i>=rows||*j>=cols){return Err(api::error("sparse matrix dimensions or indices invalid"));}Ok(matrix(MatrixData::Sparse(Csr::from_triplets(rows,cols,&triplets.0))))});
object_binding!(band_pattern,"num_analysis_jacobian_band_pattern","crate::kernels::jacobian::band_pattern","Actual native normalized band sparsity pattern with checked dimension/bandwidth budgets.",true,(n:usize,lower:usize,upper:usize)->NewObject<MatrixObject>,|_ctx:&CallCtx<'_>,n:usize,lower:usize,upper:usize|{if n==0||n>256||lower>=n||upper>=n{return Err(api::error("band pattern needs n1..256 and bandwidths<n"));}Ok(matrix(MatrixData::Sparse(j::band_pattern(n,lower,upper))))});
struct Coloring {
    pattern: Csr,
    colors: (Vec<usize>, usize),
}
impl NativeType for Coloring {
    const KIND: &'static str = "rssn/native/analysis/jacobian/coloring/v1";
    const MAX_BYTES: usize = 2 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.pattern.cols > 256
            || self.colors.0.len() != self.pattern.cols
            || self.colors.1 > 256
        {
            return Err(api::error("native coloring shape mismatch"));
        }
        Ok(())
    }
}
object_binding!(color,"num_analysis_jacobian_color_columns","crate::kernels::jacobian::color_columns","Actual native graph coloring retained with its original sparsity pattern; prevents mismatched or forged colored finite-difference descriptors.",true,(pattern:Handle<MatrixObject>)->NewObject<Coloring>,|ctx:&CallCtx<'_>,pattern:Handle<MatrixObject>|pattern.with(ctx,|p|{let MatrixData::Sparse(c)=&p.data else{return Err(api::error("coloring requires sparse pattern"));};Ok(NewObject::normalized(Coloring{pattern:c.clone(),colors:j::color_columns(c)}))}));
struct BandFactor {
    native: j::BandLu,
    n: usize,
}
impl NativeType for BandFactor {
    const KIND: &'static str = "rssn/native/analysis/jacobian/band_factor/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.n > 512 || self.n == 0 {
            Err(api::error("invalid band factor size"))
        } else {
            Ok(())
        }
    }
}
struct Factor {
    native: j::Factor,
    n: usize,
}
impl NativeType for Factor {
    const KIND: &'static str = "rssn/native/analysis/jacobian/factor/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.n > 512 || self.n == 0 {
            Err(api::error("invalid factor size"))
        } else {
            Ok(())
        }
    }
}
object_binding!(band_factor,"num_analysis_band_lu_factor","crate::kernels::jacobian::BandLu::factor","Actual native band LU; exact triplet bounds and band widths are checked before allocation; singularity remains None.",true,(n:usize,lower:usize,upper:usize,triplets:Triplets)->Maybe<NewObject<BandFactor>>,|_ctx:&CallCtx<'_>,n:usize,lower:usize,upper:usize,triplets:Triplets|{triplet_bounds(&triplets,n)?;if lower>=n||upper>=n{return Err(api::error("bandwidth must be<n"));}Ok(Maybe(j::BandLu::factor(n,lower,upper,&triplets.0).map(|native|NewObject::normalized(BandFactor{native,n}))))});
object_binding!(band_solve,"num_analysis_band_lu_solve","crate::kernels::jacobian::BandLu::solve","Actual native band-factor solve in normalized coordinates.",true,(factor:Handle<BandFactor>,rhs:Vec<f64>)->Vec<f64>,|ctx:&CallCtx<'_>,factor:Handle<BandFactor>,rhs:Vec<f64>|factor.with(ctx,|f|{if rhs.len()!=f.n{return Err(api::error("band solve RHS size mismatch"));}Ok(f.native.solve(&rhs))}));
object_binding!(factor,"num_analysis_jacobian_factor","crate::kernels::jacobian::Factor::from_triplets","Actual native dense/banded/sparse stage-block factorization. n*stages<=512; solver Dense,Banded,Sparse; all triplets are checked before native indexing.",true,(n:usize,stages:usize,triplets:Triplets,solver:String,band:(usize,usize))->Maybe<NewObject<Factor>>,|_ctx:&CallCtx<'_>,n:usize,stages:usize,triplets:Triplets,solver:String,band:(usize,usize)|{let total=n.checked_mul(stages).filter(|n|*n>0&&*n<=512).ok_or_else(||api::error("stage-block factor dimension budget"))?;triplet_bounds(&triplets,total)?;if band.0>=n||band.1>=n{return Err(api::error("factor bandwidth must be<n"));}let kind=match solver.as_str(){"Dense"=>j::SolverKind::Dense,"Banded"=>j::SolverKind::Banded,"Sparse"=>j::SolverKind::Sparse,_=>return Err(api::error("factor solver must be Dense,Banded,Sparse"))};Ok(Maybe(j::Factor::from_triplets(n,stages,&triplets.0,kind,band).map(|native|NewObject::normalized(Factor{native,n:total}))))});
object_binding!(factor_solve,"num_analysis_jacobian_factor_solve","crate::kernels::jacobian::Factor::solve","Actual native factor solve in its stage-major normalized coordinate layout.",true,(factor:Handle<Factor>,rhs:Vec<f64>)->Vec<f64>,|ctx:&CallCtx<'_>,factor:Handle<Factor>,rhs:Vec<f64>|factor.with(ctx,|f|{if rhs.len()!=f.n{return Err(api::error("factor solve RHS size mismatch"));}Ok(f.native.solve(&rhs))}));

object_binding!(mass_rows,"num_analysis_mass_rows","crate::kernels::jacobian::Mass<'_>::rows","Actual native normalized mass-matrix rows.",true,(matrix:Handle<MatrixObject>)->usize,|ctx:&CallCtx<'_>,matrix:Handle<MatrixObject>|matrix.with(ctx,|m|Ok(m.mass().rows())));

object_binding!(mass_cols,"num_analysis_mass_cols","crate::kernels::jacobian::Mass<'_>::cols","Actual native normalized mass-matrix cols.",true,(matrix:Handle<MatrixObject>)->usize,|ctx:&CallCtx<'_>,matrix:Handle<MatrixObject>|matrix.with(ctx,|m|Ok(m.mass().cols())));

object_binding!(mass_band,"num_analysis_mass_band","crate::kernels::jacobian::Mass<'_>::band","Actual native normalized mass-matrix band.",true,(matrix:Handle<MatrixObject>)->(usize,usize),|ctx:&CallCtx<'_>,matrix:Handle<MatrixObject>|matrix.with(ctx,|m|Ok(m.mass().band())));

object_binding!(mass_matvec,"num_analysis_mass_matvec","crate::kernels::jacobian::Mass<'_>::matvec","Actual native mass-matrix action in explicitly normalized coordinates; vector length is checked before native indexing.",true,(matrix:Handle<MatrixObject>,x:Vec<f64>)->Vec<f64>,|ctx:&CallCtx<'_>,matrix:Handle<MatrixObject>,x:Vec<f64>|matrix.with(ctx,|m|{if x.len()!=m.shape().1{return Err(api::error("mass action vector length mismatch"));}Ok(m.mass().matvec(&x))}));
object_binding!(color_data,"num_analysis_jacobian_coloring_data","","Inspect the actual native color assignment/count retained with its sparsity pattern.",true,(coloring:Handle<Coloring>)->(Vec<usize>,usize),|ctx:&CallCtx<'_>,coloring:Handle<Coloring>|coloring.with(ctx,|c|Ok(c.colors.clone())));
object_binding!(matrix_data,"num_analysis_jacobian_matrix_data","","Inspect normalized native matrix data and actual evaluation count. RHS-derived Jacobians additionally expose their original time/state reference quantities; normalized coefficients are never mislabeled as one common physical matrix domain.",false,(matrix:Handle<MatrixObject>)->Dynamic,|ctx:&CallCtx<'_>,matrix:Handle<MatrixObject>|matrix.with(ctx,|m|{let normalized=match &m.data{MatrixData::Dense(d)=>Output::encode(d.clone())?,MatrixData::Sparse(_)=>{let mut trip=Vec::new();m.mass().for_each(|i,j,v|trip.push((i,j,v)));Output::encode(Triplets(trip))?}};let refs=if let Some(f)=&m.frame{let state=f.state_scales.iter().zip(&f.state_dimensions).map(|(s,d)|api::quantity(*s,d.clone())).collect::<Result<Vec<_>,_>>()?;api::record(vec![("time",api::quantity(f.time_scale,f.time_dimension.clone())?),("state",Value::Tensor(crate::kernel::value::TensorValue{shape:vec![state.len()],elem:api::dynamic_type(),data:state}))])}else{api::record(vec![])};Ok(Dynamic(api::record(vec![("normalized",normalized),("evaluations",Output::encode(m.evaluations)?),("references",refs)])))}));
struct Fd {
    native: j::FdMode<RateFn>,
}
impl NativeType for Fd {
    const KIND: &'static str = "rssn/native/analysis/jacobian/fd_mode/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
object_binding!(fd_serial,"num_analysis_fd_mode_serial","crate::kernels::jacobian::FdMode<F>::serial","Instantiate the real native serial finite-difference mode for owned compiled RHS functions.",true,()->NewObject<Fd>,|_ctx:&CallCtx<'_>|Ok(NewObject::normalized(Fd{native:j::FdMode::<RateFn>::serial()})));
object_binding!(fd_parallel,"num_analysis_fd_mode_parallel","crate::kernels::jacobian::FdMode<F>::parallel","Instantiate the real native parallel finite-difference mode for immutable owned RHS functions; no host callback is shared unsafely.",true,()->NewObject<Fd>,|_ctx:&CallCtx<'_>|Ok(NewObject::normalized(Fd{native:j::FdMode::<RateFn>::parallel()})));

#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn rhs_state(
    ctx: &CallCtx<'_>,
    rhs: Handle<OdeTape>,
    y: &[f64],
    f0: &[f64],
    budget: usize,
) -> Result<(Arc<OdeTape>, Arc<Budget>, RateFn), OpError> {
    let rhs = Arc::new(rhs.with(ctx, |r| Ok(r.clone()))?);
    if y.len() != rhs.frame.state_scales.len() || f0.len() != y.len() {
        return Err(api::error("finite-difference state/rate size mismatch"));
    }
    if budget
        .checked_mul(
            rhs.tape
                .insts
                .len()
                .saturating_add(rhs.tape.outputs.len())
                .max(1),
        )
        .is_none_or(|n| n > 100_000_000)
    {
        return Err(api::error(
            "finite-difference instruction work budget exceeded",
        ));
    }
    let budget = Budget::new(budget)?;
    let f = rhs.clone().function(budget.clone());
    Ok((rhs, budget, f))
}

object_binding!(call_fd_jacobian_dense,"num_analysis_fd_jacobian_dense","crate::kernels::jacobian::fd_jacobian_dense","Actual native fd_jacobian_dense on explicitly normalized time/state/rate data; RHS retains physical references and returned matrix reports them. Parallel mode calls the genuine Sync native implementation on an owned immutable tape.",true,(rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize|{let(rhs,budget,f)=rhs_state(ctx,rhs,&state,&initial_rate,budget)?;let data=MatrixData::Dense(j::fd_jacobian_dense(&f,time,&state,&initial_rate));budget.finish()?;Ok(NewObject::new(MatrixObject{data,frame:Some(rhs.frame.clone()),evaluations:budget.used()},rhs.frame.dimensions()))});

object_binding!(call_fd_jacobian_colored,"num_analysis_fd_jacobian_colored","crate::kernels::jacobian::fd_jacobian_colored","Actual native fd_jacobian_colored on explicitly normalized time/state/rate data; RHS retains physical references and returned matrix reports them. Parallel mode calls the genuine Sync native implementation on an owned immutable tape.",true,(rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize,pattern:Handle<MatrixObject>,colors:Handle<Coloring>)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize,pattern:Handle<MatrixObject>,colors:Handle<Coloring>|{let(rhs,budget,f)=rhs_state(ctx,rhs,&state,&initial_rate,budget)?;let pattern=pattern.with(ctx,|p|{let MatrixData::Sparse(c)=&p.data else{return Err(api::error("colored FD needs sparse pattern"));};Ok(c.clone())})?;if pattern.rows!=state.len()||pattern.cols!=state.len(){return Err(api::error("colored FD pattern must match the state square dimension"));}let colors=colors.with(ctx,|c|{if c.pattern.rows!=pattern.rows||c.pattern.cols!=pattern.cols||c.pattern.indptr!=pattern.indptr||c.pattern.indices!=pattern.indices{return Err(api::error("coloring belongs to a different sparsity pattern"));}Ok(c.colors.clone())})?;let mut nfev=0usize;let data=MatrixData::Sparse(j::fd_jacobian_colored(&f,time,&state,&initial_rate,&pattern,&colors,&mut nfev));budget.finish()?;Ok(NewObject::new(MatrixObject{data,frame:Some(rhs.frame.clone()),evaluations:budget.used()},rhs.frame.dimensions()))});

object_binding!(call_fd_jacobian_dense_par,"num_analysis_fd_jacobian_dense_par","crate::kernels::jacobian::fd_jacobian_dense_par","Actual native fd_jacobian_dense_par on explicitly normalized time/state/rate data; RHS retains physical references and returned matrix reports them. Parallel mode calls the genuine Sync native implementation on an owned immutable tape.",true,(rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize|{let(rhs,budget,f)=rhs_state(ctx,rhs,&state,&initial_rate,budget)?;let data=MatrixData::Dense(j::fd_jacobian_dense_par(&f,time,&state,&initial_rate));budget.finish()?;Ok(NewObject::new(MatrixObject{data,frame:Some(rhs.frame.clone()),evaluations:budget.used()},rhs.frame.dimensions()))});

object_binding!(call_fd_jacobian_colored_par,"num_analysis_fd_jacobian_colored_par","crate::kernels::jacobian::fd_jacobian_colored_par","Actual native fd_jacobian_colored_par on explicitly normalized time/state/rate data; RHS retains physical references and returned matrix reports them. Parallel mode calls the genuine Sync native implementation on an owned immutable tape.",true,(rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize,pattern:Handle<MatrixObject>,colors:Handle<Coloring>)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,rhs:Handle<OdeTape>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,budget:usize,pattern:Handle<MatrixObject>,colors:Handle<Coloring>|{let(rhs,budget,f)=rhs_state(ctx,rhs,&state,&initial_rate,budget)?;let pattern=pattern.with(ctx,|p|{let MatrixData::Sparse(c)=&p.data else{return Err(api::error("colored FD needs sparse pattern"));};Ok(c.clone())})?;if pattern.rows!=state.len()||pattern.cols!=state.len(){return Err(api::error("colored FD pattern must match the state square dimension"));}let colors=colors.with(ctx,|c|{if c.pattern.rows!=pattern.rows||c.pattern.cols!=pattern.cols||c.pattern.indptr!=pattern.indptr||c.pattern.indices!=pattern.indices{return Err(api::error("coloring belongs to a different sparsity pattern"));}Ok(c.colors.clone())})?;let mut nfev=0usize;let data=MatrixData::Sparse(j::fd_jacobian_colored_par(&f,time,&state,&initial_rate,&pattern,&colors,&mut nfev));budget.finish()?;Ok(NewObject::new(MatrixObject{data,frame:Some(rhs.frame.clone()),evaluations:budget.used()},rhs.frame.dimensions()))});
object_binding!(block,"num_analysis_jacobian_block_triplets","crate::kernels::jacobian::block_triplets","Actual native stage-block assembly in normalized coordinates. cm/cj are immutable two-index compiled coefficient programs with dimensionless input/output; mass is None or a normalized matrix handle. n*stages<=128 bounds output triplets before allocation.",false,(n:usize,stages:usize,cm:Input,cj:Input,mass:Input,jacobian:Handle<MatrixObject>,budget:usize)->Triplets,|ctx:&CallCtx<'_>,n:usize,stages:usize,cm:Input,cj:Input,mass:Input,jacobian:Handle<MatrixObject>,budget:usize|{
 if n==0||stages==0||stages>7||n.checked_mul(stages).is_none_or(|n|n>128){return Err(api::error("block assembly dimension/stage budget exceeded"));}let jacobian=jacobian.with(ctx,|m|Ok(m.clone()))?;if jacobian.shape()!=(n,n){return Err(api::error("block Jacobian shape mismatch"));}
 let cm=compiled_snapshot(ctx,&cm.0,budget)?;let cj=compiled_snapshot(ctx,&cj.0,budget)?;for f in[&cm,&cj]{if f.input_dimensions.len()!=2{return Err(api::error("block coefficient program must have two index inputs"));}for d in f.input_dimensions.iter().chain(std::iter::once(&f.output_dimension)){api::same_dimension(d,&crate::kernel::signature::DimExpr::Dimless)?;}}
 let mass=match &mass.0{Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()=>None,Value::Opaque(_)=>{let m=<Handle<MatrixObject> as crate::rules::numerical_abi::native_object_api::NativeArg>::decode(&mass.0,4,ctx)?;Some(m.with(ctx,|m|Ok(m.clone()))?)},_=>return Err(api::error("block mass must be None or normalized matrix handle"))};
if mass.as_ref().is_some_and(|m|m.shape()!=(n,n)){return Err(api::error("block mass shape mismatch"));}
 let out=j::block_triplets(n,stages,&|i,j|cm.sample(&[i as f64,j as f64]),&|i,j|cj.sample(&[i as f64,j as f64]),mass.as_ref().map(MatrixObject::mass),&jacobian.jac());cm.finish()?;cj.finish()?;Ok(Triplets(out))
});
#[derive(Clone)]
struct Provider {
    rhs: Arc<OdeTape>,
    structure: String,
    band: (usize, usize),
    pattern: Option<Csr>,
    analytic: usize,
}
impl NativeType for Provider {
    const KIND: &'static str = "rssn/native/analysis/jacobian/provider/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.rhs.validate()?;
        if self.analytic > 2 {
            return Err(api::error("invalid Jacobian source selector"));
        }
        Ok(())
    }
}
impl Provider {
    fn structure(&self) -> j::JacStructure<'_> {
        match self.structure.as_str() {
            | "Banded" => {
                j::JacStructure::Banded {
                    kl: self.band.0,
                    ku: self.band.1,
                }
            },
            | "Sparse" => {
                self.pattern
                    .as_ref()
                    .map_or(j::JacStructure::Dense, j::JacStructure::Sparse)
            },
            | _ => j::JacStructure::Dense,
        }
    }

    fn prepare<'a>(
        &'a self,
        dense: &'a crate::kernels::ode_stiff::JacobianFn<'a>,
        sparse: &'a j::SparseJacobianFn<'a>,
    ) -> Result<j::JacobianProvider<'a>, OpError> {
        j::JacobianProvider::new(
            self.rhs.frame.state_scales.len(),
            self.structure(),
            if self.analytic == 1 {
                Some(dense)
            } else {
                None
            },
            if self.analytic == 2 {
                Some(sparse)
            } else {
                None
            },
        )
        .ok_or_else(|| api::error("invalid native Jacobian provider structure"))
    }

    fn sparse_jac(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
    ) -> Csr {
        let dense = self.rhs.jacobian(budget, t, y);
        if let Some(pattern) = &self.pattern {
            let mut sparse = pattern.clone();
            for i in 0..sparse.rows {
                for k in sparse.indptr[i]..sparse.indptr[i + 1] {
                    sparse.values[k] = dense.at(i, sparse.indices[k]);
                }
            }
            sparse
        } else {
            let mut trip = Vec::new();
            j::Mass::Dense(&dense).for_each(|i, j, v| {
                if v != 0.0 {
                    trip.push((i, j, v));
                }
            });
            Csr::from_triplets(dense.rows, dense.cols, &trip)
        }
    }
}
object_binding!(provider_new,"num_analysis_jacobian_provider_new","crate::kernels::jacobian::JacobianProvider<'a>::new","Actual native JacobianProvider builder over an owned RHS tape. Structure Dense,Banded,Sparse; analytic0 finite differences,1 exact native tape dense AD,2 exact native tape sparse AD. Sparse structure requires its normalized matrix pattern. Provider borrows are reconstructed per invocation to retain safe owned lifetimes; every evaluation uses the real native provider.",false,(rhs:Handle<OdeTape>,structure:String,band:(usize,usize),pattern:Input,analytic:usize)->NewObject<Provider>,|ctx:&CallCtx<'_>,rhs:Handle<OdeTape>,structure:String,band:(usize,usize),pattern:Input,analytic:usize|{
 let rhs=Arc::new(rhs.with(ctx,|r|Ok(r.clone()))?);let n=rhs.frame.state_scales.len();if !["Dense","Banded","Sparse"].contains(&structure.as_str())||band.0>=n||band.1>=n||analytic>2{return Err(api::error("invalid Jacobian provider option"));}
 let pattern=if structure=="Sparse"{let h=<Handle<MatrixObject> as crate::rules::numerical_abi::native_object_api::NativeArg>::decode(&pattern.0,3,ctx)?;Some(h.with(ctx,|m|{if m.shape()!=(n,n){return Err(api::error("provider pattern must match square state dimension"));}let MatrixData::Sparse(p)=&m.data else{return Err(api::error("provider pattern must be sparse"));};Ok(p.clone())})?)}else{if !matches!(&pattern.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){return Err(api::error("non-sparse provider takes None pattern"));}None};
 let provider=Provider{rhs,structure,band,pattern,analytic};let budget=Budget::new(1)?;let dense=|t,y:&[f64]|provider.rhs.jacobian(&budget,t,y);let sparse=|t,y:&[f64]|provider.sparse_jac(&budget,t,y);provider.prepare(&dense,&sparse)?;Ok(NewObject::normalized(provider))
});
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn provider_eval(
    ctx: &CallCtx<'_>,
    provider: Handle<Provider>,
    time: f64,
    state: Vec<f64>,
    initial_rate: Vec<f64>,
    nfev: usize,
    budget: usize,
    fd: Option<Handle<Fd>>,
) -> Result<NewObject<MatrixObject>, OpError> {
    let provider = provider.with(ctx, |p| Ok(p.clone()))?;
    let n = provider.rhs.frame.state_scales.len();
    if state.len() != n || initial_rate.len() != n || nfev > 1_000_000 {
        return Err(api::error("provider state/rate/evaluation-count mismatch"));
    }
    if budget
        .checked_mul(provider.rhs.tape.insts.len().saturating_mul(n).max(1))
        .is_none_or(|n| n > 100_000_000)
    {
        return Err(api::error("provider instruction work budget exceeded"));
    }
    let budget = Budget::new(budget)?;
    let f = provider.rhs.clone().function(budget.clone());
    let dense = |t, y: &[f64]| provider.rhs.jacobian(&budget, t, y);
    let sparse = |t, y: &[f64]| provider.sparse_jac(&budget, t, y);
    let native = provider.prepare(&dense, &sparse)?;
    let mut updated = nfev;
    let out = if let Some(fd) = fd {
        fd.with(ctx, |fd| {
            Ok(native.eval_with(&f, time, &state, &initial_rate, &mut updated, &fd.native))
        })?
    } else {
        native.eval(&f, time, &state, &initial_rate, &mut updated)
    };
    budget.finish()?;
    let data = match out {
        | j::Jac::Dense(m) => MatrixData::Dense(m),
        | j::Jac::Sparse(c) => MatrixData::Sparse(c),
    };
    Ok(NewObject::new(
        MatrixObject {
            data,
            frame: Some(provider.rhs.frame.clone()),
            evaluations: updated,
        },
        provider.rhs.frame.dimensions(),
    ))
}
object_binding!(provider_evaluate,"num_analysis_jacobian_provider_eval","crate::kernels::jacobian::JacobianProvider<'a>::eval","Actual native provider evaluation on explicitly normalized time/state/rate data; updated caller evaluation count and physical frame remain with the returned matrix.",true,(provider:Handle<Provider>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,nfev:usize,budget:usize)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,provider:Handle<Provider>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,nfev:usize,budget:usize|provider_eval(ctx,provider,time,state,initial_rate,nfev,budget,None));
object_binding!(provider_evaluate_with,"num_analysis_jacobian_provider_eval_with","crate::kernels::jacobian::JacobianProvider<'a>::eval_with","Actual native provider evaluation with explicit real serial/parallel FD mode; normalized frame contracts and actual updated nfev are preserved.",true,(provider:Handle<Provider>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,nfev:usize,budget:usize,fd:Handle<Fd>)->NewObject<MatrixObject>,|ctx:&CallCtx<'_>,provider:Handle<Provider>,time:f64,state:Vec<f64>,initial_rate:Vec<f64>,nfev:usize,budget:usize,fd:Handle<Fd>|provider_eval(ctx,provider,time,state,initial_rate,nfev,budget,Some(fd)));

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &dense_new::BINDING,
            &sparse_new::BINDING,
            &band_pattern::BINDING,
            &color::BINDING,
            &band_factor::BINDING,
            &band_solve::BINDING,
            &factor::BINDING,
            &factor_solve::BINDING,
            &mass_rows::BINDING,
            &mass_cols::BINDING,
            &mass_band::BINDING,
            &mass_matvec::BINDING,
            &color_data::BINDING,
            &matrix_data::BINDING,
            &fd_serial::BINDING,
            &fd_parallel::BINDING,
            &call_fd_jacobian_dense::BINDING,
            &call_fd_jacobian_colored::BINDING,
            &call_fd_jacobian_dense_par::BINDING,
            &call_fd_jacobian_colored_par::BINDING,
            &block::BINDING,
            &provider_new::BINDING,
            &provider_evaluate::BINDING,
            &provider_evaluate_with::BINDING,
        ],
        "numerical_vector",
    )
}
