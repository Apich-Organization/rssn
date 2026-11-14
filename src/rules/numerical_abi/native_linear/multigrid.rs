//! Native geometric and algebraic multigrid hierarchies with operator reference domains.
use super::super::native_object_api::Maybe;
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
use super::Output;
use super::RuleError;
use super::api;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::scale;
use super::sparse::Csr;
use super::sparse::inverse;
use super::sparse::result;
use super::sparse::vector;
use crate::kernels::amg_general as ag;
use crate::kernels::multigrid as mg;
struct Amg<T> {
    data: T,
    n: usize,
    scale: f64,
}
macro_rules! hierarchy {
    ($t:path,$kind:literal) => {
        impl NativeType for Amg<$t> {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = 16 * 1024 * 1024;

            fn validate(&self) -> Result<(), OpError> {
                if self.n == 0 || self.n > 512 || !self.scale.is_finite() || self.scale <= 0.0 {
                    return Err(err("AMG hierarchy size/reference budget exceeded"));
                }
                Ok(())
            }
        }
    };
}
hierarchy!(mg::Amg, "rssn/native/amg/v1");
hierarchy!(ag::AmgGeneral, "rssn/native/amg-general/v1");
struct Geometric {
    data: mg::GeometricMg,
    n: usize,
    dimension: usize,
    scale: f64,
}
impl NativeType for Geometric {
    const KIND: &'static str = "rssn/native/geometric-mg/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.size() > 32768 || !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(err("geometric hierarchy size/reference budget exceeded"));
        }
        Ok(())
    }
}
struct MgResult {
    data: mg::MgResult,
    scale: f64,
}
impl NativeType for MgResult {
    const KIND: &'static str = "rssn/native/multigrid-result/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 65536;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.x.len() > LIMIT
            || self.data.history.len() > 2049
            || self
                .data
                .x
                .iter()
                .chain(&self.data.history)
                .any(|v| !v.is_finite())
            || !self.scale.is_finite()
            || self.scale <= 0.0
        {
            return Err(err("invalid finite multigrid result or work budget"));
        }
        Ok(())
    }
}
fn cycle(name: &str) -> Result<mg::Cycle, OpError> {
    match name {
        | "V" => Ok(mg::Cycle::V),
        | "W" => Ok(mg::Cycle::W),
        | _ => Err(err("cycle must be V or W")),
    }
}
fn controls(
    name: &str,
    smoother: &str,
    pre: usize,
    post: usize,
    tolerance: f64,
    cycles: usize,
) -> Result<mg::MgOptions, OpError> {
    if pre > 8 || post > 8 || cycles == 0 || cycles > 2048 || tolerance <= 0.0 {
        return Err(err("multigrid control work budget exceeded"));
    }
    let smoother = match smoother {
        | "Jacobi" => mg::Smoother::Jacobi,
        | "RedBlackGaussSeidel" => mg::Smoother::RedBlackGaussSeidel,
        | _ => return Err(err("invalid native geometric smoother")),
    };
    Ok(mg::MgOptions {
        cycle: cycle(name)?,
        pre,
        post,
        smoother,
        tol: tolerance,
        max_cycles: cycles,
    })
}
fn setup(
    a: &Csr,
    theta: f64,
    coarse: usize,
    levels: usize,
    sweeps: usize,
    near: &[Vec<f64>],
) -> Result<(), OpError> {
    if a.data.rows == 0
        || a.data.rows != a.data.cols
        || a.data.rows > 512
        || !(0.0..=1.0).contains(&theta)
        || coarse == 0
        || coarse > 512
        || levels == 0
        || levels > 10
        || sweeps > 8
        || near.len() > 16
        || near.iter().any(|v| v.len() != a.data.rows)
    {
        return Err(err("AMG setup shape/options/work budget exceeded"));
    }
    Ok(())
}
object_binding!(mg_amg_build,"num_linear_amg_build","crate::kernels::multigrid::Amg::build","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Csr>,theta:f64,coarse_size:usize,max_levels:usize,sweeps:usize,cycle_name:String,near_nullspace:Vec<Vec<f64>>,smooth:bool)->NewObject<Amg<mg::Amg>>,|ctx,a,theta,coarse_size,max_levels,sweeps,cycle_name,near_nullspace,smooth|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{setup(a,theta,coarse_size,max_levels,sweeps,&near_nullspace)?;let options=mg::AmgOptions{theta,coarse_size,max_levels,sweeps,cycle:cycle(&cycle_name)?,near_nullspace,smooth_prolongator:smooth};let data=mg::Amg::build(&a.data,&options).map_err(|e|err(format!("{e:?}")))?;Ok(NewObject::new(Amg{data,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(mg_amg_general_build,"num_linear_amg_general_build","crate::kernels::amg_general::AmgGeneral::build","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Csr>,kind:String,theta:f64,coarse_size:usize,max_levels:usize,sweeps:usize,cycle_name:String,near_nullspace:Vec<Vec<f64>>)->NewObject<Amg<ag::AmgGeneral>>,|ctx,a,kind,theta,coarse_size,max_levels,sweeps,cycle_name,near_nullspace|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{setup(a,theta,coarse_size,max_levels,sweeps,&near_nullspace)?;let kind=match kind.as_str(){"RugeStuben"=>ag::AmgKind::RugeStuben,"SmoothedAggregation"=>ag::AmgKind::SmoothedAggregation,_=>return Err(err("invalid AMG coarsening kind"))};let options=ag::AmgGeneralOptions{kind,theta,coarse_size,max_levels,sweeps,cycle:cycle(&cycle_name)?,near_nullspace};let data=ag::AmgGeneral::build(&a.data,&options).map_err(|e|err(format!("{e:?}")))?;Ok(NewObject::new(Amg{data,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(mg_amg_num_levels,"num_linear_amg_num_levels","crate::kernels::multigrid::Amg::num_levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.num_levels()))});
object_binding!(mg_amg_operator_complexity,"num_linear_amg_operator_complexity","crate::kernels::multigrid::Amg::operator_complexity","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.operator_complexity()))});
object_binding!(mg_amg_precondition,"num_linear_amg_precondition","crate::kernels::multigrid::Amg::precondition","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.precondition(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_vcycle,"num_linear_amg_vcycle","crate::kernels::multigrid::Amg::vcycle","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>,initial:Vec<f64>)->NewObject<Array>,|ctx,a,rhs,initial|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.len()!=a.n{return Err(err("initial iterate length mismatch"));}let mut out=initial;a.data.vcycle(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_solve,"num_linear_amg_solve","crate::kernels::multigrid::Amg::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;if tolerance<=0.0||max_cycles==0||max_cycles>2048{return Err(err("AMG solve control/work budget invalid"));}a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.as_ref().is_some_and(|v|v.len()!=a.n){return Err(err("initial iterate length mismatch"));}Ok(NewObject::new(MgResult{data:a.data.solve(input,initial.as_deref(),tolerance,max_cycles),scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_amg_general_num_levels,"num_linear_amg_general_num_levels","crate::kernels::amg_general::AmgGeneral::num_levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.num_levels()))});
object_binding!(mg_amg_general_operator_complexity,"num_linear_amg_general_operator_complexity","crate::kernels::amg_general::AmgGeneral::operator_complexity","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.operator_complexity()))});
object_binding!(mg_amg_general_precondition,"num_linear_amg_general_precondition","crate::kernels::amg_general::AmgGeneral::precondition","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.precondition(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_general_vcycle,"num_linear_amg_general_vcycle","crate::kernels::amg_general::AmgGeneral::vcycle","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>,initial:Vec<f64>)->NewObject<Array>,|ctx,a,rhs,initial|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.len()!=a.n{return Err(err("initial iterate length mismatch"));}let mut out=initial;a.data.vcycle(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_general_solve,"num_linear_amg_general_solve","crate::kernels::amg_general::AmgGeneral::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;if tolerance<=0.0||max_cycles==0||max_cycles>2048{return Err(err("AMG solve control/work budget invalid"));}a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.as_ref().is_some_and(|v|v.len()!=a.n){return Err(err("initial iterate length mismatch"));}Ok(NewObject::new(MgResult{data:a.data.solve(input,initial.as_deref(),tolerance,max_cycles),scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_multigrid_mean_factor,"num_linear_multigrid_mean_factor","crate::kernels::multigrid::MgResult::mean_factor","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<MgResult>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.mean_factor()))});
object_binding!(mg_multigrid_result,"num_linear_multigrid_result","","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<MgResult>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::record(vec![("solution",a.data.x.clone().encode()?),("reference",api::quantity(a.scale,d)?),("history",a.data.history.clone().encode()?),("converged",a.data.converged.encode()?)]))))});
object_binding!(mg_geometric_new,"num_linear_geometric_new","crate::kernels::multigrid::GeometricMg::new","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",false,(dimension:usize,n:usize,shift:Normalized<f64>,reference:Input)->Maybe<NewObject<Geometric>>,|ctx,dimension,n,shift,reference|{let _=ctx;if dimension==0||dimension>3||n>127||n.checked_pow(dimension as u32).is_none_or(|n|n>32768){return Err(err("geometric hierarchy extent/work budget exceeded"));}let(s,d)=api::reference(&reference.0)?;Ok(Maybe(mg::GeometricMg::new(dimension,n,shift.0).map(|data|NewObject::new(Geometric{data,n,dimension,scale:s},vec![d]))))});
object_binding!(mg_geometric_levels,"num_linear_geometric_levels","crate::kernels::multigrid::GeometricMg::levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.levels()))});
object_binding!(mg_geometric_size,"num_linear_geometric_size","crate::kernels::multigrid::GeometricMg::size","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.size()))});
object_binding!(mg_geometric_apply,"num_linear_geometric_apply","crate::kernels::multigrid::GeometricMg::apply","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,level:usize,values:Handle<Array>)->NewObject<Array>,|ctx,a,level,values|{let _=ctx;let d=product(dim(ctx,&a.value)?,dim(ctx,&values.value)?)?;a.with(ctx,|a|values.with(ctx,|u|{if level>=a.data.levels(){return Err(err("multigrid level out of range"));}let n=((a.n+1)>>level)-1;let size=n.pow(a.dimension as u32);let input=vector(u,size)?;let mut out=vec![0.0;size];a.data.apply(level,input,&mut out);result(out,scale(a.scale,u.scale)?,d)}))});
object_binding!(mg_geometric_residual_norm,"num_linear_geometric_residual_norm","crate::kernels::multigrid::GeometricMg::residual_norm","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,values:Vec<f64>,rhs:Handle<Array>)->Dynamic,|ctx,a,values,rhs|{let _=ctx;let d=dim(ctx,&rhs.value)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if values.len()!=input.len(){return Err(err("initial iterate length mismatch"));}Ok(Dynamic(api::quantity(a.data.residual_norm(&values,input)*r.scale,d)?))}))});
object_binding!(mg_geometric_solve,"num_linear_geometric_solve","crate::kernels::multigrid::GeometricMg::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let data=a.data.solve(input,initial.as_deref(),&options);Ok(NewObject::new(MgResult{data,scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_geometric_fmg,"num_linear_geometric_fmg","crate::kernels::multigrid::GeometricMg::fmg","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let data=a.data.fmg(input,max_cycles,&options);Ok(NewObject::new(MgResult{data,scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_geometric_cycle_once,"num_linear_geometric_cycle_once","crate::kernels::multigrid::GeometricMg::cycle_once","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<Array>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let mut out=initial.unwrap_or_else(||vec![0.0;input.len()]);a.data.cycle_once(&mut out,input,&options);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &mg_amg_build::BINDING,
    &mg_amg_general_build::BINDING,
    &mg_amg_num_levels::BINDING,
    &mg_amg_operator_complexity::BINDING,
    &mg_amg_precondition::BINDING,
    &mg_amg_vcycle::BINDING,
    &mg_amg_solve::BINDING,
    &mg_amg_general_num_levels::BINDING,
    &mg_amg_general_operator_complexity::BINDING,
    &mg_amg_general_precondition::BINDING,
    &mg_amg_general_vcycle::BINDING,
    &mg_amg_general_solve::BINDING,
    &mg_multigrid_mean_factor::BINDING,
    &mg_multigrid_result::BINDING,
    &mg_geometric_new::BINDING,
    &mg_geometric_levels::BINDING,
    &mg_geometric_size::BINDING,
    &mg_geometric_apply::BINDING,
    &mg_geometric_residual_norm::BINDING,
    &mg_geometric_solve::BINDING,
    &mg_geometric_fmg::BINDING,
    &mg_geometric_cycle_once::BINDING,
];
