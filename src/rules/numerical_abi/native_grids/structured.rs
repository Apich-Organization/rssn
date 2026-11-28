//! Checked structured spatial operators and native continuum evolvers.
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Grid;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeArg;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::TypeExpr;
use super::Value;
use super::api;
use super::canonical;
use super::compiled_snapshot;
use super::err;
use super::object_binding;
use super::payload;
use super::structured;
pub(super) fn solver(
    name: &str,
    tolerance: f64,
) -> Result<structured::LinearSolver, OpError> {
    if tolerance <= 0.0 || tolerance > 1.0 {
        return Err(err("linear solver tolerance must be in (0,1]"));
    }
    Ok(match name {
        | "Auto" => structured::LinearSolver::Auto,
        | "Direct" => structured::LinearSolver::Direct,
        | "BiCgStab" => structured::LinearSolver::BiCgStab(tolerance),
        | "Gmres" => structured::LinearSolver::Gmres(tolerance),
        | _ => return Err(err("solver must be Auto, Direct, BiCgStab or Gmres")),
    })
}
fn average(name: &str) -> Result<structured::FaceAverage, OpError> {
    match name {
        | "Harmonic" => Ok(structured::FaceAverage::Harmonic),
        | "Arithmetic" => Ok(structured::FaceAverage::Arithmetic),
        | _ => Err(err("face average must be Harmonic or Arithmetic")),
    }
}
pub(super) fn checked_snapshot(
    ctx: &CallCtx<'_>,
    program: &Input,
    budget: usize,
    inputs: &[DimExpr],
    output: &DimExpr,
) -> Result<CompiledSnapshot, OpError> {
    let s = compiled_snapshot(ctx, &program.0, budget)?;
    if s.input_dimensions.len() != inputs.len() {
        return Err(err(
            "compiled field callback arity differs from its physical domain",
        ));
    }
    for (a, b) in s.input_dimensions.iter().zip(inputs) {
        api::same_dimension(a, b)?;
    }
    api::same_dimension(&s.output_dimension, output)?;
    Ok(s)
}
fn spatial(
    ctx: &CallCtx<'_>,
    g: &Grid,
    p: &Input,
    budget: usize,
    output: &DimExpr,
) -> Result<CompiledSnapshot, OpError> {
    checked_snapshot(ctx, p, budget, &vec![g.frame.dims[0].clone(); 3], output)
}
pub(super) fn space(
    s: &CompiledSnapshot,
    g: &Frame,
    x: [f64; 3],
    out_scale: f64,
) -> f64 {
    s.sample(&x.map(|x| x * g.length)) / out_scale
}
pub(super) fn st(
    s: &CompiledSnapshot,
    g: &Frame,
    t: f64,
    x: [f64; 3],
    out_scale: f64,
) -> f64 {
    s.sample(&[
        t * g.time,
        x[0] * g.length,
        x[1] * g.length,
        x[2] * g.length,
    ]) / out_scale
}
fn field(
    g: &Grid,
    v: &[f64],
) -> Result<(), OpError> {
    if v.len() != g.data.len() {
        return Err(err("PDE field length does not match grid"));
    }
    Ok(())
}
fn evolution(
    g: &Grid,
    steps: usize,
    t_end: f64,
) -> Result<(), OpError> {
    if steps == 0
        || steps > 8192
        || t_end <= 0.0
        || !t_end.is_finite()
        || g.data
            .len()
            .checked_mul(steps)
            .is_none_or(|n| n > 8_388_608)
        || g.data.len() > 4096
    {
        return Err(err("PDE step/work/allocation budget exceeded"));
    }
    Ok(())
}
pub(super) fn dims_forcing(
    f: &Frame,
    time_power: i64,
) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Mul(vec![
        f.dims[2].clone(),
        DimExpr::Pow(
            Box::new(f.dims[1].clone()),
            Exp::Rat(Rational::new(-time_power, 1)),
        ),
    ]))
}
fn diffusivity(f: &Frame) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Mul(vec![
        DimExpr::Pow(Box::new(f.dims[0].clone()), Exp::Rat(Rational::new(2, 1))),
        DimExpr::Inverse(Box::new(f.dims[1].clone())),
    ]))
}
fn adaptive(
    tolerance: f64,
    dt0: f64,
    dt_min: f64,
    dt_max: f64,
    max_steps: usize,
) -> Result<structured::AdaptiveOptions, OpError> {
    if tolerance <= 0.0
        || dt_min <= 0.0
        || dt0 < dt_min
        || dt_max < dt0
        || max_steps == 0
        || max_steps > 8192
    {
        return Err(err("invalid adaptive tolerance/steps or work budget"));
    }
    Ok(structured::AdaptiveOptions {
        tol: tolerance,
        dt0,
        dt_min,
        dt_max,
        max_steps,
    })
}
fn adaptive_payload(
    g: &Grid,
    out: Vec<f64>,
    stats: structured::AdaptiveStats,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("field", payload(&g.frame, out, vec![g.data.len()])?.0),
        ("accepted", stats.accepted.encode()?),
        ("rejected", stats.rejected.encode()?),
        (
            "dt_min",
            api::quantity(stats.dt_min_used * g.frame.time, g.frame.dims[1].clone())?,
        ),
        (
            "dt_max",
            api::quantity(stats.dt_max_used * g.frame.time, g.frame.dims[1].clone())?,
        ),
    ])))
}
struct Laplacian {
    data: structured::LaplacianOp,
    grid: Grid,
}
struct Diffusion {
    data: structured::VariableDiffusionOp,
    grid: Grid,
}
impl NativeType for Laplacian {
    const KIND: &'static str = "rssn/native/laplacian/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.grid.validate()?;
        if self.grid.data.len() > 4096 {
            return Err(err("native operator assembly budget exceeded"));
        }
        Ok(())
    }
}
impl NativeType for Diffusion {
    const KIND: &'static str = "rssn/native/variable-diffusion/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.grid.validate()?;
        if self.grid.data.len() > 4096 {
            return Err(err("native operator assembly budget exceeded"));
        }
        Ok(())
    }
}
object_binding!(structured_structured_sample,"num_grid_structured_sample","crate::sim::structured::Grid::sample","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,program:Input,budget:usize)->Dynamic,|ctx,grid,program,budget|{let _=ctx;grid.with(ctx,|g|{let s=spatial(ctx,g,&program,budget,&g.frame.dims[2])?;let exact=|x|space(&s,&g.frame,x,g.frame.value);let out=g.data.sample(exact);s.finish()?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_structured_max_error,"num_grid_structured_max_error","crate::sim::structured::Grid::max_error","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,program:Input,budget:usize,values:Vec<f64>)->Dynamic,|ctx,grid,program,budget,values|{let _=ctx;grid.with(ctx,|g|{let s=spatial(ctx,g,&program,budget,&g.frame.dims[2])?;let exact=|x|space(&s,&g.frame,x,g.frame.value);field(g,&values)?;let out=g.data.max_error(&values,exact);s.finish()?;Ok(Dynamic(api::quantity(out*g.frame.value,g.frame.dims[2].clone())?))})});
object_binding!(structured_structured_l2_error,"num_grid_structured_l2_error","crate::sim::structured::Grid::l2_error","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,program:Input,budget:usize,values:Vec<f64>)->Dynamic,|ctx,grid,program,budget,values|{let _=ctx;grid.with(ctx,|g|{let s=spatial(ctx,g,&program,budget,&g.frame.dims[2])?;let exact=|x|space(&s,&g.frame,x,g.frame.value);field(g,&values)?;let out=g.data.l2_error(&values,exact);s.finish()?;let d=canonical(DimExpr::Mul(vec![g.frame.dims[2].clone(),DimExpr::Pow(Box::new(g.frame.dims[0].clone()),Exp::Rat(Rational::new(g.data.dim as i64,2)))]))?;Ok(Dynamic(api::quantity(out*g.frame.value*g.frame.length.powf(g.data.dim as f64/2.0),d)?))})});
object_binding!(structured_laplacian_new,"num_grid_laplacian_new","crate::sim::structured::LaplacianOp::new","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>)->NewObject<Laplacian>,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{if g.data.len()>4096{return Err(err("operator assembly budget exceeded"));}let data=structured::LaplacianOp::new(&g.data,&g.bcs.data).map_err(err)?;Ok(NewObject::dependent(Laplacian{data,grid:g.clone()},g.frame.dims.to_vec(),vec![grid.value.clone()]))})});
object_binding!(structured_variable_diffusion_new,"num_grid_variable_diffusion_new","crate::sim::structured::VariableDiffusionOp::new","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,alpha:Vec<f64>,averaging:String)->NewObject<Diffusion>,|ctx,grid,alpha,averaging|{let _=ctx;grid.with(ctx,|g|{field(g,&alpha)?;if g.data.len()>4096{return Err(err("operator assembly budget exceeded"));}let data=structured::VariableDiffusionOp::new(&g.data,&g.bcs.data,&alpha,average(&averaging)?).map_err(err)?;Ok(NewObject::dependent(Diffusion{data,grid:g.clone()},g.frame.dims.to_vec(),vec![grid.value.clone()]))})});
object_binding!(structured_laplacian_dirichlet_count,"num_grid_laplacian_dirichlet_count","crate::sim::structured::LaplacianOp::dirichlet_count","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>)->usize,|ctx,op|{let _=ctx;op.with(ctx,|a|{Ok(a.data.dirichlet_count())})});
object_binding!(structured_laplacian_has_absorbing,"num_grid_laplacian_has_absorbing","crate::sim::structured::LaplacianOp::has_absorbing","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>)->bool,|ctx,op|{let _=ctx;op.with(ctx,|a|{Ok(a.data.has_absorbing())})});
object_binding!(structured_laplacian_is_dirichlet,"num_grid_laplacian_is_dirichlet","crate::sim::structured::LaplacianOp::is_dirichlet","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>,index:usize)->bool,|ctx,op,index|{let _=ctx;op.with(ctx,|a|{if index>=a.grid.data.len(){return Err(err("node index out of range"));}Ok(a.data.is_dirichlet(index))})});
object_binding!(structured_laplacian_apply,"num_grid_laplacian_apply","crate::sim::structured::LaplacianOp::apply","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>,values:Vec<f64>)->Dynamic,|ctx,op,values|{let _=ctx;op.with(ctx,|a|{field(&a.grid,&values)?;let mut out=vec![0.0;values.len()];a.data.apply(&values,&mut out);operator_payload(&a.grid,out,false)})});
object_binding!(structured_laplacian_forcing,"num_grid_laplacian_forcing","crate::sim::structured::LaplacianOp::forcing","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>,time:f64)->Dynamic,|ctx,op,time|{let _=ctx;op.with(ctx,|a|{let out=a.data.forcing(time);a.grid.bcs.finish_for(a.grid.data.dim)?;operator_payload(&a.grid,out,false)})});
object_binding!(structured_laplacian_set_dirichlet,"num_grid_laplacian_set_dirichlet","crate::sim::structured::LaplacianOp::set_dirichlet","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>,values:Vec<f64>,time:f64)->Dynamic,|ctx,op,values,time|{let _=ctx;op.with(ctx,|a|{field(&a.grid,&values)?;let mut out=values;a.data.set_dirichlet(&mut out,time);a.grid.bcs.finish_for(a.grid.data.dim)?;payload(&a.grid.frame,out,vec![a.grid.data.len()])})});
object_binding!(structured_laplacian_system,"num_grid_laplacian_system","crate::sim::structured::LaplacianOp::system","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Laplacian>,shift:f64,scale:f64)->Dynamic,|ctx,op,shift,scale|{let _=ctx;op.with(ctx,|a|{let m=a.data.system(shift,scale);Ok(Dynamic(api::record(vec![("rows",m.rows.encode()?),("cols",m.cols.encode()?),("indptr",m.indptr.encode()?),("indices",m.indices.encode()?),("values",m.values.encode()?),("frame",payload(&a.grid.frame,vec![],vec![0])?.0)])))})});
object_binding!(structured_variable_diffusion_dirichlet_count,"num_grid_variable_diffusion_dirichlet_count","crate::sim::structured::VariableDiffusionOp::dirichlet_count","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>)->usize,|ctx,op|{let _=ctx;op.with(ctx,|a|{Ok(a.data.dirichlet_count())})});
object_binding!(structured_variable_diffusion_has_absorbing,"num_grid_variable_diffusion_has_absorbing","crate::sim::structured::VariableDiffusionOp::has_absorbing","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>)->bool,|ctx,op|{let _=ctx;op.with(ctx,|a|{Ok(a.data.has_absorbing())})});
object_binding!(structured_variable_diffusion_is_dirichlet,"num_grid_variable_diffusion_is_dirichlet","crate::sim::structured::VariableDiffusionOp::is_dirichlet","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>,index:usize)->bool,|ctx,op,index|{let _=ctx;op.with(ctx,|a|{if index>=a.grid.data.len(){return Err(err("node index out of range"));}Ok(a.data.is_dirichlet(index))})});
object_binding!(structured_variable_diffusion_apply,"num_grid_variable_diffusion_apply","crate::sim::structured::VariableDiffusionOp::apply","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>,values:Vec<f64>)->Dynamic,|ctx,op,values|{let _=ctx;op.with(ctx,|a|{field(&a.grid,&values)?;let mut out=vec![0.0;values.len()];a.data.apply(&values,&mut out);operator_payload(&a.grid,out,true)})});
object_binding!(structured_variable_diffusion_forcing,"num_grid_variable_diffusion_forcing","crate::sim::structured::VariableDiffusionOp::forcing","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>,time:f64)->Dynamic,|ctx,op,time|{let _=ctx;op.with(ctx,|a|{let out=a.data.forcing(time);a.grid.bcs.finish_for(a.grid.data.dim)?;operator_payload(&a.grid,out,true)})});
object_binding!(structured_variable_diffusion_set_dirichlet,"num_grid_variable_diffusion_set_dirichlet","crate::sim::structured::VariableDiffusionOp::set_dirichlet","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>,values:Vec<f64>,time:f64)->Dynamic,|ctx,op,values,time|{let _=ctx;op.with(ctx,|a|{field(&a.grid,&values)?;let mut out=values;a.data.set_dirichlet(&mut out,time);a.grid.bcs.finish_for(a.grid.data.dim)?;payload(&a.grid.frame,out,vec![a.grid.data.len()])})});
object_binding!(structured_variable_diffusion_system,"num_grid_variable_diffusion_system","crate::sim::structured::VariableDiffusionOp::system","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(op:Handle<Diffusion>,shift:f64,scale:f64)->Dynamic,|ctx,op,shift,scale|{let _=ctx;op.with(ctx,|a|{let m=a.data.system(shift,scale);Ok(Dynamic(api::record(vec![("rows",m.rows.encode()?),("cols",m.cols.encode()?),("indptr",m.indptr.encode()?),("indices",m.indices.encode()?),("values",m.values.encode()?),("frame",payload(&a.grid.frame,vec![],vec![0])?.0)])))})});
object_binding!(structured_solve_helmholtz,"num_grid_solve_helmholtz","crate::sim::structured::solve_helmholtz","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,sigma:f64,program:Input,budget:usize,solver_name:String,tolerance:f64)->Dynamic,|ctx,grid,sigma,program,budget,solver_name,tolerance|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_m2_0_1=g.frame.reference([-2, 0, 1])?;if g.data.len()>4096{return Err(err("elliptic solve allocation/work budget exceeded"));}let output=canonical(DimExpr::Mul(vec![g.frame.dims[2].clone(),DimExpr::Pow(Box::new(g.frame.dims[0].clone()),Exp::Rat(Rational::new(-2,1)))]))?;let s=spatial(ctx,g,&program,budget,&output)?;let source=|x|space(&s,&g.frame,x,g_frame_reference_m2_0_1);let out=structured::solve_helmholtz(&g.data,&g.bcs.data,sigma,source,solver(&solver_name,tolerance)?).map_err(err)?;s.finish()?;g.bcs.finish_for(g.data.dim)?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_solve_variable_helmholtz,"num_grid_solve_variable_helmholtz","crate::sim::structured::solve_variable_helmholtz","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,sigma:f64,program:Input,budget:usize,solver_name:String,tolerance:f64,alpha:Vec<f64>,averaging:String)->Dynamic,|ctx,grid,sigma,program,budget,solver_name,tolerance,alpha,averaging|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_m2_0_1=g.frame.reference([-2, 0, 1])?;if g.data.len()>4096{return Err(err("elliptic solve allocation/work budget exceeded"));}let output=canonical(DimExpr::Mul(vec![g.frame.dims[2].clone(),DimExpr::Pow(Box::new(g.frame.dims[0].clone()),Exp::Rat(Rational::new(-2,1)))]))?;let s=spatial(ctx,g,&program,budget,&output)?;let source=|x|space(&s,&g.frame,x,g_frame_reference_m2_0_1);field(g,&alpha)?;let out=structured::solve_variable_helmholtz(&g.data,&g.bcs.data,&alpha,average(&averaging)?,sigma,source,solver(&solver_name,tolerance)?).map_err(err)?;s.finish()?;g.bcs.finish_for(g.data.dim)?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_heat_fixed,"num_grid_heat_fixed","crate::sim::structured::heat_fixed","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,initial:Input,source:Input,budget:usize,alpha:f64,t_end:f64,steps:usize,theta:f64)->Dynamic,|ctx,grid,initial,source,budget,alpha,t_end,steps,theta|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -1, 1])?};evolution(g,steps,t_end)?;if !(0.0..=1.0).contains(&theta){return Err(err("theta must be in [0,1]"));}let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let source=if matches!(&source.0,Value::Variant(v) if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone()],&dims_forcing(&g.frame,1)?)?)};let initial_fn=|x|space(&init,&g.frame,x,g.frame.value);let source_fn=|t,x|source.as_ref().map_or(0.0,|s|st(s,&g.frame,t,x,g_frame_reference_0_m1_1));if alpha<0.0{return Err(err("diffusivity must be nonnegative"));}let problem=structured::HeatProblem{grid:&g.data,bcs:&g.bcs.data,alpha,init:&initial_fn,source:Some(&source_fn),solver:structured::LinearSolver::Auto};let out=structured::heat_fixed(&problem,t_end,steps,theta).map_err(err)?;init.finish()?;if let Some(source)=source{source.finish()?;}g.bcs.finish_for(g.data.dim)?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_heat_adaptive,"num_grid_heat_adaptive","crate::sim::structured::heat_adaptive","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,initial:Input,source:Input,budget:usize,alpha:f64,t_end:f64,tolerance:f64,dt0:f64,dt_min:f64,dt_max:f64,max_steps:usize)->Dynamic,|ctx,grid,initial,source,budget,alpha,t_end,tolerance,dt0,dt_min,dt_max,max_steps|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -1, 1])?};let opts=adaptive(tolerance,dt0,dt_min,dt_max,max_steps)?;evolution(g,max_steps,t_end)?;let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let source=if matches!(&source.0,Value::Variant(v) if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone()],&dims_forcing(&g.frame,1)?)?)};let initial_fn=|x|space(&init,&g.frame,x,g.frame.value);let source_fn=|t,x|source.as_ref().map_or(0.0,|s|st(s,&g.frame,t,x,g_frame_reference_0_m1_1));if alpha<0.0{return Err(err("diffusivity must be nonnegative"));}let problem=structured::HeatProblem{grid:&g.data,bcs:&g.bcs.data,alpha,init:&initial_fn,source:Some(&source_fn),solver:structured::LinearSolver::Auto};let out=structured::heat_adaptive(&problem,t_end,&opts).map_err(err)?;init.finish()?;if let Some(source)=source{source.finish()?;}g.bcs.finish_for(g.data.dim)?;adaptive_payload(g,out.0,out.1)})});
object_binding!(structured_heat_variable_fixed,"num_grid_heat_variable_fixed","crate::sim::structured::heat_variable_fixed","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,initial:Input,source:Input,budget:usize,alpha:Vec<f64>,averaging:String,t_end:f64,steps:usize,theta:f64)->Dynamic,|ctx,grid,initial,source,budget,alpha,averaging,t_end,steps,theta|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -1, 1])?};evolution(g,steps,t_end)?;if !(0.0..=1.0).contains(&theta){return Err(err("theta must be in [0,1]"));}let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let source=if matches!(&source.0,Value::Variant(v) if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone()],&dims_forcing(&g.frame,1)?)?)};let initial_fn=|x|space(&init,&g.frame,x,g.frame.value);let source_fn=|t,x|source.as_ref().map_or(0.0,|s|st(s,&g.frame,t,x,g_frame_reference_0_m1_1));field(g,&alpha)?;let problem=structured::VariableHeatProblem{grid:&g.data,bcs:&g.bcs.data,alpha:&alpha,averaging:average(&averaging)?,init:&initial_fn,source:Some(&source_fn),solver:structured::LinearSolver::Auto};let out=structured::heat_variable_fixed(&problem,t_end,steps,theta).map_err(err)?;init.finish()?;if let Some(source)=source{source.finish()?;}g.bcs.finish_for(g.data.dim)?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_heat_variable_adaptive,"num_grid_heat_variable_adaptive","crate::sim::structured::heat_variable_adaptive","Native structured PDE method in the grid reference frame; normalized numeric inputs are converted to coherent callback coordinates and checked physical callback output domains. Callback and discretization work are bounded; native boundary errors propagate.",true,(grid:Handle<Grid>,initial:Input,source:Input,budget:usize,alpha:Vec<f64>,averaging:String,t_end:f64,tolerance:f64,dt0:f64,dt_min:f64,dt_max:f64,max_steps:usize)->Dynamic,|ctx,grid,initial,source,budget,alpha,averaging,t_end,tolerance,dt0,dt_min,dt_max,max_steps|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=g.frame.reference([0, -1, 1])?;let opts=adaptive(tolerance,dt0,dt_min,dt_max,max_steps)?;evolution(g,max_steps,t_end)?;let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let source=if matches!(&source.0,Value::Variant(v) if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone()],&dims_forcing(&g.frame,1)?)?)};let initial_fn=|x|space(&init,&g.frame,x,g.frame.value);let source_fn=|t,x|source.as_ref().map_or(0.0,|s|st(s,&g.frame,t,x,g_frame_reference_0_m1_1));field(g,&alpha)?;let problem=structured::VariableHeatProblem{grid:&g.data,bcs:&g.bcs.data,alpha:&alpha,averaging:average(&averaging)?,init:&initial_fn,source:Some(&source_fn),solver:structured::LinearSolver::Auto};let out=structured::heat_variable_adaptive(&problem,t_end,&opts).map_err(err)?;init.finish()?;if let Some(source)=source{source.finish()?;}g.bcs.finish_for(g.data.dim)?;adaptive_payload(g,out.0,out.1)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &structured_wave_leapfrog::BINDING,
    &structured_nonlinear_diffusion::BINDING,
    &structured_nonlinear_diffusion_adaptive::BINDING,
    &structured_linear_system_new::BINDING,
    &structured_linear_system_solve::BINDING,
    &structured_structured_sample::BINDING,
    &structured_structured_max_error::BINDING,
    &structured_structured_l2_error::BINDING,
    &structured_laplacian_new::BINDING,
    &structured_variable_diffusion_new::BINDING,
    &structured_laplacian_dirichlet_count::BINDING,
    &structured_laplacian_has_absorbing::BINDING,
    &structured_laplacian_is_dirichlet::BINDING,
    &structured_laplacian_apply::BINDING,
    &structured_laplacian_forcing::BINDING,
    &structured_laplacian_set_dirichlet::BINDING,
    &structured_laplacian_system::BINDING,
    &structured_variable_diffusion_dirichlet_count::BINDING,
    &structured_variable_diffusion_has_absorbing::BINDING,
    &structured_variable_diffusion_is_dirichlet::BINDING,
    &structured_variable_diffusion_apply::BINDING,
    &structured_variable_diffusion_forcing::BINDING,
    &structured_variable_diffusion_set_dirichlet::BINDING,
    &structured_variable_diffusion_system::BINDING,
    &structured_solve_helmholtz::BINDING,
    &structured_solve_variable_helmholtz::BINDING,
    &structured_heat_fixed::BINDING,
    &structured_heat_adaptive::BINDING,
    &structured_heat_variable_fixed::BINDING,
    &structured_heat_variable_adaptive::BINDING,
];

fn operator_payload(
    grid: &Grid,
    out: Vec<f64>,
    variable: bool,
) -> Result<Dynamic, OpError> {
    let mut f = grid.frame.clone();
    if variable {
        f.value /= f.time;
        f.dims[2] = dims_forcing(&grid.frame, 1)?;
    } else {
        f.value /= f.length * f.length;
        f.dims[2] = canonical(DimExpr::Mul(vec![
            f.dims[2].clone(),
            DimExpr::Pow(Box::new(f.dims[0].clone()), Exp::Rat(Rational::new(-2, 1))),
        ]))?;
    }
    payload(&f, out, vec![grid.data.len()])
}
object_binding!(structured_wave_leapfrog,"num_grid_wave_leapfrog","crate::sim::structured::wave_leapfrog","Native PDE in an explicit reference frame. Callback coordinates and outputs retain their physical domain; normalized numeric coefficients are checked. Native step/factor budgets, boundary errors and callback errors propagate.",true,(grid:Handle<Grid>,initial:Input,velocity:Input,source:Input,budget:usize,c:f64,t_end:f64,steps:usize)->Dynamic,|ctx,grid,initial,velocity,source,budget,c,t_end,steps|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m2_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -2, 1])?};let g_frame_reference_0_m1_1=g.frame.reference([0, -1, 1])?;evolution(g,steps,t_end)?;if c<=0.0{return Err(err("wave speed must be positive"));}let u=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let v=spatial(ctx,g,&velocity,budget,&dims_forcing(&g.frame,1)?)?;let source=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone()],&dims_forcing(&g.frame,2)?)?)};let ufn=|x|space(&u,&g.frame,x,g.frame.value);let vfn=|x|space(&v,&g.frame,x,g_frame_reference_0_m1_1);let sfn=|t,x|source.as_ref().map_or(0.0,|s|st(s,&g.frame,t,x,g_frame_reference_0_m2_1));let problem=structured::WaveProblem{grid:&g.data,bcs:&g.bcs.data,c,u0:&ufn,v0:&vfn,source:Some(&sfn)};let out=structured::wave_leapfrog(&problem,t_end,steps).map_err(err)?;u.finish()?;v.finish()?;if let Some(s)=source{s.finish()?;}g.bcs.finish_for(g.data.dim)?;payload(&g.frame,out,vec![g.data.len()])})});
object_binding!(structured_nonlinear_diffusion,"num_grid_nonlinear_diffusion","crate::sim::structured::nonlinear_diffusion","Native PDE in an explicit reference frame. Callback coordinates and outputs retain their physical domain; normalized numeric coefficients are checked. Native step/factor budgets, boundary errors and callback errors propagate.",true,(grid:Handle<Grid>,initial:Input,diffusivity_program:Input,source:Input,budget:usize,t_end:f64,steps:usize,theta:f64)->Dynamic,|ctx,grid,initial,diffusivity_program,source,budget,t_end,steps,theta|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -1, 1])?};let g_frame_reference_2_m1_0=g.frame.reference([2, -1, 0])?;evolution(g,steps,t_end)?;if !(0.0..=1.0).contains(&theta){return Err(err("theta must be in [0,1]"));}let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let diffusion=checked_snapshot(ctx,&diffusivity_program,budget,&[g.frame.dims[2].clone()],&diffusivity(&g.frame)?)?;let source=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[2].clone()],&dims_forcing(&g.frame,1)?)?)};let ufn=|x|space(&init,&g.frame,x,g.frame.value);let dfn=|u|diffusion.sample(&[u*g.frame.value])/(g_frame_reference_2_m1_0);let sfn=|t,x:[f64;3],u|source.as_ref().map_or(0.0,|s|s.sample(&[t*g.frame.time,x[0]*g.frame.length,x[1]*g.frame.length,x[2]*g.frame.length,u*g.frame.value])/(g_frame_reference_0_m1_1));let problem=structured::NonlinearDiffusion{grid:&g.data,bcs:&g.bcs.data,d:&dfn,source:Some(&sfn),init:&ufn};let out=structured::nonlinear_diffusion(&problem,t_end,steps,theta).map_err(err)?;init.finish()?;diffusion.finish()?;if let Some(s)=source{s.finish()?;}g.bcs.finish_for(g.data.dim)?;Ok(Dynamic(api::record(vec![("field",payload(&g.frame,out.0,vec![g.data.len()])?.0),("newton_iters",out.1.newton_iters.encode()?),("krylov_iters",out.1.krylov_iters.encode()?)])))})});
object_binding!(structured_nonlinear_diffusion_adaptive,"num_grid_nonlinear_diffusion_adaptive","crate::sim::structured::nonlinear_diffusion_adaptive","Native PDE in an explicit reference frame. Callback coordinates and outputs retain their physical domain; normalized numeric coefficients are checked. Native step/factor budgets, boundary errors and callback errors propagate.",true,(grid:Handle<Grid>,initial:Input,diffusivity_program:Input,source:Input,budget:usize,t_end:f64,tolerance:f64,dt0:f64,dt_min:f64,dt_max:f64,max_steps:usize)->Dynamic,|ctx,grid,initial,diffusivity_program,source,budget,t_end,tolerance,dt0,dt_min,dt_max,max_steps|{let _=ctx;grid.with(ctx,|g|{let g_frame_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{g.frame.reference([0, -1, 1])?};let g_frame_reference_2_m1_0=g.frame.reference([2, -1, 0])?;let opts=adaptive(tolerance,dt0,dt_min,dt_max,max_steps)?;evolution(g,max_steps,t_end)?;let init=spatial(ctx,g,&initial,budget,&g.frame.dims[2])?;let diffusion=checked_snapshot(ctx,&diffusivity_program,budget,&[g.frame.dims[2].clone()],&diffusivity(&g.frame)?)?;let source=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[g.frame.dims[1].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[0].clone(),g.frame.dims[2].clone()],&dims_forcing(&g.frame,1)?)?)};let ufn=|x|space(&init,&g.frame,x,g.frame.value);let dfn=|u|diffusion.sample(&[u*g.frame.value])/(g_frame_reference_2_m1_0);let sfn=|t,x:[f64;3],u|source.as_ref().map_or(0.0,|s|s.sample(&[t*g.frame.time,x[0]*g.frame.length,x[1]*g.frame.length,x[2]*g.frame.length,u*g.frame.value])/(g_frame_reference_0_m1_1));let problem=structured::NonlinearDiffusion{grid:&g.data,bcs:&g.bcs.data,d:&dfn,source:Some(&sfn),init:&ufn};let out=structured::nonlinear_diffusion_adaptive(&problem,t_end,&opts).map_err(err)?;init.finish()?;diffusion.finish()?;if let Some(s)=source{s.finish()?;}g.bcs.finish_for(g.data.dim)?;Ok(Dynamic(api::record(vec![("solution",adaptive_payload(g,out.0,out.1)?.0),("newton_iters",out.2.newton_iters.encode()?),("krylov_iters",out.2.krylov_iters.encode()?)])))})});

struct System {
    data: structured::LinearSystem,
    frame: Frame,
    n: usize,
}
impl NativeType for System {
    const KIND: &'static str = "rssn/native/structured-linear-system/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.n > 512 {
            return Err(err("factored system side budget exceeded"));
        }
        Ok(())
    }
}
struct CsrInput(crate::kernels::krylov::Csr);
impl NativeArg for CsrInput {
    fn ty() -> TypeExpr {
        super::super::iterative::csr_type()
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(v)?;
        let m = super::super::iterative::decode_csr(v)?;
        if m.rows > 512 {
            return Err(err("factored system side budget exceeded"));
        }
        Ok(Self(m))
    }
}
object_binding!(structured_linear_system_new,"num_grid_linear_system_new","crate::sim::structured::LinearSystem::new","Native PDE in an explicit reference frame. Callback coordinates and outputs retain their physical domain; normalized numeric coefficients are checked. Native step/factor budgets, boundary errors and callback errors propagate.",true,(frame:Handle<Frame>,matrix:CsrInput,solver_name:String,tolerance:f64)->NewObject<System>,|ctx,frame,matrix,solver_name,tolerance|{let _=ctx;frame.with(ctx,|f|{let n=matrix.0.rows;Ok(NewObject::new(System{data:structured::LinearSystem::new(matrix.0,solver(&solver_name,tolerance)?).map_err(err)?,frame:f.clone(),n},f.dims.to_vec()))})});
object_binding!(structured_linear_system_solve,"num_grid_linear_system_solve","crate::sim::structured::LinearSystem::solve","Native PDE in an explicit reference frame. Callback coordinates and outputs retain their physical domain; normalized numeric coefficients are checked. Native step/factor budgets, boundary errors and callback errors propagate.",true,(system:Handle<System>,rhs:Vec<f64>,initial:Option<Vec<f64>>)->Dynamic,|ctx,system,rhs,initial|{let _=ctx;system.with(ctx,|s|{if rhs.len()!=s.n||initial.as_ref().is_some_and(|v|v.len()!=s.n){return Err(err("system RHS/initial length mismatch"));}payload(&s.frame,s.data.solve(&rhs,initial.as_deref()).map_err(err)?,vec![s.n])})});
