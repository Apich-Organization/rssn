//! Native educational Poisson elements with explicit coordinate/solution reference scales.
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::OpError;
use super::Rational;
use super::RuleError;
use super::api;
use super::boundary_models::Markers;
use super::canonical;
use super::err;
use super::fem::Space;
use super::object_binding;
use super::payload;
use super::structured_ops::checked_snapshot;
use super::structured_ops::st;
use crate::sim::physics_fem as native;
fn cap(n: &[usize]) -> Result<Vec<usize>, OpError> {
    if n.is_empty() || n.len() > 3 || n.iter().any(|n| *n == 0 || *n > 512) {
        return Err(err("FEM element counts invalid"));
    }
    let shape = n.iter().map(|n| n + 1).collect::<Vec<_>>();
    if shape
        .iter()
        .try_fold(1usize, |a, b| a.checked_mul(*b))
        .is_none_or(|n| n > 4096)
    {
        return Err(err("FEM allocation/work cap exceeded"));
    }
    Ok(shape)
}
fn source_dim(f: &Frame) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Mul(vec![
        f.dims[2].clone(),
        DimExpr::Pow(Box::new(f.dims[0].clone()), Exp::Rat(Rational::new(-2, 1))),
    ]))
}
object_binding!(fem_legacy_solve1,"num_grid_fem_legacy_solve1","crate::sim::physics_fem::solve_poisson_1d","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize,length:f64)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize,length:f64|{let shape=cap(&elements)?;if elements.len()!=1{return Err(err("FEM dimension mismatch"));}
if length<=0.0{return Err(err("positive domain length required"));}frame.with(ctx,|f|{let scale=f.value/(f.length*f.length);let snap=checked_snapshot(ctx,&source,budget,&vec![f.dims[0].clone();1],&source_dim(f)?)?;let result=native::solve_poisson_1d(elements[0],length,|x|snap.sample(&[x*f.length])/scale);snap.finish()?;payload(f,result.map_err(err)?,shape)})});
object_binding!(fem_legacy_solve2,"num_grid_fem_legacy_solve2","crate::sim::physics_fem::solve_poisson_2d","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize|{let shape=cap(&elements)?;if elements.len()!=2{return Err(err("FEM dimension mismatch"));}frame.with(ctx,|f|{let scale=f.value/(f.length*f.length);let snap=checked_snapshot(ctx,&source,budget,&vec![f.dims[0].clone();2],&source_dim(f)?)?;let result=native::solve_poisson_2d(elements[0],elements[1],|x,y|snap.sample(&[x*f.length,y*f.length])/scale);snap.finish()?;payload(f,result.map_err(err)?,shape)})});
object_binding!(fem_legacy_solve3,"num_grid_fem_legacy_solve3","crate::sim::physics_fem::solve_poisson_3d","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>,elements:Vec<usize>,source:Input,budget:usize|{let shape=cap(&elements)?;if elements.len()!=3{return Err(err("FEM dimension mismatch"));}
if elements.iter().any(|n|*n!=elements[0]){return Err(err("native 3D FEM requires equal subdivisions"));}frame.with(ctx,|f|{let scale=f.value/(f.length*f.length);let snap=checked_snapshot(ctx,&source,budget,&vec![f.dims[0].clone();3],&source_dim(f)?)?;let result=native::solve_poisson_3d(elements[0],|x,y,z|snap.sample(&[x*f.length,y*f.length,z*f.length])/scale);snap.finish()?;payload(f,result.map_err(err)?,shape)})});
object_binding!(fem_legacy_scenario1,"num_grid_fem_legacy_scenario1","crate::sim::physics_fem::simulate_1d_poisson_scenario","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>|{frame.with(ctx,|f|payload(f,native::simulate_1d_poisson_scenario().map_err(err)?,vec![51usize]))});
object_binding!(fem_legacy_scenario2,"num_grid_fem_legacy_scenario2","crate::sim::physics_fem::simulate_2d_poisson_scenario","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>|{frame.with(ctx,|f|payload(f,native::simulate_2d_poisson_scenario().map_err(err)?,vec![21usize,21usize]))});
object_binding!(fem_legacy_scenario3,"num_grid_fem_legacy_scenario3","crate::sim::physics_fem::simulate_3d_poisson_scenario","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(frame:Handle<Frame>)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<Frame>|{frame.with(ctx,|f|payload(f,native::simulate_3d_poisson_scenario().map_err(err)?,vec![11usize,11usize,11usize]))});
object_binding!(fem_legacy_simplex,"num_grid_fem_legacy_simplex","crate::sim::physics_fem::solve_poisson_simplex","Native Poisson finite elements with explicit PDE reference frame. Numeric domain lengths and field arrays are normalized by Length and solution reference U; source callbacks accept measured coordinates (and measured Time for simplex) and return U/Length^2. Allocation/DOF and callback evaluation budgets are checked.",true,(space_handle:Handle<Space>,markers:Handle<Markers>,reaction:f64,source:Input,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,space_handle:Handle<Space>,markers:Handle<Markers>,reaction:f64,source:Input,budget:usize|{space_handle.with(ctx,|s|markers.with(ctx,|m|{s.frame.same(&m.frame)?;if reaction<0.0{return Err(err("nonnegative normalized reaction coefficient required"));}let snap=checked_snapshot(ctx,&source,budget,&[s.frame.dims[1].clone(),s.frame.dims[0].clone(),s.frame.dims[0].clone(),s.frame.dims[0].clone()],&source_dim(&s.frame)?)?;let scale=s.frame.value/(s.frame.length*s.frame.length);let result=native::solve_poisson_simplex(&s.data,reaction,&|t,x|st(&snap,&s.frame,t,x,scale),&m.data);snap.finish()?;for d in m.callbacks.values().flatten().chain(&m.default){d.finish()?;}payload(&s.frame,result.map_err(err)?,vec![s.data.size()])}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &fem_legacy_solve1::BINDING,
            &fem_legacy_solve2::BINDING,
            &fem_legacy_solve3::BINDING,
            &fem_legacy_scenario1::BINDING,
            &fem_legacy_scenario2::BINDING,
            &fem_legacy_scenario3::BINDING,
            &fem_legacy_simplex::BINDING,
        ],
        "numerical_vector",
    )
}
