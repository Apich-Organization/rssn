//! Native continuum multigrid with physical source and boundary references.
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Faces;
use super::Frame;
use super::Grid;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
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
use super::err;
use super::grid_size;
use super::object_binding;
use super::payload;
use super::structured_ops::checked_snapshot;
use super::structured_ops::space;
use crate::sim::multigrid as mg;
struct Options(mg::MgOptions);
impl NativeArg for Options {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("cycle".into(), TypeExpr::String),
            ("pre".into(), TypeExpr::Int),
            ("post".into(), TypeExpr::Int),
            ("max_cycles".into(), TypeExpr::Int),
            ("tol".into(), TypeExpr::Float),
            ("coarsest_nodes".into(), TypeExpr::Int),
        ])
    }

    fn decode(
        v: &Value,
        i: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(v)?;
        let Value::Record(r) = v else {
            return Err(err("multigrid options must be a record"));
        };
        if r.fields.len() != 6 {
            return Err(err("multigrid options require exactly six fields"));
        }
        let u =
            |n| <usize as crate::rules::numerical_abi::codecs::Arg>::decode(api::field(v, n)?, i);
        let name = <String as crate::rules::numerical_abi::codecs::Arg>::decode(
            api::field(v, "cycle")?,
            i,
        )?;
        let cycle = match name.as_str() {
            | "V" => mg::MgCycle::V,
            | "W" => mg::MgCycle::W,
            | "Fmg" => mg::MgCycle::Fmg,
            | _ => return Err(err("cycle must be V, W or Fmg")),
        };
        let options = mg::MgOptions {
            cycle,
            pre: u("pre")?,
            post: u("post")?,
            max_cycles: u("max_cycles")?,
            tol: <f64 as crate::rules::numerical_abi::codecs::Arg>::decode(
                api::field(v, "tol")?,
                i,
            )?,
            coarsest_nodes: u("coarsest_nodes")?,
        };
        if options.pre > 8
            || options.post > 8
            || options.max_cycles == 0
            || options.max_cycles > 512
            || options.tol <= 0.0
            || options.coarsest_nodes == 0
            || options.coarsest_nodes > 4096
        {
            return Err(err("multigrid options exceed work bounds"));
        }
        Ok(Self(options))
    }
}
struct Hierarchy {
    data: mg::MultigridSolver,
    grid: Grid,
}
impl NativeType for Hierarchy {
    const KIND: &'static str = "rssn/native/continuum-multigrid/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.grid.validate()?;
        if self.grid.data.len() > 4096 {
            return Err(err("continuum multigrid node budget exceeded"));
        }
        Ok(())
    }
}
struct Solution {
    data: mg::MgResult,
    frame: Frame,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/continuum-multigrid-result/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 16384;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.u.len() > LIMIT
            || self.data.residuals.len() > 513
            || self
                .data
                .u
                .iter()
                .chain(&self.data.residuals)
                .any(|v| !v.is_finite())
        {
            return Err(err("invalid finite multigrid result or work budget"));
        }
        Ok(())
    }
}
object_binding!(new,"num_grid_multigrid_new","crate::sim::multigrid::MultigridSolver::new","Actual native structured multigrid hierarchy with explicit Length/Time/field references and original V/W/Fmg controls. Nonzero/mixed/periodic native boundary conditions are retained; sources are normalized by field/length².",true,(grid:Handle<Grid>,sigma:f64,options:Options)->NewObject<Hierarchy>,|ctx,grid,sigma,options|{grid.with(ctx,|g|{if g.data.len()>4096||g.data.len().checked_mul(options.0.max_cycles).is_none_or(|w|w>8_388_608){return Err(err("multigrid hierarchy work budget exceeded"));}let data=mg::MultigridSolver::new(&g.data,&g.bcs.data,sigma,options.0).map_err(err)?;Ok(NewObject::dependent(Hierarchy{data,grid:g.clone()},g.frame.dims.to_vec(),vec![grid.value.clone()]))})});
object_binding!(levels,"num_grid_multigrid_num_levels","crate::sim::multigrid::MultigridSolver::num_levels","Actual number of native continuum multigrid levels.",true,(hierarchy:Handle<Hierarchy>)->usize,|ctx,hierarchy|{hierarchy.with(ctx,|h|Ok(h.data.num_levels()))});
object_binding!(solve,"num_grid_multigrid_solve","crate::sim::multigrid::MultigridSolver::solve","Native continuum multigrid solve, preserving normalized residual history and convergence status. Source arrays use field/length² reference and time is normalized by the retained Time reference; output field retains its physical amplitude.",true,(hierarchy:Handle<Hierarchy>,source:Vec<f64>,time:f64)->NewObject<Solution>,|ctx,hierarchy,source,time|{hierarchy.with(ctx,|h|{if source.len()!=h.grid.data.len(){return Err(err("multigrid source length mismatch"));}let data=h.data.solve(&source,time).map_err(err)?;h.grid.bcs.finish_for(h.grid.data.dim)?;Ok(NewObject::new(Solution{data,frame:h.grid.frame.clone()},h.grid.frame.dims.to_vec()))})});
object_binding!(factor,"num_grid_multigrid_convergence_factor","crate::sim::multigrid::MgResult::convergence_factor","Native residual reduction factor, without replacing the kernel's stopping report.",true,(solution:Handle<Solution>)->f64,|ctx,solution|{solution.with(ctx,|s|Ok(s.data.convergence_factor()))});
object_binding!(data,"num_grid_multigrid_result","","Native continuum multigrid solution with explicit physical field reference and original normalized residual history/status.",true,(solution:Handle<Solution>)->Dynamic,|ctx,solution|{solution.with(ctx,|s|Ok(Dynamic(api::record(vec![("field",payload(&s.frame,s.data.u.clone(),vec![s.data.u.len()])?.0),("residuals",s.data.residuals.clone().encode()?),("converged",s.data.converged.encode()?)]))))});
object_binding!(helmholtz,"num_grid_solve_helmholtz_multigrid","crate::sim::multigrid::solve_helmholtz_multigrid","Actual native continuum multigrid free solver with coherent compiled spatial source of field/length² units, retained physical boundary data, and original cycle controls.",true,(grid:Handle<Grid>,sigma:f64,source:Input,budget:usize,options:Options)->NewObject<Solution>,|ctx,grid,sigma,source,budget,options|{grid.with(ctx,|g|{let g_frame_reference_m2_0_1=g.frame.reference([-2, 0, 1])?;if g.data.len()>4096||g.data.len().checked_mul(options.0.max_cycles).is_none_or(|w|w>8_388_608){return Err(err("multigrid solve work budget exceeded"));}let d=canonical(DimExpr::Mul(vec![g.frame.dims[2].clone(),DimExpr::Pow(Box::new(g.frame.dims[0].clone()),Exp::Rat(Rational::new(-2,1)))]))?;let source=checked_snapshot(ctx,&source,budget,&vec![g.frame.dims[0].clone();3],&d)?;let f=|x|space(&source,&g.frame,x,g_frame_reference_m2_0_1);let data=mg::solve_helmholtz_multigrid(&g.data,&g.bcs.data,sigma,f,options.0).map_err(err)?;source.finish()?;g.bcs.finish_for(g.data.dim)?;Ok(NewObject::new(Solution{data,frame:g.frame.clone()},g.frame.dims.to_vec()))})});
object_binding!(poisson_nd,"num_grid_solve_poisson_multigrid_nd","crate::sim::physics_mtm::solve_poisson_multigrid_nd","Actual native ND multigrid wrapper on a zero-origin rectangular grid. Extents and positive spacings are normalized by the retained Length reference; source uses field/length². Original V/W/Fmg controls, boundary conditions, residuals and convergence status are preserved.",true,(bcs:Handle<Faces>,extents:Vec<usize>,spacing:Vec<f64>,source:Vec<f64>,sigma:f64,options:Options)->NewObject<Solution>,|ctx,bcs,extents,spacing,source,sigma,options|{grid_size(&extents)?;let n=extents.iter().product::<usize>();if n>4096||source.len()!=n||spacing.len()!=extents.len()||spacing.iter().any(|h|*h<=0.0)||n.checked_mul(options.0.max_cycles).is_none_or(|w|w>8_388_608){return Err(err("multigrid ND grid/source/work budget invalid"));}bcs.with(ctx,|b|{let data=crate::sim::physics_mtm::solve_poisson_multigrid_nd(&extents,&spacing,&source,sigma,&b.data,options.0).map_err(err)?;b.finish_for(extents.len())?;Ok(NewObject::new(Solution{data,frame:b.frame.clone()},b.frame.dims.to_vec()))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &new::BINDING,
            &levels::BINDING,
            &solve::BINDING,
            &factor::BINDING,
            &data::BINDING,
            &helmholtz::BINDING,
            &poisson_nd::BINDING,
        ],
        "numerical_vector",
    )
}
