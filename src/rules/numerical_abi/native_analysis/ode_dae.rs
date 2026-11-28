//! Native semi-explicit DAE objects with independent state/constraint reference domains.
use super::Arg;
use super::Budget;
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::Maybe;
use super::NativeType;
use super::NewObject;
use super::OdeFrame;
use super::OpError;
use super::RuleError;
use super::a;
use super::api;
use super::extended;
use super::object_binding;
use super::options;
use crate::kernels::ode_geom as g;
use crate::rules::numerical_abi::native_analysis_ode_common::VectorOracle;
use crate::rules::numerical_abi::native_analysis_ode_common::reference_frame;
#[derive(Clone)]
struct Model {
    y: OdeFrame,
    z: OdeFrame,
    f: VectorOracle,
    g: VectorOracle,
    index: usize,
}
impl NativeType for Model {
    const KIND: &'static str = "rssn/native/analysis/dae/model/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.y.state_scales.len() + self.z.state_scales.len() > 30
            || !matches!(self.index, 1 | 2)
        {
            return Err(api::error("DAE dimension/index budget exceeded"));
        }
        Ok(())
    }
}
fn inputs(
    t: f64,
    y: &[f64],
    z: &[f64],
) -> Vec<f64> {
    std::iter::once(t)
        .chain(y.iter().copied())
        .chain(z.iter().copied())
        .collect()
}
impl Model {
    fn rate(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
        z: &[f64],
    ) -> Vec<f64> {
        self.f.sample(budget, &inputs(t, y, z))
    }

    fn constraint(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
        z: &[f64],
    ) -> Vec<f64> {
        self.g
            .sample(budget, &inputs(t, y, if self.index == 1 { z } else { &[] }))
    }

    fn pair(
        &self,
        y: &[f64],
        z: &[f64],
    ) -> Result<Dynamic, OpError> {
        Ok(Dynamic(api::record(vec![
            ("differential", self.y.restore(y)?),
            ("algebraic", self.z.restore(z)?),
        ])))
    }
}
struct Dae {
    data: g::DaeSolution,
    y: OdeFrame,
    z: OdeFrame,
}
impl NativeType for Dae {
    const KIND: &'static str = "rssn/native/analysis/dae/solution/v1";
    const MAX_BYTES: usize = 12 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        let n = self.y.state_scales.len() + self.z.state_scales.len();
        if self.data.t.len().checked_mul(n).is_none_or(|n| n > 47900)
            || self.data.y.len() != self.data.t.len()
            || self.data.z.len() != self.data.t.len()
        {
            return Err(api::error("DAE retained sample budget/shape mismatch"));
        }
        Ok(())
    }
}
object_binding!(model_new,"num_analysis_dae_model_new","","Build a genuine semi-explicit native DAE model. f tape receives(time,y,z) and outputs physical dy/dt. Constraint tape receives(time,y,z) for index1or(time,y) for index2. Positive independent constraint references define its normalized Newton/violation metric; physical differential/algebraic units may differ by component.",false,(time_reference:Input,y_references:Many<Input>,z_references:Many<Input>,constraint_references:Many<Input>,rhs:Input,rhs_params:Input,constraint:Input,constraint_params:Input,index:usize)->NewObject<Model>,|ctx,time_reference,y_references,z_references,constraint_references,rhs,rhs_params,constraint,constraint_params,index|{let y=reference_frame(&time_reference,&y_references)?;let z=reference_frame(&time_reference,&z_references)?;let c=reference_frame(&time_reference,&constraint_references)?;if c.state_scales.len()!=z.state_scales.len()||!matches!(index,1|2){return Err(api::error("constraint count must equal algebraic count; index must be1or2"));}let mut scales=vec![y.time_scale];scales.extend(y.state_scales.clone());scales.extend(z.state_scales.clone());let mut dims=y.dimensions();dims.extend(z.state_dimensions.clone());let rate_scales=y.state_scales.iter().map(|s|api::scale_product(*s,&[],&[y.time_scale])).collect::<Result<Vec<_>,_>>()?;let rate_dims=y.state_dimensions.iter().map(|d|crate::kernel::signature::DimExpr::Mul(vec![d.clone(),crate::kernel::signature::DimExpr::Inverse(Box::new(y.time_dimension.clone()))])).collect::<Vec<_>>();let f=VectorOracle::new(ctx,&rhs,&rhs_params,&scales,&dims,&rate_scales,&rate_dims)?;if index==2{scales.truncate(1+y.state_scales.len());dims.truncate(scales.len());}let g=VectorOracle::new(ctx,&constraint,&constraint_params,&scales,&dims,&c.state_scales,&c.state_dimensions)?;let mut dimensions=y.dimensions();dimensions.extend(z.state_dimensions.clone());let model=Model{y,z,f,g,index};model.validate()?;Ok(NewObject::new(model,dimensions))});
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn consistent(
    ctx: &CallCtx<'_>,
    model: Handle<Model>,
    y: Input,
    z: Input,
    time: Input,
    tolerance: f64,
    max_iterations: usize,
    maxeval: usize,
    index: usize,
) -> Result<Dynamic, OpError> {
    model.with(ctx, |m| {
        if m.index != index
            || max_iterations == 0
            || max_iterations > 100
            || tolerance <= 0.0
            || !tolerance.is_finite()
        {
            return Err(api::error(
                "DAE consistency index/tolerance/iteration budget mismatch",
            ));
        }
        let y = m.y.state(&y)?;
        let z = m.z.state(&z)?;
        let t = m.y.time(&time)?;
        let budget = Budget::new(maxeval)?;
        let mut f = |t, y: &[f64], z: &[f64]| m.rate(&budget, t, y, z);
        let result = if index == 1 {
            let mut g = |t, y: &[f64], z: &[f64]| m.constraint(&budget, t, y, z);
            g::dae_consistent_ic_index1(&mut f, &mut g, t, &y, &z, tolerance, max_iterations)
        } else {
            let mut g = |t, y: &[f64]| m.constraint(&budget, t, y, &[]);
            g::dae_consistent_ic_index2(&mut f, &mut g, t, &y, &z, tolerance, max_iterations)
        };
        budget.finish()?;
        let (y, z) = result.map_err(api::error)?;
        m.pair(&y, &z)
    })
}
object_binding!(consistent_index1,"num_analysis_dae_consistent_index1","crate::kernels::ode_geom::dae_consistent_ic_index1","Actual native index1 DAE consistency projection with independently measured variables, explicit constraint metric and original oracle errors.",false,(model:Handle<Model>,y:Input,z:Input,time:Input,tolerance:f64,max_iterations:usize,maxeval:usize)->Dynamic,|ctx,model,y,z,time,tolerance,max_iterations,maxeval|consistent(ctx,model,y,z,time,tolerance,max_iterations,maxeval,1));
object_binding!(consistent_index2,"num_analysis_dae_consistent_index2","crate::kernels::ode_geom::dae_consistent_ic_index2","Actual native index2 DAE consistency projection with independently measured variables, explicit constraint metric and original oracle errors.",false,(model:Handle<Model>,y:Input,z:Input,time:Input,tolerance:f64,max_iterations:usize,maxeval:usize)->Dynamic,|ctx,model,y,z,time,tolerance,max_iterations,maxeval|consistent(ctx,model,y,z,time,tolerance,max_iterations,maxeval,2));
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn solve(
    ctx: &CallCtx<'_>,
    model: Handle<Model>,
    initial_y: Input,
    initial_z: Input,
    start: Input,
    end: Input,
    config: Input,
    maxeval: usize,
    index: usize,
) -> Result<NewObject<Dae>, OpError> {
    model.with(ctx, |m| {
        if m.index != index {
            return Err(api::error("DAE model index mismatch"));
        }
        let y = m.y.state(&initial_y)?;
        let z = m.z.state(&initial_z)?;
        let o = options(&config, &m.y)?;
        if o.max_steps
            .checked_mul(y.len() + z.len())
            .is_none_or(|n| n > 47900)
        {
            return Err(api::error("DAE retained state sample budget exceeded"));
        }
        let newton_tol = api::field(&config.0, "newton_tol")?;
        api::normalized(newton_tol)?;
        let newton_tol = <f64 as Arg>::decode(newton_tol, 0)?;
        let max_newton_iter = <usize as Arg>::decode(api::field(&config.0, "max_newton_iter")?, 0)?;
        if max_newton_iter == 0 || max_newton_iter > 100 {
            return Err(api::error("DAE Newton iteration budget must be1..100"));
        }
        let options = g::DaeOptions {
            rtol: o.rtol,
            atol: o.atol,
            h0: o.h0.unwrap_or(0.0),
            h_min: o.hmin,
            h_max: o.hmax,
            max_steps: o.max_steps,
            max_evaluations: maxeval,
            newton_tol,
            max_newton_iter,
        };
        let span = (m.y.time(&start)?, m.y.time(&end)?);
        let budget = Budget::new(maxeval)?;
        let mut f = |t, y: &[f64], z: &[f64]| m.rate(&budget, t, y, z);
        let result = if index == 1 {
            let mut g = |t, y: &[f64], z: &[f64]| m.constraint(&budget, t, y, z);
            g::solve_dae_index1(&mut f, &mut g, &y, &z, span, options)
        } else {
            let mut g = |t, y: &[f64]| m.constraint(&budget, t, y, &[]);
            g::solve_dae_index2(&mut f, &mut g, &y, &z, span, options)
        };
        budget.finish()?;
        let mut dims = m.y.dimensions();
        dims.extend(m.z.state_dimensions.clone());
        Ok(NewObject::new(
            Dae {
                data: result.map_err(api::error)?,
                y: m.y.clone(),
                z: m.z.clone(),
            },
            dims,
        ))
    })
}
object_binding!(solve_index1,"num_analysis_dae_solve_index1","crate::kernels::ode_geom::solve_dae_index1","Actual adaptive native index1 DAE solve with independent differential/algebraic physical domains, normalized error/constraint metrics, physical step bounds and explicit Newton/attempt/evaluation budgets.",false,(model:Handle<Model>,initial_y:Input,initial_z:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Dae>,|ctx,model,initial_y,initial_z,start,end,config,maxeval|solve(ctx,model,initial_y,initial_z,start,end,config,maxeval,1));
object_binding!(solve_index2,"num_analysis_dae_solve_index2","crate::kernels::ode_geom::solve_dae_index2","Actual adaptive native index2 DAE solve with independent differential/algebraic physical domains, normalized error/constraint metrics, physical step bounds and explicit Newton/attempt/evaluation budgets.",false,(model:Handle<Model>,initial_y:Input,initial_z:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Dae>,|ctx,model,initial_y,initial_z,start,end,config,maxeval|solve(ctx,model,initial_y,initial_z,start,end,config,maxeval,2));
object_binding!(eval,"num_analysis_dae_eval","crate::kernels::ode_geom::DaeSolution::eval","Actual native DAE dense evaluation restoring differential/algebraic component units; outside the integrated interval remains None.",false,(solution:Handle<Dae>,time:Input)->Maybe<Dynamic>,|ctx,solution,time|solution.with(ctx,|s|Ok(Maybe(s.data.eval(s.y.time(&time)?).map(|(y,z)|Ok::<_,OpError>(Dynamic(api::record(vec![("differential",s.y.restore(&y)?),("algebraic",s.z.restore(&z)?)])))).transpose()?))));
object_binding!(violation,"num_analysis_dae_max_constraint_violation","crate::kernels::ode_geom::DaeSolution::max_constraint_violation","Actual native maximum constraint violation in the explicit independently scaled normalized constraint metric. No common physical residual unit is invented for heterogeneous equations.",false,(solution:Handle<Dae>,model:Handle<Model>,maxeval:usize)->f64,|ctx,solution,model,maxeval|{let m=model.with(ctx,|m|Ok(m.clone()))?;solution.with(ctx,|s|{super::extended::same_frame(&m.y,&s.y)?;super::extended::same_frame(&m.z,&s.z)?;let budget=Budget::new(maxeval)?;let value=s.data.max_constraint_violation(|t,y,z|m.constraint(&budget,t,y,z));budget.finish()?;Ok(value)})});
object_binding!(fixed_index1,"num_analysis_dae_fixed_index1","crate::kernels::ode_adaptive::dae_index1","Actual native fixed backward-Euler index1DAE, with independently normalized differential/algebraic variables and physical time step. Returned states concatenate the original differential/algebraic order.",false,(model:Handle<Model>,initial_y:Input,initial_z:Input,start:Input,h:Input,steps:usize,maxeval:usize)->Many<Dynamic>,|ctx,model,initial_y,initial_z,start,h,steps,maxeval|model.with(ctx,|m|{if m.index!=1||steps==0||steps.checked_mul(m.y.state_scales.len()+m.z.state_scales.len()+4).is_none_or(|n|n>200000){return Err(api::error("fixed DAE index/sample budget mismatch"));}let y=m.y.state(&initial_y)?;let z=m.z.state(&initial_z)?;let t=m.y.time(&start)?;let h=m.y.time(&h)?;if h<=0.0{return Err(api::error("fixed DAE step must be positive"));}let budget=Budget::new(maxeval)?;let f=|t,y:&[f64],z:&[f64]|m.rate(&budget,t,y,z);let g=|t,y:&[f64],z:&[f64]|m.constraint(&budget,t,y,z);let result=a::dae_index1(f,g,&y,&z,t,h,steps);budget.finish()?;Ok(Many(result.map_err(api::error)?.iter().map(|r|m.pair(&r[..y.len()],&r[y.len()..])).collect::<Result<_,_>>()?))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &model_new::BINDING,
            &consistent_index1::BINDING,
            &consistent_index2::BINDING,
            &solve_index1::BINDING,
            &solve_index2::BINDING,
            &eval::BINDING,
            &violation::BINDING,
            &fixed_index1::BINDING,
        ],
        "numerical_vector",
    )
}
