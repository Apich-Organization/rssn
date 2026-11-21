//! Native state-dependent and fixed-delay solvers using scoped native history lookup.
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
use super::object_binding;
use super::options;
use crate::kernels::ode_geom as g;
use crate::rules::numerical_abi::native_analysis_ode_common::VectorOracle;
use crate::rules::numerical_abi::native_analysis_ode_common::reference_frame;
struct Dde {
    data: g::DdeSolution,
    frame: OdeFrame,
}
impl NativeType for Dde {
    const KIND: &'static str = "rssn/native/analysis/dde/solution/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.t.len() != self.data.y.len()
            || self
                .data
                .t
                .len()
                .checked_mul(self.frame.state_scales.len())
                .is_none_or(|n| n > 47900)
        {
            return Err(api::error("DDE retained sample budget exceeded"));
        }
        Ok(())
    }
}
fn state_oracle(
    ctx: &CallCtx<'_>,
    frame: &OdeFrame,
    rhs: &Input,
    rhs_params: &Input,
    lags: usize,
) -> Result<VectorOracle, OpError> {
    let mut scales = vec![frame.time_scale];
    let mut dims = vec![frame.time_dimension.clone()];
    for _ in 0..=lags {
        scales.extend(frame.state_scales.clone());
        dims.extend(frame.state_dimensions.clone());
    }
    let rates = frame
        .state_scales
        .iter()
        .map(|s| api::scale_product(*s, &[], &[frame.time_scale]))
        .collect::<Result<Vec<_>, _>>()?;
    let rate_dims = frame
        .state_dimensions
        .iter()
        .map(|d| {
            crate::kernel::signature::DimExpr::Mul(vec![
                d.clone(),
                crate::kernel::signature::DimExpr::Inverse(Box::new(frame.time_dimension.clone())),
            ])
        })
        .collect::<Vec<_>>();
    VectorOracle::new(ctx, rhs, rhs_params, &scales, &dims, &rates, &rate_dims)
}
fn history_oracle(
    ctx: &CallCtx<'_>,
    frame: &OdeFrame,
    history: &Input,
    params: &Input,
) -> Result<VectorOracle, OpError> {
    VectorOracle::new(
        ctx,
        history,
        params,
        &[frame.time_scale],
        std::slice::from_ref(&frame.time_dimension),
        &frame.state_scales,
        &frame.state_dimensions,
    )
}
object_binding!(state_dependent,"num_analysis_dde_state_dependent","crate::kernels::ode_geom::solve_dde_state_dependent","Actual native state-dependent DDE. Each query tape(time,state)->physical past time invokes the native scoped history lookup at that stage; RHS tape receives(time,state,queried states in declared order). No lookup closure escapes or is retained. Independent state/history/rate references and explicit query/RHS evaluation budgets are preserved.",false,(time_reference:Input,state_references:Many<Input>,rhs:Input,rhs_params:Input,history:Input,history_params:Input,queries:Many<Input>,query_params:Many<Input>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Dde>,|ctx,time_reference,state_references,rhs,rhs_params,history,history_params,queries,query_params,initial,start,end,config,maxeval|{let frame=reference_frame(&time_reference,&state_references)?;if queries.0.len()>8||queries.0.len()!=query_params.0.len(){return Err(api::error("DDE query/parameter count mismatch or budget exceeds8"));}let rhs=state_oracle(ctx,&frame,&rhs,&rhs_params,queries.0.len())?;let history=history_oracle(ctx,&frame,&history,&history_params)?;let mut input_scales=vec![frame.time_scale];input_scales.extend(frame.state_scales.clone());let input_dims=frame.dimensions();let queries=queries.0.iter().zip(&query_params.0).map(|(q,p)|VectorOracle::new(ctx,q,p,&input_scales,&input_dims,&[frame.time_scale],std::slice::from_ref(&frame.time_dimension))).collect::<Result<Vec<_>,_>>()?;let y=frame.state(&initial)?;let o=options(&config,&frame)?;let options=g::DdeOptions{rtol:o.rtol,atol:o.atol,h0:o.h0.unwrap_or(0.0),h_min:o.hmin,h_max:o.hmax,max_steps:o.max_steps,max_evaluations:maxeval};let span=(frame.time(&start)?,frame.time(&end)?);let budget=Budget::new(maxeval)?;let mut rate=|t,y:&[f64],lookup:&dyn Fn(f64)->Vec<f64>|{let mut input=std::iter::once(t).chain(y.iter().copied()).collect::<Vec<_>>();let query_input=input.clone();for query in &queries{let q=query.sample(&budget,&query_input);let state=lookup(q[0]);if state.len()!=y.len(){budget.fail(api::error("native DDE lookup/state size mismatch"));return vec![f64::NAN;y.len()];}input.extend(state);}rhs.sample(&budget,&input)};let mut past=|t|history.sample(&budget,&[t]);let result=g::solve_dde_state_dependent(&mut rate,&mut past,&y,span,options);budget.finish()?;Ok(NewObject::new(Dde{data:result.map_err(api::error)?,frame:frame.clone()},frame.dimensions()))});
object_binding!(eval,"num_analysis_dde_eval","crate::kernels::ode_geom::DdeSolution::eval","Actual native DDE Hermite history evaluation, restoring independent physical component units; outside the integrated interval returns None.",false,(solution:Handle<Dde>,time:Input)->Maybe<Dynamic>,|ctx,solution,time|solution.with(ctx,|s|Ok(Maybe(s.data.eval(s.frame.time(&time)?).map(|y|s.frame.restore(&y).map(Dynamic)).transpose()?))));
object_binding!(fixed_rk4,"num_analysis_dde_fixed_rk4","crate::kernels::ode_adaptive::dde_rk4","Actual native fixed-delay RK4 method of steps. History tape(time)->state sets the initial condition; RHS tape receives(time,state,delayed state). tau/h are physical times, tau>=h, and the final native step clips the endpoint. Returned samples carry true physical coordinates, not an invented dense interpolant.",false,(time_reference:Input,state_references:Many<Input>,rhs:Input,rhs_params:Input,history:Input,history_params:Input,tau:Input,start:Input,end:Input,h:Input,maxeval:usize)->Many<Dynamic>,|ctx,time_reference,state_references,rhs,rhs_params,history,history_params,tau,start,end,h,maxeval|{let frame=reference_frame(&time_reference,&state_references)?;let rhs=state_oracle(ctx,&frame,&rhs,&rhs_params,1)?;let history=history_oracle(ctx,&frame,&history,&history_params)?;let t0=frame.time(&start)?;let t1=frame.time(&end)?;let tau=frame.time(&tau)?;let h=frame.time(&h)?;if h<=0.0||tau<h||t1<t0||((t1-t0)/h).ceil()*(frame.state_scales.len()+4)as f64>200000.0{return Err(api::error("fixed DDE span/delay/step output budget invalid"));}let budget=Budget::new(maxeval)?;let rate=|t,y:&[f64],lag:&[f64]|rhs.sample(&budget,&std::iter::once(t).chain(y.iter().copied()).chain(lag.iter().copied()).collect::<Vec<_>>());let past=|t|history.sample(&budget,&[t]);let result=a::dde_rk4(rate,past,tau,t0,t1,h);budget.finish()?;let(times,states)=result.map_err(api::error)?;Ok(Many(times.into_iter().zip(states).map(|(t,y)|Ok(Dynamic(api::record(vec![("time",api::quantity(api::scale_product(t,&[frame.time_scale],&[])?,frame.time_dimension.clone())?),("state",frame.restore(&y)?)])))).collect::<Result<_,OpError>>()?))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &state_dependent::BINDING,
            &eval::BINDING,
            &fixed_rk4::BINDING,
        ],
        "numerical_vector",
    )
}
