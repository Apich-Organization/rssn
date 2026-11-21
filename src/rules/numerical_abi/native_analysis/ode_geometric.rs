//! Actual native geometric steps and schedules over physical RHS reference frames.
use super::Budget;
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::NewObject;
use super::OdeTape;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Solution;
use super::SolutionData;
use super::api;
use super::event_value;
use super::events;
use super::object_binding;
use super::options;
use crate::kernels::ode_geom as g;
fn method(value: usize) -> Result<g::GeometricMethod, OpError> {
    Ok(match value {
        | 0 => g::GeometricMethod::Rk4,
        | 1 => g::GeometricMethod::ImplicitMidpoint,
        | 2 => g::GeometricMethod::StormerVerlet,
        | 3 => g::GeometricMethod::Ruth4,
        | 4 => g::GeometricMethod::GaussLegendre4,
        | 5 => g::GeometricMethod::GaussLegendre6,
        | _ => {
            return Err(api::error(
                "geometric method must be0RK4/1midpoint/2Verlet/3Ruth4/4GL4/5GL6",
            ));
        },
    })
}
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn step(
    ctx: &CallCtx<'_>,
    rhs: Handle<OdeTape>,
    initial: Input,
    time: Input,
    h: Input,
    stages: usize,
    maxeval: usize,
    kind: usize,
) -> Result<Dynamic, OpError> {
    let m = rhs.with(ctx, |m| Ok(m.clone()))?;
    let y = m.frame.state(&initial)?;
    let t = m.frame.time(&time)?;
    let h = m.frame.time(&h)?;
    if y.len() % 2 != 0 && matches!(kind, 1 | 2) {
        return Err(api::error(
            "partitioned Verlet/Ruth require an even state split; no fallback method is selected",
        ));
    }
    let budget = Budget::new(maxeval)?;
    let mut rate = |t, y: &[f64]| {
        let mut out = vec![0.0; y.len()];
        m.rate(&budget, t, y, &mut out);
        out
    };
    let result = match kind {
        | 0 => g::rk4_step(&mut rate, t, &y, h).map(|y| (y, None)),
        | 1 => g::partitioned_verlet_step(&mut rate, t, &y, h).map(|(y, n)| (y, Some(n))),
        | 2 => g::partitioned_ruth4_step(&mut rate, t, &y, h).map(|(y, n)| (y, Some(n))),
        | 3 => g::gauss_legendre_step(&mut rate, t, &y, h, stages).map(|(y, n)| (y, Some(n))),
        | _ => return Err(api::error("unknown geometric step")),
    };
    budget.finish()?;
    let (y, n) = result.map_err(api::error)?;
    Ok(Dynamic(api::record(vec![
        ("state", m.frame.restore(&y)?),
        ("native_evaluations", Output::encode(n)?),
        ("evaluations", Output::encode(budget.used())?),
    ])))
}
object_binding!(rk4_step,"num_analysis_ode_rk4_step","crate::kernels::ode_geom::rk4_step","Actual native rk4_step. Physical time/state/rate references are retained by the RHS tape; no alternate method is substituted. stages selects the native1/2/3-stage Gauss-Legendre method where applicable.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,stages:usize,maxeval:usize)->Dynamic,|ctx,rhs,initial,time,h,stages,maxeval|step(ctx,rhs,initial,time,h,stages,maxeval,0));
object_binding!(verlet_step,"num_analysis_ode_verlet_step","crate::kernels::ode_geom::partitioned_verlet_step","Actual native partitioned_verlet_step. Physical time/state/rate references are retained by the RHS tape; no alternate method is substituted. stages selects the native1/2/3-stage Gauss-Legendre method where applicable.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,stages:usize,maxeval:usize)->Dynamic,|ctx,rhs,initial,time,h,stages,maxeval|step(ctx,rhs,initial,time,h,stages,maxeval,1));
object_binding!(ruth4_step,"num_analysis_ode_ruth4_step","crate::kernels::ode_geom::partitioned_ruth4_step","Actual native partitioned_ruth4_step. Physical time/state/rate references are retained by the RHS tape; no alternate method is substituted. stages selects the native1/2/3-stage Gauss-Legendre method where applicable.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,stages:usize,maxeval:usize)->Dynamic,|ctx,rhs,initial,time,h,stages,maxeval|step(ctx,rhs,initial,time,h,stages,maxeval,2));
object_binding!(gauss_step,"num_analysis_ode_gauss_step","crate::kernels::ode_geom::gauss_legendre_step","Actual native gauss_legendre_step. Physical time/state/rate references are retained by the RHS tape; no alternate method is substituted. stages selects the native1/2/3-stage Gauss-Legendre method where applicable.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,stages:usize,maxeval:usize)->Dynamic,|ctx,rhs,initial,time,h,stages,maxeval|step(ctx,rhs,initial,time,h,stages,maxeval,3));

#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn schedule(
    ctx: &CallCtx<'_>,
    rhs: Handle<OdeTape>,
    initial: Input,
    start: Input,
    end: Input,
    config: Input,
    controls: Vec<usize>,
    event_values: Many<Input>,
    maxeval: usize,
    adaptive: bool,
) -> Result<NewObject<Solution>, OpError> {
    let m = rhs.with(ctx, |m| Ok(m.clone()))?;
    let frame = m.frame.clone();
    let y = frame.state(&initial)?;
    let span = (frame.time(&start)?, frame.time(&end)?);
    let o = options(&config, &frame)?;
    let method = method(controls.first().copied().unwrap_or(0))?;
    if matches!(
        method,
        g::GeometricMethod::StormerVerlet | g::GeometricMethod::Ruth4
    ) && !y.len().is_multiple_of(2)
    {
        return Err(api::error(
            "partitioned geometric methods require even state split",
        ));
    }
    let steps = controls.get(1).copied().unwrap_or(1);
    if steps == 0 || steps.checked_mul(y.len()).is_none_or(|n| n > 47900) {
        return Err(api::error("geometric schedule sample budget exceeded"));
    }
    let budget = Budget::new(maxeval)?;
    let programs = events(ctx, &event_values, &frame, maxeval)?;
    let specs = programs
        .iter()
        .map(|p| {
            g::EventSpec {
                terminal: p.terminal,
                direction: p.direction,
            }
        })
        .collect::<Vec<_>>();
    let mut rate = |t, y: &[f64]| {
        let mut out = vec![0.0; y.len()];
        m.rate(&budget, t, y, &mut out);
        out
    };
    let mut event = |index, t, y: &[f64]| event_value(&programs[index], &frame, &budget, t, y);
    let invalid = std::cell::Cell::new(false);
    let result = if adaptive {
        g::solve_adaptive_geometric(
            &mut rate,
            &y,
            span,
            method,
            g::AdaptiveGeomOptions {
                rtol: o.rtol,
                atol: o.atol,
                h0: o.h0.unwrap_or(0.0),
                h_min: o.hmin,
                h_max: o.hmax,
                max_steps: o.max_steps,
                max_evaluations: maxeval,
            },
            &specs,
            &mut event,
            &invalid,
        )
    } else {
        g::solve_fixed_geometric(
            &mut rate,
            &y,
            span,
            steps,
            method,
            controls.get(2).copied(),
            &specs,
            &mut event,
            &invalid,
        )
    };
    budget.finish()?;
    for p in &programs {
        p.program.finish()?;
    }
    let data = SolutionData::Geometric(result.map_err(api::error)?);
    Ok(NewObject::new(
        Solution {
            frame: frame.clone(),
            data,
            evaluations: budget.used(),
        },
        frame.dimensions(),
    ))
}
object_binding!(fixed,"num_analysis_ode_geometric_fixed","crate::kernels::ode_geom::solve_fixed_geometric","Actual fixed geometric driver with physical state components, native Hermite event localization, explicit method/steps/optional stages controls and budgets.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|schedule(ctx,rhs,initial,start,end,config,controls,events,maxeval,false));
object_binding!(adaptive,"num_analysis_ode_geometric_adaptive","crate::kernels::ode_geom::solve_adaptive_geometric","Actual tolerance-controlled native geometric driver; method-specific full/two-half-step error control, physical step bounds, native events and explicit attempt/evaluation budgets.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|schedule(ctx,rhs,initial,start,end,config,controls,events,maxeval,true));
object_binding!(sample,"num_analysis_ode_geometric_sample","crate::kernels::ode_geom::sample_fixed_geometric","Actual native fixed-geometric requested-time sampler. Times have the RHS time dimension; mixed physical state components retain independent references.",false,(rhs:Handle<OdeTape>,initial:Input,times:Many<Input>,method_index:usize,steps:usize,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,times,method_index,steps,maxeval|{let m=rhs.with(ctx,|m|Ok(m.clone()))?;let y=m.frame.state(&initial)?;if times.0.is_empty()||times.0.len().checked_mul(y.len()).is_none_or(|n|n>47900)||steps==0||steps>47900{return Err(api::error("geometric sampler work/storage budget exceeded"));}let method=method(method_index)?;if matches!(method,g::GeometricMethod::StormerVerlet|g::GeometricMethod::Ruth4)&&!y.len().is_multiple_of(2){return Err(api::error("partitioned method requires even state split"));}let ts=times.0.iter().map(|t|m.frame.time(t)).collect::<Result<Vec<_>,_>>()?;let budget=Budget::new(maxeval)?;let mut rate=|t,y:&[f64]|{let mut out=vec![0.0;y.len()];m.rate(&budget,t,y,&mut out);out};let result=g::sample_fixed_geometric(&mut rate,&y,&ts,method,steps);budget.finish()?;Ok(Many(result.map_err(api::error)?.iter().map(|y|m.frame.restore(y).map(Dynamic)).collect::<Result<Vec<_>,_>>()?))});
object_binding!(hermite,"num_analysis_ode_hermite_segment","crate::kernels::ode_geom::hermite_segment","Native Hermite segment; endpoints are physical state vectors and endpoint rates have each state dimension divided by the RHS time dimension.",false,(rhs:Handle<OdeTape>,start:Input,end:Input,y0:Input,y1:Input,f0:Input,f1:Input,time:Input)->Dynamic,|ctx,rhs,start,end,y0,y1,f0,f1,time|{let m=rhs.with(ctx,|m|Ok(m.clone()))?;let frame=&m.frame;let t0=frame.time(&start)?;let t1=frame.time(&end)?;let y0=frame.state(&y0)?;let y1=frame.state(&y1)?;let dims=frame.state_dimensions.iter().map(|d|crate::kernel::signature::DimExpr::Mul(vec![d.clone(),crate::kernel::signature::DimExpr::Inverse(Box::new(frame.time_dimension.clone()))])).collect::<Vec<_>>();let scales=frame.state_scales.iter().map(|s|api::scale_product(*s,&[],&[frame.time_scale])).collect::<Result<Vec<_>,_>>()?;let f0=crate::rules::numerical_abi::native_analysis_newton_cg::parameters(&f0,&scales,&dims)?;let f1=crate::rules::numerical_abi::native_analysis_newton_cg::parameters(&f1,&scales,&dims)?;let t=frame.time(&time)?;Ok(Dynamic(frame.restore(&g::hermite_segment(t0,t1,&y0,&y1,&f0,&f1,t).map_err(api::error)?)?))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &rk4_step::BINDING,
            &verlet_step::BINDING,
            &ruth4_step::BINDING,
            &gauss_step::BINDING,
            &fixed::BINDING,
            &adaptive::BINDING,
            &sample::BINDING,
            &hermite::BINDING,
        ],
        "numerical_vector",
    )
}
