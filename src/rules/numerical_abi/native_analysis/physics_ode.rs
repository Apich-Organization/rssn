//! Native generic ODE-system instantiations on immutable physical tape models.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_analysis_ode_common::Budget;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeFrame;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeTape;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::Many;
use crate::rules::numerical_abi::native_object_api::Maybe;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use crate::sim::integrators as i;
use crate::sim::physics_em as e;
use crate::sim::physics_rkm as r;
use r::OdeSystem;
use std::sync::Arc;
#[path = "physics_mechanics.rs"]
mod mechanics;
pub(in crate::rules) use mechanics::Acceleration;

const HISTORY_ELEMENT_LIMIT: usize = 262_000;
const ORACLE_ELEMENT_LIMIT: usize = 1_000_000;

struct System {
    model: OdeTape,
    budget: Arc<Budget>,
}
impl r::OdeSystem for System {
    fn dim(&self) -> usize {
        self.model.frame.state_scales.len()
    }

    fn eval(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        self.model.rate(&self.budget, t, y, out);
    }
}
impl NativeType for r::DormandPrince54 {
    const KIND: &'static str = "rssn/native/analysis/rkm/dormand-prince/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
fn output(
    frame: &OdeFrame,
    rows: Vec<(f64, Vec<f64>)>,
) -> Result<Many<Dynamic>, OpError> {
    if rows
        .len()
        .checked_mul(frame.state_scales.len() + 4)
        .is_none_or(|n| n > HISTORY_ELEMENT_LIMIT)
    {
        return Err(api::error("native ODE history element budget exceeded"));
    }
    Ok(Many(
        rows.into_iter()
            .map(|(t, y)| {
                Ok(Dynamic(api::record(vec![
                    (
                        "time",
                        api::quantity(
                            api::scale_product(t, &[frame.time_scale], &[])?,
                            frame.time_dimension.clone(),
                        )?,
                    ),
                    ("state", frame.restore(&y)?),
                ])))
            })
            .collect::<Result<_, OpError>>()?,
    ))
}
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn run(
    ctx: &CallCtx<'_>,
    rhs: Handle<OdeTape>,
    initial: Input,
    start: Input,
    end: Input,
    h: Input,
    rtol: f64,
    atol: f64,
    maxeval: usize,
    kind: usize,
) -> Result<Many<Dynamic>, OpError> {
    let model = rhs.with(ctx, |m| Ok(m.clone()))?;
    let frame = model.frame.clone();
    let y = frame.state(&initial)?;
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let h = frame.time(&h)?;
    if t1 <= t0
        || h <= 0.0
        || ((t1 - t0) / h).ceil() * (y.len() + 4) as f64 > HISTORY_ELEMENT_LIMIT as f64
        || !rtol.is_finite()
        || rtol <= 0.0
        || !atol.is_finite()
        || atol <= 0.0
        || maxeval
            .checked_mul(y.len() + 4)
            .is_none_or(|n| n > ORACLE_ELEMENT_LIMIT)
    {
        return Err(api::error(
            "native ODE span/step/tolerance or worst-case output work budget invalid",
        ));
    }
    let system = System {
        model,
        budget: Budget::new(maxeval)?,
    };
    let rows = match kind {
        | 0 => r::solve_rk4(&system, &y, (t0, t1), h),
        | 1 => r::DormandPrince54::new().solve(&system, &y, (t0, t1), h, (rtol, atol)),
        | 2 => r::CashKarp45::default().solve(&system, &y, (t0, t1), h, (rtol, atol)),
        | 3 => r::BogackiShampine23::default().solve(&system, &y, (t0, t1), h, (rtol, atol)),
        | 4 => e::solve_forward_euler(&system, &y, (t0, t1), h),
        | 5 => e::solve_heun_euler(&system, &y, (t0, t1), h),
        | 6 => e::solve_midpoint_euler(&system, &y, (t0, t1), h),
        | _ => return Err(api::error("unknown native ODE solver")),
    };
    system.budget.finish()?;
    if rows.last().is_none_or(|(t, _)| {
        if kind == 1 {
            // Dormand–Prince stores the requested end exactly on acceptance.
            // An absolute tolerance can hide an incomplete short interval.
            !crate::kernels::exact_eq(*t, t1)
        } else {
            (*t - t1).abs() > 1e-12 * (1.0 + t1.abs())
        }
    }) {
        return Err(api::error(
            "native adaptive solver stopped before requested endpoint",
        ));
    }
    output(&frame, rows)
}
object_binding!(dormand_new,"num_analysis_physics_dormand_new","crate::sim::physics_rkm::DormandPrince54::new","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,()->NewObject<r::DormandPrince54>,|_ctx|Ok(NewObject::normalized(r::DormandPrince54::new())));
object_binding!(rk4,"num_analysis_physics_rk4","crate::sim::physics_rkm::solve_rk4","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,0));
object_binding!(dormand,"num_analysis_physics_dormand","crate::sim::physics_rkm::DormandPrince54::solve","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,1));
object_binding!(cash_karp,"num_analysis_physics_cash_karp","crate::sim::physics_rkm::CashKarp45::solve","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,2));
object_binding!(bogacki,"num_analysis_physics_bogacki","crate::sim::physics_rkm::BogackiShampine23::solve","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,3));
object_binding!(forward_euler,"num_analysis_physics_forward_euler","crate::sim::physics_em::solve_forward_euler","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,4));
object_binding!(heun_euler,"num_analysis_physics_heun_euler","crate::sim::physics_em::solve_heun_euler","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,5));
object_binding!(midpoint_euler,"num_analysis_physics_midpoint_euler","crate::sim::physics_em::solve_midpoint_euler","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,h:Input,rtol:f64,atol:f64,maxeval:usize)->Many<Dynamic>,|ctx,rhs,initial,start,end,h,rtol,atol,maxeval|run(ctx,rhs,initial,start,end,h,rtol,atol,maxeval,6));
object_binding!(
    explicit_euler, "num_analysis_physics_explicit_euler", "crate::sim::physics_em::solve_explicit_euler",
    "Evolves a native measured ODE tape using Forward, Midpoint or Heun with a signed time step. The last interval is shortened to reach the requested endpoint. Returns time and mixed physical state for each accepted grid point. Forward needs one tape evaluation per step and Midpoint and Heun need two.",
    false,
    (rhs: Handle<OdeTape>, initial: Input, start: Input, end: Input, step: Input, method: String, maxeval: usize) -> Many<Dynamic>,
    |ctx, rhs, initial, start, end, step, method, maxeval| {
        let method = match method.as_str() {
            | "Forward" => e::ExplicitEulerMethod::Forward,
            | "Midpoint" => e::ExplicitEulerMethod::Midpoint,
            | "Heun" => e::ExplicitEulerMethod::Heun,
            | _ => return Err(api::error("explicit Euler method must be Forward, Midpoint or Heun")),
        };
        let model = rhs.with(ctx, |model| Ok(model.clone()))?;
        let frame = model.frame.clone();
        let state = frame.state(&initial)?;
        let span = (frame.time(&start)?, frame.time(&end)?);
        let dt = frame.time(&step)?;
        let steps = e::explicit_euler_step_bound(span, dt).map_err(api::error)?;
        let row_size = state.len() + 4;
        if steps.checked_add(1).and_then(|n| n.checked_mul(row_size)).is_none_or(|n| n > HISTORY_ELEMENT_LIMIT)
            || maxeval.checked_mul(row_size).is_none_or(|n| n > ORACLE_ELEMENT_LIMIT)
        {
            return Err(api::error("explicit Euler history or oracle work budget exceeded"));
        }
        let system = System { model, budget: Budget::new(maxeval)? };
        let result = e::solve_explicit_euler(&system, &state, span, dt, method);
        system.budget.finish()?;
        output(&frame, result.map_err(api::error)?)
    }
);
object_binding!(system_dim,"num_analysis_physics_system_dim","crate::sim::physics_rkm::OdeSystem::dim","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>)->usize,|ctx,rhs|rhs.with(ctx,|m|{let system=System{model:m.clone(),budget:Budget::new(1)?};Ok(system.dim())}));
object_binding!(system_eval,"num_analysis_physics_system_eval","crate::sim::physics_rkm::OdeSystem::eval","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,time:Input,state:Input,maxeval:usize)->Dynamic,|ctx,rhs,time,state,maxeval|rhs.with(ctx,|m|{let t=m.frame.time(&time)?;let y=m.frame.state(&state)?;let system=System{model:m.clone(),budget:Budget::new(maxeval)?};let mut out=vec![0.0;y.len()];system.eval(t,&y,&mut out);system.budget.finish()?;let frame=&m.frame;let data=out.iter().zip(&frame.state_scales).zip(&frame.state_dimensions).map(|((v,s),d)|api::quantity(api::scale_product(*v,&[*s],&[frame.time_scale])?,crate::kernel::signature::DimExpr::Mul(vec![d.clone(),crate::kernel::signature::DimExpr::Inverse(Box::new(frame.time_dimension.clone()))]))).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![data.len()],elem:api::dynamic_type(),data})) )}));
object_binding!(backward_euler_step,"num_analysis_physics_backward_euler_step","crate::sim::integrators::backward_euler_step","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,maxeval:usize)->Maybe<Dynamic>,|ctx,rhs,initial,time,h,maxeval|rhs.with(ctx,|m|{let t=m.frame.time(&time)?;let h=m.frame.time(&h)?;let y=m.frame.state(&initial)?;let budget=Budget::new(maxeval)?;let f=|t,y:&[f64]|{let mut out=vec![0.0;y.len()];m.rate(&budget,t,y,&mut out);out};let result=i::backward_euler_step(&f,t,&y,h);budget.finish()?;Ok(Maybe(result.map(|y|m.frame.restore(&y).map(Dynamic)).transpose()?))}));
object_binding!(implicit_midpoint_step,"num_analysis_physics_implicit_midpoint_step","crate::sim::integrators::implicit_midpoint_step","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,time:Input,h:Input,maxeval:usize)->Maybe<Dynamic>,|ctx,rhs,initial,time,h,maxeval|rhs.with(ctx,|m|{let t=m.frame.time(&time)?;let h=m.frame.time(&h)?;let y=m.frame.state(&initial)?;let budget=Budget::new(maxeval)?;let f=|t,y:&[f64]|{let mut out=vec![0.0;y.len()];m.rate(&budget,t,y,&mut out);out};let result=i::implicit_midpoint_step(&f,t,&y,h);budget.finish()?;Ok(Maybe(result.map(|y|m.frame.restore(&y).map(Dynamic)).transpose()?))}));
object_binding!(rosenbrock_solve,"num_analysis_physics_rosenbrock_solve","crate::sim::integrators::rosenbrock_solve","Actual native generic ODE solver instantiated on an immutable tape implementing the real OdeSystem trait. Independent physical state/time references, normalized rtol/atol and a shared explicit oracle budget are preserved; incomplete native paths are errors.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,tolerance:f64,maxeval:usize)->Maybe<Many<Dynamic>>,|ctx,rhs,initial,start,end,tolerance,maxeval|rhs.with(ctx,|m|{let t0=m.frame.time(&start)?;let t1=m.frame.time(&end)?;let y=m.frame.state(&initial)?;if tolerance<=0.0||!tolerance.is_finite()||maxeval.checked_mul(y.len()+4).is_none_or(|n|n>1000000){return Err(api::error("Rosenbrock tolerance/output budget invalid"));}let budget=Budget::new(maxeval)?;let f=|t,y:&[f64]|{let mut out=vec![0.0;y.len()];m.rate(&budget,t,y,&mut out);out};let result=i::rosenbrock_solve(f,t0,&y,t1,tolerance);budget.finish()?;Ok(Maybe(result.map(|rows|output(&m.frame,rows)).transpose()?))}));
object_binding!(simulate_damped_oscillator_scenario,"num_analysis_physics_simulate_damped_oscillator_scenario","crate::sim::physics_rkm::simulate_damped_oscillator_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_rkm::simulate_damped_oscillator_scenario();if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_lorenz_attractor_scenario,"num_analysis_physics_simulate_lorenz_attractor_scenario","crate::sim::physics_rkm::simulate_lorenz_attractor_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_rkm::simulate_lorenz_attractor_scenario();if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_lotka_volterra_scenario,"num_analysis_physics_simulate_lotka_volterra_scenario","crate::sim::physics_rkm::simulate_lotka_volterra_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_rkm::simulate_lotka_volterra_scenario();if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_vanderpol_scenario,"num_analysis_physics_simulate_vanderpol_scenario","crate::sim::physics_rkm::simulate_vanderpol_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_rkm::simulate_vanderpol_scenario();if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_gravity_semi_implicit_euler_scenario,"num_analysis_physics_simulate_gravity_semi_implicit_euler_scenario","crate::sim::physics_em::simulate_gravity_semi_implicit_euler_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_em::simulate_gravity_semi_implicit_euler_scenario().map_err(api::error)?;if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_oscillator_forward_euler_scenario,"num_analysis_physics_simulate_oscillator_forward_euler_scenario","crate::sim::physics_em::simulate_oscillator_forward_euler_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_em::simulate_oscillator_forward_euler_scenario();if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
object_binding!(simulate_stiff_decay_scenario,"num_analysis_physics_simulate_stiff_decay_scenario","crate::sim::physics_em::simulate_stiff_decay_scenario","Actual native fixed built-in scenario with its literal normalized parameters. Returned normalized time/state values are explicitly dimensionless; use a physical RHS model for user-specified quantities.",false,()->Many<Dynamic>,|_ctx|{let rows=crate::sim::physics_em::simulate_stiff_decay_scenario().map_err(api::error)?;if rows.iter().map(|(_,v)|v.len()+4).sum::<usize>()>262000{return Err(api::error("built-in scenario output budget exceeded"));}Ok(Many(rows.into_iter().map(|(t,y)|Ok(Dynamic(api::record(vec![("normalized_time",Output::encode(t)?),("normalized_state",Output::encode(y)?)])))).collect::<Result<_,OpError>>()?))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    mechanics::install(i)?;
    api::install(
        i,
        &[
            &dormand_new::BINDING,
            &rk4::BINDING,
            &dormand::BINDING,
            &cash_karp::BINDING,
            &bogacki::BINDING,
            &forward_euler::BINDING,
            &heun_euler::BINDING,
            &midpoint_euler::BINDING,
            &explicit_euler::BINDING,
            &system_dim::BINDING,
            &system_eval::BINDING,
            &backward_euler_step::BINDING,
            &implicit_midpoint_step::BINDING,
            &rosenbrock_solve::BINDING,
            &simulate_damped_oscillator_scenario::BINDING,
            &simulate_lorenz_attractor_scenario::BINDING,
            &simulate_lotka_volterra_scenario::BINDING,
            &simulate_vanderpol_scenario::BINDING,
            &simulate_gravity_semi_implicit_euler_scenario::BINDING,
            &simulate_oscillator_forward_euler_scenario::BINDING,
            &simulate_stiff_decay_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
