//! Native ODE trajectories and dense solutions in independently scaled physical coordinates.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernels::ode;
use crate::kernels::ode_adaptive as a;
use crate::kernels::ode_sdirk as s;
use crate::kernels::ode_stiff as stiff;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_analysis_ode_common::Budget;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeFrame;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeTape;
use crate::rules::numerical_abi::native_backend::CompiledSnapshot;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
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
use std::sync::Arc;
type EventFunction<'a> = Box<dyn Fn(f64, &[f64]) -> f64 + 'a>;
#[path = "ode_boundary.rs"]
mod boundary;
#[path = "ode_dae.rs"]
mod dae;
#[path = "ode_dde.rs"]
mod dde;
#[path = "ode_extended.rs"]
mod extended;
#[path = "ode_geometric.rs"]
mod geometric;
#[path = "ode_hamiltonian.rs"]
mod hamiltonian;
struct Solution {
    frame: OdeFrame,
    data: SolutionData,
    evaluations: usize,
}
enum SolutionData {
    Trajectory(ode::Trajectory),
    FixedEvents(ode::Trajectory, Vec<ode::FixedEvent>),
    Geometric(crate::kernels::ode_geom::FixedGeomSolution),
    Dense(a::OdeSolution),
    Sdirk(s::SdirkSolution),
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/analysis/ode/solution/v1";
    // Bound includes native dense polynomial vectors and their per-step allocation headers.
    const MAX_BYTES: usize = 24 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        let (t, y) = self.nodes();
        if t.len() != y.len()
            || t.len()
                .checked_mul(self.frame.state_scales.len())
                .is_none_or(|n| n > 48000)
            || y.iter().any(|r| {
                r.len() != self.frame.state_scales.len() || r.iter().any(|v| !v.is_finite())
            })
        {
            return Err(api::error(
                "invalid ODE solution or retained sample budget exceeded",
            ));
        }
        Ok(())
    }
}
impl Solution {
    fn nodes(&self) -> (&[f64], &[Vec<f64>]) {
        match &self.data {
            | SolutionData::Trajectory(x) | SolutionData::FixedEvents(x, _) => (&x.t, &x.y),
            | SolutionData::Dense(x) => (&x.t, &x.y),
            | SolutionData::Sdirk(x) => (&x.t, &x.y),
            | SolutionData::Geometric(x) => (&x.t, &x.y),
        }
    }
}
fn options(
    v: &Input,
    frame: &OdeFrame,
) -> Result<a::OdeOptions, OpError> {
    let f = |name| {
        let value = api::field(&v.0, name)?;
        api::normalized(value)?;
        <f64 as Arg>::decode(value, 0)
    };
    let time = |name| frame.time(&Input(api::field(&v.0, name)?.clone()));
    let rtol = f("rtol")?;
    let atol = f("atol")?;
    let h0 = time("h0")?;
    let hmin = time("hmin")?;
    let hmax = time("hmax")?;
    let max_steps = <usize as Arg>::decode(api::field(&v.0, "max_steps")?, 0)?;
    if !rtol.is_finite()
        || rtol < 0.0
        || !atol.is_finite()
        || atol < 0.0
        || (rtol <= 0.0 && atol <= 0.0)
        || h0 < 0.0
        || hmin < 0.0
        || hmax <= 0.0
        || hmin > hmax
        || max_steps == 0
        || max_steps
            .checked_mul(frame.state_scales.len())
            .is_none_or(|n| n > 47900)
    {
        return Err(api::error(
            "ODE tolerances must be positive in the normalized state metric; physical step bounds and retained sample budget must be valid",
        ));
    }
    Ok(a::OdeOptions {
        rtol,
        atol,
        h0: if h0 > 0.0 { Some(h0) } else { None },
        hmin,
        hmax,
        max_steps,
    })
}
fn sconfig(
    v: &Input,
    o: &a::OdeOptions,
) -> Result<s::SdirkOptions, OpError> {
    if o.hmin != 0.0 {
        return Err(api::error(
            "native SDIRK does not expose hmin; specify zero",
        ));
    }
    let newton_max_iter = <usize as Arg>::decode(api::field(&v.0, "newton_max_iter")?, 0)?;
    let newton_tol = <f64 as Arg>::decode(api::field(&v.0, "newton_tol")?, 0)?;
    if newton_max_iter == 0 || newton_max_iter > 100 || !newton_tol.is_finite() || newton_tol <= 0.0
    {
        return Err(api::error("invalid SDIRK Newton tolerance/budget"));
    }
    Ok(s::SdirkOptions {
        rtol: o.rtol,
        atol: o.atol,
        h0: o.h0,
        hmax: o.hmax,
        max_steps: o.max_steps,
        newton_max_iter,
        newton_tol,
    })
}
pub(in crate::rules) struct EventProgram {
    pub(in crate::rules) program: CompiledSnapshot,
    pub(in crate::rules) terminal: bool,
    pub(in crate::rules) direction: i8,
}
pub(in crate::rules) fn events(
    ctx: &CallCtx<'_>,
    values: &Many<Input>,
    frame: &OdeFrame,
    budget: usize,
) -> Result<Vec<EventProgram>, OpError> {
    if values.0.len() > 32 {
        return Err(api::error("ODE event budget exceeds32"));
    }
    values
        .0
        .iter()
        .map(|v| {
            let program = compiled_snapshot(ctx, api::field(&v.0, "program")?, budget)?;
            if program.input_dimensions.len() != frame.state_scales.len() + 1 {
                return Err(api::error("event requires time and each state component"));
            }
            api::same_dimension(&program.input_dimensions[0], &frame.time_dimension)?;
            for (i, d) in frame.state_dimensions.iter().enumerate() {
                api::same_dimension(d, &program.input_dimensions[i + 1])?;
            }
            let direction = <i64 as Arg>::decode(api::field(&v.0, "direction")?, 0)?;
            if !(-1..=1).contains(&direction) {
                return Err(api::error("event direction must be -1,0,1"));
            }
            Ok(EventProgram {
                program,
                terminal: <bool as Arg>::decode(api::field(&v.0, "terminal")?, 0)?,
                direction: direction as i8,
            })
        })
        .collect()
}
pub(in crate::rules) fn event_value(
    event: &EventProgram,
    frame: &OdeFrame,
    budget: &Budget,
    t: f64,
    y: &[f64],
) -> f64 {
    if !budget.charge() {
        return f64::NAN;
    }
    let x = std::iter::once(api::scale_product(t, &[frame.time_scale], &[]))
        .chain(
            y.iter()
                .zip(&frame.state_scales)
                .map(|(v, s)| api::scale_product(*v, &[*s], &[])),
        )
        .collect::<Result<Vec<_>, _>>();
    match x.and_then(|x| event.program.evaluate_coherent(&x)) {
        | Ok(v) => v,
        | Err(e) => {
            budget.fail(e);
            f64::NAN
        },
    }
}
#[allow(deprecated)] // The legacy native entry remains intentionally callable.
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
    config: Input,
    controls: Vec<usize>,
    event_values: Many<Input>,
    maxeval: usize,
    kind: usize,
) -> Result<NewObject<Solution>, OpError> {
    let model = Arc::new(rhs.with(ctx, |x| Ok(x.clone()))?);
    let frame = model.frame.clone();
    let y = frame.state(&initial)?;
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let o = options(&config, &frame)?;
    // This legacy kernel has no step-bound arguments. Reject constraints it
    // cannot honor instead of silently accepting the shared option fields.
    if kind == 1 && (o.h0.is_some() || o.hmin > 0.0 || o.hmax < (t1 - t0).abs()) {
        return Err(api::error(
            "legacy adaptive driver does not expose h0/hmin/hmax: require h0=0,hmin=0,hmax covering the complete span; use the dedicated adaptive drivers for tighter bounds",
        ));
    }
    let steps = controls.first().copied().unwrap_or(1);
    if steps == 0 || steps.checked_mul(y.len()).is_none_or(|n| n > 47900) {
        return Err(api::error("ODE step/control storage budget exceeded"));
    }
    let budget = Budget::new(maxeval)?;
    if !event_values.0.is_empty() && matches!(kind, 0 | 1 | 4 | 9 | 11 | 12 | 13 | 14 | 16) {
        return Err(api::error(
            "this native entry has no event argument; select its events variant",
        ));
    }
    let programs = events(ctx, &event_values, &frame, maxeval)?;
    let functions: Vec<EventFunction<'_>> = programs
        .iter()
        .map(|event| {
            let frame = &frame;
            let budget = &budget;
            Box::new(move |t, y: &[f64]| event_value(event, frame, budget, t, y))
                as Box<dyn Fn(f64, &[f64]) -> f64>
        })
        .collect();
    let es: Vec<a::Event<'_>> = programs
        .iter()
        .zip(&functions)
        .map(|(e, f)| {
            a::Event {
                g: f.as_ref(),
                terminal: e.terminal,
                direction: e.direction,
            }
        })
        .collect();
    let rate = |t: f64, y: &[f64], out: &mut [f64]| model.rate(&budget, t, y, out);
    let jac = |t: f64, y: &[f64]| model.jacobian(&budget, t, y);
    let result: Result<SolutionData, OpError> = (|| {
        Ok(match kind {
            | 0 => {
                SolutionData::Trajectory(
                    ode::solve_fixed(
                        rate,
                        &y,
                        (t0, t1),
                        steps,
                        match controls.get(1).copied().unwrap_or(2) {
                            | 0 => ode::OdeSolverMethod::Euler,
                            | 1 => ode::OdeSolverMethod::Heun,
                            | 2 => ode::OdeSolverMethod::RungeKutta4,
                            | _ => {
                                return Err(api::error(
                                    "fixed ODE method must be0Euler/1Heun/2RK4",
                                ));
                            },
                        },
                    )
                    .map_err(api::error)?,
                )
            },
            | 1 => {
                SolutionData::Trajectory(
                    ode::solve_adaptive(rate, &y, (t0, t1), o.rtol, o.atol, o.max_steps)
                        .map_err(api::error)?,
                )
            },
            | 2 => {
                let (traj, hits) = ode::solve_fixed_events(
                    rate,
                    &y,
                    (t0, t1),
                    steps,
                    match controls.get(1).copied().unwrap_or(2) {
                        | 0 => ode::OdeSolverMethod::Euler,
                        | 1 => ode::OdeSolverMethod::Heun,
                        | 2 => ode::OdeSolverMethod::RungeKutta4,
                        | _ => return Err(api::error("fixed ODE method must be0Euler/1Heun/2RK4")),
                    },
                    &es,
                )
                .map_err(api::error)?;
                SolutionData::FixedEvents(traj, hits)
            },
            | 3 => {
                SolutionData::Dense(
                    stiff::bdf_events(rate, t0, t1, &y, &o, steps, &es).map_err(api::error)?,
                )
            },
            | 4 => {
                SolutionData::Dense(
                    stiff::bdf_jac(rate, jac, t0, t1, &y, &o, steps).map_err(api::error)?,
                )
            },
            | 5 => {
                SolutionData::Dense(
                    stiff::bdf_ex(rate, Some(&jac), t0, t1, &y, &o, steps, &es)
                        .map_err(api::error)?,
                )
            },
            | 6 | 7 => {
                let bo = stiff::BdfOptions {
                    jac: if controls.get(1) == Some(&1) {
                        Some(&jac)
                    } else {
                        None
                    },
                    max_order: steps,
                    events: &es,
                    ..stiff::BdfOptions::default()
                };
                SolutionData::Dense(
                    if kind == 6 {
                        stiff::bdf_with(rate, t0, t1, &y, &o, &bo)
                    } else {
                        stiff::bdf_with_par(rate, t0, t1, &y, &o, &bo)
                    }
                    .map_err(api::error)?,
                )
            },
            | 8 => {
                SolutionData::Dense(
                    stiff::radau5_events(rate, t0, t1, &y, &o, &es).map_err(api::error)?,
                )
            },
            | 9 => {
                SolutionData::Dense(
                    stiff::radau5_jac(rate, jac, t0, t1, &y, &o).map_err(api::error)?,
                )
            },
            | 10 => {
                let ro = stiff::RadauOptions {
                    jac: if controls.get(2) == Some(&1) {
                        Some(&jac)
                    } else {
                        None
                    },
                    min_stages: steps,
                    max_stages: controls.get(1).copied().unwrap_or(steps),
                    events: &es,
                    ..stiff::RadauOptions::default()
                };
                SolutionData::Dense(
                    stiff::radau_with(rate, t0, t1, &y, &o, &ro).map_err(api::error)?,
                )
            },
            | 11 => {
                let last = stiff::radau5_fixed(rate, t0, t1, &y, steps).map_err(api::error)?;
                SolutionData::Trajectory(ode::Trajectory {
                    t: vec![t0, t1],
                    y: vec![y.clone(), last],
                })
            },
            | 12 => {
                let last = stiff::radau_fixed(
                    rate,
                    t0,
                    t1,
                    &y,
                    steps,
                    controls.get(1).copied().unwrap_or(3),
                )
                .map_err(api::error)?;
                SolutionData::Trajectory(ode::Trajectory {
                    t: vec![t0, t1],
                    y: vec![y.clone(), last],
                })
            },
            | 13 => {
                SolutionData::Sdirk(
                    s::sdirk3(rate, jac, t0, t1, &y, &sconfig(&config, &o)?).map_err(api::error)?,
                )
            },
            | 14 => {
                SolutionData::Sdirk(
                    s::sdirk3_fd(rate, t0, t1, &y, &sconfig(&config, &o)?).map_err(api::error)?,
                )
            },
            | 15 => {
                SolutionData::Dense(
                    s::sdirk3_events(rate, jac, t0, t1, &y, &sconfig(&config, &o)?, &es)
                        .map_err(api::error)?,
                )
            },
            | 16 => {
                let last = s::sdirk3_fixed(rate, jac, t0, t1, &y, steps).map_err(api::error)?;
                SolutionData::Trajectory(ode::Trajectory {
                    t: vec![t0, t1],
                    y: vec![y.clone(), last],
                })
            },
            | 17 => {
                SolutionData::Dense(
                    s::sdirk3_fixed_events(rate, jac, t0, t1, &y, steps, &es)
                        .map_err(api::error)?,
                )
            },
            | 18 => {
                let cfg = stiff::RadauConfig {
                    mass: None,
                    min_stages: steps,
                    max_stages: controls.get(1).copied().unwrap_or(steps),
                };
                SolutionData::Dense(
                    stiff::radau(rate, Some(&jac), t0, t1, &y, &o, &es, &cfg)
                        .map_err(api::error)?,
                )
            },
            | 19 => {
                SolutionData::Dense(
                    a::rosenbrock23_events(rate, t0, t1, &y, &o, &es).map_err(api::error)?,
                )
            },
            | _ => return Err(api::error("unknown native ODE entry")),
        })
    })();
    budget.finish()?;
    for e in &programs {
        e.program.finish()?;
    }
    let data = result?;
    drop(es);
    drop(functions);
    let dimensions = frame.dimensions();
    Ok(NewObject::new(
        Solution {
            frame,
            data,
            evaluations: budget.used(),
        },
        dimensions,
    ))
}

object_binding!(fixed,"num_analysis_ode_fixed","crate::kernels::ode::solve_fixed","Native ode::solve_fixed; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,0));

object_binding!(adaptive,"num_analysis_ode_adaptive","crate::kernels::ode::solve_adaptive","Native ode::solve_adaptive: this legacy driver supports rtol/atol/max_steps; h0/hmin must be zero and hmax must cover the whole span, with unsupported requests rejected. independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,1));

object_binding!(fixed_events,"num_analysis_ode_fixed_events","crate::kernels::ode::solve_fixed_events","Native ode::solve_fixed_events; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,2));

object_binding!(bdf_events,"num_analysis_ode_bdf_events","crate::kernels::ode_stiff::bdf_events","Native ode_stiff::bdf_events; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,3));

object_binding!(bdf_jac,"num_analysis_ode_bdf_jac","crate::kernels::ode_stiff::bdf_jac","Native ode_stiff::bdf_jac; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,4));

object_binding!(bdf_ex,"num_analysis_ode_bdf_ex","crate::kernels::ode_stiff::bdf_ex","Native ode_stiff::bdf_ex; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,5));

object_binding!(bdf_with,"num_analysis_ode_bdf_with","crate::kernels::ode_stiff::bdf_with","Native ode_stiff::bdf_with; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,6));

object_binding!(bdf_parallel,"num_analysis_ode_bdf_parallel","crate::kernels::ode_stiff::bdf_with_par","Native ode_stiff::bdf_with_par; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,7));

object_binding!(radau_events,"num_analysis_ode_radau_events","crate::kernels::ode_stiff::radau5_events","Native ode_stiff::radau5_events; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,8));

object_binding!(radau_jac,"num_analysis_ode_radau_jac","crate::kernels::ode_stiff::radau5_jac","Native ode_stiff::radau5_jac; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,9));

object_binding!(radau_with,"num_analysis_ode_radau_with","crate::kernels::ode_stiff::radau_with","Native ode_stiff::radau_with; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,10));

object_binding!(radau5_fixed,"num_analysis_ode_radau5_fixed","crate::kernels::ode_stiff::radau5_fixed","Native ode_stiff::radau5_fixed; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,11));

object_binding!(radau_fixed,"num_analysis_ode_radau_fixed","crate::kernels::ode_stiff::radau_fixed","Native ode_stiff::radau_fixed; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,12));

object_binding!(sdirk,"num_analysis_ode_sdirk","crate::kernels::ode_sdirk::sdirk3","Native ode_sdirk::sdirk3; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,13));

object_binding!(sdirk_fd,"num_analysis_ode_sdirk_fd","crate::kernels::ode_sdirk::sdirk3_fd","Native ode_sdirk::sdirk3_fd; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,14));

object_binding!(sdirk_events,"num_analysis_ode_sdirk_events","crate::kernels::ode_sdirk::sdirk3_events","Native ode_sdirk::sdirk3_events; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,15));

object_binding!(sdirk_fixed,"num_analysis_ode_sdirk_fixed","crate::kernels::ode_sdirk::sdirk3_fixed","Native ode_sdirk::sdirk3_fixed; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,16));

object_binding!(sdirk_fixed_events,"num_analysis_ode_sdirk_fixed_events","crate::kernels::ode_sdirk::sdirk3_fixed_events","Native ode_sdirk::sdirk3_fixed_events; independently measured initial components and physical times are normalized by the retained RHS tape references. Config rtol/atol are explicit normalized-state tolerances; h0/hmin/hmax are physical times and max_steps is mandatory. Controls contain steps/order/stages as documented by the native function; events contain program/terminal/direction. Real native evaluation and retained sample budgets are enforced.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,17));

object_binding!(radau_legacy,"num_analysis_ode_radau","crate::kernels::ode_stiff::radau","Actual legacy native variable-stage Radau entry, including its true collocation dense output and event handling. Physical frames and explicit bounds match the modern Radau adapter.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,18));
object_binding!(rosenbrock_events,"num_analysis_ode_rosenbrock23_events","crate::kernels::ode_adaptive::rosenbrock23_events","Actual native Rosenbrock2/3 with native event localization and dense output in independently scaled physical state coordinates.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,controls:Vec<usize>,events:Many<Input>,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,controls,events,maxeval|run(ctx,rhs,initial,start,end,config,controls,events,maxeval,19));
object_binding!(view,"num_analysis_ode_view","crate::kernels::ode_adaptive::OdeSolution::view","Actual native named solution view; native component naming/padding is preserved while samples and events retain physical dimensions.",false,(solution:Handle<Solution>,components:Vec<String>)->Dynamic,|ctx,solution,components|solution.with(ctx,|s|{if components.len()>128||components.iter().any(|s|s.len()>128){return Err(api::error("component name budget exceeded"));}let SolutionData::Dense(x)=&s.data else{return Err(api::error("OdeSolution.view requires a native adaptive dense solution"));};let names=components.iter().map(String::as_str).collect::<Vec<_>>();let view=x.view(&names);let states=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![view.y.len()],elem:api::dynamic_type(),data:view.y.iter().map(|y|s.frame.restore(y)).collect::<Result<_,_>>()?});Ok(Dynamic(api::record(vec![("times",restored_times(s,&view.t)?),("states",states),("components",Output::encode(view.flatten_order)?),("events",restored_hits(s,x.events.iter().map(|h|(h.index,h.t,h.y.as_slice())))?),("nfev",Output::encode(view.nfev)?),("rejected",Output::encode(view.rejected)?)])))}));
object_binding!(last,"num_analysis_ode_last","crate::kernels::ode::Trajectory::last","Actual native trajectory last; only valid on a native discrete trajectory.",false,(solution:Handle<Solution>)->Dynamic,|ctx,solution|solution.with(ctx,|s|match &s.data{SolutionData::Trajectory(t)|SolutionData::FixedEvents(t,_)=>Ok(Dynamic(s.frame.restore(t.last())?)),_=>Err(api::error("Trajectory.last requires a discrete native trajectory"))}));
object_binding!(eval,"num_analysis_ode_eval","crate::kernels::ode_adaptive::OdeSolution::eval","Actual native dense output; no invented interpolant for a discrete trajectory. Outside the solved interval returns None.",false,(solution:Handle<Solution>,time:Input)->Maybe<Dynamic>,|ctx,solution,time|solution.with(ctx,|s|{let time=s.frame.time(&time)?;let value=match &s.data{SolutionData::Dense(x)=>x.eval(time),SolutionData::Sdirk(x)=>x.eval(time),SolutionData::Geometric(x)=>crate::kernels::ode_geom::evaluate_history(&x.t,&x.y,&x.rates,time).ok(),SolutionData::Trajectory(_)|SolutionData::FixedEvents(_,_)=>return Err(api::error("discrete trajectories have no native dense output"))};Ok(Maybe(value.map(|v|s.frame.restore(&v).map(Dynamic)).transpose()?))}));
object_binding!(sdirk_eval,"num_analysis_sdirk_eval","crate::kernels::ode_sdirk::SdirkSolution::eval","Actual native SDIRK Hermite dense output in physical state coordinates.",false,(solution:Handle<Solution>,time:Input)->Maybe<Dynamic>,|ctx,solution,time|solution.with(ctx,|s|{let time=s.frame.time(&time)?;let SolutionData::Sdirk(x)=&s.data else{return Err(api::error("SdirkSolution.eval requires native SDIRK solution"));};Ok(Maybe(x.eval(time).map(|v|s.frame.restore(&v).map(Dynamic)).transpose()?))}));
fn restored_times(
    s: &Solution,
    t: &[f64],
) -> Result<Value, OpError> {
    Ok(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![t.len()],
        elem: crate::kernel::signature::TypeExpr::Dim(
            Box::new(crate::kernel::signature::TypeExpr::Float),
            s.frame.time_dimension.clone(),
        ),
        data: t
            .iter()
            .map(|v| api::scale_product(*v, &[s.frame.time_scale], &[]).map(Value::Float))
            .collect::<Result<_, _>>()?,
    }))
}
fn restored_hits<'a>(
    s: &Solution,
    hits: impl Iterator<Item = (usize, f64, &'a [f64])>,
) -> Result<Value, OpError> {
    let data = hits
        .map(|(index, t, y)| {
            Ok(api::record(vec![
                ("index", Output::encode(index)?),
                (
                    "time",
                    api::quantity(
                        api::scale_product(t, &[s.frame.time_scale], &[])?,
                        s.frame.time_dimension.clone(),
                    )?,
                ),
                ("state", s.frame.restore(y)?),
            ]))
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![data.len()],
        elem: api::dynamic_type(),
        data,
    }))
}
object_binding!(data,"num_analysis_ode_data","","Physical accepted samples and native event/order statistics. Fixed endpoint-only solvers return only their actual initial/final states.",false,(solution:Handle<Solution>)->Dynamic,|ctx,solution|solution.with(ctx,|s|{let(t,y)=s.nodes();let mut fields=vec![("times",restored_times(s,t)?),("states",Value::Tensor(crate::kernel::value::TensorValue{shape:vec![y.len()],elem:api::dynamic_type(),data:y.iter().map(|row|s.frame.restore(row)).collect::<Result<Vec<_>,_>>()?})),("evaluations",Output::encode(s.evaluations)?)];match &s.data{SolutionData::Dense(x)=>{fields.push(("orders",Output::encode(x.orders.clone())?));fields.push(("rejected",Output::encode(x.rejected)?));fields.push(("events",restored_hits(s,x.events.iter().map(|e|(e.index,e.t,e.y.as_slice())))?));},SolutionData::FixedEvents(_,hits)=>{fields.push(("events",restored_hits(s,hits.iter().map(|e|(e.index,e.t,e.y.as_slice())))?));},SolutionData::Geometric(x)=>{fields.push(("rejected",Output::encode(x.rejected)?));fields.push(("events",restored_hits(s,x.hits.iter().map(|e|(e.index,e.t,e.y.as_slice())))?));},SolutionData::Sdirk(x)=>{fields.push(("rejected",Output::encode(x.rejected)?));fields.push(("newton_iterations",Output::encode(x.newton_iterations)?));fields.push(("jacobian_evaluations",Output::encode(x.njev)?));},SolutionData::Trajectory(_)=>{}}Ok(Dynamic(api::record(fields)))}));

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    geometric::install(i)?;
    extended::install(i)?;
    dae::install(i)?;
    dde::install(i)?;
    hamiltonian::install(i)?;
    boundary::install(i)?;
    api::install(
        i,
        &[
            &radau_legacy::BINDING,
            &rosenbrock_events::BINDING,
            &view::BINDING,
            &fixed::BINDING,
            &adaptive::BINDING,
            &fixed_events::BINDING,
            &bdf_events::BINDING,
            &bdf_jac::BINDING,
            &bdf_ex::BINDING,
            &bdf_with::BINDING,
            &bdf_parallel::BINDING,
            &radau_events::BINDING,
            &radau_jac::BINDING,
            &radau_with::BINDING,
            &radau5_fixed::BINDING,
            &radau_fixed::BINDING,
            &sdirk::BINDING,
            &sdirk_fd::BINDING,
            &sdirk_events::BINDING,
            &sdirk_fixed::BINDING,
            &sdirk_fixed_events::BINDING,
            &last::BINDING,
            &eval::BINDING,
            &sdirk_eval::BINDING,
            &data::BINDING,
        ],
        "numerical_vector",
    )
}
