//! Native coordinate-splitting and extended-Hamiltonian TTL drivers.
use super::Budget;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::TypeExpr;
use super::Value;
use super::api;
use super::c;
use super::energy;
use super::l;
use super::m;
use super::number;
use super::object_binding;
use super::order;
use super::parameters;
use super::s;
use super::scheme;
use super::t;
use super::tape_snapshot;
use super::vector;
use super::velocity;
use crate::kernels::numeric_common::trajectory_sample_count;
use crate::rules::numerical_abi::native_analysis_physics_ode::Acceleration;
use crate::rules::numerical_abi::native_object_api::Many;
use std::num::NonZeroUsize;
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn split(
    ctx: &CallCtx<'_>,
    acceleration: Handle<Acceleration>,
    initial: Input,
    h: Input,
    steps: usize,
    name: String,
    stages: Vec<[f64; 2]>,
    maxeval: usize,
    custom: bool,
) -> Result<Dynamic, OpError> {
    let a = acceleration.with(ctx, |a| Ok(a.clone()))?;
    let frame = a.frame()?;
    let state = frame.state(&initial)?;
    let n = state.len() / 2;
    let (mut q, mut p) = (state[..n].to_vec(), state[n..].to_vec());
    let h = frame.time(&h)?;
    if steps == 0 || steps > 100000 || stages.len() > 128 {
        return Err(api::error("native split stage/step budget exceeded"));
    }
    let budget = Budget::new(maxeval)?;
    let mut acceleration = |q: &[f64], out: &mut [f64]| a.evaluate(&budget, q, out);
    if custom {
        if stages.is_empty() {
            return Err(api::error("custom splitting needs at least one stage"));
        }
        let stages = stages.into_iter().map(|[a, b]| (a, b)).collect::<Vec<_>>();
        s::step(&stages, &mut q, &mut p, h, &mut acceleration);
    } else {
        s::integrate(scheme(&name)?, &mut q, &mut p, h, steps, &mut acceleration);
    }
    budget.finish()?;
    q.extend(p);
    Ok(Dynamic(frame.restore(&q)?))
}
object_binding!(split_step,"num_analysis_symplectic_split_step","crate::sim::symplectic::step","Actual native custom-stage drift/kick step. State contains physical positions then velocities and stage coefficients are explicitly normalized; a pure-position acceleration tape preserves coordinate/time² units.",false,(acceleration:Handle<Acceleration>,initial:Input,h:Input,stages:Vec<[f64;2]>,maxeval:usize)->Dynamic,|ctx,acceleration,initial,h,stages,maxeval|split(ctx,acceleration,initial,h,1,"Verlet".into(),stages,maxeval,true));
object_binding!(split_integrate,"num_analysis_symplectic_integrate","crate::sim::symplectic::integrate","Native Hamiltonian splitting with six schemes, physical coordinate, velocity and time references, and a shared evaluation budget.",false,(acceleration:Handle<Acceleration>,initial:Input,h:Input,steps:usize,name:String,maxeval:usize)->Dynamic,|ctx,acceleration,initial,h,steps,name,maxeval|split(ctx,acceleration,initial,h,steps,name,Vec::new(),maxeval,false));

struct SplitSamples {
    steps: Vec<usize>,
    states: Many<Dynamic>,
}

impl api::NativeOutput for SplitSamples {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("steps".into(), api::tensor_type(1, TypeExpr::Int)),
            ("states".into(), <Many<Dynamic> as api::NativeOutput>::ty()),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("steps", Output::encode(self.steps)?),
            ("states", api::NativeOutput::encode(self.states, ctx)?),
        ]))
    }
}

object_binding!(
    split_integrate_sampled,
    "num_analysis_symplectic_integrate_sampled",
    "crate::sim::symplectic::integrate_sampled",
    "Integrates every step with a separable symplectic scheme. Retains step zero, each positive save_every multiple, and the final step. Each state preserves the physical positions followed by velocities. Step indices are exact and the signed step size uses the acceleration tape's time dimension.",
    false,
    (
        acceleration: Handle<Acceleration>,
        initial: Input,
        h: Input,
        steps: usize,
        name: String,
        save_every: usize,
        maxeval: usize,
    ) -> SplitSamples,
    |ctx, acceleration, initial, h, steps, name, save_every, maxeval| {
        if steps > 100_000 {
            return Err(api::error("native split step budget exceeded"));
        }
        let cadence = NonZeroUsize::new(save_every)
            .ok_or_else(|| api::error("native split sample cadence must be positive"))?;
        let a = acceleration.with(ctx, |a| Ok(a.clone()))?;
        let frame = a.frame()?;
        let state = frame.state(&initial)?;
        let h = frame.time(&h)?;
        if trajectory_sample_count(steps, cadence)
            .and_then(|count| count.checked_mul(state.len()))
            .is_none_or(|n| n > 262_000)
        {
            return Err(api::error("native split sample storage exceeded"));
        }
        let (q, p) = state.split_at(state.len() / 2);
        let budget = Budget::new(maxeval)?;
        let samples = s::integrate_sampled(
            scheme(&name)?,
            q,
            p,
            h,
            steps,
            save_every,
            &mut |q, out| a.evaluate(&budget, q, out),
        );
        budget.finish()?;
        let samples = samples.map_err(api::error)?;
        let mut steps = Vec::with_capacity(samples.len());
        let states = samples.into_iter().map(|sample| {
            steps.push(sample.step);
            let mut state = sample.q;
            state.extend(sample.p);
            frame.restore(&state).map(Dynamic)
        }).collect::<Result<Vec<_>, _>>()?;
        Ok(SplitSamples { steps, states: Many(states) })
    }
);

struct Ttl {
    data: s::TtlState,
    frame: Frame,
    masses: Vec<f64>,
}
impl NativeType for Ttl {
    const KIND: &'static str = "rssn/native/analysis/symplectic/ttl/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.q.is_empty()
            || self.data.q.len() > 32
            || self.data.q.len() != self.data.p.len()
            || self.data.q.len() != self.masses.len()
            || self
                .data
                .q
                .iter()
                .chain(&self.data.p)
                .any(|v| !v.is_finite())
            || !self.data.t.is_finite()
        {
            return Err(api::error("invalid native TTL state"));
        }
        Ok(())
    }
}
object_binding!(ttl_new,"num_analysis_symplectic_ttl_new","","Build actual native TTL phase state with physical positions, momenta, masses and initial time. The fictitious integration increment has energy*time dimensions.",false,(positions:Input,momenta:Input,masses:Input,start:Input,length_ref:Input,mass_ref:Input,time_ref:Input)->NewObject<Ttl>,|_ctx,positions,momenta,masses,start,length_ref,mass_ref,time_ref|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let q=vector(&positions,frame.length,l())?;let p=vector(&momenta,api::scale_product(1.0,&[frame.mass,frame.length],&[frame.time])?,DimExpr::Mul(vec![m(),velocity()]))?;let masses=vector(&masses,frame.mass,m())?;if masses.iter().any(|m|*m<=0.0){return Err(api::error("TTL masses must be positive"));}let start_time=number(&start,frame.time,t())?;Ok(NewObject::new(Ttl{data:s::TtlState{q,p,t:start_time},frame,masses},vec![l(),m(),t()]))});
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn ttl(
    ctx: &CallCtx<'_>,
    state: Handle<Ttl>,
    binding: Input,
    fictitious_step: Input,
    degree: usize,
    potential: Input,
    params: Input,
    maxeval: usize,
    kind: usize,
) -> Result<f64, OpError> {
    let snapshot = tape_snapshot(ctx, &potential.0)?;
    let params = parameters(
        &params,
        &snapshot.parameter_scales,
        &snapshot.parameter_dimensions,
    )?;
    if snapshot.output_scales.len() != 1 {
        return Err(api::error("TTL potential tape must have one output"));
    }
    api::same_dimension(&snapshot.output_dimensions[0], &energy())?;
    for d in &snapshot.input_dimensions {
        api::same_dimension(d, &l())?;
    }
    state.with_mut(ctx, |s| {
        if snapshot.input_scales.len() != s.data.q.len() {
            return Err(api::error("TTL potential input/state mismatch"));
        }
        let b = number(&binding, s.frame.energy()?, energy())?;
        let step = number(
            &fictitious_step,
            api::scale_product(s.frame.energy()?, &[s.frame.time], &[])?,
            DimExpr::Mul(vec![energy(), t()]),
        )?;
        let kinetic = s::ttl_kinetic(&s.data, &s.masses);
        let drift_denominator = kinetic + b;
        if !drift_denominator.is_finite() || drift_denominator <= 0.0 {
            return Err(api::error(
                "TTL initial kinetic energy plus binding must be finite and positive (T+B>0)",
            ));
        }
        let budget = Budget::new(maxeval)?;
        let er = s.frame.energy()?;
        let mut potential = |q: &[f64], out: &mut [f64]| {
            out.fill(f64::NAN);
            let run = || {
                if !budget.charge() {
                    return Err(OpError::BudgetExceeded);
                }
                let input = q
                    .iter()
                    .zip(&snapshot.input_scales)
                    .map(|(v, r)| api::scale_product(*v, &[s.frame.length], &[*r]))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut value = [0.0];
                snapshot.tape.eval(&input, &params, &mut value);
                let value = api::scale_product(value[0], &[snapshot.output_scales[0]], &[er])?;
                if value >= 0.0 {
                    return Err(api::error("TTL requires a strictly negative potential"));
                }
                let gradient = snapshot.tape.gradient_reverse(0, &input, &params);
                let gradient = gradient
                    .iter()
                    .zip(&snapshot.input_scales)
                    .map(|(g, r)| {
                        api::scale_product(
                            *g,
                            &[snapshot.output_scales[0], s.frame.length],
                            &[er, *r],
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok::<_, OpError>((value, gradient))
            };
            match run() {
                | Ok((v, g)) => {
                    out.copy_from_slice(&g);
                    v
                },
                | Err(e) => {
                    budget.fail(e);
                    f64::NAN
                },
            }
        };
        let error = match kind {
            | 0 => {
                s::ttl_step(&mut s.data, &s.masses, b, step, &mut potential);
                0.0
            },
            | 1 => {
                order(degree)?;
                let result =
                    s::ttl_step_composed(&mut s.data, &s.masses, b, step, degree, &mut potential);
                budget.finish()?;
                result.map_err(api::error)?;
                0.0
            },
            | 2 => {
                if !(1..=8).contains(&degree) {
                    return Err(api::error("TTL GBS extrapolation cap must be1..8"));
                }
                c::ttl_gbs_step(&mut s.data, &s.masses, b, step, degree, &mut potential)
            },
            | _ => return Err(api::error("unknown native TTL operation")),
        };
        budget.finish()?;
        s.validate()?;
        Ok(error)
    })
}
object_binding!(ttl_step,"num_analysis_symplectic_ttl_step","crate::sim::symplectic::ttl_step","Actual native ttl_step with a true native tape reverse gradient of the physical negative potential. Binding has energy units and may be negative when initial T+B>0; fictitious step has energy*time units; the returned error is explicitly in the native normalized state metric.",false,(state:Handle<Ttl>,binding:Input,fictitious_step:Input,degree:usize,potential:Input,params:Input,maxeval:usize)->f64,|ctx,state,binding,fictitious_step,degree,potential,params,maxeval|ttl(ctx,state,binding,fictitious_step,degree,potential,params,maxeval,0));
object_binding!(ttl_step_composed,"num_analysis_symplectic_ttl_step_composed","crate::sim::symplectic::ttl_step_composed","Actual native ttl_step_composed with a true native tape reverse gradient of the physical negative potential. Binding has energy units and may be negative when initial T+B>0; fictitious step has energy*time units; the returned error is explicitly in the native normalized state metric.",false,(state:Handle<Ttl>,binding:Input,fictitious_step:Input,degree:usize,potential:Input,params:Input,maxeval:usize)->f64,|ctx,state,binding,fictitious_step,degree,potential,params,maxeval|ttl(ctx,state,binding,fictitious_step,degree,potential,params,maxeval,1));
object_binding!(ttl_gbs_step,"num_analysis_symplectic_ttl_gbs_step","crate::sim::chain::ttl_gbs_step","Actual native ttl_gbs_step with a true native tape reverse gradient of the physical negative potential. Binding has energy units and may be negative when initial T+B>0; fictitious step has energy*time units; the returned error is explicitly in the native normalized state metric.",false,(state:Handle<Ttl>,binding:Input,fictitious_step:Input,degree:usize,potential:Input,params:Input,maxeval:usize)->f64,|ctx,state,binding,fictitious_step,degree,potential,params,maxeval|ttl(ctx,state,binding,fictitious_step,degree,potential,params,maxeval,2));
object_binding!(ttl_data,"num_analysis_symplectic_ttl_data","","Physical retained native TTL time, positions and momenta.",false,(state:Handle<Ttl>)->Dynamic,|ctx,state|state.with(ctx,|s|{let frame=&s.frame;let vector=|v:&[f64],scale,dim|Ok::<_,OpError>(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![v.len()],elem:TypeExpr::Dim(Box::new(TypeExpr::Float),dim),data:v.iter().map(|v|api::scale_product(*v,&[scale],&[]).map(Value::Float)).collect::<Result<_,_>>()?}));Ok(Dynamic(api::record(vec![("positions",vector(&s.data.q,frame.length,l())?),("momenta",vector(&s.data.p,api::scale_product(1.0,&[frame.mass,frame.length],&[frame.time])?,DimExpr::Mul(vec![m(),velocity()]))?),("time",api::quantity(api::scale_product(s.data.t,&[frame.time],&[])?,t())?)])))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &split_step::BINDING,
            &split_integrate::BINDING,
            &split_integrate_sampled::BINDING,
            &ttl_new::BINDING,
            &ttl_step::BINDING,
            &ttl_step_composed::BINDING,
            &ttl_gbs_step::BINDING,
            &ttl_data::BINDING,
        ],
        "numerical_vector",
    )
}
