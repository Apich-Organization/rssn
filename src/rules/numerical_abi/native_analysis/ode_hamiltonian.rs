//! Native general Hamiltonian and separable symplectic trajectories.
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
use crate::kernel::signature::DimExpr;
use crate::kernels::ode_geom as g;
use crate::rules::numerical_abi::native_analysis_ode_common::VectorOracle;
use crate::rules::numerical_abi::native_analysis_ode_common::reference_frame;
use crate::rules::numerical_abi::native_analysis_physics_ode::Acceleration;
struct Phase {
    data: g::SymplecticSolution,
    q: OdeFrame,
    p: OdeFrame,
}
impl NativeType for Phase {
    const KIND: &'static str = "rssn/native/analysis/hamiltonian/solution/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.t.len() != self.data.q.len()
            || self.data.t.len() != self.data.p.len()
            || self
                .data
                .t
                .len()
                .checked_mul(self.q.state_scales.len() + self.p.state_scales.len())
                .is_none_or(|n| n > 47900)
        {
            return Err(api::error("Hamiltonian retained sample budget exceeded"));
        }
        Ok(())
    }
}
fn method(name: &str) -> Result<g::SymplecticMethod, OpError> {
    Ok(match name {
        | "StormerVerlet" => g::SymplecticMethod::StormerVerlet,
        | "Ruth3" => g::SymplecticMethod::Ruth3,
        | "Ruth4" => g::SymplecticMethod::Ruth4,
        | "Yoshida6" => g::SymplecticMethod::Yoshida6,
        | "Yoshida8" => g::SymplecticMethod::Yoshida8,
        | "ImplicitMidpoint" => g::SymplecticMethod::ImplicitMidpoint,
        | "GaussLegendre4" => g::SymplecticMethod::GaussLegendre4,
        | "GaussLegendre6" => g::SymplecticMethod::GaussLegendre6,
        | _ => return Err(api::error("unknown native SymplecticMethod variant")),
    })
}
fn rate_references(f: &OdeFrame) -> Result<(Vec<f64>, Vec<DimExpr>), OpError> {
    Ok((
        f.state_scales
            .iter()
            .map(|s| api::scale_product(*s, &[], &[f.time_scale]))
            .collect::<Result<_, _>>()?,
        f.state_dimensions
            .iter()
            .map(|d| {
                DimExpr::Mul(vec![
                    d.clone(),
                    DimExpr::Inverse(Box::new(f.time_dimension.clone())),
                ])
            })
            .collect(),
    ))
}
object_binding!(hamiltonian,"num_analysis_ode_integrate_hamiltonian","crate::kernels::ode_geom::integrate_hamiltonian","Actual general native Hamiltonian integrator with all8native methods. dq/dt and dp/dt tapes receive(time,q,p); q/p components carry independent physical references and rates preserve each component/time. Returned object retains actual native Hermite trajectory.",false,(time_reference:Input,q_references:Many<Input>,p_references:Many<Input>,q_rate:Input,q_params:Input,p_rate:Input,p_params:Input,initial_q:Input,initial_p:Input,start:Input,end:Input,steps:usize,name:String,maxeval:usize)->NewObject<Phase>,|ctx,time_reference,q_references,p_references,q_rate,q_params,p_rate,p_params,initial_q,initial_p,start,end,steps,name,maxeval|{let q=reference_frame(&time_reference,&q_references)?;let p=reference_frame(&time_reference,&p_references)?;if q.state_scales.len()!=p.state_scales.len()||steps==0||steps.checked_mul(2*q.state_scales.len()).is_none_or(|n|n>47900){return Err(api::error("Hamiltonian phase shape/sample budget mismatch"));}let mut input_scales=vec![q.time_scale];input_scales.extend(q.state_scales.clone());input_scales.extend(p.state_scales.clone());let mut input_dims=q.dimensions();input_dims.extend(p.state_dimensions.clone());let(qs,qd)=rate_references(&q)?;let(ps,pd)=rate_references(&p)?;let fq=VectorOracle::new(ctx,&q_rate,&q_params,&input_scales,&input_dims,&qs,&qd)?;let fp=VectorOracle::new(ctx,&p_rate,&p_params,&input_scales,&input_dims,&ps,&pd)?;let q0=q.state(&initial_q)?;let p0=p.state(&initial_p)?;let span=(q.time(&start)?,q.time(&end)?);let budget=Budget::new(maxeval)?;let mut qrate=|t,q:&[f64],p:&[f64]|fq.sample(&budget,&std::iter::once(t).chain(q.iter().copied()).chain(p.iter().copied()).collect::<Vec<_>>());let mut prate=|t,q:&[f64],p:&[f64]|fp.sample(&budget,&std::iter::once(t).chain(q.iter().copied()).chain(p.iter().copied()).collect::<Vec<_>>());let result=g::integrate_hamiltonian(&mut qrate,&mut prate,&q0,&p0,span,steps,method(&name)?);budget.finish()?;let mut dims=q.dimensions();dims.extend(p.state_dimensions.clone());Ok(NewObject::new(Phase{data:result.map_err(api::error)?,q,p},dims))});
object_binding!(symplectic,"num_analysis_ode_integrate_symplectic","crate::kernels::ode_geom::integrate_symplectic","Actual native separable symplectic integrator. Acceleration tape receives(time,position) and outputs position/time²; velocity references are derived from position/time. All8native methods and original native dense data are preserved.",false,(time_reference:Input,q_references:Many<Input>,acceleration:Input,params:Input,initial_q:Input,initial_velocity:Input,start:Input,end:Input,steps:usize,name:String,maxeval:usize)->NewObject<Phase>,|ctx,time_reference,q_references,acceleration,params,initial_q,initial_velocity,start,end,steps,name,maxeval|{let q=reference_frame(&time_reference,&q_references)?;if steps==0||steps.checked_mul(2*q.state_scales.len()).is_none_or(|n|n>47900){return Err(api::error("symplectic sample budget exceeded"));}let(ps,pd)=rate_references(&q)?;let p=OdeFrame{time_scale:q.time_scale,time_dimension:q.time_dimension.clone(),state_scales:ps,state_dimensions:pd};let(ascl,adim)=rate_references(&p)?;let mut input_scales=vec![q.time_scale];input_scales.extend(q.state_scales.clone());let oracle=VectorOracle::new(ctx,&acceleration,&params,&input_scales,&q.dimensions(),&ascl,&adim)?;let q0=q.state(&initial_q)?;let p0=p.state(&initial_velocity)?;let span=(q.time(&start)?,q.time(&end)?);let budget=Budget::new(maxeval)?;let mut acceleration=|t,q:&[f64]|oracle.sample(&budget,&std::iter::once(t).chain(q.iter().copied()).collect::<Vec<_>>());let result=g::integrate_symplectic(&mut acceleration,&q0,&p0,span,steps,method(&name)?);budget.finish()?;let mut dims=q.dimensions();dims.extend(p.state_dimensions.clone());Ok(NewObject::new(Phase{data:result.map_err(api::error)?,q,p},dims))});
object_binding!(eval,"num_analysis_ode_phase_eval","crate::kernels::ode_geom::SymplecticSolution::eval","Actual native symplectic trajectory Hermite evaluation. Independent position/momentum or velocity units are preserved; out-of-span queries remain None.",false,(solution:Handle<Phase>,time:Input)->Maybe<Dynamic>,|ctx,solution,time|solution.with(ctx,|s|Ok(Maybe(s.data.eval(s.q.time(&time)?).map(|(q,p)|Ok::<_,OpError>(Dynamic(api::record(vec![("q",s.q.restore(&q)?),("p",s.p.restore(&p)?)])))).transpose()?))));
object_binding!(energy_drift,"num_analysis_ode_max_energy_drift","crate::kernels::ode_geom::SymplecticSolution::max_energy_drift","Actual native maximum absolute energy drift on accepted phase nodes. Hamiltonian tape receives(q,p); positive energy reference must match every q_i*p_i/time domain. Native reverse/finite-difference substitution is not involved.",false,(solution:Handle<Phase>,hamiltonian:Input,params:Input,energy_reference:Input,maxeval:usize)->Dynamic,|ctx,solution,hamiltonian,params,energy_reference,maxeval|solution.with(ctx,|s|{let(energy,dim)=api::reference(&energy_reference.0)?;for(q,p)in s.q.state_dimensions.iter().zip(&s.p.state_dimensions){api::same_dimension(&dim,&DimExpr::Mul(vec![q.clone(),p.clone(),DimExpr::Inverse(Box::new(s.q.time_dimension.clone()))]))?;}let scales=[s.q.state_scales.clone(),s.p.state_scales.clone()].concat();let dims=[s.q.state_dimensions.clone(),s.p.state_dimensions.clone()].concat();let oracle=VectorOracle::new(ctx,&hamiltonian,&params,&scales,&dims,&[energy],std::slice::from_ref(&dim))?;let budget=Budget::new(maxeval)?;let value=s.data.max_energy_drift(|q,p|oracle.sample(&budget,&q.iter().chain(p).copied().collect::<Vec<_>>())[0]);budget.finish()?;Ok(Dynamic(api::quantity(api::scale_product(value,&[energy],&[])?,dim)?))}));
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn primitive(
    ctx: &CallCtx<'_>,
    acceleration: Handle<Acceleration>,
    initial: Input,
    dt: Input,
    steps: usize,
    maxeval: usize,
    fourth: bool,
) -> Result<Dynamic, OpError> {
    acceleration.with(ctx, |a| {
        let frame = a.frame()?;
        let initial = frame.state(&initial)?;
        let n = initial.len() / 2;
        let dt = frame.time(&dt)?;
        if steps == 0 || steps > 100000 {
            return Err(api::error(
                "primitive symplectic step budget must be1..100000",
            ));
        }
        let budget = Budget::new(maxeval)?;
        let f = |q: &[f64]| {
            let mut out = vec![0.0; n];
            a.evaluate(&budget, q, &mut out);
            out
        };
        let (mut q, p) = if fourth {
            a::yoshida4(f, &initial[..n], &initial[n..], dt, steps)
        } else {
            a::leapfrog(f, &initial[..n], &initial[n..], dt, steps)
        };
        budget.finish()?;
        q.extend(p);
        Ok(Dynamic(frame.restore(&q)?))
    })
}
object_binding!(leapfrog,"num_analysis_ode_leapfrog","crate::kernels::ode_adaptive::leapfrog","Actual native leapfrog over a pure-position acceleration tape with physical position/velocity/time references and explicit evaluation budget.",false,(acceleration:Handle<Acceleration>,initial:Input,dt:Input,steps:usize,maxeval:usize)->Dynamic,|ctx,acceleration,initial,dt,steps,maxeval|primitive(ctx,acceleration,initial,dt,steps,maxeval,false));
object_binding!(yoshida4,"num_analysis_ode_yoshida4","crate::kernels::ode_adaptive::yoshida4","Actual native fourth-order Yoshida over a pure-position acceleration tape, preserving physical phase/time references and original native evaluation errors.",false,(acceleration:Handle<Acceleration>,initial:Input,dt:Input,steps:usize,maxeval:usize)->Dynamic,|ctx,acceleration,initial,dt,steps,maxeval|primitive(ctx,acceleration,initial,dt,steps,maxeval,true));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &hamiltonian::BINDING,
            &symplectic::BINDING,
            &eval::BINDING,
            &energy_drift::BINDING,
            &leapfrog::BINDING,
            &yoshida4::BINDING,
        ],
        "numerical_vector",
    )
}
