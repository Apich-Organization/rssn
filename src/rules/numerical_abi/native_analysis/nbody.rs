//! Actual chain-regularized and Hermite native N-body systems with SI frames.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use crate::sim::chain as c;
use crate::sim::hermite as h;
#[path = "symplectic.rs"]
mod symplectic;
#[derive(Clone)]
struct Frame {
    length: f64,
    mass: f64,
    time: f64,
}
fn l() -> DimExpr {
    DimExpr::base("Length")
}
fn m() -> DimExpr {
    DimExpr::base("Mass")
}
fn t() -> DimExpr {
    DimExpr::base("Time")
}
fn inv(d: DimExpr) -> DimExpr {
    DimExpr::Inverse(Box::new(d))
}
fn power(
    d: DimExpr,
    n: i64,
) -> DimExpr {
    d.pow(Exp::Rat(Rational::int(n)))
}
fn velocity() -> DimExpr {
    DimExpr::Mul(vec![l(), inv(t())])
}
fn energy() -> DimExpr {
    DimExpr::Mul(vec![m(), power(l(), 2), power(t(), -2)])
}
fn angular() -> DimExpr {
    DimExpr::Mul(vec![m(), power(l(), 2), inv(t())])
}
fn gdim() -> DimExpr {
    DimExpr::Mul(vec![power(l(), 3), inv(m()), power(t(), -2)])
}
impl Frame {
    fn new(
        length: &Input,
        mass: &Input,
        time: &Input,
    ) -> Result<Self, OpError> {
        let (length, ld) = api::reference(&length.0)?;
        let (mass, md) = api::reference(&mass.0)?;
        let (time, td) = api::reference(&time.0)?;
        api::same_dimension(&ld, &l())?;
        api::same_dimension(&md, &m())?;
        api::same_dimension(&td, &t())?;
        Ok(Self { length, mass, time })
    }

    fn vel(&self) -> Result<f64, OpError> {
        api::scale_product(self.length, &[], &[self.time])
    }

    fn energy(&self) -> Result<f64, OpError> {
        api::scale_product(
            1.0,
            &[self.mass, self.length, self.length],
            &[self.time, self.time],
        )
    }

    fn angular(&self) -> Result<f64, OpError> {
        api::scale_product(1.0, &[self.mass, self.length, self.length], &[self.time])
    }

    fn g(
        &self,
        v: &Input,
    ) -> Result<f64, OpError> {
        number(
            v,
            api::scale_product(
                1.0,
                &[self.length, self.length, self.length],
                &[self.mass, self.time, self.time],
            )?,
            gdim(),
        )
    }
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn number(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    api::scale_product(
        crate::rules::numerical_abi::numeric(&v.0, 0)?,
        &[],
        &[scale],
    )
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn vector(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<Vec<f64>, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    let v = <Vec<f64> as Arg>::decode(&v.0, 0)?;
    if v.len() > 128 {
        return Err(api::error("N-body count budget exceeds128"));
    }
    v.into_iter()
        .map(|v| api::scale_product(v, &[], &[scale]))
        .collect()
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn vectors(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<Vec<[f64; 3]>, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    let mut v = <Vec<[f64; 3]> as Arg>::decode(&v.0, 0)?;
    if v.len() > 128 {
        return Err(api::error("N-body count budget exceeds128"));
    }
    for row in &mut v {
        for v in row {
            *v = api::scale_product(*v, &[], &[scale])?;
        }
    }
    Ok(v)
}
fn out(
    v: &[[f64; 3]],
    scale: f64,
    dim: DimExpr,
) -> Result<Value, OpError> {
    Ok(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![v.len(), 3],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dim),
        data: v
            .iter()
            .flatten()
            .map(|v| api::scale_product(*v, &[scale], &[]).map(Value::Float))
            .collect::<Result<_, _>>()?,
    }))
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn cartesian(
    frame: &Frame,
    pair: (Vec<[f64; 3]>, Vec<[f64; 3]>),
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("positions", out(&pair.0, frame.length, l())?),
        ("velocities", out(&pair.1, frame.vel()?, velocity())?),
    ])))
}
struct Chain {
    data: c::ChainSystem,
    frame: Frame,
}
impl NativeType for Chain {
    const KIND: &'static str = "rssn/native/analysis/nbody/chain/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.mass.len() > 128
            || !self.data.t.is_finite()
            || self
                .data
                .xc
                .iter()
                .chain(&self.data.vc)
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(api::error("invalid native chain size/state"));
        }
        Ok(())
    }
}
struct Hermite {
    data: h::HermiteSystem,
    frame: Frame,
    min_step: f64,
}
impl NativeType for Hermite {
    const KIND: &'static str = "rssn/native/analysis/nbody/hermite/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.mass.len() > 128
            || self
                .data
                .pos
                .iter()
                .chain(&self.data.vel)
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(api::error("invalid native Hermite size/state"));
        }
        Ok(())
    }
}
object_binding!(chain_new,"num_analysis_nbody_chain_new","crate::sim::chain::ChainSystem::new","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(g:Input,masses:Input,positions:Input,velocities:Input,length_ref:Input,mass_ref:Input,time_ref:Input)->NewObject<Chain>,|_ctx,g,masses,positions,velocities,length_ref,mass_ref,time_ref|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let g=frame.g(&g)?;let mass=vector(&masses,frame.mass,m())?;let pos=vectors(&positions,frame.length,l())?;let vel=vectors(&velocities,frame.vel()?,velocity())?;let data=c::ChainSystem::new(g,&mass,&pos,&vel).map_err(api::error)?;Ok(NewObject::new(Chain{data,frame},vec![l(),m(),t()]))});
object_binding!(chain_cartesian,"num_analysis_nbody_chain_cartesian","crate::sim::chain::ChainSystem::cartesian","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->Dynamic,|ctx,chain|chain.with(ctx,|c|cartesian(&c.frame,c.data.cartesian())));
object_binding!(chain_energy,"num_analysis_nbody_chain_energy","crate::sim::chain::ChainSystem::energy","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->Dynamic,|ctx,chain|chain.with(ctx,|c|Ok(Dynamic(api::quantity(api::scale_product(c.data.energy(),&[c.frame.energy()?],&[])?,energy())?))));
object_binding!(chain_kinetic,"num_analysis_nbody_chain_kinetic","crate::sim::chain::ChainSystem::kinetic","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->Dynamic,|ctx,chain|chain.with(ctx,|c|Ok(Dynamic(api::quantity(api::scale_product(c.data.kinetic(),&[c.frame.energy()?],&[])?,energy())?))));
object_binding!(chain_potential,"num_analysis_nbody_chain_potential","crate::sim::chain::ChainSystem::potential","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->Dynamic,|ctx,chain|chain.with(ctx,|c|Ok(Dynamic(api::quantity(api::scale_product(c.data.potential(),&[c.frame.energy()?],&[])?,energy())?))));
object_binding!(chain_angular_momentum,"num_analysis_nbody_chain_angular_momentum","crate::sim::chain::ChainSystem::angular_momentum","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->Dynamic,|ctx,chain|chain.with(ctx,|c|Ok(Dynamic(out(&[c.data.angular_momentum()],c.frame.angular()?,angular())?))));
object_binding!(chain_rechain,"num_analysis_nbody_chain_rechain","crate::sim::chain::ChainSystem::rechain","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>)->(),|ctx,chain|chain.with_mut(ctx,|c|{c.data.rechain();Ok(())}));
object_binding!(chain_leapfrog,"num_analysis_nbody_chain_leapfrog","crate::sim::chain::ChainSystem::leapfrog","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>,fictitious_step:Input)->(),|ctx,chain,fictitious_step|chain.with_mut(ctx,|c|{let step=number(&fictitious_step,api::scale_product(c.frame.energy()?,&[c.frame.time],&[])?,DimExpr::Mul(vec![energy(),t()]))?;c.data.leapfrog(step);c.validate()}));
object_binding!(chain_gbs_step,"num_analysis_nbody_chain_gbs_step","crate::sim::chain::ChainSystem::gbs_step","Actual native GBS step. Its returned error is explicitly in the native normalized state metric, while fictitious time has physical energy*time units.",false,(chain:Handle<Chain>,fictitious_step:Input,kmax:usize)->f64,|ctx,chain,fictitious_step,kmax|chain.with_mut(ctx,|c|{if !(1..=8).contains(&kmax){return Err(api::error("chain extrapolation cap must be1..8"));}let ds=number(&fictitious_step,api::scale_product(c.frame.energy()?,&[c.frame.time],&[])?,DimExpr::Mul(vec![energy(),t()]))?;let error=c.data.gbs_step(ds,kmax);c.validate()?;Ok(error)}));
object_binding!(chain_evolve_until,"num_analysis_nbody_chain_evolve_until","crate::sim::chain::ChainSystem::evolve_until","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(chain:Handle<Chain>,end:Input,fictitious_step:Input,kmax:usize,tolerance:f64,max_native_work:usize)->usize,|ctx,chain,end,fictitious_step,kmax,tolerance,max_native_work|chain.with_mut(ctx,|c|{let n=c.data.mass.len();if !(1..=8).contains(&kmax)||max_native_work>1000000000||5000000usize.checked_mul(n*n).and_then(|w|w.checked_mul(kmax*kmax)).is_none_or(|w|w>max_native_work)||tolerance<=0.0||!tolerance.is_finite(){return Err(api::error("native chain worst-case5000000attempts exceed explicit work bound"));}let end=number(&end,c.frame.time,t())?;let ds=number(&fictitious_step,api::scale_product(c.frame.energy()?,&[c.frame.time],&[])?,DimExpr::Mul(vec![energy(),t()]))?;if ds<=0.0{return Err(api::error("chain fictitious step must be positive"));}let steps=c.data.evolve_until(end,ds,kmax,tolerance).map_err(api::error)?;c.validate()?;Ok(steps)}));
object_binding!(hermite_new,"num_analysis_nbody_hermite_new","crate::sim::hermite::HermiteSystem::new","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(g:Input,masses:Input,positions:Input,velocities:Input,length_ref:Input,mass_ref:Input,time_ref:Input,eps:Input,start:Input,order:usize,policy:String,hmin:Input,hmax:Input,eta:f64)->NewObject<Hermite>,|_ctx,g,masses,positions,velocities,length_ref,mass_ref,time_ref,eps,start,order,policy,hmin,hmax,eta|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let g=frame.g(&g)?;let mass=vector(&masses,frame.mass,m())?;let pos=vectors(&positions,frame.length,l())?;let vel=vectors(&velocities,frame.vel()?,velocity())?;let eps=number(&eps,frame.length,l())?;let t0=number(&start,frame.time,t())?;let min_step=number(&hmin,frame.time,t())?;let max_step=number(&hmax,frame.time,t())?;if min_step<=0.0||max_step<min_step||eta<=0.0||!eta.is_finite(){return Err(api::error("invalid Hermite step/eta bounds"));}let order=match order{4=>h::HermiteOrder::Fourth,6=>h::HermiteOrder::Sixth,_=>return Err(api::error("Hermite order must be4or6"))};let policy=match policy.as_str(){"Shared"=>{if !crate::kernels::exact_eq(min_step,max_step){return Err(api::error("shared Hermite step requires equal min/max"));}h::StepPolicy::Shared(max_step)},"Block"=>h::StepPolicy::Block{eta,dt_max:max_step,dt_min:min_step},_=>return Err(api::error("Hermite policy must beShared/Block"))};let data=h::HermiteSystem::new(mass,pos,vel,&h::HermiteConfig{g,eps,order,policy,t0}).map_err(api::error)?;Ok(NewObject::new(Hermite{data,frame,min_step},vec![l(),m(),t()]))});
object_binding!(hermite_step,"num_analysis_nbody_hermite_step","crate::sim::hermite::HermiteSystem::step","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>)->(),|ctx,system|system.with_mut(ctx,|s|{s.data.step().map_err(api::error)?;s.validate()}));
object_binding!(hermite_advance,"num_analysis_nbody_hermite_advance","crate::sim::hermite::HermiteSystem::advance_to","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>,end:Input,max_work:usize)->(),|ctx,system,end,max_work|system.with_mut(ctx,|s|{let end=number(&end,s.frame.time,t())?;let n=s.data.mass.len();let first=s.data.t_last.iter().copied().fold(s.data.time,f64::min);let attempted=((end-first).max(0.0)/s.min_step).ceil()+1.0;if max_work>100000000||attempted*n.pow(3)as f64>max_work as f64{return Err(api::error("Hermite conservative block/force work bound exceeded"));}s.data.advance_to(end).map_err(api::error)?;s.validate()}));
object_binding!(hermite_next_time,"num_analysis_nbody_hermite_next_time","crate::sim::hermite::HermiteSystem::next_time","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>)->Dynamic,|ctx,system|system.with(ctx,|s|Ok(Dynamic(api::quantity(api::scale_product(s.data.next_time(),&[s.frame.time],&[])?,t())?))));
object_binding!(hermite_state_at,"num_analysis_nbody_hermite_state_at","crate::sim::hermite::HermiteSystem::state_at","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>,time:Input)->Dynamic,|ctx,system,time|system.with(ctx,|s|cartesian(&s.frame,s.data.state_at(number(&time,s.frame.time,t())?).map_err(api::error)?)));
object_binding!(hermite_energy,"num_analysis_nbody_hermite_energy","crate::sim::hermite::HermiteSystem::total_energy","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>)->Dynamic,|ctx,system|system.with(ctx,|s|Ok(Dynamic(api::quantity(api::scale_product(s.data.total_energy().map_err(api::error)?,&[s.frame.energy()?],&[])?,energy())?))));
object_binding!(hermite_enable_post_newtonian,"num_analysis_nbody_hermite_enable_post_newtonian","crate::sim::hermite::HermiteSystem::enable_post_newtonian","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>,light_speed:Input)->(),|ctx,system,light_speed|system.with_mut(ctx,|s|{let c=number(&light_speed,s.frame.vel()?,velocity())?;s.data.enable_post_newtonian(c).map_err(api::error)}));
object_binding!(hermite_add_binary,"num_analysis_nbody_hermite_add_binary","crate::sim::hermite::HermiteSystem::add_binary","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(system:Handle<Hermite>,a:usize,b:usize,ks_eta:f64)->(),|ctx,system,a,b,ks_eta|system.with_mut(ctx,|s|{if !ks_eta.is_finite()||!(0.001..=1.0).contains(&ks_eta){return Err(api::error("KS accuracy parameter must be0.001..1"));}s.data.add_binary(a,b,ks_eta).map_err(api::error)}));
object_binding!(eih_accelerations,"num_analysis_nbody_eih_accelerations","crate::sim::hermite::eih_accelerations","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(masses:Input,positions:Input,velocities:Input,newton_acceleration:Input,g:Input,light_speed:Input,length_ref:Input,mass_ref:Input,time_ref:Input)->Dynamic,|_ctx,masses,positions,velocities,newton_acceleration,g,light_speed,length_ref,mass_ref,time_ref|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let mass=vector(&masses,frame.mass,m())?;let pos=vectors(&positions,frame.length,l())?;let vel=vectors(&velocities,frame.vel()?,velocity())?;if mass.len()!=pos.len()||mass.len()!=vel.len()||mass.iter().any(|m|*m<=0.0){return Err(api::error("EIH body shape/positive mass mismatch"));}let g=frame.g(&g)?;let light_speed=number(&light_speed,frame.vel()?,velocity())?;if light_speed<=0.0{return Err(api::error("EIH speed must bepositive"));}let dim=DimExpr::Mul(vec![l(),power(t(),-2)]);let scale=api::scale_product(frame.length,&[],&[frame.time,frame.time])?;let acc=vectors(&newton_acceleration,scale,dim.clone())?;if acc.len()!=mass.len(){return Err(api::error("EIH acceleration shape mismatch"));}Ok(Dynamic(out(&h::eih_accelerations(&mass,&pos,&vel,&acc,g,light_speed),scale,dim)?))});
object_binding!(eih_energy,"num_analysis_nbody_eih_energy","crate::sim::hermite::eih_energy","Actual native N-body operation. Positive SI mass/length/time references preserve physical input/output domains. Particle, extrapolation and attempted-work bounds are checked before native computation.",false,(masses:Input,positions:Input,velocities:Input,g:Input,light_speed:Input,length_ref:Input,mass_ref:Input,time_ref:Input)->Dynamic,|_ctx,masses,positions,velocities,g,light_speed,length_ref,mass_ref,time_ref|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let mass=vector(&masses,frame.mass,m())?;let pos=vectors(&positions,frame.length,l())?;let vel=vectors(&velocities,frame.vel()?,velocity())?;if mass.len()!=pos.len()||mass.len()!=vel.len()||mass.iter().any(|m|*m<=0.0){return Err(api::error("EIH body shape/positive mass mismatch"));}let g=frame.g(&g)?;let light_speed=number(&light_speed,frame.vel()?,velocity())?;if light_speed<=0.0{return Err(api::error("EIH speed must bepositive"));}Ok(Dynamic(api::quantity(api::scale_product(h::eih_energy(&mass,&pos,&vel,g,light_speed),&[frame.energy()?],&[])?,energy())?))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    symplectic::install(i)?;
    api::install(
        i,
        &[
            &chain_new::BINDING,
            &chain_cartesian::BINDING,
            &chain_energy::BINDING,
            &chain_kinetic::BINDING,
            &chain_potential::BINDING,
            &chain_angular_momentum::BINDING,
            &chain_rechain::BINDING,
            &chain_leapfrog::BINDING,
            &chain_gbs_step::BINDING,
            &chain_evolve_until::BINDING,
            &hermite_new::BINDING,
            &hermite_step::BINDING,
            &hermite_advance::BINDING,
            &hermite_next_time::BINDING,
            &hermite_state_at::BINDING,
            &hermite_energy::BINDING,
            &hermite_enable_post_newtonian::BINDING,
            &hermite_add_binary::BINDING,
            &eih_accelerations::BINDING,
            &eih_energy::BINDING,
        ],
        "numerical_vector",
    )
}
