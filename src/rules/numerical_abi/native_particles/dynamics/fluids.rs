#![allow(clippy::wildcard_imports)] // Native dissipative/SPH fluid states retain mixed physical dimensions.
use super::*;
use crate::sim::dpd;
use crate::sim::physics_mm as mm;
use crate::sim::sph;
struct Dpd {
    native: dpd::Dpd,
    frame: Frame,
}
struct Gas<const D: usize> {
    native: sph::Sph<D>,
    frame: Frame,
}
impl NativeType for Dpd {
    const KIND: &'static str = "rssn/native/sim/dpd/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        engine::bound(self.native.len(), 512)?;
        engine::finite(
            self.native
                .pos
                .iter()
                .chain(&self.native.vel)
                .flatten()
                .copied(),
        )?;
        engine::finite([self.native.t])
    }
}
impl<const D: usize> NativeType for Gas<D> {
    const KIND: &'static str = match D {
        | 1 => "rssn/native/sim/sph/d1/v1",
        | 2 => "rssn/native/sim/sph/d2/v1",
        | _ => "rssn/native/sim/sph/d3/v1",
    };
    const MAX_BYTES: usize = 524288;

    fn validate(&self) -> Result<(), OpError> {
        let n = self.native.pos.len();
        if n == 0
            || n > 128
            || self.native.vel.len() != n
            || self.native.mass.len() != n
            || self.native.u.len() != n
            || self.native.h.len() != n
            || self.native.rho.len() != n
        {
            return Err(api::error("SPH state shape/particle bound invalid"));
        }
        if self.native.mass.iter().any(|x| *x <= 0.0)
            || self.native.h.iter().any(|x| *x <= 0.0)
            || self.native.rho.iter().any(|x| *x <= 0.0)
            || self.native.u.iter().any(|x| *x < 0.0)
        {
            return Err(api::error(
                "SPH mass/density/smoothing/energy domain invalid",
            ));
        }
        engine::finite(
            self.native
                .pos
                .iter()
                .chain(&self.native.vel)
                .flatten()
                .copied(),
        )?;
        engine::finite(
            self.native
                .mass
                .iter()
                .chain(&self.native.u)
                .chain(&self.native.rho)
                .chain(&self.native.h)
                .copied(),
        )
    }
}
fn named(
    f: &Frame,
    v: &Value,
    name: &str,
    p: [i32; 4],
    temp: bool,
) -> Result<f64, OpError> {
    f.scalar(api::field(v, name)?, p, temp)
}
fn kernel(s: &str) -> Result<sph::KernelKind, OpError> {
    match s {
        | "CubicSpline" => Ok(sph::KernelKind::CubicSpline),
        | "WendlandC2" => Ok(sph::KernelKind::WendlandC2),
        | _ => Err(api::error("SPH kernel must be CubicSpline/WendlandC2")),
    }
}
fn boundary(value: &Value) -> Result<sph::SphBoundary, OpError> {
    api::normalized(value)?;
    let text = <String as Arg>::decode(value, 0)?;
    match text.as_str() {
        | "Open" => Ok(sph::SphBoundary::Open),
        | "Periodic" => Ok(sph::SphBoundary::Periodic),
        | "Reflecting" => Ok(sph::SphBoundary::Reflecting),
        | _ => Err(api::error("SPH boundary must be Open/Periodic/Reflecting")),
    }
}
fn fluid_steps(
    n: usize,
    steps: usize,
) -> Result<(), OpError> {
    if steps > 65536
        || n.checked_mul(n)
            .and_then(|n| n.checked_mul(steps))
            .is_none_or(|n| n > 50_000_000)
    {
        return Err(api::error("fluid particle-step work budget exceeded"));
    }
    Ok(())
}
object_binding!(dpd_new,"num_native_dpd_new","crate::sim::dpd::Dpd::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,n:usize,dimension:usize,parameters:Input,seed:u64)->NewObject<Dpd>,|ctx:&CallCtx<'_>,f:Handle<Frame>,n:usize,dimension:usize,parameters:Input,seed:u64|{let f=frame(ctx,&f)?;if !(2..=512).contains(&n)||!(2..=3).contains(&dimension){return Err(api::error("DPD needs2..512particles and dimension2/3"));}let len=engine::array::<3>(&f,api::field(&parameters.0,"box_length")?,LENGTH)?;let a=named(&f,&parameters.0,"repulsion",FORCE,false)?;let gamma=named(&f,&parameters.0,"friction",[1,0,-1,0],false)?;let kt=named(&f,&parameters.0,"temperature",ENERGY,true)?;let rc=named(&f,&parameters.0,"cutoff",LENGTH,false)?;let mass=named(&f,&parameters.0,"mass",MASS,false)?;if a<0.0||gamma<0.0||kt<0.0||rc<=0.0||mass<=0.0||len[..dimension].iter().any(|l|*l<3.0*rc){return Err(api::error("DPD physical parameters outside domain"));}let cells=len[..dimension].iter().try_fold(1usize,|n,l|{let side=(l/rc).floor();if side>1000000.0{return None;}n.checked_mul(side as usize).filter(|n|*n<=1_000_000)});if cells.is_none(){return Err(api::error("DPD cell grid allocation budget exceeded"));}let native=dpd::Dpd::new(n,dpd::DpdParams{dim:dimension,box_len:len,a,gamma,kt,rc,mass},seed).map_err(api::error)?;let dims=f.dims.to_vec();Ok(NewObject::new(Dpd{native,frame:f},dims))});
object_binding!(dpd_set_walls,"num_native_dpd_set_walls","crate::sim::dpd::Dpd::set_walls","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>,axis:usize,repulsion:Input,kind:String,friction:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>,axis:usize,repulsion:Input,kind:String,friction:f64|{fluid.with_mut(ctx,|s|{let a_wall=s.frame.scalar(&repulsion.0,FORCE,false)?;let kind=match kind.as_str(){"Specular"=>dpd::WallKind::Specular,"BounceBack"=>dpd::WallKind::BounceBack,"NoSlip" if (0.0..=1.0).contains(&friction)=>dpd::WallKind::NoSlip{friction},_=>return Err(api::error("wall kind/friction invalid"))};s.native.set_walls(dpd::DpdWall{axis,a_wall,kind}).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(dpd_len,"num_native_dpd_len","crate::sim::dpd::Dpd::len","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->usize,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(s.native.len()))});
object_binding!(dpd_is_empty,"num_native_dpd_is_empty","crate::sim::dpd::Dpd::is_empty","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(s.native.is_empty()))});
object_binding!(dpd_density,"num_native_dpd_density","crate::sim::dpd::Dpd::density","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.density(),[0,-(s.native.params.dim as i32),0,0],false)?)))});
object_binding!(dpd_momentum,"num_native_dpd_momentum","crate::sim::dpd::Dpd::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?)))});
object_binding!(dpd_temperature,"num_native_dpd_temperature","crate::sim::dpd::Dpd::temperature","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.temperature(),ENERGY,true)?)))});
object_binding!(dpd_pressure,"num_native_dpd_pressure","crate::sim::dpd::Dpd::pressure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.pressure(),[1,2-s.native.params.dim as i32,-2,0],false)?)))});
object_binding!(dpd_mean_velocity,"num_native_dpd_mean_velocity","crate::sim::dpd::Dpd::mean_velocity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>,axis:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>,axis:usize|{fluid.with(ctx,|s|{if axis>=s.native.params.dim{return Err(api::error("mean velocity axis outside active dimension"));}Ok(Dynamic(s.frame.output(s.native.mean_velocity(axis),VELOCITY,false)?))})});
object_binding!(dpd_step,"num_native_dpd_step","crate::sim::dpd::Dpd::step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>,dt:Input)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>,dt:Input|{fluid.with_mut(ctx,|s|{let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0{return Err(api::error("DPD time step must be positive"));}fluid_steps(s.native.len(),1)?;s.native.step(dt);s.validate()?;Ok(true)})});
object_binding!(dpd_run,"num_native_dpd_run","crate::sim::dpd::Dpd::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>,dt:Input,steps:usize)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>,dt:Input,steps:usize|{fluid.with_mut(ctx,|s|{let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0{return Err(api::error("DPD time step must be positive"));}fluid_steps(s.native.len(),steps)?;s.native.run(dt,steps);s.validate()?;Ok(true)})});
object_binding!(dpd_state,"num_native_dpd_state","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Dpd>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Dpd>|{fluid.with(ctx,|s|Ok(Dynamic(api::record(vec![("positions",engine::packed(&s.frame,s.native.pos.clone(),LENGTH)?),("velocities",engine::packed(&s.frame,s.native.vel.clone(),VELOCITY)?),("time",s.frame.output(s.native.t,TIME,false)?)]))))});
object_binding!(sph_kernel_w,"num_native_sph_kernel_w","crate::sim::sph::KernelKind::w","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,kind:String,dimension:usize,radius:Input,smoothing_length:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,kind:String,dimension:usize,radius:Input,smoothing_length:Input|{if !(1..=3).contains(&dimension){return Err(api::error("SPH kernel dimension must be1..3"));}let f=frame(ctx,&f)?;let r=f.scalar(&radius.0,LENGTH,false)?;let h=f.scalar(&smoothing_length.0,LENGTH,false)?;if r<0.0||h<=0.0{return Err(api::error("SPH kernel radius/smoothing domain invalid"));}let out=kernel(&kind)?.w(dimension,r,h);Ok(Dynamic(f.output(out,[0,-(dimension as i32),0,0],false)?))});
object_binding!(sph_kernel_dw,"num_native_sph_kernel_dw","crate::sim::sph::KernelKind::dw","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,kind:String,dimension:usize,radius:Input,smoothing_length:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,kind:String,dimension:usize,radius:Input,smoothing_length:Input|{if !(1..=3).contains(&dimension){return Err(api::error("SPH kernel dimension must be1..3"));}let f=frame(ctx,&f)?;let r=f.scalar(&radius.0,LENGTH,false)?;let h=f.scalar(&smoothing_length.0,LENGTH,false)?;if r<0.0||h<=0.0{return Err(api::error("SPH kernel radius/smoothing domain invalid"));}let out=kernel(&kind)?.dw(dimension,r,h);Ok(Dynamic(f.output(out,[0,-(dimension as i32)-1,0,0],false)?))});
object_binding!(sph_new_d1,"num_native_sph_new_d1","crate::sim::sph::Sph<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input)->NewObject<Gas<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input|{let f=frame(ctx,&f)?;let pos=engine::rows::<1>(&f,&positions.0,LENGTH)?;let mass=f.vector(&masses.0,MASS)?;let u=f.vector(&specific_energy.0,[0,2,-2,0])?;if pos.is_empty()||pos.len()>128||mass.len()!=pos.len()||u.len()!=pos.len()||mass.iter().any(|x|*x<=0.0)||u.iter().any(|x|*x<0.0){return Err(api::error("SPH needs1..128positions and positive masses/nonnegative specificenergy"));}let lo=engine::array::<1>(&f,&lower.0,LENGTH)?;let hi=engine::array::<1>(&f,&upper.0,LENGTH)?;if lo.iter().zip(&hi).any(|(lo,hi)|*lo>=*hi){return Err(api::error("SPH domain bounds must increase"));}let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("SPH boundary vector required"));};
if t.shape!=[1]{return Err(api::error("SPH boundary count differs fromD"));}let bc=t.data.iter().map(boundary).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("SPH boundary shape"))?;let native=sph::Sph::new(pos,mass,u,lo,hi,bc);let dims=f.dims.to_vec();Ok(NewObject::new(Gas{native,frame:f},dims))});
object_binding!(sph_update_density_d1,"num_native_sph_update_density_d1","crate::sim::sph::Sph<D>::update_density","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with_mut(ctx,|s|{s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_derivatives_d1,"num_native_sph_derivatives_d1","crate::sim::sph::Sph<D>::derivatives","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with_mut(ctx,|s|{let(a,u)=s.native.derivatives();s.validate()?;Ok(Dynamic(api::record(vec![("accelerations",engine::packed(&s.frame,a,[0,1,-2,0])?),("energy_rates",s.frame.array(u,[0,2,-3,0])?)])))})});
object_binding!(sph_total_energy_d1,"num_native_sph_total_energy_d1","crate::sim::sph::Sph<D>::total_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(sph_gravitational_energy_d1,"num_native_sph_gravitational_energy_d1","crate::sim::sph::Sph<D>::gravitational_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with_mut(ctx,|s|{let v=s.native.gravitational_energy().map(|x|s.frame.output(x,ENERGY,false)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if v.is_some(){"Some"}else{"None"}.into(),payload:v.into_iter().collect()})))})});
object_binding!(sph_momentum_d1,"num_native_sph_momentum_d1","crate::sim::sph::Sph<D>::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?))})});
object_binding!(sph_pressure_d1,"num_native_sph_pressure_d1","crate::sim::sph::Sph<D>::pressure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.pressure(index),[1,1,-2,0],false)?))})});
object_binding!(sph_sound_speed_d1,"num_native_sph_sound_speed_d1","crate::sim::sph::Sph<D>::sound_speed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.sound_speed(index),VELOCITY,false)?))})});
object_binding!(sph_step_d1,"num_native_sph_step_d1","crate::sim::sph::Sph<D>::step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,dt:Input)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,dt:Input|{fluid.with_mut(ctx,|s|{let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0{return Err(api::error("SPH time step must bepositive"));}s.native.step(dt);s.validate()?;Ok(true)})});
object_binding!(sph_suggest_dt_d1,"num_native_sph_suggest_dt_d1","crate::sim::sph::Sph<D>::suggest_dt","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,cfl:f64)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,cfl:f64|{if cfl<=0.0||cfl>1.0{return Err(api::error("SPH CFL must be in(0,1]"));}fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.suggest_dt(cfl),TIME,false)?)))});
object_binding!(sph_state_d1,"num_native_sph_state_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>|{fluid.with(ctx,|s|Ok(Dynamic(api::record(vec![("positions",engine::packed(&s.frame,s.native.pos.clone(),LENGTH)?),("velocities",engine::packed(&s.frame,s.native.vel.clone(),VELOCITY)?),("specific_energy",s.frame.array(s.native.u.clone(),[0,2,-2,0])?),("density",s.frame.array(s.native.rho.clone(),[1,-1,0,0])?),("smoothing_length",s.frame.array(s.native.h.clone(),LENGTH)?)]))))});
object_binding!(sph_controls_d1,"num_native_sph_controls_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool|{if gamma<=1.0||gamma>3.0||eta<=0.0||eta>4.0||alpha<0.0||alpha>10.0||beta<0.0||beta>20.0{return Err(api::error("SPH equationofstate/smoothing/viscosity controls invalid"));}let k=kernel(&kind)?;fluid.with_mut(ctx,|s|{s.native.gamma=gamma;s.native.eta=eta;s.native.alpha=alpha;s.native.beta=beta;s.native.kernel=k;s.native.grad_h=grad_h;s.native.balsara=balsara;s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_set_gravity_d1,"num_native_sph_set_gravity_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<1>>,constant:Input,softening:Input,theta:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<1>>,constant:Input,softening:Input,theta:f64|{if theta<=0.0||theta>1.0{return Err(api::error("gravity opening angle outside(0,1]"));}fluid.with_mut(ctx,|s|{let g=s.frame.scalar(&constant.0,[-1,1,-2,0],false)?;let eps=s.frame.scalar(&softening.0,LENGTH,false)?;if g<0.0||eps<0.0{return Err(api::error("gravity coupling/softening must benonnegative"));}s.native.gravity=Some(sph::SphGravity{g,eps,theta});Ok(true)})});
object_binding!(sph_new_d2,"num_native_sph_new_d2","crate::sim::sph::Sph<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input)->NewObject<Gas<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input|{let f=frame(ctx,&f)?;let pos=engine::rows::<2>(&f,&positions.0,LENGTH)?;let mass=f.vector(&masses.0,MASS)?;let u=f.vector(&specific_energy.0,[0,2,-2,0])?;if pos.is_empty()||pos.len()>128||mass.len()!=pos.len()||u.len()!=pos.len()||mass.iter().any(|x|*x<=0.0)||u.iter().any(|x|*x<0.0){return Err(api::error("SPH needs1..128positions and positive masses/nonnegative specificenergy"));}let lo=engine::array::<2>(&f,&lower.0,LENGTH)?;let hi=engine::array::<2>(&f,&upper.0,LENGTH)?;if lo.iter().zip(&hi).any(|(lo,hi)|*lo>=*hi){return Err(api::error("SPH domain bounds must increase"));}let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("SPH boundary vector required"));};
if t.shape!=[2]{return Err(api::error("SPH boundary count differs fromD"));}let bc=t.data.iter().map(boundary).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("SPH boundary shape"))?;let native=sph::Sph::new(pos,mass,u,lo,hi,bc);let dims=f.dims.to_vec();Ok(NewObject::new(Gas{native,frame:f},dims))});
object_binding!(sph_update_density_d2,"num_native_sph_update_density_d2","crate::sim::sph::Sph<D>::update_density","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with_mut(ctx,|s|{s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_derivatives_d2,"num_native_sph_derivatives_d2","crate::sim::sph::Sph<D>::derivatives","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with_mut(ctx,|s|{let(a,u)=s.native.derivatives();s.validate()?;Ok(Dynamic(api::record(vec![("accelerations",engine::packed(&s.frame,a,[0,1,-2,0])?),("energy_rates",s.frame.array(u,[0,2,-3,0])?)])))})});
object_binding!(sph_total_energy_d2,"num_native_sph_total_energy_d2","crate::sim::sph::Sph<D>::total_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(sph_gravitational_energy_d2,"num_native_sph_gravitational_energy_d2","crate::sim::sph::Sph<D>::gravitational_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with_mut(ctx,|s|{let v=s.native.gravitational_energy().map(|x|s.frame.output(x,ENERGY,false)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if v.is_some(){"Some"}else{"None"}.into(),payload:v.into_iter().collect()})))})});
object_binding!(sph_momentum_d2,"num_native_sph_momentum_d2","crate::sim::sph::Sph<D>::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?))})});
object_binding!(sph_pressure_d2,"num_native_sph_pressure_d2","crate::sim::sph::Sph<D>::pressure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.pressure(index),[1,0,-2,0],false)?))})});
object_binding!(sph_sound_speed_d2,"num_native_sph_sound_speed_d2","crate::sim::sph::Sph<D>::sound_speed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.sound_speed(index),VELOCITY,false)?))})});
object_binding!(sph_step_d2,"num_native_sph_step_d2","crate::sim::sph::Sph<D>::step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,dt:Input)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,dt:Input|{fluid.with_mut(ctx,|s|{let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0{return Err(api::error("SPH time step must bepositive"));}s.native.step(dt);s.validate()?;Ok(true)})});
object_binding!(sph_suggest_dt_d2,"num_native_sph_suggest_dt_d2","crate::sim::sph::Sph<D>::suggest_dt","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,cfl:f64)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,cfl:f64|{if cfl<=0.0||cfl>1.0{return Err(api::error("SPH CFL must be in(0,1]"));}fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.suggest_dt(cfl),TIME,false)?)))});
object_binding!(sph_state_d2,"num_native_sph_state_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>|{fluid.with(ctx,|s|Ok(Dynamic(api::record(vec![("positions",engine::packed(&s.frame,s.native.pos.clone(),LENGTH)?),("velocities",engine::packed(&s.frame,s.native.vel.clone(),VELOCITY)?),("specific_energy",s.frame.array(s.native.u.clone(),[0,2,-2,0])?),("density",s.frame.array(s.native.rho.clone(),[1,-2,0,0])?),("smoothing_length",s.frame.array(s.native.h.clone(),LENGTH)?)]))))});
object_binding!(sph_controls_d2,"num_native_sph_controls_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool|{if gamma<=1.0||gamma>3.0||eta<=0.0||eta>4.0||alpha<0.0||alpha>10.0||beta<0.0||beta>20.0{return Err(api::error("SPH equationofstate/smoothing/viscosity controls invalid"));}let k=kernel(&kind)?;fluid.with_mut(ctx,|s|{s.native.gamma=gamma;s.native.eta=eta;s.native.alpha=alpha;s.native.beta=beta;s.native.kernel=k;s.native.grad_h=grad_h;s.native.balsara=balsara;s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_set_gravity_d2,"num_native_sph_set_gravity_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<2>>,constant:Input,softening:Input,theta:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<2>>,constant:Input,softening:Input,theta:f64|{if theta<=0.0||theta>1.0{return Err(api::error("gravity opening angle outside(0,1]"));}fluid.with_mut(ctx,|s|{let g=s.frame.scalar(&constant.0,[-1,2,-2,0],false)?;let eps=s.frame.scalar(&softening.0,LENGTH,false)?;if g<0.0||eps<0.0{return Err(api::error("gravity coupling/softening must benonnegative"));}s.native.gravity=Some(sph::SphGravity{g,eps,theta});Ok(true)})});
object_binding!(sph_new_d3,"num_native_sph_new_d3","crate::sim::sph::Sph<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input)->NewObject<Gas<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,specific_energy:Input,lower:Input,upper:Input,boundaries:Input|{let f=frame(ctx,&f)?;let pos=engine::rows::<3>(&f,&positions.0,LENGTH)?;let mass=f.vector(&masses.0,MASS)?;let u=f.vector(&specific_energy.0,[0,2,-2,0])?;if pos.is_empty()||pos.len()>128||mass.len()!=pos.len()||u.len()!=pos.len()||mass.iter().any(|x|*x<=0.0)||u.iter().any(|x|*x<0.0){return Err(api::error("SPH needs1..128positions and positive masses/nonnegative specificenergy"));}let lo=engine::array::<3>(&f,&lower.0,LENGTH)?;let hi=engine::array::<3>(&f,&upper.0,LENGTH)?;if lo.iter().zip(&hi).any(|(lo,hi)|*lo>=*hi){return Err(api::error("SPH domain bounds must increase"));}let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("SPH boundary vector required"));};
if t.shape!=[3]{return Err(api::error("SPH boundary count differs fromD"));}let bc=t.data.iter().map(boundary).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("SPH boundary shape"))?;let native=sph::Sph::new(pos,mass,u,lo,hi,bc);let dims=f.dims.to_vec();Ok(NewObject::new(Gas{native,frame:f},dims))});
object_binding!(sph_update_density_d3,"num_native_sph_update_density_d3","crate::sim::sph::Sph<D>::update_density","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with_mut(ctx,|s|{s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_derivatives_d3,"num_native_sph_derivatives_d3","crate::sim::sph::Sph<D>::derivatives","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with_mut(ctx,|s|{let(a,u)=s.native.derivatives();s.validate()?;Ok(Dynamic(api::record(vec![("accelerations",engine::packed(&s.frame,a,[0,1,-2,0])?),("energy_rates",s.frame.array(u,[0,2,-3,0])?)])))})});
object_binding!(sph_total_energy_d3,"num_native_sph_total_energy_d3","crate::sim::sph::Sph<D>::total_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(sph_gravitational_energy_d3,"num_native_sph_gravitational_energy_d3","crate::sim::sph::Sph<D>::gravitational_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with_mut(ctx,|s|{let v=s.native.gravitational_energy().map(|x|s.frame.output(x,ENERGY,false)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if v.is_some(){"Some"}else{"None"}.into(),payload:v.into_iter().collect()})))})});
object_binding!(sph_momentum_d3,"num_native_sph_momentum_d3","crate::sim::sph::Sph<D>::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with_mut(ctx,|s|{Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?))})});
object_binding!(sph_pressure_d3,"num_native_sph_pressure_d3","crate::sim::sph::Sph<D>::pressure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.pressure(index),[1,-1,-2,0],false)?))})});
object_binding!(sph_sound_speed_d3,"num_native_sph_sound_speed_d3","crate::sim::sph::Sph<D>::sound_speed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,index:usize|{fluid.with(ctx,|s|{if index>=s.native.pos.len(){return Err(api::error("SPH particle index outside state"));}Ok(Dynamic(s.frame.output(s.native.sound_speed(index),VELOCITY,false)?))})});
object_binding!(sph_step_d3,"num_native_sph_step_d3","crate::sim::sph::Sph<D>::step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,dt:Input)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,dt:Input|{fluid.with_mut(ctx,|s|{let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0{return Err(api::error("SPH time step must bepositive"));}s.native.step(dt);s.validate()?;Ok(true)})});
object_binding!(sph_suggest_dt_d3,"num_native_sph_suggest_dt_d3","crate::sim::sph::Sph<D>::suggest_dt","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,cfl:f64)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,cfl:f64|{if cfl<=0.0||cfl>1.0{return Err(api::error("SPH CFL must be in(0,1]"));}fluid.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.suggest_dt(cfl),TIME,false)?)))});
object_binding!(sph_state_d3,"num_native_sph_state_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>)->Dynamic,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>|{fluid.with(ctx,|s|Ok(Dynamic(api::record(vec![("positions",engine::packed(&s.frame,s.native.pos.clone(),LENGTH)?),("velocities",engine::packed(&s.frame,s.native.vel.clone(),VELOCITY)?),("specific_energy",s.frame.array(s.native.u.clone(),[0,2,-2,0])?),("density",s.frame.array(s.native.rho.clone(),[1,-3,0,0])?),("smoothing_length",s.frame.array(s.native.h.clone(),LENGTH)?)]))))});
object_binding!(sph_controls_d3,"num_native_sph_controls_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,gamma:f64,eta:f64,alpha:f64,beta:f64,kind:String,grad_h:bool,balsara:bool|{if gamma<=1.0||gamma>3.0||eta<=0.0||eta>4.0||alpha<0.0||alpha>10.0||beta<0.0||beta>20.0{return Err(api::error("SPH equationofstate/smoothing/viscosity controls invalid"));}let k=kernel(&kind)?;fluid.with_mut(ctx,|s|{s.native.gamma=gamma;s.native.eta=eta;s.native.alpha=alpha;s.native.beta=beta;s.native.kernel=k;s.native.grad_h=grad_h;s.native.balsara=balsara;s.native.update_density();s.validate()?;Ok(true)})});
object_binding!(sph_set_gravity_d3,"num_native_sph_set_gravity_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fluid:Handle<Gas<3>>,constant:Input,softening:Input,theta:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Gas<3>>,constant:Input,softening:Input,theta:f64|{if theta<=0.0||theta>1.0{return Err(api::error("gravity opening angle outside(0,1]"));}fluid.with_mut(ctx,|s|{let g=s.frame.scalar(&constant.0,[-1,3,-2,0],false)?;let eps=s.frame.scalar(&softening.0,LENGTH,false)?;if g<0.0||eps<0.0{return Err(api::error("gravity coupling/softening must benonnegative"));}s.native.gravity=Some(sph::SphGravity{g,eps,theta});Ok(true)})});
struct Mm {
    native: mm::SPHSystem,
}
impl NativeType for Mm {
    const KIND: &'static str = "rssn/native/sim/legacy-sph-normalized/v1";
    const MAX_BYTES: usize = 524288;

    fn validate(&self) -> Result<(), OpError> {
        engine::bound(self.native.particles.len(), 128)?;
        for p in &self.native.particles {
            engine::finite([
                p.pos.x, p.pos.y, p.vel.x, p.vel.y, p.force.x, p.force.y, p.mass, p.density,
                p.pressure,
            ])?;
        }
        Ok(())
    }
}
impl NativeType for mm::Vector2D {
    const KIND: &'static str = "rssn/native/sim/mm-vector-normalized/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        engine::finite([self.x, self.y])
    }
}
impl NativeType for mm::Poly6Kernel {
    const KIND: &'static str = "rssn/native/sim/mm-poly6-normalized/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if self.h_sq <= 0.0 {
            return Err(api::error("normalized Poly6 smoothing scale underflowed"));
        }
        engine::finite([self.h_sq, self.factor])
    }
}
impl NativeType for mm::SpikyKernel {
    const KIND: &'static str = "rssn/native/sim/mm-spiky-normalized/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if self.h <= 0.0 {
            return Err(api::error("normalized Spiky smoothing scale invalid"));
        }
        engine::finite([self.h, self.factor])
    }
}
object_binding!(mm_vector_new,"num_native_mm_vector_new","crate::sim::physics_mm::Vector2D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(x:f64,y:f64)->NewObject<mm::Vector2D>,|_ctx:&CallCtx<'_>,x:f64,y:f64|{Ok(NewObject::normalized(mm::Vector2D::new(x,y)))});
object_binding!(mm_poly6_new,"num_native_mm_poly6_new","crate::sim::physics_mm::Poly6Kernel::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(smoothing_length:f64)->NewObject<mm::Poly6Kernel>,|_ctx:&CallCtx<'_>,smoothing_length:f64|{if smoothing_length<=0.0{return Err(api::error("normalized smoothing length must bepositive"));}Ok(NewObject::normalized(mm::Poly6Kernel::new(smoothing_length)))});
object_binding!(mm_spiky_new,"num_native_mm_spiky_new","crate::sim::physics_mm::SpikyKernel::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(smoothing_length:f64)->NewObject<mm::SpikyKernel>,|_ctx:&CallCtx<'_>,smoothing_length:f64|{if smoothing_length<=0.0{return Err(api::error("normalized smoothing length must bepositive"));}Ok(NewObject::normalized(mm::SpikyKernel::new(smoothing_length)))});
object_binding!(mm_system_new,"num_native_mm_system_new","","Legacy native educational SPH uses explicitly normalized unit-mass coefficients and3D smoothing kernels on2D sites. This adapter rejects physical quantities rather than assigning inconsistent SI force dimensions; use the general Sph<D> adapters for coherent physical models.",true,(positions:Vec<[f64;2]>,smoothing_length:f64,mass:f64,gravity:[f64;2],viscosity:f64,gas_constant:f64,rest_density:f64,bounds:[f64;2])->NewObject<Mm>,|_ctx:&CallCtx<'_>,positions:Vec<[f64;2]>,smoothing_length:f64,mass:f64,gravity:[f64;2],viscosity:f64,gas_constant:f64,rest_density:f64,bounds:[f64;2]|{if positions.is_empty()||positions.len()>128||smoothing_length<=0.0||mass<=0.0||viscosity<0.0||gas_constant<=0.0||rest_density<=0.0||bounds.iter().any(|x|*x<=0.0){return Err(api::error("legacy reduced-unit SPH controls invalid"));}let particles=positions.into_iter().map(|[x,y]|mm::Particle{pos:mm::Vector2D::new(x,y),vel:mm::Vector2D::new(0.0,0.0),force:mm::Vector2D::new(0.0,0.0),density:0.0,pressure:0.0,mass}).collect();let mut native=mm::SPHSystem{particles,poly6:mm::Poly6Kernel::new(smoothing_length),spiky:mm::SpikyKernel::new(smoothing_length),gravity:mm::Vector2D::new(gravity[0],gravity[1]),viscosity,gas_const:gas_constant,rest_density,bounds:mm::Vector2D::new(bounds[0],bounds[1])};native.compute_density_pressure();Ok(NewObject::normalized(Mm{native}))});
object_binding!(mm_compute_density_pressure,"num_native_mm_compute_density_pressure","crate::sim::physics_mm::SPHSystem::compute_density_pressure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(fluid:Handle<Mm>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Mm>|{fluid.with_mut(ctx,|s|{s.native.compute_density_pressure();s.validate()?;Ok(true)})});
object_binding!(mm_compute_forces,"num_native_mm_compute_forces","crate::sim::physics_mm::SPHSystem::compute_forces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(fluid:Handle<Mm>)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Mm>|{fluid.with_mut(ctx,|s|{if s.native.particles.iter().any(|p|p.density<=0.0){return Err(api::error("legacy SPH density must becomputed andpositive"));}s.native.compute_forces();s.validate()?;Ok(true)})});
object_binding!(mm_integrate,"num_native_mm_integrate","crate::sim::physics_mm::SPHSystem::integrate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(fluid:Handle<Mm>,dt:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Mm>,dt:f64|{if dt<=0.0{return Err(api::error("legacy reducedtime step must bepositive"));}fluid.with_mut(ctx,|s|{if s.native.particles.iter().any(|p|p.density<=0.0){return Err(api::error("legacy SPH density must becomputed andpositive"));}s.native.integrate(dt);s.validate()?;Ok(true)})});
object_binding!(mm_update,"num_native_mm_update","crate::sim::physics_mm::SPHSystem::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(fluid:Handle<Mm>,dt:f64)->bool,|ctx:&CallCtx<'_>,fluid:Handle<Mm>,dt:f64|{if dt<=0.0{return Err(api::error("legacy reducedtime step must bepositive"));}fluid.with_mut(ctx,|s|{if s.native.particles.iter().any(|p|p.density<=0.0){return Err(api::error("legacy SPH density must becomputed andpositive"));}s.native.update(dt);s.validate()?;Ok(true)})});
object_binding!(mm_dam_break_scenario,"num_native_mm_dam_break_scenario","crate::sim::physics_mm::simulate_dam_break_2d_scenario","Execute the original hard-coded native demonstration in its explicitly normalized model coordinates. No physical-unit interpretation is inferred.",true,()->Vec<[f64;2]>,|_ctx:&CallCtx<'_>|{Ok(mm::simulate_dam_break_2d_scenario().into_iter().map(|(x,y)|[x,y]).collect())});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &dpd_new::BINDING,
            &dpd_set_walls::BINDING,
            &dpd_len::BINDING,
            &dpd_is_empty::BINDING,
            &dpd_density::BINDING,
            &dpd_momentum::BINDING,
            &dpd_temperature::BINDING,
            &dpd_pressure::BINDING,
            &dpd_mean_velocity::BINDING,
            &dpd_step::BINDING,
            &dpd_run::BINDING,
            &dpd_state::BINDING,
            &sph_kernel_w::BINDING,
            &sph_kernel_dw::BINDING,
            &sph_new_d1::BINDING,
            &sph_update_density_d1::BINDING,
            &sph_derivatives_d1::BINDING,
            &sph_total_energy_d1::BINDING,
            &sph_gravitational_energy_d1::BINDING,
            &sph_momentum_d1::BINDING,
            &sph_pressure_d1::BINDING,
            &sph_sound_speed_d1::BINDING,
            &sph_step_d1::BINDING,
            &sph_suggest_dt_d1::BINDING,
            &sph_state_d1::BINDING,
            &sph_controls_d1::BINDING,
            &sph_set_gravity_d1::BINDING,
            &sph_new_d2::BINDING,
            &sph_update_density_d2::BINDING,
            &sph_derivatives_d2::BINDING,
            &sph_total_energy_d2::BINDING,
            &sph_gravitational_energy_d2::BINDING,
            &sph_momentum_d2::BINDING,
            &sph_pressure_d2::BINDING,
            &sph_sound_speed_d2::BINDING,
            &sph_step_d2::BINDING,
            &sph_suggest_dt_d2::BINDING,
            &sph_state_d2::BINDING,
            &sph_controls_d2::BINDING,
            &sph_set_gravity_d2::BINDING,
            &sph_new_d3::BINDING,
            &sph_update_density_d3::BINDING,
            &sph_derivatives_d3::BINDING,
            &sph_total_energy_d3::BINDING,
            &sph_gravitational_energy_d3::BINDING,
            &sph_momentum_d3::BINDING,
            &sph_pressure_d3::BINDING,
            &sph_sound_speed_d3::BINDING,
            &sph_step_d3::BINDING,
            &sph_suggest_dt_d3::BINDING,
            &sph_state_d3::BINDING,
            &sph_controls_d3::BINDING,
            &sph_set_gravity_d3::BINDING,
            &mm_vector_new::BINDING,
            &mm_poly6_new::BINDING,
            &mm_spiky_new::BINDING,
            &mm_system_new::BINDING,
            &mm_compute_density_pressure::BINDING,
            &mm_compute_forces::BINDING,
            &mm_integrate::BINDING,
            &mm_update::BINDING,
            &mm_dam_break_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
