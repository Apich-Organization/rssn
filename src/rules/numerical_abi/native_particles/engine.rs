#![allow(clippy::wildcard_imports)] // Physical D=1/2/3 native engine specializations share the parent reference frame.
use super::*;
use crate::kernels::random::Rng;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::sim::md_engine as e;
pub(super) struct BoxState<const D: usize> {
    pub(super) native: e::SimBox<D>,
    pub(super) frame: Frame,
}
pub(super) struct System<const D: usize> {
    pub(super) native: e::MdSystem<D>,
    pub(super) frame: Frame,
}
pub(super) struct Field<const D: usize> {
    pub(super) native: Option<e::ForceField<D>>,
    pub(super) extra_work: usize,
    pub(super) frame: Frame,
}
struct List<const D: usize> {
    native: e::NeighbourList<D>,
    frame: Frame,
}
struct Simulation<const D: usize> {
    extra_work: usize,
    native: e::MdSimulation<D>,
    frame: Frame,
}
struct ForceState<const D: usize> {
    native: e::Forces<D>,
    frame: Frame,
}
struct Radial<const D: usize> {
    native: e::Rdf<D>,
    frame: Frame,
    frames: usize,
}

impl<const D: usize> Radial<D> {
    const MAX_FRAMES: usize = 1_000_000;

    fn merge(
        &mut self,
        source: &e::Rdf<D>,
        frames: usize,
    ) -> Result<bool, OpError> {
        let frames = self
            .frames
            .checked_add(frames)
            .filter(|&frames| frames <= Self::MAX_FRAMES)
            .ok_or_else(|| api::error("RDF frame-count budget exhausted"))?;
        self.native.merge(source).map_err(api::error)?;
        self.frames = frames;
        Ok(true)
    }
}

fn merge_radial<const D: usize>(
    ctx: &CallCtx<'_>,
    destination: &Handle<Radial<D>>,
    source: &Handle<Radial<D>>,
) -> Result<bool, OpError> {
    if destination.value == source.value {
        return destination.with_mut(ctx, |rdf| rdf.merge(&rdf.native.clone(), rdf.frames));
    }
    // Take a consistent snapshot before locking the destination.
    let (native, frame, frames) = source.with(ctx, |rdf| {
        Ok((rdf.native.clone(), rdf.frame.clone(), rdf.frames))
    })?;
    destination.with_mut(ctx, |rdf| {
        rdf.frame.same(&frame)?;
        rdf.merge(&native, frames)
    })
}

pub(super) struct Pair {
    pub(super) native: e::PairTerm,
    pub(super) frame: Frame,
}
struct Chain {
    native: e::NoseHooverChain,
    frame: Frame,
}
pub(super) fn bound(
    n: usize,
    max: usize,
) -> Result<(), OpError> {
    if n > max {
        Err(api::error("native engine allocation/work bound exceeded"))
    } else {
        Ok(())
    }
}
pub(super) fn finite(values: impl IntoIterator<Item = f64>) -> Result<(), OpError> {
    if values.into_iter().any(|x| !x.is_finite()) {
        Err(api::error("native engine state contains nonfinite values"))
    } else {
        Ok(())
    }
}
pub(super) fn array<const D: usize>(
    f: &Frame,
    v: &Value,
    p: [i32; 4],
) -> Result<[f64; D], OpError> {
    f.vector(v, p)?
        .try_into()
        .map_err(|_| api::error("native engine vector dimension mismatch"))
}
pub(super) fn rows<const D: usize>(
    f: &Frame,
    v: &Value,
    p: [i32; 4],
) -> Result<Vec<[f64; D]>, OpError> {
    let Value::Tensor(t) = v else {
        return Err(api::error("native engine matrix expected"));
    };
    if t.shape.len() != 2 || t.shape[1] != D || t.shape[0] > 512 || t.shape[0] * D != t.data.len() {
        return Err(api::error(
            "native engine matrix shape/particle bound exceeded",
        ));
    }
    api::same_dimension(&api::dimension(v)?, &f.dim(p)?)?;
    t.data
        .chunks(D)
        .map(|r| {
            r.iter()
                .map(|v| f.apply(crate::rules::numerical_abi::numeric(v, 0)?, p, false, false))
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|_| api::error("native engine row shape"))
        })
        .collect()
}
pub(super) fn packed<const D: usize>(
    f: &Frame,
    v: Vec<[f64; D]>,
    p: [i32; 4],
) -> Result<Value, OpError> {
    let n = v.len();
    let data = v
        .into_iter()
        .flatten()
        .map(|x| f.output(x, p, false))
        .collect::<Result<_, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape: vec![n, D],
        data,
        elem: TypeExpr::Float.with_dim(f.dim(p)?),
    }))
}
pub(super) fn matrix<const D: usize>(
    f: &Frame,
    v: &Value,
    p: [i32; 4],
) -> Result<[[f64; D]; D], OpError> {
    rows(f, v, p)?
        .try_into()
        .map_err(|_| api::error("native engine matrix must be square D-by-D"))
}
pub(super) fn check_box<const D: usize>(b: &e::SimBox<D>) -> Result<(), OpError> {
    if !(1..=3).contains(&D)
        || b.len.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || !b.volume().is_finite()
        || b.volume() <= 0.0
    {
        return Err(api::error(
            "box must have positive finite nonsingular lengths/volume",
        ));
    }
    finite(b.matrix().into_iter().flatten())
}
const PARTICLE_LIMIT: usize = 512;
const NEIGHBOUR_PAIR_LIMIT: usize = PARTICLE_LIMIT * (PARTICLE_LIMIT - 1) / 2;

pub(super) fn check_system<const D: usize>(s: &e::MdSystem<D>) -> Result<(), OpError> {
    bound(s.n(), PARTICLE_LIMIT)?;
    check_box(&s.bx)?;
    if s.n() == 0
        || s.vel.len() != s.n()
        || s.unwrapped.len() != s.n()
        || s.mass.len() != s.n()
        || s.charge.len() != s.n()
        || s.mass.iter().any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(api::error("native MD system shape/mass invalid"));
    }
    finite(
        s.pos
            .iter()
            .chain(&s.vel)
            .chain(&s.unwrapped)
            .flatten()
            .copied(),
    )?;
    finite(s.charge.iter().copied())?;
    finite([s.time])
}
macro_rules! object_type {
    ($ty:ident,$kind:literal,$bytes:expr,$value:ident,$body:expr) => {
        impl<const D: usize> NativeType for $ty<D> {
            const KIND: &'static str = match D {
                | 1 => concat!($kind, "/d1/v1"),
                | 2 => concat!($kind, "/d2/v1"),
                | _ => concat!($kind, "/d3/v1"),
            };
            const MAX_BYTES: usize = $bytes;

            fn validate(&self) -> Result<(), OpError> {
                let $value = self;
                $body
            }
        }
    };
}
object_type!(
    BoxState,
    "rssn/native/sim/engine-box",
    4096,
    v,
    check_box(&v.native)
);
object_type!(
    System,
    "rssn/native/sim/engine-system",
    262144,
    v,
    check_system(&v.native)
);
object_type!(Field, "rssn/native/sim/engine-force-field", 2_097_152, v, {
    let ff = v
        .native
        .as_ref()
        .ok_or_else(|| api::error("force field ownership has been transferred"))?;
    bound(ff.pair.len(), 64)?;
    bound(ff.bonds.len() + ff.angles.len() + ff.dihedrals.len(), 4096)?;
    ff.validate().map_err(api::error)
});
object_type!(List, "rssn/native/sim/engine-neighbours", 4_194_304, v, {
    bound(v.native.pairs.len(), NEIGHBOUR_PAIR_LIMIT)?;
    if v.native.rcut <= 0.0 || v.native.skin < 0.0 {
        return Err(api::error("neighbour cutoff/skin invalid"));
    }
    Ok(())
});
object_type!(
    Simulation,
    "rssn/native/sim/engine-simulation",
    8_388_608,
    v,
    {
        check_system(&v.native.sys)?;
        if v.native.dt <= 0.0 || !v.native.dt.is_finite() {
            return Err(api::error("simulation time step invalid"));
        }
        Ok(())
    }
);
object_type!(ForceState, "rssn/native/sim/engine-forces", 262144, v, {
    bound(v.native.f.len(), 512)?;
    finite(v.native.f.iter().flatten().copied())?;
    finite(v.native.virial.into_iter().flatten())?;
    finite([v.native.epot])
});
object_type!(Radial, "rssn/native/sim/engine-rdf", 65536, v, {
    bound(v.native.hist.len(), 4096)?;
    if v.native.hist.is_empty() || v.native.rmax <= 0.0 {
        return Err(api::error("RDF bins/radius invalid"));
    }
    Ok(())
});
impl NativeType for Pair {
    const KIND: &'static str = "rssn/native/sim/engine-pair/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
impl NativeType for Chain {
    const KIND: &'static str = "rssn/native/sim/engine-thermostat-chain/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        bound(self.native.q.len(), 16)?;
        if self.native.q.is_empty()
            || self.native.q.iter().any(|x| *x <= 0.0)
            || self.native.nf <= 0.0
        {
            return Err(api::error("thermostat chain masses/degrees invalid"));
        }
        finite(
            self.native
                .q
                .iter()
                .chain(&self.native.xi)
                .chain(&self.native.vxi)
                .copied(),
        )
    }
}
pub(super) fn system_data<const D: usize>(
    s: &e::MdSystem<D>,
    f: &Frame,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("positions", packed(f, s.pos.clone(), LENGTH)?),
        ("velocities", packed(f, s.vel.clone(), VELOCITY)?),
        ("unwrapped", packed(f, s.unwrapped.clone(), LENGTH)?),
        ("masses", f.array(s.mass.clone(), MASS)?),
        ("charges", f.array(s.charge.clone(), CHARGE)?),
        ("time", f.output(s.time, TIME, false)?),
    ]))
}
pub(super) fn force_data<const D: usize>(
    s: &e::Forces<D>,
    f: &Frame,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("forces", packed(f, s.f.clone(), FORCE)?),
        ("energy", f.output(s.epot, ENERGY, false)?),
        ("virial", packed(f, s.virial.to_vec(), ENERGY)?),
    ]))
}
fn boundary(
    f: &Frame,
    v: &Value,
) -> Result<e::Boundary, OpError> {
    match v {
        | Value::String(s) => {
            match s.as_str() {
                | "Open" => Ok(e::Boundary::Open),
                | "Periodic" => Ok(e::Boundary::Periodic),
                | "Reflecting" => Ok(e::Boundary::Reflecting),
                | _ => {
                    Err(api::error(
                        "boundary must be Open/Periodic/Reflecting or Thermal record",
                    ))
                },
            }
        },
        | _ => {
            let t = f.scalar(api::field(v, "temperature")?, ENERGY, true)?;
            if t < 0.0 {
                return Err(api::error(
                    "thermal boundary temperature must be nonnegative",
                ));
            }
            Ok(e::Boundary::Thermal(t))
        },
    }
}
pub(super) fn safe_grid<const D: usize>(
    s: &e::MdSystem<D>,
    list: &e::NeighbourList<D>,
) -> Result<(), OpError> {
    // Dense neighbour storage is limited to eight bins per particle. Sparse
    // storage contains at most one bin per occupied cell, so box volume no
    // longer determines allocation. The possible pair count is bounded by n.
    bound(s.n(), PARTICLE_LIMIT)?;
    let r = list.rcut + list.skin;
    if r <= 0.0 || !r.is_finite() {
        return Err(api::error("cutoff+skin must be positive finite"));
    }
    Ok(())
}
object_binding!(engine_pair_new,"num_native_engine_pair_new","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,kind:String,parameters:Input)->NewObject<Pair>,|ctx:&CallCtx<'_>,f:Handle<Frame>,kind:String,parameters:Input|{let f=frame(ctx,&f)?;let field=|name,p|f.scalar(api::field(&parameters.0,name)?,p,false);let native=match kind.as_str(){"LennardJones"=>{let eps=field("epsilon",ENERGY)?;let sigma=field("sigma",LENGTH)?;if eps<0.0||sigma<=0.0{return Err(api::error("LJ epsilon/sigma invalid"));}e::PairTerm::LennardJones{eps,sigma}},"Morse"=>{let d=field("depth",ENERGY)?;let alpha=field("alpha",[0,-1,0,0])?;let r0=field("r0",LENGTH)?;if d<0.0||alpha<=0.0||r0<=0.0{return Err(api::error("Morse controls invalid"));}e::PairTerm::Morse{d,alpha,r0}},"Coulomb"=>{let k=field("constant",[1,3,-2,-2])?;let alpha=field("alpha",[0,-1,0,0])?;if alpha<0.0{return Err(api::error("Ewald alpha must be nonnegative"));}e::PairTerm::Coulomb{k,alpha}},_=>return Err(api::error("pair kind must be LennardJones/Morse/Coulomb"))};let dims=f.dims.to_vec();Ok(NewObject::new(Pair{native,frame:f},dims))});
object_binding!(engine_pair_eval,"num_native_engine_pair_eval","crate::sim::md_engine::PairTerm::eval","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(pair:Handle<Pair>,distance_squared:Input,charge_product:Input)->Dynamic,|ctx:&CallCtx<'_>,pair:Handle<Pair>,distance_squared:Input,charge_product:Input|{pair.with(ctx,|p|{let r=p.frame.scalar(&distance_squared.0,[0,2,0,0],false)?;let q=p.frame.scalar(&charge_product.0,[0,0,0,2],false)?;if r<=0.0{return Err(api::error("pair squared distance must be positive"));}let(u,k)=p.native.eval(r,q);Ok(Dynamic(api::record(vec![("energy",p.frame.output(u,ENERGY,false)?),("force_over_distance",p.frame.output(k,[1,0,-2,0],false)?)])))})});
object_binding!(engine_box_periodic_d1,"num_native_engine_box_periodic_d1","crate::sim::md_engine::SimBox<D>::periodic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input)->NewObject<BoxState<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input|{let f=frame(ctx,&f)?;let lengths=array::<1>(&f,&lengths.0,LENGTH)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::periodic(lengths),frame:f},dims))});
object_binding!(engine_box_new_d1,"num_native_engine_box_new_d1","crate::sim::md_engine::SimBox<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input,boundaries:Input)->NewObject<BoxState<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input,boundaries:Input|{let f=frame(ctx,&f)?;let lengths=array::<1>(&f,&lengths.0,LENGTH)?;let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("boundaries require a rank1 vector"));};
if t.shape!=[1]{return Err(api::error("boundary count must equalD"));}let bc=t.data.iter().map(|v|boundary(&f,v)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("boundary array mismatch"))?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::new(lengths,bc),frame:f},dims))});
object_binding!(engine_box_triclinic_d1,"num_native_engine_box_triclinic_d1","crate::sim::md_engine::SimBox<D>::triclinic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cell:Input)->NewObject<BoxState<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cell:Input|{let f=frame(ctx,&f)?;let h=matrix::<1>(&f,&cell.0,LENGTH)?;let native=e::SimBox::triclinic(h).map_err(api::error)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native,frame:f},dims))});
object_binding!(engine_box_set_cell_d1,"num_native_engine_box_set_cell_d1","crate::sim::md_engine::SimBox<D>::set_cell","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,cell:Input)->bool,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,cell:Input|{b.with_mut(ctx,|b|{let h=matrix::<1>(&b.frame,&cell.0,LENGTH)?;b.native.set_cell(h).map_err(api::error)?;b.validate()?;Ok(true)})});
object_binding!(engine_box_matrix_d1,"num_native_engine_box_matrix_d1","crate::sim::md_engine::SimBox<D>::matrix","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>|{b.with(ctx,|b|Ok(Dynamic(packed(&b.frame,b.native.matrix().to_vec(),LENGTH)?)))});
object_binding!(engine_box_volume_d1,"num_native_engine_box_volume_d1","crate::sim::md_engine::SimBox<D>::volume","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>|{b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.volume(),[0,1,0,0],false)?)))});
object_binding!(engine_box_height_d1,"num_native_engine_box_height_d1","crate::sim::md_engine::SimBox<D>::height","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,axis:usize)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,axis:usize|{if axis>=1{return Err(api::error("box axis outsideD"));}b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.height(axis),LENGTH,false)?)))});
object_binding!(engine_box_to_fractional_d1,"num_native_engine_box_to_fractional_d1","crate::sim::md_engine::SimBox<D>::to_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,point:Input)->[f64;1],|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,point:Input|{b.with(ctx,|b|Ok(b.native.to_fractional(&array::<1>(&b.frame,&point.0,LENGTH)?)))});
object_binding!(engine_box_from_fractional_d1,"num_native_engine_box_from_fractional_d1","crate::sim::md_engine::SimBox<D>::from_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,point:[f64;1])->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,point:[f64;1]|{b.with(ctx,|b|Ok(Dynamic(b.frame.array(b.native.from_fractional(&point).to_vec(),LENGTH)?)))});
object_binding!(engine_box_min_image_d1,"num_native_engine_box_min_image_d1","crate::sim::md_engine::SimBox<D>::min_image","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,displacement:Input)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,displacement:Input|{b.with(ctx,|b|{let mut x=array::<1>(&b.frame,&displacement.0,LENGTH)?;b.native.min_image(&mut x);Ok(Dynamic(b.frame.array(x.to_vec(),LENGTH)?))})});
object_binding!(engine_system_new_d1,"num_native_engine_system_new_d1","crate::sim::md_engine::MdSystem<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,positions:Input)->NewObject<System<1>>,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,positions:Input|{b.with(ctx,|b|{let positions=rows::<1>(&b.frame,&positions.0,LENGTH)?;if positions.is_empty(){return Err(api::error("system needs at least one particle"));}Ok(NewObject::new(System{native:e::MdSystem::new(positions,b.native.clone()),frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(engine_system_data_d1,"num_native_engine_system_data_d1","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<1>>|{system.with(ctx,|s|Ok(Dynamic(system_data(&s.native,&s.frame)?)))});
object_binding!(engine_system_set_state_d1,"num_native_engine_system_set_state_d1","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,velocities:Input,masses:Input,charges:Input)->bool,|ctx:&CallCtx<'_>,system:Handle<System<1>>,velocities:Input,masses:Input,charges:Input|{system.with_mut(ctx,|s|{let vel=rows::<1>(&s.frame,&velocities.0,VELOCITY)?;let mass=s.frame.vector(&masses.0,MASS)?;let charge=s.frame.vector(&charges.0,CHARGE)?;if vel.len()!=s.native.n()||mass.len()!=s.native.n()||charge.len()!=s.native.n()||mass.iter().any(|x|*x<=0.0){return Err(api::error("system state field shape/mass mismatch"));}s.native.vel=vel;s.native.mass=mass;s.native.charge=charge;s.validate()?;Ok(true)})});
object_binding!(engine_system_n_d1,"num_native_engine_system_n_d1","crate::sim::md_engine::MdSystem<D>::n","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>)->usize,|ctx:&CallCtx<'_>,system:Handle<System<1>>|{system.with(ctx,|s|Ok(s.native.n()))});
object_binding!(engine_system_kinetic_energy_d1,"num_native_engine_system_kinetic_energy_d1","crate::sim::md_engine::MdSystem<D>::kinetic_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<1>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.kinetic_energy(),ENERGY,false)?)))});
object_binding!(engine_system_momentum_d1,"num_native_engine_system_momentum_d1","crate::sim::md_engine::MdSystem<D>::momentum","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<1>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?)))});
object_binding!(engine_system_temperature_d1,"num_native_engine_system_temperature_d1","crate::sim::md_engine::MdSystem<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,dof:f64)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<1>>,dof:f64|{if dof<=0.0{return Err(api::error("thermal degrees must be positive"));}system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.temperature(dof),ENERGY,true)?)))});
object_binding!(engine_system_remove_com_velocity_d1,"num_native_engine_system_remove_com_velocity_d1","crate::sim::md_engine::MdSystem<D>::remove_com_velocity","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<1>>|{system.with_mut(ctx,|s|{s.native.remove_com_velocity();s.validate()?;Ok(true)})});
object_binding!(engine_system_maxwell_boltzmann_d1,"num_native_engine_system_maxwell_boltzmann_d1","crate::sim::md_engine::MdSystem<D>::maxwell_boltzmann","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,temperature:Input,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<1>>,temperature:Input,rng:Handle<Rng>|{system.with_mut(ctx,|s|{let t=s.frame.scalar(&temperature.0,ENERGY,true)?;if t<0.0{return Err(api::error("temperature must be nonnegative"));}rng.with_mut(ctx,|r|{s.native.maxwell_boltzmann(t,r);s.validate()?;Ok(true)})})});
object_binding!(engine_system_apply_boundaries_d1,"num_native_engine_system_apply_boundaries_d1","crate::sim::md_engine::MdSystem<D>::apply_boundaries","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<1>>,rng:Handle<Rng>|{system.with_mut(ctx,|s|{if s.native.pos.iter().flatten().any(|x|x.abs()>1e9*s.native.bx.len.iter().copied().fold(1.0,f64::min)){return Err(api::error("boundary displacement exceeds traversal work bound"));}rng.with_mut(ctx,|r|{s.native.apply_boundaries(r);s.validate()?;Ok(true)})})});
object_binding!(engine_lattice_cubic_d1,"num_native_engine_lattice_cubic_d1","crate::sim::md_engine::lattice_cubic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,counts:[usize;1],spacing:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,counts:[usize;1],spacing:Input|{let f=frame(ctx,&f)?;if counts.contains(&0)||counts.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).is_none_or(|n|n>512){return Err(api::error("lattice particle allocation bound is512"));}let a=f.scalar(&spacing.0,LENGTH,false)?;if a<=0.0{return Err(api::error("lattice spacing must be positive"));}Ok(Dynamic(packed(&f,e::lattice_cubic(counts,a),LENGTH)?))});
object_binding!(engine_random_positions_d1,"num_native_engine_random_positions_d1","crate::sim::md_engine::random_positions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<1>>,n:usize,minimum_distance:Input,rng:Handle<Rng>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<1>>,n:usize,minimum_distance:Input,rng:Handle<Rng>|{bound(n,128)?;b.with(ctx,|b|{let d=b.frame.scalar(&minimum_distance.0,LENGTH,false)?;if d<0.0{return Err(api::error("minimum separation must be nonnegative"));}let out=rng.with_mut(ctx,|r|Ok(e::random_positions(n,&b.native,d,r)))?;Ok(Dynamic(match out{None=>Value::Variant(crate::kernel::value::VariantValue{tag:"None".into(),payload:vec![]}),Some(v)=>Value::Variant(crate::kernel::value::VariantValue{tag:"Some".into(),payload:vec![packed(&b.frame,v,LENGTH)?]})}))})});
object_binding!(engine_forcefield_new_d1,"num_native_engine_forcefield_new_d1","crate::sim::md_engine::ForceField<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,terms:Input,cutoff:Input)->NewObject<Field<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,terms:Input,cutoff:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;if rcut<=0.0{return Err(api::error("force cutoff must be positive"));}let Value::Tensor(t)=&terms.0 else{return Err(api::error("pair terms must be rank1handles"));};
if t.shape!=[t.data.len()]||t.data.len()>64{return Err(api::error("pair term shape/count bound exceeded"));}let pair=t.data.iter().map(|v|{let h=<Handle<Pair> as NativeArg>::decode(v,0,ctx)?;h.with(ctx,|p|{f.same(&p.frame)?;Ok(p.native)})}).collect::<Result<_,_>>()?;let dims=f.dims.to_vec();Ok(NewObject::new(Field{native:Some(e::ForceField::new(pair,rcut)),extra_work:0,frame:f},dims))});
object_binding!(engine_forcefield_validate_d1,"num_native_engine_forcefield_validate_d1","crate::sim::md_engine::ForceField<D>::validate","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<1>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<1>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.validate().map_err(api::error)?;Ok(true)})});
object_binding!(engine_forcefield_rebuild_exclusions_d1,"num_native_engine_forcefield_rebuild_exclusions_d1","crate::sim::md_engine::ForceField<D>::rebuild_exclusions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<1>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<1>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.rebuild_exclusions();Ok(true)})});
object_binding!(engine_forces_zeros_d1,"num_native_engine_forces_zeros_d1","crate::sim::md_engine::Forces<D>::zeros","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,n:usize)->NewObject<ForceState<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,n:usize|{bound(n,512)?;let f=frame(ctx,&f)?;let dims=f.dims.to_vec();Ok(NewObject::new(ForceState{native:e::Forces::<1>::zeros(n),frame:f},dims))});
object_binding!(engine_forces_data_d1,"num_native_engine_forces_data_d1","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(forces:Handle<ForceState<1>>)->Dynamic,|ctx:&CallCtx<'_>,forces:Handle<ForceState<1>>|{forces.with(ctx,|s|Ok(Dynamic(force_data(&s.native,&s.frame)?)))});
object_binding!(engine_neighbours_new_d1,"num_native_engine_neighbours_new_d1","crate::sim::md_engine::NeighbourList<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cutoff:Input,skin:Input)->NewObject<List<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cutoff:Input,skin:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;let skin=f.scalar(&skin.0,LENGTH,false)?;if rcut<=0.0||skin<0.0{return Err(api::error("cutoff must bepositive andskin nonnegative"));}let dims=f.dims.to_vec();Ok(NewObject::new(List{native:e::NeighbourList::<1>::new(rcut,skin),frame:f},dims))});
object_binding!(engine_neighbours_needs_rebuild_d1,"num_native_engine_neighbours_needs_rebuild_d1","crate::sim::md_engine::NeighbourList<D>::needs_rebuild","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<1>>,system:Handle<System<1>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<1>>,system:Handle<System<1>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;Ok(l.native.needs_rebuild(&s.native))}))});
object_binding!(engine_neighbours_build_d1,"num_native_engine_neighbours_build_d1","crate::sim::md_engine::NeighbourList<D>::build","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<1>>,system:Handle<System<1>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<1>>,system:Handle<System<1>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;l.native.build(&s.native);l.validate()?;Ok(true)}))});
object_binding!(engine_compute_forces_d1,"num_native_engine_compute_forces_d1","crate::sim::md_engine::compute_forces","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,field:Handle<Field<1>>,list:Handle<List<1>>)->NewObject<ForceState<1>>,|ctx:&CallCtx<'_>,system:Handle<System<1>>,field:Handle<Field<1>>,list:Handle<List<1>>|{system.with(ctx,|s|field.with(ctx,|ff|list.with_mut(ctx,|l|{s.frame.same(&ff.frame)?;s.frame.same(&l.frame)?;safe_grid(&s.native,&l.native)?;field_indices(ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,s.native.n())?;let native=e::compute_forces(&s.native,ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,&mut l.native).map_err(api::error)?;Ok(NewObject::new(ForceState{native,frame:s.frame.clone()},s.frame.dims.to_vec()))})))});
object_binding!(engine_simulation_new_d1,"num_native_engine_simulation_new_d1","crate::sim::md_engine::MdSimulation<D>::new","Construct the actual native driver with an owned force field transferred once. The system is cloned; the force-field handle becomes consumed. Preflight checks run before the transfer; physical frame witnesses and compiled extra-force ownership stay in the simulation.",false,(system:Handle<System<1>>,field:Handle<Field<1>>,skin:Input,dt:Input,seed:u64)->NewObject<Simulation<1>>,|ctx:&CallCtx<'_>,system:Handle<System<1>>,field:Handle<Field<1>>,skin:Input,dt:Input,seed:u64|{system.with(ctx,|s|field.with_mut(ctx,|f|{s.frame.same(&f.frame)?;let skin=s.frame.scalar(&skin.0,LENGTH,false)?;let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0||skin<0.0{return Err(api::error("time step/skin invalid"));}let ff=f.native.as_ref().ok_or_else(||api::error("force field consumed"))?;field_indices(ff,s.native.n())?;ff.validate().map_err(api::error)?;let list=e::NeighbourList::new(ff.rcut,skin);safe_grid(&s.native,&list)?;let ff=f.native.take().ok_or_else(||api::error("force field already consumed"))?;let native=e::MdSimulation::new(s.native.clone(),ff,skin,dt,seed).map_err(api::error)?;Ok(NewObject::new(Simulation{native,extra_work:f.extra_work,frame:s.frame.clone()},s.frame.dims.to_vec()))}))});
object_binding!(engine_simulation_dof_d1,"num_native_engine_simulation_dof_d1","crate::sim::md_engine::MdSimulation<D>::dof","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->f64,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{Ok(s.native.dof())})});
object_binding!(engine_simulation_temperature_d1,"num_native_engine_simulation_temperature_d1","crate::sim::md_engine::MdSimulation<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{if s.native.dof()<=0.0{return Err(api::error("simulation temperature needs positive unconstrained DOF"));}Ok(Dynamic(s.frame.output(s.native.temperature(),ENERGY,true)?))})});
object_binding!(engine_simulation_total_energy_d1,"num_native_engine_simulation_total_energy_d1","crate::sim::md_engine::MdSimulation<D>::total_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_conserved_energy_d1,"num_native_engine_simulation_conserved_energy_d1","crate::sim::md_engine::MdSimulation<D>::conserved_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.conserved_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_pressure_d1,"num_native_engine_simulation_pressure_d1","crate::sim::md_engine::MdSimulation<D>::pressure","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.pressure(),[1,1,-2,0],false)?))})});
object_binding!(engine_simulation_pressure_tensor_d1,"num_native_engine_simulation_pressure_tensor_d1","crate::sim::md_engine::MdSimulation<D>::pressure_tensor","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|{Ok(Dynamic(packed(&s.frame,s.native.pressure_tensor().to_vec(),[1,1,-2,0])?))})});
object_binding!(engine_simulation_step_d1,"num_native_engine_simulation_step_d1","crate::sim::md_engine::MdSimulation<D>::step","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with_mut(ctx,|s|{let steps=1usize;if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.step().map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_run_d1,"num_native_engine_simulation_run_d1","crate::sim::md_engine::MdSimulation<D>::run","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>,steps:usize)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>,steps:usize|{simulation.with_mut(ctx,|s|{if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.run(steps).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_state_d1,"num_native_engine_simulation_state_d1","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<1>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with(ctx,|s|Ok(Dynamic(system_data(&s.native.sys,&s.frame)?)))});
object_binding!(engine_rdf_new_d1,"num_native_engine_rdf_new_d1","crate::sim::md_engine::Rdf<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,bins:usize,radius:Input)->NewObject<Radial<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,bins:usize,radius:Input|{bound(bins,4096)?;if bins==0{return Err(api::error("RDF needs nonzero bins"));}let f=frame(ctx,&f)?;let rmax=f.scalar(&radius.0,LENGTH,false)?;if rmax<=0.0{return Err(api::error("RDF radius must be positive"));}let dims=f.dims.to_vec();Ok(NewObject::new(Radial{native:e::Rdf::<1>::new(bins,rmax),frame:f,frames:0},dims))});
object_binding!(engine_rdf_add_frame_d1,"num_native_engine_rdf_add_frame_d1","crate::sim::md_engine::Rdf<D>::add_frame","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<1>>,system:Handle<System<1>>)->bool,|ctx:&CallCtx<'_>,rdf:Handle<Radial<1>>,system:Handle<System<1>>|{rdf.with_mut(ctx,|r|system.with(ctx,|s|{r.frame.same(&s.frame)?;if s.native.n()<2||s.native.bx.len.iter().any(|l|2.0*r.native.rmax>*l){return Err(api::error("RDF needs>=2particles and radius<=halfbox"));}
if r.frames>=Radial::<1>::MAX_FRAMES{return Err(api::error("RDF frame-count budget exhausted"));}r.native.add_frame(&s.native);r.frames+=1;Ok(true)}))});
object_binding!(engine_rdf_result_d1,"num_native_engine_rdf_result_d1","crate::sim::md_engine::Rdf<D>::result","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<1>>)->Dynamic,|ctx:&CallCtx<'_>,rdf:Handle<Radial<1>>|{rdf.with(ctx,|r|{if r.frames==0{return Err(api::error("RDF needs accumulated frames"));}let(x,y):(Vec<_>,Vec<_>)=r.native.result().into_iter().unzip();Ok(Dynamic(api::record(vec![("radius",r.frame.array(x,LENGTH)?),("density_ratio",Output::encode(y)?)])))})});
object_binding!(engine_msd_d1,"num_native_engine_msd_d1","crate::sim::md_engine::msd","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<1>>,reference:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<1>>,reference:Input|{system.with(ctx,|s|{let r=rows::<1>(&s.frame,&reference.0,LENGTH)?;if r.len()!=s.native.n(){return Err(api::error("MSD reference size differs from system"));}Ok(Dynamic(s.frame.output(e::msd(&s.native,&r),[0,2,0,0],false)?))})});
object_binding!(engine_vacf_d1,"num_native_engine_vacf_d1","crate::sim::md_engine::vacf","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,frames:Input,lags:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,frames:Input,lags:usize|{let f=frame(ctx,&f)?;let Value::Tensor(t)=&frames.0 else{return Err(api::error("VACF frame tensor expected"));};
if t.shape.len()!=3||t.shape[2]!=1||t.shape[0]==0||t.shape[1]==0||lags==0||lags>t.shape[0]||t.data.len().checked_mul(lags).is_none_or(|n|n>10_000_000){return Err(api::error("VACF shape/lag/work invalid"));}api::same_dimension(&api::dimension(&frames.0)?,&f.dim(VELOCITY)?)?;let points=t.data.chunks(1).map(|r|r.iter().map(|v|f.apply(crate::rules::numerical_abi::numeric(v,0)?,VELOCITY,false,false)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("VACF row"))).collect::<Result<Vec<[f64;1]>,OpError>>()?;let rows=points.chunks(t.shape[1]).map(<[[f64;1]]>::to_vec).collect::<Vec<_>>();Ok(Dynamic(f.array(e::vacf(&rows,lags),[0,2,-2,0])?))});
object_binding!(engine_box_periodic_d2,"num_native_engine_box_periodic_d2","crate::sim::md_engine::SimBox<D>::periodic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input)->NewObject<BoxState<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input|{let f=frame(ctx,&f)?;let lengths=array::<2>(&f,&lengths.0,LENGTH)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::periodic(lengths),frame:f},dims))});
object_binding!(engine_box_new_d2,"num_native_engine_box_new_d2","crate::sim::md_engine::SimBox<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input,boundaries:Input)->NewObject<BoxState<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input,boundaries:Input|{let f=frame(ctx,&f)?;let lengths=array::<2>(&f,&lengths.0,LENGTH)?;let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("boundaries require a rank1 vector"));};
if t.shape!=[2]{return Err(api::error("boundary count must equalD"));}let bc=t.data.iter().map(|v|boundary(&f,v)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("boundary array mismatch"))?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::new(lengths,bc),frame:f},dims))});
object_binding!(engine_box_triclinic_d2,"num_native_engine_box_triclinic_d2","crate::sim::md_engine::SimBox<D>::triclinic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cell:Input)->NewObject<BoxState<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cell:Input|{let f=frame(ctx,&f)?;let h=matrix::<2>(&f,&cell.0,LENGTH)?;let native=e::SimBox::triclinic(h).map_err(api::error)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native,frame:f},dims))});
object_binding!(engine_box_set_cell_d2,"num_native_engine_box_set_cell_d2","crate::sim::md_engine::SimBox<D>::set_cell","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,cell:Input)->bool,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,cell:Input|{b.with_mut(ctx,|b|{let h=matrix::<2>(&b.frame,&cell.0,LENGTH)?;b.native.set_cell(h).map_err(api::error)?;b.validate()?;Ok(true)})});
object_binding!(engine_box_matrix_d2,"num_native_engine_box_matrix_d2","crate::sim::md_engine::SimBox<D>::matrix","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>|{b.with(ctx,|b|Ok(Dynamic(packed(&b.frame,b.native.matrix().to_vec(),LENGTH)?)))});
object_binding!(engine_box_volume_d2,"num_native_engine_box_volume_d2","crate::sim::md_engine::SimBox<D>::volume","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>|{b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.volume(),[0,2,0,0],false)?)))});
object_binding!(engine_box_height_d2,"num_native_engine_box_height_d2","crate::sim::md_engine::SimBox<D>::height","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,axis:usize)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,axis:usize|{if axis>=2{return Err(api::error("box axis outsideD"));}b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.height(axis),LENGTH,false)?)))});
object_binding!(engine_box_to_fractional_d2,"num_native_engine_box_to_fractional_d2","crate::sim::md_engine::SimBox<D>::to_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,point:Input)->[f64;2],|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,point:Input|{b.with(ctx,|b|Ok(b.native.to_fractional(&array::<2>(&b.frame,&point.0,LENGTH)?)))});
object_binding!(engine_box_from_fractional_d2,"num_native_engine_box_from_fractional_d2","crate::sim::md_engine::SimBox<D>::from_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,point:[f64;2])->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,point:[f64;2]|{b.with(ctx,|b|Ok(Dynamic(b.frame.array(b.native.from_fractional(&point).to_vec(),LENGTH)?)))});
object_binding!(engine_box_min_image_d2,"num_native_engine_box_min_image_d2","crate::sim::md_engine::SimBox<D>::min_image","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,displacement:Input)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,displacement:Input|{b.with(ctx,|b|{let mut x=array::<2>(&b.frame,&displacement.0,LENGTH)?;b.native.min_image(&mut x);Ok(Dynamic(b.frame.array(x.to_vec(),LENGTH)?))})});
object_binding!(engine_system_new_d2,"num_native_engine_system_new_d2","crate::sim::md_engine::MdSystem<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,positions:Input)->NewObject<System<2>>,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,positions:Input|{b.with(ctx,|b|{let positions=rows::<2>(&b.frame,&positions.0,LENGTH)?;if positions.is_empty(){return Err(api::error("system needs at least one particle"));}Ok(NewObject::new(System{native:e::MdSystem::new(positions,b.native.clone()),frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(engine_system_data_d2,"num_native_engine_system_data_d2","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<2>>|{system.with(ctx,|s|Ok(Dynamic(system_data(&s.native,&s.frame)?)))});
object_binding!(engine_system_set_state_d2,"num_native_engine_system_set_state_d2","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,velocities:Input,masses:Input,charges:Input)->bool,|ctx:&CallCtx<'_>,system:Handle<System<2>>,velocities:Input,masses:Input,charges:Input|{system.with_mut(ctx,|s|{let vel=rows::<2>(&s.frame,&velocities.0,VELOCITY)?;let mass=s.frame.vector(&masses.0,MASS)?;let charge=s.frame.vector(&charges.0,CHARGE)?;if vel.len()!=s.native.n()||mass.len()!=s.native.n()||charge.len()!=s.native.n()||mass.iter().any(|x|*x<=0.0){return Err(api::error("system state field shape/mass mismatch"));}s.native.vel=vel;s.native.mass=mass;s.native.charge=charge;s.validate()?;Ok(true)})});
object_binding!(engine_system_n_d2,"num_native_engine_system_n_d2","crate::sim::md_engine::MdSystem<D>::n","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>)->usize,|ctx:&CallCtx<'_>,system:Handle<System<2>>|{system.with(ctx,|s|Ok(s.native.n()))});
object_binding!(engine_system_kinetic_energy_d2,"num_native_engine_system_kinetic_energy_d2","crate::sim::md_engine::MdSystem<D>::kinetic_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<2>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.kinetic_energy(),ENERGY,false)?)))});
object_binding!(engine_system_momentum_d2,"num_native_engine_system_momentum_d2","crate::sim::md_engine::MdSystem<D>::momentum","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<2>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?)))});
object_binding!(engine_system_temperature_d2,"num_native_engine_system_temperature_d2","crate::sim::md_engine::MdSystem<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,dof:f64)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<2>>,dof:f64|{if dof<=0.0{return Err(api::error("thermal degrees must be positive"));}system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.temperature(dof),ENERGY,true)?)))});
object_binding!(engine_system_remove_com_velocity_d2,"num_native_engine_system_remove_com_velocity_d2","crate::sim::md_engine::MdSystem<D>::remove_com_velocity","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<2>>|{system.with_mut(ctx,|s|{s.native.remove_com_velocity();s.validate()?;Ok(true)})});
object_binding!(engine_system_maxwell_boltzmann_d2,"num_native_engine_system_maxwell_boltzmann_d2","crate::sim::md_engine::MdSystem<D>::maxwell_boltzmann","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,temperature:Input,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<2>>,temperature:Input,rng:Handle<Rng>|{system.with_mut(ctx,|s|{let t=s.frame.scalar(&temperature.0,ENERGY,true)?;if t<0.0{return Err(api::error("temperature must be nonnegative"));}rng.with_mut(ctx,|r|{s.native.maxwell_boltzmann(t,r);s.validate()?;Ok(true)})})});
object_binding!(engine_system_apply_boundaries_d2,"num_native_engine_system_apply_boundaries_d2","crate::sim::md_engine::MdSystem<D>::apply_boundaries","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<2>>,rng:Handle<Rng>|{system.with_mut(ctx,|s|{if s.native.pos.iter().flatten().any(|x|x.abs()>1e9*s.native.bx.len.iter().copied().fold(1.0,f64::min)){return Err(api::error("boundary displacement exceeds traversal work bound"));}rng.with_mut(ctx,|r|{s.native.apply_boundaries(r);s.validate()?;Ok(true)})})});
object_binding!(engine_lattice_cubic_d2,"num_native_engine_lattice_cubic_d2","crate::sim::md_engine::lattice_cubic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,counts:[usize;2],spacing:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,counts:[usize;2],spacing:Input|{let f=frame(ctx,&f)?;if counts.contains(&0)||counts.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).is_none_or(|n|n>512){return Err(api::error("lattice particle allocation bound is512"));}let a=f.scalar(&spacing.0,LENGTH,false)?;if a<=0.0{return Err(api::error("lattice spacing must be positive"));}Ok(Dynamic(packed(&f,e::lattice_cubic(counts,a),LENGTH)?))});
object_binding!(engine_random_positions_d2,"num_native_engine_random_positions_d2","crate::sim::md_engine::random_positions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<2>>,n:usize,minimum_distance:Input,rng:Handle<Rng>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<2>>,n:usize,minimum_distance:Input,rng:Handle<Rng>|{bound(n,128)?;b.with(ctx,|b|{let d=b.frame.scalar(&minimum_distance.0,LENGTH,false)?;if d<0.0{return Err(api::error("minimum separation must be nonnegative"));}let out=rng.with_mut(ctx,|r|Ok(e::random_positions(n,&b.native,d,r)))?;Ok(Dynamic(match out{None=>Value::Variant(crate::kernel::value::VariantValue{tag:"None".into(),payload:vec![]}),Some(v)=>Value::Variant(crate::kernel::value::VariantValue{tag:"Some".into(),payload:vec![packed(&b.frame,v,LENGTH)?]})}))})});
object_binding!(engine_forcefield_new_d2,"num_native_engine_forcefield_new_d2","crate::sim::md_engine::ForceField<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,terms:Input,cutoff:Input)->NewObject<Field<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,terms:Input,cutoff:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;if rcut<=0.0{return Err(api::error("force cutoff must be positive"));}let Value::Tensor(t)=&terms.0 else{return Err(api::error("pair terms must be rank1handles"));};
if t.shape!=[t.data.len()]||t.data.len()>64{return Err(api::error("pair term shape/count bound exceeded"));}let pair=t.data.iter().map(|v|{let h=<Handle<Pair> as NativeArg>::decode(v,0,ctx)?;h.with(ctx,|p|{f.same(&p.frame)?;Ok(p.native)})}).collect::<Result<_,_>>()?;let dims=f.dims.to_vec();Ok(NewObject::new(Field{native:Some(e::ForceField::new(pair,rcut)),extra_work:0,frame:f},dims))});
object_binding!(engine_forcefield_validate_d2,"num_native_engine_forcefield_validate_d2","crate::sim::md_engine::ForceField<D>::validate","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<2>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<2>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.validate().map_err(api::error)?;Ok(true)})});
object_binding!(engine_forcefield_rebuild_exclusions_d2,"num_native_engine_forcefield_rebuild_exclusions_d2","crate::sim::md_engine::ForceField<D>::rebuild_exclusions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<2>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<2>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.rebuild_exclusions();Ok(true)})});
object_binding!(engine_forces_zeros_d2,"num_native_engine_forces_zeros_d2","crate::sim::md_engine::Forces<D>::zeros","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,n:usize)->NewObject<ForceState<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,n:usize|{bound(n,512)?;let f=frame(ctx,&f)?;let dims=f.dims.to_vec();Ok(NewObject::new(ForceState{native:e::Forces::<2>::zeros(n),frame:f},dims))});
object_binding!(engine_forces_data_d2,"num_native_engine_forces_data_d2","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(forces:Handle<ForceState<2>>)->Dynamic,|ctx:&CallCtx<'_>,forces:Handle<ForceState<2>>|{forces.with(ctx,|s|Ok(Dynamic(force_data(&s.native,&s.frame)?)))});
object_binding!(engine_neighbours_new_d2,"num_native_engine_neighbours_new_d2","crate::sim::md_engine::NeighbourList<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cutoff:Input,skin:Input)->NewObject<List<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cutoff:Input,skin:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;let skin=f.scalar(&skin.0,LENGTH,false)?;if rcut<=0.0||skin<0.0{return Err(api::error("cutoff must bepositive andskin nonnegative"));}let dims=f.dims.to_vec();Ok(NewObject::new(List{native:e::NeighbourList::<2>::new(rcut,skin),frame:f},dims))});
object_binding!(engine_neighbours_needs_rebuild_d2,"num_native_engine_neighbours_needs_rebuild_d2","crate::sim::md_engine::NeighbourList<D>::needs_rebuild","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<2>>,system:Handle<System<2>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<2>>,system:Handle<System<2>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;Ok(l.native.needs_rebuild(&s.native))}))});
object_binding!(engine_neighbours_build_d2,"num_native_engine_neighbours_build_d2","crate::sim::md_engine::NeighbourList<D>::build","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<2>>,system:Handle<System<2>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<2>>,system:Handle<System<2>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;l.native.build(&s.native);l.validate()?;Ok(true)}))});
object_binding!(engine_compute_forces_d2,"num_native_engine_compute_forces_d2","crate::sim::md_engine::compute_forces","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,field:Handle<Field<2>>,list:Handle<List<2>>)->NewObject<ForceState<2>>,|ctx:&CallCtx<'_>,system:Handle<System<2>>,field:Handle<Field<2>>,list:Handle<List<2>>|{system.with(ctx,|s|field.with(ctx,|ff|list.with_mut(ctx,|l|{s.frame.same(&ff.frame)?;s.frame.same(&l.frame)?;safe_grid(&s.native,&l.native)?;field_indices(ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,s.native.n())?;let native=e::compute_forces(&s.native,ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,&mut l.native).map_err(api::error)?;Ok(NewObject::new(ForceState{native,frame:s.frame.clone()},s.frame.dims.to_vec()))})))});
object_binding!(engine_simulation_new_d2,"num_native_engine_simulation_new_d2","crate::sim::md_engine::MdSimulation<D>::new","Construct the actual native driver with an owned force field transferred once. The system is cloned; the force-field handle becomes consumed. Preflight checks run before the transfer; physical frame witnesses and compiled extra-force ownership stay in the simulation.",false,(system:Handle<System<2>>,field:Handle<Field<2>>,skin:Input,dt:Input,seed:u64)->NewObject<Simulation<2>>,|ctx:&CallCtx<'_>,system:Handle<System<2>>,field:Handle<Field<2>>,skin:Input,dt:Input,seed:u64|{system.with(ctx,|s|field.with_mut(ctx,|f|{s.frame.same(&f.frame)?;let skin=s.frame.scalar(&skin.0,LENGTH,false)?;let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0||skin<0.0{return Err(api::error("time step/skin invalid"));}let ff=f.native.as_ref().ok_or_else(||api::error("force field consumed"))?;field_indices(ff,s.native.n())?;ff.validate().map_err(api::error)?;let list=e::NeighbourList::new(ff.rcut,skin);safe_grid(&s.native,&list)?;let ff=f.native.take().ok_or_else(||api::error("force field already consumed"))?;let native=e::MdSimulation::new(s.native.clone(),ff,skin,dt,seed).map_err(api::error)?;Ok(NewObject::new(Simulation{native,extra_work:f.extra_work,frame:s.frame.clone()},s.frame.dims.to_vec()))}))});
object_binding!(engine_simulation_dof_d2,"num_native_engine_simulation_dof_d2","crate::sim::md_engine::MdSimulation<D>::dof","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->f64,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{Ok(s.native.dof())})});
object_binding!(engine_simulation_temperature_d2,"num_native_engine_simulation_temperature_d2","crate::sim::md_engine::MdSimulation<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{if s.native.dof()<=0.0{return Err(api::error("simulation temperature needs positive unconstrained DOF"));}Ok(Dynamic(s.frame.output(s.native.temperature(),ENERGY,true)?))})});
object_binding!(engine_simulation_total_energy_d2,"num_native_engine_simulation_total_energy_d2","crate::sim::md_engine::MdSimulation<D>::total_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_conserved_energy_d2,"num_native_engine_simulation_conserved_energy_d2","crate::sim::md_engine::MdSimulation<D>::conserved_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.conserved_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_pressure_d2,"num_native_engine_simulation_pressure_d2","crate::sim::md_engine::MdSimulation<D>::pressure","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.pressure(),[1,0,-2,0],false)?))})});
object_binding!(engine_simulation_pressure_tensor_d2,"num_native_engine_simulation_pressure_tensor_d2","crate::sim::md_engine::MdSimulation<D>::pressure_tensor","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|{Ok(Dynamic(packed(&s.frame,s.native.pressure_tensor().to_vec(),[1,0,-2,0])?))})});
object_binding!(engine_simulation_step_d2,"num_native_engine_simulation_step_d2","crate::sim::md_engine::MdSimulation<D>::step","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with_mut(ctx,|s|{let steps=1usize;if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.step().map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_run_d2,"num_native_engine_simulation_run_d2","crate::sim::md_engine::MdSimulation<D>::run","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>,steps:usize)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>,steps:usize|{simulation.with_mut(ctx,|s|{if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.run(steps).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_state_d2,"num_native_engine_simulation_state_d2","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<2>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with(ctx,|s|Ok(Dynamic(system_data(&s.native.sys,&s.frame)?)))});
object_binding!(engine_rdf_new_d2,"num_native_engine_rdf_new_d2","crate::sim::md_engine::Rdf<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,bins:usize,radius:Input)->NewObject<Radial<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,bins:usize,radius:Input|{bound(bins,4096)?;if bins==0{return Err(api::error("RDF needs nonzero bins"));}let f=frame(ctx,&f)?;let rmax=f.scalar(&radius.0,LENGTH,false)?;if rmax<=0.0{return Err(api::error("RDF radius must be positive"));}let dims=f.dims.to_vec();Ok(NewObject::new(Radial{native:e::Rdf::<2>::new(bins,rmax),frame:f,frames:0},dims))});
object_binding!(engine_rdf_add_frame_d2,"num_native_engine_rdf_add_frame_d2","crate::sim::md_engine::Rdf<D>::add_frame","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<2>>,system:Handle<System<2>>)->bool,|ctx:&CallCtx<'_>,rdf:Handle<Radial<2>>,system:Handle<System<2>>|{rdf.with_mut(ctx,|r|system.with(ctx,|s|{r.frame.same(&s.frame)?;if s.native.n()<2||s.native.bx.len.iter().any(|l|2.0*r.native.rmax>*l){return Err(api::error("RDF needs>=2particles and radius<=halfbox"));}
if r.frames>=Radial::<2>::MAX_FRAMES{return Err(api::error("RDF frame-count budget exhausted"));}r.native.add_frame(&s.native);r.frames+=1;Ok(true)}))});
object_binding!(engine_rdf_result_d2,"num_native_engine_rdf_result_d2","crate::sim::md_engine::Rdf<D>::result","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<2>>)->Dynamic,|ctx:&CallCtx<'_>,rdf:Handle<Radial<2>>|{rdf.with(ctx,|r|{if r.frames==0{return Err(api::error("RDF needs accumulated frames"));}let(x,y):(Vec<_>,Vec<_>)=r.native.result().into_iter().unzip();Ok(Dynamic(api::record(vec![("radius",r.frame.array(x,LENGTH)?),("density_ratio",Output::encode(y)?)])))})});
object_binding!(engine_msd_d2,"num_native_engine_msd_d2","crate::sim::md_engine::msd","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<2>>,reference:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<2>>,reference:Input|{system.with(ctx,|s|{let r=rows::<2>(&s.frame,&reference.0,LENGTH)?;if r.len()!=s.native.n(){return Err(api::error("MSD reference size differs from system"));}Ok(Dynamic(s.frame.output(e::msd(&s.native,&r),[0,2,0,0],false)?))})});
object_binding!(engine_vacf_d2,"num_native_engine_vacf_d2","crate::sim::md_engine::vacf","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,frames:Input,lags:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,frames:Input,lags:usize|{let f=frame(ctx,&f)?;let Value::Tensor(t)=&frames.0 else{return Err(api::error("VACF frame tensor expected"));};
if t.shape.len()!=3||t.shape[2]!=2||t.shape[0]==0||t.shape[1]==0||lags==0||lags>t.shape[0]||t.data.len().checked_mul(lags).is_none_or(|n|n>10_000_000){return Err(api::error("VACF shape/lag/work invalid"));}api::same_dimension(&api::dimension(&frames.0)?,&f.dim(VELOCITY)?)?;let points=t.data.chunks(2).map(|r|r.iter().map(|v|f.apply(crate::rules::numerical_abi::numeric(v,0)?,VELOCITY,false,false)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("VACF row"))).collect::<Result<Vec<[f64;2]>,OpError>>()?;let rows=points.chunks(t.shape[1]).map(<[[f64;2]]>::to_vec).collect::<Vec<_>>();Ok(Dynamic(f.array(e::vacf(&rows,lags),[0,2,-2,0])?))});
object_binding!(engine_box_periodic_d3,"num_native_engine_box_periodic_d3","crate::sim::md_engine::SimBox<D>::periodic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input)->NewObject<BoxState<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input|{let f=frame(ctx,&f)?;let lengths=array::<3>(&f,&lengths.0,LENGTH)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::periodic(lengths),frame:f},dims))});
object_binding!(engine_box_new_d3,"num_native_engine_box_new_d3","crate::sim::md_engine::SimBox<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,lengths:Input,boundaries:Input)->NewObject<BoxState<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,lengths:Input,boundaries:Input|{let f=frame(ctx,&f)?;let lengths=array::<3>(&f,&lengths.0,LENGTH)?;let Value::Tensor(t)=&boundaries.0 else{return Err(api::error("boundaries require a rank1 vector"));};
if t.shape!=[3]{return Err(api::error("boundary count must equalD"));}let bc=t.data.iter().map(|v|boundary(&f,v)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("boundary array mismatch"))?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native:e::SimBox::new(lengths,bc),frame:f},dims))});
object_binding!(engine_box_triclinic_d3,"num_native_engine_box_triclinic_d3","crate::sim::md_engine::SimBox<D>::triclinic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cell:Input)->NewObject<BoxState<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cell:Input|{let f=frame(ctx,&f)?;let h=matrix::<3>(&f,&cell.0,LENGTH)?;let native=e::SimBox::triclinic(h).map_err(api::error)?;let dims=f.dims.to_vec();Ok(NewObject::new(BoxState{native,frame:f},dims))});
object_binding!(engine_box_set_cell_d3,"num_native_engine_box_set_cell_d3","crate::sim::md_engine::SimBox<D>::set_cell","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,cell:Input)->bool,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,cell:Input|{b.with_mut(ctx,|b|{let h=matrix::<3>(&b.frame,&cell.0,LENGTH)?;b.native.set_cell(h).map_err(api::error)?;b.validate()?;Ok(true)})});
object_binding!(engine_box_matrix_d3,"num_native_engine_box_matrix_d3","crate::sim::md_engine::SimBox<D>::matrix","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>|{b.with(ctx,|b|Ok(Dynamic(packed(&b.frame,b.native.matrix().to_vec(),LENGTH)?)))});
object_binding!(engine_box_volume_d3,"num_native_engine_box_volume_d3","crate::sim::md_engine::SimBox<D>::volume","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>|{b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.volume(),[0,3,0,0],false)?)))});
object_binding!(engine_box_height_d3,"num_native_engine_box_height_d3","crate::sim::md_engine::SimBox<D>::height","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,axis:usize)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,axis:usize|{if axis>=3{return Err(api::error("box axis outsideD"));}b.with(ctx,|b|Ok(Dynamic(b.frame.output(b.native.height(axis),LENGTH,false)?)))});
object_binding!(engine_box_to_fractional_d3,"num_native_engine_box_to_fractional_d3","crate::sim::md_engine::SimBox<D>::to_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,point:Input)->[f64;3],|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,point:Input|{b.with(ctx,|b|Ok(b.native.to_fractional(&array::<3>(&b.frame,&point.0,LENGTH)?)))});
object_binding!(engine_box_from_fractional_d3,"num_native_engine_box_from_fractional_d3","crate::sim::md_engine::SimBox<D>::from_fractional","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,point:[f64;3])->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,point:[f64;3]|{b.with(ctx,|b|Ok(Dynamic(b.frame.array(b.native.from_fractional(&point).to_vec(),LENGTH)?)))});
object_binding!(engine_box_min_image_d3,"num_native_engine_box_min_image_d3","crate::sim::md_engine::SimBox<D>::min_image","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,displacement:Input)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,displacement:Input|{b.with(ctx,|b|{let mut x=array::<3>(&b.frame,&displacement.0,LENGTH)?;b.native.min_image(&mut x);Ok(Dynamic(b.frame.array(x.to_vec(),LENGTH)?))})});
object_binding!(engine_system_new_d3,"num_native_engine_system_new_d3","crate::sim::md_engine::MdSystem<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,positions:Input)->NewObject<System<3>>,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,positions:Input|{b.with(ctx,|b|{let positions=rows::<3>(&b.frame,&positions.0,LENGTH)?;if positions.is_empty(){return Err(api::error("system needs at least one particle"));}Ok(NewObject::new(System{native:e::MdSystem::new(positions,b.native.clone()),frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(engine_system_data_d3,"num_native_engine_system_data_d3","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>|{system.with(ctx,|s|Ok(Dynamic(system_data(&s.native,&s.frame)?)))});
object_binding!(engine_system_set_state_d3,"num_native_engine_system_set_state_d3","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,velocities:Input,masses:Input,charges:Input)->bool,|ctx:&CallCtx<'_>,system:Handle<System<3>>,velocities:Input,masses:Input,charges:Input|{system.with_mut(ctx,|s|{let vel=rows::<3>(&s.frame,&velocities.0,VELOCITY)?;let mass=s.frame.vector(&masses.0,MASS)?;let charge=s.frame.vector(&charges.0,CHARGE)?;if vel.len()!=s.native.n()||mass.len()!=s.native.n()||charge.len()!=s.native.n()||mass.iter().any(|x|*x<=0.0){return Err(api::error("system state field shape/mass mismatch"));}s.native.vel=vel;s.native.mass=mass;s.native.charge=charge;s.validate()?;Ok(true)})});
object_binding!(engine_system_n_d3,"num_native_engine_system_n_d3","crate::sim::md_engine::MdSystem<D>::n","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>)->usize,|ctx:&CallCtx<'_>,system:Handle<System<3>>|{system.with(ctx,|s|Ok(s.native.n()))});
object_binding!(engine_system_kinetic_energy_d3,"num_native_engine_system_kinetic_energy_d3","crate::sim::md_engine::MdSystem<D>::kinetic_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.kinetic_energy(),ENERGY,false)?)))});
object_binding!(engine_system_momentum_d3,"num_native_engine_system_momentum_d3","crate::sim::md_engine::MdSystem<D>::momentum","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>|{system.with(ctx,|s|Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?)))});
object_binding!(engine_system_temperature_d3,"num_native_engine_system_temperature_d3","crate::sim::md_engine::MdSystem<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,dof:f64)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,dof:f64|{if dof<=0.0{return Err(api::error("thermal degrees must be positive"));}system.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.temperature(dof),ENERGY,true)?)))});
object_binding!(engine_system_remove_com_velocity_d3,"num_native_engine_system_remove_com_velocity_d3","crate::sim::md_engine::MdSystem<D>::remove_com_velocity","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<3>>|{system.with_mut(ctx,|s|{s.native.remove_com_velocity();s.validate()?;Ok(true)})});
object_binding!(engine_system_maxwell_boltzmann_d3,"num_native_engine_system_maxwell_boltzmann_d3","crate::sim::md_engine::MdSystem<D>::maxwell_boltzmann","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,temperature:Input,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<3>>,temperature:Input,rng:Handle<Rng>|{system.with_mut(ctx,|s|{let t=s.frame.scalar(&temperature.0,ENERGY,true)?;if t<0.0{return Err(api::error("temperature must be nonnegative"));}rng.with_mut(ctx,|r|{s.native.maxwell_boltzmann(t,r);s.validate()?;Ok(true)})})});
object_binding!(engine_system_apply_boundaries_d3,"num_native_engine_system_apply_boundaries_d3","crate::sim::md_engine::MdSystem<D>::apply_boundaries","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,system:Handle<System<3>>,rng:Handle<Rng>|{system.with_mut(ctx,|s|{if s.native.pos.iter().flatten().any(|x|x.abs()>1e9*s.native.bx.len.iter().copied().fold(1.0,f64::min)){return Err(api::error("boundary displacement exceeds traversal work bound"));}rng.with_mut(ctx,|r|{s.native.apply_boundaries(r);s.validate()?;Ok(true)})})});
object_binding!(engine_lattice_cubic_d3,"num_native_engine_lattice_cubic_d3","crate::sim::md_engine::lattice_cubic","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,counts:[usize;3],spacing:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,counts:[usize;3],spacing:Input|{let f=frame(ctx,&f)?;if counts.contains(&0)||counts.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).is_none_or(|n|n>512){return Err(api::error("lattice particle allocation bound is512"));}let a=f.scalar(&spacing.0,LENGTH,false)?;if a<=0.0{return Err(api::error("lattice spacing must be positive"));}Ok(Dynamic(packed(&f,e::lattice_cubic(counts,a),LENGTH)?))});
object_binding!(engine_random_positions_d3,"num_native_engine_random_positions_d3","crate::sim::md_engine::random_positions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(b:Handle<BoxState<3>>,n:usize,minimum_distance:Input,rng:Handle<Rng>)->Dynamic,|ctx:&CallCtx<'_>,b:Handle<BoxState<3>>,n:usize,minimum_distance:Input,rng:Handle<Rng>|{bound(n,128)?;b.with(ctx,|b|{let d=b.frame.scalar(&minimum_distance.0,LENGTH,false)?;if d<0.0{return Err(api::error("minimum separation must be nonnegative"));}let out=rng.with_mut(ctx,|r|Ok(e::random_positions(n,&b.native,d,r)))?;Ok(Dynamic(match out{None=>Value::Variant(crate::kernel::value::VariantValue{tag:"None".into(),payload:vec![]}),Some(v)=>Value::Variant(crate::kernel::value::VariantValue{tag:"Some".into(),payload:vec![packed(&b.frame,v,LENGTH)?]})}))})});
object_binding!(engine_forcefield_new_d3,"num_native_engine_forcefield_new_d3","crate::sim::md_engine::ForceField<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,terms:Input,cutoff:Input)->NewObject<Field<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,terms:Input,cutoff:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;if rcut<=0.0{return Err(api::error("force cutoff must be positive"));}let Value::Tensor(t)=&terms.0 else{return Err(api::error("pair terms must be rank1handles"));};
if t.shape!=[t.data.len()]||t.data.len()>64{return Err(api::error("pair term shape/count bound exceeded"));}let pair=t.data.iter().map(|v|{let h=<Handle<Pair> as NativeArg>::decode(v,0,ctx)?;h.with(ctx,|p|{f.same(&p.frame)?;Ok(p.native)})}).collect::<Result<_,_>>()?;let dims=f.dims.to_vec();Ok(NewObject::new(Field{native:Some(e::ForceField::new(pair,rcut)),extra_work:0,frame:f},dims))});
object_binding!(engine_forcefield_validate_d3,"num_native_engine_forcefield_validate_d3","crate::sim::md_engine::ForceField<D>::validate","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<3>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<3>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.validate().map_err(api::error)?;Ok(true)})});
object_binding!(engine_forcefield_rebuild_exclusions_d3,"num_native_engine_forcefield_rebuild_exclusions_d3","crate::sim::md_engine::ForceField<D>::rebuild_exclusions","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(field:Handle<Field<3>>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<3>>|{field.with_mut(ctx,|s|{let ff=s.native.as_mut().ok_or_else(||api::error("force field consumed"))?;ff.rebuild_exclusions();Ok(true)})});
object_binding!(engine_forces_zeros_d3,"num_native_engine_forces_zeros_d3","crate::sim::md_engine::Forces<D>::zeros","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,n:usize)->NewObject<ForceState<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,n:usize|{bound(n,512)?;let f=frame(ctx,&f)?;let dims=f.dims.to_vec();Ok(NewObject::new(ForceState{native:e::Forces::<3>::zeros(n),frame:f},dims))});
object_binding!(engine_forces_data_d3,"num_native_engine_forces_data_d3","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(forces:Handle<ForceState<3>>)->Dynamic,|ctx:&CallCtx<'_>,forces:Handle<ForceState<3>>|{forces.with(ctx,|s|Ok(Dynamic(force_data(&s.native,&s.frame)?)))});
object_binding!(engine_neighbours_new_d3,"num_native_engine_neighbours_new_d3","crate::sim::md_engine::NeighbourList<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,cutoff:Input,skin:Input)->NewObject<List<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,cutoff:Input,skin:Input|{let f=frame(ctx,&f)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;let skin=f.scalar(&skin.0,LENGTH,false)?;if rcut<=0.0||skin<0.0{return Err(api::error("cutoff must bepositive andskin nonnegative"));}let dims=f.dims.to_vec();Ok(NewObject::new(List{native:e::NeighbourList::<3>::new(rcut,skin),frame:f},dims))});
object_binding!(engine_neighbours_needs_rebuild_d3,"num_native_engine_neighbours_needs_rebuild_d3","crate::sim::md_engine::NeighbourList<D>::needs_rebuild","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<3>>,system:Handle<System<3>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<3>>,system:Handle<System<3>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;Ok(l.native.needs_rebuild(&s.native))}))});
object_binding!(engine_neighbours_build_d3,"num_native_engine_neighbours_build_d3","crate::sim::md_engine::NeighbourList<D>::build","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(list:Handle<List<3>>,system:Handle<System<3>>)->bool,|ctx:&CallCtx<'_>,list:Handle<List<3>>,system:Handle<System<3>>|{list.with_mut(ctx,|l|system.with(ctx,|s|{l.frame.same(&s.frame)?;safe_grid(&s.native,&l.native)?;l.native.build(&s.native);l.validate()?;Ok(true)}))});
object_binding!(engine_compute_forces_d3,"num_native_engine_compute_forces_d3","crate::sim::md_engine::compute_forces","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,field:Handle<Field<3>>,list:Handle<List<3>>)->NewObject<ForceState<3>>,|ctx:&CallCtx<'_>,system:Handle<System<3>>,field:Handle<Field<3>>,list:Handle<List<3>>|{system.with(ctx,|s|field.with(ctx,|ff|list.with_mut(ctx,|l|{s.frame.same(&ff.frame)?;s.frame.same(&l.frame)?;safe_grid(&s.native,&l.native)?;field_indices(ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,s.native.n())?;let native=e::compute_forces(&s.native,ff.native.as_ref().ok_or_else(||api::error("force field consumed"))?,&mut l.native).map_err(api::error)?;Ok(NewObject::new(ForceState{native,frame:s.frame.clone()},s.frame.dims.to_vec()))})))});
object_binding!(engine_simulation_new_d3,"num_native_engine_simulation_new_d3","crate::sim::md_engine::MdSimulation<D>::new","Construct the actual native driver with an owned force field transferred once. The system is cloned; the force-field handle becomes consumed. Preflight checks run before the transfer; physical frame witnesses and compiled extra-force ownership stay in the simulation.",false,(system:Handle<System<3>>,field:Handle<Field<3>>,skin:Input,dt:Input,seed:u64)->NewObject<Simulation<3>>,|ctx:&CallCtx<'_>,system:Handle<System<3>>,field:Handle<Field<3>>,skin:Input,dt:Input,seed:u64|{system.with(ctx,|s|field.with_mut(ctx,|f|{s.frame.same(&f.frame)?;let skin=s.frame.scalar(&skin.0,LENGTH,false)?;let dt=s.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0||skin<0.0{return Err(api::error("time step/skin invalid"));}let ff=f.native.as_ref().ok_or_else(||api::error("force field consumed"))?;field_indices(ff,s.native.n())?;ff.validate().map_err(api::error)?;let list=e::NeighbourList::new(ff.rcut,skin);safe_grid(&s.native,&list)?;let ff=f.native.take().ok_or_else(||api::error("force field already consumed"))?;let native=e::MdSimulation::new(s.native.clone(),ff,skin,dt,seed).map_err(api::error)?;Ok(NewObject::new(Simulation{native,extra_work:f.extra_work,frame:s.frame.clone()},s.frame.dims.to_vec()))}))});
object_binding!(engine_simulation_dof_d3,"num_native_engine_simulation_dof_d3","crate::sim::md_engine::MdSimulation<D>::dof","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->f64,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{Ok(s.native.dof())})});
object_binding!(engine_simulation_temperature_d3,"num_native_engine_simulation_temperature_d3","crate::sim::md_engine::MdSimulation<D>::temperature","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{if s.native.dof()<=0.0{return Err(api::error("simulation temperature needs positive unconstrained DOF"));}Ok(Dynamic(s.frame.output(s.native.temperature(),ENERGY,true)?))})});
object_binding!(engine_simulation_total_energy_d3,"num_native_engine_simulation_total_energy_d3","crate::sim::md_engine::MdSimulation<D>::total_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_conserved_energy_d3,"num_native_engine_simulation_conserved_energy_d3","crate::sim::md_engine::MdSimulation<D>::conserved_energy","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.conserved_energy(),ENERGY,false)?))})});
object_binding!(engine_simulation_pressure_d3,"num_native_engine_simulation_pressure_d3","crate::sim::md_engine::MdSimulation<D>::pressure","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.pressure(),[1,-1,-2,0],false)?))})});
object_binding!(engine_simulation_pressure_tensor_d3,"num_native_engine_simulation_pressure_tensor_d3","crate::sim::md_engine::MdSimulation<D>::pressure_tensor","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|{Ok(Dynamic(packed(&s.frame,s.native.pressure_tensor().to_vec(),[1,-1,-2,0])?))})});
object_binding!(engine_simulation_step_d3,"num_native_engine_simulation_step_d3","crate::sim::md_engine::MdSimulation<D>::step","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with_mut(ctx,|s|{let steps=1usize;if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.step().map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_run_d3,"num_native_engine_simulation_run_d3","crate::sim::md_engine::MdSimulation<D>::run","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>,steps:usize)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>,steps:usize|{simulation.with_mut(ctx,|s|{if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100_000_000){return Err(api::error("MD pair-step work budget exceeded"));}s.native.run(steps).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(engine_simulation_state_d3,"num_native_engine_simulation_state_d3","","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(simulation:Handle<Simulation<3>>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with(ctx,|s|Ok(Dynamic(system_data(&s.native.sys,&s.frame)?)))});
object_binding!(engine_rdf_new_d3,"num_native_engine_rdf_new_d3","crate::sim::md_engine::Rdf<D>::new","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,bins:usize,radius:Input)->NewObject<Radial<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,bins:usize,radius:Input|{bound(bins,4096)?;if bins==0{return Err(api::error("RDF needs nonzero bins"));}let f=frame(ctx,&f)?;let rmax=f.scalar(&radius.0,LENGTH,false)?;if rmax<=0.0{return Err(api::error("RDF radius must be positive"));}let dims=f.dims.to_vec();Ok(NewObject::new(Radial{native:e::Rdf::<3>::new(bins,rmax),frame:f,frames:0},dims))});
object_binding!(engine_rdf_add_frame_d3,"num_native_engine_rdf_add_frame_d3","crate::sim::md_engine::Rdf<D>::add_frame","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<3>>,system:Handle<System<3>>)->bool,|ctx:&CallCtx<'_>,rdf:Handle<Radial<3>>,system:Handle<System<3>>|{rdf.with_mut(ctx,|r|system.with(ctx,|s|{r.frame.same(&s.frame)?;if s.native.n()<2||s.native.bx.len.iter().any(|l|2.0*r.native.rmax>*l){return Err(api::error("RDF needs>=2particles and radius<=halfbox"));}
if r.frames>=Radial::<3>::MAX_FRAMES{return Err(api::error("RDF frame-count budget exhausted"));}r.native.add_frame(&s.native);r.frames+=1;Ok(true)}))});
object_binding!(engine_rdf_result_d3,"num_native_engine_rdf_result_d3","crate::sim::md_engine::Rdf<D>::result","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(rdf:Handle<Radial<3>>)->Dynamic,|ctx:&CallCtx<'_>,rdf:Handle<Radial<3>>|{rdf.with(ctx,|r|{if r.frames==0{return Err(api::error("RDF needs accumulated frames"));}let(x,y):(Vec<_>,Vec<_>)=r.native.result().into_iter().unzip();Ok(Dynamic(api::record(vec![("radius",r.frame.array(x,LENGTH)?),("density_ratio",Output::encode(y)?)])))})});
object_binding!(engine_msd_d3,"num_native_engine_msd_d3","crate::sim::md_engine::msd","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(system:Handle<System<3>>,reference:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,reference:Input|{system.with(ctx,|s|{let r=rows::<3>(&s.frame,&reference.0,LENGTH)?;if r.len()!=s.native.n(){return Err(api::error("MSD reference size differs from system"));}Ok(Dynamic(s.frame.output(e::msd(&s.native,&r),[0,2,0,0],false)?))})});
object_binding!(engine_vacf_d3,"num_native_engine_vacf_d3","crate::sim::md_engine::vacf","Actual native MD method with checked physical or explicitly reduced reference frames. Shapes, finite state, object identity and allocation/work limits precede native work.",false,(f:Handle<Frame>,frames:Input,lags:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,frames:Input,lags:usize|{let f=frame(ctx,&f)?;let Value::Tensor(t)=&frames.0 else{return Err(api::error("VACF frame tensor expected"));};
if t.shape.len()!=3||t.shape[2]!=3||t.shape[0]==0||t.shape[1]==0||lags==0||lags>t.shape[0]||t.data.len().checked_mul(lags).is_none_or(|n|n>10_000_000){return Err(api::error("VACF shape/lag/work invalid"));}api::same_dimension(&api::dimension(&frames.0)?,&f.dim(VELOCITY)?)?;let points=t.data.chunks(3).map(|r|r.iter().map(|v|f.apply(crate::rules::numerical_abi::numeric(v,0)?,VELOCITY,false,false)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|api::error("VACF row"))).collect::<Result<Vec<[f64;3]>,OpError>>()?;let rows=points.chunks(t.shape[1]).map(<[[f64;3]]>::to_vec).collect::<Vec<_>>();Ok(Dynamic(f.array(e::vacf(&rows,lags),[0,2,-2,0])?))});
object_binding!(engine_lattice_fcc,"num_native_engine_lattice_fcc","crate::sim::md_engine::lattice_fcc","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,n:usize,spacing:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,n:usize,spacing:Input|{let f=frame(ctx,&f)?;if n.checked_pow(3).and_then(|v|v.checked_mul(4)).is_none_or(|n|n>512){return Err(api::error("FCC particle bound512"));}let a=f.scalar(&spacing.0,LENGTH,false)?;if a<=0.0{return Err(api::error("lattice spacing must be positive"));}Ok(Dynamic(packed(&f,e::lattice_fcc(n,a),LENGTH)?))});
object_binding!(engine_lattice_triangular,"num_native_engine_lattice_triangular","crate::sim::md_engine::lattice_triangular","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,nx:usize,ny:usize,spacing:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,nx:usize,ny:usize,spacing:Input|{let f=frame(ctx,&f)?;if nx.checked_mul(ny).is_none_or(|n|n>512){return Err(api::error("triangular particle bound512"));}let a=f.scalar(&spacing.0,LENGTH,false)?;if a<=0.0{return Err(api::error("lattice spacing must be positive"));}Ok(Dynamic(packed(&f,e::lattice_triangular(nx,ny,a),LENGTH)?))});
object_binding!(engine_energy_drift_slope,"num_native_engine_energy_drift_slope","crate::sim::md_engine::energy_drift_slope","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,times:Input,energies:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,times:Input,energies:Input|{let f=frame(ctx,&f)?;let t=f.vector(&times.0,TIME)?;let v=f.vector(&energies.0,ENERGY)?;if t.len()!=v.len()||t.len()<2||t.iter().all(|x|(*x-t[0]).abs()<=0.0){return Err(api::error("energy drift requires equal samples and distinct times"));}Ok(Dynamic(f.output(e::energy_drift_slope(&t,&v),[1,2,-3,0],false)?))});
object_binding!(engine_chain_new,"num_native_engine_chain_new","crate::sim::md_engine::NoseHooverChain::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,temperature:Input,tau:Input,length:usize,dof:f64)->NewObject<Chain>,|ctx:&CallCtx<'_>,f:Handle<Frame>,temperature:Input,tau:Input,length:usize,dof:f64|{let f=frame(ctx,&f)?;let t=f.scalar(&temperature.0,ENERGY,true)?;let tau=f.scalar(&tau.0,TIME,false)?;if length==0||length>16||t<=0.0||tau<=0.0||dof<=0.0{return Err(api::error("Nose Hoover chain controls invalid"));}let dims=f.dims.to_vec();Ok(NewObject::new(Chain{native:e::NoseHooverChain::new(t,tau,length,dof),frame:f},dims))});
object_binding!(engine_chain_propagate,"num_native_engine_chain_propagate","crate::sim::md_engine::NoseHooverChain::propagate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(chain:Handle<Chain>,twice_kinetic_energy:Input,step:Input)->f64,|ctx:&CallCtx<'_>,chain:Handle<Chain>,twice_kinetic_energy:Input,step:Input|{chain.with_mut(ctx,|c|{let k=c.frame.scalar(&twice_kinetic_energy.0,ENERGY,false)?;let h=c.frame.scalar(&step.0,TIME,false)?;if k<0.0{return Err(api::error("kinetic energy must be nonnegative"));}let scale=c.native.propagate(k,h);c.validate()?;Ok(scale)})});
object_binding!(engine_chain_energy,"num_native_engine_chain_energy","crate::sim::md_engine::NoseHooverChain::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(chain:Handle<Chain>)->Dynamic,|ctx:&CallCtx<'_>,chain:Handle<Chain>|{chain.with(ctx,|c|Ok(Dynamic(c.frame.output(c.native.energy(),ENERGY,false)?)))});
struct Flexible {
    native: e::FlexibleBarostat,
    frame: Frame,
}
impl NativeType for Flexible {
    const KIND: &'static str = "rssn/native/sim/engine-flexible-barostat/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        if !(1..=3).contains(&self.native.d)
            || self.native.p_ext.len() != self.native.d * self.native.d
            || self.native.pg.len() != self.native.d * self.native.d
            || self.native.w <= 0.0
        {
            return Err(api::error("flexible barostat dimension/state invalid"));
        }
        finite(self.native.p_ext.iter().chain(&self.native.pg).copied())
    }
}
fn cell_mode(value: &str) -> Result<e::CellMode, OpError> {
    match value {
        | "Full" => Ok(e::CellMode::Full),
        | "Orthorhombic" => Ok(e::CellMode::Orthorhombic),
        | "Isotropic" => Ok(e::CellMode::Isotropic),
        | _ => Err(api::error("cell mode must be Full/Orthorhombic/Isotropic")),
    }
}
fn records(value: &Value) -> Result<&[Value], OpError> {
    let Value::Tensor(t) = value else {
        return Err(api::error("record vector expected"));
    };
    if t.shape != [t.data.len()] || t.data.len() > 4096 {
        return Err(api::error("record vector shape/count invalid"));
    }
    Ok(&t.data)
}
fn atom(
    value: &Value,
    name: &str,
    n: usize,
) -> Result<usize, OpError> {
    let index = <usize as Arg>::decode(api::field(value, name)?, 0)?;
    if index >= n {
        return Err(api::error(
            "bonded atom index outside declared particle count",
        ));
    }
    Ok(index)
}
fn set_bonded<const D: usize>(
    f: &Frame,
    ff: &mut e::ForceField<D>,
    v: &Value,
    n: usize,
) -> Result<(), OpError> {
    let mut bonds = vec![];
    let mut angles = vec![];
    let mut dihedrals = vec![];
    for b in records(api::field(v, "bonds")?)? {
        let i = atom(b, "i", n)?;
        let j = atom(b, "j", n)?;
        let k = f.scalar(api::field(b, "k")?, [1, 0, -2, 0], false)?;
        let r0 = f.scalar(api::field(b, "r0")?, LENGTH, false)?;
        if i == j || k < 0.0 || r0 <= 0.0 {
            return Err(api::error("invalid harmonic bond"));
        }
        bonds.push(e::Bond { i, j, k, r0 });
    }
    for a in records(api::field(v, "angles")?)? {
        let i = atom(a, "i", n)?;
        let j = atom(a, "j", n)?;
        let k = atom(a, "k", n)?;
        let ktheta = f.scalar(api::field(a, "constant")?, ENERGY, false)?;
        api::normalized(api::field(a, "theta0")?)?;
        let theta0 = crate::rules::numerical_abi::numeric(api::field(a, "theta0")?, 0)?;
        if i == j
            || j == k
            || i == k
            || ktheta < 0.0
            || !(0.0..=std::f64::consts::PI).contains(&theta0)
        {
            return Err(api::error("invalid harmonic angle"));
        }
        angles.push(e::Angle {
            i,
            j,
            k,
            ktheta,
            theta0,
        });
    }
    for d in records(api::field(v, "dihedrals")?)? {
        let atoms = [
            atom(d, "i", n)?,
            atom(d, "j", n)?,
            atom(d, "k", n)?,
            atom(d, "l", n)?,
        ];
        let mut unique = atoms.to_vec();
        unique.sort_unstable();
        unique.dedup();
        let k = f.scalar(api::field(d, "constant")?, ENERGY, false)?;
        api::normalized(api::field(d, "multiplicity")?)?;
        api::normalized(api::field(d, "phase")?)?;
        let multiplicity = crate::rules::numerical_abi::numeric(api::field(d, "multiplicity")?, 0)?;
        let delta = crate::rules::numerical_abi::numeric(api::field(d, "phase")?, 0)?;
        if unique.len() != 4 || k < 0.0 || multiplicity <= 0.0 || multiplicity > 32.0 {
            return Err(api::error("invalid periodic dihedral"));
        }
        dihedrals.push(e::Dihedral {
            atoms,
            k,
            n: multiplicity,
            delta,
        });
    }
    if D < 2 && !angles.is_empty() || D != 3 && !dihedrals.is_empty() {
        return Err(api::error(
            "angle/dihedral terms unsupported in this spatial dimension",
        ));
    }
    ff.bonds = bonds;
    ff.angles = angles;
    ff.dihedrals = dihedrals;
    ff.shift = <bool as Arg>::decode(api::field(v, "shift")?, 0)?;
    ff.shift_force = <bool as Arg>::decode(api::field(v, "shift_force")?, 0)?;
    ff.exclude_bonded = <bool as Arg>::decode(api::field(v, "exclude_bonded")?, 0)?;
    ff.rebuild_exclusions();
    Ok(())
}
fn engine_thermostat(
    f: &Frame,
    v: &Value,
    dof: f64,
) -> Result<e::Thermostat, OpError> {
    if matches!(v,Value::String(s)if s=="None") {
        return Ok(e::Thermostat::None);
    }
    let kind = <String as Arg>::decode(api::field(v, "kind")?, 0)?;
    let t0 = f.scalar(api::field(v, "temperature")?, ENERGY, true)?;
    if t0 < 0.0 {
        return Err(api::error("target temperature must be nonnegative"));
    }
    match kind.as_str() {
        | "Berendsen" => {
            let tau = f.scalar(api::field(v, "tau")?, TIME, false)?;
            if tau <= 0.0 {
                return Err(api::error("relaxation time must bepositive"));
            }
            Ok(e::Thermostat::Berendsen { t0, tau })
        },
        | "Langevin" => {
            let gamma = f.scalar(api::field(v, "friction")?, [0, 0, -1, 0], false)?;
            if gamma < 0.0 {
                return Err(api::error("friction must be nonnegative"));
            }
            Ok(e::Thermostat::Langevin { t0, gamma })
        },
        | "NoseHooverChain" => {
            let tau = f.scalar(api::field(v, "tau")?, TIME, false)?;
            let m = <usize as Arg>::decode(api::field(v, "length")?, 0)?;
            if t0 <= 0.0 || tau <= 0.0 || m == 0 || m > 16 || dof <= 0.0 {
                return Err(api::error("Nose Hoover controls invalid"));
            }
            Ok(e::Thermostat::NoseHooverChain(e::NoseHooverChain::new(
                t0, tau, m, dof,
            )))
        },
        | _ => Err(api::error("unknown engine thermostat")),
    }
}
object_binding!(engine_flexible_hydrostatic,"num_native_engine_flexible_hydrostatic","crate::sim::md_engine::FlexibleBarostat::hydrostatic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dimension:usize,pressure:Input,piston_mass:Input,mode:String)->NewObject<Flexible>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dimension:usize,pressure:Input,piston_mass:Input,mode:String|{if !(1..=3).contains(&dimension){return Err(api::error("barostat dimension must be1..3"));}let f=frame(ctx,&f)?;let w=f.scalar(&piston_mass.0,[1,2,0,0],false)?;if w<=0.0{return Err(api::error("piston mass must bepositive"));}let mode=cell_mode(&mode)?;let pressure_dimension=[1,2-dimension as i32,-2,0];let p=f.scalar(&pressure.0,pressure_dimension,false)?;let native=e::FlexibleBarostat::hydrostatic(dimension,p,w,mode);let dims=f.dims.to_vec();Ok(NewObject::new(Flexible{native,frame:f},dims))});
object_binding!(engine_flexible_with_tensor,"num_native_engine_flexible_with_tensor","crate::sim::md_engine::FlexibleBarostat::with_tensor","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dimension:usize,pressure:Input,piston_mass:Input,mode:String)->NewObject<Flexible>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dimension:usize,pressure:Input,piston_mass:Input,mode:String|{if !(1..=3).contains(&dimension){return Err(api::error("barostat dimension must be1..3"));}let f=frame(ctx,&f)?;let w=f.scalar(&piston_mass.0,[1,2,0,0],false)?;if w<=0.0{return Err(api::error("piston mass must bepositive"));}let mode=cell_mode(&mode)?;let pressure_dimension=[1,2-dimension as i32,-2,0];let Value::Tensor(t)=&pressure.0 else{return Err(api::error("pressure tensor expected"));};
if t.shape!=[dimension,dimension]{return Err(api::error("pressure tensor dimensions differ"));}api::same_dimension(&api::dimension(&pressure.0)?,&f.dim(pressure_dimension)?)?;let p=t.data.iter().map(|v|f.apply(crate::rules::numerical_abi::numeric(v,0)?,pressure_dimension,false,false)).collect::<Result<Vec<_>,_>>()?;let native=e::FlexibleBarostat::with_tensor(dimension,p,w,mode);let dims=f.dims.to_vec();Ok(NewObject::new(Flexible{native,frame:f},dims))});
object_binding!(engine_flexible_shear_flow,"num_native_engine_flexible_shear_flow","crate::sim::md_engine::FlexibleBarostat::shear_flow","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dimension:usize,rate:Input)->NewObject<Flexible>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dimension:usize,rate:Input|{if !(1..=3).contains(&dimension){return Err(api::error("barostat dimension must be1..3"));}let f=frame(ctx,&f)?;if dimension<2{return Err(api::error("shear flow needs dimension>=2"));}let rate=f.scalar(&rate.0,[0,0,-1,0],false)?;let native=e::FlexibleBarostat::shear_flow(dimension,rate);let dims=f.dims.to_vec();Ok(NewObject::new(Flexible{native,frame:f},dims))});
object_binding!(engine_flexible_with_langevin_cell,"num_native_engine_flexible_with_langevin_cell","crate::sim::md_engine::FlexibleBarostat::with_langevin_cell","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(barostat:Handle<Flexible>,friction:Input)->NewObject<Flexible>,|ctx:&CallCtx<'_>,barostat:Handle<Flexible>,friction:Input|{barostat.with(ctx,|b|{let gamma=b.frame.scalar(&friction.0,[0,0,-1,0],false)?;if gamma<0.0{return Err(api::error("cell friction must be nonnegative"));}Ok(NewObject::new(Flexible{native:b.native.clone().with_langevin_cell(gamma),frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(engine_flexible_with_chain,"num_native_engine_flexible_with_chain","crate::sim::md_engine::FlexibleBarostat::with_chain","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(barostat:Handle<Flexible>,temperature:Input,tau:Input,length:usize)->NewObject<Flexible>,|ctx:&CallCtx<'_>,barostat:Handle<Flexible>,temperature:Input,tau:Input,length:usize|{barostat.with(ctx,|b|{let t=b.frame.scalar(&temperature.0,ENERGY,true)?;let tau=b.frame.scalar(&tau.0,TIME,false)?;if t<=0.0||tau<=0.0||length==0||length>16{return Err(api::error("cell chain controls invalid"));}Ok(NewObject::new(Flexible{native:b.native.clone().with_chain(t,tau,length),frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(engine_flexible_cell_dof,"num_native_engine_flexible_cell_dof","crate::sim::md_engine::FlexibleBarostat::cell_dof","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(barostat:Handle<Flexible>)->usize,|ctx:&CallCtx<'_>,barostat:Handle<Flexible>|{barostat.with(ctx,|b|Ok(b.native.cell_dof()))});
object_binding!(engine_flexible_energy,"num_native_engine_flexible_energy","crate::sim::md_engine::FlexibleBarostat::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(barostat:Handle<Flexible>,volume:Input)->Dynamic,|ctx:&CallCtx<'_>,barostat:Handle<Flexible>,volume:Input|{barostat.with(ctx,|b|{let v=b.frame.scalar(&volume.0,[0,b.native.d as i32,0,0],false)?;if v<=0.0{return Err(api::error("cell volume must bepositive"));}Ok(Dynamic(b.frame.output(b.native.energy(v),ENERGY,false)?))})});
object_binding!(engine_forcefield_set_bonded_d1,"num_native_engine_forcefield_set_bonded_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(field:Handle<Field<1>>,config:Input,particle_count:usize)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<1>>,config:Input,particle_count:usize|{bound(particle_count,512)?;field.with_mut(ctx,|f|{let ff=f.native.as_mut().ok_or_else(||api::error("force field consumed"))?;set_bonded(&f.frame,ff,&config.0,particle_count)?;f.validate()?;Ok(true)})});
object_binding!(engine_simulation_set_thermostat_d1,"num_native_engine_simulation_set_thermostat_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<1>>,control:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>,control:Input|{simulation.with_mut(ctx,|s|{s.native.thermostat=engine_thermostat(&s.frame,&control.0,s.native.dof())?;Ok(true)})});
object_binding!(engine_simulation_set_barostat_d1,"num_native_engine_simulation_set_barostat_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<1>>,barostat:Handle<Flexible>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>,barostat:Handle<Flexible>|{simulation.with_mut(ctx,|s|barostat.with(ctx,|b|{s.frame.same(&b.frame)?;if b.native.d!=1{return Err(api::error("barostat dimension differs fromsimulation"));}s.native.barostat=e::Barostat::Flexible(b.native.clone());Ok(true)}))});
object_binding!(engine_simulation_set_constraints_d1,"num_native_engine_simulation_set_constraints_d1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<1>>,constraints:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>,constraints:Input|{simulation.with_mut(ctx,|s|{let values=records(&constraints.0)?;if values.len()>64||values.len()>s.native.sys.n().saturating_mul(1).saturating_sub(2){return Err(api::error("constraint count exceedswork/DOF budget"));}let mut constraints=vec![];let mut pairs=std::collections::BTreeSet::new();for v in values{let i=atom(v,"i",s.native.sys.n())?;let j=atom(v,"j",s.native.sys.n())?;let d=s.frame.scalar(api::field(v,"distance")?,LENGTH,false)?;if i==j||d<=0.0||!pairs.insert((i.min(j),i.max(j))){return Err(api::error("constraint needs distinctatoms/positive distance/unique pair"));}constraints.push(e::Constraint{i,j,d});}s.native.constraints=constraints;s.native.refresh_constraint_virial();Ok(true)})});
object_binding!(engine_simulation_refresh_constraint_virial_d1,"num_native_engine_simulation_refresh_constraint_virial_d1","crate::sim::md_engine::MdSimulation<D>::refresh_constraint_virial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<1>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<1>>|{simulation.with_mut(ctx,|s|{s.native.refresh_constraint_virial();finite(s.native.constraint_virial.into_iter().flatten())?;Ok(true)})});
object_binding!(engine_forcefield_set_bonded_d2,"num_native_engine_forcefield_set_bonded_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(field:Handle<Field<2>>,config:Input,particle_count:usize)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<2>>,config:Input,particle_count:usize|{bound(particle_count,512)?;field.with_mut(ctx,|f|{let ff=f.native.as_mut().ok_or_else(||api::error("force field consumed"))?;set_bonded(&f.frame,ff,&config.0,particle_count)?;f.validate()?;Ok(true)})});
object_binding!(engine_simulation_set_thermostat_d2,"num_native_engine_simulation_set_thermostat_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<2>>,control:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>,control:Input|{simulation.with_mut(ctx,|s|{s.native.thermostat=engine_thermostat(&s.frame,&control.0,s.native.dof())?;Ok(true)})});
object_binding!(engine_simulation_set_barostat_d2,"num_native_engine_simulation_set_barostat_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<2>>,barostat:Handle<Flexible>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>,barostat:Handle<Flexible>|{simulation.with_mut(ctx,|s|barostat.with(ctx,|b|{s.frame.same(&b.frame)?;if b.native.d!=2{return Err(api::error("barostat dimension differs fromsimulation"));}s.native.barostat=e::Barostat::Flexible(b.native.clone());Ok(true)}))});
object_binding!(engine_simulation_set_constraints_d2,"num_native_engine_simulation_set_constraints_d2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<2>>,constraints:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>,constraints:Input|{simulation.with_mut(ctx,|s|{let values=records(&constraints.0)?;if values.len()>64||values.len()>s.native.sys.n().saturating_mul(2).saturating_sub(3){return Err(api::error("constraint count exceedswork/DOF budget"));}let mut constraints=vec![];let mut pairs=std::collections::BTreeSet::new();for v in values{let i=atom(v,"i",s.native.sys.n())?;let j=atom(v,"j",s.native.sys.n())?;let d=s.frame.scalar(api::field(v,"distance")?,LENGTH,false)?;if i==j||d<=0.0||!pairs.insert((i.min(j),i.max(j))){return Err(api::error("constraint needs distinctatoms/positive distance/unique pair"));}constraints.push(e::Constraint{i,j,d});}s.native.constraints=constraints;s.native.refresh_constraint_virial();Ok(true)})});
object_binding!(engine_simulation_refresh_constraint_virial_d2,"num_native_engine_simulation_refresh_constraint_virial_d2","crate::sim::md_engine::MdSimulation<D>::refresh_constraint_virial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<2>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<2>>|{simulation.with_mut(ctx,|s|{s.native.refresh_constraint_virial();finite(s.native.constraint_virial.into_iter().flatten())?;Ok(true)})});
object_binding!(engine_forcefield_set_bonded_d3,"num_native_engine_forcefield_set_bonded_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(field:Handle<Field<3>>,config:Input,particle_count:usize)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<3>>,config:Input,particle_count:usize|{bound(particle_count,512)?;field.with_mut(ctx,|f|{let ff=f.native.as_mut().ok_or_else(||api::error("force field consumed"))?;set_bonded(&f.frame,ff,&config.0,particle_count)?;f.validate()?;Ok(true)})});
object_binding!(engine_simulation_set_thermostat_d3,"num_native_engine_simulation_set_thermostat_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<3>>,control:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>,control:Input|{simulation.with_mut(ctx,|s|{s.native.thermostat=engine_thermostat(&s.frame,&control.0,s.native.dof())?;Ok(true)})});
object_binding!(engine_simulation_set_barostat_d3,"num_native_engine_simulation_set_barostat_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<3>>,barostat:Handle<Flexible>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>,barostat:Handle<Flexible>|{simulation.with_mut(ctx,|s|barostat.with(ctx,|b|{s.frame.same(&b.frame)?;if b.native.d!=3{return Err(api::error("barostat dimension differs fromsimulation"));}s.native.barostat=e::Barostat::Flexible(b.native.clone());Ok(true)}))});
object_binding!(engine_simulation_set_constraints_d3,"num_native_engine_simulation_set_constraints_d3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<3>>,constraints:Input)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>,constraints:Input|{simulation.with_mut(ctx,|s|{let values=records(&constraints.0)?;if values.len()>64||values.len()>s.native.sys.n().saturating_mul(3).saturating_sub(4){return Err(api::error("constraint count exceedswork/DOF budget"));}let mut constraints=vec![];let mut pairs=std::collections::BTreeSet::new();for v in values{let i=atom(v,"i",s.native.sys.n())?;let j=atom(v,"j",s.native.sys.n())?;let d=s.frame.scalar(api::field(v,"distance")?,LENGTH,false)?;if i==j||d<=0.0||!pairs.insert((i.min(j),i.max(j))){return Err(api::error("constraint needs distinctatoms/positive distance/unique pair"));}constraints.push(e::Constraint{i,j,d});}s.native.constraints=constraints;s.native.refresh_constraint_virial();Ok(true)})});
object_binding!(engine_simulation_refresh_constraint_virial_d3,"num_native_engine_simulation_refresh_constraint_virial_d3","crate::sim::md_engine::MdSimulation<D>::refresh_constraint_virial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Simulation<3>>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Simulation<3>>|{simulation.with_mut(ctx,|s|{s.native.refresh_constraint_virial();finite(s.native.constraint_virial.into_iter().flatten())?;Ok(true)})});
fn plain_matrix<const D: usize>(v: &Input) -> Result<[[f64; D]; D], OpError> {
    let rows = <Vec<Vec<f64>> as Arg>::decode(&v.0, 0)?;
    if rows.len() != D || rows.iter().any(|r| r.len() != D) {
        return Err(api::error("cell helper matrix must beD-by-D"));
    }
    rows.into_iter()
        .map(|r| {
            r.try_into()
                .map_err(|_| api::error("cell helper matrix row"))
        })
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| api::error("cell helper matrix shape"))
}
fn plain_output<const D: usize>(
    v: [[f64; D]; D],
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![D, D],
        data: v
            .into_iter()
            .flatten()
            .map(|x| api::quantity(x, d.clone()))
            .collect::<Result<_, _>>()?,
        elem: TypeExpr::Float.with_dim(d),
    })))
}
fn canonical_dimension(d: &DimExpr) -> Result<DimExpr, OpError> {
    super::super::runtime_dimensions::expression(&super::super::runtime_dimensions::canonical(d)?)
}
object_binding!(engine_mat_mul_d1,"num_native_engine_mat_mul_d1","crate::sim::md_engine::mat_mul","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<1>(&left)?;let d=api::dimension(&left.0)?;let b=plain_matrix::<1>(&right)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;plain_output(e::mat_mul(&a,&b),dim)});
object_binding!(engine_mat_vec_d1,"num_native_engine_mat_vec_d1","crate::sim::md_engine::mat_vec","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<1>(&left)?;let d=api::dimension(&left.0)?;let b=<[f64;1] as Arg>::decode(&right.0,0)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;let out=e::mat_vec(&a,&b);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![1],data:out.into_iter().map(|x|api::quantity(x,dim.clone())).collect::<Result<_,_>>()?,elem:TypeExpr::Float.with_dim(dim)})))});
object_binding!(engine_mat_det_d1,"num_native_engine_mat_det_d1","crate::sim::md_engine::mat_det","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<1>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(d.pow(Exp::Rat(Rational::int(1)))))?;Ok(Dynamic(api::quantity(e::mat_det(&a),dim)?))});
object_binding!(engine_mat_inverse_d1,"num_native_engine_mat_inverse_d1","crate::sim::md_engine::mat_inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<1>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(DimExpr::Inverse(Box::new(d))))?;let out=e::mat_inverse(&a).map(|v|plain_output(v,dim).map(|v|v.0)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if out.is_some(){"Some"}else{"None"}.into(),payload:out.into_iter().collect()})))});
object_binding!(engine_mat_exp_phi1_d1,"num_native_engine_mat_exp_phi1_d1","crate::sim::md_engine::mat_exp_phi1","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<1>(&left)?;api::normalized(&left.0)?;let(a,b)=e::mat_exp_phi1(&a);Ok(Dynamic(api::record(vec![("exponential",plain_output(a,DimExpr::Dimless)?.0),("phi1",plain_output(b,DimExpr::Dimless)?.0)])))});
object_binding!(engine_mat_mul_d2,"num_native_engine_mat_mul_d2","crate::sim::md_engine::mat_mul","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<2>(&left)?;let d=api::dimension(&left.0)?;let b=plain_matrix::<2>(&right)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;plain_output(e::mat_mul(&a,&b),dim)});
object_binding!(engine_mat_vec_d2,"num_native_engine_mat_vec_d2","crate::sim::md_engine::mat_vec","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<2>(&left)?;let d=api::dimension(&left.0)?;let b=<[f64;2] as Arg>::decode(&right.0,0)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;let out=e::mat_vec(&a,&b);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![2],data:out.into_iter().map(|x|api::quantity(x,dim.clone())).collect::<Result<_,_>>()?,elem:TypeExpr::Float.with_dim(dim)})))});
object_binding!(engine_mat_det_d2,"num_native_engine_mat_det_d2","crate::sim::md_engine::mat_det","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<2>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(d.pow(Exp::Rat(Rational::int(2)))))?;Ok(Dynamic(api::quantity(e::mat_det(&a),dim)?))});
object_binding!(engine_mat_inverse_d2,"num_native_engine_mat_inverse_d2","crate::sim::md_engine::mat_inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<2>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(DimExpr::Inverse(Box::new(d))))?;let out=e::mat_inverse(&a).map(|v|plain_output(v,dim).map(|v|v.0)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if out.is_some(){"Some"}else{"None"}.into(),payload:out.into_iter().collect()})))});
object_binding!(engine_mat_exp_phi1_d2,"num_native_engine_mat_exp_phi1_d2","crate::sim::md_engine::mat_exp_phi1","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<2>(&left)?;api::normalized(&left.0)?;let(a,b)=e::mat_exp_phi1(&a);Ok(Dynamic(api::record(vec![("exponential",plain_output(a,DimExpr::Dimless)?.0),("phi1",plain_output(b,DimExpr::Dimless)?.0)])))});
object_binding!(engine_mat_mul_d3,"num_native_engine_mat_mul_d3","crate::sim::md_engine::mat_mul","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<3>(&left)?;let d=api::dimension(&left.0)?;let b=plain_matrix::<3>(&right)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;plain_output(e::mat_mul(&a,&b),dim)});
object_binding!(engine_mat_vec_d3,"num_native_engine_mat_vec_d3","crate::sim::md_engine::mat_vec","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input,right:Input|{let a=plain_matrix::<3>(&left)?;let d=api::dimension(&left.0)?;let b=<[f64;3] as Arg>::decode(&right.0,0)?;let dim=canonical_dimension(&(DimExpr::Mul(vec![d,api::dimension(&right.0)?])))?;let out=e::mat_vec(&a,&b);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],data:out.into_iter().map(|x|api::quantity(x,dim.clone())).collect::<Result<_,_>>()?,elem:TypeExpr::Float.with_dim(dim)})))});
object_binding!(engine_mat_det_d3,"num_native_engine_mat_det_d3","crate::sim::md_engine::mat_det","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<3>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(d.pow(Exp::Rat(Rational::int(3)))))?;Ok(Dynamic(api::quantity(e::mat_det(&a),dim)?))});
object_binding!(engine_mat_inverse_d3,"num_native_engine_mat_inverse_d3","crate::sim::md_engine::mat_inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<3>(&left)?;let d=api::dimension(&left.0)?;let dim=canonical_dimension(&(DimExpr::Inverse(Box::new(d))))?;let out=e::mat_inverse(&a).map(|v|plain_output(v,dim).map(|v|v.0)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if out.is_some(){"Some"}else{"None"}.into(),payload:out.into_iter().collect()})))});
object_binding!(engine_mat_exp_phi1_d3,"num_native_engine_mat_exp_phi1_d3","crate::sim::md_engine::mat_exp_phi1","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input)->Dynamic,|_ctx:&CallCtx<'_>,left:Input|{let a=plain_matrix::<3>(&left)?;api::normalized(&left.0)?;let(a,b)=e::mat_exp_phi1(&a);Ok(Dynamic(api::record(vec![("exponential",plain_output(a,DimExpr::Dimless)?.0),("phi1",plain_output(b,DimExpr::Dimless)?.0)])))});
pub(super) fn field_indices<const D: usize>(
    ff: &e::ForceField<D>,
    n: usize,
) -> Result<(), OpError> {
    if ff.bonds.iter().any(|b| b.i >= n || b.j >= n)
        || ff.angles.iter().any(|a| a.i >= n || a.j >= n || a.k >= n)
        || ff.dihedrals.iter().any(|d| d.atoms.iter().any(|i| *i >= n))
        || ff.extra_exclusions.iter().any(|(i, j)| *i >= n || *j >= n)
    {
        return Err(api::error(
            "force field bonded/excluded indices exceed actual system size",
        ));
    }
    Ok(())
}
macro_rules! radial_merge_binding {
    ($module:ident, $name:literal, $dim:literal) => {
        object_binding!(
            $module,
            $name,
            "crate::sim::md_engine::Rdf<D>::merge",
            "Combine radial distribution frame counts and per-frame ideal-gas references in the destination. Both accumulators must use the same bins, radius and physical reference frame. The source can be an independent accumulator.",
            false,
            (destination: Handle<Radial<$dim>>, source: Handle<Radial<$dim>>) -> bool,
            |ctx, destination, source| merge_radial(ctx, &destination, &source)
        );
    };
}
radial_merge_binding!(engine_rdf_merge_d1, "num_native_engine_rdf_merge_d1", 1);
radial_merge_binding!(engine_rdf_merge_d2, "num_native_engine_rdf_merge_d2", 2);
radial_merge_binding!(engine_rdf_merge_d3, "num_native_engine_rdf_merge_d3", 3);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &engine_flexible_hydrostatic::BINDING,
            &engine_flexible_with_tensor::BINDING,
            &engine_flexible_shear_flow::BINDING,
            &engine_flexible_with_langevin_cell::BINDING,
            &engine_flexible_with_chain::BINDING,
            &engine_flexible_cell_dof::BINDING,
            &engine_flexible_energy::BINDING,
            &engine_forcefield_set_bonded_d1::BINDING,
            &engine_simulation_set_thermostat_d1::BINDING,
            &engine_simulation_set_barostat_d1::BINDING,
            &engine_simulation_set_constraints_d1::BINDING,
            &engine_simulation_refresh_constraint_virial_d1::BINDING,
            &engine_forcefield_set_bonded_d2::BINDING,
            &engine_simulation_set_thermostat_d2::BINDING,
            &engine_simulation_set_barostat_d2::BINDING,
            &engine_simulation_set_constraints_d2::BINDING,
            &engine_simulation_refresh_constraint_virial_d2::BINDING,
            &engine_forcefield_set_bonded_d3::BINDING,
            &engine_simulation_set_thermostat_d3::BINDING,
            &engine_simulation_set_barostat_d3::BINDING,
            &engine_simulation_set_constraints_d3::BINDING,
            &engine_simulation_refresh_constraint_virial_d3::BINDING,
            &engine_mat_mul_d1::BINDING,
            &engine_mat_vec_d1::BINDING,
            &engine_mat_det_d1::BINDING,
            &engine_mat_inverse_d1::BINDING,
            &engine_mat_exp_phi1_d1::BINDING,
            &engine_mat_mul_d2::BINDING,
            &engine_mat_vec_d2::BINDING,
            &engine_mat_det_d2::BINDING,
            &engine_mat_inverse_d2::BINDING,
            &engine_mat_exp_phi1_d2::BINDING,
            &engine_mat_mul_d3::BINDING,
            &engine_mat_vec_d3::BINDING,
            &engine_mat_det_d3::BINDING,
            &engine_mat_inverse_d3::BINDING,
            &engine_mat_exp_phi1_d3::BINDING,
            &engine_pair_new::BINDING,
            &engine_pair_eval::BINDING,
            &engine_box_periodic_d1::BINDING,
            &engine_box_new_d1::BINDING,
            &engine_box_triclinic_d1::BINDING,
            &engine_box_set_cell_d1::BINDING,
            &engine_box_matrix_d1::BINDING,
            &engine_box_volume_d1::BINDING,
            &engine_box_height_d1::BINDING,
            &engine_box_to_fractional_d1::BINDING,
            &engine_box_from_fractional_d1::BINDING,
            &engine_box_min_image_d1::BINDING,
            &engine_system_new_d1::BINDING,
            &engine_system_data_d1::BINDING,
            &engine_system_set_state_d1::BINDING,
            &engine_system_n_d1::BINDING,
            &engine_system_kinetic_energy_d1::BINDING,
            &engine_system_momentum_d1::BINDING,
            &engine_system_temperature_d1::BINDING,
            &engine_system_remove_com_velocity_d1::BINDING,
            &engine_system_maxwell_boltzmann_d1::BINDING,
            &engine_system_apply_boundaries_d1::BINDING,
            &engine_lattice_cubic_d1::BINDING,
            &engine_random_positions_d1::BINDING,
            &engine_forcefield_new_d1::BINDING,
            &engine_forcefield_validate_d1::BINDING,
            &engine_forcefield_rebuild_exclusions_d1::BINDING,
            &engine_forces_zeros_d1::BINDING,
            &engine_forces_data_d1::BINDING,
            &engine_neighbours_new_d1::BINDING,
            &engine_neighbours_needs_rebuild_d1::BINDING,
            &engine_neighbours_build_d1::BINDING,
            &engine_compute_forces_d1::BINDING,
            &engine_simulation_new_d1::BINDING,
            &engine_simulation_dof_d1::BINDING,
            &engine_simulation_temperature_d1::BINDING,
            &engine_simulation_total_energy_d1::BINDING,
            &engine_simulation_conserved_energy_d1::BINDING,
            &engine_simulation_pressure_d1::BINDING,
            &engine_simulation_pressure_tensor_d1::BINDING,
            &engine_simulation_step_d1::BINDING,
            &engine_simulation_run_d1::BINDING,
            &engine_simulation_state_d1::BINDING,
            &engine_rdf_new_d1::BINDING,
            &engine_rdf_add_frame_d1::BINDING,
            &engine_rdf_merge_d1::BINDING,
            &engine_rdf_result_d1::BINDING,
            &engine_msd_d1::BINDING,
            &engine_vacf_d1::BINDING,
            &engine_box_periodic_d2::BINDING,
            &engine_box_new_d2::BINDING,
            &engine_box_triclinic_d2::BINDING,
            &engine_box_set_cell_d2::BINDING,
            &engine_box_matrix_d2::BINDING,
            &engine_box_volume_d2::BINDING,
            &engine_box_height_d2::BINDING,
            &engine_box_to_fractional_d2::BINDING,
            &engine_box_from_fractional_d2::BINDING,
            &engine_box_min_image_d2::BINDING,
            &engine_system_new_d2::BINDING,
            &engine_system_data_d2::BINDING,
            &engine_system_set_state_d2::BINDING,
            &engine_system_n_d2::BINDING,
            &engine_system_kinetic_energy_d2::BINDING,
            &engine_system_momentum_d2::BINDING,
            &engine_system_temperature_d2::BINDING,
            &engine_system_remove_com_velocity_d2::BINDING,
            &engine_system_maxwell_boltzmann_d2::BINDING,
            &engine_system_apply_boundaries_d2::BINDING,
            &engine_lattice_cubic_d2::BINDING,
            &engine_random_positions_d2::BINDING,
            &engine_forcefield_new_d2::BINDING,
            &engine_forcefield_validate_d2::BINDING,
            &engine_forcefield_rebuild_exclusions_d2::BINDING,
            &engine_forces_zeros_d2::BINDING,
            &engine_forces_data_d2::BINDING,
            &engine_neighbours_new_d2::BINDING,
            &engine_neighbours_needs_rebuild_d2::BINDING,
            &engine_neighbours_build_d2::BINDING,
            &engine_compute_forces_d2::BINDING,
            &engine_simulation_new_d2::BINDING,
            &engine_simulation_dof_d2::BINDING,
            &engine_simulation_temperature_d2::BINDING,
            &engine_simulation_total_energy_d2::BINDING,
            &engine_simulation_conserved_energy_d2::BINDING,
            &engine_simulation_pressure_d2::BINDING,
            &engine_simulation_pressure_tensor_d2::BINDING,
            &engine_simulation_step_d2::BINDING,
            &engine_simulation_run_d2::BINDING,
            &engine_simulation_state_d2::BINDING,
            &engine_rdf_new_d2::BINDING,
            &engine_rdf_add_frame_d2::BINDING,
            &engine_rdf_merge_d2::BINDING,
            &engine_rdf_result_d2::BINDING,
            &engine_msd_d2::BINDING,
            &engine_vacf_d2::BINDING,
            &engine_box_periodic_d3::BINDING,
            &engine_box_new_d3::BINDING,
            &engine_box_triclinic_d3::BINDING,
            &engine_box_set_cell_d3::BINDING,
            &engine_box_matrix_d3::BINDING,
            &engine_box_volume_d3::BINDING,
            &engine_box_height_d3::BINDING,
            &engine_box_to_fractional_d3::BINDING,
            &engine_box_from_fractional_d3::BINDING,
            &engine_box_min_image_d3::BINDING,
            &engine_system_new_d3::BINDING,
            &engine_system_data_d3::BINDING,
            &engine_system_set_state_d3::BINDING,
            &engine_system_n_d3::BINDING,
            &engine_system_kinetic_energy_d3::BINDING,
            &engine_system_momentum_d3::BINDING,
            &engine_system_temperature_d3::BINDING,
            &engine_system_remove_com_velocity_d3::BINDING,
            &engine_system_maxwell_boltzmann_d3::BINDING,
            &engine_system_apply_boundaries_d3::BINDING,
            &engine_lattice_cubic_d3::BINDING,
            &engine_random_positions_d3::BINDING,
            &engine_forcefield_new_d3::BINDING,
            &engine_forcefield_validate_d3::BINDING,
            &engine_forcefield_rebuild_exclusions_d3::BINDING,
            &engine_forces_zeros_d3::BINDING,
            &engine_forces_data_d3::BINDING,
            &engine_neighbours_new_d3::BINDING,
            &engine_neighbours_needs_rebuild_d3::BINDING,
            &engine_neighbours_build_d3::BINDING,
            &engine_compute_forces_d3::BINDING,
            &engine_simulation_new_d3::BINDING,
            &engine_simulation_dof_d3::BINDING,
            &engine_simulation_temperature_d3::BINDING,
            &engine_simulation_total_energy_d3::BINDING,
            &engine_simulation_conserved_energy_d3::BINDING,
            &engine_simulation_pressure_d3::BINDING,
            &engine_simulation_pressure_tensor_d3::BINDING,
            &engine_simulation_step_d3::BINDING,
            &engine_simulation_run_d3::BINDING,
            &engine_simulation_state_d3::BINDING,
            &engine_rdf_new_d3::BINDING,
            &engine_rdf_add_frame_d3::BINDING,
            &engine_rdf_merge_d3::BINDING,
            &engine_rdf_result_d3::BINDING,
            &engine_msd_d3::BINDING,
            &engine_vacf_d3::BINDING,
            &engine_lattice_fcc::BINDING,
            &engine_lattice_triangular::BINDING,
            &engine_energy_drift_slope::BINDING,
            &engine_chain_new::BINDING,
            &engine_chain_propagate::BINDING,
            &engine_chain_energy::BINDING,
        ],
        "numerical_vector",
    )
}
