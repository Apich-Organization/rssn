#![allow(clippy::wildcard_imports)] // Site forces and molecular rigid degrees of freedom use the same native physical frame.
use super::*;
use crate::sim::md_rigid;
struct Rigid {
    native: md_rigid::RigidMdSimulation,
    frame: Frame,
    extra_work: usize,
}
impl NativeType for Rigid {
    const KIND: &'static str = "rssn/native/sim/rigid-md/v1";
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        engine::check_system(&self.native.sys)?;
        engine::bound(self.native.mols.len(), 512)?;
        engine::finite([self.native.total_energy(), self.native.max_geometry_error()])
    }
}
fn molecules(
    v: &Input,
    n: usize,
) -> Result<Vec<Vec<usize>>, OpError> {
    api::normalized(&v.0)?;
    let offsets = <Vec<usize> as Arg>::decode(api::field(&v.0, "offsets")?, 0)?;
    let indices = <Vec<usize> as Arg>::decode(api::field(&v.0, "indices")?, 0)?;
    if offsets.len() < 2
        || offsets.len() > 513
        || offsets.first() != Some(&0)
        || offsets.last() != Some(&indices.len())
        || offsets.windows(2).any(|w| w[0] >= w[1])
        || indices.len() != n
    {
        return Err(api::error(
            "molecular CSR must partition every site into nonempty rows",
        ));
    }
    let mut sorted = indices.clone();
    sorted.sort_unstable();
    if sorted != (0..n).collect::<Vec<_>>() {
        return Err(api::error(
            "molecular CSR must contain each site exactlyonce",
        ));
    }
    Ok(offsets
        .windows(2)
        .map(|w| indices[w[0]..w[1]].to_vec())
        .collect())
}
object_binding!(rigid_md_new,"num_native_rigid_md_new","crate::sim::md_rigid::RigidMdSimulation::new","Construct the actual native rigid-molecule site driver. Checked CSR partitions identify rigid molecules; the physical forcefield is transferred once. Unitful centers/sites, translational and rotational DOFs retain mass-length-time dimensions; rotational integration remains native.",false,(system:Handle<System<3>>,field:Handle<Field<3>>,skin:Input,dt:Input,partition:Input)->NewObject<Rigid>,|ctx:&CallCtx<'_>,system:Handle<System<3>>,field:Handle<Field<3>>,skin:Input,dt:Input,partition:Input|{system.with(ctx,|s|field.with_mut(ctx,|f|{s.frame.same(&f.frame)?;let skin=s.frame.scalar(&skin.0,LENGTH,false)?;let dt=s.frame.scalar(&dt.0,TIME,false)?;if skin<0.0||dt<=0.0{return Err(api::error("rigid MD time step/skin invalid"));}let molecules=molecules(&partition,s.native.n())?;if s.native.bx.bc.iter().any(|b|*b!=e::Boundary::Periodic)||s.native.bx.shear_rate!=0.0{return Err(api::error("rigid MD needs fullyperiodic nonsheared box"));}let ff=f.native.as_ref().ok_or_else(||api::error("force field consumed"))?;engine::field_indices(ff,s.native.n())?;ff.validate().map_err(api::error)?;let list=e::NeighbourList::new(ff.rcut,skin);engine::safe_grid(&s.native,&list)?;let ff=f.native.take().ok_or_else(||api::error("force field consumed"))?;let native=md_rigid::RigidMdSimulation::new(s.native.clone(),ff,skin,dt,&molecules).map_err(api::error)?;Ok(NewObject::new(Rigid{native,frame:s.frame.clone(),extra_work:f.extra_work},s.frame.dims.to_vec()))}))});
object_binding!(rigid_md_kinetic_energy,"num_native_rigid_md_kinetic_energy","crate::sim::md_rigid::RigidMdSimulation::kinetic_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.kinetic_energy(),ENERGY,false)?))})});
object_binding!(rigid_md_total_energy,"num_native_rigid_md_total_energy","crate::sim::md_rigid::RigidMdSimulation::total_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.total_energy(),ENERGY,false)?))})});
object_binding!(rigid_md_dof,"num_native_rigid_md_dof","crate::sim::md_rigid::RigidMdSimulation::dof","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->f64,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(s.native.dof())})});
object_binding!(rigid_md_temperature,"num_native_rigid_md_temperature","crate::sim::md_rigid::RigidMdSimulation::temperature","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{if s.native.dof()<=0.0{return Err(api::error("rigid temperature requires positive DOF"));}Ok(Dynamic(s.frame.output(s.native.temperature(),ENERGY,true)?))})});
object_binding!(rigid_md_momentum,"num_native_rigid_md_momentum","crate::sim::md_rigid::RigidMdSimulation::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.array(s.native.momentum().to_vec(),MOMENTUM)?))})});
object_binding!(rigid_md_pressure_tensor,"num_native_rigid_md_pressure_tensor","crate::sim::md_rigid::RigidMdSimulation::pressure_tensor","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(Dynamic(engine::packed(&s.frame,s.native.pressure_tensor().to_vec(),[1,-1,-2,0])?))})});
object_binding!(rigid_md_max_geometry_error,"num_native_rigid_md_max_geometry_error","crate::sim::md_rigid::RigidMdSimulation::max_geometry_error","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.max_geometry_error(),LENGTH,false)?))})});
object_binding!(rigid_md_sync_sites,"num_native_rigid_md_sync_sites","crate::sim::md_rigid::RigidMdSimulation::sync_sites","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with_mut(ctx,|s|{s.native.sync_sites();s.validate()?;Ok(true)})});
object_binding!(rigid_md_step,"num_native_rigid_md_step","crate::sim::md_rigid::RigidMdSimulation::step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with_mut(ctx,|s|{let steps=1usize;if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>50_000_000){return Err(api::error("rigid force-step work bound exceeded"));}s.native.step().map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(rigid_md_run,"num_native_rigid_md_run","crate::sim::md_rigid::RigidMdSimulation::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>,steps:usize)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>,steps:usize|{simulation.with_mut(ctx,|s|{if steps>65536||s.native.sys.n().checked_mul(s.native.sys.n()).and_then(|n|n.checked_add(s.extra_work)).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>50_000_000){return Err(api::error("rigid force-step work bound exceeded"));}s.native.run(steps).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(rigid_md_controls,"num_native_rigid_md_controls","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>,order:usize,top_split:bool,exact_rotor:bool)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>,order:usize,top_split:bool,exact_rotor:bool|{if !matches!(order,2|4|6|8){return Err(api::error("rotational composition order must be2/4/6/8"));}simulation.with_mut(ctx,|s|{s.native.order=order;s.native.top_split=top_split;s.native.exact_rotor=exact_rotor;Ok(true)})});
object_binding!(rigid_md_state,"num_native_rigid_md_state","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Rigid>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Rigid>|{simulation.with(ctx,|s|Ok(Dynamic(engine::system_data(&s.native.sys,&s.frame)?)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &rigid_md_new::BINDING,
            &rigid_md_kinetic_energy::BINDING,
            &rigid_md_total_energy::BINDING,
            &rigid_md_dof::BINDING,
            &rigid_md_temperature::BINDING,
            &rigid_md_momentum::BINDING,
            &rigid_md_pressure_tensor::BINDING,
            &rigid_md_max_geometry_error::BINDING,
            &rigid_md_sync_sites::BINDING,
            &rigid_md_step::BINDING,
            &rigid_md_run::BINDING,
            &rigid_md_controls::BINDING,
            &rigid_md_state::BINDING,
        ],
        "numerical_vector",
    )
}
