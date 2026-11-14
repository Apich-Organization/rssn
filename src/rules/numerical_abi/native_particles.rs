//! Physical and reduced-unit specializations of the native particle API.
mod classical;
mod dynamics;
mod engine;
mod multigrid;
mod spin;
mod wave;
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::TensorValue;
use crate::sim::physics_md as md;
const KB: f64 = 1.380_649e-23;
const MASS: [i32; 4] = [1, 0, 0, 0];
const LENGTH: [i32; 4] = [0, 1, 0, 0];
const TIME: [i32; 4] = [0, 0, 1, 0];
const CHARGE: [i32; 4] = [0, 0, 0, 1];
const VELOCITY: [i32; 4] = [0, 1, -1, 0];
const FORCE: [i32; 4] = [1, 1, -2, 0];
const ENERGY: [i32; 4] = [1, 2, -2, 0];
const MOMENTUM: [i32; 4] = [1, 1, -1, 0];
#[derive(Clone)]
struct Frame {
    values: [f64; 4],
    dims: [DimExpr; 4],
    physical: bool,
}
impl NativeType for Frame {
    const KIND: &'static str = "rssn/native/sim/particle-frame/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        if self.values.iter().all(|v| v.is_finite() && *v > 0.0) {
            Ok(())
        } else {
            Err(api::error(
                "particle reference values must be finite and positive",
            ))
        }
    }
}
impl Frame {
    fn dim(
        &self,
        p: [i32; 4],
    ) -> Result<DimExpr, OpError> {
        let parts = self
            .dims
            .iter()
            .zip(p)
            .filter(|(_, n)| *n != 0)
            .map(|(d, n)| d.clone().pow(Exp::Rat(Rational::int(i64::from(n)))))
            .collect();
        super::runtime_dimensions::expression(&super::runtime_dimensions::canonical(
            &DimExpr::Mul(parts),
        )?)
    }

    fn apply(
        &self,
        v: f64,
        p: [i32; 4],
        restore: bool,
        temp: bool,
    ) -> Result<f64, OpError> {
        let mut num = Vec::new();
        let mut den = Vec::new();
        for (i, n) in p.iter().enumerate() {
            for _ in 0..n.unsigned_abs() {
                if (*n > 0) == restore {
                    num.push(self.values[i]);
                } else {
                    den.push(self.values[i]);
                }
            }
        }
        if temp && self.physical {
            if restore {
                den.push(KB);
            } else {
                num.push(KB);
            }
        }
        api::scale_product(v, &num, &den)
    }

    fn temperature_dim(&self) -> DimExpr {
        if self.physical {
            DimExpr::base("Temperature")
        } else {
            DimExpr::Dimless
        }
    }

    fn scalar(
        &self,
        v: &Value,
        p: [i32; 4],
        temp: bool,
    ) -> Result<f64, OpError> {
        api::same_dimension(
            &api::dimension(v)?,
            &if temp {
                self.temperature_dim()
            } else {
                self.dim(p)?
            },
        )?;
        self.apply(super::numeric(v, 0)?, p, false, temp)
    }

    fn output(
        &self,
        v: f64,
        p: [i32; 4],
        temp: bool,
    ) -> Result<Value, OpError> {
        let dim = if temp {
            self.temperature_dim()
        } else {
            self.dim(p)?
        };
        if v.is_infinite() {
            return Ok(Value::Dyn(crate::kernel::value::DynValue {
                dim,
                value: Box::new(Value::Float(v)),
            }));
        }
        api::quantity(self.apply(v, p, true, temp)?, dim)
    }

    fn vector(
        &self,
        v: &Value,
        p: [i32; 4],
    ) -> Result<Vec<f64>, OpError> {
        let Value::Tensor(t) = v else {
            return Err(api::error("expected a measured vector"));
        };
        if t.shape != [t.data.len()] || t.data.len() > 4096 {
            return Err(api::error(
                "particle vector shape/allocation budget exceeded",
            ));
        }
        api::same_dimension(&api::dimension(v)?, &self.dim(p)?)?;
        t.data
            .iter()
            .map(|v| self.apply(super::numeric(v, 0)?, p, false, false))
            .collect()
    }

    fn array(
        &self,
        data: Vec<f64>,
        p: [i32; 4],
    ) -> Result<Value, OpError> {
        let dim = self.dim(p)?;
        Ok(Value::Tensor(TensorValue {
            shape: vec![data.len()],
            elem: TypeExpr::Float.with_dim(dim),
            data: data
                .into_iter()
                .map(|v| self.output(v, p, false))
                .collect::<Result<_, _>>()?,
        }))
    }

    fn same(
        &self,
        other: &Self,
    ) -> Result<(), OpError> {
        for i in 0..4 {
            api::same_dimension(&self.dims[i], &other.dims[i])?;
            if self.values[i].to_bits() != other.values[i].to_bits() {
                return Err(api::error(
                    "particle objects require identical physical reference scales",
                ));
            }
        }
        Ok(())
    }
}
struct Particle {
    native: md::Particle,
    frame: Frame,
}
struct Cloud {
    native: Vec<md::Particle>,
    frame: Frame,
    dimension: usize,
}
impl NativeType for Particle {
    const KIND: &'static str = "rssn/native/sim/physical-particle/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        validate_particle(&self.native)
    }
}
impl NativeType for Cloud {
    const KIND: &'static str = "rssn/native/sim/particle-cloud/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.len() > 512 || !(1..=8).contains(&self.dimension) {
            return Err(api::error("particle cloud dimension/count budget exceeded"));
        }
        for p in &self.native {
            validate_particle(p)?;
            if p.position.len() != self.dimension {
                return Err(api::error("particle cloud dimensions differ"));
            }
        }
        Ok(())
    }
}
impl NativeType for md::GaussianStream {
    const KIND: &'static str = "rssn/native/sim/gaussian-stream/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
fn validate_particle(p: &md::Particle) -> Result<(), OpError> {
    let d = p.position.len();
    if !(1..=8).contains(&d)
        || p.velocity.len() != d
        || p.force.len() != d
        || !p.mass.is_finite()
        || p.mass <= 0.0
        || !p.charge.is_finite()
        || p.position
            .iter()
            .chain(&p.velocity)
            .chain(&p.force)
            .any(|x| !x.is_finite())
    {
        return Err(api::error("particle mass/vectors/charge are invalid"));
    }
    Ok(())
}
fn frame(
    ctx: &CallCtx<'_>,
    f: &Handle<Frame>,
) -> Result<Frame, OpError> {
    f.with(ctx, |f| Ok(f.clone()))
}
fn particle_value(
    p: &md::Particle,
    f: &Frame,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("id", Value::int(p.id)),
        ("mass", f.output(p.mass, MASS, false)?),
        ("position", f.array(p.position.clone(), LENGTH)?),
        ("velocity", f.array(p.velocity.clone(), VELOCITY)?),
        ("force", f.array(p.force.clone(), FORCE)?),
        ("charge", f.output(p.charge, CHARGE, false)?),
        ("particle_type", Value::int(p.particle_type)),
    ]))
}
fn pair_result(
    result: Result<(f64, Vec<f64>), String>,
    f: &Frame,
) -> Result<Dynamic, OpError> {
    let (energy, force) = result.map_err(api::error)?;
    Ok(Dynamic(api::record(vec![
        ("energy", f.output(energy, ENERGY, false)?),
        ("force", f.array(force, FORCE)?),
    ])))
}
fn cloud_state(
    ps: &[md::Particle],
    f: &Frame,
    d: usize,
) -> Result<Value, OpError> {
    let array = |values: Vec<f64>, p: [i32; 4], matrix: bool| -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: if matrix {
                vec![ps.len(), d]
            } else {
                vec![ps.len()]
            },
            elem: TypeExpr::Float.with_dim(f.dim(p)?),
            data: values
                .into_iter()
                .map(|v| f.output(v, p, false))
                .collect::<Result<_, _>>()?,
        }))
    };
    Ok(api::record(vec![
        (
            "masses",
            array(ps.iter().map(|p| p.mass).collect(), MASS, false)?,
        ),
        (
            "positions",
            array(
                ps.iter().flat_map(|p| p.position.clone()).collect(),
                LENGTH,
                true,
            )?,
        ),
        (
            "velocities",
            array(
                ps.iter().flat_map(|p| p.velocity.clone()).collect(),
                VELOCITY,
                true,
            )?,
        ),
        (
            "forces",
            array(
                ps.iter().flat_map(|p| p.force.clone()).collect(),
                FORCE,
                true,
            )?,
        ),
        (
            "charges",
            array(ps.iter().map(|p| p.charge).collect(), CHARGE, false)?,
        ),
    ]))
}
struct ForceCallback(CallbackRef);
impl NativeArg for ForceCallback {
    fn ty() -> TypeExpr {
        let m = DimExpr::Var("M".into());
        let l = DimExpr::Var("L".into());
        let t = DimExpr::Var("T".into());
        let q = DimExpr::Var("Q".into());
        let v = DimExpr::Mul(vec![l.clone(), DimExpr::Inverse(Box::new(t.clone()))]);
        let force = DimExpr::Mul(vec![
            m.clone(),
            l.clone(),
            DimExpr::Inverse(Box::new(t.pow(2))),
        ]);
        let tensor = |rank, d| api::tensor_type(rank, TypeExpr::Float.with_dim(d));
        TypeExpr::Function(
            vec![TypeExpr::Record(vec![
                ("masses".into(), tensor(1, m)),
                ("positions".into(), tensor(2, l)),
                ("velocities".into(), tensor(2, v)),
                ("forces".into(), tensor(2, force.clone())),
                ("charges".into(), tensor(1, q)),
            ])],
            Box::new(tensor(2, force)),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Function(handle) = v else {
            return Err(api::error("expected a checked force callback"));
        };
        if ctx.callback_signature(*handle)?.params.len() != 1 {
            return Err(api::error(
                "force callback expects one particle-state record",
            ));
        }
        Ok(Self(*handle))
    }
}
#[allow(
    clippy::too_many_arguments,
    reason = "Native force ABI shares explicit state, physical frame, work counter and first-error state"
)]
fn force_call(
    ctx: &CallCtx<'_>,
    cb: &ForceCallback,
    ps: &mut [md::Particle],
    f: &Frame,
    d: usize,
    calls: &mut usize,
    budget: usize,
    failure: &mut Option<OpError>,
) -> Result<(), String> {
    let mut run = || -> Result<Vec<f64>, OpError> {
        if *calls >= budget {
            return Err(api::error("particle force callback budget exhausted"));
        }
        *calls += 1;
        let result = ctx.invoke_callback(cb.0, &[cloud_state(ps, f, d)?])?;
        let Value::Tensor(t) = &result else {
            return Err(api::error("force callback must return a force matrix"));
        };
        if t.shape != [ps.len(), d] || t.data.len() != ps.len() * d {
            return Err(api::error("force callback matrix shape mismatch"));
        }
        api::same_dimension(&api::dimension(&result)?, &f.dim(FORCE)?)?;
        t.data
            .iter()
            .map(|v| f.apply(super::numeric(v, 0)?, FORCE, false, false))
            .collect()
    };
    match run() {
        | Ok(forces) => {
            for (p, force) in ps.iter_mut().zip(forces.chunks(d)) {
                p.force.copy_from_slice(force);
            }
            Ok(())
        },
        | Err(e) => {
            let message = e.to_string();
            if failure.is_none() {
                *failure = Some(e);
            }
            Err(message)
        },
    }
}
object_binding!(frame,"num_native_md_frame","","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(mass:Input,length:Input,time:Input,charge:Input)->NewObject<Frame>,|ctx,mass,length,time,charge|{let _=ctx;let refs=[api::reference(&mass.0)?,api::reference(&length.0)?,api::reference(&time.0)?,api::reference(&charge.0)?];let expected=[DimExpr::base("Mass"),DimExpr::base("Length"),DimExpr::base("Time"),DimExpr::Mul(vec![DimExpr::base("Current"),DimExpr::base("Time")])];for(i,(_,d))in refs.iter().enumerate(){api::same_dimension(d,&expected[i])?;}Ok(NewObject::new(Frame{values:refs.clone().map(|r|r.0),dims:refs.map(|r|r.1),physical:true},expected.to_vec()))});
object_binding!(normalized_frame,"num_native_md_normalized_frame","","Explicit reduced-unit native frame; all numeric values are dimensionless and kB equals one.",false,()->NewObject<Frame>,|ctx,|{let _=ctx;Ok(NewObject::normalized(Frame{values:[1.0;4],dims:std::array::from_fn(|_|DimExpr::Dimless),physical:false}))});
object_binding!(particle_new,"num_native_md_particle_new","crate::sim::physics_md::Particle::new","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,id:usize,mass:Input,position:Input,velocity:Input)->NewObject<Particle>,|ctx,f,id,mass,position,velocity|{let _=ctx;let f=frame(ctx,&f)?;let mass=f.scalar(&mass.0,MASS,false)?;let position=f.vector(&position.0,LENGTH)?;let velocity=f.vector(&velocity.0,VELOCITY)?;if mass<=0.0||position.is_empty()||position.len()>8||position.len()!=velocity.len(){return Err(api::error("particle constructor dimension/mass mismatch"));}let dims=f.dims.to_vec();Ok(NewObject::new(Particle{native:md::Particle::new(id,mass,position,velocity),frame:f},dims))});
object_binding!(particle_with_charge,"num_native_md_particle_with_charge","crate::sim::physics_md::Particle::with_charge","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,id:usize,mass:Input,position:Input,velocity:Input,charge:Input)->NewObject<Particle>,|ctx,f,id,mass,position,velocity,charge|{let _=ctx;let f=frame(ctx,&f)?;let mass=f.scalar(&mass.0,MASS,false)?;let position=f.vector(&position.0,LENGTH)?;let velocity=f.vector(&velocity.0,VELOCITY)?;if mass<=0.0||position.is_empty()||position.len()>8||position.len()!=velocity.len(){return Err(api::error("particle constructor dimension/mass mismatch"));}let dims=f.dims.to_vec();Ok(NewObject::new(Particle{native:md::Particle::with_charge(id,mass,position,velocity,f.scalar(&charge.0,CHARGE,false)?),frame:f},dims))});
object_binding!(particle_data,"num_native_md_particle_data","","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(p:Handle<Particle>)->Dynamic,|ctx,p|{let _=ctx;p.with(ctx,|p|Ok(Dynamic(particle_value(&p.native,&p.frame)?)))});
object_binding!(particle_kinetic_energy,"num_native_md_particle_kinetic_energy","crate::sim::physics_md::Particle::kinetic_energy","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(p:Handle<Particle>)->Dynamic,|ctx,p|{let _=ctx;p.with(ctx,|p|Ok(Dynamic(p.frame.output(p.native.kinetic_energy(),ENERGY,false)?)))});
object_binding!(particle_speed,"num_native_md_particle_speed","crate::sim::physics_md::Particle::speed","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(p:Handle<Particle>)->Dynamic,|ctx,p|{let _=ctx;p.with(ctx,|p|Ok(Dynamic(p.frame.output(p.native.speed(),VELOCITY,false)?)))});
object_binding!(particle_momentum,"num_native_md_particle_momentum","crate::sim::physics_md::Particle::momentum","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(p:Handle<Particle>)->Dynamic,|ctx,p|{let _=ctx;p.with(ctx,|p|Ok(Dynamic(p.frame.array(p.native.momentum(),MOMENTUM)?)))});
object_binding!(particle_distance_to,"num_native_md_particle_distance_to","crate::sim::physics_md::Particle::distance_to","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>)->Dynamic,|ctx,a,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;Ok(Dynamic(a.frame.output(a.native.distance_to(&b.native).map_err(api::error)?,LENGTH,false)?))}))});
object_binding!(cloud,"num_native_md_cloud","","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,dimension:usize,particles:Input)->NewObject<Cloud>,|ctx,f,dimension,particles|{let _=ctx;let f=frame(ctx,&f)?;let Value::Tensor(t)=particles.0 else{return Err(api::error("expected particle handle vector"));};if t.shape!=[t.data.len()]||t.data.len()>512||!(1..=8).contains(&dimension){return Err(api::error("particle cloud count/dimension budget"));}let native=t.data.iter().map(|v|api::get::<Particle,_>(ctx,v,Particle::KIND,|p|{f.same(&p.frame)?;if p.native.position.len()!=dimension{return Err(api::error("particle cloud dimensional mismatch"));}Ok(p.native.clone())})).collect::<Result<Vec<_>,_>>()?;let dims=f.dims.to_vec();Ok(NewObject::new(Cloud{native,frame:f,dimension},dims))});
object_binding!(cloud_data,"num_native_md_cloud_data","","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->Dynamic,|ctx,cloud|{let _=ctx;cloud.with(ctx,|c|Ok(Dynamic(cloud_state(&c.native,&c.frame,c.dimension)?)))});
object_binding!(total_kinetic_energy,"num_native_md_total_kinetic_energy","crate::sim::physics_md::total_kinetic_energy","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->Dynamic,|ctx,cloud|{let _=ctx;cloud.with(ctx,|c|Ok(Dynamic(c.frame.output(md::total_kinetic_energy(&c.native),ENERGY,false)?)))});
object_binding!(temperature,"num_native_md_temperature","crate::sim::physics_md::temperature","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->Dynamic,|ctx,cloud|{let _=ctx;cloud.with(ctx,|c|Ok(Dynamic(c.frame.output(md::temperature(&c.native),ENERGY,true)?)))});
object_binding!(total_momentum,"num_native_md_total_momentum","crate::sim::physics_md::total_momentum","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->Dynamic,|ctx,cloud|{let _=ctx;cloud.with(ctx,|c|Ok(Dynamic(c.frame.array(md::total_momentum(&c.native).map_err(api::error)?,MOMENTUM)?)))});
object_binding!(center_of_mass,"num_native_md_center_of_mass","crate::sim::physics_md::center_of_mass","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->Dynamic,|ctx,cloud|{let _=ctx;cloud.with(ctx,|c|Ok(Dynamic(c.frame.array(md::center_of_mass(&c.native).map_err(api::error)?,LENGTH)?)))});
object_binding!(pressure,"num_native_md_pressure","crate::sim::physics_md::pressure","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,volume:Input,virial:Input)->Dynamic,|ctx,cloud,volume,virial|{let _=ctx;cloud.with(ctx,|c|{let mut vol=[0;4];vol[1]=c.dimension as i32;let v=c.frame.scalar(&volume.0,vol,false)?;if v<=0.0{return Err(api::error("particle volume must be positive"));}let virial=c.frame.scalar(&virial.0,ENERGY,false)?;let mut pressure=ENERGY;pressure[1]-=c.dimension as i32;Ok(Dynamic(c.frame.output(md::pressure(&c.native,v,virial),pressure,false)?))})});
object_binding!(create_cubic_lattice,"num_native_md_create_cubic_lattice","crate::sim::physics_md::create_cubic_lattice","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,n:usize,spacing:Input,mass:Input)->NewObject<Cloud>,|ctx,f,n,spacing,mass|{let _=ctx;let f=frame(ctx,&f)?;if n.checked_pow(3).and_then(|n|n.checked_mul(1)).is_none_or(|n|n>512){return Err(api::error("particle lattice allocation budget exceeded"));}let spacing=f.scalar(&spacing.0,LENGTH,false)?;let mass=f.scalar(&mass.0,MASS,false)?;if spacing<=0.0||mass<=0.0{return Err(api::error("lattice spacing and mass must be positive"));}let native=md::create_cubic_lattice(n,spacing,mass);let dims=f.dims.to_vec();Ok(NewObject::new(Cloud{native,frame:f,dimension:3},dims))});
object_binding!(create_fcc_lattice,"num_native_md_create_fcc_lattice","crate::sim::physics_md::create_fcc_lattice","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,n:usize,spacing:Input,mass:Input)->NewObject<Cloud>,|ctx,f,n,spacing,mass|{let _=ctx;let f=frame(ctx,&f)?;if n.checked_pow(3).and_then(|n|n.checked_mul(4)).is_none_or(|n|n>512){return Err(api::error("particle lattice allocation budget exceeded"));}let spacing=f.scalar(&spacing.0,LENGTH,false)?;let mass=f.scalar(&mass.0,MASS,false)?;if spacing<=0.0||mass<=0.0{return Err(api::error("lattice spacing and mass must be positive"));}let native=md::create_fcc_lattice(n,spacing,mass);let dims=f.dims.to_vec();Ok(NewObject::new(Cloud{native,frame:f,dimension:3},dims))});
object_binding!(apply_pbc,"num_native_md_apply_pbc","crate::sim::physics_md::apply_pbc","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,position:Input,box_size:Input)->Dynamic,|ctx,f,position,box_size|{let _=ctx;let f=frame(ctx,&f)?;let p=f.vector(&position.0,LENGTH)?;let size=f.vector(&box_size.0,LENGTH)?;if p.is_empty()||p.len()>8||p.len()!=size.len()||size.iter().any(|v|*v<=0.0){return Err(api::error("periodic position/box shape or extent invalid"));}Ok(Dynamic(f.array(md::apply_pbc(&p,&size),LENGTH)?))});
object_binding!(minimum_image_distance,"num_native_md_minimum_image_distance","crate::sim::physics_md::minimum_image_distance","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(f:Handle<Frame>,position:Input,box_size:Input)->Dynamic,|ctx,f,position,box_size|{let _=ctx;let f=frame(ctx,&f)?;let p=f.vector(&position.0,LENGTH)?;let size=f.vector(&box_size.0,LENGTH)?;if p.is_empty()||p.len()>8||p.len()!=size.len()||size.iter().any(|v|*v<=0.0){return Err(api::error("periodic position/box shape or extent invalid"));}Ok(Dynamic(f.array(md::minimum_image_distance(&p,&size),LENGTH)?))});
object_binding!(mean_square_displacement,"num_native_md_mean_square_displacement","crate::sim::physics_md::mean_square_displacement","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Cloud>,b:Handle<Cloud>)->Dynamic,|ctx,a,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;if a.dimension!=b.dimension||a.native.len()!=b.native.len(){return Err(api::error("MSD particle state shape mismatch"));}Ok(Dynamic(a.frame.output(md::mean_square_displacement(&a.native,&b.native),[0,2,0,0],false)?))}))});
object_binding!(radial_distribution_function,"num_native_md_radial_distribution_function","crate::sim::physics_md::radial_distribution_function","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,box_size:Input,bins:usize,radius:Input)->Dynamic,|ctx,cloud,box_size,bins,radius|{let _=ctx;cloud.with(ctx,|c|{if bins==0||bins>4096{return Err(api::error("RDF bin count must be1..4096"));}let size=c.frame.vector(&box_size.0,LENGTH)?;let radius=c.frame.scalar(&radius.0,LENGTH,false)?;if size.len()!=c.dimension||size.iter().any(|v|*v<=0.0)||radius<=0.0{return Err(api::error("RDF box/radius invalid"));}let(r,g)=md::radial_distribution_function(&c.native,&size,bins,radius);Ok(Dynamic(api::record(vec![("radius",c.frame.array(r,LENGTH)?),("density_ratio",c.frame.array(g,[0;4])?)])))})});
object_binding!(neighbour_pairs,"num_native_md_neighbour_pairs","crate::sim::physics_md::neighbour_pairs","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,box_length:Input,cutoff:Input)->Dynamic,|ctx,cloud,box_length,cutoff|{let _=ctx;cloud.with(ctx,|c|{let length=c.frame.scalar(&box_length.0,LENGTH,false)?;let cutoff=c.frame.scalar(&cutoff.0,LENGTH,false)?;if length<=0.0||cutoff<=0.0{return Err(api::error("neighbor box and cutoff must be positive"));}let pairs=md::neighbour_pairs(&c.native,length,cutoff);let n=pairs.len();Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![n,2],elem:TypeExpr::Int,data:pairs.into_iter().flat_map(|(a,b)|[Value::int(a),Value::int(b)]).collect()})))})});
object_binding!(lennard_jones_interaction,"num_native_md_lennard_jones_interaction","crate::sim::physics_md::lennard_jones_interaction","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>,epsilon:Input,sigma:Input)->Dynamic,|ctx,a,b,epsilon,sigma|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;let f=&a.frame;let epsilon=f.scalar(&epsilon.0,ENERGY,false)?;let sigma=f.scalar(&sigma.0,LENGTH,false)?;if epsilon<0.0||sigma<=0.0{return Err(api::error("LJ energy/length invalid"));}pair_result(md::lennard_jones_interaction(&a.native,&b.native,epsilon,sigma),f)}))});
object_binding!(morse_interaction,"num_native_md_morse_interaction","crate::sim::physics_md::morse_interaction","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>,energy:Input,inverse_length:Input,equilibrium:Input)->Dynamic,|ctx,a,b,energy,inverse_length,equilibrium|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;let f=&a.frame;let energy=f.scalar(&energy.0,ENERGY,false)?;let inverse_length=f.scalar(&inverse_length.0,[0,-1,0,0],false)?;let equilibrium=f.scalar(&equilibrium.0,LENGTH,false)?;if energy<0.0||inverse_length<=0.0||equilibrium<0.0{return Err(api::error("Morse parameters invalid"));}pair_result(md::morse_interaction(&a.native,&b.native,energy,inverse_length,equilibrium),f)}))});
object_binding!(harmonic_interaction,"num_native_md_harmonic_interaction","crate::sim::physics_md::harmonic_interaction","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>,stiffness:Input,equilibrium:Input)->Dynamic,|ctx,a,b,stiffness,equilibrium|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;let f=&a.frame;let stiffness=f.scalar(&stiffness.0,[1,0,-2,0],false)?;let equilibrium=f.scalar(&equilibrium.0,LENGTH,false)?;if stiffness<0.0||equilibrium<0.0{return Err(api::error("harmonic parameters invalid"));}pair_result(md::harmonic_interaction(&a.native,&b.native,stiffness,equilibrium),f)}))});
object_binding!(coulomb_interaction,"num_native_md_coulomb_interaction","crate::sim::physics_md::coulomb_interaction","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>,constant:Input)->Dynamic,|ctx,a,b,constant|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;let f=&a.frame;let constant=f.scalar(&constant.0,[1,3,-2,-2],false)?;pair_result(md::coulomb_interaction(&a.native,&b.native,constant),f)}))});
object_binding!(soft_sphere_interaction,"num_native_md_soft_sphere_interaction","crate::sim::physics_md::soft_sphere_interaction","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(a:Handle<Particle>,b:Handle<Particle>,energy:Input,length:Input,power:i32)->Dynamic,|ctx,a,b,energy,length,power|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;let f=&a.frame;let energy=f.scalar(&energy.0,ENERGY,false)?;let length=f.scalar(&length.0,LENGTH,false)?;if energy<0.0||length<=0.0||power<=0||power>256{return Err(api::error("soft-sphere parameters invalid"));}pair_result(md::soft_sphere_interaction(&a.native,&b.native,energy,length,power),f)}))});
object_binding!(remove_com_velocity,"num_native_md_remove_com_velocity","crate::sim::physics_md::remove_com_velocity","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>)->bool,|ctx,cloud|{let _=ctx;cloud.with_mut(ctx,|c|{md::remove_com_velocity(&mut c.native).map_err(api::error)?;c.validate()?;Ok(true)})});
object_binding!(initialize_velocities_maxwell_boltzmann,"num_native_md_initialize_velocities_maxwell_boltzmann","crate::sim::physics_md::initialize_velocities_maxwell_boltzmann","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,temperature:Input,seed:u64)->bool,|ctx,cloud,temperature,seed|{let _=ctx;cloud.with_mut(ctx,|c|{let temperature=c.frame.scalar(&temperature.0,ENERGY,true)?;if temperature<0.0{return Err(api::error("temperature must be nonnegative"));}md::initialize_velocities_maxwell_boltzmann(&mut c.native,temperature,seed);c.validate()?;Ok(true)})});
object_binding!(velocity_rescale,"num_native_md_velocity_rescale","crate::sim::physics_md::velocity_rescale","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,temperature:Input)->bool,|ctx,cloud,temperature|{let _=ctx;cloud.with_mut(ctx,|c|{let temperature=c.frame.scalar(&temperature.0,ENERGY,true)?;if temperature<0.0{return Err(api::error("temperature must be nonnegative"));}md::velocity_rescale(&mut c.native,temperature);c.validate()?;Ok(true)})});
object_binding!(berendsen_thermostat,"num_native_md_berendsen_thermostat","crate::sim::physics_md::berendsen_thermostat","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,temperature:Input,tau:Input,dt:Input)->bool,|ctx,cloud,temperature,tau,dt|{let _=ctx;cloud.with_mut(ctx,|c|{let temperature=c.frame.scalar(&temperature.0,ENERGY,true)?;let tau=c.frame.scalar(&tau.0,TIME,false)?;let dt=c.frame.scalar(&dt.0,TIME,false)?;if temperature<0.0||tau<=0.0||dt<=0.0||dt>tau{return Err(api::error("Berendsen temperature/time controls invalid"));}md::berendsen_thermostat(&mut c.native,temperature,tau,dt);c.validate()?;Ok(true)})});
object_binding!(gaussian_stream_new,"num_native_md_gaussian_stream_new","crate::sim::physics_md::GaussianStream::new","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(seed:u64)->NewObject<md::GaussianStream>,|ctx,seed|{let _=ctx;Ok(NewObject::normalized(md::GaussianStream::new(seed)))});
object_binding!(gaussian_stream_next,"num_native_md_gaussian_stream_next","crate::sim::physics_md::GaussianStream::next_normal","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(noise:Handle<md::GaussianStream>)->f64,|ctx,noise|{let _=ctx;noise.with_mut(ctx,|g|Ok(g.next_normal()))});
object_binding!(integrate_velocity_verlet,"num_native_md_integrate_velocity_verlet","crate::sim::physics_md::integrate_velocity_verlet","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,dt:Input,steps:usize,budget:usize,force:ForceCallback)->Dynamic,|ctx,cloud,dt,steps,budget,force|{let _=ctx;cloud.with_mut(ctx,|c|{let dt=c.frame.scalar(&dt.0,TIME,false)?;if dt<=0.0||steps>65536||budget==0||budget>1000000||steps.checked_add(1).and_then(|n|n.checked_mul(c.native.len().max(1))).and_then(|n|n.checked_mul(c.dimension)).is_none_or(|n|n>131072){return Err(api::error("Verlet step/output/work controls invalid"));}let mut calls=0;let mut failure=None;let f=c.frame.clone();let d=c.dimension;let history=md::integrate_velocity_verlet(&mut c.native,dt,steps,|ps|force_call(ctx,&force,ps,&f,d,&mut calls,budget,&mut failure));if let Some(e)=failure{return Err(e);}let history=history.map_err(api::error)?;c.validate()?;let data=history.iter().map(|ps|cloud_state(ps,&f,d)).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(api::record(vec![("states",Value::Tensor(TensorValue{shape:vec![data.len()],elem:api::dynamic_type(),data})),("force_evaluations",Value::int(calls))])))})});
object_binding!(integrate_langevin_baoab,"num_native_md_integrate_langevin_baoab","crate::sim::physics_md::integrate_langevin_baoab","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,noise:Handle<md::GaussianStream>,dt:Input,steps:usize,friction:Input,temperature:Input,budget:usize,force:ForceCallback)->Dynamic,|ctx,cloud,noise,dt,steps,friction,temperature,budget,force|{let _=ctx;cloud.with_mut(ctx,|c|{let dt=c.frame.scalar(&dt.0,TIME,false)?;let gamma=c.frame.scalar(&friction.0,[0,0,-1,0],false)?;let kt=c.frame.scalar(&temperature.0,ENERGY,true)?;if dt<=0.0||gamma<0.0||kt<0.0||steps>65536||budget==0||budget>1000000||steps.checked_mul(c.native.len().max(1)).is_none_or(|n|n>10000000){return Err(api::error("BAOAB controls/work invalid"));}let f=c.frame.clone();let d=c.dimension;let mut calls=0;let mut failure=None;let result=noise.with_mut(ctx,|g|md::integrate_langevin_baoab(&mut c.native,dt,steps,gamma,kt,g,|ps|force_call(ctx,&force,ps,&f,d,&mut calls,budget,&mut failure)).map_err(api::error));if let Some(e)=failure{return Err(e);}result?;c.validate()?;Ok(Dynamic(api::record(vec![("state",cloud_state(&c.native,&c.frame,c.dimension)?),("force_evaluations",Value::int(calls))])))})});
object_binding!(integrate_nose_hoover,"num_native_md_integrate_nose_hoover","crate::sim::physics_md::integrate_nose_hoover","Physical particle adapter. Coordinates, velocity, mass, charge, force, energy and thermal targets use checked reference units; count, shape and work bounds precede native work.",false,(cloud:Handle<Cloud>,dt:Input,steps:usize,temperature:Input,thermostat_mass:Input,xi:Input,budget:usize,force:ForceCallback)->Dynamic,|ctx,cloud,dt,steps,temperature,thermostat_mass,xi,budget,force|{let _=ctx;cloud.with_mut(ctx,|c|{let dt=c.frame.scalar(&dt.0,TIME,false)?;let kt=c.frame.scalar(&temperature.0,ENERGY,true)?;let mass=c.frame.scalar(&thermostat_mass.0,[1,2,0,0],false)?;let xi=c.frame.scalar(&xi.0,[0,0,-1,0],false)?;if dt<=0.0||kt<0.0||mass<=0.0||steps>65536||budget==0||budget>1000000||steps.checked_mul(c.native.len().max(1)).is_none_or(|n|n>10000000){return Err(api::error("Nose-Hoover controls/work invalid"));}let f=c.frame.clone();let d=c.dimension;let mut calls=0;let mut failure=None;let result=md::integrate_nose_hoover(&mut c.native,dt,steps,kt,mass,xi,|ps|force_call(ctx,&force,ps,&f,d,&mut calls,budget,&mut failure));if let Some(e)=failure{return Err(e);}let(xi,eta)=result.map_err(api::error)?;c.validate()?;Ok(Dynamic(api::record(vec![("state",cloud_state(&c.native,&c.frame,c.dimension)?),("xi",c.frame.output(xi,[0,0,-1,0],false)?),("eta",Value::Float(eta)),("force_evaluations",Value::int(calls))])))})});

fn thermostat(
    f: &Frame,
    v: &Value,
) -> Result<md::EngineThermostat, OpError> {
    if matches!(v,Value::String(s)if s=="None") {
        return Ok(md::EngineThermostat::None);
    }
    let Value::String(kind) = api::field(v, "kind")? else {
        return Err(api::error("thermostat kind must be text"));
    };
    let temperature = f.scalar(api::field(v, "temperature")?, ENERGY, true)?;
    if temperature < 0.0 {
        return Err(api::error("thermostat temperature must be nonnegative"));
    }
    match kind.as_str() {
        | "Langevin" => {
            let g = f.scalar(api::field(v, "friction")?, [0, 0, -1, 0], false)?;
            if g < 0.0 {
                return Err(api::error("friction must be nonnegative"));
            }
            Ok(md::EngineThermostat::Langevin(temperature, g))
        },
        | "Berendsen" => {
            let tau = f.scalar(api::field(v, "relaxation")?, TIME, false)?;
            if tau <= 0.0 {
                return Err(api::error("thermostat relaxation must be positive"));
            }
            Ok(md::EngineThermostat::Berendsen(temperature, tau))
        },
        | "NoseHooverChain" => {
            let tau = f.scalar(api::field(v, "relaxation")?, TIME, false)?;
            let n = <usize as super::codecs::Arg>::decode(api::field(v, "chain_length")?, 0)?;
            if tau <= 0.0 || n == 0 || n > 16 {
                return Err(api::error("Nose Hoover relaxation/chain bounds invalid"));
            }
            Ok(md::EngineThermostat::NoseHooverChain(temperature, tau, n))
        },
        | _ => Err(api::error("unknown thermostat kind")),
    }
}
object_binding!(run_engine_md,"num_native_md_run_engine_md","crate::sim::physics_md::run_engine_md","Run the real native MD engine. Physical thermal targets are Kelvin in SI frames and reduced kBT in normalized frames; LJ energy/length, time, box, count and work bounds are validated.",false,(cloud:Handle<Cloud>,box_size:Input,epsilon:Input,sigma:Input,cutoff:Input,dt:Input,steps:usize,control:Input)->Dynamic,|ctx,cloud,box_size,epsilon,sigma,cutoff,dt,steps,control|cloud.with_mut(ctx,|c|{if c.dimension>3||steps>65536||c.native.len().checked_mul(c.native.len()).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100000000){return Err(api::error("MD engine dimension/pair-step work budget exceeded"));}let size=c.frame.vector(&box_size.0,LENGTH)?;let epsilon=c.frame.scalar(&epsilon.0,ENERGY,false)?;let sigma=c.frame.scalar(&sigma.0,LENGTH,false)?;let cutoff=c.frame.scalar(&cutoff.0,LENGTH,false)?;let dt=c.frame.scalar(&dt.0,TIME,false)?;if size.len()!=c.dimension||size.iter().any(|v|*v<=0.0)||epsilon<0.0||sigma<=0.0||cutoff<=0.0||dt<=0.0{return Err(api::error("MD engine physical controls invalid"));}let t=thermostat(&c.frame,&control.0)?;let energy=md::run_engine_md(&mut c.native,&size,(epsilon,sigma,cutoff),dt,steps,t).map_err(api::error)?;c.validate()?;Ok(Dynamic(c.frame.output(energy,ENERGY,false)?))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    engine::install(i)?;
    dynamics::install(i)?;
    spin::install(i)?;
    multigrid::install(i)?;
    classical::install(i)?;
    wave::install(i)?;
    api::install(
        i,
        &[
            &frame::BINDING,
            &normalized_frame::BINDING,
            &particle_new::BINDING,
            &particle_with_charge::BINDING,
            &particle_data::BINDING,
            &particle_kinetic_energy::BINDING,
            &particle_speed::BINDING,
            &particle_momentum::BINDING,
            &particle_distance_to::BINDING,
            &cloud::BINDING,
            &cloud_data::BINDING,
            &total_kinetic_energy::BINDING,
            &temperature::BINDING,
            &total_momentum::BINDING,
            &center_of_mass::BINDING,
            &pressure::BINDING,
            &create_cubic_lattice::BINDING,
            &create_fcc_lattice::BINDING,
            &apply_pbc::BINDING,
            &minimum_image_distance::BINDING,
            &mean_square_displacement::BINDING,
            &radial_distribution_function::BINDING,
            &neighbour_pairs::BINDING,
            &lennard_jones_interaction::BINDING,
            &morse_interaction::BINDING,
            &harmonic_interaction::BINDING,
            &coulomb_interaction::BINDING,
            &soft_sphere_interaction::BINDING,
            &remove_com_velocity::BINDING,
            &initialize_velocities_maxwell_boltzmann::BINDING,
            &velocity_rescale::BINDING,
            &berendsen_thermostat::BINDING,
            &gaussian_stream_new::BINDING,
            &gaussian_stream_next::BINDING,
            &integrate_velocity_verlet::BINDING,
            &integrate_langevin_baoab::BINDING,
            &integrate_nose_hoover::BINDING,
            &run_engine_md::BINDING,
        ],
        "numerical_vector",
    )
}
