//! Geometrized native relativity with explicit physical length and time restoration.
mod models;
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::Maybe;
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
use crate::kernel::value::TensorValue;
use crate::sim::kerr;
use crate::sim::relativity as r;
const C: f64 = 299_792_458.;
#[derive(Clone)]
struct Frame {
    length: f64,
    physical: bool,
}
impl Frame {
    fn dim(
        &self,
        p: i32,
    ) -> DimExpr {
        if !self.physical || p == 0 {
            DimExpr::Dimless
        } else {
            DimExpr::base("Length").pow(Exp::Rat(Rational::int(i64::from(p))))
        }
    }

    fn read(
        &self,
        v: &Value,
        p: i32,
    ) -> Result<f64, OpError> {
        api::same_dimension(&api::dimension(v)?, &self.dim(p))?;
        let refs = vec![self.length; p.unsigned_abs() as usize];
        api::scale_product(
            super::numeric(v, 0)?,
            if p < 0 { &refs } else { &[] },
            if p > 0 { &refs } else { &[] },
        )
    }

    fn out(
        &self,
        x: f64,
        p: i32,
    ) -> Result<Value, OpError> {
        let refs = vec![self.length; p.unsigned_abs() as usize];
        api::quantity(
            api::scale_product(
                x,
                if p > 0 { &refs } else { &[] },
                if p < 0 { &refs } else { &[] },
            )?,
            self.dim(p),
        )
    }

    fn time(
        &self,
        x: f64,
    ) -> Result<Value, OpError> {
        if self.physical {
            api::quantity(
                api::scale_product(x, &[self.length], &[C])?,
                DimExpr::base("Time"),
            )
        } else {
            Ok(Value::Float(x))
        }
    }

    fn time_in(
        &self,
        x: &Value,
    ) -> Result<f64, OpError> {
        api::same_dimension(
            &api::dimension(x)?,
            &if self.physical {
                DimExpr::base("Time")
            } else {
                DimExpr::Dimless
            },
        )?;
        api::scale_product(
            super::numeric(x, 0)?,
            if self.physical { &[C] } else { &[1.] },
            &[self.length],
        )
    }

    fn states(
        &self,
        state: (f64, f64, f64, f64),
    ) -> Result<Dynamic, OpError> {
        Ok(Dynamic(api::record(vec![
            ("time", self.time(state.0)?),
            ("radius", self.out(state.1, 1)?),
            ("theta", api::quantity(state.2, DimExpr::Dimless)?),
            ("phi", api::quantity(state.3, DimExpr::Dimless)?),
        ])))
    }
}
struct Metric {
    native: r::KerrNewman,
    frame: Frame,
}
struct Geodesic {
    native: r::Geodesic,
    frame: Frame,
}
struct Orbit {
    native: kerr::KerrOrbit,
    frame: Frame,
}
struct Unbound {
    native: kerr::UnboundGeodesic,
    frame: Frame,
}
struct Schild {
    native: kerr::KerrSchildGeodesic,
    frame: Frame,
}
macro_rules! native {
    ($t:ty,$kind:literal,$body:expr) => {
        impl NativeType for $t {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = 16384;

            fn validate(&self) -> Result<(), OpError> {
                ($body)(self)
            }
        }
    };
}
native!(Frame, "rssn/native/relativity-frame/v1", |f: &Frame| {
    if f.length > 0.0 && f.length.is_finite() {
        Ok(())
    } else {
        Err(api::error("length reference must be finite and positive"))
    }
});
native!(Metric, "rssn/native/metric/v1", |f: &Metric| {
    if f.native.m > 0.0
        && [f.native.m, f.native.a, f.native.q]
            .iter()
            .all(|x| x.is_finite())
    {
        Ok(())
    } else {
        Err(api::error("metric parameters invalid"))
    }
});
native!(Geodesic, "rssn/native/geodesic/v1", |f: &Geodesic| {
    if [
        f.native.t,
        f.native.r,
        f.native.theta,
        f.native.phi,
        f.native.pr,
        f.native.ptheta,
        f.native.energy,
        f.native.lz,
        f.native.lambda,
    ]
    .iter()
    .all(|x| x.is_finite())
    {
        Ok(())
    } else {
        Err(api::error("geodesic state became nonfinite"))
    }
});
native!(Orbit, "rssn/native/kerr-orbit/v1", |_: &Orbit| Ok(()));
native!(Unbound, "rssn/native/unbound-geodesic/v1", |_: &Unbound| {
    Ok(())
});
native!(Schild, "rssn/native/kerr-schild/v1", |_: &Schild| Ok(()));
fn frame(
    ctx: &CallCtx<'_>,
    f: &Handle<Frame>,
) -> Result<Frame, OpError> {
    f.with(ctx, |f| Ok(f.clone()))
}
fn dims(f: &Frame) -> Vec<DimExpr> {
    vec![f.dim(1)]
}
fn geo(g: &Geodesic) -> Result<Dynamic, OpError> {
    let n = &g.native;
    let f = &g.frame;
    Ok(Dynamic(api::record(vec![
        ("time", f.time(n.t)?),
        ("radius", f.out(n.r, 1)?),
        ("theta", Value::Float(n.theta)),
        ("phi", Value::Float(n.phi)),
        ("radial_momentum", Value::Float(n.pr)),
        ("polar_momentum", f.out(n.ptheta, 1)?),
        ("energy_per_rest_mass", Value::Float(n.energy)),
        ("angular_momentum_per_rest_mass", f.out(n.lz, 1)?),
        ("proper_time", f.time(n.lambda)?),
        ("constraint", Value::Float(n.constraint())),
    ])))
}
fn field(
    v: &Input,
    name: &str,
) -> Result<f64, OpError> {
    let v = api::field(&v.0, name)?;
    api::normalized(v)?;
    super::numeric(v, 0)
}
fn ray(x: kerr::RayOutcome) -> Dynamic {
    Dynamic(match x {
        | kerr::RayOutcome::Captured => {
            api::record(vec![("status", Value::String("Captured".into()))])
        },
        | kerr::RayOutcome::Escaped { theta, phi } => {
            api::record(vec![
                ("status", Value::String("Escaped".into())),
                ("theta", Value::Float(theta)),
                ("phi", Value::Float(phi)),
            ])
        },
    })
}
fn scatter(x: kerr::ScatterOutcome) -> Dynamic {
    Dynamic(match x {
        | kerr::ScatterOutcome::Captured => {
            api::record(vec![("status", Value::String("Captured".into()))])
        },
        | kerr::ScatterOutcome::Escaped { azimuth, deflection } => {
            api::record(vec![
                ("status", Value::String("Escaped".into())),
                ("azimuth", Value::Float(azimuth)),
                ("deflection", Value::Float(deflection)),
            ])
        },
    })
}
fn curve(
    data: Vec<(f64, f64)>,
    f: &Frame,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![data.len(), 2],
        elem: TypeExpr::Float.with_dim(f.dim(1)),
        data: data
            .into_iter()
            .flat_map(|p| [p.0, p.1])
            .map(|v| f.out(v, 1))
            .collect::<Result<_, _>>()?,
    })))
}
fn tolerance(v: f64) -> Result<(), OpError> {
    if v <= 0.0 || v > 0.1 || !v.is_finite() {
        Err(api::error("relative tolerance must be finite in(0,0.1]"))
    } else {
        Ok(())
    }
}
object_binding!(frame,"num_native_relativity_frame","","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(reference:Input)->NewObject<Frame>,|ctx,reference|{let _=ctx;let(length,d)=api::reference(&reference.0)?;api::same_dimension(&d,&DimExpr::base("Length"))?;Ok(NewObject::new(Frame{length,physical:true},vec![d]))});
object_binding!(frame_normalized,"num_native_relativity_frame_normalized","","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,()->NewObject<Frame>,|ctx,|{let _=ctx;Ok(NewObject::normalized(Frame{length:1.,physical:false}))});
object_binding!(metric,"num_native_relativity_metric","","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass_geometric:Input,spin:Input,charge_geometric:Input)->NewObject<Metric>,|ctx,f,mass_geometric,spin,charge_geometric|{let _=ctx;let f=frame(ctx,&f)?;let native=r::KerrNewman{m:f.read(&mass_geometric.0,1)?,a:f.read(&spin.0,1)?,q:f.read(&charge_geometric.0,1)?};let d=dims(&f);Ok(NewObject::new(Metric{native,frame:f},d))});
object_binding!(schwarzschild,"num_native_relativity_schwarzschild","crate::sim::relativity::KerrNewman::schwarzschild","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass_geometric:Input)->NewObject<Metric>,|ctx,f,mass_geometric|{let _=ctx;let f=frame(ctx,&f)?;let native=r::KerrNewman::schwarzschild(f.read(&mass_geometric.0,1)?);let d=dims(&f);Ok(NewObject::new(Metric{native,frame:f},d))});
object_binding!(horizon,"num_native_relativity_horizon","crate::sim::relativity::KerrNewman::horizon","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>)->Maybe<Dynamic>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.native.horizon().map(|x|a.frame.out(x,1).map(Dynamic)).transpose()?)))});
object_binding!(isco,"num_native_relativity_isco","crate::sim::relativity::KerrNewman::isco","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,prograde:bool)->Maybe<Dynamic>,|ctx,a,prograde|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.native.isco(prograde).map(|x|a.frame.out(x,1).map(Dynamic)).transpose()?)))});
object_binding!(photon_orbit,"num_native_relativity_photon_orbit","crate::sim::relativity::KerrNewman::photon_orbit","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,prograde:bool)->Maybe<Dynamic>,|ctx,a,prograde|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.native.photon_orbit(prograde).map(|x|a.frame.out(x,1).map(Dynamic)).transpose()?)))});
object_binding!(circular_orbit,"num_native_relativity_circular_orbit","crate::sim::relativity::KerrNewman::circular_orbit","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,radius:Input,prograde:bool)->Maybe<Dynamic>,|ctx,a,radius,prograde|{let _=ctx;a.with(ctx,|a|{let r=a.frame.read(&radius.0,1)?;Ok(Maybe(a.native.circular_orbit(r,prograde).map(|(energy,angular,frequency)|Ok::<_,OpError>(Dynamic(api::record(vec![("energy_per_rest_mass",Value::Float(energy)),("angular_momentum_per_rest_mass",a.frame.out(angular,1)?),("angular_frequency",if a.frame.physical{api::quantity(api::scale_product(frequency,&[C],&[a.frame.length])?,DimExpr::base("Time").pow(-1))?}else{Value::Float(frequency)})])))).transpose()?))})});
object_binding!(hamiltonian,"num_native_relativity_hamiltonian","crate::sim::relativity::KerrNewman::hamiltonian","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,radius:Input,theta:f64,pr:f64,ptheta:Input,energy:f64,angular:Input)->f64,|ctx,a,radius,theta,pr,ptheta,energy,angular|{let _=ctx;a.with(ctx,|a|Ok(a.native.hamiltonian(a.frame.read(&radius.0,1)?,theta,pr,a.frame.read(&ptheta.0,1)?,energy,a.frame.read(&angular.0,1)?)))});
object_binding!(geodesic_new,"num_native_relativity_geodesic_new","crate::sim::relativity::Geodesic::new","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,radius:Input,theta:f64,phi:f64,energy:f64,angular:Input,ptheta:Input,radial_sign:f64,mu2:f64)->NewObject<Geodesic>,|ctx,a,radius,theta,phi,energy,angular,ptheta,radial_sign,mu2|{let _=ctx;a.with(ctx,|a|{let native=r::Geodesic::new(a.native,a.frame.read(&radius.0,1)?,theta,phi,energy,a.frame.read(&angular.0,1)?,a.frame.read(&ptheta.0,1)?,radial_sign,mu2).map_err(api::error)?;Ok(NewObject::new(Geodesic{native,frame:a.frame.clone()},dims(&a.frame)))})});
object_binding!(geodesic_from_init,"num_native_relativity_geodesic_from_init","crate::sim::relativity::Geodesic::from_init","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,init:Input)->NewObject<Geodesic>,|ctx,a,init|{let _=ctx;a.with(ctx,|a|{let f=&a.frame;let init=r::GeodesicInit{r:f.read(api::field(&init.0,"radius")?,1)?,theta:field(&init,"theta")?,phi:field(&init,"phi")?,energy:field(&init,"energy_per_rest_mass")?,lz:f.read(api::field(&init.0,"angular_momentum_per_rest_mass")?,1)?,ptheta:f.read(api::field(&init.0,"polar_momentum")?,1)?,radial_sign:field(&init,"radial_sign")?,mu2:field(&init,"mass_shell")?};let native=r::Geodesic::from_init(a.native,&init).map_err(api::error)?;Ok(NewObject::new(Geodesic{native,frame:f.clone()},dims(f)))})});
object_binding!(geodesic_circular,"num_native_relativity_geodesic_circular","crate::sim::relativity::Geodesic::circular","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Metric>,radius:Input,prograde:bool)->NewObject<Geodesic>,|ctx,a,radius,prograde|{let _=ctx;a.with(ctx,|a|Ok(NewObject::new(Geodesic{native:r::Geodesic::circular(a.native,a.frame.read(&radius.0,1)?,prograde).map_err(api::error)?,frame:a.frame.clone()},dims(&a.frame))))});
object_binding!(geodesic_state,"num_native_relativity_geodesic_state","","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Geodesic>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,geo)});
object_binding!(geodesic_constraint,"num_native_relativity_geodesic_constraint","crate::sim::relativity::Geodesic::constraint","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Geodesic>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.native.constraint()))});
object_binding!(geodesic_step,"num_native_relativity_geodesic_step","crate::sim::relativity::Geodesic::step","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Geodesic>,h:Input,project:bool)->Dynamic,|ctx,a,h,project|{let _=ctx;a.with_mut(ctx,|a|{let h=a.frame.time_in(&h.0)?;if h==0.0{return Err(api::error("zero proper-time step"));}a.native.step(h,project).map_err(api::error)?;a.validate()?;geo(a)})});
object_binding!(geodesic_integrate,"num_native_relativity_geodesic_integrate","crate::sim::relativity::Geodesic::integrate","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Geodesic>,h:Input,steps:usize,project:bool,r_stop:Input,r_max:Input)->Dynamic,|ctx,a,h,steps,project,r_stop,r_max|{let _=ctx;a.with_mut(ctx,|a|{if steps>65536{return Err(api::error("geodesic step budget exceeded"));}let h=a.frame.time_in(&h.0)?;let lo=a.frame.read(&r_stop.0,1)?;let hi=a.frame.read(&r_max.0,1)?;if h==0.0||lo<=0.0||hi<=lo{return Err(api::error("geodesic controls invalid"));}let status=a.native.integrate(h,steps,project,lo,hi).map_err(api::error)?;a.validate()?;let data=geo(a)?.0;Ok(Dynamic(api::record(vec![("termination",Value::String(format!("{status:?}"))),("state",data)])))})});
object_binding!(orbit_new,"num_native_relativity_orbit_new","crate::sim::kerr::KerrOrbit::new","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,energy:f64,angular:Input,carter:Input,radius:Input,theta:f64,phi:f64,radial_out:bool,polar_down:bool)->NewObject<Orbit>,|ctx,f,mass,spin,energy,angular,carter,radius,theta,phi,radial_out,polar_down|{let _=ctx;let f=frame(ctx,&f)?;let native=kerr::KerrOrbit::new(f.read(&mass.0,1)?,f.read(&spin.0,1)?,(energy,f.read(&angular.0,1)?,f.read(&carter.0,2)?),(f.read(&radius.0,1)?,theta,phi),(radial_out,polar_down)).map_err(api::error)?;let d=dims(&f);Ok(NewObject::new(Orbit{native,frame:f},d))});
object_binding!(orbit_from_apsides,"num_native_relativity_orbit_from_apsides","crate::sim::kerr::KerrOrbit::from_apsides","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,r_peri:Input,r_apo:Input,z_minus:f64,prograde:bool)->NewObject<Orbit>,|ctx,f,mass,spin,r_peri,r_apo,z_minus,prograde|{let _=ctx;let f=frame(ctx,&f)?;let native=kerr::KerrOrbit::from_apsides(f.read(&mass.0,1)?,f.read(&spin.0,1)?,(f.read(&r_peri.0,1)?,f.read(&r_apo.0,1)?),z_minus,prograde).map_err(api::error)?;let d=dims(&f);Ok(NewObject::new(Orbit{native,frame:f},d))});
object_binding!(orbit_radius,"num_native_relativity_orbit_radius","crate::sim::kerr::KerrOrbit::radius","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>,mino_time:Input)->Dynamic,|ctx,a,mino_time|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(a.frame.out(a.native.radius(a.frame.read(&mino_time.0,-1)?),1)?)))});
object_binding!(orbit_cos_theta,"num_native_relativity_orbit_cos_theta","crate::sim::kerr::KerrOrbit::cos_theta","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>,mino_time:Input)->f64,|ctx,a,mino_time|{let _=ctx;a.with(ctx,|a|Ok(a.native.cos_theta(a.frame.read(&mino_time.0,-1)?)))});
object_binding!(orbit_state,"num_native_relativity_orbit_state","crate::sim::kerr::KerrOrbit::state","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>,mino_time:Input)->Dynamic,|ctx,a,mino_time|{let _=ctx;a.with(ctx,|a|a.frame.states(a.native.state(a.frame.read(&mino_time.0,-1)?)))});
object_binding!(orbit_radial_period,"num_native_relativity_orbit_radial_period","crate::sim::kerr::KerrOrbit::radial_period","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(a.frame.out(a.native.radial_period(),-1)?)))});
object_binding!(orbit_polar_period,"num_native_relativity_orbit_polar_period","crate::sim::kerr::KerrOrbit::polar_period","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>)->Maybe<Dynamic>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.native.polar_period().map(|v|a.frame.out(v,-1).map(Dynamic)).transpose()?)))});
object_binding!(orbit_frequencies,"num_native_relativity_orbit_frequencies","crate::sim::kerr::KerrOrbit::frequencies","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Orbit>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let(u,o)=a.native.frequencies();let f=&a.frame;let frequencies=o.into_iter().map(|v|if f.physical{api::quantity(api::scale_product(v,&[C],&[f.length])?,DimExpr::base("Time").pow(-1))}else{Ok(Value::Float(v))}).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(api::record(vec![("upsilon_r",f.out(u[0],1)?),("upsilon_theta",f.out(u[1],1)?),("upsilon_phi",f.out(u[2],1)?),("gamma",f.out(u[3],2)?),("coordinate_frequencies",Value::Tensor(TensorValue{shape:vec![3],elem:api::dynamic_type(),data:frequencies}))])))})});
object_binding!(unbound_new,"num_native_relativity_unbound_new","crate::sim::kerr::UnboundGeodesic::new","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,mu2:f64,energy:f64,angular:Input,carter:Input,radius:Input,theta:f64,phi:f64,radial_out:bool,polar_down:bool)->NewObject<Unbound>,|ctx,f,mass,spin,mu2,energy,angular,carter,radius,theta,phi,radial_out,polar_down|{let _=ctx;let f=frame(ctx,&f)?;let native=kerr::UnboundGeodesic::new(f.read(&mass.0,1)?,f.read(&spin.0,1)?,mu2,(energy,f.read(&angular.0,1)?,f.read(&carter.0,2)?),(f.read(&radius.0,1)?,theta,phi),(radial_out,polar_down)).map_err(api::error)?;let d=dims(&f);Ok(NewObject::new(Unbound{native,frame:f},d))});
object_binding!(unbound_advance,"num_native_relativity_unbound_advance","crate::sim::kerr::UnboundGeodesic::advance","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>,delta:Input,tol:f64)->bool,|ctx,a,delta,tol|{let _=ctx;a.with_mut(ctx,|a|{tolerance(tol)?;let h=a.frame.read(&delta.0,-1)?;if h==0.0||h.abs()>1e6{return Err(api::error("geodesic advance size budget exceeded"));}Ok(a.native.advance(h,tol))})});
object_binding!(schild_new,"num_native_relativity_schild_new","crate::sim::kerr::KerrSchildGeodesic::from_constants","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,mu2:f64,energy:f64,angular:Input,carter:Input,radius:Input,theta:f64,phi:f64,radial_out:bool,polar_down:bool)->NewObject<Schild>,|ctx,f,mass,spin,mu2,energy,angular,carter,radius,theta,phi,radial_out,polar_down|{let _=ctx;let f=frame(ctx,&f)?;let native=kerr::KerrSchildGeodesic::from_constants(f.read(&mass.0,1)?,f.read(&spin.0,1)?,(energy,f.read(&angular.0,1)?,f.read(&carter.0,2)?),(f.read(&radius.0,1)?,theta,phi),mu2,(radial_out,polar_down)).map_err(api::error)?;let d=dims(&f);Ok(NewObject::new(Schild{native,frame:f},d))});
object_binding!(schild_advance,"num_native_relativity_schild_advance","crate::sim::kerr::KerrSchildGeodesic::advance","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>,delta:Input,tol:f64,r_stop:Input)->bool,|ctx,a,delta,tol,r_stop|{let _=ctx;a.with_mut(ctx,|a|{tolerance(tol)?;let h=a.frame.time_in(&delta.0)?;if h==0.0||h.abs()>1e6{return Err(api::error("geodesic advance size budget exceeded"));}Ok(a.native.advance(h,tol,a.frame.read(&r_stop.0,1)?))})});
object_binding!(unbound_state,"num_native_relativity_unbound_state","crate::sim::kerr::UnboundGeodesic::state","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{a.frame.states(a.native.state())})});
object_binding!(unbound_impact_parameter,"num_native_relativity_unbound_impact_parameter","crate::sim::kerr::UnboundGeodesic::impact_parameter","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>)->Maybe<Dynamic>,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(Maybe(a.native.impact_parameter().map(|v|a.frame.out(v,1).map(Dynamic)).transpose()?))})});
object_binding!(unbound_carter_constant,"num_native_relativity_unbound_carter_constant","crate::sim::kerr::UnboundGeodesic::carter_constant","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(Dynamic(a.frame.out(a.native.carter_constant(),2)?))})});
object_binding!(unbound_normalisation,"num_native_relativity_unbound_normalisation","crate::sim::kerr::UnboundGeodesic::normalisation","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(a.native.normalisation())})});
object_binding!(unbound_scatter,"num_native_relativity_unbound_scatter","crate::sim::kerr::UnboundGeodesic::scatter","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Unbound>,tol:f64)->Dynamic,|ctx,a,tol|{let _=ctx;a.with_mut(ctx,|a|{tolerance(tol)?;Ok(scatter(a.native.scatter(tol).map_err(api::error)?))})});
object_binding!(schild_radius,"num_native_relativity_schild_radius","crate::sim::kerr::KerrSchildGeodesic::radius","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(a.frame.out(a.native.radius(),1)?)))});
object_binding!(schild_angular_momentum,"num_native_relativity_schild_angular_momentum","crate::sim::kerr::KerrSchildGeodesic::angular_momentum","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(a.frame.out(a.native.angular_momentum(),1)?)))});
object_binding!(schild_carter_constant,"num_native_relativity_schild_carter_constant","crate::sim::kerr::KerrSchildGeodesic::carter_constant","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(a.frame.out(a.native.carter_constant(),2)?)))});
object_binding!(schild_hamiltonian,"num_native_relativity_schild_hamiltonian","crate::sim::kerr::KerrSchildGeodesic::hamiltonian","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.native.hamiltonian()))});
object_binding!(schild_energy,"num_native_relativity_schild_energy","crate::sim::kerr::KerrSchildGeodesic::energy","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.native.energy()))});
object_binding!(schild_boyer_lindquist,"num_native_relativity_schild_boyer_lindquist","crate::sim::kerr::KerrSchildGeodesic::boyer_lindquist","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(a:Handle<Schild>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let(r,theta,phi)=a.native.boyer_lindquist();Ok(Dynamic(api::record(vec![("radius",a.frame.out(r,1)?),("theta",Value::Float(theta)),("phi",Value::Float(phi))])))})});
object_binding!(photon_radius,"num_native_relativity_photon_radius","crate::sim::kerr::kerr_photon_orbit_radius","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,prograde:bool)->Dynamic,|ctx,f,mass,spin,prograde|{let _=ctx;let f=frame(ctx,&f)?;Ok(Dynamic(f.out(kerr::kerr_photon_orbit_radius(f.read(&mass.0,1)?,f.read(&spin.0,1)?,prograde),1)?))});
object_binding!(shadow_analytic,"num_native_relativity_shadow_analytic","crate::sim::kerr::shadow_boundary_analytic","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,theta:f64,points:usize)->Dynamic,|ctx,f,mass,spin,theta,points|{let _=ctx;if points>2048{return Err(api::error("shadow curve point budget exceeded"));}let f=frame(ctx,&f)?;curve(kerr::shadow_boundary_analytic(f.read(&mass.0,1)?,f.read(&spin.0,1)?,theta,points).map_err(api::error)?,&f)});
object_binding!(trace_ray,"num_native_relativity_trace_ray","crate::sim::kerr::trace_ray","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,theta:f64,radius:Input,xi:Input,eta:Input,polar_down:bool)->Dynamic,|ctx,f,mass,spin,theta,radius,xi,eta,polar_down|{let _=ctx;let f=frame(ctx,&f)?;Ok(ray(kerr::trace_ray(f.read(&mass.0,1)?,f.read(&spin.0,1)?,theta,f.read(&radius.0,1)?,(f.read(&xi.0,1)?,f.read(&eta.0,2)?),polar_down).map_err(api::error)?))});
object_binding!(trace_shadow,"num_native_relativity_trace_shadow","crate::sim::kerr::trace_shadow","Native geometrized relativity with positive length reference. Mass, spin and charge metric parameters use geometric lengths; coordinate/proper time restore seconds, Mino time uses inverse length, specific energy and angles are dimensionless.",false,(f:Handle<Frame>,mass:Input,spin:Input,theta:f64,radius:Input,points:usize)->Dynamic,|ctx,f,mass,spin,theta,radius,points|{let _=ctx;if points>64{return Err(api::error("numerical shadow ray work budget exceeded"));}let f=frame(ctx,&f)?;curve(kerr::trace_shadow(f.read(&mass.0,1)?,f.read(&spin.0,1)?,theta,f.read(&radius.0,1)?,points).map_err(api::error)?,&f)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    models::install(i)?;
    api::install(
        i,
        &[
            &frame::BINDING,
            &frame_normalized::BINDING,
            &metric::BINDING,
            &schwarzschild::BINDING,
            &horizon::BINDING,
            &isco::BINDING,
            &photon_orbit::BINDING,
            &circular_orbit::BINDING,
            &hamiltonian::BINDING,
            &geodesic_new::BINDING,
            &geodesic_from_init::BINDING,
            &geodesic_circular::BINDING,
            &geodesic_state::BINDING,
            &geodesic_constraint::BINDING,
            &geodesic_step::BINDING,
            &geodesic_integrate::BINDING,
            &orbit_new::BINDING,
            &orbit_from_apsides::BINDING,
            &orbit_radius::BINDING,
            &orbit_cos_theta::BINDING,
            &orbit_state::BINDING,
            &orbit_radial_period::BINDING,
            &orbit_polar_period::BINDING,
            &orbit_frequencies::BINDING,
            &unbound_new::BINDING,
            &unbound_advance::BINDING,
            &schild_new::BINDING,
            &schild_advance::BINDING,
            &unbound_state::BINDING,
            &unbound_impact_parameter::BINDING,
            &unbound_carter_constant::BINDING,
            &unbound_normalisation::BINDING,
            &unbound_scatter::BINDING,
            &schild_radius::BINDING,
            &schild_angular_momentum::BINDING,
            &schild_carter_constant::BINDING,
            &schild_hamiltonian::BINDING,
            &schild_energy::BINDING,
            &schild_boyer_lindquist::BINDING,
            &photon_radius::BINDING,
            &shadow_analytic::BINDING,
            &trace_ray::BINDING,
            &trace_shadow::BINDING,
        ],
        "numerical_vector",
    )
}
