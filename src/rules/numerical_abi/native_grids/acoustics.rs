//! Staggered acoustics with explicit dimensional reduction and retained compiled sources.
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Faces;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Value;
use super::api;
use super::err;
use super::mechanics::References;
use super::object_binding;
use super::structured_ops::checked_snapshot;
use crate::sim::acoustics as sound;
#[derive(Clone)]
struct Config {
    data: sound::AcousticsConfig,
    refs: References,
    faces: Faces,
}
impl NativeType for Config {
    const KIND: &'static str = "rssn/native/acoustics-config/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        self.refs.validate()?;
        for axis in 0..self.data.dim.min(3) {
            for max in [false, true] {
                let face = crate::sim::boundary::Face::of_axis(axis, max);
                if let crate::sim::boundary::BoundaryCondition::Neumann(value) =
                    self.faces.data.get(face)
                {
                    if self.faces.datums[face.index()]
                        .iter()
                        .any(|d| d.snapshot.is_some())
                    {
                        return Err(err(
                            "native acoustics supports homogeneous constant Neumann faces only",
                        ));
                    }
                    let value = value.eval(0.0, [0.0; 3]);
                    if !value.is_finite() || value.abs() > 0.0 {
                        return Err(err(
                            "native acoustics does not implement nonzero Neumann pressure gradients",
                        ));
                    }
                }
            }
        }
        let retained = self
            .faces
            .datums
            .iter()
            .flatten()
            .filter_map(|d| d.snapshot.as_ref())
            .try_fold(0usize, |n, s| n.checked_add(s.retained_bytes));
        if retained.is_none_or(|n| n > 524288) {
            return Err(err("acoustic boundary compiled storage cap exceeded"));
        }
        if self.data.dim == 0 || self.data.dim > 3 || self.data.n.iter().product::<usize>() > 16384
        {
            return Err(err("acoustics configuration cap exceeded"));
        }
        Ok(())
    }
}
struct Acoustic {
    data: Option<sound::Acoustics>,
    config: Config,
    sources: Vec<CompiledSnapshot>,
}
impl NativeType for Acoustic {
    const KIND: &'static str = "rssn/native/acoustics/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.config.validate()?;
        let retained = self
            .sources
            .iter()
            .try_fold(0usize, |n, s| n.checked_add(s.retained_bytes));
        if retained.is_none_or(|n| n > 2097152) {
            return Err(err("acoustic source compiled storage cap exceeded"));
        }
        if self.sources.len() > 64 {
            return Err(err("acoustic source cap exceeded"));
        }
        let s = self
            .data
            .as_ref()
            .ok_or_else(|| err("acoustic solver invalidated by a failed mutation"))?;
        if !s.t.is_finite()
            || !s.dt.is_finite()
            || s.dt <= 0.0
            || s.p
                .iter()
                .chain(s.v.iter().flatten())
                .any(|x| !x.is_finite())
        {
            return Err(err("acoustic state is nonfinite"));
        }
        Ok(())
    }
}
fn pressure(c: &Config) -> Result<(f64, DimExpr), OpError> {
    c.refs.derived([2 - c.data.dim as i64, -2, 1, 0])
}
fn field(
    c: &Config,
    values: Vec<f64>,
    powers: [i64; 4],
) -> Result<Value, OpError> {
    let (s, d) = c.refs.derived(powers)?;
    Ok(api::record(vec![
        ("normalized", values.encode()?),
        ("reference", api::quantity(s, d)?),
    ]))
}
#[derive(Debug)]
struct CallbackAbort;
fn sample(
    s: &CompiledSnapshot,
    x: &[f64],
    scale: f64,
) -> f64 {
    let v = s.sample(x) / scale;
    if !v.is_finite() {
        std::panic::resume_unwind(Box::new(CallbackAbort));
    }
    v
}
fn guard<T>(run: impl FnOnce() -> T) -> Result<T, OpError> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        | Ok(v) => Ok(v),
        | Err(e) if e.is::<CallbackAbort>() => Err(err("acoustics compiled callback failed")),
        | Err(e) => std::panic::resume_unwind(e),
    }
}
fn mutate<T>(
    a: &mut Acoustic,
    run: impl FnOnce(&mut sound::Acoustics) -> Result<T, OpError>,
) -> Result<T, OpError> {
    let mut state = a
        .data
        .take()
        .ok_or_else(|| err("acoustic solver invalidated by previous failure"))?;
    let result = guard(|| run(&mut state));
    for s in &a.sources {
        s.finish()?;
    }
    a.config.faces.finish_for(a.config.data.dim)?;
    let out = result??;
    a.data = Some(state);
    a.validate()?;
    Ok(out)
}
object_binding!(acoustics_config,"num_grid_acoustics_config","crate::sim::acoustics::AcousticsConfig::new","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(refs:Handle<References>,faces:Handle<Faces>,n:Vec<usize>,spacing:Vec<f64>,speed:f64,density:f64,cfl:f64,pml:usize,pml_order:f64)->NewObject<Config>,|ctx:&CallCtx<'_>,refs:Handle<References>,faces:Handle<Faces>,n:Vec<usize>,spacing:Vec<f64>,speed:f64,density:f64,cfl:f64,pml:usize,pml_order:f64|{refs.with(ctx,|r|faces.with(ctx,|f|{if n.is_empty()||n.len()>3||n.iter().any(|x|*x<4||*x>1024)||n.iter().try_fold(1usize,|a,b|a.checked_mul(*b)).is_none_or(|n|n>16384)||spacing.len()!=n.len()||spacing.iter().any(|x|*x<=0.0)||speed<=0.0||density<=0.0||cfl<=0.0||cfl>1.0||pml>512||pml_order<=0.0{return Err(err("acoustics grid/medium/CFL/PML invalid"));}let(ps,pd)=r.derived([2-n.len() as i64,-2,1,0])?;for(a,b)in f.frame.dims[..2].iter().zip(&r.dims[..2]){api::same_dimension(a,b)?;}api::same_dimension(&pd,&f.frame.dims[2])?;if f.frame.length.to_bits()!=r.scales[0].to_bits()||f.frame.time.to_bits()!=r.scales[1].to_bits()||f.frame.value.to_bits()!=ps.to_bits(){return Err(err("acoustic pressure boundary references do not match mechanical frame"));}let mut data=sound::AcousticsConfig::new(&n,&spacing,speed).map_err(err)?;data.rho=density;data.cfl=cfl;data.pml_cells=pml;data.pml_order=pml_order;data.bcs=f.data.clone();Ok(NewObject::new(Config{data,refs:r.clone(),faces:f.clone()},r.dims.to_vec()))}))});
object_binding!(acoustics_new,"num_grid_acoustics_new","crate::sim::acoustics::Acoustics::new","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(config:Handle<Config>)->NewObject<Acoustic>,|ctx:&CallCtx<'_>,config:Handle<Config>|{config.with(ctx,|c|Ok(NewObject::new(Acoustic{data:Some(sound::Acoustics::new(c.data.clone()).map_err(err)?),config:c.clone(),sources:vec![]},c.refs.dims.to_vec())))});
object_binding!(acoustics_set_medium,"num_grid_acoustics_set_medium","crate::sim::acoustics::Acoustics::set_medium","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>,speed:Input,density:Input,budget:usize)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>,speed:Input,density:Input,budget:usize|{solver.with_mut(ctx,|a|{let c=a.config.clone();let(vs,vd)=c.refs.derived([1,-1,0,0])?;let(rs,rd)=c.refs.derived([-(c.data.dim as i64),0,1,0])?;let inputs=vec![c.refs.dims[0].clone();3];let v=checked_snapshot(ctx,&speed,budget,&inputs,&vd)?;let r=checked_snapshot(ctx,&density,budget,&inputs,&rd)?;let result=mutate(a,|s|s.set_medium(|x|{let x=x.map(|v|v*c.refs.scales[0]);(sample(&v,&x,vs),sample(&r,&x,rs))}).map_err(err));v.finish()?;r.finish()?;result?;Ok(true)})});
object_binding!(acoustics_set_pressure,"num_grid_acoustics_set_pressure","crate::sim::acoustics::Acoustics::set_pressure","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>,pressure_field:Input,budget:usize)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>,pressure_field:Input,budget:usize|{solver.with_mut(ctx,|a|{let c=a.config.clone();let(ps,pd)=pressure(&c)?;let p=checked_snapshot(ctx,&pressure_field,budget,&vec![c.refs.dims[0].clone();3],&pd)?;let result=mutate(a,|s|{s.set_pressure(|x|sample(&p,&x.map(|x|x*c.refs.scales[0]),ps));Ok(())});p.finish()?;result?;Ok(true)})});
object_binding!(acoustics_source,"num_grid_acoustics_source","crate::sim::acoustics::Acoustics::add_source","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>,cell:usize,source:Input,budget:usize)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>,cell:usize,source:Input,budget:usize|{solver.with_mut(ctx,|a|{if cell>=a.config.data.n.iter().product::<usize>()||a.sources.len()>=64{return Err(err("acoustic source index/count invalid"));}let c=a.config.clone();let(ps,_)=pressure(&c)?;let(_,pd)=c.refs.derived([2-c.data.dim as i64,-3,1,0])?;let snapshot=checked_snapshot(ctx,&source,budget,&[c.refs.dims[1].clone()],&pd)?;if a.sources.iter().try_fold(snapshot.retained_bytes,|n,s|n.checked_add(s.retained_bytes)).is_none_or(|n|n>2097152){return Err(err("acoustic source storage cap exceeded"));}let source=snapshot.clone();a.data.as_mut().ok_or_else(||err("acoustic solver invalidated"))?.add_source(cell,move |t|sample(&source,&[t*c.refs.scales[1]],ps/c.refs.scales[1]));a.sources.push(snapshot);Ok(true)})});
object_binding!(acoustics_n,"num_grid_acoustics_n","crate::sim::acoustics::Acoustics::n","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>)->Vec<usize>,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>|{solver.with(ctx,|a|Ok(a.data.as_ref().ok_or_else(||err("acoustics invalidated"))?.n().to_vec()))});
object_binding!(acoustics_index,"num_grid_acoustics_index","crate::sim::acoustics::Acoustics::idx","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>,ijk:Vec<usize>)->usize,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>,ijk:Vec<usize>|{solver.with(ctx,|a|{if ijk.len()!=3||ijk.iter().zip(a.config.data.n).any(|(i,n)|*i>=n){return Err(err("acoustics cell index invalid"));}Ok(a.data.as_ref().ok_or_else(||err("acoustics invalidated"))?.idx(ijk[0],ijk[1],ijk[2]))})});
object_binding!(acoustics_energy,"num_grid_acoustics_energy","crate::sim::acoustics::Acoustics::energy","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>|{solver.with(ctx,|a|{let(s,d)=a.config.refs.derived([2,-2,1,0])?;Ok(Dynamic(api::quantity(a.data.as_ref().ok_or_else(||err("acoustics invalidated"))?.energy()*s,d)?))})});
object_binding!(acoustics_rest,"num_grid_acoustics_rest","crate::sim::acoustics::Acoustics::start_from_rest","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>|{solver.with_mut(ctx,|a|{mutate(a,|s|{s.start_from_rest();Ok(())})?;Ok(true)})});
object_binding!(acoustics_step,"num_grid_acoustics_step","crate::sim::acoustics::Acoustics::step","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>|{solver.with_mut(ctx,|a|{mutate(a,|s|{s.step();Ok(())})?;Ok(true)})});
object_binding!(acoustics_run,"num_grid_acoustics_run","crate::sim::acoustics::Acoustics::run","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>,steps:usize)->bool,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>,steps:usize|{solver.with_mut(ctx,|a|{if steps>8192||a.config.data.n.iter().product::<usize>().checked_mul(steps).is_none_or(|v|v>8388608){return Err(err("acoustic time/work budget exceeded"));}mutate(a,|s|{s.run(steps);Ok(())})?;Ok(true)})});
object_binding!(acoustics_fields,"num_grid_acoustics_fields","","Native staggered acoustic call in a mechanical reference frame. In dimension d, density is Mass/Length^d and pressure is Energy/Length^d, representing the explicitly reduced continuum. Velocity is Length/Time and total energy is Energy. All numeric coordinates, spacings, speed, density and arrays are normalized by these references. Compiled physical callbacks retain dimensions and finite evaluation budgets; failed mutations invalidate their solver handle.",true,(solver:Handle<Acoustic>)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Acoustic>|{solver.with(ctx,|a|{let c=&a.config;let s=a.data.as_ref().ok_or_else(||err("acoustics invalidated"))?;Ok(Dynamic(api::record(vec![("shape",s.n().to_vec().encode()?),("pressure",field(c,s.p.clone(),[2-c.data.dim as i64,-2,1,0])?),("velocity_x",field(c,s.v[0].clone(),[1,-1,0,0])?),("velocity_y",field(c,s.v[1].clone(),[1,-1,0,0])?),("velocity_z",field(c,s.v[2].clone(),[1,-1,0,0])?),("time",api::quantity(s.t*c.refs.scales[1],c.refs.dims[1].clone())?),("dt",api::quantity(s.dt*c.refs.scales[1],c.refs.dims[1].clone())?)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &acoustics_config::BINDING,
            &acoustics_new::BINDING,
            &acoustics_set_medium::BINDING,
            &acoustics_set_pressure::BINDING,
            &acoustics_source::BINDING,
            &acoustics_n::BINDING,
            &acoustics_index::BINDING,
            &acoustics_energy::BINDING,
            &acoustics_rest::BINDING,
            &acoustics_step::BINDING,
            &acoustics_run::BINDING,
            &acoustics_fields::BINDING,
        ],
        "numerical_vector",
    )
}
