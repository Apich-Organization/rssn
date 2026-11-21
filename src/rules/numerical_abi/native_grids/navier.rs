//! Native staggered Navier–Stokes state in a mechanical reference frame.
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::TypeExpr;
use super::Value;
use super::api;
use super::err;
use super::mechanics::References;
use super::object_binding;
use crate::sim::ns_mac;
pub(super) struct Ns {
    pub(super) data: ns_mac::NsSolver,
    pub(super) frame: References,
    callbacks: Vec<CompiledSnapshot>,
}
impl NativeType for Ns {
    const KIND: &'static str = "rssn/native/navier-stokes/v1";
    const MAX_BYTES: usize = 24 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        let n = self.data.n();
        if n.iter()
            .try_fold(1usize, |a, b| a.checked_mul(*b))
            .is_none_or(|n| n > 512)
            || self.data.p.len() != n.iter().product::<usize>()
            || self
                .data
                .u
                .iter()
                .flatten()
                .chain(&self.data.p)
                .any(|v| !v.is_finite())
            || !self.data.dt.is_finite()
            || self.data.dt <= 0.0
            || !self.data.t.is_finite()
        {
            return Err(err(
                "Navier–Stokes state exceeds finite shape/work/retention bounds",
            ));
        }
        Ok(())
    }
}
pub(super) fn from_native(
    data: ns_mac::NsSolver,
    frame: References,
) -> NewObject<Ns> {
    let dimensions = frame.dims.to_vec();
    NewObject::new(
        Ns {
            data,
            frame,
            callbacks: vec![],
        },
        dimensions,
    )
}

use super::super::codecs::Arg;
use super::super::native_object_api::Maybe;
use super::boundary_models::FluidFaces;
struct Config {
    data: Option<ns_mac::NsConfig>,
    frame: References,
    callbacks: Vec<CompiledSnapshot>,
}
impl NativeType for Config {
    const KIND: &'static str = "rssn/native/ns-config/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
fn decode<T: Arg>(
    v: &Value,
    name: &str,
) -> Result<T, OpError> {
    let v = api::field(v, name)?;
    api::normalized(v)?;
    T::decode(v, 0)
}
fn vec3(v: Vec<f64>) -> Result<[f64; 3], OpError> {
    v.try_into()
        .map_err(|_| err("three normalized components required"))
}
fn triple(v: Vec<usize>) -> Result<[usize; 3], OpError> {
    v.try_into().map_err(|_| err("three indices required"))
}
const fn dim(s: &Ns) -> usize {
    s.data.spatial_dimension()
}
fn finish(s: &Ns) -> Result<(), OpError> {
    for c in &s.callbacks {
        c.finish()?;
    }
    s.validate()
}
fn measured(
    s: &Ns,
    v: f64,
    p: [i64; 4],
) -> Result<Dynamic, OpError> {
    let (scale, d) = s.frame.derived(p)?;
    Ok(Dynamic(api::quantity(v * scale, d)?))
}
fn field_result(
    s: &Ns,
    v: Vec<f64>,
    p: [i64; 4],
) -> Result<Dynamic, OpError> {
    let (scale, d) = s.frame.derived(p)?;
    Ok(Dynamic(api::record(vec![
        ("normalized", v.encode()?),
        ("reference", api::quantity(scale, d)?),
    ])))
}

struct NsDivergenceSummary(Value);

impl api::NativeOutput for NsDivergenceSummary {
    fn ty() -> TypeExpr {
        let inverse_time =
            TypeExpr::Float.with_dim(DimExpr::Base("Time".into()).pow(Exp::Rat(Rational::int(-1))));
        TypeExpr::Record(vec![
            ("max_abs".into(), inverse_time.clone()),
            ("rms".into(), inverse_time),
            ("net_outflow".into(), api::dynamic_type()),
            ("fluid_volume".into(), api::dynamic_type()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(self.0)
    }
}

object_binding!(
    ns_divergence_summary,
    "num_ns_divergence_summary",
    "crate::sim::ns_mac::NsSolver::divergence_summary",
    "Reports maximum absolute divergence and a fluid-volume-weighted RMS, both in inverse time. Net outward flux has length^rank/time units, and fluid volume has length^rank units for a 2D or 3D grid. Cut-cell norms use actual fluid fractions. Fully solid domains return zero diagnostics.",
    true,
    (state: Handle<Ns>) -> NsDivergenceSummary,
    |ctx, state| state.with(ctx, |s| {
        let summary = s.data.divergence_summary().map_err(err)?;
        let rank = dim(s) as i64;
        let (rate_scale, rate_dimension) = s.frame.derived([0, -1, 0, 0])?;
        let rate = |v| api::quantity(v * rate_scale, rate_dimension.clone());
        Ok(NsDivergenceSummary(api::record(vec![
            ("max_abs", rate(summary.max_abs)?),
            ("rms", rate(summary.rms)?),
            ("net_outflow", measured(s, summary.net_outflow, [rank, -1, 0, 0])?.0),
            ("fluid_volume", measured(s, summary.fluid_volume, [rank, 0, 0, 0])?.0),
        ])))
    })
);

fn spatial(
    c: &CompiledSnapshot,
    f: &References,
    x: [f64; 3],
    out: f64,
) -> f64 {
    c.sample(&x.map(|v| v * f.scales[0])) / out
}
fn spacetime(
    c: &CompiledSnapshot,
    f: &References,
    t: f64,
    x: [f64; 3],
    out: f64,
) -> f64 {
    c.sample(&[
        t * f.scales[1],
        x[0] * f.scales[0],
        x[1] * f.scales[0],
        x[2] * f.scales[0],
    ]) / out
}
#[allow(
    clippy::too_many_arguments,
    reason = "checks one physical callback declaration with its domain, quota and retained parents"
)]
fn optional_callback(
    ctx: &CallCtx<'_>,
    v: &Value,
    name: &str,
    refs: &References,
    output: [i64; 4],
    time: bool,
    budget: usize,
    parents: &mut Vec<Value>,
) -> Result<Option<CompiledSnapshot>, OpError> {
    let v = api::field(v, name)?;
    let Value::Variant(v) = v else {
        return Err(err(
            "callback must be None or Some(native compiled program)",
        ));
    };
    if v.tag == "None" && v.payload.is_empty() {
        return Ok(None);
    }
    if v.tag != "Some" || v.payload.len() != 1 {
        return Err(err("invalid callback option envelope"));
    }
    let mut inputs = vec![refs.dims[0].clone(); 3];
    if time {
        inputs.insert(0, refs.dims[1].clone());
    }
    let c = super::structured_ops::checked_snapshot(
        ctx,
        &Input(v.payload[0].clone()),
        budget,
        &inputs,
        &refs.derived(output)?.1,
    )?;
    parents.push(v.payload[0].clone());
    Ok(Some(c))
}
object_binding!(ns_config_new,"num_ns_config_new","crate::sim::ns_mac::NsConfig::new","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",false,(frame:Handle<References>,faces:Handle<FluidFaces>,extents:Vec<usize>,lengths:Vec<f64>,viscosity:f64,options:Input,callbacks:Input,budget:usize)->NewObject<Config>,|ctx,frame,faces,extents,lengths,viscosity,options,callbacks,budget|{frame.with(ctx,|r|faces.with(ctx,|b|{
if !(2..=3).contains(&extents.len())||extents.iter().any(|n|*n<3)||extents.iter().try_fold(1usize,|a,b|a.checked_mul(*b)).is_none_or(|n|n>512)||lengths.len()!=extents.len()||lengths.iter().any(|v|*v<=0.0)||viscosity<=0.0||budget==0||budget>8_388_608{return Err(err("NS configuration work/shape/positive coefficient bounds exceeded"));}
api::same_dimension(&b.frame.dims[0],&r.dims[0])?;api::same_dimension(&b.frame.dims[1],&r.dims[1])?;api::same_dimension(&b.frame.dims[2],&r.derived([1,-1,0,0])?.1)?;
if b.frame.length.to_bits()!=r.scales[0].to_bits()||b.frame.time.to_bits()!=r.scales[1].to_bits()||b.frame.value.to_bits()!=r.derived([1,-1,0,0])?.0.to_bits(){return Err(err("fluid BC reference must equal the mechanical velocity reference"));}
let mut cfg=ns_mac::NsConfig::new(&extents,&lengths,viscosity).map_err(err)?;cfg.bcs=b.data.clone();
cfg.rotational=decode(&options.0,"rotational")?;cfg.solver=crate::rules::numerical_abi::native_grids::structured_ops::solver(&decode::<String>(&options.0,"solver")?,decode(&options.0,"solver_tolerance")?)?;
let conv=decode::<String>(&options.0,"convection")?;cfg.convection=match conv.as_str(){"Central"=>ns_mac::Convection::Central,"Quick"=>ns_mac::Convection::Quick,"Minmod"=>ns_mac::Convection::Tvd(ns_mac::ConvLimiter::Minmod),"VanLeer"=>ns_mac::Convection::Tvd(ns_mac::ConvLimiter::VanLeer),"Superbee"=>ns_mac::Convection::Tvd(ns_mac::ConvLimiter::Superbee),"Mc"=>ns_mac::Convection::Tvd(ns_mac::ConvLimiter::Mc),_=>return Err(err("unknown NS convection"))};
let lc:f64=decode(&options.0,"les_constant")?;if !(0.0..=1.0).contains(&lc){return Err(err("LES constant must be in [0,1]"));}cfg.les=match decode::<String>(&options.0,"les")?.as_str(){"None"=>ns_mac::LesModel::None,"Smagorinsky"=>ns_mac::LesModel::Smagorinsky{cs:lc},"Wale"=>ns_mac::LesModel::Wale{cw:lc},_=>return Err(err("unknown LES model"))};
let solid:Vec<bool>=decode(&options.0,"solid")?;if !solid.is_empty(){if solid.len()!=extents.iter().product::<usize>(){return Err(err("solid mask shape mismatch"));}cfg.solid=Some(solid);}
let beta:f64=decode(&options.0,"beta_min")?;if !(0.0..=1.0).contains(&beta){return Err(err("cut-cell beta_min must be in [0,1]"));}cfg.cut_cell=ns_mac::CutCellOptions{conservative:decode(&options.0,"conservative")?,beta_min:beta};
let mut parents=vec![frame.value.clone(),faces.value.clone()];let mut retained=vec![];
let mut force=vec![];for name in ["force_x","force_y","force_z"]{force.push(optional_callback(ctx,&callbacks.0,name,r,[1,-2,0,0],true,budget,&mut parents)?);}
if force.iter().any(Option::is_some){if force.iter().any(Option::is_none){return Err(err("force needs three compiled scalar component programs"));}let c:Vec<_>=force.into_iter().flatten().collect();retained.extend(c.clone());let rf=r.clone();let scale=r.derived([1,-2,0,0])?.0;cfg.force=Some(Box::new(move|t,x|std::array::from_fn(|i|spacetime(&c[i],&rf,t,x,scale))));}
let level=optional_callback(ctx,&callbacks.0,"level_set",r,[1,0,0,0],false,budget,&mut parents)?;
let mut velocity=vec![];for name in ["body_vx","body_vy","body_vz"]{velocity.push(optional_callback(ctx,&callbacks.0,name,r,[1,-1,0,0],false,budget,&mut parents)?);}
if let Some(level)=level{if velocity.iter().any(Option::is_none){return Err(err("immersed body needs all three velocity component programs"));}let vel:Vec<_>=velocity.into_iter().flatten().collect();retained.push(level.clone());retained.extend(vel.clone());let rl=r.clone();let rv=r.clone();let scale=r.derived([1,-1,0,0])?.0;cfg.body=Some(ns_mac::ImmersedBody{level_set:Box::new(move|x|spatial(&level,&rl,x,rl.scales[0])),velocity:Box::new(move|x|std::array::from_fn(|i|spatial(&vel[i],&rv,x,scale)))});}else if velocity.iter().any(Option::is_some){return Err(err("body velocity requires a level set"));}
let temp=optional_callback(ctx,&callbacks.0,"temperature",r,[0,0,0,1],false,budget,&mut parents)?;
let density=optional_callback(ctx,&callbacks.0,"density",r,[-3,0,1,0],false,budget,&mut parents)?;
if temp.is_some()&&density.is_some(){return Err(err("density and Boussinesq thermal models are exclusive"));}
if let Some(temp)=temp{let kappa:f64=decode(&options.0,"thermal_diffusivity")?;if kappa<=0.0{return Err(err("thermal diffusivity must be positive"));}let gbeta=vec3(decode(&options.0,"buoyancy")?)?;let kinds:Vec<String>=decode(&options.0,"thermal_bcs")?;let values:Vec<f64>=decode(&options.0,"thermal_values")?;if kinds.len()!=6||values.len()!=6{return Err(err("thermal faces require six BC kinds and normalized values"));}let mut bcs=[ns_mac::ScalarBc::Adiabatic;6];for i in 0..6{bcs[i]=match kinds[i].as_str(){"Adiabatic"=>ns_mac::ScalarBc::Adiabatic,"Periodic"=>ns_mac::ScalarBc::Periodic,"Dirichlet"=>ns_mac::ScalarBc::Dirichlet(values[i]),_=>return Err(err("unknown thermal BC"))};}let rf=r.clone();retained.push(temp.clone());cfg.thermal=Some(ns_mac::Thermal{kappa,gbeta,t_ref:decode(&options.0,"temperature_reference")?,bcs,init:Box::new(move|x|spatial(&temp,&rf,x,rf.scales[3]))});}
if let Some(density)=density{let rf=r.clone();let scale=r.derived([-3,0,1,0])?.0;retained.push(density.clone());cfg.density=Some(ns_mac::Density{init:Box::new(move|x|spatial(&density,&rf,x,scale)),transport:decode(&options.0,"transport_density")?,gravity:vec3(decode(&options.0,"gravity")?)?});}
Ok(NewObject::dependent(Config{data:Some(cfg),frame:r.clone(),callbacks:retained},r.dims.to_vec(),parents))
}))});
object_binding!(ns_new,"num_ns_new","crate::sim::ns_mac::NsSolver::new","Consumes one native NS configuration and constructs its actual staggered solver. Initial velocity is three compiled physical spatial programs. All callbacks retain parent handles; callback budgets and original failures are preserved.",false,(config:Handle<Config>,dt:f64,initial:Input,budget:usize)->NewObject<Ns>,|ctx,config,dt,initial,budget|{if dt<=0.0||budget==0||budget>8_388_608{return Err(err("NS dt/callback budget invalid"));}let mut parents=vec![config.value.clone()];config.with_mut(ctx,|c|{let Value::Tensor(t)=&initial.0 else{return Err(err("initial velocity needs a vector of three compiled component programs"));};
if t.shape!=[3]||t.data.len()!=3{return Err(err("initial velocity needs three components"));}let mut init=vec![];for v in &t.data{init.push(crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&Input(v.clone()),budget,&vec![c.frame.dims[0].clone();3],&c.frame.derived([1,-1,0,0])?.1)?);parents.push(v.clone());}let cfg=c.data.take().ok_or_else(||err("NS config has already been consumed"))?;let scale=c.frame.derived([1,-1,0,0])?.0;let mut callbacks=c.callbacks.clone();callbacks.extend(init.clone());let data=ns_mac::NsSolver::new(cfg,dt,&|x|std::array::from_fn(|i|spatial(&init[i],&c.frame,x,scale))).map_err(err);for p in &callbacks{p.finish()?;}let out=Ns{data:data?,frame:c.frame.clone(),callbacks};out.validate()?;Ok(NewObject::dependent(out,c.frame.dims.to_vec(),parents))})});
object_binding!(ns_step,"num_ns_step","crate::sim::ns_mac::NsSolver::step","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->bool,|ctx,state|{state.with_mut(ctx,|s|{let result=s.data.step().map_err(err);finish(s)?;result?;Ok(true)})});
object_binding!(ns_run,"num_ns_run","crate::sim::ns_mac::NsSolver::run","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,steps:usize)->bool,|ctx,state,steps|{if steps>4096{return Err(err("NS step budget exceeded"));}state.with_mut(ctx,|s|{if s.data.n().iter().product::<usize>().checked_mul(steps).is_none_or(|n|n>1_048_576){return Err(err("NS work budget exceeded"));}let result=s.data.run(steps).map_err(err);finish(s)?;result?;Ok(true)})});
object_binding!(ns_steps,"num_ns_steps","crate::sim::ns_mac::NsSolver::steps","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->usize,|ctx,state|{state.with(ctx,|s|{Ok(s.data.steps())})});
object_binding!(ns_n,"num_ns_n","crate::sim::ns_mac::NsSolver::n","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Vec<usize>,|ctx,state|{state.with(ctx,|s|{Ok(s.data.n().to_vec())})});
object_binding!(ns_divergence,"num_ns_divergence","crate::sim::ns_mac::NsSolver::divergence","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|{field_result(s,s.data.divergence(),[0,-1,0,0])})});
object_binding!(ns_eddy_viscosity,"num_ns_eddy_viscosity","crate::sim::ns_mac::NsSolver::eddy_viscosity","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|{field_result(s,s.data.eddy_viscosity(),[2,-1,0,0])})});
object_binding!(ns_max_speed,"num_ns_max_speed","crate::sim::ns_mac::NsSolver::max_speed","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|{measured(s,s.data.max_speed(),[1,-1,0,0])})});
object_binding!(ns_mass,"num_ns_mass","crate::sim::ns_mac::NsSolver::mass","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|{measured(s,s.data.mass(),[dim(s) as i64-3,0,1,0])})});
object_binding!(ns_kinetic_energy,"num_ns_kinetic_energy","crate::sim::ns_mac::NsSolver::kinetic_energy","Native integral of speed squared divided by two, without density weighting; output has dimensions Length^(dimension+2)/Time^2. For 2D this is the native integral per omitted depth.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|{measured(s,s.data.kinetic_energy(),[dim(s) as i64+2,-2,0,0])})});
object_binding!(ns_temperature,"num_ns_temperature","crate::sim::ns_mac::NsSolver::temperature","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Maybe<Dynamic>,|ctx,state|{state.with(ctx,|s|{Ok(Maybe(s.data.temperature().map(|v|field_result(s,v.to_vec(),[0,0,0,1])).transpose()?))})});
object_binding!(ns_density,"num_ns_density","crate::sim::ns_mac::NsSolver::density","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Maybe<Dynamic>,|ctx,state|{state.with(ctx,|s|{Ok(Maybe(s.data.density().map(|v|field_result(s,v.to_vec(),[-3,0,1,0])).transpose()?))})});
object_binding!(ns_cell_size,"num_ns_cell_size","crate::sim::ns_mac::NsSolver::cell_size","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,axis:usize)->Dynamic,|ctx,state,axis|{state.with(ctx,|s|{if axis>=dim(s){return Err(err("inactive axis"));}measured(s,s.data.cell_size(axis),[1,0,0,0])})});
object_binding!(ns_cell_index,"num_ns_cell_index","crate::sim::ns_mac::NsSolver::cell_index","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,index:Vec<usize>)->usize,|ctx,state,index|{state.with(ctx,|s|{let index=triple(index)?;if index.iter().zip(s.data.n()).any(|(i,n)|*i>=n){return Err(err("cell index outside grid"));}Ok(s.data.cell_index(index))})});
object_binding!(ns_volume_fraction,"num_ns_volume_fraction","crate::sim::ns_mac::NsSolver::volume_fraction","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,index:Vec<usize>)->f64,|ctx,state,index|{state.with(ctx,|s|{let index=triple(index)?;if index.iter().zip(s.data.n()).any(|(i,n)|*i>=n){return Err(err("cell index outside grid"));}Ok(s.data.volume_fraction(index))})});
object_binding!(ns_face_aperture,"num_ns_face_aperture","crate::sim::ns_mac::NsSolver::face_aperture","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,axis:usize,index:Vec<usize>)->f64,|ctx,state,axis,index|{state.with(ctx,|s|{let index=triple(index)?;if axis>=dim(s)||index.iter().enumerate().any(|(a,i)|*i>=s.data.n()[a]+usize::from(a==axis)){return Err(err("face index outside grid"));}Ok(s.data.face_aperture(axis,index))})});
object_binding!(ns_flux_through,"num_ns_flux_through","crate::sim::ns_mac::NsSolver::flux_through","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,axis:usize,index:usize)->Dynamic,|ctx,state,axis,index|{state.with(ctx,|s|{if axis>=dim(s)||index>s.data.n()[axis]{return Err(err("face plane outside grid"));}measured(s,s.data.flux_through(axis,index),[dim(s) as i64,-1,0,0])})});
object_binding!(ns_sample_pressure,"num_ns_sample_pressure","crate::sim::ns_mac::NsSolver::sample_pressure","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,position:Vec<f64>)->Dynamic,|ctx,state,position|{state.with(ctx,|s|{let x=vec3(position)?;measured(s,s.data.sample_pressure(x),[-1,-2,1,0])})});
object_binding!(ns_sample_velocity,"num_ns_sample_velocity","crate::sim::ns_mac::NsSolver::sample_velocity","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,position:Vec<f64>,component:usize)->Dynamic,|ctx,state,position,component|{state.with(ctx,|s|{let x=vec3(position)?;if component>=dim(s){return Err(err("inactive velocity component"));}measured(s,s.data.sample_velocity(component,x),[1,-1,0,0])})});
object_binding!(ns_set_density,"num_ns_set_density","crate::sim::ns_mac::NsSolver::set_density","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,values:Vec<f64>)->bool,|ctx,state,values|{state.with_mut(ctx,|s|{s.data.set_density(&values).map_err(err)?;finish(s)?;Ok(true)})});
object_binding!(ns_set_temperature,"num_ns_set_temperature","crate::sim::ns_mac::NsSolver::set_temperature","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,values:Vec<f64>)->bool,|ctx,state,values|{state.with_mut(ctx,|s|{s.data.set_temperature(&values).map_err(err)?;finish(s)?;Ok(true)})});
object_binding!(ns_force_on_region,"num_ns_force_on_region","crate::sim::ns_mac::NsSolver::force_on_region","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>,lower:Vec<f64>,upper:Vec<f64>)->Dynamic,|ctx,state,lower,upper|{state.with(ctx,|s|{let lo=vec3(lower)?;let hi=vec3(upper)?;for a in 0..dim(s){let length=s.data.cell_size(a)*s.data.n()[a]as f64;if lo[a]<0.0||hi[a]<=lo[a]||hi[a]>length{return Err(err("force box must lie inside the native domain"));}}field_result(s,s.data.force_on_region(lo,hi).to_vec(),[dim(s)as i64-2,-2,1,0])})});
object_binding!(ns_data,"num_ns_data","","Native staggered flow in explicit Length/Time/Mass/Temperature references. Fields and coordinates are normalized in that frame; physical outputs retain their dimensions. Work, indices, callback domains, session and parent lifetime are checked.",true,(state:Handle<Ns>)->Dynamic,|ctx,state|{state.with(ctx,|s|Ok(Dynamic(api::record(vec![("pressure",field_result(s,s.data.p.clone(),[-1,-2,1,0])?.0),("velocity_x",field_result(s,s.data.u[0].clone(),[1,-1,0,0])?.0),("velocity_y",field_result(s,s.data.u[1].clone(),[1,-1,0,0])?.0),("velocity_z",field_result(s,s.data.u[2].clone(),[1,-1,0,0])?.0),("time",measured(s,s.data.t,[0,1,0,0])?.0)]))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &ns_config_new::BINDING,
    &ns_new::BINDING,
    &ns_step::BINDING,
    &ns_run::BINDING,
    &ns_steps::BINDING,
    &ns_n::BINDING,
    &ns_divergence::BINDING,
    &ns_divergence_summary::BINDING,
    &ns_eddy_viscosity::BINDING,
    &ns_max_speed::BINDING,
    &ns_mass::BINDING,
    &ns_kinetic_energy::BINDING,
    &ns_temperature::BINDING,
    &ns_density::BINDING,
    &ns_cell_size::BINDING,
    &ns_cell_index::BINDING,
    &ns_volume_fraction::BINDING,
    &ns_face_aperture::BINDING,
    &ns_flux_through::BINDING,
    &ns_sample_pressure::BINDING,
    &ns_sample_velocity::BINDING,
    &ns_set_density::BINDING,
    &ns_set_temperature::BINDING,
    &ns_force_on_region::BINDING,
    &ns_data::BINDING,
];
