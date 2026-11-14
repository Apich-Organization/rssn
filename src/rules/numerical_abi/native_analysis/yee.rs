//! Native Yee FDTD with physical wave speed, impedance and field references.
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
use crate::rules::numerical_abi::native_backend::CompiledSnapshot;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::Many;
use crate::rules::numerical_abi::native_object_api::Maybe;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use crate::sim::boundary::BoundaryCondition as B;
use crate::sim::boundary::FaceBcs;
use crate::sim::yee_fdtd as y;
#[derive(Clone)]
struct Frame {
    length: f64,
    time: f64,
    electric: f64,
    magnetic: f64,
    impedance: f64,
    speed: f64,
}
fn length() -> DimExpr {
    DimExpr::base("Length")
}
fn time() -> DimExpr {
    DimExpr::base("Time")
}
fn mass() -> DimExpr {
    DimExpr::base("Mass")
}
fn charge() -> DimExpr {
    DimExpr::Mul(vec![DimExpr::base("Current"), time()])
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
fn electric() -> DimExpr {
    DimExpr::Mul(vec![mass(), length(), inv(charge()), power(time(), -2)])
}
fn magnetic() -> DimExpr {
    DimExpr::Mul(vec![charge(), inv(time()), inv(length())])
}
fn impedance() -> DimExpr {
    DimExpr::Mul(vec![
        mass(),
        power(length(), 2),
        inv(time()),
        power(charge(), -2),
    ])
}
fn energy() -> DimExpr {
    DimExpr::Mul(vec![mass(), power(length(), 2), power(time(), -2)])
}
fn conductivity() -> DimExpr {
    DimExpr::Mul(vec![inv(impedance()), inv(length())])
}
fn current() -> DimExpr {
    DimExpr::Mul(vec![charge(), inv(time()), power(length(), -2)])
}
impl Frame {
    fn new(
        l: &Input,
        speed: &Input,
        z: &Input,
        e: &Input,
    ) -> Result<Self, OpError> {
        let (length, ld) = api::reference(&l.0)?;
        let (speed, sd) = api::reference(&speed.0)?;
        let (impedance, zd) = api::reference(&z.0)?;
        let (electric, ed) = api::reference(&e.0)?;
        api::same_dimension(&ld, &self::length())?;
        api::same_dimension(&sd, &DimExpr::Mul(vec![self::length(), inv(time())]))?;
        api::same_dimension(&zd, &self::impedance())?;
        api::same_dimension(&ed, &self::electric())?;
        let time = api::scale_product(length, &[], &[speed])?;
        let magnetic = api::scale_product(electric, &[], &[impedance])?;
        Ok(Self {
            length,
            time,
            electric,
            magnetic,
            impedance,
            speed,
        })
    }

    fn energy_scale(&self) -> Result<f64, OpError> {
        api::scale_product(
            1.0,
            &[
                self.electric,
                self.electric,
                self.length,
                self.length,
                self.length,
            ],
            &[self.impedance, self.speed],
        )
    }

    fn sigma(&self) -> Result<f64, OpError> {
        api::scale_product(1.0, &[], &[self.impedance, self.length])
    }

    fn current(&self) -> Result<f64, OpError> {
        api::scale_product(self.electric, &[], &[self.impedance, self.length])
    }
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn scalar(
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
fn output(
    data: &[f64],
    scale: f64,
    dim: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![data.len()],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dim),
        data: data
            .iter()
            .map(|v| api::scale_product(*v, &[scale], &[]).map(Value::Float))
            .collect::<Result<_, _>>()?,
    })))
}
fn comp(name: &str) -> Result<y::Comp, OpError> {
    Ok(match name {
        | "Ex" => y::Comp::Ex,
        | "Ey" => y::Comp::Ey,
        | "Ez" => y::Comp::Ez,
        | "Hx" => y::Comp::Hx,
        | "Hy" => y::Comp::Hy,
        | "Hz" => y::Comp::Hz,
        | _ => return Err(api::error("unknown native Yee field component")),
    })
}
const fn electric_comp(c: y::Comp) -> bool {
    matches!(c, y::Comp::Ex | y::Comp::Ey | y::Comp::Ez)
}
#[derive(Clone)]
struct Pole {
    data: y::Pole,
    time: f64,
}
impl NativeType for Pole {
    const KIND: &'static str = "rssn/native/analysis/yee/pole/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.data.omega0, self.data.gamma, self.data.strength]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(api::error(
                "Yee pole needs finite nonnegative frequency/damping/strength",
            ));
        }
        Ok(())
    }
}
struct Config {
    data: Option<y::YeeConfig>,
    frame: Frame,
    programs: Vec<CompiledSnapshot>,
}
impl NativeType for Config {
    const KIND: &'static str = "rssn/native/analysis/yee/config/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        let cfg = self
            .data
            .as_ref()
            .ok_or_else(|| api::error("native Yee configuration was consumed by Yee.new"))?;
        if cfg.n.iter().product::<usize>() > 4096 {
            return Err(api::error("Yee grid cell budget exceeds4096"));
        }
        Ok(())
    }
}
struct Yee {
    data: y::Yee,
    frame: Frame,
    programs: Vec<CompiledSnapshot>,
    sources: usize,
    snapshot_bytes: usize,
}
impl NativeType for Yee {
    const KIND: &'static str = "rssn/native/analysis/yee/solver/v1";
    const MAX_BYTES: usize = 48 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.size() > 4096
            || self.sources > 4
            || self.snapshot_bytes > 32 * 1024 * 1024
            || self.data.f.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(api::error(
                "Yee grid/source/32MiB independent snapshot budget or finite state invalid",
            ));
        }
        Ok(())
    }
}
impl Yee {
    fn finish(&self) -> Result<(), OpError> {
        for p in &self.programs {
            p.finish()?;
        }
        self.validate()
    }

    fn retain(
        &mut self,
        p: CompiledSnapshot,
    ) -> Result<(), OpError> {
        if self.sources >= 4
            || self
                .snapshot_bytes
                .checked_add(p.retained_bytes)
                .is_none_or(|n| n > 32 * 1024 * 1024)
        {
            return Err(api::error(
                "Yee retained source count4or32MiB snapshot reservation exhausted",
            ));
        }
        self.snapshot_bytes += p.retained_bytes;
        self.programs.push(p);
        self.sources += 1;
        Ok(())
    }
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn profile(
    ctx: &CallCtx<'_>,
    value: &Input,
    out_dim: Option<DimExpr>,
    maxeval: usize,
) -> Result<CompiledSnapshot, OpError> {
    let p = compiled_snapshot(ctx, &value.0, maxeval)?;
    if p.input_dimensions.len() != 3 {
        return Err(api::error("Yee material profile needs3spatial inputs"));
    }
    for d in &p.input_dimensions {
        api::same_dimension(d, &length())?;
    }
    if let Some(out_dim) = out_dim {
        api::same_dimension(&p.output_dimension, &out_dim)?;
    }
    Ok(p)
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn spatial(
    p: CompiledSnapshot,
    frame: Frame,
    scale: f64,
) -> Box<dyn Fn([f64; 3]) -> f64 + Send + Sync> {
    Box::new(move |x| {
        let input = x.map(|x| api::scale_product(x, &[frame.length], &[]).unwrap_or(f64::NAN));
        let value = p.sample(&input);
        api::scale_product(value, &[], &[scale]).unwrap_or(f64::NAN)
    })
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn waveform(
    ctx: &CallCtx<'_>,
    handle: &Input,
    frame: &Frame,
    out: DimExpr,
    maxeval: usize,
) -> Result<CompiledSnapshot, OpError> {
    let p = compiled_snapshot(ctx, &handle.0, maxeval)?;
    if p.input_dimensions.len() != 1 {
        return Err(api::error("Yee waveform requires one physical time input"));
    }
    api::same_dimension(&p.input_dimensions[0], &time())?;
    api::same_dimension(&p.output_dimension, &out)?;
    let _ = frame;
    Ok(p)
}
fn wave(
    p: CompiledSnapshot,
    frame: &Frame,
    scale: f64,
) -> Box<dyn Fn(f64) -> f64 + Send + Sync> {
    let time_scale = frame.time;
    Box::new(move |t| {
        let t = api::scale_product(t, &[time_scale], &[]).unwrap_or(f64::NAN);
        api::scale_product(p.sample(&[t]), &[], &[scale]).unwrap_or(f64::NAN)
    })
}
object_binding!(pole_lorentz,"num_analysis_yee_pole_lorentz","crate::sim::yee_fdtd::Pole::lorentz","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(delta:f64,omega:Input,gamma:Input,time_reference:Input)->NewObject<Pole>,|_ctx,delta,omega,gamma,time_reference|{let(time,dim)=api::reference(&time_reference.0)?;api::same_dimension(&dim,&super::time())?;let gamma=scalar(&gamma,api::scale_product(1.0,&[],&[time])?,inv(super::time()))?;let omega=scalar(&omega,api::scale_product(1.0,&[],&[time])?,inv(super::time()))?;let data=y::Pole::lorentz(delta,omega,gamma);Ok(NewObject::new(Pole{data,time},vec![super::time()]))});
object_binding!(pole_drude,"num_analysis_yee_pole_drude","crate::sim::yee_fdtd::Pole::drude","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(omega:Input,gamma:Input,time_reference:Input)->NewObject<Pole>,|_ctx,omega,gamma,time_reference|{let(time,dim)=api::reference(&time_reference.0)?;api::same_dimension(&dim,&super::time())?;let gamma=scalar(&gamma,api::scale_product(1.0,&[],&[time])?,inv(super::time()))?;let omega=scalar(&omega,api::scale_product(1.0,&[],&[time])?,inv(super::time()))?;let data=y::Pole::drude(omega,gamma);Ok(NewObject::new(Pole{data,time},vec![super::time()]))});
object_binding!(susceptibility,"num_analysis_yee_susceptibility","crate::sim::yee_fdtd::Pole::susceptibility","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(pole:Handle<Pole>,omega:Input)->[f64;2],|ctx,pole,omega|pole.with(ctx,|p|{let omega=scalar(&omega,api::scale_product(1.0,&[],&[p.time])?,inv(time()))?;let(a,b)=p.data.susceptibility(omega);Ok([a,b])}));
object_binding!(config_new,"num_analysis_yee_config_new","crate::sim::yee_fdtd::YeeConfig::new","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). PEC level sets may have any scalar output domain normalized by their explicit positive program reference; only signs and scale-invariant fractions are used. Grid≤4096and source/work bounds precede execution.",false,(n:Vec<usize>,spacing:Input,length_reference:Input,speed_reference:Input,impedance_reference:Input,electric_reference:Input,controls:Input,profiles:Many<Input>,poles:Many<Handle<Pole>>,maxeval:usize)->NewObject<Config>,|ctx,n,spacing,length_reference,speed_reference,impedance_reference,electric_reference,controls,profiles,poles,maxeval|{if n.is_empty()||n.len()>3||n.iter().try_fold(1usize,|a,b|a.checked_mul(*b)).is_none_or(|p|p>4096)||profiles.0.len()>3||poles.0.len()>8{return Err(api::error("Yee cell/profile/pole budget exceeded"));}let frame=Frame::new(&length_reference,&speed_reference,&impedance_reference,&electric_reference)?;api::same_dimension(&api::dimension(&spacing.0)?,&length())?;let spacing=<Vec<f64> as Arg>::decode(&spacing.0,0)?.into_iter().map(|h|api::scale_product(h,&[],&[frame.length])).collect::<Result<Vec<_>,_>>()?;let mut data=y::YeeConfig::new(&n,&spacing).map_err(api::error)?;let f=|name|{let v=api::field(&controls.0,name)?;api::normalized(v)?;<f64 as Arg>::decode(v,0)};data.cfl=f("cfl")?;data.pml_cells=<usize as Arg>::decode(api::field(&controls.0,"pml_cells")?,0)?;if data.pml_cells>64{return Err(api::error("PML layer budget exceeds64"));}data.pml_order=f("pml_order")?;data.pml_kappa_max=f("pml_kappa_max")?;data.dm_min_area=f("dm_min_area")?;data.pml_alpha_max=scalar(&Input(api::field(&controls.0,"pml_alpha_max")?.clone()),api::scale_product(1.0,&[],&[frame.time])?,inv(time()))?;let bc=<String as Arg>::decode(api::field(&controls.0,"boundary")?,0)?;data.bcs=FaceBcs::uniform(&match bc.as_str(){"PEC"=>B::dirichlet(0.0),"PMC"=>B::neumann(0.0),"Periodic"=>B::Periodic,"Absorbing"=>B::Absorbing,_=>return Err(api::error("Yee boundary must bePEC/PMC/Periodic/Absorbing"))});let mut programs=Vec::new();let mut parents=Vec::new();let mut kinds=std::collections::BTreeSet::new();for descriptor in &profiles.0{let kind=<String as Arg>::decode(api::field(&descriptor.0,"kind")?,0)?;if !kinds.insert(kind.clone()){return Err(api::error("duplicate Yee material profile"));}let value=Input(api::field(&descriptor.0,"program")?.clone());let(dim,scale)=match kind.as_str(){"eps_r"=>(Some(DimExpr::Dimless),1.0),"sigma"=>(Some(conductivity()),frame.sigma()?),"pec_body"=>(None,1.0),_=>return Err(api::error("Yee profile kind must beeps_r/sigma/pec_body"))};let p=profile(ctx,&value,dim,maxeval)?;let scale=if kind=="pec_body"{p.output_scale}else{scale};parents.push(value.0);programs.push(p.clone());let callback=spatial(p,frame.clone(),scale);match kind.as_str(){"eps_r"=>data.eps_r=Some(callback),"sigma"=>data.sigma=Some(callback),_=>data.pec_body=Some(callback)}}let mut native_poles=Vec::new();for pole in &poles.0{native_poles.push(pole.with(ctx,|p|{let ratio=api::scale_product(frame.time,&[],&[p.time])?;Ok(y::Pole{omega0:api::scale_product(p.data.omega0,&[ratio],&[])?,gamma:api::scale_product(p.data.gamma,&[ratio],&[])?,strength:api::scale_product(p.data.strength,&[ratio,ratio],&[])?})})?);}
if !native_poles.is_empty(){data.dispersive=Some(Box::new(move|_|native_poles.clone()));}
if programs.iter().map(|p|p.retained_bytes).sum::<usize>()>32*1024*1024{return Err(api::error("Yee initial snapshot reservation exceeds32MiB"));}Ok(NewObject::dependent(Config{data:Some(data),frame,programs},vec![length(),time(),electric(),magnetic()],parents))});
object_binding!(new,"num_analysis_yee_new","crate::sim::yee_fdtd::Yee::new","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(config:Handle<Config>)->NewObject<Yee>,|ctx,config|{let(frame,data,programs)=config.with_mut(ctx,|c|Ok((c.frame.clone(),c.data.take().ok_or_else(||api::error("Yee config was already consumed"))?,std::mem::take(&mut c.programs))))?;let result=y::Yee::new(data);for p in &programs{p.finish()?;}let data=result.map_err(api::error)?;let snapshot_bytes=programs.iter().map(|p|p.retained_bytes).sum();Ok(NewObject::new(Yee{data,frame,programs,sources:0,snapshot_bytes},vec![length(),time(),electric(),magnetic()]))});
object_binding!(size,"num_analysis_yee_size","crate::sim::yee_fdtd::Yee::size","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->usize,|ctx,solver|solver.with(ctx,|s|Ok(s.data.size())));
object_binding!(n,"num_analysis_yee_n","crate::sim::yee_fdtd::Yee::n","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->[usize;3],|ctx,solver|solver.with(ctx,|s|Ok(s.data.n())));
object_binding!(idx,"num_analysis_yee_idx","crate::sim::yee_fdtd::Yee::idx","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,i:usize,j:usize,k:usize)->usize,|ctx,solver,i,j,k|solver.with(ctx,|s|{let n=s.data.n();if i>=n[0]||j>=n[1]||k>=n[2]{return Err(api::error("Yee cell index out of range"));}Ok(s.data.idx(i,j,k))}));
object_binding!(h,"num_analysis_yee_h","crate::sim::yee_fdtd::Yee::h","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,axis:usize)->Dynamic,|ctx,solver,axis|solver.with(ctx,|s|{if axis>2{return Err(api::error("Yee axis must be0..2"));}Ok(Dynamic(api::quantity(api::scale_product(s.data.h(axis),&[s.frame.length],&[])?,length())?))}));
object_binding!(energy,"num_analysis_yee_energy","crate::sim::yee_fdtd::Yee::energy","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->Dynamic,|ctx,solver|solver.with(ctx,|s|Ok(Dynamic(api::quantity(api::scale_product(s.data.energy(),&[s.frame.energy_scale()?],&[])?,energy())?))));
object_binding!(ohmic_loss,"num_analysis_yee_ohmic_loss","crate::sim::yee_fdtd::Yee::ohmic_loss","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->Dynamic,|ctx,solver|solver.with(ctx,|s|Ok(Dynamic(api::quantity(api::scale_product(s.data.ohmic_loss(),&[s.frame.energy_scale()?],&[])?,energy())?))));
object_binding!(ade_dt_limit,"num_analysis_yee_ade_dt_limit","crate::sim::yee_fdtd::Yee::ade_dt_limit","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->Maybe<Dynamic>,|ctx,solver|solver.with(ctx,|s|{let limit=s.data.ade_dt_limit();if limit.is_infinite()&&limit.is_sign_positive(){return Ok(Maybe(None));}Ok(Maybe(Some(Dynamic(api::quantity(api::scale_product(limit,&[s.frame.time],&[])?,time())?))))}));
object_binding!(field,"num_analysis_yee_field","crate::sim::yee_fdtd::Yee::field","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,component:String)->Dynamic,|ctx,solver,component|solver.with(ctx,|s|{let c=comp(&component)?;output(s.data.field(c),if electric_comp(c){s.frame.electric}else{s.frame.magnetic},if electric_comp(c){electric()}else{magnetic()})}));
object_binding!(field_mut,"num_analysis_yee_field_mut","crate::sim::yee_fdtd::Yee::field_mut","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,component:String,values:Input)->(),|ctx,solver,component,values|solver.with_mut(ctx,|s|{let c=comp(&component)?;let(scale,dim)=if electric_comp(c){(s.frame.electric,electric())}else{(s.frame.magnetic,magnetic())};api::same_dimension(&api::dimension(&values.0)?,&dim)?;let values=<Vec<f64> as Arg>::decode(&values.0,0)?;if values.len()!=s.data.size(){return Err(api::error("Yee field replacement size mismatch"));}let values=values.iter().map(|v|api::scale_product(*v,&[],&[scale])).collect::<Result<Vec<_>,_>>()?;s.data.field_mut(c).copy_from_slice(&values);s.finish()}));
object_binding!(sigma,"num_analysis_yee_sigma","crate::sim::yee_fdtd::Yee::sigma","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->Dynamic,|ctx,solver|solver.with(ctx,|s|output(s.data.sigma(),s.frame.sigma()?,conductivity())));
object_binding!(step,"num_analysis_yee_step","crate::sim::yee_fdtd::Yee::step","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->(),|ctx,solver|solver.with_mut(ctx,|s|{s.data.step();s.finish()}));
object_binding!(half_step_h,"num_analysis_yee_half_step_h","crate::sim::yee_fdtd::Yee::half_step_h","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>)->(),|ctx,solver|solver.with_mut(ctx,|s|{s.data.half_step_h();s.finish()}));
object_binding!(run,"num_analysis_yee_run","crate::sim::yee_fdtd::Yee::run","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,steps:usize)->(),|ctx,solver,steps|solver.with_mut(ctx,|s|{if steps.checked_mul(s.data.size()).is_none_or(|n|n>1000000){return Err(api::error("Yee step-times-cell work budget exceeds1000000"));}s.data.run(steps);s.finish()}));
object_binding!(add_source,"num_analysis_yee_add_source","crate::sim::yee_fdtd::Yee::add_source","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,cell:[usize;3],component:String,program:Input,maxeval:usize)->(),|ctx,solver,cell,component,program,maxeval|solver.with_mut(ctx,|s|{let c=comp(&component)?;if !electric_comp(c)||cell.iter().zip(s.data.n()).any(|(i,n)|*i>=n){return Err(api::error("current source requires an electric component and valid cell"));}let p=waveform(ctx,&program,&s.frame,current(),maxeval)?;let scale=s.frame.current()?;let callback=wave(p.clone(),&s.frame,scale);s.retain(p)?;s.data.add_source(y::CurrentSource{cell,comp:c,waveform:callback});Ok(())}));
object_binding!(set_plane_wave,"num_analysis_yee_set_plane_wave","crate::sim::yee_fdtd::Yee::set_plane_wave","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,first:usize,last:usize,psi:f64,program:Input,total_steps:usize,maxeval:usize)->(),|ctx,solver,first,last,psi,program,total_steps,maxeval|solver.with_mut(ctx,|s|{if total_steps>100000||first>=last||last>=s.data.n()[0]{return Err(api::error("TF/SF bounds/incident-history budget invalid"));}let p=waveform(ctx,&program,&s.frame,electric(),maxeval)?;if s.sources>=4||s.snapshot_bytes.checked_add(p.retained_bytes).is_none_or(|n|n>32*1024*1024){return Err(api::error("Yee source snapshot reservation exhausted"));}let waveform=wave(p.clone(),&s.frame,s.frame.electric);s.data.set_plane_wave(y::PlaneWave{ia:first,ib:last,psi,waveform},total_steps).map_err(api::error)?;s.retain(p)?;s.finish()}));
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn oblique(
    ctx: &CallCtx<'_>,
    solver: Handle<Yee>,
    lo: Input,
    hi: Input,
    direction: [f64; 3],
    polarization: [f64; 3],
    origin: Input,
    program: Input,
    rotation: [[f64; 3]; 3],
    maxeval: usize,
    rotated: bool,
) -> Result<(), OpError> {
    solver.with_mut(ctx, |s| {
        let a = |v: &Input| {
            api::same_dimension(&api::dimension(&v.0)?, &length())?;
            let mut a = <[f64; 3] as Arg>::decode(&v.0, 0)?;
            for v in &mut a {
                *v = api::scale_product(*v, &[], &[s.frame.length])?;
            }
            Ok::<_, OpError>(a)
        };
        let lo = a(&lo)?;
        let hi = a(&hi)?;
        let origin = a(&origin)?;
        let p = waveform(ctx, &program, &s.frame, electric(), maxeval)?;
        if s.sources >= 4
            || s.snapshot_bytes
                .checked_add(p.retained_bytes)
                .is_none_or(|n| n > 32 * 1024 * 1024)
        {
            return Err(api::error(
                "Yee oblique source snapshot reservation exhausted",
            ));
        }
        let w = y::ObliquePlaneWave {
            lo,
            hi,
            direction,
            polarization,
            origin,
            waveform: wave(p.clone(), &s.frame, s.frame.electric),
        };
        if rotated {
            s.data.set_oblique_plane_wave_rotated(w, rotation)
        } else {
            s.data.set_oblique_plane_wave(w)
        }
        .map_err(api::error)?;
        s.retain(p)?;
        s.finish()
    })
}
object_binding!(set_oblique_plane_wave,"num_analysis_yee_set_oblique_plane_wave","crate::sim::yee_fdtd::Yee::set_oblique_plane_wave","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,lo:Input,hi:Input,direction:[f64;3],polarization:[f64;3],origin:Input,program:Input,rotation:[[f64;3];3],maxeval:usize)->(),|ctx,solver,lo,hi,direction,polarization,origin,program,rotation,maxeval|oblique(ctx,solver,lo,hi,direction,polarization,origin,program,rotation,maxeval,false));
object_binding!(set_oblique_plane_wave_rotated,"num_analysis_yee_set_oblique_plane_wave_rotated","crate::sim::yee_fdtd::Yee::set_oblique_plane_wave_rotated","Actual native Yee FDTD operation. Explicit positive SI length, wave-speed, impedance and electric-field references define c=epsilon0=mu0=1native coordinates. H uses E/Z, time L/c, energy E²L³/(Zc), conductivity1/(ZL). Grid≤4096and source/work bounds precede execution.",false,(solver:Handle<Yee>,lo:Input,hi:Input,direction:[f64;3],polarization:[f64;3],origin:Input,program:Input,rotation:[[f64;3];3],maxeval:usize)->(),|ctx,solver,lo,hi,direction,polarization,origin,program,rotation,maxeval|oblique(ctx,solver,lo,hi,direction,polarization,origin,program,rotation,maxeval,true));
object_binding!(state,"num_analysis_yee_state","","Native Yee state. The solver reserves48MiB: up to16MiBgrid/auxiliary/history storage and32MiBfully independent compiled snapshots, including at most4mutable sources. Snapshots remain valid after their original source handles are released.",false,(solver:Handle<Yee>)->Dynamic,|ctx,solver|solver.with(ctx,|s|Ok(Dynamic(api::record(vec![("time",api::quantity(api::scale_product(s.data.t,&[s.frame.time],&[])?,time())?),("dt",api::quantity(api::scale_product(s.data.dt,&[s.frame.time],&[])?,time())?),("sources",Output::encode(s.sources)?),("independent_snapshot_bytes",Output::encode(s.snapshot_bytes)?)])))));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &pole_lorentz::BINDING,
            &pole_drude::BINDING,
            &susceptibility::BINDING,
            &config_new::BINDING,
            &new::BINDING,
            &size::BINDING,
            &n::BINDING,
            &idx::BINDING,
            &h::BINDING,
            &energy::BINDING,
            &ohmic_loss::BINDING,
            &ade_dt_limit::BINDING,
            &field::BINDING,
            &field_mut::BINDING,
            &sigma::BINDING,
            &step::BINDING,
            &half_step_h::BINDING,
            &run::BINDING,
            &add_source::BINDING,
            &set_plane_wave::BINDING,
            &set_oblique_plane_wave::BINDING,
            &set_oblique_plane_wave_rotated::BINDING,
            &state::BINDING,
        ],
        "numerical_vector",
    )
}
