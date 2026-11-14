//! Native discontinuous Galerkin laws, meshes, boundaries and measured solutions.
use super::super::codecs::Arg;
use super::super::native_object_api::Maybe;
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::canonical;
use super::compiled_snapshot;
use super::err;
use super::mechanics::References;
use super::object_binding;
use super::payload;
use crate::sim::dg;
use crate::sim::dg::DgLaw;
use crate::sim::dg3d;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
#[derive(Clone)]
struct Law {
    data: Arc<dyn DgLaw>,
    frame: Frame,
    refs: Vec<(f64, DimExpr)>,
    dim: usize,
    gamma: Option<f64>,
}
impl NativeType for Law {
    const KIND: &'static str = "rssn/native/dg-law/v1";
    const MAX_BYTES: usize = 16384;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !(1..=3).contains(&self.dim)
            || self.refs.len() != self.data.nvar()
            || self.refs.len() > 5
            || self.refs.iter().any(|r| !r.0.is_finite() || r.0 <= 0.0)
        {
            return Err(err("DG law domain/ref invariant invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Mesh {
    data: dg::DgMesh,
    frame: Frame,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/dg-mesh/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !(1..=2).contains(&self.data.dim)
            || !(1..=8).contains(&self.data.p)
            || self.data.ne.contains(&0)
            || self.data.ne.iter().any(|n| *n > 256)
            || self.data.elements() * self.data.nodes_per_element() > 8192
            || (0..self.data.dim).any(|a| {
                !self.data.lo[a].is_finite()
                    || !self.data.hi[a].is_finite()
                    || self.data.hi[a] <= self.data.lo[a]
            })
        {
            return Err(err("DG mesh geometry/DOF budget invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Mesh3 {
    data: dg3d::Dg3Mesh,
    frame: Frame,
}
impl NativeType for Mesh3 {
    const KIND: &'static str = "rssn/native/dg3-mesh/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !(1..=6).contains(&self.data.p)
            || self.data.ne.contains(&0)
            || self.data.ne.iter().any(|n| *n > 128)
            || self.data.elements() * self.data.nodes_per_element() > 8192
            || (0..3).any(|a| {
                !self.data.lo[a].is_finite()
                    || !self.data.hi[a].is_finite()
                    || self.data.hi[a] <= self.data.lo[a]
            })
        {
            return Err(err("DG3 mesh geometry/DOF budget invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Bcs {
    data: [dg::DgBc; 4],
    law: Law,
    callbacks: Vec<CompiledSnapshot>,
}
impl NativeType for Bcs {
    const KIND: &'static str = "rssn/native/dg-bcs/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        self.law.validate()
    }
}
#[derive(Clone)]
struct Bcs3 {
    data: [dg3d::Dg3Bc; 6],
    law: Law,
    callbacks: Vec<CompiledSnapshot>,
}
impl NativeType for Bcs3 {
    const KIND: &'static str = "rssn/native/dg3-bcs/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        self.law.validate()
    }
}
#[derive(Clone)]
struct Solution {
    data: dg::DgSolution,
    law: Law,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/dg-solution/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.law.validate()?;
        if self.data.q.len() > 40960
            || self.data.q.iter().any(|v| !v.is_finite())
            || self.data.nvar != self.law.refs.len()
        {
            return Err(err(
                "DG solution finite component/storage invariant invalid",
            ));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Solution3 {
    data: dg3d::Dg3Solution,
    law: Law,
}
impl NativeType for Solution3 {
    const KIND: &'static str = "rssn/native/dg3-solution/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.law.validate()?;
        if self.data.q.len() > 40960
            || self.data.q.iter().any(|v| !v.is_finite())
            || self.data.nvar != self.law.refs.len()
        {
            return Err(err(
                "DG3 solution finite component/storage invariant invalid",
            ));
        }
        Ok(())
    }
}
fn geometry(
    f: &Frame,
    g: &Frame,
) -> Result<(), OpError> {
    for i in 0..2 {
        api::same_dimension(&f.dims[i], &g.dims[i])?;
    }
    if f.length.to_bits() != g.length.to_bits() || f.time.to_bits() != g.time.to_bits() {
        return Err(err("DG geometry/time references differ"));
    }
    Ok(())
}
fn states(
    l: &Law,
    q: &[f64],
    dir: usize,
) -> Result<(), OpError> {
    if dir >= l.dim || q.len() != l.refs.len() || q.iter().any(|v| !v.is_finite()) {
        return Err(err("DG state shape/direction invalid"));
    }
    if let Some((rho, p)) = l.data.rho_p(q)
        && (!rho.is_finite() || rho <= 0.0 || !p.is_finite() || p <= 0.0)
    {
        return Err(err("Euler density/pressure must be finite positive"));
    }
    Ok(())
}
fn measured(
    v: f64,
    reference: f64,
    dim: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::quantity(v * reference, dim)?))
}
fn power(
    d: &DimExpr,
    p: i64,
    q: i64,
) -> DimExpr {
    DimExpr::Pow(Box::new(d.clone()), Exp::Rat(Rational::new(p, q)))
}
fn result(
    l: &Law,
    v: Vec<f64>,
    flux: bool,
) -> Result<Dynamic, OpError> {
    let fields = v
        .into_iter()
        .enumerate()
        .map(|(i, v)| {
            let (scale, dim) = &l.refs[i];
            let (scale, dim) = if flux {
                (
                    api::scale_product(1.0, &[*scale, l.frame.length], &[l.frame.time])?,
                    canonical(DimExpr::Mul(vec![
                        dim.clone(),
                        l.frame.dims[0].clone(),
                        DimExpr::Inverse(Box::new(l.frame.dims[1].clone())),
                    ]))?,
                )
            } else {
                (*scale, dim.clone())
            };
            Ok((format!("component{i}"), api::quantity(v * scale, dim)?))
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Dynamic(Value::Record(crate::kernel::value::RecordValue {
        fields,
    })))
}
fn field(
    l: &Law,
    values: Vec<f64>,
    var: usize,
) -> Result<Dynamic, OpError> {
    if var >= l.refs.len() {
        return Err(err("DG variable out of range"));
    }
    Ok(Dynamic(api::record(vec![
        ("normalized", values.encode()?),
        (
            "reference",
            api::quantity(l.refs[var].0, l.refs[var].1.clone())?,
        ),
    ])))
}
fn programs(
    ctx: &CallCtx<'_>,
    v: &Input,
    l: &Law,
    spatial: usize,
    time: bool,
    state: bool,
    budget: usize,
) -> Result<Vec<CompiledSnapshot>, OpError> {
    let Value::Tensor(t) = &v.0 else {
        return Err(err("DG callback needs vector of scalar compiled programs"));
    };
    if t.shape != [l.refs.len()]
        || t.data.len() != l.refs.len()
        || budget == 0
        || budget > 8_388_608
    {
        return Err(err("DG callback component/budget mismatch"));
    }
    let mut dims = vec![l.frame.dims[0].clone(); spatial];
    if time {
        dims.insert(0, l.frame.dims[1].clone());
    }
    if state {
        dims.extend(l.refs.iter().map(|r| r.1.clone()));
    }
    t.data
        .iter()
        .zip(&l.refs)
        .map(|(v, r)| {
            super::structured_ops::checked_snapshot(ctx, &Input(v.clone()), budget, &dims, &r.1)
        })
        .collect()
}
fn eval(
    c: &[CompiledSnapshot],
    l: &Law,
    t: Option<f64>,
    x: &[f64],
    q: Option<&[f64]>,
) -> Vec<f64> {
    let mut input = vec![];
    if let Some(t) = t {
        input.push(t * l.frame.time);
    }
    input.extend(x.iter().map(|v| v * l.frame.length));
    if let Some(q) = q {
        input.extend(q.iter().zip(&l.refs).map(|(v, r)| v * r.0));
    }
    c.iter()
        .zip(&l.refs)
        .map(|(c, r)| c.sample(&input) / r.0)
        .collect()
}
fn finish(c: &[CompiledSnapshot]) -> Result<(), OpError> {
    for c in c {
        c.finish()?;
    }
    Ok(())
}
fn read<T: Arg>(
    v: &Value,
    name: &str,
) -> Result<T, OpError> {
    let v = api::field(v, name)?;
    api::normalized(v)?;
    T::decode(v, 0)
}
fn filter(v: &Value) -> Result<Option<dg::ModalFilter>, OpError> {
    let enabled: bool = read(v, "filter")?;
    let alpha: f64 = read(v, "filter_alpha")?;
    let order: u32 = read(v, "filter_order")?;
    if alpha < 0.0 || order == 0 || order > 32 || !order.is_multiple_of(2) {
        return Err(err(
            "modal filter needs nonnegative alpha and even order 2..32",
        ));
    }
    Ok(enabled.then_some(dg::ModalFilter { alpha, order }))
}
fn parse_options(v: &Value) -> Result<dg::DgOptions, OpError> {
    let m: f64 = read(v, "tvb_m")?;
    if m < 0.0 {
        return Err(err("TVB M must be nonnegative"));
    }
    let shock = match read::<String>(v, "shock")?.as_str() {
        | "None" => dg::ShockLimiter::None,
        | "Tvb" => dg::ShockLimiter::Tvb { m },
        | "TvbCharacteristic" => dg::ShockLimiter::TvbCharacteristic { m },
        | "Hierarchical" => dg::ShockLimiter::Hierarchical,
        | "Weno" => dg::ShockLimiter::Weno { m },
        | _ => return Err(err("unknown DG shock limiter")),
    };
    Ok(dg::DgOptions {
        filter: filter(v)?,
        shock,
        positivity: read(v, "positivity")?,
        entropy_stable: read(v, "entropy_stable")?,
    })
}
struct Meter<'a> {
    law: &'a dyn DgLaw,
    left: AtomicUsize,
    exhausted: std::sync::atomic::AtomicBool,
}
impl Meter<'_> {
    fn debit(&self) -> bool {
        if self
            .left
            .try_update(Ordering::Relaxed, Ordering::Relaxed, |v| v.checked_sub(1))
            .is_ok()
        {
            true
        } else {
            self.exhausted.store(true, Ordering::Relaxed);
            false
        }
    }

    fn finish(&self) -> Result<(), OpError> {
        if self.exhausted.load(Ordering::Relaxed) {
            Err(err("DG native law evaluation work budget exhausted"))
        } else {
            Ok(())
        }
    }
}
impl DgLaw for Meter<'_> {
    fn nvar(&self) -> usize {
        self.law.nvar()
    }

    fn flux(
        &self,
        q: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.flux(q, dir, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn max_speed(
        &self,
        q: &[f64],
        dir: usize,
    ) -> f64 {
        if self.debit() {
            self.law.max_speed(q, dir)
        } else {
            f64::INFINITY
        }
    }

    fn rho_p(
        &self,
        q: &[f64],
    ) -> Option<(f64, f64)> {
        self.law.rho_p(q)
    }

    fn has_characteristics(&self) -> bool {
        self.law.has_characteristics()
    }

    fn to_characteristic(
        &self,
        q: &[f64],
        d: usize,
        v: &[f64],
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.to_characteristic(q, d, v, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn from_characteristic(
        &self,
        q: &[f64],
        d: usize,
        v: &[f64],
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.from_characteristic(q, d, v, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn wall_state(
        &self,
        q: &[f64],
        d: usize,
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.wall_state(q, d, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn numerical_flux(
        &self,
        l: &[f64],
        r: &[f64],
        d: usize,
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.numerical_flux(l, r, d, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn two_point_flux(
        &self,
        l: &[f64],
        r: &[f64],
        d: usize,
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.two_point_flux(l, r, d, out);
        } else {
            out.fill(f64::NAN);
        }
    }

    fn entropy_stable_flux(
        &self,
        l: &[f64],
        r: &[f64],
        d: usize,
        out: &mut [f64],
    ) {
        if self.debit() {
            self.law.entropy_stable_flux(l, r, d, out);
        } else {
            out.fill(f64::NAN);
        }
    }
}
object_binding!(dg_law_advection,"num_dg_advection","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<Frame>,velocity:Vec<f64>)->NewObject<Law>,|ctx,frame,velocity|{frame.with(ctx,|f|{if !(1..=3).contains(&velocity.len()){return Err(err("advection dimension must be 1..3"));}let dim=velocity.len();let data:Arc<dyn DgLaw>=if dim==3{Arc::new(dg3d::LinearAdvection3d{a:velocity.try_into().map_err(|_|err("advection shape"))?})}else{Arc::new(dg::LinearAdvection{a:[velocity[0],*velocity.get(1).unwrap_or(&0.0)]})};Ok(NewObject::new(Law{data,frame:f.clone(),refs:vec![(f.value,f.dims[2].clone())],dim,gamma:None},f.dims.to_vec()))})});
object_binding!(dg_law_burgers,"num_dg_burgers","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<Frame>,dimension:usize)->NewObject<Law>,|ctx,frame,dimension|{frame.with(ctx,|f|{let f_reference_1_m1_0=f.reference([1, -1, 0])?;if !(1..=2).contains(&dimension){return Err(err("Burgers supports dimension 1/2"));}let velocity=canonical(DimExpr::Mul(vec![f.dims[0].clone(),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;api::same_dimension(&f.dims[2],&velocity)?;if f.value.to_bits()!=(f_reference_1_m1_0).to_bits(){return Err(err("Burgers field reference must equal L/T"));}Ok(NewObject::new(Law{data:Arc::new(dg::Burgers),frame:f.clone(),refs:vec![(f.value,velocity)],dim:dimension,gamma:None},f.dims.to_vec()))})});
object_binding!(dg_law_euler,"num_dg_euler","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<References>,dimension:usize,gamma:f64)->NewObject<Law>,|ctx,frame,dimension,gamma|{frame.with(ctx,|r|{if !(1..=3).contains(&dimension)||gamma<=1.0||gamma>3.0{return Err(err("Euler dimension/gamma invalid"));}let mut refs=vec![r.derived([-3,0,1,0])?];for _ in 0..dimension{refs.push(r.derived([-2,-1,1,0])?);}refs.push(r.derived([-1,-2,1,0])?);let f=Frame{length:r.scales[0],time:r.scales[1],value:refs[0].0,dims:[r.dims[0].clone(),r.dims[1].clone(),refs[0].1.clone()]};let data:Arc<dyn DgLaw>=match dimension{1=>Arc::new(dg::Euler1d{gamma}),2=>Arc::new(dg::Euler2d{gamma}),_=>Arc::new(dg3d::Euler3d{gamma})};let dims=refs.iter().map(|r|r.1.clone()).collect();Ok(NewObject::new(Law{data,frame:f,refs,dim:dimension,gamma:Some(gamma)},dims))})});
object_binding!(dg_law_nvar,"num_dg_nvar","crate::sim::dg::DgLaw::nvar","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>)->usize,|ctx,law|{law.with(ctx,|l|Ok(l.data.nvar()))});
object_binding!(dg_law_has_characteristics,"num_dg_has_characteristics","crate::sim::dg::DgLaw::has_characteristics","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>)->bool,|ctx,law|{law.with(ctx,|l|Ok(l.data.has_characteristics()))});
object_binding!(dg_law_flux,"num_dg_flux","crate::sim::dg::DgLaw::flux","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize)->Dynamic,|ctx,law,state,direction|{law.with(ctx,|l|{states(l,&state,direction)?;let mut out=vec![0.0;l.refs.len()];l.data.flux(&state,direction,&mut out);result(l,out,true)})});
object_binding!(dg_law_numerical_flux,"num_dg_numerical_flux","crate::sim::dg::DgLaw::numerical_flux","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize,right:Vec<f64>)->Dynamic,|ctx,law,state,direction,right|{law.with(ctx,|l|{states(l,&state,direction)?;states(l,&right,direction)?;let mut out=vec![0.0;l.refs.len()];l.data.numerical_flux(&state,&right,direction,&mut out);result(l,out,true)})});
object_binding!(dg_law_two_point_flux,"num_dg_two_point_flux","crate::sim::dg::DgLaw::two_point_flux","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize,right:Vec<f64>)->Dynamic,|ctx,law,state,direction,right|{law.with(ctx,|l|{states(l,&state,direction)?;states(l,&right,direction)?;let mut out=vec![0.0;l.refs.len()];l.data.two_point_flux(&state,&right,direction,&mut out);result(l,out,true)})});
object_binding!(dg_law_entropy_stable_flux,"num_dg_entropy_stable_flux","crate::sim::dg::DgLaw::entropy_stable_flux","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize,right:Vec<f64>)->Dynamic,|ctx,law,state,direction,right|{law.with(ctx,|l|{states(l,&state,direction)?;states(l,&right,direction)?;let mut out=vec![0.0;l.refs.len()];l.data.entropy_stable_flux(&state,&right,direction,&mut out);result(l,out,true)})});
object_binding!(dg_law_wall_state,"num_dg_wall_state","crate::sim::dg::DgLaw::wall_state","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize)->Dynamic,|ctx,law,state,direction|{law.with(ctx,|l|{states(l,&state,direction)?;let mut out=vec![0.0;l.refs.len()];l.data.wall_state(&state,direction,&mut out);result(l,out,false)})});
object_binding!(dg_law_max_speed,"num_dg_max_speed","crate::sim::dg::DgLaw::max_speed","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize)->Dynamic,|ctx,law,state,direction|{law.with(ctx,|l|{states(l,&state,direction)?;measured(l.data.max_speed(&state,direction),l.frame.length/l.frame.time,canonical(DimExpr::Mul(vec![l.frame.dims[0].clone(),DimExpr::Inverse(Box::new(l.frame.dims[1].clone()))]))?)})});
object_binding!(dg_law_rho_p,"num_dg_rho_p","crate::sim::dg::DgLaw::rho_p","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>)->Maybe<Dynamic>,|ctx,law,state|{law.with(ctx,|l|{states(l,&state,0)?;Ok(Maybe(l.data.rho_p(&state).map(|(rho,p)|Ok::<Dynamic,OpError>(Dynamic(api::record(vec![("density",api::quantity(rho*l.refs[0].0,l.refs[0].1.clone())?),("pressure",api::quantity(p*l.refs.last().ok_or_else(||err("missing energy ref"))?.0,l.refs.last().ok_or_else(||err("missing energy ref"))?.1.clone())?)])))).transpose()?))})});
object_binding!(dg_law_to_characteristic,"num_dg_to_characteristic","crate::sim::dg::DgLaw::to_characteristic","Native characteristic transformation in the explicitly normalized conserved coordinate basis. Eigenvectors mix physical component domains, so its coordinates are dimensionless relative to the law references; inverse uses the identical basis.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize,vector:Vec<f64>)->Vec<f64>,|ctx,law,state,direction,vector|{law.with(ctx,|l|{states(l,&state,direction)?;if vector.len()!=l.refs.len(){return Err(err("characteristic vector size mismatch"));}let mut out=vec![0.0;vector.len()];l.data.to_characteristic(&state,direction,&vector,&mut out);Ok(out)})});
object_binding!(dg_law_from_characteristic,"num_dg_from_characteristic","crate::sim::dg::DgLaw::from_characteristic","Native characteristic transformation in the explicitly normalized conserved coordinate basis. Eigenvectors mix physical component domains, so its coordinates are dimensionless relative to the law references; inverse uses the identical basis.",true,(law:Handle<Law>,state:Vec<f64>,direction:usize,vector:Vec<f64>)->Vec<f64>,|ctx,law,state,direction,vector|{law.with(ctx,|l|{states(l,&state,direction)?;if vector.len()!=l.refs.len(){return Err(err("characteristic vector size mismatch"));}let mut out=vec![0.0;vector.len()];l.data.from_characteristic(&state,direction,&vector,&mut out);Ok(out)})});
object_binding!(dg_mesh_line,"num_dg_mesh_line","crate::sim::dg::DgMesh::line","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<Frame>,elements:usize,lower:f64,upper:f64,degree:usize)->NewObject<Mesh>,|ctx,frame,elements,lower,upper,degree|{frame.with(ctx,|f|{if elements==0||elements>256||degree==0||degree>8{return Err(err("DG mesh count/degree budget invalid"));}Ok(NewObject::new(Mesh{data:dg::DgMesh::line(elements,lower,upper,degree),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(dg_mesh_rect,"num_dg_mesh_rect","crate::sim::dg::DgMesh::rect","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<Frame>,elements:Vec<usize>,lower:Vec<f64>,upper:Vec<f64>,degree:usize)->NewObject<Mesh>,|ctx,frame,elements,lower,upper,degree|{frame.with(ctx,|f|{let ne:[usize;2]=elements.try_into().map_err(|_|err("DG rectangle needs two extents"))?;if ne.contains(&0)||ne.iter().any(|n|*n>256)||degree==0||degree>8{return Err(err("DG mesh count/degree budget invalid"));}Ok(NewObject::new(Mesh{data:dg::DgMesh::rect(ne,lower.try_into().map_err(|_|err("DG lower shape"))?,upper.try_into().map_err(|_|err("DG upper shape"))?,degree),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(dg3_mesh,"num_dg3_mesh","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(frame:Handle<Frame>,elements:Vec<usize>,lower:Vec<f64>,upper:Vec<f64>,degree:usize)->NewObject<Mesh3>,|ctx,frame,elements,lower,upper,degree|{frame.with(ctx,|f|{let ne:[usize;3]=elements.try_into().map_err(|_|err("DG box needs three extents"))?;if ne.contains(&0)||ne.iter().any(|n|*n>128)||degree==0||degree>6{return Err(err("DG3 mesh count/degree budget invalid"));}Ok(NewObject::new(Mesh3{data:dg3d::Dg3Mesh{ne,lo:lower.try_into().map_err(|_|err("DG3 lower shape"))?,hi:upper.try_into().map_err(|_|err("DG3 upper shape"))?,p:degree},frame:f.clone()},f.dims.to_vec()))})});
object_binding!(dg_mesh_elements,"num_dg_elements","crate::sim::dg::DgMesh::elements","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh>)->usize,|ctx,mesh|{mesh.with(ctx,|m|Ok(m.data.elements()))});
object_binding!(dg_mesh_nodes_per_element,"num_dg_nodes_per_element","crate::sim::dg::DgMesh::nodes_per_element","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh>)->usize,|ctx,mesh|{mesh.with(ctx,|m|Ok(m.data.nodes_per_element()))});
object_binding!(dg3_mesh_elements,"num_dg3_elements","crate::sim::dg3d::Dg3Mesh::elements","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh3>)->usize,|ctx,mesh|{mesh.with(ctx,|m|Ok(m.data.elements()))});
object_binding!(dg3_mesh_nodes_per_element,"num_dg3_nodes_per_element","crate::sim::dg3d::Dg3Mesh::nodes_per_element","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh3>)->usize,|ctx,mesh|{mesh.with(ctx,|m|Ok(m.data.nodes_per_element()))});
object_binding!(dg_lgl,"num_dg_lgl","crate::sim::dg::lgl","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(degree:usize)->Dynamic,|ctx,degree|{let _=ctx;if !(1..=64).contains(&degree){return Err(err("LGL degree must be 1..64"));}let(x,w)=dg::lgl(degree);Ok(Dynamic(api::record(vec![("nodes",x.encode()?),("weights",w.encode()?)])))});
object_binding!(dg_lgl_diff,"num_dg_lgl_diff","crate::sim::dg::lgl_diff_matrix","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(nodes:Vec<f64>)->Vec<Vec<f64>>,|ctx,nodes|{let _=ctx;if nodes.len()<2||nodes.len()>65||nodes.iter().enumerate().any(|(i,x)|nodes.iter().skip(i+1).any(|y|crate::kernels::exact_eq(*x,*y))){return Err(err("differentiation nodes require 2..65 distinct normalized values"));}Ok(dg::lgl_diff_matrix(&nodes))});
object_binding!(dg_bcs,"num_dg_bcs","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,kinds:Vec<String>)->NewObject<Bcs>,|ctx,law,kinds|{law.with(ctx,|l|{if l.dim>2||kinds.len()!=4{return Err(err("DG BC dimension/face count mismatch"));}let mut bcs:[dg::DgBc;4]=std::array::from_fn(|_|dg::DgBc::Outflow);for(i,k)in kinds.iter().enumerate(){bcs[i]=match k.as_str(){"Periodic"=>dg::DgBc::Periodic,"Outflow"=>dg::DgBc::Outflow,"Wall"=>dg::DgBc::Wall,_=>return Err(err("use Periodic,Outflow,Wall then install a Fixed/Custom callback face"))};}Ok(NewObject::new(Bcs{data:bcs,law:l.clone(),callbacks:vec![]},l.refs.iter().map(|r|r.1.clone()).collect()))})});
object_binding!(dg_bc_fixed,"num_dg_bc_fixed","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(bcs:Handle<Bcs>,face:usize,callbacks:Input,budget:usize)->NewObject<Bcs>,|ctx,bcs,face,callbacks,budget|{bcs.with(ctx,|b|{if face>=4{return Err(err("DG face index outside array"));}let c=programs(ctx,&callbacks,&b.law,2,true,false,budget)?;let mut out=b.clone();let l=b.law.clone();let cc=c.clone();out.data[face]=dg::DgBc::Fixed(Arc::new(move|t,x|eval(&cc,&l,Some(t),&x,None)));out.callbacks.extend(c);if out.callbacks.len()>64{return Err(err("DG retained callback count budget exceeded"));}let mut parents=vec![bcs.value.clone()];if let Value::Tensor(t)=&callbacks.0{parents.extend(t.data.clone());}Ok(NewObject::dependent(out,b.law.refs.iter().map(|r|r.1.clone()).collect(),parents))})});
object_binding!(dg_bc_custom,"num_dg_bc_custom","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(bcs:Handle<Bcs>,face:usize,callbacks:Input,budget:usize)->NewObject<Bcs>,|ctx,bcs,face,callbacks,budget|{bcs.with(ctx,|b|{if face>=4{return Err(err("DG face index outside array"));}let c=programs(ctx,&callbacks,&b.law,2,true,true,budget)?;let mut out=b.clone();let l=b.law.clone();let cc=c.clone();out.data[face]=dg::DgBc::Custom(Arc::new(move|t,x,q|eval(&cc,&l,Some(t),&x,Some(q))));out.callbacks.extend(c);if out.callbacks.len()>64{return Err(err("DG retained callback count budget exceeded"));}let mut parents=vec![bcs.value.clone()];if let Value::Tensor(t)=&callbacks.0{parents.extend(t.data.clone());}Ok(NewObject::dependent(out,b.law.refs.iter().map(|r|r.1.clone()).collect(),parents))})});
object_binding!(dg3_bcs,"num_dg3_bcs","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,kinds:Vec<String>)->NewObject<Bcs3>,|ctx,law,kinds|{law.with(ctx,|l|{if l.dim!=3||kinds.len()!=6{return Err(err("DG BC dimension/face count mismatch"));}let mut bcs:[dg3d::Dg3Bc;6]=std::array::from_fn(|_|dg3d::Dg3Bc::Outflow);for(i,k)in kinds.iter().enumerate(){bcs[i]=match k.as_str(){"Periodic"=>dg3d::Dg3Bc::Periodic,"Outflow"=>dg3d::Dg3Bc::Outflow,"Wall"=>dg3d::Dg3Bc::Wall,_=>return Err(err("use Periodic,Outflow,Wall then install a Fixed/Custom callback face"))};}Ok(NewObject::new(Bcs3{data:bcs,law:l.clone(),callbacks:vec![]},l.refs.iter().map(|r|r.1.clone()).collect()))})});
object_binding!(dg3_bc_fixed,"num_dg3_bc_fixed","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(bcs:Handle<Bcs3>,face:usize,callbacks:Input,budget:usize)->NewObject<Bcs3>,|ctx,bcs,face,callbacks,budget|{bcs.with(ctx,|b|{if face>=6{return Err(err("DG face index outside array"));}let c=programs(ctx,&callbacks,&b.law,3,true,false,budget)?;let mut out=b.clone();let l=b.law.clone();let cc=c.clone();out.data[face]=dg3d::Dg3Bc::Fixed(Arc::new(move|t,x|eval(&cc,&l,Some(t),&x,None)));out.callbacks.extend(c);if out.callbacks.len()>64{return Err(err("DG retained callback count budget exceeded"));}let mut parents=vec![bcs.value.clone()];if let Value::Tensor(t)=&callbacks.0{parents.extend(t.data.clone());}Ok(NewObject::dependent(out,b.law.refs.iter().map(|r|r.1.clone()).collect(),parents))})});
object_binding!(dg_solve_dg,"num_solve_dg","crate::sim::dg::solve_dg","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(mesh:Handle<Mesh>,bcs:Handle<Bcs>,initial:Input,t_end:f64,cfl:f64,options:Input,callback_budget:usize,work_budget:usize)->NewObject<Solution>,|ctx,mesh,bcs,initial,t_end,cfl,options,callback_budget,work_budget|{mesh.with(ctx,|m|bcs.with(ctx,|b|{geometry(&m.frame,&b.law.frame)?;if b.law.dim!=m.data.dim||t_end<0.0||cfl<=0.0||cfl>1.0||work_budget==0||work_budget>8_388_608{return Err(err("DG dimension/time/CFL/work budget invalid"));}let init=programs(ctx,&initial,&b.law,2,false,false,callback_budget)?;let meter=Meter{law:b.law.data.as_ref(),left:AtomicUsize::new(work_budget),exhausted:std::sync::atomic::AtomicBool::new(false)};let out=dg::solve_dg(&m.data,&meter,&b.data,&|x|eval(&init,&b.law,None,&x,None),t_end,cfl,filter(&options.0)?).map_err(err);finish(&init)?;finish(&b.callbacks)?;meter.finish()?;let data=out?;Ok(NewObject::new(Solution{data,law:b.law.clone()},b.law.refs.iter().map(|r|r.1.clone()).collect()))}))});
object_binding!(dg_solve_dg_limited,"num_solve_dg_limited","crate::sim::dg::solve_dg_limited","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(mesh:Handle<Mesh>,bcs:Handle<Bcs>,initial:Input,t_end:f64,cfl:f64,options:Input,callback_budget:usize,work_budget:usize)->NewObject<Solution>,|ctx,mesh,bcs,initial,t_end,cfl,options,callback_budget,work_budget|{mesh.with(ctx,|m|bcs.with(ctx,|b|{geometry(&m.frame,&b.law.frame)?;if b.law.dim!=m.data.dim||t_end<0.0||cfl<=0.0||cfl>1.0||work_budget==0||work_budget>8_388_608{return Err(err("DG dimension/time/CFL/work budget invalid"));}let init=programs(ctx,&initial,&b.law,2,false,false,callback_budget)?;let meter=Meter{law:b.law.data.as_ref(),left:AtomicUsize::new(work_budget),exhausted:std::sync::atomic::AtomicBool::new(false)};let tvb_m:f64=read(&options.0,"tvb_m")?;if tvb_m<0.0{return Err(err("TVB M must be nonnegative"));}let limiter=read::<bool>(&options.0,"limiter")?.then_some(dg::DgLimiter{tvb_m:read::<bool>(&options.0,"tvb")?.then_some(tvb_m),positivity:read(&options.0,"positivity")?});let out=dg::solve_dg_limited(&m.data,&meter,&b.data,&|x|eval(&init,&b.law,None,&x,None),t_end,cfl,filter(&options.0)?,limiter).map_err(err);finish(&init)?;finish(&b.callbacks)?;meter.finish()?;let data=out?;Ok(NewObject::new(Solution{data,law:b.law.clone()},b.law.refs.iter().map(|r|r.1.clone()).collect()))}))});
object_binding!(dg_solve_dg_with,"num_solve_dg_with","crate::sim::dg::solve_dg_with","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(mesh:Handle<Mesh>,bcs:Handle<Bcs>,initial:Input,t_end:f64,cfl:f64,options:Input,callback_budget:usize,work_budget:usize)->NewObject<Solution>,|ctx,mesh,bcs,initial,t_end,cfl,options,callback_budget,work_budget|{mesh.with(ctx,|m|bcs.with(ctx,|b|{geometry(&m.frame,&b.law.frame)?;if b.law.dim!=m.data.dim||t_end<0.0||cfl<=0.0||cfl>1.0||work_budget==0||work_budget>8_388_608{return Err(err("DG dimension/time/CFL/work budget invalid"));}let init=programs(ctx,&initial,&b.law,2,false,false,callback_budget)?;let meter=Meter{law:b.law.data.as_ref(),left:AtomicUsize::new(work_budget),exhausted:std::sync::atomic::AtomicBool::new(false)};let opt=parse_options(&options.0)?;let out=dg::solve_dg_with(&m.data,&meter,&b.data,&|x|eval(&init,&b.law,None,&x,None),t_end,cfl,&opt).map_err(err);finish(&init)?;finish(&b.callbacks)?;meter.finish()?;let data=out?;Ok(NewObject::new(Solution{data,law:b.law.clone()},b.law.refs.iter().map(|r|r.1.clone()).collect()))}))});
object_binding!(dg_solve_dg3d,"num_solve_dg3d","crate::sim::dg3d::solve_dg3d","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(mesh:Handle<Mesh3>,bcs:Handle<Bcs3>,initial:Input,t_end:f64,cfl:f64,options:Input,callback_budget:usize,work_budget:usize)->NewObject<Solution3>,|ctx,mesh,bcs,initial,t_end,cfl,options,callback_budget,work_budget|{mesh.with(ctx,|m|bcs.with(ctx,|b|{geometry(&m.frame,&b.law.frame)?;if b.law.dim!=3||t_end<0.0||cfl<=0.0||cfl>1.0||work_budget==0||work_budget>8_388_608{return Err(err("DG dimension/time/CFL/work budget invalid"));}let init=programs(ctx,&initial,&b.law,3,false,false,callback_budget)?;let meter=Meter{law:b.law.data.as_ref(),left:AtomicUsize::new(work_budget),exhausted:std::sync::atomic::AtomicBool::new(false)};let limiter=match read::<String>(&options.0,"limiter")?.as_str(){"None"=>dg3d::Dg3Limiter::None,"Minmax"=>dg3d::Dg3Limiter::Minmax,"Characteristic"=>dg3d::Dg3Limiter::Characteristic,"Hierarchical"=>dg3d::Dg3Limiter::Hierarchical,"Weno"=>dg3d::Dg3Limiter::Weno,_=>return Err(err("unknown DG3 limiter"))};let tvd:f64=read(&options.0,"tvb_m")?;let eps:f64=read(&options.0,"weno_epsilon")?;if tvd<0.0||eps<=0.0{return Err(err("DG3 limiter parameters invalid"));}let opt=dg3d::Dg3Options{entropy_stable:read(&options.0,"entropy_stable")?,positivity:read(&options.0,"positivity")?,limiter,limiter_tvd:tvd,limiter_weno_eps:eps};let out=dg3d::solve_dg3d(&m.data,&meter,&b.data,&|x|eval(&init,&b.law,None,&x,None),t_end,cfl,&opt).map_err(err);finish(&init)?;finish(&b.callbacks)?;meter.finish()?;let data=out?;Ok(NewObject::new(Solution3{data,law:b.law.clone()},b.law.refs.iter().map(|r|r.1.clone()).collect()))}))});
object_binding!(dg_advection_step,"num_dg_advection_step","crate::sim::dg::advection_step","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh>,velocity:f64,state:Vec<f64>,time:f64,dt:f64)->Dynamic,|ctx,mesh,velocity,state,time,dt|{mesh.with(ctx,|m|{let out=dg::advection_step(&m.data,velocity,&state,time,dt).map_err(err)?;payload(&m.frame,out,vec![state.len()])})});
object_binding!(dg_advection_step_nd,"num_dg_advection_step_nd","crate::sim::dg::advection_step_nd","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(mesh:Handle<Mesh>,velocity:Vec<f64>,state:Vec<f64>,time:f64,dt:f64)->Dynamic,|ctx,mesh,velocity,state,time,dt|{mesh.with(ctx,|m|{let out=dg::advection_step_nd(&m.data,velocity.try_into().map_err(|_|err("two advection velocity components required"))?,&state,time,dt).map_err(err)?;payload(&m.frame,out,vec![state.len()])})});
object_binding!(dg_value,"num_dg_value","crate::sim::dg::DgSolution::value","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution>,element:usize,node:usize,variable:usize)->Dynamic,|ctx,solution,element,node,variable|{solution.with(ctx,|s|{if element>=s.data.mesh.elements()||node>=s.data.mesh.nodes_per_element()||variable>=s.law.refs.len(){return Err(err("DG solution index outside shape"));}measured(s.data.value(element,node,variable),s.law.refs[variable].0,s.law.refs[variable].1.clone())})});
object_binding!(dg_coord,"num_dg_coord","crate::sim::dg::DgSolution::coord","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution>,element:usize,node:usize)->Dynamic,|ctx,solution,element,node|{solution.with(ctx,|s|{if element>=s.data.mesh.elements()||node>=s.data.mesh.nodes_per_element(){return Err(err("DG solution node index outside shape"));}let x=s.data.coord(element,node);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![x.len()],elem:TypeExpr::Dim(Box::new(TypeExpr::Float),s.law.frame.dims[0].clone()),data:x.into_iter().map(|v|api::quantity(v*s.law.frame.length,s.law.frame.dims[0].clone())).collect::<Result<_,_>>()?})))})});
object_binding!(dg_l2_error,"num_dg_l2_error","crate::sim::dg::DgSolution::l2_error","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(solution:Handle<Solution>,variable:usize,exact:Input,budget:usize)->Dynamic,|ctx,solution,variable,exact,budget|{solution.with(ctx,|s|{if variable>=s.law.refs.len(){return Err(err("DG component out of range"));}let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![s.law.frame.dims[0].clone();2],&s.law.refs[variable].1)?;let value=s.data.l2_error(variable,|x|c.sample(&x.map(|v|v*s.law.frame.length))/s.law.refs[variable].0);c.finish()?;let exponent=s.data.mesh.dim as i64;measured(value,s.law.refs[variable].0*s.law.frame.length.powf(exponent as f64/2.0),canonical(DimExpr::Mul(vec![s.law.refs[variable].1.clone(),power(&s.law.frame.dims[0],exponent,2)]))?)})});
object_binding!(dg_integral,"num_dg_integral","crate::sim::dg::DgSolution::integral","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution>,variable:usize)->Dynamic,|ctx,solution,variable|{solution.with(ctx,|s|{if variable>=s.law.refs.len(){return Err(err("DG component out of range"));}let exponent=s.data.mesh.dim as i64;measured(s.data.integral(variable),s.law.refs[variable].0*s.law.frame.length.powi(exponent as i32),canonical(DimExpr::Mul(vec![s.law.refs[variable].1.clone(),power(&s.law.frame.dims[0],exponent,1)]))?)})});
object_binding!(dg_data,"num_dg_data","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution>,variable:usize)->Dynamic,|ctx,solution,variable|{solution.with(ctx,|s|{if variable>=s.law.refs.len(){return Err(err("DG component out of range"));}field(&s.law,s.data.q.chunks(s.law.refs.len()).map(|v|v[variable]).collect(),variable)})});
object_binding!(dg3_value,"num_dg3_value","crate::sim::dg3d::Dg3Solution::value","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution3>,element:usize,node:usize,variable:usize)->Dynamic,|ctx,solution,element,node,variable|{solution.with(ctx,|s|{if element>=s.data.mesh.elements()||node>=s.data.mesh.nodes_per_element()||variable>=s.law.refs.len(){return Err(err("DG solution index outside shape"));}measured(s.data.value(element,node,variable),s.law.refs[variable].0,s.law.refs[variable].1.clone())})});
object_binding!(dg3_coord,"num_dg3_coord","crate::sim::dg3d::Dg3Solution::coord","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution3>,element:usize,node:usize)->Dynamic,|ctx,solution,element,node|{solution.with(ctx,|s|{if element>=s.data.mesh.elements()||node>=s.data.mesh.nodes_per_element(){return Err(err("DG solution node index outside shape"));}let x=s.data.coord(element,node);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![x.len()],elem:TypeExpr::Dim(Box::new(TypeExpr::Float),s.law.frame.dims[0].clone()),data:x.into_iter().map(|v|api::quantity(v*s.law.frame.length,s.law.frame.dims[0].clone())).collect::<Result<_,_>>()?})))})});
object_binding!(dg3_l2_error,"num_dg3_l2_error","crate::sim::dg3d::Dg3Solution::l2_error","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(solution:Handle<Solution3>,variable:usize,exact:Input,budget:usize)->Dynamic,|ctx,solution,variable,exact,budget|{solution.with(ctx,|s|{if variable>=s.law.refs.len(){return Err(err("DG component out of range"));}let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![s.law.frame.dims[0].clone();3],&s.law.refs[variable].1)?;let value=s.data.l2_error(variable,|x|c.sample(&x.map(|v|v*s.law.frame.length))/s.law.refs[variable].0);c.finish()?;let exponent=3;measured(value,s.law.refs[variable].0*s.law.frame.length.powf(exponent as f64/2.0),canonical(DimExpr::Mul(vec![s.law.refs[variable].1.clone(),power(&s.law.frame.dims[0],exponent,2)]))?)})});
object_binding!(dg3_integral,"num_dg3_integral","crate::sim::dg3d::Dg3Solution::integral","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",false,(solution:Handle<Solution3>,function:Input,budget:usize)->Dynamic,|ctx,solution,function,budget|{solution.with(ctx,|s|{let c=compiled_snapshot(ctx,&function.0,budget)?;if c.input_dimensions.len()!=s.law.refs.len(){return Err(err("DG3 state functional arity mismatch"));}for(a,b)in c.input_dimensions.iter().zip(&s.law.refs){api::same_dimension(a,&b.1)?;}let value=s.data.integral(|q|c.sample(&q.iter().zip(&s.law.refs).map(|(v,r)|v*r.0).collect::<Vec<_>>()));c.finish()?;measured(value,s.law.frame.length.powi(3),canonical(DimExpr::Mul(vec![c.output_dimension,power(&s.law.frame.dims[0],3,1)]))?)})});
object_binding!(dg3_data,"num_dg3_data","","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(solution:Handle<Solution3>,variable:usize)->Dynamic,|ctx,solution,variable|{solution.with(ctx,|s|{if variable>=s.law.refs.len(){return Err(err("DG component out of range"));}field(&s.law,s.data.q.chunks(s.law.refs.len()).map(|v|v[variable]).collect(),variable)})});
object_binding!(dg3_pressure,"num_dg3_pressure","crate::sim::dg3d::Euler3d::pressure","Actual native DG method with explicit geometry/time and per-component physical references; normalized states, separate flux units and coherent compiled callbacks preserve mixed Euler dimensions. Native algorithm and limiter identity are retained; shapes, work and parent lifetimes are checked.",true,(law:Handle<Law>,state:Vec<f64>)->Dynamic,|ctx,law,state|{law.with(ctx,|l|{states(l,&state,0)?;if l.dim!=3{return Err(err("Euler3d pressure needs a three-dimensional Euler law"));}let gamma=l.gamma.ok_or_else(||err("pressure needs an Euler law"))?;measured(dg3d::Euler3d{gamma}.pressure(&state),l.refs[4].0,l.refs[4].1.clone())})});
object_binding!(dg_entropy,"num_dg_entropy","crate::sim::dg3d::euler_entropy","Native mathematical Euler entropy in dimensionless reference coordinates. The logarithm uses normalized pressure and density; changing physical reference scales intentionally changes this entropy coordinate.",true,(gamma:f64,velocity_dimensions:usize,state:Vec<f64>)->f64,|ctx,gamma,velocity_dimensions,state|{let _=ctx;if !(1..=3).contains(&velocity_dimensions)||state.len()!=velocity_dimensions+2||gamma<=1.0||gamma>3.0||state[0]<=0.0{return Err(err("normalized Euler entropy domain invalid"));}let p=(gamma-1.0)*(state[velocity_dimensions+1]-0.5*state[1..=velocity_dimensions].iter().map(|v|v*v).sum::<f64>()/state[0]);if p<=0.0{return Err(err("normalized Euler entropy requires positive pressure"));}Ok(dg3d::euler_entropy(gamma,velocity_dimensions,&state))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &dg_law_advection::BINDING,
    &dg_law_burgers::BINDING,
    &dg_law_euler::BINDING,
    &dg_law_nvar::BINDING,
    &dg_law_has_characteristics::BINDING,
    &dg_law_flux::BINDING,
    &dg_law_numerical_flux::BINDING,
    &dg_law_two_point_flux::BINDING,
    &dg_law_entropy_stable_flux::BINDING,
    &dg_law_wall_state::BINDING,
    &dg_law_max_speed::BINDING,
    &dg_law_rho_p::BINDING,
    &dg_law_to_characteristic::BINDING,
    &dg_law_from_characteristic::BINDING,
    &dg_mesh_line::BINDING,
    &dg_mesh_rect::BINDING,
    &dg3_mesh::BINDING,
    &dg_mesh_elements::BINDING,
    &dg_mesh_nodes_per_element::BINDING,
    &dg3_mesh_elements::BINDING,
    &dg3_mesh_nodes_per_element::BINDING,
    &dg_lgl::BINDING,
    &dg_lgl_diff::BINDING,
    &dg_bcs::BINDING,
    &dg_bc_fixed::BINDING,
    &dg_bc_custom::BINDING,
    &dg3_bcs::BINDING,
    &dg3_bc_fixed::BINDING,
    &dg_solve_dg::BINDING,
    &dg_solve_dg_limited::BINDING,
    &dg_solve_dg_with::BINDING,
    &dg_solve_dg3d::BINDING,
    &dg_advection_step::BINDING,
    &dg_advection_step_nd::BINDING,
    &dg_value::BINDING,
    &dg_coord::BINDING,
    &dg_l2_error::BINDING,
    &dg_integral::BINDING,
    &dg_data::BINDING,
    &dg3_value::BINDING,
    &dg3_coord::BINDING,
    &dg3_l2_error::BINDING,
    &dg3_integral::BINDING,
    &dg3_data::BINDING,
    &dg3_pressure::BINDING,
    &dg_entropy::BINDING,
];
