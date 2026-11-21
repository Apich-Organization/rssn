//! Native compressible/shallow-water finite volumes and adaptive mesh refinement.
use super::super::codecs::Arg;
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
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::canonical;
use super::err;
use super::mechanics::References;
use super::object_binding;
use crate::sim::euler_amr as amr;
use crate::sim::euler_nd as fv;
#[derive(Clone)]
struct Grid {
    data: fv::FvGrid,
    frame: References,
}
impl NativeType for Grid {
    const KIND: &'static str = "rssn/native/fv-grid/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !(1..=3).contains(&self.data.dim)
            || self.data.n.iter().any(|n| *n == 0 || *n > 512)
            || self.data.len() > 4096
            || (0..self.data.dim).any(|a| {
                !self.data.lo[a].is_finite()
                    || !self.data.hi[a].is_finite()
                    || self.data.hi[a] <= self.data.lo[a]
            })
        {
            return Err(err("finite-volume geometry/shape/work budget invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Config {
    data: fv::FvConfig,
    frame: References,
}
impl NativeType for Config {
    const KIND: &'static str = "rssn/native/fv-config/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !self.data.cfl.is_finite() || self.data.cfl <= 0.0 || self.data.cfl >= 1.0 {
            return Err(err("FV CFL must be in (0,1)"));
        }
        match self.data.law {
            | fv::Law::Euler { gamma } if !gamma.is_finite() || gamma <= 1.0 || gamma > 3.0 => {
                Err(err("gamma must be in (1,3]"))
            },
            | fv::Law::ShallowWater { g } if !g.is_finite() || g <= 0.0 => {
                Err(err("shallow gravity must be finite positive"))
            },
            | _ => Ok(()),
        }
    }
}
#[derive(Clone)]
struct Amr {
    data: amr::AmrConfig,
    frame: References,
}
impl NativeType for Amr {
    const KIND: &'static str = "rssn/native/fv-amr-config/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        Config {
            data: self.data.fv,
            frame: self.frame.clone(),
        }
        .validate()?;
        if self.data.max_level > 3
            || self.data.regrid_interval == 0
            || self.data.regrid_interval > 4096
            || self.data.buffer > 4
            || !self.data.threshold.is_finite()
            || self.data.threshold < 0.0
            || !self.data.efficiency.is_finite()
            || self.data.efficiency <= 0.0
            || self.data.efficiency > 1.0
        {
            return Err(err("AMR level/regrid/tagging/efficiency bounds invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Solution {
    data: fv::FvSolution,
    cfg: Config,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/fv-solution/v1";
    const MAX_BYTES: usize = 512 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.cfg.validate()?;
        if self.data.q.len() != self.data.grid.len()
            || self.data.q.len() > 4096
            || self.data.q.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(err("FV solution finite shape invariant invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct AmrSolution {
    data: amr::AmrSolution,
    cfg: Config,
}
impl NativeType for AmrSolution {
    const KIND: &'static str = "rssn/native/fv-amr-solution/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.cfg.validate()?;
        if self.data.levels.len() > 4
            || self
                .data
                .levels
                .iter()
                .flat_map(|l| &l.patches)
                .try_fold(0usize, |n, p| n.checked_add(p.2.len()))
                .is_none_or(|n| n > 16384)
            || self
                .data
                .levels
                .iter()
                .flat_map(|l| &l.patches)
                .flat_map(|p| &p.2)
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(err("AMR retained hierarchy finite storage budget invalid"));
        }
        Ok(())
    }
}
fn refs(
    c: &Config,
    primitive: bool,
) -> Result<Vec<(f64, DimExpr)>, OpError> {
    match c.data.law {
        | fv::Law::Euler { .. } => {
            let r = c.frame.derived([-3, 0, 1, 0])?;
            let v = c.frame.derived(if primitive {
                [1, -1, 0, 0]
            } else {
                [-2, -1, 1, 0]
            })?;
            let p = c.frame.derived([-1, -2, 1, 0])?;
            Ok(vec![r, v.clone(), v.clone(), v, p])
        },
        | fv::Law::ShallowWater { .. } => {
            let r = c.frame.derived([1, 0, 0, 0])?;
            let v = c.frame.derived(if primitive {
                [1, -1, 0, 0]
            } else {
                [2, -1, 0, 0]
            })?;
            Ok(vec![r, v.clone(), v.clone(), v, (1.0, DimExpr::Dimless)])
        },
    }
}
fn same(
    a: &References,
    b: &References,
) -> Result<(), OpError> {
    for (i, x) in a.scales.iter().enumerate() {
        if x.to_bits() != b.scales[i].to_bits() {
            return Err(err("FV mechanical reference frames differ"));
        }
        api::same_dimension(&a.dims[i], &b.dims[i])?;
    }
    Ok(())
}
fn state(v: Vec<f64>) -> Result<fv::State, OpError> {
    v.try_into()
        .map_err(|_| err("native primitive/conserved state has five slots"))
}
fn bc(
    kind: &str,
    v: Vec<f64>,
) -> Result<fv::GasBc, OpError> {
    Ok(match kind {
        | "Transmissive" if v.is_empty() => fv::GasBc::Transmissive,
        | "Reflective" if v.is_empty() => fv::GasBc::Reflective,
        | "Periodic" if v.is_empty() => fv::GasBc::Periodic,
        | "Fixed" => fv::GasBc::Fixed(state(v)?),
        | "NonReflecting" => fv::GasBc::NonReflecting(state(v)?),
        | _ => return Err(err("unknown gas BC or unused state controls")),
    })
}
fn read<T: Arg>(
    v: &Value,
    name: &str,
) -> Result<T, OpError> {
    let v = api::field(v, name)?;
    api::normalized(v)?;
    T::decode(v, 0)
}
fn scheme(
    c: &mut fv::FvConfig,
    options: &Value,
) -> Result<(), OpError> {
    c.cfl = read(options, "cfl")?;
    c.recon = match read::<String>(options, "reconstruction")?.as_str() {
        | "FirstOrder" => fv::Reconstruction::FirstOrder,
        | "Weno5" => fv::Reconstruction::Weno5,
        | "Minmod" => fv::Reconstruction::Muscl(fv::SlopeLimiter::Minmod),
        | "VanLeer" => fv::Reconstruction::Muscl(fv::SlopeLimiter::VanLeer),
        | "Mc" => fv::Reconstruction::Muscl(fv::SlopeLimiter::Mc),
        | "Superbee" => fv::Reconstruction::Muscl(fv::SlopeLimiter::Superbee),
        | _ => return Err(err("unknown FV reconstruction")),
    };
    c.riemann = match read::<String>(options, "riemann")?.as_str() {
        | "Hllc" => fv::RiemannSolver::Hllc,
        | "Hlle" => fv::RiemannSolver::Hlle,
        | "Roe" => fv::RiemannSolver::Roe,
        | "Rusanov" => fv::RiemannSolver::Rusanov,
        | _ => return Err(err("unknown Riemann solver")),
    };
    Ok(())
}
fn programs(
    ctx: &CallCtx<'_>,
    c: &Config,
    v: &Input,
    budget: usize,
) -> Result<Vec<CompiledSnapshot>, OpError> {
    let Value::Tensor(t) = &v.0 else {
        return Err(err(
            "initial primitive state requires five compiled scalar programs",
        ));
    };
    if t.shape != [5] || t.data.len() != 5 || budget == 0 || budget > 8_388_608 {
        return Err(err("FV initial callback shape/work budget invalid"));
    }
    let r = refs(c, true)?;
    t.data
        .iter()
        .zip(r)
        .map(|(p, r)| {
            super::structured_ops::checked_snapshot(
                ctx,
                &Input(p.clone()),
                budget,
                &vec![c.frame.dims[0].clone(); 3],
                &r.1,
            )
        })
        .collect()
}
fn eval(
    c: &Config,
    p: &[CompiledSnapshot],
    x: [f64; 3],
    r: &[(f64, DimExpr)],
) -> fv::State {
    std::array::from_fn(|i| p[i].sample(&x.map(|v| v * c.frame.scales[0])) / r[i].0)
}
fn finish(p: &[CompiledSnapshot]) -> Result<(), OpError> {
    for p in p {
        p.finish()?;
    }
    Ok(())
}
fn checked_steps(
    grid: &fv::FvGrid,
    steps: usize,
    levels: usize,
) -> Result<(), OpError> {
    if steps == 0
        || steps > 4096
        || grid
            .len()
            .checked_mul(steps)
            .and_then(|n| n.checked_mul(1usize << levels))
            .is_none_or(|n| n > 4_194_304)
    {
        return Err(err("FV/AMR step/work budget exceeded"));
    }
    Ok(())
}
fn measured(
    v: f64,
    c: &Config,
    component: usize,
    primitive: bool,
    length_power: usize,
) -> Result<Dynamic, OpError> {
    let r = refs(c, primitive)?;
    let r = r
        .get(component)
        .ok_or_else(|| err("state component out of range"))?;
    let d = canonical(DimExpr::Mul(vec![
        r.1.clone(),
        DimExpr::Pow(
            Box::new(c.frame.dims[0].clone()),
            Exp::Rat(Rational::new(length_power as i64, 1)),
        ),
    ]))?;
    Ok(Dynamic(api::quantity(
        v * r.0 * c.frame.scales[0].powi(length_power as i32),
        d,
    )?))
}
fn result(
    values: fv::State,
    c: &Config,
    primitive: bool,
    dim: usize,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(Value::Record(crate::kernel::value::RecordValue {
        fields: values
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                Ok((
                    format!("component{i}"),
                    measured(v, c, i, primitive, dim)?.0,
                ))
            })
            .collect::<Result<_, OpError>>()?,
    })))
}
object_binding!(fv_grid_new,"num_fv_grid_new","crate::sim::euler_nd::FvGrid::new","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(frame:Handle<References>,extents:Vec<usize>,lower:Vec<f64>,upper:Vec<f64>)->NewObject<Grid>,|ctx,frame,extents,lower,upper|{frame.with(ctx,|f|{if extents.len()>3||extents.iter().any(|n|*n>512)||extents.iter().try_fold(1usize,|a,b|a.checked_mul(*b)).is_none_or(|n|n>4096){return Err(err("FV grid allocation budget exceeded"));}Ok(NewObject::new(Grid{data:fv::FvGrid::new(&extents,&lower,&upper).map_err(err)?,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fv_grid_len,"num_fv_grid_len","crate::sim::euler_nd::FvGrid::len","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>)->usize,|ctx,grid|{grid.with(ctx,|g|Ok(g.data.len()))});
object_binding!(fv_grid_is_empty,"num_fv_grid_is_empty","crate::sim::euler_nd::FvGrid::is_empty","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>)->bool,|ctx,grid|{grid.with(ctx,|g|Ok(g.data.is_empty()))});
object_binding!(fv_grid_h,"num_fv_grid_h","crate::sim::euler_nd::FvGrid::h","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>,axis:usize)->Dynamic,|ctx,grid,axis|{grid.with(ctx,|g|{if axis>=g.data.dim{return Err(err("inactive FV axis"));}Ok(Dynamic(api::quantity(g.data.h(axis)*g.frame.scales[0],g.frame.dims[0].clone())?))})});
object_binding!(fv_grid_volume,"num_fv_grid_cell_volume","crate::sim::euler_nd::FvGrid::cell_volume","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>)->Dynamic,|ctx,grid|{grid.with(ctx,|g|{let(scale,d)=g.frame.derived([g.data.dim as i64,0,0,0])?;Ok(Dynamic(api::quantity(g.data.cell_volume()*scale,d)?))})});
object_binding!(fv_grid_index,"num_fv_grid_index","crate::sim::euler_nd::FvGrid::index","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>,index:Vec<usize>)->usize,|ctx,grid,index|{grid.with(ctx,|g|{let i:[usize;3]=index.try_into().map_err(|_|err("FV index needs three slots"))?;if i.iter().zip(g.data.n).any(|(i,n)|*i>=n){return Err(err("FV index outside grid"));}Ok(g.data.index(i))})});
object_binding!(fv_grid_ijk,"num_fv_grid_ijk","crate::sim::euler_nd::FvGrid::ijk","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>,flat:usize)->Vec<usize>,|ctx,grid,flat|{grid.with(ctx,|g|{if flat>=g.data.len(){return Err(err("FV flat index outside grid"));}Ok(g.data.ijk(flat).to_vec())})});
object_binding!(fv_grid_centre,"num_fv_grid_centre","crate::sim::euler_nd::FvGrid::centre","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(grid:Handle<Grid>,flat:usize)->Dynamic,|ctx,grid,flat|{grid.with(ctx,|g|{if flat>=g.data.len(){return Err(err("FV flat index outside grid"));}let x=g.data.centre(flat);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],elem:TypeExpr::Dim(Box::new(TypeExpr::Float),g.frame.dims[0].clone()),data:x.into_iter().map(|v|api::quantity(v*g.frame.scales[0],g.frame.dims[0].clone())).collect::<Result<_,_>>()?})))})});
object_binding!(fv_config_euler,"num_fv_config_euler","crate::sim::euler_nd::FvConfig::euler","Creates a ideal-gas Euler configuration with the dimensionless adiabatic exponent. The options record contains cfl, reconstruction and riemann. Riemann solver names are Hllc, Hlle, Roe and Rusanov. Hlle uses two Einfeldt wave speeds.",true,(frame:Handle<References>,gamma:f64,options:Input)->NewObject<Config>,|ctx,frame,gamma,options|{frame.with(ctx,|f|{let mut data=fv::FvConfig::euler(gamma);scheme(&mut data,&options.0)?;Ok(NewObject::new(Config{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fv_config_shallow_water,"num_fv_config_shallow_water","crate::sim::euler_nd::FvConfig::shallow_water","Creates a shallow-water configuration with gravity in the frame’s acceleration reference. The options record contains cfl, reconstruction and riemann. Riemann solver names are Hllc, Hlle, Roe and Rusanov. Hlle uses two Einfeldt wave speeds.",true,(frame:Handle<References>,gravity:f64,options:Input)->NewObject<Config>,|ctx,frame,gravity,options|{frame.with(ctx,|f|{let mut data=fv::FvConfig::shallow_water(gravity);scheme(&mut data,&options.0)?;Ok(NewObject::new(Config{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fv_config_with_face,"num_fv_with_face","crate::sim::euler_nd::FvConfig::with_face","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(config:Handle<Config>,face:usize,kind:String,primitive:Vec<f64>)->NewObject<Config>,|ctx,config,face,kind,primitive|{config.with(ctx,|c|{if face>=6{return Err(err("FV BC index outside range"));}let data=c.data.with_face(face,bc(&kind,primitive)?);Ok(NewObject::new(Config{data,frame:c.frame.clone()},c.frame.dims.to_vec()))})});
object_binding!(fv_config_with_axis,"num_fv_with_axis","crate::sim::euler_nd::FvConfig::with_axis","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(config:Handle<Config>,axis:usize,kind:String,primitive:Vec<f64>)->NewObject<Config>,|ctx,config,axis,kind,primitive|{config.with(ctx,|c|{if axis>=3{return Err(err("FV BC index outside range"));}let data=c.data.with_axis(axis,bc(&kind,primitive)?);Ok(NewObject::new(Config{data,frame:c.frame.clone()},c.frame.dims.to_vec()))})});
object_binding!(fv_solve,"num_fv_solve","crate::sim::euler_nd::solve","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",false,(grid:Handle<Grid>,config:Handle<Config>,initial:Input,t_end:f64,max_steps:usize,callback_budget:usize)->NewObject<Solution>,|ctx,grid,config,initial,t_end,max_steps,callback_budget|{grid.with(ctx,|g|config.with(ctx,|c|{same(&g.frame,&c.frame)?;checked_steps(&g.data,max_steps,0)?;if t_end<0.0{return Err(err("FV final time must be nonnegative"));}let p=programs(ctx,c,&initial,callback_budget)?;let r=refs(c,true)?;let out=crate::sim::models_geodesic_relativity::with_native_budget(max_steps,||fv::solve(&g.data,&c.data,&|x|eval(c,&p,x,&r),t_end)).map_err(err)?;finish(&p)?;let data=out.map_err(err)?;Ok(NewObject::new(Solution{data,cfg:c.clone()},c.frame.dims.to_vec()))}))});
object_binding!(fv_solve_with_sources,"num_fv_solve_with_sources","crate::sim::euler_nd::solve_with_sources","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",false,(grid:Handle<Grid>,config:Handle<Config>,initial:Input,t_end:f64,max_steps:usize,callback_budget:usize,gravity:Vec<f64>,bed:Input,manning:f64)->NewObject<Solution>,|ctx,grid,config,initial,t_end,max_steps,callback_budget,gravity,bed,manning|{grid.with(ctx,|g|config.with(ctx,|c|{same(&g.frame,&c.frame)?;checked_steps(&g.data,max_steps,0)?;if t_end<0.0{return Err(err("FV final time must be nonnegative"));}let p=programs(ctx,c,&initial,callback_budget)?;let r=refs(c,true)?;let gravity:[f64;3]=gravity.try_into().map_err(|_|err("gravity requires three normalized acceleration components"))?;if manning<0.0{return Err(err("normalized Manning coefficient must be nonnegative"));}let Value::Variant(b)=&bed.0 else{return Err(err("bed callback needs exact Option envelope"));};let b=match(b.tag.as_str(),b.payload.as_slice()){("None",[])=>None,("Some",[v])=>Some(crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&Input(v.clone()),callback_budget,&vec![c.frame.dims[0].clone();3],&c.frame.dims[0])?),_=>return Err(err("malformed bed callback option"))};let bedfn=|x:[f64;3]|b.as_ref().map_or(0.0,|b|b.sample(&x.map(|v|v*c.frame.scales[0]))/c.frame.scales[0]);let source=fv::Sources{gravity,bed:b.as_ref().map(|_|&bedfn as &dyn Fn([f64;3])->f64),manning};let out=crate::sim::models_geodesic_relativity::with_native_budget(max_steps,||fv::solve_with_sources(&g.data,&c.data,&|x|eval(c,&p,x,&r),&source,t_end)).map_err(err)?;finish(&p)?;if let Some(b)=&b{b.finish()?;}let data=out.map_err(err)?;Ok(NewObject::new(Solution{data,cfg:c.clone()},c.frame.dims.to_vec()))}))});
object_binding!(fv_primitive,"num_fv_primitive","crate::sim::euler_nd::FvSolution::primitive","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(solution:Handle<Solution>,cell:usize)->Dynamic,|ctx,solution,cell|{solution.with(ctx,|s|{if cell>=s.data.grid.len(){return Err(err("FV cell outside grid"));}result(s.data.primitive(cell),&s.cfg,true,0)})});
object_binding!(fv_totals,"num_fv_totals","crate::sim::euler_nd::FvSolution::totals","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(solution:Handle<Solution>)->Dynamic,|ctx,solution|{solution.with(ctx,|s|result(s.data.totals(),&s.cfg,false,s.data.grid.dim))});
object_binding!(fv_error,"num_fv_l1_error","crate::sim::euler_nd::FvSolution::l1_error","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",false,(solution:Handle<Solution>,component:usize,exact:Input,budget:usize)->Dynamic,|ctx,solution,component,exact,budget|{solution.with(ctx,|s|{let refs=refs(&s.cfg,true)?;let r=refs.get(component).ok_or_else(||err("FV component out of range"))?;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![s.cfg.frame.dims[0].clone();3],&r.1)?;let v=s.data.l1_error(component,|x|c.sample(&x.map(|v|v*s.cfg.frame.scales[0]))/r.0);c.finish()?;measured(v,&s.cfg,component,true,s.data.grid.dim)})});
object_binding!(amr_config,"num_fv_amr_config","crate::sim::euler_amr::AmrConfig::new","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(config:Handle<Config>,max_level:usize,regrid_interval:usize,threshold:f64,buffer:usize,efficiency:f64,reflux:bool)->NewObject<Amr>,|ctx,config,max_level,regrid_interval,threshold,buffer,efficiency,reflux|{config.with(ctx,|c|{let mut data=amr::AmrConfig::new(c.data);data.max_level=max_level;data.regrid_interval=regrid_interval;data.threshold=threshold;data.buffer=buffer;data.efficiency=efficiency;data.reflux=reflux;Ok(NewObject::new(Amr{data,frame:c.frame.clone()},c.frame.dims.to_vec()))})});
object_binding!(amr_solve,"num_fv_solve_amr","crate::sim::euler_amr::solve_amr","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",false,(grid:Handle<Grid>,config:Handle<Amr>,initial:Input,t_end:f64,max_steps:usize,callback_budget:usize)->NewObject<AmrSolution>,|ctx,grid,config,initial,t_end,max_steps,callback_budget|{grid.with(ctx,|g|config.with(ctx,|a|{same(&g.frame,&a.frame)?;let factor=1usize<<(a.data.max_level*g.data.dim);if g.data.len().checked_mul(factor).is_none_or(|n|n>4096){return Err(err("AMR full tagging hierarchy exceeds allocation budget"));}checked_steps(&g.data,max_steps,a.data.max_level*g.data.dim)?;if t_end<0.0{return Err(err("AMR final time must be nonnegative"));}let c=Config{data:a.data.fv,frame:a.frame.clone()};let p=programs(ctx,&c,&initial,callback_budget)?;let r=refs(&c,true)?;let out=crate::sim::models_geodesic_relativity::with_native_budget(max_steps,||amr::solve_amr(&g.data,&a.data,&|x|eval(&c,&p,x,&r),t_end)).map_err(err)?;finish(&p)?;Ok(NewObject::new(AmrSolution{data:out.map_err(err)?,cfg:c},a.frame.dims.to_vec()))}))});
object_binding!(amr_cells,"num_fv_amr_cells_per_level","crate::sim::euler_amr::AmrSolution::cells_per_level","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(solution:Handle<AmrSolution>)->Vec<usize>,|ctx,solution|{solution.with(ctx,|s|Ok(s.data.cells_per_level()))});
object_binding!(amr_totals,"num_fv_amr_totals","crate::sim::euler_amr::AmrSolution::totals","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(solution:Handle<AmrSolution>)->Dynamic,|ctx,solution|{solution.with(ctx,|s|result(s.data.totals(),&s.cfg,false,s.data.grid.dim))});
object_binding!(amr_primitive,"num_fv_amr_primitive_on_level","crate::sim::euler_amr::AmrSolution::primitive_on_level","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",true,(solution:Handle<AmrSolution>,level:usize,component:usize)->Dynamic,|ctx,solution,level,component|{solution.with(ctx,|s|{if level>3||component>=5||s.data.grid.len().checked_mul(1usize<<(level*s.data.grid.dim)).is_none_or(|n|n>4096){return Err(err("AMR uniform projection allocation/component budget exceeded"));}let values=s.data.primitive_on_level(level);let r=refs(&s.cfg,true)?;Ok(Dynamic(api::record(vec![("normalized",values.iter().map(|v|v[component]).collect::<Vec<_>>().encode()?),("reference",api::quantity(r[component].0,r[component].1.clone())?)])))})});
object_binding!(amr_error,"num_fv_amr_composite_l1_error","crate::sim::euler_amr::AmrSolution::composite_l1_error","Actual native finite-volume/AMR method in explicit Length/Time/Mass/Temperature references. Primitive slots are density or depth, three velocities, pressure or unused zero; conserved and integrated outputs use their distinct physical units. Native reconstruction, Riemann solver, sources, boundary identity and scoped step/work limits are retained.",false,(solution:Handle<AmrSolution>,component:usize,exact:Input,budget:usize)->Dynamic,|ctx,solution,component,exact,budget|{solution.with(ctx,|s|{let r=refs(&s.cfg,true)?;let r=r.get(component).ok_or_else(||err("AMR component out of range"))?;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![s.cfg.frame.dims[0].clone();3],&r.1)?;let v=s.data.composite_l1_error(component,|x|c.sample(&x.map(|v|v*s.cfg.frame.scales[0]))/r.0);c.finish()?;measured(v,&s.cfg,component,true,s.data.grid.dim)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &fv_grid_new::BINDING,
    &fv_grid_len::BINDING,
    &fv_grid_is_empty::BINDING,
    &fv_grid_h::BINDING,
    &fv_grid_volume::BINDING,
    &fv_grid_index::BINDING,
    &fv_grid_ijk::BINDING,
    &fv_grid_centre::BINDING,
    &fv_config_euler::BINDING,
    &fv_config_shallow_water::BINDING,
    &fv_config_with_face::BINDING,
    &fv_config_with_axis::BINDING,
    &fv_solve::BINDING,
    &fv_solve_with_sources::BINDING,
    &fv_primitive::BINDING,
    &fv_totals::BINDING,
    &fv_error::BINDING,
    &amr_config::BINDING,
    &amr_solve::BINDING,
    &amr_cells::BINDING,
    &amr_totals::BINDING,
    &amr_primitive::BINDING,
    &amr_error::BINDING,
];
