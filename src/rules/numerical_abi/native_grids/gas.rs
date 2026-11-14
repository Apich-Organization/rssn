//! Mixed-dimension one-dimensional compressible gas calls.
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::OpError;
use super::Output;
use super::RuleError;
use super::api;
use super::err;
use super::mechanics::References;
use super::object_binding;
use super::structured_ops::checked_snapshot;
use crate::sim::gas_dynamics as gas;
fn state(v: &[f64]) -> Result<gas::Primitive, OpError> {
    if v.len() != 3 || v[0] <= 0.0 || v[2] <= 0.0 {
        return Err(err(
            "normalized primitive must be [positive density, velocity, positive pressure]",
        ));
    }
    Ok(gas::Primitive {
        rho: v[0],
        u: v[1],
        p: v[2],
    })
}
fn gamma(g: f64) -> Result<(), OpError> {
    if g <= 1.0 || g > 10.0 {
        return Err(err("gas gamma must be in (1,10]"));
    }
    Ok(())
}
fn limiter(s: &str) -> Result<gas::Limiter, OpError> {
    match s {
        | "None" => Ok(gas::Limiter::None),
        | "Minmod" => Ok(gas::Limiter::Minmod),
        | "VanLeer" => Ok(gas::Limiter::VanLeer),
        | _ => Err(err("limiter must be None, Minmod or VanLeer")),
    }
}
fn output(
    f: &References,
    x: Vec<f64>,
    w: Vec<gas::Primitive>,
) -> Result<Dynamic, OpError> {
    let refs = [
        f.derived([-3, 0, 1, 0])?,
        f.derived([1, -1, 0, 0])?,
        f.derived([-1, -2, 1, 0])?,
    ];
    Ok(Dynamic(api::record(vec![
        ("coordinates", x.encode()?),
        (
            "states",
            w.into_iter()
                .map(|w| vec![w.rho, w.u, w.p])
                .collect::<Vec<_>>()
                .encode()?,
        ),
        (
            "coordinate_reference",
            api::quantity(f.scales[0], f.dims[0].clone())?,
        ),
        (
            "density_reference",
            api::quantity(refs[0].0, refs[0].1.clone())?,
        ),
        (
            "velocity_reference",
            api::quantity(refs[1].0, refs[1].1.clone())?,
        ),
        (
            "pressure_reference",
            api::quantity(refs[2].0, refs[2].1.clone())?,
        ),
    ])))
}
object_binding!(gas_hllc,"num_grid_gas_hllc","crate::sim::gas_dynamics::hllc_flux","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64|{gamma(adiabatic)?;frame.with(ctx,|f|{let flux=gas::hllc_flux(state(&left)?,state(&right)?,adiabatic);let powers=[[-2,-1,1,0],[-1,-2,1,0],[0,-3,1,0]];let mut fields=Vec::new();for (i,name) in ["mass_flux","momentum_flux","energy_flux"].into_iter().enumerate(){let(s,d)=f.derived(powers[i])?;fields.push((name,api::quantity(flux[i]*s,d)?));}Ok(Dynamic(api::record(fields)))})});
object_binding!(gas_exact,"num_grid_gas_exact","crate::sim::gas_dynamics::exact_riemann","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64,similarity:f64)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64,similarity:f64|{gamma(adiabatic)?;frame.with(ctx,|f|{let w=gas::exact_riemann(state(&left)?,state(&right)?,adiabatic,similarity);output(f,vec![],vec![w])})});
object_binding!(gas_sod,"num_grid_gas_sod","crate::sim::gas_dynamics::sod_shock_tube","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,cells:usize,time:f64,reconstruction:String,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,cells:usize,time:f64,reconstruction:String,budget:usize|{if !(6..=1024).contains(&cells)||time<0.0||budget==0||budget>100000{return Err(err("gas grid/time/step budget invalid"));}let lim=limiter(&reconstruction)?;frame.with(ctx,|f|{let(x,rho,p)=crate::sim::models_geodesic_relativity::with_native_budget(budget,||gas::sod_shock_tube(cells,time,lim)).map_err(err)?;Ok(Dynamic(api::record(vec![("coordinates",x.encode()?),("density",rho.encode()?),("pressure",p.encode()?),("coordinate_reference",api::quantity(f.scales[0],f.dims[0].clone())?),("density_reference",{let(s,d)=f.derived([-3,0,1,0])?;api::quantity(s,d)?}),("pressure_reference",{let(s,d)=f.derived([-1,-2,1,0])?;api::quantity(s,d)?})])))})});
object_binding!(gas_solve,"num_grid_gas_solve","crate::sim::gas_dynamics::solve_euler_1d","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,density:Input,velocity:Input,pressure:Input,lo:f64,hi:f64,cells:usize,time:f64,adiabatic:f64,cfl:f64,reconstruction:String,boundary:String,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,density:Input,velocity:Input,pressure:Input,lo:f64,hi:f64,cells:usize,time:f64,adiabatic:f64,cfl:f64,reconstruction:String,boundary:String,budget:usize|{gamma(adiabatic)?;if !(6..=1024).contains(&cells)||hi<=lo||time<0.0||cfl<=0.0||cfl>=1.0||budget==0||budget>100000{return Err(err("gas domain/time/grid/CFL/budget invalid"));}let lim=limiter(&reconstruction)?;let bc=match boundary.as_str(){"Reflective"=>gas::Boundary::Reflective,"Transmissive"=>gas::Boundary::Transmissive,_=>return Err(err("gas boundary must be Reflective or Transmissive"))};frame.with(ctx,|f|{let(rs,rd)=f.derived([-3,0,1,0])?;let(vs,vd)=f.derived([1,-1,0,0])?;let(ps,pd)=f.derived([-1,-2,1,0])?;let inputs=[f.dims[0].clone()];let rho=checked_snapshot(ctx,&density,cells*8,&inputs,&rd)?;let vel=checked_snapshot(ctx,&velocity,cells*8,&inputs,&vd)?;let pre=checked_snapshot(ctx,&pressure,cells*8,&inputs,&pd)?;let invalid=std::sync::atomic::AtomicBool::new(false);let result=crate::sim::models_geodesic_relativity::with_native_budget(budget,||gas::solve_euler_1d(|x|{let w=gas::Primitive{rho:rho.sample(&[x*f.scales[0]])/rs,u:vel.sample(&[x*f.scales[0]])/vs,p:pre.sample(&[x*f.scales[0]])/ps};
if !w.rho.is_finite()||w.rho<=0.0||!w.u.is_finite()||!w.p.is_finite()||w.p<=0.0{invalid.store(true,std::sync::atomic::Ordering::Relaxed);}w},lo,hi,cells,time,adiabatic,cfl,lim,bc));rho.finish()?;vel.finish()?;pre.finish()?;if invalid.load(std::sync::atomic::Ordering::Relaxed){return Err(err("gas initial callback returned invalid primitive"));}let(x,w)=result.map_err(err)?;output(f,x,w)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &gas_hllc::BINDING,
            &gas_exact::BINDING,
            &gas_sod::BINDING,
            &gas_solve::BINDING,
        ],
        "numerical_vector",
    )
}
