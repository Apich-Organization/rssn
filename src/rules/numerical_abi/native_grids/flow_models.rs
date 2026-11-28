//! Actual legacy square-channel/cavity model calls in mechanical references.
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Installer;
use super::OpError;
use super::Output;
use super::RuleError;
use super::Value;
use super::api;
use super::err;
use super::mechanics::References;
use super::object_binding;
use crate::sim::models_navier_stokes_fluid as flow;
use ndarray::Array2;
fn validate(
    n: usize,
    re: f64,
    dt: f64,
    steps: usize,
) -> Result<(), OpError> {
    if !(3..=65).contains(&n)
        || re <= 0.0
        || dt <= 0.0
        || steps > 4096
        || n.checked_mul(n)
            .and_then(|n| n.checked_mul(steps))
            .is_none_or(|n| n > 8_388_608)
    {
        return Err(err("square flow shape/coefficient/work budget invalid"));
    }
    Ok(())
}
fn result(
    out: flow::NavierStokesOutput,
    r: &References,
) -> Result<Dynamic, OpError> {
    let (u, v, p) = out.map_err(err)?;
    let velocity = r.derived([1, -1, 0, 0])?;
    let pressure = r.derived([-1, -2, 1, 0])?;
    let encode = |a: Array2<f64>, reference: &(f64, DimExpr)| -> Result<Value, OpError> {
        if a.len() > 65536 || a.iter().any(|v| !v.is_finite()) {
            return Err(err(
                "flow output is not finite or exceeds allocation bounds",
            ));
        }
        Ok(api::record(vec![
            ("shape", a.shape().to_vec().encode()?),
            (
                "normalized",
                a.iter().copied().collect::<Vec<_>>().encode()?,
            ),
            (
                "reference",
                api::quantity(reference.0, reference.1.clone())?,
            ),
        ]))
    };
    Ok(Dynamic(api::record(vec![
        ("u", encode(u, &velocity)?),
        ("v", encode(v, &velocity)?),
        ("pressure", encode(p, &pressure)?),
    ])))
}
object_binding!(flow_profile,"num_flow_inflow_velocity","crate::sim::models_navier_stokes_fluid::InflowProfile::velocity","Actual native unit-height inflow profile evaluated at normalized y; amplitude is normalized by the mechanical L/T reference and result retains velocity units.",true,(frame:Handle<References>,kind:String,amplitude:f64,height:f64)->Dynamic,|ctx,frame,kind,amplitude,height|{frame.with(ctx,|f|{if !(0.0..=1.0).contains(&height){return Err(err("inflow height must be in [0,1]"));}let profile=match kind.as_str(){"Uniform"=>flow::InflowProfile::Uniform(amplitude),"Parabolic"=>flow::InflowProfile::Parabolic(amplitude),_=>return Err(err("inflow kind must be Uniform or Parabolic"))};let(scale,d)=f.derived([1,-1,0,0])?;Ok(Dynamic(api::quantity(profile.velocity(height)*scale,d)?))})});
object_binding!(flow_channel,"num_flow_run_channel","crate::sim::models_navier_stokes_fluid::run_channel_flow","Native square unit-domain channel with original free-slip walls and uniform unit inflow; normalized dt uses the explicit Time reference, output velocities and pressure have mechanical units. The native rectangular rejection is preserved.",true,(frame:Handle<References>,nx:usize,ny:usize,reynolds:f64,dt:f64,steps:usize,obstacle:Vec<bool>)->Dynamic,|ctx,frame,nx,ny,reynolds,dt,steps,obstacle|{frame.with(ctx,|f|{validate(nx,reynolds,dt,steps)?;if nx!=ny||obstacle.len()!=nx*ny{return Err(err("channel needs a square grid and matching obstacle mask"));}let mask=Array2::from_shape_vec((ny,nx),obstacle).map_err(err)?;result(flow::run_channel_flow(nx,ny,reynolds,dt,steps,&mask),f)})});
object_binding!(flow_channel_with,"num_flow_run_channel_with","crate::sim::models_navier_stokes_fluid::run_channel_flow_with","Actual native configurable channel model on a square unit domain; inflow, no-slip/slip wall, obstacle mask, Reynolds number and normalized time step are passed unchanged. Physical outputs retain independent velocity/pressure references.",true,(frame:Handle<References>,size:usize,reynolds:f64,dt:f64,steps:usize,obstacle:Vec<bool>,profile:String,amplitude:f64,walls:String)->Dynamic,|ctx,frame,size,reynolds,dt,steps,obstacle,profile,amplitude,walls|{frame.with(ctx,|f|{validate(size,reynolds,dt,steps)?;if obstacle.len()!=size*size{return Err(err("channel obstacle mask shape mismatch"));}let mask=Array2::from_shape_vec((size,size),obstacle).map_err(err)?;let inflow=match profile.as_str(){"Uniform"=>flow::InflowProfile::Uniform(amplitude),"Parabolic"=>flow::InflowProfile::Parabolic(amplitude),_=>return Err(err("unknown inflow profile"))};let walls=match walls.as_str(){"Slip"=>flow::WallCondition::Slip,"NoSlip"=>flow::WallCondition::NoSlip,_=>return Err(err("unknown channel wall condition"))};result(flow::run_channel_flow_with(&flow::ChannelFlowConfig{n:size,re:reynolds,dt,n_iter:steps,obstacle_mask:&mask,inflow,walls}),f)})});
object_binding!(flow_cavity,"num_flow_run_lid_cavity","crate::sim::models_navier_stokes_fluid::run_lid_driven_cavity","Actual native legacy staggered cavity on a square unit domain. Lid velocity and time are normalized by explicit mechanical references. Native pressure/multigrid behavior is retained and work is bounded before invocation.",true,(frame:Handle<References>,nx:usize,ny:usize,reynolds:f64,dt:f64,steps:usize,lid_velocity:f64)->Dynamic,|ctx,frame,nx,ny,reynolds,dt,steps,lid_velocity|{frame.with(ctx,|f|{validate(nx,reynolds,dt,steps)?;if nx!=ny{return Err(err("cavity requires a square grid"));}result(flow::run_lid_driven_cavity(&flow::NavierStokesParameters{nx,ny,re:reynolds,dt,n_iter:steps,lid_velocity}),f)})});
object_binding!(flow_scenario,"num_flow_save_lid_cavity_scenario","crate::sim::models_navier_stokes_fluid::simulate_lid_driven_cavity_scenario","Runs the actual native fixed 65x65, 200-step cavity scenario and writes its three native NPY files to the explicitly provided output directory. This legacy routine reports failures itself; the result only acknowledges return from the routine, not successful file creation.",true,(output_directory:String,work_budget:usize)->bool,|ctx,output_directory,work_budget|{let _=ctx;if output_directory.is_empty()||output_directory.len()>4096||!(65*65*200..=8_388_608).contains(&work_budget){return Err(err("scenario output path/work budget invalid"));}flow::simulate_lid_driven_cavity_scenario(std::path::Path::new(&output_directory));Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &flow_profile::BINDING,
            &flow_channel::BINDING,
            &flow_channel_with::BINDING,
            &flow_cavity::BINDING,
            &flow_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
