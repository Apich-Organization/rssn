#![allow(clippy::wildcard_imports)] // Geometrized model parameters use the same measured relativity frame as the native integrators.
use super::*;
use crate::sim::models_geodesic_relativity as model;
struct Parameters {
    native: model::GeodesicParameters,
    frame: Frame,
}
impl NativeType for Parameters {
    const KIND: &'static str = "rssn/native/relativity-model-parameters/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.black_hole_mass <= 0.0
            || self.native.initial_state[0] <= 0.0
            || self.native.initial_dt <= 0.0
            || self.native.proper_time_end < 0.0
            || self.native.initial_state.iter().any(|x| !x.is_finite())
        {
            return Err(api::error("geodesic model parameters invalid"));
        }
        Ok(())
    }
}
fn work(budget: usize) -> Result<(), OpError> {
    if budget == 0 || budget > 100000 {
        return Err(api::error("geodesic model step/RHS budget mustbe1..100000"));
    }
    Ok(())
}
object_binding!(geodesic_model_parameters,"num_native_geodesic_model_parameters","","Actual native geodesic model parameter state. Mass is geometricLength, radius physicalLength, radial velocityLength/Time, angular rate1/Time; proper time restores seconds through the same geometrized length reference. Every run separately charges actual implicit steps or adaptive RHS work.",false,(f:Handle<Frame>,mass_geometric:Input,radius:Input,radial_velocity:Input,phi:f64,angular_velocity:Input,proper_time_end:Input,initial_dt:Input)->NewObject<Parameters>,|ctx:&CallCtx<'_>,f:Handle<Frame>,mass_geometric:Input,radius:Input,radial_velocity:Input,phi:f64,angular_velocity:Input,proper_time_end:Input,initial_dt:Input|{let f=frame(ctx,&f)?;let mass=f.read(&mass_geometric.0,1)?;let radius=f.read(&radius.0,1)?;let velocity_dimension=if f.physical{DimExpr::Mul(vec![DimExpr::base("Length"),DimExpr::base("Time").pow(-1)])}else{DimExpr::Dimless};api::same_dimension(&api::dimension(&radial_velocity.0)?,&velocity_dimension)?;let radial=crate::rules::numerical_abi::numeric(&radial_velocity.0,0)?/if f.physical{C}else{1.0};let rate_dimension=if f.physical{DimExpr::base("Time").pow(-1)}else{DimExpr::Dimless};api::same_dimension(&api::dimension(&angular_velocity.0)?,&rate_dimension)?;let rate=api::scale_product(crate::rules::numerical_abi::numeric(&angular_velocity.0,0)?,&[f.length],&[if f.physical{C}else{1.0}])?;let end=f.time_in(&proper_time_end.0)?;let dt=f.time_in(&initial_dt.0)?;if mass<=0.0||radius<=2.0*mass||dt<=0.0||end<0.0||end/dt>100000.0{return Err(api::error("geodesic model needsoutsidehorizon positive mass/step and boundednominalsteps"));}let native=model::GeodesicParameters{black_hole_mass:mass,initial_state:[radius,radial,phi,rate],proper_time_end:end,initial_dt:dt};let dim=dims(&f);Ok(NewObject::new(Parameters{native,frame:f},dim))});
object_binding!(geodesic_model_effective_potential,"num_native_geodesic_model_effective_potential","crate::sim::models_geodesic_relativity::GeodesicParameters::effective_potential","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,radius:Input,angular_momentum_geometric:Input)->f64,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,radius:Input,angular_momentum_geometric:Input|{parameters.with(ctx,|p|{let r=p.frame.read(&radius.0,1)?;let l=p.frame.read(&angular_momentum_geometric.0,1)?;if r<=0.0{return Err(api::error("effectivepotential radius mustbepositive"));}Ok(p.native.effective_potential(r,l))})});
object_binding!(geodesic_model_run_geodesic_simulation,"num_native_geodesic_model_run_geodesic_simulation","crate::sim::models_geodesic_relativity::run_geodesic_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{let path=model::with_native_budget(step_rhs_budget,||model::run_geodesic_simulation(&p.native)).map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(geodesic_model_run_geodesic_simulation_rk45,"num_native_geodesic_model_run_geodesic_simulation_rk45","crate::sim::models_geodesic_relativity::run_geodesic_simulation_rk45","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{let path=model::with_native_budget(step_rhs_budget,||model::run_geodesic_simulation_rk45(&p.native)).map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(geodesic_model_run_kerr_newman_geodesic,"num_native_geodesic_model_run_kerr_newman_geodesic","crate::sim::models_geodesic_relativity::run_kerr_newman_geodesic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,spin:Input,charge_geometric:Input,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,spin:Input,charge_geometric:Input,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{let spin=p.frame.read(&spin.0,1)?;let charge=p.frame.read(&charge_geometric.0,1)?;let path=model::with_native_budget(step_rhs_budget,||model::run_kerr_newman_geodesic(&p.native,spin,charge)).map_err(api::error)?.map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(geodesic_model_file_scenario,"num_native_geodesic_model_file_scenario","crate::sim::models_geodesic_relativity::simulate_black_hole_orbits_scenario","IO: actual native reduced-unit Schwarzschild orbit demonstration writes CSV files in the chosen directory. The explicit shared step/RHS budget is inherited by every native Rayon orbit worker and returns a checked error on exhaustion.",true,(output_directory:String,step_rhs_budget:usize)->bool,|_ctx:&CallCtx<'_>,output_directory:String,step_rhs_budget:usize|{work(step_rhs_budget)?;if output_directory.is_empty()||output_directory.len()>4096{return Err(api::error("native orbit scenario needsnonemptyoutputdirectory"));}model::with_native_budget(step_rhs_budget,||model::simulate_black_hole_orbits_scenario(std::path::Path::new(&output_directory))).map_err(api::error)?.map_err(api::error)?;Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &geodesic_model_parameters::BINDING,
            &geodesic_model_effective_potential::BINDING,
            &geodesic_model_run_geodesic_simulation::BINDING,
            &geodesic_model_run_geodesic_simulation_rk45::BINDING,
            &geodesic_model_run_kerr_newman_geodesic::BINDING,
            &geodesic_model_file_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
