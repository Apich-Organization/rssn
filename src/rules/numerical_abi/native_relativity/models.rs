#![allow(clippy::wildcard_imports)] // Geometrized model parameters use the same measured relativity frame as the native integrators.
use super::*;
use crate::sim::models_geodesic_relativity as model;
use std::num::NonZeroUsize;

struct Parameters {
    native: model::GeodesicParameters,
    frame: Frame,
}
impl NativeType for Parameters {
    const KIND: &'static str = "rssn/native/relativity-model-parameters/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        let p = &self.native;
        if p.black_hole_mass <= 0.0
            || p.initial_state[0] <= 0.0
            || p.initial_dt <= 0.0
            || p.proper_time_end < 0.0
            || p.proper_time_end / p.initial_dt > 100_000.0
            || ![p.black_hole_mass, p.initial_dt, p.proper_time_end]
                .into_iter()
                .all(f64::is_finite)
            || p.initial_state.iter().any(|x| !x.is_finite())
        {
            return Err(api::error(
                "geodesic parameters need positive mass, radius and step, with a finite nonnegative span and at most 100,000 nominal steps",
            ));
        }
        Ok(())
    }
}
impl Parameters {
    fn metric(
        &self,
        spin: &Input,
        charge: &Input,
    ) -> Result<r::KerrNewman, OpError> {
        Ok(r::KerrNewman {
            m: self.native.black_hole_mass,
            a: self.frame.read(&spin.0, 1)?,
            q: self.frame.read(&charge.0, 1)?,
        })
    }

    fn require_schwarzschild_exterior(&self) -> Result<(), OpError> {
        if self.native.initial_state[0] <= 2.0 * self.native.black_hole_mass {
            return Err(api::error(
                "initial radius must be outside the Schwarzschild horizon",
            ));
        }
        Ok(())
    }
}

const fn energy_branch(nonpositive: bool) -> model::EquatorialEnergyBranch {
    if nonpositive {
        model::EquatorialEnergyBranch::Nonpositive
    } else {
        model::EquatorialEnergyBranch::Nonnegative
    }
}

fn work(budget: usize) -> Result<(), OpError> {
    if budget == 0 || budget > 100000 {
        return Err(api::error("geodesic model step/RHS budget mustbe1..100000"));
    }
    Ok(())
}
object_binding!(
    geodesic_model_parameters, "num_native_geodesic_model_parameters", "",
    "Prepare equatorial geodesic parameters. In a physical frame, mass and radius have length units, radial velocity has length per time, angular velocity is inverse time, and simulation times are seconds. Positive radius is checked here. The selected spacetime determines the horizon at initialization or integration.",
    false,
    (f: Handle<Frame>, mass_geometric: Input, radius: Input, radial_velocity: Input, phi: f64, angular_velocity: Input, proper_time_end: Input, initial_dt: Input) -> NewObject<Parameters>,
    |ctx:&CallCtx<'_>, f:Handle<Frame>, mass_geometric:Input, radius:Input, radial_velocity:Input, phi:f64, angular_velocity:Input, proper_time_end:Input, initial_dt:Input| {
        let f = frame(ctx, &f)?;
        let mass = f.read(&mass_geometric.0, 1)?;
        let radius = f.read(&radius.0, 1)?;
        let (velocity_dim, rate_dim, speed) = if f.physical {
            let rate = DimExpr::base("Time").pow(-1);
            (DimExpr::Mul(vec![DimExpr::base("Length"), rate.clone()]), rate, C)
        } else {
            (DimExpr::Dimless, DimExpr::Dimless, 1.0)
        };
        api::same_dimension(&api::dimension(&radial_velocity.0)?, &velocity_dim)?;
        api::same_dimension(&api::dimension(&angular_velocity.0)?, &rate_dim)?;
        let radial = crate::rules::numerical_abi::numeric(&radial_velocity.0, 0)? / speed;
        let rate = api::scale_product(
            crate::rules::numerical_abi::numeric(&angular_velocity.0, 0)?,
            &[f.length],
            &[speed],
        )?;
        let native = model::GeodesicParameters {
            black_hole_mass: mass,
            initial_state: [radius, radial, phi, rate],
            proper_time_end: f.time_in(&proper_time_end.0)?,
            initial_dt: f.time_in(&initial_dt.0)?,
        };
        let dimensions = dims(&f);
        Ok(NewObject::new(Parameters { native, frame: f }, dimensions))
    }
);
object_binding!(
    geodesic_model_effective_potential, "num_native_geodesic_model_effective_potential",
    "crate::sim::models_geodesic_relativity::GeodesicParameters::effective_potential",
    "Evaluates the equatorial Schwarzschild effective potential for the retained mass. Supply a positive radius and specific angular momentum, both as geometric lengths in a physical frame.\n\nReturns the dimensionless potential −M/r + (1/2 − M/r)(L/r)^2. The scaled calculation preserves finite ratios even when individual length powers fall outside Float.\n\nInvalid units, nonpositive radius and unrepresentable results produce errors.",
    false, (parameters:Handle<Parameters>,radius:Input,angular_momentum_geometric:Input) -> f64,
    |ctx,parameters,radius,angular_momentum_geometric| parameters.with(ctx, |p| {
        let radius = p.frame.read(&radius.0, 1)?;
        let angular = p.frame.read(&angular_momentum_geometric.0, 1)?.abs();
        if radius <= 0.0 { return Err(api::error("effective potential radius must be positive")); }
        let mass = p.native.black_hole_mass;
        let mass_ratio = api::scale_product(mass, &[], &[radius])?;
        if angular == 0.0 { return Ok(-mass_ratio); }
        let coefficient = if (0.25..=0.75).contains(&mass_ratio) {
            // Subtract before dividing near the horizon. Form r - 2M rather
            // than r/2 - M to preserve differences on a subnormal grid.
            ((-2.0_f64).mul_add(mass, radius) / radius) * 0.5
        } else {
            0.5 - mass_ratio
        };
        let orbital = api::scale_product(coefficient, &[angular, angular], &[radius, radius])?;
        Ok(orbital - mass_ratio)
    })
);
object_binding!(geodesic_model_run_geodesic_simulation,"num_native_geodesic_model_run_geodesic_simulation","crate::sim::models_geodesic_relativity::run_geodesic_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{p.require_schwarzschild_exterior()?;let path=model::with_native_budget(step_rhs_budget,||model::run_geodesic_simulation(&p.native)).map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(geodesic_model_run_geodesic_simulation_rk45,"num_native_geodesic_model_run_geodesic_simulation_rk45","crate::sim::models_geodesic_relativity::run_geodesic_simulation_rk45","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{p.require_schwarzschild_exterior()?;let path=model::with_native_budget(step_rhs_budget,||model::run_geodesic_simulation_rk45(&p.native)).map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(geodesic_model_run_kerr_newman_geodesic,"num_native_geodesic_model_run_kerr_newman_geodesic","crate::sim::models_geodesic_relativity::run_kerr_newman_geodesic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Handle<Parameters>,spin:Input,charge_geometric:Input,step_rhs_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,parameters:Handle<Parameters>,spin:Input,charge_geometric:Input,step_rhs_budget:usize|{work(step_rhs_budget)?;parameters.with(ctx,|p|{let spin=p.frame.read(&spin.0,1)?;let charge=p.frame.read(&charge_geometric.0,1)?;let path=model::with_native_budget(step_rhs_budget,||model::run_kerr_newman_geodesic(&p.native,spin,charge)).map_err(api::error)?.map_err(api::error)?;curve(path,&p.frame)})});
object_binding!(
    geodesic_model_equatorial_timelike, "num_native_geodesic_model_equatorial_timelike",
    "crate::sim::models_geodesic_relativity::equatorial_timelike_geodesic",
    "Prepare a future-directed equatorial timelike Kerr–Newman geodesic. Spin and charge are geometric lengths in a physical frame. Set nonpositive_energy for the second branch available for compatible states inside the stationary limit and outside the horizon. Use the relativity state operator to read energy, angular momentum and the initial state.",
    false,
    (parameters: Handle<Parameters>, spin: Input, charge_geometric: Input, nonpositive_energy: bool) -> NewObject<Geodesic>,
    |ctx, parameters, spin, charge_geometric, nonpositive_energy| {
        parameters.with(ctx, |p| {
            let native = model::equatorial_timelike_geodesic(
                p.metric(&spin, &charge_geometric)?,
                p.native.initial_state,
                energy_branch(nonpositive_energy),
            ).map_err(api::error)?;
            let frame = p.frame.clone();
            Ok(NewObject::new(Geodesic { native, frame }, dims(&p.frame)))
        })
    }
);
object_binding!(
    geodesic_model_run_kerr_newman_sampled, "num_native_geodesic_model_run_kerr_newman_geodesic_sampled",
    "crate::sim::models_geodesic_relativity::run_kerr_newman_geodesic_sampled",
    "Return a two-column equatorial Cartesian orbit in the frame's length units. Retain the initial point, every specified positive number of steps, and the final or termination point. Sampling leaves the implicit-step budget and spacetime termination rules unchanged. nonpositive_energy requests the second future-directed branch.",
    false,
    (parameters: Handle<Parameters>, spin: Input, charge_geometric: Input, nonpositive_energy: bool, retain_every: usize, step_rhs_budget: usize) -> Dynamic,
    |ctx, parameters, spin, charge_geometric, nonpositive_energy, retain_every, step_rhs_budget| {
        let stride = NonZeroUsize::new(retain_every)
            .ok_or_else(|| api::error("geodesic sampling stride must be positive"))?;
        work(step_rhs_budget)?;
        parameters.with(ctx, |p| {
            let metric = p.metric(&spin, &charge_geometric)?;
            let path = model::with_native_budget(step_rhs_budget, || {
                model::run_kerr_newman_geodesic_sampled(
                    &p.native, metric.a, metric.q, energy_branch(nonpositive_energy), stride,
                )
            }).map_err(api::error)?.map_err(api::error)?;
            curve(path, &p.frame)
        })
    }
);
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
            &geodesic_model_equatorial_timelike::BINDING,
            &geodesic_model_run_kerr_newman_sampled::BINDING,
            &geodesic_model_file_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
