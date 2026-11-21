//! Native method-of-lines bindings for spatial and time-dependent endpoints.
use super::BoundaryObject;
use super::CallCtx;
use super::CompiledSnapshot;
use super::FieldObject;
use super::Frame;
use super::Grid;
use super::GridObject;
use super::Handle;
use super::Input;
use super::NewObject;
use super::OpError;
use super::api;
use super::boundary_list;
use super::compiled_snapshot;
use super::object_binding;
use super::p;
use super::program;
use crate::kernels::ode_adaptive::OdeOptions;

struct MolInputs {
    end: Input,
    method: String,
    rtol: f64,
    atol: Input,
    h0: Input,
    hmin: Input,
    hmax: Input,
    max_steps: usize,
    budget: usize,
    rhs: Input,
    initial: Input,
}

#[derive(Clone, Copy)]
enum BoundaryInputs<'a> {
    Fixed(&'a Input),
    Timed(&'a [String; 2], [&'a Input; 2]),
}

struct TimedProgram {
    tape: CompiledSnapshot,
    gradient: bool,
}

impl TimedProgram {
    fn compile(
        ctx: &CallCtx<'_>,
        kind: &str,
        value: &Input,
        frame: &Frame,
        budget: usize,
    ) -> Result<Self, OpError> {
        let gradient = match kind {
            | "Dirichlet" => false,
            | "Neumann" => true,
            | _ => {
                return Err(api::error(
                    "method-of-lines endpoint must be Dirichlet or Neumann",
                ));
            },
        };
        let tape = compiled_snapshot(ctx, &value.0, budget)?;
        let [input] = tape.input_dimensions.as_slice() else {
            return Err(api::error(
                "timed boundary requires one measured Time input",
            ));
        };
        api::same_dimension(input, &frame.unit(0, 1, 0))?;
        api::same_dimension(
            &tape.output_dimension,
            &frame.unit(-i64::from(gradient), 0, 1),
        )?;
        Ok(Self { tape, gradient })
    }

    fn sample(
        &self,
        time: f64,
        frame: &Frame,
    ) -> f64 {
        let Ok(time) = api::scale_product(time, &[frame.time], &[]) else {
            return f64::NAN;
        };
        let result = self.tape.sample(&[time]);
        let scale = [frame.length];
        let numerator: &[f64] = if self.gradient { &scale } else { &[] };
        api::scale_product(result, numerator, &[frame.value]).unwrap_or(f64::NAN)
    }

    fn boundary<'a>(
        &self,
        sample: &'a dyn Fn(f64) -> f64,
    ) -> p::TimedBoundary1D<'a> {
        if self.gradient {
            p::TimedBoundary1D::Neumann(sample)
        } else {
            p::TimedBoundary1D::Dirichlet(sample)
        }
    }
}

#[allow(
    clippy::large_enum_variant,
    reason = "retain two endpoint tapes on the stack for one integration"
)]
enum Conditions {
    Fixed(Vec<BoundaryObject>),
    Timed([TimedProgram; 2]),
}

impl Conditions {
    fn finish(&self) -> Result<(), OpError> {
        match self {
            | Self::Fixed(boundaries) => {
                for snapshot in boundaries.iter().filter_map(|b| b.snapshot.as_ref()) {
                    snapshot.finish()?;
                }
            },
            | Self::Timed(boundaries) => {
                for condition in boundaries {
                    condition.tape.finish()?;
                }
            },
        }
        Ok(())
    }
}

fn solve(
    ctx: &CallCtx<'_>,
    grid: &Handle<GridObject>,
    boundary: BoundaryInputs<'_>,
    inputs: &MolInputs,
) -> Result<NewObject<FieldObject>, OpError> {
    let grid = grid.with(ctx, |grid| Ok(grid.clone()))?;
    let Grid::D1(native_grid) = grid.native else {
        return Err(api::error("method of lines requires rank1 grid"));
    };
    if inputs.max_steps == 0
        || inputs
            .max_steps
            .checked_mul(grid.count())
            .is_none_or(|n| n > 4_000_000)
    {
        return Err(api::error("method-of-lines work budget exceeded"));
    }
    let frame = &grid.frame;
    let budget = inputs.budget;
    let conditions = match boundary {
        | BoundaryInputs::Fixed(values) => Conditions::Fixed(boundary_list(ctx, values, &grid)?),
        | BoundaryInputs::Timed(kinds, values) => {
            let [left_kind, right_kind] = kinds;
            let [left, right] = values;
            Conditions::Timed([
                TimedProgram::compile(ctx, left_kind, left, frame, budget)?,
                TimedProgram::compile(ctx, right_kind, right, frame, budget)?,
            ])
        },
    };
    let tape = compiled_snapshot(ctx, &inputs.rhs.0, budget)?;
    let expected = [
        frame.unit(0, 1, 0),
        frame.unit(1, 0, 0),
        frame.unit(0, 0, 1),
        frame.unit(-1, 0, 1),
        frame.unit(-2, 0, 1),
    ];
    if tape.input_dimensions.len() != expected.len() {
        return Err(api::error(
            "method-of-lines RHS requires five physical inputs",
        ));
    }
    for (actual, expected) in tape.input_dimensions.iter().zip(&expected) {
        api::same_dimension(actual, expected)?;
    }
    api::same_dimension(&tape.output_dimension, &frame.unit(0, -1, 1))?;
    let initial = program(ctx, &inputs.initial, frame, 1, frame.unit(0, 0, 1), budget)?;
    let end = frame.measured(&inputs.end, 0, 1, 0)?;
    let atol = frame.measured(&inputs.atol, 0, 0, 1)?;
    let h0 = frame.measured(&inputs.h0, 0, 1, 0)?;
    let hmin = frame.measured(&inputs.hmin, 0, 1, 0)?;
    let hmax = frame.measured(&inputs.hmax, 0, 1, 0)?;
    if inputs.rtol < 0.0 || atol <= 0.0 || h0 <= 0.0 || hmin < 0.0 || hmax <= 0.0 || hmin > hmax {
        return Err(api::error(
            "invalid method-of-lines tolerance or time bounds",
        ));
    }
    let method = match inputs.method.as_str() {
        | "Radau5" => p::StiffMethod::Radau5,
        | "Bdf1" => p::StiffMethod::Bdf(1),
        | "Bdf2" => p::StiffMethod::Bdf(2),
        | "Bdf3" => p::StiffMethod::Bdf(3),
        | "Bdf4" => p::StiffMethod::Bdf(4),
        | "Bdf5" => p::StiffMethod::Bdf(5),
        | _ => return Err(api::error("method must be Radau5 or Bdf1..Bdf5")),
    };
    let opts = OdeOptions {
        rtol: inputs.rtol,
        atol,
        h0: Some(h0),
        hmin,
        hmax,
        max_steps: inputs.max_steps,
    };
    let l = frame.length;
    let t = frame.time;
    let v = frame.value;
    let rhs = |tn: f64, xn: f64, un: f64, ux: f64, uxx: f64| {
        tape.sample(&[tn * t, xn * l, un * v, ux * v / l, uxx * v / l / l]) / (v / t)
    };
    let u0 = |xn: f64| initial.sample(&[xn * l]) / v;
    let result = match &conditions {
        | Conditions::Fixed(boundaries) => {
            let [left, right] = boundaries.as_slice() else {
                return Err(api::error("method-of-lines boundary count"));
            };
            let bc = [left.native.clone(), right.native.clone()];
            p::method_of_lines_1d(&rhs, &native_grid, &bc, &u0, end, method, &opts)
        },
        | Conditions::Timed([left, right]) => {
            let left_sample = |time| left.sample(time, frame);
            let right_sample = |time| right.sample(time, frame);
            let bc = [left.boundary(&left_sample), right.boundary(&right_sample)];
            p::method_of_lines_1d_timed(&rhs, &native_grid, &bc, &u0, end, method, &opts)
        },
    };
    tape.finish()?;
    initial.finish()?;
    conditions.finish()?;
    let result = result.map_err(api::error)?;
    Ok(NewObject::new(
        FieldObject {
            native: p::PdeSolution::D1(result.field),
            frame: frame.clone(),
            work: Some((result.steps, result.nfev)),
        },
        frame.dims(),
    ))
}

object_binding!(
    fixed, "num_analysis_pde_method_of_lines_1d", "crate::kernels::pde::method_of_lines_1d",
    "Native stiff1D method of lines. Compiled RHS inputs have domains(time,length,value,value/length,value/length^2),output value/time. Initial program is length->value. Method Radau5 or Bdf1..Bdf5; physical atol and all step bounds retain the grid time/value domains; explicit attempted-step and per-program evaluation budgets. field_data exposes native accepted steps/evaluations.",
    false,
    (grid: Handle<GridObject>, boundaries: Input, end: Input, method: String, rtol: f64, atol: Input, h0: Input, hmin: Input, hmax: Input, max_steps: usize, budget: usize, rhs: Input, initial: Input) -> NewObject<FieldObject>,
    |ctx, grid, boundaries, end, method, rtol, atol, h0, hmin, hmax, max_steps, budget, rhs, initial| {
        solve(ctx, &grid, BoundaryInputs::Fixed(&boundaries), &MolInputs { end, method, rtol, atol, h0, hmin, hmax, max_steps, budget, rhs, initial })
    }
);
object_binding!(
    timed, "num_analysis_pde_method_of_lines_1d_timed", "crate::kernels::pde::method_of_lines_1d_timed",
    "Stiff 1D method of lines with two time-dependent compiled boundary programs in left-to-right order. Set each kind to Dirichlet for a measured field value or Neumann for a measured outward derivative. Each boundary program takes a measured Time. The returned field retains physical samples and exact prescribed Dirichlet endpoint values at the requested final time.",
    false,
    (grid: Handle<GridObject>, kinds: [String; 2], left: Input, right: Input, end: Input, method: String, rtol: f64, atol: Input, h0: Input, hmin: Input, hmax: Input, max_steps: usize, budget: usize, rhs: Input, initial: Input) -> NewObject<FieldObject>,
    |ctx, grid, kinds, left, right, end, method, rtol, atol, h0, hmin, hmax, max_steps, budget, rhs, initial| {
        solve(ctx, &grid, BoundaryInputs::Timed(&kinds, [&left, &right]), &MolInputs { end, method, rtol, atol, h0, hmin, hmax, max_steps, budget, rhs, initial })
    }
);

pub(super) const FIXED: &api::Binding = &fixed::BINDING;
pub(super) const TIMED: &api::Binding = &timed::BINDING;
