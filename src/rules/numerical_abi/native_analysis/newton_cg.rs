//! Native tape-based Newton-CG preserving mixed coordinate reference domains.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::value::TensorValue;
use crate::kernels::newton_cg::NewtonCgOptions;
use crate::kernels::newton_cg::NewtonCgWorkspace;
use crate::kernels::newton_cg::{
    self,
};
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
pub(super) fn parameters(
    v: &Input,
    scales: &[f64],
    dimensions: &[crate::kernel::signature::DimExpr],
) -> Result<Vec<f64>, OpError> {
    let Value::Tensor(t) = &v.0 else {
        return Err(api::error(
            "native tape arguments require a rank-one tensor",
        ));
    };
    if t.shape != [scales.len()] || t.data.len() != scales.len() || scales.len() > 4096 {
        return Err(api::error(
            "native tape argument size or4096component budget mismatch",
        ));
    }
    // Read the actual uniform element metadata without forcing independently
    // wrapped Value cells to share one physical dimension. Each cell is checked below.
    let uniform = api::dimension(&Value::Tensor(TensorValue {
        shape: vec![0],
        elem: t.elem.clone(),
        data: Vec::new(),
    }))?;
    t.data
        .iter()
        .zip(scales)
        .zip(dimensions)
        .map(|((v, s), d)| {
            let actual = if matches!(v, Value::Dyn(_) | Value::Quantity(_)) {
                api::dimension(v)?
            } else {
                uniform.clone()
            };
            api::same_dimension(d, &actual)?;
            api::scale_product(super::numeric(v, 0)?, &[], &[*s])
        })
        .collect()
}
fn options(v: &Input) -> Result<NewtonCgOptions, OpError> {
    api::normalized(&v.0)?;
    let n = |field| <usize as Arg>::decode(api::field(&v.0, field)?, 0);
    let f = |field| <f64 as Arg>::decode(api::field(&v.0, field)?, 0);
    Ok(NewtonCgOptions {
        gradient_tolerance: f("gradient_tolerance")?,
        max_iterations: n("max_iterations")?,
        max_cg_iterations: n("max_cg_iterations")?,
        cg_relative_tolerance: f("cg_relative_tolerance")?,
        initial_radius: f("initial_radius")?,
        max_radius: f("max_radius")?,
        min_radius: f("min_radius")?,
        acceptance_threshold: f("acceptance_threshold")?,
    })
}
struct PreparedNewtonCg {
    snapshot: TapeSnapshot,
    output: usize,
    workspace: NewtonCgWorkspace,
}

impl NativeType for PreparedNewtonCg {
    const KIND: &'static str = "rssn/native/analysis/newton-cg-tape-plan/v1";
    // Reserve the owned tape budget and the largest workspace that the
    // checked 4096-coordinate entry point can grow into across repeated calls.
    const MAX_BYTES: usize = 5 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        check_output(&self.snapshot, self.output)
    }
}

fn check_output(
    snapshot: &TapeSnapshot,
    output: usize,
) -> Result<(), OpError> {
    if output >= snapshot.output_scales.len()
        || snapshot.output_dimensions.len() != snapshot.output_scales.len()
    {
        return Err(api::error("native tape output index out of range"));
    }
    Ok(())
}

fn solve_tape(
    snapshot: &TapeSnapshot,
    output: usize,
    initial: &Input,
    params: &Input,
    config: &Input,
    workspace: Option<&mut NewtonCgWorkspace>,
) -> Result<Dynamic, OpError> {
    check_output(snapshot, output)?;
    let initial = parameters(initial, &snapshot.input_scales, &snapshot.input_dimensions)?;
    let params = parameters(
        params,
        &snapshot.parameter_scales,
        &snapshot.parameter_dimensions,
    )?;
    let options = options(config)?;
    if options.max_cg_iterations > 4096
        || options
            .max_iterations
            .checked_mul(options.max_cg_iterations.saturating_add(8))
            .and_then(|n| n.checked_mul(snapshot.tape.insts.len().max(1)))
            .is_none_or(|n| n > 100_000_000)
    {
        return Err(api::error("Newton-CG tape work budget exceeded"));
    }
    let tape = &snapshot.tape;
    let r = workspace
        .map_or_else(
            || newton_cg::minimize_tape(tape, output, &initial, &params, &options),
            |workspace| workspace.minimize_tape(tape, output, &initial, &params, &options),
        )
        .map_err(api::error)?;
    let x = Value::Tensor(TensorValue {
        shape: vec![r.x.len()],
        elem: api::dynamic_type(),
        data: r
            .x
            .into_iter()
            .zip(&snapshot.input_scales)
            .zip(&snapshot.input_dimensions)
            .map(|((v, s), d)| api::quantity(api::scale_product(v, &[*s], &[])?, d.clone()))
            .collect::<Result<_, _>>()?,
    });
    let value = api::scale_product(r.value, &[snapshot.output_scales[output]], &[])?;
    Ok(Dynamic(api::record(vec![
        ("x", x),
        (
            "value",
            api::quantity(value, snapshot.output_dimensions[output].clone())?,
        ),
        ("normalized_gradient_norm", Output::encode(r.gradient_norm)?),
        ("status", Value::String(format!("{:?}", r.status))),
        ("iterations", Output::encode(r.iterations)?),
        (
            "objective_evaluations",
            Output::encode(r.objective_evaluations)?,
        ),
        (
            "gradient_evaluations",
            Output::encode(r.gradient_evaluations)?,
        ),
        (
            "hessian_vector_evaluations",
            Output::encode(r.hessian_vector_evaluations)?,
        ),
        (
            "negative_curvature_steps",
            Output::encode(r.negative_curvature_steps)?,
        ),
        ("boundary_steps", Output::encode(r.boundary_steps)?),
        ("accepted_steps", Output::encode(r.accepted_steps)?),
    ])))
}

object_binding!(
    minimize,
    "num_analysis_newton_cg_tape",
    "crate::kernels::newton_cg::minimize_tape",
    "Actual native Tape Newton-CG. Tape input, parameter and output references preserve mixed physical coordinate units. Initial and parameter vectors may hold separately measured Value cells. Options use the normalized coordinate metric. The normalized gradient norm is dimensionless. Iteration, CG-action and instruction work is bounded before the solve.",
    false,
    (tape: Input, output: usize, initial: Input, params: Input, config: Input) -> Dynamic,
    |ctx, tape, output, initial, params, config| {
        let snapshot = tape_snapshot(ctx, &tape.0)?;
        solve_tape(&snapshot, output, &initial, &params, &config, None)
    }
);

object_binding!(
    prepare,
    "num_analysis_newton_cg_tape_plan_new",
    "",
    "Retain an independent snapshot of a native tape and select one objective output. The plan also holds the Newton and CG buffers used by repeated solves. Releasing the input tape does not invalidate the plan.",
    false,
    (tape: Input, output: usize) -> NewObject<PreparedNewtonCg>,
    |ctx, tape, output| {
        let snapshot = tape_snapshot(ctx, &tape.0)?;
        Ok(NewObject::normalized(PreparedNewtonCg {
            snapshot,
            output,
            workspace: NewtonCgWorkspace::default(),
        }))
    }
);

object_binding!(
    prepared_minimize,
    "num_analysis_newton_cg_tape_plan_minimize",
    "crate::kernels::newton_cg::NewtonCgWorkspace::minimize_tape",
    "Solve a prepared tape objective from a new measured initial point and parameters, using the plan's retained buffers. Each solve has its own normalized options and work budget. Results and counters have the same units and meaning as the one-shot tape solver.",
    false,
    (plan: Handle<PreparedNewtonCg>, initial: Input, params: Input, config: Input) -> Dynamic,
    |ctx, plan, initial, params, config| plan.with_mut(ctx, |plan| {
        solve_tape(&plan.snapshot, plan.output, &initial, &params, &config, Some(&mut plan.workspace))
    })
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &minimize::BINDING,
            &prepare::BINDING,
            &prepared_minimize::BINDING,
        ],
        "numerical_vector",
    )
}
