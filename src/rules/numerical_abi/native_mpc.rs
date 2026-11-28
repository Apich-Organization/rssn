//! Bounded normalized MPC trajectories with a reusable condensation.
use super::codecs::Arg;
use super::native_control::row_matrix;
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernels::control::mpc::CondensedMpc;
use crate::kernels::control::mpc::MpcOptions;
use crate::kernels::control::state_space::StateSpace;
use crate::kernels::dense::Mat;
use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::field;
use api::object_binding;

const MAX_STATES: usize = 8;
const MAX_INPUTS: usize = 4;
const MAX_PREDICTION_HORIZON: usize = 32;
const MAX_KKT_DIMENSION: usize = 64;
const MAX_STEPS: usize = 256;
const MAX_TRAJECTORY_WORK: usize = 1 << 20;

impl Arg for MpcOptions {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("horizon".into(), TypeExpr::Int),
            ("control_horizon".into(), TypeExpr::Int),
            ("q".into(), Vec::<f64>::ty()),
            ("r".into(), Vec::<f64>::ty()),
            ("u_min".into(), Option::<Vec<f64>>::ty()),
            ("u_max".into(), Option::<Vec<f64>>::ty()),
            ("x_min".into(), Option::<Vec<f64>>::ty()),
            ("x_max".into(), Option::<Vec<f64>>::ty()),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Record(record) = value else {
            return Err(error("expected MPC options"));
        };
        if record.fields.len() != 8 {
            return Err(error(
                "expected horizons, weights and four optional MPC bounds",
            ));
        }
        Ok(Self {
            horizon: usize::decode(field(value, "horizon")?, index)?,
            control_horizon: usize::decode(field(value, "control_horizon")?, index)?,
            q: Vec::<f64>::decode(field(value, "q")?, index)?,
            r: Vec::<f64>::decode(field(value, "r")?, index)?,
            u_min: Option::<Vec<f64>>::decode(field(value, "u_min")?, index)?,
            u_max: Option::<Vec<f64>>::decode(field(value, "u_max")?, index)?,
            x_min: Option::<Vec<f64>>::decode(field(value, "x_min")?, index)?,
            x_max: Option::<Vec<f64>>::decode(field(value, "x_max")?, index)?,
        })
    }
}

struct PreparedMpc {
    problem: CondensedMpc,
    states: usize,
    max_kkt: usize,
}

impl NativeType for PreparedMpc {
    const KIND: &'static str = "rssn/native/control/mpc/v1";
    // Includes both prediction maps, weighted Hessian, constraint rows and
    // their Vec storage at the configured shape limits.
    const MAX_BYTES: usize = 256 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

fn rows(matrix: &Mat) -> Vec<Vec<f64>> {
    matrix
        .data
        .chunks_exact(matrix.cols)
        .map(<[f64]>::to_vec)
        .collect()
}

fn mpc_error(message: String) -> OpError {
    crate::cancellation::checkpoint().map_or_else(OpError::from, |()| error(message))
}

object_binding!(
    mpc_prepare,
    "num_control_mpc_prepare",
    "",
    "Prepares a bounded normalized discrete MPC from state and input matrices, positive sample time, horizons, weights and optional ordered bounds. Each input weight must be positive.",
    true,
    (a: Mat, b: Mat, dt: f64, options: MpcOptions) -> NewObject<PreparedMpc>,
    |_, a, b, dt, options| {
        let states = a.rows;
        let inputs = b.cols;
        if !(1..=MAX_STATES).contains(&states) || !(1..=MAX_INPUTS).contains(&inputs)
            || a.cols != states || b.rows != states
            || !(1..=MAX_PREDICTION_HORIZON).contains(&options.horizon)
            || !(1..=options.horizon).contains(&options.control_horizon)
            || options.r.iter().any(|weight| *weight <= 0.0)
        {
            return Err(error("MPC shape, horizons or input weights exceed the supported range"));
        }
        let variables = inputs.saturating_mul(options.control_horizon);
        let predictions = states.saturating_mul(options.horizon);
        let count_input_bounds = usize::from(options.u_min.is_some()).saturating_add(usize::from(options.u_max.is_some()));
        let count_state_bounds = usize::from(options.x_min.is_some()).saturating_add(usize::from(options.x_max.is_some()));
        let max_kkt = variables.saturating_mul(count_input_bounds.saturating_add(1))
            .saturating_add(predictions.saturating_mul(count_state_bounds));
        if variables > 32 || max_kkt > MAX_KKT_DIMENSION {
            return Err(error("MPC constraints exceed the condensed QP work budget"));
        }
        let model = StateSpace::discrete(rows(&a), rows(&b), vec![], vec![], dt).map_err(error)?;
        let problem = CondensedMpc::new(&model, &options).map_err(mpc_error)?;
        Ok(NewObject::normalized(PreparedMpc { problem, states, max_kkt }))
    }
);

object_binding!(
    mpc_trajectory,
    "num_control_mpc_closed_loop",
    "crate::kernels::control::mpc::mpc_closed_loop",
    "Runs a normalized MPC from an initial state toward a fixed state reference, reusing the condensed QP for every step. The trajectory includes the initial state.",
    true,
    (plan: Handle<PreparedMpc>, initial: Vec<f64>, reference: Vec<f64>, steps: usize) -> Mat,
    |ctx, plan, initial, reference, steps| plan.with(ctx, |plan| {
        if steps > MAX_STEPS || steps.saturating_mul(plan.max_kkt.saturating_pow(3)) > MAX_TRAJECTORY_WORK {
            return Err(error("MPC trajectory exceeds the work budget"));
        }
        let values = plan.problem.closed_loop(&initial, &reference, steps).map_err(mpc_error)?;
        Ok(row_matrix(values, plan.states))
    })
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps MPC registrations with their domain contracts"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[&mpc_prepare::BINDING, &mpc_trajectory::BINDING],
        "numerical_vector",
        DimRule::Dimensionless,
    )
}
