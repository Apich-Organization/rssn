//! Checked normalized PID trajectories and eigensystem realisation.
use super::codecs::Arg;
use super::codecs::Output;
use super::native_object_api as api;
use super::native_transfer::TransferPlan;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernels::control::ident::EraModel;
use crate::kernels::control::ident::era;
use crate::kernels::control::ident::impulse_response;
use crate::kernels::control::pid::Pid;
use crate::kernels::control::pid::pid_response;
use crate::kernels::control::pid::simulate_pid_checked;
use crate::kernels::dense::Mat;
use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::field;
use api::object_binding;

const MAX_PID_SAMPLES: usize = 16_384;
const MAX_PID_WORK: usize = 1 << 24;
const MAX_ERA_SAMPLES: usize = 384;
const MAX_ERA_OUTPUTS: usize = 16;
const MAX_ERA_ORDER: usize = 24;
const MAX_HANKEL_CELLS: usize = 32_768;
const MAX_IMPULSE_STEPS: usize = 4_096;

impl Arg for Pid {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("kp".into(), TypeExpr::Float),
            ("ki".into(), TypeExpr::Float),
            ("kd".into(), TypeExpr::Float),
            ("filter_time".into(), <Option<f64> as Arg>::ty()),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Record(record) = value else {
            return Err(error("expected PID gain fields"));
        };
        if record.fields.len() != 4 {
            return Err(error("PID gains require kp, ki, kd, and filter_time"));
        }
        Ok(Self {
            kp: f64::decode(field(value, "kp")?, index)?,
            ki: f64::decode(field(value, "ki")?, index)?,
            kd: f64::decode(field(value, "kd")?, index)?,
            tf: Option::<f64>::decode(field(value, "filter_time")?, index)?,
        })
    }
}

struct OutputLimits {
    low: f64,
    high: f64,
}
impl Arg for OutputLimits {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("low".into(), TypeExpr::Float),
            ("high".into(), TypeExpr::Float),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Record(record) = value else {
            return Err(error("expected PID output limits"));
        };
        if record.fields.len() != 2 {
            return Err(error("PID output limits require low and high"));
        }
        Ok(Self {
            low: f64::decode(field(value, "low")?, index)?,
            high: f64::decode(field(value, "high")?, index)?,
        })
    }
}

object_binding!(
    pid_closed_loop,
    "num_control_pid_response",
    "crate::kernels::control::pid::pid_response",
    "Simulates a parallel PID against a prepared proper transfer model on a finite increasing normalized time grid from rest.",
    true,
    (pid: Pid, plant: Handle<TransferPlan>, time: Vec<f64>, setpoint: f64) -> Vec<f64>,
    |ctx, pid, plant, time, setpoint| plant.with(ctx, |plant| {
        let order = plant.model.order();
        if time.len() > MAX_PID_SAMPLES
            || time.len().saturating_mul(order.saturating_mul(order)) > MAX_PID_WORK
        {
            return Err(error("PID plant response exceeds the work budget"));
        }
        pid_response(&pid, &plant.model, &time, setpoint).map_err(error)
    })
);

object_binding!(
    pid_sampled,
    "num_control_simulate_pid",
    "crate::kernels::control::pid::simulate_pid",
    "Simulates a positional parallel PID over finite normalized error samples using a positive step, optional derivative filtering and bounded-output anti-windup.",
    true,
    (pid: Pid, errors: Vec<f64>, dt: f64, limits: Option<OutputLimits>) -> Vec<f64>,
    |_, pid, errors, dt, limits| {
        if errors.len() > MAX_PID_SAMPLES {
            return Err(error("PID sample count exceeds the work limit"));
        }
        simulate_pid_checked(&pid, &errors, dt, limits.map(|pair| (pair.low, pair.high))).map_err(error)
    }
);

impl NativeType for EraModel {
    const KIND: &'static str = "rssn/native/control/era/v1";
    const MAX_BYTES: usize = 16 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

pub(super) fn row_matrix(
    rows: Vec<Vec<f64>>,
    cols: usize,
) -> Mat {
    Mat {
        rows: rows.len(),
        cols,
        data: rows.into_iter().flatten().collect(),
    }
}

struct EraMatrices {
    a: Mat,
    b: Mat,
    c: Mat,
    order: usize,
}
impl Output for EraMatrices {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("a".into(), <Mat as Output>::ty()),
            ("b".into(), <Mat as Output>::ty()),
            ("c".into(), <Mat as Output>::ty()),
            ("order".into(), TypeExpr::Int),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("a", self.a.encode()?),
            ("b", self.b.encode()?),
            ("c", self.c.encode()?),
            ("order", self.order.encode()?),
        ]))
    }
}

object_binding!(
    era_fit,
    "num_control_era",
    "crate::kernels::control::ident::era",
    "Fits a normalized discrete one-input state-space model from finite input and output samples. Output matrix rows are sample times. The first input must be nonzero. The fitted model has no direct feedthrough.",
    true,
    (outputs: Mat, inputs: Vec<f64>, order: usize) -> NewObject<EraModel>,
    |_, outputs, inputs, order| {
        if outputs.rows != inputs.len() || outputs.rows > MAX_ERA_SAMPLES
            || !(1..=MAX_ERA_OUTPUTS).contains(&outputs.cols)
            || !(1..=MAX_ERA_ORDER).contains(&order)
            || order > outputs.rows.saturating_sub(1) / 2
        {
            return Err(error("ERA sample shape or order exceeds the supported range"));
        }
        let columns = outputs.rows.saturating_sub(order).saturating_sub(1);
        if outputs.cols.saturating_mul(order).saturating_mul(columns) > MAX_HANKEL_CELLS {
            return Err(error("ERA Hankel matrix exceeds the work limit"));
        }
        let rows = outputs.data.chunks_exact(outputs.cols).map(<[f64]>::to_vec).collect::<Vec<_>>();
        Ok(NewObject::normalized(era(&rows, &inputs, order).map_err(error)?))
    }
);

object_binding!(
    era_impulse,
    "num_control_impulse_response",
    "crate::kernels::control::ident::impulse_response",
    "Computes the first-input impulse response in sample-by-output order from a fitted ERA model. The first row is zero and a zero-step result preserves the output column count.",
    true,
    (model: Handle<EraModel>, steps: usize) -> Mat,
    |ctx, model, steps| model.with(ctx, |model| {
        if steps > MAX_IMPULSE_STEPS {
            return Err(error("ERA impulse sample count exceeds the work limit"));
        }
        Ok(row_matrix(impulse_response(model, steps).map_err(error)?, model.c.len()))
    })
);

object_binding!(
    era_matrices,
    "num_control_era_matrices",
    "",
    "Returns the fitted A, B and C matrices and realisation order from a live ERA handle.",
    true,
    (model: Handle<EraModel>) -> EraMatrices,
    |ctx, model| model.with(ctx, |model| Ok(EraMatrices {
        a: row_matrix(model.a.clone(), model.order), b: row_matrix(model.b.clone(), 1),
        c: row_matrix(model.c.clone(), model.order), order: model.order,
    }))
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps control contracts and registrations with their native adapters"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[
            &pid_closed_loop::BINDING,
            &pid_sampled::BINDING,
            &era_fit::BINDING,
            &era_impulse::BINDING,
            &era_matrices::BINDING,
        ],
        "numerical_vector",
        DimRule::Dimensionless,
    )
}
