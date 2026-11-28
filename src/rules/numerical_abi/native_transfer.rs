//! Checked, reusable SISO transfer functions for the native value ABI.
use super::codecs::Output;
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernels::control::transfer::TransferFunction;
use crate::kernels::control::transfer::feedback;
use crate::kernels::control::transfer::gain_margin_phase_margin;
use crate::kernels::control::transfer::parallel;
use crate::kernels::control::transfer::series;
use crate::kernels::control::transfer::stable_poles;
use crate::kernels::control::transfer::step_response;
use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::object_binding;
use num_complex::Complex64;
use std::cell::OnceCell;
use std::mem::size_of;

const MAX_COEFFICIENTS: usize = 129;
const MAX_STEP_SAMPLES: usize = 16_384;
const MAX_STEP_WORK: usize = 1 << 24;

type Roots = Result<Vec<Complex64>, String>;

pub(super) struct TransferPlan {
    pub(super) model: TransferFunction,
    poles: OnceCell<Roots>,
    zeros: OnceCell<Roots>,
}

impl NativeType for TransferPlan {
    const KIND: &'static str = "rssn/native/control/transfer/v1";
    // Reserve the maximum coefficient and root allocations, including caches
    // populated by later queries. Compositions can retain coefficient capacity
    // from before leading zeros were trimmed. Cached errors also fit.
    const MAX_BYTES: usize = size_of::<Self>()
        + 4 * MAX_COEFFICIENTS * size_of::<f64>()
        + 2 * MAX_COEFFICIENTS * size_of::<Complex64>()
        + 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

impl TransferPlan {
    fn prepare(model: TransferFunction) -> Result<NewObject<Self>, OpError> {
        check_coefficient_count(model.num.len(), model.den.len())?;
        Ok(NewObject::normalized(Self {
            model,
            poles: OnceCell::new(),
            zeros: OnceCell::new(),
        }))
    }
}

fn check_coefficient_count(
    numerator: usize,
    denominator: usize,
) -> Result<(), OpError> {
    if numerator.max(denominator) > MAX_COEFFICIENTS {
        return Err(error("transfer degree exceeds 128"));
    }
    Ok(())
}

struct TransferMargins {
    gain_db: Option<f64>,
    phase_degrees: Option<f64>,
}

impl Output for TransferMargins {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("gain_db".into(), Option::<f64>::ty()),
            ("phase_degrees".into(), Option::<f64>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("gain_db", self.gain_db.encode()?),
            ("phase_degrees", self.phase_degrees.encode()?),
        ]))
    }
}

macro_rules! constructor {
    ($name:ident, $op:literal, $source:literal, $doc:literal, $make:ident) => {
        object_binding!(
            $name, $op, $source, $doc, true,
            (num: Vec<f64>, den: Vec<f64>) -> NewObject<TransferPlan>,
            |_, num, den| {
                check_coefficient_count(num.len(), den.len())?;
                TransferPlan::prepare(TransferFunction::$make(num, den).map_err(error)?)
            }
        );
    };
}

constructor!(
    transfer_new,
    "num_control_transfer_new",
    "crate::kernels::control::transfer::TransferFunction::new",
    "Prepares a real SISO transfer function from finite coefficients in descending powers of normalized s. A nonzero denominator is required; improper functions are allowed.",
    new
);
constructor!(
    transfer_new_proper,
    "num_control_transfer_new_proper",
    "crate::kernels::control::transfer::TransferFunction::new_proper",
    "Prepares a proper real SISO transfer function from finite coefficients in descending powers of normalized s.",
    new_proper
);

macro_rules! inquiry {
    ($name:ident, $op:literal, $source:literal, $doc:literal, $method:ident -> $ret:ty) => {
        object_binding!(
            $name, $op, $source, $doc, true,
            (plan: Handle<TransferPlan>) -> $ret,
            |ctx, plan| plan.with(ctx, |plan| Ok(plan.model.$method()))
        );
    };
}
inquiry!(
    transfer_is_proper, "num_control_transfer_is_proper",
    "crate::kernels::control::transfer::TransferFunction::is_proper",
    "Returns whether the numerator degree is at most the denominator degree.",
    is_proper -> bool
);
inquiry!(
    transfer_is_strictly_proper, "num_control_transfer_is_strictly_proper",
    "crate::kernels::control::transfer::TransferFunction::is_strictly_proper",
    "Returns whether the numerator degree is below the denominator degree. The zero transfer is strictly proper.",
    is_strictly_proper -> bool
);
inquiry!(
    transfer_order, "num_control_transfer_order",
    "crate::kernels::control::transfer::TransferFunction::order",
    "Returns the denominator degree after trimming leading zeros.",
    order -> usize
);

object_binding!(
    transfer_is_stable,
    "num_control_transfer_is_stable",
    "crate::kernels::control::transfer::TransferFunction::is_stable",
    "Checks that all continuous-time poles have negative real parts. Returns an error if the poles cannot be computed.",
    true,
    (plan: Handle<TransferPlan>) -> bool,
    |ctx, plan| plan.with(ctx, |plan| plan.poles.get_or_init(|| plan.model.poles_checked())
        .as_ref().map(|poles| stable_poles(poles)).map_err(error))
);

object_binding!(
    transfer_zeros,
    "num_control_transfer_zeros",
    "crate::kernels::control::transfer::TransferFunction::zeros",
    "Returns the finite zeros of the transfer numerator as a complex vector. Constant and zero numerators have an empty result.",
    true,
    (plan: Handle<TransferPlan>) -> Vec<Complex64>,
    |ctx, plan| plan.with(ctx, |plan| plan.zeros.get_or_init(|| plan.model.zeros_checked())
        .as_ref().cloned().map_err(error))
);

object_binding!(
    transfer_evaluate,
    "num_control_transfer_evaluate",
    "crate::kernels::control::transfer::TransferFunction::evaluate",
    "Evaluates a transfer function at a finite normalized complex frequency. Rejects a pole or a response outside Float range.",
    true,
    (plan: Handle<TransferPlan>, s: Complex64) -> Complex64,
    |ctx, plan, s| plan.with(ctx, |plan| plan.model.evaluate_checked(s).map_err(error))
);

macro_rules! connect {
    ($name:ident, $op:literal, $source:literal, $doc:literal, $connect:ident) => {
        object_binding!(
            $name, $op, $source, $doc, true,
            (left: Handle<TransferPlan>, right: Handle<TransferPlan>) -> NewObject<TransferPlan>,
            |ctx, left, right| left.with(ctx, |left| right.with(ctx, |right| {
                TransferPlan::prepare($connect(&left.model, &right.model).map_err(error)?)
            }))
        );
    };
}
connect!(
    transfer_series,
    "num_control_transfer_series",
    "crate::kernels::control::transfer::series",
    "Prepares the series connection a b.",
    series
);
connect!(
    transfer_parallel,
    "num_control_transfer_parallel",
    "crate::kernels::control::transfer::parallel",
    "Prepares the parallel connection a + b.",
    parallel
);
connect!(
    transfer_feedback,
    "num_control_transfer_feedback",
    "crate::kernels::control::transfer::feedback",
    "Prepares the negative feedback connection a / (1 + a b).",
    feedback
);

object_binding!(
    transfer_margins,
    "num_control_transfer_gain_margin_phase_margin",
    "crate::kernels::control::transfer::gain_margin_phase_margin",
    "Returns gain margin in dB and phase margin in degrees. Each margin is None when its crossover is absent on the normalized log-frequency grid.",
    true,
    (plan: Handle<TransferPlan>) -> TransferMargins,
    |ctx, plan| plan.with(ctx, |plan| {
        let (gain_db, phase_degrees) = gain_margin_phase_margin(&plan.model).map_err(error)?;
        Ok(TransferMargins { gain_db, phase_degrees })
    })
);

object_binding!(
    transfer_step_response,
    "num_control_transfer_step_response",
    "crate::kernels::control::transfer::step_response",
    "Returns the fourth-order Runge-Kutta unit-step response of a proper model on a finite, strictly increasing normalized time grid. The initial state is zero at the first sample.",
    true,
    (plan: Handle<TransferPlan>, time: Vec<f64>) -> Vec<f64>,
    |ctx, plan, time| plan.with(ctx, |plan| {
        let order = plan.model.order();
        if time.len() > MAX_STEP_SAMPLES
            || time.len().saturating_mul(order.saturating_mul(order)) > MAX_STEP_WORK
        {
            return Err(error("transfer step response exceeds the work budget"));
        }
        step_response(&plan.model, &time).map_err(error)
    })
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps the transfer bindings with their graph dimension contract"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[
            &transfer_new::BINDING,
            &transfer_new_proper::BINDING,
            &transfer_is_proper::BINDING,
            &transfer_is_strictly_proper::BINDING,
            &transfer_order::BINDING,
            &transfer_is_stable::BINDING,
            &transfer_zeros::BINDING,
            &transfer_evaluate::BINDING,
            &transfer_series::BINDING,
            &transfer_parallel::BINDING,
            &transfer_feedback::BINDING,
            &transfer_margins::BINDING,
            &transfer_step_response::BINDING,
        ],
        "numerical_vector",
        DimRule::Dimensionless,
    )
}
