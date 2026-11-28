//! Session-owned normalized state-space models.
use std::mem::size_of;
use std::sync::OnceLock;

use num_complex::Complex64;

use super::codecs::StateSpaceShape;
use super::codecs::checked_state_matrix_payload;
use super::codecs::state_space_error;
use super::codecs::state_space_result;
use super::codecs::state_space_rows;
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernels::control::StateSpace;
use crate::kernels::dense::Mat;
use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::object_binding;

struct PreparedStateSpace {
    model: StateSpace,
    shape: StateSpaceShape,
    poles: OnceLock<Vec<Complex64>>,
    stable: OnceLock<bool>,
    controllable: OnceLock<bool>,
    observable: OnceLock<bool>,
}

impl NativeType for PreparedStateSpace {
    const KIND: &'static str = "rssn/native/control/state-space/v1";
    const MAX_BYTES: usize = 64 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        let Self { model, shape, .. } = self;
        let mut bytes =
            size_of::<Self>().saturating_add(shape.states.saturating_mul(size_of::<Complex64>()));
        for rows in [&model.a, &model.b, &model.c, &model.d] {
            bytes = bytes.saturating_add(rows.capacity().saturating_mul(size_of::<Vec<f64>>()));
            for row in rows {
                bytes = bytes.saturating_add(row.capacity().saturating_mul(size_of::<f64>()));
            }
        }
        bytes
    }
}

const fn new_state(
    model: StateSpace,
    shape: StateSpaceShape,
) -> NewObject<PreparedStateSpace> {
    NewObject::normalized(PreparedStateSpace {
        model,
        shape,
        poles: OnceLock::new(),
        stable: OnceLock::new(),
        controllable: OnceLock::new(),
        observable: OnceLock::new(),
    })
}

macro_rules! continuous_model {
    ($id:ident, $name:literal, $ctor:ident, $path:literal, $doc:literal) => {
        object_binding!(
            $id, $name, $path, $doc, true,
            (a: Mat, b: Mat, c: Mat, d: Mat) -> NewObject<PreparedStateSpace>,
            |_, a, b, c, d| {
                let ([a, b, c, d], shape) = checked_state_matrix_payload(a, b, c, d)?;
                Ok(new_state(StateSpace::$ctor(a, b, c, d).map_err(state_space_error)?, shape))
            }
        );
    };
}

continuous_model!(
    state_default,
    "num_control_state_space_new",
    new,
    "crate::kernels::control::state_space::StateSpace::new",
    "Prepares a continuous model from normalized A, B, C and D matrices. Simulation times use the model's chosen time unit."
);

continuous_model!(
    state_continuous,
    "num_control_state_space_continuous",
    continuous,
    "crate::kernels::control::state_space::StateSpace::continuous",
    "Prepares a continuous state-space model from normalized A, B, C and D matrices."
);

object_binding!(
    state_discrete, "num_control_state_space_discrete",
    "crate::kernels::control::state_space::StateSpace::discrete",
    "Prepares a discrete model from normalized A, B, C and D matrices and a positive sampling interval.",
    true,
    (a: Mat, b: Mat, c: Mat, d: Mat, dt: f64) -> NewObject<PreparedStateSpace>,
    |_, a, b, c, d, dt| {
        let ([a, b, c, d], shape) = checked_state_matrix_payload(a, b, c, d)?;
        Ok(new_state(StateSpace::discrete(a, b, c, d, dt).map_err(state_space_error)?, shape))
    }
);

object_binding!(
    state_poles, "num_control_state_space_poles",
    "crate::kernels::control::state_space::StateSpace::poles",
    "Returns finite poles of a prepared model. A static model has no poles. Eigenvalue failure returns an error.",
    true,
    (model: Handle<PreparedStateSpace>) -> Vec<Complex64>,
    |ctx, model| model.with(ctx, |m| {
        m.shape.check_diagnostic_work(1)?;
        if let Some(poles) = m.poles.get() {
            return Ok(poles.clone());
        }
        crate::cancellation::checkpoint()?;
        let poles = m.model.poles();
        if poles.len() != m.shape.states {
            return Err(state_space_error("a finite pole is required for each state".into()));
        }
        crate::cancellation::checkpoint()?;
        let _ = m.poles.set(poles.clone());
        Ok(poles)
    })
);

macro_rules! state_property {
    ($id:ident, $name:literal, $field:ident, $method:ident, $path:literal, $doc:literal, $($check:tt)*) => {
        object_binding!(
            $id, $name, $path, $doc, true,
            (model: Handle<PreparedStateSpace>) -> bool,
            |ctx, model| model.with(ctx, |m| {
                let check = $($check)*;
                check(&m.shape)?;
                if let Some(result) = m.$field.get() {
                    return Ok(*result);
                }
                crate::cancellation::checkpoint()?;
                let result = m.model.$method();
                crate::cancellation::checkpoint()?;
                let _ = m.$field.set(result);
                Ok(result)
            })
        );
    };
}

state_property!(
    state_stable,
    "num_control_state_space_is_stable",
    stable,
    is_stable,
    "crate::kernels::control::state_space::StateSpace::is_stable",
    "Checks whether all continuous poles have negative real part or all discrete poles lie inside the unit circle.",
    |shape: &StateSpaceShape| shape.check_diagnostic_work(1)
);
state_property!(
    state_controllable,
    "num_control_state_space_is_controllable",
    controllable,
    is_controllable,
    "crate::kernels::control::state_space::StateSpace::is_controllable",
    "Checks whether the normalized state-space controllability matrix has full row rank.",
    |shape: &StateSpaceShape| shape.check_diagnostic_matrix(shape.inputs)
);
state_property!(
    state_observable,
    "num_control_state_space_is_observable",
    observable,
    is_observable,
    "crate::kernels::control::state_space::StateSpace::is_observable",
    "Checks whether the normalized state-space observability matrix has full column rank.",
    |shape: &StateSpaceShape| shape.check_diagnostic_matrix(shape.outputs)
);

object_binding!(
    state_simulate, "num_control_state_space_simulate_with_dt",
    "crate::kernels::control::state_space::StateSpace::simulate_with_dt",
    "Simulates a normalized model at positive sampling interval dt. A discrete model requires its configured sample interval.",
    true,
    (model: Handle<PreparedStateSpace>, x0: Vec<f64>, input: Mat, dt: f64) -> Mat,
    |ctx, model, x0, input, dt| model.with(ctx, |m| {
        if input.cols != m.shape.inputs {
            return Err(error("the input matrix width must match the prepared model"));
        }
        if input.rows.saturating_mul(m.shape.outputs.saturating_add(1)) >= api::VALUE_LIMIT {
            return Err(error("state-space trajectory exceeds the native tensor size limit"));
        }
        if m.model.dt.is_none() {
            m.shape.check_zoh_workspace()?;
        }
        let rows = m.model.simulate_with_dt(&x0, &state_space_rows(input), dt).map_err(state_space_error)?;
        crate::cancellation::checkpoint()?;
        Ok(state_space_result(rows, m.shape.outputs))
    })
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[
            &state_default::BINDING,
            &state_continuous::BINDING,
            &state_discrete::BINDING,
            &state_poles::BINDING,
            &state_stable::BINDING,
            &state_controllable::BINDING,
            &state_observable::BINDING,
            &state_simulate::BINDING,
        ],
        "numerical_matrix",
        DimRule::Dimensionless,
    )
}
