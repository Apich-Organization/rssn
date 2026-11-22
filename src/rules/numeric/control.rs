//! Control systems and system identification.
//!
//! | operator | value |
//! |---|---|
//! | `care(a, b, q, r)` | the stabilising solution `P` of the continuous algebraic Riccati equation |
//! | `dare(a, b, q, r)` | the stabilising solution `P` of the discrete algebraic Riccati equation |
//! | `lqr(a, b, q, r)` | `list(K, P)`: the continuous optimal gain and the Riccati solution |
//! | `dlqr(a, b, q, r)` | `list(K, P)`: the discrete optimal gain and the Riccati solution |
//! | `lqg(a, b, c, q, r, w, v)` | `list(K, L)`: the LQR gain and the steady-state Kalman gain |
//! | `kalman_gain(a, c, w, v)` | the steady-state Kalman gain of `x+ = a x + w` observed through `c` with covariances `w` and `v` |
//! | `tf_state_space(num, den)` | `list(a, b, c, d)`: the controllable canonical realisation of `num(s)/den(s)` |
//! | `ss_transfer(a, b, c, d)` | `list(num, den)`: the transfer function from the first input to the first output |
//! | `tf_poles(num, den)` | the poles as `list(re, im)` pairs |
//! | `tf_bode(num, den, omega)` | `list(magnitude_db, phase_deg)` at the frequencies `omega` |
//! | `tf_nyquist(num, den, omega)` | `list(re, im)` of `G(i omega)` |
//! | `pid_tune(num, den, rule, lambda)` | `list(kp, ki, kd)` from a first-order-plus-dead-time model of the plant |
//! | `mpc_plan(a, b, c, d, x0, reference, horizon, q, r)` | `list(u, predicted_states, objective, iterations)`: the condensed-QP receding-horizon plan from `x0` |
//! | `arx_fit(u, y, na, nb, delay)` | `list(a, b, fit)`: the ARX coefficients and the fit quality |
//!
//! `num` and `den` are coefficient lists in **descending** powers of `s`;
//! `a`, `b`, `c`, `d` are the state-space matrices as lists of rows. `rule` is
//! `ziegler_nichols`, `cohen_coon`, `imc`, `amigo` or `chien_hrones_reswick`,
//! and `lambda` is the closed-loop time constant for `imc` (use `0` to let the
//! rule choose).
//!
//! A model the method cannot use — a singular `r`, a plant with no
//! first-order-plus-dead-time approximation, an uncontrollable pair where the
//! Riccati equation has no stabilising solution — is a [`DiagKind::Domain`]
//! finding, not a matrix of numbers that is not a gain.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::require;
use super::conic::matrix;
use super::conic::numbers;
use crate::kernel::diag::DiagKind;
use crate::kernels::control::StateSpace;
use crate::kernels::control::TransferFunction;
use crate::kernels::control::TuningRule;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "care",
        params: "tttt",
        eval: care,
        doc: "`care(a, b, q, r)`: the stabilising solution `P` of the continuous algebraic Riccati equation `A'P + PA - PBR^-1B'P + Q = 0`.",
    },
    Entry {
        name: "dare",
        params: "tttt",
        eval: dare,
        doc: "`dare(a, b, q, r)`: the stabilising solution `P` of the discrete algebraic Riccati equation.",
    },
    Entry {
        name: "lqr",
        params: "tttt",
        eval: lqr,
        doc: "`lqr(a, b, q, r)`: `list(K, P)`, the continuous linear-quadratic gain `K = R^-1 B' P` and the Riccati solution.",
    },
    Entry {
        name: "dlqr",
        params: "tttt",
        eval: dlqr,
        doc: "`dlqr(a, b, q, r)`: `list(K, P)`, the discrete linear-quadratic gain and the Riccati solution.",
    },
    Entry {
        name: "lqg",
        params: "ttttttt",
        eval: lqg,
        doc: "`lqg(a, b, c, q, r, w, v)`: `list(K, L)`, the LQR gain and the steady-state Kalman gain of the dual problem.",
    },
    Entry {
        name: "kalman_gain",
        params: "tttt",
        eval: kalman_gain,
        doc: "`kalman_gain(a, c, w, v)`: the steady-state Kalman gain of `x+ = a x + w` observed through `c`, with process covariance `w` and noise covariance `v`.",
    },
    Entry {
        name: "tf_state_space",
        params: "tt",
        eval: tf_state_space,
        doc: "`tf_state_space(num, den)`: `list(a, b, c, d)`, the controllable canonical realisation of the transfer function `num(s)/den(s)` (descending powers).",
    },
    Entry {
        name: "ss_transfer",
        params: "tttt",
        eval: ss_transfer,
        doc: "`ss_transfer(a, b, c, d)`: `list(num, den)`, the transfer function from the first input to the first output.",
    },
    Entry {
        name: "tf_poles",
        params: "tt",
        eval: tf_poles,
        doc: "`tf_poles(num, den)`: the poles of the transfer function as `list(re, im)` pairs.",
    },
    Entry {
        name: "tf_bode",
        params: "ttv",
        eval: tf_bode,
        doc: "`tf_bode(num, den, omega)`: `list(magnitude_db, phase_deg)` of the frequency response at `omega`.",
    },
    Entry {
        name: "tf_nyquist",
        params: "ttv",
        eval: tf_nyquist,
        doc: "`tf_nyquist(num, den, omega)`: `list(re, im)` of `G(i omega)` at `omega`.",
    },
    Entry {
        name: "pid_tune",
        params: "ttnd",
        eval: pid_tune,
        doc: "`pid_tune(num, den, rule, lambda)`: `list(kp, ki, kd)` from a first-order-plus-dead-time model of the plant; `rule` is `ziegler_nichols`, `cohen_coon`, `imc`, `amigo` or `chien_hrones_reswick`.",
    },
    Entry {
        name: "mpc_plan",
        params: "ttttttddd",
        eval: mpc_plan,
        doc: "`mpc_plan(a, b, c, d, x0, reference, horizon, q, r)`: `list(u, predicted_states, objective, iterations)` of one receding-horizon step solved as a condensed quadratic program.",
    },
    Entry {
        name: "arx_fit",
        params: "vvddd",
        eval: arx_fit,
        doc: "`arx_fit(u, y, na, nb, delay)`: `list(a, b, fit)` of an ARX model by least squares.",
    },
];

/// A matrix result.
fn matrix_out(rows: Vec<Vec<f64>>) -> Data {
    Data::Tree(rows.into_iter().map(Data::vector).collect())
}

/// A transfer function from its two coefficient lists.
fn transfer(args: &[Data]) -> Option<TransferFunction> {
    let numerator = numbers(args, 0)?;
    let denominator = numbers(args, 1)?;
    match TransferFunction::new(numerator, denominator) {
        | Ok(model) => Some(model),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

/// A state-space model from four matrix arguments.
fn state_space(args: &[Data]) -> Option<StateSpace> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let states = a.len();
    let Some(b) = matrix(args, 1, None) else {
        return None;
    };
    let Some(c) = matrix(args, 2, None) else {
        return None;
    };
    let Some(d) = matrix(args, 3, None) else {
        return None;
    };
    if a.iter().any(|row| row.len() != states) {
        return bad(0, DiagKind::Shape, format_args!("`a` must be square"));
    }
    match StateSpace::new(a, b, c, d) {
        | Ok(model) => Some(model),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn care(args: &[Data]) -> Option<Data> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let size = a.len();
    let b = matrix(args, 1, None)?;
    let q = matrix(args, 2, Some(size))?;
    let r = matrix(args, 3, None)?;
    match crate::kernels::control::care(&a, &b, &q, &r) {
        | Ok(p) => Some(matrix_out(p)),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn dare(args: &[Data]) -> Option<Data> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let size = a.len();
    let b = matrix(args, 1, None)?;
    let q = matrix(args, 2, Some(size))?;
    let r = matrix(args, 3, None)?;
    match crate::kernels::control::dare(&a, &b, &q, &r) {
        | Ok(p) => Some(matrix_out(p)),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

/// The common shape of `lqr` and `dlqr`.
fn quadratic_gain(
    args: &[Data],
    discrete: bool,
) -> Option<Data> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let size = a.len();
    let b = matrix(args, 1, None)?;
    let q = matrix(args, 2, Some(size))?;
    let r = matrix(args, 3, None)?;
    let solved = if discrete {
        crate::kernels::control::dlqr(&a, &b, &q, &r)
    } else {
        crate::kernels::control::lqr(&a, &b, &q, &r)
    };
    match solved {
        | Ok((gain, p)) => Some(Data::Tree(vec![matrix_out(gain), matrix_out(p)])),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn lqr(args: &[Data]) -> Option<Data> {
    quadratic_gain(args, false)
}

fn dlqr(args: &[Data]) -> Option<Data> {
    quadratic_gain(args, true)
}

fn lqg(args: &[Data]) -> Option<Data> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let size = a.len();
    let b = matrix(args, 1, None)?;
    let c = matrix(args, 2, None)?;
    let q = matrix(args, 3, Some(size))?;
    let r = matrix(args, 4, None)?;
    let w = matrix(args, 5, Some(size))?;
    let v = matrix(args, 6, None)?;
    match crate::kernels::control::lqg(&a, &b, &c, &q, &r, &w, &v) {
        | Ok((gain, estimator)) => Some(Data::Tree(vec![matrix_out(gain), matrix_out(estimator)])),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn kalman_gain(args: &[Data]) -> Option<Data> {
    let Some(a) = args.get(0).and_then(super::bridge::rows) else {
        return bad(0, DiagKind::Shape, format_args!("expected the matrix `a`"));
    };
    let size = a.len();
    let c = matrix(args, 1, None)?;
    let w = matrix(args, 2, Some(size))?;
    let v = matrix(args, 3, None)?;
    match crate::kernels::control::kalman_gain(&a, &c, &w, &v) {
        | Ok(gain) => Some(matrix_out(gain)),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn tf_state_space(args: &[Data]) -> Option<Data> {
    let model = transfer(args)?;
    match model.to_state_space() {
        | Ok(space) => {
            Some(Data::Tree(vec![
                matrix_out(space.a),
                matrix_out(space.b),
                matrix_out(space.c),
                matrix_out(space.d),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn ss_transfer(args: &[Data]) -> Option<Data> {
    let space = state_space(args)?;
    match crate::kernels::control::ss_to_tf(&space) {
        | Ok(models) => {
            let Some(first) = models.first() else {
                return bad(
                    0,
                    DiagKind::Domain,
                    format_args!("the model has no input/output pair"),
                );
            };
            Some(Data::Tree(vec![
                Data::vector(first.num.clone()),
                Data::vector(first.den.clone()),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn tf_poles(args: &[Data]) -> Option<Data> {
    let model = transfer(args)?;
    let poles = model
        .poles()
        .into_iter()
        .map(|pole| Data::vector(vec![pole.re, pole.im]))
        .collect::<Vec<Data>>();
    Some(Data::Tree(poles))
}

fn tf_bode(args: &[Data]) -> Option<Data> {
    let model = transfer(args)?;
    let omega = numbers(args, 2)?;
    match crate::kernels::control::bode(&model, &omega) {
        | Ok((magnitude, phase)) => {
            Some(Data::Tree(vec![
                Data::vector(magnitude),
                Data::vector(phase),
            ]))
        },
        | Err(message) => bad(2, DiagKind::Domain, format_args!("{message}")),
    }
}

fn tf_nyquist(args: &[Data]) -> Option<Data> {
    let model = transfer(args)?;
    let omega = numbers(args, 2)?;
    match crate::kernels::control::nyquist(&model, &omega) {
        | Ok(points) => {
            let real = points.iter().map(|&(re, _)| re).collect();
            let imaginary = points.iter().map(|&(_, im)| im).collect();
            Some(Data::Tree(vec![
                Data::vector(real),
                Data::vector(imaginary),
            ]))
        },
        | Err(message) => bad(2, DiagKind::Domain, format_args!("{message}")),
    }
}

fn pid_tune(args: &[Data]) -> Option<Data> {
    let model = transfer(args)?;
    let rule_text = super::bridge::name(args, 2)?;
    let rule = match rule_text {
        | "ziegler_nichols" | "zn" => TuningRule::ZieglerNichols,
        | "cohen_coon" | "cc" => TuningRule::CohenCoon,
        | "imc" => TuningRule::Imc,
        | "amigo" => TuningRule::Amigo,
        | "chien_hrones_reswick" | "chr" => TuningRule::ChienHronesReswick,
        | _ => {
            return bad(
                2,
                DiagKind::Domain,
                format_args!(
                    "`{rule_text}` is not a tuning rule (ziegler_nichols, cohen_coon, imc, \
                     amigo, chien_hrones_reswick)"
                ),
            );
        },
    };
    let lambda = super::bridge::scalar(args, 3)?;
    let chosen = if lambda > 0.0 {
        Some(lambda)
    } else {
        None
    };
    match crate::kernels::control::tune(&model, rule, chosen) {
        | Ok(pid) => Some(Data::vector(vec![pid.kp, pid.ki, pid.kd])),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn mpc_plan(args: &[Data]) -> Option<Data> {
    let space = state_space(args)?;
    let start = numbers(args, 4)?;
    let reference = numbers(args, 5)?;
    let horizon = count(args, 6, 512)?;
    let state_weights = numbers(args, 7)?;
    let input_weights = numbers(args, 8)?;
    require(
        start.len() == space.n_states(),
        Some(4),
        DiagKind::Shape,
        format_args!(
            "the initial state needs {} entries, got {}",
            space.n_states(),
            start.len()
        ),
    )?;
    require(
        state_weights.len() == space.n_states(),
        Some(7),
        DiagKind::Shape,
        format_args!(
            "the state weights need {} entries, got {}",
            space.n_states(),
            state_weights.len()
        ),
    )?;
    require(
        input_weights.len() == space.n_inputs(),
        Some(8),
        DiagKind::Shape,
        format_args!(
            "the input weights need {} entries, got {}",
            space.n_inputs(),
            input_weights.len()
        ),
    )?;
    require(
        horizon > 0,
        Some(6),
        DiagKind::Domain,
        format_args!("the horizon must be at least 1"),
    )?;
    let options = crate::kernels::control::MpcOptions {
        horizon,
        control_horizon: horizon,
        q: state_weights,
        r: input_weights,
        u_min: None,
        u_max: None,
        x_min: None,
        x_max: None,
    };
    match crate::kernels::control::mpc_step(&space, &start, &reference, &options) {
        | Ok(step) => {
            #[allow(clippy::cast_precision_loss)]
            let iterations = step.iterations as f64;
            Some(Data::Tree(vec![
                Data::vector(step.u),
                Data::Tree(
                    step.predicted_states
                        .into_iter()
                        .map(Data::vector)
                        .collect(),
                ),
                Data::Scalar(step.objective),
                Data::Scalar(iterations),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn arx_fit(args: &[Data]) -> Option<Data> {
    let input = super::bridge::vector(args, 0)?;
    let output = super::bridge::vector(args, 1)?;
    let na = count(args, 2, 32)?;
    let nb = count(args, 3, 32)?;
    let delay = count(args, 4, 32)?;
    match crate::kernels::control::arx(input, output, na, nb, delay) {
        | Ok(model) => {
            Some(Data::Tree(vec![
                Data::vector(model.a),
                Data::vector(model.b),
                Data::Scalar(model.fit),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
