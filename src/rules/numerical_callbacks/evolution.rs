//! Dimensioned vector IVPs share one callback, shape and result contract.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::canonical;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::expression;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::ode_adaptive::OdeError;
use crate::kernels::ode_adaptive::OdeOptions;
use crate::kernels::ode_adaptive::OdeSolution;
use crate::kernels::ode_adaptive::{
    self,
};
use crate::kernels::ode_stiff;

#[derive(Clone, Copy)]
enum Method {
    Dopri5,
    Rosenbrock23,
    Radau5,
    Bdf,
}
pub(super) fn tensor_type(
    rank: u64,
    dimension: DimExpr,
) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(dimensional_type(dimension)),
    }
}
// N is an operator-side existential shared by the initial state and RHS input.
// Runtime preflight additionally checks that relationship against the concrete
// callback bank signature, before any host invocation.
pub(super) fn state_type(dimension: DimExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: Some(vec![NatExpr::Var("N".into())]),
        elem: Box::new(dimensional_type(dimension)),
    }
}

pub(super) fn vector_dimension(
    ty: &TypeExpr,
    length: Option<usize>,
) -> Result<DimExpr, OpError> {
    let TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape,
        elem,
    } = ty
    else {
        return Err(OpError::Unsupported {
            backend: "numerical callback".into(),
            feature: "rank-one real vector schema".into(),
        });
    };
    if let Some(shape) = shape {
        let [NatExpr::Lit(extent)] = shape.as_slice() else {
            return Err(OpError::Unsupported {
                backend: "numerical callback".into(),
                feature: "concrete vector extent".into(),
            });
        };
        if length.is_some_and(|length| u64::try_from(length).ok() != Some(*extent)) {
            return Err(OpError::Domain(
                "callback vector extent differs from the expected length".into(),
            ));
        }
    }
    scalar_dimension(elem)
}
pub(super) fn vector(
    value: &Value,
    dimension: &DimExpr,
    length: Option<usize>,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<Vec<f64>, OpError> {
    let Value::Tensor(tensor) = value else {
        return Err(OpError::Type {
            expected: Box::new(tensor_type(1, dimension.clone())),
            got: Box::new(ctx.value_type(value)?),
            arg: index,
        });
    };
    if tensor.shape.len() != 1
        || tensor.shape[0] != tensor.data.len()
        || length.is_some_and(|length| length != tensor.data.len())
    {
        return Err(OpError::Domain(format!(
            "argument {index} must be a vector of the expected length with consistent storage"
        )));
    }
    let actual = vector_dimension(&ctx.value_type(value)?, length)?;
    if canonical(&actual)? != canonical(dimension)? {
        return Err(OpError::Type {
            expected: Box::new(tensor_type(1, dimension.clone())),
            got: Box::new(ctx.value_type(value)?),
            arg: index,
        });
    }
    tensor
        .data
        .iter()
        .map(|value| {
            // Bare coherent coefficients use the declared tensor element schema;
            // an explicit per-element unit may never contradict that schema.
            if matches!(value, Value::Dyn(_) | Value::Quantity(_)) {
                check_dimension(ctx, value, dimension, index)?;
            }
            numeric(value, index)
        })
        .collect()
}
pub(super) fn output_tensor(
    shape: Vec<usize>,
    data: Vec<f64>,
    dimension: DimExpr,
) -> Result<Value, OpError> {
    if shape
        .iter()
        .try_fold(1usize, |size, extent| size.checked_mul(*extent))
        != Some(data.len())
        || data.iter().any(|value| !value.is_finite())
    {
        return Err(OpError::Domain(
            "numerical result has malformed storage or non-finite coefficients".into(),
        ));
    }
    Ok(Value::Tensor(TensorValue {
        shape,
        data: data.into_iter().map(Value::Float).collect(),
        elem: dimensional_type(dimension),
    }))
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let rate = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let mut signature = Signature::new(
            vec![
                dimensional_type(d.clone()),
                dimensional_type(d.clone()),
                state_type(e.clone()),
                dimensional_type(e.clone()),
                TypeExpr::Float,
                dimensional_type(d.clone()),
                TypeExpr::Int,
                TypeExpr::Function(
                    vec![dimensional_type(d.clone()), state_type(e.clone())],
                    Box::new(tensor_type(1, rate)),
                ),
            ],
            TypeExpr::Record(vec![
                ("times".into(), tensor_type(1, d)),
                ("states".into(), tensor_type(2, e)),
                ("evaluations".into(), TypeExpr::Int),
                ("rejected".into(), TypeExpr::Int),
            ]),
        );
        signature.type_params = ["D", "E"]
            .into_iter()
            .map(|name| {
                TypeParam {
                    name: name.into(),
                    kind: TypeParamKind::Dim,
                    bounds: vec![],
                }
            })
            .collect();
        signature.type_params.push(TypeParam {
            name: "N".into(),
            kind: TypeParamKind::Nat,
            bounds: vec![],
        });
        super::common::effectful(signature)
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 8 {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(8),
                got: args.len(),
            });
        }
        let (callback, signature) =
            CallbackRun::new(ctx, &args[7], 7, 2, self.signature().params[7].clone())?;
        let time_dimension = scalar_dimension(&signature.params[0])?;
        let state_dimension = vector_dimension(&signature.params[1], None)?;
        let state = vector(&args[2], &state_dimension, None, 2, ctx)?;
        if state.is_empty() {
            return Err(OpError::Domain("ODE initial state must be nonempty".into()));
        }
        vector_dimension(&signature.params[1], Some(state.len()))?;
        let rate_dimension = expression(&canonical(&DimExpr::Mul(vec![
            state_dimension.clone(),
            DimExpr::Inverse(Box::new(time_dimension.clone())),
        ]))?)?;
        let actual_rate = vector_dimension(&signature.ret, Some(state.len()))?;
        if canonical(&actual_rate)? != canonical(&rate_dimension)? {
            return Err(OpError::Domain(
                "ODE callback derivative dimension must be state divided by time".into(),
            ));
        }
        for index in [0usize, 1, 5] {
            check_dimension(ctx, &args[index], &time_dimension, index)?;
        }
        let t0 = numeric(&args[0], 0)?;
        let t1 = numeric(&args[1], 1)?;
        let atol = numeric(&args[3], 3)?;
        if !crate::kernels::exact_eq(atol, 0.0) {
            check_dimension(ctx, &args[3], &state_dimension, 3)?;
        }
        check_dimension(ctx, &args[4], &DimExpr::Dimless, 4)?;
        let rtol = numeric(&args[4], 4)?;
        let hmax = numeric(&args[5], 5)?;
        let max_steps = positive_budget(&args[6], 6, "ODE attempt budget")?;
        let options = OdeOptions {
            atol,
            rtol,
            hmax,
            max_steps,
            h0: None,
        };
        let derivative = |time: f64, state: &[f64], output: &mut [f64]| {
            let arguments = [
                dimensional_float(time, time_dimension.clone()),
                Value::Tensor(TensorValue {
                    shape: vec![state.len()],
                    data: state.iter().copied().map(Value::Float).collect(),
                    elem: dimensional_type(state_dimension.clone()),
                }),
            ];
            match callback.evaluate(&arguments, |value| {
                vector(&value, &rate_dimension, Some(state.len()), 7, ctx)
            }) {
                | Some(values) => output.copy_from_slice(&values),
                | None => output.fill(f64::NAN),
            }
        };
        let stop = || callback.failed();
        let solution = match self {
            | Self::Dopri5 => {
                ode_adaptive::dopri5_controlled(
                    derivative,
                    t0,
                    t1,
                    &state,
                    &options,
                    &[],
                    Some(&stop),
                )
            },
            | Self::Rosenbrock23 => {
                ode_adaptive::rosenbrock23_events_controlled(
                    derivative,
                    t0,
                    t1,
                    &state,
                    &options,
                    &[],
                    Some(&stop),
                )
            },
            | Self::Radau5 => {
                ode_stiff::radau5_controlled(derivative, t0, t1, &state, &options, &stop)
            },
            | Self::Bdf => ode_stiff::bdf_controlled(derivative, t0, t1, &state, &options, &stop),
        };
        let evaluations = callback.evaluations();
        callback.finish()?;
        let OdeSolution { t, y, rejected, .. } = solution.map_err(|error| {
            match error {
                | OdeError::Invalid => {
                    OpError::Domain(
                        "invalid ODE times, state, tolerances, maximum step or attempt budget"
                            .into(),
                    )
                },
                | OdeError::NonFinite => {
                    OpError::Domain("ODE state or derivative became non-finite".into())
                },
                | OdeError::StepLimit => {
                    OpError::NotReduced(
                        "ODE attempt budget exhausted or time step is not representable".into(),
                    )
                },
                | OdeError::SolveFailed => {
                    OpError::NotReduced("ODE linear/nonlinear solve failed".into())
                },
            }
        })?;
        if y.len() != t.len() || y.iter().any(|row| row.len() != state.len()) {
            return Err(OpError::Domain(
                "ODE solver returned inconsistent shape or evaluation counters".into(),
            ));
        }
        let samples = t.len();
        Ok(Value::Record(RecordValue {
            fields: vec![
                (
                    "times".into(),
                    output_tensor(vec![samples], t, time_dimension)?,
                ),
                (
                    "states".into(),
                    output_tensor(
                        vec![samples, state.len()],
                        y.into_iter().flatten().collect(),
                        state_dimension,
                    )?,
                ),
                ("evaluations".into(), Value::int(evaluations)),
                ("rejected".into(), Value::int(rejected)),
            ],
        }))
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        (
            "num_ode_dopri5",
            "crate::kernels::ode_adaptive::dopri5",
            Method::Dopri5,
        ),
        (
            "num_ode_rosenbrock23",
            "crate::kernels::ode_adaptive::rosenbrock23",
            Method::Rosenbrock23,
        ),
        (
            "num_ode_radau5",
            "crate::kernels::ode_stiff::radau5",
            Method::Radau5,
        ),
        ("num_ode_bdf", "crate::kernels::ode_stiff::bdf", Method::Bdf),
    ] {
        register(
            installer,
            name,
            path,
            "Serial adaptive vector IVP: initial time, final time, initial state vector, absolute state tolerance, relative tolerance, positive finite maximal time step, positive attempt budget, RHS callback. RHS(time,state) returns a same-length vector with state/time dimensions. Returns accepted times/states/evaluations/rejected; the rank-two state tensor has one row per accepted time. Coefficients use coherent units. Callback failures propagate unchanged, rejected attempts consume the budget, and numerical error estimates are heuristic. Events, DAEs and dense-output handles use the native Rust API.",
            8,
            method,
        )?;
    }
    Ok(())
}
