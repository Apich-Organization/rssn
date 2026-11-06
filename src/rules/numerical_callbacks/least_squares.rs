//! Dimensioned vector nonlinear least squares through a scoped residual callback.
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
use super::check_dimension;
use super::common::CallbackRun;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::state_type;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::numeric;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::rootfind::LeastSquaresOptions;
use crate::kernels::rootfind::LeastSquaresStatus;
use crate::kernels::rootfind::{
    self,
};

struct LeastSquares;
impl EvalDyn for LeastSquares {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let mut signature = Signature::new(
            vec![
                state_type(d.clone()),
                dimensional_type(e.clone()),
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
                TypeExpr::Function(
                    vec![state_type(d.clone())],
                    Box::new(tensor_type(1, e.clone())),
                ),
            ],
            TypeExpr::Record(vec![
                ("x".into(), state_type(d)),
                ("residual_norm".into(), dimensional_type(e)),
                ("iterations".into(), TypeExpr::Int),
                ("evaluations".into(), TypeExpr::Int),
                ("jacobian_evaluations".into(), TypeExpr::Int),
                ("accepted_steps".into(), TypeExpr::Int),
                ("rejected_steps".into(), TypeExpr::Int),
                ("converged".into(), TypeExpr::Bool),
                ("status".into(), TypeExpr::String),
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
            CallbackRun::new(ctx, &args[7], 7, 1, self.signature().params[7].clone())?;
        let domain = vector_dimension(&signature.params[0], None)?;
        let residual_dim = vector_dimension(&signature.ret, None)?;
        let initial = vector(&args[0], &domain, None, 0, ctx)?;
        if initial.is_empty() {
            return Err(OpError::Domain(
                "least-squares initial vector must be nonempty".into(),
            ));
        }
        vector_dimension(&signature.params[0], Some(initial.len()))?;
        check_dimension(ctx, &args[1], &residual_dim, 1)?;
        for index in 2..=4 {
            check_dimension(ctx, &args[index], &DimExpr::Dimless, index)?;
        }
        let options = LeastSquaresOptions {
            residual_tolerance: numeric(&args[1], 1)?,
            gradient_tolerance: numeric(&args[2], 2)?,
            step_tolerance: numeric(&args[3], 3)?,
            initial_damping: numeric(&args[4], 4)?,
            max_iterations: positive_budget(&args[5], 5, "least-squares iteration budget")?,
            max_evaluations: positive_budget(&args[6], 6, "least-squares evaluation budget")?,
        };
        let fixed_output = match &signature.ret {
            | TypeExpr::Tensor {
                shape: Some(shape), ..
            } => {
                match shape.as_slice() {
                    | [NatExpr::Lit(length)] => {
                        Some(usize::try_from(*length).map_err(|_| {
                            OpError::Domain("residual extent does not fit usize".into())
                        })?)
                    },
                    | _ => return Err(OpError::Domain("residual extent must be concrete".into())),
                }
            },
            | _ => None,
        };
        let result = rootfind::least_squares(
            |x| {
                let arg = Value::Tensor(TensorValue {
                    shape: vec![x.len()],
                    data: x.iter().copied().map(Value::Float).collect(),
                    elem: dimensional_type(domain.clone()),
                });
                callback
                    .evaluate(&[arg], |value| {
                        vector(&value, &residual_dim, fixed_output, 7, ctx)
                    })
                    .unwrap_or_default()
            },
            None,
            &initial,
            options,
        );
        callback.finish()?;
        let result = result
            .map_err(|error| OpError::Domain(format!("least-squares solver failed: {error:?}")))?;
        let status = match result.status {
            | LeastSquaresStatus::ResidualConverged => "residual_converged",
            | LeastSquaresStatus::GradientConverged => "gradient_converged",
            | LeastSquaresStatus::IterationLimit => "iteration_limit",
            | LeastSquaresStatus::EvaluationLimit => "evaluation_limit",
            | LeastSquaresStatus::Stagnated => "stagnated",
        };
        Ok(Value::Record(RecordValue {
            fields: vec![
                (
                    "x".into(),
                    output_tensor(vec![result.x.len()], result.x, domain)?,
                ),
                (
                    "residual_norm".into(),
                    dimensional_float(result.residual_norm, residual_dim),
                ),
                ("iterations".into(), Value::int(result.iterations)),
                ("evaluations".into(), Value::int(result.evaluations)),
                (
                    "jacobian_evaluations".into(),
                    Value::int(result.jacobian_evaluations),
                ),
                ("accepted_steps".into(), Value::int(result.accepted_steps)),
                ("rejected_steps".into(), Value::int(result.rejected_steps)),
                ("converged".into(), Value::Bool(result.status.converged())),
                ("status".into(), Value::String(status.into())),
            ],
        }))
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        installer,
        "num_least_squares",
        "crate::kernels::rootfind::least_squares",
        "Damped nonlinear least squares with a scoped vector residual callback. Initial parameters and residuals may have independent dimensions. Returns the best point, residual norm, work counters, convergence and explicit stopping status; callback failures propagate.",
        8,
        LeastSquares,
    )
}
