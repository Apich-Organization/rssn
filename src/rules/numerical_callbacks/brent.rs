//! Bracketed root finding with independently dimensioned domain and residual.
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
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::rootfind;
use crate::kernels::rootfind::RootError;

struct Brent {
    with_bracket: bool,
}
impl EvalDyn for Brent {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let mut result = vec![
            ("root".into(), dimensional_type(d.clone())),
            ("evaluations".into(), TypeExpr::Int),
        ];
        if self.with_bracket {
            result.push((
                "bracket".into(),
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(1),
                    shape: Some(vec![NatExpr::Lit(2)]),
                    elem: Box::new(dimensional_type(d.clone())),
                },
            ));
        }
        let mut signature = Signature::new(
            vec![
                dimensional_type(d.clone()),
                dimensional_type(d.clone()),
                dimensional_type(d.clone()),
                TypeExpr::Int,
                TypeExpr::Function(vec![dimensional_type(d)], Box::new(dimensional_type(e))),
            ],
            TypeExpr::Record(result),
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
        super::common::effectful(signature)
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 5 {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(5),
                got: args.len(),
            });
        }
        let (callback, signature) =
            CallbackRun::new(ctx, &args[4], 4, 1, self.signature().params[4].clone())?;
        let domain = scalar_dimension(&signature.params[0])?;
        scalar_dimension(&signature.ret)?;
        for (index, value) in args[..3].iter().enumerate() {
            check_dimension(ctx, value, &domain, index)?;
        }
        let a = numeric(&args[0], 0)?;
        let b = numeric(&args[1], 1)?;
        let tolerance = numeric(&args[2], 2)?;
        if tolerance <= 0.0 {
            return Err(OpError::Domain(
                "root time/argument tolerance must be positive".into(),
            ));
        }
        let iterations = positive_budget(&args[3], 3, "root iteration budget")?;
        let evaluate = |x| {
            callback
                .evaluate(&[dimensional_float(x, domain.clone())], |value| {
                    numeric(&value, 4)
                })
                .unwrap_or(f64::NAN)
        };
        let result = if self.with_bracket {
            rootfind::brent_with_bracket(evaluate, a, b, tolerance, iterations)
                .map(|result| (result.root, Some(result.bracket)))
        } else {
            rootfind::brent(evaluate, a, b, tolerance, iterations).map(|root| (root, None))
        };
        let evaluations = callback.evaluations();
        callback.finish()?;
        let (root, bracket) = result.map_err(|error| {
            match error {
                | RootError::NoConvergence => {
                    OpError::NotReduced("Brent root iteration budget exhausted".into())
                },
                | RootError::NotBracketed => {
                    OpError::Domain("root interval does not bracket a sign change".into())
                },
                | other => OpError::Domain(format!("Brent root solver failed: {other:?}")),
            }
        })?;
        let mut fields = vec![
            ("root".into(), dimensional_float(root, domain.clone())),
            ("evaluations".into(), Value::int(evaluations)),
        ];
        if let Some(bracket) = bracket {
            fields.push((
                "bracket".into(),
                Value::Tensor(TensorValue {
                    shape: vec![2],
                    data: bracket
                        .into_iter()
                        .map(|endpoint| dimensional_float(endpoint, domain.clone()))
                        .collect(),
                    elem: dimensional_type(domain),
                }),
            ));
        }
        Ok(Value::Record(RecordValue { fields }))
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        installer,
        "num_brent",
        "crate::kernels::rootfind::brent",
        "Serial bracketed Brent root solver: lower bound, upper bound, positive argument tolerance, positive iteration budget, callback. The root and tolerance have callback-domain dimensions; residual dimensions are independent. Returns root/evaluations without an extra residual evaluation. Endpoint samples count as evaluations, callback failures propagate unchanged, and iteration exhaustion is an explicit error.",
        5,
        Brent { with_bracket: false },
    )?;
    register(
        installer,
        "num_brent_with_bracket",
        "crate::kernels::rootfind::brent_with_bracket",
        "Bracketed Brent root solver returning the root, actual evaluation count and final two endpoints in ascending order. The argument tolerance and all returned coordinates use the callback's input dimension. An exact zero repeats its coordinate at both endpoints.",
        5,
        Brent { with_bracket: true },
    )
}
