//! Higher-order numerical exports through checked, serial host callbacks.
//!
//! These evaluators deliberately have no rewrite/scalar/JIT semantics: a
//! callback may read or mutate host state. The function signature, dimensions,
//! callback errors and explicit quadrature budget remain observable.
mod brent;
mod common;
mod evolution;
mod extended;
mod least_squares;
use super::numerical_abi::numeric;
use super::numerical_abi::runtime_dimensions::canonical;
use super::numerical_abi::runtime_dimensions::expression;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Effect;
use crate::kernel::signature::EffectRow;
use crate::kernel::signature::Purity;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;
use crate::kernel::value::DynValue;
use crate::kernel::value::RecordValue;
use crate::kernel::value::Value;
use crate::kernels::quadrature::GkRule;
use crate::kernels::quadrature::integrate_gk;
use num_traits::ToPrimitive;
use std::cell::RefCell;
use std::sync::Arc;

fn dimensional_type(dimension: DimExpr) -> TypeExpr {
    if dimension == DimExpr::Dimless {
        TypeExpr::Float
    } else {
        TypeExpr::Dim(Box::new(TypeExpr::Float), dimension)
    }
}
fn dimensional_float(
    value: f64,
    dimension: DimExpr,
) -> Value {
    if dimension == DimExpr::Dimless {
        Value::Float(value)
    } else {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim: dimension,
        })
    }
}
fn scalar_dimension(ty: &TypeExpr) -> Result<DimExpr, OpError> {
    match ty {
        | TypeExpr::Float | TypeExpr::Int | TypeExpr::Real => Ok(DimExpr::Dimless),
        | TypeExpr::Dim(inner, dim)
            if matches!(
                inner.as_ref(),
                TypeExpr::Float | TypeExpr::Int | TypeExpr::Real
            ) =>
        {
            canonical(dim)?;
            Ok(dim.clone())
        },
        | _ => {
            Err(OpError::Unsupported {
                backend: "numerical callback".into(),
                feature: "real scalar callback domain and codomain".into(),
            })
        },
    }
}
fn check_dimension(
    ctx: &CallCtx<'_>,
    value: &Value,
    dimension: &DimExpr,
    index: usize,
) -> Result<(), OpError> {
    let actual = ctx.value_type(value)?;
    if canonical(&scalar_dimension(&actual)?)? == canonical(dimension)? {
        return Ok(());
    }
    Err(OpError::Type {
        expected: Box::new(dimensional_type(dimension.clone())),
        got: Box::new(actual),
        arg: index,
    })
}

struct Integrate;
impl EvalDyn for Integrate {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let output = DimExpr::Mul(vec![d.clone(), e.clone()]);
        let mut signature = Signature::new(
            vec![
                dimensional_type(d.clone()),
                dimensional_type(d.clone()),
                dimensional_type(output),
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Function(
                    vec![dimensional_type(d.clone())],
                    Box::new(dimensional_type(e.clone())),
                ),
            ],
            super::numerical_abi::codecs::structured_type::<crate::kernels::quadrature::Integral>(
                super::numerical_abi::codecs::STRUCTURED_INTEGRAL,
                &d,
                &e,
            ),
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
        signature.effects = EffectRow {
            effects: vec![Effect::Stateful, Effect::Diverges],
        };
        signature
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 6 {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(6),
                got: args.len(),
            });
        }
        let (callback, signature) =
            common::CallbackRun::new(ctx, &args[5], 5, 1, self.signature().params[5].clone())?;
        let domain = scalar_dimension(&signature.params[0])?;
        let codomain = scalar_dimension(&signature.ret)?;
        let output = expression(&canonical(&DimExpr::Mul(vec![
            domain.clone(),
            codomain.clone(),
        ]))?)?;
        check_dimension(ctx, &args[0], &domain, 0)?;
        check_dimension(ctx, &args[1], &domain, 1)?;
        let a = numeric(&args[0], 0)?;
        let b = numeric(&args[1], 1)?;
        let abs_tol = numeric(&args[2], 2)?;
        let rel_tol = numeric(&args[3], 3)?;
        // Exact zero contributes no absolute tolerance and is unit-polymorphic.
        if abs_tol != 0.0 {
            check_dimension(ctx, &args[2], &output, 2)?;
        }
        check_dimension(ctx, &args[3], &DimExpr::Dimless, 3)?;
        if abs_tol < 0.0 || rel_tol < 0.0 || (abs_tol == 0.0 && rel_tol == 0.0) {
            return Err(OpError::Domain(
                "quadrature tolerances must be nonnegative and not both zero".into(),
            ));
        }
        let panels = common::positive_budget(&args[4], 4, "quadrature panel budget")?;
        let evaluate = |x: f64| {
            callback
                .evaluate(&[dimensional_float(x, domain.clone())], |value| {
                    numeric(&value, 5)
                })
                .unwrap_or(f64::NAN)
        };
        let result = integrate_gk(evaluate, a, b, GkRule::G7K15, abs_tol, rel_tol, panels);
        callback.finish()?;
        if !result.value.is_finite() || !result.error.is_finite() {
            return Err(OpError::Domain(
                "quadrature result or error estimate is not finite".into(),
            ));
        }
        super::numerical_abi::codecs::structured_value(
            result,
            super::numerical_abi::codecs::STRUCTURED_INTEGRAL,
            &domain,
            &codomain,
        )
    }
}
fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    common::register(
        installer,
        "num_integrate",
        "crate::kernels::quadrature::integrate_gk",
        "Serial adaptive Gauss-Kronrod quadrature: a, b, absolute tolerance, relative tolerance, max panels, callback. Returns value/error/evaluations/converged. Bounds match the callback domain; value/error/absolute tolerance carry domain times codomain dimensions; zero absolute tolerance is polymorphic. Callback errors propagate unchanged, and the estimate is heuristic, not a rigorous enclosure.",
        6,
        Integrate,
    )?;
    brent::install(installer)?;
    evolution::install(installer)?;
    least_squares::install(installer)?;
    extended::install(installer)?;
    Ok(())
}
/// Serial numerical operators with checked host callbacks.
#[must_use]
pub fn numerical_callback() -> RuleSet {
    RuleSet::new("numerical_callback", install)
}
