//! Explicit, scaled callback adapters for native numerical APIs.
#[path = "native_callbacks/autodiff.rs"]
mod autodiff;
#[path = "native_callbacks/calculus.rs"]
mod calculus;
#[path = "native_callbacks/complex.rs"]
mod complex;
#[path = "native_callbacks/complex_certification.rs"]
mod complex_certification;
#[path = "native_callbacks/contour.rs"]
mod contour;
#[path = "native_callbacks/saddle.rs"]
mod saddle;
#[path = "native_callbacks/transforms.rs"]
mod transforms;

#[path = "native_callbacks/convergence.rs"]
mod convergence;
mod laurent;
#[path = "native_callbacks/optimization.rs"]
mod optimization;
#[path = "native_callbacks/riemann_hilbert.rs"]
mod riemann_hilbert;
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
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::expression;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::value::TensorValue;
use std::cell::Cell;

pub(super) struct Callbacks<'a, 'ctx> {
    pub(super) runs: Vec<CallbackRun<'a, 'ctx>>,
    pub(super) signatures: Vec<Signature>,
    pub(super) used: Cell<usize>,
    pub(super) limit: usize,
    pub(super) exhausted: Cell<bool>,
}
impl<'a, 'ctx> Callbacks<'a, 'ctx> {
    pub(super) fn new(
        ctx: &'a CallCtx<'ctx>,
        args: &[Value],
        signature: &Signature,
        start: usize,
        limit: usize,
    ) -> Result<Self, OpError> {
        arity(args, signature.params.len())?;
        let mut runs = Vec::new();
        let mut signatures = Vec::new();
        for (index, expected) in signature.params.iter().enumerate().skip(start) {
            let TypeExpr::Function(params, _) = expected else {
                return Err(domain("callback parameter is not a function"));
            };
            let (run, sig) = CallbackRun::new(
                ctx,
                &args[index],
                index,
                u8::try_from(params.len()).map_err(|_| domain("callback arity exceeds ABI"))?,
                expected.clone(),
            )?;
            runs.push(run);
            signatures.push(sig);
        }
        Ok(Self {
            runs,
            signatures,
            used: Cell::new(0),
            limit,
            exhausted: Cell::new(false),
        })
    }

    pub(super) fn call<T>(
        &self,
        index: usize,
        args: &[Value],
        decode: impl FnOnce(Value) -> Result<T, OpError>,
    ) -> Option<T> {
        if self.exhausted.get() || self.runs.iter().any(CallbackRun::failed) {
            return None;
        }
        if self.used.get() >= self.limit {
            self.exhausted.set(true);
            return None;
        }
        self.used.set(self.used.get() + 1);
        self.runs[index].evaluate(args, decode)
    }

    pub(super) fn finish(self) -> Result<(), OpError> {
        for run in self.runs {
            run.finish()?;
        }
        if self.exhausted.get() {
            Err(OpError::BudgetExceeded)
        } else {
            Ok(())
        }
    }
}
pub(super) fn domain(message: impl Into<String>) -> OpError {
    OpError::Domain(message.into())
}
const fn arity(
    args: &[Value],
    expected: usize,
) -> Result<(), OpError> {
    if args.len() == expected {
        Ok(())
    } else {
        Err(OpError::Arity {
            expected: ArityInfo::fixed(expected as u8),
            got: args.len(),
        })
    }
}
pub(super) fn typed(
    params: Vec<TypeExpr>,
    ret: TypeExpr,
    dims: &[&str],
) -> Signature {
    let mut sig = effectful(Signature::new(params, ret));
    sig.type_params = dims
        .iter()
        .map(|name| {
            TypeParam {
                name: (*name).into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            }
        })
        .collect();
    sig
}
pub(super) fn real(d: DimExpr) -> TypeExpr {
    dimensional_type(d)
}
pub(super) fn function(
    params: Vec<TypeExpr>,
    ret: TypeExpr,
) -> TypeExpr {
    TypeExpr::Function(params, Box::new(ret))
}
pub(super) fn var(name: &str) -> DimExpr {
    DimExpr::Var(name.into())
}
pub(super) fn ratio(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<DimExpr, OpError> {
    expression(&canonical(&DimExpr::Mul(vec![
        a.clone(),
        DimExpr::Inverse(Box::new(b.clone())),
    ]))?)
}
pub(super) fn equal(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<(), OpError> {
    if canonical(a)? == canonical(b)? {
        Ok(())
    } else {
        Err(domain("callback or argument dimension mismatch"))
    }
}
pub(super) fn scalar(
    v: &Value,
    d: &DimExpr,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    check_dimension(ctx, v, d, index)?;
    numeric(v, index)
}
pub(super) fn positive(
    v: &Value,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    let v = scalar(v, &DimExpr::Dimless, index, ctx)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(domain("control must be finite and positive"))
    }
}
pub(super) fn scale(
    v: &Value,
    d: &DimExpr,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    let v = scalar(v, d, index, ctx)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(domain("physical reference scale must be positive"))
    }
}
pub(super) fn count(
    v: &Value,
    index: usize,
) -> Result<usize, OpError> {
    let n = positive_budget(v, index, "finite native work count")?;
    if n > 1_000_000 {
        Err(domain("native work count exceeds one million"))
    } else {
        Ok(n)
    }
}
pub(super) fn vec_value(
    v: &[f64],
    d: &DimExpr,
) -> Value {
    vector_from_iter(v.iter().copied(), d)
}
pub(super) fn scaled_vec_value(
    coordinates: &[f64],
    reference: f64,
    dimension: &DimExpr,
) -> Value {
    vector_from_iter(coordinates.iter().map(|x| x * reference), dimension)
}
fn vector_from_iter(
    coordinates: impl ExactSizeIterator<Item = f64>,
    dimension: &DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![coordinates.len()],
        data: coordinates.map(Value::Float).collect(),
        elem: real(dimension.clone()),
    })
}
pub(super) fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields.into_iter().map(|(n, v)| (n.into(), v)).collect(),
    })
}
pub(super) fn finite(v: f64) -> Result<f64, OpError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(domain("native numerical result is not finite"))
    }
}
pub(super) fn bounded_size(
    rows: usize,
    cols: usize,
) -> Result<(), OpError> {
    if rows.checked_mul(cols).is_none_or(|n| n > 8_388_608) {
        Err(domain("native allocation exceeds checked storage budget"))
    } else {
        Ok(())
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    calculus::install(i)?;
    transforms::install(i)?;
    autodiff::install(i)?;
    complex::install(i)?;
    complex_certification::install(i)?;
    saddle::install(i)?;
    contour::install(i)?;
    laurent::install(i)?;
    convergence::install(i)?;
    riemann_hilbert::install(i)?;
    optimization::install(i)
}
