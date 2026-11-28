//! Concrete, checked specializations of the native discrete/scientific API.
use super::codecs::Arg;
use super::codecs::Output;
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeOutput;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::error;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use super::numeric;
use super::runtime_dimensions;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::VariantValue;

struct Input(Value);
impl NativeArg for Input {
    fn ty() -> TypeExpr {
        api::dynamic_type()
    }

    fn decode(
        value: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::bounded(value)?;
        Ok(Self(value.clone()))
    }
}
struct Maybe<T>(Option<T>);
impl<T: NativeOutput> NativeOutput for Maybe<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![T::ty()]),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Variant(match self.0 {
            | None => {
                VariantValue {
                    tag: "None".into(),
                    payload: vec![],
                }
            },
            | Some(value) => {
                VariantValue {
                    tag: "Some".into(),
                    payload: vec![NativeOutput::encode(value, ctx)?],
                }
            },
        }))
    }
}
struct Many<T>(Vec<T>);
impl<T: NativeOutput> NativeOutput for Many<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(T::ty()),
        }
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if self.0.len() > api::VALUE_LIMIT {
            return Err(error("native object output count exceeds budget"));
        }
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.0.len()],
            data: self
                .0
                .into_iter()
                .map(|value| NativeOutput::encode(value, ctx))
                .collect::<Result<_, _>>()?,
            elem: T::ty(),
        }))
    }
}
fn limit(
    value: usize,
    max: usize,
) -> Result<usize, OpError> {
    if value > max {
        Err(error("native work/allocation bound exceeded"))
    } else {
        Ok(value)
    }
}
fn scalar(value: &Input) -> Result<(f64, DimExpr), OpError> {
    Ok((numeric(&value.0, 0)?, api::dimension(&value.0)?))
}
fn canonical(dim: &DimExpr) -> Result<DimExpr, OpError> {
    runtime_dimensions::expression(&runtime_dimensions::canonical(dim)?)
}
fn object_dimension<T: NativeType>(
    object: &Handle<T>,
    ctx: &CallCtx<'_>,
) -> Result<DimExpr, OpError> {
    Ok(object
        .dimensions(ctx)?
        .into_iter()
        .next()
        .unwrap_or(DimExpr::Dimless))
}
fn same<T: NativeType, U: NativeType>(
    left: &Handle<T>,
    right: &Handle<U>,
    ctx: &CallCtx<'_>,
) -> Result<DimExpr, OpError> {
    let dim = object_dimension(left, ctx)?;
    api::same_dimension(&dim, &object_dimension(right, ctx)?)?;
    Ok(dim)
}
fn quantities(
    data: Vec<f64>,
    dimension: &DimExpr,
) -> Result<Dynamic, OpError> {
    let mut value = Output::encode(data)?;
    let Value::Tensor(t) = &mut value else {
        return Err(error("native vector encoding invariant"));
    };
    t.elem = TypeExpr::Float.with_dim(dimension.clone());
    t.data = t
        .data
        .drain(..)
        .map(|value| api::quantity(numeric(&value, 0)?, dimension.clone()))
        .collect::<Result<_, _>>()?;
    Ok(Dynamic(value))
}
fn quantities_matrix(
    data: Vec<Vec<f64>>,
    dimension: &DimExpr,
) -> Result<Dynamic, OpError> {
    let mut value = Output::encode(data)?;
    let Value::Tensor(t) = &mut value else {
        return Err(error("native matrix encoding invariant"));
    };
    t.elem = TypeExpr::Float.with_dim(dimension.clone());
    t.data = t
        .data
        .drain(..)
        .map(|value| api::quantity(numeric(&value, 0)?, dimension.clone()))
        .collect::<Result<_, _>>()?;
    Ok(Dynamic(value))
}

mod arithmetic;
mod ball;
mod callbacks;
mod coding;
mod distribution_values;
mod distributions;
mod epd;
mod graph;
mod graphics;
mod homology_maps;
mod homotopy;
mod lie;
mod misc;
mod parallel;
mod permgroup;
mod qmatrix;
mod stochastic;
pub(super) fn install(
    i: &mut Installer<'_>,
    namespace: &str,
) -> Result<(), RuleError> {
    if namespace != "numerical_vector" {
        return Ok(());
    }
    coding::install(i)?;
    graph::install(i)?;
    graphics::install(i)?;
    lie::install(i)?;
    parallel::install(i)?;
    callbacks::install(i)?;
    epd::install(i)?;
    homotopy::install(i)?;
    ball::install(i)?;
    permgroup::install(i)?;
    homology_maps::install(i)?;
    misc::install(i)?;
    qmatrix::install(i)?;
    arithmetic::install(i)?;
    stochastic::install(i)?;
    distributions::install(i)?;
    distribution_values::install(i)?;
    Ok(())
}

fn extended_quantity(
    value: f64,
    dim: DimExpr,
) -> Result<Dynamic, OpError> {
    if value.is_nan() {
        return Err(error("native quantity is mathematically undefined (NaN)"));
    }
    Ok(Dynamic(if dim == DimExpr::Dimless {
        Value::Float(value)
    } else {
        Value::Dyn(crate::kernel::value::DynValue {
            dim,
            value: Box::new(Value::Float(value)),
        })
    }))
}
