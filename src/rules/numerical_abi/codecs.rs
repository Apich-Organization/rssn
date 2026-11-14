//! Composable typed codecs for numerical exports.
use super::CallCtx;
use super::DimensionPolicy;
use super::EvalDyn;
use super::NatExpr;
use super::OpError;
use super::Signature;
use super::Stability;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::dense;
use super::dense_error;
use super::dimensions;
use super::matrix;
use super::numeric;
use super::output_matrix;
use super::output_scalar;
use super::output_tensor;
use super::output_vector;
use super::runtime_dimensions;
use super::special;
use super::tensor;
use super::tensor_type;
use crate::kernel::value::RecordValue;
use num_traits::ToPrimitive;

/// A value type that can be decoded from and described in the kernel ABI.
trait Arg: Sized {
    fn ty() -> TypeExpr;
    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError>;
}
/// A result type that can be described and encoded in the kernel ABI.
pub(in crate::rules) trait Output {
    fn ty() -> TypeExpr;
    fn encode(self) -> Result<Value, OpError>;
}
impl Arg for f64 {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        numeric(v, index)
    }
}
impl Output for f64 {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn encode(self) -> Result<Value, OpError> {
        output_scalar(self)
    }
}
impl Arg for u64 {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match v {
            | Value::Int(n) => {
                n.to_u64().ok_or_else(|| {
                    OpError::Domain(format!(
                        "argument {index} must be a nonnegative integer fitting u64"
                    ))
                })
            },
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(TypeExpr::Int),
                    got: Box::new(v.type_expr()),
                    arg: index,
                })
            },
        }
    }
}
impl Output for usize {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::int(self))
    }
}
impl Arg for Vec<f64> {
    fn ty() -> TypeExpr {
        tensor_type(1)
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        Ok(tensor(v, 1, index)?.1)
    }
}
impl Output for Vec<f64> {
    fn ty() -> TypeExpr {
        tensor_type(1)
    }

    fn encode(self) -> Result<Value, OpError> {
        output_vector(self)
    }
}
impl Arg for dense::Mat {
    fn ty() -> TypeExpr {
        tensor_type(2)
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        matrix(v, index)
    }
}
impl Output for dense::Mat {
    fn ty() -> TypeExpr {
        tensor_type(2)
    }

    fn encode(self) -> Result<Value, OpError> {
        output_matrix(self)
    }
}
impl<T: Output> Output for Result<T, dense::DenseError> {
    fn ty() -> TypeExpr {
        T::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        self.map_err(dense_error)?.encode()
    }
}
impl<T: Output> Output for Result<T, OpError> {
    fn ty() -> TypeExpr {
        T::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        self?.encode()
    }
}

impl Output for bool {
    fn ty() -> TypeExpr {
        TypeExpr::Bool
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Bool(self))
    }
}
impl Arg for Vec<usize> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Int),
        }
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = v else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(v.type_expr()),
                arg: index,
            });
        };
        if t.shape.len() != 1 || t.shape[0] != t.data.len() {
            return Err(OpError::Domain(format!(
                "argument {index} must be a well-formed integer vector"
            )));
        }
        t.data
            .iter()
            .map(|value| {
                let n = <u64 as Arg>::decode(value, index)?;
                usize::try_from(n)
                    .map_err(|_| OpError::Domain(format!("argument {index} exceeds platform size")))
            })
            .collect()
    }
}
impl Output for Vec<usize> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            data: self.into_iter().map(Value::int).collect(),
            elem: TypeExpr::Int,
        }))
    }
}
struct TrainCores(Vec<crate::kernels::tensor_train::TrainCore>, f64);
impl Output for TrainCores {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(tensor_type(3)),
        }
    }

    fn encode(self) -> Result<Value, OpError> {
        let mut data = Vec::with_capacity(self.0.len());
        let last = self.0.len().saturating_sub(1);
        for (index, core) in self.0.into_iter().enumerate() {
            let values = if index == last {
                core.data().iter().map(|x| x / self.1).collect()
            } else {
                core.data().to_vec()
            };
            data.push(output_tensor(core.shape().to_vec(), values)?);
        }
        Ok(Value::Tensor(TensorValue {
            shape: vec![data.len()],
            data,
            elem: tensor_type(3),
        }))
    }
}
// One record declaration controls both the signature and the value fields.
macro_rules! record {
    ($ty:path { $($field:ident : $field_ty:ty),+ $(,)? }) => {
        impl Output for $ty {
            fn ty() -> TypeExpr { TypeExpr::Record(vec![$((stringify!($field).to_string(), <$field_ty as Output>::ty())),+]) }
            fn encode(self) -> Result<Value, OpError> {
                let $ty { $($field),+ } = self;
                Ok(Value::Record(RecordValue { fields: vec![$((stringify!($field).to_string(), $field.encode()?)),+] }))
            }
        }
    };
}
// The same concrete codecs serve every source-generated binding.
macro_rules! integer_codec {
    ($($ty:ty),*) => {$(
        impl Arg for $ty {
            fn ty() -> TypeExpr { TypeExpr::Int }
            fn decode(value: &Value, index: usize) -> Result<Self, OpError> {
                let Value::Int(value) = value else { return Err(OpError::Type { expected: Box::new(TypeExpr::Int), got: Box::new(value.type_expr()), arg:index }); };
                value.to_string().parse::<Self>().map_err(|_| OpError::Domain(format!("argument {index} does not fit {}",stringify!($ty))))
            }
        }
    )*};
}
integer_codec!(usize, i64, u32, i32, i128);
macro_rules! integer_output {
    ($($ty:ty),*) => {$(impl Output for $ty {
        fn ty() -> TypeExpr { TypeExpr::Int }
        fn encode(self) -> Result<Value,OpError> { Ok(Value::int(self)) }
    })*};
}
integer_output!(u64, i64, u32, i32, i128, num_bigint::BigInt);
impl Arg for bool {
    fn ty() -> TypeExpr {
        TypeExpr::Bool
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match value {
            | Value::Bool(value) => Ok(*value),
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(TypeExpr::Bool),
                    got: Box::new(value.type_expr()),
                    arg: index,
                })
            },
        }
    }
}
impl Arg for num_bigint::BigInt {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match value {
            | Value::Int(value) => Ok(value.clone()),
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(TypeExpr::Int),
                    got: Box::new(value.type_expr()),
                    arg: index,
                })
            },
        }
    }
}
impl Output for num_complex::Complex64 {
    fn ty() -> TypeExpr {
        TypeExpr::Complex
    }

    fn encode(self) -> Result<Value, OpError> {
        if !self.re.is_finite() || !self.im.is_finite() {
            return Err(OpError::Domain("non-finite complex result".into()));
        }
        Ok(Value::Complex(crate::kernel::value::ComplexValue::Float {
            re: self.re,
            im: self.im,
        }))
    }
}
impl<T: Output> Output for Option<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![T::ty()]),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Variant(crate::kernel::value::VariantValue {
            tag: if self.is_some() {
                "Some"
            } else {
                "None"
            }
            .into(),
            payload: match self {
                | Some(value) => vec![value.encode()?],
                | None => vec![],
            },
        }))
    }
}
impl<T: Arg, const N: usize> Arg for [T; N] {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Lit(N as u64)]),
            elem: Box::new(T::ty()),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if t.shape != [N] || t.data.len() != N {
            return Err(OpError::Domain(format!(
                "argument {index} must have shape [{N}]"
            )));
        }
        t.data
            .iter()
            .map(|v| T::decode(v, index))
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| OpError::Domain("invalid fixed array length".into()))
    }
}
impl<T: Output, const N: usize> Output for [T; N] {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Lit(N as u64)]),
            elem: Box::new(T::ty()),
        }
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![N],
            data: self
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: T::ty(),
        }))
    }
}
impl<T: Arg, const N: usize> Arg for Vec<[T; N]> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(2),
            shape: None,
            elem: Box::new(T::ty()),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if t.shape.len() != 2
            || t.shape.get(1) != Some(&N)
            || t.shape.first().and_then(|n| n.checked_mul(N)) != Some(t.data.len())
            || N == 0
        {
            return Err(OpError::Domain(format!(
                "argument {index} requires a matrix with {N} columns"
            )));
        }
        t.data
            .as_chunks::<N>()
            .0
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| T::decode(v, index))
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .map_err(|_| OpError::Domain("invalid row length".into()))
            })
            .collect()
    }
}
impl<T: Output, const N: usize> Output for Vec<[T; N]> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(2),
            shape: None,
            elem: Box::new(T::ty()),
        }
    }

    fn encode(self) -> Result<Value, OpError> {
        let n = self.len();
        Ok(Value::Tensor(TensorValue {
            shape: vec![n, N],
            data: self
                .into_iter()
                .flatten()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: T::ty(),
        }))
    }
}
impl Arg for Vec<i64> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Int),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if t.shape != [t.data.len()] {
            return Err(OpError::Domain("expected integer vector".into()));
        }
        t.data
            .iter()
            .map(|v| <i64 as Arg>::decode(v, index))
            .collect()
    }
}
impl Output for Vec<i64> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            data: self.into_iter().map(Value::int).collect(),
            elem: TypeExpr::Int,
        }))
    }
}
impl Arg for Vec<Vec<f64>> {
    fn ty() -> TypeExpr {
        tensor_type(2)
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let (shape, data) = tensor(value, 2, index)?;
        if shape[1] == 0 {
            Ok(vec![vec![]; shape[0]])
        } else {
            Ok(data.chunks_exact(shape[1]).map(<[f64]>::to_vec).collect())
        }
    }
}
impl Output for Vec<Vec<f64>> {
    fn ty() -> TypeExpr {
        tensor_type(2)
    }

    fn encode(self) -> Result<Value, OpError> {
        let n = self.len();
        let m = self.first().map_or(0, Vec::len);
        if self.iter().any(|r| r.len() != m) {
            return Err(OpError::Domain("ragged matrix result".into()));
        }
        output_tensor(vec![n, m], self.into_iter().flatten().collect())
    }
}
macro_rules! tuple_output {
    ($($ty:ident:$field:ident:$index:tt),+) => {impl<$($ty:Output),+> Output for ($($ty,)+) {
        fn ty()->TypeExpr {TypeExpr::Record(vec![$((stringify!($field).into(),$ty::ty())),+])}
        fn encode(self)->Result<Value,OpError>{Ok(Value::Record(RecordValue{fields:vec![$((stringify!($field).into(),self.$index.encode()?)),+]}))}
    }};
}
tuple_output!(A:item0:0,B:item1:1);
tuple_output!(A:item0:0,B:item1:1,C:item2:2);
tuple_output!(A:item0:0,B:item1:1,C:item2:2,D:item3:3);

impl Arg for num_complex::Complex64 {
    fn ty() -> TypeExpr {
        TypeExpr::Complex
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let result = match value {
            | Value::Complex(c) => {
                let (re, im) = c.to_f64();
                Self::new(re, im)
            },
            | _ => Self::new(numeric(value, index)?, 0.0),
        };
        if !result.re.is_finite() || !result.im.is_finite() {
            return Err(OpError::Domain(format!(
                "argument {index} contains a non-finite complex value"
            )));
        }
        Ok(result)
    }
}
impl Arg for Vec<num_complex::Complex64> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Complex),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if t.shape != [t.data.len()] {
            return Err(OpError::Domain("expected a complex vector".into()));
        }
        t.data
            .iter()
            .map(|v| <num_complex::Complex64 as Arg>::decode(v, index))
            .collect()
    }
}
impl Output for Vec<num_complex::Complex64> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            data: self
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: TypeExpr::Complex,
        }))
    }
}
impl<T: Arg> Arg for Option<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![T::ty()]),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match value {
            | Value::Variant(v) if v.tag == "None" && v.payload.is_empty() => Ok(None),
            | Value::Variant(v) if v.tag == "Some" && v.payload.len() == 1 => {
                Ok(Some(T::decode(&v.payload[0], index)?))
            },
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(<Self as Arg>::ty()),
                    got: Box::new(value.type_expr()),
                    arg: index,
                })
            },
        }
    }
}

impl Arg for String {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match value {
            | Value::String(value) => Ok(value.clone()),
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(TypeExpr::String),
                    got: Box::new(value.type_expr()),
                    arg: index,
                })
            },
        }
    }
}
impl Arg for Vec<String> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::String),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if t.shape != [t.data.len()] {
            return Err(OpError::Domain("expected a vector of names".into()));
        }
        t.data
            .iter()
            .map(|v| <String as Arg>::decode(v, index))
            .collect()
    }
}

include!("structured_generated.rs");

/// Build a record's physical result schema from generated field formulas.
pub(in crate::rules) fn structured_type<T: Output>(
    fields: &[(&str, super::dimensions::FieldDimension)],
    d: &crate::kernel::signature::DimExpr,
    e: &crate::kernel::signature::DimExpr,
) -> TypeExpr {
    let mut ty = T::ty();
    super::dimensions::apply_record_fields(&mut ty, fields, d, e);
    ty
}

/// Encode a validated record and apply the same field formulas as its schema.
pub(in crate::rules) fn structured_value<T: Output>(
    value: T,
    fields: &[(&str, super::dimensions::FieldDimension)],
    d: &crate::kernel::signature::DimExpr,
    e: &crate::kernel::signature::DimExpr,
) -> Result<Value, OpError> {
    if fields
        .iter()
        .any(|(_, formula)| matches!(formula, super::dimensions::FieldDimension::Unbound))
    {
        return Err(OpError::Unsupported {
            backend: "numerical".into(),
            feature: "unbound structured dimensional policy".into(),
        });
    }
    let signature = Signature::new(Vec::new(), structured_type::<T>(fields, d, e));
    runtime_dimensions::finish(
        &signature,
        &std::collections::BTreeMap::new(),
        value.encode()?,
    )
}

struct Qr {
    q: dense::Mat,
    r: dense::Mat,
}
record!(Qr {
    q: dense::Mat,
    r: dense::Mat
});
struct SymmetricEigen {
    values: Vec<f64>,
    vectors: dense::Mat,
}
record!(SymmetricEigen { values: Vec<f64>, vectors: dense::Mat });

/// A descriptor whose arity, types, and conversion are generated together.
pub struct CodecExport {
    /// Public Rust kernel function invoked by the adapter.
    pub rust_path: &'static str,
    /// Name as exposed to hosts.
    pub name: &'static str,
    /// Rule-set namespace.
    pub namespace: &'static str,
    /// Description of domain, output and conventions.
    pub doc: &'static str,
    /// Stability of this export.
    pub stability: Stability,
    /// Shared argument/result dimensional relationship.
    pub dimension_policy: DimensionPolicy,
    base_signature: fn() -> Signature,
    call: fn(&[Value]) -> Result<Value, OpError>,
}
impl CodecExport {
    /// The exact signature derived from codecs and dimensional policy.
    #[must_use]
    pub fn signature(&self) -> Signature {
        dimensions::apply((self.base_signature)(), self.dimension_policy)
    }
}
// Compound raw-number arguments must not hide physical quantities in Some,
// records, or tensor elements. Top-level dimension inference cannot see these.
fn require_dimensionless_type(ty: &TypeExpr) -> Result<(), OpError> {
    match ty {
        | TypeExpr::Dim(inner, dim) => {
            if !runtime_dimensions::canonical(dim)?.is_empty() {
                return Err(OpError::Domain(
                    "raw numerical arguments must be dimensionless".into(),
                ));
            }
            require_dimensionless_type(inner)
        },
        | TypeExpr::Point(_) => {
            Err(OpError::Domain(
                "raw numerical arguments cannot be affine points".into(),
            ))
        },
        | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => {
            require_dimensionless_type(elem)
        },
        | TypeExpr::Record(fields) => {
            for (_, ty) in fields {
                require_dimensionless_type(ty)?;
            }
            Ok(())
        },
        | TypeExpr::Variant(variants) => {
            for (_, payload) in variants {
                for ty in payload {
                    require_dimensionless_type(ty)?;
                }
            }
            Ok(())
        },
        | _ => Ok(()),
    }
}
fn require_dimensionless(value: &Value) -> Result<(), OpError> {
    require_dimensionless_type(&value.type_expr())?;
    match value {
        | Value::Quantity(quantity) => {
            if !runtime_dimensions::canonical(&quantity.unit.dimension)?.is_empty() {
                return Err(OpError::Domain(
                    "raw numerical arguments must be dimensionless".into(),
                ));
            }
            require_dimensionless(&quantity.value)
        },
        | Value::Dyn(dynamic) => {
            if !runtime_dimensions::canonical(&dynamic.dim)?.is_empty() {
                return Err(OpError::Domain(
                    "raw numerical arguments must be dimensionless".into(),
                ));
            }
            require_dimensionless(&dynamic.value)
        },
        | Value::Tensor(tensor) => {
            for value in &tensor.data {
                require_dimensionless(value)?;
            }
            Ok(())
        },
        | Value::Record(record) => {
            for (_, value) in &record.fields {
                require_dimensionless(value)?;
            }
            Ok(())
        },
        | Value::Variant(variant) => {
            for value in &variant.payload {
                require_dimensionless(value)?;
            }
            Ok(())
        },
        | _ => Ok(()),
    }
}

impl EvalDyn for &'static CodecExport {
    fn signature(&self) -> Signature {
        CodecExport::signature(self)
    }

    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let arity = self.signature().params.len();
        if args.len() != arity {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(arity as u8),
                got: args.len(),
            });
        }
        let contract = dimensions::contract_for(self.dimension_policy);
        if contract.other_args_dimensionless {
            for (index, value) in args.iter().enumerate() {
                if !contract.inputs.iter().any(|input| input.argument == index)
                    && !contract
                        .constraints
                        .iter()
                        .any(|input| input.argument == index)
                {
                    require_dimensionless(value)?;
                }
            }
        }
        let dimensions = runtime_dimensions::prepare(self.dimension_policy, args)?;
        let result = (self.call)(args)?;
        runtime_dimensions::finish(&self.signature(), &dimensions, result)
    }
}
// A descriptor contains only the kernel call. Its codecs generate both the
// signature and all checked argument/result conversions.
macro_rules! binding {
    ($module:ident, $name:literal, $namespace:literal, $doc:literal, $rust_path:literal, $policy:expr,
     ($($arg:ident : $argty:ty),* $(,)?) -> $ret:ty, $body:expr) => {
        mod $module {
            // Descriptor bodies use different subsets of these names.
            #[allow(unused_imports)]
            use super::{Arg, Output, CodecExport, Signature, Value, OpError, Stability,
                DimensionPolicy, dense, dense_error, special, Qr, SymmetricEigen, TrainResult, TrainCores};
            fn signature() -> Signature { Signature::new(vec![$(<$argty as Arg>::ty()),*], <$ret as Output>::ty()) }
            fn call(args: &[Value]) -> Result<Value, OpError> {
                let mut index = 0usize;
                $(let $arg: $argty = { let value = <$argty as Arg>::decode(&args[index], index)?; index += 1; value };)*
                let _ = index;
                let result: $ret = ($body)($($arg),*);
                result.encode()
            }
            pub static EXPORT: CodecExport = CodecExport { rust_path: $rust_path, name: $name, namespace: $namespace,
                doc: $doc, stability: Stability::Experimental, dimension_policy: $policy, base_signature: signature, call };
        }
    };
}
binding!(svd, "num_svd", "numerical_matrix", "Thin SVD A = U diag(s) V^T; returns u, s and v.", "crate::kernels::dense::svd", DimensionPolicy::Svd,
    (a: dense::Mat) -> Result<dense::Svd, dense::DenseError>, |a: dense::Mat| dense::svd(&a));
binding!(qr, "num_qr", "numerical_matrix", "Householder QR with q square and r matching the input.", "crate::kernels::dense::qr", DimensionPolicy::Qr,
    (a: dense::Mat) -> Qr, |a: dense::Mat| { let (q,r) = dense::qr(&a); Qr {q,r} });
binding!(eigen, "num_eigen_symmetric", "numerical_matrix", "Real eigenvalues and orthonormal eigenvectors of a symmetric matrix.", "crate::kernels::dense::eigen_symmetric", DimensionPolicy::SymmetricEigen,
    (a: dense::Mat) -> Result<SymmetricEigen, dense::DenseError>, |a: dense::Mat| dense::eigen_symmetric(&a).map(|(values,vectors)| SymmetricEigen {values,vectors}));
binding!(lstsq, "num_lstsq", "numerical_matrix", "Full-column-rank least-squares solution using QR.", "crate::kernels::dense::lstsq", DimensionPolicy::LeastSquares,
    (a: dense::Mat, b: Vec<f64>) -> Result<Vec<f64>, dense::DenseError>, |a: dense::Mat, b: Vec<f64>| dense::lstsq(&a,&b));
binding!(rank, "num_rank_tol", "numerical_matrix", "SVD rank using a nonnegative absolute singular-value tolerance.", "crate::kernels::dense::rank", DimensionPolicy::RankTolerance,
(a: dense::Mat, tol: f64) -> Result<usize, OpError>, |a: dense::Mat, tol: f64| {
    if tol < 0.0 { return Err(OpError::Domain("rank tolerance must be nonnegative".into())); }
    dense::rank(&a,Some(tol)).map_err(dense_error)
});
binding!(factorial, "num_factorial", "numerical_special", "Floating-point factorial of a nonnegative u64 integer.", "crate::kernels::special::factorial", DimensionPolicy::Dimensionless,
(n: u64) -> Result<f64, OpError>, |n: u64| {
    if n == u64::MAX { return Err(OpError::Domain("factorial argument exceeds implementation range".into())); }
    Ok(special::factorial(n))
});
binding!(binomial, "num_binomial", "numerical_special", "Floating-point binomial coefficient; returns zero when k exceeds n.", "crate::kernels::special::binomial", DimensionPolicy::Dimensionless,
(n: u64, k: u64) -> Result<f64, OpError>, |n: u64, k: u64| {
    if n == u64::MAX { return Err(OpError::Domain("binomial argument exceeds implementation range".into())); }
    Ok(special::binomial(n,k))
});
struct TrainResult {
    cores: TrainCores,
    shape: Vec<usize>,
    scale: f64,
    error: f64,
    tolerance_met: bool,
}
record!(TrainResult { cores: TrainCores, shape: Vec<usize>, scale: f64, error: f64, tolerance_met: bool });
binding!(tt_svd, "num_tt_svd", "numerical_matrix", "Tensor-train SVD of row-major data; rank cap must be positive. Returns dimensionless normalized cores, shape, scale, Frobenius truncation error and tolerance status; reconstruction multiplies the core contraction by scale.", "crate::kernels::tensor_train::tt_svd", DimensionPolicy::TensorTrain,
(data: Vec<f64>, shape: Vec<usize>, tolerance: f64, max_rank: u64) -> Result<TrainResult, OpError>,
|data: Vec<f64>, shape: Vec<usize>, tolerance: f64, max_rank: u64| {
    use crate::kernels::tensor_train::{tt_svd, TrainError};
    let cap = usize::try_from(max_rank).map_err(|_| OpError::Domain("rank cap exceeds platform size".into()))?;
    let scale = data.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    let scale = if scale == 0.0 { 1.0 } else { scale };
    let result = tt_svd(&data,&shape,tolerance,Some(cap)).map_err(|error| match error {
        TrainError::NoConvergence => OpError::NotReduced(error.to_string()),
        other => OpError::Domain(other.to_string()),
    })?;
    Ok(TrainResult { cores: TrainCores(result.cores().to_vec(), scale), shape: result.shape().to_vec(), scale, error: result.truncation_error(), tolerance_met: result.tolerance_met() })
});
/// Additional exports using composable codecs and structured return values.
pub static EXPORTS: &[&CodecExport] = &[
    &svd::EXPORT,
    &qr::EXPORT,
    &eigen::EXPORT,
    &lstsq::EXPORT,
    &rank::EXPORT,
    &factorial::EXPORT,
    &binomial::EXPORT,
    &tt_svd::EXPORT,
];

include!("codec_generated.rs");
