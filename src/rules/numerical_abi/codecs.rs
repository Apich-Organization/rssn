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
use num_complex::Complex64;
use num_traits::ToPrimitive;
use std::sync::OnceLock;

/// A value type that can be decoded from and described in the kernel ABI.
pub(in crate::rules) trait Arg: Sized {
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
                Self::try_from(value).map_err(|_| OpError::Domain(format!("argument {index} does not fit {}",stringify!($ty))))
            }
        }
    )*};
}
integer_codec!(usize, i64, u32, i32, i128, u8, u128);
macro_rules! integer_output {
    ($($ty:ty),*) => {$(impl Output for $ty {
        fn ty() -> TypeExpr { TypeExpr::Int }
        fn encode(self) -> Result<Value,OpError> { Ok(Value::int(self)) }
    })*};
}
integer_output!(u64, i64, u32, i32, i128, u8, u128, num_bigint::BigInt);
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
// Keep fixed-width rows on the stack and stop on the first bad cell.
fn decode_fixed_array<T: Arg, const N: usize>(
    cells: &[Value; N],
    index: usize,
) -> Result<[T; N], OpError> {
    let mut decoded: [Option<T>; N] = std::array::from_fn(|_| None);
    for (slot, cell) in decoded.iter_mut().zip(cells) {
        *slot = Some(T::decode(cell, index)?);
    }
    Ok(decoded.map(|item| item.expect("every fixed array slot was filled")))
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
        let cells = t
            .data
            .as_slice()
            .try_into()
            .map_err(|_| OpError::Domain("invalid fixed array length".into()))?;
        decode_fixed_array(cells, index)
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
            .map(|row| decode_fixed_array(row, index))
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
            | Value::Dyn(dynamic) => return Self::decode(&dynamic.value, index),
            | Value::Quantity(quantity) => {
                let scale = super::linear_unit_scale(&quantity.unit, index)?;
                Self::decode(&quantity.value, index)? * scale
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

impl Output for Vec<String> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        if self.len() > 65536
            || self
                .iter()
                .try_fold(0usize, |n, v| n.checked_add(v.len()))
                .is_none_or(|n| n > 1_048_576)
        {
            return Err(OpError::Domain(
                "native text-vector allocation budget exceeded".into(),
            ));
        }
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            elem: TypeExpr::String,
            data: self.into_iter().map(Value::String).collect(),
        }))
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
    /// Public Rust function invoked by the adapter; empty for trait operators
    /// without a public inherent-function source path.
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
    signature: OnceLock<Signature>,
    call: fn(&[Value]) -> Result<Value, OpError>,
}
impl CodecExport {
    /// The exact signature derived from codecs and dimensional policy.
    #[must_use]
    pub fn signature(&self) -> Signature {
        self.checked_signature().clone()
    }

    fn checked_signature(&self) -> &Signature {
        self.signature
            .get_or_init(|| dimensions::apply((self.base_signature)(), self.dimension_policy))
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
        let signature = self.checked_signature();
        let arity = signature.params.len();
        if args.len() != arity {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(arity as u8),
                got: args.len(),
            });
        }
        for value in args {
            super::native_object_api::bounded(value)?;
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
        native_scalar_preflight(self.name, args)?;
        let result = (self.call)(args)?;
        runtime_dimensions::finish_shaped(&signature.ret, &dimensions, args, result)
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
                // A zero-arity export has no argument bindings.
                #[allow(unused_mut)]
                let mut index = 0usize;
                let _ = args;
                $(let $arg: $argty = { let value = <$argty as Arg>::decode(&args[index], index)?; index += 1; value };)*
                let _ = index;
                let result: $ret = ($body)($($arg),*);
                result.encode()
            }
            pub static EXPORT: CodecExport = CodecExport { rust_path: $rust_path, name: $name, namespace: $namespace,
                doc: $doc, stability: Stability::Experimental, dimension_policy: $policy, base_signature: signature, signature: std::sync::OnceLock::new(), call };
        }
    };
}
binding!(svd, "num_svd", "numerical_matrix", "Thin SVD A = U diag(s) V^T; returns u, s and v.", "crate::kernels::dense::svd", DimensionPolicy::Svd,
    (a: dense::Mat) -> Result<dense::Svd, dense::DenseError>, |a: dense::Mat| dense::svd(&a));
binding!(qr, "num_qr", "numerical_matrix", "Householder QR with q square and r matching the input.", "crate::kernels::dense::qr", DimensionPolicy::Qr,
    (a: dense::Mat) -> Qr, |a: dense::Mat| { let (q,r) = dense::qr(&a); Qr {q,r} });

binding!(
    qr_thin,
    "num_qr_thin",
    "numerical_matrix",
    "Economy Householder QR. For an m by n matrix, q is m by min(m,n) and dimensionless, and r is min(m,n) by n with the input units. Unrepresentable factors are errors.",
    "crate::kernels::dense::qr_thin",
    DimensionPolicy::Qr,
    (a: dense::Mat) -> Result<Qr, OpError>,
    |a: dense::Mat| {
        let (q, r) = dense::qr_thin(&a).map_err(dense_error)?;
        Ok(Qr { q, r })
    }
);
binding!(
    solve_linear_rows,
    "num_solve_linear_rows",
    "numerical_matrix",
    "Solves a finite dense square system A x = b. The solution has the units of b divided by the units of A.",
    "crate::kernels::dense::solve_linear_rows",
    DimensionPolicy::LeastSquares,
    (a: Vec<Vec<f64>>, b: Vec<f64>) -> Result<Vec<f64>, OpError>,
    |a: Vec<Vec<f64>>, b: Vec<f64>| {
        dense::solve_linear_rows(&a, &b).map_err(dense_error)
    }
);
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
(n: u64, k: u64) -> f64, special::binomial);
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
    &csr_transpose::EXPORT,
    &csr_dense_product::EXPORT,
    &csr_sparse_product::EXPORT,
    &csr_lsqr::EXPORT,
    &qr_thin::EXPORT,
    &solve_linear_rows::EXPORT,
    &rh_matrix_new::EXPORT,
    &rh_matrix_identity::EXPORT,
    &rh_matrix_zero::EXPORT,
    &rh_matrix_add::EXPORT,
    &rh_matrix_sub::EXPORT,
    &rh_matrix_mul::EXPORT,
    &rh_matrix_scale::EXPORT,
    &rh_matrix_inv::EXPORT,
    &rh_matrix_norm::EXPORT,
    &rh_solution_eval::EXPORT,
    &eigen::EXPORT,
    &lstsq::EXPORT,
    &rank::EXPORT,
    &complex_disk_new::EXPORT,
    &complex_disk_of::EXPORT,
    &complex_disk_point::EXPORT,
    &complex_disk_zero::EXPORT,
    &complex_disk_one::EXPORT,
    &complex_disk_i::EXPORT,
    &complex_disk_center::EXPORT,
    &complex_disk_radius::EXPORT,
    &complex_disk_contains::EXPORT,
    &complex_disk_contains_zero::EXPORT,
    &complex_disk_is_finite::EXPORT,
    &complex_disk_erf::EXPORT,
    &complex_disk_log_gamma::EXPORT,
    &complex_disk_gamma::EXPORT,
    &complex_disk_bessel_j::EXPORT,
    &complex_disk_inside_strictly::EXPORT,
    &complex_disk_overlaps::EXPORT,
    &complex_disk_modulus_lower::EXPORT,
    &complex_disk_modulus_upper::EXPORT,
    &complex_disk_neg::EXPORT,
    &complex_disk_conj::EXPORT,
    &complex_disk_add::EXPORT,
    &complex_disk_sub::EXPORT,
    &complex_disk_mul::EXPORT,
    &complex_disk_square::EXPORT,
    &complex_disk_scale::EXPORT,
    &complex_disk_recip::EXPORT,
    &complex_disk_div::EXPORT,
    &complex_disk_jet_variable::EXPORT,
    &complex_disk_jet_constant::EXPORT,
    &complex_disk_jet_add::EXPORT,
    &complex_disk_jet_sub::EXPORT,
    &complex_disk_jet_mul::EXPORT,
    &complex_disk_jet_div::EXPORT,
    &complex_disk_jet_powi::EXPORT,
    &complex_disk_matvec::EXPORT,
    &complex_disk_matmul::EXPORT,
    &complex_disk_solve::EXPORT,
    &complex_disk_jet_exp::EXPORT,
    &complex_disk_jet_sin::EXPORT,
    &complex_disk_jet_cos::EXPORT,
    &complex_disk_jet_sinh::EXPORT,
    &complex_disk_jet_cosh::EXPORT,
    &complex_disk_jet_ln::EXPORT,
    &complex_disk_jet_sqrt::EXPORT,
    &complex_disk_jet_tan::EXPORT,
    &complex_disk_jet_tanh::EXPORT,
    &complex_disk_jet_asin::EXPORT,
    &complex_disk_jet_acos::EXPORT,
    &complex_disk_jet_atan::EXPORT,
    &complex_disk_jet_asinh::EXPORT,
    &complex_disk_jet_acosh::EXPORT,
    &complex_disk_jet_atanh::EXPORT,
    &complex_disk_jet_pow::EXPORT,
    &complex_disk_sec::EXPORT,
    &complex_disk_csc::EXPORT,
    &complex_disk_cot::EXPORT,
    &complex_disk_sech::EXPORT,
    &complex_disk_csch::EXPORT,
    &complex_disk_coth::EXPORT,
    &factorial::EXPORT,
    &continuation_pade::EXPORT,
    &glm_fitted_model::EXPORT,
    &glm_model_predict::EXPORT,
    &survival_fit_cox::EXPORT,
    &survival_predict_cox::EXPORT,
    &survival_kaplan_meier::EXPORT,
    &survival_median::EXPORT,
    &glm_link_value::EXPORT,
    &glm_inverse_link_value::EXPORT,
    &glm_link_derivative::EXPORT,
    &glm_variance::EXPORT,
    &special_complex_legendre_p::EXPORT,
    &special_complex_associated_legendre_p::EXPORT,
    &special_complex_chebyshev_t::EXPORT,
    &special_complex_chebyshev_u::EXPORT,
    &special_complex_hermite_h::EXPORT,
    &special_complex_laguerre_l::EXPORT,
    &special_complex_jacobi_p::EXPORT,
    &special_complex_gegenbauer_c::EXPORT,
    &special_complex_spherical_harmonic::EXPORT,
    &special_complex_bernoulli::EXPORT,
    &special_complex_euler::EXPORT,
    &special_complex_clausen::EXPORT,
    &special_complex_log_barnes_g::EXPORT,
    &special_complex_fresnel_auxiliary::EXPORT,
    &special_complex_bessel_j_order::EXPORT,
    &special_complex_bessel_i_order::EXPORT,
    &special_complex_anger_j::EXPORT,
    &special_complex_weber_e::EXPORT,
    &special_complex_kelvin::EXPORT,
    &special_complex_exp_integral_en::EXPORT,
    &special_complex_legendre_q::EXPORT,
    &continuation_eval_pade::EXPORT,
    &continuation_wynn::EXPORT,
    &continuation_recenter::EXPORT,
    &continuation_holonomic_order::EXPORT,
    &continuation_holonomic_follow::EXPORT,
    &continuation_holonomic_monodromy::EXPORT,
    &continuation_borel_pade::EXPORT,
    &binomial::EXPORT,
    &harmonic::EXPORT,
    &log_binomial::EXPORT,
    &checked_multinomial::EXPORT,
    &complex_dynamics_mandelbrot_escape::EXPORT,
    &complex_dynamics_julia_escape::EXPORT,
    &complex_dynamics_in_main_components::EXPORT,
    &complex_dynamics_mandelbrot_distance::EXPORT,
    &complex_dynamics_quadratic_periodic_points::EXPORT,
    &complex_dynamics_bottcher::EXPORT,
    &complex_dynamics_external_ray::EXPORT,
    &complex_dynamics_quadratic_orbit::EXPORT,
    &complex_dynamics_quadratic_multiplier::EXPORT,
    &complex_dynamics_quadratic_lyapunov::EXPORT,
    &inverse_beta_reg::EXPORT,
    &tt_svd::EXPORT,
    &vector_axpy::EXPORT,
    &vector_axpby::EXPORT,
    &vector_outer_product::EXPORT,
    &signal_filtfilt::EXPORT,
    &contour_deform::EXPORT,
    &contour_keyhole::EXPORT,
    &contour_hankel::EXPORT,
    &contour_dogbone::EXPORT,
    &contour_indented_semicircle::EXPORT,
    &contour_bromwich::EXPORT,
    &contour_rectangle::EXPORT,
    &contour_sector::EXPORT,
    &timeseries_integrate::EXPORT,
    &timeseries_rolling_moments::EXPORT,
    &timeseries_kalman_filter::EXPORT,
    &timeseries_kalman_smooth::EXPORT,
    &control_c2d_zoh::EXPORT,
    &control_simulate_discrete::EXPORT,
    &control_simulate_continuous::EXPORT,
    &control_discrete_n_states::EXPORT,
    &control_discrete_n_inputs::EXPORT,
    &control_discrete_n_outputs::EXPORT,
    &control_discrete_poles::EXPORT,
    &control_discrete_is_stable::EXPORT,
    &control_discrete_controllability_matrix::EXPORT,
    &control_discrete_observability_matrix::EXPORT,
    &polynomial_barycentric_eval::EXPORT,
    &polynomial_power::EXPORT,
    &polynomial_gcd::EXPORT,
    &polynomial_durand_kerner::EXPORT,
    &polynomial_aberth_ehrlich::EXPORT,
    &polynomial_laguerre::EXPORT,
    &polynomial_isolate_verified::EXPORT,
    &cauchy_radius::EXPORT,
    &sharp_cauchy_radius::EXPORT,
    &fujiwara_radius::EXPORT,
    &enestrom_kakeya::EXPORT,
    &koebe_growth::EXPORT,
    &koebe_distortion::EXPORT,
    &koebe_rotation::EXPORT,
    &complex_coord_arg_branch::EXPORT,
    &complex_coord_log_cut::EXPORT,
    &complex_coord_sqrt_cut::EXPORT,
    &complex_coord_pow_cut::EXPORT,
    &carlson_rf_complex::EXPORT,
    &carlson_rc_complex::EXPORT,
    &carlson_rd_complex::EXPORT,
    &carlson_rg_complex::EXPORT,
    &carlson_rj_complex::EXPORT,
    &weierstrass_p_prime::EXPORT,
    &weierstrass_invariants::EXPORT,
    &weierstrass_sigma::EXPORT,
    &weierstrass_zeta::EXPORT,
    &complex_modular_discriminant::EXPORT,
    &mp_real_exp::EXPORT,
    &mp_ball_exp::EXPORT,
    &mp_real_zero::EXPORT,
    &mp_real_one::EXPORT,
    &mp_real_from_i64::EXPORT,
    &mp_real_from_bigint::EXPORT,
    &mp_real_from_f64::EXPORT,
    &mp_real_to_f64::EXPORT,
    &mp_real_sign::EXPORT,
    &mp_decimal_from_text::EXPORT,
    &mp_decimal_from_text_in::EXPORT,
    &mp_decimal_to_text::EXPORT,
    &mp_decimal_to_text_in::EXPORT,
    &mp_decimal_round_text::EXPORT,
    &mp_decimal_enclose::EXPORT,
    &mp_decimal_format_ball::EXPORT,
    &mp_decimal_exponent::EXPORT,
    &mp_decimal_significant_digits::EXPORT,
    &mp_decimal_working_bits::EXPORT,
    &mp_real_round::EXPORT,
    &mp_real_compare::EXPORT,
    &mp_real_add::EXPORT,
    &mp_real_sub::EXPORT,
    &mp_real_mul::EXPORT,
    &mp_real_div::EXPORT,
    &mp_ball_add::EXPORT,
    &mp_ball_sub::EXPORT,
    &mp_ball_mul::EXPORT,
    &mp_ball_div::EXPORT,
    &mp_real_const_pi::EXPORT,
    &mp_real_const_e::EXPORT,
    &mp_real_const_ln2::EXPORT,
    &mp_real_const_sqrt2::EXPORT,
    &mp_real_const_phi::EXPORT,
    &mp_real_const_named::EXPORT,
    &mp_real_const_enclose::EXPORT,
    &mp_real_const_names::EXPORT,
    &mp_real_const_pi_scaled::EXPORT,
    &mp_real_const_functions::EXPORT,
    &mp_real_abs::EXPORT,
    &mp_real_neg::EXPORT,
    &mp_real_is_zero::EXPORT,
    &mp_real_ln::EXPORT,
    &mp_ball_ln::EXPORT,
    &mp_real_log2::EXPORT,
    &mp_ball_log2::EXPORT,
    &mp_real_log10::EXPORT,
    &mp_ball_log10::EXPORT,
    &mp_real_exp2::EXPORT,
    &mp_ball_exp2::EXPORT,
    &mp_real_scale_decimal::EXPORT,
    &mp_real_scale_decimal_in::EXPORT,
    &mp_real_is_integer::EXPORT,
    &mp_real_floor_log2::EXPORT,
    &mp_real_trunc_to_bigint::EXPORT,
    &mp_real_with_precision::EXPORT,
    &mp_real_sqrt::EXPORT,
    &mp_real_sin::EXPORT,
    &mp_ball_sin::EXPORT,
    &mp_real_cos::EXPORT,
    &mp_ball_cos::EXPORT,
    &mp_real_tan::EXPORT,
    &mp_ball_tan::EXPORT,
    &mp_real_sin_cos::EXPORT,
    &mp_real_scale_pow2::EXPORT,
    &mp_ball_new::EXPORT,
    &mp_ball_point::EXPORT,
    &mp_ball_contains_zero::EXPORT,
    &mp_ball_width::EXPORT,
    &mp_real_mantissa::EXPORT,
    &mp_real_exponent::EXPORT,
    &mp_real_precision::EXPORT,
    &mp_ball_midpoint::EXPORT,
    &mp_ball_radius::EXPORT,
    &mp_real_limits::EXPORT,
    &mp_real_rounding_mode::EXPORT,
    &mp_real_cbrt::EXPORT,
    &mp_real_cbrt_in::EXPORT,
    &mp_real_hypot::EXPORT,
    &mp_ball_sqrt::EXPORT,
    &mp_ball_neg::EXPORT,
    &mp_real_sinh::EXPORT,
    &mp_ball_sinh::EXPORT,
    &mp_real_cosh::EXPORT,
    &mp_ball_cosh::EXPORT,
    &mp_real_tanh::EXPORT,
    &mp_ball_tanh::EXPORT,
    &mp_ball_abs::EXPORT,
    &mp_ball_endpoints::EXPORT,
    &mp_ball_lipschitz::EXPORT,
    &mp_expression_enclose::EXPORT,
    &mp_complex_new::EXPORT,
    &mp_complex_point::EXPORT,
    &mp_complex_from_real::EXPORT,
    &mp_complex_zero::EXPORT,
    &mp_complex_one::EXPORT,
    &mp_complex_i::EXPORT,
    &mp_complex_sin::EXPORT,
    &mp_complex_cos::EXPORT,
    &mp_complex_tan::EXPORT,
    &mp_complex_sinh::EXPORT,
    &mp_complex_cosh::EXPORT,
    &mp_complex_tanh::EXPORT,
    &mp_complex_asin::EXPORT,
    &mp_complex_acos::EXPORT,
    &mp_complex_atan::EXPORT,
    &mp_complex_asinh::EXPORT,
    &mp_complex_acosh::EXPORT,
    &mp_complex_atanh::EXPORT,
    &mp_complex_erf::EXPORT,
    &mp_complex_erfc::EXPORT,
    &mp_complex_bessel_j::EXPORT,
    &mp_complex_log_gamma::EXPORT,
    &mp_complex_gamma::EXPORT,
    &mp_complex_digamma::EXPORT,
    &mp_complex_riemann_zeta::EXPORT,
    &mp_complex_polylog::EXPORT,
    &mp_complex_sqrt::EXPORT,
    &mp_complex_pow::EXPORT,
    &mp_complex_abs::EXPORT,
    &mp_complex_exp::EXPORT,
    &mp_complex_ln::EXPORT,
    &mp_complex_neg::EXPORT,
    &mp_complex_conj::EXPORT,
    &mp_complex_add::EXPORT,
    &mp_complex_sub::EXPORT,
    &mp_complex_mul::EXPORT,
    &mp_complex_div::EXPORT,
    &mp_complex_recip::EXPORT,
    &mp_real_atan::EXPORT,
    &mp_ball_atan::EXPORT,
    &mp_real_atan2::EXPORT,
    &mp_ball_atan2::EXPORT,
    &mp_real_asin::EXPORT,
    &mp_ball_asin::EXPORT,
    &mp_real_acos::EXPORT,
    &mp_ball_acos::EXPORT,
    &mp_complex_arg::EXPORT,
    &mp_complex_abs_sq::EXPORT,
    &mp_real_pow::EXPORT,
    &mp_ball_pow::EXPORT,
    &mp_real_powi::EXPORT,
    &mp_ball_powi::EXPORT,
    &mp_real_pow_ratio::EXPORT,
    &mp_ball_pow_ratio::EXPORT,
    &mp_complex_powi::EXPORT,
    &finite_field_batch_inverse::EXPORT,
    &finite_field_batch_inverse_skip_zeros::EXPORT,
    &finite_field_slice_add::EXPORT,
    &finite_field_slice_sub::EXPORT,
    &finite_field_slice_scale::EXPORT,
    &finite_field_slice_axpy::EXPORT,
    &finite_field_dot::EXPORT,
    &number_theory_mul_mod::EXPORT,
    &number_theory_mul_mod_u128::EXPORT,
    &number_theory_mod_pow_u128::EXPORT,
    &extended_gcd_u64::EXPORT,
    &modular_inverse_u64::EXPORT,
    &jacobi_i128::EXPORT,
    &legendre_u64::EXPORT,
    &primitive_root_u64::EXPORT,
    &finite_field_prime_generator::EXPORT,
    &finite_field_poly_mul_ntt::EXPORT,
    &finite_field_ntt_plan::EXPORT,
    &finite_field_ntt_transform::EXPORT,
    &uhp_to_disk::EXPORT,
    &disk_to_uhp::EXPORT,
    &poincare_disk_distance::EXPORT,
    &uhp_distance::EXPORT,
    &riemann_sphere_project::EXPORT,
    &complex_coordinate_to_polar::EXPORT,
    &complex_coordinate_from_polar::EXPORT,
    &complex_coordinate_abs_squared::EXPORT,
    &complex_coordinate_sign::EXPORT,
    &complex_coordinate_unwrap_phase::EXPORT,
    &complex_coordinate_continue_root::EXPORT,
    &riemann_sphere_unproject::EXPORT,
    &joukowsky_inverse::EXPORT,
    &conformal_disk_automorphism::EXPORT,
    &conformal_mobius_class::EXPORT,
    &conformal_fixed_points::EXPORT,
    &conformal_bergman_disk::EXPORT,
    &riemann_hyperelliptic::EXPORT,
    &riemann_superelliptic::EXPORT,
    &riemann_monodromy::EXPORT,
    &conformal_bergman_basis::EXPORT,
    &conformal_mobius_three_points::EXPORT,
    &conformal_schwarz_christoffel_evaluate::EXPORT,
    &conformal_schwarz_christoffel_solve::EXPORT,
    &conformal_szego_prepare_sampled::EXPORT,
    &conformal_szego_eval::EXPORT,
    &conformal_szego_eval_map::EXPORT,
    &conformal_szego_eval_bergman::EXPORT,
    &univalent_bieberbach::EXPORT,
    &univalent_fekete_szego::EXPORT,
    &univalent_hankel::EXPORT,
    &univalent_grunsky_matrix::EXPORT,
    &univalent_grunsky_check::EXPORT,
    &univalent_becker::EXPORT,
    &univalent_loewner_step::EXPORT,
];

include!("riemann_hilbert_codecs.rs");
include!("vector_codecs.rs");
include!("csr_structural_codecs.rs");
include!("mp_constant_codecs.rs");
include!("contour_codecs.rs");
include!("signal_filter_codecs.rs");
include!("complex_disk_codecs.rs");
include!("timeseries_codecs.rs");
include!("state_space_codecs.rs");
include!("complex_polynomial_codecs.rs");
include!("real_polynomial_codecs.rs");
include!("complex_bound_codecs.rs");
include!("branch_cut_codecs.rs");
include!("weierstrass_codecs.rs");
include!("carlson_complex_codecs.rs");
include!("continuation_codecs.rs");
include!("modular_codecs.rs");
include!("prime_group_codecs.rs");
include!("ntt_plan_codecs.rs");
include!("complex_dynamics_codecs.rs");
include!("special_probability_codecs.rs");
include!("glm_link_codecs.rs");
include!("survival_model_codecs.rs");
include!("special_complex_codecs.rs");
include!("mp_exp_codecs.rs");
include!("native_scalar_codecs.rs");
include!("mp_real_codecs.rs");
include!("mp_decimal_codecs.rs");
include!("mp_ball_codecs.rs");
include!("mp_enclosure_codecs.rs");
include!("mp_complex_codecs.rs");
include!("codec_generated.rs");
include!("mp_complex_trigonometry_codecs.rs");
include!("mp_complex_inverse_codecs.rs");
include!("mp_error_function_codecs.rs");
include!("mp_bessel_codecs.rs");
include!("mp_gamma_codecs.rs");
include!("mp_dirichlet_codecs.rs");
include!("mp_complex_root_power_codecs.rs");
include!("mp_complex_exp_log_codecs.rs");

include!("functional_codecs.rs");
include!("mp_exact_root_codecs.rs");
include!("object_codecs.rs");
include!("mp_logarithm_codecs.rs");
include!("dual_codecs.rs");
include!("mp_trigonometry_codecs.rs");
include!("geometry_codecs.rs");
include!("mp_hyperbolic_codecs.rs");
include!("native_analysis/mod.rs");
include!("mp_inverse_tangent_codecs.rs");
include!("mp_inverse_sine_codecs.rs");
include!("mp_power_codecs.rs");
include!("native_linear/codecs.rs");
include!("hyperbolic_codecs.rs");
include!("complex_coordinate_codecs.rs");
include!("riemann_surface_codecs.rs");
include!("conformal_mobius_codecs.rs");
include!("conformal_kernel_codecs.rs");
include!("conformal_szego_codecs.rs");
include!("univalent_codecs.rs");

impl Output for String {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::String(self))
    }
}
// Native void methods acknowledge completion with Bool; the ABI has no unit tuple.
impl Output for () {
    fn ty() -> TypeExpr {
        TypeExpr::Bool
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Bool(true))
    }
}
