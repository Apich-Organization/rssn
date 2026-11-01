//! The dynamic value type.
//!
//! [`Value`] is what crosses a dynamic call boundary: into
//! [`Session::call_value`](crate::api::Session::call_value), out of a
//! [`Plugin`](super::Plugin), through the C ABI, and over a process boundary
//! as deterministic CBOR.
//!
//! It is a *closed* type: a host can construct every value it will ever see,
//! and rssn can encode every one of them. That is what makes it usable as an
//! ABI. Anything rssn does not model is carried as
//! [`Value::Opaque`] — a kind tag and bytes the kernel never interprets —
//! so extension never needs a new variant.
//!
//! # Relation to terms
//!
//! Terms are the engine's representation; `Value` is the boundary's. The
//! mapping is deliberately partial in one direction and total in the other:
//!
//! * every [`Term`](crate::api::Term) has a `Value`
//!   ([`Term::to_value`](crate::api::Term::to_value)); a term that is not a
//!   number, a quantity or a tensor becomes [`Value::Term`], a handle to the
//!   term itself;
//! * every `Value` that mentions no handle becomes a term
//!   ([`Session::term_of`](crate::api::Session::term_of)). A `Value::Term`
//!   from *another* session is refused, because a node id means nothing
//!   outside the graph that minted it.

use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use serde::Deserialize;
use serde::Serialize;

use crate::graph::Number;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::TypeExpr;

/// Identifies one [`Session`](crate::api::Session) for the lifetime of the
/// process.
///
/// Handles minted by one session are refused by another, which turns a
/// use-after-free of a term handle into an error instead of a wrong answer.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct SessionId(pub u64);

impl fmt::Display for SessionId {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "session#{}", self.0)
    }
}

/// An owned handle to a term of a session.
///
/// Unlike [`Term`](crate::api::Term), which borrows the session, a `TermRef`
/// is `Copy + Send + Sync` and owns nothing: it is a `(session, node)` pair.
/// It can be stored in a closure, sent to another thread and turned back
/// into a `Term` with
/// [`Session::term_ref`](crate::api::Session::term_ref) — on the session
/// that minted it.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct TermRef {
    /// The session that owns the node.
    pub session: SessionId,
    /// The node, as returned by [`Term::id`](crate::api::Term::id).
    pub node: u32,
}

/// An exact or floating complex number.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ComplexValue {
    /// Both parts exact.
    Exact {
        /// Real part.
        re: Number,
        /// Imaginary part.
        im: Number,
    },
    /// Both parts double precision.
    Float {
        /// Real part.
        re: f64,
        /// Imaginary part.
        im: f64,
    },
}

impl ComplexValue {
    /// The value as a double-precision pair, rounding exact parts.
    #[must_use]
    pub fn to_f64(&self) -> (f64, f64) {
        match self {
            | Self::Exact { re, im } => (re.to_f64(), im.to_f64()),
            | Self::Float { re, im } => (*re, *im),
        }
    }

    /// Whether both parts are exact.
    #[must_use]
    pub const fn is_exact(&self) -> bool {
        matches!(self, Self::Exact { .. })
    }
}

/// A unit definition carried by value: enough to print, convert and check a
/// quantity without a registry lookup.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnitValue {
    /// Symbolic name, as written in a term (`m`, `km`, `degC`).
    pub name: String,
    /// How it is printed (`m`, `°C`).
    pub display: String,
    /// Its dimension.
    pub dimension: DimExpr,
    /// Scale to coherent SI units, exact: `SI = value * scale + offset`.
    pub scale: ScaleValue,
    /// Additive offset for affine units, exact.
    pub offset: Option<ScaleValue>,
    /// Whether the unit is affine (a point scale such as `degC`).
    pub affine: bool,
}

/// An exact scale factor `rat * pi^pi`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScaleValue {
    /// Numerator and denominator of the rational factor.
    pub numer: BigInt,
    /// Denominator of the rational factor, always positive.
    pub denom: BigInt,
    /// Power of pi in the factor.
    pub pi: i32,
}

impl ScaleValue {
    /// The scale `1`.
    #[must_use]
    pub fn one() -> Self {
        Self {
            numer: BigInt::from(1),
            denom: BigInt::from(1),
            pi: 0,
        }
    }

    /// The exact rational value of the scale, ignoring any factor of pi.
    #[must_use]
    pub fn rational(&self) -> BigRational {
        BigRational::new(self.numer.clone(), self.denom.clone())
    }
}

impl From<BigRational> for ScaleValue {
    fn from(value: BigRational) -> Self {
        let (numer, denom) = value.into();
        Self {
            numer,
            denom,
            pi: 0,
        }
    }
}

/// A value with a unit: `3 m`, `20 degC`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct QuantityValue {
    /// The numeric part.
    pub value: Box<Value>,
    /// The unit.
    pub unit: UnitValue,
}

/// A dense tensor in row-major order.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TensorValue {
    /// Extents, outermost first. Empty for a scalar.
    pub shape: Vec<usize>,
    /// Elements, row-major. `data.len()` is the product of `shape`.
    pub data: Vec<Value>,
}

impl TensorValue {
    /// A rank-0 tensor holding `value`.
    #[must_use]
    pub fn scalar(value: Value) -> Self {
        Self {
            shape: Vec::new(),
            data: vec![value],
        }
    }

    /// A rank-1 tensor.
    #[must_use]
    pub fn vector(data: Vec<Value>) -> Self {
        let len = data.len();
        Self {
            shape: vec![len],
            data,
        }
    }

    /// The rank.
    #[must_use]
    pub const fn rank(&self) -> usize {
        self.shape.len()
    }

    /// The number of elements the shape implies.
    #[must_use]
    pub fn extent(&self) -> usize {
        self.shape.iter().product()
    }

    /// Whether `data` matches `shape`.
    #[must_use]
    pub fn is_consistent(&self) -> bool {
        self.extent() == self.data.len()
    }
}

/// How a sparse tensor stores its entries.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SparseFormat {
    /// Compressed sparse row.
    Csr,
    /// Compressed sparse column.
    Csc,
    /// Coordinate list.
    Coo,
}

/// A sparse tensor.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SparseValue {
    /// Extents, outermost first.
    pub shape: Vec<usize>,
    /// Storage format.
    pub format: SparseFormat,
    /// Flat row-major indices of the stored entries, ascending.
    pub indices: Vec<u64>,
    /// The stored entries, in the order of `indices`.
    pub data: Vec<Value>,
}

/// A record: named fields in declaration order.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RecordValue {
    /// Field names and values.
    pub fields: Vec<(String, Value)>,
}

impl RecordValue {
    /// The value of a field.
    #[must_use]
    pub fn get(
        &self,
        name: &str,
    ) -> Option<&Value> {
        self.fields
            .iter()
            .find(|(field, _)| field == name)
            .map(|(_, value)| value)
    }
}

/// A tagged union value.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct VariantValue {
    /// The tag.
    pub tag: String,
    /// The payload.
    pub payload: Vec<Value>,
}

/// A callback the host supplies for a higher-order call.
///
/// The index refers to the callback table of the
/// [`CallCtx`](super::CallCtx) of the call in progress; it is not a global
/// identifier.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct CallbackRef(pub u32);

/// A value whose dimension is only known at run time.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DynValue {
    /// The value.
    pub value: Box<Value>,
    /// Its dimension.
    pub dim: DimExpr,
}

/// A value rssn does not interpret: a host type, or a plugin's own data.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpaqueValue {
    /// What the bytes mean, by convention `vendor/type`.
    pub kind: String,
    /// The payload, opaque to the kernel.
    pub bytes: Vec<u8>,
}

/// A dynamically typed value.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Value {
    /// An exact integer.
    Int(BigInt),
    /// An exact rational, or a symbolic real when it came from a term.
    Real(Number),
    /// A double-precision float.
    Float(f64),
    /// A complex number.
    Complex(ComplexValue),
    /// A truth value.
    Bool(bool),
    /// Text.
    String(String),
    /// A free symbol of the expression graph.
    ///
    /// Distinct from [`Value::String`] because the two are not interchangeable
    /// at the boundary: a string is data an operator may concatenate, a symbol
    /// is a hole in an expression that an operator may differentiate against.
    Symbol(String),
    /// A unit, not a quantity.
    Unit(UnitValue),
    /// A value with a unit.
    Quantity(QuantityValue),
    /// A dense tensor.
    Tensor(TensorValue),
    /// A sparse tensor.
    Sparse(SparseValue),
    /// A record.
    Record(RecordValue),
    /// A tagged union.
    Variant(VariantValue),
    /// A handle to a term of the calling session.
    Term(TermRef),
    /// A host callback, for higher-order operators.
    Function(CallbackRef),
    /// A value with a run-time dimension.
    Dyn(DynValue),
    /// Host data the kernel does not interpret.
    Opaque(OpaqueValue),
}

impl Value {
    /// The value `0`.
    #[must_use]
    pub fn zero() -> Self {
        Self::Int(BigInt::from(0))
    }

    /// An exact integer.
    #[must_use]
    pub fn int(
        value: impl Into<BigInt>,
    ) -> Self {
        Self::Int(value.into())
    }

    /// An exact rational.
    #[must_use]
    pub fn rational(value: BigRational) -> Self {
        if value.is_integer() {
            Self::Int(value.to_integer())
        } else {
            Self::Real(Number::rat(value))
        }
    }

    /// The value as an exact number, when it is one.
    #[must_use]
    pub fn as_number(&self) -> Option<Number> {
        match self {
            | Self::Int(n) => Some(Number::Int(n.clone())),
            | Self::Real(n) => Some(n.clone()),
            | Self::Float(f) => Some(Number::from(*f)),
            | _ => None,
        }
    }

    /// The value as an exact rational, when it is one.
    #[must_use]
    pub fn as_rational(&self) -> Option<BigRational> {
        self.as_number().and_then(|n| n.to_rational())
    }

    /// The value as a float, when it has one.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            | Self::Int(n) => n.to_f64(),
            | Self::Real(n) => Some(n.to_f64()),
            | Self::Float(f) => Some(*f),
            | _ => None,
        }
    }

    /// Whether the value is a scalar number.
    #[must_use]
    pub const fn is_number(&self) -> bool {
        matches!(self, Self::Int(_) | Self::Real(_) | Self::Float(_))
    }

    /// Whether the value is exactly zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match self {
            | Self::Int(n) => n.sign() == num_bigint::Sign::NoSign,
            | Self::Real(n) => n.is_zero(),
            | Self::Float(f) => *f == 0.0,
            | _ => false,
        }
    }

    /// The type of the value, for error reporting.
    #[must_use]
    pub fn type_expr(&self) -> TypeExpr {
        match self {
            | Self::Int(_) => TypeExpr::Int,
            | Self::Real(_) => TypeExpr::Real,
            | Self::Float(_) => TypeExpr::Float,
            | Self::Complex(_) => TypeExpr::Complex,
            | Self::Bool(_) => TypeExpr::Bool,
            | Self::String(_) => TypeExpr::String,
            | Self::Symbol(_) => TypeExpr::Var("Symbol".to_owned()),
            | Self::Unit(u) => TypeExpr::Unit(u.dimension.clone()),
            | Self::Quantity(q) => TypeExpr::Dim(
                Box::new(TypeExpr::Real),
                q.unit.dimension.clone(),
            ),
            | Self::Tensor(t) => TypeExpr::Tensor {
                rank: NatExpr::Lit(u64::try_from(t.rank()).unwrap_or(u64::MAX)),
                shape: Some(
                    t.shape
                        .iter()
                        .map(|&n| NatExpr::Lit(u64::try_from(n).unwrap_or(u64::MAX)))
                        .collect(),
                ),
                elem: Box::new(t.data.first().map_or(TypeExpr::Real, Self::type_expr)),
            },
            | Self::Sparse(s) => TypeExpr::Sparse {
                format: match s.format {
                    | SparseFormat::Csr => "csr".to_owned(),
                    | SparseFormat::Csc => "csc".to_owned(),
                    | SparseFormat::Coo => "coo".to_owned(),
                },
                rank: NatExpr::Lit(u64::try_from(s.shape.len()).unwrap_or(u64::MAX)),
                elem: Box::new(s.data.first().map_or(TypeExpr::Real, Self::type_expr)),
            },
            | Self::Record(r) => TypeExpr::Record(
                r.fields
                    .iter()
                    .map(|(name, value)| (name.clone(), value.type_expr()))
                    .collect(),
            ),
            | Self::Variant(v) => TypeExpr::Variant(vec![(
                v.tag.clone(),
                v.payload.iter().map(Self::type_expr).collect(),
            )]),
            | Self::Term(_) => TypeExpr::Dyn(crate::kernel::signature::RepExpr {
                repr: "term".to_owned(),
                dim: None,
            }),
            | Self::Function(_) => TypeExpr::Function(Vec::new(), Box::new(TypeExpr::Never)),
            | Self::Dyn(d) => TypeExpr::Dyn(crate::kernel::signature::RepExpr {
                repr: "dyn".to_owned(),
                dim: Some(d.dim.clone()),
            }),
            | Self::Opaque(o) => TypeExpr::Opaque(o.kind.clone()),
        }
    }

    /// The dimension this value is known to carry, when it is a quantity.
    #[must_use]
    pub const fn dimension(&self) -> Option<&DimExpr> {
        match self {
            | Self::Quantity(q) => Some(&q.unit.dimension),
            | Self::Unit(u) => Some(&u.dimension),
            | Self::Dyn(d) => Some(&d.dim),
            | _ => None,
        }
    }

    /// The encoding of this value as deterministic CBOR.
    ///
    /// # Errors
    /// See [`encode_cbor`](super::codec::encode_cbor).
    pub fn to_cbor(&self) -> Result<Vec<u8>, super::codec::CodecError> {
        super::codec::encode_cbor(self)
    }

    /// Reads a value written by [`Value::to_cbor`].
    ///
    /// # Errors
    /// See [`decode_cbor`](super::codec::decode_cbor).
    pub fn from_cbor(bytes: &[u8]) -> Result<Self, super::codec::CodecError> {
        super::codec::decode_cbor(bytes)
    }
}

impl fmt::Display for Value {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Int(n) => write!(f, "{n}"),
            | Self::Real(n) => write!(f, "{n}"),
            | Self::Float(x) => write!(f, "{x}"),
            | Self::Complex(c) => {
                let (re, im) = c.to_f64();
                write!(f, "{re}{}{im}i", if im < 0.0 { "-" } else { "+" })
            },
            | Self::Bool(b) => write!(f, "{b}"),
            | Self::String(s) => write!(f, "{s:?}"),
            | Self::Symbol(s) => f.write_str(s),
            | Self::Unit(u) => f.write_str(&u.display),
            | Self::Quantity(q) => write!(f, "{} {}", q.value, q.unit.display),
            | Self::Tensor(t) => write!(f, "<tensor {:?}>", t.shape),
            | Self::Sparse(s) => write!(f, "<sparse {:?} {} entries>", s.shape, s.data.len()),
            | Self::Record(r) => {
                f.write_str("{")?;
                for (i, (name, value)) in r.fields.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    write!(f, "{name}: {value}")?;
                }
                f.write_str("}")
            },
            | Self::Variant(v) => write!(f, "{} {:?}", v.tag, v.payload.len()),
            | Self::Term(t) => write!(f, "term#{}", t.node),
            | Self::Function(c) => write!(f, "<callback {}>", c.0),
            | Self::Dyn(d) => write!(f, "{}", d.value),
            | Self::Opaque(o) => write!(f, "<{}: {} bytes>", o.kind, o.bytes.len()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_round_trip() {
        let values = [
            Value::int(7_i32),
            Value::rational(BigRational::new(BigInt::from(1), BigInt::from(3))),
            Value::Float(2.5),
            Value::Bool(true),
            Value::String("x".to_owned()),
        ];
        for value in values {
            let bytes = value.to_cbor().unwrap_or_else(|e| panic!("{e}"));
            let back = Value::from_cbor(&bytes).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(back, value);
        }
    }

    #[test]
    fn exact_arithmetic_survives_the_boundary() {
        let third = Value::rational(BigRational::new(BigInt::from(1), BigInt::from(3)));
        assert!(third.as_number().is_some_and(|n| n.is_exact()));
        assert_eq!(Value::int(6_i32).as_f64(), Some(6.0));
        assert!(Value::int(0_i32).is_zero());
    }

    #[test]
    fn tensors_report_shape_and_type() {
        let t = TensorValue {
            shape: vec![2, 2],
            data: vec![Value::int(1), Value::int(2), Value::int(3), Value::int(4)],
        };
        assert!(t.is_consistent());
        assert_eq!(t.rank(), 2);
        assert_eq!(t.extent(), 4);
        let TypeExpr::Tensor { rank, shape, .. } = Value::Tensor(t).type_expr() else {
            panic!("expected a tensor type");
        };
        assert_eq!(rank, NatExpr::Lit(2));
        assert_eq!(shape.map(|s| s.len()), Some(2));
    }

    #[test]
    fn opaque_values_carry_their_kind() {
        let value = Value::Opaque(OpaqueValue {
            kind: "acme/poly".to_owned(),
            bytes: vec![1, 2, 3],
        });
        let bytes = value.to_cbor().unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(Value::from_cbor(&bytes), Ok(value));
    }
}
