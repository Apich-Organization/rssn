//! Portable `NumPy` bytes expressed through the existing string and tensor ABI.

use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::io::npy;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::RecordValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use num_traits::ToPrimitive;
use std::sync::Arc;

const MAX_VALUES: usize = 262_144;
const MAX_HEADER_BYTES: usize = 10 * 1024 + 12;
const MAX_ENCODED_BYTES: usize = 2 * (MAX_VALUES * size_of::<f64>() + MAX_HEADER_BYTES);

#[derive(Clone, Copy)]
enum Codec {
    Encode,
    Decode,
}

fn vector_type(elem: TypeExpr) -> TypeExpr {
    elem.tensor(NatExpr::Lit(1))
}

fn vector_parts(
    value: &Value,
    limit: usize,
) -> Result<&[Value], OpError> {
    if let Value::Tensor(tensor) = value
        && tensor.shape == [tensor.data.len()]
        && tensor.data.len() <= limit
    {
        Ok(&tensor.data)
    } else {
        Err(OpError::Domain("invalid NPY vector shape or size".into()))
    }
}

fn hex_digit(nibble: u8) -> char {
    match nibble {
        | 0..=9 => char::from(b'0'.saturating_add(nibble)),
        | _ => char::from(b'a'.saturating_add(nibble.saturating_sub(10))),
    }
}

fn encode(args: &[Value]) -> Result<Value, OpError> {
    let [shape, data] = args else {
        return Err(OpError::Arity {
            expected: ArityInfo::fixed(2),
            got: args.len(),
        });
    };
    let shape = vector_parts(shape, MAX_HEADER_BYTES / 2)?
        .iter()
        .map(|extent| {
            match extent {
                | Value::Int(n) => {
                    n.to_usize().ok_or_else(|| {
                        OpError::Domain(
                            "NPY dimensions must fit nonnegative platform indices".into(),
                        )
                    })
                },
                | _ => Err(OpError::Domain("NPY dimensions must be integers".into())),
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let values = vector_parts(data, MAX_VALUES)?
        .iter()
        .map(|value| {
            match value {
                | Value::Float(value) => Ok(*value),
                | _ => {
                    Err(OpError::Domain(
                        "NPY storage must be dimensionless Float".into(),
                    ))
                },
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bytes = npy::encode(&shape, &values).map_err(OpError::Domain)?;
    let mut hex = String::with_capacity(bytes.len().saturating_mul(2));
    for byte in bytes {
        hex.push(hex_digit(byte >> 4));
        hex.push(hex_digit(byte & 0x0f));
    }
    Ok(Value::String(hex))
}

fn decode(args: &[Value]) -> Result<Value, OpError> {
    let [source] = args else {
        return Err(OpError::Arity {
            expected: ArityInfo::fixed(1),
            got: args.len(),
        });
    };
    let Value::String(source) = source else {
        return Err(OpError::Domain(
            "NPY data must be a hexadecimal string".into(),
        ));
    };
    let (pairs, remainder) = source.as_bytes().as_chunks::<2>();
    if source.len() > MAX_ENCODED_BYTES || !remainder.is_empty() {
        return Err(OpError::Domain("invalid hexadecimal NPY length".into()));
    }
    let bytes = pairs
        .iter()
        .map(|&[high, low]| {
            let invalid = || OpError::Domain("NPY data contains a nonhexadecimal digit".into());
            let high = char::from(high)
                .to_digit(16)
                .and_then(|n| u8::try_from(n).ok())
                .ok_or_else(invalid)?;
            let low = char::from(low)
                .to_digit(16)
                .and_then(|n| u8::try_from(n).ok())
                .ok_or_else(invalid)?;
            Ok((high << 4) | low)
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    let (shape, values) = npy::decode_with_limit(&bytes, MAX_VALUES).map_err(OpError::Domain)?;
    Ok(Value::Record(RecordValue {
        fields: vec![
            (
                "shape".into(),
                Value::Tensor(TensorValue {
                    shape: vec![shape.len()],
                    elem: TypeExpr::Int,
                    data: shape.into_iter().map(Value::int).collect(),
                }),
            ),
            (
                "data".into(),
                Value::Tensor(TensorValue {
                    shape: vec![values.len()],
                    elem: TypeExpr::Float,
                    data: values.into_iter().map(Value::Float).collect(),
                }),
            ),
        ],
    }))
}

impl EvalDyn for Codec {
    fn signature(&self) -> Signature {
        let shape = vector_type(TypeExpr::Int);
        let data = vector_type(TypeExpr::Float);
        match self {
            | Self::Encode => Signature::new(vec![shape, data], TypeExpr::String),
            | Self::Decode => {
                Signature::new(
                    vec![TypeExpr::String],
                    TypeExpr::Record(vec![("shape".into(), shape), ("data".into(), data)]),
                )
            },
        }
    }

    fn eval(
        &self,
        args: &[Value],
        _ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if !cfg!(feature = "npy") {
            return Err(OpError::Unsupported {
                backend: "rssn".into(),
                feature: "npy".into(),
            });
        }
        match self {
            | Self::Encode => encode(args),
            | Self::Decode => decode(args),
        }
    }
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, path, doc, codec) in [
        (
            "num_npy_encode",
            2,
            "crate::io::npy::encode",
            "Encode shape and dimensionless row-major Float storage as little-endian NumPy float64 bytes represented in lowercase hexadecimal. Accepts scalars and empty arrays, up to 262144 values.",
            Codec::Encode,
        ),
        (
            "num_npy_decode",
            1,
            "crate::io::npy::decode",
            "Decode hexadecimal NumPy float32 or float64 bytes into shape and row-major Float storage. Accepts either byte order and C/Fortran layouts, up to 262144 values.",
            Codec::Decode,
        ),
    ] {
        let op = installer.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .rust_path(path)
            .category("numerical_matrix")
            .doc(doc)
            .signature(codec.signature())
            .eval_dyn(Arc::new(codec)))?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
