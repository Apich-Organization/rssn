//! Validated portable algebraic models with explicit Float and exact Real specializations.
//! Every model carries defining values, a kind/version, and no process-local pointers.
use super::dimensions::DimContract;
use super::dimensions::DimFormula;
use super::dimensions::FieldDimension;
use super::dimensions::InputBinding;
use super::dimensions::InputConstraint;
use super::dimensions::OutputBinding;
use super::dimensions::{
    self,
};
use super::numeric;
use super::runtime_dimensions;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::RecordValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernel::value::VariantValue;
use crate::kernels::matrix::Backend;
use crate::kernels::matrix::FaerDecompositionResult;
use crate::kernels::matrix::FaerDecompositionType;
use crate::kernels::matrix::Field;
use crate::kernels::matrix::Matrix;
use crate::kernels::polynomial::Polynomial;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use num_traits::Zero;
use std::sync::Arc;

mod factors;
mod multilinear;

const LIMIT: usize = 16_384;
const SIDE: usize = 128;
fn error(message: impl std::fmt::Display) -> OpError {
    OpError::Domain(message.to_string())
}
fn field<'a>(
    value: &'a Value,
    name: &str,
    names: &[&str],
) -> Result<&'a Value, OpError> {
    let Value::Record(record) = value else {
        return Err(error("expected a portable model record"));
    };
    if record.fields.len() != names.len()
        || names
            .iter()
            .any(|name| record.fields.iter().filter(|(key, _)| key == name).count() != 1)
    {
        return Err(error(
            "portable model has missing, duplicate, or undeclared fields",
        ));
    }
    record
        .fields
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value)
        .ok_or_else(|| error("missing model field"))
}
fn kind(
    value: &Value,
    expected: &str,
) -> Result<(), OpError> {
    if matches!(value, Value::String(actual) if actual==expected) {
        Ok(())
    } else {
        Err(error("unknown portable model kind/version"))
    }
}
fn integer(value: &Value) -> Result<usize, OpError> {
    match value {
        | Value::Int(n) => {
            n.to_usize()
                .ok_or_else(|| error("expected a nonnegative machine integer"))
        },
        | _ => Err(error("expected an integer")),
    }
}
fn tensor_type(
    rank: u64,
    elem: TypeExpr,
) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(elem),
    }
}
fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    })
}
fn option_type(payload: TypeExpr) -> TypeExpr {
    TypeExpr::Variant(vec![
        ("None".into(), vec![]),
        ("Some".into(), vec![payload]),
    ])
}
fn option(value: Option<Value>) -> Value {
    Value::Variant(VariantValue {
        tag: if value.is_some() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: value.into_iter().collect(),
    })
}
trait Scalar: Field + ToPrimitive {
    const EXACT: bool;
    fn ty() -> TypeExpr;
    fn decode(value: &Value) -> Result<Self, OpError>;
    fn encode(self) -> Result<Value, OpError>;
}
impl Scalar for f64 {
    const EXACT: bool = false;

    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn decode(value: &Value) -> Result<Self, OpError> {
        numeric(value, 0)
    }

    fn encode(self) -> Result<Value, OpError> {
        super::output_scalar(self)
    }
}
impl Scalar for BigRational {
    const EXACT: bool = true;

    fn ty() -> TypeExpr {
        TypeExpr::Real
    }

    fn decode(value: &Value) -> Result<Self, OpError> {
        let result = match value {
            | Value::Int(n) | Value::Real(crate::graph::Number::Int(n)) => {
                Self::from_integer(n.clone())
            },
            | Value::Real(crate::graph::Number::Rat(n)) => n.clone(),
            | Value::Dyn(d) => return Self::decode(&d.value),
            | Value::Quantity(q) => {
                if q.unit.affine
                    || q.unit.offset.is_some()
                    || q.unit.scale.pi != 0
                    || q.unit.scale.denom.is_zero()
                {
                    return Err(error("exact matrices require rational linear unit scales"));
                }
                Self::decode(&q.value)? * q.unit.scale.rational()
            },
            | _ => {
                return Err(error(
                    "exact matrix entries require Int or Real, never Float",
                ));
            },
        };
        if result.numer().bits() > 4096 || result.denom().bits() > 4096 {
            return Err(error("exact matrix entry exceeds 4096-bit budget"));
        }
        Ok(result)
    }

    fn encode(self) -> Result<Value, OpError> {
        if self.numer().bits() > 4096 || self.denom().bits() > 4096 {
            return Err(error("exact matrix result exceeds 4096-bit budget"));
        }
        Ok(Value::Real(crate::graph::Number::rat(self)))
    }
}
fn matrix_type(exact: bool) -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        ("backend".into(), TypeExpr::Int),
        (
            "entries".into(),
            tensor_type(
                2,
                if exact {
                    TypeExpr::Real
                } else {
                    TypeExpr::Float
                },
            ),
        ),
    ])
}
fn array<T: Scalar>(
    value: &Value,
    rank: usize,
) -> Result<(Vec<usize>, Vec<T>), OpError> {
    let Value::Tensor(t) = value else {
        return Err(error("expected a dense tensor"));
    };
    if t.shape.len() != rank
        || t.shape.iter().any(|n| *n > SIDE)
        || t.shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(t.data.len())
        || t.data.len() > LIMIT
    {
        return Err(error("invalid model tensor shape or allocation budget"));
    }
    Ok((
        t.shape.clone(),
        t.data.iter().map(T::decode).collect::<Result<_, _>>()?,
    ))
}
fn encode_array<T: Scalar>(
    shape: Vec<usize>,
    data: Vec<T>,
) -> Result<Value, OpError> {
    if data.len() > LIMIT
        || shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(data.len())
    {
        return Err(error("result exceeds portable tensor budget"));
    }
    Ok(Value::Tensor(TensorValue {
        shape,
        elem: T::ty(),
        data: data.into_iter().map(T::encode).collect::<Result<_, _>>()?,
    }))
}
fn backend(
    value: &Value,
    exact: bool,
) -> Result<Backend, OpError> {
    match integer(value)? {
        | 0 => Ok(Backend::Native),
        | 1 if !exact => Ok(Backend::Faer),
        | _ => Err(error("backend must be 0=Native or 1=Faer (Float only)")),
    }
}
fn decode_matrix<T: Scalar>(value: &Value) -> Result<Matrix<T>, OpError> {
    let names = ["format", "backend", "entries"];
    kind(
        field(value, "format", &names)?,
        if T::EXACT {
            "rssn/matrix-real/v1"
        } else {
            "rssn/matrix-float/v1"
        },
    )?;
    let backend = backend(field(value, "backend", &names)?, T::EXACT)?;
    let (shape, data) = array::<T>(field(value, "entries", &names)?, 2)?;
    if T::EXACT && shape.iter().any(|n| *n > 32) {
        return Err(error("exact matrices are limited to 32 rows and columns"));
    }
    Ok(Matrix::new(shape[0], shape[1], data).with_backend(backend))
}
fn encode_matrix<T: Scalar>(matrix: Matrix<T>) -> Result<Value, OpError> {
    let backend = usize::from(matrix.backend == Backend::Faer);
    let shape = vec![matrix.rows(), matrix.cols()];
    Ok(record(vec![
        (
            "format",
            Value::String(
                if T::EXACT {
                    "rssn/matrix-real/v1"
                } else {
                    "rssn/matrix-float/v1"
                }
                .into(),
            ),
        ),
        ("backend", Value::int(backend)),
        ("entries", encode_array(shape, matrix.into_data())?),
    ]))
}
fn contract(
    inputs: Vec<InputBinding>,
    constraints: Vec<InputConstraint>,
    paths: Vec<OutputBinding>,
    result: Option<&'static DimFormula>,
) -> DimContract {
    static CONTRACTS: std::sync::OnceLock<
        std::sync::Mutex<std::collections::BTreeMap<String, DimContract>>,
    > = std::sync::OnceLock::new();
    let key = format!("{inputs:?}|{constraints:?}|{paths:?}|{result:?}");
    let mut contracts = CONTRACTS
        .get_or_init(Default::default)
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    *contracts.entry(key).or_insert_with(|| {
        DimContract {
            inputs: Box::leak(inputs.into_boxed_slice()),
            constraints: Box::leak(constraints.into_boxed_slice()),
            outputs: &[],
            output_paths: Box::leak(paths.into_boxed_slice()),
            result,
            other_args_dimensionless: true,
        }
    })
}
const D: DimFormula = DimFormula::Var("D");
const E: DimFormula = DimFormula::Var("E");
const INVERSE: DimFormula = DimFormula::Pow(&D, -1, 1);
const PRODUCT: DimFormula = DimFormula::Mul(&[D, E]);
const fn input(
    name: &'static str,
    argument: usize,
    path: &'static [&'static str],
) -> InputBinding {
    InputBinding { name, argument, path }
}
const fn constraint(
    argument: usize,
    path: &'static [&'static str],
    formula: &'static DimFormula,
) -> InputConstraint {
    InputConstraint {
        argument,
        path,
        formula,
    }
}
const fn output(
    path: &'static [&'static str],
    formula: &'static DimFormula,
) -> OutputBinding {
    OutputBinding { path, formula }
}
#[derive(Clone, Copy)]
enum MatrixMethod {
    New,
    Zeros,
    Identity,
    WithBackend,
    SetBackend,
    Get,
    Set,
    Rows,
    Cols,
    Data,
    IntoData,
    GetCols,
    Rref,
    Transpose,
    Strassen,
    Det,
    DetBlock,
    DetLu,
    Lu,
    Inverse,
    NullSpace,
    Rank,
    Trace,
    Symmetric,
    Diagonal,
    Frobenius,
    L1,
    LInf,
    IsIdentity,
    Orthogonal,
    Jacobi,
    DecomposeLu,
    DecomposeQr,
    DecomposeCholesky,
    DecomposeSvd,
    DecomposeEigen,
}
struct MatrixAdapter {
    method: MatrixMethod,
    exact: bool,
    contract: DimContract,
    signature: Signature,
}
#[allow(clippy::too_many_lines)] // Reviewed signature and dimension cases remain adjacent.
fn matrix_adapter(
    method: MatrixMethod,
    exact: bool,
) -> MatrixAdapter {
    use MatrixMethod as M;
    let scalar = if exact {
        TypeExpr::Real
    } else {
        TypeExpr::Float
    };
    let model = matrix_type(exact);
    let mut args = vec![model.clone()];
    let mut ret = model.clone();
    let mut inputs = vec![input("D", 0, &["entries"])];
    let mut checks = vec![];
    let mut outputs = vec![output(&["entries"], &D)];
    let mut result = None;
    match method {
        | M::New => {
            args = vec![tensor_type(2, scalar), TypeExpr::Int];
            inputs = vec![input("D", 0, &[])];
        },
        | M::Zeros => {
            args = vec![TypeExpr::Int, TypeExpr::Int, scalar];
            inputs = vec![input("D", 2, &[])];
        },
        | M::Identity => {
            args = vec![TypeExpr::Int];
            inputs.clear();
            outputs.clear();
        },
        | M::WithBackend | M::SetBackend => args.push(TypeExpr::Int),
        | M::Get => {
            args.extend([TypeExpr::Int, TypeExpr::Int]);
            ret = scalar;
            outputs.clear();
            result = Some(&D);
        },
        | M::Set => {
            args.extend([TypeExpr::Int, TypeExpr::Int, scalar]);
            checks.push(constraint(3, &[], &D));
        },
        | M::Rows | M::Cols | M::Rank => {
            ret = TypeExpr::Int;
            outputs.clear();
        },
        | M::Data | M::IntoData => {
            ret = tensor_type(1, scalar);
            outputs.clear();
            result = Some(&D);
        },
        | M::GetCols => {
            ret = tensor_type(2, scalar);
            outputs.clear();
            result = Some(&D);
        },
        | M::Rref => {
            ret = TypeExpr::Record(vec![
                ("matrix".into(), model),
                ("rank".into(), TypeExpr::Int),
            ]);
            outputs.clear();
        },
        | M::Strassen => {
            args.push(model);
            inputs.push(input("E", 1, &["entries"]));
            outputs = vec![output(&["entries"], &PRODUCT)];
        },
        | M::Det | M::DetBlock | M::DetLu => {
            ret = TypeExpr::Record(vec![
                ("normalized".into(), scalar.clone()),
                ("entry_unit".into(), scalar),
                ("order".into(), TypeExpr::Int),
            ]);
            outputs = vec![output(&["entry_unit"], &D)];
        },
        | M::Lu => {
            ret = TypeExpr::Record(vec![
                ("packed".into(), model),
                ("swaps".into(), TypeExpr::Int),
                ("entry_unit".into(), scalar),
            ]);
            outputs = vec![output(&["entry_unit"], &D)];
        },
        | M::Inverse => {
            ret = option_type(model);
            outputs = vec![output(&["entries"], &INVERSE)];
        },
        | M::NullSpace => outputs.clear(),
        | M::Trace => {
            ret = scalar;
            outputs.clear();
            result = Some(&D);
        },
        | M::Frobenius | M::L1 | M::LInf => {
            ret = TypeExpr::Float;
            outputs.clear();
            result = Some(&D);
        },
        | M::Symmetric | M::Diagonal => {
            ret = TypeExpr::Bool;
            outputs.clear();
        },
        | M::IsIdentity | M::Orthogonal => {
            inputs.clear();
            args.push(TypeExpr::Float);
            ret = TypeExpr::Bool;
            outputs.clear();
        },
        | M::Jacobi => {
            args.extend([TypeExpr::Int, TypeExpr::Float]);
            checks.push(constraint(2, &[], &D));
            ret = TypeExpr::Record(vec![
                ("values".into(), tensor_type(1, TypeExpr::Float)),
                ("vectors".into(), model),
            ]);
            outputs = vec![output(&["values"], &D)];
        },
        | M::DecomposeLu => {
            ret = TypeExpr::Record(vec![
                ("l".into(), model.clone()),
                ("u".into(), model),
                ("p".into(), tensor_type(1, TypeExpr::Int)),
            ]);
            outputs = vec![output(&["u", "entries"], &D)];
        },
        | M::DecomposeQr => {
            ret = TypeExpr::Record(vec![("q".into(), model.clone()), ("r".into(), model)]);
            outputs = vec![output(&["r", "entries"], &D)];
        },
        | M::DecomposeCholesky => {
            ret = TypeExpr::Record(vec![("l".into(), model)]);
            outputs = vec![output(&["l", "entries"], &DimFormula::Pow(&D, 1, 2))];
        },
        | M::DecomposeSvd => {
            ret = TypeExpr::Record(vec![
                ("u".into(), model.clone()),
                ("s".into(), tensor_type(1, TypeExpr::Float)),
                ("v".into(), model),
            ]);
            outputs = vec![output(&["s"], &D)];
        },
        | M::DecomposeEigen => {
            ret = TypeExpr::Record(vec![
                ("values".into(), tensor_type(1, TypeExpr::Float)),
                ("vectors".into(), model),
            ]);
            outputs = vec![output(&["values"], &D)];
        },
        | M::Transpose => {},
    }
    let mut contract = contract(inputs, checks, outputs, result);
    if matches!(method, M::Inverse) {
        contract.output_paths = &[];
        contract.outputs = &[
            ("format", FieldDimension::NonNumeric),
            ("backend", FieldDimension::NonNumeric),
            ("entries", FieldDimension::Formula(&INVERSE)),
        ];
    }
    let signature = dimensions::apply_contract(Signature::new(args, ret), &contract);
    MatrixAdapter {
        method,
        exact,
        contract,
        signature,
    }
}
impl MatrixAdapter {
    #[allow(clippy::too_many_lines)] // Native dispatch shared by both scalar specializations.
    fn run<T: Scalar>(
        &self,
        args: &[Value],
    ) -> Result<Value, OpError> {
        use MatrixMethod as M;
        match self.method {
            | M::New => {
                let (shape, data) = array::<T>(&args[0], 2)?;
                if T::EXACT && shape.iter().any(|n| *n > 32) {
                    return Err(error("exact matrices are limited to 32 rows and columns"));
                }
                return encode_matrix(
                    Matrix::new(shape[0], shape[1], data)
                        .with_backend(backend(&args[1], T::EXACT)?),
                );
            },
            | M::Zeros | M::Identity => {
                let rows = integer(&args[0])?;
                let cols = if matches!(self.method, M::Identity) {
                    rows
                } else {
                    integer(&args[1])?
                };
                if T::EXACT && (rows > 32 || cols > 32) {
                    return Err(error("exact matrices are limited to 32 rows and columns"));
                }
                if rows > SIDE || cols > SIDE || rows.checked_mul(cols).is_none_or(|n| n > LIMIT) {
                    return Err(error("matrix constructor exceeds allocation budget"));
                }
                if matches!(self.method, M::Zeros) && !T::decode(&args[2])?.is_one() {
                    return Err(error(
                        "zero matrix unit witness must equal one canonical unit",
                    ));
                }
                return encode_matrix(if matches!(self.method, M::Identity) {
                    Matrix::<T>::identity(rows)
                } else {
                    Matrix::<T>::zeros(rows, cols)
                });
            },
            | _ => {},
        }
        let mut a = decode_matrix::<T>(&args[0])?;
        if matches!(
            self.method,
            M::Frobenius | M::L1 | M::LInf | M::IsIdentity | M::Orthogonal
        ) && a
            .data()
            .iter()
            .any(|value| value.to_f64().is_none_or(|value| !value.is_finite()))
        {
            return Err(error(
                "approximate norm/tolerance query cannot represent an exact entry as Float",
            ));
        }
        match self.method {
            | M::WithBackend => encode_matrix(a.with_backend(backend(&args[1], T::EXACT)?)),
            | M::SetBackend => {
                a.set_backend(backend(&args[1], T::EXACT)?);
                encode_matrix(a)
            },
            | M::Get | M::Set => {
                let row = integer(&args[1])?;
                let col = integer(&args[2])?;
                if row >= a.rows() || col >= a.cols() {
                    return Err(error("matrix index outside shape"));
                }
                if matches!(self.method, M::Get) {
                    a.get(row, col).clone().encode()
                } else {
                    *a.get_mut(row, col) = T::decode(&args[3])?;
                    encode_matrix(a)
                }
            },
            | M::Rows => Ok(Value::int(a.rows())),
            | M::Cols => Ok(Value::int(a.cols())),
            | M::Data => encode_array(vec![a.data().len()], a.data().clone()),
            | M::IntoData => encode_array(vec![a.data().len()], a.into_data()),
            | M::GetCols => {
                let shape = vec![a.cols(), a.rows()];
                encode_array(shape, a.get_cols().into_iter().flatten().collect())
            },
            | M::Rref => {
                let rank = a.rref().map_err(error)?;
                Ok(record(vec![
                    ("matrix", encode_matrix(a)?),
                    ("rank", Value::int(rank)),
                ]))
            },
            | M::Transpose => encode_matrix(a.transpose()),
            | M::Strassen => {
                let b = decode_matrix::<T>(&args[1])?;
                let side = a
                    .rows()
                    .max(a.cols())
                    .max(b.rows())
                    .max(b.cols())
                    .checked_next_power_of_two()
                    .ok_or_else(|| error("Strassen padding overflows"))?;
                if side
                    .checked_mul(side)
                    .and_then(|n| n.checked_mul(12))
                    .is_none_or(|n| n > LIMIT * 12)
                {
                    return Err(error("Strassen padding exceeds work budget"));
                }
                encode_matrix(a.mul_strassen(&b).map_err(error)?)
            },
            | M::Det | M::DetBlock | M::DetLu => {
                let value = match self.method {
                    | M::Det => a.determinant(),
                    | M::DetBlock => a.determinant_block(),
                    | _ => a.determinant_lu(),
                }
                .map_err(error)?;
                Ok(record(vec![
                    ("normalized", value.encode()?),
                    ("entry_unit", T::one().encode()?),
                    ("order", Value::int(a.rows())),
                ]))
            },
            | M::Lu => {
                let (packed, swaps) = a.lu_decomposition().map_err(error)?;
                Ok(record(vec![
                    ("packed", encode_matrix(packed)?),
                    ("swaps", Value::int(swaps)),
                    ("entry_unit", T::one().encode()?),
                ]))
            },
            | M::Inverse => Ok(option(a.inverse().map(encode_matrix).transpose()?)),
            | M::NullSpace => encode_matrix(a.null_space().map_err(error)?),
            | M::Rank => Ok(Value::int(a.rank().map_err(error)?)),
            | M::Trace => a.trace().map_err(error)?.encode(),
            | M::Symmetric => Ok(Value::Bool(a.is_symmetric())),
            | M::Diagonal => Ok(Value::Bool(a.is_diagonal())),
            | M::Frobenius => super::output_scalar(a.frobenius_norm()),
            | M::L1 => super::output_scalar(a.l1_norm()),
            | M::LInf => super::output_scalar(a.linf_norm()),
            | M::IsIdentity | M::Orthogonal => {
                let tolerance = numeric(&args[1], 1)?;
                if tolerance < 0.0 {
                    return Err(error("matrix tolerance must be nonnegative"));
                }
                Ok(Value::Bool(if matches!(self.method, M::IsIdentity) {
                    a.is_identity(tolerance)
                } else {
                    a.is_orthogonal(tolerance)
                }))
            },
            | M::Jacobi
            | M::DecomposeLu
            | M::DecomposeQr
            | M::DecomposeCholesky
            | M::DecomposeSvd
            | M::DecomposeEigen => Err(error("decomposition specialization requires Float")),
            | M::New | M::Zeros | M::Identity => unreachable!(),
        }
    }
}
impl EvalDyn for MatrixAdapter {
    fn signature(&self) -> Signature {
        self.signature.clone()
    }

    #[allow(clippy::too_many_lines)] // Keep native method validation beside dispatch.
    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != self.signature.params.len() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.signature.params.len() as u8),
                got: args.len(),
            });
        }
        let dims = runtime_dimensions::prepare_contract(&self.contract, args)?;
        let result = if matches!(
            self.method,
            MatrixMethod::DecomposeLu
                | MatrixMethod::DecomposeQr
                | MatrixMethod::DecomposeCholesky
                | MatrixMethod::DecomposeSvd
                | MatrixMethod::DecomposeEigen
        ) {
            let a = decode_matrix::<f64>(&args[0])?;
            if matches!(
                self.method,
                MatrixMethod::DecomposeCholesky | MatrixMethod::DecomposeEigen
            ) && !a.is_symmetric()
            {
                return Err(error("symmetric decomposition requires symmetric input"));
            }
            let kind = match self.method {
                | MatrixMethod::DecomposeLu => FaerDecompositionType::Lu,
                | MatrixMethod::DecomposeQr => FaerDecompositionType::Qr,
                | MatrixMethod::DecomposeCholesky => FaerDecompositionType::Cholesky,
                | MatrixMethod::DecomposeSvd => FaerDecompositionType::Svd,
                | _ => FaerDecompositionType::EigenSymmetric,
            };
            match a
                .decompose(kind)
                .ok_or_else(|| error("native decomposition is unavailable or failed"))?
            {
                | FaerDecompositionResult::Lu { l, u, p } => {
                    record(vec![
                        ("l", encode_matrix(l)?),
                        ("u", encode_matrix(u)?),
                        (
                            "p",
                            Value::Tensor(TensorValue {
                                shape: vec![p.len()],
                                elem: TypeExpr::Int,
                                data: p.into_iter().map(Value::int).collect(),
                            }),
                        ),
                    ])
                },
                | FaerDecompositionResult::Qr { q, r } => {
                    record(vec![("q", encode_matrix(q)?), ("r", encode_matrix(r)?)])
                },
                | FaerDecompositionResult::Cholesky { l } => record(vec![("l", encode_matrix(l)?)]),
                | FaerDecompositionResult::Svd { u, s, v } => {
                    record(vec![
                        ("u", encode_matrix(u)?),
                        ("s", encode_array(vec![s.len()], s)?),
                        ("v", encode_matrix(v)?),
                    ])
                },
                | FaerDecompositionResult::EigenSymmetric { values, vectors } => {
                    record(vec![
                        ("values", encode_array(vec![values.len()], values)?),
                        ("vectors", encode_matrix(vectors)?),
                    ])
                },
            }
        } else if matches!(self.method, MatrixMethod::Jacobi) {
            let a = decode_matrix::<f64>(&args[0])?;
            if !a.is_symmetric() {
                return Err(error("Jacobi eigendecomposition requires symmetry"));
            }
            let sweeps = integer(&args[1])?;
            let tolerance = numeric(&args[2], 2)?;
            if sweeps == 0 || sweeps > 100000 || tolerance <= 0.0 {
                return Err(error(
                    "positive bounded sweep budget and tolerance required",
                ));
            }
            let (values, vectors) = a
                .jacobi_eigen_decomposition(sweeps, tolerance)
                .map_err(error)?;
            record(vec![
                ("values", encode_array(vec![values.len()], values)?),
                ("vectors", encode_matrix(vectors)?),
            ])
        } else if self.exact {
            self.run::<BigRational>(args)?
        } else {
            self.run::<f64>(args)?
        };
        runtime_dimensions::finish(&self.signature, &dims, result)
    }
}
pub(super) fn install(
    installer: &mut Installer<'_>,
    namespace: &str,
) -> Result<(), RuleError> {
    if namespace == "numerical_vector" {
        install_polynomials(installer)?;
        install_approximations(installer)?;
        multilinear::install(installer)?;
        return Ok(());
    }
    if namespace != "numerical_matrix" {
        return Ok(());
    }
    for exact in [false, true] {
        for (suffix, path, method) in MATRIX_METHODS {
            if exact
                && matches!(
                    method,
                    MatrixMethod::Jacobi
                        | MatrixMethod::DecomposeLu
                        | MatrixMethod::DecomposeQr
                        | MatrixMethod::DecomposeCholesky
                        | MatrixMethod::DecomposeSvd
                        | MatrixMethod::DecomposeEigen
                )
            {
                continue;
            }
            let name = format!(
                "num_matrix_{}_{}",
                if exact { "real" } else { "float" },
                suffix
            );
            let adapter = matrix_adapter(*method, exact);
            let signature = adapter.signature();
            let descriptor=OpDescriptor::new(&name,Arity::Fixed(signature.params.len()as u8)).rust_path(path).category("numerical_matrix")
                .doc("Validated portable Matrix specialization. Format/version, dense rank/shape, indices, finite values and allocation limits are checked. Real specializes exact rational arithmetic and permits only Native backend. Backend changes are immutable. Determinants return normalized magnitude, entry_unit and order, representing normalized*entry_unit^order; packed LU similarly operates on entries normalized by one canonical input unit. RREF/null-space vectors are dimensionless. No pointer or process-local handle is serialized.")
                .stability(Stability::Experimental).signature(signature).eval_dyn(Arc::new(adapter));
            let op = installer.op(descriptor)?;
            installer.dim(op, crate::graph::dim::DimRule::Container);
        }
    }
    factors::install(installer)?;
    Ok(())
}
const MATRIX_METHODS: &[(&str, &str, MatrixMethod)] = &[
    (
        "decompose_lu",
        "crate::kernels::matrix::Matrix<T>::decompose",
        MatrixMethod::DecomposeLu,
    ),
    (
        "decompose_qr",
        "crate::kernels::matrix::Matrix<T>::decompose",
        MatrixMethod::DecomposeQr,
    ),
    (
        "decompose_cholesky",
        "crate::kernels::matrix::Matrix<T>::decompose",
        MatrixMethod::DecomposeCholesky,
    ),
    (
        "decompose_svd",
        "crate::kernels::matrix::Matrix<T>::decompose",
        MatrixMethod::DecomposeSvd,
    ),
    (
        "decompose_eigen",
        "crate::kernels::matrix::Matrix<T>::decompose",
        MatrixMethod::DecomposeEigen,
    ),
    (
        "new",
        "crate::kernels::matrix::Matrix<T>::new",
        MatrixMethod::New,
    ),
    (
        "zeros",
        "crate::kernels::matrix::Matrix<T>::zeros",
        MatrixMethod::Zeros,
    ),
    (
        "identity",
        "crate::kernels::matrix::Matrix<T>::identity",
        MatrixMethod::Identity,
    ),
    (
        "with_backend",
        "crate::kernels::matrix::Matrix<T>::with_backend",
        MatrixMethod::WithBackend,
    ),
    (
        "set_backend",
        "crate::kernels::matrix::Matrix<T>::set_backend",
        MatrixMethod::SetBackend,
    ),
    (
        "get",
        "crate::kernels::matrix::Matrix<T>::get",
        MatrixMethod::Get,
    ),
    (
        "get_mut",
        "crate::kernels::matrix::Matrix<T>::get_mut",
        MatrixMethod::Set,
    ),
    (
        "rows",
        "crate::kernels::matrix::Matrix<T>::rows",
        MatrixMethod::Rows,
    ),
    (
        "cols",
        "crate::kernels::matrix::Matrix<T>::cols",
        MatrixMethod::Cols,
    ),
    (
        "data",
        "crate::kernels::matrix::Matrix<T>::data",
        MatrixMethod::Data,
    ),
    (
        "into_data",
        "crate::kernels::matrix::Matrix<T>::into_data",
        MatrixMethod::IntoData,
    ),
    (
        "get_cols",
        "crate::kernels::matrix::Matrix<T>::get_cols",
        MatrixMethod::GetCols,
    ),
    (
        "rref",
        "crate::kernels::matrix::Matrix<T>::rref",
        MatrixMethod::Rref,
    ),
    (
        "transpose",
        "crate::kernels::matrix::Matrix<T>::transpose",
        MatrixMethod::Transpose,
    ),
    (
        "mul_strassen",
        "crate::kernels::matrix::Matrix<T>::mul_strassen",
        MatrixMethod::Strassen,
    ),
    (
        "determinant",
        "crate::kernels::matrix::Matrix<T>::determinant",
        MatrixMethod::Det,
    ),
    (
        "determinant_block",
        "crate::kernels::matrix::Matrix<T>::determinant_block",
        MatrixMethod::DetBlock,
    ),
    (
        "determinant_lu",
        "crate::kernels::matrix::Matrix<T>::determinant_lu",
        MatrixMethod::DetLu,
    ),
    (
        "lu_decomposition",
        "crate::kernels::matrix::Matrix<T>::lu_decomposition",
        MatrixMethod::Lu,
    ),
    (
        "inverse",
        "crate::kernels::matrix::Matrix<T>::inverse",
        MatrixMethod::Inverse,
    ),
    (
        "null_space",
        "crate::kernels::matrix::Matrix<T>::null_space",
        MatrixMethod::NullSpace,
    ),
    (
        "rank",
        "crate::kernels::matrix::Matrix<T>::rank",
        MatrixMethod::Rank,
    ),
    (
        "trace",
        "crate::kernels::matrix::Matrix<T>::trace",
        MatrixMethod::Trace,
    ),
    (
        "is_symmetric",
        "crate::kernels::matrix::Matrix<T>::is_symmetric",
        MatrixMethod::Symmetric,
    ),
    (
        "is_diagonal",
        "crate::kernels::matrix::Matrix<T>::is_diagonal",
        MatrixMethod::Diagonal,
    ),
    (
        "frobenius_norm",
        "crate::kernels::matrix::Matrix<T>::frobenius_norm",
        MatrixMethod::Frobenius,
    ),
    (
        "l1_norm",
        "crate::kernels::matrix::Matrix<T>::l1_norm",
        MatrixMethod::L1,
    ),
    (
        "linf_norm",
        "crate::kernels::matrix::Matrix<T>::linf_norm",
        MatrixMethod::LInf,
    ),
    (
        "is_identity",
        "crate::kernels::matrix::Matrix<T>::is_identity",
        MatrixMethod::IsIdentity,
    ),
    (
        "is_orthogonal",
        "crate::kernels::matrix::Matrix<T>::is_orthogonal",
        MatrixMethod::Orthogonal,
    ),
    (
        "jacobi_eigen_decomposition",
        "crate::kernels::matrix::Matrix<f64>::jacobi_eigen_decomposition",
        MatrixMethod::Jacobi,
    ),
];

// A normalized polynomial carries ordinate coefficients and one physical axis
// scale, so c[k]*(x/scale)^k never mixes different physical powers in a tensor.
const POLY_FIELDS: &[&str] = &["format", "coefficients", "coordinate_scale"];
const QUOTIENT: DimFormula = DimFormula::Mul(&[D, DimFormula::Pow(&E, -1, 1)]);
const X: DimFormula = DimFormula::Var("X");
const DERIVATIVE: DimFormula = DimFormula::Mul(&[D, DimFormula::Pow(&X, -1, 1)]);
const INTEGRAL: DimFormula = DimFormula::Mul(&[D, X]);
fn polynomial_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        ("coefficients".into(), tensor_type(1, TypeExpr::Float)),
        ("coordinate_scale".into(), TypeExpr::Float),
    ])
}
fn decode_polynomial(value: &Value) -> Result<(Polynomial, f64), OpError> {
    kind(
        field(value, "format", POLY_FIELDS)?,
        "rssn/polynomial-normalized/v1",
    )?;
    let (_, coefficients) = array::<f64>(field(value, "coefficients", POLY_FIELDS)?, 1)?;
    if coefficients.is_empty() || coefficients.len() > 129 {
        return Err(error(
            "polynomial must contain 1..=129 descending coefficients",
        ));
    }
    let scale = numeric(field(value, "coordinate_scale", POLY_FIELDS)?, 0)?;
    if scale <= 0.0 {
        return Err(error("coordinate scale must be positive"));
    }
    Ok((Polynomial::new(coefficients), scale))
}
fn encode_polynomial(
    mut polynomial: Polynomial,
    scale: f64,
) -> Result<Value, OpError> {
    if polynomial.coeffs.is_empty() {
        polynomial.coeffs.push(0.0);
    }
    if polynomial.coeffs.len() > 129 {
        return Err(error("polynomial result exceeds degree budget"));
    }
    Ok(record(vec![
        (
            "format",
            Value::String("rssn/polynomial-normalized/v1".into()),
        ),
        (
            "coefficients",
            encode_array(vec![polynomial.coeffs.len()], polynomial.coeffs)?,
        ),
        ("coordinate_scale", super::output_scalar(scale)?),
    ]))
}
#[derive(Clone, Copy)]
enum PolynomialMethod {
    New,
    Eval,
    Roots,
    Derivative,
    Integral,
    Degree,
    IsZero,
    LongDivision,
    TryLongDivision,
    Add,
    Sub,
    Mul,
    Div,
    Scale,
    DivScalar,
}
struct PolynomialAdapter {
    method: PolynomialMethod,
    signature: Signature,
    contract: DimContract,
}
#[allow(clippy::too_many_lines)] // Reviewed signature and dimension cases remain adjacent.
fn polynomial_adapter(method: PolynomialMethod) -> PolynomialAdapter {
    use PolynomialMethod as P;
    let model = polynomial_type();
    let mut args = vec![model.clone()];
    let mut ret = model.clone();
    let mut inputs = vec![
        input("D", 0, &["coefficients"]),
        input("X", 0, &["coordinate_scale"]),
    ];
    let mut checks = vec![];
    let mut outputs = vec![
        output(&["coefficients"], &D),
        output(&["coordinate_scale"], &X),
    ];
    let mut result = None;
    match method {
        | P::New => {
            args = vec![tensor_type(1, TypeExpr::Float), TypeExpr::Float];
            inputs = vec![input("D", 0, &[]), input("X", 1, &[])];
        },
        | P::Eval => {
            args.push(TypeExpr::Float);
            checks.push(constraint(1, &[], &X));
            ret = TypeExpr::Float;
            outputs.clear();
            result = Some(&D);
        },
        | P::Roots => {
            ret = tensor_type(1, TypeExpr::Float);
            outputs.clear();
            result = Some(&X);
        },
        | P::Derivative => outputs[0] = output(&["coefficients"], &DERIVATIVE),
        | P::Integral => outputs[0] = output(&["coefficients"], &INTEGRAL),
        | P::Degree => {
            ret = TypeExpr::Int;
            outputs.clear();
        },
        | P::IsZero => {
            args.push(TypeExpr::Float);
            checks.push(constraint(1, &[], &D));
            ret = TypeExpr::Bool;
            outputs.clear();
        },
        | P::Add | P::Sub | P::Mul | P::Div | P::LongDivision | P::TryLongDivision => {
            args.push(model.clone());
            checks.push(constraint(1, &["coordinate_scale"], &X));
            if matches!(method, P::Add | P::Sub) {
                checks.push(constraint(1, &["coefficients"], &D));
            } else {
                inputs.push(input("E", 1, &["coefficients"]));
            }
            if matches!(method, P::Mul) {
                outputs[0] = output(&["coefficients"], &PRODUCT);
            }
            if matches!(method, P::Div) {
                outputs[0] = output(&["coefficients"], &QUOTIENT);
            }
            if matches!(method, P::LongDivision | P::TryLongDivision) {
                ret = TypeExpr::Record(vec![
                    ("quotient".into(), model.clone()),
                    ("remainder".into(), model),
                ]);
                outputs = vec![
                    output(&["quotient", "coefficients"], &QUOTIENT),
                    output(&["quotient", "coordinate_scale"], &X),
                    output(&["remainder", "coefficients"], &D),
                    output(&["remainder", "coordinate_scale"], &X),
                ];
            }
        },
        | P::Scale | P::DivScalar => {
            args.push(TypeExpr::Float);
            inputs.push(input("E", 1, &[]));
            outputs[0] = output(
                &["coefficients"],
                if matches!(method, P::Scale) {
                    &PRODUCT
                } else {
                    &QUOTIENT
                },
            );
        },
    }
    let contract = contract(inputs, checks, outputs, result);
    let signature = dimensions::apply_contract(Signature::new(args, ret), &contract);
    PolynomialAdapter {
        method,
        signature,
        contract,
    }
}
impl EvalDyn for PolynomialAdapter {
    fn signature(&self) -> Signature {
        self.signature.clone()
    }

    #[allow(clippy::too_many_lines)] // Keep native method validation beside dispatch.
    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use PolynomialMethod as P;
        if args.len() != self.signature.params.len() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.signature.params.len() as u8),
                got: args.len(),
            });
        }
        let dims = runtime_dimensions::prepare_contract(&self.contract, args)?;
        let result = if matches!(self.method, P::New) {
            let (_, coefficients) = array::<f64>(&args[0], 1)?;
            let scale = numeric(&args[1], 1)?;
            if coefficients.is_empty() || coefficients.len() > 65 || scale <= 0.0 {
                return Err(error(
                    "polynomial requires 1..=65 coefficients and positive coordinate scale",
                ));
            }
            encode_polynomial(Polynomial::new(coefficients), scale)?
        } else {
            let (a, scale) = decode_polynomial(&args[0])?;
            match self.method {
                | P::Eval => super::output_scalar(a.eval(numeric(&args[1], 1)? / scale))?,
                | P::Roots => {
                    if a.coeffs.iter().all(|v| crate::kernels::exact_eq(*v, 0.0)) {
                        return Err(error("zero polynomial has no finite isolated root list"));
                    }
                    let roots = a.find_roots().map_err(error)?;
                    encode_array(
                        vec![roots.len()],
                        roots.into_iter().map(|v| v * scale).collect(),
                    )?
                },
                | P::Derivative => {
                    encode_polynomial(a.derivative().div_scalar(scale).map_err(error)?, scale)?
                },
                | P::Integral => encode_polynomial(a.integral() * scale, scale)?,
                | P::Degree => Value::int(a.degree()),
                | P::IsZero => {
                    let epsilon = numeric(&args[1], 1)?;
                    if epsilon < 0.0 {
                        return Err(error("polynomial zero tolerance must be nonnegative"));
                    }
                    Value::Bool(a.is_zero(epsilon))
                },
                | P::Scale => encode_polynomial(a * numeric(&args[1], 1)?, scale)?,
                | P::DivScalar => {
                    encode_polynomial(a.div_scalar(numeric(&args[1], 1)?).map_err(error)?, scale)?
                },
                | P::Add | P::Sub | P::Mul | P::Div | P::LongDivision | P::TryLongDivision => {
                    let (b, other_scale) = decode_polynomial(&args[1])?;
                    if !crate::kernels::exact_eq(scale, other_scale) {
                        return Err(error(
                            "polynomial arithmetic requires identical physical coordinate scales",
                        ));
                    }
                    if matches!(self.method, P::Div | P::LongDivision)
                        && b.coeffs
                            .first()
                            .is_none_or(|v| crate::kernels::exact_eq(*v, 0.0))
                    {
                        return Err(error("divisor must have a nonzero leading coefficient"));
                    }
                    match self.method {
                        | P::Add => encode_polynomial(a + b, scale)?,
                        | P::Sub => encode_polynomial(a - b, scale)?,
                        | P::Mul => encode_polynomial(a * b, scale)?,
                        | P::Div => encode_polynomial(a / b, scale)?,
                        | _ => {
                            let (q, r) = if matches!(self.method, P::TryLongDivision) {
                                a.try_long_division(&b).map_err(error)?
                            } else {
                                a.long_division(&b)
                            };
                            record(vec![
                                ("quotient", encode_polynomial(q, scale)?),
                                ("remainder", encode_polynomial(r, scale)?),
                            ])
                        },
                    }
                },
                | P::New => unreachable!(),
            }
        };
        runtime_dimensions::finish(&self.signature, &dims, result)
    }
}
fn install_polynomials(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (suffix, path, method) in POLYNOMIAL_METHODS {
        let name = format!("num_polynomial_model_{suffix}");
        let adapter = polynomial_adapter(*method);
        let signature = adapter.signature();
        let doc = if matches!(method, PolynomialMethod::TryLongDivision) {
            "Divide measured descending-coefficient polynomials with the same coordinate reference. Leading zeros in the divisor are skipped. Returns quotient coefficients in dividend units divided by divisor units and remainder coefficients in dividend units. Zero divisors and nonfinite results are errors."
        } else {
            "Portable descending-coefficient Polynomial of normalized x/coordinate_scale. Coefficients retain ordinate dimensions; derivative and integral restore physical scale powers. Arithmetic requires identical coordinate scales. Kind/version, positive scale, finite values and degree/allocation budgets are validated. Native polynomial algorithms are called directly; no pointer or process-local handle is exposed."
        };
        let mut descriptor = OpDescriptor::new(&name, Arity::Fixed(signature.params.len() as u8))
            .category("numerical_vector")
            .doc(doc)
            .stability(Stability::Experimental)
            .signature(signature)
            .eval_dyn(Arc::new(adapter));
        if !path.is_empty() {
            descriptor = descriptor.rust_path(path);
        }
        let op = installer.op(descriptor)?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
const POLYNOMIAL_METHODS: &[(&str, &str, PolynomialMethod)] = &[
    (
        "new",
        "crate::kernels::polynomial::Polynomial::new",
        PolynomialMethod::New,
    ),
    (
        "eval",
        "crate::kernels::polynomial::Polynomial::eval",
        PolynomialMethod::Eval,
    ),
    (
        "find_roots",
        "crate::kernels::polynomial::Polynomial::find_roots",
        PolynomialMethod::Roots,
    ),
    (
        "derivative",
        "crate::kernels::polynomial::Polynomial::derivative",
        PolynomialMethod::Derivative,
    ),
    (
        "integral",
        "crate::kernels::polynomial::Polynomial::integral",
        PolynomialMethod::Integral,
    ),
    (
        "degree",
        "crate::kernels::polynomial::Polynomial::degree",
        PolynomialMethod::Degree,
    ),
    (
        "is_zero",
        "crate::kernels::polynomial::Polynomial::is_zero",
        PolynomialMethod::IsZero,
    ),
    (
        "long_division",
        "crate::kernels::polynomial::Polynomial::long_division",
        PolynomialMethod::LongDivision,
    ),
    (
        "try_long_division",
        "crate::kernels::polynomial::Polynomial::try_long_division",
        PolynomialMethod::TryLongDivision,
    ),
    ("add", "", PolynomialMethod::Add),
    ("sub", "", PolynomialMethod::Sub),
    ("mul", "", PolynomialMethod::Mul),
    ("div", "", PolynomialMethod::Div),
    ("scale", "", PolynomialMethod::Scale),
    (
        "div_scalar",
        "crate::kernels::polynomial::Polynomial::div_scalar",
        PolynomialMethod::DivScalar,
    ),
];

use crate::kernels::approx::Aaa;
use crate::kernels::approx::ChebSeries;
use crate::kernels::approx::CubicSpline;
use crate::kernels::approx::Pchip;
use crate::kernels::approx::SplineBc;
use crate::kernels::approx::{
    self,
};
const CURVE_FIELDS: &[&str] = &["format", "x", "y", "boundary", "left", "right"];
const CHEB_FIELDS: &[&str] = &["format", "coefficients", "a", "b"];
const AAA_FIELDS: &[&str] = &["format", "z", "f", "weights", "error"];
fn curve_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        ("x".into(), tensor_type(1, TypeExpr::Float)),
        ("y".into(), tensor_type(1, TypeExpr::Float)),
        ("boundary".into(), TypeExpr::Int),
        ("left".into(), TypeExpr::Float),
        ("right".into(), TypeExpr::Float),
    ])
}
pub(in crate::rules) fn cheb_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        ("coefficients".into(), tensor_type(1, TypeExpr::Float)),
        ("a".into(), TypeExpr::Float),
        ("b".into(), TypeExpr::Float),
    ])
}
fn aaa_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        ("z".into(), tensor_type(1, TypeExpr::Float)),
        ("f".into(), tensor_type(1, TypeExpr::Float)),
        ("weights".into(), tensor_type(1, TypeExpr::Float)),
        ("error".into(), TypeExpr::Float),
    ])
}
struct CurveData {
    x: Vec<f64>,
    y: Vec<f64>,
    bc: SplineBc,
    boundary: usize,
    left: f64,
    right: f64,
}
fn curve_data(
    value: &Value,
    pchip: bool,
) -> Result<CurveData, OpError> {
    kind(
        field(value, "format", CURVE_FIELDS)?,
        if pchip {
            "rssn/pchip/v1"
        } else {
            "rssn/cubic-spline/v1"
        },
    )?;
    let (_, x) = array::<f64>(field(value, "x", CURVE_FIELDS)?, 1)?;
    let (_, y) = array::<f64>(field(value, "y", CURVE_FIELDS)?, 1)?;
    let boundary = integer(field(value, "boundary", CURVE_FIELDS)?)?;
    let left = numeric(field(value, "left", CURVE_FIELDS)?, 0)?;
    let right = numeric(field(value, "right", CURVE_FIELDS)?, 0)?;
    let bc = match boundary {
        | 0 => SplineBc::Natural,
        | 1 => SplineBc::Clamped(left, right),
        | 2 => SplineBc::NotAKnot,
        | _ => {
            return Err(error(
                "spline boundary must be 0=natural, 1=clamped, 2=not-a-knot",
            ));
        },
    };
    if x.len() != y.len()
        || x.len() < 2
        || x.windows(2).any(|w| w[0] >= w[1])
        || pchip && boundary != 0
    {
        return Err(error("invalid interpolation grid or PCHIP boundary"));
    }
    Ok(CurveData {
        x,
        y,
        bc,
        boundary,
        left,
        right,
    })
}
fn encode_curve(
    data: CurveData,
    pchip: bool,
) -> Result<Value, OpError> {
    Ok(record(vec![
        (
            "format",
            Value::String(
                if pchip {
                    "rssn/pchip/v1"
                } else {
                    "rssn/cubic-spline/v1"
                }
                .into(),
            ),
        ),
        ("x", encode_array(vec![data.x.len()], data.x)?),
        ("y", encode_array(vec![data.y.len()], data.y)?),
        ("boundary", Value::int(data.boundary)),
        ("left", super::output_scalar(data.left)?),
        ("right", super::output_scalar(data.right)?),
    ]))
}
fn decode_cheb(value: &Value) -> Result<ChebSeries, OpError> {
    kind(field(value, "format", CHEB_FIELDS)?, "rssn/chebyshev/v1")?;
    let (_, coeffs) = array::<f64>(field(value, "coefficients", CHEB_FIELDS)?, 1)?;
    let a = numeric(field(value, "a", CHEB_FIELDS)?, 0)?;
    let b = numeric(field(value, "b", CHEB_FIELDS)?, 0)?;
    if coeffs.is_empty() || coeffs.len() > 65 || a >= b {
        return Err(error(
            "Chebyshev needs 1..=65 coefficients on an increasing finite interval",
        ));
    }
    Ok(ChebSeries { coeffs, a, b })
}
pub(in crate::rules) fn encode_cheb(series: ChebSeries) -> Result<Value, OpError> {
    Ok(record(vec![
        ("format", Value::String("rssn/chebyshev/v1".into())),
        (
            "coefficients",
            encode_array(vec![series.coeffs.len()], series.coeffs)?,
        ),
        ("a", super::output_scalar(series.a)?),
        ("b", super::output_scalar(series.b)?),
    ]))
}
fn decode_aaa(value: &Value) -> Result<Aaa, OpError> {
    kind(field(value, "format", AAA_FIELDS)?, "rssn/aaa/v1")?;
    let (_, z) = array::<f64>(field(value, "z", AAA_FIELDS)?, 1)?;
    let (_, f) = array::<f64>(field(value, "f", AAA_FIELDS)?, 1)?;
    let (_, w) = array::<f64>(field(value, "weights", AAA_FIELDS)?, 1)?;
    let error_value = numeric(field(value, "error", AAA_FIELDS)?, 0)?;
    if z.is_empty()
        || z.len() != f.len()
        || z.len() != w.len()
        || error_value < 0.0
        || w.iter().all(|v| crate::kernels::exact_eq(*v, 0.0))
        || z.iter()
            .enumerate()
            .any(|(i, v)| z[..i].iter().any(|u| crate::kernels::exact_eq(*u, *v)))
    {
        return Err(error("invalid AAA support, weights or error"));
    }
    Ok(Aaa {
        z,
        f,
        w,
        error: error_value,
    })
}
fn encode_aaa(model: Aaa) -> Result<Value, OpError> {
    Ok(record(vec![
        ("format", Value::String("rssn/aaa/v1".into())),
        ("z", encode_array(vec![model.z.len()], model.z)?),
        ("f", encode_array(vec![model.f.len()], model.f)?),
        ("weights", encode_array(vec![model.w.len()], model.w)?),
        ("error", super::output_scalar(model.error)?),
    ]))
}
#[derive(Clone, Copy)]
enum ApproxMethod {
    CubicNew,
    CubicEval,
    CubicD1,
    CubicD2,
    CubicD3,
    CubicIntegral,
    PchipNew,
    PchipEval,
    PchipDerivative,
    ChebNew,
    ChebFromValues,
    ChebEval,
    ChebDegree,
    ChebDerivative,
    ChebIntegral,
    ChebDefinite,
    ChebRoots,
    AaaNew,
    AaaEval,
    AaaDegree,
}
struct ApproxAdapter {
    method: ApproxMethod,
    signature: Signature,
    contract: DimContract,
}
const D2: DimFormula = DimFormula::Mul(&[D, DimFormula::Pow(&X, -2, 1)]);
const D3: DimFormula = DimFormula::Mul(&[D, DimFormula::Pow(&X, -3, 1)]);
#[allow(clippy::too_many_lines)] // Reviewed signature and dimension cases remain adjacent.
fn approximation_adapter(method: ApproxMethod) -> ApproxAdapter {
    use ApproxMethod as A;
    let curve = matches!(
        method,
        A::CubicNew
            | A::CubicEval
            | A::CubicD1
            | A::CubicD2
            | A::CubicD3
            | A::CubicIntegral
            | A::PchipNew
            | A::PchipEval
            | A::PchipDerivative
    );
    let aaa = matches!(method, A::AaaNew | A::AaaEval | A::AaaDegree);
    let model = if curve {
        curve_type()
    } else if aaa {
        aaa_type()
    } else {
        cheb_type()
    };
    let mut args = vec![model.clone()];
    let mut ret = model;
    let mut inputs = if curve {
        vec![input("D", 0, &["y"]), input("X", 0, &["x"])]
    } else if aaa {
        vec![input("D", 0, &["f"]), input("X", 0, &["z"])]
    } else {
        vec![input("D", 0, &["coefficients"]), input("X", 0, &["a"])]
    };
    let mut checks = if curve {
        vec![
            constraint(0, &["left"], &DERIVATIVE),
            constraint(0, &["right"], &DERIVATIVE),
        ]
    } else if aaa {
        vec![constraint(0, &["error"], &D)]
    } else {
        vec![constraint(0, &["b"], &X)]
    };
    let mut outputs = if curve {
        vec![
            output(&["x"], &X),
            output(&["y"], &D),
            output(&["left"], &DERIVATIVE),
            output(&["right"], &DERIVATIVE),
        ]
    } else if aaa {
        vec![
            output(&["z"], &X),
            output(&["f"], &D),
            output(&["error"], &D),
        ]
    } else {
        vec![
            output(&["coefficients"], &D),
            output(&["a"], &X),
            output(&["b"], &X),
        ]
    };
    let mut result = None;
    match method {
        | A::CubicNew => {
            args = vec![
                tensor_type(1, TypeExpr::Float),
                tensor_type(1, TypeExpr::Float),
                TypeExpr::Int,
                TypeExpr::Float,
                TypeExpr::Float,
            ];
            inputs = vec![input("X", 0, &[]), input("D", 1, &[])];
            checks = vec![
                constraint(3, &[], &DERIVATIVE),
                constraint(4, &[], &DERIVATIVE),
            ];
        },
        | A::PchipNew => {
            args = vec![
                tensor_type(1, TypeExpr::Float),
                tensor_type(1, TypeExpr::Float),
            ];
            inputs = vec![input("X", 0, &[]), input("D", 1, &[])];
            checks.clear();
        },
        | A::ChebNew | A::ChebFromValues => {
            args = vec![
                tensor_type(1, TypeExpr::Float),
                TypeExpr::Float,
                TypeExpr::Float,
            ];
            inputs = vec![input("D", 0, &[]), input("X", 1, &[])];
            checks = vec![constraint(2, &[], &X)];
        },
        | A::AaaNew => {
            args = vec![
                tensor_type(1, TypeExpr::Float),
                tensor_type(1, TypeExpr::Float),
                TypeExpr::Float,
                TypeExpr::Int,
            ];
            inputs = vec![input("X", 0, &[]), input("D", 1, &[])];
            checks.clear();
        },
        | A::CubicEval
        | A::PchipEval
        | A::ChebEval
        | A::AaaEval
        | A::CubicD1
        | A::CubicD2
        | A::CubicD3
        | A::PchipDerivative => {
            args.push(TypeExpr::Float);
            checks.push(constraint(1, &[], &X));
            outputs.clear();
            ret = TypeExpr::Float;
            result = Some(match method {
                | A::CubicD1 | A::PchipDerivative => &DERIVATIVE,
                | A::CubicD2 => &D2,
                | A::CubicD3 => &D3,
                | _ => &D,
            });
        },
        | A::CubicIntegral => {
            args.extend([TypeExpr::Float, TypeExpr::Float]);
            checks.extend([constraint(1, &[], &X), constraint(2, &[], &X)]);
            outputs.clear();
            ret = TypeExpr::Float;
            result = Some(&INTEGRAL);
        },
        | A::ChebDegree | A::AaaDegree => {
            ret = TypeExpr::Int;
            outputs.clear();
        },
        | A::ChebDerivative => outputs[0] = output(&["coefficients"], &DERIVATIVE),
        | A::ChebIntegral => outputs[0] = output(&["coefficients"], &INTEGRAL),
        | A::ChebDefinite => {
            ret = TypeExpr::Float;
            outputs.clear();
            result = Some(&INTEGRAL);
        },
        | A::ChebRoots => {
            ret = tensor_type(1, TypeExpr::Float);
            outputs.clear();
            result = Some(&X);
        },
    }
    let contract = contract(inputs, checks, outputs, result);
    let signature = dimensions::apply_contract(Signature::new(args, ret), &contract);
    ApproxAdapter {
        method,
        signature,
        contract,
    }
}
impl EvalDyn for ApproxAdapter {
    fn signature(&self) -> Signature {
        self.signature.clone()
    }

    #[allow(clippy::too_many_lines)] // Keep native method validation beside dispatch.
    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use ApproxMethod as A;
        if args.len() != self.signature.params.len() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.signature.params.len() as u8),
                got: args.len(),
            });
        }
        let dims = runtime_dimensions::prepare_contract(&self.contract, args)?;
        let result = match self.method {
            | A::CubicNew | A::PchipNew => {
                let (_, x) = array::<f64>(&args[0], 1)?;
                let (_, y) = array::<f64>(&args[1], 1)?;
                let pchip = matches!(self.method, A::PchipNew);
                let (boundary, left, right) = if pchip {
                    (0, 0.0, 0.0)
                } else {
                    (
                        integer(&args[2])?,
                        numeric(&args[3], 3)?,
                        numeric(&args[4], 4)?,
                    )
                };
                let bc = match boundary {
                    | 0 => SplineBc::Natural,
                    | 1 => SplineBc::Clamped(left, right),
                    | 2 => SplineBc::NotAKnot,
                    | _ => return Err(error("invalid spline boundary selector")),
                };
                if pchip {
                    Pchip::validate_samples(&x, &y).map_err(error)?;
                } else {
                    CubicSpline::new(&x, &y, bc).map_err(error)?;
                }
                encode_curve(
                    CurveData {
                        x,
                        y,
                        bc,
                        boundary,
                        left,
                        right,
                    },
                    pchip,
                )?
            },
            | A::CubicEval | A::CubicD1 | A::CubicD2 | A::CubicD3 | A::CubicIntegral => {
                let data = curve_data(&args[0], false)?;
                let spline = CubicSpline::new(&data.x, &data.y, data.bc).map_err(error)?;
                let x = numeric(&args[1], 1)?;
                let value = match self.method {
                    | A::CubicEval => spline.eval(x),
                    | A::CubicIntegral => spline.integrate(x, numeric(&args[2], 2)?),
                    | _ => {
                        spline.derivative(
                            x,
                            match self.method {
                                | A::CubicD1 => 1,
                                | A::CubicD2 => 2,
                                | _ => 3,
                            },
                        )
                    },
                };
                super::output_scalar(value)?
            },
            | A::PchipEval | A::PchipDerivative => {
                let data = curve_data(&args[0], true)?;
                let x = numeric(&args[1], 1)?;
                super::output_scalar(Pchip::query_checked_samples(
                    &data.x,
                    &data.y,
                    x,
                    matches!(self.method, A::PchipDerivative),
                ))?
            },
            | A::ChebNew | A::ChebFromValues => {
                let (_, samples) = array::<f64>(&args[0], 1)?;
                let a = numeric(&args[1], 1)?;
                let b = numeric(&args[2], 2)?;
                if samples.is_empty() || samples.len() > 65 || a >= b {
                    return Err(error(
                        "Chebyshev model requires 1 to 65 values and increasing endpoints",
                    ));
                }
                let series = if matches!(self.method, A::ChebFromValues) {
                    ChebSeries::from_values(&samples, a, b)
                } else {
                    ChebSeries {
                        coeffs: samples,
                        a,
                        b,
                    }
                };
                encode_cheb(series)?
            },
            | A::ChebEval
            | A::ChebDegree
            | A::ChebDerivative
            | A::ChebIntegral
            | A::ChebDefinite
            | A::ChebRoots => {
                let series = decode_cheb(&args[0])?;
                match self.method {
                    | A::ChebEval => super::output_scalar(series.eval(numeric(&args[1], 1)?))?,
                    | A::ChebDegree => Value::int(series.degree()),
                    | A::ChebDerivative => encode_cheb(series.derivative())?,
                    | A::ChebIntegral => encode_cheb(series.integral())?,
                    | A::ChebDefinite => super::output_scalar(series.definite_integral())?,
                    | _ => {
                        let roots = series.roots().map_err(error)?;
                        encode_array(vec![roots.len()], roots)?
                    },
                }
            },
            | A::AaaNew => {
                let (_, z) = array::<f64>(&args[0], 1)?;
                let (_, f) = array::<f64>(&args[1], 1)?;
                let tolerance = numeric(&args[2], 2)?;
                let terms = integer(&args[3])?;
                if tolerance <= 0.0
                    || terms == 0
                    || terms > 64
                    || z.iter()
                        .enumerate()
                        .any(|(i, v)| z[..i].iter().any(|u| crate::kernels::exact_eq(*u, *v)))
                {
                    return Err(error(
                        "AAA needs distinct sites, positive relative tolerance, and 1..=64 terms",
                    ));
                }
                encode_aaa(approx::aaa(&z, &f, tolerance, terms).map_err(error)?)?
            },
            | A::AaaEval | A::AaaDegree => {
                let model = decode_aaa(&args[0])?;
                if matches!(self.method, A::AaaDegree) {
                    Value::int(model.degree())
                } else {
                    super::output_scalar(model.eval(numeric(&args[1], 1)?))?
                }
            },
        };
        runtime_dimensions::finish(&self.signature, &dims, result)
    }
}
fn install_approximations(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in APPROX_METHODS {
        let adapter = approximation_adapter(*method);
        let signature = adapter.signature();
        let mut descriptor=OpDescriptor::new(name,Arity::Fixed(signature.params.len()as u8)).category("numerical_vector")
            .doc("Portable native approximation model with independent coordinate and ordinate dimensions. Every query checks kind/version, defining data lengths, finite data, monotone grids where required, interval and allocation budgets. Derivatives restore ordinate/coordinate powers; integrals restore their product. Model values contain only defining data and no process-local references.")
            .stability(Stability::Experimental).signature(signature).eval_dyn(Arc::new(adapter));
        if !path.is_empty() {
            descriptor = descriptor.rust_path(path);
        }
        let op = installer.op(descriptor)?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
const APPROX_METHODS: &[(&str, &str, ApproxMethod)] = &[
    (
        "num_native_cubic_spline_new",
        "crate::kernels::approx::CubicSpline::new",
        ApproxMethod::CubicNew,
    ),
    (
        "num_native_cubic_spline_eval",
        "crate::kernels::approx::CubicSpline::eval",
        ApproxMethod::CubicEval,
    ),
    (
        "num_native_cubic_spline_d1",
        "crate::kernels::approx::CubicSpline::derivative",
        ApproxMethod::CubicD1,
    ),
    (
        "num_native_cubic_spline_d2",
        "crate::kernels::approx::CubicSpline::derivative",
        ApproxMethod::CubicD2,
    ),
    (
        "num_native_cubic_spline_d3",
        "crate::kernels::approx::CubicSpline::derivative",
        ApproxMethod::CubicD3,
    ),
    (
        "num_native_cubic_spline_integrate",
        "crate::kernels::approx::CubicSpline::integrate",
        ApproxMethod::CubicIntegral,
    ),
    (
        "num_native_pchip_new",
        "crate::kernels::approx::Pchip::new",
        ApproxMethod::PchipNew,
    ),
    (
        "num_native_pchip_eval",
        "crate::kernels::approx::Pchip::eval",
        ApproxMethod::PchipEval,
    ),
    (
        "num_native_pchip_derivative",
        "crate::kernels::approx::Pchip::derivative",
        ApproxMethod::PchipDerivative,
    ),
    ("num_native_cheb_coefficients", "", ApproxMethod::ChebNew),
    (
        "num_native_cheb_from_values",
        "crate::kernels::approx::ChebSeries::from_values",
        ApproxMethod::ChebFromValues,
    ),
    (
        "num_native_cheb_eval",
        "crate::kernels::approx::ChebSeries::eval",
        ApproxMethod::ChebEval,
    ),
    (
        "num_native_cheb_degree",
        "crate::kernels::approx::ChebSeries::degree",
        ApproxMethod::ChebDegree,
    ),
    (
        "num_native_cheb_derivative",
        "crate::kernels::approx::ChebSeries::derivative",
        ApproxMethod::ChebDerivative,
    ),
    (
        "num_native_cheb_integral",
        "crate::kernels::approx::ChebSeries::integral",
        ApproxMethod::ChebIntegral,
    ),
    (
        "num_native_cheb_definite_integral",
        "crate::kernels::approx::ChebSeries::definite_integral",
        ApproxMethod::ChebDefinite,
    ),
    (
        "num_native_cheb_roots",
        "crate::kernels::approx::ChebSeries::roots",
        ApproxMethod::ChebRoots,
    ),
    (
        "num_native_aaa_build",
        "crate::kernels::approx::aaa",
        ApproxMethod::AaaNew,
    ),
    (
        "num_native_aaa_eval",
        "crate::kernels::approx::Aaa::eval",
        ApproxMethod::AaaEval,
    ),
    (
        "num_native_aaa_degree",
        "crate::kernels::approx::Aaa::degree",
        ApproxMethod::AaaDegree,
    ),
];
