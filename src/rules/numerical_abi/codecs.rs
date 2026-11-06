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
include!("structured_generated.rs");

/// Build a record's physical result schema from generated field formulas.
pub(in crate::rules) fn structured_type<T: Output>(fields: &[(&str, super::dimensions::FieldDimension)], d: &crate::kernel::signature::DimExpr, e: &crate::kernel::signature::DimExpr) -> TypeExpr {
    let mut ty = T::ty();
    super::dimensions::apply_record_fields(&mut ty, fields, d, e);
    ty
}

/// Encode a validated record and apply the same field formulas as its schema.
pub(in crate::rules) fn structured_value<T: Output>(value: T, fields: &[(&str, super::dimensions::FieldDimension)], d: &crate::kernel::signature::DimExpr, e: &crate::kernel::signature::DimExpr) -> Result<Value, OpError> {
    if fields.iter().any(|(_, formula)| matches!(formula, super::dimensions::FieldDimension::Unbound)) {
        return Err(OpError::Unsupported { backend: "numerical".into(), feature: "unbound structured dimensional policy".into() });
    }
    let signature = Signature::new(Vec::new(), structured_type::<T>(fields, d, e));
    runtime_dimensions::finish(&signature, &std::collections::BTreeMap::new(), value.encode()?)
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
