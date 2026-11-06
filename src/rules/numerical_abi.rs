//! Declarative host ABI for closure-free numerical kernels.
//!
//! Each entry provides its documentation, namespace and Rust function in one
//! place. `install` derives the manifest signature and dynamic dispatcher from
//! its adapter shape; adding an entry requires no host-language wrapper.
/// Structured numerical exports built from reusable argument and result codecs.
pub mod codecs;
/// Declarative dimensional relationships for numerical exports.
pub mod dimensions;
mod iterative;
pub(super) mod runtime_dimensions;
use dimensions::DimensionPolicy;

use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernels::dense;
use crate::kernels::special;
use crate::kernels::stats;
use std::sync::Arc;

type VectorBinary = fn(&[f64], &[f64]) -> Result<Vec<f64>, String>;

/// An adapter family. Each variant determines both argument validation and signature.
#[derive(Clone, Copy)]
pub enum Adapter {
    /// Scalar to scalar.
    Scalar(fn(f64) -> f64),
    /// Scalar to scalar; infinities are valid mathematical outputs, NaN is not.
    ScalarAllowInfinity(fn(f64) -> f64),
    /// Vector to scalar.
    VectorScalar(fn(&[f64]) -> f64),
    /// Fallible vector to scalar with a domain-specific constraint.
    VectorScalarFallible(fn(&[f64]) -> Result<f64, String>),
    /// A mutable vector kernel; the input is copied before invocation.
    MutableVectorScalar(fn(&mut [f64]) -> f64),
    /// Vector and scalar to scalar.
    VectorScalarScalar(fn(&[f64], f64) -> f64),
    /// Vector and scalar to vector.
    VectorScalarVector(fn(&[f64], f64) -> Vec<f64>),
    /// Vector to vector.
    VectorVector(fn(&[f64]) -> Vec<f64>),
    /// Fallible vector to vector.
    VectorVectorFallible(fn(&[f64]) -> Result<Vec<f64>, String>),
    /// Two equal-length vectors to scalar.
    VectorsScalar(fn(&[f64], &[f64]) -> Result<f64, String>),
    /// Two equal-length vectors to vector.
    VectorsVector(VectorBinary),
    /// Dense matrix to scalar.
    MatrixScalar(fn(&dense::Mat) -> Result<f64, dense::DenseError>),
    /// Matrix to scalar, allowing an infinite (but never NaN) mathematical result.
    MatrixScalarAllowInfinity(fn(&dense::Mat) -> Result<f64, dense::DenseError>),
    /// Dense matrix to matrix.
    MatrixMatrix(fn(&dense::Mat) -> Result<dense::Mat, dense::DenseError>),
    /// Two dense matrices to matrix.
    MatricesMatrix(fn(&dense::Mat, &dense::Mat) -> Result<dense::Mat, dense::DenseError>),
}

/// An auditable numerical export, consumed by manifest registration and dispatch.
pub struct Export {
    /// Public kernel function this adapter invokes.
    pub rust_path: &'static str,
    /// Unique ABI name.
    pub name: &'static str,
    /// Rule set and generated namespace.
    pub namespace: &'static str,
    /// Kernel category.
    pub category: &'static str,
    /// User-facing behavior and conventions.
    pub doc: &'static str,
    /// Stability of this interface.
    pub stability: Stability,
    /// Whether an infinite scalar result is a documented outcome. NaN is always rejected.
    pub allows_infinite_output: bool,
    /// Shared argument/result dimensional relationship.
    pub dimension_policy: DimensionPolicy,
    /// Signature and implementation family.
    pub adapter: Adapter,
}

macro_rules! export {
    ($name:literal, $namespace:literal, $doc:literal, $adapter:expr, $rust_path:literal, $policy:expr) => {
        Export {
            rust_path: $rust_path,
            name: $name,
            namespace: $namespace,
            category: $namespace,
            doc: $doc,
            stability: Stability::Experimental,
            allows_infinite_output: false,
            dimension_policy: $policy,
            adapter: $adapter,
        }
    };
}
impl Export {
    /// The exact signature derived from adapter shape and dimensional policy.
    #[must_use]
    pub fn signature(&self) -> Signature {
        dimensions::apply(self.adapter.signature(), self.dimension_policy)
    }
}
impl EvalDyn for &'static Export {
    fn signature(&self) -> Signature {
        Export::signature(self)
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let arity = self.signature().params.len();
        if args.len() != arity {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(arity as u8),
                got: args.len(),
            });
        }
        let dimensions = runtime_dimensions::prepare(self.dimension_policy, args)?;
        let result = self.adapter.eval(args, ctx)?;
        runtime_dimensions::finish(&self.signature(), &dimensions, result)
    }
}
include!("numerical_abi/generated.rs");

define_numerical_exports!(EXPORTS;
    export!(
        "num_entropy",
        "numerical_stats",
        "Shannon entropy of the supplied probabilities.",
        Adapter::VectorScalarFallible(entropy),
        "crate::kernels::stats::shannon_entropy",
        DimensionPolicy::Dimensionless
    ),
    export!(
        "num_geometric_mean",
        "numerical_stats",
        "Geometric mean of a sample.",
        Adapter::VectorScalarFallible(geometric_mean),
        "crate::kernels::stats::geometric_mean",
        DimensionPolicy::Preserve
    ),
    export!(
        "num_harmonic_mean",
        "numerical_stats",
        "Harmonic mean of a sample.",
        Adapter::VectorScalarFallible(harmonic_mean),
        "crate::kernels::stats::harmonic_mean",
        DimensionPolicy::Preserve
    ),
    export!(
        "num_covariance",
        "numerical_stats",
        "Sample covariance of equally sized samples.",
        Adapter::VectorsScalar(covariance),
        "crate::kernels::stats::covariance",
        DimensionPolicy::Product
    ),
    export!(
        "num_correlation",
        "numerical_stats",
        "Pearson correlation of equally sized samples.",
        Adapter::VectorsScalar(correlation),
        "crate::kernels::stats::correlation",
        DimensionPolicy::IndependentUnitless
    ),
    Export {
        rust_path: "crate::kernels::dense::cond",
        name: "num_cond",
        namespace: "numerical_matrix",
        category: "numerical_matrix",
        doc: "2-norm condition number; infinity for a singular matrix.",
        stability: Stability::Experimental,
        allows_infinite_output: true,
        dimension_policy: DimensionPolicy::UnitlessResult,
        adapter: Adapter::MatrixScalarAllowInfinity(dense::cond),
    },
    export!(
        "num_spectral_radius",
        "numerical_matrix",
        "Largest absolute eigenvalue.",
        Adapter::MatrixScalar(dense::spectral_radius),
        "crate::kernels::dense::spectral_radius",
        DimensionPolicy::Preserve
    ),
    export!(
        "num_pinv",
        "numerical_matrix",
        "Moore-Penrose pseudoinverse.",
        Adapter::MatrixMatrix(dense::pinv),
        "crate::kernels::dense::pinv",
        DimensionPolicy::Inverse
    ),
    export!(
        "num_cholesky",
        "numerical_matrix",
        "Lower Cholesky factor of a positive definite matrix.",
        Adapter::MatrixMatrix(dense::cholesky),
        "crate::kernels::dense::cholesky",
        DimensionPolicy::Root
    ),
    export!(
        "num_transpose",
        "numerical_matrix",
        "Transpose a dense matrix.",
        Adapter::MatrixMatrix(transpose),
        "crate::kernels::dense::Mat::transpose",
        DimensionPolicy::Preserve
    ),
    export!(
        "num_matmul",
        "numerical_matrix",
        "Dense matrix product.",
        Adapter::MatricesMatrix(dense::Mat::matmul),
        "crate::kernels::dense::Mat::matmul",
        DimensionPolicy::Product
    ),
);

define_numerical_exports!(GENERATED_EXPORTS;);

/// All descriptor-based numerical exports, including generated entries.
pub fn exports() -> impl Iterator<Item = &'static Export> {
    EXPORTS.iter()
}

fn entropy(p: &[f64]) -> Result<f64, String> {
    if p.iter().any(|v| *v < 0.0 || *v > 1.0) || (p.iter().sum::<f64>() - 1.0).abs() > 1e-10 {
        return Err("probabilities must lie in [0, 1] and sum to one".into());
    }
    Ok(stats::shannon_entropy(p))
}
fn geometric_mean(data: &[f64]) -> Result<f64, String> {
    if data.is_empty() || data.iter().any(|v| *v <= 0.0) {
        return Err("geometric mean requires a non-empty positive sample".into());
    }
    Ok(stats::geometric_mean(data))
}
fn harmonic_mean(data: &[f64]) -> Result<f64, String> {
    if data.is_empty() || data.iter().any(|v| *v <= 0.0) {
        return Err("harmonic mean requires a non-empty positive sample".into());
    }
    Ok(stats::harmonic_mean(data))
}
// These infallible functions implement fallible adapter slots shared with
// other kernels in the same family.
#[allow(clippy::unnecessary_wraps)]
fn covariance(
    a: &[f64],
    b: &[f64],
) -> Result<f64, String> {
    Ok(stats::covariance(a, b))
}
#[allow(clippy::unnecessary_wraps)]
fn correlation(
    a: &[f64],
    b: &[f64],
) -> Result<f64, String> {
    Ok(stats::correlation(a, b))
}
#[allow(clippy::unnecessary_wraps)]
fn transpose(a: &dense::Mat) -> Result<dense::Mat, dense::DenseError> {
    Ok(a.transpose())
}
fn tensor_type(rank: u64) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(TypeExpr::Float),
    }
}
impl Adapter {
    fn signature(self) -> Signature {
        use Adapter::MatricesMatrix;
        use Adapter::MatrixMatrix;
        use Adapter::MatrixScalar;
        use Adapter::MatrixScalarAllowInfinity;
        use Adapter::MutableVectorScalar;
        use Adapter::Scalar;
        use Adapter::ScalarAllowInfinity;
        use Adapter::VectorScalar;
        use Adapter::VectorScalarFallible;
        use Adapter::VectorScalarScalar;
        use Adapter::VectorScalarVector;
        use Adapter::VectorVector;
        use Adapter::VectorVectorFallible;
        use Adapter::VectorsScalar;
        use Adapter::VectorsVector;
        match self {
            | Scalar(_) | ScalarAllowInfinity(_) => {
                Signature::new(vec![TypeExpr::Float], TypeExpr::Float)
            },
            | VectorScalar(_) | VectorScalarFallible(_) | MutableVectorScalar(_) => {
                Signature::new(vec![tensor_type(1)], TypeExpr::Float)
            },
            | VectorScalarScalar(_) => {
                Signature::new(vec![tensor_type(1), TypeExpr::Float], TypeExpr::Float)
            },
            | VectorScalarVector(_) => {
                Signature::new(vec![tensor_type(1), TypeExpr::Float], tensor_type(1))
            },
            | VectorVector(_) | VectorVectorFallible(_) => {
                Signature::new(vec![tensor_type(1)], tensor_type(1))
            },
            | VectorsScalar(_) => {
                Signature::new(vec![tensor_type(1), tensor_type(1)], TypeExpr::Float)
            },
            | VectorsVector(_) => {
                Signature::new(vec![tensor_type(1), tensor_type(1)], tensor_type(1))
            },
            | MatrixScalar(_) | MatrixScalarAllowInfinity(_) => {
                Signature::new(vec![tensor_type(2)], TypeExpr::Float)
            },
            | MatrixMatrix(_) => Signature::new(vec![tensor_type(2)], tensor_type(2)),
            | MatricesMatrix(_) => {
                Signature::new(vec![tensor_type(2), tensor_type(2)], tensor_type(2))
            },
        }
    }
}
pub(super) fn numeric(
    v: &Value,
    arg: usize,
) -> Result<f64, OpError> {
    use num_traits::ToPrimitive;
    use num_traits::Zero;
    let value = match v {
        | Value::Dyn(dynamic) => numeric(&dynamic.value, arg)?,
        | Value::Quantity(quantity) => {
            if quantity.unit.affine || quantity.unit.offset.is_some() {
                return Err(OpError::Domain(format!(
                    "argument {arg} is an affine point, not a linear quantity"
                )));
            }
            if quantity.unit.scale.denom.is_zero() {
                return Err(OpError::Domain(format!(
                    "argument {arg} has an invalid unit scale"
                )));
            }
            let scale = quantity.unit.scale.rational().to_f64().ok_or_else(|| {
                OpError::Domain(format!("argument {arg} unit scale does not fit Float"))
            })? * std::f64::consts::PI.powi(quantity.unit.scale.pi);
            if !scale.is_finite() || scale == 0.0 {
                return Err(OpError::Domain(format!(
                    "argument {arg} unit scale is not representable as a finite nonzero Float"
                )));
            }
            numeric(&quantity.value, arg)? * scale
        },
        | _ => {
            v.as_f64().ok_or_else(|| {
                OpError::Type {
                    expected: Box::new(TypeExpr::Float),
                    got: Box::new(v.type_expr()),
                    arg,
                }
            })?
        },
    };
    if !value.is_finite() {
        return Err(OpError::Domain(format!(
            "argument {arg} must be finite and representable as Float"
        )));
    }
    Ok(value)
}

fn tensor(
    v: &Value,
    rank: usize,
    arg: usize,
) -> Result<(Vec<usize>, Vec<f64>), OpError> {
    let Value::Tensor(t) = v else {
        return Err(OpError::Type {
            expected: Box::new(tensor_type(rank as u64)),
            got: Box::new(v.type_expr()),
            arg,
        });
    };
    if t.shape.len() != rank {
        return Err(OpError::Domain(format!(
            "argument {arg} must have rank {rank}, got {}",
            t.shape.len()
        )));
    }
    if t.shape.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(t.data.len()) {
        return Err(OpError::Domain(format!(
            "argument {arg} tensor shape does not match its data"
        )));
    }
    let data = t
        .data
        .iter()
        .map(|x| numeric(x, arg))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((t.shape.clone(), data))
}
fn matrix(
    v: &Value,
    arg: usize,
) -> Result<dense::Mat, OpError> {
    let (shape, data) = tensor(v, 2, arg)?;
    Ok(dense::Mat {
        rows: shape[0],
        cols: shape[1],
        data,
    })
}
fn output_scalar(x: f64) -> Result<Value, OpError> {
    if !x.is_finite() {
        return Err(OpError::Domain(
            "numerical kernel returned a non-finite result".into(),
        ));
    }
    Ok(Value::Float(x))
}
fn output_tensor(
    shape: Vec<usize>,
    data: Vec<f64>,
) -> Result<Value, OpError> {
    let data = data
        .into_iter()
        .map(output_scalar)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape,
        data,
        elem: TypeExpr::Float,
    }))
}
fn output_vector(data: Vec<f64>) -> Result<Value, OpError> {
    output_tensor(vec![data.len()], data)
}
fn output_matrix(m: dense::Mat) -> Result<Value, OpError> {
    output_tensor(vec![m.rows, m.cols], m.data)
}
fn dense_error(error: dense::DenseError) -> OpError {
    match error {
        | dense::DenseError::NoConvergence => {
            OpError::NotReduced("numerical matrix algorithm did not converge".into())
        },
        | other => OpError::Domain(other.to_string()),
    }
}
impl EvalDyn for Adapter {
    fn signature(&self) -> Signature {
        (*self).signature()
    }

    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use Adapter::MatricesMatrix;
        use Adapter::MatrixMatrix;
        use Adapter::MatrixScalar;
        use Adapter::MatrixScalarAllowInfinity;
        use Adapter::MutableVectorScalar;
        use Adapter::Scalar;
        use Adapter::ScalarAllowInfinity;
        use Adapter::VectorScalar;
        use Adapter::VectorScalarFallible;
        use Adapter::VectorScalarScalar;
        use Adapter::VectorScalarVector;
        use Adapter::VectorVector;
        use Adapter::VectorVectorFallible;
        use Adapter::VectorsScalar;
        use Adapter::VectorsVector;
        let sig = self.signature();
        if args.len() != sig.params.len() {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo {
                    min: sig.params.len() as u8,
                    max: Some(sig.params.len() as u8),
                },
                got: args.len(),
            });
        }
        match self {
            | Scalar(f) => output_scalar(f(numeric(&args[0], 0)?)),
            | ScalarAllowInfinity(f) => {
                let result = f(numeric(&args[0], 0)?);
                if result.is_nan() {
                    return Err(OpError::Domain("numerical kernel returned NaN".into()));
                }
                Ok(Value::Float(result))
            },
            | VectorScalar(f) => output_scalar(f(&tensor(&args[0], 1, 0)?.1)),
            | MutableVectorScalar(f) => output_scalar(f(&mut tensor(&args[0], 1, 0)?.1)),
            | VectorScalarScalar(f) => {
                output_scalar(f(&tensor(&args[0], 1, 0)?.1, numeric(&args[1], 1)?))
            },
            | VectorScalarVector(f) => {
                output_vector(f(&tensor(&args[0], 1, 0)?.1, numeric(&args[1], 1)?))
            },
            | VectorScalarFallible(f) => {
                output_scalar(f(&tensor(&args[0], 1, 0)?.1).map_err(OpError::Domain)?)
            },
            | VectorVector(f) => output_vector(f(&tensor(&args[0], 1, 0)?.1)),
            | VectorVectorFallible(f) => {
                output_vector(f(&tensor(&args[0], 1, 0)?.1).map_err(OpError::Domain)?)
            },
            | VectorsScalar(f) => {
                let a = tensor(&args[0], 1, 0)?.1;
                let b = tensor(&args[1], 1, 1)?.1;
                if a.len() != b.len() {
                    return Err(OpError::Domain("vector lengths differ".into()));
                }
                output_scalar(f(&a, &b).map_err(OpError::Domain)?)
            },
            | VectorsVector(f) => {
                let a = tensor(&args[0], 1, 0)?.1;
                let b = tensor(&args[1], 1, 1)?.1;
                if a.len() != b.len() {
                    return Err(OpError::Domain("vector lengths differ".into()));
                }
                output_vector(f(&a, &b).map_err(OpError::Domain)?)
            },
            | MatrixScalar(f) => output_scalar(f(&matrix(&args[0], 0)?).map_err(dense_error)?),
            | MatrixScalarAllowInfinity(f) => {
                let result = f(&matrix(&args[0], 0)?).map_err(dense_error)?;
                if result.is_nan() {
                    return Err(OpError::Domain("numerical kernel returned NaN".into()));
                }
                Ok(Value::Float(result))
            },
            | MatrixMatrix(f) => output_matrix(f(&matrix(&args[0], 0)?).map_err(dense_error)?),
            | MatricesMatrix(f) => {
                let a = matrix(&args[0], 0)?;
                let b = matrix(&args[1], 1)?;
                if a.rows.checked_mul(b.cols).is_none() {
                    return Err(OpError::Domain("matrix result shape overflows".into()));
                }
                output_matrix(f(&a, &b).map_err(dense_error)?)
            },
        }
    }
}
fn install(
    i: &mut Installer<'_>,
    namespace: &str,
) -> Result<(), RuleError> {
    for export in exports().filter(|entry| entry.namespace == namespace) {
        let signature = export.signature();
        i.op(
            OpDescriptor::new(export.name, Arity::Fixed(signature.params.len() as u8))
                .rust_path(export.rust_path)
                .category(export.category)
                .doc(export.doc)
                .stability(export.stability)
                .signature(signature)
                .eval_dyn(Arc::new(export)),
        )?;
    }
    for export in codecs::EXPORTS
        .iter()
        .copied()
        .filter(|entry| entry.namespace == namespace)
    {
        let signature = export.signature();
        i.op(
            OpDescriptor::new(export.name, Arity::Fixed(signature.params.len() as u8))
                .rust_path(export.rust_path)
                .category(export.namespace)
                .doc(export.doc)
                .stability(export.stability)
                .signature(signature)
                .eval_dyn(Arc::new(export)),
        )?;
    }
    if namespace == "numerical_matrix" {
        for export in iterative::EXPORTS {
            let signature = export.signature();
            i.op(
                OpDescriptor::new(export.name, Arity::Fixed(signature.params.len() as u8))
                    .rust_path(export.rust_path)
                    .category(namespace)
                    .doc(export.doc)
                    .stability(Stability::Experimental)
                    .signature(signature)
                    .eval_dyn(Arc::new(export)),
            )?;
        }
    }
    Ok(())
}
macro_rules! rule_set {
    ($function:ident, $install:ident, $namespace:literal) => {
        fn $install(i: &mut Installer<'_>) -> Result<(), RuleError> {
            install(i, $namespace)
        }
        /// Numerical ABI rule set.
        #[must_use]
        pub fn $function() -> RuleSet {
            RuleSet::new($namespace, $install)
        }
    };
}
rule_set!(numerical_vector, install_vector, "numerical_vector");
rule_set!(numerical_stats, install_stats, "numerical_stats");
rule_set!(numerical_special, install_special, "numerical_special");
rule_set!(numerical_matrix, install_matrix, "numerical_matrix");
