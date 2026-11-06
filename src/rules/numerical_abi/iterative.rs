//! Validated structural CSR adapters for iterative solvers.
use super::codecs::Output;
use super::dimensions::DimContract;
use super::dimensions::DimFormula;
use super::dimensions::InputBinding;
use super::dimensions::InputConstraint;
use super::numeric;
use super::tensor;
use super::tensor_type;
use crate::graph::op::EvalDyn;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::Value;
use crate::kernels::krylov::Csr;
use crate::kernels::krylov::IterResult;
use crate::kernels::krylov::{self};
use num_traits::ToPrimitive;

#[derive(Clone, Copy)]
enum Method {
    Cg,
    Gmres,
    Bicgstab,
}

/// Reviewed solver semantics; the Rust kernel supplies the actual algorithm.
pub(super) struct IterativeExport {
    pub(super) name: &'static str,
    pub(super) rust_path: &'static str,
    pub(super) doc: &'static str,
    method: Method,
}
pub(super) static EXPORTS: &[IterativeExport] = &[
    IterativeExport {
        name: "num_csr_cg",
        rust_path: "crate::kernels::krylov::cg",
        doc: "Conjugate gradients with identity preconditioning; the matrix must be symmetric positive definite. Convergence reports the independently recomputed relative residual.",
        method: Method::Cg,
    },
    IterativeExport {
        name: "num_csr_gmres",
        rust_path: "crate::kernels::krylov::gmres",
        doc: "Restarted GMRES with identity preconditioning. Convergence reports the independently recomputed relative residual.",
        method: Method::Gmres,
    },
    IterativeExport {
        name: "num_csr_bicgstab",
        rust_path: "crate::kernels::krylov::bicgstab",
        doc: "BiCGSTAB with identity preconditioning. Convergence reports the independently recomputed relative residual.",
        method: Method::Bicgstab,
    },
];
static ITERATIVE_CONTRACT: DimContract = DimContract {
    inputs: &[
        InputBinding {
            name: "D",
            argument: 0,
            path: &["values"],
        },
        InputBinding {
            name: "E",
            argument: 1,
            path: &[],
        },
    ],
    constraints: &[InputConstraint {
        argument: 2,
        path: &[],
        formula: &DimFormula::One,
    }],
    outputs: super::codecs::STRUCTURED_ITERATIVE_LINEAR_SOLVE,
    output_paths: &[],
    result: None,
    other_args_dimensionless: true,
};
fn csr_type() -> TypeExpr {
    let vector = |ty| TypeExpr::Tensor {
        rank: crate::kernel::signature::NatExpr::Lit(1),
        shape: None,
        elem: Box::new(ty),
    };
    TypeExpr::Record(vec![
        ("rows".into(), TypeExpr::Int),
        ("cols".into(), TypeExpr::Int),
        ("indptr".into(), vector(TypeExpr::Int)),
        ("indices".into(), vector(TypeExpr::Int)),
        ("values".into(), vector(TypeExpr::Float)),
    ])
}
impl IterativeExport {
    pub(super) fn signature(&self) -> Signature {
        let mut params = vec![csr_type(), tensor_type(1), TypeExpr::Float, TypeExpr::Int];
        if matches!(self.method, Method::Gmres) {
            params.push(TypeExpr::Int);
        }
        super::dimensions::apply_contract(
            Signature::new(params, IterResult::ty()),
            &ITERATIVE_CONTRACT,
        )
    }
}
fn integer(
    value: &Value,
    label: &str,
) -> Result<usize, OpError> {
    match value {
        | Value::Int(n) => n.to_usize().ok_or_else(|| {
            OpError::Domain(format!("{label} must fit a nonnegative platform integer"))
        }),
        | _ => Err(OpError::Domain(format!("{label} must be an integer"))),
    }
}
fn index_vector(
    value: &Value,
    label: &str,
) -> Result<Vec<usize>, OpError> {
    let Value::Tensor(t) = value else {
        return Err(OpError::Domain(format!(
            "{label} must be an integer vector"
        )));
    };
    if t.shape.len() != 1 || t.shape[0] != t.data.len() {
        return Err(OpError::Domain(format!(
            "{label} has malformed vector shape"
        )));
    }
    t.data.iter().map(|item| integer(item, label)).collect()
}
fn decode_csr(value: &Value) -> Result<Csr, OpError> {
    let Value::Record(record) = value else {
        return Err(OpError::Domain("CSR matrix must be a record".into()));
    };
    if record.fields.len() != 5 {
        return Err(OpError::Domain(
            "CSR record must contain exactly five fields".into(),
        ));
    }
    let get = |name: &str| -> Result<&Value, OpError> {
        let mut matches = record.fields.iter().filter(|(label, _)| label == name);
        let value = matches
            .next()
            .ok_or_else(|| OpError::Domain(format!("CSR field {name} missing")))?;
        if matches.next().is_some() {
            return Err(OpError::Domain(format!("duplicate CSR field {name}")));
        }
        Ok(&value.1)
    };
    let rows = integer(get("rows")?, "rows")?;
    let cols = integer(get("cols")?, "cols")?;
    let indptr = index_vector(get("indptr")?, "indptr")?;
    let indices = index_vector(get("indices")?, "indices")?;
    let data = get("values")?;
    let (_, values) = tensor(data, 1, 0)?;
    if rows == 0
        || rows != cols
        || rows.checked_add(1) != Some(indptr.len())
        || indices.len() != values.len()
        || indptr.first() != Some(&0)
        || indptr.last() != Some(&values.len())
        || indptr
            .windows(2)
            .any(|p| p[0] > p[1] || p[1] > values.len())
    {
        return Err(OpError::Domain(
            "invalid square CSR shape or row pointers".into(),
        ));
    }
    for row in indptr.windows(2) {
        let columns = &indices[row[0]..row[1]];
        if columns.iter().any(|&column| column >= cols)
            || columns.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(OpError::Domain(
                "CSR columns must be in bounds and strictly increasing per row".into(),
            ));
        }
    }
    Ok(Csr {
        rows,
        cols,
        indptr,
        indices,
        values,
    })
}
// Reuse the solver's scale-safe original-system residual, including a positive
// subnormal sentinel when a nonzero ratio is below f64's representable range.
fn relative_residual(
    a: &Csr,
    x: &[f64],
    b: &[f64],
) -> Result<f64, OpError> {
    let mut ax = vec![0.0; a.rows];
    a.matvec(x, &mut ax);
    let ratio = krylov::residual_ratio(b, &ax)
        .ok_or_else(|| OpError::Domain("nonfinite recomputed residual".into()))?;
    if ratio.is_finite() {
        Ok(ratio)
    } else {
        Err(OpError::Domain(
            "relative residual is undefined or exceeds the finite range".into(),
        ))
    }
}

impl EvalDyn for &'static IterativeExport {
    fn signature(&self) -> Signature {
        IterativeExport::signature(self)
    }

    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = if matches!(self.method, Method::Gmres) {
            5
        } else {
            4
        };
        if args.len() != expected {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(expected as u8),
                got: args.len(),
            });
        }
        let dims = super::runtime_dimensions::prepare_contract(&ITERATIVE_CONTRACT, args)?;
        let matrix = decode_csr(&args[0])?;
        let (_, b) = tensor(&args[1], 1, 1)?;
        if b.len() != matrix.rows {
            return Err(OpError::Domain("CSR matrix and RHS shape mismatch".into()));
        }
        let tol = numeric(&args[2], 2)?;
        if tol < 0.0 {
            return Err(OpError::Domain(
                "relative tolerance must be nonnegative".into(),
            ));
        }
        let max_iter = integer(&args[3], "max_iter")?;
        if max_iter == 0 {
            return Err(OpError::Domain("max_iter must be positive".into()));
        }
        let apply = |x: &[f64], y: &mut [f64]| matrix.matvec(x, y);
        let identity = |r: &[f64], z: &mut [f64]| z.copy_from_slice(r);
        let mut result = match self.method {
            | Method::Cg => krylov::cg(apply, identity, &b, None, tol, max_iter),
            | Method::Bicgstab => krylov::bicgstab(apply, identity, &b, None, tol, max_iter),
            | Method::Gmres => {
                let restart = integer(&args[4], "restart")?;
                if restart == 0 || restart > matrix.rows {
                    return Err(OpError::Domain("restart must be in 1..=rows".into()));
                }
                krylov::gmres(apply, identity, &b, None, restart, tol, max_iter)
            },
        };
        if result.x.iter().any(|x| !x.is_finite()) {
            return Err(OpError::Domain(
                "iterative solver returned nonfinite state".into(),
            ));
        }
        result.residual = relative_residual(&matrix, &result.x, &b)?;
        result.converged = result.converged && result.residual <= tol;
        super::runtime_dimensions::finish(&self.signature(), &dims, result.encode()?)
    }
}
