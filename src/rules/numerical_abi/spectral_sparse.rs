//! Dimensioned structural sparse eigenvalue and multigrid adapters.
use super::dimensions::DimContract;
use super::dimensions::DimFormula;
use super::dimensions::FieldDimension;
use super::dimensions::InputBinding;
use super::dimensions::InputConstraint;
use super::dimensions::{
    self,
};
use super::iterative::csr_type;
use super::iterative::decode_csr;
use super::iterative::integer;
use super::iterative::relative_residual;
use super::numeric;
use super::tensor;
use super::tensor_type;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::RecordValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernels::amg_general;
use crate::kernels::eigs;
use crate::kernels::krylov;
use crate::kernels::multigrid;
use std::sync::Arc;

static EIG_FIELDS: &[(&str, FieldDimension)] = &[
    ("values_re", FieldDimension::Formula(&DimFormula::Var("D"))),
    ("values_im", FieldDimension::Formula(&DimFormula::Var("D"))),
    ("vectors_re", FieldDimension::Formula(&DimFormula::One)),
    ("vectors_im", FieldDimension::Formula(&DimFormula::One)),
    ("residuals", FieldDimension::Formula(&DimFormula::Var("D"))),
    ("converged", FieldDimension::NonNumeric),
    ("matvecs", FieldDimension::NonNumeric),
    ("restarts", FieldDimension::NonNumeric),
];
static EIG_CONTRACT: DimContract = DimContract {
    inputs: &[InputBinding {
        name: "D",
        argument: 0,
        path: &["values"],
    }],
    constraints: &[],
    outputs: EIG_FIELDS,
    output_paths: &[],
    result: None,
    other_args_dimensionless: true,
};
static SHIFT_CONTRACT: DimContract = DimContract {
    inputs: &[InputBinding {
        name: "D",
        argument: 0,
        path: &["values"],
    }],
    constraints: &[InputConstraint {
        argument: 2,
        path: &[],
        formula: &DimFormula::Var("D"),
    }],
    outputs: EIG_FIELDS,
    output_paths: &[],
    result: None,
    other_args_dimensionless: true,
};
static AMG_CONTRACT: DimContract = DimContract {
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
    constraints: &[],
    outputs: &[
        (
            "x",
            FieldDimension::Formula(&DimFormula::Mul(&[
                DimFormula::Var("E"),
                DimFormula::Pow(&DimFormula::Var("D"), -1, 1),
            ])),
        ),
        ("residual", FieldDimension::Formula(&DimFormula::One)),
        ("history", FieldDimension::Formula(&DimFormula::One)),
        ("iterations", FieldDimension::NonNumeric),
        ("converged", FieldDimension::NonNumeric),
    ],
    output_paths: &[],
    result: None,
    other_args_dimensionless: true,
};
#[derive(Clone, Copy)]
enum Method {
    Lanczos,
    Arnoldi,
    ShiftInvert,
    Amg,
    GeneralAmg,
    CgAmg,
}
fn output_type(method: Method) -> TypeExpr {
    match method {
        | Method::Lanczos | Method::Arnoldi | Method::ShiftInvert => {
            TypeExpr::Record(vec![
                ("values_re".into(), tensor_type(1)),
                ("values_im".into(), tensor_type(1)),
                ("vectors_re".into(), tensor_type(2)),
                ("vectors_im".into(), tensor_type(2)),
                ("residuals".into(), tensor_type(1)),
                ("converged".into(), TypeExpr::Bool),
                ("matvecs".into(), TypeExpr::Int),
                ("restarts".into(), TypeExpr::Int),
            ])
        },
        | _ => {
            TypeExpr::Record(vec![
                ("x".into(), tensor_type(1)),
                ("residual".into(), TypeExpr::Float),
                ("history".into(), tensor_type(1)),
                ("iterations".into(), TypeExpr::Int),
                ("converged".into(), TypeExpr::Bool),
            ])
        },
    }
}
fn contract(method: Method) -> &'static DimContract {
    match method {
        | Method::ShiftInvert => &SHIFT_CONTRACT,
        | Method::Lanczos | Method::Arnoldi => &EIG_CONTRACT,
        | _ => &AMG_CONTRACT,
    }
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let params = match self {
            | Self::Lanczos | Self::Arnoldi => {
                vec![
                    csr_type(),
                    TypeExpr::Int,
                    TypeExpr::String,
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | Self::ShiftInvert => {
                vec![
                    csr_type(),
                    TypeExpr::Int,
                    TypeExpr::Float,
                    TypeExpr::Bool,
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | _ => vec![csr_type(), tensor_type(1), TypeExpr::Float, TypeExpr::Int],
        };
        dimensions::apply_contract(Signature::new(params, output_type(*self)), contract(*self))
    }

    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        let count = signature.params.len();
        if args.len() != count {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(count as u8),
                got: args.len(),
            });
        }
        let dims = super::runtime_dimensions::prepare_contract(contract(*self), args)?;
        let a = decode_csr(&args[0])?;
        let result = match self {
            | Self::Lanczos | Self::Arnoldi | Self::ShiftInvert => {
                let k = integer(&args[1], "k")?;
                if k == 0 || k > a.rows || matches!(self, Self::Arnoldi) && k == a.rows {
                    return Err(OpError::Domain(
                        "eigenpair count is outside the method's matrix dimension range".into(),
                    ));
                }
                let (ti, mi, ni) = if matches!(self, Self::ShiftInvert) {
                    (4, 5, 6)
                } else {
                    (3, 4, 5)
                };
                let tol = numeric(&args[ti], ti)?;
                if !tol.is_finite() || tol <= 0.0 {
                    return Err(OpError::Domain(
                        "sparse eigen tolerance must be finite and positive".into(),
                    ));
                }
                let max_restarts = integer(&args[mi], "max_restarts")?;
                let ncv = integer(&args[ni], "ncv")?;
                if max_restarts == 0 || (ncv != 0 && (ncv <= k || ncv > a.rows)) {
                    return Err(OpError::Domain(
                        "positive restart budget and ncv=0 or k<ncv<=rows are required".into(),
                    ));
                }
                let opts = eigs::EigsOptions {
                    ncv,
                    tol,
                    max_restarts,
                    seed: 0x5eed,
                };
                let eigen = if matches!(self, Self::ShiftInvert) {
                    let Value::Bool(symmetric) = &args[3] else {
                        return Err(OpError::Domain(
                            "shift-invert symmetry must be Boolean".into(),
                        ));
                    };
                    if *symmetric {
                        require_symmetric(&a)?;
                    }
                    let sigma = numeric(&args[2], 2)?;
                    eigs::eigs_csr_shift_invert(&a, k, sigma, *symmetric, &opts).map_err(|e| {
                        OpError::NotReduced(format!("sparse shift-invert failed: {e:?}"))
                    })?
                } else {
                    let Value::String(which) = &args[2] else {
                        return Err(OpError::Domain(
                            "sparse eigen selection must be text".into(),
                        ));
                    };
                    let which = match which.as_str() {
                        | "largest_magnitude" => eigs::Which::LargestMagnitude,
                        | "smallest_magnitude" => eigs::Which::SmallestMagnitude,
                        | "largest_real" => eigs::Which::LargestReal,
                        | "smallest_real" => eigs::Which::SmallestReal,
                        | "largest_imag" => eigs::Which::LargestImag,
                        | _ => {
                            return Err(OpError::Domain(
                                "unknown sparse eigen spectrum selection".into(),
                            ));
                        },
                    };
                    if matches!(self, Self::Lanczos) {
                        require_symmetric(&a)?;
                        eigs::eigs_sym(|x, y| a.matvec(x, y), a.rows, k, which, &opts)
                    } else {
                        eigs::eigs(|x, y| a.matvec(x, y), a.rows, k, which, &opts)
                    }
                    .map_err(|e| OpError::NotReduced(format!("sparse eigensolver failed: {e:?}")))?
                };
                encode_eigen(&a, eigen, k, tol)?
            },
            | _ => {
                let (_, b) = tensor(&args[1], 1, 1)?;
                if b.len() != a.rows {
                    return Err(OpError::Domain("AMG matrix and RHS shapes disagree".into()));
                }
                let tol = numeric(&args[2], 2)?;
                let budget = integer(&args[3], "iteration/cycle budget")?;
                if !tol.is_finite() || tol <= 0.0 || budget == 0 {
                    return Err(OpError::Domain(
                        "AMG tolerance and cycle budget must be positive".into(),
                    ));
                }
                let (x, iterations, history, converged) = match self {
                    | Self::GeneralAmg => {
                        let hierarchy = amg_general::AmgGeneral::build(
                            &a,
                            &amg_general::AmgGeneralOptions::default(),
                        )
                        .map_err(|e| OpError::Domain(format!("AMG setup: {e:?}")))?;
                        let r = hierarchy.solve(&b, None, tol, budget);
                        (
                            r.x,
                            r.history.len().saturating_sub(1),
                            r.history,
                            r.converged,
                        )
                    },
                    | Self::Amg | Self::CgAmg => {
                        require_symmetric(&a)?;
                        let hierarchy =
                            multigrid::Amg::build(&a, &multigrid::AmgOptions::default())
                                .map_err(|e| OpError::Domain(format!("AMG setup: {e:?}")))?;
                        if matches!(self, Self::CgAmg) {
                            let r = krylov::cg(
                                |x, y| a.matvec(x, y),
                                |r, z| hierarchy.precondition(r, z),
                                &b,
                                None,
                                tol,
                                budget,
                            );
                            (r.x, r.iterations, vec![r.residual], r.converged)
                        } else {
                            let r = hierarchy.solve(&b, None, tol, budget);
                            (
                                r.x,
                                r.history.len().saturating_sub(1),
                                r.history,
                                r.converged,
                            )
                        }
                    },
                    | _ => unreachable!(),
                };
                let residual = relative_residual(&a, &x, &b)?;
                Value::Record(RecordValue {
                    fields: vec![
                        ("x".into(), output(vec![x.len()], x)?),
                        ("residual".into(), Value::Float(residual)),
                        ("history".into(), output(vec![history.len()], history)?),
                        ("iterations".into(), Value::int(iterations)),
                        (
                            "converged".into(),
                            Value::Bool(converged && residual <= tol),
                        ),
                    ],
                })
            },
        };
        super::runtime_dimensions::finish(&signature, &dims, result)
    }
}
fn require_symmetric(a: &krylov::Csr) -> Result<(), OpError> {
    for i in 0..a.rows {
        for p in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[p];
            let start = a.indptr[j];
            let end = a.indptr[j + 1];
            let other = a.indices[start..end]
                .binary_search(&i)
                .ok()
                .map_or(0.0, |offset| a.values[start + offset]);
            let value = a.values[p];
            if (value - other).abs() > 32.0 * f64::EPSILON * value.abs().max(other.abs()) {
                return Err(OpError::Domain(
                    "symmetric sparse method requires numerically symmetric coefficients".into(),
                ));
            }
        }
    }
    Ok(())
}
fn output(
    shape: Vec<usize>,
    data: Vec<f64>,
) -> Result<Value, OpError> {
    if shape.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(data.len())
        || data.iter().any(|v| !v.is_finite())
    {
        return Err(OpError::Domain(
            "sparse solver returned nonfinite or malformed result".into(),
        ));
    }
    Ok(Value::Tensor(TensorValue {
        shape,
        data: data.into_iter().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    }))
}
fn encode_eigen(
    a: &krylov::Csr,
    mut result: eigs::EigsResult,
    k: usize,
    tol: f64,
) -> Result<Value, OpError> {
    let count = result.values.len();
    if count != result.vectors_re.len() || count != result.vectors_im.len() {
        return Err(OpError::Domain(
            "sparse eigensolver returned inconsistent vector counts".into(),
        ));
    }
    let scale = a.values.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    let mut relative = Vec::new();
    let mut residuals = Vec::new();
    for i in 0..count {
        let xr = &mut result.vectors_re[i];
        let xi = &mut result.vectors_im[i];
        if xr.len() != a.rows || xi.len() != a.rows {
            return Err(OpError::Domain("malformed sparse eigenvector".into()));
        }
        let norm = xr.iter().chain(xi.iter()).fold(0.0_f64, |n, v| n.hypot(*v));
        if !norm.is_finite() || norm <= 0.0 {
            return Err(OpError::Domain(
                "sparse eigenvector has invalid norm".into(),
            ));
        }
        for x in xr.iter_mut().chain(xi.iter_mut()) {
            *x /= norm;
        }
        let mut ar = vec![0.0; a.rows];
        let mut ai = vec![0.0; a.rows];
        a.matvec(xr, &mut ar);
        a.matvec(xi, &mut ai);
        let (lr, li) = result.values[i];
        let r = (0..a.rows).fold(0.0_f64, |n, j| {
            n.hypot(ar[j] - lr * xr[j] + li * xi[j])
                .hypot(ai[j] - lr * xi[j] - li * xr[j])
        });
        if !r.is_finite() {
            return Err(OpError::Domain(
                "sparse eigen residual exceeds finite range".into(),
            ));
        }
        residuals.push(r);
        let reference = scale.max(lr.hypot(li));
        relative.push(if reference > 0.0 {
            r / reference
        } else {
            r
        });
    }
    let mut vr = Vec::new();
    let mut vi = Vec::new();
    for row in 0..a.rows {
        for column in 0..count {
            vr.push(result.vectors_re[column][row]);
            vi.push(result.vectors_im[column][row]);
        }
    }
    let values_re = result.values.iter().map(|v| v.0).collect();
    let values_im = result.values.iter().map(|v| v.1).collect();
    Ok(Value::Record(RecordValue {
        fields: vec![
            ("values_re".into(), output(vec![count], values_re)?),
            ("values_im".into(), output(vec![count], values_im)?),
            ("vectors_re".into(), output(vec![a.rows, count], vr)?),
            ("vectors_im".into(), output(vec![a.rows, count], vi)?),
            ("residuals".into(), output(vec![count], residuals)?),
            (
                "converged".into(),
                Value::Bool(result.converged && count == k && relative.iter().all(|v| *v <= tol)),
            ),
            ("matvecs".into(), Value::int(result.matvecs)),
            ("restarts".into(), Value::int(result.restarts)),
        ],
    }))
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        (
            "num_csr_eigs_sym",
            "crate::kernels::eigs::eigs_sym",
            Method::Lanczos,
        ),
        (
            "num_csr_eigs",
            "crate::kernels::eigs::eigs",
            Method::Arnoldi,
        ),
        (
            "num_csr_eigs_shift_invert",
            "crate::kernels::eigs::eigs_csr_shift_invert",
            Method::ShiftInvert,
        ),
        (
            "num_csr_amg_solve",
            "crate::kernels::multigrid::Amg::solve",
            Method::Amg,
        ),
        (
            "num_csr_amg_general_solve",
            "crate::kernels::amg_general::AmgGeneral::solve",
            Method::GeneralAmg,
        ),
        (
            "num_csr_cg_amg",
            "crate::kernels::multigrid::Amg::precondition",
            Method::CgAmg,
        ),
    ] {
        let signature = method.signature();
        let op=installer.op(OpDescriptor::new(name,Arity::Fixed(signature.params.len()as u8)).rust_path(path).category("numerical_matrix")
        .doc("Validated sorted structural CSR sparse eigen/AMG adapters. Lanczos/Arnoldi: matrix,k,selection text,positive relative tolerance,max_restarts,ncv (0 auto). Shift-invert: matrix,k,shift,symmetric,relative tolerance,max_restarts,ncv. Eigenvalues and independently recomputed original-system residual norms carry matrix dimensions; eigenvectors are normalized and dimensionless. Symmetric methods reject nonsymmetric coefficients. AMG: matrix,RHS,positive relative tolerance,cycle/iteration budget; solution has RHS/matrix dimensions, history/residual are relative and dimensionless. Smoothed aggregation requires SPD-compatible setup; general AMG uses its positive-diagonal setup contract. Convergence flags require verified original-system residuals.")
        .stability(Stability::Experimental).signature(signature).eval_dyn(Arc::new(method)))?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
