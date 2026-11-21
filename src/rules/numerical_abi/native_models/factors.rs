//! Portable defining matrices reconstruct native dense and sparse factors.
use super::Arc;
use super::Arity;
use super::ArityInfo;
use super::CallCtx;
use super::D;
use super::DimContract;
use super::DimFormula;
use super::E;
use super::EvalDyn;
use super::INVERSE;
use super::Installer;
use super::LIMIT;
use super::OpDescriptor;
use super::OpError;
use super::PRODUCT;
use super::RuleError;
use super::SIDE;
use super::Signature;
use super::Stability;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::array;
use super::contract;
use super::dimensions;
use super::encode_array;
use super::error;
use super::field;
use super::input;
use super::integer;
use super::kind;
use super::output;
use super::record;
use super::runtime_dimensions;
use super::tensor_type;
use crate::kernels::dense;
use crate::kernels::krylov;
fn sparse_error(value: krylov::SparseLuError) -> OpError {
    error(format!("sparse LU factorization failed: {value:?}"))
}
const CSR_FIELDS: &[&str] = &["rows", "cols", "indptr", "indices", "values"];
const FACTOR_FIELDS: &[&str] = &["format", "matrix", "kind"];
fn indices(values: Vec<usize>) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Int,
        data: values.into_iter().map(Value::int).collect(),
    })
}
fn decode_indices(value: &Value) -> Result<Vec<usize>, OpError> {
    let Value::Tensor(t) = value else {
        return Err(error("expected integer vector"));
    };
    if t.shape != vec![t.data.len()] || t.data.len() > LIMIT {
        return Err(error("invalid integer vector shape/budget"));
    }
    t.data.iter().map(integer).collect()
}
fn csr(value: &Value) -> Result<krylov::Csr, OpError> {
    let rows = integer(field(value, "rows", CSR_FIELDS)?)?;
    let cols = integer(field(value, "cols", CSR_FIELDS)?)?;
    if rows > SIDE || cols > SIDE {
        return Err(error("portable CSR exceeds matrix-order budget"));
    }
    let Value::Tensor(t) = field(value, "values", CSR_FIELDS)? else {
        return Err(error("CSR values must be a tensor"));
    };
    if t.data.len() > LIMIT {
        return Err(error("portable CSR exceeds nonzero budget"));
    }
    let indptr = decode_indices(field(value, "indptr", CSR_FIELDS)?)?;
    let indices = decode_indices(field(value, "indices", CSR_FIELDS)?)?;
    let (_, values) = super::super::tensor(field(value, "values", CSR_FIELDS)?, 1, 0)?;
    if rows.checked_add(1) != Some(indptr.len())
        || indices.len() != values.len()
        || indptr.first() != Some(&0)
        || indptr.last() != Some(&values.len())
        || indptr
            .windows(2)
            .any(|w| w[0] > w[1] || w[1] > values.len())
    {
        return Err(error("invalid CSR shape or row offsets"));
    }
    for pair in indptr.windows(2) {
        let row = &indices[pair[0]..pair[1]];
        if row.iter().any(|column| *column >= cols) || row.windows(2).any(|w| w[0] >= w[1]) {
            return Err(error(
                "CSR column indices must be bounded and strictly increasing per row",
            ));
        }
    }
    Ok(krylov::Csr {
        rows,
        cols,
        indptr,
        indices,
        values,
    })
}
fn encode_csr(a: krylov::Csr) -> Result<Value, OpError> {
    Ok(record(vec![
        ("rows", Value::int(a.rows)),
        ("cols", Value::int(a.cols)),
        ("indptr", indices(a.indptr)),
        ("indices", indices(a.indices)),
        ("values", encode_array(vec![a.values.len()], a.values)?),
    ]))
}
fn factor_type(sparse: bool) -> TypeExpr {
    TypeExpr::Record(vec![
        ("format".into(), TypeExpr::String),
        (
            "matrix".into(),
            if sparse {
                super::super::iterative::csr_type()
            } else {
                tensor_type(2, TypeExpr::Float)
            },
        ),
        ("kind".into(), TypeExpr::Int),
    ])
}
fn encode_factor(
    matrix: Value,
    kind: usize,
    sparse: bool,
) -> Value {
    record(vec![
        (
            "format",
            Value::String(
                if sparse {
                    "rssn/sparse-factor/v1"
                } else {
                    "rssn/dense-lu/v1"
                }
                .into(),
            ),
        ),
        ("matrix", matrix),
        ("kind", Value::int(kind)),
    ])
}
fn factor(
    value: &Value,
    sparse: bool,
) -> Result<(&Value, usize), OpError> {
    kind(
        field(value, "format", FACTOR_FIELDS)?,
        if sparse {
            "rssn/sparse-factor/v1"
        } else {
            "rssn/dense-lu/v1"
        },
    )?;
    let kind = integer(field(value, "kind", FACTOR_FIELDS)?)?;
    if kind > 4 || !sparse && kind != 0 {
        return Err(error("invalid native factor kind"));
    }
    Ok((field(value, "matrix", FACTOR_FIELDS)?, kind))
}
fn dense_matrix(value: &Value) -> Result<dense::Mat, OpError> {
    let (shape, data) = array::<f64>(value, 2)?;
    if shape[0] != shape[1] {
        return Err(error("dense LU requires square matrix"));
    }
    Ok(dense::Mat {
        rows: shape[0],
        cols: shape[1],
        data,
    })
}
#[derive(Clone, Copy)]
enum Method {
    CsrFromTriplets,
    DenseBuild,
    DenseSolve,
    DenseInverse,
    DenseDet,
    CsrMatvec,
    CsrDiagonal,
    CsrJacobi,
    CsrBandwidth,
    CsrRcm,
    CsrMinDegree,
    CsrAmd,
    CsrNested,
    CsrPermute,
    IluBuild,
    IcBuild,
    AutoBuild,
    LuBuild,
    IluApply,
    IcApply,
    AutoApply,
    AutoBorrowApply,
    AutoKind,
    LuSolve,
    LuNnz,
    LuDim,
    LuPermutation,
}
struct Adapter {
    method: Method,
    signature: Signature,
    contract: DimContract,
}
fn make(method: Method) -> Adapter {
    use Method as M;
    let sparse = !matches!(
        method,
        M::DenseBuild | M::DenseSolve | M::DenseInverse | M::DenseDet
    );
    let build = matches!(
        method,
        M::DenseBuild | M::IluBuild | M::IcBuild | M::AutoBuild | M::LuBuild
    );
    let direct = matches!(
        method,
        M::CsrMatvec
            | M::CsrDiagonal
            | M::CsrJacobi
            | M::CsrBandwidth
            | M::CsrRcm
            | M::CsrMinDegree
            | M::CsrAmd
            | M::CsrNested
            | M::CsrPermute
    );
    let source = if sparse {
        super::super::iterative::csr_type()
    } else {
        tensor_type(2, TypeExpr::Float)
    };
    let mut args = vec![if build || direct {
        source
    } else {
        factor_type(sparse)
    }];
    let path = if build || direct {
        if sparse {
            &["values"][..]
        } else {
            &[][..]
        }
    } else if sparse {
        &["matrix", "values"][..]
    } else {
        &["matrix"][..]
    };
    let mut inputs = vec![input("D", 0, path)];
    let mut outputs = vec![];
    let mut result = None;
    let ret = match method {
        | M::DenseBuild => {
            outputs.push(output(&["matrix"], &D));
            factor_type(false)
        },
        | M::IluBuild | M::IcBuild | M::AutoBuild | M::LuBuild => {
            outputs.push(output(&["matrix", "values"], &D));
            factor_type(true)
        },
        | M::DenseSolve
        | M::LuSolve
        | M::IluApply
        | M::IcApply
        | M::AutoApply
        | M::AutoBorrowApply => {
            args.push(tensor_type(1, TypeExpr::Float));
            inputs.push(input("E", 1, &[]));
            result = Some(&DimFormula::Mul(&[E, DimFormula::Pow(&D, -1, 1)]));
            tensor_type(1, TypeExpr::Float)
        },
        | M::DenseInverse => {
            result = Some(&INVERSE);
            tensor_type(2, TypeExpr::Float)
        },
        | M::DenseDet => {
            outputs.push(output(&["entry_unit"], &D));
            TypeExpr::Record(vec![
                ("normalized".into(), TypeExpr::Float),
                ("entry_unit".into(), TypeExpr::Float),
                ("order".into(), TypeExpr::Int),
            ])
        },
        | M::CsrMatvec => {
            args.push(tensor_type(1, TypeExpr::Float));
            inputs.push(input("E", 1, &[]));
            result = Some(&PRODUCT);
            tensor_type(1, TypeExpr::Float)
        },
        | M::CsrDiagonal => {
            result = Some(&D);
            tensor_type(1, TypeExpr::Float)
        },
        | M::CsrJacobi => {
            result = Some(&INVERSE);
            tensor_type(1, TypeExpr::Float)
        },
        | M::CsrBandwidth | M::LuNnz | M::LuDim => TypeExpr::Int,
        | M::AutoKind => TypeExpr::String,
        | M::CsrPermute => {
            args.push(tensor_type(1, TypeExpr::Int));
            outputs.push(output(&["values"], &D));
            super::super::iterative::csr_type()
        },
        | _ => tensor_type(1, TypeExpr::Int),
    };
    if matches!(method, M::CsrFromTriplets) {
        args = vec![
            TypeExpr::Int,
            TypeExpr::Int,
            tensor_type(1, TypeExpr::Int),
            tensor_type(1, TypeExpr::Int),
            tensor_type(1, TypeExpr::Float),
        ];
        inputs = vec![input("D", 4, &[])];
        outputs = vec![output(&["values"], &D)];
        let contract = contract(inputs, vec![], outputs, None);
        let signature = dimensions::apply_contract(
            Signature::new(args, super::super::iterative::csr_type()),
            &contract,
        );
        return Adapter {
            method,
            signature,
            contract,
        };
    }
    let contract = contract(inputs, vec![], outputs, result);
    let signature = dimensions::apply_contract(Signature::new(args, ret), &contract);
    Adapter {
        method,
        signature,
        contract,
    }
}
impl EvalDyn for Adapter {
    fn signature(&self) -> Signature {
        self.signature.clone()
    }

    #[allow(clippy::too_many_lines)] // Native factor families share validation and dimensional preflight.
    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use Method as M;
        if args.len() != self.signature.params.len() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.signature.params.len() as u8),
                got: args.len(),
            });
        }
        let dims = runtime_dimensions::prepare_contract(&self.contract, args)?;
        let value = match self.method {
            | M::CsrFromTriplets => {
                let rows = integer(&args[0])?;
                let cols = integer(&args[1])?;
                if rows > SIDE || cols > SIDE {
                    return Err(error("CSR order exceeds portable budget"));
                }
                let r = decode_indices(&args[2])?;
                let c = decode_indices(&args[3])?;
                let Value::Tensor(t) = &args[4] else {
                    return Err(error("triplet values must be a tensor"));
                };
                if t.data.len() > LIMIT {
                    return Err(error("triplet count exceeds portable budget"));
                }
                let (_, v) = super::super::tensor(&args[4], 1, 4)?;
                if r.len() != v.len()
                    || c.len() != v.len()
                    || r.iter().any(|n| *n >= rows)
                    || c.iter().any(|n| *n >= cols)
                {
                    return Err(error(
                        "triplet arrays have unequal lengths or out-of-bounds indices",
                    ));
                }
                let triplets = r
                    .into_iter()
                    .zip(c)
                    .zip(v)
                    .map(|((r, c), v)| (r, c, v))
                    .collect::<Vec<_>>();
                encode_csr(krylov::Csr::from_triplets(rows, cols, &triplets))?
            },
            | M::DenseBuild => {
                let a = dense_matrix(&args[0])?;
                dense::lu_factor(&a).map_err(error)?;
                encode_factor(args[0].clone(), 0, false)
            },
            | M::DenseSolve | M::DenseInverse | M::DenseDet => {
                let (value, _) = factor(&args[0], false)?;
                let a = dense_matrix(value)?;
                let lu = dense::lu_factor(&a).map_err(error)?;
                match self.method {
                    | M::DenseSolve => {
                        let (_, b) = array::<f64>(&args[1], 1)?;
                        if b.len() != a.rows {
                            return Err(error("LU right-hand side length mismatch"));
                        }
                        encode_array(vec![b.len()], lu.solve(&b))?
                    },
                    | M::DenseInverse => {
                        let a = lu.inverse();
                        encode_array(vec![a.rows, a.cols], a.data)?
                    },
                    | _ => {
                        record(vec![
                            ("normalized", super::super::output_scalar(lu.det())?),
                            ("entry_unit", Value::Float(1.0)),
                            ("order", Value::int(a.rows)),
                        ])
                    },
                }
            },
            | M::IluBuild | M::IcBuild | M::AutoBuild | M::LuBuild => {
                let a = csr(&args[0])?;
                if a.rows != a.cols {
                    return Err(error("sparse factor requires a square matrix"));
                }
                let kind = match self.method {
                    | M::IluBuild => {
                        a.ilu0()
                            .ok_or_else(|| error("ILU(0) factorization failed"))?;
                        1
                    },
                    | M::IcBuild => {
                        a.ic0().ok_or_else(|| error("IC(0) factorization failed"))?;
                        2
                    },
                    | M::AutoBuild => {
                        let _ = krylov::AutoPreconditioner::new(&a);
                        3
                    },
                    | _ => {
                        a.sparse_lu(krylov::SparseOrdering::Natural)
                            .map_err(sparse_error)?;
                        4
                    },
                };
                encode_factor(encode_csr(a)?, kind, true)
            },
            | M::IluApply
            | M::IcApply
            | M::AutoApply
            | M::AutoBorrowApply
            | M::AutoKind
            | M::LuSolve
            | M::LuNnz
            | M::LuDim
            | M::LuPermutation => {
                let (value, kind) = factor(&args[0], true)?;
                let a = csr(value)?;
                if a.rows != a.cols {
                    return Err(error("sparse factor requires a square matrix"));
                }
                let expected = match self.method {
                    | M::IluApply => 1,
                    | M::IcApply => 2,
                    | M::AutoApply | M::AutoBorrowApply | M::AutoKind => 3,
                    | _ => 4,
                };
                if kind != expected {
                    return Err(error("factor method does not match portable factor kind"));
                }
                if matches!(self.method, M::AutoKind) {
                    Value::String(format!("{:?}", krylov::AutoPreconditioner::new(&a).kind()))
                } else if matches!(self.method, M::LuNnz | M::LuDim | M::LuPermutation) {
                    let lu = a
                        .sparse_lu(krylov::SparseOrdering::Natural)
                        .map_err(sparse_error)?;
                    match self.method {
                        | M::LuNnz => Value::int(lu.nnz()),
                        | M::LuDim => Value::int(lu.dim()),
                        | _ => indices(lu.permutation().to_vec()),
                    }
                } else {
                    let (_, b) = array::<f64>(&args[1], 1)?;
                    if b.len() != a.rows {
                        return Err(error("factor right-hand side length mismatch"));
                    }
                    let mut z = vec![0.0; b.len()];
                    match self.method {
                        | M::IluApply => {
                            a.ilu0()
                                .ok_or_else(|| error("ILU(0) factorization failed"))?
                                .apply(&b, &mut z);
                        },
                        | M::IcApply => {
                            a.ic0()
                                .ok_or_else(|| error("IC(0) factorization failed"))?
                                .apply(&b, &mut z);
                        },
                        | M::AutoApply => krylov::AutoPreconditioner::new(&a).apply(&b, &mut z),
                        | M::AutoBorrowApply => {
                            krylov::AutoPreconditioner::new(&a)
                                .as_preconditioner()
                                .apply(&b, &mut z);
                        },
                        | _ => {
                            z = a
                                .sparse_lu(krylov::SparseOrdering::Natural)
                                .map_err(sparse_error)?
                                .solve(&b);
                        },
                    }
                    encode_array(vec![z.len()], z)?
                }
            },
            | _ => {
                let a = csr(&args[0])?;
                if a.rows != a.cols
                    && matches!(
                        self.method,
                        M::CsrRcm | M::CsrMinDegree | M::CsrAmd | M::CsrNested | M::CsrPermute
                    )
                {
                    return Err(error("graph ordering requires a square CSR matrix"));
                }
                match self.method {
                    | M::CsrMatvec => {
                        let (_, x) = array::<f64>(&args[1], 1)?;
                        if x.len() != a.cols {
                            return Err(error("CSR vector length mismatch"));
                        }
                        let mut y = vec![0.0; a.rows];
                        a.matvec(&x, &mut y);
                        encode_array(vec![y.len()], y)?
                    },
                    | M::CsrDiagonal => {
                        let d = a.diagonal();
                        encode_array(vec![d.len()], d)?
                    },
                    | M::CsrJacobi => {
                        let d = a.jacobi_inverse();
                        encode_array(vec![d.len()], d)?
                    },
                    | M::CsrBandwidth => Value::int(a.bandwidth()),
                    | M::CsrRcm => indices(a.rcm_ordering()),
                    | M::CsrMinDegree => indices(a.min_degree_ordering()),
                    | M::CsrAmd => indices(a.amd_ordering()),
                    | M::CsrNested => indices(a.nested_dissection_ordering()),
                    | M::CsrPermute => {
                        let permutation = decode_indices(&args[1])?;
                        let mut sorted = permutation.clone();
                        sorted.sort_unstable();
                        if sorted != (0..a.rows).collect::<Vec<_>>() {
                            return Err(error(
                                "CSR permutation must contain each row index exactly once",
                            ));
                        }
                        encode_csr(a.permute_symmetric(&permutation))?
                    },
                    | _ => unreachable!(),
                }
            },
        };
        runtime_dimensions::finish(&self.signature, &dims, value)
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in METHODS {
        let adapter = make(*method);
        let signature = adapter.signature();
        let op=installer.op(OpDescriptor::new(name,Arity::Fixed(signature.params.len()as u8)).rust_path(path).category("numerical_matrix")
        .doc("Validated portable native dense/sparse factor definition. Reconstructs native factors from finite defining matrices, preserving method identity and right-hand-side/matrix quotient units. Format/version/kind, CSR sorting, shapes, indices, exact permutations and allocation budgets are checked before native methods. No process-local factor pointers are serialized.")
        .stability(Stability::Experimental).signature(signature).eval_dyn(Arc::new(adapter)))?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
const METHODS: &[(&str, &str, Method)] = &[
    (
        "num_csr_model_from_triplets",
        "crate::kernels::krylov::Csr::from_triplets",
        Method::CsrFromTriplets,
    ),
    (
        "num_dense_lu_model_build",
        "crate::kernels::dense::lu_factor",
        Method::DenseBuild,
    ),
    (
        "num_dense_lu_model_solve",
        "crate::kernels::dense::Lu::solve",
        Method::DenseSolve,
    ),
    (
        "num_dense_lu_model_inverse",
        "crate::kernels::dense::Lu::inverse",
        Method::DenseInverse,
    ),
    (
        "num_dense_lu_model_det",
        "crate::kernels::dense::Lu::det",
        Method::DenseDet,
    ),
    (
        "num_csr_model_matvec",
        "crate::kernels::krylov::Csr::matvec",
        Method::CsrMatvec,
    ),
    (
        "num_csr_model_diagonal",
        "crate::kernels::krylov::Csr::diagonal",
        Method::CsrDiagonal,
    ),
    (
        "num_csr_model_jacobi_inverse",
        "crate::kernels::krylov::Csr::jacobi_inverse",
        Method::CsrJacobi,
    ),
    (
        "num_csr_model_bandwidth",
        "crate::kernels::krylov::Csr::bandwidth",
        Method::CsrBandwidth,
    ),
    (
        "num_csr_model_rcm",
        "crate::kernels::krylov::Csr::rcm_ordering",
        Method::CsrRcm,
    ),
    (
        "num_csr_model_min_degree",
        "crate::kernels::krylov::Csr::min_degree_ordering",
        Method::CsrMinDegree,
    ),
    (
        "num_csr_model_amd",
        "crate::kernels::krylov::Csr::amd_ordering",
        Method::CsrAmd,
    ),
    (
        "num_csr_model_nested_dissection",
        "crate::kernels::krylov::Csr::nested_dissection_ordering",
        Method::CsrNested,
    ),
    (
        "num_csr_model_permute",
        "crate::kernels::krylov::Csr::permute_symmetric",
        Method::CsrPermute,
    ),
    (
        "num_sparse_ilu_model_build",
        "crate::kernels::krylov::Csr::ilu0",
        Method::IluBuild,
    ),
    (
        "num_sparse_ic_model_build",
        "crate::kernels::krylov::Csr::ic0",
        Method::IcBuild,
    ),
    (
        "num_sparse_auto_model_build",
        "crate::kernels::krylov::AutoPreconditioner::new",
        Method::AutoBuild,
    ),
    (
        "num_sparse_lu_model_build",
        "crate::kernels::krylov::Csr::sparse_lu",
        Method::LuBuild,
    ),
    (
        "num_sparse_ilu_model_apply",
        "crate::kernels::krylov::Ilu0::apply",
        Method::IluApply,
    ),
    (
        "num_sparse_ic_model_apply",
        "crate::kernels::krylov::Ic0::apply",
        Method::IcApply,
    ),
    (
        "num_sparse_auto_model_apply",
        "crate::kernels::krylov::AutoPreconditioner::apply",
        Method::AutoApply,
    ),
    (
        "num_sparse_auto_model_borrow_apply",
        "crate::kernels::krylov::AutoPreconditioner::as_preconditioner",
        Method::AutoBorrowApply,
    ),
    (
        "num_sparse_auto_model_kind",
        "crate::kernels::krylov::AutoPreconditioner::kind",
        Method::AutoKind,
    ),
    (
        "num_sparse_lu_model_solve",
        "crate::kernels::krylov::SparseLu::solve",
        Method::LuSolve,
    ),
    (
        "num_sparse_lu_model_nnz",
        "crate::kernels::krylov::SparseLu::nnz",
        Method::LuNnz,
    ),
    (
        "num_sparse_lu_model_dim",
        "crate::kernels::krylov::SparseLu::dim",
        Method::LuDim,
    ),
    (
        "num_sparse_lu_model_permutation",
        "crate::kernels::krylov::SparseLu::permutation",
        Method::LuPermutation,
    ),
];
