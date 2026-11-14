use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::rules::numerical_abi::EXPORTS;
use rssn::rules::numerical_abi::GENERATED_EXPORTS;
use std::collections::BTreeSet;

fn vector(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn matrix(
    rows: usize,
    cols: usize,
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows, cols],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}
#[test]
fn descriptors_register_exact_manifest_overloads_and_unique_keys() {
    let session = Session::new();
    let manifest = session.manifest();
    let mut names = BTreeSet::new();
    assert!(!GENERATED_EXPORTS.is_empty());
    for export in EXPORTS {
        assert!(names.insert(export.name));
        let info = manifest.op(export.name).expect(export.name);
        assert_eq!(info.rule_set, export.namespace);
        assert_eq!(info.rust_path.as_deref(), Some(export.rust_path));
        assert_eq!(info.doc, export.doc);
        assert_eq!(info.stability, export.stability);
        assert_eq!(info.overloads.len(), 1);
        assert_eq!(info.keys.len(), 1);
        assert_eq!(
            info.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
    }
}
#[test]
fn numeric_families_dispatch_through_name_and_key() {
    let session = Session::new();
    let manifest = session.manifest();
    let args = [vector(&[1.0, 2.0]), vector(&[3.0, 4.0])];
    let key = &manifest.op("num_dot").unwrap().keys[0];
    assert_eq!(
        call(&session, "num_dot", &args).unwrap(),
        Value::Float(11.0)
    );
    assert_eq!(
        session.call_key(key, &args, &CallCtx::new()).unwrap(),
        Value::Float(11.0)
    );
    assert_eq!(
        call(&session, "num_mean", &[vector(&[1.0, 2.0, 3.0])]).unwrap(),
        Value::Float(2.0)
    );
    assert_eq!(
        call(&session, "num_sigmoid", &[Value::Float(0.0)]).unwrap(),
        Value::Float(0.5)
    );
    assert_eq!(
        call(
            &session,
            "num_matmul",
            &[matrix(1, 2, &[2.0, 3.0]), matrix(2, 1, &[4.0, 5.0])]
        )
        .unwrap(),
        matrix(1, 1, &[23.0])
    );
}
#[test]
fn numeric_adapters_reject_bad_shapes_and_non_finite_values() {
    let session = Session::new();
    assert!(matches!(
        call(&session, "num_dot", &[vector(&[1.0]), vector(&[1.0, 2.0])]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(&session, "num_norm", &[vector(&[f64::NAN])]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_matmul",
            &[matrix(1, 2, &[1.0, 2.0]), matrix(1, 2, &[1.0, 2.0])]
        ),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_cholesky",
            &[matrix(2, 2, &[1.0, 2.0, 2.0, 1.0])]
        ),
        Err(OpError::Domain(_))
    ));
    let malformed = matrix(2, 2, &[1.0]);
    assert!(matches!(
        call(&session, "num_transpose", &[malformed]),
        Err(OpError::Domain(_))
    ));
    let transposed = call(&session, "num_transpose", &[matrix(2, 1, &[4.0, 5.0])]).unwrap();
    assert_eq!(
        call(
            &session,
            "num_matmul",
            &[transposed, matrix(2, 1, &[2.0, 3.0])]
        )
        .unwrap(),
        matrix(1, 1, &[23.0])
    );
}

#[test]
fn structured_codecs_preserve_schema_and_decompositions() {
    let session = Session::new();
    let manifest = session.manifest();
    for (name, fields) in [
        ("num_svd", vec!["u", "s", "v"]),
        ("num_qr", vec!["q", "r"]),
        ("num_eigen_symmetric", vec!["values", "vectors"]),
    ] {
        let info = manifest.op(name).expect(name);
        let TypeExpr::Record(schema) = &info.overloads[0].ret else {
            panic!("{name} has no record signature");
        };
        assert_eq!(
            schema
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            fields
        );
        let result = call(&session, name, &[matrix(2, 2, &[3.0, 0.0, 0.0, 2.0])]).unwrap();
        let Value::Record(record) = result else {
            panic!("{name} has no record result");
        };
        assert_eq!(
            record
                .fields
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            fields
        );
        for ((_, field_ty), (_, field)) in schema.iter().zip(&record.fields) {
            assert!(rssn::kernel::overload::matches_type(field_ty, &field.type_expr()).accepted());
        }
    }
    let Value::Record(svd) =
        call(&session, "num_svd", &[matrix(2, 2, &[3.0, 0.0, 0.0, 2.0])]).unwrap()
    else {
        panic!("svd record");
    };
    let (Value::Tensor(u), Value::Tensor(s), Value::Tensor(v)) = (
        svd.get("u").unwrap(),
        svd.get("s").unwrap(),
        svd.get("v").unwrap(),
    ) else {
        panic!("svd tensors");
    };
    for i in 0..2 {
        for j in 0..2 {
            let mut actual = 0.0;
            for k in 0..2 {
                actual += u.data[i * 2 + k].as_f64().unwrap()
                    * s.data[k].as_f64().unwrap()
                    * v.data[j * 2 + k].as_f64().unwrap();
            }
            let expected = [3.0, 0.0, 0.0, 2.0][i * 2 + j];
            assert!((actual - expected).abs() < 1e-12);
        }
    }
}
#[test]
fn integer_and_multiple_argument_codecs_check_domains() {
    let session = Session::new();
    let Value::Float(fact) = call(&session, "num_factorial", &[Value::int(5)]).unwrap() else {
        panic!("factorial");
    };
    assert!((fact - 120.0).abs() < 1e-9);
    let Value::Float(comb) =
        call(&session, "num_binomial", &[Value::int(6), Value::int(2)]).unwrap()
    else {
        panic!("binomial");
    };
    assert!((comb - 15.0).abs() < 1e-9);
    assert!(matches!(
        call(&session, "num_factorial", &[Value::int(-1)]),
        Err(OpError::Domain(_))
    ));
    assert!(call(&session, "num_factorial", &[Value::Float(2.5)]).is_err());
    assert!(matches!(
        call(
            &session,
            "num_rank_tol",
            &[matrix(1, 1, &[1.0]), Value::Float(-1.0)]
        ),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_lstsq",
            &[matrix(2, 1, &[1.0, 2.0]), vector(&[1.0])]
        ),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn documented_infinity_is_not_confused_with_invalid_nan() {
    let session = Session::new();
    let condition = call(&session, "num_cond", &[matrix(2, 2, &[1.0, 0.0, 0.0, 0.0])]).unwrap();
    assert!(matches!(condition, Value::Float(value) if value == f64::INFINITY));
    assert!(matches!(
        call(&session, "num_geometric_mean", &[vector(&[-1.0, 2.0])]),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn tensor_train_export_keeps_shapes_cores_and_error_report() {
    let session = Session::new();
    let shape = Value::Tensor(TensorValue {
        shape: vec![2],
        data: vec![Value::int(2), Value::int(2)],
        elem: TypeExpr::Int,
    });
    let args = [
        vector(&[2.0, 0.0, 0.0, 2.0]),
        shape,
        Value::Float(0.0),
        Value::int(2),
    ];
    let Value::Record(result) = call(&session, "num_tt_svd", &args).unwrap() else {
        panic!("tensor train record");
    };
    assert_eq!(result.get("shape"), Some(&args[1]));
    assert_eq!(result.get("tolerance_met"), Some(&Value::Bool(true)));
    assert!(matches!(result.get("error"),Some(Value::Float(x)) if *x < 1e-12));
    let Some(Value::Tensor(cores)) = result.get("cores") else {
        panic!("core list");
    };
    assert_eq!(cores.shape, vec![2]);
    assert_eq!(result.get("scale"), Some(&Value::Float(2.0)));
    let (Value::Tensor(left), Value::Tensor(right)) = (&cores.data[0], &cores.data[1]) else {
        panic!("core tensors");
    };
    let rank = left.shape[2];
    for i in 0..2 {
        for j in 0..2 {
            let value = (0..rank)
                .map(|r| {
                    left.data[i * rank + r].as_f64().unwrap()
                        * right.data[r * 2 + j].as_f64().unwrap()
                })
                .sum::<f64>()
                * 2.0;
            assert!((value - if i == j { 2.0 } else { 0.0 }).abs() < 1e-12);
        }
    }
    assert!(
        cores
            .data
            .iter()
            .all(|v| matches!(v,Value::Tensor(t) if t.shape.len() == 3))
    );
}

#[test]
fn dimensional_policies_preserve_and_combine_witnesses() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let session = Session::new();
    let tensor = |unit: &str, xs: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![xs.len()],
            data: xs.iter().copied().map(Value::Float).collect(),
            elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base(unit.into())),
        })
    };
    let a = tensor("Length", &[1.0, 2.0]);
    let b = tensor("Time", &[3.0, 4.0]);
    let Value::Dyn(DynValue { value, dim }) =
        call(&session, "num_dot", &[a.clone(), b.clone()]).unwrap()
    else {
        panic!("dot witness");
    };
    assert_eq!(*value, Value::Float(11.0));
    assert_eq!(
        dim,
        DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Base("Time".into())
        ])
    );
    let Value::Dyn(DynValue { value, dim }) =
        call(&session, "num_mean", std::slice::from_ref(&a)).unwrap()
    else {
        panic!("mean witness");
    };
    assert_eq!(*value, Value::Float(1.5));
    assert_eq!(dim, DimExpr::Base("Length".into()));
    assert!(matches!(
        call(&session, "num_vector_add", &[a, b]),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn quantities_are_converted_to_coherent_units_before_numerical_work() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::QuantityValue;
    let session = Session::new();
    let mut unit = session
        .manifest()
        .units
        .iter()
        .find(|u| u.name == "m")
        .expect("metre")
        .clone();
    unit.scale.numer = 1000.into();
    let argument = Value::Tensor(TensorValue {
        shape: vec![1],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base("Length".into())),
        data: vec![Value::Quantity(QuantityValue {
            value: Box::new(Value::Float(2.0)),
            unit,
        })],
    });
    let Value::Dyn(DynValue { value, dim }) = call(&session, "num_mean", &[argument]).unwrap()
    else {
        panic!("quantity result witness");
    };
    assert_eq!(*value, Value::Float(2000.0));
    assert_eq!(dim, DimExpr::Base("Length".into()));
}

#[test]
fn generated_source_bindings_dispatch_and_obey_output_policies() {
    let session = Session::new();
    let product = call(
        &session,
        "num_cross_product",
        &[vector(&[1.0, 0.0, 0.0]), vector(&[0.0, 1.0, 0.0])],
    )
    .expect("3D cross product");
    let Value::Tensor(product) = product else {
        panic!("expected tensor")
    };
    assert_eq!(
        product.data,
        vec![Value::Float(0.0), Value::Float(0.0), Value::Float(1.0)]
    );
    assert!(
        call(
            &session,
            "num_cross_product",
            &[vector(&[1.0, 0.0]), vector(&[0.0, 1.0])]
        )
        .is_err()
    );
    assert_eq!(
        call(&session, "num_inverse_erfc", &[Value::Float(0.0)]).unwrap(),
        Value::Float(f64::INFINITY)
    );
    assert!(call(&session, "num_inverse_erfc", &[Value::Float(-1.0)]).is_err());
    let Value::Float(result) = call(&session, "num_erfi", &[Value::Float(1.0)]).unwrap() else {
        panic!("expected scalar")
    };
    assert!((result - 1.6504257588).abs() < 1e-8);
    let sample = vector(&[9.0, 1.0, 3.0]);
    assert_eq!(
        call(&session, "num_median", std::slice::from_ref(&sample)).unwrap(),
        Value::Float(3.0)
    );
    assert_eq!(
        sample,
        vector(&[9.0, 1.0, 3.0]),
        "mutable Rust kernels must not mutate caller tensors"
    );
    assert!(call(&session, "num_skewness", &[vector(&[1.0, 2.0])]).is_err());
    assert_eq!(
        call(
            &session,
            "num_lp_norm",
            &[vector(&[3.0, 4.0]), Value::Float(2.0)]
        )
        .unwrap(),
        Value::Float(5.0)
    );
    assert!(
        call(
            &session,
            "num_lp_norm",
            &[vector(&[3.0]), Value::Float(0.0)]
        )
        .is_err()
    );
    let Value::Tensor(scaled) = call(
        &session,
        "num_scalar_mul",
        &[vector(&[2.0, 3.0]), Value::Float(4.0)],
    )
    .unwrap() else {
        panic!("expected tensor")
    };
    assert_eq!(scaled.data, vec![Value::Float(8.0), Value::Float(12.0)]);
}

#[test]
fn generated_structured_policies_are_exhaustive_and_honest_about_unbound_types() {
    use rssn::rules::numerical_abi::codecs::STRUCTURED_POLICIES;
    use rssn::rules::numerical_abi::dimensions::FieldDimension;
    let svd = STRUCTURED_POLICIES
        .iter()
        .find(|policy| policy.rust_path == "crate::kernels::dense::Svd")
        .unwrap();
    assert!(svd.bound);
    assert_eq!(
        svd.fields.iter().map(|(name, _)| *name).collect::<Vec<_>>(),
        ["u", "s", "v"]
    );
    assert!(
        svd.fields
            .iter()
            .all(|(_, dimension)| matches!(dimension, FieldDimension::Formula(_)))
    );
    let integral = STRUCTURED_POLICIES
        .iter()
        .find(|policy| policy.rust_path == "crate::kernels::quadrature::Integral")
        .unwrap();
    assert!(integral.bound);
    assert!(matches!(integral.fields[0].1, FieldDimension::Formula(_)));
    assert!(matches!(integral.fields[1].1, FieldDimension::Formula(_)));
    assert!(matches!(integral.fields[2].1, FieldDimension::NonNumeric));
    let iterative = STRUCTURED_POLICIES
        .iter()
        .find(|policy| {
            policy.rust_path == "crate::kernels::krylov::IterResult"
                && policy.contract == "ITERATIVE_LINEAR_SOLVE"
        })
        .unwrap();
    assert!(iterative.bound);
    assert!(matches!(iterative.fields[0].1, FieldDimension::Formula(_)));
    assert!(matches!(iterative.fields[1].1, FieldDimension::Formula(_)));
    assert!(matches!(iterative.fields[2].1, FieldDimension::NonNumeric));
}

fn csr(
    indptr: &[i64],
    indices: &[i64],
    values: &[f64],
) -> Value {
    use rssn::kernel::value::RecordValue;
    let ints = |items: &[i64]| {
        Value::Tensor(TensorValue {
            shape: vec![items.len()],
            data: items.iter().copied().map(Value::int).collect(),
            elem: TypeExpr::Int,
        })
    };
    Value::Record(RecordValue {
        fields: vec![
            ("rows".into(), Value::int(2)),
            ("cols".into(), Value::int(2)),
            ("indptr".into(), ints(indptr)),
            ("indices".into(), ints(indices)),
            ("values".into(), vector(values)),
        ],
    })
}
fn result_field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected record")
    };
    &record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .unwrap()
        .1
}
#[test]
fn structural_csr_solvers_check_residuals_and_manifest_origins() {
    let session = Session::new();
    for name in ["num_csr_cg", "num_csr_gmres", "num_csr_bicgstab"] {
        let mut args = vec![
            csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
            vector(&[2.0, 8.0]),
            Value::Float(1e-10),
            Value::int(40),
        ];
        if name == "num_csr_gmres" {
            args.push(Value::int(2));
        }
        let result = call(&session, name, &args).unwrap();
        assert_eq!(result_field(&result, "converged"), &Value::Bool(true));
        assert!(result_field(&result, "residual").as_f64().unwrap() < 1e-10);
        let key = session.manifest().op(name).unwrap().keys[0].clone();
        assert_eq!(
            session.call_key(&key, &args, &CallCtx::new()).unwrap(),
            result
        );
    }
}
#[test]
fn structural_csr_rejects_bad_rows_columns_and_controls() {
    let session = Session::new();
    let good_rhs = vector(&[2.0, 8.0]);
    for matrix in [
        csr(&[0, 2, 1], &[0, 1], &[2.0, 4.0]),
        csr(&[0, 1, 2], &[2, 1], &[2.0, 4.0]),
        csr(&[0, 2, 2], &[0, 0], &[2.0, 4.0]),
        csr(&[0, 1, 3], &[0, 1], &[2.0, 4.0]),
    ] {
        assert!(
            call(
                &session,
                "num_csr_cg",
                &[matrix, good_rhs.clone(), Value::Float(1e-8), Value::int(10)]
            )
            .is_err()
        );
    }
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
                vector(&[1.0]),
                Value::Float(1e-8),
                Value::int(10)
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_csr_gmres",
            &[
                csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
                good_rhs.clone(),
                Value::Float(1e-8),
                Value::int(10),
                Value::int(0)
            ]
        )
        .is_err()
    );
    let mut oversized = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut record) = oversized else {
        unreachable!()
    };
    record.fields[0].1 = Value::int(-1);
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                oversized,
                good_rhs.clone(),
                Value::Float(1e-8),
                Value::int(10)
            ]
        )
        .is_err()
    );
    let mut duplicate = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut record) = duplicate else {
        unreachable!()
    };
    record.fields[1].0 = "rows".into();
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[duplicate, good_rhs, Value::Float(1e-8), Value::int(10)]
        )
        .is_err()
    );
}

#[test]
fn structural_csr_zero_rhs_and_iteration_limit_are_honest() {
    let session = Session::new();
    let matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let zero = call(
        &session,
        "num_csr_cg",
        &[
            matrix.clone(),
            vector(&[0.0, 0.0]),
            Value::Float(1e-10),
            Value::int(5),
        ],
    )
    .unwrap();
    assert_eq!(result_field(&zero, "converged"), &Value::Bool(true));
    assert_eq!(result_field(&zero, "residual").as_f64(), Some(0.0));
    let limited = call(
        &session,
        "num_csr_cg",
        &[
            matrix,
            vector(&[2.0, 8.0]),
            Value::Float(1e-14),
            Value::int(1),
        ],
    )
    .unwrap();
    assert_eq!(result_field(&limited, "converged"), &Value::Bool(false));
    assert!(result_field(&limited, "residual").as_f64().unwrap() > 1e-14);
}
#[test]
fn structural_csr_preserves_units_and_rejects_conflicting_element_witnesses() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let mut matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut fields) = matrix else {
        unreachable!()
    };
    let Value::Tensor(values) = &mut fields.fields[4].1 else {
        unreachable!()
    };
    values.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), length.clone());
    let mut rhs = vector(&[2.0, 8.0]);
    let Value::Tensor(rhs_values) = &mut rhs else {
        unreachable!()
    };
    rhs_values.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), time.clone());
    let session = Session::new();
    let result = call(
        &session,
        "num_csr_cg",
        &[
            matrix.clone(),
            rhs.clone(),
            Value::Float(1e-10),
            Value::int(10),
        ],
    )
    .unwrap();
    let Value::Tensor(x) = result_field(&result, "x") else {
        panic!("solution tensor")
    };
    let TypeExpr::Dim(_, dimension) = &x.elem else {
        panic!("solution unit")
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(dimension.clone()),
            &TypeExpr::Float.with_dim(DimExpr::Mul(vec![
                time.clone(),
                DimExpr::Inverse(Box::new(length.clone()))
            ]))
        )
        .accepted()
    );
    let Value::Record(ref mut fields) = matrix else {
        unreachable!()
    };
    let Value::Tensor(values) = &mut fields.fields[4].1 else {
        unreachable!()
    };
    values.data[0] = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: time,
    });
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                matrix.clone(),
                rhs.clone(),
                Value::Float(1e-10),
                Value::int(10)
            ]
        )
        .is_err()
    );
    let matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let bad_tol = Value::Dyn(DynValue {
        value: Box::new(Value::Float(1e-10)),
        dim: length,
    });
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[matrix, rhs, bad_tol, Value::int(10)]
        )
        .is_err()
    );
}

#[test]
fn structural_csr_bindings_solve_uniformly_tiny_and_huge_systems() {
    let session = Session::new();
    for scale in [f64::from_bits(1), 1e-320, 1e-300, 1e300] {
        for name in ["num_csr_cg", "num_csr_gmres", "num_csr_bicgstab"] {
            let mut args = vec![
                csr(&[0, 1, 2], &[0, 1], &[scale, scale]),
                vector(&[scale, -scale]),
                Value::Float(1e-10),
                Value::int(12),
            ];
            if name == "num_csr_gmres" {
                args.push(Value::int(2));
            }
            let result = call(&session, name, &args).unwrap();
            assert_eq!(
                result_field(&result, "converged"),
                &Value::Bool(true),
                "{name}, scale {scale:e}"
            );
            let Value::Tensor(x) = result_field(&result, "x") else {
                panic!("solution tensor")
            };
            assert!(
                (x.data[0].as_f64().unwrap() - 1.0).abs() < 1e-10,
                "{name}, scale {scale:e}"
            );
            assert!(
                (x.data[1].as_f64().unwrap() + 1.0).abs() < 1e-10,
                "{name}, scale {scale:e}"
            );
        }
    }
}

#[test]
fn generated_contract_table_preserves_dimension_rules_across_families() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let session = Session::new();
    let dimensioned = |name: &str, xs: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![xs.len()],
            data: xs.iter().copied().map(Value::Float).collect(),
            elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base(name.into())),
        })
    };
    let length = dimensioned("Length", &[1.0, 3.0]);
    let time = dimensioned("Time", &[1.0, 0.0]);
    let Value::Dyn(variance) =
        call(&session, "num_variance", std::slice::from_ref(&length)).unwrap()
    else {
        panic!("variance quantity")
    };
    let expected = DimExpr::Mul(vec![
        DimExpr::Base("Length".into()),
        DimExpr::Base("Length".into()),
    ]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(variance.dim),
            &TypeExpr::Float.with_dim(expected)
        )
        .accepted()
    );
    assert!(matches!(
        call(
            &session,
            "num_cosine_similarity",
            &[length.clone(), time.clone()]
        )
        .unwrap(),
        Value::Float(_)
    ));
    let Value::Tensor(projected) = call(&session, "num_project", &[length.clone(), time]).unwrap()
    else {
        panic!("projection")
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &projected.elem,
            &TypeExpr::Float.with_dim(DimExpr::Base("Length".into()))
        )
        .accepted()
    );
    let empty = dimensioned("Length", &[]);
    let Value::Dyn(mean) = call(&session, "num_mean", &[empty]).unwrap() else {
        panic!("empty mean witness")
    };
    assert_eq!(mean.dim, DimExpr::Base("Length".into()));
    let bad_scale = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: DimExpr::Base("Time".into()),
    });
    assert!(call(&session, "num_lp_norm", &[length, bad_scale]).is_err());
}

#[test]
fn source_codecs_reject_dimensions_hidden_inside_optional_arguments() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::VariantValue;
    let session = Session::new();
    let optional = |value| {
        Value::Variant(VariantValue {
            tag: "Some".into(),
            payload: vec![value],
        })
    };
    let args = [
        vector(&[0.0, 1.0, 2.0]),
        vector(&[1.0]),
        optional(Value::Float(0.5)),
    ];
    let valid = call(&session, "num_stats_ext_kde", &args);
    assert!(
        valid.is_ok(),
        "valid KDE optional bandwidth rejected: {valid:?}"
    );
    let none = Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    let inferred = call(
        &session,
        "num_stats_ext_kde",
        &[args[0].clone(), args[1].clone(), none],
    );
    assert!(
        inferred.is_ok(),
        "valid KDE None bandwidth rejected: {inferred:?}"
    );
    for malformed in [
        Value::Variant(VariantValue {
            tag: "None".into(),
            payload: vec![Value::Float(0.5)],
        }),
        Value::Variant(VariantValue {
            tag: "Some".into(),
            payload: vec![],
        }),
        Value::Variant(VariantValue {
            tag: "Unknown".into(),
            payload: vec![],
        }),
        optional(Value::Bool(true)),
    ] {
        assert!(
            call(
                &session,
                "num_stats_ext_kde",
                &[args[0].clone(), args[1].clone(), malformed]
            )
            .is_err()
        );
    }
    let test_result = call(
        &session,
        "num_stats_ext_t_test_one_sample",
        &[vector(&[1.0, 2.0, 3.0]), Value::Float(2.0)],
    )
    .unwrap();
    let Value::Record(test_result) = test_result else {
        panic!("expected test-result record")
    };
    assert_eq!(test_result.get("df"), Some(&optional(Value::Float(2.0))));
    let ks_result = call(
        &session,
        "num_stats_ext_ks_test_two_sample",
        &[vector(&[1.0, 2.0]), vector(&[1.0, 2.0])],
    )
    .unwrap();
    let Value::Record(ks_result) = ks_result else {
        panic!("expected test-result record")
    };
    assert_eq!(
        ks_result.get("df"),
        Some(&Value::Variant(VariantValue {
            tag: "None".into(),
            payload: vec![]
        }))
    );
    let dimensioned = Value::Dyn(DynValue {
        value: Box::new(Value::Float(0.5)),
        dim: DimExpr::Base("Length".into()),
    });
    let bad_args = [args[0].clone(), args[1].clone(), optional(dimensioned)];
    let rejected = call(&session, "num_stats_ext_kde", &bad_args);
    assert!(
        rejected.is_err(),
        "dimensioned optional bandwidth was accepted: {rejected:?}"
    );
}
