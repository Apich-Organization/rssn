//! Native tensor and block-sparse adapters preserve explicit reference scales.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Result<Value, OpError> {
    s.call_value(n, a, &CallCtx::new())
}
fn ints(a: &[usize]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![a.len()],
        elem: TypeExpr::Int,
        data: a.iter().map(|v| Value::int(*v)).collect(),
    })
}
fn floats(a: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![a.len()],
        elem: TypeExpr::Float,
        data: a.iter().copied().map(Value::Float).collect(),
    })
}
fn measured(
    x: f64,
    d: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: DimExpr::Base(d.into()),
    })
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    &r.fields.iter().find(|(k, _)| k == n).unwrap().1
}
fn scalar(v: &Value) -> f64 {
    if let Value::Dyn(d) = v {
        scalar(&d.value)
    } else {
        v.as_f64().unwrap()
    }
}
fn array(
    s: &Session,
    shape: &[usize],
    data: &[f64],
    reference: Value,
) -> Value {
    call(
        s,
        "num_linear_array_new",
        &[ints(shape), floats(data), reference],
    )
    .unwrap()
}

#[test]
fn arnoldi_handles_extreme_normalized_scales_and_keeps_the_operator_reference() {
    let s = Session::new();
    let reference = measured(2.0, "Time");
    for gain in [2.0_f64.powi(-600), 2.0_f64.powi(600)] {
        let operator = array(
            &s,
            &[2, 2],
            &[gain, 0.0, 0.0, 2.0 * gain],
            reference.clone(),
        );
        let result = call(
            &s,
            "num_linear_arnoldi",
            &[operator, floats(&[1e200, 1e200]), Value::int(2)],
        )
        .unwrap();
        assert_eq!(field(&result, "operator_reference"), &reference);
        let Value::Tensor(h) = field(&result, "hessenberg") else {
            panic!("expected a Hessenberg matrix")
        };
        assert_eq!(h.shape, [2, 2]);
        assert_eq!(h.elem, TypeExpr::Float);
        assert!((scalar(&h.data[0]) / gain - 1.5).abs() <= 1e-14);
        assert!((scalar(&h.data[2]) / gain - 0.5).abs() <= 1e-14);
    }
}

#[test]
fn randomized_svd_preserves_extreme_singular_values_and_reference_units() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let reference = measured(2.0, "Length");
    for gain in [2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
        let m = array(
            &session,
            &[4, 3],
            &[
                3.0 * gain,
                0.0,
                0.0,
                0.0,
                2.0 * gain,
                0.0,
                0.0,
                0.0,
                gain,
                0.0,
                0.0,
                0.0,
            ],
            reference.clone(),
        );
        let actual = call(
            &session,
            "num_linear_randomized_svd",
            &[
                m,
                Value::int(3),
                Value::int(1),
                Value::int(2),
                Value::int(17),
            ],
        )
        .unwrap();
        let tensor = |field_name| {
            match field(&actual, field_name) {
                | Value::Tensor(tensor) => tensor,
                | other => panic!("expected tensor, got {other:?}"),
            }
        };
        assert_eq!(tensor("u").shape, [4, 3]);
        assert_eq!(tensor("v").shape, [3, 3]);
        assert_eq!(tensor("s").shape, [3]);
        let used_reference = field(&actual, "singular_value_reference");
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(length.clone()),
                &used_reference.type_expr(),
            )
            .accepted()
        );
        assert_eq!(scalar(used_reference), 2.0);
        for (value, expected) in tensor("s").data.iter().zip([6.0, 4.0, 2.0]) {
            let scaled = scalar(value) / gain * scalar(used_reference);
            assert!((scaled - expected).abs() < 1e-12, "{scaled} != {expected}");
        }
    }
}

#[test]
fn tensor_contraction_preserves_uncontracted_axis_order_and_product_units() {
    let s = Session::new();
    let a = array(&s, &[2, 2], &[1.0, 2.0, 2.0, 4.0], measured(3.0, "Time"));
    let b = array(
        &s,
        &[2, 3],
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        measured(3.0, "Length"),
    );
    let op = "num_linear_tensordot";
    let out = call(&s, op, &[a.clone(), b.clone(), ints(&[1]), ints(&[0])]).unwrap();
    let packed = call(&s, "num_linear_array_data", &[out]).unwrap();
    assert_eq!(field(&packed, "shape"), &ints(&[2, 3]));
    assert_eq!(
        field(&packed, "data"),
        &floats(&[9.0, 12.0, 15.0, 18.0, 24.0, 30.0])
    );
    let reference = field(&packed, "reference");
    let product = DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time")]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(product),
            &reference.type_expr()
        )
        .accepted()
    );
    assert_eq!(scalar(reference), 9.0);
    assert!(call(&s, op, &[a, b, ints(&[1, 1]), ints(&[0, 0])]).is_err());
}

#[test]
fn schur_plans_keep_quotient_units_and_accept_independent_right_hand_sides() {
    let s = Session::new();
    let a = array(&s, &[1, 1], &[2.0], measured(2.0, "Time"));
    let b = array(&s, &[1, 1], &[1.0], measured(4.0, "Time"));
    let new = "num_linear_sylvester_plan_new";
    let lyapunov = "num_linear_lyapunov_plan_new";
    let solve = "num_linear_sylvester_plan_solve";
    let syl = call(&s, new, &[a.clone(), b]).unwrap();
    let lyap = call(&s, lyapunov, std::slice::from_ref(&a)).unwrap();
    let quotient = DimExpr::Mul(vec![
        DimExpr::Base("Length".into()),
        DimExpr::Inverse(Box::new(DimExpr::Base("Time".into()))),
    ]);
    let expected = TypeExpr::Float.with_dim(quotient);
    for plan in [&syl, &lyap] {
        for sample in [2.0, -4.0] {
            let rhs = array(&s, &[1, 1], &[sample], measured(6.0, "Length"));
            let result = call(&s, solve, &[plan.clone(), rhs]).unwrap();
            let packed = call(&s, "num_linear_array_data", &[result]).unwrap();
            assert_eq!(field(&packed, "shape"), &ints(&[1, 1]));
            let reference = field(&packed, "reference");
            assert!(
                rssn::kernel::overload::matches_type(&expected, &reference.type_expr()).accepted()
            );
            let Value::Tensor(data) = field(&packed, "data") else {
                panic!("expected normalized solution tensor");
            };
            assert_eq!(
                data.data[0].as_f64().unwrap() * scalar(reference),
                sample * 0.75
            );
        }
    }
    let wrong_rhs = array(&s, &[1, 2], &[1.0, 2.0], measured(1.0, "Length"));
    assert!(call(&s, solve, &[syl, wrong_rhs]).is_err());
    let wrong_coefficient = array(&s, &[1, 1], &[1.0], measured(1.0, "Length"));
    assert!(call(&s, new, &[a, wrong_coefficient]).is_err());
}

#[test]
fn amg_precondition_many_reuses_one_hierarchy_and_keeps_quotient_units() {
    let s = Session::new();
    let matrix = call(
        &s,
        "num_linear_csr_new",
        &[
            Value::int(2),
            Value::int(2),
            ints(&[0, 1]),
            ints(&[0, 1]),
            floats(&[1.0, 2.0]),
            measured(2.0, "Time"),
        ],
    )
    .unwrap();
    let hierarchy = call(
        &s,
        "num_linear_amg_general_build",
        &[
            matrix,
            Value::String("RugeStuben".into()),
            Value::Float(0.25),
            Value::int(2),
            Value::int(4),
            Value::int(1),
            Value::String("V".into()),
            Value::Tensor(TensorValue {
                shape: vec![0, 2],
                elem: TypeExpr::Float,
                data: vec![],
            }),
        ],
    )
    .unwrap();
    let size = call(
        &s,
        "num_linear_amg_general_workspace_len",
        std::slice::from_ref(&hierarchy),
    )
    .unwrap();
    assert_eq!(size, Value::int(0));
    let name = "num_linear_amg_general_precondition_many";
    let rhs = array(&s, &[2, 2], &[2.0, 4.0, 8.0, -2.0], measured(6.0, "Length"));
    let x = call(&s, name, &[hierarchy.clone(), rhs]).unwrap();
    let result = call(&s, "num_linear_array_data", &[x]).unwrap();
    assert_eq!(field(&result, "shape"), &ints(&[2, 2]));
    assert_eq!(field(&result, "data"), &floats(&[2.0, 2.0, 8.0, -1.0]));
    let Value::Dyn(reference) = field(&result, "reference") else {
        panic!("measured quotient");
    };
    let expected_dim = DimExpr::Mul(vec![
        DimExpr::Base("Length".into()),
        DimExpr::Inverse(Box::new(DimExpr::Base("Time".into()))),
    ]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(expected_dim),
            &TypeExpr::Float.with_dim(reference.dim.clone()),
        )
        .accepted()
    );
    assert_eq!(scalar(&reference.value), 3.0);

    let empty = array(&s, &[0, 2], &[], measured(6.0, "Length"));
    let out = call(&s, name, &[hierarchy.clone(), empty]).unwrap();
    let data = call(&s, "num_linear_array_data", &[out]).unwrap();
    assert_eq!(field(&data, "shape"), &ints(&[0, 2]));
    assert_eq!(field(&data, "data"), &floats(&[]));
    let wrong = array(&s, &[2, 1], &[1.0, 2.0], measured(6.0, "Length"));
    assert!(call(&s, name, &[hierarchy, wrong]).is_err());
}

#[test]
fn small_immutable_amg_hierarchies_share_the_session_budget() {
    let session = Session::new();
    let matrix = call(
        &session,
        "num_linear_csr_new",
        &[
            Value::int(2),
            Value::int(2),
            ints(&[0, 1]),
            ints(&[0, 1]),
            floats(&[1.0, 2.0]),
            measured(2.0, "Time"),
        ],
    )
    .unwrap();
    let controls = [
        Value::Float(0.25),
        Value::int(2),
        Value::int(4),
        Value::int(1),
        Value::String("V".into()),
        Value::Tensor(TensorValue {
            shape: vec![0, 2],
            elem: TypeExpr::Float,
            data: vec![],
        }),
    ];
    let mut general = vec![matrix.clone(), Value::String("RugeStuben".into())];
    general.extend_from_slice(&controls);
    let mut symmetric = vec![matrix];
    symmetric.extend_from_slice(&controls);
    symmetric.push(Value::Bool(true));
    let mut hierarchies = Vec::new();
    for _ in 0..8 {
        hierarchies.push(call(&session, "num_linear_amg_general_build", &general).unwrap());
        hierarchies.push(call(&session, "num_linear_amg_build", &symmetric).unwrap());
    }
    let rhs = array(&session, &[2], &[2.0, 4.0], measured(6.0, "Length"));
    for (index, hierarchy) in hierarchies.into_iter().enumerate() {
        let operation = if index % 2 == 0 {
            "num_linear_amg_general_precondition"
        } else {
            "num_linear_amg_precondition"
        };
        let solved = call(&session, operation, &[hierarchy, rhs.clone()]).unwrap();
        let packed = call(&session, "num_linear_array_data", &[solved]).unwrap();
        assert_eq!(field(&packed, "data"), &floats(&[2.0, 2.0]));
    }
}

#[test]
fn small_immutable_sparse_objects_charge_their_owned_buffers() {
    let session = Session::new();
    let mut matrices = Vec::new();
    for _ in 0..12 {
        let csr = call(
            &session,
            "num_linear_csr_new",
            &[
                Value::int(2),
                Value::int(2),
                ints(&[0, 1]),
                ints(&[0, 1]),
                floats(&[1.0, 2.0]),
                measured(2.0, "Time"),
            ],
        )
        .unwrap();
        let blocks = call(
            &session,
            "num_linear_bsr_from_csr",
            &[csr.clone(), Value::int(1), Value::int(1)],
        )
        .unwrap();
        matrices.push((csr, blocks));
    }
    for (csr, blocks) in matrices {
        assert_eq!(
            call(&session, "num_linear_bsr_stored_len", &[blocks]).unwrap(),
            Value::int(2)
        );
        assert_eq!(
            call(&session, "num_linear_csr_bandwidth", &[csr]).unwrap(),
            Value::int(0)
        );
    }
}

#[test]
fn interpolative_decomposition_has_a_typed_dimensionless_result() {
    let session = Session::new();
    let matrix = array(
        &session,
        &[2, 3],
        &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0],
        measured(2.0, "Length"),
    );
    let name = "num_linear_interpolative_decomposition";
    let manifest = session.manifest();
    let signature = &manifest.op(name).unwrap().overloads[0];
    assert!(matches!(&signature.ret, TypeExpr::Record(fields) if fields.len() == 2));
    let result = call(&session, name, &[matrix.clone(), Value::int(2)]).unwrap();
    assert_eq!(field(&result, "columns"), &ints(&[2]));
    let Value::Tensor(coefficients) = field(&result, "coefficients") else {
        panic!("expected interpolation coefficient tensor");
    };
    assert_eq!(coefficients.shape, [1, 3]);
    assert_eq!(coefficients.elem, TypeExpr::Float);
    for (actual, want) in coefficients
        .data
        .iter()
        .map(scalar)
        .zip([1.0 / 3.0, 2.0 / 3.0, 1.0])
    {
        assert!((actual - want).abs() < 1e-14);
    }
    assert!(call(&session, name, &[matrix, Value::int(3)]).is_err());
}

#[test]
fn train_inner_restores_a_representable_result_with_extreme_reference_scales() {
    let s = Session::new();
    let train = |row, dimension| {
        let a = array(&s, &[2], row, measured(1e200, dimension));
        call(
            &s,
            "num_linear_tt_svd",
            &[a, Value::Float(1e-12), Value::int(1)],
        )
        .unwrap()
    };
    let a = train(&[1.0, 0.0], "Length");
    let b = train(&[0.0, 1.0], "Time");
    let args = [a, b];
    assert_eq!(
        call(&s, "num_linear_train_inner_workspace_len", &args).unwrap(),
        Value::int(0)
    );
    let result = call(&s, "num_linear_train_inner", &args).unwrap();
    let Value::Dyn(result) = result else {
        panic!("measured inner product");
    };
    assert_eq!(scalar(&result.value), 0.0);
    assert_eq!(
        result.dim,
        DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Base("Time".into())
        ])
    );
}

#[test]
fn train_native_entry_roundtrip_and_inner_product_retain_units() {
    let s = Session::new();
    let a = array(&s, &[2, 2], &[1.0, 2.0, 2.0, 4.0], measured(3.0, "Length"));
    let train = call(
        &s,
        "num_linear_tt_svd",
        &[a.clone(), Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    let entry = call(
        &s,
        "num_linear_train_entry",
        &[train.clone(), ints(&[1, 1])],
    )
    .unwrap();
    assert!((scalar(&entry) - 12.0).abs() < 1e-10);
    assert!(matches!(entry, Value::Dyn(_)));
    let restored = call(
        &s,
        "num_linear_train_to_dense",
        std::slice::from_ref(&train),
    )
    .unwrap();
    let raw = call(&s, "num_linear_array_data", &[restored]).unwrap();
    let Value::Tensor(t) = field(&raw, "data") else {
        panic!("tensor")
    };
    assert!((scalar(&t.data[3]) - 4.0).abs() < 1e-10);
    assert_eq!(scalar(field(&raw, "reference")), 3.0);
    let inner = call(&s, "num_linear_train_inner", &[train.clone(), train]).unwrap();
    assert!((scalar(&inner) - 225.0).abs() < 1e-8);
    let other = array(&s, &[2, 2], &[1.0; 4], measured(3.0, "Time"));
    let b = call(
        &s,
        "num_linear_tt_svd",
        &[other, Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    let first = call(
        &s,
        "num_linear_tt_svd",
        &[a, Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    assert!(call(&s, "num_linear_train_add", &[first, b]).is_err());
}
#[test]
fn bsr_native_transpose_product_and_input_guards_are_real() {
    let s = Session::new();
    let a = call(
        &s,
        "num_linear_csr_new",
        &[
            Value::int(2),
            Value::int(2),
            ints(&[0, 0, 1]),
            ints(&[0, 1, 1]),
            floats(&[2.0, 3.0, 4.0]),
            measured(2.0, "Time"),
        ],
    )
    .unwrap();
    let b = call(
        &s,
        "num_linear_bsr_from_csr",
        &[a, Value::int(1), Value::int(1)],
    )
    .unwrap();
    let x = array(&s, &[2], &[1.0, 2.0], measured(3.0, "Length"));
    let y = call(&s, "num_linear_bsr_matvec_transpose", &[b.clone(), x]).unwrap();
    let y = call(&s, "num_linear_array_data", &[y]).unwrap();
    let Value::Tensor(t) = field(&y, "data") else {
        panic!("tensor")
    };
    assert_eq!(scalar(&t.data[0]), 2.0);
    assert_eq!(scalar(&t.data[1]), 11.0);
    assert_eq!(scalar(field(&y, "reference")), 6.0);
    assert!(
        call(
            &s,
            "num_linear_bsr_from_csr",
            &[b, Value::int(1), Value::int(1)]
        )
        .is_err()
    );
    assert!(
        call(
            &s,
            "num_linear_array_new",
            &[ints(&[usize::MAX, 2]), floats(&[]), Value::Float(1.0)]
        )
        .is_err()
    );
    let foreign = Session::new();
    let local = array(&s, &[1], &[1.0], Value::Float(1.0));
    assert!(call(&foreign, "num_linear_array_norm", &[local]).is_err());
}

#[test]
fn bsr_native_block_product_keeps_rhs_shape_and_reference() {
    let s = Session::new();
    let bsr = call(
        &s,
        "num_linear_bsr_from_triplets",
        &[
            Value::int(2),
            Value::int(2),
            Value::int(2),
            Value::int(1),
            ints(&[0, 0, 1]),
            ints(&[0, 1, 1]),
            floats(&[2.0, 3.0, 4.0]),
            measured(2.0, "Time"),
        ],
    )
    .unwrap();
    let dim = DimExpr::Mul(vec![
        DimExpr::Base("Length".into()),
        DimExpr::Base("Time".into()),
    ]);
    let op = "num_linear_bsr_matmul_dense";
    for (cols, input, output) in [
        (
            3,
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0][..],
            &[14.0, 19.0, 24.0, 16.0, 20.0, 24.0][..],
        ),
        (0, &[][..], &[][..]),
    ] {
        let x = array(&s, &[2, cols], input, measured(3.0, "Length"));
        let result = call(&s, op, &[bsr.clone(), x]).unwrap();
        let data = call(&s, "num_linear_array_data", &[result]).unwrap();
        assert_eq!(field(&data, "shape"), &ints(&[2, cols]));
        assert_eq!(field(&data, "data"), &floats(output));
        let Value::Dyn(reference) = field(&data, "reference") else {
            panic!("expected product dimension");
        };
        assert_eq!(reference.dim, dim);
        assert_eq!(scalar(&reference.value), 6.0);
    }
    for (shape, data) in [(&[2][..], &[1.0, 2.0][..]), (&[1, 2][..], &[1.0, 2.0][..])] {
        let x = array(&s, shape, data, measured(3.0, "Length"));
        assert!(call(&s, op, &[bsr.clone(), x]).is_err());
    }
}

#[test]
fn sparse_direct_batch_solves_keep_quotient_reference_and_method() {
    let s = Session::new();
    for (matrix_reference, rhs_reference, expected_scale) in
        [(2.0, 6.0, 3.0), (1e-311, 1e-311, 1.0)]
    {
        let matrix = call(
            &s,
            "num_linear_csr_new",
            &[
                Value::int(2),
                Value::int(2),
                ints(&[0, 1]),
                ints(&[0, 1]),
                floats(&[2.0, 4.0]),
                measured(matrix_reference, "Length"),
            ],
        )
        .unwrap();
        let rhs = array(
            &s,
            &[2, 2],
            &[2.0, 6.0, 4.0, 12.0],
            measured(rhs_reference, "Time"),
        );
        let args = [matrix.clone(), rhs];
        for (name, automatic) in [
            ("num_linear_sparse_native_lu_solve_many", false),
            ("num_linear_sparse_native_cholesky_solve_many", false),
            ("num_linear_sparse_native_solve_auto_many", true),
        ] {
            let result = call(&s, name, &args).unwrap();
            let solution = if automatic {
                assert_eq!(
                    field(&result, "algorithm"),
                    &Value::String("Cholesky".into())
                );
                field(&result, "solution")
            } else {
                &result
            };
            let packed = call(&s, "num_linear_array_data", std::slice::from_ref(solution)).unwrap();
            assert_eq!(field(&packed, "shape"), &ints(&[2, 2]));
            let Value::Tensor(values) = field(&packed, "data") else {
                panic!("expected the packed solution")
            };
            for (actual, expected) in values.data.iter().map(scalar).zip([1.0, 3.0, 1.0, 3.0]) {
                assert!((actual - expected).abs() < 1e-14);
            }
            let reference = field(&packed, "reference");
            let Value::Dyn(reference) = reference else {
                panic!("expected a physical quotient")
            };
            let expected = DimExpr::Mul(vec![
                DimExpr::Inverse(Box::new(DimExpr::base("Length"))),
                DimExpr::base("Time"),
            ]);
            assert!(
                rssn::kernel::overload::matches_type(
                    &TypeExpr::Float.with_dim(expected),
                    &TypeExpr::Float.with_dim(reference.dim.clone()),
                )
                .accepted()
            );
            assert_eq!(scalar(&reference.value), expected_scale);
        }
        let wrong = array(
            &s,
            &[4],
            &[2.0, 6.0, 4.0, 12.0],
            measured(rhs_reference, "Time"),
        );
        assert!(
            call(
                &s,
                "num_linear_sparse_native_lu_solve_many",
                &[matrix, wrong]
            )
            .is_err()
        );
    }
}

#[test]
fn csc_products_reuse_a_measured_destination_with_its_own_reference() {
    let s = Session::new();
    let a = call(
        &s,
        "num_linear_csc_from_triplets",
        &[
            Value::int(2),
            Value::int(3),
            ints(&[0, 1, 0, 1]),
            ints(&[0, 0, 2, 1]),
            floats(&[2.0, 1.0, -1.0, 3.0]),
            measured(2.0, "Length"),
        ],
    )
    .unwrap();
    let dimension = DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time")]);
    let reference = |value| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim: dimension.clone(),
        })
    };
    let from_rows = array(&s, &[2], &[4.0, -2.0], measured(5.0, "Time"));
    let from_cols = array(&s, &[3], &[3.0, 4.0, 5.0], measured(5.0, "Time"));
    let output_rows = array(&s, &[2], &[0.0; 2], reference(4.0));
    let output_cols = array(&s, &[3], &[0.0; 3], reference(4.0));
    let expected = call(
        &s,
        "num_linear_csc_transpose_matvec",
        &[a.clone(), from_rows.clone()],
    )
    .unwrap();
    call(
        &s,
        "num_linear_csc_transpose_matvec_into",
        &[a.clone(), from_rows, output_cols.clone()],
    )
    .unwrap();
    call(
        &s,
        "num_linear_csc_matvec_into",
        &[a, from_cols, output_rows.clone()],
    )
    .unwrap();
    for (value, scale, numbers) in [
        (expected, 10.0, &[6.0, -6.0, -4.0][..]),
        (output_cols, 4.0, &[15.0, -15.0, -10.0][..]),
        (output_rows, 4.0, &[2.5, 37.5][..]),
    ] {
        let packed = call(&s, "num_linear_array_data", &[value]).unwrap();
        let actual_ref = field(&packed, "reference");
        assert_eq!(scalar(actual_ref), scale);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(dimension.clone()),
                &actual_ref.type_expr()
            )
            .accepted()
        );
        let Value::Tensor(data) = field(&packed, "data") else {
            panic!("expected packed array data")
        };
        for (&expected, actual) in numbers.iter().zip(data.data.iter().map(scalar)) {
            assert!((actual - expected).abs() < 1e-13);
        }
    }
}

#[test]
fn csc_in_place_rejects_nonfinite_result_without_partially_overwriting() {
    let s = Session::new();
    let a = call(
        &s,
        "num_linear_csc_from_triplets",
        &[
            Value::int(2),
            Value::int(2),
            ints(&[0, 1]),
            ints(&[0, 1]),
            floats(&[1.0, f64::MAX]),
            Value::Float(1.0),
        ],
    )
    .unwrap();
    let input = array(&s, &[2], &[2.0; 2], measured(1.0, "Time"));
    let output = array(&s, &[2], &[11.0, 13.0], measured(1.0, "Time"));
    for name in [
        "num_linear_csc_matvec_into",
        "num_linear_csc_transpose_matvec_into",
    ] {
        assert!(call(&s, name, &[a.clone(), input.clone(), output.clone()]).is_err());
    }
    let packed = call(&s, "num_linear_array_data", &[output]).unwrap();
    assert_eq!(field(&packed, "data"), &floats(&[11.0, 13.0]));
}

#[test]
fn csc_norms_return_three_measured_norms_in_one_call() {
    let session = Session::new();
    let matrix = call(
        &session,
        "num_linear_csc_from_triplets",
        &[
            Value::int(2),
            Value::int(3),
            ints(&[0, 0, 1]),
            ints(&[0, 2, 2]),
            floats(&[-3.0, 4.0, 12.0]),
            measured(2.0, "Mass"),
        ],
    )
    .unwrap();
    let value = call(&session, "num_linear_csc_norms", &[matrix]).unwrap();
    let expected = TypeExpr::Float.with_dim(DimExpr::base("Mass"));
    for (name, number) in [("frobenius", 26.0), ("l1", 32.0), ("linf", 24.0)] {
        let norm = field(&value, name);
        assert!(rssn::kernel::overload::matches_type(&expected, &norm.type_expr()).accepted());
        assert!((scalar(norm) - number).abs() < 1e-13);
    }
}

#[test]
fn independent_csc_handles_charge_retained_storage_per_matrix() {
    let session = Session::new();
    let matrices = (1..=24)
        .map(|value| {
            call(
                &session,
                "num_linear_csc_from_triplets",
                &[
                    Value::int(1),
                    Value::int(1),
                    ints(&[0]),
                    ints(&[0]),
                    floats(&[f64::from(value)]),
                    measured(2.0, "Length"),
                ],
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    for (index, matrix) in matrices.into_iter().enumerate() {
        let norm = call(&session, "num_linear_csc_l1_norm", &[matrix]).unwrap();
        assert_eq!(scalar(&norm), 2.0 * (index + 1) as f64);
    }
}
