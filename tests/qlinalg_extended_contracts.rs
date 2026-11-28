//! Exact matrix operators use normalized handles and preserve integer and rational outputs.
use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

fn call(
    s: &Session,
    operation: &str,
    args: &[Value],
) -> Value {
    s.call_value(operation, args, &CallCtx::new())
        .unwrap_or_else(|error| panic!("{operation}: {error:?}"))
}

fn ints<const N: usize>(
    shape: [usize; N],
    values: impl IntoIterator<Item = impl Into<BigInt>>,
) -> Value {
    Value::Tensor(TensorValue::of(
        shape.to_vec(),
        values.into_iter().map(Value::int).collect(),
        TypeExpr::Int,
    ))
}

fn matrix(
    s: &Session,
    values: Value,
) -> Value {
    call(s, "num_native_qmatrix_new", &[values, Value::int(1)])
}

fn some(value: Value) -> Value {
    let Value::Variant(VariantValue { tag, mut payload }) = value else {
        panic!("expected an option");
    };
    assert_eq!(tag, "Some");
    assert_eq!(payload.len(), 1);
    payload.pop().unwrap()
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected a normal-form record");
    };
    &record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .unwrap()
        .1
}

fn integers2(value: &Value) -> [BigInt; 4] {
    let Value::Tensor(value) = value else {
        panic!("expected an exact integer matrix");
    };
    assert_eq!(value.shape, [2, 2]);
    assert_eq!(value.elem, TypeExpr::Int);
    let values = value.data.iter().map(|v| {
        let Value::Int(n) = v else {
            panic!("expected an exact integer coefficient");
        };
        n.clone()
    });
    values.collect::<Vec<_>>().try_into().unwrap()
}

fn product2(
    a: &[BigInt; 4],
    b: &[BigInt; 4],
) -> [BigInt; 4] {
    [
        &a[0] * &b[0] + &a[1] * &b[2],
        &a[0] * &b[1] + &a[1] * &b[3],
        &a[2] * &b[0] + &a[3] * &b[2],
        &a[2] * &b[1] + &a[3] * &b[3],
    ]
}

fn unimodular(m: &[BigInt; 4]) -> bool {
    let det = &m[0] * &m[3] - &m[1] * &m[2];
    det == BigInt::from(1) || det == BigInt::from(-1)
}

#[test]
fn exact_determinants_and_normal_forms_use_integer_arithmetic() {
    let s = Session::new();
    let wide = BigInt::from(-6) << 130_u32;
    let mut coefficients = vec![BigInt::from(0); 25];
    for (i, v) in [wide.clone(), (-3).into(), 5.into(), 7.into(), 11.into()]
        .into_iter()
        .enumerate()
    {
        coefficients[6 * i] = v;
    }
    let a = matrix(&s, ints([5, 5], coefficients));
    let operations = [
        "num_native_qlinalg_det_bareiss",
        "num_native_qlinalg_det_modular_int",
    ];
    for operation in operations {
        assert_eq!(
            call(&s, operation, std::slice::from_ref(&a)),
            Value::int(-1155_i32 * &wide)
        );
    }

    let coefficients = [4, 6, 10, 8].map(BigInt::from);
    let a = matrix(&s, ints([2, 2], coefficients.clone()));
    let hermite = call(
        &s,
        "num_native_qlinalg_hermite_normal_form",
        std::slice::from_ref(&a),
    );
    assert_eq!(hermite, ints([2, 2], [2, 10, 0, 14]));
    let h = call(
        &s,
        "num_native_qlinalg_hermite_normal_form_with_transform",
        std::slice::from_ref(&a),
    );
    let u = integers2(field(&h, "transform"));
    assert!(unimodular(&u));
    assert_eq!(field(&h, "form"), &hermite);
    assert_eq!(product2(&u, &coefficients), integers2(&hermite));

    let smith = call(
        &s,
        "num_native_qlinalg_smith_normal_form",
        std::slice::from_ref(&a),
    );
    let u = integers2(field(&smith, "left"));
    let v = integers2(field(&smith, "right"));
    assert!(unimodular(&u) && unimodular(&v));
    assert_eq!(field(&smith, "form"), &ints([2, 2], [2, 0, 0, 14]));
    assert_eq!(
        product2(&product2(&u, &coefficients), &v),
        integers2(field(&smith, "form"))
    );
    assert_eq!(
        call(&s, "num_native_qlinalg_smith_normal_form_diag", &[a]),
        ints([2], [2, 14])
    );

    let wide_empty = matrix(&s, ints([0, 2], [] as [i32; 0]));
    let smith = call(&s, "num_native_qlinalg_smith_normal_form", &[wide_empty]);
    assert_eq!(field(&smith, "form"), &ints([0, 2], [] as [i32; 0]));
    assert_eq!(field(&smith, "left"), &ints([0, 0], [] as [i32; 0]));
    assert_eq!(field(&smith, "right"), &ints([2, 2], [1, 0, 0, 1]));
}

#[test]
fn exact_solvers_preserve_residuals_and_degenerate_shapes() {
    let s = Session::new();
    let a = matrix(&s, ints([2, 2], [3, 0, 0, 2]));
    let polynomial = some(call(
        &s,
        "num_native_qlinalg_charpoly",
        std::slice::from_ref(&a),
    ));
    let Value::Tensor(polynomial) = polynomial else {
        panic!("expected coefficients")
    };
    assert_eq!(polynomial.shape, [3]);
    assert_eq!(
        polynomial.data,
        vec![Value::int(6), Value::int(-5), Value::int(1)]
    );
    let exact = Value::Tensor(TensorValue::of(
        vec![2],
        vec![
            Value::rational(BigRational::new(1.into(), 3.into())),
            Value::rational(BigRational::new((-1).into(), 2.into())),
        ],
        TypeExpr::Real,
    ));
    for operation in [
        "num_native_qlinalg_solve_dixon",
        "num_native_qlinalg_solve_unique",
    ] {
        assert_eq!(
            some(call(&s, operation, &[a.clone(), ints([2], [1, -1])])),
            exact
        );
    }
    let a = matrix(&s, ints([1, 2], [2, 4]));
    let b = matrix(&s, ints([1, 2], [6, 8]));
    let x = some(call(&s, "num_native_qlinalg_solve_mat", &[a.clone(), b]));
    let Value::Tensor(x) = call(&s, "num_native_qmatrix_export", &[x]) else {
        panic!("expected exact tensor");
    };
    assert_eq!(x.shape, [2, 2]);
    assert_eq!(
        x.data,
        vec![Value::int(3), Value::int(4), Value::int(0), Value::int(0)]
    );
    let no_solution = Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    assert_eq!(
        call(&s, "num_native_qlinalg_solve_unique", &[a, ints([1], [6])]),
        no_solution,
    );
    let singular = matrix(&s, ints([2, 2], [1, 2, 2, 4]));
    assert_eq!(
        call(
            &s,
            "num_native_qlinalg_solve_dixon",
            &[singular, ints([2], [1, 2])]
        ),
        no_solution
    );

    let a = matrix(&s, ints([0, 2], [] as [i32; 0]));
    let b = matrix(&s, ints([0, 3], [] as [i32; 0]));
    let x = some(call(&s, "num_native_qlinalg_solve_mat", &[a, b]));
    let Value::Tensor(x) = call(&s, "num_native_qmatrix_export", &[x]) else {
        panic!("expected exact tensor");
    };
    assert_eq!(x.shape, [2, 3]);
    assert_eq!(x.data, vec![Value::int(0); 6]);
}

#[test]
fn exact_matrix_operators_reject_invalid_shape_fraction_and_work() {
    let s = Session::new();
    let attempt = |operation, args: &[Value]| s.call_value(operation, args, &CallCtx::new());
    let a = matrix(&s, ints([1, 2], [1, 2]));
    assert!(
        matches!(attempt("num_native_qlinalg_charpoly", std::slice::from_ref(&a)), Err(OpError::Domain(message)) if message.contains("square"))
    );
    assert!(
        matches!(attempt("num_native_qlinalg_solve_unique", &[a.clone(), ints([2], [1,2])]), Err(OpError::Domain(message)) if message.contains("row count"))
    );
    let b = matrix(&s, ints([2, 1], [1, 2]));
    assert!(
        matches!(attempt("num_native_qlinalg_solve_mat", &[a, b]), Err(OpError::Domain(message)) if message.contains("row count"))
    );

    let fraction = Value::rational(BigRational::new(1.into(), 3.into()));
    let a = matrix(
        &s,
        Value::Tensor(TensorValue::of(vec![1, 1], vec![fraction], TypeExpr::Real)),
    );
    for operation in [
        "num_native_qlinalg_det_bareiss",
        "num_native_qlinalg_smith_normal_form_diag",
    ] {
        assert!(
            matches!(attempt(operation, std::slice::from_ref(&a)), Err(OpError::Domain(message)) if message.contains("integral"))
        );
    }

    let mut coefficients = vec![BigInt::from(0); 1024];
    coefficients[0] = BigInt::from(1) << 256_u32;
    let a = matrix(&s, ints([32, 32], coefficients));
    assert!(
        matches!(attempt("num_native_qlinalg_hermite_normal_form", &[a]), Err(OpError::Domain(message)) if message.contains("work budget"))
    );
}


#[test]
fn exact_matrix_handles_are_session_bound_and_rhs_must_be_normalized() {
    let s = Session::new();
    let a = matrix(&s, ints([1, 1], [3]));
    let rhs = Value::Tensor(TensorValue::of(
        vec![1],
        vec![Value::Dyn(DynValue {
            value: Box::new(Value::int(1)),
            dim: DimExpr::base("Length"),
        })],
        TypeExpr::Int.with_dim(DimExpr::base("Length")),
    ));
    assert!(
        s.call_value(
            "num_native_qlinalg_solve_dixon",
            &[a.clone(), rhs],
            &CallCtx::new()
        )
        .is_err()
    );
    let alien = Session::new();
    assert!(
        alien
            .call_value(
                "num_native_qlinalg_charpoly",
                std::slice::from_ref(&a),
                &CallCtx::new()
            )
            .is_err()
    );
    assert_eq!(
        some(call(&s, "num_native_qlinalg_charpoly", &[a])),
        Value::Tensor(TensorValue::of(
            vec![2],
            vec![Value::int(-3), Value::int(1)],
            TypeExpr::Real
        ))
    );
}
