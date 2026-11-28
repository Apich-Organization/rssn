use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::riemann_hilbert::Matrix2;

fn complex(
    value: f64,
    dim: &DimExpr,
) -> Value {
    let cell = Value::Complex(ComplexValue::Float { re: value, im: 0.0 });
    if *dim == DimExpr::Dimless {
        cell
    } else {
        Value::Dyn(DynValue {
            value: Box::new(cell),
            dim: dim.clone(),
        })
    }
}

fn tensor<const N: usize>(
    shape: Vec<usize>,
    entries: [f64; N],
    dim: &DimExpr,
) -> Value {
    Value::Tensor(TensorValue::of(
        shape,
        entries
            .into_iter()
            .map(|value| complex(value, dim))
            .collect(),
        TypeExpr::Complex.with_dim(dim.clone()),
    ))
}

fn matrix_entry(
    value: &Value,
    index: usize,
) -> (f64, f64) {
    let Value::Tensor(matrix) = value else {
        panic!("expected matrix")
    };
    assert_eq!(matrix.shape, [2, 2]);
    let cell = match &matrix.data[index] {
        | Value::Dyn(v) => v.value.as_ref(),
        | v => v,
    };
    let Value::Complex(cell) = cell else {
        panic!("expected complex matrix entry")
    };
    cell.to_f64()
}

fn call(
    session: &Session,
    op: &str,
    args: &[Value],
) -> Result<Value, rssn::kernel::OpError> {
    session.call_value(&format!("num_riemann_hilbert_{op}"), args, &CallCtx::new())
}

#[test]
fn fixed_matrix_arithmetic_preserves_product_and_reciprocal_units() {
    let session = Session::new();
    let dim = DimExpr::base("Length");
    let dimless = DimExpr::Dimless;
    let entries = tensor(vec![4], [2.0, 1.0, 0.0, 4.0], &dim);
    let a = call(&session, "matrix_new", &[entries]).unwrap();
    assert_eq!(matrix_entry(&a, 1), (1.0, 0.0));
    let sum = call(&session, "matrix_add", &[a.clone(), a.clone()]).unwrap();
    assert_eq!(matrix_entry(&sum, 1), (2.0, 0.0));
    let delta = call(&session, "matrix_sub", &[sum, a.clone()]).unwrap();
    assert_eq!(delta, a);
    let inverse = call(&session, "matrix_inv", std::slice::from_ref(&a)).unwrap();
    let Value::Variant(inverse) = inverse else {
        panic!("expected option")
    };
    assert_eq!(inverse.tag.as_str(), "Some");
    let product = call(
        &session,
        "matrix_mul",
        &[a.clone(), inverse.payload[0].clone()],
    )
    .unwrap();
    assert_eq!(product, call(&session, "matrix_identity", &[]).unwrap());
    let scaled = call(
        &session,
        "matrix_scale",
        &[a.clone(), complex(3.0, &DimExpr::base("Time"))],
    )
    .unwrap();
    let Value::Tensor(scaled) = scaled else {
        panic!("expected fixed matrix")
    };
    let combined = DimExpr::Mul(vec![dim.clone(), DimExpr::base("Time")]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Complex.with_dim(combined), &scaled.elem)
            .accepted()
    );
    assert!(
        call(
            &session,
            "matrix_add",
            &[a.clone(), call(&session, "matrix_zero", &[]).unwrap()]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "matrix_norm",
            &[tensor(vec![1, 4], [1.0, 0.0, 0.0, 1.0], &dimless)]
        )
        .is_err()
    );
}

#[test]
fn inverse_and_norm_stay_representable_across_binary64_scales() {
    for scale in [1e-200, 1.0, 1e200] {
        let matrix = Matrix2::identity().scale(Complex64::new(scale, 0.0));
        assert!((matrix.norm() / scale - 2.0_f64.sqrt()).abs() < 1e-15);
        let inverse = matrix.inv().unwrap();
        assert!(matrix.mul(inverse).sub(Matrix2::identity()).norm() < 1e-15);
    }
    let diagonal = Matrix2::new(1e200.into(), 0.0.into(), 0.0.into(), 1e-200.into());
    assert!(
        diagonal
            .mul(diagonal.inv().unwrap())
            .sub(Matrix2::identity())
            .norm()
            < 1e-15
    );
    assert!(
        Matrix2::new(1.0.into(), 2.0.into(), 1.0.into(), 2.0.into())
            .inv()
            .is_none()
    );
    let unrepresentable = Matrix2::identity().scale(Complex64::new(f64::from_bits(1), 0.0));
    assert!(unrepresentable.inv().is_none());
}

fn constant_density(nodes: [f64; 3]) -> Value {
    let dimless = DimExpr::Dimless;
    Value::Record(RecordValue {
        fields: vec![
            (
                "nodes".into(),
                Value::Tensor(TensorValue::vector(
                    nodes.into_iter().map(Value::Float).collect(),
                )),
            ),
            (
                "weights".into(),
                Value::Tensor(TensorValue::vector(
                    [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
                        .into_iter()
                        .map(Value::Float)
                        .collect(),
                )),
            ),
            (
                "density".into(),
                tensor(
                    vec![3, 2, 2],
                    [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    &dimless,
                ),
            ),
        ],
    })
}

#[test]
fn the_same_collocation_record_handles_near_and_far_evaluations() {
    let session = Session::new();
    let node = 0.6_f64.sqrt();
    let solution = constant_density([-node, 0.0, node]);
    let at = |y| {
        call(
            &session,
            "solution_eval",
            &[
                solution.clone(),
                Value::Complex(ComplexValue::Float { re: 0.0, im: y }),
            ],
        )
    };
    let near = at(1e-320).unwrap();
    let lower = at(-1e-320).unwrap();
    assert!((matrix_entry(&near, 1).0 - 0.5).abs() < 1e-14);
    assert!((matrix_entry(&lower, 1).0 + 0.5).abs() < 1e-14);
    let far = at(1e200).unwrap();
    assert!((matrix_entry(&far, 1).0 * 1e200 - 1.0 / std::f64::consts::PI).abs() < 1e-14);
    assert!(at(0.0).is_err());
    let malformed = constant_density([0.0, 0.0, node]);
    assert!(
        call(
            &session,
            "solution_eval",
            &[
                malformed,
                Value::Complex(ComplexValue::Float { re: 0.0, im: 1.0 })
            ]
        )
        .is_err()
    );
}

#[test]
fn callback_solver_reuses_both_rhs_and_returns_the_analytic_nilpotent_solution() {
    use rssn::kernel::NativeCallback;
    use rssn::kernel::OpError;
    use rssn::kernel::signature::NatExpr;
    use rssn::kernel::signature::Signature;
    use std::sync::Arc;
    use std::sync::atomic::AtomicUsize;
    use std::sync::atomic::Ordering::Relaxed;

    let session = Session::new();
    let mut ctx = CallCtx::new();
    let evaluations = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&evaluations);
    let jump_type = TypeExpr::Tensor {
        rank: NatExpr::Lit(2),
        shape: Some(vec![NatExpr::Lit(2); 2]),
        elem: Box::new(TypeExpr::Complex),
    };
    let jump = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                counter.fetch_add(1, Relaxed);
                let t = args[0].as_f64().unwrap();
                let u = 1.0 - t * t;
                Ok(tensor(
                    vec![2, 2],
                    [1.0 + 2.0 * u, -u, 4.0 * u, 1.0 - 2.0 * u],
                    &DimExpr::Dimless,
                ))
            },
            Signature::new(vec![TypeExpr::Float], jump_type.clone()),
        )))
        .unwrap();
    let bad_endpoints = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(tensor(vec![2, 2], [2.0, 0.0, 0.0, 2.0], &DimExpr::Dimless)),
            Signature::new(vec![TypeExpr::Float], jump_type),
        )))
        .unwrap();
    let solve = |limit, callback| {
        session.call_ho(
            "num_riemann_hilbert_solve",
            &[Value::int(16), Value::int(limit)],
            &[callback],
            &ctx,
        )
    };
    assert!(matches!(solve(17, jump), Err(OpError::BudgetExceeded)));
    assert_eq!(evaluations.load(Relaxed), 0);
    assert!(matches!(solve(18, bad_endpoints), Err(OpError::Domain(_))));
    let value = solve(18, jump).unwrap();
    let Value::Record(solved) = value else {
        panic!("expected solver result")
    };
    assert_eq!(solved.get("evaluations"), Some(&Value::int(18)));
    assert_eq!(evaluations.load(Relaxed), 18);
    let point = Value::Complex(ComplexValue::Float { re: 0.0, im: 1.0 });
    let matrix = call(
        &session,
        "solution_eval",
        &[solved.get("solution").unwrap().clone(), point],
    )
    .unwrap();
    let coefficient = (std::f64::consts::PI - 2.0) / std::f64::consts::TAU;
    let expected = [
        1.0 + 2.0 * coefficient,
        -coefficient,
        4.0 * coefficient,
        1.0 - 2.0 * coefficient,
    ];
    for (i, wanted) in expected.into_iter().enumerate() {
        let (re, im) = matrix_entry(&matrix, i);
        assert!(
            (re - wanted).abs() < 1e-9 && im.abs() < 1e-10,
            "wrong entry {i}: {re} + {im}i"
        );
    }
}
