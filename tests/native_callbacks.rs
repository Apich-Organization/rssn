use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;
fn scalar(
    v: f64,
    d: DimExpr,
) -> Value {
    if d == DimExpr::Dimless {
        Value::Float(v)
    } else {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(v)),
            dim: d,
        })
    }
}
fn dim(s: &str) -> DimExpr {
    DimExpr::Base(s.into())
}
fn ty(d: DimExpr) -> TypeExpr {
    if d == DimExpr::Dimless {
        TypeExpr::Float
    } else {
        TypeExpr::Dim(Box::new(TypeExpr::Float), d)
    }
}
fn number(v: &Value) -> f64 {
    match v {
        | Value::Dyn(v) => number(&v.value),
        | _ => v.as_f64().unwrap(),
    }
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else { panic!() };
    r.get(n).unwrap()
}
fn session() -> rssn::Session {
    rssn::Session::new()
}

#[test]
fn indefinite_recurrences_validate_callback_work_before_sampling() {
    let session = session();
    let length = dim("Length");
    let time = dim("Time");
    let sampled = Arc::new(AtomicUsize::new(0));
    let callback_counter = Arc::clone(&sampled);
    let mut ctx = CallCtx::new();
    let measured_callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |a: &[Value]| {
                let sample = callback_counter.fetch_add(1, Relaxed) % 3;
                assert_eq!(number(&a[0]), [2.0, 4.0, 6.0][sample]);
                Ok(scalar(2.0, dim("Time")))
            },
            Signature::new(vec![ty(length.clone())], ty(time.clone())),
        )))
        .unwrap();
    let config = Value::Record(RecordValue {
        fields: vec![
            ("h".into(), scalar(2.0, length.clone())),
            ("step".into(), scalar(2.0, length.clone())),
            ("tolerance".into(), Value::Float(1e-10)),
            ("shift".into(), Value::int(20)),
            ("taylor_terms".into(), Value::int(8)),
        ],
    });
    let call = |op, budget| {
        session.call_value(
            op,
            &[
                Value::Function(measured_callback),
                scalar(8.0, length.clone()),
                config.clone(),
                scalar(2.0, length.clone()),
                scalar(1.0, time.clone()),
                Value::int(budget),
            ],
            &ctx,
        )
    };
    for (op, expected, expected_dim) in [
        ("num_native_sum_indefinite_sum", 6.0, time.clone()),
        ("num_native_sum_indefinite_product", 8.0, DimExpr::Dimless),
    ] {
        let before = sampled.load(Relaxed);
        let err = call(op, 2).unwrap_err();
        assert!(err.to_string().contains("evaluation budget"), "{err}");
        assert_eq!(sampled.load(Relaxed), before);
        let result = call(op, 3).unwrap();
        assert_eq!(sampled.load(Relaxed), before + 3);
        assert!(
            rssn::kernel::overload::matches_type(&ty(expected_dim), &result.type_expr()).accepted()
        );
        assert!((number(&result) - expected).abs() < 1e-12);
    }
}

#[test]
fn scaled_derivative_and_source_failure() {
    let s = session();
    let mut ctx = CallCtx::new();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(scalar(number(&args[0]).powi(2), dim("Time"))),
            Signature::new(vec![ty(dim("Length"))], ty(dim("Time"))),
        )))
        .unwrap();
    let result = s
        .call_ho(
            "num_native_derivative_scaled",
            &[
                scalar(2., dim("Length")),
                scalar(1., dim("Length")),
                scalar(1., dim("Time")),
                Value::int(10),
            ],
            &[callback],
            &ctx,
        )
        .unwrap();
    assert!((number(&result) - 4.).abs() < 1e-9);
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("sentinel".into())),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let err = s.call_ho(
        "num_native_derivative_scaled",
        &[
            Value::Float(2.),
            Value::Float(1.),
            Value::Float(1.),
            Value::int(10),
        ],
        &[cb],
        &ctx,
    );
    assert!(matches!(err,Err(OpError::Host(ref s))if s=="sentinel"));
}
#[test]
fn scaled_series_and_complex_circle() {
    let s = session();
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| Ok(scalar(0.5f64.powf(number(&a[0])), dim("Length"))),
            Signature::new(vec![TypeExpr::Float], ty(dim("Length"))),
        )))
        .unwrap();
    let out = s
        .call_ho(
            "num_native_sum_series",
            &[
                Value::int(0),
                Value::int(40),
                scalar(1e-10, dim("Length")),
                scalar(1., dim("Length")),
                Value::int(100),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    assert!((number(field(&out, "value")) - 2.).abs() < 1e-8);
    let mut ctx = CallCtx::new();
    let zty = |d| TypeExpr::Dim(Box::new(TypeExpr::Complex), d);
    let reciprocal = DimExpr::Inverse(Box::new(dim("Length")));
    let outdim = reciprocal.clone();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |a: &[Value]| {
                let Value::Dyn(v) = &a[0] else { panic!() };
                let Value::Complex(c) = v.value.as_ref() else {
                    panic!()
                };
                let (re, im) = c.to_f64();
                let q = re * re + im * im;
                Ok(Value::Dyn(DynValue {
                    value: Box::new(Value::Complex(ComplexValue::Float {
                        re: re / q,
                        im: -im / q,
                    })),
                    dim: outdim.clone(),
                }))
            },
            Signature::new(vec![zty(dim("Length"))], zty(reciprocal.clone())),
        )))
        .unwrap();
    let center = Value::Dyn(DynValue {
        value: Box::new(Value::Complex(ComplexValue::Float { re: 0., im: 0. })),
        dim: dim("Length"),
    });
    let out = s
        .call_ho(
            "num_complex_circle_integral_scaled",
            &[
                center,
                scalar(2., dim("Length")),
                Value::int(128),
                scalar(1., dim("Length")),
                scalar(1., reciprocal),
                Value::int(200),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    let Value::Complex(c) = out else {
        panic!("{out:?}")
    };
    let (re, im) = c.to_f64();
    assert!(re.abs() < 1e-12 && (im - 2. * std::f64::consts::PI).abs() < 1e-12);
}

fn vector(
    values: &[f64],
    d: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().map(|&x| Value::Float(x)).collect(),
        elem: ty(d),
    })
}
fn vector_type(d: DimExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: rssn::kernel::signature::NatExpr::Lit(1),
        shape: None,
        elem: Box::new(ty(d)),
    }
}
fn dual(
    v: f64,
    d: f64,
) -> Value {
    Value::Record(rssn::kernel::value::RecordValue {
        fields: vec![
            ("val".into(), Value::Float(v)),
            ("deriv".into(), Value::Float(d)),
        ],
    })
}
fn dual_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("val".into(), TypeExpr::Float),
        ("deriv".into(), TypeExpr::Float),
    ])
}
#[test]
fn forward_duals_restore_physical_derivatives_and_honor_budget() {
    let s = session();
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| {
                let x = rssn::kernels::autodiff::Dual::new(
                    number(field(&a[0], "val")),
                    number(field(&a[0], "deriv")),
                );
                let y = x * x;
                Ok(dual(y.val, y.deriv))
            },
            Signature::new(vec![dual_type()], dual_type()),
        )))
        .unwrap();
    let result = s
        .call_ho(
            "num_dual_forward_derivative_scaled",
            &[
                scalar(6., dim("Length")),
                scalar(2., dim("Length")),
                scalar(5., dim("Time")),
                Value::int(1),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    assert_eq!(number(field(&result, "value")), 45.);
    assert_eq!(number(field(&result, "derivative")), 15.);
    assert_ne!(
        rssn::kernel::overload::matches_type(
            &ty(DimExpr::Mul(vec![
                dim("Time"),
                DimExpr::Inverse(Box::new(dim("Length")))
            ])),
            &ctx.value_type(field(&result, "derivative")).unwrap(),
        ),
        rssn::kernel::overload::Match::No,
    );
    let mut ctx = CallCtx::new();
    let vector_dual = TypeExpr::Record(vec![
        ("val".into(), vector_type(DimExpr::Dimless)),
        ("deriv".into(), vector_type(DimExpr::Dimless)),
    ]);
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| {
                let Value::Tensor(v) = field(&a[0], "val") else {
                    panic!()
                };
                let Value::Tensor(d) = field(&a[0], "deriv") else {
                    panic!()
                };
                let x = rssn::kernels::autodiff::Dual::new(number(&v.data[0]), number(&d.data[0]));
                let y = rssn::kernels::autodiff::Dual::new(number(&v.data[1]), number(&d.data[1]));
                let z = x * x + y * y;
                Ok(dual(z.val, z.deriv))
            },
            Signature::new(vec![vector_dual], dual_type()),
        )))
        .unwrap();
    let result = s
        .call_ho(
            "num_dual_forward_gradient_scaled",
            &[
                vector(&[2., 3.], dim("Length")),
                scalar(1., dim("Length")),
                scalar(1., dim("Time")),
                Value::int(2),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(t) = result else { panic!() };
    assert_eq!(t.data.iter().map(number).collect::<Vec<_>>(), vec![4., 6.]);
    let error = s
        .call_ho(
            "num_dual_forward_gradient_scaled",
            &[
                vector(&[2., 3.], dim("Length")),
                scalar(1., dim("Length")),
                scalar(1., dim("Time")),
                Value::int(1),
            ],
            &[cb],
            &ctx,
        )
        .unwrap_err();
    assert!(matches!(error, OpError::BudgetExceeded));
}
#[test]
fn reverse_native_handles_are_scoped_and_restore_units() {
    let s = session();
    let mut ctx = CallCtx::new();
    let var = TypeExpr::Opaque("rssn/native/ad/var/v1".into());
    let vars = TypeExpr::Tensor {
        rank: rssn::kernel::signature::NatExpr::Lit(1),
        shape: None,
        elem: Box::new(var.clone()),
    };
    let escaped = Arc::new(std::sync::Mutex::new(None));
    let save = escaped.clone();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |a: &[Value]| {
                let Value::Tensor(v) = &a[0] else { panic!() };
                *save.lock().unwrap() = Some(v.data[0].clone());
                Ok(v.data[0].clone())
            },
            Signature::new(vec![vars], var),
        )))
        .unwrap();
    let result = s
        .call_ho(
            "num_native_reverse_gradient_scaled",
            &[
                vector(&[6., 8.], dim("Length")),
                scalar(2., dim("Length")),
                scalar(5., dim("Time")),
                Value::int(1),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    assert_eq!(number(field(&result, "value")), 15.);
    let Value::Tensor(g) = field(&result, "gradient") else {
        panic!()
    };
    assert_eq!(g.data.iter().map(number).collect::<Vec<_>>(), vec![2.5, 0.]);
    let handle = escaped.lock().unwrap().clone().unwrap();
    assert!(s.call_value("num_ad_var_value", &[handle], &ctx).is_err());
}
#[test]
fn native_simplex_optimizes_physical_objective() {
    let s = session();
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| {
                let Value::Tensor(v) = &a[0] else { panic!() };
                Ok(scalar((number(&v.data[0]) - 3.).powi(2), dim("Time")))
            },
            Signature::new(vec![vector_type(dim("Length"))], ty(dim("Time"))),
        )))
        .unwrap();
    let out = s
        .call_ho(
            "num_native_nelder_mead_scaled",
            &[
                vector(&[0.], dim("Length")),
                scalar(1., dim("Length")),
                scalar(2., dim("Time")),
                scalar(1., dim("Length")),
                Value::Float(1e-10),
                Value::int(1000),
                Value::int(10000),
            ],
            &[cb],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(x) = field(&out, "x") else {
        panic!()
    };
    assert!((number(&x.data[0]) - 3.).abs() < 1e-4);
    assert!(number(field(&out, "value")) < 1e-8);
}
