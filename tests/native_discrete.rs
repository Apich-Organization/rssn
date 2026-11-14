//! Concrete source-inventory objects and callbacks execute through session APIs.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
fn call(
    s: &Session,
    n: &str,
    args: &[Value],
) -> Value {
    s.call_value(n, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{n}: {e:?}"))
}
fn dim(n: &str) -> DimExpr {
    DimExpr::Base(n.into())
}
fn scalar(
    x: f64,
    d: DimExpr,
) -> Value {
    if d == DimExpr::Dimless {
        Value::Float(x)
    } else {
        Value::Dyn(DynValue {
            dim: d,
            value: Box::new(Value::Float(x)),
        })
    }
}
fn number(x: &Value) -> f64 {
    match x {
        | Value::Dyn(x) => number(&x.value),
        | _ => x.as_f64().unwrap(),
    }
}
fn tensor(
    shape: Vec<usize>,
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue { shape, data, elem })
}
fn ints(data: &[i64]) -> Value {
    tensor(
        vec![data.len()],
        data.iter().copied().map(Value::int).collect(),
        TypeExpr::Int,
    )
}
fn floats(
    data: &[f64],
    d: DimExpr,
) -> Value {
    tensor(
        vec![data.len()],
        data.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float.with_dim(d),
    )
}
fn some(x: &Value) -> &Value {
    let Value::Variant(v) = x else {
        panic!("option")
    };
    assert_eq!(v.tag, "Some");
    &v.payload[0]
}
#[test]
fn native_group_and_exact_lie_objects_are_actual_checked_algorithms() {
    let s = Session::new();
    let generators = tensor(
        vec![2, 3],
        vec![1, 0, 2, 1, 2, 0].into_iter().map(Value::int).collect(),
        TypeExpr::Int,
    );
    let g = call(&s, "num_native_group_new", &[Value::int(3), generators]);
    assert_eq!(
        call(&s, "num_native_group_order", &[g.clone()]),
        Value::int(6)
    );
    let center = call(&s, "num_native_group_center", &[g.clone()]);
    assert_eq!(call(&s, "num_native_group_order", &[center]), Value::int(1));
    assert!(
        s.call_value(
            "num_native_group_sylow_subgroup",
            &[g, Value::int(1), Value::int(1), Value::int(10)],
            &CallCtx::new()
        )
        .is_err()
    );
    let c = call(
        &s,
        "num_native_lie_parse_type",
        &[Value::String("A2".into())],
    );
    let r = call(&s, "num_native_root_system_new", &[some(&c).clone()]);
    let r = some(&r).clone();
    assert_eq!(
        call(&s, "num_native_root_system_dimension", &[r.clone()]),
        Value::int(8)
    );
    assert_eq!(
        some(&call(
            &s,
            "num_native_root_system_weyl_dim",
            &[r.clone(), ints(&[1, 0])]
        )),
        &Value::int(3)
    );
    assert_eq!(
        call(&s, "num_native_root_system_weyl_order", &[r]),
        Value::int(6)
    );
}
#[test]
fn native_geometry_frames_and_ball_witnesses_retain_physical_dimensions() {
    let s = Session::new();
    let m = |x| scalar(x, dim("Length"));
    let t = |x| scalar(x, dim("Time"));
    let p = call(&s, "num_native_graphics_point3_new", &[m(1.), m(2.), m(3.)]);
    let transform = call(
        &s,
        "num_native_graphics_translation_matrix",
        &[m(4.), m(0.), m(0.), m(2.)],
    );
    let p = call(&s, "num_native_graphics_transform_point", &[transform, p]);
    let v = call(&s, "num_native_graphics_point3_components", &[p]);
    let Value::Tensor(v) = v else {
        panic!("tensor")
    };
    assert_eq!(v.elem, TypeExpr::Float.with_dim(dim("Length")));
    assert_eq!(
        v.data.iter().map(number).collect::<Vec<_>>(),
        vec![5., 2., 3.]
    );
    let b = call(&s, "num_native_cball_square", &[m(2.), m(0.), m(0.25)]);
    assert_eq!(
        call(
            &s,
            "num_native_cball_contains",
            &[b.clone(), scalar(2., dim("Length"))]
        ),
        Value::Bool(true)
    );
    assert!(
        s.call_value(
            "num_native_cball_contains",
            &[b.clone(), t(2.)],
            &CallCtx::new()
        )
        .is_err()
    );
    let squared = call(&s, "num_native_cball_powi", &[b, Value::int(2)]);
    let width = call(&s, "num_native_cball_half_width", &[squared]);
    let Value::Dyn(witness) = &width else {
        panic!("measured ball width must retain its dimension witness")
    };
    let expected = TypeExpr::Float.with_dim(dim("Length").pow(rssn::kernel::signature::Exp::Rat(
        rssn::kernel::signature::Rational::int(2),
    )));
    assert!(
        rssn::kernel::overload::matches_type(
            &expected,
            &TypeExpr::Float.with_dim(witness.dim.clone())
        )
        .accepted()
    );
    assert!((number(&width) - 1.125).abs() < 1e-12);
}
#[test]
fn native_polynomial_certification_and_repeated_roots_are_scale_consistent() {
    let s = Session::new();
    let system = call(
        &s,
        "num_native_polynomial_system_new",
        &[
            floats(&[1., -2.], DimExpr::Dimless),
            tensor(
                vec![2, 1],
                vec![Value::int(2), Value::int(0)],
                TypeExpr::Int,
            ),
            ints(&[0, 2]),
            scalar(2., dim("Length")),
            scalar(3., dim("Time")),
        ],
    );
    let p = floats(&[2. * 2f64.sqrt()], dim("Length"));
    let proof = call(&s, "num_native_homotopy_certify_zero", &[system, p]);
    assert!(matches!(proof,Value::Variant(ref v)if v.tag=="Some"));
    let poly = call(
        &s,
        "num_native_polynomial_normalized_new",
        &[floats(&[1., 0., -4., 0., 4.], DimExpr::Dimless)],
    );
    let roots = call(
        &s,
        "num_native_real_roots_find_roots",
        &[
            poly,
            scalar(1e-10, dim("Length")),
            scalar(1., dim("Length")),
        ],
    );
    let Value::Tensor(roots) = roots else {
        panic!("roots")
    };
    assert_eq!(roots.data.len(), 2);
    assert!((number(&roots.data[0]) + 2f64.sqrt()).abs() < 1e-9);
    assert!((number(&roots.data[1]) - 2f64.sqrt()).abs() < 1e-9);
}
#[test]
fn native_parallel_callbacks_run_on_owner_thread_preserve_units_and_failures() {
    let s = Session::new();
    let mut ctx = CallCtx::new();
    let owner = std::thread::current().id();
    let calls = Arc::new(AtomicUsize::new(0));
    let seen = calls.clone();
    let h = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                assert_eq!(std::thread::current().id(), owner);
                seen.fetch_add(1, Ordering::Relaxed);
                Ok(scalar(args[0].as_f64().unwrap() + 1., dim("Length")))
            },
            Signature::new(vec![TypeExpr::Int], TypeExpr::Float.with_dim(dim("Length"))),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_par_map_indexed",
            &[
                Value::int(100),
                Value::int(2),
                Value::Function(h),
                Value::int(100),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(out) = out else {
        panic!("tensor")
    };
    assert_eq!(out.data.len(), 100);
    assert_eq!(number(&out.data[99]), 100.);
    assert_eq!(calls.load(Ordering::Relaxed), 100);
    let err = s
        .call_value(
            "num_native_par_sum_fixed",
            &[
                Value::int(100),
                Value::int(4),
                Value::int(2),
                Value::Function(h),
                Value::int(2),
            ],
            &ctx,
        )
        .unwrap_err();
    assert!(err.to_string().contains("budget"));
    let mut failed = CallCtx::new();
    let h = failed
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("parallel-sentinel".into())),
            Signature::new(vec![TypeExpr::Int], TypeExpr::Float),
        )))
        .unwrap();
    let err = s
        .call_value(
            "num_native_par_map_indexed",
            &[
                Value::int(100),
                Value::int(1),
                Value::Function(h),
                Value::int(100),
            ],
            &failed,
        )
        .unwrap_err();
    assert!(err.to_string().contains("parallel-sentinel"));
}
#[test]
fn native_metric_and_bootstrap_callbacks_keep_independent_input_output_dimensions() {
    let s = Session::new();
    let mut ctx = CallCtx::new();
    let vector_ty = |d| {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Float.with_dim(d)),
        }
    };
    let statistic = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let Value::Tensor(t) = &args[0] else {
                    panic!("sample")
                };
                Ok(scalar(
                    t.data.iter().map(number).sum::<f64>() / t.data.len() as f64,
                    dim("Time"),
                ))
            },
            Signature::new(
                vec![vector_ty(dim("Length"))],
                TypeExpr::Float.with_dim(dim("Time")),
            ),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_robust_bootstrap_ci",
            &[
                floats(&[2., 2., 2.], dim("Length")),
                Value::Function(statistic),
                Value::int(8),
                Value::Float(0.1),
                Value::int(1),
                Value::int(8),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(out) = out else {
        panic!("CI")
    };
    assert_eq!(out.elem, TypeExpr::Float.with_dim(dim("Time")));
    assert_eq!(number(&out.data[0]), 2.);
    let metric = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| {
                Ok(tensor(
                    vec![2, 2],
                    vec![
                        Value::Float(1.),
                        Value::Float(0.),
                        Value::Float(0.),
                        Value::Float(1.),
                    ],
                    TypeExpr::Float,
                ))
            },
            Signature::new(
                vec![vector_ty(dim("Length"))],
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(2),
                    shape: None,
                    elem: Box::new(TypeExpr::Float),
                },
            ),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_geometry_ricci_scalar",
            &[
                Value::Function(metric),
                floats(&[1., 2.], dim("Length")),
                scalar(1., dim("Length")),
                Value::Float(1.),
                Value::int(1000),
            ],
            &ctx,
        )
        .unwrap();
    assert!(number(&out).abs() < 1e-9);
}
#[test]
fn native_rng_and_online_covariance_are_mutable_session_objects() {
    let s = Session::new();
    let a = call(&s, "num_native_rng_new", &[Value::int(42)]);
    let b = call(&s, "num_native_rng_new", &[Value::int(42)]);
    assert_eq!(
        call(&s, "num_native_rng_next_u64", &[a.clone()]),
        call(&s, "num_native_rng_next_u64", &[b])
    );
    assert!(
        Session::new()
            .call_value("num_native_rng_uniform", &[a], &CallCtx::new())
            .is_err()
    );
    let c = call(
        &s,
        "num_native_online_covariance_new",
        &[Value::int(2), scalar(1., dim("Length"))],
    );
    for p in [[1., 3.], [3., 7.]] {
        call(
            &s,
            "num_native_online_covariance_push",
            &[c.clone(), floats(&p, dim("Length"))],
        );
    }
    assert_eq!(
        call(&s, "num_native_online_covariance_count", &[c.clone()]),
        Value::int(2)
    );
    let cov = call(&s, "num_native_online_covariance_sample_covariance", &[c]);
    let Value::Tensor(cov) = some(&cov) else {
        panic!("covariance")
    };
    assert_eq!(
        cov.data.iter().map(number).collect::<Vec<_>>(),
        vec![2., 4., 4., 8.]
    );
}
#[test]
fn general_multivector_inverse_covers_nonhomogeneous_and_tiny_elements() {
    use rssn::kernels::geometric_algebra::Multivector3D;
    let a = Multivector3D {
        s: 2.,
        v1: 1.,
        ..Default::default()
    };
    let inverse = a.inv().unwrap();
    assert!((inverse.s - 2. / 3.).abs() < 1e-15);
    assert!((inverse.v1 + 1. / 3.).abs() < 1e-15);
    assert!(
        Multivector3D {
            s: 1e-300,
            ..Default::default()
        }
        .inv()
        .is_some()
    );
    assert!(
        Multivector3D {
            s: 1.,
            v1: 1.,
            ..Default::default()
        }
        .inv()
        .is_none()
    );
}
