//! Expected values and work semantics for checked numerical transforms.

use num_complex::Complex64 as C;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::complex_transforms as native;
use rssn::kernels::special_complex::gamma_complex;
use std::cell::Cell;
use std::sync::Arc;

fn cx(z: C) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn real(x: f64) -> Value {
    Value::Float(x)
}

fn vector(
    data: &[C],
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data: data.iter().copied().map(cx).collect(),
        elem,
    })
}

fn as_complex(value: &Value) -> C {
    match value {
        | Value::Complex(z) => {
            let (re, im) = z.to_f64();
            C::new(re, im)
        },
        | Value::Dyn(v) => as_complex(&v.value),
        | _ => C::new(value.as_f64().unwrap(), 0.0),
    }
}

fn field<'a>(
    value: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected transform status record")
    };
    record.get(key).unwrap()
}

fn near(
    actual: C,
    expected: C,
    tolerance: f64,
) {
    assert!(
        (actual - expected).norm() <= tolerance,
        "{actual} vs {expected}"
    );
}

#[test]
fn backward_borel_recovers_terms_lost_to_forward_weight_underflow() {
    let mut coefficients = vec![C::new(0.0, 0.0); 200];
    coefficients[199] = C::new(1e300, 0.0);
    let result = native::borel_transform_series(&coefficients, C::new(1.0, 0.0));
    let expected = (1e300_f64.ln() - (1..=199).map(|i| f64::from(i).ln()).sum::<f64>()).exp();
    near(result / expected, C::new(1.0, 0.0), 1e-12);
    let edge = [C::new(0.0, 0.0), C::new(1.7e308, 0.0), C::new(1.7e308, 0.0)];
    let restored = native::borel_transform_series(&edge, C::new(0.75, 0.0));
    near(restored / 1.7e308, C::new(1.03125, 0.0), 1e-14);
}

#[test]
fn pure_transform_routes_preserve_coefficients_and_cancel_large_gamma_products() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let time = DimExpr::base("Time");
    let coefficients = vector(
        &[C::new(1.0, 0.0), C::new(2.0, 0.0), C::new(6.0, 0.0)],
        TypeExpr::Complex.with_dim(time.clone()),
    );
    let borel = session
        .call_value(
            "num_transform_borel_series",
            &[coefficients, cx(C::new(0.5, 0.0))],
            &ctx,
        )
        .unwrap();
    near(as_complex(&borel), C::new(2.75, 0.0), 1e-14);
    assert!(matches!(borel, Value::Dyn(v) if v.dim == time));

    let v = |values: &[C]| vector(values, TypeExpr::Complex);
    let a = C::new(200.0, 0.0);
    let zero = C::new(0.0, 0.0);
    let z = C::new(1.3, 0.2);
    let result = session
        .call_value(
            "num_transform_mellin_barnes_checked",
            &[
                v(&[a, a, zero, a]),
                v(&[]),
                v(&[a, a, a]),
                v(&[]),
                cx(z),
                real(0.5),
                real(25.0),
                real(1e-9),
                Value::int(256),
            ],
            &ctx,
        )
        .unwrap();
    assert_eq!(field(&result, "converged"), &Value::Bool(true));
    near(as_complex(field(&result, "value")), (-z).exp(), 1e-9);
    let left_of_pole = native::mellin_barnes_integral(&[zero], &[], &[], &[], z, -0.5, 25.0, 1e-9);
    near(left_of_pole, (-z).exp() - 1.0, 1e-9);
    assert!(
        !native::mellin_barnes_integral(&[zero], &[], &[], &[], z, 0.0, 25.0, 1e-9).is_finite()
    );
    let meijer = session
        .call_value(
            "num_transform_meijer_g_checked",
            &[
                v(&[]),
                v(&[zero]),
                Value::int(1),
                Value::int(0),
                cx(z),
                real(0.5),
                real(25.0),
                real(1e-9),
                Value::int(256),
            ],
            &ctx,
        )
        .unwrap();
    near(as_complex(field(&meijer, "value")), (-z).exp(), 1e-9);
    assert!(
        session
            .call_value(
                "num_transform_meijer_g_checked",
                &[
                    v(&[zero]),
                    v(&[zero]),
                    Value::int(1),
                    Value::int(1),
                    cx(z),
                    real(2.0),
                    real(25.0),
                    real(1e-9),
                    Value::int(256)
                ],
                &ctx,
            )
            .is_err()
    );
}

#[test]
fn adaptive_routes_replay_independent_closed_forms() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let sig = Signature::new(vec![TypeExpr::Float], TypeExpr::Complex);
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(cx(C::new((-args[0].as_f64().unwrap()).exp(), 0.0))),
            sig.clone(),
        )))
        .unwrap();
    let gaussian = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                Ok(cx(C::new(
                    (-0.5 * args[0].as_f64().unwrap().powi(2)).exp(),
                    0.0,
                )))
            },
            sig.clone(),
        )))
        .unwrap();
    let constant = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(C::new(1.0, 0.0))),
            sig,
        )))
        .unwrap();
    let gamma = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(cx(gamma_complex(as_complex(&args[0])))),
            Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex),
        )))
        .unwrap();
    let linear = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(real(args[0].as_f64().unwrap())),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let s = C::new(1.5, 0.4);
    let cases = [
        (
            "mellin",
            vec![cx(s), real(20.0)],
            callback,
            gamma_complex(s),
        ),
        (
            "inverse_mellin",
            vec![cx(C::new(0.7, 0.0)), real(1.0), real(25.0)],
            gamma,
            C::new((-0.7_f64).exp(), 0.0),
        ),
        (
            "hilbert",
            vec![real(0.0), real(1.0)],
            linear,
            C::new(-2.0 / std::f64::consts::PI, 0.0),
        ),
        (
            "hankel",
            vec![real(0.0), real(1.2), real(12.0)],
            gaussian,
            C::new((-0.72_f64).exp(), 0.0),
        ),
        (
            "stieltjes",
            vec![cx(C::new(1.0, 0.0)), real(2.0)],
            constant,
            C::new(3.0_f64.ln(), 0.0),
        ),
    ];
    for (op, mut args, callback, expected) in cases {
        args.extend([real(1e-10), Value::int(256), Value::int(2000)]);
        let result = session
            .call_ho(
                &format!("num_transform_{op}_checked"),
                &args,
                &[callback],
                &ctx,
            )
            .unwrap();
        near(as_complex(field(&result, "value")), expected, 1e-8);
        assert_eq!(field(&result, "converged"), &Value::Bool(true));
    }
}

#[test]
fn laplace_variants_recover_real_and_complex_poles_with_observed_counts() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let pole = C::new(-1.25, 0.4);
    let complex = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| Ok(cx(C::new(1.0, 0.0) / (as_complex(&args[0]) - pole))),
            Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex),
        )))
        .unwrap();
    let real_cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(real(1.0 / (args[0].as_f64().unwrap() + 1.5))),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    for (op, shift, n, tol, calls) in [
        ("talbot", None, 32, 1e-6, 48),
        ("talbot_shifted", Some(1.0), 32, 1e-6, 48),
        ("parabola", None, 40, 1e-5, 62),
        ("parabola_shifted", Some(1.0), 40, 1e-5, 62),
    ] {
        let mut args = vec![real(0.75)];
        args.extend(shift.map(real));
        args.extend([real(tol), Value::int(n), Value::int(calls)]);
        let value = session
            .call_ho(
                &format!("num_transform_laplace_{op}_checked"),
                &args,
                &[complex],
                &ctx,
            )
            .unwrap();
        near(as_complex(field(&value, "value")), (pole * 0.75).exp(), tol);
        assert_eq!(field(&value, "evaluations"), &Value::int(calls));
    }
    let stehfest = session
        .call_ho(
            "num_transform_laplace_stehfest_checked",
            &[real(1.7), real(5e-4), Value::int(12), Value::int(22)],
            &[real_cb],
            &ctx,
        )
        .unwrap();
    near(
        as_complex(field(&stehfest, "value")),
        C::new((-1.5_f64 * 1.7).exp(), 0.0),
        1e-4,
    );
    assert_eq!(field(&stehfest, "evaluations"), &Value::int(22));
}

#[test]
fn complex_quadrature_samples_once_and_reports_real_host_budget() {
    let calls = Cell::new(0);
    let value = native::mellin_transform(
        |_x| {
            calls.set(calls.get() + 1);
            C::new(2.0, 1.0)
        },
        C::new(0.0, 0.0),
        1.0,
        1e-10,
    );
    near(value, C::new(4.0, 2.0), 1e-12);
    assert_eq!(calls.get(), 15);
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(C::new(2.0, 1.0))),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Complex),
        )))
        .unwrap();
    let result = |budget, s, limit| {
        session.call_ho(
            "num_transform_mellin_checked",
            &[
                cx(s),
                real(limit),
                real(1e-10),
                Value::int(1),
                Value::int(budget),
            ],
            &[cb],
            &ctx,
        )
    };
    assert!(matches!(
        result(14, C::new(0.0, 0.0), 1.0),
        Err(OpError::BudgetExceeded)
    ));
    let result = result(15, C::new(0.0, 10.0), 1.0).unwrap();
    assert_eq!(field(&result, "evaluations"), &Value::int(15));
    assert_eq!(field(&result, "converged"), &Value::Bool(false));
}
