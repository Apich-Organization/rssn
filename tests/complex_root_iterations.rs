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
use rssn::kernel::value::DynValue;
use rssn::kernels::complex_roots::delves_lyness_roots;
use rssn::kernels::complex_roots::halley_root;
use rssn::kernels::complex_roots::homotopy_track_1d;
use rssn::kernels::complex_roots::muller_root;
use std::cell::Cell;
use std::sync::Arc;

fn complex(z: C) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn measured(
    value: Value,
    dim: &DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim: dim.clone(),
    })
}

fn trace(value: &Value) -> C {
    match value {
        | Value::Complex(z) => {
            let (re, im) = z.to_f64();
            C::new(re, im)
        },
        | Value::Dyn(z) => trace(&z.value),
        | _ => panic!("expected measured complex coordinate"),
    }
}

fn field<'a>(
    result: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = result else {
        panic!("expected root report")
    };
    record.get(name).unwrap()
}

#[test]
fn muller_reuses_two_samples_and_rejects_nonroot_duplicates() {
    let calls = Cell::new(0);
    let f = |z: C| {
        calls.set(calls.get() + 1);
        z * z - 2.0
    };
    let root = muller_root(
        f,
        C::new(1.0, 0.0),
        C::new(1.2, 0.0),
        C::new(1.4, 0.0),
        1,
        1e-13,
    )
    .unwrap();
    assert!((root.re - 2.0_f64.sqrt()).abs() < 1e-13);
    assert_eq!(calls.get(), 4);
    let bad_seed = C::new(1.1, 0.1);
    assert!(muller_root(f, bad_seed, bad_seed, bad_seed, 12, 1e-13).is_err());
    assert!(
        muller_root(
            |_| C::new(f64::NAN, 0.0),
            bad_seed,
            bad_seed,
            bad_seed,
            12,
            1e-13
        )
        .is_err()
    );
}

#[test]
fn halley_scales_large_ordinate_values_and_keeps_residual_semantics() {
    let pure = |z: C| z * z - 2.0;
    let root = halley_root(
        |z| 1e200 * pure(z),
        |z| 2e200 * z,
        |_| C::new(2e200, 0.0),
        C::new(1.4, 0.1),
        12,
        1e188,
    )
    .unwrap();
    assert!(pure(root).norm() < 1e-12);
    let impossible = halley_root(
        |_| C::new(1.0, 0.0),
        |_| C::new(1e20, 0.0),
        |_| C::new(1e20, 0.0),
        C::new(0.0, 0.0),
        4,
        1e-10,
    );
    assert!(impossible.is_err());
}

#[test]
fn measured_muller_reports_exact_callback_budget_and_host_diagnostic() {
    let session = Session::new();
    let axis = DimExpr::base("Length");
    let amplitude = DimExpr::base("Time");
    let callback_type = Signature::new(
        vec![TypeExpr::Complex.with_dim(axis.clone())],
        TypeExpr::Complex.with_dim(amplitude.clone()),
    );
    let output = amplitude.clone();
    let mut ctx = CallCtx::new();
    let f = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                let z = trace(&args[0]);
                Ok(measured(complex(z * z - 2.0), &output))
            },
            callback_type.clone(),
        )))
        .unwrap();
    let throws = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| {
                Err(OpError::HostDiagnostic {
                    code: "E-test",
                    message: "root fixture failed".into(),
                })
            },
            callback_type,
        )))
        .unwrap();
    let args = |budget| {
        [
            measured(complex(C::new(1.0, 0.0)), &axis),
            measured(complex(C::new(1.2, 0.0)), &axis),
            measured(complex(C::new(1.4, 0.0)), &axis),
            measured(Value::Float(2.0), &axis),
            measured(Value::Float(5.0), &amplitude),
            Value::Float(1e-12),
            Value::int(1),
            Value::int(budget),
        ]
    };
    let result = session
        .call_ho("num_complex_root_muller", &args(4), &[f], &ctx)
        .unwrap();
    assert!((trace(field(&result, "root")).re - 2.0_f64.sqrt()).abs() < 1e-12);
    assert_eq!(field(&result, "evaluations"), &Value::int(4));
    assert!(matches!(
        session.call_ho("num_complex_root_muller", &args(3), &[f], &ctx),
        Err(OpError::BudgetExceeded)
    ));
    assert!(
        matches!(session.call_ho("num_complex_root_muller", &args(4), &[throws], &ctx), Err(OpError::HostDiagnostic { code, .. }) if code == "E-test")
    );
}

#[test]
fn homotopy_checks_both_endpoint_residuals() {
    let h = |z: C, t| z - (1.0 + 2.0 * t);
    let hz = |_, _| C::new(1.0, 0.0);
    let ht = |_, _| C::new(-2.0, 0.0);
    let end = homotopy_track_1d(h, hz, ht, C::new(1.0, 0.0), 8, 3, 1e-12).unwrap();
    assert!((end.re - 3.0).abs() < 1e-12);
    assert!(homotopy_track_1d(h, hz, ht, C::new(1.001, 0.0), 8, 3, 1e-12).is_err());
    assert!(homotopy_track_1d(h, hz, ht, C::new(1.0, 0.0), 8, 0, 1e-12).is_err());
}

#[test]
fn delves_lyness_uses_one_callback_sample_pair_per_circle_point() {
    let values = Cell::new(0);
    let slopes = Cell::new(0);
    let f = |z: C| {
        values.set(values.get() + 1);
        (z - 1.0) * (z + 1.0) * z.exp()
    };
    let df = |z: C| {
        slopes.set(slopes.get() + 1);
        (2.0 * z + z * z - 1.0) * z.exp()
    };
    let roots = delves_lyness_roots(f, df, C::new(0.0, 0.0), 1.7, 96).unwrap();
    assert_eq!(roots.len(), 2);
    for expected in [-1.0, 1.0] {
        assert!(roots.iter().any(|z| (*z - expected).norm() < 1e-11));
    }
    assert!((96..112).contains(&values.get()) && (96..112).contains(&slopes.get()));
}
