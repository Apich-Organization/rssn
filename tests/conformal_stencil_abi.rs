//! Accuracy, finite sampling and dimensional regressions for complex stencils.
use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use num_complex::Complex64;
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
use rssn::kernels::complex::check_cauchy_riemann_numeric;
use rssn::kernels::complex::complex_laplacian_numeric;
use rssn::kernels::complex::complex_step_derivative;
use rssn::kernels::complex::wirtinger_dz_numeric;
use rssn::kernels::complex::wirtinger_dzbar_numeric;
use rssn::kernels::conformal::is_conformal_numeric;

#[test]
fn holomorphic_diagnostics_reuse_four_complex_samples() {
    let count = Cell::new(0);
    let f = |z| {
        count.set(count.get() + 1);
        Complex64::new(2.0, -3.0) * z
    };
    let z = Complex64::new(1.0, 2.0);
    assert!(is_conformal_numeric(&f, z, 0.125, 1e-12));
    assert_eq!(count.replace(0), 4);
    assert!(check_cauchy_riemann_numeric(f, z, 0.125, 1e-12));
    assert_eq!(count.replace(0), 4);
    assert!(!check_cauchy_riemann_numeric(|z| z.conj(), z, 0.125, 1e-12));
    assert!(!is_conformal_numeric(
        &|_| Complex64::default(),
        z,
        0.125,
        1e-12
    ));
}

#[test]
fn extreme_spans_and_unrepresentable_perturbations_are_distinguished() {
    let z = Complex64::default();
    assert_eq!(wirtinger_dz_numeric(|z| z, z, 1e308), Complex64::from(1.0));
    assert_eq!(
        wirtinger_dzbar_numeric(|z| z, z, 1e308),
        Complex64::default()
    );
    let huge = Complex64::new(1e200, 0.0);
    let calls = Cell::new(0);
    assert!(
        !wirtinger_dz_numeric(
            |z| {
                calls.set(calls.get() + 1);
                z
            },
            huge,
            1.0,
        )
        .is_finite()
    );
    assert_eq!(calls.get(), 0);
}

#[test]
fn laplacian_and_complex_step_use_exact_work_counts() {
    let calls = Cell::new(0);
    let z = Complex64::new(1.0, 2.0);
    let laplacian = complex_laplacian_numeric(
        |z| {
            calls.set(calls.get() + 1);
            Complex64::new(
                3.0 * z.re * z.re + 2.0 * z.im * z.im,
                z.re * z.re - 5.0 * z.im * z.im,
            )
        },
        z,
        0.125,
    );
    assert_eq!(laplacian, Complex64::new(10.0, -8.0));
    assert_eq!(calls.replace(0), 5);
    let step = complex_step_derivative(
        |z| {
            calls.set(calls.get() + 1);
            z * z * z
        },
        2.0,
        1e-120,
    );
    assert!((step - 12.0).abs() < 1e-12);
    assert_eq!(calls.get(), 1);
}

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}
fn measured(
    value: Value,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim,
    })
}

#[test]
fn declared_stencils_preserve_units_and_work_bounds() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let calls = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&calls);
    let out = time.clone();
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                let [Value::Dyn(arg)] = args else {
                    panic!("expected measured position")
                };
                let Value::Complex(z) = arg.value.as_ref() else {
                    panic!("expected complex scalar")
                };
                counter.fetch_add(1, Ordering::Relaxed);
                let (re, im) = z.to_f64();
                Ok(measured(complex(2.0 * re, 2.0 * im), out.clone()))
            },
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();

    let args = |budget| {
        [
            measured(complex(1.0, 2.0), length.clone()),
            measured(Value::Float(0.125), length.clone()),
            measured(Value::Float(3.0), time.clone()),
            Value::int(budget),
        ]
    };
    let result = session
        .call_ho("num_complex_wirtinger_dz_checked", &args(4), &[cb], &ctx)
        .unwrap();
    let Value::Dyn(result) = result else {
        panic!("expected measured derivative")
    };
    let Value::Complex(z) = result.value.as_ref() else {
        panic!("expected complex derivative")
    };
    let (re, im) = z.to_f64();
    assert!((re - 2.0).abs() < 1e-12 && im.abs() < 1e-12);
    assert_eq!(calls.swap(0, Ordering::Relaxed), 4);
    assert_eq!(
        session.call_ho("num_complex_wirtinger_dz_checked", &args(3), &[cb], &ctx),
        Err(OpError::BudgetExceeded),
    );
    assert_eq!(calls.swap(0, Ordering::Relaxed), 3);
    let mut conformal = args(4).to_vec();
    conformal.insert(2, measured(Value::Float(4.0), length));
    conformal.insert(4, Value::Float(1e-12));
    assert_eq!(
        session
            .call_ho("num_conformal_numeric_checked", &conformal, &[cb], &ctx)
            .unwrap(),
        Value::Bool(true)
    );
    assert_eq!(calls.load(Ordering::Relaxed), 4);
}
