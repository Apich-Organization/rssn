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
use rssn::kernel::value::TensorValue;
use rssn::kernels::analytic_continuation as a;
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

fn coordinate(value: &Value) -> C {
    match value {
        | Value::Complex(z) => {
            let (re, im) = z.to_f64();
            C::new(re, im)
        },
        | Value::Dyn(z) => coordinate(&z.value),
        | _ => panic!("expected a complex scalar"),
    }
}

fn tensor(
    shape: Vec<usize>,
    entries: &[C],
) -> Value {
    Value::Tensor(TensorValue::of(
        shape,
        entries.iter().copied().map(complex).collect(),
        TypeExpr::Complex,
    ))
}

#[test]
fn holonomic_long_step_adapts_and_refuses_invalid_state_or_singularity() {
    let exponential = a::HolonomicOde {
        coeffs: vec![vec![-C::ONE], vec![C::ONE]],
    };
    let path = [C::ZERO, C::new(30.0, 0.0)];
    let result = a::continue_holonomic_along_path(&exponential, &path, &[C::ONE], 24).unwrap();
    let expected = 30.0_f64.exp();
    assert!((result[0].re / expected - 1.0).abs() < 5e-13);
    assert!(a::continue_holonomic_along_path(&exponential, &path, &[C::ONE, C::ONE], 24).is_none());
    let logarithmic = a::HolonomicOde {
        coeffs: vec![vec![], vec![C::ONE], vec![C::ZERO, C::ONE]],
    };
    let at_singularity = [C::ONE, C::ZERO];
    assert!(
        a::continue_holonomic_along_path(&logarithmic, &at_singularity, &[C::ZERO, C::ONE], 24)
            .is_none()
    );
}

#[test]
fn a_delayed_nonzero_taylor_term_cannot_make_a_step_look_constant() {
    let mut factor = vec![C::ZERO; 13];
    factor[12] = -C::ONE;
    let ode = a::HolonomicOde {
        coeffs: vec![factor, vec![C::ONE]],
    };
    let end = 1.2_f64;
    let path = [C::ZERO, C::new(end, 0.0)];
    let got = a::continue_holonomic_along_path(&ode, &path, &[C::ONE], 4).unwrap();
    let reference = (end.powi(13) / 13.0).exp();
    assert!((got[0].re / reference - 1.0).abs() < 2e-12);
}

#[test]
fn holonomic_monodromy_transports_basis_columns_for_logarithm() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let coefficients = tensor(
        vec![3, 2],
        &[C::ZERO, C::ZERO, C::ONE, C::ZERO, C::ZERO, C::ONE],
    );
    assert_eq!(
        session
            .call_value(
                "num_continuation_holonomic_order",
                std::slice::from_ref(&coefficients),
                &ctx
            )
            .unwrap(),
        Value::int(2)
    );
    let path = tensor(vec![4], &[C::ONE, C::I, -C::ONE, -C::I]);
    let result = session
        .call_value(
            "num_continuation_holonomic_monodromy",
            &[coefficients, path, Value::int(32)],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(matrix) = result else {
        panic!("expected monodromy matrix")
    };
    assert_eq!(matrix.shape, [2, 2]);
    for (actual, expected) in
        matrix
            .data
            .iter()
            .zip([C::ONE, C::I * std::f64::consts::TAU, C::ZERO, C::ONE])
    {
        assert!((coordinate(actual) - expected).norm() < 1e-11);
    }
}

#[test]
fn upper_borel_lateral_branch_is_selected_only_for_axial_poles() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for (sign, expected) in [
        (-1.0, C::new(0.7226572337764452, 0.0)),
        (1.0, C::new(1.3409654195801466, 0.8503366631752727)),
    ] {
        let mut factorial = 1.0;
        let sequence: Vec<_> = (0..=8)
            .map(|k| {
                if k > 0 {
                    factorial *= sign * f64::from(k);
                }
                C::new(factorial, 0.0)
            })
            .collect();
        let result = session
            .call_value(
                "num_continuation_borel_pade",
                &[
                    tensor(vec![sequence.len()], &sequence),
                    complex(C::new(0.5, 0.0)),
                    Value::int(0),
                    Value::int(1),
                ],
                &ctx,
            )
            .unwrap();
        assert!((coordinate(&result) - expected).norm() < 2e-11);
    }
}

#[test]
fn laurent_fft_retains_amplitude_units_and_the_actual_budget() {
    let session = Session::new();
    let axis = DimExpr::base("Length");
    let amplitude = DimExpr::base("Time");
    let output = amplitude.clone();
    let mut ctx = CallCtx::new();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                let z = coordinate(&args[0]) / 2.0;
                Ok(measured(complex(3.0 * (C::ONE / z + z * z)), &output))
            },
            Signature::new(
                vec![TypeExpr::Complex.with_dim(axis.clone())],
                TypeExpr::Complex.with_dim(amplitude),
            ),
        )))
        .unwrap();
    let args = |budget| {
        [
            measured(complex(C::ZERO), &axis),
            measured(Value::Float(2.0), &axis),
            Value::int(1),
            Value::int(2),
            Value::int(64),
            Value::int(budget),
        ]
    };
    assert!(matches!(
        session.call_ho("num_continuation_laurent_fft", &args(63), &[callback], &ctx),
        Err(OpError::BudgetExceeded)
    ));
    let result = session
        .call_ho("num_continuation_laurent_fft", &args(64), &[callback], &ctx)
        .unwrap();
    let Value::Record(report) = result else {
        panic!("expected Laurent report")
    };
    assert_eq!(report.get("first_power"), Some(&Value::int(-1)));
    assert_eq!(report.get("evaluations"), Some(&Value::int(64)));
    let Some(Value::Tensor(values)) = report.get("coefficients") else {
        panic!("expected Fourier coefficients")
    };
    for (value, expected) in values.data.iter().zip([3.0, 0.0, 0.0, 3.0]) {
        assert!((coordinate(value) - expected).norm() < 1e-12);
    }
}

#[test]
fn eighth_order_monodromy_is_identity_around_an_ordinary_point() {
    let order = 8;
    let mut polynomials = vec![vec![C::ZERO]; order + 1];
    polynomials[0] = vec![C::new(0.5, 0.0), C::new(-0.3, 0.01)];
    polynomials[order] = vec![C::ONE, C::new(-0.1, 0.0)];
    let ode = a::HolonomicOde { coeffs: polynomials };
    let loop_points: Vec<_> = (0..32)
        .map(|k| C::from_polar(0.3, std::f64::consts::TAU * f64::from(k) / 32.0))
        .collect();
    let monodromy = a::holonomic_monodromy_matrix(&ode, &loop_points, 24).unwrap();
    for (row, values) in monodromy.iter().enumerate() {
        for (column, &actual) in values.iter().enumerate() {
            let expected = if row == column {
                C::ONE
            } else {
                C::ZERO
            };
            assert!((actual - expected).norm() < 4e-12);
        }
    }
}
