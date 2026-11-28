//! Exact rational references for integer powers of forward dual numbers.

use num_rational::BigRational;
use num_traits::ToPrimitive;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernels::autodiff::Dual;

fn reference(
    value: f64,
    tangent: f64,
    n: i32,
) -> [f64; 3] {
    let x = BigRational::from_float(value).unwrap();
    let d = BigRational::from_float(tangent).unwrap();
    let f = BigRational::from_integer(n.into());
    let g = BigRational::from_integer((i64::from(n) - 1).into());
    [
        x.pow(n),
        &f * &d * x.pow(n - 1),
        f * g * &d * &d * x.pow(n - 2),
    ]
    .map(|v| v.to_f64().unwrap())
}

fn assert_float(
    actual: f64,
    expected: f64,
) {
    if expected == 0.0 || !expected.is_finite() {
        assert_eq!(actual, expected);
    } else if expected.is_subnormal() {
        assert!(
            actual.is_sign_negative() == expected.is_sign_negative()
                && actual.to_bits().abs_diff(expected.to_bits()) <= 2,
            "{actual:e} != {expected:e}",
        );
    } else {
        assert!(
            actual.is_finite() && (actual / expected - 1.0).abs() <= 4e-14,
            "{actual:e} != {expected:e}",
        );
    }
}

fn assert_power(
    value: f64,
    tangent: f64,
    n: i32,
) {
    let x = Dual::new(value, tangent);
    let [expected, slope, curvature] = reference(value, tangent, n);
    let first = x.powi(n);
    assert_float(first.val, expected);
    assert_float(first.deriv, slope);
    let (second, actual_curvature) = x.powi_with_second_derivative(n);
    assert_float(second.val, expected);
    assert_float(second.deriv, slope);
    assert_float(actual_curvature, curvature);
}

#[test]
fn native_integer_power_exposes_primal_slope_and_curvature() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let dual = session
        .call_value(
            "num_dual_new",
            &[Value::Float(2.0), Value::Float(3.0)],
            &ctx,
        )
        .unwrap();
    let result = session
        .call_value(
            "num_dual_powi_with_second_derivative",
            &[dual.clone(), Value::int(4)],
            &ctx,
        )
        .unwrap();
    let Value::Record(record) = result else {
        panic!("expected a flat derivative record");
    };
    assert_eq!(
        record.fields,
        [
            ("val".into(), Value::Float(16.0)),
            ("deriv".into(), Value::Float(96.0)),
            ("second".into(), Value::Float(432.0)),
        ],
    );
    assert!(
        session
            .call_value(
                "num_dual_powi_with_second_derivative",
                &[dual, Value::int(i64::from(i32::MAX) + 1)],
                &ctx,
            )
            .is_err(),
    );
}

#[test]
fn powers_recover_finite_derivatives_across_intermediate_scale_limits() {
    for (value, tangent, n) in [
        (1e-200, 1e200, 3),
        (1e200, 1e-200, 3),
        (1e-200, 1e200, 4),
        (1e200, 1e-200, 4),
        (-1e200, -1e-200, 3),
        (-1e-200, -1e200, 4),
        (1e160, 1e200, -2),
        (-1e-160, 1e-200, -2),
        (0.5, f64::MAX, 2),
        (2.0_f64.powi(-537), 1.0, 2),
        (2.0_f64.powi(-538), 1.0, 2),
        (f64::MIN_POSITIVE, f64::from_bits(1), -1),
        (1e200, -0.0, 3),
    ] {
        assert_power(value, tangent, n);
    }
}

#[test]
fn finite_integer_power_and_both_derivatives_match_exact_rationals() {
    for value in [-1.125, -0.75, 0.16, 0.99, 1.5, 2.0] {
        for n in [-70, -11, -3, -1, 2, 3, 4, 8, 17, 70] {
            assert_power(value, -0.625, n);
        }
    }
}

#[test]
fn minimum_integer_exponent_keeps_parity_for_both_derivatives() {
    for base in [-1.0, 1.0] {
        let x = Dual::new(base, 0.5);
        for exponent in [i32::MIN, i32::MIN + 1, i32::MAX] {
            let n = f64::from(exponent);
            let sign = if exponent % 2 == 0 {
                1.0
            } else {
                base
            };
            let derivative = n * 0.5 * base * sign;
            let curvature = n * (n - 1.0) * 0.25 * sign;
            assert_eq!(x.powi(exponent), Dual::new(sign, derivative));
            assert_eq!(
                x.powi_with_second_derivative(exponent),
                (Dual::new(sign, derivative), curvature),
            );
        }
    }
}

#[test]
fn constant_linear_and_ieee_boundary_cases_keep_their_conventions() {
    let x = Dual::new(-0.0, 2.0);
    assert_eq!(x.powi_with_second_derivative(0), (Dual::constant(1.0), 0.0));
    assert_eq!(x.powi_with_second_derivative(1), (x, 0.0));
    assert_eq!(x.powi_with_second_derivative(2), (Dual::constant(0.0), 8.0));
    assert!(x.powi(3).val.is_sign_negative());
    assert!(x.powi(-1).val == f64::NEG_INFINITY);
    assert!(Dual::constant(0.0).powi(-1).deriv.is_nan());
    assert!(Dual::var(f64::INFINITY).powi(-1).val == 0.0);
    assert!(Dual::var(f64::NAN).powi(5).val.is_nan());
    assert_eq!(Dual::var(f64::NAN).powi(0), Dual::constant(1.0));
    assert!(Dual::constant(1e200).powi(3).deriv == 0.0);
}
