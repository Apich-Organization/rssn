//! Numerical root bounds and unit-disk value API regressions.

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::complex_roots;
use rssn::kernels::univalent;

fn tensor(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data,
        elem,
    })
}

#[test]
fn sharp_cauchy_radius_keeps_nonzero_roots_when_the_constant_term_is_zero() {
    let coefficient = |value| Complex64::new(value, 0.0);
    let base = [0.0, -4.0, 1.0, -0.0];
    for scale in [1.0, 2.0_f64.powi(-900), 2.0_f64.powi(900)] {
        let polynomial = base.map(|a| coefficient(a * scale));
        let radius = complex_roots::cauchy_sharp_root_bound(&polynomial).unwrap();
        assert!((radius - 4.0).abs() <= 4.0 * f64::EPSILON);
    }
    let monomial = [coefficient(0.0), coefficient(0.0), coefficient(8.0)];
    assert_eq!(complex_roots::cauchy_sharp_root_bound(&monomial), Some(0.0));
    let nonfinite = [coefficient(1.0), coefficient(2.0), coefficient(f64::NAN)];
    assert_eq!(complex_roots::cauchy_sharp_root_bound(&nonfinite), None);
}

#[test]
fn root_bounds_use_finite_ratios_when_coefficient_norms_do_not_fit() {
    let extreme = Complex64::new(f64::MAX, f64::MAX);
    assert_eq!(
        complex_roots::cauchy_root_bound(&[extreme, extreme]),
        Some(2.0)
    );
    assert!(
        (complex_roots::cauchy_sharp_root_bound(&[extreme, extreme]).unwrap() - 1.0).abs() < 1e-14
    );
    assert_eq!(
        complex_roots::fujiwara_root_bound(&[extreme, extreme]),
        Some(1.0)
    );

    let small = f64::from_bits(1);
    let subnormal = [Complex64::new(small, small), Complex64::new(small, 0.0)];
    assert!(
        (complex_roots::cauchy_sharp_root_bound(&subnormal).unwrap() - std::f64::consts::SQRT_2)
            .abs()
            < 2e-14
    );

    let mut polynomial = [Complex64::ZERO; 5];
    polynomial[0].re = 2.0_f64.powi(1_022);
    polynomial[4].re = 2.0_f64.powi(-1_022);
    let expected = 2.0_f64.powi(511);
    let sharp = complex_roots::cauchy_sharp_root_bound(&polynomial).unwrap();
    assert!((sharp / expected - 1.0).abs() < 2e-13);
    let fujiwara = complex_roots::fujiwara_root_bound(&polynomial).unwrap();
    assert!((1.68..1.69).contains(&(fujiwara / expected)));
    assert_eq!(
        complex_roots::cauchy_root_bound(&polynomial),
        Some(f64::INFINITY)
    );
}

#[test]
fn annulus_handles_underflow_and_rotation_stays_finite_near_the_disk_boundary() {
    let small = f64::from_bits(1);
    assert_eq!(
        complex_roots::enestrom_kakeya_bounds(&[small, f64::MAX]),
        Some((0.0, small))
    );
    assert_eq!(
        complex_roots::enestrom_kakeya_bounds(&[f64::MAX, small]),
        Some((f64::MAX, f64::INFINITY))
    );
    let switch = std::f64::consts::FRAC_1_SQRT_2;
    for radius in [switch.next_down(), switch, switch.next_up()] {
        let bound = univalent::koebe_rotation_bound(radius).unwrap();
        assert!((bound - std::f64::consts::PI).abs() < 5e-15);
    }
    // High-precision evaluation at the exact binary64 radius 1 - 2^-53.
    let edge = univalent::koebe_rotation_bound(1.0_f64.next_down()).unwrap();
    assert!((edge - 39.18524604270695).abs() < 1e-14);
    assert_eq!(univalent::koebe_rotation_bound(1.0), None);
}

#[test]
fn value_bounds_preserve_the_dimensionless_coordinate_and_report_input_errors() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let dim = DimExpr::base("Mass");
    let annulus = "num_enestrom_kakeya_bounds";
    let coefficients = [1.0, 2.0, 4.0].into_iter().map(Value::Float).collect();
    let Value::Record(bounds) = session
        .call_value(
            annulus,
            &[tensor(coefficients, TypeExpr::Float.with_dim(dim.clone()))],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected a named pair of root bounds");
    };
    assert_eq!(
        bounds.fields,
        [
            ("lower".to_owned(), Value::Float(0.5)),
            ("upper".to_owned(), Value::Float(0.5))
        ]
    );

    let complex = |re| Value::Complex(ComplexValue::Float { re, im: 0.0 });
    let coefficients = tensor(
        vec![complex(0.0), complex(-4.0), complex(1.0)],
        TypeExpr::Complex.with_dim(dim.clone()),
    );
    let observed = session
        .call_value("num_cauchy_sharp_root_bound", &[coefficients], &ctx)
        .unwrap();
    let Value::Float(radius) = observed else {
        panic!("expected a dimensionless radius, got {observed:?}");
    };
    assert!((radius - 4.0).abs() <= 4.0 * f64::EPSILON);

    let mixed = tensor(
        vec![
            Value::Float(1.0),
            Value::Dyn(DynValue {
                value: Box::new(Value::Float(2.0)),
                dim: DimExpr::base("Time"),
            }),
        ],
        TypeExpr::Float.with_dim(dim.clone()),
    );
    assert!(session.call_value(annulus, &[mixed], &ctx).is_err());

    let dimensioned_radius = Value::Dyn(DynValue {
        value: Box::new(Value::Float(0.5)),
        dim,
    });
    assert!(
        session
            .call_value("num_koebe_growth_bounds", &[dimensioned_radius], &ctx)
            .is_err()
    );
    assert!(matches!(
        session.call_value("num_koebe_distortion_bounds", &[Value::Float(1.0)], &ctx),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        session.call_value(
            annulus,
            &[tensor(
                vec![Value::Float(1.0), Value::Float(-2.0)],
                TypeExpr::Float
            )],
            &ctx
        ),
        Err(OpError::Domain(_))
    ));
}
