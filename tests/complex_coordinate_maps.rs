//! Checked scalar, record and optional surfaces for complex coordinate maps.

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;

const PROJECT: &str = "num_complex_coords_to_riemann_sphere";
const UNPROJECT: &str = "num_complex_coords_from_riemann_sphere";
const JOUKOWSKY: &str = "num_conformal_joukowsky";
const INVERSE: &str = "num_conformal_joukowsky_inverse";
const CAYLEY: &str = "num_conformal_cayley_transform";

const fn complex(value: Complex64) -> Value {
    Value::Complex(ComplexValue::Float {
        re: value.re,
        im: value.im,
    })
}

fn decoded(value: &Value) -> Complex64 {
    let Value::Complex(value) = value else {
        panic!("expected a complex value, got {value:?}");
    };
    let (re, im) = value.to_f64();
    Complex64::new(re, im)
}

#[test]
fn projection_round_trips_even_when_the_height_rounds_to_one() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for point in [Complex64::new(0.25, -0.5), Complex64::new(1e200, -2e200)] {
        let Ok(Value::Record(record)) = session.call_value(PROJECT, &[complex(point)], &ctx) else {
            panic!("expected a sphere coordinate record");
        };
        let [(x, a), (y, b), (z, c)] = record.fields.as_slice() else {
            panic!("expected three coordinate fields");
        };
        assert_eq!([x, y, z], ["x", "y", "z"]);
        let Ok(Value::Variant(option)) =
            session.call_value(UNPROJECT, &[a.clone(), b.clone(), c.clone()], &ctx)
        else {
            panic!("expected an optional complex result");
        };
        assert_eq!(option.tag, "Some");
        let [image] = option.payload.as_slice() else {
            panic!("expected the finite preimage");
        };
        assert!((decoded(image) - point).norm() / point.norm() <= 4e-15);
    }
    let Ok(Value::Variant(option)) = session.call_value(
        UNPROJECT,
        &[Value::Float(0.0), Value::Float(0.0), Value::Float(1.0)],
        &ctx,
    ) else {
        panic!("expected the north-pole option");
    };
    assert_eq!(option.tag, "None");
    assert!(option.payload.is_empty());

    let Ok(Value::Float(small)) = session.call_value(
        "num_complex_coords_spherical_metric_factor",
        &[complex(Complex64::new(1e151, 0.0))],
        &ctx,
    ) else {
        panic!("expected a scalar metric factor");
    };
    assert!((small / 2e-302 - 1.0).abs() <= 4e-15);

    for (component, expected_bits) in [
        (6.29973938625772e161, 1),
        (3.6856537881532267e161, 1),
        (2.8396892071351655e161, 3),
    ] {
        let point = Complex64::new(component, component);
        let Ok(Value::Float(metric)) = session.call_value(
            "num_complex_coords_spherical_metric_factor",
            &[complex(point)],
            &ctx,
        ) else {
            panic!("expected a subnormal scalar metric factor");
        };
        assert_eq!(metric.to_bits(), expected_bits);
    }
}

#[test]
fn inverse_map_retains_both_branches_and_their_order() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for point in [Complex64::new(1e200, -3e199), Complex64::new(-1e200, 3e199)] {
        let Ok(Value::Record(record)) = session.call_value(INVERSE, &[complex(point)], &ctx) else {
            panic!("expected both inverse branches");
        };
        let [(plus_name, plus), (minus_name, minus)] = record.fields.as_slice() else {
            panic!("expected plus and minus fields");
        };
        assert_eq!([plus_name, minus_name], ["plus", "minus"]);
        let (plus, minus) = (decoded(plus), decoded(minus));
        assert_eq!(plus.norm() > minus.norm(), point.re > 0.0);
        assert!((plus * minus - 1.0).norm() <= 4e-15);
        for branch in [plus, minus] {
            let Ok(image) = session.call_value(JOUKOWSKY, &[complex(branch)], &ctx) else {
                panic!("expected both branches to give finite forward values");
            };
            assert!((decoded(&image) - point).norm() / point.norm() <= 4e-15);
        }
    }
    let Ok(value) = session.call_value(CAYLEY, &[complex(Complex64::new(0.0, -2.0))], &ctx) else {
        panic!("expected a finite Möbius image outside the upper chart");
    };
    assert_eq!(decoded(&value), Complex64::new(3.0, 0.0));
}

#[test]
fn inverse_map_retains_the_subnormal_root_when_the_other_overflows() {
    use rssn::kernels::conformal::joukowsky_inverse;

    for sign in [1.0, -1.0] {
        let point = Complex64::new(sign * f64::MAX, 0.0);
        let (plus, minus) = joukowsky_inverse(point);
        let (large, small) = if sign > 0.0 {
            (plus, minus)
        } else {
            (minus, plus)
        };
        assert!(large.re.is_infinite());
        assert_eq!(large.re.signum(), sign);
        assert_eq!(small, Complex64::new(sign * (0.5 / f64::MAX), 0.0));
    }
}

#[test]
fn malformed_points_poles_and_measured_coordinates_fail() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for (op, args) in [
        (UNPROJECT, vec![Value::Float(0.0); 3]),
        (PROJECT, vec![complex(Complex64::new(f64::INFINITY, 0.0))]),
        (CAYLEY, vec![complex(Complex64::new(0.0, -1.0))]),
        (JOUKOWSKY, vec![complex(Complex64::new(0.0, 0.0))]),
    ] {
        assert!(
            matches!(session.call_value(op, &args, &ctx), Err(OpError::Domain(_))),
            "{op}"
        );
    }
    let measured = Value::Dyn(DynValue {
        value: Box::new(complex(Complex64::new(0.0, 1.0))),
        dim: DimExpr::base("Length"),
    });
    assert!(session.call_value(PROJECT, &[measured], &ctx).is_err());
}

#[test]
fn named_and_optional_shapes_are_declared_dimensionless() {
    let manifest = Session::new().manifest();
    let shapes = [
        (
            PROJECT,
            "complex_coords::to_riemann_sphere",
            vec![TypeExpr::Complex],
            TypeExpr::Record(
                ["x", "y", "z"]
                    .map(|key| (key.into(), TypeExpr::Float))
                    .to_vec(),
            ),
        ),
        (
            UNPROJECT,
            "complex_coords::from_riemann_sphere",
            vec![TypeExpr::Float; 3],
            TypeExpr::Variant(vec![
                ("None".into(), vec![]),
                ("Some".into(), vec![TypeExpr::Complex]),
            ]),
        ),
        (
            INVERSE,
            "conformal::joukowsky_inverse",
            vec![TypeExpr::Complex],
            TypeExpr::Record(
                ["plus", "minus"]
                    .map(|key| (key.into(), TypeExpr::Complex))
                    .to_vec(),
            ),
        ),
        (
            CAYLEY,
            "conformal::cayley_transform",
            vec![TypeExpr::Complex],
            TypeExpr::Complex,
        ),
    ];
    for (name, source, params, result) in shapes {
        let Some(op) = manifest.op(name) else {
            panic!("missing {name}");
        };
        let path = format!("crate::kernels::{source}");
        assert_eq!(op.rust_path.as_deref(), Some(path.as_str()));
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        let [signature] = op.overloads.as_slice() else {
            panic!("expected one declared signature for {name}");
        };
        assert_eq!(signature.params, params);
        assert_eq!(signature.ret, result);
    }
}

#[test]
fn polar_coordinates_keep_the_measured_radius_and_dimensionless_angle() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let order = DimExpr::Base("Length".into());
    let value = Value::Dyn(DynValue {
        value: Box::new(complex(Complex64::new(3.0, 4.0))),
        dim: order.clone(),
    });
    let Value::Record(result) = session
        .call_value(
            "num_complex_coords_to_polar",
            std::slice::from_ref(&value),
            &ctx,
        )
        .unwrap()
    else {
        panic!("the polar calculation must return a record");
    };
    let [
        (_, angle @ Value::Float(radians)),
        (_, radius @ Value::Dyn(radius_value)),
    ] = result.fields.as_slice()
    else {
        panic!("expected the radius witness and dimensionless angle: {result:?}");
    };
    assert_eq!(radius_value.dim, order);
    assert_eq!(*radius_value.value, Value::Float(5.0));
    assert!((*radians - 4.0_f64.atan2(3.0)).abs() <= f64::EPSILON);
    let Value::Dyn(restored) = session
        .call_value(
            "num_complex_coords_from_polar",
            &[radius.clone(), angle.clone()],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected the original physical dimension");
    };
    assert_eq!(restored.dim, order);
    assert!((decoded(&restored.value) - Complex64::new(3.0, 4.0)).norm() < 2e-15);
    let Value::Dyn(squared) = session
        .call_value("num_complex_coords_abs_sq", &[value], &ctx)
        .unwrap()
    else {
        panic!("expected squared physical units");
    };
    assert_eq!(*squared.value, Value::Float(25.0));
    assert_eq!(
        squared.dim,
        order.pow(rssn::kernel::signature::Exp::Rat(
            rssn::kernel::signature::Rational::new(2, 1)
        ))
    );
    assert!(
        session
            .call_value(
                "num_complex_coords_from_polar",
                &[Value::Float(-1.0), Value::Float(0.0)],
                &ctx
            )
            .is_err()
    );
}

#[test]
fn unit_direction_and_root_radius_survive_both_exponent_extremes() {
    use rssn::kernels::complex_coords::analytic_continuation_root;
    use rssn::kernels::complex_coords::sign_complex;
    use rssn::kernels::complex_coords::to_polar;
    use std::f64::consts::FRAC_1_SQRT_2;
    use std::f64::consts::PI;
    let smallest = f64::from_bits(1);
    for scale in [smallest, f64::MAX] {
        let point = Complex64::new(scale, -scale);
        let expected = Complex64::new(FRAC_1_SQRT_2, -FRAC_1_SQRT_2);
        assert!((sign_complex(point) - expected).norm() <= 2e-15);
        let root = analytic_continuation_root(&[point], 2, 0)[0];
        let root_radius = scale.sqrt() * 2.0_f64.powf(0.25);
        assert!((root.norm() / root_radius - 1.0).abs() <= 2e-15);
    }
    assert_eq!(
        to_polar(Complex64::new(-1.0, -0.0)).1.to_bits(),
        (-PI).to_bits()
    );
    assert_eq!(
        to_polar(Complex64::new(-1.0, 0.0)).1.to_bits(),
        PI.to_bits()
    );
}

#[test]
fn phase_reconstruction_tracks_many_turns_in_both_directions() {
    use rssn::kernels::complex_coords::analytic_continuation_root;
    use rssn::kernels::complex_coords::unwrap_phase;
    use std::f64::consts::FRAC_PI_2;
    use std::f64::consts::PI;
    let corners = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    let turns = 100_000;
    for direction in [-1.0_f64, 1.0] {
        let points = (0..=4 * turns)
            .map(|step| {
                let point = corners[step % 4];
                if direction.is_sign_negative() {
                    point.conj()
                } else {
                    point
                }
            })
            .collect::<Vec<_>>();
        let phase = unwrap_phase(&points);
        for (step, radians) in phase.into_iter().enumerate().step_by(197) {
            let expected = direction * f64::from(u32::try_from(step).unwrap()) * FRAC_PI_2;
            assert!((radians - expected).abs() <= 2.0 * f64::EPSILON * expected.abs().max(1.0));
        }
    }
    let roots = analytic_continuation_root(
        &[corners[0], corners[1], corners[2], corners[3], corners[0]],
        2,
        -2,
    );
    assert!((roots[4] - Complex64::from_polar(1.0, PI)).norm() <= 1e-15);
}
