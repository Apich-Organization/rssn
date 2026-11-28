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
