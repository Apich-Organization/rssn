//! Typed value contracts for dimensionless Weierstrass kernels.

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

const PRIME: &str = "num_weierstrass_p_prime_complex";
const INVARIANTS: &str = "num_weierstrass_invariants_complex";

const fn complex(value: Complex64) -> Value {
    Value::Complex(ComplexValue::Float {
        re: value.re,
        im: value.im,
    })
}

fn decoded(value: &Value) -> Complex64 {
    let Value::Complex(value) = value else {
        panic!("expected a complex result, got {value:?}");
    };
    let (re, im) = value.to_f64();
    Complex64::new(re, im)
}

fn near(
    got: Complex64,
    expected: Complex64,
    tolerance: f64,
) {
    assert!(
        (got - expected).norm() <= tolerance * (1.0 + expected.norm()),
        "got {got}, expected {expected}"
    );
}

#[test]
fn oblique_lattice_matches_independent_references_and_period_translations() {
    let session = Session::new();
    let context = CallCtx::new();
    let omega1 = Complex64::new(1.0, 0.0);
    let omega2 = Complex64::new(0.2, 1.1);
    let z = Complex64::new(0.3, 0.4);

    // Independently computed cotangent row sums of the defining lattice.
    let expected = Complex64::new(15.117_919_150_272_616, 6.084_359_385_787_47);
    let g2 = Complex64::new(8.703_070_400_834_047, 1.856_095_769_282_975_1);
    let g3 = Complex64::new(3.820_148_695_061_241_6, -2.168_333_627_017_415_4);
    for (first, second) in [(omega1, omega2), (omega1, -omega2), (omega2, omega1)] {
        let prime = |point| {
            decoded(
                &session
                    .call_value(
                        PRIME,
                        &[complex(point), complex(first), complex(second)],
                        &context,
                    )
                    .unwrap(),
            )
        };
        near(prime(z), expected, 2e-13);
        near(prime(-z), -expected, 2e-13);
        near(prime(z + 2.0 * first + 2.0 * second), expected, 2e-13);
        near(prime(z + 16.0 * second), expected, 2e-13);
        near(prime(first), Complex64::new(0.0, 0.0), 2e-13);

        let Value::Record(record) = session
            .call_value(INVARIANTS, &[complex(first), complex(second)], &context)
            .unwrap()
        else {
            panic!("expected a named complex record");
        };
        assert_eq!(
            record
                .fields
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["g2", "g3"]
        );
        near(decoded(&record.fields[0].1), g2, 2e-13);
        near(decoded(&record.fields[1].1), g3, 2e-13);
    }
}

#[test]
fn lattice_poles_dimensions_and_invalid_period_pairs_raise_errors() {
    let session = Session::new();
    let context = CallCtx::new();
    let one = Complex64::new(1.0, 0.0);
    let zero = Complex64::new(0.0, 0.0);
    let i = Complex64::new(0.0, 1.0);
    let value = Complex64::new(0.3, 0.4);

    for pole in [zero, 2.0 * one, 2.0 * i, 2.0 * (one + i)] {
        let error = session.call_value(PRIME, &[complex(pole), complex(one), complex(i)], &context);
        assert!(matches!(error, Err(OpError::Domain(message)) if message.contains("pole")));
    }

    let epsilon = 2.0_f64.powi(-20);
    let near_pole = session
        .call_value(
            PRIME,
            &[complex(2.0 * one + epsilon), complex(one), complex(i)],
            &context,
        )
        .unwrap();
    near(
        decoded(&near_pole),
        Complex64::new(-2.0 / epsilon.powi(3), 0.0),
        2e-13,
    );

    for periods in [
        (zero, i),
        (one, zero),
        (one, 2.0 * one),
        (one, Complex64::new(f64::INFINITY, 1.0)),
    ] {
        let w1 = complex(periods.0);
        let w2 = complex(periods.1);
        assert!(matches!(
            session.call_value(INVARIANTS, &[w1.clone(), w2.clone()], &context),
            Err(OpError::Domain(_))
        ));
        assert!(matches!(
            session.call_value(PRIME, &[complex(value), w1, w2], &context),
            Err(OpError::Domain(_))
        ));
    }

    let measured = Value::Dyn(DynValue {
        value: Box::new(complex(value)),
        dim: DimExpr::base("Length"),
    });
    assert!(
        session
            .call_value(PRIME, &[measured, complex(one), complex(i)], &context)
            .is_err()
    );
    let nonfinite = complex(Complex64::new(0.0, f64::NAN));
    assert!(
        session
            .call_value(PRIME, &[nonfinite, complex(one), complex(i)], &context)
            .is_err()
    );
    // Even valid finite inputs can put the native derivative outside binary64.
    let tiny = complex(Complex64::new(1e-200, 0.0));
    assert!(matches!(
        session.call_value(PRIME, &[tiny, complex(one), complex(i)], &context),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn weierstrass_kernels_advertise_strongly_typed_dimensionless_contracts() {
    let manifest = Session::new().manifest();
    for (name, method, params, result) in [
        (
            PRIME,
            "weierstrass_p_prime_complex",
            vec![TypeExpr::Complex; 3],
            TypeExpr::Complex,
        ),
        (
            INVARIANTS,
            "weierstrass_invariants_complex",
            vec![TypeExpr::Complex; 2],
            TypeExpr::Record(vec![
                ("g2".into(), TypeExpr::Complex),
                ("g3".into(), TypeExpr::Complex),
            ]),
        ),
    ] {
        let op = manifest.op(name).unwrap();
        let path = format!("crate::kernels::special_complex_ext::elliptic_modular::{method}");
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(op.rust_path.as_deref(), Some(path.as_str()));
        let [signature] = op.overloads.as_slice() else {
            panic!("expected one Weierstrass overload for {name}");
        };
        assert_eq!(signature.params, params);
        assert_eq!(signature.ret, result);
    }
}
