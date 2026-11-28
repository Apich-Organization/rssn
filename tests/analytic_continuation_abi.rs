//! Typed complex continuation through the same API used by language bindings.
use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::analytic_continuation;

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}

fn vector(entries: &[(f64, f64)]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![entries.len()],
        elem: TypeExpr::Complex,
        data: entries.iter().map(|&(re, im)| complex(re, im)).collect(),
    })
}

#[test]
fn pade_round_trips_its_named_record_and_rejects_a_pole() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let pade = session
        .call_value(
            "num_continuation_pade",
            &[
                vector(&[(1.0, 0.0), (1.0, 0.0)]),
                Value::int(0),
                Value::int(1),
            ],
            &ctx,
        )
        .expect("Padé coefficients of 1 / (1 - z)");
    assert_eq!(
        session.call_value(
            "num_continuation_eval_pade",
            &[pade.clone(), complex(0.5, 0.0)],
            &ctx
        ),
        Ok(complex(2.0, 0.0))
    );
    assert!(matches!(
        session.call_value(
            "num_continuation_eval_pade",
            &[pade, complex(1.0, 0.0)],
            &ctx
        ),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn extreme_pade_ratios_use_scaled_complex_division() {
    for scale in [1e-300, 1e300] {
        let coefficient = Complex64::new(scale, scale);
        assert_eq!(
            analytic_continuation::eval_pade_complex(
                &[coefficient],
                &[coefficient],
                Complex64::ZERO
            ),
            Complex64::ONE
        );
    }
    let coefficients = [Complex64::ONE];
    assert!(
        analytic_continuation::pade_approximant_complex(&coefficients, usize::MAX, 1).is_none()
    );
    assert!(
        analytic_continuation::borel_pade_resum(&coefficients, Complex64::ZERO, usize::MAX, 1)
            .is_none()
    );
}

#[test]
fn taylor_recentering_and_wynn_agree_with_exact_values() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let polynomial = vector(&[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]);
    let translated = session.call_value(
        "num_continuation_recenter",
        &[polynomial, vector(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])],
        &ctx,
    );
    assert_eq!(
        translated,
        Ok(vector(&[(17.0, 0.0), (14.0, 0.0), (3.0, 0.0)]))
    );
    let accelerated = session.call_value(
        "num_continuation_wynn",
        &[vector(&[(1.0, 0.0), (1.5, 0.0), (1.75, 0.0)])],
        &ctx,
    );
    assert_eq!(accelerated, Ok(complex(2.0, 0.0)));
}

#[test]
fn declared_continuation_paths_have_explicit_dimension_contracts() {
    let session = Session::new();
    let manifest = session.manifest();
    for (operator, function) in [
        ("pade", "pade_approximant_complex"),
        ("eval_pade", "eval_pade_complex"),
        ("wynn", "wynn_epsilon_complex"),
        ("recenter", "weierstrass_chain_continue"),
    ] {
        let name = format!("num_continuation_{operator}");
        let path = format!("crate::kernels::analytic_continuation::{function}");
        let entry = manifest
            .ops
            .iter()
            .find(|op| op.name == name)
            .expect("declared operator");
        assert_eq!(entry.rust_path.as_deref(), Some(path.as_str()));
        assert_eq!(entry.keys.len(), 1);
        assert_eq!(entry.signature_source, SignatureSource::Declared);
        assert!(entry.dim_rule.is_some());
    }
}
