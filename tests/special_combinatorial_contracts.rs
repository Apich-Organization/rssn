//! Boundary and reference cases for the native numeric and value APIs.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::special;

fn integers(values: &[u64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

#[test]
fn logarithmic_coefficients_keep_the_small_count_when_n_exceeds_binary64_precision() {
    let n = u64::MAX;
    assert_eq!(special::ln_binomial(n, 0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(special::ln_binomial(n, 1), (n as f64).ln());
    assert_eq!(special::ln_binomial(n, n - 1), special::ln_binomial(n, 1));
    assert_eq!(special::ln_binomial(3, 4), f64::NEG_INFINITY);
    // References are log values of exact integer binomial coefficients. The
    // final two cases straddle the product/stirling branch threshold.
    for (n, k, expected) in [
        (1_000, 500, 689.467_261_567_851_2),
        ((1_u64 << 60) + 23, 70, 2_680.779_114_785_993),
        (u64::MAX, 32, 1_338.007_466_330_653),
        (u64::MAX, 33, 1_378.872_378_325_023),
    ] {
        let observed = special::ln_binomial(n, k);
        assert!(
            (observed - expected).abs() <= 4e-15 * expected,
            "{n}, {k}: {observed}"
        );
        assert_eq!(observed, special::ln_binomial(n, n - k));
    }
    assert!(special::binomial(1_000, 500).is_finite());
    assert_eq!(special::binomial(2_048, 1_024), f64::INFINITY);
}

#[test]
fn harmonic_matches_exact_small_values_and_large_references() {
    assert_eq!(special::harmonic(0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(special::harmonic(1), 1.0);
    // Exact fractions include H_5 and the values on both sides of the table boundary.
    for (n, expected) in [
        (5, 137.0 / 60.0),
        (15, 1_195_757.0 / 360_360.0),
        (16, 2_436_559.0 / 720_720.0),
        (17, 42_142_223.0 / 12_252_240.0),
        (1_024, 7.509_175_672_278_133_5_f64),
        ((1_u64 << 53) + 1, 37.314_016_234_578_64_f64),
        (u64::MAX, 44.938_635_220_738_036_f64),
    ] {
        let observed = special::harmonic(n);
        assert!(
            observed.to_bits().abs_diff(expected.to_bits()) <= 2,
            "H_{n}: {observed:.17e}, expected {expected:.17e}"
        );
    }
}

#[test]
fn checked_multinomial_keeps_empty_partitions_and_validates_the_entire_sum() {
    let n = u64::MAX;
    assert_eq!(special::checked_multinomial(0, &[]), Some(1.0));
    assert_eq!(special::multinomial(9, &[2, 3, 0, 4]), 1_260.0);
    assert_eq!(special::checked_multinomial(n, &[n - 1, 1]), Some(n as f64));
    assert_eq!(special::checked_multinomial(n, &[n, 1]), None);
    assert_eq!(special::checked_multinomial(2_050, &[1_024, 1_024]), None);
    assert_eq!(
        special::checked_multinomial(2_048, &[1_024, 1_024]),
        Some(f64::INFINITY)
    );
}

#[test]
fn inverse_beta_preserves_analytic_tails_and_brackets_polynomial_quantiles() {
    let tiny = 2.0_f64.powi(-54);
    let x = special::inverse_beta_reg(1.0, 2.0, tiny).unwrap();
    assert!((x / (0.5 * tiny) - 1.0).abs() <= 2.0 * f64::EPSILON);
    assert!((special::inverse_beta_reg(2.0, 1.0, 0.25).unwrap() - 0.5).abs() <= f64::EPSILON);
    // I_x(2,2) = x^2(3 - 2x).
    for (p, exact) in [(5.0 / 32.0, 0.25), (27.0 / 32.0, 0.75)] {
        assert!((special::inverse_beta_reg(2.0, 2.0, p).unwrap() - exact).abs() <= 2e-14);
    }
    assert_eq!(special::inverse_beta_reg(f64::MAX, f64::MAX, 0.5), Ok(0.5));
    assert!(special::inverse_beta_reg(f64::MAX, f64::MAX, 0.4).is_err());
    assert!(special::inverse_beta_reg(f64::INFINITY, 1.0, 0.5).is_err());
    assert!(special::inverse_beta_reg(2.0, 3.0, f64::NAN).is_err());
}

#[test]
fn scalar_adapters_keep_finite_results_and_report_domain_errors() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for (name, inputs, expected) in [
        ("num_harmonic", vec![Value::int(0)], 0.0),
        ("num_harmonic", vec![Value::int(5)], 137.0 / 60.0),
        (
            "num_binomial",
            vec![Value::int(u64::MAX), Value::int(1)],
            u64::MAX as f64,
        ),
        (
            "num_ln_binomial",
            vec![Value::int(u64::MAX), Value::int(1)],
            (u64::MAX as f64).ln(),
        ),
        (
            "num_multinomial",
            vec![Value::int(9), integers(&[2, 3, 4])],
            1_260.0,
        ),
        (
            "num_inverse_beta_reg",
            vec![Value::Float(1.0), Value::Float(1.0), Value::Float(0.75)],
            0.75,
        ),
    ] {
        assert_eq!(
            session.call_value(name, &inputs, &ctx).unwrap(),
            Value::Float(expected),
            "{name}"
        );
    }
    for (name, inputs) in [
        ("num_ln_binomial", vec![Value::int(7), Value::int(8)]),
        (
            "num_multinomial",
            vec![Value::int(u64::MAX), integers(&[u64::MAX, 1])],
        ),
        (
            "num_multinomial",
            vec![Value::int(2_048), integers(&[1_024, 1_024])],
        ),
        (
            "num_inverse_beta_reg",
            vec![Value::Float(2.0), Value::Float(3.0), Value::Float(1.1)],
        ),
    ] {
        assert!(
            matches!(
                session.call_value(name, &inputs, &ctx),
                Err(OpError::Domain(_))
            ),
            "{name}"
        );
    }
    let measured = Value::Dyn(DynValue {
        value: Box::new(Value::int(3)),
        dim: DimExpr::base("Length"),
    });
    assert!(
        session
            .call_value("num_harmonic", &[measured], &ctx)
            .is_err()
    );
}
