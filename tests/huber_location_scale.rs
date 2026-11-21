//! Independent clipped-score and scale checks for the Huber location estimates.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::robust::huber_location;
use rssn::kernels::robust::huber_location_with_scale;

const NORMAL_MAD: f64 = 1.482_602_218_505_602;

fn signed_score(
    sample: &[f64],
    x: f64,
    cutoff: f64,
) -> f64 {
    sample
        .iter()
        .map(|&value| (value - x).clamp(-cutoff, cutoff))
        .sum()
}

fn bisect(
    sample: &[f64],
    cutoff: f64,
) -> f64 {
    let mut left = sample.iter().copied().fold(f64::INFINITY, f64::min);
    let mut right = sample.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    for _ in 0..160 {
        let x = f64::midpoint(left, right);
        if x == left || x == right {
            return x;
        }
        let score = signed_score(sample, x, cutoff);
        if score > 0.0 {
            left = x;
        } else if score < 0.0 {
            right = x;
        } else {
            return x;
        }
    }
    f64::midpoint(left, right)
}

fn sorted_mad(sample: &[f64]) -> f64 {
    let mut ordered = sample.to_vec();
    ordered.sort_by(f64::total_cmp);
    let midpoint = ordered[ordered.len() / 2];
    for value in &mut ordered {
        *value = (*value - midpoint).abs();
    }
    ordered.sort_by(f64::total_cmp);
    ordered[ordered.len() / 2]
}

#[test]
fn a_known_active_set_gives_the_exact_clipped_score_root() {
    let sample = [-3.0, -1.0, 0.0, 2.0, 100.0];
    let expected = 1.0 / 3.0;
    let actual = huber_location_with_scale(&sample, 2.0, 1.0);
    assert!((actual - expected).abs() < 1e-15);
    assert!(signed_score(&sample, actual, 2.0).abs() < 1e-14);
    assert_eq!(huber_location_with_scale(&sample, 1.0, 2.0), actual);
}

#[test]
fn native_huber_preserves_the_sample_dimension_and_checks_scale_units() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let observations = [0.0, 0.0, 0.2, 0.4, 100.0];
    let samples = Value::Tensor(TensorValue {
        shape: vec![observations.len()],
        elem: TypeExpr::Float.with_dim(length.clone()),
        data: observations.into_iter().map(Value::Float).collect(),
    });
    let estimate = |dim| {
        session.call_value(
            "num_robust_huber_location_with_scale",
            &[
                samples.clone(),
                Value::Float(1.0),
                Value::Dyn(DynValue {
                    value: Box::new(Value::Float(1.0)),
                    dim,
                }),
            ],
            &ctx,
        )
    };
    let Value::Dyn(estimate_value) = estimate(length.clone()).unwrap() else {
        panic!("expected a measured location");
    };
    assert_eq!(estimate_value.dim, length);
    let Value::Float(coherent) = *estimate_value.value else {
        panic!("expected a scalar location");
    };
    assert!((coherent - 0.4).abs() < 1e-12);
    assert!(estimate(DimExpr::base("Time")).is_err());
}

#[test]
fn tuned_and_mad_based_estimates_match_independent_bisection() {
    let mut seed = 32_u64;
    for count in [5, 17, 129, 513] {
        let mut sample = Vec::new();
        for index in 0..count {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let base = ((seed >> 11) as f64 / 2.0_f64.powi(53)).powi(2);
            sample.push((base * 4.0 - 1.0) + if index % 13 == 0 { 25.0 } else { 0.0 });
        }
        let scale = sorted_mad(&sample) * NORMAL_MAD;
        for k in [0.005, 0.05, 0.4, 1.345, 3.0, 20.0, 100.0] {
            let expected = bisect(&sample, k * scale);
            let actual = huber_location(&sample, k);
            assert!(
                (actual - expected).abs() < 6e-13,
                "n={count}, k={k}, {actual} vs {expected}"
            );
            let supplied = huber_location_with_scale(&sample, k, 0.75);
            let expected_supplied = bisect(&sample, k * 0.75);
            assert!((supplied - expected_supplied).abs() < 6e-13);
        }
    }
}

#[test]
fn huge_mad_and_cutoff_products_keep_a_finite_location() {
    let sample = [-1.0, -0.875, 0.0, 0.75, 1.0].map(|x| x * f64::MAX);
    assert!((0.875 * f64::MAX * NORMAL_MAD).is_infinite());
    let actual = huber_location(&sample, 1.345);
    assert!(actual.is_finite());
    assert!((actual / f64::MAX + 0.025).abs() < 1e-15);

    let tiny_unit = f64::from_bits(1);
    let tiny = [0.0, tiny_unit, 5.0 * tiny_unit];
    assert_eq!(
        huber_location_with_scale(&tiny, f64::MAX, f64::MAX).to_bits(),
        2
    );
    assert_eq!(huber_location_with_scale(&tiny, 1.0, f64::MAX).to_bits(), 2);

    let balanced = [
        f64::NAN,
        f64::NEG_INFINITY,
        -2.0,
        0.0,
        2.0,
        3.0,
        f64::INFINITY,
    ];
    assert_eq!(
        huber_location_with_scale(&balanced, f64::MAX, f64::MAX),
        0.75
    );
}

#[test]
fn subnormal_observations_and_clipping_keep_the_resolved_location() {
    let units = [-30.0, -10.0, 0.0, 20.0, 1_000.0];
    let min_positive = f64::from_bits(1);
    let scaled = units.map(|value| value * min_positive);
    let expected = bisect(&units, 1.345 * NORMAL_MAD * sorted_mad(&units));
    let actual = huber_location(&scaled, 1.345) / min_positive;
    assert!((actual - expected).abs() <= 0.51, "{actual} vs {expected}");

    let sample = [-2.0, 1.0, 3.0];
    assert_eq!(
        huber_location_with_scale(&sample, min_positive, min_positive),
        1.0
    );
}

#[test]
fn flat_intervals_and_missing_values_have_deterministic_semantics() {
    assert_eq!(huber_location_with_scale(&[-3.0, 1.0], 0.1, 1.0), -1.0);
    assert_eq!(huber_location(&[4.0, 4.0, 4.0, 100.0], 1.345), 4.0);
    assert_eq!(
        huber_location(&[1.0, f64::INFINITY, f64::INFINITY], 1.345),
        f64::INFINITY
    );
    assert!(huber_location(&[f64::NEG_INFINITY, f64::INFINITY], 1.345).is_nan());
    assert!(huber_location(&[f64::NAN], 1.345).is_nan());
    assert!(huber_location_with_scale(&[], 1.0, 1.0).is_nan());
    assert!(huber_location_with_scale(&[1.0], 1.0, 0.0).is_nan());
    assert!(huber_location_with_scale(&[1.0], f64::INFINITY, 1.0).is_nan());
    assert!(huber_location_with_scale(&[1.0], 1.0, f64::NAN).is_nan());
}
