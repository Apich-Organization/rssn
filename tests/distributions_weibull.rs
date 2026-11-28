//! Weibull tails, inverse probabilities, and extreme scales.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::distributions::DistError;
use rssn::kernels::distributions::Distribution;
use rssn::kernels::distributions::weibull_inverse_survival;
use rssn::kernels::distributions::weibull_quantiles;

fn assert_relative(
    actual: f64,
    expected: f64,
) {
    assert!(
        actual.is_finite() && (actual / expected - 1.0).abs() <= 9e-13,
        "{actual:e} != {expected:e}",
    );
}

fn coherent_float(
    value: &Value,
    expected: &DimExpr,
) -> f64 {
    let Value::Dyn(measured) = value else {
        panic!("expected a measured Float, got {value:?}");
    };
    assert_eq!(&measured.dim, expected);
    let Value::Float(coherent) = *measured.value else {
        panic!("expected a coherent Float");
    };
    coherent
}

#[test]
fn native_weibull_quantiles_and_upper_tail_restore_scale_units() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let duration = DimExpr::base("Time");
    let scale = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: duration.clone(),
    });
    let probabilities = Value::Tensor(TensorValue {
        shape: vec![2],
        data: vec![Value::Float(1e-250), Value::Float(0.5)],
        elem: TypeExpr::Float,
    });
    let batch = session
        .call_value(
            "num_distributions_weibull_quantiles",
            &[Value::Float(2.0), scale.clone(), probabilities],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(batch) = batch else {
        panic!("expected a vector of measured quantiles");
    };
    assert_eq!(batch.elem, TypeExpr::Float.with_dim(duration.clone()));
    let [first, second] = batch.data.as_slice() else {
        panic!("expected both quantiles");
    };
    assert_relative(coherent_float(first, &duration), 2e-125);
    assert_relative(
        coherent_float(second, &duration),
        2.0 * std::f64::consts::LN_2.sqrt(),
    );
    let upper = session
        .call_value(
            "num_distributions_weibull_inverse_survival",
            &[Value::Float(2.0), scale, Value::Float(1e-250)],
            &ctx,
        )
        .unwrap();
    assert_relative(
        coherent_float(&upper, &duration),
        2.0 * (250.0 * std::f64::consts::LN_10).sqrt(),
    );
}

#[test]
fn weibull_cdf_and_density_keep_a_representable_hazard() {
    // References calculated from the binary64 inputs in 440-digit Decimal.
    for (x, shape, scale, cdf, pdf) in [
        (
            1e300,
            0.001,
            1e-300,
            0.981_334_375_438_481_1,
            7.430_918_980_800_158e-305,
        ),
        (
            1e-300,
            0.001,
            1e300,
            0.222_124_383_190_155_04,
            1.953_935_207_266_794_8e296,
        ),
        (7.5e-308, 1.0, 1e-310, 1.0, 1.901_684_963_470_690_3e-16),
    ] {
        let distribution = Distribution::Weibull { shape, scale };
        assert_relative(distribution.cdf(x), cdf);
        assert_relative(distribution.pdf(x), pdf);
    }
}

#[test]
fn steep_weibull_resolves_points_adjacent_to_the_scale() {
    let scale = 0.95_f64;
    let distribution = Distribution::Weibull { shape: 5e15, scale };
    for (x, cdf, pdf) in [
        (
            scale.next_down(),
            0.427_349_911_289_630_2,
            1_680_216_887_689_696.5,
        ),
        (
            scale.next_up(),
            0.833_670_578_686_138_5,
            1_570_311_681_944_650.8,
        ),
    ] {
        assert_relative(distribution.cdf(x), cdf);
        assert_relative(distribution.pdf(x), pdf);
    }
}

#[test]
fn lower_quantiles_keep_tiny_probabilities_and_cancelling_scale() {
    for (probability, shape, scale, expected) in [
        (1e-250, 2.0, 2.0, 2e-125),
        (f64::from_bits(1), 2.0, 3.0, 6.668_276_248_455_232e-162),
        (
            0.950_212_931_632_136,
            0.001,
            1e-300,
            1.322_070_819_480_701_9e177,
        ),
        (
            0.329_679_953_964_360_7,
            0.001,
            1e300,
            1.148_130_695_274_161_7e-98,
        ),
    ] {
        let d = Distribution::Weibull { shape, scale };
        assert_relative(d.quantile(probability).unwrap(), expected);
        assert_relative(d.cdf(d.quantile(probability).unwrap()), probability);
    }
}

#[test]
fn inverse_survival_keeps_the_upper_tail_and_the_pivot_for_small_shape() {
    let smallest = f64::from_bits(1);
    assert_relative(
        weibull_inverse_survival(2.0, 3.0, smallest).unwrap(),
        81.853_287_333_450_64,
    );
    let survival = 0.367_879_441_171_442_33;
    let probability = 1.0 - survival;
    let d = Distribution::Weibull {
        shape: 1e-17,
        scale: 2.0,
    };
    let expected = 0.068_198_115_576_778_99;
    assert_relative(d.quantile(probability).unwrap(), expected);
    assert_relative(
        weibull_inverse_survival(1e-17, 2.0, survival).unwrap(),
        expected,
    );
    let mut calls = 0;
    let sampled = d
        .sample(&mut || {
            calls += 1;
            survival
        })
        .unwrap();
    assert_eq!(calls, 1);
    assert_relative(sampled, expected);
}

#[test]
fn weibull_quantile_batch_keeps_order_and_validates_inputs() {
    let d = Distribution::Weibull {
        shape: 1.7,
        scale: 2.3,
    };
    let probabilities = [0.7, 0.0, 1.0, 1e-45, 0.7, 0.1];
    let expected = probabilities.map(|p| d.quantile(p).unwrap());
    assert_eq!(
        weibull_quantiles(1.7, 2.3, &probabilities).unwrap(),
        expected
    );
    assert!(weibull_quantiles(1.7, 2.3, &[]).unwrap().is_empty());
    assert_eq!(
        weibull_inverse_survival(1.7, 2.3, 0.0).unwrap(),
        f64::INFINITY
    );
    assert_eq!(weibull_inverse_survival(1.7, 2.3, 1.0).unwrap(), 0.0);
    for invalid in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            weibull_quantiles(1.7, 2.3, &[0.5, invalid]),
            Err(DistError::Domain(_))
        ));
        assert!(matches!(
            weibull_inverse_survival(1.7, 2.3, invalid),
            Err(DistError::Domain(_))
        ));
    }
    for (shape, scale) in [
        (0.0, 2.3),
        (1.7, 0.0),
        (f64::INFINITY, 2.3),
        (1.7, f64::NAN),
    ] {
        assert!(matches!(
            weibull_quantiles(shape, scale, &[]),
            Err(DistError::Domain(_))
        ));
        assert!(matches!(
            weibull_inverse_survival(shape, scale, 0.5),
            Err(DistError::Domain(_))
        ));
    }
}

#[test]
fn weibull_density_endpoints_keep_their_existing_conventions() {
    for (shape, at_zero) in [(0.5, f64::INFINITY), (1.0, 0.5), (2.0, 0.0)] {
        let d = Distribution::Weibull { shape, scale: 2.0 };
        assert_eq!(d.pdf(0.0), at_zero);
        assert_eq!(d.pdf(-1.0), 0.0);
        assert_eq!(d.pdf(f64::INFINITY), 0.0);
        assert_eq!(d.cdf(0.0), 0.0);
        assert_eq!(d.cdf(f64::INFINITY), 1.0);
        assert!(d.cdf(f64::NAN).is_nan());
    }
}

#[test]
fn ordinary_weibull_matches_exponential_and_rayleigh_reductions() {
    for (shape, hazard, rate) in [(1.0, 1.5_f64, 0.5), (2.0, 2.25, 1.5)] {
        let distribution = Distribution::Weibull { shape, scale: 2.0 };
        let survival = (-hazard).exp();
        assert_relative(distribution.cdf(3.0), 1.0 - survival);
        assert_relative(distribution.pdf(3.0), rate * survival);
        assert_relative(distribution.quantile(1.0 - survival).unwrap(), 3.0);
        assert_relative(weibull_inverse_survival(shape, 2.0, survival).unwrap(), 3.0);
    }
}
