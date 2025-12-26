use assert_approx_eq::assert_approx_eq;
use rssn::prelude::numerical::*;

#[test]

fn test_aitken_acceleration() {

    // Sequence converging to 1: s_n = 1 + 0.5^n
    // This is linearly convergent. Aitken should capture it perfectly or very well.
    let mut seq = Vec::new();

    for i in 0 .. 10 {

        seq.push(1.0 + 0.5f64.powi(i));
    }

    let acc =
        numerical_aitken_acceleration(
            &seq,
        );

    let last = acc.last().unwrap();

    // 1 + 0.5^9 is ~1.00195
    // Aitken on geometric series is actually exact (or very close).
    assert_approx_eq!(last, 1.0, 1e-10);
}

#[test]

fn test_richardson_extrapolation() {

    // Derivative of e^x at x=0 is 1.
    // Central difference: (e^h - e^-h)/(2h). Error O(h^2).
    // Steps: 0.8, 0.4, 0.2, 0.1
    let f = |x : f64| x.exp();

    let derivative = |h : f64| {

        (f(h) - f(-h)) / (2.0 * h)
    };

    let steps =
        vec![0.4, 0.2, 0.1, 0.05];

    let approximations : Vec<f64> =
        steps
            .iter()
            .map(|&h| derivative(h))
            .collect();

    let extrapolated = numerical_richardson_extrapolation(&approximations);

    let best = extrapolated
        .last()
        .unwrap();

    // Normal approximation for h=0.05
    let normal_err =
        (derivative(0.05) - 1.0).abs();

    // Extrapolated error should be significantly smaller
    let extrap_err = (best - 1.0).abs();

    assert!(extrap_err < normal_err);

    assert!(
        extrap_err < 1e-9,
        "Extrapolated error: {}, \
         Normal error: {}",
        extrap_err,
        normal_err
    );
}

#[test]

fn test_wynn_epsilon() {

    // Leibniz series for pi/4: 1 - 1/3 + 1/5 - 1/7 ...
    // Converges very slowly.
    let mut seq = Vec::new();

    let mut sum = 0.0;

    for i in 0 .. 10 {

        let term = if i % 2 == 0 {

            1.0
        } else {

            -1.0
        } / (2.0 * i as f64
            + 1.0);

        sum += term;

        seq.push(sum);
    }

    // Last term of raw sequence
    let raw_err = (sum
        - std::f64::consts::FRAC_PI_4)
        .abs();

    let acc =
        numerical_wynn_epsilon(&seq);

    let best = acc.last().unwrap();

    let acc_err = (best
        - std::f64::consts::FRAC_PI_4)
        .abs();

    // Wynn epsilon is very effective for alternating series
    assert!(acc_err < raw_err);

    // With 10 terms, raw error is ~1/20 = 0.05.
    // Wynn is usually much better.
    assert!(acc_err < 1e-4);
}

#[cfg(test)]

mod proptests {

    use assert_approx_eq::assert_approx_eq;
    use proptest::prelude::*;

    use super::*;

    proptest! {
        // Test that Aitken acceleration preserves the limit of an already convergent constant sequence
        #[test]
        fn test_aitken_constant_sequence(c in -100.0..100.0f64) {
            // Sequence: c + 1/2^n
            let mut s = Vec::new();
            for i in 0..10 {
                s.push(c + 0.5f64.powi(i));
            }
            let res = numerical_aitken_acceleration(&s);
            if let Some(last) = res.last() {
                assert_approx_eq!(last, c, 1e-9);
            }
        }

        // Test that Richardson extrapolation improves or maintains accuracy for a simple function
        #[test]
        fn test_richardson_consistency(val in 0.1..2.0f64) {
            // f(x) = x^2, f'(val) = 2*val.
            // Central diff is exact for quadratic. So Richardson should also be exact.

            let f = |x: f64| x * x;
            let derivative = |h: f64| (f(val + h) - f(val - h)) / (2.0 * h);

            let steps = vec![0.1, 0.05, 0.025];
            let approxs: Vec<f64> = steps.iter().map(|&h| derivative(h)).collect();
            let extraps = numerical_richardson_extrapolation(&approxs);

            if let Some(best) = extraps.last() {
                assert_approx_eq!(best, 2.0 * val, 1e-9);
            }
        }
    }
}
