//! Positive-parameter Gauss values at one across the floating-point range.

use num_complex::Complex64 as C;
use rssn::kernels::hypergeometric::HypError;
use rssn::kernels::hypergeometric::hyp_pfq_checked;
use rssn::kernels::hypergeometric::hyp2f1_cplx;

fn f(
    a: f64,
    b: f64,
    c: f64,
) -> C {
    hyp2f1_cplx(C::new(a, 0.0), C::new(b, 0.0), C::new(c, 0.0), C::ONE)
}

#[test]
fn gamma_recurrence_at_large_arguments_and_one_ulp_above_divergence() {
    // Gauss's theorem and Gamma(x + 1) = x Gamma(x) give F(1,b;c;1).
    for (b, c) in [(1.0, 200.0), (17.5, 800.0), (0.5, 1e15)] {
        let expected = (c - 1.0) / ((c - 1.0) - b);
        let actual = f(1.0, b, c);
        assert_eq!(actual.im, 0.0);
        assert!((actual.re / expected - 1.0).abs() < 8e-14);
        assert_eq!(actual, f(b, 1.0, c));
    }
    // Computing c - (a + b) first rounds away the positive excess.
    let b = 3.0 * 2.0_f64.powi(-54);
    let c = 1.0_f64.next_up();
    assert!((f(1.0, b, c).re - 4.0).abs() < 2e-14);
    assert_eq!(f(1.0, b, c), f(b, 1.0, c));
    assert_eq!(
        hyp_pfq_checked(&[C::ONE, C::new(b, 0.0)], &[C::new(c, 0.0)], C::ONE).unwrap(),
        f(1.0, b, c)
    );
}

#[test]
fn finite_independent_high_precision_gauss_ratios() {
    // mpmath log-gamma differences at 450 decimal digits, using the exact
    // binary inputs before rounding these results to f64.
    let minimum = f64::from_bits(1);
    for (a, b, c, expected) in [
        (0.3, 0.4, 2.5, 1.073_810_773_722_559_7),
        (123.5, 47.75, 300.125, 2_383_069_056_708.137_7),
        (1e8, 1e8, 1e16, 2.718_281_855_641_864),
        (1e150, 1.5e150, 1e308, 1.000_000_015_000_000_1),
        (
            0.75,
            9_007_199_254_740_992.0,
            9_007_199_254_740_994.0,
            838_037_420_390.968_5,
        ),
        (0.25, 0.5, 0.75_f64.next_up(), 1_717_582_293_339_987.8),
        (0.75, 1.125, 2.375_f64.next_down(), 2.665_813_791_000_11),
        (0.75, 1.125, 2.375_f64.next_up(), 2.665_813_791_000_106_7),
        (minimum, minimum, 4.0 * minimum, 1.125),
    ] {
        let actual = f(a, b, c);
        assert_eq!(actual.im, 0.0);
        assert!(
            (actual.re / expected - 1.0).abs() < 8e-14,
            "F({a},{b};{c};1): {actual}, expected {expected}"
        );
        assert_eq!(actual, f(b, a, c));
    }
}

#[test]
fn retains_divergence_poles_and_termination_at_unity() {
    assert_eq!(f(1.0, 1.0, 2.0), C::new(f64::INFINITY, 0.0));
    assert_eq!(f(1.0, 1.0, 2.0_f64.next_down()), C::new(f64::INFINITY, 0.0));
    assert_eq!(
        hyp_pfq_checked(&[C::ONE, C::ONE], &[C::new(-2.0, 0.0)], C::ONE),
        Err(HypError::DenominatorPole)
    );
    assert_eq!(f(-2.0, 0.5, 2.0), C::new(0.625, 0.0));
    assert_eq!(f(0.0, 1.0, 2.0), C::ONE);
}
