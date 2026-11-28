//! Branch-cut and near-origin checks for complex Lambert W.

use num_complex::Complex64;
use rssn::kernels::special_complex::lambert_w_complex;

#[test]
fn signed_zero_selects_the_correct_branch_at_the_branch_point() {
    let point = -0.367_879_441_171_442_33_f64;
    let upper = Complex64::new(point, 0.0);
    let lower = upper.conj();
    assert_eq!(lambert_w_complex(upper, -1).re, -1.0);
    assert_eq!(lambert_w_complex(lower, 1).re, -1.0);
    assert!(lambert_w_complex(lower, 0).im.is_sign_negative());

    let upper_distant = lambert_w_complex(upper, 1);
    let lower_distant = lambert_w_complex(lower, -1);
    let want = Complex64::new(-3.088_843_015_613_044, -7.461_489_285_654_254);
    assert!((lower_distant - want).norm() <= 1e-14);
    assert!((upper_distant - want.conj()).norm() <= 1e-14);
    assert!((lower_distant * lower_distant.exp() - lower).norm() <= 1e-14);

    let right = point.next_up();
    let adjacent_lower = lambert_w_complex(Complex64::new(right, -0.0), 1);
    assert!((adjacent_lower.re + 1.000_000_015_304_254_3).abs() <= 2.0 * f64::EPSILON);
    assert_eq!(adjacent_lower.im, 0.0);
}

#[test]
fn small_principal_values_keep_their_tiny_components() {
    // Rounded references from 130-digit arithmetic, checked componentwise.
    for (z, expected) in [
        (
            Complex64::new(0.0, 5.95e-8),
            Complex64::new(3.540_249_999_999_966_4e-15, 5.949_999_999_999_968e-8),
        ),
        (Complex64::new(0.0, 1e-100), Complex64::new(1e-200, 1e-100)),
        (
            Complex64::new(5.95e-8, -2.1e-8),
            Complex64::new(5.949_999_690_075_02e-8, -2.099_999_750_100_032e-8),
        ),
    ] {
        let result = lambert_w_complex(z, 0);
        for (value, reference) in [(result.re, expected.re), (result.im, expected.im)] {
            assert!((value - reference).abs() <= 4.0 * f64::EPSILON * reference.abs());
        }
    }
    let tiny = Complex64::new(f64::from_bits(1), -0.0);
    let result = lambert_w_complex(tiny, 0);
    assert_eq!(result.re.to_bits(), tiny.re.to_bits());
    assert_eq!(result.im.to_bits(), tiny.im.to_bits());
}

#[test]
fn small_principal_values_preserve_the_real_residual_under_cancellation() {
    for (z, expected_real) in [
        (
            Complex64::new(-1.000_000_000_000_000_1e-20, 1e-10),
            -2.274_382_039_513_355_5e-37,
        ),
        (
            Complex64::new(-3.530_543_912_777_516e-15, 5.941_838_026_046_743e-8),
            -1.928_202_531_690_729e-31,
        ),
    ] {
        let result = lambert_w_complex(z, 0);
        let relative_error = (result.re / expected_real - 1.0).abs();
        assert!(relative_error <= 1e-13, "{z}: got {result}");
    }
}
