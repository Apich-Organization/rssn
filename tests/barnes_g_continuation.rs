//! Numerical and branch regressions for the complex Barnes G-function.

use num_complex::Complex64;
use rssn::kernels::special_complex::gamma_complex;
use rssn::kernels::special_complex::ln_gamma_complex;
use rssn::kernels::special_complex_ext::barnes_g_complex;
use rssn::kernels::special_complex_ext::ln_barnes_g_complex;
use std::f64::consts::PI;

fn assert_close(
    actual: Complex64,
    expected: Complex64,
) {
    let tolerance = 3e-13 * (1.0 + expected.norm());
    assert!(
        (actual - expected).norm() <= tolerance,
        "got {actual}, expected {expected} within {tolerance}"
    );
}

#[test]
fn logarithmic_barnes_g_matches_reference_across_the_shift_boundary() {
    // Computed at 80 decimal digits with the Hurwitz zeta derivative and
    // Kinkelin reflection, and independently checked after exponentiation.
    for (z, expected) in [
        (
            Complex64::new(0.5, 0.0),
            Complex64::new(-0.505_433_054_489_695_3, 0.0),
        ),
        (
            Complex64::new(2.5, 1.1),
            Complex64::new(-0.340_455_400_328_172_56, -0.162_724_926_797_941_07),
        ),
        (
            Complex64::new(-7.999, 0.17),
            Complex64::new(20.593_146_874_625_443, 125.327_797_308_665_04),
        ),
        (
            Complex64::new(-8.001, 0.17),
            Complex64::new(20.615_023_971_887_42, 125.442_811_645_587_26),
        ),
        (
            Complex64::new(-128.25, 0.1),
            Complex64::new(28_017.504_466_965_955, 26_254.256_991_687_274),
        ),
        (
            Complex64::new(-1000.0, 2_f64.powi(-40)),
            Complex64::new(2_682_954.352_025_72, 1_573_939.490_244_808),
        ),
    ] {
        assert_close(ln_barnes_g_complex(z), expected);
        assert_close(ln_barnes_g_complex(z.conj()), expected.conj());
    }
}

#[test]
fn barnes_g_preserves_the_functional_equation_and_cut_sides() {
    for z in [
        Complex64::new(2.5, 1.1),
        Complex64::new(-8.25, 0.2),
        Complex64::new(-31.125, -0.25),
    ] {
        assert_close(
            ln_barnes_g_complex(z + 1.0) - ln_barnes_g_complex(z),
            ln_gamma_complex(z),
        );
    }
    let z = Complex64::new(-8.25, 0.2);
    assert_close(
        barnes_g_complex(z + 1.0),
        gamma_complex(z) * barnes_g_complex(z),
    );

    for (real, phase) in [(-3.25, 10.0 * PI), (-8.25, 45.0 * PI)] {
        // Crossing the zeros of orders 1 through n accumulates a phase of
        // pi*n*(n+1)/2, with the sign determined by the side of the cut.
        let upper = ln_barnes_g_complex(Complex64::new(real, 0.0));
        let lower = ln_barnes_g_complex(Complex64::new(real, -0.0));
        assert!((upper.im - phase).abs() <= 2e-13);
        assert_close(lower, upper.conj());
    }

    let zero = Complex64::new(0.0, 0.0);
    let one = Complex64::new(1.0, 0.0);
    assert_eq!(barnes_g_complex(one), one);
    assert_eq!(barnes_g_complex(Complex64::new(-9.0, 0.0)), zero);
    assert_eq!(
        ln_barnes_g_complex(Complex64::new(-9.0, 0.0)),
        Complex64::new(f64::NEG_INFINITY, 0.0)
    );
}

#[test]
fn large_finite_negative_arguments_use_the_logarithmic_continuation() {
    for (z, expected) in [
        (
            Complex64::new(-1e16, 0.25),
            Complex64::new(1.767_068_074_395_237e33, 1.570_796_326_794_896_1e32),
        ),
        (
            Complex64::new(-1e100, 1.0),
            Complex64::new(1.143_792_546_497_022_8e202, 1.570_796_326_794_896_8e200),
        ),
    ] {
        assert_close(ln_barnes_g_complex(z), expected);
        assert_close(ln_barnes_g_complex(z.conj()), expected.conj());
    }
}
