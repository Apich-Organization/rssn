//! Regression coverage for the real-axis continuation shared by Lerch and polylog.

use num_complex::Complex64;
use rssn::kernels::special_complex_ext::lerch_phi_complex;
use rssn::kernels::special_ext::lerch;
use rssn::kernels::special_ext::polylog;
use rssn::kernels::special_ext::polylog_complex;

#[test]
fn real_and_complex_continuations_agree_with_reference_values() {
    // Evaluated independently at 55 decimal digits with mpmath.
    for (z, s, expected) in [
        (
            -2.934_378_196_865_409_5,
            0.805_578_583_264_250_4,
            0.426_578_979_084_660_85,
        ),
        (-1.0, 0.001, 0.500_225_760_834_162_5),
        (
            -1.639_498_081_505_755_6,
            -1.049_891_850_324_448_8,
            0.132_439_858_616_325_2,
        ),
        (-1234.0, 2.0, 0.021_861_589_304_464_496),
        (-1234.0, -1.3, -0.000_016_244_157_230_436_844),
        (0.999, 0.805, 13.581_999_233_838_04),
    ] {
        let zz = Complex64::from(z);
        let ss = Complex64::from(s);
        for actual in [
            Complex64::from(lerch(z, s, 1.0)),
            lerch_phi_complex(zz, ss, Complex64::ONE),
            Complex64::from(polylog(s, z) / z),
            polylog_complex(ss, zz) / zz,
        ] {
            assert!(
                (actual - expected).norm() <= 2e-13 * expected.abs(),
                "z={z}, s={s}: {actual} should be {expected}"
            );
        }
    }
}

#[test]
fn lerch_keeps_an_independent_shift_parameter() {
    for (order, expected) in [
        (0.5, 0.136_730_464_276_828_83),
        (-1.25, 0.103_078_727_532_341_3),
    ] {
        let real = lerch(-6.5, order, 1.7);
        let complex = lerch_phi_complex((-6.5).into(), order.into(), 1.7.into());
        assert!((real - expected).abs() < 2e-14);
        assert!((complex - expected).norm() < 2e-14);
    }
}

#[test]
fn zero_order_is_a_rational_function_across_the_cut() {
    let z = Complex64::new(2.230_490_308_529_074, 0.0);
    let a = Complex64::new(-1.687_053_938_838_183_6, 0.0);
    let expected = Complex64::ONE / (Complex64::ONE - z);
    let actual = lerch_phi_complex(z, Complex64::ZERO, a);
    assert_eq!(actual, expected);
}

#[test]
fn complex_lerch_uses_the_correct_lip_of_the_cut() {
    let order = Complex64::from(2.0);
    let shift = Complex64::from(2.0);
    for im in [0.0, 0.02, -0.02] {
        let z = Complex64::new(2.0, im);
        let expected = (polylog_complex(order, z) / z - Complex64::ONE) / z;
        let actual = lerch_phi_complex(z, order, shift);
        assert!(
            (actual - expected).norm() < 6e-14,
            "z={z}: {actual} should be {expected}"
        );
        assert_eq!(actual.im.is_sign_positive(), im > 0.0);
    }
}

#[test]
fn complex_lerch_euler_integral_matches_high_precision_references() {
    // DLMF 25.14.5, evaluated independently by high-precision quadrature.
    for (z, s, a, expected) in [
        (
            (1.4, 0.2),
            (1.2, 0.3),
            (2.0, 0.1),
            (0.672_591_447_771_341_5, 1.208_044_549_182_517),
        ),
        (
            (0.95, 0.01),
            (0.001, 0.002),
            (0.3, 0.9),
            (19.217_924_072_724_835, 3.738_772_004_659_223_6),
        ),
        (
            (199.3, -0.0002),
            (3.0, 2.0),
            (1.0, 0.2),
            (0.036_245_802_981_909_6, -0.293_526_217_701_394),
        ),
    ] {
        let z = Complex64::new(z.0, z.1);
        let s = Complex64::new(s.0, s.1);
        let a = Complex64::new(a.0, a.1);
        let expected = Complex64::new(expected.0, expected.1);
        let actual = lerch_phi_complex(z, s, a);
        assert!(
            (actual - expected).norm() < 8e-14 * expected.norm().max(1.0),
            "z={z}, s={s}, a={a}: {actual} should be {expected}"
        );
    }
}
