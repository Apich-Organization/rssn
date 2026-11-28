//! Range and branch regressions for the single-pole Meijer G exponential.

#![allow(clippy::unwrap_used)]

use std::f64::consts::FRAC_PI_4;

use num_complex::Complex64 as C;
use rssn::kernels::meijer::meijer_g;
use rssn::kernels::meijer::meijer_g_real;

fn eval(
    b: C,
    z: C,
) -> Option<C> {
    meijer_g(1, 0, &[], &[b], z)
}

fn relative(
    got: f64,
    want: f64,
    tolerance: f64,
) {
    assert!(
        (got / want - 1.0).abs() <= tolerance,
        "got {got:e}, want {want:e}"
    );
}

#[test]
fn combined_logarithm_resolves_cancellation_and_compensated_range() {
    relative(
        meijer_g_real(1, 0, &[], &[0.5], 80.0).unwrap(),
        80.0_f64.sqrt() * (-80.0_f64).exp(),
        1e-13,
    );
    // The separate factor 800^120 overflows, while exp(-800) rounds to zero.
    relative(
        meijer_g_real(1, 0, &[], &[120.0], 800.0).unwrap(),
        8.614_159_649_215_166,
        2e-13,
    );

    let tiny = f64::from_bits(1);
    relative(
        meijer_g_real(1, 0, &[], &[0.5], tiny).unwrap(),
        tiny.sqrt(),
        1e-13,
    );
}

#[test]
fn exponential_rounds_gradually_through_underflow() {
    for (z, expected) in [
        (744.0, 2),
        (745.0, 1),
        (745.133_219_101_941_1, 1),
        (745.133_219_101_941_2, 0),
        (746.0, 0),
    ] {
        assert_eq!(
            meijer_g_real(1, 0, &[], &[0.0], z).unwrap().to_bits(),
            expected
        );
    }

    let b = C::new(0.0, 0.0);
    let got = eval(b, C::new(744.0, FRAC_PI_4)).unwrap();
    assert_eq!(got.re.to_bits(), 1);
    assert_eq!(got.im.to_bits(), (1_u64 << 63) | 1);

    let got = eval(b, C::new(745.0, FRAC_PI_4)).unwrap();
    assert_eq!(got.re.to_bits(), 0);
    assert_eq!(got.im.to_bits(), 1_u64 << 63);
}

#[test]
fn principal_cut_is_upper_for_both_zero_signs() {
    let b = C::new(0.5, 0.0);
    let want = 2.0_f64.sqrt() * 2.0_f64.exp();
    for zero in [0.0, -0.0] {
        let got = eval(b, C::new(-2.0, zero)).unwrap();
        relative(got.im, want, 1e-13);
        assert!(got.re.abs() <= 1e-14 * want);
    }
    let lower = eval(b, C::new(-2.0, -1e-10)).unwrap();
    assert!(lower.im < 0.0);
    assert_eq!(meijer_g_real(1, 0, &[], &[0.5], -2.0), None);

    let even_power = 4_503_599_627_370_496.0;
    for (power, sign) in [(even_power, 1.0), (even_power + 1.0, -1.0)] {
        let got = eval(C::new(power, 0.0), C::new(-1.0, -0.0)).unwrap();
        relative(got.re, sign * 1.0_f64.exp(), 1e-14);
        assert!(got.im.abs() <= 1e-14);
    }
}

#[test]
fn complex_phases_and_large_cartesian_components_match_the_reference() {
    // Independently rounded from a 100-digit evaluation of z^b exp(-z).
    let got = eval(C::new(0.75, -0.125), C::new(0.4, 0.9)).unwrap();
    relative(got.re, 0.764_970_847_402_371_8, 1e-13);
    relative(got.im, -0.025_764_239_464_560_548, 1e-13);

    let got = eval(C::new(-1.25, 0.2), C::new(0.001, 1e14)).unwrap();
    relative(got.re, 9.820_225_938_399_63e-19, 1e-13);
    relative(got.im, 2.088_026_205_996_764_7e-18, 1e-13);

    let got = eval(C::new(0.0, 0.0), C::new(-710.0, -FRAC_PI_4)).unwrap();
    relative(got.re, 1.579_672_848_288_201_5e308, 1e-13);
    relative(got.im, 1.579_672_848_288_201_3e308, 1e-13);
    assert!(eval(C::new(0.0, 0.0), C::new(-712.0, 0.0)).is_none());
}

#[test]
fn exact_zero_limits_and_an_overflowing_argument_norm() {
    let zero = C::new(0.0, 0.0);
    assert_eq!(eval(zero, zero), Some(C::new(1.0, 0.0)));
    assert_eq!(eval(C::new(0.5, 0.75), zero), Some(zero));
    assert_eq!(eval(C::new(0.0, 0.75), zero), None);
    assert_eq!(eval(C::new(-0.5, 0.0), zero), None);
    assert_eq!(
        eval(C::new(0.5, 0.0), C::new(f64::MAX, f64::MAX)),
        Some(zero)
    );
    assert_eq!(
        meijer_g(0, 0, &[], &[C::new(0.5, 0.0)], C::new(1.0, 0.0)),
        Some(zero)
    );
}
