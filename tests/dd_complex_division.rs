//! Exponent-range and cancellation regressions for complex double-double division.

use rssn::kernels::dd::Cdd;
use rssn::kernels::dd::Dd;

#[test]
fn divides_without_materializing_the_squared_norm_at_input_scale() {
    for s in [1e-200, 1e200, f64::from_bits(1), f64::MAX] {
        let z = Cdd::new(s, s);
        assert_eq!(z / Cdd::new(s, s), Cdd::new(1.0, 0.0));
        assert_eq!(Cdd::new(s, 0.0) / z, Cdd::new(0.5, -0.5));
        let axis = Cdd::new(s, 0.0);
        assert_eq!(axis / Cdd::new(s, 0.0), Cdd::new(1.0, 0.0));
    }
}

#[test]
fn each_output_component_retains_its_own_scale() {
    let large = 2.0_f64.powi(1000);
    let small = 2.0_f64.powi(-1000);
    let z = Cdd::new(large, small);
    assert_eq!(z / Cdd::new(2.0, 0.0), Cdd::new(0.5 * large, 0.5 * small));
    assert_eq!(z / Cdd::new(0.0, -2.0), Cdd::new(-0.5 * small, 0.5 * large));

    let a = 2.0_f64.powi(600);
    let b = 2.0_f64.powi(-400);
    let q = Cdd::new(0.0, a) / Cdd::new(a, b);
    assert_eq!(q, Cdd::new(small, 1.0));

    let q = Cdd::new(-large, 1.0) / Cdd::new(small, 0.0);
    assert_eq!(q, Cdd::new(f64::NEG_INFINITY, large));

    let boundary = Cdd {
        re: Dd::from_parts(f64::MAX, 2.0_f64.powi(969)),
        im: Dd::ZERO,
    };
    assert_eq!(boundary / Cdd::new(1.0, 0.0), boundary);
}

#[test]
fn cancellation_keeps_the_low_low_product() {
    let leading = 2.0_f64.powi(700);
    let tail = 2.0_f64.powi(646);
    let z = Cdd {
        re: Dd::from_parts(leading, tail),
        im: Dd::from_f64(leading),
    };
    let w = Cdd {
        re: Dd::from_f64(leading),
        im: Dd::from_parts(leading, -tail),
    };
    let q = z / w;
    // With u = 2^-54 the exact imaginary numerator is 1 - (1+u)(1-u) = u^2.
    // Its quotient below is rounded from 280-digit Decimal arithmetic.
    let expected = Dd::from_parts(
        f64::from_bits(0x3920_0000_0000_0000),
        f64::from_bits(0x35c0_0000_0000_0000),
    );
    assert!((q.im - expected).abs().hi <= expected.hi * 2e-31, "{q:?}");
    assert_eq!(q.re, Dd::from_parts(1.0, 2.0_f64.powi(-54)));
}

#[test]
fn subnormal_quotients_round_to_even_and_use_the_low_part() {
    let t = f64::from_bits(1);
    assert_eq!(
        Cdd::new(3.0 * t, 5.0 * t) / Cdd::new(2.0, 0.0),
        Cdd::new(2.0 * t, 2.0 * t)
    );
    let high = 2.0_f64.powi(-54);
    let tail = 2.0_f64.powi(-108);
    let z = Cdd {
        re: Dd::from_parts(high, tail),
        im: Dd::from_parts(high, -tail),
    };
    let q = z / Cdd::new(2.0_f64.powi(1021), 0.0);
    assert_eq!(q, Cdd::new(t, 0.0));
}

#[test]
fn reciprocal_reuses_scaled_division_and_preserves_singular_results() {
    for e in [-700, 700] {
        let s = 2.0_f64.powi(e);
        let expected = 2.0_f64.powi(-1 - e);
        assert_eq!(Cdd::new(s, s).recip(), Cdd::new(expected, -expected));
    }
    assert_eq!(
        Cdd::new(f64::from_bits(1), 0.0).recip(),
        Cdd::new(f64::INFINITY, 0.0)
    );
    let small = f64::from_bits(1 << 49);
    assert_eq!(
        Cdd::new(f64::MAX, f64::MAX).recip(),
        Cdd::new(small, -small)
    );
    let singular = Cdd::new(0.0, 0.0).recip();
    assert!(singular.re.hi.is_nan() && singular.im.hi.is_nan());
}
