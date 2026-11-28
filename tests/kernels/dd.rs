//! Tests for double-double arithmetic.

use rssn::kernels::dd::Cdd;
use rssn::kernels::dd::Dd;

fn atan_inv(k: f64) -> Dd {
    // atan(1/k) = sum (-1)^j / ((2j+1) k^(2j+1))
    let kk = Dd::from_f64(k);
    let k2 = kk * kk;
    let mut pw = Dd::ONE / kk;
    let mut sum = Dd::ZERO;
    for j in 0..60 {
        let term = pw / Dd::from_f64(f64::from(2 * j + 1));
        if j % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
        pw /= k2;
    }
    sum
}

#[test]
fn real_arithmetic_reaches_thirty_digits() {
    let third = Dd::ONE / Dd::from_f64(3.0);
    let err = (third * Dd::from_f64(3.0) - Dd::ONE).abs();
    assert!(err.hi < 1e-31, "{err:?}");
    let s2 = Dd::from_f64(2.0).sqrt();
    assert!((s2 * s2 - Dd::from_f64(2.0)).abs().hi < 1e-31);
    // Machin: pi = 16 atan(1/5) - 4 atan(1/239)
    let pi = Dd::from_f64(16.0) * atan_inv(5.0) - Dd::from_f64(4.0) * atan_inv(239.0);
    let want = Dd {
        hi: std::f64::consts::PI,
        lo: 1.224_646_799_147_353_2e-16,
    };
    assert!((pi - want).abs().hi < 2e-31, "{pi:?}");
    // e = sum 1/k!
    let mut e = Dd::ZERO;
    let mut t = Dd::ONE;
    for k in 1..40 {
        e += t;
        t /= Dd::from_f64(f64::from(k));
    }
    let want = Dd {
        hi: std::f64::consts::E,
        lo: 1.445_646_891_729_250_2e-16,
    };
    assert!((e - want).abs().hi < 2e-31, "{e:?}");
    // 1 + 2^-80 is not representable in f64 but is in Dd
    let tiny = Dd::from_f64(2.0_f64.powi(-80));
    let x = Dd::ONE + tiny;
    assert!((x - Dd::ONE - tiny).abs().hi == 0.0);
    assert!(Dd::from_f64(-4.0).sqrt().hi.is_nan());
    assert!(Dd::from_f64(-3.5).abs().hi == 3.5);
}

#[test]
fn complex_arithmetic() {
    let i = Cdd::new(0.0, 1.0);
    let m = i * i + Cdd::new(1.0, 0.0);
    assert!(m.re.hi == 0.0 && m.im.hi == 0.0);
    let a = Cdd::new(1.0 / 3.0, -2.0 / 7.0);
    let b = Cdd::new(0.7, 1.3);
    let r = (a / b) * b - a;
    assert!(r.norm() < 1e-31, "{}", r.norm());
    let mut z = a;
    z += b;
    z -= b;
    z *= b;
    z /= b;
    assert!((z - a).norm() < 1e-31);
    assert!((-a + a).norm() == 0.0);
}

#[test]
fn square_root_recovers_the_low_part_at_both_float_limits() {
    // Independently rounded high and residual from a 1800-digit square root.
    for (input, high, low) in [
        (
            f64::from_bits(3),
            0x1e6b_b67a_e858_4caa,
            0x1b0c_ec95_d0b5_c1e3,
        ),
        (f64::MAX, 0x5fef_ffff_ffff_ffff, 0x5c90_0000_0000_0000),
    ] {
        let root = Dd::from_f64(input).sqrt();
        assert_eq!(root.hi.to_bits(), high);
        assert_eq!(root.lo.to_bits(), low);
    }
    assert_eq!(
        Dd::from_f64(f64::INFINITY).sqrt(),
        Dd::from_f64(f64::INFINITY)
    );
    assert!(Dd::from_f64(f64::NAN).sqrt().hi.is_nan());
}

#[test]
fn euclidean_norm_keeps_both_parts_and_scales_its_squares() {
    let scale = Dd::from_parts(1.0, 2.0_f64.powi(-60));
    let c = Cdd {
        re: Dd::from_f64(3.0) * scale,
        im: Dd::from_f64(-4.0) * scale,
    };
    let want = Dd::from_f64(5.0) * scale;
    assert!((c.norm_dd() - want).abs().hi < 2e-31);
    assert_ne!(c.norm_dd().lo, 0.0);
    assert_eq!(Dd::ZERO.hypot(-scale), scale);
    assert_eq!(Cdd::new(-0.0, -0.0).norm_dd().hi.to_bits(), 0);

    for scale in [f64::from_bits(1), 2.0_f64.powi(1000)] {
        assert_eq!(
            Dd::from_f64(3.0 * scale).hypot(Dd::from_f64(4.0 * scale)),
            Dd::from_f64(5.0 * scale)
        );
    }
    let large = Dd::from_f64(f64::MAX);
    assert_eq!(large.hypot(large), Dd::from_f64(f64::INFINITY));
}

#[test]
fn euclidean_norm_uses_the_low_part_to_round_subnormal_midpoints() {
    // In units of the minimum subnormal, each norm is within 1e-15 of a
    // half integer. Comparing x^2 + y^2 with (x + 1/2)^2 fixes the direction.
    let y = Dd::from_f64(f64::from_bits(1 << 25));
    for (x, expected) in [
        ((1 << 50) - 2, (1 << 50) - 1),
        ((1 << 50) + 1, (1 << 50) + 1),
    ] {
        let norm = Dd::from_f64(f64::from_bits(x)).hypot(y);
        assert_eq!(norm, Dd::from_f64(f64::from_bits(expected)));
    }
}
