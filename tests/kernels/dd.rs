//! Tests for double-double arithmetic.

use rssn::kernels::dd::{Cdd, Dd};

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
    let want = Dd { hi: 3.141_592_653_589_793, lo: 1.224_646_799_147_353_2e-16 };
    assert!((pi - want).abs().hi < 2e-31, "{pi:?}");
    // e = sum 1/k!
    let mut e = Dd::ZERO;
    let mut t = Dd::ONE;
    for k in 1..40 {
        e += t;
        t /= Dd::from_f64(f64::from(k));
    }
    let want = Dd { hi: 2.718_281_828_459_045, lo: 1.445_646_891_729_250_2e-16 };
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
