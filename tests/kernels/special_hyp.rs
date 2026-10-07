//! Tests for near-integer `c - a - b` and complex-argument `2F1`.

use num_complex::Complex64;
use rssn::kernels::special_ext::{hyp2f1, hyp2f1_complex};

fn series(a: f64, b: f64, c: f64, x: f64) -> f64 {
    let (mut t, mut s) = (1.0, 1.0);
    for n in 0..200_000 {
        let nf = f64::from(n);
        t *= (a + nf) * (b + nf) / ((c + nf) * (nf + 1.0)) * x;
        s += t;
        if t.abs() < 1e-18 * s.abs() {
            break;
        }
    }
    s
}

#[test]
fn near_integer_c_minus_a_minus_b() {
    for &eps in &[1e-12, 1e-10, 1e-9, 1e-8, 1e-6, 1e-5, -1e-9, -1e-7, 5e-5] {
        for &(a, b, m, x) in &[
            (0.5, 0.75, 0.0, 0.9),
            (1.0, 1.0, 0.0, 0.95),
            (0.3, 0.6, 1.0, 0.85),
            (1.5, 0.25, 2.0, 0.8),
            (0.5, 0.5, 0.0, 0.99),
        ] {
            let c = a + b + m + eps;
            let want = series(a, b, c, x);
            let got = hyp2f1(a, b, c, x);
            assert!(
                (got - want).abs() <= 1e-12 * want.abs(),
                "a={a} b={b} m={m} eps={eps} x={x}: {got} vs {want}"
            );
        }
    }
    // negative integer c - a - b
    let (a, b, x) = (1.2, 0.9, 0.9);
    let c = a + b - 1.0 + 3e-9;
    assert!((hyp2f1(a, b, c, x) - series(a, b, c, x)).abs() < 1e-10 * series(a, b, c, x).abs());
    // continuity across the exact-integer logarithmic case
    let exact = hyp2f1(0.5, 0.5, 1.0, 0.9);
    let near = hyp2f1(0.5, 0.5, 1.0 + 1e-9, 0.9);
    assert!((exact - near).abs() < 1e-8);
}

#[test]
fn complex_argument_closed_forms() {
    let pts = [
        Complex64::new(2.0, 1.0),
        Complex64::new(-3.0, 0.5),
        Complex64::new(0.4, -2.0),
        Complex64::new(-0.9, -0.1),
        Complex64::new(10.0, 10.0),
        Complex64::new(0.3, 0.9),
        Complex64::new(1.0, 1e-3),
        Complex64::new(-40.0, 3.0),
        Complex64::new(1e4, -2.0),
    ];
    for &z in &pts {
        // 2F1(1,1;2;z) = -ln(1-z)/z
        let want = -(Complex64::new(1.0, 0.0) - z).ln() / z;
        let got = hyp2f1_complex(1.0, 1.0, 2.0, z);
        assert!((got - want).norm() <= 1e-12 * want.norm(), "z={z}: {got} vs {want}");
        // 2F1(a,b;b;z) = (1-z)^(-a)
        let want = (Complex64::new(1.0, 0.0) - z).powf(-0.7);
        let got = hyp2f1_complex(0.7, 1.3, 1.3, z);
        assert!((got - want).norm() <= 1e-12 * want.norm(), "z={z}: {got} vs {want}");
        // 2F1(1/2,1/2;3/2;z) = asin(sqrt z)/sqrt z  (principal branch away from the cut)
        if z.im.abs() > 0.05 {
            let sq = z.sqrt();
            let want = sq.asin() / sq;
            let got = hyp2f1_complex(0.5, 0.5, 1.5, z);
            assert!((got - want).norm() <= 1e-11 * want.norm(), "z={z}: {got} vs {want}");
        }
    }
}

#[test]
fn complex_matches_real_and_cut_convention() {
    for &x in &[-5.0, -0.7, 0.3, 0.8, 0.97] {
        let r = hyp2f1(0.4, 1.1, 2.3, x);
        let c = hyp2f1_complex(0.4, 1.1, 2.3, Complex64::new(x, 0.0));
        assert!((c.re - r).abs() <= 1e-12 * r.abs(), "x={x}: {c} vs {r}");
        assert!(c.im.abs() <= 1e-12 * r.abs());
    }
    // on the cut the value is the limit from above: 2F1(1,1;2;2) = i pi / 2
    let v = hyp2f1_complex(1.0, 1.0, 2.0, Complex64::new(2.0, 0.0));
    assert!((v - Complex64::new(0.0, std::f64::consts::FRAC_PI_2)).norm() < 1e-12, "{v}");
    let below = hyp2f1_complex(1.0, 1.0, 2.0, Complex64::new(2.0, -1e-12));
    assert!((below + Complex64::new(0.0, std::f64::consts::FRAC_PI_2)).norm() < 1e-9);
    // conjugate symmetry off the cut
    let z = Complex64::new(-2.5, 1.7);
    let (u, w) = (hyp2f1_complex(0.3, 0.9, 1.7, z), hyp2f1_complex(0.3, 0.9, 1.7, z.conj()));
    assert!((u - w.conj()).norm() < 1e-13 * u.norm());
    // terminating series for |z| > 1, and failure modes
    let p = hyp2f1_complex(-2.0, 3.0, 1.5, Complex64::new(1.5, 2.0));
    let z = Complex64::new(1.5, 2.0);
    let want = Complex64::new(1.0, 0.0) - z * 4.0 + z * z * 3.2;
    assert!((p - want).norm() < 1e-12);
    assert!(hyp2f1_complex(0.5, 0.5, -1.0, z).re.is_nan());
    // z = 1 delegates to Gauss's theorem
    assert!((hyp2f1_complex(1.0, 2.0, 4.0, Complex64::new(1.0, 0.0)).re - 3.0).abs() < 1e-12);
}
