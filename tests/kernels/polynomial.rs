//! Dense real polynomials (coefficients stored highest degree first). Ported
//! from `numerical_polynomial_test.rs` and `numerical/polynomial.rs`.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::polynomial::Polynomial;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn eval_and_degree() {
    let p = Polynomial {
        coeffs: vec![1.0, 2.0, 1.0],
    }; // x^2 + 2x + 1
    assert_eq!(p.eval(0.0), 1.0);
    assert_eq!(p.eval(1.0), 4.0);
    assert_eq!(p.eval(-1.0), 0.0);
    assert_eq!(p.degree(), 2);
}

#[test]
fn eval_and_derivative_high_precision() {
    // P(x) = 3x^3 - 2x^2 + 0.5x - 7 ; reference value from 50-digit arithmetic.
    let p = Polynomial {
        coeffs: vec![3.0, -2.0, 0.5, -7.0],
    };
    let x = 1.23456789f64;
    let expected = -3.7860026896706396f64;
    assert!((p.eval(x) - expected).abs() < 1e-12, "got {}", p.eval(x));

    let d = p.derivative();
    assert_eq!(d.coeffs, vec![9.0, -4.0, 0.5]);
}

#[test]
fn arithmetic() {
    let p1 = Polynomial {
        coeffs: vec![1.0, 1.0],
    }; // x + 1
    let p2 = Polynomial {
        coeffs: vec![1.0, -1.0],
    }; // x - 1
    assert_eq!((p1.clone() + p2.clone()).coeffs, vec![2.0, 0.0]);
    assert_eq!((p1.clone() * p2.clone()).coeffs, vec![1.0, 0.0, -1.0]);
    let diff = p1 - p2;
    assert!(diff.eval(5.0) - 2.0 == 0.0);
}

#[test]
fn calculus() {
    let p = Polynomial {
        coeffs: vec![1.0, 0.0, 0.0],
    }; // x^2
    assert_eq!(p.derivative().coeffs, vec![2.0, 0.0]);
    assert_eq!(p.integral().coeffs, vec![1.0 / 3.0, 0.0, 0.0, 0.0]);
}

#[test]
fn long_division_with_remainder() {
    // (x^3 - 6x^2 + 11x - 6 + 5) / (x - 1)  ->  q = x^2 - 5x + 6, r = 5
    let p = Polynomial::new(vec![1.0, -6.0, 11.0, -1.0]);
    let (q, r) = p.long_division(&Polynomial::new(vec![1.0, -1.0]));
    let want = [1.0, -5.0, 6.0];
    assert_eq!(q.coeffs.len(), 3, "quotient {:?}", q.coeffs);
    for (g, w) in q.coeffs.iter().zip(want) {
        assert_approx_eq!(*g, w, 1e-12);
    }
    assert_eq!(r.coeffs.len(), 1);
    assert_approx_eq!(r.coeffs[0], 5.0, 1e-12);
}

#[test]
fn long_division_symmetric_quotient_and_exact_division() {
    // (x^2 - 1) / (x - 1) = x + 1 exactly (quotient happens to be palindromic).
    let (q, r) =
        Polynomial::new(vec![1.0, 0.0, -1.0]).long_division(&Polynomial::new(vec![1.0, -1.0]));
    assert_eq!(&q.coeffs[..2], &[1.0, 1.0]);
    assert!(r.is_zero(1e-12));
}

#[test]
fn zero_test() {
    assert!(Polynomial::new(vec![0.0, 1e-15]).is_zero(1e-12));
    assert!(!Polynomial::new(vec![0.0, 1.0]).is_zero(1e-12));
}

#[test]
fn roots() {
    let p = Polynomial {
        coeffs: vec![1.0, 0.0, -1.0],
    }; // x^2 - 1
    let mut roots = p.find_roots().unwrap_or_else(|e| panic!("{e}"));
    roots.sort_by(f64::total_cmp);
    assert_eq!(roots.len(), 2);
    assert_approx_eq!(roots[0], -1.0, 1e-9);
    assert_approx_eq!(roots[1], 1.0, 1e-9);
}

#[test]
fn division_by_zero_scalar_is_error() {
    let p = Polynomial {
        coeffs: vec![1.0, 2.0, 3.0],
    };
    assert!(p.clone().div_scalar(0.0).is_err());
    let half = p.div_scalar(2.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(half.coeffs, vec![0.5, 1.0, 1.5]);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_add_then_sub_recovers(
        cp in prop::collection::vec(-1e6f64..1e6, 1..6),
        cq in prop::collection::vec(-1e6f64..1e6, 1..6),
        x in -1e3f64..1e3,
    ) {
        let p = Polynomial { coeffs: cp };
        let q = Polynomial { coeffs: cq };
        let back = (p.clone() + q.clone()) - q;
        let scale = p.eval(x).abs().max(1.0);
        prop_assert!((p.eval(x) - back.eval(x)).abs() < 1e-9 * scale * 1e6);
    }

    #[test]
    fn prop_scalar_mul_div_inverse(cs in prop::collection::vec(-1e6f64..1e6, 1..6), s in -1e3f64..1e3) {
        prop_assume!(s.abs() > 1e-6);
        let p = Polynomial { coeffs: cs };
        let back = (p.clone() * s).div_scalar(s).map_err(TestCaseError::fail)?;
        for x in [0.0, 0.618_033_988_7, -std::f64::consts::PI] {
            let scale = p.eval(x).abs().max(1.0);
            prop_assert!((p.eval(x) - back.eval(x)).abs() < 1e-9 * scale);
        }
    }

    #[test]
    fn prop_derivative_matches_finite_difference(cs in prop::collection::vec(-1e3f64..1e3, 2..6), x in -10.0f64..10.0) {
        let p = Polynomial { coeffs: cs };
        let h = 1e-5;
        let numeric = (p.eval(x + h) - p.eval(x - h)) / (2.0 * h);
        let analytic = p.derivative().eval(x);
        prop_assert!((analytic - numeric).abs() < 1e-4 * analytic.abs().max(1.0) * 10.0, "{analytic} vs {numeric}");
    }

    #[test]
    fn prop_eval_of_sum(c1 in prop::collection::vec(-100.0..100.0f64, 1..10),
                        c2 in prop::collection::vec(-100.0..100.0f64, 1..10), x in -10.0..10.0f64) {
        let p1 = Polynomial { coeffs: c1 };
        let p2 = Polynomial { coeffs: c2 };
        let sum = p1.clone() + p2.clone();
        let want = p1.eval(x) + p2.eval(x);
        prop_assert!((sum.eval(x) - want).abs() < 1e-9 * want.abs().max(1.0) * 1e3);
    }

    #[test]
    fn prop_eval_of_product(c1 in prop::collection::vec(-10.0..10.0f64, 1..5),
                            c2 in prop::collection::vec(-10.0..10.0f64, 1..5), x in -2.0..2.0f64) {
        let p1 = Polynomial { coeffs: c1 };
        let p2 = Polynomial { coeffs: c2 };
        let prod = p1.clone() * p2.clone();
        let want = p1.eval(x) * p2.eval(x);
        prop_assert!((prod.eval(x) - want).abs() < 1e-9 * want.abs().max(1.0) * 1e2);
    }

    #[test]
    fn prop_integral_then_derivative_is_identity(cs in prop::collection::vec(-10.0..10.0f64, 1..6), x in -3.0..3.0f64) {
        let p = Polynomial { coeffs: cs };
        let back = p.integral().derivative();
        prop_assert!((back.eval(x) - p.eval(x)).abs() < 1e-8 * p.eval(x).abs().max(1.0));
    }
}

#[test]
fn long_division_reconstructs_dividend() {
    // p = q * d + r must hold pointwise for a non-palindromic case.
    let p = Polynomial::new(vec![2.0, -3.0, 0.5, 7.0, -4.0]);
    let d = Polynomial::new(vec![1.0, 0.0, 2.0]);
    let (q, r) = p.clone().long_division(&d);
    assert_eq!(q.coeffs.len(), 3);
    for &x in &[-2.0, -0.5, 0.0, 1.3, 3.0] {
        let lhs = p.eval(x);
        let rhs = q.eval(x) * d.eval(x) + r.eval(x);
        assert_approx_eq!(lhs, rhs, 1e-9);
    }
}

#[test]
fn long_division_by_higher_degree_gives_zero_quotient() {
    let p = Polynomial::new(vec![1.0, 2.0]);
    let d = Polynomial::new(vec![1.0, 0.0, 1.0]);
    let (q, r) = p.long_division(&d);
    assert_eq!(q.coeffs, vec![0.0]);
    assert_eq!(r.coeffs, vec![1.0, 2.0]);
}
