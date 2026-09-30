//! Real-root isolation via Sturm sequences (ported from `numerical_real_roots_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::polynomial::Polynomial;
use rssn::kernels::real_roots::{find_roots, isolate_real_roots, refine_root_bisection, sturm_sequence};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn quadratic_roots() {
    // x^2 - 4
    let poly = Polynomial::new(vec![1.0, 0.0, -4.0]);
    let roots = find_roots(&poly, 1e-9).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(roots.len(), 2);
    assert_approx_eq!(roots[0], -2.0);
    assert_approx_eq!(roots[1], 2.0);
}

#[test]
fn cubic_roots() {
    // (x-1)(x-2)(x-3)
    let poly = Polynomial::new(vec![1.0, -6.0, 11.0, -6.0]);
    let roots = find_roots(&poly, 1e-9).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(roots.len(), 3);
    for (r, w) in roots.iter().zip([1.0, 2.0, 3.0]) {
        assert_approx_eq!(*r, w);
    }
}

#[test]
fn close_roots_are_separated() {
    // (x - 1)(x - 1.001)
    let poly = Polynomial::new(vec![1.0, -2.001, 1.001]);
    let roots = find_roots(&poly, 1e-9).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(roots.len(), 2);
    assert_approx_eq!(roots[0], 1.0, 1e-6);
    assert_approx_eq!(roots[1], 1.001, 1e-6);
}

#[test]
fn no_real_roots() {
    let poly = Polynomial::new(vec![1.0, 0.0, 1.0]);
    let roots = find_roots(&poly, 1e-9).unwrap_or_else(|e| panic!("{e}"));
    assert!(roots.is_empty(), "x^2 + 1 has no real roots, got {roots:?}");
}

#[test]
fn sturm_sequence_of_cubic_starts_with_poly_and_derivative() {
    let poly = Polynomial::new(vec![1.0, -6.0, 11.0, -6.0]);
    let seq = sturm_sequence(&poly);
    assert!(seq.len() >= 3, "a cubic with 3 simple roots has a full chain, got {}", seq.len());
    assert_eq!(seq[0].coeffs, poly.coeffs);
    assert_eq!(seq[1].coeffs, poly.derivative().coeffs);
    assert!(sturm_sequence(&Polynomial::new(vec![])).is_empty());
}

#[test]
fn isolation_returns_one_interval_per_root_each_containing_it() {
    let poly = Polynomial::new(vec![1.0, -6.0, 11.0, -6.0]);
    let mut intervals = isolate_real_roots(&poly, 1e-6).unwrap_or_else(|e| panic!("{e}"));
    intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert_eq!(intervals.len(), 3);
    for ((lo, hi), root) in intervals.iter().zip([1.0, 2.0, 3.0]) {
        assert!(*lo <= root + 1e-6 && root - 1e-6 <= *hi, "[{lo}, {hi}] misses {root}");
    }
}

#[test]
fn bisection_refines_root() {
    let poly = Polynomial::new(vec![1.0, 0.0, -2.0]);
    let r = refine_root_bisection(&poly, (1.0, 2.0), 1e-12);
    assert_approx_eq!(r, 2f64.sqrt(), 1e-10);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_quadratic_roots_match_formula(a in 0.5..10.0f64, b in -10.0..10.0f64, c in -10.0..10.0f64) {
        let disc = b * b - 4.0 * a * c;
        prop_assume!(disc.abs() > 1e-3);
        let poly = Polynomial::new(vec![a, b, c]);
        let roots = find_roots(&poly, 1e-9).map_err(TestCaseError::fail)?;
        if disc > 0.0 {
            prop_assert_eq!(roots.len(), 2, "roots {:?}", roots);
            let lo = (-b - disc.sqrt()) / (2.0 * a);
            let hi = (-b + disc.sqrt()) / (2.0 * a);
            prop_assert!((roots[0] - lo).abs() < 1e-5 && (roots[1] - hi).abs() < 1e-5, "{roots:?} vs {lo}, {hi}");
            for r in roots {
                prop_assert!(poly.eval(r).abs() < 1e-4);
            }
        } else {
            prop_assert!(roots.is_empty(), "expected no roots, got {:?}", roots);
        }
    }

    #[test]
    fn prop_planted_cubic_roots_are_found(r1 in -5.0..-1.0f64, r2 in -0.5..0.5f64, r3 in 1.0..5.0f64) {
        // (x - r1)(x - r2)(x - r3)
        let s1 = r1 + r2 + r3;
        let s2 = r1 * r2 + r1 * r3 + r2 * r3;
        let s3 = r1 * r2 * r3;
        let poly = Polynomial::new(vec![1.0, -s1, s2, -s3]);
        let roots = find_roots(&poly, 1e-10).map_err(TestCaseError::fail)?;
        prop_assert_eq!(roots.len(), 3);
        for (got, want) in roots.iter().zip([r1, r2, r3]) {
            prop_assert!((got - want).abs() < 1e-6, "{got} vs {want}");
        }
    }
}
