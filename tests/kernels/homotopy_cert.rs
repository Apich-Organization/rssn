//! Tests for Smale alpha-theory certificates.

use rssn::kernels::homotopy::{HomotopyOptions, Precision};
use rssn::kernels::homotopy_cert::{ALPHA_0, certify_endpoints, homotopy_solve_certified, smale_alpha};
use rssn::kernels::homotopy::homotopy_solve;
use rssn::kernels::rootfind::PolyTerm;

fn t(coef: f64, exps: &[u32]) -> PolyTerm {
    PolyTerm { coef, exps: exps.to_vec() }
}

#[test]
fn alpha_threshold_constant() {
    let exact = (13.0 - 3.0 * 17.0_f64.sqrt()) / 4.0;
    assert!((ALPHA_0 - exact).abs() < 1e-15);
}

#[test]
fn invariants_of_a_scalar_quadratic() {
    // f = x^2 - 2: beta = |f/f'|, gamma = |f''/2 / f'| = 1 / |2x|
    let sys = vec![vec![t(1.0, &[2]), t(-2.0, &[0])]];
    let c = smale_alpha(&sys, &[(1.5, 0.0)]).unwrap();
    assert!((c.beta - 0.25 / 3.0).abs() < 1e-14);
    assert!((c.gamma - 1.0 / 3.0).abs() < 1e-14);
    assert!((c.alpha - 0.25 / 9.0).abs() < 1e-14);
    assert!(c.certified);
    // the associated zero lies within 2 beta
    assert!((1.5 - 2.0_f64.sqrt()).abs() <= c.radius);
    // the test is sufficient, not necessary: x = 1 converges under Newton but has alpha = 1/4
    let c = smale_alpha(&sys, &[(1.0, 0.0)]).unwrap();
    assert!((c.alpha - 0.25).abs() < 1e-14 && !c.certified && c.radius.is_infinite());
    // near a double root the point is not certified
    let dbl = vec![vec![t(1.0, &[2]), t(-2.0, &[1]), t(1.0, &[0])]];
    let c = smale_alpha(&dbl, &[(1.1, 0.0)]).unwrap();
    assert!(!c.certified, "alpha {}", c.alpha);
    // singular Jacobian
    assert!(smale_alpha(&dbl, &[(1.0, 0.0)]).is_none());
    assert!(smale_alpha(&dbl, &[(1.0, 0.0), (2.0, 0.0)]).is_none());
}

#[test]
fn higher_degree_invariants_match_a_direct_computation() {
    // f = x^3 - 1 at x = 1 + 0.1 i: f' = 3x^2, f''/2 = 3x, f'''/6 = 1
    let sys = vec![vec![t(1.0, &[3]), t(-1.0, &[0])]];
    let x = num_complex::Complex64::new(1.0, 0.1);
    let c = smale_alpha(&sys, &[(x.re, x.im)]).unwrap();
    let fp = 3.0 * x * x;
    let beta = ((x * x * x - 1.0) / fp).norm();
    let g2 = (3.0 * x / fp).norm();
    let g3 = (1.0 / fp).norm().sqrt();
    assert!((c.beta - beta).abs() < 1e-13);
    assert!((c.gamma - g2.max(g3)).abs() < 1e-13, "{} vs {}", c.gamma, g2.max(g3));
    // two variables: x^2 + y^2 - 1, x y - 1/4 at an approximate solution (the Frobenius bound is an upper bound)
    let sys = vec![vec![t(1.0, &[2, 0]), t(1.0, &[0, 2]), t(-1.0, &[0, 0])], vec![t(1.0, &[1, 1]), t(-0.25, &[0, 0])]];
    let c = smale_alpha(&sys, &[(0.97, 0.0), (0.26, 0.0)]).unwrap();
    assert!(c.certified && c.beta < 0.05, "{c:?}");
}

#[test]
fn certified_complete_solution_lists() {
    // cyclic 3-roots: six simple solutions, all certified and separated, so the list is complete
    let sys = vec![
        vec![t(1.0, &[1, 0, 0]), t(1.0, &[0, 1, 0]), t(1.0, &[0, 0, 1])],
        vec![t(1.0, &[1, 1, 0]), t(1.0, &[0, 1, 1]), t(1.0, &[1, 0, 1])],
        vec![t(1.0, &[1, 1, 1]), t(-1.0, &[0, 0, 0])],
    ];
    let r = homotopy_solve_certified(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.report.paths, 6);
    assert_eq!(r.report.certified, 6);
    assert_eq!(r.report.distinct, 6);
    assert!(r.report.complete);
    for c in r.report.certificates.iter().flatten() {
        assert!(c.certified && c.alpha < 1e-6, "{c:?}");
    }
    // a system with solutions at infinity is not complete in the certified sense
    let sys = vec![vec![t(1.0, &[2, 0]), t(-1.0, &[0, 0])], vec![t(1.0, &[1, 1]), t(-1.0, &[0, 0])]];
    let r = homotopy_solve_certified(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.report.certified, 2);
    assert_eq!(r.report.distinct, 2);
    assert!(!r.report.complete);
    // duplicated endpoints are not "distinct": feed the same endpoint twice
    let mut ends = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    let dup = ends.iter().find(|p| !p.at_infinity).unwrap().clone();
    ends.push(dup);
    let rep = certify_endpoints(&sys, &ends, 1.0);
    assert!(rep.distinct <= 1, "duplicates must not be counted as distinct: {}", rep.distinct);
}

#[test]
fn certification_of_tight_clusters_in_double_double() {
    let d = 2.0_f64.powi(-27);
    let sys = vec![vec![t(1.0, &[3]), t(-4.0 * d, &[2]), t(3.0 * d * d, &[1]), t(0.0, &[0])]];
    let opts = HomotopyOptions { precision: Precision::DoubleDouble, ..HomotopyOptions::default() };
    let r = homotopy_solve_certified(&sys, &opts).unwrap();
    assert!(r.report.complete, "{:?}", r.report);
    assert_eq!(r.report.distinct, 3);
    // the default entry point escalates to extended precision by itself
    let r = homotopy_solve_certified(&sys, &HomotopyOptions::default()).unwrap();
    assert!(r.report.complete, "{:?}", r.report);
    assert!(r.ends.iter().any(|p| p.extended));
}
