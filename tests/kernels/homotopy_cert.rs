//! Tests for Smale alpha-theory certificates and the rigorous
//! Krawczyk certificates.

use rssn::kernels::homotopy::HomotopyOptions;
use rssn::kernels::homotopy::Precision;
use rssn::kernels::homotopy::homotopy_solve;
use rssn::kernels::homotopy_cert::ALPHA_0;
use rssn::kernels::homotopy_cert::certify_endpoints;
use rssn::kernels::homotopy_cert::certify_endpoints_rigorous;
use rssn::kernels::homotopy_cert::certify_homotopy_node;
use rssn::kernels::homotopy_cert::certify_zero;
use rssn::kernels::homotopy_cert::disjoint;
use rssn::kernels::homotopy_cert::homotopy_solve_certified;
use rssn::kernels::homotopy_cert::homotopy_solve_verified;
use rssn::kernels::homotopy_cert::smale_alpha;
use rssn::kernels::homotopy_cert::verify_zero;
use rssn::kernels::rootfind::PolyTerm;

fn t(
    coef: f64,
    exps: &[u32],
) -> PolyTerm {
    PolyTerm {
        coef,
        exps: exps.to_vec(),
    }
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
    assert!(
        (c.gamma - g2.max(g3)).abs() < 1e-13,
        "{} vs {}",
        c.gamma,
        g2.max(g3)
    );
    // two variables: x^2 + y^2 - 1, x y - 1/4 at an approximate solution (the Frobenius bound is an upper bound)
    let sys = vec![
        vec![t(1.0, &[2, 0]), t(1.0, &[0, 2]), t(-1.0, &[0, 0])],
        vec![t(1.0, &[1, 1]), t(-0.25, &[0, 0])],
    ];
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
    let sys = vec![
        vec![t(1.0, &[2, 0]), t(-1.0, &[0, 0])],
        vec![t(1.0, &[1, 1]), t(-1.0, &[0, 0])],
    ];
    let r = homotopy_solve_certified(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.report.certified, 2);
    assert_eq!(r.report.distinct, 2);
    assert!(!r.report.complete);
    // duplicated endpoints are not "distinct": feed the same endpoint twice
    let mut ends = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    let dup = ends.iter().find(|p| !p.at_infinity).unwrap().clone();
    ends.push(dup);
    let rep = certify_endpoints(&sys, &ends, 1.0);
    assert!(
        rep.distinct <= 1,
        "duplicates must not be counted as distinct: {}",
        rep.distinct
    );
}

#[test]
fn certification_of_tight_clusters_in_double_double() {
    let d = 2.0_f64.powi(-27);
    let sys = vec![vec![
        t(1.0, &[3]),
        t(-4.0 * d, &[2]),
        t(3.0 * d * d, &[1]),
        t(0.0, &[0]),
    ]];
    let opts = HomotopyOptions {
        precision: Precision::DoubleDouble,
        ..HomotopyOptions::default()
    };
    let r = homotopy_solve_certified(&sys, &opts).unwrap();
    assert!(r.report.complete, "{:?}", r.report);
    assert_eq!(r.report.distinct, 3);
    // the default entry point escalates to extended precision by itself
    let r = homotopy_solve_certified(&sys, &HomotopyOptions::default()).unwrap();
    assert!(r.report.complete, "{:?}", r.report);
    assert!(r.ends.iter().any(|p| p.extended));
}

// ----------------------------------------------------------------------
// Rigorous (interval) certificates
// ----------------------------------------------------------------------

#[test]
fn krawczyk_proves_a_simple_zero() {
    // x^2 - 2 at the double-precision root: a box around it is certified,
    // and the box contains the exact root.
    let sys = vec![vec![t(1.0, &[2]), t(-2.0, &[0])]];
    let root = 2.0_f64.sqrt();
    let ball = certify_zero(&sys, &[(root, 0.0)]).unwrap();
    assert!(ball.contains(root, 0.0), "{ball:?}");
    assert!(ball.radius < 1e-12, "radius {}", ball.radius);
    // a modestly over-wide initial box is contracted to the root's ball
    let wide = verify_zero(&sys, &[(root, 0.0)], 1e-6).unwrap();
    assert!(wide.contains(root, 0.0), "{wide:?}");
    assert!(wide.radius < 1e-13, "contracted radius {}", wide.radius);
    // a point with no root nearby cannot be certified
    assert!(certify_zero(&sys, &[(0.0, 0.0)]).is_none());
    // a double root has a singular Jacobian: no uniqueness certificate
    let double = vec![vec![t(1.0, &[2]), t(-2.0, &[1]), t(1.0, &[0])]];
    assert!(certify_zero(&double, &[(1.0, 0.0)]).is_none());
    assert!(certify_zero(&double, &[(1.05, 0.0)]).is_none());
    // a linear system: the box is a single point up to rounding
    let linear = vec![vec![t(2.0, &[1]), t(-1.0, &[0])]];
    let ball = certify_zero(&linear, &[(0.5, 0.0)]).unwrap();
    assert!(ball.contains(0.5, 0.0) && ball.radius < 1e-13);
}

#[test]
fn rigorous_endpoint_certification_of_a_system() {
    // x^2 + y^2 - 1, x y - 1/4: four simple solutions.
    let sys = vec![
        vec![t(1.0, &[2, 0]), t(1.0, &[0, 2]), t(-1.0, &[0, 0])],
        vec![t(1.0, &[1, 1]), t(-0.25, &[0, 0])],
    ];
    let ends = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    let report = certify_endpoints_rigorous(&sys, &ends);
    assert_eq!(report.paths, 4);
    assert_eq!(report.certified, 4, "{report:?}");
    assert_eq!(report.distinct, 4);
    assert!(report.complete);
    // every certified ball is centred on a point that nearly solves the system
    for ball in report.certificates.iter().flatten() {
        let z = [
            num_complex::Complex64::new(ball.centre[0].0, ball.centre[0].1),
            num_complex::Complex64::new(ball.centre[1].0, ball.centre[1].1),
        ];
        let f0 = z[0] * z[0] + z[1] * z[1] - 1.0;
        let f1 = z[0] * z[1] - 0.25;
        assert!(f0.norm() < 1e-8 && f1.norm() < 1e-8, "{f0} {f1}");
    }
    // two balls around different roots are disjoint and the predicate agrees
    let balls: Vec<_> = report.certificates.iter().flatten().cloned().collect();
    for (i, a) in balls.iter().enumerate() {
        for (j, b) in balls.iter().enumerate() {
            if i != j {
                assert!(disjoint(a, b));
            }
        }
    }
    // a duplicated endpoint is not counted as distinct: the two copies are
    // the same zero, so neither is separated from all the others and only
    // the remaining three zeros count
    let mut duplicated = ends.clone();
    duplicated.push(ends[0].clone());
    let report = certify_endpoints_rigorous(&sys, &duplicated);
    assert_eq!(report.certified, 5);
    assert_eq!(report.distinct, 3);
    assert!(!report.complete);
}

#[test]
fn verified_solve_is_complete_and_certifies_steps() {
    let sys = vec![
        vec![t(1.0, &[1, 0, 0]), t(1.0, &[0, 1, 0]), t(1.0, &[0, 0, 1])],
        vec![t(1.0, &[1, 1, 0]), t(1.0, &[0, 1, 1]), t(1.0, &[1, 0, 1])],
        vec![t(1.0, &[1, 1, 1]), t(-1.0, &[0, 0, 0])],
    ];
    let result = homotopy_solve_verified(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(result.rigorous.paths, 6);
    assert_eq!(result.rigorous.distinct, 6);
    assert!(result.rigorous.complete);
    assert_eq!(result.floating.certified, 6);
    // every regular path's nodes were interval-verified
    for end in &result.ends {
        assert!(end.verified_steps > 0, "{end:?}");
        assert_eq!(end.unverified_steps, 0, "{end:?}");
        assert!(end.steps_certified, "{end:?}");
    }
}

#[test]
fn solving_without_step_certification_reports_no_verified_steps() {
    let sys = vec![vec![t(1.0, &[2]), t(-2.0, &[0])]];
    let ends = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    assert!(
        ends.iter()
            .all(|e| e.verified_steps == 0 && !e.steps_certified)
    );
}

#[test]
fn homotopy_node_certificate() {
    // x^2 - 2 homogenised to degree 2: X1^2 - 2 X0^2, patch X1 = 1.
    let homog = vec![vec![t(1.0, &[0, 2]), t(-2.0, &[2, 0])]];
    let root = 2.0_f64.sqrt();
    let node = [(1.0 / root, 0.0), (1.0, 0.0)];
    let ball = certify_homotopy_node(
        &homog,
        &[2],
        num_complex::Complex64::new(1.0, 0.0),
        &[
            num_complex::Complex64::new(0.0, 0.0),
            num_complex::Complex64::new(1.0, 0.0),
        ],
        &node,
        (1.0, 0.0),
    );
    assert!(
        ball.is_some(),
        "the patch point of a simple root is certified"
    );
    // a point off the patch is not a zero of the homogenised system
    let off = [(0.5, 0.0), (1.0, 0.0)];
    assert!(
        certify_homotopy_node(
            &homog,
            &[2],
            num_complex::Complex64::new(1.0, 0.0),
            &[
                num_complex::Complex64::new(0.0, 0.0),
                num_complex::Complex64::new(1.0, 0.0),
            ],
            &off,
            (1.0, 0.0)
        )
        .is_none()
    );
}

#[test]
fn extended_precision_and_rigorous_certificates_agree() {
    let d = 2.0_f64.powi(-27);
    let sys = vec![vec![
        t(1.0, &[3]),
        t(-4.0 * d, &[2]),
        t(3.0 * d * d, &[1]),
        t(0.0, &[0]),
    ]];
    let opts = HomotopyOptions {
        precision: Precision::DoubleDouble,
        ..HomotopyOptions::default()
    };
    let result = homotopy_solve_verified(&sys, &opts).unwrap();
    assert!(result.rigorous.complete, "{:?}", result.rigorous);
    // the roots are 0, d, 3d: the certified boxes separate them
    let balls: Vec<_> = result
        .rigorous
        .certificates
        .iter()
        .flatten()
        .cloned()
        .collect();
    assert_eq!(balls.len(), 3);
    assert!(balls.iter().any(|b| b.contains(0.0, 0.0)));
    assert!(balls.iter().any(|b| b.contains(d, 0.0)));
    assert!(balls.iter().any(|b| b.contains(3.0 * d, 0.0)));
    for i in 0..balls.len() {
        for j in i + 1..balls.len() {
            assert!(disjoint(&balls[i], &balls[j]));
        }
    }
}
