//! Boundary, multiple-root and large-interval tests for numerical Sturm search.

use rssn::kernels::polynomial::Polynomial;
use rssn::kernels::real_roots::find_roots;
use rssn::kernels::real_roots::find_roots_in;
use rssn::kernels::real_roots::isolate_real_roots;
use rssn::kernels::real_roots::refine_root_bisection;

#[test]
fn split_root_at_zero_does_not_replace_the_next_root() {
    for coefficients in [vec![1.0, -0.25, 0.0], vec![1.0, -0.25, 0.0, 0.0]] {
        let polynomial = Polynomial::new(coefficients);
        let roots = find_roots(&polynomial, 1e-13).unwrap();
        assert_eq!(roots.len(), 2);
        assert_eq!(roots[0], 0.0);
        assert!((roots[1] - 0.25).abs() < 1e-13);

        let intervals = isolate_real_roots(&polynomial, 1e-13).unwrap();
        assert_eq!(intervals.len(), 2);
        assert!(intervals[0].0 <= 0.0 && intervals[0].1 == 0.0);
        assert!(intervals[1].0 > 0.0 && intervals[1].0 <= 0.25 && intervals[1].1 >= 0.25);
        assert_eq!(refine_root_bisection(&polynomial, intervals[0], 1e-13), 0.0);
        let positive = refine_root_bisection(&polynomial, intervals[1], 1e-13);
        assert!((positive - 0.25).abs() < 1e-13);
    }
}

#[test]
fn nonzero_split_root_stays_in_one_closed_interval() {
    // The Cauchy bound is 32, so isolation splits at -3 exactly. Dividing
    // the coefficients by 31 makes floating Horner nonzero at this split.
    let polynomial = Polynomial::new(vec![1.0, 10.0, 31.0, 30.0]);
    let roots = [-5.0, -3.0, -2.0];
    let intervals = isolate_real_roots(&polynomial, 1e-12).unwrap();
    assert_eq!(intervals.len(), roots.len());
    for (&expected, &(lo, hi)) in roots.iter().zip(&intervals) {
        assert!((lo..=hi).contains(&expected));
        assert_eq!(
            roots
                .iter()
                .filter(|root| (lo..=hi).contains(*root))
                .count(),
            1
        );
        assert!((refine_root_bisection(&polynomial, (lo, hi), 1e-12) - expected).abs() < 1e-12);
    }
}

#[test]
fn bounded_sturm_search_includes_both_endpoints_and_singletons() {
    let polynomial = Polynomial::new(vec![1.0, 0.0, -1.0, 0.0]);
    assert_eq!(
        find_roots_in(&polynomial, (-1.0, 1.0), 1e-12).unwrap(),
        [-1.0, 0.0, 1.0]
    );
    assert_eq!(
        find_roots_in(&polynomial, (0.0, 0.0), 1e-12).unwrap(),
        [0.0]
    );
    assert_eq!(
        find_roots_in(&polynomial, (1.0, 2.0), 1e-12).unwrap(),
        [1.0]
    );
    assert!(
        find_roots_in(&polynomial, (0.1, 0.9), 1e-12)
            .unwrap()
            .is_empty()
    );
    assert!(
        find_roots_in(&polynomial, (0.5, 0.5), 1e-12)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn bounded_search_uses_original_coefficients_for_endpoint_roots() {
    // Scaling x(x - 1/4)(x - 1) by its largest coefficient perturbs the
    // polynomial's exact zero at one.
    let polynomial = Polynomial::new(vec![1.0, -1.25, 0.25, 0.0]);
    assert_eq!(
        find_roots_in(&polynomial, (1.0, 1.0), 1e-12).unwrap(),
        [1.0]
    );
    assert_eq!(
        find_roots_in(&polynomial, (0.5, 1.0), 1e-12).unwrap(),
        [1.0]
    );

    // Even Horner with FMA loses the linear term during cancellation.
    let large = 9_007_199_254_740_992.0;
    let polynomial = Polynomial::new(vec![large, 1.0, -large, -1.0]);
    assert_eq!(
        find_roots_in(&polynomial, (1.0, 1.0), 1e-12).unwrap(),
        [1.0]
    );
    assert_eq!(refine_root_bisection(&polynomial, (1.0, 1.0), 1e-12), 1.0);

    let polynomial = Polynomial::new(vec![f64::from_bits(1), -f64::MAX, 0.0]);
    assert_eq!(
        find_roots_in(&polynomial, (0.0, 0.0), 1e-12).unwrap(),
        [0.0]
    );
    assert_eq!(refine_root_bisection(&polynomial, (0.0, 0.0), 1e-12), 0.0);
}

#[test]
fn bounded_search_avoids_an_unrepresentable_global_root_bound() {
    let leading = f64::from_bits(1_u64 << 49);
    let polynomial = Polynomial::new(vec![leading, 0.0, -1.0]);
    let tolerance = 1e141;
    assert!(find_roots(&polynomial, tolerance).is_err());
    let roots = find_roots_in(&polynomial, (0.0, f64::MAX), tolerance).unwrap();
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 1.0 / leading.sqrt()).abs() < tolerance);
}

#[test]
fn global_search_includes_roots_at_the_float_limit() {
    for root in [-f64::MAX, f64::MAX] {
        let polynomial = Polynomial::new(vec![1.0, -root]);
        assert_eq!(find_roots(&polynomial, 1e292).unwrap(), [root]);
        assert_eq!(isolate_real_roots(&polynomial, 1e292).unwrap().len(), 1);
    }
}

#[test]
fn bounded_search_refines_an_even_root_across_a_huge_interval() {
    // (x^2 - 2)^2 has two distinct roots and never changes residual sign.
    let polynomial = Polynomial::new(vec![1.0, 0.0, -4.0, 0.0, 4.0]);
    let tolerance = 1e-12;
    let roots = find_roots_in(&polynomial, (0.0, f64::MAX), tolerance).unwrap();
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 2.0_f64.sqrt()).abs() < tolerance);
    let global = find_roots(&polynomial, tolerance).unwrap();
    assert_eq!(global.len(), 2);
    assert!((global[0] + 2.0_f64.sqrt()).abs() < tolerance);
}

#[test]
fn sturm_search_reports_invalid_coefficients_or_bounds() {
    let valid = Polynomial::new(vec![1.0, -1.0]);
    assert!(find_roots_in(&valid, (2.0, 1.0), 1e-6).is_err());
    assert!(find_roots_in(&valid, (0.0, f64::INFINITY), 1e-6).is_err());
    assert!(find_roots_in(&valid, (0.0, 2.0), f64::NAN).is_err());
    for coefficients in [
        vec![0.0, 1.0],
        vec![1.0, f64::NAN],
        vec![1.0, f64::INFINITY],
    ] {
        let invalid = Polynomial::new(coefficients);
        assert!(find_roots(&invalid, 1e-6).is_err());
        assert!(isolate_real_roots(&invalid, 1e-6).is_err());
        assert!(find_roots_in(&invalid, (0.0, 2.0), 1e-6).is_err());
        assert!(refine_root_bisection(&invalid, (0.0, 2.0), 1e-6).is_nan());
    }
}
