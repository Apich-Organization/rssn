use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;

#[test]

fn test_find_roots_quadratic() {

    // x^2 - 4 = 0 => roots -2, 2
    let poly =
        numerical_Polynomial::new(
            vec![1.0, 0.0, -4.0],
        );

    let roots = numerical_find_roots(
        &poly, 1e-9,
    )
    .unwrap();

    assert_eq!(roots.len(), 2);

    assert_approx_eq!(roots[0], -2.0);

    assert_approx_eq!(roots[1], 2.0);
}

#[test]

fn test_find_roots_cubic() {

    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    let poly =
        numerical_Polynomial::new(
            vec![
                1.0, -6.0, 11.0, -6.0,
            ],
        );

    let roots = numerical_find_roots(
        &poly, 1e-9,
    )
    .unwrap();

    assert_eq!(roots.len(), 3);

    assert_approx_eq!(roots[0], 1.0);

    assert_approx_eq!(roots[1], 2.0);

    assert_approx_eq!(roots[2], 3.0);
}

#[test]

fn test_find_roots_close() {

    // Roots at 1.0 and 1.001
    // (x-1)(x-1.001) = x^2 - 2.001x + 1.001
    let poly =
        numerical_Polynomial::new(
            vec![1.0, -2.001, 1.001],
        );

    let roots = numerical_find_roots(
        &poly, 1e-9,
    )
    .unwrap();

    assert_eq!(roots.len(), 2);

    assert_approx_eq!(
        roots[0], 1.0, 1e-6
    );

    assert_approx_eq!(
        roots[1], 1.001, 1e-6
    );
}

#[cfg(test)]

mod proptests {

    use super::*;

    proptest! {
        #[test]
        fn prop_find_roots_quadratic(a in -10.0..10.0f64, b in -10.0..10.0f64, c in -10.0..10.0f64) {
            // ax^2 + bx + c = 0
            // Roots: (-b +/- sqrt(b^2 - 4ac)) / 2a
            if a.abs() > 1e-3 {
                let disc = b*b - 4.0*a*c;
                let poly = numerical_Polynomial::new(vec![a, b, c]);
                let roots_res = numerical_find_roots(&poly, 1e-9);

                if disc > 1e-6 {
                    // Two real roots
                    if let Ok(roots) = roots_res {
                         assert!(roots.len() == 2, "Expected 2 roots for disc included, got {}", roots.len());
                         // Verify roots
                         for r in roots {
                             let val = poly.eval(r);
                             // Quadratic evaluation can be numerically unstable for large coeffs, but here reasonable range
                             assert!(val.abs() < 1e-4, "Root {}, val {}", r, val);
                         }
                    }
                } else if disc < -1e-6 {
                    // No real roots
                    if let Ok(roots) = roots_res {
                        assert!(roots.is_empty(), "Expected 0 roots, got {}", roots.len());
                    }
                }
            }
        }
    }
}
