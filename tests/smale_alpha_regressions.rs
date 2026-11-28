//! Analytic checks for Smale derivative scaling and factor reuse.

use rssn::kernels::homotopy_cert::smale_alpha;
use rssn::kernels::rootfind::PolyTerm;

fn term(
    coef: f64,
    exps: &[u32],
) -> PolyTerm {
    PolyTerm {
        coef,
        exps: exps.to_vec(),
    }
}

#[test]
fn factorial_overflow_does_not_erase_a_high_derivative() {
    // At x=0 the Jacobian is 1 and the only higher Taylor coefficient is
    // the coefficient of x^171. Taking its 170th root gives gamma.
    for (coef, gamma) in [(1.0, 1.0), (2.0_f64.powi(340), 4.0)] {
        let system = vec![vec![term(coef, &[171]), term(1.0, &[1]), term(-1.0, &[0])]];
        let cert = smale_alpha(&system, &[(0.0, 0.0)]).unwrap();
        assert_eq!(cert.beta, 1.0);
        assert!((cert.gamma - gamma).abs() < 2e-15);
        assert!((cert.alpha - gamma).abs() < 2e-15);
        assert!(!cert.certified && cert.radius.is_infinite());
    }
    // At x=1 the largest normalized derivative of x^d occurs at k=2,
    // giving gamma=(d-1)/2 and beta=1/d.
    let cert = smale_alpha(&[vec![term(1.0, &[171])]], &[(1.0, 0.0)]).unwrap();
    assert!((cert.beta - 1.0 / 171.0).abs() < 1e-17);
    assert!((cert.gamma - 85.0).abs() < 2e-13);
    assert!(!cert.certified);

    // The mixed derivative x^170 y occurs in 171 ordered columns.
    // A coefficient of 171 gives a normalized derivative of one in each.
    let mixed = vec![
        vec![
            term(171.0, &[170, 1]),
            term(1.0, &[1, 0]),
            term(-0.5, &[0, 0]),
        ],
        vec![term(1.0, &[0, 1])],
    ];
    let cert = smale_alpha(&mixed, &[(0.0, 0.0); 2]).unwrap();
    assert_eq!(cert.beta, 0.5);
    assert!((cert.gamma - 171.0_f64.powf(1.0 / 340.0)).abs() < 3e-15);
    assert!(!cert.certified);
}

#[test]
fn tiny_and_large_steps_have_representable_norms_and_radii() {
    for magnitude in [2.0_f64.powi(-600), 2.0_f64.powi(600)] {
        let linear = vec![vec![term(1.0, &[1]), term(-magnitude, &[0])]];
        let cert = smale_alpha(&linear, &[(0.0, 0.0)]).unwrap();
        assert!(cert.certified);
        assert_eq!(cert.beta, magnitude);
        assert_eq!(cert.radius, 2.0 * magnitude);
        assert_eq!((cert.gamma, cert.alpha), (0.0, 0.0));

        let scaled = vec![vec![term(magnitude, &[1]), term(-0.5 * magnitude, &[0])]];
        let cert = smale_alpha(&scaled, &[(0.0, 0.0)]).unwrap();
        assert!(cert.certified);
        assert_eq!(cert.beta, 0.5);
        assert_eq!(cert.radius, 1.0);
    }
    let tiny = 2.0_f64.powi(-600);
    let quadratic = vec![vec![term(tiny, &[2]), term(1.0, &[1]), term(-0.25, &[0])]];
    let cert = smale_alpha(&quadratic, &[(0.0, 0.0)]).unwrap();
    assert_eq!(cert.gamma, tiny);
    assert_eq!(cert.alpha, 0.25 * tiny);

    let outside = vec![vec![term(tiny.recip(), &[1]), term(-tiny, &[0])]];
    assert!(smale_alpha(&outside, &[(0.0, 0.0)]).is_none());
}

#[test]
fn pivoted_jacobian_solves_the_step_and_mixed_derivative_columns() {
    // J swaps the coordinates. The Taylor tensor has one x^2 entry of
    // size 2 and two mixed entries of size 1/2 on the other component.
    let system = vec![
        vec![
            term(1.0, &[0, 1]),
            term(1.0, &[1, 1]),
            term(-0.125, &[0, 0]),
        ],
        vec![term(1.0, &[1, 0]), term(2.0, &[2, 0]), term(-0.25, &[0, 0])],
    ];
    let cert = smale_alpha(&system, &[(0.0, 0.0); 2]).unwrap();
    let beta = 0.125_f64.hypot(0.25);
    let gamma = 4.5_f64.sqrt();
    assert!((cert.beta - beta).abs() < 2e-16);
    assert!((cert.gamma - gamma).abs() < 2e-15);
    assert!((cert.alpha - beta * gamma).abs() < 2e-15);
    assert!(!cert.certified);
}

#[test]
fn malformed_degree_and_nonfinite_input_return_no_certificate() {
    let overflow = vec![
        vec![term(1.0, &[u32::MAX, 2]), term(1.0, &[1, 0])],
        vec![term(1.0, &[0, 1])],
    ];
    assert!(smale_alpha(&overflow, &[(1.0, 0.0); 2]).is_none());
    assert!(smale_alpha(&[vec![term(1.0, &[1])]], &[(f64::INFINITY, 0.0)]).is_none());
    assert!(smale_alpha(&[vec![term(f64::NAN, &[1])]], &[(1.0, 0.0)]).is_none());
    assert!(smale_alpha(&[vec![term(1.0, &[2])]], &[(0.0, 0.0)]).is_none());
}
