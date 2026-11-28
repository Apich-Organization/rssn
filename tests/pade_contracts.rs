//! Native Padé entry points agree on coefficients and retain their error contracts.

use rssn::kernels::approx::ApproxError;
use rssn::kernels::approx::{
    self,
};
use rssn::kernels::interp_ext::InterpError;
use rssn::kernels::interp_ext::{
    self,
};

#[test]
fn both_entry_points_reconstruct_the_taylor_polynomial() {
    let coeffs = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0];
    let (p, q) = approx::pade(&coeffs, 2, 2).unwrap();
    assert_eq!(
        (p.clone(), q.clone()),
        interp_ext::pade(&coeffs, 2, 2).unwrap()
    );
    assert_eq!(q[0], 1.0);
    for (k, expected) in p.iter().copied().map(Some).chain([None, None]).enumerate() {
        let product = q
            .iter()
            .zip(coeffs[..=k].iter().rev())
            .map(|(&denominator, &coefficient)| denominator * coefficient)
            .sum::<f64>();
        assert!((product - expected.unwrap_or(0.0)).abs() < 1.0e-15);
    }

    // Here the denominator has larger degree than the numerator, so the
    // triangular system includes zero coefficients at negative indices.
    let coeffs = [1.0, 0.5, 0.25, 0.125];
    let expected = (vec![1.0], vec![1.0, -0.5, 0.0, 0.0]);
    assert_eq!(approx::pade(&coeffs, 0, 3).unwrap(), expected);
    assert_eq!(interp_ext::pade(&coeffs, 0, 3).unwrap(), expected);
}

#[test]
fn overflow_and_singular_inputs_keep_their_own_error_channels() {
    for (m, n) in [(usize::MAX, 0), (1, usize::MAX), (0, usize::MAX)] {
        assert_eq!(approx::pade(&[1.0], m, n), Err(ApproxError::Invalid));
        let Err(InterpError::Domain(detail)) = interp_ext::pade(&[1.0], m, n) else {
            panic!("expected a coefficient-count error");
        };
        assert!(detail.contains("coefficient count exceeds usize"));
    }
    assert_eq!(
        interp_ext::pade(&[1.0], 2, 2),
        Err(InterpError::Domain(
            "Padé [2/2] needs at least 5 coefficients, got 1".into()
        )),
    );
    assert_eq!(approx::pade(&[0.0; 3], 1, 1), Err(ApproxError::Singular));
    assert_eq!(
        interp_ext::pade(&[0.0; 3], 1, 1),
        Err(InterpError::Singular)
    );
    assert_eq!(
        approx::pade(&[1.0, 1.0, f64::NAN], 1, 1),
        Err(ApproxError::Singular)
    );
    assert_eq!(
        interp_ext::pade(&[1.0, 1.0, f64::NAN], 1, 1),
        Err(InterpError::Domain("dense solve failed: NonFinite".into())),
    );
}

#[test]
fn constant_approximants_keep_each_entry_point_signed_zero() {
    let (numerator, denominator) = approx::pade(&[-0.0], 0, 0).unwrap();
    assert_eq!(numerator[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(denominator, [1.0]);
    let (numerator, denominator) = interp_ext::pade(&[-0.0], 0, 0).unwrap();
    assert_eq!(numerator[0].to_bits(), 0.0_f64.to_bits());
    assert_eq!(denominator, [1.0]);
}
