//! Untrusted portable descriptors must fail at the kernel preflight, never panic.
use rssn::kernels::interp_ext::bspline_eval;
#[test]
fn invalid_bspline_descriptors_fail_before_indexing() {
    assert!(bspline_eval(&[0.0, 0.0, 1.0], &[1.0], 1, 0.5).is_err());
    assert!(bspline_eval(&[0.0, 1.0], &[1.0], usize::MAX, 0.5).is_err());
    assert!(bspline_eval(&[0.0, 0.0], &[1.0], 0, 0.0).is_err());
    assert!(bspline_eval(&[0.0, f64::NAN], &[1.0], 0, 0.0).is_err());
    assert!(bspline_eval(&[0.0, 1.0], &[f64::INFINITY], 0, 0.5).is_err());
    assert!(bspline_eval(&[0.0, 1.0], &[1.0], 0, f64::NAN).is_err());
    assert_eq!(bspline_eval(&[0.0, 1.0], &[3.0], 0, 0.5).unwrap(), 3.0);
}
