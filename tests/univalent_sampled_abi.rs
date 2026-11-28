//! Checked sampled univalence diagnostics and the shared angular grid.
use num_complex::Complex64;
use rssn::kernels::univalent::check_nehari_univalence;
use rssn::kernels::univalent::radius_of_convexity_numeric;
use rssn::kernels::univalent::radius_of_starlikeness_numeric;
use rssn::kernels::univalent::schwarzian_derivative_numeric;
use std::cell::Cell;

#[test]
fn schwarzian_quotient_uses_one_cauchy_circle() {
    let calls = Cell::new(0);
    let square = |z: Complex64| {
        calls.set(calls.get() + 1);
        z * z * 1e150
    };
    let result = schwarzian_derivative_numeric(&square, Complex64::from(1.0), 0.25).unwrap();
    assert!((result - Complex64::from(-1.5)).norm() < 1e-11);
    assert_eq!(calls.get(), 64);
    assert!(
        schwarzian_derivative_numeric(&|_| Complex64::from(2.0), Complex64::ZERO, 0.125).is_err()
    );
}

#[test]
fn sampling_rejects_nonfinite_functions_and_out_of_disk_circles() {
    let bad = |_| Complex64::new(f64::NAN, 0.0);
    assert!(check_nehari_univalence(&bad, &[0.125], 8, 0.125).is_err());
    assert!(check_nehari_univalence(&|z| z, &[0.99], 8, 0.125).is_err());
    assert!(check_nehari_univalence(&|z| z, &[], 8, 0.125).is_err());
    assert!(radius_of_convexity_numeric(&bad, &|_| Complex64::ZERO, 0.5, 8).is_nan());
}

#[test]
fn cached_grid_preserves_nested_radius_and_exact_work() {
    let f_calls = Cell::new(0);
    let d_calls = Cell::new(0);
    let map = |z| {
        f_calls.set(f_calls.get() + 1);
        z
    };
    let derivative = |_| {
        d_calls.set(d_calls.get() + 1);
        Complex64::from(1.0)
    };
    for angles in [8, 11] {
        let radius = radius_of_starlikeness_numeric(&map, &derivative, 1.0, angles);
        assert!((radius - 1.0).abs() < 1e-11);
        assert_eq!(f_calls.replace(0), angles * 40);
        assert_eq!(d_calls.replace(0), angles * 40);
    }
}
