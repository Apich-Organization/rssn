//! Regression tests for clamped bilinear interpolation at finite extremes.

use rssn::kernels::interp_ext::interp2;

#[test]
fn clamped_bilinear_interpolation_scales_extreme_intervals() {
    let extreme = [-f64::MAX, f64::MAX];
    let regular = [0.0, 4.0];
    let values = [0.0, 2.0, 6.0, 8.0];
    // The plane has value 2*u + 6*v in normalized cell coordinates.
    for (xs, ys, x, y, expected) in [
        (extreme, regular, 0.0, 2.0, 4.0),
        (regular, extreme, 2.0, f64::MAX / 2.0, 5.5),
        (extreme, extreme, -f64::MAX / 2.0, 0.0, 3.5),
        (extreme, regular, f64::MAX / 2.0, -1.0, 1.5),
        (extreme, regular, f64::MAX, 5.0, 8.0),
    ] {
        assert_eq!(interp2(&xs, &ys, &values, x, y).unwrap(), expected);
    }
}
