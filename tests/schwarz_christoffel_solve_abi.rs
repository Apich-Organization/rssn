//! Prepared polygon solver scale and topology regressions.
use num_complex::Complex64;
use rssn::kernels::conformal::solve_sc_uhp;

#[test]
fn square_and_triangle_are_consistent_under_scaling() {
    let square = [
        Complex64::ZERO,
        Complex64::ONE,
        Complex64::new(1.0, 1.0),
        Complex64::I,
    ];
    let point = Complex64::new(0.5, 0.5);
    for polygon in [&square[..], &square[..3]] {
        let control = solve_sc_uhp(polygon, 40, 1e-7)
            .unwrap()
            .eval(point)
            .unwrap();
        assert!(control.re.is_finite() && control.im.is_finite());
        for extent in [1e-120, 1e120] {
            let p: Vec<_> = polygon.iter().map(|&v| v * extent).collect();
            let actual = solve_sc_uhp(&p, 40, 1e-7).unwrap().eval(point).unwrap() / extent;
            assert!((actual - control).norm() < 1e-11);
        }
    }
}

#[test]
fn clockwise_degenerate_and_self_intersecting_polygons_are_rejected() {
    let cw = [
        Complex64::ZERO,
        Complex64::I,
        Complex64::new(1.0, 1.0),
        Complex64::ONE,
    ];
    assert!(solve_sc_uhp(&cw, 40, 1e-7).is_err());
    let cross = [cw[0], cw[2], cw[1], cw[3]];
    assert!(solve_sc_uhp(&cross, 40, 1e-7).is_err());
    let collinear = [Complex64::ZERO, Complex64::ONE, Complex64::from(2.0)];
    assert!(solve_sc_uhp(&collinear, 40, 1e-7).is_err());
}
