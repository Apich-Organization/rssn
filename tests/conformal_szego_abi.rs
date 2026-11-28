//! Prepared Szegő checks on the unit disk and invalid sample geometry.
use num_complex::Complex64;
use rssn::kernels::conformal::SzegoKerzmanSteinSolver;
use rssn::kernels::conformal::szego_kerzman_stein_map;
use std::f64::consts::PI;

fn circle(t: f64) -> Complex64 {
    Complex64::from_polar(1.0, t)
}
fn derivative(t: f64) -> Complex64 {
    Complex64::I * circle(t)
}

#[test]
fn prepared_unit_disk_recovers_map_and_both_kernels() {
    let center = Complex64::from(0.2);
    let point = Complex64::new(0.4, 0.3);
    let prepared = SzegoKerzmanSteinSolver::solve(circle, derivative, 64, center).unwrap();
    let exact = (point - center) / (Complex64::ONE - center.conj() * point);
    let szego = Complex64::from(1.0 / (2.0 * PI)) / (Complex64::ONE - center.conj() * point);
    let (image, s) = prepared.eval_map(point);
    assert!((image - exact).norm() < 1e-10);
    assert!((s - szego).norm() < 1e-10);
    assert!((prepared.eval_bergman(point) - 4.0 * PI * szego * szego).norm() < 1e-10);
    assert!((prepared.eval_szego(point) - s).norm() < 1e-15);
    let (again, boundary) = szego_kerzman_stein_map(circle, derivative, 64, center, point).unwrap();
    assert!((again - image).norm() < 1e-15 && (boundary - s).norm() < 1e-15);
}

#[test]
fn clockwise_collapsed_and_outside_geometry_fails() {
    let zero = Complex64::ZERO;
    assert!(SzegoKerzmanSteinSolver::solve(|t| circle(-t), |t| -derivative(-t), 32, zero).is_err());
    assert!(SzegoKerzmanSteinSolver::solve(circle, |_| zero, 32, zero).is_err());
    assert!(SzegoKerzmanSteinSolver::solve(circle, derivative, 32, Complex64::from(2.0)).is_err());
    let p = SzegoKerzmanSteinSolver::solve(circle, derivative, 32, zero).unwrap();
    assert!(!p.eval_map(Complex64::from(2.0)).0.is_finite());
}
