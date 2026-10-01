//! Christoffel symbols and curvature from closure metrics.

use rssn::kernels::differential_geometry::{
    christoffel_symbols, metric_tensor_at_point, ricci_scalar, ricci_tensor, riemann_tensor,
};

fn polar(p: &[f64]) -> Vec<Vec<f64>> {
    vec![vec![1.0, 0.0], vec![0.0, p[0] * p[0]]]
}

fn sphere(a: f64) -> impl Fn(&[f64]) -> Vec<Vec<f64>> {
    move |p: &[f64]| vec![vec![a * a, 0.0], vec![0.0, a * a * p[0].sin().powi(2)]]
}

#[test]
fn metric_is_evaluated_and_checked() {
    let g = metric_tensor_at_point(polar, &[2.0, 0.3]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(g, vec![vec![1.0, 0.0], vec![0.0, 4.0]]);
    assert!(metric_tensor_at_point(polar, &[2.0, 0.3, 1.0]).is_err());
}

#[test]
fn polar_christoffel_symbols() {
    let r = 1.7;
    let gamma = christoffel_symbols(polar, &[r, 0.4]).unwrap_or_else(|e| panic!("{e}"));
    assert!((gamma[0][1][1] + r).abs() < 1e-8, "{}", gamma[0][1][1]);
    assert!((gamma[1][0][1] - 1.0 / r).abs() < 1e-8);
    assert!((gamma[1][1][0] - 1.0 / r).abs() < 1e-8);
    assert!(gamma[0][0][0].abs() < 1e-8 && gamma[1][1][1].abs() < 1e-8);
}

#[test]
fn flat_polar_plane_has_no_curvature() {
    let p = [1.3, 0.5];
    let riemann = riemann_tensor(polar, &p).unwrap_or_else(|e| panic!("{e}"));
    for value in riemann.iter().flatten().flatten().flatten() {
        assert!(value.abs() < 1e-6, "{value}");
    }
    assert!(ricci_scalar(polar, &p).unwrap_or_else(|e| panic!("{e}")).abs() < 1e-6);
}

#[test]
fn sphere_ricci_scalar_is_two_over_a_squared() {
    for a in [1.0, 2.5] {
        let r = ricci_scalar(sphere(a), &[1.1, 0.3]).unwrap_or_else(|e| panic!("{e}"));
        assert!((r - 2.0 / (a * a)).abs() < 1e-6, "a={a}: {r}");
    }
    // The Ricci tensor of a 2-sphere is (R/2) g.
    let a = 1.5;
    let ricci = ricci_tensor(sphere(a), &[0.9, 0.0]).unwrap_or_else(|e| panic!("{e}"));
    assert!((ricci[0][0] - 1.0).abs() < 1e-6);
    assert!((ricci[1][1] - 0.9_f64.sin().powi(2)).abs() < 1e-6);
    assert!(ricci[0][1].abs() < 1e-6);
}

#[test]
fn singular_metric_is_an_error() {
    let degenerate = |_: &[f64]| vec![vec![1.0, 1.0], vec![1.0, 1.0]];
    assert!(christoffel_symbols(degenerate, &[0.0, 0.0]).is_err());
}
