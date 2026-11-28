//! Christoffel symbols and curvature from closure metrics.

use std::cell::Cell;
use std::cell::RefCell;
use std::collections::HashSet;

use rssn::kernels::differential_geometry::christoffel_symbols;
use rssn::kernels::differential_geometry::metric_tensor_at_point;
use rssn::kernels::differential_geometry::ricci_curvature;
use rssn::kernels::differential_geometry::ricci_scalar;
use rssn::kernels::differential_geometry::ricci_tensor;
use rssn::kernels::differential_geometry::riemann_tensor;

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
    assert!(
        ricci_scalar(polar, &p)
            .unwrap_or_else(|e| panic!("{e}"))
            .abs()
            < 1e-6
    );
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

#[test]
fn ricci_reuses_metric_samples_and_matches_a_full_riemann_contraction() {
    let point = [0.9, 0.3];
    let metric = sphere(1.5);
    let seen = RefCell::new(HashSet::new());
    let curvature = ricci_curvature(
        |point| {
            let key = point
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>();
            assert!(
                seen.borrow_mut().insert(key),
                "duplicate metric evaluation at {point:?}"
            );
            metric(point)
        },
        &point,
    )
    .unwrap_or_else(|error| panic!("{error}"));
    assert!(
        seen.borrow().len() < 81,
        "{} metric evaluations",
        seen.borrow().len()
    );
    assert_eq!(curvature.tensor, ricci_tensor(&metric, &point).unwrap());
    assert_eq!(curvature.scalar, ricci_scalar(&metric, &point).unwrap());

    let full = riemann_tensor(&metric, &point).unwrap();
    for (sigma, row) in curvature.tensor.iter().enumerate() {
        for (nu, &value) in row.iter().enumerate() {
            let expected: f64 = (0..point.len()).map(|mu| full[mu][sigma][mu][nu]).sum();
            assert!(
                (value - expected).abs() < 1e-12,
                "{value} versus {expected}"
            );
        }
    }
}

#[test]
fn uniform_metric_scale_changes_only_the_scalar_contraction() {
    let point = [1.1, 0.3];
    let base = ricci_curvature(sphere(1.0), &point).unwrap();
    for scale in [1e-200, 1e200, 1e308] {
        let metric =
            |point: &[f64]| vec![vec![scale, 0.0], vec![0.0, scale * point[0].sin().powi(2)]];
        let result = ricci_curvature(metric, &point).unwrap();
        assert!(
            (result.scalar * scale - 2.0).abs() < 1e-6,
            "scale {scale}, scalar {}",
            result.scalar
        );
        for (value, expected) in result
            .tensor
            .iter()
            .flatten()
            .zip(base.tensor.iter().flatten())
        {
            assert!((value - expected).abs() < 1e-8);
        }
    }

    // A symmetric indefinite metric may have a zero diagonal and a valid pivot.
    let lorentz = |_: &[f64]| vec![vec![0.0, 1e-200], vec![1e-200, 0.0]];
    let result = ricci_curvature(lorentz, &[0.0, 0.0]).unwrap();
    assert_eq!(result.scalar, 0.0);
    assert_eq!(result.tensor, vec![vec![0.0, 0.0], vec![0.0, 0.0]]);
}

#[test]
fn a_three_dimensional_conformal_metric_has_the_expected_cross_term() {
    let gradient = [0.31, -0.27, 0.0];
    let point = [0.4, -0.3, 0.6];
    let scale_at =
        |point: &[f64]| (2.0 * gradient.iter().zip(point).map(|(a, b)| a * b).sum::<f64>()).exp();
    let metric = |point: &[f64]| {
        let scale = scale_at(point);
        vec![
            vec![scale, 0.0, 0.0],
            vec![0.0, scale, 0.0],
            vec![0.0, 0.0, scale],
        ]
    };
    let curvature = ricci_curvature(metric, &point).unwrap();
    let gradient_square: f64 = gradient.iter().map(|value| value * value).sum();
    for (i, row) in curvature.tensor.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            let expected = gradient[i] * gradient[j]
                - if i == j {
                    gradient_square
                } else {
                    0.0
                };
            assert!((value - expected).abs() < 1e-8, "{value} versus {expected}");
        }
    }
    let expected_scalar = -2.0 * gradient_square / scale_at(&point);
    assert!(
        (curvature.scalar - expected_scalar).abs() < 1e-8,
        "{} versus {expected_scalar}",
        curvature.scalar
    );
}

#[test]
fn every_stencil_sample_must_be_finite_symmetric_and_square() {
    let point = [0.0, 0.0];
    let small_asymmetric = |_: &[f64]| vec![vec![1e-200, 1e-201], vec![0.0, 1e-200]];
    assert!(
        metric_tensor_at_point(small_asymmetric, &point)
            .unwrap_err()
            .contains("symmetric")
    );
    let wrong_shape = |point: &[f64]| {
        if point[0] == 0.0 {
            vec![vec![1.0, 0.0], vec![0.0, 1.0]]
        } else {
            vec![vec![1.0]]
        }
    };
    assert!(
        christoffel_symbols(wrong_shape, &point)
            .unwrap_err()
            .contains("2x2")
    );
    assert!(
        ricci_curvature(wrong_shape, &point)
            .unwrap_err()
            .contains("2x2")
    );
    let asymmetric_stencil = |point: &[f64]| vec![vec![1.0, point[0]], vec![0.0, 1.0]];
    assert!(
        riemann_tensor(asymmetric_stencil, &point)
            .unwrap_err()
            .contains("symmetric")
    );
    let nonfinite_stencil = |point: &[f64]| {
        let first = if point[0] == 0.0 {
            1.0
        } else {
            f64::NAN
        };
        vec![vec![first, 0.0], vec![0.0, 1.0]]
    };
    assert!(
        ricci_tensor(nonfinite_stencil, &point)
            .unwrap_err()
            .contains("finite")
    );
}

#[test]
fn nonfinite_coordinates_never_reach_the_callback() {
    let calls = Cell::new(0_usize);
    let metric = |point: &[f64]| {
        calls.set(calls.get() + 1);
        assert!(point.iter().all(|value| value.is_finite()));
        vec![vec![1.0]]
    };
    assert!(metric_tensor_at_point(metric, &[f64::NAN]).is_err());
    assert_eq!(calls.get(), 0);
    assert!(christoffel_symbols(metric, &[f64::MAX]).is_err());
    assert_eq!(calls.get(), 1);
}
