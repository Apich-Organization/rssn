//! Ordered contractions and callback behavior for finite-difference curvature.

use std::cell::RefCell;

use rayon::ThreadPoolBuilder;
use rssn::kernels::differential_geometry::christoffel_symbols;
use rssn::kernels::differential_geometry::riemann_tensor;

fn product_of_spheres(point: &[f64]) -> Vec<Vec<f64>> {
    let mut metric = vec![vec![0.0; point.len()]; point.len()];
    for (axis, row) in metric.iter_mut().enumerate() {
        row[axis] = if axis % 2 == 0 {
            1.0
        } else {
            point[axis - 1].sin().powi(2)
        };
    }
    metric
}

#[test]
fn sphere_product_curvature_has_exact_last_pair_antisymmetry() {
    let point = [0.7, 0.1, 1.1, 0.2];
    let tensor = riemann_tensor(product_of_spheres, &point).unwrap();
    assert!((tensor[0][1][0][1] - 0.7_f64.sin().powi(2)).abs() < 1e-8);
    assert!((tensor[1][0][0][1] + 1.0).abs() < 1e-8);
    assert!((tensor[2][3][2][3] - 1.1_f64.sin().powi(2)).abs() < 1e-8);
    assert!(tensor[0][2][0][2].abs() < 1e-12);
    for plane in tensor.iter().flatten() {
        for (mu, row) in plane.iter().enumerate() {
            assert_eq!(row[mu].to_bits(), 0.0_f64.to_bits());
            for (nu, value) in row.iter().enumerate().skip(mu + 1) {
                if *value == 0.0 {
                    assert_eq!(plane[nu][mu], 0.0);
                } else {
                    assert_eq!(plane[nu][mu].to_bits(), (-value).to_bits());
                }
            }
        }
    }
}

#[test]
fn workers_leave_metric_callback_local_and_preserve_each_output_bit() {
    let point = [0.9; 12];
    let evaluate = |threads| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let caller = std::thread::current().id();
                // The callback can hold non-Sync state and its order is visible.
                let samples = RefCell::new(Vec::new());
                let result = riemann_tensor(
                    |p| {
                        assert_eq!(std::thread::current().id(), caller);
                        samples.borrow_mut().push(p.to_vec());
                        product_of_spheres(p)
                    },
                    &point,
                )
                .unwrap();
                let bits = result
                    .into_iter()
                    .flatten()
                    .flatten()
                    .flatten()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>();
                (samples.into_inner(), bits)
            })
    };
    assert_eq!(evaluate(1), evaluate(4));
}

#[test]
fn empty_charts_and_one_dimensional_curvature_keep_their_shapes() {
    let metric = |p: &[f64]| vec![vec![1.0 + p[0] * p[0]]];
    let first = christoffel_symbols(metric, &[0.4]).unwrap();
    assert!((first[0][0][0] - 0.4 / 1.16).abs() < 1e-12);
    let curvature = riemann_tensor(metric, &[0.4]).unwrap();
    assert_eq!(curvature, vec![vec![vec![vec![0.0]]]]);
    let empty = |_: &[f64]| Vec::new();
    assert!(christoffel_symbols(empty, &[]).unwrap().is_empty());
    assert!(riemann_tensor(empty, &[]).unwrap().is_empty());
}

#[test]
fn diagonal_shortcut_preserves_errors_for_overflowing_connection_products() {
    let metric = |p: &[f64]| vec![vec![1e-155 + p[0]]];
    let connection = christoffel_symbols(metric, &[0.0]).unwrap()[0][0][0];
    assert!(connection.is_finite() && connection.abs() > 1e154);
    assert_eq!(
        riemann_tensor(metric, &[0.0]).unwrap_err(),
        "Riemann tensor is outside the representable range"
    );
}
