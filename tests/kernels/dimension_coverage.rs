//! Space-dimension coverage of the kernels: n-D dot products, tensor-product
//! interpolation (1D/2D/3D) and planar slices / streamlines of 3D fields.

use std::f64::consts::PI;
use std::path::PathBuf;

use rssn::io::plot::SliceAxis;
use rssn::io::plot::contour_lines_slice_3d;
use rssn::io::plot::plot_contour_slice_3d;
use rssn::io::plot::plot_heatmap_slice_3d;
use rssn::io::plot::plot_quiver_slice_3d;
use rssn::io::plot::plot_streamlines_3d;
use rssn::io::plot::slice_scalar_field_3d;
use rssn::io::plot::streamline_paths_3d;
use rssn::kernels::computer_graphics::Vector2D;
use rssn::kernels::computer_graphics::Vector3D;
use rssn::kernels::computer_graphics::cross_product;
use rssn::kernels::computer_graphics::cross_product_2d;
use rssn::kernels::computer_graphics::dot;
use rssn::kernels::computer_graphics::dot_product;
use rssn::kernels::computer_graphics::dot_product_2d;
use rssn::kernels::computer_graphics::dot_product_3d;
use rssn::kernels::computer_graphics::dot_product_nd;
use rssn::kernels::interpolate::SplineEnd;
use rssn::kernels::interpolate::TensorOrder;
use rssn::kernels::interpolate::bilinear_interpolation;
use rssn::kernels::interpolate::cubic_spline_with;
use rssn::kernels::interpolate::interpolate_nd;
use rssn::kernels::interpolate::multilinear_interpolation;
use rssn::kernels::interpolate::tensor_cubic_spline_interpolation;
use rssn::kernels::interpolate::trilinear_interpolation;

// ---------------------------------------------------------------------------
// dot products
// ---------------------------------------------------------------------------

#[test]
fn dot_products_agree_across_dimensions() {
    let a2 = Vector2D { x: 1.5, y: -2.0 };
    let b2 = Vector2D { x: 0.5, y: 4.0 };
    assert!((dot_product_2d(&a2, &b2) - (-7.25)).abs() < 1e-15);
    assert!((dot(&a2, &b2) - dot_product_2d(&a2, &b2)).abs() < 1e-15);

    let a3 = Vector3D {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let b3 = Vector3D {
        x: -4.0,
        y: 0.5,
        z: 2.0,
    };
    assert!((dot_product_3d(&a3, &b3) - dot_product(&a3, &b3)).abs() < 1e-15);
    assert!((dot(&a3, &b3) - 3.0).abs() < 1e-15);
    // the n-D routine is the same function on slices of any length
    assert!((dot_product_nd(&[1.5, -2.0], &[0.5, 4.0]).unwrap() - (-7.25)).abs() < 1e-15);
    assert!((dot_product_nd(&[1.0, 2.0, 3.0], &[-4.0, 0.5, 2.0]).unwrap() - 3.0).abs() < 1e-15);
    assert!((dot_product_nd(&[1.0; 7], &[2.0; 7]).unwrap() - 14.0).abs() < 1e-15);
    assert_eq!(dot_product_nd(&[], &[]).unwrap(), 0.0);
    assert!(dot_product_nd(&[1.0, 2.0], &[1.0]).is_err());
}

#[test]
fn cross_product_2d_is_the_z_component_of_the_3d_cross_product() {
    let a = Vector2D { x: 1.0, y: 2.0 };
    let b = Vector2D { x: -3.0, y: 0.5 };
    let c3 = cross_product(
        &Vector3D {
            x: a.x,
            y: a.y,
            z: 0.0,
        },
        &Vector3D {
            x: b.x,
            y: b.y,
            z: 0.0,
        },
    );
    assert!((cross_product_2d(&a, &b) - c3.z).abs() < 1e-15);
    assert!((cross_product_2d(&a, &a)).abs() < 1e-15);
}

// ---------------------------------------------------------------------------
// interpolation
// ---------------------------------------------------------------------------

fn axis(
    n: usize,
    lo: f64,
    hi: f64,
) -> Vec<f64> {
    (0..n)
        .map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64)
        .collect()
}

/// Values on a rectilinear grid, first axis fastest.
fn grid_values(
    axes: &[Vec<f64>],
    f: impl Fn(&[f64]) -> f64,
) -> Vec<f64> {
    let total: usize = axes.iter().map(Vec::len).product();
    (0..total)
        .map(|mut p| {
            let x: Vec<f64> = axes
                .iter()
                .map(|a| {
                    let v = a[p % a.len()];
                    p /= a.len();
                    v
                })
                .collect();
            f(&x)
        })
        .collect()
}

#[test]
fn multilinear_is_exact_for_functions_affine_in_each_variable() {
    // non-uniform axes on purpose
    let ax = vec![0.0, 0.3, 1.0, 1.1, 2.0];
    let ay = vec![-1.0, -0.2, 0.5, 2.0];
    let az = vec![0.0, 1.0, 1.5];
    let f = |x: &[f64]| {
        let (a, b, c) = (x[0], x[1], x[2]);
        1.0 + 2.0 * a - 3.0 * b + 0.5 * c + a * b - 2.0 * b * c + 0.7 * a * c + 1.5 * a * b * c
    };
    let axes = [ax, ay, az];
    let v = grid_values(&axes, f);
    for p in [
        [0.1, -0.7, 0.2],
        [1.05, 1.0, 1.4],
        [2.0, 2.0, 1.5],
        [0.0, -1.0, 0.0],
        [1.7, 0.1, 0.9],
    ] {
        let got = multilinear_interpolation(&axes, &v, &p).unwrap();
        assert!((got - f(&p)).abs() < 1e-12, "{p:?}");
        let tri =
            trilinear_interpolation(&axes[0], &axes[1], &axes[2], &v, p[0], p[1], p[2]).unwrap();
        assert!((tri - got).abs() < 1e-15);
    }
    // the 2D name is the 2D case
    let f2 = |x: &[f64]| 2.0 - x[0] + 3.0 * x[1] + 0.25 * x[0] * x[1];
    let v2 = grid_values(&[axes[0].clone(), axes[1].clone()], f2);
    let b = bilinear_interpolation(&axes[0], &axes[1], &v2, 0.77, 0.31).unwrap();
    assert!((b - f2(&[0.77, 0.31])).abs() < 1e-12);
    // nodes are reproduced exactly
    assert!(
        (bilinear_interpolation(&axes[0], &axes[1], &v2, 1.1, 0.5).unwrap() - f2(&[1.1, 0.5]))
            .abs()
            < 1e-13
    );
}

#[test]
fn multilinear_one_axis_is_piecewise_linear_interpolation_and_second_order() {
    let ax = vec![0.0, 1.0, 2.5, 3.0];
    let v = [0.0, 2.0, 1.0, 4.0];
    let got = multilinear_interpolation(std::slice::from_ref(&ax), &v, &[1.75]).unwrap();
    assert!((got - 1.5).abs() < 1e-15);
    // smooth 2D data: halving h divides the error by ~4
    let err = |n: usize| {
        let axes = [axis(n, 0.0, 2.0), axis(n, 0.0, 2.0)];
        let f = |x: &[f64]| (x[0]).sin() * (x[1]).cos();
        let v = grid_values(&axes, f);
        let mut worst = 0.0_f64;
        for i in 0..20 {
            for j in 0..20 {
                let p = [0.05 + 0.097 * i as f64, 0.03 + 0.098 * j as f64];
                worst =
                    worst.max((multilinear_interpolation(&axes, &v, &p).unwrap() - f(&p)).abs());
            }
        }
        worst
    };
    let order = (err(9) / err(17)).log2();
    assert!((1.7..2.3).contains(&order), "observed order {order}");
}

#[test]
fn tensor_cubic_spline_is_exact_for_cubics_and_fourth_order() {
    let axes = [axis(6, 0.0, 1.0), axis(5, -1.0, 1.0), axis(7, 0.0, 2.0)];
    let cubic = |x: &[f64]| {
        let (a, b, c) = (x[0], x[1], x[2]);
        1.0 + a - 2.0 * b * b + 0.5 * c * c * c + a * a * a - a * b * c * c + 3.0 * b * b * b
    };
    let v = grid_values(&axes, cubic);
    for p in [[0.13, -0.77, 1.31], [0.99, 0.45, 0.05], [0.5, 0.0, 1.0]] {
        let got = tensor_cubic_spline_interpolation(&axes, &v, &p).unwrap();
        assert!((got - cubic(&p)).abs() < 1e-10, "{p:?}");
        let via = interpolate_nd(&axes, &v, &p, TensorOrder::CubicSpline).unwrap();
        assert_eq!(got, via);
    }
    // a single axis is the not-a-knot spline of the 1D module
    let ax = axis(9, 0.0, PI);
    let ys: Vec<f64> = ax.iter().map(|x| x.sin()).collect();
    let pts: Vec<(f64, f64)> = ax.iter().copied().zip(ys.iter().copied()).collect();
    let spline = cubic_spline_with(&pts, SplineEnd::NotAKnot).unwrap();
    let one = tensor_cubic_spline_interpolation(std::slice::from_ref(&ax), &ys, &[1.234]).unwrap();
    assert!((one - spline(1.234)).abs() < 1e-14);
    // smooth 2D data: refinement h -> h/2 reduces the error by roughly 2^4
    let err = |n: usize| {
        let axes = [axis(n, 0.0, 2.0), axis(n, 0.0, 2.0)];
        let f = |x: &[f64]| (x[0]).sin() * (x[1]).cos();
        let v = grid_values(&axes, f);
        let mut worst = 0.0_f64;
        for i in 0..9 {
            for j in 0..9 {
                let p = [0.3 + 0.17 * i as f64, 0.2 + 0.18 * j as f64];
                worst = worst
                    .max((tensor_cubic_spline_interpolation(&axes, &v, &p).unwrap() - f(&p)).abs());
            }
        }
        worst
    };
    let order = (err(9) / err(17)).log2();
    assert!((3.0..5.0).contains(&order), "observed order {order}");
}

#[test]
fn tensor_interpolation_rejects_bad_input() {
    let axes = [vec![0.0, 1.0, 2.0], vec![0.0, 1.0]];
    let v = vec![0.0; 6];
    assert!(
        multilinear_interpolation(&axes, &v, &[0.5]).is_err(),
        "point dimension"
    );
    assert!(
        multilinear_interpolation(&axes, &v[..5], &[0.5, 0.5]).is_err(),
        "value count"
    );
    assert!(
        multilinear_interpolation(&axes, &v, &[2.5, 0.5]).is_err(),
        "outside"
    );
    assert!(
        multilinear_interpolation(&[vec![0.0, 0.0], vec![0.0, 1.0]], &[0.0; 4], &[0.0, 0.5])
            .is_err(),
        "repeated node"
    );
    assert!(multilinear_interpolation(&[], &[], &[]).is_err());
    assert!(
        tensor_cubic_spline_interpolation(&axes, &v, &[0.5, 0.5]).is_err(),
        "too few nodes"
    );
}

// ---------------------------------------------------------------------------
// 3D plots
// ---------------------------------------------------------------------------

struct TempFile(PathBuf);

impl TempFile {
    fn new(tag: &str) -> Self {
        Self(std::env::temp_dir().join(format!("rssn_dimcov_{}_{tag}.svg", std::process::id())))
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

const BOX: [(f64, f64); 3] = [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)];

#[test]
fn slice_sampling_places_the_in_plane_axes_in_order() {
    let f = |x: f64, y: f64, z: f64| 100.0 * x + 10.0 * y + z;
    // z = const: columns follow x, rows follow y
    let s = slice_scalar_field_3d(f, BOX, SliceAxis::Z, 0.5, 5).unwrap();
    assert_eq!(s.dim(), (5, 5));
    assert!((s[[0, 0]] - (-100.0 - 10.0 + 0.5)).abs() < 1e-12);
    assert!((s[[0, 4]] - (100.0 - 10.0 + 0.5)).abs() < 1e-12);
    assert!((s[[4, 0]] - (-100.0 + 10.0 + 0.5)).abs() < 1e-12);
    // y = const: in-plane axes are (x, z)
    let s = slice_scalar_field_3d(f, BOX, SliceAxis::Y, -0.5, 3).unwrap();
    assert!((s[[2, 1]] - (0.0 - 5.0 + 1.0)).abs() < 1e-12);
    // x = const: in-plane axes are (y, z)
    let s = slice_scalar_field_3d(f, BOX, SliceAxis::X, 1.0, 3).unwrap();
    assert!((s[[2, 0]] - (100.0 - 10.0 + 1.0)).abs() < 1e-12);
    assert!(slice_scalar_field_3d(f, BOX, SliceAxis::Z, 1.5, 5).is_err());
    assert!(slice_scalar_field_3d(f, BOX, SliceAxis::Z, 0.0, 1).is_err());
}

#[test]
fn contour_of_a_sphere_slice_is_a_circle_of_the_analytic_radius() {
    let r2 = |x: f64, y: f64, z: f64| x * x + y * y + z * z;
    for (axis, position) in [
        (SliceAxis::Z, 0.0),
        (SliceAxis::Z, 0.3),
        (SliceAxis::X, -0.4),
        (SliceAxis::Y, 0.2),
    ] {
        let lines = contour_lines_slice_3d(r2, BOX, axis, position, 81, &[0.25]).unwrap();
        assert!(!lines.is_empty(), "{axis:?}");
        let want = (0.25_f64 - position * position).sqrt();
        for line in &lines {
            assert!(line.points.len() > 20);
            for &(u, v) in &line.points {
                assert!(
                    (u.hypot(v) - want).abs() < 2e-3,
                    "{axis:?} at {position}: {}",
                    u.hypot(v)
                );
            }
        }
    }
    // above the pole there is no contour at this level
    let lines = contour_lines_slice_3d(r2, BOX, SliceAxis::Z, 0.9, 41, &[0.25]).unwrap();
    assert!(lines.iter().all(|l| l.points.is_empty()));
}

#[test]
fn slice_plots_write_svg_files() {
    let f = |x: f64, y: f64, z: f64| (x * y).sin() + z;
    let heat = TempFile::new("heat");
    plot_heatmap_slice_3d(f, BOX, SliceAxis::Z, 0.0, 16, &heat.0).unwrap();
    let text = std::fs::read_to_string(&heat.0).unwrap();
    assert!(text.starts_with("<svg") || text.contains("<svg"));
    assert!(text.matches("<rect").count() >= 256);

    let contour = TempFile::new("contour");
    plot_contour_slice_3d(f, BOX, SliceAxis::X, 0.2, 31, &[-0.3, 0.0, 0.3], &contour.0).unwrap();
    assert!(
        std::fs::read_to_string(&contour.0)
            .unwrap()
            .contains("<polyline")
    );

    let quiver = TempFile::new("quiver");
    plot_quiver_slice_3d(
        |x, y, z| (-y, x, z),
        BOX,
        SliceAxis::Z,
        0.0,
        7,
        0.8,
        &quiver.0,
    )
    .unwrap();
    assert!(
        std::fs::read_to_string(&quiver.0)
            .unwrap()
            .contains("<line")
    );
    assert!(
        plot_quiver_slice_3d(
            |x, y, z| (x, y, z),
            BOX,
            SliceAxis::Z,
            0.0,
            7,
            -1.0,
            &quiver.0
        )
        .is_err()
    );
}

#[test]
fn streamlines_3d_follow_a_helix_with_conserved_radius_and_linear_climb() {
    let ranges = [(-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)];
    let field = |x: f64, y: f64, _z: f64| (-y, x, 0.1);
    let seeds = [[1.0, 0.0, -1.0], [0.0, 0.5, 0.0]];
    let paths = streamline_paths_3d(field, ranges, &seeds, 0.05, 100).unwrap();
    assert_eq!(paths.len(), 2);
    for (path, seed) in paths.iter().zip(&seeds) {
        assert_eq!(path.len(), 101);
        let r0 = seed[0].hypot(seed[1]);
        for (k, p) in path.iter().enumerate() {
            assert!((p[0].hypot(p[1]) - r0).abs() < 1e-6);
            assert!((p[2] - (seed[2] + 0.1 * 0.05 * k as f64)).abs() < 1e-12);
        }
    }
    // paths stop at the box, at a stagnation point, and bad input is rejected
    let escaping =
        streamline_paths_3d(|_, _, _| (1.0, 0.0, 0.0), ranges, &[[0.0; 3]], 0.25, 100).unwrap();
    assert_eq!(escaping[0].len(), 9);
    let still =
        streamline_paths_3d(|_, _, _| (0.0, 0.0, 0.0), ranges, &[[0.0; 3]], 0.25, 10).unwrap();
    assert_eq!(still[0].len(), 1);
    assert!(streamline_paths_3d(field, ranges, &[[3.0, 0.0, 0.0]], 0.05, 10).is_err());
    assert!(streamline_paths_3d(field, ranges, &[], 0.05, 10).is_err());
    assert!(streamline_paths_3d(field, ranges, &seeds, 0.0, 10).is_err());

    let out = TempFile::new("stream");
    plot_streamlines_3d(field, ranges, &seeds, 0.05, 100, &out.0).unwrap();
    assert!(
        std::fs::read_to_string(&out.0)
            .unwrap()
            .contains("<polyline")
    );
}
