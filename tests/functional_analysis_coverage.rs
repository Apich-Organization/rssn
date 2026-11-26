//! Coverage tests for `rssn::kernels::functional_analysis`.

use rssn::kernels::functional_analysis::L2Grid;
use rssn::kernels::functional_analysis::gram_schmidt;
use rssn::kernels::functional_analysis::gram_schmidt_orthonormal;
use rssn::kernels::functional_analysis::infinity_norm;
use rssn::kernels::functional_analysis::inner_product;
use rssn::kernels::functional_analysis::l1_norm;
use rssn::kernels::functional_analysis::l2_norm;
use rssn::kernels::functional_analysis::normalize;
use rssn::kernels::functional_analysis::project;

#[test]
fn norms_on_analytical_and_edge_inputs() {
    assert_eq!(l1_norm(&[]), 0.0);
    assert_eq!(l2_norm(&[]), 0.0);
    assert_eq!(infinity_norm(&[]), 0.0);
    assert_eq!(l1_norm(&[(0.0, 5.0)]), 0.0);
    assert_eq!(l2_norm(&[(0.0, 5.0)]), 0.0);
    assert_eq!(infinity_norm(&[(0.0, -5.0)]), 5.0);

    // f(x) = x on [0, 2]: ∫|x|dx = 2, (∫x^2 dx)^(1/2) = sqrt(8/3), ||f||_∞ = 2.
    let n = 401;
    let linear: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let x = 2.0 * i as f64 / (n - 1) as f64;
            (x, x)
        })
        .collect();
    assert!((l1_norm(&linear) - 2.0).abs() < 1e-4);
    assert!((l2_norm(&linear) - (8.0_f64 / 3.0).sqrt()).abs() < 1e-4);
    assert!((infinity_norm(&linear) - 2.0).abs() < 1e-12);

    // Repeated coordinates give zero-width panels and leave norms finite.
    let repeated = [(0.0, 1.0), (0.0, 2.0), (1.0, 2.0), (1.0, 3.0)];
    assert!((l1_norm(&repeated) - 2.0).abs() < 1e-12);
    assert!((l2_norm(&repeated) - 2.0).abs() < 1e-12);
}

#[test]
fn l2_grid_validates_and_normalizes_extreme_scales() {
    assert!(L2Grid::new(&[0.0, 2.0, 1.0]).is_err());
    assert!(L2Grid::new(&[0.0, f64::NAN]).is_err());

    let grid = L2Grid::new(&[0.0, 0.5, 1.0]).unwrap();
    assert!(grid.norm(&[1.0, 2.0]).is_err());
    assert!(grid.norm(&[1.0, f64::INFINITY, 1.0]).is_err());

    // Zero-measure function stays unchanged on normalize.
    let (z_norm, z_vec) = grid.norm_and_normalize(&[0.0, 0.0, 0.0]).unwrap();
    assert_eq!(z_norm, 0.0);
    assert_eq!(z_vec, vec![0.0, 0.0, 0.0]);

    // Extreme magnitudes avoid overflow/underflow during norm and normalization.
    for scale in [1e200, 1e-200] {
        let (norm, unit) = grid.norm_and_normalize(&[scale, scale, scale]).unwrap();
        assert!((norm / scale - 1.0).abs() < 1e-12, "scale {scale}");
        let check = grid.norm(&unit).unwrap();
        assert!((check - 1.0).abs() < 1e-12, "normalized norm at {scale}");
    }

    // Standalone normalize also preserves unit norm at extreme scales.
    let pts = [(0.0, 1e200), (1.0, 1e200)];
    let norm_pts = normalize(&pts);
    assert!((l2_norm(&norm_pts) - 1.0).abs() < 1e-12);
}

#[test]
fn inner_product_projection_and_gram_schmidt() {
    let f = [(0.0, 1.0), (0.5, 1.0), (1.0, 1.0)];
    let g = [(0.0, -1.0), (0.5, 0.0), (1.0, 1.0)];
    assert!(inner_product(&f, &[(0.0, 1.0)]).is_err());
    assert!(inner_product(&f, &[(0.0, 1.0), (0.6, 1.0), (1.0, 1.0)]).is_err());

    // <1, 2x - 1> on [0, 1] = 0.
    let ip = inner_product(&f, &g).unwrap();
    assert!(ip.abs() < 1e-12);

    let proj = project(&g, &f).unwrap();
    for &(_, y) in &proj {
        assert!(y.abs() < 1e-12);
    }

    // Zero target projects to zero.
    let zero = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0)];
    let proj_zero = project(&f, &zero).unwrap();
    assert!(proj_zero.iter().all(|&(_, y)| y == 0.0));

    let ortho = gram_schmidt(&[f.to_vec(), vec![(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)]]).unwrap();
    assert!(inner_product(&ortho[0], &ortho[1]).unwrap().abs() < 1e-12);

    let orthonorm =
        gram_schmidt_orthonormal(&[f.to_vec(), vec![(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)]]).unwrap();
    assert!((l2_norm(&orthonorm[0]) - 1.0).abs() < 1e-12);
    assert!((l2_norm(&orthonorm[1]) - 1.0).abs() < 1e-12);
}
