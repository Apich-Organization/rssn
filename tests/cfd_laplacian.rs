//! Centered two-dimensional CFD Laplacian and preallocated output.

use rssn::kernels::matrix::Matrix;
use rssn::sim::physics_cfd::compute_laplacian;
use rssn::sim::physics_cfd::compute_laplacian_into;

#[test]
fn rectangular_axes_and_row_major_output_match_a_quadratic_field() {
    let (rows, cols) = (7, 12);
    let (dx, dy) = (0.125, 0.25);
    let values = (0..rows)
        .flat_map(|i| {
            let x = i as f64 * dx;
            (0..cols).map(move |j| {
                let y = j as f64 * dy;
                3.0 * x * x - 4.0 * y * y + 2.0 * x * y + 5.0
            })
        })
        .collect();
    let field = Matrix::new(rows, cols, values);
    let mut buffer = vec![f64::NAN; rows * cols];
    compute_laplacian_into(&field, dx, dy, &mut buffer).unwrap();
    let result = compute_laplacian(&field, -dx, -dy);
    assert_eq!((result.rows(), result.cols()), (rows, cols));
    for (i, row) in buffer.chunks_exact(cols).enumerate() {
        for (j, value) in row.iter().enumerate() {
            let expected: f64 = if i == 0 || i + 1 == rows || j == 0 || j + 1 == cols {
                0.0
            } else {
                -2.0
            };
            assert_eq!(value.to_bits(), expected.to_bits());
            assert_eq!(result.get(i, j).to_bits(), expected.to_bits());
        }
    }
    let mut wrong = [9.0; 3];
    assert!(compute_laplacian_into(&field, dx, dy, &mut wrong).is_err());
    assert_eq!(wrong.map(f64::to_bits), [9.0_f64.to_bits(); 3]);
}

#[test]
fn extreme_spacing_and_intermediate_overflow_keep_representable_results() {
    for dx in [1e-200, 1e-155, 1.0, 1e200, 1e308] {
        let field = Matrix::new(3, 4, [[dx; 4], [0.0; 4], [dx; 4]].concat());
        let result = compute_laplacian(&field, dx, 1.25);
        assert!((result.get(1, 2) * dx - 2.0).abs() < 5e-15);
        let result = compute_laplacian(&field, -dx, 1.25);
        assert!((result.get(1, 2) * dx - 2.0).abs() < 5e-15);
    }

    let m = f64::MAX;
    let field = Matrix::new(3, 3, [[m; 3], [m; 3], [0.0; 3]].concat());
    let result = compute_laplacian(&field, 1.0, 1.0);
    assert_eq!(result.get(1, 1).to_bits(), (-m).to_bits());
}

#[test]
fn empty_and_edge_only_grids_and_signed_zero_spacing() {
    for (rows, cols) in [(0, 0), (0, 4), (4, 0), (1, 7), (7, 1), (2, 5), (5, 2)] {
        let field = Matrix::new(rows, cols, vec![3.0; rows * cols]);
        let result = compute_laplacian(&field, 0.0, -0.0);
        assert_eq!((rows, cols), (result.rows(), result.cols()));
        assert!(result.data().iter().all(|value| value.to_bits() == 0));
        let mut output = vec![f64::NAN; rows * cols];
        compute_laplacian_into(&field, 0.0, -0.0, &mut output).unwrap();
        assert!(output.iter().all(|value| value.to_bits() == 0));
    }

    let field = Matrix::new(3, 3, [[1.0; 3], [0.0; 3], [1.0; 3]].concat());
    for spacing in [0.0, -0.0] {
        let result = compute_laplacian(&field, spacing, 1.0);
        assert_eq!(result.get(1, 1).to_bits(), f64::INFINITY.to_bits());
        let zero = compute_laplacian(&Matrix::zeros(3, 3), 1.0, spacing);
        assert!(zero.get(1, 1).is_nan());
    }
}

#[test]
fn serial_parallel_and_preallocated_rows_preserve_each_value() {
    let (rows, cols) = (149, 267);
    let field = Matrix::new(
        rows,
        cols,
        (0..rows * cols).map(|i| (i as f64 * 0.331).sin()).collect(),
    );
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let parallel = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let expected = serial.install(|| compute_laplacian(&field, 0.815, 1.125));
    let mut buffer = vec![42.0; rows * cols];
    parallel.install(|| compute_laplacian_into(&field, 0.815, 1.125, &mut buffer).unwrap());
    for (a, b) in expected.data().iter().zip(&buffer) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
