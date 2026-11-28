//! Symmetric RBF fitting with finite-scale coordinate normalization.

use rssn::kernels::interp_ext::InterpError;
use rssn::kernels::interp_ext::RbfKernel;
use rssn::kernels::interp_ext::rbf_interp;

fn close(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}

#[test]
fn gaussian_fits_nonzero_distances_at_finite_scale_extremes() {
    let midpoint = 4.0 * (-0.25_f64).exp() / (1.0 + (-1.0_f64).exp());
    for scale in [1e-200, 1.0, 1e200] {
        close(
            rbf_interp(
                &[&[0.0], &[scale]],
                &[2.0, 2.0],
                &[0.5 * scale],
                RbfKernel::Gaussian(scale),
            )
            .unwrap(),
            midpoint,
            3e-14,
        );
    }
    let expected = 4.0 * (-1.0_f64).exp() / (1.0 + (-4.0_f64).exp());
    close(
        rbf_interp(
            &[&[-1e308], &[1e308]],
            &[2.0, 2.0],
            &[0.0],
            RbfKernel::Gaussian(1e308),
        )
        .unwrap(),
        expected,
        3e-14,
    );
}

#[test]
fn multiquadric_weights_keep_finite_large_and_tiny_values() {
    let large = rbf_interp(&[&[0.0]], &[2.0], &[1e200], RbfKernel::Multiquadric(1.0)).unwrap();
    close(large / 1e200, 2.0, 2e-15);
    let tiny = rbf_interp(
        &[&[0.0]],
        &[2.0],
        &[1e308],
        RbfKernel::InverseMultiquadric(0.25),
    )
    .unwrap();
    close(tiny / 5e-309, 1.0, 4e-15);
}

#[test]
fn all_kernels_interpolate_sample_sites() {
    let points: [&[f64]; 4] = [&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1.2, 0.8]];
    let values = [1.75, -1.0, 0.6, 1.4];
    let kernels = [
        RbfKernel::Gaussian(0.9),
        RbfKernel::Multiquadric(0.9),
        RbfKernel::InverseMultiquadric(0.9),
        RbfKernel::ThinPlateSpline,
        RbfKernel::Cubic,
        RbfKernel::Linear,
    ];
    for kernel in kernels {
        for (&point, expected) in points.iter().zip(values) {
            close(
                rbf_interp(&points, &values, point, kernel).unwrap(),
                expected,
                4e-13,
            );
        }
    }
}

#[test]
fn invalid_kernel_controls_and_singular_sites_are_reported() {
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rbf_interp(&[&[0.0]], &[1.0], &[1.0], RbfKernel::Gaussian(scale)),
            Err(InterpError::Domain(_))
        ));
    }
    assert!(matches!(
        rbf_interp(&[&[f64::NAN]], &[1.0], &[0.0], RbfKernel::Linear),
        Err(InterpError::Domain(_))
    ));
    assert!(matches!(
        rbf_interp(&[&[0.0, 1.0]], &[1.0], &[0.0], RbfKernel::Linear),
        Err(InterpError::ShapeMismatch(_))
    ));
    assert!(matches!(
        rbf_interp(
            &[&[0.0], &[0.0]],
            &[1.0, 2.0],
            &[0.0],
            RbfKernel::Gaussian(1.0)
        ),
        Err(InterpError::Singular)
    ));
}

#[test]
fn symmetric_matrix_rows_are_worker_independent() {
    let points: Vec<Vec<_>> = (0..128)
        .map(|i| {
            (0..7)
                .map(|axis| f64::from((i * (axis + 23) + 111 * axis) % 1279) / 1279.0)
                .collect()
        })
        .collect();
    let refs: Vec<_> = points.iter().map(Vec::as_slice).collect();
    let values: Vec<_> = (0..128).map(|i| (f64::from(i) * 0.51).sin()).collect();
    let query = [0.21, 0.53, 0.38, 0.73, 0.96, 0.59, 0.12];
    let evaluate = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| rbf_interp(&refs, &values, &query, RbfKernel::Gaussian(0.23)).unwrap())
    };
    let single = evaluate(1);
    assert!(single.is_finite());
    assert_eq!(single.to_bits(), evaluate(4).to_bits());
}
