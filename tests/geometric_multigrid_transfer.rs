use rssn::kernels::multigrid::GeometricMg;
use rssn::kernels::multigrid::GeometricTransferError;
use rssn::kernels::multigrid::GeometricTransferWorkspace;

#[test]
fn geometric_hierarchy_rejects_non_finite_shifts() {
    for shift in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
        assert!(GeometricMg::new(2, 7, shift).is_none());
    }
    assert!(GeometricMg::new(2, 7, -0.0).is_some());
}

fn hat_weight(
    mut fine: usize,
    mut coarse: usize,
    nf: usize,
    nc: usize,
    dim: usize,
) -> f64 {
    let mut weight = 1.0;
    for _ in 0..dim {
        let f = fine % nf + 1;
        let c = 2 * (coarse % nc + 1);
        weight *= match f.abs_diff(c) {
            | 0 => 1.0,
            | 1 => 0.5,
            | _ => 0.0,
        };
        fine /= nf;
        coarse /= nc;
    }
    weight
}

#[test]
fn prepared_transfers_match_dense_tensor_weights_and_the_adjoint_identity() {
    let mut workspace = GeometricTransferWorkspace::default();
    for dim in (1..=3).rev() {
        let mg = GeometricMg::new(dim, 7, 0.125).unwrap();
        for (level, nf, nc) in [(0, 7, 3), (1, 3, 1)] {
            let fine: Vec<_> = (0..mg.level_size(level).unwrap())
                .map(|i| ((i * 31 + 7) % 19) as f64 - 9.0)
                .collect();
            let coarse: Vec<_> = (0..mg.level_size(level + 1).unwrap())
                .map(|i| ((i * 13 + 11) % 17) as f64 - 8.0)
                .collect();
            let mut restricted = vec![f64::NAN; coarse.len()];
            let mut prolonged = vec![f64::NAN; fine.len()];
            mg.restrict_into(level, &fine, &mut restricted, &mut workspace)
                .unwrap();
            mg.prolong_into(level, &coarse, &mut prolonged, &mut workspace)
                .unwrap();

            let normalization = 2.0_f64.powi(-(dim as i32));
            for (j, value) in restricted.iter().enumerate() {
                let expected: f64 = fine
                    .iter()
                    .enumerate()
                    .map(|(i, value)| hat_weight(i, j, nf, nc, dim) * value * normalization)
                    .sum();
                assert_eq!(
                    *value, expected,
                    "restriction dimension {dim}, level {level}, point {j}"
                );
            }
            for (i, value) in prolonged.iter().enumerate() {
                let expected: f64 = coarse
                    .iter()
                    .enumerate()
                    .map(|(j, value)| hat_weight(i, j, nf, nc, dim) * value)
                    .sum();
                assert_eq!(
                    *value, expected,
                    "interpolation dimension {dim}, level {level}, point {i}"
                );
            }
            let lhs: f64 = coarse.iter().zip(restricted).map(|(c, r)| c * r).sum();
            let rhs: f64 = fine
                .iter()
                .zip(prolonged)
                .map(|(f, p)| f * p * normalization)
                .sum();
            assert_eq!(lhs, rhs);
        }
    }
}

#[test]
fn restriction_preserves_tiny_uniform_weights_and_rounds_halfway_values() {
    let mut workspace = GeometricTransferWorkspace::default();
    let small = f64::from_bits(1);
    for dim in 1..=3 {
        let mg = GeometricMg::new(dim, 7, 0.0).unwrap();
        let mut target = vec![f64::NAN; mg.level_size(1).unwrap()];
        for constant in [small, -small, f64::MIN_POSITIVE, f64::MAX] {
            mg.restrict_into(0, &vec![constant; mg.size()], &mut target, &mut workspace)
                .unwrap();
            assert!(target.iter().all(|v| v.to_bits() == constant.to_bits()));
        }
    }

    let mg = GeometricMg::new(1, 3, 0.0).unwrap();
    for (input, expected) in [
        ([0.0, small, 0.0], 0.0),
        ([0.0, -small, 0.0], -0.0),
        ([small, 2.0 * small, small], 2.0 * small),
        ([-small, -2.0 * small, -small], -2.0 * small),
    ] {
        let mut target = [f64::NAN];
        mg.restrict_into(0, &input, &mut target, &mut workspace)
            .unwrap();
        assert_eq!(target[0].to_bits(), expected.to_bits());
    }
}

#[test]
fn prolongation_uses_finite_midpoints_and_zero_dirichlet_boundaries() {
    let mg = GeometricMg::new(1, 7, 0.0).unwrap();
    let mut workspace = GeometricTransferWorkspace::default();
    let large = f64::MAX;
    let mut fine = [f64::NAN; 7];
    mg.prolong_into(0, &[large; 3], &mut fine, &mut workspace)
        .unwrap();
    assert_eq!(
        fine,
        [0.5 * large, large, large, large, large, large, 0.5 * large]
    );

    mg.prolong_into(0, &[large, -large, 0.0], &mut fine, &mut workspace)
        .unwrap();
    assert_eq!(
        fine,
        [0.5 * large, large, 0.0, -large, -0.5 * large, 0.0, 0.0]
    );
}

#[test]
fn invalid_grid_shapes_and_levels_leave_outputs_untouched() {
    let mg = GeometricMg::new(2, 3, 0.0).unwrap();
    let mut workspace = GeometricTransferWorkspace::default();
    assert_eq!(mg.level_size(0), Some(mg.size()));
    assert_eq!(mg.level_size(1), Some(1));
    assert_eq!(mg.level_size(usize::MAX), None);

    let mut coarse = [42.0];
    for level in [1, usize::MAX] {
        assert_eq!(
            mg.restrict_into(level, &[0.0; 9], &mut coarse, &mut workspace),
            Err(GeometricTransferError::Level)
        );
    }
    assert_eq!(
        mg.restrict_into(0, &[0.0; 8], &mut coarse, &mut workspace),
        Err(GeometricTransferError::Shape)
    );
    assert_eq!(coarse, [42.0]);
    let mut fine = [42.0; 8];
    assert_eq!(
        mg.prolong_into(0, &coarse, &mut fine, &mut workspace),
        Err(GeometricTransferError::Shape)
    );
    assert_eq!(fine, [42.0; 8]);

    let n = (1_usize << (usize::BITS - 2)) - 1;
    assert!(GeometricMg::new(1, usize::MAX, 0.0).is_none());
    assert!(GeometricMg::new(3, n, 0.0).is_none());
}
