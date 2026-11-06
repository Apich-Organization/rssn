//! Public threshold incomplete-LU integration and validation.
use rssn::kernels::krylov::Csr;
use rssn::kernels::krylov::IlutError;
use rssn::kernels::krylov::IlutOptions;
use rssn::kernels::krylov::gmres;

fn opts(cap: usize) -> IlutOptions {
    IlutOptions {
        drop_tolerance: 0.0,
        max_entries_per_side: cap,
        pivot_tolerance: 0.0,
    }
}
fn mat() -> Csr {
    Csr::from_triplets(
        4,
        4,
        &[
            (0, 0, 4.0),
            (0, 1, 1.0),
            (0, 3, 0.5),
            (1, 0, 2.0),
            (1, 1, 5.0),
            (1, 2, -1.0),
            (2, 1, 3.0),
            (2, 2, 6.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
            (3, 2, 2.0),
            (3, 3, 7.0),
        ],
    )
}
#[test]
fn unlimited_fill_agrees_with_original_system_and_cap_is_enforced() {
    let a = mat();
    let factor = a.ilut(opts(4)).unwrap();
    let truth = [1.0, -2.0, 3.0, 0.5];
    let mut rhs = [0.0; 4];
    a.matvec(&truth, &mut rhs);
    let mut actual = [0.0; 4];
    factor.apply(&rhs, &mut actual).unwrap();
    for (a, b) in actual.iter().zip(truth) {
        assert!((a - b).abs() < 1e-12);
    }
    let mut check = [0.0; 4];
    a.matvec(&actual, &mut check);
    for (a, b) in check.iter().zip(rhs) {
        assert!((a - b).abs() < 1e-12);
    }
    let diagonal = a.ilut(opts(0)).unwrap();
    let stats = diagonal.stats();
    assert_eq!(stats.lower_entries, 0);
    assert_eq!(stats.upper_entries, 4);
    assert!(stats.dropped_entries > 0);
    let limited = a.ilut(opts(1)).unwrap();
    let stats = limited.stats();
    assert!(stats.lower_entries <= 4 && stats.upper_entries <= 8 && stats.dropped_entries > 0);
    assert_eq!(
        format!("{limited:?}"),
        format!("{:?}", a.ilut(opts(1)).unwrap())
    );
}
#[test]
fn malformed_matrices_and_bad_pivots_fail_without_panics() {
    let mut a = mat();
    a.cols = 3;
    assert!(matches!(a.ilut(opts(4)), Err(IlutError::InvalidInput)));
    assert!(a.ilu0().is_none());
    let mut a = mat();
    a.indptr[1] = usize::MAX;
    assert!(matches!(a.ilut(opts(4)), Err(IlutError::InvalidInput)));
    assert!(a.ilu0().is_none());
    let mut a = mat();
    a.indices[0] = usize::MAX;
    assert!(matches!(a.ilut(opts(4)), Err(IlutError::InvalidInput)));
    assert!(a.ilu0().is_none());
    let mut a = mat();
    a.values[0] = f64::NAN;
    assert!(matches!(a.ilut(opts(4)), Err(IlutError::NonFiniteInput)));
    assert!(a.ilu0().is_none());
    let a = Csr::from_triplets(2, 2, &[(0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]);
    assert!(matches!(
        a.ilut(opts(2)),
        Err(IlutError::UnsafePivot { row: 0 })
    ));
    let a = Csr::from_triplets(2, 2, &[(0, 0, 1e-16), (0, 1, 1.0), (1, 1, 1.0)]);
    assert!(matches!(
        a.ilut(IlutOptions {
            pivot_tolerance: 1e-12,
            ..opts(2)
        }),
        Err(IlutError::UnsafePivot { row: 0 })
    ));
}
#[test]
fn apply_is_transactional_and_handles_uniform_extreme_scales() {
    for scale in [1e-200, 1.0, 1e200] {
        let a = Csr::from_triplets(
            2,
            2,
            &[
                (0, 0, 4.0 * scale),
                (0, 1, scale),
                (1, 0, scale),
                (1, 1, 3.0 * scale),
            ],
        );
        let factor = a.ilut(opts(2)).unwrap();
        let mut rhs = [0.0; 2];
        a.matvec(&[2.0, -1.0], &mut rhs);
        let mut output = [7.0, 8.0];
        factor.apply(&rhs, &mut output).unwrap();
        assert!(
            (output[0] - 2.0).abs() < 1e-13 && (output[1] + 1.0).abs() < 1e-13,
            "{scale}: {output:?}"
        );
        let mut unchanged = [7.0, 8.0];
        assert!(factor.apply(&[f64::NAN, 0.0], &mut unchanged).is_err());
        assert_eq!(unchanged, [7.0, 8.0]);
        assert!(factor.apply(&[1.0], &mut unchanged).is_err());
        assert_eq!(unchanged, [7.0, 8.0]);
    }
}
#[test]
fn bounded_fill_preconditioner_accelerates_nonsymmetric_gmres() {
    let side = 12;
    let n = side * side;
    let mut trip = Vec::new();
    for i in 0..n {
        trip.push((i, i, 4.2));
        if i % side > 0 {
            trip.push((i, i - 1, -1.6));
        }
        if i % side + 1 < side {
            trip.push((i, i + 1, -0.4));
        }
        if i >= side {
            trip.push((i, i - side, -1.2));
        }
        if i + side < n {
            trip.push((i, i + side, -0.8));
        }
    }
    let a = Csr::from_triplets(n, n, &trip);
    let factor = a
        .ilut(IlutOptions {
            drop_tolerance: 1e-5,
            max_entries_per_side: 12,
            pivot_tolerance: 1e-14,
        })
        .unwrap();
    let rhs: Vec<_> = (0..n).map(|i| (i as f64 * 0.17).sin()).collect();
    let plain = gmres(
        |x, y| a.matvec(x, y),
        |r, z| z.copy_from_slice(r),
        &rhs,
        None,
        25,
        1e-9,
        500,
    );
    let pre = gmres(
        |x, y| a.matvec(x, y),
        |r, z| factor.apply(r, z).unwrap(),
        &rhs,
        None,
        25,
        1e-9,
        500,
    );
    let zero = a.ilu0().unwrap();
    let pre0 = gmres(
        |x, y| a.matvec(x, y),
        |r, z| zero.apply(r, z),
        &rhs,
        None,
        25,
        1e-9,
        500,
    );
    assert!(
        plain.converged && pre0.converged && pre.converged,
        "{plain:?} {pre0:?} {pre:?}"
    );
    assert!(
        pre.iterations < pre0.iterations,
        "ILU(0)={} ILUT={}",
        pre0.iterations,
        pre.iterations
    );
    assert!(
        pre.iterations < plain.iterations,
        "plain={} pre={}",
        plain.iterations,
        pre.iterations
    );
    let mut ax = vec![0.0; n];
    a.matvec(&pre.x, &mut ax);
    let nr = rhs
        .iter()
        .zip(&ax)
        .fold(0.0_f64, |acc, (b, a)| acc.hypot(b - a));
    let nb = rhs.iter().fold(0.0_f64, |acc, b| acc.hypot(*b));
    assert!(nr / nb <= 1e-9, "{}", nr / nb);
}

#[test]
fn tiny_triangular_numerator_is_divided_before_underflow() {
    let tiny = f64::from_bits(1);
    let matrix = Csr::from_triplets(2, 2, &[(0, 0, tiny), (0, 1, tiny), (1, 1, 1.0)]);
    let factor = matrix.ilut(opts(2)).unwrap();
    let mut out = [0.0; 2];
    factor.apply(&[tiny, 0.5], &mut out).unwrap();
    assert_eq!(out, [0.5, 0.5]);
}

#[test]
fn arithmetic_breakdowns_and_invalid_options_are_explicit() {
    let a = mat();
    for tolerance in [f64::NAN, -0.1, 1.1] {
        assert!(matches!(
            a.ilut(IlutOptions {
                drop_tolerance: tolerance,
                ..opts(4)
            }),
            Err(IlutError::InvalidInput)
        ));
        assert!(matches!(
            a.ilut(IlutOptions {
                pivot_tolerance: tolerance,
                ..opts(4)
            }),
            Err(IlutError::InvalidInput)
        ));
    }
    let tiny = f64::from_bits(1);
    let a = Csr::from_triplets(1, 1, &[(0, 0, tiny)]);
    let f = a.ilut(opts(1)).unwrap();
    let mut out = [42.0];
    assert!(matches!(
        f.apply(&[f64::MAX], &mut out),
        Err(IlutError::NumericalBreakdown { row: 0 })
    ));
    assert_eq!(out, [42.0]);
    let a = Csr::from_triplets(
        3,
        3,
        &[
            (0, 0, 2.0),
            (0, 2, tiny),
            (1, 0, 1.0),
            (1, 1, 1.0),
            (2, 2, 1.0),
        ],
    );
    assert!(matches!(
        a.ilut(opts(2)),
        Err(IlutError::NumericalBreakdown { row: 1 })
    ));
    let a = Csr::from_triplets(2, 2, &[(0, 0, f64::MAX), (1, 0, tiny), (1, 1, 1.0)]);
    assert!(a.ilu0().is_none());
}

#[test]
fn zero_drop_full_fill_matches_random_diagonally_dominant_systems() {
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    for _ in 0..80 {
        let n = 6;
        let mut trip = Vec::new();
        for row in 0..n {
            let mut sum = 0.0;
            for col in 0..n {
                if row == col {
                    continue;
                }
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                if seed % 3 != 0 {
                    let value = ((seed % 2001) as f64 - 1000.0) / 1000.0;
                    sum += value.abs();
                    trip.push((row, col, value));
                }
            }
            trip.push((row, row, sum + 1.0));
        }
        let a = Csr::from_triplets(n, n, &trip);
        let factor = a.ilut(opts(n)).unwrap();
        let expected: Vec<_> = (0..n).map(|i| (i as f64 * 0.41).sin()).collect();
        let mut rhs = vec![0.0; n];
        a.matvec(&expected, &mut rhs);
        let mut actual = vec![0.0; n];
        factor.apply(&rhs, &mut actual).unwrap();
        let mut check = vec![0.0; n];
        a.matvec(&actual, &mut check);
        let error = check
            .iter()
            .zip(&rhs)
            .fold(0.0_f64, |acc, (a, b)| acc.max((a - b).abs()));
        assert!(error < 1e-12, "error {error}");
    }
}
