use rssn::kernels::krylov::Csr;

#[test]
fn separable_square_scale_recovers_both_representable_spd_pivots() {
    for (d, a, b) in [(-900, -200, 501), (900, 200, -499)] {
        for sign in [-1.0, 1.0] {
            let d = 2.0_f64.powi(d);
            let a = sign * 2.0_f64.powi(a);
            let b = 2.0_f64.powi(b);
            let matrix = Csr::from_triplets(2, 2, &[(0, 0, d), (1, 0, a), (1, 1, b)]);
            let factor = matrix.ic0().unwrap();
            let mut solution = [-1.0; 2];
            // The full two-dimensional lower triangle makes this an exact
            // LDLᵀ solve. Its RHS is the second column of A.
            factor.apply(&[a, b], &mut solution);
            assert_eq!(solution.map(f64::to_bits), [0.0_f64, 1.0].map(f64::to_bits));
        }
    }
}

#[test]
fn wide_lower_factors_share_the_correct_pivot_in_later_rows() {
    // Construct A from the exact dyadic factor data
    // D = (2^-900, 2^500, 2^498), L10=2^700, L20=2^699, L21=1/2.
    let a20 = 2.0_f64.powi(-201);
    let a21 = 2.0_f64.powi(500);
    let a22 = 3.0 * 2.0_f64.powi(498);
    let matrix = Csr::from_triplets(
        3,
        3,
        &[
            (0, 0, 2.0_f64.powi(-900)),
            (1, 0, 2.0_f64.powi(-200)),
            (1, 1, 2.0_f64.powi(501)),
            (2, 0, a20),
            (2, 1, a21),
            (2, 2, a22),
        ],
    );
    let factor = matrix.ic0().unwrap();
    let mut solution = [-1.0; 3];
    factor.apply(&[a20, a21, a22], &mut solution);
    assert_eq!(
        solution.map(f64::to_bits),
        [0.0_f64, 0.0, 1.0].map(f64::to_bits)
    );
}

#[test]
fn malformed_csr_fails_without_indexing_past_storage() {
    for (indptr, indices) in [
        (vec![0, 1], vec![0, 1]),
        (vec![0, 3, 2], vec![0, 1]),
        (vec![0, 1, 2], vec![0, 2]),
        (vec![0, 2, 3], vec![0, 0, 1]),
        (vec![0, 1, 3], vec![0, 1, 0]),
    ] {
        let matrix = Csr {
            rows: 2,
            cols: 2,
            values: vec![1.0; indices.len()],
            indptr,
            indices,
        };
        assert!(matrix.ic0().is_none());
    }
    let empty = Csr::from_triplets(0, 0, &[]);
    empty.ic0().unwrap().apply(&[], &mut []);
}

#[test]
fn upper_triangle_values_are_ignored_and_unrepresentable_factors_are_rejected() {
    let matrix = Csr::from_triplets(
        2,
        2,
        &[(0, 0, 2.0), (0, 1, f64::NAN), (1, 0, 1.0), (1, 1, 2.5)],
    );
    let mut solution = [0.0; 2];
    matrix.ic0().unwrap().apply(&[1.0, 2.5], &mut solution);
    assert_eq!(solution.map(f64::to_bits), [0.0_f64, 1.0].map(f64::to_bits));

    // The positive matrix can be represented as CSR, but its L10 = 2^1034
    // cannot be stored in a Float.
    let unsupported = Csr::from_triplets(
        2,
        2,
        &[
            (0, 0, f64::from_bits(1)),
            (1, 0, 2.0_f64.powi(-40)),
            (1, 1, 2.0_f64.powi(995)),
        ],
    );
    assert!(unsupported.ic0().is_none());
}
