use rssn::kernels::krylov::Csr;
use rssn::kernels::sparse_native;

#[test]
fn sparse_matvec_survives_intermediate_overflow_and_cancellation() {
    let matrix = Csr::from_triplets(
        2,
        2,
        &[(0, 0, 1e308), (0, 1, 9e307), (1, 0, 9e307), (1, 1, 1e308)],
    );
    let mut out = [0.0; 2];
    matrix.matvec(&[2.0, -2.0], &mut out);
    assert!((out[0] / 2e307 - 1.0).abs() < 1e-15);
    assert!((out[1] / -2e307 - 1.0).abs() < 1e-15);
}

#[test]
fn sparse_matvec_preserves_small_cross_terms() {
    let matrix = Csr::from_triplets(1, 2, &[(0, 0, 1e300), (0, 1, 1e-300)]);
    let mut out = [0.0];
    matrix.matvec(&[1e-300, 1e300], &mut out);
    assert_eq!(out[0], 2.0);
}

#[test]
fn wide_csr_keeps_exact_bits_between_worker_counts_and_storage_views() {
    let (rows, cols) = (128, 512);
    let triplets: Vec<_> = (0..rows)
        .flat_map(|row| {
            (0..cols).map(move |col| {
                let weight = 1.0_f64 / ((1 + row + col) as f64);
                let value = if row == 0 && col < 2 {
                    if col == 0 { 1e308 } else { -1e308 }
                } else if (row + col) % 3 == 0 {
                    -weight
                } else {
                    weight
                };
                (row, col, value)
            })
        })
        .collect();
    let matrix = Csr::from_triplets(rows, cols, &triplets);
    let borrowed = sparse_native::Csr {
        rows,
        cols,
        row_ptr: &matrix.indptr,
        col_idx: &matrix.indices,
        values: &matrix.values,
    };
    let input = vec![2.0; cols];
    let compute = |workers| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(|| {
                let mut owned = vec![0.0; rows];
                matrix.matvec(&input, &mut owned);
                let checked = borrowed.matvec(&input).unwrap();
                let bits =
                    |values: Vec<f64>| values.into_iter().map(f64::to_bits).collect::<Vec<_>>();
                assert!(owned.iter().all(|v| v.is_finite()));
                let owned = bits(owned);
                assert_eq!(owned, bits(checked));
                owned
            })
    };
    assert_eq!(compute(1), compute(4));
}
