//! Shape and reproducibility contracts for packed Householder columns.

use rayon::ThreadPoolBuilder;
use rssn::kernels::dense::Mat;
use rssn::kernels::lowrank::Gaussian;
use rssn::kernels::lowrank::LowRankError;
use rssn::kernels::lowrank::interpolative_decomposition;
use rssn::kernels::lowrank::qr_pivoted;

#[test]
fn empty_and_rank_zero_factors_keep_the_original_extents() {
    for (rows, cols) in [(0, 0), (0, 7), (6, 0), (7, 1), (1, 7)] {
        let matrix = Mat::zeros(rows, cols);
        let factor = qr_pivoted(&matrix).unwrap();
        let thin = rows.min(cols);
        assert_eq!((factor.q.rows, factor.q.cols), (rows, thin));
        assert_eq!((factor.r.rows, factor.r.cols), (thin, cols));
        assert_eq!(factor.rank(None), 0);
        let id = interpolative_decomposition(&matrix, 0).unwrap();
        assert!(id.columns.is_empty());
        assert_eq!((id.coefficients.rows, id.coefficients.cols), (0, cols));
    }
    let malformed = Mat {
        rows: 3,
        cols: 2,
        data: vec![1.0; 5],
    };
    assert!(matches!(
        qr_pivoted(&malformed),
        Err(LowRankError::Shape(_))
    ));
    assert!(matches!(
        interpolative_decomposition(&malformed, 0),
        Err(LowRankError::Shape(_))
    ));
}

#[test]
fn exact_equal_norm_ties_retain_the_existing_pivot_order() {
    let identity = Mat::identity(3);
    assert_eq!(qr_pivoted(&identity).unwrap().perm, [2, 0, 1]);
    assert_eq!(
        interpolative_decomposition(&identity, 2).unwrap().columns,
        [2, 0]
    );
}

#[test]
fn seeded_full_and_truncated_factors_are_exact_across_worker_counts() {
    let (rows, cols) = (340, 615);
    let mut generator = Gaussian::new(837);
    let a = Mat {
        rows,
        cols,
        data: (0..rows * cols).map(|_| generator.sample()).collect(),
    };
    let snapshot = |threads| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let qr = qr_pivoted(&a).unwrap();
                let id = interpolative_decomposition(&a, 48).unwrap();
                let bits =
                    qr.q.data
                        .into_iter()
                        .chain(qr.r.data)
                        .chain(id.coefficients.data)
                        .map(f64::to_bits)
                        .collect::<Vec<_>>();
                (qr.perm, id.columns, bits)
            })
    };
    let (pivots, columns, bits) = snapshot(1);
    assert_eq!(pivots[..48], columns);
    assert_eq!(bits.len(), rows * rows + rows * cols + 48 * cols);
    assert_eq!((pivots, columns, bits), snapshot(4));
}
