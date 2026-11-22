use super::*;
use crate::kernels::krylov::cg;
use crate::kernels::random::Rng;

fn random_csr(
    rng: &mut Rng,
    rows: usize,
    cols: usize,
    density: f64,
) -> Csr {
    let mut trip = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            if rng.uniform() < density {
                trip.push((r, c, rng.normal() + 0.1));
            }
        }
    }
    Csr::from_triplets(rows, cols, &trip)
}

fn dense(csr: &Csr) -> Vec<f64> {
    let mut out = vec![0.0; csr.rows * csr.cols];
    for r in 0..csr.rows {
        for k in csr.indptr[r]..csr.indptr[r + 1] {
            out[r * csr.cols + csr.indices[k]] = csr.values[k];
        }
    }
    out
}

#[test]
fn roundtrip_csr_bsr_csr() {
    let mut rng = Rng::new(1);
    for (br, bc) in [(1, 1), (2, 3), (3, 2), (4, 3), (12, 18)] {
        let csr = random_csr(&mut rng, 12, 18, 0.2);
        let bsr = Bsr::from_csr(&csr, br, bc).unwrap();
        assert_eq!(bsr.validate(), Ok(()));
        assert_eq!(bsr.to_csr(), csr);
        assert_eq!(bsr.nnz(), csr.values.len());
    }
}

#[test]
fn matvec_equals_csr_bit_for_bit() {
    let mut rng = Rng::new(2);
    let csr = random_csr(&mut rng, 24, 30, 0.15);
    let x: Vec<f64> = (0..30).map(|_| rng.normal()).collect();
    let mut want = vec![0.0; 24];
    csr.matvec(&x, &mut want);
    for (br, bc) in [(1, 1), (2, 3), (4, 5), (3, 6)] {
        let bsr = Bsr::from_csr(&csr, br, bc).unwrap();
        let mut got = vec![9.0; 24];
        bsr.matvec(&x, &mut got);
        assert_eq!(got, want, "block {br}x{bc}");
    }
}

#[test]
fn matvec_parallel_path_is_identical() {
    let mut rng = Rng::new(3);
    let n = 600;
    let csr = random_csr(&mut rng, n, n, 0.2);
    let bsr = Bsr::from_csr(&csr, 4, 4).unwrap();
    assert!(bsr.stored_len() >= BSR_PAR_NNZ);
    let x: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let (mut a, mut b) = (vec![0.0; n], vec![0.0; n]);
    csr.matvec(&x, &mut a);
    bsr.matvec(&x, &mut b);
    assert_eq!(a, b);
}

#[test]
fn transpose_and_matvec_transpose() {
    let mut rng = Rng::new(4);
    let csr = random_csr(&mut rng, 12, 18, 0.25);
    let bsr = Bsr::from_csr(&csr, 3, 2).unwrap();
    let t = bsr.transpose();
    assert_eq!((t.br, t.bc, t.nrows, t.ncols), (2, 3, 18, 12));
    assert_eq!(t.validate(), Ok(()));
    let d = dense(&csr);
    let dt = dense(&t.to_csr());
    for r in 0..12 {
        for c in 0..18 {
            assert_eq!(d[r * 18 + c], dt[c * 12 + r]);
        }
    }
    assert_eq!(t.transpose(), bsr);
    let x: Vec<f64> = (0..12).map(|_| rng.normal()).collect();
    let (mut y1, mut y2) = (vec![0.0; 18], vec![0.0; 18]);
    bsr.matvec_transpose(&x, &mut y1);
    t.matvec(&x, &mut y2);
    for (a, b) in y1.iter().zip(&y2) {
        assert!((a - b).abs() < 1e-12);
    }
}

#[test]
fn matmul_matches_dense() {
    let mut rng = Rng::new(5);
    let a = random_csr(&mut rng, 8, 12, 0.3);
    let b = random_csr(&mut rng, 12, 6, 0.3);
    let c = Bsr::from_csr(&a, 2, 3)
        .unwrap()
        .matmul(&Bsr::from_csr(&b, 3, 2).unwrap())
        .unwrap();
    assert_eq!(c.validate(), Ok(()));
    assert_eq!((c.br, c.bc), (2, 2));
    let (da, db, dc) = (dense(&a), dense(&b), dense(&c.to_csr()));
    for i in 0..8 {
        for j in 0..6 {
            let want: f64 = (0..12).map(|k| da[i * 12 + k] * db[k * 6 + j]).sum();
            assert!((dc[i * 6 + j] - want).abs() < 1e-12);
        }
    }
    let bad = Bsr::from_csr(&b, 1, 2).unwrap();
    assert!(matches!(
        Bsr::from_csr(&a, 2, 3).unwrap().matmul(&bad),
        Err(BsrError::ShapeMismatch(_))
    ));
}

/// A block-tridiagonal matrix of `n` block rows with random `b x b` blocks.
fn block_tridiagonal(
    rng: &mut Rng,
    n: usize,
    b: usize,
) -> Csr {
    let mut trip = Vec::new();
    for i in 0..n {
        for j in [i.wrapping_sub(1), i, i + 1] {
            if j >= n {
                continue;
            }
            for r in 0..b {
                for c in 0..b {
                    let mut v = rng.normal();
                    if i == j && r == c {
                        v += 6.0;
                    }
                    trip.push((i * b + r, j * b + c, v));
                }
            }
        }
    }
    Csr::from_triplets(n * b, n * b, &trip)
}

#[test]
fn bilu0_is_exact_without_fill() {
    let mut rng = Rng::new(6);
    let (n, b) = (7, 3);
    let csr = block_tridiagonal(&mut rng, n, b);
    let bsr = Bsr::from_csr(&csr, b, b).unwrap();
    let ilu = bsr.bilu0().unwrap();
    let x: Vec<f64> = (0..n * b).map(|_| rng.normal()).collect();
    let mut rhs = vec![0.0; n * b];
    bsr.matvec(&x, &mut rhs);
    let mut z = vec![0.0; n * b];
    ilu.apply(&rhs, &mut z);
    for (a, e) in z.iter().zip(&x) {
        assert!((a - e).abs() < 1e-10, "{a} vs {e}");
    }
}

#[test]
fn block_jacobi_inverts_block_diagonal_matrix() {
    let mut rng = Rng::new(7);
    let csr = block_tridiagonal(&mut rng, 1, 4);
    let bsr = Bsr::from_csr(&csr, 4, 4).unwrap();
    let pre = bsr.block_jacobi().unwrap();
    let x: Vec<f64> = (0..4).map(|_| rng.normal()).collect();
    let mut rhs = vec![0.0; 4];
    bsr.matvec(&x, &mut rhs);
    let mut z = vec![0.0; 4];
    pre.apply(&rhs, &mut z);
    for (a, e) in z.iter().zip(&x) {
        assert!((a - e).abs() < 1e-10);
    }
    let singular = Bsr::from_triplets(
        2,
        2,
        2,
        2,
        &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 4.0)],
    )
    .unwrap();
    assert_eq!(
        singular.block_jacobi(),
        Err(BsrError::SingularBlock { block_row: 0 })
    );
    let missing = Bsr::from_triplets(
        4,
        4,
        2,
        2,
        &[(0, 0, 1.0), (1, 1, 1.0), (0, 2, 1.0), (2, 0, 1.0)],
    )
    .unwrap();
    assert_eq!(
        missing.bilu0(),
        Err(BsrError::MissingDiagonal { block_row: 1 })
    );
}

#[test]
fn shape_validation() {
    let csr = Csr::from_triplets(5, 4, &[(0, 0, 1.0)]);
    assert_eq!(
        Bsr::from_csr(&csr, 2, 2),
        Err(BsrError::NotDivisible { extent: 5, block: 2 })
    );
    assert_eq!(Bsr::from_csr(&csr, 0, 2), Err(BsrError::ZeroBlock));
    let mut bsr = Bsr::from_csr(&Csr::from_triplets(4, 4, &[(0, 0, 1.0)]), 2, 2).unwrap();
    assert_eq!(bsr.validate(), Ok(()));
    bsr.indices[0] = 9;
    assert!(bsr.validate().is_err());
}

#[test]
fn validation_rejects_malformed_offsets_and_oversized_storage() {
    let mut a = Bsr::from_triplets(4, 4, 2, 2, &[(0, 0, 1.0), (3, 3, 2.0)]).unwrap();
    a.indptr = vec![0, usize::MAX, 2];
    assert!(a.validate().is_err());
    assert!(matches!(a.bilu0(), Err(BsrError::Invalid(_))));
    assert!(matches!(a.diagonal_blocks(), Err(BsrError::Invalid(_))));
    assert!(matches!(a.matmul(&a), Err(BsrError::Invalid(_))));
    assert!(matches!(
        a.matmul_dense(&[0.0; 4], 1),
        Err(BsrError::Invalid(_))
    ));

    let mut csr = Csr::from_triplets(2, 2, &[(0, 0, 1.0)]);
    csr.indptr[1] = usize::MAX;
    assert!(matches!(
        Bsr::from_csr(&csr, 1, 1),
        Err(BsrError::Invalid(_))
    ));
    assert_eq!(choose_block_size(&csr, &[(2, 2)]), (1, 1));
    csr.indptr[1] = 1;
    csr.indices[0] = 2;
    assert!(Bsr::from_csr(&csr, 1, 1).is_err());
    assert!(Bsr::zeros(usize::MAX, 0, 1, 1).is_err());
    assert!(Bsr::zeros(0, 0, usize::MAX, 2).is_err());
    assert!(Bsr::zeros(0, 0, 1, usize::MAX).is_ok());
    let a = Bsr::zeros(0, 0, usize::MAX, 1).unwrap();
    let b = Bsr::zeros(0, 0, 1, 2).unwrap();
    assert!(matches!(a.matmul(&b), Err(BsrError::Invalid(_))));
}

#[test]
fn sparse_storage_scales_with_touched_blocks_on_wide_shapes() {
    let cols = usize::MAX - 1;
    let trip = [
        (0, cols - 2, 2.0),
        (0, cols - 1, 3.0),
        (1, cols - 2, 5.0),
        (1, cols - 1, 7.0),
    ];
    let csr = Csr::from_triplets(2, cols, &trip);
    assert_eq!(choose_block_size(&csr, &[(2, 2), (usize::MAX, 2)]), (2, 2));
    let right = Bsr::from_csr(&csr, 2, 2).unwrap();
    assert_eq!(right, Bsr::from_triplets(2, cols, 2, 2, &trip).unwrap());
    assert_eq!(right.to_csr(), csr);

    let left = Bsr::from_triplets(2, 2, 2, 2, &[(0, 0, 1.0), (1, 0, 1.0), (1, 1, 1.0)]).unwrap();
    let c = left.matmul(&right).unwrap();
    assert_eq!(c.validate(), Ok(()));
    assert_eq!(c.indptr, [0, 1]);
    assert_eq!(c.indices, [cols / 2 - 1]);
    assert_eq!(c.values, [2.0, 3.0, 7.0, 10.0]);

    let tall = Bsr::from_triplets(usize::MAX, 1, usize::MAX, 1, &[]).unwrap();
    assert_eq!(tall.validate(), Ok(()));
    assert_eq!(tall.block_count(), 0);
}

#[test]
fn triplets_preserve_duplicate_order_and_drop_cancellations() {
    let trip = [
        (2, 2, 8.0),
        (0, 0, 1e16),
        (0, 2, 4.0),
        (0, 0, 1.0),
        (0, 0, -1e16),
        (0, 2, -5.0),
        (2, 3, 3.0),
        (2, 2, -8.0),
        (7, 8, 42.0),
    ];
    let via_csr = Bsr::from_csr(&Csr::from_triplets(4, 4, &trip), 2, 2).unwrap();
    let direct = Bsr::from_triplets(4, 4, 2, 2, &trip).unwrap();
    assert_eq!(direct, via_csr);
    assert_eq!(direct.indptr, [0, 1, 2]);
    assert_eq!(direct.indices, [1, 1]);
    assert_eq!(direct.values, [-1.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0]);
}

#[test]
fn padded_zero_blocks_do_not_touch_infinite_or_nan_rhs_entries() {
    let a = Bsr::from_triplets(2, 2, 2, 2, &[(0, 0, 1.0), (1, 1, 2.0)]).unwrap();
    let mut product = [0.0, 0.0, 91.0];
    let mut reference = [0.0; 2];
    let x = [f64::INFINITY, 3.0];
    a.matvec(&x, &mut product);
    a.to_csr().matvec(&x, &mut reference);
    assert_eq!(product[..2], reference);
    assert_eq!(product[2], 91.0);

    a.matvec_transpose(&[f64::NAN, 3.0], &mut product);
    assert!(product[0].is_nan());
    assert_eq!(product[1..], [6.0, 91.0]);
    let batch = a
        .matmul_dense(&[f64::INFINITY, 2.0, 3.0, f64::NAN], 2)
        .unwrap();
    assert_eq!(batch[..3], [f64::INFINITY, 2.0, 6.0]);
    assert!(batch[3].is_nan());
}

#[test]
fn dense_rhs_covers_exact_cancellation_shape_errors_and_empty_extents() {
    let a = Bsr::from_triplets(2, 3, 2, 3, &[(0, 0, 1e16), (0, 1, 1.0), (0, 2, -1e16)]).unwrap();
    let x = [1.0, 2.0, 1.0, 3.0, 1.0, 2.0];
    assert_eq!(a.matmul_dense(&x, 2).unwrap(), [1.0, 3.0, 0.0, 0.0]);
    assert!(matches!(
        a.matmul_dense(&x, 3),
        Err(BsrError::ShapeMismatch(_))
    ));
    assert!(a.matmul_dense(&[], 0).unwrap().is_empty());
    let no_input_rows = Bsr::zeros(2, 0, 2, 1).unwrap();
    assert_eq!(no_input_rows.matmul_dense(&[], 3).unwrap(), [0.0; 6]);
    assert!(no_input_rows.matmul_dense(&[], usize::MAX).is_err());
    let no_output_rows = Bsr::zeros(0, 3, 2, 3).unwrap();
    assert!(no_output_rows.matmul_dense(&x, 2).unwrap().is_empty());
}

#[test]
fn sparse_and_dense_batches_match_across_worker_counts() {
    let csr = block_tridiagonal(&mut Rng::new(202), 360, 8);
    let a = Bsr::from_csr(&csr, 8, 8).unwrap();
    let mut rng = Rng::new(21);
    let rhs_cols = 5;
    let x: Vec<_> = (0..a.ncols * rhs_cols).map(|_| rng.normal()).collect();
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| (a.matmul(&a).unwrap(), a.matmul_dense(&x, rhs_cols).unwrap()))
    };
    let (seq, dense) = run(1);
    let (par, parallel_dense) = run(4);
    assert_eq!(seq.indptr, par.indptr);
    assert_eq!(seq.indices, par.indices);
    assert!(
        seq.values
            .iter()
            .zip(&par.values)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    assert!(
        dense
            .iter()
            .zip(&parallel_dense)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );

    let mut product = vec![0.0; a.nrows];
    for rhs in 0..rhs_cols {
        let column: Vec<_> = x.iter().skip(rhs).step_by(rhs_cols).copied().collect();
        csr.matvec(&column, &mut product);
        for (row, expected) in product.iter().enumerate() {
            assert_eq!(dense[row * rhs_cols + rhs].to_bits(), expected.to_bits());
        }
    }
}

#[test]
fn factorization_workspace_resets_pivots_between_blocks() {
    let a = Bsr::from_triplets(
        4,
        4,
        2,
        2,
        &[
            (0, 1, 1.0),
            (1, 0, 1.0),
            (2, 2, 2.0),
            (2, 3, 3.0),
            (3, 2, 2.0),
            (3, 3, -1.0),
        ],
    )
    .unwrap();
    let rhs = [2.0, 1.0, 18.0, 2.0];
    let mut solved = [0.0; 4];
    a.block_jacobi().unwrap().apply(&rhs, &mut solved);
    assert_eq!(solved, [1.0, 2.0, 3.0, 4.0]);
    solved.fill(0.0);
    a.bilu0().unwrap().apply(&rhs, &mut solved);
    assert_eq!(solved, [1.0, 2.0, 3.0, 4.0]);

    let empty = Bsr::zeros(0, 0, 2, 2).unwrap();
    empty.block_jacobi().unwrap().apply(&[], &mut []);
    empty.bilu0().unwrap().apply(&[], &mut []);
}

#[test]
fn inverse_handles_finite_blocks_near_float_limits() {
    for scale in [1e-308, 1e308] {
        let a = Bsr::from_triplets(
            2,
            2,
            2,
            2,
            &[(0, 0, scale), (0, 1, scale), (1, 0, -scale), (1, 1, scale)],
        )
        .unwrap();
        let x = [0.125, -0.25];
        let mut rhs = [0.0; 2];
        a.matvec(&x, &mut rhs);
        let mut out = [0.0; 2];
        a.block_jacobi().unwrap().apply(&rhs, &mut out);
        assert!(
            out.iter()
                .zip(x)
                .all(|(actual, want)| (actual - want).abs() < 1e-15)
        );
        a.bilu0().unwrap().apply(&rhs, &mut out);
        assert!(
            out.iter()
                .zip(x)
                .all(|(actual, want)| (actual - want).abs() < 1e-15)
        );
    }
    let bad = Bsr::from_triplets(1, 1, 1, 1, &[(0, 0, f64::NAN)]).unwrap();
    assert!(matches!(
        bad.block_jacobi(),
        Err(BsrError::SingularBlock { block_row: 0 })
    ));
    assert!(matches!(
        bad.bilu0(),
        Err(BsrError::SingularBlock { block_row: 0 })
    ));
}

#[test]
fn sparse_multiply_drops_blocks_after_exact_cancellation() {
    let a = Bsr::from_triplets(2, 4, 2, 2, &[(0, 0, 1.0), (0, 2, 1.0), (1, 1, 2.0)]).unwrap();
    let b = Bsr::from_triplets(4, 2, 2, 2, &[(0, 0, 3.0), (2, 0, -3.0)]).unwrap();
    assert_eq!(a.matmul(&b).unwrap(), Bsr::zeros(2, 2, 2, 2).unwrap());
}

#[test]
fn empty_matrix() {
    let empty = Csr::from_triplets(0, 0, &[]);
    let bsr = Bsr::from_csr(&empty, 2, 2).unwrap();
    assert_eq!(bsr.block_count(), 0);
    assert_eq!(bsr.validate(), Ok(()));
    assert_eq!(bsr.to_csr(), empty);
    let zero = Bsr::from_csr(&Csr::from_triplets(4, 6, &[]), 2, 3).unwrap();
    assert_eq!(zero.block_count(), 0);
    assert!((zero.fill_ratio() - 1.0).abs() < f64::EPSILON);
    let mut y = vec![1.0; 4];
    zero.matvec(&[0.0; 6], &mut y);
    assert_eq!(y, vec![0.0; 4]);
    assert_eq!(zero.transpose().validate(), Ok(()));
}

#[test]
fn unit_blocks_are_csr() {
    let mut rng = Rng::new(8);
    let csr = random_csr(&mut rng, 10, 10, 0.3);
    let bsr = Bsr::from_csr(&csr, 1, 1).unwrap();
    assert_eq!(bsr.indptr, csr.indptr);
    assert_eq!(bsr.indices, csr.indices);
    assert_eq!(bsr.values, csr.values);
    assert!((bsr.fill_ratio() - 1.0).abs() < f64::EPSILON);
}

#[test]
fn block_size_selection() {
    // 2x2 nodal blocks on a chain: dense blocks, ratio 1 (4x4 would store 1.86x).
    let mut rng = Rng::new(9);
    let csr = block_tridiagonal(&mut rng, 10, 2);
    let cands = [(4, 4), (3, 3), (2, 2), (1, 1)];
    assert_eq!(choose_block_size(&csr, &cands), (2, 2));
    let bsr = Bsr::from_csr(&csr, 2, 2).unwrap();
    assert!((bsr.fill_ratio() - 1.0).abs() < 1e-12);
    // A scalar tridiagonal matrix gains nothing from blocks.
    let trip: Vec<_> = (0..12)
        .flat_map(|i| [(i, i, 2.0), (i, (i + 1) % 12, -1.0)])
        .collect();
    let scalar = Csr::from_triplets(12, 12, &trip);
    assert_eq!(choose_block_size(&scalar, &cands), (1, 1));
    assert_eq!(
        choose_block_size(&Csr::from_triplets(4, 4, &[]), &cands),
        (1, 1)
    );
}

/// `P (x) S` with the 2-D five-point Laplacian `P` on an `m x m` grid and a
/// dense SPD `s x s` nodal coupling `S`: scalar Jacobi cannot see `S`.
fn nodal_system(
    m: usize,
    s: usize,
    coupling: &[f64],
) -> Csr {
    let mut trip = Vec::new();
    for gy in 0..m {
        for gx in 0..m {
            let node = gy * m + gx;
            let mut nbrs = vec![(node, 4.0)];
            if gx > 0 {
                nbrs.push((node - 1, -1.0));
            }
            if gx + 1 < m {
                nbrs.push((node + 1, -1.0));
            }
            if gy > 0 {
                nbrs.push((node - m, -1.0));
            }
            if gy + 1 < m {
                nbrs.push((node + m, -1.0));
            }
            for (other, w) in nbrs {
                for a in 0..s {
                    for b in 0..s {
                        trip.push((node * s + a, other * s + b, w * coupling[a * s + b]));
                    }
                }
            }
        }
    }
    Csr::from_triplets(m * m * s, m * m * s, &trip)
}

fn cg_iterations(
    csr: &Csr,
    pre: &dyn Fn(&[f64], &mut [f64]),
) -> usize {
    let n = csr.rows;
    let b: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64).collect();
    let result = cg(|x, y| csr.matvec(x, y), pre, &b, None, 1e-8, 5000);
    assert!(result.converged);
    result.iterations
}

#[test]
fn block_jacobi_beats_scalar_jacobi_in_cg() {
    let couplings: [(usize, Vec<f64>); 2] = [
        (2, vec![1.0, 0.95, 0.95, 1.0]),
        (3, vec![1.0, 0.9, 0.9, 0.9, 1.0, 0.9, 0.9, 0.9, 1.0]),
    ];
    for (s, coupling) in couplings {
        let csr = nodal_system(12, s, &coupling);
        let inv = csr.jacobi_inverse();
        let scalar = cg_iterations(&csr, &|r, z| {
            for i in 0..r.len() {
                z[i] = r[i] * inv[i];
            }
        });
        let bsr = Bsr::from_csr(&csr, s, s).unwrap();
        assert!((bsr.fill_ratio() - 1.0).abs() < 1e-12);
        let pre = bsr.block_jacobi().unwrap();
        let block = cg_iterations(&csr, &|r, z| pre.apply(r, z));
        assert!(block < scalar, "s={s}: block {block} vs scalar {scalar}");
        let ilu = bsr.bilu0().unwrap();
        let bilu = cg_iterations(&csr, &|r, z| ilu.apply(r, z));
        assert!(bilu < block, "s={s}: bilu {bilu} vs block {block}");
    }
}

#[test]
fn block_jacobi_setup_and_apply_are_thread_count_independent() {
    // 1100 block rows of 4x4 blocks: the set-up and the apply both cross the
    // parallel work threshold (apply work is b^2 per block).
    let csr = block_tridiagonal(&mut Rng::new(404), 1100, 4);
    let a = Bsr::from_csr(&csr, 4, 4).unwrap();
    let rhs: Vec<f64> = (0..a.nrows).map(|i| 0.5 + (i % 11) as f64 * 0.25).collect();
    let run = |threads: usize| {
        crate::execution::with_scoped_pool(Some(threads), || {
            let pre = a.block_jacobi().unwrap();
            let mut z = vec![0.0; a.nrows];
            pre.apply(&rhs, &mut z);
            (pre, z)
        })
    };
    let (serial_pre, serial_z) = run(1);
    for threads in [2usize, 4] {
        let (parallel_pre, parallel_z) = run(threads);
        assert_eq!(serial_pre, parallel_pre);
        assert!(
            serial_z
                .iter()
                .zip(&parallel_z)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}
