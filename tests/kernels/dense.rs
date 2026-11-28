//! Tests for the dense linear algebra kernels.

use rssn::kernels::dense::DenseError;
use rssn::kernels::dense::Mat;
use rssn::kernels::dense::cholesky;
use rssn::kernels::dense::cholesky_solve;
use rssn::kernels::dense::cond;
use rssn::kernels::dense::eigen_symmetric;
use rssn::kernels::dense::eigenvalues;
use rssn::kernels::dense::lstsq;
use rssn::kernels::dense::lstsq_svd;
use rssn::kernels::dense::lu_factor;
use rssn::kernels::dense::pinv;
use rssn::kernels::dense::qr;
use rssn::kernels::dense::rank;
use rssn::kernels::dense::solve;
use rssn::kernels::dense::spectral_radius;
use rssn::kernels::dense::svd;

fn m(rows: &[&[f64]]) -> Mat {
    Mat::from_rows(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap()
}

#[test]
fn lu_solve_det_inverse() {
    let a = m(&[&[2.0, 1.0, 1.0], &[4.0, -6.0, 0.0], &[-2.0, 7.0, 2.0]]);
    let lu = lu_factor(&a).unwrap();
    assert!((lu.det() - (-16.0)).abs() < 1e-12);
    let x = lu.solve(&[5.0, -2.0, 9.0]);
    for (v, e) in x.iter().zip([1.0, 1.0, 2.0]) {
        assert!((v - e).abs() < 1e-12);
    }
    let prod = a.matmul(&lu.inverse()).unwrap();
    assert!((prod.at(0, 0) - 1.0).abs() < 1e-12 && prod.at(0, 1).abs() < 1e-12);
    assert!(lu_factor(&m(&[&[1.0, 2.0], &[2.0, 4.0]])).is_err());
}

#[test]
fn lu_reuses_caller_storage_and_validates_before_writing() {
    let a = m(&[&[2.0, 1.0, 1.0], &[4.0, -6.0, 0.0], &[-2.0, 7.0, 2.0]]);
    let lu = lu_factor(&a).unwrap();
    let b = [5.0, -2.0, 9.0];
    let mut out = [99.0; 3];
    lu.solve_into(&b, &mut out).unwrap();
    assert_eq!(out.as_slice(), lu.solve(&b));
    let expected = out;
    assert_eq!(lu.solve_into(&b[..2], &mut out), Err(DenseError::Shape));
    assert_eq!(
        lu.solve_into(&[f64::NAN, 0.0, 0.0], &mut out),
        Err(DenseError::NonFinite)
    );
    assert_eq!(out, expected);
    assert!(lu.solve(&b[..2]).is_empty());
    assert_eq!(
        solve(&a, &[0.0, f64::INFINITY, 0.0]),
        Err(DenseError::NonFinite)
    );
}

#[test]
fn dense_shapes_and_full_exponent_range() {
    let malformed = Mat {
        rows: 2,
        cols: 2,
        data: vec![1.0],
    };
    assert!(matches!(lu_factor(&malformed), Err(DenseError::Shape)));
    assert_eq!(malformed.matmul(&Mat::identity(2)), Err(DenseError::Shape));
    assert!(matches!(
        eigen_symmetric(&malformed),
        Err(DenseError::Shape)
    ));
    assert_eq!(lstsq_svd(&Mat::identity(2), &[1.0]), Err(DenseError::Shape));

    for factor in [f64::MAX, f64::from_bits(1)] {
        let a = m(&[&[factor]]);
        assert_eq!(a.norm_fro(), factor);
        assert_eq!(eigen_symmetric(&a).unwrap().0, [factor]);
        assert_eq!(eigenvalues(&a).unwrap(), [(factor, 0.0)]);
    }
    let a = m(&[&[1e308], &[1e308]]);
    let solution = lstsq_svd(&a, &[f64::MAX, f64::MAX]).unwrap();
    assert!((solution[0] - f64::MAX / 1e308).abs() < 1e-14);
}

#[test]
fn frobenius_norm_covers_vector_tails_and_rescaling() {
    for scale in [1.0, 1e300, 1e-300] {
        let mut entries = [0.0; 19];
        entries[0] = 3.0 * scale;
        entries[1] = 4.0 * scale;
        entries[18] = 12.0 * scale;
        assert!((m(&[&entries]).norm_fro() / scale - 13.0).abs() < 4e-15);
    }
    assert!(m(&[&[f64::NAN, f64::INFINITY]]).norm_fro().is_nan());
}

#[test]
fn cholesky_factor() {
    let a = m(&[
        &[4.0, 12.0, -16.0],
        &[12.0, 37.0, -43.0],
        &[-16.0, -43.0, 98.0],
    ]);
    let l = cholesky(&a).unwrap();
    assert_eq!((l.at(0, 0), l.at(1, 0), l.at(2, 0)), (2.0, 6.0, -8.0));
    assert!((l.at(2, 2) - 3.0).abs() < 1e-12);
    let x = cholesky_solve(&a, &[1.0, 2.0, 3.0]).unwrap();
    let back = a.matvec(&x);
    assert!((back[0] - 1.0).abs() < 1e-9 && (back[2] - 3.0).abs() < 1e-9);
    assert!(cholesky(&m(&[&[1.0, 2.0], &[2.0, 1.0]])).is_err());
}

#[test]
fn qr_and_least_squares() {
    let a = m(&[&[1.0, 1.0], &[1.0, 2.0], &[1.0, 3.0], &[1.0, 4.0]]);
    let (q, r) = qr(&a);
    let back = q.matmul(&r).unwrap();
    for i in 0..4 {
        for j in 0..2 {
            assert!((back.at(i, j) - a.at(i, j)).abs() < 1e-12);
        }
    }
    let y: Vec<f64> = (1..=4).map(|x| 2.0 + 3.0 * f64::from(x)).collect();
    let c = lstsq(&a, &y).unwrap();
    assert!((c[0] - 2.0).abs() < 1e-12 && (c[1] - 3.0).abs() < 1e-12);
    let c2 = lstsq_svd(&a, &y).unwrap();
    assert!((c2[0] - 2.0).abs() < 1e-10 && (c2[1] - 3.0).abs() < 1e-10);
}

#[test]
fn symmetric_eigen() {
    let a = m(&[&[2.0, -1.0, 0.0], &[-1.0, 2.0, -1.0], &[0.0, -1.0, 2.0]]);
    let (vals, vecs) = eigen_symmetric(&a).unwrap();
    let s2 = 2.0_f64.sqrt();
    for (v, e) in vals.iter().zip([2.0 - s2, 2.0, 2.0 + s2]) {
        assert!((v - e).abs() < 1e-12);
    }
    for k in 0..3 {
        let v: Vec<f64> = (0..3).map(|i| vecs.at(i, k)).collect();
        let av = a.matvec(&v);
        for i in 0..3 {
            assert!((av[i] - vals[k] * v[i]).abs() < 1e-12);
        }
    }
}

#[test]
fn svd_cond_rank_pinv() {
    let a = m(&[&[3.0, 0.0], &[0.0, 1.0]]);
    let s = svd(&a).unwrap();
    assert!((s.s[0] - 3.0).abs() < 1e-12 && (s.s[1] - 1.0).abs() < 1e-12);
    assert!((cond(&a).unwrap() - 3.0).abs() < 1e-12);
    let h = Mat {
        rows: 4,
        cols: 4,
        data: (0..16)
            .map(|k| 1.0 / f64::from(k / 4 + k % 4 + 1))
            .collect(),
    };
    assert!((cond(&h).unwrap() - 15_513.738_738).abs() < 1e-2);
    let r = m(&[&[1.0, 2.0], &[2.0, 4.0], &[3.0, 6.0]]);
    assert_eq!(rank(&r, None).unwrap(), 1);
    let p = pinv(&r).unwrap();
    let rp = r.matmul(&p).unwrap().matmul(&r).unwrap();
    for (x, y) in rp.data.iter().zip(&r.data) {
        assert!((x - y).abs() < 1e-12);
    }
    let w = m(&[&[1.0, 0.0, 0.0], &[0.0, 2.0, 0.0]]);
    let sw = svd(&w).unwrap();
    assert!((sw.s[0] - 2.0).abs() < 1e-12);
}

#[test]
fn nonsymmetric_eigenvalues() {
    let rot = m(&[&[0.0, -1.0], &[1.0, 0.0]]);
    let mut ev = eigenvalues(&rot).unwrap();
    ev.sort_by(|a, b| a.1.total_cmp(&b.1));
    assert!(ev[0].0.abs() < 1e-14 && (ev[0].1 + 1.0).abs() < 1e-14);
    assert!((ev[1].1 - 1.0).abs() < 1e-14);
    let a = m(&[
        &[1.0, 2.0, 3.0, 4.0],
        &[0.0, 2.0, 5.0, 6.0],
        &[0.0, 0.0, 3.0, 7.0],
        &[0.0, 0.0, 0.0, 4.0],
    ]);
    let mut re: Vec<f64> = eigenvalues(&a).unwrap().iter().map(|e| e.0).collect();
    re.sort_by(f64::total_cmp);
    for (v, e) in re.iter().zip([1.0, 2.0, 3.0, 4.0]) {
        assert!((v - e).abs() < 1e-10);
    }
    let c = m(&[
        &[0.0, 0.0, 0.0, 0.0, 120.0],
        &[1.0, 0.0, 0.0, 0.0, -274.0],
        &[0.0, 1.0, 0.0, 0.0, 225.0],
        &[0.0, 0.0, 1.0, 0.0, -85.0],
        &[0.0, 0.0, 0.0, 1.0, 15.0],
    ]);
    let mut re: Vec<f64> = eigenvalues(&c).unwrap().iter().map(|e| e.0).collect();
    re.sort_by(f64::total_cmp);
    for (v, e) in re.iter().zip([1.0, 2.0, 3.0, 4.0, 5.0]) {
        assert!((v - e).abs() < 1e-8, "{re:?}");
    }
    assert!((spectral_radius(&c).unwrap() - 5.0).abs() < 1e-8);
}

#[test]
fn solve_shape_error() {
    let a = m(&[&[1.0, 0.0], &[0.0, 1.0]]);
    assert!(solve(&a, &[1.0]).is_err());
}

#[test]
fn generalized_symmetric_eigenproblem_with_positive_and_indefinite_b() {
    use rssn::kernels::dense::eigen_generalized_symmetric;
    // SPD B: pencil (A, B) with A = diag(2, 3, 5) in B-coordinates, A = S^T D S, B = S^T S.
    let s = Mat::from_rows(&[
        vec![1.0, 0.5, 0.0],
        vec![0.0, 1.0, 0.25],
        vec![0.3, 0.0, 1.0],
    ])
    .unwrap();
    let d = [2.0, 3.0, 5.0];
    let (st, mut dm) = (s.transpose(), Mat::zeros(3, 3));
    for i in 0..3 {
        dm.set(i, i, d[i]);
    }
    let b = st.matmul(&s).unwrap();
    let a = st.matmul(&dm).unwrap().matmul(&s).unwrap();
    let r = eigen_generalized_symmetric(&a, &b).unwrap();
    assert!(r.b_positive_definite);
    for (got, want) in r.values.iter().zip(d) {
        assert!(
            (got.0 - want).abs() < 1e-10 && got.1 == 0.0,
            "{got:?} {want}"
        );
    }
    let x = r.vectors.unwrap();
    let xtbx = x.transpose().matmul(&b).unwrap().matmul(&x).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((xtbx.at(i, j) - want).abs() < 1e-10);
        }
    }
    // indefinite B, real spectrum: A = diag(2, 3), B = diag(1, -1) -> {2, -3}
    let a = Mat::from_rows(&[vec![2.0, 0.0], vec![0.0, 3.0]]).unwrap();
    let b = Mat::from_rows(&[vec![1.0, 0.0], vec![0.0, -1.0]]).unwrap();
    let r = eigen_generalized_symmetric(&a, &b).unwrap();
    assert!(!r.b_positive_definite && r.vectors.is_none());
    let mut v: Vec<f64> = r.values.iter().map(|e| e.0).collect();
    v.sort_by(f64::total_cmp);
    assert!((v[0] + 3.0).abs() < 1e-12 && (v[1] - 2.0).abs() < 1e-12);
    // indefinite B, complex pair: A = [[0, 1], [1, 0]], B = diag(1, -1): B^-1 A = [[0, 1], [-1, 0]]
    let a = Mat::from_rows(&[vec![0.0, 1.0], vec![1.0, 0.0]]).unwrap();
    let r = eigen_generalized_symmetric(&a, &b).unwrap();
    let mut im: Vec<f64> = r.values.iter().map(|e| e.1).collect();
    im.sort_by(f64::total_cmp);
    assert!(r.values.iter().all(|e| e.0.abs() < 1e-12));
    assert!((im[0] + 1.0).abs() < 1e-12 && (im[1] - 1.0).abs() < 1e-12);
    // shapes
    assert!(eigen_generalized_symmetric(&a, &Mat::identity(3)).is_err());
}

#[test]
fn optimized_matmul_preserves_rectangular_layout_and_validates_empty_products() {
    use rssn::kernels::dense::optimized::matmul;

    let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let b = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
    assert_eq!(matmul(2, 3, 2, &a, &b), Ok(vec![58.0, 64.0, 139.0, 154.0]));
    assert_eq!(matmul(3, 0, 2, &[], &[]), Ok(vec![0.0; 6]));
    assert_eq!(
        matmul(0, 2, 1, &[], &[1.0, f64::NAN]),
        Err(DenseError::NonFinite)
    );
    assert_eq!(matmul(2, 3, 2, &a[..5], &b), Err(DenseError::Shape));
    assert_eq!(matmul(usize::MAX, 0, 1, &[], &[]), Err(DenseError::Shape));
}

#[test]
fn optimized_matmul_uses_the_active_pool_without_changing_dot_products() {
    use rssn::kernels::dense::optimized::matmul;

    let (rows, inner, cols) = (160, 210, 144);
    let lhs: Vec<_> = (0..rows * inner)
        .map(|i| (i % 31) as f64 * 0.13 - 1.0)
        .collect();
    let rhs: Vec<_> = (0..inner * cols)
        .map(|i| (i % 27) as f64 * 0.41 - 3.0)
        .collect();
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| matmul(rows, inner, cols, &lhs, &rhs).unwrap())
    };
    let sequential = run(1);
    let parallel = run(4);
    assert_eq!(sequential.len(), rows * cols);
    assert!(
        sequential
            .iter()
            .zip(&parallel)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

mod source_leaf_contracts {
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::One;
    use num_traits::Zero;
    use rssn::kernels::dense;
    use rssn::kernels::dense::exact;
    use rssn::kernels::dense::optimized;
    use rssn::kernels::dense::reference;

    fn rat(
        num: i64,
        den: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(num), BigInt::from(den))
    }

    fn multiply(
        a: &[BigRational],
        b: &[BigRational],
        n: usize,
    ) -> Vec<BigRational> {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| (0..n).map(|k| &a[i * n + k] * &b[k * n + j]).sum()))
            .collect()
    }

    #[test]
    fn exact_lu_reconstructs_singular_rational_matrices() {
        let entries = [
            rat(0, 1),
            rat(1, 3),
            rat(2, 3),
            rat(2, 5),
            rat(3, 7),
            rat(5, 9),
            rat(0, 1),
            rat(1, 3),
            rat(2, 3),
        ];
        let (p, l, u) = exact::lu_exact(3, &entries).unwrap();
        assert_eq!(multiply(&p, &entries, 3), multiply(&l, &u, 3));
        for i in 0..3 {
            assert_eq!(l[i * 3 + i], BigRational::one());
            for j in 0..i {
                assert!(u[i * 3 + j].is_zero());
                assert!(l[j * 3 + i].is_zero());
            }
        }
        assert!(u[8].is_zero());
        assert_eq!(exact::lu_exact(0, &[]).unwrap(), (vec![], vec![], vec![]));
        assert_eq!(
            exact::lu_exact(usize::MAX, &[]),
            Err(exact::ExactError::Shape)
        );
    }

    #[test]
    fn exact_lu_shares_independent_rational_rows_between_workers() {
        let n = 96;
        let a: Vec<_> = (0..n)
            .flat_map(|row| {
                (0..n).map(move |col| {
                    if row == col {
                        BigRational::one()
                    } else if col == 0 || (row == 0 && col >= n - 5) {
                        rat(1, 3)
                    } else {
                        BigRational::zero()
                    }
                })
            })
            .collect();
        let factor = |workers| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap()
                .install(|| exact::lu_exact(n, &a).unwrap())
        };
        let serial = factor(1);
        let parallel = factor(4);
        assert_eq!(serial, parallel);
        let (p, l, u) = parallel;
        for (row, col) in [(0, 0), (15, 0), (15, n - 1), (n - 1, n - 1)] {
            let pa: BigRational = (0..n).map(|i| &p[row * n + i] * &a[i * n + col]).sum();
            let lu: BigRational = (0..n).map(|i| &l[row * n + i] * &u[i * n + col]).sum();
            assert_eq!(pa, lu);
        }
    }

    #[test]
    fn exact_qr_reconstructs_rank_deficiency_and_reports_irrational_norms() {
        let a = [3, 6, 0, 4, 8, 0, 0, 0, 5].map(|value| rat(value, 1));
        let (q, r) = exact::qr_exact(3, 3, &a).unwrap();
        assert_eq!(multiply(&q, &r, 3), a);
        let flattened: Vec<_> = (0..9).map(|idx| q[idx % 3 * 3 + idx / 3].clone()).collect();
        assert_eq!(
            multiply(&flattened, &q, 3),
            [1, 0, 0, 0, 1, 0, 0, 0, 0].map(|value| rat(value, 1))
        );
        assert_eq!(
            exact::qr_exact(2, 1, &[rat(1, 1), rat(1, 1)]),
            Err(exact::ExactError::IrrationalNorm)
        );
        assert_eq!(exact::exact_sqrt(&rat(49, 225)), Some(rat(7, 15)));
        assert_eq!(exact::exact_sqrt(&rat(-1, 1)), None);
        assert_eq!(exact::exact_sqrt(&rat(13, 7)), None);
    }

    #[test]
    fn exact_qr_zero_rows_do_not_scan_columns() {
        assert_eq!(
            exact::qr_exact(0, usize::MAX, &[]).unwrap(),
            (vec![], vec![])
        );
        assert_eq!(
            exact::qr_exact(2, 0, &[]).unwrap(),
            (vec![BigRational::zero(); 4], vec![])
        );
    }

    #[test]
    fn reference_and_faer_solve_and_det_preserve_domain_and_row_layout() {
        let a = [0.0, 2.0, 3.0, -1.0, 4.0, 1.0, 5.0, 5.0, 1.0];
        let rhs = [5.0, -6.0, -4.0];
        for solve in [reference::solve, optimized::solve] {
            let x = solve(3, 3, &a, &rhs).unwrap();
            for (row, want) in a.as_chunks::<3>().0.iter().zip(rhs) {
                let got = row.iter().zip(&x).map(|(a, b)| a * b).sum::<f64>();
                assert!((got - want).abs() < 1e-12);
            }
            assert_eq!(
                solve(1, 1, &[0.0], &[1.0]),
                Err(dense::DenseError::Singular)
            );
            assert_eq!(
                solve(1, 1, &[1.0], &[f64::NAN]),
                Err(dense::DenseError::NonFinite)
            );
        }
        let d = reference::det(3, 3, &a).unwrap();
        assert!((optimized::det(3, 3, &a).unwrap() - d).abs() < 1e-12);
        assert_eq!(reference::det(0, 0, &[]), Ok(1.0));
        assert_eq!(optimized::det(0, 0, &[]), Ok(1.0));
        assert_eq!(
            reference::det(1, 2, &[1.0, 2.0]),
            Err(dense::DenseError::Shape)
        );
        assert_eq!(
            optimized::det(1, 2, &[1.0, 2.0]),
            Err(dense::DenseError::Shape)
        );
    }

    #[test]
    fn both_qr_interfaces_reconstruct_a_tall_input_on_their_shared_columns() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0];
        type Factor = fn(usize, usize, &[f64]) -> Result<(Vec<f64>, Vec<f64>), dense::DenseError>;
        let methods: [(Factor, usize); 2] = [(reference::qr, 3), (optimized::qr, 2)];
        for (qr, r_rows) in methods {
            let (q, r) = qr(3, 2, &a).unwrap();
            assert_eq!(r.len(), r_rows * 2);
            let q_prefix: Vec<_> = q
                .as_chunks::<3>()
                .0
                .iter()
                .flat_map(|row| row[..r_rows].iter().copied())
                .collect();
            let rebuilt = reference::matmul(3, r_rows, 2, &q_prefix, &r);
            assert!(
                rebuilt
                    .iter()
                    .zip(a)
                    .all(|(got, want)| (got - want).abs() < 1e-12)
            );
            for i in 0..2 {
                let norm = (0..3)
                    .map(|row| q[row * 3 + i] * q[row * 3 + i])
                    .sum::<f64>();
                assert!((norm - 1.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn reference_products_keep_scalar_order_and_single_product_zero_sign() {
        let a = [1e16, 1.0, -1e16, -0.0, -0.0, -0.0];
        let b = [1.0, 5.0, 1.0, 6.0, 1.0, 5.0];
        let result = reference::matmul(2, 3, 2, &a, &b);
        let expected = [0.0_f64, 8.0, -0.0, -0.0];
        for (&got, &want) in result.iter().zip(&expected) {
            assert_eq!(got.to_bits(), want.to_bits());
        }
        let wide_b = [
            1.0, 5.0, 1.0, 5.0, 1.0, 5.0, 1.0, 6.0, 1.0, 6.0, 1.0, 6.0, 1.0, 5.0, 1.0, 5.0, 1.0,
            5.0,
        ];
        let wide = reference::matmul(2, 3, 6, &a, &wide_b);
        for (got, want) in wide.into_iter().zip([
            0.0_f64, 8.0, 0.0, 8.0, 0.0, 8.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
        ]) {
            assert_eq!(got.to_bits(), want.to_bits());
        }
        let one = reference::matmul(1, 1, 2, &[-0.0], &[1.0, -1.0]);
        assert_eq!(one[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(one[1].to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            reference::matmul(2, 0, 2, &[], &[])
                .into_iter()
                .map(f64::to_bits)
                .collect::<Vec<_>>(),
            [0_u64; 4]
        );
    }
}
