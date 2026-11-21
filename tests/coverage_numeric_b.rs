//! Coverage tests: block sparse rows (BSR), exact sparse linear algebra and
//! the dense rational helpers. Public API only; independent implementations
//! (dense rational, sparse rational, modular, black-box) are cross-checked.

#![allow(clippy::needless_range_loop)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::bsr::Bsr;
use rssn::kernels::bsr::BsrError;
use rssn::kernels::bsr::choose_block_size;
use rssn::kernels::krylov::Csr;
use rssn::kernels::qlinalg::Q;
use rssn::kernels::qlinalg::QMat;
use rssn::kernels::qlinalg::q;
use rssn::kernels::qlinalg::{
    self,
};
use rssn::kernels::sparse_exact::Gfp;
use rssn::kernels::sparse_exact::Rationals;
use rssn::kernels::sparse_exact::Solution;
use rssn::kernels::sparse_exact::SparseMat;
use rssn::kernels::sparse_exact::berlekamp_massey;
use rssn::kernels::sparse_exact::det;
use rssn::kernels::sparse_exact::det_rational;
use rssn::kernels::sparse_exact::is_prime_u64;
use rssn::kernels::sparse_exact::minpoly_wiedemann;
use rssn::kernels::sparse_exact::nullspace;
use rssn::kernels::sparse_exact::rank;
use rssn::kernels::sparse_exact::rational_reconstruct;
use rssn::kernels::sparse_exact::solve;
use rssn::kernels::sparse_exact::solve_dixon;
use rssn::kernels::sparse_exact::solve_wiedemann;
use rssn::kernels::sparse_exact::word_primes;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn unit(&mut self) -> f64 {
        self.next() as f64 / (1u64 << 31) as f64
    }

    fn int(
        &mut self,
        lo: i64,
        hi: i64,
    ) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }
}

// --------------------------------------------------------------------- BSR

fn csr_dense(c: &Csr) -> Vec<Vec<f64>> {
    let mut d = vec![vec![0.0; c.cols]; c.rows];
    for i in 0..c.rows {
        for k in c.indptr[i]..c.indptr[i + 1] {
            d[i][c.indices[k]] += c.values[k];
        }
    }
    d
}

/// Random matrix whose nonzeros come in full `b x b` blocks.
fn block_triplets(
    nb: usize,
    b: usize,
    density: f64,
    seed: u64,
    dominant: bool,
) -> Vec<(usize, usize, f64)> {
    let mut rng = Lcg(seed);
    let mut t = Vec::new();
    for ib in 0..nb {
        for jb in 0..nb {
            let on_band = ib.abs_diff(jb) <= 1;
            let present = ib == jb || (on_band && rng.unit() < density);
            if !present {
                continue;
            }
            for r in 0..b {
                for c in 0..b {
                    let mut v = rng.unit() * 2.0 - 1.0 + 0.1;
                    if dominant && ib == jb && r == c {
                        v += 8.0;
                    }
                    t.push((ib * b + r, jb * b + c, v));
                }
            }
        }
    }
    t
}

#[test]
fn bsr_round_trip_and_matvec_are_exact() {
    let t = block_triplets(6, 3, 0.7, 1, false);
    let csr = Csr::from_triplets(18, 18, &t);
    let mut rng = Lcg(2);
    let x: Vec<f64> = (0..18).map(|_| rng.unit() - 0.5).collect();
    let mut y_ref = vec![0.0; 18];
    csr.matvec(&x, &mut y_ref);
    for (br, bc) in [(1, 1), (3, 3), (2, 3), (3, 2), (6, 9), (18, 18)] {
        let b = Bsr::from_csr(&csr, br, bc).unwrap();
        b.validate().unwrap();
        assert_eq!(b.to_csr(), csr, "block {br}x{bc}");
        assert!(b.fill_ratio() >= 1.0);
        let mut y = vec![f64::NAN; 18];
        b.matvec(&x, &mut y);
        for i in 0..18 {
            assert_eq!(
                y[i].to_bits(),
                y_ref[i].to_bits(),
                "block {br}x{bc} row {i}"
            );
        }
        // A^T x against the dense matrix
        let d = csr_dense(&csr);
        let mut yt = vec![f64::NAN; 18];
        b.matvec_transpose(&x, &mut yt);
        for j in 0..18 {
            let want: f64 = (0..18).map(|i| d[i][j] * x[i]).sum();
            assert!((yt[j] - want).abs() < 1e-12, "block {br}x{bc} col {j}");
        }
    }
    // full blocks: no padding at the natural block size
    let natural = Bsr::from_csr(&csr, 3, 3).unwrap();
    assert!((natural.fill_ratio() - 1.0).abs() < 1e-12);
    assert_eq!(natural.nnz(), csr.values.len());
    // from_triplets sums duplicates
    let dup = Bsr::from_triplets(4, 4, 2, 2, &[(0, 0, 1.0), (0, 0, 2.0), (3, 3, 5.0)]).unwrap();
    let mut y = vec![0.0; 4];
    dup.matvec(&[1.0, 1.0, 1.0, 1.0], &mut y);
    assert_eq!(y, vec![3.0, 0.0, 0.0, 5.0]);
}

#[test]
fn bsr_transpose_and_matmul_match_dense() {
    let a = Bsr::from_triplets(12, 8, 3, 2, &{
        let mut rng = Lcg(5);
        let mut t = Vec::new();
        for i in 0..12 {
            for j in 0..8 {
                if rng.unit() < 0.4 {
                    t.push((i, j, rng.unit() + 0.25));
                }
            }
        }
        t
    })
    .unwrap();
    let at = a.transpose();
    assert_eq!((at.nrows, at.ncols, at.br, at.bc), (8, 12, 2, 3));
    at.validate().unwrap();
    assert_eq!(at.transpose(), a);
    let da = csr_dense(&a.to_csr());
    let dat = csr_dense(&at.to_csr());
    for i in 0..12 {
        for j in 0..8 {
            assert_eq!(da[i][j], dat[j][i]);
        }
    }
    // A^T A (8x8) vs dense
    let prod = at.matmul(&a).unwrap();
    prod.validate().unwrap();
    let dp = csr_dense(&prod.to_csr());
    for i in 0..8 {
        for j in 0..8 {
            let want: f64 = (0..12).map(|k| da[k][i] * da[k][j]).sum();
            assert!((dp[i][j] - want).abs() < 1e-12, "({i},{j})");
        }
    }
    // incompatible inner extents / block sizes
    assert!(matches!(a.matmul(&a), Err(BsrError::ShapeMismatch(_))));
    let other_block = Bsr::zeros(8, 8, 4, 4).unwrap();
    assert!(matches!(
        a.matmul(&other_block),
        Err(BsrError::ShapeMismatch(_))
    ));
}

#[test]
fn bsr_constructor_errors() {
    assert_eq!(Bsr::zeros(4, 4, 0, 2).unwrap_err(), BsrError::ZeroBlock);
    assert_eq!(Bsr::zeros(4, 4, 2, 0).unwrap_err(), BsrError::ZeroBlock);
    assert!(matches!(
        Bsr::zeros(5, 4, 2, 2),
        Err(BsrError::NotDivisible { extent: 5, block: 2 })
    ));
    assert!(matches!(
        Bsr::zeros(4, 5, 2, 2),
        Err(BsrError::NotDivisible { .. })
    ));
    let z = Bsr::zeros(6, 6, 2, 3).unwrap();
    assert_eq!(z.nnz(), 0);
    let mut y = vec![9.0; 6];
    z.matvec(&[1.0; 6], &mut y);
    assert_eq!(y, vec![0.0; 6]);
    let csr = Csr::from_triplets(5, 5, &[(0, 0, 1.0)]);
    assert!(Bsr::from_csr(&csr, 2, 2).is_err());
    // 0 x 0 does not panic
    let _ = Bsr::zeros(0, 0, 1, 1);
    // display of every error is non-empty
    for e in [
        BsrError::ZeroBlock,
        BsrError::NotSquare,
        BsrError::MissingDiagonal { block_row: 1 },
        BsrError::SingularBlock { block_row: 2 },
        BsrError::ShapeMismatch("x".into()),
        BsrError::Invalid("y".into()),
    ] {
        assert!(!e.to_string().is_empty());
    }
}

#[test]
fn bsr_block_preconditioners_are_exact_on_their_pattern() {
    let (nb, b) = (7, 3);
    let n = nb * b;
    let mut rng = Lcg(11);
    let r: Vec<f64> = (0..n).map(|_| rng.unit() - 0.5).collect();
    // block diagonal: block Jacobi is the exact inverse
    let diag_only: Vec<_> = block_triplets(nb, b, 0.0, 3, true);
    let bd = Bsr::from_triplets(n, n, b, b, &diag_only).unwrap();
    let bj = bd.block_jacobi().unwrap();
    assert_eq!(bj.block_size(), b);
    let mut z = vec![0.0; n];
    bj.apply(&r, &mut z);
    let mut back = vec![0.0; n];
    bd.matvec(&z, &mut back);
    for i in 0..n {
        assert!((back[i] - r[i]).abs() < 1e-12);
    }
    // diagonal blocks accessor
    let blocks = bd.diagonal_blocks().unwrap();
    assert_eq!(blocks.len(), nb);
    assert!(blocks.iter().all(|blk| blk.len() == b * b));
    // block tridiagonal: BILU(0) has no dropped fill, so it is an exact solver
    let tri = block_triplets(nb, b, 1.0, 4, true);
    let a = Bsr::from_triplets(n, n, b, b, &tri).unwrap();
    let ilu = a.bilu0().unwrap();
    let mut z = vec![0.0; n];
    ilu.apply(&r, &mut z);
    let mut back = vec![0.0; n];
    a.matvec(&z, &mut back);
    for i in 0..n {
        assert!(
            (back[i] - r[i]).abs() < 1e-10,
            "row {i}: {} vs {}",
            back[i],
            r[i]
        );
    }
}

#[test]
fn bsr_preconditioner_errors() {
    // non-square
    let rect = Bsr::zeros(4, 6, 2, 2).unwrap();
    assert_eq!(rect.block_jacobi().unwrap_err(), BsrError::NotSquare);
    assert_eq!(rect.bilu0().unwrap_err(), BsrError::NotSquare);
    // non-square blocks
    let nsq = Bsr::zeros(6, 6, 2, 3).unwrap();
    assert_eq!(nsq.bilu0().unwrap_err(), BsrError::NotSquare);
    // block (1,1) has no diagonal block
    let missing = Bsr::from_triplets(
        4,
        4,
        2,
        2,
        &[
            (0, 0, 1.0),
            (0, 1, 1.0),
            (0, 3, 1.0),
            (2, 0, 1.0),
            (1, 1, 2.0),
            (3, 3, 0.0),
        ],
    );
    if let Ok(m) = missing {
        assert!(matches!(m.bilu0(), Err(BsrError::MissingDiagonal { .. })));
    }
    // singular diagonal block [[1,2],[2,4]]
    let sing = Bsr::from_triplets(
        4,
        4,
        2,
        2,
        &[
            (0, 0, 1.0),
            (0, 1, 2.0),
            (1, 0, 2.0),
            (1, 1, 4.0),
            (2, 2, 1.0),
            (3, 3, 1.0),
        ],
    )
    .unwrap();
    assert!(matches!(
        sing.block_jacobi(),
        Err(BsrError::SingularBlock { block_row: 0 })
    ));
    assert!(sing.bilu0().is_err());
}

#[test]
fn choose_block_size_picks_the_natural_block() {
    let t = block_triplets(8, 4, 0.5, 9, false);
    let csr = Csr::from_triplets(32, 32, &t);
    let cand = [(1, 1), (2, 2), (4, 4), (8, 8), (3, 3), (0, 4)];
    let (br, bc) = choose_block_size(&csr, &cand);
    assert_eq!((br, bc), (4, 4));
    // a pure diagonal matrix gains nothing from blocks
    let diag: Vec<_> = (0..32).map(|i| (i, i, 1.0 + i as f64)).collect();
    assert_eq!(
        choose_block_size(&Csr::from_triplets(32, 32, &diag), &cand),
        (1, 1)
    );
    // no nonzeros, or no candidates
    assert_eq!(
        choose_block_size(&Csr::from_triplets(32, 32, &[]), &cand),
        (1, 1)
    );
    assert_eq!(choose_block_size(&csr, &[]), (1, 1));
    // candidates that do not divide the shape are ignored
    assert_eq!(choose_block_size(&csr, &[(5, 5), (7, 3)]), (1, 1));
}

// ------------------------------------------------------------- qlinalg

fn qmat(rows: &[&[i64]]) -> QMat {
    rows.iter()
        .map(|r| r.iter().map(|&v| q(v)).collect())
        .collect()
}

#[test]
fn qlinalg_basic_identities() {
    let a = qmat(&[&[2, 1, 0], &[1, 3, 1], &[0, 1, 4]]);
    assert_eq!(qlinalg::det(&a), q(2 * 11 - 4));
    let inv = qlinalg::inverse(&a).unwrap();
    assert_eq!(qlinalg::matmul(&a, &inv), qlinalg::identity(3));
    assert_eq!(qlinalg::matmul(&inv, &a), qlinalg::identity(3));
    // det multiplicative, transpose invariant
    let b = qmat(&[&[1, 2, 3], &[0, 1, 4], &[5, 6, 0]]);
    assert_eq!(
        qlinalg::det(&qlinalg::matmul(&a, &b)),
        qlinalg::det(&a) * qlinalg::det(&b)
    );
    assert_eq!(qlinalg::det(&qlinalg::transpose(&b, 3)), qlinalg::det(&b));
    assert_eq!(qlinalg::det(&b), q(1));
    // singular
    let s = qmat(&[&[1, 2], &[2, 4]]);
    assert_eq!(qlinalg::det(&s), q(0));
    assert!(qlinalg::inverse(&s).is_none());
    assert_eq!(qlinalg::det(&[]), q(1));
    // matvec
    assert_eq!(
        qlinalg::matvec(&a, &[q(1), q(0), q(-1)]),
        vec![q(2), q(0), q(-4)]
    );
    assert!(qlinalg::is_zero_vec(&qlinalg::zeros(4)));
    assert!(!qlinalg::is_zero_vec(&[q(0), q(1)]));
}

#[test]
fn qlinalg_rref_rank_nullspace_solve() {
    let mut rng = Lcg(77);
    for _ in 0..30 {
        let (m, n) = (rng.int(1, 6) as usize, rng.int(1, 6) as usize);
        let a: QMat = (0..m)
            .map(|_| (0..n).map(|_| q(rng.int(-2, 2))).collect())
            .collect();
        let r = qlinalg::rank(&a, n);
        let ns = qlinalg::nullspace(&a, n);
        assert_eq!(ns.len(), n - r, "rank-nullity");
        for v in &ns {
            assert!(qlinalg::is_zero_vec(&qlinalg::matvec(&a, v)));
        }
        // rref pivots and rank agree
        let mut red = a.clone();
        assert_eq!(qlinalg::rref(&mut red, n).len(), r);
        assert_eq!(red.len(), r, "zero rows dropped");
        // consistent system A x = A x0
        let x0: Vec<Q> = (0..n).map(|_| q(rng.int(-3, 3))).collect();
        let b = qlinalg::matvec(&a, &x0);
        let x = qlinalg::solve(&a, n, &b).expect("consistent");
        assert_eq!(qlinalg::matvec(&a, &x), b);
        // span basis and reduce_mod: every row reduces to zero
        let mut basis = a.clone();
        let pivots = qlinalg::rref(&mut basis, n);
        assert_eq!(qlinalg::span_basis(&a, n), basis);
        for row in &a {
            assert!(qlinalg::is_zero_vec(&qlinalg::reduce_mod(
                row, &basis, &pivots
            )));
        }
    }
    // inconsistent
    let a = qmat(&[&[1, 1], &[1, 1]]);
    assert!(qlinalg::solve(&a, 2, &[q(1), q(2)]).is_none());
}

// ------------------------------------------------------------ sparse_exact

fn to_sparse(a: &QMat) -> SparseMat<BigRational> {
    SparseMat::from_dense(&mut Rationals, a)
}

fn random_int_mat(
    n: usize,
    density: u64,
    seed: u64,
    dominant: bool,
) -> QMat {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j && dominant {
                        q(n as i64 * 3 + rng.int(0, 3))
                    } else if rng.next() % 100 < density {
                        q(rng.int(-3, 3))
                    } else {
                        q(0)
                    }
                })
                .collect()
        })
        .collect()
}

#[test]
fn sparse_exact_determinants_agree_across_methods() {
    for seed in 0..25u64 {
        let n = 2 + (seed % 6) as usize;
        let a = random_int_mat(n, 45, seed, false);
        let want = qlinalg::det(&a);
        let sp = to_sparse(&a);
        assert_eq!(det(&mut Rationals, &sp), Some(want.clone()), "seed {seed}");
        assert_eq!(
            det_rational(&sp),
            Some(want.clone()),
            "modular, seed {seed}"
        );
        assert_eq!(
            rank(&mut Rationals, &sp),
            qlinalg::rank(&a, n),
            "seed {seed}"
        );
        // integer determinant reduced mod a prime equals the GF(p) elimination
        let p = Gfp(1_000_000_007);
        let rows: Vec<Vec<u64>> = a
            .iter()
            .map(|r| r.iter().map(|v| p.reduce(&v.to_integer())).collect())
            .collect();
        let spp = SparseMat::from_dense(&mut p.clone(), &rows);
        assert_eq!(
            det(&mut p.clone(), &spp),
            Some(p.reduce(&want.to_integer())),
            "mod p, seed {seed}"
        );
    }
    // non-square and empty
    let rect = SparseMat::<BigRational>::zeros(2, 3);
    assert!(det(&mut Rationals, &rect).is_none());
    assert!(det_rational(&rect).is_none());
    assert_eq!(
        det_rational(&SparseMat::<BigRational>::zeros(0, 0)),
        Some(BigRational::one())
    );
    assert_eq!(
        det_rational(&SparseMat::<BigRational>::zeros(3, 3)),
        Some(BigRational::zero())
    );
    // rational (non-integer) entries
    let h = vec![
        vec![
            BigRational::new(1.into(), 2.into()),
            BigRational::new(1.into(), 3.into()),
        ],
        vec![
            BigRational::new(1.into(), 4.into()),
            BigRational::new(1.into(), 5.into()),
        ],
    ];
    let want = BigRational::new(1.into(), 10.into()) - BigRational::new(1.into(), 12.into());
    assert_eq!(det_rational(&to_sparse(&h)), Some(want.clone()));
    assert_eq!(qlinalg::det(&h), want);
}

#[test]
fn sparse_exact_solvers_agree() {
    for seed in 100..115u64 {
        let n = 3 + (seed % 5) as usize;
        let a = random_int_mat(n, 35, seed, true);
        let sp = to_sparse(&a);
        let mut rng = Lcg(seed + 1);
        let b: Vec<Q> = (0..n).map(|_| q(rng.int(-9, 9))).collect();
        let dixon = solve_dixon(&sp, &b).expect("non-singular");
        assert_eq!(
            sp.matvec(&mut Rationals, &dixon).unwrap(),
            b,
            "Dixon residual, seed {seed}"
        );
        let inv = qlinalg::inverse(&a).unwrap();
        assert_eq!(
            qlinalg::matvec(&inv, &b),
            dixon,
            "vs dense inverse, seed {seed}"
        );
        match solve(&mut Rationals, &sp, &b).unwrap() {
            | Solution::Family {
                particular,
                nullspace,
            } => {
                assert_eq!(particular, dixon);
                assert!(nullspace.is_empty());
            },
            | Solution::Inconsistent => panic!("seed {seed}: nonsingular reported inconsistent"),
        }
    }
    // shape errors
    let sp = to_sparse(&qmat(&[&[1, 2], &[3, 4]]));
    assert!(solve_dixon(&sp, &[q(1)]).is_none());
    assert!(solve(&mut Rationals, &sp, &[q(1)]).is_none());
    assert!(solve_dixon(&SparseMat::zeros(2, 3), &[q(1), q(1)]).is_none());
    assert_eq!(
        solve_dixon(&SparseMat::<BigRational>::zeros(0, 0), &[]),
        Some(vec![])
    );
    // singular: Dixon declines, solve reports a family or inconsistency
    let sing = to_sparse(&qmat(&[&[1, 1], &[1, 1]]));
    assert!(solve_dixon(&sing, &[q(2), q(2)]).is_none());
    assert_eq!(
        solve(&mut Rationals, &sing, &[q(1), q(2)]),
        Some(Solution::Inconsistent)
    );
    match solve(&mut Rationals, &sing, &[q(2), q(2)]).unwrap() {
        | Solution::Family {
            particular,
            nullspace,
        } => {
            assert_eq!(
                sing.matvec(&mut Rationals, &particular).unwrap(),
                vec![q(2), q(2)]
            );
            assert_eq!(nullspace.len(), 1);
            assert_eq!(
                sing.matvec(&mut Rationals, &nullspace[0]).unwrap(),
                vec![q(0), q(0)]
            );
        },
        | Solution::Inconsistent => panic!("consistent system reported inconsistent"),
    }
}

#[test]
fn sparse_exact_nullspace_and_rank_nullity() {
    for seed in 200..220u64 {
        let n = 3 + (seed % 5) as usize;
        let a = random_int_mat(n, 30, seed, false);
        let sp = to_sparse(&a);
        let ns = nullspace(&mut Rationals, &sp);
        let r = rank(&mut Rationals, &sp);
        assert_eq!(ns.len() + r, n, "seed {seed}");
        for v in &ns {
            let zero = sp.matvec(&mut Rationals, v).unwrap();
            assert!(zero.iter().all(Zero::is_zero), "seed {seed}");
        }
    }
}

#[test]
fn sparse_matrix_algebra_matches_dense() {
    let a = random_int_mat(5, 50, 31, false);
    let b = random_int_mat(5, 50, 32, false);
    let (sa, sb) = (to_sparse(&a), to_sparse(&b));
    let mut d = Rationals;
    assert_eq!(
        sa.matmul(&mut d, &sb).unwrap().to_dense(&mut d),
        qlinalg::matmul(&a, &b)
    );
    let sum = sa.add(&mut d, &sb).unwrap().to_dense(&mut d);
    for i in 0..5 {
        for j in 0..5 {
            assert_eq!(sum[i][j], a[i][j].clone() + b[i][j].clone());
        }
    }
    assert_eq!(sa.transpose().to_dense(&mut d), qlinalg::transpose(&a, 5));
    let scaled = sa.scale(&mut d, &q(-3)).to_dense(&mut d);
    for i in 0..5 {
        for j in 0..5 {
            assert_eq!(scaled[i][j], a[i][j].clone() * q(-3));
        }
    }
    // scaling by zero empties the matrix
    assert_eq!(sa.scale(&mut d, &q(0)).nnz(), 0);
    // a + (-a) cancels exactly and stores no zeros
    let neg = sa.scale(&mut d, &q(-1));
    assert_eq!(sa.add(&mut d, &neg).unwrap().nnz(), 0);
    // shapes
    assert!(sa.add(&mut d, &SparseMat::zeros(4, 5)).is_none());
    assert!(sa.matmul(&mut d, &SparseMat::zeros(4, 4)).is_none());
    assert!(sa.matvec(&mut d, &[q(1)]).is_none());
    // triplets: duplicates summed, cancelling pairs vanish, out of range is an error
    let t = SparseMat::from_triplets(
        &mut d,
        2,
        2,
        [(0, 0, q(1)), (0, 0, q(2)), (1, 1, q(5)), (1, 1, q(-5))],
    )
    .unwrap();
    assert_eq!(t.nnz(), 1);
    assert_eq!(
        t.triplets()
            .map(|(i, j, v)| (i, j, v.clone()))
            .collect::<Vec<_>>(),
        vec![(0, 0, q(3))]
    );
    assert!(SparseMat::from_triplets(&mut d, 2, 2, [(2, 0, q(1))]).is_err());
    assert!(SparseMat::from_triplets(&mut d, 2, 2, [(0, 2, q(1))]).is_err());
}

#[test]
fn prime_field_helpers() {
    // Miller-Rabin against a sieve
    let limit = 5000usize;
    let mut sieve = vec![true; limit + 1];
    sieve[0] = false;
    sieve[1] = false;
    for i in 2..=limit {
        if sieve[i] {
            let mut j = i * i;
            while j <= limit {
                sieve[j] = false;
                j += i;
            }
        }
    }
    for n in 0..=limit {
        assert_eq!(is_prime_u64(n as u64), sieve[n], "n = {n}");
    }
    for p in [
        2_305_843_009_213_693_951_u64,
        1_000_000_007,
        18_446_744_073_709_551_557,
    ] {
        assert!(is_prime_u64(p), "{p}");
    }
    // Carmichael numbers, strong pseudoprimes to small bases, squares of primes
    for c in [
        561_u64,
        41041,
        3_215_031_751,
        3_825_123_056_546_413_051,
        1_000_000_007 * 998_244_353,
        4_294_967_291 * 4_294_967_279,
    ] {
        assert!(!is_prime_u64(c), "{c}");
    }
    assert!(!is_prime_u64(u64::MAX));
    // word primes: descending, below 2^62, prime
    let ps: Vec<u64> = word_primes().take(5).collect();
    assert_eq!(ps.len(), 5);
    assert!(ps.windows(2).all(|w| w[0] > w[1]));
    assert!(ps.iter().all(|&p| p < 1 << 62 && is_prime_u64(p)));
    // Gfp arithmetic
    let p = Gfp(2_305_843_009_213_693_951);
    let mut rng = Lcg(5);
    for _ in 0..100 {
        let (a, b) = (rng.next() << 20 | rng.next(), rng.next() << 25 | rng.next());
        assert_eq!(
            p.mul_mod(a, b) as u128,
            (a as u128 * b as u128) % p.0 as u128
        );
        let a = a % p.0;
        if a != 0 {
            assert_eq!(p.mul_mod(a, p.inv(a).unwrap()), 1);
            assert_eq!(p.pow_mod(a, p.0 - 1), 1, "Fermat");
        }
    }
    assert!(p.inv(0).is_none());
    assert!(p.inv(p.0).is_none());
    assert_eq!(p.pow_mod(7, 0), 1);
    assert_eq!(p.reduce(&BigInt::from(-1)), p.0 - 1);
}

#[test]
fn rational_reconstruction_round_trip() {
    let m = BigInt::from(2_305_843_009_213_693_951_u64) * BigInt::from(1_000_000_007_u64);
    let mut rng = Lcg(3);
    for _ in 0..50 {
        let n = rng.int(-1_000_000, 1_000_000);
        let d = rng.int(1, 1_000_000);
        // a = n * d^{-1} (mod m)
        let egcd = BigInt::from(d).extended_gcd_for_test(&m);
        let a = (BigInt::from(n) * egcd) % &m;
        let got = rational_reconstruct(&a, &m).expect("small fraction reconstructs");
        assert_eq!(got, BigRational::new(n.into(), d.into()), "{n}/{d}");
    }
    // a fraction too large for the modulus has no reconstruction or a different one: never panics
    let _ = rational_reconstruct(
        &BigInt::from(123_456_789_012_345_u64),
        &BigInt::from(1_000_003),
    );
    let _ = rational_reconstruct(&BigInt::zero(), &BigInt::one());
}

/// Modular inverse via the extended Euclidean algorithm.
trait ModInv {
    fn extended_gcd_for_test(
        &self,
        m: &BigInt,
    ) -> BigInt;
}
impl ModInv for BigInt {
    fn extended_gcd_for_test(
        &self,
        m: &BigInt,
    ) -> BigInt {
        let (mut r0, mut r1) = (m.clone(), self.clone());
        let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
        while !r1.is_zero() {
            let qt = &r0 / &r1;
            let r2 = &r0 - &qt * &r1;
            let t2 = &t0 - &qt * &t1;
            r0 = r1;
            r1 = r2;
            t0 = t1;
            t1 = t2;
        }
        assert!(r0.is_one(), "not invertible");
        ((t0 % m) + m) % m
    }
}

#[test]
fn black_box_methods_over_a_prime_field() {
    let p = Gfp(1_000_003);
    // Berlekamp-Massey on Fibonacci numbers mod p: order-2 recurrence
    let mut fib = vec![1_u64, 1];
    for i in 2..30 {
        let next = (fib[i - 1] + fib[i - 2]) % p.0;
        fib.push(next);
    }
    let c = berlekamp_massey(&fib, p);
    assert_eq!(c.len(), 3, "recurrence order 2, got {c:?}");
    assert_eq!(*c.last().unwrap(), 1, "monic");
    for i in 0..fib.len() - 2 {
        let s: u64 = (0..3)
            .map(|k| p.mul_mod(c[k], fib[i + k]))
            .fold(0, |a, b| (a + b) % p.0);
        assert_eq!(s, 0, "recurrence fails at {i}: {c:?}");
    }
    // LFSR of order 4 with random coefficients
    let coeffs = [3_u64, 5, 7, 11];
    let mut seq = vec![1_u64, 2, 3, 4];
    for i in 4..60 {
        let mut v = 0u64;
        for (k, &ck) in coeffs.iter().enumerate() {
            v = (v + p.mul_mod(ck, seq[i - 1 - k])) % p.0;
        }
        seq.push(v);
    }
    let c = berlekamp_massey(&seq, p);
    assert_eq!(c.len(), 5);
    for i in 0..seq.len() - 4 {
        let s: u64 = (0..5)
            .map(|k| p.mul_mod(c[k], seq[i + k]))
            .fold(0, |a, b| (a + b) % p.0);
        assert_eq!(s, 0, "order 4, index {i}");
    }
    // empty and constant-zero sequences do not panic
    let _ = berlekamp_massey(&[], p);
    let _ = berlekamp_massey(&[0, 0, 0, 0], p);

    // minimal polynomial of diag(1,2,3) is (x-1)(x-2)(x-3)
    let mut dom = p;
    let diag = SparseMat::from_triplets(&mut dom, 3, 3, [(0, 0, 1), (1, 1, 2), (2, 2, 3)]).unwrap();
    let minpoly = minpoly_wiedemann(&diag, p, 3, 42).unwrap();
    let want = vec![p.0 - 6, 11, p.0 - 6, 1];
    assert_eq!(minpoly, want);
    // repeated eigenvalues lower the degree
    let rep =
        SparseMat::from_triplets(&mut dom, 4, 4, [(0, 0, 2), (1, 1, 2), (2, 2, 2), (3, 3, 5)])
            .unwrap();
    assert_eq!(minpoly_wiedemann(&rep, p, 3, 7).unwrap().len(), 3);
    assert!(minpoly_wiedemann(&SparseMat::<u64>::zeros(2, 3), p, 1, 1).is_none());
    // determinism
    assert_eq!(
        minpoly_wiedemann(&diag, p, 2, 9),
        minpoly_wiedemann(&diag, p, 2, 9)
    );
}

#[test]
fn wiedemann_solves_sparse_systems() {
    let p = Gfp(1_000_000_007);
    let n = 12;
    let a = random_int_mat(n, 25, 4242, true);
    let rows: Vec<Vec<u64>> = a
        .iter()
        .map(|r| r.iter().map(|v| p.reduce(&v.to_integer())).collect())
        .collect();
    let mut dom = p;
    let sp = SparseMat::from_dense(&mut dom, &rows);
    let mut rng = Lcg(8);
    let b: Vec<u64> = (0..n).map(|_| rng.next() % p.0).collect();
    let x = solve_wiedemann(&sp, &b, p, 5, 1).expect("solvable");
    assert_eq!(sp.matvec(&mut dom, &x).unwrap(), b);
    // wrong shapes
    assert!(solve_wiedemann(&sp, &b[..3], p, 1, 1).is_none());
    assert!(solve_wiedemann(&SparseMat::<u64>::zeros(2, 3), &[1, 1], p, 1, 1).is_none());
    // singular: either None or a verified solution
    let sing = SparseMat::from_dense(&mut dom, &[vec![1, 1], vec![1, 1]]);
    if let Some(x) = solve_wiedemann(&sing, &[1, 1], p, 3, 1) {
        assert_eq!(sing.matvec(&mut dom, &x).unwrap(), vec![1, 1]);
    }
}
