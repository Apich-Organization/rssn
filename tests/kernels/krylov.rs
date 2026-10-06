//! Tests for the sparse matrix and Krylov solver kernels.

use rssn::kernels::krylov::{
    Csr, SparseLuError, SparseOrdering, bicgstab, cg, gmres, identity_precond,
};

fn laplacian(n: usize) -> Csr {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i > 0 {
            t.push((i, i - 1, -1.0));
        }
        if i + 1 < n {
            t.push((i, i + 1, -1.0));
        }
    }
    Csr::from_triplets(n, n, &t)
}

fn convection(n: usize) -> Csr {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 3.0));
        if i > 0 {
            t.push((i, i - 1, -2.0));
        }
        if i + 1 < n {
            t.push((i, i + 1, -0.5));
        }
    }
    Csr::from_triplets(n, n, &t)
}

fn check(a: &Csr, x: &[f64], b: &[f64], tol: f64) {
    let mut ax = vec![0.0; b.len()];
    a.matvec(x, &mut ax);
    let r: f64 = ax.iter().zip(b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt();
    assert!(r < tol, "residual {r}");
}

#[test]
fn csr_basics() {
    let a = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 0, 2.0), (1, 1, 4.0), (0, 1, 5.0)]);
    assert_eq!(a.diagonal(), vec![3.0, 4.0]);
    let mut y = [0.0; 2];
    a.matvec(&[1.0, 1.0], &mut y);
    assert_eq!(y, [8.0, 4.0]);
}

#[test]
fn cg_with_preconditioners() {
    let n = 100;
    let a = laplacian(n);
    let b: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let plain = cg(|x, y| a.matvec(x, y), identity_precond, &b, None, 1e-10, 1000);
    assert!(plain.converged);
    check(&a, &plain.x, &b, 1e-8);
    let dinv = a.jacobi_inverse();
    let jac = cg(|x, y| a.matvec(x, y), |r, z| { for i in 0..r.len() { z[i] = dinv[i] * r[i]; } }, &b, None, 1e-10, 1000);
    assert!(jac.converged);
    let ilu = a.ilu0().unwrap();
    let pre = cg(|x, y| a.matvec(x, y), |r, z| ilu.apply(r, z), &b, None, 1e-10, 1000);
    assert!(pre.converged);
    // tridiagonal ILU(0) is the exact LU: one iteration
    assert!(pre.iterations <= 2, "{}", pre.iterations);
    check(&a, &pre.x, &b, 1e-8);
}

#[test]
fn gmres_and_bicgstab_nonsymmetric() {
    let n = 200;
    let a = convection(n);
    let b = vec![1.0; n];
    let g = gmres(|x, y| a.matvec(x, y), identity_precond, &b, None, 30, 1e-10, 2000);
    assert!(g.converged, "{g:?}");
    check(&a, &g.x, &b, 1e-8);
    let ilu = a.ilu0().unwrap();
    let gp = gmres(|x, y| a.matvec(x, y), |r, z| ilu.apply(r, z), &b, None, 10, 1e-12, 100);
    assert!(gp.converged && gp.iterations <= 3, "{gp:?}");
    let bi = bicgstab(|x, y| a.matvec(x, y), identity_precond, &b, None, 1e-10, 1000);
    assert!(bi.converged, "{bi:?}");
    check(&a, &bi.x, &b, 1e-8);
    let bp = bicgstab(|x, y| a.matvec(x, y), |r, z| ilu.apply(r, z), &b, None, 1e-12, 100);
    assert!(bp.converged);
}

#[test]
fn ilu_requires_diagonal() {
    let a = Csr::from_triplets(2, 2, &[(0, 1, 1.0), (1, 0, 1.0)]);
    assert!(a.ilu0().is_none());
}

fn poisson2d(m: usize) -> Csr {
    let n = m * m;
    let mut t = Vec::new();
    for i in 0..m {
        for j in 0..m {
            let r = i * m + j;
            t.push((r, r, 4.0));
            if i > 0 {
                t.push((r, r - m, -1.0));
            }
            if i + 1 < m {
                t.push((r, r + m, -1.0));
            }
            if j > 0 {
                t.push((r, r - 1, -1.0));
            }
            if j + 1 < m {
                t.push((r, r + 1, -1.0));
            }
        }
    }
    Csr::from_triplets(n, n, &t)
}

fn rhs(n: usize) -> Vec<f64> {
    (0..n).map(|i| ((i * 7919 % 101) as f64) / 50.0 - 1.0).collect()
}

#[test]
fn sparse_lu_poisson_matches_cg() {
    let a = poisson2d(24);
    let n = a.rows;
    let b = rhs(n);
    let cgres = cg(|x, y| a.matvec(x, y), identity_precond, &b, None, 1e-13, 5000);
    assert!(cgres.converged);
    for ord in [SparseOrdering::Natural, SparseOrdering::Rcm, SparseOrdering::MinDegree] {
        let lu = a.sparse_lu(ord).unwrap();
        let x = lu.solve(&b);
        check(&a, &x, &b, 1e-10);
        let d: f64 = x.iter().zip(&cgres.x).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        assert!(d < 1e-8, "{ord:?}: {d}");
        assert_eq!(lu.dim(), n);
    }
}

#[test]
fn sparse_lu_ordering_reduces_fill() {
    let a = poisson2d(30);
    let n = a.rows;
    // scramble the natural ordering with a fixed pseudo-random permutation
    let mut perm: Vec<usize> = (0..n).collect();
    let mut s = 12345u64;
    for i in (1..n).rev() {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let j = (s >> 33) as usize % (i + 1);
        perm.swap(i, j);
    }
    let scrambled = a.permute_symmetric(&perm);
    assert!(scrambled.bandwidth() > 10 * a.bandwidth());
    let nat = scrambled.sparse_lu(SparseOrdering::Natural).unwrap();
    let rcm = scrambled.sparse_lu(SparseOrdering::Rcm).unwrap();
    let md = scrambled.sparse_lu(SparseOrdering::MinDegree).unwrap();
    assert!(rcm.nnz() * 3 < nat.nnz(), "rcm {} natural {}", rcm.nnz(), nat.nnz());
    assert!(md.nnz() * 4 < nat.nnz(), "md {} natural {}", md.nnz(), nat.nnz());
    // minimum degree beats RCM on 2D grids
    assert!(md.nnz() < rcm.nnz());
    // RCM restores a banded matrix
    let p = scrambled.rcm_ordering();
    let banded = scrambled.permute_symmetric(&p);
    assert!(banded.bandwidth() <= 2 * a.bandwidth(), "bw {}", banded.bandwidth());
    // all orderings agree on the solution
    let b = rhs(n);
    for lu in [&nat, &rcm, &md] {
        check(&scrambled, &lu.solve(&b), &b, 1e-9);
    }
}

#[test]
fn sparse_lu_nonsymmetric_and_pivoting() {
    let a = convection(200);
    let b = rhs(200);
    let x = a.sparse_lu(SparseOrdering::Rcm).unwrap().solve(&b);
    check(&a, &x, &b, 1e-10);
    // zero diagonal forces row pivoting: [[0,1,0],[2,0,1],[0,3,1]]
    let a = Csr::from_triplets(3, 3, &[(0, 1, 1.0), (1, 0, 2.0), (1, 2, 1.0), (2, 1, 3.0), (2, 2, 1.0)]);
    let b = [1.0, 2.0, 3.0];
    for ord in [SparseOrdering::Natural, SparseOrdering::Rcm, SparseOrdering::MinDegree] {
        let x = a.sparse_lu_with(ord, 1.0).unwrap().solve(&b);
        check(&a, &x, &b, 1e-12);
    }
    // random sparse unsymmetric diagonally dominant matrix
    let n = 300;
    let mut t = Vec::new();
    let mut s = 99u64;
    for i in 0..n {
        t.push((i, i, 10.0));
        for _ in 0..3 {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let j = (s >> 33) as usize % n;
            t.push((i, j, ((s >> 20) % 100) as f64 / 100.0 - 0.5));
        }
    }
    let a = Csr::from_triplets(n, n, &t);
    let b = rhs(n);
    for ord in [SparseOrdering::Natural, SparseOrdering::Rcm, SparseOrdering::MinDegree] {
        check(&a, &a.sparse_lu(ord).unwrap().solve(&b), &b, 1e-10);
    }
}

#[test]
fn sparse_lu_errors() {
    let r = Csr::from_triplets(2, 3, &[(0, 0, 1.0)]);
    assert_eq!(r.sparse_lu(SparseOrdering::Natural).unwrap_err(), SparseLuError::NotSquare);
    let sing = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 4.0)]);
    assert!(matches!(sing.sparse_lu(SparseOrdering::Natural), Err(SparseLuError::Singular(_))));
    // empty pattern column
    let sing = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (1, 0, 1.0)]);
    assert!(sing.sparse_lu(SparseOrdering::Rcm).is_err());
}
