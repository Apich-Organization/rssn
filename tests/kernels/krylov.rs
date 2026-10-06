//! Tests for the sparse matrix and Krylov solver kernels.

use rssn::kernels::krylov::{Csr, bicgstab, cg, gmres, identity_precond};

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
