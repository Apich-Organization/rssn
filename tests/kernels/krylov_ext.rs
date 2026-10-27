//! Tests for MINRES, BiCGSTAB(l), IDR(s), FGMRES and the Lanczos/Arnoldi eigensolvers.

use rssn::kernels::eigs::{EigsOptions, Which, arnoldi, eigs, eigs_shift_invert, eigs_sym, eigs_sym_shift_invert};
use rssn::kernels::krylov::{Csr, SparseOrdering, gmres, identity_precond};
use rssn::kernels::krylov_ext::{IdrOptions, bicgstab_l, fgmres, idr_s_with, minres};

fn laplace2d(nx: usize, shift: f64) -> Csr {
    let n = nx * nx;
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let v = j * nx + i;
            t.push((v, v, 4.0 - shift));
            if i + 1 < nx {
                t.push((v, v + 1, -1.0));
                t.push((v + 1, v, -1.0));
            }
            if j + 1 < nx {
                t.push((v, v + nx, -1.0));
                t.push((v + nx, v, -1.0));
            }
        }
    }
    Csr::from_triplets(n, n, &t)
}

fn conv_diff(nx: usize, peclet: f64) -> Csr {
    let n = nx * nx;
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let v = j * nx + i;
            t.push((v, v, 4.0));
            if i + 1 < nx {
                t.push((v, v + 1, -1.0 + peclet));
                t.push((v + 1, v, -1.0 - peclet));
            }
            if j + 1 < nx {
                t.push((v, v + nx, -1.0 + 0.5 * peclet));
                t.push((v + nx, v, -1.0 - 0.5 * peclet));
            }
        }
    }
    Csr::from_triplets(n, n, &t)
}

fn rhs(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1.0 + (i as f64 * 0.31).sin()).collect()
}

fn true_res(a: &Csr, x: &[f64], b: &[f64]) -> f64 {
    let mut ax = vec![0.0; b.len()];
    a.matvec(x, &mut ax);
    let r: f64 = ax.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum::<f64>().sqrt();
    r / b.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn minres_spd_and_indefinite() {
    let b = rhs(400);
    let a = laplace2d(20, 0.0);
    let r = minres(|x, y| a.matvec(x, y), identity_precond, &b, None, 1e-10, 2000);
    assert!(r.converged && true_res(&a, &r.x, &b) < 1e-9, "{r:?}");
    // indefinite: shift inside the spectrum
    let a = laplace2d(20, 2.013);
    let r = minres(|x, y| a.matvec(x, y), identity_precond, &b, None, 1e-10, 4000);
    assert!(true_res(&a, &r.x, &b) < 1e-8, "{}", true_res(&a, &r.x, &b));
    // Jacobi preconditioner (SPD) on the SPD matrix
    let a = laplace2d(20, 0.0);
    let dinv = a.jacobi_inverse();
    let r2 = minres(
        |x, y| a.matvec(x, y),
        |r, z| {
            for i in 0..r.len() {
                z[i] = dinv[i] * r[i];
            }
        },
        &b,
        None,
        1e-10,
        2000,
    );
    assert!(true_res(&a, &r2.x, &b) < 1e-8);
}

#[test]
fn bicgstab_l_variants() {
    let a = conv_diff(24, 0.4);
    let b = rhs(a.rows);
    let mut iters = Vec::new();
    for l in [1usize, 2, 4] {
        let r = bicgstab_l(|x, y| a.matvec(x, y), identity_precond, l, &b, None, 1e-10, 2000);
        assert!(true_res(&a, &r.x, &b) < 1e-8, "l={l}: {}", true_res(&a, &r.x, &b));
        iters.push(r.iterations);
    }
    // with ILU(0) preconditioning
    let ilu = a.ilu0().unwrap();
    let r = bicgstab_l(|x, y| a.matvec(x, y), |r, z| ilu.apply(r, z), 2, &b, None, 1e-12, 500);
    assert!(true_res(&a, &r.x, &b) < 1e-10 && r.iterations < iters[1]);
}

#[test]
fn idr_s_variants() {
    let a = conv_diff(24, 0.4);
    let b = rhs(a.rows);
    let g = gmres(|x, y| a.matvec(x, y), identity_precond, &b, None, 200, 1e-10, 2000);
    assert!(g.converged);
    for s in [1usize, 2, 4, 8] {
        let r = idr_s_with(|x, y| a.matvec(x, y), identity_precond, &b, None, &IdrOptions { s, tol: 1e-10, max_iter: 3000, seed: 7 });
        assert!(true_res(&a, &r.x, &b) < 1e-8, "s={s}: {}", true_res(&a, &r.x, &b));
        if s >= 4 {
            assert!(r.iterations < 3 * g.iterations + 50, "s={s}: {} vs gmres {}", r.iterations, g.iterations);
        }
    }
    let ilu = a.ilu0().unwrap();
    let r = idr_s_with(|x, y| a.matvec(x, y), |r, z| ilu.apply(r, z), &b, None, &IdrOptions { s: 4, tol: 1e-12, max_iter: 500, seed: 3 });
    assert!(true_res(&a, &r.x, &b) < 1e-10);
}

#[test]
fn fgmres_variable_preconditioner() {
    let a = conv_diff(24, 0.4);
    let b = rhs(a.rows);
    let ilu = a.ilu0().unwrap();
    let calls = std::cell::Cell::new(0usize);
    // preconditioner that changes every call: ILU followed by a varying number of Jacobi sweeps
    let dinv = a.jacobi_inverse();
    let pre = |r: &[f64], z: &mut [f64]| {
        calls.set(calls.get() + 1);
        ilu.apply(r, z);
        let sweeps = calls.get() % 3;
        let mut az = vec![0.0; r.len()];
        for _ in 0..sweeps {
            a.matvec(z, &mut az);
            for i in 0..r.len() {
                z[i] += dinv[i] * (r[i] - az[i]);
            }
        }
    };
    let r = fgmres(|x, y| a.matvec(x, y), pre, &b, None, 30, 1e-11, 300);
    assert!(r.converged && true_res(&a, &r.x, &b) < 1e-9, "{r:?}");
    // identity preconditioner reproduces plain GMRES behaviour
    let r0 = fgmres(|x, y| a.matvec(x, y), identity_precond, &b, None, 60, 1e-10, 1000);
    assert!(true_res(&a, &r0.x, &b) < 1e-8);
}

#[test]
fn lanczos_symmetric_laplacian() {
    let n = 400;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i + 1 < n {
            t.push((i, i + 1, -1.0));
            t.push((i + 1, i, -1.0));
        }
    }
    let a = Csr::from_triplets(n, n, &t);
    let exact = |k: usize| 2.0 - 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
    let opts = EigsOptions { tol: 1e-10, ..EigsOptions::default() };
    let r = eigs_sym(|x, y| a.matvec(x, y), n, 4, Which::LargestReal, &opts).unwrap();
    assert!(r.converged, "{r:?}");
    for q in 0..4 {
        assert!((r.values[q].0 - exact(n - q)).abs() < 1e-8, "{} vs {}", r.values[q].0, exact(n - q));
        // residual check with the returned vector
        let mut ax = vec![0.0; n];
        a.matvec(&r.vectors_re[q], &mut ax);
        let res: f64 = ax
            .iter()
            .zip(&r.vectors_re[q])
            .map(|(p, v)| (p - r.values[q].0 * v).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(res < 1e-7, "residual {res}");
    }
    // shift-invert around sigma = 0.5
    let sigma = 0.5;
    let mut shifted = t.clone();
    for i in 0..n {
        shifted.push((i, i, -sigma));
    }
    let lu = Csr::from_triplets(n, n, &shifted).sparse_lu(SparseOrdering::Rcm).unwrap();
    let r = eigs_sym_shift_invert(|x, y| y.copy_from_slice(&lu.solve(x)), n, 3, sigma, &opts).unwrap();
    let mut all: Vec<f64> = (1..=n).map(exact).collect();
    all.sort_by(|p, q| (p - sigma).abs().total_cmp(&(q - sigma).abs()));
    for q in 0..3 {
        assert!((r.values[q].0 - all[q]).abs() < 1e-8, "{} vs {}", r.values[q].0, all[q]);
    }
    // smallest algebraic without inversion
    let r = eigs_sym(|x, y| a.matvec(x, y), n, 2, Which::SmallestReal, &EigsOptions { tol: 1e-8, max_restarts: 2000, ..opts }).unwrap();
    assert!((r.values[0].0 - exact(1)).abs() < 1e-7 && (r.values[1].0 - exact(2)).abs() < 1e-7);
    assert!(eigs_sym(|x, y| a.matvec(x, y), n, 0, Which::LargestReal, &opts).is_err());
}

#[test]
fn arnoldi_nonsymmetric_real_spectrum() {
    // tridiagonal Toeplitz (sub -1, diag 2, super -0.5): lambda_k = 2 + 2 sqrt(0.3) cos(k pi / (n + 1))
    let n = 100;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i + 1 < n {
            t.push((i, i + 1, -0.5));
            t.push((i + 1, i, -0.6));
        }
    }
    let a = Csr::from_triplets(n, n, &t);
    let exact = |k: usize| 2.0 + 2.0 * 0.3_f64.sqrt() * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
    let opts = EigsOptions { tol: 1e-9, ncv: 40, max_restarts: 3000, ..EigsOptions::default() };
    let r = eigs(|x, y| a.matvec(x, y), n, 3, Which::LargestReal, &opts).unwrap();
    assert!(r.converged, "restarts {}", r.restarts);
    for q in 0..3 {
        assert!((r.values[q].0 - exact(1 + q)).abs() < 1e-6, "{:?} vs {}", r.values[q], exact(1 + q));
        assert!(r.values[q].1.abs() < 1e-8);
    }
    // shift-invert near 2.0
    let sigma = 2.05;
    let mut sh = t.clone();
    for i in 0..n {
        sh.push((i, i, -sigma));
    }
    let lu = Csr::from_triplets(n, n, &sh).sparse_lu(SparseOrdering::Rcm).unwrap();
    let r = eigs_shift_invert(|x, y| y.copy_from_slice(&lu.solve(x)), n, 2, sigma, &opts).unwrap();
    let mut all: Vec<f64> = (1..=n).map(exact).collect();
    all.sort_by(|p, q| (p - sigma).abs().total_cmp(&(q - sigma).abs()));
    for q in 0..2 {
        assert!((r.values[q].0 - all[q]).abs() < 1e-6, "{:?} vs {}", r.values[q], all[q]);
    }
}

#[test]
fn arnoldi_complex_pairs() {
    // block upper-triangular matrix: 2x2 rotation-scaling blocks on the diagonal,
    // arbitrary coupling above it; eigenvalues alpha_i +- i beta_i
    let nb = 40;
    let n = 2 * nb;
    let mut t = Vec::new();
    let mut want: Vec<(f64, f64)> = Vec::new();
    for b in 0..nb {
        let alpha = 0.05 * b as f64;
        let beta = 1.0 + 0.15 * b as f64;
        t.push((2 * b, 2 * b, alpha));
        t.push((2 * b + 1, 2 * b + 1, alpha));
        t.push((2 * b, 2 * b + 1, -beta));
        t.push((2 * b + 1, 2 * b, beta));
        want.push((alpha, beta));
        want.push((alpha, -beta));
        if b + 1 < nb {
            t.push((2 * b, 2 * b + 2, 0.3));
            t.push((2 * b + 1, 2 * b + 3, -0.2));
        }
    }
    let a = Csr::from_triplets(n, n, &t);
    want.sort_by(|p, q| q.0.hypot(q.1).total_cmp(&p.0.hypot(p.1)).then(q.1.total_cmp(&p.1)));
    let opts = EigsOptions { tol: 1e-9, ncv: 40, max_restarts: 2000, ..EigsOptions::default() };
    let r = eigs(|x, y| a.matvec(x, y), n, 4, Which::LargestMagnitude, &opts).unwrap();
    assert!(r.converged);
    assert_eq!(r.values.len(), 4);
    let mut got = r.values.clone();
    got.sort_by(|p, q| q.0.hypot(q.1).total_cmp(&p.0.hypot(p.1)).then(q.1.total_cmp(&p.1)));
    for q in 0..4 {
        assert!((got[q].0 - want[q].0).abs() < 1e-6 && (got[q].1 - want[q].1).abs() < 1e-6, "{:?} vs {:?}", got[q], want[q]);
    }
    // complex eigenvector residual: A (xr + i xi) = (lr + i li)(xr + i xi)
    for q in 0..4 {
        let (lr, li) = r.values[q];
        let (xr, xi) = (&r.vectors_re[q], &r.vectors_im[q]);
        let (mut axr, mut axi) = (vec![0.0; n], vec![0.0; n]);
        a.matvec(xr, &mut axr);
        a.matvec(xi, &mut axi);
        let mut res = 0.0;
        for i in 0..n {
            res += (axr[i] - (lr * xr[i] - li * xi[i])).powi(2) + (axi[i] - (lr * xi[i] + li * xr[i])).powi(2);
        }
        assert!(res.sqrt() < 1e-6, "residual {}", res.sqrt());
    }
}

#[test]
fn arnoldi_factorisation_relation() {
    let a = conv_diff(8, 0.3);
    let v0 = rhs(a.rows);
    let m = 12;
    let (v, h, beta) = arnoldi(|x, y| a.matvec(x, y), &v0, m);
    assert_eq!(v.len(), m + 1);
    for j in 0..m {
        let mut av = vec![0.0; a.rows];
        a.matvec(&v[j], &mut av);
        for e in 0..a.rows {
            let mut rhs = 0.0;
            for i in 0..m {
                rhs += v[i][e] * h.at(i, j);
            }
            if j == m - 1 {
                rhs += beta * v[m][e];
            }
            assert!((av[e] - rhs).abs() < 1e-10, "col {j}");
        }
    }
}
