//! Tests for Krylov-Schur, block Lanczos, Jacobi-Davidson and generalised eigensolvers.

use num_complex::Complex64;
use rssn::kernels::dense::{self, Mat};
use rssn::kernels::eigs::{EigsOptions, Which};
use rssn::kernels::eigs_ext::{
    GenEigsOptions, JdOptions, complex_schur, eigs_gen_shift_invert, eigs_gen_sym_shift_invert,
    eigs_jd_sym, eigs_krylov_schur, eigs_sym_block, reorder_schur,
};
use rssn::kernels::krylov::{Csr, SparseOrdering};

type C = Complex64;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 33) as f64) / f64::from(1u32 << 31) - 0.5
}

fn cmatmul(a: &[Vec<C>], b: &[Vec<C>]) -> Vec<Vec<C>> {
    let n = a.len();
    let m = b[0].len();
    (0..n)
        .map(|i| (0..m).map(|j| (0..b.len()).map(|k| a[i][k] * b[k][j]).sum()).collect())
        .collect()
}

fn adjoint(a: &[Vec<C>]) -> Vec<Vec<C>> {
    let n = a.len();
    (0..n).map(|i| (0..n).map(|j| a[j][i].conj()).collect()).collect()
}

fn max_diff(a: &[Vec<C>], b: &[Vec<C>]) -> f64 {
    let mut m = 0.0f64;
    for i in 0..a.len() {
        for j in 0..a[i].len() {
            m = m.max((a[i][j] - b[i][j]).norm());
        }
    }
    m
}

#[test]
fn complex_schur_reconstructs_and_reorders() {
    let n = 14;
    let mut seed = 42u64;
    let a: Vec<Vec<C>> =
        (0..n).map(|_| (0..n).map(|_| C::new(lcg(&mut seed), lcg(&mut seed))).collect()).collect();
    let (mut t, mut q) = complex_schur(&a).unwrap();
    let recon = cmatmul(&cmatmul(&q, &t), &adjoint(&q));
    assert!(max_diff(&recon, &a) < 1e-12);
    for i in 0..n {
        for j in 0..i {
            assert!(t[i][j].norm() < 1e-13);
        }
    }
    let eye: Vec<Vec<C>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { C::new(1.0, 0.0) } else { C::new(0.0, 0.0) }).collect())
        .collect();
    assert!(max_diff(&cmatmul(&adjoint(&q), &q), &eye) < 1e-13);
    // sort by decreasing real part
    reorder_schur(&mut t, &mut q, n, &|z: C| -z.re);
    for i in 0..n - 1 {
        assert!(t[i][i].re >= t[i + 1][i + 1].re - 1e-12);
    }
    let recon = cmatmul(&cmatmul(&q, &t), &adjoint(&q));
    assert!(max_diff(&recon, &a) < 1e-11);
}

#[test]
fn krylov_schur_non_normal_convection_diffusion() {
    // tridiag(-1-c, 2, -1+c): eigenvalues 2 + 2 sqrt(1 - c^2) cos(j pi/(n+1)), non-normal (eigenvalue condition ~1e8)
    let n = 60;
    let c = 0.3;
    let mut trip = Vec::new();
    for i in 0..n {
        trip.push((i, i, 2.0));
        if i > 0 {
            trip.push((i, i - 1, -1.0 - c));
        }
        if i + 1 < n {
            trip.push((i, i + 1, -1.0 + c));
        }
    }
    let a = Csr::from_triplets(n, n, &trip);
    let opts = EigsOptions { tol: 1e-10, max_restarts: 3000, ..EigsOptions::default() };
    let r = eigs_krylov_schur(|x, y| a.matvec(x, y), n, 4, Which::LargestReal, &opts).unwrap();
    assert!(r.converged, "restarts {}", r.restarts);
    let s = 2.0 * (1.0 - c * c).sqrt();
    for (j, v) in r.values.iter().enumerate() {
        let want = 2.0 + s * (std::f64::consts::PI * (j + 1) as f64 / (n as f64 + 1.0)).cos();
        // the eigenvalue condition number is huge, so allow for it
        assert!((v.0 - want).abs() < 1e-5, "{j}: {} vs {want}", v.0);
        assert!(v.1.abs() < 1e-5);
    }
    // true residual of the returned vectors
    for i in 0..4 {
        let mut ax = vec![0.0; n];
        a.matvec(&r.vectors_re[i], &mut ax);
        let res: f64 = ax
            .iter()
            .zip(&r.vectors_re[i])
            .map(|(p, q)| (p - r.values[i].0 * q).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(res < 1e-6, "residual {res}");
    }
}

#[test]
fn krylov_schur_complex_pairs() {
    // block upper triangular with 2x2 rotation blocks: eigenvalues a_j +- i b_j
    let nb = 25;
    let n = 2 * nb;
    let mut m = Mat::zeros(n, n);
    let mut seed = 7u64;
    let mut want = Vec::new();
    for j in 0..nb {
        let (a, b) = (0.1 * j as f64, 1.0 + 0.05 * j as f64);
        m.set(2 * j, 2 * j, a);
        m.set(2 * j, 2 * j + 1, -b);
        m.set(2 * j + 1, 2 * j, b);
        m.set(2 * j + 1, 2 * j + 1, a);
        want.push((a, b));
    }
    for i in 0..n {
        for j in 2 * (i / 2) + 2..n {
            m.set(i, j, 0.4 * lcg(&mut seed));
        }
    }
    let opts = EigsOptions { tol: 1e-10, max_restarts: 2000, ..EigsOptions::default() };
    let r = eigs_krylov_schur(|x, y| y.copy_from_slice(&m.matvec(x)), n, 6, Which::LargestReal, &opts).unwrap();
    assert!(r.converged);
    // the three largest-real-part pairs: j = 24, 23, 22
    for p in 0..3 {
        let (a, b) = want[nb - 1 - p];
        let found: Vec<_> = r.values.iter().filter(|v| (v.0 - a).abs() < 1e-7).collect();
        assert_eq!(found.len(), 2, "pair {p}: {:?}", r.values);
        assert!(found.iter().all(|v| (v.1.abs() - b).abs() < 1e-7));
        assert!(found[0].1 * found[1].1 < 0.0);
    }
}

fn laplacian2d(nx: usize) -> Csr {
    let n = nx * nx;
    let mut t = Vec::new();
    for i in 0..nx {
        for j in 0..nx {
            let p = i * nx + j;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - nx, -1.0));
            }
            if i + 1 < nx {
                t.push((p, p + nx, -1.0));
            }
            if j > 0 {
                t.push((p, p - 1, -1.0));
            }
            if j + 1 < nx {
                t.push((p, p + 1, -1.0));
            }
        }
    }
    Csr::from_triplets(n, n, &t)
}

#[test]
fn block_lanczos_finds_multiple_eigenvalues() {
    let nx = 30;
    let a = laplacian2d(nx);
    let n = nx * nx;
    let opts = EigsOptions { tol: 1e-9, max_restarts: 2000, ..EigsOptions::default() };
    let k = 8;
    let r = eigs_sym_block(|x, y| a.matvec(x, y), n, k, Which::SmallestReal, 2, &opts).unwrap();
    assert!(r.converged, "restarts {}", r.restarts);
    // exact: 4 - 2cos(p pi/(nx+1)) - 2cos(q pi/(nx+1))
    let th = std::f64::consts::PI / (nx as f64 + 1.0);
    let mut exact: Vec<f64> = Vec::new();
    for p in 1..=nx {
        for q in 1..=nx {
            exact.push(4.0 - 2.0 * (p as f64 * th).cos() - 2.0 * (q as f64 * th).cos());
        }
    }
    exact.sort_by(f64::total_cmp);
    for i in 0..k {
        assert!((r.values[i].0 - exact[i]).abs() < 1e-8, "{i}: {} vs {}", r.values[i].0, exact[i]);
    }
    // both members of every double eigenvalue are present
    assert!((r.values[1].0 - r.values[2].0).abs() < 1e-8);
    assert!((r.values[4].0 - r.values[5].0).abs() < 1e-8);
    assert!((r.values[6].0 - r.values[7].0).abs() < 1e-8);
}

#[test]
fn block_lanczos_tight_cluster() {
    // exact multiplicity 3 and a cluster with gaps 1e-7 at the top of a diagonal operator
    let n = 600;
    let mut d: Vec<f64> = (0..n).map(|i| 10.0 * i as f64 / n as f64).collect();
    d[10] = 100.0;
    d[20] = 100.0;
    d[30] = 100.0;
    d[40] = 90.0;
    d[50] = 90.0 + 1e-7;
    d[60] = 90.0 + 2e-7;
    let opts = EigsOptions { tol: 1e-10, max_restarts: 2000, ..EigsOptions::default() };
    let r = eigs_sym_block(|x, y| {
        for i in 0..n {
            y[i] = d[i] * x[i];
        }
    }, n, 6, Which::LargestReal, 3, &opts)
    .unwrap();
    assert!(r.converged);
    let want = [100.0, 100.0, 100.0, 90.0 + 2e-7, 90.0 + 1e-7, 90.0];
    for i in 0..6 {
        assert!((r.values[i].0 - want[i]).abs() < 1e-9, "{i}: {}", r.values[i].0);
    }
    // invalid arguments
    assert!(eigs_sym_block(|_, _| {}, n, 0, Which::LargestReal, 3, &opts).is_err());
    assert!(eigs_sym_block(|_, _| {}, n, 2, Which::LargestReal, 0, &opts).is_err());
}

fn laplacian1d(n: usize) -> Csr {
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

#[test]
fn jacobi_davidson_interior_eigenvalues() {
    let n = 400;
    let a = laplacian1d(n);
    let target = 2.0;
    let mut exact: Vec<f64> = (1..=n)
        .map(|j| 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos())
        .collect();
    exact.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    // preconditioner: sparse LU of A - (target + 0.003) I
    let mut shifted = a.clone();
    for i in 0..n {
        for p in shifted.indptr[i]..shifted.indptr[i + 1] {
            if shifted.indices[p] == i {
                shifted.values[p] -= target + 0.003;
            }
        }
    }
    let lu = shifted.sparse_lu(SparseOrdering::Rcm).unwrap();
    let pc = |x: &[f64], y: &mut [f64]| y.copy_from_slice(&lu.solve(x));
    let opts = JdOptions { tol: 1e-9, ..JdOptions::default() };
    let r = eigs_jd_sym(|x, y| a.matvec(x, y), n, 3, target, &opts, Some(&pc)).unwrap();
    assert!(r.converged, "iterations {}", r.restarts);
    let mut got: Vec<f64> = r.values.iter().map(|v| v.0).collect();
    got.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    let mut want: Vec<f64> = exact[..3].to_vec();
    want.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    for i in 0..3 {
        assert!((got[i] - want[i]).abs() < 1e-8, "{i}: {} vs {}", got[i], want[i]);
    }
    for i in 0..3 {
        assert!(r.residuals[i] < 1e-8);
    }
}

#[test]
fn jacobi_davidson_without_preconditioner() {
    let n = 100;
    let a = laplacian1d(n);
    let target = 1.0;
    let mut exact: Vec<f64> = (1..=n)
        .map(|j| 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos())
        .collect();
    exact.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    let opts = JdOptions { tol: 1e-8, inner_iters: 60, max_iter: 1500, ..JdOptions::default() };
    let r = eigs_jd_sym(|x, y| a.matvec(x, y), n, 2, target, &opts, None).unwrap();
    assert!(r.converged);
    let mut got: Vec<f64> = r.values.iter().map(|v| v.0).collect();
    got.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    let mut want = exact[..2].to_vec();
    want.sort_by(|x, y| (x - target).abs().total_cmp(&(y - target).abs()));
    for i in 0..2 {
        assert!((got[i] - want[i]).abs() < 1e-7, "{i}: {} vs {}", got[i], want[i]);
    }
    assert!(eigs_jd_sym(|_, _| {}, n, 0, target, &opts, None).is_err());
}

fn fe_matrices(n: usize) -> (Csr, Csr, f64) {
    let h = 1.0 / (n as f64 + 1.0);
    let (mut k, mut m) = (Vec::new(), Vec::new());
    for i in 0..n {
        k.push((i, i, 2.0 / h));
        m.push((i, i, 4.0 * h / 6.0));
        if i > 0 {
            k.push((i, i - 1, -1.0 / h));
            m.push((i, i - 1, h / 6.0));
        }
        if i + 1 < n {
            k.push((i, i + 1, -1.0 / h));
            m.push((i, i + 1, h / 6.0));
        }
    }
    (Csr::from_triplets(n, n, &k), Csr::from_triplets(n, n, &m), h)
}

#[test]
fn generalised_symmetric_shift_invert() {
    let n = 300;
    let (k, m, h) = fe_matrices(n);
    let exact = |j: usize| {
        let c = (j as f64 * std::f64::consts::PI * h).cos();
        6.0 * (1.0 - c) / (h * h * (2.0 + c))
    };
    let sigma = 3000.0;
    let mut sh_trip = Vec::new();
    for i in 0..n {
        for (mat, f) in [(&k, 1.0), (&m, -sigma)] {
            for p in mat.indptr[i]..mat.indptr[i + 1] {
                sh_trip.push((i, mat.indices[p], f * mat.values[p]));
            }
        }
    }
    let shifted = Csr::from_triplets(n, n, &sh_trip);
    let lu = shifted.sparse_lu(SparseOrdering::Amd).unwrap();
    let opts = GenEigsOptions {
        sigma,
        block: 1,
        eigs: EigsOptions { tol: 1e-10, ..EigsOptions::default() },
    };
    let r = eigs_gen_sym_shift_invert(
        |x, y| k.matvec(x, y),
        |x, y| m.matvec(x, y),
        |x, y| y.copy_from_slice(&lu.solve(x)),
        n,
        4,
        &opts,
    )
    .unwrap();
    let mut all: Vec<f64> = (1..=n).map(exact).collect();
    all.sort_by(|a, b| (a - sigma).abs().total_cmp(&(b - sigma).abs()));
    let mut got: Vec<f64> = r.values.iter().map(|v| v.0).collect();
    got.sort_by(|a, b| (a - sigma).abs().total_cmp(&(b - sigma).abs()));
    for i in 0..4 {
        assert!((got[i] - all[i]).abs() < 1e-6 * all[i], "{i}: {} vs {}", got[i], all[i]);
        assert!(r.residuals[i] < 1e-5 * all[i], "residual {}", r.residuals[i]);
    }
    // B-orthonormal eigenvectors
    for i in 0..4 {
        let mut bx = vec![0.0; n];
        m.matvec(&r.vectors_re[i], &mut bx);
        let nb: f64 = r.vectors_re[i].iter().zip(&bx).map(|(a, b)| a * b).sum();
        assert!((nb - 1.0).abs() < 1e-8);
    }
}

#[test]
fn generalised_nonsymmetric_matches_dense() {
    // A = K + convection, B = FE mass: compare with the dense eigenvalues of B^{-1} A
    let n = 40;
    let (k, m, _h) = fe_matrices(n);
    let mut a = Mat::zeros(n, n);
    let mut bm = Mat::zeros(n, n);
    for i in 0..n {
        for p in k.indptr[i]..k.indptr[i + 1] {
            a.set(i, k.indices[p], k.values[p]);
        }
        for p in m.indptr[i]..m.indptr[i + 1] {
            bm.set(i, m.indices[p], m.values[p]);
        }
        if i + 1 < n {
            a.set(i, i + 1, a.at(i, i + 1) + 3.0);
        }
        if i > 0 {
            a.set(i, i - 1, a.at(i, i - 1) - 3.0);
        }
    }
    let binv = dense::lu_factor(&bm).unwrap().inverse();
    let c = binv.matmul(&a).unwrap();
    let mut ev = dense::eigenvalues(&c).unwrap();
    let sigma = 500.0;
    ev.sort_by(|p, q| {
        let dp = (p.0 - sigma).hypot(p.1);
        let dq = (q.0 - sigma).hypot(q.1);
        dp.total_cmp(&dq)
    });
    let mut s = a.clone();
    for i in 0..n {
        for j in 0..n {
            s.set(i, j, a.at(i, j) - sigma * bm.at(i, j));
        }
    }
    let slu = dense::lu_factor(&s).unwrap();
    let opts = GenEigsOptions {
        sigma,
        block: 1,
        eigs: EigsOptions { tol: 1e-11, ncv: 30, ..EigsOptions::default() },
    };
    let r = eigs_gen_shift_invert(
        |x, y| y.copy_from_slice(&a.matvec(x)),
        |x, y| y.copy_from_slice(&bm.matvec(x)),
        |x, y| y.copy_from_slice(&slu.solve(x)),
        n,
        3,
        &opts,
    )
    .unwrap();
    for v in &r.values {
        let dmin = ev.iter().map(|e| (e.0 - v.0).hypot(e.1 - v.1)).fold(f64::INFINITY, f64::min);
        assert!(dmin < 1e-6 * v.0.abs().max(1.0), "{v:?} not in the dense spectrum (dist {dmin})");
    }
    for i in 0..3 {
        assert!(r.residuals[i] < 1e-5 * r.values[i].0.abs().max(1.0), "{}", r.residuals[i]);
    }
}
