//! Tests for expm, sqrtm, logm, phi functions and exp(tA) v.

use rssn::kernels::dense::Mat;
use rssn::kernels::dense::{self};
use rssn::kernels::krylov::Csr;
use rssn::kernels::matfun::expm;
use rssn::kernels::matfun::expm_krylov;
use rssn::kernels::matfun::expm_multiply;
use rssn::kernels::matfun::logm;
use rssn::kernels::matfun::phi_matrices;
use rssn::kernels::matfun::sqrtm;

fn maxdiff(
    a: &Mat,
    b: &Mat,
) -> f64 {
    a.data
        .iter()
        .zip(&b.data)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn pseudo(
    n: usize,
    seed: u64,
    scale: f64,
) -> Mat {
    let mut s = seed;
    let mut m = Mat::zeros(n, n);
    for v in &mut m.data {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *v = scale * (((s >> 33) as f64) / f64::from(1u32 << 31) - 0.5);
    }
    m
}

#[test]
fn expm_closed_forms() {
    // rotation generator
    let th = 0.9;
    let a = Mat::from_rows(&[vec![0.0, -th], vec![th, 0.0]]).unwrap();
    let e = expm(&a).unwrap();
    assert!(
        maxdiff(
            &e,
            &Mat::from_rows(&[vec![th.cos(), -th.sin()], vec![th.sin(), th.cos()]]).unwrap()
        ) < 1e-14
    );
    // large rotation exercises scaling and squaring
    let a = Mat::from_rows(&[vec![0.0, -40.0], vec![40.0, 0.0]]).unwrap();
    let e = expm(&a).unwrap();
    assert!((e.at(0, 0) - 40f64.cos()).abs() < 1e-11 && (e.at(1, 0) - 40f64.sin()).abs() < 1e-11);
    // nilpotent: exp(N) = I + N + N^2/2
    let n = Mat::from_rows(&[
        vec![0.0, 1.0, 2.0],
        vec![0.0, 0.0, 3.0],
        vec![0.0, 0.0, 0.0],
    ])
    .unwrap();
    let e = expm(&n).unwrap();
    let want = Mat::from_rows(&[
        vec![1.0, 1.0, 3.5],
        vec![0.0, 1.0, 3.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    assert!(maxdiff(&e, &want) < 1e-14);
    // diagonal
    let d = Mat::from_rows(&[vec![-2.0, 0.0], vec![0.0, 3.5]]).unwrap();
    let e = expm(&d).unwrap();
    assert!(
        (e.at(0, 0) - (-2f64).exp()).abs() < 1e-15 && (e.at(1, 1) - 3.5f64.exp()).abs() < 1e-12
    );
    // exp(A) exp(-A) = I for every Pade branch
    for scale in [0.005, 0.1, 0.5, 1.5, 4.0, 30.0] {
        let a = pseudo(6, 3, scale);
        let ea = expm(&a).unwrap();
        let mut neg = a.clone();
        neg.data.iter_mut().for_each(|x| *x = -*x);
        let prod = ea.matmul(&expm(&neg).unwrap()).unwrap();
        assert!(
            maxdiff(&prod, &Mat::identity(6)) < 1e-9 * (1.0 + scale.exp().min(1e9)),
            "scale {scale}"
        );
        // semigroup property exp(A) = exp(A/2)^2
        let mut half = a.clone();
        half.data.iter_mut().for_each(|x| *x *= 0.5);
        let eh = expm(&half).unwrap();
        let sq = eh.matmul(&eh).unwrap();
        assert!(
            maxdiff(&sq, &ea) < 1e-10 * ea.data.iter().fold(1.0_f64, |m, x| m.max(x.abs())),
            "scale {scale}"
        );
    }
    // symmetric: compare with the spectral decomposition
    let s0 = pseudo(5, 11, 2.0);
    let sym = lin_sym(&s0);
    let (vals, vecs) = dense::eigen_symmetric(&sym).unwrap();
    let mut spec = Mat::zeros(5, 5);
    for i in 0..5 {
        for j in 0..5 {
            spec.set(
                i,
                j,
                (0..5)
                    .map(|k| vecs.at(i, k) * vals[k].exp() * vecs.at(j, k))
                    .sum(),
            );
        }
    }
    assert!(maxdiff(&expm(&sym).unwrap(), &spec) < 1e-12);
    assert!(expm(&Mat::zeros(2, 3)).is_err());
}

fn lin_sym(a: &Mat) -> Mat {
    let t = a.transpose();
    let mut r = a.clone();
    for (x, y) in r.data.iter_mut().zip(&t.data) {
        *x = 0.5 * (*x + y);
    }
    r
}

#[test]
fn sqrtm_and_logm() {
    // SPD: S^2 = A
    let b = pseudo(5, 5, 1.0);
    let spd = {
        let g = b.matmul(&b.transpose()).unwrap();
        let mut r = g;
        for i in 0..5 {
            r.set(i, i, r.at(i, i) + 0.5);
        }
        r
    };
    let s = sqrtm(&spd).unwrap();
    assert!(maxdiff(&s.matmul(&s).unwrap(), &spd) < 1e-12);
    // rotation: sqrt of R(theta) = R(theta / 2), log = generator
    let th: f64 = 1.3;
    let r = Mat::from_rows(&[vec![th.cos(), -th.sin()], vec![th.sin(), th.cos()]]).unwrap();
    let half = Mat::from_rows(&[
        vec![(th / 2.0).cos(), -(th / 2.0).sin()],
        vec![(th / 2.0).sin(), (th / 2.0).cos()],
    ])
    .unwrap();
    assert!(maxdiff(&sqrtm(&r).unwrap(), &half) < 1e-13);
    let l = logm(&r).unwrap();
    assert!(
        maxdiff(
            &l,
            &Mat::from_rows(&[vec![0.0, -th], vec![th, 0.0]]).unwrap()
        ) < 1e-12
    );
    // logm(expm(A)) = A for a nonsymmetric A with small spectrum
    let a = pseudo(6, 9, 1.2);
    let la = logm(&expm(&a).unwrap()).unwrap();
    assert!(maxdiff(&la, &a) < 1e-10, "{}", maxdiff(&la, &a));
    // diagonal matrix with a wide spread of eigenvalues
    let d = Mat::from_rows(&[
        vec![1e-3, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1e4],
    ])
    .unwrap();
    let ld = logm(&d).unwrap();
    assert!(
        (ld.at(0, 0) - 1e-3f64.ln()).abs() < 1e-11 && (ld.at(2, 2) - 1e4f64.ln()).abs() < 1e-11
    );
    // exp(log(A)) = A for SPD
    assert!(maxdiff(&expm(&logm(&spd).unwrap()).unwrap(), &spd) < 1e-10);
    // eigenvalue on the negative axis is rejected
    let neg = Mat::from_rows(&[vec![-1.0, 0.0], vec![0.0, 1.0]]).unwrap();
    assert!(sqrtm(&neg).is_err() || logm(&neg).is_err());
    assert!(logm(&Mat::zeros(2, 2)).is_err());
}

#[test]
fn phi_functions() {
    // scalar: phi_1(z) = (e^z - 1) / z, phi_2(z) = (e^z - 1 - z) / z^2, phi_3(z) = (e^z - 1 - z - z^2/2) / z^3
    for &z in &[0.7_f64, -1.9, 2.5] {
        let a = Mat::from_rows(&[vec![z]]).unwrap();
        let p = phi_matrices(&a, 3).unwrap();
        assert!((p[0].at(0, 0) - z.exp()).abs() < 1e-13);
        assert!((p[1].at(0, 0) - (z.exp() - 1.0) / z).abs() < 1e-13);
        assert!((p[2].at(0, 0) - (z.exp() - 1.0 - z) / (z * z)).abs() < 1e-13);
        assert!((p[3].at(0, 0) - (z.exp() - 1.0 - z - z * z / 2.0) / z.powi(3)).abs() < 1e-12);
    }
    // phi_k(0) = 1 / k!
    let p = phi_matrices(&Mat::zeros(2, 2), 3).unwrap();
    assert!(
        (p[1].at(0, 0) - 1.0).abs() < 1e-14
            && (p[2].at(1, 1) - 0.5).abs() < 1e-14
            && (p[3].at(0, 0) - 1.0 / 6.0).abs() < 1e-14
    );
    // phi_1(A) A = exp(A) - I
    let a = pseudo(4, 2, 1.5);
    let p = phi_matrices(&a, 1).unwrap();
    let lhs = p[1].matmul(&a).unwrap();
    let mut rhs = expm(&a).unwrap();
    for i in 0..4 {
        rhs.set(i, i, rhs.at(i, i) - 1.0);
    }
    assert!(maxdiff(&lhs, &rhs) < 1e-12);
}

fn laplacian1d(n: usize) -> Csr {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, -2.0 * (n as f64 + 1.0).powi(2)));
        if i + 1 < n {
            t.push((i, i + 1, (n as f64 + 1.0).powi(2)));
            t.push((i + 1, i, (n as f64 + 1.0).powi(2)));
        }
    }
    Csr::from_triplets(n, n, &t)
}

#[test]
fn expm_action_matches_dense_and_heat_kernel() {
    let n = 30;
    let a = laplacian1d(n);
    // dense reference
    let mut ad = Mat::zeros(n, n);
    for i in 0..n {
        for k in a.indptr[i]..a.indptr[i + 1] {
            ad.set(i, a.indices[k], a.values[k]);
        }
    }
    let v: Vec<f64> = (0..n)
        .map(|i| {
            ((i + 1) as f64 * std::f64::consts::PI / (n as f64 + 1.0)).sin()
                + 0.3 * (i as f64).cos()
        })
        .collect();
    for &t in &[1e-4, 2e-3, 0.05] {
        let mut at = ad.clone();
        at.data.iter_mut().for_each(|x| *x *= t);
        let want = expm(&at).unwrap().matvec(&v);
        let got = expm_multiply(&a, &v, t);
        let err = got
            .iter()
            .zip(&want)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-11, "t={t}: {err}");
        let kr = expm_krylov(|x, y| a.matvec(x, y), &v, t, 20, 1e-13).unwrap();
        let err = kr
            .iter()
            .zip(&want)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-9, "krylov t={t}: {err}");
    }
    // first sine mode decays at its exact rate
    let mode: Vec<f64> = (0..n)
        .map(|i| ((i + 1) as f64 * std::f64::consts::PI / (n as f64 + 1.0)).sin())
        .collect();
    let lam =
        -2.0 * (n as f64 + 1.0).powi(2) * (1.0 - (std::f64::consts::PI / (n as f64 + 1.0)).cos());
    let w = expm_multiply(&a, &mode, 0.02);
    for i in 0..n {
        assert!((w[i] - (lam * 0.02).exp() * mode[i]).abs() < 1e-12);
    }
    // backward in time with a skew operator: exact rotation
    let mut t = Vec::new();
    for i in 0..4 {
        t.push((2 * i, 2 * i + 1, -1.0 - i as f64));
        t.push((2 * i + 1, 2 * i, 1.0 + i as f64));
    }
    let s = Csr::from_triplets(8, 8, &t);
    let v0 = vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
    let w = expm_multiply(&s, &v0, -2.5);
    for i in 0..4 {
        let th = -2.5 * (1.0 + i as f64);
        assert!((w[2 * i] - th.cos()).abs() < 1e-12 && (w[2 * i + 1] - th.sin()).abs() < 1e-12);
    }
    let w = expm_krylov(|x, y| s.matvec(x, y), &v0, 3.0, 8, 1e-12).unwrap();
    for i in 0..4 {
        let th = 3.0 * (1.0 + i as f64);
        assert!((w[2 * i] - th.cos()).abs() < 1e-8 && (w[2 * i + 1] - th.sin()).abs() < 1e-8);
    }
}
