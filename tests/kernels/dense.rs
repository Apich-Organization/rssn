//! Tests for the dense linear algebra kernels.

use rssn::kernels::dense::{
    Mat, cholesky, cholesky_solve, cond, eigen_symmetric, eigenvalues, lstsq, lstsq_svd,
    lu_factor, pinv, qr, rank, solve, spectral_radius, svd,
};

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
fn cholesky_factor() {
    let a = m(&[&[4.0, 12.0, -16.0], &[12.0, 37.0, -43.0], &[-16.0, -43.0, 98.0]]);
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
        data: (0..16).map(|k| 1.0 / f64::from(k / 4 + k % 4 + 1)).collect(),
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
