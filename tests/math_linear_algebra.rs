//! Dense linear algebra on the inputs numerical libraries usually get wrong:
//! empty and 1x1 matrices, non-finite entries, rank deficiency, extreme
//! scales, ill-conditioning, defective and complex-spectrum matrices, and
//! the matrix functions on singular and symmetric input.

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::dense::*;
use rssn::kernels::matfun::expm;
use rssn::kernels::matfun::logm;
use rssn::kernels::matfun::sqrtm;

fn m(rows: &[Vec<f64>]) -> Mat {
    Mat::from_rows(rows).unwrap()
}

fn residual(
    a: &Mat,
    b: &Mat,
) -> f64 {
    a.data
        .iter()
        .zip(&b.data)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn hilbert(n: usize) -> Mat {
    Mat::from_rows(
        &(0..n)
            .map(|i| (0..n).map(|j| 1.0 / (i + j + 1) as f64).collect())
            .collect::<Vec<Vec<f64>>>(),
    )
    .unwrap()
}

#[test]
fn non_finite_input_is_an_error_everywhere() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let a = m(&[vec![bad, 1.0], vec![1.0, 2.0]]);
        let b = [1.0, 1.0];
        assert_eq!(lu_factor(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(solve(&a, &b).err(), Some(DenseError::NonFinite));
        assert_eq!(cholesky(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(cholesky_solve(&a, &b).err(), Some(DenseError::NonFinite));
        assert_eq!(lstsq(&a, &b).err(), Some(DenseError::NonFinite));
        assert_eq!(svd(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(pinv(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(rank(&a, None).err(), Some(DenseError::NonFinite));
        assert_eq!(cond(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(eigen_symmetric(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(eigenvalues(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(eigen(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(null_space(&a, None).err(), Some(DenseError::NonFinite));
        assert_eq!(expm(&a).err().is_some(), true);
        assert_eq!(sqrtm(&a).err(), Some(DenseError::NonFinite));
        assert_eq!(logm(&a).err(), Some(DenseError::NonFinite));
        assert!(eigen_generalized_symmetric(&a, &Mat::identity(2)).is_err());
        // A non-finite right-hand side too.
        let ok = m(&[vec![2.0, 0.0], vec![0.0, 2.0]]);
        assert_eq!(lstsq(&ok, &[bad, 1.0]).err(), Some(DenseError::NonFinite));
    }
}

#[test]
fn wrong_lengths_are_errors_not_panics() {
    let a = Mat::identity(3);
    assert_eq!(solve(&a, &[1.0]).err(), Some(DenseError::Shape));
    assert_eq!(
        cholesky_solve(&a, &[1.0, 2.0]).err(),
        Some(DenseError::Shape)
    );
    assert_eq!(
        cholesky_solve(&a, &[1.0, 2.0, 3.0, 4.0]).err(),
        Some(DenseError::Shape)
    );
    assert_eq!(lstsq(&a, &[1.0]).err(), Some(DenseError::Shape));
    assert_eq!(lu_factor(&Mat::zeros(2, 3)).err(), Some(DenseError::Shape));
    assert_eq!(cholesky(&Mat::zeros(2, 3)).err(), Some(DenseError::Shape));
    assert_eq!(sqrtm(&Mat::zeros(2, 3)).err(), Some(DenseError::Shape));
    assert!(Mat::from_rows(&[vec![1.0], vec![1.0, 2.0]]).is_err());
}

#[test]
fn empty_and_one_by_one() {
    let e = Mat::zeros(0, 0);
    assert!(lu_factor(&e).is_ok());
    assert!(cholesky(&e).is_ok());
    assert_eq!(eigenvalues(&e).unwrap(), vec![]);
    assert_eq!(rank(&e, None).unwrap(), 0);
    assert_eq!(pinv(&Mat::zeros(0, 3)).unwrap().rows, 3);
    assert_eq!(pinv(&Mat::zeros(0, 3)).unwrap().cols, 0);
    assert_eq!(
        null_space(&Mat::zeros(0, 3), None).unwrap(),
        Mat::identity(3)
    );
    assert_eq!(null_space(&Mat::zeros(3, 0), None).unwrap().cols, 0);
    let (q, r) = qr(&Mat::zeros(0, 4));
    assert_eq!((q.rows, r.rows, r.cols), (0, 0, 4));
    let one = m(&[vec![-3.0]]);
    assert_eq!(solve(&one, &[6.0]).unwrap(), vec![-2.0]);
    assert_eq!(lu_factor(&one).unwrap().det(), -3.0);
    assert!(cholesky(&one).is_err());
    assert_eq!(svd(&one).unwrap().s, vec![3.0]);
    assert_eq!(eigenvalues(&one).unwrap(), vec![(-3.0, 0.0)]);
    assert_eq!(pinv(&m(&[vec![0.0]])).unwrap().data, vec![0.0]);
    assert!(sqrtm(&one).is_err(), "no real square root of -3");
    assert!(logm(&one).is_err());
    let e1 = expm(&one).unwrap();
    assert!((e1.at(0, 0) - (-3f64).exp()).abs() < 1e-15);
}

#[test]
fn singular_systems_are_detected() {
    let a = m(&[vec![1.0, 2.0], vec![2.0, 4.0]]);
    assert_eq!(solve(&a, &[1.0, 1.0]).err(), Some(DenseError::Singular));
    assert_eq!(
        lstsq(
            &m(&[vec![1.0, 2.0], vec![2.0, 4.0], vec![3.0, 6.0]]),
            &[1.0, 2.0, 3.0]
        )
        .err(),
        Some(DenseError::Singular)
    );
    assert_eq!(
        cholesky(&m(&[vec![1.0, 1.0], vec![1.0, 1.0]])).err(),
        Some(DenseError::Singular)
    );
    assert!(cond(&a).unwrap().is_infinite() || cond(&a).unwrap() > 1e15);
}

#[test]
fn rank_deficient_pseudo_inverse_satisfies_the_penrose_conditions() {
    let a = m(&[
        vec![1.0, 2.0, 3.0],
        vec![2.0, 4.0, 6.0],
        vec![1.0, 0.0, 1.0],
        vec![0.0, 1.0, 1.0],
    ]);
    assert_eq!(rank(&a, None).unwrap(), 2);
    let p = pinv(&a).unwrap();
    let apa = a.matmul(&p).unwrap().matmul(&a).unwrap();
    let pap = p.matmul(&a).unwrap().matmul(&p).unwrap();
    assert!(residual(&apa, &a) < 1e-12);
    assert!(residual(&pap, &p) < 1e-12);
    let ap = a.matmul(&p).unwrap();
    assert!(residual(&ap, &ap.transpose()) < 1e-12);
    let pa = p.matmul(&a).unwrap();
    assert!(residual(&pa, &pa.transpose()) < 1e-12);
    // Minimum-norm solution of a consistent rank-deficient system.
    let b = [1.0, 2.0, 3.0];
    let line = m(&[vec![1.0, 2.0], vec![2.0, 4.0], vec![3.0, 6.0]]);
    let x = lstsq_svd(&line, &b).unwrap();
    assert!(
        (x[0] - 0.2).abs() < 1e-12 && (x[1] - 0.4).abs() < 1e-12,
        "{x:?}"
    );
}

#[test]
fn null_space_and_range() {
    // Wide, rank 2 of 4 columns.
    let a = m(&[vec![1.0, 1.0, 0.0, 2.0], vec![0.0, 1.0, 1.0, 1.0]]);
    let n = null_space(&a, None).unwrap();
    assert_eq!((n.rows, n.cols), (4, 2));
    assert!(a.matmul(&n).unwrap().data.iter().all(|v| v.abs() < 1e-12));
    let gram = n.transpose().matmul(&n).unwrap();
    assert!(
        residual(&gram, &Mat::identity(2)) < 1e-12,
        "orthonormal basis"
    );
    // Tall, full column rank: trivial null space.
    let t = m(&[vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0]]);
    assert_eq!(null_space(&t, None).unwrap().cols, 0);
    // Zero matrix: everything is null.
    assert_eq!(null_space(&Mat::zeros(2, 3), None).unwrap().cols, 3);
    // Range of a rank-1 matrix is one-dimensional and contains the columns.
    let r1 = m(&[vec![1.0, 2.0], vec![2.0, 4.0], vec![3.0, 6.0]]);
    let q = orth(&r1, None).unwrap();
    assert_eq!(q.cols, 1);
    let proj = q.matmul(&q.transpose()).unwrap().matmul(&r1).unwrap();
    assert!(residual(&proj, &r1) < 1e-12);
    // A tolerance that is too big drops real directions.
    assert_eq!(null_space(&Mat::identity(3), Some(2.0)).unwrap().cols, 3);
}

#[test]
fn extreme_scales_keep_their_rank_and_values() {
    for scale in [1e-300, 1e-150, 1.0, 1e150, 1e300] {
        let a = m(&[
            vec![3.0 * scale, 1.0 * scale],
            vec![1.0 * scale, 2.0 * scale],
        ]);
        let s = svd(&a).unwrap().s;
        let want = [
            (5.0 + 5f64.sqrt()) / 2.0 * scale,
            (5.0 - 5f64.sqrt()) / 2.0 * scale,
        ];
        for (g, w) in s.iter().zip(want) {
            assert!((g - w).abs() <= 1e-13 * w, "scale {scale}: {g} vs {w}");
        }
        assert_eq!(rank(&a, None).unwrap(), 2, "scale {scale}");
        let ev = eigenvalues(&a).unwrap();
        assert!(
            ev.iter().all(|&(re, im)| im == 0.0 && re > 0.0),
            "scale {scale}: {ev:?}"
        );
    }
}

#[test]
fn ill_conditioned_hilbert_matrices() {
    let h = hilbert(12);
    let c = cond(&h).unwrap();
    assert!(c > 1e15 && c < 1e19, "cond(H12) = {c:e}");
    // Solve H x = H * ones: the error is bounded by cond * eps.
    let ones = vec![1.0; 12];
    let b = h.matvec(&ones);
    let x = solve(&h, &b).unwrap();
    let err = x.iter().map(|v| (v - 1.0).abs()).fold(0.0, f64::max);
    assert!(err.is_finite() && err < 1.0, "forward error {err}"); // cond ~ 1e16: only the order of magnitude survives
    // The pseudo-inverse with the default tolerance truncates instead of
    // blowing up.
    let p = pinv(&h).unwrap();
    assert!(p.data.iter().all(|v| v.is_finite()));
    // H5 is solved to the expected accuracy.
    let h5 = hilbert(5);
    let b5 = h5.matvec(&vec![1.0; 5]);
    let x5 = solve(&h5, &b5).unwrap();
    assert!(x5.iter().all(|v| (v - 1.0).abs() < 1e-9));
}

#[test]
fn eigenvalues_of_hard_matrices() {
    // Rotation: purely imaginary pair.
    let rot = m(&[vec![0.0, -1.0], vec![1.0, 0.0]]);
    let ev = eigenvalues(&rot).unwrap();
    assert!(
        ev.iter()
            .all(|&(re, im)| re.abs() < 1e-14 && (im.abs() - 1.0).abs() < 1e-14),
        "{ev:?}"
    );
    // Jordan block: double eigenvalue, one eigenvector.
    let jordan = m(&[vec![1.0, 1.0], vec![0.0, 1.0]]);
    let ev = eigenvalues(&jordan).unwrap();
    assert!(
        ev.iter()
            .all(|&(re, im)| (re - 1.0).abs() < 1e-7 && im.abs() < 1e-7),
        "{ev:?}"
    );
    // Badly scaled: eigenvalues 0 and 2 (balancing needed).
    let bad = m(&[vec![1.0, 1e12], vec![1e-12, 1.0]]);
    let mut re: Vec<f64> = eigenvalues(&bad).unwrap().iter().map(|e| e.0).collect();
    re.sort_by(f64::total_cmp);
    assert!(re[0].abs() < 1e-8 && (re[1] - 2.0).abs() < 1e-8, "{re:?}");
    // Companion matrix of (x - 1)(x - 2)(x - 3)(x - 4)(x - 5).
    let c = m(&[
        vec![15.0, -85.0, 225.0, -274.0, 120.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0],
    ]);
    let mut re: Vec<f64> = eigenvalues(&c).unwrap().iter().map(|e| e.0).collect();
    re.sort_by(f64::total_cmp);
    for (k, v) in re.iter().enumerate() {
        assert!((v - (k + 1) as f64).abs() < 1e-8, "{re:?}");
    }
    // Zero matrix and identity.
    assert!(
        eigenvalues(&Mat::zeros(3, 3))
            .unwrap()
            .iter()
            .all(|e| *e == (0.0, 0.0))
    );
    assert!(
        eigenvalues(&Mat::identity(4))
            .unwrap()
            .iter()
            .all(|e| (e.0 - 1.0).abs() < 1e-15)
    );
    // Symmetric: ascending and orthonormal vectors.
    let s = m(&[
        vec![2.0, -1.0, 0.0],
        vec![-1.0, 2.0, -1.0],
        vec![0.0, -1.0, 2.0],
    ]);
    let (vals, vecs) = eigen_symmetric(&s).unwrap();
    assert!(vals.windows(2).all(|w| w[0] <= w[1]));
    assert!((vals[0] - (2.0 - 2f64.sqrt())).abs() < 1e-13);
    assert!(residual(&vecs.transpose().matmul(&vecs).unwrap(), &Mat::identity(3)) < 1e-13);
}

#[test]
fn matrix_functions_on_singular_symmetric_and_nilpotent_input() {
    // sqrtm of a singular PSD matrix (the Denman-Beavers iteration cannot
    // factor it).
    let psd = m(&[vec![1.0, 0.0], vec![0.0, 0.0]]);
    let r = sqrtm(&psd).unwrap();
    assert!(residual(&r, &psd) < 1e-14);
    let psd2 = m(&[vec![2.0, 2.0], vec![2.0, 2.0]]);
    let r = sqrtm(&psd2).unwrap();
    assert!(residual(&r.matmul(&r).unwrap(), &psd2) < 1e-13);
    assert!(
        sqrtm(&Mat::zeros(3, 3))
            .unwrap()
            .data
            .iter()
            .all(|v| v.abs() < 1e-15)
    );
    // Negative eigenvalue: no real principal root.
    assert!(sqrtm(&m(&[vec![-4.0, 0.0], vec![0.0, 1.0]])).is_err());
    // Non-symmetric still works: sqrt([[4, 1], [0, 9]])^2.
    let u = m(&[vec![4.0, 1.0], vec![0.0, 9.0]]);
    let r = sqrtm(&u).unwrap();
    assert!(residual(&r.matmul(&r).unwrap(), &u) < 1e-12);
    // logm of a symmetric positive definite matrix is exact; of a singular
    // or indefinite one it is an error.
    let spd = m(&[vec![2.0, 1.0], vec![1.0, 2.0]]);
    let l = logm(&spd).unwrap();
    assert!(residual(&expm(&l).unwrap(), &spd) < 1e-12);
    assert!(logm(&psd).is_err());
    assert!(logm(&m(&[vec![-1.0, 0.0], vec![0.0, 1.0]])).is_err());
    // expm of a nilpotent matrix is a finite series: exp(N) = I + N + N^2/2.
    let n = m(&[
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![0.0, 0.0, 0.0],
    ]);
    let e = expm(&n).unwrap();
    let want = m(&[
        vec![1.0, 1.0, 0.5],
        vec![0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0],
    ]);
    assert!(residual(&e, &want) < 1e-14);
    // expm(0) = I; large norm is scaled and squared.
    assert!(residual(&expm(&Mat::zeros(3, 3)).unwrap(), &Mat::identity(3)) < 1e-15);
    let big = m(&[vec![0.0, 40.0], vec![-40.0, 0.0]]);
    let e = expm(&big).unwrap();
    assert!((e.at(0, 0) - 40f64.cos()).abs() < 1e-10 && (e.at(0, 1) - 40f64.sin()).abs() < 1e-10);
    // exp(A) exp(-A) = I for a non-normal matrix.
    let a = m(&[vec![1.0, 100.0], vec![0.0, -1.0]]);
    let prod = expm(&a)
        .unwrap()
        .matmul(&expm(&m(&[vec![-1.0, -100.0], vec![0.0, 1.0]])).unwrap())
        .unwrap();
    assert!(residual(&prod, &Mat::identity(2)) < 1e-9);
    // An overflowing exponential is an error, not Inf entries.
    assert!(expm(&m(&[vec![1e4, 0.0], vec![0.0, 1.0]])).is_err());
}

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x1157),
        failure_persistence: None,
        cases: 200,
        ..ProptestConfig::default()
    }
}

fn entry() -> impl Strategy<Value = f64> {
    prop_oneof![
        8 => -10.0..10.0f64,
        1 => Just(0.0),
        1 => Just(f64::NAN),
        1 => Just(f64::INFINITY),
        1 => Just(1e300),
        1 => Just(-1e-300),
    ]
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn dense_routines_never_panic(rows in 0usize..5, cols in 0usize..5, data in prop::collection::vec(entry(), 25), b in prop::collection::vec(entry(), 0..6)) {
        let a = Mat { rows, cols, data: data[..rows * cols].to_vec() };
        let _ = lu_factor(&a);
        let _ = solve(&a, &b);
        let _ = cholesky(&a);
        let _ = cholesky_solve(&a, &b);
        let _ = qr(&a);
        let _ = lstsq(&a, &b);
        let _ = lstsq_svd(&a, &b);
        let _ = svd(&a);
        let _ = pinv(&a);
        let _ = rank(&a, None);
        let _ = cond(&a);
        let _ = null_space(&a, None);
        let _ = orth(&a, None);
        let _ = eigen_symmetric(&a);
        let _ = eigenvalues(&a);
        let _ = eigen(&a);
        let _ = spectral_radius(&a);
        let _ = expm(&a);
        let _ = sqrtm(&a);
        let _ = logm(&a);
        let _ = eigen_generalized_symmetric(&a, &a);
    }

    #[test]
    fn pinv_solves_consistent_systems(rows in 1usize..5, cols in 1usize..5, seed in prop::collection::vec(-5.0..5.0f64, 25), x in prop::collection::vec(-5.0..5.0f64, 5)) {
        let a = Mat { rows, cols, data: seed[..rows * cols].to_vec() };
        let b = a.matvec(&x[..cols]);
        let sol = lstsq_svd(&a, &b).unwrap();
        // The residual of the least-squares solution must vanish for a
        // consistent system, regardless of rank.
        let back = a.matvec(&sol);
        let scale = 1.0 + b.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let cn = cond(&a).unwrap_or(f64::INFINITY);
        if cn < 1e8 || rank(&a, None).unwrap() < rows.min(cols) {
            for (u, v) in back.iter().zip(&b) {
                prop_assert!((u - v).abs() <= 1e-6 * scale, "{u} vs {v}");
            }
        }
    }
}
