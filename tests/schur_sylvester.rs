//! Residual, scale, and reuse tests for the Sylvester and Lyapunov solvers.

use rssn::kernels::dense::DenseError;
use rssn::kernels::dense::Mat;
use rssn::kernels::schur::SylvesterPlan;
use rssn::kernels::schur::lyapunov;
use rssn::kernels::schur::sylvester;

fn matrix<const N: usize, const M: usize>(rows: [[f64; M]; N]) -> Mat {
    Mat {
        rows: N,
        cols: M,
        data: rows.into_iter().flatten().collect(),
    }
}

fn scale(
    a: &Mat,
    factor: f64,
) -> Mat {
    Mat {
        rows: a.rows,
        cols: a.cols,
        data: a.data.iter().map(|x| x * factor).collect(),
    }
}

fn apply(
    a: &Mat,
    x: &Mat,
    b: &Mat,
) -> Mat {
    let mut rhs = Mat::zeros(a.rows, b.cols);
    for i in 0..a.rows {
        for j in 0..b.cols {
            let left: f64 = (0..a.cols).map(|k| a.at(i, k) * x.at(k, j)).sum();
            let right: f64 = (0..b.rows).map(|k| x.at(i, k) * b.at(k, j)).sum();
            rhs.set(i, j, left + right);
        }
    }
    rhs
}

fn assert_relative(
    actual: &Mat,
    expected: &Mat,
    tolerance: f64,
) {
    assert_eq!((actual.rows, actual.cols), (expected.rows, expected.cols));
    let scale = expected.data.iter().fold(0.0_f64, |m, &x| m.max(x.abs()));
    let error = actual
        .data
        .iter()
        .zip(&expected.data)
        .fold(0.0_f64, |m, (&x, &y)| m.max((x - y).abs()));
    assert!(
        actual.is_finite() && error <= tolerance * scale,
        "relative error {}",
        error / scale
    );
}

fn dense_nonsymmetric_a() -> Mat {
    let t = matrix([[-2.0, -0.7, 0.3], [0.7, -2.0, -0.25], [0.0, 0.0, -4.0]]);
    let q = matrix([[0.8, 0.48, 0.36], [-0.6, 0.64, 0.48], [0.0, -0.6, 0.8]]);
    q.matmul(&t).unwrap().matmul(&q.transpose()).unwrap()
}

#[test]
fn reusable_sylvester_solves_multiple_rectangular_right_hand_sides() {
    let a = dense_nonsymmetric_a();
    let b = matrix([[2.0, -0.4], [0.4, 2.0]]);
    let plan = SylvesterPlan::new(&a, &b).unwrap();
    let malformed = Mat {
        rows: 3,
        cols: 2,
        data: vec![],
    };
    assert!(matches!(plan.solve(&malformed), Err(DenseError::Shape)));
    for expected in [
        matrix([[1.0, -2.0], [-0.5, 0.125], [2.25, -3.25]]),
        matrix([[0.2, 1.0], [1.0, 1.0], [0.3, -0.4]]),
    ] {
        let c = apply(&a, &expected, &b);
        let x = plan.solve(&c).unwrap();
        assert_relative(&x, &expected, 5e-12);
        assert_relative(&apply(&a, &x, &b), &c, 5e-12);
        assert_eq!(x, sylvester(&a, &b, &c).unwrap());
    }
}

#[test]
fn actual_schur_blocks_are_retained_at_small_and_extreme_scales() {
    let a0 = dense_nonsymmetric_a();
    let b0 = matrix([[3.0, -0.75], [0.75, 3.0]]);
    let expected = matrix([[1.0, -2.0], [0.5, 3.0], [0.75, 0.25]]);
    for factor in [1e-250, 1e-24, 1e24, 1e250] {
        let (a, b) = (scale(&a0, factor), scale(&b0, factor));
        let c = apply(&a, &expected, &b);
        let x = sylvester(&a, &b, &c).unwrap();
        assert_relative(&x, &expected, 5e-11);
        assert_relative(&apply(&a, &x, &b), &c, 5e-11);
    }

    let a = matrix([[2e-20, -1e-20], [1e-20, 2e-20]]);
    let b = matrix([[3e-20]]);
    let expected = matrix([[1.0], [-2.0]]);
    let c = apply(&a, &expected, &b);
    assert_relative(&sylvester(&a, &b, &c).unwrap(), &expected, 1e-14);
}

#[test]
fn representable_solutions_survive_large_coefficient_sums_and_pivot_growth() {
    let s = matrix([[1e308]]);
    assert_eq!(sylvester(&s, &s, &s).unwrap(), matrix([[0.5]]));

    let a = matrix([[1e308, 1e308], [-1e308, 1e308]]);
    let b = matrix([[0.0]]);
    let expected = matrix([[0.5], [0.2]]);
    let c = apply(&a, &expected, &b);
    let x = sylvester(&a, &b, &c).unwrap();
    assert_relative(&x, &expected, 1e-14);
    assert_relative(&apply(&a, &x, &b), &c, 1e-14);
}

#[test]
fn lyapunov_reuses_one_schur_decomposition_with_the_existing_sign() {
    let a = dense_nonsymmetric_a();
    let at = a.transpose();
    let plan = SylvesterPlan::for_lyapunov(&a).unwrap();
    for expected in [
        matrix([[2.0, 0.2, -0.4], [0.2, 1.0, 0.25], [-0.4, 0.25, 1.5]]),
        matrix([[1.0, 0.3, 0.5], [-0.25, 2.0, 0.7], [0.4, -0.1, 3.0]]),
    ] {
        let q = apply(&a, &expected, &at);
        let x = plan.solve(&q).unwrap();
        assert_relative(&x, &expected, 5e-12);
        assert_relative(&apply(&a, &x, &at), &q, 5e-12);
        assert_eq!(x, lyapunov(&a, &q).unwrap());
    }
}

#[test]
fn shape_finiteness_singularity_and_empty_systems_are_checked() {
    let malformed = Mat {
        rows: 2,
        cols: 2,
        data: vec![1.0],
    };
    let rhs = Mat::identity(2);
    assert!(matches!(lyapunov(&malformed, &rhs), Err(DenseError::Shape)));
    assert!(matches!(
        SylvesterPlan::new(&rhs, &malformed),
        Err(DenseError::Shape)
    ));
    assert!(matches!(
        SylvesterPlan::for_lyapunov(&matrix([[f64::NAN]])),
        Err(DenseError::NonFinite)
    ));
    let a = matrix([[1.0]]);
    let b = matrix([[-1.0]]);
    let singular = SylvesterPlan::new(&a, &b).unwrap();
    assert!(matches!(singular.solve(&a), Err(DenseError::Singular)));
    assert!(matches!(
        singular.solve(&matrix([[f64::INFINITY]])),
        Err(DenseError::NonFinite)
    ));
    let tiny = matrix([[f64::from_bits(1)]]);
    assert!(matches!(
        sylvester(&tiny, &matrix([[0.0]]), &a),
        Err(DenseError::NonFinite)
    ));

    let empty = Mat::zeros(0, 0);
    let empty_rhs = Mat::zeros(0, 1);
    let plan = SylvesterPlan::new(&empty, &a).unwrap();
    assert_eq!(plan.solve(&empty_rhs).unwrap(), empty_rhs);
    assert_eq!(sylvester(&empty, &a, &empty_rhs).unwrap(), empty_rhs);
    assert_eq!(lyapunov(&empty, &empty).unwrap(), empty);
}

#[test]
fn cancelling_a_solve_leaves_the_reduction_reusable() {
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    let a = matrix([[2.0]]);
    let b = matrix([[3.0]]);
    let c = matrix([[10.0]]);
    let plan = SylvesterPlan::new(&a, &b).unwrap();
    {
        let _scope = rssn::cancellation::install(Arc::new(AtomicBool::new(true)));
        assert!(matches!(plan.solve(&c), Err(DenseError::Cancelled)));
        assert!(matches!(
            SylvesterPlan::new(&a, &b),
            Err(DenseError::Cancelled)
        ));
    }
    assert_eq!(plan.solve(&c).unwrap(), a);
}
