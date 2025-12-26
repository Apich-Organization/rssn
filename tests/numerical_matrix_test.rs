use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;
use rssn::prelude::*;

#[test]

fn test_matrix_basic() {

    let data = vec![1.0, 2.0, 3.0, 4.0];

    let m = numerical_Matrix::new(2, 2, data.clone());

    assert_eq!(m.rows(), 2);

    assert_eq!(m.cols(), 2);

    assert_eq!(m.data(), &data);

    assert_eq!(*m.get(0, 1), 2.0);
}

#[test]

fn test_matrix_arithmetic() {

    let m1 = numerical_Matrix::new(2, 2, vec![1.0, 2.0, 3.0, 4.0]);

    let m2 = numerical_Matrix::new(2, 2, vec![5.0, 6.0, 7.0, 8.0]);

    // Add
    let m3 = m1.clone() + m2.clone();

    assert_eq!(m3.data(), &vec![6.0, 8.0, 10.0, 12.0]);

    // Sub
    let m4 = m2.clone() - m1.clone();

    assert_eq!(m4.data(), &vec![4.0, 4.0, 4.0, 4.0]);

    // Mul
    let m5 = m1.clone() * m2.clone();

    // [1 2] [5 6] = [1*5+2*7 1*6+2*8] = [19 22]
    // [3 4] [7 8]   [3*5+4*7 3*6+4*8] = [43 50]
    assert_eq!(m5.data(), &vec![19.0, 22.0, 43.0, 50.0]);
}

#[test]

fn test_matrix_transpose() {

    let m = numerical_Matrix::new(
        2,
        3,
        vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        ],
    );

    let mt = m.transpose();

    assert_eq!(mt.rows(), 3);

    assert_eq!(mt.cols(), 2);

    assert_eq!(mt.data(), &vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
}

#[test]

fn test_matrix_determinant() {

    let m = numerical_Matrix::new(2, 2, vec![1.0, 2.0, 3.0, 4.0]);

    assert_approx_eq!(
        m.determinant()
            .unwrap(),
        -2.0
    );

    let m3 = numerical_Matrix::new(
        3,
        3,
        vec![
            1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0,
        ],
    );

    // det = 1*(1*0-4*6) - 2*(0*0-4*5) + 3*(0*6-1*5)
    // det = -24 + 40 - 15 = 1
    assert_approx_eq!(
        m3.determinant()
            .unwrap(),
        1.0
    );
}

#[test]

fn test_matrix_inverse() {

    let m = numerical_Matrix::new(2, 2, vec![4.0, 7.0, 2.0, 6.0]);

    // det = 24 - 14 = 10
    // inv = 1/10 * [6 -7; -2 4] = [0.6 -0.7; -0.2 0.4]
    let inv = m.inverse().unwrap();

    assert_approx_eq!(inv.get(0, 0), 0.6);

    assert_approx_eq!(inv.get(0, 1), -0.7);

    assert_approx_eq!(inv.get(1, 0), -0.2);

    assert_approx_eq!(inv.get(1, 1), 0.4);

    let identity = m * inv;

    assert_approx_eq!(identity.get(0, 0), 1.0);

    assert_approx_eq!(identity.get(0, 1), 0.0);

    assert_approx_eq!(identity.get(1, 0), 0.0);

    assert_approx_eq!(identity.get(1, 1), 1.0);
}

#[test]

fn test_matrix_norms() {

    let m = numerical_Matrix::new(2, 2, vec![1.0, -2.0, -3.0, 4.0]);

    assert_approx_eq!(m.frobenius_norm(), (1.0 + 4.0 + 9.0 + 16.0f64).sqrt());

    assert_eq!(m.l1_norm(), 6.0);
}

#[test]

fn test_matrix_trace() {

    let m = numerical_Matrix::new(
        3,
        3,
        vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
        ],
    );

    assert_eq!(m.trace().unwrap(), 15.0);
}

#[test]

fn test_matrix_rank() {

    let m = numerical_Matrix::new(
        3,
        3,
        vec![
            1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.0, 1.0, 1.0,
        ],
    );

    assert_eq!(m.rank().unwrap(), 2);
}

#[test]

fn test_matrix_identity_orthogonal() {

    let m = numerical_Matrix::<f64>::identity(3);

    assert!(m.is_identity(1e-9));

    assert!(m.is_orthogonal(1e-9));

    let m2 = numerical_Matrix::new(2, 2, vec![0.0, 1.0, 1.0, 0.0]);

    assert!(!m2.is_identity(1e-9));

    assert!(m2.is_orthogonal(1e-9));
}

proptest! {
    #[test]
    fn prop_matrix_transpose_transpose(rows in 1..10usize, cols in 1..10usize) {
        let data: Vec<f64> = vec![0.0; rows * cols]; // Simplified for prop
        let m = numerical_Matrix::new(rows, cols, data);
        assert_eq!(m.transpose().transpose(), m);
    }
}
