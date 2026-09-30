use rssn::compute::{
    charpoly, compute, det, eigenvalues, matrix_inv, ComputeConfig,
};
use rssn::symbolic::core::Expr;
use num_bigint::BigInt;

#[test]
fn test_symbolic_matrix_inversion_2x2() {
    let a = Expr::new_variable("a");
    let b = Expr::new_variable("b");
    let c = Expr::new_variable("c");
    let d = Expr::new_variable("d");

    let mat = Expr::Matrix(vec![
        vec![a.clone(), b.clone()],
        vec![c.clone(), d.clone()],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_symbolic();

    let inv_op = matrix_inv(mat);
    let result = compute(&inv_op, &config);

    eprintln!("Symbolic matrix inversion result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].len(), 2);
            assert_eq!(rows[1].len(), 2);
        }
        other => panic!("Expected symbolic Expr::Matrix, got: {:?}", other),
    }
}

#[test]
fn test_symbolic_matrix_determinant_2x2() {
    let a = Expr::new_variable("a");
    let b = Expr::new_variable("b");
    let c = Expr::new_variable("c");
    let d = Expr::new_variable("d");

    let mat = Expr::Matrix(vec![
        vec![a.clone(), b.clone()],
        vec![c.clone(), d.clone()],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_symbolic();

    let det_op = det(mat);
    let result = compute(&det_op, &config);

    eprintln!("Symbolic determinant result: {:?}", result);
    // Determinant of [[a, b], [c, d]] is a*d - b*c
    match result {
        Expr::Sub(lhs, rhs) => {
            eprintln!("LHS: {:?}, RHS: {:?}", lhs, rhs);
            assert!(matches!(*lhs, Expr::Mul(..) | Expr::Variable(..)));
            assert!(matches!(*rhs, Expr::Mul(..) | Expr::Variable(..)));
        }
        Expr::Add(..) => {}
        other => panic!("Expected Expr::Sub/Add for symbolic determinant, got: {:?}", other),
    }
}

#[test]
fn test_symbolic_characteristic_polynomial() {
    let a = Expr::new_variable("a");
    let b = Expr::new_variable("b");
    let c = Expr::new_variable("c");
    let d = Expr::new_variable("d");

    let mat = Expr::Matrix(vec![
        vec![a.clone(), b.clone()],
        vec![c.clone(), d.clone()],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_symbolic();

    let cp_op = charpoly(mat, "lambda");
    let result = compute(&cp_op, &config);

    eprintln!("Symbolic characteristic polynomial result: {:?}", result);
    // Should be an algebraic expression representing det(A - lambda*I)
    assert!(!matches!(result, Expr::BinaryList(..)));
}

#[test]
fn test_symbolic_eigenvalues_diagonal() {
    let mat = Expr::Matrix(vec![
        vec![Expr::new_bigint(BigInt::from(2)), Expr::new_bigint(BigInt::from(0))],
        vec![Expr::new_bigint(BigInt::from(0)), Expr::new_bigint(BigInt::from(5))],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_symbolic();

    let eig_op = eigenvalues(mat);
    let result = compute(&eig_op, &config);

    eprintln!("Symbolic eigenvalues result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            let mut vals: Vec<i64> = rows
                .iter()
                .filter_map(|r| match &r[0] {
                    Expr::BigInt(b) => num_traits::ToPrimitive::to_i64(b),
                    Expr::Constant(c) => Some(*c as i64),
                    _ => None,
                })
                .collect();
            vals.sort();
            assert_eq!(vals, vec![2, 5]);
        }
        other => panic!("Expected Expr::Matrix for eigenvalues, got: {:?}", other),
    }
}

#[test]
fn test_numerical_matrix_inversion_2x2() {
    // A = [[4.0, 7.0], [2.0, 6.0]]
    // det(A) = 24 - 14 = 10
    // A^-1 = [[0.6, -0.7], [-0.2, 0.4]]
    let mat = Expr::Matrix(vec![
        vec![Expr::Constant(4.0), Expr::Constant(7.0)],
        vec![Expr::Constant(2.0), Expr::Constant(6.0)],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_numerical(1e-7);

    let inv_op = matrix_inv(mat);
    let result = compute(&inv_op, &config);

    eprintln!("Numerical matrix inversion result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].len(), 2);

            let get_val = |r: usize, c: usize| match &rows[r][c] {
                Expr::Constant(v) => *v,
                other => panic!("Expected Constant at ({},{}), got {:?}", r, c, other),
            };

            assert!((get_val(0, 0) - 0.6).abs() < 1e-6);
            assert!((get_val(0, 1) - (-0.7)).abs() < 1e-6);
            assert!((get_val(1, 0) - (-0.2)).abs() < 1e-6);
            assert!((get_val(1, 1) - 0.4).abs() < 1e-6);
        }
        other => panic!("Expected numerical Expr::Matrix, got: {:?}", other),
    }
}

#[test]
fn test_numerical_matrix_inversion_3x3() {
    // Invertible 3x3 matrix
    // A = [[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [5.0, 6.0, 0.0]]
    // det(A) = 1*(0 - 24) - 2*(0 - 20) + 3*(0 - 5) = -24 + 40 - 15 = 1
    let original_data = vec![
        vec![1.0, 2.0, 3.0],
        vec![0.0, 1.0, 4.0],
        vec![5.0, 6.0, 0.0],
    ];

    let mat = Expr::Matrix(
        original_data
            .iter()
            .map(|r| r.iter().map(|&v| Expr::Constant(v)).collect())
            .collect(),
    );

    let config = ComputeConfig::default()
        .with_matrix()
        .target_numerical(1e-7);

    let inv_op = matrix_inv(mat);
    let result = compute(&inv_op, &config);

    eprintln!("Numerical 3x3 matrix inversion result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 3);
            assert_eq!(rows[0].len(), 3);

            let mut inv_mat = vec![vec![0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    match &rows[i][j] {
                        Expr::Constant(v) => inv_mat[i][j] = *v,
                        other => panic!("Expected Constant at ({},{}), got {:?}", i, j, other),
                    }
                }
            }

            // Verify A * A^-1 == Identity(3)
            for i in 0..3 {
                for j in 0..3 {
                    let mut dot = 0.0;
                    for k in 0..3 {
                        dot += original_data[i][k] * inv_mat[k][j];
                    }
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (dot - expected).abs() < 1e-5,
                        "Mismatch at ({}, {}): got {}, expected {}",
                        i,
                        j,
                        dot,
                        expected
                    );
                }
            }
        }
        other => panic!("Expected numerical Expr::Matrix, got: {:?}", other),
    }
}

#[test]
fn test_numerical_matrix_determinant() {
    // A = [[4.0, 7.0], [2.0, 6.0]], det = 10.0
    let mat = Expr::Matrix(vec![
        vec![Expr::Constant(4.0), Expr::Constant(7.0)],
        vec![Expr::Constant(2.0), Expr::Constant(6.0)],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_numerical(1e-7);

    let det_op = det(mat);
    let result = compute(&det_op, &config);

    eprintln!("Numerical determinant result: {:?}", result);
    match result {
        Expr::Constant(val) => {
            assert!((val - 10.0).abs() < 1e-6, "Expected det = 10.0, got {}", val);
        }
        other => panic!("Expected numerical Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_eigenvalues_symmetric() {
    // Symmetric matrix [[2.0, 1.0], [1.0, 2.0]]
    // Eigenvalues are 1.0 and 3.0
    let mat = Expr::Matrix(vec![
        vec![Expr::Constant(2.0), Expr::Constant(1.0)],
        vec![Expr::Constant(1.0), Expr::Constant(2.0)],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_numerical(1e-7);

    let eig_op = eigenvalues(mat);
    let result = compute(&eig_op, &config);

    eprintln!("Numerical symmetric eigenvalues result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            let mut eigs: Vec<f64> = rows
                .iter()
                .map(|r| match &r[0] {
                    Expr::Constant(v) => *v,
                    other => panic!("Expected Constant eigenvalue, got {:?}", other),
                })
                .collect();
            eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());

            assert!((eigs[0] - 1.0).abs() < 1e-5, "Expected 1.0, got {}", eigs[0]);
            assert!((eigs[1] - 3.0).abs() < 1e-5, "Expected 3.0, got {}", eigs[1]);
        }
        other => panic!("Expected Expr::Matrix for eigenvalues, got: {:?}", other),
    }
}

#[test]
fn test_numerical_eigenvalues_general() {
    // Non-symmetric real matrix [[1.0, 2.0], [3.0, 4.0]]
    // Trace = 5.0, Det = -2.0
    // Characteristic poly: lambda^2 - 5*lambda - 2 = 0
    // Eigenvalues: (5 +- sqrt(33)) / 2 => ~5.3722813 and ~-0.3722813
    let mat = Expr::Matrix(vec![
        vec![Expr::Constant(1.0), Expr::Constant(2.0)],
        vec![Expr::Constant(3.0), Expr::Constant(4.0)],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .target_numerical(1e-7);

    let eig_op = eigenvalues(mat);
    let result = compute(&eig_op, &config);

    eprintln!("Numerical general eigenvalues result: {:?}", result);
    match result {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            let mut eigs: Vec<f64> = rows
                .iter()
                .map(|r| match &r[0] {
                    Expr::Constant(v) => *v,
                    other => panic!("Expected Constant eigenvalue, got {:?}", other),
                })
                .collect();
            eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let expected_1 = (5.0 - (33.0_f64).sqrt()) / 2.0;
            let expected_2 = (5.0 + (33.0_f64).sqrt()) / 2.0;

            assert!(
                (eigs[0] - expected_1).abs() < 1e-5,
                "Expected {}, got {}",
                expected_1,
                eigs[0]
            );
            assert!(
                (eigs[1] - expected_2).abs() < 1e-5,
                "Expected {}, got {}",
                expected_2,
                eigs[1]
            );
        }
        other => panic!("Expected Expr::Matrix for eigenvalues, got: {:?}", other),
    }
}

#[test]
fn test_matrix_with_bindings_numerical() {
    // Matrix with variables x and y:
    // [[x, 1.0], [1.0, y]] with x = 2.0, y = 2.0
    // Eigenvalues should be 1.0 and 3.0
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");

    let mat = Expr::Matrix(vec![
        vec![x, Expr::Constant(1.0)],
        vec![Expr::Constant(1.0), y],
    ]);

    let config = ComputeConfig::default()
        .with_matrix()
        .bind("x", 2.0)
        .bind("y", 2.0)
        .target_numerical(1e-7);

    let eig_op = eigenvalues(mat.clone());
    let eig_res = compute(&eig_op, &config);

    eprintln!("Bound variables eigenvalues result: {:?}", eig_res);
    match eig_res {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            let mut eigs: Vec<f64> = rows
                .iter()
                .map(|r| match &r[0] {
                    Expr::Constant(v) => *v,
                    other => panic!("Expected Constant eigenvalue, got {:?}", other),
                })
                .collect();
            eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert!((eigs[0] - 1.0).abs() < 1e-5);
            assert!((eigs[1] - 3.0).abs() < 1e-5);
        }
        other => panic!("Expected Expr::Matrix for eigenvalues, got: {:?}", other),
    }

    // Inversion with bindings
    let inv_op = matrix_inv(mat);
    let inv_res = compute(&inv_op, &config);

    eprintln!("Bound variables inversion result: {:?}", inv_res);
    match inv_res {
        Expr::Matrix(rows) => {
            assert_eq!(rows.len(), 2);
            // Inverse of [[2, 1], [1, 2]] is 1/3 * [[2, -1], [-1, 2]] = [[2/3, -1/3], [-1/3, 2/3]]
            let get_val = |r: usize, c: usize| match &rows[r][c] {
                Expr::Constant(v) => *v,
                other => panic!("Expected Constant at ({},{}), got {:?}", r, c, other),
            };
            assert!((get_val(0, 0) - (2.0 / 3.0)).abs() < 1e-5);
            assert!((get_val(0, 1) - (-1.0 / 3.0)).abs() < 1e-5);
            assert!((get_val(1, 0) - (-1.0 / 3.0)).abs() < 1e-5);
            assert!((get_val(1, 1) - (2.0 / 3.0)).abs() < 1e-5);
        }
        other => panic!("Expected Expr::Matrix for inversion, got: {:?}", other),
    }
}
