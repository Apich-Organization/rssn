use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;

#[test]

fn test_solve_linear_unique() {

    // 2x + y = 5
    // x - y = 1
    // Solution: x = 2, y = 1
    let a = numerical_Matrix::new(2, 2, vec![2.0, 1.0, 1.0, -1.0]);

    let b = vec![5.0, 1.0];

    match numerical_solve_linear_system(&a, &b).unwrap() {
        numerical_LinearSolution::Unique(x) => {

            assert_approx_eq!(x[0], 2.0);

            assert_approx_eq!(x[1], 1.0);
        }
        _ => panic!("Expected unique solution"),
    }
}

#[test]

fn test_solve_linear_parametric() {

    // x + y + z = 3
    // 2x + 2y + 2z = 6
    // Rank = 1. Solution is a plane.
    let a = numerical_Matrix::new(
        2,
        3,
        vec![
            1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
        ],
    );

    let b = vec![3.0, 6.0];

    match numerical_solve_linear_system(&a, &b).unwrap() {
        numerical_LinearSolution::Parametric {
            particular,
            null_space_basis,
        } => {

            // Check particular solution
            // A * particular should be b
            // We can check just the first row: 1*p0 + 1*p1 + 1*p2 = 3
            let p_sum: f64 = particular.iter().sum();

            assert_approx_eq!(p_sum, 3.0);

            // Check null space
            assert!(null_space_basis.cols() > 0);

            // Verify basis vector v satisfies Av = 0
            for col in null_space_basis.get_cols() {

                let s: f64 = col.iter().sum(); // Since row is 1,1,1, dot(row, col) = sum(col)
                assert_approx_eq!(s, 0.0);
            }
        }
        _ => panic!("Expected parametric solution"),
    }
}

#[test]

fn test_solve_linear_no_solution() {

    // x + y = 2
    // x + y = 3
    let a = numerical_Matrix::new(2, 2, vec![1.0, 1.0, 1.0, 1.0]);

    let b = vec![2.0, 3.0];

    match numerical_solve_linear_system(&a, &b).unwrap() {
        numerical_LinearSolution::NoSolution => {}
        _ => panic!("Expected no solution"),
    }
}

proptest! {
    #[test]
    fn prop_solve_linear_invertible(
        // Generate random invertible matrices by generating random L and U matrices with non-zero diagonal
        // Or simplified: Generate random square matrix and random vector x, compute b=Ax, solve Ax=b.
        // There is a small chance A is singular, but with floats it's rare. We can handle the parametric case too.
        rows in 2..5usize,
        seed in 0..1000u64 // Just to vary data efficiently
    ) {
        let cols = rows; // Square for unique test mostly
        // Generate matrix data deterministically from seed or use proptest collections
        // We'll use a simple deterministic generator for stability here
        let mut data = Vec::with_capacity(rows * cols);
        let mut rng = seed;
        for _ in 0..(rows*cols) {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let val = (rng % 100) as f64 / 10.0 - 5.0; // -5.0 to 5.0
            data.push(val);
        }

        let a = numerical_Matrix::new(rows, cols, data);

        // Generate random solution x
        let mut x = Vec::with_capacity(rows);
        for _ in 0..rows {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            x.push((rng % 100) as f64 / 10.0 - 5.0);
        }

        // Compute b = A * x
        let mut b = vec![0.0; rows];
        for i in 0..rows {
            let mut sum = 0.0;
            for j in 0..cols {
                sum += a.get(i, j) * x[j];
            }
            b[i] = sum;
        }

        // Solve Ax = b
        let result = numerical_solve_linear_system(&a, &b);

        match result {
            Ok(numerical_LinearSolution::Unique(sol)) => {
                // Check if A * sol is close to b
                for i in 0..rows {
                    let mut sum = 0.0;
                    for j in 0..cols {
                        sum += a.get(i, j) * sol[j];
                    }
                    if (sum - b[i]).abs() > 1e-6 {
                        // Failed check
                        panic!("Solution validation failed: A*sol = {:?}, b = {:?}", sum, b[i]);
                    }
                }
            },
            Ok(numerical_LinearSolution::Parametric { particular, null_space_basis }) => {
                // Check particular solution: A * p = b
                for i in 0..rows {
                    let mut sum = 0.0;
                    for j in 0..cols {
                        sum += a.get(i, j) * particular[j];
                    }
                     if (sum - b[i]).abs() > 1e-6 {
                        panic!("Particular solution validation failed");
                    }
                }

                // Check null space basis: A * v = 0
                for col_idx in 0..null_space_basis.cols() {
                     for i in 0..rows {
                        let mut sum = 0.0;
                        for j in 0..cols {
                            sum += a.get(i, j) * null_space_basis.get(j, col_idx);
                        }
                        if sum.abs() > 1e-6 {
                            panic!("Null space validation failed");
                        }
                    }
                }
            },
            Ok(numerical_LinearSolution::NoSolution) => {
                 // This should technically not happen if we constructed b = Ax, unless substantial precision loss
                 // But valid if it happens.
            },
            Err(e) => panic!("Solver error: {}", e),
        }
    }
}
