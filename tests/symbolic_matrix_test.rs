use num_traits::ToPrimitive;
use rssn::symbolic::core::Expr;
use rssn::symbolic::matrix::*;
use rssn::symbolic::simplify_dag::simplify;

#[test]

fn test_add_matrices() {

    let m1 = Expr::Matrix(vec![
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
        ],
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
        ],
    ]);

    let m2 = Expr::Matrix(vec![
        vec![
            Expr::new_constant(5.0),
            Expr::new_constant(6.0),
        ],
        vec![
            Expr::new_constant(7.0),
            Expr::new_constant(8.0),
        ],
    ]);

    let sum = add_matrices(&m1, &m2);

    if let Expr::Matrix(rows) = sum {

        assert_eq!(rows[0][0], Expr::new_constant(6.0));

        assert_eq!(rows[0][1], Expr::new_constant(8.0));

        assert_eq!(rows[1][0], Expr::new_constant(10.0));

        assert_eq!(rows[1][1], Expr::new_constant(12.0));
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_mul_matrices() {

    let m1 = Expr::Matrix(vec![
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
        ],
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
        ],
    ]);

    let m2 = Expr::Matrix(vec![
        vec![
            Expr::new_constant(2.0),
            Expr::new_constant(0.0),
        ],
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
        ],
    ]);

    // [1 2] * [2 0] = [1*2+2*1 1*0+2*2] = [4 4]
    // [3 4]   [1 2]   [3*2+4*1 3*0+4*2]   [10 8]
    let prod = mul_matrices(&m1, &m2);

    if let Expr::Matrix(rows) = prod {

        assert_eq!(rows[0][0], Expr::new_constant(4.0));

        assert_eq!(rows[0][1], Expr::new_constant(4.0));

        assert_eq!(rows[1][0], Expr::new_constant(10.0));

        assert_eq!(rows[1][1], Expr::new_constant(8.0));
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_transpose_matrix() {

    let m = Expr::Matrix(vec![
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
        ],
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
        ],
    ]);

    let t = transpose_matrix(&m);

    if let Expr::Matrix(rows) = t {

        assert_eq!(rows[0][0], Expr::new_constant(1.0));

        assert_eq!(rows[0][1], Expr::new_constant(3.0));

        assert_eq!(rows[1][0], Expr::new_constant(2.0));

        assert_eq!(rows[1][1], Expr::new_constant(4.0));
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_determinant() {

    let m = Expr::Matrix(vec![
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
        ],
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
        ],
    ]);

    // det = 1*4 - 2*3 = 4 - 6 = -2
    let det = determinant(&m);

    assert_eq!(det, Expr::new_constant(-2.0));
}

#[test]

fn test_inverse_matrix() {

    let m = Expr::Matrix(vec![
        vec![
            Expr::new_constant(4.0),
            Expr::new_constant(7.0),
        ],
        vec![
            Expr::new_constant(2.0),
            Expr::new_constant(6.0),
        ],
    ]);

    // det = 24 - 14 = 10
    // inv = 1/10 * [6 -7] = [0.6 -0.7]
    //              [-2 4]   [-0.2 0.4]
    let inv = inverse_matrix(&m);

    if let Expr::Matrix(rows) = inv {

        assert_eq!(rows[0][0], Expr::new_constant(0.6));

        assert_eq!(rows[0][1], Expr::new_constant(-0.7));

        assert_eq!(rows[1][0], Expr::new_constant(-0.2));

        assert_eq!(rows[1][1], Expr::new_constant(0.4));
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_solve_linear_system() {

    // 4x + 7y = 5
    // 2x + 6y = -2
    let a = Expr::Matrix(vec![
        vec![
            Expr::new_constant(4.0),
            Expr::new_constant(7.0),
        ],
        vec![
            Expr::new_constant(2.0),
            Expr::new_constant(6.0),
        ],
    ]);

    let b = Expr::Matrix(vec![
        vec![Expr::new_constant(
            5.0,
        )],
        vec![Expr::new_constant(
            -2.0,
        )],
    ]);

    let sol = solve_linear_system(&a, &b).unwrap();

    if let Expr::Matrix(rows) = sol {

        // x = 4.4, y = -1.8
        let x_expr = &rows[0][0];

        let y_expr = &rows[1][0];

        let x = x_expr.to_ast().unwrap_or(x_expr.clone());

        let y = y_expr.to_ast().unwrap_or(y_expr.clone());

        println!("x = {:?}", x);

        println!("y = {:?}", y);

        // Check x approx 4.4
        if let Expr::Constant(val) = x {

            assert!((val - 4.4).abs() < 1e-10);
        } else if let Expr::Rational(r) = x {

            let val = r.to_f64().unwrap();

            assert!((val - 4.4).abs() < 1e-10);
        } else {

            panic!("Expected constant for x, got {:?}", x);
        }

        // Check y approx -1.8
        if let Expr::Constant(val) = y {

            assert!((val - (-1.8)).abs() < 1e-10);
        } else if let Expr::Rational(r) = y {

            let val = r.to_f64().unwrap();

            assert!((val - (-1.8)).abs() < 1e-10);
        } else {

            panic!("Expected constant for y, got {:?}", y);
        }
    } else {

        panic!("Expected matrix solution");
    }
}
