use std::sync::Arc;

use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::number_theory::chinese_remainder;
use rssn::symbolic::number_theory::extended_gcd;
use rssn::symbolic::number_theory::is_prime;
use rssn::symbolic::number_theory::solve_diophantine;

#[test]

fn test_solve_linear_diophantine() {

    // 3x + 5y = 1
    // Solution: x = 2 + 5t, y = -1 - 3t (or similar)
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let eq = Expr::Eq(
        Arc::new(Expr::new_add(
            Expr::new_mul(
                Expr::new_constant(3.0),
                x.clone(),
            ),
            Expr::new_mul(
                Expr::new_constant(5.0),
                y.clone(),
            ),
        )),
        Arc::new(Expr::new_constant(
            1.0,
        )),
    );

    let vars = vec!["x", "y"];

    let result =
        solve_diophantine(&eq, &vars);

    assert!(result.is_ok());

    let solutions = result.unwrap();

    assert_eq!(solutions.len(), 2);
}

#[test]

fn test_solve_pell() {

    // x^2 - 2y^2 = 1
    // Fundamental solution: (3, 2) -> 3^2 - 2*2^2 = 9 - 8 = 1
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let eq = Expr::Eq(
        Arc::new(Expr::new_sub(
            Expr::new_pow(
                x.clone(),
                Expr::new_constant(2.0),
            ),
            Expr::new_mul(
                Expr::new_constant(2.0),
                Expr::new_pow(
                    y.clone(),
                    Expr::new_constant(
                        2.0,
                    ),
                ),
            ),
        )),
        Arc::new(Expr::new_constant(
            1.0,
        )),
    );

    let vars = vec!["x", "y"];

    let result =
        solve_diophantine(&eq, &vars);

    assert!(result.is_ok());

    let solutions = result.unwrap();

    assert_eq!(solutions.len(), 2);

    // Check if solution is (3, 2)
    // Note: solve_pell returns BigInts wrapped in Expr
    assert_eq!(
        solutions[0],
        Expr::BigInt(BigInt::from(3))
    );

    assert_eq!(
        solutions[1],
        Expr::BigInt(BigInt::from(2))
    );
}

#[test]

fn test_solve_pythagorean() {

    // x^2 + y^2 = z^2
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let z = Expr::new_variable("z");

    let eq = Expr::Eq(
        Arc::new(Expr::new_add(
            Expr::new_pow(
                x.clone(),
                Expr::new_constant(2.0),
            ),
            Expr::new_pow(
                y.clone(),
                Expr::new_constant(2.0),
            ),
        )),
        Arc::new(Expr::new_pow(
            z.clone(),
            Expr::new_constant(2.0),
        )),
    );

    let vars = vec!["x", "y", "z"];

    let result =
        solve_diophantine(&eq, &vars);

    assert!(result.is_ok());

    let solutions = result.unwrap();

    assert_eq!(solutions.len(), 3);
}

#[test]

fn test_extended_gcd() {

    // gcd(10, 6) = 2
    // 10x + 6y = 2 -> 10(-1) + 6(2) = -10 + 12 = 2
    let a =
        Expr::BigInt(BigInt::from(10));

    let b =
        Expr::BigInt(BigInt::from(6));

    let (g, x, y) =
        extended_gcd(&a, &b);

    assert_eq!(
        g,
        Expr::BigInt(BigInt::from(2))
    );

    assert_eq!(
        x,
        Expr::BigInt(BigInt::from(-1))
    );

    assert_eq!(
        y,
        Expr::BigInt(BigInt::from(2))
    );
}

#[test]

fn test_chinese_remainder() {

    // x = 2 mod 3
    // x = 3 mod 5
    // x = 2 mod 7
    // Solution: 23
    // 23 = 2 mod 3 (21+2)
    // 23 = 3 mod 5 (20+3)
    // 23 = 2 mod 7 (21+2)

    let congruences = vec![
        (
            Expr::BigInt(BigInt::from(
                2,
            )),
            Expr::BigInt(BigInt::from(
                3,
            )),
        ),
        (
            Expr::BigInt(BigInt::from(
                3,
            )),
            Expr::BigInt(BigInt::from(
                5,
            )),
        ),
        (
            Expr::BigInt(BigInt::from(
                2,
            )),
            Expr::BigInt(BigInt::from(
                7,
            )),
        ),
    ];

    let result =
        chinese_remainder(&congruences);

    assert!(result.is_some());

    // Result is modulo 3*5*7 = 105
    // 23 mod 105
    // The function returns Expr::Mod(23, 105)
    let res = result.unwrap();

    let (val, modulus) =
        if let Expr::Dag(node) = res {

            match node
                .to_expr()
                .unwrap()
            {
                | Expr::Mod(v, m) => {
                    (v, m)
                },
                | _ => {

                    panic!(
                        "Expected Mod \
                         expression \
                         inside Dag"
                    );
                },
            }
        } else if let Expr::Mod(v, m) =
            res
        {

            (v, m)
        } else {

            panic!(
                "Expected Mod \
                 expression, got {:?}",
                res
            );
        };

    match &*modulus {
        | Expr::BigInt(i) => {

            println!(
                "DEBUG: modulus is \
                 BigInt({})",
                i
            )
        },
        | Expr::Constant(c) => {

            println!(
                "DEBUG: modulus is \
                 Constant({})",
                c
            )
        },
        | Expr::Dag(_) => {

            println!(
                "DEBUG: modulus is Dag"
            )
        },
        | _ => {

            println!(
                "DEBUG: modulus is \
                 {:?}",
                *modulus
            )
        },
    }

    assert_eq!(
        *val,
        Expr::Constant(23.0)
    );

    assert_eq!(
        *modulus,
        Expr::BigInt(BigInt::from(105))
    );
}

#[test]

fn test_is_prime() {

    assert_eq!(
        is_prime(&Expr::BigInt(
            BigInt::from(2)
        )),
        Expr::Boolean(true)
    );

    assert_eq!(
        is_prime(&Expr::BigInt(
            BigInt::from(3)
        )),
        Expr::Boolean(true)
    );

    assert_eq!(
        is_prime(&Expr::BigInt(
            BigInt::from(4)
        )),
        Expr::Boolean(false)
    );

    assert_eq!(
        is_prime(&Expr::BigInt(
            BigInt::from(17)
        )),
        Expr::Boolean(true)
    );

    assert_eq!(
        is_prime(&Expr::BigInt(
            BigInt::from(100)
        )),
        Expr::Boolean(false)
    );

    // Symbolic
    let x = Expr::new_variable("x");

    assert!(matches!(
        is_prime(&x),
        Expr::IsPrime(_)
    ));
}
