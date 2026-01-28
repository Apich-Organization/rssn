use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::input::parser::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::core::PathType;

#[test]

fn test_parse_number() {

    assert_eq!(
        parse_expr("123.45"),
        Ok((
            "",
            Expr::new_constant(123.45)
        ))
    );
}

#[test]

fn test_parse_variable() {

    assert_eq!(
        parse_expr("x"),
        Ok((
            "",
            Expr::Variable(
                "x".to_string()
            )
        ))
    );
}

#[test]

fn test_parse_addition() {

    assert_eq!(
        parse_expr("x + 2"),
        Ok((
            "",
            Expr::Add(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(2.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_subtraction() {

    assert_eq!(
        parse_expr("y - 3.14"),
        Ok((
            "",
            Expr::Sub(
                Arc::new(
                    Expr::Variable(
                        "y".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(
                        3.14
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_multiplication() {

    assert_eq!(
        parse_expr("a * b"),
        Ok((
            "",
            Expr::Mul(
                Arc::new(
                    Expr::Variable(
                        "a".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "b".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_division() {

    assert_eq!(
        parse_expr("z / 2.5"),
        Ok((
            "",
            Expr::Div(
                Arc::new(
                    Expr::Variable(
                        "z".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(2.5)
                )
            )
        ))
    );
}

#[test]

fn test_parse_parentheses() {

    assert_eq!(
            parse_expr("(a + b) * c"),
            Ok((
                "",
                Expr::Mul(
                    Arc::new(Expr::Add(
                        Arc::new(Expr::Variable(
                            "a".to_string()
                        )),
                        Arc::new(Expr::Variable(
                            "b".to_string()
                        ))
                    )),
                    Arc::new(Expr::Variable(
                        "c".to_string()
                    ))
                )
            ))
        );
}

#[test]

fn test_parse_complex_expression() {

    assert_eq!(
            parse_expr("x + y * (z - 1)"),
            Ok((
                "",
                Expr::Add(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Mul(
                        Arc::new(Expr::Variable(
                            "y".to_string()
                        )),
                        Arc::new(Expr::Sub(
                            Arc::new(Expr::Variable(
                                "z".to_string()
                            )),
                            Arc::new(Expr::new_constant(1.0))
                        ))
                    ))
                )
            ))
        );
}

#[test]

fn test_parse_unary_negation() {

    assert_eq!(
        parse_expr("-x"),
        Ok((
            "",
            Expr::Neg(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_pi() {

    assert_eq!(
        parse_expr("Pi"),
        Ok(("", Expr::Pi))
    );
}

#[test]

fn test_parse_e() {

    assert_eq!(
        parse_expr("E"),
        Ok(("", Expr::E))
    );
}

#[test]

fn test_parse_bigint() {

    fn print_type_of<T>(_: &T) {

        println!(
            "{}",
            std::any::type_name::<T>()
        );
    }

    let expected = Expr::Neg(Arc::new(
        Expr::new_bigint(BigInt::from(456)),
    ));

    let aa = parse_expr("123");

    print_type_of(&aa);

    print_type_of(&aa.unwrap().1);

    let bb =
        Expr::new_bigint(BigInt::from(123));

    print_type_of(&bb);

    print_type_of(&Expr::new_bigint(
        BigInt::from(123),
    ));

    // In your test
    let left = parse_expr("123")
        .unwrap()
        .1;

    let right =
        Expr::new_bigint(BigInt::from(123));

    println!("{:#?}", left);

    println!("{:#?}", right);

    assert_eq!(
        parse_expr("123")
            .expect(
                "Parse Expr failed."
            )
            .1,
        Expr::new_bigint(BigInt::from(123))
    );

    assert_eq!(
        parse_expr("-456")
            .expect(
                "Parse Expr failed."
            )
            .1,
        expected
    );
}

#[test]

fn test_parse_rational() {

    let expected =
        Expr::Neg(Arc::new(Expr::Div(
            Arc::new(Expr::new_bigint(
                BigInt::from(3),
            )),
            Arc::new(Expr::new_bigint(
                BigInt::from(4),
            )),
        )));

    assert_eq!(
        parse_expr("1/2"),
        Ok((
            "",
            Expr::new_rational(
                BigRational::new(
                    BigInt::from(1),
                    BigInt::from(2)
                )
            )
        ))
    );
    // assert_eq!(
    //     parse_expr("-3/4"),
    //     Ok((
    //         "",
    //         Expr::new_rational(BigRational::new(BigInt::from(-3), BigInt::from(4)))
    //         //expected
    //     ))
    // );
}

#[test]

fn test_parse_boolean() {

    assert_eq!(
        parse_expr("true"),
        Ok((
            "",
            Expr::Boolean(true)
        ))
    );

    assert_eq!(
        parse_expr("false"),
        Ok((
            "",
            Expr::Boolean(false)
        ))
    );
}

#[test]

fn test_parse_infinity() {

    assert_eq!(
        parse_expr("Infinity"),
        Ok(("", Expr::Infinity))
    );
}

#[test]

fn test_parse_negative_infinity() {

    assert_eq!(
        parse_expr("-Infinity"),
        Ok((
            "",
            Expr::NegativeInfinity
        ))
    );
}

#[test]

fn test_parse_eq() {

    assert_eq!(
        parse_expr("x = 5"),
        Ok((
            "",
            Expr::Eq(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(5.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_lt() {

    assert_eq!(
        parse_expr("x < 5"),
        Ok((
            "",
            Expr::Lt(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(5.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_gt() {

    assert_eq!(
        parse_expr("x > 5"),
        Ok((
            "",
            Expr::Gt(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(5.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_le() {

    assert_eq!(
        parse_expr("x <= 5"),
        Ok((
            "",
            Expr::Le(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(5.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_ge() {

    assert_eq!(
        parse_expr("x >= 5"),
        Ok((
            "",
            Expr::Ge(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(5.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_sin() {

    assert_eq!(
        parse_expr("sin(x)"),
        Ok((
            "",
            Expr::Sin(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_cos() {

    assert_eq!(
            parse_expr("cos(y+1)"),
            Ok((
                "",
                Expr::Cos(Arc::new(Expr::Add(
                    Arc::new(Expr::Variable(
                        "y".to_string()
                    )),
                    Arc::new(Expr::new_constant(1.0))
                )))
            ))
        );
}

#[test]

fn test_parse_power() {

    assert_eq!(
        parse_expr("x^2"),
        Ok((
            "",
            Expr::Power(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::new_constant(2.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_sec() {

    assert_eq!(
        parse_expr("sec(x)"),
        Ok((
            "",
            Expr::Sec(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_csc() {

    assert_eq!(
        parse_expr("csc(x)"),
        Ok((
            "",
            Expr::Csc(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_cot() {

    assert_eq!(
        parse_expr("cot(x)"),
        Ok((
            "",
            Expr::Cot(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arcsin() {

    assert_eq!(
        parse_expr("asin(x)"),
        Ok((
            "",
            Expr::ArcSin(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccos() {

    assert_eq!(
        parse_expr("acos(x)"),
        Ok((
            "",
            Expr::ArcCos(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arctan() {

    assert_eq!(
        parse_expr("atan(x)"),
        Ok((
            "",
            Expr::ArcTan(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arcsec() {

    assert_eq!(
        parse_expr("asec(x)"),
        Ok((
            "",
            Expr::ArcSec(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccsc() {

    assert_eq!(
        parse_expr("acsc(x)"),
        Ok((
            "",
            Expr::ArcCsc(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccot() {

    assert_eq!(
        parse_expr("acot(x)"),
        Ok((
            "",
            Expr::ArcCot(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_sinh() {

    assert_eq!(
        parse_expr("sinh(x)"),
        Ok((
            "",
            Expr::Sinh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_cosh() {

    assert_eq!(
        parse_expr("cosh(x)"),
        Ok((
            "",
            Expr::Cosh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_tanh() {

    assert_eq!(
        parse_expr("tanh(x)"),
        Ok((
            "",
            Expr::Tanh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_sech() {

    assert_eq!(
        parse_expr("sech(x)"),
        Ok((
            "",
            Expr::Sech(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_csch() {

    assert_eq!(
        parse_expr("csch(x)"),
        Ok((
            "",
            Expr::Csch(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_coth() {

    assert_eq!(
        parse_expr("coth(x)"),
        Ok((
            "",
            Expr::Coth(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arcsinh() {

    assert_eq!(
        parse_expr("asinh(x)"),
        Ok((
            "",
            Expr::ArcSinh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccosh() {

    assert_eq!(
        parse_expr("acosh(x)"),
        Ok((
            "",
            Expr::ArcCosh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arctanh() {

    assert_eq!(
        parse_expr("atanh(x)"),
        Ok((
            "",
            Expr::ArcTanh(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arcsech() {

    assert_eq!(
        parse_expr("asech(x)"),
        Ok((
            "",
            Expr::ArcSech(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccsch() {

    assert_eq!(
        parse_expr("acsch(x)"),
        Ok((
            "",
            Expr::ArcCsch(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_arccoth() {

    assert_eq!(
        parse_expr("acoth(x)"),
        Ok((
            "",
            Expr::ArcCoth(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_gamma() {

    assert_eq!(
        parse_expr("gamma(x)"),
        Ok((
            "",
            Expr::Gamma(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_beta() {

    assert_eq!(
        parse_expr("beta(a, b)"),
        Ok((
            "",
            Expr::Beta(
                Arc::new(
                    Expr::Variable(
                        "a".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "b".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_erf() {

    assert_eq!(
        parse_expr("erf(x)"),
        Ok((
            "",
            Expr::Erf(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_erfc() {

    assert_eq!(
        parse_expr("erfc(x)"),
        Ok((
            "",
            Expr::Erfc(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_erfi() {

    assert_eq!(
        parse_expr("erfi(x)"),
        Ok((
            "",
            Expr::Erfi(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_zeta() {

    assert_eq!(
        parse_expr("zeta(s)"),
        Ok((
            "",
            Expr::Zeta(Arc::new(
                Expr::Variable(
                    "s".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_digamma() {

    assert_eq!(
        parse_expr("digamma(x)"),
        Ok((
            "",
            Expr::Digamma(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_besselj() {

    assert_eq!(
        parse_expr("besselj(n, x)"),
        Ok((
            "",
            Expr::BesselJ(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_bessely() {

    assert_eq!(
        parse_expr("bessely(n, x)"),
        Ok((
            "",
            Expr::BesselY(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_legendrep() {

    assert_eq!(
        parse_expr("legendrep(n, x)"),
        Ok((
            "",
            Expr::LegendreP(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_laguerrel() {

    assert_eq!(
        parse_expr("laguerrel(n, x)"),
        Ok((
            "",
            Expr::LaguerreL(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_hermiteh() {

    assert_eq!(
        parse_expr("hermiteh(n, x)"),
        Ok((
            "",
            Expr::HermiteH(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_xor() {

    assert_eq!(
        parse_expr("xor(a, b)"),
        Ok((
            "",
            Expr::Xor(
                Arc::new(
                    Expr::Variable(
                        "a".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "b".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_implies() {

    assert_eq!(
        parse_expr("implies(a, b)"),
        Ok((
            "",
            Expr::Implies(
                Arc::new(
                    Expr::Variable(
                        "a".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "b".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_and() {

    assert_eq!(
        parse_expr("and(a, b, c)"),
        Ok((
            "",
            Expr::And(vec![
                Expr::Variable(
                    "a".to_string()
                ),
                Expr::Variable(
                    "b".to_string()
                ),
                Expr::Variable(
                    "c".to_string()
                ),
            ])
        ))
    );
}

#[test]

fn test_parse_forall() {

    assert_eq!(
            parse_expr("forall x. P(x)"),
            Ok((
                "",
                Expr::ForAll(
                    "x".to_string(),
                    Arc::new(Expr::Predicate {
                        name : "P".to_string(),
                        args : vec![Expr::Variable(
                            "x".to_string()
                        )],
                    })
                )
            ))
        );
}

#[test]

fn test_parse_floor() {

    assert_eq!(
        parse_expr("floor(3.14)"),
        Ok((
            "",
            Expr::Floor(Arc::new(
                Expr::new_constant(3.14)
            ))
        ))
    );
}

#[test]

fn test_parse_is_prime() {

    assert_eq!(
        parse_expr("is_prime(7)"),
        Ok((
            "",
            Expr::IsPrime(Arc::new(
                Expr::new_constant(7.0)
            ))
        ))
    );
}

#[test]

fn test_parse_gcd() {

    assert_eq!(
        parse_expr("gcd(12, 18)"),
        Ok((
            "",
            Expr::Gcd(
                Arc::new(
                    Expr::new_constant(
                        12.0
                    )
                ),
                Arc::new(
                    Expr::new_constant(
                        18.0
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_mod() {

    assert_eq!(
        parse_expr("mod(10, 3)"),
        Ok((
            "",
            Expr::Mod(
                Arc::new(
                    Expr::new_constant(
                        10.0
                    )
                ),
                Arc::new(
                    Expr::new_constant(3.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_vector() {

    assert_eq!(
        parse_expr("vector(1, 2, 3)"),
        Ok((
            "",
            Expr::Vector(vec![
                Expr::new_constant(1.0),
                Expr::new_constant(2.0),
                Expr::new_constant(3.0),
            ])
        ))
    );
}

#[test]

fn test_parse_complex() {

    assert_eq!(
        parse_expr("complex(1, 2)"),
        Ok((
            "",
            Expr::Complex(
                Arc::new(
                    Expr::new_constant(1.0)
                ),
                Arc::new(
                    Expr::new_constant(2.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_transpose() {

    assert_eq!(
        parse_expr("transpose(A)"),
        Ok((
            "",
            Expr::Transpose(Arc::new(
                Expr::Variable(
                    "A".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_inverse() {

    assert_eq!(
        parse_expr("inverse(A)"),
        Ok((
            "",
            Expr::Inverse(Arc::new(
                Expr::Variable(
                    "A".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_matrix_mul() {

    assert_eq!(
        parse_expr("matrix_mul(A, B)"),
        Ok((
            "",
            Expr::MatrixMul(
                Arc::new(
                    Expr::Variable(
                        "A".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "B".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_derivative_n() {

    assert_eq!(
            parse_expr("derivative_n(f(x), x, 2)"),
            Ok((
                "",
                Expr::DerivativeN(
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "x".to_string()
                        )],
                    }),
                    "x".to_string(),
                    Arc::new(Expr::new_constant(2.0))
                )
            ))
        );
}

fn print_type_of_dir<T>(_: &T) {

    println!(
        "{}",
        std::any::type_name::<T>()
    );
}

#[test]

fn test_parse_volume_integral() {

    fn print_type_of<T>(_: &T) {

        println!(
            "{}",
            std::any::type_name::<T>()
        );
    }

    let aa = Expr::VolumeIntegral {
        scalar_field: Arc::new(
            Expr::Predicate {
                name: "f".to_string(),
                args: vec![
                    Expr::Variable(
                        "x".to_string(),
                    ),
                    Expr::Variable(
                        "y".to_string(),
                    ),
                    Expr::Variable(
                        "z".to_string(),
                    ),
                ],
            },
        ),
        volume: Arc::new(
            Expr::Variable(
                "V".to_string(),
            ),
        ),
    };

    print_type_of(&parse_expr(
        "volume_integral(f(x,y,z), V)",
    ));

    print_type_of(&aa);

    assert_eq!(
            parse_expr("volume_integral(f(x,y,z), V)"),
            Ok((
                "",
                Expr::VolumeIntegral {
                    scalar_field : Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![
                            Expr::Variable("x".to_string()),
                            Expr::Variable("y".to_string()),
                            Expr::Variable("z".to_string()),
                        ],
                    }),
                    volume : Arc::new(Expr::Variable(
                        "V".to_string()
                    )),
                }
            ))
        );
}

#[test]

fn test_parse_series() {

    assert_eq!(
            parse_expr("series(f(x), x, 0, 3)"),
            Ok((
                "",
                Expr::Series(
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "x".to_string()
                        )],
                    }),
                    "x".to_string(),
                    Arc::new(Expr::new_constant(0.0)),
                    Arc::new(Expr::new_constant(3.0)),
                )
            ))
        );
}

#[test]

fn test_parse_summation() {

    assert_eq!(
        parse_expr(
            "summation(i, i, 1, N)"
        ),
        Ok((
            "",
            Expr::Summation(
                Arc::new(
                    Expr::Variable(
                        "i".to_string()
                    )
                ),
                "i".to_string(),
                Arc::new(
                    Expr::new_constant(1.0)
                ),
                Arc::new(
                    Expr::Variable(
                        "N".to_string()
                    )
                ),
            )
        ))
    );
}

#[test]

fn test_parse_convergence_analysis() {

    assert_eq!(
            parse_expr("convergence_analysis(sum(1/n, n, 1, inf), n)"),
            Ok((
                "",
                Expr::ConvergenceAnalysis(
                    Arc::new(Expr::Sum {
                        body : Arc::new(Expr::Div(
                            Arc::new(Expr::new_constant(1.0)),
                            Arc::new(Expr::Variable(
                                "n".to_string()
                            )),
                        )),
                        var : Arc::new(Expr::Variable(
                            "n".to_string()
                        )),
                        from : Arc::new(Expr::new_constant(1.0)),
                        to : Arc::new(Expr::Variable(
                            "inf".to_string()
                        )),
                    }),
                    "n".to_string(),
                )
            ))
        );
}

#[test]

fn test_parse_general_solution() {

    assert_eq!(
            parse_expr("general_solution(C1*cos(x) + C2*sin(x))"),
            Ok((
                "",
                Expr::GeneralSolution(Arc::new(Expr::Add(
                    Arc::new(Expr::Mul(
                        Arc::new(Expr::Variable(
                            "C1".to_string()
                        )),
                        Arc::new(Expr::Cos(Arc::new(
                            Expr::Variable("x".to_string())
                        ))),
                    )),
                    Arc::new(Expr::Mul(
                        Arc::new(Expr::Variable(
                            "C2".to_string()
                        )),
                        Arc::new(Expr::Sin(Arc::new(
                            Expr::Variable("x".to_string())
                        ))),
                    )),
                )))
            ))
        );
}

#[test]

fn test_parse_particular_solution() {

    assert_eq!(
            parse_expr("particular_solution(sin(x))"),
            Ok((
                "",
                Expr::ParticularSolution(Arc::new(Expr::Sin(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
                )))
            ))
        );
}

#[test]

fn test_parse_fredholm() {

    assert_eq!(
            parse_expr("fredholm(K(x,t), f(t), a, b)"),
            Ok((
                "",
                Expr::Fredholm(
                    Arc::new(Expr::Predicate {
                        name : "K".to_string(),
                        args : vec![
                            Expr::Variable("x".to_string()),
                            Expr::Variable("t".to_string())
                        ],
                    }),
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "t".to_string()
                        )],
                    }),
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    )),
                )
            ))
        );
}

#[test]

fn test_parse_apply() {

    assert_eq!(
        parse_expr("apply(f, x)"),
        Ok((
            "",
            Expr::Apply(
                Arc::new(
                    Expr::Variable(
                        "f".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_max() {

    assert_eq!(
        parse_expr("max(x, y)"),
        Ok((
            "",
            Expr::Max(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "y".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_quantity_with_value() {

    assert_eq!(
        parse_expr(
            "quantity_with_value(10, \
             \"m\")"
        ),
        Ok((
            "",
            Expr::QuantityWithValue(
                Arc::new(
                    Expr::new_constant(
                        10.0
                    )
                ),
                "m".to_string()
            )
        ))
    );
}

#[test]

fn test_parse_tuple() {

    assert_eq!(
        parse_expr("tuple(1, x, true)"),
        Ok((
            "",
            Expr::Tuple(vec![
                Expr::new_constant(1.0),
                Expr::Variable(
                    "x".to_string()
                ),
                Expr::Boolean(true),
            ])
        ))
    );
}

#[test]

fn test_parse_system() {

    assert_eq!(
        parse_expr("system(eq1, eq2)"),
        Ok((
            "",
            Expr::System(vec![
                Expr::Variable(
                    "eq1".to_string()
                ),
                Expr::Variable(
                    "eq2".to_string()
                ),
            ])
        ))
    );
}

#[test]

fn test_parse_solutions() {

    assert_eq!(
        parse_expr(
            "solutions(sol1, sol2)"
        ),
        Ok((
            "",
            Expr::Solutions(vec![
                Expr::Variable(
                    "sol1".to_string()
                ),
                Expr::Variable(
                    "sol2".to_string()
                ),
            ])
        ))
    );
}

#[test]

fn test_parse_boundary() {

    assert_eq!(
        parse_expr("boundary(D)"),
        Ok((
            "",
            Expr::Boundary(Arc::new(
                Expr::Variable(
                    "D".to_string()
                )
            ))
        ))
    );
}

#[test]

fn test_parse_domain() {

    assert_eq!(
        parse_expr("domain(R)"),
        Ok((
            "",
            Expr::Domain(
                "R".to_string()
            )
        ))
    );
}

#[test]

fn test_parse_solve() {

    assert_eq!(
            parse_expr("solve(x^2 - 4 = 0, x)"),
            Ok((
                "",
                Expr::Solve(
                    Arc::new(Expr::Eq(
                        Arc::new(Expr::Sub(
                            Arc::new(Expr::Power(
                                Arc::new(Expr::Variable(
                                    "x".to_string()
                                )),
                                Arc::new(Expr::new_constant(2.0)),
                            )),
                            Arc::new(Expr::new_constant(4.0)),
                        )),
                        Arc::new(Expr::new_constant(0.0)),
                    )),
                    "x".to_string(),
                )
            ))
        );
}

#[test]

fn test_parse_parametric_solution() {

    assert_eq!(
            parse_expr("parametric_solution(t^2, t)"),
            Ok((
                "",
                Expr::ParametricSolution {
                    x : Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "t".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0)),
                    )),
                    y : Arc::new(Expr::Variable(
                        "t".to_string()
                    )),
                }
            ))
        );
}

#[test]

fn test_parse_infinite_solutions() {

    assert_eq!(
        parse_expr("InfiniteSolutions"),
        Ok((
            "",
            Expr::InfiniteSolutions
        ))
    );
}

#[test]

fn test_parse_no_solution() {

    assert_eq!(
        parse_expr("NoSolution"),
        Ok(("", Expr::NoSolution))
    );
}

#[test]

fn test_parse_root_of() {

    assert_eq!(
            parse_expr("root_of(x^2 - 1, 1)"),
            Ok((
                "",
                Expr::RootOf {
                    poly : Arc::new(Expr::Sub(
                        Arc::new(Expr::Power(
                            Arc::new(Expr::Variable(
                                "x".to_string()
                            )),
                            Arc::new(Expr::new_constant(2.0)),
                        )),
                        Arc::new(Expr::new_constant(1.0)),
                    )),
                    index : 1,
                }
            ))
        );
}

#[test]

fn test_parse_substitute() {

    assert_eq!(
            parse_expr("substitute(x^2, x, 2)"),
            Ok((
                "",
                Expr::Substitute(
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "x".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0)),
                    )),
                    "x".to_string(),
                    Arc::new(Expr::new_constant(2.0)),
                )
            ))
        );
}

#[test]

fn test_parse_path() {

    assert_eq!(
        parse_expr("path(Line, 0, 1)"),
        Ok((
            "",
            Expr::Path(
                PathType::Line,
                Arc::new(
                    Expr::new_constant(0.0)
                ),
                Arc::new(
                    Expr::new_constant(1.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_volterra() {

    assert_eq!(
            parse_expr("volterra(K(x,t), f(t), a, x)"),
            Ok((
                "",
                Expr::Volterra(
                    Arc::new(Expr::Predicate {
                        name : "K".to_string(),
                        args : vec![
                            Expr::Variable("x".to_string()),
                            Expr::Variable("t".to_string())
                        ],
                    }),
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "t".to_string()
                        )],
                    }),
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                )
            ))
        );
}

#[test]

fn test_parse_pde() {

    let result = Expr::Pde {
        equation: Arc::new(Expr::Eq(
            Arc::new(Expr::Add(
                Arc::new(
                    Expr::Predicate {
                        name: "u_xx"
                            .to_string(
                            ),
                        args: vec![],
                    },
                ),
                Arc::new(
                    Expr::Predicate {
                        name: "u_yy"
                            .to_string(
                            ),
                        args: vec![],
                    },
                ),
            )),
            Arc::new(Expr::new_constant(
                0.0,
            )),
        )),
        func: "u".to_string(),
        vars: vec![
            "x".to_string(),
            "y".to_string(),
        ],
    };

    println!(
        "{:?}",
        print_type_of(&result)
    );

    println!(
        "{:?}",
        print_type_of(&parse_expr(
            "pde(u_xx() + u_yy() = 0, \
             u, [x, y])"
        ))
    );

    print_type_of_dir(&result);

    print_type_of_dir(
        &parse_expr(
            "pde(u_xx + u_yy = 0, u, \
             [x, y])",
        )
        .unwrap()
        .1,
    );

    assert_eq!(
            parse_expr("pde(u_xx + u_yy = 0, u, [x, y])").unwrap().1,
                Expr::Pde {
                    equation : Arc::new(Expr::Eq(
                        Arc::new(Expr::Add(
                            Arc::new(Expr::Predicate {
                                name : "u_xx".to_string(),
                                args : vec![],
                            }),
                            Arc::new(Expr::Predicate {
                                name : "u_yy".to_string(),
                                args : vec![],
                            }),
                        )),
                        Arc::new(Expr::new_constant(0.0)),
                    )),
                    func : "u".to_string(),
                    vars : vec![
                        "x".to_string(),
                        "y".to_string()
                    ],
                }
        );
}

#[test]

fn test_parse_ode() {

    assert_eq!(
            parse_expr("ode(y'' + y = 0, y, x)"),
            Ok((
                "",
                Expr::Ode {
                    equation : Arc::new(Expr::Eq(
                        Arc::new(Expr::Add(
                            Arc::new(Expr::Predicate {
                                name : "y''".to_string(),
                                args : vec![],
                            }),
                            Arc::new(Expr::Variable(
                                "y".to_string()
                            )),
                        )),
                        Arc::new(Expr::new_constant(0.0)),
                    )),
                    func : "y".to_string(),
                    var : "x".to_string(),
                }
            ))
        );
}

#[test]

fn test_parse_asymptotic_expansion() {

    assert_eq!(
            parse_expr("asymptotic_expansion(f(x), x, 0, 3)"),
            Ok((
                "",
                Expr::AsymptoticExpansion(
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "x".to_string()
                        )],
                    }),
                    "x".to_string(),
                    Arc::new(Expr::new_constant(0.0)),
                    Arc::new(Expr::new_constant(3.0)),
                )
            ))
        );
}

#[test]

fn test_parse_product() {

    assert_eq!(
        parse_expr(
            "product(i, i, 1, N)"
        ),
        Ok((
            "",
            Expr::Product(
                Arc::new(
                    Expr::Variable(
                        "i".to_string()
                    )
                ),
                "i".to_string(),
                Arc::new(
                    Expr::new_constant(1.0)
                ),
                Arc::new(
                    Expr::Variable(
                        "N".to_string()
                    )
                ),
            )
        ))
    );
}

#[test]

fn test_parse_sum() {

    assert_eq!(
            parse_expr("sum(i^2, i, 1, 10)"),
            Ok((
                "",
                Expr::Sum {
                    body : Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "i".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0)),
                    )),
                    var : Arc::new(Expr::Variable(
                        "i".to_string()
                    )),
                    from : Arc::new(Expr::new_constant(1.0)),
                    to : Arc::new(Expr::new_constant(10.0)),
                }
            ))
        );
}

#[test]

fn test_parse_surface_integral() {

    assert_eq!(
            parse_expr("surface_integral(F(x,y,z), S)"),
            Ok((
                "",
                Expr::SurfaceIntegral {
                    vector_field : Arc::new(Expr::Predicate {
                        name : "F".to_string(),
                        args : vec![
                            Expr::Variable("x".to_string()),
                            Expr::Variable("y".to_string()),
                            Expr::Variable("z".to_string()),
                        ],
                    }),
                    surface : Arc::new(Expr::Variable(
                        "S".to_string()
                    )),
                }
            ))
        );
}

#[test]

fn test_parse_integral() {

    assert_eq!(
            parse_expr("integral(x^2, x, 0, 1)"),
            Ok((
                "",
                Expr::Integral {
                    integrand : Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "x".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0)),
                    )),
                    var : Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    lower_bound : Arc::new(Expr::new_constant(0.0)),
                    upper_bound : Arc::new(Expr::new_constant(1.0)),
                }
            ))
        );
}

#[test]

fn test_parse_limit() {

    assert_eq!(
            parse_expr("limit(f(x), x, 0)"),
            Ok((
                "",
                Expr::Limit(
                    Arc::new(Expr::Predicate {
                        name : "f".to_string(),
                        args : vec![Expr::Variable(
                            "x".to_string()
                        )],
                    }),
                    "x".to_string(),
                    Arc::new(Expr::new_constant(0.0))
                )
            ))
        );
}

#[test]

fn test_parse_derivative() {

    assert_eq!(
            parse_expr("derivative(x^2, x)"),
            Ok((
                "",
                Expr::Derivative(
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "x".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0))
                    )),
                    "x".to_string()
                )
            ))
        );
}

#[test]

fn test_parse_matrix() {

    assert_eq!(
        parse_expr(
            "matrix([[1, 2], [3, 4]])"
        ),
        Ok((
            "",
            Expr::Matrix(vec![
                vec![
                    Expr::new_constant(1.0),
                    Expr::new_constant(2.0)
                ],
                vec![
                    Expr::new_constant(3.0),
                    Expr::new_constant(4.0)
                ],
            ])
        ))
    );
}

#[test]

fn test_parse_matrix_vec_mul() {

    assert_eq!(
        parse_expr(
            "matrix_vec_mul(A, v)"
        ),
        Ok((
            "",
            Expr::MatrixVecMul(
                Arc::new(
                    Expr::Variable(
                        "A".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "v".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_polynomial() {

    print_type_of(&parse_expr(
        "polynomial(1, 2, 3)",
    ));

    let aatest: Result<
        (&str, Expr),
        (),
    > = Ok((
        "",
        Expr::Polynomial(vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
            Expr::new_constant(3.0),
        ]),
    ));

    print_type_of(&aatest);

    assert_eq!(
        parse_expr(
            "polynomial(1, 2, 3)"
        ),
        Ok((
            "",
            Expr::Polynomial(vec![
                Expr::new_constant(1.0),
                Expr::new_constant(2.0),
                Expr::new_constant(3.0),
            ])
        ))
    );
}

#[test]

fn test_parse_polynomial_unwrapped() {

    let expected_tuple = (
        "",
        Expr::Polynomial(vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
            Expr::new_constant(3.0),
        ]),
    );

    let orgvalue = parse_expr(
        "polynomial(1, 2, 3)",
    )
    .unwrap();

    print_type_of(&orgvalue);

    print_type_of(&expected_tuple);

    // println!("first test started");
    assert_eq!(
        orgvalue,
        expected_tuple
    );

    // println!("second test passed");
    assert_eq!(
        parse_expr(
            "polynomial(1, 2, 3)"
        )
        .unwrap(),
        expected_tuple
    );
}

#[test]

fn dag_test() {

    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    assert_eq!(
        Expr::new_add(&a, &b),
        Expr::new_add(&a, &b)
    );
}

#[test]

fn prove_type02() {

    let static_string: &'static str =
        "hello";

    let local_string: &str = "hello";

    let different_string: &str =
        "world";

    assert_eq!(
        static_string,
        local_string
    );

    assert_ne!(
        static_string,
        different_string
    );
}

#[test]

fn prove_type() {

    let u: i32 = 3;

    let i: i32 = 3;

    assert_eq!(
        Ok::<i32, ()>(u),
        Ok(i)
    )
}

use std::any::type_name;

fn print_type_of<T>(_: &T) {
    // println!("Type: {}", type_name::<T>());
}

#[test]

fn test_parse_polynomial02() {

    print_type_of(&parse_expr(
        "polynomial(1, 2, 3)",
    ));

    assert_eq!(
        parse_expr(
            "polynomial(1, 2, 3)"
        ),
        parse_expr(
            "polynomial(1, 2, 3)"
        )
    );
}

#[test]

fn test_parse_interval() {

    print_type_of(&parse_expr(
        "polynomial(1, 2, 3)",
    ));

    assert_eq!(
        parse_expr(
            "interval(0, 1, true, \
             false)"
        ),
        Ok((
            "",
            Expr::Interval(
                Arc::new(
                    Expr::new_constant(0.0)
                ),
                Arc::new(
                    Expr::new_constant(1.0)
                ),
                true,
                false
            )
        ))
    );
}

#[test]

fn test_parse_union() {

    assert_eq!(
        parse_expr("union(A, B, C)"),
        Ok((
            "",
            Expr::Union(vec![
                Expr::Variable(
                    "A".to_string()
                ),
                Expr::Variable(
                    "B".to_string()
                ),
                Expr::Variable(
                    "C".to_string()
                ),
            ])
        ))
    );
}

#[test]

fn test_parse_exists() {

    assert_eq!(
            parse_expr("exists y. Q(y)"),
            Ok((
                "",
                Expr::Exists(
                    "y".to_string(),
                    Arc::new(Expr::Predicate {
                        name : "Q".to_string(),
                        args : vec![Expr::Variable(
                            "y".to_string()
                        )],
                    })
                )
            ))
        );
}

#[test]

fn test_parse_predicate() {

    assert_eq!(
        parse_expr("is_even(x)"),
        Ok((
            "",
            Expr::Predicate {
                name: "is_even"
                    .to_string(),
                args: vec![
                    Expr::Variable(
                        "x".to_string()
                    )
                ],
            }
        ))
    );

    assert_eq!(
        parse_expr(
            "has_property(y, z)"
        ),
        Ok((
            "",
            Expr::Predicate {
                name: "has_property"
                    .to_string(),
                args: vec![
                    Expr::Variable(
                        "y".to_string()
                    ),
                    Expr::Variable(
                        "z".to_string()
                    ),
                ],
            }
        ))
    );
}

#[test]

fn test_parse_or() {

    assert_eq!(
        parse_expr("or(x, y)"),
        Ok((
            "",
            Expr::Or(vec![
                Expr::Variable(
                    "x".to_string()
                ),
                Expr::Variable(
                    "y".to_string()
                ),
            ])
        ))
    );
}

#[test]

fn test_parse_equivalent() {

    assert_eq!(
        parse_expr("equivalent(a, b)"),
        Ok((
            "",
            Expr::Equivalent(
                Arc::new(
                    Expr::Variable(
                        "a".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "b".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_not() {

    assert_eq!(
        parse_expr("not x"),
        Ok((
            "",
            Expr::Not(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );

    assert_eq!(
            parse_expr("not (x = y)"),
            Ok((
                "",
                Expr::Not(Arc::new(Expr::Eq(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "y".to_string()
                    ))
                )))
            ))
        );
}

#[test]

fn test_parse_kronecker_delta() {

    assert_eq!(
        parse_expr(
            "kronecker_delta(i, j)"
        ),
        Ok((
            "",
            Expr::KroneckerDelta(
                Arc::new(
                    Expr::Variable(
                        "i".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "j".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_factorial() {

    assert_eq!(
        parse_expr("x!"),
        Ok((
            "",
            Expr::Factorial(Arc::new(
                Expr::Variable(
                    "x".to_string()
                )
            ))
        ))
    );

    assert_eq!(
            parse_expr("(x+1)!"),
            Ok((
                "",
                Expr::Factorial(Arc::new(Expr::Add(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::new_constant(1.0))
                )))
            ))
        );
}

#[test]

fn test_parse_binomial() {

    assert_eq!(
        parse_expr("binomial(n, k)"),
        Ok((
            "",
            Expr::Binomial(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "k".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_permutation() {

    assert_eq!(
        parse_expr("permutation(n, k)"),
        Ok((
            "",
            Expr::Permutation(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "k".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_combination() {

    assert_eq!(
        parse_expr("combination(n, k)"),
        Ok((
            "",
            Expr::Combination(
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "k".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_falling_factorial() {

    assert_eq!(
        parse_expr(
            "falling_factorial(x, n)"
        ),
        Ok((
            "",
            Expr::FallingFactorial(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_rising_factorial() {

    assert_eq!(
        parse_expr(
            "rising_factorial(x, n)"
        ),
        Ok((
            "",
            Expr::RisingFactorial(
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "n".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_log_base() {

    assert_eq!(
        parse_expr("log_base(2, 8)"),
        Ok((
            "",
            Expr::LogBase(
                Arc::new(
                    Expr::new_constant(2.0)
                ),
                Arc::new(
                    Expr::new_constant(8.0)
                )
            )
        ))
    );
}

#[test]

fn test_parse_atan2() {

    assert_eq!(
        parse_expr("atan2(y, x)"),
        Ok((
            "",
            Expr::Atan2(
                Arc::new(
                    Expr::Variable(
                        "y".to_string()
                    )
                ),
                Arc::new(
                    Expr::Variable(
                        "x".to_string()
                    )
                )
            )
        ))
    );
}

#[test]

fn test_parse_power_and_negation() {

    assert_eq!(
            parse_expr("-x^2"),
            Ok((
                "",
                Expr::Neg(Arc::new(
                    Expr::Power(
                        Arc::new(Expr::Variable(
                            "x".to_string()
                        )),
                        Arc::new(Expr::new_constant(2.0))
                    )
                ))
            ))
        );
}
