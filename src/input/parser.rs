#![allow(
    clippy::wildcard_in_or_patterns
)]

use std::sync::Arc;

use nom::IResult;
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::alpha1;
use nom::character::complete::char;
use nom::character::complete::digit1;
use nom::character::complete::i64 as nom_i64;
use nom::character::complete::multispace0;
use nom::character::complete::multispace1;
use nom::combinator::map;
use nom::combinator::map_res;
use nom::combinator::opt;
use nom::combinator::recognize;
use nom::multi::fold_many0;
use nom::multi::separated_list1;
use nom::sequence::delimited;
use nom::sequence::pair;
use nom::sequence::preceded;
use num_bigint::BigInt;
use num_rational::BigRational;

use crate::symbolic::core::Expr;
use crate::symbolic::core::PathType;

/// Checks if a character is a valid identifier character.

fn is_identifier_char(c: char) -> bool {

    c.is_alphanumeric()
        || c == '_'
        || c == '\''
}

pub(crate) fn identifier_name(
    input: &str
) -> IResult<&str, &str> {

    recognize(pair(
        alpha1,
        // nom::multi::many0(nom::character::complete::alphanumeric1)
        nom::bytes::complete::take_while(is_identifier_char),
    ))(input)
}

pub(crate) fn parse_rational_structure(
    input: &str
) -> IResult<&str, ()> {

    let (input, _) = nom_i64(input)?;

    let (input, _) = char('/')(input)?;

    let (input, _) = nom_i64(input)?;

    Ok((input, ()))
}

// Entry point for parsing an expression
/// Parses a string into a symbolic expression.
///
/// # Arguments
/// * `input` - The string to parse.
///
/// # Returns
/// A `Result` containing the remaining input and the parsed `Expr`.
///
/// # Errors
/// Returns a `nom::Err` if the input cannot be parsed as a valid expression.

pub fn parse_expr(
    input: &str
) -> IResult<&str, Expr> {

    expr(input)
}

// expr = comparison_expr
pub(crate) fn expr(
    input: &str
) -> IResult<&str, Expr> {

    comparison_expr(input)
}

// comparison_expr = additive_expr { ("=" | "<" | ">" | "<=" | ">=") additive_expr }
pub(crate) fn comparison_expr(
    input: &str
) -> IResult<&str, Expr> {

    let (input, init) =
        additive_expr(input)?;

    fold_many0(
        pair(
            alt((
                tag("="),
                tag("<="),
                tag(">="),
                tag("<"),
                tag(">"),
            )),
            additive_expr,
        ),
        move || init.clone(),
        |acc, (op, val)| {

            match op {
                | "=" => {
                    Expr::Eq(
                        Arc::new(acc),
                        Arc::new(val),
                    )
                },
                | "<" => {
                    Expr::Lt(
                        Arc::new(acc),
                        Arc::new(val),
                    )
                },
                | ">" => {
                    Expr::Gt(
                        Arc::new(acc),
                        Arc::new(val),
                    )
                },
                | "<=" => {
                    Expr::Le(
                        Arc::new(acc),
                        Arc::new(val),
                    )
                },
                | ">=" => {
                    Expr::Ge(
                        Arc::new(acc),
                        Arc::new(val),
                    )
                },
                | _ => unreachable!(),
            }
        },
    )(input)
}

// additive_expr = term { ("+" | "-") term }
pub(crate) fn additive_expr(
    input: &str
) -> IResult<&str, Expr> {

    let (input, init) = term(input)?;

    fold_many0(
        pair(
            alt((char('+'), char('-'))),
            term,
        ),
        move || init.clone(),
        |acc, (op, val)| {
            if op == '+' {

                Expr::Add(
                    Arc::new(acc),
                    Arc::new(val),
                )
            } else {

                Expr::Sub(
                    Arc::new(acc),
                    Arc::new(val),
                )
            }
        },
    )(input)
}

// term = factor { ("*" | "/") factor }
pub(crate) fn term(
    input: &str
) -> IResult<&str, Expr> {

    let (input, init) = factor(input)?;

    fold_many0(
        pair(
            alt((char('*'), char('/'))),
            factor,
        ),
        move || init.clone(),
        |acc, (op, val)| {
            if op == '*' {

                Expr::Mul(
                    Arc::new(acc),
                    Arc::new(val),
                )
            } else {

                Expr::Div(
                    Arc::new(acc),
                    Arc::new(val),
                )
            }
        },
    )(input)
}

// factor = unary | "(" expr ")"
pub(crate) fn factor(
    input: &str
) -> IResult<&str, Expr> {

    delimited(
        multispace0,
        alt((
            unary,
            parenthesized_expr,
        )),
        multispace0,
    )(input)
}

// unary = ["-"] ["not"] power
pub(crate) fn unary(
    input: &str
) -> IResult<&str, Expr> {

    // println!("in unary staring");
    // println!("{}",input);
    let original_input = input;

    let (input, neg) =
        opt(char('-'))(input)?;

    let (input, not_op) =
        opt(preceded(
            tag("not"),
            multispace1,
        ))(input)?;

    // println!("{}",input);
    ////println!("{}",neg);
    if neg.is_some() {

        if let Ok((_, ())) = nom::combinator::peek(parse_rational_structure)(input) {

            // If it looks like a rational number, we need to parse it differently
            // to handle negative rationals properly
            return parse_rational(original_input);
        }
    }

    let (input, mut expr) =
        power(input)?;

    // println!("in unary start again");
    if neg.is_some() {

        // Special case: -Infinity should be NegativeInfinity, not Neg(Infinity)
        if matches!(
            expr,
            Expr::Infinity
        ) {

            expr =
                Expr::NegativeInfinity;
        } else {

            expr = Expr::Neg(Arc::new(
                expr,
            ));
        }
    }

    if not_op.is_some() {

        expr =
            Expr::Not(Arc::new(expr));
    }

    Ok((input, expr))
}

// power = atom [ "^" unary ] [ "!" ]
pub(crate) fn power(
    input: &str
) -> IResult<&str, Expr> {

    let (input, base) = atom(input)?;

    let (input, power_expr) = opt(
        preceded(char('^'), unary),
    )(
        input
    )?;

    let (input, factorial_op) =
        opt(char('!'))(input)?;

    let mut result = base;

    if let Some(p) = power_expr {

        result = Expr::Power(
            Arc::new(result),
            Arc::new(p),
        );
    }

    if factorial_op.is_some() {

        result = Expr::Factorial(
            Arc::new(result),
        );
    }

    Ok((input, result))
}

pub(crate) fn parse_bigint(
    input: &str
) -> IResult<&str, Expr> {

    map(nom_i64, |n| {

        Expr::BigInt(BigInt::from(n))
    })(input)
}

// pub(crate) fn parse_rational(input: &str) -> IResult<&str, Expr> {
//     //println!("parse_rational");
//     map(pair(nom_i64, preceded(char('/'), nom_i64)), |(num, den)| {
//         Expr::Rational(BigRational::new(BigInt::from(num), BigInt::from(den)))
//     })(input)
// }
// ----------------------------------------------------

pub(crate) fn parse_rational(
    input: &str
) -> IResult<&str, Expr> {

    // println!("parse_rational start");
    let (input, sign) =
        opt(char('-'))(input)?;

    // println!("{}", input);
    let (input, numerator) =
        nom_i64(input)?;

    // println!("{}", numerator);
    let (input, denominator) =
        preceded(char('/'), nom_i64)(
            input,
        )?;

    // println!("{}", denominator);
    let final_numerator =
        if sign.is_some() {

            // println!("matched -");
            // println!("{}", -numerator);
            BigInt::from(-numerator)
        } else {

            // println!("not matched -");
            BigInt::from(numerator)
        };

    // println!("{}", final_numerator);
    // println!("parse_rational end");

    if denominator <= 0 {

        return nom::combinator::fail(
            input,
        );
    }

    // println!("input:");
    // println!("{}", input);
    // println!("output:");
    // println!("{}", Expr::Rational(BigRational::new(final_numerator.clone(), BigInt::from(denominator))));
    Ok((
        input,
        Expr::Rational(
            BigRational::new(
                final_numerator,
                BigInt::from(
                    denominator,
                ),
            ),
        ),
        // Expr::Rational(BigRational::new(BigInt::from(-3), BigInt::from(4))),
    ))
}

// ----------------------------------------------------

pub(crate) fn parse_boolean(
    input: &str
) -> IResult<&str, Expr> {

    alt((
        map(tag("true"), |_| {

            Expr::Boolean(true)
        }),
        map(tag("false"), |_| {

            Expr::Boolean(false)
        }),
    ))(input)
}

pub(crate) fn parse_infinity(
    input: &str
) -> IResult<&str, Expr> {

    map(
        tag("Infinity"),
        |_| Expr::Infinity,
    )(input)
}

pub(crate) fn parse_negative_infinity(
    input: &str
) -> IResult<&str, Expr> {

    map(
        tag("-Infinity"),
        |_| Expr::NegativeInfinity,
    )(input)
}

pub(crate) fn parse_string_literal(
    input: &str
) -> IResult<&str, Expr> {

    map(
        delimited(
            char('"'),
            nom::bytes::complete::take_while(|c| c != '"'),
            char('"'),
        ),
        |s : &str| Expr::Variable(s.to_string()),
    )(input)
}

pub(crate) fn parse_numeric_literals(
    input: &str
) -> IResult<&str, Expr> {

    // println!("in numeric");
    alt((
        parse_rational,
        parse_number,
        parse_bigint,
    ))(input)
    // alt((parse_bigint, parse_number))(input)
}

pub(crate) fn parse_boolean_and_infinities(
    input: &str
) -> IResult<&str, Expr> {

    alt((
        parse_boolean,
        parse_infinity,
        parse_negative_infinity,
    ))(input)
}

#[allow(clippy::too_many_lines)]

pub(crate) fn parse_function_call(
    input: &str
) -> IResult<&str, Expr> {

    // println!("parse function call");
    // let (input, func_name) = alpha1(input)?;
    let (input, func_name) =
        identifier_name(input)?;

    // println!("Parsed function name: {}", func_name);
    let (input, args) = delimited(
        char('('),
        separated_list1(
            delimited(
                multispace0,
                char(','),
                multispace0,
            ),
            expr,
        ),
        char(')'),
    )(input)?;

    // println!("Remaining input after args: '{}'", input);
    match func_name {
        // Two-argument functions with specific names
        | "log_base" => {
            Ok((
                input,
                Expr::LogBase(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "atan2" => {
            Ok((
                input,
                Expr::Atan2(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "falling_factorial" => {
            Ok((
                input,
                Expr::FallingFactorial(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "rising_factorial" => {
            Ok((
                input,
                Expr::RisingFactorial(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "kronecker_delta" => {
            Ok((
                input,
                Expr::KroneckerDelta(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "matrix_mul" => {
            Ok((
                input,
                Expr::MatrixMul(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "matrix_vec_mul" => {
            Ok((
                input,
                Expr::MatrixVecMul(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "quantity_with_value" => {
            Ok((
                input,
                Expr::QuantityWithValue(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                ),
            ))
        },

        // Single-argument functions
        | "sin" => {
            Ok((
                input,
                Expr::Sin(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "cos" => {
            Ok((
                input,
                Expr::Cos(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "tan" => {
            Ok((
                input,
                Expr::Tan(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "log" => {
            Ok((
                input,
                Expr::Log(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "exp" => {
            Ok((
                input,
                Expr::Exp(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "sqrt" => {
            Ok((
                input,
                Expr::Sqrt(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "abs" => {
            Ok((
                input,
                Expr::Abs(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "sec" => {
            Ok((
                input,
                Expr::Sec(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "csc" => {
            Ok((
                input,
                Expr::Csc(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "cot" => {
            Ok((
                input,
                Expr::Cot(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "asin" => {
            Ok((
                input,
                Expr::ArcSin(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acos" => {
            Ok((
                input,
                Expr::ArcCos(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "atan" => {
            Ok((
                input,
                Expr::ArcTan(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "asec" => {
            Ok((
                input,
                Expr::ArcSec(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acsc" => {
            Ok((
                input,
                Expr::ArcCsc(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acot" => {
            Ok((
                input,
                Expr::ArcCot(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "sinh" => {
            Ok((
                input,
                Expr::Sinh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "cosh" => {
            Ok((
                input,
                Expr::Cosh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "tanh" => {
            Ok((
                input,
                Expr::Tanh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "sech" => {
            Ok((
                input,
                Expr::Sech(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "csch" => {
            Ok((
                input,
                Expr::Csch(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "coth" => {
            Ok((
                input,
                Expr::Coth(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "asinh" => {
            Ok((
                input,
                Expr::ArcSinh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acosh" => {
            Ok((
                input,
                Expr::ArcCosh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "atanh" => {
            Ok((
                input,
                Expr::ArcTanh(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "asech" => {
            Ok((
                input,
                Expr::ArcSech(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acsch" => {
            Ok((
                input,
                Expr::ArcCsch(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "acoth" => {
            Ok((
                input,
                Expr::ArcCoth(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "gamma" => {
            Ok((
                input,
                Expr::Gamma(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "erf" => {
            Ok((
                input,
                Expr::Erf(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "erfc" => {
            Ok((
                input,
                Expr::Erfc(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "erfi" => {
            Ok((
                input,
                Expr::Erfi(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "zeta" => {
            Ok((
                input,
                Expr::Zeta(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "digamma" => {
            Ok((
                input,
                Expr::Digamma(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "floor" => {
            Ok((
                input,
                Expr::Floor(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "is_prime" => {
            Ok((
                input,
                Expr::IsPrime(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "transpose" => {
            Ok((
                input,
                Expr::Transpose(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "inverse" => {
            Ok((
                input,
                Expr::Inverse(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "general_solution" => {
            Ok((
                input,
                Expr::GeneralSolution(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "particular_solution" => {
            Ok((
                input,
                Expr::ParticularSolution(Arc::new(
                    args[0].clone(),
                )),
            ))
        },
        | "boundary" => {
            Ok((
                input,
                Expr::Boundary(Arc::new(
                    args[0].clone(),
                )),
            ))
        },

        // Two-argument functions
        | "binomial" => {
            Ok((
                input,
                Expr::Binomial(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "permutation" => {
            Ok((
                input,
                Expr::Permutation(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "combination" => {
            Ok((
                input,
                Expr::Combination(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "beta" => {
            Ok((
                input,
                Expr::Beta(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "besselj" => {
            Ok((
                input,
                Expr::BesselJ(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "bessely" => {
            Ok((
                input,
                Expr::BesselY(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "legendrep" => {
            Ok((
                input,
                Expr::LegendreP(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "laguerrel" => {
            Ok((
                input,
                Expr::LaguerreL(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "hermiteh" => {
            Ok((
                input,
                Expr::HermiteH(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "xor" => {
            Ok((
                input,
                Expr::Xor(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "implies" => {
            Ok((
                input,
                Expr::Implies(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "equivalent" => {
            Ok((
                input,
                Expr::Equivalent(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "gcd" => {
            Ok((
                input,
                Expr::Gcd(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "mod" => {
            Ok((
                input,
                Expr::Mod(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "complex" => {
            Ok((
                input,
                Expr::Complex(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "apply" => {
            Ok((
                input,
                Expr::Apply(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "max" => {
            Ok((
                input,
                Expr::Max(
                    Arc::new(args[0].clone()),
                    Arc::new(args[1].clone()),
                ),
            ))
        },
        | "parametric_solution" => {
            Ok((
                input,
                Expr::ParametricSolution {
                    x : Arc::new(args[0].clone()),
                    y : Arc::new(args[1].clone()),
                },
            ))
        },

        // N-ary functions
        | "and" => {
            Ok((
                input,
                Expr::And(args),
            ))
        },
        | "or" => {
            Ok((
                input,
                Expr::Or(args),
            ))
        },
        | "union" => {
            Ok((
                input,
                Expr::Union(args),
            ))
        },
        | "polynomial" => {
            Ok((
                input,
                Expr::Polynomial(args),
            ))
        },
        | "vector" => {
            Ok((
                input,
                Expr::Vector(args),
            ))
        },
        | "tuple" => {
            Ok((
                input,
                Expr::Tuple(args),
            ))
        },
        | "system" => {
            Ok((
                input,
                Expr::System(args),
            ))
        },
        | "solutions" => {
            Ok((
                input,
                Expr::Solutions(args),
            ))
        },

        // Functions with special parsing
        | "derivative" => {
            Ok((
                input,
                Expr::Derivative(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                ),
            ))
        },
        | "convergence_analysis" => {
            Ok((
                input,
                Expr::ConvergenceAnalysis(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                ),
            ))
        },
        | "solve" => {
            Ok((
                input,
                Expr::Solve(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                ),
            ))
        },
        | "derivative_n" => {
            Ok((
                input,
                Expr::DerivativeN(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                    Arc::new(args[2].clone()),
                ),
            ))
        },
        | "limit" => {
            Ok((
                input,
                Expr::Limit(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                    Arc::new(args[2].clone()),
                ),
            ))
        },
        | "substitute" => {
            Ok((
                input,
                Expr::Substitute(
                    Arc::new(args[0].clone()),
                    match &args[1] {
                        | Expr::Variable(s) => s.clone(),
                        | _ => return nom::combinator::fail(input),
                    },
                    Arc::new(args[2].clone()),
                ),
            ))
        },
        | "predicate" | _ => {
            Ok((
                input,
                Expr::Predicate {
                    name : func_name.to_string(),
                    args,
                },
            ))
        },
    }
}

pub(crate) fn parse_matrix(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("matrix")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, rows) = delimited(
        char('['),
        separated_list1(
            delimited(
                multispace0,
                char(','),
                multispace0,
            ),
            delimited(
                char('['),
                separated_list1(
                    delimited(
                        multispace0,
                        char(','),
                        multispace0,
                    ),
                    expr,
                ),
                char(']'),
            ),
        ),
        char(']'),
    )(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Matrix(rows),
    ))
}

pub(crate) fn parse_pde(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) = tag("pde")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, equation) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, func_name) =
        alpha1(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, vars_list) =
        delimited(
            char('['),
            separated_list1(
                delimited(
                    multispace0,
                    char(','),
                    multispace0,
                ),
                alpha1,
            ),
            char(']'),
        )(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Pde {
            equation: Arc::new(
                equation,
            ),
            func: func_name.to_string(),
            vars: vars_list
                .iter()
                .map(|s| {

                    (*s).to_string()
                })
                .collect(),
        },
    ))
}

// atom = numeric_literals | boolean_and_infinities | matrix | function_call | variable | constant
pub(crate) fn parse_path_type(
    input: &str
) -> IResult<&str, PathType> {

    alt((
        map(tag("Line"), |_| {

            PathType::Line
        }),
        map(
            tag("Circle"),
            |_| PathType::Circle,
        ),
        map(
            tag("Rectangle"),
            |_| PathType::Rectangle,
        ),
    ))(input)
}

pub(crate) fn parse_path(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("path")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, path_type) =
        parse_path_type(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg1) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg2) = expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Path(
            path_type,
            Arc::new(arg1),
            Arc::new(arg2),
        ),
    ))
}

pub(crate) fn parse_interval(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("interval")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, lower_bound) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, upper_bound) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, incl_lower) =
        parse_boolean(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, incl_upper) =
        parse_boolean(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Interval(
            Arc::new(lower_bound),
            Arc::new(upper_bound),
            match incl_lower {
                | Expr::Boolean(b) => b,
                | _ => return nom::combinator::fail(input),
            },
            match incl_upper {
                | Expr::Boolean(b) => b,
                | _ => return nom::combinator::fail(input),
            },
        ),
    ))
}

pub(crate) fn parse_quantifier(
    input: &str
) -> IResult<&str, Expr> {

    let (input, quantifier_type) =
        alt((
            tag("forall"),
            tag("exists"),
        ))(input)?;

    let (input, _) =
        multispace1(input)?;

    let (input, var_name) =
        alpha1(input)?;

    let (input, _) = delimited(
        multispace0,
        char('.'),
        multispace0,
    )(input)?;

    let (input, body) = expr(input)?;

    match quantifier_type {
        | "forall" => {
            Ok((
                input,
                Expr::ForAll(
                    var_name
                        .to_string(),
                    Arc::new(body),
                ),
            ))
        },
        | "exists" => {
            Ok((
                input,
                Expr::Exists(
                    var_name
                        .to_string(),
                    Arc::new(body),
                ),
            ))
        },
        | _ => unreachable!(),
    }
}

pub(crate) fn parse_domain(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("domain")(input)?;

    let (input, domain_name) =
        delimited(
            char('('),
            alpha1,
            char(')'),
        )(input)?;

    Ok((
        input,
        Expr::Domain(
            domain_name.to_string(),
        ),
    ))
}

pub(crate) fn parse_ode(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) = tag("ode")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, equation) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, func_name) =
        alpha1(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, var_name) =
        alpha1(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Ode {
            equation: Arc::new(
                equation,
            ),
            func: func_name.to_string(),
            var: var_name.to_string(),
        },
    ))
}

pub(crate) fn parse_sum(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) = tag("sum")(input)?;

    // println!("step1");
    // println!("{}", input);
    let (input, _) = char('(')(input)?;

    // println!("step2");
    // println!("{}", input);
    let (input, body) = expr(input)?;

    // println!("step3");
    // println!("{}", input);
    // println!("{}", body);
    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    // println!("step4");
    // println!("{}", input);
    let (input, var) = expr(input)?;

    // println!("step5");
    // println!("{}", input);
    // println!("{}", input);
    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, from) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, to) = expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Sum {
            body: Arc::new(body),
            var: Arc::new(var),
            from: Arc::new(from),
            to: Arc::new(to),
        },
    ))
}

pub(crate) fn parse_integral(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("integral")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, integrand) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, var) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, lower_bound) =
        expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, upper_bound) =
        expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Integral {
            integrand: Arc::new(
                integrand,
            ),
            var: Arc::new(var),
            lower_bound: Arc::new(
                lower_bound,
            ),
            upper_bound: Arc::new(
                upper_bound,
            ),
        },
    ))
}

pub(crate) fn parse_series_like_function(
    input: &str
) -> IResult<&str, Expr> {

    let (input, func_name) =
        alt((
            tag("series"),
            tag("summation"),
            tag("product"),
        ))(input)?;

    let (input, _) = char('(')(input)?;

    let (input, arg1) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, var_name_expr) =
        expr(input)?;

    let var_name = match var_name_expr {
        | Expr::Variable(s) => s,
        | _ => return nom::combinator::fail(input),
    };

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg3) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg4) = expr(input)?;

    let (input, _) = char(')')(input)?;

    match func_name {
        | "series" => {
            Ok((
                input,
                Expr::Series(
                    Arc::new(arg1),
                    var_name,
                    Arc::new(arg3),
                    Arc::new(arg4),
                ),
            ))
        },
        | "summation" => {
            Ok((
                input,
                Expr::Summation(
                    Arc::new(arg1),
                    var_name,
                    Arc::new(arg3),
                    Arc::new(arg4),
                ),
            ))
        },
        | "product" => {
            Ok((
                input,
                Expr::Product(
                    Arc::new(arg1),
                    var_name,
                    Arc::new(arg3),
                    Arc::new(arg4),
                ),
            ))
        },
        | _ => unreachable!(),
    }
}

pub(crate) fn parse_asymptotic_expansion(
    input: &str
) -> IResult<&str, Expr> {

    // println!("asymptotic_expansion started");
    let (input, _) = tag(
        "asymptotic_expansion",
    )(input)?;

    // println!("step2");
    // println!("{}", input);
    let (input, _) = char('(')(input)?;

    // println!("step3");
    // println!("{}", input);
    let (input, arg1) = expr(input)?;

    // println!("step4");
    // println!("{}", arg1);
    // println!("{}", input);
    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    // println!("step5");
    let (input, var_name_expr) =
        expr(input)?;

    // println!("Nom started");
    let var_name = match var_name_expr {
        | Expr::Variable(s) => s,
        | _ => return nom::combinator::fail(input),
    };

    // println!("Nom end");
    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg3) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg4) = expr(input)?;

    let (input, _) = char(')')(input)?;

    // println!("{}", input);
    // println!("{}", arg1);
    // println!("{}", arg3);
    // println!("{}", arg4);
    Ok((
        input,
        Expr::AsymptoticExpansion(
            Arc::new(arg1),
            var_name,
            Arc::new(arg3),
            Arc::new(arg4),
        ),
    ))
}

pub(crate) fn parse_fredholm(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("fredholm")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, arg1) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg2) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg3) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg4) = expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Fredholm(
            Arc::new(arg1),
            Arc::new(arg2),
            Arc::new(arg3),
            Arc::new(arg4),
        ),
    ))
}

pub(crate) fn parse_volterra(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("volterra")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, arg1) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg2) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg3) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, arg4) = expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::Volterra(
            Arc::new(arg1),
            Arc::new(arg2),
            Arc::new(arg3),
            Arc::new(arg4),
        ),
    ))
}

pub(crate) fn parse_parametric_solution(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) = tag(
        "parametric_solution",
    )(input)?;

    let (input, _) = char('(')(input)?;

    let (input, x_expr) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, y_expr) = expr(input)?;

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::ParametricSolution {
            x: Arc::new(x_expr),
            y: Arc::new(y_expr),
        },
    ))
}

pub(crate) fn parse_root_of(
    input: &str
) -> IResult<&str, Expr> {

    let (input, _) =
        tag("root_of")(input)?;

    let (input, _) = char('(')(input)?;

    let (input, poly) = expr(input)?;

    let (input, _) = delimited(
        multispace0,
        char(','),
        multispace0,
    )(input)?;

    let (input, index_expr) =
        expr(input)?;

    let index = match index_expr {
        | Expr::Constant(c) => (c as i64).try_into().unwrap_or(0),
        | _ => return nom::combinator::fail(input),
    };

    let (input, _) = char(')')(input)?;

    Ok((
        input,
        Expr::RootOf {
            poly: Arc::new(poly),
            index,
        },
    ))
}

pub(crate) fn atom(
    input: &str
) -> IResult<&str, Expr> {

    alt((
        parse_numeric_literals,
        parse_series_like_function,
        parse_asymptotic_expansion,
        parse_matrix,
        parse_pde,
        parse_path,
        parse_interval,
        parse_quantifier,
        parse_domain,
        parse_ode,
        parse_sum,
        parse_integral,
        parse_fredholm,
        parse_volterra,
        parse_parametric_solution,
        parse_root_of,
        parse_function_call, /* New unified function parser */
        parse_boolean_and_infinities,
        parse_constant,
        alt((
            parenthesized_expr,
            parse_string_literal,
            parse_variable,
        )),
    ))(input)
}

pub(crate) fn parse_float(
    input: &str
) -> IResult<&str, f64> {

    map_res(
        recognize(pair(
            opt(char('-')),
            alt((
                // Leading decimal point (e.g., ".45")
                recognize(pair(
                    char('.'),
                    digit1,
                )),
                // Digits with decimal point and more digits (e.g., "123.45")
                recognize(pair(
                    digit1,
                    pair(
                        char('.'),
                        digit1,
                    ),
                )),
                // Just digits (e.g., "123")
                digit1,
            )),
        )),
        |s: &str| s.parse::<f64>(),
    )(input)
}

// Parses a floating-point number
pub(crate) fn parse_number(
    input: &str
) -> IResult<&str, Expr> {

    map(
        parse_float,
        Expr::Constant,
    )(input)
}

// Parses a mathematical constant
pub(crate) fn parse_constant(
    input: &str
) -> IResult<&str, Expr> {

    alt((
        map(tag("Pi"), |_| {

            Expr::Pi
        }),
        map(tag("E"), |_| {

            Expr::E
        }),
        map(
            tag("InfiniteSolutions"),
            |_| Expr::InfiniteSolutions,
        ),
        map(
            tag("NoSolution"),
            |_| Expr::NoSolution,
        ),
    ))(input)
}

// Parses a variable (a sequence of alphabetic characters)
pub(crate) fn parse_variable(
    input: &str
) -> IResult<&str, Expr> {

    map(
        identifier_name,
        |s: &str| {

            // If the identifier contains a quote, treat it as a Predicate (e.g., y'', y')
            if s.contains('\'') {

                Expr::Predicate {
                    name: s.to_string(),
                    args: vec![],
                }
            } else {

                Expr::Variable(
                    s.to_string(),
                )
            }
        },
    )(input)
}

// Parses an expression enclosed in parentheses
pub(crate) fn parenthesized_expr(
    input: &str
) -> IResult<&str, Expr> {

    delimited(
        char('('),
        expr,
        char(')'),
    )(input)
}

#[cfg(test)]

mod tests {

    use std::sync::Arc;

    use super::*;
    use crate::symbolic::core::Expr;

    #[test]

    fn test_parse_number() {

        assert_eq!(
            parse_expr("123.45"),
            Ok((
                "",
                Expr::Constant(123.45)
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(2.0))
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
                    Arc::new(Expr::Variable(
                        "y".to_string()
                    )),
                    Arc::new(Expr::Constant(3.14))
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
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "z".to_string()
                    )),
                    Arc::new(Expr::Constant(2.5))
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
                            Arc::new(Expr::Constant(1.0))
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
            println!("{}", std::any::type_name::<T>());
        }
        let expected = Expr::Neg(
            Arc::new(Expr::BigInt(
                BigInt::from(456),
            )),
        );
        let aa = parse_expr("123");
        print_type_of(&aa);
        print_type_of(&aa.unwrap().1);
        let bb = Expr::BigInt(BigInt::from(123));
        print_type_of(&bb);
        print_type_of(&Expr::BigInt(BigInt::from(123)));
        assert_eq!(
            parse_expr("123").expect("Parse Expr failed.").1,
            Expr::BigInt(
                BigInt::from(123)
            )
        );

        assert_eq!(
            parse_expr("-456").expect("Parse Expr failed.").1,
            expected
        );
    }

    #[test]

    fn test_parse_rational() {

        let expected = Expr::Neg(
            Arc::new(Expr::Div(
                Arc::new(Expr::BigInt(
                    BigInt::from(3),
                )),
                Arc::new(Expr::BigInt(
                    BigInt::from(4),
                )),
            )),
        );

        assert_eq!(
            parse_expr("1/2"),
            Ok((
                "",
                Expr::Rational(
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
        //         Expr::Rational(BigRational::new(BigInt::from(-3), BigInt::from(4)))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(5.0))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(5.0))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(5.0))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(5.0))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(5.0))
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
                    Arc::new(Expr::Constant(1.0))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Constant(2.0))
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
                    Expr::Variable("x".to_string())
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
                    Expr::Variable("x".to_string())
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
                    Expr::Variable("x".to_string())
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
                    Expr::Variable("x".to_string())
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
                    Expr::Variable("x".to_string())
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
                    Expr::Variable("x".to_string())
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
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    ))
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
                    Expr::Variable("x".to_string())
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    ))
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
                    Expr::Constant(
                        3.14
                    )
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
                Expr::IsPrime(
                    Arc::new(
                        Expr::Constant(
                            7.0
                        )
                    )
                )
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
                        Expr::Constant(
                            12.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
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
                        Expr::Constant(
                            10.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
                            3.0
                        )
                    )
                )
            ))
        );
    }

    #[test]

    fn test_parse_vector() {

        assert_eq!(
            parse_expr(
                "vector(1, 2, 3)"
            ),
            Ok((
                "",
                Expr::Vector(vec![
                    Expr::Constant(1.0),
                    Expr::Constant(2.0),
                    Expr::Constant(3.0),
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
                        Expr::Constant(
                            1.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
                            2.0
                        )
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
                    Expr::Variable("A".to_string())
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
                    Expr::Variable("A".to_string())
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
                    Arc::new(Expr::Variable(
                        "A".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "B".to_string()
                    ))
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
                    Arc::new(Expr::Constant(2.0))
                )
            ))
        );
    }

    #[test]

    fn test_parse_volume_integral() {

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
                    Arc::new(Expr::Constant(0.0)),
                    Arc::new(Expr::Constant(3.0)),
                )
            ))
        );
    }

    #[test]

    fn test_parse_summation() {

        assert_eq!(
            parse_expr("summation(i, i, 1, N)"),
            Ok((
                "",
                Expr::Summation(
                    Arc::new(Expr::Variable(
                        "i".to_string()
                    )),
                    "i".to_string(),
                    Arc::new(Expr::Constant(1.0)),
                    Arc::new(Expr::Variable(
                        "N".to_string()
                    )),
                )
            ))
        );
    }

    #[test]

    fn test_parse_convergence_analysis()
    {

        assert_eq!(
            parse_expr("convergence_analysis(sum(1/n, n, 1, inf), n)"),
            Ok((
                "",
                Expr::ConvergenceAnalysis(
                    Arc::new(Expr::Sum {
                        body : Arc::new(Expr::Div(
                            Arc::new(Expr::Constant(1.0)),
                            Arc::new(Expr::Variable(
                                "n".to_string()
                            )),
                        )),
                        var : Arc::new(Expr::Variable(
                            "n".to_string()
                        )),
                        from : Arc::new(Expr::Constant(1.0)),
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

    fn test_parse_particular_solution()
    {

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
                    Arc::new(Expr::Variable(
                        "f".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "y".to_string()
                    ))
                )
            ))
        );
    }

    #[test]

    fn test_parse_quantity_with_value()
    {

        assert_eq!(
            parse_expr(
                "quantity_with_value(10, \"m\")"
            ),
            Ok((
                "",
                Expr::QuantityWithValue(
                    Arc::new(
                        Expr::Constant(
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
            parse_expr(
                "tuple(1, x, true)"
            ),
            Ok((
                "",
                Expr::Tuple(vec![
                    Expr::Constant(1.0),
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
            parse_expr(
                "system(eq1, eq2)"
            ),
            Ok((
                "",
                Expr::System(vec![
                    Expr::Variable(
                        "eq1"
                            .to_string(
                            )
                    ),
                    Expr::Variable(
                        "eq2"
                            .to_string(
                            )
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
                        "sol1"
                            .to_string(
                            )
                    ),
                    Expr::Variable(
                        "sol2"
                            .to_string(
                            )
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
                    Expr::Variable("D".to_string())
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
                                Arc::new(Expr::Constant(2.0)),
                            )),
                            Arc::new(Expr::Constant(4.0)),
                        )),
                        Arc::new(Expr::Constant(0.0)),
                    )),
                    "x".to_string(),
                )
            ))
        );
    }

    #[test]

    fn test_parse_parametric_solution()
    {

        assert_eq!(
            parse_expr("parametric_solution(t^2, t)"),
            Ok((
                "",
                Expr::ParametricSolution {
                    x : Arc::new(Expr::Power(
                        Arc::new(Expr::Variable(
                            "t".to_string()
                        )),
                        Arc::new(Expr::Constant(2.0)),
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
            parse_expr(
                "InfiniteSolutions"
            ),
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
                            Arc::new(Expr::Constant(2.0)),
                        )),
                        Arc::new(Expr::Constant(1.0)),
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
                        Arc::new(Expr::Constant(2.0)),
                    )),
                    "x".to_string(),
                    Arc::new(Expr::Constant(2.0)),
                )
            ))
        );
    }

    #[test]

    fn test_parse_path() {

        assert_eq!(
            parse_expr(
                "path(Line, 0, 1)"
            ),
            Ok((
                "",
                Expr::Path(
                    PathType::Line,
                    Arc::new(
                        Expr::Constant(
                            0.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
                            1.0
                        )
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

        assert_eq!(
            parse_expr("pde(u_xx + u_yy = 0, u, [x, y])"),
            Ok((
                "",
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
                        Arc::new(Expr::Constant(0.0)),
                    )),
                    func : "u".to_string(),
                    vars : vec![
                        "x".to_string(),
                        "y".to_string()
                    ],
                }
            ))
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
                        Arc::new(Expr::Constant(0.0)),
                    )),
                    func : "y".to_string(),
                    var : "x".to_string(),
                }
            ))
        );
    }

    #[test]

    fn test_parse_asymptotic_expansion()
    {

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
                    Arc::new(Expr::Constant(0.0)),
                    Arc::new(Expr::Constant(3.0)),
                )
            ))
        );
    }

    #[test]

    fn test_parse_product() {

        assert_eq!(
            parse_expr("product(i, i, 1, N)"),
            Ok((
                "",
                Expr::Product(
                    Arc::new(Expr::Variable(
                        "i".to_string()
                    )),
                    "i".to_string(),
                    Arc::new(Expr::Constant(1.0)),
                    Arc::new(Expr::Variable(
                        "N".to_string()
                    )),
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
                        Arc::new(Expr::Constant(2.0)),
                    )),
                    var : Arc::new(Expr::Variable(
                        "i".to_string()
                    )),
                    from : Arc::new(Expr::Constant(1.0)),
                    to : Arc::new(Expr::Constant(10.0)),
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
                        Arc::new(Expr::Constant(2.0)),
                    )),
                    var : Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    lower_bound : Arc::new(Expr::Constant(0.0)),
                    upper_bound : Arc::new(Expr::Constant(1.0)),
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
                    Arc::new(Expr::Constant(0.0))
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
                        Arc::new(Expr::Constant(2.0))
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
                "matrix([[1, 2], [3, \
                 4]])"
            ),
            Ok((
                "",
                Expr::Matrix(vec![
                    vec![
                        Expr::Constant(
                            1.0
                        ),
                        Expr::Constant(
                            2.0
                        )
                    ],
                    vec![
                        Expr::Constant(
                            3.0
                        ),
                        Expr::Constant(
                            4.0
                        )
                    ],
                ])
            ))
        );
    }

    #[test]

    fn test_parse_matrix_vec_mul() {

        assert_eq!(
            parse_expr("matrix_vec_mul(A, v)"),
            Ok((
                "",
                Expr::MatrixVecMul(
                    Arc::new(Expr::Variable(
                        "A".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "v".to_string()
                    ))
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
                Expr::Constant(1.0),
                Expr::Constant(2.0),
                Expr::Constant(3.0),
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
                    Expr::Constant(1.0),
                    Expr::Constant(2.0),
                    Expr::Constant(3.0),
                ])
            ))
        );
    }

    #[test]

    fn test_parse_polynomial_unwrapped()
    {

        let expected_tuple = (
            "",
            Expr::Polynomial(vec![
                Expr::Constant(1.0),
                Expr::Constant(2.0),
                Expr::Constant(3.0),
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

        let static_string : &'static str = "hello";

        let local_string: &str =
            "hello";

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
                        Expr::Constant(
                            0.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
                            1.0
                        )
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
            parse_expr(
                "union(A, B, C)"
            ),
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
            parse_expr("is_prime(x)"),
            Ok((
                "",
                Expr::Predicate {
                    name : "is_prime".to_string(),
                    args : vec![Expr::Variable(
                        "x".to_string()
                    )],
                }
            ))
        );

        assert_eq!(
            parse_expr("has_property(y, z)"),
            Ok((
                "",
                Expr::Predicate {
                    name : "has_property".to_string(),
                    args : vec![
                        Expr::Variable("y".to_string()),
                        Expr::Variable("z".to_string()),
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
                    Arc::new(Expr::Variable(
                        "a".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "b".to_string()
                    ))
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
            parse_expr("kronecker_delta(i, j)"),
            Ok((
                "",
                Expr::KroneckerDelta(
                    Arc::new(Expr::Variable(
                        "i".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "j".to_string()
                    ))
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
                    Expr::Variable("x".to_string())
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
                    Arc::new(Expr::Constant(1.0))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "k".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "k".to_string()
                    ))
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
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "k".to_string()
                    ))
                )
            ))
        );
    }

    #[test]

    fn test_parse_falling_factorial() {

        assert_eq!(
            parse_expr("falling_factorial(x, n)"),
            Ok((
                "",
                Expr::FallingFactorial(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    ))
                )
            ))
        );
    }

    #[test]

    fn test_parse_rising_factorial() {

        assert_eq!(
            parse_expr("rising_factorial(x, n)"),
            Ok((
                "",
                Expr::RisingFactorial(
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "n".to_string()
                    ))
                )
            ))
        );
    }

    #[test]

    fn test_parse_log_base() {

        assert_eq!(
            parse_expr(
                "log_base(2, 8)"
            ),
            Ok((
                "",
                Expr::LogBase(
                    Arc::new(
                        Expr::Constant(
                            2.0
                        )
                    ),
                    Arc::new(
                        Expr::Constant(
                            8.0
                        )
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
                    Arc::new(Expr::Variable(
                        "y".to_string()
                    )),
                    Arc::new(Expr::Variable(
                        "x".to_string()
                    ))
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
                        Arc::new(Expr::Constant(2.0))
                    )
                ))
            ))
        );
    }
}
