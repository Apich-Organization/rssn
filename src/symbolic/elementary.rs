//! # Elementary Functions and Expression Manipulation
//!
//! This module provides constructor functions for creating elementary mathematical expressions
//! (like trigonometric, exponential, and power functions) and tools for manipulating
//! these expressions, such as `expand` using algebraic and trigonometric identities.

use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

/// Creates a sine expression: `sin(expr)`.
#[must_use]

pub fn sin(expr: Expr) -> Expr {

    Expr::new_sin(expr)
}

/// Creates a cosine expression: `cos(expr)`.
#[must_use]

pub fn cos(expr: Expr) -> Expr {

    Expr::new_cos(expr)
}

/// Creates a tangent expression: `tan(expr)`.
#[must_use]

pub fn tan(expr: Expr) -> Expr {

    Expr::new_tan(expr)
}

/// Creates a hyperbolic sine expression: `sinh(expr)`.
#[must_use]

pub fn sinh(expr: Expr) -> Expr {

    Expr::new_sinh(expr)
}

/// Creates a hyperbolic cosine expression: `cosh(expr)`.
#[must_use]

pub fn cosh(expr: Expr) -> Expr {

    Expr::new_cosh(expr)
}

/// Creates a hyperbolic tangent expression: `tanh(expr)`.
#[must_use]

pub fn tanh(expr: Expr) -> Expr {

    Expr::new_tanh(expr)
}

/// Creates a natural logarithm expression: `ln(expr)`.
#[must_use]

pub fn ln(expr: Expr) -> Expr {

    Expr::new_log(expr)
}

/// Creates an exponential expression: `e^(expr)`.
#[must_use]

pub fn exp(expr: Expr) -> Expr {

    Expr::new_exp(expr)
}

/// Creates a square root expression: `sqrt(expr)`.
#[must_use]

pub fn sqrt(expr: Expr) -> Expr {

    Expr::new_sqrt(expr)
}

/// Creates a power expression: `base^exp`.
#[must_use]

pub fn pow(
    base: Expr,
    exp: Expr,
) -> Expr {

    Expr::new_pow(base, exp)
}

/// Returns the symbolic representation of positive infinity.
#[must_use]

pub const fn infinity() -> Expr {

    Expr::Infinity
}

/// Returns the symbolic representation of negative infinity.
#[must_use]

pub const fn negative_infinity() -> Expr
{

    Expr::NegativeInfinity
}

/// Creates a logarithm expression with a specified base: `log_base(expr)`.
#[must_use]

pub fn log_base(
    base: Expr,
    expr: Expr,
) -> Expr {

    Expr::new_log_base(base, expr)
}

/// Creates a cotangent expression: `cot(expr)`.
#[must_use]

pub fn cot(expr: Expr) -> Expr {

    Expr::new_cot(expr)
}

/// Creates a secant expression: `sec(expr)`.
#[must_use]

pub fn sec(expr: Expr) -> Expr {

    Expr::new_sec(expr)
}

/// Creates a cosecant expression: `csc(expr)`.
#[must_use]

pub fn csc(expr: Expr) -> Expr {

    Expr::new_csc(expr)
}

/// Creates an inverse cotangent expression: `acot(expr)`.
#[must_use]

pub fn acot(expr: Expr) -> Expr {

    Expr::new_arccot(expr)
}

/// Creates an inverse secant expression: `asec(expr)`.
#[must_use]

pub fn asec(expr: Expr) -> Expr {

    Expr::new_arcsec(expr)
}

/// Creates an inverse cosecant expression: `acsc(expr)`.
#[must_use]

pub fn acsc(expr: Expr) -> Expr {

    Expr::new_arccsc(expr)
}

/// Creates a hyperbolic cotangent expression: `coth(expr)`.
#[must_use]

pub fn coth(expr: Expr) -> Expr {

    Expr::new_coth(expr)
}

/// Creates a hyperbolic secant expression: `sech(expr)`.
#[must_use]

pub fn sech(expr: Expr) -> Expr {

    Expr::new_sech(expr)
}

/// Creates a hyperbolic cosecant expression: `csch(expr)`.
#[must_use]

pub fn csch(expr: Expr) -> Expr {

    Expr::new_csch(expr)
}

/// Creates an inverse hyperbolic sine expression: `asinh(expr)`.
#[must_use]

pub fn asinh(expr: Expr) -> Expr {

    Expr::new_arcsinh(expr)
}

/// Creates an inverse hyperbolic cosine expression: `acosh(expr)`.
#[must_use]

pub fn acosh(expr: Expr) -> Expr {

    Expr::new_arccosh(expr)
}

/// Creates an inverse hyperbolic tangent expression: `atanh(expr)`.
#[must_use]

pub fn atanh(expr: Expr) -> Expr {

    Expr::new_arctanh(expr)
}

/// Creates an inverse hyperbolic cotangent expression: `acoth(expr)`.
#[must_use]

pub fn acoth(expr: Expr) -> Expr {

    Expr::new_arccoth(expr)
}

/// Creates an inverse hyperbolic secant expression: `asech(expr)`.
#[must_use]

pub fn asech(expr: Expr) -> Expr {

    Expr::new_arcsech(expr)
}

/// Creates an inverse hyperbolic cosecant expression: `acsch(expr)`.
#[must_use]

pub fn acsch(expr: Expr) -> Expr {

    Expr::new_arccsch(expr)
}

/// Creates a 2-argument inverse tangent expression: `atan2(y, x)`.
#[must_use]

pub fn atan2(
    y: Expr,
    x: Expr,
) -> Expr {

    Expr::new_atan2(y, x)
}

/// Returns the symbolic representation of Pi.
#[must_use]

pub const fn pi() -> Expr {

    Expr::Pi
}

/// Returns the symbolic representation of Euler's number (e).
#[must_use]

pub const fn e() -> Expr {

    Expr::E
}

/// Expands a symbolic expression by applying distributive, power, and trigonometric identities.
///
/// This is often the reverse of simplification and can reveal the underlying structure of an expression.
/// The expansion is applied recursively to all parts of the expression.
///
/// # Arguments
/// * `expr` - The expression to expand.
///
/// # Returns
/// A new, expanded `Expr`.
#[must_use]

pub fn expand(expr: Expr) -> Expr {

    simplify(&expand_internal(
        expr,
    ))
}

fn expand_internal(expr: Expr) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            expand_internal(
                node.to_expr()
                    .expect("Expand"),
            )
        },
        | Expr::Complex(re, im) => {

            // Convert Complex(re, im) to Add(re, Mul(im, i)) for expansion
            let i = Expr::Variable(
                "i".to_string(),
            );

            let expanded_re =
                expand_internal(
                    (*re).clone(),
                );

            let expanded_im =
                expand_internal(
                    (*im).clone(),
                );

            Expr::Add(
                Arc::new(expanded_re),
                Arc::new(Expr::Mul(
                    Arc::new(
                        expanded_im,
                    ),
                    Arc::new(i),
                )),
            )
        },
        | Expr::Add(a, b) => {
            Expr::Add(
                Arc::new(
                    expand_internal(
                        (*a).clone(),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        (*b).clone(),
                    ),
                ),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::Sub(
                Arc::new(
                    expand_internal(
                        (*a).clone(),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        (*b).clone(),
                    ),
                ),
            )
        },
        | Expr::Mul(a, b) => {
            expand_mul(
                (*a).clone(),
                (*b).clone(),
            )
        },
        | Expr::Div(a, b) => {
            Expr::Div(
                Arc::new(
                    expand_internal(
                        (*a).clone(),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        (*b).clone(),
                    ),
                ),
            )
        },
        | Expr::Power(b, e) => {
            expand_power(&b, &e)
        },
        | Expr::Log(arg) => {
            expand_log(&arg)
        },
        | Expr::Sin(arg) => {
            expand_sin(&arg)
        },
        | Expr::Cos(arg) => {
            expand_cos(&arg)
        },
        | Expr::Exp(arg) => {
            expand_exp(&arg)
        },
        | _ => expr,
    }
}

/// Expands multiplication over addition: `a*(b+c) -> a*b + a*c`.

pub(crate) fn expand_mul(
    a: Expr,
    b: Expr,
) -> Expr {

    let a_exp = expand_internal(a);

    let b_exp = expand_internal(b);

    match (a_exp, b_exp) {
        | (l, Expr::Add(m, n)) => {
            Expr::Add(
                Arc::new(expand_internal(
                    Expr::Mul(
                        Arc::new(l.clone()),
                        m,
                    ),
                )),
                Arc::new(expand_internal(
                    Expr::Mul(Arc::new(l), n),
                )),
            )
        },
        | (Expr::Add(m, n), r) => {
            Expr::Add(
                Arc::new(expand_internal(
                    Expr::Mul(
                        m,
                        Arc::new(r.clone()),
                    ),
                )),
                Arc::new(expand_internal(
                    Expr::Mul(n, Arc::new(r)),
                )),
            )
        },
        | (l, r) => {
            Expr::Mul(
                Arc::new(l),
                Arc::new(r),
            )
        },
    }
}

/// Expands powers, e.g., `(a*b)^c -> a^c * b^c` and `(a+b)^n -> a^n + ...` (binomial expansion).

pub(crate) fn expand_power(
    base: &Arc<Expr>,
    exp: &Arc<Expr>,
) -> Expr {

    let b_exp = expand_internal(
        base.as_ref()
            .clone(),
    );

    let e_exp = expand_internal(
        exp.as_ref().clone(),
    );

    match (b_exp, e_exp) {
        | (Expr::Mul(f, g), e) => {
            Expr::Mul(
                Arc::new(expand_internal(
                    Expr::Power(
                        f,
                        Arc::new(e.clone()),
                    ),
                )),
                Arc::new(expand_internal(
                    Expr::Power(g, Arc::new(e)),
                )),
            )
        },
        | (Expr::Add(a, b), Expr::BigInt(n)) => {
            if let Some(n_usize) = n.to_usize() {

                expand_binomial(a, b, n_usize)
            } else {

                Expr::Power(
                    Arc::new(Expr::Add(a, b)),
                    Arc::new(Expr::BigInt(n)),
                )
            }
        },
        | (Expr::Add(a, b), Expr::Constant(c)) if c.fract() == 0.0 && c >= 0.0 => {

            let n_usize = c as usize;

            expand_binomial(a, b, n_usize)
        },
        | (b, e) => {
            Expr::Power(
                Arc::new(b),
                Arc::new(e),
            )
        },
    }
}

fn expand_binomial(
    a: Arc<Expr>,
    b: Arc<Expr>,
    n: usize,
) -> Expr {

    let mut sum =
        Expr::BigInt(BigInt::zero());

    for k in 0..=n {

        let bin_coeff = Expr::BigInt(
            binomial_coefficient(n, k),
        );

        let term1 = Expr::Power(
            a.clone(),
            Arc::new(Expr::BigInt(
                BigInt::from(k),
            )),
        );

        let term2 = Expr::Power(
            b.clone(),
            Arc::new(Expr::BigInt(
                BigInt::from(n - k),
            )),
        );

        let term = Expr::Mul(
            Arc::new(bin_coeff),
            Arc::new(Expr::Mul(
                Arc::new(term1),
                Arc::new(term2),
            )),
        );

        sum = Expr::Add(
            Arc::new(sum),
            Arc::new(expand_internal(
                term,
            )),
        );
    }

    sum
}

/// Expands logarithms using identities like `log(a*b) -> log(a) + log(b)`.

pub(crate) fn expand_log(
    arg: &Arc<Expr>
) -> Expr {

    let arg_exp = expand_internal(
        arg.as_ref().clone(),
    );

    match arg_exp {
        | Expr::Mul(a, b) => {
            Expr::Add(
                Arc::new(
                    expand_internal(
                        Expr::Log(a),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        Expr::Log(b),
                    ),
                ),
            )
        },
        | Expr::Div(a, b) => {
            Expr::Sub(
                Arc::new(
                    expand_internal(
                        Expr::Log(a),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        Expr::Log(b),
                    ),
                ),
            )
        },
        | Expr::Power(b, e) => {
            Expr::Mul(
                e,
                Arc::new(
                    expand_internal(
                        Expr::Log(b),
                    ),
                ),
            )
        },
        | a => Expr::Log(Arc::new(a)),
    }
}

/// Expands `sin` using sum-angle identities, e.g., `sin(a+b)`.

pub(crate) fn expand_sin(
    arg: &Arc<Expr>
) -> Expr {

    let arg_exp = expand_internal(
        arg.as_ref().clone(),
    );

    match arg_exp {
        | Expr::Add(a, b) => {
            Expr::Add(
                Arc::new(Expr::Mul(
                    Arc::new(sin(a
                        .as_ref()
                        .clone())),
                    Arc::new(cos(b
                        .as_ref()
                        .clone())),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(cos(a
                        .as_ref()
                        .clone())),
                    Arc::new(sin(b
                        .as_ref()
                        .clone())),
                )),
            )
        },
        | a => Expr::Sin(Arc::new(a)),
    }
}

/// Expands `cos` using sum-angle identities, e.g., `cos(a+b)`.

pub(crate) fn expand_cos(
    arg: &Arc<Expr>
) -> Expr {

    let arg_exp = expand_internal(
        arg.as_ref().clone(),
    );

    match arg_exp {
        | Expr::Add(a, b) => {
            Expr::Sub(
                Arc::new(Expr::Mul(
                    Arc::new(cos(a
                        .as_ref()
                        .clone())),
                    Arc::new(cos(b
                        .as_ref()
                        .clone())),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(sin(a
                        .as_ref()
                        .clone())),
                    Arc::new(sin(b
                        .as_ref()
                        .clone())),
                )),
            )
        },
        | a => Expr::Cos(Arc::new(a)),
    }
}

/// Expands `exp` using identities like `exp(a+b) -> exp(a) * exp(b)`.

pub(crate) fn expand_exp(
    arg: &Arc<Expr>
) -> Expr {

    let arg_exp = expand_internal(
        arg.as_ref().clone(),
    );

    match arg_exp {
        | Expr::Add(a, b) => {
            Expr::Mul(
                Arc::new(
                    expand_internal(
                        Expr::Exp(a),
                    ),
                ),
                Arc::new(
                    expand_internal(
                        Expr::Exp(b),
                    ),
                ),
            )
        },
        | Expr::Mul(a, b) => {

            let is_i = |e: &Expr| matches!(e, Expr::Variable(name) if name == "i");

            if is_i(&a) {

                // exp(i * b) = cos(b) + i * sin(b)
                let cos_b = Expr::Cos(
                    b.clone(),
                );

                let sin_b =
                    Expr::Sin(b);

                Expr::Add(
                    Arc::new(cos_b),
                    Arc::new(
                        Expr::Mul(
                            a,
                            Arc::new(
                                sin_b,
                            ),
                        ),
                    ),
                )
            } else if is_i(&b) {

                // exp(a * i) = cos(a) + i * sin(a)
                let cos_a = Expr::Cos(
                    a.clone(),
                );

                let sin_a =
                    Expr::Sin(a);

                Expr::Add(
                    Arc::new(cos_a),
                    Arc::new(
                        Expr::Mul(
                            b,
                            Arc::new(
                                sin_a,
                            ),
                        ),
                    ),
                )
            } else {

                Expr::Exp(Arc::new(
                    Expr::Mul(a, b),
                ))
            }
        },
        | a => Expr::Exp(Arc::new(a)),
    }
}

/// Helper to compute binomial coefficients C(n, k) = n! / (k! * (n-k)!).
#[must_use]

pub fn binomial_coefficient(
    n: usize,
    k: usize,
) -> BigInt {

    if k > n {

        return BigInt::zero();
    }

    let mut res = BigInt::one();

    for i in 0..k {

        res = (res * (n - i)) / (i + 1);
    }

    res
}
