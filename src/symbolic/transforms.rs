//! # Symbolic Integral Transforms
//!
//! This module provides functions for performing symbolic integral transforms,
//! including the Fourier, Laplace, and Z-transforms, as well as their inverses.
//! It also includes implementations of key transform properties and theorems,
//! such as the convolution theorem.

use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use crate::symbolic::calculus::definite_integrate;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::path_integrate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::solve::solve;

pub(crate) fn i_complex() -> Expr {

    Expr::new_complex(
        Expr::BigInt(BigInt::zero()),
        Expr::BigInt(BigInt::one()),
    )
}

/// Applies the time-shift property of the Fourier Transform.
///
/// If `F(ω)` is the Fourier Transform of `f(t)`, then the Fourier Transform
/// of `f(t - a)` is `e^(-jωa) * F(ω)`.
///
/// # Arguments
/// * `f_omega` - The Fourier Transform of the original function `f(t)`.
/// * `a` - The time shift amount.
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the Fourier Transform of the time-shifted function.
#[must_use]

pub fn fourier_time_shift(
    f_omega: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_exp(Expr::new_mul(
            Expr::new_mul(
                i_complex(),
                Expr::Variable(
                    out_var.to_string(),
                ),
            ),
            Expr::new_neg(a.clone()),
        )),
        f_omega.clone(),
    ))
}

/// Applies the frequency-shift property of the Fourier Transform.
///
/// If `F(ω)` is the Fourier Transform of `f(t)`, then the Fourier Transform
/// of `e^(jat) * f(t)` is `F(ω - a)`.
///
/// # Arguments
/// * `f_omega` - The Fourier Transform of the original function `f(t)`.
/// * `a` - The frequency shift amount.
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the Fourier Transform of the frequency-shifted function.
#[must_use]

pub fn fourier_frequency_shift(
    f_omega: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::Substitute(
        Arc::new(f_omega.clone()),
        out_var.to_string(),
        Arc::new(Expr::new_sub(
            Expr::Variable(
                out_var.to_string(),
            ),
            a.clone(),
        )),
    ))
}

/// Applies the scaling property of the Fourier Transform.
///
/// If `F(ω)` is the Fourier Transform of `f(t)`, then the Fourier Transform
/// of `f(at)` is `(1/|a|) * F(ω/a)`.
///
/// # Arguments
/// * `f_omega` - The Fourier Transform of the original function `f(t)`.
/// * `a` - The scaling factor.
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the Fourier Transform of the scaled function.
#[must_use]

pub fn fourier_scaling(
    f_omega: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_div(
            Expr::BigInt(BigInt::one()),
            Expr::new_abs(a.clone()),
        ),
        Expr::Substitute(
            Arc::new(f_omega.clone()),
            out_var.to_string(),
            Arc::new(Expr::new_div(
                Expr::Variable(
                    out_var.to_string(),
                ),
                a.clone(),
            )),
        ),
    ))
}

/// Applies the differentiation property of the Fourier Transform.
///
/// If `F(ω)` is the Fourier Transform of `f(t)`, then the Fourier Transform
/// of `df(t)/dt` is `jω * F(ω)`.
///
/// # Arguments
/// * `f_omega` - The Fourier Transform of the original function `f(t)`.
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the Fourier Transform of the differentiated function.
#[must_use]

pub fn fourier_differentiation(
    f_omega: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_mul(
            i_complex(),
            Expr::Variable(
                out_var.to_string(),
            ),
        ),
        f_omega.clone(),
    ))
}

/// Applies the time-shift property of the Laplace Transform.
///
/// If `F(s)` is the Laplace Transform of `f(t)`, then the Laplace Transform
/// of `f(t - a)u(t - a)` is `e^(-as) * F(s)`, where `u(t - a)` is the unit step function.
///
/// # Arguments
/// * `f_s` - The Laplace Transform of the original function `f(t)`.
/// * `a` - The time shift amount.
/// * `out_var` - The output complex frequency variable (e.g., "s").
///
/// # Returns
/// An `Expr` representing the Laplace Transform of the time-shifted function.
#[must_use]

pub fn laplace_time_shift(
    f_s: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_exp(Expr::new_mul(
            Expr::new_neg(a.clone()),
            Expr::Variable(
                out_var.to_string(),
            ),
        )),
        f_s.clone(),
    ))
}

/// Applies the differentiation property of the Laplace Transform.
///
/// If `F(s)` is the Laplace Transform of `f(t)`, then the Laplace Transform
/// of `df(t)/dt` is `sF(s) - f(0)`.
///
/// # Arguments
/// * `f_s` - The Laplace Transform of the original function `f(t)`.
/// * `out_var` - The output complex frequency variable (e.g., "s").
/// * `f_zero` - The value of the function `f(t)` at `t=0`.
///
/// # Returns
/// An `Expr` representing the Laplace Transform of the differentiated function.
#[must_use]

pub fn laplace_differentiation(
    f_s: &Expr,
    out_var: &str,
    f_zero: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        Expr::new_mul(
            Expr::Variable(
                out_var.to_string(),
            ),
            f_s.clone(),
        ),
        f_zero.clone(),
    ))
}

/// Applies the frequency-shift property of the Laplace Transform.
///
/// If `F(s)` is the Laplace Transform of `f(t)`, then the Laplace Transform
/// of `e^(at)f(t)` is `F(s - a)`.
///
/// # Arguments
/// * `f_s` - The Laplace Transform of the original function.
/// * `a` - The shift amount.
/// * `out_var` - The output variable "s".
///
/// # Returns
/// An `Expr` representing `F(s - a)`.
#[must_use]

pub fn laplace_frequency_shift(
    f_s: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::Substitute(
        Arc::new(f_s.clone()),
        out_var.to_string(),
        Arc::new(Expr::new_sub(
            Expr::Variable(
                out_var.to_string(),
            ),
            a.clone(),
        )),
    ))
}

/// Applies the scaling property of the Laplace Transform.
///
/// If `F(s)` is the Laplace Transform of `f(t)`, then the Laplace Transform
/// of `f(at)` is `(1/a)F(s/a)`.
///
/// # Arguments
/// * `f_s` - The Laplace Transform of the original function.
/// * `a` - The scaling factor.
/// * `out_var` - The output variable "s".
///
/// # Returns
/// An `Expr` representing `(1/a)F(s/a)`.
#[must_use]

pub fn laplace_scaling(
    f_s: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_div(
            Expr::BigInt(BigInt::one()),
            a.clone(),
        ),
        Expr::Substitute(
            Arc::new(f_s.clone()),
            out_var.to_string(),
            Arc::new(Expr::new_div(
                Expr::Variable(
                    out_var.to_string(),
                ),
                a.clone(),
            )),
        ),
    ))
}

/// Applies the integration property of the Laplace Transform.
///
/// L{∫f(t)dt} = F(s)/s
#[must_use]

pub fn laplace_integration(
    f_s: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_div(
        f_s.clone(),
        Expr::Variable(
            out_var.to_string(),
        ),
    ))
}

/// Applies the time-shift property of the Z-Transform.
///
/// Z{x[n-k]} = z^(-k)X(z)
#[must_use]

pub fn z_time_shift(
    f_z: &Expr,
    k: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_pow(
            Expr::Variable(
                out_var.to_string(),
            ),
            Expr::new_neg(k.clone()),
        ),
        f_z.clone(),
    ))
}

/// Applies the scaling property of the Z-Transform.
///
/// Z{a^n x[n]} = X(z/a)
#[must_use]

pub fn z_scaling(
    f_z: &Expr,
    a: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::Substitute(
        Arc::new(f_z.clone()),
        out_var.to_string(),
        Arc::new(Expr::new_div(
            Expr::Variable(
                out_var.to_string(),
            ),
            a.clone(),
        )),
    ))
}

/// Applies the differentiation property of the Z-Transform.
///
/// Z{n x[n]} = -z dX/dz
#[must_use]

pub fn z_differentiation(
    f_z: &Expr,
    out_var: &str,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_neg(Expr::Variable(
            out_var.to_string(),
        )),
        differentiate(f_z, out_var),
    ))
}

/// Computes the partial fraction decomposition of a rational expression.
/// Computes the continuous Fourier Transform of an expression.
///
/// The Fourier Transform `F(ω)` of a function `f(t)` is defined as:
/// `F(ω) = ∫(-∞ to ∞) f(t) * e^(-jωt) dt`.
///
/// # Arguments
/// * `expr` - The expression `f(t)` to transform.
/// * `in_var` - The input time variable (e.g., "t").
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the symbolic Fourier Transform.
#[must_use]

pub fn fourier_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let integrand = Expr::new_mul(
        expr.clone(),
        Expr::new_exp(Expr::new_neg(
            Expr::new_mul(
                i_complex(),
                Expr::new_mul(
                    Expr::Variable(
                        out_var
                            .to_string(
                            ),
                    ),
                    Expr::Variable(
                        in_var
                            .to_string(
                            ),
                    ),
                ),
            ),
        )),
    );

    definite_integrate(
        &integrand,
        in_var,
        &Expr::NegativeInfinity,
        &Expr::Infinity,
    )
}

/// Computes the inverse continuous Fourier Transform of an expression.
///
/// The inverse Fourier Transform `f(t)` of `F(ω)` is defined as:
/// `f(t) = (1/(2π)) * ∫(-∞ to ∞) F(ω) * e^(jωt) dω`.
///
/// # Arguments
/// * `expr` - The expression `F(ω)` to inverse transform.
/// * `in_var` - The input frequency variable (e.g., "omega").
/// * `out_var` - The output time variable (e.g., "t").
///
/// # Returns
/// An `Expr` representing the symbolic inverse Fourier Transform.
#[must_use]

pub fn inverse_fourier_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let factor = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_mul(
            Expr::BigInt(BigInt::from(
                2,
            )),
            Expr::Variable(
                "pi".to_string(),
            ),
        ),
    );

    let integrand = Expr::new_mul(
        expr.clone(),
        Expr::new_exp(Expr::new_mul(
            i_complex(),
            Expr::new_mul(
                Expr::Variable(
                    in_var.to_string(),
                ),
                Expr::Variable(
                    out_var.to_string(),
                ),
            ),
        )),
    );

    let integral = definite_integrate(
        &integrand,
        in_var,
        &Expr::NegativeInfinity,
        &Expr::Infinity,
    );

    simplify(&Expr::new_mul(
        factor,
        integral,
    ))
}

/// Computes the unilateral Laplace Transform of an expression.
///
/// The Laplace Transform `F(s)` of a function `f(t)` is defined as:
/// `F(s) = ∫(0 to ∞) f(t) * e^(-st) dt`.
///
/// # Arguments
/// * `expr` - The expression `f(t)` to transform.
/// * `in_var` - The input time variable (e.g., "t").
/// * `out_var` - The output complex frequency variable (e.g., "s").
///
/// # Returns
/// An `Expr` representing the symbolic Laplace Transform.
#[must_use]

pub fn laplace_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let integrand = Expr::new_mul(
        expr.clone(),
        Expr::new_exp(Expr::new_neg(
            Expr::new_mul(
                Expr::Variable(
                    out_var.to_string(),
                ),
                Expr::Variable(
                    in_var.to_string(),
                ),
            ),
        )),
    );

    definite_integrate(
        &integrand,
        in_var,
        &Expr::BigInt(BigInt::zero()),
        &Expr::Infinity,
    )
}

/// Computes the inverse Laplace Transform of an expression.
///
/// The inverse Laplace Transform `f(t)` of `F(s)` is defined by the Bromwich integral:
/// `f(t) = (1/(2πj)) * ∫(c-j∞ to c+j∞) F(s) * e^(st) ds`.
///
/// This function attempts to use lookup tables and partial fraction decomposition first.
/// If these methods are insufficient, it falls back to the Bromwich integral representation
/// as a path integral.
///
/// # Arguments
/// * `expr` - The expression `F(s)` to inverse transform.
/// * `in_var` - The input complex frequency variable (e.g., "s").
/// * `out_var` - The output time variable (e.g., "t").
///
/// # Returns
/// An `Expr` representing the symbolic inverse Laplace Transform.
#[must_use]

pub fn inverse_laplace_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    if let Some(result) =
        lookup_inverse_laplace(
            expr,
            in_var,
            out_var,
        )
    {

        return result;
    }

    if let Some(terms) =
        partial_fraction_decomposition(
            expr,
            in_var,
        )
    {

        let mut result_expr =
            Expr::BigInt(BigInt::zero());

        for term in terms {

            result_expr = simplify(&Expr::new_add(
                result_expr,
                inverse_laplace_transform(
                    &term,
                    in_var,
                    out_var,
                ),
            ));
        }

        return result_expr;
    }

    let c =
        Expr::Variable("c".to_string());

    let integrand = Expr::new_mul(
        expr.clone(),
        Expr::new_exp(Expr::new_mul(
            Expr::Variable(
                in_var.to_string(),
            ),
            Expr::Variable(
                out_var.to_string(),
            ),
        )),
    );

    let factor = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_mul(
            Expr::new_mul(
                Expr::BigInt(
                    BigInt::from(2),
                ),
                Expr::Variable(
                    "pi".to_string(),
                ),
            ),
            i_complex(),
        ),
    );

    let integral = path_integrate(
        &integrand,
        in_var,
        &Expr::Path(
            crate::symbolic::core::PathType::Line,
            Arc::new(Expr::new_sub(
                c.clone(),
                Expr::Infinity,
            )),
            Arc::new(Expr::new_add(
                c,
                Expr::Infinity,
            )),
        ),
    );

    simplify(&Expr::new_mul(
        factor,
        integral,
    ))
}

/// Computes the unilateral Z-Transform of a discrete-time signal.
///
/// The Z-Transform `X(z)` of a discrete-time signal `x[n]` is defined as:
/// `X(z) = Σ(n=0 to ∞) x[n] * z^(-n)`.
///
/// # Arguments
/// * `expr` - The expression `x[n]` representing the discrete-time signal.
/// * `in_var` - The input discrete time variable (e.g., "n").
/// * `out_var` - The output complex frequency variable (e.g., "z").
///
/// # Returns
/// An `Expr` representing the symbolic Z-Transform.
#[must_use]

pub fn z_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let term = Expr::new_mul(
        expr.clone(),
        Expr::new_pow(
            Expr::Variable(
                out_var.to_string(),
            ),
            Expr::new_neg(
                Expr::Variable(
                    in_var.to_string(),
                ),
            ),
        ),
    );

    simplify(&Expr::Summation(
        Arc::new(term),
        in_var.to_string(),
        Arc::new(
            Expr::NegativeInfinity,
        ),
        Arc::new(Expr::Infinity),
    ))
}

/// Computes the inverse Z-Transform of an expression.
///
/// The inverse Z-Transform `x[n]` of `X(z)` is defined by the contour integral:
/// `x[n] = (1/(2πj)) * ∮(C) X(z) * z^(n-1) dz`.
///
/// # Arguments
/// * `expr` - The expression `X(z)` to inverse transform.
/// * `in_var` - The input complex frequency variable (e.g., "z").
/// * `out_var` - The output discrete time variable (e.g., "n").
///
/// # Returns
/// An `Expr` representing the symbolic inverse Z-Transform.
#[must_use]

pub fn inverse_z_transform(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let factor = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_mul(
            Expr::new_mul(
                Expr::BigInt(
                    BigInt::from(2),
                ),
                Expr::Variable(
                    "pi".to_string(),
                ),
            ),
            i_complex(),
        ),
    );

    let integrand = Expr::new_mul(
        expr.clone(),
        Expr::new_pow(
            Expr::Variable(
                in_var.to_string(),
            ),
            Expr::new_sub(
                Expr::Variable(
                    out_var.to_string(),
                ),
                Expr::BigInt(
                    BigInt::one(),
                ),
            ),
        ),
    );

    let integral = path_integrate(
        &integrand,
        in_var,
        &Expr::Path(
            crate::symbolic::core::PathType::Circle,
            Arc::new(Expr::BigInt(
                BigInt::zero(),
            )),
            Arc::new(Expr::Variable(
                "R".to_string(),
            )),
        ),
    );

    simplify(&Expr::new_mul(
        factor,
        integral,
    ))
}

#[must_use]

/// Performs partial fraction decomposition on a rational expression.
///
/// This function decomposes a rational function (a ratio of polynomials) into a sum of
/// simpler fractions. It handles cases with both distinct and repeated roots in the denominator.
///
/// # Arguments
/// * `expr` - The rational expression to decompose.
/// * `var` - The variable with respect to which the decomposition is performed.
///
/// # Returns
/// An `Option<Vec<Expr>>` containing the terms of the decomposition if successful, or `None` otherwise.

pub fn partial_fraction_decomposition(
    expr: &Expr,
    var: &str,
) -> Option<Vec<Expr>> {

    if let Expr::Dag(node) = expr {

        return node
            .to_expr()
            .ok()
            .and_then(|e| partial_fraction_decomposition(&e, var));
    }

    if let Expr::Div(num, den) = expr {

        let roots: Vec<Expr> =
            solve(den, var)
                .into_iter()
                .map(|r| simplify(&r))
                .collect();

        if roots.is_empty()
            || roots
                .iter()
                .any(|r| {

                    matches!(
                        r,
                        Expr::Solve(
                            _,
                            _
                        )
                    )
                })
        {

            return None;
        }

        let mut terms = Vec::new();

        let mut temp_den = den.clone();

        for root in roots {

            let factor = Expr::new_sub(
                Expr::Variable(
                    var.to_string(),
                ),
                root.clone(),
            );

            let mut multiplicity = 0;

            while is_zero(&simplify(
                &crate::symbolic::calculus::evaluate_at_point(
                    &temp_den,
                    var,
                    &root,
                ),
            )) {

                multiplicity += 1;

                let (quotient, _) = crate::symbolic::polynomial::polynomial_long_division(
                    &temp_den,
                    &factor,
                    var,
                );

                temp_den = Arc::new(quotient);
            }

            for k in 1 ..= multiplicity
            {

                let mut g = simplify(
                    &Expr::new_div(
                        num.clone(),
                        temp_den
                            .as_ref()
                            .clone(),
                    ),
                );

                for _ in 0
                    .. (multiplicity
                        - k)
                {

                    g = differentiate(
                        &g, var,
                    );
                }

                let c = simplify(&Expr::new_div(
                    crate::symbolic::calculus::evaluate_at_point(&g, var, &root),
                    Expr::Constant(crate::symbolic::calculus::factorial(multiplicity - k)),
                ));

                terms.push(simplify(
                    &Expr::new_div(
                        c,
                        Expr::new_pow(
                            factor.clone(),
                            Expr::BigInt(BigInt::from(k)),
                        ),
                    ),
                ));
            }
        }

        return Some(terms);
    }

    None
}

pub(crate) fn lookup_inverse_laplace(
    expr: &Expr,
    in_var: &str,
    out_var: &str,
) -> Option<Expr> {

    match expr {
        | Expr::Dag(node) => {
            lookup_inverse_laplace(
                &node
                    .to_expr()
                    .expect(
                        "Dag Inverse",
                    ),
                in_var,
                out_var,
            )
        },
        | Expr::Div(num, den) => {
            match (&**num, &**den) {
                | (
                    Expr::BigInt(n),
                    Expr::Variable(v),
                ) if n.is_one()
                    && v == in_var =>
                {
                    Some(Expr::BigInt(
                        BigInt::one(),
                    ))
                },
                | (
                    Expr::BigInt(n),
                    Expr::Sub(
                        s_var,
                        a_const,
                    ),
                ) if n.is_one() => {

                    if let (
                        Expr::Variable(
                            v,
                        ),
                        Expr::Constant(
                            a,
                        ),
                    ) = (
                        &**s_var,
                        &**a_const,
                    ) {

                        if v == in_var {

                            return Some(Expr::new_exp(
                                Expr::new_mul(
                                    Expr::Constant(*a),
                                    Expr::Variable(out_var.to_string()),
                                ),
                            ));
                        }
                    }

                    None
                },
                | (
                    Expr::Constant(w),
                    Expr::Add(
                        s_sq,
                        w_sq,
                    ),
                ) => {

                    if let (
                        Expr::Power(
                            s_var,
                            s_exp,
                        ),
                        Expr::Power(
                            w_const,
                            _w_exp,
                        ),
                    ) = (
                        &**s_sq,
                        &**w_sq,
                    ) {

                        if let (Expr::Variable(v), s_exp_expr) = (
                            &**s_var,
                            s_exp.clone(),
                        ) {

                            if let Expr::BigInt(s_exp_val) = &*s_exp_expr {

                                if s_exp_val == &BigInt::from(2)
                                    && v == in_var
                                    && if let Expr::Constant(val) = **w_const {

                                        val
                                    } else {

                                        return None;
                                    } == *w
                                {

                                    return Some(Expr::new_sin(
                                        Expr::new_mul(
                                            w_const.clone(),
                                            Expr::Variable(out_var.to_string()),
                                        ),
                                    ));
                                }
                            }
                        }
                    }

                    None
                },
                | (
                    Expr::Variable(v),
                    Expr::Add(
                        s_sq,
                        w_sq,
                    ),
                ) if v == in_var => {

                    if let (
                        Expr::Power(
                            s_var,
                            s_exp,
                        ),
                        Expr::Power(
                            w_const,
                            _w_exp,
                        ),
                    ) = (
                        &**s_sq,
                        &**w_sq,
                    ) {

                        if let (Expr::Variable(s), s_exp_expr) = (
                            &**s_var,
                            s_exp.clone(),
                        ) {

                            if let Expr::BigInt(s_exp_val) = &*s_exp_expr {

                                if s_exp_val == &BigInt::from(2) && s == in_var {

                                    return Some(Expr::new_cos(
                                        Expr::new_mul(
                                            w_const.clone(),
                                            Expr::Variable(out_var.to_string()),
                                        ),
                                    ));
                                }
                            }
                        }
                    }

                    None
                },
                | _ => None,
            }
        },
        | _ => None,
    }
}

/// Applies the Convolution Theorem for Fourier Transforms.
///
/// The Convolution Theorem states that the Fourier Transform of a convolution
/// of two functions is the product of their individual Fourier Transforms:
/// `FT{f(t) * g(t)} = F(ω) * G(ω)`.
///
/// # Arguments
/// * `f` - The first function in the time domain.
/// * `g` - The second function in the time domain.
/// * `in_var` - The input time variable (e.g., "t").
/// * `out_var` - The output frequency variable (e.g., "omega").
///
/// # Returns
/// An `Expr` representing the Fourier Transform of the convolution.
#[must_use]

pub fn convolution_fourier(
    f: &Expr,
    g: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let ft_f = fourier_transform(
        f,
        in_var,
        out_var,
    );

    let ft_g = fourier_transform(
        g,
        in_var,
        out_var,
    );

    simplify(&Expr::new_mul(
        ft_f, ft_g,
    ))
}

/// Applies the Convolution Theorem for Laplace Transforms.
///
/// The Convolution Theorem states that the Laplace Transform of a convolution
/// of two functions is the product of their individual Laplace Transforms:
/// `LT{f(t) * g(t)} = F(s) * G(s)`.
///
/// # Arguments
/// * `f` - The first function in the time domain.
/// * `g` - The second function in the time domain.
/// * `in_var` - The input time variable (e.g., "t").
/// * `out_var` - The output complex frequency variable (e.g., "s").
///
/// # Returns
/// An `Expr` representing the Laplace Transform of the convolution.
#[must_use]

pub fn convolution_laplace(
    f: &Expr,
    g: &Expr,
    in_var: &str,
    out_var: &str,
) -> Expr {

    let lt_f = laplace_transform(
        f,
        in_var,
        out_var,
    );

    let lt_g = laplace_transform(
        g,
        in_var,
        out_var,
    );

    simplify(&Expr::new_mul(
        lt_f, lt_g,
    ))
}
