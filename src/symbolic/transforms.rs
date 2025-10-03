//! # Symbolic Integral Transforms
//!
//! This module provides functions for performing symbolic integral transforms,
//! including the Fourier, Laplace, and Z-transforms, as well as their inverses.
//! It also includes implementations of key transform properties and theorems,
//! such as the convolution theorem.

use crate::symbolic::calculus::{definite_integrate, differentiate, path_integrate};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::{is_zero, simplify};
use crate::symbolic::solve::solve;
use num_bigint::BigInt;
use num_traits::{One, Zero};

pub(crate) fn i_complex() -> Expr {
    Expr::Complex(
        Box::new(Expr::BigInt(BigInt::zero())),
        Box::new(Expr::BigInt(BigInt::one())),
    )
}

// =====================================================================================
// region: Transform Properties
// =====================================================================================

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
pub fn fourier_time_shift(f_omega: &Expr, a: &Expr, out_var: &str) -> Expr {
    simplify(Expr::Mul(
        Box::new(Expr::Exp(Box::new(Expr::Mul(
            Box::new(Expr::Mul(
                Box::new(i_complex()),
                Box::new(Expr::Variable(out_var.to_string())),
            )),
            Box::new(Expr::Neg(Box::new(a.clone()))),
        )))),
        Box::new(f_omega.clone()),
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
pub fn fourier_frequency_shift(f_omega: &Expr, a: &Expr, out_var: &str) -> Expr {
    simplify(Expr::Substitute(
        Box::new(f_omega.clone()),
        out_var.to_string(),
        Box::new(Expr::Sub(
            Box::new(Expr::Variable(out_var.to_string())),
            Box::new(a.clone()),
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
pub fn fourier_scaling(f_omega: &Expr, a: &Expr, out_var: &str) -> Expr {
    simplify(Expr::Mul(
        Box::new(Expr::Div(
            Box::new(Expr::BigInt(BigInt::one())),
            Box::new(Expr::Abs(Box::new(a.clone()))),
        )),
        Box::new(Expr::Substitute(
            Box::new(f_omega.clone()),
            out_var.to_string(),
            Box::new(Expr::Div(
                Box::new(Expr::Variable(out_var.to_string())),
                Box::new(a.clone()),
            )),
        )),
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
pub fn fourier_differentiation(f_omega: &Expr, out_var: &str) -> Expr {
    simplify(Expr::Mul(
        Box::new(Expr::Mul(
            Box::new(i_complex()),
            Box::new(Expr::Variable(out_var.to_string())),
        )),
        Box::new(f_omega.clone()),
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
pub fn laplace_time_shift(f_s: &Expr, a: &Expr, out_var: &str) -> Expr {
    simplify(Expr::Mul(
        Box::new(Expr::Exp(Box::new(Expr::Mul(
            Box::new(Expr::Neg(Box::new(a.clone()))),
            Box::new(Expr::Variable(out_var.to_string())),
        )))),
        Box::new(f_s.clone()),
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
pub fn laplace_differentiation(f_s: &Expr, out_var: &str, f_zero: &Expr) -> Expr {
    simplify(Expr::Sub(
        Box::new(Expr::Mul(
            Box::new(Expr::Variable(out_var.to_string())),
            Box::new(f_s.clone()),
        )),
        Box::new(f_zero.clone()),
    ))
}

// endregion

// =====================================================================================
// region: Main Transform Functions
// =====================================================================================

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
pub fn fourier_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    let integrand = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
            Box::new(i_complex()),
            Box::new(Expr::Mul(
                Box::new(Expr::Variable(out_var.to_string())),
                Box::new(Expr::Variable(in_var.to_string())),
            )),
        )))))),
    );
    definite_integrate(&integrand, in_var, &Expr::NegativeInfinity, &Expr::Infinity)
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
pub fn inverse_fourier_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    let factor = Expr::Div(
        Box::new(Expr::BigInt(BigInt::one())),
        Box::new(Expr::Mul(
            Box::new(Expr::BigInt(BigInt::from(2))),
            Box::new(Expr::Variable("pi".to_string())),
        )),
    );
    let integrand = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Exp(Box::new(Expr::Mul(
            Box::new(i_complex()),
            Box::new(Expr::Mul(
                Box::new(Expr::Variable(in_var.to_string())),
                Box::new(Expr::Variable(out_var.to_string())),
            )),
        )))),
    );
    let integral = definite_integrate(&integrand, in_var, &Expr::NegativeInfinity, &Expr::Infinity);
    simplify(Expr::Mul(Box::new(factor), Box::new(integral)))
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
pub fn laplace_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    let integrand = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
            Box::new(Expr::Variable(out_var.to_string())),
            Box::new(Expr::Variable(in_var.to_string())),
        )))))),
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
pub fn inverse_laplace_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    if let Some(result) = lookup_inverse_laplace(expr, in_var, out_var) {
        return result;
    }
    if let Some(terms) = partial_fraction_decomposition(expr, in_var) {
        let mut result_expr = Expr::BigInt(BigInt::zero());
        for term in terms {
            result_expr = simplify(Expr::Add(
                Box::new(result_expr),
                Box::new(inverse_laplace_transform(&term, in_var, out_var)),
            ));
        }
        return result_expr;
    }
    let c = Expr::Variable("c".to_string());
    let integrand = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Exp(Box::new(Expr::Mul(
            Box::new(Expr::Variable(in_var.to_string())),
            Box::new(Expr::Variable(out_var.to_string())),
        )))),
    );
    let factor = Expr::Div(
        Box::new(Expr::BigInt(BigInt::one())),
        Box::new(Expr::Mul(
            Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Variable("pi".to_string())),
            )),
            Box::new(i_complex()),
        )),
    );
    let integral = path_integrate(
        &integrand,
        in_var,
        &Expr::Path(
            crate::symbolic::core::PathType::Line,
            Box::new(Expr::Sub(Box::new(c.clone()), Box::new(Expr::Infinity))),
            Box::new(Expr::Add(Box::new(c), Box::new(Expr::Infinity))),
        ),
    );
    simplify(Expr::Mul(Box::new(factor), Box::new(integral)))
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
pub fn z_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    let term = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Power(
            Box::new(Expr::Variable(out_var.to_string())),
            Box::new(Expr::Neg(Box::new(Expr::Variable(in_var.to_string())))),
        )),
    );
    simplify(Expr::Summation(
        Box::new(term),
        in_var.to_string(),
        Box::new(Expr::NegativeInfinity),
        Box::new(Expr::Infinity),
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
pub fn inverse_z_transform(expr: &Expr, in_var: &str, out_var: &str) -> Expr {
    let factor = Expr::Div(
        Box::new(Expr::BigInt(BigInt::one())),
        Box::new(Expr::Mul(
            Box::new(Expr::Mul(
                Box::new(Expr::BigInt(BigInt::from(2))),
                Box::new(Expr::Variable("pi".to_string())),
            )),
            Box::new(i_complex()),
        )),
    );
    let integrand = Expr::Mul(
        Box::new(expr.clone()),
        Box::new(Expr::Power(
            Box::new(Expr::Variable(in_var.to_string())),
            Box::new(Expr::Sub(
                Box::new(Expr::Variable(out_var.to_string())),
                Box::new(Expr::BigInt(BigInt::one())),
            )),
        )),
    );
    let integral = path_integrate(
        &integrand,
        in_var,
        &Expr::Path(
            crate::symbolic::core::PathType::Circle,
            Box::new(Expr::BigInt(BigInt::zero())),
            Box::new(Expr::Variable("R".to_string())),
        ),
    );
    simplify(Expr::Mul(Box::new(factor), Box::new(integral)))
}

// endregion

// =====================================================================================
// region: Helpers and Theorems
// =====================================================================================

pub(crate) fn partial_fraction_decomposition(expr: &Expr, var: &str) -> Option<Vec<Expr>> {
    if let Expr::Div(num, den) = expr {
        let roots = solve(den, var);
        if roots.is_empty() || roots.iter().any(|r| matches!(r, Expr::Solve(_, _))) {
            return None; // Can't solve for roots
        }

        let mut terms = Vec::new();
        let mut temp_den = den.clone();

        for root in roots {
            let factor = Expr::Sub(
                Box::new(Expr::Variable(var.to_string())),
                Box::new(root.clone()),
            );
            let mut multiplicity = 0;
            while is_zero(&simplify(crate::symbolic::calculus::evaluate_at_point(
                &temp_den, var, &root,
            ))) {
                multiplicity += 1;
                temp_den = Box::new(simplify(Expr::Div(
                    Box::new(*temp_den),
                    Box::new(factor.clone()),
                )));
            }

            for k in 1..=multiplicity {
                let mut g = simplify(Expr::Div(num.clone(), Box::new(*temp_den.clone())));
                for _ in 0..(multiplicity - k) {
                    g = differentiate(&g, var);
                }
                let c = simplify(Expr::Div(
                    Box::new(crate::symbolic::calculus::evaluate_at_point(&g, var, &root)),
                    Box::new(Expr::Constant(crate::symbolic::calculus::factorial(
                        multiplicity - k,
                    ))),
                ));
                terms.push(simplify(Expr::Div(
                    Box::new(c),
                    Box::new(Expr::Power(
                        Box::new(factor.clone()),
                        Box::new(Expr::BigInt(BigInt::from(k))),
                    )),
                )));
            }
        }
        return Some(terms);
    }
    None
}

pub(crate) fn lookup_inverse_laplace(expr: &Expr, in_var: &str, out_var: &str) -> Option<Expr> {
    match expr {
        Expr::Div(num, den) => match (&**num, &**den) {
            (Expr::BigInt(n), Expr::Variable(v)) if n.is_one() && v == in_var => {
                Some(Expr::BigInt(BigInt::one()))
            }
            (Expr::BigInt(n), Expr::Sub(s_var, a_const)) if n.is_one() => {
                if let (Expr::Variable(v), Expr::Constant(a)) = (&**s_var, &**a_const) {
                    if v == in_var {
                        return Some(Expr::Exp(Box::new(Expr::Mul(
                            Box::new(Expr::Constant(*a)),
                            Box::new(Expr::Variable(out_var.to_string())),
                        ))));
                    }
                }
                None
            }
            (Expr::Constant(w), Expr::Add(s_sq, w_sq)) => {
                if let (Expr::Power(s_var, s_exp), Expr::Power(w_const, _w_exp)) =
                    (&**s_sq, &**w_sq)
                {
                    if let (Expr::Variable(v), s_exp_expr) = (&**s_var, s_exp.clone()) {
                        if let Expr::BigInt(s_exp_val) = &*s_exp_expr {
                            if s_exp_val == &BigInt::from(2)
                                && v == in_var
                                && if let Expr::Constant(val) = **w_const {
                                    val
                                } else {
                                    return None;
                                } == *w
                            {
                                return Some(Expr::Sin(Box::new(Expr::Mul(
                                    w_const.clone(),
                                    Box::new(Expr::Variable(out_var.to_string())),
                                ))));
                            }
                        }
                    }
                }
                None
            }
            (Expr::Variable(v), Expr::Add(s_sq, w_sq)) if v == in_var => {
                if let (Expr::Power(s_var, s_exp), Expr::Power(w_const, _w_exp)) =
                    (&**s_sq, &**w_sq)
                {
                    if let (Expr::Variable(s), s_exp_expr) = (&**s_var, s_exp.clone()) {
                        if let Expr::BigInt(s_exp_val) = &*s_exp_expr {
                            if s_exp_val == &BigInt::from(2) && s == in_var {
                                return Some(Expr::Cos(Box::new(Expr::Mul(
                                    w_const.clone(),
                                    Box::new(Expr::Variable(out_var.to_string())),
                                ))));
                            }
                        }
                    }
                }
                None
            }
            _ => None,
        },
        _ => None,
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
pub fn convolution_fourier(f: &Expr, g: &Expr, in_var: &str, out_var: &str) -> Expr {
    let ft_f = fourier_transform(f, in_var, out_var);
    let ft_g = fourier_transform(g, in_var, out_var);
    simplify(Expr::Mul(Box::new(ft_f), Box::new(ft_g)))
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
pub fn convolution_laplace(f: &Expr, g: &Expr, in_var: &str, out_var: &str) -> Expr {
    let lt_f = laplace_transform(f, in_var, out_var);
    let lt_g = laplace_transform(g, in_var, out_var);
    simplify(Expr::Mul(Box::new(lt_f), Box::new(lt_g)))
}

// endregion
