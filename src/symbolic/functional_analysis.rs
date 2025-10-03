//! # Functional Analysis
//!
//! This module provides structures and functions for computations in functional analysis.
//! Functional analysis is a branch of mathematics that studies vector spaces endowed with
//! some kind of limit-related structure (e.g., an inner product, a norm) and the linear
//! functions that act upon these spaces. It includes implementations for Hilbert and Banach
//! spaces, linear operators, inner products, and various norms.

use crate::symbolic::calculus::{definite_integrate, differentiate};
use crate::symbolic::core::Expr;
use crate::symbolic::elementary::sqrt;
use crate::symbolic::simplify::{is_zero, simplify};
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Represents a Hilbert space, a complete inner product space.
/// This implementation specifically models L^2([a, b]), the space of square-integrable
/// complex-valued functions on an interval [a, b].
#[derive(Clone, Debug, PartialEq)]
pub struct HilbertSpace {
    /// The variable of the functions in this space, e.g., "x".
    pub var: String,
    /// The lower bound of the integration interval.
    pub lower_bound: Expr,
    /// The upper bound of the integration interval.
    pub upper_bound: Expr,
}

impl HilbertSpace {
    /// Creates a new L^2 space on the interval `[a, b]`.
    ///
    /// # Arguments
    /// * `var` - The variable of the functions in this space, e.g., "x".
    /// * `lower_bound` - The lower bound of the integration interval.
    /// * `upper_bound` - The upper bound of the integration interval.
    ///
    /// # Returns
    /// A new `HilbertSpace` instance.
    pub fn new(var: &str, lower_bound: Expr, upper_bound: Expr) -> Self {
        Self {
            var: var.to_string(),
            lower_bound,
            upper_bound,
        }
    }
}

/// Represents a Banach space, a complete normed vector space.
/// This implementation specifically models L^p([a, b]), the space of functions for which
/// the p-th power of their absolute value is Lebesgue integrable.
#[derive(Clone, Debug, PartialEq)]
pub struct BanachSpace {
    /// The variable of the functions in this space, e.g., "x".
    pub var: String,
    /// The lower bound of the integration interval.
    pub lower_bound: Expr,
    /// The upper bound of the integration interval.
    pub upper_bound: Expr,
    /// The p-value for the L^p norm, where p >= 1.
    pub p: Expr,
}

impl BanachSpace {
    /// Creates a new L^p space on the interval `[a, b]`.
    ///
    /// # Arguments
    /// * `var` - The variable of the functions in this space, e.g., "x".
    /// * `lower_bound` - The lower bound of the integration interval.
    /// * `upper_bound` - The upper bound of the integration interval.
    /// * `p` - The p-value for the L^p norm, where `p >= 1`.
    ///
    /// # Returns
    /// A new `BanachSpace` instance.
    pub fn new(var: &str, lower_bound: Expr, upper_bound: Expr, p: Expr) -> Self {
        Self {
            var: var.to_string(),
            lower_bound,
            upper_bound,
            p,
        }
    }
}

/// Represents common linear operators that act on functions in a vector space.
#[derive(Clone, Debug, PartialEq)]
pub enum LinearOperator {
    /// The derivative operator d/dx.
    Derivative(String),
    /// An integral operator ∫_a^x, where a is the lower bound.
    Integral(Expr, String),
}

impl LinearOperator {
    /// Applies the operator to a given expression (function).
    ///
    /// # Arguments
    /// * `expr` - The expression to apply the operator to.
    ///
    /// # Returns
    /// A new `Expr` representing the result of the operation.
    pub fn apply(&self, expr: &Expr) -> Expr {
        match self {
            LinearOperator::Derivative(var) => differentiate(expr, var),
            LinearOperator::Integral(lower_bound, var) => {
                let x = Expr::Variable(var.clone());
                definite_integrate(expr, var, lower_bound, &x)
            }
        }
    }
}

/// Computes the inner product of two functions, `f` and `g`, in a given Hilbert space.
///
/// For the L^2([a, b]) space, the inner product is defined as:
/// `<f, g> = ∫_a^b f(x)g*(x) dx`.
/// For simplicity with real functions, this implementation computes `∫_a^b f(x)g(x) dx`.
///
/// # Arguments
/// * `space` - The `HilbertSpace` defining the integration interval.
/// * `f` - The first function.
/// * `g` - The second function.
///
/// # Returns
/// An `Expr` representing the symbolic result of the inner product integral.
pub fn inner_product(space: &HilbertSpace, f: &Expr, g: &Expr) -> Expr {
    let integrand = simplify(Expr::Mul(Box::new(f.clone()), Box::new(g.clone())));
    definite_integrate(
        &integrand,
        &space.var,
        &space.lower_bound,
        &space.upper_bound,
    )
}

/// Computes the norm of a function `f` in a given Hilbert space.
///
/// The norm is a measure of the "length" of the function and is induced by the inner product:
/// `||f|| = sqrt(<f, f>)`.
///
/// # Arguments
/// * `space` - The `HilbertSpace`.
/// * `f` - The function.
///
/// # Returns
/// An `Expr` representing the norm of the function.
pub fn norm(space: &HilbertSpace, f: &Expr) -> Expr {
    let inner_product_f_f = inner_product(space, f, f);
    sqrt(inner_product_f_f)
}

/// Computes the L^p norm of a function `f` in a given Banach space.
///
/// The L^p norm is defined as: `||f||_p = (∫_a^b |f(x)|^p dx)^(1/p)`.
///
/// # Arguments
/// * `space` - The `BanachSpace` defining the interval and p-value.
/// * `f` - The function.
///
/// # Returns
/// An `Expr` representing the L^p norm of the function.
pub fn banach_norm(space: &BanachSpace, f: &Expr) -> Expr {
    // Integrand is |f(x)|^p
    let integrand = Expr::Power(
        Box::new(Expr::Abs(Box::new(f.clone()))),
        Box::new(space.p.clone()),
    );

    // Integral part: ∫_a^b |f(x)|^p dx
    let integral = definite_integrate(
        &integrand,
        &space.var,
        &space.lower_bound,
        &space.upper_bound,
    );

    // Final result: ( ... )^(1/p)
    let one_over_p = Expr::Div(
        Box::new(Expr::BigInt(BigInt::one())),
        Box::new(space.p.clone()),
    );

    simplify(Expr::Power(Box::new(integral), Box::new(one_over_p)))
}

/// Checks if two functions are orthogonal in a given Hilbert space.
///
/// Two functions are orthogonal if their inner product is zero.
///
/// # Arguments
/// * `space` - The `HilbertSpace`.
/// * `f` - The first function.
/// * `g` - The second function.
///
/// # Returns
/// `true` if the functions are orthogonal, `false` otherwise.
pub fn are_orthogonal(space: &HilbertSpace, f: &Expr, g: &Expr) -> bool {
    let prod = simplify(inner_product(space, f, g));
    is_zero(&prod)
}

/// Computes the projection of function `f` onto function `g` in a given Hilbert space.
///
/// The projection of `f` onto `g` finds the component of `f` that lies in the direction of `g`.
/// Formula: `proj_g(f) = (<f, g> / <g, g>) * g`.
///
/// # Arguments
/// * `space` - The `HilbertSpace`.
/// * `f` - The function to project.
/// * `g` - The function to project onto.
///
/// # Returns
/// An `Expr` representing the projected function. Returns a zero expression if `g` is the zero vector.
pub fn project(space: &HilbertSpace, f: &Expr, g: &Expr) -> Expr {
    let inner_product_f_g = inner_product(space, f, g);
    let inner_product_g_g = inner_product(space, g, g);

    // If the norm of g is zero, g is the zero vector, so the projection is zero.
    if is_zero(&simplify(inner_product_g_g.clone())) {
        return Expr::BigInt(num_bigint::BigInt::zero());
    }

    let coefficient = simplify(Expr::Div(
        Box::new(inner_product_f_g),
        Box::new(inner_product_g_g),
    ));

    simplify(Expr::Mul(Box::new(coefficient), Box::new(g.clone())))
}
