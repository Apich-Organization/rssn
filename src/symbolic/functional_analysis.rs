//! # Functional Analysis
//!
//! This module provides structures and functions for computations in functional analysis.
//! Functional analysis is a branch of mathematics that studies vector spaces endowed with
//! some kind of limit-related structure (e.g., an inner product, a norm) and the linear
//! functions that act upon these spaces. It includes implementations for Hilbert and Banach
//! spaces, linear operators, inner products, and various norms.

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::calculus::definite_integrate;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::elementary::sqrt;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;

/// Represents a Hilbert space, a complete inner product space.
/// This implementation specifically models L^2([a, b]), the space of square-integrable
/// complex-valued functions on an interval [a, b].
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

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
    #[must_use]

    pub fn new(
        var: &str,
        lower_bound: Expr,
        upper_bound: Expr,
    ) -> Self {

        Self {
            var: var.to_string(),
            lower_bound,
            upper_bound,
        }
    }
}

/// Represents a Banach space, a complete normed vector space.
///
/// This implementation specifically models L^p([a, b]), the space of functions for which
/// the p-th power of their absolute value is Lebesgue integrable.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

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
    #[must_use]

    pub fn new(
        var: &str,
        lower_bound: Expr,
        upper_bound: Expr,
        p: Expr,
    ) -> Self {

        Self {
            var: var.to_string(),
            lower_bound,
            upper_bound,
            p,
        }
    }
}

/// Represents common linear operators that act on functions in a vector space.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub enum LinearOperator {
    /// The identity operator I(f) = f.
    Identity,
    /// The derivative operator d/dx.
    Derivative(String),
    /// An integral operator ∫_a^x, where a is the lower bound.
    Integral(Expr, String),
    /// A multiplication operator `M_g(f)` = g * f.
    Multiplication(Expr),
    /// A composition of two operators (A ∘ B).
    Composition(
        Box<LinearOperator>,
        Box<LinearOperator>,
    ),
}

impl LinearOperator {
    /// Applies the operator to a given expression (function).
    #[must_use]

    pub fn apply(
        &self,
        expr: &Expr,
    ) -> Expr {

        match self {
            | Self::Identity => {
                expr.clone()
            },
            | Self::Derivative(var) => {
                differentiate(expr, var)
            },
            | Self::Integral(
                lower_bound,
                var,
            ) => {

                let x = Expr::Variable(
                    var.clone(),
                );

                definite_integrate(
                    expr,
                    var,
                    lower_bound,
                    &x,
                )
            },
            | Self::Multiplication(
                g,
            ) => {
                simplify(
                    &Expr::new_mul(
                        g.clone(),
                        expr.clone(),
                    ),
                )
            },
            | Self::Composition(
                op1,
                op2,
            ) => {
                op1.apply(
                    &op2.apply(expr),
                )
            },
        }
    }
}

/// Computes the inner product of two functions, `f` and `g`, in a given Hilbert space.
///
/// For the L^2([a, b]) space, the inner product is defined as:
/// `<f, g> = ∫_a^b f(x)g*(x) dx`.
/// For simplicity with real functions, this implementation computes `∫_a^b f(x)g(x) dx`.
#[must_use]

pub fn inner_product(
    space: &HilbertSpace,
    f: &Expr,
    g: &Expr,
) -> Expr {

    let integrand =
        simplify(&Expr::new_mul(
            f.clone(),
            g.clone(),
        ));

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
#[must_use]

pub fn norm(
    space: &HilbertSpace,
    f: &Expr,
) -> Expr {

    let inner_product_f_f =
        inner_product(space, f, f);

    sqrt(inner_product_f_f)
}

/// Computes the L^p norm of a function `f` in a given Banach space.
///
/// The L^p norm is defined as: `||f||_p = (∫_a^b |f(x)|^p dx)^(1/p)`.
#[must_use]

pub fn banach_norm(
    space: &BanachSpace,
    f: &Expr,
) -> Expr {

    let integrand = Expr::new_pow(
        Expr::new_abs(f.clone()),
        space.p.clone(),
    );

    let integral = definite_integrate(
        &integrand,
        &space.var,
        &space.lower_bound,
        &space.upper_bound,
    );

    let one_over_p = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        space.p.clone(),
    );

    simplify(&Expr::new_pow(
        integral,
        one_over_p,
    ))
}

/// Checks if two functions are orthogonal in a given Hilbert space.
///
/// Two functions are orthogonal if their inner product is zero.
#[must_use]

pub fn are_orthogonal(
    space: &HilbertSpace,
    f: &Expr,
    g: &Expr,
) -> bool {

    let prod = simplify(
        &inner_product(space, f, g),
    );

    is_zero(&prod)
}

/// Computes the projection of function `f` onto function `g` in a given Hilbert space.
///
/// The projection of `f` onto `g` finds the component of `f` that lies in the direction of `g`.
/// Formula: `proj_g(f) = (<f, g> / <g, g>) * g`.
#[must_use]

pub fn project(
    space: &HilbertSpace,
    f: &Expr,
    g: &Expr,
) -> Expr {

    let inner_product_f_g =
        inner_product(space, f, g);

    let inner_product_g_g =
        inner_product(space, g, g);

    if is_zero(&simplify(
        &inner_product_g_g,
    )) {

        return Expr::BigInt(
            num_bigint::BigInt::zero(),
        );
    }

    let coefficient =
        simplify(&Expr::new_div(
            inner_product_f_g,
            inner_product_g_g,
        ));

    simplify(&Expr::new_mul(
        coefficient,
        g.clone(),
    ))
}

/// Performs the Gram-Schmidt process to orthogonalize a set of functions.
///
/// Given a set of linearly independent functions {`u_1`, ..., `u_n`}, this function produces
/// a set of orthogonal functions {`v_1`, ..., `v_n`} that span the same subspace.
///
/// `v_1` = `u_1`
/// `v_k` = `u_k` - sum_{j=1}^{k-1} proj_{`v_j}(u_k)`
#[must_use]

pub fn gram_schmidt(
    space: &HilbertSpace,
    basis: &[Expr],
) -> Vec<Expr> {

    let mut orthogonal_basis =
        Vec::new();

    for b in basis {

        let mut v = b.clone();

        for u in &orthogonal_basis {

            let proj = project(
                space,
                b,
                u,
            );

            v = Expr::new_sub(v, proj);
        }

        orthogonal_basis
            .push(simplify(&v));
    }

    orthogonal_basis
}

/// Performs the Gram-Schmidt process to orthonormalize a set of functions.
///
/// This is similar to `gram_schmidt`, but each resulting vector is normalized to have length 1.
#[must_use]

pub fn gram_schmidt_orthonormal(
    space: &HilbertSpace,
    basis: &[Expr],
) -> Vec<Expr> {

    let orthogonal_basis =
        gram_schmidt(space, basis);

    let mut orthonormal_basis =
        Vec::new();

    for v in orthogonal_basis {

        let n = norm(space, &v);

        if is_zero(&n) {

            // Should not happen if basis is linearly independent
            orthonormal_basis.push(v);
        } else {

            orthonormal_basis.push(
                simplify(
                    &Expr::new_div(
                        v, n,
                    ),
                ),
            );
        }
    }

    orthonormal_basis
}
