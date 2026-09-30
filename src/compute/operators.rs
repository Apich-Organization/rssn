//! # RSSN Compute DSL Operators
//!
//! Provides clean, zero-cost DSL constructors for calculus, linear algebra,
//! and continuous analysis identity transformation operators.

use std::sync::Arc;
use crate::symbolic::core::Expr;

/// Constructs an unevaluated derivative operator `d(expr)/d(var)` in $O(1)$ memory.
#[must_use]
pub fn d(expr: Expr, var: impl Into<String>) -> Expr {
    Expr::new_derivative(expr, var.into())
}

/// Constructs an unevaluated integral operator $\int_{a}^{b} f(x) dx$ in $O(1)$ memory.
///
/// Supports definite integration with `bounds = Some((lower, upper))` and
/// indefinite integration with `bounds = None`.
#[must_use]
pub fn integral<L: Into<Expr>, U: Into<Expr>>(
    expr: Expr,
    var: impl Into<String>,
    bounds: Option<(L, U)>,
) -> Expr {
    let var_str = var.into();
    let (lb, ub) = match bounds {
        Some((l, u)) => (l.into(), u.into()),
        None => (Expr::Variable("a".to_string()), Expr::Variable("b".to_string())),
    };
    Expr::Integral {
        integrand: Arc::new(expr),
        var: Arc::new(Expr::Variable(var_str)),
        lower_bound: Arc::new(lb),
        upper_bound: Arc::new(ub),
    }
}

/// Constructs an unevaluated definite integral operator with explicit lower and upper bounds in $O(1)$ memory.
#[must_use]
pub fn definite_integral(
    expr: Expr,
    var: impl Into<String>,
    lower: impl Into<Expr>,
    upper: impl Into<Expr>,
) -> Expr {
    integral(expr, var, Some((lower.into(), upper.into())))
}

/// Constructs an unevaluated indefinite integral operator in $O(1)$ memory.
#[must_use]
pub fn indefinite_integral(
    expr: Expr,
    var: impl Into<String>,
) -> Expr {
    integral::<Expr, Expr>(expr, var, None)
}

/// Constructs an unevaluated limit operator $\lim_{var \to to} f(var)$ in $O(1)$ memory.
#[must_use]
pub fn limit(
    expr: Expr,
    var: impl Into<String>,
    to: impl Into<Expr>,
) -> Expr {
    Expr::Limit(
        Arc::new(expr),
        var.into(),
        Arc::new(to.into()),
    )
}

/// Constructs an unevaluated gradient operator $\nabla f$ in $O(1)$ memory.
#[must_use]
pub fn gradient(
    expr: Expr,
    vars: &[impl AsRef<str>],
) -> Expr {
    let mut args = Vec::with_capacity(vars.len() + 1);
    args.push(expr);
    for v in vars {
        args.push(Expr::Variable(v.as_ref().to_string()));
    }
    Expr::NaryList("gradient".to_string(), args)
}

/// Constructs a symbolic variable expression in $O(1)$ memory.
#[must_use]
pub fn var(name: impl Into<String>) -> Expr {
    let s = name.into();
    Expr::new_variable(&s)
}

/// Constructs an unevaluated ordinary differential equation (ODE) operator in $O(1)$ memory.
#[must_use]
pub fn ode(
    equation: Expr,
    func: impl Into<String>,
    var: impl Into<String>,
) -> Expr {
    Expr::Ode {
        equation: Arc::new(equation),
        func: func.into(),
        var: var.into(),
    }
}

/// Constructs an unevaluated equation solving operator in $O(1)$ memory.
#[must_use]
pub fn solve(equation: Expr, target: impl Into<String>) -> Expr {
    Expr::Solve(Arc::new(equation), target.into())
}

/// Constructs an equation operator `lhs = rhs` in $O(1)$ memory.
#[must_use]
pub fn eq(lhs: Expr, rhs: Expr) -> Expr {
    Expr::Eq(Arc::new(lhs), Arc::new(rhs))
}

/// Constructs an unevaluated matrix inverse operator $A^{-1}$ in $O(1)$ memory.
#[must_use]
pub fn matrix_inv(matrix: Expr) -> Expr {
    Expr::Inverse(Arc::new(matrix))
}

/// Constructs an unevaluated matrix determinant operator $\det(A)$ in $O(1)$ memory.
#[must_use]
pub fn det(matrix: Expr) -> Expr {
    Expr::UnaryList("det".to_string(), Arc::new(matrix))
}

/// Constructs an unevaluated matrix eigenvalues operator in $O(1)$ memory.
#[must_use]
pub fn eigenvalues(matrix: Expr) -> Expr {
    Expr::UnaryList("eigenvalues".to_string(), Arc::new(matrix))
}

/// Constructs an unevaluated characteristic polynomial operator $\det(A - \lambda I)$ in $O(1)$ memory.
#[must_use]
pub fn charpoly(matrix: Expr, var: impl Into<String>) -> Expr {
    Expr::BinaryList(
        "characteristic_polynomial".to_string(),
        Arc::new(matrix),
        Arc::new(Expr::Variable(var.into())),
    )
}

// Re-export common elementary functions for calculus DSL convenience
pub use crate::symbolic::elementary::{cos, cosh, exp, ln, sin, sinh, tan, tanh};
