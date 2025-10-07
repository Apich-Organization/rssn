//! # Radical Expression Simplification
//!
//! This module provides functions for simplifying radical expressions, particularly
//! focusing on the denesting of nested square roots of the form `sqrt(A + B*sqrt(C))`.

use std::sync::Arc;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;

/// Attempts to denest a nested square root of the form `sqrt(A + B*sqrt(C))`.
///
/// This function applies the denesting formula: `sqrt(X ± sqrt(Y)) = sqrt((X+sqrt(X^2-Y))/2) ± sqrt((X-sqrt(X^2-Y))/2)`.
/// It looks for a pattern `A + B*sqrt(C)` inside the outermost square root.
///
/// # Arguments
/// * `expr` - The expression containing the nested square root.
///
/// # Returns
/// The simplified expression if denesting is successful, or the original expression if no simplification is found.
pub fn denest_sqrt(expr: &Expr) -> Expr {
    if let Expr::Sqrt(inner) = expr {
        if let Some((a, b, c)) = match_nested_sqrt_pattern(inner) {
            // We are trying to simplify sqrt(a + b*sqrt(c))
            // We look for a solution of the form sqrt(x) + sqrt(y)
            // This leads to solving the quadratic z^2 - a*z + (b^2*c)/4 = 0
            let discriminant = simplify(Expr::Sub(
                Arc::new(Expr::Power(
                    Arc::new(a.clone()),
                    Arc::new(Expr::Constant(2.0)),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(b.clone()),
                    Arc::new(Expr::Power(
                        Arc::new(c.clone()),
                        Arc::new(Expr::Constant(2.0)),
                    )),
                )), // This should be b^2*c
            ));

            if let Some(alpha) = is_perfect_square(&discriminant) {
                // The roots are (a ± alpha) / 2
                let two = Expr::Constant(2.0);
                let x = simplify(Expr::Div(
                    Arc::new(Expr::Add(Arc::new(a.clone()), Arc::new(alpha.clone()))),
                    Arc::new(two.clone()),
                ));
                let y = simplify(Expr::Div(
                    Arc::new(Expr::Sub(Arc::new(a), Arc::new(alpha))),
                    Arc::new(two),
                ));

                return simplify(Expr::Add(
                    Arc::new(Expr::Sqrt(Arc::new(x))),
                    Arc::new(Expr::Sqrt(Arc::new(y))),
                ));
            }
        }
    }
    expr.clone()
}

/// Matches an expression of the form A + B*sqrt(C).
pub(crate) fn match_nested_sqrt_pattern(expr: &Expr) -> Option<(Expr, Expr, Expr)> {
    if let Expr::Add(a, term_b) = expr {
        if let Expr::Mul(b, sqrt_c) = &**term_b {
            if let Expr::Sqrt(c) = &**sqrt_c {
                return Some((a.as_ref().clone(), b.as_ref().clone(), c.as_ref().clone()));
            }
        }
    }
    None
}

/// Checks if an expression is a perfect square and returns its root if so.
pub(crate) fn is_perfect_square(expr: &Expr) -> Option<Expr> {
    if let Expr::Constant(c) = expr {
        if *c >= 0.0 {
            let root = c.sqrt();
            if root.fract() == 0.0 {
                return Some(Expr::Constant(root));
            }
        }
    }
    // More cases for BigInt, Rational, and symbolic squares can be added here.
    None
}
