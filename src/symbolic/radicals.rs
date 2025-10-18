//! # Radical Expression Simplification
//!
//! This module provides functions for simplifying radical expressions, particularly
//! focusing on the denesting of nested square roots of the form `sqrt(A + B*sqrt(C))`.
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use std::sync::Arc;
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
            let discriminant = simplify(Expr::new_sub(
                Expr::new_pow(a.clone(), Expr::Constant(2.0)),
                Expr::new_mul(b.clone(), Expr::new_pow(c.clone(), Expr::Constant(2.0))),
            ));
            if let Some(alpha) = is_perfect_square(&discriminant) {
                let two = Expr::Constant(2.0);
                let x = simplify(Expr::new_div(
                    Expr::new_add(a.clone(), alpha.clone()),
                    two.clone(),
                ));
                let y = simplify(Expr::new_div(Expr::new_sub(a, alpha), two));
                return simplify(Expr::new_add(Expr::new_sqrt(x), Expr::new_sqrt(y)));
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
    None
}
