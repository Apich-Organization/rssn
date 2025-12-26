//! # Radical Expression Simplification
//!
//! This module provides functions for simplifying radical expressions, particularly
//! focusing on the denesting of nested square roots of the form `sqrt(A + B*sqrt(C))`.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

/// Recursively simplifies radical expressions in the given expression tree.
///
/// This function traverses the expression and applies denesting algorithms to
/// any square root sub-expressions it encounters.
///
/// # Arguments
/// * `expr` - The expression to simplify.
///
/// # Returns
/// The simplified expression with denested radicals where possible.
#[must_use]

pub fn simplify_radicals(
    expr : &Expr
) -> Expr {

    match expr {
        | Expr::Sqrt(inner) => {

            let simplified_inner =
                simplify_radicals(
                    inner,
                );

            denest_sqrt(
                &Expr::new_sqrt(
                    simplified_inner,
                ),
            )
        },
        | Expr::Power(base, exp) => {

            let simplified_base =
                simplify_radicals(base);

            let simplified_exp =
                simplify_radicals(exp);

            // Check if this is a square root (power of 1/2)
            if let Expr::Constant(c) =
                &simplified_exp
            {

                if (c - 0.5).abs()
                    < f64::EPSILON
                {

                    return denest_sqrt(&Expr::new_sqrt(
                        simplified_base,
                    ));
                }
            }

            Expr::new_pow(
                simplified_base,
                simplified_exp,
            )
        },
        | Expr::Add(a, b) => {
            Expr::new_add(
                simplify_radicals(a),
                simplify_radicals(b),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                simplify_radicals(a),
                simplify_radicals(b),
            )
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                simplify_radicals(a),
                simplify_radicals(b),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                simplify_radicals(a),
                simplify_radicals(b),
            )
        },
        | Expr::Neg(a) => {
            Expr::new_neg(
                simplify_radicals(a),
            )
        },
        | Expr::Dag(node) => {
            simplify_radicals(
                &node
                    .to_expr()
                    .unwrap_or_else(
                        |_| {
                            expr.clone()
                        },
                    ),
            )
        },
        | _ => expr.clone(),
    }
}

/// Attempts to denest a nested square root of the form `sqrt(A ± B*sqrt(C))`.
///
/// This function applies the denesting formula: `sqrt(X ± sqrt(Y)) = sqrt((X+sqrt(X^2-Y))/2) ± sqrt((X-sqrt(X^2-Y))/2)`.
/// It looks for patterns `A + B*sqrt(C)` and `A - B*sqrt(C)` inside the outermost square root.
///
/// # Arguments
/// * `expr` - The expression containing the nested square root.
///
/// # Returns
/// The simplified expression if denesting is successful, or the original expression if no simplification is found.
#[must_use]

pub fn denest_sqrt(
    expr : &Expr
) -> Expr {

    let expr_resolved =
        if let Expr::Dag(node) = expr {

            node.to_expr()
                .unwrap_or_else(|_| {
                    expr.clone()
                })
        } else {

            expr.clone()
        };

    if let Expr::Sqrt(inner) =
        &expr_resolved
    {

        // Handle A + B*sqrt(C)
        if let Some((a, b, c)) =
            match_nested_sqrt_pattern(
                inner,
            )
        {

            if let Some(res) =
                apply_denesting(
                    a, b, c, true,
                )
            {

                return res;
            }
        }

        // Handle A - B*sqrt(C)
        if let Some((a, b, c)) = match_nested_sqrt_sub_pattern(inner) {

            if let Some(res) = apply_denesting(a, b, c, false) {

                return res;
            }
        }

        // Handle simple A + sqrt(C) where B=1
        if let Expr::Add(a, term_b) =
            inner.as_ref()
        {

            if let Expr::Sqrt(c) =
                term_b.as_ref()
            {

                if let Some(res) = apply_denesting(
                    a.as_ref().clone(),
                    Expr::new_constant(1.0),
                    c.as_ref().clone(),
                    true,
                ) {

                    return res;
                }
            }
        }

        // Handle simple A - sqrt(C) where B=1
        if let Expr::Sub(a, term_b) =
            inner.as_ref()
        {

            if let Expr::Sqrt(c) =
                term_b.as_ref()
            {

                if let Some(res) = apply_denesting(
                    a.as_ref().clone(),
                    Expr::new_constant(1.0),
                    c.as_ref().clone(),
                    false,
                ) {

                    return res;
                }
            }
        }
    }

    expr.clone()
}

fn apply_denesting(
    a : Expr,
    b : Expr,
    c : Expr,
    is_add : bool,
) -> Option<Expr> {

    // We have sqrt(A ± B*sqrt(C))
    // This is equivalent to sqrt(A ± sqrt(B^2 * C))
    // Let X = A, Y = B^2 * C
    // Formula: sqrt((X + sqrt(X^2 - Y))/2) ± sqrt((X - sqrt(X^2 - Y))/2)

    let x = a;

    let y = simplify(&Expr::new_mul(
        Expr::new_pow(
            b,
            Expr::new_constant(2.0),
        ),
        c,
    ));

    let discriminant =
        simplify(&Expr::new_sub(
            Expr::new_pow(
                x.clone(),
                Expr::new_constant(2.0),
            ),
            y,
        ));

    if let Some(alpha) =
        is_perfect_square(&discriminant)
    {

        let two =
            Expr::new_constant(2.0);

        let term1 =
            simplify(&Expr::new_div(
                Expr::new_add(
                    x.clone(),
                    alpha.clone(),
                ),
                two.clone(),
            ));

        let term2 =
            simplify(&Expr::new_div(
                Expr::new_sub(x, alpha),
                two,
            ));

        let sqrt_term1 =
            Expr::new_sqrt(term1);

        let sqrt_term2 =
            Expr::new_sqrt(term2);

        return Some(simplify(
            &if is_add {

                Expr::new_add(
                    sqrt_term1,
                    sqrt_term2,
                )
            } else {

                Expr::new_sub(
                    sqrt_term1,
                    sqrt_term2,
                )
            },
        ));
    }

    None
}

/// Matches an expression of the form A + B*sqrt(C).

pub(crate) fn match_nested_sqrt_pattern(
    expr : &Expr
) -> Option<(Expr, Expr, Expr)> {

    if let Expr::Add(a, term_b) = expr {

        if let Expr::Mul(b, sqrt_c) =
            &**term_b
        {

            if let Expr::Sqrt(c) =
                &**sqrt_c
            {

                return Some((
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                ));
            }
        }

        // Handle commutative case: B*sqrt(C) + A
        if let Expr::Mul(b, sqrt_c) =
            &**a
        {

            if let Expr::Sqrt(c) =
                &**sqrt_c
            {

                return Some((
                    term_b
                        .as_ref()
                        .clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                ));
            }
        }
    }

    None
}

/// Matches an expression of the form A - B*sqrt(C).

pub(crate) fn match_nested_sqrt_sub_pattern(
    expr : &Expr
) -> Option<(Expr, Expr, Expr)> {

    if let Expr::Sub(a, term_b) = expr {

        if let Expr::Mul(b, sqrt_c) =
            &**term_b
        {

            if let Expr::Sqrt(c) =
                &**sqrt_c
            {

                return Some((
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                ));
            }
        }
    }

    None
}

/// Checks if an expression is a perfect square and returns its root if so.

pub(crate) fn is_perfect_square(
    expr : &Expr
) -> Option<Expr> {

    let expr =
        if let Expr::Dag(node) = expr {

            node.to_expr()
                .unwrap_or_else(|_| {
                    expr.clone()
                })
        } else {

            expr.clone()
        };

    match expr {
        | Expr::Constant(c) => {

            if c >= 0.0 {

                let root = c.sqrt();

                if (root - root.round())
                    .abs()
                    < f64::EPSILON
                {

                    return Some(
                        Expr::Constant(
                            root,
                        ),
                    );
                }
            }
        },
        | Expr::BigInt(n) => {

            use num_traits::ToPrimitive;

            if let Some(f) = n.to_f64()
            {

                if f >= 0.0 {

                    let root = f.sqrt();

                    if (root
                        - root.round())
                    .abs()
                        < f64::EPSILON
                    {

                        return Some(Expr::Constant(root));
                    }
                }
            }
        },
        | _ => {},
    }

    // TODO: Add symbolic perfect square check
    None
}
