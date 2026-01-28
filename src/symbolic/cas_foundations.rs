use std::collections::HashMap;

use crate::symbolic::core::Expr;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::MonomialOrder;
use crate::symbolic::grobner::{
    self,
};
use crate::symbolic::polynomial::expr_to_sparse_poly;
use crate::symbolic::polynomial::sparse_poly_to_expr;

const ERROR_MARGIN: f64 = 1e-9;

/// Breaks a single term (like `2*x^2*y`) into a map of its base factors and their counts.
///
/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.
#[must_use]

pub fn get_term_factors(
    expr: &Expr
) -> HashMap<Expr, i32> {

    let mut factors = HashMap::new();

    match expr {
        | Expr::Dag(node) => {

            return get_term_factors(
                &node
                    .to_expr()
                    .expect(
                    "Get term factors",
                ),
            );
        },
        | Expr::Mul(a, b) => {

            let af =
                get_term_factors(a);

            let bf =
                get_term_factors(b);

            for (k, v) in af {

                *factors
                    .entry(k)
                    .or_insert(0) += v;
            }

            for (k, v) in bf {

                *factors
                    .entry(k)
                    .or_insert(0) += v;
            }
        },
        | Expr::Power(base, exp) => {
            if let Some(n) =
                exp.to_f64()
            {

                if n.fract() == 0.0 {

                    *factors
                        .entry(
                            base.as_ref()
                                .clone(),
                        )
                        .or_insert(0) += n as i32;
                }
            } else {

                factors.insert(
                    expr.clone(),
                    1,
                );
            }
        },
        | Expr::Neg(a) => {

            factors.insert(
                Expr::new_constant(-1.0),
                1,
            );

            let af =
                get_term_factors(a);

            for (k, v) in af {

                *factors
                    .entry(k)
                    .or_insert(0) += v;
            }
        },
        | _ => {

            factors.insert(
                expr.clone(),
                1,
            );
        },
    }

    factors
}

/// Reconstructs an expression from a map of factors and their counts.

pub fn build_expr_from_factors<
    S: ::std::hash::BuildHasher,
>(
    factors: HashMap<Expr, i32, S>
) -> Expr {

    if factors.is_empty() {

        return Expr::new_constant(1.0);
    }

    let mut terms: Vec<Expr> = factors
        .into_iter()
        .filter(|(_, count)| {

            *count != 0
        })
        .map(|(base, count)| {
            if count == 1 {

                base
            } else {

                Expr::new_pow(
                    base,
                    Expr::new_constant(
                        f64::from(
                            count,
                        ),
                    ),
                )
            }
        })
        .collect();

    if terms.is_empty() {

        return Expr::new_constant(1.0);
    }

    terms.sort();

    let mut tree = terms.remove(0);

    for term in terms {

        tree =
            Expr::new_mul(tree, term);
    }

    tree
}

/// Flattens a nested chain of `Add` expressions into a vector of terms.

/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.

pub(crate) fn flatten_sum(
    expr: Expr,
    terms: &mut Vec<Expr>,
) {

    match expr {
        | Expr::Dag(node) => {

            flatten_sum(
                node.to_expr()
                    .expect(
                        "Flatten Sum",
                    ),
                terms,
            );
        },
        | Expr::Add(a, b) => {

            flatten_sum(
                (*a).clone(),
                terms,
            );

            flatten_sum(
                (*b).clone(),
                terms,
            );
        },
        | _ => terms.push(expr),
    }
}

/// Flattens a nested chain of `Mul` expressions into two vectors: numeric and other factors.

/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.

pub(crate) fn flatten_product(
    expr: Expr,
    numeric_factors: &mut Vec<f64>,
    other_factors: &mut Vec<Expr>,
) {

    match expr {
        | Expr::Dag(node) => {

            flatten_product(
                node.to_expr()
                    .expect(
                    "Flatten product",
                ),
                numeric_factors,
                other_factors,
            );
        },
        | Expr::Mul(a, b) => {

            flatten_product(
                (*a).clone(),
                numeric_factors,
                other_factors,
            );

            flatten_product(
                (*b).clone(),
                numeric_factors,
                other_factors,
            );
        },
        | Expr::new_constant(n) => {

            numeric_factors.push(n);
        },
        | _ => other_factors.push(expr),
    }
}

/// Normalizes an expression to a canonical form.
///
/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.

pub fn normalize(expr: Expr) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            normalize(
                node.to_expr()
                    .expect("Normalize"),
            )
        },
        | Expr::Add(a, b) => {

            let mut terms = Vec::new();

            flatten_sum(
                Expr::new_add(
                    normalize((*a).clone()),
                    normalize((*b).clone()),
                ),
                &mut terms,
            );

            if terms.len() == 1 {

                return match terms.pop() {
                    | Some(t) => t,
                    | _none => unreachable!(),
                };
            }

            terms.sort();

            build_sum_from_vec(terms)
        },
        | Expr::Mul(a, b) => {

            let mut numeric_factors = Vec::new();

            let mut other_factors = Vec::new();

            flatten_product(
                Expr::new_mul(
                    normalize((*a).clone()),
                    normalize((*b).clone()),
                ),
                &mut numeric_factors,
                &mut other_factors,
            );

            if numeric_factors.is_empty() && other_factors.len() == 1 {

                return match other_factors.pop() {
                    | Some(f) => f,
                    | _none => unreachable!(),
                };
            }

            other_factors.sort();

            build_product_from_vecs(
                &numeric_factors,
                other_factors,
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                normalize((*a).clone()),
                normalize((*b).clone()),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                normalize((*a).clone()),
                normalize((*b).clone()),
            )
        },
        | Expr::Power(a, b) => {
            Expr::new_pow(
                normalize((*a).clone()),
                normalize((*b).clone()),
            )
        },
        | Expr::Neg(a) => {
            Expr::new_neg(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Sin(a) => {
            Expr::new_sin(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Cos(a) => {
            Expr::new_cos(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Tan(a) => {
            Expr::new_tan(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Exp(a) => {
            Expr::new_exp(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Log(a) => {
            Expr::new_log(normalize(
                (*a).clone(),
            ))
        },
        | Expr::Vector(v) => {
            Expr::Vector(
                v.into_iter()
                    .map(normalize)
                    .collect(),
            )
        },
        | Expr::Matrix(m) => {
            Expr::Matrix(
                m.into_iter()
                    .map(|row| {

                        row.into_iter()
                            .map(normalize)
                            .collect()
                    })
                    .collect(),
            )
        },
        | e => e,
    }
}

/// Expands expressions by applying the distributive property and expanding powers.
///
/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.

pub fn expand(expr: Expr) -> Expr {

    let expanded_expr = match expr {
        | Expr::Dag(node) => {

            return expand(
                node.to_expr()
                    .expect("Expand"),
            );
        },
        | Expr::Mul(a, b) => {

            let exp_a = expand(
                a.as_ref().clone(),
            );

            let exp_b = expand(
                b.as_ref().clone(),
            );

            let exp_a_unwrapped =
                match exp_a {
                    | Expr::Dag(
                        ref node,
                    ) => {
                        node.to_expr()
                            .unwrap_or_else(|_| exp_a.clone())
                    },
                    | _ => exp_a,
                };

            let exp_b_unwrapped =
                match exp_b {
                    | Expr::Dag(
                        ref node,
                    ) => {
                        node.to_expr()
                            .unwrap_or_else(|_| exp_b.clone())
                    },
                    | _ => exp_b,
                };

            match (
                exp_a_unwrapped,
                exp_b_unwrapped,
            ) {
                | (
                    Expr::Add(a1, a2),
                    b_expr,
                ) => {

                    let t1 =
                        Expr::new_mul(
                            a1.as_ref()
                                .clone(
                                ),
                            b_expr
                                .clone(
                                ),
                        );

                    let t2 =
                        Expr::new_mul(
                            a2.as_ref()
                                .clone(
                                ),
                            b_expr,
                        );

                    expand(
                        Expr::new_add(
                            t1, t2,
                        ),
                    )
                },
                | (
                    Expr::AddList(
                        terms,
                    ),
                    ref b_expr,
                ) => {

                    let new_terms : Vec<Expr> = terms
                        .iter()
                        .map(|t| {

                            Expr::new_mul(
                                t.clone(),
                                b_expr.clone(),
                            )
                        })
                        .collect();

                    expand(build_sum_from_vec(
                        new_terms,
                    ))
                },
                | (
                    a_expr,
                    Expr::Add(b1, b2),
                ) => {

                    let t1 =
                        Expr::new_mul(
                            a_expr
                                .clone(
                                ),
                            b1.as_ref()
                                .clone(
                                ),
                        );

                    let t2 =
                        Expr::new_mul(
                            a_expr,
                            b2.as_ref()
                                .clone(
                                ),
                        );

                    expand(
                        Expr::new_add(
                            t1, t2,
                        ),
                    )
                },
                | (
                    ref a_expr,
                    Expr::AddList(
                        terms,
                    ),
                ) => {

                    let new_terms : Vec<Expr> = terms
                        .iter()
                        .map(|t| {

                            Expr::new_mul(
                                a_expr.clone(),
                                t.clone(),
                            )
                        })
                        .collect();

                    expand(build_sum_from_vec(
                        new_terms,
                    ))
                },
                | (
                    exp_a_u,
                    exp_b_u,
                ) => {
                    Expr::new_mul(
                        exp_a_u,
                        exp_b_u,
                    )
                },
            }
        },
        | Expr::Power(base, exp) => {

            let exp_base = expand(
                base.as_ref()
                    .clone(),
            );

            let exp_exp = expand(
                exp.as_ref().clone(),
            );

            if let Some(n) =
                exp_exp.to_f64()
            {

                if n.fract() == 0.0
                    && n > 1.0
                {

                    let n_us = (n
                        as i64)
                        .try_into()
                        .unwrap_or(0);

                    let mut result =
                        exp_base
                            .clone();

                    for _ in 1 .. n_us {

                        result = Expr::new_mul(
                            result,
                            exp_base.clone(),
                        );
                    }

                    return expand(
                        result,
                    );
                }
            }

            Expr::new_pow(
                exp_base,
                exp_exp,
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                expand(
                    a.as_ref().clone(),
                ),
                expand(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Neg(a) => {

            match expand(
                a.as_ref().clone(),
            ) {
                | Expr::Add(b, c) => {
                    Expr::new_add(
                        Expr::new_neg(
                            b,
                        ),
                        Expr::new_neg(
                            c,
                        ),
                    )
                },
                | Expr::Neg(b) => {
                    b.as_ref().clone()
                },
                | expanded_a => {
                    Expr::new_neg(
                        expanded_a,
                    )
                },
            }
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                expand(
                    a.as_ref().clone(),
                ),
                expand(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Sin(a) => {
            Expr::new_sin(expand(
                a.as_ref().clone(),
            ))
        },
        | Expr::Cos(a) => {
            Expr::new_cos(expand(
                a.as_ref().clone(),
            ))
        },
        | Expr::Tan(a) => {
            Expr::new_tan(expand(
                a.as_ref().clone(),
            ))
        },
        | Expr::Exp(a) => {
            Expr::new_exp(expand(
                a.as_ref().clone(),
            ))
        },
        | Expr::Log(a) => {
            Expr::new_log(expand(
                a.as_ref().clone(),
            ))
        },
        | Expr::Vector(v) => {
            Expr::Vector(
                v.into_iter()
                    .map(expand)
                    .collect(),
            )
        },
        | Expr::Matrix(m) => {
            Expr::Matrix(
                m.into_iter()
                    .map(|row| {

                        row.into_iter()
                            .map(expand)
                            .collect()
                    })
                    .collect(),
            )
        },
        | _ => expr,
    };

    normalize(expanded_expr)
}

/// Factorizes an expression by extracting common factors from sums.
#[must_use]

pub fn factorize(expr: Expr) -> Expr {

    let expanded = expand(expr);

    let expanded_unwrapped =
        match &expanded {
            | Expr::Dag(node) => {
                node.to_expr()
                    .unwrap_or_else(
                        |_| {

                            expanded
                                .clone()
                        },
                    )
            },
            | _ => expanded.clone(),
        };

    match expanded_unwrapped {
        | Expr::Add(a, b) => {

            let mut terms = Vec::new();

            flatten_sum(
                Expr::new_add(a, b),
                &mut terms,
            );

            factorize_terms(terms)
        },
        | Expr::AddList(terms) => {

            let mut flat_terms =
                Vec::new();

            for term in terms {

                flatten_sum(
                    term,
                    &mut flat_terms,
                );
            }

            factorize_terms(flat_terms)
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                factorize(
                    a.as_ref().clone(),
                ),
                factorize(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::MulList(factors) => {

            let new_factors: Vec<Expr> =
                factors
                    .iter()
                    .map(|f| {

                        factorize(
                            f.clone(),
                        )
                    })
                    .collect();

            // Reconstruct MulList or Mul
            // For simplicity, just fold
            if new_factors.is_empty() {

                Expr::new_constant(1.0)
            } else {

                let mut res =
                    new_factors[0]
                        .clone();

                for f in new_factors
                    .iter()
                    .skip(1)
                {

                    res = Expr::new_mul(
                        res,
                        f.clone(),
                    );
                }

                res
            }
        },
        | Expr::Power(a, b) => {
            Expr::new_pow(
                factorize(
                    a.as_ref().clone(),
                ),
                factorize(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Neg(a) => {
            Expr::new_neg(factorize(
                a.as_ref().clone(),
            ))
        },
        | e => e,
    }
}

pub(crate) fn factorize_terms(
    terms: Vec<Expr>
) -> Expr {

    let term_factors: Vec<
        HashMap<Expr, i32>,
    > = terms
        .iter()
        .map(get_term_factors)
        .collect();

    let mut gcd_factors = term_factors
        .first()
        .cloned()
        .unwrap_or_default();

    for next_term_map in term_factors
        .iter()
        .skip(1)
    {

        let mut next_gcd =
            HashMap::new();

        for (base, count) in
            &gcd_factors
        {

            if let Some(next_count) =
                next_term_map.get(base)
            {

                next_gcd.insert(
                    base.clone(),
                    (*count).min(
                        *next_count,
                    ),
                );
            }
        }

        gcd_factors = next_gcd;
    }

    if gcd_factors.is_empty()
        || (gcd_factors.len() == 1
            && gcd_factors
                .keys()
                .next()
                == Some(
                    &Expr::new_constant(
                        1.0,
                    ),
                ))
    {

        return build_sum_from_vec(
            terms,
        );
    }

    let gcd_expr =
        build_expr_from_factors(
            gcd_factors.clone(),
        );

    let mut new_terms = Vec::new();

    for term_map in &term_factors {

        let mut remaining_factors =
            term_map.clone();

        for (base, gcd_count) in
            &gcd_factors
        {

            if let Some(term_count) =
                remaining_factors
                    .get_mut(base)
            {

                *term_count -=
                    gcd_count;
            }
        }

        new_terms.push(
            build_expr_from_factors(
                remaining_factors,
            ),
        );
    }

    let remaining_sum =
        build_sum_from_vec(new_terms);

    Expr::new_mul(
        gcd_expr,
        remaining_sum,
    )
}

/// Helper to build a normalized sum from a vector of expressions.

pub(crate) fn build_sum_from_vec(
    mut terms: Vec<Expr>
) -> Expr {

    if terms.is_empty() {

        return Expr::new_constant(0.0);
    }

    if terms.len() == 1 {

        return match terms.pop() {
            | Some(t) => t,
            | _none => unreachable!(),
        };
    }

    terms.sort();

    let mut tree = terms.remove(0);

    for term in terms {

        tree =
            Expr::new_add(tree, term);
    }

    tree
}

/// Helper to build a normalized product from vectors of numeric and other factors.

pub(crate) fn build_product_from_vecs(
    numeric_factors: &[f64],
    other_factors: Vec<Expr>,
) -> Expr {

    let numeric_product: f64 =
        numeric_factors
            .iter()
            .product();

    let has_numeric_term =
        (numeric_product - 1.0).abs()
            > ERROR_MARGIN
            || other_factors.is_empty();

    let mut tree: Option<Expr> =
        if has_numeric_term {

            Some(Expr::new_constant(
                numeric_product,
            ))
        } else {

            None
        };

    for factor in other_factors {

        tree = Some(
            if let Some(t) = tree {

                Expr::new_mul(t, factor)
            } else {

                factor
            },
        );
    }

    tree.unwrap_or(Expr::new_constant(1.0))
}

/// Placeholder for Risch algorithm for symbolic integration.
///
/// **Note:** This function is deprecated.
#[deprecated(
    since = "0.1.9",
    note = "Please use \
            symbolic/integrate \
            instead."
)]
#[must_use]

pub fn risch_integrate(
    expr: &Expr,
    var: &str,
) -> Expr {

    Expr::Variable(format!(
        "RischIntegrate({expr}, {var})"
    ))
}

/// Placeholder for Gröbner Basis computation for solving polynomial systems.
///
/// **Note:** This function is deprecated.
#[deprecated(
    since = "0.1.9",
    note = "Please use \
            symbolic/grobner instead."
)]
#[must_use]

pub fn grobner_basis(
    _polynomials: Vec<Expr>,
    _variables: Vec<String>,
) -> Vec<Expr> {

    vec![Expr::Variable(
        "GröbnerBasis(system)"
            .to_string(),
    )]
}

/// Placeholder for Cylindrical Algebraic Decomposition (CAD) for real algebraic geometry.
///
/// **Note:** This function is deprecated.
#[deprecated(
    since = "0.1.9",
    note = "Please use symbolic/cad \
            instead."
)]
#[must_use]

pub fn cylindrical_algebraic_decomposition(
    _polynomials: Vec<Expr>,
    _variables: Vec<String>,
) -> Expr {

    Expr::Variable("CylindricalAlgebraicDecomposition(system)".to_string())
}

/// Simplifies an expression using a set of polynomial side-relations.
///
/// This function provides a powerful simplification mechanism by calculating the normal form
/// of the expression with respect to a polynomial ideal defined by the given relations.
/// It computes the Gröbner basis of the relations and uses it to reduce the expression.
///
/// # Arguments
/// * `expr` - The expression to simplify.
/// * `relations` - A slice of expressions representing the polynomial relations (e.g., `x^2 + y^2 - 1`).
/// * `vars` - A slice of variable names to establish an ordering for the polynomial ring.
/// * `order` - The monomial ordering to use for the Gröbner basis computation.
///
/// # Returns
/// A new, simplified `Expr`. Returns the original expression if any step fails.
#[must_use]

pub fn simplify_with_relations(
    expr: &Expr,
    relations: &[Expr],
    vars: &[&str],
    order: MonomialOrder,
) -> Expr {

    // 1. Convert the expression and relations to sparse polynomials.
    let p =
        expr_to_sparse_poly(expr, vars);

    let relation_polys: Vec<
        SparsePolynomial,
    > = relations
        .iter()
        .map(|r| {

            expr_to_sparse_poly(r, vars)
        })
        .collect();

    // 2. Compute the Gröbner basis of the relations.
    let g_basis =
        match grobner::buchberger(
            &relation_polys,
            order,
        ) {
            | Ok(basis) => basis,
            | Err(_) => {

                return expr.clone();
            }, /* Return original on failure */
        };

    // 3. Compute the normal form (remainder) of the expression w.r.t. the basis.
    match grobner::poly_division_multivariate(&p, &g_basis, order) {
        | Ok((_, remainder)) => {

            // 4. Convert the remainder polynomial back to an expression.
            sparse_poly_to_expr(&remainder)
        },
        | Err(_) => expr.clone(), // Return original on failure
    }
}

/// Normalizes an expression to a canonical form using a set of polynomial side-relations.
///
/// This is an alias for `simplify_with_relations`, as the normal form with respect to a
/// Gröbner basis is a canonical representation for polynomials within the ideal.
#[must_use]

pub fn normalize_with_relations(
    expr: &Expr,
    relations: &[Expr],
    vars: &[&str],
    order: MonomialOrder,
) -> Expr {

    simplify_with_relations(
        expr,
        relations,
        vars,
        order,
    )
}
