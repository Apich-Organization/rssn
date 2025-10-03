use crate::symbolic::core::Expr;
use std::collections::HashMap;

// --- Helper Functions for Factorization ---

/// Breaks a single term (like `2*x^2*y`) into a map of its base factors and their counts.
pub(crate) fn get_term_factors(expr: &Expr) -> HashMap<Expr, i32> {
    let mut factors = HashMap::new();
    match expr {
        Expr::Mul(a, b) => {
            let af = get_term_factors(a);
            let bf = get_term_factors(b);
            for (k, v) in af {
                *factors.entry(k).or_insert(0) += v;
            }
            for (k, v) in bf {
                *factors.entry(k).or_insert(0) += v;
            }
        }
        Expr::Power(base, exp) => {
            if let Some(n) = exp.to_f64() {
                if n.fract() == 0.0 {
                    *factors.entry(*base.clone()).or_insert(0) += n as i32;
                }
            } else {
                // Treat non-constant powers as a single factor
                factors.insert(expr.clone(), 1);
            }
        }
        Expr::Neg(a) => {
            factors.insert(Expr::Constant(-1.0), 1);
            let af = get_term_factors(a);
            for (k, v) in af {
                *factors.entry(k).or_insert(0) += v;
            }
        }
        _ => {
            factors.insert(expr.clone(), 1);
        }
    }
    factors
}

/// Reconstructs an expression from a map of factors and their counts.
pub(crate) fn build_expr_from_factors(factors: HashMap<Expr, i32>) -> Expr {
    if factors.is_empty() {
        return Expr::Constant(1.0);
    }

    let mut terms: Vec<Expr> = factors
        .into_iter()
        .filter(|(_, count)| *count != 0) // Filter out factors with a zero count
        .map(|(base, count)| {
            if count == 1 {
                base
            } else {
                Expr::Power(Box::new(base), Box::new(Expr::Constant(count as f64)))
            }
        })
        .collect();

    if terms.is_empty() {
        return Expr::Constant(1.0);
    }

    terms.sort_unstable();

    let mut tree = terms.remove(0);
    for term in terms {
        tree = Expr::Mul(Box::new(tree), Box::new(term));
    }
    tree
}

// --- Main CAS Functions ---

/// Flattens a nested chain of `Add` expressions into a vector of terms.
pub(crate) fn flatten_sum(expr: Expr, terms: &mut Vec<Expr>) {
    match expr {
        Expr::Add(a, b) => {
            flatten_sum(*a, terms);
            flatten_sum(*b, terms);
        }
        _ => terms.push(expr),
    }
}

/// Flattens a nested chain of `Mul` expressions into two vectors: numeric and other factors.
pub(crate) fn flatten_product(
    expr: Expr,
    numeric_factors: &mut Vec<f64>,
    other_factors: &mut Vec<Expr>,
) {
    match expr {
        Expr::Mul(a, b) => {
            flatten_product(*a, numeric_factors, other_factors);
            flatten_product(*b, numeric_factors, other_factors);
        }
        Expr::Constant(n) => numeric_factors.push(n),
        _ => other_factors.push(expr),
    }
}

/// Normalizes an expression to a canonical form.
pub fn normalize(expr: Expr) -> Expr {
    match expr {
        Expr::Add(a, b) => {
            let mut terms = Vec::new();
            flatten_sum(
                Expr::Add(Box::new(normalize(*a)), Box::new(normalize(*b))),
                &mut terms,
            );
            if terms.len() == 1 {
                return terms.pop().unwrap();
            }
            terms.sort_unstable();
            build_sum_from_vec(terms)
        }
        Expr::Mul(a, b) => {
            let mut numeric_factors = Vec::new();
            let mut other_factors = Vec::new();
            flatten_product(
                Expr::Mul(Box::new(normalize(*a)), Box::new(normalize(*b))),
                &mut numeric_factors,
                &mut other_factors,
            );
            if numeric_factors.is_empty() && other_factors.len() == 1 {
                return other_factors.pop().unwrap();
            }
            other_factors.sort_unstable();
            build_product_from_vecs(numeric_factors, other_factors)
        }
        Expr::Sub(a, b) => Expr::Sub(Box::new(normalize(*a)), Box::new(normalize(*b))),
        Expr::Div(a, b) => Expr::Div(Box::new(normalize(*a)), Box::new(normalize(*b))),
        Expr::Power(a, b) => Expr::Power(Box::new(normalize(*a)), Box::new(normalize(*b))),
        Expr::Neg(a) => Expr::Neg(Box::new(normalize(*a))),
        Expr::Sin(a) => Expr::Sin(Box::new(normalize(*a))),
        Expr::Cos(a) => Expr::Cos(Box::new(normalize(*a))),
        Expr::Tan(a) => Expr::Tan(Box::new(normalize(*a))),
        Expr::Exp(a) => Expr::Exp(Box::new(normalize(*a))),
        Expr::Log(a) => Expr::Log(Box::new(normalize(*a))),
        Expr::Vector(v) => Expr::Vector(v.into_iter().map(normalize).collect()),
        Expr::Matrix(m) => Expr::Matrix(
            m.into_iter()
                .map(|row| row.into_iter().map(normalize).collect())
                .collect(),
        ),
        e => e,
    }
}

/// Expands expressions by applying the distributive property and expanding powers.
pub fn expand(expr: Expr) -> Expr {
    let expanded_expr = match expr {
        Expr::Mul(a, b) => {
            let exp_a = expand(*a);
            let exp_b = expand(*b);
            match (exp_a, exp_b) {
                (Expr::Add(a1, a2), b_expr) => {
                    let term1 = expand(Expr::Mul(a1, Box::new(b_expr.clone())));
                    let term2 = expand(Expr::Mul(a2, Box::new(b_expr)));
                    Expr::Add(Box::new(term1), Box::new(term2))
                }
                (a_expr, Expr::Add(b1, b2)) => {
                    let term1 = expand(Expr::Mul(Box::new(a_expr.clone()), b1));
                    let term2 = expand(Expr::Mul(Box::new(a_expr), b2));
                    Expr::Add(Box::new(term1), Box::new(term2))
                }
                (exp_a, exp_b) => Expr::Mul(Box::new(exp_a), Box::new(exp_b)),
            }
        }
        Expr::Power(base, exp) => {
            let exp_base = expand(*base);
            let exp_exp = expand(*exp);
            if let Some(n) = exp_exp.to_f64() {
                if n.fract() == 0.0 && n > 1.0 {
                    let n_us = n as usize;
                    let mut result = exp_base.clone();
                    for _ in 1..n_us {
                        result = Expr::Mul(Box::new(result), Box::new(exp_base.clone()));
                    }
                    return expand(result);
                }
            }
            Expr::Power(Box::new(exp_base), Box::new(exp_exp))
        }
        Expr::Sub(a, b) => Expr::Sub(Box::new(expand(*a)), Box::new(expand(*b))),
        Expr::Neg(a) => match expand(*a) {
            Expr::Add(b, c) => Expr::Add(Box::new(Expr::Neg(b)), Box::new(Expr::Neg(c))),
            Expr::Neg(b) => *b,
            expanded_a => Expr::Neg(Box::new(expanded_a)),
        },
        Expr::Div(a, b) => Expr::Div(Box::new(expand(*a)), Box::new(expand(*b))),
        Expr::Sin(a) => Expr::Sin(Box::new(expand(*a))),
        Expr::Cos(a) => Expr::Cos(Box::new(expand(*a))),
        Expr::Tan(a) => Expr::Tan(Box::new(expand(*a))),
        Expr::Exp(a) => Expr::Exp(Box::new(expand(*a))),
        Expr::Log(a) => Expr::Log(Box::new(expand(*a))),
        Expr::Vector(v) => Expr::Vector(v.into_iter().map(expand).collect()),
        Expr::Matrix(m) => Expr::Matrix(
            m.into_iter()
                .map(|row| row.into_iter().map(expand).collect())
                .collect(),
        ),
        e => e,
    };
    normalize(expanded_expr)
}

/// Factorizes an expression by extracting common factors from sums.
pub fn factorize(expr: Expr) -> Expr {
    let expanded = expand(expr);
    match expanded {
        Expr::Add(a, b) => {
            let mut terms = Vec::new();
            flatten_sum(Expr::Add(a, b), &mut terms);
            let term_factors: Vec<HashMap<Expr, i32>> =
                terms.iter().map(get_term_factors).collect();
            let mut gcd_factors = term_factors.first().cloned().unwrap_or_default();
            for next_term_map in term_factors.iter().skip(1) {
                let mut next_gcd = HashMap::new();
                for (base, count) in &gcd_factors {
                    if let Some(next_count) = next_term_map.get(base) {
                        next_gcd.insert(base.clone(), (*count).min(*next_count));
                    }
                }
                gcd_factors = next_gcd;
            }
            if gcd_factors.is_empty()
                || (gcd_factors.len() == 1
                    && gcd_factors.keys().next() == Some(&Expr::Constant(1.0)))
            {
                return build_sum_from_vec(terms);
            }
            let gcd_expr = build_expr_from_factors(gcd_factors.clone());
            let mut new_terms = Vec::new();
            for term_map in &term_factors {
                let mut remaining_factors = term_map.clone();
                for (base, gcd_count) in &gcd_factors {
                    if let Some(term_count) = remaining_factors.get_mut(base) {
                        *term_count -= gcd_count;
                    }
                }
                new_terms.push(build_expr_from_factors(remaining_factors));
            }
            let remaining_sum = build_sum_from_vec(new_terms);
            Expr::Mul(Box::new(gcd_expr), Box::new(remaining_sum))
        }
        Expr::Mul(a, b) => Expr::Mul(Box::new(factorize(*a)), Box::new(factorize(*b))),
        Expr::Power(a, b) => Expr::Power(Box::new(factorize(*a)), Box::new(factorize(*b))),
        Expr::Neg(a) => Expr::Neg(Box::new(factorize(*a))),
        e => e,
    }
}

/// Helper to build a normalized sum from a vector of expressions.
pub(crate) fn build_sum_from_vec(mut terms: Vec<Expr>) -> Expr {
    if terms.is_empty() {
        return Expr::Constant(0.0);
    }
    if terms.len() == 1 {
        return terms.pop().unwrap();
    }
    terms.sort_unstable();
    let mut tree = terms.remove(0);
    for term in terms {
        tree = Expr::Add(Box::new(tree), Box::new(term));
    }
    tree
}

/// Helper to build a normalized product from vectors of numeric and other factors.
pub(crate) fn build_product_from_vecs(numeric_factors: Vec<f64>, other_factors: Vec<Expr>) -> Expr {
    let numeric_product: f64 = numeric_factors.iter().product();
    let has_numeric_term = numeric_product != 1.0 || other_factors.is_empty();

    let mut tree: Option<Expr> = None;

    if has_numeric_term {
        tree = Some(Expr::Constant(numeric_product));
    }

    for factor in other_factors {
        if let Some(t) = tree {
            tree = Some(Expr::Mul(Box::new(t), Box::new(factor)));
        } else {
            tree = Some(factor);
        }
    }
    tree.unwrap_or(Expr::Constant(1.0))
}

/// Placeholder for Risch algorithm for symbolic integration.
pub fn risch_integrate(expr: Expr, var: &str) -> Expr {
    Expr::Variable(format!("RischIntegrate({}, {})", expr, var))
}

/// Placeholder for Gröbner Basis computation for solving polynomial systems.
pub fn grobner_basis(_polynomials: Vec<Expr>, _variables: Vec<String>) -> Vec<Expr> {
    vec![Expr::Variable("GröbnerBasis(system)".to_string())]
}

/// Placeholder for Cylindrical Algebraic Decomposition (CAD) for real algebraic geometry.
pub fn cylindrical_algebraic_decomposition(
    _polynomials: Vec<Expr>,
    _variables: Vec<String>,
) -> Expr {
    Expr::Variable("CylindricalAlgebraicDecomposition(system)".to_string())
}
