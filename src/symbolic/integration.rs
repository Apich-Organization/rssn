//! # Advanced Symbolic Integration Techniques
//!
//! This module provides advanced symbolic integration techniques, particularly focusing
//! on the Risch-Norman algorithm for integrating elementary functions. It includes
//! implementations for integrating rational functions (Hermite-Ostrogradsky method)
//! and handling transcendental extensions (logarithmic and exponential cases).

use crate::symbolic::calculus::{
    differentiate,
    integrate,
    substitute,
};
use crate::symbolic::core::{
    Expr,
    Monomial,
    SparsePolynomial,
};
use crate::symbolic::matrix::determinant;
use crate::symbolic::number_theory::expr_to_sparse_poly;
use crate::symbolic::polynomial::gcd;
use crate::symbolic::polynomial::poly_mul_scalar_expr;
use crate::symbolic::polynomial::{
    contains_var,
    differentiate_poly,
    sparse_poly_to_expr,
};
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::solve::solve;
use crate::symbolic::solve::solve_system;
use std::collections::{
    BTreeMap,
    HashMap,
};
use std::sync::Arc;

/// Integrates a rational function `P(x)/Q(x)` using the Hermite-Ostrogradsky method.
///
/// This method decomposes the integral of a rational function into a rational part
/// and a transcendental (logarithmic) part. It involves polynomial long division,
/// square-free factorization of the denominator, and solving a system of linear equations.
///
/// # Arguments
/// * `p` - The numerator polynomial as a `SparsePolynomial`.
/// * `q` - The denominator polynomial as a `SparsePolynomial`.
/// * `x` - The variable of integration.
///
/// # Returns
/// A `Result` containing an `Expr` representing the integral, or an error string if computation fails.

pub fn integrate_rational_function(
    p: &SparsePolynomial,
    q: &SparsePolynomial,
    x: &str,
) -> Result<Expr, String> {

    let (quotient, remainder) = p
        .clone()
        .long_division(q.clone(), x);

    let integral_of_quotient =
        poly_integrate(&quotient, x);

    if remainder
        .terms
        .is_empty()
    {

        return Ok(
            integral_of_quotient,
        );
    }

    let q_prime =
        differentiate_poly(q, x);

    let d = gcd(
        q.clone(),
        q_prime.clone(),
        x,
    );

    let b = q
        .clone()
        .long_division(d.clone(), x)
        .0;

    let (a_poly, c_poly) =
        build_and_solve_hermite_system(
            &remainder, &b, &d,
            &q_prime, x,
        )?;

    let rational_part = Expr::new_div(
        sparse_poly_to_expr(&c_poly),
        sparse_poly_to_expr(&d),
    );

    let integral_of_transcendental_part = integrate_square_free_rational_part(&a_poly, &b, x)?;

    Ok(simplify(
        &Expr::new_add(
            integral_of_quotient,
            Expr::new_add(
                rational_part,
                integral_of_transcendental_part,
            ),
        ),
    ))
}

/// Constructs and solves the linear system for coefficients in Hermite integration.

pub(crate) fn build_and_solve_hermite_system(
    p: &SparsePolynomial,
    b: &SparsePolynomial,
    d: &SparsePolynomial,
    q_prime: &SparsePolynomial,
    x: &str,
) -> Result<
    (
        SparsePolynomial,
        SparsePolynomial,
    ),
    String,
> {

    let deg_d = d.degree(x) as usize;

    let deg_b = b.degree(x) as usize;

    let a_coeffs: Vec<_> = (0..deg_b)
        .map(|i| {
            Expr::Variable(format!(
                "a{i}"
            ))
        })
        .collect();

    let c_coeffs: Vec<_> = (0..deg_d)
        .map(|i| {
            Expr::Variable(format!(
                "c{i}"
            ))
        })
        .collect();

    let a_sym =
        poly_from_coeffs(&a_coeffs, x);

    let c_sym =
        poly_from_coeffs(&c_coeffs, x);

    let c_prime_sym =
        differentiate_poly(&c_sym, x);

    let t = (b.clone()
        * q_prime.clone())
    .long_division(d.clone(), x)
    .0;

    let term1 = b.clone() * c_prime_sym;

    let term2 = t * c_sym;

    let term3 = d.clone() * a_sym;

    let rhs_poly =
        (term1 - term2) + term3;

    let mut equations = Vec::new();

    let num_unknowns = deg_b + deg_d;

    for i in 0..=num_unknowns {

        let p_coeff = p
            .get_coeff_for_power(x, i)
            .unwrap_or_else(|| {
                Expr::Constant(0.0)
            });

        let rhs_coeff = rhs_poly
            .get_coeff_for_power(x, i)
            .unwrap_or_else(|| {
                Expr::Constant(0.0)
            });

        equations.push(simplify(
            &Expr::Eq(
                Arc::new(p_coeff),
                Arc::new(rhs_coeff),
            ),
        ));
    }

    let mut unknown_vars_str: Vec<String> = a_coeffs
        .iter()
        .map(std::string::ToString::to_string)
        .collect();

    unknown_vars_str.extend(
        c_coeffs
            .iter()
            .map(std::string::ToString::to_string),
    );

    let unknown_vars: Vec<&str> = unknown_vars_str
        .iter()
        .map(std::string::String::as_str)
        .collect();

    let solutions = solve_system(
        &equations,
        &unknown_vars,
    )
    .ok_or(
        "Failed to solve linear \
         system for coefficients.",
    )?;

    let sol_map: HashMap<_, _> =
        solutions
            .into_iter()
            .collect();

    let final_a_coeffs: Result<
        Vec<Expr>,
        _,
    > = a_coeffs
        .iter()
        .map(|v| {

            sol_map
                .get(v)
                .cloned()
                .ok_or_else(|| {
                    format!(
                        "Solver did \
                         not return a \
                         solution for \
                         coefficient \
                         {v}"
                    )
                })
        })
        .collect();

    let final_a_coeffs =
        final_a_coeffs?;

    let final_c_coeffs: Result<
        Vec<Expr>,
        _,
    > = c_coeffs
        .iter()
        .map(|v| {

            sol_map
                .get(v)
                .cloned()
                .ok_or_else(|| {
                    format!(
                        "Solver did \
                         not return a \
                         solution for \
                         coefficient \
                         {v}"
                    )
                })
        })
        .collect();

    let final_c_coeffs =
        final_c_coeffs?;

    Ok((
        poly_from_coeffs(
            &final_a_coeffs,
            x,
        ),
        poly_from_coeffs(
            &final_c_coeffs,
            x,
        ),
    ))
}

/// Main entry point for Risch-Norman style integration.
#[must_use]

pub fn risch_norman_integrate(
    expr: &Expr,
    x: &str,
) -> Expr {

    if let Some(t) =
        find_outermost_transcendental(
            expr, x,
        )
    {

        if let Ok((a_t, d_t)) =
            expr_to_rational_poly(
                expr, &t, x,
            )
        {

            let (p_t, r_t) = a_t
                .long_division(
                    d_t.clone(),
                    x,
                );

            let poly_integral = match t
            {
                | Expr::Exp(_) => {
                    integrate_poly_exp(
                        &p_t, &t, x,
                    )
                },
                | Expr::Log(_) => {
                    integrate_poly_log(
                        &p_t, &t, x,
                    )
                },
                | _ => Err(
                    "Unsupported \
                     transcendental \
                     type"
                        .to_string(),
                ),
            };

            let rational_integral =
                if r_t.terms.is_empty()
                {

                    Ok(Expr::Constant(
                        0.0,
                    ))
                } else {

                    hermite_integrate_rational(
                    &r_t,
                    &d_t,
                    &t.to_string(),
                )
                };

            if let (Ok(pi), Ok(ri)) = (
                poly_integral,
                rational_integral,
            ) {

                return simplify(
                    &Expr::new_add(
                        pi, ri,
                    ),
                );
            }
        }
    }

    integrate_rational_function_expr(
        expr, x,
    )
    .unwrap_or_else(|_| {
        integrate(expr, x, None, None)
    })
}

/// Integrates the polynomial part of a transcendental function extension F(t) for the logarithmic case.

pub(crate) fn integrate_poly_log(
    p_in_t: &SparsePolynomial,
    t: &Expr,
    x: &str,
) -> Result<Expr, String> {

    let t_var = "t_var";

    if p_in_t.degree(t_var) < 0 {

        return Ok(Expr::Constant(0.0));
    }

    let n =
        p_in_t.degree(t_var) as usize;

    let p_coeffs =
        p_in_t.get_coeffs_as_vec(t_var);

    let p_n = p_coeffs[0].clone();

    let q_n =
        risch_norman_integrate(&p_n, x);

    if let Expr::Integral { .. } = q_n {

        return Err(
            "Recursive integration of \
             leading coefficient \
             failed."
                .to_string(),
        );
    }

    let t_pow_n = SparsePolynomial {
        terms: BTreeMap::from([(
            Monomial(BTreeMap::from([
                (
                    t_var.to_string(),
                    n as u32,
                ),
            ])),
            Expr::Constant(1.0),
        )]),
    };

    let _q_poly_term =
        poly_mul_scalar_expr(
            &t_pow_n, &q_n,
        );

    // Compute d/dx(q(x) * t^n) = q'(x) * t^n + q(x) * n * t^(n-1) * (dt/dx)
    // We need to differentiate q(x) with respect to x, then multiply by t^n
    // Plus q(x) * n * t^(n-1) * (dt/dx)

    // First term: q'(x) * t^n
    let q_n_deriv =
        differentiate(&q_n, x);

    let term1 = poly_mul_scalar_expr(
        &t_pow_n, &q_n_deriv,
    );

    // Second term: q(x) * n * t^(n-1) * (dt/dx)
    // For t = ln(x), dt/dx = 1/x
    let dt_dx = if let Expr::Log(arg) =
        t
    {

        // dt/dx for ln(f(x)) = f'(x)/f(x)
        let f_prime =
            differentiate(arg, x);

        simplify(&Expr::new_div(
            f_prime,
            (**arg).clone(),
        ))
    } else {

        return Err(
            "Only logarithmic case is \
             currently supported"
                .to_string(),
        );
    };

    let mut deriv = term1;

    if n > 0 {

        let t_pow_n_minus_1 = SparsePolynomial {
            terms: BTreeMap::from([(
                Monomial(BTreeMap::from([(
                    t_var.to_string(),
                    (n - 1) as u32,
                )])),
                Expr::Constant(1.0),
            )]),
        };

        let coeff =
            simplify(&Expr::new_mul(
                Expr::new_mul(
                    q_n.clone(),
                    Expr::Constant(
                        n as f64,
                    ),
                ),
                dt_dx,
            ));

        let term2 =
            poly_mul_scalar_expr(
                &t_pow_n_minus_1,
                &coeff,
            );

        deriv = deriv + term2;
    }

    let mut p_star =
        (*p_in_t).clone() - deriv;

    p_star.prune_zeros(); // Remove zero coefficients to ensure degree decreases
    let recursive_integral =
        integrate_poly_log(
            &p_star, t, x,
        )?;

    let q_term_expr = Expr::new_mul(
        q_n,
        Expr::new_pow(
            t.clone(),
            Expr::Constant(n as f64),
        ),
    );

    Ok(simplify(
        &Expr::new_add(
            q_term_expr,
            recursive_integral,
        ),
    ))
}

pub(crate) fn find_outermost_transcendental(
    expr: &Expr,
    x: &str,
) -> Option<Expr> {

    let mut found_exp = None;

    let mut found_log = None;

    expr.pre_order_walk(&mut |e| {

        if let Expr::Exp(arg) = e {

            if contains_var(arg, x) {

                found_exp =
                    Some(e.clone());
            }
        }

        if let Expr::Log(arg) = e {

            if contains_var(arg, x) {

                found_log =
                    Some(e.clone());
            }
        }
    });

    found_exp.or(found_log)
}

/// Integrates the polynomial part of a transcendental function extension F(t).
///
/// This implementation handles the exponential case, where t = exp(g(x)).
/// Integrates the polynomial part of a transcendental function extension F(t).
/// This implementation handles the exponential case, where t = exp(g(x)).

pub fn integrate_poly_exp(
    p_in_t: &SparsePolynomial,
    t: &Expr,
    x: &str,
) -> Result<Expr, String> {

    let g = if let Expr::Exp(inner) = t
    {

        &**inner
    } else {

        return Err(
            "t is not exponential"
                .to_string(),
        );
    };

    let g_prime = differentiate(g, x);

    let p_coeffs =
        p_in_t.get_coeffs_as_vec(x);

    let n = p_in_t.degree(x) as usize;

    let mut q_coeffs =
        vec![
            Expr::Constant(0.0);
            n + 1
        ];

    for i in (0..=n).rev() {

        let p_i = p_coeffs
            .get(i)
            .cloned()
            .unwrap_or_else(|| {
                Expr::Constant(0.0)
            });

        let rhs = if i < n {

            let q_i_plus_1 =
                q_coeffs[i + 1].clone();

            let factor = Expr::new_mul(
                Expr::Constant(
                    (i + 1) as f64,
                ),
                g_prime.clone(),
            );

            simplify(&Expr::new_sub(
                p_i,
                Expr::new_mul(
                    factor, q_i_plus_1,
                ),
            ))
        } else {

            p_i
        };

        let q_i_var = format!("q_{i}");

        let q_i_expr = Expr::Variable(
            q_i_var.clone(),
        );

        let q_i_prime =
            differentiate(&q_i_expr, x);

        let ode_p_term =
            simplify(&Expr::new_mul(
                Expr::Constant(
                    i as f64,
                ),
                g_prime.clone(),
            ));

        let ode = simplify(&Expr::Eq(
            Arc::new(Expr::new_add(
                q_i_prime,
                Expr::new_mul(
                    ode_p_term,
                    q_i_expr,
                ),
            )),
            Arc::new(rhs),
        ));

        let sol_eq = crate::symbolic::ode::solve_ode(
            &ode, &q_i_var, x, None,
        );

        if let Expr::Eq(_, sol) = sol_eq
        {

            q_coeffs[i] =
                sol.as_ref().clone();
        } else {

            return Err(format!(
                "Failed to solve ODE \
                 for coefficient q_{i}"
            ));
        }
    }

    let q_poly =
        poly_from_coeffs(&q_coeffs, x);

    Ok(substitute(
        &sparse_poly_to_expr(&q_poly),
        x,
        t,
    ))
}

/// Helper to create a `SparsePolynomial` from a dense vector of coefficients.
#[must_use]

pub fn poly_from_coeffs(
    coeffs: &[Expr],
    var: &str,
) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    let n = coeffs.len() - 1;

    for (i, coeff) in coeffs
        .iter()
        .enumerate()
    {

        if !is_zero(&simplify(
            &coeff.clone(),
        )) {

            let mut mono_map =
                BTreeMap::new();

            let power = (n - i) as u32;

            if power > 0 {

                mono_map.insert(
                    var.to_string(),
                    power,
                );
            }

            terms.insert(
                Monomial(mono_map),
                coeff.clone(),
            );
        }
    }

    SparsePolynomial { terms }
}

/// Integrates a proper rational function A/B where B is square-free, using the Rothstein-Trager method.

pub fn partial_fraction_integrate(
    a: &SparsePolynomial,
    b: &SparsePolynomial,
    x: &str,
) -> Result<Expr, String> {

    let z =
        Expr::Variable("z".to_string());

    let b_prime =
        differentiate_poly(b, x);

    let r_poly_sym = a.clone()
        - (b_prime.clone()
            * poly_from_coeffs(
                &[z],
                x,
            ));

    let sylvester_mat =
        sylvester_matrix(
            &r_poly_sym,
            b,
            x,
        );

    let resultant =
        determinant(&sylvester_mat);

    let roots_c =
        solve(&resultant, "z");

    if roots_c.is_empty() {

        return Ok(Expr::Constant(0.0));
    }

    let mut total_log_sum =
        Expr::Constant(0.0);

    for c_i in roots_c {

        let a_minus_ci_b_prime = a.clone()
            - (b_prime.clone()
                * poly_from_coeffs(
                    std::slice::from_ref(&c_i),
                    x,
                ));

        let v_i = gcd(
            a_minus_ci_b_prime,
            b.clone(),
            x,
        );

        let log_term = Expr::new_log(
            sparse_poly_to_expr(&v_i),
        );

        let term =
            simplify(&Expr::new_mul(
                c_i, log_term,
            ));

        total_log_sum =
            simplify(&Expr::new_add(
                total_log_sum,
                term,
            ));
    }

    Ok(total_log_sum)
}

/// Constructs the Sylvester matrix of two polynomials.

pub(crate) fn sylvester_matrix(
    p: &SparsePolynomial,
    q: &SparsePolynomial,
    x: &str,
) -> Expr {

    let n = p.degree(x) as usize;

    let m = q.degree(x) as usize;

    let mut matrix =
        vec![
            vec![
                Expr::Constant(0.0);
                n + m
            ];
            n + m
        ];

    let p_coeffs =
        p.get_coeffs_as_vec(x);

    let q_coeffs =
        q.get_coeffs_as_vec(x);

    for i in 0..m {

        for j in 0..=n {

            matrix[i][i + j] = p_coeffs
                .get(j)
                .cloned()
                .unwrap_or_else(|| {
                    Expr::Constant(0.0)
                });
        }
    }

    for i in 0..n {

        for j in 0..=m {

            matrix[i + m][i + j] = q_coeffs
                .get(j)
                .cloned()
                .unwrap_or_else(|| Expr::Constant(0.0));
        }
    }

    Expr::Matrix(matrix)
}

/// Helper to integrate a simple polynomial.

pub(crate) fn poly_integrate(
    p: &SparsePolynomial,
    x: &str,
) -> Expr {

    let mut integral_expr =
        Expr::Constant(0.0);

    if p.terms.is_empty() {

        return integral_expr;
    }

    for (mono, coeff) in &p.terms {

        let exp = f64::from(
            mono.0
                .get(x)
                .copied()
                .unwrap_or(0),
        );

        let new_exp = exp + 1.0;

        let new_coeff =
            simplify(&Expr::new_div(
                coeff.clone(),
                Expr::Constant(new_exp),
            ));

        let term = Expr::new_mul(
            new_coeff,
            Expr::new_pow(
                Expr::Variable(
                    x.to_string(),
                ),
                Expr::Constant(new_exp),
            ),
        );

        integral_expr =
            simplify(&Expr::new_add(
                integral_expr,
                term,
            ));
    }

    integral_expr
}

pub fn hermite_integrate_rational(
    p: &SparsePolynomial,
    q: &SparsePolynomial,
    x: &str,
) -> Result<Expr, String> {

    /// Integrates a rational function `P(x)/Q(x)` using the Hermite-Ostrogradsky method.
    ///
    /// This function is a specialized version of `integrate_rational_function` that directly
    /// applies the Hermite-Ostrogradsky decomposition to a proper rational function.
    ///
    /// # Arguments
    /// * `p` - The numerator polynomial as a `SparsePolynomial`.
    /// * `q` - The denominator polynomial as a `SparsePolynomial`.
    /// * `x` - The variable of integration.
    ///
    /// # Returns
    /// A `Result` containing an `Expr` representing the integral, or an error string if computation fails.
    let (quotient, remainder) = p
        .clone()
        .long_division(q.clone(), x);

    let integral_of_quotient =
        poly_integrate(&quotient, x);

    if remainder
        .terms
        .is_empty()
    {

        return Ok(
            integral_of_quotient,
        );
    }

    let q_prime =
        differentiate_poly(q, x);

    let d = gcd(
        q.clone(),
        q_prime.clone(),
        x,
    );

    let b = q
        .clone()
        .long_division(d.clone(), x)
        .0;

    let (a_poly, c_poly) =
        build_and_solve_hermite_system(
            &remainder, &b, &d,
            &q_prime, x,
        )?;

    let rational_part = Expr::new_div(
        sparse_poly_to_expr(&c_poly),
        sparse_poly_to_expr(&d),
    );

    let integral_of_transcendental_part = integrate_square_free_rational_part(&a_poly, &b, x)?;

    Ok(simplify(
        &Expr::new_add(
            integral_of_quotient,
            Expr::new_add(
                rational_part,
                integral_of_transcendental_part,
            ),
        ),
    ))
}

/// Integrates a rational function A/B where B is square-free, using the Rothstein-Trager method.

pub(crate) fn integrate_square_free_rational_part(
    a: &SparsePolynomial,
    b: &SparsePolynomial,
    x: &str,
) -> Result<Expr, String> {

    let z =
        Expr::Variable("z".to_string());

    let b_prime =
        differentiate_poly(b, x);

    let r_poly_sym = a.clone()
        - (b_prime.clone()
            * expr_to_sparse_poly(&z));

    let sylvester_mat =
        sylvester_matrix(
            &r_poly_sym,
            b,
            x,
        );

    let resultant =
        determinant(&sylvester_mat);

    let roots_c =
        solve(&resultant, "z");

    if roots_c.is_empty() {

        return Ok(Expr::Constant(0.0));
    }

    let mut total_log_sum =
        Expr::Constant(0.0);

    for c_i in roots_c {

        let a_minus_ci_b_prime = a
            .clone()
            - (b_prime.clone()
                * expr_to_sparse_poly(
                    &c_i,
                ));

        let v_i = gcd(
            a_minus_ci_b_prime,
            b.clone(),
            x,
        );

        let log_term = Expr::new_log(
            sparse_poly_to_expr(&v_i),
        );

        let term =
            simplify(&Expr::new_mul(
                c_i, log_term,
            ));

        total_log_sum =
            simplify(&Expr::new_add(
                total_log_sum,
                term,
            ));
    }

    Ok(total_log_sum)
}

/// Converts an expression into a rational function A(t)/D(t) of a transcendental element t.

pub(crate) fn expr_to_rational_poly(
    expr: &Expr,
    t: &Expr,
    _x: &str,
) -> Result<
    (
        SparsePolynomial,
        SparsePolynomial,
    ),
    String,
> {

    // Substitute t with a variable "t_var" to convert to polynomial
    // First, we need to replace occurrences of t in expr with a variable
    let t_var_name = "t_var";

    let expr_with_t_var =
        substitute_expr_for_var(
            expr, t, t_var_name,
        );

    let poly = crate::symbolic::polynomial::expr_to_sparse_poly(
        &expr_with_t_var,
        &[t_var_name],
    );

    let one_poly = SparsePolynomial {
        terms: BTreeMap::from([(
            Monomial(BTreeMap::new()),
            Expr::Constant(1.0),
        )]),
    };

    Ok((poly, one_poly))
}

/// Helper to substitute an expression with a variable name

fn substitute_expr_for_var(
    expr: &Expr,
    target: &Expr,
    var_name: &str,
) -> Expr {

    if expr == target {

        return Expr::Variable(
            var_name.to_string(),
        );
    }

    match expr {
        | Expr::Add(a, b) => {
            Expr::new_add(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
                substitute_expr_for_var(
                    b, target, var_name,
                ),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
                substitute_expr_for_var(
                    b, target, var_name,
                ),
            )
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
                substitute_expr_for_var(
                    b, target, var_name,
                ),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
                substitute_expr_for_var(
                    b, target, var_name,
                ),
            )
        },
        | Expr::Power(a, b) => {
            Expr::new_pow(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
                substitute_expr_for_var(
                    b, target, var_name,
                ),
            )
        },
        | Expr::Neg(a) => {
            Expr::new_neg(
                substitute_expr_for_var(
                    a, target, var_name,
                ),
            )
        },
        | _ => expr.clone(),
    }
}

pub fn integrate_rational_function_expr(
    expr: &Expr,
    x: &str,
) -> Result<Expr, String> {

    let p = expr_to_sparse_poly(expr);

    let q = SparsePolynomial {
        terms: BTreeMap::from([(
            Monomial(BTreeMap::new()),
            Expr::Constant(1.0),
        )]),
    };

    integrate_rational_function(
        &p, &q, x,
    )
}

#[must_use]

pub fn poly_derivative_symbolic(
    p: &SparsePolynomial,
    x: &str,
) -> SparsePolynomial {

    differentiate_poly(p, x)
}
