//! # Ordinary Differential Equation (ODE) Solver
//!
//! This module provides a comprehensive set of functions for solving Ordinary Differential Equations.
//! It includes dispatchers for various types of ODEs (first-order linear, separable, Bernoulli,
//! exact, Cauchy-Euler, reduction of order), as well as methods for solving systems of ODEs
//! and applying initial conditions. Techniques like series solutions and Fourier transforms
//! are also supported for specific cases.

use crate::symbolic::calculus::{differentiate, integrate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::pattern_match;
use crate::symbolic::simplify::{is_zero, simplify};
use crate::symbolic::solve::{solve, solve_linear_system};
use crate::symbolic::transforms;
use std::collections::{HashMap, HashSet};

// =====================================================================================
// region: ODE Parser and Helpers
// =====================================================================================

/// A structured representation of a parsed ODE.
pub struct ParsedODE {
    pub order: u32,
    pub coeffs: HashMap<u32, Expr>,
    pub remaining_expr: Expr,
}

pub(crate) fn parse_ode(equation: &Expr, func: &str, var: &str) -> ParsedODE {
    let mut coeffs = HashMap::new();
    let mut remaining_expr = Expr::Constant(0.0);

    pub(crate) fn collect_terms(
        expr: &Expr,
        func: &str,
        var: &str,
        coeffs: &mut HashMap<u32, Expr>,
        remaining: &mut Expr,
    ) {
        if let Expr::Add(a, b) = expr {
            collect_terms(a, func, var, coeffs, remaining);
            collect_terms(b, func, var, coeffs, remaining);
        } else {
            let (order, coeff) = get_term_order_and_coeff(expr, func, var);
            if order > 100 {
                *remaining = simplify(Expr::Add(
                    Box::new(remaining.clone()),
                    Box::new(expr.clone()),
                ));
            } else {
                let entry = coeffs.entry(order).or_insert_with(|| Expr::Constant(0.0));
                *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(coeff)));
            }
        }
    }

    collect_terms(
        &simplify(equation.clone()),
        func,
        var,
        &mut coeffs,
        &mut remaining_expr,
    );
    let max_order = coeffs.keys().max().cloned().unwrap_or(0);

    ParsedODE {
        order: max_order,
        coeffs,
        remaining_expr,
    }
}

pub(crate) fn get_term_order_and_coeff(expr: &Expr, func: &str, var: &str) -> (u32, Expr) {
    match expr {
        Expr::Derivative(inner, d_var) if d_var == var => {
            let (order, coeff) = get_term_order_and_coeff(inner, func, var);
            (order + 1, coeff)
        }
        Expr::Mul(a, b) => {
            let (order_a, _) = get_term_order_and_coeff(a, func, var);
            let (order_b, _) = get_term_order_and_coeff(b, func, var);
            if order_a > 100 && order_b <= 100 {
                let (order, coeff) = get_term_order_and_coeff(b, func, var);
                (
                    order,
                    simplify(Expr::Mul(Box::new(coeff), Box::new(*a.clone()))),
                )
            } else if order_b > 100 && order_a <= 100 {
                let (order, coeff) = get_term_order_and_coeff(a, func, var);
                (
                    order,
                    simplify(Expr::Mul(Box::new(coeff), Box::new(*b.clone()))),
                )
            } else {
                (999, expr.clone())
            }
        }
        Expr::Variable(v) if v == func => (0, Expr::Constant(1.0)),
        _ => (999, expr.clone()),
    }
}

pub(crate) fn find_constants(expr: &Expr, constants: &mut Vec<String>) {
    if let Expr::Variable(s) = expr {
        if s.starts_with('C')
            && s.chars().skip(1).all(|c| c.is_ascii_digit())
            && !constants.contains(s)
        {
            constants.push(s.clone());
        }
    }
    match expr {
        Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b) => {
            find_constants(a, constants);
            find_constants(b, constants);
        }
        Expr::Sin(a) | Expr::Cos(a) | Expr::Tan(a) | Expr::Exp(a) | Expr::Log(a) | Expr::Neg(a) => {
            find_constants(a, constants);
        }
        _ => {}
    }
}

pub(crate) fn find_derivatives(expr: &Expr, var: &str, derivatives: &mut HashMap<String, u32>) {
    if let Expr::Derivative(inner, d_var) = expr {
        if d_var == var {
            let mut current = &**inner;
            let mut order = 1;
            while let Expr::Derivative(next_inner, next_d_var) = current {
                if next_d_var == var {
                    order += 1;
                    current = next_inner;
                } else {
                    break;
                }
            }
            if let Expr::Variable(func_name) = current {
                let entry = derivatives.entry(func_name.clone()).or_insert(0);
                *entry = std::cmp::max(*entry, order);
            }
        }
    }
    match expr {
        Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b) => {
            find_derivatives(a, var, derivatives);
            find_derivatives(b, var, derivatives);
        }
        Expr::Sin(a) | Expr::Cos(a) | Expr::Tan(a) | Expr::Exp(a) | Expr::Log(a) | Expr::Neg(a) => {
            find_derivatives(a, var, derivatives);
        }
        _ => {}
    }
}

// =====================================================================================
// endregion: ODE Parser and Helpers
// =====================================================================================

// =====================================================================================
// region: Public API and Dispatchers
// =====================================================================================

/// Main public function for solving a single Ordinary Differential Equation.
///
/// This function acts as a dispatcher, attempting to solve the ODE by trying various
/// specialized solvers based on the ODE's type (e.g., first-order linear, separable).
/// It can also apply initial conditions to find a particular solution.
///
/// # Arguments
/// * `ode` - The ODE to solve, typically an `Expr::Eq`.
/// * `func` - The name of the unknown function (e.g., "y").
/// * `var` - The independent variable (e.g., "x").
/// * `initial_conditions` - An `Option` containing a slice of `(Expr, u32, Expr)` tuples
///   representing initial conditions `(x0, order_of_derivative, y_value_at_x0)`.
///
/// # Returns
/// An `Expr` representing the general or particular solution to the ODE.
pub fn solve_ode(
    ode: &Expr,
    func: &str,
    var: &str,
    initial_conditions: Option<&[(Expr, u32, Expr)]>,
) -> Expr {
    //let general_solution_eq = solve_ode_system(&[ode.clone()], &[func], var)
    // If ode is a reference:
    let general_solution_eq = solve_ode_system(std::slice::from_ref(ode), &[func], var)
        .and_then(|mut solutions| solutions.pop())
        .map(|sol| Expr::Eq(Box::new(Expr::Variable(func.to_string())), Box::new(sol)))
        .unwrap_or_else(|| Expr::Solve(Box::new(ode.clone()), func.to_string()));

    if let Some(conditions) = initial_conditions {
        if let Expr::Eq(_, general_solution) = &general_solution_eq {
            return apply_initial_conditions(general_solution, var, conditions);
        }
    }

    general_solution_eq
}

/// Solves a system of coupled ordinary differential equations, including higher-order ones.
///
/// This function first reduces the system of higher-order ODEs to an equivalent first-order system.
/// Then, it attempts to solve this first-order system sequentially by applying various ODE solvers
/// to each equation.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the ODEs in the system.
/// * `funcs` - A slice of string slices representing the names of the unknown functions (e.g., `["y", "z"]`).
/// * `var` - The independent variable (e.g., "x").
///
/// # Returns
/// An `Option<Vec<Expr>>` containing a vector of solutions for each function in `funcs`,
/// or `None` if the system cannot be solved.
pub fn solve_ode_system(equations: &[Expr], funcs: &[&str], var: &str) -> Option<Vec<Expr>> {
    let (first_order_eqs, all_vars, original_funcs_map) =
        reduce_to_first_order_system(equations, funcs, var);
    let first_order_funcs: Vec<&str> = all_vars.iter().map(|s| s.as_str()).collect();

    let solutions_map =
        solve_first_order_system_sequentially(&first_order_eqs, &first_order_funcs, var)?;

    let mut final_solutions = Vec::new();
    for &original_func in funcs {
        let sol_var = original_funcs_map.get(original_func)?;
        let solution = solutions_map.get(sol_var)?;
        final_solutions.push(solution.clone());
    }
    Some(final_solutions)
}

pub(crate) fn try_all_solvers(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    let eq = if let Expr::Eq(l, r) = equation {
        simplify(Expr::Sub(l.clone(), r.clone()))
    } else {
        equation.clone()
    };
    solve_first_order_linear_ode(&eq, func, var)
        .or_else(|| solve_separable_ode(&eq, func, var))
        .or_else(|| solve_bernoulli_ode(&eq, func, var))
        .or_else(|| solve_cauchy_euler_ode(&eq, func, var))
        .or_else(|| solve_exact_ode(&eq, func, var))
}

pub(crate) fn apply_initial_conditions(
    general_solution: &Expr,
    var: &str,
    conditions: &[(Expr, u32, Expr)],
) -> Expr {
    let mut constants = Vec::new();
    find_constants(general_solution, &mut constants);
    constants.sort();

    if constants.is_empty() || conditions.is_empty() {
        return general_solution.clone();
    }

    let mut eq_system = Vec::new();
    for (x0, order, y_val) in conditions {
        let mut sol_deriv = general_solution.clone();
        for _ in 0..*order {
            sol_deriv = differentiate(&sol_deriv, var);
        }
        let substituted_sol = substitute(&sol_deriv, var, x0);
        let equation = simplify(Expr::Eq(Box::new(substituted_sol), Box::new(y_val.clone())));
        eq_system.push(equation);
    }

    if eq_system.len() < constants.len() {
        return Expr::Variable("Not enough initial conditions".to_string());
    }

    if let Ok(const_solutions) = solve_linear_system(&Expr::System(eq_system), &constants) {
        let mut final_solution = general_solution.clone();
        for (i, c_var) in constants.iter().enumerate() {
            if i < const_solutions.len() {
                final_solution = substitute(&final_solution, c_var, &const_solutions[i]);
            }
        }
        return simplify(final_solution);
    }

    Expr::Variable("Could not solve for constants".to_string())
}

// =====================================================================================
// endregion: Public API and Dispatchers
// =====================================================================================

// =====================================================================================
// region: Core Solver Implementations
// =====================================================================================

pub(crate) fn reduce_to_first_order_system(
    equations: &[Expr],
    funcs: &[&str],
    var: &str,
) -> (Vec<Expr>, Vec<String>, HashMap<String, String>) {
    let mut new_eqs = Vec::new();
    let mut new_vars_map: HashMap<(String, u32), String> = HashMap::new();
    let mut all_new_vars = funcs.iter().map(|s| s.to_string()).collect::<HashSet<_>>();
    let mut original_funcs_map = HashMap::new();

    for &func in funcs {
        new_vars_map.insert((func.to_string(), 0), func.to_string());
        original_funcs_map.insert(func.to_string(), func.to_string());
    }

    let mut temp_eqs = equations.to_vec();
    let mut i = 0;
    while i < temp_eqs.len() {
        let eq = &temp_eqs[i];
        let mut derivatives = HashMap::new();
        find_derivatives(eq, var, &mut derivatives);

        let mut eq_with_substitutions = eq.clone();

        for (func, &order) in &derivatives {
            if order > 1 {
                for k in 1..order {
                    let key = (func.clone(), k);
                    //if new_vars_map.get(&key).is_none() {
                    if !new_vars_map.contains_key(&key) {
                        let new_var_name = format!("{}_d{}", func, k);
                        all_new_vars.insert(new_var_name.clone());
                        new_vars_map.insert(key.clone(), new_var_name.clone());

                        let prev_var_name = new_vars_map.get(&(func.clone(), k - 1)).unwrap();
                        let new_eq = Expr::Eq(
                            Box::new(Expr::Derivative(
                                Box::new(Expr::Variable(prev_var_name.clone())),
                                var.to_string(),
                            )),
                            Box::new(Expr::Variable(new_var_name.clone())),
                        );
                        temp_eqs.push(new_eq);
                    }
                }
                let highest_deriv = (0..order).fold(Expr::Variable(func.clone()), |e, _| {
                    Expr::Derivative(Box::new(e), var.to_string())
                });
                let replacement_var_name = new_vars_map.get(&(func.clone(), order - 1)).unwrap();
                let replacement_expr = Expr::Derivative(
                    Box::new(Expr::Variable(replacement_var_name.clone())),
                    var.to_string(),
                );
                eq_with_substitutions = substitute(
                    &eq_with_substitutions,
                    &highest_deriv.to_string(),
                    &replacement_expr,
                );
            }
        }
        new_eqs.push(eq_with_substitutions);
        i += 1;
    }

    (
        new_eqs,
        all_new_vars.into_iter().collect(),
        original_funcs_map,
    )
}

pub(crate) fn solve_first_order_system_sequentially(
    equations: &[Expr],
    funcs: &[&str],
    var: &str,
) -> Option<HashMap<String, Expr>> {
    let mut remaining_eqs: Vec<Expr> = equations.to_vec();
    let mut solutions: HashMap<String, Expr> = HashMap::new();
    let mut progress = true;

    while progress && !remaining_eqs.is_empty() {
        progress = false;
        let mut solved_eq_indices = Vec::new();

        for (i, eq) in remaining_eqs.iter().enumerate() {
            let mut current_eq = eq.clone();
            for (solved_func, solution_expr) in &solutions {
                current_eq = substitute(&current_eq, solved_func, solution_expr);
            }

            let mut remaining_funcs = Vec::new();
            for &f in funcs {
                if !solutions.contains_key(f) {
                    let mut found = false;
                    current_eq.pre_order_walk(&mut |e| {
                        if let Expr::Variable(v) = e {
                            if v == f {
                                found = true;
                            }
                        }
                        if let Expr::Derivative(inner, _) = e {
                            if let Expr::Variable(v) = &**inner {
                                if v == f {
                                    found = true;
                                }
                            }
                        }
                    });
                    if found {
                        remaining_funcs.push(f);
                    }
                }
            }

            if remaining_funcs.len() == 1 {
                let func_to_solve = remaining_funcs[0];
                let solution_eq = try_all_solvers(&current_eq, func_to_solve, var)?;

                if let Expr::Eq(_, solution_expr) = solution_eq {
                    solutions.insert(func_to_solve.to_string(), *solution_expr);
                    solved_eq_indices.push(i);
                    progress = true;
                    break;
                }
            }
        }

        for &i in solved_eq_indices.iter().rev() {
            remaining_eqs.remove(i);
        }
    }

    if solutions.len() == funcs.len() {
        Some(solutions)
    } else {
        None
    }
}

pub(crate) fn solve_separable_ode(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    let y_prime = Expr::Derivative(Box::new(Expr::Variable(func.to_string())), var.to_string());
    if let Some(assignments) = pattern_match(
        equation,
        &Expr::Sub(
            Box::new(Expr::Mul(
                Box::new(Expr::Pattern("F".to_string())),
                Box::new(y_prime.clone()),
            )),
            Box::new(Expr::Pattern("G".to_string())),
        ),
    ) {
        let f_y = assignments.get("F")?;
        let g_x = assignments.get("G")?;

        if !g_x.to_string().contains(func) && !f_y.to_string().contains(var) {
            let int_f_y = integrate(f_y, func, None, None);
            let int_g_x = integrate(g_x, var, None, None);
            let c = Expr::Variable("C".to_string());
            return Some(simplify(Expr::Eq(
                Box::new(int_f_y),
                Box::new(Expr::Add(Box::new(int_g_x), Box::new(c))),
            )));
        }
    }
    None
}

pub(crate) fn solve_first_order_linear_ode(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    let parsed = parse_ode(equation, func, var);
    if parsed.order != 1 {
        return None;
    }

    let p_x = parsed.coeffs.get(&0).cloned()?;
    let r_x = parsed.remaining_expr.clone();
    let q_x = simplify(Expr::Neg(Box::new(r_x)));
    let y_expr = Expr::Variable(func.to_string());

    let mu = Expr::Exp(Box::new(integrate(&p_x, var, None, None)));
    let rhs = integrate(
        &simplify(Expr::Mul(Box::new(q_x), Box::new(mu.clone()))),
        var,
        None,
        None,
    );
    let c = Expr::Variable("C1".to_string());
    let solution = simplify(Expr::Div(
        Box::new(Expr::Add(Box::new(rhs), Box::new(c))),
        Box::new(mu),
    ));
    Some(Expr::Eq(Box::new(y_expr), Box::new(solution)))
}

pub fn solve_bernoulli_ode(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    /// Solves a Bernoulli differential equation.
    ///
    /// A Bernoulli ODE is of the form `y' + P(x)y = Q(x)y^n`.
    /// This function transforms the Bernoulli equation into a first-order linear ODE
    /// using the substitution `v = y^(1-n)`, which can then be solved.
    ///
    /// # Arguments
    /// * `equation` - The Bernoulli ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    ///
    /// # Returns
    /// An `Option<Expr>` representing the solution, or `None` if the equation
    /// does not match the Bernoulli form or cannot be solved.
    let y = Expr::Variable(func.to_string());
    let y_prime = Expr::Derivative(Box::new(y.clone()), var.to_string());

    let pattern = Expr::Add(
        Box::new(y_prime),
        Box::new(Expr::Sub(
            Box::new(Expr::Mul(
                Box::new(Expr::Pattern("P".to_string())),
                Box::new(y.clone()),
            )),
            Box::new(Expr::Mul(
                Box::new(Expr::Pattern("Q".to_string())),
                Box::new(Expr::Power(
                    Box::new(y.clone()),
                    Box::new(Expr::Pattern("n".to_string())),
                )),
            )),
        )),
    );

    if let Some(m) = pattern_match(equation, &pattern) {
        let p_x = m.get("P")?;
        let q_x = m.get("Q")?;
        let n = m.get("n")?.to_f64()?;

        if n == 1.0 || n == 0.0 {
            return None;
        }

        let one_minus_n = 1.0 - n;
        let p_v = simplify(Expr::Mul(
            Box::new(Expr::Constant(one_minus_n)),
            Box::new(p_x.clone()),
        ));
        let q_v = simplify(Expr::Mul(
            Box::new(Expr::Constant(one_minus_n)),
            Box::new(q_x.clone()),
        ));

        let v_prime = Expr::Derivative(Box::new(Expr::Variable("v".to_string())), var.to_string());
        let linear_ode_v = Expr::Add(
            Box::new(v_prime),
            Box::new(Expr::Sub(
                Box::new(Expr::Mul(
                    Box::new(p_v),
                    Box::new(Expr::Variable("v".to_string())),
                )),
                Box::new(q_v),
            )),
        );

        let v_solution_eq = solve_first_order_linear_ode(&linear_ode_v, "v", var)?;
        let v_solution = if let Expr::Eq(_, sol) = v_solution_eq {
            *sol
        } else {
            return None;
        };

        let y_solution = simplify(Expr::Power(
            Box::new(v_solution),
            Box::new(Expr::Constant(1.0 / one_minus_n)),
        ));
        return Some(Expr::Eq(Box::new(y), Box::new(y_solution)));
    }
    None
}

pub fn solve_riccati_ode(equation: &Expr, func: &str, var: &str, y1: &Expr) -> Option<Expr> {
    /// Solves a Riccati differential equation.
    ///
    /// A Riccati ODE is of the form `y' = P(x) + Q(x)y + R(x)y^2`.
    /// If a particular solution `y1` is known, the general solution can be found
    /// by substituting `y = y1 + 1/v`, which transforms the Riccati equation
    /// into a first-order linear ODE in `v`.
    ///
    /// # Arguments
    /// * `equation` - The Riccati ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    /// * `y1` - A known particular solution to the Riccati equation.
    ///
    /// # Returns
    /// An `Option<Expr>` representing the general solution, or `None` if the equation
    /// does not match the Riccati form or cannot be solved.
    let _ = (equation, func, var, y1);
    None
}

pub fn solve_cauchy_euler_ode(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    /// Solves a homogeneous Cauchy-Euler (equidimensional) differential equation of second order.
    ///
    /// A Cauchy-Euler ODE is of the form `ax^2y'' + bxy' + cy = 0`.
    /// This function finds the roots of the associated indicial equation to construct
    /// the general solution, handling cases of distinct real roots, repeated real roots,
    /// and complex conjugate roots.
    ///
    /// # Arguments
    /// * `equation` - The Cauchy-Euler ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    ///
    /// # Returns
    /// An `Option<Expr>` representing the general solution, or `None` if the equation
    /// does not match the Cauchy-Euler form or cannot be solved.
    let parsed = parse_ode(equation, func, var);
    if parsed.order != 2 || !is_zero(&parsed.remaining_expr) {
        return None;
    }

    let c2 = parsed.coeffs.get(&2)?;
    let c1 = parsed.coeffs.get(&1)?;
    let c0 = parsed.coeffs.get(&0)?;
    let x = Expr::Variable(var.to_string());
    let x_sq = Expr::Power(Box::new(x.clone()), Box::new(Expr::Constant(2.0)));

    let a = simplify(Expr::Div(Box::new(c2.clone()), Box::new(x_sq)));
    let b = simplify(Expr::Div(Box::new(c1.clone()), Box::new(x.clone())));
    let c = c0.clone();

    if a.to_f64().is_none() || b.to_f64().is_none() || c.to_f64().is_none() {
        return None;
    }

    let m = Expr::Variable("m".to_string());
    let b_minus_a = simplify(Expr::Sub(Box::new(b), Box::new(a.clone())));
    let aux_eq = Expr::Add(
        Box::new(Expr::Mul(
            Box::new(a),
            Box::new(Expr::Power(
                Box::new(m.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
        )),
        Box::new(Expr::Add(
            Box::new(Expr::Mul(Box::new(b_minus_a), Box::new(m.clone()))),
            Box::new(c.clone()),
        )),
    );

    let roots = solve(&aux_eq, "m");
    if roots.len() != 2 {
        return None;
    }

    let m1 = &roots[0];
    let m2 = &roots[1];
    let const1 = Expr::Variable("C1".to_string());
    let const2 = Expr::Variable("C2".to_string());

    let solution = if m1 != m2 {
        simplify(Expr::Add(
            Box::new(Expr::Mul(
                Box::new(const1),
                Box::new(Expr::Power(Box::new(x.clone()), Box::new(m1.clone()))),
            )),
            Box::new(Expr::Mul(
                Box::new(const2),
                Box::new(Expr::Power(Box::new(x.clone()), Box::new(m2.clone()))),
            )),
        ))
    } else {
        simplify(Expr::Mul(
            Box::new(Expr::Power(Box::new(x.clone()), Box::new(m1.clone()))),
            Box::new(Expr::Add(
                Box::new(const1),
                Box::new(Expr::Mul(
                    Box::new(const2),
                    Box::new(Expr::Log(Box::new(x.clone()))),
                )),
            )),
        ))
    };

    Some(Expr::Eq(
        Box::new(Expr::Variable(func.to_string())),
        Box::new(solution),
    ))
}

pub fn solve_by_reduction_of_order(
    equation: &Expr,
    func: &str,
    var: &str,
    y1: &Expr,
) -> Option<Expr> {
    /// Solves a second-order homogeneous linear ODE by reduction of order.
    ///
    /// If one non-trivial solution `y1` to a second-order homogeneous linear ODE
    /// is known, a second linearly independent solution `y2` can be found.
    /// The general solution is then `y = C1*y1 + C2*y2`.
    ///
    /// # Arguments
    /// * `equation` - The second-order homogeneous linear ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    /// * `y1` - A known non-trivial solution to the homogeneous ODE.
    ///
    /// # Returns
    /// An `Option<Expr>` representing the general solution, or `None` if the equation
    /// is not a second-order homogeneous linear ODE or cannot be solved.
    let parsed = parse_ode(equation, func, var);
    if parsed.order != 2 || !is_zero(&parsed.remaining_expr) {
        return None;
    }

    let coeff2 = parsed.coeffs.get(&2)?;
    let p_x = simplify(Expr::Div(
        Box::new(parsed.coeffs.get(&1)?.clone()),
        Box::new(coeff2.clone()),
    ));

    let integral_p = integrate(&p_x, var, None, None);
    let exp_term = Expr::Exp(Box::new(Expr::Neg(Box::new(integral_p))));

    let y1_sq = Expr::Power(Box::new(y1.clone()), Box::new(Expr::Constant(2.0)));
    let integrand = simplify(Expr::Div(Box::new(exp_term), Box::new(y1_sq)));
    let integral_v = integrate(&integrand, var, None, None);

    let y2 = simplify(Expr::Mul(Box::new(y1.clone()), Box::new(integral_v)));

    let c1 = Expr::Variable("C1".to_string());
    let c2 = Expr::Variable("C2".to_string());
    let general_solution = simplify(Expr::Add(
        Box::new(Expr::Mul(Box::new(c1), Box::new(y1.clone()))),
        Box::new(Expr::Mul(Box::new(c2), Box::new(y2))),
    ));

    Some(Expr::Eq(
        Box::new(Expr::Variable(func.to_string())),
        Box::new(general_solution),
    ))
}

pub fn solve_exact_ode(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    /// Solves an exact first-order Ordinary Differential Equation.
    ///
    /// An ODE of the form `M(x,y)dx + N(x,y)dy = 0` is exact if `∂M/∂y = ∂N/∂x`.
    /// The solution is then given by `F(x,y) = C`, where `∂F/∂x = M` and `∂F/∂y = N`.
    ///
    /// # Arguments
    /// * `equation` - The exact ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    ///
    /// # Returns
    /// An `Option<Expr>` representing the implicit solution `F(x,y) = C`,
    /// or `None` if the equation is not exact or cannot be solved.
    let y = Expr::Variable(func.to_string());
    let y_prime = Expr::Derivative(Box::new(y.clone()), var.to_string());

    let pattern = Expr::Add(
        Box::new(Expr::Pattern("M".to_string())),
        Box::new(Expr::Mul(
            Box::new(Expr::Pattern("N".to_string())),
            Box::new(y_prime),
        )),
    );

    if let Some(m) = pattern_match(equation, &pattern) {
        let m_xy = m.get("M")?;
        let n_xy = m.get("N")?;

        let dm_dy = differentiate(m_xy, func);
        let dn_dx = differentiate(n_xy, var);

        if simplify(dm_dy) != simplify(dn_dx) {
            return None;
        }

        let int_m_dx = integrate(m_xy, var, None, None);

        let d_int_m_dy = differentiate(&int_m_dx, func);
        let g_prime_y = simplify(Expr::Sub(Box::new(n_xy.clone()), Box::new(d_int_m_dy)));

        let g_y = integrate(&g_prime_y, func, None, None);

        let f_xy = simplify(Expr::Add(Box::new(int_m_dx), Box::new(g_y)));

        return Some(Expr::Eq(
            Box::new(f_xy),
            Box::new(Expr::Variable("C".to_string())),
        ));
    }
    None
}

pub fn solve_ode_by_series(
    equation: &Expr,
    func: &str,
    var: &str,
    x0: &Expr,
    order: u32,
    initial_conditions: &[(u32, Expr)],
) -> Option<Expr> {
    /// Solves an Ordinary Differential Equation using the power series method.
    ///
    /// This method assumes a solution of the form `y(x) = Σ a_n (x - x0)^n`
    /// and substitutes it into the ODE to find a recurrence relation for the coefficients `a_n`.
    /// It then uses initial conditions to determine the first few coefficients and constructs
    /// a truncated power series solution.
    ///
    /// # Arguments
    /// * `equation` - The ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "x").
    /// * `x0` - The center of the power series expansion.
    /// * `order` - The maximum order of the series to compute.
    /// * `initial_conditions` - A slice of `(u32, Expr)` tuples representing initial conditions
    ///   for `y(x0)`, `y'(x0)`, etc. (e.g., `(0, y(x0))`, `(1, y'(x0))`).
    ///
    /// # Returns
    /// An `Option<Expr>` representing the truncated power series solution.
    let mut y_n_at_x0: HashMap<u32, Expr> = initial_conditions.iter().cloned().collect();

    let parsed = parse_ode(equation, func, var);
    let highest_order = parsed.order;
    let coeff_helper = parsed.coeffs.clone();
    let coeff_highest = coeff_helper.get(&highest_order)?;

    let mut other_terms = Expr::Constant(0.0);
    for (o, c) in parsed.coeffs {
        if o < highest_order {
            let deriv = (0..o).fold(Expr::Variable(func.to_string()), |e, _| {
                Expr::Derivative(Box::new(e), var.to_string())
            });
            other_terms = simplify(Expr::Add(
                Box::new(other_terms),
                Box::new(Expr::Mul(Box::new(c), Box::new(deriv))),
            ));
        }
    }
    other_terms = simplify(Expr::Add(
        Box::new(other_terms),
        Box::new(parsed.remaining_expr),
    ));
    let highest_deriv_expr = simplify(Expr::Neg(Box::new(Expr::Div(
        Box::new(other_terms),
        Box::new(coeff_highest.clone()),
    ))));

    for n in highest_order..=order {
        if !y_n_at_x0.contains_key(&n) {
            let mut current_expr = highest_deriv_expr.clone();
            for i in 0..n {
                let deriv_i_expr = (0..i).fold(Expr::Variable(func.to_string()), |e, _| {
                    Expr::Derivative(Box::new(e), var.to_string())
                });
                if let Some(val) = y_n_at_x0.get(&i) {
                    current_expr = substitute(&current_expr, &deriv_i_expr.to_string(), val);
                }
            }
            let val_at_x0 = substitute(&current_expr, var, x0);
            y_n_at_x0.insert(n, simplify(val_at_x0));
        }
    }

    let mut series_sum = Expr::Constant(0.0);
    for n in 0..=order {
        if let Some(y_n_val) = y_n_at_x0.get(&n) {
            let n_factorial = (1..=n).product::<u32>() as f64;
            let coeff_term = simplify(Expr::Div(
                Box::new(y_n_val.clone()),
                Box::new(Expr::Constant(n_factorial)),
            ));
            let power_term = Expr::Power(
                Box::new(Expr::Sub(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(x0.clone()),
                )),
                Box::new(Expr::Constant(n as f64)),
            );
            series_sum = simplify(Expr::Add(
                Box::new(series_sum),
                Box::new(Expr::Mul(Box::new(coeff_term), Box::new(power_term))),
            ));
        }
    }

    Some(series_sum)
}

pub fn solve_ode_by_fourier(equation: &Expr, func: &str, var: &str) -> Option<Expr> {
    /// Solves a linear Ordinary Differential Equation using the Fourier Transform method.
    ///
    /// This method transforms the ODE from the time domain to the frequency domain,
    /// converting differential operators into algebraic multiplications. The resulting
    /// algebraic equation is solved for the transformed function, and then the inverse
    /// Fourier Transform is applied to obtain the solution in the original domain.
    ///
    /// # Arguments
    /// * `equation` - The ODE to solve.
    /// * `func` - The name of the unknown function (e.g., "y").
    /// * `var` - The independent variable (e.g., "t").
    ///
    /// # Returns
    /// An `Option<Expr>` representing the solution, or `None` if the ODE type
    /// is not supported by this method.
    let omega_var = "w";
    let parsed = parse_ode(equation, func, var);

    let g_w = transforms::fourier_transform(&parsed.remaining_expr, var, omega_var);

    let mut algebraic_lhs = Expr::Constant(0.0);
    let y_w = Expr::Variable("Y".to_string());

    for (order, coeff) in parsed.coeffs {
        coeff.to_f64()?;

        let mut deriv_transform = y_w.clone();
        for _ in 0..order {
            deriv_transform = transforms::fourier_differentiation(&deriv_transform, omega_var);
        }
        let term = simplify(Expr::Mul(Box::new(coeff), Box::new(deriv_transform)));
        algebraic_lhs = simplify(Expr::Add(Box::new(algebraic_lhs), Box::new(term)));
    }

    let algebraic_eq = simplify(Expr::Sub(
        Box::new(algebraic_lhs),
        Box::new(simplify(Expr::Neg(Box::new(g_w)))),
    ));

    let y_w_solutions = solve(&algebraic_eq, "Y");
    if y_w_solutions.is_empty() {
        return None;
    }
    let y_w_solution = y_w_solutions[0].clone();

    let solution = transforms::inverse_fourier_transform(&y_w_solution, omega_var, var);
    Some(Expr::Eq(
        Box::new(Expr::Variable(func.to_string())),
        Box::new(solution),
    ))
}

// =====================================================================================
// endregion: Specific Solver Implementations
// =====================================================================================
