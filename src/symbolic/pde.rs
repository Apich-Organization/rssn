//! # Partial Differential Equation (PDE) Solver
//!
//! This module provides functions for solving various types of Partial Differential Equations.
//! It includes strategies for first-order PDEs (method of characteristics), second-order PDEs
//! (separation of variables, D'Alembert's formula for wave equation), and techniques like
//! Green's functions and Fourier transforms for specific PDE types.

use crate::symbolic::calculus::{differentiate, integrate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::ode::solve_ode;
use crate::symbolic::simplify::{collect_and_order_terms, is_zero, pattern_match, simplify};
use crate::symbolic::transforms;
use std::collections::HashMap;

// =====================================================================================
// region: Public API and Dispatchers
// =====================================================================================

/// Main dispatcher for solving Partial Differential Equations.
///
/// This function attempts to solve a given PDE by trying various specialized solvers
/// based on the PDE's type, order, and provided boundary/initial conditions.
///
/// # Arguments
/// * `pde` - The PDE to solve, typically an `Expr::Eq`.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["x", "t"]`).
/// * `conditions` - An `Option` containing a slice of `Expr` representing boundary or initial conditions.
///
/// # Returns
/// An `Expr` representing the solution to the PDE, or an unevaluated `Expr::Solve` if no solution is found.
pub fn solve_pde(pde: &Expr, func: &str, vars: &[&str], conditions: Option<&[Expr]>) -> Expr {
    let equation = if let Expr::Eq(lhs, rhs) = pde {
        simplify(Expr::Sub(lhs.clone(), rhs.clone()))
    } else {
        pde.clone()
    };

    solve_pde_dispatch(&equation, func, vars, conditions)
        .unwrap_or_else(|| Expr::Solve(Box::new(pde.clone()), func.to_string()))
}

/// Internal dispatcher that attempts various solving strategies.
pub(crate) fn solve_pde_dispatch(
    equation: &Expr,
    func: &str,
    vars: &[&str],
    conditions: Option<&[Expr]>,
) -> Option<Expr> {
    if let Some(conds) = conditions {
        if let Some(solution) = solve_pde_by_separation_of_variables(equation, func, vars, conds) {
            return Some(solution);
        }
    }

    let order = get_pde_order(equation, func, vars);
    match order {
        1 => solve_pde_by_characteristics(equation, func, vars)
            .or_else(|| solve_burgers_equation(equation, func, vars, conditions)),
        2 => solve_second_order_pde(equation, func, vars)
            .or_else(|| solve_pde_by_greens_function(equation, func, vars))
            .or_else(|| solve_with_fourier_transform(equation, func, vars, conditions)),
        _ => None,
    }
}

// =====================================================================================
// endregion: Public API and Dispatchers
// =====================================================================================

// =====================================================================================
// region: Solver Implementations
// =====================================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
enum BoundaryConditionType {
    Dirichlet,
    Neumann,
}

#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    at_zero: BoundaryConditionType,
    at_l: BoundaryConditionType,
    l: Expr,
    initial_cond: Expr,
    initial_cond_deriv: Option<Expr>,
}

/// Solves 1D linear, homogeneous PDEs with homogeneous boundary conditions using Separation of Variables.
///
/// This method assumes a solution of the form `u(x,t) = X(x)T(t)` and separates the PDE
/// into two ordinary differential equations. It then applies boundary conditions to solve
/// the spatial (Sturm-Liouville) problem and initial conditions to solve the temporal problem.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["x", "t"]`).
/// * `conditions` - A slice of `Expr` representing boundary and initial conditions.
///
/// # Returns
/// An `Option<Expr>` representing the series solution, or `None` if the PDE type
/// or conditions are not supported.
pub fn solve_pde_by_separation_of_variables(
    equation: &Expr,
    func: &str,
    vars: &[&str],
    conditions: &[Expr],
) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    let x_var = vars[0];
    let t_var = vars[1];

    let bc = parse_conditions(conditions, func, x_var, t_var)?;

    let u = Expr::Variable(func.to_string());
    let u_t = differentiate(&u, t_var);
    let u_tt = differentiate(&u_t, t_var);
    let u_xx = differentiate(&differentiate(&u, x_var), x_var);

    let n = Expr::Variable("n".to_string());
    let x = Expr::Variable(x_var.to_string());
    let l = bc.l.clone();

    // 1. Solve the spatial Sturm-Liouville problem based on BCs
    let (lambda_n_sq, x_n) = match (bc.at_zero, bc.at_l) {
        (BoundaryConditionType::Dirichlet, BoundaryConditionType::Dirichlet) => {
            let lambda_n = Expr::Div(
                Box::new(Expr::Mul(Box::new(n.clone()), Box::new(Expr::Pi))),
                Box::new(l.clone()),
            );
            let x_n = Expr::Sin(Box::new(Expr::Mul(
                Box::new(lambda_n.clone()),
                Box::new(x.clone()),
            )));
            (
                Expr::Power(Box::new(lambda_n), Box::new(Expr::Constant(2.0))),
                x_n,
            )
        }
        (BoundaryConditionType::Neumann, BoundaryConditionType::Neumann) => {
            let lambda_n = Expr::Div(
                Box::new(Expr::Mul(Box::new(n.clone()), Box::new(Expr::Pi))),
                Box::new(l.clone()),
            );
            let x_n = Expr::Cos(Box::new(Expr::Mul(
                Box::new(lambda_n.clone()),
                Box::new(x.clone()),
            )));
            (
                Expr::Power(Box::new(lambda_n), Box::new(Expr::Constant(2.0))),
                x_n,
            )
        }
        (BoundaryConditionType::Dirichlet, BoundaryConditionType::Neumann) => {
            let lambda_n = Expr::Div(
                Box::new(Expr::Mul(
                    Box::new(Expr::Add(
                        Box::new(n.clone()),
                        Box::new(Expr::Constant(0.5)),
                    )),
                    Box::new(Expr::Pi),
                )),
                Box::new(l.clone()),
            );
            let x_n = Expr::Sin(Box::new(Expr::Mul(
                Box::new(lambda_n.clone()),
                Box::new(x.clone()),
            )));
            (
                Expr::Power(Box::new(lambda_n), Box::new(Expr::Constant(2.0))),
                x_n,
            )
        }
        (BoundaryConditionType::Neumann, BoundaryConditionType::Dirichlet) => {
            let lambda_n = Expr::Div(
                Box::new(Expr::Mul(
                    Box::new(Expr::Add(
                        Box::new(n.clone()),
                        Box::new(Expr::Constant(0.5)),
                    )),
                    Box::new(Expr::Pi),
                )),
                Box::new(l.clone()),
            );
            let x_n = Expr::Cos(Box::new(Expr::Mul(
                Box::new(lambda_n.clone()),
                Box::new(x.clone()),
            )));
            (
                Expr::Power(Box::new(lambda_n), Box::new(Expr::Constant(2.0))),
                x_n,
            )
        }
    };

    // 2. Identify PDE type and solve temporal part
    // Case A: Heat Equation (u_t = alpha * u_xx)
    let heat_pattern = Expr::Sub(
        Box::new(u_t.clone()),
        Box::new(Expr::Mul(
            Box::new(Expr::Pattern("alpha".to_string())),
            Box::new(u_xx.clone()),
        )),
    );
    if let Some(m) = pattern_match(equation, &heat_pattern) {
        let alpha = m.get("alpha")?;
        let t_n = Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
            Box::new(alpha.clone()),
            Box::new(Expr::Mul(
                Box::new(lambda_n_sq),
                Box::new(Expr::Variable(t_var.to_string())),
            )),
        )))));
        let cn_integrand = Expr::Mul(Box::new(bc.initial_cond.clone()), Box::new(x_n.clone()));
        let cn_integral = integrate(&cn_integrand, x_var, Some(&Expr::Constant(0.0)), Some(&l));
        let cn = Expr::Mul(
            Box::new(Expr::Div(
                Box::new(Expr::Constant(2.0)),
                Box::new(l.clone()),
            )),
            Box::new(cn_integral),
        );
        let series_term = Expr::Mul(
            Box::new(cn),
            Box::new(Expr::Mul(Box::new(t_n), Box::new(x_n))),
        );
        let solution = Expr::Summation(
            Box::new(series_term),
            "n".to_string(),
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Infinity),
        );
        return Some(Expr::Eq(Box::new(u), Box::new(solution)));
    }

    // Case B: Wave Equation (u_tt = c^2 * u_xx)
    let wave_pattern = Expr::Sub(
        Box::new(u_tt.clone()),
        Box::new(Expr::Mul(
            Box::new(Expr::Power(
                Box::new(Expr::Pattern("c".to_string())),
                Box::new(Expr::Constant(2.0)),
            )),
            Box::new(u_xx.clone()),
        )),
    );
    if let Some(m) = pattern_match(equation, &wave_pattern) {
        let c = m.get("c")?;
        let lambda_n = Expr::Sqrt(Box::new(lambda_n_sq));
        let omega_n = simplify(Expr::Mul(Box::new(c.clone()), Box::new(lambda_n)));
        let f_x = bc.initial_cond;
        let g_x = bc.initial_cond_deriv?;

        let an_integrand = Expr::Mul(Box::new(f_x.clone()), Box::new(x_n.clone()));
        let an_integral = integrate(&an_integrand, x_var, Some(&Expr::Constant(0.0)), Some(&l));
        let an = Expr::Mul(
            Box::new(Expr::Div(
                Box::new(Expr::Constant(2.0)),
                Box::new(l.clone()),
            )),
            Box::new(an_integral),
        );

        let bn_integrand = Expr::Mul(Box::new(g_x.clone()), Box::new(x_n.clone()));
        let bn_integral = integrate(&bn_integrand, x_var, Some(&Expr::Constant(0.0)), Some(&l));
        let bn = Expr::Mul(
            Box::new(Expr::Div(
                Box::new(Expr::Constant(2.0)),
                Box::new(Expr::Mul(Box::new(l.clone()), Box::new(omega_n.clone()))),
            )),
            Box::new(bn_integral),
        );

        let t_n = Expr::Add(
            Box::new(Expr::Mul(
                Box::new(an),
                Box::new(Expr::Cos(Box::new(Expr::Mul(
                    Box::new(omega_n.clone()),
                    Box::new(Expr::Variable(t_var.to_string())),
                )))),
            )),
            Box::new(Expr::Mul(
                Box::new(bn),
                Box::new(Expr::Sin(Box::new(Expr::Mul(
                    Box::new(omega_n),
                    Box::new(Expr::Variable(t_var.to_string())),
                )))),
            )),
        );
        let series_term = Expr::Mul(Box::new(t_n), Box::new(x_n));
        let solution = Expr::Summation(
            Box::new(series_term),
            "n".to_string(),
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Infinity),
        );
        return Some(Expr::Eq(Box::new(u), Box::new(solution)));
    }

    None
}

// ... [rest of the functions from previous state] ...

/// Solves first-order Partial Differential Equations using the method of characteristics.
///
/// This method transforms a PDE into a system of Ordinary Differential Equations (ODEs)
/// along characteristic curves. It is particularly effective for linear and quasi-linear
/// first-order PDEs.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["x", "y"]`).
///
/// # Returns
/// An `Option<Expr>` representing the solution, or `None` if the PDE does not match
/// a recognizable first-order linear/quasi-linear form.
pub fn solve_pde_by_characteristics(equation: &Expr, func: &str, vars: &[&str]) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    let x_var = vars[0];
    let y_var = vars[1];
    let u_func = Expr::Variable(func.to_string());
    let u_x = differentiate(&u_func, x_var);
    let u_y = differentiate(&u_func, y_var);

    let pattern_a = Expr::Mul(Box::new(Expr::Pattern("a".to_string())), Box::new(u_x));
    let pattern_b = Expr::Mul(Box::new(Expr::Pattern("b".to_string())), Box::new(u_y));
    let pattern_c = Expr::Pattern("c".to_string());
    let pattern = Expr::Sub(
        Box::new(Expr::Add(Box::new(pattern_a), Box::new(pattern_b))),
        Box::new(pattern_c),
    );

    if let Some(m) = pattern_match(equation, &pattern) {
        let a = m.get("a")?;
        let b = m.get("b")?;
        let c = m.get("c")?;

        let t = "t";
        let ode1 = Expr::Eq(
            Box::new(differentiate(&Expr::Variable(vars[0].to_string()), t)),
            Box::new(a.clone()),
        );
        let ode2 = Expr::Eq(
            Box::new(differentiate(&Expr::Variable(vars[1].to_string()), t)),
            Box::new(b.clone()),
        );
        let ode3 = Expr::Eq(
            Box::new(differentiate(&Expr::Variable(func.to_string()), t)),
            Box::new(c.clone()),
        );

        return Some(Expr::System(vec![
            Expr::Variable("Solve the following system of characteristic ODEs:".to_string()),
            ode1, ode2, ode3,
            Expr::Variable(format!("Then, find the relationship between the constants of integration to get the solution for u({},{}).", x_var, y_var))
        ]));
    }
    None
}

/// Solves a Partial Differential Equation using Green's functions.
///
/// Green's functions are used to solve inhomogeneous linear differential equations
/// with boundary conditions. The solution is expressed as an integral of the Green's
/// function multiplied by the source term.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables.
///
/// # Returns
/// An `Option<Expr>` representing the integral solution, or `None` if the differential
/// operator is not recognized or its Green's function is not implemented.
pub fn solve_pde_by_greens_function(equation: &Expr, func: &str, vars: &[&str]) -> Option<Expr> {
    let (lhs, rhs) = if let Expr::Eq(l, r) = equation {
        (&**l, &**r)
    } else {
        (equation, &Expr::Constant(0.0))
    };
    let f = if is_zero(rhs) {
        simplify(Expr::Neg(Box::new(lhs.clone())))
    } else {
        rhs.clone()
    };
    let operator_expr = if is_zero(rhs) {
        lhs.clone()
    } else {
        simplify(Expr::Sub(Box::new(lhs.clone()), Box::new(f.clone())))
    };

    let operator = identify_differential_operator(&operator_expr, func, vars);
    if operator == "Unknown_Operator" {
        return None;
    }

    let (green_function, integration_vars) = match operator.as_str() {
        "Laplacian_2D" => {
            let x_p = Expr::Variable(format!("{}_p", vars[0]));
            let y_p = Expr::Variable(format!("{}_p", vars[1]));
            let r_sq = Expr::Add(
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Variable(vars[0].to_string())),
                        Box::new(x_p.clone()),
                    )),
                    Box::new(Expr::Constant(2.0)),
                )),
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Variable(vars[1].to_string())),
                        Box::new(y_p.clone()),
                    )),
                    Box::new(Expr::Constant(2.0)),
                )),
            );
            let g = Expr::Mul(
                Box::new(Expr::Div(
                    Box::new(Expr::Constant(1.0)),
                    Box::new(Expr::Mul(Box::new(Expr::Constant(2.0)), Box::new(Expr::Pi))),
                )),
                Box::new(Expr::Log(Box::new(Expr::Sqrt(Box::new(r_sq))))),
            );
            (g, vec![x_p, y_p])
        }
        "Wave_1D" => {
            let x_p = Expr::Variable(format!("{}_p", vars[0]));
            let t_p = Expr::Variable(format!("{}_p", vars[1]));
            let c = Expr::Variable("c".to_string());
            let term = Expr::Sub(
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Variable(vars[1].to_string())),
                        Box::new(t_p.clone()),
                    )),
                    Box::new(Expr::Constant(2.0)),
                )),
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Variable(vars[0].to_string())),
                        Box::new(x_p.clone()),
                    )),
                    Box::new(Expr::Constant(2.0)),
                )),
            );
            let heaviside = Expr::Apply(Box::new(Expr::Variable("H".to_string())), Box::new(term));
            let g = Expr::Mul(
                Box::new(Expr::Div(
                    Box::new(Expr::Constant(1.0)),
                    Box::new(Expr::Mul(Box::new(Expr::Constant(2.0)), Box::new(c))),
                )),
                Box::new(heaviside),
            );
            (g, vec![x_p, t_p])
        }
        _ => return None,
    };

    let mut f_prime = f.clone();
    for (i, var) in vars.iter().enumerate() {
        f_prime = substitute(&f_prime, var, &integration_vars[i]);
    }

    let integrand = simplify(Expr::Mul(Box::new(green_function), Box::new(f_prime)));
    let mut final_integral = integrand;
    for var in integration_vars.into_iter().rev() {
        final_integral = Expr::Integral {
            integrand: Box::new(final_integral),
            var: Box::new(var),
            lower_bound: Box::new(Expr::NegativeInfinity),
            upper_bound: Box::new(Expr::Infinity),
        };
    }
    Some(final_integral)
}

/// Solves a second-order Partial Differential Equation.
///
/// This function acts as a dispatcher for various second-order PDE types.
/// It first classifies the PDE (hyperbolic, parabolic, elliptic) and then
/// attempts to apply the appropriate solution method.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables.
///
/// # Returns
/// An `Option<Expr>` representing the solution, or `None` if the PDE type
/// is not supported or cannot be solved.
pub fn solve_second_order_pde(equation: &Expr, func: &str, vars: &[&str]) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    if let Some((_a, _b, _c, pde_type)) = classify_second_order_pde(equation, func, vars) {
        match pde_type.as_str() {
            "Hyperbolic" => solve_wave_equation_1d_dalembert(equation, func, vars),
            _ => None,
        }
    } else {
        None
    }
}

/// Solves the 1D homogeneous wave equation `u_tt = c^2 * u_xx` using D'Alembert's formula.
///
/// D'Alembert's formula provides a general solution for the 1D wave equation
/// in terms of two arbitrary functions `F` and `G`:
/// `u(x,t) = F(x + ct) + G(x - ct)`.
/// Initial conditions are typically used to determine `F` and `G`.
///
/// # Arguments
/// * `equation` - The wave equation to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["t", "x"]`).
///
/// # Returns
/// An `Option<Expr>` representing the general solution, or `None` if the equation
/// does not match the 1D wave equation pattern.
pub fn solve_wave_equation_1d_dalembert(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    let t_var = vars[0];
    let x_var = vars[1];
    let u = Expr::Variable(func.to_string());

    let u_tt = differentiate(&differentiate(&u, t_var), t_var);
    let u_xx = differentiate(&differentiate(&u, x_var), x_var);

    let pattern = Expr::Sub(
        Box::new(u_tt.clone()),
        Box::new(Expr::Mul(
            Box::new(Expr::Power(
                Box::new(Expr::Pattern("c".to_string())),
                Box::new(Expr::Constant(2.0)),
            )),
            Box::new(u_xx.clone()),
        )),
    );

    if let Some(assignments) = pattern_match(equation, &pattern) {
        let c = assignments.get("c").unwrap();
        let f = Expr::Variable("F".to_string());
        let g = Expr::Variable("G".to_string());
        let x = Expr::Variable(x_var.to_string());
        let t = Expr::Variable(t_var.to_string());

        let term1 = substitute(
            &f,
            &x.to_string(),
            &simplify(Expr::Add(
                Box::new(x.clone()),
                Box::new(Expr::Mul(Box::new(c.clone()), Box::new(t.clone()))),
            )),
        );
        let term2 = substitute(
            &g,
            &x.to_string(),
            &simplify(Expr::Sub(
                Box::new(x.clone()),
                Box::new(Expr::Mul(Box::new(c.clone()), Box::new(t.clone()))),
            )),
        );

        let solution = simplify(Expr::Add(Box::new(term1), Box::new(term2)));
        return Some(Expr::Eq(Box::new(u), Box::new(solution)));
    }
    None
}

/// Solves the 1D Burgers' equation `u_t + u*u_x = 0`.
///
/// Burgers' equation is a fundamental partial differential equation occurring in various areas
/// of applied mathematics, including fluid mechanics, nonlinear acoustics, and traffic flow.
/// This function provides an implicit solution based on initial conditions.
///
/// # Arguments
/// * `equation` - The Burgers' equation to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["t", "x"]`).
/// * `initial_conditions` - An `Option` containing a slice of `Expr` representing initial conditions.
///
/// # Returns
/// An `Option<Expr>` representing the implicit solution, or `None` if the equation
/// does not match the Burgers' equation pattern or initial conditions are missing.
pub fn solve_burgers_equation(
    equation: &Expr,
    func: &str,
    vars: &[&str],
    initial_conditions: Option<&[Expr]>,
) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    let t_var = vars[0];
    let x_var = vars[1];
    let u = Expr::Variable(func.to_string());

    let u_t = differentiate(&u, t_var);
    let u_x = differentiate(&u, x_var);

    let pattern = Expr::Add(
        Box::new(u_t),
        Box::new(Expr::Mul(Box::new(u.clone()), Box::new(u_x))),
    );
    if simplify(equation.clone()) != simplify(pattern) {
        return None;
    }

    let f_of_x = initial_conditions?
        .iter()
        .find(|cond| matches!(cond, Expr::Eq(lhs, _) if **lhs == u))?;

    if let Expr::Eq(_, initial_func) = f_of_x {
        let x_minus_ut = Expr::Sub(
            Box::new(Expr::Variable(x_var.to_string())),
            Box::new(Expr::Mul(
                Box::new(u.clone()),
                Box::new(Expr::Variable(t_var.to_string())),
            )),
        );
        let implicit_solution = substitute(initial_func, x_var, &x_minus_ut);
        return Some(Expr::Eq(Box::new(u), Box::new(implicit_solution)));
    }
    None
}

/// Solves a Partial Differential Equation using the Fourier Transform method.
///
/// This method transforms the PDE from the spatial domain to the frequency domain,
/// often converting it into a simpler Ordinary Differential Equation (ODE).
/// The ODE is then solved, and the inverse Fourier Transform is applied to obtain
/// the solution in the original domain.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["t", "x"]`).
/// * `initial_conditions` - An `Option` containing a slice of `Expr` representing initial conditions.
///
/// # Returns
/// An `Option<Expr>` representing the solution, or `None` if the PDE type
/// or conditions are not supported by this method.
pub fn solve_with_fourier_transform(
    equation: &Expr,
    func: &str,
    vars: &[&str],
    initial_conditions: Option<&[Expr]>,
) -> Option<Expr> {
    if vars.len() != 2 {
        return None;
    }
    let t_var = vars[0];
    let x_var = vars[1];
    let k_var = "k";

    let initial_cond = initial_conditions?.iter().find(
        |cond| matches!(cond, Expr::Eq(lhs, _) if **lhs == Expr::Variable(func.to_string())),
    )?;
    let f_x = if let Expr::Eq(_, ic) = initial_cond {
        ic
    } else {
        return None;
    };

    let u_k_0 = transforms::fourier_transform(f_x, x_var, k_var);

    let u = Expr::Variable(func.to_string());
    let u_t = differentiate(&u, t_var);
    let u_xx = differentiate(&u_t, x_var);

    let pattern = Expr::Sub(
        Box::new(u_t),
        Box::new(Expr::Mul(
            Box::new(Expr::Pattern("alpha".to_string())),
            Box::new(u_xx),
        )),
    );
    if let Some(m) = pattern_match(equation, &pattern) {
        let alpha = m.get("alpha")?;
        let k = Expr::Variable(k_var.to_string());
        let _t = Expr::Variable(t_var.to_string());

        let neg_alpha_k_sq = Expr::Neg(Box::new(Expr::Mul(
            Box::new(alpha.clone()),
            Box::new(Expr::Power(
                Box::new(k.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
        )));
        let ode_in_t = Expr::Sub(
            Box::new(differentiate(&Expr::Variable("U".to_string()), t_var)),
            Box::new(Expr::Mul(
                Box::new(neg_alpha_k_sq),
                Box::new(Expr::Variable("U".to_string())),
            )),
        );

        let u_k_t_sol = solve_ode(&ode_in_t, "U", t_var, None);
        if let Expr::Eq(_, general_sol) = u_k_t_sol {
            let c1 = Expr::Variable("C1".to_string());
            let u_k_t = substitute(&general_sol, &c1.to_string(), &u_k_0);

            let solution = transforms::inverse_fourier_transform(&u_k_t, k_var, x_var);
            return Some(Expr::Eq(Box::new(u.clone()), Box::new(solution)));
        }
    }
    None
}

// =====================================================================================
// endregion: Core Solver Implementations
// =====================================================================================

// =====================================================================================
// region: Helpers
// =====================================================================================

pub(crate) fn get_pde_order(expr: &Expr, _func: &str, vars: &[&str]) -> usize {
    let mut max_order = 0;
    expr.pre_order_walk(&mut |sub_expr| {
        if let Expr::Derivative(inner_expr, deriv_var) = sub_expr {
            if vars.contains(&deriv_var.as_str()) {
                let mut current_order = 1;
                let mut current_inner = inner_expr.clone();
                while let Expr::Derivative(next_inner, next_deriv_var) = &*current_inner {
                    if vars.contains(&next_deriv_var.as_str()) {
                        current_order += 1;
                        current_inner = next_inner.clone();
                    } else {
                        break;
                    }
                }
                if current_order > max_order {
                    max_order = current_order;
                }
            }
        }
    });
    max_order
}

pub(crate) fn classify_second_order_pde(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<(Expr, Expr, Expr, String)> {
    let x = &vars[0];
    let y = &vars[1];
    let u = Expr::Variable(func.to_string());

    let u_x = differentiate(&u, x);
    let u_y = differentiate(&u, y);
    let u_xx = differentiate(&u_x, x);
    let u_yy = differentiate(&u_y, y);
    let u_xy = differentiate(&u_x, y);

    let (_, terms) = collect_and_order_terms(equation);
    let mut coeffs = HashMap::new();

    for (term, coeff) in &terms {
        if *term == u_xx {
            coeffs.insert("A", coeff.clone());
        } else if *term == u_xy {
            coeffs.insert("B", coeff.clone());
        } else if *term == u_yy {
            coeffs.insert("C", coeff.clone());
        }
    }

    let a = coeffs
        .get("A")
        .cloned()
        .unwrap_or_else(|| Expr::Constant(0.0));
    let b = coeffs
        .get("B")
        .cloned()
        .unwrap_or_else(|| Expr::Constant(0.0));
    let c = coeffs
        .get("C")
        .cloned()
        .unwrap_or_else(|| Expr::Constant(0.0));

    let discriminant = simplify(Expr::Sub(
        Box::new(Expr::Power(
            Box::new(b.clone()),
            Box::new(Expr::Constant(2.0)),
        )),
        Box::new(Expr::Mul(
            Box::new(Expr::Constant(4.0)),
            Box::new(Expr::Mul(Box::new(a.clone()), Box::new(c.clone()))),
        )),
    ));

    let pde_type = if let Some(d) = discriminant.to_f64() {
        if d > 0.0 {
            "Hyperbolic".to_string()
        } else if d == 0.0 {
            "Parabolic".to_string()
        } else {
            "Elliptic".to_string()
        }
    } else {
        format!("Undetermined (discriminant is symbolic: {})", discriminant)
    };

    Some((a, b, c, pde_type))
}

pub(crate) fn identify_differential_operator(lhs: &Expr, func: &str, vars: &[&str]) -> String {
    let u = Expr::Variable(func.to_string());
    if vars.len() == 2 {
        let x = vars[0];
        let y = vars[1];
        let u_xx = differentiate(&differentiate(&u, x), x);
        let u_yy = differentiate(&differentiate(&u, y), y);

        if *lhs == simplify(Expr::Add(Box::new(u_xx.clone()), Box::new(u_yy.clone()))) {
            return "Laplacian_2D".to_string();
        }

        let pattern = Expr::Sub(
            Box::new(u_yy),
            Box::new(Expr::Mul(
                Box::new(Expr::Pattern("c_sq".to_string())),
                Box::new(u_xx),
            )),
        );
        if pattern_match(lhs, &pattern).is_some() {
            return "Wave_1D".to_string();
        }
    }
    "Unknown_Operator".to_string()
}

pub(crate) fn parse_conditions(
    conditions: &[Expr],
    func: &str,
    x_var: &str,
    t_var: &str,
) -> Option<BoundaryConditions> {
    let mut at_zero: Option<BoundaryConditionType> = None;
    let mut at_l: Option<BoundaryConditionType> = None;
    let mut l: Option<Expr> = None;
    let mut initial_cond: Option<Expr> = None;
    let mut initial_cond_deriv: Option<Expr> = None;

    let u = Expr::Variable(func.to_string());
    let u_x = differentiate(&u, x_var);

    for cond in conditions {
        if let Expr::Eq(lhs, rhs) = cond {
            // Initial conditions u(x,0) or u_t(x,0)
            if let Some(val) = get_value_at_point(lhs, t_var, &Expr::Constant(0.0)) {
                if val == u {
                    initial_cond = Some(*rhs.clone());
                }
                if val == differentiate(&u, t_var) {
                    initial_cond_deriv = Some(*rhs.clone());
                }
            }
            // Boundary conditions u(0,t)=0, u(L,t)=0, u_x(0,t)=0, etc.
            else if is_zero(rhs) {
                if let Some(val) = get_value_at_point(lhs, x_var, &Expr::Constant(0.0)) {
                    if val == u {
                        at_zero = Some(BoundaryConditionType::Dirichlet);
                    }
                    if val == u_x {
                        at_zero = Some(BoundaryConditionType::Neumann);
                    }
                }
                if let Some(val) = get_value_at_point(lhs, x_var, &Expr::Variable("L".to_string()))
                {
                    if val == u {
                        at_l = Some(BoundaryConditionType::Dirichlet);
                        l = Some(Expr::Variable("L".to_string()));
                    }
                    if val == u_x {
                        at_l = Some(BoundaryConditionType::Neumann);
                        l = Some(Expr::Variable("L".to_string()));
                    }
                }
            }
        }
    }

    Some(BoundaryConditions {
        at_zero: at_zero?,
        at_l: at_l?,
        l: l?,
        initial_cond: initial_cond?,
        initial_cond_deriv,
    })
}

pub(crate) fn get_value_at_point(expr: &Expr, var: &str, point: &Expr) -> Option<Expr> {
    // A simplified way to check for expressions like `u(x,0)`
    // This needs a more robust implementation based on the expression structure.
    if let Expr::Variable(s) = expr {
        if s.contains(var) && s.contains(&point.to_string()) {
            // This is not a robust way to do this.
            // A proper implementation would parse the function call notation.
            return Some(Expr::Variable(
                s.replace(&format!("({},{})", var, point), ""),
            ));
        }
    }
    None
}

// =====================================================================================
// endregion: Helpers
// =====================================================================================
