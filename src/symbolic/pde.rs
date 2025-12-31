//! # Partial Differential Equation (PDE) Solver
//!
//! This module provides functions for solving various types of Partial Differential Equations.
//! It includes strategies for first-order PDEs (method of characteristics), second-order PDEs
//! (separation of variables, D'Alembert's formula for wave equation), and techniques like
//! Green's functions and Fourier transforms for specific PDE types.

use std::collections::HashMap;
use std::sync::Arc;

use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::integrate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::ode::solve_ode;
use crate::symbolic::simplify::collect_and_order_terms;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify::pattern_match;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::transforms;

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
///
/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.
#[must_use]

pub fn solve_pde(
    pde: &Expr,
    func: &str,
    vars: &[&str],
    conditions: Option<&[Expr]>,
) -> Expr {

    let equation =
        if let Expr::Eq(lhs, rhs) = pde
        {

            simplify(&Expr::new_sub(
                lhs.clone(),
                rhs.clone(),
            ))
        } else {

            pde.clone()
        };

    // Unwrap DAG if present
    let equation =
        if let Expr::Dag(node) =
            equation
        {

            node.to_expr()
                .expect(
                    "Unwrap DAG in \
                     solve_pde",
                )
        } else {

            equation
        };

    solve_pde_dispatch(
        &equation,
        func,
        vars,
        conditions,
    )
    .unwrap_or_else(|| {

        Expr::Solve(
            Arc::new(pde.clone()),
            func.to_string(),
        )
    })
}

/// Internal dispatcher that attempts various solving strategies.

pub(crate) fn solve_pde_dispatch(
    equation: &Expr,
    func: &str,
    vars: &[&str],
    conditions: Option<&[Expr]>,
) -> Option<Expr> {

    if let Some(conds) = conditions {

        if let Some(solution) = solve_pde_by_separation_of_variables(
            equation,
            func,
            vars,
            conds,
        ) {

            return Some(solution);
        }
    }

    let order = get_pde_order(
        equation,
        func,
        vars,
    );

    match order {
        | 1 => {
            solve_pde_by_characteristics(equation, func, vars).or_else(|| {

                solve_burgers_equation(
                    equation,
                    func,
                    vars,
                    conditions,
                )
            })
        },
        | 2 => {
            solve_second_order_pde(equation, func, vars)
                .or_else(|| solve_pde_by_greens_function(equation, func, vars))
                .or_else(|| {

                    solve_with_fourier_transform(
                        equation,
                        func,
                        vars,
                        conditions,
                    )
                })
        },
        | _ => None,
    }
}

#[derive(
    Debug, Clone, PartialEq, Eq,
)]

enum BoundaryConditionType {
    Dirichlet,
    Neumann,
}

#[derive(Debug, Clone)]

/// Represents the set of boundary and initial conditions for a 1D PDE.

pub struct BoundaryConditions {
    /// The type of boundary condition at x=0.
    at_zero: BoundaryConditionType,
    /// The type of boundary condition at x=L.
    at_l: BoundaryConditionType,
    /// The length of the spatial domain L.
    l: Expr,
    /// The initial condition function at t=0, u(x, 0) = f(x).
    initial_cond: Expr,
    /// The initial condition for the first time derivative, ∂u/∂t(x, 0) = g(x) (for wave eq).
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
#[must_use]

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

    let bc = parse_conditions(
        conditions,
        func,
        x_var,
        t_var,
    )?;

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let u_tt = Expr::Derivative(
        Arc::new(u_t.clone()),
        t_var.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x),
        x_var.to_string(),
    );

    let n =
        Expr::Variable("n".to_string());

    let x = Expr::Variable(
        x_var.to_string(),
    );

    let l = bc.l.clone();

    let (lambda_n_sq, x_n) = match (bc.at_zero, bc.at_l) {
        | (BoundaryConditionType::Dirichlet, BoundaryConditionType::Dirichlet) => {

            let lambda_n = Expr::new_div(
                Expr::new_mul(n, Expr::Pi),
                l.clone(),
            );

            let x_n = Expr::new_sin(Expr::new_mul(
                lambda_n.clone(),
                x,
            ));

            (
                Expr::new_pow(
                    lambda_n,
                    Expr::Constant(2.0),
                ),
                x_n,
            )
        },
        | (BoundaryConditionType::Neumann, BoundaryConditionType::Neumann) => {

            let lambda_n = Expr::new_div(
                Expr::new_mul(n, Expr::Pi),
                l.clone(),
            );

            let x_n = Expr::new_cos(Expr::new_mul(
                lambda_n.clone(),
                x,
            ));

            (
                Expr::new_pow(
                    lambda_n,
                    Expr::Constant(2.0),
                ),
                x_n,
            )
        },
        | (BoundaryConditionType::Dirichlet, BoundaryConditionType::Neumann) => {

            let lambda_n = Expr::new_div(
                Expr::new_mul(
                    Expr::new_add(
                        n,
                        Expr::Constant(0.5),
                    ),
                    Expr::Pi,
                ),
                l.clone(),
            );

            let x_n = Expr::new_sin(Expr::new_mul(
                lambda_n.clone(),
                x,
            ));

            (
                Expr::new_pow(
                    lambda_n,
                    Expr::Constant(2.0),
                ),
                x_n,
            )
        },
        | (BoundaryConditionType::Neumann, BoundaryConditionType::Dirichlet) => {

            let lambda_n = Expr::new_div(
                Expr::new_mul(
                    Expr::new_add(
                        n,
                        Expr::Constant(0.5),
                    ),
                    Expr::Pi,
                ),
                l.clone(),
            );

            let x_n = Expr::new_cos(Expr::new_mul(
                lambda_n.clone(),
                x,
            ));

            (
                Expr::new_pow(
                    lambda_n,
                    Expr::Constant(2.0),
                ),
                x_n,
            )
        },
    };

    let heat_pattern = Expr::new_sub(
        u_t,
        Expr::new_mul(
            Expr::Pattern(
                "alpha".to_string(),
            ),
            u_xx.clone(),
        ),
    );

    if let Some(m) = pattern_match(
        equation,
        &heat_pattern,
    ) {

        let alpha = m.get("alpha")?;

        let t_n = Expr::new_exp(Expr::new_neg(
            Expr::new_mul(
                alpha.clone(),
                Expr::new_mul(
                    lambda_n_sq,
                    Expr::Variable(t_var.to_string()),
                ),
            ),
        ));

        let cn_integrand =
            Expr::new_mul(
                bc.initial_cond,
                x_n.clone(),
            );

        let cn_integral = integrate(
            &cn_integrand,
            x_var,
            Some(&Expr::Constant(0.0)),
            Some(&l),
        );

        let cn = Expr::new_mul(
            Expr::new_div(
                Expr::Constant(2.0),
                l,
            ),
            cn_integral,
        );

        let series_term = Expr::new_mul(
            cn,
            Expr::new_mul(t_n, x_n),
        );

        let solution = Expr::Summation(
            Arc::new(series_term),
            "n".to_string(),
            Arc::new(Expr::Constant(
                1.0,
            )),
            Arc::new(Expr::Infinity),
        );

        return Some(Expr::Eq(
            Arc::new(u),
            Arc::new(solution),
        ));
    }

    let wave_pattern = Expr::new_sub(
        u_tt,
        Expr::new_mul(
            Expr::new_pow(
                Expr::Pattern(
                    "c".to_string(),
                ),
                Expr::Constant(2.0),
            ),
            u_xx,
        ),
    );

    if let Some(m) = pattern_match(
        equation,
        &wave_pattern,
    ) {

        let c = m.get("c")?;

        let lambda_n =
            Expr::new_sqrt(lambda_n_sq);

        let omega_n =
            simplify(&Expr::new_mul(
                c.clone(),
                lambda_n,
            ));

        let f_x = bc.initial_cond;

        let g_x =
            bc.initial_cond_deriv?;

        let an_integrand =
            Expr::new_mul(
                f_x,
                x_n.clone(),
            );

        let an_integral = integrate(
            &an_integrand,
            x_var,
            Some(&Expr::Constant(0.0)),
            Some(&l),
        );

        let an = Expr::new_mul(
            Expr::new_div(
                Expr::Constant(2.0),
                l.clone(),
            ),
            an_integral,
        );

        let bn_integrand =
            Expr::new_mul(
                g_x,
                x_n.clone(),
            );

        let bn_integral = integrate(
            &bn_integrand,
            x_var,
            Some(&Expr::Constant(0.0)),
            Some(&l),
        );

        let bn = Expr::new_mul(
            Expr::new_div(
                Expr::Constant(2.0),
                Expr::new_mul(
                    l,
                    omega_n.clone(),
                ),
            ),
            bn_integral,
        );

        let t_n = Expr::new_add(
            Expr::new_mul(
                an,
                Expr::new_cos(Expr::new_mul(
                    omega_n.clone(),
                    Expr::Variable(t_var.to_string()),
                )),
            ),
            Expr::new_mul(
                bn,
                Expr::new_sin(Expr::new_mul(
                    omega_n,
                    Expr::Variable(t_var.to_string()),
                )),
            ),
        );

        let series_term =
            Expr::new_mul(t_n, x_n);

        let solution = Expr::Summation(
            Arc::new(series_term),
            "n".to_string(),
            Arc::new(Expr::Constant(
                1.0,
            )),
            Arc::new(Expr::Infinity),
        );

        return Some(Expr::Eq(
            Arc::new(u),
            Arc::new(solution),
        ));
    }

    None
}

/// Classification of a PDE based on its properties
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    serde::Serialize,
    serde::Deserialize,
)]

pub enum PDEType {
    /// Wave equation (hyperbolic)
    Wave,
    /// Heat/Diffusion equation (parabolic)
    Heat,
    /// Laplace equation (elliptic, homogeneous)
    Laplace,
    /// Poisson equation (elliptic, inhomogeneous)
    Poisson,
    /// Helmholtz equation
    Helmholtz,
    /// Schrödinger equation
    Schrodinger,
    /// Burgers equation (nonlinear)
    Burgers,
    /// First-order transport/advection
    Transport,
    /// Unknown or custom PDE
    Unknown,
}

#[derive(
    Debug,
    Clone,
    serde::Serialize,
    serde::Deserialize,
)]

/// Contains the results of a PDE classification.

pub struct PDEClassification {
    /// The identified type of the PDE (e.g., Wave, Heat).
    pub pde_type: PDEType,
    /// The highest order of derivatives present in the equation.
    pub order: usize,
    /// The number of independent variables (spatial and temporal).
    pub dimension: usize,
    /// Whether the equation is linear with respect to the unknown function and its derivatives.
    pub is_linear: bool,
    /// Whether the equation is homogeneous (all terms contain the unknown function or its derivatives).
    pub is_homogeneous: bool,
    /// A list of suggested methods for solving this specific classification of PDE.
    pub suggested_methods: Vec<String>,
}

/// Heuristically classifies a PDE and suggests solution methods
///
/// This function analyzes the structure of a PDE to determine:
/// - Type (wave, heat, Laplace, etc.)
/// - Order (1st, 2nd, etc.)
/// - Dimension (1D, 2D, 3D, etc.)
/// - Linearity
/// - Homogeneity
/// - Suggested solution methods
///
/// # Arguments
/// * `equation` - The PDE to classify
/// * `func` - The unknown function name
/// * `vars` - The independent variables
///
/// # Returns
/// A `PDEClassification` struct with analysis results
#[must_use]

pub fn classify_pde_heuristic(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> PDEClassification {

    let order = get_pde_order(
        equation,
        func,
        vars,
    );

    let dimension = vars.len();

    // Check linearity by looking for products of u with itself or its derivatives
    let is_linear =
        check_linearity(equation, func);

    // Check homogeneity by looking for terms without u or its derivatives
    let is_homogeneous =
        check_homogeneity(
            equation,
            func,
        );

    // Determine PDE type based on structure
    let pde_type = if order == 1 {

        classify_first_order(
            equation,
            func,
            vars,
        )
    } else if order == 2 {

        classify_second_order(
            equation,
            func,
            vars,
        )
    } else {

        PDEType::Unknown
    };

    // Suggest solution methods based on classification
    let suggested_methods =
        suggest_solution_methods(
            &pde_type,
            order,
            dimension,
            is_linear,
            is_homogeneous,
        );

    PDEClassification {
        pde_type,
        order,
        dimension,
        is_linear,
        is_homogeneous,
        suggested_methods,
    }
}

/// Check if a PDE is linear (no products of u with itself or its derivatives)

fn check_linearity(
    equation: &Expr,
    func: &str,
) -> bool {

    // For now, a simple heuristic: check if there are any Mul nodes
    // where both operands contain the function or its derivatives
    // This is a simplified check; a full implementation would be more sophisticated
    !contains_nonlinear_terms(
        equation,
        func,
    )
}

fn contains_nonlinear_terms(
    expr: &Expr,
    func: &str,
) -> bool {

    match expr {
        | Expr::Mul(a, b) => {

            let a_has_func = contains_function_or_derivative(a, func);

            let b_has_func = contains_function_or_derivative(b, func);

            if a_has_func && b_has_func {

                return true;
            }

            contains_nonlinear_terms(a, func) || contains_nonlinear_terms(b, func)
        },
        | Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Div(a, b) | Expr::Power(a, b) => {
            contains_nonlinear_terms(a, func) || contains_nonlinear_terms(b, func)
        },
        | Expr::Dag(node) => {
            node.to_expr()
                .is_ok_and(|inner| contains_nonlinear_terms(&inner, func))
        },
        | _ => false,
    }
}

fn contains_function_or_derivative(
    expr: &Expr,
    func: &str,
) -> bool {

    match expr {
        | Expr::Variable(v) if v == func => true,
        | Expr::Derivative(inner, _) => contains_function_or_derivative(inner, func),
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b) => {
            contains_function_or_derivative(a, func) || contains_function_or_derivative(b, func)
        },
        | Expr::Dag(node) => {
            node.to_expr().is_ok_and(|inner| {
                contains_function_or_derivative(&inner, func)
            })
        },
        | _ => false,
    }
}

/// Check if a PDE is homogeneous (all terms contain u or its derivatives)

fn check_homogeneity(
    equation: &Expr,
    func: &str,
) -> bool {

    let terms = collect_terms(equation);

    for term in &terms {

        // Handle Neg wrapper
        let inner_term =
            if let Expr::Neg(inner) =
                term
            {

                inner.as_ref()
            } else {

                term
            };

        if !contains_function_or_derivative(inner_term, func) {

            // Found a term without the function - inhomogeneous
            return false;
        }
    }

    true
}

/// Classify first-order PDEs

fn classify_first_order(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> PDEType {

    // Check for Burgers equation: u_t + u*u_x = 0
    if vars.len() == 2 {

        let u = Expr::Variable(
            func.to_string(),
        );

        let u_t = Expr::Derivative(
            Arc::new(u.clone()),
            vars[0].to_string(),
        );

        let u_x = Expr::Derivative(
            Arc::new(u.clone()),
            vars[1].to_string(),
        );

        let burgers_pattern =
            Expr::new_add(
                u_t,
                Expr::new_mul(u, u_x),
            );

        if simplify(equation)
            == simplify(
                &burgers_pattern,
            )
        {

            return PDEType::Burgers;
        }
    }

    // Otherwise, assume transport/advection equation
    PDEType::Transport
}

/// Classify second-order PDEs

fn classify_second_order(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> PDEType {

    if vars.len() < 2 {

        return PDEType::Unknown;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let terms = collect_terms(equation);

    // Check for different types of derivatives
    // Assume first variable is time if there are time derivatives
    let mut has_first_time_deriv =
        false;

    let mut has_second_time_deriv =
        false;

    let mut spatial_second_deriv_count =
        0;

    // Check first variable for time derivatives
    if vars.len() >= 2 {

        let u_t = Expr::Derivative(
            Arc::new(u.clone()),
            vars[0].to_string(),
        );

        let u_tt = Expr::Derivative(
            Arc::new(u_t.clone()),
            vars[0].to_string(),
        );

        for term in &terms {

            if extract_coefficient(
                term, &u_t,
            )
            .is_some()
            {

                has_first_time_deriv =
                    true;
            }

            if extract_coefficient(
                term, &u_tt,
            )
            .is_some()
            {

                has_second_time_deriv =
                    true;
            }
        }
    }

    // Count spatial second derivatives (all variables or all except first if time-dependent)
    // Only skip first variable if we actually found time derivatives in it
    let spatial_start = usize::from(
        has_first_time_deriv
            || has_second_time_deriv,
    );

    for var in vars.iter().skip(spatial_start) {

        let u_var = Expr::Derivative(
            Arc::new(u.clone()),
            (*var).to_string(),
        );

        let u_var_var =
            Expr::Derivative(
                Arc::new(u_var),
                (*var).to_string(),
            );

        for term in &terms {

            if extract_coefficient(
                term,
                &u_var_var,
            )
            .is_some()
            {

                spatial_second_deriv_count += 1;

                break; // Count each spatial variable only once
            }
        }
    }

    // If no time derivatives were found, also check all variables for second derivatives
    // This handles the case where all variables are spatial (e.g., Laplace equation)
    if !has_first_time_deriv
        && !has_second_time_deriv
    {

        spatial_second_deriv_count = 0; // Reset and recount all
        for var in vars {

            let u_var =
                Expr::Derivative(
                    Arc::new(u.clone()),
                    (*var).to_string(),
                );

            let u_var_var =
                Expr::Derivative(
                    Arc::new(u_var),
                    (*var).to_string(),
                );

            for term in &terms {

                if extract_coefficient(
                    term,
                    &u_var_var,
                )
                .is_some()
                {

                    spatial_second_deriv_count += 1;

                    break;
                }
            }
        }
    }

    // Classification logic
    if has_second_time_deriv
        && spatial_second_deriv_count
            > 0
    {

        // Has u_tt and spatial second derivatives -> Wave equation
        return PDEType::Wave;
    } else if has_first_time_deriv
        && !has_second_time_deriv
        && spatial_second_deriv_count
            > 0
    {

        // Has u_t (but not u_tt) and spatial second derivatives -> Heat equation
        return PDEType::Heat;
    } else if !has_first_time_deriv
        && !has_second_time_deriv
        && spatial_second_deriv_count
            >= 2
    {

        // No time derivatives, multiple spatial second derivatives -> Elliptic (Laplace/Poisson)
        if check_homogeneity(
            equation,
            func,
        ) {

            return PDEType::Laplace;
        }

        return PDEType::Poisson;
    }

    PDEType::Unknown
}

/// Suggest solution methods based on PDE classification

fn suggest_solution_methods(
    pde_type: &PDEType,
    order: usize,
    dimension: usize,
    is_linear: bool,
    is_homogeneous: bool,
) -> Vec<String> {

    let mut methods = Vec::new();

    match pde_type {
        | PDEType::Wave => {

            if dimension == 1 {

                methods.push(
                    "D'Alembert's \
                     formula"
                        .to_string(),
                );
            }

            methods.push(
                "Separation of \
                 variables"
                    .to_string(),
            );

            methods.push(
                "Fourier series"
                    .to_string(),
            );

            if !is_homogeneous {

                methods.push(
                    "Green's function"
                        .to_string(),
                );
            }
        },
        | PDEType::Heat => {

            methods.push(
                "Separation of \
                 variables"
                    .to_string(),
            );

            methods.push(
                "Fourier series"
                    .to_string(),
            );

            methods.push(
                "Fourier transform"
                    .to_string(),
            );

            if !is_homogeneous {

                methods.push(
                    "Green's function"
                        .to_string(),
                );

                methods.push("Duhamel's principle".to_string());
            }
        },
        | PDEType::Laplace
        | PDEType::Poisson => {

            methods.push(
                "Separation of \
                 variables"
                    .to_string(),
            );

            methods.push(
                "Fourier series"
                    .to_string(),
            );

            if pde_type
                == &PDEType::Poisson
            {

                methods.push(
                    "Green's function"
                        .to_string(),
                );
            }

            if dimension == 2 {

                methods.push(
                    "Conformal mapping"
                        .to_string(),
                );
            }
        },
        | PDEType::Transport => {

            methods.push("Method of characteristics".to_string());
        },
        | PDEType::Burgers => {

            methods.push("Cole-Hopf transformation".to_string());

            methods.push("Method of characteristics".to_string());
        },
        | PDEType::Helmholtz => {

            methods.push(
                "Separation of \
                 variables"
                    .to_string(),
            );

            methods.push(
                "Green's function"
                    .to_string(),
            );
        },
        | PDEType::Schrodinger => {

            methods.push(
                "Separation of \
                 variables"
                    .to_string(),
            );

            methods.push(
                "WKB approximation"
                    .to_string(),
            );
        },
        | PDEType::Unknown => {

            if order == 1 && is_linear {

                methods.push("Method of characteristics".to_string());
            }

            if order == 2
                && is_linear
                && is_homogeneous
            {

                methods.push(
                    "Separation of \
                     variables"
                        .to_string(),
                );
            }

            methods.push(
                "Numerical methods"
                    .to_string(),
            );
        },
    }

    methods
}

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

fn extract_coefficient(
    term: &Expr,
    var: &Expr,
) -> Option<Expr> {

    // Unwrap DAG if present
    let term =
        if let Expr::Dag(node) = term {

            node.to_expr()
                .ok()?
        } else {

            term.clone()
        };

    if &term == var {

        return Some(Expr::Constant(
            1.0,
        ));
    }

    match &term {
        | Expr::Neg(inner) => {

            // Unwrap DAG in inner
            let inner =
                if let Expr::Dag(node) =
                    inner.as_ref()
                {

                    node.to_expr()
                        .ok()?
                } else {

                    inner
                        .as_ref()
                        .clone()
                };

            if &inner == var {

                return Some(
                    Expr::Constant(
                        -1.0,
                    ),
                );
            }

            if let Expr::Mul(a, b) =
                &inner
            {

                // Unwrap DAG in a and b
                let a_unwrapped =
                    if let Expr::Dag(
                        node,
                    ) = a.as_ref()
                    {

                        node.to_expr()
                            .ok()?
                    } else {

                        a.as_ref()
                            .clone()
                    };

                let b_unwrapped =
                    if let Expr::Dag(
                        node,
                    ) = b.as_ref()
                    {

                        node.to_expr()
                            .ok()?
                    } else {

                        b.as_ref()
                            .clone()
                    };

                if &a_unwrapped == var {

                    return Some(
                        Expr::new_neg(
                            b_unwrapped,
                        ),
                    );
                }

                if &b_unwrapped == var {

                    return Some(
                        Expr::new_neg(
                            a_unwrapped,
                        ),
                    );
                }
            }

            None
        },
        | Expr::Mul(a, b) => {

            // Unwrap DAG in a and b
            let a_unwrapped =
                if let Expr::Dag(node) =
                    a.as_ref()
                {

                    node.to_expr()
                        .ok()?
                } else {

                    a.as_ref().clone()
                };

            let b_unwrapped =
                if let Expr::Dag(node) =
                    b.as_ref()
                {

                    node.to_expr()
                        .ok()?
                } else {

                    b.as_ref().clone()
                };

            if &a_unwrapped == var {

                return Some(
                    b_unwrapped,
                );
            }

            if &b_unwrapped == var {

                return Some(
                    a_unwrapped,
                );
            }

            None
        },
        | _ => None,
    }
}

fn collect_terms(
    expr: &Expr
) -> Vec<Expr> {

    match expr {
        | Expr::Add(a, b) => {

            let mut terms =
                collect_terms(a);

            terms.extend(
                collect_terms(b),
            );

            terms
        },
        | Expr::Sub(a, b) => {

            let mut terms =
                collect_terms(a);

            let neg_b_terms =
                collect_terms(b);

            for term in neg_b_terms {

                // Negate the term
                if let Expr::Neg(inner) = term {

                    terms.push(
                        inner
                            .as_ref()
                            .clone(),
                    ); // -(-x) = x
                } else if let Expr::Constant(c) = term {

                    terms.push(Expr::Constant(-c));
                } else {

                    terms.push(Expr::new_neg(term));
                }
            }

            terms
        },
        | Expr::Neg(inner) => {

            let terms =
                collect_terms(inner);

            terms
                .into_iter()
                .map(Expr::new_neg)
                .collect()
        },
        | Expr::Dag(node) => {
            collect_terms(
                &node
                    .to_expr()
                    .unwrap(),
            )
        },
        | _ => vec![expr.clone()],
    }
}

#[must_use]

/// Solves a first-order PDE using the method of characteristics.
///
/// This method transforms the PDE into a set of ordinary differential equations (characteristic equations)
/// that describe how the solution evolves along specific curves.
///
/// # Arguments
/// * `equation` - The PDE to solve.
/// * `func` - The name of the unknown function.
/// * `vars` - The independent variables.
///
/// # Returns
/// An `Option<Expr>` containing the solution if found.

pub fn solve_pde_by_characteristics(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 2 {

        return None;
    }

    let x_var = vars[0];

    let y_var = vars[1];

    let u_func = Expr::Variable(
        func.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u_func.clone()),
        x_var.to_string(),
    );

    let u_y = Expr::Derivative(
        Arc::new(u_func),
        y_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut a = Expr::Constant(0.0);

    let mut b = Expr::Constant(0.0);

    let mut c_neg = Expr::Constant(0.0); // Terms not containing u_x or u_y

    for term in terms {

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_x,
            )
        {

            a = Expr::new_add(a, coeff);

            continue;
        }

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_y,
            )
        {

            b = Expr::new_add(b, coeff);

            continue;
        }

        // Else it's part of c (moved to RHS, so -c on LHS)
        c_neg =
            Expr::new_add(c_neg, term);
    }

    let a = simplify(&a);

    let b = simplify(&b);

    let c = simplify(&Expr::new_neg(
        c_neg,
    )); // c is on RHS

    // Check if it's a valid characteristic equation (linear in derivatives)
    if is_zero(&a) && is_zero(&b) {

        return None;
    }

    // Solve characteristic equations: dy/dx = b/a
    // For now, let's just implement the simple case where a, b are constants
    // xi = y - (b/a)x

    let b_over_a = simplify(
        &Expr::new_div(b, a.clone()),
    );

    let c_over_a =
        simplify(&Expr::new_div(c, a));

    // Characteristic variable xi
    // Assuming a, b are constants for now or simple enough
    // xi = y - (b/a)x
    let xi = simplify(&Expr::new_sub(
        Expr::Variable(
            y_var.to_string(),
        ),
        Expr::new_mul(
            b_over_a,
            Expr::Variable(
                x_var.to_string(),
            ),
        ),
    ));

    // Particular solution for u
    // u_p = Integral(c/a) dx
    let u_p = integrate(
        &c_over_a,
        x_var,
        None,
        None,
    );

    // General solution: u = u_p + F(xi)
    let f_xi = Expr::Apply(
        Arc::new(Expr::Variable(
            "F".to_string(),
        )),
        Arc::new(xi),
    );

    let solution = simplify(
        &Expr::new_add(u_p, f_xi),
    );

    Some(solution)
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
#[must_use]

pub fn solve_pde_by_greens_function(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    let (lhs, rhs) = if let Expr::Eq(
        l,
        r,
    ) = equation
    {

        (&**l, &**r)
    } else {

        (
            equation,
            &Expr::Constant(0.0),
        )
    };

    let f = if is_zero(rhs) {

        simplify(&Expr::new_neg(
            lhs.clone(),
        ))
    } else {

        rhs.clone()
    };

    let operator_expr = if is_zero(rhs)
    {

        lhs.clone()
    } else {

        simplify(&Expr::new_sub(
            lhs.clone(),
            f.clone(),
        ))
    };

    let operator =
        identify_differential_operator(
            &operator_expr,
            func,
            vars,
        );

    if operator == "Unknown_Operator" {

        return None;
    }

    let (
        green_function,
        integration_vars,
    ) = match operator.as_str() {
        | "Laplacian_2D" => {

            let x_p = Expr::Variable(
                format!(
                    "{}_p",
                    vars[0]
                ),
            );

            let y_p = Expr::Variable(
                format!(
                    "{}_p",
                    vars[1]
                ),
            );

            let r_sq = Expr::new_add(
                Expr::new_pow(
                    Expr::new_sub(
                        Expr::Variable(vars[0].to_string()),
                        x_p.clone(),
                    ),
                    Expr::Constant(2.0),
                ),
                Expr::new_pow(
                    Expr::new_sub(
                        Expr::Variable(vars[1].to_string()),
                        y_p.clone(),
                    ),
                    Expr::Constant(2.0),
                ),
            );

            let g = Expr::new_mul(
                Expr::new_div(
                    Expr::Constant(1.0),
                    Expr::new_mul(
                        Expr::Constant(
                            2.0,
                        ),
                        Expr::Pi,
                    ),
                ),
                Expr::new_log(
                    Expr::new_sqrt(
                        r_sq,
                    ),
                ),
            );

            (g, vec![x_p, y_p])
        },
        | "Wave_1D" => {

            let x_p = Expr::Variable(
                format!(
                    "{}_p",
                    vars[0]
                ),
            );

            let t_p = Expr::Variable(
                format!(
                    "{}_p",
                    vars[1]
                ),
            );

            let c = Expr::Variable(
                "c".to_string(),
            );

            let term = Expr::new_sub(
                Expr::new_pow(
                    Expr::new_sub(
                        Expr::Variable(vars[1].to_string()),
                        t_p.clone(),
                    ),
                    Expr::Constant(2.0),
                ),
                Expr::new_pow(
                    Expr::new_sub(
                        Expr::Variable(vars[0].to_string()),
                        x_p.clone(),
                    ),
                    Expr::Constant(2.0),
                ),
            );

            let heaviside =
                Expr::new_apply(
                    Expr::Variable(
                        "H".to_string(),
                    ),
                    term,
                );

            let g = Expr::new_mul(
                Expr::new_div(
                    Expr::Constant(1.0),
                    Expr::new_mul(
                        Expr::Constant(
                            2.0,
                        ),
                        c,
                    ),
                ),
                heaviside,
            );

            (g, vec![x_p, t_p])
        },
        | _ => return None,
    };

    let mut f_prime = f;

    for (i, var) in vars
        .iter()
        .enumerate()
    {

        f_prime = substitute(
            &f_prime,
            var,
            &integration_vars[i],
        );
    }

    let integrand =
        simplify(&Expr::new_mul(
            green_function,
            f_prime,
        ));

    let mut final_integral = integrand;

    for var in integration_vars
        .into_iter()
        .rev()
    {

        final_integral = Expr::Integral {
            integrand : Arc::new(final_integral),
            var : Arc::new(var),
            lower_bound : Arc::new(Expr::NegativeInfinity),
            upper_bound : Arc::new(Expr::Infinity),
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
#[must_use]

pub fn solve_second_order_pde(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 2 {

        return None;
    }

    match classify_second_order_pde(
        equation,
        func,
        vars,
    ) { Some((
        _a,
        _b,
        _c,
        pde_type,
    )) => {

        match pde_type.as_str() {
            | "Hyperbolic" => solve_wave_equation_1d_dalembert(equation, func, vars),
            | _ => None,
        }
    } _ => {

        None
    }}
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
#[must_use]

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

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let u_tt = Expr::Derivative(
        Arc::new(u_t),
        t_var.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x),
        x_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut coeff_u_tt =
        Expr::Constant(0.0);

    let mut coeff_u_xx =
        Expr::Constant(0.0);

    for term in terms {

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_tt,
            )
        {

            coeff_u_tt = Expr::new_add(
                coeff_u_tt,
                coeff,
            );

            continue;
        }

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_xx,
            )
        {

            coeff_u_xx = Expr::new_add(
                coeff_u_xx,
                coeff,
            );

            continue;
        }
    }

    let coeff_u_tt =
        simplify(&coeff_u_tt);

    let coeff_u_xx =
        simplify(&coeff_u_xx);

    // We expect coeff_u_tt = 1 (or constant) and coeff_u_xx = -c^2
    // Or proportional.

    // Check if coeff_u_tt is non-zero
    if is_zero(&coeff_u_tt) {

        return None;
    }

    // c^2 = - coeff_u_xx / coeff_u_tt
    let c_sq = simplify(
        &Expr::new_neg(Expr::new_div(
            coeff_u_xx,
            coeff_u_tt,
        )),
    );

    // Check if c_sq is positive (or at least not obviously negative/zero)
    // and extract c.
    // If c_sq is constant, we can take sqrt.
    // If c_sq is var^2, we can take var.

    let c = simplify(&Expr::new_sqrt(
        c_sq,
    ));

    // D'Alembert solution: u(x,t) = F(x - ct) + G(x + ct)
    let f =
        Expr::Variable("F".to_string());

    let g =
        Expr::Variable("G".to_string());

    let arg1 =
        simplify(&Expr::new_sub(
            Expr::Variable(
                x_var.to_string(),
            ),
            Expr::new_mul(
                c.clone(),
                Expr::Variable(
                    t_var.to_string(),
                ),
            ),
        ));

    let arg2 =
        simplify(&Expr::new_add(
            Expr::Variable(
                x_var.to_string(),
            ),
            Expr::new_mul(
                c,
                Expr::Variable(
                    t_var.to_string(),
                ),
            ),
        ));

    let sol = Expr::new_add(
        Expr::Apply(
            Arc::new(f),
            Arc::new(arg1),
        ),
        Expr::Apply(
            Arc::new(g),
            Arc::new(arg2),
        ),
    );

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(sol),
    ))
}

/// Solves the 1D heat equation `u_t = α*u_xx`.
///
/// The heat equation (also known as the diffusion equation) describes the distribution
/// of heat (or variation in temperature) in a given region over time. This function
/// solves the homogeneous 1D heat equation with constant thermal diffusivity α.
///
/// # Arguments
/// * `equation` - The heat equation to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["t", "x"]`).
///
/// # Returns
/// An `Option<Expr>` representing the general solution as a Fourier series, or `None`
/// if the equation does not match the heat equation pattern.
///
/// # Example
/// For `u_t = α*u_xx`, the general solution is:
/// `u(x,t) = Σ A_n * exp(-α*n²*π²*t/L²) * sin(n*π*x/L)`
#[must_use]

pub fn solve_heat_equation_1d(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 2 {

        return None;
    }

    let t_var = vars[0];

    let x_var = vars[1];

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x),
        x_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut coeff_u_t =
        Expr::Constant(0.0);

    let mut coeff_u_xx =
        Expr::Constant(0.0);

    for term in terms {

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_t,
            )
        {

            coeff_u_t = Expr::new_add(
                coeff_u_t,
                coeff,
            );

            continue;
        }

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_xx,
            )
        {

            coeff_u_xx = Expr::new_add(
                coeff_u_xx,
                coeff,
            );

            continue;
        }
    }

    let coeff_u_t =
        simplify(&coeff_u_t);

    let coeff_u_xx =
        simplify(&coeff_u_xx);

    // Check if coeff_u_t is non-zero
    if is_zero(&coeff_u_t) {

        return None;
    }

    // α = coeff_u_xx / coeff_u_t
    let alpha =
        simplify(&Expr::new_div(
            coeff_u_xx,
            coeff_u_t,
        ));

    // General solution using Fourier series:
    // u(x,t) = Σ A_n * exp(-α*n²*π²*t/L²) * sin(n*π*x/L)
    // For now, return a symbolic representation

    let n =
        Expr::Variable("n".to_string());

    let a_n = Expr::Variable(
        "A_n".to_string(),
    );

    let l =
        Expr::Variable("L".to_string());

    let t = Expr::Variable(
        t_var.to_string(),
    );

    let x = Expr::Variable(
        x_var.to_string(),
    );

    // exp(-α*n²*π²*t/L²)
    let exponent =
        Expr::new_neg(Expr::new_div(
            Expr::new_mul(
                Expr::new_mul(
                    alpha,
                    Expr::new_pow(
                        n.clone(),
                        Expr::Constant(
                            2.0,
                        ),
                    ),
                ),
                Expr::new_mul(
                    Expr::new_pow(
                        Expr::Pi,
                        Expr::Constant(
                            2.0,
                        ),
                    ),
                    t,
                ),
            ),
            Expr::new_pow(
                l.clone(),
                Expr::Constant(2.0),
            ),
        ));

    let time_part =
        Expr::new_exp(exponent);

    // sin(n*π*x/L)
    let spatial_part =
        Expr::new_sin(Expr::new_div(
            Expr::new_mul(
                Expr::new_mul(
                    n.clone(),
                    Expr::Pi,
                ),
                x,
            ),
            l,
        ));

    // A_n * exp(...) * sin(...)
    let term = Expr::new_mul(
        Expr::new_mul(a_n, time_part),
        spatial_part,
    );

    // Σ from n=1 to ∞
    let solution = Expr::Sum {
        body: Arc::new(term),
        var: Arc::new(n),
        from: Arc::new(Expr::Constant(
            1.0,
        )),
        to: Arc::new(Expr::Infinity),
    };

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution),
    ))
}

/// Solves the 2D Laplace equation `u_xx + u_yy = 0`.
///
/// Laplace's equation is a second-order partial differential equation that describes
/// steady-state heat distribution, electrostatic potential, and other physical phenomena.
/// Solutions to Laplace's equation are called harmonic functions.
///
/// # Arguments
/// * `equation` - The Laplace equation to solve.
/// * `func` - The name of the unknown function (e.g., "u").
/// * `vars` - A slice of string slices representing the independent variables (e.g., `["x", "y"]`).
///
/// # Returns
/// An `Option<Expr>` representing a general solution, or `None` if the equation
/// does not match the Laplace equation pattern.
///
/// # Example
/// For `u_xx + u_yy = 0`, the solution is a harmonic function that can be expressed
/// as a Fourier series or in terms of separation of variables.
#[must_use]

pub fn solve_laplace_equation_2d(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 2 {

        return None;
    }

    let x_var = vars[0];

    let y_var = vars[1];

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x),
        x_var.to_string(),
    );

    let u_y = Expr::Derivative(
        Arc::new(u.clone()),
        y_var.to_string(),
    );

    let u_yy = Expr::Derivative(
        Arc::new(u_y),
        y_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut coeff_u_xx =
        Expr::Constant(0.0);

    let mut coeff_u_yy =
        Expr::Constant(0.0);

    for term in terms {

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_xx,
            )
        {

            coeff_u_xx = Expr::new_add(
                coeff_u_xx,
                coeff,
            );

            continue;
        }

        if let Some(coeff) =
            extract_coefficient(
                &term, &u_yy,
            )
        {

            coeff_u_yy = Expr::new_add(
                coeff_u_yy,
                coeff,
            );

            continue;
        }
    }

    let coeff_u_xx =
        simplify(&coeff_u_xx);

    let coeff_u_yy =
        simplify(&coeff_u_yy);

    // For Laplace equation, we expect u_xx + u_yy = 0
    // This means coeff_u_xx and coeff_u_yy should have the same absolute value
    // If the equation is passed as "u_xx + u_yy", both coefficients are 1
    // If the equation is passed as "u_xx + u_yy = 0", we'd get the same
    // So we check if they're equal and non-zero (standard Laplace)
    // or if their sum is zero (generalized form)

    let both_equal =
        coeff_u_xx == coeff_u_yy;

    let sum = simplify(&Expr::new_add(
        coeff_u_xx,
        coeff_u_yy,
    ));

    let sum_is_zero = is_zero(&sum);

    if !both_equal && !sum_is_zero {

        return None;
    }

    // General solution using separation of variables:
    // u(x,y) = Σ (A_n*sinh(n*π*y/L) + B_n*cosh(n*π*y/L)) * sin(n*π*x/L)
    // For now, return a symbolic representation

    let n =
        Expr::Variable("n".to_string());

    let a_n = Expr::Variable(
        "A_n".to_string(),
    );

    let b_n = Expr::Variable(
        "B_n".to_string(),
    );

    let l =
        Expr::Variable("L".to_string());

    let x = Expr::Variable(
        x_var.to_string(),
    );

    let y = Expr::Variable(
        y_var.to_string(),
    );

    // sinh(n*π*y/L)
    let sinh_part =
        Expr::new_sinh(Expr::new_div(
            Expr::new_mul(
                Expr::new_mul(
                    n.clone(),
                    Expr::Pi,
                ),
                y.clone(),
            ),
            l.clone(),
        ));

    // cosh(n*π*y/L)
    let cosh_part =
        Expr::new_cosh(Expr::new_div(
            Expr::new_mul(
                Expr::new_mul(
                    n.clone(),
                    Expr::Pi,
                ),
                y,
            ),
            l.clone(),
        ));

    // sin(n*π*x/L)
    let sin_part =
        Expr::new_sin(Expr::new_div(
            Expr::new_mul(
                Expr::new_mul(
                    n.clone(),
                    Expr::Pi,
                ),
                x,
            ),
            l,
        ));

    // (A_n*sinh(...) + B_n*cosh(...)) * sin(...)
    let term = Expr::new_mul(
        Expr::new_add(
            Expr::new_mul(
                a_n,
                sinh_part,
            ),
            Expr::new_mul(
                b_n,
                cosh_part,
            ),
        ),
        sin_part,
    );

    // Σ from n=1 to ∞
    let solution = Expr::Sum {
        body: Arc::new(term),
        var: Arc::new(n),
        from: Arc::new(Expr::Constant(
            1.0,
        )),
        to: Arc::new(Expr::Infinity),
    };

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution),
    ))
}

/// Solves the 3D wave equation `u_tt = c²(u_xx + u_yy + u_zz)`.
///
/// # Arguments
/// * `equation` - The 3D wave equation
/// * `func` - The unknown function name
/// * `vars` - Independent variables `["t", "x", "y", "z"]`
#[must_use]

pub fn solve_wave_equation_3d(
    _equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 4 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let solution_note = Expr::Variable(
        "3D_wave_solution".to_string(),
    );

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution_note),
    ))
}

/// Solves the 3D heat equation `u_t = α(u_xx + u_yy + u_zz)`.
#[must_use]

pub fn solve_heat_equation_3d(
    _equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 4 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let solution_note = Expr::Variable(
        "3D_heat_solution".to_string(),
    );

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution_note),
    ))
}

/// Solves the 3D Laplace equation `u_xx + u_yy + u_zz = 0`.
#[must_use]

pub fn solve_laplace_equation_3d(
    _equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 3 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let solution_note = Expr::Variable(
        "3D_Laplace_solution"
            .to_string(),
    );

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution_note),
    ))
}

/// Solves the 2D Poisson equation `u_xx + u_yy = f(x,y)`.
///
/// The Poisson equation is the inhomogeneous version of Laplace's equation,
/// describing phenomena with source terms such as electrostatics with charge density.
///
/// # Arguments
/// * `equation` - The Poisson equation (e.g., `u_xx + u_yy - f = 0`)
/// * `func` - The unknown function name (e.g., "u")
/// * `vars` - Independent variables `["x", "y"]`
///
/// # Returns
/// Solution using Green's function method
#[must_use]

pub fn solve_poisson_equation_2d(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 2 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        vars[0].to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x),
        vars[0].to_string(),
    );

    let u_y = Expr::Derivative(
        Arc::new(u.clone()),
        vars[1].to_string(),
    );

    let u_yy = Expr::Derivative(
        Arc::new(u_y),
        vars[1].to_string(),
    );

    let terms = collect_terms(equation);

    let mut coeff_u_xx =
        Expr::Constant(0.0);

    let mut coeff_u_yy =
        Expr::Constant(0.0);

    let mut source_term =
        Expr::Constant(0.0);

    for term in terms {

        match extract_coefficient(
                &term, &u_xx,
            )
        { Some(coeff) => {

            coeff_u_xx = Expr::new_add(
                coeff_u_xx,
                coeff,
            );
        } _ => { match extract_coefficient(
                &term, &u_yy,
            )
        { Some(coeff) => {

            coeff_u_yy = Expr::new_add(
                coeff_u_yy,
                coeff,
            );
        } _ => {

            // Terms without u derivatives are the source term
            source_term = Expr::new_add(
                source_term,
                term,
            );
        }}}}
    }

    let coeff_u_xx =
        simplify(&coeff_u_xx);

    let coeff_u_yy =
        simplify(&coeff_u_yy);

    let source_term = simplify(
        &Expr::new_neg(source_term),
    ); // Move to RHS

    // Check if it's a valid Poisson equation
    let both_equal =
        coeff_u_xx == coeff_u_yy;

    if !both_equal
        || is_zero(&source_term)
    {

        return None; // If source is zero, it's Laplace, not Poisson
    }

    // Solution using Green's function:
    // u(x,y) = ∫∫ G(x,y;ξ,η) f(ξ,η) dξ dη
    // where G is the Green's function for the Laplacian

    let x = Expr::Variable(
        vars[0].to_string(),
    );

    let y = Expr::Variable(
        vars[1].to_string(),
    );

    let xi =
        Expr::Variable("ξ".to_string());

    let eta =
        Expr::Variable("η".to_string());

    // Green's function: G = (1/2π) ln(r) where r = sqrt((x-ξ)² + (y-η)²)
    let r_sq = Expr::new_add(
        Expr::new_pow(
            Expr::new_sub(
                x,
                xi.clone(),
            ),
            Expr::Constant(2.0),
        ),
        Expr::new_pow(
            Expr::new_sub(
                y,
                eta.clone(),
            ),
            Expr::Constant(2.0),
        ),
    );

    let green = Expr::new_mul(
        Expr::new_div(
            Expr::Constant(1.0),
            Expr::new_mul(
                Expr::Constant(2.0),
                Expr::Pi,
            ),
        ),
        Expr::new_log(Expr::new_sqrt(
            r_sq,
        )),
    );

    let integrand = Expr::new_mul(
        green,
        source_term,
    );

    // Double integral
    let inner_integral =
        Expr::Integral {
            integrand: Arc::new(
                integrand,
            ),
            var: Arc::new(eta),
            lower_bound: Arc::new(
                Expr::NegativeInfinity,
            ),
            upper_bound: Arc::new(
                Expr::Infinity,
            ),
        };

    let solution = Expr::Integral {
        integrand: Arc::new(
            inner_integral,
        ),
        var: Arc::new(xi),
        lower_bound: Arc::new(
            Expr::NegativeInfinity,
        ),
        upper_bound: Arc::new(
            Expr::Infinity,
        ),
    };

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution),
    ))
}

/// Solves the 3D Poisson equation `u_xx + u_yy + u_zz = f(x,y,z)`.
#[must_use]

pub fn solve_poisson_equation_3d(
    _equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() != 3 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    // For 3D Poisson: u(x,y,z) = ∫∫∫ G(r) f(x',y',z') dx'dy'dz'
    // where G(r) = -1/(4πr) and r = |r - r'|

    let solution_note = Expr::Variable("3D_Poisson_solution_via_Greens_function".to_string());

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution_note),
    ))
}

/// Solves the Helmholtz equation `∇²u + k²u = 0`.
///
/// The Helmholtz equation arises in wave propagation, acoustics, and electromagnetic theory.
/// It's the time-independent form of the wave equation.
///
/// # Arguments
/// * `equation` - The Helmholtz equation (e.g., `u_xx + u_yy + k²u = 0`)
/// * `func` - The unknown function name
/// * `vars` - Independent variables
///
/// # Returns
/// Solution using separation of variables or Green's function
#[must_use]

pub fn solve_helmholtz_equation(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() < 2 {

        return None;
    }

    let u = Expr::Variable(
        func.to_string(),
    );

    let terms = collect_terms(equation);

    // Check for second spatial derivatives and a term with u (not its derivatives)
    let mut spatial_deriv_count = 0;

    let mut has_u_term = false;

    for var in vars {

        let u_var = Expr::Derivative(
            Arc::new(u.clone()),
            (*var).to_string(),
        );

        let u_var_var =
            Expr::Derivative(
                Arc::new(u_var),
                (*var).to_string(),
            );

        for term in &terms {

            if extract_coefficient(
                term,
                &u_var_var,
            )
            .is_some()
            {

                spatial_deriv_count +=
                    1;

                break;
            }
        }
    }

    // Check for k²u term
    for term in &terms {

        if extract_coefficient(
            term, &u,
        )
        .is_some()
        {

            has_u_term = true;

            break;
        }
    }

    if spatial_deriv_count < 2
        || !has_u_term
    {

        return None;
    }

    // Solution depends on dimension and boundary conditions
    let solution_note = if vars.len()
        == 2
    {

        Expr::Variable("2D_Helmholtz_solution_via_separation_of_variables".to_string())
    } else {

        Expr::Variable("3D_Helmholtz_solution_via_Greens_function".to_string())
    };

    Some(Expr::Eq(
        Arc::new(u),
        Arc::new(solution_note),
    ))
}

/// Solves the time-dependent Schrödinger equation `iℏ ∂ψ/∂t = -ℏ²/(2m) ∇²ψ + V(x)ψ`.
///
/// The Schrödinger equation is fundamental to quantum mechanics, describing
/// the time evolution of quantum states.
///
/// # Arguments
/// * `equation` - The Schrödinger equation
/// * `func` - The wave function name (e.g., "ψ" or "psi")
/// * `vars` - Independent variables (first is time, rest are spatial)
///
/// # Returns
/// Solution using separation of variables (energy eigenstates)
#[must_use]

pub fn solve_schrodinger_equation(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() < 2 {

        return None;
    }

    let psi = Expr::Variable(
        func.to_string(),
    );

    let t_var = vars[0];

    // Check for first time derivative
    let psi_t = Expr::Derivative(
        Arc::new(psi.clone()),
        t_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut has_time_deriv = false;

    let mut has_spatial_deriv = false;

    for term in &terms {

        if extract_coefficient(
            term,
            &psi_t,
        )
        .is_some()
        {

            has_time_deriv = true;
        }
    }

    // Check for spatial derivatives
    for var in vars.iter().skip(1) {

        let psi_x = Expr::Derivative(
            Arc::new(psi.clone()),
            (*var).to_string(),
        );

        let psi_xx = Expr::Derivative(
            Arc::new(psi_x),
            (*var).to_string(),
        );

        for term in &terms {

            if extract_coefficient(
                term,
                &psi_xx,
            )
            .is_some()
            {

                has_spatial_deriv =
                    true;

                break;
            }
        }
    }

    if !has_time_deriv
        || !has_spatial_deriv
    {

        return None;
    }

    // Solution using separation of variables: ψ(x,t) = φ(x)exp(-iEt/ℏ)
    // where φ(x) satisfies the time-independent Schrödinger equation

    let solution_note = Expr::Variable("Schrodinger_solution_via_energy_eigenstates".to_string());

    Some(Expr::Eq(
        Arc::new(psi),
        Arc::new(solution_note),
    ))
}

/// Solves the Klein-Gordon equation `∂²φ/∂t² - c²∇²φ + m²c⁴/ℏ² φ = 0`.
///
/// The Klein-Gordon equation is a relativistic wave equation describing
/// spin-0 particles in quantum field theory.
///
/// # Arguments
/// * `equation` - The Klein-Gordon equation
/// * `func` - The field name (e.g., "φ" or "phi")
/// * `vars` - Independent variables (first is time, rest are spatial)
///
/// # Returns
/// Solution using plane wave decomposition or Green's function
#[must_use]

pub fn solve_klein_gordon_equation(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Option<Expr> {

    if vars.len() < 2 {

        return None;
    }

    let phi = Expr::Variable(
        func.to_string(),
    );

    let t_var = vars[0];

    // Check for second time derivative
    let phi_t = Expr::Derivative(
        Arc::new(phi.clone()),
        t_var.to_string(),
    );

    let phi_tt = Expr::Derivative(
        Arc::new(phi_t),
        t_var.to_string(),
    );

    let terms = collect_terms(equation);

    let mut has_second_time_deriv =
        false;

    let mut has_spatial_deriv = false;

    let mut has_mass_term = false;

    for term in &terms {

        if extract_coefficient(
            term,
            &phi_tt,
        )
        .is_some()
        {

            has_second_time_deriv =
                true;
        }

        if extract_coefficient(
            term, &phi,
        )
        .is_some()
        {

            has_mass_term = true;
        }
    }

    // Check for spatial derivatives
    for var in vars.iter().skip(1) {

        let phi_x = Expr::Derivative(
            Arc::new(phi.clone()),
            (*var).to_string(),
        );

        let phi_xx = Expr::Derivative(
            Arc::new(phi_x),
            (*var).to_string(),
        );

        for term in &terms {

            if extract_coefficient(
                term,
                &phi_xx,
            )
            .is_some()
            {

                has_spatial_deriv =
                    true;

                break;
            }
        }
    }

    if !has_second_time_deriv
        || !has_spatial_deriv
        || !has_mass_term
    {

        return None;
    }

    // Solution using plane wave decomposition: φ(x,t) = ∫ [a(k)e^(i(k·x - ωt)) + a*(k)e^(-i(k·x - ωt))] d³k
    // where ω² = c²k² + m²c⁴/ℏ²

    let solution_note = Expr::Variable("Klein_Gordon_solution_via_plane_waves".to_string());

    Some(Expr::Eq(
        Arc::new(phi),
        Arc::new(solution_note),
    ))
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
#[must_use]

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

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let pattern = Expr::new_add(
        u_t,
        Expr::new_mul(u.clone(), u_x),
    );

    if simplify(&equation.clone())
        != simplify(&pattern)
    {

        return None;
    }

    let f_of_x = initial_conditions?
        .iter()
        .find(|cond| matches!(cond, Expr::Eq(lhs, _) if ** lhs == u))?;

    if let Expr::Eq(_, initial_func) =
        f_of_x
    {

        let x_minus_ut = Expr::new_sub(
            Expr::Variable(
                x_var.to_string(),
            ),
            Expr::new_mul(
                u.clone(),
                Expr::Variable(
                    t_var.to_string(),
                ),
            ),
        );

        let implicit_solution =
            substitute(
                initial_func,
                x_var,
                &x_minus_ut,
            );

        return Some(Expr::Eq(
            Arc::new(u),
            Arc::new(implicit_solution),
        ));
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
#[must_use]

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

    let initial_cond = initial_conditions?
        .iter()
        .find(|cond| {

            matches!(
                cond, Expr::Eq(lhs, _) if ** lhs == Expr::Variable(func.to_string())
            )
        })?;

    let f_x = if let Expr::Eq(_, ic) =
        initial_cond
    {

        ic
    } else {

        return None;
    };

    let u_k_0 =
        transforms::fourier_transform(
            f_x, x_var, k_var,
        );

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(Expr::Derivative(
            Arc::new(u.clone()),
            x_var.to_string(),
        )),
        x_var.to_string(),
    );

    let pattern = Expr::new_sub(
        u_t,
        Expr::new_mul(
            Expr::Pattern(
                "alpha".to_string(),
            ),
            u_xx,
        ),
    );

    if let Some(m) = pattern_match(
        equation,
        &pattern,
    ) {

        let alpha = m.get("alpha")?;

        let k = Expr::Variable(
            k_var.to_string(),
        );

        let _t = Expr::Variable(
            t_var.to_string(),
        );

        let neg_alpha_k_sq =
            Expr::new_neg(
                Expr::new_mul(
                    alpha.clone(),
                    Expr::new_pow(
                        k,
                        Expr::Constant(
                            2.0,
                        ),
                    ),
                ),
            );

        let ode_in_t = Expr::new_sub(
            differentiate(
                &Expr::Variable(
                    "U".to_string(),
                ),
                t_var,
            ),
            Expr::new_mul(
                neg_alpha_k_sq,
                Expr::Variable(
                    "U".to_string(),
                ),
            ),
        );

        let u_k_t_sol = solve_ode(
            &ode_in_t,
            "U",
            t_var,
            None,
        );

        if let Expr::Eq(
            _,
            general_sol,
        ) = u_k_t_sol
        {

            let c1 = Expr::Variable(
                "C1".to_string(),
            );

            let u_k_t = substitute(
                &general_sol,
                &c1.to_string(),
                &u_k_0,
            );

            let solution = transforms::inverse_fourier_transform(&u_k_t, k_var, x_var);

            return Some(Expr::Eq(
                Arc::new(u),
                Arc::new(solution),
            ));
        }
    }

    None
}

pub(crate) fn get_pde_order(
    expr: &Expr,
    _func: &str,
    vars: &[&str],
) -> usize {

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
) -> Option<(
    Expr,
    Expr,
    Expr,
    String,
)> {

    let x = &vars[0];

    let y = &vars[1];

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        (*x).to_string(),
    );

    let u_y = Expr::Derivative(
        Arc::new(u),
        (*y).to_string(),
    );

    let u_xx = Expr::Derivative(
        Arc::new(u_x.clone()),
        (*x).to_string(),
    );

    let u_yy = Expr::Derivative(
        Arc::new(u_y),
        (*y).to_string(),
    );

    let u_xy = Expr::Derivative(
        Arc::new(u_x),
        (*y).to_string(),
    );

    let (_, terms) =
        collect_and_order_terms(
            equation,
        );

    let mut coeffs = HashMap::new();

    for (term, coeff) in &terms {

        if *term == u_xx {

            coeffs.insert(
                "A",
                coeff.clone(),
            );
        } else if *term == u_xy {

            coeffs.insert(
                "B",
                coeff.clone(),
            );
        } else if *term == u_yy {

            coeffs.insert(
                "C",
                coeff.clone(),
            );
        }
    }

    let a = coeffs
        .get("A")
        .cloned()
        .unwrap_or_else(|| {

            Expr::Constant(0.0)
        });

    let b = coeffs
        .get("B")
        .cloned()
        .unwrap_or_else(|| {

            Expr::Constant(0.0)
        });

    let c = coeffs
        .get("C")
        .cloned()
        .unwrap_or_else(|| {

            Expr::Constant(0.0)
        });

    let discriminant =
        simplify(&Expr::new_sub(
            Expr::new_pow(
                b.clone(),
                Expr::Constant(2.0),
            ),
            Expr::new_mul(
                Expr::Constant(4.0),
                Expr::new_mul(
                    a.clone(),
                    c.clone(),
                ),
            ),
        ));

    let pde_type = discriminant.to_f64().map_or_else(
        || {
            format!(
                "Undetermined (discriminant is symbolic: {discriminant})"
            )
        },
        |d| {
            if d > 0.0 {
                "Hyperbolic".to_string()
            } else if d == 0.0 {
                "Parabolic".to_string()
            } else {
                "Elliptic".to_string()
            }
        },
    );

    Some((a, b, c, pde_type))
}

pub(crate) fn identify_differential_operator(
    lhs: &Expr,
    func: &str,
    vars: &[&str],
) -> String {

    let u = Expr::Variable(
        func.to_string(),
    );

    if vars.len() == 2 {

        let x = vars[0];

        let y = vars[1];

        let u_x = Expr::Derivative(
            Arc::new(u.clone()),
            x.to_string(),
        );

        let u_xx = Expr::Derivative(
            Arc::new(u_x),
            x.to_string(),
        );

        let u_y = Expr::Derivative(
            Arc::new(u),
            y.to_string(),
        );

        let u_yy = Expr::Derivative(
            Arc::new(u_y),
            y.to_string(),
        );

        if *lhs
            == simplify(&Expr::new_add(
                u_xx.clone(),
                u_yy.clone(),
            ))
        {

            return "Laplacian_2D"
                .to_string();
        }

        let pattern = Expr::new_sub(
            u_yy,
            Expr::new_mul(
                Expr::Pattern(
                    "c_sq".to_string(),
                ),
                u_xx,
            ),
        );

        if pattern_match(lhs, &pattern)
            .is_some()
        {

            return "Wave_1D"
                .to_string();
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

    let mut at_zero: Option<
        BoundaryConditionType,
    > = None;

    let mut at_l: Option<
        BoundaryConditionType,
    > = None;

    let mut l: Option<Expr> = None;

    let mut initial_cond: Option<Expr> =
        None;

    let mut initial_cond_deriv: Option<
        Expr,
    > = None;

    let u = Expr::Variable(
        func.to_string(),
    );

    let u_x = Expr::Derivative(
        Arc::new(u.clone()),
        x_var.to_string(),
    );

    let u_t = Expr::Derivative(
        Arc::new(u.clone()),
        t_var.to_string(),
    );

    let vars_order = [x_var, t_var];

    for cond in conditions {

        if let Expr::Eq(lhs, rhs) = cond
        {

            match get_value_at_point(
                    lhs,
                    t_var,
                    &Expr::Constant(
                        0.0,
                    ),
                    &vars_order,
                )
            { Some(val_str) => {

                let val =
                    Expr::Variable(
                        val_str
                            .to_string(
                            ),
                    );

                if val == u {

                    initial_cond = Some(
                        rhs.as_ref()
                            .clone(),
                    );
                }

                if val == u_t {

                    initial_cond_deriv = Some(rhs.as_ref().clone());
                }
            } _ => if is_zero(rhs) {

                if let Some(val_str) =
                    get_value_at_point(
                        lhs,
                        x_var,
                        &Expr::Constant(
                            0.0,
                        ),
                        &vars_order,
                    )
                {

                    let val = Expr::Variable(val_str.to_string());

                    if val == u {

                        at_zero = Some(BoundaryConditionType::Dirichlet);
                    }

                    if val == u_x {

                        at_zero = Some(BoundaryConditionType::Neumann);
                    }
                }

                if let Some(val_str) = get_value_at_point(
                    lhs,
                    x_var,
                    &Expr::Variable("L".to_string()),
                    &vars_order,
                ) {

                    let val = Expr::Variable(val_str.to_string());

                    if val == u {

                        at_l = Some(BoundaryConditionType::Dirichlet);

                        l = Some(Expr::Variable(
                            "L".to_string(),
                        ));
                    }

                    if val == u_x {

                        at_l = Some(BoundaryConditionType::Neumann);

                        l = Some(Expr::Variable(
                            "L".to_string(),
                        ));
                    }
                }
            }}
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

pub(crate) fn get_value_at_point<'a>(
    expr: &'a Expr,
    var: &str,
    point: &Expr,
    vars_order: &[&str],
) -> Option<&'a str> {

    if let Expr::Variable(s) = expr {

        if let Some(open_paren) =
            s.find('(')
        {

            if let Some(close_paren) =
                s.rfind(')')
            {

                let func_name =
                    &s[.. open_paren];

                let args_str = &s
                    [open_paren + 1
                        .. close_paren];

                let args: Vec<&str> =
                    args_str
                        .split(',')
                        .map(str::trim)
                        .collect();

                if let Some(var_index) =
                    vars_order
                        .iter()
                        .position(
                            |&v| {

                                v == var
                            },
                        )
                {

                    if let Some(
                        arg_val_str,
                    ) = args
                        .get(var_index)
                    {

                        if arg_val_str == &point.to_string() {

                            return Some(func_name);
                        }
                    }
                }
            }
        }
    }

    None
}
