//! # Heuristic E-Graph Symbolic Computation Engine
//!
//! This module provides a state-of-the-art equality saturation and expression rewriting
//! engine built directly upon Directed Acyclic Graph (DAG) equivalence classes.
//!
//! ## Architecture Highlights
//! - **DAG-Native Expression Sharing**: Subexpressions share underlying nodes via `Arc<Expr>`,
//!   eliminating memory explosion during large-scale saturation.
//! - **Heuristic Priority Cocooning**: High-weight operators (`Derivative`, `Integral`, `Ode`,
//!   `Solve`) are prioritized for reduction before structural rewrites are triggered.
//! - **Algorithmic Oracle Identities**: Integrates advanced domain algorithms (ODE solving,
//!   Gröbner bases, polynomial factorization) as first-class rewrite rules.
//! - **Chained Typestate Machine**: Type-safe pipeline builder ensuring valid rule set
//!   composition and execution state transitions.

#![allow(missing_docs)]
#![allow(clippy::module_name_repetitions)]

use std::sync::{Arc, LazyLock};

pub mod cost;
pub mod eclass;
pub mod egraph;
pub mod enode;
pub mod heuristics;
pub mod id;
pub mod rules;
pub mod typestate;
pub mod union_find;

#[cfg(test)]
mod tests;

pub use cost::Extractor;
pub use eclass::EClass;
pub use egraph::EGraph;
pub use enode::ENode;
pub use heuristics::HeuristicBudget;
pub use id::Id;
pub use typestate::{EGraphPipeline, PipelineBuilder, Ready, Unconfigured};

use crate::symbolic::core::Expr;

/// Global standard pipeline instance with full rule sets and oracles enabled.
static STANDARD_PIPELINE: LazyLock<EGraphPipeline<Ready>> = LazyLock::new(EGraphPipeline::standard);

/// Simplifies a symbolic expression using the standard heuristic E-Graph pipeline.
#[must_use]
pub fn simplify(expr: &Expr) -> Expr {
    STANDARD_PIPELINE.simplify(expr)
}

/// Differentiates an expression with respect to `var` by formulating `Derivative(expr, var)`
/// and simplifying it through the heuristic E-Graph engine.
#[must_use]
pub fn diff(
    expr: &Expr,
    var: &str,
) -> Expr {
    let deriv = Expr::Derivative(Arc::new(expr.clone()), var.to_string());
    simplify(&deriv)
}

/// Solves an equation or system by formulating `Solve(expr, var)` and resolving it through the E-Graph.
#[must_use]
pub fn solve_equation(
    expr: &Expr,
    var: &str,
) -> Expr {
    let solve_expr = Expr::Solve(Arc::new(expr.clone()), var.to_string());
    simplify(&solve_expr)
}

/// Integrates an expression with respect to `var` by formulating `Expr::Integral`
/// and resolving it through the heuristic E-Graph engine.
#[must_use]
pub fn integrate(
    expr: &Expr,
    var: &str,
    lower_bound: Option<&Expr>,
    upper_bound: Option<&Expr>,
) -> Expr {
    let int_expr = match (lower_bound, upper_bound) {
        (Some(lb), Some(ub)) => Expr::Integral {
            integrand: Arc::new(expr.clone()),
            var: Arc::new(Expr::Variable(var.to_string())),
            lower_bound: Arc::new((*lb).clone()),
            upper_bound: Arc::new((*ub).clone()),
        },
        _ => Expr::Integral {
            integrand: Arc::new(expr.clone()),
            var: Arc::new(Expr::Variable(var.to_string())),
            lower_bound: Arc::new(Expr::Variable("a".to_string())),
            upper_bound: Arc::new(Expr::Variable("b".to_string())),
        },
    };
    simplify(&int_expr)
}

/// Computes the limit of an expression by formulating `Expr::Limit`
/// and resolving it through the heuristic E-Graph engine.
#[must_use]
pub fn limit(
    expr: &Expr,
    var: &str,
    to: &Expr,
) -> Expr {
    let lim_expr = Expr::Limit(
        Arc::new(expr.clone()),
        var.to_string(),
        Arc::new(to.clone()),
    );
    simplify(&lim_expr)
}

/// Computes series expansion by formulating `Expr::Series` through the E-Graph.
#[must_use]
pub fn series(
    expr: &Expr,
    var: &str,
    point: &Expr,
    order: usize,
) -> Expr {
    let series_expr = Expr::Series(
        Arc::new(expr.clone()),
        var.to_string(),
        Arc::new(point.clone()),
        Arc::new(Expr::Constant(order as f64)),
    );
    simplify(&series_expr)
}

/// Solves an ODE by formulating `Expr::Ode` through the E-Graph.
#[must_use]
pub fn solve_ode(
    equation: &Expr,
    func: &str,
    var: &str,
) -> Expr {
    let ode_expr = Expr::Ode {
        equation: Arc::new(equation.clone()),
        func: func.to_string(),
        var: var.to_string(),
    };
    simplify(&ode_expr)
}

/// Solves a PDE by formulating `Expr::Pde` through the E-Graph.
#[must_use]
pub fn solve_pde(
    equation: &Expr,
    func: &str,
    vars: &[&str],
) -> Expr {
    let pde_expr = Expr::Pde {
        equation: Arc::new(equation.clone()),
        func: func.to_string(),
        vars: vars.iter().map(|s| (*s).to_string()).collect(),
    };
    simplify(&pde_expr)
}

/// Factorizes an expression by formulating `UnaryList("factor", ...)` and resolving through the E-Graph.
#[must_use]
pub fn factor(expr: &Expr) -> Expr {
    let factor_expr = Expr::UnaryList("factor".to_string(), Arc::new(expr.clone()));
    simplify(&factor_expr)
}
