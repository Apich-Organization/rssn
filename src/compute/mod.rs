//! # RSSN Unified Compute Engine
//!
//! Provides the primary computation entry point [`compute`] and configuration
//! [`ComputeConfig`] for transforming mathematical expressions into symbolic,
//! numerical, or exact result representations.

pub mod config;
pub mod operators;

pub use config::{ComputeConfig, TargetRepresentation};
pub use operators::*;

use ordered_float::OrderedFloat;

use crate::symbolic::core::Expr;
use crate::symbolic::egraph::cost::Extractor;
use crate::symbolic::egraph::egraph::EGraph;
use crate::symbolic::egraph::enode::ENode;
use crate::symbolic::egraph::rules::algebra::{ConstantFoldingRule, IdentityReductionRule};
use crate::symbolic::egraph::rules::calculus::DifferentiationRule;
use crate::symbolic::egraph::rules::oracles::{
    CombinatoricsOracleRule, ElementaryOracleRule, EulerLagrangeOracleRule, FactorOracleRule,
    FunctionalAnalysisOracleRule, GrobnerOracleRule, IndefiniteSumOracleRule,
    IntegralEquationsOracleRule, IntegralOracleRule, LieAlgebraOracleRule, LimitOracleRule,
    LogicOracleRule, MatrixOracleRule, NumberTheoryOracleRule, OdeOracleRule,
    OptimizationOracleRule, PdeOracleRule, PolynomialOracleRule, ProductOracleRule,
    QuantumOracleRule, RadicalsOracleRule, ResidueOracleRule, SeriesOracleRule, SolveOracleRule,
    SpecialFunctionsOracleRule, StatsOracleRule, SumOracleRule, TransformOracleRule,
    VectorCalculusOracleRule,
};
use crate::symbolic::egraph::rules::trig::TrigIdentitiesRule;
use crate::symbolic::egraph::rules::Rule;

/// Unified mathematical computation and equivalence reduction function.
///
/// Treats mathematical operators as identity transformation operators and reduces
/// them to the desired target representation kind (symbolic closed-form, numerical scalar/tensor,
/// or exact discrete form) according to the provided [`ComputeConfig`].
///
/// # Examples
///
/// ```rust
/// use rssn::compute::{compute, ComputeConfig};
/// use rssn::symbolic::core::Expr;
/// use num_bigint::BigInt;
///
/// let x = Expr::new_variable("x");
/// let x_sq = Expr::new_pow(x.clone(), Expr::new_bigint(2.into()));
/// let d_expr = Expr::new_derivative(x_sq, "x");
///
/// // 1. Symbolic reduction: d(x^2)/dx -> 2*x
/// let config_sym = ComputeConfig::default()
///     .with_calculus()
///     .target_symbolic();
/// let res_sym = compute(&d_expr, &config_sym);
/// assert!(matches!(res_sym, Expr::Mul(_, _)));
///
/// // 2. Numerical reduction with bound variable: d(x^2)/dx @ x=3 -> 6.0
/// let config_num = ComputeConfig::default()
///     .with_calculus()
///     .bind("x", 3.0)
///     .target_numerical(1e-7);
/// let res_num = compute(&d_expr, &config_num);
/// match res_num {
///     Expr::Constant(val) => assert!((val - 6.0).abs() < 1e-6),
///     Expr::BigInt(val) => assert_eq!(val, BigInt::from(6)),
///     other => panic!("Expected numerical constant, got: {:?}", other),
/// }
/// ```
#[must_use]
pub fn compute(expr: &Expr, config: &ComputeConfig) -> Expr {
    let mut egraph = EGraph::new();
    let root = egraph.add_expr(expr);

    egraph.rebuild();

    // Assemble active transformation rules based on ComputeConfig
    let mut rules: Vec<Box<dyn Rule>> = Vec::new();

    if config.enable_algebra {
        rules.push(Box::new(ConstantFoldingRule));
        rules.push(Box::new(IdentityReductionRule));
        rules.push(Box::new(TrigIdentitiesRule));
        rules.push(Box::new(ElementaryOracleRule));
        rules.push(Box::new(PolynomialOracleRule));
        rules.push(Box::new(GrobnerOracleRule));
        rules.push(Box::new(FactorOracleRule));
        rules.push(Box::new(RadicalsOracleRule));
    }

    if config.enable_calculus {
        rules.push(Box::new(DifferentiationRule));
        rules.push(Box::new(IntegralOracleRule::from_config(config)));
        rules.push(Box::new(LimitOracleRule));
        rules.push(Box::new(SeriesOracleRule));
        rules.push(Box::new(SumOracleRule));
        rules.push(Box::new(ProductOracleRule));
        rules.push(Box::new(IndefiniteSumOracleRule));
        rules.push(Box::new(VectorCalculusOracleRule));
        rules.push(Box::new(EulerLagrangeOracleRule));
        rules.push(Box::new(ResidueOracleRule));
    }

    if config.enable_ode {
        rules.push(Box::new(OdeOracleRule::with_config(config.clone())));
    }

    if config.enable_pde {
        rules.push(Box::new(PdeOracleRule));
    }

    if config.enable_matrix {
        rules.push(Box::new(MatrixOracleRule::from_config(config)));
    }

    if config.enable_transforms {
        rules.push(Box::new(TransformOracleRule));
    }

    if config.enable_special_functions {
        rules.push(Box::new(SpecialFunctionsOracleRule));
    }

    if config.enable_optimization {
        rules.push(Box::new(OptimizationOracleRule));
    }

    if config.enable_stats {
        rules.push(Box::new(StatsOracleRule));
    }

    if config.enable_physics {
        rules.push(Box::new(QuantumOracleRule));
    }

    // Always include Solve if calculus or algebra is on
    if config.enable_algebra || config.enable_ode || config.enable_calculus {
        rules.push(Box::new(SolveOracleRule::with_config(config.clone())));
        rules.push(Box::new(CombinatoricsOracleRule));
        rules.push(Box::new(LogicOracleRule));
        rules.push(Box::new(NumberTheoryOracleRule));
        rules.push(Box::new(FunctionalAnalysisOracleRule));
        rules.push(Box::new(LieAlgebraOracleRule));
        rules.push(Box::new(IntegralEquationsOracleRule));
    }

    // Separate rules by priority tiers
    let tier0_rules: Vec<&dyn Rule> = rules.iter().filter(|r| r.tier() == 0).map(AsRef::as_ref).collect();
    let tier1_rules: Vec<&dyn Rule> = rules.iter().filter(|r| r.tier() == 1).map(AsRef::as_ref).collect();
    let other_rules: Vec<&dyn Rule> = rules.iter().filter(|r| r.tier() >= 2).map(AsRef::as_ref).collect();

    let budget = &config.budget;
    let mut bindings_applied = false;

    for _iter in 0..budget.max_iterations {
        let mut iteration_changes = 0;

        // Phase 1: Fully resolve all high-weight operators (Derivatives, Integrals, ODEs, Solves)
        if budget.prioritize_high_weight_ops {
            for _ in 0..20 {
                let mut tier0_changes = 0;
                for rule in &tier0_rules {
                    tier0_changes += rule.apply(&mut egraph);
                }
                egraph.rebuild();
                iteration_changes += tier0_changes;
                if tier0_changes == 0 {
                    break;
                }
            }
        }

        // Apply bindings after operator de-cocooning so d(f(x), x) evaluates before x is bound to a constant
        if !bindings_applied && !config.bindings.is_empty() {
            for (var_name, &val) in &config.bindings {
                let var_id = egraph.add_node(ENode::Variable(var_name.clone()));
                let val_id = egraph.add_node(ENode::Constant(OrderedFloat(val)));
                egraph.union(var_id, val_id);
            }
            egraph.rebuild();
            bindings_applied = true;
            iteration_changes += 1;
        }

        // Phase 2: Constant folding & Identity reduction
        for rule in &tier1_rules {
            iteration_changes += rule.apply(&mut egraph);
        }
        egraph.rebuild();

        // Phase 3: Structural transformations
        if egraph.total_nodes() < budget.max_nodes && egraph.total_classes() < budget.max_classes {
            for rule in &other_rules {
                iteration_changes += rule.apply(&mut egraph);
            }
            egraph.rebuild();
        }

        // Fixpoint reached
        if iteration_changes == 0 {
            break;
        }

        // Node / Class budget bounds
        if egraph.total_nodes() >= budget.max_nodes || egraph.total_classes() >= budget.max_classes {
            break;
        }
    }

    let extractor = Extractor::new(&egraph);
    let extracted = extractor.extract(&egraph, root);
    apply_target_representation(extracted, config)
}

/// Transforms the extracted equivalence class root expression into the desired
/// representation domain specified by [`ComputeConfig::target`].
fn apply_target_representation(expr: Expr, config: &ComputeConfig) -> Expr {
    match &config.target {
        TargetRepresentation::Symbolic => expr,
        TargetRepresentation::Exact => expr,
        TargetRepresentation::Numerical { tolerance: _, max_iterations: _ } => {
            // 1. Direct floating-point conversion if available
            if let Some(val) = expr.to_f64() {
                return Expr::Constant(val);
            }
            // 2. Expression evaluation with bound variables
            if let Ok(val) = crate::numerical::elementary::eval_expr(&expr, &config.bindings) {
                return Expr::Constant(val);
            }
            // 3. Fallback numerical quadrature for unevaluated Integrals
            if let Expr::Integral { integrand, var, lower_bound, upper_bound } = &expr {
                let var_str = match &**var {
                    Expr::Variable(name) => name.clone(),
                    _ => var.to_string(),
                };
                let a_opt = lower_bound.to_f64().or_else(|| crate::numerical::elementary::eval_expr(lower_bound, &config.bindings).ok());
                let b_opt = upper_bound.to_f64().or_else(|| crate::numerical::elementary::eval_expr(upper_bound, &config.bindings).ok());
                if let (Some(a), Some(b)) = (a_opt, b_opt) {
                    if a.is_finite() && b.is_finite() {
                        let quad_res = crate::numerical::integrate::adaptive_quadrature(
                            |x: f64| -> f64 {
                                let mut vars = config.bindings.clone();
                                vars.insert(var_str.clone(), x);
                                crate::numerical::elementary::eval_expr(integrand, &vars).unwrap_or(f64::NAN)
                            },
                            (a, b),
                            1e-7,
                        );
                        if !quad_res.is_nan() {
                            return Expr::Constant(quad_res);
                        }
                    }
                }
            }
            // 4. Fallback numerical derivative for unevaluated Derivatives
            if let Expr::Derivative(body, var) = &expr {
                if let Some(&x_val) = config.bindings.get(var) {
                    if let Ok(val) = crate::numerical::calculus::partial_derivative(body, var, x_val) {
                        return Expr::Constant(val);
                    }
                }
            }
            // 5. Evaluate vector components if applicable (e.g. gradient)
            if let Expr::Vector(comps) = expr {
                let new_comps: Vec<Expr> = comps
                    .into_iter()
                    .map(|c| apply_target_representation(c, config))
                    .collect();
                return Expr::Vector(new_comps);
            }
            expr
        }
        TargetRepresentation::Auto => {
            if !config.bindings.is_empty() {
                if let Ok(val) = crate::numerical::elementary::eval_expr(&expr, &config.bindings) {
                    return Expr::Constant(val);
                }
            }
            expr
        }
    }
}
