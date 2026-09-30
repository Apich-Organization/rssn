use std::sync::Arc;
use std::collections::HashMap;
use num_traits::{One, Zero, Signed, ToPrimitive};
use num_bigint::BigInt;
use super::Rule;
use crate::compute::config::{ComputeConfig, TargetRepresentation};
use crate::symbolic::cas_foundations::{factorize_internal, simplify_with_relations_internal};
use crate::symbolic::core::Expr;
use crate::symbolic::egraph::cost::Extractor;
use crate::symbolic::egraph::egraph::EGraph;
use crate::symbolic::egraph::enode::ENode;
use crate::symbolic::egraph::id::Id;
use crate::symbolic::egraph::simplify;
use crate::symbolic::grobner::MonomialOrder;
use crate::symbolic::ode::{
    parse_ode,
    solve_ode_internal, solve_ode_system_internal, solve_separable_ode_internal,
    solve_first_order_linear_ode_internal, solve_bernoulli_ode_internal,
    solve_riccati_ode_internal, solve_cauchy_euler_ode_internal,
    solve_by_reduction_of_order_internal, solve_exact_ode_internal,
    solve_ode_by_series_internal, solve_ode_by_fourier_internal,
};
use crate::symbolic::solve::solve_internal;
use crate::symbolic::pde::{
    solve_pde_internal,
    solve_pde_by_separation_of_variables_internal,
    solve_pde_by_characteristics_internal,
    solve_pde_by_greens_function_internal,
    solve_second_order_pde_internal,
    solve_wave_equation_1d_dalembert_internal,
    solve_heat_equation_1d_internal,
    solve_laplace_equation_2d_internal,
    solve_wave_equation_3d_internal,
    solve_heat_equation_3d_internal,
    solve_laplace_equation_3d_internal,
    solve_poisson_equation_2d_internal,
    solve_poisson_equation_3d_internal,
    solve_helmholtz_equation_internal,
    solve_schrodinger_equation_internal,
    solve_klein_gordon_equation_internal,
    solve_burgers_equation_internal,
    solve_with_fourier_transform_internal,
};

fn extract_ode_rhs(eq: &Expr, func: &str, var: &str) -> Option<Expr> {
    let y_prime = Expr::Derivative(Arc::new(Expr::Variable(func.to_string())), var.to_string());
    match eq {
        Expr::Eq(l, r) => {
            if **l == y_prime {
                return Some(r.as_ref().clone());
            }
            if **r == y_prime {
                return Some(l.as_ref().clone());
            }
            if let Expr::Mul(c, inner) = l.as_ref() {
                if **inner == y_prime {
                    return Some(simplify(&Expr::new_div(r.as_ref().clone(), c.as_ref().clone())));
                }
            }
            let sub = simplify(&Expr::new_sub(l.as_ref().clone(), r.as_ref().clone()));
            extract_ode_rhs(&sub, func, var)
        }
        Expr::Sub(l, r) => {
            if **l == y_prime {
                return Some(r.as_ref().clone());
            }
            let parsed = parse_ode(eq, func, var);
            if parsed.order == 1 {
                let c1 = parsed.coeffs.get(&1).cloned().unwrap_or(Expr::Constant(1.0));
                let c0 = parsed.coeffs.get(&0).cloned().unwrap_or(Expr::Constant(0.0));
                let rem = parsed.remaining_expr;
                let numerator = simplify(&Expr::new_neg(Expr::new_add(Expr::new_mul(c0, Expr::Variable(func.to_string())), rem)));
                return Some(simplify(&Expr::new_div(numerator, c1)));
            }
            None
        }
        _ => {
            let parsed = parse_ode(eq, func, var);
            if parsed.order == 1 {
                let c1 = parsed.coeffs.get(&1).cloned().unwrap_or(Expr::Constant(1.0));
                let c0 = parsed.coeffs.get(&0).cloned().unwrap_or(Expr::Constant(0.0));
                let rem = parsed.remaining_expr;
                let numerator = simplify(&Expr::new_neg(Expr::new_add(Expr::new_mul(c0, Expr::Variable(func.to_string())), rem)));
                return Some(simplify(&Expr::new_div(numerator, c1)));
            }
            None
        }
    }
}

/// Ordinary Differential Equation Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Ode` and `ENode::OracleCall` with ODE solvers.
#[derive(Debug, Clone, Default)]
pub struct OdeOracleRule {
    pub config: Option<ComputeConfig>,
}

impl OdeOracleRule {
    #[must_use]
    pub fn new() -> Self {
        Self { config: None }
    }

    #[must_use]
    pub fn with_config(config: ComputeConfig) -> Self {
        Self { config: Some(config) }
    }
}

impl Rule for OdeOracleRule {
    fn name(&self) -> &str {
        "oracles::ode_solver"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut ode_targets = Vec::new();
        let mut oracle_call_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Ode { equation, func, var } => {
                        ode_targets.push((class_id, *equation, func.clone(), var.clone()));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "solve_ode_system"
                            | "solve_separable_ode"
                            | "solve_first_order_linear_ode"
                            | "solve_bernoulli_ode"
                            | "solve_riccati_ode"
                            | "solve_cauchy_euler_ode"
                            | "solve_by_reduction_of_order"
                            | "solve_exact_ode"
                            | "solve_ode_by_series"
                            | "solve_ode_by_fourier" => {
                                oracle_call_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if ode_targets.is_empty() && oracle_call_targets.is_empty() {
            return 0;
        }

        let is_numerical = match &self.config {
            Some(cfg) => match cfg.target {
                TargetRepresentation::Numerical { .. } => true,
                _ => false,
            },
            None => false,
        };

        let extractor = Extractor::new(egraph);
        let mut solutions = Vec::new();

        for (class_id, eq_id, func, var) in ode_targets {
            let eq_expr = extractor.extract(egraph, eq_id);
            if is_numerical {
                if let Some(cfg) = &self.config {
                    if let Some(rhs) = extract_ode_rhs(&eq_expr, &func, &var) {
                        let y0 = cfg.ode_initial_conditions.get(&func)
                            .or_else(|| cfg.bindings.get(&format!("{func}0")))
                            .or_else(|| cfg.bindings.get(&func))
                            .copied()
                            .unwrap_or(1.0);
                        let x0 = cfg.ode_range.map(|r| r.0)
                            .or_else(|| cfg.bindings.get(&format!("{var}0")).copied())
                            .unwrap_or(0.0);
                        let x_end = cfg.ode_range.map(|r| r.1)
                            .or_else(|| cfg.bindings.get(&var).copied())
                            .unwrap_or(1.0);
                        let num_steps = cfg.ode_steps.unwrap_or(100);

                        if let Ok(trajectory) = crate::numerical::ode::solve_ode_system_rk4_named(
                            &[rhs],
                            &var,
                            &[&func],
                            &[y0],
                            (x0, x_end),
                            num_steps,
                        ) {
                            if let Some(last_pt) = trajectory.last() {
                                if !last_pt.is_empty() {
                                    solutions.push((class_id, Expr::Constant(last_pt[0])));
                                    continue;
                                }
                            }
                        }
                    }
                }
            }

            let ics: Option<Vec<(Expr, u32, Expr)>> = self.config.as_ref().and_then(|cfg| {
                if !cfg.ode_initial_conditions.is_empty() {
                    let x0 = cfg.ode_range.map(|r| r.0).unwrap_or(0.0);
                    let list = cfg.ode_initial_conditions.iter().map(|(_f, &val)| {
                        (Expr::Constant(x0), 0, Expr::Constant(val))
                    }).collect();
                    Some(list)
                } else {
                    None
                }
            });
            let sol_expr = solve_ode_internal(&eq_expr, &func, &var, ics.as_deref());
            solutions.push((class_id, sol_expr));
        }

        for (class_id, name, args) in oracle_call_targets {
            match name.as_str() {
                "solve_ode_system" if args.len() >= 3 => {
                    let eqs_expr = extractor.extract(egraph, args[0]);
                    let fns_expr = extractor.extract(egraph, args[1]);
                    let var_expr = extractor.extract(egraph, args[2]);
                    let var_str = match &var_expr {
                        Expr::Variable(v) => v.clone(),
                        _ => format!("{var_expr}"),
                    };
                    let eqs = match eqs_expr {
                        Expr::Tuple(list) => list,
                        other => vec![other],
                    };
                    let fns_strings: Vec<String> = match fns_expr {
                        Expr::Tuple(list) => list
                            .into_iter()
                            .map(|e| match e {
                                Expr::Variable(v) => v,
                                other => format!("{other}"),
                            })
                            .collect(),
                        Expr::Variable(v) => vec![v],
                        other => vec![format!("{other}")],
                    };
                    let fn_refs: Vec<&str> = fns_strings.iter().map(String::as_str).collect();
                    if let Some(sols) = solve_ode_system_internal(&eqs, &fn_refs, &var_str) {
                        solutions.push((class_id, Expr::Solutions(sols)));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_separable_ode" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_separable_ode_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_first_order_linear_ode" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_first_order_linear_ode_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_bernoulli_ode" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_bernoulli_ode_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_riccati_ode" if args.len() >= 4 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    let y1 = extractor.extract(egraph, args[3]);
                    if let Some(sol) = solve_riccati_ode_internal(&eq, &func, &var, &y1) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_cauchy_euler_ode" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_cauchy_euler_ode_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_by_reduction_of_order" if args.len() >= 4 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    let known = extractor.extract(egraph, args[3]);
                    if let Some(sol) = solve_by_reduction_of_order_internal(&eq, &func, &var, &known) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_exact_ode" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_exact_ode_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_ode_by_series" if args.len() >= 6 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    let x0 = extractor.extract(egraph, args[3]);
                    let order_expr = extractor.extract(egraph, args[4]);
                    let order = order_expr.to_f64().unwrap_or(5.0) as u32;
                    let ics_expr = extractor.extract(egraph, args[5]);
                    let mut ics = Vec::new();
                    if let Expr::Tuple(list) = ics_expr {
                        for item in list {
                            if let Expr::Tuple(pair) = item {
                                if pair.len() == 2 {
                                    let k = pair[0].to_f64().unwrap_or(0.0) as u32;
                                    ics.push((k, pair[1].clone()));
                                }
                            }
                        }
                    }
                    if let Some(sol) = solve_ode_by_series_internal(&eq, &func, &var, &x0, order, &ics) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                "solve_ode_by_fourier" if args.len() >= 3 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let func = format!("{}", extractor.extract(egraph, args[1]));
                    let var = format!("{}", extractor.extract(egraph, args[2]));
                    if let Some(sol) = solve_ode_by_fourier_internal(&eq, &func, &var) {
                        solutions.push((class_id, sol));
                    } else {
                        solutions.push((class_id, Expr::NoSolution));
                    }
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, sol_expr) in solutions {
            let sol_id = egraph.add_expr(&sol_expr);
            if egraph.union(class_id, sol_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Algebraic Equation Solver Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Solve` with `rssn::symbolic::solve::solve` and numerical root solvers.
#[derive(Debug, Clone, Default)]
pub struct SolveOracleRule {
    pub config: Option<ComputeConfig>,
}

impl SolveOracleRule {
    #[must_use]
    pub fn new() -> Self {
        Self { config: None }
    }

    #[must_use]
    pub fn with_config(config: ComputeConfig) -> Self {
        Self { config: Some(config) }
    }
}

impl Rule for SolveOracleRule {
    fn name(&self) -> &str {
        "oracles::equation_solver"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Solve(equation, var) = node {
                    targets.push((class_id, *equation, var.clone()));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let is_numerical = match &self.config {
            Some(cfg) => match cfg.target {
                TargetRepresentation::Numerical { .. } => true,
                TargetRepresentation::Auto => !cfg.bindings.is_empty(),
                _ => false,
            },
            None => false,
        };

        let extractor = Extractor::new(egraph);
        let mut solutions = Vec::new();

        for (class_id, eq_id, var) in targets {
            let eq_expr = extractor.extract(egraph, eq_id);
            if is_numerical {
                if let Some(cfg) = &self.config {
                    let normalized = if let Expr::Eq(l, r) = &eq_expr {
                        simplify(&Expr::new_sub(l.as_ref().clone(), r.as_ref().clone()))
                    } else {
                        eq_expr.clone()
                    };
                    let (tol, max_iter) = match cfg.target {
                        TargetRepresentation::Numerical { tolerance, max_iterations } => (tolerance, max_iterations),
                        _ => (1e-7, 1000),
                    };
                    let start_guess = cfg.bindings.get(&var).copied();
                    if let Ok(root_val) = crate::numerical::solve::solve_root(&normalized, &var, start_guess, tol, max_iter) {
                        solutions.push((class_id, vec![Expr::Constant(root_val)]));
                        continue;
                    }
                }
            }

            let sols = solve_internal(&eq_expr, &var);
            solutions.push((class_id, sols));
        }

        let mut applied = 0;
        for (class_id, sols) in solutions {
            if sols.is_empty() {
                let no_sol = egraph.add_node(ENode::NoSolution);
                if egraph.union(class_id, no_sol) {
                    applied += 1;
                }
            } else if sols.len() == 1 {
                let sol_id = egraph.add_expr(&sols[0]);
                if egraph.union(class_id, sol_id) {
                    applied += 1;
                }
            } else {
                let sol_ids: Vec<Id> = sols.iter().map(|s| egraph.add_expr(s)).collect();
                let sols_node = egraph.add_node(ENode::Solutions(sol_ids));
                if egraph.union(class_id, sols_node) {
                    applied += 1;
                }
            }
        }
        applied
    }
}

/// Gröbner Basis Reduction Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::SimplifyWithRelations` with `rssn::symbolic::cas_foundations::simplify_with_relations`.
#[derive(Debug)]
pub struct GrobnerOracleRule;

impl Rule for GrobnerOracleRule {
    fn name(&self) -> &str {
        "oracles::grobner_relations"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::SimplifyWithRelations { expr, relations, vars } => {
                        targets.push((class_id, *expr, relations.clone(), vars.clone()));
                    }
                    ENode::OracleCall(name, args) if name == "simplify_with_relations" && args.len() >= 2 => {
                        oracle_targets.push((class_id, args.clone()));
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, expr_id, rel_ids, vars) in targets {
            let expr = extractor.extract(egraph, expr_id);
            let relations: Vec<Expr> = rel_ids.iter().map(|&id| extractor.extract(egraph, id)).collect();
            let var_refs: Vec<&str> = vars.iter().map(String::as_str).collect();

            let simplified = simplify_with_relations_internal(
                &expr,
                &relations,
                &var_refs,
                MonomialOrder::Lexicographical,
            );
            results.push((class_id, simplified));
        }

        for (class_id, args) in oracle_targets {
            let expr = extractor.extract(egraph, args[0]);
            let rel_extracted = extractor.extract(egraph, args[1]);
            let relations: Vec<Expr> = match rel_extracted {
                Expr::Vector(vec) => vec,
                Expr::Tuple(vec) => vec,
                other => vec![other],
            };
            let mut vars = Vec::new();
            for &arg_id in &args[2..] {
                let v_expr = extractor.extract(egraph, arg_id);
                match v_expr {
                    Expr::Variable(name) => vars.push(name),
                    other => vars.push(format!("{}", other)),
                }
            }
            let var_refs: Vec<&str> = vars.iter().map(String::as_str).collect();

            let simplified = simplify_with_relations_internal(
                &expr,
                &relations,
                &var_refs,
                MonomialOrder::Lexicographical,
            );
            results.push((class_id, simplified));
        }

        let mut applied = 0;
        for (class_id, simplified_expr) in results {
            let res_id = egraph.add_expr(&simplified_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Symbolic & Numerical Integration Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Integral` with `rssn::symbolic::calculus::integrate_internal`
/// and numerical quadrature via `rssn::numerical::integrate::adaptive_quadrature`.
#[derive(Debug, Clone, Default)]
pub struct IntegralOracleRule {
    pub target: Option<crate::compute::config::TargetRepresentation>,
    pub bindings: HashMap<String, f64>,
}

impl IntegralOracleRule {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn from_config(config: &crate::compute::config::ComputeConfig) -> Self {
        Self {
            target: Some(config.target.clone()),
            bindings: config.bindings.clone(),
        }
    }

    pub fn is_numerical_target(&self) -> bool {
        match &self.target {
            Some(crate::compute::config::TargetRepresentation::Numerical { .. }) => true,
            _ => false,
        }
    }

    pub fn numerical_tolerance(&self) -> f64 {
        match &self.target {
            Some(crate::compute::config::TargetRepresentation::Numerical { tolerance, .. }) => *tolerance,
            _ => 1e-6,
        }
    }
}

impl Rule for IntegralOracleRule {
    fn name(&self) -> &str {
        "oracles::integral_solver"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Integral { integrand, var, lower_bound, upper_bound } = node {
                    targets.push((class_id, *integrand, *var, *lower_bound, *upper_bound));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, int_id, var_id, lb_id, ub_id) in targets {
            let integrand = extractor.extract(egraph, int_id);
            let var_expr = extractor.extract(egraph, var_id);
            let var_str = match &var_expr {
                Expr::Variable(name) => name.clone(),
                _ => var_expr.to_string(),
            };
            let lb = extractor.extract(egraph, lb_id);
            let ub = extractor.extract(egraph, ub_id);

            let is_indefinite = match (&lb, &ub) {
                (Expr::Variable(a), Expr::Variable(b)) if a == "a" && b == "b" => true,
                _ => false,
            };

            let sym_res = if is_indefinite {
                crate::symbolic::calculus::integrate_internal(&integrand, &var_str, None, None)
            } else {
                crate::symbolic::calculus::integrate_internal(&integrand, &var_str, Some(&lb), Some(&ub))
            };

            let sym_success = !matches!(sym_res, Expr::Integral { .. });

            // If definite integral with numeric bounds, evaluate numerical quadrature
            let lb_num = lb.to_f64().or_else(|| crate::numerical::elementary::eval_expr(&lb, &self.bindings).ok());
            let ub_num = ub.to_f64().or_else(|| crate::numerical::elementary::eval_expr(&ub, &self.bindings).ok());

            let mut quad_val: Option<f64> = None;
            if !is_indefinite {
                if let (Some(a), Some(b)) = (lb_num, ub_num) {
                    if a.is_finite() && b.is_finite() {
                        let tol = self.numerical_tolerance();
                        let q = crate::numerical::integrate::adaptive_quadrature(
                            |x: f64| -> f64 {
                                let mut vars: HashMap<String, f64> = self.bindings.clone();
                                vars.insert(var_str.clone(), x);
                                crate::numerical::elementary::eval_expr(&integrand, &vars).unwrap_or(f64::NAN)
                            },
                            (a, b),
                            tol,
                        );
                        if !q.is_nan() {
                            quad_val = Some(q);
                        }
                    }
                }
            }

            if self.is_numerical_target() {
                if let Some(val) = quad_val {
                    results.push((class_id, Expr::Constant(val)));
                } else if sym_success {
                    results.push((class_id, sym_res));
                }
            } else {
                if sym_success {
                    results.push((class_id, sym_res));
                    if let Some(val) = quad_val {
                        results.push((class_id, Expr::Constant(val)));
                    }
                } else if let Some(val) = quad_val {
                    results.push((class_id, Expr::Constant(val)));
                }
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Limit Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Limit` with `rssn::symbolic::calculus::limit_internal`.
#[derive(Debug)]
pub struct LimitOracleRule;

impl Rule for LimitOracleRule {
    fn name(&self) -> &str {
        "oracles::limit_solver"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Limit(body, var, point) = node {
                    targets.push((class_id, *body, var.clone(), *point));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, body_id, var, point_id) in targets {
            let body = extractor.extract(egraph, body_id);
            let point = extractor.extract(egraph, point_id);
            let res = crate::symbolic::calculus::limit_internal(&body, &var, &point, 0);
            if !matches!(res, Expr::Limit(..)) {
                results.push((class_id, res));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Partial Differential Equation Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Pde` with `rssn::symbolic::pde::solve_pde`.
#[derive(Debug)]
pub struct PdeOracleRule;

impl Rule for PdeOracleRule {
    fn name(&self) -> &str {
        "oracles::pde_solver"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_call_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Pde { equation, func, vars } => {
                        targets.push((class_id, *equation, func.clone(), vars.clone()));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "solve_pde_by_separation_of_variables"
                            | "solve_pde_by_characteristics"
                            | "solve_pde_by_greens_function"
                            | "solve_second_order_pde"
                            | "solve_wave_equation_1d_dalembert"
                            | "solve_heat_equation_1d"
                            | "solve_laplace_equation_2d"
                            | "solve_wave_equation_3d"
                            | "solve_heat_equation_3d"
                            | "solve_laplace_equation_3d"
                            | "solve_poisson_equation_2d"
                            | "solve_poisson_equation_3d"
                            | "solve_helmholtz_equation"
                            | "solve_schrodinger_equation"
                            | "solve_klein_gordon_equation"
                            | "solve_burgers_equation"
                            | "solve_with_fourier_transform" => {
                                oracle_call_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_call_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, eq_id, func, vars) in targets {
            let eq = extractor.extract(egraph, eq_id);
            let var_refs: Vec<&str> = vars.iter().map(String::as_str).collect();
            let res = solve_pde_internal(&eq, &func, &var_refs, None);
            if !matches!(res, Expr::Pde { .. }) {
                results.push((class_id, res));
            }
        }

        for (class_id, name, args) in oracle_call_targets {
            if args.len() < 3 {
                continue;
            }
            let eq = extractor.extract(egraph, args[0]);
            let func = match &extractor.extract(egraph, args[1]) {
                Expr::Variable(v) => v.clone(),
                other => format!("{other}"),
            };
            let vars_expr = extractor.extract(egraph, args[2]);
            let vars_strings: Vec<String> = match vars_expr {
                Expr::Tuple(list) => list
                    .into_iter()
                    .map(|e| match e {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    })
                    .collect(),
                Expr::Variable(v) => vec![v],
                other => vec![format!("{other}")],
            };
            let var_refs: Vec<&str> = vars_strings.iter().map(String::as_str).collect();

            let conds_vec: Option<Vec<Expr>> = if args.len() > 3 {
                match extractor.extract(egraph, args[3]) {
                    Expr::Tuple(list) => Some(list),
                    other => Some(vec![other]),
                }
            } else {
                None
            };
            let conds_slice = conds_vec.as_deref();

            let maybe_sol = match name.as_str() {
                "solve_pde_by_separation_of_variables" => {
                    solve_pde_by_separation_of_variables_internal(&eq, &func, &var_refs, conds_slice.unwrap_or(&[]))
                }
                "solve_pde_by_characteristics" => {
                    solve_pde_by_characteristics_internal(&eq, &func, &var_refs)
                }
                "solve_pde_by_greens_function" => {
                    solve_pde_by_greens_function_internal(&eq, &func, &var_refs)
                }
                "solve_second_order_pde" => {
                    solve_second_order_pde_internal(&eq, &func, &var_refs)
                }
                "solve_wave_equation_1d_dalembert" => {
                    solve_wave_equation_1d_dalembert_internal(&eq, &func, &var_refs)
                }
                "solve_heat_equation_1d" => {
                    solve_heat_equation_1d_internal(&eq, &func, &var_refs)
                }
                "solve_laplace_equation_2d" => {
                    solve_laplace_equation_2d_internal(&eq, &func, &var_refs)
                }
                "solve_wave_equation_3d" => {
                    solve_wave_equation_3d_internal(&eq, &func, &var_refs)
                }
                "solve_heat_equation_3d" => {
                    solve_heat_equation_3d_internal(&eq, &func, &var_refs)
                }
                "solve_laplace_equation_3d" => {
                    solve_laplace_equation_3d_internal(&eq, &func, &var_refs)
                }
                "solve_poisson_equation_2d" => {
                    solve_poisson_equation_2d_internal(&eq, &func, &var_refs)
                }
                "solve_poisson_equation_3d" => {
                    solve_poisson_equation_3d_internal(&eq, &func, &var_refs)
                }
                "solve_helmholtz_equation" => {
                    solve_helmholtz_equation_internal(&eq, &func, &var_refs)
                }
                "solve_schrodinger_equation" => {
                    solve_schrodinger_equation_internal(&eq, &func, &var_refs)
                }
                "solve_klein_gordon_equation" => {
                    solve_klein_gordon_equation_internal(&eq, &func, &var_refs)
                }
                "solve_burgers_equation" => {
                    solve_burgers_equation_internal(&eq, &func, &var_refs, conds_slice)
                }
                "solve_with_fourier_transform" => {
                    solve_with_fourier_transform_internal(&eq, &func, &var_refs, conds_slice)
                }
                _ => None,
            };

            if let Some(sol) = maybe_sol {
                results.push((class_id, sol));
            } else {
                results.push((class_id, Expr::NoSolution));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Polynomial Factorization Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Factor` with `rssn::symbolic::cas_foundations::factorize`.
#[derive(Debug)]
pub struct FactorOracleRule;

impl Rule for FactorOracleRule {
    fn name(&self) -> &str {
        "oracles::factorizer"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Factor(inner) = node {
                    targets.push((class_id, *inner));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, inner_id) in targets {
            let inner = extractor.extract(egraph, inner_id);
            let res = factorize_internal(inner);
            results.push((class_id, res));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Series Expansion Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Series` with `rssn::symbolic::series::taylor_series`.
#[derive(Debug)]
pub struct SeriesOracleRule;

impl Rule for SeriesOracleRule {
    fn name(&self) -> &str {
        "oracles::series_expansion"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Series(body, var, point, order) => {
                        targets.push((class_id, *body, var.clone(), *point, *order));
                    }
                    ENode::OracleCall(name, args) => {
                        if (name == "laurent_series" || name == "fourier_series") && args.len() >= 4 {
                            oracle_targets.push((class_id, name.clone(), args.clone()));
                        }
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, body_id, var, point_id, order_id) in targets {
            let body = extractor.extract(egraph, body_id);
            let point = extractor.extract(egraph, point_id);
            let order_expr = extractor.extract(egraph, order_id);
            let order = order_expr.to_f64().unwrap_or(5.0) as usize;
            let res = crate::symbolic::series::taylor_series_internal(&body, &var, &point, order);
            results.push((class_id, res));
        }

        for (class_id, name, args) in oracle_targets {
            let expr = extractor.extract(egraph, args[0]);
            let var_expr = extractor.extract(egraph, args[1]);
            let var = match var_expr {
                Expr::Variable(v) => v,
                other => format!("{}", other),
            };
            let p_or_c = extractor.extract(egraph, args[2]);
            let order_expr = extractor.extract(egraph, args[3]);
            let order = order_expr.to_f64().unwrap_or(5.0) as usize;

            if name == "laurent_series" {
                let res = crate::symbolic::series::laurent_series_internal(&expr, &var, &p_or_c, order);
                results.push((class_id, res));
            } else if name == "fourier_series" {
                let res = crate::symbolic::series::fourier_series_internal(&expr, &var, &p_or_c, order);
                results.push((class_id, res));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Summation Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Sum` with `rssn::symbolic::series::summation`.
#[derive(Debug)]
pub struct SumOracleRule;

impl Rule for SumOracleRule {
    fn name(&self) -> &str {
        "oracles::summation"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Sum { body, var, from, to } = node {
                    targets.push((class_id, *body, *var, *from, *to));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, body_id, var_id, from_id, to_id) in targets {
            let body = extractor.extract(egraph, body_id);
            let var_expr = extractor.extract(egraph, var_id);
            let var_str = match &var_expr {
                Expr::Variable(name) => name.clone(),
                _ => var_expr.to_string(),
            };
            let from = extractor.extract(egraph, from_id);
            let to = extractor.extract(egraph, to_id);
            let res = crate::symbolic::series::summation_internal(&body, &var_str, &from, &to);
            results.push((class_id, res));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Product Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Product` with `rssn::symbolic::series::product`.
#[derive(Debug)]
pub struct ProductOracleRule;

impl Rule for ProductOracleRule {
    fn name(&self) -> &str {
        "oracles::product"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Product(body, var, from, to) = node {
                    targets.push((class_id, *body, var.clone(), *from, *to));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, body_id, var, from_id, to_id) in targets {
            let body = extractor.extract(egraph, body_id);
            let from = extractor.extract(egraph, from_id);
            let to = extractor.extract(egraph, to_id);
            let res = crate::symbolic::series::product_internal(&body, &var, &from, &to);
            results.push((class_id, res));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Indefinite Summation Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::IndefiniteSum` with closed-form summation in `crate::numerical::indefinite_sum`.
#[derive(Debug)]
pub struct IndefiniteSumOracleRule;

impl Rule for IndefiniteSumOracleRule {
    fn name(&self) -> &str {
        "oracles::indefinite_sum"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::IndefiniteSum { body, var, step } = node {
                    targets.push((class_id, *body, var.clone(), *step));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, body_id, var, step_id) in targets {
            let body = extractor.extract(egraph, body_id);
            let step = extractor.extract(egraph, step_id);

            let is_step_one = match &step {
                Expr::Constant(c) => (*c - 1.0).abs() < 1e-9,
                Expr::BigInt(i) => i.is_one(),
                Expr::Rational(r) => r.is_one(),
                _ => false,
            };

            if is_step_one {
                if let Some(closed_form) = crate::numerical::indefinite_sum::try_closed_form_sum(&body, &var) {
                    results.push((class_id, closed_form));
                }
            } else {
                let u_var = "_u";
                let u_expr = Expr::new_variable(u_var);
                let replacement = Expr::new_mul(u_expr, step.clone());
                let substituted_body = crate::symbolic::calculus::substitute(&body, &var, &replacement);

                if let Some(closed_form_u) = crate::numerical::indefinite_sum::try_closed_form_sum(&substituted_body, u_var) {
                    let u_back = Expr::new_div(Expr::new_variable(&var), step);
                    let final_expr = crate::symbolic::calculus::substitute(&closed_form_u, u_var, &u_back);
                    results.push((class_id, final_expr));
                }
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Radicals Denesting Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges `ENode::Sqrt` with `rssn::symbolic::radicals::denest_sqrt`.
#[derive(Debug)]
pub struct RadicalsOracleRule;

impl Rule for RadicalsOracleRule {
    fn name(&self) -> &str {
        "oracles::radicals_denesting"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Sqrt(inner) = node {
                    targets.push((class_id, *inner));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, inner_id) in targets {
            let inner = extractor.extract(egraph, inner_id);
            let sqrt_expr = Expr::new_sqrt(inner);
            let denested = crate::symbolic::radicals::denest_sqrt(&sqrt_expr);
            if denested != sqrt_expr {
                results.push((class_id, denested));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Matrix Operations Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges E-Graph matrix operations (Inverse, Determinant, Trace, Transpose, MatrixMul)
/// with `rssn::symbolic::matrix`.
#[derive(Debug, Clone)]
pub struct MatrixOracleRule {
    pub target: Option<crate::compute::config::TargetRepresentation>,
    pub bindings: HashMap<String, f64>,
}

impl Default for MatrixOracleRule {
    fn default() -> Self {
        Self {
            target: None,
            bindings: HashMap::new(),
        }
    }
}

impl MatrixOracleRule {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn from_config(config: &crate::compute::config::ComputeConfig) -> Self {
        Self {
            target: Some(config.target.clone()),
            bindings: config.bindings.clone(),
        }
    }

    pub fn is_numerical_target(&self) -> bool {
        match &self.target {
            Some(crate::compute::config::TargetRepresentation::Numerical { .. }) => true,
            _ => false,
        }
    }

    pub fn numerical_tolerance(&self) -> f64 {
        match &self.target {
            Some(crate::compute::config::TargetRepresentation::Numerical { tolerance, .. }) => *tolerance,
            _ => 1e-10,
        }
    }
}

fn expr_to_f64_val(expr: &Expr, bindings: &HashMap<String, f64>) -> Option<f64> {
    match expr {
        Expr::Constant(c) => Some(*c),
        Expr::BigInt(b) => b.to_f64(),
        Expr::Rational(r) => r.to_f64(),
        Expr::Variable(v) => bindings.get(v).copied(),
        Expr::Neg(inner) => expr_to_f64_val(inner, bindings).map(|x| -x),
        Expr::Add(a, b) => Some(expr_to_f64_val(a, bindings)? + expr_to_f64_val(b, bindings)?),
        Expr::Sub(a, b) => Some(expr_to_f64_val(a, bindings)? - expr_to_f64_val(b, bindings)?),
        Expr::Mul(a, b) => Some(expr_to_f64_val(a, bindings)? * expr_to_f64_val(b, bindings)?),
        Expr::Div(a, b) => {
            let denom = expr_to_f64_val(b, bindings)?;
            if denom.abs() > 1e-15 {
                Some(expr_to_f64_val(a, bindings)? / denom)
            } else {
                None
            }
        }
        _ => None,
    }
}

fn extract_f64_matrix(
    expr: &Expr,
    bindings: &HashMap<String, f64>,
) -> Option<(usize, usize, Vec<f64>)> {
    if let Expr::Matrix(rows) = expr {
        if rows.is_empty() {
            return Some((0, 0, Vec::new()));
        }
        let num_rows = rows.len();
        let num_cols = rows[0].len();
        if !rows.iter().all(|r| r.len() == num_cols) {
            return None;
        }
        let mut data = Vec::with_capacity(num_rows * num_cols);
        for row in rows {
            for elem in row {
                data.push(expr_to_f64_val(elem, bindings)?);
            }
        }
        Some((num_rows, num_cols, data))
    } else {
        None
    }
}

fn numerical_matrix_inverse(rows: usize, cols: usize, data: &[f64]) -> Option<Expr> {
    if rows != cols || rows == 0 {
        return None;
    }
    // 1. Try Faer backend via crate::numerical::matrix::Matrix
    let faer_mat = crate::numerical::matrix::Matrix::new(rows, cols, data.to_vec())
        .with_backend(crate::numerical::matrix::Backend::Faer);
    if let Some(inv) = faer_mat.inverse() {
        let inv_data = inv.into_data();
        let mut res_rows = Vec::with_capacity(rows);
        for i in 0..rows {
            let mut row = Vec::with_capacity(cols);
            for j in 0..cols {
                row.push(Expr::Constant(inv_data[i * cols + j]));
            }
            res_rows.push(row);
        }
        return Some(Expr::Matrix(res_rows));
    }

    // 2. Try nalgebra as fallback
    let dmat = nalgebra::DMatrix::from_row_slice(rows, cols, data);
    if let Some(inv_dmat) = dmat.try_inverse() {
        let mut res_rows = Vec::with_capacity(rows);
        for i in 0..rows {
            let mut row = Vec::with_capacity(cols);
            for j in 0..cols {
                row.push(Expr::Constant(inv_dmat[(i, j)]));
            }
            res_rows.push(row);
        }
        return Some(Expr::Matrix(res_rows));
    }

    None
}

fn numerical_determinant(rows: usize, cols: usize, data: &[f64]) -> Option<f64> {
    if rows != cols || rows == 0 {
        return None;
    }
    // 1. Try Faer/Matrix LU determinant
    let mat = crate::numerical::matrix::Matrix::new(rows, cols, data.to_vec());
    if let Ok(det_val) = mat.determinant_lu() {
        return Some(det_val);
    }
    // 2. Try nalgebra
    let dmat = nalgebra::DMatrix::from_row_slice(rows, cols, data);
    Some(dmat.determinant())
}

fn numerical_eigenvalues(
    rows: usize,
    cols: usize,
    data: &[f64],
    tolerance: f64,
) -> Option<Expr> {
    if rows != cols || rows == 0 {
        return None;
    }

    // Check if matrix is symmetric
    let mut is_symmetric = true;
    for i in 0..rows {
        for j in (i + 1)..cols {
            if (data[i * cols + j] - data[j * cols + i]).abs() > 1e-12 {
                is_symmetric = false;
                break;
            }
        }
        if !is_symmetric {
            break;
        }
    }

    if is_symmetric {
        // Option A: Try Faer self-adjoint eigen decomposition
        let mat = crate::numerical::matrix::Matrix::new(rows, cols, data.to_vec())
            .with_backend(crate::numerical::matrix::Backend::Faer);
        if let Some(crate::numerical::matrix::FaerDecompositionResult::EigenSymmetric { values, .. }) =
            mat.decompose(crate::numerical::matrix::FaerDecompositionType::EigenSymmetric)
        {
            let eig_rows: Vec<Vec<Expr>> = values.into_iter().map(|v| vec![Expr::Constant(v)]).collect();
            return Some(Expr::Matrix(eig_rows));
        }

        // Option B: Try Jacobi eigen decomposition from Matrix<f64>
        if let Ok((values, _)) = mat.jacobi_eigen_decomposition(2000, tolerance.max(1e-12)) {
            let eig_rows: Vec<Vec<Expr>> = values.into_iter().map(|v| vec![Expr::Constant(v)]).collect();
            return Some(Expr::Matrix(eig_rows));
        }

        // Option C: Try nalgebra symmetric eigenvalues
        let dmat = nalgebra::DMatrix::from_row_slice(rows, cols, data);
        let values = dmat.symmetric_eigenvalues();
        let eig_rows: Vec<Vec<Expr>> = values.iter().map(|&v| vec![Expr::Constant(v)]).collect();
        return Some(Expr::Matrix(eig_rows));
    } else {
        // Non-symmetric real matrix: use nalgebra complex_eigenvalues
        let dmat = nalgebra::DMatrix::from_row_slice(rows, cols, data);
        let complex_eigs = dmat.complex_eigenvalues();
        let mut eig_rows = Vec::with_capacity(rows);
        for eig in complex_eigs.iter() {
            let elem = if eig.im.abs() <= 1e-12 {
                Expr::Constant(eig.re)
            } else {
                Expr::new_complex(Expr::Constant(eig.re), Expr::Constant(eig.im))
            };
            eig_rows.push(vec![elem]);
        }
        return Some(Expr::Matrix(eig_rows));
    }
}

impl Rule for MatrixOracleRule {
    fn name(&self) -> &str {
        "oracles::matrix_operations"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut inv_targets = Vec::new();
        let mut det_targets = Vec::new();
        let mut trace_targets = Vec::new();
        let mut transpose_targets = Vec::new();
        let mut mul_targets = Vec::new();
        let mut oracle_call_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Inverse(mat) => inv_targets.push((class_id, *mat)),
                    ENode::Determinant(mat) => det_targets.push((class_id, *mat)),
                    ENode::Trace(mat) => trace_targets.push((class_id, *mat)),
                    ENode::Transpose(mat) => transpose_targets.push((class_id, *mat)),
                    ENode::MatrixMul(m1, m2) => mul_targets.push((class_id, *m1, *m2)),
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "solve_linear_system"
                            | "characteristic_polynomial"
                            | "charpoly"
                            | "eigenvalues"
                            | "rref"
                            | "null_space" => {
                                oracle_call_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if inv_targets.is_empty()
            && det_targets.is_empty()
            && trace_targets.is_empty()
            && transpose_targets.is_empty()
            && mul_targets.is_empty()
            && oracle_call_targets.is_empty()
        {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, mat_id) in inv_targets {
            let mat_expr = extractor.extract(egraph, mat_id);
            if matches!(mat_expr, Expr::Matrix(..)) {
                if self.is_numerical_target() {
                    if let Some((rows, cols, data)) = extract_f64_matrix(&mat_expr, &self.bindings) {
                        if let Some(inv) = numerical_matrix_inverse(rows, cols, &data) {
                            results.push((class_id, inv));
                            continue;
                        }
                    }
                }
                let inv = crate::symbolic::matrix::inverse_internal(&mat_expr);
                if !matches!(inv, Expr::Inverse(..)) {
                    results.push((class_id, inv));
                }
            }
        }

        for (class_id, mat_id) in det_targets {
            let mat_expr = extractor.extract(egraph, mat_id);
            if matches!(mat_expr, Expr::Matrix(..)) {
                if self.is_numerical_target() {
                    if let Some((rows, cols, data)) = extract_f64_matrix(&mat_expr, &self.bindings) {
                        if let Some(det_val) = numerical_determinant(rows, cols, &data) {
                            results.push((class_id, Expr::Constant(det_val)));
                            continue;
                        }
                    }
                }
                let det = crate::symbolic::matrix::determinant_internal(&mat_expr);
                results.push((class_id, det));
            }
        }

        for (class_id, mat_id) in trace_targets {
            let mat_expr = extractor.extract(egraph, mat_id);
            if matches!(mat_expr, Expr::Matrix(..)) {
                if let Ok(tr) = crate::symbolic::matrix::trace_internal(&mat_expr) {
                    results.push((class_id, tr));
                }
            }
        }

        for (class_id, mat_id) in transpose_targets {
            let mat_expr = extractor.extract(egraph, mat_id);
            if matches!(mat_expr, Expr::Matrix(..)) {
                let tr = crate::symbolic::matrix::transpose_internal(&mat_expr);
                results.push((class_id, tr));
            }
        }

        for (class_id, m1_id, m2_id) in mul_targets {
            let m1_expr = extractor.extract(egraph, m1_id);
            let m2_expr = extractor.extract(egraph, m2_id);
            if matches!(m1_expr, Expr::Matrix(..)) && matches!(m2_expr, Expr::Matrix(..)) {
                let prod = crate::symbolic::matrix::mul_matrices_internal(&m1_expr, &m2_expr);
                results.push((class_id, prod));
            }
        }

        for (class_id, name, args) in oracle_call_targets {
            match name.as_str() {
                "solve_linear_system" if args.len() >= 2 => {
                    let a = extractor.extract(egraph, args[0]);
                    let b = extractor.extract(egraph, args[1]);
                    if let Ok(sol) = crate::symbolic::matrix::solve_linear_system_internal(&a, &b) {
                        results.push((class_id, sol));
                    } else {
                        results.push((class_id, Expr::NoSolution));
                    }
                }
                "characteristic_polynomial" | "charpoly" if !args.is_empty() => {
                    let m = extractor.extract(egraph, args[0]);
                    let var = if args.len() >= 2 {
                        match extractor.extract(egraph, args[1]) {
                            Expr::Variable(v) => v,
                            other => format!("{other}"),
                        }
                    } else {
                        "lambda".to_string()
                    };
                    if let Ok(poly) = crate::symbolic::matrix::characteristic_polynomial_internal(&m, &var) {
                        results.push((class_id, poly));
                    } else {
                        results.push((class_id, Expr::NoSolution));
                    }
                }
                "eigenvalues" if !args.is_empty() => {
                    let m = extractor.extract(egraph, args[0]);
                    if matches!(m, Expr::Matrix(..)) {
                        if self.is_numerical_target() {
                            if let Some((rows, cols, data)) = extract_f64_matrix(&m, &self.bindings) {
                                if let Some(eigs) = numerical_eigenvalues(rows, cols, &data, self.numerical_tolerance()) {
                                    results.push((class_id, eigs));
                                    continue;
                                }
                            }
                        }
                        // Symbolic fallback / target
                        if let Ok((eigs, _)) = crate::symbolic::matrix::eigen_decomposition(&m) {
                            results.push((class_id, eigs));
                        } else if let Some((rows, cols, data)) = extract_f64_matrix(&m, &self.bindings) {
                            // If symbolic failed, fallback to numerical
                            if let Some(eigs) = numerical_eigenvalues(rows, cols, &data, self.numerical_tolerance()) {
                                results.push((class_id, eigs));
                            } else {
                                results.push((class_id, Expr::NoSolution));
                            }
                        } else {
                            results.push((class_id, Expr::NoSolution));
                        }
                    }
                }
                "rref" if !args.len().is_zero() => {
                    let m = extractor.extract(egraph, args[0]);
                    if let Ok(r) = crate::symbolic::matrix::rref_internal(&m) {
                        results.push((class_id, r));
                    } else {
                        results.push((class_id, Expr::NoSolution));
                    }
                }
                "null_space" if !args.len().is_zero() => {
                    let m = extractor.extract(egraph, args[0]);
                    if let Ok(ns) = crate::symbolic::matrix::null_space_internal(&m) {
                        results.push((class_id, ns));
                    } else {
                        results.push((class_id, Expr::NoSolution));
                    }
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Integral Transforms Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges E-Graph transform operations (Laplace, Fourier, Z-Transform)
/// with `rssn::symbolic::transforms`.
#[derive(Debug)]
pub struct TransformOracleRule;

impl Rule for TransformOracleRule {
    fn name(&self) -> &str {
        "oracles::integral_transforms"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_call_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Laplace { expr, var, target } => {
                        targets.push((class_id, 0, *expr, var.clone(), target.clone()));
                    }
                    ENode::InverseLaplace { expr, var, target } => {
                        targets.push((class_id, 1, *expr, var.clone(), target.clone()));
                    }
                    ENode::Fourier { expr, var, target } => {
                        targets.push((class_id, 2, *expr, var.clone(), target.clone()));
                    }
                    ENode::InverseFourier { expr, var, target } => {
                        targets.push((class_id, 3, *expr, var.clone(), target.clone()));
                    }
                    ENode::ZTransform { expr, var, target } => {
                        targets.push((class_id, 4, *expr, var.clone(), target.clone()));
                    }
                    ENode::InverseZTransform { expr, var, target } => {
                        targets.push((class_id, 5, *expr, var.clone(), target.clone()));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "partial_fraction_decomposition" if args.len() >= 2 => {
                                oracle_call_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_call_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, op_kind, expr_id, var, target) in targets {
            let expr = extractor.extract(egraph, expr_id);
            let res = match op_kind {
                0 => crate::symbolic::transforms::laplace_transform_internal(&expr, &var, &target),
                1 => crate::symbolic::transforms::inverse_laplace_transform_internal(&expr, &var, &target),
                2 => crate::symbolic::transforms::fourier_transform_internal(&expr, &var, &target),
                3 => crate::symbolic::transforms::inverse_fourier_transform_internal(&expr, &var, &target),
                4 => crate::symbolic::transforms::z_transform_internal(&expr, &var, &target),
                5 => crate::symbolic::transforms::inverse_z_transform_internal(&expr, &var, &target),
                _ => continue,
            };
            results.push((class_id, res));
        }

        for (class_id, name, args) in oracle_call_targets {
            match name.as_str() {
                "partial_fraction_decomposition" if args.len() >= 2 => {
                    let expr = extractor.extract(egraph, args[0]);
                    let var = match extractor.extract(egraph, args[1]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    if let Some(terms) = crate::symbolic::transforms::partial_fraction_decomposition_internal(&expr, &var) {
                        results.push((class_id, Expr::Tuple(terms)));
                    } else {
                        results.push((class_id, Expr::NoSolution));
                    }
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Vector Calculus Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges Gradient, Divergence, Curl, Laplacian with symbolic differentiation.
#[derive(Debug)]
pub struct VectorCalculusOracleRule;

impl Rule for VectorCalculusOracleRule {
    fn name(&self) -> &str {
        "oracles::vector_calculus"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Gradient { expr, vars } => {
                        targets.push((class_id, 0, *expr, vars.clone()));
                    }
                    ENode::Divergence { expr, vars } => {
                        targets.push((class_id, 1, *expr, vars.clone()));
                    }
                    ENode::Curl { expr, vars } => {
                        targets.push((class_id, 2, *expr, vars.clone()));
                    }
                    ENode::Laplacian { expr, vars } => {
                        targets.push((class_id, 3, *expr, vars.clone()));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "line_integral_scalar" if args.len() >= 5 => {
                                oracle_targets.push((class_id, name.clone(), args.clone()));
                            }
                            "line_integral_vector" if args.len() >= 5 => {
                                oracle_targets.push((class_id, name.clone(), args.clone()));
                            }
                            "surface_integral" if args.len() >= 8 => {
                                oracle_targets.push((class_id, name.clone(), args.clone()));
                            }
                            "volume_integral" if args.len() >= 10 => {
                                oracle_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, op_kind, expr_id, vars) in targets {
            let expr = extractor.extract(egraph, expr_id);
            match op_kind {
                0 => {
                    let components: Vec<Expr> = vars
                        .iter()
                        .map(|v| crate::symbolic::egraph::diff(&expr, v))
                        .collect();
                    results.push((class_id, Expr::Vector(components)));
                }
                1 => {
                    if let Expr::Vector(comps) = &expr {
                        let mut sum = Expr::BigInt(BigInt::zero());
                        for (i, v) in vars.iter().enumerate() {
                            if i < comps.len() {
                                let df = crate::symbolic::egraph::diff(&comps[i], v);
                                sum = Expr::new_add(sum, df);
                            }
                        }
                        results.push((class_id, sum));
                    }
                }
                2 => {
                    if let Expr::Vector(comps) = &expr {
                        if comps.len() >= 3 && vars.len() >= 3 {
                            let fx = &comps[0];
                            let fy = &comps[1];
                            let fz = &comps[2];
                            let x_comp = Expr::new_sub(
                                crate::symbolic::egraph::diff(fz, &vars[1]),
                                crate::symbolic::egraph::diff(fy, &vars[2]),
                            );
                            let y_comp = Expr::new_sub(
                                crate::symbolic::egraph::diff(fx, &vars[2]),
                                crate::symbolic::egraph::diff(fz, &vars[0]),
                            );
                            let z_comp = Expr::new_sub(
                                crate::symbolic::egraph::diff(fy, &vars[0]),
                                crate::symbolic::egraph::diff(fx, &vars[1]),
                            );
                            results.push((class_id, Expr::Vector(vec![x_comp, y_comp, z_comp])));
                        }
                    }
                }
                3 => {
                    let mut sum = Expr::BigInt(BigInt::zero());
                    for v in &vars {
                        let d1 = crate::symbolic::egraph::diff(&expr, v);
                        let d2 = crate::symbolic::egraph::diff(&d1, v);
                        sum = Expr::new_add(sum, d2);
                    }
                    results.push((class_id, sum));
                }
                _ => {}
            }
        }

        for (class_id, name, args) in oracle_targets {
            match name.as_str() {
                "line_integral_scalar" => {
                    let field = extractor.extract(egraph, args[0]);
                    let r_expr = extractor.extract(egraph, args[1]);
                    let r_vec = match r_expr {
                        Expr::Vector(v) => v,
                        _ => vec![Expr::Constant(0.0), Expr::Constant(0.0), Expr::Constant(0.0)],
                    };
                    let rx = r_vec.get(0).cloned().unwrap_or(Expr::Constant(0.0));
                    let ry = r_vec.get(1).cloned().unwrap_or(Expr::Constant(0.0));
                    let rz = r_vec.get(2).cloned().unwrap_or(Expr::Constant(0.0));
                    let t_var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let t0 = extractor.extract(egraph, args[3]);
                    let t1 = extractor.extract(egraph, args[4]);
                    let curve = crate::symbolic::vector_calculus::ParametricCurve {
                        r: crate::symbolic::vector::Vector::new(rx, ry, rz),
                        t_var,
                        t_bounds: (t0, t1),
                    };
                    let res = crate::symbolic::vector_calculus::line_integral_scalar_internal(&field, &curve);
                    results.push((class_id, res));
                }
                "line_integral_vector" => {
                    let f_expr = extractor.extract(egraph, args[0]);
                    let f_vec = match f_expr {
                        Expr::Vector(v) => v,
                        _ => vec![Expr::Constant(0.0), Expr::Constant(0.0), Expr::Constant(0.0)],
                    };
                    let fx = f_vec.get(0).cloned().unwrap_or(Expr::Constant(0.0));
                    let fy = f_vec.get(1).cloned().unwrap_or(Expr::Constant(0.0));
                    let fz = f_vec.get(2).cloned().unwrap_or(Expr::Constant(0.0));
                    let field = crate::symbolic::vector::Vector::new(fx, fy, fz);

                    let r_expr = extractor.extract(egraph, args[1]);
                    let r_vec = match r_expr {
                        Expr::Vector(v) => v,
                        _ => vec![Expr::Constant(0.0), Expr::Constant(0.0), Expr::Constant(0.0)],
                    };
                    let rx = r_vec.get(0).cloned().unwrap_or(Expr::Constant(0.0));
                    let ry = r_vec.get(1).cloned().unwrap_or(Expr::Constant(0.0));
                    let rz = r_vec.get(2).cloned().unwrap_or(Expr::Constant(0.0));
                    let t_var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let t0 = extractor.extract(egraph, args[3]);
                    let t1 = extractor.extract(egraph, args[4]);
                    let curve = crate::symbolic::vector_calculus::ParametricCurve {
                        r: crate::symbolic::vector::Vector::new(rx, ry, rz),
                        t_var,
                        t_bounds: (t0, t1),
                    };
                    let res = crate::symbolic::vector_calculus::line_integral_vector_internal(&field, &curve);
                    results.push((class_id, res));
                }
                "surface_integral" => {
                    let f_expr = extractor.extract(egraph, args[0]);
                    let f_vec = match f_expr {
                        Expr::Vector(v) => v,
                        _ => vec![Expr::Constant(0.0), Expr::Constant(0.0), Expr::Constant(0.0)],
                    };
                    let fx = f_vec.get(0).cloned().unwrap_or(Expr::Constant(0.0));
                    let fy = f_vec.get(1).cloned().unwrap_or(Expr::Constant(0.0));
                    let fz = f_vec.get(2).cloned().unwrap_or(Expr::Constant(0.0));
                    let field = crate::symbolic::vector::Vector::new(fx, fy, fz);

                    let r_expr = extractor.extract(egraph, args[1]);
                    let r_vec = match r_expr {
                        Expr::Vector(v) => v,
                        _ => vec![Expr::Constant(0.0), Expr::Constant(0.0), Expr::Constant(0.0)],
                    };
                    let rx = r_vec.get(0).cloned().unwrap_or(Expr::Constant(0.0));
                    let ry = r_vec.get(1).cloned().unwrap_or(Expr::Constant(0.0));
                    let rz = r_vec.get(2).cloned().unwrap_or(Expr::Constant(0.0));

                    let u_var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let u0 = extractor.extract(egraph, args[3]);
                    let u1 = extractor.extract(egraph, args[4]);
                    let v_var = match extractor.extract(egraph, args[5]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let v0 = extractor.extract(egraph, args[6]);
                    let v1 = extractor.extract(egraph, args[7]);

                    let surface = crate::symbolic::vector_calculus::ParametricSurface {
                        r: crate::symbolic::vector::Vector::new(rx, ry, rz),
                        u_var,
                        u_bounds: (u0, u1),
                        v_var,
                        v_bounds: (v0, v1),
                    };
                    let res = crate::symbolic::vector_calculus::surface_integral_internal(&field, &surface);
                    results.push((class_id, res));
                }
                "volume_integral" => {
                    let field = extractor.extract(egraph, args[0]);
                    let x_var = match extractor.extract(egraph, args[1]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let y_var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let z_var = match extractor.extract(egraph, args[3]) {
                        Expr::Variable(v) => v,
                        other => format!("{other}"),
                    };
                    let x0 = extractor.extract(egraph, args[4]);
                    let x1 = extractor.extract(egraph, args[5]);
                    let y0 = extractor.extract(egraph, args[6]);
                    let y1 = extractor.extract(egraph, args[7]);
                    let z0 = extractor.extract(egraph, args[8]);
                    let z1 = extractor.extract(egraph, args[9]);

                    let volume = crate::symbolic::vector_calculus::Volume {
                        vars: (x_var, y_var, z_var),
                        x_bounds: (x0, x1),
                        y_bounds: (y0, y1),
                        z_bounds: (z0, z1),
                    };
                    let res = crate::symbolic::vector_calculus::volume_integral_internal(&field, &volume);
                    results.push((class_id, res));
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Calculus of Variations Euler-Lagrange Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges EulerLagrange node with symbolic derivation: d/dt(dL/dq') - dL/dq.
#[derive(Debug)]
pub struct EulerLagrangeOracleRule;

impl Rule for EulerLagrangeOracleRule {
    fn name(&self) -> &str {
        "oracles::euler_lagrange"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut solve_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::EulerLagrange { lagrangian, func, var } => {
                        targets.push((class_id, *lagrangian, func.clone(), var.clone()));
                    }
                    ENode::OracleCall(name, args) if name == "euler_lagrange" && args.len() >= 3 => {
                        let extractor = Extractor::new(egraph);
                        let func_expr = extractor.extract(egraph, args[1]);
                        let var_expr = extractor.extract(egraph, args[2]);
                        let func = match func_expr { Expr::Variable(s) => s, _ => func_expr.to_string() };
                        let var = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
                        targets.push((class_id, args[0], func, var));
                    }
                    ENode::OracleCall(name, args) if name == "solve_euler_lagrange" && args.len() >= 3 => {
                        solve_targets.push((class_id, args[0], args[1], args[2]));
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && solve_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, lagrangian_id, func, var) in targets {
            let l_expr = extractor.extract(egraph, lagrangian_id);
            let result = crate::symbolic::calculus_of_variations::euler_lagrange_internal(&l_expr, &func, &var);
            results.push((class_id, result));
        }

        for (class_id, l_id, func_id, var_id) in solve_targets {
            let l_expr = extractor.extract(egraph, l_id);
            let func_expr = extractor.extract(egraph, func_id);
            let var_expr = extractor.extract(egraph, var_id);
            let func_str = match func_expr { Expr::Variable(s) => s, _ => func_expr.to_string() };
            let var_str = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
            let result = crate::symbolic::calculus_of_variations::solve_euler_lagrange_internal(&l_expr, &func_str, &var_str);
            results.push((class_id, result));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Quantum Mechanics Operator Commutator Oracle Rule (Tier 0).
/// Rewrites [A, B] <=> A*B - B*A and {A, B} <=> A*B + B*A.
#[derive(Debug)]
pub struct QuantumOracleRule;

impl Rule for QuantumOracleRule {
    fn name(&self) -> &str {
        "oracles::quantum_operators"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_calls = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Commutator(a, b) => targets.push((class_id, true, *a, *b)),
                    ENode::Anticommutator(a, b) => targets.push((class_id, false, *a, *b)),
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "quantum_commutator" | "quantum_bra_ket" | "quantum_expectation_value" => {
                                oracle_calls.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut applied = 0;
        for (class_id, is_comm, a, b) in targets {
            let ab = egraph.add_node(ENode::Mul(a, b));
            let ba = egraph.add_node(ENode::Mul(b, a));
            let res = if is_comm {
                egraph.add_node(ENode::Sub(ab, ba))
            } else {
                egraph.add_node(ENode::Add(ab, ba))
            };
            if egraph.union(class_id, res) {
                applied += 1;
            }
        }

        if !oracle_calls.is_empty() {
            let extractor = Extractor::new(egraph);
            let mut results = Vec::new();
            for (class_id, name, args) in oracle_calls {
                match name.as_str() {
                    "quantum_commutator" if args.len() >= 3 => {
                        let a_expr = extractor.extract(egraph, args[0]);
                        let b_expr = extractor.extract(egraph, args[1]);
                        let ket_expr = extractor.extract(egraph, args[2]);
                        let a_op = crate::symbolic::quantum_mechanics::Operator::new(a_expr);
                        let b_op = crate::symbolic::quantum_mechanics::Operator::new(b_expr);
                        let ket = crate::symbolic::quantum_mechanics::Ket { state: ket_expr };
                        let res = crate::symbolic::quantum_mechanics::commutator_internal(&a_op, &b_op, &ket);
                        results.push((class_id, res));
                    }
                    "quantum_bra_ket" if args.len() >= 2 => {
                        let bra_expr = extractor.extract(egraph, args[0]);
                        let ket_expr = extractor.extract(egraph, args[1]);
                        let bra = crate::symbolic::quantum_mechanics::Bra { state: bra_expr };
                        let ket = crate::symbolic::quantum_mechanics::Ket { state: ket_expr };
                        let res = crate::symbolic::quantum_mechanics::bra_ket_internal(&bra, &ket);
                        results.push((class_id, res));
                    }
                    "quantum_expectation_value" if args.len() >= 2 => {
                        let op_expr = extractor.extract(egraph, args[0]);
                        let psi_expr = extractor.extract(egraph, args[1]);
                        let op = crate::symbolic::quantum_mechanics::Operator::new(op_expr);
                        let psi = crate::symbolic::quantum_mechanics::Ket { state: psi_expr };
                        let res = crate::symbolic::quantum_mechanics::expectation_value_internal(&op, &psi);
                        results.push((class_id, res));
                    }
                    _ => {}
                }
            }
            for (class_id, res_expr) in results {
                let res_id = egraph.add_expr(&res_expr);
                if egraph.union(class_id, res_id) {
                    applied += 1;
                }
            }
        }

        applied
    }
}

/// Complex Analysis Residue Oracle Rule (Tier 0: High-weight De-cocooning).
/// Bridges Residue(f, z, z0) with lim_{z -> z0} (z - z0)*f(z).
#[derive(Debug)]
pub struct ResidueOracleRule;

impl Rule for ResidueOracleRule {
    fn name(&self) -> &str {
        "oracles::residue"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();
        let mut oracle_calls = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Residue { expr, var, point } => {
                        targets.push((class_id, *expr, var.clone(), *point));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "residue"
                            | "contour_integral_residue_theorem"
                            | "cauchy_integral_formula"
                            | "cauchy_derivative_formula" => {
                                oracle_calls.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        if targets.is_empty() && oracle_calls.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, expr_id, var, point_id) in targets {
            let f = extractor.extract(egraph, expr_id);
            let pt = extractor.extract(egraph, point_id);
            let z = Expr::Variable(var.clone());
            let factor = Expr::new_sub(z, pt.clone());
            let product = Expr::new_mul(factor, f);
            let res = crate::symbolic::egraph::limit(&product, &var, &pt);
            results.push((class_id, res));
        }

        for (class_id, name, args) in oracle_calls {
            match name.as_str() {
                "residue" if args.len() >= 3 => {
                    let f = extractor.extract(egraph, args[0]);
                    let var_expr = extractor.extract(egraph, args[1]);
                    let var = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
                    let pt = extractor.extract(egraph, args[2]);
                    let res = crate::symbolic::complex_analysis::calculate_residue_internal(&f, &var, &pt);
                    results.push((class_id, res));
                }
                "contour_integral_residue_theorem" if args.len() >= 2 => {
                    let f = extractor.extract(egraph, args[0]);
                    let var_expr = extractor.extract(egraph, args[1]);
                    let var = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
                    let singularities: Vec<Expr> = args[2..].iter().map(|&id| extractor.extract(egraph, id)).collect();
                    let res = crate::symbolic::complex_analysis::contour_integral_residue_theorem_internal(&f, &var, &singularities);
                    results.push((class_id, res));
                }
                "cauchy_integral_formula" if args.len() >= 3 => {
                    let f = extractor.extract(egraph, args[0]);
                    let var_expr = extractor.extract(egraph, args[1]);
                    let var = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
                    let z0 = extractor.extract(egraph, args[2]);
                    let res = crate::symbolic::complex_analysis::cauchy_integral_formula_internal(&f, &var, &z0);
                    results.push((class_id, res));
                }
                "cauchy_derivative_formula" if args.len() >= 4 => {
                    let f = extractor.extract(egraph, args[0]);
                    let var_expr = extractor.extract(egraph, args[1]);
                    let var = match var_expr { Expr::Variable(s) => s, _ => var_expr.to_string() };
                    let z0 = extractor.extract(egraph, args[2]);
                    let n = extractor.extract(egraph, args[3]).to_f64().unwrap_or(1.0) as usize;
                    let res = crate::symbolic::complex_analysis::cauchy_derivative_formula_internal(&f, &var, &z0, n);
                    results.push((class_id, res));
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

fn gcd_bigint(mut a: BigInt, mut b: BigInt) -> BigInt {
    while !b.is_zero() {
        let r = &a % &b;
        a = b;
        b = r;
    }
    if a < BigInt::zero() { -a } else { a }
}

fn factorial_bigint(n: u64) -> BigInt {
    let mut res = BigInt::one();
    for i in 2..=n {
        res *= BigInt::from(i);
    }
    res
}

fn binomial_bigint(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    if k == 0 || k == n {
        return BigInt::one();
    }
    let k = k.min(n - k);
    let mut num = BigInt::one();
    let mut den = BigInt::one();
    for i in 0..k {
        num *= BigInt::from(n - i);
        den *= BigInt::from(i + 1);
    }
    num / den
}

fn permutation_bigint(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let mut res = BigInt::one();
    for i in 0..k {
        res *= BigInt::from(n - i);
    }
    res
}

/// Combinatorics & Arithmetic Oracle Rule (Tier 0: Numeric and Algebraic Reductions).
/// Evaluates Factorial, Binomial, Permutation, Combination, Gcd, Lcm, Mod, Floor, Max.
#[derive(Debug)]
pub struct CombinatoricsOracleRule;

impl Rule for CombinatoricsOracleRule {
    fn name(&self) -> &str {
        "oracles::combinatorics"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut fact_targets = Vec::new();
        let mut binom_targets = Vec::new();
        let mut perm_targets = Vec::new();
        let mut comb_targets = Vec::new();
        let mut gcd_targets = Vec::new();
        let mut lcm_targets = Vec::new();
        let mut mod_targets = Vec::new();
        let mut floor_targets = Vec::new();
        let mut max_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Factorial(n) => fact_targets.push((class_id, *n)),
                    ENode::Binomial(n, k) => binom_targets.push((class_id, *n, *k)),
                    ENode::Permutation(n, k) => perm_targets.push((class_id, *n, *k)),
                    ENode::Combination(n, k) => comb_targets.push((class_id, *n, *k)),
                    ENode::Gcd(a, b) => gcd_targets.push((class_id, *a, *b)),
                    ENode::Lcm(a, b) => lcm_targets.push((class_id, *a, *b)),
                    ENode::Mod(a, b) => mod_targets.push((class_id, *a, *b)),
                    ENode::Floor(a) => floor_targets.push((class_id, *a)),
                    ENode::Max(a, b) => max_targets.push((class_id, *a, *b)),
                    _ => {}
                }
            }
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, n_id) in fact_targets {
            let n_expr = extractor.extract(egraph, n_id);
            if let Some(n_val) = n_expr.to_f64() {
                if n_val >= 0.0 && n_val.fract() == 0.0 && n_val <= 30.0 {
                    let k = n_val as u64;
                    let f = factorial_bigint(k);
                    results.push((class_id, Expr::BigInt(f)));
                }
            }
        }

        for (class_id, n_id, k_id) in binom_targets {
            let n_expr = extractor.extract(egraph, n_id);
            let k_expr = extractor.extract(egraph, k_id);
            if let (Some(n_val), Some(k_val)) = (n_expr.to_f64(), k_expr.to_f64()) {
                if n_val >= 0.0 && k_val >= 0.0 && n_val.fract() == 0.0 && k_val.fract() == 0.0 && n_val <= 100.0 {
                    let b = binomial_bigint(n_val as u64, k_val as u64);
                    results.push((class_id, Expr::BigInt(b)));
                }
            }
        }

        let mut symbolic_perms = Vec::new();
        let mut symbolic_combs = Vec::new();

        for (class_id, n_id, k_id) in perm_targets {
            let n_expr = extractor.extract(egraph, n_id);
            let k_expr = extractor.extract(egraph, k_id);
            if let (Some(n_val), Some(k_val)) = (n_expr.to_f64(), k_expr.to_f64()) {
                if n_val >= 0.0 && k_val >= 0.0 && n_val.fract() == 0.0 && k_val.fract() == 0.0 && n_val <= 30.0 {
                    let p = permutation_bigint(n_val as u64, k_val as u64);
                    results.push((class_id, Expr::BigInt(p)));
                    continue;
                }
            }
            symbolic_perms.push((class_id, n_id, k_id));
        }

        for (class_id, n_id, k_id) in comb_targets {
            let n_expr = extractor.extract(egraph, n_id);
            let k_expr = extractor.extract(egraph, k_id);
            if let (Some(n_val), Some(k_val)) = (n_expr.to_f64(), k_expr.to_f64()) {
                if n_val >= 0.0 && k_val >= 0.0 && n_val.fract() == 0.0 && k_val.fract() == 0.0 && n_val <= 100.0 {
                    let b = binomial_bigint(n_val as u64, k_val as u64);
                    results.push((class_id, Expr::BigInt(b)));
                    continue;
                }
            }
            symbolic_combs.push((class_id, n_id, k_id));
        }

        for (class_id, a_id, b_id) in gcd_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);
            if let (Some(a_val), Some(b_val)) = (a_expr.to_f64(), b_expr.to_f64()) {
                if a_val.fract() == 0.0 && b_val.fract() == 0.0 {
                    let g = gcd_bigint(BigInt::from(a_val as i64), BigInt::from(b_val as i64));
                    results.push((class_id, Expr::BigInt(g)));
                }
            }
        }

        for (class_id, a_id, b_id) in lcm_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);
            if let (Some(a_val), Some(b_val)) = (a_expr.to_f64(), b_expr.to_f64()) {
                if a_val.fract() == 0.0 && b_val.fract() == 0.0 {
                    let a_bi = BigInt::from(a_val as i64);
                    let b_bi = BigInt::from(b_val as i64);
                    let g = gcd_bigint(a_bi.clone(), b_bi.clone());
                    if !g.is_zero() {
                        let lcm = (a_bi * b_bi).abs() / g;
                        results.push((class_id, Expr::BigInt(lcm)));
                    }
                }
            }
        }

        for (class_id, a_id, b_id) in mod_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);
            if let (Some(a_val), Some(b_val)) = (a_expr.to_f64(), b_expr.to_f64()) {
                if a_val.fract() == 0.0 && b_val.fract() == 0.0 && b_val != 0.0 {
                    let a_i = a_val as i64;
                    let b_i = b_val as i64;
                    let m = a_i.rem_euclid(b_i);
                    if a_i != m {
                        results.push((class_id, Expr::Mod(Arc::new(Expr::BigInt(BigInt::from(m))), Arc::new(b_expr))));
                    }
                }
            }
        }

        for (class_id, a_id) in floor_targets {
            let a_expr = extractor.extract(egraph, a_id);
            if let Some(a_val) = a_expr.to_f64() {
                results.push((class_id, Expr::BigInt(BigInt::from(a_val.floor() as i64))));
            }
        }

        for (class_id, a_id, b_id) in max_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);
            if let (Some(a_val), Some(b_val)) = (a_expr.to_f64(), b_expr.to_f64()) {
                let m = if a_val >= b_val { a_expr } else { b_expr };
                results.push((class_id, m));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }

        for (class_id, n_id, k_id) in symbolic_perms {
            let n_fact = egraph.add_node(ENode::Factorial(n_id));
            let sub = egraph.add_node(ENode::Sub(n_id, k_id));
            let sub_fact = egraph.add_node(ENode::Factorial(sub));
            let div = egraph.add_node(ENode::Div(n_fact, sub_fact));
            if egraph.union(class_id, div) {
                applied += 1;
            }
        }

        for (class_id, n_id, k_id) in symbolic_combs {
            let n_fact = egraph.add_node(ENode::Factorial(n_id));
            let k_fact = egraph.add_node(ENode::Factorial(k_id));
            let sub = egraph.add_node(ENode::Sub(n_id, k_id));
            let sub_fact = egraph.add_node(ENode::Factorial(sub));
            let denom = egraph.add_node(ENode::Mul(k_fact, sub_fact));
            let div = egraph.add_node(ENode::Div(n_fact, denom));
            if egraph.union(class_id, div) {
                applied += 1;
            }
        }

        applied
    }
}

/// Boolean Logic Oracle Rule (Tier 0: Propositional Simplifications).
/// Implication, equivalence, XOR, De Morgan, double negation, excluded middle, contradiction, and truth table reduction.
#[derive(Debug)]
pub struct LogicOracleRule;

impl Rule for LogicOracleRule {
    fn name(&self) -> &str {
        "oracles::logic"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut not_targets = Vec::new();
        let mut and_targets = Vec::new();
        let mut or_targets = Vec::new();
        let mut implies_targets = Vec::new();
        let mut equiv_targets = Vec::new();
        let mut xor_targets = Vec::new();
        let mut forall_targets = Vec::new();
        let mut exists_targets = Vec::new();
        let mut oracle_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Not(a) => not_targets.push((class_id, *a)),
                    ENode::And(list) => and_targets.push((class_id, list.clone())),
                    ENode::Or(list) => or_targets.push((class_id, list.clone())),
                    ENode::Implies(a, b) => implies_targets.push((class_id, *a, *b)),
                    ENode::Equivalent(a, b) => equiv_targets.push((class_id, *a, *b)),
                    ENode::Xor(a, b) => xor_targets.push((class_id, *a, *b)),
                    ENode::ForAll(var, body) => forall_targets.push((class_id, var.clone(), *body)),
                    ENode::Exists(var, body) => exists_targets.push((class_id, var.clone(), *body)),
                    ENode::OracleCall(name, args) => {
                        if (name == "to_cnf" || name == "to_dnf") && !args.is_empty() {
                            oracle_targets.push((class_id, name.clone(), args[0]));
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut applied = 0;

        // Implication: A => B  -->  Not(A) Or B
        for (class_id, a, b) in implies_targets {
            let not_a = egraph.add_node(ENode::Not(a));
            let or_node = egraph.add_node(ENode::Or(vec![not_a, b]));
            if egraph.union(class_id, or_node) {
                applied += 1;
            }
        }

        // Equivalence: A <=> B  -->  (Not(A) Or B) And (Not(B) Or A)
        for (class_id, a, b) in equiv_targets {
            let not_a = egraph.add_node(ENode::Not(a));
            let not_b = egraph.add_node(ENode::Not(b));
            let or1 = egraph.add_node(ENode::Or(vec![not_a, b]));
            let or2 = egraph.add_node(ENode::Or(vec![not_b, a]));
            let and_node = egraph.add_node(ENode::And(vec![or1, or2]));
            if egraph.union(class_id, and_node) {
                applied += 1;
            }
        }

        // Xor: A ^ B  -->  (A Or B) And Not(A And B)
        for (class_id, a, b) in xor_targets {
            let or_node = egraph.add_node(ENode::Or(vec![a, b]));
            let and_node = egraph.add_node(ENode::And(vec![a, b]));
            let not_and = egraph.add_node(ENode::Not(and_node));
            let res = egraph.add_node(ENode::And(vec![or_node, not_and]));
            if egraph.union(class_id, res) {
                applied += 1;
            }
        }

        // Quantifier reduction: ForAll(x, P) -> P if x not free in P
        for (class_id, var, body) in forall_targets {
            let can_body = egraph.union_find.find(body);
            if let Some(c) = egraph.get_class(can_body) {
                if !c.free_vars.contains(&var) {
                    if egraph.union(class_id, can_body) {
                        applied += 1;
                    }
                }
            }
        }

        // Quantifier reduction: Exists(x, P) -> P if x not free in P
        for (class_id, var, body) in exists_targets {
            let can_body = egraph.union_find.find(body);
            if let Some(c) = egraph.get_class(can_body) {
                if !c.free_vars.contains(&var) {
                    if egraph.union(class_id, can_body) {
                        applied += 1;
                    }
                }
            }
        }

        // Not rewrites
        for (class_id, child_id) in not_targets {
            let can_child = egraph.union_find.find(child_id);
            let child_class = match egraph.get_class(can_child) {
                Some(c) => c,
                None => continue,
            };
            for node in child_class.nodes.clone() {
                match node {
                    ENode::Boolean(b) => {
                        let opp = egraph.add_node(ENode::Boolean(!b));
                        if egraph.union(class_id, opp) {
                            applied += 1;
                        }
                    }
                    ENode::Not(inner) => {
                        if egraph.union(class_id, inner) {
                            applied += 1;
                        }
                    }
                    ENode::ForAll(var, body) => {
                        let not_body = egraph.add_node(ENode::Not(body));
                        let exists_node = egraph.add_node(ENode::Exists(var, not_body));
                        if egraph.union(class_id, exists_node) {
                            applied += 1;
                        }
                    }
                    ENode::Exists(var, body) => {
                        let not_body = egraph.add_node(ENode::Not(body));
                        let forall_node = egraph.add_node(ENode::ForAll(var, not_body));
                        if egraph.union(class_id, forall_node) {
                            applied += 1;
                        }
                    }
                    ENode::And(list) => {
                        let not_terms: Vec<Id> = list.iter().map(|&x| egraph.add_node(ENode::Not(x))).collect();
                        let or_node = egraph.add_node(ENode::Or(not_terms));
                        if egraph.union(class_id, or_node) {
                            applied += 1;
                        }
                    }
                    ENode::Or(list) => {
                        let not_terms: Vec<Id> = list.iter().map(|&x| egraph.add_node(ENode::Not(x))).collect();
                        let and_node = egraph.add_node(ENode::And(not_terms));
                        if egraph.union(class_id, and_node) {
                            applied += 1;
                        }
                    }
                    _ => {}
                }
            }
        }

        // And rewrites: flattening, false propagation, true elimination, contradiction detection
        for (class_id, list) in and_targets {
            let mut flattened = Vec::new();
            for id in &list {
                let can_id = egraph.union_find.find(*id);
                let sub_opt = if let Some(c) = egraph.get_class(can_id) {
                    c.nodes.iter().find_map(|node| {
                        if let ENode::And(sub) = node {
                            Some(sub.clone())
                        } else {
                            None
                        }
                    })
                } else {
                    None
                };

                if let Some(sub) = sub_opt {
                    flattened.extend(sub.into_iter().map(|x| egraph.union_find.find(x)));
                } else {
                    flattened.push(can_id);
                }
            }

            let mut has_false = false;
            let mut filtered = Vec::new();
            for id in flattened {
                let can_id = egraph.union_find.find(id);
                let mut is_true = false;
                if let Some(c) = egraph.get_class(can_id) {
                    for node in &c.nodes {
                        match node {
                            ENode::Boolean(false) => {
                                has_false = true;
                                break;
                            }
                            ENode::Boolean(true) => {
                                is_true = true;
                            }
                            _ => {}
                        }
                    }
                }
                if has_false {
                    break;
                }
                if !is_true && !filtered.contains(&can_id) {
                    filtered.push(can_id);
                }
            }

            if has_false {
                let f = egraph.add_node(ENode::Boolean(false));
                if egraph.union(class_id, f) {
                    applied += 1;
                }
                continue;
            }

            // Check contradiction: P and Not(P)
            let mut has_contradiction = false;
            for &id in &filtered {
                let not_inners: Vec<Id> = if let Some(c) = egraph.get_class(id) {
                    c.nodes
                        .iter()
                        .filter_map(|node| if let ENode::Not(inner) = node { Some(*inner) } else { None })
                        .collect()
                } else {
                    Vec::new()
                };

                for inner in not_inners {
                    let can_inner = egraph.union_find.find(inner);
                    if filtered.contains(&can_inner) {
                        has_contradiction = true;
                        break;
                    }
                }
                if has_contradiction {
                    break;
                }
            }

            if has_contradiction {
                let f = egraph.add_node(ENode::Boolean(false));
                if egraph.union(class_id, f) {
                    applied += 1;
                }
                continue;
            }

            if filtered.is_empty() {
                let t = egraph.add_node(ENode::Boolean(true));
                if egraph.union(class_id, t) {
                    applied += 1;
                }
            } else if filtered.len() == 1 {
                if egraph.union(class_id, filtered[0]) {
                    applied += 1;
                }
            } else if filtered.len() != list.len() {
                let new_and = egraph.add_node(ENode::And(filtered));
                if egraph.union(class_id, new_and) {
                    applied += 1;
                }
            }
        }

        // Or rewrites: flattening, true propagation, false elimination, excluded middle (tautology)
        for (class_id, list) in or_targets {
            let mut flattened = Vec::new();
            for id in &list {
                let can_id = egraph.union_find.find(*id);
                let sub_opt = if let Some(c) = egraph.get_class(can_id) {
                    c.nodes.iter().find_map(|node| {
                        if let ENode::Or(sub) = node {
                            Some(sub.clone())
                        } else {
                            None
                        }
                    })
                } else {
                    None
                };

                if let Some(sub) = sub_opt {
                    flattened.extend(sub.into_iter().map(|x| egraph.union_find.find(x)));
                } else {
                    flattened.push(can_id);
                }
            }

            let mut has_true = false;
            let mut filtered = Vec::new();
            for id in flattened {
                let can_id = egraph.union_find.find(id);
                let mut is_false = false;
                if let Some(c) = egraph.get_class(can_id) {
                    for node in &c.nodes {
                        match node {
                            ENode::Boolean(true) => {
                                has_true = true;
                                break;
                            }
                            ENode::Boolean(false) => {
                                is_false = true;
                            }
                            _ => {}
                        }
                    }
                }
                if has_true {
                    break;
                }
                if !is_false && !filtered.contains(&can_id) {
                    filtered.push(can_id);
                }
            }

            if has_true {
                let t = egraph.add_node(ENode::Boolean(true));
                if egraph.union(class_id, t) {
                    applied += 1;
                }
                continue;
            }

            // Check tautology (law of excluded middle): P or Not(P)
            let mut has_tautology = false;
            for &id in &filtered {
                let not_inners: Vec<Id> = if let Some(c) = egraph.get_class(id) {
                    c.nodes
                        .iter()
                        .filter_map(|node| if let ENode::Not(inner) = node { Some(*inner) } else { None })
                        .collect()
                } else {
                    Vec::new()
                };

                for inner in not_inners {
                    let can_inner = egraph.union_find.find(inner);
                    if filtered.contains(&can_inner) {
                        has_tautology = true;
                        break;
                    }
                }
                if has_tautology {
                    break;
                }
            }

            if has_tautology {
                let t = egraph.add_node(ENode::Boolean(true));
                if egraph.union(class_id, t) {
                    applied += 1;
                }
                continue;
            }

            if filtered.is_empty() {
                let f = egraph.add_node(ENode::Boolean(false));
                if egraph.union(class_id, f) {
                    applied += 1;
                }
            } else if filtered.len() == 1 {
                if egraph.union(class_id, filtered[0]) {
                    applied += 1;
                }
            } else if filtered.len() != list.len() {
                let new_or = egraph.add_node(ENode::Or(filtered));
                if egraph.union(class_id, new_or) {
                    applied += 1;
                }
            }
        }

        if !oracle_targets.is_empty() {
            let extractor = Extractor::new(egraph);
            let mut oracle_results = Vec::new();
            for (class_id, name, arg_id) in oracle_targets {
                let expr = extractor.extract(egraph, arg_id);
                if name == "to_cnf" {
                    let res = crate::symbolic::logic::to_cnf_internal(&expr);
                    oracle_results.push((class_id, res));
                } else if name == "to_dnf" {
                    let res = crate::symbolic::logic::to_dnf_internal(&expr);
                    oracle_results.push((class_id, res));
                }
            }
            for (class_id, res_expr) in oracle_results {
                let res_id = egraph.add_expr(&res_expr);
                if egraph.union(class_id, res_id) {
                    applied += 1;
                }
            }
        }

        applied
    }
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

fn lcm_i64(a: i64, b: i64) -> i64 {
    if a == 0 || b == 0 {
        0
    } else {
        (a / gcd_i64(a, b) * b).abs()
    }
}

/// Number Theory Oracle Rule (Tier 0).
/// Handles Mod, Gcd, Lcm, IsPrime, and number theoretic utilities.
#[derive(Debug)]
pub struct NumberTheoryOracleRule;

impl Rule for NumberTheoryOracleRule {
    fn name(&self) -> &str {
        "oracles::number_theory"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut mod_targets = Vec::new();
        let mut gcd_targets = Vec::new();
        let mut lcm_targets = Vec::new();
        let mut prime_targets = Vec::new();
        let mut oracle_call_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Mod(a, b) => mod_targets.push((class_id, *a, *b)),
                    ENode::Gcd(a, b) => gcd_targets.push((class_id, *a, *b)),
                    ENode::Lcm(a, b) => lcm_targets.push((class_id, *a, *b)),
                    ENode::IsPrime(arg) => {
                        prime_targets.push((class_id, *arg));
                    }
                    ENode::OracleCall(name, args) => {
                        match name.as_str() {
                            "extended_gcd" | "solve_diophantine" | "chinese_remainder" => {
                                oracle_call_targets.push((class_id, name.clone(), args.clone()));
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut applied = 0;
        let extractor = Extractor::new(egraph);

        for (class_id, a_id, b_id) in mod_targets {
            let a_canon = egraph.find(a_id);
            let b_canon = egraph.find(b_id);
            if a_canon == b_canon {
                let zero = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(0.0)));
                if egraph.union(class_id, zero) {
                    applied += 1;
                }
                continue;
            }

            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);

            if let (Some(a_val), Some(b_val)) = (crate::symbolic::simplify::as_f64(&a_expr), crate::symbolic::simplify::as_f64(&b_expr)) {
                if b_val != 0.0 && a_val.fract() == 0.0 && b_val.fract() == 0.0 {
                    let a_int = a_val as i64;
                    let b_int = b_val as i64;
                    let rem = ((a_int % b_int) + b_int) % b_int;
                    let rem_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(rem as f64)));
                    if egraph.union(class_id, rem_node) {
                        applied += 1;
                    }
                }
            }
        }

        for (class_id, a_id, b_id) in gcd_targets {
            let a_canon = egraph.find(a_id);
            let b_canon = egraph.find(b_id);
            if a_canon == b_canon {
                if egraph.union(class_id, a_canon) {
                    applied += 1;
                }
                continue;
            }

            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);

            if let (Some(a_val), Some(b_val)) = (crate::symbolic::simplify::as_f64(&a_expr), crate::symbolic::simplify::as_f64(&b_expr)) {
                if a_val.fract() == 0.0 && b_val.fract() == 0.0 {
                    let a_int = (a_val as i64).abs();
                    let b_int = (b_val as i64).abs();
                    let g = gcd_i64(a_int, b_int);
                    let g_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(g as f64)));
                    if egraph.union(class_id, g_node) {
                        applied += 1;
                    }
                }
            }
        }

        for (class_id, a_id, b_id) in lcm_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);

            if let (Some(a_val), Some(b_val)) = (crate::symbolic::simplify::as_f64(&a_expr), crate::symbolic::simplify::as_f64(&b_expr)) {
                if a_val.fract() == 0.0 && b_val.fract() == 0.0 {
                    let a_int = (a_val as i64).abs();
                    let b_int = (b_val as i64).abs();
                    let l = lcm_i64(a_int, b_int);
                    let l_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(l as f64)));
                    if egraph.union(class_id, l_node) {
                        applied += 1;
                    }
                }
            }
        }

        for (class_id, arg_id) in prime_targets {
            let arg_expr = extractor.extract(egraph, arg_id);
            let prime_res = crate::symbolic::number_theory::is_prime_internal(&arg_expr);
            if matches!(prime_res, Expr::Boolean(_)) {
                let res_id = egraph.add_expr(&prime_res);
                if egraph.union(class_id, res_id) {
                    applied += 1;
                }
            }
        }

        for (class_id, name, args) in oracle_call_targets {
            match name.as_str() {
                "extended_gcd" if args.len() >= 2 => {
                    let a = extractor.extract(egraph, args[0]);
                    let b = extractor.extract(egraph, args[1]);
                    let (g, x, y) = crate::symbolic::number_theory::extended_gcd_internal(&a, &b);
                    let tuple_expr = Expr::Tuple(vec![g, x, y]);
                    let tuple_id = egraph.add_expr(&tuple_expr);
                    if egraph.union(class_id, tuple_id) {
                        applied += 1;
                    }
                }
                "solve_diophantine" if args.len() >= 2 => {
                    let eq = extractor.extract(egraph, args[0]);
                    let vars_expr = extractor.extract(egraph, args[1]);
                    let var_strings: Vec<String> = match vars_expr {
                        Expr::Tuple(list) => list
                            .into_iter()
                            .map(|e| match e {
                                Expr::Variable(v) => v,
                                other => format!("{other}"),
                            })
                            .collect(),
                        Expr::Variable(v) => vec![v],
                        other => vec![format!("{other}")],
                    };
                    let var_refs: Vec<&str> = var_strings.iter().map(String::as_str).collect();
                    if let Ok(sols) = crate::symbolic::number_theory::solve_diophantine_internal(&eq, &var_refs) {
                        let sols_expr = Expr::Solutions(sols);
                        let sols_id = egraph.add_expr(&sols_expr);
                        if egraph.union(class_id, sols_id) {
                            applied += 1;
                        }
                    } else {
                        let no_sol_id = egraph.add_expr(&Expr::NoSolution);
                        if egraph.union(class_id, no_sol_id) {
                            applied += 1;
                        }
                    }
                }
                "chinese_remainder" if !args.len().is_zero() => {
                    let pairs_expr = extractor.extract(egraph, args[0]);
                    let mut congruences = Vec::new();
                    if let Expr::Tuple(list) = pairs_expr {
                        for item in list {
                            if let Expr::Tuple(pair) = item {
                                if pair.len() == 2 {
                                    congruences.push((pair[0].clone(), pair[1].clone()));
                                }
                            }
                        }
                    }
                    if let Some(sol) = crate::symbolic::number_theory::chinese_remainder_internal(&congruences) {
                        let sol_id = egraph.add_expr(&sol);
                        if egraph.union(class_id, sol_id) {
                            applied += 1;
                        }
                    } else {
                        let no_sol_id = egraph.add_expr(&Expr::NoSolution);
                        if egraph.union(class_id, no_sol_id) {
                            applied += 1;
                        }
                    }
                }
                _ => {}
            }
        }

        applied
    }
}

/// Elementary Algebra Oracle Rule (Tier 0).
/// Bridges UnaryList("expand", expr) with `rssn::symbolic::elementary::expand_internal`.
#[derive(Debug)]
pub struct ElementaryOracleRule;

impl Rule for ElementaryOracleRule {
    fn name(&self) -> &str {
        "oracles::elementary_algebra"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Expand(arg) = node {
                    targets.push((class_id, *arg));
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, arg_id) in targets {
            let arg_expr = extractor.extract(egraph, arg_id);
            let expanded = crate::symbolic::elementary::expand_internal(arg_expr);
            results.push((class_id, expanded));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Special Functions Oracle Rule (Tier 0).
/// Evaluates Gamma, Zeta, Erf for special known values.
#[derive(Debug)]
pub struct SpecialFunctionsOracleRule;

impl Rule for SpecialFunctionsOracleRule {
    fn name(&self) -> &str {
        "oracles::special_functions"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut gamma_targets = Vec::new();
        let mut zeta_targets = Vec::new();
        let mut erf_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Gamma(arg) => gamma_targets.push((class_id, *arg)),
                    ENode::Zeta(arg) => zeta_targets.push((class_id, *arg)),
                    ENode::Erf(arg) => erf_targets.push((class_id, *arg)),
                    _ => {}
                }
            }
        }

        let mut applied = 0;
        let extractor = Extractor::new(egraph);

        for (class_id, arg_id) in gamma_targets {
            let arg_expr = extractor.extract(egraph, arg_id);
            if let Some(val) = crate::symbolic::simplify::as_f64(&arg_expr) {
                if val > 0.0 && val.fract() == 0.0 && val <= 20.0 {
                    let n = val as u64;
                    let mut fact = 1u64;
                    for k in 1..n {
                        fact *= k;
                    }
                    let res_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(fact as f64)));
                    if egraph.union(class_id, res_node) {
                        applied += 1;
                    }
                }
            }
        }

        for (class_id, arg_id) in zeta_targets {
            let arg_expr = extractor.extract(egraph, arg_id);
            if let Some(val) = crate::symbolic::simplify::as_f64(&arg_expr) {
                if val == 0.0 {
                    let res_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(-0.5)));
                    if egraph.union(class_id, res_node) {
                        applied += 1;
                    }
                }
            }
        }

        for (class_id, arg_id) in erf_targets {
            let arg_expr = extractor.extract(egraph, arg_id);
            if let Some(val) = crate::symbolic::simplify::as_f64(&arg_expr) {
                if val == 0.0 {
                    let res_node = egraph.add_node(ENode::Constant(ordered_float::OrderedFloat(0.0)));
                    if egraph.union(class_id, res_node) {
                        applied += 1;
                    }
                }
            }
        }

        applied
    }
}

/// Polynomial Algebra Oracle Rule (Tier 0).
/// Handles polynomial GCD, division, and canonical polynomial simplification.
#[derive(Debug)]
pub struct PolynomialOracleRule;

impl Rule for PolynomialOracleRule {
    fn name(&self) -> &str {
        "oracles::polynomial_ops"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut gcd_targets = Vec::new();
        let mut div_targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    ENode::Gcd(a, b) => {
                        gcd_targets.push((class_id, *a, *b));
                    }
                    ENode::OracleCall(name, args) if name == "polynomial_long_division" && args.len() >= 3 => {
                        div_targets.push((class_id, args[0], args[1], args[2]));
                    }
                    _ => {}
                }
            }
        }

        if gcd_targets.is_empty() && div_targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, a_id, b_id) in gcd_targets {
            let a_expr = extractor.extract(egraph, a_id);
            let b_expr = extractor.extract(egraph, b_id);

            let mut vars = std::collections::BTreeSet::new();
            a_expr.pre_order_walk(&mut |e| {
                if let Expr::Variable(v) = e {
                    vars.insert(v.clone());
                }
            });
            b_expr.pre_order_walk(&mut |e| {
                if let Expr::Variable(v) = e {
                    vars.insert(v.clone());
                }
            });

            if let Some(var) = vars.iter().next() {
                let var_slice = [var.as_str()];
                let p1 = crate::symbolic::polynomial::expr_to_sparse_poly(&a_expr, &var_slice);
                let p2 = crate::symbolic::polynomial::expr_to_sparse_poly(&b_expr, &var_slice);
                let gcd_poly = crate::symbolic::polynomial::gcd(p1, p2, var);
                let gcd_expr = crate::symbolic::polynomial::sparse_poly_to_expr(&gcd_poly);
                results.push((class_id, gcd_expr));
            }
        }

        for (class_id, n_id, d_id, var_id) in div_targets {
            let n = extractor.extract(egraph, n_id);
            let d = extractor.extract(egraph, d_id);
            let var = match extractor.extract(egraph, var_id) {
                Expr::Variable(v) => v,
                other => format!("{other}"),
            };
            let (q, r) = crate::symbolic::polynomial::polynomial_long_division_internal(&n, &d, &var);
            results.push((class_id, Expr::Tuple(vec![q, r])));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Optimization Oracle Rule (Tier 0).
/// Handles Hessian matrix calculation.
#[derive(Debug)]
pub struct OptimizationOracleRule;

impl Rule for OptimizationOracleRule {
    fn name(&self) -> &str {
        "oracles::optimization"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::OracleCall(name, args) = node {
                    if name == "hessian_matrix" && !args.is_empty() {
                        targets.push((class_id, args.clone()));
                    }
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, args) in targets {
            let f = extractor.extract(egraph, args[0]);
            let mut vars = Vec::new();
            for &arg_id in &args[1..] {
                let v_expr = extractor.extract(egraph, arg_id);
                match v_expr {
                    Expr::Variable(name) => vars.push(name),
                    other => vars.push(format!("{}", other)),
                }
            }
            let var_refs: Vec<&str> = vars.iter().map(String::as_str).collect();
            let res = crate::symbolic::optimize::hessian_matrix_internal(&f, &var_refs);
            results.push((class_id, res));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Functional Analysis Oracle Rule (Tier 0).
/// Bridges Hilbert space inner products, norms, projections, and Banach space norms.
#[derive(Debug)]
pub struct FunctionalAnalysisOracleRule;

impl Rule for FunctionalAnalysisOracleRule {
    fn name(&self) -> &str {
        "oracles::functional_analysis"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::OracleCall(name, args) = node {
                    match name.as_str() {
                        "hilbert_inner_product" if args.len() >= 5 => {
                            targets.push((class_id, name.clone(), args.clone()));
                        }
                        "hilbert_norm" if args.len() >= 4 => {
                            targets.push((class_id, name.clone(), args.clone()));
                        }
                        "banach_norm" if args.len() >= 5 => {
                            targets.push((class_id, name.clone(), args.clone()));
                        }
                        "hilbert_project" if args.len() >= 5 => {
                            targets.push((class_id, name.clone(), args.clone()));
                        }
                        _ => {}
                    }
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, name, args) in targets {
            match name.as_str() {
                "hilbert_inner_product" => {
                    let f = extractor.extract(egraph, args[0]);
                    let g = extractor.extract(egraph, args[1]);
                    let var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{}", other),
                    };
                    let lower = extractor.extract(egraph, args[3]);
                    let upper = extractor.extract(egraph, args[4]);
                    let space = crate::symbolic::functional_analysis::HilbertSpace {
                        var,
                        lower_bound: lower,
                        upper_bound: upper,
                    };
                    let res = crate::symbolic::functional_analysis::inner_product_internal(&space, &f, &g);
                    results.push((class_id, res));
                }
                "hilbert_norm" => {
                    let f = extractor.extract(egraph, args[0]);
                    let var = match extractor.extract(egraph, args[1]) {
                        Expr::Variable(v) => v,
                        other => format!("{}", other),
                    };
                    let lower = extractor.extract(egraph, args[2]);
                    let upper = extractor.extract(egraph, args[3]);
                    let space = crate::symbolic::functional_analysis::HilbertSpace {
                        var,
                        lower_bound: lower,
                        upper_bound: upper,
                    };
                    let res = crate::symbolic::functional_analysis::norm_internal(&space, &f);
                    results.push((class_id, res));
                }
                "banach_norm" => {
                    let f = extractor.extract(egraph, args[0]);
                    let var = match extractor.extract(egraph, args[1]) {
                        Expr::Variable(v) => v,
                        other => format!("{}", other),
                    };
                    let lower = extractor.extract(egraph, args[2]);
                    let upper = extractor.extract(egraph, args[3]);
                    let p = extractor.extract(egraph, args[4]);
                    let space = crate::symbolic::functional_analysis::BanachSpace {
                        var,
                        lower_bound: lower,
                        upper_bound: upper,
                        p,
                    };
                    let res = crate::symbolic::functional_analysis::banach_norm_internal(&space, &f);
                    results.push((class_id, res));
                }
                "hilbert_project" => {
                    let f = extractor.extract(egraph, args[0]);
                    let g = extractor.extract(egraph, args[1]);
                    let var = match extractor.extract(egraph, args[2]) {
                        Expr::Variable(v) => v,
                        other => format!("{}", other),
                    };
                    let lower = extractor.extract(egraph, args[3]);
                    let upper = extractor.extract(egraph, args[4]);
                    let space = crate::symbolic::functional_analysis::HilbertSpace {
                        var,
                        lower_bound: lower,
                        upper_bound: upper,
                    };
                    let res = crate::symbolic::functional_analysis::project_internal(&space, &f, &g);
                    results.push((class_id, res));
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Lie Algebra Oracle Rule (Tier 0).
/// Bridges Lie bracket commutator [X, Y] = XY - YX.
#[derive(Debug)]
pub struct LieAlgebraOracleRule;

impl Rule for LieAlgebraOracleRule {
    fn name(&self) -> &str {
        "oracles::lie_algebra"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::OracleCall(name, args) = node {
                    if name == "lie_bracket" && args.len() >= 2 {
                        targets.push((class_id, args[0], args[1]));
                    }
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, x_id, y_id) in targets {
            let x = extractor.extract(egraph, x_id);
            let y = extractor.extract(egraph, y_id);
            if let Ok(res) = crate::symbolic::lie_groups_and_algebras::lie_bracket_internal(&x, &y) {
                results.push((class_id, res));
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Integral Equations Oracle Rule (Tier 0).
/// Solves singular integral equations such as the airfoil equation.
#[derive(Debug)]
pub struct IntegralEquationsOracleRule;

impl Rule for IntegralEquationsOracleRule {
    fn name(&self) -> &str {
        "oracles::integral_equations"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::OracleCall(name, args) = node {
                    if name == "solve_airfoil_equation" && args.len() >= 3 {
                        targets.push((class_id, args[0], args[1], args[2]));
                    }
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, f_id, vx_id, vt_id) in targets {
            let f_x = extractor.extract(egraph, f_id);
            let var_x = match extractor.extract(egraph, vx_id) {
                Expr::Variable(v) => v,
                other => format!("{}", other),
            };
            let var_t = match extractor.extract(egraph, vt_id) {
                Expr::Variable(v) => v,
                other => format!("{}", other),
            };
            let res = crate::symbolic::integral_equations::solve_airfoil_equation_internal(&f_x, &var_x, &var_t);
            results.push((class_id, res));
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Statistics and Information Theory Oracle Rule (Tier 0).
/// Bridges mean, variance, std_dev, covariance, correlation, shannon_entropy, and gini_impurity.
#[derive(Debug)]
pub struct StatsOracleRule;

impl Rule for StatsOracleRule {
    fn name(&self) -> &str {
        "oracles::stats"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut targets = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::OracleCall(name, args) = node {
                    match name.as_str() {
                        "stats_mean"
                        | "stats_variance"
                        | "stats_std_dev"
                        | "stats_covariance"
                        | "stats_correlation"
                        | "shannon_entropy"
                        | "gini_impurity" => {
                            targets.push((class_id, name.clone(), args.clone()));
                        }
                        _ => {}
                    }
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let extractor = Extractor::new(egraph);
        let mut results = Vec::new();

        for (class_id, name, args) in targets {
            let extracted_args: Vec<Expr> = args.iter().map(|&id| extractor.extract(egraph, id)).collect();
            match name.as_str() {
                "stats_mean" => {
                    let res = crate::symbolic::stats::mean_internal(&extracted_args);
                    results.push((class_id, res));
                }
                "stats_variance" => {
                    let res = crate::symbolic::stats::variance_internal(&extracted_args);
                    results.push((class_id, res));
                }
                "stats_std_dev" => {
                    let res = crate::symbolic::stats::std_dev_internal(&extracted_args);
                    results.push((class_id, res));
                }
                "stats_covariance" => {
                    if !extracted_args.is_empty() {
                        let len = extracted_args[0].to_f64().unwrap_or(0.0) as usize;
                        if extracted_args.len() >= 1 + 2 * len {
                            let data1 = &extracted_args[1..1 + len];
                            let data2 = &extracted_args[1 + len..1 + 2 * len];
                            let res = crate::symbolic::stats::covariance_internal(data1, data2);
                            results.push((class_id, res));
                        }
                    }
                }
                "stats_correlation" => {
                    if !extracted_args.is_empty() {
                        let len = extracted_args[0].to_f64().unwrap_or(0.0) as usize;
                        if extracted_args.len() >= 1 + 2 * len {
                            let data1 = &extracted_args[1..1 + len];
                            let data2 = &extracted_args[1 + len..1 + 2 * len];
                            let res = crate::symbolic::stats::correlation_internal(data1, data2);
                            results.push((class_id, res));
                        }
                    }
                }
                "shannon_entropy" => {
                    let res = crate::symbolic::stats_information_theory::shannon_entropy_internal(&extracted_args);
                    results.push((class_id, res));
                }
                "gini_impurity" => {
                    let res = crate::symbolic::stats_information_theory::gini_impurity_internal(&extracted_args);
                    results.push((class_id, res));
                }
                _ => {}
            }
        }

        let mut applied = 0;
        for (class_id, res_expr) in results {
            let res_id = egraph.add_expr(&res_expr);
            if egraph.union(class_id, res_id) {
                applied += 1;
            }
        }
        applied
    }
}


