use std::collections::HashMap;
use std::sync::Arc;

use super::egraph::EGraph;
use super::enode::ENode;
use super::id::Id;
use crate::symbolic::core::Expr;

/// Extracts the minimum-cost canonical `Expr` from an `EGraph`, preserving DAG structural
/// sharing via `Arc<Expr>` and guaranteeing cycle-freedom.
#[derive(Clone, Debug)]
pub struct Extractor {
    costs: HashMap<Id, (u64, ENode)>,
}

impl Extractor {
    /// Computes the optimal (lowest-cost) representation for every reachable E-Class
    /// in the E-Graph using a fixpoint relaxation algorithm.
    pub fn new(egraph: &EGraph) -> Self {
        let mut costs: HashMap<Id, (u64, ENode)> = HashMap::new();
        let mut changed = true;

        // Iterate until costs converge to a fixpoint
        while changed {
            changed = false;

            for (idx, opt_class) in egraph.classes.iter().enumerate() {
                let class = match opt_class {
                    Some(c) => c,
                    None => continue,
                };
                let class_id = egraph.find_immut(Id::from_usize(idx));

                for node in &class.nodes {
                    let mut node_cost = node.complexity_weight() as u64;

                    // Heuristic penalty for unnormalized negation (Not wrapping compound expressions or quantifiers)
                    if let ENode::Not(child) = node {
                        let child_root = egraph.find_immut(*child);
                        if let Some(c) = egraph.get_class(child_root) {
                            let is_compound = c.nodes.iter().any(|n| {
                                matches!(
                                    n,
                                    ENode::ForAll(..)
                                        | ENode::Exists(..)
                                        | ENode::And(..)
                                        | ENode::Or(..)
                                        | ENode::Not(..)
                                        | ENode::Implies(..)
                                        | ENode::Equivalent(..)
                                        | ENode::Xor(..)
                                )
                            });
                            if is_compound {
                                node_cost = node_cost.saturating_add(50);
                            }
                        }
                    }

                    let mut all_children_have_cost = true;

                    for child_id in node.children() {
                        let child_root = egraph.find_immut(child_id);
                        if let Some(&(child_cost, _)) = costs.get(&child_root) {
                            node_cost = node_cost.saturating_add(child_cost);
                        } else {
                            all_children_have_cost = false;
                            break;
                        }
                    }

                    if all_children_have_cost {
                        let is_better = match costs.get(&class_id) {
                            Some(&(existing_cost, _)) => node_cost < existing_cost,
                            None => true,
                        };

                        if is_better {
                            costs.insert(class_id, (node_cost, node.clone()));
                            changed = true;
                        }
                    }
                }
            }
        }

        Self { costs }
    }

    /// Extracts the optimal `Expr` for the given class `root`, caching subexpressions
    /// in `memo` to guarantee maximal DAG structural sharing.
    pub fn extract(&self, egraph: &EGraph, root: Id) -> Expr {
        let mut memo: HashMap<Id, Arc<Expr>> = HashMap::new();
        let canonical_root = egraph.find_immut(root);
        (*self.extract_shared(egraph, canonical_root, &mut memo)).clone()
    }

    /// Recursively extracts an `Arc<Expr>`, sharing identical subexpression pointers in $O(1)$.
    fn extract_shared(
        &self,
        egraph: &EGraph,
        class_id: Id,
        memo: &mut HashMap<Id, Arc<Expr>>,
    ) -> Arc<Expr> {
        let canonical_id = egraph.find_immut(class_id);

        if let Some(existing) = memo.get(&canonical_id) {
            return existing.clone();
        }

        let best_node = match self.costs.get(&canonical_id) {
            Some((_, node)) => node.clone(),
            None => {
                // Fallback: take the first node in class if disconnected
                if let Some(class) = egraph.get_class(canonical_id) {
                    class.nodes.first().cloned().unwrap_or(ENode::Variable("undefined".to_string()))
                } else {
                    ENode::Variable("undefined".to_string())
                }
            }
        };

        let expr = match best_node {
            ENode::Constant(c) => Expr::Constant(c.0),
            ENode::BigInt(b) => Expr::BigInt(b),
            ENode::Rational(r) => Expr::Rational(r),
            ENode::Boolean(b) => Expr::Boolean(b),
            ENode::Variable(name) => Expr::Variable(name),
            ENode::Pattern(name) => Expr::Pattern(name),
            ENode::Pi => Expr::Pi,
            ENode::E => Expr::E,
            ENode::Infinity => Expr::Infinity,
            ENode::NegativeInfinity => Expr::NegativeInfinity,

            ENode::Add(a, b) => {
                let ea = self.extract_shared(egraph, a, memo);
                let eb = self.extract_shared(egraph, b, memo);
                Expr::Add(ea, eb)
            }
            ENode::Sub(a, b) => {
                let ea = self.extract_shared(egraph, a, memo);
                let eb = self.extract_shared(egraph, b, memo);
                Expr::Sub(ea, eb)
            }
            ENode::Mul(a, b) => {
                let ea = self.extract_shared(egraph, a, memo);
                let eb = self.extract_shared(egraph, b, memo);
                Expr::Mul(ea, eb)
            }
            ENode::Div(a, b) => {
                let ea = self.extract_shared(egraph, a, memo);
                let eb = self.extract_shared(egraph, b, memo);
                Expr::Div(ea, eb)
            }
            ENode::Power(a, b) => {
                let ea = self.extract_shared(egraph, a, memo);
                let eb = self.extract_shared(egraph, b, memo);
                Expr::Power(ea, eb)
            }
            ENode::Neg(a) => {
                let ea = self.extract_shared(egraph, a, memo);
                Expr::Neg(ea)
            }
            ENode::Complex(re, im) => {
                let ere = self.extract_shared(egraph, re, memo);
                let eim = self.extract_shared(egraph, im, memo);
                Expr::Complex(ere, eim)
            }

            ENode::AddList(list) => {
                let exprs = list.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::AddList(exprs)
            }
            ENode::MulList(list) => {
                let exprs = list.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::MulList(exprs)
            }

            ENode::Sin(a) => Expr::Sin(self.extract_shared(egraph, a, memo)),
            ENode::Cos(a) => Expr::Cos(self.extract_shared(egraph, a, memo)),
            ENode::Tan(a) => Expr::Tan(self.extract_shared(egraph, a, memo)),
            ENode::Sec(a) => Expr::Sec(self.extract_shared(egraph, a, memo)),
            ENode::Csc(a) => Expr::Csc(self.extract_shared(egraph, a, memo)),
            ENode::Cot(a) => Expr::Cot(self.extract_shared(egraph, a, memo)),
            ENode::ArcSin(a) => Expr::ArcSin(self.extract_shared(egraph, a, memo)),
            ENode::ArcCos(a) => Expr::ArcCos(self.extract_shared(egraph, a, memo)),
            ENode::ArcTan(a) => Expr::ArcTan(self.extract_shared(egraph, a, memo)),
            ENode::Sinh(a) => Expr::Sinh(self.extract_shared(egraph, a, memo)),
            ENode::Cosh(a) => Expr::Cosh(self.extract_shared(egraph, a, memo)),
            ENode::Tanh(a) => Expr::Tanh(self.extract_shared(egraph, a, memo)),
            ENode::Exp(a) => Expr::Exp(self.extract_shared(egraph, a, memo)),
            ENode::Log(a) => Expr::Log(self.extract_shared(egraph, a, memo)),
            ENode::Abs(a) => Expr::Abs(self.extract_shared(egraph, a, memo)),
            ENode::Sqrt(a) => Expr::Sqrt(self.extract_shared(egraph, a, memo)),

            ENode::Derivative(body, var) => {
                let ebody = self.extract_shared(egraph, body, memo);
                Expr::Derivative(ebody, var)
            }
            ENode::DerivativeN(body, var, n) => {
                let ebody = self.extract_shared(egraph, body, memo);
                let en = self.extract_shared(egraph, n, memo);
                Expr::DerivativeN(ebody, var, en)
            }
            ENode::Integral { integrand, var, lower_bound, upper_bound } => {
                let ei = self.extract_shared(egraph, integrand, memo);
                let ev = self.extract_shared(egraph, var, memo);
                let elb = self.extract_shared(egraph, lower_bound, memo);
                let eub = self.extract_shared(egraph, upper_bound, memo);
                Expr::Integral {
                    integrand: ei,
                    var: ev,
                    lower_bound: elb,
                    upper_bound: eub,
                }
            }
            ENode::Limit(body, var, target) => {
                let eb = self.extract_shared(egraph, body, memo);
                let et = self.extract_shared(egraph, target, memo);
                Expr::Limit(eb, var, et)
            }
            ENode::Sum { body, var, from, to } => {
                let eb = self.extract_shared(egraph, body, memo);
                let ev = self.extract_shared(egraph, var, memo);
                let ef = self.extract_shared(egraph, from, memo);
                let et = self.extract_shared(egraph, to, memo);
                Expr::Sum {
                    body: eb,
                    var: ev,
                    from: ef,
                    to: et,
                }
            }
            ENode::Product(body, var, from, to) => {
                let eb = self.extract_shared(egraph, body, memo);
                let ef = self.extract_shared(egraph, from, memo);
                let et = self.extract_shared(egraph, to, memo);
                Expr::Product(eb, var, ef, et)
            }
            ENode::Series(body, var, pt, order) => {
                let eb = self.extract_shared(egraph, body, memo);
                let ep = self.extract_shared(egraph, pt, memo);
                let eo = self.extract_shared(egraph, order, memo);
                Expr::Series(eb, var, ep, eo)
            }
            ENode::Solve(eq, var) => {
                let ee = self.extract_shared(egraph, eq, memo);
                Expr::Solve(ee, var)
            }
            ENode::RootOf(poly, index) => {
                let ep = self.extract_shared(egraph, poly, memo);
                Expr::RootOf { poly: ep, index }
            }
            ENode::Ode { equation, func, var } => {
                let ee = self.extract_shared(egraph, equation, memo);
                Expr::Ode {
                    equation: ee,
                    func,
                    var,
                }
            }
            ENode::Pde { equation, func, vars } => {
                let ee = self.extract_shared(egraph, equation, memo);
                Expr::Pde {
                    equation: ee,
                    func,
                    vars,
                }
            }

            ENode::Matrix(rows) => {
                let row_exprs = rows.iter().map(|row| {
                    row.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect()
                }).collect();
                Expr::Matrix(row_exprs)
            }
            ENode::Vector(v) => {
                let exprs = v.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::Vector(exprs)
            }
            ENode::Transpose(a) => Expr::Transpose(self.extract_shared(egraph, a, memo)),
            ENode::MatrixMul(a, b) => Expr::MatrixMul(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Inverse(a) => Expr::Inverse(self.extract_shared(egraph, a, memo)),

            ENode::Eq(a, b) => Expr::Eq(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Lt(a, b) => Expr::Lt(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Gt(a, b) => Expr::Gt(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Le(a, b) => Expr::Le(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Ge(a, b) => Expr::Ge(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),

            ENode::Mod(a, b) => Expr::Mod(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Floor(a) => Expr::Floor(self.extract_shared(egraph, a, memo)),
            ENode::Factorial(a) => Expr::Factorial(self.extract_shared(egraph, a, memo)),
            ENode::Gcd(a, b) => Expr::Gcd(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Max(a, b) => Expr::Max(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Expand(a) => Expr::UnaryList("expand".to_string(), self.extract_shared(egraph, a, memo)),
            ENode::IsPrime(a) => Expr::IsPrime(self.extract_shared(egraph, a, memo)),
            ENode::IndefiniteSum { body, var, step } => {
                let eb = self.extract_shared(egraph, body, memo);
                let es = self.extract_shared(egraph, step, memo);
                Expr::IndefiniteSum { body: eb, var, step: es }
            }
            ENode::IndefiniteProduct { body, var, step } => {
                let eb = self.extract_shared(egraph, body, memo);
                let es = self.extract_shared(egraph, step, memo);
                Expr::IndefiniteProduct { body: eb, var, step: es }
            }

            ENode::ArcSec(a) => Expr::ArcSec(self.extract_shared(egraph, a, memo)),
            ENode::ArcCsc(a) => Expr::ArcCsc(self.extract_shared(egraph, a, memo)),
            ENode::ArcCot(a) => Expr::ArcCot(self.extract_shared(egraph, a, memo)),
            ENode::Sech(a) => Expr::Sech(self.extract_shared(egraph, a, memo)),
            ENode::Csch(a) => Expr::Csch(self.extract_shared(egraph, a, memo)),
            ENode::Coth(a) => Expr::Coth(self.extract_shared(egraph, a, memo)),
            ENode::ArcSinh(a) => Expr::ArcSinh(self.extract_shared(egraph, a, memo)),
            ENode::ArcCosh(a) => Expr::ArcCosh(self.extract_shared(egraph, a, memo)),
            ENode::ArcTanh(a) => Expr::ArcTanh(self.extract_shared(egraph, a, memo)),
            ENode::ArcSech(a) => Expr::ArcSech(self.extract_shared(egraph, a, memo)),
            ENode::ArcCsch(a) => Expr::ArcCsch(self.extract_shared(egraph, a, memo)),
            ENode::ArcCoth(a) => Expr::ArcCoth(self.extract_shared(egraph, a, memo)),
            ENode::LogBase(a, b) => Expr::LogBase(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Atan2(a, b) => Expr::Atan2(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),

            ENode::Gamma(a) => Expr::Gamma(self.extract_shared(egraph, a, memo)),
            ENode::Beta(a, b) => Expr::Beta(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Erf(a) => Expr::Erf(self.extract_shared(egraph, a, memo)),
            ENode::Erfc(a) => Expr::Erfc(self.extract_shared(egraph, a, memo)),
            ENode::Erfi(a) => Expr::Erfi(self.extract_shared(egraph, a, memo)),
            ENode::Zeta(a) => Expr::Zeta(self.extract_shared(egraph, a, memo)),
            ENode::Digamma(a) => Expr::Digamma(self.extract_shared(egraph, a, memo)),
            ENode::BesselJ(a, b) => Expr::BesselJ(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::BesselY(a, b) => Expr::BesselY(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::LegendreP(a, b) => Expr::LegendreP(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::LaguerreL(a, b) => Expr::LaguerreL(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::HermiteH(a, b) => Expr::HermiteH(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),

            ENode::And(v) => {
                let exprs = v.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::And(exprs)
            }
            ENode::Or(v) => {
                let exprs = v.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::Or(exprs)
            }
            ENode::Not(a) => Expr::Not(self.extract_shared(egraph, a, memo)),
            ENode::Xor(a, b) => Expr::Xor(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Implies(a, b) => Expr::Implies(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Equivalent(a, b) => Expr::Equivalent(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Predicate { name, args } => {
                let exprs = args.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::Predicate { name, args: exprs }
            }
            ENode::ForAll(var, body) => Expr::ForAll(var, self.extract_shared(egraph, body, memo)),
            ENode::Exists(var, body) => Expr::Exists(var, self.extract_shared(egraph, body, memo)),
            ENode::Apply(a, b) => Expr::Apply(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Tuple(v) => {
                let exprs = v.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::Tuple(exprs)
            }
            ENode::Solutions(v) => {
                let exprs = v.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                Expr::Solutions(exprs)
            }
            ENode::InfiniteSolutions => Expr::InfiniteSolutions,
            ENode::NoSolution => Expr::NoSolution,

            // Factor and SimplifyWithRelations are E-Graph-internal oracle ops.
            // Extract them as their inner expression (the oracle will have already
            // resolved them to concrete results in the E-Class).
            ENode::Factor(a) => (*self.extract_shared(egraph, a, memo)).clone(),
            ENode::SimplifyWithRelations { expr, .. } => {
                (*self.extract_shared(egraph, expr, memo)).clone()
            }

            ENode::Determinant(a) => Expr::UnaryList("det".to_string(), self.extract_shared(egraph, a, memo)),
            ENode::Trace(a) => Expr::UnaryList("trace".to_string(), self.extract_shared(egraph, a, memo)),
            ENode::Laplace { expr, var, target } => Expr::NaryList(
                "laplace".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::InverseLaplace { expr, var, target } => Expr::NaryList(
                "inv_laplace".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::Fourier { expr, var, target } => Expr::NaryList(
                "fourier".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::InverseFourier { expr, var, target } => Expr::NaryList(
                "inv_fourier".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::ZTransform { expr, var, target } => Expr::NaryList(
                "z_transform".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::InverseZTransform { expr, var, target } => Expr::NaryList(
                "inv_z_transform".to_string(),
                vec![(*self.extract_shared(egraph, expr, memo)).clone(), Expr::Variable(var), Expr::Variable(target)],
            ),
            ENode::Gradient { expr, vars } => Expr::NaryList(
                "gradient".to_string(),
                std::iter::once((*self.extract_shared(egraph, expr, memo)).clone())
                    .chain(vars.into_iter().map(Expr::Variable))
                    .collect(),
            ),
            ENode::Divergence { expr, vars } => Expr::NaryList(
                "divergence".to_string(),
                std::iter::once((*self.extract_shared(egraph, expr, memo)).clone())
                    .chain(vars.into_iter().map(Expr::Variable))
                    .collect(),
            ),
            ENode::Curl { expr, vars } => Expr::NaryList(
                "curl".to_string(),
                std::iter::once((*self.extract_shared(egraph, expr, memo)).clone())
                    .chain(vars.into_iter().map(Expr::Variable))
                    .collect(),
            ),
            ENode::Laplacian { expr, vars } => Expr::NaryList(
                "laplacian".to_string(),
                std::iter::once((*self.extract_shared(egraph, expr, memo)).clone())
                    .chain(vars.into_iter().map(Expr::Variable))
                    .collect(),
            ),
            ENode::EulerLagrange { lagrangian, func, var } => Expr::NaryList(
                "euler_lagrange".to_string(),
                vec![(*self.extract_shared(egraph, lagrangian, memo)).clone(), Expr::Variable(func), Expr::Variable(var)],
            ),
            ENode::Commutator(a, b) => Expr::BinaryList(
                "commutator".to_string(),
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Anticommutator(a, b) => Expr::BinaryList(
                "anticommutator".to_string(),
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Residue { expr, var, point } => Expr::NaryList(
                "residue".to_string(),
                vec![
                    (*self.extract_shared(egraph, expr, memo)).clone(),
                    Expr::Variable(var),
                    (*self.extract_shared(egraph, point, memo)).clone(),
                ],
            ),
            ENode::Binomial(a, b) => Expr::Binomial(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Permutation(a, b) => Expr::Permutation(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Combination(a, b) => Expr::Combination(
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::Lcm(a, b) => Expr::BinaryList(
                "lcm".to_string(),
                self.extract_shared(egraph, a, memo),
                self.extract_shared(egraph, b, memo),
            ),
            ENode::OracleCall(name, args) => {
                let exprs: Vec<Expr> = args.iter().map(|id| (*self.extract_shared(egraph, *id, memo)).clone()).collect();
                if exprs.len() == 1 {
                    Expr::UnaryList(name, Arc::new(exprs.into_iter().next().unwrap()))
                } else if exprs.len() == 2 {
                    let mut it = exprs.into_iter();
                    let a = it.next().unwrap();
                    let b = it.next().unwrap();
                    Expr::BinaryList(name, Arc::new(a), Arc::new(b))
                } else {
                    Expr::NaryList(name, exprs)
                }
            }
        };

        let shared = Arc::new(expr);
        memo.insert(canonical_id, shared.clone());
        shared
    }
}
