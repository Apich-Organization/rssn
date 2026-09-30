use std::collections::HashMap;
use ordered_float::OrderedFloat;

use super::eclass::EClass;
use super::enode::ENode;
use super::id::Id;
use super::union_find::UnionFind;
use crate::symbolic::core::Expr;

/// The central Equivalence Graph (E-Graph) maintaining equivalence classes of expressions.
#[derive(Clone, Debug)]
pub struct EGraph {
    /// Disjoint-set union-find for maintaining equivalence relations.
    pub union_find: UnionFind,
    /// Table of equivalence classes, indexed by canonical `Id`.
    pub classes: Vec<Option<EClass>>,
    /// Hashconsing memo table mapping canonical `ENode`s to their containing `Id`.
    pub memo: HashMap<ENode, Id>,
    /// Queue of pending union operations to be processed during rebuild.
    pub pending_unions: Vec<(Id, Id)>,
}

impl Default for EGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl EGraph {
    /// Creates a new, empty `EGraph`.
    #[must_use]
    pub fn new() -> Self {
        Self {
            union_find: UnionFind::new(),
            classes: Vec::new(),
            memo: HashMap::new(),
            pending_unions: Vec::new(),
        }
    }

    /// Returns the canonical representative `Id` for `id`.
    #[inline(always)]
    pub fn find(&mut self, id: Id) -> Id {
        self.union_find.find(id)
    }

    /// Returns the canonical representative `Id` for `id` immutably.
    #[inline(always)]
    #[must_use]
    pub fn find_immut(&self, id: Id) -> Id {
        self.union_find.find_immut(id)
    }

    /// Retrieves an immutable reference to the canonical `EClass` for `id`.
    #[must_use]
    pub fn get_class(&self, id: Id) -> Option<&EClass> {
        let root = self.find_immut(id);
        self.classes.get(root.as_usize()).and_then(Option::as_ref)
    }

    /// Retrieves a mutable reference to the canonical `EClass` for `id`.
    pub fn get_class_mut(&mut self, id: Id) -> Option<&mut EClass> {
        let root = self.find(id);
        self.classes.get_mut(root.as_usize()).and_then(Option::as_mut)
    }

    /// Adds an `Expr` (AST / shared DAG) into the E-Graph and returns its canonical `Id`.
    pub fn add_expr(&mut self, expr: &Expr) -> Id {
        match expr {
            Expr::Constant(c) => self.add_node(ENode::Constant(OrderedFloat(*c))),
            Expr::BigInt(b) => self.add_node(ENode::BigInt(b.clone())),
            Expr::Rational(r) => self.add_node(ENode::Rational(r.clone())),
            Expr::Boolean(b) => self.add_node(ENode::Boolean(*b)),
            Expr::Variable(name) => self.add_node(ENode::Variable(name.clone())),
            Expr::Pattern(name) => self.add_node(ENode::Pattern(name.clone())),
            Expr::Pi => self.add_node(ENode::Pi),
            Expr::E => self.add_node(ENode::E),
            Expr::Infinity => self.add_node(ENode::Infinity),
            Expr::NegativeInfinity => self.add_node(ENode::NegativeInfinity),

            Expr::Add(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Add(id_a, id_b))
            }
            Expr::Sub(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Sub(id_a, id_b))
            }
            Expr::Mul(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Mul(id_a, id_b))
            }
            Expr::Div(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Div(id_a, id_b))
            }
            Expr::Power(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Power(id_a, id_b))
            }
            Expr::Neg(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Neg(id_a))
            }
            Expr::Complex(re, im) => {
                let id_re = self.add_expr(re);
                let id_im = self.add_expr(im);
                self.add_node(ENode::Complex(id_re, id_im))
            }

            Expr::AddList(list) => {
                let ids: Vec<Id> = list.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::AddList(ids))
            }
            Expr::MulList(list) => {
                let ids: Vec<Id> = list.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::MulList(ids))
            }

            Expr::Sin(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Sin(id_a))
            }
            Expr::Cos(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Cos(id_a))
            }
            Expr::Tan(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Tan(id_a))
            }
            Expr::Exp(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Exp(id_a))
            }
            Expr::Log(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Log(id_a))
            }
            Expr::Abs(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Abs(id_a))
            }
            Expr::Sqrt(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Sqrt(id_a))
            }

            // High-order operations
            Expr::Derivative(body, var) => {
                let id_body = self.add_expr(body);
                self.add_node(ENode::Derivative(id_body, var.clone()))
            }
            Expr::DerivativeN(body, var, n) => {
                let id_body = self.add_expr(body);
                let id_n = self.add_expr(n);
                self.add_node(ENode::DerivativeN(id_body, var.clone(), id_n))
            }
            Expr::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => {
                let id_int = self.add_expr(integrand);
                let id_var = self.add_expr(var);
                let id_lb = self.add_expr(lower_bound);
                let id_ub = self.add_expr(upper_bound);
                self.add_node(ENode::Integral {
                    integrand: id_int,
                    var: id_var,
                    lower_bound: id_lb,
                    upper_bound: id_ub,
                })
            }
            Expr::Limit(body, var, target) => {
                let id_body = self.add_expr(body);
                let id_target = self.add_expr(target);
                self.add_node(ENode::Limit(id_body, var.clone(), id_target))
            }
            Expr::Sum { body, var, from, to } => {
                let id_body = self.add_expr(body);
                let id_var = self.add_expr(var);
                let id_from = self.add_expr(from);
                let id_to = self.add_expr(to);
                self.add_node(ENode::Sum {
                    body: id_body,
                    var: id_var,
                    from: id_from,
                    to: id_to,
                })
            }
            Expr::Product(body, var, from, to) => {
                let id_body = self.add_expr(body);
                let id_from = self.add_expr(from);
                let id_to = self.add_expr(to);
                self.add_node(ENode::Product(id_body, var.clone(), id_from, id_to))
            }
            Expr::Series(body, var, pt, order) => {
                let id_body = self.add_expr(body);
                let id_pt = self.add_expr(pt);
                let id_order = self.add_expr(order);
                self.add_node(ENode::Series(id_body, var.clone(), id_pt, id_order))
            }
            Expr::Solve(eq, var) => {
                let id_eq = self.add_expr(eq);
                self.add_node(ENode::Solve(id_eq, var.clone()))
            }
            Expr::RootOf { poly, index } => {
                let id_poly = self.add_expr(poly);
                self.add_node(ENode::RootOf(id_poly, *index))
            }
            Expr::Ode { equation, func, var } => {
                let id_eq = self.add_expr(equation);
                self.add_node(ENode::Ode {
                    equation: id_eq,
                    func: func.clone(),
                    var: var.clone(),
                })
            }
            Expr::Pde { equation, func, vars } => {
                let id_eq = self.add_expr(equation);
                self.add_node(ENode::Pde {
                    equation: id_eq,
                    func: func.clone(),
                    vars: vars.clone(),
                })
            }

            Expr::Sec(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Sec(id))
            }
            Expr::Csc(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Csc(id))
            }
            Expr::Cot(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Cot(id))
            }
            Expr::ArcSin(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::ArcSin(id))
            }
            Expr::ArcCos(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::ArcCos(id))
            }
            Expr::ArcTan(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::ArcTan(id))
            }
            Expr::Sinh(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Sinh(id))
            }
            Expr::Cosh(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Cosh(id))
            }
            Expr::Tanh(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Tanh(id))
            }

            Expr::Eq(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Eq(id_a, id_b))
            }
            Expr::Lt(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Lt(id_a, id_b))
            }
            Expr::Gt(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Gt(id_a, id_b))
            }
            Expr::Le(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Le(id_a, id_b))
            }
            Expr::Ge(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Ge(id_a, id_b))
            }

            Expr::Matrix(rows) => {
                let row_ids: Vec<Vec<Id>> = rows
                    .iter()
                    .map(|r| r.iter().map(|e| self.add_expr(e)).collect())
                    .collect();
                self.add_node(ENode::Matrix(row_ids))
            }
            Expr::Vector(v) => {
                let ids = v.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::Vector(ids))
            }
            Expr::Transpose(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Transpose(id))
            }
            Expr::MatrixMul(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::MatrixMul(id_a, id_b))
            }
            Expr::Inverse(a) => {
                let id = self.add_expr(a);
                self.add_node(ENode::Inverse(id))
            }

            // Fallback for legacy Dag wrapper if encountered during migration
            #[allow(deprecated)]
            Expr::Dag(dag_node) => {
                if let Ok(ast) = dag_node.to_expr() {
                    self.add_expr(&ast)
                } else {
                    self.add_node(ENode::Variable("err".to_string()))
                }
            }

            // Number theory / arithmetic
            Expr::Mod(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Mod(id_a, id_b))
            }
            Expr::Floor(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Floor(id_a))
            }
            Expr::IsPrime(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::IsPrime(id_a))
            }
            Expr::Factorial(a) => {
                let id_a = self.add_expr(a);
                self.add_node(ENode::Factorial(id_a))
            }
            Expr::Gcd(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Gcd(id_a, id_b))
            }
            Expr::Max(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Max(id_a, id_b))
            }
            Expr::Atan2(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Atan2(id_a, id_b))
            }

            // Indefinite summation / product
            Expr::IndefiniteSum { body, var, step } => {
                let id_body = self.add_expr(body);
                let id_step = self.add_expr(step);
                self.add_node(ENode::IndefiniteSum { body: id_body, var: var.clone(), step: id_step })
            }
            Expr::IndefiniteProduct { body, var, step } => {
                let id_body = self.add_expr(body);
                let id_step = self.add_expr(step);
                self.add_node(ENode::IndefiniteProduct { body: id_body, var: var.clone(), step: id_step })
            }

            // Additional trig / hyp not yet mapped
            Expr::Sech(a) => { let id = self.add_expr(a); self.add_node(ENode::Sech(id)) }
            Expr::Csch(a) => { let id = self.add_expr(a); self.add_node(ENode::Csch(id)) }
            Expr::Coth(a) => { let id = self.add_expr(a); self.add_node(ENode::Coth(id)) }
            Expr::ArcSinh(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcSinh(id)) }
            Expr::ArcCosh(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcCosh(id)) }
            Expr::ArcTanh(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcTanh(id)) }
            Expr::ArcSec(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcSec(id)) }
            Expr::ArcCsc(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcCsc(id)) }
            Expr::ArcCot(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcCot(id)) }
            Expr::ArcSech(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcSech(id)) }
            Expr::ArcCsch(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcCsch(id)) }
            Expr::ArcCoth(a) => { let id = self.add_expr(a); self.add_node(ENode::ArcCoth(id)) }
            Expr::LogBase(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::LogBase(id_a, id_b))
            }

            // Special functions
            Expr::Gamma(a) => { let id = self.add_expr(a); self.add_node(ENode::Gamma(id)) }
            Expr::Beta(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::Beta(id_a, id_b))
            }
            Expr::Erf(a) => { let id = self.add_expr(a); self.add_node(ENode::Erf(id)) }
            Expr::Erfc(a) => { let id = self.add_expr(a); self.add_node(ENode::Erfc(id)) }
            Expr::Erfi(a) => { let id = self.add_expr(a); self.add_node(ENode::Erfi(id)) }
            Expr::Zeta(a) => { let id = self.add_expr(a); self.add_node(ENode::Zeta(id)) }
            Expr::Digamma(a) => { let id = self.add_expr(a); self.add_node(ENode::Digamma(id)) }
            Expr::BesselJ(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::BesselJ(id_a, id_b))
            }
            Expr::BesselY(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::BesselY(id_a, id_b))
            }
            Expr::LegendreP(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::LegendreP(id_a, id_b))
            }
            Expr::LaguerreL(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::LaguerreL(id_a, id_b))
            }
            Expr::HermiteH(a, b) => {
                let id_a = self.add_expr(a); let id_b = self.add_expr(b);
                self.add_node(ENode::HermiteH(id_a, id_b))
            }

            // Logic
            Expr::And(items) => {
                let ids: Vec<Id> = items.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::And(ids))
            }
            Expr::Or(items) => {
                let ids: Vec<Id> = items.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::Or(ids))
            }
            Expr::Not(a) => { let id = self.add_expr(a); self.add_node(ENode::Not(id)) }
            Expr::Xor(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Xor(id_a, id_b))
            }
            Expr::Implies(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Implies(id_a, id_b))
            }
            Expr::Equivalent(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Equivalent(id_a, id_b))
            }
            Expr::Predicate { name, args } => {
                let ids: Vec<Id> = args.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::Predicate { name: name.clone(), args: ids })
            }
            Expr::ForAll(var, body) => {
                let id_body = self.add_expr(body);
                self.add_node(ENode::ForAll(var.clone(), id_body))
            }
            Expr::Exists(var, body) => {
                let id_body = self.add_expr(body);
                self.add_node(ENode::Exists(var.clone(), id_body))
            }

            // Solutions / Apply / Tuple
            Expr::Apply(a, b) => {
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::Apply(id_a, id_b))
            }
            Expr::Tuple(items) => {
                let ids: Vec<Id> = items.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::Tuple(ids))
            }
            Expr::Solutions(items) => {
                let ids: Vec<Id> = items.iter().map(|e| self.add_expr(e)).collect();
                self.add_node(ENode::Solutions(ids))
            }
            Expr::InfiniteSolutions => self.add_node(ENode::InfiniteSolutions),
            Expr::NoSolution => self.add_node(ENode::NoSolution),


            Expr::Summation(body, var, from, to) => {
                let id_body = self.add_expr(body);
                let id_var = self.add_node(ENode::Variable(var.clone()));
                let id_from = self.add_expr(from);
                let id_to = self.add_expr(to);
                self.add_node(ENode::Sum { body: id_body, var: id_var, from: id_from, to: id_to })
            }
            Expr::MatrixVecMul(a, b) => {
                // Treat as regular matrix multiplication for E-Graph purposes
                let id_a = self.add_expr(a);
                let id_b = self.add_expr(b);
                self.add_node(ENode::MatrixMul(id_a, id_b))
            }

            Expr::Binomial(n, k) => {
                let id_n = self.add_expr(n);
                let id_k = self.add_expr(k);
                self.add_node(ENode::Binomial(id_n, id_k))
            }
            Expr::Permutation(n, k) => {
                let id_n = self.add_expr(n);
                let id_k = self.add_expr(k);
                self.add_node(ENode::Permutation(id_n, id_k))
            }
            Expr::Combination(n, k) => {
                let id_n = self.add_expr(n);
                let id_k = self.add_expr(k);
                self.add_node(ENode::Combination(id_n, id_k))
            }
            Expr::UnaryList(op, inner) => match op.as_str() {
                "det" | "determinant" => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::Determinant(id))
                }
                "trace" => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::Trace(id))
                }
                "expand" => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::Expand(id))
                }
                "is_prime" => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::IsPrime(id))
                }
                "factor" => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::Factor(id))
                }
                _ => {
                    let id = self.add_expr(inner);
                    self.add_node(ENode::OracleCall(op.clone(), vec![id]))
                }
            },
            Expr::BinaryList(op, a, b) => match op.as_str() {
                "commutator" => {
                    let id_a = self.add_expr(a);
                    let id_b = self.add_expr(b);
                    self.add_node(ENode::Commutator(id_a, id_b))
                }
                "anticommutator" => {
                    let id_a = self.add_expr(a);
                    let id_b = self.add_expr(b);
                    self.add_node(ENode::Anticommutator(id_a, id_b))
                }
                "lcm" => {
                    let id_a = self.add_expr(a);
                    let id_b = self.add_expr(b);
                    self.add_node(ENode::Lcm(id_a, id_b))
                }
                "gcd" | "polynomial_gcd" => {
                    let id_a = self.add_expr(a);
                    let id_b = self.add_expr(b);
                    self.add_node(ENode::Gcd(id_a, id_b))
                }
                _ => {
                    let id_a = self.add_expr(a);
                    let id_b = self.add_expr(b);
                    self.add_node(ENode::OracleCall(op.clone(), vec![id_a, id_b]))
                }
            },
            Expr::NaryList(op, args) => match op.as_str() {
                "laplace" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::Laplace { expr: id_f, var, target })
                }
                "inv_laplace" | "inverse_laplace" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::InverseLaplace { expr: id_f, var, target })
                }
                "fourier" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::Fourier { expr: id_f, var, target })
                }
                "inv_fourier" | "inverse_fourier" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::InverseFourier { expr: id_f, var, target })
                }
                "z_transform" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::ZTransform { expr: id_f, var, target })
                }
                "inv_z_transform" | "inverse_z_transform" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let target = format!("{}", args[2]);
                    self.add_node(ENode::InverseZTransform { expr: id_f, var, target })
                }
                "gradient" if !args.is_empty() => {
                    let id_f = self.add_expr(&args[0]);
                    let vars = args[1..].iter().map(|e| format!("{}", e)).collect();
                    self.add_node(ENode::Gradient { expr: id_f, vars })
                }
                "divergence" if !args.is_empty() => {
                    let id_f = self.add_expr(&args[0]);
                    let vars = args[1..].iter().map(|e| format!("{}", e)).collect();
                    self.add_node(ENode::Divergence { expr: id_f, vars })
                }
                "curl" if !args.is_empty() => {
                    let id_f = self.add_expr(&args[0]);
                    let vars = args[1..].iter().map(|e| format!("{}", e)).collect();
                    self.add_node(ENode::Curl { expr: id_f, vars })
                }
                "laplacian" if !args.is_empty() => {
                    let id_f = self.add_expr(&args[0]);
                    let vars = args[1..].iter().map(|e| format!("{}", e)).collect();
                    self.add_node(ENode::Laplacian { expr: id_f, vars })
                }
                "euler_lagrange" if args.len() >= 3 => {
                    let id_lagrangian = self.add_expr(&args[0]);
                    let func = format!("{}", args[1]);
                    let var = format!("{}", args[2]);
                    self.add_node(ENode::EulerLagrange { lagrangian: id_lagrangian, func, var })
                }
                "residue" if args.len() >= 3 => {
                    let id_f = self.add_expr(&args[0]);
                    let var = format!("{}", args[1]);
                    let id_pt = self.add_expr(&args[2]);
                    self.add_node(ENode::Residue { expr: id_f, var, point: id_pt })
                }
                _ => {
                    let ids: Vec<Id> = args.iter().map(|e| self.add_expr(e)).collect();
                    self.add_node(ENode::OracleCall(op.clone(), ids))
                }
            },

            // Any remaining expressions map to Variable or default node
            _ => self.add_node(ENode::Variable(format!("{:?}", expr))),
        }
    }

    /// Adds a single canonicalized `ENode` into the E-Graph.
    /// If an equivalent node already exists, returns its canonical `Id`.
    /// Otherwise, creates a new `EClass` and registers parent edges.
    pub fn add_node(&mut self, mut node: ENode) -> Id {
        node.canonicalize(&mut self.union_find);

        if let Some(&id) = self.memo.get(&node) {
            return self.union_find.find(id);
        }

        let new_id = self.union_find.make_set();
        let mut class = EClass::new(new_id, node.clone());

        // Register parent pointers in children and propagate free_vars
        for child_id in node.children() {
            let child_root = self.union_find.find(child_id);
            if let Some(Some(child_class)) = self.classes.get_mut(child_root.as_usize()) {
                class.free_vars.extend(child_class.free_vars.iter().cloned());
                child_class.parents.push((node.clone(), new_id));
            }
        }

        // For Derivative/DerivativeN nodes, the variable string is metadata
        // (not a child Id), but the expression functionally depends on it.
        // Add it to free_vars so contains_var correctly detects dependency.
        match &node {
            ENode::Derivative(_, var)
            | ENode::DerivativeN(_, var, _)
            | ENode::Limit(_, var, _)
            | ENode::Solve(_, var)
            | ENode::Product(_, var, _, _)
            | ENode::Series(_, var, _, _) => {
                class.free_vars.insert(var.clone());
            }
            ENode::Ode { var, func, .. } => {
                class.free_vars.insert(var.clone());
                class.free_vars.insert(func.clone());
            }
            ENode::Pde { vars, func, .. } => {
                for v in vars {
                    class.free_vars.insert(v.clone());
                }
                class.free_vars.insert(func.clone());
            }
            ENode::Laplace { var, target, .. }
            | ENode::InverseLaplace { var, target, .. }
            | ENode::Fourier { var, target, .. }
            | ENode::InverseFourier { var, target, .. }
            | ENode::ZTransform { var, target, .. }
            | ENode::InverseZTransform { var, target, .. } => {
                class.free_vars.insert(var.clone());
                class.free_vars.insert(target.clone());
            }
            ENode::Gradient { vars, .. }
            | ENode::Divergence { vars, .. }
            | ENode::Curl { vars, .. }
            | ENode::Laplacian { vars, .. } => {
                for v in vars {
                    class.free_vars.insert(v.clone());
                }
            }
            ENode::EulerLagrange { func, var, .. } => {
                class.free_vars.insert(func.clone());
                class.free_vars.insert(var.clone());
            }
            ENode::Residue { var, .. } => {
                class.free_vars.insert(var.clone());
            }
            ENode::ForAll(var, _) | ENode::Exists(var, _) => {
                class.free_vars.remove(var);
            }
            _ => {}
        }

        // Store new class
        if new_id.as_usize() >= self.classes.len() {
            self.classes.resize_with(new_id.as_usize() + 1, || None);
        }
        self.classes[new_id.as_usize()] = Some(class);

        self.memo.insert(node, new_id);
        new_id
    }

    /// Merges the equivalence classes containing `id1` and `id2`.
    /// Returns `true` if this merged two previously distinct classes.
    pub fn union(&mut self, id1: Id, id2: Id) -> bool {
        let root1 = self.find(id1);
        let root2 = self.find(id2);
        if root1 == root2 {
            return false;
        }

        let (winner, changed) = self.union_find.union(root1, root2);
        if !changed {
            return false;
        }

        let loser = if winner == root1 { root2 } else { root1 };

        // Take loser class and merge into winner class
        if let Some(loser_class) = self.classes[loser.as_usize()].take() {
            if let Some(Some(winner_class)) = self.classes.get_mut(winner.as_usize()) {
                winner_class.merge(loser_class);
            }
        }

        self.pending_unions.push((winner, loser));
        true
    }

    /// Restores the Congruence Closure invariant across all E-Classes.
    /// Re-canonicalizes parent nodes and merges classes that have become structurally equivalent.
    /// Returns the number of class unions performed during rebuilding.
    pub fn rebuild(&mut self) -> usize {
        let mut total_unions = 0;

        while !self.pending_unions.is_empty() {
            let batch = std::mem::take(&mut self.pending_unions);

            for (winner, _) in batch {
                let canonical_winner = self.union_find.find(winner);
                let parents = if let Some(Some(class)) = self.classes.get_mut(canonical_winner.as_usize()) {
                    std::mem::take(&mut class.parents)
                } else {
                    continue;
                };

                let mut new_parents = Vec::with_capacity(parents.len());

                for (mut parent_node, parent_id) in parents {
                    self.memo.remove(&parent_node);
                    parent_node.canonicalize(&mut self.union_find);
                    let canonical_parent_id = self.union_find.find(parent_id);

                    if let Some(&existing_id) = self.memo.get(&parent_node) {
                        let canonical_existing = self.union_find.find(existing_id);
                        if canonical_existing != canonical_parent_id {
                            if self.union(canonical_existing, canonical_parent_id) {
                                total_unions += 1;
                            }
                        }
                    }

                    // Propagate child free_vars up to parent class
                    let mut child_vars = Vec::new();
                    for child_id in parent_node.children() {
                        let cr = self.union_find.find(child_id);
                        if let Some(Some(child_class)) = self.classes.get(cr.as_usize()) {
                            child_vars.extend(child_class.free_vars.iter().cloned());
                        }
                    }
                    if let Some(Some(parent_class)) = self.classes.get_mut(canonical_parent_id.as_usize()) {
                        parent_class.free_vars.extend(child_vars);
                    }

                    self.memo.insert(parent_node.clone(), canonical_parent_id);
                    new_parents.push((parent_node, canonical_parent_id));
                }

                if let Some(Some(class)) = self.classes.get_mut(canonical_winner.as_usize()) {
                    class.parents = new_parents;
                }
            }
        }

        // Final deduplication and canonicalization pass over memo keys
        let old_memo = std::mem::take(&mut self.memo);
        for (mut node, id) in old_memo {
            node.canonicalize(&mut self.union_find);
            let canonical_id = self.union_find.find(id);
            self.memo.insert(node, canonical_id);
        }

        total_unions
    }

    /// Returns the total count of active E-Classes.
    #[must_use]
    pub fn total_classes(&self) -> usize {
        self.classes.iter().filter(|c| c.is_some()).count()
    }

    /// Returns the total count of E-Nodes across all active classes.
    #[must_use]
    pub fn total_nodes(&self) -> usize {
        self.classes
            .iter()
            .filter_map(Option::as_ref)
            .map(EClass::len)
            .sum()
    }
}
