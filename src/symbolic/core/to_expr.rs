#![allow(deprecated)]


use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;

use super::dag_mgr::DagNode;
use super::dag_mgr::DagOp;
use super::expr::Distribution;
use super::expr::Expr;
use crate::symbolic::unit_unification::UnitQuantity;

impl DagNode {
    /// Converts a DAG node into an `Expr` (AST) structure.
    ///
    /// This is an iterative implementation to avoid stack overflow on deep trees.
    ///
    /// # Errors
    /// Returns an error if the node limit is exceeded or if any child cannot be converted.

    pub fn to_expr(
        &self
    ) -> Result<Expr, String> {

        use std::collections::HashMap;

        // Iterative implementation using explicit stack to prevent stack overflow
        // This uses a post-order (bottom-up) traversal strategy

        const MAX_NODES: usize =
            100_000;

        const MAX_CHILDREN: usize =
            10000;

        // Memoization: maps node hash to its converted Expr
        let mut memo: HashMap<
            u64,
            Expr,
        > = HashMap::new();

        // Work stack: nodes to process
        let mut work_stack: Vec<
            Arc<Self>,
        > = vec![Arc::new(
            self.clone(),
        )];

        // Track which nodes we've pushed to avoid cycles
        let mut visited: HashMap<
            u64,
            bool,
        > = HashMap::new();

        let mut nodes_processed = 0;

        while let Some(node) =
            work_stack.pop()
        {

            // Safety check: prevent processing too many nodes
            nodes_processed += 1;

            if nodes_processed
                > MAX_NODES
            {

                return Err(format!(
                    "Exceeded maximum \
                     node limit of \
                     {MAX_NODES}"
                ));
            }

            // If already converted, skip
            if memo.contains_key(
                &node.hash,
            ) {

                continue;
            }

            // Safety check: limit children count
            if node.children.len()
                > MAX_CHILDREN
            {

                return Err(format!(
                    "Node has too \
                     many children \
                     ({}), exceeds \
                     limit of {}",
                    node.children.len(),
                    MAX_CHILDREN
                ));
            }

            // Check if all children are already converted
            let children_ready = node
                .children
                .iter()
                .all(|child| {

                    memo.contains_key(
                        &child.hash,
                    )
                });

            if children_ready {

                // All children converted, now convert this node
                let children_exprs : Vec<Expr> = node
                    .children
                    .iter()
                    .filter_map(|child| {

                        memo.get(&child.hash)
                            .cloned()
                    })
                    .collect();

                // Helper macro to create Arc from children_exprs
                macro_rules! arc {
                    ($idx:expr_2021) => {

                        if $idx < children_exprs.len() {

                            Arc::new(children_exprs[$idx].clone())
                        } else {

                            return Err(format!(
                                "Index {} out of bounds for children array of length {}",
                                $idx,
                                children_exprs.len()
                            ));
                        }
                    };
                }

                let expr = match &node.op {
                    // --- Leaf Nodes ---
                    | DagOp::Constant(c) => Expr::Constant(c.into_inner()),
                    | DagOp::BigInt(i) => Expr::BigInt(i.clone()),
                    | DagOp::Rational(r) => Expr::Rational(r.clone()),
                    | DagOp::Boolean(b) => Expr::Boolean(*b),
                    | DagOp::Variable(s) => Expr::Variable(s.clone()),
                    | DagOp::Pattern(s) => Expr::Pattern(s.clone()),
                    | DagOp::Domain(s) => Expr::Domain(s.clone()),
                    | DagOp::Pi => Expr::Pi,
                    | DagOp::E => Expr::E,
                    | DagOp::Infinity => Expr::Infinity,
                    | DagOp::NegativeInfinity => Expr::NegativeInfinity,
                    | DagOp::InfiniteSolutions => Expr::InfiniteSolutions,
                    | DagOp::NoSolution => Expr::NoSolution,

                    // --- Operators with associated data ---
                    | DagOp::Derivative(s) => {

                        if children_exprs.is_empty() {

                            return Err(
                                "Derivative operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::Derivative(arc!(0), s.clone())
                    },
                    | DagOp::DerivativeN(s) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "DerivativeN operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::DerivativeN(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                        )
                    },
                    | DagOp::Limit(s) => {

                        if children_exprs.len() < 2 {

                            return Err("Limit operator requires at least 2 children".to_string());
                        }

                        Expr::Limit(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                        )
                    },
                    | DagOp::Solve(s) => {

                        if children_exprs.is_empty() {

                            return Err("Solve operator requires at least 1 child".to_string());
                        }

                        Expr::Solve(arc!(0), s.clone())
                    },
                    | DagOp::ConvergenceAnalysis(s) => {

                        if children_exprs.is_empty() {

                            return Err(
                                "ConvergenceAnalysis operator requires at least 1 child"
                                    .to_string(),
                            );
                        }

                        Expr::ConvergenceAnalysis(arc!(0), s.clone())
                    },
                    | DagOp::ForAll(s) => {

                        if children_exprs.is_empty() {

                            return Err("ForAll operator requires at least 1 child".to_string());
                        }

                        Expr::ForAll(s.clone(), arc!(0))
                    },
                    | DagOp::Exists(s) => {

                        if children_exprs.is_empty() {

                            return Err("Exists operator requires at least 1 child".to_string());
                        }

                        Expr::Exists(s.clone(), arc!(0))
                    },
                    | DagOp::Substitute(s) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Substitute operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Substitute(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                        )
                    },
                    | DagOp::Ode {
                        func,
                        var,
                    } => {

                        if children_exprs.is_empty() {

                            return Err("Ode operator requires at least 1 child".to_string());
                        }

                        Expr::Ode {
                            equation : arc!(0),
                            func : func.clone(),
                            var : var.clone(),
                        }
                    },
                    | DagOp::Pde {
                        func,
                        vars,
                    } => {

                        if children_exprs.is_empty() {

                            return Err("Pde operator requires at least 1 child".to_string());
                        }

                        Expr::Pde {
                            equation : arc!(0),
                            func : func.clone(),
                            vars : vars.clone(),
                        }
                    },
                    | DagOp::Predicate {
                        name,
                    } => {
                        Expr::Predicate {
                            name : name.clone(),
                            args : children_exprs.clone(),
                        }
                    },
                    | DagOp::Path(pt) => {

                        if children_exprs.len() < 2 {

                            return Err("Path operator requires at least 2 children".to_string());
                        }

                        Expr::Path(
                            pt.clone(),
                            arc!(0),
                            arc!(1),
                        )
                    },
                    | DagOp::Interval(incl_lower, incl_upper) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Interval operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Interval(
                            arc!(0),
                            arc!(1),
                            *incl_lower,
                            *incl_upper,
                        )
                    },
                    | DagOp::RootOf {
                        index,
                    } => {

                        if children_exprs.is_empty() {

                            return Err("RootOf operator requires at least 1 child".to_string());
                        }

                        Expr::RootOf {
                            poly : arc!(0),
                            index : *index,
                        }
                    },
                    | DagOp::SparsePolynomial(p) => Expr::SparsePolynomial(p.clone()),
                    | DagOp::QuantityWithValue(u) => {

                        if children_exprs.is_empty() {

                            return Err(
                                "QuantityWithValue operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::QuantityWithValue(arc!(0), u.clone())
                    },

                    // --- Binary operators ---
                    | DagOp::Add => {

                        if children_exprs.len() < 2 {

                            return Err("Add operator requires at least 2 children".to_string());
                        }

                        if children_exprs.len() == 2 {

                            Expr::Add(arc!(0), arc!(1))
                        } else {

                            Expr::AddList(children_exprs.clone())
                        }
                    },
                    | DagOp::Sub => {

                        if children_exprs.len() < 2 {

                            return Err("Sub operator requires at least 2 children".to_string());
                        }

                        Expr::Sub(arc!(0), arc!(1))
                    },
                    | DagOp::Mul => {

                        if children_exprs.len() < 2 {

                            return Err("Mul operator requires at least 2 children".to_string());
                        }

                        if children_exprs.len() == 2 {

                            Expr::Mul(arc!(0), arc!(1))
                        } else {

                            Expr::MulList(children_exprs.clone())
                        }
                    },
                    | DagOp::Div => {

                        if children_exprs.len() < 2 {

                            return Err("Div operator requires at least 2 children".to_string());
                        }

                        Expr::Div(arc!(0), arc!(1))
                    },
                    | DagOp::Eq => {

                        if children_exprs.len() < 2 {

                            return Err("Eq operator requires at least 2 children".to_string());
                        }

                        Expr::Eq(arc!(0), arc!(1))
                    },
                    | DagOp::Lt => {

                        if children_exprs.len() < 2 {

                            return Err("Lt operator requires at least 2 children".to_string());
                        }

                        Expr::Lt(arc!(0), arc!(1))
                    },
                    | DagOp::Gt => {

                        if children_exprs.len() < 2 {

                            return Err("Gt operator requires at least 2 children".to_string());
                        }

                        Expr::Gt(arc!(0), arc!(1))
                    },
                    | DagOp::Le => {

                        if children_exprs.len() < 2 {

                            return Err("Le operator requires at least 2 children".to_string());
                        }

                        Expr::Le(arc!(0), arc!(1))
                    },
                    | DagOp::Ge => {

                        if children_exprs.len() < 2 {

                            return Err("Ge operator requires at least 2 children".to_string());
                        }

                        Expr::Ge(arc!(0), arc!(1))
                    },
                    | DagOp::LogBase => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "LogBase operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::LogBase(arc!(0), arc!(1))
                    },
                    | DagOp::Atan2 => {

                        if children_exprs.len() < 2 {

                            return Err("Atan2 operator requires at least 2 children".to_string());
                        }

                        Expr::Atan2(arc!(0), arc!(1))
                    },
                    | DagOp::Binomial => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Binomial operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Binomial(arc!(0), arc!(1))
                    },
                    | DagOp::Beta => {

                        if children_exprs.len() < 2 {

                            return Err("Beta operator requires at least 2 children".to_string());
                        }

                        Expr::Beta(arc!(0), arc!(1))
                    },
                    | DagOp::BesselJ => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "BesselJ operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::BesselJ(arc!(0), arc!(1))
                    },
                    | DagOp::BesselY => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "BesselY operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::BesselY(arc!(0), arc!(1))
                    },
                    | DagOp::LegendreP => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "LegendreP operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::LegendreP(arc!(0), arc!(1))
                    },
                    | DagOp::LaguerreL => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "LaguerreL operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::LaguerreL(arc!(0), arc!(1))
                    },
                    | DagOp::HermiteH => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "HermiteH operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::HermiteH(arc!(0), arc!(1))
                    },
                    | DagOp::KroneckerDelta => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "KroneckerDelta operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::KroneckerDelta(arc!(0), arc!(1))
                    },
                    | DagOp::Permutation => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Permutation operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Permutation(arc!(0), arc!(1))
                    },
                    | DagOp::Combination => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Combination operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Combination(arc!(0), arc!(1))
                    },
                    | DagOp::FallingFactorial => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "FallingFactorial operator requires at least 2 children"
                                    .to_string(),
                            );
                        }

                        Expr::FallingFactorial(arc!(0), arc!(1))
                    },
                    | DagOp::RisingFactorial => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "RisingFactorial operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::RisingFactorial(arc!(0), arc!(1))
                    },
                    | DagOp::Xor => {

                        if children_exprs.len() < 2 {

                            return Err("Xor operator requires at least 2 children".to_string());
                        }

                        Expr::Xor(arc!(0), arc!(1))
                    },
                    | DagOp::Implies => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Implies operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Implies(arc!(0), arc!(1))
                    },
                    | DagOp::Equivalent => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Equivalent operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Equivalent(arc!(0), arc!(1))
                    },
                    | DagOp::Gcd => {

                        if children_exprs.len() < 2 {

                            return Err("Gcd operator requires at least 2 children".to_string());
                        }

                        Expr::Gcd(arc!(0), arc!(1))
                    },
                    | DagOp::Mod => {

                        if children_exprs.len() < 2 {

                            return Err("Mod operator requires at least 2 children".to_string());
                        }

                        Expr::Mod(arc!(0), arc!(1))
                    },
                    | DagOp::Max => {

                        if children_exprs.len() < 2 {

                            return Err("Max operator requires at least 2 children".to_string());
                        }

                        Expr::Max(arc!(0), arc!(1))
                    },
                    | DagOp::MatrixMul => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "MatrixMul operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::MatrixMul(arc!(0), arc!(1))
                    },
                    | DagOp::MatrixVecMul => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "MatrixVecMul operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::MatrixVecMul(arc!(0), arc!(1))
                    },
                    | DagOp::Apply => {

                        if children_exprs.len() < 2 {

                            return Err("Apply operator requires at least 2 children".to_string());
                        }

                        Expr::Apply(arc!(0), arc!(1))
                    },

                    // --- Unary operators ---
                    | DagOp::Neg => {

                        if children_exprs.is_empty() {

                            return Err("Neg operator requires at least 1 child".to_string());
                        }

                        Expr::Neg(arc!(0))
                    },
                    | DagOp::Power => {

                        if children_exprs.len() < 2 {

                            return Err("Power operator requires at least 2 children".to_string());
                        }

                        Expr::Power(arc!(0), arc!(1))
                    },
                    | DagOp::Sin => {

                        if children_exprs.is_empty() {

                            return Err("Sin operator requires at least 1 child".to_string());
                        }

                        Expr::Sin(arc!(0))
                    },
                    | DagOp::Cos => {

                        if children_exprs.is_empty() {

                            return Err("Cos operator requires at least 1 child".to_string());
                        }

                        Expr::Cos(arc!(0))
                    },
                    | DagOp::Tan => {

                        if children_exprs.is_empty() {

                            return Err("Tan operator requires at least 1 child".to_string());
                        }

                        Expr::Tan(arc!(0))
                    },
                    | DagOp::Exp => {

                        if children_exprs.is_empty() {

                            return Err("Exp operator requires at least 1 child".to_string());
                        }

                        Expr::Exp(arc!(0))
                    },
                    | DagOp::Log => {

                        if children_exprs.is_empty() {

                            return Err("Log operator requires at least 1 child".to_string());
                        }

                        Expr::Log(arc!(0))
                    },
                    | DagOp::Abs => {

                        if children_exprs.is_empty() {

                            return Err("Abs operator requires at least 1 child".to_string());
                        }

                        Expr::Abs(arc!(0))
                    },
                    | DagOp::Sqrt => {

                        if children_exprs.is_empty() {

                            return Err("Sqrt operator requires at least 1 child".to_string());
                        }

                        Expr::Sqrt(arc!(0))
                    },
                    | DagOp::Transpose => {

                        if children_exprs.is_empty() {

                            return Err("Transpose operator requires at least 1 child".to_string());
                        }

                        Expr::Transpose(arc!(0))
                    },
                    | DagOp::Inverse => {

                        if children_exprs.is_empty() {

                            return Err("Inverse operator requires at least 1 child".to_string());
                        }

                        Expr::Inverse(arc!(0))
                    },
                    | DagOp::Sec => {

                        if children_exprs.is_empty() {

                            return Err("Sec operator requires at least 1 child".to_string());
                        }

                        Expr::Sec(arc!(0))
                    },
                    | DagOp::Csc => {

                        if children_exprs.is_empty() {

                            return Err("Csc operator requires at least 1 child".to_string());
                        }

                        Expr::Csc(arc!(0))
                    },
                    | DagOp::Cot => {

                        if children_exprs.is_empty() {

                            return Err("Cot operator requires at least 1 child".to_string());
                        }

                        Expr::Cot(arc!(0))
                    },
                    | DagOp::ArcSin => {

                        if children_exprs.is_empty() {

                            return Err("ArcSin operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSin(arc!(0))
                    },
                    | DagOp::ArcCos => {

                        if children_exprs.is_empty() {

                            return Err("ArcCos operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCos(arc!(0))
                    },
                    | DagOp::ArcTan => {

                        if children_exprs.is_empty() {

                            return Err("ArcTan operator requires at least 1 child".to_string());
                        }

                        Expr::ArcTan(arc!(0))
                    },
                    | DagOp::ArcSec => {

                        if children_exprs.is_empty() {

                            return Err("ArcSec operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSec(arc!(0))
                    },
                    | DagOp::ArcCsc => {

                        if children_exprs.is_empty() {

                            return Err("ArcCsc operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCsc(arc!(0))
                    },
                    | DagOp::ArcCot => {

                        if children_exprs.is_empty() {

                            return Err("ArcCot operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCot(arc!(0))
                    },
                    | DagOp::Sinh => {

                        if children_exprs.is_empty() {

                            return Err("Sinh operator requires at least 1 child".to_string());
                        }

                        Expr::Sinh(arc!(0))
                    },
                    | DagOp::Cosh => {

                        if children_exprs.is_empty() {

                            return Err("Cosh operator requires at least 1 child".to_string());
                        }

                        Expr::Cosh(arc!(0))
                    },
                    | DagOp::Tanh => {

                        if children_exprs.is_empty() {

                            return Err("Tanh operator requires at least 1 child".to_string());
                        }

                        Expr::Tanh(arc!(0))
                    },
                    | DagOp::Sech => {

                        if children_exprs.is_empty() {

                            return Err("Sech operator requires at least 1 child".to_string());
                        }

                        Expr::Sech(arc!(0))
                    },
                    | DagOp::Csch => {

                        if children_exprs.is_empty() {

                            return Err("Csch operator requires at least 1 child".to_string());
                        }

                        Expr::Csch(arc!(0))
                    },
                    | DagOp::Coth => {

                        if children_exprs.is_empty() {

                            return Err("Coth operator requires at least 1 child".to_string());
                        }

                        Expr::Coth(arc!(0))
                    },
                    | DagOp::ArcSinh => {

                        if children_exprs.is_empty() {

                            return Err("ArcSinh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSinh(arc!(0))
                    },
                    | DagOp::ArcCosh => {

                        if children_exprs.is_empty() {

                            return Err("ArcCosh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCosh(arc!(0))
                    },
                    | DagOp::ArcTanh => {

                        if children_exprs.is_empty() {

                            return Err("ArcTanh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcTanh(arc!(0))
                    },
                    | DagOp::ArcSech => {

                        if children_exprs.is_empty() {

                            return Err("ArcSech operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSech(arc!(0))
                    },
                    | DagOp::ArcCsch => {

                        if children_exprs.is_empty() {

                            return Err("ArcCsch operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCsch(arc!(0))
                    },
                    | DagOp::ArcCoth => {

                        if children_exprs.is_empty() {

                            return Err("ArcCoth operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCoth(arc!(0))
                    },
                    | DagOp::Factorial => {

                        if children_exprs.is_empty() {

                            return Err("Factorial operator requires at least 1 child".to_string());
                        }

                        Expr::Factorial(arc!(0))
                    },
                    | DagOp::Boundary => {

                        if children_exprs.is_empty() {

                            return Err("Boundary operator requires at least 1 child".to_string());
                        }

                        Expr::Boundary(arc!(0))
                    },
                    | DagOp::Gamma => {

                        if children_exprs.is_empty() {

                            return Err("Gamma operator requires at least 1 child".to_string());
                        }

                        Expr::Gamma(arc!(0))
                    },
                    | DagOp::Erf => {

                        if children_exprs.is_empty() {

                            return Err("Erf operator requires at least 1 child".to_string());
                        }

                        Expr::Erf(arc!(0))
                    },
                    | DagOp::Erfc => {

                        if children_exprs.is_empty() {

                            return Err("Erfc operator requires at least 1 child".to_string());
                        }

                        Expr::Erfc(arc!(0))
                    },
                    | DagOp::Erfi => {

                        if children_exprs.is_empty() {

                            return Err("Erfi operator requires at least 1 child".to_string());
                        }

                        Expr::Erfi(arc!(0))
                    },
                    | DagOp::Zeta => {

                        if children_exprs.is_empty() {

                            return Err("Zeta operator requires at least 1 child".to_string());
                        }

                        Expr::Zeta(arc!(0))
                    },
                    | DagOp::Digamma => {

                        if children_exprs.is_empty() {

                            return Err("Digamma operator requires at least 1 child".to_string());
                        }

                        Expr::Digamma(arc!(0))
                    },
                    | DagOp::Not => {

                        if children_exprs.is_empty() {

                            return Err("Not operator requires at least 1 child".to_string());
                        }

                        Expr::Not(arc!(0))
                    },
                    | DagOp::Floor => {

                        if children_exprs.is_empty() {

                            return Err("Floor operator requires at least 1 child".to_string());
                        }

                        Expr::Floor(arc!(0))
                    },
                    | DagOp::IsPrime => {

                        if children_exprs.is_empty() {

                            return Err("IsPrime operator requires at least 1 child".to_string());
                        }

                        Expr::IsPrime(arc!(0))
                    },
                    | DagOp::GeneralSolution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "GeneralSolution operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::GeneralSolution(arc!(0))
                    },
                    | DagOp::ParticularSolution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "ParticularSolution operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::ParticularSolution(arc!(0))
                    },

                    // --- Complex structures ---
                    | DagOp::Matrix {
                        rows: _,
                        cols,
                    } => {

                        if children_exprs
                            .len()
                            .is_multiple_of(*cols) {

                            let reconstructed_matrix : Vec<Vec<Expr>> = children_exprs
                                .chunks(*cols)
                                .map(<[Expr]>::to_vec)
                                .collect();

                            Expr::Matrix(reconstructed_matrix)
                        } else {

                            let complete_rows = (children_exprs.len() / cols) * cols;

                            let reconstructed_matrix : Vec<Vec<Expr>> = children_exprs
                                .iter()
                                .take(complete_rows)
                                .cloned()
                                .collect::<Vec<_>>()
                                .chunks(*cols)
                                .map(<[Expr]>::to_vec)
                                .collect();

                            Expr::Matrix(reconstructed_matrix)
                        }
                    },
                    | DagOp::Complex => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Complex operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::Complex(arc!(0), arc!(1))
                    },
                    | DagOp::Integral => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Integral operator requires at least 4 children".to_string(),
                            );
                        }

                        Expr::Integral {
                            integrand : arc!(0),
                            var : arc!(1),
                            lower_bound : arc!(2),
                            upper_bound : arc!(3),
                        }
                    },
                    | DagOp::VolumeIntegral => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "VolumeIntegral operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::VolumeIntegral {
                            scalar_field : arc!(0),
                            volume : arc!(1),
                        }
                    },
                    | DagOp::SurfaceIntegral => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "SurfaceIntegral operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::SurfaceIntegral {
                            vector_field : arc!(0),
                            surface : arc!(1),
                        }
                    },
                    | DagOp::Sum => {

                        if children_exprs.len() < 4 {

                            return Err("Sum operator requires at least 4 children".to_string());
                        }

                        Expr::Sum {
                            body : arc!(0),
                            var : arc!(1),
                            from : arc!(2),
                            to : arc!(3),
                        }
                    },
                    | DagOp::Series(s) => {

                        if children_exprs.len() < 3 {

                            return Err("Series operator requires at least 3 children".to_string());
                        }

                        Expr::Series(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                            arc!(2),
                        )
                    },
                    | DagOp::Summation(s) => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "Summation operator requires at least 3 children".to_string(),
                            );
                        }

                        Expr::Summation(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                            arc!(2),
                        )
                    },
                    | DagOp::Product(s) => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "Product operator requires at least 3 children".to_string(),
                            );
                        }

                        Expr::Product(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                            arc!(2),
                        )
                    },
                    | DagOp::AsymptoticExpansion(s) => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "AsymptoticExpansion operator requires at least 3 children"
                                    .to_string(),
                            );
                        }

                        Expr::AsymptoticExpansion(
                            arc!(0),
                            s.clone(),
                            arc!(1),
                            arc!(2),
                        )
                    },
                    | DagOp::ParametricSolution => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "ParametricSolution operator requires at least 2 children"
                                    .to_string(),
                            );
                        }

                        Expr::ParametricSolution {
                            x : arc!(0),
                            y : arc!(1),
                        }
                    },
                    | DagOp::Fredholm => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Fredholm operator requires at least 4 children".to_string(),
                            );
                        }

                        Expr::Fredholm(
                            arc!(0),
                            arc!(1),
                            arc!(2),
                            arc!(3),
                        )
                    },
                    | DagOp::Volterra => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Volterra operator requires at least 4 children".to_string(),
                            );
                        }

                        Expr::Volterra(
                            arc!(0),
                            arc!(1),
                            arc!(2),
                            arc!(3),
                        )
                    },
                    | DagOp::Distribution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "Distribution operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::Distribution(children_exprs[0].clone_box_dist()?)
                    },
                    | DagOp::Quantity => {

                        if children_exprs.is_empty() {

                            return Err("Quantity operator requires at least 1 child".to_string());
                        }

                        Expr::Quantity(children_exprs[0].clone_box_quant()?)
                    },

                    // --- List operators ---
                    | DagOp::Vector => Expr::Vector(children_exprs.clone()),
                    | DagOp::And => Expr::And(children_exprs.clone()),
                    | DagOp::Or => Expr::Or(children_exprs.clone()),
                    | DagOp::Union => Expr::Union(children_exprs.clone()),
                    | DagOp::Polynomial => Expr::Polynomial(children_exprs.clone()),
                    | DagOp::System => Expr::System(children_exprs.clone()),
                    | DagOp::Solutions => Expr::Solutions(children_exprs.clone()),
                    | DagOp::Tuple => Expr::Tuple(children_exprs.clone()),

                    // --- Custom ---
                    | DagOp::CustomZero => Expr::CustomZero,
                    | DagOp::CustomString(s) => Expr::CustomString(s.clone()),
                    | DagOp::CustomArcOne => {

                        if children_exprs.is_empty() {

                            return Err(
                                "CustomArcOne operator requires at least 1 child".to_string(),
                            );
                        }

                        Expr::CustomArcOne(arc!(0))
                    },
                    | DagOp::CustomArcTwo => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "CustomArcTwo operator requires at least 2 children".to_string(),
                            );
                        }

                        Expr::CustomArcTwo(arc!(0), arc!(1))
                    },
                    | DagOp::CustomArcThree => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "CustomArcThree operator requires at least 3 children".to_string(),
                            );
                        }

                        Expr::CustomArcThree(
                            arc!(0),
                            arc!(1),
                            arc!(2),
                        )
                    },
                    | DagOp::CustomArcFour => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "CustomArcFour operator requires at least 4 children".to_string(),
                            );
                        }

                        Expr::CustomArcFour(
                            arc!(0),
                            arc!(1),
                            arc!(2),
                            arc!(3),
                        )
                    },
                    | DagOp::CustomArcFive => {

                        if children_exprs.len() < 5 {

                            return Err(
                                "CustomArcFive operator requires at least 5 children".to_string(),
                            );
                        }

                        Expr::CustomArcFive(
                            arc!(0),
                            arc!(1),
                            arc!(2),
                            arc!(3),
                            arc!(4),
                        )
                    },
                    | DagOp::CustomVecOne => Expr::CustomVecOne(children_exprs.clone()),
                    | DagOp::CustomVecTwo => {
                        return Err("CustomVecTwo to_expr is ambiguous".to_string())
                    },
                    | DagOp::CustomVecThree => {
                        return Err("CustomVecThree to_expr is ambiguous".to_string())
                    },
                    | DagOp::CustomVecFour => {
                        return Err("CustomVecFour to_expr is ambiguous".to_string())
                    },
                    | DagOp::CustomVecFive => {
                        return Err("CustomVecFive to_expr is ambiguous".to_string())
                    },

                    | DagOp::UnaryList(s) => {

                        if children_exprs.is_empty() {

                            return Err(format!(
                                "UnaryList operator {s} requires at least 1 child"
                            ));
                        }

                        Expr::UnaryList(s.clone(), arc!(0))
                    },
                    | DagOp::BinaryList(s) => {

                        if children_exprs.len() < 2 {

                            return Err(format!(
                                "BinaryList operator {s} requires at least 2 children"
                            ));
                        }

                        Expr::BinaryList(
                            s.clone(),
                            arc!(0),
                            arc!(1),
                        )
                    },
                    | DagOp::NaryList(s) => {
                        Expr::NaryList(
                            s.clone(),
                            children_exprs.clone(),
                        )
                    },
                };

                // Store the converted expression
                memo.insert(
                    node.hash,
                    expr,
                );
            } else {

                // Not all children ready, push node back and push children
                work_stack
                    .push(node.clone());

                // Push children in reverse order (so they're processed in correct order)
                if visited
                    .insert(
                        node.hash,
                        true,
                    )
                    .is_none()
                {

                    for child in node
                        .children
                        .iter()
                        .rev()
                    {

                        if !memo.contains_key(&child.hash) {

                            work_stack.push(child.clone());
                        }
                    }
                }
            }
        }

        // Return the converted expression for the root node
        memo.get(&self.hash)
            .cloned()
            .ok_or_else(|| {

                "Failed to convert \
                 root node"
                    .to_string()
            })
    }

    #[must_use]

    /// Creates a new `DagNode` with the given operation and children.
    ///
    /// Automatically handles hashing and enforces safety limits on the number of children.

    pub fn new(
        op: DagOp,
        children: Vec<Arc<Self>>,
    ) -> Arc<Self> {

        // Safety check: limit number of children to prevent excessive memory allocation
        const MAX_CHILDREN: usize =
            10000;

        if children.len() > MAX_CHILDREN
        {

            // This should not happen in normal usage, but we handle it gracefully
            // by truncating the children list - this is a defensive programming approach
            let safe_children: Vec<_> =
                children
                    .into_iter()
                    .take(MAX_CHILDREN)
                    .collect();

            let mut hasher = std::collections::hash_map::DefaultHasher::new();

            op.hash(&mut hasher);

            safe_children
                .hash(&mut hasher);

            let hash = hasher.finish();

            return Arc::new(Self {
                op,
                children: safe_children,
                hash,
            });
        }

        let mut hasher = std::collections::hash_map::DefaultHasher::new();

        op.hash(&mut hasher);

        children.hash(&mut hasher);

        let hash = hasher.finish();

        Arc::new(Self {
            op,
            children,
            hash,
        })
    }
}

impl Expr {
    /// Clones the distribution if the expression represents a probability distribution.
    ///
    /// # Errors
    /// Returns an error if the expression is not a distribution.

    pub fn clone_box_dist(
        &self
    ) -> Result<
        Arc<dyn Distribution>,
        String,
    > {

        if let Self::Distribution(d) =
            self
        {

            Ok(d.clone_box())
        } else {

            Err("Cannot clone into \
                 Distribution"
                .to_string())
        }
    }

    /// Clones the unit quantity if the expression represents a quantity.
    ///
    /// # Errors
    /// Returns an error if the expression is not a quantity.

    pub fn clone_box_quant(
        &self
    ) -> Result<Arc<UnitQuantity>, String>
    {

        if let Self::Quantity(q) = self
        {

            Ok(q.clone())
        } else {

            Err("Cannot clone into \
                 UnitQuantity"
                .to_string())
        }
    }
}
