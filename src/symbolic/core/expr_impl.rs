#![allow(deprecated)]
#![allow(clippy::match_same_arms)]

use std::cmp::Ordering;
use std::convert::AsRef;
use std::fmt::Debug;
use std::fmt::Write;
use std::fmt::{
    self,
};
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;

use num_bigint::Sign;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;

use super::api::get_dynamic_op_properties;
use super::dag_mgr::DagOp;
use super::expr::Expr;

impl PartialEq for Expr {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {

        if let (
            Self::Dag(n1),
            Self::Dag(n2),
        ) = (self, other)
        {

            if Arc::ptr_eq(n1, n2) {

                return true;
            }
        }

        let op1 = self.op();
        let op2 = other.op();

        if op1 != op2 {
            // Allow numeric cross-comparison even if the operations differ,
            // as long as both sides represent equivalent numerical values.
            // We use a safe version that avoids recursive DAG conversion to prevent stack overflow (panic).
            if let (Some(f1), Some(f2)) = (get_numeric_value_efficient(self), get_numeric_value_efficient(other)) {
                return (f1 - f2).abs() < f64::EPSILON;
            }

            return false;
        }

        match (self, other) {
            // --- COMMUTATIVE OPERATORS (A+B == B+A) ---
            | (Self::Add(l1, r1), Self::Add(l2, r2)) | (Self::Mul(l1, r1), Self::Mul(l2, r2)) => {

                // Check for (l1 == l2 AND r1 == r2) OR (l1 == r2 AND r1 == l2)
                let standard_order_match = l1
                    .as_ref()
                    .eq(l2.as_ref())
                    && r1
                        .as_ref()
                        .eq(r2.as_ref());

                let inverse_order_match = l1
                    .as_ref()
                    .eq(r2.as_ref())
                    && r1
                        .as_ref()
                        .eq(l2.as_ref());

                return standard_order_match || inverse_order_match;
            },

            // --- NON-COMMUTATIVE OPERATORS (A-B != B-A) ---
            | (Self::Sub(l1, r1), Self::Sub(l2, r2))
            | (Self::Div(l1, r1), Self::Div(l2, r2))
            | (Self::Power(l1, r1), Self::Power(l2, r2)) => {

                // Positional comparison is required for non-commutative ops
                return l1
                    .as_ref()
                    .eq(l2.as_ref())
                    && r1
                        .as_ref()
                        .eq(r2.as_ref());
            },

            // Special handling for Derivative to compare both expression and variable
            | (Self::Derivative(e1, v1), Self::Derivative(e2, v2)) => {
                return v1 == v2
                    && e1
                        .as_ref()
                        .eq(e2.as_ref())
            },

            // Special handling for other variants with String parameters
            | (Self::Solve(e1, v1), Self::Solve(e2, v2))
            | (Self::ConvergenceAnalysis(e1, v1), Self::ConvergenceAnalysis(e2, v2))
            | (Self::ForAll(v1, e1), Self::ForAll(v2, e2))
            | (Self::Exists(v1, e1), Self::Exists(v2, e2)) => {
                return v1 == v2
                    && e1
                        .as_ref()
                        .eq(e2.as_ref())
            },

            | (Self::Constant(f1), Self::Constant(f2)) => return (f1 - f2).abs() < f64::EPSILON,
            | (Self::BigInt(b1), Self::BigInt(b2)) => return b1 == b2,
            | (Self::Rational(r1), Self::Rational(r2)) => return r1 == r2,

            // BigInt <=> Rational
            | (Self::BigInt(b), Self::Rational(r)) | (Self::Rational(r), Self::BigInt(b)) => {

                let temp_rational = BigRational::from(b.clone());

                return r == &temp_rational;
            },

            // BigInt / Rational <=> Constant(f64)
            | (Self::Constant(f), Self::Rational(r)) | (Self::Rational(r), Self::Constant(f)) => {
                match r.to_f64() {
                    | Some(r_f64) => return (f - r_f64).abs() < f64::EPSILON,
                    | None => return false,
                }
            },

            | (Self::Constant(f), Self::BigInt(b)) | (Self::BigInt(b), Self::Constant(f)) => {

                if f.fract().abs() < f64::EPSILON {

                    match b.to_f64() {
                        | Some(b_f64) => return (f - b_f64).abs() < f64::EPSILON,
                        | None => return false,
                    }
                }

                return false;
            },

            | _ => { /* Ignore, Enter Next Step */ },
        }

        let self_children =
            self.children();

        let other_children =
            other.children();

        if self_children.len()
            != other_children.len()
        {

            return false;
        }

        self_children
            .iter()
            .zip(other_children.iter())
            .all(
                |(
                    l_child_expr,
                    r_child_expr,
                )| {

                    l_child_expr.eq(
                        r_child_expr,
                    )
                },
            )
    }
}

impl Eq for Expr {
}

impl Hash for Expr {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {

        // Use the unified view
        let op = self.op();

        op.hash(state);

        let mut children =
            self.children();

        match op {
            | DagOp::Add
            | DagOp::Mul
            | DagOp::And
            | DagOp::Or
            | DagOp::Xor
            | DagOp::Equivalent
            | DagOp::Eq => {

                // Commutative: hash children in a specified order (sorted)
                // to ensure the hash is canonical.
                children.sort(); // Relies on Ord
                for child in children {

                    child.hash(state);
                }
            },
            | _ => {

                // Non-commutative: hash children in order.
                for child in children {

                    child.hash(state);
                }
            },
        }
    }
}

impl PartialOrd for Expr {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {

        Some(self.cmp(other))
    }
}

impl Ord for Expr {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {

        // Fast path for identical DAG nodes.
        if let (
            Self::Dag(n1),
            Self::Dag(n2),
        ) = (self, other)
        {

            if Arc::ptr_eq(n1, n2) {

                return Ordering::Equal;
            }
        }

        // Compare by operator.
        let op_ordering = self
            .op()
            .cmp(&other.op());

        if op_ordering
            != Ordering::Equal
        {

            return op_ordering;
        }

        // For canonical DAG nodes, children are already sorted, so we can compare directly.
        // For non-DAG nodes or mixed comparisons, this relies on the slow sorting path in the
        // Ord implementation of the children expressions.
        self.children()
            .cmp(&other.children())
    }
}

/// Custom error type for symbolic operations.
#[derive(Debug)]

pub enum SymbolicError {
    /// Error message variant.
    Msg(String),
}

impl fmt::Display for SymbolicError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {

        match self {
            | Self::Msg(s) => {

                write!(f, "{s}")
            },
        }
    }
}

impl From<String> for SymbolicError {
    fn from(s: String) -> Self {

        Self::Msg(s)
    }
}

impl From<&str> for SymbolicError {
    fn from(s: &str) -> Self {

        Self::Msg(s.to_string())
    }
}

impl Expr {
    /// Performs a pre-order traversal of the expression tree.
    /// It visits the current node first, then its children.
    ///
    /// # Arguments
    /// * `f` - A mutable function that takes a reference to an `Expr` and is applied to each node during traversal

    pub fn pre_order_walk<F>(
        &self,
        f: &mut F,
    ) where
        F: FnMut(&Self),
    {

        f(self); // Visit parent
        match self {
            // Binary operators
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Eq(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Binomial(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::KroneckerDelta(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b)
            | Self::FallingFactorial(a, b)
            | Self::RisingFactorial(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Gcd(a, b)
            | Self::Mod(a, b)
            | Self::Max(a, b)
            | Self::MatrixMul(a, b)
            | Self::MatrixVecMul(a, b)
            | Self::Apply(a, b)
            | Self::Path(_, a, b) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            },
            // Unary operators
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Neg(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Sec(a)
            | Self::Csc(a)
            | Self::Cot(a)
            | Self::ArcSin(a)
            | Self::ArcCos(a)
            | Self::ArcTan(a)
            | Self::ArcSec(a)
            | Self::ArcCsc(a)
            | Self::ArcCot(a)
            | Self::Sinh(a)
            | Self::Cosh(a)
            | Self::Tanh(a)
            | Self::Sech(a)
            | Self::Csch(a)
            | Self::Coth(a)
            | Self::ArcSinh(a)
            | Self::ArcCosh(a)
            | Self::ArcTanh(a)
            | Self::ArcSech(a)
            | Self::ArcCsch(a)
            | Self::ArcCoth(a)
            | Self::Boundary(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Not(a)
            | Self::Floor(a)
            | Self::IsPrime(a)
            | Self::Factorial(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::GeneralSolution(a)
            | Self::ParticularSolution(a)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::ConvergenceAnalysis(a, _)
            | Self::ForAll(_, a)
            | Self::Exists(_, a) => {

                a.pre_order_walk(f);
            },
            // N-ary operators
            | Self::Matrix(m) => {
                for row in m {
                    for e in row {
                        e.pre_order_walk(f);
                    }
                }
            },
            | Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v)
            | Self::AddList(v)
            | Self::MulList(v) => {
                for e in v {
                    e.pre_order_walk(f);
                }
            },
            | Self::Predicate {
                args,
                ..
            } => {
                for e in args {
                    e.pre_order_walk(f);
                }
            },
            | Self::SparsePolynomial(p) => {
                for c in p.terms.values() {
                    c.pre_order_walk(f);
                }
            },
            // More complex operators
            | Self::Sum {
                body,
                var,
                from,
                to,
            } => {

                f(self);

                body.pre_order_walk(f);

                var.pre_order_walk(f);

                from.pre_order_walk(f);

                to.pre_order_walk(f);
            },
            | Self::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {

                integrand.pre_order_walk(f);

                lower_bound.pre_order_walk(f);

                upper_bound.pre_order_walk(f);
            },
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                scalar_field.pre_order_walk(f);

                volume.pre_order_walk(f);
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                vector_field.pre_order_walk(f);

                surface.pre_order_walk(f);
            },
            | Self::DerivativeN(e, _, n) => {

                e.pre_order_walk(f);

                n.pre_order_walk(f);
            },
            | Self::Series(a, _, c, d)
            | Self::Summation(a, _, c, d)
            | Self::Product(a, _, c, d) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            },
            | Self::AsymptoticExpansion(a, _, c, d) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            },
            | Self::Interval(a, b, _, _) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            },
            | Self::Substitute(a, _, c) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);
            },
            | Self::Limit(a, _, c) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);
            },
            | Self::Ode {
                equation,
                ..
            } => equation.pre_order_walk(f),
            | Self::Pde {
                equation,
                ..
            } => equation.pre_order_walk(f),
            | Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            },
            | Self::ParametricSolution {
                x,
                y,
            } => {

                x.pre_order_walk(f);

                y.pre_order_walk(f);
            },
            | Self::RootOf {
                poly,
                ..
            } => poly.pre_order_walk(f),
            | Self::QuantityWithValue(v, _) => v.pre_order_walk(f),

            | Self::CustomArcOne(a) => {

                a.pre_order_walk(f);
            },
            | Self::CustomArcTwo(a, b) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            },
            | Self::CustomArcThree(a, b, c) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);
            },
            | Self::CustomArcFour(a, b, c, d) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            },
            | Self::CustomArcFive(a, b, c, d, e) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);

                e.pre_order_walk(f);
            },
            | Self::CustomVecOne(v) => {
                for e in v {
                    e.pre_order_walk(f);
                }
            },
            | Self::CustomVecTwo(v1, v2) => {

                for e in v1 {

                    e.pre_order_walk(f);
                }

                for e in v2 {

                    e.pre_order_walk(f);
                }
            },
            | Self::CustomVecThree(v1, v2, v3) => {

                for e in v1 {

                    e.pre_order_walk(f);
                }

                for e in v2 {

                    e.pre_order_walk(f);
                }

                for e in v3 {

                    e.pre_order_walk(f);
                }
            },
            | Self::CustomVecFour(v1, v2, v3, v4) => {

                for e in v1 {

                    e.pre_order_walk(f);
                }

                for e in v2 {

                    e.pre_order_walk(f);
                }

                for e in v3 {

                    e.pre_order_walk(f);
                }

                for e in v4 {

                    e.pre_order_walk(f);
                }
            },
            | Self::CustomVecFive(v1, v2, v3, v4, v5) => {

                for e in v1 {

                    e.pre_order_walk(f);
                }

                for e in v2 {

                    e.pre_order_walk(f);
                }

                for e in v3 {

                    e.pre_order_walk(f);
                }

                for e in v4 {

                    e.pre_order_walk(f);
                }

                for e in v5 {

                    e.pre_order_walk(f);
                }
            },
            | Self::UnaryList(_, a) => a.pre_order_walk(f),
            | Self::BinaryList(_, a, b) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            },
            | Self::NaryList(_, v) => {
                for e in v {

                    e.pre_order_walk(f);
                }
            },

            | Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.pre_order_walk(f);
                }
            },
            // Leaf nodes
            | Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Domain(_)
            | Self::Pi
            | Self::Quantity(_)
            | Self::E
            | Self::Infinity
            | Self::NegativeInfinity
            | Self::InfiniteSolutions
            | Self::NoSolution
            | Self::CustomZero
            | Self::CustomString(_)
            | Self::Distribution(_) => {},
        }
    }

    /// Performs a post-order traversal of the expression tree.
    /// It visits the children first, then the current node.
    ///
    /// # Arguments
    /// * `f` - A mutable function that takes a reference to an `Expr` and is applied to each node during traversal

    pub fn post_order_walk<F>(
        &self,
        f: &mut F,
    ) where
        F: FnMut(&Self),
    {

        match self {
            // Binary operators
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Eq(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Binomial(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::KroneckerDelta(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b)
            | Self::FallingFactorial(a, b)
            | Self::RisingFactorial(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Gcd(a, b)
            | Self::Mod(a, b)
            | Self::Max(a, b)
            | Self::MatrixMul(a, b)
            | Self::MatrixVecMul(a, b)
            | Self::Apply(a, b)
            | Self::Path(_, a, b) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            },
            // Unary operators
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Neg(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Sec(a)
            | Self::Csc(a)
            | Self::Cot(a)
            | Self::ArcSin(a)
            | Self::ArcCos(a)
            | Self::ArcTan(a)
            | Self::ArcSec(a)
            | Self::ArcCsc(a)
            | Self::ArcCot(a)
            | Self::Sinh(a)
            | Self::Cosh(a)
            | Self::Tanh(a)
            | Self::Sech(a)
            | Self::Csch(a)
            | Self::Coth(a)
            | Self::ArcSinh(a)
            | Self::ArcCosh(a)
            | Self::ArcTanh(a)
            | Self::ArcSech(a)
            | Self::ArcCsch(a)
            | Self::ArcCoth(a)
            | Self::Boundary(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Not(a)
            | Self::Floor(a)
            | Self::IsPrime(a)
            | Self::Factorial(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::GeneralSolution(a)
            | Self::ParticularSolution(a)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::ConvergenceAnalysis(a, _)
            | Self::ForAll(_, a)
            | Self::Exists(_, a) => {

                a.post_order_walk(f);
            },
            // N-ary operators
            | Self::Matrix(m) => {
                for row in m {
                    for e in row {
                        e.post_order_walk(f);
                    }
                }
            },
            | Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v)
            | Self::AddList(v)
            | Self::MulList(v) => {
                for e in v {
                    e.post_order_walk(f);
                }
            },
            | Self::Predicate {
                args,
                ..
            } => {
                for e in args {
                    e.post_order_walk(f);
                }
            },
            | Self::SparsePolynomial(p) => {
                for c in p.terms.values() {
                    c.post_order_walk(f);
                }
            },
            // More complex operators
            | Self::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {

                integrand.post_order_walk(f);

                lower_bound.post_order_walk(f);

                upper_bound.post_order_walk(f);
            },
            | Self::Sum {
                body,
                var,
                from,
                to,
            } => {

                body.post_order_walk(f);

                var.post_order_walk(f);

                from.post_order_walk(f);

                to.post_order_walk(f);
            },
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                scalar_field.post_order_walk(f);

                volume.post_order_walk(f);
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                vector_field.post_order_walk(f);

                surface.post_order_walk(f);
            },
            | Self::DerivativeN(e, _, n) => {

                e.post_order_walk(f);

                n.post_order_walk(f);
            },
            | Self::Series(a, _, c, d)
            | Self::Summation(a, _, c, d)
            | Self::Product(a, _, c, d) => {

                a.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            },
            | Self::AsymptoticExpansion(a, _, c, d) => {

                a.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            },
            | Self::Interval(a, b, _, _) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            },
            | Self::Substitute(a, _, c) => {

                a.post_order_walk(f);

                c.post_order_walk(f);
            },
            | Self::Limit(a, _, c) => {

                a.post_order_walk(f);

                c.post_order_walk(f);
            },
            | Self::Ode {
                equation,
                ..
            } => equation.post_order_walk(f),
            | Self::Pde {
                equation,
                ..
            } => equation.post_order_walk(f),
            | Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            },
            | Self::ParametricSolution {
                x,
                y,
            } => {

                x.post_order_walk(f);

                y.post_order_walk(f);
            },
            | Self::QuantityWithValue(v, _) => v.post_order_walk(f),
            | Self::RootOf {
                poly,
                ..
            } => poly.post_order_walk(f),

            | Self::CustomArcOne(a) => {

                a.post_order_walk(f);
            },
            | Self::CustomArcTwo(a, b) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            },
            | Self::CustomArcThree(a, b, c) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);
            },
            | Self::CustomArcFour(a, b, c, d) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            },
            | Self::CustomArcFive(a, b, c, d, e) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);

                e.post_order_walk(f);
            },
            | Self::CustomVecOne(v)
            | Self::CustomVecTwo(v, _)
            | Self::CustomVecThree(v, _, _)
            | Self::CustomVecFour(v, _, _, _)
            | Self::CustomVecFive(v, _, _, _, _) => {
                for e in v {

                    e.post_order_walk(f);
                }
            },
            | Self::UnaryList(_, a) => a.post_order_walk(f),
            | Self::BinaryList(_, a, b) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            },
            | Self::NaryList(_, v) => {
                for e in v {

                    e.post_order_walk(f);
                }
            },
            | Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.post_order_walk(f);
                }
            },
            // Leaf nodes
            | Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Domain(_)
            | Self::Pi
            | Self::Quantity(_)
            | Self::E
            | Self::Infinity
            | Self::NegativeInfinity
            | Self::InfiniteSolutions
            | Self::NoSolution
            | Self::CustomZero
            | Self::CustomString(_)
            | Self::Distribution(_) => {},
        }

        f(self); // Visit parent
    }

    /// Performs an in-order traversal of the expression tree.
    /// For binary operators, it visits the left child, the node itself, then the right child.
    /// For other nodes, the behavior is adapted as it's not strictly defined.
    ///
    /// # Arguments
    /// * `f` - A mutable function that takes a reference to an `Expr` and is applied to each node during traversal

    pub fn in_order_walk<F>(
        &self,
        f: &mut F,
    ) where
        F: FnMut(&Self),
    {

        match self {
            // Binary operators
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Eq(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Binomial(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::KroneckerDelta(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b)
            | Self::FallingFactorial(a, b)
            | Self::RisingFactorial(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Gcd(a, b)
            | Self::Mod(a, b)
            | Self::Max(a, b)
            | Self::MatrixMul(a, b)
            | Self::MatrixVecMul(a, b)
            | Self::Apply(a, b)
            | Self::Path(_, a, b) => {

                a.in_order_walk(f);

                f(self);

                b.in_order_walk(f);
            },
            // Unary operators (treat as pre-order)
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Neg(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Sec(a)
            | Self::Csc(a)
            | Self::Cot(a)
            | Self::ArcSin(a)
            | Self::ArcCos(a)
            | Self::ArcTan(a)
            | Self::ArcSec(a)
            | Self::ArcCsc(a)
            | Self::ArcCot(a)
            | Self::Sinh(a)
            | Self::Cosh(a)
            | Self::Tanh(a)
            | Self::Sech(a)
            | Self::Csch(a)
            | Self::Coth(a)
            | Self::ArcSinh(a)
            | Self::ArcCosh(a)
            | Self::ArcTanh(a)
            | Self::ArcSech(a)
            | Self::ArcCsch(a)
            | Self::ArcCoth(a)
            | Self::Boundary(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Not(a)
            | Self::Floor(a)
            | Self::IsPrime(a)
            | Self::Factorial(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::GeneralSolution(a)
            | Self::ParticularSolution(a)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::ConvergenceAnalysis(a, _)
            | Self::ForAll(_, a)
            | Self::Exists(_, a) => {

                f(self);

                a.in_order_walk(f);
            },
            // N-ary operators (visit self, then children)
            | Self::Matrix(m) => {

                f(self);

                m.iter()
                    .flatten()
                    .for_each(|e| e.in_order_walk(f));
            },
            | Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v)
            | Self::AddList(v)
            | Self::MulList(v) => {

                f(self);

                for e in v {

                    e.in_order_walk(f);
                }
            },
            | Self::Predicate {
                args,
                ..
            } => {

                f(self);

                for e in args {

                    e.in_order_walk(f);
                }
            },
            | Self::SparsePolynomial(p) => {

                f(self);

                p.terms
                    .values()
                    .for_each(|c| c.in_order_walk(f));
            },
            // More complex operators (visit self, then children)
            | Self::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {

                f(self);

                integrand.in_order_walk(f);

                lower_bound.in_order_walk(f);

                upper_bound.in_order_walk(f);
            },
            | Self::Sum {
                body,
                var,
                from,
                to,
            } => {

                f(self);

                body.in_order_walk(f);

                var.in_order_walk(f);

                from.in_order_walk(f);

                to.in_order_walk(f);
            },
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                f(self);

                scalar_field.in_order_walk(f);

                volume.in_order_walk(f);
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                f(self);

                vector_field.in_order_walk(f);

                surface.in_order_walk(f);
            },
            | Self::DerivativeN(e, _, n) => {

                f(self);

                e.in_order_walk(f);

                n.in_order_walk(f);
            },
            | Self::Series(a, _, c, d)
            | Self::Summation(a, _, c, d)
            | Self::Product(a, _, c, d) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);

                d.in_order_walk(f);
            },
            | Self::AsymptoticExpansion(a, _, c, _d) => {

                f(self);

                a.in_order_walk(f);

                c.pre_order_walk(f);
            },
            | Self::Interval(a, b, _, _) => {

                f(self);

                a.in_order_walk(f);

                b.in_order_walk(f);
            },
            | Self::Substitute(a, _, c) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);
            },
            | Self::Limit(a, _, c) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);
            },
            | Self::Ode {
                equation,
                ..
            } => {

                f(self);

                equation.in_order_walk(f);
            },
            | Self::Pde {
                equation,
                ..
            } => {

                f(self);

                equation.in_order_walk(f);
            },
            | Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                f(self);

                a.in_order_walk(f);

                b.in_order_walk(f);

                c.in_order_walk(f);

                d.pre_order_walk(f);
            },
            | Self::ParametricSolution {
                x,
                y,
            } => {

                f(self);

                x.in_order_walk(f);

                y.in_order_walk(f);
            },
            | Self::RootOf {
                poly,
                ..
            } => {

                f(self);

                poly.in_order_walk(f);
            },
            | Self::QuantityWithValue(v, _) => v.in_order_walk(f),

            | Self::CustomArcOne(a) => {

                f(self);

                a.in_order_walk(f);
            },
            | Self::CustomArcTwo(a, b) => {

                a.in_order_walk(f);

                f(self);

                b.in_order_walk(f);
            },
            | Self::CustomArcThree(a, b, c) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);
            },
            | Self::CustomArcFour(a, b, c, d) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);

                d.in_order_walk(f);
            },
            | Self::CustomArcFive(a, b, c, d, e) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);

                d.in_order_walk(f);

                e.in_order_walk(f);
            },
            | Self::CustomVecOne(v) => {

                f(self);

                for e in v {

                    e.in_order_walk(f);
                }
            },
            | Self::CustomVecTwo(v1, v2) => {

                f(self);

                for e in v1 {

                    e.in_order_walk(f);
                }

                for e in v2 {

                    e.in_order_walk(f);
                }
            },
            | Self::CustomVecThree(v1, v2, v3) => {

                f(self);

                for e in v1 {

                    e.in_order_walk(f);
                }

                for e in v2 {

                    e.in_order_walk(f);
                }

                for e in v3 {

                    e.in_order_walk(f);
                }
            },
            | Self::CustomVecFour(v1, v2, v3, v4) => {

                f(self);

                for e in v1 {

                    e.in_order_walk(f);
                }

                for e in v2 {

                    e.in_order_walk(f);
                }

                for e in v3 {

                    e.in_order_walk(f);
                }

                for e in v4 {

                    e.in_order_walk(f);
                }
            },
            | Self::CustomVecFive(v1, v2, v3, v4, v5) => {

                f(self);

                for e in v1 {

                    e.in_order_walk(f);
                }

                for e in v2 {

                    e.in_order_walk(f);
                }

                for e in v3 {

                    e.in_order_walk(f);
                }

                for e in v4 {

                    e.in_order_walk(f);
                }

                for e in v5 {

                    e.in_order_walk(f);
                }
            },
            | Self::UnaryList(_, a) => {

                f(self);

                a.in_order_walk(f);
            },
            | Self::BinaryList(_, a, b) => {

                a.in_order_walk(f);

                f(self);

                b.in_order_walk(f);
            },
            | Self::NaryList(_, v) => {

                f(self);

                for e in v {

                    e.in_order_walk(f);
                }
            },

            // Leaf nodes
            | Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Domain(_)
            | Self::Pi
            | Self::Quantity(_)
            | Self::E
            | Self::Infinity
            | Self::NegativeInfinity
            | Self::InfiniteSolutions
            | Self::NoSolution => {},
            | Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.in_order_walk(f);
                }
            },
            | Self::CustomZero | Self::CustomString(_) | Self::Distribution(_) => {},
        }

        f(self); // Visit parent
    }

    /// Returns the direct children of this expression in a vector.
    ///
    /// For leaf nodes, returns an empty vector.
    /// For unary operations, returns a vector with one element.
    /// For binary operations, returns a vector with two elements, and so on.
    ///
    /// # Returns
    /// * `Vec<Expr>` - A vector containing the direct children of this expression

    pub(crate) fn get_children_internal(
        &self
    ) -> Vec<Self> {

        match self {
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Eq(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Binomial(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::KroneckerDelta(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b)
            | Self::FallingFactorial(a, b)
            | Self::RisingFactorial(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Gcd(a, b)
            | Self::Mod(a, b)
            | Self::Max(a, b)
            | Self::MatrixMul(a, b)
            | Self::MatrixVecMul(a, b)
            | Self::Apply(a, b) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ]
            },
            | Self::AddList(v) | Self::MulList(v) => v.clone(),
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Neg(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Sec(a)
            | Self::Csc(a)
            | Self::Cot(a)
            | Self::ArcSin(a)
            | Self::ArcCos(a)
            | Self::ArcTan(a)
            | Self::ArcSec(a)
            | Self::ArcCsc(a)
            | Self::ArcCot(a)
            | Self::Sinh(a)
            | Self::Cosh(a)
            | Self::Tanh(a)
            | Self::Sech(a)
            | Self::Csch(a)
            | Self::Coth(a)
            | Self::ArcSinh(a)
            | Self::ArcCosh(a)
            | Self::ArcTanh(a)
            | Self::ArcSech(a)
            | Self::ArcCsch(a)
            | Self::ArcCoth(a)
            | Self::Boundary(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Not(a)
            | Self::Floor(a)
            | Self::IsPrime(a)
            | Self::Factorial(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::GeneralSolution(a)
            | Self::ParticularSolution(a)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::ConvergenceAnalysis(a, _)
            | Self::ForAll(_, a)
            | Self::Exists(_, a) => vec![a.as_ref().clone()],
            | Self::Matrix(m) => {
                m.iter()
                    .flatten()
                    .cloned()
                    .collect()
            },
            | Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v) => v.clone(),
            | Self::Predicate {
                args,
                ..
            } => args.clone(),
            | Self::SparsePolynomial(p) => {
                p.terms
                    .values()
                    .cloned()
                    .collect()
            },
            | Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => {

                vec![
                    integrand
                        .as_ref()
                        .clone(),
                    var.as_ref().clone(),
                    lower_bound
                        .as_ref()
                        .clone(),
                    upper_bound
                        .as_ref()
                        .clone(),
                ]
            },
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                vec![
                    scalar_field
                        .as_ref()
                        .clone(),
                    volume
                        .as_ref()
                        .clone(),
                ]
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                vec![
                    vector_field
                        .as_ref()
                        .clone(),
                    surface
                        .as_ref()
                        .clone(),
                ]
            },
            | Self::DerivativeN(e, _, n) => {

                vec![
                    e.as_ref().clone(),
                    n.as_ref().clone(),
                ]
            },
            | Self::Series(a, _, c, d)
            | Self::Summation(a, _, c, d)
            | Self::Product(a, _, c, d) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            },
            | Self::AsymptoticExpansion(a, _, c, d) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            },
            | Self::Interval(a, b, _, _) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ]
            },
            | Self::Substitute(a, _, c) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                ]
            },
            | Self::Limit(a, _, c) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                ]
            },
            | Self::Ode {
                equation,
                ..
            } => {

                vec![equation
                    .as_ref()
                    .clone()]
            },
            | Self::Pde {
                equation,
                ..
            } => {

                vec![equation
                    .as_ref()
                    .clone()]
            },
            | Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            },
            | Self::ParametricSolution {
                x,
                y,
            } => {

                vec![
                    x.as_ref().clone(),
                    y.as_ref().clone(),
                ]
            },
            | Self::RootOf {
                poly,
                ..
            } => {

                vec![poly
                    .as_ref()
                    .clone()]
            },
            | Self::QuantityWithValue(v, _) => vec![v.as_ref().clone()],
            | Self::CustomArcOne(a) => vec![a.as_ref().clone()],
            | Self::CustomArcTwo(a, b) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ]
            },
            | Self::CustomArcThree(a, b, c) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                ]
            },
            | Self::CustomArcFour(a, b, c, d) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            },
            | Self::CustomArcFive(a, b, c, d, e) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                    e.as_ref().clone(),
                ]
            },
            | Self::CustomVecOne(v) => v.clone(),
            | Self::CustomVecTwo(v1, v2) => {
                v1.iter()
                    .chain(v2.iter())
                    .cloned()
                    .collect()
            },
            | Self::CustomVecThree(v1, v2, v3) => {
                v1.iter()
                    .chain(v2.iter())
                    .chain(v3.iter())
                    .cloned()
                    .collect()
            },
            | Self::CustomVecFour(v1, v2, v3, v4) => {
                v1.iter()
                    .chain(v2.iter())
                    .chain(v3.iter())
                    .chain(v4.iter())
                    .cloned()
                    .collect()
            },
            | Self::CustomVecFive(v1, v2, v3, v4, v5) => {
                v1.iter()
                    .chain(v2.iter())
                    .chain(v3.iter())
                    .chain(v4.iter())
                    .chain(v5.iter())
                    .cloned()
                    .collect()
            },
            | Self::UnaryList(_, a) => vec![a.as_ref().clone()],
            | Self::BinaryList(_, a, b) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ]
            },
            | Self::NaryList(_, v) => v.clone(),
            | _ => vec![],
        }
    }

    #[must_use]

    /// Normalizes the expression by sorting sub-expressions of commutative operators.
    ///
    /// This helps in identifying identical expressions that differ only in terms of operand order.

    pub fn normalize(&self) -> Self {

        match self {
            | Self::Add(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Add(
                    Arc::new(
                        children[0]
                            .clone(),
                    ),
                    Arc::new(
                        children[1]
                            .clone(),
                    ),
                )
            },
            | Self::AddList(list) => {

                let mut children =
                    Vec::new();

                for child in list {

                    if let Self::AddList(sub_list) = child {

                        children.extend(sub_list.clone());
                    } else {

                        children.push(child.clone());
                    }
                }

                children.sort();

                Self::AddList(children)
            },
            | Self::Mul(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Mul(
                    Arc::new(
                        children[0]
                            .clone(),
                    ),
                    Arc::new(
                        children[1]
                            .clone(),
                    ),
                )
            },
            | Self::MulList(list) => {

                let mut children =
                    Vec::new();

                for child in list {

                    if let Self::MulList(sub_list) = child {

                        children.extend(sub_list.clone());
                    } else {

                        children.push(child.clone());
                    }
                }

                children.sort();

                Self::MulList(children)
            },
            | Self::Sub(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Sub(
                    Arc::new(
                        children[0]
                            .clone(),
                    ),
                    Arc::new(
                        children[1]
                            .clone(),
                    ),
                )
            },
            | Self::Div(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Div(
                    Arc::new(
                        children[0]
                            .clone(),
                    ),
                    Arc::new(
                        children[1]
                            .clone(),
                    ),
                )
            },
            | Self::UnaryList(s, a) => {
                Self::UnaryList(
                    s.clone(),
                    Arc::new(
                        a.normalize(),
                    ),
                )
            },
            | Self::BinaryList(
                s,
                a,
                b,
            ) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                if let Some(props) = get_dynamic_op_properties(s) {

                    if props.is_commutative {

                        children.sort();
                    }
                }

                Self::BinaryList(
                    s.clone(),
                    Arc::new(
                        children[0]
                            .clone(),
                    ),
                    Arc::new(
                        children[1]
                            .clone(),
                    ),
                )
            },
            | Self::NaryList(
                s,
                list,
            ) => {

                let mut children =
                    list.clone();

                if let Some(props) = get_dynamic_op_properties(s) {

                    if props.is_commutative {

                        children.sort();
                    }
                }

                Self::NaryList(
                    s.clone(),
                    children,
                )
            },
            | _ => self.clone(),
        }
    }

    /// Converts this expression to its corresponding DAG operation.
    ///
    /// This method extracts the operation from an expression without its children,
    /// mapping it to the appropriate `DagOp` variant for use in the DAG system.
    ///
    /// # Returns
    /// * `Result<DagOp, String>` - The corresponding DAG operation or an error if conversion fails

    pub(crate) fn to_dag_op_internal(
        &self
    ) -> Result<DagOp, String> {

        match self {
            | Self::Constant(c) => {
                Ok(DagOp::Constant(
                    OrderedFloat(*c),
                ))
            },
            | Self::BigInt(i) => {
                Ok(DagOp::BigInt(
                    i.clone(),
                ))
            },
            | Self::Rational(r) => {
                Ok(DagOp::Rational(
                    r.clone(),
                ))
            },
            | Self::Boolean(b) => Ok(DagOp::Boolean(*b)),
            | Self::Variable(s) => {
                Ok(DagOp::Variable(
                    s.clone(),
                ))
            },
            | Self::Pattern(s) => {
                Ok(DagOp::Pattern(
                    s.clone(),
                ))
            },
            | Self::Domain(s) => {
                Ok(DagOp::Domain(
                    s.clone(),
                ))
            },
            | Self::Pi => Ok(DagOp::Pi),
            | Self::E => Ok(DagOp::E),
            | Self::Infinity => Ok(DagOp::Infinity),
            | Self::NegativeInfinity => Ok(DagOp::NegativeInfinity),
            | Self::InfiniteSolutions => Ok(DagOp::InfiniteSolutions),
            | Self::NoSolution => Ok(DagOp::NoSolution),

            | Self::Derivative(_, s) => {
                Ok(DagOp::Derivative(
                    s.clone(),
                ))
            },
            | Self::DerivativeN(_, s, _) => {
                Ok(DagOp::DerivativeN(
                    s.clone(),
                ))
            },
            | Self::Limit(_, s, _) => {
                Ok(DagOp::Limit(
                    s.clone(),
                ))
            },
            | Self::Solve(_, s) => {
                Ok(DagOp::Solve(
                    s.clone(),
                ))
            },
            | Self::ConvergenceAnalysis(_, s) => Ok(DagOp::ConvergenceAnalysis(s.clone())),
            | Self::ForAll(s, _) => {
                Ok(DagOp::ForAll(
                    s.clone(),
                ))
            },
            | Self::Exists(s, _) => {
                Ok(DagOp::Exists(
                    s.clone(),
                ))
            },
            | Self::Substitute(_, s, _) => {
                Ok(DagOp::Substitute(
                    s.clone(),
                ))
            },
            | Self::Ode {
                func,
                var,
                ..
            } => {
                Ok(DagOp::Ode {
                    func : func.clone(),
                    var : var.clone(),
                })
            },
            | Self::Pde {
                func,
                vars,
                ..
            } => {
                Ok(DagOp::Pde {
                    func : func.clone(),
                    vars : vars.clone(),
                })
            },
            | Self::Predicate {
                name,
                ..
            } => {
                Ok(DagOp::Predicate {
                    name : name.clone(),
                })
            },
            | Self::Path(pt, _, _) => {
                Ok(DagOp::Path(
                    pt.clone(),
                ))
            },
            | Self::Interval(_, _, incl_lower, incl_upper) => {
                Ok(DagOp::Interval(
                    *incl_lower,
                    *incl_upper,
                ))
            },
            | Self::RootOf {
                index,
                ..
            } => {
                Ok(DagOp::RootOf {
                    index : *index,
                })
            },
            | Self::SparsePolynomial(p) => Ok(DagOp::SparsePolynomial(p.clone())),
            | Self::QuantityWithValue(_, u) => Ok(DagOp::QuantityWithValue(u.clone())),

            | Self::Add(_, _) => Ok(DagOp::Add),
            | Self::AddList(_) => Ok(DagOp::Add),
            | Self::Sub(_, _) => Ok(DagOp::Sub),
            | Self::Mul(_, _) => Ok(DagOp::Mul),
            | Self::MulList(_) => Ok(DagOp::Mul),
            | Self::Div(_, _) => Ok(DagOp::Div),
            | Self::Neg(_) => Ok(DagOp::Neg),
            | Self::Power(_, _) => Ok(DagOp::Power),
            | Self::Sin(_) => Ok(DagOp::Sin),
            | Self::Cos(_) => Ok(DagOp::Cos),
            | Self::Tan(_) => Ok(DagOp::Tan),
            | Self::Exp(_) => Ok(DagOp::Exp),
            | Self::Log(_) => Ok(DagOp::Log),
            | Self::Abs(_) => Ok(DagOp::Abs),
            | Self::Sqrt(_) => Ok(DagOp::Sqrt),
            | Self::Eq(_, _) => Ok(DagOp::Eq),
            | Self::Lt(_, _) => Ok(DagOp::Lt),
            | Self::Gt(_, _) => Ok(DagOp::Gt),
            | Self::Le(_, _) => Ok(DagOp::Le),
            | Self::Ge(_, _) => Ok(DagOp::Ge),
            | Self::Matrix(m) => {

                let rows = m.len();

                let cols = if rows > 0 {

                    m[0].len()
                } else {

                    0
                };

                Ok(DagOp::Matrix {
                    rows,
                    cols,
                })
            },
            | Self::Vector(_) => Ok(DagOp::Vector),
            | Self::Complex(_, _) => Ok(DagOp::Complex),
            | Self::Transpose(_) => Ok(DagOp::Transpose),
            | Self::MatrixMul(_, _) => Ok(DagOp::MatrixMul),
            | Self::MatrixVecMul(_, _) => Ok(DagOp::MatrixVecMul),
            | Self::Inverse(_) => Ok(DagOp::Inverse),
            | Self::Integral {
                ..
            } => Ok(DagOp::Integral),
            | Self::VolumeIntegral {
                ..
            } => Ok(DagOp::VolumeIntegral),
            | Self::SurfaceIntegral {
                ..
            } => Ok(DagOp::SurfaceIntegral),
            | Self::Sum {
                ..
            } => Ok(DagOp::Sum),
            | Self::Series(_, s, _, _) => {
                Ok(DagOp::Series(
                    s.clone(),
                ))
            },
            | Self::Summation(_, s, _, _) => {
                Ok(DagOp::Summation(
                    s.clone(),
                ))
            },
            | Self::Product(_, s, _, _) => {
                Ok(DagOp::Product(
                    s.clone(),
                ))
            },
            | Self::AsymptoticExpansion(_, s, _, _) => Ok(DagOp::AsymptoticExpansion(s.clone())),
            | Self::Sec(_) => Ok(DagOp::Sec),
            | Self::Csc(_) => Ok(DagOp::Csc),
            | Self::Cot(_) => Ok(DagOp::Cot),
            | Self::ArcSin(_) => Ok(DagOp::ArcSin),
            | Self::ArcCos(_) => Ok(DagOp::ArcCos),
            | Self::ArcTan(_) => Ok(DagOp::ArcTan),
            | Self::ArcSec(_) => Ok(DagOp::ArcSec),
            | Self::ArcCsc(_) => Ok(DagOp::ArcCsc),
            | Self::ArcCot(_) => Ok(DagOp::ArcCot),
            | Self::Sinh(_) => Ok(DagOp::Sinh),
            | Self::Cosh(_) => Ok(DagOp::Cosh),
            | Self::Tanh(_) => Ok(DagOp::Tanh),
            | Self::Sech(_) => Ok(DagOp::Sech),
            | Self::Csch(_) => Ok(DagOp::Csch),
            | Self::Coth(_) => Ok(DagOp::Coth),
            | Self::ArcSinh(_) => Ok(DagOp::ArcSinh),
            | Self::ArcCosh(_) => Ok(DagOp::ArcCosh),
            | Self::ArcTanh(_) => Ok(DagOp::ArcTanh),
            | Self::ArcSech(_) => Ok(DagOp::ArcSech),
            | Self::ArcCsch(_) => Ok(DagOp::ArcCsch),
            | Self::ArcCoth(_) => Ok(DagOp::ArcCoth),
            | Self::LogBase(_, _) => Ok(DagOp::LogBase),
            | Self::Atan2(_, _) => Ok(DagOp::Atan2),
            | Self::Binomial(_, _) => Ok(DagOp::Binomial),
            | Self::Factorial(_) => Ok(DagOp::Factorial),
            | Self::Permutation(_, _) => Ok(DagOp::Permutation),
            | Self::Combination(_, _) => Ok(DagOp::Combination),
            | Self::FallingFactorial(_, _) => Ok(DagOp::FallingFactorial),
            | Self::RisingFactorial(_, _) => Ok(DagOp::RisingFactorial),
            | Self::Boundary(_) => Ok(DagOp::Boundary),
            | Self::Gamma(_) => Ok(DagOp::Gamma),
            | Self::Beta(_, _) => Ok(DagOp::Beta),
            | Self::Erf(_) => Ok(DagOp::Erf),
            | Self::Erfc(_) => Ok(DagOp::Erfc),
            | Self::Erfi(_) => Ok(DagOp::Erfi),
            | Self::Zeta(_) => Ok(DagOp::Zeta),
            | Self::BesselJ(_, _) => Ok(DagOp::BesselJ),
            | Self::BesselY(_, _) => Ok(DagOp::BesselY),
            | Self::LegendreP(_, _) => Ok(DagOp::LegendreP),
            | Self::LaguerreL(_, _) => Ok(DagOp::LaguerreL),
            | Self::HermiteH(_, _) => Ok(DagOp::HermiteH),
            | Self::Digamma(_) => Ok(DagOp::Digamma),
            | Self::KroneckerDelta(_, _) => Ok(DagOp::KroneckerDelta),
            | Self::And(_) => Ok(DagOp::And),
            | Self::Or(_) => Ok(DagOp::Or),
            | Self::Not(_) => Ok(DagOp::Not),
            | Self::Xor(_, _) => Ok(DagOp::Xor),
            | Self::Implies(_, _) => Ok(DagOp::Implies),
            | Self::Equivalent(_, _) => Ok(DagOp::Equivalent),
            | Self::Union(_) => Ok(DagOp::Union),
            | Self::Polynomial(_) => Ok(DagOp::Polynomial),
            | Self::Floor(_) => Ok(DagOp::Floor),
            | Self::IsPrime(_) => Ok(DagOp::IsPrime),
            | Self::Gcd(_, _) => Ok(DagOp::Gcd),
            | Self::Mod(_, _) => Ok(DagOp::Mod),
            | Self::System(_) => Ok(DagOp::System),
            | Self::Solutions(_) => Ok(DagOp::Solutions),
            | Self::ParametricSolution {
                ..
            } => Ok(DagOp::ParametricSolution),
            | Self::GeneralSolution(_) => Ok(DagOp::GeneralSolution),
            | Self::ParticularSolution(_) => Ok(DagOp::ParticularSolution),
            | Self::Fredholm(_, _, _, _) => Ok(DagOp::Fredholm),
            | Self::Volterra(_, _, _, _) => Ok(DagOp::Volterra),
            | Self::Apply(_, _) => Ok(DagOp::Apply),
            | Self::Tuple(_) => Ok(DagOp::Tuple),
            | Self::Distribution(_) => Ok(DagOp::Distribution),
            | Self::Max(_, _) => Ok(DagOp::Max),
            | Self::Quantity(_) => Ok(DagOp::Quantity),
            | Self::Dag(_) => Err("Cannot convert Dag to DagOp".to_string()),

            | Self::CustomZero => Ok(DagOp::CustomZero),
            | Self::CustomString(s) => {
                Ok(DagOp::CustomString(
                    s.clone(),
                ))
            },
            | Self::CustomArcOne(_) => Ok(DagOp::CustomArcOne),
            | Self::CustomArcTwo(_, _) => Ok(DagOp::CustomArcTwo),
            | Self::CustomArcThree(_, _, _) => Ok(DagOp::CustomArcThree),
            | Self::CustomArcFour(_, _, _, _) => Ok(DagOp::CustomArcFour),
            | Self::CustomArcFive(_, _, _, _, _) => Ok(DagOp::CustomArcFive),
            | Self::CustomVecOne(_) => Ok(DagOp::CustomVecOne),
            | Self::CustomVecTwo(_, _) => Ok(DagOp::CustomVecTwo),
            | Self::CustomVecThree(_, _, _) => Ok(DagOp::CustomVecThree),
            | Self::CustomVecFour(_, _, _, _) => Ok(DagOp::CustomVecFour),
            | Self::CustomVecFive(_, _, _, _, _) => Ok(DagOp::CustomVecFive),
            | Self::UnaryList(s, _) => {
                Ok(DagOp::UnaryList(
                    s.clone(),
                ))
            },
            | Self::BinaryList(s, _, _) => {
                Ok(DagOp::BinaryList(
                    s.clone(),
                ))
            },
            | Self::NaryList(s, _) => {
                Ok(DagOp::NaryList(
                    s.clone(),
                ))
            },
        }
    }
}

/// Helper function to extract a numeric float value from an Expr reasonably efficiently.
/// Handles Dag nodes without full conversion to Expr (no recursion) and adds safety
/// checks for extremely large BigInt values to prevent potential panics during conversion.
fn get_numeric_value_efficient(e: &Expr) -> Option<f64> {
    match e {
        Expr::Constant(f) => Some(*f),
        Expr::BigInt(b) => {
            // Safety check for extremely large integers beyond f64 range.
            // A BigInt with > 1024 bits is definitely infinity in f64 (~1.8e308).
            if b.bits() > 1024 {
                return Some(if b.sign() == Sign::Minus {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                });
            }
            b.to_f64()
        }
        Expr::Rational(r) => r.to_f64(),
        Expr::Pi => Some(std::f64::consts::PI),
        Expr::E => Some(std::f64::consts::E),
        Expr::Dag(node) => {
            // Directly inspect the DagOp without recursive to_expr() call.
            match &node.op {
                DagOp::Constant(f) => Some(f.into_inner()),
                DagOp::BigInt(b) => {
                    if b.bits() > 1024 {
                        return Some(if b.sign() == Sign::Minus {
                            f64::NEG_INFINITY
                        } else {
                            f64::INFINITY
                        });
                    }
                    b.to_f64()
                }
                DagOp::Rational(r) => r.to_f64(),
                DagOp::Pi => Some(std::f64::consts::PI),
                DagOp::E => Some(std::f64::consts::E),
                _ => None,
            }
        }
        _ => None,
    }
}
