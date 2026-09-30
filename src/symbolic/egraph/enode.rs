use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use ordered_float::OrderedFloat;

use super::id::Id;
use super::union_find::UnionFind;

/// An equivalence node representing a mathematical operation whose children are E-Class IDs.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ENode {
    // --- Leaf / Terminal Nodes ---
    Constant(OrderedFloat<f64>),
    BigInt(BigInt),
    Rational(BigRational),
    Boolean(bool),
    Variable(String),
    Pattern(String),
    Pi,
    E,
    Infinity,
    NegativeInfinity,

    // --- Basic Arithmetic ---
    Add(Id, Id),
    Sub(Id, Id),
    Mul(Id, Id),
    Div(Id, Id),
    Power(Id, Id),
    Neg(Id),
    Complex(Id, Id),
    AddList(Vec<Id>),
    MulList(Vec<Id>),

    // --- Elementary & Transcendental Functions ---
    Sin(Id),
    Cos(Id),
    Tan(Id),
    Sec(Id),
    Csc(Id),
    Cot(Id),
    ArcSin(Id),
    ArcCos(Id),
    ArcTan(Id),
    ArcSec(Id),
    ArcCsc(Id),
    ArcCot(Id),
    Sinh(Id),
    Cosh(Id),
    Tanh(Id),
    Sech(Id),
    Csch(Id),
    Coth(Id),
    ArcSinh(Id),
    ArcCosh(Id),
    ArcTanh(Id),
    ArcSech(Id),
    ArcCsch(Id),
    ArcCoth(Id),
    Exp(Id),
    Log(Id),
    LogBase(Id, Id),
    Abs(Id),
    Sqrt(Id),
    Atan2(Id, Id),

    // --- High-Order Mathematical Operations (Primary targets for heuristic de-cocooning) ---
    /// First-order derivative: d(body)/d(var)
    Derivative(Id, String),
    /// N-th order derivative: d^n(body)/d(var)^n
    DerivativeN(Id, String, Id),
    /// Definite or Indefinite Integral
    Integral {
        integrand: Id,
        var: Id,
        lower_bound: Id,
        upper_bound: Id,
    },
    /// Limit of body as var approaches target
    Limit(Id, String, Id),
    /// Summation
    Sum {
        body: Id,
        var: Id,
        from: Id,
        to: Id,
    },
    /// Product
    Product(Id, String, Id, Id),
    /// Series expansion
    Series(Id, String, Id, Id),
    /// Algebraic Equation Solver
    Solve(Id, String),
    /// Root of a polynomial
    RootOf(Id, u32),
    /// Ordinary Differential Equation
    Ode {
        equation: Id,
        func: String,
        var: String,
    },
    /// Partial Differential Equation
    Pde {
        equation: Id,
        func: String,
        vars: Vec<String>,
    },
    /// Polynomial factorization request
    Factor(Id),
    /// Simplification with polynomial side-relations (Gröbner basis reduction)
    SimplifyWithRelations {
        expr: Id,
        relations: Vec<Id>,
        vars: Vec<String>,
    },

    // --- Special Functions ---
    Gamma(Id),
    Beta(Id, Id),
    Erf(Id),
    Erfc(Id),
    Erfi(Id),
    Zeta(Id),
    BesselJ(Id, Id),
    BesselY(Id, Id),
    LegendreP(Id, Id),
    LaguerreL(Id, Id),
    HermiteH(Id, Id),
    Digamma(Id),

    // --- Linear Algebra ---
    Matrix(Vec<Vec<Id>>),
    Vector(Vec<Id>),
    Transpose(Id),
    MatrixMul(Id, Id),
    Inverse(Id),
    Determinant(Id),
    Trace(Id),

    // --- Transforms ---
    Laplace { expr: Id, var: String, target: String },
    InverseLaplace { expr: Id, var: String, target: String },
    Fourier { expr: Id, var: String, target: String },
    InverseFourier { expr: Id, var: String, target: String },
    ZTransform { expr: Id, var: String, target: String },
    InverseZTransform { expr: Id, var: String, target: String },

    // --- Vector Calculus ---
    Gradient { expr: Id, vars: Vec<String> },
    Divergence { expr: Id, vars: Vec<String> },
    Curl { expr: Id, vars: Vec<String> },
    Laplacian { expr: Id, vars: Vec<String> },

    // --- Calculus of Variations ---
    EulerLagrange { lagrangian: Id, func: String, var: String },

    // --- Quantum Mechanics ---
    Commutator(Id, Id),
    Anticommutator(Id, Id),

    // --- Complex Analysis ---
    Residue { expr: Id, var: String, point: Id },

    // --- Logic & Relations ---
    Eq(Id, Id),
    Lt(Id, Id),
    Gt(Id, Id),
    Le(Id, Id),
    Ge(Id, Id),
    And(Vec<Id>),
    Or(Vec<Id>),
    Not(Id),
    Xor(Id, Id),
    Implies(Id, Id),
    Equivalent(Id, Id),
    Predicate { name: String, args: Vec<Id> },
    ForAll(String, Id),
    Exists(String, Id),

    // --- General Function Application / Other ---
    Apply(Id, Id),
    Tuple(Vec<Id>),
    Solutions(Vec<Id>),
    InfiniteSolutions,
    NoSolution,

    // --- Number Theory / Arithmetic ---
    Mod(Id, Id),
    Floor(Id),
    Factorial(Id),
    Gcd(Id, Id),
    Lcm(Id, Id),
    Max(Id, Id),
    Binomial(Id, Id),
    Permutation(Id, Id),
    Combination(Id, Id),

    // --- Indefinite Summation ---
    IndefiniteSum { body: Id, var: String, step: Id },
    IndefiniteProduct { body: Id, var: String, step: Id },

    // --- Algebraic Expansion & Number Theory ---
    Expand(Id),
    IsPrime(Id),

    // --- High-Order Algorithmic Oracle Call ---
    OracleCall(String, Vec<Id>),
}

impl ENode {
    /// Returns the slice or vector of child IDs for this node.
    pub fn children(&self) -> Vec<Id> {
        match self {
            Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Pi
            | Self::E
            | Self::Infinity
            | Self::NegativeInfinity
            | Self::InfiniteSolutions
            | Self::NoSolution => Vec::new(),

            Self::Neg(a)
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
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
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::Determinant(a)
            | Self::Trace(a)
            | Self::Laplace { expr: a, .. }
            | Self::InverseLaplace { expr: a, .. }
            | Self::Fourier { expr: a, .. }
            | Self::InverseFourier { expr: a, .. }
            | Self::ZTransform { expr: a, .. }
            | Self::InverseZTransform { expr: a, .. }
            | Self::Gradient { expr: a, .. }
            | Self::Divergence { expr: a, .. }
            | Self::Curl { expr: a, .. }
            | Self::Laplacian { expr: a, .. }
            | Self::EulerLagrange { lagrangian: a, .. }
            | Self::Not(a)
            | Self::ForAll(_, a)
            | Self::Exists(_, a)
            | Self::Factor(a)
            | Self::RootOf(a, _)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::Expand(a)
            | Self::IsPrime(a)
            | Self::Ode { equation: a, .. } => vec![*a],

            Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::MatrixMul(a, b)
            | Self::Commutator(a, b)
            | Self::Anticommutator(a, b)
            | Self::Eq(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Apply(a, b)
            | Self::DerivativeN(a, _, b)
            | Self::Limit(a, _, b)
            | Self::Residue { expr: a, point: b, .. } => vec![*a, *b],

            Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => vec![*integrand, *var, *lower_bound, *upper_bound],

            Self::Sum { body, var, from, to } => vec![*body, *var, *from, *to],

            Self::Product(a, _, b, c) | Self::Series(a, _, b, c) => vec![*a, *b, *c],

            Self::Pde { equation, .. } => vec![*equation],

            Self::SimplifyWithRelations { expr, relations, .. } => {
                let mut c = Vec::with_capacity(relations.len() + 1);
                c.push(*expr);
                c.extend_from_slice(relations);
                c
            }

            Self::AddList(list)
            | Self::MulList(list)
            | Self::Vector(list)
            | Self::And(list)
            | Self::Or(list)
            | Self::Tuple(list)
            | Self::Solutions(list)
            | Self::Predicate { args: list, .. } => list.clone(),

            Self::Matrix(rows) => {
                let mut c = Vec::new();
                for row in rows {
                    c.extend_from_slice(row);
                }
                c
            }

            Self::Mod(a, b)
            | Self::Gcd(a, b)
            | Self::Lcm(a, b)
            | Self::Max(a, b)
            | Self::Binomial(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b) => vec![*a, *b],
            Self::Floor(a) | Self::Factorial(a) => vec![*a],
            Self::IndefiniteSum { body, step, .. }
            | Self::IndefiniteProduct { body, step, .. } => vec![*body, *step],
            Self::OracleCall(_, args) => args.clone(),
        }
    }

    /// Canonicalizes the node by replacing all child IDs with their representative in `uf`.
    pub fn canonicalize(&mut self, uf: &mut UnionFind) {
        match self {
            Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Pi
            | Self::E
            | Self::Infinity
            | Self::NegativeInfinity
            | Self::InfiniteSolutions
            | Self::NoSolution => {}

            Self::Neg(a)
            | Self::Sin(a)
            | Self::Cos(a)
            | Self::Tan(a)
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
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Abs(a)
            | Self::Sqrt(a)
            | Self::Gamma(a)
            | Self::Erf(a)
            | Self::Erfc(a)
            | Self::Erfi(a)
            | Self::Zeta(a)
            | Self::Digamma(a)
            | Self::Transpose(a)
            | Self::Inverse(a)
            | Self::Determinant(a)
            | Self::Trace(a)
            | Self::Laplace { expr: a, .. }
            | Self::InverseLaplace { expr: a, .. }
            | Self::Fourier { expr: a, .. }
            | Self::InverseFourier { expr: a, .. }
            | Self::ZTransform { expr: a, .. }
            | Self::InverseZTransform { expr: a, .. }
            | Self::Gradient { expr: a, .. }
            | Self::Divergence { expr: a, .. }
            | Self::Curl { expr: a, .. }
            | Self::Laplacian { expr: a, .. }
            | Self::EulerLagrange { lagrangian: a, .. }
            | Self::Not(a)
            | Self::ForAll(_, a)
            | Self::Exists(_, a)
            | Self::Factor(a)
            | Self::Floor(a)
            | Self::Factorial(a)
            | Self::RootOf(a, _)
            | Self::Derivative(a, _)
            | Self::Solve(a, _)
            | Self::Expand(a)
            | Self::IsPrime(a)
            | Self::Ode { equation: a, .. } => {
                *a = uf.find(*a);
            }

            Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::Power(a, b)
            | Self::Complex(a, b)
            | Self::LogBase(a, b)
            | Self::Atan2(a, b)
            | Self::Beta(a, b)
            | Self::BesselJ(a, b)
            | Self::BesselY(a, b)
            | Self::LegendreP(a, b)
            | Self::LaguerreL(a, b)
            | Self::HermiteH(a, b)
            | Self::MatrixMul(a, b)
            | Self::Commutator(a, b)
            | Self::Anticommutator(a, b)
            | Self::Eq(a, b)
            | Self::Lt(a, b)
            | Self::Gt(a, b)
            | Self::Le(a, b)
            | Self::Ge(a, b)
            | Self::Xor(a, b)
            | Self::Implies(a, b)
            | Self::Equivalent(a, b)
            | Self::Apply(a, b)
            | Self::Mod(a, b)
            | Self::Gcd(a, b)
            | Self::Lcm(a, b)
            | Self::Max(a, b)
            | Self::Binomial(a, b)
            | Self::Permutation(a, b)
            | Self::Combination(a, b)
            | Self::DerivativeN(a, _, b)
            | Self::Limit(a, _, b)
            | Self::Residue { expr: a, point: b, .. } => {
                *a = uf.find(*a);
                *b = uf.find(*b);
            }

            Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => {
                *integrand = uf.find(*integrand);
                *var = uf.find(*var);
                *lower_bound = uf.find(*lower_bound);
                *upper_bound = uf.find(*upper_bound);
            }

            Self::Sum { body, var, from, to } => {
                *body = uf.find(*body);
                *var = uf.find(*var);
                *from = uf.find(*from);
                *to = uf.find(*to);
            }

            Self::Product(a, _, b, c) | Self::Series(a, _, b, c) => {
                *a = uf.find(*a);
                *b = uf.find(*b);
                *c = uf.find(*c);
            }

            Self::Pde { equation, .. } => {
                *equation = uf.find(*equation);
            }

            Self::SimplifyWithRelations { expr, relations, .. } => {
                *expr = uf.find(*expr);
                for r in relations.iter_mut() {
                    *r = uf.find(*r);
                }
            }

            Self::AddList(list)
            | Self::MulList(list)
            | Self::Vector(list)
            | Self::And(list)
            | Self::Or(list)
            | Self::Tuple(list)
            | Self::Solutions(list)
            | Self::Predicate { args: list, .. }
            | Self::OracleCall(_, list) => {
                for item in list.iter_mut() {
                    *item = uf.find(*item);
                }
            }

            Self::Matrix(rows) => {
                for row in rows.iter_mut() {
                    for item in row.iter_mut() {
                        *item = uf.find(*item);
                    }
                }
            }

            Self::IndefiniteSum { body, step, .. }
            | Self::IndefiniteProduct { body, step, .. } => {
                *body = uf.find(*body);
                *step = uf.find(*step);
            }
        }
    }

    /// Evaluates the intrinsic complexity weight of this operator.
    /// Higher weight means more complex; high-weight operations are prioritized
    /// for de-cocooning (化茧) during heuristic equality saturation.
    #[must_use]
    pub fn complexity_weight(&self) -> u32 {
        match self {
            Self::Constant(_)
            | Self::BigInt(_)
            | Self::Rational(_)
            | Self::Boolean(_)
            | Self::Variable(_)
            | Self::Pattern(_)
            | Self::Pi
            | Self::E => 1,

            Self::Neg(_) | Self::Mul(_, _) => 2,
            Self::Add(_, _) | Self::Sub(_, _) => 3,
            Self::AddList(list) => ((list.len() * 2) as u32).max(3),
            Self::MulList(list) => (list.len() as u32).max(2),

            Self::Div(_, _) | Self::Power(_, _) | Self::Complex(_, _) => 8,

            Self::Sin(_)
            | Self::Cos(_)
            | Self::Tan(_)
            | Self::Sec(_)
            | Self::Csc(_)
            | Self::Cot(_)
            | Self::ArcSin(_)
            | Self::ArcCos(_)
            | Self::ArcTan(_)
            | Self::Sinh(_)
            | Self::Cosh(_)
            | Self::Tanh(_)
            | Self::Exp(_)
            | Self::Log(_)
            | Self::Abs(_)
            | Self::Sqrt(_) => 12,

            Self::ArcSec(_)
            | Self::ArcCsc(_)
            | Self::ArcCot(_)
            | Self::Sech(_)
            | Self::Csch(_)
            | Self::Coth(_)
            | Self::ArcSinh(_)
            | Self::ArcCosh(_)
            | Self::ArcTanh(_)
            | Self::ArcSech(_)
            | Self::ArcCsch(_)
            | Self::ArcCoth(_)
            | Self::LogBase(_, _)
            | Self::Atan2(_, _)
            | Self::Gamma(_)
            | Self::Beta(_, _)
            | Self::Erf(_)
            | Self::Erfc(_)
            | Self::Erfi(_)
            | Self::Zeta(_)
            | Self::BesselJ(_, _)
            | Self::BesselY(_, _)
            | Self::LegendreP(_, _)
            | Self::LaguerreL(_, _)
            | Self::HermiteH(_, _)
            | Self::Digamma(_) => 25,

            Self::Matrix(rows) => {
                let total_elements: usize = rows.iter().map(Vec::len).sum();
                (total_elements as u32).max(10)
            }
            Self::Vector(v) => (v.len() as u32).max(5),
            Self::Transpose(_) | Self::MatrixMul(_, _) => 6_000,
            Self::Inverse(_) => 7_000,
            Self::Determinant(_) => 8_000,
            Self::Trace(_) => 500,

            // --- High-Order Operators (Target of Heuristic De-cocooning) ---
            Self::Derivative(_, _) => 10_000,
            Self::DerivativeN(_, _, _) => 12_000,
            Self::Integral { .. } => 13_000,
            Self::Limit(_, _, _) => 9_000,
            Self::Sum { .. } | Self::Product(_, _, _, _) => 7_000,
            Self::Series(_, _, _, _) => 8_000,
            Self::Factor(_) => 8_000,
            Self::SimplifyWithRelations { .. } => 11_000,
            Self::Solve(_, _) => 15_000,
            Self::RootOf(_, _) => 12_000,
            Self::Ode { .. } => 16_000,
            Self::Pde { .. } => 18_000,
            Self::Laplace { .. }
            | Self::InverseLaplace { .. }
            | Self::Fourier { .. }
            | Self::InverseFourier { .. }
            | Self::ZTransform { .. }
            | Self::InverseZTransform { .. } => 14_000,
            Self::Gradient { .. }
            | Self::Divergence { .. }
            | Self::Curl { .. }
            | Self::Laplacian { .. } => 9_000,
            Self::EulerLagrange { .. } => 15_000,
            Self::Commutator(_, _) | Self::Anticommutator(_, _) => 100,
            Self::Residue { .. } => 9_000,

            Self::Infinity | Self::NegativeInfinity => 2,
            Self::Eq(_, _) | Self::Lt(_, _) | Self::Gt(_, _) | Self::Le(_, _) | Self::Ge(_, _) => 4,
            Self::And(v) | Self::Or(v) => (v.len() as u32).max(4),
            Self::Not(_) => 3,
            Self::Xor(_, _) | Self::Implies(_, _) | Self::Equivalent(_, _) => 8,
            Self::Predicate { args, .. } => (args.len() as u32).max(5),
            Self::ForAll(_, _) | Self::Exists(_, _) => 12,
            Self::Apply(_, _) => 40,
            Self::Tuple(v) | Self::Solutions(v) => (v.len() as u32).max(5),
            Self::InfiniteSolutions | Self::NoSolution => 5,

            Self::Mod(_, _) | Self::Gcd(_, _) | Self::Lcm(_, _) | Self::Max(_, _) => 6,
            Self::Floor(_) | Self::Factorial(_) => 5,
            Self::Binomial(_, _) | Self::Permutation(_, _) | Self::Combination(_, _) => 20,
            Self::IndefiniteSum { .. } | Self::IndefiniteProduct { .. } => 12_000,
            Self::Expand(_) => 200,
            Self::IsPrime(_) => 20,
            Self::OracleCall(..) => 15_000,
        }
    }

    /// Returns whether this node is a high-weight "complex" operator that should be
    /// prioritized for de-cocooning (化茧).
    #[inline(always)]
    #[must_use]
    pub fn is_high_weight_op(&self) -> bool {
        self.complexity_weight() >= 70
    }

    /// Returns whether this node represents numeric zero.
    #[inline]
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match self {
            Self::BigInt(b) => b.is_zero(),
            Self::Rational(r) => r.is_zero(),
            Self::Constant(c) => c.0 == 0.0,
            _ => false,
        }
    }

    /// Returns whether this node represents numeric one.
    #[inline]
    #[must_use]
    pub fn is_one(&self) -> bool {
        match self {
            Self::BigInt(b) => b.is_one(),
            Self::Rational(r) => r.is_one(),
            Self::Constant(c) => c.0 == 1.0,
            _ => false,
        }
    }

    /// Returns whether this node represents numeric -1.
    #[inline]
    #[must_use]
    pub fn is_neg_one(&self) -> bool {
        match self {
            Self::BigInt(b) => *b == BigInt::from(-1),
            Self::Rational(r) => *r == BigRational::from(BigInt::from(-1)),
            Self::Constant(c) => c.0 == -1.0,
            _ => false,
        }
    }
}
