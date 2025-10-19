use std::convert::AsRef;
use std::sync::Arc;

//use std::collections::{BTreeMap, HashMap};
use crate::symbolic::unit_unification::UnitQuantity;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{self, Debug, Write};
use std::hash::{Hash, Hasher};
use std::sync::Mutex;

use lazy_static::lazy_static;

lazy_static! {
    pub static ref DAG_MANAGER: DagManager = DagManager::new();
}

// --- Distribution Trait ---
// Moved here to break circular dependency
pub trait Distribution: Debug + Send + Sync {
    fn pdf(&self, x: &Expr) -> Expr;
    fn cdf(&self, x: &Expr) -> Expr;
    fn expectation(&self) -> Expr;
    fn variance(&self) -> Expr;
    fn mgf(&self, t: &Expr) -> Expr;
    fn clone_box(&self) -> Arc<dyn Distribution>;
}

// --- End Distribution Trait ---

/// `PathType` enum
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
pub enum PathType {
    Line,

    Circle,

    Rectangle,
}

/// Represents a single term in a multivariate polynomial.
///
/// A monomial is a product of variables raised to non-negative integer powers,
/// such as `x^2*y^3`. This struct stores it as a map from variable names (String)
/// to their exponents (u32).
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
pub struct Monomial(pub BTreeMap<String, u32>);

/// Represents a sparse multivariate polynomial.
///
/// A sparse polynomial is stored as a map from `Monomial`s to their `Expr` coefficients.
/// This representation is highly efficient for polynomials with a small number of non-zero
/// terms relative to the degree, such as `x^1000 + 1`.
#[derive(
    Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize, PartialOrd, Ord,
)]
pub struct SparsePolynomial {
    /// The terms of the polynomial, mapping each monomial to its coefficient.
    pub terms: BTreeMap<Monomial, Expr>,
}

/// The central enum representing a mathematical expression in the symbolic system.
///
/// `Expr` is an Abstract Syntax Tree (AST) that can represent a wide variety of
/// mathematical objects and operations. Manual implementations for `Debug`, `Clone`,
/// `PartialEq`, `Eq`, and `Hash` are provided to handle variants containing types
/// that do not derive these traits automatically (e.g., `f64`, `Arc<dyn Distribution>`).
#[derive(serde::Serialize, serde::Deserialize)]
pub enum Expr {
    // --- Basic & Numeric Types ---
    /// A floating-point constant (64-bit).
    Constant(f64),
    /// An arbitrary-precision integer.
    BigInt(BigInt),
    /// An arbitrary-precision rational number.
    Rational(BigRational),
    /// A boolean value (`true` or `false`).
    Boolean(bool),
    /// A symbolic variable, represented by a string name.
    Variable(String),
    /// A pattern variable, used for rule-matching systems.
    Pattern(String),

    // --- Arithmetic Operations ---
    /// Addition of two expressions.
    Add(Arc<Expr>, Arc<Expr>),
    /// Subtraction of two expressions.
    Sub(Arc<Expr>, Arc<Expr>),
    /// Multiplication of two expressions.
    Mul(Arc<Expr>, Arc<Expr>),
    /// Division of two expressions.
    Div(Arc<Expr>, Arc<Expr>),
    /// Exponentiation (`base` raised to the power of `exponent`).
    Power(Arc<Expr>, Arc<Expr>),
    /// Negation of an expression.
    Neg(Arc<Expr>),

    // --- Basic Mathematical Functions ---
    /// The sine function.
    Sin(Arc<Expr>),
    /// The cosine function.
    Cos(Arc<Expr>),
    /// The tangent function.
    Tan(Arc<Expr>),
    /// The natural exponential function, `e^x`.
    Exp(Arc<Expr>),
    /// The natural logarithm, `ln(x)`.
    Log(Arc<Expr>),
    /// The absolute value function, `|x|`.
    Abs(Arc<Expr>),
    /// The square root function.
    Sqrt(Arc<Expr>),

    // --- Equations and Relations ---
    /// Represents an equation (`left = right`).
    Eq(Arc<Expr>, Arc<Expr>),
    /// Less than (`<`).
    Lt(Arc<Expr>, Arc<Expr>),
    /// Greater than (`>`).
    Gt(Arc<Expr>, Arc<Expr>),
    /// Less than or equal to (`<=`).
    Le(Arc<Expr>, Arc<Expr>),
    /// Greater than or equal to (`>=`).
    Ge(Arc<Expr>, Arc<Expr>),

    // --- Linear Algebra ---
    /// A matrix, represented as a vector of row vectors.
    Matrix(Vec<Vec<Expr>>),
    /// A vector (or column matrix).
    Vector(Vec<Expr>),
    /// A complex number with real and imaginary parts.
    Complex(Arc<Expr>, Arc<Expr>),
    /// Matrix transpose.
    Transpose(Arc<Expr>),
    /// Matrix-matrix multiplication.
    MatrixMul(Arc<Expr>, Arc<Expr>),
    /// Matrix-vector multiplication.
    MatrixVecMul(Arc<Expr>, Arc<Expr>),
    /// Matrix inverse.
    Inverse(Arc<Expr>),

    // --- Calculus ---
    /// The derivative of an expression with respect to a variable.
    Derivative(Arc<Expr>, String),
    /// The N-th derivative of an expression.
    DerivativeN(Arc<Expr>, String, Arc<Expr>),
    /// A definite integral of `integrand` with respect to `var` from `lower_bound` to `upper_bound`.
    Integral {
        integrand: Arc<Expr>,
        var: Arc<Expr>,
        lower_bound: Arc<Expr>,
        upper_bound: Arc<Expr>,
    },
    /// A volume integral of a scalar field over a specified volume.
    VolumeIntegral {
        scalar_field: Arc<Expr>,
        volume: Arc<Expr>,
    },
    /// A surface integral of a vector field over a specified surface.
    SurfaceIntegral {
        vector_field: Arc<Expr>,
        surface: Arc<Expr>,
    },
    /// A limit of an expression as a variable approaches a point.
    Limit(Arc<Expr>, String, Arc<Expr>),

    // --- Series and Summations ---
    /// A summation of `body` with `var` from `from` to `to`.
    Sum {
        body: Arc<Expr>,
        var: Arc<Expr>,
        from: Arc<Expr>,
        to: Arc<Expr>,
    },
    /// A finite or infinite series expansion.
    Series(Arc<Expr>, String, Arc<Expr>, Arc<Expr>),
    /// A summation over a range (similar to `Sum`).
    Summation(Arc<Expr>, String, Arc<Expr>, Arc<Expr>),
    /// A product of terms over a range.
    Product(Arc<Expr>, String, Arc<Expr>, Arc<Expr>),
    /// Represents a convergence analysis for a series.
    ConvergenceAnalysis(Arc<Expr>, String),
    /// An asymptotic expansion of a function.
    AsymptoticExpansion(Arc<Expr>, String, Arc<Expr>, Arc<Expr>),

    // --- Trigonometric & Hyperbolic Functions (Extended) ---
    /// Secant function.
    Sec(Arc<Expr>),
    /// Cosecant function.
    Csc(Arc<Expr>),
    /// Cotangent function.
    Cot(Arc<Expr>),
    /// Arcsine (inverse sine).
    ArcSin(Arc<Expr>),
    /// Arccosine (inverse cosine).
    ArcCos(Arc<Expr>),
    /// Arctangent (inverse tangent).
    ArcTan(Arc<Expr>),
    /// Arcsecant (inverse secant).
    ArcSec(Arc<Expr>),
    /// Arccosecant (inverse cosecant).
    ArcCsc(Arc<Expr>),
    /// Arccotangent (inverse cotangent).
    ArcCot(Arc<Expr>),
    /// Hyperbolic sine.
    Sinh(Arc<Expr>),
    /// Hyperbolic cosine.
    Cosh(Arc<Expr>),
    /// Hyperbolic tangent.
    Tanh(Arc<Expr>),
    /// Hyperbolic secant.
    Sech(Arc<Expr>),
    /// Hyperbolic cosecant.
    Csch(Arc<Expr>),
    /// Hyperbolic cotangent.
    Coth(Arc<Expr>),
    /// Inverse hyperbolic sine.
    ArcSinh(Arc<Expr>),
    /// Inverse hyperbolic cosine.
    ArcCosh(Arc<Expr>),
    /// Inverse hyperbolic tangent.
    ArcTanh(Arc<Expr>),
    /// Inverse hyperbolic secant.
    ArcSech(Arc<Expr>),
    /// Inverse hyperbolic cosecant.
    ArcCsch(Arc<Expr>),
    /// Inverse hyperbolic cotangent.
    ArcCoth(Arc<Expr>),
    /// Logarithm with a specified base.
    LogBase(Arc<Expr>, Arc<Expr>),
    /// Two-argument arctangent.
    Atan2(Arc<Expr>, Arc<Expr>),

    // --- Combinatorics ---
    /// Binomial coefficient, "n choose k".
    Binomial(Arc<Expr>, Arc<Expr>),
    /// Factorial, `n!`.
    Factorial(Arc<Expr>),
    /// Permutations, `P(n, k)`.
    Permutation(Arc<Expr>, Arc<Expr>),
    /// Combinations, `C(n, k)`.
    Combination(Arc<Expr>, Arc<Expr>),
    /// Falling factorial.
    FallingFactorial(Arc<Expr>, Arc<Expr>),
    /// Rising factorial.
    RisingFactorial(Arc<Expr>, Arc<Expr>),

    // --- Geometry & Vector Calculus ---
    /// A path for path integrals (e.g., line, circle).
    Path(PathType, Arc<Expr>, Arc<Expr>),
    /// Represents the boundary of a domain.
    Boundary(Arc<Expr>),
    /// Represents a named domain (e.g., for integrals).
    Domain(String),

    // --- Mathematical Constants ---
    /// The mathematical constant Pi (π).
    Pi,
    /// The mathematical constant E (e, Euler's number).
    E,
    /// Represents infinity.
    Infinity,
    /// Represents negative infinity.
    NegativeInfinity,

    // --- Special Functions ---
    /// The Gamma function.
    Gamma(Arc<Expr>),
    /// The Beta function.
    Beta(Arc<Expr>, Arc<Expr>),
    /// The error function.
    Erf(Arc<Expr>),
    /// The complementary error function.
    Erfc(Arc<Expr>),
    /// The imaginary error function.
    Erfi(Arc<Expr>),
    /// The Riemann Zeta function.
    Zeta(Arc<Expr>),
    /// Bessel function of the first kind.
    BesselJ(Arc<Expr>, Arc<Expr>),
    /// Bessel function of the second kind.
    BesselY(Arc<Expr>, Arc<Expr>),
    /// Legendre polynomial.
    LegendreP(Arc<Expr>, Arc<Expr>),
    /// Laguerre polynomial.
    LaguerreL(Arc<Expr>, Arc<Expr>),
    /// Hermite polynomial.
    HermiteH(Arc<Expr>, Arc<Expr>),
    /// The digamma function (psi function).
    Digamma(Arc<Expr>),
    /// The Kronecker delta function.
    KroneckerDelta(Arc<Expr>, Arc<Expr>),

    // --- Logic & Sets ---
    /// Logical AND of a vector of expressions.
    And(Vec<Expr>),
    /// Logical OR of a vector of expressions.
    Or(Vec<Expr>),
    /// Logical NOT.
    Not(Arc<Expr>),
    /// Logical XOR (exclusive OR).
    Xor(Arc<Expr>, Arc<Expr>),
    /// Logical implication (`A => B`).
    Implies(Arc<Expr>, Arc<Expr>),
    /// Logical equivalence (`A <=> B`).
    Equivalent(Arc<Expr>, Arc<Expr>),
    /// A predicate with a name and arguments.
    Predicate {
        name: String,
        args: Vec<Expr>,
    },
    /// Universal quantifier ("for all").
    ForAll(String, Arc<Expr>),
    /// Existential quantifier ("there exists").
    Exists(String, Arc<Expr>),
    /// A union of sets or intervals.
    Union(Vec<Expr>),
    /// An interval with a lower and upper bound, and flags for inclusion.
    Interval(Arc<Expr>, Arc<Expr>, bool, bool),

    // --- Polynomials & Number Theory ---
    /// A dense polynomial represented by its coefficients.
    Polynomial(Vec<Expr>),
    /// A sparse polynomial.
    SparsePolynomial(SparsePolynomial),
    /// The floor function.
    Floor(Arc<Expr>),
    /// A predicate to check if a number is prime.
    IsPrime(Arc<Expr>),
    /// Greatest Common Divisor (GCD).
    Gcd(Arc<Expr>, Arc<Expr>),
    /// Modulo operation.
    Mod(Arc<Expr>, Arc<Expr>),

    // --- Solving & Substitution ---
    /// Represents the action of solving an equation for a variable.
    Solve(Arc<Expr>, String),
    /// Represents the substitution of a variable in an expression with another expression.
    Substitute(Arc<Expr>, String, Arc<Expr>),
    /// Represents a system of equations to be solved.
    System(Vec<Expr>),
    /// Represents the set of solutions to an equation or system.
    Solutions(Vec<Expr>),
    /// A parametric solution, e.g., for a system of ODEs.
    ParametricSolution {
        x: Arc<Expr>,
        y: Arc<Expr>,
    },
    /// Represents the `i`-th root of a polynomial.
    RootOf {
        poly: Arc<Expr>,
        index: u32,
    },
    /// Represents infinite solutions.
    InfiniteSolutions,
    /// Represents that no solution exists.
    NoSolution,

    // --- Differential Equations ---
    /// An ordinary differential equation (ODE).
    Ode {
        equation: Arc<Expr>,
        func: String,
        var: String,
    },
    /// A partial differential equation (PDE).
    Pde {
        equation: Arc<Expr>,
        func: String,
        vars: Vec<String>,
    },
    /// The general solution to a differential equation.
    GeneralSolution(Arc<Expr>),
    /// A particular solution to a differential equation.
    ParticularSolution(Arc<Expr>),

    // --- Integral Equations ---
    /// A Fredholm integral equation.
    Fredholm(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),
    /// A Volterra integral equation.
    Volterra(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),

    // --- Miscellaneous ---
    /// Application of a function to an argument.
    Apply(Arc<Expr>, Arc<Expr>),
    /// A tuple of expressions.
    Tuple(Vec<Expr>),
    /// A node in a Directed Acyclic Graph (DAG) for expression sharing.
    #[serde(skip_serializing, skip_deserializing)]
    Dag(Arc<DagNode>),
    /// A probability distribution.
    #[serde(skip_serializing, skip_deserializing)]
    Distribution(Arc<dyn Distribution>),
    /// Maximum of two expressions.
    Max(Arc<Expr>, Arc<Expr>),
    /// A unified quantity with its value and unit string.
    Quantity(Arc<UnitQuantity>),
    /// A temporary representation of a value with a unit string, before unification.
    QuantityWithValue(Arc<Expr>, String),

    // --- Custom ---
    CustomZero,
    CustomString(String),

    CustomArcOne(Arc<Expr>),
    CustomArcTwo(Arc<Expr>, Arc<Expr>),
    CustomArcThree(Arc<Expr>, Arc<Expr>, Arc<Expr>),
    CustomArcFour(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),
    CustomArcFive(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),

    CustomVecOne(Vec<Expr>),
    CustomVecTwo(Vec<Expr>, Vec<Expr>),
    CustomVecThree(Vec<Expr>, Vec<Expr>, Vec<Expr>),
    CustomVecFour(Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>),
    CustomVecFive(Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>),
}

impl Clone for Expr {
    fn clone(&self) -> Self {
        match self {
            Expr::Constant(c) => Expr::Constant(*c),
            Expr::BigInt(i) => Expr::BigInt(i.clone()),
            Expr::Rational(r) => Expr::Rational(r.clone()),
            Expr::Boolean(b) => Expr::Boolean(*b),
            Expr::Variable(s) => Expr::Variable(s.clone()),
            Expr::Pattern(s) => Expr::Pattern(s.clone()),
            Expr::Add(a, b) => Expr::Add(a.clone(), b.clone()),
            Expr::Sub(a, b) => Expr::Sub(a.clone(), b.clone()),
            Expr::Mul(a, b) => Expr::Mul(a.clone(), b.clone()),
            Expr::Div(a, b) => Expr::Div(a.clone(), b.clone()),
            Expr::Power(a, b) => Expr::Power(a.clone(), b.clone()),
            Expr::Sin(a) => Expr::Sin(a.clone()),
            Expr::Cos(a) => Expr::Cos(a.clone()),
            Expr::Tan(a) => Expr::Tan(a.clone()),
            Expr::Exp(a) => Expr::Exp(a.clone()),
            Expr::Log(a) => Expr::Log(a.clone()),
            Expr::Neg(a) => Expr::Neg(a.clone()),
            Expr::Eq(a, b) => Expr::Eq(a.clone(), b.clone()),
            Expr::Matrix(m) => Expr::Matrix(m.clone()),
            Expr::Vector(v) => Expr::Vector(v.clone()),
            Expr::Complex(re, im) => Expr::Complex(re.clone(), im.clone()),
            Expr::Derivative(e, s) => Expr::Derivative(e.clone(), s.clone()),
            Expr::Sum {
                body,
                var,
                from,
                to,
            } => Expr::Sum {
                body: body.clone(),
                var: var.clone(),
                from: from.clone(),
                to: to.clone(),
            },
            Expr::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => Expr::Integral {
                integrand: integrand.clone(),
                var: var.clone(),
                lower_bound: lower_bound.clone(),
                upper_bound: upper_bound.clone(),
            },
            Expr::Path(pt, p1, p2) => Expr::Path(pt.clone(), p1.clone(), p2.clone()),
            Expr::Abs(a) => Expr::Abs(a.clone()),
            Expr::Sqrt(a) => Expr::Sqrt(a.clone()),
            Expr::Sec(a) => Expr::Sec(a.clone()),
            Expr::Csc(a) => Expr::Csc(a.clone()),
            Expr::Cot(a) => Expr::Cot(a.clone()),
            Expr::ArcSin(a) => Expr::ArcSin(a.clone()),
            Expr::ArcCos(a) => Expr::ArcCos(a.clone()),
            Expr::ArcTan(a) => Expr::ArcTan(a.clone()),
            Expr::ArcSec(a) => Expr::ArcSec(a.clone()),
            Expr::ArcCsc(a) => Expr::ArcCsc(a.clone()),
            Expr::ArcCot(a) => Expr::ArcCot(a.clone()),
            Expr::Sinh(a) => Expr::Sinh(a.clone()),
            Expr::Cosh(a) => Expr::Cosh(a.clone()),
            Expr::Tanh(a) => Expr::Tanh(a.clone()),
            Expr::Sech(a) => Expr::Sech(a.clone()),
            Expr::Csch(a) => Expr::Csch(a.clone()),
            Expr::Coth(a) => Expr::Coth(a.clone()),
            Expr::ArcSinh(a) => Expr::ArcSinh(a.clone()),
            Expr::ArcCosh(a) => Expr::ArcCosh(a.clone()),
            Expr::ArcTanh(a) => Expr::ArcTanh(a.clone()),
            Expr::ArcSech(a) => Expr::ArcSech(a.clone()),
            Expr::ArcCsch(a) => Expr::ArcCsch(a.clone()),
            Expr::ArcCoth(a) => Expr::ArcCoth(a.clone()),
            Expr::LogBase(b, a) => Expr::LogBase(b.clone(), a.clone()),
            Expr::Atan2(y, x) => Expr::Atan2(y.clone(), x.clone()),
            Expr::Binomial(n, k) => Expr::Binomial(n.clone(), k.clone()),
            Expr::Boundary(e) => Expr::Boundary(e.clone()),
            Expr::Domain(s) => Expr::Domain(s.clone()),
            Expr::VolumeIntegral {
                scalar_field,
                volume,
            } => Expr::VolumeIntegral {
                scalar_field: scalar_field.clone(),
                volume: volume.clone(),
            },
            Expr::SurfaceIntegral {
                vector_field,
                surface,
            } => Expr::SurfaceIntegral {
                vector_field: vector_field.clone(),
                surface: surface.clone(),
            },
            Expr::Pi => Expr::Pi,
            Expr::E => Expr::E,
            Expr::Infinity => Expr::Infinity,
            Expr::NegativeInfinity => Expr::NegativeInfinity,
            Expr::Apply(a, b) => Expr::Apply(a.clone(), b.clone()),
            Expr::Tuple(v) => Expr::Tuple(v.clone()),
            Expr::Gamma(a) => Expr::Gamma(a.clone()),
            Expr::Beta(a, b) => Expr::Beta(a.clone(), b.clone()),
            Expr::Erf(a) => Expr::Erf(a.clone()),
            Expr::Erfc(a) => Expr::Erfc(a.clone()),
            Expr::Erfi(a) => Expr::Erfi(a.clone()),
            Expr::Zeta(a) => Expr::Zeta(a.clone()),
            Expr::BesselJ(a, b) => Expr::BesselJ(a.clone(), b.clone()),
            Expr::BesselY(a, b) => Expr::BesselY(a.clone(), b.clone()),
            Expr::LegendreP(a, b) => Expr::LegendreP(a.clone(), b.clone()),
            Expr::LaguerreL(a, b) => Expr::LaguerreL(a.clone(), b.clone()),
            Expr::HermiteH(a, b) => Expr::HermiteH(a.clone(), b.clone()),
            Expr::Digamma(a) => Expr::Digamma(a.clone()),
            Expr::KroneckerDelta(a, b) => Expr::KroneckerDelta(a.clone(), b.clone()),
            Expr::DerivativeN(e, s, n) => Expr::DerivativeN(e.clone(), s.clone(), n.clone()),
            Expr::Series(a, b, c, d) => Expr::Series(a.clone(), b.clone(), c.clone(), d.clone()),
            Expr::Summation(a, b, c, d) => {
                Expr::Summation(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::Product(a, b, c, d) => Expr::Product(a.clone(), b.clone(), c.clone(), d.clone()),
            Expr::ConvergenceAnalysis(e, s) => Expr::ConvergenceAnalysis(e.clone(), s.clone()),
            Expr::AsymptoticExpansion(a, b, c, d) => {
                Expr::AsymptoticExpansion(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::Lt(a, b) => Expr::Lt(a.clone(), b.clone()),
            Expr::Gt(a, b) => Expr::Gt(a.clone(), b.clone()),
            Expr::Le(a, b) => Expr::Le(a.clone(), b.clone()),
            Expr::Ge(a, b) => Expr::Ge(a.clone(), b.clone()),
            Expr::Union(v) => Expr::Union(v.clone()),
            Expr::Interval(a, b, c, d) => Expr::Interval(a.clone(), b.clone(), *c, *d),
            Expr::Solve(e, s) => Expr::Solve(e.clone(), s.clone()),
            Expr::Substitute(a, b, c) => Expr::Substitute(a.clone(), b.clone(), c.clone()),
            Expr::Limit(a, b, c) => Expr::Limit(a.clone(), b.clone(), c.clone()),
            Expr::InfiniteSolutions => Expr::InfiniteSolutions,
            Expr::NoSolution => Expr::NoSolution,
            Expr::Dag(n) => Expr::Dag(n.clone()),
            Expr::Factorial(a) => Expr::Factorial(a.clone()),
            Expr::Permutation(a, b) => Expr::Permutation(a.clone(), b.clone()),
            Expr::Combination(a, b) => Expr::Combination(a.clone(), b.clone()),
            Expr::FallingFactorial(a, b) => Expr::FallingFactorial(a.clone(), b.clone()),
            Expr::RisingFactorial(a, b) => Expr::RisingFactorial(a.clone(), b.clone()),
            Expr::Ode {
                equation,
                func,
                var,
            } => Expr::Ode {
                equation: equation.clone(),
                func: func.clone(),
                var: var.clone(),
            },
            Expr::Pde {
                equation,
                func,
                vars,
            } => Expr::Pde {
                equation: equation.clone(),
                func: func.clone(),
                vars: vars.clone(),
            },
            Expr::GeneralSolution(e) => Expr::GeneralSolution(e.clone()),
            Expr::ParticularSolution(e) => Expr::ParticularSolution(e.clone()),
            Expr::Fredholm(a, b, c, d) => {
                Expr::Fredholm(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::Volterra(a, b, c, d) => {
                Expr::Volterra(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::And(v) => Expr::And(v.clone()),
            Expr::Or(v) => Expr::Or(v.clone()),
            Expr::Not(a) => Expr::Not(a.clone()),
            Expr::Xor(a, b) => Expr::Xor(a.clone(), b.clone()),
            Expr::Implies(a, b) => Expr::Implies(a.clone(), b.clone()),
            Expr::Equivalent(a, b) => Expr::Equivalent(a.clone(), b.clone()),
            Expr::Predicate { name, args } => Expr::Predicate {
                name: name.clone(),
                args: args.clone(),
            },
            Expr::ForAll(s, e) => Expr::ForAll(s.clone(), e.clone()),
            Expr::Exists(s, e) => Expr::Exists(s.clone(), e.clone()),
            Expr::Polynomial(c) => Expr::Polynomial(c.clone()),
            Expr::SparsePolynomial(p) => Expr::SparsePolynomial(p.clone()),
            Expr::Floor(a) => Expr::Floor(a.clone()),
            Expr::IsPrime(a) => Expr::IsPrime(a.clone()),
            Expr::Gcd(a, b) => Expr::Gcd(a.clone(), b.clone()),
            Expr::Distribution(d) => Expr::Distribution(d.clone()),
            Expr::Mod(a, b) => Expr::Mod(a.clone(), b.clone()),
            Expr::Max(a, b) => Expr::Max(a.clone(), b.clone()),
            Expr::Quantity(q) => Expr::Quantity(q.clone()),
            Expr::QuantityWithValue(v, u) => Expr::QuantityWithValue(v.clone(), u.clone()),
            Expr::Transpose(a) => Expr::Transpose(a.clone()),
            Expr::MatrixMul(a, b) => Expr::MatrixMul(a.clone(), b.clone()),
            Expr::MatrixVecMul(a, b) => Expr::MatrixVecMul(a.clone(), b.clone()),
            Expr::Inverse(a) => Expr::Inverse(a.clone()),
            Expr::System(v) => Expr::System(v.clone()),
            Expr::Solutions(v) => Expr::Solutions(v.clone()),
            Expr::ParametricSolution { x, y } => Expr::ParametricSolution {
                x: x.clone(),
                y: y.clone(),
            },
            Expr::RootOf { poly, index } => Expr::RootOf {
                poly: poly.clone(),
                index: *index,
            },

            Expr::CustomZero => Expr::CustomZero,
            Expr::CustomString(a) => Expr::CustomString(a.clone()),
            Expr::CustomArcOne(a) => Expr::CustomArcOne(a.clone()),
            Expr::CustomArcTwo(a, b) => Expr::CustomArcTwo(a.clone(), b.clone()),
            Expr::CustomArcThree(a, b, c) => Expr::CustomArcThree(a.clone(), b.clone(), c.clone()),
            Expr::CustomArcFour(a, b, c, d) => {
                Expr::CustomArcFour(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::CustomArcFive(a, b, c, d, e) => {
                Expr::CustomArcFive(a.clone(), b.clone(), c.clone(), d.clone(), e.clone())
            }
            Expr::CustomVecOne(a) => Expr::CustomVecOne(a.clone()),
            Expr::CustomVecTwo(a, b) => Expr::CustomVecTwo(a.clone(), b.clone()),
            Expr::CustomVecThree(a, b, c) => Expr::CustomVecThree(a.clone(), b.clone(), c.clone()),
            Expr::CustomVecFour(a, b, c, d) => {
                Expr::CustomVecFour(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Expr::CustomVecFive(a, b, c, d, e) => {
                Expr::CustomVecFive(a.clone(), b.clone(), c.clone(), d.clone(), e.clone())
            }
        }
    }
}

impl Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Use Display for a more compact representation in debug outputs
        write!(f, "{}", self)
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Expr::Dag(node) => match node.to_expr() {
                Ok(expr) => write!(f, "{}", expr),
                Err(e) => write!(f, "<Error converting DAG to Expr: {}>", e),
            },
            Expr::Constant(c) => write!(f, "{}", c),
            Expr::BigInt(i) => write!(f, "{}", i),
            Expr::Rational(r) => write!(f, "{}", r),
            Expr::Boolean(b) => write!(f, "{}", b),
            Expr::Variable(s) => write!(f, "{}", s),
            Expr::Pattern(s) => write!(f, "{}", s),
            Expr::Add(a, b) => write!(f, "({} + {})", a, b),
            Expr::Sub(a, b) => write!(f, "({} - {})", a, b),
            Expr::Mul(a, b) => write!(f, "({} * {})", a, b),
            Expr::Div(a, b) => write!(f, "({} / {})", a, b),
            Expr::Power(a, b) => write!(f, "({}^({}))", a, b),
            Expr::Sin(a) => write!(f, "sin({})", a),
            Expr::Cos(a) => write!(f, "cos({})", a),
            Expr::Tan(a) => write!(f, "tan({})", a),
            Expr::Exp(a) => write!(f, "exp({})", a),
            Expr::Log(a) => write!(f, "ln({})", a),
            Expr::Neg(a) => write!(f, "-({})", a),
            Expr::Eq(a, b) => write!(f, "{} = {}", a, b),
            Expr::Matrix(m) => {
                let mut s = String::new();
                s.push('[');
                for (i, row) in m.iter().enumerate() {
                    if i > 0 {
                        s.push_str("; ");
                    }
                    s.push('[');
                    for (j, val) in row.iter().enumerate() {
                        if j > 0 {
                            s.push_str(", ");
                        }
                        s.push_str(&val.to_string());
                    }
                    s.push(']');
                }
                s.push(']');
                write!(f, "{}", s)
            }
            Expr::Vector(v) => {
                let mut s = String::new();
                s.push('[');
                for (i, val) in v.iter().enumerate() {
                    if i > 0 {
                        s.push_str(", ");
                    }
                    s.push_str(&val.to_string());
                }
                s.push(']');
                write!(f, "{}", s)
            }
            Expr::Complex(re, im) => write!(f, "({} + {}i)", re, im),
            Expr::Derivative(expr, var) => write!(f, "d/d{}({})", var, expr),
            Expr::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => write!(
                f,
                "∫ from {} to {} of {} dx",
                lower_bound, upper_bound, integrand
            ),
            Expr::Sum {
                body,
                var,
                from,
                to,
            } => {
                write!(f, "sum({}, {}, {}, {})", body, var, from, to)
            }
            Expr::Path(path_type, p1, p2) => write!(f, "Path({:?}, {}, {})", path_type, p1, p2),
            Expr::Abs(a) => write!(f, "|{}|", a),
            Expr::Sqrt(a) => write!(f, "sqrt({})", a),
            Expr::Sec(a) => write!(f, "sec({})", a),
            Expr::Csc(a) => write!(f, "csc({})", a),
            Expr::Cot(a) => write!(f, "cot({})", a),
            Expr::ArcSin(a) => write!(f, "asin({})", a),
            Expr::ArcCos(a) => write!(f, "acos({})", a),
            Expr::ArcTan(a) => write!(f, "atan({})", a),
            Expr::ArcSec(a) => write!(f, "asec({})", a),
            Expr::ArcCsc(a) => write!(f, "acsc({})", a),
            Expr::ArcCot(a) => write!(f, "acot({})", a),
            Expr::Sinh(a) => write!(f, "sinh({})", a),
            Expr::Cosh(a) => write!(f, "cosh({})", a),
            Expr::Tanh(a) => write!(f, "tanh({})", a),
            Expr::Sech(a) => write!(f, "sech({})", a),
            Expr::Csch(a) => write!(f, "csch({})", a),
            Expr::Coth(a) => write!(f, "coth({})", a),
            Expr::ArcSinh(a) => write!(f, "asinh({})", a),
            Expr::ArcCosh(a) => write!(f, "acosh({})", a),
            Expr::ArcTanh(a) => write!(f, "atanh({})", a),
            Expr::ArcSech(a) => write!(f, "asech({})", a),
            Expr::ArcCsch(a) => write!(f, "acsch({})", a),
            Expr::ArcCoth(a) => write!(f, "acoth({})", a),
            Expr::LogBase(b, a) => write!(f, "log_{}({})", b, a),
            Expr::Atan2(y, x) => write!(f, "atan2({}, {})", y, x),
            Expr::Pi => write!(f, "pi"),
            Expr::E => write!(f, "e"),
            Expr::Infinity => write!(f, "oo"),
            Expr::NegativeInfinity => write!(f, "-oo"),
            Expr::Ode {
                equation,
                func,
                var,
            } => write!(f, "ODE({}, {}, {})", equation, func, var),
            Expr::Pde {
                equation,
                func,
                vars,
            } => write!(f, "PDE({}, {}, {:?})", equation, func, vars),
            Expr::Fredholm(a, b, c, d) => write!(f, "Fredholm({}, {}, {}, {})", a, b, c, d),
            Expr::Volterra(a, b, c, d) => write!(f, "Volterra({}, {}, {}, {})", a, b, c, d),
            Expr::And(v) => write!(
                f,
                "({})",
                v.iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<String>>()
                    .join(" && ")
            ),
            Expr::Or(v) => write!(
                f,
                "({})",
                v.iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<String>>()
                    .join(" || ")
            ),
            Expr::Not(a) => write!(f, "!({})", a),
            Expr::Xor(a, b) => write!(f, "({} ^ {})", a, b),
            Expr::Implies(a, b) => write!(f, "({} => {})", a, b),
            Expr::Equivalent(a, b) => write!(f, "({} <=> {})", a, b),
            Expr::Predicate { name, args } => {
                let args_str = args
                    .iter()
                    .map(|arg| arg.to_string())
                    .collect::<Vec<_>>()
                    .join(", ");
                write!(f, "{}({})", name, args_str)
            }
            Expr::ForAll(s, e) => write!(f, "∀{}. ({})", s, e),
            Expr::Exists(s, e) => write!(f, "∃{}. ({})", s, e),
            Expr::Polynomial(coeffs) => {
                let mut s = String::new();
                for (i, coeff) in coeffs.iter().enumerate().rev() {
                    if !s.is_empty() {
                        s.push_str(" + ");
                    }
                    let _ = write!(s, "{}*x^{}", coeff, i);
                }
                write!(f, "{}", s)
            }
            Expr::SparsePolynomial(p) => {
                let mut s = String::new();
                for (monomial, coeff) in &p.terms {
                    if !s.is_empty() {
                        s.push_str(" + ");
                    }
                    let _ = write!(s, "({})", coeff);
                    for (var, exp) in &monomial.0 {
                        let _ = write!(s, "*{}^{}", var, exp);
                    }
                }
                write!(f, "{}", s)
            }
            Expr::Floor(a) => write!(f, "⌊{}⌋", a),
            Expr::IsPrime(a) => write!(f, "IsPrime({})", a),
            Expr::Gcd(a, b) => write!(f, "gcd({}, {})", a, b),
            Expr::Factorial(a) => write!(f, "{}!", a),
            Expr::Distribution(d) => write!(f, "{:?}", d),
            Expr::Mod(a, b) => write!(f, "({} mod {})", a, b),
            Expr::Max(a, b) => write!(f, "max({}, {})", a, b),
            Expr::System(v) => write!(f, "System({:?})", v),
            Expr::Solutions(v) => write!(f, "Solutions({:?})", v),
            Expr::ParametricSolution { x, y } => write!(f, "(x(t)={}, y(t)={})", x, y),
            Expr::RootOf { poly, index } => write!(f, "RootOf({}, {})", poly, index),
            Expr::Erfc(a) => write!(f, "erfc({})", a),
            Expr::Erfi(a) => write!(f, "erfi({})", a),
            Expr::Zeta(a) => write!(f, "zeta({})", a),

            Expr::CustomZero => write!(f, "CustomZero"),
            Expr::CustomString(s) => write!(f, "CustomString({})", s),
            Expr::CustomArcOne(a) => write!(f, "CustomArcOne({})", a),
            Expr::CustomArcTwo(a, b) => write!(f, "CustomArcTwo({}, {})", a, b),
            Expr::CustomArcThree(a, b, c) => write!(f, "CustomArcThree({}, {}, {})", a, b, c),
            Expr::CustomArcFour(a, b, c, d) => {
                write!(f, "CustomArcFour({}, {}, {}, {})", a, b, c, d)
            }
            Expr::CustomArcFive(a, b, c, d, e) => {
                write!(f, "CustomArcFive({}, {}, {}, {}, {})", a, b, c, d, e)
            }
            Expr::CustomVecOne(v) => write!(f, "CustomVecOne({:?})", v),
            Expr::CustomVecTwo(v1, v2) => write!(f, "CustomVecTwo({:?}, {:?})", v1, v2),
            Expr::CustomVecThree(v1, v2, v3) => {
                write!(f, "CustomVecThree({:?}, {:?}, {:?})", v1, v2, v3)
            }
            Expr::CustomVecFour(v1, v2, v3, v4) => {
                write!(f, "CustomVecFour({:?}, {:?}, {:?}, {:?})", v1, v2, v3, v4)
            }
            Expr::CustomVecFive(v1, v2, v3, v4, v5) => write!(
                f,
                "CustomVecFive({:?}, {:?}, {:?}, {:?}, {:?})",
                v1, v2, v3, v4, v5
            ),

            _ => write!(f, "Unimplemented Display for Expr variant"),
        }
    }
}

impl Expr {
    #[must_use]
    #[inline]
    pub fn re(&self) -> Self {
        if let Expr::Complex(re, _) = self {
            re.as_ref().clone()
        } else {
            self.clone()
        }
    }

    #[must_use]
    #[inline]
    pub fn im(&self) -> Self {
        if let Expr::Complex(_, im) = self {
            im.as_ref().clone()
        } else {
            Expr::Constant(0.0)
        }
    }

    #[inline]
    pub fn to_f64(&self) -> Option<f64> {
        match self {
            Expr::Constant(val) => Some(*val),
            Expr::BigInt(val) => val.to_f64(),
            Expr::Rational(val) => val.to_f64(),
            Expr::Pi => Some(std::f64::consts::PI),
            Expr::E => Some(std::f64::consts::E),
            _ => None,
        }
    }

    /// Returns the operation type of the expression in a unified way.
    pub fn op(&self) -> DagOp {
        match self {
            Expr::Dag(node) => node.op.clone(),
            _ => self.to_dag_op_internal().expect(
                "Failed to convert Expr to DagOp; this should be impossible for any valid Expr",
            ),
        }
    }

    /// Returns the children of the expression in a unified way.
    pub fn children(&self) -> Vec<Expr> {
        match self {
            Expr::Dag(node) => node.children.iter().map(|n| Expr::Dag(n.clone())).collect(),
            _ => self.get_children_internal(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn variant_order(&self) -> i32 {
        match self {
            Expr::Constant(_) => 0,
            Expr::BigInt(_) => 1,
            Expr::Rational(_) => 2,
            Expr::Boolean(_) => 3,
            Expr::Variable(_) => 4,
            Expr::Pattern(_) => 5,
            Expr::Add(_, _) => 6,
            Expr::Sub(_, _) => 7,
            Expr::Mul(_, _) => 8,
            Expr::Div(_, _) => 9,
            Expr::Power(_, _) => 10,
            Expr::Sin(_) => 11,
            Expr::Cos(_) => 12,
            Expr::Tan(_) => 13,
            Expr::Exp(_) => 14,
            Expr::Log(_) => 15,
            Expr::Neg(_) => 16,
            Expr::Eq(_, _) => 17,
            Expr::Matrix(_) => 18,
            Expr::Vector(_) => 19,
            Expr::Complex(_, _) => 20,
            Expr::Derivative(_, _) => 21,
            Expr::Integral { .. } => 22,
            Expr::Sum { .. } => 22, // Assign same order as Integral for now
            Expr::Path(_, _, _) => 23,
            Expr::Abs(_) => 24,
            Expr::Sqrt(_) => 25,
            Expr::Sec(_) => 26,
            Expr::Csc(_) => 27,
            Expr::Cot(_) => 28,
            Expr::ArcSin(_) => 29,
            Expr::ArcCos(_) => 30,
            Expr::ArcTan(_) => 31,
            Expr::ArcSec(_) => 32,
            Expr::ArcCsc(_) => 33,
            Expr::ArcCot(_) => 34,
            Expr::Sinh(_) => 35,
            Expr::Cosh(_) => 36,
            Expr::Tanh(_) => 37,
            Expr::Sech(_) => 38,
            Expr::Csch(_) => 39,
            Expr::Coth(_) => 40,
            Expr::ArcSinh(_) => 41,
            Expr::ArcCosh(_) => 42,
            Expr::ArcTanh(_) => 43,
            Expr::ArcSech(_) => 44,
            Expr::ArcCsch(_) => 45,
            Expr::ArcCoth(_) => 46,
            Expr::LogBase(_, _) => 47,
            Expr::Atan2(_, _) => 48,
            Expr::Binomial(_, _) => 49,
            Expr::Boundary(_) => 50,
            Expr::Domain(_) => 51,
            Expr::VolumeIntegral { .. } => 52,
            Expr::SurfaceIntegral { .. } => 53,
            Expr::Pi => 54,
            Expr::E => 55,
            Expr::Infinity => 56,
            Expr::NegativeInfinity => 57,
            Expr::Apply(_, _) => 58,
            Expr::Tuple(_) => 59,
            Expr::Gamma(_) => 60,
            Expr::Beta(_, _) => 61,
            Expr::Erf(_) => 62,
            Expr::Erfc(_) => 63,
            Expr::Erfi(_) => 64,
            Expr::Zeta(_) => 65,
            Expr::BesselJ(_, _) => 66,
            Expr::BesselY(_, _) => 67,
            Expr::LegendreP(_, _) => 68,
            Expr::LaguerreL(_, _) => 69,
            Expr::HermiteH(_, _) => 70,
            Expr::Digamma(_) => 71,
            Expr::KroneckerDelta(_, _) => 72,
            Expr::DerivativeN(_, _, _) => 73,
            Expr::Series(_, _, _, _) => 74,
            Expr::Summation(_, _, _, _) => 75,
            Expr::Product(_, _, _, _) => 76,
            Expr::ConvergenceAnalysis(_, _) => 77,
            Expr::AsymptoticExpansion(_, _, _, _) => 78,
            Expr::Lt(_, _) => 79,
            Expr::Gt(_, _) => 80,
            Expr::Le(_, _) => 81,
            Expr::Ge(_, _) => 82,
            Expr::Union(_) => 83,
            Expr::Interval(_, _, _, _) => 84,
            Expr::Solve(_, _) => 85,
            Expr::Substitute(_, _, _) => 86,
            Expr::Limit(_, _, _) => 87,
            Expr::InfiniteSolutions => 88,
            Expr::NoSolution => 89,
            Expr::Dag(_) => 90,
            Expr::Factorial(_) => 91,
            Expr::Permutation(_, _) => 92,
            Expr::Combination(_, _) => 93,
            Expr::FallingFactorial(_, _) => 94,
            Expr::RisingFactorial(_, _) => 95,
            Expr::Ode { .. } => 96,
            Expr::Pde { .. } => 97,
            Expr::GeneralSolution(_) => 98,
            Expr::ParticularSolution(_) => 99,
            Expr::Fredholm(_, _, _, _) => 100,
            Expr::Volterra(_, _, _, _) => 101,
            Expr::And(_) => 102,
            Expr::Or(_) => 103,
            Expr::Not(_) => 104,
            Expr::Xor(_, _) => 105,
            Expr::Implies(_, _) => 106,
            Expr::Equivalent(_, _) => 107,
            Expr::Predicate { .. } => 108,
            Expr::ForAll(_, _) => 109,
            Expr::Exists(_, _) => 110,
            Expr::Polynomial(_) => 111,
            Expr::SparsePolynomial(_) => 112,
            Expr::Floor(_) => 113,
            Expr::IsPrime(_) => 114,
            Expr::Gcd(_, _) => 115,
            Expr::Distribution(_) => 116,
            Expr::Mod(_, _) => 117,
            Expr::Max(_, _) => 118,
            Expr::Transpose(_) => 119,
            Expr::MatrixMul(_, _) => 120,
            Expr::MatrixVecMul(_, _) => 121,
            Expr::Inverse(_) => 122,
            Expr::System(_) => 123,
            Expr::Solutions(_) => 124,
            Expr::ParametricSolution { .. } => 125,
            Expr::RootOf { .. } => 126,
            Expr::Quantity(_) => 127,
            Expr::QuantityWithValue(_, _) => 128,
            Expr::CustomZero => 129,
            Expr::CustomString(_) => 130,
            Expr::CustomArcOne(_) => 131,
            Expr::CustomArcTwo(_, _) => 132,
            Expr::CustomArcThree(_, _, _) => 133,
            Expr::CustomArcFour(_, _, _, _) => 134,
            Expr::CustomArcFive(_, _, _, _, _) => 135,
            Expr::CustomVecOne(_) => 136,
            Expr::CustomVecTwo(_, _) => 137,
            Expr::CustomVecThree(_, _, _) => 138,
            Expr::CustomVecFour(_, _, _, _) => 139,
            Expr::CustomVecFive(_, _, _, _, _) => 140,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DagNode {
    pub op: DagOp,
    pub children: Vec<Arc<DagNode>>,
    pub hash: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum DagOp {
    // --- Leaf Nodes ---
    Constant(OrderedFloat<f64>),
    BigInt(BigInt),
    Rational(BigRational),
    Boolean(bool),
    Variable(String),
    Pattern(String),
    Domain(String),
    Pi,
    E,
    Infinity,
    NegativeInfinity,
    InfiniteSolutions,
    NoSolution,

    // --- Operators with associated data ---
    Derivative(String),
    DerivativeN(String),
    Limit(String),
    Solve(String),
    ConvergenceAnalysis(String),
    ForAll(String),
    Exists(String),
    Substitute(String),
    Ode { func: String, var: String },
    Pde { func: String, vars: Vec<String> },
    Predicate { name: String },
    Path(PathType),
    Interval(bool, bool),

    // --- Operators without associated data (children only) ---
    Add,
    Sub,
    Mul,
    Div,
    Neg,
    Power,
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
    Abs,
    Sqrt,
    Eq,
    Lt,
    Gt,
    Le,
    Ge,
    Matrix { rows: usize, cols: usize },
    Vector,
    Complex,
    Transpose,
    MatrixMul,
    MatrixVecMul,
    Inverse,
    Integral,
    VolumeIntegral,
    SurfaceIntegral,
    Sum,
    Series(String),
    Summation(String),
    Product(String),
    AsymptoticExpansion(String),
    Sec,
    Csc,
    Cot,
    ArcSin,
    ArcCos,
    ArcTan,
    ArcSec,
    ArcCsc,
    ArcCot,
    Sinh,
    Cosh,
    Tanh,
    Sech,
    Csch,
    Coth,
    ArcSinh,
    ArcCosh,
    ArcTanh,
    ArcSech,
    ArcCsch,
    ArcCoth,
    LogBase,
    Atan2,
    Binomial,
    Factorial,
    Permutation,
    Combination,
    FallingFactorial,
    RisingFactorial,
    Boundary,
    Gamma,
    Beta,
    Erf,
    Erfc,
    Erfi,
    Zeta,
    BesselJ,
    BesselY,
    LegendreP,
    LaguerreL,
    HermiteH,
    Digamma,
    KroneckerDelta,
    And,
    Or,
    Not,
    Xor,
    Implies,
    Equivalent,
    Union,
    Polynomial,
    SparsePolynomial(SparsePolynomial), // Note: Storing whole struct for simplicity
    Floor,
    IsPrime,
    Gcd,
    Mod,
    System,
    Solutions,
    ParametricSolution,
    RootOf { index: u32 },
    GeneralSolution,
    ParticularSolution,
    Fredholm,
    Volterra,
    Apply,
    Tuple,
    Distribution, // Trait objects are handled separately
    Max,
    Quantity, // Handled separately
    QuantityWithValue(String),

    // --- Custom ---
    CustomZero,
    CustomString(String),
    CustomArcOne,
    CustomArcTwo,
    CustomArcThree,
    CustomArcFour,
    CustomArcFive,
    CustomVecOne,
    CustomVecTwo,
    CustomVecThree,
    CustomVecFour,
    CustomVecFive,
}

impl PartialEq for DagNode {
    fn eq(&self, other: &Self) -> bool {
        self.op == other.op && self.children == other.children
    }
}

impl Eq for DagNode {}

impl Hash for DagNode {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.op.hash(state);
        self.children.hash(state);
    }
}

impl From<DagNode> for Expr {
    fn from(node: DagNode) -> Self {
        node.to_expr().expect("Cannot convert DagNode to Expr.")
    }
}

impl DagNode {
    pub fn to_expr(&self) -> Result<Expr, String> {
        let children_exprs: Result<Vec<Expr>, String> =
            self.children.iter().map(|child| child.to_expr()).collect();
        let children_exprs = children_exprs?;

        macro_rules! arc {
            ($idx:expr) => {
                Arc::new(children_exprs[$idx].clone())
            };
        }

        match &self.op {
            // --- Leaf Nodes ---
            DagOp::Constant(c) => Ok(Expr::Constant(c.into_inner())),
            DagOp::BigInt(i) => Ok(Expr::BigInt(i.clone())),
            DagOp::Rational(r) => Ok(Expr::Rational(r.clone())),
            DagOp::Boolean(b) => Ok(Expr::Boolean(*b)),
            DagOp::Variable(s) => Ok(Expr::Variable(s.clone())),
            DagOp::Pattern(s) => Ok(Expr::Pattern(s.clone())),
            DagOp::Domain(s) => Ok(Expr::Domain(s.clone())),
            DagOp::Pi => Ok(Expr::Pi),
            DagOp::E => Ok(Expr::E),
            DagOp::Infinity => Ok(Expr::Infinity),
            DagOp::NegativeInfinity => Ok(Expr::NegativeInfinity),
            DagOp::InfiniteSolutions => Ok(Expr::InfiniteSolutions),
            DagOp::NoSolution => Ok(Expr::NoSolution),

            // --- Operators with associated data ---
            DagOp::Derivative(s) => Ok(Expr::Derivative(arc!(0), s.clone())),
            DagOp::DerivativeN(s) => Ok(Expr::DerivativeN(arc!(0), s.clone(), arc!(1))),
            DagOp::Limit(s) => Ok(Expr::Limit(arc!(0), s.clone(), arc!(1))),
            DagOp::Solve(s) => Ok(Expr::Solve(arc!(0), s.clone())),
            DagOp::ConvergenceAnalysis(s) => Ok(Expr::ConvergenceAnalysis(arc!(0), s.clone())),
            DagOp::ForAll(s) => Ok(Expr::ForAll(s.clone(), arc!(0))),
            DagOp::Exists(s) => Ok(Expr::Exists(s.clone(), arc!(0))),
            DagOp::Substitute(s) => Ok(Expr::Substitute(arc!(0), s.clone(), arc!(1))),
            DagOp::Ode { func, var } => Ok(Expr::Ode {
                equation: arc!(0),
                func: func.clone(),
                var: var.clone(),
            }),
            DagOp::Pde { func, vars } => Ok(Expr::Pde {
                equation: arc!(0),
                func: func.clone(),
                vars: vars.clone(),
            }),
            DagOp::Predicate { name } => Ok(Expr::Predicate {
                name: name.clone(),
                args: children_exprs,
            }),
            DagOp::Path(pt) => Ok(Expr::Path(pt.clone(), arc!(0), arc!(1))),
            DagOp::Interval(incl_lower, incl_upper) => {
                Ok(Expr::Interval(arc!(0), arc!(1), *incl_lower, *incl_upper))
            }
            DagOp::RootOf { index } => Ok(Expr::RootOf {
                poly: arc!(0),
                index: *index,
            }),
            DagOp::SparsePolynomial(p) => Ok(Expr::SparsePolynomial(p.clone())),
            DagOp::QuantityWithValue(u) => Ok(Expr::QuantityWithValue(arc!(0), u.clone())),

            // --- Operators without associated data (children only) ---
            DagOp::Add => Ok(Expr::Add(arc!(0), arc!(1))),
            DagOp::Sub => Ok(Expr::Sub(arc!(0), arc!(1))),
            DagOp::Mul => Ok(Expr::Mul(arc!(0), arc!(1))),
            DagOp::Div => Ok(Expr::Div(arc!(0), arc!(1))),
            DagOp::Neg => Ok(Expr::Neg(arc!(0))),
            DagOp::Power => Ok(Expr::Power(arc!(0), arc!(1))),
            DagOp::Sin => Ok(Expr::Sin(arc!(0))),
            DagOp::Cos => Ok(Expr::Cos(arc!(0))),
            DagOp::Tan => Ok(Expr::Tan(arc!(0))),
            DagOp::Exp => Ok(Expr::Exp(arc!(0))),
            DagOp::Log => Ok(Expr::Log(arc!(0))),
            DagOp::Abs => Ok(Expr::Abs(arc!(0))),
            DagOp::Sqrt => Ok(Expr::Sqrt(arc!(0))),
            DagOp::Eq => Ok(Expr::Eq(arc!(0), arc!(1))),
            DagOp::Lt => Ok(Expr::Lt(arc!(0), arc!(1))),
            DagOp::Gt => Ok(Expr::Gt(arc!(0), arc!(1))),
            DagOp::Le => Ok(Expr::Le(arc!(0), arc!(1))),
            DagOp::Ge => Ok(Expr::Ge(arc!(0), arc!(1))),
            DagOp::Matrix { rows: _, cols } => {
                let reconstructed_matrix: Vec<Vec<Expr>> = children_exprs
                    .chunks(*cols)
                    .map(|chunk| chunk.to_vec())
                    .collect();
                Ok(Expr::Matrix(reconstructed_matrix))
            }
            DagOp::Vector => Ok(Expr::Vector(children_exprs)),
            DagOp::Complex => Ok(Expr::Complex(arc!(0), arc!(1))),
            DagOp::Transpose => Ok(Expr::Transpose(arc!(0))),
            DagOp::MatrixMul => Ok(Expr::MatrixMul(arc!(0), arc!(1))),
            DagOp::MatrixVecMul => Ok(Expr::MatrixVecMul(arc!(0), arc!(1))),
            DagOp::Inverse => Ok(Expr::Inverse(arc!(0))),
            DagOp::Integral => Ok(Expr::Integral {
                integrand: arc!(0),
                var: arc!(1),
                lower_bound: arc!(2),
                upper_bound: arc!(3),
            }),
            DagOp::VolumeIntegral => Ok(Expr::VolumeIntegral {
                scalar_field: arc!(0),
                volume: arc!(1),
            }),
            DagOp::SurfaceIntegral => Ok(Expr::SurfaceIntegral {
                vector_field: arc!(0),
                surface: arc!(1),
            }),
            DagOp::Sum => Ok(Expr::Sum {
                body: arc!(0),
                var: arc!(1),
                from: arc!(2),
                to: arc!(3),
            }),
            DagOp::Series(s) => Ok(Expr::Series(arc!(0), s.clone(), arc!(1), arc!(2))),
            DagOp::Summation(s) => Ok(Expr::Summation(arc!(0), s.clone(), arc!(1), arc!(2))),
            DagOp::Product(s) => Ok(Expr::Product(arc!(0), s.clone(), arc!(1), arc!(2))),
            DagOp::AsymptoticExpansion(s) => Ok(Expr::AsymptoticExpansion(
                arc!(0),
                s.clone(),
                arc!(1),
                arc!(2),
            )),
            DagOp::Sec => Ok(Expr::Sec(arc!(0))),
            DagOp::Csc => Ok(Expr::Csc(arc!(0))),
            DagOp::Cot => Ok(Expr::Cot(arc!(0))),
            DagOp::ArcSin => Ok(Expr::ArcSin(arc!(0))),
            DagOp::ArcCos => Ok(Expr::ArcCos(arc!(0))),
            DagOp::ArcTan => Ok(Expr::ArcTan(arc!(0))),
            DagOp::ArcSec => Ok(Expr::ArcSec(arc!(0))),
            DagOp::ArcCsc => Ok(Expr::ArcCsc(arc!(0))),
            DagOp::ArcCot => Ok(Expr::ArcCot(arc!(0))),
            DagOp::Sinh => Ok(Expr::Sinh(arc!(0))),
            DagOp::Cosh => Ok(Expr::Cosh(arc!(0))),
            DagOp::Tanh => Ok(Expr::Tanh(arc!(0))),
            DagOp::Sech => Ok(Expr::Sech(arc!(0))),
            DagOp::Csch => Ok(Expr::Csch(arc!(0))),
            DagOp::Coth => Ok(Expr::Coth(arc!(0))),
            DagOp::ArcSinh => Ok(Expr::ArcSinh(arc!(0))),
            DagOp::ArcCosh => Ok(Expr::ArcCosh(arc!(0))),
            DagOp::ArcTanh => Ok(Expr::ArcTanh(arc!(0))),
            DagOp::ArcSech => Ok(Expr::ArcSech(arc!(0))),
            DagOp::ArcCsch => Ok(Expr::ArcCsch(arc!(0))),
            DagOp::ArcCoth => Ok(Expr::ArcCoth(arc!(0))),
            DagOp::LogBase => Ok(Expr::LogBase(arc!(0), arc!(1))),
            DagOp::Atan2 => Ok(Expr::Atan2(arc!(0), arc!(1))),
            DagOp::Binomial => Ok(Expr::Binomial(arc!(0), arc!(1))),
            DagOp::Factorial => Ok(Expr::Factorial(arc!(0))),
            DagOp::Permutation => Ok(Expr::Permutation(arc!(0), arc!(1))),
            DagOp::Combination => Ok(Expr::Combination(arc!(0), arc!(1))),
            DagOp::FallingFactorial => Ok(Expr::FallingFactorial(arc!(0), arc!(1))),
            DagOp::RisingFactorial => Ok(Expr::RisingFactorial(arc!(0), arc!(1))),
            DagOp::Boundary => Ok(Expr::Boundary(arc!(0))),
            DagOp::Gamma => Ok(Expr::Gamma(arc!(0))),
            DagOp::Beta => Ok(Expr::Beta(arc!(0), arc!(1))),
            DagOp::Erf => Ok(Expr::Erf(arc!(0))),
            DagOp::Erfc => Ok(Expr::Erfc(arc!(0))),
            DagOp::Erfi => Ok(Expr::Erfi(arc!(0))),
            DagOp::Zeta => Ok(Expr::Zeta(arc!(0))),
            DagOp::BesselJ => Ok(Expr::BesselJ(arc!(0), arc!(1))),
            DagOp::BesselY => Ok(Expr::BesselY(arc!(0), arc!(1))),
            DagOp::LegendreP => Ok(Expr::LegendreP(arc!(0), arc!(1))),
            DagOp::LaguerreL => Ok(Expr::LaguerreL(arc!(0), arc!(1))),
            DagOp::HermiteH => Ok(Expr::HermiteH(arc!(0), arc!(1))),
            DagOp::Digamma => Ok(Expr::Digamma(arc!(0))),
            DagOp::KroneckerDelta => Ok(Expr::KroneckerDelta(arc!(0), arc!(1))),
            DagOp::And => Ok(Expr::And(children_exprs)),
            DagOp::Or => Ok(Expr::Or(children_exprs)),
            DagOp::Not => Ok(Expr::Not(arc!(0))),
            DagOp::Xor => Ok(Expr::Xor(arc!(0), arc!(1))),
            DagOp::Implies => Ok(Expr::Implies(arc!(0), arc!(1))),
            DagOp::Equivalent => Ok(Expr::Equivalent(arc!(0), arc!(1))),
            DagOp::Union => Ok(Expr::Union(children_exprs)),
            DagOp::Polynomial => Ok(Expr::Polynomial(children_exprs)),
            DagOp::Floor => Ok(Expr::Floor(arc!(0))),
            DagOp::IsPrime => Ok(Expr::IsPrime(arc!(0))),
            DagOp::Gcd => Ok(Expr::Gcd(arc!(0), arc!(1))),
            DagOp::Mod => Ok(Expr::Mod(arc!(0), arc!(1))),
            DagOp::System => Ok(Expr::System(children_exprs)),
            DagOp::Solutions => Ok(Expr::Solutions(children_exprs)),
            DagOp::ParametricSolution => Ok(Expr::ParametricSolution {
                x: arc!(0),
                y: arc!(1),
            }),
            DagOp::GeneralSolution => Ok(Expr::GeneralSolution(arc!(0))),
            DagOp::ParticularSolution => Ok(Expr::ParticularSolution(arc!(0))),
            DagOp::Fredholm => Ok(Expr::Fredholm(arc!(0), arc!(1), arc!(2), arc!(3))),
            DagOp::Volterra => Ok(Expr::Volterra(arc!(0), arc!(1), arc!(2), arc!(3))),
            DagOp::Apply => Ok(Expr::Apply(arc!(0), arc!(1))),
            DagOp::Tuple => Ok(Expr::Tuple(children_exprs)),
            DagOp::Distribution => Ok(Expr::Distribution(children_exprs[0].clone_box_dist()?)),
            DagOp::Max => Ok(Expr::Max(arc!(0), arc!(1))),
            DagOp::Quantity => Ok(Expr::Quantity(children_exprs[0].clone_box_quant()?)),

            // --- Custom ---
            DagOp::CustomZero => Ok(Expr::CustomZero),
            DagOp::CustomString(s) => Ok(Expr::CustomString(s.clone())),
            DagOp::CustomArcOne => Ok(Expr::CustomArcOne(arc!(0))),
            DagOp::CustomArcTwo => Ok(Expr::CustomArcTwo(arc!(0), arc!(1))),
            DagOp::CustomArcThree => Ok(Expr::CustomArcThree(arc!(0), arc!(1), arc!(2))),
            DagOp::CustomArcFour => Ok(Expr::CustomArcFour(arc!(0), arc!(1), arc!(2), arc!(3))),
            DagOp::CustomArcFive => Ok(Expr::CustomArcFive(
                arc!(0),
                arc!(1),
                arc!(2),
                arc!(3),
                arc!(4),
            )),
            DagOp::CustomVecOne => Ok(Expr::CustomVecOne(children_exprs)),
            DagOp::CustomVecTwo => Err("CustomVecTwo to_expr is ambiguous".to_string()),
            DagOp::CustomVecThree => Err("CustomVecThree to_expr is ambiguous".to_string()),
            DagOp::CustomVecFour => Err("CustomVecFour to_expr is ambiguous".to_string()),
            DagOp::CustomVecFive => Err("CustomVecFive to_expr is ambiguous".to_string()),
        }
    }

    pub fn new(op: DagOp, children: Vec<Arc<DagNode>>) -> Arc<Self> {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        op.hash(&mut hasher);
        children.hash(&mut hasher);
        let hash = hasher.finish();
        Arc::new(DagNode { op, children, hash })
    }
}

impl Expr {
    pub fn clone_box_dist(&self) -> Result<Arc<dyn Distribution>, String> {
        if let Expr::Distribution(d) = self {
            Ok(d.clone_box())
        } else {
            Err("Cannot clone into Distribution".to_string())
        }
    }
    pub fn clone_box_quant(&self) -> Result<Arc<UnitQuantity>, String> {
        if let Expr::Quantity(q) = self {
            Ok(q.clone())
        } else {
            Err("Cannot clone into UnitQuantity".to_string())
        }
    }
}

pub struct DagManager {
    nodes: Mutex<HashMap<u64, Vec<Arc<DagNode>>>>,
}

impl Default for DagManager {
    fn default() -> Self {
        Self::new()
    }
}

impl DagManager {
    #[inline]
    pub fn new() -> Self {
        DagManager {
            nodes: Mutex::new(HashMap::new()),
        }
    }

    /// Get an existing node identical to (op, children) if present; otherwise create and insert.
    ///
    /// This implementation avoids returning a node solely based on the `u64` hash.
    /// When a hash bucket is found, we iterate the bucket and compare structural equality
    /// (op + children count + children's hashes). Only when no equal node is found do we insert.
    #[inline]
    pub fn get_or_create_normalized(
        &self,
        op: DagOp,
        mut children: Vec<Arc<DagNode>>,
    ) -> Result<Arc<DagNode>, String> {
        match op {
            DagOp::Add | DagOp::Mul => {
                children.sort_by(|a, b| a.hash.cmp(&b.hash));
            }
            _ => {}
        }

        // Compute 64-bit hash key
        let mut hasher = ahash::AHasher::default();
        op.hash(&mut hasher);
        for c in &children {
            // Use stored hash if present to avoid recursing
            Self::c_hash_for_hasher(c, &mut hasher);
        }
        let hash = hasher.finish();

        // Acquire lock safely: handle PoisonError by recovering the inner guard.
        let mut nodes_guard = match self.nodes.lock() {
            Ok(g) => g,
            Err(pe) => {
                // If a thread panicked previously, recover the poisoned lock's inner data.
                // We prefer to continue with a best-effort recovery instead of panicking.
                pe.into_inner()
            }
        };

        // Ensure the bucket is a vector of candidates to support collision buckets.
        // nodes: HashMap<u64, Vec<Arc<DagNode>>>
        match nodes_guard.entry(hash) {
            Entry::Occupied(mut occ) => {
                // occ.get_mut() is a Vec<Arc<DagNode>>
                let bucket = occ.get_mut();
                // Build a temporary DagNode candidate for structural comparison.
                // We avoid allocating the Arc until we know it's needed.
                for cand in bucket.iter() {
                    if Self::dag_nodes_structurally_equal(cand, &op, &children) {
                        // Found exact structural match; return shared instance.
                        return Ok(cand.clone());
                    }
                }
                // No structural match found in bucket: create new node and push.
                let node = Arc::new(DagNode { op, children, hash });
                bucket.push(node.clone());
                Ok(node)
            }
            Entry::Vacant(vac) => {
                // No bucket yet: create a new vector with the node.
                let node = Arc::new(DagNode { op, children, hash });
                vac.insert(vec![node.clone()]);
                Ok(node)
            }
        }
    }

    /// Helper: compare candidate node with the provided op + children for structural equality.
    /// Uses hash + op + child count + child hashes to decide equality (cheap check),
    /// and falls back to recursive compare of children ops if necessary.
    pub(crate) fn dag_nodes_structurally_equal(
        cand: &Arc<DagNode>,
        op: &DagOp,
        children: &Vec<Arc<DagNode>>,
    ) -> bool {
        // Quick checks: hash, op discrimination, length
        if cand.hash != Self::compute_op_children_hash(op, children) {
            return false;
        }
        if &cand.op != op {
            return false;
        }
        if cand.children.len() != children.len() {
            return false;
        }
        // Compare children's hashes to avoid deep recursion in most cases.
        for (a, b) in cand.children.iter().zip(children.iter()) {
            if a.hash != b.hash {
                return false;
            }
        }
        // All quick checks passed: consider structurally equal.
        true
    }

    /// Compute the same hash that we use as bucket key for an op+children.
    pub(crate) fn compute_op_children_hash(op: &DagOp, children: &Vec<Arc<DagNode>>) -> u64 {
        let mut hasher = ahash::AHasher::default();
        op.hash(&mut hasher);
        for c in children {
            Self::c_hash_for_hasher(c, &mut hasher);
        }
        hasher.finish()
    }

    /// Helper to feed a child's hash into hasher; uses stored hash field if available.
    pub(crate) fn c_hash_for_hasher(c: &Arc<DagNode>, hasher: &mut ahash::AHasher) {
        // Use the child's precomputed hash to avoid deep recursion.
        hasher.write_u64(c.hash);
    }

    #[inline]
    pub fn get_or_create(&self, expr: &Expr) -> Result<Arc<DagNode>, String> {
        if let Expr::Dag(node) = expr {
            return Ok(node.clone());
        }
        let op = expr.to_dag_op_internal()?;
        let children_exprs = expr.get_children_internal();
        let children_nodes = children_exprs
            .iter()
            .map(|child| self.get_or_create(child))
            .collect::<Result<Vec<_>, _>>()?;
        self.get_or_create_normalized(op, children_nodes)
    }
}

impl PartialEq for Expr {
    fn eq(&self, other: &Self) -> bool {
        // Shortcut for identical DAG nodes
        if let (Expr::Dag(n1), Expr::Dag(n2)) = (self, other) {
            if Arc::ptr_eq(n1, n2) {
                return true;
            }
        }

        // Use the unified view
        let self_op = self.op();
        let other_op = other.op();

        if self_op != other_op {
            return false;
        }

        let self_children = self.children();
        let other_children = other.children();

        match self_op {
            DagOp::Add
            | DagOp::Mul
            | DagOp::And
            | DagOp::Or
            | DagOp::Xor
            | DagOp::Equivalent
            | DagOp::Eq => {
                // Eq is commutative
                // Commutative operations: compare children as multisets.
                if self_children.len() != other_children.len() {
                    return false;
                }
                let mut sorted_self = self_children;
                let mut sorted_other = other_children;
                sorted_self.sort();
                sorted_other.sort();
                sorted_self == sorted_other
            }
            _ => {
                // Non-commutative operations: compare children in order.
                self_children == other_children
            }
        }
    }
}

impl Eq for Expr {}

impl Hash for Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Use the unified view
        let op = self.op();
        op.hash(state);

        let mut children = self.children();

        match op {
            DagOp::Add
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
            }
            _ => {
                // Non-commutative: hash children in order.
                for child in children {
                    child.hash(state);
                }
            }
        }
    }
}

impl PartialOrd for Expr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Expr {
    fn cmp(&self, other: &Self) -> Ordering {
        // Shortcut for identical DAG nodes
        if let (Expr::Dag(n1), Expr::Dag(n2)) = (self, other) {
            if Arc::ptr_eq(n1, n2) {
                return Ordering::Equal;
            }
        }

        // Use the unified view. First, compare by operator.
        let op_ordering = self.op().cmp(&other.op());
        if op_ordering != Ordering::Equal {
            return op_ordering;
        }

        // If operators are the same, compare children.
        let mut self_children = self.children();
        let mut other_children = other.children();

        match self.op() {
            DagOp::Add
            | DagOp::Mul
            | DagOp::And
            | DagOp::Or
            | DagOp::Xor
            | DagOp::Equivalent
            | DagOp::Eq => {
                // Commutative: sort children before comparing to establish a canonical order.
                self_children.sort();
                other_children.sort();
                self_children.cmp(&other_children)
            }
            _ => {
                // Non-commutative: compare children in their natural order.
                self_children.cmp(&other_children)
            }
        }
    }
}
#[derive(Debug)]
pub enum SymbolicError {
    Msg(String),
}

impl fmt::Display for SymbolicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SymbolicError::Msg(s) => write!(f, "{}", s),
        }
    }
}

impl From<String> for SymbolicError {
    fn from(s: String) -> Self {
        SymbolicError::Msg(s)
    }
}

impl From<&str> for SymbolicError {
    fn from(s: &str) -> Self {
        SymbolicError::Msg(s.to_string())
    }
}

impl Expr {
    /// Performs a pre-order traversal of the expression tree.
    /// It visits the current node first, then its children.
    pub fn pre_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Expr),
    {
        f(self); // Visit parent
        match self {
            // Binary operators
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b)
            | Expr::Eq(a, b)
            | Expr::Complex(a, b)
            | Expr::LogBase(a, b)
            | Expr::Atan2(a, b)
            | Expr::Binomial(a, b)
            | Expr::Beta(a, b)
            | Expr::BesselJ(a, b)
            | Expr::BesselY(a, b)
            | Expr::LegendreP(a, b)
            | Expr::LaguerreL(a, b)
            | Expr::HermiteH(a, b)
            | Expr::KroneckerDelta(a, b)
            | Expr::Lt(a, b)
            | Expr::Gt(a, b)
            | Expr::Le(a, b)
            | Expr::Ge(a, b)
            | Expr::Permutation(a, b)
            | Expr::Combination(a, b)
            | Expr::FallingFactorial(a, b)
            | Expr::RisingFactorial(a, b)
            | Expr::Xor(a, b)
            | Expr::Implies(a, b)
            | Expr::Equivalent(a, b)
            | Expr::Gcd(a, b)
            | Expr::Mod(a, b)
            | Expr::Max(a, b)
            | Expr::MatrixMul(a, b)
            | Expr::MatrixVecMul(a, b)
            | Expr::Apply(a, b)
            | Expr::Path(_, a, b) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
            }
            // Unary operators
            Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Neg(a)
            | Expr::Abs(a)
            | Expr::Sqrt(a)
            | Expr::Sec(a)
            | Expr::Csc(a)
            | Expr::Cot(a)
            | Expr::ArcSin(a)
            | Expr::ArcCos(a)
            | Expr::ArcTan(a)
            | Expr::ArcSec(a)
            | Expr::ArcCsc(a)
            | Expr::ArcCot(a)
            | Expr::Sinh(a)
            | Expr::Cosh(a)
            | Expr::Tanh(a)
            | Expr::Sech(a)
            | Expr::Csch(a)
            | Expr::Coth(a)
            | Expr::ArcSinh(a)
            | Expr::ArcCosh(a)
            | Expr::ArcTanh(a)
            | Expr::ArcSech(a)
            | Expr::ArcCsch(a)
            | Expr::ArcCoth(a)
            | Expr::Boundary(a)
            | Expr::Gamma(a)
            | Expr::Erf(a)
            | Expr::Erfc(a)
            | Expr::Erfi(a)
            | Expr::Zeta(a)
            | Expr::Digamma(a)
            | Expr::Not(a)
            | Expr::Floor(a)
            | Expr::IsPrime(a)
            | Expr::Factorial(a)
            | Expr::Transpose(a)
            | Expr::Inverse(a)
            | Expr::GeneralSolution(a)
            | Expr::ParticularSolution(a)
            | Expr::Derivative(a, _)
            | Expr::Solve(a, _)
            | Expr::ConvergenceAnalysis(a, _)
            | Expr::ForAll(_, a)
            | Expr::Exists(_, a) => {
                a.pre_order_walk(f);
            }
            // N-ary operators
            Expr::Matrix(m) => m.iter().flatten().for_each(|e| e.pre_order_walk(f)),
            Expr::Vector(v)
            | Expr::Tuple(v)
            | Expr::Polynomial(v)
            | Expr::And(v)
            | Expr::Or(v)
            | Expr::Union(v)
            | Expr::System(v)
            | Expr::Solutions(v) => v.iter().for_each(|e| e.pre_order_walk(f)),
            Expr::Predicate { args, .. } => args.iter().for_each(|e| e.pre_order_walk(f)),
            Expr::SparsePolynomial(p) => p.terms.values().for_each(|c| c.pre_order_walk(f)),
            // More complex operators
            Expr::Sum {
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
            }
            Expr::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {
                integrand.pre_order_walk(f);
                lower_bound.pre_order_walk(f);
                upper_bound.pre_order_walk(f);
            }
            Expr::VolumeIntegral {
                scalar_field,
                volume,
            } => {
                scalar_field.pre_order_walk(f);
                volume.pre_order_walk(f);
            }
            Expr::SurfaceIntegral {
                vector_field,
                surface,
            } => {
                vector_field.pre_order_walk(f);
                surface.pre_order_walk(f);
            }
            Expr::DerivativeN(e, _, n) => {
                e.pre_order_walk(f);
                n.pre_order_walk(f);
            }
            Expr::Series(a, _, c, d) | Expr::Summation(a, _, c, d) | Expr::Product(a, _, c, d) => {
                a.pre_order_walk(f);
                c.pre_order_walk(f);
                d.pre_order_walk(f);
            }
            Expr::AsymptoticExpansion(a, _, c, d) => {
                a.pre_order_walk(f);
                c.pre_order_walk(f);
                d.pre_order_walk(f);
            }
            Expr::Interval(a, b, _, _) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
            }
            Expr::Substitute(a, _, c) => {
                a.pre_order_walk(f);
                c.pre_order_walk(f);
            }
            Expr::Limit(a, _, c) => {
                a.pre_order_walk(f);
                c.pre_order_walk(f);
            }
            Expr::Ode { equation, .. } => equation.pre_order_walk(f),
            Expr::Pde { equation, .. } => equation.pre_order_walk(f),
            Expr::Fredholm(a, b, c, d) | Expr::Volterra(a, b, c, d) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
                c.pre_order_walk(f);
                d.pre_order_walk(f);
            }
            Expr::ParametricSolution { x, y } => {
                x.pre_order_walk(f);
                y.pre_order_walk(f);
            }
            Expr::RootOf { poly, .. } => poly.pre_order_walk(f),
            Expr::QuantityWithValue(v, _) => v.pre_order_walk(f),

            Expr::CustomArcOne(a) => {
                a.pre_order_walk(f);
            }
            Expr::CustomArcTwo(a, b) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
            }
            Expr::CustomArcThree(a, b, c) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
                c.pre_order_walk(f);
            }
            Expr::CustomArcFour(a, b, c, d) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
                c.pre_order_walk(f);
                d.pre_order_walk(f);
            }
            Expr::CustomArcFive(a, b, c, d, e) => {
                a.pre_order_walk(f);
                b.pre_order_walk(f);
                c.pre_order_walk(f);
                d.pre_order_walk(f);
                e.pre_order_walk(f);
            }
            Expr::CustomVecOne(v) => v.iter().for_each(|e| e.pre_order_walk(f)),
            Expr::CustomVecTwo(v1, v2) => {
                for e in v1 {
                    e.pre_order_walk(f);
                }
                for e in v2 {
                    e.pre_order_walk(f);
                }
            }
            Expr::CustomVecThree(v1, v2, v3) => {
                for e in v1 {
                    e.pre_order_walk(f);
                }
                for e in v2 {
                    e.pre_order_walk(f);
                }
                for e in v3 {
                    e.pre_order_walk(f);
                }
            }
            Expr::CustomVecFour(v1, v2, v3, v4) => {
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
            }
            Expr::CustomVecFive(v1, v2, v3, v4, v5) => {
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
            }

            // Leaf nodes
            Expr::Constant(_)
            | Expr::BigInt(_)
            | Expr::Rational(_)
            | Expr::Boolean(_)
            | Expr::Variable(_)
            | Expr::Pattern(_)
            | Expr::Domain(_)
            | Expr::Pi
            | Expr::Quantity(_)
            | Expr::E
            | Expr::Infinity
            | Expr::NegativeInfinity
            | Expr::InfiniteSolutions
            | Expr::NoSolution
            | Expr::Dag(_)
            | Expr::CustomZero
            | Expr::CustomString(_)
            | Expr::Distribution(_) => {}
        }
    }

    /// Performs a post-order traversal of the expression tree.
    /// It visits the children first, then the current node.
    pub fn post_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Expr),
    {
        match self {
            // Binary operators
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b)
            | Expr::Eq(a, b)
            | Expr::Complex(a, b)
            | Expr::LogBase(a, b)
            | Expr::Atan2(a, b)
            | Expr::Binomial(a, b)
            | Expr::Beta(a, b)
            | Expr::BesselJ(a, b)
            | Expr::BesselY(a, b)
            | Expr::LegendreP(a, b)
            | Expr::LaguerreL(a, b)
            | Expr::HermiteH(a, b)
            | Expr::KroneckerDelta(a, b)
            | Expr::Lt(a, b)
            | Expr::Gt(a, b)
            | Expr::Le(a, b)
            | Expr::Ge(a, b)
            | Expr::Permutation(a, b)
            | Expr::Combination(a, b)
            | Expr::FallingFactorial(a, b)
            | Expr::RisingFactorial(a, b)
            | Expr::Xor(a, b)
            | Expr::Implies(a, b)
            | Expr::Equivalent(a, b)
            | Expr::Gcd(a, b)
            | Expr::Mod(a, b)
            | Expr::Max(a, b)
            | Expr::MatrixMul(a, b)
            | Expr::MatrixVecMul(a, b)
            | Expr::Apply(a, b)
            | Expr::Path(_, a, b) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
            }
            // Unary operators
            Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Neg(a)
            | Expr::Abs(a)
            | Expr::Sqrt(a)
            | Expr::Sec(a)
            | Expr::Csc(a)
            | Expr::Cot(a)
            | Expr::ArcSin(a)
            | Expr::ArcCos(a)
            | Expr::ArcTan(a)
            | Expr::ArcSec(a)
            | Expr::ArcCsc(a)
            | Expr::ArcCot(a)
            | Expr::Sinh(a)
            | Expr::Cosh(a)
            | Expr::Tanh(a)
            | Expr::Sech(a)
            | Expr::Csch(a)
            | Expr::Coth(a)
            | Expr::ArcSinh(a)
            | Expr::ArcCosh(a)
            | Expr::ArcTanh(a)
            | Expr::ArcSech(a)
            | Expr::ArcCsch(a)
            | Expr::ArcCoth(a)
            | Expr::Boundary(a)
            | Expr::Gamma(a)
            | Expr::Erf(a)
            | Expr::Erfc(a)
            | Expr::Erfi(a)
            | Expr::Zeta(a)
            | Expr::Digamma(a)
            | Expr::Not(a)
            | Expr::Floor(a)
            | Expr::IsPrime(a)
            | Expr::Factorial(a)
            | Expr::Transpose(a)
            | Expr::Inverse(a)
            | Expr::GeneralSolution(a)
            | Expr::ParticularSolution(a)
            | Expr::Derivative(a, _)
            | Expr::Solve(a, _)
            | Expr::ConvergenceAnalysis(a, _)
            | Expr::ForAll(_, a)
            | Expr::Exists(_, a) => {
                a.post_order_walk(f);
            }
            // N-ary operators
            Expr::Matrix(m) => m.iter().flatten().for_each(|e| e.post_order_walk(f)),
            Expr::Vector(v)
            | Expr::Tuple(v)
            | Expr::Polynomial(v)
            | Expr::And(v)
            | Expr::Or(v)
            | Expr::Union(v)
            | Expr::System(v)
            | Expr::Solutions(v) => v.iter().for_each(|e| e.post_order_walk(f)),
            Expr::Predicate { args, .. } => args.iter().for_each(|e| e.post_order_walk(f)),
            Expr::SparsePolynomial(p) => p.terms.values().for_each(|c| c.post_order_walk(f)),
            // More complex operators
            Expr::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {
                integrand.post_order_walk(f);
                lower_bound.post_order_walk(f);
                upper_bound.post_order_walk(f);
            }
            Expr::Sum {
                body,
                var,
                from,
                to,
            } => {
                body.post_order_walk(f);
                var.post_order_walk(f);
                from.post_order_walk(f);
                to.post_order_walk(f);
            }
            Expr::VolumeIntegral {
                scalar_field,
                volume,
            } => {
                scalar_field.post_order_walk(f);
                volume.post_order_walk(f);
            }
            Expr::SurfaceIntegral {
                vector_field,
                surface,
            } => {
                vector_field.post_order_walk(f);
                surface.post_order_walk(f);
            }
            Expr::DerivativeN(e, _, n) => {
                e.post_order_walk(f);
                n.post_order_walk(f);
            }
            Expr::Series(a, _, c, d) | Expr::Summation(a, _, c, d) | Expr::Product(a, _, c, d) => {
                a.post_order_walk(f);
                c.post_order_walk(f);
                d.post_order_walk(f);
            }
            Expr::AsymptoticExpansion(a, _, c, d) => {
                a.post_order_walk(f);
                c.post_order_walk(f);
                d.post_order_walk(f);
            }
            Expr::Interval(a, b, _, _) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
            }
            Expr::Substitute(a, _, c) => {
                a.post_order_walk(f);
                c.post_order_walk(f);
            }
            Expr::Limit(a, _, c) => {
                a.post_order_walk(f);
                c.post_order_walk(f);
            }
            Expr::Ode { equation, .. } => equation.post_order_walk(f),
            Expr::Pde { equation, .. } => equation.post_order_walk(f),
            Expr::Fredholm(a, b, c, d) | Expr::Volterra(a, b, c, d) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
                c.post_order_walk(f);
                d.post_order_walk(f);
            }
            Expr::ParametricSolution { x, y } => {
                x.post_order_walk(f);
                y.post_order_walk(f);
            }
            Expr::QuantityWithValue(v, _) => v.post_order_walk(f),
            Expr::RootOf { poly, .. } => poly.post_order_walk(f),

            Expr::CustomArcOne(a) => {
                a.post_order_walk(f);
            }
            Expr::CustomArcTwo(a, b) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
            }
            Expr::CustomArcThree(a, b, c) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
                c.post_order_walk(f);
            }
            Expr::CustomArcFour(a, b, c, d) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
                c.post_order_walk(f);
                d.post_order_walk(f);
            }
            Expr::CustomArcFive(a, b, c, d, e) => {
                a.post_order_walk(f);
                b.post_order_walk(f);
                c.post_order_walk(f);
                d.post_order_walk(f);
                e.post_order_walk(f);
            }
            Expr::CustomVecOne(v)
            | Expr::CustomVecTwo(v, _)
            | Expr::CustomVecThree(v, _, _)
            | Expr::CustomVecFour(v, _, _, _)
            | Expr::CustomVecFive(v, _, _, _, _) => {
                for e in v {
                    e.post_order_walk(f);
                }
            }

            // Leaf nodes
            Expr::Constant(_)
            | Expr::BigInt(_)
            | Expr::Rational(_)
            | Expr::Boolean(_)
            | Expr::Variable(_)
            | Expr::Pattern(_)
            | Expr::Domain(_)
            | Expr::Pi
            | Expr::Quantity(_)
            | Expr::E
            | Expr::Infinity
            | Expr::NegativeInfinity
            | Expr::InfiniteSolutions
            | Expr::NoSolution
            | Expr::Dag(_)
            | Expr::CustomZero
            | Expr::CustomString(_)
            | Expr::Distribution(_) => {}
        }
        f(self); // Visit parent
    }

    /// Performs an in-order traversal of the expression tree.
    /// For binary operators, it visits the left child, the node itself, then the right child.
    /// For other nodes, the behavior is adapted as it's not strictly defined.
    pub fn in_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Expr),
    {
        match self {
            // Binary operators
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b)
            | Expr::Eq(a, b)
            | Expr::Complex(a, b)
            | Expr::LogBase(a, b)
            | Expr::Atan2(a, b)
            | Expr::Binomial(a, b)
            | Expr::Beta(a, b)
            | Expr::BesselJ(a, b)
            | Expr::BesselY(a, b)
            | Expr::LegendreP(a, b)
            | Expr::LaguerreL(a, b)
            | Expr::HermiteH(a, b)
            | Expr::KroneckerDelta(a, b)
            | Expr::Lt(a, b)
            | Expr::Gt(a, b)
            | Expr::Le(a, b)
            | Expr::Ge(a, b)
            | Expr::Permutation(a, b)
            | Expr::Combination(a, b)
            | Expr::FallingFactorial(a, b)
            | Expr::RisingFactorial(a, b)
            | Expr::Xor(a, b)
            | Expr::Implies(a, b)
            | Expr::Equivalent(a, b)
            | Expr::Gcd(a, b)
            | Expr::Mod(a, b)
            | Expr::Max(a, b)
            | Expr::MatrixMul(a, b)
            | Expr::MatrixVecMul(a, b)
            | Expr::Apply(a, b)
            | Expr::Path(_, a, b) => {
                a.in_order_walk(f);
                f(self);
                b.in_order_walk(f);
            }
            // Unary operators (treat as pre-order)
            Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Neg(a)
            | Expr::Abs(a)
            | Expr::Sqrt(a)
            | Expr::Sec(a)
            | Expr::Csc(a)
            | Expr::Cot(a)
            | Expr::ArcSin(a)
            | Expr::ArcCos(a)
            | Expr::ArcTan(a)
            | Expr::ArcSec(a)
            | Expr::ArcCsc(a)
            | Expr::ArcCot(a)
            | Expr::Sinh(a)
            | Expr::Cosh(a)
            | Expr::Tanh(a)
            | Expr::Sech(a)
            | Expr::Csch(a)
            | Expr::Coth(a)
            | Expr::ArcSinh(a)
            | Expr::ArcCosh(a)
            | Expr::ArcTanh(a)
            | Expr::ArcSech(a)
            | Expr::ArcCsch(a)
            | Expr::ArcCoth(a)
            | Expr::Boundary(a)
            | Expr::Gamma(a)
            | Expr::Erf(a)
            | Expr::Erfc(a)
            | Expr::Erfi(a)
            | Expr::Zeta(a)
            | Expr::Digamma(a)
            | Expr::Not(a)
            | Expr::Floor(a)
            | Expr::IsPrime(a)
            | Expr::Factorial(a)
            | Expr::Transpose(a)
            | Expr::Inverse(a)
            | Expr::GeneralSolution(a)
            | Expr::ParticularSolution(a)
            | Expr::Derivative(a, _)
            | Expr::Solve(a, _)
            | Expr::ConvergenceAnalysis(a, _)
            | Expr::ForAll(_, a)
            | Expr::Exists(_, a) => {
                f(self);
                a.in_order_walk(f);
            }
            // N-ary operators (visit self, then children)
            Expr::Matrix(m) => {
                f(self);
                m.iter().flatten().for_each(|e| e.in_order_walk(f));
            }
            Expr::Vector(v)
            | Expr::Tuple(v)
            | Expr::Polynomial(v)
            | Expr::And(v)
            | Expr::Or(v)
            | Expr::Union(v)
            | Expr::System(v)
            | Expr::Solutions(v) => {
                f(self);
                for e in v {
                    e.in_order_walk(f);
                }
            }
            Expr::Predicate { args, .. } => {
                f(self);
                for e in args {
                    e.in_order_walk(f);
                }
            }
            Expr::SparsePolynomial(p) => {
                f(self);
                p.terms.values().for_each(|c| c.in_order_walk(f));
            }
            // More complex operators (visit self, then children)
            Expr::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {
                f(self);
                integrand.in_order_walk(f);
                lower_bound.in_order_walk(f);
                upper_bound.in_order_walk(f);
            }
            Expr::Sum {
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
            }
            Expr::VolumeIntegral {
                scalar_field,
                volume,
            } => {
                f(self);
                scalar_field.in_order_walk(f);
                volume.in_order_walk(f);
            }
            Expr::SurfaceIntegral {
                vector_field,
                surface,
            } => {
                f(self);
                vector_field.in_order_walk(f);
                surface.in_order_walk(f);
            }
            Expr::DerivativeN(e, _, n) => {
                f(self);
                e.in_order_walk(f);
                n.in_order_walk(f);
            }
            Expr::Series(a, _, c, d) | Expr::Summation(a, _, c, d) | Expr::Product(a, _, c, d) => {
                f(self);
                a.in_order_walk(f);
                c.in_order_walk(f);
                d.in_order_walk(f);
            }
            Expr::AsymptoticExpansion(a, _, c, _d) => {
                f(self);
                a.in_order_walk(f);
                c.pre_order_walk(f);
            }
            Expr::Interval(a, b, _, _) => {
                f(self);
                a.in_order_walk(f);
                b.in_order_walk(f);
            }
            Expr::Substitute(a, _, c) => {
                f(self);
                a.in_order_walk(f);
                c.in_order_walk(f);
            }
            Expr::Limit(a, _, c) => {
                f(self);
                a.in_order_walk(f);
                c.in_order_walk(f);
            }
            Expr::Ode { equation, .. } => {
                f(self);
                equation.in_order_walk(f);
            }
            Expr::Pde { equation, .. } => {
                f(self);
                equation.in_order_walk(f);
            }
            Expr::Fredholm(a, b, c, d) | Expr::Volterra(a, b, c, d) => {
                f(self);
                a.in_order_walk(f);
                b.in_order_walk(f);
                c.in_order_walk(f);
                d.pre_order_walk(f);
            }
            Expr::ParametricSolution { x, y } => {
                f(self);
                x.in_order_walk(f);
                y.in_order_walk(f);
            }
            Expr::RootOf { poly, .. } => {
                f(self);
                poly.in_order_walk(f);
            }
            Expr::QuantityWithValue(v, _) => v.in_order_walk(f),

            Expr::CustomArcOne(a) => {
                f(self);
                a.in_order_walk(f);
            }
            Expr::CustomArcTwo(a, b) => {
                a.in_order_walk(f);
                f(self);
                b.in_order_walk(f);
            }
            Expr::CustomArcThree(a, b, c) => {
                a.in_order_walk(f);
                b.in_order_walk(f);
                f(self);
                c.in_order_walk(f);
            }
            Expr::CustomArcFour(a, b, c, d) => {
                a.in_order_walk(f);
                b.in_order_walk(f);
                f(self);
                c.in_order_walk(f);
                d.in_order_walk(f);
            }
            Expr::CustomArcFive(a, b, c, d, e) => {
                a.in_order_walk(f);
                b.in_order_walk(f);
                f(self);
                c.in_order_walk(f);
                d.in_order_walk(f);
                e.in_order_walk(f);
            }
            Expr::CustomVecOne(v) => {
                f(self);
                for e in v {
                    e.in_order_walk(f);
                }
            }
            Expr::CustomVecTwo(v1, v2) => {
                f(self);
                for e in v1 {
                    e.in_order_walk(f);
                }
                for e in v2 {
                    e.in_order_walk(f);
                }
            }
            Expr::CustomVecThree(v1, v2, v3) => {
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
            }
            Expr::CustomVecFour(v1, v2, v3, v4) => {
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
            }
            Expr::CustomVecFive(v1, v2, v3, v4, v5) => {
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
            }

            // Leaf nodes
            Expr::Constant(_)
            | Expr::BigInt(_)
            | Expr::Rational(_)
            | Expr::Boolean(_)
            | Expr::Variable(_)
            | Expr::Pattern(_)
            | Expr::Domain(_)
            | Expr::Pi
            | Expr::Quantity(_)
            | Expr::E
            | Expr::Infinity
            | Expr::NegativeInfinity
            | Expr::InfiniteSolutions
            | Expr::NoSolution
            | Expr::Dag(_)
            | Expr::CustomZero
            | Expr::CustomString(_)
            | Expr::Distribution(_) => {}
        }
        f(self); // Visit parent
    }

    pub(crate) fn get_children_internal(&self) -> Vec<Expr> {
        match self {
            Expr::Add(a, b)
            | Expr::Sub(a, b)
            | Expr::Mul(a, b)
            | Expr::Div(a, b)
            | Expr::Power(a, b)
            | Expr::Eq(a, b)
            | Expr::Lt(a, b)
            | Expr::Gt(a, b)
            | Expr::Le(a, b)
            | Expr::Ge(a, b)
            | Expr::Complex(a, b)
            | Expr::LogBase(a, b)
            | Expr::Atan2(a, b)
            | Expr::Binomial(a, b)
            | Expr::Beta(a, b)
            | Expr::BesselJ(a, b)
            | Expr::BesselY(a, b)
            | Expr::LegendreP(a, b)
            | Expr::LaguerreL(a, b)
            | Expr::HermiteH(a, b)
            | Expr::KroneckerDelta(a, b)
            | Expr::Permutation(a, b)
            | Expr::Combination(a, b)
            | Expr::FallingFactorial(a, b)
            | Expr::RisingFactorial(a, b)
            | Expr::Xor(a, b)
            | Expr::Implies(a, b)
            | Expr::Equivalent(a, b)
            | Expr::Gcd(a, b)
            | Expr::Mod(a, b)
            | Expr::Max(a, b)
            | Expr::MatrixMul(a, b)
            | Expr::MatrixVecMul(a, b)
            | Expr::Apply(a, b) => vec![a.as_ref().clone(), b.as_ref().clone()],
            Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Neg(a)
            | Expr::Abs(a)
            | Expr::Sqrt(a)
            | Expr::Sec(a)
            | Expr::Csc(a)
            | Expr::Cot(a)
            | Expr::ArcSin(a)
            | Expr::ArcCos(a)
            | Expr::ArcTan(a)
            | Expr::ArcSec(a)
            | Expr::ArcCsc(a)
            | Expr::ArcCot(a)
            | Expr::Sinh(a)
            | Expr::Cosh(a)
            | Expr::Tanh(a)
            | Expr::Sech(a)
            | Expr::Csch(a)
            | Expr::Coth(a)
            | Expr::ArcSinh(a)
            | Expr::ArcCosh(a)
            | Expr::ArcTanh(a)
            | Expr::ArcSech(a)
            | Expr::ArcCsch(a)
            | Expr::ArcCoth(a)
            | Expr::Boundary(a)
            | Expr::Gamma(a)
            | Expr::Erf(a)
            | Expr::Erfc(a)
            | Expr::Erfi(a)
            | Expr::Zeta(a)
            | Expr::Digamma(a)
            | Expr::Not(a)
            | Expr::Floor(a)
            | Expr::IsPrime(a)
            | Expr::Factorial(a)
            | Expr::Transpose(a)
            | Expr::Inverse(a)
            | Expr::GeneralSolution(a)
            | Expr::ParticularSolution(a) => vec![a.as_ref().clone()],
            Expr::Matrix(m) => m.iter().flatten().cloned().collect(),
            Expr::Vector(v)
            | Expr::Tuple(v)
            | Expr::Polynomial(v)
            | Expr::And(v)
            | Expr::Or(v)
            | Expr::Union(v)
            | Expr::System(v)
            | Expr::Solutions(v) => v.clone(),
            Expr::Predicate { args, .. } => args.clone(),
            Expr::SparsePolynomial(p) => p.terms.values().cloned().collect(),
            Expr::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => vec![
                integrand.as_ref().clone(),
                var.as_ref().clone(),
                lower_bound.as_ref().clone(),
                upper_bound.as_ref().clone(),
            ],
            Expr::VolumeIntegral {
                scalar_field,
                volume,
            } => vec![scalar_field.as_ref().clone(), volume.as_ref().clone()],
            Expr::SurfaceIntegral {
                vector_field,
                surface,
            } => vec![vector_field.as_ref().clone(), surface.as_ref().clone()],
            Expr::DerivativeN(e, _, n) => vec![e.as_ref().clone(), n.as_ref().clone()],
            Expr::Series(a, _, c, d) | Expr::Summation(a, _, c, d) | Expr::Product(a, _, c, d) => {
                vec![a.as_ref().clone(), c.as_ref().clone(), d.as_ref().clone()]
            }
            Expr::AsymptoticExpansion(a, _, c, d) => {
                vec![a.as_ref().clone(), c.as_ref().clone(), d.as_ref().clone()]
            }
            Expr::Interval(a, b, _, _) => vec![a.as_ref().clone(), b.as_ref().clone()],
            Expr::Substitute(a, _, c) => vec![a.as_ref().clone(), c.as_ref().clone()],
            Expr::Limit(a, _, c) => vec![a.as_ref().clone(), c.as_ref().clone()],
            Expr::Ode { equation, .. } => vec![equation.as_ref().clone()],
            Expr::Pde { equation, .. } => vec![equation.as_ref().clone()],
            Expr::Fredholm(a, b, c, d) | Expr::Volterra(a, b, c, d) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
            ],
            Expr::ParametricSolution { x, y } => vec![x.as_ref().clone(), y.as_ref().clone()],
            Expr::RootOf { poly, .. } => vec![poly.as_ref().clone()],
            Expr::QuantityWithValue(v, _) => vec![v.as_ref().clone()],
            Expr::CustomArcOne(a) => vec![a.as_ref().clone()],
            Expr::CustomArcTwo(a, b) => vec![a.as_ref().clone(), b.as_ref().clone()],
            Expr::CustomArcThree(a, b, c) => {
                vec![a.as_ref().clone(), b.as_ref().clone(), c.as_ref().clone()]
            }
            Expr::CustomArcFour(a, b, c, d) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
            ],
            Expr::CustomArcFive(a, b, c, d, e) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
                e.as_ref().clone(),
            ],
            Expr::CustomVecOne(v) => v.clone(),
            Expr::CustomVecTwo(v1, v2) => v1.iter().chain(v2.iter()).cloned().collect(),
            Expr::CustomVecThree(v1, v2, v3) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .cloned()
                .collect(),
            Expr::CustomVecFour(v1, v2, v3, v4) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .chain(v4.iter())
                .cloned()
                .collect(),
            Expr::CustomVecFive(v1, v2, v3, v4, v5) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .chain(v4.iter())
                .chain(v5.iter())
                .cloned()
                .collect(),
            _ => vec![],
        }
    }

    #[must_use]
    pub fn normalize(&self) -> Expr {
        match self {
            Expr::Add(a, b) => {
                let mut children = [a.as_ref().clone(), b.as_ref().clone()];
                children.sort();
                Expr::Add(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Expr::Mul(a, b) => {
                let mut children = [a.as_ref().clone(), b.as_ref().clone()];
                children.sort();
                Expr::Mul(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Expr::Sub(a, b) => {
                let mut children = [a.as_ref().clone(), b.as_ref().clone()];
                children.sort();
                Expr::Sub(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Expr::Div(a, b) => {
                let mut children = [a.as_ref().clone(), b.as_ref().clone()];
                children.sort();
                Expr::Div(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            _ => self.clone(),
        }
    }

    pub(crate) fn to_dag_op_internal(&self) -> Result<DagOp, String> {
        match self {
            Expr::Constant(c) => Ok(DagOp::Constant(OrderedFloat(*c))),
            Expr::BigInt(i) => Ok(DagOp::BigInt(i.clone())),
            Expr::Rational(r) => Ok(DagOp::Rational(r.clone())),
            Expr::Boolean(b) => Ok(DagOp::Boolean(*b)),
            Expr::Variable(s) => Ok(DagOp::Variable(s.clone())),
            Expr::Pattern(s) => Ok(DagOp::Pattern(s.clone())),
            Expr::Domain(s) => Ok(DagOp::Domain(s.clone())),
            Expr::Pi => Ok(DagOp::Pi),
            Expr::E => Ok(DagOp::E),
            Expr::Infinity => Ok(DagOp::Infinity),
            Expr::NegativeInfinity => Ok(DagOp::NegativeInfinity),
            Expr::InfiniteSolutions => Ok(DagOp::InfiniteSolutions),
            Expr::NoSolution => Ok(DagOp::NoSolution),

            Expr::Derivative(_, s) => Ok(DagOp::Derivative(s.clone())),
            Expr::DerivativeN(_, s, _) => Ok(DagOp::DerivativeN(s.clone())),
            Expr::Limit(_, s, _) => Ok(DagOp::Limit(s.clone())),
            Expr::Solve(_, s) => Ok(DagOp::Solve(s.clone())),
            Expr::ConvergenceAnalysis(_, s) => Ok(DagOp::ConvergenceAnalysis(s.clone())),
            Expr::ForAll(s, _) => Ok(DagOp::ForAll(s.clone())),
            Expr::Exists(s, _) => Ok(DagOp::Exists(s.clone())),
            Expr::Substitute(_, s, _) => Ok(DagOp::Substitute(s.clone())),
            Expr::Ode { func, var, .. } => Ok(DagOp::Ode {
                func: func.clone(),
                var: var.clone(),
            }),
            Expr::Pde { func, vars, .. } => Ok(DagOp::Pde {
                func: func.clone(),
                vars: vars.clone(),
            }),
            Expr::Predicate { name, .. } => Ok(DagOp::Predicate { name: name.clone() }),
            Expr::Path(pt, _, _) => Ok(DagOp::Path(pt.clone())),
            Expr::Interval(_, _, incl_lower, incl_upper) => {
                Ok(DagOp::Interval(*incl_lower, *incl_upper))
            }
            Expr::RootOf { index, .. } => Ok(DagOp::RootOf { index: *index }),
            Expr::SparsePolynomial(p) => Ok(DagOp::SparsePolynomial(p.clone())),
            Expr::QuantityWithValue(_, u) => Ok(DagOp::QuantityWithValue(u.clone())),

            Expr::Add(_, _) => Ok(DagOp::Add),
            Expr::Sub(_, _) => Ok(DagOp::Sub),
            Expr::Mul(_, _) => Ok(DagOp::Mul),
            Expr::Div(_, _) => Ok(DagOp::Div),
            Expr::Neg(_) => Ok(DagOp::Neg),
            Expr::Power(_, _) => Ok(DagOp::Power),
            Expr::Sin(_) => Ok(DagOp::Sin),
            Expr::Cos(_) => Ok(DagOp::Cos),
            Expr::Tan(_) => Ok(DagOp::Tan),
            Expr::Exp(_) => Ok(DagOp::Exp),
            Expr::Log(_) => Ok(DagOp::Log),
            Expr::Abs(_) => Ok(DagOp::Abs),
            Expr::Sqrt(_) => Ok(DagOp::Sqrt),
            Expr::Eq(_, _) => Ok(DagOp::Eq),
            Expr::Lt(_, _) => Ok(DagOp::Lt),
            Expr::Gt(_, _) => Ok(DagOp::Gt),
            Expr::Le(_, _) => Ok(DagOp::Le),
            Expr::Ge(_, _) => Ok(DagOp::Ge),
            Expr::Matrix(m) => {
                let rows = m.len();
                let cols = if rows > 0 { m[0].len() } else { 0 };
                Ok(DagOp::Matrix { rows, cols })
            }
            Expr::Vector(_) => Ok(DagOp::Vector),
            Expr::Complex(_, _) => Ok(DagOp::Complex),
            Expr::Transpose(_) => Ok(DagOp::Transpose),
            Expr::MatrixMul(_, _) => Ok(DagOp::MatrixMul),
            Expr::MatrixVecMul(_, _) => Ok(DagOp::MatrixVecMul),
            Expr::Inverse(_) => Ok(DagOp::Inverse),
            Expr::Integral { .. } => Ok(DagOp::Integral),
            Expr::VolumeIntegral { .. } => Ok(DagOp::VolumeIntegral),
            Expr::SurfaceIntegral { .. } => Ok(DagOp::SurfaceIntegral),
            Expr::Sum { .. } => Ok(DagOp::Sum),
            Expr::Series(_, s, _, _) => Ok(DagOp::Series(s.clone())),
            Expr::Summation(_, s, _, _) => Ok(DagOp::Summation(s.clone())),
            Expr::Product(_, s, _, _) => Ok(DagOp::Product(s.clone())),
            Expr::AsymptoticExpansion(_, s, _, _) => Ok(DagOp::AsymptoticExpansion(s.clone())),
            Expr::Sec(_) => Ok(DagOp::Sec),
            Expr::Csc(_) => Ok(DagOp::Csc),
            Expr::Cot(_) => Ok(DagOp::Cot),
            Expr::ArcSin(_) => Ok(DagOp::ArcSin),
            Expr::ArcCos(_) => Ok(DagOp::ArcCos),
            Expr::ArcTan(_) => Ok(DagOp::ArcTan),
            Expr::ArcSec(_) => Ok(DagOp::ArcSec),
            Expr::ArcCsc(_) => Ok(DagOp::ArcCsc),
            Expr::ArcCot(_) => Ok(DagOp::ArcCot),
            Expr::Sinh(_) => Ok(DagOp::Sinh),
            Expr::Cosh(_) => Ok(DagOp::Cosh),
            Expr::Tanh(_) => Ok(DagOp::Tanh),
            Expr::Sech(_) => Ok(DagOp::Sech),
            Expr::Csch(_) => Ok(DagOp::Csch),
            Expr::Coth(_) => Ok(DagOp::Coth),
            Expr::ArcSinh(_) => Ok(DagOp::ArcSinh),
            Expr::ArcCosh(_) => Ok(DagOp::ArcCosh),
            Expr::ArcTanh(_) => Ok(DagOp::ArcTanh),
            Expr::ArcSech(_) => Ok(DagOp::ArcSech),
            Expr::ArcCsch(_) => Ok(DagOp::ArcCsch),
            Expr::ArcCoth(_) => Ok(DagOp::ArcCoth),
            Expr::LogBase(_, _) => Ok(DagOp::LogBase),
            Expr::Atan2(_, _) => Ok(DagOp::Atan2),
            Expr::Binomial(_, _) => Ok(DagOp::Binomial),
            Expr::Factorial(_) => Ok(DagOp::Factorial),
            Expr::Permutation(_, _) => Ok(DagOp::Permutation),
            Expr::Combination(_, _) => Ok(DagOp::Combination),
            Expr::FallingFactorial(_, _) => Ok(DagOp::FallingFactorial),
            Expr::RisingFactorial(_, _) => Ok(DagOp::RisingFactorial),
            Expr::Boundary(_) => Ok(DagOp::Boundary),
            Expr::Gamma(_) => Ok(DagOp::Gamma),
            Expr::Beta(_, _) => Ok(DagOp::Beta),
            Expr::Erf(_) => Ok(DagOp::Erf),
            Expr::Erfc(_) => Ok(DagOp::Erfc),
            Expr::Erfi(_) => Ok(DagOp::Erfi),
            Expr::Zeta(_) => Ok(DagOp::Zeta),
            Expr::BesselJ(_, _) => Ok(DagOp::BesselJ),
            Expr::BesselY(_, _) => Ok(DagOp::BesselY),
            Expr::LegendreP(_, _) => Ok(DagOp::LegendreP),
            Expr::LaguerreL(_, _) => Ok(DagOp::LaguerreL),
            Expr::HermiteH(_, _) => Ok(DagOp::HermiteH),
            Expr::Digamma(_) => Ok(DagOp::Digamma),
            Expr::KroneckerDelta(_, _) => Ok(DagOp::KroneckerDelta),
            Expr::And(_) => Ok(DagOp::And),
            Expr::Or(_) => Ok(DagOp::Or),
            Expr::Not(_) => Ok(DagOp::Not),
            Expr::Xor(_, _) => Ok(DagOp::Xor),
            Expr::Implies(_, _) => Ok(DagOp::Implies),
            Expr::Equivalent(_, _) => Ok(DagOp::Equivalent),
            Expr::Union(_) => Ok(DagOp::Union),
            Expr::Polynomial(_) => Ok(DagOp::Polynomial),
            Expr::Floor(_) => Ok(DagOp::Floor),
            Expr::IsPrime(_) => Ok(DagOp::IsPrime),
            Expr::Gcd(_, _) => Ok(DagOp::Gcd),
            Expr::Mod(_, _) => Ok(DagOp::Mod),
            Expr::System(_) => Ok(DagOp::System),
            Expr::Solutions(_) => Ok(DagOp::Solutions),
            Expr::ParametricSolution { .. } => Ok(DagOp::ParametricSolution),
            Expr::GeneralSolution(_) => Ok(DagOp::GeneralSolution),
            Expr::ParticularSolution(_) => Ok(DagOp::ParticularSolution),
            Expr::Fredholm(_, _, _, _) => Ok(DagOp::Fredholm),
            Expr::Volterra(_, _, _, _) => Ok(DagOp::Volterra),
            Expr::Apply(_, _) => Ok(DagOp::Apply),
            Expr::Tuple(_) => Ok(DagOp::Tuple),
            Expr::Distribution(_) => Ok(DagOp::Distribution),
            Expr::Max(_, _) => Ok(DagOp::Max),
            Expr::Quantity(_) => Ok(DagOp::Quantity),
            Expr::Dag(_) => Err("Cannot convert Dag to DagOp".to_string()),

            Expr::CustomZero => Ok(DagOp::CustomZero),
            Expr::CustomString(s) => Ok(DagOp::CustomString(s.clone())),
            Expr::CustomArcOne(_) => Ok(DagOp::CustomArcOne),
            Expr::CustomArcTwo(_, _) => Ok(DagOp::CustomArcTwo),
            Expr::CustomArcThree(_, _, _) => Ok(DagOp::CustomArcThree),
            Expr::CustomArcFour(_, _, _, _) => Ok(DagOp::CustomArcFour),
            Expr::CustomArcFive(_, _, _, _, _) => Ok(DagOp::CustomArcFive),
            Expr::CustomVecOne(_) => Ok(DagOp::CustomVecOne),
            Expr::CustomVecTwo(_, _) => Ok(DagOp::CustomVecTwo),
            Expr::CustomVecThree(_, _, _) => Ok(DagOp::CustomVecThree),
            Expr::CustomVecFour(_, _, _, _) => Ok(DagOp::CustomVecFour),
            Expr::CustomVecFive(_, _, _, _, _) => Ok(DagOp::CustomVecFive),
        }
    }
}

impl AsRef<Expr> for Expr {
    fn as_ref(&self) -> &Expr {
        self
    }
}

// --- Helper Macros ---
macro_rules! unary_constructor {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        pub fn $name<A>(a: A) -> Expr
        where
            A: AsRef<Expr>,
        {
            let dag_a = DAG_MANAGER.get_or_create(a.as_ref()).unwrap();
            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a])
                .unwrap();
            Expr::Dag(node)
        }
    };
}
macro_rules! binary_constructor {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        pub fn $name<A, B>(a: A, b: B) -> Expr
        where
            A: AsRef<Expr>,
            B: AsRef<Expr>,
        {
            let dag_a = DAG_MANAGER.get_or_create(a.as_ref()).unwrap();
            let dag_b = DAG_MANAGER.get_or_create(b.as_ref()).unwrap();
            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a, dag_b])
                .unwrap();
            Expr::Dag(node)
        }
    };
}
macro_rules! n_ary_constructor {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        pub fn $name<I, T>(elements: I) -> Expr
        where
            I: IntoIterator<Item = T>,
            T: AsRef<Expr>,
        {
            let children_nodes = elements
                .into_iter()
                .map(|child| DAG_MANAGER.get_or_create(child.as_ref()).unwrap())
                .collect::<Vec<_>>();
            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, children_nodes)
                .unwrap();
            Expr::Dag(node)
        }
    };
}

// --- Smart Constructors ---
impl Expr {
    // --- Leaf Node Constructors ---
    /// Creates a new Constant expression, managed by the DAG.
    pub fn new_constant(c: f64) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Constant(OrderedFloat(c)), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new Variable expression, managed by the DAG.
    pub fn new_variable(name: &str) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Variable(name.to_string()), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new BigInt expression, managed by the DAG.
    pub fn new_bigint(i: BigInt) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::BigInt(i), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new Rational expression, managed by the DAG.
    pub fn new_rational(r: BigRational) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Rational(r), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }
    /// Creates a new Pi expression, managed by the DAG.
    pub fn new_pi() -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Pi, vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new E expression, managed by the DAG.
    pub fn new_e() -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::E, vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new Infinity expression, managed by the DAG.
    pub fn new_infinity() -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Infinity, vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    /// Creates a new NegativeInfinity expression, managed by the DAG.
    pub fn new_negative_infinity() -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::NegativeInfinity, vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    // --- Unary Operator Constructors ---
    unary_constructor!(new_sin, Sin);
    unary_constructor!(new_cos, Cos);
    unary_constructor!(new_tan, Tan);
    unary_constructor!(new_exp, Exp);
    unary_constructor!(new_log, Log);
    unary_constructor!(new_neg, Neg);
    unary_constructor!(new_abs, Abs);
    unary_constructor!(new_sqrt, Sqrt);
    unary_constructor!(new_transpose, Transpose);
    unary_constructor!(new_inverse, Inverse);
    unary_constructor!(new_sec, Sec);
    unary_constructor!(new_csc, Csc);
    unary_constructor!(new_cot, Cot);
    unary_constructor!(new_arcsin, ArcSin);
    unary_constructor!(new_arccos, ArcCos);
    unary_constructor!(new_arctan, ArcTan);
    unary_constructor!(new_arcsec, ArcSec);
    unary_constructor!(new_arccsc, ArcCsc);
    unary_constructor!(new_arccot, ArcCot);
    unary_constructor!(new_sinh, Sinh);
    unary_constructor!(new_cosh, Cosh);
    unary_constructor!(new_tanh, Tanh);
    unary_constructor!(new_sech, Sech);
    unary_constructor!(new_csch, Csch);
    unary_constructor!(new_coth, Coth);
    unary_constructor!(new_arcsinh, ArcSinh);
    unary_constructor!(new_arccosh, ArcCosh);
    unary_constructor!(new_arctanh, ArcTanh);
    unary_constructor!(new_arcsech, ArcSech);
    unary_constructor!(new_arccsch, ArcCsch);
    unary_constructor!(new_arccoth, ArcCoth);
    unary_constructor!(new_not, Not);
    unary_constructor!(new_floor, Floor);
    unary_constructor!(new_gamma, Gamma);
    unary_constructor!(new_erf, Erf);
    unary_constructor!(new_erfc, Erfc);
    unary_constructor!(new_erfi, Erfi);
    unary_constructor!(new_zeta, Zeta);
    unary_constructor!(new_digamma, Digamma);

    // --- Binary Operator Constructors ---
    binary_constructor!(new_add, Add);
    binary_constructor!(new_sub, Sub);
    binary_constructor!(new_mul, Mul);
    binary_constructor!(new_div, Div);
    binary_constructor!(new_pow, Power);
    binary_constructor!(new_complex, Complex);
    binary_constructor!(new_matrix_mul, MatrixMul);
    binary_constructor!(new_matrix_vec_mul, MatrixVecMul);
    binary_constructor!(new_log_base, LogBase);
    binary_constructor!(new_atan2, Atan2);
    binary_constructor!(new_xor, Xor);
    binary_constructor!(new_implies, Implies);
    binary_constructor!(new_equivalent, Equivalent);
    binary_constructor!(new_beta, Beta);
    binary_constructor!(new_bessel_j, BesselJ);
    binary_constructor!(new_bessel_y, BesselY);
    binary_constructor!(new_legendre_p, LegendreP);
    binary_constructor!(new_laguerre_l, LaguerreL);
    binary_constructor!(new_hermite_h, HermiteH);
    binary_constructor!(new_kronecker_delta, KroneckerDelta);
    binary_constructor!(new_apply, Apply);

    // --- N-ary Constructors ---
    n_ary_constructor!(new_vector, Vector);
    n_ary_constructor!(new_and, And);
    n_ary_constructor!(new_or, Or);
    n_ary_constructor!(new_union, Union);
    n_ary_constructor!(new_polynomial, Polynomial);
    n_ary_constructor!(new_tuple, Tuple);

    // --- Special Constructors ---
    /// Creates a new Matrix expression, managed by the DAG.
    pub fn new_matrix<I, J, T>(elements: I) -> Expr
    where
        I: IntoIterator<Item = J>,
        J: IntoIterator<Item = T>,
        T: AsRef<Expr>,
    {
        let mut flat_children_nodes = Vec::new();
        let mut rows = 0;
        let mut cols = 0;

        for row_iter in elements {
            rows += 1;
            let mut current_cols = 0;
            for element in row_iter {
                let node = DAG_MANAGER
                    .get_or_create(element.as_ref())
                    .expect("Value is valid");
                flat_children_nodes.push(node);
                current_cols += 1;
            }

            if cols == 0 {
                cols = current_cols;
            } else if current_cols != cols {
                panic!("Matrix rows must have consistent length");
            }
        }

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Matrix { rows, cols }, flat_children_nodes)
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_predicate<I, T>(name: &str, args: I) -> Expr
    where
        I: IntoIterator<Item = T>,
        T: AsRef<Expr>,
    {
        let children_nodes = args
            .into_iter()
            .map(|child| {
                DAG_MANAGER
                    .get_or_create(child.as_ref())
                    .expect("Value is valid")
            })
            .collect::<Vec<_>>();
        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Predicate {
                    name: name.to_string(),
                },
                children_nodes,
            )
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_forall<A>(var: &str, expr: A) -> Expr
    where
        A: AsRef<Expr>,
    {
        let child_node = DAG_MANAGER
            .get_or_create(expr.as_ref())
            .expect("Value is valid");
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::ForAll(var.to_string()), vec![child_node])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_exists<A>(var: &str, expr: A) -> Expr
    where
        A: AsRef<Expr>,
    {
        let child_node = DAG_MANAGER
            .get_or_create(expr.as_ref())
            .expect("Value is valid");
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Exists(var.to_string()), vec![child_node])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_interval<A, B>(lower: A, upper: B, incl_lower: bool, incl_upper: bool) -> Expr
    where
        A: AsRef<Expr>,
        B: AsRef<Expr>,
    {
        let dag_lower = DAG_MANAGER
            .get_or_create(lower.as_ref())
            .expect("Value is valid");
        let dag_upper = DAG_MANAGER
            .get_or_create(upper.as_ref())
            .expect("Value is valid");
        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Interval(incl_lower, incl_upper),
                vec![dag_lower, dag_upper],
            )
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_sparse_polynomial(p: SparsePolynomial) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::SparsePolynomial(p), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    // --- Custom Constructors ---
    pub fn new_custom_zero() -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomZero, vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_custom_string(s: &str) -> Expr {
        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomString(s.to_string()), vec![])
            .expect("Value is valid");
        Expr::Dag(node)
    }

    unary_constructor!(new_custom_arc_one, CustomArcOne);
    binary_constructor!(new_custom_arc_two, CustomArcTwo);

    pub fn new_custom_arc_three<A, B, C>(a: A, b: B, c: C) -> Expr
    where
        A: AsRef<Expr>,
        B: AsRef<Expr>,
        C: AsRef<Expr>,
    {
        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");
        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");
        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");

        let children = vec![dag_a, dag_b, dag_c];

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomArcThree, children)
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_custom_arc_four<A, B, C, D>(a: A, b: B, c: C, d: D) -> Expr
    where
        A: AsRef<Expr>,
        B: AsRef<Expr>,
        C: AsRef<Expr>,
        D: AsRef<Expr>,
    {
        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");
        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");
        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");
        let dag_d = DAG_MANAGER
            .get_or_create(d.as_ref())
            .expect("Value is valid");

        let children = vec![dag_a, dag_b, dag_c, dag_d];

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomArcFour, children)
            .expect("Value is valid");
        Expr::Dag(node)
    }

    pub fn new_custom_arc_five<A, B, C, D, E>(a: A, b: B, c: C, d: D, e: E) -> Expr
    where
        A: AsRef<Expr>,
        B: AsRef<Expr>,
        C: AsRef<Expr>,
        D: AsRef<Expr>,
        E: AsRef<Expr>,
    {
        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");
        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");
        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");
        let dag_d = DAG_MANAGER
            .get_or_create(d.as_ref())
            .expect("Value is valid");
        let dag_e = DAG_MANAGER
            .get_or_create(e.as_ref())
            .expect("Value is valid");

        let children = vec![dag_a, dag_b, dag_c, dag_d, dag_e];

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomArcFive, children)
            .expect("Value is valid");
        Expr::Dag(node)
    }

    n_ary_constructor!(new_custom_vec_one, CustomVecOne);
    n_ary_constructor!(new_custom_vec_two, CustomVecTwo);
    n_ary_constructor!(new_custom_vec_three, CustomVecThree);
    n_ary_constructor!(new_custom_vec_four, CustomVecFour);
    n_ary_constructor!(new_custom_vec_five, CustomVecFive);
}
