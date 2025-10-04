//use std::collections::{BTreeMap, HashMap};
use crate::symbolic::unit_unification::UnitQuantity;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{self, Debug};
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex};

// --- Distribution Trait ---
// Moved here to break circular dependency
pub trait Distribution: Debug + Send + Sync {
    fn pdf(&self, x: &Expr) -> Expr;
    fn cdf(&self, x: &Expr) -> Expr;
    fn expectation(&self) -> Expr;
    fn variance(&self) -> Expr;
    fn mgf(&self, t: &Expr) -> Expr;
    fn clone_box(&self) -> Box<dyn Distribution>;
}

impl Clone for Box<dyn Distribution> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
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
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct SparsePolynomial {
    /// The terms of the polynomial, mapping each monomial to its coefficient.
    pub terms: BTreeMap<Monomial, Expr>,
}

/// The central enum representing a mathematical expression in the symbolic system.
///
/// `Expr` is an Abstract Syntax Tree (AST) that can represent a wide variety of
/// mathematical objects and operations. Manual implementations for `Debug`, `Clone`,
/// `PartialEq`, `Eq`, and `Hash` are provided to handle variants containing types
/// that do not derive these traits automatically (e.g., `f64`, `Box<dyn Distribution>`).
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
    Add(Box<Expr>, Box<Expr>),
    /// Subtraction of two expressions.
    Sub(Box<Expr>, Box<Expr>),
    /// Multiplication of two expressions.
    Mul(Box<Expr>, Box<Expr>),
    /// Division of two expressions.
    Div(Box<Expr>, Box<Expr>),
    /// Exponentiation (`base` raised to the power of `exponent`).
    Power(Box<Expr>, Box<Expr>),
    /// Negation of an expression.
    Neg(Box<Expr>),

    // --- Basic Mathematical Functions ---
    /// The sine function.
    Sin(Box<Expr>),
    /// The cosine function.
    Cos(Box<Expr>),
    /// The tangent function.
    Tan(Box<Expr>),
    /// The natural exponential function, `e^x`.
    Exp(Box<Expr>),
    /// The natural logarithm, `ln(x)`.
    Log(Box<Expr>),
    /// The absolute value function, `|x|`.
    Abs(Box<Expr>),
    /// The square root function.
    Sqrt(Box<Expr>),

    // --- Equations and Relations ---
    /// Represents an equation (`left = right`).
    Eq(Box<Expr>, Box<Expr>),
    /// Less than (`<`).
    Lt(Box<Expr>, Box<Expr>),
    /// Greater than (`>`).
    Gt(Box<Expr>, Box<Expr>),
    /// Less than or equal to (`<=`).
    Le(Box<Expr>, Box<Expr>),
    /// Greater than or equal to (`>=`).
    Ge(Box<Expr>, Box<Expr>),

    // --- Linear Algebra ---
    /// A matrix, represented as a vector of row vectors.
    Matrix(Vec<Vec<Expr>>),
    /// A vector (or column matrix).
    Vector(Vec<Expr>),
    /// A complex number with real and imaginary parts.
    Complex(Box<Expr>, Box<Expr>),
    /// Matrix transpose.
    Transpose(Box<Expr>),
    /// Matrix-matrix multiplication.
    MatrixMul(Box<Expr>, Box<Expr>),
    /// Matrix-vector multiplication.
    MatrixVecMul(Box<Expr>, Box<Expr>),
    /// Matrix inverse.
    Inverse(Box<Expr>),

    // --- Calculus ---
    /// The derivative of an expression with respect to a variable.
    Derivative(Box<Expr>, String),
    /// The N-th derivative of an expression.
    DerivativeN(Box<Expr>, String, Box<Expr>),
    /// A definite integral of `integrand` with respect to `var` from `lower_bound` to `upper_bound`.
    Integral {
        integrand: Box<Expr>,
        var: Box<Expr>,
        lower_bound: Box<Expr>,
        upper_bound: Box<Expr>,
    },
    /// A volume integral of a scalar field over a specified volume.
    VolumeIntegral {
        scalar_field: Box<Expr>,
        volume: Box<Expr>,
    },
    /// A surface integral of a vector field over a specified surface.
    SurfaceIntegral {
        vector_field: Box<Expr>,
        surface: Box<Expr>,
    },
    /// A limit of an expression as a variable approaches a point.
    Limit(Box<Expr>, String, Box<Expr>),

    // --- Series and Summations ---
    /// A summation of `body` with `var` from `from` to `to`.
    Sum {
        body: Box<Expr>,
        var: Box<Expr>,
        from: Box<Expr>,
        to: Box<Expr>,
    },
    /// A finite or infinite series expansion.
    Series(Box<Expr>, String, Box<Expr>, Box<Expr>),
    /// A summation over a range (similar to `Sum`).
    Summation(Box<Expr>, String, Box<Expr>, Box<Expr>),
    /// A product of terms over a range.
    Product(Box<Expr>, String, Box<Expr>, Box<Expr>),
    /// Represents a convergence analysis for a series.
    ConvergenceAnalysis(Box<Expr>, String),
    /// An asymptotic expansion of a function.
    AsymptoticExpansion(Box<Expr>, String, Box<Expr>, Box<Expr>),

    // --- Trigonometric & Hyperbolic Functions (Extended) ---
    /// Secant function.
    Sec(Box<Expr>),
    /// Cosecant function.
    Csc(Box<Expr>),
    /// Cotangent function.
    Cot(Box<Expr>),
    /// Arcsine (inverse sine).
    ArcSin(Box<Expr>),
    /// Arccosine (inverse cosine).
    ArcCos(Box<Expr>),
    /// Arctangent (inverse tangent).
    ArcTan(Box<Expr>),
    /// Arcsecant (inverse secant).
    ArcSec(Box<Expr>),
    /// Arccosecant (inverse cosecant).
    ArcCsc(Box<Expr>),
    /// Arccotangent (inverse cotangent).
    ArcCot(Box<Expr>),
    /// Hyperbolic sine.
    Sinh(Box<Expr>),
    /// Hyperbolic cosine.
    Cosh(Box<Expr>),
    /// Hyperbolic tangent.
    Tanh(Box<Expr>),
    /// Hyperbolic secant.
    Sech(Box<Expr>),
    /// Hyperbolic cosecant.
    Csch(Box<Expr>),
    /// Hyperbolic cotangent.
    Coth(Box<Expr>),
    /// Inverse hyperbolic sine.
    ArcSinh(Box<Expr>),
    /// Inverse hyperbolic cosine.
    ArcCosh(Box<Expr>),
    /// Inverse hyperbolic tangent.
    ArcTanh(Box<Expr>),
    /// Inverse hyperbolic secant.
    ArcSech(Box<Expr>),
    /// Inverse hyperbolic cosecant.
    ArcCsch(Box<Expr>),
    /// Inverse hyperbolic cotangent.
    ArcCoth(Box<Expr>),
    /// Logarithm with a specified base.
    LogBase(Box<Expr>, Box<Expr>),
    /// Two-argument arctangent.
    Atan2(Box<Expr>, Box<Expr>),

    // --- Combinatorics ---
    /// Binomial coefficient, "n choose k".
    Binomial(Box<Expr>, Box<Expr>),
    /// Factorial, `n!`.
    Factorial(Box<Expr>),
    /// Permutations, `P(n, k)`.
    Permutation(Box<Expr>, Box<Expr>),
    /// Combinations, `C(n, k)`.
    Combination(Box<Expr>, Box<Expr>),
    /// Falling factorial.
    FallingFactorial(Box<Expr>, Box<Expr>),
    /// Rising factorial.
    RisingFactorial(Box<Expr>, Box<Expr>),

    // --- Geometry & Vector Calculus ---
    /// A path for path integrals (e.g., line, circle).
    Path(PathType, Box<Expr>, Box<Expr>),
    /// Represents the boundary of a domain.
    Boundary(Box<Expr>),
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
    Gamma(Box<Expr>),
    /// The Beta function.
    Beta(Box<Expr>, Box<Expr>),
    /// The error function.
    Erf(Box<Expr>),
    /// The complementary error function.
    Erfc(Box<Expr>),
    /// The imaginary error function.
    Erfi(Box<Expr>),
    /// The Riemann Zeta function.
    Zeta(Box<Expr>),
    /// Bessel function of the first kind.
    BesselJ(Box<Expr>, Box<Expr>),
    /// Bessel function of the second kind.
    BesselY(Box<Expr>, Box<Expr>),
    /// Legendre polynomial.
    LegendreP(Box<Expr>, Box<Expr>),
    /// Laguerre polynomial.
    LaguerreL(Box<Expr>, Box<Expr>),
    /// Hermite polynomial.
    HermiteH(Box<Expr>, Box<Expr>),
    /// The digamma function (psi function).
    Digamma(Box<Expr>),
    /// The Kronecker delta function.
    KroneckerDelta(Box<Expr>, Box<Expr>),

    // --- Logic & Sets ---
    /// Logical AND of a vector of expressions.
    And(Vec<Expr>),
    /// Logical OR of a vector of expressions.
    Or(Vec<Expr>),
    /// Logical NOT.
    Not(Box<Expr>),
    /// Logical XOR (exclusive OR).
    Xor(Box<Expr>, Box<Expr>),
    /// Logical implication (`A => B`).
    Implies(Box<Expr>, Box<Expr>),
    /// Logical equivalence (`A <=> B`).
    Equivalent(Box<Expr>, Box<Expr>),
    /// A predicate with a name and arguments.
    Predicate { name: String, args: Vec<Expr> },
    /// Universal quantifier ("for all").
    ForAll(String, Box<Expr>),
    /// Existential quantifier ("there exists").
    Exists(String, Box<Expr>),
    /// A union of sets or intervals.
    Union(Vec<Expr>),
    /// An interval with a lower and upper bound, and flags for inclusion.
    Interval(Box<Expr>, Box<Expr>, bool, bool),

    // --- Polynomials & Number Theory ---
    /// A dense polynomial represented by its coefficients.
    Polynomial(Vec<Expr>),
    /// A sparse polynomial.
    SparsePolynomial(SparsePolynomial),
    /// The floor function.
    Floor(Box<Expr>),
    /// A predicate to check if a number is prime.
    IsPrime(Box<Expr>),
    /// Greatest Common Divisor (GCD).
    Gcd(Box<Expr>, Box<Expr>),
    /// Modulo operation.
    Mod(Box<Expr>, Box<Expr>),

    // --- Solving & Substitution ---
    /// Represents the action of solving an equation for a variable.
    Solve(Box<Expr>, String),
    /// Represents the substitution of a variable in an expression with another expression.
    Substitute(Box<Expr>, String, Box<Expr>),
    /// Represents a system of equations to be solved.
    System(Vec<Expr>),
    /// Represents the set of solutions to an equation or system.
    Solutions(Vec<Expr>),
    /// A parametric solution, e.g., for a system of ODEs.
    ParametricSolution { x: Box<Expr>, y: Box<Expr> },
    /// Represents the `i`-th root of a polynomial.
    RootOf { poly: Box<Expr>, index: u32 },
    /// Represents infinite solutions.
    InfiniteSolutions,
    /// Represents that no solution exists.
    NoSolution,

    // --- Differential Equations ---
    /// An ordinary differential equation (ODE).
    Ode {
        equation: Box<Expr>,
        func: String,
        var: String,
    },
    /// A partial differential equation (PDE).
    Pde {
        equation: Box<Expr>,
        func: String,
        vars: Vec<String>,
    },
    /// The general solution to a differential equation.
    GeneralSolution(Box<Expr>),
    /// A particular solution to a differential equation.
    ParticularSolution(Box<Expr>),

    // --- Integral Equations ---
    /// A Fredholm integral equation.
    Fredholm(Box<Expr>, Box<Expr>, Box<Expr>, Box<Expr>),
    /// A Volterra integral equation.
    Volterra(Box<Expr>, Box<Expr>, Box<Expr>, Box<Expr>),

    // --- Miscellaneous ---
    /// Application of a function to an argument.
    Apply(Box<Expr>, Box<Expr>),
    /// A tuple of expressions.
    Tuple(Vec<Expr>),
    /// A node in a Directed Acyclic Graph (DAG) for expression sharing.
    #[serde(skip_serializing, skip_deserializing)]
    Dag(Arc<DagNode>),
    /// A probability distribution.
    #[serde(skip_serializing, skip_deserializing)]
    Distribution(Box<dyn Distribution>),
    /// Maximum of two expressions.
    Max(Box<Expr>, Box<Expr>),
    /// A unified quantity with its value and unit string.
    Quantity(Box<UnitQuantity>),
    /// A temporary representation of a value with a unit string, before unification.
    QuantityWithValue(Box<Expr>, String),
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
                    s.push_str(&format!("{}*x^{}", coeff, i));
                }
                write!(f, "{}", s)
            }
            Expr::SparsePolynomial(p) => {
                let mut s = String::new();
                for (monomial, coeff) in &p.terms {
                    if !s.is_empty() {
                        s.push_str(" + ");
                    }
                    s.push_str(&format!("({})", coeff));
                    for (var, exp) in &monomial.0 {
                        s.push_str(&format!("*{}^{}", var, exp));
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
            _ => write!(f, "Unimplemented Display for Expr variant"),
        }
    }
}

impl Expr {
    pub fn re(&self) -> Self {
        if let Expr::Complex(re, _) = self {
            *re.clone()
        } else {
            self.clone()
        }
    }

    pub fn im(&self) -> Self {
        if let Expr::Complex(_, im) = self {
            *im.clone()
        } else {
            Expr::Constant(0.0)
        }
    }

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
    Constant(OrderedFloat<f64>),
    BigInt(BigInt),
    Rational(BigRational),
    Boolean(bool),
    Variable(String),
    Pattern(String),
    Add,
    Sub,
    Mul,
    Div,
    Power,
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
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

impl DagNode {
    pub fn new(op: DagOp, children: Vec<Arc<DagNode>>) -> Arc<Self> {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        op.hash(&mut hasher);
        children.hash(&mut hasher);
        let hash = hasher.finish();
        Arc::new(DagNode { op, children, hash })
    }
}

pub struct DagManager {
    nodes: Mutex<HashMap<u64, Arc<DagNode>>>,
}

impl Default for DagManager {
    fn default() -> Self {
        Self::new()
    }
}

impl DagManager {
    pub fn new() -> Self {
        DagManager {
            nodes: Mutex::new(HashMap::new()),
        }
    }

    pub fn get_or_create(&self, op: DagOp, children: Vec<Arc<DagNode>>) -> Arc<DagNode> {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        op.hash(&mut hasher);
        children.hash(&mut hasher);
        let hash = hasher.finish();

        let mut nodes = self.nodes.lock().unwrap();
        if let Some(node) = nodes.get(&hash) {
            return node.clone();
        }

        let new_node = Arc::new(DagNode { op, children, hash });
        nodes.insert(hash, new_node.clone());
        new_node
    }
}

impl PartialEq for Expr {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Expr::Constant(a), Expr::Constant(b)) => a == b,
            (Expr::BigInt(a), Expr::BigInt(b)) => a == b,
            (Expr::Rational(a), Expr::Rational(b)) => a == b,
            (Expr::Boolean(a), Expr::Boolean(b)) => a == b,
            (Expr::Variable(a), Expr::Variable(b)) => a == b,
            (Expr::Pattern(a), Expr::Pattern(b)) => a == b,
            (Expr::Pi, Expr::Pi) | (Expr::E, Expr::E) => true,
            (Expr::Add(a1, b1), Expr::Add(a2, b2)) => {
                (a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2)
            }
            (Expr::Sub(a1, b1), Expr::Sub(a2, b2)) => a1 == a2 && b1 == b2,
            (Expr::Mul(a1, b1), Expr::Mul(a2, b2)) => {
                (a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2)
            }
            (Expr::Div(a1, b1), Expr::Div(a2, b2)) => a1 == a2 && b1 == b2,
            (Expr::Power(a1, b1), Expr::Power(a2, b2)) => a1 == a2 && b1 == b2,
            (Expr::Sin(a1), Expr::Sin(a2)) => a1 == a2,
            (Expr::Cos(a1), Expr::Cos(a2)) => a1 == a2,
            (Expr::Tan(a1), Expr::Tan(a2)) => a1 == a2,
            (Expr::Exp(a1), Expr::Exp(a2)) => a1 == a2,
            (Expr::Log(a1), Expr::Log(a2)) => a1 == a2,
            (Expr::Dag(a), Expr::Dag(b)) => Arc::ptr_eq(a, b),
            (Expr::And(v1), Expr::And(v2)) => {
                v1.iter().all(|item| v2.contains(item)) && v2.iter().all(|item| v1.contains(item))
            }
            (Expr::Or(v1), Expr::Or(v2)) => {
                v1.iter().all(|item| v2.contains(item)) && v2.iter().all(|item| v1.contains(item))
            }
            (Expr::Not(a1), Expr::Not(a2)) => a1 == a2,
            (Expr::Equivalent(a1, b1), Expr::Equivalent(a2, b2)) => a1 == a2 && b1 == b2,
            (Expr::Predicate { name: n1, args: a1 }, Expr::Predicate { name: n2, args: a2 }) => {
                n1 == n2 && a1 == a2
            }
            (Expr::ForAll(s1, e1), Expr::ForAll(s2, e2)) => s1 == s2 && e1 == e2,
            (Expr::Exists(s1, e1), Expr::Exists(s2, e2)) => s1 == s2 && e1 == e2,
            (Expr::Distribution(_), Expr::Distribution(_)) => false, // Cannot compare distributions
            (Expr::System(v1), Expr::System(v2)) => v1 == v2,
            (Expr::Solutions(v1), Expr::Solutions(v2)) => v1 == v2,
            (
                Expr::ParametricSolution { x: x1, y: y1 },
                Expr::ParametricSolution { x: x2, y: y2 },
            ) => x1 == x2 && y1 == y2,
            (
                Expr::RootOf {
                    poly: p1,
                    index: i1,
                },
                Expr::RootOf {
                    poly: p2,
                    index: i2,
                },
            ) => p1 == p2 && i1 == i2,
            (
                Expr::Sum {
                    body: b1,
                    var: v1,
                    from: f1,
                    to: t1,
                },
                Expr::Sum {
                    body: b2,
                    var: v2,
                    from: f2,
                    to: t2,
                },
            ) => b1 == b2 && v1 == v2 && f1 == f2 && t1 == t2,
            (Expr::Erfc(a1), Expr::Erfc(a2)) => a1 == a2,
            (Expr::Erfi(a1), Expr::Erfi(a2)) => a1 == a2,
            (Expr::Zeta(a1), Expr::Zeta(a2)) => a1 == a2,
            (Expr::SparsePolynomial(p1), Expr::SparsePolynomial(p2)) => p1 == p2,
            (Expr::Quantity(q1), Expr::Quantity(q2)) => q1 == q2,
            (Expr::QuantityWithValue(v1, u1), Expr::QuantityWithValue(v2, u2)) => {
                v1 == v2 && u1 == u2
            }
            _ => false,
        }
    }
}

impl Eq for Expr {}

impl Hash for Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match self {
            Expr::Constant(f) => OrderedFloat(*f).hash(state),
            Expr::BigInt(i) => i.hash(state),
            Expr::Rational(r) => r.hash(state),
            Expr::Boolean(b) => b.hash(state),
            Expr::Variable(s) | Expr::Pattern(s) => s.hash(state),
            Expr::Pi => "pi".hash(state),
            Expr::E => "e".hash(state),
            Expr::Add(a, b)
            | Expr::Mul(a, b)
            | Expr::Xor(a, b)
            | Expr::Implies(a, b)
            | Expr::Equivalent(a, b)
            | Expr::Gcd(a, b) => {
                // Hash in a commutative way
                let mut h1 = std::collections::hash_map::DefaultHasher::new();
                let mut h2 = std::collections::hash_map::DefaultHasher::new();
                a.hash(&mut h1);
                b.hash(&mut h2);
                if h1.finish() < h2.finish() {
                    a.hash(state);
                    b.hash(state);
                } else {
                    b.hash(state);
                    a.hash(state);
                }
            }
            Expr::Predicate { name, args } => {
                name.hash(state);
                args.hash(state);
            }
            Expr::Sub(a, b) | Expr::Div(a, b) | Expr::Power(a, b) | Expr::Eq(a, b) => {
                a.hash(state);
                b.hash(state);
            }
            Expr::ForAll(s, e) | Expr::Exists(s, e) => {
                s.hash(state);
                e.hash(state);
            }
            Expr::Sin(a)
            | Expr::Cos(a)
            | Expr::Tan(a)
            | Expr::Exp(a)
            | Expr::Log(a)
            | Expr::Neg(a)
            | Expr::Not(a)
            | Expr::Factorial(a)
            | Expr::Floor(a)
            | Expr::IsPrime(a) => a.hash(state),
            Expr::And(v) | Expr::Or(v) => {
                let mut hashes: Vec<_> = v
                    .iter()
                    .map(|e| {
                        let mut h = std::collections::hash_map::DefaultHasher::new();
                        e.hash(&mut h);
                        h.finish()
                    })
                    .collect();
                hashes.sort();
                hashes.hash(state);
            }
            Expr::Dag(node) => node.hash.hash(state),
            Expr::System(v) => v.hash(state),
            Expr::Solutions(v) => v.hash(state),
            Expr::ParametricSolution { x, y } => {
                x.hash(state);
                y.hash(state);
            }
            Expr::RootOf { poly, index } => {
                poly.hash(state);
                index.hash(state);
            }
            Expr::Sum {
                body,
                var,
                from,
                to,
            } => {
                body.hash(state);
                var.hash(state);
                from.hash(state);
                to.hash(state);
            }
            /*
            Expr::Sum { body, var, from, to } => {
                body.hash(state);
                var.hash(state);
                from.hash(state);
                to.hash(state);
            }
            */
            Expr::Erfc(a) => a.hash(state),
            Expr::Erfi(a) => a.hash(state),
            Expr::Zeta(a) => a.hash(state),
            Expr::SparsePolynomial(p) => {
                for (monomial, coeff) in &p.terms {
                    monomial.hash(state);
                    coeff.hash(state);
                }
            }
            Expr::Quantity(q) => q.hash(state),
            Expr::QuantityWithValue(v, u) => {
                v.hash(state);
                u.hash(state);
            }
            // Note: Hashing for many variants is omitted for brevity, but should be implemented
            _ => {}
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
        let self_order = self.variant_order();
        let other_order = other.variant_order();

        if self_order != other_order {
            return self_order.cmp(&other_order);
        }

        match (self, other) {
            (Expr::Constant(a), Expr::Constant(b)) => OrderedFloat(*a).cmp(&OrderedFloat(*b)),
            (Expr::BigInt(a), Expr::BigInt(b)) => a.cmp(b),
            (Expr::Rational(a), Expr::Rational(b)) => a.cmp(b),
            (Expr::Boolean(a), Expr::Boolean(b)) => a.cmp(b),
            (Expr::Variable(a), Expr::Variable(b)) => a.cmp(b),
            (Expr::Pattern(a), Expr::Pattern(b)) => a.cmp(b),
            (Expr::Pi, Expr::Pi) | (Expr::E, Expr::E) => Ordering::Equal,
            (Expr::Add(a1, b1), Expr::Add(a2, b2)) | (Expr::Mul(a1, b1), Expr::Mul(a2, b2)) => {
                let (s1, s2) = if a1.cmp(b1) < Ordering::Equal {
                    (a1, b1)
                } else {
                    (b1, a1)
                };
                let (o1, o2) = if a2.cmp(b2) < Ordering::Equal {
                    (a2, b2)
                } else {
                    (b2, a2)
                };
                s1.cmp(o1).then_with(|| s2.cmp(o2))
            }
            (Expr::Sub(a1, b1), Expr::Sub(a2, b2))
            | (Expr::Div(a1, b1), Expr::Div(a2, b2))
            | (Expr::Power(a1, b1), Expr::Power(a2, b2))
            | (Expr::Eq(a1, b1), Expr::Eq(a2, b2))
            | (Expr::Complex(a1, b1), Expr::Complex(a2, b2))
            | (Expr::LogBase(a1, b1), Expr::LogBase(a2, b2))
            | (Expr::Atan2(a1, b1), Expr::Atan2(a2, b2))
            | (Expr::Binomial(a1, b1), Expr::Binomial(a2, b2))
            | (Expr::Beta(a1, b1), Expr::Beta(a2, b2))
            | (Expr::BesselJ(a1, b1), Expr::BesselJ(a2, b2))
            | (Expr::BesselY(a1, b1), Expr::BesselY(a2, b2))
            | (Expr::LegendreP(a1, b1), Expr::LegendreP(a2, b2))
            | (Expr::LaguerreL(a1, b1), Expr::LaguerreL(a2, b2))
            | (Expr::HermiteH(a1, b1), Expr::HermiteH(a2, b2))
            | (Expr::KroneckerDelta(a1, b1), Expr::KroneckerDelta(a2, b2))
            | (Expr::Lt(a1, b1), Expr::Lt(a2, b2))
            | (Expr::Gt(a1, b1), Expr::Gt(a2, b2))
            | (Expr::Le(a1, b1), Expr::Le(a2, b2))
            | (Expr::Ge(a1, b1), Expr::Ge(a2, b2))
            | (Expr::Permutation(a1, b1), Expr::Permutation(a2, b2))
            | (Expr::Combination(a1, b1), Expr::Combination(a2, b2))
            | (Expr::FallingFactorial(a1, b1), Expr::FallingFactorial(a2, b2))
            | (Expr::RisingFactorial(a1, b1), Expr::RisingFactorial(a2, b2))
            | (Expr::Xor(a1, b1), Expr::Xor(a2, b2))
            | (Expr::Implies(a1, b1), Expr::Implies(a2, b2))
            | (Expr::Equivalent(a1, b1), Expr::Equivalent(a2, b2))
            | (Expr::Gcd(a1, b1), Expr::Gcd(a2, b2))
            | (Expr::Mod(a1, b1), Expr::Mod(a2, b2))
            | (Expr::Max(a1, b1), Expr::Max(a2, b2))
            | (Expr::MatrixMul(a1, b1), Expr::MatrixMul(a2, b2))
            | (Expr::MatrixVecMul(a1, b1), Expr::MatrixVecMul(a2, b2)) => {
                a1.cmp(a2).then_with(|| b1.cmp(b2))
            }

            (Expr::Predicate { name: n1, args: a1 }, Expr::Predicate { name: n2, args: a2 }) => {
                n1.cmp(n2).then_with(|| a1.cmp(a2))
            }

            (Expr::ForAll(s1, e1), Expr::ForAll(s2, e2)) => s1.cmp(s2).then_with(|| e1.cmp(e2)),
            (Expr::Exists(s1, e1), Expr::Exists(s2, e2)) => s1.cmp(s2).then_with(|| e1.cmp(e2)),

            (Expr::System(v1), Expr::System(v2)) => v1.cmp(v2),
            (Expr::Solutions(v1), Expr::Solutions(v2)) => v1.cmp(v2),
            (
                Expr::ParametricSolution { x: x1, y: y1 },
                Expr::ParametricSolution { x: x2, y: y2 },
            ) => x1.cmp(x2).then_with(|| y1.cmp(y2)),

            (
                Expr::RootOf {
                    poly: p1,
                    index: i1,
                },
                Expr::RootOf {
                    poly: p2,
                    index: i2,
                },
            ) => p1.cmp(p2).then_with(|| i1.cmp(i2)),
            (
                Expr::Sum {
                    body: b1,
                    var: v1,
                    from: f1,
                    to: t1,
                },
                Expr::Sum {
                    body: b2,
                    var: v2,
                    from: f2,
                    to: t2,
                },
            ) => b1
                .cmp(b2)
                .then_with(|| v1.cmp(v2))
                .then_with(|| f1.cmp(f2))
                .then_with(|| t1.cmp(t2)),

            (Expr::Erfc(a1), Expr::Erfc(a2)) => a1.cmp(a2),
            (Expr::Erfi(a1), Expr::Erfi(a2)) => a1.cmp(a2),
            (Expr::Zeta(a1), Expr::Zeta(a2)) => a1.cmp(a2),

            (Expr::Sin(a1), Expr::Sin(a2))
            | (Expr::Cos(a1), Expr::Cos(a2))
            | (Expr::Tan(a1), Expr::Tan(a2))
            | (Expr::Exp(a1), Expr::Exp(a2))
            | (Expr::Log(a1), Expr::Log(a2))
            | (Expr::Neg(a1), Expr::Neg(a2))
            | (Expr::Abs(a1), Expr::Abs(a2))
            | (Expr::Sqrt(a1), Expr::Sqrt(a2))
            | (Expr::Sec(a1), Expr::Sec(a2))
            | (Expr::Csc(a1), Expr::Csc(a2))
            | (Expr::Cot(a1), Expr::Cot(a2))
            | (Expr::ArcSin(a1), Expr::ArcSin(a2))
            | (Expr::ArcCos(a1), Expr::ArcCos(a2))
            | (Expr::ArcTan(a1), Expr::ArcTan(a2))
            | (Expr::ArcSec(a1), Expr::ArcSec(a2))
            | (Expr::ArcCsc(a1), Expr::ArcCsc(a2))
            | (Expr::ArcCot(a1), Expr::ArcCot(a2))
            | (Expr::Sinh(a1), Expr::Sinh(a2))
            | (Expr::Cosh(a1), Expr::Cosh(a2))
            | (Expr::Tanh(a1), Expr::Tanh(a2))
            | (Expr::Sech(a1), Expr::Sech(a2))
            | (Expr::Csch(a1), Expr::Csch(a2))
            | (Expr::Coth(a1), Expr::Coth(a2))
            | (Expr::ArcSinh(a1), Expr::ArcSinh(a2))
            | (Expr::ArcCosh(a1), Expr::ArcCosh(a2))
            | (Expr::ArcTanh(a1), Expr::ArcTanh(a2))
            | (Expr::ArcSech(a1), Expr::ArcSech(a2))
            | (Expr::ArcCsch(a1), Expr::ArcCsch(a2))
            | (Expr::ArcCoth(a1), Expr::ArcCoth(a2))
            | (Expr::Boundary(a1), Expr::Boundary(a2))
            | (Expr::Gamma(a1), Expr::Gamma(a2))
            | (Expr::Erf(a1), Expr::Erf(a2))
            | (Expr::Digamma(a1), Expr::Digamma(a2))
            | (Expr::Not(a1), Expr::Not(a2))
            | (Expr::Floor(a1), Expr::Floor(a2))
            | (Expr::IsPrime(a1), Expr::IsPrime(a2))
            | (Expr::Factorial(a1), Expr::Factorial(a2))
            | (Expr::Transpose(a1), Expr::Transpose(a2))
            | (Expr::Inverse(a1), Expr::Inverse(a2))
            | (Expr::GeneralSolution(a1), Expr::GeneralSolution(a2))
            | (Expr::ParticularSolution(a1), Expr::ParticularSolution(a2)) => a1.cmp(a2),
            //(Expr::Solve(a1, _), Expr::Solve(a2, _)) => a1.cmp(a2),
            //(Expr::ConvergenceAnalysis(a1, _), Expr::ConvergenceAnalysis(a2, _)) => a1.cmp(a2),
            (Expr::Matrix(m1), Expr::Matrix(m2)) => m1.cmp(m2),
            (Expr::Vector(v1), Expr::Vector(v2)) => v1.cmp(v2),
            (Expr::Tuple(t1), Expr::Tuple(t2)) => t1.cmp(t2),
            (Expr::Polynomial(p1), Expr::Polynomial(p2)) => p1.cmp(p2),
            (Expr::SparsePolynomial(p1), Expr::SparsePolynomial(p2)) => {
                p1.terms.iter().cmp(p2.terms.iter())
            }
            (Expr::And(v1), Expr::And(v2))
            | (Expr::Or(v1), Expr::Or(v2))
            | (Expr::Union(v1), Expr::Union(v2)) => {
                let mut sorted1 = v1.clone();
                let mut sorted2 = v2.clone();
                sorted1.sort();
                sorted2.sort();
                sorted1.cmp(&sorted2)
            }
            (Expr::Derivative(e1, s1), Expr::Derivative(e2, s2)) => {
                e1.cmp(e2).then_with(|| s1.cmp(s2))
            }
            (Expr::Quantity(_), Expr::Quantity(_)) => std::cmp::Ordering::Equal, // Cannot be ordered
            (Expr::QuantityWithValue(v1, u1), Expr::QuantityWithValue(v2, u2)) => {
                v1.cmp(v2).then_with(|| u1.cmp(u2))
            }
            //(Expr::Derivative(a1, _), Expr::Derivative(a2, _)) => a1.cmp(a2),
            (
                Expr::Integral {
                    integrand: i1,
                    lower_bound: l1,
                    upper_bound: u1,
                    ..
                },
                Expr::Integral {
                    integrand: i2,
                    lower_bound: l2,
                    upper_bound: u2,
                    ..
                },
            ) => i1.cmp(i2).then_with(|| l1.cmp(l2)).then_with(|| u1.cmp(u2)),
            (Expr::Path(pt1, p1_1, p1_2), Expr::Path(pt2, p2_1, p2_2)) => pt1
                .cmp(pt2)
                .then_with(|| p1_1.cmp(p2_1))
                .then_with(|| p1_2.cmp(p2_2)),
            (Expr::Domain(s1), Expr::Domain(s2)) => s1.cmp(s2),
            (
                Expr::VolumeIntegral {
                    scalar_field: s1,
                    volume: v1,
                },
                Expr::VolumeIntegral {
                    scalar_field: s2,
                    volume: v2,
                },
            ) => s1.cmp(s2).then_with(|| v1.cmp(v2)),
            (
                Expr::SurfaceIntegral {
                    vector_field: v1,
                    surface: s1,
                },
                Expr::SurfaceIntegral {
                    vector_field: v2,
                    surface: s2,
                },
            ) => v1.cmp(v2).then_with(|| s1.cmp(s2)),
            (Expr::Infinity, Expr::Infinity)
            | (Expr::NegativeInfinity, Expr::NegativeInfinity)
            | (Expr::InfiniteSolutions, Expr::InfiniteSolutions)
            | (Expr::NoSolution, Expr::NoSolution) => Ordering::Equal,
            (Expr::Apply(a1, b1), Expr::Apply(a2, b2)) => a1.cmp(a2).then_with(|| b1.cmp(b2)),
            (Expr::DerivativeN(e1, s1, n1), Expr::DerivativeN(e2, s2, n2)) => {
                e1.cmp(e2).then_with(|| s1.cmp(s2)).then_with(|| n1.cmp(n2))
            }
            (Expr::Series(a1, b1, c1, d1), Expr::Series(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::Summation(a1, b1, c1, d1), Expr::Summation(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::Product(a1, b1, c1, d1), Expr::Product(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (
                Expr::AsymptoticExpansion(a1, b1, c1, d1),
                Expr::AsymptoticExpansion(a2, b2, c2, d2),
            ) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::Fredholm(a1, b1, c1, d1), Expr::Fredholm(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::Volterra(a1, b1, c1, d1), Expr::Volterra(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::ConvergenceAnalysis(e1, s1), Expr::ConvergenceAnalysis(e2, s2)) => {
                e1.cmp(e2).then_with(|| s1.cmp(s2))
            }
            (Expr::Interval(a1, b1, c1, d1), Expr::Interval(a2, b2, c2, d2)) => a1
                .cmp(a2)
                .then_with(|| b1.cmp(b2))
                .then_with(|| c1.cmp(c2))
                .then_with(|| d1.cmp(d2)),
            (Expr::Solve(e1, s1), Expr::Solve(e2, s2)) => e1.cmp(e2).then_with(|| s1.cmp(s2)),
            (Expr::Substitute(a1, b1, c1), Expr::Substitute(a2, b2, c2)) => {
                a1.cmp(a2).then_with(|| b1.cmp(b2)).then_with(|| c1.cmp(c2))
            }
            (Expr::Limit(a1, b1, c1), Expr::Limit(a2, b2, c2)) => {
                a1.cmp(a2).then_with(|| b1.cmp(b2)).then_with(|| c1.cmp(c2))
            }
            (Expr::Dag(n1), Expr::Dag(n2)) => n1.hash.cmp(&n2.hash),
            (
                Expr::Ode {
                    equation: e1,
                    func: f1,
                    var: v1,
                },
                Expr::Ode {
                    equation: e2,
                    func: f2,
                    var: v2,
                },
            ) => e1.cmp(e2).then_with(|| f1.cmp(f2)).then_with(|| v1.cmp(v2)),
            (
                Expr::Pde {
                    equation: e1,
                    func: f1,
                    vars: v1,
                },
                Expr::Pde {
                    equation: e2,
                    func: f2,
                    vars: v2,
                },
            ) => e1.cmp(e2).then_with(|| f1.cmp(f2)).then_with(|| v1.cmp(v2)),
            (Expr::Distribution(_), Expr::Distribution(_)) => Ordering::Equal, // Cannot be ordered
            _ => Ordering::Equal,
        }
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
                v.iter().for_each(|e| e.in_order_walk(f));
            }
            Expr::Predicate { args, .. } => {
                f(self);
                args.iter().for_each(|e| e.in_order_walk(f));
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
            /*
            Expr::Sum { body, var, from, to } => {
                f(self);
                body.in_order_walk(f);
                var.in_order_walk(f);
                from.in_order_walk(f);
                to.in_order_walk(f);
            }
            */
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
            | Expr::Distribution(_) => {}
        }
        f(self); // Visit parent
    }
}
