//! # Symbolic Expression Core Module
//!
//! This module provides the foundational data structures and operations for symbolic
//! mathematics in the RSSN library. It implements a hybrid AST/DAG (Abstract Syntax Tree /
//! Directed Acyclic Graph) representation system for mathematical expressions.
//!
//! ## Architecture Overview
//!
//! The core of this module is the [`Expr`] enum, which represents symbolic mathematical
//! expressions. The system is designed around a dual representation strategy:
//!
//! 1. **AST Representation**: Traditional tree-based expression structure (legacy)
//! 2. **DAG Representation**: Graph-based structure with shared subexpressions (modern)
//!
//! The DAG representation is managed by the [`DagManager`], which provides:
//! - Automatic deduplication of identical subexpressions
//! - Hash-based node lookup for O(1) retrieval
//! - Canonical normalization of expressions
//! - Memory-efficient storage through structural sharing
//!
//! ## Key Components
//!
//! ### Expression Types ([`Expr`])
//!
//! The [`Expr`] enum supports a comprehensive set of mathematical operations:
//!
//! - **Atomic Values**: Constants, variables, patterns, special constants (π, e, ∞)
//! - **Arithmetic**: Addition, subtraction, multiplication, division, power, negation
//! - **Trigonometric**: sin, cos, tan, sec, csc, cot and their inverses
//! - **Hyperbolic**: sinh, cosh, tanh, sech, csch, coth and their inverses
//! - **Special Functions**: Gamma, Beta, Bessel, Legendre, Laguerre, Hermite, etc.
//! - **Calculus**: Derivatives, integrals, limits, series, summations
//! - **Linear Algebra**: Matrices, vectors, transpose, inverse, matrix multiplication
//! - **Logic**: Boolean operations, predicates, quantifiers
//! - **Advanced**: ODEs, PDEs, distributions, complex analysis
//!
//! ### N-ary Operations
//!
//! The module includes efficient n-ary operation variants:
//!
//! - [`Expr::AddList`]: Sum of multiple terms in a single operation
//! - [`Expr::MulList`]: Product of multiple factors in a single operation
//!
//! These variants improve performance by reducing tree depth and enabling better
//! optimization opportunities during simplification.
//!
//! ### Dynamic Operations
//!
//! The module supports runtime-extensible operations through:
//!
//! - [`Expr::UnaryList`]: Custom unary operations
//! - [`Expr::BinaryList`]: Custom binary operations  
//! - [`Expr::NaryList`]: Custom n-ary operations
//!
//! These are registered in the [`DYNAMIC_OP_REGISTRY`] with properties like
//! associativity and commutativity, enabling plugin systems and domain-specific
//! extensions without modifying the core enum.
//!
//! ### DAG Management
//!
//! The [`DagManager`] provides centralized management of DAG nodes:
//!
//! ```rust
//! use rssn::symbolic::core::{Expr, DAG_MANAGER};
//!
//! // Create expressions using smart constructors
//! let x = Expr::new_variable("x");
//! let two = Expr::new_constant(2.0);
//! let expr = Expr::new_add(x, two);
//!
//! // The DAG_MANAGER automatically deduplicates identical subexpressions
//! ```
//!
//! ### Smart Constructors
//!
//! All operations have corresponding smart constructors (e.g., `new_add`, `new_mul`)
//! that automatically:
//! - Convert to DAG representation
//! - Normalize the expression
//! - Deduplicate subexpressions
//! - Apply basic simplifications
//!
//! ## AST to DAG Migration
//!
//! The module is undergoing a gradual migration from AST to DAG representation:
//!
//! - **Legacy AST forms** (e.g., `Expr::Add(Arc<Expr>, Arc<Expr>)`) remain for compatibility
//! - **Modern DAG forms** use `Expr::Dag(Arc<DagNode>)` wrapper
//! - **Smart constructors** automatically create DAG forms
//! - **Conversion utilities** (`to_dag()`, `to_ast()`) enable interoperability
//!
//! This hybrid approach ensures backward compatibility while enabling new optimizations.
//!
//! ## Expression Traversal
//!
//! The module provides multiple traversal methods:
//!
//! - [`Expr::pre_order_walk`]: Visit parent before children
//! - [`Expr::post_order_walk`]: Visit children before parent
//! - [`Expr::in_order_walk`]: Visit left child, parent, then right child
//!
//! ## Normalization and Canonicalization
//!
//! The [`Expr::normalize`] method provides canonical forms:
//! - Sorts commutative operation children
//! - Flattens nested associative operations
//! - Applies consistent ordering for hashing
//!
//! ## Examples
//!
//! ### Basic Expression Creation
//!
//! ```rust
//! use rssn::symbolic::core::Expr;
//!
//! // Using smart constructors (recommended)
//! let x = Expr::new_variable("x");
//! let y = Expr::new_variable("y");
//! let sum = Expr::new_add(x.clone(), y.clone());
//! let product = Expr::new_mul(x, y);
//! ```
//!
//! ### N-ary Operations
//!
//! ```rust
//! use rssn::symbolic::core::Expr;
//!
//! // Efficient multi-term addition
//! let sum = Expr::AddList(vec![
//!     Expr::Variable("a".to_string()),
//!     Expr::Variable("b".to_string()),
//!     Expr::Variable("c".to_string()),
//!     Expr::Variable("d".to_string()),
//! ]);
//! ```
//!
//! ### Dynamic Operations
//!
//! ```rust
//! use rssn::symbolic::core::{Expr, register_dynamic_op, DynamicOpProperties};
//! use std::sync::Arc;
//!
//! // Register a custom operation
//! register_dynamic_op("custom_func", DynamicOpProperties {
//!     name: "custom_func".to_string(),
//!     description: "My custom function".to_string(),
//!     is_associative: false,
//!     is_commutative: false,
//! });
//!
//! // Use it
//! let expr = Expr::UnaryList(
//!     "custom_func".to_string(),
//!     Arc::new(Expr::Variable("x".to_string()))
//! );
//! ```
//!
//! ## Performance Considerations
//!
//! - **DAG representation** reduces memory usage through structural sharing
//! - **Hash-based deduplication** provides O(1) lookup for common subexpressions
//! - **N-ary operations** reduce tree depth and improve cache locality
//! - **Lazy evaluation** defers expensive operations until needed
//!
//! ## Thread Safety
//!
//! - The [`DAG_MANAGER`] uses internal locking for thread-safe access
//! - The [`DYNAMIC_OP_REGISTRY`] uses `RwLock` for concurrent reads
//! - Individual [`Expr`] values are immutable and can be shared across threads
//!
//! ## See Also
//!
//! - [`simplify_dag`](crate::symbolic::simplify_dag) - Modern DAG-based simplification
//! - [`simplify`](crate::symbolic::simplify) - Legacy AST-based simplification (deprecated)
//! - [`calculus`](crate::symbolic::calculus) - Symbolic differentiation and integration
//! - [`elementary`](crate::symbolic::elementary) - Elementary function transformations

#![allow(deprecated)]

use std::convert::AsRef;

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
use std::sync::{Arc, Mutex, RwLock};

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

    // --- N-ary Arithmetic Operations ---
    // These variants support N-ary operations, which are more efficient for
    // associative operations like addition and multiplication.
    // They are intended to eventually replace nested binary Add/Mul chains.
    /// N-ary Addition (Sum of a list of expressions).
    ///
    /// This variant represents the sum of multiple expressions in a single operation,
    /// which is more efficient than nested binary `Add` operations for associative operations.
    /// It's part of the AST to DAG migration strategy, allowing for more flexible
    /// expression representation without breaking existing code.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Representing a + b + c + d as a single n-ary operation
    /// let sum = Expr::AddList(vec![
    ///     Expr::Variable("a".to_string()),
    ///     Expr::Variable("b".to_string()),
    ///     Expr::Variable("c".to_string()),
    ///     Expr::Variable("d".to_string()),
    /// ]);
    /// ```
    AddList(Vec<Expr>),
    /// N-ary Multiplication (Product of a list of expressions).
    ///
    /// This variant represents the product of multiple expressions in a single operation,
    /// which is more efficient than nested binary `Mul` operations for associative operations.
    /// It's part of the AST to DAG migration strategy, allowing for more flexible
    /// expression representation without breaking existing code.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Representing a * b * c * d as a single n-ary operation
    /// let product = Expr::MulList(vec![
    ///     Expr::Variable("a".to_string()),
    ///     Expr::Variable("b".to_string()),
    ///     Expr::Variable("c".to_string()),
    ///     Expr::Variable("d".to_string()),
    /// ]);
    /// ```
    MulList(Vec<Expr>),

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
    Predicate { name: String, args: Vec<Expr> },
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
    ParametricSolution { x: Arc<Expr>, y: Arc<Expr> },
    /// Represents the `i`-th root of a polynomial.
    RootOf { poly: Arc<Expr>, index: u32 },
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
    ///
    /// This is now the preferred representation for all expressions.
    /// When serialized, the DAG structure is preserved.
    Dag(Arc<DagNode>),
    /// A probability distribution.
    #[serde(
        skip_serializing,
        skip_deserializing
    )]
    Distribution(Arc<dyn Distribution>),
    /// Maximum of two expressions.
    Max(Arc<Expr>, Arc<Expr>),
    /// A unified quantity with its value and unit string.
    Quantity(Arc<UnitQuantity>),
    /// A temporary representation of a value with a unit string, before unification.
    QuantityWithValue(Arc<Expr>, String),

    // --- Custom Variants (Old and Deprecated)---
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    CustomZero,
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    CustomString(String),

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    CustomArcOne(Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'BinaryList' variant instead."
    )]
    CustomArcTwo(Arc<Expr>, Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomArcThree(Arc<Expr>, Arc<Expr>, Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomArcFour(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomArcFive(Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>, Arc<Expr>),

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    CustomVecOne(Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'BinaryList' variant instead."
    )]
    CustomVecTwo(Vec<Expr>, Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomVecThree(Vec<Expr>, Vec<Expr>, Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomVecFour(Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]
    CustomVecFive(Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>, Vec<Expr>),

    // --- Dynamic/Generic Operations ---
    /// Generic unary operation identified by a name.
    ///
    /// This variant allows for extensible unary operations without modifying the `Expr` enum.
    /// Operations are registered in the `DYNAMIC_OP_REGISTRY` with their properties
    /// (associativity, commutativity, etc.). This is part of the plugin system and
    /// future-proofing strategy.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::{Expr, register_dynamic_op, DynamicOpProperties};
    /// use std::sync::Arc;
    ///
    /// // Register a custom operation
    /// register_dynamic_op("my_func", DynamicOpProperties {
    ///     name: "my_func".to_string(),
    ///     description: "My custom function".to_string(),
    ///     is_associative: false,
    ///     is_commutative: false,
    /// });
    ///
    /// // Use it
    /// let expr = Expr::UnaryList(
    ///     "my_func".to_string(),
    ///     Arc::new(Expr::Variable("x".to_string()))
    /// );
    /// ```
    UnaryList(String, Arc<Expr>),
    /// Generic binary operation identified by a name.
    ///
    /// This variant allows for extensible binary operations without modifying the `Expr` enum.
    /// Operations are registered in the `DYNAMIC_OP_REGISTRY` with their properties.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::{Expr, register_dynamic_op, DynamicOpProperties};
    /// use std::sync::Arc;
    ///
    /// // Register a custom binary operation
    /// register_dynamic_op("my_binop", DynamicOpProperties {
    ///     name: "my_binop".to_string(),
    ///     description: "My custom binary operation".to_string(),
    ///     is_associative: true,
    ///     is_commutative: true,
    /// });
    ///
    /// // Use it
    /// let expr = Expr::BinaryList(
    ///     "my_binop".to_string(),
    ///     Arc::new(Expr::Variable("x".to_string())),
    ///     Arc::new(Expr::Variable("y".to_string()))
    /// );
    /// ```
    BinaryList(String, Arc<Expr>, Arc<Expr>),
    /// Generic n-ary operation identified by a name.
    ///
    /// This variant allows for extensible n-ary operations without modifying the `Expr` enum.
    /// Operations are registered in the `DYNAMIC_OP_REGISTRY` with their properties.
    /// This is particularly useful for operations that can take a variable number of arguments.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::{Expr, register_dynamic_op, DynamicOpProperties};
    ///
    /// // Register a custom n-ary operation
    /// register_dynamic_op("my_nary", DynamicOpProperties {
    ///     name: "my_nary".to_string(),
    ///     description: "My custom n-ary operation".to_string(),
    ///     is_associative: true,
    ///     is_commutative: false,
    /// });
    ///
    /// // Use it
    /// let expr = Expr::NaryList(
    ///     "my_nary".to_string(),
    ///     vec![
    ///         Expr::Variable("a".to_string()),
    ///         Expr::Variable("b".to_string()),
    ///         Expr::Variable("c".to_string()),
    ///     ]
    /// );
    /// ```
    NaryList(String, Vec<Expr>),
}

impl Clone for Expr {
    fn clone(&self) -> Self {

        match self {
            Self::Constant(c) => Self::Constant(*c),
            Self::BigInt(i) => Self::BigInt(i.clone()),
            Self::Rational(r) => Self::Rational(r.clone()),
            Self::Boolean(b) => Self::Boolean(*b),
            Self::Variable(s) => Self::Variable(s.clone()),
            Self::Pattern(s) => Self::Pattern(s.clone()),
            Self::Add(a, b) => Self::Add(a.clone(), b.clone()),
            Self::AddList(list) => Self::AddList(list.clone()),
            Self::Sub(a, b) => Self::Sub(a.clone(), b.clone()),
            Self::Mul(a, b) => Self::Mul(a.clone(), b.clone()),
            Self::MulList(list) => Self::MulList(list.clone()),
            Self::Div(a, b) => Self::Div(a.clone(), b.clone()),
            Self::Power(a, b) => Self::Power(a.clone(), b.clone()),
            Self::Sin(a) => Self::Sin(a.clone()),
            Self::Cos(a) => Self::Cos(a.clone()),
            Self::Tan(a) => Self::Tan(a.clone()),
            Self::Exp(a) => Self::Exp(a.clone()),
            Self::Log(a) => Self::Log(a.clone()),
            Self::Neg(a) => Self::Neg(a.clone()),
            Self::Eq(a, b) => Self::Eq(a.clone(), b.clone()),
            Self::Matrix(m) => Self::Matrix(m.clone()),
            Self::Vector(v) => Self::Vector(v.clone()),
            Self::Complex(re, im) => Self::Complex(re.clone(), im.clone()),
            Self::Derivative(e, s) => Self::Derivative(e.clone(), s.clone()),
            Self::Sum {
                body,
                var,
                from,
                to,
            } => Self::Sum {
                body: body.clone(),
                var: var.clone(),
                from: from.clone(),
                to: to.clone(),
            },
            Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => Self::Integral {
                integrand: integrand.clone(),
                var: var.clone(),
                lower_bound: lower_bound.clone(),
                upper_bound: upper_bound.clone(),
            },
            Self::Path(pt, p1, p2) => Self::Path(pt.clone(), p1.clone(), p2.clone()),
            Self::Abs(a) => Self::Abs(a.clone()),
            Self::Sqrt(a) => Self::Sqrt(a.clone()),
            Self::Sec(a) => Self::Sec(a.clone()),
            Self::Csc(a) => Self::Csc(a.clone()),
            Self::Cot(a) => Self::Cot(a.clone()),
            Self::ArcSin(a) => Self::ArcSin(a.clone()),
            Self::ArcCos(a) => Self::ArcCos(a.clone()),
            Self::ArcTan(a) => Self::ArcTan(a.clone()),
            Self::ArcSec(a) => Self::ArcSec(a.clone()),
            Self::ArcCsc(a) => Self::ArcCsc(a.clone()),
            Self::ArcCot(a) => Self::ArcCot(a.clone()),
            Self::Sinh(a) => Self::Sinh(a.clone()),
            Self::Cosh(a) => Self::Cosh(a.clone()),
            Self::Tanh(a) => Self::Tanh(a.clone()),
            Self::Sech(a) => Self::Sech(a.clone()),
            Self::Csch(a) => Self::Csch(a.clone()),
            Self::Coth(a) => Self::Coth(a.clone()),
            Self::ArcSinh(a) => Self::ArcSinh(a.clone()),
            Self::ArcCosh(a) => Self::ArcCosh(a.clone()),
            Self::ArcTanh(a) => Self::ArcTanh(a.clone()),
            Self::ArcSech(a) => Self::ArcSech(a.clone()),
            Self::ArcCsch(a) => Self::ArcCsch(a.clone()),
            Self::ArcCoth(a) => Self::ArcCoth(a.clone()),
            Self::LogBase(b, a) => Self::LogBase(b.clone(), a.clone()),
            Self::Atan2(y, x) => Self::Atan2(y.clone(), x.clone()),
            Self::Binomial(n, k) => Self::Binomial(n.clone(), k.clone()),
            Self::Boundary(e) => Self::Boundary(e.clone()),
            Self::Domain(s) => Self::Domain(s.clone()),
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => Self::VolumeIntegral {
                scalar_field: scalar_field.clone(),
                volume: volume.clone(),
            },
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => Self::SurfaceIntegral {
                vector_field: vector_field.clone(),
                surface: surface.clone(),
            },
            Self::Pi => Self::Pi,
            Self::E => Self::E,
            Self::Infinity => Self::Infinity,
            Self::NegativeInfinity => Self::NegativeInfinity,
            Self::Apply(a, b) => Self::Apply(a.clone(), b.clone()),
            Self::Tuple(v) => Self::Tuple(v.clone()),
            Self::Gamma(a) => Self::Gamma(a.clone()),
            Self::Beta(a, b) => Self::Beta(a.clone(), b.clone()),
            Self::Erf(a) => Self::Erf(a.clone()),
            Self::Erfc(a) => Self::Erfc(a.clone()),
            Self::Erfi(a) => Self::Erfi(a.clone()),
            Self::Zeta(a) => Self::Zeta(a.clone()),
            Self::BesselJ(a, b) => Self::BesselJ(a.clone(), b.clone()),
            Self::BesselY(a, b) => Self::BesselY(a.clone(), b.clone()),
            Self::LegendreP(a, b) => Self::LegendreP(a.clone(), b.clone()),
            Self::LaguerreL(a, b) => Self::LaguerreL(a.clone(), b.clone()),
            Self::HermiteH(a, b) => Self::HermiteH(a.clone(), b.clone()),
            Self::Digamma(a) => Self::Digamma(a.clone()),
            Self::KroneckerDelta(a, b) => Self::KroneckerDelta(a.clone(), b.clone()),
            Self::DerivativeN(e, s, n) => Self::DerivativeN(e.clone(), s.clone(), n.clone()),
            Self::Series(a, b, c, d) => Self::Series(a.clone(), b.clone(), c.clone(), d.clone()),
            Self::Summation(a, b, c, d) => {
                Self::Summation(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::Product(a, b, c, d) => Self::Product(a.clone(), b.clone(), c.clone(), d.clone()),
            Self::ConvergenceAnalysis(e, s) => Self::ConvergenceAnalysis(e.clone(), s.clone()),
            Self::AsymptoticExpansion(a, b, c, d) => {
                Self::AsymptoticExpansion(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::Lt(a, b) => Self::Lt(a.clone(), b.clone()),
            Self::Gt(a, b) => Self::Gt(a.clone(), b.clone()),
            Self::Le(a, b) => Self::Le(a.clone(), b.clone()),
            Self::Ge(a, b) => Self::Ge(a.clone(), b.clone()),
            Self::Union(v) => Self::Union(v.clone()),
            Self::Interval(a, b, c, d) => Self::Interval(a.clone(), b.clone(), *c, *d),
            Self::Solve(e, s) => Self::Solve(e.clone(), s.clone()),
            Self::Substitute(a, b, c) => Self::Substitute(a.clone(), b.clone(), c.clone()),
            Self::Limit(a, b, c) => Self::Limit(a.clone(), b.clone(), c.clone()),
            Self::InfiniteSolutions => Self::InfiniteSolutions,
            Self::NoSolution => Self::NoSolution,
            Self::Dag(n) => Self::Dag(n.clone()),
            Self::Factorial(a) => Self::Factorial(a.clone()),
            Self::Permutation(a, b) => Self::Permutation(a.clone(), b.clone()),
            Self::Combination(a, b) => Self::Combination(a.clone(), b.clone()),
            Self::FallingFactorial(a, b) => Self::FallingFactorial(a.clone(), b.clone()),
            Self::RisingFactorial(a, b) => Self::RisingFactorial(a.clone(), b.clone()),
            Self::Ode {
                equation,
                func,
                var,
            } => Self::Ode {
                equation: equation.clone(),
                func: func.clone(),
                var: var.clone(),
            },
            Self::Pde {
                equation,
                func,
                vars,
            } => Self::Pde {
                equation: equation.clone(),
                func: func.clone(),
                vars: vars.clone(),
            },
            Self::GeneralSolution(e) => Self::GeneralSolution(e.clone()),
            Self::ParticularSolution(e) => Self::ParticularSolution(e.clone()),
            Self::Fredholm(a, b, c, d) => {
                Self::Fredholm(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::Volterra(a, b, c, d) => {
                Self::Volterra(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::And(v) => Self::And(v.clone()),
            Self::Or(v) => Self::Or(v.clone()),
            Self::Not(a) => Self::Not(a.clone()),
            Self::Xor(a, b) => Self::Xor(a.clone(), b.clone()),
            Self::Implies(a, b) => Self::Implies(a.clone(), b.clone()),
            Self::Equivalent(a, b) => Self::Equivalent(a.clone(), b.clone()),
            Self::Predicate { name, args } => Self::Predicate {
                name: name.clone(),
                args: args.clone(),
            },
            Self::ForAll(s, e) => Self::ForAll(s.clone(), e.clone()),
            Self::Exists(s, e) => Self::Exists(s.clone(), e.clone()),
            Self::Polynomial(c) => Self::Polynomial(c.clone()),
            Self::SparsePolynomial(p) => Self::SparsePolynomial(p.clone()),
            Self::Floor(a) => Self::Floor(a.clone()),
            Self::IsPrime(a) => Self::IsPrime(a.clone()),
            Self::Gcd(a, b) => Self::Gcd(a.clone(), b.clone()),
            Self::Distribution(d) => Self::Distribution(d.clone()),
            Self::Mod(a, b) => Self::Mod(a.clone(), b.clone()),
            Self::Max(a, b) => Self::Max(a.clone(), b.clone()),
            Self::Quantity(q) => Self::Quantity(q.clone()),
            Self::QuantityWithValue(v, u) => Self::QuantityWithValue(v.clone(), u.clone()),
            Self::Transpose(a) => Self::Transpose(a.clone()),
            Self::MatrixMul(a, b) => Self::MatrixMul(a.clone(), b.clone()),
            Self::MatrixVecMul(a, b) => Self::MatrixVecMul(a.clone(), b.clone()),
            Self::Inverse(a) => Self::Inverse(a.clone()),
            Self::System(v) => Self::System(v.clone()),
            Self::Solutions(v) => Self::Solutions(v.clone()),
            Self::ParametricSolution { x, y } => Self::ParametricSolution {
                x: x.clone(),
                y: y.clone(),
            },
            Self::RootOf { poly, index } => Self::RootOf {
                poly: poly.clone(),
                index: *index,
            },

            Self::CustomZero => Self::CustomZero,
            Self::CustomString(a) => Self::CustomString(a.clone()),
            Self::CustomArcOne(a) => Self::CustomArcOne(a.clone()),
            Self::CustomArcTwo(a, b) => Self::CustomArcTwo(a.clone(), b.clone()),
            Self::CustomArcThree(a, b, c) => Self::CustomArcThree(a.clone(), b.clone(), c.clone()),
            Self::CustomArcFour(a, b, c, d) => {
                Self::CustomArcFour(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::CustomArcFive(a, b, c, d, e) => {
                Self::CustomArcFive(a.clone(), b.clone(), c.clone(), d.clone(), e.clone())
            }
            Self::CustomVecOne(a) => Self::CustomVecOne(a.clone()),
            Self::CustomVecTwo(a, b) => Self::CustomVecTwo(a.clone(), b.clone()),
            Self::CustomVecThree(a, b, c) => Self::CustomVecThree(a.clone(), b.clone(), c.clone()),
            Self::CustomVecFour(a, b, c, d) => {
                Self::CustomVecFour(a.clone(), b.clone(), c.clone(), d.clone())
            }
            Self::CustomVecFive(a, b, c, d, e) => {
                Self::CustomVecFive(a.clone(), b.clone(), c.clone(), d.clone(), e.clone())
            }
            Self::UnaryList(s, a) => Self::UnaryList(s.clone(), a.clone()),
            Self::BinaryList(s, a, b) => Self::BinaryList(s.clone(), a.clone(), b.clone()),
            Self::NaryList(s, v) => Self::NaryList(s.clone(), v.clone()),
        }
    }
}

impl Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        // Use Display for a more compact representation in debug outputs
        write!(f, "{self}")
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        match self {
            Self::Dag(node) => match node.to_expr() {
                Ok(expr) => write!(f, "{expr}"),
                Err(e) => write!(f, "<Error converting DAG to Expr: {e}>"),
            },
            Self::Constant(c) => write!(f, "{c}"),
            Self::BigInt(i) => write!(f, "{i}"),
            Self::Rational(r) => write!(f, "{r}"),
            Self::Boolean(b) => write!(f, "{b}"),
            Self::Variable(s) => write!(f, "{s}"),
            Self::Pattern(s) => write!(f, "{s}"),
            Self::Add(a, b) => write!(f, "({a} + {b})"),
            Self::AddList(list) => {

                write!(f, "(")?;

                for (i, item) in list
                    .iter()
                    .enumerate()
                {

                    if i > 0 {

                        write!(f, " + ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            }
            Self::Sub(a, b) => write!(f, "({a} - {b})"),
            Self::Mul(a, b) => write!(f, "({a} * {b})"),
            Self::MulList(list) => {

                write!(f, "(")?;

                for (i, item) in list
                    .iter()
                    .enumerate()
                {

                    if i > 0 {

                        write!(f, " * ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            }
            Self::Div(a, b) => write!(f, "({a} / {b})"),
            Self::Power(a, b) => write!(f, "({a}^({b}))"),
            Self::Sin(a) => write!(f, "sin({a})"),
            Self::Cos(a) => write!(f, "cos({a})"),
            Self::Tan(a) => write!(f, "tan({a})"),
            Self::Exp(a) => write!(f, "exp({a})"),
            Self::Log(a) => write!(f, "ln({a})"),
            Self::Neg(a) => write!(f, "-({a})"),
            Self::Eq(a, b) => write!(f, "{a} = {b}"),
            Self::Matrix(m) => {

                let mut s = String::new();

                s.push('[');

                for (i, row) in m.iter().enumerate() {

                    if i > 0 {

                        s.push_str("; ");
                    }

                    s.push('[');

                    for (j, val) in row
                        .iter()
                        .enumerate()
                    {

                        if j > 0 {

                            s.push_str(", ");
                        }

                        s.push_str(&val.to_string());
                    }

                    s.push(']');
                }

                s.push(']');

                write!(f, "{s}")
            }
            Self::Vector(v) => {

                let mut s = String::new();

                s.push('[');

                for (i, val) in v.iter().enumerate() {

                    if i > 0 {

                        s.push_str(", ");
                    }

                    s.push_str(&val.to_string());
                }

                s.push(']');

                write!(f, "{s}")
            }
            Self::Complex(re, im) => write!(f, "({re} + {im}i)"),
            Self::Derivative(expr, var) => write!(f, "d/d{var}({expr})"),
            Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => write!(
                f,
                "integral({integrand}, {var}, {lower_bound}, {upper_bound})"
            ),
            Self::Sum {
                body,
                var,
                from,
                to,
            } => {

                write!(f, "sum({body}, {var}, {from}, {to})")
            }
            Self::Path(path_type, p1, p2) => write!(f, "path({path_type:?}, {p1}, {p2})"),
            Self::Abs(a) => write!(f, "|{a}|"),
            Self::Sqrt(a) => write!(f, "sqrt({a})"),
            Self::Sec(a) => write!(f, "sec({a})"),
            Self::Csc(a) => write!(f, "csc({a})"),
            Self::Cot(a) => write!(f, "cot({a})"),
            Self::ArcSin(a) => write!(f, "asin({a})"),
            Self::ArcCos(a) => write!(f, "acos({a})"),
            Self::ArcTan(a) => write!(f, "atan({a})"),
            Self::ArcSec(a) => write!(f, "asec({a})"),
            Self::ArcCsc(a) => write!(f, "acsc({a})"),
            Self::ArcCot(a) => write!(f, "acot({a})"),
            Self::Sinh(a) => write!(f, "sinh({a})"),
            Self::Cosh(a) => write!(f, "cosh({a})"),
            Self::Tanh(a) => write!(f, "tanh({a})"),
            Self::Sech(a) => write!(f, "sech({a})"),
            Self::Csch(a) => write!(f, "csch({a})"),
            Self::Coth(a) => write!(f, "coth({a})"),
            Self::ArcSinh(a) => write!(f, "asinh({a})"),
            Self::ArcCosh(a) => write!(f, "acosh({a})"),
            Self::ArcTanh(a) => write!(f, "atanh({a})"),
            Self::ArcSech(a) => write!(f, "asech({a})"),
            Self::ArcCsch(a) => write!(f, "acsch({a})"),
            Self::ArcCoth(a) => write!(f, "acoth({a})"),
            Self::LogBase(b, a) => write!(f, "log_base({b}, {a})"),
            Self::Atan2(y, x) => write!(f, "atan2({y}, {x})"),
            Self::Pi => write!(f, "Pi"),
            Self::E => write!(f, "E"),
            Self::Infinity => write!(f, "Infinity"),
            Self::NegativeInfinity => write!(f, "-Infinity"),
            Self::Ode {
                equation,
                func,
                var,
            } => write!(f, "ode({equation}, {func}, {var})"),
            Self::Pde {
                equation,
                func,
                vars,
            } => write!(f, "pde({equation}, {func}, {vars:?})"),
            Self::Fredholm(a, b, c, d) => write!(f, "fredholm({a}, {b}, {c}, {d})"),
            Self::Volterra(a, b, c, d) => write!(f, "volterra({a}, {b}, {c}, {d})"),
            Self::And(v) => write!(
                f,
                "({})",
                v.iter()
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<String>>()
                    .join(" && ")
            ),
            Self::Or(v) => write!(
                f,
                "({})",
                v.iter()
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<String>>()
                    .join(" || ")
            ),
            Self::Not(a) => write!(f, "!({a})"),
            Self::Xor(a, b) => write!(f, "({a} ^ {b})"),
            Self::Implies(a, b) => write!(f, "({a} => {b})"),
            Self::Equivalent(a, b) => write!(f, "({a} <=> {b})"),
            Self::Predicate { name, args } => {

                let args_str = args
                    .iter()
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(", ");

                write!(f, "{name}({args_str})")
            }
            Self::ForAll(s, e) => write!(f, "∀{s}. ({e})"),
            Self::Exists(s, e) => write!(f, "∃{s}. ({e})"),
            Self::Polynomial(coeffs) => {

                let mut s = String::new();

                for (i, coeff) in coeffs
                    .iter()
                    .enumerate()
                    .rev()
                {

                    if !s.is_empty() {

                        s.push_str(" + ");
                    }

                    let _ = write!(s, "{coeff}*x^{i}");
                }

                write!(f, "{s}")
            }
            Self::SparsePolynomial(p) => {

                let mut s = String::new();

                for (monomial, coeff) in &p.terms {

                    if !s.is_empty() {

                        s.push_str(" + ");
                    }

                    let _ = write!(s, "({coeff})");

                    for (var, exp) in &monomial.0 {

                        let _ = write!(s, "*{var}^{exp}");
                    }
                }

                write!(f, "{s}")
            }
            Self::Floor(a) => write!(f, "floor({a})"),
            Self::IsPrime(a) => write!(f, "is_prime({a})"),
            Self::Gcd(a, b) => write!(f, "gcd({a}, {b})"),
            Self::Factorial(a) => write!(f, "factorial({a})"),
            Self::Distribution(d) => write!(f, "{d:?}"),
            Self::Mod(a, b) => write!(f, "({a} mod {b})"),
            Self::Max(a, b) => write!(f, "max({a}, {b})"),
            Self::System(v) => write!(f, "system({v:?})"),
            Self::Solutions(v) => write!(f, "solutions({v:?})"),
            Self::ParametricSolution { x, y } => write!(f, "parametric_solution({x}, {y})"),
            Self::RootOf { poly, index } => write!(f, "root_of({poly}, {index})"),
            Self::Erfc(a) => write!(f, "erfc({a})"),
            Self::Erfi(a) => write!(f, "erfi({a})"),
            Self::Zeta(a) => write!(f, "zeta({a})"),

            Self::CustomZero => write!(f, "CustomZero"),
            Self::CustomString(s) => write!(f, "CustomString({s})"),
            Self::CustomArcOne(a) => write!(f, "CustomArcOne({a})"),
            Self::CustomArcTwo(a, b) => write!(f, "CustomArcTwo({a}, {b})"),
            Self::CustomArcThree(a, b, c) => write!(f, "CustomArcThree({a}, {b}, {c})"),
            Self::CustomArcFour(a, b, c, d) => {

                write!(f, "CustomArcFour({a}, {b}, {c}, {d})")
            }
            Self::CustomArcFive(a, b, c, d, e) => {

                write!(f, "CustomArcFive({a}, {b}, {c}, {d}, {e})")
            }
            Self::CustomVecOne(v) => write!(f, "CustomVecOne({v:?})"),
            Self::CustomVecTwo(v1, v2) => write!(f, "CustomVecTwo({v1:?}, {v2:?})"),
            Self::CustomVecThree(v1, v2, v3) => {

                write!(f, "CustomVecThree({v1:?}, {v2:?}, {v3:?})")
            }
            Self::CustomVecFour(v1, v2, v3, v4) => {

                write!(f, "CustomVecFour({v1:?}, {v2:?}, {v3:?}, {v4:?})")
            }
            Self::CustomVecFive(v1, v2, v3, v4, v5) => {

                write!(f, "CustomVecFive({v1:?}, {v2:?}, {v3:?}, {v4:?}, {v5:?})")
            }

            Self::UnaryList(s, a) => write!(f, "{s}({a})"),
            Self::BinaryList(s, a, b) => write!(f, "{s}({a}, {b})"),
            Self::NaryList(s, v) => {

                write!(f, "{s}(")?;

                for (i, item) in v.iter().enumerate() {

                    if i > 0 {

                        write!(f, ", ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            }

            Self::GeneralSolution(e) => write!(f, "general_solution({e})"),
            Self::ParticularSolution(e) => write!(f, "particular_solution({e})"),
            Self::InfiniteSolutions => write!(f, "InfiniteSolutions"),
            Self::NoSolution => write!(f, "NoSolution"),
            Self::KroneckerDelta(a, b) => write!(f, "kronecker_delta({a}, {b})"),
            Self::MatrixMul(a, b) => write!(f, "matrix_mul({a}, {b})"),
            Self::MatrixVecMul(a, b) => write!(f, "matrix_vec_mul({a}, {b})"),
            Self::QuantityWithValue(v, u) => write!(f, "quantity_with_value({v}, \"{u}\")"),
            Self::Tuple(v) => write!(f, "tuple({v:?})"),
            Self::Interval(a, b, c, d) => write!(f, "interval({a}, {b}, {c}, {d})"),
            Self::Domain(s) => write!(f, "domain({s})"),
            Self::AsymptoticExpansion(a, b, c, d) => {

                write!(f, "asymptotic_expansion({a}, {b}, {c}, {d})")
            }
            Self::ConvergenceAnalysis(e, s) => write!(f, "convergence_analysis({e}, {s})"),
            Self::DerivativeN(e, s, n) => write!(f, "derivative_n({e}, {s}, {n})"),
            Self::FallingFactorial(a, b) => write!(f, "falling_factorial({a}, {b})"),
            Self::RisingFactorial(a, b) => write!(f, "rising_factorial({a}, {b})"),
            Self::Product(a, b, c, d) => write!(f, "product({a}, {b}, {c}, {d})"),
            Self::Series(a, b, c, d) => write!(f, "series({a}, {b}, {c}, {d})"),
            Self::Summation(a, b, c, d) => write!(f, "summation({a}, {b}, {c}, {d})"),
            Self::Lt(a, b) => write!(f, "({a} < {b})"),
            Self::Gt(a, b) => write!(f, "({a} > {b})"),
            Self::Le(a, b) => write!(f, "({a} <= {b})"),
            Self::Ge(a, b) => write!(f, "({a} >= {b})"),
            Self::Transpose(a) => write!(f, "transpose({a})"),
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => write!(f, "volume_integral({scalar_field}, {volume})"),
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => write!(f, "surface_integral({vector_field}, {surface})"),
            Self::Union(v) => write!(f, "union({v:?})"),
            Self::Solve(e, s) => write!(f, "solve({e}, {s})"),
            Self::Apply(a, b) => write!(f, "apply({a}, {b})"),
            Self::Quantity(q) => write!(f, "{q:?}"),
            Self::Inverse(a) => write!(f, "inverse({a})"),
            Self::Limit(a, b, c) => write!(f, "limit({a}, {b}, {c})"),
            Self::Binomial(a, b) => write!(f, "binomial({a}, {b})"),
            Self::Permutation(a, b) => write!(f, "permutation({a}, {b})"),
            Self::Combination(a, b) => write!(f, "combination({a}, {b})"),
            Self::Boundary(a) => write!(f, "boundary({a})"),
            Self::Gamma(a) => write!(f, "gamma({a})"),
            Self::Beta(a, b) => write!(f, "beta({a}, {b})"),
            Self::Erf(a) => write!(f, "erf({a})"),
            Self::BesselJ(a, b) => write!(f, "BesselJ({a}, {b})"),
            Self::BesselY(a, b) => write!(f, "BesselY({a}, {b})"),
            Self::LegendreP(a, b) => write!(f, "LegendreP({a}, {b})"),
            Self::LaguerreL(a, b) => write!(f, "LaguerreL({a}, {b})"),
            Self::HermiteH(a, b) => write!(f, "HermiteH({a}, {b})"),
            Self::Digamma(a) => write!(f, "Digamma({a})"),
            Self::Substitute(a, b, c) => write!(f, "substitute({a}, {b}, {c})"),
        }
    }
}

impl Expr {
    #[must_use]
    #[inline]

    pub fn re(&self) -> Self {

        if let Self::Complex(re, _) = self {

            re.as_ref().clone()
        } else {

            self.clone()
        }
    }

    #[must_use]
    #[inline]

    pub fn im(&self) -> Self {

        if let Self::Complex(_, im) = self {

            im.as_ref().clone()
        } else {

            Self::Constant(0.0)
        }
    }

    #[inline]
    #[must_use]

    pub fn to_f64(&self) -> Option<f64> {

        match self {
            Self::Constant(val) => Some(*val),
            Self::BigInt(val) => val.to_f64(),
            Self::Rational(val) => val.to_f64(),
            Self::Pi => Some(std::f64::consts::PI),
            Self::E => Some(std::f64::consts::E),
            Self::Dag(node) => node
                .to_expr()
                .ok()?
                .to_f64(),
            _ => None,
        }
    }

    /// Returns the operation type of the expression in a unified way.
    ///
    /// This method handles both regular expressions and DAG nodes, returning
    /// the operation regardless of internal representation.
    ///
    /// # Returns
    /// * `DagOp` - The operation type corresponding to this expression
    ///
    #[must_use]

    pub fn op(&self) -> DagOp {

        match self {
            Self::Dag(node) => node.op.clone(),
            _ => self
                .to_dag_op_internal()
                .expect(
                    "Failed to convert Expr to DagOp; this should be impossible for any valid Expr",
                ),
        }
    }

    /// Returns the children of the expression in a unified way.
    ///
    /// This method handles both regular expressions and DAG nodes, returning
    /// the direct child expressions regardless of internal representation.
    ///
    /// # Returns
    /// * `Vec<Expr>` - A vector containing the direct children of this expression
    ///
    #[must_use]

    pub fn children(&self) -> Vec<Self> {

        match self {
            Self::Dag(node) => node
                .children
                .iter()
                .map(|n| Self::Dag(n.clone()))
                .collect(),
            _ => self.get_children_internal(),
        }
    }

    #[allow(dead_code)]

    pub(crate) const fn variant_order(&self) -> i32 {

        match self {
            Self::Constant(_) => 0,
            Self::BigInt(_) => 1,
            Self::Rational(_) => 2,
            Self::Boolean(_) => 3,
            Self::Variable(_) => 4,
            Self::Pattern(_) => 5,
            Self::Add(_, _) => 6,
            Self::AddList(_) => 6, // Same order as Add
            Self::Sub(_, _) => 7,
            Self::Mul(_, _) => 8,
            Self::MulList(_) => 8, // Same order as Mul
            Self::Div(_, _) => 9,
            Self::Power(_, _) => 10,
            Self::Sin(_) => 11,
            Self::Cos(_) => 12,
            Self::Tan(_) => 13,
            Self::Exp(_) => 14,
            Self::Log(_) => 15,
            Self::Neg(_) => 16,
            Self::Eq(_, _) => 17,
            Self::Matrix(_) => 18,
            Self::Vector(_) => 19,
            Self::Complex(_, _) => 20,
            Self::Derivative(_, _) => 21,
            Self::Integral { .. } => 22,
            Self::Sum { .. } => 22, // Assign same order as Integral for now
            Self::Path(_, _, _) => 23,
            Self::Abs(_) => 24,
            Self::Sqrt(_) => 25,
            Self::Sec(_) => 26,
            Self::Csc(_) => 27,
            Self::Cot(_) => 28,
            Self::ArcSin(_) => 29,
            Self::ArcCos(_) => 30,
            Self::ArcTan(_) => 31,
            Self::ArcSec(_) => 32,
            Self::ArcCsc(_) => 33,
            Self::ArcCot(_) => 34,
            Self::Sinh(_) => 35,
            Self::Cosh(_) => 36,
            Self::Tanh(_) => 37,
            Self::Sech(_) => 38,
            Self::Csch(_) => 39,
            Self::Coth(_) => 40,
            Self::ArcSinh(_) => 41,
            Self::ArcCosh(_) => 42,
            Self::ArcTanh(_) => 43,
            Self::ArcSech(_) => 44,
            Self::ArcCsch(_) => 45,
            Self::ArcCoth(_) => 46,
            Self::LogBase(_, _) => 47,
            Self::Atan2(_, _) => 48,
            Self::Binomial(_, _) => 49,
            Self::Boundary(_) => 50,
            Self::Domain(_) => 51,
            Self::VolumeIntegral { .. } => 52,
            Self::SurfaceIntegral { .. } => 53,
            Self::Pi => 54,
            Self::E => 55,
            Self::Infinity => 56,
            Self::NegativeInfinity => 57,
            Self::Apply(_, _) => 58,
            Self::Tuple(_) => 59,
            Self::Gamma(_) => 60,
            Self::Beta(_, _) => 61,
            Self::Erf(_) => 62,
            Self::Erfc(_) => 63,
            Self::Erfi(_) => 64,
            Self::Zeta(_) => 65,
            Self::BesselJ(_, _) => 66,
            Self::BesselY(_, _) => 67,
            Self::LegendreP(_, _) => 68,
            Self::LaguerreL(_, _) => 69,
            Self::HermiteH(_, _) => 70,
            Self::Digamma(_) => 71,
            Self::KroneckerDelta(_, _) => 72,
            Self::DerivativeN(_, _, _) => 73,
            Self::Series(_, _, _, _) => 74,
            Self::Summation(_, _, _, _) => 75,
            Self::Product(_, _, _, _) => 76,
            Self::ConvergenceAnalysis(_, _) => 77,
            Self::AsymptoticExpansion(_, _, _, _) => 78,
            Self::Lt(_, _) => 79,
            Self::Gt(_, _) => 80,
            Self::Le(_, _) => 81,
            Self::Ge(_, _) => 82,
            Self::Union(_) => 83,
            Self::Interval(_, _, _, _) => 84,
            Self::Solve(_, _) => 85,
            Self::Substitute(_, _, _) => 86,
            Self::Limit(_, _, _) => 87,
            Self::InfiniteSolutions => 88,
            Self::NoSolution => 89,
            Self::Dag(_) => 90,
            Self::Factorial(_) => 91,
            Self::Permutation(_, _) => 92,
            Self::Combination(_, _) => 93,
            Self::FallingFactorial(_, _) => 94,
            Self::RisingFactorial(_, _) => 95,
            Self::Ode { .. } => 96,
            Self::Pde { .. } => 97,
            Self::GeneralSolution(_) => 98,
            Self::ParticularSolution(_) => 99,
            Self::Fredholm(_, _, _, _) => 100,
            Self::Volterra(_, _, _, _) => 101,
            Self::And(_) => 102,
            Self::Or(_) => 103,
            Self::Not(_) => 104,
            Self::Xor(_, _) => 105,
            Self::Implies(_, _) => 106,
            Self::Equivalent(_, _) => 107,
            Self::Predicate { .. } => 108,
            Self::ForAll(_, _) => 109,
            Self::Exists(_, _) => 110,
            Self::Polynomial(_) => 111,
            Self::SparsePolynomial(_) => 112,
            Self::Floor(_) => 113,
            Self::IsPrime(_) => 114,
            Self::Gcd(_, _) => 115,
            Self::Distribution(_) => 116,
            Self::Mod(_, _) => 117,
            Self::Max(_, _) => 118,
            Self::Transpose(_) => 119,
            Self::MatrixMul(_, _) => 120,
            Self::MatrixVecMul(_, _) => 121,
            Self::Inverse(_) => 122,
            Self::System(_) => 123,
            Self::Solutions(_) => 124,
            Self::ParametricSolution { .. } => 125,
            Self::RootOf { .. } => 126,
            Self::Quantity(_) => 127,
            Self::QuantityWithValue(_, _) => 128,
            Self::CustomZero => 129,
            Self::CustomString(_) => 130,
            Self::CustomArcOne(_) => 131,
            Self::CustomArcTwo(_, _) => 132,
            Self::CustomArcThree(_, _, _) => 133,
            Self::CustomArcFour(_, _, _, _) => 134,
            Self::CustomArcFive(_, _, _, _, _) => 135,
            Self::CustomVecOne(_) => 136,
            Self::CustomVecTwo(_, _) => 137,
            Self::CustomVecThree(_, _, _) => 138,
            Self::CustomVecFour(_, _, _, _) => 139,
            Self::CustomVecFive(_, _, _, _, _) => 140,
            Self::UnaryList(_, _) => 141,
            Self::BinaryList(_, _, _) => 142,
            Self::NaryList(_, _) => 143,
        }
    }
}

#[derive(Debug, Clone, serde::Serialize)]

pub struct DagNode {
    pub op: DagOp,
    pub children: Vec<Arc<DagNode>>,
    #[serde(skip)]
    pub hash: u64,
}

impl<'de> serde::Deserialize<'de> for DagNode {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {

        #[derive(serde::Deserialize)]

        struct DagNodeHelper {
            op: DagOp,
            children: Vec<Arc<DagNode>>,
        }

        let helper = DagNodeHelper::deserialize(deserializer)?;

        // Recompute hash after deserialization
        let mut hasher = std::collections::hash_map::DefaultHasher::new();

        helper
            .op
            .hash(&mut hasher);

        for child in &helper.children {

            child
                .hash
                .hash(&mut hasher);
        }

        let hash = hasher.finish();

        Ok(Self {
            op: helper.op,
            children: helper.children,
            hash,
        })
    }
}

#[derive(
    Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]

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

    // --- Dynamic/Generic Operations ---
    UnaryList(String),
    BinaryList(String),
    NaryList(String),
}

impl PartialEq for DagNode {
    fn eq(&self, other: &Self) -> bool {

        // 1. Check the operation
        if self.op != other.op {

            return false;
        }

        // 2. Check the length
        if self.children.len() != other.children.len() {

            return false;
        }

        // 3. Recursively compare the CONTENT of the children (O(N))
        self.children
            .iter()
            .zip(
                other
                    .children
                    .iter(),
            )
            .all(|(l_child_arc, r_child_arc)| {

                // This calls PartialEq recursively on the DagNode contents
                l_child_arc
                    .as_ref()
                    .eq(r_child_arc.as_ref())
            })

        // NOTE: The hash field is IMPLICITLY IGNORED by this manual implementation.
    }
}

impl Eq for DagNode {}

impl PartialOrd for DagNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {

        Some(self.cmp(other))
    }
}

impl Ord for DagNode {
    fn cmp(&self, other: &Self) -> Ordering {

        // A stable, structural comparison for canonical sorting.
        // Compare by operation type first, then recursively by children.
        self.op
            .cmp(&other.op)
            .then_with(|| {
                self.children
                    .cmp(&other.children)
            })
    }
}

impl Hash for DagNode {
    fn hash<H: Hasher>(&self, state: &mut H) {

        self.op.hash(state);

        self.children
            .hash(state);
    }
}

impl From<DagNode> for Expr {
    fn from(node: DagNode) -> Self {

        node.to_expr()
            .expect("Cannot convert DagNode to Expr.")
    }
}

impl DagNode {
    pub fn to_expr(&self) -> Result<Expr, String> {

        use std::collections::HashMap;

        // Iterative implementation using explicit stack to prevent stack overflow
        // This uses a post-order (bottom-up) traversal strategy

        const MAX_NODES: usize = 100000;

        const MAX_CHILDREN: usize = 10000;

        // Memoization: maps node hash to its converted Expr
        let mut memo: HashMap<u64, Expr> = HashMap::new();

        // Work stack: nodes to process
        let mut work_stack: Vec<Arc<Self>> = vec![Arc::new(
            self.clone(),
        )];

        // Track which nodes we've pushed to avoid cycles
        let mut visited: HashMap<u64, bool> = HashMap::new();

        let mut nodes_processed = 0;

        while let Some(node) = work_stack.pop() {

            // Safety check: prevent processing too many nodes
            nodes_processed += 1;

            if nodes_processed > MAX_NODES {

                return Err(format!("Exceeded maximum node limit of {MAX_NODES}"));
            }

            // If already converted, skip
            if memo.contains_key(&node.hash) {

                continue;
            }

            // Safety check: limit children count
            if node.children.len() > MAX_CHILDREN {

                return Err(format!(
                    "Node has too many children ({}), exceeds limit of {}",
                    node.children.len(),
                    MAX_CHILDREN
                ));
            }

            // Check if all children are already converted
            let children_ready = node
                .children
                .iter()
                .all(|child| memo.contains_key(&child.hash));

            if children_ready {

                // All children converted, now convert this node
                let children_exprs: Vec<Expr> = node
                    .children
                    .iter()
                    .filter_map(|child| {
                        memo.get(&child.hash)
                            .cloned()
                    })
                    .collect();

                // Helper macro to create Arc from children_exprs
                macro_rules! arc {
                    ($idx:expr) => {

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
                    DagOp::Constant(c) => Expr::Constant(c.into_inner()),
                    DagOp::BigInt(i) => Expr::BigInt(i.clone()),
                    DagOp::Rational(r) => Expr::Rational(r.clone()),
                    DagOp::Boolean(b) => Expr::Boolean(*b),
                    DagOp::Variable(s) => Expr::Variable(s.clone()),
                    DagOp::Pattern(s) => Expr::Pattern(s.clone()),
                    DagOp::Domain(s) => Expr::Domain(s.clone()),
                    DagOp::Pi => Expr::Pi,
                    DagOp::E => Expr::E,
                    DagOp::Infinity => Expr::Infinity,
                    DagOp::NegativeInfinity => Expr::NegativeInfinity,
                    DagOp::InfiniteSolutions => Expr::InfiniteSolutions,
                    DagOp::NoSolution => Expr::NoSolution,

                    // --- Operators with associated data ---
                    DagOp::Derivative(s) => {

                        if children_exprs.is_empty() {

                            return Err("Derivative operator requires at least 1 child".to_string());
                        }

                        Expr::Derivative(arc!(0), s.clone())
                    }
                    DagOp::DerivativeN(s) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "DerivativeN operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::DerivativeN(arc!(0), s.clone(), arc!(1))
                    }
                    DagOp::Limit(s) => {

                        if children_exprs.len() < 2 {

                            return Err("Limit operator requires at least 2 children".to_string());
                        }

                        Expr::Limit(arc!(0), s.clone(), arc!(1))
                    }
                    DagOp::Solve(s) => {

                        if children_exprs.is_empty() {

                            return Err("Solve operator requires at least 1 child".to_string());
                        }

                        Expr::Solve(arc!(0), s.clone())
                    }
                    DagOp::ConvergenceAnalysis(s) => {

                        if children_exprs.is_empty() {

                            return Err("ConvergenceAnalysis operator requires at least 1 child"
                                .to_string());
                        }

                        Expr::ConvergenceAnalysis(arc!(0), s.clone())
                    }
                    DagOp::ForAll(s) => {

                        if children_exprs.is_empty() {

                            return Err("ForAll operator requires at least 1 child".to_string());
                        }

                        Expr::ForAll(s.clone(), arc!(0))
                    }
                    DagOp::Exists(s) => {

                        if children_exprs.is_empty() {

                            return Err("Exists operator requires at least 1 child".to_string());
                        }

                        Expr::Exists(s.clone(), arc!(0))
                    }
                    DagOp::Substitute(s) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Substitute operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Substitute(arc!(0), s.clone(), arc!(1))
                    }
                    DagOp::Ode { func, var } => {

                        if children_exprs.is_empty() {

                            return Err("Ode operator requires at least 1 child".to_string());
                        }

                        Expr::Ode {
                            equation: arc!(0),
                            func: func.clone(),
                            var: var.clone(),
                        }
                    }
                    DagOp::Pde { func, vars } => {

                        if children_exprs.is_empty() {

                            return Err("Pde operator requires at least 1 child".to_string());
                        }

                        Expr::Pde {
                            equation: arc!(0),
                            func: func.clone(),
                            vars: vars.clone(),
                        }
                    }
                    DagOp::Predicate { name } => Expr::Predicate {
                        name: name.clone(),
                        args: children_exprs.clone(),
                    },
                    DagOp::Path(pt) => {

                        if children_exprs.len() < 2 {

                            return Err("Path operator requires at least 2 children".to_string());
                        }

                        Expr::Path(pt.clone(), arc!(0), arc!(1))
                    }
                    DagOp::Interval(incl_lower, incl_upper) => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Interval operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Interval(arc!(0), arc!(1), *incl_lower, *incl_upper)
                    }
                    DagOp::RootOf { index } => {

                        if children_exprs.is_empty() {

                            return Err("RootOf operator requires at least 1 child".to_string());
                        }

                        Expr::RootOf {
                            poly: arc!(0),
                            index: *index,
                        }
                    }
                    DagOp::SparsePolynomial(p) => Expr::SparsePolynomial(p.clone()),
                    DagOp::QuantityWithValue(u) => {

                        if children_exprs.is_empty() {

                            return Err(
                                "QuantityWithValue operator requires at least 1 child".to_string()
                            );
                        }

                        Expr::QuantityWithValue(arc!(0), u.clone())
                    }

                    // --- Binary operators ---
                    DagOp::Add => {

                        if children_exprs.len() < 2 {

                            return Err("Add operator requires at least 2 children".to_string());
                        }

                        if children_exprs.len() == 2 {

                            Expr::Add(arc!(0), arc!(1))
                        } else {

                            Expr::AddList(children_exprs.clone())
                        }
                    }
                    DagOp::Sub => {

                        if children_exprs.len() < 2 {

                            return Err("Sub operator requires at least 2 children".to_string());
                        }

                        Expr::Sub(arc!(0), arc!(1))
                    }
                    DagOp::Mul => {

                        if children_exprs.len() < 2 {

                            return Err("Mul operator requires at least 2 children".to_string());
                        }

                        if children_exprs.len() == 2 {

                            Expr::Mul(arc!(0), arc!(1))
                        } else {

                            Expr::MulList(children_exprs.clone())
                        }
                    }
                    DagOp::Div => {

                        if children_exprs.len() < 2 {

                            return Err("Div operator requires at least 2 children".to_string());
                        }

                        Expr::Div(arc!(0), arc!(1))
                    }
                    DagOp::Eq => {

                        if children_exprs.len() < 2 {

                            return Err("Eq operator requires at least 2 children".to_string());
                        }

                        Expr::Eq(arc!(0), arc!(1))
                    }
                    DagOp::Lt => {

                        if children_exprs.len() < 2 {

                            return Err("Lt operator requires at least 2 children".to_string());
                        }

                        Expr::Lt(arc!(0), arc!(1))
                    }
                    DagOp::Gt => {

                        if children_exprs.len() < 2 {

                            return Err("Gt operator requires at least 2 children".to_string());
                        }

                        Expr::Gt(arc!(0), arc!(1))
                    }
                    DagOp::Le => {

                        if children_exprs.len() < 2 {

                            return Err("Le operator requires at least 2 children".to_string());
                        }

                        Expr::Le(arc!(0), arc!(1))
                    }
                    DagOp::Ge => {

                        if children_exprs.len() < 2 {

                            return Err("Ge operator requires at least 2 children".to_string());
                        }

                        Expr::Ge(arc!(0), arc!(1))
                    }
                    DagOp::LogBase => {

                        if children_exprs.len() < 2 {

                            return Err("LogBase operator requires at least 2 children".to_string());
                        }

                        Expr::LogBase(arc!(0), arc!(1))
                    }
                    DagOp::Atan2 => {

                        if children_exprs.len() < 2 {

                            return Err("Atan2 operator requires at least 2 children".to_string());
                        }

                        Expr::Atan2(arc!(0), arc!(1))
                    }
                    DagOp::Binomial => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Binomial operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Binomial(arc!(0), arc!(1))
                    }
                    DagOp::Beta => {

                        if children_exprs.len() < 2 {

                            return Err("Beta operator requires at least 2 children".to_string());
                        }

                        Expr::Beta(arc!(0), arc!(1))
                    }
                    DagOp::BesselJ => {

                        if children_exprs.len() < 2 {

                            return Err("BesselJ operator requires at least 2 children".to_string());
                        }

                        Expr::BesselJ(arc!(0), arc!(1))
                    }
                    DagOp::BesselY => {

                        if children_exprs.len() < 2 {

                            return Err("BesselY operator requires at least 2 children".to_string());
                        }

                        Expr::BesselY(arc!(0), arc!(1))
                    }
                    DagOp::LegendreP => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "LegendreP operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::LegendreP(arc!(0), arc!(1))
                    }
                    DagOp::LaguerreL => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "LaguerreL operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::LaguerreL(arc!(0), arc!(1))
                    }
                    DagOp::HermiteH => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "HermiteH operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::HermiteH(arc!(0), arc!(1))
                    }
                    DagOp::KroneckerDelta => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "KroneckerDelta operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::KroneckerDelta(arc!(0), arc!(1))
                    }
                    DagOp::Permutation => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Permutation operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Permutation(arc!(0), arc!(1))
                    }
                    DagOp::Combination => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Combination operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Combination(arc!(0), arc!(1))
                    }
                    DagOp::FallingFactorial => {

                        if children_exprs.len() < 2 {

                            return Err("FallingFactorial operator requires at least 2 children"
                                .to_string());
                        }

                        Expr::FallingFactorial(arc!(0), arc!(1))
                    }
                    DagOp::RisingFactorial => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "RisingFactorial operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::RisingFactorial(arc!(0), arc!(1))
                    }
                    DagOp::Xor => {

                        if children_exprs.len() < 2 {

                            return Err("Xor operator requires at least 2 children".to_string());
                        }

                        Expr::Xor(arc!(0), arc!(1))
                    }
                    DagOp::Implies => {

                        if children_exprs.len() < 2 {

                            return Err("Implies operator requires at least 2 children".to_string());
                        }

                        Expr::Implies(arc!(0), arc!(1))
                    }
                    DagOp::Equivalent => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "Equivalent operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::Equivalent(arc!(0), arc!(1))
                    }
                    DagOp::Gcd => {

                        if children_exprs.len() < 2 {

                            return Err("Gcd operator requires at least 2 children".to_string());
                        }

                        Expr::Gcd(arc!(0), arc!(1))
                    }
                    DagOp::Mod => {

                        if children_exprs.len() < 2 {

                            return Err("Mod operator requires at least 2 children".to_string());
                        }

                        Expr::Mod(arc!(0), arc!(1))
                    }
                    DagOp::Max => {

                        if children_exprs.len() < 2 {

                            return Err("Max operator requires at least 2 children".to_string());
                        }

                        Expr::Max(arc!(0), arc!(1))
                    }
                    DagOp::MatrixMul => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "MatrixMul operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::MatrixMul(arc!(0), arc!(1))
                    }
                    DagOp::MatrixVecMul => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "MatrixVecMul operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::MatrixVecMul(arc!(0), arc!(1))
                    }
                    DagOp::Apply => {

                        if children_exprs.len() < 2 {

                            return Err("Apply operator requires at least 2 children".to_string());
                        }

                        Expr::Apply(arc!(0), arc!(1))
                    }

                    // --- Unary operators ---
                    DagOp::Neg => {

                        if children_exprs.is_empty() {

                            return Err("Neg operator requires at least 1 child".to_string());
                        }

                        Expr::Neg(arc!(0))
                    }
                    DagOp::Power => {

                        if children_exprs.len() < 2 {

                            return Err("Power operator requires at least 2 children".to_string());
                        }

                        Expr::Power(arc!(0), arc!(1))
                    }
                    DagOp::Sin => {

                        if children_exprs.is_empty() {

                            return Err("Sin operator requires at least 1 child".to_string());
                        }

                        Expr::Sin(arc!(0))
                    }
                    DagOp::Cos => {

                        if children_exprs.is_empty() {

                            return Err("Cos operator requires at least 1 child".to_string());
                        }

                        Expr::Cos(arc!(0))
                    }
                    DagOp::Tan => {

                        if children_exprs.is_empty() {

                            return Err("Tan operator requires at least 1 child".to_string());
                        }

                        Expr::Tan(arc!(0))
                    }
                    DagOp::Exp => {

                        if children_exprs.is_empty() {

                            return Err("Exp operator requires at least 1 child".to_string());
                        }

                        Expr::Exp(arc!(0))
                    }
                    DagOp::Log => {

                        if children_exprs.is_empty() {

                            return Err("Log operator requires at least 1 child".to_string());
                        }

                        Expr::Log(arc!(0))
                    }
                    DagOp::Abs => {

                        if children_exprs.is_empty() {

                            return Err("Abs operator requires at least 1 child".to_string());
                        }

                        Expr::Abs(arc!(0))
                    }
                    DagOp::Sqrt => {

                        if children_exprs.is_empty() {

                            return Err("Sqrt operator requires at least 1 child".to_string());
                        }

                        Expr::Sqrt(arc!(0))
                    }
                    DagOp::Transpose => {

                        if children_exprs.is_empty() {

                            return Err("Transpose operator requires at least 1 child".to_string());
                        }

                        Expr::Transpose(arc!(0))
                    }
                    DagOp::Inverse => {

                        if children_exprs.is_empty() {

                            return Err("Inverse operator requires at least 1 child".to_string());
                        }

                        Expr::Inverse(arc!(0))
                    }
                    DagOp::Sec => {

                        if children_exprs.is_empty() {

                            return Err("Sec operator requires at least 1 child".to_string());
                        }

                        Expr::Sec(arc!(0))
                    }
                    DagOp::Csc => {

                        if children_exprs.is_empty() {

                            return Err("Csc operator requires at least 1 child".to_string());
                        }

                        Expr::Csc(arc!(0))
                    }
                    DagOp::Cot => {

                        if children_exprs.is_empty() {

                            return Err("Cot operator requires at least 1 child".to_string());
                        }

                        Expr::Cot(arc!(0))
                    }
                    DagOp::ArcSin => {

                        if children_exprs.is_empty() {

                            return Err("ArcSin operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSin(arc!(0))
                    }
                    DagOp::ArcCos => {

                        if children_exprs.is_empty() {

                            return Err("ArcCos operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCos(arc!(0))
                    }
                    DagOp::ArcTan => {

                        if children_exprs.is_empty() {

                            return Err("ArcTan operator requires at least 1 child".to_string());
                        }

                        Expr::ArcTan(arc!(0))
                    }
                    DagOp::ArcSec => {

                        if children_exprs.is_empty() {

                            return Err("ArcSec operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSec(arc!(0))
                    }
                    DagOp::ArcCsc => {

                        if children_exprs.is_empty() {

                            return Err("ArcCsc operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCsc(arc!(0))
                    }
                    DagOp::ArcCot => {

                        if children_exprs.is_empty() {

                            return Err("ArcCot operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCot(arc!(0))
                    }
                    DagOp::Sinh => {

                        if children_exprs.is_empty() {

                            return Err("Sinh operator requires at least 1 child".to_string());
                        }

                        Expr::Sinh(arc!(0))
                    }
                    DagOp::Cosh => {

                        if children_exprs.is_empty() {

                            return Err("Cosh operator requires at least 1 child".to_string());
                        }

                        Expr::Cosh(arc!(0))
                    }
                    DagOp::Tanh => {

                        if children_exprs.is_empty() {

                            return Err("Tanh operator requires at least 1 child".to_string());
                        }

                        Expr::Tanh(arc!(0))
                    }
                    DagOp::Sech => {

                        if children_exprs.is_empty() {

                            return Err("Sech operator requires at least 1 child".to_string());
                        }

                        Expr::Sech(arc!(0))
                    }
                    DagOp::Csch => {

                        if children_exprs.is_empty() {

                            return Err("Csch operator requires at least 1 child".to_string());
                        }

                        Expr::Csch(arc!(0))
                    }
                    DagOp::Coth => {

                        if children_exprs.is_empty() {

                            return Err("Coth operator requires at least 1 child".to_string());
                        }

                        Expr::Coth(arc!(0))
                    }
                    DagOp::ArcSinh => {

                        if children_exprs.is_empty() {

                            return Err("ArcSinh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSinh(arc!(0))
                    }
                    DagOp::ArcCosh => {

                        if children_exprs.is_empty() {

                            return Err("ArcCosh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCosh(arc!(0))
                    }
                    DagOp::ArcTanh => {

                        if children_exprs.is_empty() {

                            return Err("ArcTanh operator requires at least 1 child".to_string());
                        }

                        Expr::ArcTanh(arc!(0))
                    }
                    DagOp::ArcSech => {

                        if children_exprs.is_empty() {

                            return Err("ArcSech operator requires at least 1 child".to_string());
                        }

                        Expr::ArcSech(arc!(0))
                    }
                    DagOp::ArcCsch => {

                        if children_exprs.is_empty() {

                            return Err("ArcCsch operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCsch(arc!(0))
                    }
                    DagOp::ArcCoth => {

                        if children_exprs.is_empty() {

                            return Err("ArcCoth operator requires at least 1 child".to_string());
                        }

                        Expr::ArcCoth(arc!(0))
                    }
                    DagOp::Factorial => {

                        if children_exprs.is_empty() {

                            return Err("Factorial operator requires at least 1 child".to_string());
                        }

                        Expr::Factorial(arc!(0))
                    }
                    DagOp::Boundary => {

                        if children_exprs.is_empty() {

                            return Err("Boundary operator requires at least 1 child".to_string());
                        }

                        Expr::Boundary(arc!(0))
                    }
                    DagOp::Gamma => {

                        if children_exprs.is_empty() {

                            return Err("Gamma operator requires at least 1 child".to_string());
                        }

                        Expr::Gamma(arc!(0))
                    }
                    DagOp::Erf => {

                        if children_exprs.is_empty() {

                            return Err("Erf operator requires at least 1 child".to_string());
                        }

                        Expr::Erf(arc!(0))
                    }
                    DagOp::Erfc => {

                        if children_exprs.is_empty() {

                            return Err("Erfc operator requires at least 1 child".to_string());
                        }

                        Expr::Erfc(arc!(0))
                    }
                    DagOp::Erfi => {

                        if children_exprs.is_empty() {

                            return Err("Erfi operator requires at least 1 child".to_string());
                        }

                        Expr::Erfi(arc!(0))
                    }
                    DagOp::Zeta => {

                        if children_exprs.is_empty() {

                            return Err("Zeta operator requires at least 1 child".to_string());
                        }

                        Expr::Zeta(arc!(0))
                    }
                    DagOp::Digamma => {

                        if children_exprs.is_empty() {

                            return Err("Digamma operator requires at least 1 child".to_string());
                        }

                        Expr::Digamma(arc!(0))
                    }
                    DagOp::Not => {

                        if children_exprs.is_empty() {

                            return Err("Not operator requires at least 1 child".to_string());
                        }

                        Expr::Not(arc!(0))
                    }
                    DagOp::Floor => {

                        if children_exprs.is_empty() {

                            return Err("Floor operator requires at least 1 child".to_string());
                        }

                        Expr::Floor(arc!(0))
                    }
                    DagOp::IsPrime => {

                        if children_exprs.is_empty() {

                            return Err("IsPrime operator requires at least 1 child".to_string());
                        }

                        Expr::IsPrime(arc!(0))
                    }
                    DagOp::GeneralSolution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "GeneralSolution operator requires at least 1 child".to_string()
                            );
                        }

                        Expr::GeneralSolution(arc!(0))
                    }
                    DagOp::ParticularSolution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "ParticularSolution operator requires at least 1 child".to_string()
                            );
                        }

                        Expr::ParticularSolution(arc!(0))
                    }

                    // --- Complex structures ---
                    DagOp::Matrix { rows: _, cols } => {

                        if !children_exprs
                            .len()
                            .is_multiple_of(*cols)
                        {

                            let complete_rows = (children_exprs.len() / cols) * cols;

                            let reconstructed_matrix: Vec<Vec<Expr>> = children_exprs
                                .iter()
                                .take(complete_rows)
                                .cloned()
                                .collect::<Vec<_>>()
                                .chunks(*cols)
                                .map(<[Expr]>::to_vec)
                                .collect();

                            Expr::Matrix(reconstructed_matrix)
                        } else {

                            let reconstructed_matrix: Vec<Vec<Expr>> = children_exprs
                                .chunks(*cols)
                                .map(<[Expr]>::to_vec)
                                .collect();

                            Expr::Matrix(reconstructed_matrix)
                        }
                    }
                    DagOp::Complex => {

                        if children_exprs.len() < 2 {

                            return Err("Complex operator requires at least 2 children".to_string());
                        }

                        Expr::Complex(arc!(0), arc!(1))
                    }
                    DagOp::Integral => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Integral operator requires at least 4 children".to_string()
                            );
                        }

                        Expr::Integral {
                            integrand: arc!(0),
                            var: arc!(1),
                            lower_bound: arc!(2),
                            upper_bound: arc!(3),
                        }
                    }
                    DagOp::VolumeIntegral => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "VolumeIntegral operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::VolumeIntegral {
                            scalar_field: arc!(0),
                            volume: arc!(1),
                        }
                    }
                    DagOp::SurfaceIntegral => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "SurfaceIntegral operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::SurfaceIntegral {
                            vector_field: arc!(0),
                            surface: arc!(1),
                        }
                    }
                    DagOp::Sum => {

                        if children_exprs.len() < 4 {

                            return Err("Sum operator requires at least 4 children".to_string());
                        }

                        Expr::Sum {
                            body: arc!(0),
                            var: arc!(1),
                            from: arc!(2),
                            to: arc!(3),
                        }
                    }
                    DagOp::Series(s) => {

                        if children_exprs.len() < 3 {

                            return Err("Series operator requires at least 3 children".to_string());
                        }

                        Expr::Series(arc!(0), s.clone(), arc!(1), arc!(2))
                    }
                    DagOp::Summation(s) => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "Summation operator requires at least 3 children".to_string()
                            );
                        }

                        Expr::Summation(arc!(0), s.clone(), arc!(1), arc!(2))
                    }
                    DagOp::Product(s) => {

                        if children_exprs.len() < 3 {

                            return Err("Product operator requires at least 3 children".to_string());
                        }

                        Expr::Product(arc!(0), s.clone(), arc!(1), arc!(2))
                    }
                    DagOp::AsymptoticExpansion(s) => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "AsymptoticExpansion operator requires at least 3 children"
                                    .to_string(),
                            );
                        }

                        Expr::AsymptoticExpansion(arc!(0), s.clone(), arc!(1), arc!(2))
                    }
                    DagOp::ParametricSolution => {

                        if children_exprs.len() < 2 {

                            return Err("ParametricSolution operator requires at least 2 children"
                                .to_string());
                        }

                        Expr::ParametricSolution {
                            x: arc!(0),
                            y: arc!(1),
                        }
                    }
                    DagOp::Fredholm => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Fredholm operator requires at least 4 children".to_string()
                            );
                        }

                        Expr::Fredholm(arc!(0), arc!(1), arc!(2), arc!(3))
                    }
                    DagOp::Volterra => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "Volterra operator requires at least 4 children".to_string()
                            );
                        }

                        Expr::Volterra(arc!(0), arc!(1), arc!(2), arc!(3))
                    }
                    DagOp::Distribution => {

                        if children_exprs.is_empty() {

                            return Err(
                                "Distribution operator requires at least 1 child".to_string()
                            );
                        }

                        Expr::Distribution(children_exprs[0].clone_box_dist()?)
                    }
                    DagOp::Quantity => {

                        if children_exprs.is_empty() {

                            return Err("Quantity operator requires at least 1 child".to_string());
                        }

                        Expr::Quantity(children_exprs[0].clone_box_quant()?)
                    }

                    // --- List operators ---
                    DagOp::Vector => Expr::Vector(children_exprs.clone()),
                    DagOp::And => Expr::And(children_exprs.clone()),
                    DagOp::Or => Expr::Or(children_exprs.clone()),
                    DagOp::Union => Expr::Union(children_exprs.clone()),
                    DagOp::Polynomial => Expr::Polynomial(children_exprs.clone()),
                    DagOp::System => Expr::System(children_exprs.clone()),
                    DagOp::Solutions => Expr::Solutions(children_exprs.clone()),
                    DagOp::Tuple => Expr::Tuple(children_exprs.clone()),

                    // --- Custom ---
                    DagOp::CustomZero => Expr::CustomZero,
                    DagOp::CustomString(s) => Expr::CustomString(s.clone()),
                    DagOp::CustomArcOne => {

                        if children_exprs.is_empty() {

                            return Err(
                                "CustomArcOne operator requires at least 1 child".to_string()
                            );
                        }

                        Expr::CustomArcOne(arc!(0))
                    }
                    DagOp::CustomArcTwo => {

                        if children_exprs.len() < 2 {

                            return Err(
                                "CustomArcTwo operator requires at least 2 children".to_string()
                            );
                        }

                        Expr::CustomArcTwo(arc!(0), arc!(1))
                    }
                    DagOp::CustomArcThree => {

                        if children_exprs.len() < 3 {

                            return Err(
                                "CustomArcThree operator requires at least 3 children".to_string()
                            );
                        }

                        Expr::CustomArcThree(arc!(0), arc!(1), arc!(2))
                    }
                    DagOp::CustomArcFour => {

                        if children_exprs.len() < 4 {

                            return Err(
                                "CustomArcFour operator requires at least 4 children".to_string()
                            );
                        }

                        Expr::CustomArcFour(arc!(0), arc!(1), arc!(2), arc!(3))
                    }
                    DagOp::CustomArcFive => {

                        if children_exprs.len() < 5 {

                            return Err(
                                "CustomArcFive operator requires at least 5 children".to_string()
                            );
                        }

                        Expr::CustomArcFive(arc!(0), arc!(1), arc!(2), arc!(3), arc!(4))
                    }
                    DagOp::CustomVecOne => Expr::CustomVecOne(children_exprs.clone()),
                    DagOp::CustomVecTwo => {
                        return Err("CustomVecTwo to_expr is ambiguous".to_string())
                    }
                    DagOp::CustomVecThree => {
                        return Err("CustomVecThree to_expr is ambiguous".to_string())
                    }
                    DagOp::CustomVecFour => {
                        return Err("CustomVecFour to_expr is ambiguous".to_string())
                    }
                    DagOp::CustomVecFive => {
                        return Err("CustomVecFive to_expr is ambiguous".to_string())
                    }

                    DagOp::UnaryList(s) => {

                        if children_exprs.is_empty() {

                            return Err(format!(
                                "UnaryList operator {s} requires at least 1 child"
                            ));
                        }

                        Expr::UnaryList(s.clone(), arc!(0))
                    }
                    DagOp::BinaryList(s) => {

                        if children_exprs.len() < 2 {

                            return Err(format!(
                                "BinaryList operator {s} requires at least 2 children"
                            ));
                        }

                        Expr::BinaryList(s.clone(), arc!(0), arc!(1))
                    }
                    DagOp::NaryList(s) => Expr::NaryList(s.clone(), children_exprs.clone()),
                };

                // Store the converted expression
                memo.insert(node.hash, expr);
            } else {

                // Not all children ready, push node back and push children
                work_stack.push(node.clone());

                // Push children in reverse order (so they're processed in correct order)
                if visited
                    .insert(node.hash, true)
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
            .ok_or_else(|| "Failed to convert root node".to_string())
    }

    #[must_use]

    pub fn new(op: DagOp, children: Vec<Arc<Self>>) -> Arc<Self> {

        // Safety check: limit number of children to prevent excessive memory allocation
        const MAX_CHILDREN: usize = 10000;

        if children.len() > MAX_CHILDREN {

            // This should not happen in normal usage, but we handle it gracefully
            // by truncating the children list - this is a defensive programming approach
            let safe_children: Vec<_> = children
                .into_iter()
                .take(MAX_CHILDREN)
                .collect();

            let mut hasher = std::collections::hash_map::DefaultHasher::new();

            op.hash(&mut hasher);

            safe_children.hash(&mut hasher);

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

        Arc::new(Self { op, children, hash })
    }
}

impl Expr {
    pub fn clone_box_dist(&self) -> Result<Arc<dyn Distribution>, String> {

        if let Self::Distribution(d) = self {

            Ok(d.clone_box())
        } else {

            Err("Cannot clone into Distribution".to_string())
        }
    }

    pub fn clone_box_quant(&self) -> Result<Arc<UnitQuantity>, String> {

        if let Self::Quantity(q) = self {

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
    #[must_use]

    pub fn new() -> Self {

        Self {
            nodes: Mutex::new(HashMap::new()),
        }
    }

    /// Get an existing node identical to (op, children) if present; otherwise create and insert.
    ///
    /// This implementation avoids returning a node solely based on the `u64` hash.
    /// When a hash bucket is found, we iterate the bucket and compare structural equality
    /// (op + children count + children's hashes). Only when no equal node is found do we insert.
    ///
    /// # Safety
    /// Limits the number of children to prevent excessive memory allocation
    /// and potential stack overflows during recursive operations.
    #[inline]

    pub fn get_or_create_normalized(
        &self,
        op: DagOp,
        mut children: Vec<Arc<DagNode>>,
    ) -> Result<Arc<DagNode>, String> {

        // Safety check: limit number of children to prevent excessive memory usage
        const MAX_CHILDREN: usize = 10000;

        if children.len() > MAX_CHILDREN {

            return Err(format!(
                "Too many children in node ({}), exceeds limit of {}",
                children.len(),
                MAX_CHILDREN
            ));
        }

        match op {
            DagOp::Add | DagOp::Mul => {

                // Use stable sort to ensure deterministic ordering across test runs.
                // This is critical for reproducible hashing and test stability.
                children.sort();
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

        // Prevent excessive memory usage by limiting bucket size
        const MAX_BUCKET_SIZE: usize = 1000;

        // Ensure the bucket is a vector of candidates to support collision buckets.
        // nodes: HashMap<u64, Vec<Arc<DagNode>>>
        match nodes_guard.entry(hash) {
            Entry::Occupied(mut occ) => {

                // occ.get_mut() is a Vec<Arc<DagNode>>
                let bucket = occ.get_mut();

                // Check bucket size to prevent excessive memory usage
                if bucket.len() > MAX_BUCKET_SIZE {

                    // If the bucket is too large, we just create a new node without searching
                    // This maintains correctness while limiting memory usage
                    let node = Arc::new(DagNode { op, children, hash });

                    return Ok(node);
                }

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
        for (a, b) in cand
            .children
            .iter()
            .zip(children.iter())
        {

            if a.hash != b.hash {

                return false;
            }
        }

        // Hash equality is not enough; we must verify structural equality.
        // Since we already checked length and hashes, we now do a full check.
        // Note: This relies on DagNode::eq which performs deep comparison.
        cand.children == *children
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

    /// Gets an existing DAG node for the expression if it exists, or creates a new one.
    ///
    /// If the expression is already a DAG node, returns it directly.
    /// Otherwise, constructs the DAG representation by converting the expression
    /// to its operation and children, then creating the appropriate DAG node.
    ///
    /// This method uses the shared DAG manager to ensure canonical representation
    /// and avoid duplicate nodes for identical expressions.
    ///
    /// # Arguments
    /// * `expr` - The expression to convert to or retrieve as a DAG node
    ///
    /// # Returns
    /// * `Result<Arc<DagNode>, String>` - The DAG node representation of the expression or an error
    ///
    #[inline]

    pub fn get_or_create(&self, expr: &Expr) -> Result<Arc<DagNode>, String> {

        if let Expr::Dag(node) = expr {

            return Ok(node.clone());
        }

        // Safety check: limit recursion depth to prevent stack overflow
        // We can't directly implement a recursion depth counter here since this is a single function,
        // but we can add checks for extremely complex expressions
        let op = expr.to_dag_op_internal()?;

        let children_exprs = expr.get_children_internal();

        // Limit the number of children to prevent excessive memory allocation
        const MAX_CHILDREN_PER_NODE: usize = 10000;

        if children_exprs.len() > MAX_CHILDREN_PER_NODE {

            return Err(format!(
                "Expression has too many children ({}), exceeds limit of {}",
                children_exprs.len(),
                MAX_CHILDREN_PER_NODE
            ));
        }

        let children_nodes = children_exprs
            .iter()
            .map(|child| self.get_or_create(child))
            .collect::<Result<Vec<_>, _>>()?;

        self.get_or_create_normalized(op, children_nodes)
    }
}

impl PartialEq for Expr {
    fn eq(&self, other: &Self) -> bool {

        if let (Self::Dag(n1), Self::Dag(n2)) = (self, other) {

            if Arc::ptr_eq(n1, n2) {

                return true;
            }
        }

        if self.op() != other.op() {

            return false;
        }

        match (self, other) {
            // --- COMMUTATIVE OPERATORS (A+B == B+A) ---
            (Self::Add(l1, r1), Self::Add(l2, r2)) | (Self::Mul(l1, r1), Self::Mul(l2, r2)) => {

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
            }

            // --- NON-COMMUTATIVE OPERATORS (A-B != B-A) ---
            (Self::Sub(l1, r1), Self::Sub(l2, r2))
            | (Self::Div(l1, r1), Self::Div(l2, r2))
            | (Self::Power(l1, r1), Self::Power(l2, r2)) => {

                // Positional comparison is required for non-commutative ops
                return l1
                    .as_ref()
                    .eq(l2.as_ref())
                    && r1
                        .as_ref()
                        .eq(r2.as_ref());
            }

            // Special handling for Derivative to compare both expression and variable
            (Self::Derivative(e1, v1), Self::Derivative(e2, v2)) => {
                return v1 == v2
                    && e1
                        .as_ref()
                        .eq(e2.as_ref())
            }

            // Special handling for other variants with String parameters
            (Self::Solve(e1, v1), Self::Solve(e2, v2))
            | (Self::ConvergenceAnalysis(e1, v1), Self::ConvergenceAnalysis(e2, v2))
            | (Self::ForAll(v1, e1), Self::ForAll(v2, e2))
            | (Self::Exists(v1, e1), Self::Exists(v2, e2)) => {
                return v1 == v2
                    && e1
                        .as_ref()
                        .eq(e2.as_ref())
            }

            (Self::Constant(f1), Self::Constant(f2)) => return (f1 - f2).abs() < f64::EPSILON,
            (Self::BigInt(b1), Self::BigInt(b2)) => return b1 == b2,
            (Self::Rational(r1), Self::Rational(r2)) => return r1 == r2,

            // BigInt <=> Rational
            (Self::BigInt(b), Self::Rational(r)) | (Self::Rational(r), Self::BigInt(b)) => {

                let temp_rational = BigRational::from(b.clone());

                return r == &temp_rational;
            }

            // BigInt / Rational <=> Constant(f64)
            (Self::Constant(f), Self::Rational(r)) | (Self::Rational(r), Self::Constant(f)) => {
                match r.to_f64() {
                    Some(r_f64) => return (f - r_f64).abs() < f64::EPSILON,
                    None => return false,
                }
            }

            (Self::Constant(f), Self::BigInt(b)) | (Self::BigInt(b), Self::Constant(f)) => {

                if f.fract().abs() < f64::EPSILON {

                    match b.to_f64() {
                        Some(b_f64) => return (f - b_f64).abs() < f64::EPSILON,
                        None => return false,
                    }
                }

                return false;
            }

            _ => { /* Ignore, Enter Next Step */ }
        }

        let self_children = self.children();

        let other_children = other.children();

        if self_children.len() != other_children.len() {

            return false;
        }

        self_children
            .iter()
            .zip(other_children.iter())
            .all(|(l_child_expr, r_child_expr)| l_child_expr.eq(r_child_expr))
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

        // Fast path for identical DAG nodes.
        if let (Self::Dag(n1), Self::Dag(n2)) = (self, other) {

            if Arc::ptr_eq(n1, n2) {

                return Ordering::Equal;
            }
        }

        // Compare by operator.
        let op_ordering = self
            .op()
            .cmp(&other.op());

        if op_ordering != Ordering::Equal {

            return op_ordering;
        }

        // For canonical DAG nodes, children are already sorted, so we can compare directly.
        // For non-DAG nodes or mixed comparisons, this relies on the slow sorting path in the
        // Ord implementation of the children expressions.
        self.children()
            .cmp(&other.children())
    }
}

#[derive(Debug)]

pub enum SymbolicError {
    Msg(String),
}

impl fmt::Display for SymbolicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        match self {
            Self::Msg(s) => write!(f, "{s}"),
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
    ///

    pub fn pre_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Self),
    {

        f(self); // Visit parent
        match self {
            // Binary operators
            Self::Add(a, b)
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
            }
            // Unary operators
            Self::Sin(a)
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
            }
            // N-ary operators
            Self::Matrix(m) => m
                .iter()
                .flatten()
                .for_each(|e| e.pre_order_walk(f)),
            Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v)
            | Self::AddList(v)
            | Self::MulList(v) => v
                .iter()
                .for_each(|e| e.pre_order_walk(f)),
            Self::Predicate { args, .. } => args
                .iter()
                .for_each(|e| e.pre_order_walk(f)),
            Self::SparsePolynomial(p) => p
                .terms
                .values()
                .for_each(|c| c.pre_order_walk(f)),
            // More complex operators
            Self::Sum {
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
            Self::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {

                integrand.pre_order_walk(f);

                lower_bound.pre_order_walk(f);

                upper_bound.pre_order_walk(f);
            }
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                scalar_field.pre_order_walk(f);

                volume.pre_order_walk(f);
            }
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                vector_field.pre_order_walk(f);

                surface.pre_order_walk(f);
            }
            Self::DerivativeN(e, _, n) => {

                e.pre_order_walk(f);

                n.pre_order_walk(f);
            }
            Self::Series(a, _, c, d) | Self::Summation(a, _, c, d) | Self::Product(a, _, c, d) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            }
            Self::AsymptoticExpansion(a, _, c, d) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            }
            Self::Interval(a, b, _, _) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            }
            Self::Substitute(a, _, c) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);
            }
            Self::Limit(a, _, c) => {

                a.pre_order_walk(f);

                c.pre_order_walk(f);
            }
            Self::Ode { equation, .. } => equation.pre_order_walk(f),
            Self::Pde { equation, .. } => equation.pre_order_walk(f),
            Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            }
            Self::ParametricSolution { x, y } => {

                x.pre_order_walk(f);

                y.pre_order_walk(f);
            }
            Self::RootOf { poly, .. } => poly.pre_order_walk(f),
            Self::QuantityWithValue(v, _) => v.pre_order_walk(f),

            Self::CustomArcOne(a) => {

                a.pre_order_walk(f);
            }
            Self::CustomArcTwo(a, b) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            }
            Self::CustomArcThree(a, b, c) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);
            }
            Self::CustomArcFour(a, b, c, d) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);
            }
            Self::CustomArcFive(a, b, c, d, e) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);

                c.pre_order_walk(f);

                d.pre_order_walk(f);

                e.pre_order_walk(f);
            }
            Self::CustomVecOne(v) => v
                .iter()
                .for_each(|e| e.pre_order_walk(f)),
            Self::CustomVecTwo(v1, v2) => {

                for e in v1 {

                    e.pre_order_walk(f);
                }

                for e in v2 {

                    e.pre_order_walk(f);
                }
            }
            Self::CustomVecThree(v1, v2, v3) => {

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
            Self::CustomVecFour(v1, v2, v3, v4) => {

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
            Self::CustomVecFive(v1, v2, v3, v4, v5) => {

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
            Self::UnaryList(_, a) => a.pre_order_walk(f),
            Self::BinaryList(_, a, b) => {

                a.pre_order_walk(f);

                b.pre_order_walk(f);
            }
            Self::NaryList(_, v) => {
                for e in v {

                    e.pre_order_walk(f);
                }
            }

            Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.pre_order_walk(f);
                }
            }
            // Leaf nodes
            Self::Constant(_)
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
            | Self::Distribution(_) => {}
        }
    }

    /// Performs a post-order traversal of the expression tree.
    /// It visits the children first, then the current node.
    ///
    /// # Arguments
    /// * `f` - A mutable function that takes a reference to an `Expr` and is applied to each node during traversal
    ///

    pub fn post_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Self),
    {

        match self {
            // Binary operators
            Self::Add(a, b)
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
            }
            // Unary operators
            Self::Sin(a)
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
            }
            // N-ary operators
            Self::Matrix(m) => m
                .iter()
                .flatten()
                .for_each(|e| e.post_order_walk(f)),
            Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v)
            | Self::AddList(v)
            | Self::MulList(v) => v
                .iter()
                .for_each(|e| e.post_order_walk(f)),
            Self::Predicate { args, .. } => args
                .iter()
                .for_each(|e| e.post_order_walk(f)),
            Self::SparsePolynomial(p) => p
                .terms
                .values()
                .for_each(|c| c.post_order_walk(f)),
            // More complex operators
            Self::Integral {
                integrand,
                var: _,
                lower_bound,
                upper_bound,
            } => {

                integrand.post_order_walk(f);

                lower_bound.post_order_walk(f);

                upper_bound.post_order_walk(f);
            }
            Self::Sum {
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
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                scalar_field.post_order_walk(f);

                volume.post_order_walk(f);
            }
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                vector_field.post_order_walk(f);

                surface.post_order_walk(f);
            }
            Self::DerivativeN(e, _, n) => {

                e.post_order_walk(f);

                n.post_order_walk(f);
            }
            Self::Series(a, _, c, d) | Self::Summation(a, _, c, d) | Self::Product(a, _, c, d) => {

                a.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            }
            Self::AsymptoticExpansion(a, _, c, d) => {

                a.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            }
            Self::Interval(a, b, _, _) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            }
            Self::Substitute(a, _, c) => {

                a.post_order_walk(f);

                c.post_order_walk(f);
            }
            Self::Limit(a, _, c) => {

                a.post_order_walk(f);

                c.post_order_walk(f);
            }
            Self::Ode { equation, .. } => equation.post_order_walk(f),
            Self::Pde { equation, .. } => equation.post_order_walk(f),
            Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            }
            Self::ParametricSolution { x, y } => {

                x.post_order_walk(f);

                y.post_order_walk(f);
            }
            Self::QuantityWithValue(v, _) => v.post_order_walk(f),
            Self::RootOf { poly, .. } => poly.post_order_walk(f),

            Self::CustomArcOne(a) => {

                a.post_order_walk(f);
            }
            Self::CustomArcTwo(a, b) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            }
            Self::CustomArcThree(a, b, c) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);
            }
            Self::CustomArcFour(a, b, c, d) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);
            }
            Self::CustomArcFive(a, b, c, d, e) => {

                a.post_order_walk(f);

                b.post_order_walk(f);

                c.post_order_walk(f);

                d.post_order_walk(f);

                e.post_order_walk(f);
            }
            Self::CustomVecOne(v)
            | Self::CustomVecTwo(v, _)
            | Self::CustomVecThree(v, _, _)
            | Self::CustomVecFour(v, _, _, _)
            | Self::CustomVecFive(v, _, _, _, _) => {
                for e in v {

                    e.post_order_walk(f);
                }
            }
            Self::UnaryList(_, a) => a.post_order_walk(f),
            Self::BinaryList(_, a, b) => {

                a.post_order_walk(f);

                b.post_order_walk(f);
            }
            Self::NaryList(_, v) => {
                for e in v {

                    e.post_order_walk(f);
                }
            }
            Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.post_order_walk(f);
                }
            }
            // Leaf nodes
            Self::Constant(_)
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
            | Self::Distribution(_) => {}
        }

        f(self); // Visit parent
    }

    /// Performs an in-order traversal of the expression tree.
    /// For binary operators, it visits the left child, the node itself, then the right child.
    /// For other nodes, the behavior is adapted as it's not strictly defined.
    ///
    /// # Arguments
    /// * `f` - A mutable function that takes a reference to an `Expr` and is applied to each node during traversal
    ///

    pub fn in_order_walk<F>(&self, f: &mut F)
    where
        F: FnMut(&Self),
    {

        match self {
            // Binary operators
            Self::Add(a, b)
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
            }
            // Unary operators (treat as pre-order)
            Self::Sin(a)
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
            }
            // N-ary operators (visit self, then children)
            Self::Matrix(m) => {

                f(self);

                m.iter()
                    .flatten()
                    .for_each(|e| e.in_order_walk(f));
            }
            Self::Vector(v)
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
            }
            Self::Predicate { args, .. } => {

                f(self);

                for e in args {

                    e.in_order_walk(f);
                }
            }
            Self::SparsePolynomial(p) => {

                f(self);

                p.terms
                    .values()
                    .for_each(|c| c.in_order_walk(f));
            }
            // More complex operators (visit self, then children)
            Self::Integral {
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
            Self::Sum {
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
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                f(self);

                scalar_field.in_order_walk(f);

                volume.in_order_walk(f);
            }
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                f(self);

                vector_field.in_order_walk(f);

                surface.in_order_walk(f);
            }
            Self::DerivativeN(e, _, n) => {

                f(self);

                e.in_order_walk(f);

                n.in_order_walk(f);
            }
            Self::Series(a, _, c, d) | Self::Summation(a, _, c, d) | Self::Product(a, _, c, d) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);

                d.in_order_walk(f);
            }
            Self::AsymptoticExpansion(a, _, c, _d) => {

                f(self);

                a.in_order_walk(f);

                c.pre_order_walk(f);
            }
            Self::Interval(a, b, _, _) => {

                f(self);

                a.in_order_walk(f);

                b.in_order_walk(f);
            }
            Self::Substitute(a, _, c) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);
            }
            Self::Limit(a, _, c) => {

                f(self);

                a.in_order_walk(f);

                c.in_order_walk(f);
            }
            Self::Ode { equation, .. } => {

                f(self);

                equation.in_order_walk(f);
            }
            Self::Pde { equation, .. } => {

                f(self);

                equation.in_order_walk(f);
            }
            Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => {

                f(self);

                a.in_order_walk(f);

                b.in_order_walk(f);

                c.in_order_walk(f);

                d.pre_order_walk(f);
            }
            Self::ParametricSolution { x, y } => {

                f(self);

                x.in_order_walk(f);

                y.in_order_walk(f);
            }
            Self::RootOf { poly, .. } => {

                f(self);

                poly.in_order_walk(f);
            }
            Self::QuantityWithValue(v, _) => v.in_order_walk(f),

            Self::CustomArcOne(a) => {

                f(self);

                a.in_order_walk(f);
            }
            Self::CustomArcTwo(a, b) => {

                a.in_order_walk(f);

                f(self);

                b.in_order_walk(f);
            }
            Self::CustomArcThree(a, b, c) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);
            }
            Self::CustomArcFour(a, b, c, d) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);

                d.in_order_walk(f);
            }
            Self::CustomArcFive(a, b, c, d, e) => {

                a.in_order_walk(f);

                b.in_order_walk(f);

                f(self);

                c.in_order_walk(f);

                d.in_order_walk(f);

                e.in_order_walk(f);
            }
            Self::CustomVecOne(v) => {

                f(self);

                for e in v {

                    e.in_order_walk(f);
                }
            }
            Self::CustomVecTwo(v1, v2) => {

                f(self);

                for e in v1 {

                    e.in_order_walk(f);
                }

                for e in v2 {

                    e.in_order_walk(f);
                }
            }
            Self::CustomVecThree(v1, v2, v3) => {

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
            Self::CustomVecFour(v1, v2, v3, v4) => {

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
            Self::CustomVecFive(v1, v2, v3, v4, v5) => {

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
            Self::UnaryList(_, a) => {

                f(self);

                a.in_order_walk(f);
            }
            Self::BinaryList(_, a, b) => {

                a.in_order_walk(f);

                f(self);

                b.in_order_walk(f);
            }
            Self::NaryList(_, v) => {

                f(self);

                for e in v {

                    e.in_order_walk(f);
                }
            }

            // Leaf nodes
            Self::Constant(_)
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
            | Self::NoSolution => {}
            Self::Dag(node) => {

                // Convert DAG to AST and walk that to properly expose all nodes
                if let Ok(ast_expr) = node.to_expr() {

                    ast_expr.in_order_walk(f);
                }
            }
            Self::CustomZero | Self::CustomString(_) | Self::Distribution(_) => {}
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
    ///

    pub(crate) fn get_children_internal(&self) -> Vec<Self> {

        match self {
            Self::Add(a, b)
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
            | Self::Apply(a, b) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
            ],
            Self::AddList(v) | Self::MulList(v) => v.clone(),
            Self::Sin(a)
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
            Self::Matrix(m) => m
                .iter()
                .flatten()
                .cloned()
                .collect(),
            Self::Vector(v)
            | Self::Tuple(v)
            | Self::Polynomial(v)
            | Self::And(v)
            | Self::Or(v)
            | Self::Union(v)
            | Self::System(v)
            | Self::Solutions(v) => v.clone(),
            Self::Predicate { args, .. } => args.clone(),
            Self::SparsePolynomial(p) => p
                .terms
                .values()
                .cloned()
                .collect(),
            Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => vec![
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
            ],
            Self::VolumeIntegral {
                scalar_field,
                volume,
            } => vec![
                scalar_field
                    .as_ref()
                    .clone(),
                volume
                    .as_ref()
                    .clone(),
            ],
            Self::SurfaceIntegral {
                vector_field,
                surface,
            } => vec![
                vector_field
                    .as_ref()
                    .clone(),
                surface
                    .as_ref()
                    .clone(),
            ],
            Self::DerivativeN(e, _, n) => vec![
                e.as_ref().clone(),
                n.as_ref().clone(),
            ],
            Self::Series(a, _, c, d) | Self::Summation(a, _, c, d) | Self::Product(a, _, c, d) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            }
            Self::AsymptoticExpansion(a, _, c, d) => {

                vec![
                    a.as_ref().clone(),
                    c.as_ref().clone(),
                    d.as_ref().clone(),
                ]
            }
            Self::Interval(a, b, _, _) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
            ],
            Self::Substitute(a, _, c) => vec![
                a.as_ref().clone(),
                c.as_ref().clone(),
            ],
            Self::Limit(a, _, c) => vec![
                a.as_ref().clone(),
                c.as_ref().clone(),
            ],
            Self::Ode { equation, .. } => vec![equation
                .as_ref()
                .clone()],
            Self::Pde { equation, .. } => vec![equation
                .as_ref()
                .clone()],
            Self::Fredholm(a, b, c, d) | Self::Volterra(a, b, c, d) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
            ],
            Self::ParametricSolution { x, y } => vec![
                x.as_ref().clone(),
                y.as_ref().clone(),
            ],
            Self::RootOf { poly, .. } => vec![poly
                .as_ref()
                .clone()],
            Self::QuantityWithValue(v, _) => vec![v.as_ref().clone()],
            Self::CustomArcOne(a) => vec![a.as_ref().clone()],
            Self::CustomArcTwo(a, b) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
            ],
            Self::CustomArcThree(a, b, c) => {

                vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                    c.as_ref().clone(),
                ]
            }
            Self::CustomArcFour(a, b, c, d) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
            ],
            Self::CustomArcFive(a, b, c, d, e) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
                c.as_ref().clone(),
                d.as_ref().clone(),
                e.as_ref().clone(),
            ],
            Self::CustomVecOne(v) => v.clone(),
            Self::CustomVecTwo(v1, v2) => v1
                .iter()
                .chain(v2.iter())
                .cloned()
                .collect(),
            Self::CustomVecThree(v1, v2, v3) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .cloned()
                .collect(),
            Self::CustomVecFour(v1, v2, v3, v4) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .chain(v4.iter())
                .cloned()
                .collect(),
            Self::CustomVecFive(v1, v2, v3, v4, v5) => v1
                .iter()
                .chain(v2.iter())
                .chain(v3.iter())
                .chain(v4.iter())
                .chain(v5.iter())
                .cloned()
                .collect(),
            Self::UnaryList(_, a) => vec![a.as_ref().clone()],
            Self::BinaryList(_, a, b) => vec![
                a.as_ref().clone(),
                b.as_ref().clone(),
            ],
            Self::NaryList(_, v) => v.clone(),
            _ => vec![],
        }
    }

    #[must_use]

    pub fn normalize(&self) -> Self {

        match self {
            Self::Add(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Add(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Self::AddList(list) => {

                let mut children = Vec::new();

                for child in list {

                    if let Self::AddList(sub_list) = child {

                        children.extend(sub_list.clone());
                    } else {

                        children.push(child.clone());
                    }
                }

                children.sort();

                Self::AddList(children)
            }
            Self::Mul(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Mul(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Self::MulList(list) => {

                let mut children = Vec::new();

                for child in list {

                    if let Self::MulList(sub_list) = child {

                        children.extend(sub_list.clone());
                    } else {

                        children.push(child.clone());
                    }
                }

                children.sort();

                Self::MulList(children)
            }
            Self::Sub(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Sub(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Self::Div(a, b) => {

                let mut children = [
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ];

                children.sort();

                Self::Div(Arc::new(children[0].clone()), Arc::new(children[1].clone()))
            }
            Self::UnaryList(s, a) => Self::UnaryList(s.clone(), Arc::new(a.normalize())),
            Self::BinaryList(s, a, b) => {

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
                    Arc::new(children[0].clone()),
                    Arc::new(children[1].clone()),
                )
            }
            Self::NaryList(s, list) => {

                let mut children = list.clone();

                if let Some(props) = get_dynamic_op_properties(s) {

                    if props.is_commutative {

                        children.sort();
                    }
                }

                Self::NaryList(s.clone(), children)
            }
            _ => self.clone(),
        }
    }

    /// Converts this expression to its corresponding DAG operation.
    ///
    /// This method extracts the operation from an expression without its children,
    /// mapping it to the appropriate `DagOp` variant for use in the DAG system.
    ///
    /// # Returns
    /// * `Result<DagOp, String>` - The corresponding DAG operation or an error if conversion fails
    ///

    pub(crate) fn to_dag_op_internal(&self) -> Result<DagOp, String> {

        match self {
            Self::Constant(c) => Ok(DagOp::Constant(OrderedFloat(*c))),
            Self::BigInt(i) => Ok(DagOp::BigInt(i.clone())),
            Self::Rational(r) => Ok(DagOp::Rational(r.clone())),
            Self::Boolean(b) => Ok(DagOp::Boolean(*b)),
            Self::Variable(s) => Ok(DagOp::Variable(s.clone())),
            Self::Pattern(s) => Ok(DagOp::Pattern(s.clone())),
            Self::Domain(s) => Ok(DagOp::Domain(s.clone())),
            Self::Pi => Ok(DagOp::Pi),
            Self::E => Ok(DagOp::E),
            Self::Infinity => Ok(DagOp::Infinity),
            Self::NegativeInfinity => Ok(DagOp::NegativeInfinity),
            Self::InfiniteSolutions => Ok(DagOp::InfiniteSolutions),
            Self::NoSolution => Ok(DagOp::NoSolution),

            Self::Derivative(_, s) => Ok(DagOp::Derivative(s.clone())),
            Self::DerivativeN(_, s, _) => Ok(DagOp::DerivativeN(s.clone())),
            Self::Limit(_, s, _) => Ok(DagOp::Limit(s.clone())),
            Self::Solve(_, s) => Ok(DagOp::Solve(s.clone())),
            Self::ConvergenceAnalysis(_, s) => Ok(DagOp::ConvergenceAnalysis(s.clone())),
            Self::ForAll(s, _) => Ok(DagOp::ForAll(s.clone())),
            Self::Exists(s, _) => Ok(DagOp::Exists(s.clone())),
            Self::Substitute(_, s, _) => Ok(DagOp::Substitute(s.clone())),
            Self::Ode { func, var, .. } => Ok(DagOp::Ode {
                func: func.clone(),
                var: var.clone(),
            }),
            Self::Pde { func, vars, .. } => Ok(DagOp::Pde {
                func: func.clone(),
                vars: vars.clone(),
            }),
            Self::Predicate { name, .. } => Ok(DagOp::Predicate { name: name.clone() }),
            Self::Path(pt, _, _) => Ok(DagOp::Path(pt.clone())),
            Self::Interval(_, _, incl_lower, incl_upper) => {
                Ok(DagOp::Interval(*incl_lower, *incl_upper))
            }
            Self::RootOf { index, .. } => Ok(DagOp::RootOf { index: *index }),
            Self::SparsePolynomial(p) => Ok(DagOp::SparsePolynomial(p.clone())),
            Self::QuantityWithValue(_, u) => Ok(DagOp::QuantityWithValue(u.clone())),

            Self::Add(_, _) => Ok(DagOp::Add),
            Self::AddList(_) => Ok(DagOp::Add),
            Self::Sub(_, _) => Ok(DagOp::Sub),
            Self::Mul(_, _) => Ok(DagOp::Mul),
            Self::MulList(_) => Ok(DagOp::Mul),
            Self::Div(_, _) => Ok(DagOp::Div),
            Self::Neg(_) => Ok(DagOp::Neg),
            Self::Power(_, _) => Ok(DagOp::Power),
            Self::Sin(_) => Ok(DagOp::Sin),
            Self::Cos(_) => Ok(DagOp::Cos),
            Self::Tan(_) => Ok(DagOp::Tan),
            Self::Exp(_) => Ok(DagOp::Exp),
            Self::Log(_) => Ok(DagOp::Log),
            Self::Abs(_) => Ok(DagOp::Abs),
            Self::Sqrt(_) => Ok(DagOp::Sqrt),
            Self::Eq(_, _) => Ok(DagOp::Eq),
            Self::Lt(_, _) => Ok(DagOp::Lt),
            Self::Gt(_, _) => Ok(DagOp::Gt),
            Self::Le(_, _) => Ok(DagOp::Le),
            Self::Ge(_, _) => Ok(DagOp::Ge),
            Self::Matrix(m) => {

                let rows = m.len();

                let cols = if rows > 0 { m[0].len() } else { 0 };

                Ok(DagOp::Matrix { rows, cols })
            }
            Self::Vector(_) => Ok(DagOp::Vector),
            Self::Complex(_, _) => Ok(DagOp::Complex),
            Self::Transpose(_) => Ok(DagOp::Transpose),
            Self::MatrixMul(_, _) => Ok(DagOp::MatrixMul),
            Self::MatrixVecMul(_, _) => Ok(DagOp::MatrixVecMul),
            Self::Inverse(_) => Ok(DagOp::Inverse),
            Self::Integral { .. } => Ok(DagOp::Integral),
            Self::VolumeIntegral { .. } => Ok(DagOp::VolumeIntegral),
            Self::SurfaceIntegral { .. } => Ok(DagOp::SurfaceIntegral),
            Self::Sum { .. } => Ok(DagOp::Sum),
            Self::Series(_, s, _, _) => Ok(DagOp::Series(s.clone())),
            Self::Summation(_, s, _, _) => Ok(DagOp::Summation(s.clone())),
            Self::Product(_, s, _, _) => Ok(DagOp::Product(s.clone())),
            Self::AsymptoticExpansion(_, s, _, _) => Ok(DagOp::AsymptoticExpansion(s.clone())),
            Self::Sec(_) => Ok(DagOp::Sec),
            Self::Csc(_) => Ok(DagOp::Csc),
            Self::Cot(_) => Ok(DagOp::Cot),
            Self::ArcSin(_) => Ok(DagOp::ArcSin),
            Self::ArcCos(_) => Ok(DagOp::ArcCos),
            Self::ArcTan(_) => Ok(DagOp::ArcTan),
            Self::ArcSec(_) => Ok(DagOp::ArcSec),
            Self::ArcCsc(_) => Ok(DagOp::ArcCsc),
            Self::ArcCot(_) => Ok(DagOp::ArcCot),
            Self::Sinh(_) => Ok(DagOp::Sinh),
            Self::Cosh(_) => Ok(DagOp::Cosh),
            Self::Tanh(_) => Ok(DagOp::Tanh),
            Self::Sech(_) => Ok(DagOp::Sech),
            Self::Csch(_) => Ok(DagOp::Csch),
            Self::Coth(_) => Ok(DagOp::Coth),
            Self::ArcSinh(_) => Ok(DagOp::ArcSinh),
            Self::ArcCosh(_) => Ok(DagOp::ArcCosh),
            Self::ArcTanh(_) => Ok(DagOp::ArcTanh),
            Self::ArcSech(_) => Ok(DagOp::ArcSech),
            Self::ArcCsch(_) => Ok(DagOp::ArcCsch),
            Self::ArcCoth(_) => Ok(DagOp::ArcCoth),
            Self::LogBase(_, _) => Ok(DagOp::LogBase),
            Self::Atan2(_, _) => Ok(DagOp::Atan2),
            Self::Binomial(_, _) => Ok(DagOp::Binomial),
            Self::Factorial(_) => Ok(DagOp::Factorial),
            Self::Permutation(_, _) => Ok(DagOp::Permutation),
            Self::Combination(_, _) => Ok(DagOp::Combination),
            Self::FallingFactorial(_, _) => Ok(DagOp::FallingFactorial),
            Self::RisingFactorial(_, _) => Ok(DagOp::RisingFactorial),
            Self::Boundary(_) => Ok(DagOp::Boundary),
            Self::Gamma(_) => Ok(DagOp::Gamma),
            Self::Beta(_, _) => Ok(DagOp::Beta),
            Self::Erf(_) => Ok(DagOp::Erf),
            Self::Erfc(_) => Ok(DagOp::Erfc),
            Self::Erfi(_) => Ok(DagOp::Erfi),
            Self::Zeta(_) => Ok(DagOp::Zeta),
            Self::BesselJ(_, _) => Ok(DagOp::BesselJ),
            Self::BesselY(_, _) => Ok(DagOp::BesselY),
            Self::LegendreP(_, _) => Ok(DagOp::LegendreP),
            Self::LaguerreL(_, _) => Ok(DagOp::LaguerreL),
            Self::HermiteH(_, _) => Ok(DagOp::HermiteH),
            Self::Digamma(_) => Ok(DagOp::Digamma),
            Self::KroneckerDelta(_, _) => Ok(DagOp::KroneckerDelta),
            Self::And(_) => Ok(DagOp::And),
            Self::Or(_) => Ok(DagOp::Or),
            Self::Not(_) => Ok(DagOp::Not),
            Self::Xor(_, _) => Ok(DagOp::Xor),
            Self::Implies(_, _) => Ok(DagOp::Implies),
            Self::Equivalent(_, _) => Ok(DagOp::Equivalent),
            Self::Union(_) => Ok(DagOp::Union),
            Self::Polynomial(_) => Ok(DagOp::Polynomial),
            Self::Floor(_) => Ok(DagOp::Floor),
            Self::IsPrime(_) => Ok(DagOp::IsPrime),
            Self::Gcd(_, _) => Ok(DagOp::Gcd),
            Self::Mod(_, _) => Ok(DagOp::Mod),
            Self::System(_) => Ok(DagOp::System),
            Self::Solutions(_) => Ok(DagOp::Solutions),
            Self::ParametricSolution { .. } => Ok(DagOp::ParametricSolution),
            Self::GeneralSolution(_) => Ok(DagOp::GeneralSolution),
            Self::ParticularSolution(_) => Ok(DagOp::ParticularSolution),
            Self::Fredholm(_, _, _, _) => Ok(DagOp::Fredholm),
            Self::Volterra(_, _, _, _) => Ok(DagOp::Volterra),
            Self::Apply(_, _) => Ok(DagOp::Apply),
            Self::Tuple(_) => Ok(DagOp::Tuple),
            Self::Distribution(_) => Ok(DagOp::Distribution),
            Self::Max(_, _) => Ok(DagOp::Max),
            Self::Quantity(_) => Ok(DagOp::Quantity),
            Self::Dag(_) => Err("Cannot convert Dag to DagOp".to_string()),

            Self::CustomZero => Ok(DagOp::CustomZero),
            Self::CustomString(s) => Ok(DagOp::CustomString(s.clone())),
            Self::CustomArcOne(_) => Ok(DagOp::CustomArcOne),
            Self::CustomArcTwo(_, _) => Ok(DagOp::CustomArcTwo),
            Self::CustomArcThree(_, _, _) => Ok(DagOp::CustomArcThree),
            Self::CustomArcFour(_, _, _, _) => Ok(DagOp::CustomArcFour),
            Self::CustomArcFive(_, _, _, _, _) => Ok(DagOp::CustomArcFive),
            Self::CustomVecOne(_) => Ok(DagOp::CustomVecOne),
            Self::CustomVecTwo(_, _) => Ok(DagOp::CustomVecTwo),
            Self::CustomVecThree(_, _, _) => Ok(DagOp::CustomVecThree),
            Self::CustomVecFour(_, _, _, _) => Ok(DagOp::CustomVecFour),
            Self::CustomVecFive(_, _, _, _, _) => Ok(DagOp::CustomVecFive),
            Self::UnaryList(s, _) => Ok(DagOp::UnaryList(s.clone())),
            Self::BinaryList(s, _, _) => Ok(DagOp::BinaryList(s.clone())),
            Self::NaryList(s, _) => Ok(DagOp::NaryList(s.clone())),
        }
    }
}

impl AsRef<Self> for Expr {
    fn as_ref(&self) -> &Self {

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

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a])
                .expect("DAG manager get_or_create_normalized failed");

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

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let dag_b = DAG_MANAGER
                .get_or_create(b.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a, dag_b])
                .expect("DAG manager get_or_create_normalized failed");

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
                .map(|child| {

                    DAG_MANAGER
                        .get_or_create(child.as_ref())
                        .expect("DAG manager get_or_create failed")
                })
                .collect::<Vec<_>>();

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, children_nodes)
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the 'UnaryList' variant instead."
)]

macro_rules! unary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'UnaryList' variant instead."
        )]

        pub fn $name<A>(a: A) -> Expr
        where
            A: AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a])
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the 'BinaryList' variant instead."
)]

macro_rules! binary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'BinaryList' variant instead."
        )]

        pub fn $name<A, B>(a: A, b: B) -> Expr
        where
            A: AsRef<Expr>,
            B: AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let dag_b = DAG_MANAGER
                .get_or_create(b.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, vec![dag_a, dag_b])
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the 'NaryList' variant instead."
)]

macro_rules! n_ary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        #[doc = "Creates a new "]
        #[doc = stringify!($op)]
        #[doc = " expression, managed by the DAG."]
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'NaryList' variant instead."
        )]

        pub fn $name<I, T>(elements: I) -> Expr
        where
            I: IntoIterator<Item = T>,
            T: AsRef<Expr>,
        {

            let children_nodes = elements
                .into_iter()
                .map(|child| {

                    DAG_MANAGER
                        .get_or_create(child.as_ref())
                        .expect("DAG manager get_or_create failed")
                })
                .collect::<Vec<_>>();

            let node = DAG_MANAGER
                .get_or_create_normalized(DagOp::$op, children_nodes)
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

// --- Smart Constructors ---
impl Expr {
    // --- Leaf Node Constructors ---
    /// Creates a new Constant expression, managed by the DAG.
    #[must_use]

    pub fn new_constant(c: f64) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Constant(OrderedFloat(c)), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Variable expression, managed by the DAG.
    #[must_use]

    pub fn new_variable(name: &str) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Variable(name.to_string()), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new `BigInt` expression, managed by the DAG.
    #[must_use]

    pub fn new_bigint(i: BigInt) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::BigInt(i), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Rational expression, managed by the DAG.
    #[must_use]

    pub fn new_rational(r: BigRational) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Rational(r), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Pi expression, managed by the DAG.
    #[must_use]

    pub fn new_pi() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Pi, vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new E expression, managed by the DAG.
    #[must_use]

    pub fn new_e() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::E, vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Infinity expression, managed by the DAG.
    #[must_use]

    pub fn new_infinity() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Infinity, vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new `NegativeInfinity` expression, managed by the DAG.
    #[must_use]

    pub fn new_negative_infinity() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::NegativeInfinity, vec![])
            .expect("Value is valid");

        Self::Dag(node)
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

    pub fn new_matrix<I, J, T>(elements: I) -> Self
    where
        I: IntoIterator<Item = J>,
        J: IntoIterator<Item = T>,
        T: AsRef<Self>,
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

        Self::Dag(node)
    }

    pub fn new_predicate<I, T>(name: &str, args: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<Self>,
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

        Self::Dag(node)
    }

    pub fn new_forall<A>(var: &str, expr: A) -> Self
    where
        A: AsRef<Self>,
    {

        let child_node = DAG_MANAGER
            .get_or_create(expr.as_ref())
            .expect("Value is valid");

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::ForAll(var.to_string()), vec![child_node])
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_exists<A>(var: &str, expr: A) -> Self
    where
        A: AsRef<Self>,
    {

        let child_node = DAG_MANAGER
            .get_or_create(expr.as_ref())
            .expect("Value is valid");

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::Exists(var.to_string()), vec![child_node])
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_interval<A, B>(lower: A, upper: B, incl_lower: bool, incl_upper: bool) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
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

        Self::Dag(node)
    }

    #[must_use]

    pub fn new_sparse_polynomial(p: SparsePolynomial) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::SparsePolynomial(p), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    // --- Custom Constructors ---
    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    #[must_use]

    pub fn new_custom_zero() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomZero, vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'UnaryList' variant instead."
    )]
    #[must_use]

    pub fn new_custom_string(s: &str) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomString(s.to_string()), vec![])
            .expect("Value is valid");

        Self::Dag(node)
    }

    unary_constructor_deprecated!(new_custom_arc_one, CustomArcOne);

    binary_constructor_deprecated!(new_custom_arc_two, CustomArcTwo);

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]

    pub fn new_custom_arc_three<A, B, C>(a: A, b: B, c: C) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
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

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]

    pub fn new_custom_arc_four<A, B, C, D>(a: A, b: B, c: C, d: D) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
        D: AsRef<Self>,
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

        let children = vec![
            dag_a, dag_b, dag_c, dag_d,
        ];

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomArcFour, children)
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the 'NaryList' variant instead."
    )]

    pub fn new_custom_arc_five<A, B, C, D, E>(a: A, b: B, c: C, d: D, e: E) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
        D: AsRef<Self>,
        E: AsRef<Self>,
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

        let children = vec![
            dag_a, dag_b, dag_c, dag_d, dag_e,
        ];

        let node = DAG_MANAGER
            .get_or_create_normalized(DagOp::CustomArcFive, children)
            .expect("Value is valid");

        Self::Dag(node)
    }

    n_ary_constructor_deprecated!(new_custom_vec_one, CustomVecOne);

    n_ary_constructor_deprecated!(new_custom_vec_two, CustomVecTwo);

    n_ary_constructor_deprecated!(new_custom_vec_three, CustomVecThree);

    n_ary_constructor_deprecated!(new_custom_vec_four, CustomVecFour);

    n_ary_constructor_deprecated!(new_custom_vec_five, CustomVecFive);

    // --- AST to DAG Migration Utilities ---

    /// Checks if this expression is in DAG form.
    ///
    /// # Returns
    /// * `true` if the expression is `Expr::Dag`, `false` otherwise
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::Expr;
    ///
    /// let dag_expr = Expr::new_variable("x");
    /// assert!(dag_expr.is_dag());
    ///
    /// let ast_expr = Expr::Constant(1.0);
    /// assert!(!ast_expr.is_dag());
    /// ```
    #[inline]
    #[must_use]

    pub const fn is_dag(&self) -> bool {

        matches!(self, Self::Dag(_))
    }

    /// Converts this expression to DAG form if not already.
    ///
    /// This is a key function for the AST→DAG migration. It ensures that any
    /// expression, whether in old AST form or new DAG form, is converted to
    /// a canonical DAG representation.
    ///
    /// # Returns
    /// * `Ok(Expr::Dag)` - The expression in DAG form
    /// * `Err(String)` - If conversion fails
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::Expr;
    /// use std::sync::Arc;
    ///
    /// // Old AST form
    /// let ast = Expr::Add(
    ///     Arc::new(Expr::Variable("x".to_string())),
    ///     Arc::new(Expr::Constant(1.0))
    /// );
    ///
    /// // Convert to DAG
    /// let dag = ast.to_dag().unwrap();
    /// assert!(dag.is_dag());
    /// ```

    pub fn to_dag(&self) -> Result<Self, String> {

        match self {
            // Already in DAG form, just clone
            Self::Dag(_) => Ok(self.clone()),

            // Convert AST to DAG
            _ => {

                let dag_node = DAG_MANAGER.get_or_create(self)?;

                Ok(Self::Dag(dag_node))
            }
        }
    }

    /// Converts this expression to DAG form in-place.
    ///
    /// This is a convenience method that converts the expression to DAG form
    /// and replaces the current value.
    ///
    /// # Examples
    /// ```
    /// use rssn::symbolic::core::Expr;
    ///
    /// let mut expr = Expr::Constant(1.0);
    /// assert!(!expr.is_dag());
    ///
    /// expr.to_dag_form();
    /// assert!(expr.is_dag());
    /// ```

    pub fn to_dag_form(&mut self) {

        if let Ok(dag) = self.to_dag() {

            *self = dag;
        }
    }

    /// Converts a DAG expression back to AST form (for serialization).
    ///
    /// This is used internally for backward-compatible serialization.
    /// Note: This may fail for expressions that cannot be represented in AST form.
    ///
    /// # Returns
    /// * `Ok(Expr)` - The expression in AST form
    /// * `Err(String)` - If conversion fails

    pub fn to_ast(&self) -> Result<Self, String> {

        match self {
            Self::Dag(node) => node.to_expr(),
            _ => Ok(self.clone()),
        }
    }
}

// --- Dynamic Operation Registry ---
//
// This registry allows for runtime registration of custom operations without
// modifying the core `Expr` enum. It's designed to support:
// 1. Plugin systems that need to add new operations
// 2. Domain-specific extensions without core changes
// 3. Future-proofing the expression system
//
// Operations registered here can be used with UnaryList, BinaryList, and NaryList
// variants, and their properties (associativity, commutativity) are used during
// normalization and simplification.

/// Properties of a dynamically registered operation.
///
/// This struct defines the characteristics of a custom operation that can be
/// registered at runtime. These properties are used by the simplification engine
/// to apply appropriate transformations.
#[derive(Debug, Clone, Default)]

pub struct DynamicOpProperties {
    /// The name of the operation (must be unique).
    pub name: String,
    /// A human-readable description of what the operation does.
    pub description: String,
    /// Whether the operation is associative: f(f(a,b),c) = f(a,f(b,c))
    pub is_associative: bool,
    /// Whether the operation is commutative: f(a,b) = f(b,a)
    pub is_commutative: bool,
}

lazy_static! {
    /// Global registry for dynamically registered operations.
    ///
    /// This registry maps operation names to their properties. It's thread-safe
    /// and can be accessed from multiple threads simultaneously.
    pub static ref DYNAMIC_OP_REGISTRY: RwLock<HashMap<String, DynamicOpProperties>> = RwLock::new(HashMap::new());
}

/// Registers a dynamic operation with the global registry.
///
/// This function allows you to add custom operations at runtime without modifying
/// the core `Expr` enum. Once registered, operations can be used with `UnaryList`,
/// `BinaryList`, or `NaryList` variants.
///
/// # Arguments
/// * `name` - The unique name of the operation
/// * `props` - The properties of the operation (associativity, commutativity, etc.)
///
/// # Examples
/// ```
/// use rssn::symbolic::core::{register_dynamic_op, DynamicOpProperties};
///
/// register_dynamic_op("my_custom_op", DynamicOpProperties {
///     name: "my_custom_op".to_string(),
///     description: "A custom commutative operation".to_string(),
///     is_associative: true,
///     is_commutative: true,
/// });
/// ```

pub fn register_dynamic_op(name: &str, props: DynamicOpProperties) {

    let mut registry = DYNAMIC_OP_REGISTRY
        .write()
        .unwrap();

    registry.insert(name.to_string(), props);
}

/// Retrieves the properties of a dynamically registered operation.
///
/// This function looks up an operation in the global registry and returns its
/// properties if found.
///
/// # Arguments
/// * `name` - The name of the operation to look up
///
/// # Returns
/// * `Some(DynamicOpProperties)` if the operation is registered
/// * `None` if the operation is not found
///
/// # Examples
/// ```
/// use rssn::symbolic::core::{get_dynamic_op_properties, register_dynamic_op, DynamicOpProperties};
///
/// register_dynamic_op("test_op", DynamicOpProperties {
///     name: "test_op".to_string(),
///     description: "Test operation".to_string(),
///     is_associative: false,
///     is_commutative: true,
/// });
///
/// let props = get_dynamic_op_properties("test_op");
/// assert!(props.is_some());
/// assert!(props.unwrap().is_commutative);
/// ```
#[must_use]

pub fn get_dynamic_op_properties(name: &str) -> Option<DynamicOpProperties> {

    let registry = DYNAMIC_OP_REGISTRY
        .read()
        .unwrap();

    registry
        .get(name)
        .cloned()
}
