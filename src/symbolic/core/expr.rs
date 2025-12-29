#![allow(deprecated)]

use std::collections::BTreeMap;

use std::fmt::Debug;
use std::fmt::Write;
use std::fmt;
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;


use num_bigint::BigInt;
use num_rational::BigRational;


use super::dag_mgr::*;

use crate::symbolic::unit_unification::UnitQuantity;

// --- Distribution Trait ---
// Moved here to break circular dependency
/// Trait representing a probability distribution.

pub trait Distribution:
    Debug + Send + Sync
{
    /// Probability Density Function (PDF).

    fn pdf(
        &self,
        x: &Expr,
    ) -> Expr;

    /// Cumulative Distribution Function (CDF).

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr;

    /// Calculates the expected value of the distribution.

    fn expectation(&self) -> Expr;

    /// Calculates the variance of the distribution.

    fn variance(&self) -> Expr;

    /// Moment Generating Function (MGF).

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr;

    /// Creates a boxed clone of the distribution.

    fn clone_box(
        &self
    ) -> Arc<dyn Distribution>;
}

// --- End Distribution Trait ---

/// `PathType` enum
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]

pub enum PathType {
    /// A straight line path.
    Line,

    /// A circular path.
    Circle,

    /// A rectangular path.
    Rectangle,
}

/// Represents a single term in a multivariate polynomial.
///
/// A monomial is a product of variables raised to non-negative integer powers,
/// such as `x^2*y^3`. This struct stores it as a map from variable names (String)
/// to their exponents (u32).
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]

pub struct Monomial(
    pub BTreeMap<String, u32>,
);

/// Represents a sparse multivariate polynomial.
///
/// A sparse polynomial is stored as a map from `Monomial`s to their `Expr` coefficients.
/// This representation is highly efficient for polynomials with a small number of non-zero
/// terms relative to the degree, such as `x^1000 + 1`.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Hash,
    serde::Serialize,
    serde::Deserialize,
    PartialOrd,
    Ord,
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
#[derive(
    serde::Serialize, serde::Deserialize,
)]

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
    /// 
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
    /// 
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
    DerivativeN(
        Arc<Expr>,
        String,
        Arc<Expr>,
    ),
    /// A definite integral of `integrand` with respect to `var` from `lower_bound` to `upper_bound`.
    Integral {
        /// The expression to be integrated.
        integrand: Arc<Expr>,
        /// The variable of integration.
        var: Arc<Expr>,
        /// The lower limit of integration.
        lower_bound: Arc<Expr>,
        /// The upper limit of integration.
        upper_bound: Arc<Expr>,
    },
    /// A volume integral of a scalar field over a specified volume.
    VolumeIntegral {
        /// The scalar field to be integrated over the volume.
        scalar_field: Arc<Expr>,
        /// The volume domain.
        volume: Arc<Expr>,
    },
    /// A surface integral of a vector field over a specified surface.
    SurfaceIntegral {
        /// The vector field to be integrated over the surface.
        vector_field: Arc<Expr>,
        /// The surface domain.
        surface: Arc<Expr>,
    },
    /// A limit of an expression as a variable approaches a point.
    Limit(
        Arc<Expr>,
        String,
        Arc<Expr>,
    ),

    // --- Series and Summations ---
    /// A summation of `body` with `var` from `from` to `to`.
    Sum {
        /// The expression to be summed.
        body: Arc<Expr>,
        /// The summation variable.
        var: Arc<Expr>,
        /// The starting value of the summation variable.
        from: Arc<Expr>,
        /// The ending value of the summation variable.
        to: Arc<Expr>,
    },
    /// A finite or infinite series expansion.
    Series(
        Arc<Expr>,
        String,
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// A summation over a range (similar to `Sum`).
    Summation(
        Arc<Expr>,
        String,
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// A product of terms over a range.
    Product(
        Arc<Expr>,
        String,
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// Represents a convergence analysis for a series.
    ConvergenceAnalysis(
        Arc<Expr>,
        String,
    ),
    /// An asymptotic expansion of a function.
    AsymptoticExpansion(
        Arc<Expr>,
        String,
        Arc<Expr>,
        Arc<Expr>,
    ),

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
    FallingFactorial(
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// Rising factorial.
    RisingFactorial(
        Arc<Expr>,
        Arc<Expr>,
    ),

    // --- Geometry & Vector Calculus ---
    /// A path for path integrals (e.g., line, circle).
    Path(
        PathType,
        Arc<Expr>,
        Arc<Expr>,
    ),
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
    KroneckerDelta(
        Arc<Expr>,
        Arc<Expr>,
    ),

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
        /// The name of the predicate.
        name: String,
        /// The arguments to the predicate.
        args: Vec<Expr>,
    },
    /// Universal quantifier ("for all").
    ForAll(String, Arc<Expr>),
    /// Existential quantifier ("there exists").
    Exists(String, Arc<Expr>),
    /// A union of sets or intervals.
    Union(Vec<Expr>),
    /// An interval with a lower and upper bound, and flags for inclusion.
    Interval(
        Arc<Expr>,
        Arc<Expr>,
        bool,
        bool,
    ),

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
    Substitute(
        Arc<Expr>,
        String,
        Arc<Expr>,
    ),
    /// Represents a system of equations to be solved.
    System(Vec<Expr>),
    /// Represents the set of solutions to an equation or system.
    Solutions(Vec<Expr>),
    /// A parametric solution, e.g., for a system of ODEs.
    ParametricSolution {
        /// The x-component of the parametric solution.
        x: Arc<Expr>,
        /// The y-component of the parametric solution.
        y: Arc<Expr>,
    },
    /// Represents the `i`-th root of a polynomial.
    RootOf {
        /// The polynomial to find the root of.
        poly: Arc<Expr>,
        /// The index of the root (for polynomials with multiple roots).
        index: u32,
    },
    /// Represents infinite solutions.
    InfiniteSolutions,
    /// Represents that no solution exists.
    NoSolution,

    // --- Differential Equations ---
    /// An ordinary differential equation (ODE).
    Ode {
        /// The differential equation.
        equation: Arc<Expr>,
        /// The unknown function name.
        func: String,
        /// The independent variable name.
        var: String,
    },
    /// A partial differential equation (PDE).
    Pde {
        /// The partial differential equation.
        equation: Arc<Expr>,
        /// The unknown function name.
        func: String,
        /// The independent variable names.
        vars: Vec<String>,
    },
    /// The general solution to a differential equation.
    GeneralSolution(Arc<Expr>),
    /// A particular solution to a differential equation.
    ParticularSolution(Arc<Expr>),

    // --- Integral Equations ---
    /// A Fredholm integral equation.
    Fredholm(
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// A Volterra integral equation.
    Volterra(
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
    ),

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
    QuantityWithValue(
        Arc<Expr>,
        String,
    ),

    // --- Custom Variants (Old and Deprecated)---
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    /// Deprecated: Zero-argument custom operation.
    CustomZero,
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    /// Deprecated: String-argument custom operation.
    CustomString(String),

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    /// Deprecated: One-arc-argument custom operation.
    CustomArcOne(Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'BinaryList' variant \
                instead."
    )]
    /// Deprecated: Two-arc-argument custom operation.
    CustomArcTwo(Arc<Expr>, Arc<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Three-arc-argument custom operation.
    CustomArcThree(
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
    ),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Four-arc-argument custom operation.
    CustomArcFour(
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
    ),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Five-arc-argument custom operation.
    CustomArcFive(
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
        Arc<Expr>,
    ),

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    /// Deprecated: One-vector-argument custom operation.
    CustomVecOne(Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'BinaryList' variant \
                instead."
    )]
    /// Deprecated: Two-vector-argument custom operation.
    CustomVecTwo(Vec<Expr>, Vec<Expr>),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Three-vector-argument custom operation.
    CustomVecThree(
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
    ),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Four-vector-argument custom operation.
    CustomVecFour(
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
    ),
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]
    /// Deprecated: Five-vector-argument custom operation.
    CustomVecFive(
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
        Vec<Expr>,
    ),

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
    /// 
    /// use std::sync::Arc;
    ///
    /// use rssn::symbolic::core::register_dynamic_op;
    /// use rssn::symbolic::core::DynamicOpProperties;
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Register a custom operation
    /// register_dynamic_op(
    ///     "my_func",
    ///     DynamicOpProperties {
    ///         name : "my_func".to_string(),
    ///         description : "My custom function".to_string(),
    ///         is_associative : false,
    ///         is_commutative : false,
    ///     },
    /// );
    ///
    /// // Use it
    /// let expr = Expr::UnaryList(
    ///     "my_func".to_string(),
    ///     Arc::new(Expr::Variable(
    ///         "x".to_string(),
    ///     )),
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
    /// 
    /// use std::sync::Arc;
    ///
    /// use rssn::symbolic::core::register_dynamic_op;
    /// use rssn::symbolic::core::DynamicOpProperties;
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Register a custom binary operation
    /// register_dynamic_op(
    ///     "my_binop",
    ///     DynamicOpProperties {
    ///         name : "my_binop".to_string(),
    ///         description : "My custom binary operation"
    ///             .to_string(),
    ///         is_associative : true,
    ///         is_commutative : true,
    ///     },
    /// );
    ///
    /// // Use it
    /// let expr = Expr::BinaryList(
    ///     "my_binop".to_string(),
    ///     Arc::new(Expr::Variable(
    ///         "x".to_string(),
    ///     )),
    ///     Arc::new(Expr::Variable(
    ///         "y".to_string(),
    ///     )),
    /// );
    /// ```
    BinaryList(
        String,
        Arc<Expr>,
        Arc<Expr>,
    ),
    /// Generic n-ary operation identified by a name.
    ///
    /// This variant allows for extensible n-ary operations without modifying the `Expr` enum.
    /// Operations are registered in the `DYNAMIC_OP_REGISTRY` with their properties.
    /// This is particularly useful for operations that can take a variable number of arguments.
    ///
    /// # Examples
    /// ```
    /// 
    /// use rssn::symbolic::core::register_dynamic_op;
    /// use rssn::symbolic::core::DynamicOpProperties;
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Register a custom n-ary operation
    /// register_dynamic_op(
    ///     "my_nary",
    ///     DynamicOpProperties {
    ///         name : "my_nary".to_string(),
    ///         description : "My custom n-ary operation"
    ///             .to_string(),
    ///         is_associative : true,
    ///         is_commutative : false,
    ///     },
    /// );
    ///
    /// // Use it
    /// let expr = Expr::NaryList(
    ///     "my_nary".to_string(),
    ///     vec![
    ///         Expr::Variable("a".to_string()),
    ///         Expr::Variable("b".to_string()),
    ///         Expr::Variable("c".to_string()),
    ///     ],
    /// );
    /// ```
    NaryList(String, Vec<Expr>),
}
