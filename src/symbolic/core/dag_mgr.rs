#![allow(deprecated)]

use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::convert::AsRef;
use std::fmt::Debug;
use std::fmt::Write;
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::sync::LazyLock;
use dashmap::DashMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use ordered_float::OrderedFloat;

use super::expr::Expr;
use super::expr::PathType;
use super::expr::SparsePolynomial;

/// Global singleton instance of the DagManager.

pub static DAG_MANAGER: LazyLock<
    DagManager,
> = LazyLock::new(DagManager::new);


/// Represents a node in a Direct Acyclic Graph (DAG) for expression deduplication.
#[derive(
    Debug, Clone, serde::Serialize,
)]

pub struct DagNode {
    /// The operation performed at this node.
    pub op: DagOp,
    /// The children expressions of this node.
    pub children: Vec<Arc<Self>>,
    /// Precomputed 64-bit hash of the operation and children.
    #[serde(skip)]
    pub hash: u64,
    /// Unique ID for this node, assigned by the DagManager.
    #[serde(skip)]
    pub id: u64,
}

impl<'de> serde::Deserialize<'de>
    for DagNode
{
    fn deserialize<D>(
        deserializer: D
    ) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {

        #[derive(
            serde::Deserialize,
        )]

        struct DagNodeHelper {
            op: DagOp,
            children: Vec<Arc<DagNode>>,
        }

        let helper =
            DagNodeHelper::deserialize(
                deserializer,
            )?;

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

        // Assign a unique ID from the global manager to ensure consistency across the session.
        let id = DAG_MANAGER.next_id.fetch_add(1, AtomicOrdering::Relaxed);

        Ok(Self {
            op: helper.op,
            children: helper.children,
            hash,
            id,
        })
    }
}

#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    serde::Serialize,
    serde::Deserialize,
)]

/// Represents the different types of operations and leaf nodes in a Directed Acyclic Graph (DAG).

pub enum DagOp {
    // --- Leaf Nodes ---
    /// A double-precision floating-point constant.
    Constant(OrderedFloat<f64>),
    /// A large integer constant.
    BigInt(BigInt),
    /// A rational number constant.
    Rational(BigRational),
    /// A boolean constant.
    Boolean(bool),
    /// A variable identified by its name.
    Variable(String),
    /// A pattern variable for matching.
    Pattern(String),
    /// A domain definition.
    Domain(String),
    /// The mathematical constant pi (π).
    Pi,
    /// The mathematical constant e.
    E,
    /// Positive infinity.
    Infinity,
    /// Negative infinity.
    NegativeInfinity,
    /// Represents that a system has infinite solutions.
    InfiniteSolutions,
    /// Represents that a system has no solution.
    NoSolution,

    // --- Operators with associated data ---
    /// A derivative operation with respect to a variable.
    Derivative(String),
    /// An n-th order derivative operation.
    DerivativeN(String),
    /// A limit operation as a variable approaches a point.
    Limit(String),
    /// Represents solving an equation for a variable.
    Solve(String),
    /// Analysis of series convergence.
    ConvergenceAnalysis(String),
    /// Universal quantifier.
    ForAll(String),
    /// Existential quantifier.
    Exists(String),
    /// Substitution operation.
    Substitute(String),
    /// Ordinary Differential Equation.
    Ode {
        /// The name of the function being solved for.
        func: String,
        /// The independent variable.
        var: String,
    },
    /// Partial Differential Equation.
    Pde {
        /// The name of the function being solved for.
        func: String,
        /// The independent variables.
        vars: Vec<String>,
    },
    /// A logical or mathematical predicate.
    Predicate {
        /// The name of the predicate.
        name: String,
    },
    /// A geometric path.
    Path(PathType),
    /// An interval on the real line.
    Interval(bool, bool),

    // --- Operators without associated data (children only) ---
    /// Addition.
    Add,
    /// Subtraction.
    Sub,
    /// Multiplication.
    Mul,
    /// Division.
    Div,
    /// Negation.
    Neg,
    /// Exponentiation.
    Power,
    /// Sine function.
    Sin,
    /// Cosine function.
    Cos,
    /// Tangent function.
    Tan,
    /// Natural exponential function.
    Exp,
    /// Natural logarithm.
    Log,
    /// Absolute value.
    Abs,
    /// Square root.
    Sqrt,
    /// Equality comparison.
    Eq,
    /// Less than comparison.
    Lt,
    /// Greater than comparison.
    Gt,
    /// Less than or equal comparison.
    Le,
    /// Greater than or equal comparison.
    Ge,
    /// A matrix of a specific size.
    Matrix {
        /// Number of rows.
        rows: usize,
        /// Number of columns.
        cols: usize,
    },
    /// A vector.
    Vector,
    /// A complex number.
    Complex,
    /// Matrix transpose.
    Transpose,
    /// Matrix multiplication.
    MatrixMul,
    /// Matrix-vector multiplication.
    MatrixVecMul,
    /// Matrix inverse.
    Inverse,
    /// Definite integral.
    Integral,
    /// Volume integral.
    VolumeIntegral,
    /// Surface integral.
    SurfaceIntegral,
    /// Summation.
    Sum,
    /// Series expansion.
    Series(String),
    /// Summation over a range.
    Summation(String),
    /// Product over a range.
    Product(String),
    /// Asymptotic expansion.
    AsymptoticExpansion(String),
    /// Secant function.
    Sec,
    /// Cosecant function.
    Csc,
    /// Cotangent function.
    Cot,
    /// Inverse sine function.
    ArcSin,
    /// Inverse cosine function.
    ArcCos,
    /// Inverse tangent function.
    ArcTan,
    /// Inverse secant function.
    ArcSec,
    /// Inverse cosecant function.
    ArcCsc,
    /// Inverse cotangent function.
    ArcCot,
    /// Hyperbolic sine function.
    Sinh,
    /// Hyperbolic cosine function.
    Cosh,
    /// Hyperbolic tangent function.
    Tanh,
    /// Hyperbolic secant function.
    Sech,
    /// Hyperbolic cosecant function.
    Csch,
    /// Hyperbolic cotangent function.
    Coth,
    /// Inverse hyperbolic sine function.
    ArcSinh,
    /// Inverse hyperbolic cosine function.
    ArcCosh,
    /// Inverse hyperbolic tangent function.
    ArcTanh,
    /// Inverse hyperbolic secant function.
    ArcSech,
    /// Inverse hyperbolic cosecant function.
    ArcCsch,
    /// Inverse hyperbolic cotangent function.
    ArcCoth,
    /// Logarithm with a specific base.
    LogBase,
    /// 2-argument arctangent.
    Atan2,
    /// Binomial coefficient.
    Binomial,
    /// Factorial function.
    Factorial,
    /// Permutations.
    Permutation,
    /// Combinations.
    Combination,
    /// Falling factorial.
    FallingFactorial,
    /// Rising factorial.
    RisingFactorial,
    /// Domain boundary operator.
    Boundary,
    /// Gamma function.
    Gamma,
    /// Beta function.
    Beta,
    /// Error function.
    Erf,
    /// Complementary error function.
    Erfc,
    /// Imaginary error function.
    Erfi,
    /// Riemann zeta function.
    Zeta,
    /// Bessel function of the first kind.
    BesselJ,
    /// Bessel function of the second kind.
    BesselY,
    /// Legendre polynomial.
    LegendreP,
    /// Laguerre polynomial.
    LaguerreL,
    /// Hermite polynomial.
    HermiteH,
    /// Digamma function.
    Digamma,
    /// Kronecker delta function.
    KroneckerDelta,
    /// Logical AND.
    And,
    /// Logical OR.
    Or,
    /// Logical NOT.
    Not,
    /// Logical XOR.
    Xor,
    /// Logical implication.
    Implies,
    /// Logical equivalence.
    Equivalent,
    /// Set union.
    Union,
    /// Polynomial expression.
    Polynomial,
    /// Sparse polynomial.
    SparsePolynomial(SparsePolynomial), /* Note: Storing whole struct for simplicity */
    /// Floor function.
    Floor,
    /// Primality test.
    IsPrime,
    /// Greatest common divisor.
    Gcd,
    /// Modulo operation.
    Mod,
    /// System of equations.
    System,
    /// Solution set.
    Solutions,
    /// Parametric solution.
    ParametricSolution,
    /// The i-th root of a polynomial.
    RootOf {
        /// The index of the root.
        index: u32,
    },
    /// General solution to an ODE/PDE.
    GeneralSolution,
    /// Particular solution to an ODE/PDE.
    ParticularSolution,
    /// Fredholm integral equation.
    Fredholm,
    /// Volterra integral equation.
    Volterra,
    /// Function application.
    Apply,
    /// Tuple of expressions.
    Tuple,
    /// Probability distribution.
    Distribution, /* Trait objects are handled separately */
    /// Maximum of multiple values.
    Max,
    /// Physical quantity.
    Quantity, // Handled separately
    /// Physical quantity with a specific value.
    QuantityWithValue(String),

    // --- Custom ---
    /// Custom variant with zero children (deprecated).
    CustomZero,
    /// Custom variant with a string label (deprecated).
    CustomString(String),
    /// Custom variant with one child (deprecated).
    CustomArcOne,
    /// Custom variant with two children (deprecated).
    CustomArcTwo,
    /// Custom variant with three children (deprecated).
    CustomArcThree,
    /// Custom variant with four children (deprecated).
    CustomArcFour,
    /// Custom variant with five children (deprecated).
    CustomArcFive,
    /// Custom variant with one vector child (deprecated).
    CustomVecOne,
    /// Custom variant with two vector children (deprecated).
    CustomVecTwo,
    /// Custom variant with three vector children (deprecated).
    CustomVecThree,
    /// Custom variant with four vector children (deprecated).
    CustomVecFour,
    /// Custom variant with five vector children (deprecated).
    CustomVecFive,

    // --- Dynamic/Generic Operations ---
    /// List of children for unary operations.
    UnaryList(String),
    /// List of children for binary operations.
    BinaryList(String),
    /// List of children for n-ary operations.
    NaryList(String),
}

impl PartialEq for DagNode {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        // Use unique ID for O(1) comparison if available
        if self.id != 0 && other.id != 0 {
            return self.id == other.id;
        }

        // Fallback to structural equality for nodes without IDs (e.g., during deserialization)
        // This is now iterative to prevent stack overflow.
        if self.op != other.op || self.children.len() != other.children.len() {
            return false;
        }

        let mut stack = Vec::with_capacity(16);
        for (l, r) in self.children.iter().zip(other.children.iter()) {
            stack.push((l.as_ref(), r.as_ref()));
        }

        while let Some((l_node, r_node)) = stack.pop() {
            if std::ptr::eq(l_node, r_node) {
                continue;
            }
            if l_node.id != 0 && r_node.id != 0 {
                if l_node.id != r_node.id {
                    return false;
                }
                continue;
            }
            if l_node.op != r_node.op || l_node.children.len() != r_node.children.len() {
                return false;
            }
            for (lc, rc) in l_node.children.iter().zip(r_node.children.iter()) {
                stack.push((lc, rc));
            }
        }

        true
    }
}

impl Eq for DagNode {
}

impl PartialOrd for DagNode {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {

        Some(self.cmp(other))
    }
}

impl Ord for DagNode {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        // Use unique ID for comparison if available and we don't need structural ordering.
        // However, for canonicalization (Add/Mul sorting), we need a stable structural ordering.
        // IDs are not stable across runs/managers.
        
        // Iterative structural comparison
        let mut stack = Vec::with_capacity(16);
        stack.push((self, other));

        while let Some((l, r)) = stack.pop() {
            if std::ptr::eq(l, r) || (l.id != 0 && r.id != 0 && l.id == r.id) {
                continue;
            }

            match l.op.cmp(&r.op) {
                Ordering::Equal => {
                    match l.children.len().cmp(&r.children.len()) {
                        Ordering::Equal => {
                            for (lc, rc) in l.children.iter().zip(r.children.iter()).rev() {
                                stack.push((lc.as_ref(), rc.as_ref()));
                            }
                        }
                        ord => return ord,
                    }
                }
                ord => return ord,
            }
        }
        Ordering::Equal
    }
}

impl Hash for DagNode {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {
        if self.id != 0 {
            state.write_u64(self.id);
        } else {
            // Fallback for nodes without IDs
            self.op.hash(state);
            state.write_u64(self.hash);
        }
    }
}

impl From<DagNode> for Expr {
    fn from(node: DagNode) -> Self {

        node.to_expr()
            .expect(
                "Cannot convert \
                 DagNode to Expr.",
            )
    }
}

/// Manager for Directed Acyclic Graph (DAG) nodes.
///
/// This manager maintains a cache of `DagNode` instances to ensure that
/// identical expressions are represented by the same shared memory.

pub struct DagManager {
    nodes: DashMap<u64, Vec<Arc<DagNode>>>,
    next_id: AtomicU64,
}

impl Default for DagManager {
    fn default() -> Self {

        Self::new()
    }
}


impl DagManager {
    /// Creates a new `DagManager`.
    #[inline]
    #[must_use]

    pub fn new() -> Self {

        Self {
            nodes: DashMap::new(),
            next_id: AtomicU64::new(1),
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
    ///
    /// # Errors
    /// Returns an error if the number of children exceeds the safety limit.
    #[inline]

    pub fn get_or_create_normalized(
        &self,
        op: DagOp,
        mut children: Vec<Arc<DagNode>>,
    ) -> Result<Arc<DagNode>, String>
    {

        // Safety check: limit number of children to prevent excessive memory usage
        const MAX_CHILDREN: usize =
            10000;

        // Prevent excessive memory usage by limiting bucket size
        const MAX_BUCKET_SIZE: usize =
            1000;

        if children.len() > MAX_CHILDREN
        {

            return Err(format!(
                "Too many children in \
                 node ({}), exceeds \
                 limit of {}",
                children.len(),
                MAX_CHILDREN
            ));
        }

        match op {
            | DagOp::Add
            | DagOp::Mul => {

                // Use stable sort to ensure deterministic ordering across test runs.
                // This is critical for reproducible hashing and test stability.
                children.sort();
            },
            | _ => {},
        }

        // Compute 64-bit hash key
        let mut hasher =
            ahash::AHasher::default();

        op.hash(&mut hasher);

        for c in &children {

            // Use stored hash if present to avoid recursing
            Self::c_hash_for_hasher(
                c,
                &mut hasher,
            );
        }

        let hash = hasher.finish();

        // Use DashMap entry API for sharded, thread-safe access
        let mut entry = self.nodes.entry(hash).or_insert_with(Vec::new);
        let bucket = entry.value_mut();

        // Check bucket size to prevent excessive memory usage
        if bucket.len() > MAX_BUCKET_SIZE {
            let id = self.next_id.fetch_add(1, AtomicOrdering::Relaxed);
            let node = Arc::new(DagNode {
                op,
                children,
                hash,
                id,
            });
            return Ok(node);
        }

        // Build a temporary DagNode candidate for structural comparison.
        for cand in bucket.iter() {
            if Self::dag_nodes_structurally_equal(cand, &op, &children) {
                return Ok(cand.clone());
            }
        }

        // No structural match found in bucket: create new node and push.
        let id = self.next_id.fetch_add(1, AtomicOrdering::Relaxed);
        let node = Arc::new(DagNode {
            op,
            children,
            hash,
            id,
        });

        bucket.push(node.clone());
        Ok(node)
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

        if cand.children.len()
            != children.len()
        {

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

    pub(crate) fn compute_op_children_hash(
        op: &DagOp,
        children: &Vec<Arc<DagNode>>,
    ) -> u64 {

        let mut hasher =
            ahash::AHasher::default();

        op.hash(&mut hasher);

        for c in children {

            Self::c_hash_for_hasher(
                c,
                &mut hasher,
            );
        }

        hasher.finish()
    }

    /// Helper to feed a child's hash into hasher; uses stored hash field if available.

    pub(crate) fn c_hash_for_hasher(
        c: &Arc<DagNode>,
        hasher: &mut ahash::AHasher,
    ) {

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
    /// # Errors
    /// Returns an error if conversion from AST to DAG fails.
    #[inline]

    pub fn get_or_create(
        &self,
        expr: &Expr,
    ) -> Result<Arc<DagNode>, String>
    {
        if let Expr::Dag(node) = expr {
            return Ok(node.clone());
        }

        // Iterative post-order traversal using two stacks
        let mut stack = vec![(expr.clone(), false)];
        let mut result_stack: Vec<Arc<DagNode>> = Vec::new();

        while let Some((curr_expr, visited)) = stack.pop() {
            if let Expr::Dag(node) = &curr_expr {
                result_stack.push(node.clone());
                continue;
            }

            if visited {
                let op = curr_expr.to_dag_op_internal()?;
                let children_exprs = curr_expr.get_children_internal();
                let children_count = children_exprs.len();
                
                let mut children_nodes = Vec::with_capacity(children_count);
                for _ in 0..children_count {
                    children_nodes.push(result_stack.pop().ok_or("Result stack underflow")?);
                }
                children_nodes.reverse();
                
                let node = self.get_or_create_normalized(op, children_nodes)?;
                result_stack.push(node);
            } else {
                let children = curr_expr.get_children_internal();
                if children.is_empty() {
                    // Optimized path: process leaf nodes immediately without extra stack cycles
                    let op = curr_expr.to_dag_op_internal()?;
                    let node = self.get_or_create_normalized(op, Vec::new())?;
                    result_stack.push(node);
                } else {
                    stack.push((curr_expr, true));
                    for child in children.into_iter().rev() {
                        stack.push((child, false));
                    }
                }
            }
        }

        result_stack.pop().ok_or_else(|| "Root node not found".to_string())
    }
}
