#![allow(deprecated)]

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::convert::AsRef;
use std::fmt::Debug;
use std::fmt::Write;
use std::fmt::{
    self,
};
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

use lazy_static::lazy_static;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;

use super::api::*;
use super::ast_impl::*;
use super::expr::*;
use super::expr_impl::*;
use super::to_expr::*;
use crate::symbolic::unit_unification::UnitQuantity;

lazy_static! {
    pub static ref DAG_MANAGER: DagManager =
        DagManager::new();
}

#[derive(
    Debug, Clone, serde::Serialize,
)]

pub struct DagNode {
    pub op: DagOp,
    pub children: Vec<Arc<DagNode>>,
    #[serde(skip)]
    pub hash: u64,
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

        Ok(Self {
            op: helper.op,
            children: helper.children,
            hash,
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
    Ode {
        func: String,
        var: String,
    },
    Pde {
        func: String,
        vars: Vec<String>,
    },
    Predicate {
        name: String,
    },
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
    Matrix {
        rows: usize,
        cols: usize,
    },
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
    SparsePolynomial(SparsePolynomial), /* Note: Storing whole struct for simplicity */
    Floor,
    IsPrime,
    Gcd,
    Mod,
    System,
    Solutions,
    ParametricSolution,
    RootOf {
        index: u32,
    },
    GeneralSolution,
    ParticularSolution,
    Fredholm,
    Volterra,
    Apply,
    Tuple,
    Distribution, /* Trait objects are handled separately */
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
    fn eq(
        &self,
        other: &Self,
    ) -> bool {

        // 1. Check the operation
        if self.op != other.op {

            return false;
        }

        // 2. Check the length
        if self.children.len()
            != other.children.len()
        {

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
            .all(
                |(
                    l_child_arc,
                    r_child_arc,
                )| {

                    // This calls PartialEq recursively on the DagNode contents
                    l_child_arc
                        .as_ref()
                        .eq(r_child_arc
                            .as_ref())
                },
            )

        // NOTE: The hash field is IMPLICITLY IGNORED by this manual implementation.
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

        // A stable, structural comparison for canonical sorting.
        // Compare by operation type first, then recursively by children.
        self.op
            .cmp(&other.op)
            .then_with(|| {

                self.children.cmp(
                    &other.children,
                )
            })
    }
}

impl Hash for DagNode {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {

        self.op.hash(state);

        self.children
            .hash(state);
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

pub struct DagManager {
    nodes: Mutex<
        HashMap<u64, Vec<Arc<DagNode>>>,
    >,
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
            nodes: Mutex::new(
                HashMap::new(),
            ),
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
    ) -> Result<Arc<DagNode>, String>
    {

        // Safety check: limit number of children to prevent excessive memory usage
        const MAX_CHILDREN: usize =
            10000;

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

        // Acquire lock safely: handle PoisonError by recovering the inner guard.
        let mut nodes_guard =
            match self.nodes.lock() {
                | Ok(g) => g,
                | Err(pe) => {

                    // If a thread panicked previously, recover the poisoned lock's inner data.
                    // We prefer to continue with a best-effort recovery instead of panicking.
                    pe.into_inner()
                },
            };

        // Prevent excessive memory usage by limiting bucket size
        const MAX_BUCKET_SIZE: usize =
            1000;

        // Ensure the bucket is a vector of candidates to support collision buckets.
        // nodes: HashMap<u64, Vec<Arc<DagNode>>>
        match nodes_guard.entry(hash) {
            | Entry::Occupied(
                mut occ,
            ) => {

                // occ.get_mut() is a Vec<Arc<DagNode>>
                let bucket =
                    occ.get_mut();

                // Check bucket size to prevent excessive memory usage
                if bucket.len()
                    > MAX_BUCKET_SIZE
                {

                    // If the bucket is too large, we just create a new node without searching
                    // This maintains correctness while limiting memory usage
                    let node = Arc::new(
                        DagNode {
                            op,
                            children,
                            hash,
                        },
                    );

                    return Ok(node);
                }

                // Build a temporary DagNode candidate for structural comparison.
                // We avoid allocating the Arc until we know it's needed.
                for cand in
                    bucket.iter()
                {

                    if Self::dag_nodes_structurally_equal(cand, &op, &children) {

                        // Found exact structural match; return shared instance.
                        return Ok(cand.clone());
                    }
                }

                // No structural match found in bucket: create new node and push.
                let node =
                    Arc::new(DagNode {
                        op,
                        children,
                        hash,
                    });

                bucket
                    .push(node.clone());

                Ok(node)
            },
            | Entry::Vacant(vac) => {

                // No bucket yet: create a new vector with the node.
                let node =
                    Arc::new(DagNode {
                        op,
                        children,
                        hash,
                    });

                vac.insert(vec![
                    node.clone()
                ]);

                Ok(node)
            },
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
    #[inline]

    pub fn get_or_create(
        &self,
        expr: &Expr,
    ) -> Result<Arc<DagNode>, String>
    {

        if let Expr::Dag(node) = expr {

            return Ok(node.clone());
        }

        // Safety check: limit recursion depth to prevent stack overflow
        // We can't directly implement a recursion depth counter here since this is a single function,
        // but we can add checks for extremely complex expressions
        let op =
            expr.to_dag_op_internal()?;

        let children_exprs = expr
            .get_children_internal();

        // Limit the number of children to prevent excessive memory allocation
        const MAX_CHILDREN_PER_NODE:
            usize = 10000;

        if children_exprs.len()
            > MAX_CHILDREN_PER_NODE
        {

            return Err(format!(
                "Expression has too \
                 many children ({}), \
                 exceeds limit of {}",
                children_exprs.len(),
                MAX_CHILDREN_PER_NODE
            ));
        }

        let children_nodes = children_exprs
            .iter()
            .map(|child| self.get_or_create(child))
            .collect::<Result<Vec<_>, _>>()?;

        self.get_or_create_normalized(
            op,
            children_nodes,
        )
    }
}
