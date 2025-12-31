//! # Iterative DAG-based Symbolic Simplification
//!
//! This module provides a robust, iterative, and stack-safe symbolic simplification engine
//! based on a Directed Acyclic Graph (DAG) rewriting strategy. It is designed to work
//! from scratch without external libraries like `egg`.
//!
//! ## Core Strategy
//!
//! The simplification process is built on these key principles:
//!
//! 1.  **Iterative, Bottom-Up Traversal**: To prevent stack overflows on deeply nested
//!     expressions, the simplifier uses an explicit work stack to perform a post-order
//!     (bottom-up) traversal of the expression DAG. A node is only simplified after all
//!     of its children have been simplified.
//!
//! 2.  **Fixpoint Iteration**: Simplification rules are applied repeatedly across the entire
//!     expression until no further changes can be made. This ensures that simplifications
//!     are exhaustively applied and the expression reaches a stable, simplified form.
//!
//! 3.  **Canonicalization and Hashing**: All newly created expression nodes are managed
//!     through the central `DagManager`. This ensures that structurally identical
//!     expressions are represented by a single, canonical node in memory (content-addressing),
//!     maximizing sharing and efficiency.
//!
//! 4.  **Memoization**: During each pass of the fixpoint iteration, a memoization table
//!     (cache) is used to store the simplified result of each node. This avoids redundant
//!     work within a single pass.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use ordered_float::OrderedFloat;

#[allow(unused_imports)]
use super::core::DAG_MANAGER;
#[allow(unused_imports)]
use super::core::DagManager;
#[allow(unused_imports)]
use super::core::DagNode;
#[allow(unused_imports)]
use super::core::DagOp;
#[allow(unused_imports)]
use super::core::Expr;

/// The main entry point for the iterative DAG-based simplification.
///
/// This function orchestrates the entire simplification process. It sets up a fixpoint
/// loop that repeatedly applies a bottom-up simplification pass until the expression
/// converges to a stable form.
///
/// # Arguments
/// * `expr` - The expression to simplify.
///
/// # Returns
/// A new, simplified `Expr`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::simplify_dag::simplify;
///
/// let x = Expr::new_variable("x");
///
/// // x + x -> 2*x
/// let expr = Expr::new_add(&x, &x);
///
/// let simplified = simplify(&expr);
///
/// assert_eq!(
///     simplified,
///     Expr::new_mul(
///         Expr::new_constant(2.0),
///         &x
///     )
/// );
/// ```
#[must_use]

pub fn simplify(expr: &Expr) -> Expr {

    const MAX_ITERATIONS: usize = 1000; // Prevent infinite loops

    // Get the initial root node of the DAG from the input expression.
    let mut root_node =
        match DAG_MANAGER
            .get_or_create(expr)
        {
            | Ok(node) => node,
            | Err(_) => {

                return expr.clone();
            }, /* If creation fails, return the original expression. */
        };

    // --- Fixpoint Iteration Loop ---
    // This loop continues until a full pass over the DAG results in no changes.
    // We limit iterations to prevent infinite loops in case of bugs in simplification rules
    let mut iterations = 0;

    loop {

        let (simplified_root, changed) =
            bottom_up_simplify_pass(
                root_node.clone(),
            );

        if !changed {

            // If no changes were made in the last pass, the expression is fully simplified.
            break Expr::Dag(
                simplified_root,
            );
        }

        // Update the root node for the next iteration.
        root_node = simplified_root;

        iterations += 1;

        if iterations >= MAX_ITERATIONS
        {

            // If we've reached the maximum iterations, return the current simplified expression
            // to prevent infinite loops
            break Expr::Dag(root_node);
        }
    }
}

/// Performs a single, bottom-up simplification pass over the DAG.
///
/// This function uses an explicit stack (`work_stack`) to traverse the graph in a
/// post-order fashion. It ensures that a node's children are simplified before the
/// node itself is processed.
///
/// # Arguments
/// * `root` - The root `Arc<DagNode>` of the expression to simplify in this pass.
///
/// # Returns
/// A tuple containing:
/// * The new, simplified root node of the DAG.
/// * A boolean flag indicating whether any changes were made during the pass.

pub(crate) fn bottom_up_simplify_pass(
    root: Arc<DagNode>
) -> (Arc<DagNode>, bool) {

    const MAX_NODES_PER_PASS: usize =
        10000;

    // `memo` stores the simplified version of each node encountered in this pass.
    // Key: hash of the original node, Value: the simplified node.
    let mut memo: HashMap<
        u64,
        Arc<DagNode>,
    > = HashMap::new();

    // `work_stack` manages the nodes to be visited.
    let mut work_stack: Vec<
        Arc<DagNode>,
    > = vec![root.clone()];

    // `visited` keeps track of nodes pushed to the stack to avoid cycles and redundant work.
    let mut visited: HashMap<
        u64,
        bool,
    > = HashMap::new();

    let mut changed_in_pass = false;

    // Limit the number of nodes to prevent infinite loops in case of issues
    let mut processed_nodes = 0;

    while let Some(node) =
        work_stack.pop()
    {

        // If the node is already simplified in this pass, skip it.
        if memo.contains_key(&node.hash)
        {

            continue;
        }

        processed_nodes += 1;

        if processed_nodes
            >= MAX_NODES_PER_PASS
        {

            // If we've processed too many nodes, return early to prevent infinite processing
            break;
        }

        let children_simplified = node
            .children
            .iter()
            .all(|child| {

                memo.contains_key(
                    &child.hash,
                )
            });

        if children_simplified {

            // --- All children are simplified, so we can now process this node ---

            // 1. Rebuild the node with its (already simplified) children.
            let new_children: Vec<
                Arc<DagNode>,
            > =
                node.children
                    .iter()
                    .map(|child| {

                        // Use try_get to safely handle cases where child is not in memo
                        match memo.get(
                        &child.hash,
                    ) {
                        | Some(
                            child_node,
                        ) => child_node
                            .clone(),
                        | None => {
                            // This shouldn't happen if children_simplified is true,
                            // but we handle it for safety
                            child
                                .clone()
                        },
                    }
                    })
                    .collect();

            let rebuilt_node = match DAG_MANAGER.get_or_create_normalized(
                node.op.clone(),
                new_children,
            ) {
                | Ok(node) => node,
                | Err(_) => {

                    // If normalization fails, return the original node to avoid panics
                    continue;
                },
            };

            // 2. Apply simplification rules to the rebuilt node.
            let simplified_node =
                apply_rules(
                    &rebuilt_node,
                );

            // 3. Check if simplification changed the node.
            if simplified_node.hash
                != node.hash
            {

                changed_in_pass = true;
            }

            // 4. Memoize the result for the original node's hash.
            memo.insert(
                node.hash,
                simplified_node,
            );
        } else {

            // --- Not all children are simplified yet ---

            // 1. Push the current node back onto the stack to be processed later.
            work_stack
                .push(node.clone());

            // 2. Push un-simplified children onto the stack.
            if visited
                .insert(node.hash, true)
                .is_none()
            {

                for child in node
                    .children
                    .iter()
                    .rev()
                {

                    work_stack.push(
                        child.clone(),
                    );
                }
            }
        }
    }

    // The final simplified root is the one corresponding to the original root's hash.
    // Use get to safely handle cases where root wasn't processed
    let new_root = match memo
        .get(&root.hash)
    {
        | Some(node) => node.clone(),
        | None => root, /* If root wasn't processed, return the original */
    };

    (
        new_root,
        changed_in_pass,
    )
}

/// Applies a comprehensive set of algebraic and trigonometric simplification rules.
/// This function is the core of the simplification engine, performing pattern matching
/// on a single `DagNode` whose children are assumed to be already simplified.

pub(crate) fn apply_rules(
    node: &Arc<DagNode>
) -> Arc<DagNode> {

    // --- Constant Folding ---
    if let Some(folded) =
        fold_constants(node)
    {

        return folded;
    }

    match &node.op {
        // --- Arithmetic Rules ---
        | DagOp::Add => {

            if let Some(value) =
                apply_rules_add(node)
            {

                return value;
            }
        },
        | DagOp::Sub => {

            if let Some(value) =
                apply_rules_sub(node)
            {

                return value;
            }
        },
        | DagOp::Mul => {

            if let Some(value) =
                apply_rules_mul(node)
            {

                return value;
            }
        },
        | DagOp::Div => {

            if let Some(value) =
                apply_rules_div(node)
            {

                return value;
            }
        },
        | DagOp::Neg => {

            if let Some(value) =
                apply_rules_neg(node)
            {

                return value;
            }
        },

        // --- Power & Log/Exp Rules ---
        | DagOp::Power => {

            if let Some(value) =
                apply_rules_power(node)
            {

                return value;
            }
        },
        | DagOp::Log => {

            if let Some(value) =
                apply_rules_log(node)
            {

                return value;
            }
        },
        | DagOp::LogBase => {
            if let Some(value) =
                apply_rules_logbase(
                    node,
                )
            {

                return value;
            }
        },
        | DagOp::Exp => {

            if let Some(value) =
                apply_rules_exp(node)
            {

                return value;
            }
        },

        // --- Trigonometric Rules ---
        | DagOp::Sin => {

            if let Some(value) =
                apply_rules_sin(node)
            {

                return value;
            }
        },
        | DagOp::Cos => {

            if let Some(value) =
                apply_rules_cos(node)
            {

                return value;
            }
        },
        | DagOp::Tan => {

            if let Some(value) =
                apply_rules_tan(node)
            {

                return value;
            }
        },
        | DagOp::Sec => {

            if let Some(value) =
                apply_rules_sec(node)
            {

                return value;
            }
        },
        | DagOp::Csc => {

            if let Some(value) =
                apply_rules_csc(node)
            {

                return value;
            }
        },
        | DagOp::Cot => {

            if let Some(value) =
                apply_rules_cot(node)
            {

                return value;
            }
        },

        | _ => {}, /* No rule matched for this operator */
    }

    node.clone()
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_cot(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for cot operation
    }

    let arg = &node.children[0];

    // cot(x) -> cos(x)/sin(x)

    match DAG_MANAGER
        .get_or_create_normalized(
            DagOp::Cos,
            vec![arg.clone()],
        ) {
        | Ok(cos_x) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Sin,
                vec![arg.clone()],
            ) {
                | Ok(sin_x) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Div,
                        vec![cos_x, sin_x],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if division fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if sin(x) fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if cos(x) fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_csc(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for csc operation
    }

    let arg = &node.children[0];

    // csc(x) -> 1/sin(x)

    match DAG_MANAGER.get_or_create(
        &Expr::Constant(1.0),
    ) {
        | Ok(one) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Sin,
                vec![arg.clone()],
            ) {
                | Ok(sin_x) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Div,
                        vec![one, sin_x],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if division fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if sin(x) fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if constant creation fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_sec(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for sec operation
    }

    let arg = &node.children[0];

    // sec(x) -> 1/cos(x)

    match DAG_MANAGER.get_or_create(
        &Expr::Constant(1.0),
    ) {
        | Ok(one) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Cos,
                vec![arg.clone()],
            ) {
                | Ok(cos_x) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Div,
                        vec![one, cos_x],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if division fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if cos(x) fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if constant creation fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_tan(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for tan operation
    }

    let arg = &node.children[0];

    if is_zero_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // tan(0) -> 0

    // tan(x) -> sin(x)/cos(x)
    match DAG_MANAGER
        .get_or_create_normalized(
            DagOp::Sin,
            vec![arg.clone()],
        ) {
        | Ok(sin_x) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Cos,
                vec![arg.clone()],
            ) {
                | Ok(cos_x) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Div,
                        vec![sin_x, cos_x],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if division fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if cos(x) fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if sin(x) fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_cos(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for cos operation
    }

    let arg = &node.children[0];

    // cos(0) -> 1

    if is_zero_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // cos(pi) -> -1

    if is_pi_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(-1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // cos(-x) -> cos(x)

    if matches!(&arg.op, DagOp::Neg) {

        if arg
            .children
            .is_empty()
        {

            return Some(node.clone()); // Malformed negation, return original
        }

        return Some(match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Cos,
                vec![
                    arg.children[0]
                        .clone(),
                ],
            ) {
            | Ok(result) => result,
            | Err(_) => node.clone(), /* Return original if cos operation fails */
        });
    }

    // Sum/Difference and Induction formulas for cos
    if matches!(&arg.op, DagOp::Add) {

        if arg.children.len() >= 2 {

            let a = &arg.children[0];

            let b = &arg.children[1];

            // cos(x + pi) -> -cos(x)
            if is_pi_node(b) {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Cos,
                    vec![a.clone()],
                ) {
                    | Ok(cos_a) => {

                        match DAG_MANAGER.get_or_create_normalized(
                            DagOp::Neg,
                            vec![cos_a],
                        ) {
                            | Ok(result) => return Some(result),
                            | Err(_) => return Some(node.clone()), /* Return original if negation fails */
                        }
                    },
                    | Err(_) => return Some(node.clone()), /* Return original if cos operation fails */
                }
            }

            // cos(a+b) -> cos(a)cos(b) - sin(a)sin(b)
            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Cos,
                vec![a.clone()],
            ) {
                | Ok(cos_a) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Cos,
                        vec![b.clone()],
                    ) {
                        | Ok(cos_b) => {

                            match DAG_MANAGER.get_or_create_normalized(
                                DagOp::Sin,
                                vec![a.clone()],
                            ) {
                                | Ok(sin_a) => {

                                    match DAG_MANAGER.get_or_create_normalized(
                                        DagOp::Sin,
                                        vec![b.clone()],
                                    ) {
                                        | Ok(sin_b) => {

                                            match DAG_MANAGER.get_or_create_normalized(
                                                DagOp::Mul,
                                                vec![cos_a, cos_b],
                                            ) {
                                                | Ok(term1) => {

                                                    match DAG_MANAGER
                                                        .get_or_create_normalized(
                                                            DagOp::Mul,
                                                            vec![sin_a, sin_b],
                                                        ) {
                                                        | Ok(term2) => {

                                                            match DAG_MANAGER
                                                                .get_or_create_normalized(
                                                                    DagOp::Sub,
                                                                    vec![term1, term2],
                                                                ) {
                                                                | Ok(result) => {
                                                                    return Some(result)
                                                                },
                                                                | Err(_) => {
                                                                    return Some(node.clone())
                                                                }, /* Return original if subtraction fails */
                                                            }
                                                        },
                                                        | Err(_) => return Some(node.clone()), /* Return original if sin(b) fails */
                                                    }
                                                },
                                                | Err(_) => return Some(node.clone()), /* Return original if mul fails */
                                            }
                                        },
                                        | Err(_) => return Some(node.clone()), /* Return original if sin(b) fails */
                                    }
                                },
                                | Err(_) => return Some(node.clone()), /* Return original if sin(a) fails */
                            }
                        },
                        | Err(_) => return Some(node.clone()), // Return original if cos(b) fails
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if cos(a) fails
            }
        }

        return Some(node.clone()); // Malformed add, return original
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_sin(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for sin operation
    }

    let arg = &node.children[0];

    // sin(0) -> 0

    if is_zero_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // sin(pi) -> 0

    if is_pi_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // sin(-x) -> -sin(x)

    if matches!(&arg.op, DagOp::Neg) {

        if arg
            .children
            .is_empty()
        {

            return Some(node.clone()); // Malformed negation, return original
        }

        match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Sin,
                vec![
                    arg.children[0]
                        .clone(),
                ],
            ) {
            | Ok(new_sin) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Neg,
                    vec![new_sin],
                ) {
                    | Ok(result) => return Some(result),
                    | Err(_) => return Some(node.clone()), // Return original if negation fails
                }
            },
            | Err(_) => {

                return Some(
                    node.clone(),
                );
            }, /* Return original if sin operation fails */
        }
    }

    // Sum/Difference and Induction formulas for sin
    if matches!(&arg.op, DagOp::Add) {

        if arg.children.len() >= 2 {

            let a = &arg.children[0];

            let b = &arg.children[1];

            // sin(x + pi) -> -sin(x)
            if is_pi_node(b) {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Sin,
                    vec![a.clone()],
                ) {
                    | Ok(sin_a) => {

                        match DAG_MANAGER.get_or_create_normalized(
                            DagOp::Neg,
                            vec![sin_a],
                        ) {
                            | Ok(result) => return Some(result),
                            | Err(_) => return Some(node.clone()), /* Return original if negation fails */
                        }
                    },
                    | Err(_) => return Some(node.clone()), // Return original if sin(a) fails
                }
            }

            // sin(a+b) -> sin(a)cos(b) + cos(a)sin(b)
            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Sin,
                vec![a.clone()],
            ) {
                | Ok(sin_a) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Cos,
                        vec![b.clone()],
                    ) {
                        | Ok(cos_b) => {

                            match DAG_MANAGER.get_or_create_normalized(
                                DagOp::Cos,
                                vec![a.clone()],
                            ) {
                                | Ok(cos_a) => {

                                    match DAG_MANAGER.get_or_create_normalized(
                                        DagOp::Sin,
                                        vec![b.clone()],
                                    ) {
                                        | Ok(sin_b) => {

                                            match DAG_MANAGER.get_or_create_normalized(
                                                DagOp::Mul,
                                                vec![sin_a, cos_b],
                                            ) {
                                                | Ok(term1) => {

                                                    match DAG_MANAGER
                                                        .get_or_create_normalized(
                                                            DagOp::Mul,
                                                            vec![cos_a, sin_b],
                                                        ) {
                                                        | Ok(term2) => {

                                                            match DAG_MANAGER
                                                                .get_or_create_normalized(
                                                                    DagOp::Add,
                                                                    vec![term1, term2],
                                                                ) {
                                                                | Ok(result) => {
                                                                    return Some(result)
                                                                },
                                                                | Err(_) => {
                                                                    return Some(node.clone())
                                                                }, /* Return original if addition fails */
                                                            }
                                                        },
                                                        | Err(_) => return Some(node.clone()), /* Return original if mul fails */
                                                    }
                                                },
                                                | Err(_) => return Some(node.clone()), /* Return original if mul fails */
                                            }
                                        },
                                        | Err(_) => return Some(node.clone()), /* Return original if sin(b) fails */
                                    }
                                },
                                | Err(_) => return Some(node.clone()), /* Return original if cos(a) fails */
                            }
                        },
                        | Err(_) => return Some(node.clone()), // Return original if cos(b) fails
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if sin(a) fails
            }
        }

        return Some(node.clone()); // Malformed add, return original
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_exp(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for exp operation
    }

    let arg = &node.children[0];

    // exp(0) -> 1

    if is_zero_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // exp(log(x)) -> x
    if matches!(&arg.op, DagOp::Log) {

        if arg
            .children
            .is_empty()
        {

            return Some(node.clone()); // Malformed log, return original
        }

        return Some(
            arg.children[0].clone(),
        );
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_logbase(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for logbase operation
    }

    let base = &node.children[0];

    let arg = &node.children[1];

    // log_b(b) -> 1

    if base.hash == arg.hash {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // log_b(a) -> log(a) / log(b)
    match DAG_MANAGER
        .get_or_create_normalized(
            DagOp::Log,
            vec![arg.clone()],
        ) {
        | Ok(log_a) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Log,
                vec![base.clone()],
            ) {
                | Ok(log_b) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Div,
                        vec![log_a, log_b],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if division fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if log(base) fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if log(arg) fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_log(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for log operation
    }

    let arg = &node.children[0];

    // log(1) -> 0

    if is_one_node(arg) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // log(e) -> 1

    if matches!(&arg.op, DagOp::E) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // log(exp(x)) -> x

    if matches!(&arg.op, DagOp::Exp) {

        if arg
            .children
            .is_empty()
        {

            return Some(node.clone()); // Malformed exp, return original
        }

        return Some(
            arg.children[0].clone(),
        );
    }

    // log(a*b) -> log(a) + log(b)

    if matches!(&arg.op, DagOp::Mul) {

        if arg.children.len() >= 2 {

            let a = &arg.children[0];

            let b = &arg.children[1];

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Log,
                vec![a.clone()],
            ) {
                | Ok(log_a) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Log,
                        vec![b.clone()],
                    ) {
                        | Ok(log_b) => {

                            match DAG_MANAGER.get_or_create_normalized(
                                DagOp::Add,
                                vec![log_a, log_b],
                            ) {
                                | Ok(result) => return Some(result),
                                | Err(_) => return Some(node.clone()), /* Return original if addition fails */
                            }
                        },
                        | Err(_) => return Some(node.clone()), // Return original if log(b) fails
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if log(a) fails
            }
        }

        return Some(node.clone()); // Malformed mul, return original
    }

    // log(a^b) -> b*log(a)
    if matches!(
        &arg.op,
        DagOp::Power
    ) {

        if arg.children.len() >= 2 {

            let b =
                arg.children[1].clone();

            let log_a =
                &arg.children[0];

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Log,
                vec![log_a.clone()],
            ) {
                | Ok(log_a_node) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Mul,
                        vec![b, log_a_node],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if multiplication fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if log(a) fails
            }
        }

        return Some(node.clone()); // Malformed power, return original
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_power(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for power operation
    }

    let base = &node.children[0];

    let exp = &node.children[1];

    // x ^ 1 -> x

    if is_one_node(exp) {

        return Some(base.clone());
    }

    // x ^ 0 -> 1

    if is_zero_node(exp) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // 1 ^ x -> 1

    if is_one_node(base) {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // 0 ^ -x -> Infinity

    if is_zero_node(base) {

        match &exp.op {
            | DagOp::Constant(c) if c.0 < 0.0 => {

                return Some(match DAG_MANAGER.get_or_create(&Expr::Infinity) {
                    | Ok(node) => node,
                    | Err(_) => node.clone(),
                });
            },
            | DagOp::BigInt(b) if *b < BigInt::zero() => {

                return Some(match DAG_MANAGER.get_or_create(&Expr::Infinity) {
                    | Ok(node) => node,
                    | Err(_) => node.clone(),
                });
            },
            | DagOp::Rational(r) if *r < BigRational::zero() => {

                return Some(match DAG_MANAGER.get_or_create(&Expr::Infinity) {
                    | Ok(node) => node,
                    | Err(_) => node.clone(),
                });
            },
            | _ => {},
        }
    }

    // i^2 -> -1 (imaginary unit)

    if matches!(&base.op, DagOp::Variable(name) if name == "i")
        && (is_const_node(exp, 2.0)
            || matches!(&exp.op, DagOp::BigInt(b) if *b == BigInt::from(2)))
    {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(-1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(),
        });
    }

    // Sqrt(x) ^ 2 -> x

    if matches!(
        &base.op,
        DagOp::Sqrt
    ) && (is_const_node(exp, 2.0)
        || matches!(&exp.op, DagOp::BigInt(b) if *b == BigInt::from(2)))
        && !base
            .children
            .is_empty()
    {

        return Some(
            base.children[0].clone(),
        );
    }

    // (x^a)^b -> x^(a*b)
    if matches!(
        &base.op,
        DagOp::Power
    ) && base.children.len() >= 2
    {

        match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Mul,
                vec![
                    base.children[1]
                        .clone(),
                    exp.clone(),
                ],
            ) {
            | Ok(new_exp) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Power,
                    vec![
                        base.children[0].clone(),
                        new_exp,
                    ],
                ) {
                    | Ok(result) => return Some(result),
                    | Err(_) => return Some(node.clone()), /* Return original if power operation fails */
                }
            },
            | Err(_) => {}, /* Continue with other simplifications if this fails */
        }
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_neg(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node
        .children
        .is_empty()
    {

        return Some(node.clone()); // Not enough children for neg operation
    }

    let inner = &node.children[0];

    // --x -> x
    if matches!(
        &inner.op,
        DagOp::Neg
    ) {

        if inner
            .children
            .is_empty()
        {

            return Some(node.clone()); // Malformed negation, return original
        }

        return Some(
            inner.children[0].clone(),
        );
    }

    None
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_div(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for div operation
    }

    let lhs = &node.children[0];

    let rhs = &node.children[1];

    // x / 1 -> x

    if is_one_node(rhs) {

        return Some(lhs.clone());
    }

    // x / x -> 1 (if x != 0)

    if lhs.hash == rhs.hash {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(1.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // 0 / x -> 0 (if x != 0)

    if is_zero_node(lhs) {

        if is_zero_node(rhs) {

            return Some(node.clone());
        }

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(),
        });
    }

    match DAG_MANAGER.get_or_create(
        &Expr::BigInt(BigInt::from(-1)),
    ) {
        | Ok(neg_one) => {

            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Power,
                vec![rhs.clone(), neg_one],
            ) {
                | Ok(rhs_pow_neg_one) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Mul,
                        vec![
                            lhs.clone(),
                            rhs_pow_neg_one,
                        ],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => return Some(node.clone()), /* Return original if multiplication fails */
                    }
                },
                | Err(_) => return Some(node.clone()), // Return original if power operation fails
            }
        },
        | Err(_) => {

            return Some(node.clone());
        }, /* Return original if neg_one creation fails */
    }
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_mul(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for mul operation
    }

    let lhs = &node.children[0];

    let rhs = &node.children[1];

    // x * 1 -> x

    if is_one_node(rhs) {

        return Some(lhs.clone());
    }

    if is_one_node(lhs) {

        return Some(rhs.clone());
    }

    // x * 0 -> 0

    if is_zero_node(rhs)
        || is_zero_node(lhs)
    {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // x * -1 -> -x

    if is_neg_one_node(rhs) {

        return Some(match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Neg,
                vec![lhs.clone()],
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    if is_neg_one_node(lhs) {

        return Some(match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Neg,
                vec![rhs.clone()],
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // x * x -> x^2

    if lhs.hash == rhs.hash {

        match DAG_MANAGER.get_or_create(
            &Expr::Constant(2.0),
        ) {
            | Ok(two) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Power,
                    vec![lhs.clone(), two],
                ) {
                    | Ok(result) => return Some(result),
                    | Err(_) => {}, // Continue with other simplifications if this fails
                }
            },
            | Err(_) => {}, /* Continue with other simplifications if this fails */
        }
    }

    // Distributivity: a * (b + c) -> a*b + a*c

    if matches!(&rhs.op, DagOp::Add)
        && rhs.children.len() >= 2
    {

        let a = lhs;

        let b = &rhs.children[0];

        let c = &rhs.children[1];

        match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Mul,
                vec![
                    a.clone(),
                    b.clone(),
                ],
            ) {
            | Ok(ab) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Mul,
                    vec![a.clone(), c.clone()],
                ) {
                    | Ok(ac) => {

                        return Some(match DAG_MANAGER.get_or_create_normalized(
                            DagOp::Add,
                            vec![ab, ac],
                        ) {
                            | Ok(result) => result,
                            | Err(_) => node.clone(), // Return original if addition fails
                        });
                    },
                    | Err(_) => {}, // Continue with other simplifications if this fails
                }
            },
            | Err(_) => {}, /* Continue with other simplifications if this fails */
        }
    }

    if matches!(&lhs.op, DagOp::Add)
        && lhs.children.len() >= 2
    {

        let a = &lhs.children[0];

        let b = &lhs.children[1];

        let c = rhs;

        match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Mul,
                vec![
                    a.clone(),
                    c.clone(),
                ],
            ) {
            | Ok(ac) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Mul,
                    vec![b.clone(), c.clone()],
                ) {
                    | Ok(bc) => {

                        return Some(match DAG_MANAGER.get_or_create_normalized(
                            DagOp::Add,
                            vec![ac, bc],
                        ) {
                            | Ok(result) => result,
                            | Err(_) => node.clone(), // Return original if addition fails
                        });
                    },
                    | Err(_) => {}, // Continue with other simplifications if this fails
                }
            },
            | Err(_) => {}, /* Continue with other simplifications if this fails */
        }
    }

    Some(simplify_mul(node))
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_sub(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for sub operation
    }

    let lhs = &node.children[0];

    let rhs = &node.children[1];

    // x - 0 -> x

    if is_zero_node(rhs) {

        return Some(lhs.clone());
    }

    // a - (-b) -> a + b

    if matches!(&rhs.op, DagOp::Neg)
        && !rhs
            .children
            .is_empty()
    {

        let b = &rhs.children[0];

        return Some(match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Add,
                vec![
                    lhs.clone(),
                    b.clone(),
                ],
            ) {
            | Ok(result) => result,
            | Err(_) => node.clone(),
        });
    }

    // x - x -> 0

    if lhs.hash == rhs.hash {

        return Some(match DAG_MANAGER
            .get_or_create(
                &Expr::Constant(0.0),
            ) {
            | Ok(node) => node,
            | Err(_) => node.clone(), /* Return original if creation fails */
        });
    }

    // a * x - b * x -> (a-b)*x

    if matches!(&lhs.op, DagOp::Mul)
        && matches!(&rhs.op, DagOp::Mul)
    {

        let mut terms_lhs = Vec::new();

        flatten_terms(
            lhs,
            &mut terms_lhs,
        );

        let mut terms_rhs = Vec::new();

        flatten_terms(
            rhs,
            &mut terms_rhs,
        );

        let one_node_a = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Constant(
                    OrderedFloat(1.0),
                ), // Argument 1: The DagOp, correctly constructed
                vec![], // Argument 2: The children, empty for a constant
            )
            .unwrap_or_else(|_| {
                node.clone()
            });

        let (a, b) = if lhs
            .children
            .len()
            < 2
            || terms_lhs.len() < 2
        {

            (
                one_node_a,
                terms_lhs[0].clone(),
            )
        } else {

            (
                terms_lhs[0].clone(),
                terms_lhs[1].clone(),
            )
        };

        let one_node_c = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Constant(
                    OrderedFloat(1.0),
                ), // Argument 1: The DagOp
                vec![], // Argument 2: The children
            )
            .unwrap_or_else(|_| {
                node.clone()
            });

        let (c, d) = if rhs
            .children
            .len()
            < 2
            || terms_rhs.len() < 2
        {

            (
                one_node_c,
                terms_rhs[0].clone(),
            )
        } else {

            (
                terms_rhs[0].clone(),
                terms_rhs[1].clone(),
            )
        };

        if b.hash == d.hash {

            return Some(match DAG_MANAGER.get_or_create_normalized(
                DagOp::Sub,
                vec![a, c],
            ) {
                | Ok(result) => result,
                | Err(_) => node.clone(),
            });
        }
    }

    Some(simplify_add(node))
}

#[inline(always)]
#[allow(clippy::inline_always)]
#[allow(clippy::unnecessary_wraps)]

pub(crate) fn apply_rules_add(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    if node.children.len() < 2 {

        return Some(node.clone()); // Not enough children for add operation
    }

    let lhs = &node.children[0];

    let rhs = &node.children[1];

    // x + 0 -> x

    if is_zero_node(rhs) {

        return Some(lhs.clone());
    }

    if is_zero_node(lhs) {

        return Some(rhs.clone());
    }

    // x + x -> 2*x

    if lhs.hash == rhs.hash {

        match DAG_MANAGER.get_or_create(
            &Expr::Constant(2.0),
        ) {
            | Ok(two) => {

                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Mul,
                    vec![two, lhs.clone()],
                ) {
                    | Ok(result) => return Some(result),
                    | Err(_) => {}, // Continue with other simplifications if this fails
                }
            },
            | Err(_) => {}, /* Continue with other simplifications if this fails */
        }
    }

    // Coefficient Collection: ax + bx -> (a+b)x

    if matches!(
        (&lhs.op, &rhs.op),
        (
            DagOp::Mul,
            DagOp::Mul
        )
    ) && lhs.children.len() >= 2
        && rhs.children.len() >= 2
    {

        let a = &lhs.children[0];

        let x1 = &lhs.children[1];

        let b = &rhs.children[0];

        let x2 = &rhs.children[1];

        if x1.hash == x2.hash {

            // a*x + b*x
            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Add,
                vec![a.clone(), b.clone()],
            ) {
                | Ok(a_plus_b) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Mul,
                        vec![a_plus_b, x1.clone()],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => {}, // Continue with other simplifications if this fails
                    }
                },
                | Err(_) => {}, // Continue with other simplifications if this fails
            }
        }

        if a.hash == b.hash {

            // x*a + y*a -> (x+y)*a
            match DAG_MANAGER.get_or_create_normalized(
                DagOp::Add,
                vec![
                    x1.clone(),
                    x2.clone(),
                ],
            ) {
                | Ok(x_plus_y) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Mul,
                        vec![x_plus_y, a.clone()],
                    ) {
                        | Ok(result) => return Some(result),
                        | Err(_) => {}, // Continue with other simplifications if this fails
                    }
                },
                | Err(_) => {}, // Continue with other simplifications if this fails
            }
        }
    }

    // sin(x)^2 + cos(x)^2 -> 1

    if matches!(
        (&lhs.op, &rhs.op),
        (
            DagOp::Power,
            DagOp::Power
        )
    ) && lhs.children.len() >= 2
        && rhs.children.len() >= 2
        && is_const_node(
            &lhs.children[1],
            2.0,
        )
        && is_const_node(
            &rhs.children[1],
            2.0,
        )
    {

        // Both are squared
        if matches!(
            (
                &lhs.children[0].op,
                &rhs.children[0].op
            ),
            (
                DagOp::Sin,
                DagOp::Cos
            )
        ) && !lhs.children[0]
            .children
            .is_empty()
            && !rhs.children[0]
                .children
                .is_empty()
            && lhs.children[0].children
                [0]
            .hash
                == rhs.children[0]
                    .children[0]
                    .hash
        {

            // sin(arg) and cos(arg) with same arg
            return Some(
                match DAG_MANAGER
                    .get_or_create(
                        &Expr::Constant(
                            1.0,
                        ),
                    ) {
                    | Ok(node) => node,
                    | Err(_) => {
                        node.clone()
                    }, /* Return original if creation fails */
                },
            );
        }

        if matches!(
            (
                &lhs.children[0].op,
                &rhs.children[0].op
            ),
            (
                DagOp::Cos,
                DagOp::Sin
            )
        ) && !lhs.children[0]
            .children
            .is_empty()
            && !rhs.children[0]
                .children
                .is_empty()
            && lhs.children[0].children
                [0]
            .hash
                == rhs.children[0]
                    .children[0]
                    .hash
        {

            // cos(arg) and sin(arg) with same arg
            return Some(
                match DAG_MANAGER
                    .get_or_create(
                        &Expr::Constant(
                            1.0,
                        ),
                    ) {
                    | Ok(node) => node,
                    | Err(_) => {
                        node.clone()
                    }, /* Return original if creation fails */
                },
            );
        }
    }

    Some(simplify_add(node))
}

/// Performs constant folding on a `DagNode`.
/// If the node is an operation on constant children, it computes the result and returns a new constant node.

pub(crate) fn fold_constants(
    node: &Arc<DagNode>
) -> Option<Arc<DagNode>> {

    let children_values: Option<
        Vec<Expr>,
    > = node
        .children
        .iter()
        .map(get_numeric_value)
        .collect();

    if let Some(values) =
        children_values
    {

        let result = match (
            &node.op,
            values.as_slice(),
        ) {
            | (DagOp::Add, [a, b]) => Some(add_em(a, b)),
            | (DagOp::Sub, [a, b]) => Some(sub_em(a, b)),
            | (DagOp::Mul, [a, b]) => Some(mul_em(a, b)),
            | (DagOp::Div, [a, b]) => div_em(a, b),
            | (DagOp::Power, [Expr::Constant(a), Expr::Constant(b)]) => {
                Some(Expr::Constant(
                    a.powf(*b),
                ))
            },
            | (DagOp::Neg, [a]) => Some(neg_em(a)),
            | (DagOp::Sqrt, [a]) => {
                match a.to_f64() {
                    | Some(val) if val >= 0.0 => {
                        Some(Expr::Constant(
                            val.sqrt(),
                        ))
                    },
                    | _ => None,
                }
            },
            | _ => None,
        };

        if let Some(value) = result {

            return match DAG_MANAGER
                .get_or_create(&value)
            {
                | Ok(node) => {
                    Some(node)
                },
                | Err(_) => None, /* If creation fails, return None to continue with original */
            };
        }
    }

    None
}

// --- Numeric Helper Functions ---

#[inline]

pub(crate) fn get_numeric_value(
    node: &Arc<DagNode>
) -> Option<Expr> {

    match &node.op {
        | DagOp::Constant(c) => {
            Some(Expr::Constant(
                c.into_inner(),
            ))
        },
        | DagOp::BigInt(i) => {
            Some(Expr::BigInt(
                i.clone(),
            ))
        },
        | DagOp::Rational(r) => {
            Some(Expr::Rational(
                r.clone(),
            ))
        },
        | _ => None,
    }
}

#[inline]

pub(crate) fn add_em(
    a: &Expr,
    b: &Expr,
) -> Expr {

    match (a, b) {
        | (
            Expr::Constant(va),
            Expr::Constant(vb),
        ) => {

            let result = va + vb;

            if result.is_infinite()
                || result.is_nan()
            {

                // Handle overflow/invalid results gracefully
                Expr::Constant(*va) // Return a as a fallback
            } else {

                Expr::Constant(result)
            }
        },
        | (
            Expr::BigInt(ia),
            Expr::BigInt(ib),
        ) => Expr::BigInt(ia + ib),
        | (
            Expr::Rational(ra),
            Expr::Rational(rb),
        ) => Expr::Rational(ra + rb),
        // Promote to Rational or Constant - with error handling
        | _ => {

            match (
                a.to_f64(),
                b.to_f64(),
            ) {
                | (
                    Some(va),
                    Some(vb),
                ) => {

                    let result =
                        va + vb;

                    if result
                        .is_infinite()
                        || result
                            .is_nan()
                    {

                        Expr::new_add(
                            a, b,
                        ) // Return original expression if result is invalid
                    } else {

                        Expr::Constant(
                            result,
                        )
                    }
                },
                | _ => {
                    Expr::new_add(a, b)
                }, /* Return original expression if conversion fails */
            }
        },
    }
}

#[inline]

pub(crate) fn sub_em(
    a: &Expr,
    b: &Expr,
) -> Expr {

    match (a, b) {
        | (
            Expr::Constant(va),
            Expr::Constant(vb),
        ) => {

            let result = va - vb;

            if result.is_infinite()
                || result.is_nan()
            {

                // Handle overflow/invalid results gracefully
                Expr::Constant(*va) // Return a as a fallback
            } else {

                Expr::Constant(result)
            }
        },
        | (
            Expr::BigInt(ia),
            Expr::BigInt(ib),
        ) => Expr::BigInt(ia - ib),
        | (
            Expr::Rational(ra),
            Expr::Rational(rb),
        ) => Expr::Rational(ra - rb),
        | _ => {

            match (
                a.to_f64(),
                b.to_f64(),
            ) {
                | (
                    Some(va),
                    Some(vb),
                ) => {

                    let result =
                        va - vb;

                    if result
                        .is_infinite()
                        || result
                            .is_nan()
                    {

                        Expr::new_sub(
                            a, b,
                        ) // Return original expression if result is invalid
                    } else {

                        Expr::Constant(
                            result,
                        )
                    }
                },
                | _ => {
                    Expr::new_sub(a, b)
                }, /* Return original expression if conversion fails */
            }
        },
    }
}

#[inline]

pub(crate) fn mul_em(
    a: &Expr,
    b: &Expr,
) -> Expr {

    match (a, b) {
        | (
            Expr::Constant(va),
            Expr::Constant(vb),
        ) => {

            let result = va * vb;

            if result.is_infinite()
                || result.is_nan()
            {

                // Handle overflow/invalid results gracefully
                if (va.is_infinite()
                    && is_zero_expr(b))
                    || (vb
                        .is_infinite()
                        && is_zero_expr(
                            a,
                        ))
                {

                    Expr::Constant(0.0) // 0 * inf = 0 (though mathematically indeterminate, for calculation purposes)
                } else {

                    Expr::Constant(*va) // Return a as a fallback
                }
            } else {

                Expr::Constant(result)
            }
        },
        | (
            Expr::BigInt(ia),
            Expr::BigInt(ib),
        ) => Expr::BigInt(ia * ib),
        | (
            Expr::Rational(ra),
            Expr::Rational(rb),
        ) => Expr::Rational(ra * rb),
        | _ => {

            match (
                a.to_f64(),
                b.to_f64(),
            ) {
                | (
                    Some(va),
                    Some(vb),
                ) => {

                    let result =
                        va * vb;

                    if result
                        .is_infinite()
                        || result
                            .is_nan()
                    {

                        Expr::new_mul(
                            a, b,
                        ) // Return original expression if result is invalid
                    } else {

                        Expr::Constant(
                            result,
                        )
                    }
                },
                | _ => {
                    Expr::new_mul(a, b)
                }, /* Return original expression if conversion fails */
            }
        },
    }
}

#[inline]

pub(crate) fn div_em(
    a: &Expr,
    b: &Expr,
) -> Option<Expr> {

    if is_zero_expr(b) {

        // Division by zero - return appropriate representation
        if is_zero_expr(a) {

            // 0/0 is indeterminate
            return None;
        }

        // finite/0 approaches infinity
        return Some(Expr::Infinity);
    }

    match (a, b) {
        | (
            Expr::Constant(va),
            Expr::Constant(vb),
        ) => {

            let result = va / vb;

            if result.is_infinite() {

                Some(Expr::Infinity) // Use proper infinity representation
            } else if result.is_nan() {

                None // Undefined result like 0/0
            } else {

                Some(Expr::Constant(
                    result,
                ))
            }
        },
        // For integers, create a rational
        | (
            Expr::BigInt(ia),
            Expr::BigInt(ib),
        ) => {
            Some(Expr::Rational(
                BigRational::new(
                    ia.clone(),
                    ib.clone(),
                ),
            ))
        },
        | (
            Expr::Rational(ra),
            Expr::Rational(rb),
        ) => {
            Some(Expr::Rational(
                ra / rb,
            ))
        },
        | _ => {

            match (
                a.to_f64(),
                b.to_f64(),
            ) {
                | (
                    Some(va),
                    Some(vb),
                ) => {

                    let result =
                        va / vb;

                    if result
                        .is_infinite()
                    {

                        Some(Expr::Infinity)
                    // Use proper infinity representation
                    } else if result
                        .is_nan()
                    {

                        None // Undefined result like 0/0
                    } else {

                        Some(Expr::Constant(
                            result,
                        ))
                    }
                },
                | _ => {
                    Some(Expr::new_div(
                        a, b,
                    ))
                }, /* Return original expression if conversion fails */
            }
        },
    }
}

#[inline]

pub(crate) fn neg_em(a: &Expr) -> Expr {

    match a {
        | Expr::Constant(v) => {
            Expr::Constant(-v)
        },
        | Expr::BigInt(i) => {
            Expr::BigInt(-i)
        },
        | Expr::Rational(r) => {
            Expr::Rational(-r)
        },
        | _ => unreachable!(),
    }
}

// --- Helper Functions for Node Inspection ---

#[inline]

pub(crate) fn is_numeric_node(
    node: &Arc<DagNode>
) -> bool {

    matches!(
        &node.op,
        DagOp::Constant(_)
            | DagOp::BigInt(_)
            | DagOp::Rational(_)
    )
}

#[inline]

pub(crate) fn is_zero_expr(
    expr: &Expr
) -> bool {

    match expr {
        | Expr::Constant(c)
            if *c == 0.0 =>
        {
            true
        },
        | Expr::BigInt(i)
            if i.is_zero() =>
        {
            true
        },
        | Expr::Rational(r)
            if r.is_zero() =>
        {
            true
        },
        | _ => false, /* Default case returns false */
    }
}

#[inline]

pub(crate) fn is_one_expr(
    expr: &Expr
) -> bool {

    match expr {
        | Expr::Constant(c)
            if (*c - 1.0).abs()
                < f64::EPSILON =>
        {
            true
        },
        | Expr::BigInt(i)
            if i.is_one() =>
        {
            true
        },
        | Expr::Rational(r)
            if r.is_one() =>
        {
            true
        },
        | _ => false, /* Default case returns false */
    }
}

#[inline]

pub(crate) fn zero_node() -> Arc<DagNode>
{

    match DAG_MANAGER.get_or_create(
        &Expr::BigInt(BigInt::zero()),
    ) {
        | Ok(node) => node,
        | Err(_) => {

            // Fallback: create a constant 0 node if BigInt creation fails
            DAG_MANAGER
                .get_or_create(&Expr::Constant(0.0))
                .unwrap_or_else(|_| panic!("Failed to create zero node"))
        },
    }
}

#[inline]
#[allow(dead_code)]

pub(crate) fn one_node() -> Arc<DagNode>
{

    match DAG_MANAGER.get_or_create(
        &Expr::BigInt(BigInt::one()),
    ) {
        | Ok(node) => node,
        | Err(_) => {

            // Fallback: create a constant 1 node if BigInt creation fails
            DAG_MANAGER
                .get_or_create(&Expr::Constant(1.0))
                .unwrap_or_else(|_| panic!("Failed to create one node"))
        },
    }
}

#[inline]
/// Checks if a `DagNode` is a specific constant value.

pub(crate) fn is_const_node(
    node: &Arc<DagNode>,
    val: f64,
) -> bool {

    matches!(&node.op, DagOp::Constant(c) if (c.into_inner() - val).abs() < f64::EPSILON)
}

#[inline]
/// Checks if a `DagNode` represents the constant 0.

pub(crate) fn is_zero_node(
    node: &Arc<DagNode>
) -> bool {

    // Replaced matches! with a full match expression
    match &node.op {
        | DagOp::Constant(c)
            if c.is_zero() =>
        {
            true
        },
        | DagOp::BigInt(i)
            if i.is_zero() =>
        {
            true
        },
        | DagOp::Rational(r)
            if r.is_zero() =>
        {
            true
        },
        | _ => false, /* Default case returns false */
    }
}

#[inline]
/// Checks if a `DagNode` represents the constant 1.

pub(crate) fn is_one_node(
    node: &Arc<DagNode>
) -> bool {

    // Replaced matches! with a full match expression
    match &node.op {
        | DagOp::Constant(c)
            if c.is_one() =>
        {
            true
        },
        | DagOp::BigInt(i)
            if i.is_one() =>
        {
            true
        },
        | DagOp::Rational(r)
            if r.is_one() =>
        {
            true
        },
        | _ => false, /* Default case returns false */
    }
}

#[inline]
/// Checks if a `DagNode` represents the constant -1.

pub(crate) fn is_neg_one_node(
    node: &Arc<DagNode>
) -> bool {

    matches!(&node.op, DagOp::Constant(c) if (c.into_inner() + 1.0).abs() < f64::EPSILON)
}

#[inline]
/// Checks if a `DagNode` represents the constant Pi.

pub(crate) fn is_pi_node(
    node: &Arc<DagNode>
) -> bool {

    matches!(&node.op, DagOp::Pi)
}

#[inline]
/// Flattens nested Mul operations into a single list of factors.

pub(crate) fn flatten_mul_terms(
    node: &Arc<DagNode>,
    terms: &mut Vec<Arc<DagNode>>,
) {

    if matches!(&node.op, DagOp::Mul) {

        flatten_mul_terms(
            &node.children[0],
            terms,
        );

        flatten_mul_terms(
            &node.children[1],
            terms,
        );
    } else {

        terms.push(node.clone());
    }
}

/// Simplifies a Mul operation by flattening, collecting exponents, and rebuilding.

pub(crate) fn simplify_mul(
    node: &Arc<DagNode>
) -> Arc<DagNode> {

    // 1. Flatten the nested multiplications
    let mut factors = Vec::new();

    flatten_mul_terms(
        node,
        &mut factors,
    );

    // 2. Collect exponents and constant factor
    let mut exponents: BTreeMap<
        Arc<DagNode>,
        Expr,
    > = BTreeMap::new(); // base_node -> total_exponent_expr
    let mut constant =
        Expr::BigInt(BigInt::one());

    for factor in factors {

        if let Some(val) =
            get_numeric_value(&factor)
        {

            constant =
                mul_em(&constant, &val);

            continue;
        }

        let (base_node, exponent_expr) =
            if matches!(
                &factor.op,
                DagOp::Power
            ) {

                if factor
                    .children
                    .len()
                    < 2
                {

                    // Safety check: Power node should have 2 children
                    continue; // Skip malformed nodes
                }

                // Factor is Power(base, exp)
                (
                factor.children[0].clone(),
                factor.children[1]
                    .to_expr()
                    .unwrap_or(Expr::BigInt(
                        BigInt::one(),
                    )),
            )
            } else {

                // Factor is a variable or other expression, treat as factor^1
                (
                    factor.clone(),
                    Expr::BigInt(
                        BigInt::one(),
                    ),
                )
            };

        let entry = exponents
            .entry(base_node)
            .or_insert(Expr::BigInt(
                BigInt::zero(),
            ));

        *entry = add_em(
            entry,
            &exponent_expr,
        );
    }

    // 3. Rebuild the expression
    let mut new_factors = Vec::new();

    for (base, exponent) in exponents {

        if is_zero_expr(&exponent) {

            continue; // Skip terms with a zero exponent (x^0 = 1)
        }

        if is_one_expr(&exponent) {

            new_factors
                .push(base.clone()); // x^1 -> x
        } else {

            match DAG_MANAGER
                .get_or_create(
                    &exponent,
                ) {
                | Ok(exp_node) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Power,
                        vec![
                            base.clone(),
                            exp_node,
                        ],
                    ) {
                        | Ok(power_node) => new_factors.push(power_node),
                        | Err(_) => {

                            // If creating the power fails, just add the base without exponent
                            new_factors.push(base.clone());
                        },
                    }
                },
                | Err(_) => {

                    // If creating the exponent fails, just add the base without exponent
                    new_factors.push(
                        base.clone(),
                    );
                },
            }
        }
    }

    if is_zero_expr(&constant) {

        return zero_node(); // Multiplication by zero
    }

    if !is_one_expr(&constant) {

        if let Ok(constant_node) =
            DAG_MANAGER.get_or_create(
                &constant,
            )
        {

            new_factors.insert(
                0,
                constant_node,
            );
        }
    }

    if new_factors.is_empty() {

        return one_node();
    }

    // Build the final expression tree from the simplified factors
    new_factors.sort_by_key(|n| n.hash);

    let mut result =
        new_factors[0].clone();

    for factor in new_factors
        .iter()
        .skip(1)
    {

        result = match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Mul,
                vec![
                    result.clone(),
                    factor.clone(),
                ],
            ) {
            | Ok(node) => node,
            | Err(_) => {

                // If creating the multiplication fails, return the left operand
                break;
            },
        };
    }

    result
}

#[inline]
/// Flattens nested Add and Sub operations into a single list of terms.
/// Sub(a, b) is converted to Add(a, Neg(b)) for proper coefficient collection.

pub(crate) fn flatten_terms(
    node: &Arc<DagNode>,
    terms: &mut Vec<Arc<DagNode>>,
) {

    match &node.op {
        | DagOp::Add => {

            flatten_terms(
                &node.children[0],
                terms,
            );

            flatten_terms(
                &node.children[1],
                terms,
            );
        },
        | DagOp::Sub => {

            // a - b becomes a + (-b)
            if node.children.len() >= 2
            {

                flatten_terms(
                    &node.children[0],
                    terms,
                );

                // Add the negation of the second term
                match DAG_MANAGER.get_or_create_normalized(
                    DagOp::Neg,
                    vec![node.children[1].clone()],
                ) {
                    | Ok(neg_node) => terms.push(neg_node),
                    | Err(_) => {

                        // If negation fails, just push the original Sub node
                        terms.push(node.clone());
                    },
                }
            } else {

                terms
                    .push(node.clone());
            }
        },
        | _ => {

            terms.push(node.clone());
        },
    }
}

/// Simplifies an Add operation by flattening, collecting coefficients, and rebuilding.

pub(crate) fn simplify_add(
    node: &Arc<DagNode>
) -> Arc<DagNode> {

    // 1. Flatten the nested additions
    let mut terms = Vec::new();

    flatten_terms(node, &mut terms);

    // 2. Collect coefficients and constants
    let mut coeffs: BTreeMap<
        Arc<DagNode>,
        Expr,
    > = BTreeMap::new(); // base_node -> total_coeff_expr
    let mut constant =
        Expr::BigInt(BigInt::zero());

    for term in terms {

        if let Some(val) =
            get_numeric_value(&term)
        {

            constant =
                add_em(&constant, &val);

            continue;
        }

        // Simplify Mul nodes first to flatten nested multiplications
        let simplified_term = if matches!(
            &term.op,
            DagOp::Mul
        ) {

            simplify_mul(&term)
        } else {

            term.clone()
        };

        let (coeff_expr, base_node) =
            if matches!(
                &simplified_term.op,
                DagOp::Neg
            ) {

                // Neg(x) is treated as -1 * x
                if simplified_term
                    .children
                    .is_empty()
                {

                    (
                        Expr::BigInt(
                            BigInt::one(
                            ),
                        ),
                        simplified_term
                            .clone(),
                    )
                } else {

                    let child = &simplified_term.children[0];

                    // Check if child is Mul(c, x)
                    if matches!(
                        &child.op,
                        DagOp::Mul
                    ) && child
                        .children
                        .len()
                        >= 2
                    {

                        let c = &child
                            .children
                            [0];

                        let b = &child
                            .children
                            [1];

                        if is_numeric_node(c) {

                        get_numeric_value(c).map_or_else(
                            || (Expr::Constant(-1.0), child.clone()),
                            |val| (neg_em(&val), b.clone()),
                        )
                    } else {

                        (
                            Expr::Constant(-1.0),
                            child.clone(),
                        )
                    }
                    } else {

                        (
                        Expr::Constant(-1.0),
                        child.clone(),
                    )
                    }
                }
            } else if matches!(
                &simplified_term.op,
                DagOp::Mul
            ) {

                if simplified_term
                    .children
                    .len()
                    < 2
                {

                    (
                        Expr::BigInt(
                            BigInt::one(
                            ),
                        ),
                        simplified_term
                            .clone(),
                    )
                } else {

                    let c = &simplified_term.children[0];

                    let b = &simplified_term.children[1];

                    if is_numeric_node(c) {

                    get_numeric_value(c).map_or_else(
                        || (Expr::BigInt(BigInt::one()), simplified_term.clone()),
                        |val| (val, b.clone()),
                    )
                } else if is_numeric_node(b) {

                    get_numeric_value(b).map_or_else(
                        || (Expr::BigInt(BigInt::one()), simplified_term.clone()),
                        |val| (val, c.clone()),
                    )
                } else {

                    (
                        Expr::BigInt(BigInt::one()),
                        simplified_term.clone(),
                    )
                }
                }
            } else {

                (
                    Expr::BigInt(
                        BigInt::one(),
                    ),
                    simplified_term
                        .clone(),
                )
            };

        let entry = coeffs
            .entry(base_node)
            .or_insert(Expr::BigInt(
                BigInt::zero(),
            ));

        *entry =
            add_em(entry, &coeff_expr);
    }

    // 3. Rebuild the expression
    let mut new_terms = Vec::new();

    for (base, coeff) in coeffs {

        if is_zero_expr(&coeff) {

            continue; // Skip terms with a zero coefficient
        }

        if is_one_expr(&coeff) {

            new_terms
                .push(base.clone()); // 1*x -> x
        } else {

            match DAG_MANAGER
                .get_or_create(&coeff)
            {
                | Ok(coeff_node) => {

                    match DAG_MANAGER.get_or_create_normalized(
                        DagOp::Mul,
                        vec![
                            base.clone(),
                            coeff_node,
                        ],
                    ) {
                        | Ok(mul_node) => new_terms.push(mul_node),
                        | Err(_) => {

                            // If creating the multiplication fails, just add the base
                            new_terms.push(base.clone());
                        },
                    }
                },
                | Err(_) => {

                    // If creating the coefficient fails, just add the base
                    new_terms.push(
                        base.clone(),
                    );
                },
            }
        }
    }

    if !is_zero_expr(&constant) {

        if let Ok(constant_node) =
            DAG_MANAGER.get_or_create(
                &constant,
            )
        {

            new_terms
                .push(constant_node);
        }
    }

    if new_terms.is_empty() {

        return zero_node();
    }

    // Build the final expression tree from the simplified terms
    new_terms.sort_by_key(|n| n.hash);

    let mut result =
        new_terms[0].clone();

    for term in new_terms
        .iter()
        .skip(1)
    {

        result = match DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Add,
                vec![
                    result.clone(),
                    term.clone(),
                ],
            ) {
            | Ok(node) => node,
            | Err(_) => {

                // If creating the addition fails, return the left operand
                break;
            },
        };
    }

    result
}

// // Helper function to consolidate the complex logic of extracting a coefficient and its base.
// fn extract_coeff_and_base(node: &Arc<DagNode>) -> (Expr, Arc<DagNode>) {
// We expect the node here to be a term (Mul, Neg, or a plain variable/function).
// NOTE: The node passed here is already expected to be canonicalized by the caller
// (simplify_add), so we don't call simplify_mul here anymore.
//
// 1. Handle Negation (Term: -X or -(A*X))
// if matches!(&node.op, DagOp::Neg) {
// if node.children.is_empty() {
// Malformed Neg node.
// return (Expr::BigInt(BigInt::one()), node.clone());
// }
//
// let negated_child = &node.children[0];
//
// If we are negating a multiplication, we expect it to be in canonical form
// Mul(Constant, Base) or Mul(Term, Term)
// if matches!(&negated_child.op, DagOp::Mul) {
//
// Try to extract coefficient from the canonicalized Mul node
// if negated_child.children.len() >= 2 {
// let c = &negated_child.children[0];
// let b = &negated_child.children[1];
//
// Check if the first child is the numeric coefficient (canonical form assumption)
// if is_numeric_node(c) {
// if let Some(val) = get_numeric_value(c) {
// Negate the coefficient and use the other child as the base
// let negated_coeff = neg_em(&val);
// return (negated_coeff, b.clone());
// }
// }
// }
//
// Fallback: It's a Mul with no constant or a complex Mul, treat as Neg(X)
// Coeff: -1, Base: the Mul node
// return (Expr::Constant(-1.0), negated_child.clone());
//
// } else {
// Standard Neg(x) -> Coeff: -1, Base: x
// return (Expr::Constant(-1.0), negated_child.clone());
// }
// }
//
// 2. Handle Multiplication (Term: A*X or X*A)
// else if matches!(&node.op, DagOp::Mul) {
// if node.children.len() < 2 {
// Malformed Mul node, treat as is
// return (Expr::BigInt(BigInt::one()), node.clone());
// }
//
// let c = &node.children[0]; // Child 0
// let b = &node.children[1]; // Child 1
//
// The node is canonicalized, so the constant must be the first child (c)
// if is_numeric_node(c) {
// if let Some(val) = get_numeric_value(c) {
// Coeff: A, Base: X
// return (val, b.clone());
// }
// }
//
// If the first term is not numeric, then it's a complex product X*Y.
// Fallthrough to default case is intended.
// }
//
// 3. Default Case: Term is a plain variable (X, sin(y), exp(x)*cos(y), etc.)
// Coeff: 1, Base: the node itself (already canonicalized by caller)
// (Expr::BigInt(BigInt::one()), node.clone())
// }
//
// Simplifies an Add operation by flattening, collecting coefficients, and rebuilding.
// pub(crate) fn simplify_add(node: &Arc<DagNode>) -> Arc<DagNode> {
// 1. Flatten the nested additions
// let mut terms = Vec::new();
// flatten_terms(node, &mut terms);
//
// 2. Collect coefficients and constants
// let mut coeffs: BTreeMap<u64, (Arc<DagNode>, Expr)> = BTreeMap::new(); // base_hash -> (base_node, total_coeff_expr)
// let mut constant = Expr::BigInt(BigInt::zero());
//
// for term in terms {
// if let Some(val) = get_numeric_value(&term) {
// constant = add_em(&constant, &val);
// continue;
// }
//
// --- NEW: Canonicalize the term BEFORE extraction ---
// This is the critical step to ensure that (e^x * cos(y)) always results in the same hash.
// let canonical_term = simplify_mul(&term);
// --- END NEW ---
//
// --- Use the new helper function for extraction ---
// The node passed here is now guaranteed to be in its canonical form.
// let (coeff_expr, base_node) = extract_coeff_and_base(&canonical_term);
// --- End helper function call ---
//
// let entry = coeffs
// .entry(base_node.hash)
// .or_insert((base_node, Expr::BigInt(BigInt::zero())));
// entry.1 = add_em(&entry.1, &coeff_expr);
// }
//
// 3. Rebuild the expression (The rest of this function remains unchanged)
// let mut new_terms = Vec::new();
// for (_, (base, coeff)) in coeffs {
// if is_zero_expr(&coeff) {
// continue; // Skip terms with a zero coefficient
// }
// if is_one_expr(&coeff) {
// new_terms.push(base.clone()); // 1*x -> x
// } else {
// match DAG_MANAGER.get_or_create(&coeff) {
// Ok(coeff_node) => {
// match DAG_MANAGER
// The order here matters for canonical form: Base * Coeff, but simplify_mul
// should handle re-ordering to Coeff * Base if needed.
// .get_or_create_normalized(DagOp::Mul, vec![base.clone(), coeff_node])
// {
// Ok(mul_node) => new_terms.push(mul_node),
// Err(_) => {
// If creating the multiplication fails, just add the base
// new_terms.push(base.clone());
// }
// }
// }
// Err(_) => {
// If creating the coefficient fails, just add the base
// new_terms.push(base.clone());
// }
// }
// }
// }
//
// if !is_zero_expr(&constant) {
// match DAG_MANAGER.get_or_create(&constant) {
// Ok(constant_node) => new_terms.push(constant_node),
// Err(_) => {
// If creating the constant fails, skip it (equivalent to adding 0)
// }
// }
// }
//
// if new_terms.is_empty() {
// return zero_node();
// }
//
// Build the final expression tree from the simplified terms
// new_terms.sort_by_key(|n| n.hash);
// let mut result = new_terms[0].clone();
// for i in 1..new_terms.len() {
// result = match DAG_MANAGER
// .get_or_create_normalized(DagOp::Add, vec![result.clone(), new_terms[i].clone()])
// {
// Ok(node) => node,
// Err(_) => {
// If creating the addition fails, return the left operand
// break;
// }
// };
// }
//
// result
// }

// --- Pattern Matching and Substitution ---

/// Attempts to match an expression against a pattern.
///
/// If a match is found, it returns a `HashMap` containing the assignments
/// for the pattern variables. Pattern variables are represented by `Expr::Pattern(name)`.
///
/// This function handles `Expr::Dag` nodes by automatically unwrapping them to their
/// underlying structure for matching.
#[must_use]

pub fn pattern_match(
    expr: &Expr,
    pattern: &Expr,
) -> Option<HashMap<String, Expr>> {

    let mut assignments =
        HashMap::new();

    if pattern_match_recursive(
        expr,
        pattern,
        &mut assignments,
    ) {

        Some(assignments)
    } else {

        None
    }
}

/// DEBT: Rewrite it with Iteration
/// Recursively attempts to match an expression against a pattern.

pub(crate) fn pattern_match_recursive(
    expr: &Expr,
    pattern: &Expr,
    assignments: &mut HashMap<
        String,
        Expr,
    >,
) -> bool {

    // Unwrap DAG nodes for structural matching
    let expr_unwrapped = match expr {
        | Expr::Dag(node) => {
            node.to_expr()
                .unwrap_or_else(|_| {

                    expr.clone()
                })
        },
        | _ => expr.clone(),
    };

    let pattern_unwrapped =
        match pattern {
            | Expr::Dag(node) => {
                node.to_expr()
                    .unwrap_or_else(
                        |_| {

                            pattern
                                .clone()
                        },
                    )
            },
            | _ => pattern.clone(),
        };

    match (
        &expr_unwrapped,
        &pattern_unwrapped,
    ) {
        | (_, Expr::Pattern(name)) => {

            if let Some(existing) =
                assignments.get(name)
            {

                return existing
                    == expr;
            }

            assignments.insert(
                name.clone(),
                expr.clone(),
            );

            true
        },
        | (
            Expr::Add(e1, e2),
            Expr::Add(p1, p2),
        )
        | (
            Expr::Mul(e1, e2),
            Expr::Mul(p1, p2),
        ) => {

            let original_assignments =
                assignments.clone();

            if pattern_match_recursive(e1, p1, assignments)
                && pattern_match_recursive(e2, p2, assignments)
            {

                return true;
            }

            *assignments =
                original_assignments;

            pattern_match_recursive(
                e1,
                p2,
                assignments,
            ) && pattern_match_recursive(
                e2,
                p1,
                assignments,
            )
        },
        | (
            Expr::Sub(e1, e2),
            Expr::Sub(p1, p2),
        )
        | (
            Expr::Div(e1, e2),
            Expr::Div(p1, p2),
        )
        | (
            Expr::Power(e1, e2),
            Expr::Power(p1, p2),
        ) => {
            pattern_match_recursive(
                e1,
                p1,
                assignments,
            ) && pattern_match_recursive(
                e2,
                p2,
                assignments,
            )
        },
        | (
            Expr::Sin(e),
            Expr::Sin(p),
        )
        | (
            Expr::Cos(e),
            Expr::Cos(p),
        )
        | (
            Expr::Tan(e),
            Expr::Tan(p),
        )
        | (
            Expr::Exp(e),
            Expr::Exp(p),
        )
        | (
            Expr::Log(e),
            Expr::Log(p),
        )
        | (
            Expr::Neg(e),
            Expr::Neg(p),
        )
        | (
            Expr::Abs(e),
            Expr::Abs(p),
        )
        | (
            Expr::Sqrt(e),
            Expr::Sqrt(p),
        ) => {
            pattern_match_recursive(
                e,
                p,
                assignments,
            )
        },

        | (
            Expr::NaryList(s1, v1),
            Expr::NaryList(s2, v2),
        ) => {

            if s1 != s2
                || v1.len() != v2.len()
            {

                return false;
            }

            let original_assignments =
                assignments.clone();

            for (e, p) in v1
                .iter()
                .zip(v2.iter())
            {

                if !pattern_match_recursive(e, p, assignments) {

                    *assignments = original_assignments;

                    return false;
                }
            }

            true
        },
        | (
            Expr::UnaryList(s1, e1),
            Expr::UnaryList(s2, p1),
        ) => s1 == s2
            && pattern_match_recursive(
                e1,
                p1,
                assignments,
            ),
        | (
            Expr::BinaryList(
                s1,
                e1a,
                e1b,
            ),
            Expr::BinaryList(
                s2,
                p1a,
                p1b,
            ),
        ) => s1 == s2
            && pattern_match_recursive(
                e1a,
                p1a,
                assignments,
            )
            && pattern_match_recursive(
                e1b,
                p1b,
                assignments,
            ),

        | _ => {
            expr_unwrapped
                == pattern_unwrapped
        },
    }
}

#[must_use]

/// Substitutes pattern variables in a template expression with assigned expressions (DAG version).
///
/// This version handles `Expr::Dag` nodes by automatically unwrapping them for matching
/// and substitution, making it suitable for use with expressions stored in a Directed Acyclic Graph.
///
/// # Arguments
/// * `template` - The expression containing patterns to be replaced.
/// * `assignments` - A map from pattern names to their replacement expressions.
///
/// # Returns
/// A new `Expr` with patterns substituted.

pub fn substitute_patterns<
    S: std::hash::BuildHasher,
>(
    template: &Expr,
    assignments: &HashMap<
        String,
        Expr,
        S,
    >,
) -> Expr {

    let template_unwrapped =
        match template {
            | Expr::Dag(node) => {
                node.to_expr()
                    .unwrap_or_else(
                        |_| {

                            template
                                .clone()
                        },
                    )
            },
            | _ => template.clone(),
        };

    match template_unwrapped {
        | Expr::Pattern(name) => {
            assignments
                .get(&name)
                .cloned()
                .unwrap_or_else(|| template.clone())
        },
        | Expr::Add(a, b) => {
            Expr::new_add(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                substitute_patterns(&a, assignments),
                substitute_patterns(&b, assignments),
            )
        },
        | Expr::Power(b, e) => {
            Expr::new_pow(
                substitute_patterns(&b, assignments),
                substitute_patterns(&e, assignments),
            )
        },
        | Expr::Sin(a) => {
            Expr::new_sin(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Cos(a) => {
            Expr::new_cos(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Tan(a) => {
            Expr::new_tan(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Exp(a) => {
            Expr::new_exp(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Log(a) => {
            Expr::new_log(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Neg(a) => {
            Expr::new_neg(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Abs(a) => {
            Expr::new_abs(substitute_patterns(
                &a,
                assignments,
            ))
        },
        | Expr::Sqrt(a) => {
            Expr::new_sqrt(substitute_patterns(
                &a,
                assignments,
            ))
        },

        | Expr::NaryList(s, v) => {
            Expr::NaryList(
                s,
                v.iter()
                    .map(|e| substitute_patterns(e, assignments))
                    .collect(),
            )
        },
        | Expr::UnaryList(s, e) => {
            Expr::UnaryList(
                s,
                Arc::new(substitute_patterns(
                    &e,
                    assignments,
                )),
            )
        },
        | Expr::BinaryList(s, a, b) => {
            Expr::BinaryList(
                s,
                Arc::new(substitute_patterns(
                    &a,
                    assignments,
                )),
                Arc::new(substitute_patterns(
                    &b,
                    assignments,
                )),
            )
        },
        | _ => template.clone(),
    }
}
