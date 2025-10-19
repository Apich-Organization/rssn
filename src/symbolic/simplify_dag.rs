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

use std::collections::HashMap;
use std::sync::Arc;

use super::core::{DagManager, DagNode, DagOp, Expr, DAG_MANAGER};
use num_rational::BigRational;
use num_bigint::BigInt;
use std::collections::BTreeMap;
use num_traits::{One, Zero};

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
pub fn simplify(expr: &Expr) -> Expr {
    // Get the initial root node of the DAG from the input expression.
    let mut root_node = match DAG_MANAGER.get_or_create(expr) {
        Ok(node) => node,
        Err(_) => return expr.clone(), // If creation fails, return the original expression.
    };

    // --- Fixpoint Iteration Loop ---
    // This loop continues until a full pass over the DAG results in no changes.
    loop {
        let (simplified_root, changed) = bottom_up_simplify_pass(root_node.clone());

        if !changed {
            // If no changes were made in the last pass, the expression is fully simplified.
            break Expr::Dag(simplified_root);
        }

        // Update the root node for the next iteration.
        root_node = simplified_root;
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
pub(crate) fn bottom_up_simplify_pass(root: Arc<DagNode>) -> (Arc<DagNode>, bool) {
    // `memo` stores the simplified version of each node encountered in this pass.
    // Key: hash of the original node, Value: the simplified node.
    let mut memo: HashMap<u64, Arc<DagNode>> = HashMap::new();
    // `work_stack` manages the nodes to be visited.
    let mut work_stack: Vec<Arc<DagNode>> = vec![root.clone()];
    // `visited` keeps track of nodes pushed to the stack to avoid cycles and redundant work.
    let mut visited: HashMap<u64, bool> = HashMap::new();

    let mut changed_in_pass = false;

    while let Some(node) = work_stack.pop() {
        // If the node is already simplified in this pass, skip it.
        if memo.contains_key(&node.hash) {
            continue;
        }

        let children_simplified = node
            .children
            .iter()
            .all(|child| memo.contains_key(&child.hash));

        if children_simplified {
            // --- All children are simplified, so we can now process this node ---

            // 1. Rebuild the node with its (already simplified) children.
            let new_children: Vec<Arc<DagNode>> = node
                .children
                .iter()
                .map(|child| memo.get(&child.hash).expect("DAG_MANAGER.get_or_create_normalized failed.").clone())
                .collect();

            let rebuilt_node =
                DAG_MANAGER.get_or_create_normalized(node.op.clone(), new_children).expect("DAG_MANAGER.get_or_create_normalized failed.");

            // 2. Apply simplification rules to the rebuilt node.
            let simplified_node = apply_rules(&rebuilt_node);

            // 3. Check if simplification changed the node.
            if simplified_node.hash != node.hash {
                changed_in_pass = true;
            }

            // 4. Memoize the result for the original node's hash.
            memo.insert(node.hash, simplified_node);
        } else {
            // --- Not all children are simplified yet ---

            // 1. Push the current node back onto the stack to be processed later.
            work_stack.push(node.clone());

            // 2. Push un-simplified children onto the stack.
            if visited.insert(node.hash, true).is_none() {
                for child in node.children.iter().rev() {
                     work_stack.push(child.clone());
                }
            }
        }
    }

    // The final simplified root is the one corresponding to the original root's hash.
    let new_root = memo.get(&root.hash).expect("DAG_MANAGER.get_or_create_normalized failed.").clone();
    (new_root, changed_in_pass)
}

/// Applies a comprehensive set of algebraic and trigonometric simplification rules.
/// This function is the core of the simplification engine, performing pattern matching
/// on a single `DagNode` whose children are assumed to be already simplified.
pub(crate) fn apply_rules(node: &Arc<DagNode>) -> Arc<DagNode> {
    // --- Constant Folding ---
    if let Some(folded) = fold_constants(node) {
        return folded;
    }

    match &node.op {
        // --- Arithmetic Rules ---
        DagOp::Add => {
            let lhs = &node.children[0];
            let rhs = &node.children[1];
            // x + 0 -> x
            if is_zero_node(rhs) { return lhs.clone(); }
            if is_zero_node(lhs) { return rhs.clone(); }
            // x + x -> 2*x
            if lhs.hash == rhs.hash { 
                let two = DAG_MANAGER.get_or_create(&Expr::Constant(2.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![two, lhs.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }

            // Coefficient Collection: ax + bx -> (a+b)x
            if let (DagOp::Mul, DagOp::Mul) = (&lhs.op, &rhs.op) {
                let a = &lhs.children[0];
                let x1 = &lhs.children[1];
                let b = &rhs.children[0];
                let x2 = &rhs.children[1];

                if x1.hash == x2.hash { // a*x + b*x
                    let a_plus_b = DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![a.clone(), b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                    return DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![a_plus_b, x1.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                }
                if a.hash == b.hash { // x*a + y*a -> (x+y)*a
                    let x_plus_y = DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![x1.clone(), x2.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                    return DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![x_plus_y, a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                }
            }

            // sin(x)^2 + cos(x)^2 -> 1
            if let (DagOp::Power, DagOp::Power) = (&lhs.op, &rhs.op) {
                if is_const_node(&lhs.children[1], 2.0) && is_const_node(&rhs.children[1], 2.0) { // Both are squared
                    if let (DagOp::Sin, DagOp::Cos) = (&lhs.children[0].op, &rhs.children[0].op) {
                        if lhs.children[0].children[0].hash == rhs.children[0].children[0].hash { // sin(arg) and cos(arg) with same arg
                            return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
                        }
                    }
                    if let (DagOp::Cos, DagOp::Sin) = (&lhs.children[0].op, &rhs.children[0].op) {
                        if lhs.children[0].children[0].hash == rhs.children[0].children[0].hash { // cos(arg) and sin(arg) with same arg
                            return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
                        }
                    }
                }
            }
			
			return simplify_add(node);
        }
        DagOp::Sub => {
            let lhs = &node.children[0];
            let rhs = &node.children[1];
            // x - 0 -> x
            if is_zero_node(rhs) { return lhs.clone(); }
            // x - x -> 0
            if lhs.hash == rhs.hash { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
        }
        DagOp::Mul => {
            let lhs = &node.children[0];
            let rhs = &node.children[1];
            // x * 1 -> x
            if is_one_node(rhs) { return lhs.clone(); }
            if is_one_node(lhs) { return rhs.clone(); }
            // x * 0 -> 0
            if is_zero_node(rhs) || is_zero_node(lhs) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // x * -1 -> -x
            if is_neg_one_node(rhs) { return DAG_MANAGER.get_or_create_normalized(DagOp::Neg, vec![lhs.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            if is_neg_one_node(lhs) { return DAG_MANAGER.get_or_create_normalized(DagOp::Neg, vec![rhs.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // x * x -> x^2
            if lhs.hash == rhs.hash {
                let two = DAG_MANAGER.get_or_create(&Expr::Constant(2.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Power, vec![lhs.clone(), two]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
            // Distributivity: a * (b + c) -> a*b + a*c
            if let DagOp::Add = &rhs.op {
                let a = lhs;
                let b = &rhs.children[0];
                let c = &rhs.children[1];
                let ab = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![a.clone(), b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let ac = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![a.clone(), c.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![ab, ac]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
            if let DagOp::Add = &lhs.op {
                let a = &lhs.children[0];
                let b = &lhs.children[1];
                let c = rhs;
                let ac = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![a.clone(), c.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let bc = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![b.clone(), c.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![ac, bc]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
			return simplify_mul(node);
        }
        DagOp::Div => {
            let lhs = &node.children[0];
            let rhs = &node.children[1];
            // x / 1 -> x
            if is_one_node(rhs) { return lhs.clone(); }
            // x / x -> 1 (if x != 0)
            if lhs.hash == rhs.hash { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // 0 / x -> 0 (if x != 0)
            if is_zero_node(lhs) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            let neg_one = DAG_MANAGER.get_or_create(&Expr::BigInt(BigInt::from(-1))).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let rhs_pow_neg_one = DAG_MANAGER.get_or_create_normalized(DagOp::Power, vec![rhs.clone(), neg_one]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![lhs.clone(), rhs_pow_neg_one]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }
        DagOp::Neg => {
            let inner = &node.children[0];
            // --x -> x
            if let DagOp::Neg = &inner.op {
                return inner.children[0].clone();
            }
        }

        // --- Power & Log/Exp Rules ---
        DagOp::Power => {
            let base = &node.children[0];
            let exp = &node.children[1];
            // x ^ 1 -> x
            if is_one_node(exp) { return base.clone(); }
            // x ^ 0 -> 1
            if is_zero_node(exp) { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // 1 ^ x -> 1
            if is_one_node(base) { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // (x^a)^b -> x^(a*b)
            if let DagOp::Power = &base.op {
                let new_exp = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![base.children[1].clone(), exp.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Power, vec![base.children[0].clone(), new_exp]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
        }
        DagOp::Log => {
            let arg = &node.children[0];
            // log(1) -> 0
            if is_one_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // log(e) -> 1
            if let DagOp::E = &arg.op { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // log(exp(x)) -> x
            if let DagOp::Exp = &arg.op { return arg.children[0].clone(); }
            // log(a*b) -> log(a) + log(b)
            if let DagOp::Mul = &arg.op {
                let log_a = DAG_MANAGER.get_or_create_normalized(DagOp::Log, vec![arg.children[0].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let log_b = DAG_MANAGER.get_or_create_normalized(DagOp::Log, vec![arg.children[1].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![log_a, log_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
            // log(a^b) -> b*log(a)
            if let DagOp::Power = &arg.op {
                let b = arg.children[1].clone();
                let log_a = DAG_MANAGER.get_or_create_normalized(DagOp::Log, vec![arg.children[0].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![b, log_a]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
        }
        DagOp::LogBase => {
            let base = &node.children[0];
            let arg = &node.children[1];
            // log_b(b) -> 1
            if base.hash == arg.hash { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // log_b(a) -> log(a) / log(b)
            let log_a = DAG_MANAGER.get_or_create_normalized(DagOp::Log, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let log_b = DAG_MANAGER.get_or_create_normalized(DagOp::Log, vec![base.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Div, vec![log_a, log_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }
        DagOp::Exp => {
            let arg = &node.children[0];
            // exp(0) -> 1
            if is_zero_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // exp(log(x)) -> x
            if let DagOp::Log = &arg.op { return arg.children[0].clone(); }
        }

        // --- Trigonometric Rules ---
        DagOp::Sin => {
            let arg = &node.children[0];
            // sin(0) -> 0
            if is_zero_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // sin(pi) -> 0
            if is_pi_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // sin(-x) -> -sin(x)
            if let DagOp::Neg = &arg.op {
                let new_sin = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![arg.children[0].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Neg, vec![new_sin]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
            // Sum/Difference and Induction formulas for sin
            if let DagOp::Add = &arg.op {
                let a = &arg.children[0];
                let b = &arg.children[1];
                // sin(x + pi) -> -sin(x)
                if is_pi_node(b) {
                    let sin_a = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                    return DAG_MANAGER.get_or_create_normalized(DagOp::Neg, vec![sin_a]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                }
                // sin(a+b) -> sin(a)cos(b) + cos(a)sin(b)
                let sin_a = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let cos_b = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let cos_a = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let sin_b = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let term1 = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![sin_a, cos_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let term2 = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![cos_a, sin_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![term1, term2]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
        }
        DagOp::Cos => {
            let arg = &node.children[0];
            // cos(0) -> 1
            if is_zero_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // cos(pi) -> -1
            if is_pi_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(-1.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // cos(-x) -> cos(x)
            if let DagOp::Neg = &arg.op {
                return DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![arg.children[0].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
            // Sum/Difference and Induction formulas for cos
            if let DagOp::Add = &arg.op {
                let a = &arg.children[0];
                let b = &arg.children[1];
                // cos(x + pi) -> -cos(x)
                if is_pi_node(b) {
                    let cos_a = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                    return DAG_MANAGER.get_or_create_normalized(DagOp::Neg, vec![cos_a]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                }
                // cos(a+b) -> cos(a)cos(b) - sin(a)sin(b)
                let cos_a = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let cos_b = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let sin_a = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![a.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let sin_b = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![b.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let term1 = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![cos_a, cos_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                let term2 = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![sin_a, sin_b]).expect("DAG_MANAGER.get_or_create_normalized failed.");
                return DAG_MANAGER.get_or_create_normalized(DagOp::Sub, vec![term1, term2]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            }
        }
        DagOp::Tan => {
            let arg = &node.children[0];
            // tan(0) -> 0
            if is_zero_node(arg) { return DAG_MANAGER.get_or_create(&Expr::Constant(0.0)).expect("DAG_MANAGER.get_or_create_normalized failed."); }
            // tan(x) -> sin(x)/cos(x)
            let sin_x = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let cos_x = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Div, vec![sin_x, cos_x]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }
        DagOp::Sec => { // sec(x) -> 1/cos(x)
            let arg = &node.children[0];
            let one = DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let cos_x = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Div, vec![one, cos_x]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }
        DagOp::Csc => { // csc(x) -> 1/sin(x)
            let arg = &node.children[0];
            let one = DAG_MANAGER.get_or_create(&Expr::Constant(1.0)).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let sin_x = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Div, vec![one, sin_x]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }
        DagOp::Cot => { // cot(x) -> cos(x)/sin(x)
            let arg = &node.children[0];
            let cos_x = DAG_MANAGER.get_or_create_normalized(DagOp::Cos, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            let sin_x = DAG_MANAGER.get_or_create_normalized(DagOp::Sin, vec![arg.clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
            return DAG_MANAGER.get_or_create_normalized(DagOp::Div, vec![cos_x, sin_x]).expect("DAG_MANAGER.get_or_create_normalized failed.");
        }

        _ => {} // No rule matched for this operator
    }

    node.clone()
}

/// Performs constant folding on a `DagNode`.
/// If the node is an operation on constant children, it computes the result and returns a new constant node.
pub(crate) fn fold_constants(node: &Arc<DagNode>) -> Option<Arc<DagNode>> {
    let children_values: Option<Vec<Expr>> = node.children.iter().map(get_numeric_value).collect();

    if let Some(values) = children_values {
        let result = match (&node.op, values.as_slice()) {
            (DagOp::Add, [a, b]) => Some(add_em(a, b)),
            (DagOp::Sub, [a, b]) => Some(sub_em(a, b)),
            (DagOp::Mul, [a, b]) => Some(mul_em(a, b)),
            (DagOp::Div, [a, b]) => div_em(a, b),
            (DagOp::Power, [Expr::Constant(a), Expr::Constant(b)]) => Some(Expr::Constant(a.powf(*b))),
            (DagOp::Neg, [a]) => Some(neg_em(a)),
            _ => None,
        };

        if let Some(value) = result {
            return Some(DAG_MANAGER.get_or_create(&value).expect("DAG_MANAGER.get_or_create_normalized failed."));
        }
    }

    None
}

// --- Numeric Helper Functions ---

#[inline]
pub(crate) fn get_numeric_value(node: &Arc<DagNode>) -> Option<Expr> {
    match &node.op {
        DagOp::Constant(c) => Some(Expr::Constant(c.into_inner())),
        DagOp::BigInt(i) => Some(Expr::BigInt(i.clone())),
        DagOp::Rational(r) => Some(Expr::Rational(r.clone())),
        _ => None,
    }
}
#[inline]
pub(crate) fn add_em(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        (Expr::Constant(va), Expr::Constant(vb)) => Expr::Constant(va + vb),
        (Expr::BigInt(ia), Expr::BigInt(ib)) => Expr::BigInt(ia + ib),
        (Expr::Rational(ra), Expr::Rational(rb)) => Expr::Rational(ra + rb),
        // Promote to Rational or Constant
        _ => Expr::Constant(a.to_f64().unwrap_or(0.0) + b.to_f64().unwrap_or(0.0))
    }
}
#[inline]
pub(crate) fn sub_em(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        (Expr::Constant(va), Expr::Constant(vb)) => Expr::Constant(va - vb),
        (Expr::BigInt(ia), Expr::BigInt(ib)) => Expr::BigInt(ia - ib),
        (Expr::Rational(ra), Expr::Rational(rb)) => Expr::Rational(ra - rb),
        _ => Expr::Constant(a.to_f64().unwrap_or(0.0) - b.to_f64().unwrap_or(0.0))
    }
}
#[inline]
pub(crate) fn mul_em(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        (Expr::Constant(va), Expr::Constant(vb)) => Expr::Constant(va * vb),
        (Expr::BigInt(ia), Expr::BigInt(ib)) => Expr::BigInt(ia * ib),
        (Expr::Rational(ra), Expr::Rational(rb)) => Expr::Rational(ra * rb),
        _ => Expr::Constant(a.to_f64().unwrap_or(1.0) * b.to_f64().unwrap_or(1.0))
    }
}
#[inline]
pub(crate) fn div_em(a: &Expr, b: &Expr) -> Option<Expr> {
    if is_zero_expr(b) { return None; }
    match (a, b) {
        (Expr::Constant(va), Expr::Constant(vb)) => Some(Expr::Constant(va / vb)),
        // For integers, create a rational
        (Expr::BigInt(ia), Expr::BigInt(ib)) => Some(Expr::Rational(BigRational::new(ia.clone(), ib.clone()))),
        (Expr::Rational(ra), Expr::Rational(rb)) => Some(Expr::Rational(ra / rb)),
        _ => Some(Expr::Constant(a.to_f64().unwrap_or(0.0) / b.to_f64().unwrap_or(1.0)))
    }
}
#[inline]
pub(crate) fn neg_em(a: &Expr) -> Expr {
    match a {
        Expr::Constant(v) => Expr::Constant(-v),
        Expr::BigInt(i) => Expr::BigInt(-i),
        Expr::Rational(r) => Expr::Rational(-r),
        _ => unreachable!()
    }
}

// --- Helper Functions for Node Inspection ---

#[inline]
pub(crate) fn is_numeric_node(node: &Arc<DagNode>) -> bool {
    matches!(&node.op, DagOp::Constant(_) | DagOp::BigInt(_) | DagOp::Rational(_))
}
#[inline]
pub(crate) fn is_zero_expr(expr: &Expr) -> bool {
    match expr {
        Expr::Constant(c) if *c == 0.0 => true,
        Expr::BigInt(i) if i.is_zero() => true,
        Expr::Rational(r) if r.is_zero() => true,
        _ => false, // Default case returns false
    }
}
#[inline]
pub(crate) fn is_one_expr(expr: &Expr) -> bool {
    match expr {
        Expr::Constant(c) if *c == 1.0 => true,
        Expr::BigInt(i) if i.is_one() => true,
        Expr::Rational(r) if r.is_one() => true,
        _ => false, // Default case returns false
    }
}
#[inline]
pub(crate) fn zero_node() -> Arc<DagNode> {
    DAG_MANAGER.get_or_create(&Expr::BigInt(BigInt::zero())).expect("DAG_MANAGER.get_or_create_normalized failed.")
}
#[inline]
#[allow(dead_code)]
pub(crate) fn one_node() -> Arc<DagNode> {
    DAG_MANAGER.get_or_create(&Expr::BigInt(BigInt::one())).expect("DAG_MANAGER.get_or_create_normalized failed.")
}
#[inline]
/// Checks if a `DagNode` is a specific constant value.
pub(crate) fn is_const_node(node: &Arc<DagNode>, val: f64) -> bool {
    matches!(&node.op, DagOp::Constant(c) if c.into_inner() == val)
}
#[inline]
/// Checks if a `DagNode` represents the constant 0.
pub(crate) fn is_zero_node(node: &Arc<DagNode>) -> bool {
    // Replaced matches! with a full match expression
    match &node.op {
        DagOp::Constant(c) if c.is_zero() => true,
        DagOp::BigInt(i) if i.is_zero() => true,
        DagOp::Rational(r) if r.is_zero() => true,
        _ => false, // Default case returns false
    }
}
#[inline]
/// Checks if a `DagNode` represents the constant 1.
pub(crate) fn is_one_node(node: &Arc<DagNode>) -> bool {
    // Replaced matches! with a full match expression
    match &node.op {
        DagOp::Constant(c) if c.is_one() => true,
        DagOp::BigInt(i) if i.is_one() => true,
        DagOp::Rational(r) if r.is_one() => true,
        _ => false, // Default case returns false
    }
}
#[inline]
/// Checks if a `DagNode` represents the constant -1.
pub(crate) fn is_neg_one_node(node: &Arc<DagNode>) -> bool {
    matches!(&node.op, DagOp::Constant(c) if c.into_inner() == -1.0)
}
#[inline]
/// Checks if a `DagNode` represents the constant Pi.
pub(crate) fn is_pi_node(node: &Arc<DagNode>) -> bool {
    matches!(&node.op, DagOp::Pi)
}
#[inline]
/// Flattens nested Mul operations into a single list of factors.
pub(crate) fn flatten_mul_terms(node: &Arc<DagNode>, terms: &mut Vec<Arc<DagNode>>) {
    if let DagOp::Mul = &node.op {
        flatten_mul_terms(&node.children[0], terms);
        flatten_mul_terms(&node.children[1], terms);
    } else {
        terms.push(node.clone());
    }
}

/// Simplifies a Mul operation by flattening, collecting exponents, and rebuilding.
pub(crate) fn simplify_mul(node: &Arc<DagNode>) -> Arc<DagNode> {
    // 1. Flatten the nested multiplications
    let mut factors = Vec::new();
    flatten_mul_terms(node, &mut factors);

    // 2. Collect exponents and constant factor
    let mut exponents: BTreeMap<u64, (Arc<DagNode>, Expr)> = BTreeMap::new(); // base_hash -> (base_node, total_exponent_expr)
    let mut constant = Expr::BigInt(BigInt::one());

    for factor in factors {
        if let Some(val) = get_numeric_value(&factor) {
            constant = mul_em(&constant, &val);
            continue;
        }

        let (base_node, exponent_expr) = if let DagOp::Power = &factor.op {
            // Factor is Power(base, exp)
            (factor.children[0].clone(), get_numeric_value(&factor.children[1]).unwrap_or(Expr::BigInt(BigInt::one())))
        } else {
            // Factor is a variable or other expression, treat as factor^1
            (factor.clone(), Expr::BigInt(BigInt::one()))
        };

        let entry = exponents.entry(base_node.hash).or_insert((base_node, Expr::BigInt(BigInt::zero())));
        entry.1 = add_em(&entry.1, &exponent_expr);
    }

    // 3. Rebuild the expression
    let mut new_factors = Vec::new();
    for (_, (base, exponent)) in exponents {
        if is_zero_expr(&exponent) {
            continue; // Skip terms with a zero exponent (x^0 = 1)
        }
        if is_one_expr(&exponent) {
            new_factors.push(base); // x^1 -> x
        } else {
            let exp_node = DAG_MANAGER.get_or_create(&exponent).expect("DAG_MANAGER.get_or_create_normalized failed.");
            new_factors.push(DAG_MANAGER.get_or_create_normalized(DagOp::Power, vec![base, exp_node]).expect("DAG_MANAGER.get_or_create_normalized failed."));
        }
    }

    if is_zero_expr(&constant) {
        return zero_node(); // Multiplication by zero
    }

    if !is_one_expr(&constant) {
        new_factors.insert(0, DAG_MANAGER.get_or_create(&constant).expect("DAG_MANAGER.get_or_create_normalized failed."));
    }

    if new_factors.is_empty() {
        return one_node();
    }

    // Build the final expression tree from the simplified factors
    new_factors.sort_by_key(|n| n.hash);
    let mut result = new_factors[0].clone();
    for i in 1..new_factors.len() {
        result = DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![result, new_factors[i].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
    }
    
    result
}
#[inline]
/// Flattens nested Add operations into a single list of terms.
pub(crate) fn flatten_terms(node: &Arc<DagNode>, terms: &mut Vec<Arc<DagNode>>) {
    if let DagOp::Add = &node.op {
        flatten_terms(&node.children[0], terms);
        flatten_terms(&node.children[1], terms);
    } else {
        terms.push(node.clone());
    }
}

/// Simplifies an Add operation by flattening, collecting coefficients, and rebuilding.
pub(crate) fn simplify_add(node: &Arc<DagNode>) -> Arc<DagNode> {
    // 1. Flatten the nested additions
    let mut terms = Vec::new();
    flatten_terms(node, &mut terms);

    // 2. Collect coefficients and constants
    let mut coeffs: BTreeMap<u64, (Arc<DagNode>, Expr)> = BTreeMap::new(); // base_hash -> (base_node, total_coeff_expr)
    let mut constant = Expr::BigInt(BigInt::zero());

    for term in terms {
        if let Some(val) = get_numeric_value(&term) {
            constant = add_em(&constant, &val);
            continue;
        }

        let (coeff_expr, base_node) = if let DagOp::Mul = &term.op {
            let c = &term.children[0];
            let b = &term.children[1];
            if is_numeric_node(c) {
                (get_numeric_value(c).expect("DAG_MANAGER.get_or_create_normalized failed."), b.clone())
            } else if is_numeric_node(b) {
                (get_numeric_value(b).expect("DAG_MANAGER.get_or_create_normalized failed."), c.clone())
            } else {
                (Expr::BigInt(BigInt::one()), term.clone())
            }
        } else {
            (Expr::BigInt(BigInt::one()), term.clone())
        };

        let entry = coeffs.entry(base_node.hash).or_insert((base_node, Expr::BigInt(BigInt::zero())));
        entry.1 = add_em(&entry.1, &coeff_expr);
    }

    // 3. Rebuild the expression
    let mut new_terms = Vec::new();
    for (_, (base, coeff)) in coeffs {
        if is_zero_expr(&coeff) {
            continue; // Skip terms with a zero coefficient
        }
        if is_one_expr(&coeff) {
            new_terms.push(base); // 1*x -> x
        } else {
            let coeff_node = DAG_MANAGER.get_or_create(&coeff).expect("DAG_MANAGER.get_or_create_normalized failed.");
            // To match test (x*2), put coeff second.
            new_terms.push(DAG_MANAGER.get_or_create_normalized(DagOp::Mul, vec![base, coeff_node]).expect("DAG_MANAGER.get_or_create_normalized failed."));
        }
    }

    if !is_zero_expr(&constant) {
        new_terms.push(DAG_MANAGER.get_or_create(&constant).expect("DAG_MANAGER.get_or_create_normalized failed."));
    }

    if new_terms.is_empty() {
        return zero_node();
    }

    // Build the final expression tree from the simplified terms
    new_terms.sort_by_key(|n| n.hash);
    let mut result = new_terms[0].clone();
    for i in 1..new_terms.len() {
        result = DAG_MANAGER.get_or_create_normalized(DagOp::Add, vec![result, new_terms[i].clone()]).expect("DAG_MANAGER.get_or_create_normalized failed.");
    }
    
    result
}