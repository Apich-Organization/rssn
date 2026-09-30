use std::collections::HashSet;

use super::enode::ENode;
use super::id::Id;
use crate::symbolic::core::Expr;

/// An equivalence class containing equivalent `ENode`s and upward parent pointers.
#[derive(Clone, Debug)]
pub struct EClass {
    /// The unique canonical ID of this equivalence class.
    pub id: Id,
    /// All equivalent nodes that belong to this class.
    pub nodes: Vec<ENode>,
    /// Upward parent edges: (parent_node, parent_eclass_id).
    /// Used by congruence closure to notify parents when child classes are unioned.
    pub parents: Vec<(ENode, Id)>,
    /// Cached free variables appearing in expressions of this class.
    pub free_vars: HashSet<String>,
    /// Flag indicating whether this class contains high-weight operations
    /// (Derivative, Integral, Solve, Ode, etc.) awaiting de-cocooning.
    pub has_high_weight_op: bool,
    /// Cached lowest cost representation.
    pub best_cost: u64,
    /// Best node in this class according to the complexity cost function.
    pub best_node: Option<ENode>,
    /// Optional known constant value for constant propagation / folding.
    pub constant_value: Option<Expr>,
}

impl EClass {
    /// Creates a new `EClass` with an initial node.
    pub fn new(id: Id, initial_node: ENode) -> Self {
        let is_high = initial_node.is_high_weight_op();
        let cost = initial_node.complexity_weight() as u64;

        let mut free_vars = HashSet::new();
        if let ENode::Variable(name) = &initial_node {
            free_vars.insert(name.clone());
        }

        Self {
            id,
            nodes: vec![initial_node.clone()],
            parents: Vec::new(),
            free_vars,
            has_high_weight_op: is_high,
            best_cost: cost,
            best_node: Some(initial_node),
            constant_value: None,
        }
    }

    /// Adds a node to the class if not already present.
    /// Returns `true` if a new node was added.
    pub fn add_node(&mut self, node: ENode) -> bool {
        if self.nodes.contains(&node) {
            return false;
        }

        if node.is_high_weight_op() {
            self.has_high_weight_op = true;
        }

        if let ENode::Variable(name) = &node {
            self.free_vars.insert(name.clone());
        }

        let cost = node.complexity_weight() as u64;
        if cost < self.best_cost || self.best_node.is_none() {
            self.best_cost = cost;
            self.best_node = Some(node.clone());
        }

        self.nodes.push(node);
        true
    }

    /// Merges another `EClass` into this one.
    pub fn merge(&mut self, other: Self) {
        for node in other.nodes {
            self.add_node(node);
        }
        for parent in other.parents {
            if !self.parents.contains(&parent) {
                self.parents.push(parent);
            }
        }
        self.free_vars.extend(other.free_vars);
        if other.has_high_weight_op {
            self.has_high_weight_op = true;
        }
        if other.best_cost < self.best_cost {
            self.best_cost = other.best_cost;
            self.best_node = other.best_node;
        }
        if self.constant_value.is_none() && other.constant_value.is_some() {
            self.constant_value = other.constant_value;
        }
    }

    /// Checks if this class contains a variable by name.
    #[must_use]
    pub fn contains_var(&self, var: &str) -> bool {
        self.free_vars.contains(var)
    }

    /// Returns the number of equivalent nodes in this class.
    #[inline(always)]
    #[must_use]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Returns whether this class is empty.
    #[inline(always)]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Returns whether this class contains a representation of zero.
    #[inline]
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.nodes.iter().any(ENode::is_zero)
    }

    /// Returns whether this class contains a representation of one.
    #[inline]
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.nodes.iter().any(ENode::is_one)
    }

    /// Returns whether this class contains a representation of -1.
    #[inline]
    #[must_use]
    pub fn is_neg_one(&self) -> bool {
        self.nodes.iter().any(ENode::is_neg_one)
    }
}
