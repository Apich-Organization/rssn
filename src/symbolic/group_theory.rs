//! # Group Theory
//!
//! This module provides structures for representing groups and their representations.
//! It includes definitions for `GroupElement` and `Group`, along with methods for
//! group multiplication and inverse. It also supports `Representation`s of groups
//! as matrices and character computations.

use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Represents a group element.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GroupElement(pub Expr);

/// Represents a group with its multiplication table.
#[derive(Debug, Clone)]
pub struct Group {
    pub elements: Vec<GroupElement>,
    // Multiplication table: (element1, element2) -> result
    pub multiplication_table: HashMap<(GroupElement, GroupElement), GroupElement>,
    pub identity: GroupElement,
}

impl Group {
    /// Creates a new group.
    ///
    /// # Arguments
    /// * `elements` - A `Vec<GroupElement>` containing all elements of the group.
    /// * `multiplication_table` - A `HashMap` defining the group's binary operation.
    /// * `identity` - The identity element of the group.
    ///
    /// # Returns
    /// A new `Group` instance.
    pub fn new(
        elements: Vec<GroupElement>,
        multiplication_table: HashMap<(GroupElement, GroupElement), GroupElement>,
        identity: GroupElement,
    ) -> Self {
        Group {
            elements,
            multiplication_table,
            identity,
        }
    }

    /// Multiplies two group elements.
    ///
    /// # Arguments
    /// * `a` - The first `GroupElement`.
    /// * `b` - The second `GroupElement`.
    ///
    /// # Returns
    /// An `Option<GroupElement>` representing the product `a * b`, or `None` if the product is not defined in the table.
    pub fn multiply(&self, a: &GroupElement, b: &GroupElement) -> Option<GroupElement> {
        self.multiplication_table
            .get(&(a.clone(), b.clone()))
            .cloned()
    }

    /// Computes the inverse of a group element.
    ///
    /// # Arguments
    /// * `a` - The `GroupElement` to find the inverse of.
    ///
    /// # Returns
    /// An `Option<GroupElement>` representing the inverse of `a`, or `None` if not found.
    pub fn inverse(&self, a: &GroupElement) -> Option<GroupElement> {
        for x in &self.elements {
            if let Some(product) = self.multiply(a, x) {
                if product == self.identity {
                    return Some(x.clone());
                }
            }
        }
        None
    }
}

/// Represents a group representation.
#[derive(Debug, Clone)]
pub struct Representation {
    pub group_elements: Vec<GroupElement>,
    pub matrices: HashMap<GroupElement, Expr>, // Each element maps to a matrix (Expr::Matrix)
}

impl Representation {
    /// Creates a new representation.
    ///
    /// # Arguments
    /// * `group_elements` - A `Vec<GroupElement>` of the elements in the group.
    /// * `matrices` - A `HashMap` mapping each `GroupElement` to its corresponding matrix (`Expr::Matrix`).
    ///
    /// # Returns
    /// A new `Representation` instance.
    pub fn new(group_elements: Vec<GroupElement>, matrices: HashMap<GroupElement, Expr>) -> Self {
        Representation {
            group_elements,
            matrices,
        }
    }

    /// Checks if the representation is valid (homomorphism property).
    ///
    /// A representation `ρ` is valid if it preserves the group operation:
    /// `ρ(g1 * g2) = ρ(g1) * ρ(g2)` for all `g1, g2` in the group.
    ///
    /// # Arguments
    /// * `group` - The `Group` that this is a representation of.
    ///
    /// # Returns
    /// `true` if the homomorphism property holds, `false` otherwise.
    pub fn is_valid(&self, group: &Group) -> bool {
        for g1 in &self.group_elements {
            for g2 in &self.group_elements {
                if let (Some(m1), Some(m2), Some(g1g2)) = (
                    self.matrices.get(g1),
                    self.matrices.get(g2),
                    group.multiply(g1, g2),
                ) {
                    if let Some(m_g1g2) = self.matrices.get(&g1g2) {
                        let m1m2 = crate::symbolic::matrix::mul_matrices(m1, m2);
                        if m1m2 != *m_g1g2 {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }
}

/// Computes the character of a representation.
///
/// The character `χ(g)` of a group element `g` in a representation `ρ` is the trace
/// of the matrix `ρ(g)`. It is a class function (constant on conjugacy classes).
///
/// # Arguments
/// * `representation` - The `Representation` to compute the character for.
///
/// # Returns
/// A `HashMap` mapping each `GroupElement` to its character value (`Expr`).
pub fn character(representation: &Representation) -> HashMap<GroupElement, Expr> {
    let mut chars = HashMap::new();
    for (element, matrix) in &representation.matrices {
        if let Expr::Matrix(rows) = matrix {
            let mut trace_val = Expr::Constant(0.0);
            //for i in 0..rows.len() {
            for (i, _item) in rows.iter().enumerate() {
                if let Some(diag_element) = rows[i].get(i) {
                    trace_val = Expr::Add(Box::new(trace_val), Box::new(diag_element.clone()));
                }
            }
            chars.insert(element.clone(), trace_val);
        }
    }
    chars
}
