//! # Group Theory
//!
//! This module provides structures for representing groups and their representations.
//! It includes definitions for `GroupElement` and `Group`, along with methods for
//! group multiplication and inverse. It also supports `Representation`s of groups
//! as matrices and character computations.

use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Represents a group element.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]

pub struct GroupElement(pub Expr);

/// Represents a group with its multiplication table.
#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct Group {
    pub elements: Vec<GroupElement>,
    pub multiplication_table: HashMap<(GroupElement, GroupElement), GroupElement>,
    pub identity: GroupElement,
}

impl Group {
    /// Creates a new group.
    #[must_use]

    pub const fn new(
        elements: Vec<GroupElement>,
        multiplication_table: HashMap<(GroupElement, GroupElement), GroupElement>,
        identity: GroupElement,
    ) -> Self {

        Self {
            elements,
            multiplication_table,
            identity,
        }
    }

    /// Multiplies two group elements.
    #[must_use]

    pub fn multiply(&self, a: &GroupElement, b: &GroupElement) -> Option<GroupElement> {

        self.multiplication_table
            .get(&(a.clone(), b.clone()))
            .cloned()
    }

    /// Computes the inverse of a group element.
    #[must_use]

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

    /// Checks if the group is abelian (commutative).
    #[must_use]

    pub fn is_abelian(&self) -> bool {

        for a in &self.elements {

            for b in &self.elements {

                let ab = self.multiply(a, b);

                let ba = self.multiply(b, a);

                if ab != ba {

                    return false;
                }
            }
        }

        true
    }

    /// Computes the order of an element g (smallest k such that g^k = e).
    #[must_use]

    pub fn element_order(&self, g: &GroupElement) -> Option<usize> {

        let mut current = g.clone();

        // The order of an element must divide the group order (Lagrange's theorem),
        // so we check up to group size.
        for k in 1..=self.elements.len() {

            if current == self.identity {

                return Some(k);
            }

            current = self.multiply(&current, g)?;
        }

        None
    }

    /// Finds the conjugacy classes of the group.
    #[must_use]

    pub fn conjugacy_classes(&self) -> Vec<Vec<GroupElement>> {

        let mut classes: Vec<Vec<GroupElement>> = Vec::new();

        let mut visited: Vec<GroupElement> = Vec::new();

        for x in &self.elements {

            if visited.contains(x) {

                continue;
            }

            let mut class = Vec::new();

            for g in &self.elements {

                // g * x * g^-1
                if let Some(g_inv) = self.inverse(g) {

                    if let Some(gx) = self.multiply(g, x) {

                        if let Some(conjugate) = self.multiply(&gx, &g_inv) {

                            if !class.contains(&conjugate) {

                                class.push(conjugate.clone());
                            }
                        }
                    }
                }
            }

            for c in &class {

                visited.push(c.clone());
            }

            classes.push(class);
        }

        classes
    }

    /// Finds the center of the group Z(G) = {z in G | zg = gz for all g in G}.
    #[must_use]

    pub fn center(&self) -> Vec<GroupElement> {

        let mut center_elements = Vec::new();

        for z in &self.elements {

            let mut commutes_with_all = true;

            for g in &self.elements {

                let zg = self.multiply(z, g);

                let gz = self.multiply(g, z);

                if zg != gz {

                    commutes_with_all = false;

                    break;
                }
            }

            if commutes_with_all {

                center_elements.push(z.clone());
            }
        }

        center_elements
    }
}

/// Represents a group representation.
#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct Representation {
    pub group_elements: Vec<GroupElement>,
    pub matrices: HashMap<GroupElement, Expr>,
}

impl Representation {
    /// Creates a new representation.
    #[must_use]

    pub const fn new(
        group_elements: Vec<GroupElement>,
        matrices: HashMap<GroupElement, Expr>,
    ) -> Self {

        Self {
            group_elements,
            matrices,
        }
    }

    /// Checks if the representation is valid (homomorphism property).
    #[must_use]

    pub fn is_valid(&self, group: &Group) -> bool {

        for g1 in &self.group_elements {

            for g2 in &self.group_elements {

                if let (Some(m1), Some(m2), Some(g1g2)) = (
                    self.matrices
                        .get(g1),
                    self.matrices
                        .get(g2),
                    group.multiply(g1, g2),
                ) {

                    if let Some(m_g1g2) = self
                        .matrices
                        .get(&g1g2)
                    {

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

use crate::symbolic::simplify_dag::simplify;

/// Computes the character of a representation.
#[must_use]

pub fn character(representation: &Representation) -> HashMap<GroupElement, Expr> {

    let mut chars = HashMap::new();

    for (element, matrix) in &representation.matrices {

        if let Expr::Matrix(rows) = matrix {

            let mut trace_val = Expr::Constant(0.0);

            for (i, _item) in rows
                .iter()
                .enumerate()
            {

                if let Some(diag_element) = rows[i].get(i) {

                    trace_val = simplify(&Expr::new_add(trace_val, diag_element.clone()));
                }
            }

            chars.insert(element.clone(), trace_val);
        }
    }

    chars
}
