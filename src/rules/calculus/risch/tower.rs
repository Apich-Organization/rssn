//! A tower of monomial extensions.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use super::elem::Elem;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// The variable `x` itself: `theta' = 1`.
    Primitive,
    /// `theta = log(u)`: `theta' = u'/u`.
    Log,
    /// `theta = exp(u)`: `theta' = u' * theta`.
    Exp,
}

/// Level `i` (1-based) is `theta_i`; level 0 is `Q`.
#[derive(Clone, Debug, Default)]
pub struct Tower {
    kinds: Vec<Kind>,
    /// `theta_i' = a_i * theta_i + b_i`, both of depth `i - 1`.
    a: Vec<Elem>,
    b: Vec<Elem>,
}

impl Tower {
    /// A tower holding only the variable `x`.
    pub fn with_variable() -> Self {
        Self {
            kinds: vec![Kind::Primitive],
            a: vec![Elem::zero_at(0)],
            b: vec![Elem::one_at(0)],
        }
    }

    pub fn push(
        &mut self,
        kind: Kind,
        a: Elem,
        b: Elem,
    ) {
        self.kinds.push(kind);
        self.a.push(a);
        self.b.push(b);
    }

    pub fn depth(&self) -> usize {
        self.kinds.len()
    }

    pub fn kind(
        &self,
        level: usize,
    ) -> Kind {
        self.kinds
            .get(level - 1)
            .copied()
            .unwrap_or(Kind::Primitive)
    }

    /// `(a, b)` with `theta_level' = a*theta_level + b`.
    pub fn monomial(
        &self,
        level: usize,
    ) -> (Elem, Elem) {
        (
            self.a
                .get(level - 1)
                .cloned()
                .unwrap_or_else(|| Elem::zero_at(0)),
            self.b
                .get(level - 1)
                .cloned()
                .unwrap_or_else(|| Elem::zero_at(0)),
        )
    }
}
