//! Deterministic infix rendering of concrete terms.
//!
//! The output re-sugars the canonical `add` / `mul` / `pow` forms into
//! subtraction and division and is parseable by [`Graph::parse`]. Children
//! of commutative operators are printed in a canonical order that does not
//! depend on node creation order, so the text is stable across runs.

use std::cmp::Ordering;
use std::fmt::Write;

use super::id::NodeId;
use super::number::Number;
use super::op::OpFlags;
use super::op::core;
use super::op::Arity;
use super::payload::Payload;
use super::store::Graph;

pub(crate) const PREC_EQ: u8 = 0;
pub(crate) const PREC_ADD: u8 = 1;
pub(crate) const PREC_MUL: u8 = 2;
pub(crate) const PREC_NEG: u8 = 3;
pub(crate) const PREC_POW: u8 = 4;
pub(crate) const PREC_ATOM: u8 = 5;

impl Graph {
    /// A total order on concrete terms that is independent of node ids:
    /// numbers, then symbols by name, then compound terms by operator name
    /// and children.
    #[must_use]
    pub fn term_cmp(
        &self,
        a: NodeId,
        b: NodeId,
    ) -> Ordering {
        if a == b {
            return Ordering::Equal;
        }
        let rank = |n: NodeId| match self.payload(n) {
            | Some(Payload::Num(_)) => 0_u8,
            | Some(Payload::Bool(_)) => 1,
            | Some(Payload::Str(_)) => 2,
            | Some(Payload::Sym(_)) => 3,
            | Some(Payload::Blob(_)) => 4,
            | None => 5,
        };
        rank(a)
            .cmp(&rank(b))
            .then_with(|| match (self.payload(a), self.payload(b)) {
                | (Some(Payload::Num(x)), Some(Payload::Num(y))) => x.total_cmp(y),
                | (Some(Payload::Bool(x)), Some(Payload::Bool(y))) => x.cmp(y),
                | (Some(Payload::Str(x)), Some(Payload::Str(y))) => x.cmp(y),
                | (Some(Payload::Sym(x)), Some(Payload::Sym(y))) => self
                    .interner()
                    .symbol_name(*x)
                    .cmp(self.interner().symbol_name(*y)),
                | (Some(_), Some(_)) => a.cmp(&b),
                | _ => {
                    let (ca, cb) = (self.ordered_children(a), self.ordered_children(b));
                    self.ops()
                        .get(self.op(a))
                        .name
                        .cmp(&self.ops().get(self.op(b)).name)
                        .then_with(|| {
                            ca.iter()
                                .zip(&cb)
                                .map(|(&x, &y)| self.term_cmp(x, y))
                                .find(|o| o.is_ne())
                                .unwrap_or(Ordering::Equal)
                        })
                        .then_with(|| ca.len().cmp(&cb.len()))
                },
            })
    }

    /// Children of `node` in canonical print order.
    pub(crate) fn ordered_children(
        &self,
        node: NodeId,
    ) -> Vec<NodeId> {
        let mut children = self.children(node).to_vec();
        if self.op(node) == core::MUL {
            // Coefficient first, then factors by base: `3*a^2*b`.
            children.sort_by(|&a, &b| {
                let (base_a, exp_a) = self.base_and_exponent(a);
                let (base_b, exp_b) = self.base_and_exponent(b);
                self.as_number(b)
                    .is_some()
                    .cmp(&self.as_number(a).is_some())
                    .then_with(|| self.term_cmp(base_a, base_b))
                    .then_with(|| exp_a.total_cmp(&exp_b))
                    .then_with(|| self.term_cmp(a, b))
            });
        } else if self
            .ops()
            .get(self.op(node))
            .flags
            .has(OpFlags::COMMUTATIVE)
        {
            children.sort_by(|&a, &b| self.term_cmp(a, b));
        }
        children
    }

    /// Splits `b^e` with a literal exponent into `(b, e)`; anything else is
    /// its own base with exponent one.
    pub(crate) fn base_and_exponent(
        &self,
        node: NodeId,
    ) -> (NodeId, f64) {
        match *self.children(node) {
            | [base, exp] if self.op(node) == core::POW => match self.as_number(exp) {
                | Some(n) => (base, n.to_f64()),
                | None => (node, 1.0),
            },
            | _ => (node, 1.0),
        }
    }

    /// Renders a concrete term as infix text.
    #[must_use]
    pub fn display(
        &self,
        node: NodeId,
    ) -> String {
        let mut out = String::new();
        self.write_term(node, PREC_EQ, &mut out);
        out
    }

    /// Splits a product into `(coefficient, other factors)` when it has a
    /// literal numeric factor.
    pub(crate) fn split_coefficient(
        &self,
        node: NodeId,
    ) -> Option<(Number, Vec<NodeId>)> {
        if self.op(node) != core::MUL {
            return None;
        }
        let children = self.ordered_children(node);
        let (first, rest) = children.split_first()?;
        let coeff = self.as_number(*first)?.clone();
        Some((coeff, rest.to_vec()))
    }

    fn write_term(
        &self,
        node: NodeId,
        min_prec: u8,
        out: &mut String,
    ) {
        let prec = self.precedence(node);
        let parens = prec < min_prec;
        if parens {
            out.push('(');
        }
        match self.op(node) {
            | core::LIT | core::SYM => self.write_leaf(node, out),
            | core::EQ => {
                if let [lhs, rhs] = *self.children(node) {
                    self.write_term(lhs, PREC_ADD, out);
                    out.push_str(" = ");
                    self.write_term(rhs, PREC_ADD, out);
                }
            },
            | core::ADD => self.write_sum(node, out),
            | core::MUL => self.write_product(node, out),
            | core::POW if prec == PREC_MUL => self.write_factors(&Number::from(1), &[node], out),
            | core::POW => {
                if let [base, exp] = *self.children(node) {
                    self.write_term(base, PREC_ATOM, out);
                    out.push('^');
                    self.write_term(exp, PREC_ATOM, out);
                }
            },
            | core::APPLY if self.children(node).first().is_some_and(|&f| self.as_symbol(f).is_some()) => {
                // An unknown function applied: `y(t)`, as it is parsed.
                let children = self.children(node);
                self.write_leaf(children[0], out);
                out.push('(');
                for (i, &child) in children[1..].iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    self.write_term(child, PREC_EQ, out);
                }
                out.push(')');
            },
            | op if self.ops().get(op).arity == Arity::Fixed(0) => {
                out.push_str(&self.ops().get(op).name);
            },
            | op => {
                out.push_str(&self.ops().get(op).name);
                out.push('(');
                for (i, &child) in self.ordered_children(node).iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    self.write_term(child, PREC_EQ, out);
                }
                out.push(')');
            },
        }
        if parens {
            out.push(')');
        }
    }

    pub(crate) fn precedence(
        &self,
        node: NodeId,
    ) -> u8 {
        match self.op(node) {
            | core::EQ => PREC_EQ,
            | core::ADD => PREC_ADD,
            | core::MUL => match self.split_coefficient(node) {
                | Some((c, _)) if c.is_negative() => PREC_ADD,
                | _ => PREC_MUL,
            },
            | core::POW => match self.children(node).get(1).and_then(|&e| self.as_number(e)) {
                // Printed as a quotient: `1/x^2`.
                | Some(e) if e.is_negative() => PREC_MUL,
                | _ => PREC_POW,
            },
            | core::LIT => match self.as_number(node) {
                | Some(n) if n.is_negative() => PREC_ADD,
                | Some(Number::Rat(_)) => PREC_MUL,
                | _ => PREC_ATOM,
            },
            | _ => PREC_ATOM,
        }
    }

    fn write_leaf(
        &self,
        node: NodeId,
        out: &mut String,
    ) {
        // Writing to a String cannot fail.
        let _ = match self.payload(node) {
            | Some(Payload::Num(n)) => write!(out, "{n}"),
            | Some(Payload::Sym(s)) => write!(out, "{}", self.interner().symbol_name(*s)),
            | Some(Payload::Bool(b)) => write!(out, "{b}"),
            | Some(Payload::Str(s)) => write!(out, "{s:?}"),
            | Some(Payload::Blob(b)) => write!(out, "{b:?}"),
            | None => write!(out, "?"),
        };
    }

    fn write_sum(
        &self,
        node: NodeId,
        out: &mut String,
    ) {
        let terms = self.ordered_terms(node);
        for (i, &term) in terms.iter().enumerate() {
            let negated = self.negated(term);
            match (i, &negated) {
                | (0, Some(_)) => out.push('-'),
                | (0, None) => {},
                | (_, Some(_)) => out.push_str(" - "),
                | (_, None) => out.push_str(" + "),
            }
            match negated {
                | Some(Negated::Number(n)) => {
                    let _ = write!(out, "{n}");
                },
                | Some(Negated::Product(coeff, factors)) => {
                    self.write_factors(&coeff, &factors, out);
                },
                | None => self.write_term(term, PREC_MUL, out),
            }
        }
    }

    /// The terms of a sum in canonical print order: descending degree,
    /// leading with a positive term when there is one.
    pub(crate) fn ordered_terms(
        &self,
        node: NodeId,
    ) -> Vec<NodeId> {
        // Descending degree, the constant term last: `x^2 + x + 1`.
        let mut terms = self.children(node).to_vec();
        terms.sort_by(|&a, &b| {
            self.degree(b)
                .total_cmp(&self.degree(a))
                .then_with(|| {
                    // Graded lexicographic: `a^3 + a^2*b + a*b^2 + b^3`.
                    let (ca, cb) = (self.core_factors(a), self.core_factors(b));
                    ca.iter()
                        .zip(&cb)
                        .map(|(&x, &y)| {
                            let (base_x, exp_x) = self.base_and_exponent(x);
                            let (base_y, exp_y) = self.base_and_exponent(y);
                            self.term_cmp(base_x, base_y)
                                .then_with(|| exp_y.total_cmp(&exp_x))
                        })
                        .find(|o| o.is_ne())
                        .unwrap_or_else(|| ca.len().cmp(&cb.len()))
                })
                .then_with(|| self.term_cmp(a, b))
        });
        // Lead with a positive term when there is one: `1 - x^2`.
        if terms.first().is_some_and(|&t| self.negated(t).is_some()) {
            if let Some(i) = terms.iter().position(|&t| self.negated(t).is_none()) {
                let positive = terms.remove(i);
                terms.insert(0, positive);
            }
        }
        terms
    }

    /// A rough polynomial degree used only to order the terms of a sum.
    fn degree(
        &self,
        node: NodeId,
    ) -> f64 {
        if self.as_number(node).is_some() {
            return 0.0;
        }
        match (self.op(node), self.children(node)) {
            | (core::MUL, factors) => factors.iter().map(|&f| self.degree(f)).sum(),
            | (core::ADD, terms) => terms.iter().map(|&t| self.degree(t)).fold(0.0, f64::max),
            | (core::POW, &[base, exp]) => match self.as_number(exp) {
                | Some(n) => self.degree(base) * n.to_f64(),
                | None => 1.0,
            },
            | _ => 1.0,
        }
    }

    /// The factors of a term without its numeric coefficient.
    fn core_factors(
        &self,
        node: NodeId,
    ) -> Vec<NodeId> {
        if self.op(node) == core::MUL {
            self.ordered_children(node)
                .into_iter()
                .filter(|&f| self.as_number(f).is_none())
                .collect()
        } else {
            vec![node]
        }
    }

    /// If `term` is a negative number or has a negative coefficient, returns
    /// its negation in a printable form.
    pub(crate) fn negated(
        &self,
        term: NodeId,
    ) -> Option<Negated> {
        if let Some(n) = self.as_number(term) {
            return n.is_negative().then(|| Negated::Number(n.neg()));
        }
        let (coeff, factors) = self.split_coefficient(term)?;
        coeff
            .is_negative()
            .then(|| Negated::Product(coeff.neg(), factors))
    }

    fn write_product(
        &self,
        node: NodeId,
        out: &mut String,
    ) {
        match self.split_coefficient(node) {
            | Some((coeff, factors)) if coeff.is_negative() => {
                out.push('-');
                self.write_factors(&coeff.neg(), &factors, out);
            },
            | Some((coeff, factors)) => self.write_factors(&coeff, &factors, out),
            | None => self.write_factors(&Number::from(1), &self.ordered_children(node), out),
        }
    }

    /// Writes `coeff * factors` as `numerator/denominator`, moving factors
    /// with a negative literal exponent below the fraction bar.
    fn write_factors(
        &self,
        coeff: &Number,
        factors: &[NodeId],
        out: &mut String,
    ) {
        let mut numer: Vec<NodeId> = Vec::new();
        let mut denom: Vec<(NodeId, Number)> = Vec::new();
        for &factor in factors {
            let inverse = match *self.children(factor) {
                | [base, exp] if self.op(factor) == core::POW => self
                    .as_number(exp)
                    .filter(|e| e.is_negative())
                    .map(|e| (base, e.neg())),
                | _ => None,
            };
            match inverse {
                | Some(pair) => denom.push(pair),
                | None => numer.push(factor),
            }
        }
        let show_coeff = !coeff.is_one() || numer.is_empty();
        if show_coeff {
            let _ = write!(out, "{coeff}");
        }
        let mut first = !show_coeff;
        for &factor in &numer {
            if !first {
                out.push('*');
            }
            first = false;
            self.write_term(factor, PREC_NEG, out);
        }
        if denom.is_empty() {
            return;
        }
        out.push('/');
        let group = denom.len() > 1;
        if group {
            out.push('(');
        }
        for (i, (base, exp)) in denom.iter().enumerate() {
            if i > 0 {
                out.push('*');
            }
            if exp.is_one() {
                self.write_term(*base, if group { PREC_NEG } else { PREC_POW }, out);
            } else {
                self.write_term(*base, PREC_ATOM, out);
                if exp.is_integer() {
                    let _ = write!(out, "^{exp}");
                } else {
                    let _ = write!(out, "^({exp})");
                }
            }
        }
        if group {
            out.push(')');
        }
    }
}

/// A negative term in printable form.
pub(crate) enum Negated {
    Number(Number),
    Product(Number, Vec<NodeId>),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;

    fn round_trip(src: &str) -> String {
        let mut g = Graph::new();
        assert!(
            g.ops_mut()
                .register(OpDescriptor::new("sin", Arity::Fixed(1)))
                .is_ok()
        );
        let node = g.parse(src).unwrap_or(NodeId::NONE);
        let text = g.display(node);
        assert_eq!(
            g.parse(&text).ok(),
            Some(node),
            "`{text}` must parse back to the same node"
        );
        text
    }

    #[test]
    fn sums_and_differences() {
        assert_eq!(round_trip("1 + x + x^2"), "x^2 + x + 1");
        assert_eq!(round_trip("a - b"), "a - b");
        assert_eq!(round_trip("-a - b"), "-a - b");
        assert_eq!(round_trip("a - 2*b*c"), "a - 2*b*c");
        assert_eq!(round_trip("1 - x^2"), "1 - x^2");
        assert_eq!(round_trip("-x - y"), "-x - y");
        assert_eq!(round_trip("x - 1"), "x - 1");
    }

    #[test]
    fn products_and_quotients() {
        assert_eq!(round_trip("y * 2 * x"), "2*x*y");
        assert_eq!(round_trip("a / b"), "a/b");
        assert_eq!(round_trip("a / (b * c)"), "a/(b*c)");
        assert_eq!(round_trip("1 / x^2"), "1/x^2");
        assert_eq!(round_trip("(a + b) * c"), "c*(a + b)");
        assert_eq!(round_trip("x * (1/2)"), "1/2*x");
    }

    #[test]
    fn powers_and_negation() {
        assert_eq!(round_trip("(a + b)^2"), "(a + b)^2");
        assert_eq!(round_trip("a^(b + 1)"), "a^(b + 1)");
        assert_eq!(round_trip("(-2)^x"), "(-2)^x");
        assert_eq!(round_trip("x^(-2)"), "1/x^2");
        assert_eq!(round_trip("(x + 1)^(-1)"), "1/(x + 1)");
        assert_eq!(round_trip("2^(-x)"), "2^(-x)");
        assert_eq!(round_trip("2/pi^(1/2)"), "2/pi^(1/2)");
        assert_eq!(round_trip("x^(-3/2)"), "1/x^(3/2)");
        assert_eq!(round_trip("x^(1/2)"), "x^(1/2)");
        assert_eq!(round_trip("-(a + b)"), "-(a + b)");
        assert_eq!(round_trip("a^b^c"), "a^(b^c)");
    }

    #[test]
    fn equations() {
        assert_eq!(round_trip("x^2 = 4"), "x^2 = 4");
        assert_eq!(round_trip("sin(a = b)"), "sin(a = b)");
        assert_eq!(round_trip("(a = b) + 1"), "(a = b) + 1");
        assert_eq!(round_trip("eq(x + 1, y)"), "x + 1 = y");
    }

    #[test]
    fn functions() {
        assert_eq!(round_trip("sin(x + 1)"), "sin(x + 1)");
        assert_eq!(round_trip("f(x, y)"), "f(x, y)");
    }

    #[test]
    fn order_is_independent_of_creation_order() {
        let mut g1 = Graph::new();
        let mut g2 = Graph::new();
        let a = g1.parse("c + b + a").unwrap_or(NodeId::NONE);
        let b = g2.parse("a + b + c").unwrap_or(NodeId::NONE);
        assert_eq!(g1.display(a), g2.display(b));
    }
}
