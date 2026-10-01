//! The term walker shared by the LaTeX and Typst printers.
//!
//! The walker owns everything structural: precedence, re-sugaring of
//! `add` / `mul` / `pow` into differences, quotients and roots, and the
//! recognition of derivative chains, binders and matrices. A [`Dialect`]
//! only supplies the concrete syntax of each construct.

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::Payload;
use crate::graph::op::core;
use crate::graph::print::Negated;
use crate::graph::print::PREC_ADD;
use crate::graph::print::PREC_ATOM;
use crate::graph::print::PREC_EQ;
use crate::graph::print::PREC_MUL;
use crate::graph::print::PREC_NEG;
use crate::graph::print::PREC_POW;

/// Greek letter names and their LaTeX commands. Capitals that have no
/// distinct LaTeX glyph map to the Latin letter.
pub(crate) const GREEK: &[(&str, &str)] = &[
    ("alpha", r"\alpha"),
    ("beta", r"\beta"),
    ("gamma", r"\gamma"),
    ("delta", r"\delta"),
    ("epsilon", r"\epsilon"),
    ("zeta", r"\zeta"),
    ("eta", r"\eta"),
    ("theta", r"\theta"),
    ("iota", r"\iota"),
    ("kappa", r"\kappa"),
    ("lambda", r"\lambda"),
    ("mu", r"\mu"),
    ("nu", r"\nu"),
    ("xi", r"\xi"),
    ("pi", r"\pi"),
    ("rho", r"\rho"),
    ("sigma", r"\sigma"),
    ("tau", r"\tau"),
    ("upsilon", r"\upsilon"),
    ("phi", r"\phi"),
    ("chi", r"\chi"),
    ("psi", r"\psi"),
    ("omega", r"\omega"),
    ("Alpha", "A"),
    ("Beta", "B"),
    ("Gamma", r"\Gamma"),
    ("Delta", r"\Delta"),
    ("Epsilon", "E"),
    ("Zeta", "Z"),
    ("Eta", "H"),
    ("Theta", r"\Theta"),
    ("Iota", "I"),
    ("Kappa", "K"),
    ("Lambda", r"\Lambda"),
    ("Mu", "M"),
    ("Nu", "N"),
    ("Xi", r"\Xi"),
    ("Pi", r"\Pi"),
    ("Rho", "P"),
    ("Sigma", r"\Sigma"),
    ("Tau", "T"),
    ("Upsilon", r"\Upsilon"),
    ("Phi", r"\Phi"),
    ("Chi", "X"),
    ("Psi", r"\Psi"),
    ("Omega", r"\Omega"),
];

/// A `sum` or `product` binder.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum BigOp {
    Sum,
    Product,
}

/// A chain of derivatives, already rendered.
pub(crate) struct Deriv {
    /// Partial (`∂`) rather than total (`d`) derivatives.
    pub partial: bool,
    /// Total order.
    pub total: u32,
    /// Variables with their orders, outermost first.
    pub vars: Vec<(String, u32)>,
    /// The differentiated term.
    pub body: String,
    /// Whether `body` is a bare function name that may sit in the numerator.
    pub inline: bool,
}

/// The concrete syntax of a target markup language.
pub(crate) trait Dialect {
    fn symbol(
        &self,
        name: &str,
    ) -> String;
    fn text(
        &self,
        s: &str,
    ) -> String;
    fn constant(
        &self,
        name: &str,
    ) -> Option<&'static str>;
    fn frac(
        &self,
        numer: &str,
        denom: &str,
    ) -> String;
    fn sqrt(
        &self,
        x: &str,
        index: Option<&str>,
    ) -> String;
    fn pow(
        &self,
        base: &str,
        exp: &str,
    ) -> String;
    fn paren(
        &self,
        s: &str,
    ) -> String;
    fn mul_sep(
        &self,
        prev: &str,
        next: &str,
    ) -> &'static str;
    fn function(
        &self,
        name: &str,
        args: &[String],
    ) -> String;
    fn abs(
        &self,
        x: &str,
    ) -> String;
    fn factorial(
        &self,
        x: &str,
    ) -> String;
    fn binom(
        &self,
        n: &str,
        k: &str,
    ) -> String;
    fn conj(
        &self,
        x: &str,
    ) -> String;
    fn bessel(
        &self,
        kind: char,
        order: &str,
        x: &str,
    ) -> String;
    fn deriv(
        &self,
        d: &Deriv,
    ) -> String;
    fn integral(
        &self,
        f: &str,
        var: &str,
        bounds: Option<(&str, &str)>,
    ) -> String;
    fn big_op(
        &self,
        kind: BigOp,
        f: &str,
        var: &str,
        lo: &str,
        hi: &str,
    ) -> String;
    fn limit(
        &self,
        f: &str,
        var: &str,
        to: &str,
    ) -> String;
    fn matrix(
        &self,
        rows: &[Vec<String>],
    ) -> String;
    fn list(
        &self,
        items: &[String],
    ) -> String;
    fn infix(
        &self,
        name: &str,
    ) -> Option<&'static str>;
    fn not(
        &self,
        x: &str,
    ) -> String;
}

/// Splits `name` into a base and an optional subscript at the first `_`.
pub(crate) fn split_subscript(name: &str) -> (&str, Option<&str>) {
    match name.split_once('_') {
        | Some((base, sub)) if !base.is_empty() && !sub.is_empty() => (base, Some(sub)),
        | _ => (name, None),
    }
}

/// Renders `node` in the syntax of `dialect`.
pub(crate) fn render<D: Dialect>(
    graph: &Graph,
    dialect: &D,
    node: NodeId,
) -> String {
    Walker { g: graph, d: dialect }.term(node, PREC_EQ)
}

struct Walker<'a, D> {
    g: &'a Graph,
    d: &'a D,
}

impl<D: Dialect> Walker<'_, D> {
    fn name(
        &self,
        node: NodeId,
    ) -> &str {
        &self.g.ops().get(self.g.op(node)).name
    }

    fn symbol_name(
        &self,
        node: NodeId,
    ) -> Option<&str> {
        self.g
            .as_symbol(node)
            .map(|s| self.g.interner().symbol_name(s))
    }

    fn number(
        &self,
        n: &Number,
    ) -> String {
        let magnitude = if n.is_negative() { n.neg() } else { n.clone() };
        let body = match &magnitude {
            | Number::Rat(_) => {
                let text = magnitude.to_string();
                text.split_once('/').map_or_else(
                    || text.clone(),
                    |(p, q)| self.d.frac(p, q),
                )
            },
            | Number::Float(f) if f.is_infinite() => {
                self.d.constant("oo").unwrap_or("inf").to_owned()
            },
            | _ => magnitude.to_string(),
        };
        if n.is_negative() { format!("-{body}") } else { body }
    }

    /// Prec of the node as printed by this walker.
    fn prec(
        &self,
        node: NodeId,
    ) -> u8 {
        let op = self.g.op(node);
        if op == core::EQ || self.d.infix(self.name(node)).is_some() {
            return PREC_EQ;
        }
        match self.name(node) {
            | "not" => PREC_NEG,
            | "exp" => PREC_POW,
            | "sqrt" => PREC_ATOM,
            | "diff" | "diffn" | "integral" | "defint" | "sum" | "product" | "limit" => PREC_MUL,
            | "factorial" => PREC_POW,
            | _ => self.g.precedence(node),
        }
    }

    fn term(
        &self,
        node: NodeId,
        min_prec: u8,
    ) -> String {
        let text = self.term_bare(node);
        if self.prec(node) < min_prec { self.d.paren(&text) } else { text }
    }

    fn args(
        &self,
        nodes: &[NodeId],
    ) -> Vec<String> {
        nodes.iter().map(|&c| self.term(c, PREC_EQ)).collect()
    }

    #[allow(clippy::too_many_lines)]
    fn term_bare(
        &self,
        node: NodeId,
    ) -> String {
        let children = self.g.children(node);
        match self.g.op(node) {
            | core::LIT => match self.g.payload(node) {
                | Some(Payload::Num(n)) => self.number(n),
                | Some(Payload::Bool(b)) => self
                    .d
                    .constant(if *b { "true" } else { "false" })
                    .map_or_else(|| b.to_string(), str::to_owned),
                | Some(Payload::Str(s)) => self.d.text(s),
                | _ => "?".to_owned(),
            },
            | core::SYM => self
                .symbol_name(node)
                .map_or_else(|| "?".to_owned(), |n| self.d.symbol(n)),
            | core::EQ => match *children {
                | [lhs, rhs] => format!(
                    "{} = {}",
                    self.term(lhs, PREC_ADD),
                    self.term(rhs, PREC_ADD)
                ),
                | _ => self.call(node),
            },
            | core::ADD => self.sum(node),
            | core::MUL => match self.g.split_coefficient(node) {
                | Some((coeff, factors)) if coeff.is_negative() => {
                    format!("-{}", self.product(&coeff.neg(), &factors))
                },
                | Some((coeff, factors)) => self.product(&coeff, &factors),
                | None => self.product(&Number::from(1), &self.g.ordered_children(node)),
            },
            | core::POW => {
                if self.prec(node) == PREC_MUL {
                    return self.product(&Number::from(1), &[node]);
                }
                match *children {
                    | [base, exp] => self.pow_node(base, exp),
                    | _ => self.call(node),
                }
            },
            | core::LIST => self.list(children),
            | core::APPLY => match children.split_first() {
                | Some((&head, rest)) if self.symbol_name(head).is_some() => {
                    format!("{}({})", self.term(head, PREC_ATOM), self.args(rest).join(", "))
                },
                | _ => self.call(node),
            },
            | _ => self.named(node),
        }
    }

    fn list(
        &self,
        children: &[NodeId],
    ) -> String {
        let rows: Vec<&[NodeId]> = children
            .iter()
            .filter(|&&c| self.g.op(c) == core::LIST)
            .map(|&c| self.g.children(c))
            .collect();
        let is_matrix = !children.is_empty()
            && rows.len() == children.len()
            && rows.first().is_some_and(|r| !r.is_empty())
            && rows.iter().all(|r| Some(r.len()) == rows.first().map(|f| f.len()));
        if is_matrix {
            let cells: Vec<Vec<String>> = rows.iter().map(|r| self.args(r)).collect();
            self.d.matrix(&cells)
        } else {
            self.d.list(&self.args(children))
        }
    }

    /// Known operators with a dedicated notation.
    fn named(
        &self,
        node: NodeId,
    ) -> String {
        let children = self.g.children(node);
        let name = self.name(node);
        if children.is_empty() {
            return self
                .d
                .constant(name)
                .map_or_else(|| self.d.symbol(name), str::to_owned);
        }
        if let Some(op) = self.d.infix(name) {
            if children.len() > 1 {
                let parts: Vec<String> = children.iter().map(|&c| self.term(c, PREC_ADD)).collect();
                return parts.join(op);
            }
        }
        match (name, children) {
            | ("not", &[x]) => self.d.not(&self.term(x, PREC_NEG)),
            | ("exp", &[x]) => self.d.pow(&self.d.constant("E").unwrap_or("e").to_owned(), &self.term(x, PREC_EQ)),
            | ("sqrt", &[x]) => self.d.sqrt(&self.term(x, PREC_EQ), None),
            | ("abs", &[x]) => self.d.abs(&self.term(x, PREC_EQ)),
            | ("conj", &[x]) => self.d.conj(&self.term(x, PREC_EQ)),
            | ("factorial", &[x]) => self.d.factorial(&self.term(x, PREC_ATOM)),
            | ("binomial", &[n, k]) => self.d.binom(&self.term(n, PREC_EQ), &self.term(k, PREC_EQ)),
            | ("besselj" | "bessely" | "besseli" | "besselk", &[order, x]) => {
                let kind = match name {
                    | "besselj" => 'J',
                    | "bessely" => 'Y',
                    | "besseli" => 'I',
                    | _ => 'K',
                };
                self.d.bessel(kind, &self.term(order, PREC_EQ), &self.term(x, PREC_EQ))
            },
            | ("diff" | "diffn", _) => self.derivative(node),
            | ("integral", &[f, x]) => self.d.integral(&self.term(f, PREC_EQ), &self.term(x, PREC_EQ), None),
            | ("integral" | "defint", &[f, x, a, b]) => self.d.integral(
                &self.term(f, PREC_EQ),
                &self.term(x, PREC_EQ),
                Some((&self.term(a, PREC_EQ), &self.term(b, PREC_EQ))),
            ),
            | ("sum" | "product", &[f, n, a, b]) => self.d.big_op(
                if name == "sum" { BigOp::Sum } else { BigOp::Product },
                &self.term(f, PREC_MUL),
                &self.term(n, PREC_EQ),
                &self.term(a, PREC_EQ),
                &self.term(b, PREC_EQ),
            ),
            | ("limit", &[f, x, a, ..]) => self.d.limit(
                &self.term(f, PREC_MUL),
                &self.term(x, PREC_EQ),
                &self.term(a, PREC_EQ),
            ),
            | _ => self.call(node),
        }
    }

    fn call(
        &self,
        node: NodeId,
    ) -> String {
        self.d
            .function(self.name(node), &self.args(self.g.children(node)))
    }

    fn derivative(
        &self,
        node: NodeId,
    ) -> String {
        let mut vars: Vec<(NodeId, u32)> = Vec::new();
        let mut cur = node;
        loop {
            let (var, order, inner) = match (self.name(cur), self.g.children(cur)) {
                | ("diff", &[f, x]) => (x, 1, f),
                | ("diffn", &[f, x, n]) => match self.g.as_number(n).and_then(Number::to_i64) {
                    | Some(k) if k >= 1 => (x, u32::try_from(k).unwrap_or(1), f),
                    | _ => break,
                },
                | _ => break,
            };
            match vars.iter_mut().find(|(v, _)| *v == var) {
                | Some((_, o)) => *o = o.saturating_add(order),
                | None => vars.push((var, order)),
            }
            cur = inner;
        }
        if cur == node {
            return self.call(node);
        }
        let head = match self.g.children(cur).split_first() {
            | Some((&h, rest)) if self.g.op(cur) == core::APPLY => {
                self.symbol_name(h).map(|_| (h, rest.len()))
            },
            | _ => None,
        };
        let partial = vars.len() > 1 || head.is_some_and(|(_, n)| n > 1);
        let inline = partial && head.is_some();
        let body = match head {
            | Some((h, _)) => self.term(h, PREC_ATOM),
            | None => self.term(cur, PREC_MUL),
        };
        self.d.deriv(&Deriv {
            partial,
            total: vars.iter().map(|&(_, o)| o).sum(),
            vars: vars
                .iter()
                .map(|&(v, o)| (self.term(v, PREC_ATOM), o))
                .collect(),
            body,
            inline,
        })
    }

    fn sum(
        &self,
        node: NodeId,
    ) -> String {
        let mut out = String::new();
        for (i, &term) in self.g.ordered_terms(node).iter().enumerate() {
            let negated = self.g.negated(term);
            match (i, negated.is_some()) {
                | (0, true) => out.push('-'),
                | (0, false) => {},
                | (_, true) => out.push_str(" - "),
                | (_, false) => out.push_str(" + "),
            }
            match negated {
                | Some(Negated::Number(n)) => out.push_str(&self.number(&n)),
                | Some(Negated::Product(coeff, factors)) => {
                    out.push_str(&self.product(&coeff, &factors));
                },
                | None => out.push_str(&self.term(term, PREC_MUL)),
            }
        }
        out
    }

    /// `base^exp` for a literal exponent.
    fn pow_number(
        &self,
        base: NodeId,
        exp: &Number,
    ) -> String {
        if let Number::Rat(_) = exp {
            let text = exp.to_string();
            if let Some(("1", q)) = text.split_once('/') {
                let inner = self.term(base, PREC_EQ);
                return self.d.sqrt(&inner, (q != "2").then_some(q));
            }
        }
        self.d
            .pow(&self.term(base, PREC_ATOM), &self.number(exp))
    }

    fn pow_node(
        &self,
        base: NodeId,
        exp: NodeId,
    ) -> String {
        if let Some(n) = self.g.as_number(exp) {
            return self.pow_number(base, n);
        }
        // x^(1/n): the exponent is `n^(-1)`.
        if let [index, minus_one] = *self.g.children(exp) {
            if self.g.op(exp) == core::POW
                && self.g.as_number(minus_one).is_some_and(|m| *m == Number::from(-1))
            {
                return self
                    .d
                    .sqrt(&self.term(base, PREC_EQ), Some(&self.term(index, PREC_EQ)));
            }
        }
        self.d
            .pow(&self.term(base, PREC_ATOM), &self.term(exp, PREC_EQ))
    }

    /// `coeff * factors` with `coeff > 0`, as a quotient when some factor
    /// has a negative literal exponent or the coefficient is a fraction.
    fn product(
        &self,
        coeff: &Number,
        factors: &[NodeId],
    ) -> String {
        let mut numer: Vec<NodeId> = Vec::new();
        let mut denom: Vec<String> = Vec::new();
        let mut inverses: Vec<(NodeId, Number)> = Vec::new();
        for &factor in factors {
            let inverse = match *self.g.children(factor) {
                | [base, exp] if self.g.op(factor) == core::POW => self
                    .g
                    .as_number(exp)
                    .filter(|e| e.is_negative())
                    .map(|e| (base, e.neg())),
                | _ => None,
            };
            match inverse {
                | Some(pair) => inverses.push(pair),
                | None => numer.push(factor),
            }
        }
        let text = coeff.to_string();
        let (coeff_numer, coeff_denom) = match (coeff, text.split_once('/')) {
            | (Number::Rat(_), Some((p, q))) => (p.to_owned(), Some(q.to_owned())),
            | _ => (self.number(coeff), None),
        };
        denom.extend(coeff_denom);
        let group = denom.len() + inverses.len() > 1;
        for (base, exp) in &inverses {
            denom.push(if exp.is_one() {
                self.term(*base, if group { PREC_NEG } else { PREC_EQ })
            } else {
                self.pow_number(*base, exp)
            });
        }
        let mut parts: Vec<String> = Vec::new();
        let show_coeff = coeff_numer != "1" || numer.is_empty();
        if show_coeff {
            parts.push(coeff_numer);
        }
        // A lone numerator sits inside the fraction bar and needs no parentheses.
        let numer_prec = if !denom.is_empty() && !show_coeff && numer.len() == 1 { PREC_EQ } else { PREC_NEG };
        parts.extend(numer.iter().map(|&f| self.term(f, numer_prec)));
        let join = |parts: &[String]| {
            let mut out = String::new();
            for (i, part) in parts.iter().enumerate() {
                if let Some(prev) = i.checked_sub(1).and_then(|j| parts.get(j)) {
                    out.push_str(self.d.mul_sep(prev, part));
                }
                out.push_str(part);
            }
            out
        };
        let top = join(&parts);
        if denom.is_empty() {
            top
        } else {
            self.d.frac(&top, &join(&denom))
        }
    }
}
