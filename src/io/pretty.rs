//! Two-dimensional Unicode layout of terms of a [`Graph`]: fractions with
//! bars, superscripts, radicals, stacked sums and tall brackets.

use super::markup::split_subscript;
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

const GREEK: &[(&str, char)] = &[
    ("alpha", 'α'),
    ("beta", 'β'),
    ("gamma", 'γ'),
    ("delta", 'δ'),
    ("epsilon", 'ε'),
    ("zeta", 'ζ'),
    ("eta", 'η'),
    ("theta", 'θ'),
    ("iota", 'ι'),
    ("kappa", 'κ'),
    ("lambda", 'λ'),
    ("mu", 'μ'),
    ("nu", 'ν'),
    ("xi", 'ξ'),
    ("pi", 'π'),
    ("rho", 'ρ'),
    ("sigma", 'σ'),
    ("tau", 'τ'),
    ("upsilon", 'υ'),
    ("phi", 'φ'),
    ("chi", 'χ'),
    ("psi", 'ψ'),
    ("omega", 'ω'),
    ("Gamma", 'Γ'),
    ("Delta", 'Δ'),
    ("Theta", 'Θ'),
    ("Lambda", 'Λ'),
    ("Xi", 'Ξ'),
    ("Pi", 'Π'),
    ("Sigma", 'Σ'),
    ("Phi", 'Φ'),
    ("Psi", 'Ψ'),
    ("Omega", 'Ω'),
];

/// Renders the term `node` as a multi-line Unicode drawing.
#[must_use]
pub fn to_pretty(
    graph: &Graph,
    node: NodeId,
) -> String {
    Printer { g: graph }
        .term(node, PREC_EQ)
        .lines
        .iter()
        .map(|l| l.trim_end())
        .collect::<Vec<_>>()
        .join("\n")
}

/// A rectangle of text with a baseline row.
#[derive(Clone)]
struct Block {
    lines: Vec<String>,
    width: usize,
    base: usize,
}

fn chars(s: &str) -> usize {
    s.chars().count()
}

impl Block {
    fn atom(s: &str) -> Self {
        Self {
            lines: vec![s.to_owned()],
            width: chars(s),
            base: 0,
        }
    }

    const fn height(&self) -> usize {
        self.lines.len()
    }

    const fn below(&self) -> usize {
        self.height().saturating_sub(self.base)
    }

    /// Side by side, baselines aligned.
    fn hcat(parts: &[Self]) -> Self {
        let base = parts.iter().map(|p| p.base).max().unwrap_or(0);
        let below = parts.iter().map(Self::below).max().unwrap_or(1);
        let height = base.saturating_add(below);
        let width = parts.iter().map(|p| p.width).sum();
        let mut lines = vec![String::new(); height];
        for part in parts {
            let top = base.saturating_sub(part.base);
            for (row, line) in lines.iter_mut().enumerate() {
                let cell = row
                    .checked_sub(top)
                    .and_then(|r| part.lines.get(r))
                    .map_or("", String::as_str);
                line.push_str(cell);
                line.extend(std::iter::repeat_n(' ', part.width.saturating_sub(chars(cell))));
            }
        }
        Self { lines, width, base }
    }

    fn text_row(
        s: &str,
        width: usize,
    ) -> String {
        let pad = width.saturating_sub(chars(s));
        let left = pad / 2;
        format!("{}{s}{}", " ".repeat(left), " ".repeat(pad.saturating_sub(left)))
    }

    /// Blocks stacked vertically and centred; the baseline is that of `base_part`.
    fn vcat(
        parts: &[Self],
        base_part: usize,
    ) -> Self {
        let width = parts.iter().map(|p| p.width).max().unwrap_or(0);
        let mut lines = Vec::new();
        let mut base = 0;
        for (i, part) in parts.iter().enumerate() {
            if i == base_part {
                base = lines.len().saturating_add(part.base);
            }
            for line in &part.lines {
                lines.push(Self::text_row(line, width));
            }
        }
        Self { lines, width, base }
    }

    fn frac(
        numer: &Self,
        denom: &Self,
    ) -> Self {
        let bar = Self::atom(&"─".repeat(numer.width.max(denom.width).saturating_add(2)));
        Self::vcat(&[numer.clone(), bar, denom.clone()], 1)
    }

    fn paren(inner: &Self) -> Self {
        Self::bracket(inner, ['(', ')'], ['⎛', '⎜', '⎝'], ['⎞', '⎟', '⎠'])
    }

    fn square(inner: &Self) -> Self {
        Self::bracket(inner, ['[', ']'], ['⎡', '⎢', '⎣'], ['⎤', '⎥', '⎦'])
    }

    fn bracket(
        inner: &Self,
        flat: [char; 2],
        left: [char; 3],
        right: [char; 3],
    ) -> Self {
        let h = inner.height();
        let side = |flat: char, c: [char; 3]| {
            let mut lines = Vec::with_capacity(h);
            for row in 0..h {
                lines.push(
                    match (h, row) {
                        | (1, _) => flat,
                        | (_, 0) => c[0],
                        | (_, r) if r + 1 == h => c[2],
                        | _ => c[1],
                    }
                    .to_string(),
                );
            }
            Self {
                lines,
                width: 1,
                base: inner.base,
            }
        };
        Self::hcat(&[side(flat[0], left), inner.clone(), side(flat[1], right)])
    }

    fn bars(inner: &Self) -> Self {
        let side = Self {
            lines: vec!["│".to_owned(); inner.height()],
            width: 1,
            base: inner.base,
        };
        Self::hcat(&[side.clone(), inner.clone(), side])
    }

    /// `exp` raised to the upper right of `base`.
    fn sup(
        base: &Self,
        exp: &Self,
    ) -> Self {
        let mut lines: Vec<String> = exp
            .lines
            .iter()
            .map(|l| format!("{}{l}", " ".repeat(base.width)))
            .collect();
        lines.extend(
            base.lines
                .iter()
                .map(|l| format!("{l}{}", " ".repeat(exp.width))),
        );
        Self {
            lines,
            width: base.width.saturating_add(exp.width),
            base: exp.height().saturating_add(base.base),
        }
    }

    fn padded(
        &self,
        width: usize,
    ) -> Self {
        Self {
            lines: self.lines.iter().map(|l| Self::text_row(l, width)).collect(),
            width,
            base: self.base,
        }
    }

    fn sqrt(
        inner: &Self,
        index: Option<&str>,
    ) -> Self {
        let root = index.map_or_else(|| "√".to_owned(), |n| format!("{}√", superscript(n).unwrap_or_else(|| format!("({n})"))));
        let top = format!("{}{}", " ".repeat(chars(&root)), "‾".repeat(inner.width));
        let body = Self::hcat(&[Self::atom(&root), inner.clone()]);
        Self {
            lines: std::iter::once(top).chain(body.lines).collect(),
            width: body.width,
            base: body.base.saturating_add(1),
        }
    }

    fn append(
        &self,
        text: &str,
    ) -> Self {
        Self::hcat(&[self.clone(), Self::atom(text)])
    }
}

fn superscript(s: &str) -> Option<String> {
    s.chars()
        .map(|c| match c {
            | '0' => Some('⁰'),
            | '1' => Some('¹'),
            | '2' => Some('²'),
            | '3' => Some('³'),
            | '4' => Some('⁴'),
            | '5' => Some('⁵'),
            | '6' => Some('⁶'),
            | '7' => Some('⁷'),
            | '8' => Some('⁸'),
            | '9' => Some('⁹'),
            | '-' => Some('⁻'),
            | 'n' => Some('ⁿ'),
            | _ => None,
        })
        .collect()
}

fn subscript(s: &str) -> Option<String> {
    s.chars()
        .map(|c| match c {
            | '0' => Some('₀'),
            | '1' => Some('₁'),
            | '2' => Some('₂'),
            | '3' => Some('₃'),
            | '4' => Some('₄'),
            | '5' => Some('₅'),
            | '6' => Some('₆'),
            | '7' => Some('₇'),
            | '8' => Some('₈'),
            | '9' => Some('₉'),
            | _ => None,
        })
        .collect()
}

struct Printer<'a> {
    g: &'a Graph,
}

impl Printer<'_> {
    fn name(
        &self,
        node: NodeId,
    ) -> &str {
        &self.g.ops().get(self.g.op(node)).name
    }

    fn symbol(name: &str) -> String {
        let (base, sub) = split_subscript(name);
        let head = GREEK
            .iter()
            .find(|(n, _)| *n == base)
            .map_or_else(|| base.to_owned(), |(_, c)| c.to_string());
        match sub {
            | Some(s) => format!("{head}{}", subscript(s).unwrap_or_else(|| format!("_{s}"))),
            | None => head,
        }
    }

    fn prec(
        &self,
        node: NodeId,
    ) -> u8 {
        match self.name(node) {
            | "exp" => PREC_POW,
            | "sqrt" | "abs" | "factorial" => PREC_ATOM,
            | "and" | "or" | "xor" | "implies" | "iff" | "lt" | "gt" | "le" | "ge" | "ne" => PREC_EQ,
            | _ => self.g.precedence(node),
        }
    }

    fn term(
        &self,
        node: NodeId,
        min_prec: u8,
    ) -> Block {
        let block = self.bare(node);
        if self.prec(node) < min_prec { Block::paren(&block) } else { block }
    }

    fn number(n: &Number) -> Block {
        let magnitude = if n.is_negative() { n.neg() } else { n.clone() };
        let text = magnitude.to_string();
        let body = match (&magnitude, text.split_once('/')) {
            | (Number::Rat(_), Some((p, q))) => Block::frac(&Block::atom(p), &Block::atom(q)),
            | (Number::Float(f), _) if f.is_infinite() => Block::atom("∞"),
            | _ => Block::atom(&text),
        };
        if n.is_negative() { Block::atom("-").pipe_hcat(&body) } else { body }
    }

    fn join(
        parts: &[Block],
        sep: &str,
    ) -> Block {
        let mut out: Vec<Block> = Vec::new();
        for (i, part) in parts.iter().enumerate() {
            if i > 0 {
                out.push(Block::atom(sep));
            }
            out.push(part.clone());
        }
        Block::hcat(&out)
    }

    fn args(
        &self,
        nodes: &[NodeId],
    ) -> Block {
        let parts: Vec<Block> = nodes.iter().map(|&c| self.term(c, PREC_EQ)).collect();
        Self::join(&parts, ", ")
    }

    fn call(
        &self,
        head: &str,
        nodes: &[NodeId],
    ) -> Block {
        Block::hcat(&[Block::atom(head), Block::paren(&self.args(nodes))])
    }

    #[allow(clippy::too_many_lines)]
    fn bare(
        &self,
        node: NodeId,
    ) -> Block {
        let children = self.g.children(node);
        match self.g.op(node) {
            | core::LIT => match self.g.payload(node) {
                | Some(Payload::Num(n)) => Self::number(n),
                | Some(Payload::Bool(b)) => Block::atom(if *b { "⊤" } else { "⊥" }),
                | Some(Payload::Str(s)) => Block::atom(&format!("\"{s}\"")),
                | _ => Block::atom("?"),
            },
            | core::SYM => self.g.as_symbol(node).map_or_else(
                || Block::atom("?"),
                |s| Block::atom(&Self::symbol(self.g.interner().symbol_name(s))),
            ),
            | core::EQ => match *children {
                | [l, r] => Block::hcat(&[self.term(l, PREC_ADD), Block::atom(" = "), self.term(r, PREC_ADD)]),
                | _ => self.call("eq", children),
            },
            | core::ADD => self.sum(node),
            | core::MUL => match self.g.split_coefficient(node) {
                | Some((c, f)) if c.is_negative() => Block::atom("-").pipe_hcat(&self.product(&c.neg(), &f)),
                | Some((c, f)) => self.product(&c, &f),
                | None => self.product(&Number::from(1), &self.g.ordered_children(node)),
            },
            | core::POW => {
                if self.prec(node) == PREC_MUL {
                    return self.product(&Number::from(1), &[node]);
                }
                match *children {
                    | [b, e] => self.pow_node(b, e),
                    | _ => self.call("pow", children),
                }
            },
            | core::LIST => self.list(children),
            | core::APPLY => match children.split_first() {
                | Some((&h, rest)) if self.g.as_symbol(h).is_some() => {
                    Block::hcat(&[self.term(h, PREC_ATOM), Block::paren(&self.args(rest))])
                },
                | _ => self.call("apply", children),
            },
            | _ => self.named(node),
        }
    }

    fn list(
        &self,
        children: &[NodeId],
    ) -> Block {
        let rows: Vec<&[NodeId]> = children
            .iter()
            .filter(|&&c| self.g.op(c) == core::LIST)
            .map(|&c| self.g.children(c))
            .collect();
        let width = rows.first().map(|r| r.len());
        if !children.is_empty() && rows.len() == children.len() && width.is_some_and(|w| w > 0) && rows.iter().all(|r| Some(r.len()) == width) {
            let cells: Vec<Vec<Block>> = rows
                .iter()
                .map(|r| r.iter().map(|&c| self.term(c, PREC_EQ)).collect())
                .collect();
            let widths: Vec<usize> = (0..width.unwrap_or(0))
                .map(|j| cells.iter().filter_map(|r| r.get(j)).map(|b| b.width).max().unwrap_or(0))
                .collect();
            let row_blocks: Vec<Block> = cells
                .iter()
                .map(|row| {
                    let mut parts = Vec::new();
                    for (j, cell) in row.iter().enumerate() {
                        if j > 0 {
                            parts.push(Block::atom("  "));
                        }
                        parts.push(cell.padded(widths.get(j).copied().unwrap_or(0)));
                    }
                    Block::hcat(&parts)
                })
                .collect();
            Block::paren(&Block::vcat(&row_blocks, row_blocks.len() / 2))
        } else {
            Block::square(&self.args(children))
        }
    }

    fn named(
        &self,
        node: NodeId,
    ) -> Block {
        let children = self.g.children(node);
        let name = self.name(node);
        if children.is_empty() {
            return Block::atom(match name {
                | "pi" => "π",
                | "oo" => "∞",
                | "I" => "i",
                | "E" => "e",
                | "true" => "⊤",
                | "false" => "⊥",
                | other => return Block::atom(&Self::symbol(other)),
            });
        }
        let infix = match name {
            | "and" => Some(" ∧ "),
            | "or" => Some(" ∨ "),
            | "xor" => Some(" ⊕ "),
            | "implies" => Some(" ⇒ "),
            | "iff" => Some(" ⇔ "),
            | "lt" => Some(" < "),
            | "gt" => Some(" > "),
            | "le" => Some(" ≤ "),
            | "ge" => Some(" ≥ "),
            | "ne" => Some(" ≠ "),
            | _ => None,
        };
        if let (Some(op), true) = (infix, children.len() > 1) {
            let parts: Vec<Block> = children.iter().map(|&c| self.term(c, PREC_ADD)).collect();
            return Self::join(&parts, op);
        }
        match (name, children) {
            | ("not", &[x]) => Block::atom("¬").pipe_hcat(&self.term(x, PREC_NEG)),
            | ("exp", &[x]) => Block::sup(&Block::atom("e"), &self.term(x, PREC_EQ)),
            | ("sqrt", &[x]) => Block::sqrt(&self.term(x, PREC_EQ), None),
            | ("abs", &[x]) => Block::bars(&self.term(x, PREC_EQ)),
            | ("conj", &[x]) => {
                let inner = self.term(x, PREC_ATOM);
                Block::vcat(&[Block::atom(&"‾".repeat(inner.width)), inner], 1)
            },
            | ("factorial", &[x]) => self.term(x, PREC_ATOM).append("!"),
            | ("binomial", &[n, k]) => Block::paren(&Block::vcat(&[self.term(n, PREC_EQ), self.term(k, PREC_EQ)], 0)),
            | ("diff" | "diffn", _) => self.derivative(node),
            | ("integral", &[f, x]) => Self::join(&[Block::atom("∫"), self.term(f, PREC_EQ), Block::atom(&format!("d{}", self.term(x, PREC_ATOM).lines.join("")))], " "),
            | ("integral" | "defint", &[f, x, a, b]) => {
                let sign = Block::vcat(&[self.term(b, PREC_EQ), Block::atom("∫"), self.term(a, PREC_EQ)], 1);
                Self::join(&[sign, self.term(f, PREC_EQ), Block::atom(&format!("d{}", self.term(x, PREC_ATOM).lines.join("")))], " ")
            },
            | ("sum" | "product", &[f, n, a, b]) => {
                let sign = if name == "sum" { "Σ" } else { "∏" };
                let lower = Block::hcat(&[self.term(n, PREC_EQ), Block::atom("="), self.term(a, PREC_EQ)]);
                let stack = Block::vcat(&[self.term(b, PREC_EQ), Block::atom(sign), lower], 1);
                Self::join(&[stack, self.term(f, PREC_MUL)], " ")
            },
            | ("limit", &[f, x, a, ..]) => {
                let under = Block::hcat(&[self.term(x, PREC_EQ), Block::atom(" → "), self.term(a, PREC_EQ)]);
                let stack = Block::vcat(&[Block::atom("lim"), under], 0);
                Self::join(&[stack, self.term(f, PREC_MUL)], " ")
            },
            | _ => self.call(name, children),
        }
    }

    fn derivative(
        &self,
        node: NodeId,
    ) -> Block {
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
            return self.call(self.name(node), self.g.children(node));
        }
        let partial = vars.len() > 1
            || (self.g.op(cur) == core::APPLY && self.g.children(cur).len() > 2);
        let d = if partial { "∂" } else { "d" };
        let total: u32 = vars.iter().map(|&(_, o)| o).sum();
        let lead = if total == 1 { d.to_owned() } else { format!("{d}{}", superscript(&total.to_string()).unwrap_or_default()) };
        let mut denom: Vec<Block> = Vec::new();
        for (i, &(v, o)) in vars.iter().enumerate() {
            if i > 0 {
                denom.push(Block::atom(" "));
            }
            let name = self.term(v, PREC_ATOM);
            let order = if o == 1 { String::new() } else { superscript(&o.to_string()).unwrap_or_default() };
            denom.push(Block::hcat(&[Block::atom(d), name]).append(&order));
        }
        let frac = Block::frac(&Block::atom(&lead), &Block::hcat(&denom));
        Self::join(&[frac, self.term(cur, PREC_MUL)], " ")
    }

    fn sum(
        &self,
        node: NodeId,
    ) -> Block {
        let mut parts: Vec<Block> = Vec::new();
        for (i, &term) in self.g.ordered_terms(node).iter().enumerate() {
            let negated = self.g.negated(term);
            parts.push(Block::atom(match (i, negated.is_some()) {
                | (0, true) => "-",
                | (0, false) => "",
                | (_, true) => " - ",
                | (_, false) => " + ",
            }));
            parts.push(match negated {
                | Some(Negated::Number(n)) => Self::number(&n),
                | Some(Negated::Product(c, f)) => self.product(&c, &f),
                | None => self.term(term, PREC_MUL),
            });
        }
        Block::hcat(&parts)
    }

    fn pow_number(
        &self,
        base: NodeId,
        exp: &Number,
    ) -> Block {
        let text = exp.to_string();
        if let (Number::Rat(_), Some(("1", q))) = (exp, text.split_once('/')) {
            return Block::sqrt(&self.term(base, PREC_EQ), (q != "2").then_some(q));
        }
        let b = self.term(base, PREC_ATOM);
        match (exp.is_negative(), exp) {
            | (false, Number::Int(_)) => match superscript(&text) {
                | Some(s) => b.append(&s),
                | None => Block::sup(&b, &Self::number(exp)),
            },
            | _ => Block::sup(&b, &Self::number(exp)),
        }
    }

    fn pow_node(
        &self,
        base: NodeId,
        exp: NodeId,
    ) -> Block {
        if let Some(n) = self.g.as_number(exp) {
            return self.pow_number(base, n);
        }
        if let [index, minus_one] = *self.g.children(exp) {
            if self.g.op(exp) == core::POW
                && self.g.as_number(minus_one).is_some_and(|m| *m == Number::from(-1))
            {
                let idx = self.term(index, PREC_EQ).lines.join("");
                return Block::sqrt(&self.term(base, PREC_EQ), Some(&idx));
            }
        }
        Block::sup(&self.term(base, PREC_ATOM), &self.term(exp, PREC_EQ))
    }

    fn product(
        &self,
        coeff: &Number,
        factors: &[NodeId],
    ) -> Block {
        let mut numer: Vec<NodeId> = Vec::new();
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
            | _ => (text.clone(), None),
        };
        let mut denom: Vec<Block> = coeff_denom.iter().map(|q| Block::atom(q)).collect();
        let group = denom.len() + inverses.len() > 1;
        for (base, exp) in &inverses {
            denom.push(if exp.is_one() {
                self.term(*base, if group { PREC_NEG } else { PREC_EQ })
            } else {
                self.pow_number(*base, exp)
            });
        }
        let mut parts: Vec<(Block, bool)> = Vec::new();
        if coeff_numer != "1" || numer.is_empty() {
            parts.push((Block::atom(&coeff_numer), true));
        }
        parts.extend(numer.iter().map(|&f| (self.term(f, PREC_NEG), self.g.as_number(f).is_some())));
        let lead_coeff = parts.first().is_some_and(|p| p.1) && coeff_numer != "1";
        let join = |parts: &[(Block, bool)], lead_coeff: bool| {
            let mut out: Vec<Block> = Vec::new();
            for (i, (block, numeric)) in parts.iter().enumerate() {
                if i > 0 && (*numeric || !(lead_coeff && i == 1)) {
                    out.push(Block::atom("·"));
                }
                out.push(block.clone());
            }
            Block::hcat(&out)
        };
        let top = join(&parts, lead_coeff);
        if denom.is_empty() {
            top
        } else {
            let den: Vec<(Block, bool)> = denom.into_iter().map(|b| (b, false)).collect();
            Block::frac(&top, &join(&den, false))
        }
    }
}

impl Block {
    /// `self` followed by `next`, baselines aligned.
    fn pipe_hcat(
        &self,
        next: &Self,
    ) -> Self {
        Self::hcat(&[self.clone(), next.clone()])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pretty(src: &str) -> String {
        let mut g = Graph::new();
        let _e = crate::graph::Engine::install(&mut g, &crate::rules::standard()).unwrap();
        let n = g.parse(src).unwrap();
        to_pretty(&g, n)
    }

    #[test]
    fn fraction_has_bar() {
        let out = pretty("a/b");
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].contains('a'));
        assert!(lines[1].contains('─'));
        assert!(lines[2].contains('b'));
    }

    #[test]
    fn superscripts_and_subscripts() {
        assert_eq!(pretty("sin(x)^2"), "sin(x)²");
        assert_eq!(pretty("x_1^2"), "x₁²");
    }

    #[test]
    fn root_and_greek() {
        assert!(pretty("sqrt(x+1)").contains('√'));
        assert_eq!(pretty("alpha*beta"), "α·β");
    }

    #[test]
    fn lists_and_matrices() {
        assert_eq!(pretty("list(1, 2, 3)"), "[1, 2, 3]");
        let m = pretty("list(list(1, 2), list(3, 4))");
        assert_eq!(m.lines().count(), 2);
        assert!(m.contains('⎛') && m.contains('⎝'));
    }

    #[test]
    fn big_operators() {
        let s = pretty("sum(n^2, n, 1, 10)");
        assert!(s.contains('Σ') && s.contains("n=1") && s.contains("10"));
        let i = pretty("defint(f(x), x, 0, 1)");
        assert!(i.contains('∫') && i.contains("dx"));
        let l = pretty("limit(sin(x)/x, x, 0)");
        assert!(l.contains("lim") && l.contains('→'));
    }

    #[test]
    fn derivatives() {
        assert!(pretty("diff(diff(u(x, t), x), x)").contains('∂'));
        assert!(pretty("diff(f(x), x)").contains("dx"));
    }
}
