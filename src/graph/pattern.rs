//! Patterns: the declarative half of the rule language.
//!
//! A [`Pat`] is a term with holes. It can be parsed from text
//! (`"sin(?x)^2 + cos(?x)^2"`), matched against an e-node of a [`Graph`]
//! (e-matching, modulo the equalities the graph knows) and instantiated
//! back into the graph under a substitution.
//!
//! # Matching modulo associativity and commutativity
//!
//! Commutative operators are matched under every assignment of pattern
//! children to node children. When the *root* of a pattern is an associative
//! commutative operator and the node has more children than the pattern,
//! the surplus children are collected as the **rest** of the match:
//! `?a + ?a` matches `x + y + x` with `?a = x`, rest `= [y]`. Instantiating
//! the right-hand side re-attaches the rest, which is sound because
//! `f(p…, r…) = f(f(p…), r…)` for associative `f`.
//!
//! Below the root there is nowhere to re-attach a rest, so a *nested*
//! associative-commutative pattern must account for every operand. To make
//! that practical, when the node has more operands than the pattern and the
//! pattern's **last** operand is a bare variable, that variable absorbs all
//! the operands the others did not take: in
//! `integral(?c * ?f, ?x) => ?c * integral(?f, ?x)` the pattern `?c * ?f`
//! matches `2 * x * sin(x)` with `?c = 2`, `?f = x * sin(x)` (and with the
//! other choices of `?c`).

use std::fmt;

use super::id::NodeId;
use super::id::OpId;
use super::number::Number;
use super::op::Arity;
use super::op::OpFlags;
use super::op::core;
use super::store::Graph;

/// A term with holes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Pat {
    /// A hole, identified by its index in the substitution.
    Var(u32),
    /// Matches any class known to equal this number.
    Num(Number),
    /// Matches the class of this concrete node (symbols, other literals).
    Exact(NodeId),
    /// An operator applied to sub-patterns.
    Node(OpId, Vec<Self>),
}

/// A successful match of a pattern against an e-node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Match {
    /// The matched e-node.
    pub root: NodeId,
    /// Binding of each pattern variable to a concrete node of the matched
    /// class. Unbound variables (those not occurring in the pattern) are
    /// [`NodeId::NONE`].
    pub subst: Vec<NodeId>,
    /// Children of the root left over by an associative-commutative match.
    pub rest: Vec<NodeId>,
    /// Variables bound to *several* operands of a nested associative
    /// operator: `(variable, operator, operands)`. Their entry in `subst`
    /// is [`NodeId::NONE`] until [`Match::materialize`] builds the node.
    pub multi: Vec<(u32, OpId, Vec<NodeId>)>,
}

impl Match {
    /// The complete substitution, building `operator(operands...)` for
    /// every variable that absorbed several operands.
    pub fn materialize(
        &self,
        graph: &mut Graph,
    ) -> Vec<NodeId> {
        let mut subst = self.subst.clone();
        for (var, op, operands) in &self.multi {
            if let Some(slot) = subst.get_mut(*var as usize) {
                *slot = graph.node(*op, operands);
            }
        }
        subst
    }
}

/// Error produced while parsing pattern or rule text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseError {
    /// What went wrong.
    pub message: String,
    /// Byte offset into the source text.
    pub offset: usize,
}

impl fmt::Display for ParseError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "{} at offset {}", self.message, self.offset)
    }
}

impl std::error::Error for ParseError {}

/// Names of the variables of a parsed pattern, in index order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct VarNames(Vec<String>);

impl VarNames {
    /// Index of `name`, allocating a new variable if needed.
    pub fn index(
        &mut self,
        name: &str,
    ) -> u32 {
        if let Some(i) = self.0.iter().position(|n| n == name) {
            return u32::try_from(i).unwrap_or(u32::MAX);
        }
        self.0.push(name.to_owned());
        u32::try_from(self.0.len().saturating_sub(1)).unwrap_or(u32::MAX)
    }

    /// Index of `name` if it has been seen.
    #[must_use]
    pub fn get(
        &self,
        name: &str,
    ) -> Option<u32> {
        self.0
            .iter()
            .position(|n| n == name)
            .and_then(|i| u32::try_from(i).ok())
    }

    /// Number of variables.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether no variable has been seen.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

// ----------------------------------------------------------------------
// Parsing
// ----------------------------------------------------------------------

/// Recursive-descent parser for the infix term syntax.
///
/// ```text
/// expr  := sum ('=' sum)?
/// sum   := term (('+' | '-') term)*
/// term  := unary (('*' | '/') unary)*
/// unary := '-' unary | power
/// power := atom ('^' unary)?
/// atom  := NUMBER | '?' IDENT | IDENT '(' expr,* ')' | IDENT | '(' expr ')'
/// ```
///
/// Subtraction, negation and division are sugar for the canonical
/// `add` / `mul` / `pow` forms, so no rule ever has to mention them.
pub(crate) struct Parser<'a> {
    src: &'a str,
    pos: usize,
    graph: &'a mut Graph,
    vars: &'a mut VarNames,
    /// Treat `f(x)` with an unregistered `f` as `apply(f, x)` instead of an
    /// error. Wanted for user input, not for rule definitions.
    lenient: bool,
}

impl<'a> Parser<'a> {
    pub(crate) const fn new(
        src: &'a str,
        graph: &'a mut Graph,
        vars: &'a mut VarNames,
        lenient: bool,
    ) -> Self {
        Self {
            src,
            pos: 0,
            graph,
            vars,
            lenient,
        }
    }

    pub(crate) const fn offset(&self) -> usize {
        self.pos
    }

    fn err<T>(
        &self,
        message: impl Into<String>,
    ) -> Result<T, ParseError> {
        Err(ParseError {
            message: message.into(),
            offset: self.pos,
        })
    }

    fn rest(&self) -> &'a str {
        self.src.get(self.pos..).unwrap_or("")
    }

    pub(crate) fn skip_ws(&mut self) {
        let trimmed = self.rest().trim_start();
        self.pos = self.src.len().saturating_sub(trimmed.len());
    }

    fn peek(&mut self) -> Option<char> {
        self.skip_ws();
        self.rest().chars().next()
    }

    /// Consumes `token` if it is next.
    pub(crate) fn eat(
        &mut self,
        token: &str,
    ) -> bool {
        self.skip_ws();
        if self.rest().starts_with(token) {
            self.pos = self.pos.saturating_add(token.len());
            true
        } else {
            false
        }
    }

    pub(crate) fn at_end(&mut self) -> bool {
        self.peek().is_none()
    }

    pub(crate) fn ident(&mut self) -> Option<&'a str> {
        self.skip_ws();
        let rest = self.rest();
        let first = rest.chars().next()?;
        if !(first.is_alphabetic() || first == '_') {
            return None;
        }
        let len = rest
            .find(|c: char| !(c.is_alphanumeric() || c == '_'))
            .unwrap_or(rest.len());
        self.pos = self.pos.saturating_add(len);
        rest.get(..len)
    }

    fn number(&mut self) -> Result<Number, ParseError> {
        let rest = self.rest();
        let mut len = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        let mut is_float = false;
        let after = rest.get(len..).unwrap_or("");
        let mut chars = after.chars();
        if chars.next() == Some('.') && chars.next().is_some_and(|c| c.is_ascii_digit()) {
            is_float = true;
            let frac = after.get(1..).unwrap_or("");
            let frac_len = frac
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(frac.len());
            len = len.saturating_add(1).saturating_add(frac_len);
        }
        let after = rest.get(len..).unwrap_or("");
        if after.starts_with(['e', 'E']) {
            let exp = after.get(1..).unwrap_or("");
            let sign = usize::from(exp.starts_with(['+', '-']));
            let digits = exp.get(sign..).unwrap_or("");
            let digit_len = digits
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(digits.len());
            if digit_len > 0 {
                is_float = true;
                len = len
                    .saturating_add(1)
                    .saturating_add(sign)
                    .saturating_add(digit_len);
            }
        }
        let text = rest.get(..len).unwrap_or("");
        let parsed = if is_float {
            text.parse::<f64>().ok().map(Number::Float)
        } else {
            text.parse::<num_bigint::BigInt>().ok().map(Number::Int)
        };
        match parsed {
            | Some(n) => {
                self.pos = self.pos.saturating_add(len);
                Ok(n)
            },
            | None => self.err(format!("malformed number `{text}`")),
        }
    }

    fn mk(
        &self,
        op: OpId,
        children: Vec<Pat>,
    ) -> Pat {
        // Fold literal arithmetic so that `1/2` is one rational literal and
        // `-3` one integer literal.
        if let [Pat::Num(a), Pat::Num(b)] = children.as_slice() {
            let folded = match op {
                | core::ADD => Some(a.add(b)),
                | core::MUL => Some(a.mul(b)),
                | core::POW if a.is_exact() => a.pow(b),
                | _ => None,
            };
            if let Some(n) = folded {
                return Pat::Num(n);
            }
        }
        if !self.graph.ops().get(op).flags.has(OpFlags::ASSOCIATIVE) {
            return Pat::Node(op, children);
        }
        let mut flat = Vec::with_capacity(children.len());
        for child in children {
            match child {
                | Pat::Node(inner, grand) if inner == op => flat.extend(grand),
                | other => flat.push(other),
            }
        }
        if self.lenient && (op == core::ADD || op == core::MUL) {
            // A concrete term: merge its literals and drop the identity
            // element, so that `a - 2*b` and `a + (-2)*b` are one term.
            // Rule patterns are left exactly as written.
            let identity = Number::from(i64::from(op == core::MUL));
            let mut literal = identity.clone();
            flat.retain(|child| match child {
                | Pat::Num(n) => {
                    literal = if op == core::ADD {
                        literal.add(n)
                    } else {
                        literal.mul(n)
                    };
                    false
                },
                | _ => true,
            });
            if literal != identity || flat.is_empty() {
                flat.insert(0, Pat::Num(literal));
            }
            if flat.len() == 1 {
                return flat.swap_remove(0);
            }
        }
        Pat::Node(op, flat)
    }

    /// A sum, optionally followed by `= sum` to form an equation.
    pub(crate) fn expr(&mut self) -> Result<Pat, ParseError> {
        let lhs = self.sum()?;
        self.skip_ws();
        // `=>` and `<=>` belong to the rule syntax, not to the term.
        if self.rest().starts_with('=') && !self.rest().starts_with("=>") {
            self.eat("=");
            let rhs = self.sum()?;
            return Ok(Pat::Node(core::EQ, vec![lhs, rhs]));
        }
        Ok(lhs)
    }

    fn sum(&mut self) -> Result<Pat, ParseError> {
        let mut lhs = self.term()?;
        loop {
            if self.eat("+") {
                let rhs = self.term()?;
                lhs = self.mk(core::ADD, vec![lhs, rhs]);
            } else if self.peek() == Some('-') {
                self.eat("-");
                let rhs = self.term()?;
                let neg = self.mk(core::MUL, vec![Pat::Num(Number::from(-1)), rhs]);
                lhs = self.mk(core::ADD, vec![lhs, neg]);
            } else {
                return Ok(lhs);
            }
        }
    }

    fn term(&mut self) -> Result<Pat, ParseError> {
        let mut lhs = self.unary()?;
        loop {
            if self.eat("*") {
                let rhs = self.unary()?;
                lhs = self.mk(core::MUL, vec![lhs, rhs]);
            } else if self.eat("/") {
                let inv = match self.unary()? {
                    // 1/b^e is b^(-e) for a literal exponent.
                    | Pat::Node(core::POW, mut parts)
                        if matches!(parts.as_slice(), [_, Pat::Num(_)]) =>
                    {
                        if let Some(Pat::Num(e)) = parts.last_mut() {
                            *e = e.neg();
                        }
                        Pat::Node(core::POW, parts)
                    },
                    | rhs => self.mk(core::POW, vec![rhs, Pat::Num(Number::from(-1))]),
                };
                lhs = self.mk(core::MUL, vec![lhs, inv]);
            } else {
                return Ok(lhs);
            }
        }
    }

    fn unary(&mut self) -> Result<Pat, ParseError> {
        if self.peek() == Some('-') {
            self.eat("-");
            let inner = self.unary()?;
            return Ok(self.mk(core::MUL, vec![Pat::Num(Number::from(-1)), inner]));
        }
        let base = self.atom()?;
        if self.eat("^") {
            let exp = self.unary()?;
            return Ok(self.mk(core::POW, vec![base, exp]));
        }
        Ok(base)
    }

    fn atom(&mut self) -> Result<Pat, ParseError> {
        match self.peek() {
            | Some(c) if c.is_ascii_digit() => self.number().map(Pat::Num),
            | Some('(') => {
                self.eat("(");
                let inner = self.expr()?;
                if self.eat(")") {
                    Ok(inner)
                } else {
                    self.err("expected `)`")
                }
            },
            | Some('?') => {
                self.eat("?");
                match self.ident() {
                    | Some(name) => Ok(Pat::Var(self.vars.index(name))),
                    | None => self.err("expected a variable name after `?`"),
                }
            },
            | Some(_) => {
                let start = self.pos;
                let Some(name) = self.ident() else {
                    return self.err("expected a term");
                };
                if self.eat("(") {
                    let mut args = Vec::new();
                    if !self.eat(")") {
                        loop {
                            args.push(self.expr()?);
                            if self.eat(")") {
                                break;
                            }
                            if !self.eat(",") {
                                return self.err("expected `,` or `)`");
                            }
                        }
                    }
                    self.call(name, args, start)
                } else {
                    match self.graph.ops().lookup(name) {
                        | Some(op) if self.graph.ops().get(op).arity == Arity::Fixed(0) => {
                            Ok(Pat::Node(op, Vec::new()))
                        },
                        | _ => Ok(Pat::Exact(self.graph.sym(name))),
                    }
                }
            },
            | None => self.err("unexpected end of input"),
        }
    }

    fn call(
        &mut self,
        name: &str,
        args: Vec<Pat>,
        start: usize,
    ) -> Result<Pat, ParseError> {
        match self.graph.ops().lookup(name) {
            | Some(op) => {
                let desc = self.graph.ops().get(op);
                if !desc.arity.accepts(args.len()) || desc.flags.has(OpFlags::LEAF) {
                    return Err(ParseError {
                        message: format!("operator `{name}` cannot take {} arguments", args.len()),
                        offset: start,
                    });
                }
                Ok(self.mk(op, args))
            },
            | None if self.lenient => {
                let mut children = vec![Pat::Exact(self.graph.sym(name))];
                children.extend(args);
                Ok(Pat::Node(core::APPLY, children))
            },
            | None => Err(ParseError {
                message: format!("unknown operator `{name}`"),
                offset: start,
            }),
        }
    }
}

impl Pat {
    /// Parses a pattern. Unknown function names are errors.
    ///
    /// # Errors
    /// Returns a [`ParseError`] describing the first problem in `src`.
    pub fn parse(
        src: &str,
        graph: &mut Graph,
        vars: &mut VarNames,
    ) -> Result<Self, ParseError> {
        let mut parser = Parser::new(src, graph, vars, false);
        let pat = parser.expr()?;
        if parser.at_end() {
            Ok(pat)
        } else {
            parser.err("unexpected trailing input")
        }
    }

    /// Operator at the root, if the pattern is an application.
    #[must_use]
    pub const fn root_op(&self) -> Option<OpId> {
        match self {
            | Self::Node(op, _) => Some(*op),
            | _ => None,
        }
    }

    /// Collects the variables occurring in the pattern into `out`.
    pub fn vars(
        &self,
        out: &mut Vec<u32>,
    ) {
        match self {
            | Self::Var(v) => {
                if !out.contains(v) {
                    out.push(*v);
                }
            },
            | Self::Node(_, children) => children.iter().for_each(|c| c.vars(out)),
            | Self::Num(_) | Self::Exact(_) => {},
        }
    }

    /// Builds the term obtained by replacing every variable with its
    /// binding. Returns `None` if a variable is unbound.
    pub fn instantiate(
        &self,
        graph: &mut Graph,
        subst: &[NodeId],
    ) -> Option<NodeId> {
        match self {
            | Self::Var(v) => subst.get(*v as usize).copied().filter(|n| !n.is_none()),
            | Self::Num(n) => Some(graph.num(n.clone())),
            | Self::Exact(node) => Some(*node),
            | Self::Node(op, children) => {
                let mut built = Vec::with_capacity(children.len());
                for child in children {
                    built.push(child.instantiate(graph, subst)?);
                }
                graph.try_node(*op, &built)
            },
        }
    }

    /// Finds every way this pattern matches the e-node `root`, up to
    /// `limit` matches.
    #[must_use]
    pub fn matches(
        &self,
        graph: &Graph,
        root: NodeId,
        nvars: usize,
        limit: usize,
    ) -> Vec<Match> {
        let mut out = Vec::new();
        let Self::Node(op, pats) = self else {
            return out;
        };
        if graph.op(root) != *op {
            return out;
        }
        let mut matcher = Matcher {
            graph,
            limit,
            out: &mut out,
            root,
            multi: Vec::new(),
        };
        let mut subst = vec![NodeId::NONE; nvars];
        let flags = graph.ops().get(*op).flags;
        let children = graph.children(root);
        if flags.has(OpFlags::COMMUTATIVE) {
            let allow_rest = flags.has(OpFlags::ASSOCIATIVE);
            if pats.len() == children.len() || (allow_rest && pats.len() < children.len()) {
                let mut used = vec![false; children.len()];
                matcher.assign(
                    pats,
                    children,
                    0,
                    &mut used,
                    &mut Vec::new(),
                    &mut subst,
                    true,
                );
            }
        } else if pats.len() == children.len() {
            let mut goals: Vec<(&Self, NodeId)> =
                pats.iter().zip(children.iter().copied()).rev().collect();
            matcher.solve(&mut goals, &mut subst, &[]);
        }
        out
    }
}

struct Matcher<'g, 'o> {
    graph: &'g Graph,
    limit: usize,
    out: &'o mut Vec<Match>,
    root: NodeId,
    /// Variables currently bound to several operands of a nested node.
    multi: Vec<(u32, OpId, Vec<NodeId>)>,
}

impl Matcher<'_, '_> {
    /// Tries every injection of `pats[index..]` into the unused `children`.
    #[allow(clippy::too_many_arguments)]
    fn assign<'p>(
        &mut self,
        pats: &'p [Pat],
        children: &[NodeId],
        index: usize,
        used: &mut Vec<bool>,
        goals: &mut Vec<(&'p Pat, NodeId)>,
        subst: &mut Vec<NodeId>,
        is_root: bool,
    ) {
        if self.out.len() >= self.limit {
            return;
        }
        let Some(pat) = pats.get(index) else {
            let rest: Vec<NodeId> = if is_root {
                children
                    .iter()
                    .zip(used.iter())
                    .filter(|(_, u)| !**u)
                    .map(|(c, _)| *c)
                    .collect()
            } else {
                Vec::new()
            };
            let mut pending = goals.clone();
            self.solve(&mut pending, subst, &rest);
            return;
        };
        let mut tried: Vec<NodeId> = Vec::new();
        for (i, &child) in children.iter().enumerate() {
            if used.get(i).copied().unwrap_or(true) {
                continue;
            }
            // Children in the same class are interchangeable: try one.
            let class_rep = NodeId(self.graph.find(child).raw());
            if tried.contains(&class_rep) {
                continue;
            }
            tried.push(class_rep);
            if let Some(slot) = used.get_mut(i) {
                *slot = true;
            }
            goals.push((pat, child));
            self.assign(
                pats,
                children,
                index.saturating_add(1),
                used,
                goals,
                subst,
                is_root,
            );
            goals.pop();
            if let Some(slot) = used.get_mut(i) {
                *slot = false;
            }
        }
    }

    /// Discharges `goals` (a stack) under `subst`, emitting a match when
    /// none are left.
    fn solve(
        &mut self,
        goals: &mut Vec<(&Pat, NodeId)>,
        subst: &mut Vec<NodeId>,
        rest: &[NodeId],
    ) {
        if self.out.len() >= self.limit {
            return;
        }
        let Some((pat, node)) = goals.pop() else {
            self.out.push(Match {
                root: self.root,
                subst: subst.clone(),
                rest: rest.to_vec(),
                multi: self.multi.clone(),
            });
            return;
        };
        let graph = self.graph;
        let class = graph.find(node);
        match pat {
            | Pat::Var(v) => {
                let index = *v as usize;
                let bound = subst.get(index).copied().unwrap_or(NodeId::NONE);
                if let Some((_, op, operands)) = self.multi.iter().find(|(var, _, _)| var == v) {
                    // Bound to several operands elsewhere: this occurrence
                    // must be a node made of the same operands.
                    let mut want: Vec<u32> =
                        operands.iter().map(|&o| graph.find(o).raw()).collect();
                    want.sort_unstable();
                    let op = *op;
                    let same = graph.enodes(class).any(|enode| {
                        let mut have: Vec<u32> = graph
                            .children(enode)
                            .iter()
                            .map(|&c| graph.find(c).raw())
                            .collect();
                        have.sort_unstable();
                        graph.op(enode) == op && have == want
                    });
                    if same {
                        self.solve(goals, subst, rest);
                    }
                } else if bound.is_none() {
                    if let Some(slot) = subst.get_mut(index) {
                        *slot = node;
                    }
                    self.solve(goals, subst, rest);
                    if let Some(slot) = subst.get_mut(index) {
                        *slot = NodeId::NONE;
                    }
                } else if graph.find(bound) == class {
                    self.solve(goals, subst, rest);
                }
            },
            | Pat::Num(n) => {
                if graph.class_number(class) == Some(n) {
                    self.solve(goals, subst, rest);
                }
            },
            | Pat::Exact(target) => {
                if graph.find(*target) == class {
                    self.solve(goals, subst, rest);
                }
            },
            | Pat::Node(op, pats) => {
                let commutative = graph.ops().get(*op).flags.has(OpFlags::COMMUTATIVE);
                for enode in graph.enodes(class) {
                    if graph.op(enode) != *op {
                        continue;
                    }
                    let children = graph.children(enode);
                    if commutative && children.len() > pats.len() {
                        // The last operand, if a free variable, absorbs the
                        // operands the others leave over.
                        let absorber = match pats.split_last() {
                            | Some((Pat::Var(v), front))
                                if !front.is_empty()
                                    && graph.ops().get(*op).flags.has(OpFlags::ASSOCIATIVE)
                                    && subst.get(*v as usize).is_some_and(|b| b.is_none())
                                    && !self.multi.iter().any(|(var, _, _)| var == v) =>
                            {
                                Some((*v, front))
                            },
                            | _ => None,
                        };
                        if let Some((var, front)) = absorber {
                            let mut used = vec![false; children.len()];
                            let mut inner = goals.clone();
                            self.assign_nested(
                                front,
                                children,
                                0,
                                &mut used,
                                &mut inner,
                                subst,
                                rest,
                                Some((var, *op)),
                            );
                        }
                        continue;
                    }
                    if children.len() != pats.len() {
                        continue;
                    }
                    if commutative {
                        let mut used = vec![false; children.len()];
                        let mut inner = goals.clone();
                        self.assign_nested(
                            pats, children, 0, &mut used, &mut inner, subst, rest, None,
                        );
                    } else {
                        let mut inner = goals.clone();
                        inner.extend(pats.iter().zip(children.iter().copied()).rev());
                        self.solve(&mut inner, subst, rest);
                    }
                }
            },
        }
        goals.push((pat, node));
    }

    /// Like [`Self::assign`] for a nested commutative node: every child
    /// must be consumed — by a pattern or, with `absorb`, by the variable
    /// that takes what is left — and the outer `rest` is threaded through.
    #[allow(clippy::too_many_arguments)]
    fn assign_nested<'p>(
        &mut self,
        pats: &'p [Pat],
        children: &[NodeId],
        index: usize,
        used: &mut Vec<bool>,
        goals: &mut Vec<(&'p Pat, NodeId)>,
        subst: &mut Vec<NodeId>,
        rest: &[NodeId],
        absorb: Option<(u32, OpId)>,
    ) {
        if self.out.len() >= self.limit {
            return;
        }
        let Some(pat) = pats.get(index) else {
            let mut pending = goals.clone();
            if let Some((var, op)) = absorb {
                let left: Vec<NodeId> = children
                    .iter()
                    .zip(used.iter())
                    .filter(|(_, u)| !**u)
                    .map(|(c, _)| *c)
                    .collect();
                self.multi.push((var, op, left));
                self.solve(&mut pending, subst, rest);
                self.multi.pop();
            } else {
                self.solve(&mut pending, subst, rest);
            }
            return;
        };
        let mut tried: Vec<NodeId> = Vec::new();
        for (i, &child) in children.iter().enumerate() {
            if used.get(i).copied().unwrap_or(true) {
                continue;
            }
            let class_rep = NodeId(self.graph.find(child).raw());
            if tried.contains(&class_rep) {
                continue;
            }
            tried.push(class_rep);
            if let Some(slot) = used.get_mut(i) {
                *slot = true;
            }
            goals.push((pat, child));
            self.assign_nested(
                pats,
                children,
                index.saturating_add(1),
                used,
                goals,
                subst,
                rest,
                absorb,
            );
            goals.pop();
            if let Some(slot) = used.get_mut(i) {
                *slot = false;
            }
        }
    }
}

impl Graph {
    /// Parses a concrete term (no pattern variables) and interns it.
    ///
    /// Unknown function names become applications of an undetermined
    /// function: `f(x)` parses as `apply(f, x)`.
    ///
    /// # Errors
    /// Returns a [`ParseError`] for malformed input or if the text contains
    /// pattern variables.
    pub fn parse(
        &mut self,
        src: &str,
    ) -> Result<NodeId, ParseError> {
        let mut vars = VarNames::default();
        let pat = {
            let mut parser = Parser::new(src, self, &mut vars, true);
            let pat = parser.expr()?;
            if !parser.at_end() {
                return parser.err("unexpected trailing input");
            }
            pat
        };
        if !vars.is_empty() {
            return Err(ParseError {
                message: "pattern variables are not allowed in a term".to_owned(),
                offset: 0,
            });
        }
        pat.instantiate(self, &[]).ok_or_else(|| ParseError {
            message: "term could not be built".to_owned(),
            offset: 0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::OpDescriptor;

    fn setup() -> Graph {
        let mut g = Graph::new();
        for name in ["sin", "cos"] {
            assert!(
                g.ops_mut()
                    .register(OpDescriptor::new(name, Arity::Fixed(1)))
                    .is_ok()
            );
        }
        g
    }

    fn pat(
        g: &mut Graph,
        vars: &mut VarNames,
        src: &str,
    ) -> Pat {
        Pat::parse(src, g, vars).unwrap_or(Pat::Var(u32::MAX))
    }

    #[test]
    fn sugar_desugars_to_canonical_forms() {
        let mut g = setup();
        let a = g.parse("x - y").ok();
        let b = g.parse("x + (-1)*y").ok();
        assert!(a.is_some());
        assert_eq!(a, b);
        let c = g.parse("x / y").ok();
        let d = g.parse("x * y^(-1)").ok();
        assert_eq!(c, d);
        let half = g.parse("1/2").ok();
        let expect = Number::fraction(1, 2).map(|n| g.num(n));
        assert_eq!(half, expect);
        assert_eq!(g.parse("-3").ok(), Some(g.int(-3)));
        assert_eq!(g.parse("2.5e1").ok(), Some(g.float(25.0)));
    }

    #[test]
    fn precedence() {
        let mut g = setup();
        assert_eq!(
            g.parse("a + b * c ^ d").ok(),
            g.parse("a + (b * (c ^ d))").ok()
        );
        assert_eq!(g.parse("-x^2").ok(), g.parse("-(x^2)").ok());
        assert_eq!(
            g.parse("2^-1").ok(),
            Number::fraction(1, 2).map(|n| g.num(n))
        );
        assert_eq!(g.parse("a^b^c").ok(), g.parse("a^(b^c)").ok());
    }

    #[test]
    fn errors_are_reported() {
        let mut g = setup();
        let mut vars = VarNames::default();
        assert!(Pat::parse("foo(?x)", &mut g, &mut vars).is_err());
        assert!(Pat::parse("sin(?x, ?y)", &mut g, &mut vars).is_err());
        assert!(Pat::parse("sin(?x", &mut g, &mut vars).is_err());
        assert!(g.parse("?x + 1").is_err());
        assert!(g.parse("1 +").is_err());
        assert!(g.parse("f(x)").is_ok(), "lenient mode builds apply(f, x)");
    }

    #[test]
    fn nonlinear_match() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "sin(?x)^2 + cos(?x)^2");
        let good = g.parse("sin(a)^2 + cos(a)^2").unwrap_or(NodeId::NONE);
        let bad = g.parse("sin(a)^2 + cos(b)^2").unwrap_or(NodeId::NONE);
        let a = g.sym("a");
        let m = p.matches(&g, good, vars.len(), 16);
        assert_eq!(m.len(), 1);
        assert_eq!(m.first().map(|m| m.subst.clone()), Some(vec![a]));
        assert!(p.matches(&g, bad, vars.len(), 16).is_empty());
    }

    #[test]
    fn match_modulo_equalities() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "sin(?x)^2 + cos(?x)^2");
        let node = g.parse("sin(a)^2 + cos(b)^2").unwrap_or(NodeId::NONE);
        assert!(p.matches(&g, node, vars.len(), 16).is_empty());
        let (a, b) = (g.sym("a"), g.sym("b"));
        g.union(a, b);
        g.rebuild();
        assert_eq!(p.matches(&g, node, vars.len(), 16).len(), 1);
    }

    #[test]
    fn ac_match_collects_rest() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "sin(?x)^2 + cos(?x)^2");
        let node = g
            .parse("u + sin(a)^2 + v + cos(a)^2")
            .unwrap_or(NodeId::NONE);
        let m = p.matches(&g, node, vars.len(), 16);
        assert_eq!(m.len(), 1);
        let mut rest = m.first().map(|m| m.rest.clone()).unwrap_or_default();
        rest.sort_unstable();
        let mut expect = vec![g.sym("u"), g.sym("v")];
        expect.sort_unstable();
        assert_eq!(rest, expect);
    }

    #[test]
    fn commutative_match_enumerates_assignments() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "?a * ?b");
        let node = g.parse("x * y").unwrap_or(NodeId::NONE);
        assert_eq!(p.matches(&g, node, vars.len(), 16).len(), 2);
        let square = g.parse("x * x").unwrap_or(NodeId::NONE);
        assert_eq!(
            p.matches(&g, square, vars.len(), 16).len(),
            1,
            "equal children are tried once"
        );
    }

    #[test]
    fn number_patterns_use_class_constants() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "?x ^ 2");
        let node = g.parse("y ^ k").unwrap_or(NodeId::NONE);
        assert!(p.matches(&g, node, vars.len(), 4).is_empty());
        let (k, two) = (g.sym("k"), g.int(2));
        g.union(k, two);
        g.rebuild();
        assert_eq!(p.matches(&g, node, vars.len(), 4).len(), 1);
    }

    #[test]
    fn nested_ac_patterns_absorb_into_the_last_variable() {
        let mut g = setup();
        let mut vars = VarNames::default();
        // f(?c * ?r): ?c takes one factor, ?r the others.
        let p = pat(&mut g, &mut vars, "sin(?c * ?r)");
        let node = g.parse("sin(2 * x * y)").unwrap_or(NodeId::NONE);
        let found = p.matches(&g, node, vars.len(), 16);
        assert_eq!(found.len(), 3, "each factor in turn is ?c");
        let mut seen = Vec::new();
        for m in &found {
            let subst = m.materialize(&mut g);
            seen.push(format!("{} | {}", g.display(subst[0]), g.display(subst[1])));
        }
        seen.sort();
        assert_eq!(seen, vec!["2 | x*y", "x | 2*y", "y | 2*x"]);
        // Exact arity still binds single operands.
        let node = g.parse("sin(2 * x)").unwrap_or(NodeId::NONE);
        assert_eq!(p.matches(&g, node, vars.len(), 16).len(), 2);
        // A non-variable last operand cannot absorb.
        let mut vars = VarNames::default();
        let strict = pat(&mut g, &mut vars, "sin(?c * cos(?r))");
        let node = g.parse("sin(2 * x * cos(y))").unwrap_or(NodeId::NONE);
        assert!(strict.matches(&g, node, vars.len(), 16).is_empty());
    }

    #[test]
    fn absorbed_variables_compare_as_products() {
        let mut g = setup();
        let mut vars = VarNames::default();
        // ?r occurs twice: once absorbing, once as a whole node.
        let p = pat(&mut g, &mut vars, "sin(?c * ?r) ^ cos(?r)");
        let good = g
            .parse("sin(2 * x * y) ^ cos(x * y)")
            .unwrap_or(NodeId::NONE);
        let bad = g
            .parse("sin(2 * x * y) ^ cos(x * z)")
            .unwrap_or(NodeId::NONE);
        assert_eq!(p.matches(&g, good, vars.len(), 16).len(), 1);
        assert!(p.matches(&g, bad, vars.len(), 16).is_empty());
    }

    #[test]
    fn instantiate_round_trip() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "sin(?x) * ?y");
        let (a, b) = (g.sym("a"), g.sym("b"));
        let built = p.instantiate(&mut g, &[a, b]);
        assert_eq!(built, g.parse("sin(a) * b").ok());
        assert_eq!(p.instantiate(&mut g, &[a, NodeId::NONE]), None);
    }

    #[test]
    fn match_limit_is_respected() {
        let mut g = setup();
        let mut vars = VarNames::default();
        let p = pat(&mut g, &mut vars, "?a + ?b");
        let node = g.parse("p + q + r + s + t").unwrap_or(NodeId::NONE);
        assert_eq!(p.matches(&g, node, vars.len(), 5).len(), 5);
    }
}
