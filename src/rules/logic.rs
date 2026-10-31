//! Propositional logic and comparisons.
//!
//! Boolean operators evaluate to exactly `1.0` or `0.0` (non-zero counts as
//! true on input). Identities that would leave a bare operand standing
//! alone (`and(a, true) = a`) are a kernel, because a random real sample of
//! `a` would refute them; the rest are declarative rewrites.
//!
//! The *variables* of a formula are its maximal non-logical subterms:
//! symbols, comparisons, anything else. Each distinct one is an atom.
//! Normal-form requests (`cnf`, `dnf`, `nnf`, `simplify_logic`) and
//! decision requests (`satisfiable`, `tautology`, `truth_table`) are
//! reduced by one exact kernel.
//!
//! The module also compiles formulas into binary decision diagrams for
//! `count_models`, `sat_model`, `all_models`, `bdd`, `bdd_size`, `exists`,
//! `forall`, `restrict`, `entails`, `equivalent` and `prime_implicants`
//! (see its table).

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::graph::Arity;
use crate::graph::ClassId;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimTerm;
use crate::graph::op::EvalFn;
use crate::graph::op::core;
use crate::graph::rule::Installer;


use super::arith::arith;

/// The propositional-logic rule set.
#[must_use]
pub fn logic() -> RuleSet {
    RuleSet::new("logic", install).needs(arith())
}

/// Deepest formula the normal-form kernels will look into.
const MAX_DEPTH: usize = 256;
/// Most variables a truth table or a satisfiability check enumerates.
const MAX_ENUMERATED: usize = 20;
/// Most variables `simplify_logic` and `truth_table` accept.
const MAX_MINIMISED: usize = 12;
/// Largest normal form (in nodes or terms) the distribution steps build.
const MAX_SIZE: usize = 20_000;
/// Work bound of the exact set-cover search.
const MAX_SEARCH_STEPS: usize = 200_000;

fn is_true(x: f64) -> bool {
    x != 0.0 && !x.is_nan()
}

const fn bit(b: bool) -> f64 {
    if b { 1.0 } else { 0.0 }
}

fn compare2(
    args: &[f64],
    f: fn(f64, f64) -> bool,
) -> f64 {
    match args {
        | [a, b] => bit(f(*a, *b)),
        | _ => f64::NAN,
    }
}

/// The operators the kernels need to recognise.
#[derive(Copy, Clone)]
struct Ops {
    and: OpId,
    or: OpId,
    xor: OpId,
    not: OpId,
    implies: OpId,
    iff: OpId,
    tru: OpId,
    fals: OpId,
}

/// A comparison operator.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Ne,
    Eq,
}

impl Cmp {
    /// The same comparison with its operands swapped.
    const fn flip(self) -> Self {
        match self {
            | Self::Lt => Self::Gt,
            | Self::Le => Self::Ge,
            | Self::Gt => Self::Lt,
            | Self::Ge => Self::Le,
            | other => other,
        }
    }

    const fn holds(
        self,
        order: Ordering,
    ) -> bool {
        match self {
            | Self::Lt => order.is_lt(),
            | Self::Le => order.is_le(),
            | Self::Gt => order.is_gt(),
            | Self::Ge => order.is_ge(),
            | Self::Ne => order.is_ne(),
            | Self::Eq => order.is_eq(),
        }
    }

    /// The answer to `self(a, 0)` when only the facts of `a` are known.
    fn against_zero(
        self,
        facts: Facts,
    ) -> Option<bool> {
        let (yes, no) = match self {
            | Self::Gt => (Facts::POSITIVE, Facts::NONPOSITIVE),
            | Self::Ge => (Facts::NONNEGATIVE, Facts::NEGATIVE),
            | Self::Lt => (Facts::NEGATIVE, Facts::NONNEGATIVE),
            | Self::Le => (Facts::NONPOSITIVE, Facts::POSITIVE),
            | Self::Ne => return facts.has(Facts::NONZERO).then_some(true),
            | Self::Eq => return facts.has(Facts::NONZERO).then_some(false),
        };
        if facts.has(yes) {
            Some(true)
        } else if facts.has(no) {
            Some(false)
        } else {
            None
        }
    }
}

/// A count (`count_models`, `bdd_size`) is a pure number. The argument is a
/// formula or a diagram and keeps whatever dimension its Boolean structure
/// has, so it is deliberately unconstrained.
#[allow(clippy::missing_const_for_fn, clippy::unnecessary_wraps)] // `DimRule::Custom` signature
fn dimensionless_value(
    _ctx: &DimCtx<'_>,
    _args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    Ok(DimTerm::dimensionless())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let ac = OpFlags::COMMUTATIVE
        .with(OpFlags::ASSOCIATIVE)
        .with(OpFlags::PREDICATE);
    let and = i.op(OpDescriptor::new("and", Arity::Variadic)
        .flags(ac)
        .eval(|a| bit(a.iter().all(|&x| is_true(x)))))?;
    let or = i.op(OpDescriptor::new("or", Arity::Variadic)
        .flags(ac)
        .eval(|a| bit(a.iter().any(|&x| is_true(x)))))?;
    let xor = i.op(OpDescriptor::new("xor", Arity::Variadic)
        .flags(ac)
        .eval(|a| bit(a.iter().filter(|&&x| is_true(x)).count() % 2 == 1)))?;
    let not = i.op(OpDescriptor::new("not", Arity::Fixed(1))
        .flags(OpFlags::PREDICATE)
        .eval(|a| a.first().map_or(f64::NAN, |&x| bit(!is_true(x)))))?;
    // Costlier than their expansions, so extraction prefers and/or/not.
    let implies = i.op(OpDescriptor::new("implies", Arity::Fixed(2))
        .flags(OpFlags::PREDICATE)
        .cost(12)
        .eval(|a| {
            match a {
                | [p, q] => bit(!is_true(*p) || is_true(*q)),
                | _ => f64::NAN,
            }
        }))?;
    let iff = i.op(OpDescriptor::new("iff", Arity::Fixed(2))
        .flags(OpFlags::PREDICATE)
        .cost(30)
        .eval(|a| {
            match a {
                | [p, q] => bit(is_true(*p) == is_true(*q)),
                | _ => f64::NAN,
            }
        }))?;
    let tru = i.op(OpDescriptor::new("true", Arity::Fixed(0)).eval(|_| 1.0))?;
    let fals = i.op(OpDescriptor::new("false", Arity::Fixed(0)).eval(|_| 0.0))?;
    let ops = Ops {
        and,
        or,
        xor,
        not,
        implies,
        iff,
        tru,
        fals,
    };

    let lt: EvalFn = |a| compare2(a, |x, y| x < y);
    let le: EvalFn = |a| compare2(a, |x, y| x <= y);
    let gt: EvalFn = |a| compare2(a, |x, y| x > y);
    let ge: EvalFn = |a| compare2(a, |x, y| x >= y);
    let ne: EvalFn = |a| compare2(a, |x, y| x.total_cmp(&y).is_ne());
    let mut comparisons = Vec::new();
    for (name, eval, cmp) in [
        ("lt", lt, Cmp::Lt),
        ("le", le, Cmp::Le),
        ("gt", gt, Cmp::Gt),
        ("ge", ge, Cmp::Ge),
        ("ne", ne, Cmp::Ne),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(2))
            .flags(OpFlags::PREDICATE)
            .eval(eval))?;
        comparisons.push((op, cmp));
    }
    comparisons.push((core::EQ, Cmp::Eq));
    i.kernel("logic/compare", Tier::Normalize, Compare { comparisons });

    let request = |name: &str| {
        OpDescriptor::new(name, Arity::Fixed(1))
            .flags(OpFlags::HEAVY)
            .cost(100)
    };
    let identity: EvalFn = |a| a.first().copied().unwrap_or(f64::NAN);
    let mut forms = Vec::new();
    for (name, form) in [
        ("cnf", Form::Cnf),
        ("dnf", Form::Dnf),
        ("nnf", Form::Nnf),
        ("simplify_logic", Form::Minimal),
    ] {
        let op = i.op(request(name).eval(identity))?;
        // A formula in another normal form is still a formula.
        i.dim(op, DimRule::NonMetric);
        forms.push((op, form));
    }
    for (name, form) in [
        ("satisfiable", Form::Satisfiable),
        ("tautology", Form::Tautology),
        ("truth_table", Form::Table),
    ] {
        let op = i.op(request(name))?;
        // A decision, a table of truth values or a model is not a scalar.
        i.dim(op, DimRule::NonMetric);
        forms.push((op, form));
    }

    i.kernel("logic/constants", Tier::Normalize, Constants { ops });
    i.kernel("logic/boolean", Tier::Normalize, Boolean { ops });
    i.kernel("logic/forms", Tier::Reduce, Forms { ops, forms });
    install_ext(i, ops)?;

    i.rewrites(
        Tier::Normalize,
        &[
            "logic/and-false: and(?a, false) => false",
            "logic/or-true: or(?a, true) => true",
            "logic/and-complement: and(?a, not(?a)) => false",
            "logic/or-complement: or(?a, not(?a)) => true",
            "logic/xor-self: xor(?a, ?a) => false",
            "logic/implies: implies(?a, ?b) => or(not(?a), ?b)",
            "logic/iff: iff(?a, ?b) => or(and(?a, ?b), and(not(?a), not(?b)))",
            "logic/not-lt: not(lt(?a, ?b)) => ge(?a, ?b)",
            "logic/not-le: not(le(?a, ?b)) => gt(?a, ?b)",
            "logic/not-gt: not(gt(?a, ?b)) => le(?a, ?b)",
            "logic/not-ge: not(ge(?a, ?b)) => lt(?a, ?b)",
        ],
    )?;
    i.rewrites(
        Tier::Explore,
        &[
            "logic/de-morgan-and: not(and(?a, ?b)) => or(not(?a), not(?b))",
            "logic/de-morgan-or: not(or(?a, ?b)) => and(not(?a), not(?b))",
            "logic/de-morgan-and-rev: or(not(?a), not(?b)) => not(and(?a, ?b))",
            "logic/de-morgan-or-rev: and(not(?a), not(?b)) => not(or(?a, ?b))",
        ],
    )
}

fn truth_literal(
    graph: &mut Graph,
    value: bool,
) -> NodeId {
    graph.lit(Payload::Bool(value))
}

/// The truth value `node` is known to equal, if any.
fn constant(
    graph: &Graph,
    ops: &Ops,
    node: NodeId,
) -> Option<bool> {
    let of = |n: NodeId| {
        match graph.payload(n) {
            | Some(Payload::Bool(b)) if graph.op(n) == core::LIT => Some(*b),
            | _ if graph.op(n) == ops.tru => Some(true),
            | _ if graph.op(n) == ops.fals => Some(false),
            | _ => None,
        }
    };
    of(node).or_else(|| graph.enodes(graph.find(node)).find_map(of))
}

/// An application of `op` that `node` is known to equal.
fn application(
    graph: &Graph,
    op: OpId,
    node: NodeId,
) -> Option<NodeId> {
    if graph.op(node) == op {
        return Some(node);
    }
    graph.enodes(graph.find(node)).find(|&n| graph.op(n) == op)
}

/// `tru` and `fals` are the literals `true` and `false` in disguise.
struct Constants {
    ops: Ops,
}

impl Kernel for Constants {
    fn ops(&self) -> Vec<OpId> {
        vec![self.ops.tru, self.ops.fals]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let value = cx.graph.op(node) == self.ops.tru;
        Outcome::Equal(truth_literal(cx.graph, value))
    }
}

/// Boolean simplification on the concrete term: units, idempotence,
/// absorption, double negation. Operands are assumed boolean-valued.
struct Boolean {
    ops: Ops,
}

impl Boolean {
    /// `and` (`conjunction`) or `or` of its operands with units removed,
    /// duplicates merged and absorbed operands dropped.
    fn lattice(
        &self,
        graph: &mut Graph,
        node: NodeId,
        conjunction: bool,
    ) -> Outcome {
        let ops = &self.ops;
        let (op, dual) = if conjunction {
            (ops.and, ops.or)
        } else {
            (ops.or, ops.and)
        };
        let kids = graph.children(node).to_vec();
        let mut kept: Vec<NodeId> = Vec::with_capacity(kids.len());
        for &kid in &kids {
            match constant(graph, ops, kid) {
                // `and(.., false)` and `or(.., true)`.
                | Some(value) if value != conjunction => {
                    return Outcome::Equal(truth_literal(graph, value));
                },
                | Some(_) => {},
                | None if kept.iter().any(|&k| graph.same(k, kid)) => {},
                | None => kept.push(kid),
            }
        }
        // `a and (a or b) = a`, `a or (a and b) = a`.
        let absorbed: Vec<bool> = kept
            .iter()
            .map(|&k| {
                application(graph, dual, k).is_some_and(|inner| {
                    kept.iter().any(|&other| {
                        other != k && graph.children(inner).iter().any(|&c| graph.same(c, other))
                    })
                })
            })
            .collect();
        let kept: Vec<NodeId> = kept
            .into_iter()
            .zip(absorbed)
            .filter_map(|(k, gone)| (!gone).then_some(k))
            .collect();
        if kept.len() == kids.len() {
            return Outcome::Pass;
        }
        Outcome::Equal(match kept.as_slice() {
            | [] => truth_literal(graph, conjunction),
            | [only] => *only,
            | _ => graph.node(op, &kept),
        })
    }

    /// `xor` with units and pairs cancelled: `a xor a = false`,
    /// `a xor true = not(a)`.
    fn exclusive(
        &self,
        graph: &mut Graph,
        node: NodeId,
    ) -> Outcome {
        let kids = graph.children(node).to_vec();
        let mut kept: Vec<NodeId> = Vec::with_capacity(kids.len());
        let mut negated = false;
        for &kid in &kids {
            match constant(graph, &self.ops, kid) {
                | Some(value) => negated ^= value,
                | None => {
                    match kept.iter().position(|&k| graph.same(k, kid)) {
                        | Some(at) => {
                            kept.remove(at);
                        },
                        | None => kept.push(kid),
                    }
                },
            }
        }
        if kept.len() == kids.len() && !negated {
            return Outcome::Pass;
        }
        let body = match kept.as_slice() {
            | [] => truth_literal(graph, false),
            | [only] => *only,
            | _ => graph.node(self.ops.xor, &kept),
        };
        Outcome::Equal(if negated {
            graph.node(self.ops.not, &[body])
        } else {
            body
        })
    }
}

impl Kernel for Boolean {
    fn ops(&self) -> Vec<OpId> {
        vec![self.ops.and, self.ops.or, self.ops.xor, self.ops.not]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        if op == self.ops.and {
            self.lattice(graph, node, true)
        } else if op == self.ops.or {
            self.lattice(graph, node, false)
        } else if op == self.ops.xor {
            self.exclusive(graph, node)
        } else {
            let Some(&arg) = graph.children(node).first() else {
                return Outcome::Pass;
            };
            if let Some(value) = constant(graph, &self.ops, arg) {
                return Outcome::Equal(truth_literal(graph, !value));
            }
            match application(graph, self.ops.not, arg)
                .and_then(|n| graph.children(n).first().copied())
            {
                | Some(inner) => Outcome::Equal(inner),
                | None => Outcome::Pass,
            }
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Decides comparisons of literal numbers, and of a term with zero when its
/// sign is known.
struct Compare {
    comparisons: Vec<(OpId, Cmp)>,
}

/// Order of two literals; exact when both are exact.
fn order(
    a: &Number,
    b: &Number,
) -> Option<Ordering> {
    match (a.to_rational(), b.to_rational()) {
        | (Some(x), Some(y)) => Some(x.cmp(&y)),
        | _ => a.to_f64().partial_cmp(&b.to_f64()),
    }
}

impl Kernel for Compare {
    fn ops(&self) -> Vec<OpId> {
        self.comparisons.iter().map(|c| c.0).collect()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        let (Some(&(_, cmp)), &[a, b]) = (
            self.comparisons.iter().find(|c| c.0 == op),
            graph.children(node),
        ) else {
            return Outcome::Pass;
        };
        let decided = match (graph.number_of(a), graph.number_of(b)) {
            | (Some(x), Some(y)) => order(x, y).map(|o| cmp.holds(o)),
            | (None, Some(y)) if y.is_zero() => cmp.against_zero(graph.facts(a)),
            | (Some(x), None) if x.is_zero() => cmp.flip().against_zero(graph.facts(b)),
            | _ => None,
        };
        decided.map_or(Outcome::Pass, |value| {
            Outcome::Equal(truth_literal(graph, value))
        })
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// A formula over numbered atoms.
#[derive(Clone, Debug)]
enum Formula {
    Const(bool),
    Atom(usize),
    Not(Box<Self>),
    And(Vec<Self>),
    Or(Vec<Self>),
    Xor(Vec<Self>),
    Implies(Box<Self>, Box<Self>),
    Iff(Box<Self>, Box<Self>),
}

/// A formula in negation normal form.
#[derive(Clone, Debug)]
enum Nnf {
    Const(bool),
    Lit(usize, bool),
    And(Vec<Self>),
    Or(Vec<Self>),
}

/// A conjunction or disjunction of literals `(atom, positive)`.
type Term = BTreeSet<(usize, bool)>;

/// Reads a term as a [`Formula`], collecting its atoms.
struct Reader<'a> {
    graph: &'a Graph,
    ops: &'a Ops,
    classes: Vec<ClassId>,
    nodes: Vec<NodeId>,
}

impl Reader<'_> {
    fn read(
        &mut self,
        node: NodeId,
        depth: usize,
    ) -> Option<Formula> {
        if depth > MAX_DEPTH {
            return None;
        }
        if let Some(value) = constant(self.graph, self.ops, node) {
            return Some(Formula::Const(value));
        }
        let (op, ops) = (self.graph.op(node), self.ops);
        let kids = self.graph.children(node);
        let all = |this: &mut Self| -> Option<Vec<Formula>> {
            kids.iter()
                .map(|&k| this.read(k, depth.saturating_add(1)))
                .collect()
        };
        Some(if op == ops.and {
            Formula::And(all(self)?)
        } else if op == ops.or {
            Formula::Or(all(self)?)
        } else if op == ops.xor {
            Formula::Xor(all(self)?)
        } else if let (true, &[a]) = (op == ops.not, kids) {
            Formula::Not(Box::new(self.read(a, depth.saturating_add(1))?))
        } else if let (true, &[a, b]) = (op == ops.implies, kids) {
            let (a, b) = (
                self.read(a, depth.saturating_add(1))?,
                self.read(b, depth.saturating_add(1))?,
            );
            Formula::Implies(Box::new(a), Box::new(b))
        } else if let (true, &[a, b]) = (op == ops.iff, kids) {
            let (a, b) = (
                self.read(a, depth.saturating_add(1))?,
                self.read(b, depth.saturating_add(1))?,
            );
            Formula::Iff(Box::new(a), Box::new(b))
        } else {
            let class = self.graph.find(node);
            let at = self
                .classes
                .iter()
                .position(|&c| c == class)
                .unwrap_or_else(|| {
                    self.classes.push(class);
                    self.nodes.push(node);
                    self.classes.len().saturating_sub(1)
                });
            Formula::Atom(at)
        })
    }
}

impl Formula {
    fn remap(
        &mut self,
        rank: &[usize],
    ) {
        match self {
            | Self::Const(_) => {},
            | Self::Atom(v) => *v = rank.get(*v).copied().unwrap_or(*v),
            | Self::Not(x) => x.remap(rank),
            | Self::And(xs) | Self::Or(xs) | Self::Xor(xs) => {
                for x in xs {
                    x.remap(rank);
                }
            },
            | Self::Implies(a, b) | Self::Iff(a, b) => {
                a.remap(rank);
                b.remap(rank);
            },
        }
    }

    /// Truth value when atom `v` is `mask >> v & 1`.
    fn holds(
        &self,
        mask: u32,
    ) -> bool {
        match self {
            | Self::Const(b) => *b,
            | Self::Atom(v) => mask >> v & 1 == 1,
            | Self::Not(x) => !x.holds(mask),
            | Self::And(xs) => xs.iter().all(|x| x.holds(mask)),
            | Self::Or(xs) => xs.iter().any(|x| x.holds(mask)),
            | Self::Xor(xs) => xs.iter().filter(|x| x.holds(mask)).count() % 2 == 1,
            | Self::Implies(a, b) => !a.holds(mask) || b.holds(mask),
            | Self::Iff(a, b) => a.holds(mask) == b.holds(mask),
        }
    }

    /// Negation normal form of the formula (`negate`: of its negation).
    /// `budget` bounds the number of nodes built.
    fn nnf(
        &self,
        negate: bool,
        budget: &mut usize,
    ) -> Option<Nnf> {
        *budget = budget.checked_sub(1)?;
        let junction = |xs: &[Self], conjunction: bool, budget: &mut usize| {
            let parts: Option<Vec<Nnf>> = xs.iter().map(|x| x.nnf(negate, budget)).collect();
            parts.map(|p| {
                if conjunction == negate {
                    Nnf::or(p)
                } else {
                    Nnf::and(p)
                }
            })
        };
        match self {
            | Self::Const(b) => Some(Nnf::Const(*b != negate)),
            | Self::Atom(v) => Some(Nnf::Lit(*v, !negate)),
            | Self::Not(x) => x.nnf(!negate, budget),
            | Self::And(xs) => junction(xs, true, budget),
            | Self::Or(xs) => junction(xs, false, budget),
            | Self::Implies(a, b) => {
                Some(if negate {
                    Nnf::and(vec![a.nnf(false, budget)?, b.nnf(true, budget)?])
                } else {
                    Nnf::or(vec![a.nnf(true, budget)?, b.nnf(false, budget)?])
                })
            },
            | Self::Iff(a, b) => {
                let (pa, na) = (a.nnf(false, budget)?, a.nnf(true, budget)?);
                let (pb, nb) = (b.nnf(false, budget)?, b.nnf(true, budget)?);
                // a <-> b is (!a | b) & (a | !b); its negation (a | b) & (!a | !b).
                Some(if negate {
                    Nnf::and(vec![Nnf::or(vec![pa, pb]), Nnf::or(vec![na, nb])])
                } else {
                    Nnf::and(vec![Nnf::or(vec![na, pb]), Nnf::or(vec![pa, nb])])
                })
            },
            | Self::Xor(xs) => {
                match xs.as_slice() {
                    | [] => Some(Nnf::Const(negate)),
                    | [only] => only.nnf(negate, budget),
                    // a ^ rest = !(a <-> rest)
                    | [first, rest @ ..] => {
                        Self::Iff(Box::new(first.clone()), Box::new(Self::Xor(rest.to_vec())))
                            .nnf(!negate, budget)
                    },
                }
            },
        }
    }
}

impl Nnf {
    /// A conjunction with constants folded in.
    fn and(parts: Vec<Self>) -> Self {
        Self::fold(parts, true)
    }

    /// A disjunction with constants folded in.
    fn or(parts: Vec<Self>) -> Self {
        Self::fold(parts, false)
    }

    fn fold(
        parts: Vec<Self>,
        conjunction: bool,
    ) -> Self {
        let mut kept = Vec::with_capacity(parts.len());
        for part in parts {
            match part {
                | Self::Const(b) if b == conjunction => {},
                | Self::Const(b) => return Self::Const(b),
                | other => kept.push(other),
            }
        }
        match kept.len() {
            | 0 => Self::Const(conjunction),
            | 1 => kept.pop().unwrap_or(Self::Const(conjunction)),
            | _ if conjunction => Self::And(kept),
            | _ => Self::Or(kept),
        }
    }

    /// Sum-of-products terms, without contradictory or subsumed terms.
    fn dnf(&self) -> Option<Vec<Term>> {
        let terms = match self {
            | Self::Const(true) => vec![Term::new()],
            | Self::Const(false) => Vec::new(),
            | Self::Lit(v, positive) => vec![Term::from([(*v, *positive)])],
            | Self::Or(parts) => {
                let mut all = Vec::new();
                for part in parts {
                    all.extend(part.dnf()?);
                }
                all
            },
            | Self::And(parts) => {
                let mut acc = vec![Term::new()];
                for part in parts {
                    let right = part.dnf()?;
                    let mut next = Vec::with_capacity(acc.len().saturating_mul(right.len()));
                    for left in &acc {
                        for r in &right {
                            let both: Term = left.union(r).copied().collect();
                            if !both.iter().any(|&(v, p)| both.contains(&(v, !p))) {
                                next.push(both);
                            }
                        }
                    }
                    if next.len() > MAX_SIZE {
                        return None;
                    }
                    acc = absorb(next);
                }
                acc
            },
        };
        (terms.len() <= MAX_SIZE).then(|| absorb(terms))
    }
}

/// Removes duplicate terms and terms that contain another term, and sorts
/// by size for a deterministic result.
fn absorb(mut terms: Vec<Term>) -> Vec<Term> {
    terms.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    let mut kept: Vec<Term> = Vec::with_capacity(terms.len());
    for term in terms {
        if !kept.iter().any(|k| k.is_subset(&term)) {
            kept.push(term);
        }
    }
    kept
}

/// Builds terms back into the graph.
struct Builder<'a> {
    graph: &'a mut Graph,
    ops: &'a Ops,
    atoms: &'a [NodeId],
}

impl Builder<'_> {
    fn literal(
        &mut self,
        atom: usize,
        positive: bool,
    ) -> Option<NodeId> {
        let node = *self.atoms.get(atom)?;
        Some(if positive {
            node
        } else {
            self.graph.node(self.ops.not, &[node])
        })
    }

    /// `and` (or `or`) of `items`; the unit for none, the item for one.
    fn join(
        &mut self,
        conjunction: bool,
        items: &[NodeId],
    ) -> NodeId {
        match items {
            | [] => truth_literal(self.graph, conjunction),
            | [only] => *only,
            | _ => {
                self.graph.node(
                    if conjunction {
                        self.ops.and
                    } else {
                        self.ops.or
                    },
                    items,
                )
            },
        }
    }

    /// A junction of junctions of literals: `outer(inner(literals)...)`,
    /// with every literal negated if `negate`.
    fn normal_form(
        &mut self,
        terms: &[Term],
        conjunction_outside: bool,
        negate: bool,
    ) -> Option<NodeId> {
        let mut groups = Vec::with_capacity(terms.len());
        for term in terms {
            let mut lits = Vec::with_capacity(term.len());
            for &(v, positive) in term {
                lits.push(self.literal(v, positive != negate)?);
            }
            groups.push(self.join(!conjunction_outside, &lits));
        }
        Some(self.join(conjunction_outside, &groups))
    }

    fn nnf(
        &mut self,
        n: &Nnf,
    ) -> Option<NodeId> {
        Some(match n {
            | Nnf::Const(b) => truth_literal(self.graph, *b),
            | Nnf::Lit(v, positive) => self.literal(*v, *positive)?,
            | Nnf::And(parts) | Nnf::Or(parts) => {
                let mut items = Vec::with_capacity(parts.len());
                for part in parts {
                    items.push(self.nnf(part)?);
                }
                self.join(matches!(n, Nnf::And(_)), &items)
            },
        })
    }
}

/// A prime implicant: `bits` on the variables not in `mask`.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Implicant {
    mask: u32,
    bits: u32,
}

impl Implicant {
    const fn covers(
        self,
        minterm: u32,
    ) -> bool {
        minterm & !self.mask == self.bits
    }

    const fn literals(self) -> u32 {
        (!self.mask).count_ones()
    }
}

/// Quine-McCluskey combination step on `n` variables.
fn prime_implicants(
    minterms: &[u32],
    n: usize,
) -> Vec<Implicant> {
    let mut current: BTreeSet<Implicant> = minterms
        .iter()
        .map(|&bits| Implicant { mask: 0, bits })
        .collect();
    let mut primes = Vec::new();
    while !current.is_empty() {
        let mut used = BTreeSet::new();
        let mut next = BTreeSet::new();
        for &imp in &current {
            for v in 0..n {
                let flag = 1_u32 << v;
                if imp.mask & flag != 0 || imp.bits & flag != 0 {
                    continue;
                }
                let partner = Implicant {
                    mask: imp.mask,
                    bits: imp.bits | flag,
                };
                if current.contains(&partner) {
                    used.insert(imp);
                    used.insert(partner);
                    next.insert(Implicant {
                        mask: imp.mask | flag,
                        bits: imp.bits,
                    });
                }
            }
        }
        primes.extend(current.iter().filter(|p| !used.contains(p)));
        current = next;
    }
    primes
}

/// Branch-and-bound state of the minimum cover search.
struct Cover<'a> {
    primes: &'a [Implicant],
    best: Vec<usize>,
    steps: usize,
}

impl Cover<'_> {
    fn cost(
        &self,
        chosen: &[usize],
    ) -> (usize, u32) {
        (
            chosen.len(),
            chosen
                .iter()
                .filter_map(|&c| self.primes.get(c))
                .map(|p| p.literals())
                .sum(),
        )
    }

    fn search(
        &mut self,
        uncovered: &[u32],
        chosen: &mut Vec<usize>,
    ) {
        self.steps = self.steps.saturating_add(1);
        let Some(&(pivot, _)) = uncovered
            .iter()
            .map(|&m| (m, self.primes.iter().filter(|p| p.covers(m)).count()))
            .collect::<Vec<_>>()
            .iter()
            .min_by_key(|(_, k)| *k)
        else {
            if self.cost(chosen) < self.cost(&self.best) {
                self.best = chosen.clone();
            }
            return;
        };
        let (count, lits) = self.cost(chosen);
        let (best_count, best_lits) = self.cost(&self.best);
        if count.saturating_add(1) > best_count
            || (count.saturating_add(1) == best_count && lits >= best_lits)
            || self.steps > MAX_SEARCH_STEPS
        {
            return;
        }
        let mut options: Vec<usize> = (0..self.primes.len())
            .filter(|&p| self.primes.get(p).is_some_and(|p| p.covers(pivot)))
            .collect();
        options.sort_by_key(|&p| {
            let prime = self.primes.get(p).copied();
            let gain = prime.map_or(0, |prime| {
                uncovered.iter().filter(|&&m| prime.covers(m)).count()
            });
            (
                std::cmp::Reverse(gain),
                prime.map_or(0, Implicant::literals),
                p,
            )
        });
        for p in options {
            let Some(&prime) = self.primes.get(p) else {
                continue;
            };
            let rest: Vec<u32> = uncovered
                .iter()
                .copied()
                .filter(|&m| !prime.covers(m))
                .collect();
            chosen.push(p);
            self.search(&rest, chosen);
            chosen.pop();
        }
    }
}

/// A minimal sum of products for the function with the given minterms.
///
/// Exact (minimum number of terms, then literals) unless the search bound
/// is hit, in which case the greedy cover found first is returned.
fn minimal_terms(
    minterms: &[u32],
    n: usize,
) -> Vec<Term> {
    let primes = prime_implicants(minterms, n);
    // Greedy cover as the initial bound.
    let mut greedy = Vec::new();
    let mut uncovered = minterms.to_vec();
    while !uncovered.is_empty() {
        let pick = (0..primes.len()).max_by_key(|&p| {
            let prime = primes.get(p).copied();
            let gain = prime.map_or(0, |prime| {
                uncovered.iter().filter(|&&m| prime.covers(m)).count()
            });
            (
                gain,
                std::cmp::Reverse(prime.map_or(0, Implicant::literals)),
                std::cmp::Reverse(p),
            )
        });
        let Some(prime) = pick.and_then(|p| primes.get(p).copied()) else {
            break;
        };
        uncovered.retain(|&m| !prime.covers(m));
        greedy.push(pick.unwrap_or(0));
    }
    let mut cover = Cover {
        primes: &primes,
        best: greedy,
        steps: 0,
    };
    cover.search(minterms, &mut Vec::new());
    let mut terms: Vec<Term> = cover
        .best
        .iter()
        .filter_map(|&p| primes.get(p))
        .map(|p| {
            (0..n)
                .filter(|&v| p.mask >> v & 1 == 0)
                .map(|v| (v, p.bits >> v & 1 == 1))
                .collect()
        })
        .collect();
    terms.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    terms
}

/// What a request asks for.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Form {
    Cnf,
    Dnf,
    Nnf,
    Minimal,
    Satisfiable,
    Tautology,
    Table,
}

/// Reduces normal-form and decision requests on propositional formulas.
struct Forms {
    ops: Ops,
    forms: Vec<(OpId, Form)>,
}

impl Forms {
    fn respond(
        graph: &mut Graph,
        ops: &Ops,
        form: Form,
        formula: &Formula,
        atoms: &[NodeId],
    ) -> Option<Outcome> {
        let n = atoms.len();
        let rows = |limit: usize| (n <= limit).then(|| 1_u32 << n);
        let mut budget = MAX_SIZE;
        let mut builder = Builder { graph, ops, atoms };
        match form {
            | Form::Nnf => {
                let nnf = formula.nnf(false, &mut budget)?;
                Some(Outcome::Pinned(builder.nnf(&nnf)?))
            },
            | Form::Dnf => {
                let terms = formula.nnf(false, &mut budget)?.dnf()?;
                Some(Outcome::Pinned(builder.normal_form(&terms, false, false)?))
            },
            // The clauses of the CNF are the negated terms of the DNF of the
            // negation.
            | Form::Cnf => {
                let terms = formula.nnf(true, &mut budget)?.dnf()?;
                Some(Outcome::Pinned(builder.normal_form(&terms, true, true)?))
            },
            | Form::Minimal => {
                let total = rows(MAX_MINIMISED)?;
                let minterms: Vec<u32> = (0..total).filter(|&m| formula.holds(m)).collect();
                let terms = minimal_terms(&minterms, n);
                Some(Outcome::Pinned(builder.normal_form(&terms, false, false)?))
            },
            | Form::Satisfiable => {
                let total = rows(MAX_ENUMERATED)?;
                let value = (0..total).any(|m| formula.holds(m));
                Some(Outcome::Equal(truth_literal(builder.graph, value)))
            },
            | Form::Tautology => {
                let total = rows(MAX_ENUMERATED)?;
                let value = (0..total).all(|m| formula.holds(m));
                Some(Outcome::Equal(truth_literal(builder.graph, value)))
            },
            | Form::Table => {
                let total = rows(MAX_MINIMISED)?;
                let mut table = Vec::with_capacity(1_usize << n);
                for row in 0..total {
                    // The first atom is the most significant digit.
                    let mask = (0..n).fold(0_u32, |m, v| {
                        m | (row >> n.saturating_sub(v).saturating_sub(1) & 1) << v
                    });
                    let mut cells: Vec<NodeId> = (0..n)
                        .map(|v| truth_literal(builder.graph, mask >> v & 1 == 1))
                        .collect();
                    cells.push(truth_literal(builder.graph, formula.holds(mask)));
                    table.push(builder.graph.node(core::LIST, &cells));
                }
                Some(Outcome::Pinned(builder.graph.node(core::LIST, &table)))
            },
        }
    }
}

impl Kernel for Forms {
    fn ops(&self) -> Vec<OpId> {
        self.forms.iter().map(|f| f.0).collect()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        let (Some(&(_, form)), &[arg]) =
            (self.forms.iter().find(|f| f.0 == op), graph.children(node))
        else {
            return Outcome::Pass;
        };
        let mut reader = Reader {
            graph,
            ops: &self.ops,
            classes: Vec::new(),
            nodes: Vec::new(),
        };
        let Some(mut formula) = reader.read(arg, 0) else {
            return Outcome::Pass;
        };
        let Reader { nodes, .. } = reader;
        // Atoms are ordered by name, so results do not depend on the order
        // the term happened to be traversed in.
        let mut order: Vec<usize> = (0..nodes.len()).collect();
        order.sort_by_key(|&k| nodes.get(k).map(|&n| graph.display(n)));
        let mut rank = vec![0; nodes.len()];
        for (new, &old) in order.iter().enumerate() {
            if let Some(slot) = rank.get_mut(old) {
                *slot = new;
            }
        }
        formula.remap(&rank);
        let sorted: Vec<NodeId> = order
            .iter()
            .filter_map(|&k| nodes.get(k).copied())
            .collect();
        Self::respond(graph, &self.ops, form, &formula, &sorted).unwrap_or(Outcome::Pass)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[logic()], src)
    }

    fn ev(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> f64 {
        eval(&[logic()], src, bindings)
    }

    /// Runs `check` on every assignment of `names` to 0/1 with `a` and `b`
    /// evaluated in one graph.
    fn for_all_assignments(
        a: &str,
        b: &str,
        names: &[&str],
        check: impl Fn(f64, f64, u32),
    ) {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[logic()]).is_ok());
        let ra = g.parse(a).unwrap_or_else(|e| panic!("{a}: {e}"));
        let rb = g.parse(b).unwrap_or_else(|e| panic!("{b}: {e}"));
        let symbols: Vec<_> = names.iter().map(|n| g.interner_mut().symbol(n)).collect();
        for mask in 0..1_u32 << names.len() {
            let mut env = Env::numeric(0.0);
            for (k, &sym) in symbols.iter().enumerate() {
                env.bind(sym, f64::from(mask >> k & 1));
            }
            let (x, y) = (g.eval(ra, &env), g.eval(rb, &env));
            check(x.unwrap_or(f64::NAN), y.unwrap_or(f64::NAN), mask);
        }
    }

    fn assert_equivalent(
        a: &str,
        b: &str,
        names: &[&str],
    ) {
        for_all_assignments(a, b, names, |x, y, mask| {
            assert!(
                (x - y).abs() < 1e-12,
                "`{a}` and `{b}` differ at {mask:b}: {x} vs {y}"
            );
        });
    }

    fn is_logical(name: &str) -> bool {
        matches!(name, "and" | "or" | "not" | "xor" | "implies" | "iff")
    }

    /// Whether `node` is a literal: an atom or the negation of one.
    fn is_literal(
        g: &Graph,
        node: NodeId,
    ) -> bool {
        let name = g.ops().get(g.op(node)).name.clone();
        match &*name {
            | "not" => {
                g.children(node)
                    .first()
                    .is_some_and(|&c| !is_logical(&g.ops().get(g.op(c)).name))
            },
            | n => !is_logical(n),
        }
    }

    /// Whether `text` is a normal form with `outer` at the root and `inner`
    /// below it, literals at the leaves.
    fn has_shape(
        text: &str,
        outer: &str,
        inner: &str,
    ) -> bool {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[logic()]).is_ok());
        let root = g.parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
        let name = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).name.to_string();
        let flat = |g: &Graph, n: NodeId| {
            is_literal(g, n)
                || (name(g, n) == inner && g.children(n).iter().all(|&c| is_literal(g, c)))
        };
        if name(&g, root) == outer {
            g.children(root).iter().all(|&c| flat(&g, c))
        } else {
            flat(&g, root)
        }
    }

    fn no_negated_compound(text: &str) -> bool {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[logic()]).is_ok());
        let root = g.parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
        let mut stack = vec![root];
        while let Some(n) = stack.pop() {
            let name = g.ops().get(g.op(n)).name.to_string();
            if matches!(name.as_str(), "implies" | "iff" | "xor") {
                return false;
            }
            if name == "not"
                && !g
                    .children(n)
                    .iter()
                    .all(|&c| !is_logical(&g.ops().get(g.op(c)).name))
            {
                return false;
            }
            stack.extend_from_slice(g.children(n));
        }
        true
    }

    const FORMULAS: [(&str, &[&str]); 9] = [
        ("implies(a, b)", &["a", "b"]),
        ("iff(a, b)", &["a", "b"]),
        ("xor(a, b, c)", &["a", "b", "c"]),
        ("not(and(a, or(b, not(c))))", &["a", "b", "c"]),
        ("or(and(a, b), and(not(a), c), and(b, c))", &["a", "b", "c"]),
        ("implies(iff(a, b), xor(c, d))", &["a", "b", "c", "d"]),
        (
            "and(implies(a, b), implies(b, c), implies(c, a))",
            &["a", "b", "c"],
        ),
        ("not(iff(xor(a, b), and(c, not(d))))", &["a", "b", "c", "d"]),
        ("or(a, not(a))", &["a"]),
    ];

    #[test]
    fn operators_on_literals() {
        let t = |src: &str| ev(src, &[]);
        assert!((t("and(true, true)") - 1.0).abs() < 1e-9);
        assert!((t("and(true, false)") - 0.0).abs() < 1e-9);
        assert!((t("and()") - 1.0).abs() < 1e-9);
        assert!((t("or(false, false)") - 0.0).abs() < 1e-9);
        assert!((t("or(false, true, false)") - 1.0).abs() < 1e-9);
        assert!((t("or()") - 0.0).abs() < 1e-9);
        assert!((t("xor(true, true)") - 0.0).abs() < 1e-9);
        assert!((t("xor(true, true, true)") - 1.0).abs() < 1e-9);
        assert!((t("not(true)") - 0.0).abs() < 1e-9);
        assert!((t("not(false)") - 1.0).abs() < 1e-9);
        assert!((t("implies(true, false)") - 0.0).abs() < 1e-9);
        assert!((t("implies(false, false)") - 1.0).abs() < 1e-9);
        assert!((t("implies(true, true)") - 1.0).abs() < 1e-9);
        assert!((t("iff(false, false)") - 1.0).abs() < 1e-9);
        assert!((t("iff(true, false)") - 0.0).abs() < 1e-9);
        assert!((t("true") - 1.0).abs() < 1e-9);
        assert!((t("false") - 0.0).abs() < 1e-9);
        // Any non-zero value counts as true, and results are exactly 0 or 1.
        assert!((ev("and(x, y)", &[("x", 2.5), ("y", -3.0)]) - 1.0).abs() < 1e-9);
        assert!((ev("not(x)", &[("x", 2.5)]) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn comparison_operators_evaluate() {
        let b = [("x", 1.0), ("y", 2.0)];
        assert!((ev("lt(x, y)", &b) - 1.0).abs() < 1e-9);
        assert!((ev("le(x, x)", &b) - 1.0).abs() < 1e-9);
        assert!((ev("gt(x, y)", &b) - 0.0).abs() < 1e-9);
        assert!((ev("ge(y, x)", &b) - 1.0).abs() < 1e-9);
        assert!((ev("ne(x, y)", &b) - 1.0).abs() < 1e-9);
        assert!((ev("ne(x, x)", &b) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn constants_are_the_boolean_literals() {
        assert_eq!(s("true"), "true");
        assert_eq!(s("and(true, not(false))"), "true");
        assert_eq!(s("or(false, false)"), "false");
        // The literal that number theory produces combines with the operators.
        let sets = [logic(), crate::rules::number_theory()];
        assert_eq!(simplify(&sets, "and(isprime(7), not(isprime(8)))"), "true");
        assert_eq!(simplify(&sets, "or(isprime(8), isprime(9))"), "false");
        assert_eq!(simplify(&sets, "implies(isprime(7), isprime(9))"), "false");
    }

    #[test]
    fn kernel_identities() {
        // Identity elements.
        assert_eq!(s("and(p, true)"), "p");
        assert_eq!(s("or(p, false)"), "p");
        assert_eq!(s("and(p, q, true)"), "and(p, q)");
        assert_eq!(s("or(p, false, q)"), "or(p, q)");
        assert_eq!(s("and(true, true)"), "true");
        assert_eq!(s("or(false, false)"), "false");
        // Idempotence.
        assert_eq!(s("and(p, p)"), "p");
        assert_eq!(s("or(p, q, p)"), "or(p, q)");
        // Absorption.
        assert_eq!(s("or(p, and(p, q))"), "p");
        assert_eq!(s("and(p, or(p, q))"), "p");
        assert_eq!(s("or(p, and(p, q, r), s)"), "or(p, s)");
        // Double negation and constants.
        assert_eq!(s("not(not(p))"), "p");
        assert_eq!(s("not(not(not(p)))"), "not(p)");
        assert_eq!(s("not(true)"), "false");
        assert_eq!(s("not(false)"), "true");
        // Exclusive or.
        assert_eq!(s("xor(p, false)"), "p");
        assert_eq!(s("xor(p, true)"), "not(p)");
        assert_eq!(s("xor(p, p)"), "false");
        assert_eq!(s("xor(p, q, p)"), "q");
        assert_eq!(s("xor(true, true)"), "false");
    }

    #[test]
    fn rewrites_fire() {
        assert_eq!(s("and(p, false)"), "false");
        assert_eq!(s("and(p, q, false)"), "false");
        assert_eq!(s("or(p, true)"), "true");
        assert_eq!(s("and(p, not(p))"), "false");
        assert_eq!(s("and(p, q, not(p))"), "false");
        assert_eq!(s("or(p, not(p))"), "true");
        assert_eq!(s("or(q, p, not(p))"), "true");
        assert_eq!(s("xor(p, p, q)"), "q");
        assert_eq!(s("implies(p, q)"), "or(q, not(p))");
        assert_eq!(s("implies(p, p)"), "true");
        assert_eq!(s("iff(p, p)"), "true");
        assert_eq!(s("iff(p, q) "), s("or(and(p, q), and(not(p), not(q)))"));
        assert_eq!(s("not(lt(x, y))"), "ge(x, y)");
        assert_eq!(s("not(le(x, y))"), "gt(x, y)");
        assert_eq!(s("not(gt(x, y))"), "le(x, y)");
        assert_eq!(s("not(ge(x, y))"), "lt(x, y)");
    }

    #[test]
    fn de_morgan_is_available_in_exploration() {
        // Both directions are structure-changing, so the answer is whichever
        // is smaller; check the two spellings are identified.
        let (a, _) = reduce_with(&[logic()], "not(and(p, q))", &[]);
        let (b, _) = reduce_with(&[logic()], "or(not(p), not(q))", &[]);
        assert_eq!(a, b);
        let (c, _) = reduce_with(&[logic()], "not(or(p, q))", &[]);
        let (d, _) = reduce_with(&[logic()], "and(not(p), not(q))", &[]);
        assert_eq!(c, d);
    }

    #[test]
    fn comparisons_of_literals() {
        assert_eq!(s("lt(1, 2)"), "true");
        assert_eq!(s("lt(2, 2)"), "false");
        assert_eq!(s("le(2, 2)"), "true");
        assert_eq!(s("gt(3, 2)"), "true");
        assert_eq!(s("ge(1, 2)"), "false");
        assert_eq!(s("ne(1, 2)"), "true");
        assert_eq!(s("ne(2, 2)"), "false");
        assert_eq!(s("eq(2, 2)"), "true");
        assert_eq!(s("2 = 3"), "false");
        assert_eq!(s("lt(1/3, 1/2)"), "true");
        assert_eq!(s("gt(-1/3, -1/2)"), "true");
        assert_eq!(s("lt(2^100, 2^100 + 1)"), "true");
        assert_eq!(s("and(lt(1, 2), gt(3, 2))"), "true");
        assert_eq!(s("lt(x, 2)"), "lt(x, 2)");
    }

    #[test]
    fn comparisons_from_assumptions() {
        let run = |src: &str, facts: Facts| reduce_with(&[logic()], src, &[("x", facts)]).0;
        assert_eq!(run("gt(x^2 + 1, 0)", Facts::REAL), "true");
        assert_eq!(run("gt(x, 0)", Facts::POSITIVE), "true");
        assert_eq!(run("lt(x, 0)", Facts::POSITIVE), "false");
        assert_eq!(run("ge(x, 0)", Facts::NONNEGATIVE), "true");
        assert_eq!(run("lt(0, x)", Facts::POSITIVE), "true");
        assert_eq!(run("gt(0, x)", Facts::POSITIVE), "false");
        assert_eq!(run("le(0, x)", Facts::NONNEGATIVE), "true");
        assert_eq!(run("lt(x, 0)", Facts::NEGATIVE), "true");
        assert_eq!(run("le(x, 0)", Facts::NEGATIVE), "true");
        assert_eq!(run("ge(x, 0)", Facts::NEGATIVE), "false");
        assert_eq!(run("ne(x, 0)", Facts::POSITIVE), "true");
        assert_eq!(run("eq(x, 0)", Facts::POSITIVE), "false");
        assert_eq!(run("gt(x, 0)", Facts::REAL), "gt(x, 0)");
        assert_eq!(run("gt(x, 1)", Facts::POSITIVE), "gt(x, 1)");
    }

    #[test]
    fn negation_normal_form() {
        assert_eq!(s("nnf(not(and(p, q)))"), "or(not(p), not(q))");
        assert_eq!(s("nnf(not(or(p, not(q))))"), "and(q, not(p))");
        assert_eq!(s("nnf(not(not(p)))"), "p");
        assert_eq!(s("nnf(implies(p, q))"), "or(q, not(p))");
        assert_eq!(s("nnf(not(implies(p, q)))"), "and(p, not(q))");
        for (f, names) in FORMULAS {
            let out = s(&format!("nnf({f})"));
            assert!(no_negated_compound(&out), "{f} -> {out}");
            assert_equivalent(f, &out, names);
        }
    }

    #[test]
    fn conjunctive_normal_form() {
        assert_eq!(s("cnf(or(and(p, q), r))"), "and(or(p, r), or(q, r))");
        assert_eq!(s("cnf(implies(p, q))"), "or(q, not(p))");
        assert_eq!(s("cnf(and(p, or(q, r)))"), "and(p, or(q, r))");
        assert_eq!(s("cnf(p)"), "p");
        assert_eq!(s("cnf(or(p, not(p)))"), "true");
        // A contradiction in CNF is its two unit clauses.
        assert_eq!(s("cnf(and(p, not(p)))"), "and(p, not(p))");
        for (f, names) in FORMULAS {
            let out = s(&format!("cnf({f})"));
            assert!(
                has_shape(&out, "and", "or") || out == "true" || out == "false",
                "{f} -> {out}"
            );
            assert_equivalent(f, &out, names);
        }
    }

    #[test]
    fn disjunctive_normal_form() {
        assert_eq!(s("dnf(and(or(p, q), r))"), "or(and(p, r), and(q, r))");
        assert_eq!(s("dnf(implies(p, q))"), "or(q, not(p))");
        assert_eq!(s("dnf(p)"), "p");
        assert_eq!(s("dnf(and(p, not(p)))"), "false");
        for (f, names) in FORMULAS {
            let out = s(&format!("dnf({f})"));
            assert!(
                has_shape(&out, "or", "and") || out == "true" || out == "false",
                "{f} -> {out}"
            );
            assert_equivalent(f, &out, names);
        }
    }

    #[test]
    fn minimal_sum_of_products() {
        assert_eq!(s("simplify_logic(or(and(p, q), and(p, not(q))))"), "p");
        assert_eq!(s("simplify_logic(and(p, or(p, q)))"), "p");
        assert_eq!(s("simplify_logic(or(p, not(p)))"), "true");
        assert_eq!(s("simplify_logic(and(p, not(p)))"), "false");
        assert_eq!(s("simplify_logic(implies(p, q))"), "or(q, not(p))");
        // Consensus: the middle term is redundant.
        assert_eq!(
            s("simplify_logic(or(and(a, b), and(not(a), c), and(b, c)))"),
            "or(and(a, b), and(c, not(a)))"
        );
        // Classic four-variable example with don't-care structure.
        let f = "or(and(not(a), not(b), not(c), not(d)), and(not(a), not(b), not(c), d), and(not(a), b, not(c), not(d)), \
                 and(a, not(b), not(c), not(d)), and(a, not(b), not(c), d), and(a, not(b), c, not(d)), and(a, b, not(c), not(d)))";
        let out = s(&format!("simplify_logic({f})"));
        assert!(has_shape(&out, "or", "and"), "{out}");
        assert_equivalent(f, &out, &["a", "b", "c", "d"]);
        // Minimal: 4 terms cover this function.
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[logic()]).is_ok());
        let root = g.parse(&out).unwrap_or_else(|e| panic!("{e}"));
        assert!(g.children(root).len() <= 4, "{out}");
        for (f, names) in FORMULAS {
            let out = s(&format!("simplify_logic({f})"));
            assert!(
                has_shape(&out, "or", "and") || out == "true" || out == "false",
                "{f} -> {out}"
            );
            assert_equivalent(f, &out, names);
        }
    }

    #[test]
    fn satisfiability() {
        assert_eq!(s("satisfiable(and(p, q))"), "true");
        assert_eq!(s("satisfiable(and(p, not(p)))"), "false");
        assert_eq!(s("satisfiable(true)"), "true");
        assert_eq!(s("satisfiable(false)"), "false");
        assert_eq!(s("satisfiable(and(or(p, q), not(p), not(q)))"), "false");
        assert_eq!(s("satisfiable(and(lt(x, 2), gt(x, 1)))"), "true");
        // Three pigeons, two holes: pigeon i sits in hole j iff `pij`.
        let pigeons = "and(or(p11, p12), or(p21, p22), or(p31, p32), \
             or(not(p11), not(p21)), or(not(p11), not(p31)), or(not(p21), not(p31)), \
             or(not(p12), not(p22)), or(not(p12), not(p32)), or(not(p22), not(p32)))";
        assert_eq!(s(&format!("satisfiable({pigeons})")), "false");
        // Two pigeons fit.
        let two = "and(or(p11, p12), or(p21, p22), or(not(p11), not(p21)), or(not(p12), not(p22)))";
        assert_eq!(s(&format!("satisfiable({two})")), "true");
    }

    #[test]
    fn tautologies() {
        assert_eq!(s("tautology(or(p, not(p)))"), "true");
        assert_eq!(s("tautology(p)"), "false");
        assert_eq!(s("tautology(implies(and(p, implies(p, q)), q))"), "true");
        assert_eq!(
            s("tautology(iff(not(and(p, q)), or(not(p), not(q))))"),
            "true"
        );
        assert_eq!(s("tautology(implies(p, q))"), "false");
        assert_eq!(s("tautology(true)"), "true");
    }

    #[test]
    fn truth_tables() {
        assert_eq!(
            s("truth_table(and(p, q))"),
            "list(list(false, false, false), list(false, true, false), list(true, false, false), list(true, true, true))"
        );
        assert_eq!(
            s("truth_table(not(p))"),
            "list(list(false, true), list(true, false))"
        );
        // Variables are ordered by name, not by appearance.
        assert_eq!(
            s("truth_table(implies(q, p))"),
            "list(list(false, false, true), list(false, true, false), list(true, false, true), list(true, true, true))"
        );
        assert_eq!(s("truth_table(true)"), "list(list(true))");
    }

    #[test]
    fn too_many_variables_stay_requests() {
        let vars: Vec<String> = (0..13).map(|k| format!("v{k:02}")).collect();
        let f = format!("or({})", vars.join(", "));
        let (text, reduced) = reduce_with(&[logic()], &format!("truth_table({f})"), &[]);
        assert!(!reduced, "{text}");
        let (text, reduced) = reduce_with(&[logic()], &format!("simplify_logic({f})"), &[]);
        assert!(!reduced, "{text}");
        // Normal forms and decisions still work for that size.
        assert_eq!(s(&format!("satisfiable({f})")), "true");
        assert_eq!(s(&format!("tautology({f})")), "false");
        let cnf = s(&format!("cnf({f})"));
        assert_eq!(cnf, f);
    }

    #[test]
    fn requests_evaluate_as_identity() {
        assert!((ev("cnf(x)", &[("x", 1.0)]) - 1.0).abs() < 1e-9);
        assert!((ev("simplify_logic(x)", &[("x", 0.0)]) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn atoms_may_be_comparisons() {
        assert_eq!(s("cnf(or(and(lt(x, 1), gt(y, 2)), lt(x, 1)))"), "lt(x, 1)");
        assert_eq!(
            s("simplify_logic(or(lt(x, 1), and(lt(x, 1), gt(y, 2))))"),
            "lt(x, 1)"
        );
        assert_eq!(s("tautology(or(lt(x, y), not(lt(x, y))))"), "true");
    }

    #[test]
    fn pinned_results_survive_extraction() {
        // The smaller spelling of the request's answer is not substituted for
        // the requested form.
        let sets = [logic()];
        let (text, reduced) = reduce_with(&sets, "dnf(and(or(p, q), or(r, t)))", &[]);
        assert!(reduced);
        assert_eq!(text, "or(and(p, r), and(p, t), and(q, r), and(q, t))");
        let _ = Budget::default();
    }
}

// Decision diagrams: model counting, model enumeration, quantifier
// elimination, entailment and equivalence on propositional formulas.
//
// Formulas are compiled into reduced ordered binary decision diagrams over
// their atoms (ordered by name), so the answers are canonical and scale
// with the size of the diagram rather than with `2^n`: up to 64 atoms and
// a million diagram nodes.
//
// | operator | value |
// |---|---|
// | `count_models(f)` | the number of assignments of the atoms of `f` that satisfy it |
// | `sat_model(f)` | a satisfying assignment as the conjunction of one literal per atom (the lexicographically least with `false < true`), or `false` when unsatisfiable |
// | `all_models(f)` | `list(cube, ...)` of every satisfying assignment over the atoms of `f` (at most 256) |
// | `bdd(f)` | the diagram as a nested Shannon expansion `or(and(x, hi), and(not(x), lo))` |
// | `bdd_size(f)` | the number of internal nodes of the diagram |
// | `exists(x, f)`, `forall(x, f)` | quantifier elimination over a Boolean atom or a `list` of atoms: `f[x:=0] or f[x:=1]` and `f[x:=0] and f[x:=1]` |
// | `restrict(f, x, v)` | the cofactor of `f` with the atom `x` fixed to the truth value `v` |
// | `entails(f, g)`, `equivalent(f, g)` | whether `f -> g` is a tautology, and whether `f <-> g` is |
// | `prime_implicants(f)` | `list(and(...), ...)`: the prime implicants (Quine–McCluskey, at most 12 atoms) |

use std::collections::HashMap;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;


/// Most atoms a diagram is built over.
const MAX_ATOMS: usize = 64;
/// Most nodes a diagram may grow to.
const MAX_NODES: usize = 1_000_000;
/// Most models `all_models` lists.
const MAX_MODELS: usize = 256;

const FALSE: usize = 0;
const TRUE: usize = 1;

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
enum Gate {
    And,
    Or,
    Xor,
    Implies,
    Iff,
}

impl Gate {
    const fn on(
        self,
        a: bool,
        b: bool,
    ) -> bool {
        match self {
            | Self::And => a && b,
            | Self::Or => a || b,
            | Self::Xor => a != b,
            | Self::Implies => !a || b,
            | Self::Iff => a == b,
        }
    }
}

/// A reduced ordered binary decision diagram.
struct Bdd {
    /// `(variable, low, high)`; entries 0 and 1 are the terminals.
    nodes: Vec<(usize, usize, usize)>,
    unique: HashMap<(usize, usize, usize), usize>,
    cache: HashMap<(Gate, usize, usize), usize>,
}

impl Bdd {
    fn new() -> Self {
        Self {
            nodes: vec![(usize::MAX, 0, 0), (usize::MAX, 1, 1)],
            unique: HashMap::new(),
            cache: HashMap::new(),
        }
    }

    fn variable(
        &self,
        node: usize,
    ) -> usize {
        self.nodes.get(node).map_or(usize::MAX, |n| n.0)
    }

    fn low(
        &self,
        node: usize,
    ) -> usize {
        self.nodes.get(node).map_or(0, |n| n.1)
    }

    fn high(
        &self,
        node: usize,
    ) -> usize {
        self.nodes.get(node).map_or(0, |n| n.2)
    }

    fn make(
        &mut self,
        variable: usize,
        low: usize,
        high: usize,
    ) -> Option<usize> {
        if low == high {
            return Some(low);
        }
        if let Some(&node) = self.unique.get(&(variable, low, high)) {
            return Some(node);
        }
        if self.nodes.len() >= MAX_NODES {
            return None;
        }
        self.nodes.push((variable, low, high));
        let id = self.nodes.len() - 1;
        self.unique.insert((variable, low, high), id);
        Some(id)
    }

    fn apply(
        &mut self,
        gate: Gate,
        a: usize,
        b: usize,
    ) -> Option<usize> {
        if a <= TRUE && b <= TRUE {
            return Some(usize::from(gate.on(a == TRUE, b == TRUE)));
        }
        if let Some(&hit) = self.cache.get(&(gate, a, b)) {
            return Some(hit);
        }
        let variable = self.variable(a).min(self.variable(b));
        let (a_low, a_high) = if self.variable(a) == variable {
            (self.low(a), self.high(a))
        } else {
            (a, a)
        };
        let (b_low, b_high) = if self.variable(b) == variable {
            (self.low(b), self.high(b))
        } else {
            (b, b)
        };
        let low = self.apply(gate, a_low, b_low)?;
        let high = self.apply(gate, a_high, b_high)?;
        let node = self.make(variable, low, high)?;
        self.cache.insert((gate, a, b), node);
        Some(node)
    }

    fn negate(
        &mut self,
        a: usize,
    ) -> Option<usize> {
        self.apply(Gate::Xor, a, TRUE)
    }

    fn build(
        &mut self,
        formula: &Formula,
    ) -> Option<usize> {
        match formula {
            | Formula::Const(b) => Some(usize::from(*b)),
            | Formula::Atom(v) => self.make(*v, FALSE, TRUE),
            | Formula::Not(x) => {
                let inner = self.build(x)?;
                self.negate(inner)
            },
            | Formula::And(xs) => self.fold(Gate::And, TRUE, xs),
            | Formula::Or(xs) => self.fold(Gate::Or, FALSE, xs),
            | Formula::Xor(xs) => self.fold(Gate::Xor, FALSE, xs),
            | Formula::Implies(a, b) => {
                let (a, b) = (self.build(a)?, self.build(b)?);
                self.apply(Gate::Implies, a, b)
            },
            | Formula::Iff(a, b) => {
                let (a, b) = (self.build(a)?, self.build(b)?);
                self.apply(Gate::Iff, a, b)
            },
        }
    }

    fn fold(
        &mut self,
        gate: Gate,
        unit: usize,
        items: &[Formula],
    ) -> Option<usize> {
        let mut acc = unit;
        for item in items {
            let next = self.build(item)?;
            acc = self.apply(gate, acc, next)?;
        }
        Some(acc)
    }

    fn restrict(
        &mut self,
        node: usize,
        variable: usize,
        value: bool,
        memo: &mut HashMap<usize, usize>,
    ) -> Option<usize> {
        let here = self.variable(node);
        if node <= TRUE || here > variable {
            return Some(node);
        }
        if let Some(&hit) = memo.get(&node) {
            return Some(hit);
        }
        let result = if here == variable {
            if value {
                self.high(node)
            } else {
                self.low(node)
            }
        } else {
            let (low, high) = (self.low(node), self.high(node));
            let low = self.restrict(low, variable, value, memo)?;
            let high = self.restrict(high, variable, value, memo)?;
            self.make(here, low, high)?
        };
        memo.insert(node, result);
        Some(result)
    }

    fn quantify(
        &mut self,
        node: usize,
        variable: usize,
        universal: bool,
    ) -> Option<usize> {
        let low = self.restrict(node, variable, false, &mut HashMap::new())?;
        let high = self.restrict(node, variable, true, &mut HashMap::new())?;
        self.apply(
            if universal {
                Gate::And
            } else {
                Gate::Or
            },
            low,
            high,
        )
    }

    /// Satisfying assignments over `n` variables.
    fn count(
        &self,
        root: usize,
        n: usize,
    ) -> BigInt {
        fn level(
            bdd: &Bdd,
            node: usize,
            n: usize,
        ) -> usize {
            if node <= TRUE {
                n
            } else {
                bdd.variable(node)
            }
        }
        fn go(
            bdd: &Bdd,
            node: usize,
            n: usize,
            memo: &mut HashMap<usize, BigInt>,
        ) -> BigInt {
            if node == FALSE {
                return BigInt::zero();
            }
            if node == TRUE {
                return BigInt::one();
            }
            if let Some(hit) = memo.get(&node) {
                return hit.clone();
            }
            let here = bdd.variable(node);
            let (low, high) = (bdd.low(node), bdd.high(node));
            let low_count = go(bdd, low, n, memo) << (level(bdd, low, n) - here - 1);
            let high_count = go(bdd, high, n, memo) << (level(bdd, high, n) - here - 1);
            let total = low_count + high_count;
            memo.insert(node, total.clone());
            total
        }
        go(self, root, n, &mut HashMap::new()) << level(self, root, n)
    }

    /// The least satisfying assignment (`false < true`) over `n` variables.
    fn first_model(
        &self,
        root: usize,
        n: usize,
    ) -> Option<Vec<bool>> {
        if root == FALSE {
            return None;
        }
        let mut model = vec![false; n];
        let mut node = root;
        while node > TRUE {
            let (low, high) = (self.low(node), self.high(node));
            let take_high = low == FALSE;
            if let Some(slot) = model.get_mut(self.variable(node)) {
                *slot = take_high;
            }
            node = if take_high { high } else { low };
        }
        Some(model)
    }

    /// Every satisfying assignment over `n` variables, up to `limit`.
    fn models(
        &self,
        root: usize,
        n: usize,
        limit: usize,
    ) -> Option<Vec<Vec<bool>>> {
        let mut out = Vec::new();
        let mut partial = vec![false; n];
        self.walk(root, 0, n, &mut partial, &mut out, limit)?;
        Some(out)
    }

    fn walk(
        &self,
        node: usize,
        level: usize,
        n: usize,
        partial: &mut Vec<bool>,
        out: &mut Vec<Vec<bool>>,
        limit: usize,
    ) -> Option<()> {
        if node == FALSE {
            return Some(());
        }
        if level == n {
            if out.len() >= limit {
                return None;
            }
            out.push(partial.clone());
            return Some(());
        }
        let here = if node == TRUE {
            usize::MAX
        } else {
            self.variable(node)
        };
        for value in [false, true] {
            if let Some(slot) = partial.get_mut(level) {
                *slot = value;
            }
            let next = if here == level {
                if value {
                    self.high(node)
                } else {
                    self.low(node)
                }
            } else {
                node
            };
            self.walk(next, level + 1, n, partial, out, limit)?;
        }
        Some(())
    }

    fn size(
        &self,
        root: usize,
    ) -> usize {
        let mut seen = BTreeSet::new();
        let mut stack = vec![root];
        while let Some(node) = stack.pop() {
            if node > TRUE && seen.insert(node) {
                stack.push(self.low(node));
                stack.push(self.high(node));
            }
        }
        seen.len()
    }

    fn term(
        &self,
        builder: &mut Builder<'_>,
        node: usize,
        memo: &mut HashMap<usize, NodeId>,
    ) -> Option<NodeId> {
        if node <= TRUE {
            return Some(truth_literal(builder.graph, node == TRUE));
        }
        if let Some(&hit) = memo.get(&node) {
            return Some(hit);
        }
        let (low, high) = (self.low(node), self.high(node));
        let positive = builder.literal(self.variable(node), true)?;
        let negative = builder.literal(self.variable(node), false)?;
        let built = match (low, high) {
            | (FALSE, high) => {
                let high = self.term(builder, high, memo)?;
                if high_is_true(high, builder) {
                    positive
                } else {
                    builder.join(true, &[positive, high])
                }
            },
            | (low, FALSE) => {
                let low = self.term(builder, low, memo)?;
                if high_is_true(low, builder) {
                    negative
                } else {
                    builder.join(true, &[negative, low])
                }
            },
            | (TRUE, high) => {
                let high = self.term(builder, high, memo)?;
                builder.join(false, &[negative, high])
            },
            | (low, TRUE) => {
                let low = self.term(builder, low, memo)?;
                builder.join(false, &[positive, low])
            },
            | (low, high) => {
                let (low, high) = (
                    self.term(builder, low, memo)?,
                    self.term(builder, high, memo)?,
                );
                let then = builder.join(true, &[positive, high]);
                let otherwise = builder.join(true, &[negative, low]);
                builder.join(false, &[then, otherwise])
            },
        };
        memo.insert(node, built);
        Some(built)
    }
}

/// Whether a built term is the literal `true`.
fn high_is_true(
    node: NodeId,
    builder: &Builder<'_>,
) -> bool {
    constant(builder.graph, builder.ops, node) == Some(true)
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum Kind {
    Count,
    Model,
    AllModels,
    Diagram,
    Size,
    Exists,
    Forall,
    Restrict,
    Entails,
    Equivalent,
    Primes,
}

struct Extension {
    ops: Ops,
    kinds: Vec<(OpId, Kind)>,
}

/// One compiled request.
struct Problem {
    formulas: Vec<Formula>,
    atoms: Vec<NodeId>,
}

impl Extension {
    /// Reads the arguments as formulas over a shared, name-ordered set of
    /// atoms.
    fn read(
        &self,
        graph: &Graph,
        args: &[NodeId],
    ) -> Option<Problem> {
        let mut reader = Reader {
            graph,
            ops: &self.ops,
            classes: Vec::new(),
            nodes: Vec::new(),
        };
        let mut formulas = Vec::with_capacity(args.len());
        for &arg in args {
            formulas.push(reader.read(arg, 0)?);
        }
        let Reader { nodes, .. } = reader;
        if nodes.len() > MAX_ATOMS {
            return None;
        }
        let mut order: Vec<usize> = (0..nodes.len()).collect();
        order.sort_by_key(|&k| nodes.get(k).map(|&n| graph.display(n)));
        let mut rank = vec![0; nodes.len()];
        for (new, &old) in order.iter().enumerate() {
            if let Some(slot) = rank.get_mut(old) {
                *slot = new;
            }
        }
        for formula in &mut formulas {
            formula.remap(&rank);
        }
        let atoms = order
            .iter()
            .filter_map(|&k| nodes.get(k).copied())
            .collect();
        Some(Problem { formulas, atoms })
    }

    #[allow(clippy::too_many_lines)]
    fn answer(
        &self,
        graph: &mut Graph,
        kind: Kind,
        args: &[NodeId],
    ) -> Option<Outcome> {
        // A quantifier prefix is a single atom or a list of atoms.
        let flat: Vec<NodeId> = match (kind, args) {
            | (Kind::Exists | Kind::Forall, &[prefix, body]) => {
                let mut items = if graph.op(prefix) == core::LIST {
                    graph.children(prefix).to_vec()
                } else {
                    vec![prefix]
                };
                items.push(body);
                items
            },
            | _ => args.to_vec(),
        };
        // Real polynomial formulas are decided by real quantifier
        // elimination; only genuinely Boolean bodies reach the
        // decision-diagram path below.
        if matches!(kind, Kind::Exists | Kind::Forall) {
            let &[prefix, real_body] = args else {
                return None;
            };
            match crate::rules::solve::qe::try_quantifier(
                graph,
                prefix,
                real_body,
                kind == Kind::Forall,
            ) {
                | crate::rules::solve::qe::Routed::NotReal => {},
                | crate::rules::solve::qe::Routed::Result(term) => {
                    return Some(Outcome::Equal(term));
                },
                | crate::rules::solve::qe::Routed::Unevaluated => return None,
            }
        }
        let Problem { formulas, atoms } = self.read(graph, &flat)?;
        let n = atoms.len();
        let mut bdd = Bdd::new();
        let ops = &self.ops;
        let cube = |graph: &mut Graph, model: &[bool]| -> Option<NodeId> {
            let mut builder = Builder {
                graph,
                ops,
                atoms: &atoms,
            };
            let mut lits = Vec::with_capacity(model.len());
            for (v, &value) in model.iter().enumerate() {
                lits.push(builder.literal(v, value)?);
            }
            Some(builder.join(true, &lits))
        };
        match kind {
            | Kind::Count => {
                let root = bdd.build(formulas.first()?)?;
                let count = bdd.count(root, n);
                Some(Outcome::Equal(graph.num(crate::graph::Number::Int(count))))
            },
            | Kind::Size => {
                let root = bdd.build(formulas.first()?)?;
                let size = i64::try_from(bdd.size(root)).ok()?;
                Some(Outcome::Equal(graph.int(size)))
            },
            | Kind::Model => {
                let root = bdd.build(formulas.first()?)?;
                match bdd.first_model(root, n) {
                    | None => Some(Outcome::Equal(truth_literal(graph, false))),
                    | Some(model) => Some(Outcome::Pinned(cube(graph, &model)?)),
                }
            },
            | Kind::AllModels => {
                let root = bdd.build(formulas.first()?)?;
                let models = bdd.models(root, n, MAX_MODELS)?;
                let mut items = Vec::with_capacity(models.len());
                for model in &models {
                    items.push(cube(graph, model)?);
                }
                Some(Outcome::Pinned(graph.node(core::LIST, &items)))
            },
            | Kind::Diagram => {
                let root = bdd.build(formulas.first()?)?;
                let mut builder = Builder {
                    graph,
                    ops,
                    atoms: &atoms,
                };
                let term = bdd.term(&mut builder, root, &mut HashMap::new())?;
                Some(Outcome::Pinned(term))
            },
            | Kind::Exists | Kind::Forall => {
                let (body, quantified) = formulas.split_last()?;
                let mut root = bdd.build(body)?;
                for formula in quantified {
                    let Formula::Atom(variable) = formula else {
                        return None;
                    };
                    root = bdd.quantify(root, *variable, kind == Kind::Forall)?;
                }
                let mut builder = Builder {
                    graph,
                    ops,
                    atoms: &atoms,
                };
                Some(Outcome::Pinned(bdd.term(
                    &mut builder,
                    root,
                    &mut HashMap::new(),
                )?))
            },
            | Kind::Restrict => {
                let [body, Formula::Atom(variable), Formula::Const(value)] = formulas.as_slice()
                else {
                    return None;
                };
                let root = bdd.build(body)?;
                let restricted = bdd.restrict(root, *variable, *value, &mut HashMap::new())?;
                let mut builder = Builder {
                    graph,
                    ops,
                    atoms: &atoms,
                };
                Some(Outcome::Pinned(bdd.term(
                    &mut builder,
                    restricted,
                    &mut HashMap::new(),
                )?))
            },
            | Kind::Entails | Kind::Equivalent => {
                let [a, b] = formulas.as_slice() else {
                    return None;
                };
                let (a, b) = (bdd.build(a)?, bdd.build(b)?);
                let holds = if kind == Kind::Equivalent {
                    a == b
                } else {
                    bdd.apply(Gate::Implies, a, b)? == TRUE
                };
                Some(Outcome::Equal(truth_literal(graph, holds)))
            },
            | Kind::Primes => {
                let formula = formulas.first()?;
                if n > MAX_MINIMISED {
                    return None;
                }
                let minterms: Vec<u32> = (0..1_u32 << n).filter(|&m| formula.holds(m)).collect();
                let mut builder = Builder {
                    graph,
                    ops,
                    atoms: &atoms,
                };
                let mut items = Vec::new();
                for implicant in prime_implicants(&minterms, n) {
                    let mut lits = Vec::new();
                    for v in 0..n {
                        if implicant.mask >> v & 1 == 0 {
                            lits.push(builder.literal(v, implicant.bits >> v & 1 == 1)?);
                        }
                    }
                    items.push(builder.join(true, &lits));
                }
                Some(Outcome::Pinned(graph.node(core::LIST, &items)))
            },
        }
    }
}

impl Kernel for Extension {
    fn ops(&self) -> Vec<OpId> {
        self.kinds.iter().map(|k| k.0).collect()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        let Some(&(_, kind)) = self.kinds.iter().find(|k| k.0 == op) else {
            return Outcome::Pass;
        };
        let args = graph.children(node).to_vec();
        self.answer(graph, kind, &args).unwrap_or(Outcome::Pass)
    }
}

fn install_ext(
    i: &mut Installer<'_>,
    ops: Ops,
) -> Result<(), RuleError> {
    let mut kinds = Vec::new();
    for (name, arity, kind) in [
        ("count_models", 1, Kind::Count),
        ("sat_model", 1, Kind::Model),
        ("all_models", 1, Kind::AllModels),
        ("bdd", 1, Kind::Diagram),
        ("bdd_size", 1, Kind::Size),
        ("exists", 2, Kind::Exists),
        ("forall", 2, Kind::Forall),
        ("restrict", 3, Kind::Restrict),
        ("entails", 2, Kind::Entails),
        ("equivalent", 2, Kind::Equivalent),
        ("prime_implicants", 1, Kind::Primes),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        // Counting is a pure number; every other extension returns a formula,
        // a model or a list of them.
        let rule = match kind {
            | Kind::Count | Kind::Size => DimRule::Custom(dimensionless_value),
            | _ => DimRule::NonMetric,
        };
        i.dim(op, rule);
        kinds.push((op, kind));
    }
    i.kernel("logic/diagrams", Tier::Reduce, Extension { ops, kinds });
    Ok(())
}

#[cfg(test)]
mod ext_tests {
    use crate::rules::logic::logic;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[logic()], src)
    }

    #[test]
    fn counting_and_models() {
        assert_eq!(s("count_models(or(a, b))"), "3");
        assert_eq!(s("count_models(xor(a, b, c))"), "4");
        assert_eq!(s("count_models(and(a, not(a)))"), "0");
        let parity: Vec<String> = (0..30).map(|k| format!("v{k}")).collect();
        assert_eq!(
            s(&format!("count_models(xor({}))", parity.join(", "))),
            "536870912"
        );
        assert_eq!(s("sat_model(and(a, not(b)))"), "and(a, not(b))");
        assert_eq!(s("sat_model(or(a, b))"), "and(b, not(a))");
        assert_eq!(s("sat_model(and(a, not(a)))"), "false");
        assert_eq!(
            s("all_models(iff(a, b))"),
            "list(and(not(a), not(b)), and(a, b))"
        );
        assert_eq!(s("bdd_size(xor(a, b, c))"), "5");
        assert_eq!(s("bdd(and(a, b))"), "and(a, b)");
    }

    #[test]
    fn quantifier_elimination_and_restriction() {
        assert_eq!(s("exists(y, and(x, y))"), "x");
        assert_eq!(s("forall(x, or(x, y))"), "y");
        assert_eq!(s("forall(y, and(x, y))"), "false");
        assert_eq!(s("exists(list(x, y), and(x, y))"), "true");
        assert_eq!(s("exists(x, and(x, not(x)))"), "false");
        assert_eq!(s("restrict(or(and(x, y), z), x, true)"), "or(y, z)");
        assert_eq!(s("restrict(or(and(x, y), z), x, false)"), "z");
    }

    #[test]
    fn entailment_equivalence_and_implicants() {
        assert_eq!(s("entails(and(a, b), a)"), "true");
        assert_eq!(s("entails(a, and(a, b))"), "false");
        assert_eq!(s("equivalent(implies(a, b), or(not(a), b))"), "true");
        assert_eq!(s("equivalent(xor(a, b), iff(a, b))"), "false");
        assert_eq!(
            s("prime_implicants(or(and(a, b), and(not(a), c)))"),
            "list(and(b, c), and(c, not(a)), and(a, b))"
        );
    }
}
