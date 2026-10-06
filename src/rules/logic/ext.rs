//! Decision diagrams: model counting, model enumeration, quantifier
//! elimination, entailment and equivalence on propositional formulas.
//!
//! Formulas are compiled into reduced ordered binary decision diagrams over
//! their atoms (ordered by name), so the answers are canonical and scale
//! with the size of the diagram rather than with `2^n`: up to 64 atoms and
//! a million diagram nodes.
//!
//! | operator | value |
//! |---|---|
//! | `count_models(f)` | the number of assignments of the atoms of `f` that satisfy it |
//! | `sat_model(f)` | a satisfying assignment as the conjunction of one literal per atom (the lexicographically least with `false < true`), or `false` when unsatisfiable |
//! | `all_models(f)` | `list(cube, ...)` of every satisfying assignment over the atoms of `f` (at most 256) |
//! | `bdd(f)` | the diagram as a nested Shannon expansion `or(and(x, hi), and(not(x), lo))` |
//! | `bdd_size(f)` | the number of internal nodes of the diagram |
//! | `exists(x, f)`, `forall(x, f)` | quantifier elimination over a Boolean atom or a `list` of atoms: `f[x:=0] or f[x:=1]` and `f[x:=0] and f[x:=1]` |
//! | `restrict(f, x, v)` | the cofactor of `f` with the atom `x` fixed to the truth value `v` |
//! | `entails(f, g)`, `equivalent(f, g)` | whether `f -> g` is a tautology, and whether `f <-> g` is |
//! | `prime_implicants(f)` | `list(and(...), ...)`: the prime implicants (Quine–McCluskey, at most 12 atoms) |

use std::collections::BTreeSet;
use std::collections::HashMap;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use super::prime_implicants;
use super::truth_literal;
use super::Builder;
use super::Formula;
use super::Ops;
use super::Reader;
use super::MAX_MINIMISED;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;

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
        Self { nodes: vec![(usize::MAX, 0, 0), (usize::MAX, 1, 1)], unique: HashMap::new(), cache: HashMap::new() }
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
        let (a_low, a_high) = if self.variable(a) == variable { (self.low(a), self.high(a)) } else { (a, a) };
        let (b_low, b_high) = if self.variable(b) == variable { (self.low(b), self.high(b)) } else { (b, b) };
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
            if value { self.high(node) } else { self.low(node) }
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
        self.apply(if universal { Gate::And } else { Gate::Or }, low, high)
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
            if node <= TRUE { n } else { bdd.variable(node) }
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
        let here = if node == TRUE { usize::MAX } else { self.variable(node) };
        for value in [false, true] {
            if let Some(slot) = partial.get_mut(level) {
                *slot = value;
            }
            let next = if here == level {
                if value { self.high(node) } else { self.low(node) }
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
                if high_is_true(high, builder) { positive } else { builder.join(true, &[positive, high]) }
            },
            | (low, FALSE) => {
                let low = self.term(builder, low, memo)?;
                if high_is_true(low, builder) { negative } else { builder.join(true, &[negative, low]) }
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
                let (low, high) = (self.term(builder, low, memo)?, self.term(builder, high, memo)?);
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
    super::constant(builder.graph, builder.ops, node) == Some(true)
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
        graph: &mut Graph,
        args: &[NodeId],
    ) -> Option<Problem> {
        let mut reader = Reader { graph, ops: &self.ops, classes: Vec::new(), nodes: Vec::new() };
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
        let atoms = order.iter().filter_map(|&k| nodes.get(k).copied()).collect();
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
                let mut items = if graph.op(prefix) == core::LIST { graph.children(prefix).to_vec() } else { vec![prefix] };
                items.push(body);
                items
            },
            | _ => args.to_vec(),
        };
        let Problem { formulas, atoms } = self.read(graph, &flat)?;
        let n = atoms.len();
        let mut bdd = Bdd::new();
        let ops = &self.ops;
        let cube = |graph: &mut Graph, model: &[bool]| -> Option<NodeId> {
            let mut builder = Builder { graph, ops, atoms: &atoms };
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
                let mut builder = Builder { graph, ops, atoms: &atoms };
                let term = bdd.term(&mut builder, root, &mut HashMap::new())?;
                Some(Outcome::Pinned(term))
            },
            | Kind::Exists | Kind::Forall => {
                let (body, quantified) = formulas.split_last()?;
                let mut root = bdd.build(body)?;
                for formula in quantified {
                    let Formula::Atom(variable) = formula else { return None };
                    root = bdd.quantify(root, *variable, kind == Kind::Forall)?;
                }
                let mut builder = Builder { graph, ops, atoms: &atoms };
                Some(Outcome::Pinned(bdd.term(&mut builder, root, &mut HashMap::new())?))
            },
            | Kind::Restrict => {
                let [body, Formula::Atom(variable), Formula::Const(value)] = formulas.as_slice() else {
                    return None;
                };
                let root = bdd.build(body)?;
                let restricted = bdd.restrict(root, *variable, *value, &mut HashMap::new())?;
                let mut builder = Builder { graph, ops, atoms: &atoms };
                Some(Outcome::Pinned(bdd.term(&mut builder, restricted, &mut HashMap::new())?))
            },
            | Kind::Entails | Kind::Equivalent => {
                let [a, b] = formulas.as_slice() else { return None };
                let (a, b) = (bdd.build(a)?, bdd.build(b)?);
                let holds = if kind == Kind::Equivalent { a == b } else { bdd.apply(Gate::Implies, a, b)? == TRUE };
                Some(Outcome::Equal(truth_literal(graph, holds)))
            },
            | Kind::Primes => {
                let formula = formulas.first()?;
                if n > MAX_MINIMISED {
                    return None;
                }
                let minterms: Vec<u32> = (0..1_u32 << n).filter(|&m| formula.holds(m)).collect();
                let mut builder = Builder { graph, ops, atoms: &atoms };
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

pub(super) fn install(
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
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        kinds.push((op, kind));
    }
    i.kernel("logic/diagrams", Tier::Reduce, Extension { ops, kinds });
    Ok(())
}

#[cfg(test)]
mod tests {
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
        assert_eq!(s(&format!("count_models(xor({}))", parity.join(", "))), "536870912");
        assert_eq!(s("sat_model(and(a, not(b)))"), "and(a, not(b))");
        assert_eq!(s("sat_model(or(a, b))"), "and(not(a), b)");
        assert_eq!(s("sat_model(and(a, not(a)))"), "false");
        assert_eq!(s("all_models(iff(a, b))"), "list(and(not(a), not(b)), and(a, b))");
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
        assert_eq!(s("prime_implicants(or(and(a, b), and(not(a), c)))"), "list(and(b, c), and(not(a), c), and(a, b))");
    }
}
