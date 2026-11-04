//! Extraction: choosing one concrete term out of an equivalence class.
//!
//! A [`CostModel`] prices individual e-nodes; the [`Extractor`] finds, for
//! every class reachable from the roots, the member whose total cost is
//! minimal and rebuilds that term as a concrete DAG node. Because the result
//! is hash-consed like everything else, shared subterms of the answer are
//! shared in memory.
//!
//! A model may return `None` to *forbid* an e-node. That is how phases are
//! expressed: [`ClosedForm`] forbids heavy operators, so a class only has a
//! closed-form extraction once some rule has actually reduced the request.

use std::cell::RefCell;
use std::collections::HashMap;
use std::collections::HashSet;

use super::id::ClassId;
use super::id::NodeId;
use super::id::SymbolId;
use super::op::OpFlags;
use super::op::core;
use super::store::Graph;
use crate::backend::Intrinsic;
use crate::backend::Lowering;
use crate::backend::tape::Inst;

/// Prices e-nodes for extraction.
pub trait CostModel {
    /// Cost of `node` itself, excluding its children. `None` forbids it.
    fn cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64>;

    /// Factor applied to the summed cost of `node`'s children.
    ///
    /// The default of `1` gives plain additive term size. A larger factor
    /// on an operator makes extraction prefer terms in which that operator
    /// is applied to *smaller* arguments.
    fn scale(
        &self,
        _graph: &Graph,
        _node: NodeId,
    ) -> u64 {
        1
    }
}

/// Term size weighted by the operators' declared costs. Accepts everything.
///
/// The arguments of heavy operators count eightfold, so among terms that
/// are not fully reduced the one where the unevaluated requests have been
/// pushed furthest towards the leaves wins: `f + x*diff(f, x)` is preferred
/// over `diff(x*f, x)` although it is longer.
#[derive(Copy, Clone, Debug, Default)]
pub struct SizeCost;

impl CostModel for SizeCost {
    fn cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        Some(base_cost(graph, node))
    }

    fn scale(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> u64 {
        if graph.ops().get(graph.op(node)).flags.has(OpFlags::HEAVY) {
            8
        } else {
            1
        }
    }
}

/// Penalty for arithmetic applied to a list: `2 * list(a, b)` equals
/// `list(2*a, 2*b)`, and a vector should read as one.
const LIST_ARITHMETIC: u64 = 1_000;

/// Cost of a constant or a symbol load.
const LEAF: u64 = 1;
/// Cost of an add, subtract, multiply, negate, compare or select.
const ARITHMETIC: u64 = 1;
/// Cost of a square root or cube root (a hardware instruction).
const ROOT: u64 = 4;
/// Cost of a division or an integer power (square-and-multiply).
const DIVISION: u64 = 4;
/// Cost of an integer power.
const POW_INT: u64 = 4;
/// Cost of a real power (`powf`).
const POW_REAL: u64 = 6;
/// Cost of a transcendental function (`exp`, `ln`, trigonometry, ...).
const TRANSCENDENTAL: u64 = 10;
/// Cost of `atan2` (two reductions and a division).
const ATAN2: u64 = 10;
/// Cost of a scalar call without a dedicated instruction.
const OPAQUE_CALL: u64 = 8;
/// Cost of an external C function: an ABI call, opaque to the count.
const EXTERN: u64 = 20;
/// Cost of a tensor gather (index arithmetic plus a load).
const GATHER: u64 = 3;
/// Cost of an unevaluated heavy request (`diff`, `integrate`, ...). Such a
/// term cannot be evaluated numerically at all, so it is priced far above any
/// closed form to keep extraction moving towards a reducible spelling.
const HEAVY: u64 = 100;

/// Evaluation-cost model: prices an e-node by the tape instruction its
/// operator lowers to, and a whole class by the weighted instruction count of
/// its cheapest representative.
///
/// Unlike [`SizeCost`], which counts term shape, this estimates the work a
/// numerical evaluation would do. The weights mirror the tape lowering:
///
/// | instruction | weight |
/// |---|---|
/// | constant, input, parameter | 1 |
/// | add, subtract, multiply, negate, compare, select | 1 |
/// | division, integer power | 4 |
/// | real power | 6 |
/// | `sqrt`, `cbrt` | 4 |
/// | `abs`, `floor`, `ceil` | 1 |
/// | other elementary functions, `atan2` | 10 |
/// | scalar call through the interpreter | 8 |
/// | external C function | 20 |
/// | tensor gather | 3 |
/// | heavy/unevaluatable operator | 100 |
///
/// The [`CostModel`] price of one e-node is the weight of its own operator;
/// the extractor sums these over the term, which approximates
/// [`EvalCost::term_cost`]. That method is the exact version — it lowers the
/// concrete term to the tape IR and sums the weighted instructions — and a
/// class is priced by [`EvalCost::class_cost`] from its cheapest member.
///
/// Both caches assume the graph does not change between calls; create a fresh
/// model (or call [`EvalCost::clear`]) for a graph that was unioned or
/// rewritten since.
#[derive(Debug, Default)]
pub struct EvalCost {
    /// Weighted tape cost of a concrete term, by node.
    terms: RefCell<HashMap<NodeId, Option<u64>>>,
    /// Price of a class, by class (lowest over its members).
    classes: RefCell<HashMap<ClassId, Option<u64>>>,
    /// Own-node instruction price used by the [`CostModel`] implementation.
    nodes: RefCell<HashMap<NodeId, u64>>,
}

impl EvalCost {
    /// A model with empty caches.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Weighted instruction count of evaluating the concrete term `node`.
    ///
    /// The term is lowered to the tape IR with its free symbols as
    /// parameters, optimised (constant folding and dead-code elimination), and
    /// every remaining instruction is priced by its class. `None` when the
    /// term cannot be lowered (an operator without scalar semantics). The
    /// result is cached per node.
    #[must_use]
    pub fn term_cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        if let Some(cached) = self.terms.borrow().get(&node).copied() {
            return cached;
        }
        let symbols = symbols_of(graph, node);
        let computed = crate::backend::tape::lower(graph, &[node], &[], &symbols)
            .ok()
            .map(|tape| tape.optimise())
            .map(|tape| {
                tape.insts
                    .iter()
                    .map(instruction_cost)
                    .fold(0_u64, u64::saturating_add)
            });
        self.terms.borrow_mut().insert(node, computed);
        computed
    }

    /// Weighted instruction count of the cheapest representative of `node`'s
    /// class.
    ///
    /// The result is cached per class, so repeated extraction over an
    /// unchanged graph pays for the lowering once.
    #[must_use]
    pub fn class_cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        let class = graph.find(node);
        if let Some(cached) = self.classes.borrow().get(&class).copied() {
            return cached;
        }
        let computed = graph
            .enodes(class)
            .filter_map(|member| self.term_cost(graph, member))
            .min();
        self.classes.borrow_mut().insert(class, computed);
        computed
    }

    /// Drops both caches, for use after the graph changed.
    pub fn clear(&mut self) {
        self.terms.get_mut().clear();
        self.classes.get_mut().clear();
        self.nodes.get_mut().clear();
    }

    /// Number of concrete terms in the per-node cache (diagnostics).
    #[must_use]
    pub fn cached_terms(&self) -> usize {
        self.terms.borrow().len()
    }

    /// Number of classes in the per-class cache (diagnostics).
    #[must_use]
    pub fn cached_classes(&self) -> usize {
        self.classes.borrow().len()
    }

    /// Number of e-nodes whose own price is cached (diagnostics).
    #[must_use]
    pub fn cached_nodes(&self) -> usize {
        self.nodes.borrow().len()
    }
}

impl CostModel for EvalCost {
    fn cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        if let Some(&cached) = self.nodes.borrow().get(&node) {
            return Some(cached);
        }
        let value = own_cost(graph, node);
        self.nodes.borrow_mut().insert(node, value);
        Some(value)
    }
}

/// Weight of the instruction a tape interpreter would run for `node` itself.
fn own_cost(
    graph: &Graph,
    node: NodeId,
) -> u64 {
    let op = graph.op(node);
    if op == core::LIT || op == core::SYM {
        return LEAF;
    }
    if op == core::ADD || op == core::MUL || op == core::EQ {
        return ARITHMETIC;
    }
    if op == core::POW {
        // An exact literal exponent lowers to square-and-multiply (or a root
        // followed by it); a general exponent is a `powf`.
        let exact_exponent = graph
            .children(node)
            .get(1)
            .and_then(|&e| graph.number_of(e))
            .is_some_and(|n| n.to_rational().is_some());
        return if exact_exponent {
            POW_INT
        } else {
            POW_REAL
        };
    }
    let desc = graph.ops().get(op);
    if desc.flags.has(OpFlags::HEAVY) {
        return HEAVY;
    }
    if let Some(lowering) = graph.ops().attr::<Lowering>(op).copied().or(desc.jit) {
        return match lowering {
            | Lowering::Intrinsic(kind) => intrinsic_cost(kind),
            | Lowering::Extern1(_) | Lowering::Extern2(_) => EXTERN,
            | Lowering::Gather => GATHER,
        };
    }
    if let Some(cost) = named_cost(&desc.name) {
        return cost;
    }
    if desc.eval.is_some() {
        // A scalar semantics with no dedicated instruction: the interpreter
        // calls it through a `Call`.
        return if graph.children(node).is_empty() {
            LEAF
        } else {
            OPAQUE_CALL
        };
    }
    // No scalar semantics at all: the term cannot be evaluated numerically.
    HEAVY
}

/// Weight of a named operator that the tape lowers by name.
fn named_cost(name: &str) -> Option<u64> {
    Some(match name {
        | "exp" | "ln" | "sin" | "cos" | "tan" | "asin" | "acos" | "atan" | "sinh" | "cosh"
        | "tanh" | "asinh" | "acosh" | "atanh" => TRANSCENDENTAL,
        | "sqrt" | "cbrt" => ROOT,
        | "atan2" => ATAN2,
        | "abs" | "floor" | "ceil" | "lt" | "le" | "gt" | "ge" | "ne" => ARITHMETIC,
        | _ => return None,
    })
}

/// Weight of an elementary intrinsic.
const fn intrinsic_cost(kind: Intrinsic) -> u64 {
    match kind {
        | Intrinsic::Sqrt | Intrinsic::Cbrt => ROOT,
        | Intrinsic::Abs | Intrinsic::Floor | Intrinsic::Ceil => ARITHMETIC,
        | Intrinsic::Exp
        | Intrinsic::Ln
        | Intrinsic::Sin
        | Intrinsic::Cos
        | Intrinsic::Tan
        | Intrinsic::Asin
        | Intrinsic::Acos
        | Intrinsic::Atan
        | Intrinsic::Sinh
        | Intrinsic::Cosh
        | Intrinsic::Tanh
        | Intrinsic::Asinh
        | Intrinsic::Acosh
        | Intrinsic::Atanh => TRANSCENDENTAL,
    }
}

/// Weight of one tape instruction.
const fn instruction_cost(inst: &Inst) -> u64 {
    match *inst {
        | Inst::Const(_) | Inst::Input(_) | Inst::Param(_) => LEAF,
        | Inst::Add(..)
        | Inst::Sub(..)
        | Inst::Mul(..)
        | Inst::Neg(..)
        | Inst::Cmp(..)
        | Inst::Select(..) => ARITHMETIC,
        | Inst::Div(..) | Inst::PowI(..) => DIVISION,
        | Inst::PowF(..) => POW_REAL,
        | Inst::Unary(kind, _) => intrinsic_cost(kind),
        | Inst::Atan2(..) => ATAN2,
        | Inst::Call(..) => OPAQUE_CALL,
        | Inst::Extern1(..) | Inst::Extern2(..) => EXTERN,
        | Inst::Gather(..) => GATHER,
    }
}

/// Sorted distinct symbols of the concrete term `node`.
fn symbols_of(
    graph: &Graph,
    node: NodeId,
) -> Vec<SymbolId> {
    let mut seen: HashSet<NodeId> = HashSet::new();
    let mut symbols = Vec::new();
    let mut stack = vec![node];
    while let Some(cur) = stack.pop() {
        if !seen.insert(cur) {
            continue;
        }
        if let Some(symbol) = graph.as_symbol(cur) {
            symbols.push(symbol);
        } else {
            stack.extend(graph.children(cur));
        }
    }
    symbols.sort_unstable();
    symbols.dedup();
    symbols
}

/// The declared cost of `node`'s operator, plus [`LIST_ARITHMETIC`] for a
/// sum or product with a list operand.
fn base_cost(
    graph: &Graph,
    node: NodeId,
) -> u64 {
    let op = graph.op(node);
    let mut cost = u64::from(graph.ops().get(op).cost);
    // A fraction is a little more complex than an integer: between two
    // spellings of one size, `-2^(1/2)` beats `-1/2*2^(3/2)`.
    if op == super::op::core::LIT && graph.as_number(node).is_some_and(|n| !n.is_integer()) {
        cost = cost.saturating_add(1);
    }
    if (op == super::op::core::ADD || op == super::op::core::MUL)
        && graph.children(node).iter().any(|&c| {
            graph
                .enodes(graph.find(c))
                .any(|n| graph.op(n) == super::op::core::LIST)
        })
    {
        cost = cost.saturating_add(LIST_ARITHMETIC);
    }
    cost
}

/// Like [`SizeCost`] but refuses heavy operators: only fully reduced terms
/// can be extracted.
#[derive(Copy, Clone, Debug, Default)]
pub struct ClosedForm;

impl CostModel for ClosedForm {
    fn cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        let desc = graph.ops().get(graph.op(node));
        (!desc.flags.has(OpFlags::HEAVY) || opaque_on_apply(graph, node, 8))
            .then(|| base_cost(graph, node))
    }
}

/// Whether `node` is an [`OpFlags::OPAQUE_ON_APPLY`] request applied to an
/// undetermined function (directly or through a chain of such requests).
fn opaque_on_apply(
    graph: &Graph,
    node: NodeId,
    depth: u8,
) -> bool {
    if depth == 0
        || !graph
            .ops()
            .get(graph.op(node))
            .flags
            .has(OpFlags::OPAQUE_ON_APPLY)
    {
        return false;
    }
    let Some(&first) = graph.children(node).first() else {
        return false;
    };
    graph.enodes(graph.find(first)).any(|n| {
        graph.op(n) == super::op::core::APPLY || opaque_on_apply(graph, n, depth.saturating_sub(1))
    })
}

/// Classes reachable from `roots` through e-node children, in discovery
/// order (roots first).
///
/// Class ids are dense arena indices, so membership is a bit-vector lookup
/// rather than a hash: reachability runs on every scheduler iteration and
/// inside every extraction.
#[must_use]
pub fn reachable(
    graph: &Graph,
    roots: &[NodeId],
) -> Vec<ClassId> {
    let mut seen = vec![false; graph.len()];
    let mut order = Vec::new();
    let mut stack: Vec<ClassId> = roots.iter().rev().map(|&r| graph.find(r)).collect();
    while let Some(class) = stack.pop() {
        let Some(flag) = seen.get_mut(class.index()) else {
            continue;
        };
        if *flag {
            continue;
        }
        *flag = true;
        order.push(class);
        // A class holding a literal is fully known: its other members
        // (`sin(pi)`, `x - x`, a reduced request ...) lead nowhere useful,
        // and following them would turn every constant into a hub that
        // makes the whole graph reachable from any term mentioning it.
        if graph
            .enodes(class)
            .any(|n| graph.op(n) == super::op::core::LIT)
        {
            continue;
        }
        for enode in graph.enodes(class) {
            for &child in graph.children(enode) {
                let child_class = graph.find(child);
                if !seen.get(child_class.index()).copied().unwrap_or(true) {
                    stack.push(child_class);
                }
            }
        }
    }
    order
}

/// The cheapest member of every class reachable from a set of roots.
#[derive(Clone, Debug)]
pub struct Extractor {
    /// Best `(cost, e-node)` per class, indexed by class id. Class ids are
    /// dense arena indices and this table is rebuilt on every scheduler
    /// iteration, so a direct vector beats hashing.
    best: Vec<Option<(u64, NodeId)>>,
}

impl Extractor {
    /// Computes best members under `model` for everything reachable from
    /// `roots`. The graph must be rebuilt (no pending unions).
    #[must_use]
    pub fn new(
        graph: &Graph,
        roots: &[NodeId],
        model: &dyn CostModel,
    ) -> Self {
        let classes = reachable(graph, roots);
        let mut best: Vec<Option<(u64, NodeId)>> = vec![None; graph.len()];
        // Bellman-Ford style relaxation. Processing in reverse discovery
        // order visits most children before their parents, so this usually
        // converges in two or three sweeps.
        for &class in &classes {
            if let Some(pin) = graph.pinned(class)
                && let Some(slot) = best.get_mut(class.index())
            {
                *slot = Some((1, pin));
            }
        }
        let mut changed = true;
        while changed {
            changed = false;
            for &class in classes.iter().rev() {
                if graph.pinned(class).is_some() {
                    continue;
                }
                for enode in graph.enodes(class) {
                    let Some(own) = model.cost(graph, enode) else {
                        continue;
                    };
                    let below = graph.children(enode).iter().try_fold(0_u64, |acc, &child| {
                        best.get(graph.find(child).index())
                            .copied()
                            .flatten()
                            .map(|(c, _)| acc.saturating_add(c))
                    });
                    let Some(below) = below else {
                        continue;
                    };
                    let total = own.saturating_add(below.saturating_mul(model.scale(graph, enode)));
                    // On a tie the younger node wins: simplification results
                    // are created after the terms they simplify.
                    let old = best.get(class.index()).copied().flatten();
                    let better = old.is_none_or(|(old_cost, old_node)| {
                        total < old_cost || (total == old_cost && enode > old_node)
                    });
                    if better {
                        if let Some(slot) = best.get_mut(class.index()) {
                            *slot = Some((total, enode));
                        }
                        changed = true;
                    }
                }
            }
        }
        Self { best }
    }

    /// Total cost of the best term of `node`'s class, or `None` when the
    /// model forbids every member.
    #[must_use]
    pub fn cost(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u64> {
        self.best
            .get(graph.find(node).index())
            .copied()
            .flatten()
            .map(|(c, _)| c)
    }

    /// Builds the best term of `node`'s class as a concrete node.
    ///
    /// Returns `None` when the class has no admissible term, or when the
    /// graph was changed by unions since the extractor was created in a way
    /// that makes the recorded choices cyclic.
    pub fn build(
        &self,
        graph: &mut Graph,
        node: NodeId,
    ) -> Option<NodeId> {
        let root = graph.find(node);
        let mut built: HashMap<ClassId, NodeId> = HashMap::new();
        let mut open: HashSet<ClassId> = HashSet::new();
        let mut stack = vec![(root, false)];
        while let Some((class, expanded)) = stack.pop() {
            if built.contains_key(&class) {
                continue;
            }
            if !expanded && !open.insert(class) {
                // The class is its own descendant: stale choices.
                return None;
            }
            let (_, enode) = self.best.get(class.index()).copied().flatten()?;
            if graph.pinned(class) == Some(enode) {
                built.insert(class, enode);
                continue;
            }
            let child_classes: Vec<ClassId> = graph
                .children(enode)
                .iter()
                .map(|&c| graph.find(c))
                .collect();
            if child_classes.is_empty() {
                built.insert(class, enode);
            } else if expanded {
                let mut children = Vec::with_capacity(child_classes.len());
                for child in &child_classes {
                    children.push(*built.get(child)?);
                }
                let term = graph.try_node(graph.op(enode), &children)?;
                built.insert(class, term);
            } else {
                stack.push((class, true));
                for child in child_classes {
                    if !built.contains_key(&child) {
                        stack.push((child, false));
                    }
                }
            }
        }
        built.get(&root).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;
    use crate::graph::op::core;

    #[test]
    fn picks_the_cheapest_member() {
        let mut g = Graph::new();
        let big = g.parse("a*b + a*c").unwrap_or(NodeId::NONE);
        let small = g.parse("a*(b + c)").unwrap_or(NodeId::NONE);
        g.union(big, small);
        g.rebuild();
        let ex = Extractor::new(&g, &[big], &SizeCost);
        assert_eq!(ex.build(&mut g, big), Some(small));
        assert_eq!(ex.cost(&g, big), ex.cost(&g, small));
    }

    #[test]
    fn improvements_propagate_upwards() {
        let mut g = Graph::new();
        assert!(
            g.ops_mut()
                .register(OpDescriptor::new("f", Arity::Fixed(1)))
                .is_ok()
        );
        let inner = g.parse("x + y + z").unwrap_or(NodeId::NONE);
        let outer = g.parse("f(f(x + y + z))").unwrap_or(NodeId::NONE);
        let w = g.sym("w");
        g.union(inner, w);
        g.rebuild();
        let ex = Extractor::new(&g, &[outer], &SizeCost);
        let out = ex.build(&mut g, outer);
        assert_eq!(out.map(|n| g.display(n)).as_deref(), Some("f(f(w))"));
    }

    #[test]
    fn closed_form_refuses_heavy_ops_until_reduced() {
        let mut g = Graph::new();
        let heavy = OpDescriptor::new("diff", Arity::Fixed(2)).flags(OpFlags::HEAVY);
        assert!(g.ops_mut().register(heavy).is_ok());
        let request = g.parse("diff(x^2, x) + 1").unwrap_or(NodeId::NONE);
        let ex = Extractor::new(&g, &[request], &ClosedForm);
        assert_eq!(ex.cost(&g, request), None);
        assert_eq!(ex.build(&mut g, request), None);
        assert!(
            Extractor::new(&g, &[request], &SizeCost)
                .cost(&g, request)
                .is_some()
        );

        let d = g.parse("diff(x^2, x)").unwrap_or(NodeId::NONE);
        let answer = g.parse("2*x").unwrap_or(NodeId::NONE);
        g.union(d, answer);
        g.rebuild();
        let ex = Extractor::new(&g, &[request], &ClosedForm);
        let out = ex.build(&mut g, request);
        assert_eq!(out.map(|n| g.display(n)).as_deref(), Some("2*x + 1"));
    }

    #[test]
    fn pinned_terms_are_returned_verbatim() {
        let mut g = Graph::new();
        assert!(
            g.ops_mut()
                .register(OpDescriptor::new("f", Arity::Fixed(1)))
                .is_ok()
        );
        let factored = g.parse("(a + b)^2").unwrap_or(NodeId::NONE);
        let expanded = g.parse("a^2 + 2*a*b + b^2").unwrap_or(NodeId::NONE);
        g.union(factored, expanded);
        g.pin(expanded);
        g.rebuild();
        let outer = g.parse("f((a + b)^2)").unwrap_or(NodeId::NONE);
        let ex = Extractor::new(&g, &[outer], &SizeCost);
        let out = ex.build(&mut g, outer);
        assert_eq!(
            out.map(|n| g.display(n)).as_deref(),
            Some("f(a^2 + 2*a*b + b^2)")
        );
        // A later pin does not displace the first.
        g.pin(factored);
        assert_eq!(g.pinned(g.find(factored)), Some(expanded));
    }

    #[test]
    fn cyclic_classes_terminate() {
        let mut g = Graph::new();
        assert!(
            g.ops_mut()
                .register(OpDescriptor::new("f", Arity::Fixed(1)))
                .is_ok()
        );
        let x = g.sym("x");
        let fx = g.parse("f(x)").unwrap_or(NodeId::NONE);
        // x = f(x): the class now contains a node whose child is the class.
        g.union(x, fx);
        g.rebuild();
        let ex = Extractor::new(&g, &[fx], &SizeCost);
        assert_eq!(ex.build(&mut g, fx), Some(x));
    }

    #[test]
    fn stale_extractor_fails_instead_of_looping() {
        let mut g = Graph::new();
        let term = g.parse("x * y").unwrap_or(NodeId::NONE);
        let ex = Extractor::new(&g, &[term], &SizeCost);
        // Merge the product with its own factor *after* extraction.
        let x = g.sym("x");
        g.union(term, x);
        g.rebuild();
        let out = ex.build(&mut g, term);
        assert!(out.is_none() || out == Some(x) || out == Some(term));
    }

    #[test]
    fn deep_terms_do_not_overflow_the_stack() {
        let mut g = Graph::new();
        let f = g
            .ops_mut()
            .register(OpDescriptor::new("f", Arity::Fixed(1)))
            .unwrap_or(crate::graph::OpId::NONE);
        let mut node = g.sym("x");
        for _ in 0..50_000 {
            node = g.node(f, &[node]);
        }
        let ex = Extractor::new(&g, &[node], &SizeCost);
        assert_eq!(ex.build(&mut g, node), Some(node));
        assert_eq!(ex.cost(&g, node), Some(150_003));
    }

    #[test]
    fn reachable_lists_roots_first_and_each_class_once() {
        let mut g = Graph::new();
        let root = g.parse("(a + b) * (a + b) * c").unwrap_or(NodeId::NONE);
        let classes = reachable(&g, &[root]);
        assert_eq!(classes.first(), Some(&g.find(root)));
        // root, a+b, a, b, c
        assert_eq!(classes.len(), 5);
    }

    /// Deterministic layered DAG with a few thousand distinct classes, used
    /// by the extraction benchmark below. A registered non-associative binary
    /// operator keeps the nodes from being flattened into one another, which
    /// is what a real mixed workload looks like.
    fn benchmark_graph() -> (Graph, NodeId) {
        let mut g = Graph::new();
        let op = g
            .ops_mut()
            .register(OpDescriptor::new("g", Arity::Fixed(2)))
            .unwrap_or(crate::graph::OpId::NONE);
        let mut nodes: Vec<NodeId> = (0..32).map(|i| g.sym(&format!("v{i}"))).collect();
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let pick = |state: &mut u64| {
            *state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (*state >> 33) as usize
        };
        let mut last = nodes[0];
        for _ in 0..2_000 {
            // A chain of nested operators with a shared random subterm keeps
            // the graph deep *and* wide, like a simplification run.
            let b = nodes[pick(&mut state) % nodes.len()];
            last = g.node(op, &[last, b]);
            nodes.push(last);
        }
        (g, last)
    }

    /// Reports deterministic timings of the reachability and extraction hot
    /// paths. Timing is printed, never asserted on.
    ///
    /// Measured on a development machine in the `test` profile (opt-level 1),
    /// 100 runs over the 2032-class DAG below: before the class-indexed
    /// tables `Extractor::new` took 19.2 ms/run and `reachable` 449 µs/run;
    /// after them, 8.5 ms/run (2.25×) and 149 µs/run (3.0×).
    #[test]
    fn extraction_hot_path_timing_report() {
        use std::hint::black_box;
        use std::time::Instant;

        let (g, root) = benchmark_graph();
        let classes = reachable(&g, &[root]);
        let (nodes, reachable_classes) = (g.len(), classes.len());
        println!("benchmark graph: {nodes} nodes, {reachable_classes} reachable classes");
        assert!(classes.len() > 1_500, "{} classes", classes.len());

        let runs = 100_u32;
        let start = Instant::now();
        let mut cost = None;
        for _ in 0..runs {
            let extractor = Extractor::new(&g, &[root], &SizeCost);
            cost = extractor.cost(&g, root);
            black_box(&extractor);
        }
        let per_extract = start.elapsed() / runs;
        let start = Instant::now();
        let mut count = 0;
        for _ in 0..runs {
            count += black_box(reachable(&g, &[root])).len();
        }
        let per_reachable = start.elapsed() / runs;
        assert!(count > 0 && cost.is_some());
        let reachable_count = classes.len();
        println!(
            "hot path: {reachable_count} classes; Extractor::new {per_extract:?}/run; \
             reachable {per_reachable:?}/run"
        );
    }

    /// Registers a real-valued unary `sin` without the rules, for weight
    /// tests that must not depend on a rule set.
    fn register_sin(g: &mut Graph) -> crate::graph::OpId {
        g.ops_mut()
            .register(
                OpDescriptor::new("sin", Arity::Fixed(1))
                    .eval(|a| a.first().map_or(f64::NAN, |x| x.sin())),
            )
            .unwrap_or(crate::graph::OpId::NONE)
    }

    #[test]
    fn eval_cost_term_weights_follow_the_tape() {
        let mut g = Graph::new();
        let sin = register_sin(&mut g);
        let x = g.sym("x");
        let one = g.int(1);
        let two = g.int(2);
        let sin_x = g.node(sin, &[x]);
        let square = g.node(core::POW, &[x, two]);
        let sum = g.node(core::ADD, &[x, one]);
        let model = EvalCost::new();
        assert_eq!(model.term_cost(&g, x), Some(1), "parameter load");
        assert_eq!(
            model.term_cost(&g, sin_x),
            Some(11),
            "load plus a transcendental"
        );
        assert_eq!(
            model.term_cost(&g, square),
            Some(5),
            "load plus an integer power"
        );
        assert_eq!(model.term_cost(&g, sum), Some(3), "two loads plus an add");
        assert_eq!(model.cached_terms(), 4);
        // The extractor's per-node prices are the same instruction classes.
        assert_eq!(model.cost(&g, sum), Some(1));
        assert_eq!(model.cost(&g, square), Some(POW_INT));
        assert_eq!(model.cost(&g, sin_x), Some(TRANSCENDENTAL));
        assert_eq!(model.cost(&g, x), Some(LEAF));
        assert_eq!(model.cached_nodes(), 4);
    }

    #[test]
    fn eval_cost_prefers_a_cheaper_operation_where_size_ties() {
        let mut g = Graph::new();
        let x = g.sym("x");
        let two = g.int(2);
        // Both spellings have the same term size (operator + symbol + literal),
        // and only the evaluation cost separates a multiplication from a
        // power. The product is created first so that size alone — which
        // breaks ties towards the younger node — would pick the power.
        let product = g.node(core::MUL, &[x, two]);
        let square = g.node(core::POW, &[x, two]);
        g.union(product, square);
        g.rebuild();
        let by_size = Extractor::new(&g, &[product], &SizeCost);
        let by_eval = Extractor::new(&g, &[product], &EvalCost::new());
        assert_eq!(by_size.cost(&g, product), Some(5));
        assert_eq!(by_eval.cost(&g, product), Some(3));
        let size_term = by_size.build(&mut g, product).map(|n| g.display(n));
        let eval_term = by_eval.build(&mut g, product).map(|n| g.display(n));
        assert_eq!(
            size_term.as_deref(),
            Some("x^2"),
            "a size tie goes to the younger node"
        );
        assert_eq!(
            eval_term.as_deref(),
            Some("2*x"),
            "a multiplication is cheaper to evaluate than a power"
        );
    }

    #[test]
    fn eval_cost_agrees_with_size_cost_on_a_trivial_family() {
        let mut g = Graph::new();
        let (a, b, c) = (g.sym("a"), g.sym("b"), g.sym("c"));
        let small = g.node(core::ADD, &[a, b]);
        let large = g.node(core::ADD, &[a, b, c]);
        g.union(small, large);
        g.rebuild();
        let by_size = Extractor::new(&g, &[small], &SizeCost);
        let by_eval = Extractor::new(&g, &[small], &EvalCost::new());
        assert_eq!(
            by_size
                .build(&mut g, small)
                .map(|n| g.display(n))
                .as_deref(),
            Some("a + b")
        );
        assert_eq!(
            by_eval
                .build(&mut g, small)
                .map(|n| g.display(n))
                .as_deref(),
            Some("a + b")
        );
    }

    #[test]
    fn cached_eval_cost_extractions_are_identical() {
        let mut g = Graph::new();
        let x = g.sym("x");
        let y = g.sym("y");
        let product = g.node(core::MUL, &[x, y]);
        let two = g.int(2);
        let square = g.node(core::POW, &[x, two]);
        let model = EvalCost::new();
        let first = Extractor::new(&g, &[product], &model);
        let second = Extractor::new(&g, &[product], &model);
        assert!(model.cached_nodes() > 0, "the model caches node prices");
        let first_term = first.build(&mut g, product).map(|n| g.display(n));
        let second_term = second.build(&mut g, product).map(|n| g.display(n));
        assert_eq!(first_term, second_term);
        assert_eq!(first.cost(&g, product), second.cost(&g, product));
        // The class-level view uses the cached per-term lowerings.
        g.union(product, square);
        g.rebuild();
        assert!(model.class_cost(&g, product).is_some());
        let cached = model.class_cost(&g, product);
        assert_eq!(cached, model.class_cost(&g, product));
        assert!(model.cached_classes() >= 1 && model.cached_terms() >= 1);
    }

    #[test]
    fn eval_cost_benchmark_report() {
        use std::time::Instant;

        // A large term built from real arithmetic: a chain of products, sums
        // and powers with shared operands.
        let mut g = Graph::new();
        let mut nodes: Vec<NodeId> = (0..16).map(|i| g.sym(&format!("v{i}"))).collect();
        let mut state: u64 = 0x243f_6a88_85a3_08d3;
        let pick = |state: &mut u64| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            (*state >> 11) as usize
        };
        let two = g.int(2);
        let three = g.int(3);
        let mut last = nodes[0];
        for i in 0..1_500 {
            let b = nodes[pick(&mut state) % nodes.len()];
            last = match i % 4 {
                | 0 => g.node(core::ADD, &[last, b]),
                | 1 => g.node(core::MUL, &[last, b]),
                | 2 => g.node(core::POW, &[last, two]),
                | _ => {
                    let scaled = g.node(core::MUL, &[b, three]);
                    g.node(core::ADD, &[last, scaled])
                },
            };
            nodes.push(last);
        }
        let model = EvalCost::new();
        let start = Instant::now();
        let cost = model.class_cost(&g, last);
        let first = start.elapsed();
        let start = Instant::now();
        let again = model.class_cost(&g, last);
        let cached = start.elapsed();
        assert_eq!(cost, again);
        let nodes = g.len();
        println!("EvalCost: {nodes} nodes, term cost {cost:?}; first {first:?}, cached {cached:?}");
    }
}
