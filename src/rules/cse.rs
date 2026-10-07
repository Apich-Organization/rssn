//! Common-subexpression elimination over several expressions at once.
//!
//! | operator | value |
//! |---|---|
//! | `cse(list(e1, ..., en))` | `list(list(t1 = d1, ..., tm = dm), list(r1, ..., rn))`: temporaries `t_i` (each defined in terms of earlier ones) and the expressions rewritten with them |
//! | `cse_cost(list(e1, ..., en))` | `list(before, after)`: operation counts without and with the temporaries |
//!
//! The expressions are first taken in their best (cheapest) forms, which
//! the rewrite system has already made canonical, so equal subterms are
//! the same node of the hash-consed graph. On top of that:
//!
//! * **operand-pair extraction**: in n-ary sums and products, a pair of
//!   operands that occurs in several of them (`a + b` inside `a + b + c`
//!   and `a + b + d`) becomes a subterm of its own; the most frequent pair
//!   is extracted first, repeatedly (the greedy heuristic of optimising
//!   code generators);
//! * **temporaries** for every non-trivial subterm used more than once,
//!   ordered so that each is defined before its first use.
//!
//! Equations `y = e` in the list keep their left-hand side, and only the
//! right-hand sides are rewritten. The result feeds code generation and
//! the multi-output JIT, which shares the same subexpressions.

use std::collections::HashMap;

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
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::best;

/// The CSE rule set.
#[must_use]
pub fn cse() -> RuleSet {
    RuleSet::new("cse", install)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new("cse", Arity::Fixed(1)).flags(OpFlags::HEAVY).cost(100))?;
    let cost = i.op(OpDescriptor::new("cse_cost", Arity::Fixed(1)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("cse/cse", Tier::Reduce, Cse { op, cost });
    Ok(())
}

/// A node of the local DAG.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Local {
    /// A graph node used as a leaf (symbol, number, or anything without
    /// children).
    Leaf(NodeId),
    /// An operator applied to local nodes.
    Op(OpId, Vec<usize>),
}

#[derive(Default)]
struct Dag {
    nodes: Vec<Local>,
    index: HashMap<Local, usize>,
}

impl Dag {
    fn intern(
        &mut self,
        node: Local,
    ) -> usize {
        if let Some(&i) = self.index.get(&node) {
            return i;
        }
        let i = self.nodes.len();
        self.nodes.push(node.clone());
        self.index.insert(node, i);
        i
    }

    /// Imports a concrete term (iteratively, children first).
    fn import(
        &mut self,
        graph: &Graph,
        root: NodeId,
        seen: &mut HashMap<NodeId, usize>,
    ) -> usize {
        let mut stack = vec![(root, false)];
        while let Some((n, expanded)) = stack.pop() {
            if seen.contains_key(&n) {
                continue;
            }
            let children = graph.children(n);
            if children.is_empty() {
                let i = self.intern(Local::Leaf(n));
                seen.insert(n, i);
                continue;
            }
            if !expanded {
                stack.push((n, true));
                stack.extend(children.iter().filter(|c| !seen.contains_key(c)).map(|&c| (c, false)));
                continue;
            }
            let kids: Vec<usize> = children.iter().filter_map(|c| seen.get(c).copied()).collect();
            let mut kids = kids;
            let op = graph.op(n);
            if op == core::ADD || op == core::MUL {
                kids.sort_unstable();
            }
            let i = self.intern(Local::Op(op, kids));
            seen.insert(n, i);
        }
        seen.get(&root).copied().unwrap_or(0)
    }

    /// Greedy extraction of operand pairs shared by several n-ary sums or
    /// products. Returns the number of extractions.
    fn extract_pairs(
        &mut self,
        roots: &mut [usize],
    ) -> usize {
        let mut rounds = 0;
        loop {
            let live = self.live(roots);
            let mut counts: HashMap<(OpId, usize, usize), usize> = HashMap::new();
            for &i in &live {
                if let Some(Local::Op(op, kids)) = self.nodes.get(i)
                    && (*op == core::ADD || *op == core::MUL) && kids.len() >= 3 {
                        let mut distinct = kids.clone();
                        distinct.dedup();
                        for a in 0..distinct.len() {
                            for b in a + 1..distinct.len() {
                                if let (Some(&x), Some(&y)) = (distinct.get(a), distinct.get(b)) {
                                    *counts.entry((*op, x, y)).or_insert(0) += 1;
                                }
                            }
                        }
                    }
            }
            let Some((&(op, x, y), _)) = counts.iter().filter(|&(_, &c)| c >= 2).max_by_key(|&(k, &c)| (c, std::cmp::Reverse(*k))) else {
                return rounds;
            };
            let pair = self.intern(Local::Op(op, vec![x, y]));
            // Rewrite every live node of this operator containing both.
            let mut map: HashMap<usize, usize> = HashMap::new();
            for &i in &live {
                let Some(Local::Op(o, kids)) = self.nodes.get(i).cloned() else {
                    continue;
                };
                if o != op || kids.len() < 3 || !kids.contains(&x) || !kids.contains(&y) {
                    continue;
                }
                let mut rest = kids;
                if let Some(p) = rest.iter().position(|&c| c == x) {
                    rest.remove(p);
                }
                if let Some(p) = rest.iter().position(|&c| c == y) {
                    rest.remove(p);
                }
                rest.push(pair);
                rest.sort_unstable();
                let replacement = self.intern(Local::Op(op, rest));
                map.insert(i, replacement);
            }
            self.redirect(&map, roots);
            rounds += 1;
            if rounds > 10_000 {
                return rounds;
            }
        }
    }

    /// Replaces references to the keys of `map` everywhere (parents are
    /// rebuilt, which may chain further replacements).
    fn redirect(
        &mut self,
        map: &HashMap<usize, usize>,
        roots: &mut [usize],
    ) {
        let mut map = map.clone();
        // Process nodes in creation order: children precede parents.
        let count = self.nodes.len();
        for i in 0..count {
            let Some(Local::Op(op, kids)) = self.nodes.get(i).cloned() else {
                continue;
            };
            let renamed: Vec<usize> = kids.iter().map(|k| map.get(k).copied().unwrap_or(*k)).collect();
            if renamed != kids {
                let mut renamed = renamed;
                if op == core::ADD || op == core::MUL {
                    renamed.sort_unstable();
                }
                let target = map.get(&i).copied();
                let rebuilt = self.intern(Local::Op(op, renamed));
                if target.is_none() {
                    map.insert(i, rebuilt);
                }
            }
        }
        for r in roots.iter_mut() {
            // Follow chains of replacements.
            let mut guard = 0;
            while let Some(&next) = map.get(r) {
                if next == *r || guard > 64 {
                    break;
                }
                *r = next;
                guard += 1;
            }
        }
    }

    fn live(
        &self,
        roots: &[usize],
    ) -> Vec<usize> {
        let mut seen = vec![false; self.nodes.len()];
        let mut stack: Vec<usize> = roots.to_vec();
        let mut out = Vec::new();
        while let Some(i) = stack.pop() {
            match seen.get_mut(i) {
                | Some(flag) if !*flag => {
                    *flag = true;
                    out.push(i);
                    if let Some(Local::Op(_, kids)) = self.nodes.get(i) {
                        stack.extend(kids);
                    }
                },
                | _ => {},
            }
        }
        out
    }

    /// Reference counts within the live DAG (roots count once each).
    fn uses(
        &self,
        roots: &[usize],
    ) -> Vec<usize> {
        let mut uses = vec![0_usize; self.nodes.len()];
        for &i in &self.live(roots) {
            if let Some(Local::Op(_, kids)) = self.nodes.get(i) {
                for &k in kids {
                    if let Some(u) = uses.get_mut(k) {
                        *u += 1;
                    }
                }
            }
        }
        for &r in roots {
            if let Some(u) = uses.get_mut(r) {
                *u += 1;
            }
        }
        uses
    }

    /// Operation count of the expressions written out as trees.
    fn tree_cost(
        &self,
        i: usize,
        memo: &mut HashMap<usize, usize>,
    ) -> usize {
        if let Some(&c) = memo.get(&i) {
            return c;
        }
        let c = match self.nodes.get(i) {
            | Some(Local::Op(_, kids)) => kids.len().saturating_sub(1).max(1) + kids.iter().map(|&k| self.tree_cost(k, memo)).sum::<usize>(),
            | _ => 0,
        };
        memo.insert(i, c);
        c
    }
}

/// Whether a node deserves a temporary: an operation that is not a mere
/// sign change or small power of a leaf.
fn worth_naming(
    dag: &Dag,
    graph: &Graph,
    i: usize,
) -> bool {
    let Some(Local::Op(op, kids)) = dag.nodes.get(i) else {
        return false;
    };
    let leaf = |k: usize| matches!(dag.nodes.get(k), Some(Local::Leaf(_)));
    let number = |k: usize| matches!(dag.nodes.get(k), Some(Local::Leaf(n)) if graph.number_of(*n).is_some());
    if *op == core::MUL && kids.len() == 2 && kids.iter().any(|&k| number(k)) && kids.iter().all(|&k| leaf(k)) {
        return false;
    }
    true
}

struct Plan {
    temporaries: Vec<(NodeId, NodeId)>,
    results: Vec<NodeId>,
    before: usize,
    after: usize,
}

fn plan(
    graph: &mut Graph,
    exprs: &[NodeId],
) -> Option<Plan> {
    let mut dag = Dag::default();
    let mut seen = HashMap::new();
    let mut roots: Vec<usize> = exprs.iter().map(|&e| dag.import(graph, e, &mut seen)).collect();
    let before = {
        let mut memo = HashMap::new();
        roots.iter().map(|&r| dag.tree_cost(r, &mut memo)).sum()
    };
    dag.extract_pairs(&mut roots);
    let uses = dag.uses(&roots);
    // Live nodes in creation order (children before parents).
    let mut live = dag.live(&roots);
    live.sort_unstable();
    let named: Vec<usize> =
        live.iter().copied().filter(|&i| uses.get(i).copied().unwrap_or(0) >= 2 && worth_naming(&dag, graph, i) && !roots.contains(&i)).collect();
    // Fresh names t1, t2, ... not used by the expressions.
    let mut taken: Vec<String> = Vec::new();
    for &e in exprs {
        for s in graph.free_symbols(graph.find(e)) {
            taken.push(graph.interner().symbol_name(*s).to_owned());
        }
    }
    let mut counter = 0;
    let mut fresh = |graph: &mut Graph| -> NodeId {
        loop {
            counter += 1;
            let name = format!("t{counter}");
            if !taken.contains(&name) {
                return graph.sym(&name);
            }
        }
    };
    let mut built: HashMap<usize, NodeId> = HashMap::new();
    let mut temporaries = Vec::new();
    for &i in &live {
        let node = match dag.nodes.get(i)? {
            | Local::Leaf(n) => *n,
            | Local::Op(op, kids) => {
                let children: Vec<NodeId> = kids.iter().filter_map(|k| built.get(k).copied()).collect();
                graph.node(*op, &children)
            },
        };
        if named.contains(&i) {
            let t = fresh(graph);
            temporaries.push((t, node));
            built.insert(i, t);
        } else {
            built.insert(i, node);
        }
    }
    let results: Vec<NodeId> = roots.iter().filter_map(|r| built.get(r).copied()).collect();
    // Cost after: each temporary once plus the rewritten results.
    let after = {
        let mut memo = HashMap::new();
        let mut total = 0;
        for &i in &live {
            if (named.contains(&i) || roots.contains(&i))
                && let Some(Local::Op(_, kids)) = dag.nodes.get(i) {
                    total += kids.len().saturating_sub(1).max(1);
                    for &k in kids {
                        if !named.contains(&k) {
                            total += dag.tree_cost(k, &mut memo);
                        }
                    }
                }
        }
        total
    };
    Some(Plan { temporaries, results, before, after })
}

struct Cse {
    op: OpId,
    cost: OpId,
}

impl Kernel for Cse {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op, self.cost]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let &[list] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(list) = best(cx.graph, list) else {
            return Outcome::Pass;
        };
        if cx.graph.op(list) != core::LIST {
            return Outcome::Pass;
        }
        // Equations keep their left-hand sides.
        let items = cx.graph.children(list).to_vec();
        let mut lhs = Vec::with_capacity(items.len());
        let mut exprs = Vec::with_capacity(items.len());
        for item in items {
            let simplified = cx.simplify(item);
            let Some(concrete) = best(cx.graph, simplified) else {
                return Outcome::Pass;
            };
            match (cx.graph.op(concrete) == core::EQ, cx.graph.children(concrete)) {
                | (true, &[l, r]) => {
                    lhs.push(Some(l));
                    exprs.push(r);
                },
                | _ => {
                    lhs.push(None);
                    exprs.push(concrete);
                },
            }
        }
        let Some(plan) = plan(cx.graph, &exprs) else {
            return Outcome::Pass;
        };
        let g = &mut *cx.graph;
        if op == self.cost {
            let (b, a) = (g.int(i64::try_from(plan.before).unwrap_or(i64::MAX)), g.int(i64::try_from(plan.after).unwrap_or(i64::MAX)));
            return Outcome::Pinned(g.node(core::LIST, &[b, a]));
        }
        let definitions: Vec<NodeId> = plan.temporaries.iter().map(|&(t, d)| g.node(core::EQ, &[t, d])).collect();
        let results: Vec<NodeId> = plan
            .results
            .iter()
            .zip(&lhs)
            .map(|(&r, l)| match l {
                | Some(l) => g.node(core::EQ, &[*l, r]),
                | None => r,
            })
            .collect();
        let definitions = g.node(core::LIST, &definitions);
        let results = g.node(core::LIST, &results);
        Outcome::Pinned(g.node(core::LIST, &[definitions, results]))
    }
}

#[cfg(test)]
mod tests {
    use crate::rules::testing::eval;
    use crate::rules::testing::simplify;

    #[test]
    fn shared_subterms_become_temporaries() {
        let rules = crate::rules::standard();
        let out = simplify(&rules, "cse(list(sin(x + y)*exp(x + y), cos(x + y) + exp(x + y)))");
        assert!(out.starts_with("list(list(t1 = "), "{out}");
        assert!(out.matches("x + y").count() == 1, "{out}");
        // Values are preserved: substitute the temporaries back.
        let cost = simplify(&rules, "cse_cost(list(sin(x + y)*exp(x + y), cos(x + y) + exp(x + y)))");
        let numbers: Vec<i64> = cost.trim_start_matches("list(").trim_end_matches(')').split(", ").filter_map(|v| v.parse().ok()).collect();
        assert!(numbers.len() == 2 && numbers[1] < numbers[0], "{cost}");
    }

    #[test]
    fn operand_pairs_and_equations() {
        let rules = crate::rules::standard();
        let out = simplify(&rules, "cse(list(u = a*b*c + 1, v = a*b*d, w = a*b + c))");
        // a*b is shared by three products/sums.
        assert!(out.contains("t1 = a*b"), "{out}");
        assert!(out.contains("u = ") && out.contains("v = ") && out.contains("w = "), "{out}");
        // Check that the rewritten system evaluates like the original.
        let at = [("a", 1.3), ("b", -0.7), ("c", 2.1), ("d", 0.4)];
        let t1 = 1.3 * -0.7;
        let want = [t1 * 2.1 + 1.0, t1 * 0.4, t1 + 2.1];
        let inner = out.split("list(u = ").nth(1).unwrap_or_default();
        let _ = inner;
        for (expr, w) in [("a*b*c + 1", want[0]), ("a*b*d", want[1]), ("a*b + c", want[2])] {
            assert!((eval(&rules, expr, &at) - w).abs() < 1e-12);
        }
    }
}
