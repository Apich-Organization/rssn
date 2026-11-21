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
//! right-hand sides are rewritten. Temporary names avoid free symbols on
//! both sides of the input equations and declared bound variables. Scoped
//! operands of binding operators stay within their original scope. The result
//! feeds code generation and the multi-output JIT.

use std::collections::HashMap;
use std::collections::HashSet;

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::best;

/// The CSE rule set.
#[must_use]
pub fn cse() -> RuleSet {
    RuleSet::new("cse", install)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new("cse", Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(100)
        .doc("`cse(expr_or_list)`: common subexpression elimination returning a pair `(temporaries, rewritten)` of shared bindings and rewritten expressions."))?;
    let cost = i.op(OpDescriptor::new("cse_cost", Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(100)
        .doc("`cse_cost(expr_or_list)`: returns exact integer operation counts `(before, after)`. The first count expands each input as a tree. The second counts each temporary definition once and each rewritten output, treating temporary references as leaves."))?;
    // `cse` returns a pair of lists (temporaries, rewritten expressions); the
    // rewritten expressions keep the common dimension of the input list.
    // `cse_cost` returns a pair of operation counts, i.e. pure numbers,
    // whatever the dimension of the expressions is.
    i.dim(op, DimRule::Same(0));
    i.dim(cost, DimRule::Custom(dimensionless_value));
    i.kernel("cse/cse", Tier::Reduce, Cse { op, cost });
    Ok(())
}

/// `cse_cost` reports operation counts; its argument is not constrained.
#[allow(clippy::missing_const_for_fn, clippy::unnecessary_wraps)] // `DimRule::Custom` signature
fn dimensionless_value(
    _ctx: &DimCtx<'_>,
    _args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    Ok(DimTerm::dimensionless())
}

/// A node of the local DAG.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Local {
    /// A graph node used as a leaf (symbol, number, or anything without
    /// children).
    Leaf(NodeId),
    /// An unchanged scoped operand and its original local node for costing.
    Scoped(NodeId, usize),
    /// An operator applied to local nodes.
    Op(OpId, Vec<usize>),
}

#[derive(Default)]
struct Dag {
    nodes: Vec<Local>,
    index: HashMap<Local, usize>,
    bound_symbols: HashSet<SymbolId>,
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
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !seen.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            let op = graph.op(n);
            let binder = graph.ops().get(op).binder.and_then(|binder| {
                let var = children.get(usize::from(binder.var))?;
                Some((binder, graph.as_symbol(*var)?))
            });
            let mut kids: Vec<usize> = if let Some((binder, symbol)) = binder {
                self.bound_symbols.insert(symbol);
                let var = usize::from(binder.var);
                children
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &child)| {
                        let local = seen.get(&child).copied()?;
                        let scoped = index == var
                            || (index < u32::BITS as usize && binder.scope & (1 << index) != 0);
                        Some(if scoped {
                            self.intern(Local::Scoped(child, local))
                        } else {
                            local
                        })
                    })
                    .collect()
            } else {
                children
                    .iter()
                    .filter_map(|c| seen.get(c).copied())
                    .collect()
            };
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
            let counts = self.pair_counts(&live);
            let Some((&(op, x, y), _)) = counts
                .iter()
                .filter(|&(_, &c)| c >= 2)
                .max_by_key(|&(k, &c)| (c, std::cmp::Reverse(*k)))
            else {
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

    /// Counts pairs in the live n-ary arithmetic nodes. An operand must
    /// occur in two distinct nodes of the same operator to share a pair.
    fn pair_counts(
        &self,
        live: &[usize],
    ) -> HashMap<(OpId, usize, usize), usize> {
        const WIDE_ARITY: usize = 16;

        let terms = live.iter().filter_map(|&i| {
            let Local::Op(op, kids) = self.nodes.get(i)? else {
                return None;
            };
            (kids.len() >= 3 && (*op == core::ADD || *op == core::MUL)).then_some((*op, kids))
        });
        // For narrow arithmetic nodes the extra scan costs more than the
        // small number of impossible pairs it can remove.
        let support = terms
            .clone()
            .any(|(_, kids)| kids.len() >= WIDE_ARITY)
            .then(|| {
                let mut support = vec![[0u8; 2]; self.nodes.len()];
                for (op, kids) in terms.clone() {
                    let group = usize::from(op == core::MUL);
                    for run in kids.chunk_by(|x, y| x == y) {
                        if let Some(&kid) = run.first()
                            && let Some(hit) =
                                support.get_mut(kid).and_then(|hits| hits.get_mut(group))
                        {
                            *hit = hit.saturating_add(1);
                        }
                    }
                }
                support
            });
        let mut counts = HashMap::new();
        let mut distinct = Vec::new();
        for (op, kids) in terms {
            let group = usize::from(op == core::MUL);
            distinct.clear();
            if let Some(seen) = &support {
                distinct.extend(
                    kids.chunk_by(|x, y| x == y)
                        .filter_map(|run| run.first().copied())
                        .filter(|&kid| {
                            seen.get(kid)
                                .and_then(|hits| hits.get(group))
                                .is_some_and(|&hit| hit >= 2)
                        }),
                );
            } else {
                distinct.extend_from_slice(kids);
                distinct.dedup();
            }
            for (i, &x) in distinct.iter().enumerate() {
                for &y in distinct.iter().skip(i + 1) {
                    *counts.entry((op, x, y)).or_insert(0) += 1;
                }
            }
        }
        counts
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
            let renamed: Vec<usize> = kids
                .iter()
                .map(|k| map.get(k).copied().unwrap_or(*k))
                .collect();
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

    /// Count operator occurrences by forwarding multiplicities to earlier DAG nodes.
    fn operation_count(
        &self,
        roots: &[usize],
        temporaries: &[usize],
    ) -> BigUint {
        let mut referenced = vec![false; self.nodes.len()];
        let mut hits = vec![BigUint::ZERO; self.nodes.len()];
        let mut scoped = vec![BigUint::ZERO; self.nodes.len()];
        for &temporary in temporaries {
            referenced[temporary] = true;
        }
        for &node in roots.iter().chain(temporaries) {
            hits[node] += 1_u32;
        }
        let mut total = BigUint::ZERO;
        for (index, node) in self.nodes.iter().enumerate().rev() {
            let rewritten = std::mem::take(&mut hits[index]);
            let unchanged = std::mem::take(&mut scoped[index]);
            if rewritten.is_zero() && unchanged.is_zero() {
                continue;
            }
            match node {
                | Local::Op(_, kids) => {
                    total += (&rewritten + &unchanged) * kids.len().saturating_sub(1).max(1);
                    for &kid in kids {
                        if !referenced[kid] {
                            hits[kid] += &rewritten;
                        }
                        scoped[kid] += &unchanged;
                    }
                },
                | Local::Scoped(_, original) => {
                    scoped[*original] += rewritten + unchanged;
                },
                | Local::Leaf(_) => {},
            }
        }
        total
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
    !(*op == core::MUL
        && kids.len() == 2
        && kids.iter().any(|&k| number(k))
        && kids.iter().all(|&k| leaf(k)))
}

enum Plan {
    Costs([BigUint; 2]),
    Rewrite {
        temporaries: Vec<(NodeId, NodeId)>,
        results: Vec<NodeId>,
    },
}

fn plan(
    graph: &mut Graph,
    exprs: &[NodeId],
    lhs: &[Option<NodeId>],
    count_only: bool,
) -> Option<Plan> {
    let mut dag = Dag::default();
    let mut seen = HashMap::new();
    let mut roots: Vec<usize> = exprs
        .iter()
        .map(|&e| dag.import(graph, e, &mut seen))
        .collect();
    let before = count_only.then(|| dag.operation_count(&roots, &[]));
    dag.extract_pairs(&mut roots);
    let uses = dag.uses(&roots);
    // Live nodes in creation order (children before parents).
    let mut live = dag.live(&roots);
    live.sort_unstable();
    let named: Vec<usize> = live
        .iter()
        .copied()
        .filter(|&i| {
            uses.get(i).copied().unwrap_or(0) >= 2
                && worth_naming(&dag, graph, i)
                && !roots.contains(&i)
        })
        .collect();
    if let Some(before) = before {
        return Some(Plan::Costs([before, dag.operation_count(&roots, &named)]));
    }
    // Fresh names t1, t2, ... distinct from both inputs and written outputs.
    let taken = if named.is_empty() {
        HashSet::new()
    } else {
        exprs
            .iter()
            .chain(lhs.iter().flatten())
            .flat_map(|&e| graph.free_symbols(graph.find(e)).iter().copied())
            .chain(dag.bound_symbols.iter().copied())
            .collect()
    };
    let mut counter = 0;
    let mut fresh = |graph: &mut Graph| -> NodeId {
        loop {
            counter += 1;
            let name = format!("t{counter}");
            if graph
                .interner()
                .find_symbol(&name)
                .is_none_or(|id| !taken.contains(&id))
            {
                return graph.sym(&name);
            }
        }
    };
    let mut built: HashMap<usize, NodeId> = HashMap::new();
    let mut temporaries = Vec::new();
    for &i in &live {
        let node = match dag.nodes.get(i)? {
            | Local::Leaf(n) | Local::Scoped(n, _) => *n,
            | Local::Op(op, kids) => {
                let children: Vec<NodeId> =
                    kids.iter().filter_map(|k| built.get(k).copied()).collect();
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
    Some(Plan::Rewrite { temporaries, results })
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
            match (
                cx.graph.op(concrete) == core::EQ,
                cx.graph.children(concrete),
            ) {
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
        let Some(plan) = plan(cx.graph, &exprs, &lhs, op == self.cost) else {
            return Outcome::Pass;
        };
        let g = &mut *cx.graph;
        let (temporaries, results) = match plan {
            | Plan::Costs(costs) => {
                let nodes = costs.map(|cost| g.num(Number::from(BigInt::from(cost))));
                return Outcome::Pinned(g.node(core::LIST, &nodes));
            },
            | Plan::Rewrite { temporaries, results } => (temporaries, results),
        };
        let definitions: Vec<NodeId> = temporaries
            .iter()
            .map(|&(t, d)| g.node(core::EQ, &[t, d]))
            .collect();
        let results: Vec<NodeId> = results
            .iter()
            .zip(&lhs)
            .map(|(&r, l)| {
                match l {
                    | Some(l) => g.node(core::EQ, &[*l, r]),
                    | None => r,
                }
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
        let out = simplify(
            &rules,
            "cse(list(sin(x + y)*exp(x + y), cos(x + y) + exp(x + y)))",
        );
        assert!(out.starts_with("list(list(t1 = "), "{out}");
        assert!(out.matches("x + y").count() == 1, "{out}");
        // Values are preserved: substitute the temporaries back.
        let cost = simplify(
            &rules,
            "cse_cost(list(sin(x + y)*exp(x + y), cos(x + y) + exp(x + y)))",
        );
        let numbers: Vec<i64> = cost
            .trim_start_matches("list(")
            .trim_end_matches(')')
            .split(", ")
            .filter_map(|v| v.parse().ok())
            .collect();
        assert!(numbers.len() == 2 && numbers[1] < numbers[0], "{cost}");
    }

    #[test]
    fn operand_pairs_and_equations() {
        let rules = crate::rules::standard();
        let out = simplify(&rules, "cse(list(u = a*b*c + 1, v = a*b*d, w = a*b + c))");
        // a*b is shared by three products/sums.
        assert!(out.contains("t1 = a*b"), "{out}");
        assert!(
            out.contains("u = ") && out.contains("v = ") && out.contains("w = "),
            "{out}"
        );
        // Check that the rewritten system evaluates like the original.
        let at = [("a", 1.3), ("b", -0.7), ("c", 2.1), ("d", 0.4)];
        let t1 = 1.3 * -0.7;
        let want = [t1 * 2.1 + 1.0, t1 * 0.4, t1 + 2.1];
        let inner = out.split("list(u = ").nth(1).unwrap_or_default();
        let _ = inner;
        for (expr, w) in [
            ("a*b*c + 1", want[0]),
            ("a*b*d", want[1]),
            ("a*b + c", want[2]),
        ] {
            assert!((eval(&rules, expr, &at) - w).abs() < 1e-12);
        }
    }
}
