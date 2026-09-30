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

use std::collections::HashMap;
use std::collections::HashSet;

use super::id::ClassId;
use super::id::NodeId;
use super::op::OpFlags;
use super::store::Graph;

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
        Some(u64::from(graph.ops().get(graph.op(node)).cost))
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
        (!desc.flags.has(OpFlags::HEAVY)).then_some(u64::from(desc.cost))
    }
}

/// Classes reachable from `roots` through e-node children, in discovery
/// order (roots first).
#[must_use]
pub fn reachable(
    graph: &Graph,
    roots: &[NodeId],
) -> Vec<ClassId> {
    let mut seen: HashSet<ClassId> = HashSet::new();
    let mut order = Vec::new();
    let mut stack: Vec<ClassId> = roots.iter().rev().map(|&r| graph.find(r)).collect();
    while let Some(class) = stack.pop() {
        if !seen.insert(class) {
            continue;
        }
        order.push(class);
        for enode in graph.enodes(class) {
            for &child in graph.children(enode) {
                let child_class = graph.find(child);
                if !seen.contains(&child_class) {
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
    best: HashMap<ClassId, (u64, NodeId)>,
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
        let mut best: HashMap<ClassId, (u64, NodeId)> = HashMap::with_capacity(classes.len());
        // Bellman-Ford style relaxation. Processing in reverse discovery
        // order visits most children before their parents, so this usually
        // converges in two or three sweeps.
        for &class in &classes {
            if let Some(pin) = graph.pinned(class) {
                best.insert(class, (1, pin));
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
                        best.get(&graph.find(child))
                            .map(|&(c, _)| acc.saturating_add(c))
                    });
                    let Some(below) = below else {
                        continue;
                    };
                    let total = own.saturating_add(below.saturating_mul(model.scale(graph, enode)));
                    // On a tie the younger node wins: simplification results
                    // are created after the terms they simplify.
                    let better = best.get(&class).is_none_or(|&(old, old_node)| {
                        total < old || (total == old && enode > old_node)
                    });
                    if better {
                        best.insert(class, (total, enode));
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
        self.best.get(&graph.find(node)).map(|&(c, _)| c)
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
            let &(_, enode) = self.best.get(&class)?;
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
}
