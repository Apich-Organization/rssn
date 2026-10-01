//! Substitution on concrete terms.

use std::collections::HashMap;

use super::id::NodeId;
use super::store::Graph;

impl Graph {
    /// Replaces every free occurrence of the symbol `from` in the concrete
    /// term `term` by `to`.
    ///
    /// Occurrences bound by a binding operator (the integration variable
    /// inside an integral's body, say) are left alone. Shared subterms are
    /// rewritten once. `from` must be a symbol leaf; anything else returns
    /// `term` unchanged.
    pub fn substitute(
        &mut self,
        term: NodeId,
        from: NodeId,
        to: NodeId,
    ) -> NodeId {
        let Some(symbol) = self.as_symbol(from) else {
            return term;
        };
        let mut done: HashMap<NodeId, NodeId> = HashMap::new();
        let mut stack = vec![(term, false)];
        while let Some((node, expanded)) = stack.pop() {
            if done.contains_key(&node) {
                continue;
            }
            if self.as_symbol(node) == Some(symbol) {
                done.insert(node, to);
                continue;
            }
            let children = self.children(node).to_vec();
            if children.is_empty() {
                done.insert(node, node);
                continue;
            }
            // Children in which `from` is bound by this node are kept.
            let shielded = |graph: &Self, index: usize| {
                graph.ops().get(graph.op(node)).binder.is_some_and(|b| {
                    let binds = children.get(usize::from(b.var)).and_then(|&v| graph.as_symbol(v)) == Some(symbol);
                    binds && (index == usize::from(b.var) || (index < 32 && b.scope & (1 << index) != 0))
                })
            };
            if !expanded {
                stack.push((node, true));
                for (index, &child) in children.iter().enumerate() {
                    if !shielded(self, index) && !done.contains_key(&child) {
                        stack.push((child, false));
                    }
                }
                continue;
            }
            let rebuilt: Vec<NodeId> = children
                .iter()
                .enumerate()
                .map(|(index, &child)| {
                    if shielded(self, index) { child } else { done.get(&child).copied().unwrap_or(child) }
                })
                .collect();
            let new = if rebuilt == children {
                node
            } else {
                self.try_node(self.op(node), &rebuilt).unwrap_or(node)
            };
            done.insert(node, new);
        }
        done.get(&term).copied().unwrap_or(term)
    }
}

impl Graph {
    /// Replaces every occurrence of the concrete subterm `from` in `term`
    /// by `to`, purely structurally (no regard for binders). Used to turn a
    /// compound subterm into a fresh symbol and back.
    pub fn replace_subterm(
        &mut self,
        term: NodeId,
        from: NodeId,
        to: NodeId,
    ) -> NodeId {
        let mut done: HashMap<NodeId, NodeId> = HashMap::new();
        done.insert(from, to);
        let mut stack = vec![(term, false)];
        while let Some((node, expanded)) = stack.pop() {
            if done.contains_key(&node) {
                continue;
            }
            let children = self.children(node).to_vec();
            if children.is_empty() {
                done.insert(node, node);
                continue;
            }
            if !expanded {
                stack.push((node, true));
                stack.extend(children.iter().filter(|c| !done.contains_key(c)).map(|&c| (c, false)));
                continue;
            }
            let rebuilt: Vec<NodeId> = children.iter().map(|c| done.get(c).copied().unwrap_or(*c)).collect();
            let new =
                if rebuilt == children { node } else { self.try_node(self.op(node), &rebuilt).unwrap_or(node) };
            done.insert(node, new);
        }
        done.get(&term).copied().unwrap_or(term)
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::Arity;
    use crate::graph::Graph;
    use crate::graph::NodeId;
    use crate::graph::OpDescriptor;

    fn subst(
        g: &mut Graph,
        term: &str,
        from: &str,
        to: &str,
    ) -> String {
        let term = g.parse(term).unwrap_or(NodeId::NONE);
        let from = g.sym(from);
        let to = g.parse(to).unwrap_or(NodeId::NONE);
        let out = g.substitute(term, from, to);
        g.display(out)
    }

    #[test]
    fn replaces_free_occurrences() {
        let mut g = Graph::new();
        assert_eq!(subst(&mut g, "x^2 + f(x, y)", "x", "a + 1"), "(a + 1)^2 + f(a + 1, y)");
        assert_eq!(subst(&mut g, "x + y", "z", "1"), "x + y");
        // Building a product merges its literal operands; nothing else is
        // simplified.
        assert_eq!(subst(&mut g, "x * x", "x", "3"), "9");
        assert_eq!(subst(&mut g, "x^2 * x", "x", "3"), "3*3^2", "substitution does not simplify");
    }

    #[test]
    fn respects_binders() {
        let mut g = Graph::new();
        // integral(body, var, lo, hi) binds var inside body only.
        assert!(g.ops_mut().register(OpDescriptor::new("integral", Arity::Fixed(4)).binder(1, 0b1)).is_ok());
        assert_eq!(
            subst(&mut g, "integral(x*t, x, 0, x) + x", "x", "5"),
            "integral(t*x, x, 0, 5) + 5",
            "the bound x and the binder stay; the upper limit and the outer x change"
        );
        assert_eq!(subst(&mut g, "integral(x*t, x, 0, 1)", "t", "x"), "integral(x*x, x, 0, 1)");
    }

    #[test]
    fn subterm_replacement() {
        let mut g = Graph::new();
        let term = g.parse("f(x^2 + 1) * (x^2 + 1)^3 + x").unwrap_or(NodeId::NONE);
        let inner = g.parse("x^2 + 1").unwrap_or(NodeId::NONE);
        let t = g.sym("t");
        let out = g.replace_subterm(term, inner, t);
        assert_eq!(g.display(out), "t^3*f(t) + x");
        assert_eq!(g.replace_subterm(out, t, inner), term);
    }

    #[test]
    fn non_symbol_targets_are_ignored() {
        let mut g = Graph::new();
        let term = g.parse("x + 1").unwrap_or(NodeId::NONE);
        let not_a_symbol = g.int(1);
        let y = g.sym("y");
        assert_eq!(g.substitute(term, not_a_symbol, y), term);
    }

    #[test]
    fn deep_shared_terms() {
        let mut g = Graph::new();
        let mut node = g.sym("x");
        for _ in 0..40 {
            node = g.node(crate::graph::op::core::POW, &[node, node]);
        }
        let (x, y) = (g.sym("x"), g.sym("y"));
        let out = g.substitute(node, x, y);
        assert_ne!(out, node);
        assert_eq!(g.substitute(out, y, x), node);
    }
}
