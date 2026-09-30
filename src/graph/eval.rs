//! Reference numeric evaluation of concrete terms.
//!
//! This is the semantics every rule is tested against and the baseline a
//! compiled backend must agree with: walk the term bottom-up and apply each
//! operator's [`EvalFn`](super::op::EvalFn).

use std::collections::HashMap;

use super::id::NodeId;
use super::rule::Env;
use super::store::Graph;

impl Graph {
    /// Evaluates the concrete term `node` to a float under the bindings of
    /// `env`.
    ///
    /// Returns `None` when the term contains an unbound symbol, an operator
    /// without scalar semantics, or a non-numeric literal. Shared subterms
    /// are evaluated once.
    #[must_use]
    pub fn eval(
        &self,
        node: NodeId,
        env: &Env,
    ) -> Option<f64> {
        let mut values: HashMap<NodeId, f64> = HashMap::new();
        let mut stack = vec![(node, false)];
        let mut args: Vec<f64> = Vec::new();
        while let Some((cur, expanded)) = stack.pop() {
            if values.contains_key(&cur) {
                continue;
            }
            let children = self.children(cur);
            if children.is_empty() {
                let value = match (self.as_number(cur), self.as_symbol(cur)) {
                    | (Some(n), _) => n.to_f64(),
                    | (None, Some(s)) => env.value(s)?,
                    | (None, None) => (self.ops().get(self.op(cur)).eval?)(&[]),
                };
                values.insert(cur, value);
            } else if expanded {
                let eval = self.ops().get(self.op(cur)).eval?;
                args.clear();
                for child in children {
                    args.push(*values.get(child)?);
                }
                values.insert(cur, eval(&args));
            } else {
                stack.push((cur, true));
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !values.contains_key(c))
                        .map(|&c| (c, false)),
                );
            }
        }
        values.get(&node).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates_arithmetic() {
        let mut g = Graph::new();
        let node = g.parse("(x + 1)^2 / y - 3").unwrap_or(NodeId::NONE);
        let (x, y) = (g.interner_mut().symbol("x"), g.interner_mut().symbol("y"));
        let mut env = Env::numeric(0.0);
        env.bind(x, 2.0);
        assert_eq!(g.eval(node, &env), None, "y is unbound");
        env.bind(y, 3.0);
        assert_eq!(g.eval(node, &env), Some(0.0));
    }

    #[test]
    fn operators_without_semantics_do_not_evaluate() {
        let mut g = Graph::new();
        let node = g.parse("f(1)").unwrap_or(NodeId::NONE);
        assert_eq!(g.eval(node, &Env::numeric(0.0)), None);
    }

    #[test]
    fn shared_subterms_are_evaluated_once() {
        let mut g = Graph::new();
        // 2^40 paths through the DAG, 41 nodes.
        let mut node = g.sym("x");
        for _ in 0..40 {
            node = g.node(crate::graph::op::core::MUL, &[node, node]);
        }
        let x = g.interner_mut().symbol("x");
        let mut env = Env::numeric(0.0);
        env.bind(x, 1.0);
        assert_eq!(g.eval(node, &env), Some(1.0));
    }
}
