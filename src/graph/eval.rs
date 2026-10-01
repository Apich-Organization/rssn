//! Reference numeric evaluation of concrete terms.
//!
//! This is the semantics every rule is tested against and the baseline a
//! compiled backend must agree with: walk the term bottom-up and apply each
//! operator's [`EvalFn`](super::op::EvalFn).

use std::collections::HashMap;

use num_complex::Complex64;

use super::id::NodeId;
use super::id::SymbolId;
use super::op::core;
use super::rule::Env;
use super::store::Graph;

/// Operator attribute: the operator's value over the complex numbers.
/// Operators without it are evaluated by [`Graph::eval_complex`] with their
/// real evaluator when all arguments are real.
#[derive(Copy, Clone, Debug)]
pub struct ComplexEval(pub fn(&[Complex64]) -> Option<Complex64>);

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

impl Graph {
    /// Evaluates the concrete term `node` over the complex numbers.
    ///
    /// Operators without a [`ComplexEval`] attribute are evaluated with their
    /// real evaluator when all their arguments are real.
    #[must_use]
    pub fn eval_complex(
        &self,
        node: NodeId,
        bindings: &HashMap<SymbolId, Complex64>,
    ) -> Option<Complex64> {
        let mut values: HashMap<NodeId, Complex64> = HashMap::new();
        let mut stack = vec![(node, false)];
        while let Some((cur, expanded)) = stack.pop() {
            if values.contains_key(&cur) {
                continue;
            }
            let children = self.children(cur);
            if !expanded && !children.is_empty() {
                stack.push((cur, true));
                stack.extend(children.iter().map(|&c| (c, false)));
                continue;
            }
            let args: Vec<Complex64> = children.iter().map(|c| values.get(c).copied()).collect::<Option<_>>()?;
            let value = if let Some(n) = self.as_number(cur) {
                Complex64::new(n.to_f64(), 0.0)
            } else if let Some(s) = self.as_symbol(cur) {
                *bindings.get(&s)?
            } else {
                match self.op(cur) {
                    | core::ADD => args.iter().sum(),
                    | core::MUL => args.iter().product(),
                    | core::POW => {
                        let (base, exponent) = (*args.first()?, *args.get(1)?);
                        if exponent.im == 0.0 && exponent.re.fract() == 0.0 && exponent.re.abs() <= 1024.0 {
                            #[allow(clippy::cast_possible_truncation)]
                            base.powi(exponent.re as i32)
                        } else if base.im == 0.0 && base.re > 0.0 && exponent.im == 0.0 {
                            Complex64::new(base.re.powf(exponent.re), 0.0)
                        } else {
                            base.powc(exponent)
                        }
                    },
                    | op => {
                        if let Some(eval) = self.ops().attr::<ComplexEval>(op) {
                            (eval.0)(&args)?
                        } else if args.iter().all(|a| a.im == 0.0) {
                            let eval = self.ops().get(op).eval?;
                            let reals: Vec<f64> = args.iter().map(|a| a.re).collect();
                            Complex64::new(eval(&reals), 0.0)
                        } else {
                            return None;
                        }
                    },
                }
            };
            // A negative zero imaginary part would select the lower side
            // of branch cuts; values on a cut take the principal side.
            let value = if value.im == 0.0 { Complex64::new(value.re, 0.0) } else { value };
            values.insert(cur, value);
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
