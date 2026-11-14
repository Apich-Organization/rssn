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
use super::store::Ball;
use super::store::Graph;
use crate::backend::tape::Interval;

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

    /// A rigorous enclosure of `node` over the input boxes `boxes`.
    ///
    /// Every symbol of `node` must appear in `boxes`; anything else makes the
    /// term non-numeric and returns `None`. The enclosure is computed by
    /// lowering the term to the tape IR and interpreting it over
    /// outward-rounded [`Interval`]s, so for every assignment of each symbol
    /// to a point of its box the exact value of the term lies in the returned
    /// ball. It is the interval image of the *float* semantics: an exact
    /// literal is enclosed at its nearest `f64` (widened by the usual
    /// one-ulp outward rounding), not at its rational value.
    ///
    /// Operators without interval rules widen their result to the whole line,
    /// so a ball may be [`Ball::ENTIRE`] and [`Ball::sign`] may return `None`.
    #[must_use]
    pub fn approx_interval(
        &self,
        node: NodeId,
        boxes: &[(SymbolId, Interval)],
    ) -> Option<Ball> {
        let inputs: Vec<SymbolId> = boxes.iter().map(|&(symbol, _)| symbol).collect();
        let values: Vec<Interval> = boxes.iter().map(|&(_, interval)| interval).collect();
        let tape = crate::backend::tape::lower(self, &[node], &inputs, &[]).ok()?;
        let first = tape.enclose(&values, &[]).first().copied()?;
        Ball::from_interval(first.lo, first.hi)
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
            let args: Vec<Complex64> = children
                .iter()
                .map(|c| values.get(c).copied())
                .collect::<Option<_>>()?;
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
                        if exponent.im == 0.0
                            && exponent.re.fract() == 0.0
                            && exponent.re.abs() <= 1024.0
                        {
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
            let value = if value.im == 0.0 {
                Complex64::new(value.re, 0.0)
            } else {
                value
            };
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
        let node = g
            .parse("(x + 1)^2 / y - 3")
            .unwrap_or_else(|e| panic!("{e}"));
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
        let node = g.parse("f(1)").unwrap_or_else(|e| panic!("{e}"));
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

    /// Deterministic sample of the box `[lo, hi]`.
    fn sample(
        lo: f64,
        hi: f64,
        i: usize,
        n: usize,
    ) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64 / n.saturating_sub(1).max(1) as f64;
        lo + t * (hi - lo)
    }

    #[test]
    fn interval_enclosure_contains_every_sampled_value() {
        let mut g = Graph::new();
        let node = g.parse("(x - 1)^2 + 1/3").unwrap_or_else(|e| panic!("{e}"));
        let x = g
            .interner()
            .find_symbol("x")
            .unwrap_or(crate::graph::SymbolId::NONE);
        let ball = g
            .approx_interval(node, &[(x, Interval { lo: 0.9, hi: 1.1 })])
            .unwrap_or_else(|| panic!("enclosure"));
        let (lo, hi) = ball.bounds();
        assert!(lo > 0.0, "the square plus 1/3 is positive on the box");
        let mut env = Env::numeric(0.0);
        env.bind(x, 0.0);
        for i in 0..65 {
            env.bind(x, sample(0.9, 1.1, i, 65));
            let value = g.eval(node, &env).unwrap_or(f64::NAN);
            assert!(
                lo <= value && value <= hi,
                "sample {value} outside [{lo}, {hi}]"
            );
        }
    }

    #[test]
    fn enclosure_respects_the_terms_structure() {
        let mut g = Graph::new();
        let factored = g.parse("(x - 1)^2").unwrap_or_else(|e| panic!("{e}"));
        let expanded = g.parse("x^2 - 2*x + 1").unwrap_or_else(|e| panic!("{e}"));
        let x = g
            .interner()
            .find_symbol("x")
            .unwrap_or(crate::graph::SymbolId::NONE);
        let box_x = [(x, Interval { lo: 0.9, hi: 1.1 })];
        // The same polynomial in a factored and an expanded spelling: the
        // factored one encloses the value tightly, the expanded one loses the
        // correlation between `x^2` and `2*x` and is much wider.
        let tight = g
            .approx_interval(factored, &box_x)
            .unwrap_or_else(|| panic!("enclosure"));
        let wide = g
            .approx_interval(expanded, &box_x)
            .unwrap_or_else(|| panic!("enclosure"));
        let (lo, hi) = tight.bounds();
        assert!(lo > -1e-9 && hi < 0.0101, "factored enclosure [{lo}, {hi}]");
        assert!(
            wide.bounds().1 > hi,
            "expanded form must be wider: {:?} vs {:?}",
            wide.bounds(),
            (lo, hi)
        );
        let mut env = Env::numeric(0.0);
        env.bind(x, 0.0);
        let mut values = Vec::new();
        for i in 0..65 {
            env.bind(x, sample(0.9, 1.1, i, 65));
            values.push(g.eval(factored, &env).unwrap_or(f64::NAN));
        }
        for value in values {
            assert!(lo <= value && value <= wide.bounds().1);
        }
    }

    #[test]
    fn an_operator_without_interval_rules_gives_the_whole_line() {
        let mut g = Graph::new();
        assert!(crate::graph::Engine::install(&mut g, &crate::rules::standard()).is_ok());
        let x = g.interner_mut().symbol("x");
        let node = g.parse("gamma(x) - 1/2").unwrap_or_else(|e| panic!("{e}"));
        let ball = g
            .approx_interval(node, &[(x, Interval { lo: 0.1, hi: 0.9 })])
            .unwrap_or_else(|| panic!("enclosure"));
        assert_eq!(
            ball,
            Ball::ENTIRE,
            "no interval rule for gamma on non-point box"
        );
        assert_eq!(ball.sign(), None);
    }

    #[test]
    fn unbound_symbols_are_reported_not_guessed() {
        let mut g = Graph::new();
        let node = g.parse("x + 1").unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(g.approx_interval(node, &[]), None);
    }
}
