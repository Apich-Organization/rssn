//! Differentiation kernels.

use std::collections::HashMap;

use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Extractor;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Pat;
use crate::graph::SizeCost;
use crate::graph::SymbolId;
use crate::graph::op::core;

/// Operator attribute: the partial derivative with respect to each
/// argument, as a pattern over the arguments (`?a` is argument 0, `?b`
/// argument 1, ...).
///
/// Attaching this to an operator is all it takes to teach the
/// differentiation kernel about it.
#[derive(Clone, Debug)]
pub struct Partials(pub Vec<Pat>);

/// Structural differentiation: sum rule, n-ary product rule, and the chain
/// rule through every operator that has [`Partials`].
///
/// Whatever cannot be differentiated — an undetermined function, an
/// operator without partials — is left as an inner `diff` node, so the
/// result is always an identity and partial progress is kept.
pub(super) struct Differentiate {
    pub(super) diff: OpId,
}

impl Kernel for Differentiate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.diff]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[body, var] = graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(x) = graph.symbol_of(var) else {
            return Outcome::Pass;
        };
        // Differentiate the best known form of the body, not the spelling
        // the request happened to be built with.
        let Some(term) = Extractor::new(graph, &[body], &SizeCost).build(graph, body) else {
            return Outcome::Pass;
        };
        let derivative = self.derive(graph, term, x, var);
        if graph.same(derivative, node) && graph.op(derivative) == self.diff {
            Outcome::Pass
        } else {
            Outcome::Equal(derivative)
        }
    }

    fn revisit(&self) -> bool {
        // An inner request may get reduced later and unblock this one.
        true
    }
}

impl Differentiate {
    /// Derivative of the concrete term `root` with respect to `x`.
    pub(super) fn derive(
        &self,
        graph: &mut Graph,
        root: NodeId,
        x: SymbolId,
        var: NodeId,
    ) -> NodeId {
        let zero = graph.int(0);
        let one = graph.int(1);
        let mut done: HashMap<NodeId, NodeId> = HashMap::new();
        let mut stack = vec![(root, false)];
        while let Some((node, expanded)) = stack.pop() {
            if done.contains_key(&node) {
                continue;
            }
            if !graph.depends_on(graph.find(node), x) {
                done.insert(node, zero);
                continue;
            }
            if graph.as_symbol(node) == Some(x) {
                done.insert(node, one);
                continue;
            }
            let op = graph.op(node);
            let structural = op == core::ADD
                || op == core::MUL
                || graph.ops().attr::<super::Partials>(op).is_some();
            if !structural {
                let unevaluated = graph.node(self.diff, &[node, var]);
                done.insert(node, unevaluated);
                continue;
            }
            let children = graph.children(node).to_vec();
            if !expanded {
                stack.push((node, true));
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !done.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            let derivs: Vec<NodeId> = children
                .iter()
                .map(|c| done.get(c).copied().unwrap_or(zero))
                .collect();
            let mut terms: Vec<NodeId> = Vec::new();
            if op == core::ADD {
                terms.extend(derivs.iter().copied().filter(|&d| d != zero));
            } else if op == core::MUL {
                for (i, &d) in derivs.iter().enumerate() {
                    if d == zero {
                        continue;
                    }
                    let mut factors: Vec<NodeId> = children
                        .iter()
                        .enumerate()
                        .filter(|&(j, _)| j != i)
                        .map(|(_, &c)| c)
                        .collect();
                    if d != one {
                        factors.push(d);
                    }
                    terms.push(product(graph, &factors, one));
                }
            } else {
                let partials = graph
                    .ops()
                    .attr::<super::Partials>(op)
                    .map(|p| p.0.clone())
                    .unwrap_or_default();
                for (i, &d) in derivs.iter().enumerate() {
                    if d == zero {
                        continue;
                    }
                    let partial = partials
                        .get(i)
                        .and_then(|p| p.instantiate(graph, &children));
                    let term = match partial {
                        | Some(p) if d == one => p,
                        | Some(p) => graph.node(core::MUL, &[p, d]),
                        // No rule for this argument: leave the whole node.
                        | None => graph.node(self.diff, &[node, var]),
                    };
                    terms.push(term);
                    if partial.is_none() {
                        break;
                    }
                }
                if partials.len() < derivs.len() {
                    terms = vec![graph.node(self.diff, &[node, var])];
                }
            }
            let result = match terms.as_slice() {
                | [] => zero,
                | [only] => *only,
                | _ => graph.node(core::ADD, &terms),
            };
            done.insert(node, result);
        }
        done.get(&root).copied().unwrap_or(zero)
    }
}

fn product(
    graph: &mut Graph,
    factors: &[NodeId],
    one: NodeId,
) -> NodeId {
    match factors {
        | [] => one,
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

/// Numeric differentiation by Richardson-extrapolated central differences.
///
/// Fires only when a numeric answer is wanted, the variable is bound, the
/// request has no closed form, and the body can be evaluated — typically a
/// body built from operators that have scalar semantics but no partials.
pub(super) struct FiniteDifference {
    pub(super) diff: OpId,
}

impl Kernel for FiniteDifference {
    fn ops(&self) -> Vec<OpId> {
        vec![self.diff]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let &[body, var] = graph.children(node) else {
            return Outcome::Pass;
        };
        let class = graph.find(node);
        let heavy = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).flags.has(OpFlags::HEAVY);
        if graph.enodes(class).any(|n| !heavy(graph, n)) {
            // The structural kernel got there first.
            return Outcome::Pass;
        }
        let (Some(x), Some(term)) = (
            graph.symbol_of(var),
            Extractor::new(graph, &[body], &SizeCost).build(graph, body),
        ) else {
            return Outcome::Pass;
        };
        let Some(at) = cx.env.value(x) else {
            return Outcome::Pass;
        };
        let f = |t: f64| {
            let mut env: Env = cx.env.clone();
            env.bind(x, t);
            graph.eval(term, &env)
        };
        richardson(f, at).map_or(Outcome::Pass, Outcome::Approx)
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Central differences with two Richardson extrapolation steps; the error
/// estimate is the difference between the last two extrapolants.
fn richardson(
    f: impl Fn(f64) -> Option<f64>,
    at: f64,
) -> Option<Ball> {
    let h0 = 1e-2 * at.abs().max(1.0);
    let central = |h: f64| Some((f(at + h)? - f(at - h)?) / (2.0 * h));
    let (d1, d2, d3) = (central(h0)?, central(h0 / 2.0)?, central(h0 / 4.0)?);
    let (e1, e2) = ((4.0 * d2 - d1) / 3.0, (4.0 * d3 - d2) / 3.0);
    let best = (16.0 * e2 - e1) / 15.0;
    best.is_finite().then(|| Ball {
        mid: best,
        rad: (best - e2).abs().max(f64::EPSILON * best.abs()),
    })
}
