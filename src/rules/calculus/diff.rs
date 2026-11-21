//! Differentiation kernels.

#[cfg(test)]
mod tests;

use crate::graph::Ball;
use crate::graph::Cx;
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
use crate::graph::id::FastMap;
use crate::graph::op::core;

// Ordinary products keep their expanded shape for downstream rewrites.
const PRODUCT_TREE_THRESHOLD: usize = 64;
const PRODUCT_TREE_MIN_ACTIVE: usize = 5;
const PRODUCT_LEAF_SIZE: usize = 16;

/// Operator attribute: the partial derivatives of an operator.
///
/// One entry per argument, as a pattern over the arguments (`?a` is argument 0, `?b`
/// argument 1, ...); `None` for an argument the operator cannot be
/// differentiated in (the order of a Bessel function), which is fine as
/// long as that argument does not depend on the variable.
///
/// Attaching this to an ordinary operator teaches the differentiation kernel
/// its chain rule. Binding operators require separate differentiation rules
/// that account for their scope and the meaning of the bound expression.
#[derive(Clone, Debug)]
pub struct Partials(pub Vec<Option<Pat>>);

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
        let mut done = FastMap::default();
        // Pat::instantiate needs a mutable graph. Copy the pattern definitions
        // once per operator so identical subcalls can share them.
        let mut partials = FastMap::default();
        let mut children = Vec::new();
        let mut derivs = Vec::new();
        let mut stack = vec![(root, false)];
        while let Some((node, expanded)) = stack.pop() {
            if done.contains_key(&node) {
                continue;
            }
            // Lists and equations are differentiated elementwise, keeping
            // their shape even where an element is constant.
            let op = graph.op(node);
            let elementwise = matches!(op, core::LIST | core::EQ);
            if !elementwise && !graph.depends_on(graph.find(node), x) {
                done.insert(node, zero);
                continue;
            }
            if graph.as_symbol(node) == Some(x) {
                done.insert(node, one);
                continue;
            }
            let chain = if elementwise || matches!(op, core::ADD | core::MUL) {
                None
            } else {
                let attr = partials.entry(op).or_insert_with(|| {
                    let ops = graph.ops();
                    if ops.get(op).binder.is_some() {
                        None
                    } else {
                        ops.attr::<Partials>(op).cloned()
                    }
                });
                let Some(attr) = attr else {
                    done.insert(node, graph.node(self.diff, &[node, var]));
                    continue;
                };
                Some(attr.0.as_slice())
            };
            if !expanded {
                stack.push((node, true));
                stack.extend(
                    graph
                        .children(node)
                        .iter()
                        .filter(|c| !done.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            children.clear();
            children.extend_from_slice(graph.children(node));
            derivs.clear();
            derivs.extend(
                children
                    .iter()
                    .map(|c| done.get(c).copied().unwrap_or(zero)),
            );
            let result = match op {
                | core::LIST | core::EQ => graph.node(op, &derivs),
                | core::ADD => {
                    derivs.retain(|&d| d != zero);
                    sum(graph, &derivs, zero)
                },
                | core::MUL
                    if children.len() > PRODUCT_TREE_THRESHOLD
                        && derivs
                            .iter()
                            .filter(|&&d| d != zero)
                            .take(PRODUCT_TREE_MIN_ACTIVE)
                            .count()
                            == PRODUCT_TREE_MIN_ACTIVE =>
                {
                    product_tree(graph, &children, &derivs, zero, one).1
                },
                | core::MUL => expanded_product_rule(graph, &children, &derivs, zero, one),
                | _ => {
                    chain_rule(graph, &children, &derivs, chain.unwrap_or(&[]), zero, one)
                        .unwrap_or_else(|| graph.node(self.diff, &[node, var]))
                },
            };
            done.insert(node, result);
        }
        done.get(&root).copied().unwrap_or(zero)
    }
}

fn sum(
    graph: &mut Graph,
    terms: &[NodeId],
    zero: NodeId,
) -> NodeId {
    match terms {
        | [] => zero,
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "chain-rule fallback keeps the whole derivative atomic"
)]
fn chain_rule(
    graph: &mut Graph,
    children: &[NodeId],
    derivs: &[NodeId],
    partials: &[Option<Pat>],
    zero: NodeId,
    one: NodeId,
) -> Option<NodeId> {
    // A missing partial for a dependent argument blocks the full chain rule.
    // Adding a known contribution to diff(whole, x) would count it twice.
    let active = || derivs.iter().enumerate().filter(|&(_, &d)| d != zero);
    if active().any(|(i, _)| partials.get(i).is_none_or(Option::is_none)) {
        return None;
    }
    let mut terms = Vec::new();
    for (i, &d) in active() {
        let partial = partials.get(i)?.as_ref()?.instantiate(graph, children)?;
        terms.push(if d == one {
            partial
        } else {
            graph.node(core::MUL, &[partial, d])
        });
    }
    Some(sum(graph, &terms, zero))
}

fn expanded_product_rule(
    graph: &mut Graph,
    children: &[NodeId],
    derivs: &[NodeId],
    zero: NodeId,
    one: NodeId,
) -> NodeId {
    let mut terms = Vec::new();
    let mut factors = Vec::with_capacity(children.len());
    for (i, &d) in derivs.iter().enumerate() {
        if d == zero {
            continue;
        }
        factors.clear();
        factors.extend(
            children
                .iter()
                .enumerate()
                .filter_map(|(j, &child)| (i != j).then_some(child)),
        );
        if d != one {
            factors.push(d);
        }
        terms.push(product(graph, &factors, one));
    }
    sum(graph, &terms, zero)
}

/// Compute a product and its derivative together. Splitting the product saves
/// copying every other factor for each derivative. The graph can flatten small
/// subproducts while larger subproducts remain shared nodes.
fn product_tree(
    graph: &mut Graph,
    children: &[NodeId],
    derivs: &[NodeId],
    zero: NodeId,
    one: NodeId,
) -> (NodeId, NodeId) {
    if children.len() <= PRODUCT_LEAF_SIZE {
        return (
            product(graph, children, one),
            expanded_product_rule(graph, children, derivs, zero, one),
        );
    }
    let mid = children.len() / 2;
    let (left, right) = children.split_at(mid);
    let (left_derivs, right_derivs) = derivs.split_at(mid);
    let (a, da) = product_tree(graph, left, left_derivs, zero, one);
    let (b, db) = product_tree(graph, right, right_derivs, zero, one);
    let derivative = match (da == zero, db == zero) {
        | (true, true) => zero,
        | (false, true) => graph.node(core::MUL, &[da, b]),
        | (true, false) => graph.node(core::MUL, &[a, db]),
        | (false, false) => {
            let first = graph.node(core::MUL, &[da, b]);
            let second = graph.node(core::MUL, &[a, db]);
            graph.node(core::ADD, &[first, second])
        },
    };
    (graph.node(core::MUL, &[a, b]), derivative)
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
        // Already evaluated in this run.
        if cx.graph.approx(cx.graph.find(node)).is_some() {
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
        let mut env = cx.env.clone();
        let f = |t: f64| {
            env.bind(x, t);
            graph.eval(term, &env)
        };
        richardson(f, at).map_or(Outcome::Pass, Outcome::Approx)
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Central differences with two Richardson extrapolation steps. Halving the
/// initial stencil finds a two-sided neighborhood for domain-limited functions.
fn richardson(
    mut f: impl FnMut(f64) -> Option<f64>,
    at: f64,
) -> Option<Ball> {
    if !at.is_finite() {
        return None;
    }
    let mut h = 1e-2 * at.abs().max(1.0);
    let mut difference = None;
    let mut extrapolation = None;
    for _ in 0..32 {
        let upper = at + h;
        let lower = at - h;
        if !(lower < at && at < upper) {
            return None;
        }
        if upper.is_finite() && lower.is_finite() {
            // None means the graph cannot evaluate the term at all. A
            // nonfinite value may just mean this point is outside its domain.
            let right = f(upper)?;
            let left = if right.is_finite() {
                f(lower)?
            } else {
                f64::NAN
            };
            let width = upper - lower;
            let delta = right - left;
            let current = if delta.is_finite() {
                delta / width
            } else {
                right / width - left / width
            };
            if current.is_finite() {
                if let Some(coarse) = difference.replace(current) {
                    // These forms avoid overflowing on a large finite slope.
                    let first = current + (current - coarse) / 3.0;
                    if let Some(coarse) = extrapolation.replace(first) {
                        let best = first + (first - coarse) / 15.0;
                        let error = (best - first).abs();
                        if best.is_finite() && error.is_finite() {
                            // The derivative requires the function to be
                            // defined at the point, as well as around it.
                            return f(at)?.is_finite().then_some(Ball {
                                mid: best,
                                rad: error.max(f64::EPSILON * best.abs()),
                            });
                        }
                    }
                }
                h *= 0.5;
                continue;
            }
        }
        difference = None;
        extrapolation = None;
        h *= 0.5;
    }
    None
}

/// `diffn(f, x, n)`: `n` nested `diff` requests for a literal `n ≥ 0`.
pub(super) struct Repeated {
    pub(super) diffn: OpId,
    pub(super) diff: OpId,
}

impl Kernel for Repeated {
    fn ops(&self) -> Vec<OpId> {
        vec![self.diffn]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[f, x, n] = graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(n) = graph
            .number_of(n)
            .and_then(crate::graph::Number::to_i64)
            .filter(|n| (0..=64).contains(n))
        else {
            return Outcome::Pass;
        };
        let mut term = f;
        for _ in 0..n {
            term = graph.node(self.diff, &[term, x]);
        }
        Outcome::Equal(term)
    }
}
