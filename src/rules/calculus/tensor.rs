//! The tensor overload of `sum`.
//!
//! A tensor is a nested `list` in the term graph, so `sum(t)` folds every
//! element and `sum(t, axis)` folds along one axis. The axis is kept with
//! extent 1 (`sum([[1, 2], [3, 4]], 0)` is `[[4, 6]]`), which is the shape
//! convention the host's own reduction uses; dropping the axis would make the
//! rank of the answer disagree with the rank the host's checker expects.
//!
//! The empty cases are the identity of the fold: `sum(list())` is `0`, and a
//! reduction along an axis of extent 0 gives `0` in every surviving cell
//! (`sum(list(), 0)` is `list(0)`). Symbolic elements are summed as terms, so
//! `sum(list(a, b))` is `a + b` and the ordinary arithmetic rules take it from
//! there. A non-tensor argument is left alone: the operator is defined for
//! tensors, and `sum(x)` on a bare symbol is a request, not `x`.
//!
//! The operator is declared by this rule set because `sum` already names the
//! symbolic sum `sum(f, k, a, b)`; the manifest keeps one registration per
//! name (a second one would be emitted as a second declaration by a host that
//! generates bindings from it), so the two meanings are told apart by arity.

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::op::core;

/// The tensor branch of `sum`.
pub(super) struct SumTensor {
    pub(super) sum: OpId,
}

impl Kernel for SumTensor {
    fn ops(&self) -> Vec<OpId> {
        vec![self.sum]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        reduce(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn reduce(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    match args {
        | &[tensor] => {
            let (_, items) = tensor_of(cx.graph, tensor)?;
            Some(add(cx.graph, &items))
        },
        | &[tensor, axis] => {
            let axis = cx.graph.number_of(axis).and_then(Number::to_i64)?;
            let axis = usize::try_from(axis).ok()?;
            let (shape, items) = tensor_of(cx.graph, tensor)?;
            if axis >= shape.len() {
                return None;
            }
            let (out_shape, folded) = fold_axis(cx.graph, &shape, &items, axis);
            let mut cursor = 0;
            Some(nested(cx.graph, &out_shape, &folded, &mut cursor))
        },
        // `sum(f, k, a, b)` belongs to the series kernel.
        | _ => None,
    }
}

/// The full shape of a nested `list` and its elements in row-major order.
///
/// `None` for a ragged argument: `list(1, list(2))` has no shape, and
/// padding it would be a wrong answer rather than a conservative one. A
/// non-list has no shape either — the operator is a tensor reduction.
fn tensor_of(
    graph: &Graph,
    node: NodeId,
) -> Option<(Vec<usize>, Vec<NodeId>)> {
    fn shape_of(
        graph: &Graph,
        node: NodeId,
    ) -> Option<Vec<usize>> {
        if graph.op(node) != core::LIST {
            return Some(Vec::new());
        }
        let items = graph.children(node).to_vec();
        let Some((&first, rest)) = items.split_first() else {
            return Some(vec![0]);
        };
        let tail = shape_of(graph, first)?;
        for &item in rest {
            if shape_of(graph, item)? != tail {
                return None;
            }
        }
        let mut shape = vec![items.len()];
        shape.extend(tail);
        Some(shape)
    }
    fn flatten(
        graph: &Graph,
        node: NodeId,
        items: &mut Vec<NodeId>,
    ) {
        if graph.op(node) != core::LIST {
            items.push(node);
            return;
        }
        for &item in &graph.children(node).to_vec() {
            flatten(graph, item, items);
        }
    }
    if graph.op(node) != core::LIST {
        return None;
    }
    let shape = shape_of(graph, node)?;
    let mut items = Vec::new();
    flatten(graph, node, &mut items);
    (items.len() == shape.iter().product::<usize>()).then_some((shape, items))
}

/// The sum of `terms`, with `0` for none and no wrapper for one.
fn add(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

/// Row-major strides of `shape`.
fn strides(shape: &[usize]) -> Vec<usize> {
    let mut out = vec![1_usize; shape.len()];
    for index in (0..shape.len().saturating_sub(1)).rev() {
        out[index] = out[index + 1].saturating_mul(shape[index + 1]);
    }
    out
}

/// Folds `items` along `axis`, returning the surviving shape and its cells.
fn fold_axis(
    graph: &mut Graph,
    shape: &[usize],
    items: &[NodeId],
    axis: usize,
) -> (Vec<usize>, Vec<NodeId>) {
    let strides = strides(shape);
    let mut out_shape = shape.to_vec();
    out_shape[axis] = 1;
    let count = out_shape.iter().product::<usize>();
    let mut folded = Vec::with_capacity(count);
    for index in 0..count {
        let mut rest = index;
        let mut base = 0_usize;
        // The coordinates of this output cell in the *input* shape: the
        // reduced axis contributes 0, and every other axis its own index.
        let mut coordinates = vec![0_usize; shape.len()];
        for dim in (0..shape.len()).rev() {
            let extent = out_shape[dim].max(1);
            let coordinate = rest % extent;
            rest /= extent;
            coordinates[dim] = coordinate;
            base += coordinate.saturating_mul(strides[dim]);
        }
        let terms: Vec<NodeId> = (0..shape[axis])
            .filter_map(|step| {
                let at = base.checked_add(step.checked_mul(strides[axis])?)?;
                items.get(at).copied()
            })
            .collect();
        folded.push(add(graph, &terms));
    }
    (out_shape, folded)
}

/// Builds the nested `list` of `shape` from row-major `items`.
fn nested(
    graph: &mut Graph,
    shape: &[usize],
    items: &[NodeId],
    cursor: &mut usize,
) -> NodeId {
    let Some((&extent, rest)) = shape.split_first() else {
        let node = items.get(*cursor).copied().unwrap_or_else(|| graph.int(0));
        *cursor += 1;
        return node;
    };
    let children: Vec<NodeId> = (0..extent)
        .map(|_| nested(graph, rest, items, cursor))
        .collect();
    graph.node(core::LIST, &children)
}

#[cfg(test)]
mod tests {
    use crate::graph::Engine;
    use crate::graph::Graph;
    use crate::graph::RuleSet;
    use crate::kernel::SignatureSource;
    use crate::kernel::infer::signatures;
    use crate::kernel::overload::Best;
    use crate::kernel::overload::best_types;
    use crate::kernel::signature::NatExpr;
    use crate::kernel::signature::TypeExpr;
    use crate::rules::calculus::calculus;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[calculus()], src)
    }

    fn tensor(rank: u64) -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(rank),
            shape: None,
            elem: Box::new(TypeExpr::Real),
        }
    }

    #[test]
    fn a_tensor_is_folded_over_every_element() {
        assert_eq!(run("sum(list(1, 2, 3))"), "6");
        assert_eq!(run("sum(list(list(1, 2), list(3, 4)))"), "10");
    }

    #[test]
    fn an_axis_is_folded_and_kept_with_extent_one() {
        assert_eq!(
            run("sum(list(list(1, 2), list(3, 4)), 0)"),
            "list(list(4, 6))"
        );
        assert_eq!(
            run("sum(list(list(1, 2), list(3, 4)), 1)"),
            "list(list(3), list(7))"
        );
    }

    #[test]
    fn the_empty_tensor_sums_to_the_identity() {
        assert_eq!(run("sum(list())"), "0");
        assert_eq!(run("sum(list(), 0)"), "list(0)");
        // An axis of extent 0 folds to 0 in every cell that survives.
        assert_eq!(
            run("sum(list(list(), list()), 1)"),
            "list(list(0), list(0))"
        );
    }

    #[test]
    fn symbolic_elements_are_summed_as_terms() {
        assert_eq!(run("sum(list(a, b))"), "a + b");
        assert_eq!(
            run("sum(list(list(a, b), list(c, d)), 0)"),
            "list(list(a + c, b + d))"
        );
        // A bare symbol is not a tensor: the request stands.
        let (text, reduced) = reduce_with(&[calculus()], "sum(a)", &[]);
        assert!(!reduced, "{text}");
        assert_eq!(text, "sum(a)");
    }

    #[test]
    fn an_out_of_range_axis_is_left_as_a_request() {
        let (text, reduced) = reduce_with(&[calculus()], "sum(list(1, 2), 3)", &[]);
        assert!(!reduced, "{text}");
        assert_eq!(text, "sum(list(1, 2), 3)");
    }

    /// The tensor overload must be the one a tensor argument resolves to, not
    /// the generic `Term` fallback that would accept anything.
    #[test]
    fn a_tensor_argument_resolves_to_the_tensor_overload() {
        let mut graph = Graph::new();
        let sets: Vec<RuleSet> = vec![calculus()];
        assert!(Engine::install(&mut graph, &sets).is_ok());
        let op = graph
            .ops()
            .lookup("sum")
            .unwrap_or_else(|| panic!("sum is registered"));
        let (overloads, source) = signatures(graph.ops().get(op));
        assert_eq!(source, SignatureSource::Overlay);
        let Best::Unique(chosen) = best_types(&overloads, &[tensor(1)]) else {
            panic!("a rank-1 tensor did not resolve to one overload");
        };
        assert!(
            is_tensor(&overloads[chosen].params[0]),
            "{:?}",
            overloads[chosen]
        );
        assert_eq!(overloads[chosen].ret, TypeExpr::Real);

        let Best::Unique(chosen) = best_types(&overloads, &[tensor(2), TypeExpr::Int]) else {
            panic!("a rank-2 tensor and an axis did not resolve to one overload");
        };
        assert!(
            is_tensor(&overloads[chosen].params[0]),
            "{:?}",
            overloads[chosen]
        );
        assert_eq!(overloads[chosen].params[1], TypeExpr::Int);
        assert!(is_tensor(&overloads[chosen].ret), "{:?}", overloads[chosen]);
    }

    /// A tensor parameter with an unpinned rank, which is what the overlay
    /// declares; a literal rank would refuse a rank-2 argument.
    fn is_tensor(ty: &TypeExpr) -> bool {
        matches!(
            ty,
            TypeExpr::Tensor {
                rank: NatExpr::Var(_),
                shape: None,
                elem,
            } if **elem == TypeExpr::Real
        )
    }

    /// The path a host takes: resolve the overload from values, run it, and
    /// read the answer back as a value.
    #[test]
    fn the_session_can_delegate_a_tensor_sum() {
        use crate::api::Session;
        use crate::kernel::CallCtx;
        use crate::kernel::value::TensorValue;
        use crate::kernel::value::Value;

        let session = Session::new();
        let matrix = Value::Tensor(TensorValue {
            shape: vec![2, 2],
            elem: crate::kernel::signature::TypeExpr::Int,
            data: vec![
                Value::int(1_i32),
                Value::int(2_i32),
                Value::int(3_i32),
                Value::int(4_i32),
            ],
        });
        let total = session
            .call_value("sum", std::slice::from_ref(&matrix), &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(total, Value::int(10_i32));

        let columns = session
            .call_value("sum", &[matrix, Value::int(0_i32)], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(columns) = columns else {
            panic!("expected a tensor, got {columns:?}");
        };
        assert_eq!(columns.shape, vec![1, 2]);
        assert_eq!(columns.data, vec![Value::int(4_i32), Value::int(6_i32)]);
    }
}
