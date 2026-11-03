//! Numeric kernels that have no symbolic rule set of their own.
//!
//! These are the routines of [`crate::kernels`] that a host otherwise has to
//! reimplement: they take numbers and give numbers, so they belong in the
//! manifest as operators rather than staying library functions.
//!
//! | operator | value |
//! |---|---|
//! | `fftconv(a, b)` | the linear convolution of two sequences, `c[k] = sum_i a[i] b[k - i]`, of length `len(a) + len(b) - 1` |
//! | `interpolate(xs, ys, x)` | the value at `x` of the **Lagrange polynomial** through the points `(xs[i], ys[i])` |
//!
//! Both compute only from literal numbers. A symbolic argument leaves the
//! call as a request: there is no symbolic form of a convolution or of an
//! interpolating polynomial's value, and inventing one (or returning zeros)
//! would be a wrong answer. `interpolate` is polynomial interpolation, not
//! the piecewise-linear interpolant some standard libraries call `interp1`;
//! the two agree only at the nodes, so a host must not substitute one for the
//! other.
//!
//! The other numeric kernels a host may expect are deliberately **not** here:
//! `fft` returns complex values, and a complex tensor has no term encoding
//! ([`crate::api`] refuses a `Value::Complex` with a nonzero imaginary part),
//! so an `fft` operator could not hand its answer back; `odeint` is already an
//! operator of the `ode` rule set.

use crate::graph::Arity;
use crate::graph::Cx;
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
use crate::graph::Tier;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::fftconv;
use crate::kernels::interpolate::lagrange_interpolation;

/// The numeric-kernel rule set.
#[must_use]
pub fn numeric() -> RuleSet {
    RuleSet::new("numeric", install)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let request = |name: &str, arity: u8| {
        OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100)
    };
    let fftconv = i.op(request("fftconv", 2))?;
    let interpolate = i.op(request("interpolate", 3))?;
    // A convolution multiplies the two signals' dimensions; interpolation
    // carries the ordinate's dimension through. Without these the operators
    // would be dimensionally polymorphic by omission, which the coverage test
    // refuses — correctly, because a silently dimensionless convolution is a
    // wrong answer rather than a missing one.
    i.dim(fftconv, DimRule::Product);
    i.dim(interpolate, DimRule::Same(1));
    i.kernel("numeric/fftconv", Tier::Reduce, FftConv { op: fftconv });
    i.kernel(
        "numeric/interpolate",
        Tier::Reduce,
        Interpolate { op: interpolate },
    );
    Ok(())
}

/// The numbers a `list` term holds, in order.
fn numbers(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<f64>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    graph
        .children(node)
        .iter()
        .map(|&item| graph.number_of(item).map(Number::to_f64))
        .collect()
}

/// The single number a term holds.
fn number(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.number_of(node).map(Number::to_f64)
}

/// `fftconv(a, b)`: linear convolution of two real sequences.
struct FftConv {
    op: OpId,
}

impl Kernel for FftConv {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[a, b] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(a), Some(b)) = (numbers(cx.graph, a), numbers(cx.graph, b)) else {
            return Outcome::Pass;
        };
        let out = fftconv::convolve(&a, &b);
        let items: Vec<NodeId> = out.iter().map(|&value| cx.graph.float(value)).collect();
        Outcome::Equal(cx.graph.node(core::LIST, &items))
    }
}

/// `interpolate(xs, ys, x)`: the Lagrange polynomial through the points,
/// evaluated at `x`.
struct Interpolate {
    op: OpId,
}

impl Kernel for Interpolate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[xs, ys, x] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(xs), Some(ys), Some(x)) = (
            numbers(cx.graph, xs),
            numbers(cx.graph, ys),
            number(cx.graph, x),
        ) else {
            return Outcome::Pass;
        };
        if xs.len() != ys.len() || xs.is_empty() {
            return Outcome::Pass;
        }
        let points: Vec<(f64, f64)> = xs.iter().copied().zip(ys).collect();
        // Duplicate abscissae make the interpolant undefined; the call stays
        // a request rather than picking one of the two ordinates.
        let Ok(polynomial) = lagrange_interpolation(&points) else {
            return Outcome::Pass;
        };
        Outcome::Equal(cx.graph.float(polynomial.eval(x)))
    }
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
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    use super::numeric;

    fn run(src: &str) -> String {
        simplify(&[numeric()], src)
    }

    fn vector() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Real),
        }
    }

    #[test]
    fn convolution_is_the_linear_one() {
        // [1, 2, 3] * [0, 1, 1/2] = [0, 1, 5/2, 4, 3/2].
        assert_eq!(
            run("fftconv(list(1, 2, 3), list(0, 1, 0.5))"),
            "list(0, 1, 2.5, 4, 1.5)"
        );
        // The length is len(a) + len(b) - 1, not either input's.
        assert_eq!(run("fftconv(list(1), list(2))"), "list(2)");
        // The empty case is empty, not zeros of some length.
        assert_eq!(run("fftconv(list(), list(1, 2))"), "list()");
    }

    #[test]
    fn lagrange_interpolation_hits_the_nodes_and_the_midpoint() {
        // The parabola through (0, 0), (1, 1), (2, 4) is x^2.
        assert_eq!(
            run("interpolate(list(0, 1, 2), list(0, 1, 4), 3/2)"),
            "2.25"
        );
        assert_eq!(run("interpolate(list(0, 1, 2), list(0, 1, 4), 1)"), "1");
    }

    #[test]
    fn a_symbolic_argument_is_left_as_a_request() {
        let (text, reduced) = reduce_with(
            &[numeric()],
            "interpolate(list(0, 1), list(a, b), 1/2)",
            &[],
        );
        assert!(!reduced, "{text}");
        assert_eq!(text, "interpolate(list(0, 1), list(a, b), 1/2)");
    }

    /// The tensor overload must win over the generic `Term` fallback.
    #[test]
    fn the_tensor_signature_is_the_one_that_resolves() {
        let mut graph = Graph::new();
        let sets: Vec<RuleSet> = vec![numeric()];
        assert!(Engine::install(&mut graph, &sets).is_ok());
        for (name, args) in [
            ("fftconv", vec![vector(), vector()]),
            ("interpolate", vec![vector(), vector(), TypeExpr::Real]),
        ] {
            let op = graph
                .ops()
                .lookup(name)
                .unwrap_or_else(|| panic!("{name} is registered"));
            let (overloads, source) = signatures(graph.ops().get(op));
            assert_eq!(source, SignatureSource::Overlay, "{name}");
            let Best::Unique(chosen) = best_types(&overloads, &args) else {
                panic!("{name}: {args:?} did not resolve to one overload");
            };
            assert_eq!(overloads[chosen].params, args, "{name}");
        }
    }

    /// The path a host takes: resolve the overload from values, run it, and
    /// read the answer back as a value.
    #[test]
    fn the_session_can_delegate_the_numeric_kernels() {
        use crate::api::Session;
        use crate::kernel::CallCtx;
        use crate::kernel::value::TensorValue;
        use crate::kernel::value::Value;

        let session = Session::new();
        let a = Value::Tensor(TensorValue::vector(vec![
            Value::int(1_i32),
            Value::int(2_i32),
            Value::int(3_i32),
        ]));
        let b = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::Float(0.5),
        ]));
        let out = session
            .call_value("fftconv", &[a, b], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(convolution) = out else {
            panic!("expected a tensor, got {out:?}");
        };
        assert_eq!(convolution.shape, vec![5]);
        let values: Vec<f64> = convolution
            .data
            .iter()
            .map(|value| {
                match value {
                    | Value::Int(n) => n.to_string().parse().unwrap_or(f64::NAN),
                    | Value::Real(number) => number.to_f64(),
                    | other => panic!("expected a number, got {other:?}"),
                }
            })
            .collect();
        assert_eq!(values, vec![0.0, 1.0, 2.5, 4.0, 1.5]);

        let xs = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::int(2_i32),
        ]));
        let ys = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::int(4_i32),
        ]));
        let out = session
            .call_value("interpolate", &[xs, ys, Value::Float(1.5)], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Real(answer) = out else {
            panic!("expected a number, got {out:?}");
        };
        assert_eq!(answer.to_f64(), 2.25);
    }
}
