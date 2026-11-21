use super::Differentiate;
use super::FiniteDifference;
use super::Partials;
use super::richardson;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::Outcome;
use crate::graph::Pat;
use crate::graph::id::FastMap;
use crate::graph::op::core;

fn setup() -> (Graph, Differentiate) {
    let mut graph = Graph::new();
    let diff = graph
        .ops_mut()
        .register(OpDescriptor::new("diff", Arity::Fixed(2)).flags(OpFlags::HEAVY))
        .unwrap();
    (graph, Differentiate { diff })
}

#[test]
fn missing_partial_keeps_the_entire_request_unevaluated() {
    let (mut graph, kernel) = setup();
    let op = graph
        .ops_mut()
        .register(OpDescriptor::new("partially_known", Arity::Fixed(2)))
        .unwrap();
    graph
        .ops_mut()
        .set_attr(op, Partials(vec![Some(Pat::Var(1)), None]));
    let x = graph.sym("x");
    let y = graph.sym("y");
    let symbol = graph.as_symbol(x).unwrap();
    let body = graph.node(op, &[x, x]);
    let request = graph.node(kernel.diff, &[body, x]);
    assert_eq!(kernel.derive(&mut graph, body, symbol, x), request);
    let body = graph.node(op, &[x, y]);
    assert_eq!(kernel.derive(&mut graph, body, symbol, x), y);

    // Leaving the trailing slot out is equivalent to marking it unknown.
    graph
        .ops_mut()
        .set_attr(op, Partials(vec![Some(Pat::Var(1))]));
    assert_eq!(kernel.derive(&mut graph, body, symbol, x), y);
}

#[test]
fn a_binding_operator_does_not_take_the_ordinary_scalar_chain_rule() {
    let (mut graph, kernel) = setup();
    let op = graph
        .ops_mut()
        .register(OpDescriptor::new("scoped", Arity::Fixed(3)).binder(1, 0b1))
        .unwrap();
    let partial = Some(Pat::Num(crate::graph::Number::from(1)));
    graph
        .ops_mut()
        .set_attr(op, Partials(vec![partial.clone(), None, partial]));
    let x = graph.sym("x");
    let y = graph.sym("y");
    let symbol = graph.as_symbol(x).unwrap();
    let bound_only = graph.node(op, &[x, x, y]);
    let zero = graph.int(0);
    assert_eq!(kernel.derive(&mut graph, bound_only, symbol, x), zero);
    let free_and_bound = graph.node(op, &[x, x, x]);
    let request = graph.node(kernel.diff, &[free_and_bound, x]);
    assert_eq!(
        kernel.derive(&mut graph, free_and_bound, symbol, x),
        request
    );
}

fn dag_edges(
    graph: &Graph,
    root: NodeId,
) -> usize {
    let mut seen = FastMap::default();
    let mut stack = vec![root];
    let mut edges = 0;
    while let Some(node) = stack.pop() {
        if seen.insert(node, ()).is_none() {
            let children = graph.children(node);
            edges += children.len();
            stack.extend_from_slice(children);
        }
    }
    edges
}

#[test]
fn wide_product_derivative_keeps_subproducts_shared() {
    let (mut graph, kernel) = setup();
    let x = graph.sym("x");
    let symbol = graph.as_symbol(x).unwrap();
    let constants: Vec<_> = (0..512).map(|i| 0.5 + f64::from(i) / 16384.0).collect();
    let mut factors: Vec<_> = constants
        .iter()
        .map(|&value| {
            let c = graph.float(value);
            graph.node(core::ADD, &[x, c])
        })
        .collect();
    let body = graph.node(core::MUL, &factors);
    let derivative = kernel.derive(&mut graph, body, symbol, x);
    assert!(dag_edges(&graph, derivative) < 32 * factors.len());

    let mut env = Env::default();
    env.bind(symbol, 0.5);
    let values = || constants.iter().map(|v| v + 0.5);
    let want = values().product::<f64>() * values().map(f64::recip).sum::<f64>();
    let found = graph.eval(derivative, &env).unwrap();
    assert!((found / want - 1.0).abs() < 1e-12);

    // The product rule is defined when factors vanish and must retain
    // multiplicity when the same factor has multiple edges in the graph.
    let minus_half = graph.float(-0.5);
    let zero_factor = graph.node(core::ADD, &[x, minus_half]);
    factors.push(zero_factor);
    let body = graph.node(core::MUL, &factors);
    let derivative = kernel.derive(&mut graph, body, symbol, x);
    let found = graph.eval(derivative, &env).unwrap();
    assert!((found / values().product::<f64>() - 1.0).abs() < 1e-12);
    factors.push(zero_factor);
    let body = graph.node(core::MUL, &factors);
    let derivative = kernel.derive(&mut graph, body, symbol, x);
    assert!(graph.eval(derivative, &env).unwrap().abs() < 1e-12);
}

#[test]
fn richardson_shrinks_to_a_finite_symmetric_stencil() {
    let square_in_domain = |x: f64| Some(narrow_square(&[x]));
    let result = richardson(square_in_domain, 4e-9).unwrap();
    assert!((result.mid / 8e-9 - 1.0).abs() < 1e-12);
    assert!(richardson(square_in_domain, 0.0).is_none());

    let slope = 8e307;
    let result = richardson(|x| Some(slope * x), 0.0).unwrap();
    assert!((result.mid / slope - 1.0).abs() < 1e-12);
    assert!(richardson(Some, f64::MAX).is_none());
    let result = richardson(Some, f64::MAX * 0.999_999).unwrap();
    assert!((result.mid - 1.0).abs() < 1e-12);
    assert!(richardson(|x| Some(if x == 0.0 { f64::NAN } else { x * x }), 0.0).is_none());
}

fn narrow_square(args: &[f64]) -> f64 {
    match args {
        | &[x] if (0.0..1e-8).contains(&x) => x * x,
        | _ => f64::NAN,
    }
}

#[test]
fn numeric_kernel_reuses_the_binding_inside_a_narrow_domain() {
    let (mut graph, kernel) = setup();
    let op = graph
        .ops_mut()
        .register(OpDescriptor::new("narrow_square", Arity::Fixed(1)).eval(narrow_square))
        .unwrap();
    let x = graph.sym("x");
    let body = graph.node(op, &[x]);
    let request = graph.node(kernel.diff, &[body, x]);
    let mut env = Env::numeric(1e-16);
    let symbol = graph.as_symbol(x).unwrap();
    env.bind(symbol, 4e-9);
    let mut context = Cx {
        graph: &mut graph,
        env: &env,
        engine: &Engine::default(),
    };
    assert!(matches!(
        kernel.reduce(&mut context, request),
        Outcome::Pass
    ));
    let Outcome::Approx(answer) =
        (FiniteDifference { diff: kernel.diff }).reduce(&mut context, request)
    else {
        panic!("expected a numeric derivative")
    };
    assert!((answer.mid / 8e-9 - 1.0).abs() < 1e-12);
    assert_eq!(env.value(symbol), Some(4e-9));
}
