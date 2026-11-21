//! Stable values and derivatives of symbolic AAA approximants.

use rssn::api::Session;
use rssn::backend;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use std::collections::HashSet;

fn approximate(source: &str) -> (Graph, NodeId, NodeId) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::special()]).unwrap();
    let query = graph.parse(source).unwrap();
    engine.run(
        &mut graph,
        &[query],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    let form = Extractor::new(&graph, &[query], &ClosedForm)
        .build(&mut graph, query)
        .unwrap();
    let x = graph.sym("x");
    (graph, form, x)
}

fn support_sites(
    graph: &Graph,
    form: NodeId,
    x: NodeId,
) -> Vec<f64> {
    let mut pending = vec![form];
    let mut visited = HashSet::new();
    let mut sites = Vec::new();
    while let Some(node) = pending.pop() {
        if !visited.insert(node) {
            continue;
        }
        let children = graph.children(node);
        pending.extend_from_slice(children);
        if graph.op(node) == core::ADD
            && let &[a, b] = children
            && (a == x || b == x)
            && let Some(offset) = graph.eval(if a == x { b } else { a }, &Env::numeric(0.0))
        {
            sites.push(-offset);
        }
    }
    sites.sort_by(f64::total_cmp);
    sites.dedup_by(|a, b| a.to_bits() == b.to_bits());
    sites
}

#[test]
fn constant_rational_approximant_is_defined_everywhere() {
    let (graph, form, x) = approximate("rational_approx(7, x, -1, 1)");
    let compiled = backend::compile(&graph, form, &[graph.symbol_of(x).unwrap()]).unwrap();
    let last_site = -(std::f64::consts::PI * 3999.0 / 4000.0).cos();
    for point in [-1.0, last_site, 1.0] {
        assert_eq!(compiled.call(&[point]).to_bits(), 7.0_f64.to_bits());
    }
}

#[test]
fn constant_rational_approximant_returns_a_scalar_through_the_value_api() {
    let session = Session::with_rules(&[rssn::rules::special()]).unwrap();
    let x = session.parse("x").unwrap().to_value().unwrap();
    let result = session
        .call_value(
            "rational_approx",
            &[Value::Float(7.0), x, Value::Float(-1.0), Value::Float(1.0)],
            &CallCtx::new(),
        )
        .unwrap();
    assert_eq!(result.as_f64().unwrap().to_bits(), 7.0_f64.to_bits());
}

#[test]
fn rational_values_and_gradients_remain_finite_at_support_sites() {
    let (graph, form, x) = approximate("rational_approx(3+2*x/(1+x^2), x, -1, 1)");
    let x_symbol = graph.symbol_of(x).unwrap();
    let sites = support_sites(&graph, form, x);
    assert!(sites.len() >= 3);
    let compiled = backend::compile(&graph, form, &[x_symbol]).unwrap();
    let tape = backend::tape::lower(&graph, &[form], &[x_symbol], &[]).unwrap();
    let expected = |x: f64| 3.0 + 2.0 * x / (1.0 + x * x);
    let gradient = |x: f64| 2.0 * (1.0 - x * x) / (1.0 + x * x).powi(2);
    let mut env = Env::numeric(0.0);
    let values = sites
        .into_iter()
        .flat_map(|x| [x.next_down(), x, x.next_up()]);
    for point in values.chain([-1.0, -0.5, 0.0, 0.75, 1.0]) {
        env.bind(x_symbol, point);
        let reference = expected(point);
        let graph_value = graph.eval(form, &env).unwrap();
        assert!(
            (graph_value - reference).abs() < 3e-12,
            "{point}: {graph_value}"
        );
        let compiled_value = compiled.call(&[point]);
        assert!(
            (compiled_value - reference).abs() < 3e-12,
            "{point}: {compiled_value}"
        );
        let derivative = tape.gradient_reverse(0, &[point], &[])[0];
        assert!(
            (derivative - gradient(point)).abs() < 3e-10,
            "{point}: {derivative}"
        );
    }
}

#[test]
fn higher_degree_rational_approximation_is_continuous_at_supports() {
    let (graph, form, x) = approximate("rational_approx(tanh(20*x)+1/(x^2+1/100), x, -1, 1)");
    let x_symbol = graph.symbol_of(x).unwrap();
    let sites = support_sites(&graph, form, x);
    assert!(sites.len() >= 20);
    let mut env = Env::numeric(0.0);
    let points = sites
        .into_iter()
        .flat_map(|x| [x.next_down(), x, x.next_up()]);
    for point in points.chain((0..=128).map(|i| f64::from(i) / 64.0 - 1.0)) {
        env.bind(x_symbol, point);
        let got = graph.eval(form, &env).unwrap();
        let want = (20.0 * point).tanh() + 1.0 / (point * point + 0.01);
        assert!(
            (got - want).abs() < 1e-8 * want.abs().max(1.0),
            "{point}: {got} vs {want}"
        );
    }
}
