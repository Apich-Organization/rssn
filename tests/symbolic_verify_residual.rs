//! Exact residual proofs and sampled verification controls.

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::OpDescriptor;
use rssn::graph::Saturate;
use rssn::rules::verify::verify;

fn verdict(expression: &str) -> Option<String> {
    let mut graph = Graph::new();
    run_verifier(&mut graph, expression)
}

fn run_verifier(
    graph: &mut Graph,
    expression: &str,
) -> Option<String> {
    let engine = Engine::install(graph, &[verify()]).unwrap();
    let root = graph.parse(expression).unwrap();
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{expression}");
    Extractor::new(graph, &[root], &ClosedForm)
        .build(graph, root)
        .map(|node| graph.display(node).to_string())
}

#[test]
fn symbolic_verification_refutes_every_nonzero_exact_derivative_residual() {
    for claim in [
        "verify_derivative(x,x,1+1/100000000000000)",
        "verify_derivative(x,x,1-1/100000000000000)",
        "verify_derivative(x,x,1+1/10^400)",
    ] {
        assert_eq!(verdict(claim).as_deref(), Some("false"), "{claim}");
    }
    assert_eq!(verdict("verify_derivative(x,x,1)").as_deref(), Some("true"));
    assert_eq!(
        verdict("verify_derivative(x/10^400,x,1/10^400)").as_deref(),
        Some("true")
    );
}

#[test]
fn symbolic_verification_checks_exact_solution_and_matrix_residuals() {
    for claim in [
        "verify_solution(x=0,list(x=1/100000000000000))",
        "verify_integral(1,x,x+x/100000000000000)",
        "verify_inverse(list(list(1)),list(list(1+1/100000000000000)))",
    ] {
        assert_eq!(verdict(claim).as_deref(), Some("false"), "{claim}");
    }
    assert_eq!(
        verdict("verify_inverse(list(list(2)),list(list(1/2)))").as_deref(),
        Some("true")
    );
}

#[test]
fn symbolic_verification_keeps_numeric_literal_and_sampling_tolerances() {
    fn check(
        expression: &str,
        expected: &str,
    ) {
        assert_eq!(
            verdict(expression).as_deref(),
            Some(expected),
            "{expression}"
        );
    }
    check("verify_derivative(x,x,1.000000000000001)", "true");
    check("verify_derivative(x,x,1.00001)", "false");
    check("verify_derivative(exp(x),x,exp(x)+1e-14*x)", "true");
    check("verify_derivative(exp(x),x,exp(x)+1e-6*x)", "false");
    check("verify_solution(y(x),list())", "false");
}

#[test]
fn symbolic_verification_short_circuits_serial_refutations() {
    static EVALUATIONS: AtomicUsize = AtomicUsize::new(0);
    fn refute(_: &[f64]) -> f64 {
        EVALUATIONS.fetch_add(1, Ordering::Relaxed);
        1.0
    }

    let depth = 48;
    let expression = format!(
        "verify_solution({}x{},list())",
        "probe(".repeat(depth),
        ")".repeat(depth)
    );
    for (workers, samples) in [(1, 1), (4, 12)] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        EVALUATIONS.store(0, Ordering::Relaxed);
        let verdict = pool.install(|| {
            let mut graph = Graph::new();
            graph
                .ops_mut()
                .register(OpDescriptor::new("probe", Arity::Fixed(1)).eval(refute))
                .unwrap();
            run_verifier(&mut graph, &expression)
        });
        assert_eq!(verdict.as_deref(), Some("false"));
        assert_eq!(EVALUATIONS.load(Ordering::Relaxed), depth * samples);
    }
}

#[test]
fn symbolic_verification_descriptors_show_the_kernel_argument_order() {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &[verify()]).unwrap();
    for (name, signature) in [
        ("verify_derivative", "`verify_derivative(f, x, g)`"),
        ("verify_integral", "`verify_integral(f, x, F)`"),
        (
            "verify_ode_solution",
            "`verify_ode_solution(equation, y(x), solution)`",
        ),
        ("verify_inverse", "`verify_inverse(A, B)`"),
    ] {
        let op = graph.ops().lookup(name).unwrap();
        assert!(graph.ops().get(op).doc.starts_with(signature));
    }
}
