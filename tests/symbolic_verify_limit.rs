//! Exact, signed infinite and one-sided claims through the symbolic limit API.

use rssn::Session;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::kernel::CallCtx;
use rssn::kernel::TypeExpr;
use rssn::kernel::Value;
use rssn::rules::verify::verify;

fn extract(expression: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[verify()]).unwrap();
    let root = graph.parse(expression).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{expression}");
    Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .map(|node| graph.display(node))
}

fn claims(
    inputs: &[&str],
    expected: &str,
) {
    for &input in inputs {
        assert_eq!(extract(input).as_deref(), Some(expected), "{input}");
    }
}

#[test]
fn verifies_signed_divergence_at_infinite_and_finite_points() {
    claims(
        &[
            "verify_limit(x,x,oo,oo)",
            "verify_limit(x,x,-oo,-oo)",
            "verify_limit(-x,x,-oo,oo)",
            "verify_limit(1/x^2,x,0,oo)",
            "verify_limit(-1/x^2,x,0,-oo)",
            "verify_limit(1/x,x,0,oo,plus)",
            "verify_limit(1/x,x,0,-oo,minus)",
            "verify_limit(1/(x-10^400),x,10^400,oo,plus)",
            "verify_limit(1/(x-10^400),x,10^400,-oo,minus)",
            "verify_limit(1/(x+10^400)^2,x,-10^400,oo)",
        ],
        "true",
    );
    claims(
        &[
            "verify_limit(x,x,oo,-oo)",
            "verify_limit(-1/x^2,x,0,oo)",
            "verify_limit(1/x,x,0,oo)",
            "verify_limit(1/x,x,0,oo,minus)",
            "verify_limit(10^400,x,0,oo)",
            "verify_limit(x,x,oo,10^400)",
            "verify_limit(x,x,10^400,oo)",
            "verify_limit(1/(x-10^400),x,10^400,oo)",
        ],
        "false",
    );
}

#[test]
fn verifies_finite_one_sided_limits_only_on_the_requested_side() {
    claims(
        &[
            "verify_limit(abs(x)/x,x,0,1,plus)",
            "verify_limit(abs(x)/x,x,0,-1,minus)",
            "verify_limit(sqrt(x),x,0,0,plus)",
            "verify_limit(sqrt(x),x,0,0)",
            "verify_limit(sqrt(x),x,1/10^12,1/10^6,minus)",
        ],
        "true",
    );
    claims(
        &[
            "verify_limit(abs(x)/x,x,0,-1,plus)",
            "verify_limit(abs(x)/x,x,0,1)",
            "verify_limit(abs(x)/x,x,0,1,invalid)",
            "verify_limit(sqrt(x),x,0,0,minus)",
        ],
        "false",
    );
    assert_eq!(extract("limit(abs(x)/x,x,0,plus)").as_deref(), Some("1"));
    assert_eq!(extract("limit(abs(x)/x,x,0,minus)").as_deref(), Some("-1"));
    assert_eq!(extract("limit(sqrt(x),x,0,minus)"), None);
    assert_eq!(
        extract("not(verify_limit(x,x,0,1))").as_deref(),
        Some("true")
    );
}

#[test]
fn verifies_exact_finite_values_without_tolerating_nearby_exact_claims() {
    claims(
        &[
            "verify_limit(x,x,1/3,1/3)",
            "verify_limit(x,x,1/3,0.3333333333333333)",
            "verify_limit(x,x,a,a)",
            "verify_limit(10^400,x,0,10^400)",
            "verify_limit(x,x,10^400,10^400)",
            "verify_limit(x/(x+1),x,10^400,10^400/(10^400+1))",
            "verify_limit(sin(x)/x,x,0,1)",
            "verify_limit((1+1/x)^x,x,oo,E)",
        ],
        "true",
    );
    claims(
        &[
            "verify_limit(x,x,0,1/1000)",
            "verify_limit(x,x,1/3,1/3+1/10^400)",
            "verify_limit(10^400,x,0,10^400+1)",
            "verify_limit(x/(x+1),x,10^400,1)",
        ],
        "false",
    );
}

#[test]
fn limit_verification_declares_and_dispatches_boolean_overloads() {
    let session = Session::with_rules(&[verify()]).unwrap();
    let overloads = session.overloads("verify_limit");
    assert_eq!(overloads.len(), 2);
    for (index, expressions) in [
        &["x", "x", "0", "0"][..],
        &["abs(x)/x", "x", "0", "-1", "minus"][..],
    ]
    .into_iter()
    .enumerate()
    {
        let (key, signature) = &overloads[index];
        assert_eq!(signature.params.len(), expressions.len());
        assert_eq!(signature.params[1], TypeExpr::Symbol);
        assert_eq!(signature.ret, TypeExpr::Bool);
        let arguments: Vec<_> = expressions
            .iter()
            .map(|text| session.parse(text).unwrap().to_value().unwrap())
            .collect();
        assert_eq!(
            session
                .resolve_overload("verify_limit", &arguments)
                .unwrap(),
            *key
        );
        assert_eq!(
            session.call_key(key, &arguments, &CallCtx::new()).unwrap(),
            Value::Bool(true)
        );
    }
    let arguments =
        ["x", "x", "1+1/10^14"].map(|text| session.parse(text).unwrap().to_value().unwrap());
    assert_eq!(
        session
            .call_value("verify_derivative", &arguments, &CallCtx::new())
            .unwrap(),
        Value::Bool(false)
    );
}
