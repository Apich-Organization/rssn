//! Exact quadratic, floating-point and unknown symbolic convexity cases.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Payload;
use rssn::graph::Saturate;
use rssn::rules::optimize::optimize;

fn decide(expression: &str) -> Option<bool> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[optimize()]).unwrap();
    let root = graph.parse(expression).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{expression}");
    let value = Extractor::new(&graph, &[root], &ClosedForm).build(&mut graph, root)?;
    let Some(Payload::Bool(result)) = graph.payload(value) else {
        panic!("expected a Boolean decision for {expression}");
    };
    Some(*result)
}

#[test]
fn symbolic_convexity_proves_a_rank_deficient_exact_gram_hessian() {
    let expr = "(a+b+c)^2+(d-f)^2";
    let vars = "list(a,b,c,d,e,f)";
    assert_eq!(decide(&format!("is_convex({expr},{vars})")), Some(true));
    assert_eq!(decide(&format!("is_concave(-({expr}),{vars})")), Some(true));
    assert_eq!(decide(&format!("is_concave({expr},{vars})")), Some(false));
}

#[test]
fn symbolic_convexity_distinguishes_exact_definite_indefinite_and_zero_hessians() {
    // Rational positive diagonal and its negation separate the two decisions.
    assert_eq!(decide("is_convex(x^2/3+y^2/5,list(x,y))"), Some(true));
    assert_eq!(decide("is_concave(-x^2/3-y^2/5,list(x,y))"), Some(true));
    assert_eq!(decide("is_convex(-x^2/3-y^2/5,list(x,y))"), Some(false));
    // x*y has Hessian eigenvalues +1 and -1.
    assert_eq!(decide("is_convex(x*y,list(x,y))"), Some(false));
    assert_eq!(decide("is_concave(x*y,list(x,y))"), Some(false));
    assert_eq!(decide("is_convex(2*x-y,list(x,y))"), Some(true));
    assert_eq!(decide("is_concave(2*x-y,list(x,y))"), Some(true));
}

#[test]
fn symbolic_convexity_preserves_floating_and_nonconstant_behavior() {
    assert_eq!(decide("is_convex(0.25*x^2+y^2,list(x,y))"), Some(true));
    assert_eq!(decide("is_convex(x^3,list(x))"), Some(false));
    assert_eq!(decide("is_convex(x^4,list(x))"), None);
    assert_eq!(decide("is_convex(a*x^2,list(x))"), None);
}

#[test]
fn symbolic_convexity_keeps_variable_list_boundaries() {
    assert_eq!(decide("is_convex(5,list())"), Some(true));
    assert_eq!(decide("is_convex(x^2,list(x,x))"), Some(true));
    assert_eq!(decide("is_convex(x^2,list(x,1))"), None);
    assert_eq!(decide("is_convex(x^2,list(a,b,c,d,e,f,x))"), None);
}
