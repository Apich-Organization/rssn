//! Exact determinants across the symbolic and rational interfaces.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use rssn::kernels::qlinalg;
use rssn::kernels::qlinalg::q;
use rssn::rules::linalg::linalg;

fn evaluate(
    graph: &mut Graph,
    source: &str,
) -> NodeId {
    let engine = Engine::install(graph, &[linalg()]).unwrap();
    let root = graph.parse(source).unwrap();
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    Extractor::new(graph, &[root], &ClosedForm)
        .build(graph, root)
        .unwrap()
}

#[test]
fn symbolic_determinants_keep_rational_row_swap_sign_and_cancel_symbolic_minors() {
    let mut graph = Graph::new();
    let exact = evaluate(&mut graph, "det(list(list(0,2/3),list(3/5,4/11)))");
    assert_eq!(
        graph.number_of(exact).unwrap().to_rational(),
        Some(q(-2) / q(5))
    );

    let mut graph = Graph::new();
    let symbolic = evaluate(&mut graph, "det(list(list(x,1,2),list(0,y,3),list(0,0,z)))");
    let expected = graph.parse("x*y*z").unwrap();
    assert_eq!(graph.display(symbolic), graph.display(expected));

    let mut graph = Graph::new();
    let singular = evaluate(&mut graph, "det(list(list(x,1,2),list(x,1,2),list(0,0,z)))");
    assert!(graph.number_of(singular).unwrap().is_zero());
}

#[test]
fn exact_dense_elimination_keeps_the_rank_one_determinant_in_each_pool() {
    let n = 67;
    let m: Vec<Vec<_>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| q(((i % 7 + 1) * (j % 13 + 1) + usize::from(i == j)) as i64))
                .collect()
        })
        .collect();
    // det(I + u v^T) = 1 + v^T u.
    let expected = q(1
        + (0..n)
            .map(|i| ((i % 7 + 1) * (i % 13 + 1)) as i64)
            .sum::<i64>());
    for workers in [1, 2, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        assert_eq!(pool.install(|| qlinalg::det(&m)), expected);
    }
}
