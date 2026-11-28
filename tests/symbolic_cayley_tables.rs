//! Deterministic Cayley tables through the symbolic permutation-group operators.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::discrete::discrete;

#[test]
fn symmetric_group_has_identical_public_cayley_tables_across_worker_counts() {
    let build = || {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &[discrete()]).unwrap();
        let source = "group_from_perms(list(list(2,3,4,5,1),list(2,1,3,4,5)))";
        let root = graph.parse(source).unwrap();
        engine.run(
            &mut graph,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(graph.conflicts().is_empty());
        let result = Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap();
        assert_eq!(graph.ops().get(graph.op(result)).name.as_ref(), "group");
        let &[elements, table] = graph.children(result) else {
            panic!("missing elements or table");
        };
        assert_eq!(graph.children(elements).len(), 120);
        let rows = graph.children(table);
        assert_eq!(rows.len(), 120);
        assert!(rows.iter().all(|&row| graph.children(row).len() == 120));
        graph.display(result)
    };
    let run = |workers| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(build)
    };
    let serial = run(1);
    assert!(serial == run(4));
}
