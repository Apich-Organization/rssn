//! Canonical witnesses and exhaustive finite-group automorphism counts.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::discrete::discrete;

#[test]
fn group_isomorphism_uses_table_order_and_preserves_count() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[discrete()]).unwrap();
    let relabeled_c3 = concat!(
        "group(list(20,30,10),",
        "list(list(30,10,20),list(10,20,30),list(20,30,10)))"
    );
    let witness = format!("group_isomorphism(cyclic_group(3),{relabeled_c3})");
    let quaternion = "group_from_perms(list(list(2,5,4,7,6,1,8,3),list(3,8,5,2,7,4,1,6)))";
    let nonisomorphic = format!("group_is_isomorphic(dihedral_group(4),{quaternion})");
    let quaternion_order = format!("group_order({quaternion})");
    for (source, expected) in [
        (witness.as_str(), "list(10, 20, 30)"),
        (
            "group_isomorphism(cyclic_group(1),cyclic_group(1))",
            "list(0)",
        ),
        (
            "group_is_isomorphic(cyclic_group(4),klein_four_group())",
            "false",
        ),
        (quaternion_order.as_str(), "8"),
        (nonisomorphic.as_str(), "false"),
        (
            "group_automorphism_count(group_direct_product(klein_four_group(),klein_four_group()))",
            "20160",
        ),
    ] {
        let root = graph.parse(source).unwrap();
        engine.run(
            &mut graph,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(graph.conflicts().is_empty(), "{source}");
        let value = Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap_or_else(|| panic!("unresolved {source}"));
        assert_eq!(graph.display(value), expected, "{source}");
    }
}
