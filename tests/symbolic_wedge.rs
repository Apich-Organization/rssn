//! Alternating-product and exact orientation checks for symbolic differential forms.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::geometry::geometry;

fn wedge(
    left: &str,
    right: &str,
) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    let root = graph.parse(&format!("wedge({left},{right})")).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{left} and {right}");
    Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .map(|value| graph.display(value).to_string())
}

fn basis(
    coefficient: i64,
    indices: &[i64],
) -> String {
    let indices = indices
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    format!("form(list(list({coefficient}, list({indices}))))")
}

#[test]
fn symbolic_wedge_orientation_matches_pairwise_inversion_count() {
    let left = [5, 11, 3, 1, 13, 19, 23, 15];
    let right = [22, 16, 8, 10, 26, 18, 4, 2, 12];
    let mut indices = [left.as_slice(), right.as_slice()].concat();
    let count = indices
        .iter()
        .enumerate()
        .map(|(i, index)| indices[..i].iter().filter(|prior| *prior > index).count())
        .sum::<usize>();
    assert_eq!(count % 2, 1);
    indices.sort_unstable();
    assert_eq!(
        wedge(&basis(2, &left), &basis(3, &right)),
        Some(basis(-6, &indices))
    );
}

#[test]
fn symbolic_wedge_discards_repeated_indices_and_alternates_odd_forms() {
    let a = "form(list(list(2,list(3,1,3)),list(5,list(2))))";
    let b = "form(list(list(7,list(2,4)),list(11,list(1))))";
    assert_eq!(wedge(a, b), Some(basis(-55, &[1, 2])));
    let one_form = "form(list(list(2,list(3)),list(3,list(1)),list(-7,list(2))))";
    assert_eq!(wedge(one_form, one_form).as_deref(), Some("form(list())"));
}

#[test]
fn symbolic_wedge_combines_even_form_products_in_canonical_order() {
    let even = "form(list(list(1,list(4,3)),list(1,list(1,2))))";
    assert_eq!(wedge(even, even), Some(basis(-2, &[1, 2, 3, 4])));
    let left = "form(list(list(1,list(4)),list(2,list(1))))";
    let right = "form(list(list(3,list(3)),list(5,list(2))))";
    assert_eq!(
        wedge(left, right).as_deref(),
        Some(
            "form(list(list(10, list(1, 2)), list(6, list(1, 3)), list(-5, list(2, 4)), list(-3, list(3, 4))))"
        )
    );
}

#[test]
fn symbolic_wedge_accepts_zero_forms_and_checks_the_one_based_index_range() {
    let high = i64::try_from(usize::MAX).unwrap_or(i64::MAX);
    let scalar_form = basis(2, &[]);
    let high_form = basis(-3, &[high, 1]);
    assert_eq!(wedge(&scalar_form, &high_form), Some(basis(6, &[1, high])));
    assert_eq!(wedge(&high_form, &scalar_form), Some(basis(6, &[1, high])));
    assert_eq!(
        wedge("form(list())", &high_form).as_deref(),
        Some("form(list())")
    );
    assert_eq!(wedge(&basis(1, &[0]), &high_form), None);
}
