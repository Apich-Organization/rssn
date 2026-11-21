//! Entry-only symbolic operations on sparse matrices with many empty rows.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::sparse_symbolic::sparse_symbolic;

fn reduce(source: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[sparse_symbolic()]).unwrap();
    let root = graph.parse(source).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{source}");
    assert_eq!(graph.validate(), Ok(()));
    let answer = Extractor::new(&graph, &[root], &ClosedForm).build(&mut graph, root)?;
    Some(graph.display(answer))
}

#[test]
fn sparse_symbolic_rows_handles_an_empty_matrix_at_the_index_limit() {
    let rows = usize::MAX;
    let source = format!("sparse({rows},2,list())");
    assert_eq!(
        reduce(&format!("sp_dims({source})")),
        Some(format!("list({rows}, 2)"))
    );
    assert_eq!(reduce(&format!("sp_nnz({source})")).as_deref(), Some("0"));
}

#[test]
fn sparse_symbolic_rows_preserves_large_columns_and_exact_merged_entries() {
    let max = usize::MAX;
    for rows in [3, max] {
        let source = format!(
            "sparse({rows},{max},list(list({rows},{max},1/3),list(2,7,5),list(1,{max},3/2),list({rows},{max},1/6),list(2,7,-5)))"
        );
        let canonical =
            format!("sparse({rows}, {max}, list(list(1, {max}, 3/2), list({rows}, {max}, 1/2)))");
        assert_eq!(reduce(&source), Some(canonical));
        assert_eq!(
            reduce(&format!("sp_dims({source})")),
            Some(format!("list({rows}, {max})"))
        );
        assert_eq!(reduce(&format!("sp_nnz({source})")).as_deref(), Some("2"));
    }
}

#[test]
fn sparse_symbolic_rows_preserves_symbolic_sums_and_logical_row_order() {
    let rows = usize::MAX;
    let source = format!(
        "sparse({rows},4,list(list({rows},1,x),list(17,2,y),list(1,2,z),list({rows},1,2*x),list(1,1,1/3),list(17,2,-y),list(1,1,1/6)))"
    );
    let canonical =
        format!("sparse({rows}, 4, list(list(1, 1, 1/2), list(1, 2, z), list({rows}, 1, 3*x)))");
    assert_eq!(reduce(&source), Some(canonical));
    assert_eq!(reduce(&format!("sp_nnz({source})")).as_deref(), Some("3"));
}

#[test]
fn sparse_symbolic_rows_retains_input_index_validation() {
    let max = usize::MAX;
    for source in [
        format!("sp_dims(sparse({max},2,list(list(0,1,1))))"),
        format!("sp_nnz(sparse({max},2,list(list({max},3,1))))"),
        format!("sp_dims(sparse(2,{max},list(list(3,{max},1))))"),
    ] {
        assert!(reduce(&source).is_none(), "{source}");
    }
}
