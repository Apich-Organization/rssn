//! Component traversal and deterministic nesting regressions.

use rssn::kernels::ordering::amd_order;
use rssn::kernels::ordering::cholesky_fill;
use rssn::kernels::ordering_ml::NdOptions;
use rssn::kernels::ordering_ml::hybrid_nested_dissection_order;
use rssn::kernels::ordering_ml::multilevel_nested_dissection_with;

fn options(leaf: usize) -> NdOptions {
    NdOptions {
        leaf,
        imbalance: 0.1,
        fm_passes: 4,
    }
}

#[test]
fn isolated_vertices_keep_seed_order_for_every_leaf_policy() {
    for leaf in [0, 1, 16] {
        assert!(multilevel_nested_dissection_with(&[], &options(leaf)).is_empty());
        let isolated = vec![Vec::new(); 17];
        assert_eq!(
            multilevel_nested_dissection_with(&isolated, &options(leaf)),
            (0..17).collect::<Vec<_>>()
        );
    }
    let large = vec![Vec::new(); 16_384];
    let perm = multilevel_nested_dissection_with(&large, &options(32));
    assert_eq!(perm.len(), large.len());
    assert!(perm.into_iter().eq(0..large.len()));
}

#[test]
fn disconnected_parts_remain_contiguous_in_seed_order_with_repeated_neighbours() {
    let graph = [
        vec![6, 4, 6, 0],
        vec![3, 3],
        vec![2],
        vec![1, 1],
        vec![0],
        vec![9, 7],
        vec![0, 0],
        vec![8, 5],
        vec![7, 9, 8],
        vec![5, 8],
    ];
    let expected = [vec![0, 4, 6], vec![1, 3], vec![2], vec![5, 7, 8, 9]];
    for leaf in [0, 1, 3, 64] {
        let perm = multilevel_nested_dissection_with(&graph, &options(leaf));
        assert_eq!(
            perm,
            multilevel_nested_dissection_with(&graph, &options(leaf))
        );
        let mut start = 0;
        for component in &expected {
            let end = start + component.len();
            let mut block = perm[start..end].to_vec();
            block.sort_unstable();
            assert_eq!(&block, component);
            start = end;
        }
        assert_eq!(start, perm.len());
    }
}

#[test]
fn connected_grid_keeps_a_reproducible_permutation_and_fill_score() {
    const WIDTH: usize = 10;
    const HEIGHT: usize = 9;
    let graph: Vec<_> = (0..WIDTH * HEIGHT)
        .map(|index| {
            let (x, y) = (index % WIDTH, index / WIDTH);
            let mut neighbours = Vec::new();
            if y + 1 < HEIGHT {
                neighbours.push(index + WIDTH);
            }
            if x > 0 {
                neighbours.push(index - 1);
            }
            if y > 0 {
                neighbours.push(index - WIDTH);
            }
            if x + 1 < WIDTH {
                neighbours.push(index + 1);
            }
            neighbours
        })
        .collect();
    let perm = multilevel_nested_dissection_with(&graph, &options(8));
    let mut sorted = perm.clone();
    sorted.sort_unstable();
    assert!(sorted.into_iter().eq(0..graph.len()));
    assert_eq!(perm, multilevel_nested_dissection_with(&graph, &options(8)));
    let hybrid = hybrid_nested_dissection_order(&graph);
    assert!(cholesky_fill(&graph, &hybrid) <= cholesky_fill(&graph, &amd_order(&graph)));
}
