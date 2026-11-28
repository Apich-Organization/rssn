//! Deterministic reverse Cuthill-McKee orderings of sparse patterns.

use rssn::kernels::krylov::Csr;
use rssn::kernels::krylov::SparseOrdering;

#[test]
fn rcm_preserves_bfs_and_degree_ties_across_components() {
    // From vertex 0, the equal-degree farthest nodes are visited as 4, 3.
    // In the second component, the remaining equal-degree leaves are 8, 9.
    // Vertices 5 and 10 are isolated.
    let matrix = Csr::from_triplets(
        11,
        11,
        &[
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 4, 1.0),
            (2, 3, 1.0),
            (7, 6, 1.0),
            (6, 8, 1.0),
            (9, 6, 1.0),
        ],
    );
    assert_eq!(matrix.rcm_ordering(), [10, 9, 8, 6, 7, 5, 3, 2, 0, 1, 4],);
}

#[test]
fn rcm_uses_the_common_square_pattern_and_keeps_isolates() {
    assert!(Csr::from_triplets(0, 3, &[]).rcm_ordering().is_empty());
    assert_eq!(
        Csr::from_triplets(5, 5, &[]).rcm_ordering(),
        [4, 3, 2, 1, 0],
    );
    assert_eq!(
        Csr::from_triplets(4, 6, &[(0, 5, 2.0), (1, 1, 3.0)]).rcm_ordering(),
        [3, 2, 1, 0],
    );
    assert_eq!(
        Csr::from_triplets(4, 6, &[(0, 5, 2.0), (3, 0, 0.0), (3, 3, 3.0)]).rcm_ordering(),
        [2, 1, 3, 0],
    );
    assert_eq!(
        Csr::from_triplets(6, 4, &[(5, 0, 2.0), (0, 3, 1.0), (0, 3, -1.0)]).rcm_ordering(),
        [2, 1, 3, 0],
    );
}

#[test]
fn sparse_lu_uses_the_same_rcm_permutation() {
    let matrix = Csr::from_triplets(
        5,
        5,
        &[
            (0, 0, 3.0),
            (1, 1, 3.0),
            (2, 2, 3.0),
            (3, 3, 3.0),
            (4, 4, 3.0),
            (0, 1, -1.0),
            (0, 2, -1.0),
            (1, 4, -1.0),
            (2, 3, -1.0),
        ],
    );
    let factor = matrix.sparse_lu(SparseOrdering::Rcm).unwrap();
    assert_eq!(factor.permutation(), matrix.rcm_ordering());
    let mut rhs = [0.0; 5];
    matrix.matvec(&[1.0, 2.0, 3.0, 4.0, 5.0], &mut rhs);
    for (actual, expected) in factor.solve(&rhs).iter().zip([1.0, 2.0, 3.0, 4.0, 5.0]) {
        assert!((actual - expected).abs() < 1e-14);
    }
}
