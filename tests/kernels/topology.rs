//! Computational topology kernels.

use rssn::kernels::topology::*;

fn square() -> Vec<Vec<f64>> {
    vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]]
}

#[test]
fn distance_and_components() {
    assert!((euclidean_distance(&[0.0, 0.0], &[3.0, 4.0]) - 5.0).abs() < 1e-12);
    let adj = vec![vec![1], vec![0], vec![3], vec![2], vec![]];
    assert_eq!(find_connected_components(&adj), vec![vec![0, 1], vec![2, 3], vec![4]]);
}

#[test]
fn rank_is_exact() {
    assert_eq!(integer_rank(&[vec![1, 2], vec![2, 4]]), 1);
    assert_eq!(integer_rank(&[vec![1, 0, 1], vec![0, 1, 1], vec![1, 1, 2]]), 2);
    assert_eq!(integer_rank(&[]), 0);
}

#[test]
fn closure_and_boundary() {
    let k = close_complex(&[vec![2, 0, 1]]);
    assert_eq!(k.len(), 7);
    assert_eq!(k[0], vec![0]);
    assert_eq!(k[6], vec![0, 1, 2]);
    assert_eq!(complex_dimension(&k), Some(2));
    assert_eq!(euler_characteristic(&k), 1);
    assert!(verify_boundary_property(&k));
    assert_eq!(simplex_boundary(&[0, 1, 2]).len(), 3);
}

#[test]
fn betti_numbers() {
    let hollow = close_complex(&[vec![0, 1], vec![1, 2], vec![0, 2]]);
    assert_eq!((betti_number(&hollow, 0), betti_number(&hollow, 1)), (1, 1));
    assert_eq!(cohomology_betti_number(&hollow, 1), 1);
    let filled = close_complex(&[vec![0, 1, 2]]);
    assert_eq!(betti_number(&filled, 1), 0);
    let torus = torus_complex(3, 3);
    assert_eq!(euler_characteristic(&torus), 0);
    let b: Vec<usize> = (0..3).map(|k| betti_number(&torus, k)).collect();
    assert_eq!(b, vec![1, 2, 1]);
    let c: Vec<usize> = (0..3).map(|k| cohomology_betti_number(&torus, k)).collect();
    assert_eq!(c, vec![1, 2, 1]);
    assert!(verify_boundary_property(&torus));
    assert!(verify_coboundary_property(&torus));
    let grid = grid_complex(2, 2);
    assert_eq!(betti_number(&grid, 0), 1);
    assert_eq!(betti_number(&grid, 1), 0);
}

#[test]
fn vietoris_rips() {
    assert_eq!(vietoris_rips_complex(&square(), 1.1, 2).len(), 8);
    assert_eq!(betti_numbers_at_radius(&square(), 1.1, 1), vec![1, 1]);
    assert_eq!(betti_numbers_at_radius(&square(), 1.5, 2), vec![1, 0, 1]);
    assert_eq!(betti_numbers_at_radius(&square(), 1.5, 3), vec![1, 0, 0, 0]);
    assert_eq!(betti_numbers_at_radius(&square(), 0.5, 0), vec![4]);
    let f = vietoris_rips_filtration(&square(), 2.0, 2);
    assert_eq!(f.len(), 3);
    assert_eq!(f[0].1.len(), 4);
}

#[test]
fn persistence_of_a_square() {
    let d = compute_persistence(&square(), 1.5, 3, 2);
    assert_eq!(d.len(), 3);
    assert_eq!(d[0].intervals.len(), 4);
    assert_eq!(d[1].intervals.len(), 1);
    assert!((d[1].intervals[0].birth - 1.0).abs() < 1e-12);
    assert!((d[1].intervals[0].death - 1.5).abs() < 1e-12);
}
