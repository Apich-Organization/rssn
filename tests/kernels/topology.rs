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
    // exact death: the diagonal length (the old grid sampling reported 1.5)
    assert!((d[1].intervals[0].death - 2.0_f64.sqrt()).abs() < 1e-12);
}

mod exact_persistence {
    use super::*;
    use rssn::kernels::homology::{distance_matrix, persistence_bars, rips_filtration};

    const INF: f64 = f64::INFINITY;

    fn ivs(d: &PersistenceDiagram) -> Vec<(f64, f64)> {
        d.intervals.iter().map(|i| (i.birth, i.death)).collect()
    }

    fn approx(a: &[(f64, f64)], b: &[(f64, f64)]) -> bool {
        a.len() == b.len()
            && a.iter().zip(b).all(|(x, y)| {
                (x.0 - y.0).abs() < 1e-12 && ((x.1 - y.1).abs() < 1e-12 || (x.1 == INF && y.1 == INF))
            })
    }

    fn circle(n: usize) -> Vec<Vec<f64>> {
        (0..n)
            .map(|k| {
                let a = std::f64::consts::TAU * k as f64 / n as f64;
                vec![a.cos(), a.sin()]
            })
            .collect()
    }

    #[test]
    fn hand_checked_hollow_triangle_then_filled() {
        // vertices at 0; edges at 1, 1, 2; the 2-cell at 3.
        let f = vec![
            (0.0, vec![0]),
            (0.0, vec![1]),
            (0.0, vec![2]),
            (1.0, vec![0, 1]),
            (1.0, vec![1, 2]),
            (2.0, vec![0, 2]),
            (3.0, vec![0, 1, 2]),
        ];
        let d = persistent_homology(&f, 2).unwrap();
        assert!(approx(&ivs(&d[0]), &[(0.0, 1.0), (0.0, 1.0), (0.0, INF)]));
        assert!(approx(&ivs(&d[1]), &[(2.0, 3.0)]));
        assert!(d[2].intervals.is_empty());
    }

    #[test]
    fn hollow_triangle_is_essential_h1() {
        let f = vec![
            (0.0, vec![0]),
            (0.0, vec![1]),
            (0.0, vec![2]),
            (1.0, vec![0, 1]),
            (1.0, vec![1, 2]),
            (2.0, vec![0, 2]),
        ];
        let d = persistent_homology(&f, 1).unwrap();
        assert!(approx(&ivs(&d[1]), &[(2.0, INF)]));
    }

    #[test]
    fn boundary_of_tetrahedron_has_an_essential_h2() {
        let tet = close_complex(&[vec![0, 1, 2, 3]]);
        let mut f: Vec<(f64, Vec<usize>)> = tet.iter().filter(|s| s.len() < 4).map(|s| (s.len() as f64, s.clone())).collect();
        let d = persistent_homology(&f, 2).unwrap();
        assert!(approx(&ivs(&d[2]), &[(3.0, INF)]));
        // every H1 class (an edge's cycle) is killed by the triangles
        assert!(d[1].intervals.iter().all(|i| i.death == 3.0));
        // adding the solid tetrahedron kills it
        f.push((10.0, vec![0, 1, 2, 3]));
        let d = persistent_homology(&f, 2).unwrap();
        assert!(approx(&ivs(&d[2]), &[(3.0, 10.0)]));
    }

    #[test]
    fn two_points() {
        let d = persistent_homology_rips(&[vec![0.0], vec![3.0]], 5.0, 1, 1000).unwrap();
        assert!(approx(&ivs(&d[0]), &[(0.0, 3.0), (0.0, INF)]));
        assert!(d[1].intervals.is_empty());
        // beyond max_epsilon the components never merge
        let d = persistent_homology_rips(&[vec![0.0], vec![3.0]], 2.0, 0, 1000).unwrap();
        assert!(approx(&ivs(&d[0]), &[(0.0, INF), (0.0, INF)]));
    }

    #[test]
    fn square_h1_bar_is_exact() {
        let d = persistent_homology_rips(&square(), 2.0, 2, 10_000).unwrap();
        assert!(approx(&ivs(&d[0]), &[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, INF)]));
        assert!(approx(&ivs(&d[1]), &[(1.0, 2.0_f64.sqrt())]));
        assert!(d[2].intervals.is_empty());
    }

    #[test]
    fn octahedron_has_exact_h2_bar() {
        let pts: Vec<Vec<f64>> = (0..3)
            .flat_map(|a| {
                [1.0, -1.0].map(|sgn| {
                    let mut p = vec![0.0; 3];
                    p[a] = sgn;
                    p
                })
            })
            .collect();
        let d = persistent_homology_rips(&pts, 2.5, 2, 100_000).unwrap();
        let r2 = 2.0_f64.sqrt();
        assert!(approx(&ivs(&d[0]), &[(0.0, r2), (0.0, r2), (0.0, r2), (0.0, r2), (0.0, r2), (0.0, INF)]));
        assert!(d[1].intervals.is_empty(), "triangles enter with the edges: {:?}", d[1]);
        assert!(approx(&ivs(&d[2]), &[(r2, 2.0)]));
    }

    #[test]
    fn circle_sample_has_one_long_h1_bar() {
        let n = 24;
        let d = persistent_homology_rips(&circle(n), 2.0, 1, 1_000_000).unwrap();
        let long: Vec<_> = d[1].intervals.iter().filter(|i| i.death - i.birth > 0.5).collect();
        assert_eq!(long.len(), 1, "{:?}", d[1]);
        let pi = std::f64::consts::PI;
        // born with the nearest-neighbour edges, dead when chords span a third of the circle
        assert!((long[0].birth - 2.0 * (pi / n as f64).sin()).abs() < 1e-12);
        assert!((long[0].death - 2.0 * (pi * (n / 3) as f64 / n as f64).sin()).abs() < 1e-12);
        assert_eq!(d[0].intervals.iter().filter(|i| i.death.is_infinite()).count(), 1);
    }

    #[test]
    fn sphere_sample_has_one_long_h2_bar() {
        // 24 points of a Fibonacci spiral on the unit sphere
        let n = 24usize;
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        let pts: Vec<Vec<f64>> = (0..n)
            .map(|k| {
                let z = 1.0 - 2.0 * (k as f64 + 0.5) / n as f64;
                let r = (1.0 - z * z).sqrt();
                let a = golden * k as f64;
                vec![r * a.cos(), r * a.sin(), z]
            })
            .collect();
        let d = persistent_homology_rips(&pts, 1.9, 2, 10_000_000).unwrap();
        let long: Vec<_> = d[2].intervals.iter().filter(|i| i.death - i.birth > 0.3).collect();
        assert_eq!(long.len(), 1, "{:?}", d[2]);
        assert!(long[0].birth < 1.3 && long[0].death > 1.5);
        assert!(d[1].intervals.iter().all(|i| i.death - i.birth < 0.35), "{:?}", d[1]);
        assert_eq!(d[0].intervals.iter().filter(|i| i.death.is_infinite()).count(), 1);
    }

    #[test]
    fn matches_the_unoptimised_reduction_on_random_clouds() {
        // cross-check clearing against the plain reduction in kernels::homology
        let mut seed = 0x2545_F491_4F6C_DD1D_u64;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for dim in [2usize, 3] {
            let pts: Vec<Vec<f64>> = (0..16).map(|_| (0..dim).map(|_| rnd()).collect()).collect();
            let dist = distance_matrix(&pts);
            let filt = rips_filtration(&dist, 0.9, 3, 1_000_000).unwrap();
            let plain = persistence_bars(&filt, 2).unwrap();
            let fast = persistent_homology(&filt, 2).unwrap();
            for k in 0..=2 {
                let mut want: Vec<(f64, f64)> =
                    plain.iter().filter(|b| b.dim == k).map(|b| (b.birth, b.death)).collect();
                want.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
                assert!(approx(&ivs(&fast[k]), &want), "dim {k}: {:?} vs {want:?}", ivs(&fast[k]));
            }
        }
    }

    #[test]
    fn invalid_filtrations_are_rejected() {
        // edge enters before its vertex
        let f = vec![(0.0, vec![0]), (1.0, vec![1]), (0.5, vec![0, 1])];
        assert_eq!(persistent_homology(&f, 1), Err(PersistenceError::InvalidFiltration));
        // missing face
        let f = vec![(0.0, vec![0]), (1.0, vec![0, 1])];
        assert_eq!(persistent_homology(&f, 1), Err(PersistenceError::InvalidFiltration));
        // too many simplices
        let pts = circle(30);
        assert_eq!(persistent_homology_rips(&pts, 2.0, 2, 100), Err(PersistenceError::TooManySimplices));
    }
}
