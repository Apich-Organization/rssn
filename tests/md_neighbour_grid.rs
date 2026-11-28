//! Regression tests for dense and sparse Verlet neighbour bins.

use rssn::sim::md_engine::Boundary;
use rssn::sim::md_engine::ForceField;
use rssn::sim::md_engine::MdSystem;
use rssn::sim::md_engine::NeighbourList;
use rssn::sim::md_engine::PairTerm;
use rssn::sim::md_engine::SimBox;
use rssn::sim::md_engine::compute_forces;

fn brute_orthorhombic<const D: usize>(
    sys: &MdSystem<D>,
    radius: f64,
) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    for (i, p) in sys.pos.iter().enumerate() {
        for (j, q) in sys.pos.iter().enumerate().skip(i + 1) {
            let distance2 = (0..D)
                .map(|a| {
                    let absolute = (p[a] - q[a]).abs();
                    let d = if sys.bx.bc[a] == Boundary::Periodic {
                        let l = sys.bx.len[a];
                        let wrapped = absolute % l;
                        wrapped.min(l - wrapped)
                    } else {
                        absolute
                    };
                    d * d
                })
                .sum::<f64>();
            if distance2 < radius * radius {
                result.push((i, j));
            }
        }
    }
    result
}

#[test]
fn md_neighbour_grid_huge_periodic_boxes_do_not_allocate_empty_cells() {
    for length in [1.0e7, 2.0_f64.powi(80)] {
        let positions = vec![
            [0.125, 0.25, 0.375],
            [length - 0.125, 0.25, 0.375],
            [0.375, 0.25, 0.375],
            [1.125, 0.25, 0.375],
            [length / 2.0, 0.25, 0.375],
        ];
        let sys = MdSystem::new(positions, SimBox::periodic([length; 3]));
        let mut list = NeighbourList::new(1.0, 0.0);
        list.build(&sys);
        assert!(list.pairs.contains(&(0, 2)));
        assert!(!list.pairs.contains(&(0, 3)), "strict cutoff");
        list.pairs.sort_unstable();
        assert_eq!(list.pairs, brute_orthorhombic(&sys, 1.0));
    }
}

#[test]
fn md_neighbour_grid_wraps_translated_positions_and_preserves_forces() {
    for lengths in [[10.0, 2.0, 2.0], [10.0; 3]] {
        let wrapped = vec![
            [5.125, 0.25, 0.25],
            [5.375, 0.25, 0.25],
            [0.125, 0.25, 0.25],
            [9.625, 0.25, 0.25],
            [5.125, 0.75, 0.25],
        ];
        let mut translated = wrapped.clone();
        translated[0][0] += 10.0;
        translated[3][0] -= 20.0;
        translated[4][1] += lengths[1];
        let a = MdSystem::new(wrapped, SimBox::periodic(lengths));
        let b = MdSystem::new(translated.clone(), SimBox::periodic(lengths));
        let ff = ForceField::new(
            vec![PairTerm::Morse {
                d: 2.0,
                alpha: 1.0,
                r0: 0.5,
            }],
            0.75,
        );
        let mut list_a = NeighbourList::new(0.75, 0.25);
        let mut list_b = NeighbourList::new(0.75, 0.25);
        let fa = compute_forces(&a, &ff, &mut list_a).unwrap();
        let fb = compute_forces(&b, &ff, &mut list_b).unwrap();
        assert_eq!(list_a.pairs, list_b.pairs);
        assert_eq!(fa.epot.to_bits(), fb.epot.to_bits());
        assert_eq!(fa.f, fb.f);
        assert_eq!(fa.virial, fb.virial);
        assert!(fb.epot < 0.0);
        assert_eq!(b.pos, translated);
        list_b.pairs.sort_unstable();
        assert_eq!(list_b.pairs, brute_orthorhombic(&b, 1.0));
    }
}

#[test]
fn md_neighbour_grid_mixed_boundaries_distinguish_periodic_and_open_axes() {
    let bx = SimBox::new(
        [12.0, 9.0, 1.0],
        [Boundary::Periodic, Boundary::Reflecting, Boundary::Open],
    );
    let sys = MdSystem::new(
        vec![
            [0.25, 0.25, -500.5],
            [11.75, 0.25, -500.5],
            [11.75, 8.75, -500.5],
            [0.25, 0.25, 500.5],
            [0.25, 0.25, 501.0],
            [0.25, 0.25, -499.5],
        ],
        bx,
    );
    let mut list = NeighbourList::new(0.75, 0.25);
    list.build(&sys);
    list.pairs.sort_unstable();
    assert_eq!(list.pairs, vec![(0, 1), (3, 4)]);
    assert_eq!(list.pairs, brute_orthorhombic(&sys, 1.0));
}

#[test]
fn md_neighbour_grid_keeps_the_same_cell_and_particle_order() {
    let positions = vec![
        [6.25, 0.25],
        [0.25, 0.25],
        [0.75, 0.25],
        [6.75, 0.25],
        [0.25, 1.125],
        [5.75, 0.25],
        [0.25, 2.0],
    ];
    // There are 28 cells in the dense case and 105 in the sparse case.
    // Both visit the occupied cells and their wrapping stencil in the same order.
    for height in [4.0, 15.0] {
        let sys = MdSystem::new(positions.clone(), SimBox::periodic([7.0, height]));
        let mut list = NeighbourList::new(1.0, 0.0);
        list.build(&sys);
        assert_eq!(
            list.pairs,
            vec![(1, 3), (1, 2), (1, 4), (0, 5), (0, 3), (4, 6)]
        );
        list.pairs.sort_unstable();
        assert_eq!(list.pairs, brute_orthorhombic(&sys, 1.0));
    }
}

#[test]
fn md_neighbour_grid_sparse_triclinic_matches_explicit_lattice_images() {
    let cell = [[19.0, 0.0, 0.0], [3.0, 17.0, 0.0], [2.0, 1.0, 23.0]];
    let bx = SimBox::triclinic(cell).unwrap();
    let fractions = [
        [0.01, 0.02, 0.03],
        [0.99, 0.02, 0.03],
        [0.01, 0.98, 0.03],
        [0.01, 0.02, 0.99],
        [0.5, 0.5, 0.5],
        [0.5, 0.52, 0.55],
    ];
    let sys = MdSystem::new(
        fractions.iter().map(|s| bx.from_fractional(s)).collect(),
        bx,
    );
    let mut list = NeighbourList::new(0.875, 0.375);
    list.build(&sys);
    let mut expected = Vec::new();
    for (i, p) in sys.pos.iter().enumerate() {
        for (j, q) in sys.pos.iter().enumerate().skip(i + 1) {
            let mut best2 = f64::INFINITY;
            for a in -2..=2 {
                for b in -2..=2 {
                    for c in -2..=2 {
                        let d2 = (0..3)
                            .map(|k| {
                                let d = p[k] - q[k]
                                    + f64::from(a) * cell[0][k]
                                    + f64::from(b) * cell[1][k]
                                    + f64::from(c) * cell[2][k];
                                d * d
                            })
                            .sum::<f64>();
                        best2 = best2.min(d2);
                    }
                }
            }
            if best2 < 1.5625 {
                expected.push((i, j));
            }
        }
    }
    assert!(expected.len() >= 4);
    list.pairs.sort_unstable();
    assert_eq!(list.pairs, expected);
}

#[test]
fn md_neighbour_grid_keeps_verlet_reference_on_empty_and_singleton_builds() {
    let mut sys = MdSystem::new(Vec::new(), SimBox::periodic([1.0e7; 3]));
    let mut list = NeighbourList::new(0.75, 0.25);
    assert!(list.needs_rebuild(&sys));
    list.build(&sys);
    assert!(!list.needs_rebuild(&sys));
    sys.pos.push([1.0; 3]);
    assert!(list.needs_rebuild(&sys));
    list.build(&sys);
    assert_eq!(list.builds, 2);
    assert!(list.pairs.is_empty());
    assert!(!list.needs_rebuild(&sys));
    sys.pos[0][0] += 0.0625;
    assert!(!list.needs_rebuild(&sys));
    sys.pos[0][0] += 0.125;
    assert!(list.needs_rebuild(&sys));
}
