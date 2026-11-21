//! Regression coverage for filtered simplicial events and binary reduction.

#![allow(clippy::unwrap_used)]

use std::collections::BTreeMap;

use rssn::kernels::topology::PersistenceError;
use rssn::kernels::topology::Simplex;
use rssn::kernels::topology::persistent_homology;

const INF: f64 = f64::INFINITY;

fn bars(
    filtration: &[(f64, Simplex)],
    max_dim: usize,
) -> Vec<Vec<(f64, f64)>> {
    persistent_homology(filtration, max_dim)
        .unwrap()
        .into_iter()
        .map(|diagram| {
            diagram
                .intervals
                .into_iter()
                .map(|bar| (bar.birth, bar.death))
                .collect()
        })
        .collect()
}

#[test]
fn repeated_canonical_simplices_enter_at_the_earliest_value() {
    let filtration = [
        (5.0, vec![2, 1, 0]),
        (2.0, vec![0, 2]),
        (1.0, vec![2, 1, 1]),
        (4.0, vec![1, 0]),
        (3.0, vec![1, 0, 2]),
        (1.0, vec![0, 1]),
        (4.0, vec![0]),
        (0.0, vec![2]),
        (0.0, vec![1]),
        (0.0, vec![0]),
        (1.0, vec![0, 1]),
        (2.0, vec![0]),
    ];
    assert_eq!(
        bars(&filtration, 2),
        [
            vec![(0.0, 1.0), (0.0, 1.0), (0.0, INF)],
            vec![(2.0, 3.0)],
            vec![]
        ]
    );
}

#[test]
fn signed_zeros_share_a_level_and_extended_real_values_are_ordered() {
    let filtration = [(0.0, vec![0]), (0.0, vec![1]), (-0.0, vec![0, 1])];
    assert_eq!(bars(&filtration, 1), [vec![(0.0, INF)], vec![]]);

    let extended = [(-INF, vec![0]), (0.0, vec![1]), (INF, vec![0, 1])];
    assert_eq!(bars(&extended, 0), [vec![(-INF, INF), (0.0, INF)]]);
}

#[test]
fn invalid_values_and_genuinely_late_faces_are_rejected() {
    for filtration in [
        vec![(0.0, vec![0]), (f64::NAN, vec![1])],
        vec![
            (2.0, vec![0]),
            (0.0, vec![1]),
            (1.0, vec![0, 1]),
            (3.0, vec![0]),
        ],
    ] {
        assert_eq!(
            persistent_homology(&filtration, 1),
            Err(PersistenceError::InvalidFiltration)
        );
    }
}

// A plain left-to-right F2 reduction, with a dense bit column and no clearing.
fn dense_bars(
    filtration: &[(f64, Simplex)],
    max_dim: usize,
) -> Vec<Vec<(f64, f64)>> {
    let mut earliest = BTreeMap::new();
    for (time, s) in filtration {
        let mut simplex = s.clone();
        simplex.sort_unstable();
        simplex.dedup();
        let value = earliest.entry(simplex).or_insert(*time);
        *value = value.min(*time);
    }
    let mut items: Vec<_> = earliest.into_iter().collect();
    items.sort_by(|(a, t), (b, u)| {
        t.partial_cmp(u)
            .unwrap()
            .then(a.len().cmp(&b.len()))
            .then(a.cmp(b))
    });
    assert!(items.len() <= u64::BITS as usize);

    let mut owners: Vec<Option<usize>> = vec![None; items.len()];
    let mut reduced = vec![0_u64; items.len()];
    let mut births = vec![false; items.len()];
    let mut out = vec![Vec::new(); max_dim + 1];
    for (column, (simplex, time)) in items.iter().enumerate() {
        let mut bits = 0_u64;
        for (row, (face, _)) in items[..column].iter().enumerate() {
            if face.len() + 1 == simplex.len() && face.iter().all(|v| simplex.contains(v)) {
                bits |= 1 << row;
            }
        }
        while bits != 0 {
            let low = (u64::BITS - 1 - bits.leading_zeros()) as usize;
            if let Some(owner) = owners[low] {
                bits ^= reduced[owner];
            } else {
                owners[low] = Some(column);
                let (birth, dim) = (items[low].1, items[low].0.len() - 1);
                if *time > birth && dim <= max_dim {
                    out[dim].push((birth, *time));
                }
                break;
            }
        }
        reduced[column] = bits;
        births[column] = bits == 0;
    }
    for (i, (s, birth)) in items.iter().enumerate() {
        let dim = s.len() - 1;
        if dim <= max_dim && births[i] && owners[i].is_none() {
            out[dim].push((*birth, INF));
        }
    }
    for intervals in &mut out {
        intervals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }
    out
}

#[test]
fn workspace_reduction_matches_dense_binary_columns() {
    let mut seed = 168707_u64;
    let mut next = || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        seed >> 32
    };
    for _ in 0..32 {
        let mut radii = [0_u32; 32];
        let mut filtration = Vec::new();
        for mask in 1_usize..32 {
            let vertices: Vec<_> = (0..5).filter(|&v| mask & (1 << v) != 0).collect();
            let oldest_face = vertices
                .iter()
                .map(|&v| radii[mask ^ (1 << v)])
                .max()
                .unwrap();
            let radius = oldest_face + u32::try_from(next() % 4).unwrap();
            radii[mask] = radius;
            if next() % 3 == 0 {
                let redundant = vertices
                    .iter()
                    .rev()
                    .copied()
                    .chain(vertices.first().copied())
                    .collect();
                filtration.push((f64::from(radius + 2), redundant));
            }
            filtration.push((f64::from(radius), vertices));
        }
        filtration.reverse();
        assert_eq!(bars(&filtration, 3), dense_bars(&filtration, 3));
    }
}
