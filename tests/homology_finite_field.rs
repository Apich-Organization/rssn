//! Independent finite-field rank and Betti number checks.

use rssn::kernels::homology::ChainComplex;
use rssn::kernels::homology::complex_chains;
use rssn::kernels::homology::is_prime;
use rssn::kernels::homology::projective_plane;
use rssn::kernels::homology::rank_mod_p;
use std::collections::HashSet;

const WIDE_PRIME: u64 = 18_446_744_073_709_551_557;

#[test]
fn full_word_prime_and_signed_coefficients_never_overflow() {
    assert!(is_prime(WIDE_PRIME));
    assert!(!is_prime(u64::MAX));
    assert!(!is_prime(3_215_031_751));
    assert!(!is_prime(3_825_123_056_546_413_051));
    assert_eq!(rank_mod_p(&[vec![-1, -1], vec![-1, 1]], WIDE_PRIME), 2);
    assert_eq!(
        rank_mod_p(&[vec![i64::MIN, -1], vec![i64::MIN, -1]], WIDE_PRIME),
        1
    );
    assert_eq!(rank_mod_p(&[vec![0; 4], vec![0; 4]], WIDE_PRIME), 0);
}

#[test]
fn default_binary_rows_cross_each_machine_word_boundary() {
    for cols in [63, 64, 65, 127, 128, 129] {
        let mut rows = vec![vec![i64::MIN; cols]; 3];
        rows[0][cols - 1] = -3;
        rows[1][0] = i64::MAX;
        rows[2][0] = -1;
        rows[2][cols - 1] = 9;
        assert_eq!(rank_mod_p(&rows, 2), 2);
        rows.push(vec![0; cols]);
        assert_eq!(rank_mod_p(&rows, 2), 2);
    }
    assert_eq!(rank_mod_p(&[vec![], vec![]], 2), 0);
}

fn image_size_rank(
    m: &[Vec<i64>],
    prime: i64,
) -> usize {
    let cols = m[0].len();
    let combinations = prime.pow(u32::try_from(cols).unwrap());
    let mut images = HashSet::new();
    for mut counter in 0..combinations {
        let digits = (0..cols)
            .map(|_| {
                let digit = counter % prime;
                counter /= prime;
                digit
            })
            .collect::<Vec<_>>();
        let image = m
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&digits)
                    .map(|(x, y)| x * y)
                    .sum::<i64>()
                    .rem_euclid(prime)
            })
            .collect::<Vec<_>>();
        images.insert(image);
    }
    let mut size = images.len() as i64;
    let mut dimension = 0;
    while size > 1 {
        assert_eq!(size % prime, 0);
        size /= prime;
        dimension += 1;
    }
    dimension
}

#[test]
fn row_reduction_matches_exhaustive_linear_images() {
    for prime in [2_i64, 3, 5] {
        for rows in 1..=4 {
            for cols in 1..=4 {
                for seed in 0..24 {
                    let mut matrix = vec![vec![0; cols]; rows];
                    let mut state = 19_991_u64 + seed;
                    for x in matrix.iter_mut().flatten() {
                        state = state
                            .wrapping_mul(2_862_933_555_777_941_757)
                            .wrapping_add(3_037_000_493);
                        *x = ((state >> 32) % 17) as i64 - 8;
                    }
                    assert_eq!(
                        rank_mod_p(&matrix, prime as u64),
                        image_size_rank(&matrix, prime)
                    );
                }
            }
        }
    }
}

#[test]
fn batched_betti_detects_projective_plane_two_torsion() {
    let projective = complex_chains(&projective_plane(0));
    assert_eq!(projective.betti_mod_p(2), [1, 1, 1]);
    for prime in [3, WIDE_PRIME] {
        assert_eq!(projective.betti_mod_p(prime), [1, 0, 0]);
        for (degree, dimension) in projective.betti_mod_p(prime).into_iter().enumerate() {
            assert_eq!(projective.homology_mod_p(degree, prime), dimension);
        }
    }
    assert_eq!(projective.homology_mod_p(usize::MAX, 2), 0);
    assert_eq!(projective.homology(usize::MAX).rank, 0);
    let empty = ChainComplex {
        dims: Vec::new(),
        d: Vec::new(),
    };
    assert!(empty.betti_mod_p(2).is_empty());
}
