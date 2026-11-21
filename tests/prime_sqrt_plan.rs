//! Independent enumeration and exact word-scale tests for modular square roots.

use rssn::kernels::number_theory::PrimeSqrtPlan;
use rssn::kernels::number_theory::sqrt_mod_prime;

fn trial_prime(n: u64) -> bool {
    n >= 2 && (2..).take_while(|&d| d <= n / d).all(|d| n % d != 0)
}

fn reference_square(
    n: u64,
    p: u64,
) -> u64 {
    let n = u128::from(n);
    (n * n % u128::from(p)) as u64
}

// Evaluate from the leading binary digit. This oracle shares no routine
// with the candidate or with the primality test it uses internally.
fn reference_power(
    a: u64,
    exponent: u64,
    modulus: u64,
) -> u64 {
    let modulus = u128::from(modulus);
    let mut residue = 1_u128;
    for bit in (0..u64::BITS).rev() {
        residue = residue * residue % modulus;
        if (exponent >> bit) & 1 == 1 {
            residue = residue * u128::from(a) % modulus;
        }
    }
    residue as u64
}

#[test]
fn exhaustive_small_moduli_match_the_least_root_by_enumeration() {
    for p in 0_u64..=263 {
        if !trial_prime(p) {
            assert!(PrimeSqrtPlan::new(p).is_none());
            assert_eq!(sqrt_mod_prime(0, p), None);
            assert_eq!(sqrt_mod_prime(1, p), None);
            continue;
        }

        let plan = PrimeSqrtPlan::new(p).unwrap();
        assert_eq!(plan.modulus(), p);
        let mut reference = vec![None; p as usize];
        for root in 0..p {
            reference[reference_square(root, p) as usize].get_or_insert(root);
        }
        for n in 0..2 * p {
            let expected = reference[(n % p) as usize];
            assert_eq!(plan.sqrt(n), expected, "modulus {p}, value {n}");
            assert_eq!(sqrt_mod_prime(n, p), expected, "modulus {p}, value {n}");
        }
    }
}

#[test]
fn largest_word_scale_primes_recover_canonical_planted_roots() {
    const MERSENNE: u64 = (1 << 61) - 1;
    const LARGEST: u64 = u64::MAX - 58;
    const GOLDILOCKS: u64 = u64::MAX - (1 << 32) + 2;

    // Factors of p-1 certify that 7 has multiplicative order p-1; this
    // proves primality here and supplies an independent non-residue.
    let factors = [2, 3, 5, 17, 257, 65_537];
    assert!(factors.into_iter().all(trial_prime));
    assert_eq!(
        factors[1..].iter().product::<u64>() * (1 << 32),
        GOLDILOCKS - 1
    );
    assert_eq!(reference_power(7, GOLDILOCKS - 1, GOLDILOCKS), 1);
    for q in factors {
        assert_ne!(reference_power(7, (GOLDILOCKS - 1) / q, GOLDILOCKS), 1);
    }

    for (p, non_residue) in [(MERSENNE, MERSENNE - 1), (LARGEST, 2), (GOLDILOCKS, 7)] {
        let plan = PrimeSqrtPlan::new(p).unwrap();
        assert_eq!(reference_power(non_residue, (p - 1) / 2, p), p - 1);
        assert_eq!(plan.sqrt(non_residue), None);
        assert_eq!(sqrt_mod_prime(non_residue, p), None);
        assert_eq!(plan.sqrt(u64::MAX), plan.sqrt(u64::MAX % p));
        for x in [
            0,
            1,
            2,
            p / 4,
            p / 2,
            p / 2 + 1,
            p - 1,
            p,
            u64::MAX - 1,
            u64::MAX,
        ] {
            let reduced = x % p;
            let expected = Some(reduced.min(p - reduced));
            let n = reference_square(x, p);
            assert_eq!(plan.sqrt(n), expected, "modulus {p}, root {x}");
            assert_eq!(sqrt_mod_prime(n, p), expected, "modulus {p}, root {x}");
        }
        let mut bits = p;
        for _ in 0..256 {
            bits = bits.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let x = bits % p;
            let square = reference_square(x, p);
            assert_eq!(plan.sqrt(square), Some(x.min(p - x)));
            if square != 0 {
                let non_square =
                    (u128::from(square) * u128::from(non_residue) % u128::from(p)) as u64;
                assert_eq!(plan.sqrt(non_square), None);
            }
        }
    }
}

#[test]
fn degenerate_residues_and_word_boundary_composites() {
    let binary = PrimeSqrtPlan::new(2).unwrap();
    for n in [0, 1, 2, u64::MAX - 1, u64::MAX] {
        assert_eq!(binary.sqrt(n), Some(n & 1));
        assert_eq!(sqrt_mod_prime(n, 2), Some(n & 1));
    }
    for n in [
        0,
        1,
        561,
        3_215_031_751,
        3_825_123_056_546_413_051,
        u64::MAX - 2,
        u64::MAX,
    ] {
        assert!(PrimeSqrtPlan::new(n).is_none());
        assert_eq!(sqrt_mod_prime(u64::MAX, n), None);
    }
}
