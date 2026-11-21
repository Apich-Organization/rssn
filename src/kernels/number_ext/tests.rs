use super::*;

#[test]
fn bezout_identity_holds_at_the_signed_boundaries() {
    for a in [i64::MIN, i64::MIN + 1, -89, -1, 0, 1, 55, i64::MAX] {
        for b in [i64::MIN, -55, -2, 0, 2, 89, i64::MAX] {
            let (g, x, y) = xgcd(a, b);
            assert_eq!(
                g,
                i128::from(number_theory::gcd(a.unsigned_abs(), b.unsigned_abs()))
            );
            assert_eq!(i128::from(a) * x + i128::from(b) * y, g);
        }
    }
}

#[test]
fn modular_power_and_inverse_use_signed_canonical_residues() {
    for m in [1, 2, 97, 1 << 62, i64::MAX] {
        for base in [i64::MIN, -42, -1, 0, 1, 42, i64::MAX] {
            for exponent in [0, 1, 17, u64::MAX] {
                let result = mod_pow(base, exponent, m).unwrap();
                let expected = BigInt::from(base).modpow(&BigInt::from(exponent), &BigInt::from(m));
                assert_eq!(BigInt::from(result), expected);
                assert!((0..m).contains(&result));
            }
            match mod_inv(base, m) {
                | Ok(inverse) => {
                    assert!((0..m).contains(&inverse));
                    assert_eq!(
                        (i128::from(base) * i128::from(inverse)).rem_euclid(i128::from(m)),
                        i128::from(1 % m)
                    );
                },
                | Err(_) => assert_ne!(number_theory::gcd(base.unsigned_abs(), m as u64), 1),
            }
        }
    }
    for m in [i64::MIN, -7, 0] {
        assert!(mod_pow(2, 3, m).is_err());
        assert!(mod_inv(2, m).is_err());
    }
}

#[test]
fn crt_matches_exhaustive_residues_for_non_coprime_moduli() {
    for m in 1..=12 {
        for n in 1..=12 {
            let period = (1..=m * n).find(|x| x % m == 0 && x % n == 0).unwrap();
            for r in -m..m {
                for s in -n..n {
                    let expected = (0..period).find(|x| (x - r) % m == 0 && (x - s) % n == 0);
                    assert_eq!(crt(&[r, s], &[m, n]).ok(), expected.map(|x| (x, period)));
                }
            }
        }
    }
    assert_eq!(crt(&[], &[]), Err(NumberExtError::Empty));
    assert!(matches!(
        crt(&[0], &[]),
        Err(NumberExtError::ShapeMismatch(_))
    ));
    assert!(crt(&[0, 0], &[1, 0]).is_err());
    assert!(crt(&[0], &[i64::MIN]).is_err());
    assert_eq!(
        crt(&[i64::MIN, 0], &[i64::MAX, 1]),
        Ok((i64::MAX - 1, i64::MAX))
    );
    assert_eq!(
        crt(&[-1, -1, -1], &[i64::MAX; 3]),
        Ok((i64::MAX - 1, i64::MAX))
    );
    assert!(crt(&[0, 0], &[i64::MAX, i64::MAX - 1]).is_err());
}

#[test]
fn linear_congruence_gives_every_solution_and_only_solutions() {
    for modulus in 1..=16_u64 {
        for a in -17..=17_i64 {
            for b in -17..=17_i64 {
                let solved = solve_linear_congruence(a, b, modulus).unwrap();
                let matches: Vec<_> = (0..modulus)
                    .filter(|&x| (a * x as i64 - b).rem_euclid(modulus as i64) == 0)
                    .collect();
                if let Some((first, period)) = solved {
                    assert!(first < period);
                    assert_eq!(modulus % period, 0);
                    assert_eq!(
                        (first..modulus)
                            .step_by(period as usize)
                            .collect::<Vec<_>>(),
                        matches
                    );
                } else {
                    assert_eq!(matches, [] as [u64; 0]);
                }
            }
        }
    }
    assert!(solve_linear_congruence(1, 1, 0).is_err());
    assert_eq!(solve_linear_congruence(0, 0, u64::MAX), Ok(Some((0, 1))));
    assert_eq!(solve_linear_congruence(0, i64::MIN, u64::MAX), Ok(None));
    assert_eq!(
        solve_linear_congruence(-3, -6, u64::MAX),
        Ok(Some((2, u64::MAX / 3)))
    );
    for a in [i64::MIN, -1, 1, i64::MAX] {
        for b in [i64::MIN, -1, 0, i64::MAX] {
            let modulus = u64::MAX - 58;
            let (first, period) = solve_linear_congruence(a, b, modulus).unwrap().unwrap();
            assert_eq!(period, modulus);
            assert!(first < period);
            assert_eq!(
                (i128::from(a) * i128::from(first) - i128::from(b)).rem_euclid(i128::from(modulus)),
                0
            );
        }
    }
}

#[test]
fn mobius_sieve_agrees_with_divisor_sums_and_scalar_factorization() {
    assert_eq!(mobius_sieve(0).unwrap(), [0]);
    assert_eq!(mobius_sieve(1).unwrap(), [0, 1]);
    assert!(mobius_sieve(usize::MAX).is_err());
    assert!(mobius_sieve(usize::MAX / 2).is_err());
    let values = mobius_sieve(4_096).unwrap();
    let mut sums = vec![0_i64; values.len()];
    for (n, mu) in values.into_iter().enumerate().skip(1) {
        assert_eq!(mu, mobius(n as i64), "{n}");
        for index in (n..sums.len()).step_by(n) {
            sums[index] += mu;
        }
    }
    assert_eq!(&sums[..2], [0, 1]);
    assert!(sums[2..].iter().all(|&x| x == 0));
    assert_eq!(mobius(i64::MIN), 0);
    assert_eq!(mobius(i64::MAX), 0); // Divisible by 7 squared.
    let prime = 2_147_483_647_i64;
    assert_eq!(mobius(prime), -1);
    assert_eq!(mobius(prime * 2_147_483_629), 1);
    assert_eq!(mobius(prime * prime), 0);
}

#[test]
fn partition_count_matches_an_independent_coin_change_recurrence() {
    let mut counts = vec![BigInt::default(); 101];
    counts[0] = BigInt::from(1u8);
    for coin in 1..counts.len() {
        for total in coin..counts.len() {
            let addend = counts[total - coin].clone();
            counts[total] += addend;
        }
    }
    for (n, count) in counts.into_iter().enumerate() {
        assert_eq!(partition_count(n as i64), count);
    }
    assert_eq!(partition_count(i64::MIN), BigInt::default());
    assert_eq!(
        partition_count(1_000).to_string(),
        "24061467864032622473692149727991"
    );
}

#[test]
fn stirling_kernels_match_polynomial_coefficients_and_surjections() {
    let mut coefficients = vec![BigInt::from(1u8)];
    for n in 0..=18 {
        for (k, coefficient) in coefficients.iter().enumerate() {
            assert_eq!(stirling_first_signed(n, k as i64), *coefficient);
            assert_eq!(
                stirling_first(n, k as i64),
                BigInt::from(coefficient.magnitude().clone())
            );

            let mut over = BigInt::from(1u8);
            let mut sum = BigInt::default();
            for j in 0..=k {
                let term = &over * BigInt::from(j).pow(n as u32);
                if (k - j).is_multiple_of(2) {
                    sum += term;
                } else {
                    sum -= term;
                }
                over = over * (k - j) / (j + 1);
            }
            let factorial = (1..=k as u64).product::<u64>();
            assert_eq!(stirling_second_exact(n, k as i64), sum / factorial);
        }
        let mut next = vec![BigInt::default(); coefficients.len() + 1];
        for (index, coefficient) in coefficients.into_iter().enumerate() {
            next[index] -= &coefficient * n;
            next[index + 1] += coefficient;
        }
        coefficients = next;
    }
    for stirling in [stirling_first, stirling_first_signed, stirling_second_exact] {
        assert_eq!(stirling(i64::MIN, i64::MAX), BigInt::default());
        assert_eq!(stirling(5, -1), BigInt::default());
        assert_eq!(stirling(5, 6), BigInt::default());
        assert_eq!(stirling(i64::MAX, i64::MAX), BigInt::from(1u8));
    }
    let n = i64::MAX;
    let pair = BigInt::from(n) * (n - 1) / 2;
    assert_eq!(stirling_first(n, n - 1), pair);
    assert_eq!(stirling_second_exact(n, n - 1), pair);
    assert_eq!(stirling_second_exact(n, 1), BigInt::from(1u8));
}

#[test]
fn catalan_exact_matches_the_convolution_recurrence() {
    let mut counts = vec![BigInt::from(1u8)];
    for n in 1..=25 {
        let count = (0..n)
            .map(|k| &counts[k] * &counts[n - k - 1])
            .sum::<BigInt>();
        assert_eq!(catalan_exact(n as i64), count);
        counts.push(count);
    }
    assert_eq!(catalan_exact(0), counts[0]);
    assert_eq!(catalan_exact(i64::MIN), BigInt::default());
}
