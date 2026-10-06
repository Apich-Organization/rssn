//! More counting sequences and combinatorial structures: Lah, Eulerian,
//! Narayana, Motzkin and Schröder numbers, restricted partitions and
//! compositions, necklaces, Bernoulli and Euler numbers, and permutation
//! statistics (rank, unrank, inversions, cycle type) and Young tableaux.
//!
//! | operator | value |
//! |---|---|
//! | `lah(n, k)` | Lah number `C(n-1, k-1) n!/k!`, the ways to split `n` items into `k` ordered lists |
//! | `eulerian(n, k)` | Eulerian number: permutations of `n` with `k` ascents |
//! | `narayana(n, k)` | `C(n, k) C(n, k-1) / n` |
//! | `motzkin(n)`, `schroeder(n)` | Motzkin and large Schröder numbers |
//! | `fubini(n)` | the ordered Bell number `Σ k! S(n, k)` |
//! | `central_binomial(n)` | `C(2n, n)` |
//! | `ballot(a, b)` | lattice paths with `a` up and `b` down steps staying non-negative |
//! | `involutions(n)` | permutations with `p² = id` |
//! | `rencontres(n, k)` | permutations of `n` with exactly `k` fixed points |
//! | `partitions_k(n, k)` | partitions of `n` into exactly `k` parts |
//! | `partitions_max(n, k)` | partitions of `n` with every part at most `k` |
//! | `partitions_distinct(n)` | partitions of `n` into distinct parts |
//! | `partition_list(n)` | the partitions themselves, in reverse lexicographic order (`n <= 40`) |
//! | `compositions(n, k)`, `weak_compositions(n, k)` | ordered sums of `k` positive (non-negative) integers equal to `n` |
//! | `necklaces(n, k)`, `bracelets(n, k)` | `k`-colourings of `n` beads up to rotation (and reflection), by Burnside's lemma |
//! | `bernoulli_number(n)` | the Bernoulli number `B_n` (`B_1 = -1/2`) |
//! | `euler_number(n)`, `zigzag(n)` | the secant numbers `E_n` (signed, zero for odd `n`) and the up-down numbers |
//! | `permutation_rank(p)`, `permutation_unrank(n, r)` | lexicographic rank (from 0) of a permutation of `1..n` and its inverse |
//! | `inversions(p)`, `permutation_sign(p)`, `cycle_type(p)`, `permutation_order(p)` | inversion count, sign, sorted cycle lengths and order |
//! | `hook_lengths(lambda)`, `syt_count(lambda)` | hook lengths of a Young diagram and the number of standard tableaux (hook length formula) |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::choose;
use super::factorial;
use super::nat;
use super::stirling_row;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::rules::number_theory::exact;
use crate::rules::number_theory::exact_lists;
use crate::rules::number_theory::Value;

const SET: &str = "combinatorics";
const MAX_N: u64 = 2000;

fn big(v: u64) -> BigInt {
    BigInt::from(v)
}

fn binom(
    n: u64,
    k: u64,
) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    choose(&big(n), k).unwrap_or_default()
}

fn eulerian(
    n: u64,
    k: u64,
) -> BigInt {
    let mut row = vec![BigInt::one()];
    for m in 2..=n {
        let mut next = vec![BigInt::zero(); row.len() + 1];
        for (j, slot) in next.iter_mut().enumerate() {
            let j64 = u64::try_from(j).unwrap_or(0);
            let above = row.get(j).cloned().unwrap_or_default();
            let left = if j == 0 { BigInt::zero() } else { row.get(j - 1).cloned().unwrap_or_default() };
            *slot = above * (j64 + 1) + left * (m - j64);
        }
        row = next;
    }
    usize::try_from(k).ok().and_then(|k| row.get(k).cloned()).unwrap_or_default()
}

fn motzkin(n: u64) -> BigInt {
    let (mut previous, mut current) = (BigInt::one(), BigInt::one());
    if n == 0 {
        return previous;
    }
    for m in 2..=n {
        let next = (&current * (2 * m + 1) + &previous * (3 * m - 3)) / (m + 2);
        previous = current;
        current = next;
    }
    current
}

fn schroeder(n: u64) -> BigInt {
    let (mut previous, mut current) = (BigInt::one(), BigInt::from(2));
    if n == 0 {
        return previous;
    }
    for m in 2..=n {
        let next = (&current * (6 * m - 3) - &previous * (m - 2)) / (m + 1);
        previous = current;
        current = next;
    }
    current
}

fn fubini(n: u64) -> BigInt {
    let row = stirling_row(n, |_, j| j);
    row.iter().enumerate().map(|(k, s)| s * factorial(u64::try_from(k).unwrap_or(0))).sum()
}

/// Partitions of `n` into exactly `k` parts, by `p(n, k) = p(n-1, k-1) + p(n-k, k)`.
fn partitions_exact(
    n: usize,
    k: usize,
) -> BigInt {
    let mut table = vec![vec![BigInt::zero(); k + 1]; n + 1];
    table[0][0] = BigInt::one();
    for m in 1..=n {
        for j in 1..=k.min(m) {
            let a = table[m - 1][j - 1].clone();
            let b = if m >= j { table[m - j][j].clone() } else { BigInt::zero() };
            table[m][j] = a + b;
        }
    }
    table[n][k].clone()
}

fn partitions_distinct(n: usize) -> BigInt {
    let mut ways = vec![BigInt::zero(); n + 1];
    ways[0] = BigInt::one();
    for part in 1..=n {
        for total in (part..=n).rev() {
            let add = ways[total - part].clone();
            ways[total] += add;
        }
    }
    ways[n].clone()
}

fn partition_list(n: u64) -> Vec<Vec<BigInt>> {
    fn go(
        remaining: u64,
        bound: u64,
        prefix: &mut Vec<u64>,
        out: &mut Vec<Vec<BigInt>>,
    ) {
        if remaining == 0 {
            out.push(prefix.iter().map(|&p| BigInt::from(p)).collect());
            return;
        }
        for part in (1..=remaining.min(bound)).rev() {
            prefix.push(part);
            go(remaining - part, part, prefix, out);
            prefix.pop();
        }
    }
    let mut out = Vec::new();
    go(n, n, &mut Vec::new(), &mut out);
    out
}

fn necklaces(
    n: u64,
    k: &BigInt,
) -> Option<BigInt> {
    if n == 0 {
        return Some(BigInt::one());
    }
    let mut total = BigInt::zero();
    for d in 1..=n {
        if n % d == 0 {
            let phi = (1..=d).filter(|&x| gcd_u64(x, d) == 1).count();
            total += BigInt::from(phi) * num_traits::pow::pow(k.clone(), usize::try_from(n / d).ok()?);
        }
    }
    Some(total / n)
}

const fn gcd_u64(
    mut a: u64,
    mut b: u64,
) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn bracelets(
    n: u64,
    k: &BigInt,
) -> Option<BigInt> {
    if n == 0 {
        return Some(BigInt::one());
    }
    let rotations = necklaces(n, k)?;
    let exponent = |e: u64| num_traits::pow::pow(k.clone(), usize::try_from(e).ok().unwrap_or(0));
    let reflections = if n % 2 == 1 {
        BigInt::from(n) * exponent(n.div_ceil(2))
    } else {
        BigInt::from(n / 2) * (exponent(n / 2 + 1) + exponent(n / 2))
    };
    // (rotations * n + reflections) / 2n
    Some((rotations * n + reflections) / (2 * n))
}

/// Bernoulli numbers by the Akiyama–Tanigawa transform (`B_1 = -1/2`).
fn bernoulli(n: usize) -> BigRational {
    let mut a: Vec<BigRational> = Vec::with_capacity(n + 1);
    for m in 0..=n {
        a.push(BigRational::new(BigInt::one(), BigInt::from(m + 1)));
        for j in (1..=m).rev() {
            let value = (a[j - 1].clone() - a[j].clone()) * BigRational::from_integer(BigInt::from(j));
            a[j - 1] = value;
        }
    }
    let value = a.first().cloned().unwrap_or_default();
    if n == 1 { -value } else { value }
}

/// Up-down (Euler zigzag) numbers by the Seidel–Entringer triangle.
fn zigzag(n: usize) -> BigInt {
    let mut row = vec![BigInt::one()];
    for m in 1..=n {
        let mut next = vec![BigInt::zero()];
        for k in 1..=m {
            let value = next[k - 1].clone() + row[m - k].clone();
            next.push(value);
        }
        row = next;
    }
    row.last().cloned().unwrap_or_default()
}

/// A permutation of `1..n` read from a list of integers.
fn permutation(items: &[BigInt]) -> Option<Vec<usize>> {
    let n = items.len();
    let mut seen = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for v in items {
        let v = usize::try_from(v.to_u64()?).ok()?;
        if v == 0 || v > n || seen[v - 1] {
            return None;
        }
        seen[v - 1] = true;
        out.push(v - 1);
    }
    Some(out)
}

fn cycle_lengths(perm: &[usize]) -> Vec<u64> {
    let mut seen = vec![false; perm.len()];
    let mut lengths = Vec::new();
    for start in 0..perm.len() {
        if seen[start] {
            continue;
        }
        let (mut at, mut length) = (start, 0_u64);
        while !seen[at] {
            seen[at] = true;
            at = perm[at];
            length += 1;
        }
        lengths.push(length);
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    lengths
}

fn inversions(perm: &[usize]) -> u64 {
    let mut total = 0_u64;
    for (i, &a) in perm.iter().enumerate() {
        total += u64::try_from(perm[i + 1..].iter().filter(|&&b| b < a).count()).unwrap_or(0);
    }
    total
}

fn partition_shape(items: &[BigInt]) -> Option<Vec<usize>> {
    let shape: Vec<usize> = items.iter().map(|v| v.to_u64().and_then(|v| usize::try_from(v).ok())).collect::<Option<_>>()?;
    let valid = shape.windows(2).all(|w| w[0] >= w[1]) && shape.iter().all(|&p| p > 0);
    (valid && shape.len() <= 60 && shape.iter().sum::<usize>() <= 400).then_some(shape)
}

fn hooks(shape: &[usize]) -> Vec<Vec<u64>> {
    shape
        .iter()
        .enumerate()
        .map(|(row, &length)| {
            (0..length)
                .map(|col| {
                    let arm = length - col - 1;
                    let leg = shape[row + 1..].iter().filter(|&&below| below > col).count();
                    u64::try_from(arm + leg + 1).unwrap_or(0)
                })
                .collect()
        })
        .collect()
}

#[allow(clippy::too_many_lines)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));

    exact(i, SET, binary("lah"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if n == 0 || k == 0 || k > n {
            return Some(Value::Int(if n == 0 && k == 0 { BigInt::one() } else { BigInt::zero() }));
        }
        Some(Value::Int(binom(n - 1, k - 1) * factorial(n) / factorial(k)))
    })?;
    exact(i, SET, binary("eulerian"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 500)?, nat(k, 500)?);
        if n == 0 {
            return Some(Value::Int(if k == 0 { BigInt::one() } else { BigInt::zero() }));
        }
        Some(Value::Int(eulerian(n, k)))
    })?;
    exact(i, SET, binary("narayana"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if n == 0 || k == 0 || k > n {
            return Some(Value::Int(BigInt::zero()));
        }
        Some(Value::Int(binom(n, k) * binom(n, k - 1) / n))
    })?;
    exact(i, SET, unary("motzkin"), |a| match a {
        | [n] => Some(Value::Int(motzkin(nat(n, MAX_N)?))),
        | _ => None,
    })?;
    exact(i, SET, unary("schroeder"), |a| match a {
        | [n] => Some(Value::Int(schroeder(nat(n, MAX_N)?))),
        | _ => None,
    })?;
    exact(i, SET, unary("fubini"), |a| match a {
        | [n] => Some(Value::Int(fubini(nat(n, 400)?))),
        | _ => None,
    })?;
    exact(i, SET, unary("central_binomial"), |a| match a {
        | [n] => Some(Value::Int(binom(2 * nat(n, MAX_N)?, nat(n, MAX_N)?))),
        | _ => None,
    })?;
    exact(i, SET, binary("ballot"), |a| {
        let [up, down] = a else { return None };
        let (up, down) = (nat(up, MAX_N)?, nat(down, MAX_N)?);
        if down > up {
            return Some(Value::Int(BigInt::zero()));
        }
        Some(Value::Int(binom(up + down, down) * (up - down + 1) / (up + 1)))
    })?;
    exact(i, SET, unary("involutions"), |a| {
        let [n] = a else { return None };
        let n = nat(n, MAX_N)?;
        let (mut previous, mut current) = (BigInt::one(), BigInt::one());
        for m in 2..=n {
            let next = &current + &previous * (m - 1);
            previous = current;
            current = next;
        }
        Some(Value::Int(if n == 0 { previous } else { current }))
    })?;
    exact(i, SET, binary("rencontres"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        if k > n {
            return Some(Value::Int(BigInt::zero()));
        }
        let m = n - k;
        let (mut previous, mut current) = (BigInt::one(), BigInt::zero());
        let derangements = match m {
            | 0 => previous.clone(),
            | _ => {
                for j in 2..=m {
                    let next = (&previous + &current) * (j - 1);
                    previous = current;
                    current = next;
                }
                current
            },
        };
        Some(Value::Int(binom(n, k) * derangements))
    })?;
    exact(i, SET, binary("partitions_k"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 600)?, nat(k, 600)?);
        Some(Value::Int(partitions_exact(usize::try_from(n).ok()?, usize::try_from(k).ok()?)))
    })?;
    exact(i, SET, binary("partitions_max"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, 600)?, nat(k, 600)?);
        // Conjugation: parts at most k <-> at most k parts.
        let (n, k) = (usize::try_from(n).ok()?, usize::try_from(k).ok()?);
        let total: BigInt = (0..=k).map(|j| partitions_exact(n, j)).sum();
        Some(Value::Int(total))
    })?;
    exact(i, SET, unary("partitions_distinct"), |a| match a {
        | [n] => Some(Value::Int(partitions_distinct(usize::try_from(nat(n, 2000)?).ok()?))),
        | _ => None,
    })?;
    exact(i, SET, unary("partition_list"), |a| match a {
        | [n] => Some(Value::Nested(partition_list(nat(n, 40)?))),
        | _ => None,
    })?;
    exact(i, SET, binary("compositions"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        Some(Value::Int(if n == 0 || k == 0 {
            if n == 0 && k == 0 { BigInt::one() } else { BigInt::zero() }
        } else {
            binom(n - 1, k - 1)
        }))
    })?;
    exact(i, SET, binary("weak_compositions"), |a| {
        let [n, k] = a else { return None };
        let (n, k) = (nat(n, MAX_N)?, nat(k, MAX_N)?);
        Some(Value::Int(if k == 0 {
            if n == 0 { BigInt::one() } else { BigInt::zero() }
        } else {
            binom(n + k - 1, k - 1)
        }))
    })?;
    exact(i, SET, binary("necklaces"), |a| {
        let [n, k] = a else { return None };
        necklaces(nat(n, 400)?, k).map(Value::Int)
    })?;
    exact(i, SET, binary("bracelets"), |a| {
        let [n, k] = a else { return None };
        bracelets(nat(n, 400)?, k).map(Value::Int)
    })?;
    exact(i, SET, unary("bernoulli_number"), |a| {
        let [n] = a else { return None };
        Some(Value::Rat(bernoulli(usize::try_from(nat(n, 300)?).ok()?)))
    })?;
    exact(i, SET, unary("zigzag"), |a| {
        let [n] = a else { return None };
        Some(Value::Int(zigzag(usize::try_from(nat(n, 1000)?).ok()?)))
    })?;
    exact(i, SET, unary("euler_number"), |a| {
        let [n] = a else { return None };
        let n = nat(n, 1000)?;
        if n % 2 == 1 {
            return Some(Value::Int(BigInt::zero()));
        }
        let z = zigzag(usize::try_from(n).ok()?);
        Some(Value::Int(if n % 4 == 2 { -z } else { z }))
    })?;
    exact_lists(i, SET, unary("permutation_rank"), |a| {
        let [items] = a else { return None };
        let perm = permutation(items)?;
        let n = perm.len();
        let mut rank = BigInt::zero();
        for (idx, &value) in perm.iter().enumerate() {
            let smaller = perm[idx + 1..].iter().filter(|&&later| later < value).count();
            rank += BigInt::from(smaller) * factorial(u64::try_from(n - idx - 1).ok()?);
        }
        Some(Value::Int(rank))
    })?;
    exact(i, SET, binary("permutation_unrank"), |a| {
        let [n, rank] = a else { return None };
        let n = nat(n, 500)?;
        if rank.is_negative() || *rank >= factorial(n) {
            return None;
        }
        let mut pool: Vec<u64> = (1..=n).collect();
        let mut rank = rank.clone();
        let mut out = Vec::new();
        for place in (0..n).rev() {
            let base = factorial(place);
            let index = (&rank / &base).to_usize()?;
            rank %= base;
            out.push(big(pool.remove(index)));
        }
        Some(Value::Flat(out))
    })?;
    exact_lists(i, SET, unary("inversions"), |a| {
        let [items] = a else { return None };
        Some(Value::Int(big(inversions(&permutation(items)?))))
    })?;
    exact_lists(i, SET, unary("permutation_sign"), |a| {
        let [items] = a else { return None };
        let even = inversions(&permutation(items)?) % 2 == 0;
        Some(Value::Int(BigInt::from(if even { 1 } else { -1 })))
    })?;
    exact_lists(i, SET, unary("cycle_type"), |a| {
        let [items] = a else { return None };
        Some(Value::Flat(cycle_lengths(&permutation(items)?).into_iter().map(BigInt::from).collect()))
    })?;
    exact_lists(i, SET, unary("permutation_order"), |a| {
        let [items] = a else { return None };
        let lcm = cycle_lengths(&permutation(items)?)
            .into_iter()
            .fold(BigInt::one(), |acc, l| &acc / num_integer::gcd(acc.clone(), BigInt::from(l)) * l);
        Some(Value::Int(lcm))
    })?;
    exact_lists(i, SET, unary("hook_lengths"), |a| {
        let [items] = a else { return None };
        let shape = partition_shape(items)?;
        Some(Value::Nested(hooks(&shape).into_iter().map(|row| row.into_iter().map(BigInt::from).collect()).collect()))
    })?;
    exact_lists(i, SET, unary("syt_count"), |a| {
        let [items] = a else { return None };
        let shape = partition_shape(items)?;
        let cells: usize = shape.iter().sum();
        let denominator: BigInt = hooks(&shape).into_iter().flatten().map(BigInt::from).product();
        Some(Value::Int(factorial(u64::try_from(cells).ok()?) / denominator))
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::rules::combinatorics::combinatorics;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[combinatorics()], src)
    }

    #[test]
    fn classical_triangles_and_sequences() {
        assert_eq!(s("lah(4, 2)"), "36");
        assert_eq!(s("lah(5, 1)"), "120");
        assert_eq!(s("eulerian(4, 1)"), "11");
        assert_eq!(s("eulerian(5, 2)"), "66");
        assert_eq!(s("narayana(4, 2)"), "6");
        assert_eq!(s("motzkin(6)"), "51");
        assert_eq!(s("schroeder(4)"), "90");
        assert_eq!(s("fubini(4)"), "75");
        assert_eq!(s("central_binomial(5)"), "252");
        assert_eq!(s("ballot(4, 2)"), "9");
        assert_eq!(s("involutions(5)"), "26");
        assert_eq!(s("rencontres(5, 2)"), "20");
    }

    #[test]
    fn restricted_partitions_and_compositions() {
        assert_eq!(s("partitions_k(7, 3)"), "4");
        assert_eq!(s("partitions_max(7, 3)"), "8");
        assert_eq!(s("partitions_distinct(10)"), "10");
        assert_eq!(s("partition_list(4)"), "list(list(4), list(3, 1), list(2, 2), list(2, 1, 1), list(1, 1, 1, 1))");
        assert_eq!(s("compositions(5, 3)"), "6");
        assert_eq!(s("weak_compositions(5, 3)"), "21");
        assert_eq!(s("necklaces(6, 2)"), "14");
        assert_eq!(s("bracelets(6, 2)"), "13");
    }

    #[test]
    fn bernoulli_and_euler_numbers() {
        assert_eq!(s("bernoulli_number(0)"), "1");
        assert_eq!(s("bernoulli_number(1)"), "-1/2");
        assert_eq!(s("bernoulli_number(2)"), "1/6");
        assert_eq!(s("bernoulli_number(12)"), "-691/2730");
        assert_eq!(s("bernoulli_number(13)"), "0");
        assert_eq!(s("zigzag(6)"), "61");
        assert_eq!(s("euler_number(4)"), "5");
        assert_eq!(s("euler_number(6)"), "-61");
    }

    #[test]
    fn permutation_statistics_and_tableaux() {
        assert_eq!(s("permutation_rank(list(2, 3, 1))"), "3");
        assert_eq!(s("permutation_unrank(3, 3)"), "list(2, 3, 1)");
        assert_eq!(s("inversions(list(3, 1, 2))"), "2");
        assert_eq!(s("permutation_sign(list(2, 1, 3))"), "-1");
        assert_eq!(s("cycle_type(list(2, 1, 4, 5, 3))"), "list(3, 2)");
        assert_eq!(s("permutation_order(list(2, 1, 4, 5, 3))"), "6");
        assert_eq!(s("hook_lengths(list(3, 2))"), "list(list(4, 3, 1), list(2, 1))");
        assert_eq!(s("syt_count(list(3, 2))"), "5");
        assert_eq!(s("syt_count(list(3, 2, 1))"), "16");
    }
}
