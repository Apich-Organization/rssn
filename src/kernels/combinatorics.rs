//! # Numerical Combinatorics
//!
//! This module provides numerical implementations for combinatorial functions.
//! It includes functions for calculating factorials, permutations, and combinations,
//! as well as set-partition counts and a numerical solver for linear recurrence relations.

/// Computes the factorial of `n` (`n!`) as an `f64`.
///
/// # Arguments
/// * `n` - The non-negative integer for which to compute the factorial.
///
/// # Returns
/// The factorial of `n` as an `f64`. Returns `f64::INFINITY` if `n` is too large to fit in `f64`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn factorial(n: u64) -> f64 {
    if n > 170 {
        return f64::INFINITY;
    }

    (1..=n).map(|i| i as f64).product()
}

/// Computes the number of permutations `P(n, k) = n! / (n-k)!`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of permutations as an `f64`. Returns `0.0` if `k > n`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn permutations(
    n: u64,
    k: u64,
) -> f64 {
    if k > n {
        return 0.0;
    }

    (n - k + 1..=n).map(|i| i as f64).product()
}

/// Computes the number of combinations `C(n, k) = n! / (k! * (n-k)!)`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of combinations as an `f64`. Returns `0.0` if `k > n`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn combinations(
    n: u64,
    k: u64,
) -> f64 {
    if k > n {
        return 0.0;
    }

    if k == 0 || k == n {
        return 1.0;
    }

    if k > n / 2 {
        return combinations(n, n - k);
    }

    let mut res = 1.0;

    for i in 1..=k {
        res = res * (n - i + 1) as f64 / i as f64;
    }

    res
}

/// Numerically solves a linear recurrence relation by unrolling it.
///
/// The recurrence relation is assumed to be of the form:
/// `a_n = coeffs[0]*a_{n-1} + coeffs[1]*a_{n-2} + ... + coeffs[order-1]*a_{n-order}`.
///
/// # Arguments
/// * `coeffs` - A slice of `f64` representing the coefficients of the recurrence relation.
/// * `initial_conditions` - A slice of `f64` representing the initial values `a_0, a_1, ..., a_{order-1}`.
/// * `target_n` - The index `n` for which to compute `a_n`.
///
/// # Returns
/// A `Result` containing the numerical value of `a_n`.
///
/// # Errors
///
/// Returns an error if:
/// - The number of initial conditions does not match the order of the recurrence.
/// - The internal storage fails to compute values (e.g., empty result vector).
pub fn solve_recurrence_numerical(
    coeffs: &[f64],
    initial_conditions: &[f64],
    target_n: usize,
) -> Result<f64, String> {
    let order = coeffs.len();

    if initial_conditions.len() != order {
        return Err("Number of initial \
             conditions must match \
             the order of the \
             recurrence."
            .to_string());
    }

    if target_n < order {
        return Ok(initial_conditions[target_n]);
    }

    let mut values = initial_conditions.to_vec();

    for n in order..=target_n {
        let mut next_val = 0.0;

        for i in 0..order {
            next_val += coeffs[i] * values[n - 1 - i];
        }

        values.push(next_val);
    }

    match values.last() {
        | Some(v) => Ok(*v),
        | None => {
            Err("Failed to compute \
                 the recurrence \
                 relation, values \
                 vector was empty."
                .to_string())
        },
    }
}

// B(218) fits in f64 and B(219) is larger than f64::MAX. Adding a
// singleton injects the partitions of n elements into those of n + 1.
const MAX_BELL_INDEX: usize = 218;
const SMALL_STIRLING_SIZE: usize = 41;

const fn stirling_row_start(n: usize) -> usize {
    n * (n + 1) / 2
}

#[allow(clippy::arithmetic_side_effects, clippy::indexing_slicing)]
static PARTITION_TABLES: (
    [f64; MAX_BELL_INDEX + 1],
    [f64; stirling_row_start(SMALL_STIRLING_SIZE)],
) = {
    let mut bell = [0.0; MAX_BELL_INDEX + 1];
    let mut stirling = [0.0; stirling_row_start(SMALL_STIRLING_SIZE)];
    let mut row = [0.0; MAX_BELL_INDEX + 1];
    row[0] = 1.0;
    bell[0] = 1.0;
    stirling[0] = 1.0;
    let mut n = 1;
    while n <= MAX_BELL_INDEX {
        let mut k = n;
        while k > 0 {
            row[k] = k as f64 * row[k] + row[k - 1];
            bell[n] += row[k];
            if n < SMALL_STIRLING_SIZE {
                stirling[stirling_row_start(n) + k] = row[k];
            }
            k -= 1;
        }
        row[0] = 0.0;
        n += 1;
    }
    (bell, stirling)
};

/// Computes the Stirling number of the second kind `S(n, k)` as an `f64`.
///
/// Counts partitions of `n` elements into `k` nonempty subsets. The empty
/// set has one partition into zero subsets. Returns zero for `k > n` and
/// for `k = 0 < n`, and positive infinity when the count overflows `f64`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn stirling_second(
    n: u64,
    k: u64,
) -> f64 {
    if k > n {
        return 0.0;
    }

    if k == 0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }

    if k == n || k == 1 {
        return 1.0;
    }

    if n < SMALL_STIRLING_SIZE as u64 {
        return PARTITION_TABLES.1[stirling_row_start(n as usize) + k as usize];
    }

    let excess = n - k;
    // Fix the first k elements in different blocks and freely assign the
    // rest to them. This injects k^(n-k) choices into the set of partitions.
    if excess.saturating_mul(u64::from(k.ilog2())) >= 1024 {
        return f64::INFINITY;
    }
    if excess == 1 {
        return (n as f64) * ((n - 1) as f64) * 0.5;
    }
    if k == 2 {
        return 2.0_f64.powi((n - 1) as i32) - 1.0;
    }
    if k <= 32 && n >= 4 * k {
        return stirling_few_blocks(n, k, excess);
    }
    if k > 2 * excess {
        return stirling_near_diagonal(n, excess as usize);
    }

    // Only intermediate partitions with at most k blocks and at most
    // n-k excess elements can contribute. Separate rows permit independent
    // updates within this band and one scratch allocation serves all rows.
    let k = k as usize;
    let excess = excess as usize;
    let stride = k + 1;
    let mut scratch = vec![0.0; 3 * stride];
    let (factors, rows) = scratch.split_at_mut(stride);
    let (mut old, mut next) = rows.split_at_mut(stride);
    for (blocks, scale) in factors.iter_mut().enumerate() {
        *scale = blocks as f64;
    }
    old[0] = 1.0;
    for elements in 1..=k + excess {
        let first = elements.saturating_sub(excess).max(1);
        let last = elements.min(k);
        for (((output, &scale), &same), &left) in next[first..=last]
            .iter_mut()
            .zip(&factors[first..=last])
            .zip(&old[first..=last])
            .zip(&old[first - 1..last])
        {
            *output = scale * same + left;
        }
        old[0] = 0.0;
        (old, next) = (next, old);
    }
    old[k]
}

#[allow(clippy::single_call_fn, clippy::cast_precision_loss)]
fn stirling_few_blocks(
    n: u64,
    k: u64,
    excess: u64,
) -> f64 {
    // Inclusion-exclusion on the probability that all blocks are used. The
    // probability is at least 1 - k*exp(-n/k) > 0.4 in this branch.
    let base = k as f64;
    let exponent = n as f64;
    let mut probability = 1.0;
    let mut choose = 1;
    for missed in 1..k {
        choose = choose * (k - missed + 1) / missed;
        let term = (choose as f64) * (((k - missed) as f64) / base).powf(exponent);
        // Terms decrease because every successive ratio is at most
        // k*exp(-n/k). The alternating tail is bounded by the first omission.
        if term <= f64::EPSILON / 8.0 {
            break;
        }
        if missed % 2 == 0 {
            probability += term;
        } else {
            probability -= term;
        }
    }

    // S(n,k) = probability * k^n/k!. Multiply by factors >= 1 only after
    // applying probability, so intermediate magnitudes stay bounded.
    let mut count = probability * base.powf(excess as f64);
    for divisor in 1..k {
        count *= base / divisor as f64;
    }
    count
}

#[allow(clippy::single_call_fn)]
fn stirling_near_diagonal(
    n: u64,
    excess: usize,
) -> f64 {
    // Since k > 2*excess, excess >= 128 would already have overflowed
    // the lower bound k^excess. At most 127 nonsingleton blocks are needed.
    let mut rows = [[0.0; 128]; 3];
    let [mut older, mut previous, mut row] = rows.each_mut();
    older[0] = 1.0;
    let exact_limit = 1_u128 << f64::MANTISSA_DIGITS;
    let mut exact_choose = (u128::from(n) <= exact_limit).then_some(u128::from(n));
    let mut choose = n as f64;
    let mut total = 0.0;

    for elements in 2..=2 * excess {
        let rest = elements - 1;
        let available = n - rest as u64;
        if let Some(exact) = exact_choose {
            // Keep binomial coefficients exact throughout the integer range
            // of f64. An exact numerator here fits in u128.
            let next = exact * u128::from(available) / elements as u128;
            exact_choose = (next <= exact_limit).then_some(next);
            choose = next as f64;
        } else {
            choose *= available as f64 / elements as f64;
        }
        // A partition without singletons gains its last element by adding
        // it to an existing nonsingleton block, or forming a new pair.
        row[0] = 0.0;
        for blocks in 1..=elements / 2 {
            row[blocks] = blocks as f64 * previous[blocks] + rest as f64 * older[blocks - 1];
        }
        if let Some(blocks) = elements.checked_sub(excess).filter(|&blocks| blocks > 0) {
            // Select the elements in nonsingleton blocks. All other blocks
            // are fixed singletons, so the contribution is nonnegative.
            total += choose * row[blocks];
            if total.is_infinite() {
                return f64::INFINITY;
            }
        }
        (older, previous, row) = (previous, row, older);
    }
    total
}

/// Computes the Bell number `B(n)` as an `f64`.
///
/// Counts all partitions of a set of `n` elements, including the one
/// partition of the empty set. Returns positive infinity when the count
/// overflows `f64` (`n >= 219`).
#[must_use]
pub fn bell(n: u64) -> f64 {
    if n <= MAX_BELL_INDEX as u64 {
        PARTITION_TABLES.0[n as usize]
    } else {
        f64::INFINITY
    }
}

/// Computes [`bell`] for each index using the shared table of finite counts.
///
/// Accepts unsorted indices and preserves their order and repetitions.
/// An empty input yields an empty output.
#[must_use]
pub fn bell_batch(indices: &[u64]) -> Vec<f64> {
    indices.iter().copied().map(bell).collect()
}

/// Computes the nth Catalan number `C_n`.
/// `C_n` = (1 / (n + 1)) * C(2n, n)
#[must_use]
pub fn catalan(n: u64) -> f64 {
    combinations(2 * n, n) / ((n + 1) as f64)
}

/// Computes the rising factorial (Pochhammer symbol) x^(n) = x(x+1)...(x+n-1).
#[must_use]
// precision loss is allowed here because we need to introduce bigint otherwise --- that will require a lot of work and results in breaking change and loss of performance. So we will have to handle it later.
#[allow(clippy::cast_precision_loss)]
pub fn rising_factorial(
    x: f64,
    n: u64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    let mut res = 1.0;

    for i in 0..n {
        res *= x + (i as f64);
    }

    res
}

/// Computes the falling factorial (x)_n = x(x-1)...(x-n+1).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn falling_factorial(
    x: f64,
    n: u64,
) -> f64 {
    if n == 0 {
        return 1.0;
    }

    let mut res = 1.0;

    for i in 0..n {
        res *= x - (i as f64);
    }

    res
}
