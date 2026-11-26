//! Exact linear algebra over the rationals.
//!
//! # Exact linear algebra over the rationals
//!
//! Dense matrices of [`BigRational`] stored as rows: reduced row echelon
//! forms, ranks, null spaces, solutions of linear systems, inverses,
//! determinants, characteristic polynomials, and integer normal forms (Hermite
//! and Smith invariant factors). Used by the Lie-algebra structure theory, by
//! the homology bases of simplicial complexes, and by symbolic matrix rules.
//!
//! ## Algorithms & Performance
//!
//! * **Sparse-Pivot Gauss-Jordan (`rref`)**: Pivot rows only normalize entries
//!   after the pivot column, extract non-zero `(col, val)` pairs, and update
//!   target rows over non-zero pivot columns only, with adaptive row
//!   parallelism via [`crate::execution::should_split_work`].
//! * **Denominator Clearing & Content Extraction**: Rational matrices are
//!   converted row-by-row to primitive integer matrices by clearing row LCMs
//!   and factoring out row GCDs before exact determinant or $p$-adic routines.
//! * **Exact Determinants (`det`, `det_bareiss`, `det_modular_int`)**:
//!   - $n \le 4$: closed-form / direct expansion.
//!   - $5 \le n \le 14$: fraction-free Bareiss elimination over $\mathbb{Z}$
//!     using exact integer division (zero rational GCD reductions).
//!   - $n > 14$: multi-prime 62-bit Montgomery `GF(p)` Gaussian elimination up
//!     to Hadamard's bound with parallel prime evaluation and Garner Chinese
//!     Remainder Theorem (CRT) reconstruction.
//! * **Dixon $p$-adic Lifting (`solve_dixon`, `solve`, `solve_mat`, `inverse`)**:
//!   For nonsingular square systems, inverts $A \bmod p$ once in 62-bit
//!   Montgomery arithmetic, lifts $p$-adic digits in $O(k n^2)$ modular
//!   matrix-vector products, reconstructs exact rationals via Wang's
//!   Euclidean rational reconstruction, and certifies $A X = B$ exactly.
//! * **Characteristic Polynomial (`charpoly`)**: Uses integer content-reduced
//!   Faddeev-LeVerrier for small orders ($n \le 10$) and parallel multi-prime
//!   62-bit Montgomery upper-Hessenberg similarity reduction + Garner CRT for
//!   larger matrices ($n > 10$).

#![allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    reason = "dense exact linear algebra and Montgomery modular word arithmetic enforce explicit matrix shape and bit-width bounds"
)]

use std::alloc::Layout;

use num_bigint::BigInt;
use num_bigint::Sign;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

/// An exact rational number.
pub type Q = BigRational;

/// A matrix stored as a list of rows.
pub type QMat = Vec<Vec<Q>>;

/// The integer `n` as a rational.
#[must_use]
pub fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

/// The zero vector of length `n`.
#[must_use]
pub fn zeros(n: usize) -> Vec<Q> {
    vec![Q::zero(); n]
}

/// The identity matrix of size `n`.
#[must_use]
pub fn identity(n: usize) -> QMat {
    (0..n)
        .map(|i| {
            let mut row = zeros(n);
            row[i] = Q::one();
            row
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 62-bit Montgomery modular arithmetic for multi-prime CRT & Dixon lifting
// ---------------------------------------------------------------------------

/// Verified 62-bit primes $< 2^{62}$ for Montgomery modular linear algebra.
const MONT_PRIMES_62: [u64; 32] = [
    0x3FFF_FFFF_FFFF_FFC7,
    0x3FFF_FFFF_FFFF_FFA9,
    0x3FFF_FFFF_FFFF_FF8B,
    0x3FFF_FFFF_FFFF_FF71,
    0x3FFF_FFFF_FFFF_FF67,
    0x3FFF_FFFF_FFFF_FF59,
    0x3FFF_FFFF_FFFF_FF55,
    0x3FFF_FFFF_FFFF_FF3D,
    0x3FFF_FFFF_FFFF_FF35,
    0x3FFF_FFFF_FFFF_FEEF,
    0x3FFF_FFFF_FFFF_FEE1,
    0x3FFF_FFFF_FFFF_FEC3,
    0x3FFF_FFFF_FFFF_FE45,
    0x3FFF_FFFF_FFFF_FE1D,
    0x3FFF_FFFF_FFFF_FE11,
    0x3FFF_FFFF_FFFF_FDC1,
    0x3FFF_FFFF_FFFF_FDBB,
    0x3FFF_FFFF_FFFF_FDA5,
    0x3FFF_FFFF_FFFF_FD87,
    0x3FFF_FFFF_FFFF_FD69,
    0x3FFF_FFFF_FFFF_FD03,
    0x3FFF_FFFF_FFFF_FCFB,
    0x3FFF_FFFF_FFFF_FCF7,
    0x3FFF_FFFF_FFFF_FCE9,
    0x3FFF_FFFF_FFFF_FCD3,
    0x3FFF_FFFF_FFFF_FCC1,
    0x3FFF_FFFF_FFFF_FC65,
    0x3FFF_FFFF_FFFF_FC2B,
    0x3FFF_FFFF_FFFF_FC1F,
    0x3FFF_FFFF_FFFF_FC17,
    0x3FFF_FFFF_FFFF_FC11,
    0x3FFF_FFFF_FFFF_FC07,
];

/// Montgomery arithmetic for a modulus the caller has certified as an odd
/// prime; the implementation is [`crate::kernels::modular::Montgomery`].
use crate::kernels::modular::Montgomery;

/// `a^(-1)` in Montgomery form.
///
/// Fermat's little theorem with the prime exponent `p - 2`: cheaper than an
/// extended GCD here because the base is already in the Montgomery domain, so
/// the whole exponentiation is `mul_mont` and needs no conversion.
#[inline]
#[must_use]
fn mont_inverse(
    mont: Montgomery,
    a_mont: u64,
) -> u64 {
    mont.pow_mont(a_mont, mont.modulus() - 2)
}

/// Deterministic Miller-Rabin primality test; see
/// [`crate::kernels::number_theory::is_prime_miller_rabin`].
#[must_use]
fn is_prime_u64(n: u64) -> bool {
    crate::kernels::number_theory::is_prime_miller_rabin(n)
}

/// Generates `count` distinct 62-bit primes.
#[must_use]
fn select_primes_62(count: usize) -> Vec<u64> {
    let mut out = Vec::with_capacity(count);
    let take_fixed = count.min(MONT_PRIMES_62.len());
    out.extend_from_slice(&MONT_PRIMES_62[..take_fixed]);
    if out.len() < count {
        let mut candidate = MONT_PRIMES_62[MONT_PRIMES_62.len() - 1] - 2;
        while out.len() < count {
            if is_prime_u64(candidate) {
                out.push(candidate);
            }
            candidate -= 2;
        }
    }
    out
}

/// Reduces a signed [`BigInt`] modulo `p` into `[0, p)`.
#[inline]
#[must_use]
fn bigint_mod_u64(
    val: &BigInt,
    p: u64,
) -> u64 {
    if val.is_zero() {
        return 0;
    }
    if let Some(v) = val.to_i64() {
        let r = v.rem_euclid(p as i64);
        return r as u64;
    }
    let (sign, digits) = val.to_u64_digits();
    let p128 = p as u128;
    let mut rem = 0_u128;
    for &limb in digits.iter().rev() {
        rem = ((rem << 64) | (limb as u128)) % p128;
    }
    let r = rem as u64;
    if sign == Sign::Minus && r != 0 {
        p - r
    } else {
        r
    }
}

/// Modular inverse of `a` modulo `m`.
///
/// Delegates to [`crate::kernels::modular::mod_inv_u64`], with the historical
/// convention that `m == 1` has the unique residue `0` as its inverse.
#[inline]
#[must_use]
const fn mod_inv_u64(
    a: u64,
    m: u64,
) -> Option<u64> {
    if m == 1 {
        return Some(0);
    }
    crate::kernels::modular::mod_inv_u64(a, m)
}

/// Gaussian elimination in `GF(p)` using Montgomery arithmetic on a flat $n \times n$
/// buffer `a_mod` (values in `[0, p)` cast to `i64`, already in Montgomery domain).
/// Returns `det(A) mod p` in standard `[0, p)` representation.
#[must_use]
fn det_mod_montgomery(
    a_mod: &mut [i64],
    n: usize,
    mont: Montgomery,
) -> u64 {
    if n == 0 {
        return 1;
    }
    let mut det_mont = mont.one();
    for c in 0..n {
        let Some(p_row) = (c..n).find(|&r| a_mod[r * n + c] != 0) else {
            return 0;
        };
        if p_row != c {
            for j in c..n {
                a_mod.swap(c * n + j, p_row * n + j);
            }
            det_mont = mont.neg(det_mont);
        }
        let pivot = a_mod[c * n + c].cast_unsigned();
        det_mont = mont.mul_mont(det_mont, pivot);
        let inv_pivot = mont_inverse(mont, pivot);
        let (head, tail) = a_mod.split_at_mut((c + 1) * n);
        let pivot_row = &head[c * n + c + 1..(c + 1) * n];
        for row in tail.chunks_exact_mut(n) {
            let entry = row[c].cast_unsigned();
            if entry == 0 {
                continue;
            }
            let factor = mont.mul_mont(entry, inv_pivot);
            row[c] = 0;
            for (x, &y) in row[c + 1..].iter_mut().zip(pivot_row) {
                if y != 0 {
                    *x = mont
                        .sub(x.cast_unsigned(), mont.mul_mont(factor, y.cast_unsigned()))
                        .cast_signed();
                }
            }
        }
    }
    mont.from_mont(det_mont)
}

/// Rank of an $m \times n$ matrix in `GF(p)` using Montgomery arithmetic on `a_mod`
/// (entries in `[0, p)` stored as `i64` in Montgomery domain).
#[must_use]
fn rank_mod_montgomery(
    a_mod: &mut [i64],
    rows: usize,
    cols: usize,
    mont: Montgomery,
) -> usize {
    if rows == 0 || cols == 0 {
        return 0;
    }
    let mut r = 0_usize;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p_row) = (r..rows).find(|&i| a_mod[i * cols + c] != 0) else {
            continue;
        };
        if p_row != r {
            for j in c..cols {
                a_mod.swap(r * cols + j, p_row * cols + j);
            }
        }
        let pivot = a_mod[r * cols + c].cast_unsigned();
        let inv_pivot = mont_inverse(mont, pivot);
        let (head, tail) = a_mod.split_at_mut((r + 1) * cols);
        let pivot_row = &head[r * cols + c + 1..(r + 1) * cols];
        for row in tail.chunks_exact_mut(cols) {
            let entry = row[c].cast_unsigned();
            if entry == 0 {
                continue;
            }
            let factor = mont.mul_mont(entry, inv_pivot);
            row[c] = 0;
            for (x, &y) in row[c + 1..].iter_mut().zip(pivot_row) {
                if y != 0 {
                    *x = mont
                        .sub(x.cast_unsigned(), mont.mul_mont(factor, y.cast_unsigned()))
                        .cast_signed();
                }
            }
        }
        r += 1;
    }
    r
}

/// Inverts an $n \times n$ integer matrix modulo `mont.modulus()` in Montgomery arithmetic.
/// Returns the inverse matrix in standard `[0, p)` representation as a flat $n \times n$ vector,
/// or `None` if `A` is singular modulo `p`.
#[must_use]
fn inv_mod_montgomery(
    a: &[Vec<BigInt>],
    n: usize,
    mont: Montgomery,
) -> Option<Vec<u64>> {
    let width = 2 * n;
    let mut aug = crate::kernels::reuse::global_pool().scratch_i64(n * width);
    aug.resize(n * width, 0);
    for (i, row) in a.iter().enumerate() {
        for (j, val) in row.iter().enumerate() {
            aug[i * width + j] = mont
                .to_mont(bigint_mod_u64(val, mont.modulus()))
                .cast_signed();
        }
        aug[i * width + n + i] = mont.one().cast_signed();
    }
    let mut pivot_tail = crate::kernels::reuse::global_pool().scratch_i64(width);
    for c in 0..n {
        let p_row = (c..n).find(|&r| aug[r * width + c] != 0)?;
        if p_row != c {
            for j in c..width {
                aug.swap(c * width + j, p_row * width + j);
            }
        }
        let inv_pivot = mont_inverse(mont, aug[c * width + c].cast_unsigned());
        aug[c * width + c] = mont.one().cast_signed();
        for j in c + 1..width {
            let v = mont.mul_mont(aug[c * width + j].cast_unsigned(), inv_pivot);
            aug[c * width + j] = v.cast_signed();
        }
        pivot_tail.clear();
        pivot_tail.extend_from_slice(&aug[c * width + c + 1..(c + 1) * width]);
        for r in 0..n {
            if r == c {
                continue;
            }
            let factor = aug[r * width + c].cast_unsigned();
            if factor == 0 {
                continue;
            }
            aug[r * width + c] = 0;
            let row_slice = &mut aug[r * width + c + 1..(r + 1) * width];
            for (x, &y) in row_slice.iter_mut().zip(pivot_tail.iter()) {
                if y != 0 {
                    *x = mont
                        .sub(x.cast_unsigned(), mont.mul_mont(factor, y.cast_unsigned()))
                        .cast_signed();
                }
            }
        }
    }
    let mut inv = vec![0_u64; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = mont.from_mont(aug[i * width + n + j].cast_unsigned());
        }
    }
    Some(inv)
}

/// Computes the characteristic polynomial $\det(\lambda I - A) \bmod p$ of an $n \times n$
/// matrix stored in `a_mod` (Montgomery domain) via similarity reduction to upper
/// Hessenberg form in $O(n^3)$ `u64` operations. Returns `[c_0, ..., c_{n-1}]` in `[0, p)`.
#[must_use]
fn charpoly_hessenberg_mod(
    a_mod: &mut [i64],
    n: usize,
    mont: Montgomery,
) -> Vec<u64> {
    if n == 0 {
        return Vec::new();
    }
    for c in 0..n.saturating_sub(2) {
        let Some(r) = (c + 1..n).find(|&i| a_mod[i * n + c] != 0) else {
            continue;
        };
        if a_mod[(c + 1) * n + c] == 0 {
            for j in c..n {
                a_mod.swap((c + 1) * n + j, r * n + j);
            }
            for row in 0..n {
                a_mod.swap(row * n + c + 1, row * n + r);
            }
        }
        let piv = a_mod[(c + 1) * n + c].cast_unsigned();
        let inv_piv = mont_inverse(mont, piv);
        for i in c + 2..n {
            let entry = a_mod[i * n + c].cast_unsigned();
            if entry == 0 {
                continue;
            }
            let alpha = mont.mul_mont(entry, inv_piv);
            a_mod[i * n + c] = 0;
            for j in c + 1..n {
                let h_piv = a_mod[(c + 1) * n + j].cast_unsigned();
                if h_piv != 0 {
                    let cur = a_mod[i * n + j].cast_unsigned();
                    a_mod[i * n + j] = mont.sub(cur, mont.mul_mont(alpha, h_piv)).cast_signed();
                }
            }
            for row in 0..n {
                let h_ri = a_mod[row * n + i].cast_unsigned();
                if h_ri != 0 {
                    let cur = a_mod[row * n + c + 1].cast_unsigned();
                    a_mod[row * n + c + 1] =
                        mont.add(cur, mont.mul_mont(alpha, h_ri)).cast_signed();
                }
            }
        }
    }

    let mut polys: Vec<Vec<u64>> = Vec::with_capacity(n + 1);
    polys.push(vec![mont.one()]);
    for m in 0..n {
        let mut next_p = vec![0_u64; m + 2];
        {
            let pm = &polys[m];
            next_p[1..].copy_from_slice(pm);
            let hmm = a_mod[m * n + m].cast_unsigned();
            if hmm != 0 {
                for (deg, &coeff) in pm.iter().enumerate() {
                    if coeff != 0 {
                        next_p[deg] = mont.sub(next_p[deg], mont.mul_mont(hmm, coeff));
                    }
                }
            }
        }
        let mut prod = mont.one();
        for i in (0..m).rev() {
            let sub_entry = a_mod[(i + 1) * n + i].cast_unsigned();
            if sub_entry == 0 {
                break;
            }
            prod = mont.mul_mont(prod, sub_entry);
            let him = a_mod[i * n + m].cast_unsigned();
            if him == 0 {
                continue;
            }
            let factor = mont.mul_mont(him, prod);
            for (deg, &coeff) in polys[i].iter().enumerate() {
                if coeff != 0 {
                    next_p[deg] = mont.sub(next_p[deg], mont.mul_mont(factor, coeff));
                }
            }
        }
        polys.push(next_p);
    }
    polys[n][..n].iter().map(|&c| mont.from_mont(c)).collect()
}

/// Reconstructs a unique signed [`BigInt`] in `(-M/2, M/2)` from its residues
/// `residues[i] = x mod primes[i]` using Garner's mixed-radix CRT algorithm.
#[must_use]
fn garner_crt_signed(
    primes: &[u64],
    residues: &[u64],
) -> BigInt {
    let k = primes.len();
    if k == 0 {
        return BigInt::zero();
    }
    // Mixed-radix digits v[i]
    let mut v = vec![0_u64; k];
    for i in 0..k {
        let pi = primes[i];
        let mut acc = residues[i];
        let p128 = pi as u128;
        for j in 0..i {
            let pj_mod_pi = primes[j] % pi;
            let inv = mod_inv_u64(pj_mod_pi, pi).unwrap_or(1);
            let diff = if acc >= v[j] % pi {
                acc - (v[j] % pi)
            } else {
                acc + pi - (v[j] % pi)
            };
            acc = (((diff as u128) * (inv as u128)) % p128) as u64;
        }
        v[i] = acc;
    }
    let mut x = BigInt::zero();
    let mut prod = BigInt::one();
    for i in 0..k {
        if v[i] != 0 {
            x += &prod * BigInt::from(v[i]);
        }
        prod *= BigInt::from(primes[i]);
    }
    let half = &prod >> 1_usize;
    if x > half { x - prod } else { x }
}

// ---------------------------------------------------------------------------
// Denominator clearing & Row content extraction
// ---------------------------------------------------------------------------

/// Clears denominators of each row of `m`, returning `(int_rows, row_scales)`
/// where `m[i][j] = int_rows[i][j] * row_scales[i]`.
/// Also factors out the positive integer GCD of each non-zero row so the integer
/// entries have minimum bit-length.
#[must_use]
fn clear_rows_to_primitive_ints(m: &[Vec<Q>]) -> (Vec<Vec<BigInt>>, Vec<Q>) {
    let mut int_rows = Vec::with_capacity(m.len());
    let mut scales = Vec::with_capacity(m.len());
    for row in m {
        let mut lcm_den = BigInt::one();
        let mut all_int = true;
        for entry in row {
            if !entry.is_zero() && !entry.denom().is_one() {
                all_int = false;
                lcm_den = lcm_den.lcm(entry.denom());
            }
        }
        let mut int_row: Vec<BigInt> = if all_int {
            row.iter().map(|e| e.numer().clone()).collect()
        } else {
            row.iter()
                .map(|e| {
                    if e.is_zero() {
                        BigInt::zero()
                    } else {
                        e.numer() * (&lcm_den / e.denom())
                    }
                })
                .collect()
        };
        let mut row_gcd = BigInt::zero();
        for val in &int_row {
            if !val.is_zero() {
                row_gcd = if row_gcd.is_zero() {
                    val.abs()
                } else {
                    row_gcd.gcd(val)
                };
                if row_gcd.is_one() {
                    break;
                }
            }
        }
        if row_gcd > BigInt::one() {
            for val in &mut int_row {
                if !val.is_zero() {
                    *val /= &row_gcd;
                }
            }
            scales.push(Q::new(row_gcd, lcm_den));
        } else {
            scales.push(Q::new(BigInt::one(), lcm_den));
        }
        int_rows.push(int_row);
    }
    (int_rows, scales)
}

// ---------------------------------------------------------------------------
// Sparse-Pivot Gauss-Jordan RREF
// ---------------------------------------------------------------------------

/// Adaptive grain size for row-parallel elimination over `rows` rows with
/// `active_work` operations per row.
#[inline]
#[must_use]
fn adaptive_row_grain(
    rows: usize,
    active_work: usize,
) -> usize {
    if crate::execution::should_split_work(rows, active_work.max(1).saturating_mul(64))
        || crate::execution::should_split(rows, 48)
    {
        crate::execution::choose_chunk_size(rows, 8)
    } else {
        usize::MAX
    }
}

/// Brings `m` (rows of length `cols`) to reduced row echelon form in place,
/// drops the zero rows and returns the pivot columns.
pub fn rref(
    m: &mut QMat,
    cols: usize,
) -> Vec<usize> {
    let rows = m.len();
    let mut pivots = Vec::with_capacity(rows.min(cols));
    let mut r = 0;
    let mut sparse_pivot: Vec<(usize, Q)> = Vec::new();
    let minus_one = q(-1);
    for c in 0..cols {
        if r == m.len() {
            break;
        }
        // Choose the pivot row in r..m.len() with non-zero column c that has the
        // fewest non-zero entries in c + 1..cols (Markowitz-style sparsity preference
        // while keeping deterministic earliest tie-break).
        let mut best_p: Option<(usize, usize)> = None;
        for (i, row) in m.iter().enumerate().skip(r) {
            if !row[c].is_zero() {
                let nnz = row[c + 1..cols].iter().filter(|x| !x.is_zero()).count();
                if best_p.is_none_or(|(_, best_nnz)| nnz < best_nnz) {
                    best_p = Some((i, nnz));
                    if nnz == 0 {
                        break;
                    }
                }
            }
        }
        let Some((p, _)) = best_p else {
            continue;
        };
        m.swap(r, p);
        let pivot_val = std::mem::replace(&mut m[r][c], Q::one());
        sparse_pivot.clear();
        if pivot_val.is_one() {
            for (j, x) in m[r][c + 1..cols].iter().enumerate() {
                if !x.is_zero() {
                    sparse_pivot.push((c + 1 + j, x.clone()));
                }
            }
        } else if pivot_val == minus_one {
            for (j, x) in m[r][c + 1..cols].iter_mut().enumerate() {
                if !x.is_zero() {
                    *x = -std::mem::take(x);
                    sparse_pivot.push((c + 1 + j, x.clone()));
                }
            }
        } else {
            let inv = pivot_val.recip();
            for (j, x) in m[r][c + 1..cols].iter_mut().enumerate() {
                if !x.is_zero() {
                    *x = &*x * &inv;
                    sparse_pivot.push((c + 1 + j, x.clone()));
                }
            }
        }
        if sparse_pivot.is_empty() {
            for (i, row) in m.iter_mut().enumerate() {
                if i != r && !row[c].is_zero() {
                    row[c] = Q::zero();
                }
            }
        } else {
            let grain = adaptive_row_grain(m.len(), sparse_pivot.len());
            let pivot_pairs = sparse_pivot.as_slice();
            super::par::for_each_mut(m, grain, |i, row| {
                if i == r || row[c].is_zero() {
                    return;
                }
                let f = std::mem::replace(&mut row[c], Q::zero());
                if f.is_one() {
                    for &(col_idx, ref y) in pivot_pairs {
                        row[col_idx] -= y;
                    }
                } else if f == minus_one {
                    for &(col_idx, ref y) in pivot_pairs {
                        row[col_idx] += y;
                    }
                } else {
                    for &(col_idx, ref y) in pivot_pairs {
                        row[col_idx] -= &f * y;
                    }
                }
            });
        }
        pivots.push(c);
        r += 1;
    }
    m.truncate(r);
    pivots
}

/// The rank of a matrix whose rows have length `cols`.
#[must_use]
pub fn rank(
    m: &[Vec<Q>],
    cols: usize,
) -> usize {
    let rows = m.len();
    if rows == 0 || cols == 0 {
        return 0;
    }
    // Fast exact rank certification for medium/large full-rank matrices via 2 modular primes
    let min_dim = rows.min(cols);
    if min_dim >= 16 && m.iter().all(|r| r.len() >= cols) {
        let (int_rows, _) = clear_rows_to_primitive_ints(m);
        let p0 = MONT_PRIMES_62[0];
        let p1 = MONT_PRIMES_62[1];
        let mont0 = Montgomery::new_prime(p0);
        let mont1 = Montgomery::new_prime(p1);
        let mut buf = crate::kernels::reuse::global_pool().scratch_i64(rows * cols);
        buf.resize(rows * cols, 0);
        for (i, row) in int_rows.iter().enumerate() {
            for (j, val) in row.iter().take(cols).enumerate() {
                buf[i * cols + j] = mont0.to_mont(bigint_mod_u64(val, p0)).cast_signed();
            }
        }
        let r0 = rank_mod_montgomery(&mut buf, rows, cols, mont0);
        if r0 == min_dim {
            return min_dim;
        }
        for (i, row) in int_rows.iter().enumerate() {
            for (j, val) in row.iter().take(cols).enumerate() {
                buf[i * cols + j] = mont1.to_mont(bigint_mod_u64(val, p1)).cast_signed();
            }
        }
        let r1 = rank_mod_montgomery(&mut buf, rows, cols, mont1);
        if r1 == min_dim {
            return min_dim;
        }
        if r0 == 0
            && r1 == 0
            && int_rows
                .iter()
                .all(|r| r.iter().take(cols).all(Zero::is_zero))
        {
            return 0;
        }
    }
    let mut a = m.to_vec();
    rref(&mut a, cols).len()
}

/// A basis (reduced echelon rows) of the span of `vectors` of length `dim`.
#[must_use]
pub fn span_basis(
    vectors: &[Vec<Q>],
    dim: usize,
) -> QMat {
    let mut a = vectors.to_vec();
    rref(&mut a, dim);
    a
}

/// A basis of `{x : m x = 0}` for `m` with `cols` columns.
#[must_use]
pub fn nullspace(
    m: &[Vec<Q>],
    cols: usize,
) -> QMat {
    let mut a = m.to_vec();
    let pivots = rref(&mut a, cols);
    let mut is_pivot = vec![false; cols];
    for &p in &pivots {
        if p < cols {
            is_pivot[p] = true;
        }
    }
    let mut out = Vec::with_capacity(cols.saturating_sub(pivots.len()));
    for (free, &piv) in is_pivot.iter().enumerate() {
        if piv {
            continue;
        }
        let mut x = zeros(cols);
        x[free] = Q::one();
        for (row, &p) in a.iter().zip(&pivots) {
            if !row[free].is_zero() {
                x[p] = -row[free].clone();
            }
        }
        out.push(x);
    }
    out
}

// ---------------------------------------------------------------------------
// Wang's Rational Reconstruction & Dixon p-adic Lifting
// ---------------------------------------------------------------------------

/// Reconstructs a canonical rational `n / d` from `u = n * d^{-1} mod m`
/// using Wang's Euclidean rational reconstruction algorithm with bounds
/// `|n| <= n_bound` and `1 <= d <= d_bound`.
#[must_use]
fn reconstruct_rational_bounded(
    u: &BigInt,
    m: &BigInt,
    n_bound: &BigInt,
    d_bound: &BigInt,
) -> Option<Q> {
    let u_rem = u.mod_floor(m);
    if u_rem.is_zero() {
        return Some(Q::zero());
    }
    let mut r0 = m.clone();
    let mut r1 = u_rem;
    let mut t0 = BigInt::zero();
    let mut t1 = BigInt::one();

    while &r1 > n_bound {
        let (q_div, r_next) = r0.div_rem(&r1);
        r0 = r1;
        r1 = r_next;
        let t_next = t0 - &q_div * &t1;
        t0 = t1;
        t1 = t_next;
        if t1.abs() > *d_bound {
            return None;
        }
    }
    if t1.is_zero() || t1.abs() > *d_bound {
        return None;
    }
    let (num, den) = if t1.is_negative() {
        (-r1, -t1)
    } else {
        (r1, t1)
    };
    if num.gcd(&den).is_one() {
        Some(Q::new_raw(num, den))
    } else {
        Some(Q::new(num, den))
    }
}

/// Solves a nonsingular square system $A x = b$ over $\mathbb{Q}$ using Dixon's
/// $p$-adic lifting over a 62-bit prime $p$ and exact residual certification.
///
/// Returns `None` if `a` is not square, dimensions mismatch, or `a` is singular.
#[must_use]
pub fn solve_dixon(
    a: &[Vec<Q>],
    b: &[Q],
) -> Option<Vec<Q>> {
    let n = a.len();
    if n != b.len() || a.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }
    if b.iter().all(Zero::is_zero) {
        // Verify A is nonsingular before returning 0
        return (!det(a).is_zero()).then(|| zeros(n));
    }
    let cols_b: Vec<Vec<Q>> = b.iter().map(|bi| vec![bi.clone()]).collect();
    let sol_mat = solve_dixon_dense_multi(a, &cols_b)?;
    Some(
        sol_mat
            .into_iter()
            .map(|mut row| row.pop().unwrap_or_else(Q::zero))
            .collect(),
    )
}

/// Multi-RHS Dixon $p$-adic lifting solver for $A X = B$ where $A \in \mathbb{Q}^{n \times n}$
/// and $B \in \mathbb{Q}^{n \times r}$.
///
/// Clears each row $[A_i \mid B_i]$ by its common denominator LCM to obtain an integer
/// system $\tilde{A} X = \tilde{B}$, inverts $\tilde{A} \bmod p$ once in $O(n^3)$
/// 62-bit Montgomery arithmetic, lifts $p$-adic digits, reconstructs $X \in \mathbb{Q}^{n \times r}$,
/// and certifies $\tilde{A} X = \tilde{B}$ exactly.
#[must_use]
fn solve_dixon_dense_multi(
    a: &[Vec<Q>],
    b: &[Vec<Q>],
) -> Option<QMat> {
    let n = a.len();
    let rhs_cols = b.first().map_or(0, Vec::len);
    if n == 0
        || b.len() != n
        || a.iter().any(|r| r.len() != n)
        || b.iter().any(|r| r.len() != rhs_cols)
    {
        return None;
    }
    if rhs_cols == 0 {
        return Some(vec![Vec::new(); n]);
    }

    // Clear denominators of each row [A_i | B_i] jointly so A_int * X = B_int has the exact same X!
    let mut a_int = Vec::with_capacity(n);
    let mut b_int = Vec::with_capacity(n);
    let mut max_row_bits = 1_usize;
    for (row_a, row_b) in a.iter().zip(b) {
        let mut lcm_den = BigInt::one();
        for x in row_a.iter().chain(row_b.iter()) {
            if !x.is_zero() && !x.denom().is_one() {
                lcm_den = lcm_den.lcm(x.denom());
            }
        }
        let scale_entry = |x: &Q| -> BigInt {
            if x.is_zero() {
                BigInt::zero()
            } else if lcm_den.is_one() {
                x.numer().clone()
            } else {
                x.numer() * (&lcm_den / x.denom())
            }
        };
        let mut ra: Vec<BigInt> = row_a.iter().map(scale_entry).collect();
        let mut rb: Vec<BigInt> = row_b.iter().map(scale_entry).collect();
        // Cancel common row GCD across [ra | rb]
        let mut g = BigInt::zero();
        for val in ra.iter().chain(rb.iter()) {
            if !val.is_zero() {
                g = if g.is_zero() {
                    val.abs()
                } else {
                    g.gcd(val)
                };
                if g.is_one() {
                    break;
                }
            }
        }
        if g > BigInt::one() {
            for val in ra.iter_mut().chain(rb.iter_mut()) {
                if !val.is_zero() {
                    *val /= &g;
                }
            }
        }
        for val in ra.iter().chain(rb.iter()) {
            max_row_bits = max_row_bits.max(val.bits() as usize);
        }
        a_int.push(ra);
        b_int.push(rb);
    }

    // Pick a 62-bit prime for which A_int mod p is invertible.
    // Test p0 and p1 via inversion; if both fail, probe remaining primes with a
    // cheaper determinant check before invoking full matrix inversion.
    let mut chosen: Option<(Montgomery, Vec<u64>)> = None;
    for &p in &MONT_PRIMES_62[..2] {
        let mont = Montgomery::new_prime(p);
        if let Some(inv_mod_p) = inv_mod_montgomery(&a_int, n, mont) {
            chosen = Some((mont, inv_mod_p));
            break;
        }
    }
    if chosen.is_none() {
        let mut det_buf = crate::kernels::reuse::global_pool().scratch_i64(n * n);
        det_buf.resize(n * n, 0);
        for &p in &MONT_PRIMES_62[2..8] {
            let mont = Montgomery::new_prime(p);
            for (i, row) in a_int.iter().enumerate() {
                let base = i * n;
                for (j, val) in row.iter().enumerate() {
                    det_buf[base + j] = mont.to_mont(bigint_mod_u64(val, p)).cast_signed();
                }
            }
            if det_mod_montgomery(&mut det_buf, n, mont) != 0
                && let Some(inv_mod_p) = inv_mod_montgomery(&a_int, n, mont)
            {
                chosen = Some((mont, inv_mod_p));
                break;
            }
        }
    }
    let (mont, a_inv_mod_p) = chosen?;
    let p = mont.modulus();
    let p_big = BigInt::from(p);

    // Hadamard bound on bits of Cramér determinants: 2 * n * (max_row_bits + log2(n)/2)
    let log2_n = (usize::BITS - n.leading_zeros()) as usize;
    let cramer_bits = 2 * n * (max_row_bits + log2_n + 2) + 64;
    let max_steps = cramer_bits.div_ceil(61).max(2);

    // Transpose B_int into columns for cache-friendly lifting: residuals[col][row]
    let mut residuals: Vec<Vec<BigInt>> = (0..rhs_cols)
        .map(|col| (0..n).map(|row| b_int[row][col].clone()).collect())
        .collect();
    let mut x_acc: Vec<Vec<BigInt>> = (0..rhs_cols).map(|_| vec![BigInt::zero(); n]).collect();
    let mut p_pow = BigInt::one();

    // Pre-extract sparse/compact i64 view of A_int if entries fit in i64 for ultra-fast residual update
    let fits_i64 = max_row_bits <= 62;
    let a_i64: Vec<Vec<(usize, i64)>> = if fits_i64 {
        a_int
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter_map(|(j, val)| {
                        if val.is_zero() {
                            None
                        } else {
                            val.to_i64().map(|v| (j, v))
                        }
                    })
                    .collect()
            })
            .collect()
    } else {
        Vec::new()
    };

    let check_interval = if n <= 16 {
        1
    } else if n <= 48 {
        2
    } else {
        4
    };

    for step in 0..max_steps {
        let mut all_zero_residual = true;
        // For each RHS column, compute d_k = A_inv_mod_p * (r_k mod p) mod p,
        // update x_acc += d_k * p^k, and r_{k+1} = (r_k - A_int * d_k) / p.
        let grain = if crate::execution::should_split_work(rhs_cols, n * n) {
            crate::execution::choose_chunk_size(rhs_cols, 1)
        } else {
            usize::MAX
        };
        let p_pow_ref = &p_pow;
        let a_inv_ref = &a_inv_mod_p;
        let a_int_ref = &a_int;
        let a_i64_ref = &a_i64;
        let p_big_ref = &p_big;

        let mut col_pairs: Vec<(Vec<BigInt>, Vec<BigInt>)> =
            residuals.into_iter().zip(x_acc).collect();

        super::par::for_each_mut(&mut col_pairs, grain, |_, (res_col, x_col)| {
            if res_col.iter().all(Zero::is_zero) {
                return;
            }
            let mut work_buf = crate::kernels::reuse::global_pool().scratch_i64(2 * n);
            work_buf.resize(2 * n, 0);
            let (r_mod_p, d_k) = work_buf.split_at_mut(n);
            for (dst, v) in r_mod_p.iter_mut().zip(res_col.iter()) {
                *dst = bigint_mod_u64(v, p).cast_signed();
            }
            for i in 0..n {
                let inv_row = &a_inv_ref[i * n..(i + 1) * n];
                let mut sum = 0_u128;
                for (&a_ij, &r_j) in inv_row.iter().zip(r_mod_p.iter()) {
                    if a_ij != 0 && r_j != 0 {
                        sum += (a_ij as u128) * (r_j.cast_unsigned() as u128);
                        // Prevent u128 overflow when n > 4
                        if sum >= (1_u128 << 124) {
                            sum %= p as u128;
                        }
                    }
                }
                d_k[i] = ((sum % (p as u128)) as u64).cast_signed();
            }
            for (xi, &di_signed) in x_col.iter_mut().zip(d_k.iter()) {
                let di = di_signed.cast_unsigned();
                if di != 0 {
                    *xi += p_pow_ref * BigInt::from(di);
                }
            }
            if fits_i64 {
                for (i, row_sparse) in a_i64_ref.iter().enumerate() {
                    let mut prod = BigInt::zero();
                    let mut acc_i128 = 0_i128;
                    let mut count = 0_usize;
                    for &(j, a_ij) in row_sparse {
                        let dj = d_k[j];
                        if dj != 0 {
                            acc_i128 += (a_ij as i128) * (dj as i128);
                            count += 1;
                            // 4 * 2^124 = 2^126 < i128::MAX
                            if count == 4 {
                                prod += BigInt::from(acc_i128);
                                acc_i128 = 0;
                                count = 0;
                            }
                        }
                    }
                    if acc_i128 != 0 {
                        prod += BigInt::from(acc_i128);
                    }
                    res_col[i] = (&res_col[i] - prod) / p_big_ref;
                }
            } else {
                for (i, row_a) in a_int_ref.iter().enumerate() {
                    let mut prod = BigInt::zero();
                    for (a_ij, &dj_signed) in row_a.iter().zip(d_k.iter()) {
                        let dj = dj_signed.cast_unsigned();
                        if dj != 0 && !a_ij.is_zero() {
                            prod += a_ij * BigInt::from(dj);
                        }
                    }
                    res_col[i] = (&res_col[i] - prod) / p_big_ref;
                }
            }
        });

        residuals = Vec::with_capacity(rhs_cols);
        x_acc = Vec::with_capacity(rhs_cols);
        for (r_col, x_col) in col_pairs {
            if !r_col.iter().all(Zero::is_zero) {
                all_zero_residual = false;
            }
            residuals.push(r_col);
            x_acc.push(x_col);
        }

        p_pow *= &p_big;

        let should_try_reconstruct =
            all_zero_residual || ((step + 1) % check_interval == 0) || step + 1 == max_steps;

        if should_try_reconstruct
            && let Some(certified) = try_reconstruct_and_certify(
                &a_int, &b_int, &a_i64, fits_i64, &x_acc, &p_pow, n, rhs_cols,
            )
        {
            return Some(certified);
        }
    }
    None
}

/// Attempts Wang rational reconstruction on every column of `x_acc` modulo `p_pow`
/// and certifies $A_{\text{int}} X = B_{\text{int}}$ exactly by clearing the reconstructed
/// column denominator.
#[allow(clippy::too_many_arguments)]
#[must_use]
fn try_reconstruct_and_certify(
    a_int: &[Vec<BigInt>],
    b_int: &[Vec<BigInt>],
    a_i64: &[Vec<(usize, i64)>],
    fits_i64: bool,
    x_acc: &[Vec<BigInt>],
    p_pow: &BigInt,
    n: usize,
    rhs_cols: usize,
) -> Option<QMat> {
    let bound = (p_pow >> 1_usize).sqrt();
    let mut x_cols: Vec<Vec<Q>> = Vec::with_capacity(rhs_cols);
    for col_idx in 0..rhs_cols {
        let mut q_col = Vec::with_capacity(n);
        let mut col_lcm = BigInt::one();
        for i in 0..n {
            let qi = reconstruct_rational_bounded(&x_acc[col_idx][i], p_pow, &bound, &bound)?;
            if !qi.is_zero() && !qi.denom().is_one() {
                col_lcm = col_lcm.lcm(qi.denom());
            }
            q_col.push(qi);
        }
        // Clear column denominator: x_int[i] = q_col[i] * col_lcm, then verify A_int * x_int == B_int[:, col] * col_lcm
        let x_int: Vec<BigInt> = q_col
            .iter()
            .map(|qi| {
                if qi.is_zero() {
                    BigInt::zero()
                } else if col_lcm.is_one() {
                    qi.numer().clone()
                } else {
                    qi.numer() * (&col_lcm / qi.denom())
                }
            })
            .collect();

        for i in 0..n {
            let mut lhs = BigInt::zero();
            if fits_i64 {
                for &(j, a_ij) in &a_i64[i] {
                    if !x_int[j].is_zero() {
                        lhs += &x_int[j] * BigInt::from(a_ij);
                    }
                }
            } else {
                for (a_ij, xj) in a_int[i].iter().zip(&x_int) {
                    if !a_ij.is_zero() && !xj.is_zero() {
                        lhs += a_ij * xj;
                    }
                }
            }
            let rhs = if b_int[i][col_idx].is_zero() {
                BigInt::zero()
            } else if col_lcm.is_one() {
                b_int[i][col_idx].clone()
            } else {
                &b_int[i][col_idx] * &col_lcm
            };
            if lhs != rhs {
                return None;
            }
        }
        x_cols.push(q_col);
    }

    // Transpose x_cols (rhs_cols x n) into row-major QMat (n x rhs_cols)
    let mut out = vec![vec![Q::zero(); rhs_cols]; n];
    for (col_idx, col) in x_cols.into_iter().enumerate() {
        for (row_idx, val) in col.into_iter().enumerate() {
            out[row_idx][col_idx] = val;
        }
    }
    Some(out)
}

/// Solves a matrix linear system $A X = B$ where $A$ has `cols` columns and `B`
/// has matching row count. Free variables in underdetermined systems are set to zero.
///
/// Uses Dixon $p$-adic lifting for square nonsingular systems of order $n \ge 6$
/// and sparse-pivot elimination otherwise.
#[must_use]
pub fn solve_mat(
    m: &[Vec<Q>],
    cols: usize,
    b: &[Vec<Q>],
) -> Option<QMat> {
    let rows = m.len();
    if rows != b.len() || m.iter().any(|row| row.len() != cols) {
        return None;
    }
    let rhs_cols = b.first().map_or(0, Vec::len);
    if b.iter().any(|row| row.len() != rhs_cols) {
        return None;
    }
    Layout::array::<Q>(cols.checked_add(rhs_cols)?).ok()?;
    if rows == 0 || b.iter().all(|r| r.iter().all(Zero::is_zero)) {
        return Some(vec![zeros(rhs_cols); cols]);
    }
    if rows == cols
        && rows >= 6
        && let Some(sol) = solve_dixon_dense_multi(m, b)
    {
        return Some(sol);
    }

    let total_cols = cols + rhs_cols;
    let mut aug: QMat = m
        .iter()
        .zip(b)
        .map(|(row, brow)| row.iter().chain(brow.iter()).cloned().collect())
        .collect();

    let mut pivots = Vec::with_capacity(cols.min(aug.len()));
    let mut sparse_pivot: Vec<(usize, Q)> = Vec::new();
    for c in 0..cols {
        let rank = pivots.len();
        if rank == aug.len() {
            break;
        }
        let Some(p) = (rank..aug.len()).find(|&i| !aug[i][c].is_zero()) else {
            continue;
        };
        aug.swap(rank, p);
        let (head, tail) = aug.split_at_mut(rank + 1);
        let pivot = &head[rank];
        let pivot_c = &pivot[c];
        sparse_pivot.clear();
        for (j, val) in pivot[c + 1..total_cols].iter().enumerate() {
            if !val.is_zero() {
                sparse_pivot.push((c + 1 + j, val.clone()));
            }
        }
        let grain = adaptive_row_grain(tail.len(), sparse_pivot.len());
        let pivot_pairs = sparse_pivot.as_slice();
        super::par::for_each_mut(tail, grain, |_, row| {
            if row[c].is_zero() {
                return;
            }
            let factor = std::mem::take(&mut row[c]) / pivot_c;
            for &(col_idx, ref y) in pivot_pairs {
                row[col_idx] -= &factor * y;
            }
        });
        pivots.push(c);
    }

    if aug[pivots.len()..]
        .iter()
        .any(|row| row[cols..total_cols].iter().any(|v| !v.is_zero()))
    {
        return None;
    }

    let mut x = vec![zeros(rhs_cols); cols];
    for (row, p) in aug.into_iter().zip(pivots).rev() {
        let pivot_val = &row[p];
        let mut rhs_row = row[cols..total_cols].to_vec();
        for k in p + 1..cols {
            let coeff = &row[k];
            if !coeff.is_zero() {
                for (r_val, x_val) in rhs_row.iter_mut().zip(&x[k]) {
                    if !x_val.is_zero() {
                        *r_val -= coeff * x_val;
                    }
                }
            }
        }
        for r_val in &mut rhs_row {
            if !r_val.is_zero() {
                *r_val = &*r_val / pivot_val;
            }
        }
        x[p] = rhs_row;
    }
    Some(x)
}

/// A particular solution of `m x = b`, with free variables set to zero.
///
/// Returns `None` for an inconsistent system, a mismatched right-hand
/// side, or an invalid matrix shape. Pivots are chosen from left to right.
#[must_use]
pub fn solve(
    m: &[Vec<Q>],
    cols: usize,
    b: &[Q],
) -> Option<Vec<Q>> {
    solve_impl(m, cols, b, false)
}

/// Solves `m * x = b` over `Q` when the system has a unique solution (full
/// column rank `cols`), returning `None` if inconsistent or underdetermined.
#[must_use]
pub fn solve_unique(
    m: &[Vec<Q>],
    cols: usize,
    b: &[Q],
) -> Option<Vec<Q>> {
    solve_impl(m, cols, b, true)
}

fn solve_impl(
    m: &[Vec<Q>],
    cols: usize,
    b: &[Q],
    require_unique: bool,
) -> Option<Vec<Q>> {
    if m.len() != b.len() || m.iter().any(|row| row.len() != cols) {
        return None;
    }
    Layout::array::<Q>(cols.checked_add(1)?).ok()?;
    if !require_unique && b.iter().all(Zero::is_zero) {
        return Some(zeros(cols));
    }
    // Fast Dixon p-adic lifting path for square full-rank systems of dimension >= 8
    if m.len() == cols
        && cols >= 8
        && let Some(sol) = solve_dixon(m, b)
    {
        return Some(sol);
    }
    let mut aug: QMat = m
        .iter()
        .zip(b)
        .map(|(row, bi)| {
            row.iter()
                .cloned()
                .chain(std::iter::once(bi.clone()))
                .collect()
        })
        .collect();

    // Only rows below the pivot can affect later pivot choices. Keep already
    // selected rows in echelon form for back substitution.
    let mut pivots = Vec::with_capacity(cols.min(aug.len()));
    let mut sparse_pivot: Vec<(usize, Q)> = Vec::new();
    for c in 0..cols {
        let rank = pivots.len();
        if rank == aug.len() {
            break;
        }
        let Some(p) = (rank..aug.len()).find(|&i| !aug[i][c].is_zero()) else {
            continue;
        };
        aug.swap(rank, p);
        let (head, tail) = aug.split_at_mut(rank + 1);
        let pivot = &head[rank];
        let pivot_c = &pivot[c];
        sparse_pivot.clear();
        for (j, y) in pivot[c + 1..].iter().enumerate() {
            if !y.is_zero() {
                sparse_pivot.push((c + 1 + j, y.clone()));
            }
        }
        let grain = adaptive_row_grain(tail.len(), sparse_pivot.len());
        let pivot_pairs = sparse_pivot.as_slice();
        super::par::for_each_mut(tail, grain, |_, row| {
            if row[c].is_zero() {
                return;
            }
            let factor = std::mem::take(&mut row[c]) / pivot_c;
            for &(col_idx, ref y) in pivot_pairs {
                row[col_idx] -= &factor * y;
            }
        });
        pivots.push(c);
    }
    if require_unique && pivots.len() != cols {
        return None;
    }
    // A dependent row has a zero coefficient block and must have zero RHS.
    if aug[pivots.len()..].iter().any(|row| !row[cols].is_zero()) {
        return None;
    }

    let mut x = zeros(cols);
    for (mut row, p) in aug.into_iter().zip(pivots).rev() {
        let mut residual = row.pop()?;
        for (coefficient, value) in row[p + 1..].iter().zip(&x[p + 1..]) {
            if !coefficient.is_zero() && !value.is_zero() {
                residual -= coefficient * value;
            }
        }
        x[p] = residual / &row[p];
    }
    Some(x)
}

/// The product of two matrices (`a` is `r x m`, `b` is `m x c`).
#[must_use]
pub fn matmul(
    a: &[Vec<Q>],
    b: &[Vec<Q>],
) -> QMat {
    let rows = a.len();
    let cols = b.first().map_or(0, Vec::len);
    if rows == 0 || cols == 0 || b.is_empty() {
        return vec![zeros(cols); rows];
    }
    let grain = adaptive_row_grain(rows, cols.saturating_mul(b.len()));
    super::par::map_slice(a, grain, |row| {
        let mut out_row = zeros(cols);
        for (x, brow) in row.iter().zip(b) {
            if x.is_zero() {
                continue;
            }
            if x.is_one() {
                for (dst, y) in out_row.iter_mut().zip(brow) {
                    if !y.is_zero() {
                        *dst += y;
                    }
                }
            } else {
                for (dst, y) in out_row.iter_mut().zip(brow) {
                    if !y.is_zero() {
                        *dst += x * y;
                    }
                }
            }
        }
        out_row
    })
}

/// Alias for [`matmul`]: the product of two exact rational matrices.
#[must_use]
pub fn mul(
    a: &[Vec<Q>],
    b: &[Vec<Q>],
) -> QMat {
    matmul(a, b)
}

/// The product of a matrix and a vector.
#[must_use]
pub fn matvec(
    a: &[Vec<Q>],
    v: &[Q],
) -> Vec<Q> {
    let grain = adaptive_row_grain(a.len(), v.len());
    super::par::map_slice(a, grain, |row| {
        let mut s = Q::zero();
        for (x, y) in row.iter().zip(v) {
            if !x.is_zero() && !y.is_zero() {
                s += x * y;
            }
        }
        s
    })
}

/// The transpose of an `r x c` matrix.
#[must_use]
pub fn transpose(
    a: &[Vec<Q>],
    cols: usize,
) -> QMat {
    (0..cols)
        .map(|j| a.iter().map(|row| row[j].clone()).collect())
        .collect()
}

/// The inverse of a square matrix.
#[must_use]
pub fn inverse(m: &[Vec<Q>]) -> Option<QMat> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }
    if n >= 8
        && let Some(inv) = solve_dixon_dense_multi(m, &identity(n))
    {
        return Some(inv);
    }
    let mut aug: QMat = m
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .cloned()
                .chain((0..n).map(|j| q(i64::from(i == j))))
                .collect()
        })
        .collect();
    let pivots = rref(&mut aug, 2 * n);
    if pivots.len() != n || pivots.last().is_some_and(|&p| p >= n) {
        return None;
    }
    Some(aug.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// Alias for [`inverse`]: the inverse of a square matrix.
#[must_use]
pub fn inv(m: &[Vec<Q>]) -> Option<QMat> {
    inverse(m)
}

// ---------------------------------------------------------------------------
// Exact Determinants: Closed-form (n <= 4), Bareiss (n <= 14), Modular CRT (n > 14)
// ---------------------------------------------------------------------------

/// Closed-form integer determinant for $n \le 4$.
#[must_use]
fn det_small_int(a: &[Vec<BigInt>]) -> BigInt {
    match a.len() {
        | 0 => BigInt::one(),
        | 1 => a[0][0].clone(),
        | 2 => &a[0][0] * &a[1][1] - &a[0][1] * &a[1][0],
        | 3 => {
            let m01 = &a[1][1] * &a[2][2] - &a[1][2] * &a[2][1];
            let m02 = &a[1][0] * &a[2][2] - &a[1][2] * &a[2][0];
            let m03 = &a[1][0] * &a[2][1] - &a[1][1] * &a[2][0];
            &a[0][0] * m01 - &a[0][1] * m02 + &a[0][2] * m03
        },
        | 4 => {
            // 2x2 minors of rows 2 and 3
            let s01 = &a[2][0] * &a[3][1] - &a[2][1] * &a[3][0];
            let s02 = &a[2][0] * &a[3][2] - &a[2][2] * &a[3][0];
            let s03 = &a[2][0] * &a[3][3] - &a[2][3] * &a[3][0];
            let s12 = &a[2][1] * &a[3][2] - &a[2][2] * &a[3][1];
            let s13 = &a[2][1] * &a[3][3] - &a[2][3] * &a[3][1];
            let s23 = &a[2][2] * &a[3][3] - &a[2][3] * &a[3][2];

            let c0 = &a[1][1] * &s23 - &a[1][2] * &s13 + &a[1][3] * &s12;
            let c1 = &a[1][0] * &s23 - &a[1][2] * &s03 + &a[1][3] * &s02;
            let c2 = &a[1][0] * &s13 - &a[1][1] * &s03 + &a[1][3] * &s01;
            let c3 = &a[1][0] * &s12 - &a[1][1] * &s02 + &a[1][2] * &s01;

            &a[0][0] * c0 - &a[0][1] * c1 + &a[0][2] * c2 - &a[0][3] * c3
        },
        | _ => unreachable!(),
    }
}

/// Fraction-free Bareiss determinant of a square integer matrix over $\mathbb{Z}$.
///
/// Uses exact integer division by the previous step's pivot (`(a[i][j] * pivot - a[i][k] * a[k][j]) / prev_pivot`)
/// without any rational GCD reductions.
#[must_use]
pub fn det_bareiss(m: &[Vec<BigInt>]) -> BigInt {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return BigInt::zero();
    }
    if n <= 4 {
        return det_small_int(m);
    }
    let mut a = m.to_vec();
    let mut sign = 1_i32;
    let mut prev_pivot = BigInt::one();

    for k in 0..n - 1 {
        let Some(p_row) = (k..n).find(|&i| !a[i][k].is_zero()) else {
            return BigInt::zero();
        };
        if p_row != k {
            a.swap(k, p_row);
            sign = -sign;
        }
        let pivot = a[k][k].clone();
        let (head, tail) = a.split_at_mut(k + 1);
        let pivot_row = &head[k];
        let prev = &prev_pivot;
        let grain = adaptive_row_grain(tail.len(), n - k);
        super::par::for_each_mut(tail, grain, |_, row| {
            let row_k = std::mem::replace(&mut row[k], BigInt::zero());
            if row_k.is_zero() {
                if prev.is_one() {
                    for j in k + 1..n {
                        if !row[j].is_zero() {
                            row[j] *= &pivot;
                        }
                    }
                } else {
                    for j in k + 1..n {
                        if !row[j].is_zero() {
                            row[j] = (&row[j] * &pivot) / prev;
                        }
                    }
                }
                return;
            }
            for (x, y) in row[k + 1..n].iter_mut().zip(&pivot_row[k + 1..n]) {
                let mut num = &*x * &pivot;
                if !y.is_zero() {
                    num -= &row_k * y;
                }
                *x = if prev.is_one() {
                    num
                } else {
                    num / prev
                };
            }
        });
        prev_pivot = pivot;
    }
    let last = a[n - 1][n - 1].clone();
    if sign < 0 { -last } else { last }
}

/// Computes an upper bound on $\log_2(2 |\det A| + 2)$ using Hadamard's inequality:
/// $|\det A| \le \prod_{i=1}^n \|A_{i,*}\|_2$.
#[must_use]
fn hadamard_bit_bound(m: &[Vec<BigInt>]) -> usize {
    let n = m.len();
    if n == 0 {
        return 1;
    }
    let log2_n = (usize::BITS - n.leading_zeros()) as usize;
    let mut total_bits = 2_usize;
    for row in m {
        let mut max_bits = 0_usize;
        let mut nnz = 0_usize;
        for val in row {
            if !val.is_zero() {
                nnz += 1;
                max_bits = max_bits.max(val.bits() as usize);
            }
        }
        if nnz == 0 {
            return 0;
        }
        let log2_nnz = (usize::BITS - nnz.leading_zeros()) as usize;
        // log2(||row||_2) <= max_bits + ceil(log2(nnz) / 2)
        total_bits = total_bits.saturating_add(max_bits + log2_nnz.div_ceil(2).min(log2_n));
    }
    total_bits + 4
}

/// Multi-prime modular determinant of a square integer matrix over $\mathbb{Z}$
/// using 62-bit Montgomery `GF(p)` Gaussian elimination and Garner CRT reconstruction.
#[must_use]
pub fn det_modular_int(m: &[Vec<BigInt>]) -> BigInt {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return BigInt::zero();
    }
    if n <= 4 {
        return det_small_int(m);
    }
    let bits_needed = hadamard_bit_bound(m);
    if bits_needed == 0 {
        return BigInt::zero();
    }
    let num_primes = bits_needed.div_ceil(61).max(1);
    // If Hadamard bound requires disproportionately many primes compared to n, Bareiss is faster
    if num_primes > n.saturating_mul(6).max(64) {
        return det_bareiss(m);
    }
    let primes = select_primes_62(num_primes);

    let eval_prime = |p: u64| -> u64 {
        let mont = Montgomery::new_prime(p);
        let mut buf = crate::kernels::reuse::global_pool().scratch_i64(n * n);
        buf.resize(n * n, 0);
        for (i, row) in m.iter().enumerate() {
            let base = i * n;
            for (j, val) in row.iter().enumerate() {
                buf[base + j] = mont.to_mont(bigint_mod_u64(val, p)).cast_signed();
            }
        }
        det_mod_montgomery(&mut buf, n, mont)
    };

    // Per-prime work is O(n^3) u64 operations
    let grain = if crate::execution::should_split_work(num_primes, n * n * n) {
        crate::execution::choose_chunk_size(num_primes, 1)
    } else {
        usize::MAX
    };
    let residues = super::par::map_indexed(num_primes, grain, |idx| eval_prime(primes[idx]));
    if residues.iter().all(|&r| r == 0) {
        return BigInt::zero();
    }
    garner_crt_signed(&primes, &residues)
}

/// The determinant of a square matrix.
#[must_use]
pub fn det(m: &[Vec<Q>]) -> Q {
    let n = m.len();
    if n == 0 {
        return Q::one();
    }
    if m.iter().any(|r| r.len() != n) || m.iter().any(|r| r.iter().all(Zero::is_zero)) {
        return Q::zero();
    }
    let (int_rows, row_scales) = clear_rows_to_primitive_ints(m);
    let mut scale_num = BigInt::one();
    let mut scale_den = BigInt::one();
    for s in row_scales {
        if s.is_zero() {
            return Q::zero();
        }
        scale_num *= s.numer();
        scale_den *= s.denom();
    }
    let int_det = if n <= 4 {
        det_small_int(&int_rows)
    } else if n <= 14 {
        det_bareiss(&int_rows)
    } else {
        det_modular_int(&int_rows)
    };
    if int_det.is_zero() {
        return Q::zero();
    }
    Q::new(int_det * scale_num, scale_den)
}

/// Determinant of a square matrix whose rows can be overwritten.
pub(crate) fn det_in_place(a: &mut [Vec<Q>]) -> Q {
    det(a)
}

// ---------------------------------------------------------------------------
// Characteristic Polynomial (Faddeev-LeVerrier & Modular Hessenberg CRT)
// ---------------------------------------------------------------------------

/// Faddeev-LeVerrier over $\mathbb{Z}$ for `a_int`, returning `[a_0, ..., a_{n-1}]`
/// where $\det(\mu I - A_{\text{int}}) = \mu^n + \sum_{k=0}^{n-1} a_k \mu^k$.
#[must_use]
fn charpoly_faddeev_int(a_int: &[Vec<BigInt>]) -> Vec<BigInt> {
    let n = a_int.len();
    let mut a_coeffs = vec![BigInt::zero(); n];
    let mut b_mat = vec![vec![BigInt::zero(); n]; n];
    let grain = adaptive_row_grain(n, n * n);

    for k in 1..=n {
        if k == 1 {
            for (i, row) in b_mat.iter_mut().enumerate() {
                row[i] = BigInt::one();
            }
        } else {
            let prev_b = &b_mat;
            let c_prev = &a_coeffs[n - k + 1];
            b_mat = super::par::map_indexed(n, grain, |i| {
                let mut row = vec![BigInt::zero(); n];
                let a_row = &a_int[i];
                for (l, a_il) in a_row.iter().enumerate() {
                    if !a_il.is_zero() {
                        for (j, b_lj) in prev_b[l].iter().enumerate() {
                            if !b_lj.is_zero() {
                                row[j] += a_il * b_lj;
                            }
                        }
                    }
                }
                row[i] += c_prev;
                row
            });
        }
        let mut tr = BigInt::zero();
        for i in 0..n {
            for l in 0..n {
                if !a_int[i][l].is_zero() && !b_mat[l][i].is_zero() {
                    tr += &a_int[i][l] * &b_mat[l][i];
                }
            }
        }
        a_coeffs[n - k] = -tr / BigInt::from(k);
    }
    a_coeffs
}

/// Computes an upper bound on $\log_2(2 |a_k| + 2)$ for every coefficient $a_k$ of
/// $\det(\mu I - A_{\text{int}})$ using Hadamard bounds on principal minors.
#[must_use]
fn charpoly_hadamard_bit_bound(a_int: &[Vec<BigInt>]) -> usize {
    let n = a_int.len();
    let log2_n = (usize::BITS - n.leading_zeros()) as usize;
    let mut total_bits = n.saturating_add(6);
    for row in a_int {
        let mut max_bits = 0_usize;
        let mut nnz = 0_usize;
        for val in row {
            if !val.is_zero() {
                nnz += 1;
                max_bits = max_bits.max(val.bits() as usize);
            }
        }
        if nnz > 0 {
            let log2_nnz = (usize::BITS - nnz.leading_zeros()) as usize;
            total_bits = total_bits.saturating_add(max_bits + log2_nnz.div_ceil(2).min(log2_n));
        }
    }
    total_bits
}

/// Characteristic polynomial $\det(\lambda I - A)$ of a square matrix $A \in \mathbb{Q}^{n \times n}$,
/// returned as coefficients $[c_0, c_1, \dots, c_n]$ from constant term $c_0$ to leading
/// coefficient $c_n = 1$ ($\sum_{k=0}^n c_k \lambda^k$).
#[must_use]
pub fn charpoly(m: &[Vec<Q>]) -> Option<Vec<Q>> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(vec![Q::one()]);
    }
    if n == 1 {
        return Some(vec![-m[0][0].clone(), Q::one()]);
    }
    if n == 2 {
        let tr = &m[0][0] + &m[1][1];
        let d = &m[0][0] * &m[1][1] - &m[0][1] * &m[1][0];
        return Some(vec![d, -tr, Q::one()]);
    }

    // Clear the whole matrix A to (scale_num / scale_den) * A_int so coefficients have minimum bit-length.
    let mut common_den = BigInt::one();
    for row in m {
        for x in row {
            if !x.is_zero() && !x.denom().is_one() {
                common_den = common_den.lcm(x.denom());
            }
        }
    }
    let mut a_int: Vec<Vec<BigInt>> = m
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| {
                    if x.is_zero() {
                        BigInt::zero()
                    } else if common_den.is_one() {
                        x.numer().clone()
                    } else {
                        x.numer() * (&common_den / x.denom())
                    }
                })
                .collect()
        })
        .collect();

    let mut content = BigInt::zero();
    'gcd_scan: for row in &a_int {
        for val in row {
            if !val.is_zero() {
                content = if content.is_zero() {
                    val.abs()
                } else {
                    content.gcd(val)
                };
                if content.is_one() {
                    break 'gcd_scan;
                }
            }
        }
    }
    if content.is_zero() {
        let mut coeffs = vec![Q::zero(); n + 1];
        coeffs[n] = Q::one();
        return Some(coeffs);
    }
    let (scale_num, scale_den) = if content > BigInt::one() {
        for row in &mut a_int {
            for val in row {
                if !val.is_zero() {
                    *val /= &content;
                }
            }
        }
        let g = content.gcd(&common_den);
        (content / &g, common_den / g)
    } else {
        (BigInt::one(), common_den)
    };

    // For n <= 10 use integer Faddeev-LeVerrier; for larger matrices use multi-prime
    // modular upper-Hessenberg reduction + Garner CRT when prime count is bounded.
    let bits_needed = charpoly_hadamard_bit_bound(&a_int);
    let num_primes = bits_needed.div_ceil(61).max(1);
    let a_coeffs = if n > 10 && num_primes <= n.saturating_mul(8).max(64) {
        let primes = select_primes_62(num_primes);
        let grain = if crate::execution::should_split_work(num_primes, n * n * n) {
            crate::execution::choose_chunk_size(num_primes, 1)
        } else {
            usize::MAX
        };
        let a_ref = &a_int;
        let mod_polys = super::par::map_indexed(num_primes, grain, |idx| {
            let p = primes[idx];
            let mont = Montgomery::new_prime(p);
            let mut buf = crate::kernels::reuse::global_pool().scratch_i64(n * n);
            buf.resize(n * n, 0);
            for (i, row) in a_ref.iter().enumerate() {
                let base = i * n;
                for (j, val) in row.iter().enumerate() {
                    buf[base + j] = mont.to_mont(bigint_mod_u64(val, p)).cast_signed();
                }
            }
            charpoly_hessenberg_mod(&mut buf, n, mont)
        });
        let mut out = Vec::with_capacity(n);
        let mut residues = Vec::with_capacity(num_primes);
        for coeff_idx in 0..n {
            residues.clear();
            for poly_mod in &mod_polys {
                residues.push(poly_mod[coeff_idx]);
            }
            out.push(garner_crt_signed(&primes, &residues));
        }
        out
    } else {
        charpoly_faddeev_int(&a_int)
    };

    let mut coeffs = vec![Q::zero(); n + 1];
    coeffs[n] = Q::one();
    let mut num_pow = BigInt::one();
    let mut den_pow = BigInt::one();
    for k in 1..=n {
        num_pow *= &scale_num;
        den_pow *= &scale_den;
        if !a_coeffs[n - k].is_zero() {
            coeffs[n - k] = Q::new(&a_coeffs[n - k] * &num_pow, den_pow.clone());
        }
    }
    Some(coeffs)
}

// ---------------------------------------------------------------------------
// Integer Normal Forms: Hermite Normal Form & Smith Normal Form
// ---------------------------------------------------------------------------

#[must_use]
fn int_identity(n: usize) -> Vec<Vec<BigInt>> {
    (0..n)
        .map(|i| {
            let mut row = vec![BigInt::zero(); n];
            row[i] = BigInt::one();
            row
        })
        .collect()
}

fn hermite_normal_form_impl(
    m: &[Vec<BigInt>],
    track_u: bool,
) -> (Vec<Vec<BigInt>>, Option<Vec<Vec<BigInt>>>) {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    if m.iter().any(|r| r.len() != cols) {
        return (Vec::new(), track_u.then(|| int_identity(rows)));
    }
    if rows == 0 || cols == 0 {
        return (m.to_vec(), track_u.then(|| int_identity(rows)));
    }
    let mut h = m.to_vec();
    let mut u = track_u.then(|| int_identity(rows));
    let mut lead = 0_usize;
    let mut pivots = Vec::with_capacity(rows.min(cols));

    for c in 0..cols {
        if lead >= rows {
            break;
        }
        loop {
            let mut best_row: Option<usize> = None;
            for r in lead..rows {
                if !h[r][c].is_zero() && best_row.is_none_or(|b| h[r][c].abs() < h[b][c].abs()) {
                    best_row = Some(r);
                }
            }
            let Some(smallest) = best_row else {
                break;
            };
            if smallest != lead {
                h.swap(lead, smallest);
                if let Some(u_mat) = &mut u {
                    u_mat.swap(lead, smallest);
                }
            }
            let pivot = h[lead][c].clone();
            let mut done = true;
            for r in lead + 1..rows {
                if !h[r][c].is_zero() {
                    let q_div = h[r][c].div_floor(&pivot);
                    if !q_div.is_zero() {
                        let (head, tail) = h.split_at_mut(r);
                        let pivot_row = &head[lead];
                        let target_row = &mut tail[0];
                        for (val, s) in target_row[c..cols].iter_mut().zip(&pivot_row[c..cols]) {
                            if !s.is_zero() {
                                *val -= &q_div * s;
                            }
                        }
                        if let Some(u_mat) = &mut u {
                            let (u_head, u_tail) = u_mat.split_at_mut(r);
                            for (val, s) in u_tail[0].iter_mut().zip(&u_head[lead]) {
                                if !s.is_zero() {
                                    *val -= &q_div * s;
                                }
                            }
                        }
                    }
                    if !h[r][c].is_zero() {
                        done = false;
                    }
                }
            }
            if done {
                break;
            }
        }
        if h[lead][c].is_zero() {
            continue;
        }
        if h[lead][c].is_negative() {
            for val in &mut h[lead][c..cols] {
                *val = -std::mem::take(val);
            }
            if let Some(u_mat) = &mut u {
                for val in &mut u_mat[lead] {
                    *val = -std::mem::take(val);
                }
            }
        }
        pivots.push((lead, c));
        lead += 1;
    }

    for &(r, c) in &pivots {
        let pivot = h[r][c].clone();
        // Rows `0..r` are the ones being reduced and row `r` is the pivot row
        // they are reduced against; splitting the slice keeps both borrows
        // disjoint instead of holding an immutable borrow of the whole matrix.
        let (upper, lower) = h.split_at_mut(r);
        let pivot_row = &lower[0];
        for (above, target) in upper.iter_mut().enumerate() {
            let q_div = target[c].div_floor(&pivot);
            if !q_div.is_zero() {
                for j in c..cols {
                    let s = &pivot_row[j];
                    if !s.is_zero() {
                        target[j] -= &q_div * s;
                    }
                }
                if let Some(u_mat) = &mut u {
                    let (u_upper, u_lower) = u_mat.split_at_mut(r);
                    let u_pivot = &u_lower[0];
                    if let Some(u_target) = u_upper.get_mut(above) {
                        for j in 0..rows {
                            let s = &u_pivot[j];
                            if !s.is_zero() {
                                u_target[j] -= &q_div * s;
                            }
                        }
                    }
                }
            }
        }
    }
    (h, u)
}

/// Row-style Hermite Normal Form (HNF) of an integer matrix over $\mathbb{Z}$.
///
/// Brings `m` by unimodular row operations into upper-echelon Hermite normal
/// form $H$ where every leading non-zero pivot $H_{r, c} > 0$ and all entries
/// above a pivot satisfy $0 \le H_{i, c} < H_{r, c}$.
#[must_use]
pub fn hermite_normal_form(m: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    hermite_normal_form_impl(m, false).0
}

/// Row-style Hermite Normal Form $(H, U)$ of an integer matrix over $\mathbb{Z}$
/// satisfying $U A = H$ with $U$ unimodular ($\det U = \pm 1$).
#[must_use]
pub fn hermite_normal_form_with_transform(
    m: &[Vec<BigInt>]
) -> (Vec<Vec<BigInt>>, Vec<Vec<BigInt>>) {
    let rows = m.len();
    let (h, u) = hermite_normal_form_impl(m, true);
    (h, u.unwrap_or_else(|| int_identity(rows)))
}

fn smith_normal_form_impl(
    m: &[Vec<BigInt>],
    track_uv: bool,
) -> (
    Option<Vec<Vec<BigInt>>>,
    Vec<Vec<BigInt>>,
    Option<Vec<Vec<BigInt>>>,
) {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    if m.iter().any(|r| r.len() != cols) {
        return (
            track_uv.then(|| int_identity(rows)),
            Vec::new(),
            track_uv.then(|| int_identity(cols)),
        );
    }
    if rows == 0 || cols == 0 {
        return (
            track_uv.then(|| int_identity(rows)),
            m.to_vec(),
            track_uv.then(|| int_identity(cols)),
        );
    }
    let mut d = m.to_vec();
    let mut u = track_uv.then(|| int_identity(rows));
    let mut v = track_uv.then(|| int_identity(cols));
    let limit = rows.min(cols);

    for t in 0..limit {
        loop {
            let mut best: Option<(usize, usize)> = None;
            for (i, row) in d.iter().enumerate().skip(t) {
                for (j, val) in row.iter().enumerate().skip(t) {
                    if !val.is_zero() && best.is_none_or(|(bi, bj)| val.abs() < d[bi][bj].abs()) {
                        best = Some((i, j));
                    }
                }
            }
            let Some((pi, pj)) = best else {
                break;
            };
            if pi != t {
                d.swap(t, pi);
                if let Some(u_mat) = &mut u {
                    u_mat.swap(t, pi);
                }
            }
            if pj != t {
                for row in &mut d {
                    row.swap(t, pj);
                }
                if let Some(v_mat) = &mut v {
                    for row in v_mat.iter_mut() {
                        row.swap(t, pj);
                    }
                }
            }
            let pivot = d[t][t].clone();
            let mut clean = true;
            for i in t + 1..rows {
                if !d[i][t].is_zero() {
                    let q_div = d[i][t].div_floor(&pivot);
                    if !q_div.is_zero() {
                        let (head, tail) = d.split_at_mut(i);
                        let pivot_row = &head[t];
                        let target_row = &mut tail[0];
                        for (val, s) in target_row[t..cols].iter_mut().zip(&pivot_row[t..cols]) {
                            if !s.is_zero() {
                                *val -= &q_div * s;
                            }
                        }
                        if let Some(u_mat) = &mut u {
                            let (u_head, u_tail) = u_mat.split_at_mut(i);
                            for (val, s) in u_tail[0].iter_mut().zip(&u_head[t]) {
                                if !s.is_zero() {
                                    *val -= &q_div * s;
                                }
                            }
                        }
                    }
                    if !d[i][t].is_zero() {
                        clean = false;
                    }
                }
            }
            for j in t + 1..cols {
                if !d[t][j].is_zero() {
                    let q_div = d[t][j].div_floor(&pivot);
                    if !q_div.is_zero() {
                        for row in d.iter_mut().skip(t) {
                            let s = row[t].clone();
                            if !s.is_zero() {
                                row[j] -= &q_div * s;
                            }
                        }
                        if let Some(v_mat) = &mut v {
                            for row in v_mat.iter_mut() {
                                let s = row[t].clone();
                                if !s.is_zero() {
                                    row[j] -= &q_div * s;
                                }
                            }
                        }
                    }
                    if !d[t][j].is_zero() {
                        clean = false;
                    }
                }
            }
            if !clean {
                continue;
            }
            let mut offender = None;
            'search: for (i, row) in d.iter().enumerate().skip(t + 1) {
                for entry in row.iter().skip(t + 1) {
                    if !(entry % &pivot).is_zero() {
                        offender = Some(i);
                        break 'search;
                    }
                }
            }
            if let Some(i) = offender {
                let (head, tail) = d.split_at_mut(i);
                for (val, s) in head[t][t..cols].iter_mut().zip(&tail[0][t..cols]) {
                    if !s.is_zero() {
                        *val += s;
                    }
                }
                if let Some(u_mat) = &mut u {
                    let (u_head, u_tail) = u_mat.split_at_mut(i);
                    for (val, s) in u_head[t].iter_mut().zip(&u_tail[0]) {
                        if !s.is_zero() {
                            *val += s;
                        }
                    }
                }
                continue;
            }
            break;
        }
        if d[t][t].is_negative() {
            for val in &mut d[t][t..cols] {
                *val = -std::mem::take(val);
            }
            if let Some(u_mat) = &mut u {
                for val in &mut u_mat[t] {
                    *val = -std::mem::take(val);
                }
            }
        }
    }

    (u, d, v)
}

/// Computes the non-zero diagonal invariant factors $[d_1, d_2, \dots, d_r]$
/// of the Smith Normal Form of an integer matrix over $\mathbb{Z}$, with
/// $d_i > 0$ and $d_1 \mid d_2 \mid \dots \mid d_r$.
#[must_use]
pub fn smith_normal_form_diag(m: &[Vec<BigInt>]) -> Vec<BigInt> {
    let (_, d, _) = smith_normal_form_impl(m, false);
    let rows = d.len();
    let cols = d.first().map_or(0, Vec::len);
    let limit = rows.min(cols);
    (0..limit)
        .filter(|&i| !d[i][i].is_zero())
        .map(|i| d[i][i].clone())
        .collect()
}

/// Computes the Smith Normal Form $(U, D, V)$ of an integer matrix $A$ over $\mathbb{Z}$,
/// satisfying $U A V = D$ with $U, V$ unimodular and non-negative diagonal invariant factors
/// $d_1 \mid d_2 \mid \dots \mid d_r$ in $D$.
#[must_use]
pub fn smith_normal_form(
    m: &[Vec<BigInt>]
) -> (Vec<Vec<BigInt>>, Vec<Vec<BigInt>>, Vec<Vec<BigInt>>) {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let (u, d, v) = smith_normal_form_impl(m, true);
    (
        u.unwrap_or_else(|| int_identity(rows)),
        d,
        v.unwrap_or_else(|| int_identity(cols)),
    )
}

/// Reduces `v` modulo the span of the echelon rows `basis` (with their
/// `pivots`, as returned by [`rref`]): the remainder vanishes exactly on
/// the span.
#[must_use]
pub fn reduce_mod(
    v: &[Q],
    basis: &[Vec<Q>],
    pivots: &[usize],
) -> Vec<Q> {
    let mut v = v.to_vec();
    for (row, &p) in basis.iter().zip(pivots) {
        if v[p].is_zero() {
            continue;
        }
        let f = std::mem::replace(&mut v[p], Q::zero());
        for (x, y) in v.iter_mut().zip(row).skip(p + 1) {
            if !y.is_zero() {
                *x -= &f * y;
            }
        }
    }
    v
}

/// Whether a vector is zero.
#[must_use]
pub fn is_zero_vec(v: &[Q]) -> bool {
    v.iter().all(Zero::is_zero)
}

/// Helper to Clear Rational Matrix to Integer Matrix for Testing/Internal Use
#[cfg(test)]
fn to_int_mat(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|r| r.iter().map(|&x| BigInt::from(x)).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: &[&[i64]]) -> QMat {
        rows.iter()
            .map(|r| r.iter().map(|&x| q(x)).collect())
            .collect()
    }

    fn hilbert_q(n: usize) -> QMat {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| Q::new(BigInt::one(), BigInt::from((i + j + 1) as i64)))
                    .collect()
            })
            .collect()
    }

    fn vandermonde_int(nodes: &[i64]) -> Vec<Vec<BigInt>> {
        let n = nodes.len();
        nodes
            .iter()
            .map(|&x| {
                let mut row = Vec::with_capacity(n);
                let mut acc = BigInt::one();
                let bx = BigInt::from(x);
                for _ in 0..n {
                    row.push(acc.clone());
                    acc *= &bx;
                }
                row
            })
            .collect()
    }

    fn int_matmul(
        a: &[Vec<BigInt>],
        b: &[Vec<BigInt>],
    ) -> Vec<Vec<BigInt>> {
        let rows = a.len();
        let cols = b.first().map_or(0, Vec::len);
        let mut out = vec![vec![BigInt::zero(); cols]; rows];
        for i in 0..rows {
            for (k, aik) in a[i].iter().enumerate() {
                if !aik.is_zero() {
                    for (dst, bkj) in out[i].iter_mut().zip(&b[k]) {
                        *dst += aik * bkj;
                    }
                }
            }
        }
        out
    }

    #[test]
    fn montgomery_primes_and_inverses_are_exact() {
        for &p in &MONT_PRIMES_62 {
            assert!(p < (1_u64 << 62));
            assert!(is_prime_u64(p), "prime {p:#x} failed primality");
            let mont = Montgomery::new_prime(p);
            assert_eq!(
                p.wrapping_mul(mont.n_prime().wrapping_neg()),
                1_u64,
                "n_prime mismatch for {p:#x}"
            );
            for sample in [1_u64, 2, 7, 1_000_003, p - 1] {
                let xm = mont.to_mont(sample);
                assert_eq!(mont.from_mont(xm), sample);
                let inv = mont_inverse(mont, xm);
                assert_eq!(mont.mul_mont(xm, inv), mont.one());
                let neg_xm = mont.neg(xm);
                assert_eq!(mont.add(xm, neg_xm), 0);
            }
        }
        let extra = select_primes_62(40);
        assert_eq!(extra.len(), 40);
        for &p in &extra {
            assert!(is_prime_u64(p));
        }
    }

    #[test]
    fn linear_algebra() {
        let a = m(&[&[1, 2, 3], &[2, 4, 6], &[1, 0, 1]]);
        assert_eq!(rank(&a, 3), 2);
        assert_eq!(nullspace(&a, 3).len(), 1);
        assert_eq!(det(&m(&[&[2, 1], &[1, 1]])), q(1));
        let Some(inv_a) = inverse(&m(&[&[2, 1], &[1, 1]])) else {
            panic!("invertible");
        };
        assert_eq!(inv_a, m(&[&[1, -1], &[-1, 2]]));
        assert_eq!(inv(&m(&[&[2, 1], &[1, 1]])), Some(inv_a));
        assert!(inverse(&a).is_none());
        let Some(x) = solve(&m(&[&[1, 1], &[1, -1]]), 2, &[q(3), q(1)]) else {
            panic!("solvable");
        };
        assert_eq!(x, vec![q(2), q(1)]);
        assert!(solve(&m(&[&[1, 1], &[1, 1]]), 2, &[q(1), q(2)]).is_none());
        assert_eq!(matmul(&m(&[&[1, 2]]), &m(&[&[3], &[4]])), m(&[&[11]]));
        assert_eq!(mul(&m(&[&[1, 2]]), &m(&[&[3], &[4]])), m(&[&[11]]));
        let mut b = m(&[&[0, 2], &[1, 1]]);
        let pivots = rref(&mut b, 2);
        assert_eq!(pivots, vec![0, 1]);
        assert!(is_zero_vec(&reduce_mod(&[q(5), q(7)], &b, &pivots)));
    }

    #[test]
    fn bareiss_and_modular_crt_agree_on_vandermonde_and_hilbert() {
        // Exact Vandermonde determinant: prod_{i < j} (x_j - x_i)
        let nodes: Vec<i64> = (1..=16).map(|k| k * 2 - 17).collect();
        for n in [3_usize, 4, 7, 12, 16] {
            let sub_nodes = &nodes[..n];
            let v_int = vandermonde_int(sub_nodes);
            let mut expected = BigInt::one();
            for i in 0..n {
                for j in i + 1..n {
                    expected *= BigInt::from(sub_nodes[j] - sub_nodes[i]);
                }
            }
            let d_bareiss = det_bareiss(&v_int);
            let d_mod = det_modular_int(&v_int);
            assert_eq!(d_bareiss, expected, "Bareiss Vandermonde n={n}");
            assert_eq!(d_mod, expected, "Modular CRT Vandermonde n={n}");
            let v_q: QMat = v_int
                .iter()
                .map(|r| r.iter().cloned().map(Q::from_integer).collect())
                .collect();
            assert_eq!(
                det(&v_q),
                Q::from_integer(expected),
                "det(Q) Vandermonde n={n}"
            );
        }

        // Singular 16x16 matrix gives exact zero determinant
        let mut singular_16 = vandermonde_int(&nodes);
        singular_16[15] = singular_16[0].clone();
        assert_eq!(det_modular_int(&singular_16), BigInt::zero());
        assert_eq!(det_bareiss(&singular_16), BigInt::zero());

        // Hilbert matrix determinants and inverses for n = 4, 6, 8
        // det(H_4) = 1 / 6048000
        let h4 = hilbert_q(4);
        assert_eq!(det(&h4), Q::new(BigInt::one(), BigInt::from(6_048_000_i64)));
        let h6 = hilbert_q(6);
        let Some(h6_inv) = inverse(&h6) else {
            panic!("H_6 is invertible");
        };
        assert_eq!(matmul(&h6, &h6_inv), identity(6));

        let h8 = hilbert_q(8);
        let Some(h8_inv) = inverse(&h8) else {
            panic!("H_8 is invertible via Dixon lifting");
        };
        assert_eq!(matmul(&h8, &h8_inv), identity(8));
    }

    #[test]
    fn dixon_lifting_and_solve_mat_handle_rational_and_multi_rhs_systems() {
        for n in [4_usize, 8, 12, 18] {
            let a: QMat = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            if i == j {
                                Q::new(BigInt::from((n * 3 + i + 1) as i64), BigInt::from(2))
                            } else {
                                Q::new(
                                    BigInt::from(((i * 7 + j * 11) % 9) as i64 - 4),
                                    BigInt::from(((i + j) % 5 + 1) as i64),
                                )
                            }
                        })
                        .collect()
                })
                .collect();
            let b: Vec<Q> = (0..n)
                .map(|i| {
                    Q::new(
                        BigInt::from((i as i64) * 3 - 5),
                        BigInt::from((i % 4 + 1) as i64),
                    )
                })
                .collect();
            let Some(x_dixon) = solve_dixon(&a, &b) else {
                panic!("Dixon solve succeeds");
            };
            assert_eq!(matvec(&a, &x_dixon), b, "Dixon residual n={n}");
            let Some(x_solve) = solve(&a, n, &b) else {
                panic!("solve succeeds");
            };
            assert_eq!(x_dixon, x_solve, "Dixon matches solve n={n}");

            // Multi-RHS solve_mat
            let b_mat: QMat = (0..n)
                .map(|i| {
                    vec![
                        b[i].clone(),
                        q((i as i64) + 1),
                        Q::new(BigInt::from(-1), BigInt::from((i + 2) as i64)),
                    ]
                })
                .collect();
            let Some(x_mat) = solve_mat(&a, n, &b_mat) else {
                panic!("solve_mat succeeds");
            };
            assert_eq!(matmul(&a, &x_mat), b_mat, "solve_mat residual n={n}");
        }

        // Large 62-bit integer entries exercising i128 accumulator chunking in Dixon residual update
        let n = 18;
        let big_base = (1_i64 << 61) - 12345;
        let a_big: QMat = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let sign = if (i + j) % 2 == 0 {
                            1_i64
                        } else {
                            -1_i64
                        };
                        if i == j {
                            q(big_base + (i as i64) * 97)
                        } else {
                            q(sign * (big_base / 17 + ((i * 13 + j * 7) as i64)))
                        }
                    })
                    .collect()
            })
            .collect();
        let x_expected: Vec<Q> = (0..n)
            .map(|i| {
                Q::new(
                    BigInt::from((i as i64) - 5),
                    BigInt::from((i % 3 + 2) as i64),
                )
            })
            .collect();
        let b_big = matvec(&a_big, &x_expected);
        assert_eq!(solve_dixon(&a_big, &b_big), Some(x_expected));
    }

    #[test]
    fn hermite_and_smith_normal_forms_match_invariants() {
        let a = to_int_mat(&[&[2, 4, 4], &[-6, 6, 12], &[10, 4, 16]]);
        assert_eq!(
            smith_normal_form_diag(&a),
            vec![BigInt::from(2), BigInt::from(2), BigInt::from(156)]
        );
        let (u, d, v) = smith_normal_form(&a);
        assert_eq!(det_bareiss(&u).abs(), BigInt::one());
        assert_eq!(det_bareiss(&v).abs(), BigInt::one());
        assert_eq!(int_matmul(&int_matmul(&u, &a), &v), d);

        let b = to_int_mat(&[&[2, 3, 6, 2], &[5, 6, 1, 6], &[8, 3, 1, 1]]);
        let h = hermite_normal_form(&b);
        assert_eq!(
            h,
            to_int_mat(&[&[1, 0, 50, -11], &[0, 3, 28, -2], &[0, 0, 61, -13]])
        );
        let (h_tr, u_hnf) = hermite_normal_form_with_transform(&b);
        assert_eq!(h_tr, h);
        assert_eq!(det_bareiss(&u_hnf).abs(), BigInt::one());
        assert_eq!(int_matmul(&u_hnf, &b), h);
    }

    #[test]
    fn charpoly_and_parallel_determinism_across_worker_counts() {
        let a = m(&[&[1, 2, 0], &[0, 3, 1], &[2, -1, 4]]);
        // Characteristic polynomial of a: det(l*I - a)
        let Some(cp) = charpoly(&a) else {
            panic!("square");
        };
        assert_eq!(cp.len(), 4);
        assert_eq!(cp[3], q(1));
        assert_eq!(cp[2], q(-(1 + 3 + 4)));
        assert_eq!(cp[0], -det(&a));

        // Verify Hessenberg CRT path (n = 14) matches det(lambda * I - A) and Cayley-Hamilton
        let n_cp = 14_usize;
        let m14: QMat = (0..n_cp)
            .map(|i| {
                (0..n_cp)
                    .map(|j| {
                        Q::new(
                            BigInt::from(((i * 7 + j * 5) % 11) as i64 - 5),
                            BigInt::from(((i + j) % 3 + 1) as i64),
                        )
                    })
                    .collect()
            })
            .collect();
        let Some(cp14) = charpoly(&m14) else {
            panic!("14x14 charpoly");
        };
        assert_eq!(cp14.len(), n_cp + 1);
        assert_eq!(cp14[n_cp], q(1));
        for lambda in [q(0), q(1), Q::new(BigInt::from(3), BigInt::from(5))] {
            let mut eval_poly = Q::zero();
            for c in cp14.iter().rev() {
                eval_poly = &eval_poly * &lambda + c;
            }
            let mut shifted = m14
                .iter()
                .map(|r| r.iter().map(|x| -x).collect::<Vec<Q>>())
                .collect::<QMat>();
            for (i, row) in shifted.iter_mut().enumerate() {
                row[i] += &lambda;
            }
            assert_eq!(eval_poly, det(&shifted), "charpoly eval at {lambda}");
        }

        // 24x24 dense matrix tested across 1, 2, and 4 workers
        let n = 24;
        let big: QMat = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            q(17 + (i as i64))
                        } else {
                            Q::new(
                                BigInt::from(((i * 5 + j * 3) % 7) as i64 - 3),
                                BigInt::from(((i + j) % 3 + 1) as i64),
                            )
                        }
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<Q> = (0..n).map(|i| q((i as i64) - 4)).collect();

        let (det_1, sol_1, rref_1, cp_1) = crate::execution::with_scoped_pool(Some(1), || {
            let mut r = big.clone();
            let piv = rref(&mut r, n);
            (det(&big), solve(&big, n, &rhs), (piv, r), charpoly(&big))
        });
        for workers in [2_usize, 4] {
            let (det_w, sol_w, rref_w, cp_w) =
                crate::execution::with_scoped_pool(Some(workers), || {
                    let mut r = big.clone();
                    let piv = rref(&mut r, n);
                    (det(&big), solve(&big, n, &rhs), (piv, r), charpoly(&big))
                });
            assert_eq!(det_1, det_w, "det worker determinism ({workers})");
            assert_eq!(sol_1, sol_w, "solve worker determinism ({workers})");
            assert_eq!(rref_1, rref_w, "rref worker determinism ({workers})");
            assert_eq!(cp_1, cp_w, "charpoly worker determinism ({workers})");
        }
    }
}
