//! Factorisation and greatest common divisors of multivariate and univariate
//! polynomials over `Q` and `Z`, backed by:
//!
//! * **Fast modular polynomial arithmetic & Half-GCD over $\mathbb{F}_p\[x\]$**
//!   (`hgcd_fp`, `gcd_fp`, `ext_gcd_fp`), achieving $O(M(n)\log n)$ greatest
//!   common divisors and Bézout coefficients for machine-word primes $p < 2^{31}$
//!   with scratch arena reuse from [`crate::kernels::reuse::global_pool`].
//! * **Brown's Modular GCD + Chinese Remainder Theorem over $\mathbb{Z}\[x\]$ &
//!   $\mathbb{Q}\[x\]$** (`gcd_z`, `gcd_q`), eliminating rational PRS coefficient
//!   swell via leading-coefficient normalization, unlucky-prime detection,
//!   incremental Garner CRT reconstruction, Landau–Mignotte stability checks,
//!   and parallel multi-prime evaluation.
//! * **Zippel Evaluation–Interpolation Multivariate Modular GCD** (`poly_gcd`),
//!   evaluating non-main variables at distinct prime/integer points over
//!   $\mathbb{F}_p$, interpolating univariate GCD fibres via Newton/Vandermonde
//!   interpolation, and reconstructing over $\mathbb{Z}$ via CRT with exact
//!   trial-division verification and subresultant PRS fallback.
//! * **Kronecker & Pruned Multivariate Factorisation** (`factor`), augmented
//!   with degree-vector bounds, leading/trailing coefficient divisibility filters,
//!   multi-point evaluation pruning at $0, \pm 1, 2$, and parallel candidate
//!   batch testing via [`crate::kernels::par`].

use std::collections::HashSet;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::repr::Mono;
use super::repr::Poly;
use super::univariate;
use super::univariate::QPoly;
use super::univariate::ZPoly;
use crate::execution;
use crate::graph::Number;
use crate::kernels::par;
use crate::kernels::reuse::global_pool;

const CAP: usize = 200_000;
/// The largest univariate image degree attempted in Kronecker reduction.
const MAX_IMAGE_DEGREE: u64 = 3000;
/// The largest number of univariate factors recombined.
const MAX_PIECES: usize = 18;
/// The largest number of candidate subsets tried per split.
const MAX_SUBSETS: usize = 100_000;
/// Threshold below which Euclidean `gcd_fp` is used instead of recursive Half-GCD.
const HGCD_THRESHOLD: usize = 24;
/// Parallel batch size for Kronecker subset trial division.
const FACTOR_BATCH_SIZE: usize = 64;

/// Large 31-bit primes `< 2^31` for collision-free `u128` accumulation.
pub(crate) const WORD_PRIMES: &[u64] = &[
    2_147_483_647,
    2_147_483_629,
    2_147_483_587,
    2_147_483_579,
    2_147_483_563,
    2_147_483_549,
    2_147_483_543,
    2_147_483_497,
    2_147_483_489,
    2_147_483_477,
    2_147_483_423,
    2_147_483_399,
    2_147_483_353,
    2_147_483_323,
    2_147_483_269,
    2_147_483_249,
    2_147_483_237,
    2_147_483_179,
    2_147_483_171,
    2_147_483_137,
    2_147_483_123,
    2_147_483_077,
    2_147_483_069,
    2_147_483_059,
    2_147_483_053,
    2_147_483_033,
    2_147_483_029,
    2_147_482_951,
    2_147_482_949,
    2_147_482_943,
    2_147_482_937,
    2_147_482_921,
    2_147_482_877,
    2_147_482_873,
    2_147_482_867,
    2_147_482_859,
    2_147_482_819,
    2_147_482_817,
    2_147_482_811,
    2_147_482_801,
    2_147_482_763,
    2_147_482_739,
    2_147_482_697,
    2_147_482_693,
    2_147_482_681,
    2_147_482_663,
    2_147_482_661,
    2_147_482_621,
    2_147_482_591,
    2_147_482_583,
    2_147_482_577,
    2_147_482_507,
    2_147_482_501,
    2_147_482_481,
    2_147_482_417,
    2_147_482_409,
    2_147_482_367,
    2_147_482_361,
    2_147_482_349,
    2_147_482_343,
    2_147_482_327,
    2_147_482_291,
    2_147_482_273,
    2_147_482_237,
];

/// Returns the `index`-th large prime `< 2^31`, falling back to deterministic
/// downward search past [`WORD_PRIMES`].
#[must_use]
pub(crate) fn nth_word_prime(index: usize) -> u64 {
    if let Some(&p) = WORD_PRIMES.get(index) {
        return p;
    }
    let mut candidate = WORD_PRIMES.last().copied().unwrap_or(2_147_482_237) - 2;
    let mut extra = index - WORD_PRIMES.len();
    loop {
        if is_prime_u64(candidate) {
            if extra == 0 {
                return candidate;
            }
            extra -= 1;
        }
        candidate = candidate.saturating_sub(2);
    }
}

const fn is_prime_u64(n: u64) -> bool {
    crate::kernels::number_theory::is_prime_miller_rabin(n)
}

// ---------------------------------------------------------------------------
// Fast F_p[x] arithmetic and Half-GCD (O(M(n) log n))
// ---------------------------------------------------------------------------

pub(crate) fn trim_fp(p: &mut Vec<u64>) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

#[inline]
pub(crate) const fn fp_add(
    a: u64,
    b: u64,
    p: u64,
) -> u64 {
    let s = a + b;
    if s >= p { s - p } else { s }
}

#[inline]
pub(crate) const fn fp_sub(
    a: u64,
    b: u64,
    p: u64,
) -> u64 {
    if a >= b { a - b } else { a + p - b }
}

#[inline]
pub(crate) const fn fp_mul(
    a: u64,
    b: u64,
    p: u64,
) -> u64 {
    ((a as u128 * b as u128) % (p as u128)) as u64
}

#[inline]
pub(crate) const fn fp_pow(
    mut base: u64,
    mut exp: u64,
    p: u64,
) -> u64 {
    let mut acc = 1_u64;
    base %= p;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = fp_mul(acc, base, p);
        }
        base = fp_mul(base, base, p);
        exp >>= 1;
    }
    acc
}

/// Fast modular inverse via binary/Euclidean extended GCD over `i64`.
#[inline]
pub(crate) const fn fp_inv(
    a: u64,
    p: u64,
) -> u64 {
    let m = p as i64;
    let (mut t, mut new_t) = (0_i64, 1_i64);
    let (mut r, mut new_r) = (m, (a % p) as i64);
    while new_r != 0 {
        let q = r / new_r;
        let next_t = t - q * new_t;
        t = new_t;
        new_t = next_t;
        let next_r = r - q * new_r;
        r = new_r;
        new_r = next_r;
    }
    if r > 1 {
        return 0;
    }
    if t < 0 {
        (t + m) as u64
    } else {
        t as u64
    }
}

/// Reduces an integer polynomial modulo `p`.
#[must_use]
pub(crate) fn fp_reduce_z(
    poly: &[BigInt],
    p: u64,
) -> Vec<u64> {
    let modulus = BigInt::from(p);
    let mut out: Vec<u64> = poly
        .iter()
        .map(|c| {
            c.to_i64().map_or_else(
                || c.mod_floor(&modulus).to_u64().unwrap_or(0),
                |v| v.rem_euclid(p as i64) as u64,
            )
        })
        .collect();
    trim_fp(&mut out);
    out
}

#[must_use]
pub(crate) fn fp_poly_add(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut out = vec![0_u64; n];
    for (i, &x) in a.iter().enumerate() {
        out[i] = x;
    }
    for (i, &y) in b.iter().enumerate() {
        out[i] = fp_add(out[i], y, p);
    }
    trim_fp(&mut out);
    out
}

#[must_use]
pub(crate) fn fp_poly_sub(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut out = vec![0_u64; n];
    for (i, &x) in a.iter().enumerate() {
        out[i] = x;
    }
    for (i, &y) in b.iter().enumerate() {
        out[i] = fp_sub(out[i], y, p);
    }
    trim_fp(&mut out);
    out
}

fn fp_poly_mul_schoolbook(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    // For p < 2^31 on 64-bit platforms, each product is < 2^62; accumulating up
    // to 3 nonzero rows in `usize` before reducing `% p` stays < 3 * 2^62 + p < 2^64.
    if usize::BITS >= 64 && p < (1_u64 << 31) {
        let mut scratch = global_pool().scratch_usize(out_len);
        scratch.resize(out_len, 0);
        let p_usize = p as usize;
        let mut active_rows = 0_usize;
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            let x_usize = x as usize;
            for (j, &y) in b.iter().enumerate() {
                scratch[i + j] += x_usize * (y as usize);
            }
            active_rows += 1;
            if active_rows == 3 {
                for slot in scratch.iter_mut() {
                    *slot %= p_usize;
                }
                active_rows = 0;
            }
        }
        let mut out: Vec<u64> = scratch.iter().map(|&v| (v as u64) % p).collect();
        trim_fp(&mut out);
        return out;
    }
    let mut acc = vec![0_u128; out_len];
    let p128 = u128::from(p);
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        let x128 = u128::from(x);
        for (j, &y) in b.iter().enumerate() {
            acc[i + j] += x128 * u128::from(y);
        }
    }
    let mut out: Vec<u64> = acc.into_iter().map(|v| (v % p128) as u64).collect();
    trim_fp(&mut out);
    out
}

/// Multiplies two polynomials in $\mathbb{F}_p[x]$ using pooled scratch
/// accumulation via [`global_pool`] for small degrees and Karatsuba above
/// degree 32.
#[must_use]
pub(crate) fn fp_poly_mul(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() <= 32 || b.len() <= 32 {
        return fp_poly_mul_schoolbook(a, b, p);
    }
    // Karatsuba above degree 32
    let mut out = fp_poly_mul_karatsuba(a, b, p);
    trim_fp(&mut out);
    out
}

fn fp_poly_mul_karatsuba(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() <= 32 || b.len() <= 32 {
        return fp_poly_mul_schoolbook(a, b, p);
    }
    let m = a.len().max(b.len()) / 2;
    let (a0, a1) = a.split_at(m.min(a.len()));
    let (b0, b1) = b.split_at(m.min(b.len()));
    let z0 = fp_poly_mul_karatsuba(a0, b0, p);
    let z2 = fp_poly_mul_karatsuba(a1, b1, p);
    let a_sum = fp_poly_add(a0, a1, p);
    let b_sum = fp_poly_add(b0, b1, p);
    let z1_full = fp_poly_mul_karatsuba(&a_sum, &b_sum, p);
    let z1 = fp_poly_sub(&fp_poly_sub(&z1_full, &z0, p), &z2, p);

    let out_len = a.len() + b.len() - 1;
    let mut out = vec![0_u64; out_len];
    for (i, &c) in z0.iter().enumerate() {
        out[i] = fp_add(out[i], c, p);
    }
    for (i, &c) in z1.iter().enumerate() {
        if i + m < out_len {
            out[i + m] = fp_add(out[i + m], c, p);
        }
    }
    for (i, &c) in z2.iter().enumerate() {
        if i + 2 * m < out_len {
            out[i + 2 * m] = fp_add(out[i + 2 * m], c, p);
        }
    }
    trim_fp(&mut out);
    out
}

/// Quotient and remainder of `a` divided by `b` in $\mathbb{F}_p[x]$, leasing
/// scratch space from [`global_pool`].
#[must_use]
pub(crate) fn fp_poly_divrem(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> (Vec<u64>, Vec<u64>) {
    let mut b_len = b.len();
    while b_len > 0 && b[b_len - 1] == 0 {
        b_len -= 1;
    }
    if b_len == 0 {
        return (Vec::new(), a.to_vec());
    }
    let b = &b[..b_len];
    let mut a_len = a.len();
    while a_len > 0 && a[a_len - 1] == 0 {
        a_len -= 1;
    }
    if a_len < b_len {
        return (Vec::new(), a[..a_len].to_vec());
    }
    let mut scratch = global_pool().scratch_usize(a_len);
    scratch.clear();
    scratch.extend(a[..a_len].iter().map(|&c| c as usize));
    let inv_lead = fp_inv(b[b_len - 1], p);
    let q_len = a_len - b_len + 1;
    let mut quo = vec![0_u64; q_len];
    let p128 = u128::from(p);
    for k in (0..q_len).rev() {
        let top = scratch[k + b_len - 1] as u64;
        if top == 0 {
            continue;
        }
        let factor = fp_mul(top, inv_lead, p);
        quo[k] = factor;
        let f128 = u128::from(factor);
        for j in 0..(b_len - 1) {
            let sub = ((f128 * u128::from(b[j])) % p128) as u64;
            let cur = scratch[k + j] as u64;
            scratch[k + j] = fp_sub(cur, sub, p) as usize;
        }
        scratch[k + b_len - 1] = 0;
    }
    let mut rem: Vec<u64> = scratch[..(b_len - 1)].iter().map(|&c| c as u64).collect();
    trim_fp(&mut quo);
    trim_fp(&mut rem);
    (quo, rem)
}

#[must_use]
pub(crate) fn fp_poly_rem(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    let mut b_len = b.len();
    while b_len > 0 && b[b_len - 1] == 0 {
        b_len -= 1;
    }
    if b_len == 0 {
        return a.to_vec();
    }
    let b = &b[..b_len];
    let mut a_len = a.len();
    while a_len > 0 && a[a_len - 1] == 0 {
        a_len -= 1;
    }
    if a_len < b_len {
        return a[..a_len].to_vec();
    }
    let mut scratch = global_pool().scratch_usize(a_len);
    scratch.clear();
    scratch.extend(a[..a_len].iter().map(|&c| c as usize));
    let inv_lead = fp_inv(b[b_len - 1], p);
    let q_len = a_len - b_len + 1;
    let p128 = u128::from(p);
    for k in (0..q_len).rev() {
        let top = scratch[k + b_len - 1] as u64;
        if top == 0 {
            continue;
        }
        let factor = fp_mul(top, inv_lead, p);
        let f128 = u128::from(factor);
        for j in 0..(b_len - 1) {
            let sub = ((f128 * u128::from(b[j])) % p128) as u64;
            let cur = scratch[k + j] as u64;
            scratch[k + j] = fp_sub(cur, sub, p) as usize;
        }
        scratch[k + b_len - 1] = 0;
    }
    let mut rem: Vec<u64> = scratch[..(b_len - 1)].iter().map(|&c| c as u64).collect();
    trim_fp(&mut rem);
    rem
}

#[must_use]
pub(crate) fn fp_poly_monic(
    a: &[u64],
    p: u64,
) -> Vec<u64> {
    let Some(&lead) = a.last() else {
        return Vec::new();
    };
    if lead == 1 {
        return a.to_vec();
    }
    let inv = fp_inv(lead, p);
    a.iter().map(|&c| fp_mul(c, inv, p)).collect()
}

/// $2 \times 2$ polynomial transition matrix over $\mathbb{F}_p[x]$ with an
/// explicit determinant sign (`det = +1` when `!neg_det`, `-1` when `neg_det`).
///
/// Represents $\begin{pmatrix} m_{00} & m_{01} \\ m_{10} & m_{11} \end{pmatrix}$
/// such that $\begin{pmatrix} A \\ B \end{pmatrix} = M \begin{pmatrix} A' \\ B' \end{pmatrix}$.
#[derive(Clone, Debug)]
struct PolyMatrixFp {
    m00: Vec<u64>,
    m01: Vec<u64>,
    m10: Vec<u64>,
    m11: Vec<u64>,
    neg_det: bool,
}

impl PolyMatrixFp {
    fn identity() -> Self {
        Self {
            m00: vec![1],
            m01: Vec::new(),
            m10: Vec::new(),
            m11: vec![1],
            neg_det: false,
        }
    }

    /// Multiplies `self * rhs` in $\mathbb{F}_p[x]^{2 \times 2}$.
    fn mul(
        &self,
        rhs: &Self,
        p: u64,
    ) -> Self {
        let m00 = fp_poly_add(
            &fp_poly_mul(&self.m00, &rhs.m00, p),
            &fp_poly_mul(&self.m01, &rhs.m10, p),
            p,
        );
        let m01 = fp_poly_add(
            &fp_poly_mul(&self.m00, &rhs.m01, p),
            &fp_poly_mul(&self.m01, &rhs.m11, p),
            p,
        );
        let m10 = fp_poly_add(
            &fp_poly_mul(&self.m10, &rhs.m00, p),
            &fp_poly_mul(&self.m11, &rhs.m10, p),
            p,
        );
        let m11 = fp_poly_add(
            &fp_poly_mul(&self.m10, &rhs.m01, p),
            &fp_poly_mul(&self.m11, &rhs.m11, p),
            p,
        );
        Self {
            m00,
            m01,
            m10,
            m11,
            neg_det: self.neg_det ^ rhs.neg_det,
        }
    }

    /// Right-multiplies `self` by $\begin{pmatrix} q & 1 \\ 1 & 0 \end{pmatrix}$.
    fn step_quotient(
        &mut self,
        q: &[u64],
        p: u64,
    ) {
        let next_m00 = fp_poly_add(&fp_poly_mul(&self.m00, q, p), &self.m01, p);
        let next_m10 = fp_poly_add(&fp_poly_mul(&self.m10, q, p), &self.m11, p);
        self.m01 = std::mem::replace(&mut self.m00, next_m00);
        self.m11 = std::mem::replace(&mut self.m10, next_m10);
        self.neg_det = !self.neg_det;
    }

    /// Applies $M^{-1}$ to $(a, b)^T$.
    /// Since $\det(M) = \pm 1$:
    /// $M^{-1} = \det(M) \begin{pmatrix} m_{11} & -m_{01} \\ -m_{10} & m_{00} \end{pmatrix}$.
    fn apply_inverse(
        &self,
        a: &[u64],
        b: &[u64],
        p: u64,
    ) -> (Vec<u64>, Vec<u64>) {
        let t0 = fp_poly_sub(
            &fp_poly_mul(&self.m11, a, p),
            &fp_poly_mul(&self.m01, b, p),
            p,
        );
        let t1 = fp_poly_sub(
            &fp_poly_mul(&self.m00, b, p),
            &fp_poly_mul(&self.m10, a, p),
            p,
        );
        if self.neg_det {
            (fp_poly_sub(&[], &t0, p), fp_poly_sub(&[], &t1, p))
        } else {
            (t0, t1)
        }
    }
}

/// Recursive Half-GCD (`HGCD`) over $\mathbb{F}_p[x]$.
///
/// Given `a, b` with $\deg(a) > \deg(b)$, computes a unimodular matrix $M$
/// reducing $(a, b)$ to $(a', b')$ with $\deg(b') \le \lfloor \deg(a)/2 \rfloor < \deg(a')$.
fn hgcd_fp(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> PolyMatrixFp {
    let deg_a = a.len().saturating_sub(1);
    let m = deg_a / 2;
    if b.is_empty() || b.len() - 1 <= m {
        return PolyMatrixFp::identity();
    }
    if deg_a < HGCD_THRESHOLD {
        let mut mat = PolyMatrixFp::identity();
        let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
        while !r1.is_empty() && r1.len() - 1 > m {
            let (q, r2) = fp_poly_divrem(&r0, &r1, p);
            mat.step_quotient(&q, p);
            r0 = r1;
            r1 = r2;
        }
        return mat;
    }

    let a_hi = &a[m..];
    let b_hi = if b.len() > m { &b[m..] } else { &[] };
    let mut r_mat = hgcd_fp(a_hi, b_hi, p);
    let (mut a1, mut b1) = r_mat.apply_inverse(a, b, p);
    if b1.is_empty() || b1.len() - 1 <= m {
        return r_mat;
    }

    let (q, r) = fp_poly_divrem(&a1, &b1, p);
    r_mat.step_quotient(&q, p);
    a1 = b1;
    b1 = r;
    if b1.is_empty() || b1.len() - 1 <= m {
        return r_mat;
    }

    let deg_a1 = a1.len() - 1;
    let k = (2 * m).saturating_sub(deg_a1);
    let a1_hi = if a1.len() > k {
        &a1[k..]
    } else {
        &[]
    };
    let b1_hi = if b1.len() > k {
        &b1[k..]
    } else {
        &[]
    };
    let s_mat = hgcd_fp(a1_hi, b1_hi, p);
    r_mat.mul(&s_mat, p)
}

/// Classical Euclidean GCD over $\mathbb{F}_p[x]$ (monic).
#[must_use]
pub(crate) fn gcd_fp_euclid(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    trim_fp(&mut x);
    trim_fp(&mut y);
    while !y.is_empty() {
        let r = fp_poly_rem(&x, &y, p);
        x = y;
        y = r;
    }
    fp_poly_monic(&x, p)
}

/// Fast $O(M(n)\log n)$ Half-GCD + Euclidean polynomial GCD over $\mathbb{F}_p[x]$ (monic).
#[must_use]
pub(crate) fn gcd_fp(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> Vec<u64> {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    trim_fp(&mut x);
    trim_fp(&mut y);
    if x.len() < y.len() {
        std::mem::swap(&mut x, &mut y);
    }
    while !y.is_empty() && x.len() > HGCD_THRESHOLD {
        if x.len() - y.len() > x.len() / 2 {
            let r = fp_poly_rem(&x, &y, p);
            x = y;
            y = r;
            continue;
        }
        let mat = hgcd_fp(&x, &y, p);
        let (nx, ny) = mat.apply_inverse(&x, &y, p);
        x = nx;
        y = ny;
        if !y.is_empty() {
            let r = fp_poly_rem(&x, &y, p);
            x = y;
            y = r;
        }
    }
    if y.is_empty() {
        fp_poly_monic(&x, p)
    } else {
        gcd_fp_euclid(&x, &y, p)
    }
}

/// Computes $\mathrm{res}(a, b) \pmod p$ for univariate polynomials in $\mathbb{F}_p[x]$
/// of exact degrees `deg_a = a.len() - 1` and `deg_b = b.len() - 1`.
#[must_use]
pub(crate) fn resultant_fp(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> u64 {
    let (mut left, mut right) = (a.to_vec(), b.to_vec());
    trim_fp(&mut left);
    trim_fp(&mut right);
    if left.is_empty() || right.is_empty() {
        return 0;
    }
    let mut value = 1_u64;
    while right.len() > 1 {
        let m = left.len() - 1;
        let n = right.len() - 1;
        if m % 2 == 1 && n % 2 == 1 {
            value = fp_sub(0, value, p);
        }
        if m < n {
            std::mem::swap(&mut left, &mut right);
            continue;
        }
        let remainder = fp_poly_rem(&left, &right, p);
        if remainder.is_empty() {
            return 0;
        }
        let r = remainder.len() - 1;
        let lead = right[n];
        let factor = fp_pow(lead, (m - r) as u64, p);
        value = fp_mul(value, factor, p);
        left = right;
        right = remainder;
    }
    let m = left.len() - 1;
    let c = right[0];
    fp_mul(value, fp_pow(c, m as u64, p), p)
}

/// Hadamard bound in bits on $|\mathrm{res}(A, B)|$ for $A, B \in \mathbb{Z}[x]$:
/// $|\mathrm{res}(A, B)| \le \|A\|_2^{\deg B} \|B\|_2^{\deg A}$.
#[must_use]
pub(crate) fn hadamard_resultant_bound_bits(
    a: &[BigInt],
    b: &[BigInt],
) -> usize {
    let m = a.len().saturating_sub(1);
    let n = b.len().saturating_sub(1);
    let max_a_bits = a.iter().map(BigInt::bits).max().unwrap_or(1) as usize;
    let max_b_bits = b.iter().map(BigInt::bits).max().unwrap_or(1) as usize;
    let log_len_a = (a.len().max(1) as u64).ilog2() as usize + 2;
    let log_len_b = (b.len().max(1) as u64).ilog2() as usize + 2;
    n.saturating_mul(max_a_bits + log_len_a)
        .saturating_add(m.saturating_mul(max_b_bits + log_len_b))
        .saturating_add(16)
}

/// Extended GCD over $\mathbb{F}_p[x]$ accelerated with Half-GCD: returns
/// `(g, s, t)` with `g` monic and `s * a + t * b = g (mod p)`.
#[must_use]
pub(crate) fn ext_gcd_fp(
    a: &[u64],
    b: &[u64],
    p: u64,
) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    trim_fp(&mut r0);
    trim_fp(&mut r1);
    let mut total = PolyMatrixFp::identity();
    while !r1.is_empty() {
        if r0.len() >= r1.len() && r0.len() > HGCD_THRESHOLD && r0.len() - r1.len() <= r0.len() / 2
        {
            let m = hgcd_fp(&r0, &r1, p);
            let (nr0, nr1) = m.apply_inverse(&r0, &r1, p);
            total = total.mul(&m, p);
            r0 = nr0;
            r1 = nr1;
            if r1.is_empty() {
                break;
            }
        }
        let (q, r2) = fp_poly_divrem(&r0, &r1, p);
        total.step_quotient(&q, p);
        r0 = r1;
        r1 = r2;
    }
    // Since (r0, 0)^T = total^{-1} (a, b)^T, the top row of total^{-1} is
    // det(total) * (m11, -m01).
    let (mut s, mut t) = if total.neg_det {
        (fp_poly_sub(&[], &total.m11, p), total.m01)
    } else {
        (total.m11, fp_poly_sub(&[], &total.m01, p))
    };
    if let Some(&lead) = r0.last()
        && lead != 1
    {
        let inv = fp_inv(lead, p);
        r0 = r0.into_iter().map(|c| fp_mul(c, inv, p)).collect();
        s = s.into_iter().map(|c| fp_mul(c, inv, p)).collect();
        t = t.into_iter().map(|c| fp_mul(c, inv, p)).collect();
    }
    trim_fp(&mut s);
    trim_fp(&mut t);
    (r0, s, t)
}

// ---------------------------------------------------------------------------
// Brown's Modular GCD + Chinese Remainder Theorem over Z[x] and Q[x]
// ---------------------------------------------------------------------------

fn trim_z(p: &mut ZPoly) {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
}

fn z_content_and_primitive(p: &[BigInt]) -> (BigInt, ZPoly) {
    let mut trimmed = p.to_vec();
    trim_z(&mut trimmed);
    if trimmed.is_empty() {
        return (BigInt::zero(), Vec::new());
    }
    let mut c = trimmed.iter().fold(BigInt::zero(), |acc, x| acc.gcd(x));
    if trimmed.last().is_some_and(Signed::is_negative) {
        c = -c;
    }
    if c.is_one() {
        return (c, trimmed);
    }
    let prim = trimmed.iter().map(|x| x / &c).collect();
    (c, prim)
}

/// Centered symmetric remainder in `(-m/2, m/2]`.
pub(crate) fn symmetric_mod(
    r: &BigInt,
    m: &BigInt,
    half_m: &BigInt,
) -> BigInt {
    let pos = r.mod_floor(m);
    if &pos > half_m {
        pos - m
    } else {
        pos
    }
}

/// Incremental Chinese Remainder step: given `r1 mod m1` (with `r1 in [0, m1)`)
/// and `r2 mod p`, updates `(r1, m1)` to `r mod (m1 * p)`.
/// Returns `true` if the symmetric representative changed.
pub(crate) fn crt_step_poly(
    acc_residues: &mut Vec<BigInt>,
    modulus: &mut BigInt,
    new_residues: &[u64],
    p: u64,
) -> bool {
    let p_big = BigInt::from(p);
    if modulus.is_zero() {
        *acc_residues = new_residues.iter().map(|&c| BigInt::from(c)).collect();
        *modulus = p_big;
        return true;
    }
    let m_mod_p = modulus.mod_floor(&p_big).to_u64().unwrap_or(0);
    let m_inv = fp_inv(m_mod_p, p);
    let new_modulus = &*modulus * &p_big;
    let half_old = &*modulus >> 1;
    let half_new = &new_modulus >> 1;

    if acc_residues.len() < new_residues.len() {
        acc_residues.resize(new_residues.len(), BigInt::zero());
    }
    let mut changed = false;
    for i in 0..acc_residues.len() {
        let r1 = &acc_residues[i];
        let old_sym = symmetric_mod(r1, modulus, &half_old);
        let r2 = new_residues.get(i).copied().unwrap_or(0);
        let r1_mod_p = r1.mod_floor(&p_big).to_u64().unwrap_or(0);
        let diff = fp_sub(r2, r1_mod_p, p);
        let k = fp_mul(diff, m_inv, p);
        let updated = if k == 0 {
            r1.clone()
        } else {
            r1 + &*modulus * BigInt::from(k)
        };
        let new_sym = symmetric_mod(&updated, &new_modulus, &half_new);
        if new_sym != old_sym {
            changed = true;
        }
        acc_residues[i] = updated;
    }
    *modulus = new_modulus;
    changed
}

/// Converts positive residues modulo `modulus` into centered symmetric coefficients
/// in `(-modulus/2, modulus/2]`.
fn centered_poly(
    residues: &[BigInt],
    modulus: &BigInt,
) -> ZPoly {
    let half = modulus >> 1;
    let mut out: ZPoly = residues
        .iter()
        .map(|r| symmetric_mod(r, modulus, &half))
        .collect();
    trim_z(&mut out);
    out
}

/// Landau–Mignotte bound on the coefficients of any factor of `a` and `b`
/// scaled by `gamma = gcd(lc(a), lc(b))`.
fn mignotte_bound_bits(
    a: &[BigInt],
    b: &[BigInt],
    gamma: &BigInt,
) -> u64 {
    let norm_bits = |p: &[BigInt]| -> u64 {
        p.iter().map(BigInt::bits).max().unwrap_or(1)
            + u64::from((p.len().max(1) as u64).ilog2())
            + 2
    };
    let deg = a.len().min(b.len()).saturating_sub(1) as u64;
    let min_norm = norm_bits(a).min(norm_bits(b));
    gamma.bits() + deg + min_norm + 4
}

/// Primitive greatest common divisor of two polynomials over $\mathbb{Z}\[x\]$
/// using Brown's modular GCD algorithm with Chinese Remainder reconstruction,
/// unlucky-prime detection, and trial division.
#[must_use]
pub fn gcd_z(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    let (ca, mut pa) = z_content_and_primitive(a);
    let (cb, mut pb) = z_content_and_primitive(b);
    if pa.is_empty() {
        return pb.into_iter().map(|x| x * &cb).collect();
    }
    if pb.is_empty() {
        return pa.into_iter().map(|x| x * &ca).collect();
    }
    let c_gcd = ca.gcd(&cb);

    // Strip common trailing zeros (powers of x).
    let shift_a = pa.iter().take_while(|c| c.is_zero()).count();
    let shift_b = pb.iter().take_while(|c| c.is_zero()).count();
    let common_shift = shift_a.min(shift_b);
    if shift_a > 0 {
        pa.drain(0..shift_a);
    }
    if shift_b > 0 {
        pb.drain(0..shift_b);
    }
    let (_, pa) = z_content_and_primitive(&pa);
    let (_, pb) = z_content_and_primitive(&pb);

    let finish = |mut prim: ZPoly| -> ZPoly {
        if prim.is_empty() {
            prim.push(BigInt::one());
        }
        if common_shift > 0 {
            let mut shifted = vec![BigInt::zero(); common_shift];
            shifted.extend(prim);
            prim = shifted;
        }
        if !c_gcd.is_one() {
            for c in &mut prim {
                *c *= &c_gcd;
            }
        }
        prim
    };

    if pa.len() <= 1 || pb.len() <= 1 {
        return finish(vec![BigInt::one()]);
    }
    if pa == pb {
        return finish(pa);
    }
    // Quick divisibility check when one polynomial can divide the other.
    if pa.len() >= pb.len() && univariate::exact_div(&pa, &pb).is_some() {
        return finish(pb);
    }
    if pb.len() > pa.len() && univariate::exact_div(&pb, &pa).is_some() {
        return finish(pa);
    }

    let lc_a = pa.last().cloned().unwrap_or_else(BigInt::one);
    let lc_b = pb.last().cloned().unwrap_or_else(BigInt::one);
    let gamma = lc_a.gcd(&lc_b);
    let bound_bits = mignotte_bound_bits(&pa, &pb, &gamma);

    let mut best_deg = pa.len().min(pb.len()) - 1;
    let mut crt_residues: Vec<BigInt> = Vec::new();
    let mut crt_modulus = BigInt::zero();
    let mut prime_idx = 0_usize;
    let mut stable_steps = 0_u32;

    let max_deg = pa.len().max(pb.len());
    let batch_size = if execution::should_split_work(4, max_deg * max_deg) {
        4
    } else {
        1
    };

    loop {
        // Collect a batch of admissible primes not dividing lc(A) or lc(B).
        let mut primes = Vec::with_capacity(batch_size);
        while primes.len() < batch_size {
            let p = nth_word_prime(prime_idx);
            prime_idx += 1;
            let p_big = BigInt::from(p);
            if !lc_a.mod_floor(&p_big).is_zero() && !lc_b.mod_floor(&p_big).is_zero() {
                primes.push(p);
            }
        }

        let eval_prime = |&p: &u64| -> (u64, Vec<u64>) {
            let ap = fp_reduce_z(&pa, p);
            let bp = fp_reduce_z(&pb, p);
            let gp = gcd_fp(&ap, &bp, p);
            let gamma_p = gamma.mod_floor(&BigInt::from(p)).to_u64().unwrap_or(0);
            let scaled: Vec<u64> = gp.into_iter().map(|c| fp_mul(c, gamma_p, p)).collect();
            (p, scaled)
        };

        let images: Vec<(u64, Vec<u64>)> = if primes.len() > 1 {
            par::map_slice(&primes, 1, eval_prime)
        } else {
            primes.iter().map(eval_prime).collect()
        };

        for (p, gp) in images {
            let deg = gp.len().saturating_sub(1);
            if deg == 0 {
                return finish(vec![BigInt::one()]);
            }
            if deg > best_deg {
                // Unlucky prime: modular GCD has inflated degree.
                continue;
            }
            if deg < best_deg {
                // All previous primes were unlucky; restart CRT at this smaller degree.
                best_deg = deg;
                crt_residues.clear();
                crt_modulus = BigInt::zero();
                stable_steps = 0;
            }
            let changed = crt_step_poly(&mut crt_residues, &mut crt_modulus, &gp, p);
            if changed {
                stable_steps = 0;
            } else {
                stable_steps += 1;
            }

            if stable_steps >= 1 || crt_modulus.bits() > bound_bits {
                let candidate = centered_poly(&crt_residues, &crt_modulus);
                let (_, prim_cand) = z_content_and_primitive(&candidate);
                if prim_cand.len() == best_deg + 1
                    && univariate::exact_div(&pa, &prim_cand).is_some()
                    && univariate::exact_div(&pb, &prim_cand).is_some()
                {
                    return finish(prim_cand);
                }
                if crt_modulus.bits() > bound_bits + 64 {
                    // Restart if a rare leading-coefficient factor collision occurred.
                    crt_residues.clear();
                    crt_modulus = BigInt::zero();
                    stable_steps = 0;
                }
            }
        }
    }
}

/// Monic greatest common divisor of two univariate rational polynomials
/// over $\mathbb{Q}\[x\]$ via Brown's modular GCD.
#[must_use]
pub fn gcd_q(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let mut ta = a.to_vec();
    let mut tb = b.to_vec();
    while ta.last().is_some_and(Zero::is_zero) {
        ta.pop();
    }
    while tb.last().is_some_and(Zero::is_zero) {
        tb.pop();
    }
    if ta.is_empty() {
        return univariate::monic(&tb);
    }
    if tb.is_empty() {
        return univariate::monic(&ta);
    }
    if ta.len() <= 1 || tb.len() <= 1 {
        return vec![BigRational::one()];
    }
    // For very small degree polynomials with small coefficients, check if modular is needed.
    let (_, za) = univariate::primitive(&ta);
    let (_, zb) = univariate::primitive(&tb);
    let gz = gcd_z(&za, &zb);
    let q: QPoly = gz.into_iter().map(BigRational::from_integer).collect();
    univariate::monic(&q)
}

// ---------------------------------------------------------------------------
// Zippel Sparse / Evaluation-Interpolation Multivariate Modular GCD
// ---------------------------------------------------------------------------

/// Builds a [`Poly`] from an iterator of `(mono, coeff)` terms via balanced
/// pairwise summation (`O(N log N)` instead of `O(N^2)` clones).
fn poly_from_terms<I>(terms: I) -> Poly
where
    I: IntoIterator<Item = (Mono, Number)>,
{
    let mut nodes: Vec<Poly> = terms
        .into_iter()
        .filter(|(_, c)| !c.is_zero())
        .map(|(m, c)| Poly::monomial(m, c))
        .collect();
    if nodes.is_empty() {
        return Poly::zero();
    }
    while nodes.len() > 1 {
        let mut next = Vec::with_capacity(nodes.len().div_ceil(2));
        let mut iter = nodes.into_iter();
        while let Some(a) = iter.next() {
            if let Some(b) = iter.next() {
                next.push(a.add(&b));
            } else {
                next.push(a);
            }
        }
        nodes = next;
    }
    nodes.pop().unwrap_or_else(Poly::zero)
}

/// Clears denominators of a `Poly`, returning `(rational_content, primitive_z_terms)`
/// where `poly = rational_content * sum(c_i * mono_i)` with integer `c_i` having
/// unit GCD and positive leading coefficient.
pub(crate) fn poly_primitive_z(poly: &Poly) -> Option<(BigRational, Poly)> {
    if poly.is_zero() {
        return Some((BigRational::zero(), Poly::zero()));
    }
    let mut lcm_den = BigInt::one();
    let mut rats: Vec<(Mono, BigRational)> = Vec::with_capacity(poly.len());
    for (mono, coeff) in poly.terms() {
        let r = coeff.to_rational()?;
        lcm_den = lcm_den.lcm(r.denom());
        rats.push((mono.clone(), r));
    }
    let mut gcd_num = BigInt::zero();
    let mut ints: Vec<(Mono, BigInt)> = Vec::with_capacity(rats.len());
    for (mono, r) in rats {
        let scaled = (r * &lcm_den).to_integer();
        gcd_num = gcd_num.gcd(&scaled);
        ints.push((mono, scaled));
    }
    if gcd_num.is_zero() {
        return Some((BigRational::zero(), Poly::zero()));
    }
    if poly
        .leading()
        .and_then(|(_, c)| c.to_rational())
        .is_some_and(|r| r.is_negative())
    {
        gcd_num = -gcd_num;
    }
    let out = poly_from_terms(ints.into_iter().map(|(mono, c)| {
        let prim_c = c / &gcd_num;
        (mono, Number::rat(BigRational::from_integer(prim_c)))
    }));
    Some((BigRational::new(gcd_num, lcm_den), out))
}

/// Evaluates all variables except `main_var` at the given $\mathbb{F}_p$ values,
/// returning a dense univariate polynomial in `main_var` over $\mathbb{F}_p$.
fn eval_non_main_fp(
    poly: &Poly,
    main_var: u32,
    point: &[(u32, u64)],
    p: u64,
) -> Option<Vec<u64>> {
    let deg = poly.degree_in(main_var) as usize;
    let mut dense = vec![0_u64; deg + 1];
    let modulus = BigInt::from(p);
    for (mono, coeff) in poly.terms() {
        let r = coeff.to_rational()?;
        if !r.is_integer() {
            return None;
        }
        let n = r.numer();
        let mut c = n.to_i64().map_or_else(
            || n.mod_floor(&modulus).to_u64().unwrap_or(0),
            |v| v.rem_euclid(p as i64) as u64,
        );
        let mut main_exp = 0_usize;
        for &(g, e) in mono {
            if g == main_var {
                main_exp = e as usize;
            } else {
                let val = point
                    .iter()
                    .find(|&&(v, _)| v == g)
                    .map_or(1, |&(_, val)| val);
                c = fp_mul(c, fp_pow(val, u64::from(e), p), p);
            }
        }
        dense[main_exp] = fp_add(dense[main_exp], c, p);
    }
    trim_fp(&mut dense);
    Some(dense)
}

/// Evaluates a multivariate integer `Poly` over $\mathbb{F}_p$ at a full point.
fn eval_full_fp(
    poly: &Poly,
    point: &[(u32, u64)],
    p: u64,
) -> Option<u64> {
    let modulus = BigInt::from(p);
    let mut total = 0_u64;
    for (mono, coeff) in poly.terms() {
        let r = coeff.to_rational()?;
        if !r.is_integer() {
            return None;
        }
        let n = r.numer();
        let mut c = n.to_i64().map_or_else(
            || n.mod_floor(&modulus).to_u64().unwrap_or(0),
            |v| v.rem_euclid(p as i64) as u64,
        );
        for &(g, e) in mono {
            let val = point
                .iter()
                .find(|&&(v, _)| v == g)
                .map_or(1, |&(_, val)| val);
            c = fp_mul(c, fp_pow(val, u64::from(e), p), p);
        }
        total = fp_add(total, c, p);
    }
    Some(total)
}

/// Newton univariate interpolation over $\mathbb{F}_p$: given distinct nodes
/// `xs[0..d]` and values `ys[0..d]`, returns the unique polynomial of degree `< d`
/// in ascending degree order.
fn newton_interpolate_fp(
    xs: &[u64],
    ys: &[u64],
    p: u64,
) -> Option<Vec<u64>> {
    let n = xs.len();
    if n != ys.len() || n == 0 {
        return None;
    }
    let mut dd = ys.to_vec();
    for j in 1..n {
        for i in (j..n).rev() {
            let num = fp_sub(dd[i], dd[i - 1], p);
            let den = fp_sub(xs[i], xs[i - j], p);
            if den == 0 {
                return None;
            }
            dd[i] = fp_mul(num, fp_inv(den, p), p);
        }
    }
    // Convert Newton basis sum_{j} dd[j] prod_{k < j} (x - xs[k]) to standard basis.
    let mut poly = vec![dd[n - 1]];
    for j in (0..n - 1).rev() {
        // poly = poly * (x - xs[j]) + dd[j]
        let neg_x = fp_sub(0, xs[j], p);
        let mut next = vec![0_u64; poly.len() + 1];
        for (i, &c) in poly.iter().enumerate() {
            next[i + 1] = fp_add(next[i + 1], c, p);
            next[i] = fp_add(next[i], fp_mul(c, neg_x, p), p);
        }
        next[0] = fp_add(next[0], dd[j], p);
        poly = next;
    }
    trim_fp(&mut poly);
    Some(poly)
}

enum ModularGcdFp {
    ZeroDeg,
    Terms(Vec<(Mono, u64)>),
}

/// Dense / sparse recursive evaluation-interpolation of $\gcd(A, B)$ over
/// $\mathbb{F}_p[`x_0`, \dots, `x_k`]$, normalized by `gamma = gcd(lc(A), lc(B))`.
///
/// `main_var` is the primary variable; `aux_vars` are the remaining variables
/// interpolated sequentially via Newton interpolation, with immediate
/// degree-zero short-circuiting.
#[allow(clippy::too_many_arguments)]
fn gcd_multivariate_fp(
    a: &Poly,
    b: &Poly,
    lc_a: &Poly,
    lc_b: &Poly,
    gamma: &Poly,
    main_var: u32,
    aux_vars: &[u32],
    deg_bounds: &[u32],
    p: u64,
    seed: u64,
) -> Option<ModularGcdFp> {
    if aux_vars.is_empty() {
        let ap = eval_non_main_fp(a, main_var, &[], p)?;
        let bp = eval_non_main_fp(b, main_var, &[], p)?;
        let gp = gcd_fp(&ap, &bp, p);
        if gp.len() <= 1 {
            return Some(ModularGcdFp::ZeroDeg);
        }
        let gamma_p = eval_full_fp(gamma, &[], p)?;
        if gamma_p == 0 {
            return None;
        }
        let mut terms = Vec::new();
        for (e, c) in gp.into_iter().enumerate() {
            let scaled = fp_mul(c, gamma_p, p);
            if scaled != 0 {
                let mono = if e == 0 {
                    Vec::new()
                } else {
                    vec![(main_var, u32::try_from(e).ok()?)]
                };
                terms.push((mono, scaled));
            }
        }
        return Some(ModularGcdFp::Terms(terms));
    }

    gcd_interp_rec_fp(
        a,
        b,
        lc_a,
        lc_b,
        gamma,
        main_var,
        aux_vars,
        deg_bounds,
        &[],
        p,
        seed,
    )
}

#[allow(clippy::too_many_arguments)]
fn gcd_interp_rec_fp(
    a: &Poly,
    b: &Poly,
    lc_a: &Poly,
    lc_b: &Poly,
    gamma: &Poly,
    main_var: u32,
    rem_vars: &[u32],
    rem_bounds: &[u32],
    fixed: &[(u32, u64)],
    p: u64,
    seed: u64,
) -> Option<ModularGcdFp> {
    let Some((&var, tail_vars)) = rem_vars.split_first() else {
        if eval_full_fp(lc_a, fixed, p)? == 0 || eval_full_fp(lc_b, fixed, p)? == 0 {
            return None;
        }
        let gamma_val = eval_full_fp(gamma, fixed, p)?;
        if gamma_val == 0 {
            return None;
        }
        let ap = eval_non_main_fp(a, main_var, fixed, p)?;
        let bp = eval_non_main_fp(b, main_var, fixed, p)?;
        let gp = gcd_fp(&ap, &bp, p);
        if gp.len() <= 1 {
            return Some(ModularGcdFp::ZeroDeg);
        }
        let mut out = Vec::new();
        for (e, c) in gp.into_iter().enumerate() {
            let scaled = fp_mul(c, gamma_val, p);
            if scaled != 0 {
                let mono = if e == 0 {
                    Vec::new()
                } else {
                    vec![(main_var, u32::try_from(e).ok()?)]
                };
                out.push((mono, scaled));
            }
        }
        return Some(ModularGcdFp::Terms(out));
    };

    let bound = rem_bounds.first().copied().unwrap_or(0) as usize;
    let tail_bounds = &rem_bounds[1..];
    let needed = bound + 1;

    let mut xs = Vec::with_capacity(needed);
    let mut fibres: Vec<Vec<(Mono, u64)>> = Vec::with_capacity(needed);
    let mut target_main_deg: Option<u32> = None;
    let mut candidate_x =
        (seed.wrapping_mul(131).wrapping_add(u64::from(var) * 17 + 2)) % (p - 3) + 2;
    let mut attempts = 0_usize;

    while xs.len() < needed {
        attempts += 1;
        if attempts > needed * 8 + 16 {
            return None;
        }
        let x_val = candidate_x;
        candidate_x = (candidate_x.wrapping_add(7919)) % (p - 3) + 2;
        if xs.contains(&x_val) {
            continue;
        }
        let mut next_fixed = fixed.to_vec();
        next_fixed.push((var, x_val));
        let Some(fibre_res) = gcd_interp_rec_fp(
            a,
            b,
            lc_a,
            lc_b,
            gamma,
            main_var,
            tail_vars,
            tail_bounds,
            &next_fixed,
            p,
            seed.wrapping_add(attempts as u64),
        ) else {
            continue;
        };
        let fibre = match fibre_res {
            | ModularGcdFp::ZeroDeg => return Some(ModularGcdFp::ZeroDeg),
            | ModularGcdFp::Terms(terms) => terms,
        };
        let main_deg = fibre
            .iter()
            .filter_map(|(m, _)| m.iter().find(|&&(g, _)| g == main_var).map(|&(_, e)| e))
            .max()
            .unwrap_or(0);
        match target_main_deg {
            | None => {
                target_main_deg = Some(main_deg);
                xs.push(x_val);
                fibres.push(fibre);
            },
            | Some(d) if main_deg == d => {
                xs.push(x_val);
                fibres.push(fibre);
            },
            | Some(d) if main_deg < d => {
                // Previous evaluations were unlucky and inflated main_var degree.
                target_main_deg = Some(main_deg);
                xs.clear();
                fibres.clear();
                xs.push(x_val);
                fibres.push(fibre);
            },
            | Some(_) => {
                // Unlucky evaluation point; skip.
            },
        }
    }

    // Collect the union of all monomials present in the fibres.
    let mut all_monos: Vec<Mono> = fibres
        .iter()
        .flat_map(|f| f.iter().map(|(m, _)| m.clone()))
        .collect();
    all_monos.sort_unstable();
    all_monos.dedup();

    let mut result = Vec::new();
    for base_mono in all_monos {
        let ys: Vec<u64> = fibres
            .iter()
            .map(|f| {
                f.iter()
                    .find(|(m, _)| *m == base_mono)
                    .map_or(0, |&(_, c)| c)
            })
            .collect();
        if ys.iter().all(|&y| y == 0) {
            continue;
        }
        let interp = newton_interpolate_fp(&xs, &ys, p)?;
        for (e, c) in interp.into_iter().enumerate() {
            if c == 0 {
                continue;
            }
            let mut full_mono = base_mono.clone();
            if e > 0 {
                full_mono.push((var, u32::try_from(e).ok()?));
                full_mono.sort_unstable();
            }
            result.push((full_mono, c));
        }
    }
    Some(ModularGcdFp::Terms(result))
}

/// Fast probabilistic coprimality certificate: if at one admissible modular
/// evaluation point where neither leading coefficient vanishes the univariate
/// GCD in `main_var` is `1`, then `deg_main(gcd(prim_a, prim_b)) == 0`.
fn quick_coprime_probe(
    prim_a: &Poly,
    prim_b: &Poly,
    lc_a: &Poly,
    lc_b: &Poly,
    main_var: u32,
    aux_vars: &[u32],
) -> bool {
    let p = WORD_PRIMES[0];
    for attempt in 0_u64..4 {
        let fixed: Vec<(u32, u64)> = aux_vars
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                let val =
                    ((attempt + 1) * 1009 + u64::from(v) * 131 + (i as u64) * 37 + 5) % (p - 3) + 2;
                (v, val)
            })
            .collect();
        let Some(la) = eval_full_fp(lc_a, &fixed, p) else {
            return false;
        };
        let Some(lb) = eval_full_fp(lc_b, &fixed, p) else {
            return false;
        };
        if la == 0 || lb == 0 {
            continue;
        }
        let Some(ap) = eval_non_main_fp(prim_a, main_var, &fixed, p) else {
            return false;
        };
        let Some(bp) = eval_non_main_fp(prim_b, main_var, &fixed, p) else {
            return false;
        };
        if gcd_fp(&ap, &bp, p).len() <= 1 {
            return true;
        }
        return false;
    }
    false
}

/// Removes content in `other_vars` from a vector of coefficient polynomials in
/// `main_var`.
fn upoly_primitive_content(u: &[Poly]) -> Option<Vec<Poly>> {
    let mut trimmed = u.to_vec();
    while trimmed.last().is_some_and(Poly::is_zero) {
        trimmed.pop();
    }
    if trimmed.is_empty() {
        return Some(Vec::new());
    }
    if trimmed.len() == 1 {
        return Some(vec![Poly::constant(Number::from(1))]);
    }
    let mut content = Poly::zero();
    for c in trimmed.iter().filter(|c| !c.is_zero()) {
        if content.is_zero() {
            content = c.clone();
        } else {
            content = poly_gcd(&content, c)?;
            if content.as_constant().is_some_and(|n| n.is_one()) {
                return Some(trimmed);
            }
        }
    }
    if content.is_zero() || content.as_constant().is_some_and(|n| n.is_one()) {
        return Some(trimmed);
    }
    let (_, content) = poly_primitive_z(&content)?;
    if content.as_constant().is_some_and(|n| n.is_one()) {
        return Some(trimmed);
    }
    trimmed
        .into_iter()
        .map(|c| {
            if c.is_zero() {
                Some(Poly::zero())
            } else {
                c.div_exact(&content, CAP)
            }
        })
        .collect()
}

/// Fraction-free pseudo-remainder in `main_var` over `Z[aux_vars]`.
fn upoly_prem_z(
    u: &[Poly],
    v: &[Poly],
) -> Option<Vec<Poly>> {
    let lead = v.last()?;
    if u.len() < v.len() {
        return Some(u.to_vec());
    }
    let mut rem = u.to_vec();
    let q_len = u.len() - v.len() + 1;
    let one_const = |p: &Poly| p.as_constant().is_some_and(|n| n.is_one());
    for k in (0..q_len).rev() {
        let top = rem[k + v.len() - 1].clone();
        if top.is_zero() {
            continue;
        }
        let g = if one_const(lead) || one_const(&top) {
            Poly::constant(Number::from(1))
        } else {
            poly_gcd(&top, lead)?
        };
        let (step, factor) = if one_const(&g) {
            (lead.clone(), top)
        } else {
            (lead.div_exact(&g, CAP)?, top.div_exact(&g, CAP)?)
        };
        if !one_const(&step) {
            for slot in &mut rem[..k + v.len() - 1] {
                if !slot.is_zero() {
                    *slot = slot.mul(&step, CAP)?;
                }
            }
        }
        rem[k + v.len() - 1] = Poly::zero();
        for j in 0..(v.len() - 1) {
            if !v[j].is_zero() {
                let prod = factor.mul(&v[j], CAP)?;
                rem[k + j] = rem[k + j].sub(&prod);
            }
        }
        if rem.iter().map(Poly::len).sum::<usize>() > CAP {
            return None;
        }
    }
    while rem.last().is_some_and(Poly::is_zero) {
        rem.pop();
    }
    Some(rem)
}

/// Primitive PRS fallback for multivariate GCD when modular grid or CRT bounds
/// are exceeded. Both `a` and `b` must be primitive in `main_var`.
fn poly_gcd_prs_main_var(
    a: &Poly,
    b: &Poly,
    main_var: u32,
) -> Option<Poly> {
    let trim_u = |mut v: Vec<Poly>| {
        while v.last().is_some_and(Poly::is_zero) {
            v.pop();
        }
        v
    };
    let mut u = trim_u(a.coefficients_in(main_var));
    let mut v = trim_u(b.coefficients_in(main_var));
    if u.len() < v.len() {
        std::mem::swap(&mut u, &mut v);
    }
    while v.len() > 1 {
        let rem = upoly_prem_z(&u, &v)?;
        if rem.is_empty() {
            break;
        }
        if rem.len() == 1 {
            return Some(Poly::constant(Number::from(1)));
        }
        u = v;
        v = upoly_primitive_content(&rem)?;
    }
    if v.len() <= 1 {
        return Some(Poly::constant(Number::from(1)));
    }
    let mut terms = Vec::new();
    for (power, coeff) in v.into_iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let exp = u32::try_from(power).ok()?;
        for (mono, num) in coeff.terms() {
            let mut full = mono.clone();
            if exp > 0 {
                full.push((main_var, exp));
                full.sort_unstable();
            }
            terms.push((full, num.clone()));
        }
    }
    let (_, prim) = poly_primitive_z(&poly_from_terms(terms))?;
    Some(prim)
}

/// Removes the content in `other_vars` (viewed as coefficients of `main_var`)
/// from a multivariate integer polynomial `p`.
fn primitive_in_main_var(
    p: &Poly,
    main_var: u32,
) -> Option<(Poly, Poly)> {
    let coeffs = p.coefficients_in(main_var);
    let mut content = Poly::zero();
    for c in coeffs.iter().filter(|c| !c.is_zero()) {
        if content.is_zero() {
            content = c.clone();
        } else {
            content = poly_gcd(&content, c)?;
            if content.as_constant().is_some_and(|n| n.is_one()) {
                break;
            }
        }
    }
    if content.is_zero() {
        return Some((Poly::constant(Number::from(1)), p.clone()));
    }
    let (_, content) = poly_primitive_z(&content)?;
    if content.as_constant().is_some_and(|n| n.is_one()) {
        return Some((content, p.clone()));
    }
    let prim = p.div_exact(&content, CAP)?;
    Some((content, prim))
}

/// Computes the exact greatest common divisor of multivariate rational polynomials.
///
/// Uses Zippel/Brown modular evaluation-interpolation and Chinese remainder reconstruction.
///
/// The returned GCD has primitive integer coefficients (or `1` when coprime)
/// with positive leading coefficient.
#[must_use]
pub fn poly_gcd(
    a: &Poly,
    b: &Poly,
) -> Option<Poly> {
    if a.is_zero() {
        return poly_primitive_z(b).map(|(_, p)| normalise(p));
    }
    if b.is_zero() {
        return poly_primitive_z(a).map(|(_, p)| normalise(p));
    }
    if !a.is_exact() || !b.is_exact() {
        return None;
    }

    // Extract monomial contents first.
    let (ma, mb) = (a.monomial_content(), b.monomial_content());
    let common_mono: Mono = ma
        .iter()
        .filter_map(|&(g, ea)| {
            mb.iter()
                .find(|&&(h, _)| h == g)
                .map(|&(_, eb)| (g, ea.min(eb)))
        })
        .collect();
    let a_red = if ma.is_empty() {
        a.clone()
    } else {
        a.div_monomial(&ma)
    };
    let b_red = if mb.is_empty() {
        b.clone()
    } else {
        b.div_monomial(&mb)
    };

    let (_, za) = poly_primitive_z(&a_red)?;
    let (_, zb) = poly_primitive_z(&b_red)?;

    let finish = |g: Poly| -> Poly {
        let g = normalise(g);
        if common_mono.is_empty() {
            g
        } else {
            g.mul(&Poly::monomial(common_mono.clone(), Number::from(1)), CAP)
                .unwrap_or(g)
        }
    };

    if za.as_constant().is_some() || zb.as_constant().is_some() {
        return Some(finish(Poly::constant(Number::from(1))));
    }
    if za == zb {
        return Some(finish(za));
    }
    if za.div_exact(&zb, CAP).is_some() {
        return Some(finish(zb));
    }
    if zb.div_exact(&za, CAP).is_some() {
        return Some(finish(za));
    }

    let supp_a = za.support();
    let supp_b = zb.support();
    let shared_vars: Vec<u32> = supp_a
        .iter()
        .copied()
        .filter(|v| supp_b.contains(v))
        .collect();
    if shared_vars.is_empty() {
        return Some(finish(Poly::constant(Number::from(1))));
    }

    // If any variable appears only in za or only in zb, factor out the content
    // with respect to that variable first!
    for &v in &supp_a {
        if !supp_b.contains(&v) {
            let coeffs = za.coefficients_in(v);
            let mut g = zb;
            for c in coeffs.into_iter().filter(|c| !c.is_zero()) {
                g = poly_gcd(&g, &c)?;
                if g.as_constant().is_some_and(|n| n.is_one()) {
                    return Some(finish(g));
                }
            }
            return Some(finish(g));
        }
    }
    for &v in &supp_b {
        if !supp_a.contains(&v) {
            let coeffs = zb.coefficients_in(v);
            let mut g = za;
            for c in coeffs.into_iter().filter(|c| !c.is_zero()) {
                g = poly_gcd(&g, &c)?;
                if g.as_constant().is_some_and(|n| n.is_one()) {
                    return Some(finish(g));
                }
            }
            return Some(finish(g));
        }
    }

    // Univariate fast path via Brown's modular GCD.
    if let &[only] = shared_vars.as_slice() {
        let ua: ZPoly = za
            .univariate_in(only)?
            .iter()
            .map(|n| Some(n.to_rational()?.to_integer()))
            .collect::<Option<_>>()?;
        let ub: ZPoly = zb
            .univariate_in(only)?
            .iter()
            .map(|n| Some(n.to_rational()?.to_integer()))
            .collect::<Option<_>>()?;
        let gz = gcd_z(&ua, &ub);
        let nums: Vec<Number> = gz
            .into_iter()
            .map(|c| Number::rat(BigRational::from_integer(c)))
            .collect();
        return Some(finish(Poly::from_univariate(only, &nums)));
    }

    // Choose main_var as the variable of highest maximum degree, and extract
    // multivariate contents in the remaining variables.
    let main_var = shared_vars
        .iter()
        .copied()
        .max_by_key(|&v| za.degree_in(v).max(zb.degree_in(v)))?;
    let (cont_a, prim_a) = primitive_in_main_var(&za, main_var)?;
    let (cont_b, prim_b) = primitive_in_main_var(&zb, main_var)?;
    let content_gcd = poly_gcd(&cont_a, &cont_b)?;

    let lc_a = prim_a.coefficients_in(main_var).pop()?;
    let lc_b = prim_b.coefficients_in(main_var).pop()?;
    let gamma = poly_gcd(&lc_a, &lc_b)?;

    let aux_vars: Vec<u32> = shared_vars
        .iter()
        .copied()
        .filter(|&v| v != main_var)
        .collect();
    let deg_bounds: Vec<u32> = aux_vars
        .iter()
        .map(|&v| {
            prim_a
                .degree_in(v)
                .min(prim_b.degree_in(v))
                .saturating_add(gamma.degree_in(v))
        })
        .collect();

    // Fast coprimality check at a single modular evaluation point before interpolation.
    if quick_coprime_probe(&prim_a, &prim_b, &lc_a, &lc_b, main_var, &aux_vars) {
        return Some(finish(content_gcd));
    }

    let fallback_prs = || -> Option<Poly> {
        let g = poly_gcd_prs_main_var(&prim_a, &prim_b, main_var)?;
        let full = g.mul(&content_gcd, CAP)?;
        let (_, full_prim) = poly_primitive_z(&full)?;
        Some(finish(full_prim))
    };

    // Total interpolation grid size; fall back to primitive PRS if modular grid blows up.
    let grid_points: usize = deg_bounds
        .iter()
        .try_fold(1_usize, |acc, &d| acc.checked_mul(d as usize + 1))
        .unwrap_or(usize::MAX);
    if grid_points > 4096 {
        return fallback_prs();
    }

    let mut best_main_deg = prim_a.degree_in(main_var).min(prim_b.degree_in(main_var));
    let mut crt_monos: Vec<Mono> = Vec::new();
    let mut crt_residues: Vec<BigInt> = Vec::new();
    let mut crt_modulus = BigInt::zero();

    for prime_idx in 0..32 {
        let p = nth_word_prime(prime_idx);
        let Some(res_fp) = gcd_multivariate_fp(
            &prim_a,
            &prim_b,
            &lc_a,
            &lc_b,
            &gamma,
            main_var,
            &aux_vars,
            &deg_bounds,
            p,
            (prime_idx as u64) + 1,
        ) else {
            continue;
        };
        let terms_fp = match res_fp {
            | ModularGcdFp::ZeroDeg => return Some(finish(content_gcd)),
            | ModularGcdFp::Terms(terms) => terms,
        };
        let main_deg = terms_fp
            .iter()
            .filter_map(|(m, _)| m.iter().find(|&&(g, _)| g == main_var).map(|&(_, e)| e))
            .max()
            .unwrap_or(0);
        if main_deg == 0 {
            return Some(finish(content_gcd));
        }
        if main_deg > best_main_deg {
            continue;
        }
        if main_deg < best_main_deg {
            best_main_deg = main_deg;
            crt_monos.clear();
            crt_residues.clear();
            crt_modulus = BigInt::zero();
        }

        // Align monomials between CRT accumulator and new modular image.
        let mut all_monos = crt_monos.clone();
        for (m, _) in &terms_fp {
            if !all_monos.contains(m) {
                all_monos.push(m.clone());
            }
        }
        all_monos.sort_unstable();
        let mut aligned_acc = vec![BigInt::zero(); all_monos.len()];
        for (m, r) in crt_monos.iter().zip(crt_residues.iter()) {
            if let Ok(idx) = all_monos.binary_search(m) {
                aligned_acc[idx].clone_from(r);
            }
        }
        let mut aligned_new = vec![0_u64; all_monos.len()];
        for (m, c) in &terms_fp {
            if let Ok(idx) = all_monos.binary_search(m) {
                aligned_new[idx] = *c;
            }
        }
        crt_monos = all_monos;
        crt_residues = aligned_acc;
        let changed = crt_step_poly(&mut crt_residues, &mut crt_modulus, &aligned_new, p);

        if !changed || prime_idx % 4 == 3 {
            let half = &crt_modulus >> 1;
            let candidate = poly_from_terms(crt_monos.iter().zip(crt_residues.iter()).filter_map(
                |(m, r)| {
                    let c = symmetric_mod(r, &crt_modulus, &half);
                    (!c.is_zero()).then(|| (m.clone(), Number::rat(BigRational::from_integer(c))))
                },
            ));
            if let Some((_, prim_cand)) = primitive_in_main_var(&candidate, main_var) {
                let (_, prim_cand) = poly_primitive_z(&prim_cand)?;
                if prim_cand.degree_in(main_var) == best_main_deg
                    && prim_a.div_exact(&prim_cand, CAP).is_some()
                    && prim_b.div_exact(&prim_cand, CAP).is_some()
                {
                    let full = prim_cand.mul(&content_gcd, CAP)?;
                    let (_, full_prim) = poly_primitive_z(&full)?;
                    return Some(finish(full_prim));
                }
            }
        }
    }

    fallback_prs()
}

// ---------------------------------------------------------------------------
// Multivariate Kronecker Factorization with Multi-Point Pruning & Parallelism
// ---------------------------------------------------------------------------

/// Mixed-radix weights for the variables `vars` of `f`.
fn weights(
    f: &Poly,
    vars: &[u32],
) -> Option<Vec<u64>> {
    let mut out = Vec::with_capacity(vars.len());
    let mut w: u64 = 1;
    for &v in vars {
        out.push(w);
        w = w.checked_mul(u64::from(f.degree_in(v)) + 1)?;
    }
    (w <= MAX_IMAGE_DEGREE + 1).then_some(out)
}

/// `f(t^(w_1), ..., t^(w_n))` as dense rational coefficients.
fn image(
    f: &Poly,
    vars: &[u32],
    w: &[u64],
) -> Option<Vec<BigRational>> {
    let mut dense: Vec<BigRational> = Vec::new();
    for (mono, c) in f.terms() {
        let mut e: u64 = 0;
        for &(g, k) in mono {
            let i = vars.iter().position(|&v| v == g)?;
            e += u64::from(k) * w[i];
        }
        let e = usize::try_from(e).ok()?;
        if dense.len() <= e {
            dense.resize(e + 1, BigRational::zero());
        }
        dense[e] += c.to_rational()?;
    }
    Some(dense)
}

/// Recovers the multivariate polynomial whose Kronecker image is `z`, returning
/// `None` immediately if any individual variable degree exceeds `max_degs`.
fn preimage_checked(
    z: &[BigInt],
    vars: &[u32],
    w: &[u64],
    max_degs: &[u32],
) -> Option<Poly> {
    let mut terms = Vec::new();
    for (e, c) in z.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut rest = e as u64;
        let mut mono = Vec::new();
        for i in (0..vars.len()).rev() {
            let k = rest / w[i];
            rest %= w[i];
            let k_u32 = u32::try_from(k).ok()?;
            if k_u32 > max_degs[i] {
                return None;
            }
            if k_u32 > 0 {
                mono.push((vars[i], k_u32));
            }
        }
        mono.sort_unstable();
        terms.push((mono, Number::rat(BigRational::from_integer(c.clone()))));
    }
    Some(poly_from_terms(terms))
}

fn z_mul(
    a: &[BigInt],
    b: &[BigInt],
) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

/// Evaluates an integer univariate polynomial `z(t)` at small integer `t`.
fn eval_z_at(
    z: &[BigInt],
    t: i64,
) -> BigInt {
    match t {
        | 0 => z.first().cloned().unwrap_or_else(BigInt::zero),
        | 1 => z.iter().fold(BigInt::zero(), |acc, c| acc + c),
        | -1 => {
            z.iter().enumerate().fold(BigInt::zero(), |acc, (i, c)| {
                if i % 2 == 0 {
                    acc + c
                } else {
                    acc - c
                }
            })
        },
        | _ => {
            let t_big = BigInt::from(t);
            z.iter()
                .rev()
                .fold(BigInt::zero(), |acc, c| acc * &t_big + c)
        },
    }
}

/// `g` with a positive leading coefficient.
fn normalise(g: Poly) -> Poly {
    match g.leading() {
        | Some((_, c)) if c.to_rational().is_some_and(|r| r.is_negative()) => g.neg(),
        | _ => g,
    }
}

/// Advances `c` to the next `c.len()`-subset of `0..n` in lexicographic
/// order; `false` after the last one.
fn next_combination(
    c: &mut [usize],
    n: usize,
) -> bool {
    let k = c.len();
    let mut i = k;
    while i > 0 {
        i -= 1;
        if c[i] < n - k + i {
            c[i] += 1;
            for j in i + 1..k {
                c[j] = c[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// Precomputed integer divisibility invariants of `f`'s univariate Kronecker image
/// for fast candidate pruning before multivariate `preimage` and `div_exact`.
struct PruningInvariants {
    max_degs: Vec<u32>,
    lead_z: BigInt,
    trail_nz: BigInt,
    val_1: BigInt,
    val_neg1: BigInt,
    val_2: BigInt,
}

impl PruningInvariants {
    fn from_poly_and_image(
        f: &Poly,
        vars: &[u32],
        dense_q: &[BigRational],
    ) -> Self {
        let (_, z_img) = univariate::primitive(dense_q);
        let max_degs: Vec<u32> = vars.iter().map(|&v| f.degree_in(v)).collect();
        let lead_z = z_img.last().cloned().unwrap_or_else(BigInt::one);
        let trail_nz = z_img
            .iter()
            .find(|c| !c.is_zero())
            .cloned()
            .unwrap_or_else(BigInt::one);
        let val_1 = eval_z_at(&z_img, 1);
        let val_neg1 = eval_z_at(&z_img, -1);
        let val_2 = eval_z_at(&z_img, 2);
        Self {
            max_degs,
            lead_z,
            trail_nz,
            val_1,
            val_neg1,
            val_2,
        }
    }

    /// Checks whether a candidate univariate factor product can possibly lift
    /// to a true divisor of `f`.
    fn accepts_product(
        &self,
        prod_lead: &BigInt,
        prod_trail: &BigInt,
        prod_v1: &BigInt,
        prod_vm1: &BigInt,
        prod_v2: &BigInt,
    ) -> bool {
        if prod_lead.is_zero() || !(&self.lead_z % prod_lead).is_zero() {
            return false;
        }
        if prod_trail.is_zero() || !(&self.trail_nz % prod_trail).is_zero() {
            return false;
        }
        if !self.val_1.is_zero() && (prod_v1.is_zero() || !(&self.val_1 % prod_v1).is_zero()) {
            return false;
        }
        if !self.val_neg1.is_zero()
            && (prod_vm1.is_zero() || !(&self.val_neg1 % prod_vm1).is_zero())
        {
            return false;
        }
        if !self.val_2.is_zero() && (prod_v2.is_zero() || !(&self.val_2 % prod_v2).is_zero()) {
            return false;
        }
        true
    }
}

/// Per-piece precomputed evaluations at `1, -1, 2`, leading, and trailing nonzero coefficient.
struct PieceMeta<'a> {
    coeffs: &'a [BigInt],
    lead: BigInt,
    trail_nz: BigInt,
    v1: BigInt,
    vm1: BigInt,
    v2: BigInt,
}

impl<'a> PieceMeta<'a> {
    fn new(coeffs: &'a [BigInt]) -> Self {
        Self {
            coeffs,
            lead: coeffs.last().cloned().unwrap_or_else(BigInt::one),
            trail_nz: coeffs
                .iter()
                .find(|c| !c.is_zero())
                .cloned()
                .unwrap_or_else(BigInt::one),
            v1: eval_z_at(coeffs, 1),
            vm1: eval_z_at(coeffs, -1),
            v2: eval_z_at(coeffs, 2),
        }
    }
}

/// One irreducible factor of the non-constant `f`, or `None` when `f` is
/// irreducible (or too large to split).
fn split_off(
    f: &Poly,
    vars: &[u32],
) -> Option<Poly> {
    let w = weights(f, vars)?;
    let dense = image(f, vars, &w)?;
    let (_, factors) = univariate::factor(&dense);
    let pieces: Vec<PieceMeta<'_>> = factors
        .iter()
        .flat_map(|(z, m)| std::iter::repeat_n(z.as_slice(), *m as usize))
        .map(PieceMeta::new)
        .collect();
    if pieces.len() <= 1 || pieces.len() > MAX_PIECES {
        return None;
    }

    let inv = PruningInvariants::from_poly_and_image(f, vars, &dense);
    let n = pieces.len();
    let mut tried = 0_usize;
    let mut seen: HashSet<Vec<BigInt>> = HashSet::new();

    for size in 1..=n / 2 {
        let mut choice: Vec<usize> = (0..size).collect();
        let mut batch: Vec<Vec<usize>> = Vec::with_capacity(FACTOR_BATCH_SIZE);
        let mut exhausted = false;

        while !exhausted {
            batch.clear();
            while batch.len() < FACTOR_BATCH_SIZE {
                tried += 1;
                if tried > MAX_SUBSETS {
                    exhausted = true;
                    break;
                }
                // Fast scalar pruning before enqueueing subset.
                let mut lead = BigInt::one();
                let mut trail = BigInt::one();
                let mut v1 = BigInt::one();
                let mut vm1 = BigInt::one();
                let mut v2 = BigInt::one();
                for &i in &choice {
                    let pm = &pieces[i];
                    lead *= &pm.lead;
                    trail *= &pm.trail_nz;
                    v1 *= &pm.v1;
                    vm1 *= &pm.vm1;
                    v2 *= &pm.v2;
                }
                if inv.accepts_product(&lead, &trail, &v1, &vm1, &v2) {
                    batch.push(choice.clone());
                }
                if !next_combination(&mut choice, n) {
                    exhausted = true;
                    break;
                }
            }

            if batch.is_empty() {
                if tried > MAX_SUBSETS {
                    return None;
                }
                continue;
            }

            // Deduplicate products in deterministic subset order.
            let mut candidates: Vec<Vec<BigInt>> = Vec::with_capacity(batch.len());
            for subset in &batch {
                let product = subset
                    .iter()
                    .fold(vec![BigInt::one()], |acc, &i| z_mul(&acc, pieces[i].coeffs));
                if !seen.contains(&product) && !candidates.contains(&product) {
                    candidates.push(product);
                }
            }

            let test_candidate = |product: &[BigInt]| -> Option<Poly> {
                let g = preimage_checked(product, vars, &w, &inv.max_degs)?;
                if g.as_constant().is_none() && f.div_exact(&g, CAP).is_some() {
                    Some(normalise(g))
                } else {
                    None
                }
            };

            if candidates.len() >= 4 && execution::should_split_work(candidates.len(), f.len() * 16)
            {
                let results = par::map_slice(&candidates, 2, |p| test_candidate(p));
                for (prod, res) in candidates.into_iter().zip(results) {
                    if let Some(hit) = res {
                        return Some(hit);
                    }
                    seen.insert(prod);
                }
            } else {
                for prod in candidates {
                    if let Some(hit) = test_candidate(&prod) {
                        return Some(hit);
                    }
                    seen.insert(prod);
                }
            }

            if tried > MAX_SUBSETS {
                return None;
            }
        }
    }
    None
}

/// `f = unit * prod g_i^m_i` with irreducible `g_i` of positive leading
/// coefficient, for a polynomial with exact coefficients in the
/// generators `vars`.
///
/// `None` when the polynomial is too large for the
/// reduction (the caller keeps it unfactored).
#[must_use]
pub fn factor(
    f: &Poly,
    vars: &[u32],
) -> Option<(Number, Vec<(Poly, u32)>)> {
    let mut work = f.clone();
    let mut out: Vec<(Poly, u32)> = Vec::new();
    for _ in 0..64 {
        if let Some(c) = work.as_constant() {
            out.sort_by_key(|(g, _)| (g.total_degree(), g.len()));
            return Some((c, out));
        }
        weights(&work, vars)?;
        let g = split_off(&work, vars).unwrap_or_else(|| normalise(work.clone()));
        let mut m = 0;
        while let Some(q) = work.div_exact(&g, CAP) {
            work = q;
            m += 1;
            if work.as_constant().is_some() {
                break;
            }
        }
        if m == 0 {
            return None;
        }
        out.push((g, m));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zpoly(coeffs: &[i64]) -> ZPoly {
        coeffs.iter().map(|&c| BigInt::from(c)).collect()
    }

    #[test]
    fn word_primes_are_all_prime_and_distinct() {
        let mut seen = HashSet::new();
        for &p in WORD_PRIMES {
            assert!(p < (1_u64 << 31));
            assert!(is_prime_u64(p), "{p} is composite");
            assert!(seen.insert(p), "duplicate prime {p}");
        }
        for idx in WORD_PRIMES.len()..WORD_PRIMES.len() + 8 {
            let p = nth_word_prime(idx);
            assert!(is_prime_u64(p), "fallback prime {p} is composite");
            assert!(seen.insert(p), "duplicate fallback prime {p}");
        }
    }

    #[test]
    fn hgcd_fp_matches_euclid_across_degrees() {
        let p = 2_147_483_647_u64;
        // Construct common factor g of degree 15 and two cofactors of degree 25 and 22.
        let g: Vec<u64> = (0..=15).map(|i| (i * 37 + 5) % p).collect();
        let u: Vec<u64> = (0..=25).map(|i| (i * i * 19 + 3).max(1) % p).collect();
        let v: Vec<u64> = (0..=22).map(|i| (i * 53 + 11).max(1) % p).collect();
        let a = fp_poly_mul(&g, &u, p);
        let b = fp_poly_mul(&g, &v, p);

        let fast = gcd_fp(&a, &b, p);
        let euclid = gcd_fp_euclid(&a, &b, p);
        assert_eq!(fast, euclid);
        assert_eq!(fast, fp_poly_monic(&g, p));

        // High degree Karatsuba + multi-level HGCD (deg 70 and deg 70).
        let g_hi: Vec<u64> = (0..=25).map(|i| (i * 97 + 1) % p).collect();
        let u_hi: Vec<u64> = (0..=45).map(|i| (i * 61 + 7) % p).collect();
        let v_hi: Vec<u64> = (0..=45).map(|i| (i * 83 + 13) % p).collect();
        let a_hi = fp_poly_mul(&g_hi, &u_hi, p);
        let b_hi = fp_poly_mul(&g_hi, &v_hi, p);
        assert_eq!(gcd_fp(&a_hi, &b_hi, p), gcd_fp_euclid(&a_hi, &b_hi, p));
    }

    #[test]
    fn ext_gcd_fp_satisfies_bezout_identity() {
        let p = 2_147_483_629_u64;
        for (da, db) in [(30_u64, 28_u64), (24, 24), (20, 45), (55, 55)] {
            let a: Vec<u64> = (0..=da).map(|i| ((i * 71 + 13) % p).max(1)).collect();
            let b: Vec<u64> = (0..=db).map(|i| ((i * 43 + 7) % p).max(1)).collect();
            let (g, s, t) = ext_gcd_fp(&a, &b, p);
            let lhs = fp_poly_add(&fp_poly_mul(&s, &a, p), &fp_poly_mul(&t, &b, p), p);
            assert_eq!(lhs, g, "deg ({da}, {db})");
            assert_eq!(g, gcd_fp_euclid(&a, &b, p), "deg ({da}, {db})");
        }
        // Multiple of each other at high degree.
        let a: Vec<u64> = (0..=36).map(|i| ((i * 29 + 3) % p).max(1)).collect();
        let b: Vec<u64> = a.iter().map(|&c| fp_mul(c, 5, p)).collect();
        let (g, s, t) = ext_gcd_fp(&a, &b, p);
        assert_eq!(
            fp_poly_add(&fp_poly_mul(&s, &a, p), &fp_poly_mul(&t, &b, p), p),
            g
        );
        assert_eq!(g, fp_poly_monic(&a, p));
    }

    #[test]
    fn brown_modular_gcd_handles_non_monic_and_large_leading_factors() {
        // Common factor with large non-monic leading coefficient: 120*x^4 - 77*x^2 + 35*x - 18
        let g = zpoly(&[-18, 35, -77, 0, 120]);
        let u = zpoly(&[7, -13, 29, -5, 42]);
        let v = zpoly(&[-11, 23, 0, -19, 30]);
        let a = z_mul(&g, &u);
        let b = z_mul(&g, &v);

        let computed = gcd_z(&a, &b);
        assert_eq!(computed, g);

        // High-degree polynomials where rational PRS would explode
        let mut g_hi = vec![BigInt::from(1)];
        for k in 1..=16 {
            g_hi = z_mul(&g_hi, &zpoly(&[-k, 2 * k + 1]));
        }
        let u_hi = z_mul(&zpoly(&[3, -7, 5, 11]), &zpoly(&[-5, 2, -9, 4]));
        let v_hi = z_mul(&zpoly(&[13, 2, -3, 7]), &zpoly(&[1, -4, 8, 3]));
        let a_hi = z_mul(&g_hi, &u_hi);
        let b_hi = z_mul(&g_hi, &v_hi);
        let (_, expected_prim) = z_content_and_primitive(&g_hi);
        assert_eq!(gcd_z(&a_hi, &b_hi), expected_prim);
    }

    #[test]
    fn zippel_multivariate_gcd_recovers_common_factor() {
        // g = 2*x^2*y - 3*x*z + y^2 + 5
        // u = x*y + 2*z - 1
        // v = 3*x - y*z + 2
        let x = Poly::generator(0);
        let y = Poly::generator(1);
        let z = Poly::generator(2);
        let c = |n: i64| Poly::constant(Number::from(n));

        let g = c(2)
            .mul(&x.pow(2, CAP).unwrap(), CAP)
            .unwrap()
            .mul(&y, CAP)
            .unwrap()
            .sub(&c(3).mul(&x, CAP).unwrap().mul(&z, CAP).unwrap())
            .add(&y.pow(2, CAP).unwrap())
            .add(&c(5));
        let u = x
            .mul(&y, CAP)
            .unwrap()
            .add(&c(2).mul(&z, CAP).unwrap())
            .sub(&c(1));
        let v = c(3)
            .mul(&x, CAP)
            .unwrap()
            .sub(&y.mul(&z, CAP).unwrap())
            .add(&c(2));

        let a = g.mul(&u, CAP).unwrap();
        let b = g.mul(&v, CAP).unwrap();
        let computed = poly_gcd(&a, &b).unwrap();
        assert_eq!(computed, normalise(g.clone()));

        // Direct primitive PRS fallback check and coprime check.
        let prs = poly_gcd_prs_main_var(&a, &b, 0).unwrap();
        assert_eq!(prs, normalise(g.clone()));
        assert_eq!(poly_gcd(&u, &v).unwrap(), c(1));

        // Common content in y, z times factor in x, y, z.
        let cont = y.add(&z).add(&c(1));
        let ac = a.mul(&cont, CAP).unwrap();
        let bc = b.mul(&cont, CAP).unwrap();
        let expected = normalise(g.mul(&cont, CAP).unwrap());
        assert_eq!(poly_gcd(&ac, &bc).unwrap(), expected);
    }

    #[test]
    fn multivariate_factor_is_deterministic_across_thread_counts() {
        let x = Poly::generator(0);
        let y = Poly::generator(1);
        let z = Poly::generator(2);
        let c = |n: i64| Poly::constant(Number::from(n));

        let f1 = x.add(&y).add(&z);
        let f2 = c(2).mul(&x, CAP).unwrap().sub(&y).add(&c(3));
        let f3 = x.mul(&z, CAP).unwrap().sub(&y).add(&c(1));
        let f4 = y.pow(2, CAP).unwrap().add(&x).sub(&z);
        let prod = f1
            .mul(&f2, CAP)
            .unwrap()
            .mul(&f3, CAP)
            .unwrap()
            .mul(&f4, CAP)
            .unwrap();

        let baseline = execution::with_scoped_pool(Some(1), || factor(&prod, &[0, 1, 2]));
        let par_2 = execution::with_scoped_pool(Some(2), || factor(&prod, &[0, 1, 2]));
        let par_4 = execution::with_scoped_pool(Some(4), || factor(&prod, &[0, 1, 2]));
        assert!(baseline.is_some());
        assert_eq!(baseline, par_2);
        assert_eq!(baseline, par_4);
    }
}
