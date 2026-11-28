//! # Finite Field Arithmetic
//!
//! This module provides numerical implementations for arithmetic in finite fields.
//! It includes support for prime fields GF(p) (with division-free Montgomery and
//! Barrett engines, Montgomery's trick batch inversion, 4-way unrolled and parallel
//! slice kernels, and cached Number-Theoretic Transform plans) and arithmetic for
//! GF(2^8) in the polynomial basis modulo `0x11d`. Use [`Gf256Multiplier`] to reuse
//! a fixed coefficient across byte slices.

use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;
use std::sync::Arc;

use num_traits::One;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

/// Represents an element in a prime field GF(p), where p is the modulus.
///
/// The value is stored as a `u64`, and all arithmetic operations are performed
/// modulo the specified `modulus`.
///
/// # Modulus-agnostic elements
/// A `modulus` of `0` marks a *modulus-agnostic* element: the integer
/// `value` (a non-negative multiple of the field's `1`) that has not yet
/// been bound to a field. [`Zero::zero`] and [`One::one`] return such
/// elements, so they can be created without knowing the field. In a binary
/// operation the result adopts the modulus of the concrete operand, and the
/// agnostic value is reduced modulo it; `x + zero() == x` and
/// `x * one() == x` hold for every concrete `x`. Between two agnostic
/// elements the operations act on non-negative integers (addition and
/// multiplication saturate at `u64::MAX`, subtraction saturates at `0`);
/// negating a non-zero agnostic element is not representable and leaves it
/// unchanged. When both operands are concrete with different moduli the
/// left operand's modulus is used.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrimeFieldElement {
    /// The value of the field element.
    pub value: u64,
    /// The modulus (prime number) of the field.
    pub modulus: u64,
}

impl PrimeFieldElement {
    /// Creates a new `PrimeFieldElement`.
    ///
    /// The initial value is reduced modulo the `modulus`.
    ///
    /// # Arguments
    /// * `value` - The initial value of the element.
    /// * `modulus` - The prime modulus of the field.
    #[must_use]
    pub const fn new(
        value: u64,
        modulus: u64,
    ) -> Self {
        Self {
            value: if modulus == 0 {
                value
            } else {
                value % modulus
            },
            modulus,
        }
    }

    /// Returns `true` for a modulus-agnostic element (see the type docs).
    #[must_use]
    pub const fn is_agnostic(&self) -> bool {
        self.modulus == 0
    }

    /// Modulus the result of a binary operation adopts: the left operand's
    /// unless it is agnostic.
    const fn result_modulus(
        &self,
        rhs: &Self,
    ) -> u64 {
        if self.modulus == 0 {
            rhs.modulus
        } else {
            self.modulus
        }
    }

    /// The value reduced modulo `m` (unchanged for `m == 0`).
    const fn reduced(
        &self,
        m: u64,
    ) -> u64 {
        if m == 0 {
            self.value
        } else {
            self.value % m
        }
    }

    /// Computes the multiplicative inverse of the element.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod modulus)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(PrimeFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        if self.modulus == 0 {
            // Only 1 is invertible among agnostic integers.
            return (self.value == 1).then_some(*self);
        }

        mod_inv_u64(self.value, self.modulus).map(|inv| Self::new(inv, self.modulus))
    }

    /// Computes (self^exp) modulo modulus.
    #[must_use]
    pub fn pow(
        &self,
        mut exp: u64,
    ) -> Self {
        if self.modulus == 0 {
            // Agnostic integers: saturating power.
            return Self::new(
                self.value
                    .saturating_pow(u32::try_from(exp).unwrap_or(u32::MAX)),
                0,
            );
        }

        if let Some(mont) = MontgomeryPrimeField::new(self.modulus) {
            return Self::new(mont.pow(self.value, exp), self.modulus);
        }

        let mut res = 1u128 % u128::from(self.modulus);
        let mut base = u128::from(self.value % self.modulus);
        let m = u128::from(self.modulus);

        while exp > 0 {
            if exp & 1 == 1 {
                res = (res * base) % m;
            }
            base = (base * base) % m;
            exp >>= 1;
        }

        Self::new(res as u64, self.modulus)
    }
}

impl Zero for PrimeFieldElement {
    fn zero() -> Self {
        // Modulus-agnostic zero: adopts the modulus of the other operand.
        Self { value: 0, modulus: 0 }
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }
}

impl One for PrimeFieldElement {
    fn one() -> Self {
        // Modulus-agnostic one: adopts the modulus of the other operand.
        Self { value: 1, modulus: 0 }
    }
}

impl Neg for PrimeFieldElement {
    type Output = Self;

    fn neg(self) -> Self {
        if self.value == 0 || self.modulus == 0 {
            self
        } else {
            Self::new(
                self.modulus.saturating_sub(self.value % self.modulus),
                self.modulus,
            )
        }
    }
}

// The reduction contexts and the inverse primitive live in `kernels::modular`.
// They are re-exported here because every existing path names them through
// this module.
pub use crate::kernels::modular::Barrett as BarrettMod;
pub use crate::kernels::modular::Montgomery as MontgomeryPrimeField;
pub(crate) use crate::kernels::modular::add_mod_raw;
pub use crate::kernels::modular::mod_inv_u64;
pub(crate) use crate::kernels::modular::sub_mod_raw;

/// Inverts a slice of field elements modulo `modulus` in place into `out` using
/// Montgomery's prefix-product trick (1 extended-GCD inversion + `3(n - 1)`
/// modular multiplications).
///
/// Leases temporary prefix storage from [`crate::kernels::reuse::global_pool`].
/// Returns `false` (leaving `out` empty) if `modulus < 2` or any element is not
/// invertible modulo `modulus`.
#[must_use]
#[expect(
    clippy::indexing_slicing,
    reason = "prefix and output buffers are sized to values.len() and loop indices stay in 0..n"
)]
pub fn batch_inverse_into(
    values: &[u64],
    modulus: u64,
    out: &mut Vec<u64>,
) -> bool {
    out.clear();
    let Some(barrett) = BarrettMod::new(modulus) else {
        return false;
    };
    let n = values.len();
    if n == 0 {
        return true;
    }
    out.resize(n, 0);
    let mut prefix = crate::kernels::reuse::global_pool().scratch_i64(n);
    prefix.resize(n, 0);

    let mut acc = 1u64;
    for i in 0..n {
        let v = barrett.reduce_u64(values[i]);
        if v == 0 {
            out.clear();
            return false;
        }
        prefix[i] = acc as i64;
        acc = barrett.mul(acc, v);
    }

    let Some(mut inv_acc) = barrett.inv(acc) else {
        out.clear();
        return false;
    };

    for i in (0..n).rev() {
        let v = barrett.reduce_u64(values[i]);
        out[i] = barrett.mul(inv_acc, prefix[i] as u64);
        inv_acc = barrett.mul(inv_acc, v);
    }
    true
}

/// Computes the modular inverses of all elements in `values` modulo `modulus`
/// using Montgomery's batch-inversion trick.
///
/// Returns `None` if `modulus < 2` or any element is not coprime to `modulus`.
#[must_use]
pub fn batch_inverse_mod(
    values: &[u64],
    modulus: u64,
) -> Option<Vec<u64>> {
    let mut out = Vec::with_capacity(values.len());
    batch_inverse_into(values, modulus, &mut out).then_some(out)
}

/// Computes the modular inverses of all invertible elements in `values` modulo
/// `modulus` in a single batch (mapping `0` and any non-invertible elements to `0`).
///
/// When `modulus` is prime, every non-zero element is inverted via a single
/// extended-GCD call plus `3(k - 1)` multiplications while zeros pass through as `0`.
#[must_use]
#[expect(
    clippy::indexing_slicing,
    reason = "prefix buffer and output vector have length values.len() and indices are in 0..n"
)]
pub fn batch_inverse_skip_zeros(
    values: &[u64],
    modulus: u64,
) -> Vec<u64> {
    let n = values.len();
    let mut out = vec![0u64; n];
    let Some(barrett) = BarrettMod::new(modulus) else {
        return out;
    };
    if n == 0 {
        return out;
    }

    let mut prefix = crate::kernels::reuse::global_pool().scratch_i64(n);
    prefix.resize(n, 0);

    let mut acc = 1u64;
    let mut any_nonzero = false;
    for i in 0..n {
        let v = barrett.reduce_u64(values[i]);
        if v == 0 {
            prefix[i] = 0;
        } else {
            prefix[i] = acc as i64;
            acc = barrett.mul(acc, v);
            any_nonzero = true;
        }
    }

    if !any_nonzero {
        return out;
    }

    if let Some(mut inv_acc) = barrett.inv(acc) {
        for i in (0..n).rev() {
            let v = barrett.reduce_u64(values[i]);
            if v != 0 {
                out[i] = barrett.mul(inv_acc, prefix[i] as u64);
                inv_acc = barrett.mul(inv_acc, v);
            }
        }
    } else {
        // Composite modulus where some non-zero elements share a factor with modulus:
        // invert each element individually so invertible entries still succeed.
        for i in 0..n {
            out[i] = barrett.inv(values[i]).unwrap_or(0);
        }
    }
    out
}

const SLICE_PAR_MIN_CHUNK: usize = 4_096;

#[expect(
    clippy::indexing_slicing,
    reason = "chunks_exact(4) guarantees 4 elements per chunk"
)]
fn slice_add_mod_seq(
    a: &[u64],
    b: &[u64],
    barrett: BarrettMod,
    out: &mut [u64],
) {
    let n = a.len().min(b.len()).min(out.len());
    let (a, b, out) = (&a[..n], &b[..n], &mut out[..n]);
    let (a_chunks, a_tail) = a.as_chunks::<4>();
    let (b_chunks, b_tail) = b.as_chunks::<4>();
    let (out_chunks, out_tail) = out.as_chunks_mut::<4>();
    for ((xa, xb), yo) in a_chunks.iter().zip(b_chunks).zip(out_chunks) {
        yo[0] = barrett.add(xa[0], xb[0]);
        yo[1] = barrett.add(xa[1], xb[1]);
        yo[2] = barrett.add(xa[2], xb[2]);
        yo[3] = barrett.add(xa[3], xb[3]);
    }
    for ((&xa, &xb), yo) in a_tail.iter().zip(b_tail).zip(out_tail) {
        *yo = barrett.add(xa, xb);
    }
}

/// Elementwise modular addition `out[i] = (a[i] + b[i]) mod modulus` over the
/// common prefix of `a`, `b`, and `out`.
///
/// Uses 4-way loop unrolling and splits across Rayon workers when
/// [`crate::execution::should_split_work`] permits.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunk offsets are bounded by n = a.len().min(b.len()).min(out.len())"
)]
pub fn slice_add_mod(
    a: &[u64],
    b: &[u64],
    modulus: u64,
    out: &mut [u64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    let Some(barrett) = BarrettMod::new(modulus) else {
        out[..n].fill(0);
        return;
    };
    let (a, b, out) = (&a[..n], &b[..n], &mut out[..n]);
    if crate::execution::should_split_work(n, 2) {
        let chunk = crate::execution::choose_chunk_size(n, SLICE_PAR_MIN_CHUNK);
        crate::kernels::par::for_chunks_mut(out, chunk, 2, |idx, out_chunk| {
            let start = idx * chunk;
            let end = start + out_chunk.len();
            slice_add_mod_seq(&a[start..end], &b[start..end], barrett, out_chunk);
        });
    } else {
        slice_add_mod_seq(a, b, barrett, out);
    }
}

#[expect(
    clippy::indexing_slicing,
    reason = "chunks_exact(4) guarantees 4 elements per chunk"
)]
fn slice_sub_mod_seq(
    a: &[u64],
    b: &[u64],
    barrett: BarrettMod,
    out: &mut [u64],
) {
    let n = a.len().min(b.len()).min(out.len());
    let (a, b, out) = (&a[..n], &b[..n], &mut out[..n]);
    let (a_chunks, a_tail) = a.as_chunks::<4>();
    let (b_chunks, b_tail) = b.as_chunks::<4>();
    let (out_chunks, out_tail) = out.as_chunks_mut::<4>();
    for ((xa, xb), yo) in a_chunks.iter().zip(b_chunks).zip(out_chunks) {
        yo[0] = barrett.sub(xa[0], xb[0]);
        yo[1] = barrett.sub(xa[1], xb[1]);
        yo[2] = barrett.sub(xa[2], xb[2]);
        yo[3] = barrett.sub(xa[3], xb[3]);
    }
    for ((&xa, &xb), yo) in a_tail.iter().zip(b_tail).zip(out_tail) {
        *yo = barrett.sub(xa, xb);
    }
}

/// Elementwise modular subtraction `out[i] = (a[i] - b[i]) mod modulus` over the
/// common prefix of `a`, `b`, and `out`.
///
/// Uses 4-way loop unrolling and splits across Rayon workers when
/// [`crate::execution::should_split_work`] permits.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunk offsets are bounded by n = a.len().min(b.len()).min(out.len())"
)]
pub fn slice_sub_mod(
    a: &[u64],
    b: &[u64],
    modulus: u64,
    out: &mut [u64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    let Some(barrett) = BarrettMod::new(modulus) else {
        out[..n].fill(0);
        return;
    };
    let (a, b, out) = (&a[..n], &b[..n], &mut out[..n]);
    if crate::execution::should_split_work(n, 2) {
        let chunk = crate::execution::choose_chunk_size(n, SLICE_PAR_MIN_CHUNK);
        crate::kernels::par::for_chunks_mut(out, chunk, 2, |idx, out_chunk| {
            let start = idx * chunk;
            let end = start + out_chunk.len();
            slice_sub_mod_seq(&a[start..end], &b[start..end], barrett, out_chunk);
        });
    } else {
        slice_sub_mod_seq(a, b, barrett, out);
    }
}

#[expect(
    clippy::indexing_slicing,
    reason = "chunks_exact(4) guarantees 4 elements per chunk"
)]
fn slice_scale_mod_seq(
    a: &[u64],
    scalar: u64,
    barrett: BarrettMod,
    out: &mut [u64],
) {
    let n = a.len().min(out.len());
    let (a, out) = (&a[..n], &mut out[..n]);
    let s = barrett.reduce_u64(scalar);
    if s == 0 {
        out.fill(0);
        return;
    }
    if s == 1 {
        let (a_chunks, a_tail) = a.as_chunks::<4>();
        let (out_chunks, out_tail) = out.as_chunks_mut::<4>();
        for (xa, yo) in a_chunks.iter().zip(out_chunks) {
            yo[0] = barrett.reduce_u64(xa[0]);
            yo[1] = barrett.reduce_u64(xa[1]);
            yo[2] = barrett.reduce_u64(xa[2]);
            yo[3] = barrett.reduce_u64(xa[3]);
        }
        for (&xa, yo) in a_tail.iter().zip(out_tail) {
            *yo = barrett.reduce_u64(xa);
        }
        return;
    }
    let (a_chunks, a_tail) = a.as_chunks::<4>();
    let (out_chunks, out_tail) = out.as_chunks_mut::<4>();
    for (xa, yo) in a_chunks.iter().zip(out_chunks) {
        yo[0] = barrett.mul(xa[0], s);
        yo[1] = barrett.mul(xa[1], s);
        yo[2] = barrett.mul(xa[2], s);
        yo[3] = barrett.mul(xa[3], s);
    }
    for (&xa, yo) in a_tail.iter().zip(out_tail) {
        *yo = barrett.mul(xa, s);
    }
}

/// Elementwise modular scaling `out[i] = (scalar * a[i]) mod modulus` over the
/// common prefix of `a` and `out`.
///
/// Uses 4-way loop unrolling and splits across Rayon workers when
/// [`crate::execution::should_split_work`] permits.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunk offsets are bounded by n = a.len().min(out.len())"
)]
pub fn slice_scale_mod(
    a: &[u64],
    scalar: u64,
    modulus: u64,
    out: &mut [u64],
) {
    let n = a.len().min(out.len());
    if n == 0 {
        return;
    }
    let Some(barrett) = BarrettMod::new(modulus) else {
        out[..n].fill(0);
        return;
    };
    let (a, out) = (&a[..n], &mut out[..n]);
    if crate::execution::should_split_work(n, 4) {
        let chunk = crate::execution::choose_chunk_size(n, SLICE_PAR_MIN_CHUNK);
        crate::kernels::par::for_chunks_mut(out, chunk, 2, |idx, out_chunk| {
            let start = idx * chunk;
            let end = start + out_chunk.len();
            slice_scale_mod_seq(&a[start..end], scalar, barrett, out_chunk);
        });
    } else {
        slice_scale_mod_seq(a, scalar, barrett, out);
    }
}

#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunks_exact(4) guarantees 4 elements per chunk and 128-bit fused products stay within u128"
)]
fn slice_axpy_mod_seq(
    alpha_red: u64,
    x: &[u64],
    y: &mut [u64],
    barrett: BarrettMod,
) {
    let n = x.len().min(y.len());
    let (x, y) = (&x[..n], &mut y[..n]);
    if alpha_red == 0 {
        for yi in y {
            *yi = barrett.reduce_u64(*yi);
        }
        return;
    }
    if alpha_red == 1 {
        let (x_chunks, x_tail) = x.as_chunks::<4>();
        let (y_chunks, y_tail) = y.as_chunks_mut::<4>();
        for (xx, yy) in x_chunks.iter().zip(y_chunks) {
            yy[0] = barrett.add(yy[0], xx[0]);
            yy[1] = barrett.add(yy[1], xx[1]);
            yy[2] = barrett.add(yy[2], xx[2]);
            yy[3] = barrett.add(yy[3], xx[3]);
        }
        for (&xx, yy) in x_tail.iter().zip(y_tail) {
            *yy = barrett.add(*yy, xx);
        }
        return;
    }
    let a128 = u128::from(alpha_red);
    let (x_chunks, x_tail) = x.as_chunks::<4>();
    let (y_chunks, y_tail) = y.as_chunks_mut::<4>();
    for (xx, yy) in x_chunks.iter().zip(y_chunks) {
        let y0 = u128::from(barrett.reduce_u64(yy[0]));
        let y1 = u128::from(barrett.reduce_u64(yy[1]));
        let y2 = u128::from(barrett.reduce_u64(yy[2]));
        let y3 = u128::from(barrett.reduce_u64(yy[3]));
        let x0 = u128::from(barrett.reduce_u64(xx[0]));
        let x1 = u128::from(barrett.reduce_u64(xx[1]));
        let x2 = u128::from(barrett.reduce_u64(xx[2]));
        let x3 = u128::from(barrett.reduce_u64(xx[3]));
        yy[0] = barrett.reduce_u128(y0 + a128 * x0);
        yy[1] = barrett.reduce_u128(y1 + a128 * x1);
        yy[2] = barrett.reduce_u128(y2 + a128 * x2);
        yy[3] = barrett.reduce_u128(y3 + a128 * x3);
    }
    for (&xx, yy) in x_tail.iter().zip(y_tail) {
        let yv = u128::from(barrett.reduce_u64(*yy));
        let xv = u128::from(barrett.reduce_u64(xx));
        *yy = barrett.reduce_u128(yv + a128 * xv);
    }
}

/// Modular AXPY `y[i] = (y[i] + alpha * x[i]) mod modulus` over the common
/// prefix of `x` and `y`.
///
/// Forms the core inner loop of modular Gaussian elimination, F4 Macaulay
/// matrix reduction, and Euclidean polynomial division. Uses 4-way unrolled
/// fused `u128` multiply-add with [`BarrettMod`] reduction and splits across
/// Rayon workers on large vectors.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunk offsets are bounded by n = x.len().min(y.len())"
)]
pub fn slice_axpy_mod(
    alpha: u64,
    x: &[u64],
    y: &mut [u64],
    modulus: u64,
) {
    let n = x.len().min(y.len());
    if n == 0 {
        return;
    }
    let Some(barrett) = BarrettMod::new(modulus) else {
        y[..n].fill(0);
        return;
    };
    let (x, y) = (&x[..n], &mut y[..n]);
    let alpha_red = barrett.reduce_u64(alpha);
    if crate::execution::should_split_work(n, 4) {
        let chunk = crate::execution::choose_chunk_size(n, SLICE_PAR_MIN_CHUNK);
        crate::kernels::par::for_chunks_mut(y, chunk, 2, |idx, y_chunk| {
            let start = idx * chunk;
            let end = start + y_chunk.len();
            slice_axpy_mod_seq(alpha_red, &x[start..end], y_chunk, barrett);
        });
    } else {
        slice_axpy_mod_seq(alpha_red, x, y, barrett);
    }
}

#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunks_exact(4) guarantees 4 elements and 4 products of <2^62 integers fit in u128 without overflow"
)]
fn dot_mod_block(
    a: &[u64],
    b: &[u64],
    barrett: BarrettMod,
) -> u64 {
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let can_fuse_four = barrett.modulus() < (1u64 << 62);
    let mut acc = 0u64;
    let (a_chunks, a_tail) = a.as_chunks::<4>();
    let (b_chunks, b_tail) = b.as_chunks::<4>();
    if can_fuse_four {
        for (xa, xb) in a_chunks.iter().zip(b_chunks) {
            let a0 = u128::from(barrett.reduce_u64(xa[0]));
            let a1 = u128::from(barrett.reduce_u64(xa[1]));
            let a2 = u128::from(barrett.reduce_u64(xa[2]));
            let a3 = u128::from(barrett.reduce_u64(xa[3]));
            let b0 = u128::from(barrett.reduce_u64(xb[0]));
            let b1 = u128::from(barrett.reduce_u64(xb[1]));
            let b2 = u128::from(barrett.reduce_u64(xb[2]));
            let b3 = u128::from(barrett.reduce_u64(xb[3]));
            let sum4 = (a0 * b0 + a1 * b1) + (a2 * b2 + a3 * b3);
            acc = barrett.add(acc, barrett.reduce_u128(sum4));
        }
    } else {
        for (xa, xb) in a_chunks.iter().zip(b_chunks) {
            let p0 = barrett.mul(xa[0], xb[0]);
            let p1 = barrett.mul(xa[1], xb[1]);
            let p2 = barrett.mul(xa[2], xb[2]);
            let p3 = barrett.mul(xa[3], xb[3]);
            acc = barrett.add(acc, barrett.add(barrett.add(p0, p1), barrett.add(p2, p3)));
        }
    }
    for (&xa, &xb) in a_tail.iter().zip(b_tail) {
        acc = barrett.add(acc, barrett.mul(xa, xb));
    }
    acc
}

/// Deterministic modular dot product `sum_i (a[i] * b[i]) mod modulus` over
/// the common prefix of `a` and `b`.
///
/// Accumulates 4-way unrolled `u128` products in fixed-size blocks before
/// reducing modulo `modulus`, and folds block residues deterministically in
/// index order when split across Rayon workers.
#[must_use]
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "block bounds lo..hi are clamped to n = a.len().min(b.len())"
)]
pub fn dot_mod(
    a: &[u64],
    b: &[u64],
    modulus: u64,
) -> u64 {
    const BLOCK: usize = 4_096;

    let n = a.len().min(b.len());
    if n == 0 {
        return 0;
    }
    let Some(barrett) = BarrettMod::new(modulus) else {
        return 0;
    };
    let (a, b) = (&a[..n], &b[..n]);
    if crate::execution::should_split_work(n, 4) {
        let nb = n.div_ceil(BLOCK);
        let block_sums = crate::kernels::par::map_indexed(nb, 2, |blk| {
            let lo = blk * BLOCK;
            let hi = (lo + BLOCK).min(n);
            dot_mod_block(&a[lo..hi], &b[lo..hi], barrett)
        });
        block_sums
            .into_iter()
            .fold(0u64, |acc, s| barrett.add(acc, s))
    } else {
        dot_mod_block(a, b, barrett)
    }
}

/// Checks if `p` is prime, by the crate-wide deterministic Miller–Rabin test.
///
/// The implementation is
/// [`crate::kernels::number_theory::is_prime_miller_rabin`]; the name is kept
/// because callers in this module were written against it.
pub use crate::kernels::number_theory::is_prime_miller_rabin as is_prime_u64;

/// Finds the smallest positive primitive root modulo a prime `modulus`.
///
/// The root and Montgomery context for an odd prime are cached by modulus.
/// Returns `None` if the modulus is composite, zero or one.
#[must_use]
pub fn primitive_root_prime(modulus: u64) -> Option<u64> {
    if modulus == 2 {
        return Some(1);
    }
    super::number_theory::odd_prime_generator(modulus).map(|group| group.root)
}

/// Reusable, immutable Number-Theoretic Transform (NTT) plan for power-of-two
/// transform length `len` over an NTT-friendly prime `modulus` (`len | (modulus - 1)`).
///
/// All twiddle factors (`roots_mont`, `inv_roots_mont`) and the length inverse
/// `n_inv_mont` are stored in Montgomery representation for [`MontgomeryPrimeField`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NttPrimePlan {
    /// The NTT-friendly prime modulus `p`.
    pub modulus: u64,
    /// Precomputed Montgomery arithmetic engine for `modulus`.
    pub mont: MontgomeryPrimeField,
    /// Power-of-two transform length `n`.
    pub len: usize,
    /// Base-2 logarithm of `len` (`len == 1 << log_len`).
    pub log_len: u32,
    /// Forward twiddle powers `[1, omega, omega^2, ..., omega^(n/2 - 1)]` in Montgomery form.
    pub roots_mont: Vec<u64>,
    /// Inverse twiddle powers `[1, omega^(-1), omega^(-2), ..., omega^(-(n/2 - 1))]` in Montgomery form.
    pub inv_roots_mont: Vec<u64>,
    /// Multiplicative inverse of `len` modulo `modulus` in Montgomery form.
    pub n_inv_mont: u64,
    /// Bit-reversal permutation indices of length `len`.
    pub bit_rev: Vec<usize>,
}

#[allow(
    clippy::single_call_fn,
    reason = "isolates NTT plan construction from the cache lookup"
)]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "len is a non-zero power of two dividing modulus - 1 and bit_rev has length len"
)]
fn build_ntt_plan(
    modulus: u64,
    len: usize,
) -> Option<NttPrimePlan> {
    if len == 0 || !len.is_power_of_two() || modulus <= 2 {
        return None;
    }
    let len_u64 = u64::try_from(len).ok()?;
    if !(modulus - 1).is_multiple_of(len_u64) {
        return None;
    }
    let group = super::number_theory::odd_prime_generator(modulus)?;
    let mont = group.mont;
    let log_len = len.trailing_zeros();

    let g_mont = mont.to_mont(group.root);
    let omega_mont = mont.pow_mont(g_mont, (modulus - 1) / len_u64);
    let omega_inv = mont.inv(mont.from_mont(omega_mont))?;
    let omega_inv_mont = mont.to_mont(omega_inv);

    let half = (len / 2).max(1);
    let mut roots_mont = Vec::with_capacity(half);
    let mut inv_roots_mont = Vec::with_capacity(half);
    let mut cur = mont.one_mont();
    let mut cur_inv = mont.one_mont();
    for _ in 0..half {
        roots_mont.push(cur);
        inv_roots_mont.push(cur_inv);
        cur = mont.mul_mont(cur, omega_mont);
        cur_inv = mont.mul_mont(cur_inv, omega_inv_mont);
    }

    let n_inv = mont.inv(len_u64)?;
    let n_inv_mont = mont.to_mont(n_inv);

    let mut bit_rev = vec![0usize; len];
    if log_len > 0 {
        for i in 0..len {
            bit_rev[i] = i.reverse_bits() >> (usize::BITS - log_len);
        }
    }

    Some(NttPrimePlan {
        modulus,
        mont,
        len,
        log_len,
        roots_mont,
        inv_roots_mont,
        n_inv_mont,
        bit_rev,
    })
}

/// Looks up or constructs a cached [`NttPrimePlan`] for `(modulus, len)` in
/// [`crate::kernels::reuse::global_pool`]'s named plan cache `"finite_field_ntt_plan"`.
///
/// Returns `None` if `modulus` is not an odd prime, `len` is not a power of two,
/// or `len` does not divide `modulus - 1`.
#[must_use]
pub fn ntt_plan(
    modulus: u64,
    len: usize,
) -> Option<Arc<NttPrimePlan>> {
    let cache = crate::kernels::reuse::global_pool()
        .named_plan_cache::<(u64, usize), NttPrimePlan>("finite_field_ntt_plan", 64);
    cache
        .try_get_or_insert_with((modulus, len), || {
            build_ntt_plan(modulus, len).map(Arc::new).ok_or(())
        })
        .ok()
}

/// Executes an in-place Number-Theoretic Transform on Montgomery-domain coefficients
/// `values_mont` (length `plan.len`).
///
/// When `inverse` is `true`, computes the inverse NTT including final scaling by
/// `n^(-1)` in Montgomery form.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "butterfly indices and twiddle indices stay strictly within plan.len and plan.len / 2"
)]
pub fn ntt_in_place_mont(
    values_mont: &mut [u64],
    plan: &NttPrimePlan,
    inverse: bool,
) {
    let n = plan.len;
    if values_mont.len() < n || n <= 1 {
        return;
    }
    let values = &mut values_mont[..n];
    for i in 0..n {
        let j = plan.bit_rev[i];
        if i < j {
            values.swap(i, j);
        }
    }
    let mont = plan.mont;
    let roots = if inverse {
        &plan.inv_roots_mont
    } else {
        &plan.roots_mont
    };
    let mut half_step = 1usize;
    while half_step < n {
        let step = half_step * 2;
        let stride = (n / 2) / half_step;
        let num_blocks = n / step;
        if n >= 16_384 && num_blocks >= 4 {
            crate::kernels::par::for_chunks_mut(values, step, 4, |_, block| {
                for j in 0..half_step {
                    let w = roots[j * stride];
                    let u = block[j];
                    let v = mont.mul_mont(block[j + half_step], w);
                    block[j] = mont.add_mont(u, v);
                    block[j + half_step] = mont.sub_mont(u, v);
                }
            });
        } else {
            let mut base = 0usize;
            while base < n {
                for j in 0..half_step {
                    let w = roots[j * stride];
                    let u = values[base + j];
                    let v = mont.mul_mont(values[base + j + half_step], w);
                    values[base + j] = mont.add_mont(u, v);
                    values[base + j + half_step] = mont.sub_mont(u, v);
                }
                base += step;
            }
        }
        half_step = step;
    }
    if inverse {
        let n_inv_mont = plan.n_inv_mont;
        if n >= 16_384 {
            crate::kernels::par::for_chunks_mut(values, 4096, 4, |_, chunk| {
                for slot in chunk {
                    *slot = mont.mul_mont(*slot, n_inv_mont);
                }
            });
        } else {
            for slot in values {
                *slot = mont.mul_mont(*slot, n_inv_mont);
            }
        }
    }
}

/// Multiplies two polynomials `a` and `b` modulo an NTT-friendly prime `modulus`
/// using cached [`NttPrimePlan`] tables from [`crate::kernels::reuse::global_pool`].
///
/// Inputs are in the same coefficient order (ascending or descending) and
/// the product follows that order. Its length is `a.len() + b.len() - 1`,
/// or zero if either input is empty. A nonempty product requires an odd
/// prime modulus and the next power-of-two transform length to divide
/// `modulus - 1`. Returns `None` if these conditions fail.
#[must_use]
#[expect(
    clippy::indexing_slicing,
    reason = "buffers fa and fb have length n >= out_len"
)]
pub fn poly_mul_ntt_mod(
    a: &[u64],
    b: &[u64],
    modulus: u64,
) -> Option<Vec<u64>> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let out_len = a.len().checked_add(b.len())?.checked_sub(1)?;
    let n = out_len.checked_next_power_of_two()?;
    let plan = ntt_plan(modulus, n)?;
    let mont = plan.mont;

    let mut fa = vec![0u64; n];
    let mut fb = vec![0u64; n];
    for (i, &c) in a.iter().enumerate() {
        fa[i] = mont.to_mont(c);
    }
    for (i, &c) in b.iter().enumerate() {
        fb[i] = mont.to_mont(c);
    }

    if crate::execution::should_split(n, 4096) {
        rayon::join(
            || ntt_in_place_mont(&mut fa, &plan, false),
            || ntt_in_place_mont(&mut fb, &plan, false),
        );
    } else {
        ntt_in_place_mont(&mut fa, &plan, false);
        ntt_in_place_mont(&mut fb, &plan, false);
    }

    super::par::for_each_mut(&mut fa, 16_384, |i, coefficient| {
        *coefficient = mont.mul_mont(*coefficient, fb[i]);
    });
    drop(fb);

    ntt_in_place_mont(&mut fa, &plan, true);
    fa.truncate(out_len);
    super::par::for_each_mut(&mut fa, 16_384, |_, coefficient| {
        *coefficient = mont.from_mont(*coefficient);
    });
    Some(fa)
}

impl Add for PrimeFieldElement {
    type Output = Self;

    /// Performs addition in the prime field.
    fn add(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        if m == 0 {
            let sum = u128::from(self.value) + u128::from(rhs.value);
            return Self::new(u64::try_from(sum).unwrap_or(u64::MAX), 0);
        }

        Self::new(add_mod_raw(self.value, rhs.value, m), m)
    }
}

impl Sub for PrimeFieldElement {
    type Output = Self;

    /// Performs subtraction in the prime field.
    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        if m == 0 {
            return Self::new(self.value.saturating_sub(rhs.value), 0);
        }

        Self::new(sub_mod_raw(self.value, rhs.value, m), m)
    }
}

impl Mul for PrimeFieldElement {
    type Output = Self;

    /// Performs multiplication in the prime field.
    ///
    /// Uses `u128` for the intermediate multiplication to prevent overflow
    /// before the modulo operation.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "m != 0 and product of two u64 values fits in u128"
    )]
    fn mul(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        let prod = u128::from(self.reduced(m)) * u128::from(rhs.reduced(m));

        if m == 0 {
            return Self::new(u64::try_from(prod).unwrap_or(u64::MAX), 0);
        }

        Self::new((prod % u128::from(m)) as u64, m)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for PrimeFieldElement {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        let rhs = Self::new(rhs.value, m);

        let Some(inv_rhs) = rhs.inverse() else {
            return Self::new(0, m);
        };

        self * inv_rhs
    }
}

const GF256_GENERATOR_POLY: u16 = 0x11d;
const GF256_SIZE: usize = 256;
const GF256_ORDER: usize = GF256_SIZE - 1;

struct Gf256Tables {
    log: [u8; GF256_SIZE],
    // Repeat the cycle to index by either a sum or an offset difference of logs.
    // Two spare entries give an array large enough for any two byte indices.
    exp: [u8; 2 * GF256_SIZE],
}

#[expect(
    clippy::indexing_slicing,
    reason = "the loop bounds the exponent index and polynomial reduction leaves a byte"
)]
static GF256_TABLES: Gf256Tables = {
    let mut tables = Gf256Tables {
        log: [0; GF256_SIZE],
        exp: [0; 2 * GF256_SIZE],
    };
    let mut i = 0;
    let mut x: u16 = 1;
    while i < tables.exp.len() {
        tables.exp[i] = x as u8;
        if i < GF256_ORDER {
            tables.log[x as usize] = i as u8;
        }
        x <<= 1;
        if x >= GF256_SIZE as u16 {
            x ^= GF256_GENERATOR_POLY;
        }
        i += 1;
    }
    tables
};

/// Performs addition in GF(2^8).
///
/// In GF(2^8), addition is equivalent to a bitwise XOR operation.
#[inline]
#[must_use]
pub const fn gf256_add(
    a: u8,
    b: u8,
) -> u8 {
    a ^ b
}

/// Performs multiplication in GF(2^8) using lookup tables.
///
/// The compile-time logarithm and exponentiation tables use the `0x11d` polynomial basis.
#[inline]
#[must_use]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "byte log indices fit their table and the sum of two bytes is at most 510"
)]
pub const fn gf256_mul(
    a: u8,
    b: u8,
) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }

    let log_a = GF256_TABLES.log[a as usize] as usize;
    let log_b = GF256_TABLES.log[b as usize] as usize;
    GF256_TABLES.exp[log_a + log_b]
}

/// Computes the multiplicative inverse in GF(2^8).
///
/// # Errors
/// * Returns an error if `a` is 0, as 0 has no multiplicative inverse.
#[inline]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "subtracting a byte log from 255 gives an exponent in 0..=255"
)]
pub fn gf256_inv(a: u8) -> Result<u8, String> {
    if a == 0 {
        return Err("Cannot invert 0".to_string());
    }

    Ok(GF256_TABLES.exp[GF256_ORDER - usize::from(GF256_TABLES.log[usize::from(a)])])
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the inverse of the divisor.
///
/// # Errors
/// * Returns an error if the divisor `b` is 0.
#[inline]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "adding 255 to the log bounds the nonnegative difference by 510"
)]
pub fn gf256_div(
    a: u8,
    b: u8,
) -> Result<u8, String> {
    if b == 0 {
        return Err("Division by zero".to_string());
    }

    if a == 0 {
        return Ok(0);
    }

    let log_a = usize::from(GF256_TABLES.log[usize::from(a)]);
    let log_b = usize::from(GF256_TABLES.log[usize::from(b)]);
    Ok(GF256_TABLES.exp[log_a + GF256_ORDER - log_b])
}

/// Performs exponentiation in GF(2^8).
///
/// Exponent zero returns `1`, including for a zero base.
#[must_use]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "the reduced exponent product fits u16 and the residue is at most 254"
)]
pub fn gf256_pow(
    a: u8,
    exp: u64,
) -> u8 {
    if exp == 0 {
        return 1;
    }

    if a == 0 {
        return 0;
    }

    let log_a = u16::from(GF256_TABLES.log[usize::from(a)]);
    // Nonzero elements have order dividing 255, so reducing first keeps the product in u16.
    let reduced_exp = (exp % GF256_ORDER as u64) as u16;
    let new_log = (log_a * reduced_exp) % GF256_ORDER as u16;
    GF256_TABLES.exp[usize::from(new_log)]
}

/// A reusable multiplier for a fixed coefficient in GF(2^8).
///
/// It owns a 256-byte lookup table and needs no heap allocation. Reuse it for operations
/// on several slices with the same coefficient, such as erasure-coding stripes.
///
/// ```
/// use rssn::kernels::finite_field::Gf256Multiplier;
///
/// let multiplier = Gf256Multiplier::new(2);
/// let mut parity = [0x10, 0x20];
/// multiplier.mul_add_assign(&[0x80, 1], &mut parity)?;
/// assert_eq!(parity, [0x0d, 0x22]);
/// # Ok::<(), String>(())
/// ```
#[derive(Clone, Debug)]
pub struct Gf256Multiplier {
    coefficient: u8,
    products: [u8; GF256_SIZE],
}

#[expect(
    clippy::indexing_slicing,
    reason = "the constructor and byte lookups index a 256-entry table"
)]
impl Gf256Multiplier {
    /// Builds a multiplier using the same polynomial basis as [`gf256_mul`].
    #[must_use]
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the table loop increments only while the index is below 256"
    )]
    pub const fn new(coefficient: u8) -> Self {
        let mut products = [0; GF256_SIZE];
        if coefficient > 1 {
            let mut value = 1;
            while value < GF256_SIZE {
                products[value] = gf256_mul(coefficient, value as u8);
                value += 1;
            }
        }
        Self {
            coefficient,
            products,
        }
    }

    /// Multiplies each byte of the slice by the coefficient in place.
    pub fn scale_assign(
        &self,
        values: &mut [u8],
    ) {
        match self.coefficient {
            | 0 => values.fill(0),
            | 1 => {},
            | _ => {
                for value in values {
                    *value = self.products[usize::from(*value)];
                }
            },
        }
    }

    /// Adds the scaled source into the destination using field addition (XOR).
    ///
    /// # Errors
    /// Returns an error and leaves the destination unchanged when the lengths differ.
    pub fn mul_add_assign(
        &self,
        source: &[u8],
        destination: &mut [u8],
    ) -> Result<(), String> {
        if source.len() != destination.len() {
            return Err("GF(256) source and destination lengths must match".into());
        }
        match self.coefficient {
            | 0 => {},
            | 1 => {
                for (&src, dst) in source.iter().zip(destination) {
                    *dst ^= src;
                }
            },
            | _ => {
                for (&src, dst) in source.iter().zip(destination) {
                    *dst ^= self.products[usize::from(src)];
                }
            },
        }
        Ok(())
    }
}

impl AddAssign for PrimeFieldElement {
    fn add_assign(
        &mut self,
        rhs: Self,
    ) {
        *self = *self + rhs;
    }
}

impl SubAssign for PrimeFieldElement {
    fn sub_assign(
        &mut self,
        rhs: Self,
    ) {
        *self = *self - rhs;
    }
}

impl MulAssign for PrimeFieldElement {
    fn mul_assign(
        &mut self,
        rhs: Self,
    ) {
        *self = *self * rhs;
    }
}

impl DivAssign for PrimeFieldElement {
    fn div_assign(
        &mut self,
        rhs: Self,
    ) {
        *self = *self / rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn montgomery_prime_field_matches_exact_modular_arithmetic() {
        assert!(MontgomeryPrimeField::new(0).is_none());
        assert!(MontgomeryPrimeField::new(1).is_none());
        assert!(MontgomeryPrimeField::new(2).is_none());
        assert!(MontgomeryPrimeField::new(100).is_none());

        let primes = [
            3u64,
            5,
            7,
            17,
            257,
            65_537,
            998_244_353,
            2_147_483_647,
            (1u64 << 61) - 1,
            u64::MAX - 58, // 2^64 - 59 is prime
        ];
        for p in primes {
            let mont = MontgomeryPrimeField::new(p).unwrap();
            assert_eq!(mont.modulus(), p);
            assert_eq!(p.wrapping_mul(mont.n_prime()), u64::MAX);
            assert_eq!(mont.from_mont(mont.one_mont()), 1);

            let samples = [
                0u64,
                1,
                2,
                p / 3,
                p / 2,
                p.saturating_sub(2),
                p.saturating_sub(1),
            ];
            for &a in &samples {
                let am = mont.to_mont(a);
                assert_eq!(mont.from_mont(am), a % p);
                assert_eq!(mont.from_mont(mont.neg_mont(am)), (p - (a % p)) % p);

                if a % p != 0 {
                    let inv = mont.inv(a).unwrap();
                    assert_eq!(mont.mul(a, inv), 1);
                    assert_eq!(mont.pow(a, p - 2), inv);
                } else {
                    assert_eq!(mont.inv(a), None);
                }

                for &b in &samples {
                    let bm = mont.to_mont(b);
                    let p128 = u128::from(p);
                    let expected_add = ((u128::from(a) + u128::from(b)) % p128) as u64;
                    let expected_sub = ((u128::from(a) + p128 - u128::from(b % p)) % p128) as u64;
                    let expected_mul = ((u128::from(a) * u128::from(b)) % p128) as u64;

                    assert_eq!(mont.from_mont(mont.add_mont(am, bm)), expected_add);
                    assert_eq!(mont.from_mont(mont.sub_mont(am, bm)), expected_sub);
                    assert_eq!(mont.from_mont(mont.mul_mont(am, bm)), expected_mul);
                    assert_eq!(mont.add(a, b), expected_add);
                    assert_eq!(mont.sub(a, b), expected_sub);
                    assert_eq!(mont.mul(a, b), expected_mul);

                    if b % p != 0 {
                        let q = mont.div(a, b).unwrap();
                        assert_eq!(mont.mul(q, b), a % p);
                    }
                }
            }
        }
    }

    #[test]
    fn barrett_mod_handles_even_odd_and_large_moduli() {
        assert!(BarrettMod::new(0).is_none());
        assert!(BarrettMod::new(1).is_none());

        let moduli = [
            2u64,
            3,
            4,
            7,
            10,
            256,
            65_537,
            998_244_353,
            1_000_000_007,
            (1u64 << 61) - 1,
            (1u64 << 62) - 57,
            (1u64 << 63) - 25,
        ];
        for m in moduli {
            let barrett = BarrettMod::new(m).unwrap();
            assert_eq!(barrett.modulus(), m);
            let samples = [
                0u64,
                1,
                2,
                3,
                m / 2,
                m - 1,
                m,
                m + 1,
                m.saturating_mul(2).saturating_sub(1),
                u64::MAX,
            ];
            for &a in &samples {
                assert_eq!(barrett.reduce_u64(a), a % m);
                for &b in &samples {
                    let m128 = u128::from(m);
                    let a_red = u128::from(a % m);
                    let b_red = u128::from(b % m);
                    assert_eq!(barrett.add(a, b), ((a_red + b_red) % m128) as u64);
                    assert_eq!(barrett.sub(a, b), ((a_red + m128 - b_red) % m128) as u64);
                    assert_eq!(barrett.mul(a, b), ((a_red * b_red) % m128) as u64);
                }
            }
            assert_eq!(barrett.pow(3, 19), {
                let mut acc = 1u128 % u128::from(m);
                for _ in 0..19 {
                    acc = (acc * 3) % u128::from(m);
                }
                acc as u64
            });
        }
    }

    #[test]
    fn batch_inversion_and_skip_zeros_match_elementwise_inverses() {
        let p = 998_244_353u64;
        let barrett = BarrettMod::new(p).unwrap();
        let values: Vec<u64> = (1..=257).map(|i| (i * 123_457) % p).collect();
        let invs = batch_inverse_mod(&values, p).unwrap();
        assert_eq!(invs.len(), values.len());
        for (&v, &inv) in values.iter().zip(&invs) {
            assert_eq!(barrett.mul(v, inv), 1);
            assert_eq!(Some(inv), barrett.inv(v));
        }

        // A zero element causes batch_inverse_mod to return None.
        let mut with_zero = values;
        with_zero[17] = 0;
        with_zero[99] = p;
        assert!(batch_inverse_mod(&with_zero, p).is_none());

        // batch_inverse_skip_zeros preserves 0 at zero positions and inverts all others.
        let skipped = batch_inverse_skip_zeros(&with_zero, p);
        assert_eq!(skipped.len(), with_zero.len());
        for (i, (&v, &inv)) in with_zero.iter().zip(&skipped).enumerate() {
            if i == 17 || i == 99 {
                assert_eq!(inv, 0);
            } else {
                assert_eq!(barrett.mul(v, inv), 1);
            }
        }
    }

    #[test]
    fn slice_ops_and_dot_mod_are_consistent_across_sizes() {
        let p = (1u64 << 61) - 1;
        let barrett = BarrettMod::new(p).unwrap();
        for n in [0usize, 1, 3, 4, 7, 64, 1_025, 20_000] {
            let a: Vec<u64> = (0..n)
                .map(|i| ((i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)) % p)
                .collect();
            let b: Vec<u64> = (0..n)
                .map(|i| ((i as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9) + 7) % p)
                .collect();
            let alpha = 123_456_789_012_345u64 % p;

            let mut out_add = vec![0u64; n];
            let mut out_sub = vec![0u64; n];
            let mut out_scale = vec![0u64; n];
            let mut out_axpy = b.clone();

            slice_add_mod(&a, &b, p, &mut out_add);
            slice_sub_mod(&a, &b, p, &mut out_sub);
            slice_scale_mod(&a, alpha, p, &mut out_scale);
            slice_axpy_mod(alpha, &a, &mut out_axpy, p);

            let mut expected_dot = 0u64;
            for i in 0..n {
                assert_eq!(out_add[i], barrett.add(a[i], b[i]));
                assert_eq!(out_sub[i], barrett.sub(a[i], b[i]));
                assert_eq!(out_scale[i], barrett.mul(a[i], alpha));
                assert_eq!(out_axpy[i], barrett.add(b[i], barrett.mul(alpha, a[i])));
                expected_dot = barrett.add(expected_dot, barrett.mul(a[i], b[i]));
            }
            assert_eq!(dot_mod(&a, &b, p), expected_dot);
        }
    }

    #[test]
    fn cached_ntt_plan_and_poly_mul_ntt_match_schoolbook() {
        let p = 998_244_353u64; // 119 * 2^23 + 1, primitive root 3
        assert_eq!(primitive_root_prime(p), Some(3));
        assert_eq!(primitive_root_prime(65_537), Some(3));

        let plan1 = ntt_plan(p, 64).unwrap();
        let plan2 = ntt_plan(p, 64).unwrap();
        assert!(Arc::ptr_eq(&plan1, &plan2));

        let a = vec![3u64, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        let b = vec![2u64, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5, 9];
        let prod = poly_mul_ntt_mod(&a, &b, p).unwrap();

        let barrett = BarrettMod::new(p).unwrap();
        let mut expected = vec![0u64; a.len() + b.len() - 1];
        for (i, &ca) in a.iter().enumerate() {
            for (j, &cb) in b.iter().enumerate() {
                expected[i + j] = barrett.add(expected[i + j], barrett.mul(ca, cb));
            }
        }
        assert_eq!(prod, expected);

        // Larger degree (len >= 64 -> 128-point NTT) and non-NTT prime rejection.
        let large_a: Vec<u64> = (0..45).map(|i| (i * 1_000_003 + 7) % p).collect();
        let large_b: Vec<u64> = (0..50).map(|i| (i * 999_983 + 11) % p).collect();
        let large_prod = poly_mul_ntt_mod(&large_a, &large_b, p).unwrap();
        let mut large_expected = vec![0u64; large_a.len() + large_b.len() - 1];
        for (i, &ca) in large_a.iter().enumerate() {
            for (j, &cb) in large_b.iter().enumerate() {
                large_expected[i + j] = barrett.add(large_expected[i + j], barrett.mul(ca, cb));
            }
        }
        assert_eq!(large_prod, large_expected);
        assert!(poly_mul_ntt_mod(&large_a, &large_b, 1_000_000_007).is_none());
        assert_eq!(poly_mul_ntt_mod(&[], &large_b, p), Some(Vec::new()));

        // Verify BarrettMod::reduce_u128 directly on 126-bit values.
        let p61 = (1u64 << 61) - 1;
        let b61 = BarrettMod::new(p61).unwrap();
        let x126 = u128::from(p61 - 1) * u128::from(p61 - 2) * 3 + 123_456_789;
        assert_eq!(b61.reduce_u128(x126), (x126 % u128::from(p61)) as u64);
    }
}
