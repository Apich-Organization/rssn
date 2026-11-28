//! Division-free modular arithmetic for moduli that fit a machine word.
//!
//! Two contexts cover every modulus representable in a `u64`:
//!
//! - [`Montgomery`] for an odd modulus `p > 1`, where a residue `x` may be
//!   carried as `x * R mod p` with `R = 2^64`. One REDC step replaces the
//!   hardware division of `(a * b) mod p`, which is the single most expensive
//!   instruction in a modular elimination.
//! - [`Barrett`] for any modulus `m >= 2`, including even moduli and composite
//!   moduli where no Montgomery inverse exists. It stores a 128-bit reciprocal
//!   `mu = floor(2^128 / m)` and needs no representation change, so a caller
//!   that already holds standard residues can multiply in place.
//!
//! This module is the *only* implementation of either algorithm in the crate.
//! The historical copies are now aliases:
//!
//! | former name | now |
//! |---|---|
//! | `finite_field::MontgomeryPrimeField` | [`Montgomery`] |
//! | `number_theory::MontU64` | [`Montgomery`] |
//! | `sparse_exact::MontGfp` | [`Montgomery`] |
//! | `qlinalg`'s private context | [`Montgomery`] |
//! | `finite_field::BarrettMod` | [`Barrett`] |
//! | `sparse_exact::BarrettGfp` | [`Barrett`] |
//!
//! Keeping one copy is not tidiness for its own sake: the copies had already
//! drifted (one reduced its inputs before adding, one did not; one handled the
//! carry out of the REDC sum, one relied on a modulus bound), so a contract bug
//! could be fixed in one place and survive in three others.
//!
//! Every function is `const`, so a context built from a literal modulus is
//! folded at compile time.
//!
//! [`Montgomery`]: crate::kernels::modular::Montgomery
//! [`Barrett`]: crate::kernels::modular::Barrett

#![allow(
    clippy::arithmetic_side_effects,
    reason = "this module is the reduction arithmetic itself; every operation is bounded by the documented range of its inputs"
)]

/// The Montgomery radix: residues are represented as `x * R mod p` with
/// `R = 2^64`.
const RADIX_BITS: u32 = 64;

/// `(a + b) mod modulus` for any inputs, reducing each operand first.
///
/// The difference form is used because `a + b` can overflow a `u64` for a
/// modulus near `2^64`: `modulus - b` is in `[1, modulus]`, so the comparison
/// and the subtraction both stay in range.
#[inline]
#[must_use]
pub(crate) const fn add_mod_raw(
    a: u64,
    b: u64,
    modulus: u64,
) -> u64 {
    let a_red = if a >= modulus {
        a % modulus
    } else {
        a
    };
    let b_red = if b >= modulus {
        b % modulus
    } else {
        b
    };
    let diff = modulus - b_red;
    if a_red >= diff {
        a_red - diff
    } else {
        a_red + b_red
    }
}

/// `(a - b) mod modulus` for any inputs, reducing each operand first.
#[inline]
#[must_use]
pub(crate) const fn sub_mod_raw(
    a: u64,
    b: u64,
    modulus: u64,
) -> u64 {
    let a_red = if a >= modulus {
        a % modulus
    } else {
        a
    };
    let b_red = if b >= modulus {
        b % modulus
    } else {
        b
    };
    if a_red >= b_red {
        a_red - b_red
    } else {
        modulus - (b_red - a_red)
    }
}

/// Montgomery arithmetic for an odd modulus `p > 1`.
///
/// # Representation
///
/// A value `x` in `[0, p)` has the Montgomery form `x * R mod p`. Addition and
/// subtraction are representation-independent; a product needs one REDC step,
/// which replaces a 64-bit hardware division (30–40 cycles, more when the
/// product needs `u128`) with two multiplies and a shift.
///
/// Methods ending in `_mont` take and return Montgomery-form residues and
/// require both inputs in `[0, p)`. The unsuffixed methods take and return
/// standard residues and accept any `u64`, reducing as needed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Montgomery {
    modulus: u64,
    n_prime: u64,
    r_mod: u64,
    r2_mod: u64,
}

impl Montgomery {
    /// Builds a context for an odd `modulus > 1`, or `None` for any other
    /// modulus (including every even modulus, where `2^64` is not invertible).
    #[must_use]
    pub const fn new(modulus: u64) -> Option<Self> {
        if modulus <= 1 || modulus.is_multiple_of(2) {
            return None;
        }
        Some(Self::build(modulus))
    }

    /// Builds a context for a modulus the caller has already established as an
    /// odd integer greater than one.
    ///
    /// The only callers are primality tests and prime-field constructors that
    /// have certified the modulus, so the branch is a formality; it is checked
    /// rather than assumed because a silent wrong answer in modular arithmetic
    /// is far more expensive than one compare.
    ///
    /// # Panics
    ///
    /// Panics if `modulus <= 1` or `modulus` is even.
    #[must_use]
    pub const fn new_prime(modulus: u64) -> Self {
        assert!(
            modulus > 1 && modulus % 2 == 1,
            "Montgomery modulus must be odd and greater than 1"
        );
        Self::build(modulus)
    }

    /// The shared constructor body; `modulus` must be odd and greater than one.
    const fn build(modulus: u64) -> Self {
        // Newton–Raphson for `modulus * inv == 1 (mod 2^64)`. Starting from
        // `inv = 1`, which is correct modulo 2, each step doubles the number of
        // correct low bits: 1, 2, 4, 8, 16, 32, 64.
        let mut inv = 1_u64;
        let mut round = 0;
        while round < 6 {
            inv = inv.wrapping_mul(2_u64.wrapping_sub(modulus.wrapping_mul(inv)));
            round += 1;
        }
        let n_prime = inv.wrapping_neg();
        let m128 = modulus as u128;
        // `2^64 mod p` written so that it is valid for every `p > 1`:
        // `(2^64 - 1) mod p` is in `[0, p)` and one more cannot reach `2p`.
        let r_mod_128 = (u64::MAX as u128) % m128 + 1;
        let r_mod = if r_mod_128 == m128 {
            0
        } else {
            r_mod_128
        } as u64;
        let r2_mod = ((r_mod as u128) * (r_mod as u128) % m128) as u64;
        Self {
            modulus,
            n_prime,
            r_mod,
            r2_mod,
        }
    }

    /// The modulus `p`.
    #[inline]
    #[must_use]
    pub const fn modulus(self) -> u64 {
        self.modulus
    }

    /// `n' = -p^(-1) mod 2^64`, the REDC multiplier.
    #[inline]
    #[must_use]
    pub const fn n_prime(self) -> u64 {
        self.n_prime
    }

    /// `R mod p = 2^64 mod p`, the Montgomery form of `1`.
    #[inline]
    #[must_use]
    pub const fn r_mod(self) -> u64 {
        self.r_mod
    }

    /// `R^2 mod p = 2^128 mod p`, used to convert into Montgomery form.
    #[inline]
    #[must_use]
    pub const fn r2_mod(self) -> u64 {
        self.r2_mod
    }

    /// `R mod p`, the Montgomery form of `1`.
    #[inline]
    #[must_use]
    pub const fn one(self) -> u64 {
        self.r_mod
    }

    /// `R mod p`, the Montgomery form of `1` (alias of [`Self::one`] kept for
    /// the call sites that name the representation explicitly).
    #[inline]
    #[must_use]
    pub const fn one_mont(self) -> u64 {
        self.r_mod
    }

    /// Montgomery reduction `t * R^(-1) mod p` for `0 <= t < p * 2^64`.
    ///
    /// The sum `t + m * p` can carry out of 128 bits for a modulus near `2^64`;
    /// `overflowing_add` reports that carry, and the wrap-around subtraction
    /// `u.wrapping_sub(p)` then yields `u + 2^64 - p`, which is the reduced
    /// value.
    #[inline]
    #[must_use]
    pub const fn redc(
        self,
        t: u128,
    ) -> u64 {
        let m = (t as u64).wrapping_mul(self.n_prime);
        let mn = (m as u128) * (self.modulus as u128);
        let (sum, carry) = t.overflowing_add(mn);
        let u = (sum >> RADIX_BITS) as u64;
        if carry || u >= self.modulus {
            u.wrapping_sub(self.modulus)
        } else {
            u
        }
    }

    /// Converts a standard residue `a` to `a * R mod p`.
    #[inline]
    #[must_use]
    pub const fn to_mont(
        self,
        a: u64,
    ) -> u64 {
        let a_red = if a >= self.modulus {
            a % self.modulus
        } else {
            a
        };
        self.redc((a_red as u128) * (self.r2_mod as u128))
    }

    /// Converts a `u128` value to its Montgomery form, reducing first.
    #[inline]
    #[must_use]
    pub const fn to_mont_u128(
        self,
        a: u128,
    ) -> u64 {
        let m128 = self.modulus as u128;
        let a_red = if a >= m128 {
            (a % m128) as u64
        } else {
            a as u64
        };
        self.redc((a_red as u128) * (self.r2_mod as u128))
    }

    /// Converts a Montgomery residue back to `[0, p)`.
    #[inline]
    #[must_use]
    pub const fn from_mont(
        self,
        a_mont: u64,
    ) -> u64 {
        self.redc(a_mont as u128)
    }

    /// `(a + b) mod p` for Montgomery residues in `[0, p)`.
    #[inline]
    #[must_use]
    pub const fn add_mont(
        self,
        a_mont: u64,
        b_mont: u64,
    ) -> u64 {
        // Overflow-free form: `p - b` cannot underflow, and the difference is
        // taken only when it does not wrap.
        let diff = self.modulus - b_mont;
        if a_mont >= diff {
            a_mont - diff
        } else {
            a_mont + b_mont
        }
    }

    /// `(a - b) mod p` for Montgomery residues in `[0, p)`.
    #[inline]
    #[must_use]
    pub const fn sub_mont(
        self,
        a_mont: u64,
        b_mont: u64,
    ) -> u64 {
        if a_mont >= b_mont {
            a_mont - b_mont
        } else {
            self.modulus - (b_mont - a_mont)
        }
    }

    /// `(-a) mod p` for a Montgomery residue in `[0, p)`.
    #[inline]
    #[must_use]
    pub const fn neg_mont(
        self,
        a_mont: u64,
    ) -> u64 {
        if a_mont == 0 {
            0
        } else {
            self.modulus - a_mont
        }
    }

    /// `(-a) mod p`.
    ///
    /// Negation commutes with the Montgomery representation, so this serves
    /// both domains; `neg_mont` is the same computation without the defensive
    /// reduction.
    #[inline]
    #[must_use]
    pub const fn neg(
        self,
        a: u64,
    ) -> u64 {
        let a_red = if a >= self.modulus {
            a % self.modulus
        } else {
            a
        };
        if a_red == 0 {
            0
        } else {
            self.modulus - a_red
        }
    }

    /// `a * b * R^(-1) mod p` for Montgomery residues in `[0, p)`.
    #[inline]
    #[must_use]
    pub const fn mul_mont(
        self,
        a_mont: u64,
        b_mont: u64,
    ) -> u64 {
        self.redc((a_mont as u128) * (b_mont as u128))
    }

    /// `acc + a * b * R^(-1) mod p`, the fused form used by inner products.
    #[inline]
    #[must_use]
    pub const fn fma_mont(
        self,
        a_mont: u64,
        b_mont: u64,
        acc_mont: u64,
    ) -> u64 {
        self.add_mont(acc_mont, self.mul_mont(a_mont, b_mont))
    }

    /// Halves a residue modulo the odd modulus.
    ///
    /// Halving commutes with the Montgomery representation (`(x * R) / 2` is
    /// the representation of `x / 2`), so this is valid in either domain.
    #[inline]
    #[must_use]
    pub const fn half(
        self,
        x: u64,
    ) -> u64 {
        if x.is_multiple_of(2) {
            x >> 1
        } else {
            (x >> 1) + (self.modulus >> 1) + 1
        }
    }

    /// `base^exp * R mod p` from a Montgomery base.
    #[inline]
    #[must_use]
    pub const fn pow_mont(
        self,
        mut base_mont: u64,
        mut exp: u64,
    ) -> u64 {
        let mut acc = self.r_mod;
        while exp > 0 {
            if exp % 2 == 1 {
                acc = self.mul_mont(acc, base_mont);
            }
            base_mont = self.mul_mont(base_mont, base_mont);
            exp >>= 1;
        }
        acc
    }

    /// `(a + b) mod p` in the standard domain; inputs may exceed `p`.
    #[inline]
    #[must_use]
    pub const fn add(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        add_mod_raw(a, b, self.modulus)
    }

    /// `(a - b) mod p` in the standard domain; inputs may exceed `p`.
    #[inline]
    #[must_use]
    pub const fn sub(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        sub_mod_raw(a, b, self.modulus)
    }

    /// `(a * b) mod p` in the standard domain.
    ///
    /// One conversion plus one REDC: `a * (b * R) * R^(-1) = a * b`, so only
    /// `b` needs to enter the Montgomery domain.
    #[inline]
    #[must_use]
    pub const fn mul(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        let a_red = if a >= self.modulus {
            a % self.modulus
        } else {
            a
        };
        let b_mont = self.to_mont(b);
        self.redc((a_red as u128) * (b_mont as u128))
    }

    /// `a^exp mod p` in the standard domain.
    #[inline]
    #[must_use]
    pub const fn pow(
        self,
        a: u64,
        exp: u64,
    ) -> u64 {
        self.from_mont(self.pow_mont(self.to_mont(a), exp))
    }

    /// `a^(-1) mod p`, or `None` when `a` is a multiple of `p`.
    #[inline]
    #[must_use]
    pub const fn inv(
        self,
        a: u64,
    ) -> Option<u64> {
        mod_inv_u64(a, self.modulus)
    }

    /// `a * b^(-1) mod p`, or `None` when `b` is not invertible.
    #[inline]
    #[must_use]
    pub const fn div(
        self,
        a: u64,
        b: u64,
    ) -> Option<u64> {
        match self.inv(b) {
            | Some(b_inv) => Some(self.mul(a, b_inv)),
            | None => None,
        }
    }
}

/// Barrett arithmetic for any modulus `m >= 2`.
///
/// The context stores `mu = floor(2^128 / m)` as a `u128` split into two words,
/// which is enough to reduce any product of two residues without a hardware
/// division. Standard residues are used throughout: unlike Montgomery form
/// there is no conversion, so a caller that already has residues (a matrix of
/// coefficients, say) pays nothing to opt in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Barrett {
    modulus: u64,
    mu_hi: u64,
    mu_lo: u64,
}

impl Barrett {
    /// Builds a reducer for `modulus >= 2`, or `None` for `0` and `1`.
    #[must_use]
    pub const fn new(modulus: u64) -> Option<Self> {
        if modulus < 2 {
            return None;
        }
        Some(Self::build(modulus))
    }

    /// Builds a reducer for a modulus the caller has already established as at
    /// least two.
    ///
    /// # Panics
    ///
    /// Panics if `modulus < 2`.
    #[must_use]
    pub const fn new_prime(modulus: u64) -> Self {
        assert!(modulus >= 2, "Barrett modulus must be at least 2");
        Self::build(modulus)
    }

    /// The shared constructor body; `modulus` must be at least two.
    const fn build(modulus: u64) -> Self {
        let m = modulus as u128;
        // floor(2^128 / m) = floor((2^128 - 1) / m) plus one when the remainder
        // of that division is `m - 1`, i.e. when 2^128 is an exact multiple.
        let quotient = u128::MAX / m;
        let remainder = u128::MAX % m;
        let mu = quotient + if remainder + 1 == m { 1 } else { 0 };
        Self {
            modulus,
            mu_hi: (mu >> RADIX_BITS) as u64,
            mu_lo: mu as u64,
        }
    }

    /// The modulus `m`.
    #[inline]
    #[must_use]
    pub const fn modulus(self) -> u64 {
        self.modulus
    }

    /// The high word of `mu = floor(2^128 / m)`, which is `floor(2^64 / m)`.
    #[inline]
    #[must_use]
    pub const fn mu(self) -> u64 {
        self.mu_hi
    }

    /// Reduces a `u64` modulo `m` without hardware division.
    #[inline]
    #[must_use]
    pub const fn reduce_u64(
        self,
        x: u64,
    ) -> u64 {
        if x < self.modulus {
            return x;
        }
        let lo_hi = (((x as u128) * (self.mu_lo as u128)) >> RADIX_BITS) as u64;
        let q = ((((x as u128) * (self.mu_hi as u128)) + (lo_hi as u128)) >> RADIX_BITS) as u64;
        let mut r = x.wrapping_sub(q.wrapping_mul(self.modulus));
        if r >= self.modulus {
            r -= self.modulus;
        }
        r
    }

    /// Reduces a `u128` modulo `m`.
    ///
    /// For `x < 2^126` — every product of two residues below `2^63` — the
    /// reduction is pure multiply-and-shift with at most two conditional
    /// subtractions. Above that the high word no longer fits the bound the
    /// estimate relies on, so the hardware division is used; that path only
    /// runs for moduli above `2^63`.
    #[inline]
    #[must_use]
    pub const fn reduce_u128(
        self,
        x: u128,
    ) -> u64 {
        let m = self.modulus as u128;
        if x < m {
            return x as u64;
        }
        let x_hi = (x >> RADIX_BITS) as u64;
        if x_hi >= self.modulus {
            return (x % m) as u64;
        }
        let x_lo = x as u64;
        let p_hh = (x_hi as u128) * (self.mu_hi as u128);
        let p_hl = (x_hi as u128) * (self.mu_lo as u128);
        let p_lh = (x_lo as u128) * (self.mu_hi as u128);
        let p_ll_hi = ((x_lo as u128) * (self.mu_lo as u128)) >> RADIX_BITS;
        let mid = (p_hl & (u64::MAX as u128)) + (p_lh & (u64::MAX as u128)) + p_ll_hi;
        let q = p_hh + (p_hl >> RADIX_BITS) + (p_lh >> RADIX_BITS) + (mid >> RADIX_BITS);
        let mut r = x.wrapping_sub(q.wrapping_mul(m));
        if r >= m {
            r -= m;
            if r >= m {
                r %= m;
            }
        }
        r as u64
    }

    /// `(a + b) mod m`.
    #[inline]
    #[must_use]
    pub const fn add(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        let a_red = self.reduce_u64(a);
        let b_red = self.reduce_u64(b);
        let diff = self.modulus - b_red;
        if a_red >= diff {
            a_red - diff
        } else {
            a_red + b_red
        }
    }

    /// `(a - b) mod m`.
    #[inline]
    #[must_use]
    pub const fn sub(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        let a_red = self.reduce_u64(a);
        let b_red = self.reduce_u64(b);
        if a_red >= b_red {
            a_red - b_red
        } else {
            self.modulus - (b_red - a_red)
        }
    }

    /// `(-a) mod m`.
    #[inline]
    #[must_use]
    pub const fn neg(
        self,
        a: u64,
    ) -> u64 {
        let a_red = self.reduce_u64(a);
        if a_red == 0 {
            0
        } else {
            self.modulus - a_red
        }
    }

    /// `(a * b) mod m`.
    #[inline]
    #[must_use]
    pub const fn mul(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        let a_red = self.reduce_u64(a);
        let b_red = self.reduce_u64(b);
        self.reduce_u128((a_red as u128) * (b_red as u128))
    }

    /// `a^exp mod m`.
    #[inline]
    #[must_use]
    pub const fn pow(
        self,
        a: u64,
        mut exp: u64,
    ) -> u64 {
        let mut base = self.reduce_u64(a);
        let mut acc = 1 % self.modulus;
        while exp > 0 {
            if exp % 2 == 1 {
                acc = self.mul(acc, base);
            }
            base = self.mul(base, base);
            exp >>= 1;
        }
        acc
    }

    /// `a^(-1) mod m`, or `None` when `gcd(a, m) != 1`.
    #[inline]
    #[must_use]
    pub const fn inv(
        self,
        a: u64,
    ) -> Option<u64> {
        mod_inv_u64(self.reduce_u64(a), self.modulus)
    }

    /// `a * b^(-1) mod m`, or `None` when `b` is not invertible.
    #[inline]
    #[must_use]
    pub const fn div(
        self,
        a: u64,
        b: u64,
    ) -> Option<u64> {
        match self.inv(b) {
            | Some(b_inv) => Some(self.mul(a, b_inv)),
            | None => None,
        }
    }
}

/// The extended Euclidean algorithm on `u64` inputs.
///
/// Returns `(g, x, y)` with `g = gcd(a, b)` and `a * x + b * y = g`. The
/// Bézout coefficients of `u64` inputs stay within `i128`.
#[must_use]
pub const fn extended_gcd_u64(
    a: u64,
    b: u64,
) -> (u64, i128, i128) {
    euclid_with_coefficients::<true>(a, b)
}

// The inverse uses only the coefficient of the reduced input. Discard
// the other recurrence at compile time, without duplicating Euclid.
const fn euclid_with_coefficients<const TRACK_FIRST: bool>(
    a: u64,
    b: u64,
) -> (u64, i128, i128) {
    let mut r0 = a;
    let mut r1 = b;
    let mut s0: i128 = 1;
    let mut s1: i128 = 0;
    let mut t0: i128 = 0;
    let mut t1: i128 = 1;
    while r1 != 0 {
        let q = r0 / r1;
        let r2 = r0 % r1;
        let q_i128 = q as i128;
        if TRACK_FIRST {
            let s2 = s0 - q_i128 * s1;
            s0 = s1;
            s1 = s2;
        }
        let t2 = t0 - q_i128 * t1;
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    (r0, s0, t0)
}

/// `a^(-1) mod modulus`, or `None` when `gcd(a, modulus) != 1` or
/// `modulus < 2`.
///
/// The extended Euclidean algorithm terminates in `O(log modulus)` steps with
/// small quotients, which is faster than Fermat exponentiation and works for
/// composite moduli.
#[inline]
#[must_use]
pub const fn mod_inv_u64(
    a: u64,
    modulus: u64,
) -> Option<u64> {
    if modulus <= 1 {
        return None;
    }
    let a_red = if a >= modulus {
        a % modulus
    } else {
        a
    };
    if a_red == 0 {
        return None;
    }
    let (g, _, x) = euclid_with_coefficients::<false>(modulus, a_red);
    if g != 1 {
        return None;
    }
    // For 0 < a_red < modulus, Euclid gives |x| < modulus.
    let inverse = if x < 0 {
        x + modulus as i128
    } else {
        x
    };
    Some(inverse as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Moduli that exercise the edges: tiny, a Mersenne prime, a prime just
    /// below `2^63`, the largest 64-bit prime, and an even composite.
    const MODULI: [u64; 8] = [
        3,
        5,
        97,
        65_537,
        1_000_000_007,
        (1 << 61) - 1,
        0x7FFF_FFFF_FFFF_FFC5,
        0xFFFF_FFFF_FFFF_FFC5,
    ];

    fn sample_values(m: u64) -> [u64; 6] {
        [0, 1, 2, 7 % m, m / 2, m - 1]
    }

    #[test]
    fn montgomery_rejects_even_and_trivial_moduli() {
        assert!(Montgomery::new(0).is_none());
        assert!(Montgomery::new(1).is_none());
        assert!(Montgomery::new(2).is_none());
        assert!(Montgomery::new(1_000_000_006).is_none());
        assert!(Montgomery::new(3).is_some());
        assert!(Montgomery::new(u64::MAX).is_some());
    }

    #[test]
    fn montgomery_constants_are_consistent() {
        for p in MODULI {
            let Some(m) = Montgomery::new(p) else {
                continue;
            };
            // n' is the negated inverse of p modulo 2^64.
            assert_eq!(p.wrapping_mul(m.n_prime()), u64::MAX, "p = {p}");
            // R and R^2 are the residues the conversions rely on.
            let r_want = ((1_u128 << 64) % u128::from(p)) as u64;
            assert_eq!(m.r_mod(), r_want, "p = {p}");
            let r = u128::from(m.r_mod());
            assert_eq!(m.r2_mod(), (r * r % u128::from(p)) as u64, "p = {p}");
        }
    }

    #[test]
    fn montgomery_round_trip_and_products_match_the_reference() {
        for p in MODULI {
            let Some(m) = Montgomery::new(p) else {
                continue;
            };
            for a in sample_values(p) {
                let a_m = m.to_mont(a);
                assert!(a_m < p);
                assert_eq!(m.from_mont(a_m), a, "p = {p}, a = {a}");
                for b in sample_values(p) {
                    let b_m = m.to_mont(b);
                    let want = (u128::from(a) * u128::from(b) % u128::from(p)) as u64;
                    assert_eq!(m.mul_mont(a_m, b_m), m.to_mont(want), "p = {p}, {a} * {b}");
                    assert_eq!(m.mul(a, b), want, "standard domain p = {p}, {a} * {b}");
                    // The `_mont` results are Montgomery-form, so the oracle is
                    // the Montgomery form of the standard-domain answer.
                    let sum = (u128::from(a) + u128::from(b)) % u128::from(p);
                    assert_eq!(m.add_mont(a_m, b_m), m.to_mont(sum as u64));
                    assert_eq!(m.add(a, b), sum as u64);
                    let diff = (u128::from(a) + u128::from(p) - u128::from(b)) % u128::from(p);
                    assert_eq!(m.sub_mont(a_m, b_m), m.to_mont(diff as u64));
                    assert_eq!(m.sub(a, b), diff as u64);
                    let neg = (u128::from(p) - u128::from(a)) % u128::from(p);
                    assert_eq!(m.neg_mont(a_m), m.to_mont(neg as u64));
                    let fma = (u128::from(want) + u128::from(b)) % u128::from(p);
                    assert_eq!(m.fma_mont(a_m, b_m, b_m), m.to_mont(fma as u64));
                }
                // Halving is exact for odd p.
                let doubled = m.add(a, a);
                assert_eq!(m.half(doubled), a, "p = {p}, a = {a}");
                // The `u128` entry point reduces before converting.
                assert_eq!(m.from_mont(m.to_mont_u128(u128::from(a))), a);
                let wide = u128::from(p) * 3 + u128::from(a);
                assert_eq!(m.from_mont(m.to_mont_u128(wide)), a);
            }
        }
    }

    #[test]
    fn montgomery_pow_and_inverse_match_the_reference() {
        for p in MODULI {
            let Some(m) = Montgomery::new(p) else {
                continue;
            };
            for a in sample_values(p) {
                for e in [0_u64, 1, 2, 3, 17, 1000] {
                    let mut want = 1_u128 % u128::from(p);
                    let base = u128::from(a);
                    for _ in 0..e {
                        want = want * base % u128::from(p);
                    }
                    assert_eq!(m.pow(a, e), want as u64, "p = {p}, {a}^{e}");
                }
                if a != 0 {
                    let inv = m.inv(a).expect("a is coprime to a prime modulus");
                    assert_eq!(m.mul(a, inv), 1, "p = {p}, a = {a}");
                    assert_eq!(m.div(1, a), Some(inv));
                }
            }
        }
        // A composite modulus has non-invertible elements.
        let m = Montgomery::new(15).expect("15 is odd");
        assert_eq!(m.inv(5), None);
        assert_eq!(m.div(1, 3), None);
        assert_eq!(m.div(1, 2), Some(8));
    }

    #[test]
    fn barrett_agrees_with_the_reference() {
        for p in [2_u64, 3, 4, 97, 1_000_000_007, (1 << 61) - 1, u64::MAX] {
            let b = Barrett::new(p).expect("modulus is at least 2");
            assert_eq!(b.modulus(), p);
            for a in sample_values(p) {
                assert_eq!(b.reduce_u64(a), a % p, "p = {p}, a = {a}");
                for c in sample_values(p) {
                    let want = (u128::from(a) * u128::from(c) % u128::from(p)) as u64;
                    assert_eq!(b.mul(a, c), want, "p = {p}, {a} * {c}");
                    assert_eq!(
                        b.add(a, c),
                        ((u128::from(a) + u128::from(c)) % u128::from(p)) as u64
                    );
                    assert_eq!(
                        b.sub(a, c),
                        ((u128::from(a) + u128::from(p) - u128::from(c)) % u128::from(p)) as u64
                    );
                }
                assert_eq!(b.neg(a), (u128::from(p) - u128::from(a)) as u64 % p);
                for e in [0_u64, 1, 5, 4096] {
                    let mut want = 1_u128 % u128::from(p);
                    for _ in 0..e {
                        want = want * u128::from(a) % u128::from(p);
                    }
                    assert_eq!(b.pow(a, e), want as u64, "p = {p}, {a}^{e}");
                }
            }
        }
        assert!(Barrett::new(0).is_none());
        assert!(Barrett::new(1).is_none());
    }

    #[test]
    fn barrett_reduces_products_of_large_residues() {
        // 2^63 < p < 2^64 forces the wide path, where a naive `mu` estimate
        // would be off by more than one.
        for p in [0x8000_0000_0000_0001_u64, 0xFFFF_FFFF_FFFF_FFC5, u64::MAX] {
            let b = Barrett::new(p).expect("modulus is at least 2");
            for a in [p - 1, p - 2, p / 3, p / 2 + 1] {
                for c in [p - 1, p / 7, 3] {
                    let want = (u128::from(a) * u128::from(c) % u128::from(p)) as u64;
                    assert_eq!(b.mul(a, c), want, "p = {p}, {a} * {c}");
                }
            }
        }
    }

    #[test]
    fn barrett_and_montgomery_agree_where_both_apply() {
        for p in MODULI {
            let Some(m) = Montgomery::new(p) else {
                continue;
            };
            let b = Barrett::new(p).expect("modulus is at least 2");
            for a in sample_values(p) {
                for c in sample_values(p) {
                    assert_eq!(m.mul(a, c), b.mul(a, c), "p = {p}, {a} * {c}");
                }
            }
        }
    }

    #[test]
    fn extended_gcd_reports_bezout_coefficients() {
        for (a, b) in [(3_u64, 5_u64), (240, 46), (0, 7), (7, 0), (17, 17)] {
            let (g, x, y) = extended_gcd_u64(a, b);
            assert_eq!(g, crate::kernels::number_theory::gcd(a, b));
            assert_eq!(i128::from(a) * x + i128::from(b) * y, i128::from(g));
        }
    }
}
