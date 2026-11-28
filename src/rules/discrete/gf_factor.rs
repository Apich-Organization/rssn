//! Factorisation of polynomials over GF(p).
//!
//! Polynomials use the representation of [`finite_field`](super::finite_field):
//! coefficient lists with the **highest degree first**, no leading zeros,
//! the zero polynomial being `list()`. The modulus `p` must be prime.
//!
//! | operator | value |
//! |---|---|
//! | `gfp_monic(f, p)`, `gfp_derivative(f, p)` | `f` divided by its leading coefficient; `f'` |
//! | `gfp_powmod(f, e, m, p)` | `f^e mod m` by repeated squaring (`e >= 0`) |
//! | `gfp_invmod(f, m, p)` | the inverse of `f` modulo `m` (`gcd(f, m) = 1`) |
//! | `gfp_isirreducible(f, p)` | Rabin's test: `x^(p^n) = x` and `gcd(x^(p^(n/q)) - x, f) = 1` for primes `q` dividing `n` (in `finite_field`) |
//! | `gfp_squarefree(f, p)` | `list(list(g1, m1), ...)`: monic, squarefree, pairwise coprime `g_i` with `f = lc(f) prod g_i^m_i` (Yun's algorithm, with the `p`-th root step in characteristic `p`) |
//! | `gfp_ddf(f, p)` | distinct-degree factorisation of a squarefree `f`: `list(list(g_d, d), ...)`, `g_d` the product of the irreducible factors of degree `d` |
//! | `gfp_edf(f, d, p)` | equal-degree factorisation (Cantor–Zassenhaus): the irreducible factors of a squarefree `f` all of degree `d` |
//! | `gfp_factor(f, p)` | `list(lc, list(g1, m1), ...)`: the complete factorisation, factors ordered by degree |
//! | `gfp_berlekamp(f, p)` | the same factorisation by Berlekamp's algorithm (small `p` only) |
//! | `factor_mod(f, x, p)` | like `gfp_factor`, but `f` may be an expression in `x` with integer coefficients, in which case the answer is the product `lc * g1^m1 * ...` of expressions with coefficients in `[0, p)` |
//! | `is_irreducible_mod(f, x, p)` | `gfp_isirreducible` for an expression (or list) |
//!
//! The equal-degree step draws its random polynomials from a fixed
//! deterministic generator, so results are reproducible.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::V;
use super::big;
use super::def_doc;
use super::finite_field::Poly;
use super::finite_field::PolyU64;
use super::finite_field::fast_barrett;
use super::finite_field::monic as monic_poly;
use super::finite_field::monic_u64;
use super::finite_field::norm;
use super::finite_field::norm_u64_in_place;
use super::finite_field::padd_u64;
use super::finite_field::pdivmod;
use super::finite_field::pdivmod_u64;
use super::finite_field::pegcd;
use super::finite_field::pgcd_u64;
use super::finite_field::pmul;
use super::finite_field::pmul_u64;
use super::finite_field::poly_from_u64;
use super::finite_field::poly_to_u64;
use super::finite_field::poly_value;
use super::finite_field::psub;
use super::finite_field::psub_u64;
use super::finite_field::read_poly;
use super::finite_field::reduce;
use super::finite_field::reduce_u64;
use super::finite_field::xinv;
use super::items;
use super::mod_inverse;
use super::modulo;
use super::pow;
use super::prod;
use super::sum;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::finite_field::BarrettMod;
use crate::kernels::finite_field::MontgomeryPrimeField;
use crate::kernels::finite_field::slice_axpy_mod;
use crate::kernels::par;
use crate::kernels::reuse::global_pool;
use crate::rules::number_theory::is_prime;
use crate::rules::poly::best;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::from_term;

/// A prime modulus.
fn prime(
    cx: &Cx<'_>,
    n: NodeId,
) -> Option<BigInt> {
    big(cx.graph, n).filter(is_prime)
}

fn one() -> Poly {
    vec![BigInt::one()]
}

fn x_poly() -> Poly {
    vec![BigInt::one(), BigInt::zero()]
}

/// A non-zero constant (after normalisation).
const fn is_constant(f: &Poly) -> bool {
    f.len() == 1
}

/// A non-zero constant in `u64` representation (after normalisation).
const fn is_constant_u64(f: &[u64]) -> bool {
    f.len() == 1
}

fn monic(
    f: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    monic_poly(f, p)
}

fn derivative_bigint(
    f: &Poly,
    p: &BigInt,
) -> Poly {
    let n = f.len();
    norm(
        (0..n.saturating_sub(1))
            .map(|j| &f[j] * BigInt::from(n - 1 - j))
            .collect(),
        p,
    )
}

#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "loop index j runs over 0..n - 1 where n = f.len() >= 2"
)]
fn derivative_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> PolyU64 {
    let n = f.len();
    if n <= 1 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n - 1);
    for (j, &c) in f[..n - 1].iter().enumerate() {
        let power = barrett.reduce_u64((n - 1 - j) as u64);
        out.push(barrett.mul(c, power));
    }
    norm_u64_in_place(&mut out, barrett);
    out
}

fn derivative(
    f: &Poly,
    p: &BigInt,
) -> Poly {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        return poly_from_u64(&derivative_u64(&fu, barrett));
    }
    derivative_bigint(f, p)
}

/// The monic gcd.
fn gcd(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    let (d, _, _) = pegcd(f, g, p)?;
    monic(&d, p)
}

fn quotient(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    Some(pdivmod(f, g, p)?.0)
}

fn quotient_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    Some(pdivmod_u64(f, g, barrett)?.0)
}

/// `f^e mod m` (arbitrary-precision fallback).
fn powmod_bigint(
    f: &Poly,
    e: &BigInt,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if e.is_negative() {
        return None;
    }
    let mut base = reduce(f, m, p)?;
    let mut acc = reduce(&one(), m, p)?;
    let mut e = e.clone();
    while !e.is_zero() {
        if e.bit(0) {
            acc = reduce(&pmul(&acc, &base, p), m, p)?;
        }
        base = reduce(&pmul(&base, &base, p), m, p)?;
        e >>= 1;
    }
    Some(acc)
}

/// `f^e mod m` over `u64` residues using little-endian 64-bit limbs of `e`.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "last_idx < limbs_le.len() and bit_idx / 64 <= last_idx"
)]
fn powmod_u64_limbs(
    f: &[u64],
    limbs_le: &[u64],
    m: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    let mut base = reduce_u64(f, m, barrett)?;
    let mut acc = reduce_u64(&[1u64], m, barrett)?;
    let Some(last_idx) = limbs_le.iter().rposition(|&w| w != 0) else {
        return Some(acc);
    };
    let total_bits = last_idx * 64 + (64 - limbs_le[last_idx].leading_zeros() as usize);
    for bit_idx in 0..total_bits {
        if ((limbs_le[bit_idx / 64] >> (bit_idx % 64)) & 1) == 1 {
            acc = reduce_u64(&pmul_u64(&acc, &base, barrett), m, barrett)?;
        }
        if bit_idx + 1 < total_bits {
            base = reduce_u64(&pmul_u64(&base, &base, barrett), m, barrett)?;
        }
    }
    Some(acc)
}

#[allow(clippy::single_call_fn)]
fn powmod_u64(
    f: &[u64],
    e: &BigInt,
    m: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    if e.is_negative() {
        return None;
    }
    if let Some(eu) = e.to_u64() {
        return powmod_u64_limbs(f, &[eu], m, barrett);
    }
    let (_, limbs) = e.to_u64_digits();
    powmod_u64_limbs(f, &limbs, m, barrett)
}

/// `f^e mod m`.
fn powmod(
    f: &Poly,
    e: &BigInt,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return powmod_u64(&fu, e, &mu, barrett).map(|res| poly_from_u64(&res));
    }
    powmod_bigint(f, e, m, p)
}

/// The polynomial `g` with `g(x^p) = c(x)` when only degrees divisible by
/// `p` occur in `c`.
#[allow(clippy::single_call_fn)]
fn p_th_root(
    c: &Poly,
    p: usize,
) -> Poly {
    let degree = c.len() - 1;
    (0..=degree / p)
        .rev()
        .map(|k| c[degree - k * p].clone())
        .collect()
}

/// The `u64` polynomial `g` with `g(x^p) = c(x)`.
#[allow(
    clippy::single_call_fn,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "k <= degree / p ensures degree - k * p < c.len()"
)]
fn p_th_root_u64(
    c: &[u64],
    p: usize,
) -> PolyU64 {
    let degree = c.len() - 1;
    (0..=degree / p).rev().map(|k| c[degree - k * p]).collect()
}

/// Yun's squarefree decomposition of a monic polynomial (`BigInt` reference).
fn squarefree_bigint(
    f: &Poly,
    p: &BigInt,
) -> Option<Vec<(Poly, usize)>> {
    let mut out: Vec<(Poly, usize)> = Vec::new();
    if is_constant(f) || f.is_empty() {
        return Some(out);
    }
    let df = derivative_bigint(f, p);
    let mut c = gcd(f, &df, p)?;
    let mut w = quotient(f, &c, p)?;
    let mut i = 1;
    while !is_constant(&w) {
        let y = gcd(&w, &c, p)?;
        let factor = quotient(&w, &y, p)?;
        if !is_constant(&factor) {
            out.push((factor, i));
        }
        w = y;
        c = quotient(&c, &w, p)?;
        i += 1;
    }
    if !is_constant(&c) {
        // c = g(x)^p
        let pu = p.to_usize()?;
        let root = p_th_root(&c, pu);
        for (g, m) in squarefree_bigint(&monic(&root, p)?, p)? {
            out.push((g, m * pu));
        }
    }
    // Merge equal multiplicities.
    out.sort_by_key(|(_, m)| *m);
    let mut merged: Vec<(Poly, usize)> = Vec::new();
    for (g, m) in out {
        match merged.last_mut() {
            | Some((h, k)) if *k == m => *h = pmul(h, &g, p),
            | _ => merged.push((g, m)),
        }
    }
    Some(merged)
}

/// Yun's squarefree decomposition of a monic polynomial over `u64` residues.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "multiplicity products stay bounded by polynomial degree"
)]
fn squarefree_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> Option<Vec<(PolyU64, usize)>> {
    let mut out: Vec<(PolyU64, usize)> = Vec::new();
    if is_constant_u64(f) || f.is_empty() {
        return Some(out);
    }
    let df = derivative_u64(f, barrett);
    let mut c = pgcd_u64(f, &df, barrett)?;
    let mut w = quotient_u64(f, &c, barrett)?;
    let mut i = 1usize;
    while !is_constant_u64(&w) {
        let y = pgcd_u64(&w, &c, barrett)?;
        let factor = quotient_u64(&w, &y, barrett)?;
        if !is_constant_u64(&factor) {
            out.push((factor, i));
        }
        w = y;
        c = quotient_u64(&c, &w, barrett)?;
        i += 1;
    }
    if !is_constant_u64(&c) {
        let pu = usize::try_from(barrett.modulus()).ok()?;
        let root = p_th_root_u64(&c, pu);
        for (g, m) in squarefree_u64(&monic_u64(&root, barrett)?, barrett)? {
            out.push((g, m * pu));
        }
    }
    out.sort_by_key(|(_, m)| *m);
    let mut merged: Vec<(PolyU64, usize)> = Vec::new();
    for (g, m) in out {
        match merged.last_mut() {
            | Some((h, k)) if *k == m => *h = pmul_u64(h, &g, barrett),
            | _ => merged.push((g, m)),
        }
    }
    Some(merged)
}

/// Distinct-degree factorisation of a monic squarefree polynomial (`BigInt` reference).
fn ddf_bigint(
    f: &Poly,
    p: &BigInt,
) -> Option<Vec<(Poly, usize)>> {
    let mut out = Vec::new();
    let mut h = x_poly();
    let mut rest = f.clone();
    let mut d = 1;
    while rest.len() > 2 * d {
        h = powmod_bigint(&h, p, &rest, p)?;
        let g = gcd(&psub(&h, &x_poly(), p), &rest, p)?;
        if !is_constant(&g) {
            rest = quotient(&rest, &g, p)?;
            h = reduce(&h, &rest, p)?;
            out.push((g, d));
        }
        d += 1;
    }
    if rest.len() > 1 {
        let degree = rest.len() - 1;
        out.push((rest, degree));
    }
    Some(out)
}

/// Distinct-degree factorisation of a monic squarefree polynomial over `u64` residues.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "rest.len() > 1 before computing rest.len() - 1"
)]
fn ddf_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> Option<Vec<(PolyU64, usize)>> {
    let mut out = Vec::new();
    let x = [1u64, 0u64];
    let mut h = x.to_vec();
    let mut rest = f.to_vec();
    let mut d = 1usize;
    let p_limb = [barrett.modulus()];
    while rest.len() > 2 * d {
        h = powmod_u64_limbs(&h, &p_limb, &rest, barrett)?;
        let g = pgcd_u64(&psub_u64(&h, &x, barrett), &rest, barrett)?;
        if !is_constant_u64(&g) {
            rest = quotient_u64(&rest, &g, barrett)?;
            h = reduce_u64(&h, &rest, barrett)?;
            out.push((g, d));
        }
        d += 1;
    }
    if rest.len() > 1 {
        let degree = rest.len() - 1;
        out.push((rest, degree));
    }
    Some(out)
}

/// A deterministic generator of residues.
struct Rng(u64);

impl Rng {
    const fn word(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn below(
        &mut self,
        bound: &BigInt,
    ) -> BigInt {
        let words = bound.bits() / 64 + 2;
        let mut value = BigInt::zero();
        for _ in 0..words {
            value = (value << 64) + BigInt::from(self.word());
        }
        value % bound
    }

    /// Produces the exact same residue modulo `barrett.modulus() < 2^62` as
    /// [`Self::below`], consuming two 64-bit RNG words per draw.
    #[allow(clippy::single_call_fn)]
    fn below_u64(
        &mut self,
        barrett: BarrettMod,
    ) -> u64 {
        let hi = self.word();
        let lo = self.word();
        barrett.reduce_u128((u128::from(hi) << 64) | u128::from(lo))
    }
}

/// Cantor–Zassenhaus (`BigInt` reference): the irreducible factors of the
/// monic squarefree `g`, all of degree `d`.
fn edf_bigint(
    g: &Poly,
    d: usize,
    p: &BigInt,
    rng: &mut Rng,
    out: &mut Vec<Poly>,
) -> Option<()> {
    let n = g.len().checked_sub(1)?;
    if n == d {
        out.push(g.clone());
        return Some(());
    }
    let two = BigInt::from(2);
    let exponent = (p.pow(u32::try_from(d).ok()?) - BigInt::one()) / &two;
    for _ in 0..500 {
        let a = norm((0..n).map(|_| rng.below(p)).collect(), p);
        let t = if *p == two {
            // the trace map a + a^2 + a^4 + ... + a^(2^(d-1))
            let mut s = reduce(&a, g, p)?;
            let mut t = s.clone();
            for _ in 1..d {
                s = reduce(&pmul(&s, &s, p), g, p)?;
                t = super::finite_field::padd(&t, &s, p);
            }
            t
        } else {
            psub(&powmod_bigint(&a, &exponent, g, p)?, &one(), p)
        };
        let h = gcd(&t, g, p)?;
        if !is_constant(&h) && h.len() != g.len() {
            let rest = quotient(g, &h, p)?;
            edf_bigint(&h, d, p, rng, out)?;
            return edf_bigint(&rest, d, p, rng, out);
        }
    }
    None
}

/// Precomputed Cantor–Zassenhaus exponent `(p^d - 1) / 2` in little-endian `u64` limbs.
enum EdfExponent {
    CharTwo,
    Word(u64),
    U128([u64; 2]),
    Big(Vec<u64>),
}

impl EdfExponent {
    #[allow(
        clippy::single_call_fn,
        clippy::arithmetic_side_effects,
        reason = "p >= 3 ensures pd >= 3 so pd - 1 never underflows"
    )]
    fn new(
        p: u64,
        d: usize,
    ) -> Option<Self> {
        if p == 2 {
            return Some(Self::CharTwo);
        }
        let d_u32 = u32::try_from(d).ok()?;
        if d_u32 == 1 {
            return Some(Self::Word((p - 1) / 2));
        }
        if let Some(pd) = u128::from(p).checked_pow(d_u32) {
            let exp = (pd - 1) / 2;
            return Some(Self::U128([exp as u64, (exp >> 64) as u64]));
        }
        let exp: BigInt = (BigInt::from(p).pow(d_u32) - BigInt::one()) >> 1;
        let (_, limbs) = exp.to_u64_digits();
        Some(Self::Big(limbs))
    }

    #[allow(clippy::single_call_fn)]
    fn limbs(&self) -> &[u64] {
        match self {
            | Self::CharTwo => &[],
            | Self::Word(w) => std::slice::from_ref(w),
            | Self::U128(pair) => pair,
            | Self::Big(limbs) => limbs,
        }
    }
}

fn edf_u64_rec(
    g: &[u64],
    d: usize,
    barrett: BarrettMod,
    exp: &EdfExponent,
    rng: &mut Rng,
    out: &mut Vec<PolyU64>,
) -> Option<()> {
    let n = g.len().checked_sub(1)?;
    if n == d {
        out.push(g.to_vec());
        return Some(());
    }
    for _ in 0..500 {
        let mut a: PolyU64 = (0..n).map(|_| rng.below_u64(barrett)).collect();
        norm_u64_in_place(&mut a, barrett);
        let t = if matches!(exp, EdfExponent::CharTwo) {
            let mut s = reduce_u64(&a, g, barrett)?;
            let mut trace = s.clone();
            for _ in 1..d {
                s = reduce_u64(&pmul_u64(&s, &s, barrett), g, barrett)?;
                trace = padd_u64(&trace, &s, barrett);
            }
            trace
        } else {
            psub_u64(
                &powmod_u64_limbs(&a, exp.limbs(), g, barrett)?,
                &[1u64],
                barrett,
            )
        };
        let h = pgcd_u64(&t, g, barrett)?;
        if !is_constant_u64(&h) && h.len() != g.len() {
            let rest = quotient_u64(g, &h, barrett)?;
            edf_u64_rec(&h, d, barrett, exp, rng, out)?;
            return edf_u64_rec(&rest, d, barrett, exp, rng, out);
        }
    }
    None
}

/// Cantor–Zassenhaus over `u64` residues, preserving the exact RNG draw order.
fn edf_u64(
    g: &[u64],
    d: usize,
    barrett: BarrettMod,
    rng: &mut Rng,
    out: &mut Vec<PolyU64>,
) -> Option<()> {
    let exp = EdfExponent::new(barrett.modulus(), d)?;
    edf_u64_rec(g, d, barrett, &exp, rng, out)
}

/// Rabin's irreducibility test (`p` prime, `BigInt` reference).
fn is_irreducible_bigint(
    f: &Poly,
    p: &BigInt,
) -> Option<bool> {
    let Some(degree) = f.len().checked_sub(1) else {
        return Some(false);
    };
    if degree == 0 {
        return Some(false);
    }
    let f = monic(f, p)?;
    let x = x_poly();
    // x^(p^k) mod f for the k needed, built incrementally.
    let mut primes = Vec::new();
    let mut m = degree;
    let mut q = 2;
    while q * q <= m {
        if m % q == 0 {
            primes.push(q);
            while m % q == 0 {
                m /= q;
            }
        }
        q += 1;
    }
    if m > 1 {
        primes.push(m);
    }
    let checkpoints: Vec<usize> = primes.iter().map(|q| degree / q).collect();
    let mut h = reduce(&x, &f, p)?;
    for k in 1..=degree {
        h = powmod_bigint(&h, p, &f, p)?;
        if checkpoints.contains(&k) {
            let g = gcd(&psub(&h, &x, p), &f, p)?;
            if !is_constant(&g) {
                return Some(false);
            }
        }
    }
    Some(reduce(&psub(&h, &x, p), &f, p)?.is_empty())
}

/// Rabin's irreducibility test (`p` prime) over `u64` residues.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "q >= 2 divides m and degree / q is exact"
)]
fn is_irreducible_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> Option<bool> {
    let Some(degree) = f.len().checked_sub(1) else {
        return Some(false);
    };
    if degree == 0 {
        return Some(false);
    }
    let f = monic_u64(f, barrett)?;
    let x = [1u64, 0u64];
    let mut primes = Vec::new();
    let mut m = degree;
    let mut q = 2usize;
    while q * q <= m {
        if m.is_multiple_of(q) {
            primes.push(q);
            while m.is_multiple_of(q) {
                m /= q;
            }
        }
        q += 1;
    }
    if m > 1 {
        primes.push(m);
    }
    let checkpoints: Vec<usize> = primes.iter().map(|q| degree / q).collect();
    let mut h = reduce_u64(&x, &f, barrett)?;
    let p_limb = [barrett.modulus()];
    for k in 1..=degree {
        h = powmod_u64_limbs(&h, &p_limb, &f, barrett)?;
        if checkpoints.contains(&k) {
            let g = pgcd_u64(&psub_u64(&h, &x, barrett), &f, barrett)?;
            if !is_constant_u64(&g) {
                return Some(false);
            }
        }
    }
    Some(reduce_u64(&psub_u64(&h, &x, barrett), &f, barrett)?.is_empty())
}

/// Rabin's irreducibility test (`p` prime).
pub(super) fn is_irreducible(
    f: &Poly,
    p: &BigInt,
) -> Option<bool> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        return is_irreducible_u64(&fu, barrett);
    }
    is_irreducible_bigint(f, p)
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum SplitMethod {
    CantorZassenhaus,
    Berlekamp,
}

/// Complete factorisation over `u64` residues: leading coefficient and
/// `(factor, multiplicity)` pairs ordered by degree then lexicographically.
#[allow(clippy::single_call_fn)]
fn factor_u64(
    f: &[u64],
    barrett: BarrettMod,
    method: SplitMethod,
) -> Option<(u64, Vec<(PolyU64, usize)>)> {
    let lead = *f.first()?;
    let g = monic_u64(f, barrett)?;
    let mut rng = Rng(0x5eed_1234_abcd_9876);
    let mut out: Vec<(PolyU64, usize)> = Vec::new();
    for (sq, m) in squarefree_u64(&g, barrett)? {
        for (gd, d) in ddf_u64(&sq, barrett)? {
            let parts = match method {
                | SplitMethod::CantorZassenhaus => {
                    let mut parts = Vec::new();
                    edf_u64(&gd, d, barrett, &mut rng, &mut parts)?;
                    parts
                },
                | SplitMethod::Berlekamp => berlekamp_split_u64(&gd, barrett)?,
            };
            for q in parts {
                out.push((q, m));
            }
        }
    }
    out.sort_by(|a, b| {
        a.0.len()
            .cmp(&b.0.len())
            .then_with(|| a.0.cmp(&b.0))
            .then(a.1.cmp(&b.1))
    });
    Some((lead, out))
}

/// Complete factorisation (`BigInt` reference fallback for `p >= 2^62`).
#[allow(clippy::single_call_fn)]
fn factor_bigint(
    f: &Poly,
    p: &BigInt,
    method: SplitMethod,
) -> Option<(BigInt, Vec<(Poly, usize)>)> {
    let lead = f.first()?.clone();
    let g = monic(f, p)?;
    let mut rng = Rng(0x5eed_1234_abcd_9876);
    let mut out: Vec<(Poly, usize)> = Vec::new();
    for (sq, m) in squarefree_bigint(&g, p)? {
        for (gd, d) in ddf_bigint(&sq, p)? {
            let parts = match method {
                | SplitMethod::CantorZassenhaus => {
                    let mut parts = Vec::new();
                    edf_bigint(&gd, d, p, &mut rng, &mut parts)?;
                    parts
                },
                | SplitMethod::Berlekamp => berlekamp_split_bigint(&gd, p)?,
            };
            for q in parts {
                out.push((q, m));
            }
        }
    }
    out.sort_by(|a, b| {
        a.0.len()
            .cmp(&b.0.len())
            .then_with(|| a.0.cmp(&b.0))
            .then(a.1.cmp(&b.1))
    });
    Some((lead, out))
}

/// The full factorisation: the leading coefficient and `(factor,
/// multiplicity)` pairs, ordered by degree.
fn factor(
    f: &Poly,
    p: &BigInt,
    method: SplitMethod,
) -> Option<(BigInt, Vec<(Poly, usize)>)> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let (lead, out) = factor_u64(&fu, barrett, method)?;
        return Some((
            BigInt::from(lead),
            out.into_iter()
                .map(|(q, m)| (poly_from_u64(&q), m))
                .collect(),
        ));
    }
    factor_bigint(f, p, method)
}

// ---------------- Berlekamp ----------------

/// A basis of `{v : M v = 0}` over GF(p) for an `n x n` matrix (`BigInt` reference).
#[allow(
    clippy::needless_range_loop,
    clippy::single_call_fn,
    reason = "rows and columns are indexed together"
)]
fn nullspace_bigint(
    mut m: Vec<Vec<BigInt>>,
    p: &BigInt,
) -> Option<Vec<Vec<BigInt>>> {
    let n = m.len();
    let mut pivots = Vec::new();
    let mut row = 0;
    for col in 0..n {
        let Some(found) = (row..n).find(|&r| !m[r][col].is_zero()) else {
            continue;
        };
        m.swap(row, found);
        let inv = mod_inverse(&m[row][col], p)?;
        for c in 0..n {
            m[row][c] = modulo(&(&m[row][c] * &inv), p);
        }
        for r in 0..n {
            if r != row && !m[r][col].is_zero() {
                let factor = m[r][col].clone();
                for c in 0..n {
                    let delta = &factor * &m[row][c];
                    m[r][c] = modulo(&(&m[r][c] - delta), p);
                }
            }
        }
        pivots.push(col);
        row += 1;
        if row == n {
            break;
        }
    }
    let mut basis = Vec::new();
    for free in (0..n).filter(|c| !pivots.contains(c)) {
        let mut v = vec![BigInt::zero(); n];
        v[free] = BigInt::one();
        for (r, &pc) in pivots.iter().enumerate() {
            v[pc] = modulo(&-&m[r][free], p);
        }
        basis.push(v);
    }
    Some(basis)
}

#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "chunks of stride n keep col..n within each row and pivot_tail has length n - col"
)]
fn eliminate_rows_u64(
    rows: &mut [u64],
    n: usize,
    col: usize,
    pivot_tail: &[u64],
    barrett: BarrettMod,
    mont: Option<MontgomeryPrimeField>,
) {
    if rows.is_empty() {
        return;
    }
    let p = barrett.modulus();
    par::for_chunks_mut(rows, n, 64, |_, target_row| {
        let factor = target_row[col];
        if factor == 0 {
            return;
        }
        let neg_factor = barrett.neg(factor);
        let target_tail = &mut target_row[col..n];
        if target_tail.len() >= 256 || mont.is_none() {
            slice_axpy_mod(neg_factor, pivot_tail, target_tail, p);
        } else if let Some(mont_field) = mont {
            let nf_mont = mont_field.to_mont(neg_factor);
            for (dst, &pv) in target_tail.iter_mut().zip(pivot_tail) {
                if pv != 0 {
                    let prod = mont_field.redc(u128::from(pv) * u128::from(nf_mont));
                    *dst = mont_field.add_mont(*dst, prod);
                }
            }
        }
    });
}

/// Nullspace basis of a flat row-major `n x n` matrix over `GF(p)` using
/// [`MontgomeryPrimeField`], [`slice_axpy_mod`], [`global_pool`], and [`par`].
#[allow(
    clippy::single_call_fn,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "m has length n * n and all row/col indices stay in 0..n"
)]
fn nullspace_u64(
    mut m: Vec<u64>,
    n: usize,
    barrett: BarrettMod,
) -> Option<Vec<Vec<u64>>> {
    let p = barrett.modulus();
    let mont = MontgomeryPrimeField::new(p);
    let pool = global_pool();
    let mut pivots = pool.scratch_usize(n);
    let mut is_pivot = pool.scratch_bytes(n);
    is_pivot.resize(n, 0);

    let mut row = 0usize;
    for col in 0..n {
        let Some(found) = (row..n).find(|&r| m[r * n + col] != 0) else {
            continue;
        };
        if found != row {
            let (top, bot) = m.split_at_mut(found * n);
            top[row * n..(row + 1) * n].swap_with_slice(&mut bot[..n]);
        }
        let inv = barrett.inv(m[row * n + col])?;
        if inv != 1 {
            if let Some(mont_field) = mont {
                let inv_mont = mont_field.to_mont(inv);
                for slot in &mut m[row * n + col..(row + 1) * n] {
                    *slot = mont_field.redc(u128::from(*slot) * u128::from(inv_mont));
                }
            } else {
                for slot in &mut m[row * n + col..(row + 1) * n] {
                    *slot = barrett.mul(*slot, inv);
                }
            }
        }
        let (before_pivot, pivot_and_after) = m.split_at_mut(row * n);
        let (pivot_row, after_pivot) = pivot_and_after.split_at_mut(n);
        let pivot_tail = &pivot_row[col..n];
        eliminate_rows_u64(before_pivot, n, col, pivot_tail, barrett, mont);
        eliminate_rows_u64(after_pivot, n, col, pivot_tail, barrett, mont);

        pivots.push(col);
        is_pivot[col] = 1;
        row += 1;
        if row == n {
            break;
        }
    }

    let mut basis = Vec::with_capacity(n - pivots.len());
    for free in 0..n {
        if is_pivot[free] != 0 {
            continue;
        }
        let mut v = vec![0u64; n];
        v[free] = 1;
        for (r, &pc) in pivots.iter().enumerate() {
            v[pc] = barrett.neg(m[r * n + free]);
        }
        basis.push(v);
    }
    Some(basis)
}

/// Ascending coefficient vector of length `n`.
fn ascending(
    f: &Poly,
    n: usize,
) -> Vec<BigInt> {
    let mut v: Vec<BigInt> = f.iter().rev().cloned().collect();
    v.resize(n, BigInt::zero());
    v
}

/// Berlekamp's algorithm for a monic squarefree polynomial (`BigInt` reference).
#[allow(clippy::single_call_fn)]
fn berlekamp_split_bigint(
    g: &Poly,
    p: &BigInt,
) -> Option<Vec<Poly>> {
    let n = g.len().checked_sub(1)?;
    if n <= 1 {
        return Some(vec![g.clone()]);
    }
    let pu = p.to_usize().filter(|&v| v <= 4096)?;
    let xp = powmod_bigint(&x_poly(), p, g, p)?;
    // Row i of Q: x^(p i) mod g. The sought h = sum v_i x^i satisfy
    // v (Q - I) = 0.
    let mut q = vec![ascending(&reduce(&one(), g, p)?, n)];
    let mut current = reduce(&one(), g, p)?;
    for _ in 1..n {
        current = reduce(&pmul(&current, &xp, p), g, p)?;
        q.push(ascending(&current, n));
    }
    let mut system = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let delta = if i == j {
                BigInt::one()
            } else {
                BigInt::zero()
            };
            system[j][i] = modulo(&(&q[i][j] - delta), p);
        }
    }
    let basis = nullspace_bigint(system, p)?;
    let count = basis.len();
    let mut factors = vec![g.clone()];
    for v in &basis {
        if factors.len() == count {
            break;
        }
        let h: Poly = norm(v.iter().rev().cloned().collect(), p);
        if h.len() <= 1 {
            continue;
        }
        let mut next = Vec::new();
        for f in &factors {
            let mut rest = f.clone();
            for s in 0..pu {
                let shifted = psub(&h, &vec![BigInt::from(s)], p);
                let d = gcd(&rest, &shifted, p)?;
                if !is_constant(&d) {
                    rest = quotient(&rest, &d, p)?;
                    next.push(d);
                }
            }
            if !is_constant(&rest) {
                next.push(rest);
            }
        }
        factors = next;
    }
    (factors.len() == count).then_some(factors)
}

/// Berlekamp's algorithm for a monic squarefree polynomial over `u64` residues.
#[allow(
    clippy::single_call_fn,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "n = g.len() - 1 >= 2, current.len() <= n, and system has n * n entries"
)]
fn berlekamp_split_u64(
    g: &[u64],
    barrett: BarrettMod,
) -> Option<Vec<PolyU64>> {
    let n = g.len().checked_sub(1)?;
    if n <= 1 {
        return Some(vec![g.to_vec()]);
    }
    let p_u64 = barrett.modulus();
    if p_u64 > 4096 {
        return None;
    }
    let x = [1u64, 0u64];
    let xp = powmod_u64_limbs(&x, &[p_u64], g, barrett)?;
    let mut system = vec![0u64; n * n];
    let mut current = reduce_u64(&[1u64], g, barrett)?;
    for i in 1..n {
        current = reduce_u64(&pmul_u64(&current, &xp, barrett), g, barrett)?;
        let clen = current.len();
        for j in 0..clen {
            system[j * n + i] = current[clen - 1 - j];
        }
        system[i * n + i] = barrett.sub(system[i * n + i], 1);
    }

    let basis = nullspace_u64(system, n, barrett)?;
    let count = basis.len();
    let mut factors = vec![g.to_vec()];
    for v in &basis {
        if factors.len() == count {
            break;
        }
        let mut h: PolyU64 = v.iter().rev().copied().collect();
        norm_u64_in_place(&mut h, barrett);
        if h.len() <= 1 {
            continue;
        }
        let mut next = Vec::new();
        for f in &factors {
            let mut rest = f.clone();
            for s in 0..p_u64 {
                if rest.len() <= 1 {
                    break;
                }
                let mut shifted = h.clone();
                let last = shifted.len() - 1;
                shifted[last] = barrett.sub(shifted[last], s);
                let d = pgcd_u64(&rest, &shifted, barrett)?;
                if !is_constant_u64(&d) {
                    rest = quotient_u64(&rest, &d, barrett)?;
                    next.push(d);
                }
            }
            if !is_constant_u64(&rest) {
                next.push(rest);
            }
        }
        factors = next;
    }
    (factors.len() == count).then_some(factors)
}

// ---------------- operators ----------------

fn factorisation_value(
    lead: BigInt,
    factors: Vec<(Poly, usize)>,
) -> V {
    let mut items = vec![V::Int(lead)];
    items.extend(
        factors
            .into_iter()
            .map(|(g, m)| V::List(vec![poly_value(g), V::uint(m)])),
    );
    V::List(items)
}

fn pairs_value(pairs: Vec<(Poly, usize)>) -> V {
    V::List(
        pairs
            .into_iter()
            .map(|(g, m)| V::List(vec![poly_value(g), V::uint(m)]))
            .collect(),
    )
}

fn read_prime_poly(
    cx: &Cx<'_>,
    f: NodeId,
    p: NodeId,
) -> Option<(Poly, BigInt)> {
    let p = prime(cx, p)?;
    Some((read_poly(cx, f, &p)?, p))
}

fn gfp_monic_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    Some(poly_value(monic(&f, &p)?))
}

fn gfp_derivative_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let p = big(cx.graph, *p).filter(|p| *p > BigInt::one())?;
    let f = read_poly(cx, *f, &p)?;
    Some(poly_value(derivative(&f, &p)))
}

fn gfp_powmod_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, e, m, p] = a else { return None };
    let p = big(cx.graph, *p).filter(|p| *p > BigInt::one())?;
    let (f, m) = (read_poly(cx, *f, &p)?, read_poly(cx, *m, &p)?);
    let e = big(cx.graph, *e)?;
    if m.is_empty() {
        return None;
    }
    Some(poly_value(powmod(&f, &e, &m, &p)?))
}

fn gfp_invmod_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, m, p] = a else { return None };
    let p = big(cx.graph, *p).filter(|p| *p > BigInt::one())?;
    let (f, m) = (read_poly(cx, *f, &p)?, read_poly(cx, *m, &p)?);
    if m.is_empty() {
        return None;
    }
    Some(poly_value(xinv(&f, &m, &p)?))
}

fn gfp_squarefree_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    if let Some(barrett) = fast_barrett(&p) {
        let fu = monic_u64(&poly_to_u64(&f, &p, barrett), barrett)?;
        let pairs = squarefree_u64(&fu, barrett)?
            .into_iter()
            .map(|(g, m)| (poly_from_u64(&g), m))
            .collect();
        return Some(pairs_value(pairs));
    }
    Some(pairs_value(squarefree_bigint(&monic(&f, &p)?, &p)?))
}

fn gfp_ddf_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    if let Some(barrett) = fast_barrett(&p) {
        let fu = monic_u64(&poly_to_u64(&f, &p, barrett), barrett)?;
        let df = derivative_u64(&fu, barrett);
        if !is_constant_u64(&pgcd_u64(&fu, &df, barrett)?) || is_constant_u64(&fu) {
            return None;
        }
        let pairs = ddf_u64(&fu, barrett)?
            .into_iter()
            .map(|(g, d)| (poly_from_u64(&g), d))
            .collect();
        return Some(pairs_value(pairs));
    }
    let f = monic(&f, &p)?;
    // The input must be squarefree.
    if !is_constant(&gcd(&f, &derivative(&f, &p), &p)?) || is_constant(&f) {
        return None;
    }
    Some(pairs_value(ddf_bigint(&f, &p)?))
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "fu.len() > 1 and q.len() >= d + 1 >= 2 inside length checks"
)]
fn gfp_edf_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, d, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let d = usize::try_from(big(cx.graph, *d)?)
        .ok()
        .filter(|&d| d >= 1)?;
    if let Some(barrett) = fast_barrett(&p) {
        let fu = monic_u64(&poly_to_u64(&f, &p, barrett), barrett)?;
        if fu.len() <= 1
            || (fu.len() - 1) % d != 0
            || !is_constant_u64(&pgcd_u64(&fu, &derivative_u64(&fu, barrett), barrett)?)
        {
            return None;
        }
        let mut parts = Vec::new();
        edf_u64(&fu, d, barrett, &mut Rng(0x5eed_1234_abcd_9876), &mut parts)?;
        let invalid = par::map_slice(&parts, 4, |q| {
            q.len() - 1 != d || is_irreducible_u64(q, barrett) != Some(true)
        });
        if invalid.into_iter().any(|bad| bad) {
            return None;
        }
        parts.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        return Some(V::List(
            parts
                .into_iter()
                .map(|q| poly_value(poly_from_u64(&q)))
                .collect(),
        ));
    }
    let f = monic(&f, &p)?;
    if f.len() <= 1 || (f.len() - 1) % d != 0 || !is_constant(&gcd(&f, &derivative(&f, &p), &p)?) {
        return None;
    }
    let mut parts = Vec::new();
    edf_bigint(&f, d, &p, &mut Rng(0x5eed_1234_abcd_9876), &mut parts)?;
    // Verify the premise: every part must be irreducible of degree d.
    if parts
        .iter()
        .any(|q| q.len() - 1 != d || is_irreducible_bigint(q, &p) != Some(true))
    {
        return None;
    }
    parts.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    Some(V::List(parts.into_iter().map(poly_value).collect()))
}

fn gfp_factor_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let (lead, factors) = factor(&f, &p, SplitMethod::CantorZassenhaus)?;
    Some(factorisation_value(lead, factors))
}

fn gfp_berlekamp_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let (lead, factors) = factor(&f, &p, SplitMethod::Berlekamp)?;
    Some(factorisation_value(lead, factors))
}

// ---------------- expressions in x ----------------

/// The coefficients (highest degree first) of an expression in `x` with
/// rational coefficients whose denominators are invertible mod `p`.
fn poly_of_expression(
    graph: &mut Graph,
    expr: NodeId,
    x: NodeId,
    p: &BigInt,
) -> Option<Poly> {
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let poly = from_term(graph, &mut gens, expr, Limits::default())?;
    if gens.len() != 1 {
        return None;
    }
    let mut coefficients = Vec::new();
    for c in poly.univariate_in(gx)? {
        let r: BigRational = c.to_rational()?;
        let value = modulo(r.numer(), p) * mod_inverse(r.denom(), p)?;
        coefficients.push(modulo(&value, p));
    }
    coefficients.reverse();
    Some(norm(coefficients, p))
}

/// `sum c_i x^i` from a highest-first coefficient list.
fn expression_of_poly(
    graph: &mut Graph,
    f: &Poly,
    x: NodeId,
) -> NodeId {
    let n = f.len();
    let mut terms = Vec::new();
    for (j, c) in f.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let degree = n - 1 - j;
        let coefficient = graph.num(Number::Int(c.clone()));
        let term = match degree {
            | 0 => coefficient,
            | 1 => prod(graph, &[coefficient, x]),
            | _ => {
                let e = graph.int(i64::try_from(degree).unwrap_or(i64::MAX));
                let power = pow(graph, x, e);
                prod(graph, &[coefficient, power])
            },
        };
        terms.push(term);
    }
    terms.reverse();
    sum(graph, &terms)
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    FactorMod,
    IsIrreducibleMod,
}

struct ModRequest {
    op: OpId,
    request: Request,
}

impl Kernel for ModRequest {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, x, p] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(p) = prime(cx, p) else {
            return Outcome::Pass;
        };
        let Some(f) = best(cx.graph, f) else {
            return Outcome::Pass;
        };
        let as_list = items(cx.graph, f).is_some();
        let poly = if as_list {
            read_poly(cx, f, &p)
        } else {
            let Some(x) = best(cx.graph, x) else {
                return Outcome::Pass;
            };
            poly_of_expression(cx.graph, f, x, &p)
        };
        let Some(poly) = poly else {
            return Outcome::Pass;
        };
        if self.request == Request::IsIrreducibleMod {
            return is_irreducible(&poly, &p).map_or(Outcome::Pass, |b| {
                Outcome::Equal(V::Bool(b).build(cx.graph))
            });
        }
        let Some((lead, factors)) = factor(&poly, &p, SplitMethod::CantorZassenhaus) else {
            return Outcome::Pass;
        };
        if as_list {
            return Outcome::Pinned(factorisation_value(lead, factors).build(cx.graph));
        }
        let Some(x) = best(cx.graph, x) else {
            return Outcome::Pass;
        };
        let mut parts = vec![cx.graph.num(Number::Int(lead))];
        for (g, m) in factors {
            let base = expression_of_poly(cx.graph, &g, x);
            let term = if m == 1 {
                base
            } else {
                let e = cx.graph.int(i64::try_from(m).unwrap_or(i64::MAX));
                cx.graph.node(core::POW, &[base, e])
            };
            parts.push(term);
        }
        Outcome::Pinned(prod(cx.graph, &parts))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def_doc(
        i,
        "gfp_monic",
        Arity::Fixed(2),
        gfp_monic_op,
        "Divides polynomial `f` by its leading coefficient modulo prime `p` to produce the monic polynomial.",
    )?;
    def_doc(
        i,
        "gfp_derivative",
        Arity::Fixed(2),
        gfp_derivative_op,
        "Computes the formal derivative `f'` of polynomial `f` over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_powmod",
        Arity::Fixed(4),
        gfp_powmod_op,
        "Computes `f^e mod m` for polynomial `f`, non-negative integer `e`, and modulus polynomial `m` over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_invmod",
        Arity::Fixed(3),
        gfp_invmod_op,
        "Computes the inverse of polynomial `f` modulo `m` over GF(p) when `gcd(f, m) = 1`.",
    )?;
    def_doc(
        i,
        "gfp_squarefree",
        Arity::Fixed(2),
        gfp_squarefree_op,
        "Computes the squarefree factorization of polynomial `f` over GF(p) as `list(list(g1, m1), ...)` using Yun's algorithm.",
    )?;
    def_doc(
        i,
        "gfp_ddf",
        Arity::Fixed(2),
        gfp_ddf_op,
        "Computes the distinct-degree factorization of a squarefree polynomial `f` over GF(p) as `list(list(g_d, d), ...)`.",
    )?;
    def_doc(
        i,
        "gfp_edf",
        Arity::Fixed(3),
        gfp_edf_op,
        "Computes the equal-degree factorization (Cantor-Zassenhaus) of a squarefree polynomial `f` of degree `k*d` over GF(p) into irreducible factors of degree `d`.",
    )?;
    def_doc(
        i,
        "gfp_factor",
        Arity::Fixed(2),
        gfp_factor_op,
        "Computes the complete factorization of polynomial `f` over GF(p) as `list(lc, list(g1, m1), ...)`, ordered by degree.",
    )?;
    def_doc(
        i,
        "gfp_berlekamp",
        Arity::Fixed(2),
        gfp_berlekamp_op,
        "Computes the complete factorization of polynomial `f` over a small prime field GF(p) using Berlekamp's algorithm.",
    )?;
    for (name, request, doc) in [
        (
            "factor_mod",
            Request::FactorMod,
            "Factors a polynomial expression `f` in variable `x` modulo prime `p` into a product `lc * g1^m1 * ...` with coefficients in `[0, p)`.",
        ),
        (
            "is_irreducible_mod",
            Request::IsIrreducibleMod,
            "Tests whether a polynomial expression `f` in variable `x` is irreducible modulo prime `p`.",
        ),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc(doc))?;
        // Factoring rewrites the polynomial expression, so it keeps the
        // dimension of the polynomial; irreducibility is a predicate.
        i.dim(
            op,
            if name == "factor_mod" {
                crate::graph::dim::DimRule::Same(0)
            } else {
                crate::graph::dim::DimRule::NonMetric
            },
        );
        i.kernel(
            &format!("discrete/{name}"),
            Tier::Reduce,
            ModRequest { op, request },
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::Lcg;
    use super::super::test_util::s;

    /// A parsed `list(...)` of integers and lists.
    #[derive(Debug, Clone, PartialEq)]
    enum Tree {
        Num(i64),
        List(Vec<Self>),
    }

    fn parse(text: &str) -> Tree {
        fn go(
            chars: &[char],
            at: &mut usize,
        ) -> Tree {
            let rest: String = chars[*at..].iter().collect();
            if rest.starts_with("list(") {
                *at += 5;
                let mut items = Vec::new();
                loop {
                    while chars[*at] == ' ' || chars[*at] == ',' {
                        *at += 1;
                    }
                    if chars[*at] == ')' {
                        *at += 1;
                        return Tree::List(items);
                    }
                    items.push(go(chars, at));
                }
            }
            let start = *at;
            while *at < chars.len() && (chars[*at].is_ascii_digit() || chars[*at] == '-') {
                *at += 1;
            }
            let number: String = chars[start..*at].iter().collect();
            Tree::Num(
                number
                    .parse()
                    .unwrap_or_else(|_| panic!("cannot parse `{number}`")),
            )
        }
        go(&text.chars().collect::<Vec<_>>(), &mut 0)
    }

    fn coefficients(tree: &Tree) -> Vec<i64> {
        match tree {
            | Tree::List(items) => {
                items
                    .iter()
                    .map(|t| {
                        if let Tree::Num(n) = t {
                            *n
                        } else {
                            panic!("nested")
                        }
                    })
                    .collect()
            },
            | Tree::Num(_) => panic!("not a polynomial"),
        }
    }

    fn poly(c: &[i64]) -> String {
        format!(
            "list({})",
            c.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    fn mul(
        f: &[i64],
        g: &[i64],
        p: i64,
    ) -> Vec<i64> {
        coefficients(&parse(&s(&format!(
            "gfp_mul({}, {}, {p})",
            poly(f),
            poly(g)
        ))))
    }

    /// `(lead, [(factor, multiplicity)])` from a `gfp_factor` answer.
    fn factorisation(text: &str) -> (i64, Vec<(Vec<i64>, usize)>) {
        let Tree::List(items) = parse(text) else {
            panic!("{text}")
        };
        let Tree::Num(lead) = items[0] else {
            panic!("{text}")
        };
        let factors = items[1..]
            .iter()
            .map(|t| {
                let Tree::List(pair) = t else {
                    panic!("{text}")
                };
                let Tree::Num(m) = pair[1] else {
                    panic!("{text}")
                };
                (coefficients(&pair[0]), usize::try_from(m).unwrap_or(0))
            })
            .collect();
        (lead, factors)
    }

    /// Whether `g` (monic) divides `f` over GF(p), by schoolbook division.
    fn divides(
        f: &[i64],
        g: &[i64],
        p: i64,
    ) -> bool {
        let mut r: Vec<i64> = f.iter().map(|c| c.rem_euclid(p)).collect();
        while r.len() >= g.len() {
            let lead = r[0];
            for (k, gc) in g.iter().enumerate() {
                r[k] = (r[k] - lead * gc).rem_euclid(p);
            }
            r.remove(0);
        }
        r.iter().all(|&c| c == 0)
    }

    /// Brute-force irreducibility: no monic divisor of degree `1..=n/2`.
    fn brute_irreducible(
        f: &[i64],
        p: i64,
    ) -> bool {
        let n = f.len() - 1;
        for d in 1..=n / 2 {
            let count = p.pow(u32::try_from(d).unwrap_or(1));
            for code in 0..count {
                let mut g = vec![1];
                let mut c = code;
                let mut tail = Vec::new();
                for _ in 0..d {
                    tail.push(c % p);
                    c /= p;
                }
                tail.reverse();
                g.extend(tail);
                if divides(f, &g, p) {
                    return false;
                }
            }
        }
        true
    }

    fn rebuild(
        lead: i64,
        factors: &[(Vec<i64>, usize)],
        p: i64,
    ) -> Vec<i64> {
        let mut acc = vec![lead];
        for (g, m) in factors {
            for _ in 0..*m {
                acc = mul(&acc, g, p);
            }
        }
        acc
    }

    #[test]
    fn squarefree_decomposition() {
        // (x + 1)^2 (x^2 + x + 1) over GF(2)
        let f = mul(&mul(&[1, 1], &[1, 1], 2), &[1, 1, 1], 2);
        assert_eq!(
            s(&format!("gfp_squarefree({}, 2)", poly(&f))),
            "list(list(list(1, 1, 1), 1), list(list(1, 1), 2))"
        );
        // x^p - x^0 in characteristic p: x^3 + 2 = (x + 2)^3 over GF(3)
        assert_eq!(
            s("gfp_squarefree(list(1, 0, 0, 2), 3)"),
            "list(list(list(1, 2), 3))"
        );
        // x^6 + 1 = (x^2 + 1)^3 over GF(3)
        assert_eq!(
            s("gfp_squarefree(list(1, 0, 0, 0, 0, 0, 1), 3)"),
            "list(list(list(1, 0, 1), 3))"
        );
        // squarefree input is its own decomposition (made monic)
        assert_eq!(
            s("gfp_squarefree(list(2, 0, 2), 5)"),
            "list(list(list(1, 0, 1), 1))"
        );
        assert_eq!(s("gfp_squarefree(list(3), 5)"), "list()");
        // mixed multiplicities 1, 2, 3 and p-th powers
        let f = mul(
            &mul(&[1, 2], &mul(&[1, 3], &[1, 3], 7), 7),
            &mul(&[1, 5], &mul(&[1, 5], &[1, 5], 7), 7),
            7,
        );
        assert_eq!(
            s(&format!("gfp_squarefree({}, 7)", poly(&f))),
            "list(list(list(1, 2), 1), list(list(1, 3), 2), list(list(1, 5), 3))"
        );
    }

    #[test]
    fn distinct_degree_factorisation() {
        // x^8 - x = x^8 + x over GF(2): all irreducibles of degree 1 and 3
        // x (x+1) (x^3 + x + 1)(x^3 + x^2 + 1)
        assert_eq!(
            s("gfp_ddf(list(1, 0, 0, 0, 0, 0, 0, 1, 0), 2)"),
            "list(list(list(1, 1, 0), 1), list(list(1, 1, 1, 1, 1, 1, 1), 3))"
        );
        // x^2 + 1 over GF(5) = (x + 2)(x + 3): one class of degree one
        assert_eq!(
            s("gfp_ddf(list(1, 0, 1), 5)"),
            "list(list(list(1, 0, 1), 1))"
        );
        // x^2 + 1 over GF(3) is irreducible
        assert_eq!(
            s("gfp_ddf(list(1, 0, 1), 3)"),
            "list(list(list(1, 0, 1), 2))"
        );
        // not squarefree: stays unreduced
        assert_eq!(s("gfp_ddf(list(1, 2, 1), 5)"), "gfp_ddf(list(1, 2, 1), 5)");
    }

    #[test]
    fn equal_degree_factorisation() {
        // the two cubics of GF(2)
        assert_eq!(
            s("gfp_edf(list(1, 1, 1, 1, 1, 1, 1), 3, 2)"),
            "list(list(1, 0, 1, 1), list(1, 1, 0, 1))"
        );
        // x^2 + 1 over GF(5): (x + 2)(x + 3)
        assert_eq!(
            s("gfp_edf(list(1, 0, 1), 1, 5)"),
            "list(list(1, 2), list(1, 3))"
        );
        // wrong premise: x^2 + 1 over GF(3) has no factor of degree 1
        assert_eq!(
            s("gfp_edf(list(1, 0, 1), 1, 3)"),
            "gfp_edf(list(1, 0, 1), 1, 3)"
        );
        // large prime: (x - 3)(x - 5)(x - 11) mod 1000003
        let p = 1_000_003_i64;
        let f = mul(&mul(&[1, p - 3], &[1, p - 5], p), &[1, p - 11], p);
        assert_eq!(
            s(&format!("gfp_edf({}, 1, {p})", poly(&f))),
            format!(
                "list(list(1, {}), list(1, {}), list(1, {}))",
                p - 11,
                p - 5,
                p - 3
            )
        );
    }

    #[test]
    fn complete_factorisation() {
        assert_eq!(
            s("gfp_factor(list(1, 0, 0, 0, 1), 2)"),
            "list(1, list(list(1, 1), 4))"
        );
        assert_eq!(
            s("gfp_factor(list(3, 0, 3), 5)"),
            "list(3, list(list(1, 2), 1), list(list(1, 3), 1))"
        );
        assert_eq!(s("gfp_factor(list(4), 7)"), "list(4)");
        assert_eq!(s("gfp_factor(list(), 7)"), "gfp_factor(list(), 7)");
        // composite modulus is rejected
        assert_eq!(
            s("gfp_factor(list(1, 0, 1), 4)"),
            "gfp_factor(list(1, 0, 1), 4)"
        );
        // x^4 + 1 over GF(3) = (x^2 + x + 2)(x^2 + 2x + 2)
        assert_eq!(
            s("gfp_factor(list(1, 0, 0, 0, 1), 3)"),
            "list(1, list(list(1, 1, 2), 1), list(list(1, 2, 2), 1))"
        );
        assert_eq!(
            s("gfp_berlekamp(list(1, 0, 0, 0, 1), 3)"),
            "list(1, list(list(1, 1, 2), 1), list(list(1, 2, 2), 1))"
        );
        // x^8 - 1 over GF(17) splits completely into 8 linear factors
        let (lead, factors) = factorisation(&s("gfp_factor(list(1, 0, 0, 0, 0, 0, 0, 0, 16), 17)"));
        assert_eq!((lead, factors.len()), (1, 8));
        assert!(factors.iter().all(|(g, m)| g.len() == 2 && *m == 1));
    }

    /// Random polynomials: the product of the factors is the polynomial,
    /// every factor is irreducible, and both algorithms agree.
    #[test]
    fn random_factorisations_are_consistent() {
        let mut rng = Lcg(2024);
        for p in [2_i64, 3, 5, 7, 11] {
            for _ in 0..6 {
                let degree = 2 + usize::try_from(rng.next(9)).unwrap_or(0);
                let mut f: Vec<i64> = (0..=degree)
                    .map(|_| i64::try_from(rng.next(u64::try_from(p).unwrap_or(2))).unwrap_or(0))
                    .collect();
                if f[0] == 0 {
                    f[0] = 1;
                }
                let want = coefficients(&parse(&s(&format!("gfp_norm({}, {p})", poly(&f)))));
                let (lead, factors) = factorisation(&s(&format!("gfp_factor({}, {p})", poly(&f))));
                assert_eq!(rebuild(lead, &factors, p), want, "{f:?} mod {p}");
                for (g, _) in &factors {
                    assert!(
                        brute_irreducible(g, p),
                        "{g:?} reducible mod {p} (from {f:?})"
                    );
                    assert_eq!(s(&format!("gfp_isirreducible({}, {p})", poly(g))), "true");
                }
                let (lead2, factors2) =
                    factorisation(&s(&format!("gfp_berlekamp({}, {p})", poly(&f))));
                assert_eq!((lead, &factors), (lead2, &factors2), "{f:?} mod {p}");
            }
        }
    }

    #[test]
    fn irreducibility_agrees_with_brute_force() {
        for p in [2_i64, 3, 5] {
            for degree in 1..=6_u32 {
                let count = p.pow(degree);
                for code in 0..count.min(300) {
                    let mut f = vec![1];
                    let mut c = code;
                    let mut tail = Vec::new();
                    for _ in 0..degree {
                        tail.push(c % p);
                        c /= p;
                    }
                    tail.reverse();
                    f.extend(tail);
                    assert_eq!(
                        s(&format!("gfp_isirreducible({}, {p})", poly(&f))),
                        brute_irreducible(&f, p).to_string(),
                        "{f:?} mod {p}"
                    );
                }
            }
        }
        assert_eq!(s("gfp_isirreducible(list(3), 5)"), "false");
        assert_eq!(s("gfp_isirreducible(list(2, 1), 5)"), "true");
        // a degree-12 irreducible over GF(2): x^12 + x^6 + x^4 + x + 1
        assert_eq!(
            s("gfp_isirreducible(list(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1), 2)"),
            "true"
        );
    }

    #[test]
    fn powmod_invmod_and_derivative() {
        // x^10 mod (x^2 + 1) over GF(5): x^2 = -1, x^10 = -1 = 4
        assert_eq!(s("gfp_powmod(list(1, 0), 10, list(1, 0, 1), 5)"), "list(4)");
        assert_eq!(s("gfp_powmod(list(1, 1), 0, list(1, 0, 1), 5)"), "list(1)");
        // Frobenius: (x + 1)^7 = x^7 + 1 over GF(7)
        assert_eq!(
            s("gfp_powmod(list(1, 1), 7, list(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), 7)"),
            "list(1, 0, 0, 0, 0, 0, 0, 1)"
        );
        assert_eq!(
            s("gfp_powmod(list(1, 1), -1, list(1, 0, 1), 5)"),
            "gfp_powmod(list(1, 1), -1, list(1, 0, 1), 5)"
        );
        // x (x + 1) = x^2 + x = 1 + x - 1... check f * f^-1 = 1 mod m
        let m = "list(1, 0, 1, 1)"; // x^3 + x + 1 over GF(2): a field
        for f in ["list(1, 0)", "list(1, 1)", "list(1, 0, 1)", "list(1, 1, 1)"] {
            let inv = s(&format!("gfp_invmod({f}, {m}, 2)"));
            let product = s(&format!("gfp_mul({f}, {inv}, 2)"));
            assert_eq!(
                s(&format!("gfp_divmod({product}, {m}, 2)"))
                    .rsplit_once(", ")
                    .map(|(_, r)| r.to_string()),
                Some("list(1))".to_string()),
                "{f}"
            );
        }
        assert_eq!(
            s("gfp_invmod(list(1, 1), list(1, 0, 1), 2)"),
            "gfp_invmod(list(1, 1), list(1, 0, 1), 2)"
        );
        assert_eq!(
            s("gfp_derivative(list(1, 0, 0, 5, 3), 7)"),
            "list(4, 0, 0, 5)"
        );
        assert_eq!(s("gfp_derivative(list(1, 0, 0), 2)"), "list()");
        assert_eq!(s("gfp_monic(list(2, 3), 5)"), "list(1, 4)");
    }

    #[test]
    fn expressions_modulo_p() {
        assert_eq!(s("factor_mod(x^4 + 1, x, 2)"), "(x + 1)^4");
        assert_eq!(s("factor_mod(x^2 + 1, x, 5)"), "(x + 2)*(x + 3)");
        assert_eq!(s("factor_mod(x^2 + 1, x, 3)"), "x^2 + 1");
        assert_eq!(s("factor_mod(3*x^2 + 3, x, 5)"), "3*(x + 2)*(x + 3)");
        assert_eq!(s("factor_mod(x^2 - 1, x, 7)"), "(x + 1)*(x + 6)");
        assert_eq!(s("factor_mod(x^3 + 2*x + 1, x, 3)"), "x^3 + 2*x + 1");
        assert_eq!(s("factor_mod(x^3 - x, x, 3)"), "x*(x + 1)*(x + 2)");
        // rational coefficients with denominators prime to p
        assert_eq!(s("factor_mod(x^2/2 + 1/2, x, 5)"), "3*(x + 2)*(x + 3)");
        // the list form
        assert_eq!(
            s("factor_mod(list(1, 0, 0, 0, 1), x, 2)"),
            "list(1, list(list(1, 1), 4))"
        );
        // irreducibility
        assert_eq!(s("is_irreducible_mod(x^2 + x + 1, x, 2)"), "true");
        assert_eq!(s("is_irreducible_mod(x^2 + 1, x, 2)"), "false");
        assert_eq!(s("is_irreducible_mod(list(1, 1, 1), x, 2)"), "true");
        // other variables or non-prime moduli stay unreduced
        for src in ["factor_mod(x^2 + y, x, 5)", "factor_mod(x^2 + 1, x, 6)"] {
            let (text, reduced) =
                crate::rules::testing::reduce_with(&[crate::rules::discrete::discrete()], src, &[]);
            assert_eq!((text.as_str(), reduced), (src, false));
        }
    }

    #[test]
    fn differential_parity_u64_vs_bigint_across_primes() {
        use num_bigint::BigInt;
        use num_traits::Zero;

        use super::super::finite_field::Poly;
        use super::super::finite_field::fast_barrett;
        use super::super::finite_field::monic_u64;
        use super::super::finite_field::norm;
        use super::super::finite_field::pmul;
        use super::super::finite_field::poly_from_u64;
        use super::super::finite_field::poly_to_u64;
        use super::Rng;
        use super::SplitMethod;
        use super::berlekamp_split_bigint;
        use super::berlekamp_split_u64;
        use super::ddf_bigint;
        use super::ddf_u64;
        use super::derivative_bigint;
        use super::derivative_u64;
        use super::edf_bigint;
        use super::edf_u64;
        use super::factor_bigint;
        use super::factor_u64;
        use super::is_irreducible_bigint;
        use super::is_irreducible_u64;
        use super::monic;
        use super::powmod_bigint;
        use super::powmod_u64;
        use super::squarefree_bigint;
        use super::squarefree_u64;

        let mut rng = Lcg(0xCAFE_BABE_1234_5678);
        let huge_exp: BigInt = (BigInt::from(1) << 130) + BigInt::from(12_345);

        for p_u64 in [2u64, 3, 5, 17, 101, 4093] {
            let p = BigInt::from(p_u64);
            let barrett =
                fast_barrett(&p).unwrap_or_else(|| panic!("prime {p_u64} must fit fast_barrett"));

            let mut cases: Vec<Poly> = Vec::new();
            // Repeated factors and p-th powers: lc * (x + 1)^p * (x^2 + 1)^2 * (x + a)
            let lin1 = norm(vec![BigInt::from(1), BigInt::from(1)], &p);
            let lin_a = norm(vec![BigInt::from(1), BigInt::from(p_u64 - 1)], &p);
            let quad = norm(vec![BigInt::from(1), BigInt::from(1), BigInt::from(1)], &p);
            let mut p_power = vec![BigInt::from(1)];
            for _ in 0..p_u64.min(7) {
                p_power = pmul(&p_power, &lin1, &p);
            }
            let lc_poly = vec![BigInt::from((p_u64 - 1).max(1))];
            let rep = pmul(
                &lc_poly,
                &pmul(&p_power, &pmul(&pmul(&quad, &quad, &p), &lin_a, &p), &p),
                &p,
            );
            cases.push(rep);
            cases.push(vec![BigInt::from(1)]);
            cases.push(vec![BigInt::from(1), BigInt::from(0)]);

            for _ in 0..8 {
                let deg = 2 + usize::try_from(rng.next(11)).unwrap_or(0);
                let mut coeffs: Poly = (0..=deg).map(|_| BigInt::from(rng.next(p_u64))).collect();
                if coeffs[0].is_zero() {
                    coeffs[0] = BigInt::from((1 + rng.next(p_u64 - 1)).max(1));
                }
                cases.push(norm(coeffs, &p));
            }

            for f in cases {
                let fu = poly_to_u64(&f, &p, barrett);
                assert_eq!(poly_from_u64(&fu), f);

                // derivative & irreducibility parity
                assert_eq!(
                    poly_from_u64(&derivative_u64(&fu, barrett)),
                    derivative_bigint(&f, &p)
                );
                assert_eq!(
                    is_irreducible_u64(&fu, barrett),
                    is_irreducible_bigint(&f, &p)
                );

                // powmod parity (e = 0, 1, p, and > 2^128)
                if f.len() >= 2 {
                    let base = vec![BigInt::from(1), BigInt::from(1), BigInt::from(0)];
                    let base_u = poly_to_u64(&base, &p, barrett);
                    for exp in [
                        BigInt::from(0),
                        BigInt::from(1),
                        p.clone(),
                        huge_exp.clone(),
                    ] {
                        assert_eq!(
                            powmod_u64(&base_u, &exp, &fu, barrett).map(|r| poly_from_u64(&r)),
                            powmod_bigint(&base, &exp, &f, &p),
                            "powmod mismatch for p={p_u64}, exp={exp}"
                        );
                    }
                }

                // squarefree parity
                let mf = monic(&f, &p).unwrap_or_else(|| panic!("monic failed"));
                let mfu = monic_u64(&fu, barrett).unwrap_or_else(|| panic!("monic_u64 failed"));
                let sq_big = squarefree_bigint(&mf, &p)
                    .unwrap_or_else(|| panic!("squarefree_bigint failed"));
                let sq_u64: Vec<(Poly, usize)> = squarefree_u64(&mfu, barrett)
                    .unwrap_or_else(|| panic!("squarefree_u64 failed"))
                    .into_iter()
                    .map(|(g, m)| (poly_from_u64(&g), m))
                    .collect();
                assert_eq!(sq_u64, sq_big, "squarefree mismatch mod {p_u64}");

                // ddf, edf, and berlekamp parity on squarefree components
                for (sq, _) in &sq_big {
                    let squ = poly_to_u64(sq, &p, barrett);
                    let ddf_big = ddf_bigint(sq, &p).unwrap_or_else(|| panic!("ddf_bigint failed"));
                    let ddf_fast: Vec<(Poly, usize)> = ddf_u64(&squ, barrett)
                        .unwrap_or_else(|| panic!("ddf_u64 failed"))
                        .into_iter()
                        .map(|(g, d)| (poly_from_u64(&g), d))
                        .collect();
                    assert_eq!(ddf_fast, ddf_big, "ddf mismatch mod {p_u64}");

                    for (gd, d) in ddf_big {
                        let gdu = poly_to_u64(&gd, &p, barrett);
                        let mut parts_big = Vec::new();
                        let mut parts_u64 = Vec::new();
                        let mut rng_big = Rng(0x5eed_1234_abcd_9876);
                        let mut rng_u64 = Rng(0x5eed_1234_abcd_9876);
                        edf_bigint(&gd, d, &p, &mut rng_big, &mut parts_big)
                            .unwrap_or_else(|| panic!("edf_bigint failed"));
                        edf_u64(&gdu, d, barrett, &mut rng_u64, &mut parts_u64)
                            .unwrap_or_else(|| panic!("edf_u64 failed"));
                        let parts_u64_big: Vec<Poly> =
                            parts_u64.iter().map(|q| poly_from_u64(q)).collect();
                        assert_eq!(
                            parts_u64_big, parts_big,
                            "edf RNG/order mismatch mod {p_u64}"
                        );

                        let berl_big = berlekamp_split_bigint(&gd, &p)
                            .unwrap_or_else(|| panic!("berlekamp_split_bigint failed"));
                        let berl_u64: Vec<Poly> = berlekamp_split_u64(&gdu, barrett)
                            .unwrap_or_else(|| panic!("berlekamp_split_u64 failed"))
                            .into_iter()
                            .map(|q| poly_from_u64(&q))
                            .collect();
                        assert_eq!(berl_u64, berl_big, "berlekamp mismatch mod {p_u64}");
                    }
                }

                // full gfp_factor & gfp_berlekamp parity
                let cz_big = factor_bigint(&f, &p, SplitMethod::CantorZassenhaus)
                    .unwrap_or_else(|| panic!("factor_bigint CZ failed"));
                let cz_u64_raw = factor_u64(&fu, barrett, SplitMethod::CantorZassenhaus)
                    .unwrap_or_else(|| panic!("factor_u64 CZ failed"));
                let cz_u64 = (
                    BigInt::from(cz_u64_raw.0),
                    cz_u64_raw
                        .1
                        .into_iter()
                        .map(|(q, m)| (poly_from_u64(&q), m))
                        .collect::<Vec<_>>(),
                );
                assert_eq!(cz_u64, cz_big, "gfp_factor parity mismatch mod {p_u64}");

                let bk_big = factor_bigint(&f, &p, SplitMethod::Berlekamp)
                    .unwrap_or_else(|| panic!("factor_bigint Berlekamp failed"));
                let bk_u64_raw = factor_u64(&fu, barrett, SplitMethod::Berlekamp)
                    .unwrap_or_else(|| panic!("factor_u64 Berlekamp failed"));
                let bk_u64 = (
                    BigInt::from(bk_u64_raw.0),
                    bk_u64_raw
                        .1
                        .into_iter()
                        .map(|(q, m)| (poly_from_u64(&q), m))
                        .collect::<Vec<_>>(),
                );
                assert_eq!(bk_u64, bk_big, "gfp_berlekamp parity mismatch mod {p_u64}");
                assert_eq!(cz_u64, bk_u64, "CZ vs Berlekamp mismatch mod {p_u64}");
            }
        }
    }

    #[test]
    fn large_61bit_and_ntt_primes_differential_parity() {
        use num_bigint::BigInt;

        use super::super::finite_field::Poly;
        use super::super::finite_field::fast_barrett;
        use super::super::finite_field::pmul;
        use super::super::finite_field::poly_from_u64;
        use super::super::finite_field::poly_to_u64;
        use super::Rng;
        use super::SplitMethod;
        use super::edf_bigint;
        use super::edf_u64;
        use super::factor_bigint;
        use super::factor_u64;
        use super::is_irreducible_bigint;
        use super::is_irreducible_u64;

        // 1. 30-bit NTT prime 998_244_353 with degree > 32 (exercises poly_mul_ntt_mod inside pmul_u64).
        let p_ntt_u64 = 998_244_353u64;
        let p_ntt = BigInt::from(p_ntt_u64);
        let b_ntt = fast_barrett(&p_ntt).unwrap_or_else(|| panic!("barrett"));
        let mut f_ntt = vec![BigInt::from(7)];
        let linear_roots = [1u64, 5, 19, 42, 1001, 7777];
        for (idx, &r) in linear_roots.iter().enumerate() {
            let lin = vec![BigInt::from(1), BigInt::from(p_ntt_u64 - r)];
            let mult = if idx < 2 {
                3
            } else if idx < 4 {
                2
            } else {
                1
            };
            for _ in 0..mult {
                f_ntt = pmul(&f_ntt, &lin, &p_ntt);
            }
        }
        // Multiply by irreducible quadratics to push degree >= 34.
        let mut irred_quads = Vec::new();
        for c in 1u64..100 {
            let q = vec![BigInt::from(1), BigInt::from(1), BigInt::from(c)];
            if is_irreducible_u64(&poly_to_u64(&q, &p_ntt, b_ntt), b_ntt) == Some(true) {
                irred_quads.push(q);
                if irred_quads.len() == 11 {
                    break;
                }
            }
        }
        assert_eq!(irred_quads.len(), 11);
        for q in &irred_quads {
            f_ntt = pmul(&f_ntt, q, &p_ntt);
        }
        let f_ntt_u = poly_to_u64(&f_ntt, &p_ntt, b_ntt);
        assert!(f_ntt_u.len() >= 34);

        let res_u64 = factor_u64(&f_ntt_u, b_ntt, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("factor_u64 NTT failed"));
        let res_big = factor_bigint(&f_ntt, &p_ntt, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("factor_bigint NTT failed"));
        let res_u64_big = (
            BigInt::from(res_u64.0),
            res_u64
                .1
                .into_iter()
                .map(|(q, m)| (poly_from_u64(&q), m))
                .collect::<Vec<_>>(),
        );
        assert_eq!(res_u64_big, res_big);
        assert_eq!(res_big.1.len(), linear_roots.len() + irred_quads.len());

        // 2. 61-bit Mersenne prime 2^61 - 1 (< 2^62).
        let p61_u64 = (1u64 << 61) - 1;
        let p61 = BigInt::from(p61_u64);
        let b61 =
            fast_barrett(&p61).unwrap_or_else(|| panic!("61-bit prime must use fast_barrett"));

        // Find two distinct irreducible quadratics and two distinct irreducible cubics mod p61.
        let mut quads61 = Vec::new();
        for c in 1u64..100 {
            let q = vec![BigInt::from(1), BigInt::from(0), BigInt::from(c)];
            let qu = poly_to_u64(&q, &p61, b61);
            let iu = is_irreducible_u64(&qu, b61);
            let ib = is_irreducible_bigint(&q, &p61);
            assert_eq!(iu, ib);
            if iu == Some(true) {
                quads61.push(q);
                if quads61.len() == 2 {
                    break;
                }
            }
        }
        assert_eq!(quads61.len(), 2);

        let mut cubics61 = Vec::new();
        for c in 1u64..200 {
            let cu_poly = vec![
                BigInt::from(1),
                BigInt::from(0),
                BigInt::from(1),
                BigInt::from(c),
            ];
            let cu = poly_to_u64(&cu_poly, &p61, b61);
            let iu = is_irreducible_u64(&cu, b61);
            let ib = is_irreducible_bigint(&cu_poly, &p61);
            assert_eq!(iu, ib);
            if iu == Some(true) {
                cubics61.push(cu_poly);
                if cubics61.len() == 2 {
                    break;
                }
            }
        }
        assert_eq!(cubics61.len(), 2);

        // Test gfp_edf for d = 1 (Word), d = 2 (U128), and d = 3 (Big exponent p^3 > 2^128).
        let lin1 = vec![BigInt::from(1), BigInt::from(p61_u64 - 17)];
        let lin2 = vec![BigInt::from(1), BigInt::from(p61_u64 - 999_999)];
        let lin_prod = pmul(&lin1, &lin2, &p61);
        let quad_prod = pmul(&quads61[0], &quads61[1], &p61);
        let cubic_prod = pmul(&cubics61[0], &cubics61[1], &p61);
        for (d, prod_d) in [(1usize, &lin_prod), (2, &quad_prod), (3, &cubic_prod)] {
            let prod_u = poly_to_u64(prod_d, &p61, b61);
            let mut rng_b = Rng(0x5eed_1234_abcd_9876);
            let mut rng_u = Rng(0x5eed_1234_abcd_9876);
            let mut parts_b = Vec::new();
            let mut parts_u = Vec::new();
            edf_bigint(prod_d, d, &p61, &mut rng_b, &mut parts_b)
                .unwrap_or_else(|| panic!("edf_bigint 61-bit failed for d={d}"));
            edf_u64(&prod_u, d, b61, &mut rng_u, &mut parts_u)
                .unwrap_or_else(|| panic!("edf_u64 61-bit failed for d={d}"));
            let parts_u_big: Vec<Poly> = parts_u.iter().map(|q| poly_from_u64(q)).collect();
            assert_eq!(parts_u_big, parts_b, "61-bit edf parity failed for d={d}");
        }

        // Test full gfp_factor with repeated factors and non-monic lead over GF(2^61 - 1).
        let lc61 = vec![BigInt::from(p61_u64 - 5)];
        let f61 = pmul(
            &lc61,
            &pmul(&pmul(&lin_prod, &lin_prod, &p61), &quad_prod, &p61),
            &p61,
        );
        let f61_u = poly_to_u64(&f61, &p61, b61);
        let fac61_u = factor_u64(&f61_u, b61, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("61-bit factor_u64 failed"));
        let fac61_b = factor_bigint(&f61, &p61, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("61-bit factor_bigint failed"));
        let fac61_u_big = (
            BigInt::from(fac61_u.0),
            fac61_u
                .1
                .into_iter()
                .map(|(q, m)| (poly_from_u64(&q), m))
                .collect::<Vec<_>>(),
        );
        assert_eq!(fac61_u_big, fac61_b);
        assert_eq!(fac61_u_big.0, BigInt::from(p61_u64 - 5));

        // Verify end-to-end graph operators on 61-bit prime.
        let p61_i64 = i64::try_from(p61_u64).unwrap_or_else(|_| panic!("61-bit fits i64"));
        let f61_str = poly(
            &f61.iter()
                .map(|c| i64::try_from(c).unwrap_or(0))
                .collect::<Vec<_>>(),
        );
        let (lead_s, facs_s) = factorisation(&s(&format!("gfp_factor({f61_str}, {p61_i64})")));
        assert_eq!(lead_s, p61_i64 - 5);
        assert_eq!(facs_s.len(), 4);
        assert_eq!(
            s(&format!(
                "gfp_isirreducible({}, {p61_i64})",
                poly(
                    &quads61[0]
                        .iter()
                        .map(|c| i64::try_from(c).unwrap_or(0))
                        .collect::<Vec<_>>()
                )
            )),
            "true"
        );
        assert_eq!(
            s(&format!("gfp_isirreducible({f61_str}, {p61_i64})")),
            "false"
        );
    }

    #[test]
    fn primes_above_62_bits_use_bigint_fallback() {
        use num_bigint::BigInt;

        use super::super::finite_field::fast_barrett;
        use super::super::finite_field::pmul;
        use super::SplitMethod;
        use super::ddf_bigint;
        use super::factor;
        use super::is_irreducible;

        // 2^64 - 59 fits in u64 but >= 2^62, so fast_barrett returns None.
        let p64 = BigInt::from(u64::MAX - 58);
        assert!(fast_barrett(&p64).is_none());
        let lin1 = vec![BigInt::from(1), &p64 - BigInt::from(3)];
        let lin2 = vec![BigInt::from(1), &p64 - BigInt::from(7)];
        let lc64 = vec![BigInt::from(5)];
        let f64 = pmul(&lc64, &pmul(&lin1, &lin1, &p64), &p64);
        let f64 = pmul(&f64, &lin2, &p64);
        let (lead, facs) = factor(&f64, &p64, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("64-bit fallback factor failed"));
        assert_eq!(lead, BigInt::from(5));
        assert_eq!(facs, vec![(lin2, 1), (lin1, 2)]);

        // 2^89 - 1 exceeds u64 entirely.
        let p89 = (BigInt::from(1) << 89) - BigInt::from(1);
        assert!(fast_barrett(&p89).is_none());
        let g1 = vec![BigInt::from(1), &p89 - BigInt::from(11)];
        let g2 = vec![BigInt::from(1), &p89 - BigInt::from(29)];
        let prod89 = pmul(&g1, &g2, &p89);
        assert_eq!(is_irreducible(&g1, &p89), Some(true));
        assert_eq!(is_irreducible(&prod89, &p89), Some(false));
        let ddf89 = ddf_bigint(&prod89, &p89).unwrap_or_else(|| panic!("ddf 89-bit failed"));
        assert_eq!(ddf89, vec![(prod89.clone(), 1)]);
        let (lead89, facs89) = factor(&prod89, &p89, SplitMethod::CantorZassenhaus)
            .unwrap_or_else(|| panic!("89-bit fallback factor failed"));
        assert_eq!(lead89, BigInt::from(1));
        assert_eq!(facs89, vec![(g2, 1), (g1, 1)]);
        assert_eq!(
            s("gfp_powmod(list(1, 0), 2, list(1, 0, 1), 2^89 - 1)"),
            format!("list({})", &p89 - BigInt::from(1))
        );
    }

    #[test]
    fn berlekamp_nullspace_parallel_and_deterministic_edge_cases() {
        use num_bigint::BigInt;

        use super::super::finite_field::Poly;
        use super::super::finite_field::fast_barrett;
        use super::super::finite_field::pmul;
        use super::super::finite_field::poly_from_u64;
        use super::super::finite_field::poly_to_u64;
        use super::berlekamp_split_bigint;
        use super::berlekamp_split_u64;
        use super::is_irreducible;
        use super::nullspace_bigint;
        use super::nullspace_u64;

        // 1. Differential test of nullspace_u64 vs nullspace_bigint on n = 75 (>= 64 parallel threshold)
        // across GF(2), GF(3), and GF(101).
        let n = 75usize;
        let mut rng = Lcg(0xDEAD_BEEF_9988_7766);
        for p_u64 in [2u64, 3, 101] {
            let p = BigInt::from(p_u64);
            let barrett = fast_barrett(&p).unwrap_or_else(|| panic!("barrett"));
            let mut mat_big = vec![vec![BigInt::from(0); n]; n];
            let mut mat_flat = vec![0u64; n * n];
            // Build rank-deficient matrix (last 3 rows are linear combinations of earlier rows).
            for r in 0..n - 3 {
                for c in 0..n {
                    let val = rng.next(p_u64);
                    mat_big[r][c] = BigInt::from(val);
                    mat_flat[r * n + c] = val;
                }
            }
            for c in 0..n {
                let v1 = (mat_flat[c] + mat_flat[n + c]) % p_u64;
                let v2 = (2 * mat_flat[2 * n + c] + mat_flat[3 * n + c]) % p_u64;
                let v3 = mat_flat[4 * n + c];
                for (row_offset, v) in [(n - 3, v1), (n - 2, v2), (n - 1, v3)] {
                    mat_big[row_offset][c] = BigInt::from(v);
                    mat_flat[row_offset * n + c] = v;
                }
            }
            let basis_big =
                nullspace_bigint(mat_big, &p).unwrap_or_else(|| panic!("nullspace_bigint failed"));
            let basis_u64: Vec<Poly> = nullspace_u64(mat_flat, n, barrett)
                .unwrap_or_else(|| panic!("nullspace_u64 failed"))
                .into_iter()
                .map(|v| poly_from_u64(&v))
                .collect();
            assert_eq!(basis_u64, basis_big, "nullspace mismatch for p={p_u64}");
            assert!(basis_u64.len() >= 3);
        }

        // 2. High-degree Berlekamp parity over GF(2) (degree 12 product of all 3 degree-4 irreducibles).
        let p2 = BigInt::from(2);
        let b2 = fast_barrett(&p2).unwrap_or_else(|| panic!("barrett"));
        let mut irreds_deg4 = Vec::new();
        for code in 0..16i64 {
            let mut f = vec![BigInt::from(1)];
            for bit in (0..4).rev() {
                f.push(BigInt::from((code >> bit) & 1));
            }
            if is_irreducible(&f, &p2) == Some(true) {
                irreds_deg4.push(f);
            }
        }
        assert_eq!(irreds_deg4.len(), 3);
        let mut prod2 = vec![BigInt::from(1)];
        for f in &irreds_deg4 {
            prod2 = pmul(&prod2, f, &p2);
        }
        let prod2_u = poly_to_u64(&prod2, &p2, b2);
        let berl2_u: Vec<Poly> = berlekamp_split_u64(&prod2_u, b2)
            .unwrap_or_else(|| panic!("berlekamp_split_u64 GF(2) failed"))
            .into_iter()
            .map(|q| poly_from_u64(&q))
            .collect();
        let berl2_b = berlekamp_split_bigint(&prod2, &p2)
            .unwrap_or_else(|| panic!("berlekamp_split_bigint GF(2) failed"));
        assert_eq!(berl2_u, berl2_b);

        // 3. Deterministic repeated evaluations and edge cases.
        let f_expr = "gfp_factor(list(2, 0, 1, 4, 3, 2, 1), 7)";
        assert_eq!(s(f_expr), s(f_expr));
        let b_expr = "gfp_berlekamp(list(2, 0, 1, 4, 3, 2, 1), 7)";
        assert_eq!(s(b_expr), s(f_expr));
        let edf_expr = "gfp_edf(list(1, 1, 1, 1, 1, 1, 1), 3, 2)";
        assert_eq!(s(edf_expr), s(edf_expr));

        // Empty and constant edge cases do not panic.
        assert_eq!(s("gfp_edf(list(), 1, 5)"), "gfp_edf(list(), 1, 5)");
        assert_eq!(s("gfp_isirreducible(list(), 5)"), "false");
        assert_eq!(s("gfp_ddf(list(1), 5)"), "gfp_ddf(list(1), 5)");
    }
}
