//! Finite fields with exact (arbitrary precision) integers: prime fields
//! GF(p), polynomials over GF(p) and extension fields GF(p)\[x\]/(m).
//!
//! Polynomials are coefficient lists with the **highest degree first** and
//! no leading zeros; the zero polynomial is `list()`. An element of an
//! extension field is a polynomial of degree below that of the modulus
//! polynomial `m`, which the caller supplies (irreducibility is not
//! checked, but a non-invertible element has no inverse and the request
//! stays unreduced).
//!
//! | operator | value |
//! |---|---|
//! | `gf_add(a, b, p)`, `gf_sub`, `gf_mul`, `gf_neg(a, p)` | arithmetic in GF(p) |
//! | `gf_inv(a, p)`, `gf_div(a, b, p)`, `gf_pow(a, e, p)` | inverse (when `gcd(a, p) = 1`), quotient, power (negative `e` allowed) |
//! | `gf_is_zero(a, p)`, `gf_is_one(a, p)` | tests |
//! | `gfp_norm(f, p)`, `gfp_degree(f, p)` | reduced coefficients without leading zeros; degree (`-1` for zero) |
//! | `gfp_add(f, g, p)`, `gfp_sub`, `gfp_mul` | polynomial arithmetic |
//! | `gfp_divmod(f, g, p)` | `list(quotient, remainder)` |
//! | `gfp_gcd(f, g, p)`, `gfp_egcd(f, g, p)` | monic gcd; `list(d, s, t)` with `s f + t g = d` |
//! | `gfp_eval(f, x, p)`, `gfp_isirreducible(f, p)` | value at `x`; irreducibility of `f` for prime `p` (Rabin's test) |
//! | `gfp_powmod`, `gfp_invmod`, `gfp_squarefree`, `gfp_ddf`, `gfp_edf`, `gfp_factor`, `gfp_berlekamp`, `factor_mod` | see [`gf_factor`](super::gf_factor) |
//! | `gfx_reduce(a, m, p)`, `gfx_add(a, b, m, p)`, `gfx_sub`, `gfx_mul`, `gfx_neg(a, m, p)` | extension field arithmetic |
//! | `gfx_inv(a, m, p)`, `gfx_div(a, b, m, p)`, `gfx_pow(a, e, m, p)` | inverse, quotient, power |

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::V;
use super::big;
use super::bigs;
use super::def_doc;
use super::mod_inverse;
use super::modulo;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernels::finite_field::BarrettMod;
use crate::kernels::finite_field::poly_mul_ntt_mod;

pub(super) type Poly = Vec<BigInt>;
pub(super) type PolyU64 = Vec<u64>;

/// Returns `Some(BarrettMod)` when `p` fits in the fast 62-bit `u64` path (`1 < p < 2^62`).
#[inline]
pub(super) fn fast_barrett(p: &BigInt) -> Option<BarrettMod> {
    let m = p.to_u64().filter(|&m| m > 1 && m < (1u64 << 62))?;
    BarrettMod::new(m)
}

/// Reduces a `BigInt` coefficient modulo `barrett.modulus()`.
#[inline]
pub(super) fn bigint_to_u64_mod(
    c: &BigInt,
    p: &BigInt,
    barrett: BarrettMod,
) -> u64 {
    if let Some(v) = c.to_u64() {
        return barrett.reduce_u64(v);
    }
    if let Some(v) = c.to_i64() {
        let mag = barrett.reduce_u64(v.unsigned_abs());
        return if v < 0 {
            barrett.neg(mag)
        } else {
            mag
        };
    }
    modulo(c, p).to_u64().unwrap_or(0)
}

/// Converts a highest-degree-first `Poly` (`Vec<BigInt>`) into a normalized `PolyU64`.
pub(super) fn poly_to_u64(
    f: &Poly,
    p: &BigInt,
    barrett: BarrettMod,
) -> PolyU64 {
    let mut out: Vec<u64> = f.iter().map(|c| bigint_to_u64_mod(c, p, barrett)).collect();
    norm_u64_in_place(&mut out, barrett);
    out
}

/// Converts a normalized `PolyU64` back into `Poly` (`Vec<BigInt>`).
pub(super) fn poly_from_u64(f: &[u64]) -> Poly {
    f.iter().map(|&c| BigInt::from(c)).collect()
}

/// Reduces `f` modulo `barrett.modulus()` and strips leading zeros in place.
pub(super) fn norm_u64_in_place(
    f: &mut PolyU64,
    barrett: BarrettMod,
) {
    for c in f.iter_mut() {
        *c = barrett.reduce_u64(*c);
    }
    let first = f.iter().position(|&c| c != 0).unwrap_or(f.len());
    if first > 0 {
        f.drain(..first);
    }
}

/// Returns a normalized copy of `f` with leading zeros stripped.
pub(super) fn norm_u64(
    mut f: PolyU64,
    barrett: BarrettMod,
) -> PolyU64 {
    norm_u64_in_place(&mut f, barrett);
    f
}

/// Highest-degree-first polynomial addition over `u64` residues.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "out has length n = f.len().max(g.len()) and reverse index i < f.len(), g.len()"
)]
pub(super) fn padd_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> PolyU64 {
    let n = f.len().max(g.len());
    let mut out = vec![0u64; n];
    for (i, &c) in f.iter().rev().enumerate() {
        out[n - 1 - i] = barrett.reduce_u64(c);
    }
    for (i, &c) in g.iter().rev().enumerate() {
        let idx = n - 1 - i;
        out[idx] = barrett.add(out[idx], c);
    }
    norm_u64(out, barrett)
}

/// Highest-degree-first polynomial negation over `u64` residues.
pub(super) fn pneg_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> PolyU64 {
    let out: Vec<u64> = f.iter().map(|&c| barrett.neg(c)).collect();
    norm_u64(out, barrett)
}

/// Highest-degree-first polynomial subtraction over `u64` residues.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "out has length n = f.len().max(g.len()) and reverse index i < f.len(), g.len()"
)]
pub(super) fn psub_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> PolyU64 {
    let n = f.len().max(g.len());
    let mut out = vec![0u64; n];
    for (i, &c) in f.iter().rev().enumerate() {
        out[n - 1 - i] = barrett.reduce_u64(c);
    }
    for (i, &c) in g.iter().rev().enumerate() {
        let idx = n - 1 - i;
        out[idx] = barrett.sub(out[idx], c);
    }
    norm_u64(out, barrett)
}

/// Highest-degree-first polynomial multiplication over `u64` residues.
///
/// Uses cached NTT when the degree is large and the modulus is NTT-friendly,
/// otherwise uses a 4-term fused `u128` accumulation kernel with [`BarrettMod`]
/// and pooled scratch memory from [`crate::kernels::reuse::global_pool`].
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "out has length f.len() + g.len() - 1 and i + j < out.len()"
)]
pub(super) fn pmul_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> PolyU64 {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    if f.len() >= 32
        && g.len() >= 32
        && let Some(prod) = poly_mul_ntt_mod(f, g, barrett.modulus())
    {
        return norm_u64(prod, barrett);
    }
    let out_len = f.len() + g.len() - 1;
    let mut out = vec![0u64; out_len];
    for (i, &a) in f.iter().enumerate() {
        let a_red = barrett.reduce_u64(a);
        if a_red == 0 {
            continue;
        }
        let a128 = u128::from(a_red);
        let dst = &mut out[i..i + g.len()];
        let mut cg = g.chunks_exact(4);
        let mut cd = dst.chunks_exact_mut(4);
        for (xg, yd) in cg.by_ref().zip(cd.by_ref()) {
            yd[0] = barrett.reduce_u128(u128::from(yd[0]) + a128 * u128::from(xg[0]));
            yd[1] = barrett.reduce_u128(u128::from(yd[1]) + a128 * u128::from(xg[1]));
            yd[2] = barrett.reduce_u128(u128::from(yd[2]) + a128 * u128::from(xg[2]));
            yd[3] = barrett.reduce_u128(u128::from(yd[3]) + a128 * u128::from(xg[3]));
        }
        for (&xg, yd) in cg.remainder().iter().zip(cd.into_remainder()) {
            *yd = barrett.reduce_u128(u128::from(*yd) + a128 * u128::from(xg));
        }
    }
    norm_u64(out, barrett)
}

/// Euclidean polynomial division `(quotient, remainder)` on highest-degree-first
/// `u64` slices using 4-way unrolled [`crate::kernels::finite_field::slice_axpy_mod`].
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "k + g.len() <= rem.len() for all k in 0..quot_len"
)]
pub(super) fn pdivmod_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> Option<(PolyU64, PolyU64)> {
    let lead_inv = barrett.inv(*g.first()?)?;
    let mut rem = f.to_vec();
    norm_u64_in_place(&mut rem, barrett);
    if rem.len() < g.len() {
        return Some((Vec::new(), rem));
    }
    let quot_len = rem.len() - g.len() + 1;
    let mut quot = vec![0u64; quot_len];
    let g_tail = &g[1..];
    for k in 0..quot_len {
        let coeff = barrett.mul(rem[k], lead_inv);
        quot[k] = coeff;
        rem[k] = 0;
        if coeff != 0 && !g_tail.is_empty() {
            let neg_coeff = barrett.neg(coeff);
            let window = &mut rem[k + 1..k + g.len()];
            let mut cg = g_tail.chunks_exact(4);
            let mut cw = window.chunks_exact_mut(4);
            let a128 = u128::from(neg_coeff);
            for (xg, yw) in cg.by_ref().zip(cw.by_ref()) {
                yw[0] = barrett.reduce_u128(u128::from(yw[0]) + a128 * u128::from(xg[0]));
                yw[1] = barrett.reduce_u128(u128::from(yw[1]) + a128 * u128::from(xg[1]));
                yw[2] = barrett.reduce_u128(u128::from(yw[2]) + a128 * u128::from(xg[2]));
                yw[3] = barrett.reduce_u128(u128::from(yw[3]) + a128 * u128::from(xg[3]));
            }
            for (&xg, yw) in cg.remainder().iter().zip(cw.into_remainder()) {
                *yw = barrett.reduce_u128(u128::from(*yw) + a128 * u128::from(xg));
            }
        }
    }
    Some((norm_u64(quot, barrett), norm_u64(rem, barrett)))
}

/// Polynomial remainder `f mod m` on highest-degree-first `u64` slices.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "k + m.len() <= rem.len() for all k in 0..quot_len"
)]
pub(super) fn reduce_u64(
    f: &[u64],
    m: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    let lead_inv = barrett.inv(*m.first()?)?;
    let mut rem = f.to_vec();
    norm_u64_in_place(&mut rem, barrett);
    if rem.len() < m.len() {
        return Some(rem);
    }
    let quot_len = rem.len() - m.len() + 1;
    let m_tail = &m[1..];
    for k in 0..quot_len {
        let coeff = barrett.mul(rem[k], lead_inv);
        rem[k] = 0;
        if coeff != 0 && !m_tail.is_empty() {
            let neg_coeff = barrett.neg(coeff);
            let window = &mut rem[k + 1..k + m.len()];
            let mut cm = m_tail.chunks_exact(4);
            let mut cw = window.chunks_exact_mut(4);
            let a128 = u128::from(neg_coeff);
            for (xm, yw) in cm.by_ref().zip(cw.by_ref()) {
                yw[0] = barrett.reduce_u128(u128::from(yw[0]) + a128 * u128::from(xm[0]));
                yw[1] = barrett.reduce_u128(u128::from(yw[1]) + a128 * u128::from(xm[1]));
                yw[2] = barrett.reduce_u128(u128::from(yw[2]) + a128 * u128::from(xm[2]));
                yw[3] = barrett.reduce_u128(u128::from(yw[3]) + a128 * u128::from(xm[3]));
            }
            for (&xm, yw) in cm.remainder().iter().zip(cw.into_remainder()) {
                *yw = barrett.reduce_u128(u128::from(*yw) + a128 * u128::from(xm));
            }
        }
    }
    Some(norm_u64(rem, barrett))
}

/// Extended Euclidean algorithm `(d, s, t)` on `u64` polynomials with `s*f + t*g = d`.
pub(super) fn pegcd_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> Option<(PolyU64, PolyU64, PolyU64)> {
    let mut r0 = norm_u64(f.to_vec(), barrett);
    let mut r1 = norm_u64(g.to_vec(), barrett);
    let (mut s0, mut s1) = (vec![1u64], Vec::new());
    let (mut t0, mut t1) = (Vec::new(), vec![1u64]);
    while !r1.is_empty() {
        let (q, r2) = pdivmod_u64(&r0, &r1, barrett)?;
        let qs1 = pmul_u64(&q, &s1, barrett);
        let qt1 = pmul_u64(&q, &t1, barrett);
        let s2 = psub_u64(&s0, &qs1, barrett);
        let t2 = psub_u64(&t0, &qt1, barrett);
        (r0, r1) = (r1, r2);
        (s0, s1) = (s1, s2);
        (t0, t1) = (t1, t2);
    }
    Some((r0, s0, t0))
}

/// Makes `f` monic over `u64` residues.
pub(super) fn monic_u64(
    f: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    let Some(&lead) = f.first() else {
        return Some(Vec::new());
    };
    let inv = barrett.inv(lead)?;
    if inv == 1 {
        return Some(norm_u64(f.to_vec(), barrett));
    }
    let mut out = vec![0u64; f.len()];
    crate::kernels::finite_field::slice_scale_mod(f, inv, barrett.modulus(), &mut out);
    Some(norm_u64(out, barrett))
}

/// Monic GCD of `f` and `g` on `u64` slices (without computing Bézout coefficients).
pub(super) fn pgcd_u64(
    f: &[u64],
    g: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    let mut r0 = norm_u64(f.to_vec(), barrett);
    let mut r1 = norm_u64(g.to_vec(), barrett);
    while !r1.is_empty() {
        let r2 = reduce_u64(&r0, &r1, barrett)?;
        r0 = r1;
        r1 = r2;
    }
    monic_u64(&r0, barrett)
}

/// Inverse of `x` modulo `m` over `GF(p)` on `u64` slices.
#[expect(
    clippy::indexing_slicing,
    reason = "d.len() == 1 guarantees d[0] is valid"
)]
pub(super) fn xinv_u64(
    x: &[u64],
    m: &[u64],
    barrett: BarrettMod,
) -> Option<PolyU64> {
    let x_red = reduce_u64(x, m, barrett)?;
    let (d, s, _) = pegcd_u64(&x_red, m, barrett)?;
    if d.len() != 1 {
        return None;
    }
    let d_inv = barrett.inv(d[0])?;
    let mut scaled = vec![0u64; s.len()];
    crate::kernels::finite_field::slice_scale_mod(&s, d_inv, barrett.modulus(), &mut scaled);
    reduce_u64(&norm_u64(scaled, barrett), m, barrett)
}

/// Reduces the coefficients modulo `p` and strips leading zeros.
pub(super) fn norm(
    mut f: Poly,
    p: &BigInt,
) -> Poly {
    if let Some(barrett) = fast_barrett(p) {
        return poly_from_u64(&poly_to_u64(&f, p, barrett));
    }
    for c in &mut f {
        *c = modulo(c, p);
    }
    let first = f.iter().position(|c| !c.is_zero()).unwrap_or(f.len());
    f.drain(..first);
    f
}

pub(super) fn padd(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Poly {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let gu = poly_to_u64(g, p, barrett);
        return poly_from_u64(&padd_u64(&fu, &gu, barrett));
    }
    let n = f.len().max(g.len());
    let mut out = vec![BigInt::zero(); n];
    for (i, c) in f.iter().rev().enumerate() {
        out[n - 1 - i] += c;
    }
    for (i, c) in g.iter().rev().enumerate() {
        out[n - 1 - i] += c;
    }
    norm(out, p)
}

fn pneg(
    f: &Poly,
    p: &BigInt,
) -> Poly {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        return poly_from_u64(&pneg_u64(&fu, barrett));
    }
    norm(f.iter().map(|c| -c).collect(), p)
}

pub(super) fn psub(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Poly {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let gu = poly_to_u64(g, p, barrett);
        return poly_from_u64(&psub_u64(&fu, &gu, barrett));
    }
    padd(f, &pneg(g, p), p)
}

pub(super) fn pmul(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Poly {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let gu = poly_to_u64(g, p, barrett);
        return poly_from_u64(&pmul_u64(&fu, &gu, barrett));
    }
    let mut out = vec![BigInt::zero(); f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        for (j, b) in g.iter().enumerate() {
            out[i + j] += a * b;
        }
    }
    norm(out, p)
}

/// `(quotient, remainder)`; `None` for a zero divisor or one whose leading
/// coefficient is not invertible.
pub(super) fn pdivmod(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Option<(Poly, Poly)> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let gu = poly_to_u64(g, p, barrett);
        let (q, r) = pdivmod_u64(&fu, &gu, barrett)?;
        return Some((poly_from_u64(&q), poly_from_u64(&r)));
    }
    let lead_inv = mod_inverse(g.first()?, p)?;
    let mut rem = f.clone();
    if rem.len() < g.len() {
        return Some((Vec::new(), rem));
    }
    let mut quot = vec![BigInt::zero(); rem.len() - g.len() + 1];
    for k in 0..quot.len() {
        let coeff = modulo(&(&rem[k] * &lead_inv), p);
        for (j, gc) in g.iter().enumerate() {
            rem[k + j] = modulo(&(&rem[k + j] - &coeff * gc), p);
        }
        quot[k] = coeff;
    }
    Some((norm(quot, p), norm(rem, p)))
}

/// `(d, s, t)` with `s f + t g = d`, `d` a gcd (not normalised).
pub(super) fn pegcd(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Option<(Poly, Poly, Poly)> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let gu = poly_to_u64(g, p, barrett);
        let (d, s, t) = pegcd_u64(&fu, &gu, barrett)?;
        return Some((poly_from_u64(&d), poly_from_u64(&s), poly_from_u64(&t)));
    }
    let (mut r0, mut r1) = (f.clone(), g.clone());
    let (mut s0, mut s1) = (vec![BigInt::one()], Vec::new());
    let (mut t0, mut t1) = (Vec::new(), vec![BigInt::one()]);
    while !r1.is_empty() {
        let (q, r2) = pdivmod(&r0, &r1, p)?;
        let s2 = psub(&s0, &pmul(&q, &s1, p), p);
        let t2 = psub(&t0, &pmul(&q, &t1, p), p);
        (r0, r1) = (r1, r2);
        (s0, s1) = (s1, s2);
        (t0, t1) = (t1, t2);
    }
    Some((r0, s0, t0))
}

pub(super) fn monic(
    f: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        return Some(poly_from_u64(&monic_u64(&fu, barrett)?));
    }
    let Some(lead) = f.first() else {
        return Some(Vec::new());
    };
    let inv = mod_inverse(lead, p)?;
    Some(norm(f.iter().map(|c| c * &inv).collect(), p))
}

fn peval(
    f: &Poly,
    x: &BigInt,
    p: &BigInt,
) -> BigInt {
    if let Some(barrett) = fast_barrett(p) {
        let xu = bigint_to_u64_mod(x, p, barrett);
        let acc = f.iter().fold(0u64, |acc, c| {
            let cu = bigint_to_u64_mod(c, p, barrett);
            barrett.add(barrett.mul(acc, xu), cu)
        });
        return BigInt::from(acc);
    }
    f.iter()
        .fold(BigInt::zero(), |acc, c| modulo(&(acc * x + c), p))
}

fn modulus(
    cx: &Cx<'_>,
    n: NodeId,
) -> Option<BigInt> {
    big(cx.graph, n).filter(|p| *p > BigInt::one())
}

pub(super) fn poly_value(f: Poly) -> V {
    V::ints(f)
}

pub(super) fn read_poly(
    cx: &Cx<'_>,
    n: NodeId,
    p: &BigInt,
) -> Option<Poly> {
    Some(norm(bigs(cx.graph, n)?, p))
}

// ---------------- prime field ----------------

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn scalar2(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    op: fn(&BigInt, &BigInt, &BigInt) -> Option<BigInt>,
) -> Option<V> {
    let [x, y, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (x, y) = (big(cx.graph, *x)?, big(cx.graph, *y)?);
    op(&modulo(&x, &p), &modulo(&y, &p), &p).map(V::Int)
}

fn gf_inv(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, p] = a else { return None };
    let p = modulus(cx, *p)?;
    mod_inverse(&big(cx.graph, *x)?, &p).map(V::Int)
}

fn gf_neg(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, p] = a else { return None };
    let p = modulus(cx, *p)?;
    Some(V::Int(modulo(&-big(cx.graph, *x)?, &p)))
}

fn gf_pow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, e, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (x, e) = (modulo(&big(cx.graph, *x)?, &p), big(cx.graph, *e)?);
    let base = if e.is_negative() {
        mod_inverse(&x, &p)?
    } else {
        x
    };
    Some(V::Int(base.modpow(&e.abs(), &p)))
}

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn gf_is(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    one: bool,
) -> Option<V> {
    let [x, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let x = modulo(&big(cx.graph, *x)?, &p);
    Some(V::Bool(if one {
        x.is_one()
    } else {
        x.is_zero()
    }))
}

// ---------------- polynomials ----------------

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn poly2(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    op: fn(&Poly, &Poly, &BigInt) -> Poly,
) -> Option<V> {
    let [f, g, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (f, g) = (read_poly(cx, *f, &p)?, read_poly(cx, *g, &p)?);
    Some(poly_value(op(&f, &g, &p)))
}

fn gfp_norm(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let p = modulus(cx, *p)?;
    Some(poly_value(read_poly(cx, *f, &p)?))
}

fn gfp_degree(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let f = read_poly(cx, *f, &p)?;
    Some(V::Int(BigInt::from(f.len()) - 1))
}

fn gfp_divmod(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, g, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (f, g) = (read_poly(cx, *f, &p)?, read_poly(cx, *g, &p)?);
    let (q, r) = pdivmod(&f, &g, &p)?;
    Some(V::List(vec![poly_value(q), poly_value(r)]))
}

fn gfp_gcd(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, g, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (f, g) = (read_poly(cx, *f, &p)?, read_poly(cx, *g, &p)?);
    if let Some(barrett) = fast_barrett(&p) {
        let fu = poly_to_u64(&f, &p, barrett);
        let gu = poly_to_u64(&g, &p, barrett);
        return Some(poly_value(poly_from_u64(&pgcd_u64(&fu, &gu, barrett)?)));
    }
    let (d, _, _) = pegcd(&f, &g, &p)?;
    Some(poly_value(monic(&d, &p)?))
}

fn gfp_egcd(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, g, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let (f, g) = (read_poly(cx, *f, &p)?, read_poly(cx, *g, &p)?);
    let (d, s, t) = pegcd(&f, &g, &p)?;
    Some(V::List(vec![poly_value(d), poly_value(s), poly_value(t)]))
}

fn gfp_eval(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, x, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let f = read_poly(cx, *f, &p)?;
    Some(V::Int(peval(&f, &big(cx.graph, *x)?, &p)))
}

fn gfp_isirreducible(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let p = modulus(cx, *p)?;
    let f = read_poly(cx, *f, &p)?;
    if !crate::rules::number_theory::is_prime(&p) {
        return None;
    }
    super::gf_factor::is_irreducible(&f, &p).map(V::Bool)
}

// ---------------- extension fields ----------------

fn field(
    cx: &Cx<'_>,
    m: NodeId,
    p: NodeId,
) -> Option<(Poly, BigInt)> {
    let p = modulus(cx, p)?;
    let m = read_poly(cx, m, &p)?;
    (m.len() >= 2).then_some((m, p))
}

pub(super) fn reduce(
    f: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let fu = poly_to_u64(f, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return Some(poly_from_u64(&reduce_u64(&fu, &mu, barrett)?));
    }
    Some(pdivmod(f, m, p)?.1)
}

fn gfx_reduce(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, m, p] = a else { return None };
    let (m, p) = field(cx, *m, *p)?;
    let f = read_poly(cx, *f, &p)?;
    Some(poly_value(reduce(&f, &m, &p)?))
}

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn gfx2(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    op: fn(&Poly, &Poly, &Poly, &BigInt) -> Option<Poly>,
) -> Option<V> {
    let [x, y, m, p] = a else { return None };
    let (m, p) = field(cx, *m, *p)?;
    let (x, y) = (read_poly(cx, *x, &p)?, read_poly(cx, *y, &p)?);
    Some(poly_value(op(&x, &y, &m, &p)?))
}

fn xadd(
    x: &Poly,
    y: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let xu = poly_to_u64(x, p, barrett);
        let yu = poly_to_u64(y, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return Some(poly_from_u64(&reduce_u64(
            &padd_u64(&xu, &yu, barrett),
            &mu,
            barrett,
        )?));
    }
    reduce(&padd(x, y, p), m, p)
}

fn xsub(
    x: &Poly,
    y: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let xu = poly_to_u64(x, p, barrett);
        let yu = poly_to_u64(y, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return Some(poly_from_u64(&reduce_u64(
            &psub_u64(&xu, &yu, barrett),
            &mu,
            barrett,
        )?));
    }
    reduce(&psub(x, y, p), m, p)
}

fn xmul(
    x: &Poly,
    y: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let xu = poly_to_u64(x, p, barrett);
        let yu = poly_to_u64(y, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return Some(poly_from_u64(&reduce_u64(
            &pmul_u64(&xu, &yu, barrett),
            &mu,
            barrett,
        )?));
    }
    reduce(&pmul(x, y, p), m, p)
}

/// The inverse of `x` modulo `m` over GF(p): `s x + t m = d` with `d` a
/// non-zero constant.
pub(super) fn xinv(
    x: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let xu = poly_to_u64(x, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        return Some(poly_from_u64(&xinv_u64(&xu, &mu, barrett)?));
    }
    let x = reduce(x, m, p)?;
    let (d, s, _) = pegcd(&x, m, p)?;
    if d.len() != 1 {
        return None;
    }
    let d_inv = mod_inverse(&d[0], p)?;
    reduce(&norm(s.iter().map(|c| c * &d_inv).collect(), p), m, p)
}

fn xdiv(
    x: &Poly,
    y: &Poly,
    m: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    if let Some(barrett) = fast_barrett(p) {
        let xu = poly_to_u64(x, p, barrett);
        let yu = poly_to_u64(y, p, barrett);
        let mu = poly_to_u64(m, p, barrett);
        let inv_y = xinv_u64(&yu, &mu, barrett)?;
        return Some(poly_from_u64(&reduce_u64(
            &pmul_u64(&xu, &inv_y, barrett),
            &mu,
            barrett,
        )?));
    }
    xmul(x, &xinv(y, m, p)?, m, p)
}

fn gfx_neg(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, m, p] = a else { return None };
    let (m, p) = field(cx, *m, *p)?;
    let x = read_poly(cx, *x, &p)?;
    Some(poly_value(reduce(&pneg(&x, &p), &m, &p)?))
}

fn gfx_inv(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, m, p] = a else { return None };
    let (m, p) = field(cx, *m, *p)?;
    let x = read_poly(cx, *x, &p)?;
    Some(poly_value(xinv(&x, &m, &p)?))
}

fn gfx_pow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, e, m, p] = a else { return None };
    let (m, p) = field(cx, *m, *p)?;
    let x = read_poly(cx, *x, &p)?;
    let e = big(cx.graph, *e)?;
    if let Some(barrett) = fast_barrett(&p) {
        let xu = poly_to_u64(&x, &p, barrett);
        let mu = poly_to_u64(&m, &p, barrett);
        let mut base = if e.is_negative() {
            xinv_u64(&xu, &mu, barrett)?
        } else {
            reduce_u64(&xu, &mu, barrett)?
        };
        let mut exp = e.abs();
        let mut acc = reduce_u64(&[1u64], &mu, barrett)?;
        while !exp.is_zero() {
            if exp.bit(0) {
                acc = reduce_u64(&pmul_u64(&acc, &base, barrett), &mu, barrett)?;
            }
            base = reduce_u64(&pmul_u64(&base, &base, barrett), &mu, barrett)?;
            exp >>= 1;
        }
        return Some(poly_value(poly_from_u64(&acc)));
    }
    let mut base = if e.is_negative() {
        xinv(&x, &m, &p)?
    } else {
        reduce(&x, &m, &p)?
    };
    let mut e = e.abs();
    let mut acc = reduce(&vec![BigInt::one()], &m, &p)?;
    while !e.is_zero() {
        if e.bit(0) {
            acc = xmul(&acc, &base, &m, &p)?;
        }
        base = xmul(&base, &base, &m, &p)?;
        e >>= 1;
    }
    Some(poly_value(acc))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def_doc(
        i,
        "gf_add",
        Arity::Fixed(3),
        |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x + y), p))),
        "Computes `(a + b) mod p` in the prime field GF(p).",
    )?;
    def_doc(
        i,
        "gf_sub",
        Arity::Fixed(3),
        |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x - y), p))),
        "Computes `(a - b) mod p` in the prime field GF(p).",
    )?;
    def_doc(
        i,
        "gf_mul",
        Arity::Fixed(3),
        |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x * y), p))),
        "Computes `(a * b) mod p` in the prime field GF(p).",
    )?;
    def_doc(
        i,
        "gf_div",
        Arity::Fixed(3),
        |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x * mod_inverse(y, p)?), p))),
        "Computes `a * b^(-1) mod p` in GF(p) when `gcd(b, p) = 1`.",
    )?;
    def_doc(
        i,
        "gf_neg",
        Arity::Fixed(2),
        gf_neg,
        "Computes `(-a) mod p` in the prime field GF(p).",
    )?;
    def_doc(
        i,
        "gf_inv",
        Arity::Fixed(2),
        gf_inv,
        "Computes the modular multiplicative inverse `a^(-1) mod p` when `gcd(a, p) = 1`.",
    )?;
    def_doc(
        i,
        "gf_pow",
        Arity::Fixed(3),
        gf_pow,
        "Computes `a^e mod p` in GF(p) for integer exponent `e` (negative exponents use the modular inverse).",
    )?;
    def_doc(
        i,
        "gf_is_zero",
        Arity::Fixed(2),
        |cx, a| gf_is(cx, a, false),
        "Returns `true` if `a = 0 (mod p)` in GF(p).",
    )?;
    def_doc(
        i,
        "gf_is_one",
        Arity::Fixed(2),
        |cx, a| gf_is(cx, a, true),
        "Returns `true` if `a = 1 (mod p)` in GF(p).",
    )?;

    def_doc(
        i,
        "gfp_norm",
        Arity::Fixed(2),
        gfp_norm,
        "Normalizes a highest-degree-first polynomial `f` over GF(p) by reducing coefficients modulo `p` and stripping leading zeros.",
    )?;
    def_doc(
        i,
        "gfp_degree",
        Arity::Fixed(2),
        gfp_degree,
        "Returns the degree of polynomial `f` over GF(p), or `-1` for the zero polynomial.",
    )?;
    def_doc(
        i,
        "gfp_add",
        Arity::Fixed(3),
        |cx, a| poly2(cx, a, padd),
        "Computes the sum `f + g` of two polynomials over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_sub",
        Arity::Fixed(3),
        |cx, a| poly2(cx, a, psub),
        "Computes the difference `f - g` of two polynomials over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_mul",
        Arity::Fixed(3),
        |cx, a| poly2(cx, a, pmul),
        "Computes the product `f * g` of two polynomials over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_divmod",
        Arity::Fixed(3),
        gfp_divmod,
        "Computes Euclidean polynomial division `f / g` over GF(p), returning `list(quotient, remainder)`.",
    )?;
    def_doc(
        i,
        "gfp_gcd",
        Arity::Fixed(3),
        gfp_gcd,
        "Computes the monic greatest common divisor of polynomials `f` and `g` over GF(p).",
    )?;
    def_doc(
        i,
        "gfp_egcd",
        Arity::Fixed(3),
        gfp_egcd,
        "Computes the extended Euclidean algorithm on polynomials `f` and `g` over GF(p), returning `list(d, s, t)` with `s*f + t*g = d`.",
    )?;
    def_doc(
        i,
        "gfp_eval",
        Arity::Fixed(3),
        gfp_eval,
        "Evaluates polynomial `f` at scalar `x` in GF(p).",
    )?;
    def_doc(
        i,
        "gfp_isirreducible",
        Arity::Fixed(2),
        gfp_isirreducible,
        "Tests whether polynomial `f` is irreducible over GF(p) using Rabin's irreducibility test.",
    )?;

    def_doc(
        i,
        "gfx_reduce",
        Arity::Fixed(3),
        gfx_reduce,
        "Reduces polynomial `a` modulo `m` over GF(p) in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_add",
        Arity::Fixed(4),
        |cx, a| gfx2(cx, a, xadd),
        "Computes `(a + b) mod m` in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_sub",
        Arity::Fixed(4),
        |cx, a| gfx2(cx, a, xsub),
        "Computes `(a - b) mod m` in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_mul",
        Arity::Fixed(4),
        |cx, a| gfx2(cx, a, xmul),
        "Computes `(a * b) mod m` in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_div",
        Arity::Fixed(4),
        |cx, a| gfx2(cx, a, xdiv),
        "Computes `a * b^(-1) mod m` in the extension field GF(p)[x]/(m) when `b` is invertible modulo `m`.",
    )?;
    def_doc(
        i,
        "gfx_neg",
        Arity::Fixed(3),
        gfx_neg,
        "Computes `(-a) mod m` in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_inv",
        Arity::Fixed(3),
        gfx_inv,
        "Computes the multiplicative inverse `a^(-1) mod m` in the extension field GF(p)[x]/(m).",
    )?;
    def_doc(
        i,
        "gfx_pow",
        Arity::Fixed(4),
        gfx_pow,
        "Computes `a^e mod m` in the extension field GF(p)[x]/(m) for integer exponent `e`.",
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::nested;
    use super::super::test_util::nums;
    use super::super::test_util::s;

    fn poly(c: &[i64]) -> String {
        format!(
            "list({})",
            c.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    #[test]
    fn prime_field_arithmetic() {
        assert_eq!(s("gf_add(5, 4, 7)"), "2");
        assert_eq!(s("gf_sub(2, 5, 7)"), "4");
        assert_eq!(s("gf_mul(-3, 5, 7)"), "6");
        assert_eq!(s("gf_neg(3, 7)"), "4");
        assert_eq!(s("gf_neg(0, 7)"), "0");
        assert_eq!(s("gf_inv(3, 7)"), "5");
        assert_eq!(s("gf_div(2, 3, 7)"), "3");
        assert_eq!(s("gf_pow(3, 6, 7)"), "1");
        assert_eq!(s("gf_pow(3, -1, 7)"), "5");
        assert_eq!(s("gf_is_zero(14, 7)"), "true");
        assert_eq!(s("gf_is_one(8, 7)"), "true");
        assert_eq!(s("gf_is_one(9, 7)"), "false");
        // zero has no inverse, a composite modulus has zero divisors
        assert_eq!(s("gf_inv(0, 7)"), "gf_inv(0, 7)");
        assert_eq!(s("gf_div(1, 4, 8)"), "gf_div(1, 4, 8)");
        // Fermat and inverses by brute force
        for p in [2_i64, 3, 5, 7, 11, 13] {
            for a in 1..p {
                assert_eq!(s(&format!("gf_pow({a}, {}, {p})", p - 1)), "1");
                let inv: i64 = s(&format!("gf_inv({a}, {p})")).parse().unwrap_or(-1);
                assert_eq!(a * inv % p, 1);
                assert_eq!(s(&format!("gf_mul({a}, {inv}, {p})")), "1");
            }
        }
        // big modulus
        assert_eq!(
            s("gf_mul(2^100, 3, 2^127 - 1)"),
            format!("{}", (3_u128 << 100))
        );
    }

    #[test]
    fn polynomials_over_a_prime_field() {
        assert_eq!(s("gfp_norm(list(0, 0, 8, -1, 3), 5)"), "list(3, 4, 3)");
        assert_eq!(s("gfp_norm(list(5, 10), 5)"), "list()");
        assert_eq!(s("gfp_degree(list(1, 0, 2), 5)"), "2");
        assert_eq!(s("gfp_degree(list(), 5)"), "-1");
        assert_eq!(s("gfp_add(list(1, 2, 3), list(4, 4), 5)"), "list(1, 1, 2)");
        assert_eq!(s("gfp_add(list(1, 2), list(4, 3), 5)"), "list()");
        assert_eq!(s("gfp_sub(list(1, 2, 3), list(1, 2, 3), 5)"), "list()");
        // (x + 1)^2 = x^2 + 2x + 1, (x + 1)(x - 1) = x^2 - 1
        assert_eq!(s("gfp_mul(list(1, 1), list(1, 1), 7)"), "list(1, 2, 1)");
        assert_eq!(s("gfp_mul(list(1, 1), list(1, -1), 7)"), "list(1, 0, 6)");
        assert_eq!(s("gfp_mul(list(), list(1, 1), 7)"), "list()");
        assert_eq!(s("gfp_eval(list(1, 0, 6), 3, 7)"), "1");
    }

    #[test]
    fn division_satisfies_its_definition() {
        let p = 7_i64;
        let dividends = [
            vec![3, 0, 5, 1, 2],
            vec![1, 2, 3],
            vec![5],
            vec![1, 0, 0, 0, 0, 1],
            vec![2, 3, 4, 5, 6, 1],
        ];
        let divisors = [vec![1, 1], vec![3, 0, 2], vec![5], vec![2, 1, 4, 1]];
        for f in &dividends {
            for g in &divisors {
                let dm = nested(&s(&format!("gfp_divmod({}, {}, {p})", poly(f), poly(g))));
                assert_eq!(dm.len(), 2, "{f:?} / {g:?}");
                let (q, r) = (&dm[0], &dm[1]);
                // f = q g + r and deg r < deg g
                let qg = s(&format!("gfp_mul({}, {}, {p})", poly(q), poly(g)));
                let back = s(&format!("gfp_add({qg}, {}, {p})", poly(r)));
                assert_eq!(
                    back,
                    s(&format!("gfp_norm({}, {p})", poly(f))),
                    "{f:?} / {g:?}"
                );
                assert!(r.len() < g.len());
            }
        }
        // division by zero or by a non-invertible leading coefficient stays
        assert_eq!(
            s("gfp_divmod(list(1, 1), list(), 7)"),
            "gfp_divmod(list(1, 1), list(), 7)"
        );
        assert_eq!(
            s("gfp_divmod(list(1, 1), list(2, 1), 4)"),
            "gfp_divmod(list(1, 1), list(2, 1), 4)"
        );
        // the quotient is read in the right orientation: x^2 / x = x
        assert_eq!(
            s("gfp_divmod(list(1, 0, 0), list(1, 0), 5)"),
            "list(list(1, 0), list())"
        );
    }

    #[test]
    fn gcd_and_bezout() {
        // (x+1)(x+2) and (x+1)(x+3) over GF(7)
        let f = s("gfp_mul(list(1, 1), list(1, 2), 7)");
        let g = s("gfp_mul(list(1, 1), list(1, 3), 7)");
        assert_eq!(s(&format!("gfp_gcd({f}, {g}, 7)")), "list(1, 1)");
        // 3(x+1) normalised
        let f3 = s("gfp_mul(list(3), list(1, 1), 7)");
        assert_eq!(s(&format!("gfp_gcd({f3}, {g}, 7)")), "list(1, 1)");
        assert_eq!(s("gfp_gcd(list(1, 0, 1), list(1, 1), 7)"), "list(1)");
        for (a, b) in [
            (f.as_str(), g.as_str()),
            ("list(1, 0, 1)", "list(1, 1)"),
            ("list(1, 2, 3, 4)", "list(2, 5, 1)"),
        ] {
            let e = nested(&s(&format!("gfp_egcd({a}, {b}, 7)")));
            assert_eq!(e.len(), 3);
            let sa = s(&format!("gfp_mul({}, {a}, 7)", poly(&e[1])));
            let tb = s(&format!("gfp_mul({}, {b}, 7)", poly(&e[2])));
            assert_eq!(
                s(&format!("gfp_add({sa}, {tb}, 7)")),
                poly(&e[0]),
                "{a} {b}"
            );
            // d divides both
            for x in [a, b] {
                let r = nested(&s(&format!("gfp_divmod({x}, {}, 7)", poly(&e[0]))));
                assert!(r[1].is_empty(), "not empty");
            }
        }
    }

    #[test]
    fn irreducibility_by_trial_division() {
        assert_eq!(s("gfp_isirreducible(list(1, 1, 1), 2)"), "true");
        assert_eq!(s("gfp_isirreducible(list(1, 0, 1), 2)"), "false");
        assert_eq!(s("gfp_isirreducible(list(1, 1, 0, 1), 2)"), "true");
        assert_eq!(s("gfp_isirreducible(list(1, 0, 1), 3)"), "true");
        assert_eq!(s("gfp_isirreducible(list(1, 0, 1), 5)"), "false");
        assert_eq!(s("gfp_isirreducible(list(1, 0, 0, 1, 1), 2)"), "true");
    }

    /// Checks the field axioms of GF(p)[x]/(m) over every element.
    fn check_field(
        p: i64,
        m: &[i64],
    ) {
        let deg = m.len() - 1;
        let order = p.pow(u32::try_from(deg).unwrap_or(1));
        let element = |mut code: i64| {
            let mut c = Vec::new();
            for _ in 0..deg {
                c.push(code % p);
                code /= p;
            }
            c.reverse();
            c
        };
        let mp = poly(m);
        let mul =
            |a: &[i64], b: &[i64]| s(&format!("gfx_mul({}, {}, {mp}, {p})", poly(a), poly(b)));
        let one = s(&format!("gfx_reduce(list(1), {mp}, {p})"));
        for code in 1..order {
            let a = element(code);
            let inv = s(&format!("gfx_inv({}, {mp}, {p})", poly(&a)));
            assert_eq!(mul(&a, &nums(&inv)), one, "{a:?} in GF({p}^{deg})");
            // a^(q-1) = 1
            assert_eq!(
                s(&format!("gfx_pow({}, {}, {mp}, {p})", poly(&a), order - 1)),
                one
            );
            // a^-1 through gfx_pow and division
            assert_eq!(s(&format!("gfx_pow({}, -1, {mp}, {p})", poly(&a))), inv);
            assert_eq!(
                s(&format!("gfx_div({}, {}, {mp}, {p})", poly(&a), poly(&a))),
                one
            );
        }
        // distributivity and commutativity on a sample
        for x in (0..order).step_by(2) {
            for y in (1..order).step_by(3) {
                for z in (0..order).step_by(5) {
                    let (a, b, c) = (element(x), element(y), element(z));
                    let bc = s(&format!("gfx_add({}, {}, {mp}, {p})", poly(&b), poly(&c)));
                    let lhs = mul(&a, &nums(&bc));
                    let ab = mul(&a, &b);
                    let ac = mul(&a, &c);
                    let rhs = s(&format!("gfx_add({ab}, {ac}, {mp}, {p})"));
                    assert_eq!(lhs, rhs);
                    assert_eq!(mul(&a, &b), mul(&b, &a));
                    let sum = s(&format!("gfx_add({}, {}, {mp}, {p})", poly(&a), poly(&b)));
                    assert_eq!(
                        s(&format!("gfx_sub({sum}, {}, {mp}, {p})", poly(&b))),
                        s(&format!("gfx_reduce({}, {mp}, {p})", poly(&a)))
                    );
                    let n = s(&format!("gfx_neg({}, {mp}, {p})", poly(&a)));
                    assert_eq!(
                        s(&format!("gfx_add({}, {n}, {mp}, {p})", poly(&a))),
                        "list()"
                    );
                }
            }
        }
    }

    #[test]
    fn extension_fields_are_fields() {
        check_field(2, &[1, 1, 1]); // GF(4)
        check_field(2, &[1, 1, 0, 1]); // GF(8)
        check_field(3, &[1, 0, 1]); // GF(9)
        check_field(2, &[1, 0, 0, 1, 1]); // GF(16)
        check_field(5, &[1, 0, 2]); // GF(25): x^2 + 2 is irreducible mod 5
    }

    #[test]
    fn extension_field_specifics() {
        // x^2 = x + 1 in GF(4) = GF(2)[x]/(x^2 + x + 1)
        assert_eq!(
            s("gfx_mul(list(1, 0), list(1, 0), list(1, 1, 1), 2)"),
            "list(1, 1)"
        );
        assert_eq!(
            s("gfx_reduce(list(1, 0, 0), list(1, 1, 1), 2)"),
            "list(1, 1)"
        );
        // with a reducible modulus a zero divisor has no inverse
        assert_eq!(
            s("gfx_inv(list(1, 1), list(1, 0, 1), 2)"),
            "gfx_inv(list(1, 1), list(1, 0, 1), 2)"
        );
        // zero has no inverse
        assert_eq!(
            s("gfx_inv(list(), list(1, 1, 1), 2)"),
            "gfx_inv(list(), list(1, 1, 1), 2)"
        );
        // the multiplicative group of GF(8) is cyclic of order 7
        assert_eq!(s("gfx_pow(list(1, 0), 7, list(1, 1, 0, 1), 2)"), "list(1)");
    }

    #[test]
    fn u64_ntt_and_large_prime_paths_agree_with_reference() {
        use num_bigint::BigInt;

        // 1. NTT-friendly prime (998244353) with degree >= 35 exercises poly_mul_ntt_mod inside pmul_u64.
        let p_ntt = BigInt::from(998_244_353u64);
        let b_ntt = super::fast_barrett(&p_ntt)
            .unwrap_or_else(|| panic!("998244353 must use fast_barrett"));
        let f_ntt: super::Poly = (0..36)
            .map(|i| BigInt::from((i * 12_345 + 7) % 998_244_353u64))
            .collect();
        let g_ntt: super::Poly = (0..38)
            .map(|i| BigInt::from((i * 54_321 + 13) % 998_244_353u64))
            .collect();
        let prod_ntt = super::pmul(&f_ntt, &g_ntt, &p_ntt);
        // Schoolbook BigInt reference
        let mut expected = vec![BigInt::from(0); f_ntt.len() + g_ntt.len() - 1];
        for (i, a) in f_ntt.iter().enumerate() {
            for (j, b) in g_ntt.iter().enumerate() {
                expected[i + j] += a * b;
            }
        }
        let expected = super::norm(expected, &p_ntt);
        assert_eq!(prod_ntt, expected);
        let (q_ntt, r_ntt) =
            super::pdivmod(&prod_ntt, &g_ntt, &p_ntt).unwrap_or_else(|| panic!("pdivmod failed"));
        assert_eq!(q_ntt, f_ntt);
        assert!(r_ntt.is_empty());

        // Verify pgcd_u64 matches gfp_gcd and pegcd for shared factor.
        let common = vec![BigInt::from(1), BigInt::from(17), BigInt::from(99)];
        let f_common = super::pmul(&f_ntt, &common, &p_ntt);
        let g_common = super::pmul(&g_ntt, &common, &p_ntt);
        let gcd_u64 = super::poly_from_u64(
            &super::pgcd_u64(
                &super::poly_to_u64(&f_common, &p_ntt, b_ntt),
                &super::poly_to_u64(&g_common, &p_ntt, b_ntt),
                b_ntt,
            )
            .unwrap_or_else(|| panic!("pgcd_u64 failed")),
        );
        let (d_egcd, _, _) =
            super::pegcd(&f_common, &g_common, &p_ntt).unwrap_or_else(|| panic!("pegcd failed"));
        assert_eq!(
            gcd_u64,
            super::monic(&d_egcd, &p_ntt).unwrap_or_else(|| panic!("monic failed"))
        );

        // 2. 61-bit Mersenne prime 2^61 - 1 (< 2^62) exercises Barrett u64 fast path near the boundary.
        let p61_u64 = (1u64 << 61) - 1;
        let p61 = BigInt::from(p61_u64);
        assert!(super::fast_barrett(&p61).is_some());
        let a61 = vec![
            BigInt::from(p61_u64 - 3),
            BigInt::from(p61_u64 - 19),
            BigInt::from(42),
        ];
        let b61 = vec![BigInt::from(5), BigInt::from(p61_u64 - 7)];
        let m61 = vec![BigInt::from(1), BigInt::from(0), BigInt::from(1)];
        let inv_a = super::xinv(&a61, &m61, &p61).unwrap_or_else(|| panic!("xinv failed"));
        let rem_one = super::reduce(&super::pmul(&a61, &inv_a, &p61), &m61, &p61)
            .unwrap_or_else(|| panic!("reduce failed"));
        assert_eq!(rem_one, vec![BigInt::from(1)]);
        let prod61 = super::pmul(&a61, &b61, &p61);
        let (q61, r61) =
            super::pdivmod(&prod61, &b61, &p61).unwrap_or_else(|| panic!("pdivmod61 failed"));
        assert_eq!(q61, a61);
        assert!(r61.is_empty());
        assert_eq!(
            s(&format!(
                "gfp_gcd(list(1, 0, {}), list(1, {}), {p61_u64})",
                p61_u64 - 9,
                p61_u64 - 3
            )),
            format!("list(1, {})", p61_u64 - 3)
        );

        // 3. 127-bit Mersenne prime 2^127 - 1 (> 2^62) exercises BigInt fallback cleanly.
        let p127 = (BigInt::from(1) << 127) - 1;
        assert!(super::fast_barrett(&p127).is_none());
        let a127 = vec![&p127 - BigInt::from(5), BigInt::from(11)];
        let b127 = vec![BigInt::from(3), &p127 - BigInt::from(2)];
        let prod127 = super::pmul(&a127, &b127, &p127);
        let (q127, r127) =
            super::pdivmod(&prod127, &a127, &p127).unwrap_or_else(|| panic!("pdivmod127 failed"));
        assert_eq!(q127, b127);
        assert!(r127.is_empty());
        assert_eq!(
            s("gfx_mul(list(1, 0), list(1, 1), list(1, 0, 1), 2^127 - 1)"),
            format!("list(1, {})", (BigInt::from(1) << 127) - 2)
        );
    }
}
