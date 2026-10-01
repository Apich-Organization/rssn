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
use num_traits::Zero;

use super::big;
use super::bigs;
use super::def;
use super::mod_inverse;
use super::modulo;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::RuleError;

pub(super) type Poly = Vec<BigInt>;

/// Reduces the coefficients modulo `p` and strips leading zeros.
pub(super) fn norm(
    mut f: Poly,
    p: &BigInt,
) -> Poly {
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
    norm(f.iter().map(|c| -c).collect(), p)
}

pub(super) fn psub(
    f: &Poly,
    g: &Poly,
    p: &BigInt,
) -> Poly {
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
    f.iter().fold(BigInt::zero(), |acc, c| modulo(&(acc * x + c), p))
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
    let base = if e.is_negative() { mod_inverse(&x, &p)? } else { x };
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
    Some(V::Bool(if one { x.is_one() } else { x.is_zero() }))
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

fn xadd(x: &Poly, y: &Poly, m: &Poly, p: &BigInt) -> Option<Poly> {
    reduce(&padd(x, y, p), m, p)
}

fn xsub(x: &Poly, y: &Poly, m: &Poly, p: &BigInt) -> Option<Poly> {
    reduce(&psub(x, y, p), m, p)
}

fn xmul(x: &Poly, y: &Poly, m: &Poly, p: &BigInt) -> Option<Poly> {
    reduce(&pmul(x, y, p), m, p)
}

/// The inverse of `x` modulo `m` over GF(p): `s x + t m = d` with `d` a
/// non-zero constant.
pub(super) fn xinv(x: &Poly, m: &Poly, p: &BigInt) -> Option<Poly> {
    let x = reduce(x, m, p)?;
    let (d, s, _) = pegcd(&x, m, p)?;
    if d.len() != 1 {
        return None;
    }
    let d_inv = mod_inverse(&d[0], p)?;
    reduce(&norm(s.iter().map(|c| c * &d_inv).collect(), p), m, p)
}

fn xdiv(x: &Poly, y: &Poly, m: &Poly, p: &BigInt) -> Option<Poly> {
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
    let mut base = if e.is_negative() { xinv(&x, &m, &p)? } else { reduce(&x, &m, &p)? };
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
    def(i, "gf_add", Arity::Fixed(3), |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x + y), p))))?;
    def(i, "gf_sub", Arity::Fixed(3), |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x - y), p))))?;
    def(i, "gf_mul", Arity::Fixed(3), |cx, a| scalar2(cx, a, |x, y, p| Some(modulo(&(x * y), p))))?;
    def(i, "gf_div", Arity::Fixed(3), |cx, a| {
        scalar2(cx, a, |x, y, p| Some(modulo(&(x * mod_inverse(y, p)?), p)))
    })?;
    def(i, "gf_neg", Arity::Fixed(2), gf_neg)?;
    def(i, "gf_inv", Arity::Fixed(2), gf_inv)?;
    def(i, "gf_pow", Arity::Fixed(3), gf_pow)?;
    def(i, "gf_is_zero", Arity::Fixed(2), |cx, a| gf_is(cx, a, false))?;
    def(i, "gf_is_one", Arity::Fixed(2), |cx, a| gf_is(cx, a, true))?;

    def(i, "gfp_norm", Arity::Fixed(2), gfp_norm)?;
    def(i, "gfp_degree", Arity::Fixed(2), gfp_degree)?;
    def(i, "gfp_add", Arity::Fixed(3), |cx, a| poly2(cx, a, padd))?;
    def(i, "gfp_sub", Arity::Fixed(3), |cx, a| poly2(cx, a, psub))?;
    def(i, "gfp_mul", Arity::Fixed(3), |cx, a| poly2(cx, a, pmul))?;
    def(i, "gfp_divmod", Arity::Fixed(3), gfp_divmod)?;
    def(i, "gfp_gcd", Arity::Fixed(3), gfp_gcd)?;
    def(i, "gfp_egcd", Arity::Fixed(3), gfp_egcd)?;
    def(i, "gfp_eval", Arity::Fixed(3), gfp_eval)?;
    def(i, "gfp_isirreducible", Arity::Fixed(2), gfp_isirreducible)?;

    def(i, "gfx_reduce", Arity::Fixed(3), gfx_reduce)?;
    def(i, "gfx_add", Arity::Fixed(4), |cx, a| gfx2(cx, a, xadd))?;
    def(i, "gfx_sub", Arity::Fixed(4), |cx, a| gfx2(cx, a, xsub))?;
    def(i, "gfx_mul", Arity::Fixed(4), |cx, a| gfx2(cx, a, xmul))?;
    def(i, "gfx_div", Arity::Fixed(4), |cx, a| gfx2(cx, a, xdiv))?;
    def(i, "gfx_neg", Arity::Fixed(3), gfx_neg)?;
    def(i, "gfx_inv", Arity::Fixed(3), gfx_inv)?;
    def(i, "gfx_pow", Arity::Fixed(4), gfx_pow)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::nested;
    use super::super::test_util::nums;
    use super::super::test_util::s;

    fn poly(c: &[i64]) -> String {
        format!("list({})", c.iter().map(ToString::to_string).collect::<Vec<_>>().join(", "))
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
        assert_eq!(s("gf_mul(2^100, 3, 2^127 - 1)"), format!("{}", (3_u128 << 100)));
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
        let dividends = [vec![3, 0, 5, 1, 2], vec![1, 2, 3], vec![5], vec![1, 0, 0, 0, 0, 1], vec![2, 3, 4, 5, 6, 1]];
        let divisors = [vec![1, 1], vec![3, 0, 2], vec![5], vec![2, 1, 4, 1]];
        for f in &dividends {
            for g in &divisors {
                let dm = nested(&s(&format!("gfp_divmod({}, {}, {p})", poly(f), poly(g))));
                assert_eq!(dm.len(), 2, "{f:?} / {g:?}");
                let (q, r) = (&dm[0], &dm[1]);
                // f = q g + r and deg r < deg g
                let qg = s(&format!("gfp_mul({}, {}, {p})", poly(q), poly(g)));
                let back = s(&format!("gfp_add({qg}, {}, {p})", poly(r)));
                assert_eq!(back, s(&format!("gfp_norm({}, {p})", poly(f))), "{f:?} / {g:?}");
                assert!(r.len() < g.len());
            }
        }
        // division by zero or by a non-invertible leading coefficient stays
        assert_eq!(s("gfp_divmod(list(1, 1), list(), 7)"), "gfp_divmod(list(1, 1), list(), 7)");
        assert_eq!(s("gfp_divmod(list(1, 1), list(2, 1), 4)"), "gfp_divmod(list(1, 1), list(2, 1), 4)");
        // the quotient is read in the right orientation: x^2 / x = x
        assert_eq!(s("gfp_divmod(list(1, 0, 0), list(1, 0), 5)"), "list(list(1, 0), list())");
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
        for (a, b) in [(f.as_str(), g.as_str()), ("list(1, 0, 1)", "list(1, 1)"), ("list(1, 2, 3, 4)", "list(2, 5, 1)")] {
            let e = nested(&s(&format!("gfp_egcd({a}, {b}, 7)")));
            assert_eq!(e.len(), 3);
            let sa = s(&format!("gfp_mul({}, {a}, 7)", poly(&e[1])));
            let tb = s(&format!("gfp_mul({}, {b}, 7)", poly(&e[2])));
            assert_eq!(s(&format!("gfp_add({sa}, {tb}, 7)")), poly(&e[0]), "{a} {b}");
            // d divides both
            for x in [a, b] {
                let r = nested(&s(&format!("gfp_divmod({x}, {}, 7)", poly(&e[0]))));
                assert!(r[1].is_empty());
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
        let mul = |a: &[i64], b: &[i64]| s(&format!("gfx_mul({}, {}, {mp}, {p})", poly(a), poly(b)));
        let one = s(&format!("gfx_reduce(list(1), {mp}, {p})"));
        for code in 1..order {
            let a = element(code);
            let inv = s(&format!("gfx_inv({}, {mp}, {p})", poly(&a)));
            assert_eq!(mul(&a, &nums(&inv)), one, "{a:?} in GF({p}^{deg})");
            // a^(q-1) = 1
            assert_eq!(s(&format!("gfx_pow({}, {}, {mp}, {p})", poly(&a), order - 1)), one);
            // a^-1 through gfx_pow and division
            assert_eq!(s(&format!("gfx_pow({}, -1, {mp}, {p})", poly(&a))), inv);
            assert_eq!(s(&format!("gfx_div({}, {}, {mp}, {p})", poly(&a), poly(&a))), one);
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
                    assert_eq!(s(&format!("gfx_sub({sum}, {}, {mp}, {p})", poly(&b))), s(&format!("gfx_reduce({}, {mp}, {p})", poly(&a))));
                    let n = s(&format!("gfx_neg({}, {mp}, {p})", poly(&a)));
                    assert_eq!(s(&format!("gfx_add({}, {n}, {mp}, {p})", poly(&a))), "list()");
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
        assert_eq!(s("gfx_mul(list(1, 0), list(1, 0), list(1, 1, 1), 2)"), "list(1, 1)");
        assert_eq!(s("gfx_reduce(list(1, 0, 0), list(1, 1, 1), 2)"), "list(1, 1)");
        // with a reducible modulus a zero divisor has no inverse
        assert_eq!(s("gfx_inv(list(1, 1), list(1, 0, 1), 2)"), "gfx_inv(list(1, 1), list(1, 0, 1), 2)");
        // zero has no inverse
        assert_eq!(s("gfx_inv(list(), list(1, 1, 1), 2)"), "gfx_inv(list(), list(1, 1, 1), 2)");
        // the multiplicative group of GF(8) is cyclic of order 7
        assert_eq!(s("gfx_pow(list(1, 0), 7, list(1, 1, 0, 1), 2)"), "list(1)");
    }
}
