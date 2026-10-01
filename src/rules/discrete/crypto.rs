//! Public-key cryptography on exact integers: elliptic curves over prime
//! fields (point arithmetic, ECDH, ECDSA, point compression) and RSA with
//! caller-supplied primes.
//!
//! A curve `y^2 = x^3 + a x + b (mod p)` is `list(a, b, p)` with an odd
//! modulus `p` (`ec_curve(a, b, p)` reduces `a` and `b`); an affine point is
//! `list(x, y)` and the point at infinity is `list()`.
//!
//! | operator | value |
//! |---|---|
//! | `ec_curve(a, b, p)` | the curve term with reduced coefficients |
//! | `ec_on_curve(C, P)`, `ec_is_infinity(P)` | booleans |
//! | `ec_x(P)`, `ec_y(P)` | the coordinates of an affine point |
//! | `ec_neg(C, P)`, `ec_double(C, P)`, `ec_add(C, P, Q)`, `ec_mul(C, k, P)` | the group law (`k` may be negative) |
//! | `ec_order(C, P)` | the order of `P` (found by repeated addition, at most 10^6 steps) |
//! | `ecdh_public(C, G, d)`, `ecdh_shared(C, d, Q)` | `d G` and `d Q` |
//! | `ec_compress(P)`, `ec_decompress(C, x, odd)` | `list(x, odd)` with `odd` 0 or 1; the point with that `x` and parity of `y` (any odd prime `p`, by Tonelli-Shanks) |
//! | `ecdsa_sign(C, G, n, d, h, k)` | `list(r, s)` for the hash `h`, private key `d`, nonce `k` (explicit, not random); stays unreduced when `r` or `s` is 0 |
//! | `ecdsa_verify(C, G, n, Q, h, list(r, s))` | boolean |
//! | `rsa_keygen(p, q, e)` | `list(n, e, d)` with `d e = 1 (mod (p-1)(q-1))` |
//! | `rsa_encrypt(m, e, n)`, `rsa_decrypt(c, d, n)` | modular powers |

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::big;
use super::bigs;
use super::def;
use super::items;
use super::mod_inverse;
use super::modulo;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::RuleError;

struct Curve {
    a: BigInt,
    b: BigInt,
    p: BigInt,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Pt {
    Inf,
    Aff(BigInt, BigInt),
}

fn curve(
    g: &Graph,
    n: NodeId,
) -> Option<Curve> {
    let [a, b, p] = bigs(g, n)?.try_into().ok()?;
    if p < BigInt::from(3) || !p.bit(0) {
        return None;
    }
    Some(Curve {
        a: modulo(&a, &p),
        b: modulo(&b, &p),
        p,
    })
}

fn point(
    g: &Graph,
    n: NodeId,
    c: &Curve,
) -> Option<Pt> {
    let coords = bigs(g, n)?;
    match coords.as_slice() {
        | [] => Some(Pt::Inf),
        | [x, y] => Some(Pt::Aff(modulo(x, &c.p), modulo(y, &c.p))),
        | _ => None,
    }
}

fn pt_value(p: Pt) -> V {
    match p {
        | Pt::Inf => V::List(Vec::new()),
        | Pt::Aff(x, y) => V::ints([x, y]),
    }
}

fn on_curve(
    c: &Curve,
    p: &Pt,
) -> bool {
    match p {
        | Pt::Inf => true,
        | Pt::Aff(x, y) => modulo(&(y * y), &c.p) == modulo(&(x * x * x + &c.a * x + &c.b), &c.p),
    }
}

fn neg(
    c: &Curve,
    p: &Pt,
) -> Pt {
    match p {
        | Pt::Inf => Pt::Inf,
        | Pt::Aff(x, y) => Pt::Aff(x.clone(), modulo(&-y, &c.p)),
    }
}

fn double(
    c: &Curve,
    p: &Pt,
) -> Option<Pt> {
    let Pt::Aff(x, y) = p else {
        return Some(Pt::Inf);
    };
    if y.is_zero() {
        return Some(Pt::Inf);
    }
    let m = modulo(
        &((BigInt::from(3) * x * x + &c.a) * mod_inverse(&(BigInt::from(2) * y), &c.p)?),
        &c.p,
    );
    let x3 = modulo(&(&m * &m - x - x), &c.p);
    let y3 = modulo(&(&m * (x - &x3) - y), &c.p);
    Some(Pt::Aff(x3, y3))
}

fn add(
    c: &Curve,
    p: &Pt,
    q: &Pt,
) -> Option<Pt> {
    match (p, q) {
        | (Pt::Inf, r) | (r, Pt::Inf) => Some(r.clone()),
        | (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
            if x1 == x2 {
                if y1 != y2 {
                    return Some(Pt::Inf);
                }
                return double(c, p);
            }
            let m = modulo(&((y2 - y1) * mod_inverse(&modulo(&(x2 - x1), &c.p), &c.p)?), &c.p);
            let x3 = modulo(&(&m * &m - x1 - x2), &c.p);
            let y3 = modulo(&(&m * (x1 - &x3) - y1), &c.p);
            Some(Pt::Aff(x3, y3))
        },
    }
}

fn mul(
    c: &Curve,
    k: &BigInt,
    p: &Pt,
) -> Option<Pt> {
    let base = if k.is_negative() { neg(c, p) } else { p.clone() };
    let mut k = k.abs();
    let mut acc = Pt::Inf;
    let mut step = base;
    while !k.is_zero() {
        if k.bit(0) {
            acc = add(c, &acc, &step)?;
        }
        step = double(c, &step)?;
        k >>= 1;
    }
    Some(acc)
}

/// A square root of `a` modulo the odd prime `p`, by Tonelli-Shanks.
fn sqrt_mod(
    a: &BigInt,
    p: &BigInt,
) -> Option<BigInt> {
    let a = modulo(a, p);
    if a.is_zero() {
        return Some(a);
    }
    let one = BigInt::one();
    let half = (p - &one) >> 1;
    if a.modpow(&half, p) != one {
        return None;
    }
    if (p % 4u32) == BigInt::from(3) {
        return Some(a.modpow(&((p + &one) >> 2), p));
    }
    let s = (p - &one).trailing_zeros()?;
    let q = (p - &one) >> s;
    let mut z = BigInt::from(2);
    while z.modpow(&half, p) == one {
        z += 1;
    }
    let mut m = s;
    let mut cc = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + &one) >> 1), p);
    while t != one {
        let mut i = 0;
        let mut tt = t.clone();
        while tt != one {
            tt = &tt * &tt % p;
            i += 1;
            if i == m {
                return None;
            }
        }
        let b = cc.modpow(&(BigInt::one() << (m - i - 1)), p);
        m = i;
        cc = &b * &b % p;
        t = t * &cc % p;
        r = r * &b % p;
    }
    Some(r)
}

fn ec_curve(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [ca, cb, cp] = a else { return None };
    let term = cx.graph.node(crate::graph::op::core::LIST, &[*ca, *cb, *cp]);
    let c = curve(cx.graph, term)?;
    Some(V::ints([c.a, c.b, c.p]))
}

fn with_curve(
    cx: &Cx<'_>,
    args: &[NodeId],
) -> Option<(Curve, Vec<NodeId>)> {
    let (first, rest) = args.split_first()?;
    Some((curve(cx.graph, *first)?, rest.to_vec()))
}

fn ec_on_curve(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let p = point(cx.graph, *rest.first()?, &c)?;
    Some(V::Bool(on_curve(&c, &p)))
}

fn ec_neg(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let p = point(cx.graph, *rest.first()?, &c)?;
    Some(pt_value(neg(&c, &p)))
}

fn ec_double(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let p = point(cx.graph, *rest.first()?, &c)?;
    double(&c, &p).map(pt_value)
}

fn ec_add(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [p, q] = rest.as_slice() else { return None };
    let (p, q) = (point(cx.graph, *p, &c)?, point(cx.graph, *q, &c)?);
    add(&c, &p, &q).map(pt_value)
}

fn ec_mul(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [k, p] = rest.as_slice() else { return None };
    let k = big(cx.graph, *k)?;
    let p = point(cx.graph, *p, &c)?;
    mul(&c, &k, &p).map(pt_value)
}

fn ec_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let p = point(cx.graph, *rest.first()?, &c)?;
    if p == Pt::Inf {
        return Some(V::int(1));
    }
    let mut acc = p.clone();
    for k in 1..=1_000_000_u32 {
        if acc == Pt::Inf {
            return Some(V::Int(BigInt::from(k)));
        }
        acc = add(&c, &acc, &p)?;
    }
    None
}

fn ec_is_infinity(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let list = items(cx.graph, *a.first()?)?;
    Some(V::Bool(list.is_empty()))
}

fn ec_coord(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    which: usize,
) -> Option<V> {
    let list = items(cx.graph, *a.first()?)?;
    (list.len() == 2).then(|| V::Node(list[which]))
}

fn ecdh_public(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [g, d] = rest.as_slice() else { return None };
    let g = point(cx.graph, *g, &c)?;
    let d = big(cx.graph, *d)?;
    mul(&c, &d, &g).map(pt_value)
}

fn ecdh_shared(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [d, q] = rest.as_slice() else { return None };
    let d = big(cx.graph, *d)?;
    let q = point(cx.graph, *q, &c)?;
    mul(&c, &d, &q).map(pt_value)
}

fn ec_compress(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let coords = bigs(cx.graph, *a.first()?)?;
    let [x, y] = coords.as_slice() else { return None };
    Some(V::ints([x.clone(), BigInt::from(u8::from(y.bit(0)))]))
}

fn ec_decompress(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [x, odd] = rest.as_slice() else { return None };
    let x = modulo(&big(cx.graph, *x)?, &c.p);
    let odd = big(cx.graph, *odd)?;
    let odd = if odd.is_zero() {
        false
    } else if odd.is_one() {
        true
    } else {
        return None;
    };
    let rhs = &x * &x * &x + &c.a * &x + &c.b;
    let y = sqrt_mod(&rhs, &c.p)?;
    let y = if y.bit(0) == odd { y } else { modulo(&-y, &c.p) };
    Some(V::ints([x, y]))
}

fn ecdsa_sign(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [g, n, d, h, k] = rest.as_slice() else { return None };
    let g = point(cx.graph, *g, &c)?;
    let (n, d, h, k) = (
        big(cx.graph, *n)?,
        big(cx.graph, *d)?,
        big(cx.graph, *h)?,
        big(cx.graph, *k)?,
    );
    if !n.is_positive() || !k.is_positive() {
        return None;
    }
    let Pt::Aff(rx, _) = mul(&c, &k, &g)? else {
        return None;
    };
    let r = modulo(&rx, &n);
    if r.is_zero() {
        return None;
    }
    let s = modulo(&(mod_inverse(&k, &n)? * (&h + &r * &d)), &n);
    if s.is_zero() {
        return None;
    }
    Some(V::ints([r, s]))
}

fn ecdsa_verify(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (c, rest) = with_curve(cx, a)?;
    let [g, n, q, h, sig] = rest.as_slice() else { return None };
    let g = point(cx.graph, *g, &c)?;
    let n = big(cx.graph, *n)?;
    let q = point(cx.graph, *q, &c)?;
    let h = big(cx.graph, *h)?;
    let [r, s] = bigs(cx.graph, *sig)?.try_into().ok()?;
    if !n.is_positive() {
        return None;
    }
    let range = |v: &BigInt| v.is_positive() && *v < n;
    if !range(&r) || !range(&s) {
        return Some(V::Bool(false));
    }
    let Some(w) = mod_inverse(&s, &n) else {
        return Some(V::Bool(false));
    };
    let u1 = modulo(&(&h * &w), &n);
    let u2 = modulo(&(&r * &w), &n);
    let point = add(&c, &mul(&c, &u1, &g)?, &mul(&c, &u2, &q)?)?;
    Some(V::Bool(match point {
        | Pt::Inf => false,
        | Pt::Aff(x, _) => modulo(&x, &n) == r,
    }))
}

fn rsa_keygen(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, q, e] = a else { return None };
    let (p, q, e) = (big(cx.graph, *p)?, big(cx.graph, *q)?, big(cx.graph, *e)?);
    if p < BigInt::from(2) || q < BigInt::from(2) || p == q {
        return None;
    }
    let phi = (&p - 1) * (&q - 1);
    let d = mod_inverse(&e, &phi)?;
    Some(V::ints([&p * &q, e, d]))
}

fn rsa_power(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [m, e, n] = a else { return None };
    let (m, e, n) = (big(cx.graph, *m)?, big(cx.graph, *e)?, big(cx.graph, *n)?);
    if e.is_negative() || !n.is_positive() {
        return None;
    }
    Some(V::Int(modulo(&m, &n).modpow(&e, &n)))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "ec_curve", Arity::Fixed(3), ec_curve)?;
    def(i, "ec_on_curve", Arity::Fixed(2), ec_on_curve)?;
    def(i, "ec_is_infinity", Arity::Fixed(1), ec_is_infinity)?;
    def(i, "ec_x", Arity::Fixed(1), |cx, a| ec_coord(cx, a, 0))?;
    def(i, "ec_y", Arity::Fixed(1), |cx, a| ec_coord(cx, a, 1))?;
    def(i, "ec_neg", Arity::Fixed(2), ec_neg)?;
    def(i, "ec_double", Arity::Fixed(2), ec_double)?;
    def(i, "ec_add", Arity::Fixed(3), ec_add)?;
    def(i, "ec_mul", Arity::Fixed(3), ec_mul)?;
    def(i, "ec_order", Arity::Fixed(2), ec_order)?;
    def(i, "ecdh_public", Arity::Fixed(3), ecdh_public)?;
    def(i, "ecdh_shared", Arity::Fixed(3), ecdh_shared)?;
    def(i, "ec_compress", Arity::Fixed(1), ec_compress)?;
    def(i, "ec_decompress", Arity::Fixed(3), ec_decompress)?;
    def(i, "ecdsa_sign", Arity::Fixed(6), ecdsa_sign)?;
    def(i, "ecdsa_verify", Arity::Fixed(6), ecdsa_verify)?;
    def(i, "rsa_keygen", Arity::Fixed(3), rsa_keygen)?;
    def(i, "rsa_encrypt", Arity::Fixed(3), rsa_power)?;
    def(i, "rsa_decrypt", Arity::Fixed(3), rsa_power)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::discrete;
    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[discrete()], src)
    }

    /// Independent affine arithmetic on `i64` for a small prime field.
    fn naive_add(
        (a, p): (i64, i64),
        p1: Option<(i64, i64)>,
        p2: Option<(i64, i64)>,
    ) -> Option<(i64, i64)> {
        let inv = |v: i64| (1..p).find(|x| (v.rem_euclid(p) * x) % p == 1).unwrap_or(0);
        let (Some((x1, y1)), Some((x2, y2))) = (p1, p2) else {
            return p1.or(p2);
        };
        let m = if x1 == x2 {
            if (y1 + y2) % p == 0 {
                return None;
            }
            (3 * x1 * x1 + a) % p * inv(2 * y1) % p
        } else {
            (y2 - y1).rem_euclid(p) * inv(x2 - x1) % p
        };
        let x3 = (m * m - x1 - x2).rem_euclid(p);
        let y3 = (m * (x1 - x3) - y1).rem_euclid(p);
        Some((x3, y3))
    }

    fn pt(p: Option<(i64, i64)>) -> String {
        p.map_or("list()".to_string(), |(x, y)| format!("list({x}, {y})"))
    }

    #[test]
    fn group_law_matches_naive_arithmetic() {
        let (a, b, p) = (2_i64, 3_i64, 97_i64);
        let points: Vec<(i64, i64)> = (0..p)
            .flat_map(|x| (0..p).map(move |y| (x, y)))
            .filter(|&(x, y)| (y * y - x * x * x - a * x - b).rem_euclid(p) == 0)
            .collect();
        assert!(points.len() > 50);
        let c = format!("list({a}, {b}, {p})");
        for (n, &p1) in points.iter().enumerate().step_by(7) {
            assert_eq!(s(&format!("ec_on_curve({c}, list({}, {}))", p1.0, p1.1)), "true");
            for &p2 in points.iter().skip(n % 5).step_by(11) {
                let want = pt(naive_add((a, p), Some(p1), Some(p2)));
                assert_eq!(
                    s(&format!("ec_add({c}, list({}, {}), list({}, {}))", p1.0, p1.1, p2.0, p2.1)),
                    want
                );
            }
            let want = pt(naive_add((a, p), Some(p1), Some(p1)));
            assert_eq!(s(&format!("ec_double({c}, list({}, {}))", p1.0, p1.1)), want);
        }
        // off the curve
        assert_eq!(s(&format!("ec_on_curve({c}, list(1, 1))")), "false");
        assert_eq!(s(&format!("ec_on_curve({c}, list())")), "true");
    }

    #[test]
    fn negation_infinity_and_coordinates() {
        assert_eq!(s("ec_neg(list(2, 3, 97), list(3, 6))"), "list(3, 91)");
        assert_eq!(s("ec_neg(list(2, 3, 97), list())"), "list()");
        assert_eq!(s("ec_add(list(2, 3, 97), list(3, 6), list(3, 91))"), "list()");
        assert_eq!(s("ec_add(list(2, 3, 97), list(), list(3, 6))"), "list(3, 6)");
        // a point with y = 0 has order two
        assert_eq!(s("ec_double(list(2, 3, 97), list(0, 0))"), "list()");
        assert_eq!(s("ec_is_infinity(list())"), "true");
        assert_eq!(s("ec_is_infinity(list(3, 6))"), "false");
        assert_eq!(s("ec_x(list(3, 6))"), "3");
        assert_eq!(s("ec_y(list(3, 6))"), "6");
        assert_eq!(s("ec_curve(-1, 100, 97)"), "list(96, 3, 97)");
        // even or tiny moduli are not curves
        let (text, _) = reduce_with(&[discrete()], "ec_curve(1, 1, 8)", &[]);
        assert_eq!(text, "ec_curve(1, 1, 8)");
    }

    #[test]
    fn scalar_multiplication_by_definition() {
        let c = "list(2, 2, 17)";
        let g = "list(5, 1)";
        let mut acc = "list()".to_string();
        for k in 0..=40 {
            assert_eq!(s(&format!("ec_mul({c}, {k}, {g})")), acc, "k = {k}");
            acc = s(&format!("ec_add({c}, {acc}, {g})"));
        }
        // negative scalars
        assert_eq!(s(&format!("ec_mul({c}, -3, {g})")), s(&format!("ec_neg({c}, ec_mul({c}, 3, {g}))")));
        // the textbook order is 19
        assert_eq!(s(&format!("ec_order({c}, {g})")), "19");
        assert_eq!(s(&format!("ec_mul({c}, 19, {g})")), "list()");
        assert_eq!(s(&format!("ec_order({c}, list())")), "1");
        // associativity of the group law
        let p = s(&format!("ec_mul({c}, 4, {g})"));
        let q = s(&format!("ec_mul({c}, 7, {g})"));
        let lhs = s(&format!("ec_add({c}, ec_add({c}, {g}, {p}), {q})"));
        let rhs = s(&format!("ec_add({c}, {g}, ec_add({c}, {p}, {q}))"));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn diffie_hellman_agrees() {
        let (c, g) = ("list(2, 2, 17)", "list(5, 1)");
        for (a, b) in [(3, 7), (11, 5), (2, 18)] {
            let qa = s(&format!("ecdh_public({c}, {g}, {a})"));
            let qb = s(&format!("ecdh_public({c}, {g}, {b})"));
            let sa = s(&format!("ecdh_shared({c}, {a}, {qb})"));
            let sb = s(&format!("ecdh_shared({c}, {b}, {qa})"));
            assert_eq!(sa, sb);
            assert_eq!(sa, s(&format!("ec_mul({c}, {}, {g})", a * b)));
        }
    }

    #[test]
    fn point_compression_round_trips() {
        // p = 97 is 1 mod 4 (Tonelli-Shanks), p = 103 is 3 mod 4, p = 5 is tiny.
        for (a, b, p) in [(2_i64, 3_i64, 97_i64), (2, 3, 103), (1, 1, 5)] {
            let c = format!("list({a}, {b}, {p})");
            let mut found = 0;
            for x in 0..p {
                let has_point = (0..p).any(|y| (y * y - x * x * x - a * x - b).rem_euclid(p) == 0);
                let probe = s(&format!("ec_decompress({c}, {x}, 0)"));
                assert_eq!(!probe.starts_with("ec_decompress"), has_point, "x = {x}");
                if !has_point {
                    continue;
                }
                found += 1;
                for odd in 0..2 {
                    let pt = s(&format!("ec_decompress({c}, {x}, {odd})"));
                    assert_eq!(s(&format!("ec_on_curve({c}, {pt})")), "true");
                    // a point with y = 0 is its own negative, both parities give it
                    let expected = if pt.ends_with(", 0)") { 0 } else { odd };
                    assert_eq!(s(&format!("ec_compress({pt})")), format!("list({x}, {expected})"));
                }
            }
            assert!(found > 0);
        }
        assert_eq!(s("ec_decompress(list(2, 3, 97), 3, 0)"), "list(3, 6)");
        assert_eq!(s("ec_decompress(list(2, 3, 97), 3, 1)"), "list(3, 91)");
    }

    #[test]
    fn ecdsa_signatures_verify() {
        let (c, g, n) = ("list(2, 2, 17)", "list(5, 1)", 19);
        let d = 7;
        let q = s(&format!("ecdh_public({c}, {g}, {d})"));
        for (h, k) in [(5, 3), (11, 8), (1, 18), (0, 2)] {
            let sig = s(&format!("ecdsa_sign({c}, {g}, {n}, {d}, {h}, {k})"));
            if sig.starts_with("ecdsa_sign") {
                continue;
            }
            assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {q}, {h}, {sig})")), "true", "h={h} k={k}");
            // a different message, a different key and a mangled signature fail
            let wrong_h = (h + 1) % n;
            let _ = wrong_h;
        }
        let sig = s(&format!("ecdsa_sign({c}, {g}, {n}, {d}, 5, 3)"));
        let rs: Vec<i64> = sig
            .trim_start_matches("list(")
            .trim_end_matches(')')
            .split(", ")
            .filter_map(|t| t.parse().ok())
            .collect();
        assert_eq!(rs.len(), 2);
        let bad_s = rs[1] % 18 + 1;
        assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {q}, 5, list({}, {bad_s}))", rs[0])), "false");
        let other = s(&format!("ecdh_public({c}, {g}, 5)"));
        assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {other}, 5, {sig})")), "false");
        assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {q}, 6, {sig})")), "false");
        // out-of-range components
        assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {q}, 5, list(0, 3))")), "false");
        assert_eq!(s(&format!("ecdsa_verify({c}, {g}, {n}, {q}, 5, list(3, 19))")), "false");
    }

    #[test]
    fn rsa_with_given_primes() {
        assert_eq!(s("rsa_keygen(61, 53, 17)"), "list(3233, 17, 2753)");
        assert_eq!(s("rsa_encrypt(65, 17, 3233)"), "2790");
        assert_eq!(s("rsa_decrypt(2790, 2753, 3233)"), "65");
        for m in [0, 1, 2, 1000, 3232] {
            let c = s(&format!("rsa_encrypt({m}, 17, 3233)"));
            assert_eq!(s(&format!("rsa_decrypt({c}, 2753, 3233)")), m.to_string());
        }
        // e not coprime to phi: no key
        assert_eq!(s("rsa_keygen(61, 53, 4)"), "rsa_keygen(61, 53, 4)");
        // big primes use exact integers
        assert_eq!(
            s("rsa_decrypt(rsa_encrypt(123456789, 65537, 1000000007 * 998244353), \
               invmod(65537, 1000000006 * 998244352), 1000000007 * 998244353)"),
            "123456789"
        );
    }

    #[test]
    fn symbolic_arguments_stay() {
        let sets = [discrete()];
        for src in ["ec_add(list(2, 3, 97), list(x, 6), list(3, 6))", "ec_mul(list(2, 3, p), 3, list(3, 6))", "rsa_encrypt(m, 17, 3233)"] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(reduced);
            assert_eq!(text, src);
        }
    }
}
