//! Error-correcting codes and checksums on lists of integers.
//!
//! Words are `list(...)` of bytes (bits for Hamming and BCH words). The
//! work is done by [`crate::kernels::error_correction`] (Hamming(7,4),
//! Reed-Solomon over GF(2^8) with generator `x - alpha^i`, `i = 0..n-1`,
//! BCH, CRC-8/16/32) and [`crate::kernels::finite_field`] (GF(2^8) with the
//! polynomial `0x11d`). Polynomials over GF(256) are byte lists with the
//! **highest degree first**.
//!
//! | operator | value |
//! |---|---|
//! | `hamming_distance(a, b)`, `hamming_weight(a)` | number of differing positions (equal lengths) / non-zero symbols |
//! | `hamming_encode(d)`, `hamming_check(c)`, `hamming_decode(c)` | Hamming(7,4) with parity at positions 1, 2, 4; `decode` is `list(data, position)` with the 1-based position of the corrected bit, `0` when there was no error |
//! | `rs_encode(data, n)`, `rs_check(c, n)`, `rs_decode(c, n)`, `rs_error_count(c, n)` | systematic Reed-Solomon with `n` parity symbols; `decode` returns the data part and stays unreduced when more than `n/2` symbols are wrong; `error_count` is the degree of the Berlekamp-Massey locator |
//! | `bch_encode(data, t)`, `bch_decode(c, t)` | the kernel's simplified BCH with `2 t` parity bits |
//! | `crc32(data)`, `crc32_verify(data, crc)`, `crc32_update(crc, data)`, `crc32_finalize(crc)`, `crc16(data)`, `crc8(data)` | checksums; `crc32(d) = crc32_finalize(crc32_update(4294967295, d))` |
//! | `interleave(data, depth)`, `deinterleave(data, depth)`, `conv_encode(data)` | block interleaving, rate 1/2 convolutional code |
//! | `code_min_distance(list(c1, c2, ...))`, `code_rate(k, n)` | minimum pairwise distance; `k / n` |
//! | `gf256_add(a, b)`, `gf256_mul`, `gf256_div`, `gf256_inv(a)`, `gf256_pow(a, e)` | GF(2^8) arithmetic |
//! | `gf256_exp(k)`, `gf256_log(a)` | powers of the generator 2 and discrete logarithm |
//! | `gf256_poly_eval(p, x)`, `gf256_poly_add(p, q)`, `gf256_poly_mul`, `gf256_poly_scale(p, c)`, `gf256_poly_derivative(p)` | polynomials over GF(2^8) |
//! | `gf256_poly_mod(p, d)`, `gf256_poly_divmod(p, d)`, `gf256_poly_gcd(p, q)` | remainder (with exactly `deg d` coefficients), `list(quotient, remainder)`, monic gcd |

use num_bigint::BigInt;

use super::big;
use super::byte_list;
use super::bytes;
use super::def;
use super::idx;
use super::items;
use super::small;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Number;
use crate::graph::NodeId;
use crate::graph::RuleError;
use crate::kernels::error_correction as ec;
use crate::kernels::finite_field as ff;

fn two_bytes(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<(Vec<u8>, Vec<u8>)> {
    let [x, y] = a else { return None };
    Some((bytes(cx.graph, *x)?, bytes(cx.graph, *y)?))
}

fn word_count(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<(Vec<u8>, usize)> {
    let [x, n] = a else { return None };
    Some((bytes(cx.graph, *x)?, idx(cx.graph, *n)?))
}

fn hamming_distance(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (x, y) = two_bytes(cx, a)?;
    ec::hamming_distance_numerical(&x, &y).map(V::uint)
}

fn hamming_weight(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = bytes(cx.graph, *a.first()?)?;
    Some(V::uint(ec::hamming_weight_numerical(&x)))
}

fn hamming_encode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = bytes(cx.graph, *a.first()?)?;
    ec::hamming_encode_numerical(&x).map(|c| byte_list(&c))
}

fn hamming_check(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = bytes(cx.graph, *a.first()?)?;
    Some(V::Bool(ec::hamming_check_numerical(&x)))
}

fn hamming_decode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = bytes(cx.graph, *a.first()?)?;
    let (data, pos) = ec::hamming_decode_numerical(&x).ok()?;
    Some(V::List(vec![byte_list(&data), V::uint(pos.unwrap_or(0))]))
}

fn rs_encode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (data, n) = word_count(cx, a)?;
    ec::reed_solomon_encode(&data, n).ok().map(|c| byte_list(&c))
}

fn rs_check(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (word, n) = word_count(cx, a)?;
    (n <= word.len()).then(|| V::Bool(ec::reed_solomon_check(&word, n)))
}

fn rs_decode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (mut word, n) = word_count(cx, a)?;
    if n > word.len() || word.len() > 255 {
        return None;
    }
    ec::reed_solomon_decode(&mut word, n).ok()?;
    Some(byte_list(&word[..word.len() - n]))
}

fn rs_error_count(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (word, n) = word_count(cx, a)?;
    if n > word.len() {
        return None;
    }
    let syndromes = ec::calculate_syndromes(&word, n);
    if syndromes.iter().all(|&s| s == 0) {
        return Some(V::int(0));
    }
    let sigma = ec::berlekamp_massey(&syndromes);
    Some(V::uint(sigma.iter().rposition(|&c| c != 0).unwrap_or(0)))
}

fn bch_encode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (data, t) = word_count(cx, a)?;
    Some(byte_list(&ec::bch_encode(&data, t)))
}

fn bch_decode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (word, t) = word_count(cx, a)?;
    ec::bch_decode(&word, t).ok().map(|d| byte_list(&d))
}

fn crc_value(
    cx: &Cx<'_>,
    n: NodeId,
) -> Option<u32> {
    u32::try_from(big(cx.graph, n)?).ok()
}

fn crc32(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let data = bytes(cx.graph, *a.first()?)?;
    Some(V::Int(BigInt::from(ec::crc32_compute_numerical(&data))))
}

fn crc32_verify(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [d, c] = a else { return None };
    let data = bytes(cx.graph, *d)?;
    Some(V::Bool(ec::crc32_verify_numerical(&data, crc_value(cx, *c)?)))
}

fn crc32_update(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [c, d] = a else { return None };
    let data = bytes(cx.graph, *d)?;
    Some(V::Int(BigInt::from(ec::crc32_update_numerical(crc_value(cx, *c)?, &data))))
}

fn crc32_finalize(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Int(BigInt::from(!crc_value(cx, *a.first()?)?)))
}

fn crc16(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let data = bytes(cx.graph, *a.first()?)?;
    Some(V::Int(BigInt::from(ec::crc16_compute(&data))))
}

fn crc8(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let data = bytes(cx.graph, *a.first()?)?;
    Some(V::Int(BigInt::from(ec::crc8_compute(&data))))
}

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn interleave(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    forward: bool,
) -> Option<V> {
    let (data, depth) = word_count(cx, a)?;
    Some(byte_list(&if forward { ec::interleave(&data, depth) } else { ec::deinterleave(&data, depth) }))
}

fn conv_encode(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let data = bytes(cx.graph, *a.first()?)?;
    Some(byte_list(&ec::convolutional_encode(&data)))
}

fn code_min_distance(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let words: Option<Vec<Vec<u8>>> = items(cx.graph, *a.first()?)?
        .into_iter()
        .map(|w| bytes(cx.graph, w))
        .collect();
    ec::minimum_distance(&words?).map(V::uint)
}

fn code_rate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [k, n] = a else { return None };
    let (k, n) = (small(cx.graph, *k)?, small(cx.graph, *n)?);
    let rate = Number::fraction(k, n)?;
    Some(V::Node(cx.graph.num(rate)))
}

// ---------------- GF(256) ----------------

fn byte(
    cx: &Cx<'_>,
    n: NodeId,
) -> Option<u8> {
    u8::try_from(small(cx.graph, n)?).ok()
}

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn gf256_binary(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    op: fn(u8, u8) -> Option<u8>,
) -> Option<V> {
    let [x, y] = a else { return None };
    op(byte(cx, *x)?, byte(cx, *y)?).map(|v| V::int(i64::from(v)))
}

fn gf256_inv(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = byte(cx, *a.first()?)?;
    ff::gf256_inv(x).ok().map(|v| V::int(i64::from(v)))
}

fn gf256_pow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, e] = a else { return None };
    let (x, e) = (byte(cx, *x)?, u64::try_from(big(cx.graph, *e)?).ok()?);
    Some(V::int(i64::from(ff::gf256_pow(x, e))))
}

fn gf256_exp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let k = u64::try_from(big(cx.graph, *a.first()?)?).ok()?;
    Some(V::int(i64::from(ff::gf256_pow(2, k))))
}

fn gf256_log(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let x = byte(cx, *a.first()?)?;
    (1..255_u64).find(|&k| ff::gf256_pow(2, k) == x).or_else(|| (x == 1).then_some(0)).map(|k| V::Int(BigInt::from(k)))
}

type Bytes = Vec<u8>;

fn strip(mut p: Bytes) -> Bytes {
    let first = p.iter().position(|&c| c != 0).unwrap_or(p.len());
    p.drain(..first);
    p
}

fn padd(
    p: &[u8],
    q: &[u8],
) -> Bytes {
    let n = p.len().max(q.len());
    let mut out = vec![0_u8; n];
    for (i, &c) in p.iter().rev().enumerate() {
        out[n - 1 - i] ^= c;
    }
    for (i, &c) in q.iter().rev().enumerate() {
        out[n - 1 - i] ^= c;
    }
    out
}

fn pmul(
    p: &[u8],
    q: &[u8],
) -> Bytes {
    if p.is_empty() || q.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0_u8; p.len() + q.len() - 1];
    for (i, &a) in p.iter().enumerate() {
        for (j, &b) in q.iter().enumerate() {
            out[i + j] ^= ff::gf256_mul(a, b);
        }
    }
    out
}

/// `(quotient, remainder)` with the remainder padded to `len(d) - 1`
/// coefficients (the shape of a long division loop); the divisor's leading
/// coefficient must be non-zero.
fn pdivmod(
    p: &[u8],
    d: &[u8],
) -> Option<(Bytes, Bytes)> {
    let lead = ff::gf256_inv(*d.first()?).ok()?;
    let mut rem = p.to_vec();
    let mut quot = Vec::new();
    while rem.len() >= d.len() {
        let coeff = ff::gf256_mul(rem[0], lead);
        for (j, &c) in d.iter().enumerate() {
            rem[j] ^= ff::gf256_mul(coeff, c);
        }
        quot.push(coeff);
        rem.remove(0);
    }
    Some((quot, rem))
}

fn pgcd(
    p: &[u8],
    q: &[u8],
) -> Bytes {
    let (mut a, mut b) = (strip(p.to_vec()), strip(q.to_vec()));
    while !b.is_empty() {
        let Some((_, r)) = pdivmod(&a, &b) else { break };
        a = b;
        b = strip(r);
    }
    if let Some(&lead) = a.first()
        && let Ok(inv) = ff::gf256_inv(lead) {
            a = a.iter().map(|&c| ff::gf256_mul(c, inv)).collect();
        }
    a
}

#[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
fn poly_pair(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    op: fn(&[u8], &[u8]) -> Option<V>,
) -> Option<V> {
    let (p, q) = two_bytes(cx, a)?;
    op(&p, &q)
}

fn poly_eval(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, x] = a else { return None };
    let (p, x) = (bytes(cx.graph, *p)?, byte(cx, *x)?);
    let y = p.iter().fold(0_u8, |acc, &c| ff::gf256_mul(acc, x) ^ c);
    Some(V::int(i64::from(y)))
}

fn poly_scale(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, c] = a else { return None };
    let (p, c) = (bytes(cx.graph, *p)?, byte(cx, *c)?);
    Some(byte_list(&p.iter().map(|&x| ff::gf256_mul(x, c)).collect::<Vec<_>>()))
}

fn poly_derivative(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let p = bytes(cx.graph, *a.first()?)?;
    if p.len() <= 1 {
        return Some(byte_list(&[0]));
    }
    let n = p.len() - 1;
    // In characteristic 2 only odd powers survive.
    let mut out: Vec<u8> = p.iter().enumerate().take(n).map(|(i, &c)| if (n - i) % 2 == 1 { c } else { 0 }).collect();
    while out.len() > 1 && out[0] == 0 {
        out.remove(0);
    }
    Some(byte_list(&out))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "hamming_distance", Arity::Fixed(2), hamming_distance)?;
    def(i, "hamming_weight", Arity::Fixed(1), hamming_weight)?;
    def(i, "hamming_encode", Arity::Fixed(1), hamming_encode)?;
    def(i, "hamming_check", Arity::Fixed(1), hamming_check)?;
    def(i, "hamming_decode", Arity::Fixed(1), hamming_decode)?;
    def(i, "rs_encode", Arity::Fixed(2), rs_encode)?;
    def(i, "rs_check", Arity::Fixed(2), rs_check)?;
    def(i, "rs_decode", Arity::Fixed(2), rs_decode)?;
    def(i, "rs_error_count", Arity::Fixed(2), rs_error_count)?;
    def(i, "bch_encode", Arity::Fixed(2), bch_encode)?;
    def(i, "bch_decode", Arity::Fixed(2), bch_decode)?;
    def(i, "crc32", Arity::Fixed(1), crc32)?;
    def(i, "crc32_verify", Arity::Fixed(2), crc32_verify)?;
    def(i, "crc32_update", Arity::Fixed(2), crc32_update)?;
    def(i, "crc32_finalize", Arity::Fixed(1), crc32_finalize)?;
    def(i, "crc16", Arity::Fixed(1), crc16)?;
    def(i, "crc8", Arity::Fixed(1), crc8)?;
    def(i, "interleave", Arity::Fixed(2), |cx, a| interleave(cx, a, true))?;
    def(i, "deinterleave", Arity::Fixed(2), |cx, a| interleave(cx, a, false))?;
    def(i, "conv_encode", Arity::Fixed(1), conv_encode)?;
    def(i, "code_min_distance", Arity::Fixed(1), code_min_distance)?;
    def(i, "code_rate", Arity::Fixed(2), code_rate)?;

    def(i, "gf256_add", Arity::Fixed(2), |cx, a| gf256_binary(cx, a, |x, y| Some(ff::gf256_add(x, y))))?;
    def(i, "gf256_mul", Arity::Fixed(2), |cx, a| gf256_binary(cx, a, |x, y| Some(ff::gf256_mul(x, y))))?;
    def(i, "gf256_div", Arity::Fixed(2), |cx, a| gf256_binary(cx, a, |x, y| ff::gf256_div(x, y).ok()))?;
    def(i, "gf256_inv", Arity::Fixed(1), gf256_inv)?;
    def(i, "gf256_pow", Arity::Fixed(2), gf256_pow)?;
    def(i, "gf256_exp", Arity::Fixed(1), gf256_exp)?;
    def(i, "gf256_log", Arity::Fixed(1), gf256_log)?;
    def(i, "gf256_poly_eval", Arity::Fixed(2), poly_eval)?;
    def(i, "gf256_poly_add", Arity::Fixed(2), |cx, a| poly_pair(cx, a, |p, q| Some(byte_list(&padd(p, q)))))?;
    def(i, "gf256_poly_mul", Arity::Fixed(2), |cx, a| poly_pair(cx, a, |p, q| Some(byte_list(&pmul(p, q)))))?;
    def(i, "gf256_poly_scale", Arity::Fixed(2), poly_scale)?;
    def(i, "gf256_poly_derivative", Arity::Fixed(1), poly_derivative)?;
    def(i, "gf256_poly_mod", Arity::Fixed(2), |cx, a| {
        poly_pair(cx, a, |p, d| pdivmod(p, d).map(|(_, r)| byte_list(&r)))
    })?;
    def(i, "gf256_poly_divmod", Arity::Fixed(2), |cx, a| {
        poly_pair(cx, a, |p, d| pdivmod(p, d).map(|(q, r)| V::List(vec![byte_list(&q), byte_list(&r)])))
    })?;
    def(i, "gf256_poly_gcd", Arity::Fixed(2), |cx, a| poly_pair(cx, a, |p, q| Some(byte_list(&pgcd(p, q)))))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::nested;
    use super::super::test_util::nums;
    use super::super::test_util::s;

    fn word(w: &[i64]) -> String {
        format!("list({})", w.iter().map(ToString::to_string).collect::<Vec<_>>().join(", "))
    }

    /// A small deterministic generator.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self, bound: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % bound
        }
    }

    #[test]
    fn hamming_distance_and_weight() {
        assert_eq!(s("hamming_distance(list(1, 0, 1, 1), list(1, 1, 1, 0))"), "2");
        assert_eq!(s("hamming_distance(list(), list())"), "0");
        assert_eq!(s("hamming_distance(list(1), list(1, 0))"), "hamming_distance(list(1), list(1, 0))");
        assert_eq!(s("hamming_weight(list(1, 0, 1, 1))"), "3");
        assert_eq!(s("hamming_weight(list(0, 0, 0))"), "0");
        assert_eq!(s("hamming_weight(list(1, x))"), "hamming_weight(list(1, x))");
    }

    #[test]
    fn hamming_code_corrects_every_single_error() {
        for m in 0..16_i64 {
            let data = [(m >> 3) & 1, (m >> 2) & 1, (m >> 1) & 1, m & 1];
            let cw = nums(&s(&format!("hamming_encode({})", word(&data))));
            assert_eq!(cw.len(), 7);
            assert_eq!(s(&format!("hamming_check({})", word(&cw))), "true");
            assert_eq!(s(&format!("hamming_decode({})", word(&cw))), format!("list({}, 0)", word(&data)));
            for pos in 0..7 {
                let mut bad = cw.clone();
                bad[pos] ^= 1;
                assert_eq!(s(&format!("hamming_check({})", word(&bad))), "false");
                assert_eq!(
                    s(&format!("hamming_decode({})", word(&bad))),
                    format!("list({}, {})", word(&data), pos + 1),
                    "m = {m}, pos = {pos}"
                );
            }
        }
        // minimum distance 3: all 16 codewords
        let words: Vec<String> = (0..16_i64)
            .map(|m| s(&format!("hamming_encode({})", word(&[(m >> 3) & 1, (m >> 2) & 1, (m >> 1) & 1, m & 1]))))
            .collect();
        assert_eq!(s(&format!("code_min_distance(list({}))", words.join(", "))), "3");
        assert_eq!(s("hamming_encode(list(1, 0, 1))"), "hamming_encode(list(1, 0, 1))");
    }

    #[test]
    fn reed_solomon_round_trips_with_up_to_t_errors() {
        let mut rng = Lcg(7);
        for (len, nsym) in [(3_usize, 4_usize), (5, 6), (10, 8), (20, 4), (1, 2)] {
            for _ in 0..6 {
                let data: Vec<i64> = (0..len).map(|_| rng.next(256) as i64).collect();
                let cw = nums(&s(&format!("rs_encode({}, {nsym})", word(&data))));
                assert_eq!(cw.len(), len + nsym);
                assert_eq!(&cw[..len], data.as_slice(), "systematic");
                assert_eq!(s(&format!("rs_check({}, {nsym})", word(&cw))), "true");
                assert_eq!(s(&format!("rs_error_count({}, {nsym})", word(&cw))), "0");
                for errors in 1..=nsym / 2 {
                    let mut bad = cw.clone();
                    let mut positions = Vec::new();
                    while positions.len() < errors {
                        let p = rng.next(bad.len() as u64) as usize;
                        if !positions.contains(&p) {
                            positions.push(p);
                        }
                    }
                    for &p in &positions {
                        bad[p] ^= 1 + rng.next(255) as i64;
                    }
                    assert_eq!(s(&format!("rs_check({}, {nsym})", word(&bad))), "false");
                    assert_eq!(s(&format!("rs_error_count({}, {nsym})", word(&bad))), errors.to_string());
                    assert_eq!(
                        s(&format!("rs_decode({}, {nsym})", word(&bad))),
                        word(&data),
                        "len {len}, nsym {nsym}, errors at {positions:?}"
                    );
                }
            }
        }
        // every single-symbol error of one codeword, at every position and value class
        let data = [72, 101, 108, 108, 111];
        let cw = nums(&s(&format!("rs_encode({}, 4)", word(&data))));
        for p in 0..cw.len() {
            for e in [1, 2, 0x80, 0xff] {
                let mut bad = cw.clone();
                bad[p] ^= e;
                assert_eq!(s(&format!("rs_decode({}, 4)", word(&bad))), word(&data));
            }
        }
    }

    #[test]
    fn reed_solomon_refuses_to_miscorrect() {
        // Three errors with four parity symbols: either stays unreduced or
        // (by chance) decodes to another valid codeword's data.
        let data = [1, 2, 3, 4, 5, 6];
        let cw = nums(&s(&format!("rs_encode({}, 4)", word(&data))));
        let mut bad = cw;
        bad[0] ^= 5;
        bad[3] ^= 9;
        bad[7] ^= 0x4D;
        let out = s(&format!("rs_decode({}, 4)", word(&bad)));
        assert!(out.starts_with("rs_decode") || out != word(&data));
        assert_eq!(s("rs_encode(list(1), 255)"), "rs_encode(list(1), 255)");
        assert_eq!(s("rs_decode(list(1, 2), 3)"), "rs_decode(list(1, 2), 3)");
        assert_eq!(s("rs_encode(list(1, 2, 256), 2)"), "rs_encode(list(1, 2, 256), 2)");
    }

    #[test]
    fn bch_round_trip() {
        let data = [1, 0, 1, 1, 0, 0, 1, 0];
        let cw = nums(&s(&format!("bch_encode({}, 2)", word(&data))));
        assert_eq!(cw.len(), data.len() + 4);
        assert_eq!(&cw[..8], &data);
        assert_eq!(s(&format!("bch_decode({}, 2)", word(&cw))), word(&data));
        // a single data-bit error is corrected
        for p in 0..data.len() {
            let mut bad = cw.clone();
            bad[p] ^= 1;
            assert_eq!(s(&format!("bch_decode({}, 2)", word(&bad))), word(&data), "bit {p}");
        }
    }

    /// Bitwise reference CRC-32 (IEEE 802.3).
    fn crc32_reference(data: &[u8]) -> u32 {
        let mut crc = 0xFFFF_FFFF_u32;
        for &b in data {
            crc ^= u32::from(b);
            for _ in 0..8 {
                crc = if crc & 1 == 1 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
            }
        }
        !crc
    }

    #[test]
    fn crc_checksums() {
        // the standard check value of CRC-32
        let digits: Vec<i64> = b"123456789".iter().map(|&b| i64::from(b)).collect();
        assert_eq!(s(&format!("crc32({})", word(&digits))), "3421780262");
        assert_eq!(s(&format!("crc32_verify({}, 3421780262)", word(&digits))), "true");
        assert_eq!(s(&format!("crc32_verify({}, 3421780263)", word(&digits))), "false");
        assert_eq!(s("crc32(list())"), "0");
        let mut rng = Lcg(3);
        for len in [1_usize, 2, 7, 33] {
            let data: Vec<u8> = (0..len).map(|_| rng.next(256) as u8).collect();
            let list: Vec<i64> = data.iter().map(|&b| i64::from(b)).collect();
            assert_eq!(s(&format!("crc32({})", word(&list))), crc32_reference(&data).to_string());
            // streaming: update over two chunks then finalize
            let cut = len / 2;
            let first = s(&format!("crc32_update(4294967295, {})", word(&list[..cut])));
            let second = s(&format!("crc32_update({first}, {})", word(&list[cut..])));
            assert_eq!(s(&format!("crc32_finalize({second})")), crc32_reference(&data).to_string());
        }
        assert_eq!(s("crc32_finalize(0)"), "4294967295");
        // CRC-16 and CRC-8 are deterministic and sensitive to every byte
        let base = s("crc16(list(1, 2, 3, 4))");
        assert_ne!(base, s("crc16(list(1, 2, 3, 5))"));
        assert_eq!(base, s("crc16(list(1, 2, 3, 4))"));
        let base = s("crc8(list(1, 2, 3, 4))");
        assert_ne!(base, s("crc8(list(1, 2, 3, 5))"));
        assert_eq!(s("crc32(list(300))"), "crc32(list(300))");
    }

    #[test]
    fn interleaving_and_convolutional_code() {
        let data = [1, 2, 3, 4, 5, 6, 7, 8, 9];
        let il = s(&format!("interleave({}, 3)", word(&data)));
        assert_ne!(il, word(&data));
        assert_eq!(s(&format!("deinterleave({il}, 3)")), word(&data));
        assert_eq!(nums(&s("conv_encode(list(1, 0, 1, 1))")).len() % 2, 0);
        assert_eq!(s("code_rate(4, 7)"), "4/7");
        assert_eq!(s("code_rate(1, 2)"), "1/2");
    }

    #[test]
    fn gf256_field_axioms() {
        assert_eq!(s("gf256_add(53, 202)"), "255");
        assert_eq!(s("gf256_mul(2, 128)"), "29"); // x^8 = x^4 + x^3 + x^2 + 1
        assert_eq!(s("gf256_exp(0)"), "1");
        assert_eq!(s("gf256_exp(1)"), "2");
        assert_eq!(s("gf256_exp(8)"), "29");
        assert_eq!(s("gf256_exp(255)"), "1");
        assert_eq!(s("gf256_log(2)"), "1");
        assert_eq!(s("gf256_log(1)"), "0");
        assert_eq!(s("gf256_log(0)"), "gf256_log(0)");
        assert_eq!(s("gf256_inv(0)"), "gf256_inv(0)");
        assert_eq!(s("gf256_div(5, 0)"), "gf256_div(5, 0)");
        for a in 1..=255_i64 {
            let inv = s(&format!("gf256_inv({a})"));
            assert_eq!(s(&format!("gf256_mul({a}, {inv})")), "1");
            assert_eq!(s(&format!("gf256_div({a}, {a})")), "1");
            assert_eq!(s(&format!("gf256_pow({a}, 255)")), "1");
            let log = s(&format!("gf256_log({a})"));
            assert_eq!(s(&format!("gf256_exp({log})")), a.to_string());
        }
        // shift-and-xor reference multiplication
        let reference = |mut a: u32, mut b: u32| {
            let mut r = 0;
            while b > 0 {
                if b & 1 == 1 {
                    r ^= a;
                }
                a <<= 1;
                if a & 0x100 != 0 {
                    a ^= 0x11d;
                }
                b >>= 1;
            }
            r
        };
        let mut rng = Lcg(11);
        for _ in 0..40 {
            let (a, b) = (rng.next(256) as u32, rng.next(256) as u32);
            assert_eq!(s(&format!("gf256_mul({a}, {b})")), reference(a, b).to_string());
        }
        assert_eq!(s("gf256_pow(0, 0)"), "1");
        assert_eq!(s("gf256_pow(0, 3)"), "0");
    }

    #[test]
    fn gf256_polynomials() {
        // (x + 1)(x + 2) = x^2 + 3x + 2
        assert_eq!(s("gf256_poly_mul(list(1, 1), list(1, 2))"), "list(1, 3, 2)");
        assert_eq!(s("gf256_poly_add(list(1, 2, 3), list(3))"), "list(1, 2, 0)");
        assert_eq!(s("gf256_poly_add(list(1, 2), list(1, 2))"), "list(0, 0)");
        assert_eq!(s("gf256_poly_scale(list(1, 2, 3), 2)"), "list(2, 4, 6)");
        assert_eq!(s("gf256_poly_eval(list(1, 3, 2), 1)"), "0");
        assert_eq!(s("gf256_poly_eval(list(1, 3, 2), 2)"), "0");
        // x^3 + x^2 + x + 1: derivative keeps only odd powers: x^2 + 1
        assert_eq!(s("gf256_poly_derivative(list(1, 1, 1, 1))"), "list(1, 0, 1)");
        assert_eq!(s("gf256_poly_derivative(list(5))"), "list(0)");
        // division: f = q g + r
        let mut rng = Lcg(5);
        for _ in 0..20 {
            let f: Vec<i64> = (0..6).map(|_| rng.next(256) as i64).collect();
            let g: Vec<i64> = std::iter::once(1 + rng.next(255) as i64).chain((0..2).map(|_| rng.next(256) as i64)).collect();
            let qr = nested(&s(&format!("gf256_poly_divmod({}, {})", word(&f), word(&g))));
            let qg = s(&format!("gf256_poly_mul({}, {})", word(&qr[0]), word(&g)));
            let sum = nums(&s(&format!("gf256_poly_add({qg}, {})", word(&qr[1]))));
            assert_eq!(sum, f);
            assert_eq!(qr[1].len(), g.len() - 1);
            assert_eq!(nums(&s(&format!("gf256_poly_mod({}, {})", word(&f), word(&g)))), qr[1]);
        }
        // gcd of (x+1)(x+2) and (x+1)(x+4) is x + 1
        let a = s("gf256_poly_mul(list(1, 1), list(1, 2))");
        let b = s("gf256_poly_mul(list(1, 1), list(1, 4))");
        assert_eq!(s(&format!("gf256_poly_gcd({a}, {b})")), "list(1, 1)");
        assert_eq!(s("gf256_poly_gcd(list(0, 0, 3, 3), list(2, 2))"), "list(1, 1)");
        // the Reed-Solomon generator has roots alpha^0 .. alpha^(n-1)
        let mut g = "list(1)".to_string();
        for k in 0..4 {
            let root = s(&format!("gf256_exp({k})"));
            g = s(&format!("gf256_poly_mul({g}, list(1, {root}))"));
        }
        for k in 0..4 {
            let root = s(&format!("gf256_exp({k})"));
            assert_eq!(s(&format!("gf256_poly_eval({g}, {root})")), "0");
        }
    }
}
