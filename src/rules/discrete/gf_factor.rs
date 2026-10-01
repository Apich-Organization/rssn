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

use super::big;
use super::def;
use super::finite_field::Poly;
use super::finite_field::norm;
use super::finite_field::pdivmod;
use super::finite_field::pegcd;
use super::finite_field::pmul;
use super::finite_field::poly_value;
use super::finite_field::psub;
use super::finite_field::read_poly;
use super::finite_field::reduce;
use super::finite_field::xinv;
use super::items;
use super::mod_inverse;
use super::modulo;
use super::pow;
use super::prod;
use super::sum;
use super::V;
use crate::graph::op::core;
use crate::graph::rule::Installer;
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
fn is_constant(f: &Poly) -> bool {
    f.len() == 1
}

fn monic(
    f: &Poly,
    p: &BigInt,
) -> Option<Poly> {
    super::finite_field::monic(f, p)
}

fn derivative(
    f: &Poly,
    p: &BigInt,
) -> Poly {
    let n = f.len();
    norm((0..n.saturating_sub(1)).map(|j| &f[j] * BigInt::from(n - 1 - j)).collect(), p)
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

/// `f^e mod m`.
fn powmod(
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

/// The polynomial `g` with `g(x^p) = c(x)` when only degrees divisible by
/// `p` occur in `c`.
fn p_th_root(
    c: &Poly,
    p: usize,
) -> Poly {
    let degree = c.len() - 1;
    (0..=degree / p).rev().map(|k| c[degree - k * p].clone()).collect()
}

/// Yun's squarefree decomposition of a monic polynomial: pairs (factor,
/// multiplicity), factors monic and pairwise coprime, sorted by
/// multiplicity.
fn squarefree(
    f: &Poly,
    p: &BigInt,
) -> Option<Vec<(Poly, usize)>> {
    let mut out: Vec<(Poly, usize)> = Vec::new();
    if is_constant(f) || f.is_empty() {
        return Some(out);
    }
    let df = derivative(f, p);
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
        for (g, m) in squarefree(&monic(&root, p)?, p)? {
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

/// Distinct-degree factorisation of a monic squarefree polynomial.
fn ddf(
    f: &Poly,
    p: &BigInt,
) -> Option<Vec<(Poly, usize)>> {
    let mut out = Vec::new();
    let mut h = x_poly();
    let mut rest = f.clone();
    let mut d = 1;
    while rest.len() > 2 * d {
        h = powmod(&h, p, &rest, p)?;
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

/// A deterministic generator of residues.
struct Rng(u64);

impl Rng {
    fn word(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
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
}

/// Cantor–Zassenhaus: the irreducible factors of the monic squarefree `g`,
/// all of degree `d`.
fn edf(
    g: &Poly,
    d: usize,
    p: &BigInt,
    rng: &mut Rng,
    out: &mut Vec<Poly>,
) -> Option<()> {
    let n = g.len() - 1;
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
            psub(&powmod(&a, &exponent, g, p)?, &one(), p)
        };
        let h = gcd(&t, g, p)?;
        if !is_constant(&h) && h.len() != g.len() {
            let rest = quotient(g, &h, p)?;
            edf(&h, d, p, rng, out)?;
            return edf(&rest, d, p, rng, out);
        }
    }
    None
}

/// Rabin's irreducibility test (`p` prime).
pub(super) fn is_irreducible(
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
        h = powmod(&h, p, &f, p)?;
        if checkpoints.contains(&k) {
            let g = gcd(&psub(&h, &x, p), &f, p)?;
            if !is_constant(&g) {
                return Some(false);
            }
        }
    }
    Some(reduce(&psub(&h, &x, p), &f, p)?.is_empty())
}

/// The full factorisation: the leading coefficient and (factor,
/// multiplicity) pairs, ordered by degree.
fn factor(
    f: &Poly,
    p: &BigInt,
    split: fn(&Poly, usize, &BigInt, &mut Rng) -> Option<Vec<Poly>>,
) -> Option<(BigInt, Vec<(Poly, usize)>)> {
    let lead = f.first()?.clone();
    let g = monic(f, p)?;
    let mut rng = Rng(0x5eed_1234_abcd_9876);
    let mut out: Vec<(Poly, usize)> = Vec::new();
    for (sq, m) in squarefree(&g, p)? {
        for (gd, d) in ddf(&sq, p)? {
            for q in split(&gd, d, p, &mut rng)? {
                out.push((q, m));
            }
        }
    }
    out.sort_by(|a, b| a.0.len().cmp(&b.0.len()).then_with(|| a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
    Some((lead, out))
}

fn split_cantor_zassenhaus(
    g: &Poly,
    d: usize,
    p: &BigInt,
    rng: &mut Rng,
) -> Option<Vec<Poly>> {
    let mut parts = Vec::new();
    edf(g, d, p, rng, &mut parts)?;
    Some(parts)
}

// ---------------- Berlekamp ----------------

/// A basis of `{v : M v = 0}` over GF(p) for an `n x n` matrix.
fn nullspace(
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

/// Ascending coefficient vector of length `n`.
fn ascending(
    f: &Poly,
    n: usize,
) -> Vec<BigInt> {
    let mut v: Vec<BigInt> = f.iter().rev().cloned().collect();
    v.resize(n, BigInt::zero());
    v
}

/// Berlekamp's algorithm for a monic squarefree polynomial.
fn berlekamp_split(
    g: &Poly,
    p: &BigInt,
) -> Option<Vec<Poly>> {
    let n = g.len() - 1;
    if n <= 1 {
        return Some(vec![g.clone()]);
    }
    let pu = p.to_usize().filter(|&v| v <= 4096)?;
    let xp = powmod(&x_poly(), p, g, p)?;
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
            let delta = if i == j { BigInt::one() } else { BigInt::zero() };
            system[j][i] = modulo(&(&q[i][j] - delta), p);
        }
    }
    let basis = nullspace(system, p)?;
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

fn split_berlekamp(
    g: &Poly,
    _d: usize,
    p: &BigInt,
    _rng: &mut Rng,
) -> Option<Vec<Poly>> {
    berlekamp_split(g, p)
}

// ---------------- operators ----------------

fn factorisation_value(
    lead: BigInt,
    factors: Vec<(Poly, usize)>,
) -> V {
    let mut items = vec![V::Int(lead)];
    items.extend(factors.into_iter().map(|(g, m)| V::List(vec![poly_value(g), V::uint(m)])));
    V::List(items)
}

fn pairs_value(pairs: Vec<(Poly, usize)>) -> V {
    V::List(pairs.into_iter().map(|(g, m)| V::List(vec![poly_value(g), V::uint(m)])).collect())
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
    Some(pairs_value(squarefree(&monic(&f, &p)?, &p)?))
}

fn gfp_ddf_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let f = monic(&f, &p)?;
    // The input must be squarefree.
    if !is_constant(&gcd(&f, &derivative(&f, &p), &p)?) || is_constant(&f) {
        return None;
    }
    Some(pairs_value(ddf(&f, &p)?))
}

fn gfp_edf_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, d, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let d = usize::try_from(big(cx.graph, *d)?).ok().filter(|&d| d >= 1)?;
    let f = monic(&f, &p)?;
    if is_constant(&f) || (f.len() - 1) % d != 0 || !is_constant(&gcd(&f, &derivative(&f, &p), &p)?) {
        return None;
    }
    let mut parts = Vec::new();
    edf(&f, d, &p, &mut Rng(0x5eed_1234_abcd_9876), &mut parts)?;
    // Verify the premise: every part must be irreducible of degree d.
    if parts.iter().any(|q| q.len() - 1 != d || is_irreducible(q, &p) != Some(true)) {
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
    let (lead, factors) = factor(&f, &p, split_cantor_zassenhaus)?;
    Some(factorisation_value(lead, factors))
}

fn gfp_berlekamp_op(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [f, p] = a else { return None };
    let (f, p) = read_prime_poly(cx, *f, *p)?;
    let (lead, factors) = factor(&f, &p, split_berlekamp)?;
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
            return is_irreducible(&poly, &p).map_or(Outcome::Pass, |b| Outcome::Equal(V::Bool(b).build(cx.graph)));
        }
        let Some((lead, factors)) = factor(&poly, &p, split_cantor_zassenhaus) else {
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
    def(i, "gfp_monic", Arity::Fixed(2), gfp_monic_op)?;
    def(i, "gfp_derivative", Arity::Fixed(2), gfp_derivative_op)?;
    def(i, "gfp_powmod", Arity::Fixed(4), gfp_powmod_op)?;
    def(i, "gfp_invmod", Arity::Fixed(3), gfp_invmod_op)?;
    def(i, "gfp_squarefree", Arity::Fixed(2), gfp_squarefree_op)?;
    def(i, "gfp_ddf", Arity::Fixed(2), gfp_ddf_op)?;
    def(i, "gfp_edf", Arity::Fixed(3), gfp_edf_op)?;
    def(i, "gfp_factor", Arity::Fixed(2), gfp_factor_op)?;
    def(i, "gfp_berlekamp", Arity::Fixed(2), gfp_berlekamp_op)?;
    for (name, request) in [("factor_mod", Request::FactorMod), ("is_irreducible_mod", Request::IsIrreducibleMod)] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("discrete/{name}"), Tier::Reduce, ModRequest { op, request });
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
        List(Vec<Tree>),
    }

    fn parse(text: &str) -> Tree {
        fn go(chars: &[char], at: &mut usize) -> Tree {
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
            Tree::Num(number.parse().unwrap_or_else(|_| panic!("cannot parse `{number}`")))
        }
        go(&text.chars().collect::<Vec<_>>(), &mut 0)
    }

    fn coefficients(tree: &Tree) -> Vec<i64> {
        match tree {
            | Tree::List(items) => items.iter().map(|t| if let Tree::Num(n) = t { *n } else { panic!("nested") }).collect(),
            | Tree::Num(_) => panic!("not a polynomial"),
        }
    }

    fn poly(c: &[i64]) -> String {
        format!("list({})", c.iter().map(ToString::to_string).collect::<Vec<_>>().join(", "))
    }

    fn mul(
        f: &[i64],
        g: &[i64],
        p: i64,
    ) -> Vec<i64> {
        coefficients(&parse(&s(&format!("gfp_mul({}, {}, {p})", poly(f), poly(g)))))
    }

    /// `(lead, [(factor, multiplicity)])` from a `gfp_factor` answer.
    fn factorisation(text: &str) -> (i64, Vec<(Vec<i64>, usize)>) {
        let Tree::List(items) = parse(text) else { panic!("{text}") };
        let Tree::Num(lead) = items[0] else { panic!("{text}") };
        let factors = items[1..]
            .iter()
            .map(|t| {
                let Tree::List(pair) = t else { panic!("{text}") };
                let Tree::Num(m) = pair[1] else { panic!("{text}") };
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
        assert_eq!(s(&format!("gfp_squarefree({}, 2)", poly(&f))), "list(list(list(1, 1, 1), 1), list(list(1, 1), 2))");
        // x^p - x^0 in characteristic p: x^3 + 2 = (x + 2)^3 over GF(3)
        assert_eq!(s("gfp_squarefree(list(1, 0, 0, 2), 3)"), "list(list(list(1, 2), 3))");
        // x^6 + 1 = (x^2 + 1)^3 over GF(3)
        assert_eq!(s("gfp_squarefree(list(1, 0, 0, 0, 0, 0, 1), 3)"), "list(list(list(1, 0, 1), 3))");
        // squarefree input is its own decomposition (made monic)
        assert_eq!(s("gfp_squarefree(list(2, 0, 2), 5)"), "list(list(list(1, 0, 1), 1))");
        assert_eq!(s("gfp_squarefree(list(3), 5)"), "list()");
        // mixed multiplicities 1, 2, 3 and p-th powers
        let f = mul(&mul(&[1, 2], &mul(&[1, 3], &[1, 3], 7), 7), &mul(&[1, 5], &mul(&[1, 5], &[1, 5], 7), 7), 7);
        assert_eq!(
            s(&format!("gfp_squarefree({}, 7)", poly(&f))),
            "list(list(list(1, 2), 1), list(list(1, 3), 2), list(list(1, 5), 3))"
        );
    }

    #[test]
    fn distinct_degree_factorisation() {
        // x^8 - x = x^8 + x over GF(2): all irreducibles of degree 1 and 3
        // x (x+1) (x^3 + x + 1)(x^3 + x^2 + 1)
        assert_eq!(s("gfp_ddf(list(1, 0, 0, 0, 0, 0, 0, 1, 0), 2)"), "list(list(list(1, 1, 0), 1), list(list(1, 1, 1, 1, 1, 1, 1), 3))");
        // x^2 + 1 over GF(5) = (x + 2)(x + 3): one class of degree one
        assert_eq!(s("gfp_ddf(list(1, 0, 1), 5)"), "list(list(list(1, 0, 1), 1))");
        // x^2 + 1 over GF(3) is irreducible
        assert_eq!(s("gfp_ddf(list(1, 0, 1), 3)"), "list(list(list(1, 0, 1), 2))");
        // not squarefree: stays unreduced
        assert_eq!(s("gfp_ddf(list(1, 2, 1), 5)"), "gfp_ddf(list(1, 2, 1), 5)");
    }

    #[test]
    fn equal_degree_factorisation() {
        // the two cubics of GF(2)
        assert_eq!(s("gfp_edf(list(1, 1, 1, 1, 1, 1, 1), 3, 2)"), "list(list(1, 0, 1, 1), list(1, 1, 0, 1))");
        // x^2 + 1 over GF(5): (x + 2)(x + 3)
        assert_eq!(s("gfp_edf(list(1, 0, 1), 1, 5)"), "list(list(1, 2), list(1, 3))");
        // wrong premise: x^2 + 1 over GF(3) has no factor of degree 1
        assert_eq!(s("gfp_edf(list(1, 0, 1), 1, 3)"), "gfp_edf(list(1, 0, 1), 1, 3)");
        // large prime: (x - 3)(x - 5)(x - 11) mod 1000003
        let p = 1_000_003_i64;
        let f = mul(&mul(&[1, p - 3], &[1, p - 5], p), &[1, p - 11], p);
        assert_eq!(
            s(&format!("gfp_edf({}, 1, {p})", poly(&f))),
            format!("list(list(1, {}), list(1, {}), list(1, {}))", p - 11, p - 5, p - 3)
        );
    }

    #[test]
    fn complete_factorisation() {
        assert_eq!(s("gfp_factor(list(1, 0, 0, 0, 1), 2)"), "list(1, list(list(1, 1), 4))");
        assert_eq!(s("gfp_factor(list(3, 0, 3), 5)"), "list(3, list(list(1, 2), 1), list(list(1, 3), 1))");
        assert_eq!(s("gfp_factor(list(4), 7)"), "list(4)");
        assert_eq!(s("gfp_factor(list(), 7)"), "gfp_factor(list(), 7)");
        // composite modulus is rejected
        assert_eq!(s("gfp_factor(list(1, 0, 1), 4)"), "gfp_factor(list(1, 0, 1), 4)");
        // x^4 + 1 over GF(3) = (x^2 + x + 2)(x^2 + 2x + 2)
        assert_eq!(s("gfp_factor(list(1, 0, 0, 0, 1), 3)"), "list(1, list(list(1, 1, 2), 1), list(list(1, 2, 2), 1))");
        assert_eq!(s("gfp_berlekamp(list(1, 0, 0, 0, 1), 3)"), "list(1, list(list(1, 1, 2), 1), list(list(1, 2, 2), 1))");
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
                let mut f: Vec<i64> = (0..=degree).map(|_| i64::try_from(rng.next(u64::try_from(p).unwrap_or(2))).unwrap_or(0)).collect();
                if f[0] == 0 {
                    f[0] = 1;
                }
                let want = coefficients(&parse(&s(&format!("gfp_norm({}, {p})", poly(&f)))));
                let (lead, factors) = factorisation(&s(&format!("gfp_factor({}, {p})", poly(&f))));
                assert_eq!(rebuild(lead, &factors, p), want, "{f:?} mod {p}");
                for (g, _) in &factors {
                    assert!(brute_irreducible(g, p), "{g:?} reducible mod {p} (from {f:?})");
                    assert_eq!(s(&format!("gfp_isirreducible({}, {p})", poly(g))), "true");
                }
                let (lead2, factors2) = factorisation(&s(&format!("gfp_berlekamp({}, {p})", poly(&f))));
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
                    assert_eq!(s(&format!("gfp_isirreducible({}, {p})", poly(&f))), brute_irreducible(&f, p).to_string(), "{f:?} mod {p}");
                }
            }
        }
        assert_eq!(s("gfp_isirreducible(list(3), 5)"), "false");
        assert_eq!(s("gfp_isirreducible(list(2, 1), 5)"), "true");
        // a degree-12 irreducible over GF(2): x^12 + x^6 + x^4 + x + 1
        assert_eq!(s("gfp_isirreducible(list(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1), 2)"), "true");
    }

    #[test]
    fn powmod_invmod_and_derivative() {
        // x^10 mod (x^2 + 1) over GF(5): x^2 = -1, x^10 = -1 = 4
        assert_eq!(s("gfp_powmod(list(1, 0), 10, list(1, 0, 1), 5)"), "list(4)");
        assert_eq!(s("gfp_powmod(list(1, 1), 0, list(1, 0, 1), 5)"), "list(1)");
        // Frobenius: (x + 1)^7 = x^7 + 1 over GF(7)
        assert_eq!(s("gfp_powmod(list(1, 1), 7, list(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), 7)"), "list(1, 0, 0, 0, 0, 0, 0, 1)");
        assert_eq!(s("gfp_powmod(list(1, 1), -1, list(1, 0, 1), 5)"), "gfp_powmod(list(1, 1), -1, list(1, 0, 1), 5)");
        // x (x + 1) = x^2 + x = 1 + x - 1... check f * f^-1 = 1 mod m
        let m = "list(1, 0, 1, 1)"; // x^3 + x + 1 over GF(2): a field
        for f in ["list(1, 0)", "list(1, 1)", "list(1, 0, 1)", "list(1, 1, 1)"] {
            let inv = s(&format!("gfp_invmod({f}, {m}, 2)"));
            let product = s(&format!("gfp_mul({f}, {inv}, 2)"));
            assert_eq!(s(&format!("gfp_divmod({product}, {m}, 2)")).rsplit_once(", ").map(|(_, r)| r.to_string()), Some("list(1))".to_string()), "{f}");
        }
        assert_eq!(s("gfp_invmod(list(1, 1), list(1, 0, 1), 2)"), "gfp_invmod(list(1, 1), list(1, 0, 1), 2)");
        assert_eq!(s("gfp_derivative(list(1, 0, 0, 5, 3), 7)"), "list(4, 0, 0, 5)");
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
        assert_eq!(s("factor_mod(list(1, 0, 0, 0, 1), x, 2)"), "list(1, list(list(1, 1), 4))");
        // irreducibility
        assert_eq!(s("is_irreducible_mod(x^2 + x + 1, x, 2)"), "true");
        assert_eq!(s("is_irreducible_mod(x^2 + 1, x, 2)"), "false");
        assert_eq!(s("is_irreducible_mod(list(1, 1, 1), x, 2)"), "true");
        // other variables or non-prime moduli stay unreduced
        for src in ["factor_mod(x^2 + y, x, 5)", "factor_mod(x^2 + 1, x, 6)"] {
            let (text, reduced) = crate::rules::testing::reduce_with(&[crate::rules::discrete::discrete()], src, &[]);
            assert_eq!((text.as_str(), reduced), (src, false));
        }
    }
}
