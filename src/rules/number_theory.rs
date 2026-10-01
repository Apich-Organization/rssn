//! Number theory: divisibility, modular arithmetic, primality and
//! factorisation of exact integers.
//!
//! Every operator is reduced by an exact kernel that fires only when all of
//! its arguments are literal integers; anything else stays symbolic. Results
//! are computed with arbitrary-precision integers, so nothing overflows.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

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
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

use super::arith::arith;
use super::poly::best;
use super::poly::repr::Gens;
use super::poly::repr::Limits;
use super::poly::repr::Poly;
use super::poly::repr::from_term;

/// The number-theory rule set.
#[must_use]
pub fn number_theory() -> RuleSet {
    RuleSet::new("number_theory", install).needs(arith())
}

/// What an exact kernel computed.
pub(crate) enum Value {
    /// An integer.
    Int(BigInt),
    /// An exact rational.
    Rat(BigRational),
    /// A truth value.
    Bool(bool),
    /// `list(a, b, ...)` of integers.
    Flat(Vec<BigInt>),
    /// `list(list(a, b), list(c, d), ...)` of integers.
    Nested(Vec<Vec<BigInt>>),
}

impl Value {
    /// Lists are larger than the request that produced them, so extraction
    /// would otherwise prefer the request.
    const fn is_list(&self) -> bool {
        matches!(self, Self::Flat(_) | Self::Nested(_))
    }

    fn build(
        self,
        graph: &mut Graph,
    ) -> NodeId {
        match self {
            | Self::Int(v) => graph.num(Number::Int(v)),
            | Self::Rat(v) => graph.num(Number::rat(v)),
            | Self::Bool(b) => graph.lit(Payload::Bool(b)),
            | Self::Flat(items) => list_node(graph, items),
            | Self::Nested(rows) => {
                let mut nodes = Vec::with_capacity(rows.len());
                for row in rows {
                    nodes.push(list_node(graph, row));
                }
                graph.node(core::LIST, &nodes)
            },
        }
    }
}

fn list_node(
    graph: &mut Graph,
    items: Vec<BigInt>,
) -> NodeId {
    let mut nodes = Vec::with_capacity(items.len());
    for v in items {
        nodes.push(graph.num(Number::Int(v)));
    }
    graph.node(core::LIST, &nodes)
}

/// An exact computation on integer arguments; `None` when it declines.
pub(crate) type Run = fn(&[BigInt]) -> Option<Value>;

/// An exact computation on lists of integers; `None` when it declines.
pub(crate) type ListRun = fn(&[Vec<BigInt>]) -> Option<Value>;

fn int_of(
    graph: &Graph,
    node: NodeId,
) -> Option<BigInt> {
    match graph.number_of(node)? {
        | Number::Int(v) => Some(v.clone()),
        | _ => None,
    }
}

fn list_of(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<BigInt>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    graph
        .children(node)
        .iter()
        .map(|&c| int_of(graph, c))
        .collect()
}

/// Reduces `op` when every argument is a literal integer.
pub(crate) struct Exact {
    op: OpId,
    run: Run,
}

impl Kernel for Exact {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args: Option<Vec<BigInt>> = cx
            .graph
            .children(node)
            .iter()
            .map(|&c| int_of(cx.graph, c))
            .collect();
        let Some(value) = args.and_then(|a| (self.run)(&a)) else {
            return Outcome::Pass;
        };
        let pinned = value.is_list();
        let node = value.build(cx.graph);
        if pinned {
            Outcome::Pinned(node)
        } else {
            Outcome::Equal(node)
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Reduces `op` when every argument is a `list` of literal integers.
pub(crate) struct ExactLists {
    op: OpId,
    run: ListRun,
}

impl Kernel for ExactLists {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args: Option<Vec<Vec<BigInt>>> = cx
            .graph
            .children(node)
            .iter()
            .map(|&c| list_of(cx.graph, c))
            .collect();
        let Some(value) = args.and_then(|a| (self.run)(&a)) else {
            return Outcome::Pass;
        };
        let pinned = value.is_list();
        let node = value.build(cx.graph);
        if pinned {
            Outcome::Pinned(node)
        } else {
            Outcome::Equal(node)
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Registers an operator together with its exact integer kernel.
pub(crate) fn exact(
    i: &mut Installer<'_>,
    set: &str,
    desc: OpDescriptor,
    run: Run,
) -> Result<(), RuleError> {
    let name = format!("{set}/{}", desc.name);
    let op = i.op(desc)?;
    i.kernel(&name, Tier::Normalize, Exact { op, run });
    Ok(())
}

/// Registers an operator together with its exact kernel on integer lists.
pub(crate) fn exact_lists(
    i: &mut Installer<'_>,
    set: &str,
    desc: OpDescriptor,
    run: ListRun,
) -> Result<(), RuleError> {
    let name = format!("{set}/{}", desc.name);
    let op = i.op(desc)?;
    i.kernel(&name, Tier::Normalize, ExactLists { op, run });
    Ok(())
}

/// Rounds a float argument to an integer for the numeric semantics.
fn round_arg(x: f64) -> Option<i64> {
    x.round().to_i64().filter(|v| v.unsigned_abs() < 1 << 62)
}

/// `gcd` on floats: arguments are rounded, since it is defined on integers.
fn gcd_eval(args: &[f64]) -> f64 {
    let (Some(&a), Some(&b)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    let (Some(a), Some(b)) = (round_arg(a), round_arg(b)) else {
        return f64::NAN;
    };
    let (mut x, mut y) = (a.unsigned_abs(), b.unsigned_abs());
    while y != 0 {
        (x, y) = (y, x % y);
    }
    x.to_f64().unwrap_or(f64::NAN)
}

fn lcm_eval(args: &[f64]) -> f64 {
    let g = gcd_eval(args);
    let (Some(&a), Some(&b)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if g == 0.0 {
        return 0.0;
    }
    (a.round() / g * b.round()).abs()
}

fn mod_eval(args: &[f64]) -> f64 {
    let (Some(&a), Some(&m)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    let (Some(a), Some(m)) = (round_arg(a), round_arg(m)) else {
        return f64::NAN;
    };
    if m == 0 {
        return f64::NAN;
    }
    a.rem_euclid(m.abs()).to_f64().unwrap_or(f64::NAN)
}

/// Greatest common divisor, non-negative.
pub(crate) fn gcd(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    let (mut x, mut y) = (a.abs(), b.abs());
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

/// Remainder in `[0, |m|)`; `None` for `m == 0`.
fn modulo(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigInt> {
    if m.is_zero() {
        return None;
    }
    let m = m.abs();
    let r = a % &m;
    Some(if r.is_negative() { r + m } else { r })
}

/// `(g, x, y)` with `a*x + b*y = g = gcd(a, b) >= 0`.
fn egcd(
    a: &BigInt,
    b: &BigInt,
) -> (BigInt, BigInt, BigInt) {
    let (mut r0, mut r1) = (a.abs(), b.abs());
    let (mut s0, mut s1) = (BigInt::one(), BigInt::zero());
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    while !r1.is_zero() {
        let q = &r0 / &r1;
        let r2 = &r0 - &q * &r1;
        let s2 = &s0 - &q * &s1;
        let t2 = &t0 - &q * &t1;
        (r0, r1) = (r1, r2);
        (s0, s1) = (s1, s2);
        (t0, t1) = (t1, t2);
    }
    if a.is_negative() {
        s0 = -s0;
    }
    if b.is_negative() {
        t0 = -t0;
    }
    (r0, s0, t0)
}

/// The inverse of `a` modulo `m > 0`, if `gcd(a, m) = 1`.
fn invmod(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigInt> {
    if !m.is_positive() {
        return None;
    }
    let (g, x, _) = egcd(&modulo(a, m)?, m);
    if g.is_one() {
        modulo(&x, m)
    } else if m.is_one() {
        Some(BigInt::zero())
    } else {
        None
    }
}

const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

fn mulmod(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    let wide = u128::from(a) * u128::from(b) % u128::from(m);
    u64::try_from(wide).unwrap_or(0)
}

fn powmod_u64(
    mut base: u64,
    mut exp: u64,
    m: u64,
) -> u64 {
    let mut acc = 1 % m;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mulmod(acc, base, m);
        }
        base = mulmod(base, base, m);
        exp >>= 1;
    }
    acc
}

/// Deterministic Miller-Rabin: the first 12 primes are witnesses for every
/// `n < 3.3e24`, so all of `u64`.
fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'witness: for a in BASES {
        let mut x = powmod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

fn is_prime(n: &BigInt) -> bool {
    if let Some(small) = n.to_u64() {
        return is_prime_u64(small);
    }
    if n.is_negative() {
        return false;
    }
    let one = BigInt::one();
    for p in BASES {
        if (n % p).is_zero() {
            return false;
        }
    }
    let n1 = n - &one;
    let s = n1.trailing_zeros().unwrap_or(0);
    let d = &n1 >> s;
    'witness: for a in BASES {
        let mut x = BigInt::from(a).modpow(&d, n);
        if x == one || x == n1 {
            continue;
        }
        for _ in 1..s {
            x = &x * &x % n;
            if x == n1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

const fn gcd_u64(
    mut a: u64,
    mut b: u64,
) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A non-trivial factor of the odd composite `n` by Pollard's rho.
fn rho(n: u64) -> Option<u64> {
    let step = |x: u64, c: u64| {
        let wide = (u128::from(x) * u128::from(x) + u128::from(c)) % u128::from(n);
        u64::try_from(wide).unwrap_or(0)
    };
    for c in 1..64 {
        let (mut x, mut y, mut d) = (2_u64, 2_u64, 1_u64);
        while d == 1 {
            x = step(x, c);
            y = step(step(y, c), c);
            d = gcd_u64(x.abs_diff(y), n);
        }
        if d != n {
            return Some(d);
        }
    }
    None
}

/// Appends the prime factors of `n` (with multiplicity) to `out`.
fn split(
    n: u64,
    out: &mut Vec<u64>,
) -> bool {
    if n == 1 {
        return true;
    }
    if is_prime_u64(n) {
        out.push(n);
        return true;
    }
    // Rho cycles too early on prime squares.
    let root = n.isqrt();
    if root * root == n {
        return split(root, out) && split(root, out);
    }
    rho(n).is_some_and(|d| split(d, out) && split(n / d, out))
}

/// Prime factorisation as ascending `(prime, exponent)` pairs.
pub(crate) fn factor(n: &BigInt) -> Option<Vec<(u64, u32)>> {
    let mut n = n.to_u64().filter(|&v| v >= 1)?;
    let mut primes = Vec::new();
    for d in 2..=1000_u64 {
        if d * d > n {
            break;
        }
        while n % d == 0 {
            primes.push(d);
            n /= d;
        }
    }
    if !split(n, &mut primes) {
        return None;
    }
    primes.sort_unstable();
    let mut out: Vec<(u64, u32)> = Vec::new();
    for p in primes {
        match out.last_mut() {
            | Some((q, e)) if *q == p => *e += 1,
            | _ => out.push((p, 1)),
        }
    }
    Some(out)
}

fn next_prime(n: &BigInt) -> Option<BigInt> {
    let two = BigInt::from(2);
    if *n < two {
        return Some(two);
    }
    let mut candidate: BigInt = n + 1;
    if !candidate.bit(0) {
        candidate += 1;
    }
    // Prime gaps are tiny; the bound only guards against absurd inputs.
    for _ in 0..1_000_000 {
        if is_prime(&candidate) {
            return Some(candidate);
        }
        candidate += 2;
    }
    None
}

fn jacobi(
    a: &BigInt,
    n: &BigInt,
) -> Option<i8> {
    if !n.is_positive() || !n.bit(0) {
        return None;
    }
    let (mut a, mut n) = (modulo(a, n)?, n.clone());
    let mut sign = 1_i8;
    let low = |x: &BigInt, m: u32| (x % m).to_u32().unwrap_or(0);
    while !a.is_zero() {
        while !a.bit(0) {
            a >>= 1;
            if matches!(low(&n, 8), 3 | 5) {
                sign = -sign;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if low(&a, 4) == 3 && low(&n, 4) == 3 {
            sign = -sign;
        }
        a = modulo(&a, &n)?;
    }
    Some(if n.is_one() { sign } else { 0 })
}

fn crt(args: &[Vec<BigInt>]) -> Option<Value> {
    let [residues, moduli] = args else {
        return None;
    };
    if residues.len() != moduli.len() || residues.is_empty() {
        return None;
    }
    let (mut x, mut product) = (BigInt::zero(), BigInt::one());
    for (r, m) in residues.iter().zip(moduli) {
        // Solve x + product*t = r (mod m); needs product invertible mod m,
        // i.e. the moduli are pairwise coprime.
        let inverse = invmod(&product, m)?;
        let t = modulo(&(modulo(&(r - &x), m)? * inverse), m)?;
        x += &product * t;
        product *= m;
    }
    modulo(&x, &product).map(Value::Int)
}

/// A kernel that works on the argument nodes directly, for operators whose
/// arguments are not plain integers or lists of integers.
pub(crate) type Build = fn(&mut Cx<'_>, &[NodeId]) -> Outcome;

/// Reduces `op` with a [`Build`] function.
pub(crate) struct Procedure {
    op: OpId,
    run: Build,
}

impl Kernel for Procedure {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        (self.run)(cx, &args)
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Registers an operator together with a [`Build`] kernel of the given tier.
pub(crate) fn procedure(
    i: &mut Installer<'_>,
    set: &str,
    desc: OpDescriptor,
    tier: Tier,
    run: Build,
) -> Result<OpId, RuleError> {
    let name = format!("{set}/{}", desc.name);
    let op = i.op(desc)?;
    i.kernel(&name, tier, Procedure { op, run });
    Ok(op)
}

/// `lhs - rhs` for an equation, the term itself otherwise.
pub(crate) fn residual(
    graph: &mut Graph,
    equation: NodeId,
) -> NodeId {
    match *graph.children(equation) {
        | [lhs, rhs] if graph.op(equation) == core::EQ => {
            let minus_one = graph.int(-1);
            let negated = graph.node(core::MUL, &[minus_one, rhs]);
            graph.node(core::ADD, &[lhs, negated])
        },
        | _ => equation,
    }
}

/// Longest period of a continued fraction that is searched for.
const MAX_PERIOD: usize = 200_000;
/// Most partial quotients or convergents produced by one request.
const MAX_TERMS: u64 = 10_000;
/// Largest limit of the prime sieve (9592 primes).
const MAX_SIEVE: u64 = 100_000;

/// `sqrt(n) = [a0; period, period, ...]` for a positive non-square `n`.
fn sqrt_expansion(n: &BigInt) -> Option<(BigInt, Vec<BigInt>)> {
    if !n.is_positive() {
        return None;
    }
    let a0 = n.sqrt();
    if &a0 * &a0 == *n {
        return None;
    }
    let two_a0 = &a0 * 2;
    let (mut m, mut d, mut a) = (BigInt::zero(), BigInt::one(), a0.clone());
    let mut period = Vec::new();
    while period.len() < MAX_PERIOD {
        m = &d * &a - &m;
        d = (n - &m * &m) / &d;
        a = (&a0 + &m) / &d;
        period.push(a.clone());
        if a == two_a0 {
            return Some((a0, period));
        }
    }
    None
}

/// The partial quotients of a rational, by the Euclidean algorithm.
fn rational_expansion(r: &BigRational) -> Vec<BigInt> {
    let mut out = Vec::new();
    let mut x = r.clone();
    loop {
        let a = x.floor().to_integer();
        let fraction = &x - BigRational::from_integer(a.clone());
        out.push(a);
        if fraction.is_zero() {
            return out;
        }
        x = fraction.recip();
    }
}

/// The convergents `h/k` of the continued fraction with these quotients.
fn convergent_pairs(quotients: &[BigInt]) -> Vec<(BigInt, BigInt)> {
    let (mut h1, mut h2) = (BigInt::one(), BigInt::zero());
    let (mut k1, mut k2) = (BigInt::zero(), BigInt::one());
    let mut out = Vec::with_capacity(quotients.len());
    for a in quotients {
        let h = a * &h1 + &h2;
        let k = a * &k1 + &k2;
        (h2, h1) = (h1, h.clone());
        (k2, k1) = (k1, k.clone());
        out.push((h, k));
    }
    out
}

/// Outcome of a Pell equation.
enum Pell {
    Solved(BigInt, BigInt),
    Unsolvable,
}

/// The smallest positive solution of `x^2 - n y^2 = 1` (or `= -1` when
/// `negative`). `None` when `n` is not a positive non-square (or the
/// search is too long), `Some(Pell::Unsolvable)` when there is none.
fn pell_solution(
    n: &BigInt,
    negative: bool,
) -> Option<Pell> {
    let (a0, period) = sqrt_expansion(n)?;
    let m = period.len();
    let index = if negative {
        if m % 2 == 0 {
            return Some(Pell::Unsolvable);
        }
        m - 1
    } else if m % 2 == 0 {
        m - 1
    } else {
        2 * m - 1
    };
    let mut quotients = vec![a0];
    while quotients.len() <= index {
        quotients.push(period.get((quotients.len() - 1) % m)?.clone());
    }
    convergent_pairs(&quotients)
        .pop()
        .map(|(x, y)| Pell::Solved(x, y))
}

/// A continued fraction: finite, or `[a0; period...]` of a square root.
enum Expansion {
    Finite(Vec<BigInt>),
    Surd(BigInt, Vec<BigInt>),
}

impl Expansion {
    /// The first `count` partial quotients.
    fn terms(
        &self,
        count: usize,
    ) -> Vec<BigInt> {
        match self {
            | Self::Finite(q) => q.iter().take(count).cloned().collect(),
            | Self::Surd(a0, period) => std::iter::once(a0.clone())
                .chain(period.iter().cycle().cloned())
                .take(count)
                .collect(),
        }
    }
}

/// `n` if `term` is `sqrt(n)` or `n^(1/2)` for a positive integer `n`.
fn surd_radicand(
    graph: &Graph,
    term: NodeId,
) -> Option<BigInt> {
    let half = Number::fraction(1, 2)?;
    let radicand = match (graph.op(term), graph.children(term)) {
        | (core::POW, [base, exp]) if graph.number_of(*exp).is_some_and(|e| *e == half) => *base,
        | (op, [arg]) if graph.ops().lookup("sqrt") == Some(op) => *arg,
        | _ => return None,
    };
    int_of(graph, radicand).filter(Signed::is_positive)
}

/// The continued fraction of a rational, of `sqrt(n)` and — when
/// `allow_list` — of a `list` of partial quotients.
fn expansion(
    graph: &mut Graph,
    node: NodeId,
    allow_list: bool,
) -> Option<Expansion> {
    let term = best(graph, node)?;
    if let Some(r) = graph.number_of(term).and_then(Number::to_rational) {
        return Some(Expansion::Finite(rational_expansion(&r)));
    }
    if let Some(n) = surd_radicand(graph, term) {
        let root = n.sqrt();
        if &root * &root == n {
            return Some(Expansion::Finite(vec![root]));
        }
        return sqrt_expansion(&n).map(|(a0, period)| Expansion::Surd(a0, period));
    }
    if allow_list {
        return list_of(graph, term).map(Expansion::Finite);
    }
    None
}

fn cfrac(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[x] = args else {
        return Outcome::Pass;
    };
    let graph = &mut *cx.graph;
    match expansion(graph, x, false) {
        | Some(Expansion::Finite(quotients)) => Outcome::Pinned(list_node(graph, quotients)),
        | Some(Expansion::Surd(a0, period)) => {
            let first = graph.num(Number::Int(a0));
            let period = list_node(graph, period);
            Outcome::Pinned(graph.node(core::LIST, &[first, period]))
        },
        | None => Outcome::Pass,
    }
}

fn convergents(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[x, k] = args else {
        return Outcome::Pass;
    };
    let graph = &mut *cx.graph;
    let Some(k) = int_of(graph, k)
        .and_then(|k| k.to_u64())
        .filter(|&k| k < MAX_TERMS)
    else {
        return Outcome::Pass;
    };
    let Some(expansion) = expansion(graph, x, true) else {
        return Outcome::Pass;
    };
    let count = usize::try_from(k + 1).unwrap_or(0);
    let mut nodes = Vec::new();
    for (h, k) in convergent_pairs(&expansion.terms(count)) {
        nodes.push(graph.num(Number::rat(BigRational::new(h, k))));
    }
    Outcome::Pinned(graph.node(core::LIST, &nodes))
}

/// A quadratic polynomial without mixed terms, by variable.
struct Shape {
    constant: BigRational,
    linear: Vec<BigRational>,
    square: Vec<BigRational>,
}

/// The coefficients of `poly` if it only has constant, linear and square
/// terms in the generators `vars`, all with exact coefficients.
fn shape(
    poly: &Poly,
    vars: &[u32],
) -> Option<Shape> {
    let mut out = Shape {
        constant: BigRational::zero(),
        linear: vec![BigRational::zero(); vars.len()],
        square: vec![BigRational::zero(); vars.len()],
    };
    for (mono, coeff) in poly.terms() {
        let c = coeff.to_rational()?;
        match mono.as_slice() {
            | [] => out.constant = c,
            | &[(g, e)] => {
                let index = vars.iter().position(|&v| v == g)?;
                let slot = match e {
                    | 1 => out.linear.get_mut(index)?,
                    | 2 => out.square.get_mut(index)?,
                    | _ => return None,
                };
                *slot = c;
            },
            | _ => return None,
        }
    }
    Some(out)
}

/// The values scaled by the least common multiple of their denominators.
fn integers(values: &[&BigRational]) -> Vec<BigInt> {
    let scale = values.iter().fold(BigInt::one(), |acc, v| {
        &acc / gcd(&acc, v.denom()) * v.denom()
    });
    values
        .iter()
        .map(|v| (*v * BigRational::from_integer(scale.clone())).to_integer())
        .collect()
}

/// A name for the free parameter that none of `vars` uses.
fn parameter(
    graph: &mut Graph,
    vars: &[NodeId],
    candidates: &[&str],
) -> Option<Vec<NodeId>> {
    let taken: Vec<String> = vars.iter().map(|&v| graph.display(v)).collect();
    if candidates.iter().any(|c| taken.iter().any(|t| t == c)) {
        return None;
    }
    Some(candidates.iter().map(|c| graph.sym(c)).collect())
}

/// `a x + b y = c`: all solutions with one integer parameter `t`.
fn linear_diophantine(
    graph: &mut Graph,
    shape: &Shape,
    vars: &[NodeId],
) -> Option<NodeId> {
    let (a, b) = (shape.linear.first()?, shape.linear.get(1)?);
    let c = -&shape.constant;
    let [a, b, c] = <[BigInt; 3]>::try_from(integers(&[a, b, &c])).ok()?;
    if a.is_zero() && b.is_zero() {
        // `0 = c`: no solution, or every pair (which a list cannot say).
        return (!c.is_zero()).then(|| graph.node(core::LIST, &[]));
    }
    let (g, x0, y0) = egcd(&a, &b);
    if !(&c % &g).is_zero() {
        return Some(graph.node(core::LIST, &[]));
    }
    let f = &c / &g;
    let t = *parameter(graph, vars, &["t"])?.first()?;
    let mut solution = Vec::new();
    for (particular, step) in [(x0 * &f, &b / &g), (y0 * &f, -(&a / &g))] {
        let (particular, step) = (
            graph.num(Number::Int(particular)),
            graph.num(Number::Int(step)),
        );
        let moved = graph.node(core::MUL, &[step, t]);
        solution.push(graph.node(core::ADD, &[particular, moved]));
    }
    Some(graph.node(core::LIST, &solution))
}

/// `u^2 - n v^2 = ±1` in the two variables, in either order.
fn pell_diophantine(
    graph: &mut Graph,
    shape: &Shape,
) -> Option<NodeId> {
    if shape.linear.iter().any(|c| !c.is_zero()) {
        return None;
    }
    let (s0, s1) = (shape.square.first()?, shape.square.get(1)?);
    // The variable with the positive square coefficient plays `u`.
    let (u_first, s_u, s_v) = if s0.is_positive() && s1.is_negative() {
        (true, s0, s1)
    } else if s1.is_positive() && s0.is_negative() {
        (false, s1, s0)
    } else {
        return None;
    };
    let n = -s_v / s_u;
    let rhs = -&shape.constant / s_u;
    let negative = if rhs.is_one() {
        false
    } else if rhs == -BigRational::one() {
        true
    } else {
        return None;
    };
    if !n.is_integer() {
        return None;
    }
    let (u, v) = match pell_solution(&n.to_integer(), negative)? {
        | Pell::Solved(x, y) => (x, y),
        | Pell::Unsolvable => return Some(graph.node(core::LIST, &[])),
    };
    let pair = if u_first {
        vec![u, v]
    } else {
        vec![v, u]
    };
    Some(list_node(graph, pair))
}

/// `x^2 + y^2 = z^2` (in any arrangement): the Pythagorean triples, with
/// integer parameters `m`, `n`, `k`.
fn pythagorean_diophantine(
    graph: &mut Graph,
    shape: &Shape,
    vars: &[NodeId],
) -> Option<NodeId> {
    let minus_one = -BigRational::one();
    let negatives = shape.square.iter().filter(|c| **c == minus_one).count();
    let positives = shape.square.iter().filter(|c| c.is_one()).count();
    if negatives != 1
        || positives != 2
        || !shape.constant.is_zero()
        || shape.linear.iter().any(|c| !c.is_zero())
    {
        return None;
    }
    let names = parameter(graph, vars, &["m", "n", "k"])?;
    let (&m, &n, &k) = (names.first()?, names.get(1)?, names.get(2)?);
    let two = graph.int(2);
    let (m2, n2) = (
        graph.node(core::POW, &[m, two]),
        graph.node(core::POW, &[n, two]),
    );
    let minus_one_node = graph.int(-1);
    let minus_n2 = graph.node(core::MUL, &[minus_one_node, n2]);
    let leg = {
        let difference = graph.node(core::ADD, &[m2, minus_n2]);
        graph.node(core::MUL, &[k, difference])
    };
    let other = graph.node(core::MUL, &[two, k, m, n]);
    let hypotenuse = {
        let sum = graph.node(core::ADD, &[m2, n2]);
        graph.node(core::MUL, &[k, sum])
    };
    let mut legs = [leg, other].into_iter();
    let solution: Vec<NodeId> = shape
        .square
        .iter()
        .map(|c| {
            if *c == minus_one {
                hypotenuse
            } else {
                legs.next().unwrap_or(leg)
            }
        })
        .collect();
    Some(graph.node(core::LIST, &solution))
}

fn diophantine(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let graph = &mut *cx.graph;
    let &[equation, unknowns] = args else {
        return Outcome::Pass;
    };
    if graph.op(unknowns) != core::LIST {
        return Outcome::Pass;
    }
    let vars = graph.children(unknowns).to_vec();
    if vars.iter().any(|&v| graph.as_symbol(v).is_none()) {
        return Outcome::Pass;
    }
    let expr = residual(graph, equation);
    let Some(term) = best(graph, expr) else {
        return Outcome::Pass;
    };
    let mut gens = Gens::default();
    let indices: Vec<u32> = vars.iter().map(|&v| gens.index(graph, v)).collect();
    let Some(poly) = from_term(graph, &mut gens, term, Limits::default()) else {
        return Outcome::Pass;
    };
    let Some(shape) = shape(&poly, &indices) else {
        return Outcome::Pass;
    };
    let quadratic = shape.square.iter().any(|c| !c.is_zero());
    let solution = match (vars.len(), quadratic) {
        | (2, false) => linear_diophantine(graph, &shape, &vars),
        | (2, true) => pell_diophantine(graph, &shape),
        | (3, true) => pythagorean_diophantine(graph, &shape, &vars),
        | _ => None,
    };
    solution.map_or(Outcome::Pass, Outcome::Equal)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let commutative = OpFlags::COMMUTATIVE;
    let set = "number_theory";

    exact(
        i,
        set,
        binary("gcd").flags(commutative).eval(gcd_eval),
        |a| match a {
            | [x, y] => Some(Value::Int(gcd(x, y))),
            | _ => None,
        },
    )?;
    exact(
        i,
        set,
        binary("lcm").flags(commutative).eval(lcm_eval),
        |a| match a {
            | [x, y] => {
                let g = gcd(x, y);
                if g.is_zero() {
                    return Some(Value::Int(g));
                }
                Some(Value::Int((x / g * y).abs()))
            },
            | _ => None,
        },
    )?;
    exact(i, set, binary("mod").eval(mod_eval), |a| match a {
        | [x, m] => modulo(x, m).map(Value::Int),
        | _ => None,
    })?;
    exact(
        i,
        set,
        OpDescriptor::new("powmod", Arity::Fixed(3)),
        |a| match a {
            | [b, e, m] if !e.is_negative() && m.is_positive() => {
                Some(Value::Int(modulo(b, m)?.modpow(e, m)))
            },
            | _ => None,
        },
    )?;
    exact(i, set, binary("invmod"), |a| match a {
        | [x, m] => invmod(x, m).map(Value::Int),
        | _ => None,
    })?;
    exact(i, set, unary("isprime"), |a| match a {
        | [n] => Some(Value::Bool(is_prime(n))),
        | _ => None,
    })?;
    exact(i, set, unary("nextprime"), |a| match a {
        | [n] => next_prime(n).map(Value::Int),
        | _ => None,
    })?;
    exact(i, set, unary("totient"), |a| {
        let [n] = a else {
            return None;
        };
        let mut phi = BigInt::one();
        for (p, e) in factor(n)? {
            phi *= BigInt::from(p - 1) * BigInt::from(p).pow(e - 1);
        }
        Some(Value::Int(phi))
    })?;
    exact(i, set, unary("divisor_count"), |a| {
        let [n] = a else {
            return None;
        };
        let count = factor(n)?
            .iter()
            .fold(BigInt::one(), |acc, &(_, e)| acc * (e + 1));
        Some(Value::Int(count))
    })?;
    exact(i, set, unary("divisor_sum"), |a| {
        let [n] = a else {
            return None;
        };
        let mut sigma = BigInt::one();
        for (p, e) in factor(n)? {
            // 1 + p + ... + p^e
            let p = BigInt::from(p);
            sigma *= (0..=e).fold(BigInt::zero(), |acc, k| acc + p.pow(k));
        }
        Some(Value::Int(sigma))
    })?;
    exact(i, set, unary("factorint"), |a| {
        let [n] = a else {
            return None;
        };
        let rows = factor(n)?
            .into_iter()
            .map(|(p, e)| vec![BigInt::from(p), BigInt::from(e)])
            .collect();
        Some(Value::Nested(rows))
    })?;
    exact(i, set, binary("jacobi"), |a| match a {
        | [x, n] => jacobi(x, n).map(|j| Value::Int(BigInt::from(j))),
        | _ => None,
    })?;
    exact(i, set, binary("egcd"), |a| match a {
        | [x, y] => {
            let (g, s, t) = egcd(x, y);
            Some(Value::Flat(vec![g, s, t]))
        },
        | _ => None,
    })?;
    exact_lists(i, set, binary("crt"), crt)?;
    exact(i, set, unary("pell"), |a| match a {
        | [n] => match pell_solution(n, false)? {
            | Pell::Solved(x, y) => Some(Value::Flat(vec![x, y])),
            | Pell::Unsolvable => None,
        },
        | _ => None,
    })?;
    exact(i, set, unary("primes"), |a| match a {
        | [n] => {
            let limit = usize::try_from(n.to_u64().filter(|&v| v <= MAX_SIEVE)?).ok()?;
            let primes = crate::kernels::number_theory::primes_sieve(limit);
            Some(Value::Flat(primes.into_iter().map(BigInt::from).collect()))
        },
        | _ => None,
    })?;
    procedure(i, set, unary("cfrac"), Tier::Normalize, cfrac)?;
    procedure(i, set, binary("convergents"), Tier::Normalize, convergents)?;
    procedure(
        i,
        set,
        binary("diophantine").flags(OpFlags::HEAVY).cost(100),
        Tier::Reduce,
        diophantine,
    )?;

    // Values on floats are rounded integers, so only identities that hold
    // for every rounding are stated for symbolic arguments. The others are
    // limited to literal integers, which is what they can be checked on.
    i.rewrites(
        Tier::Normalize,
        &[
            "number_theory/gcd-0: gcd(?a, 0) => ?a if integer(?a), nonnegative(?a)",
            "number_theory/gcd-1: gcd(?a, 1) => 1",
            "number_theory/gcd-gcd: gcd(gcd(?a, ?b), ?b) => gcd(?a, ?b)",
            "number_theory/lcm-0: lcm(?a, 0) => 0",
            "number_theory/lcm-1: lcm(?a, 1) => ?a if integer(?a), nonnegative(?a)",
            "number_theory/lcm-lcm: lcm(lcm(?a, ?b), ?b) => lcm(?a, ?b)",
            "number_theory/gcd-lcm: gcd(?a, ?b) * lcm(?a, ?b) => ?a * ?b \
             if integer(?a), integer(?b), positive(?a), positive(?b)",
            "number_theory/mod-1: mod(?a, 1) => 0",
            "number_theory/mod-mod: mod(mod(?a, ?m), ?m) => mod(?a, ?m)",
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Facts;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[number_theory()], src)
    }

    #[test]
    fn gcd_and_lcm() {
        assert_eq!(s("gcd(12, 18)"), "6");
        assert_eq!(s("gcd(-12, 18)"), "6");
        assert_eq!(s("gcd(0, 5)"), "5");
        assert_eq!(s("gcd(0, 0)"), "0");
        assert_eq!(s("gcd(17, 5)"), "1");
        assert_eq!(s("gcd(2^100, 6^50)"), "1125899906842624");
        assert_eq!(s("lcm(4, 6)"), "12");
        assert_eq!(s("lcm(-4, 6)"), "12");
        assert_eq!(s("lcm(0, 6)"), "0");
        assert_eq!(s("lcm(21, 6)"), "42");
        assert_eq!(s("gcd(gcd(12, 18), 8)"), "2");
    }

    #[test]
    fn mod_is_non_negative() {
        assert_eq!(s("mod(17, 5)"), "2");
        assert_eq!(s("mod(-17, 5)"), "3");
        assert_eq!(s("mod(17, -5)"), "2");
        assert_eq!(s("mod(0, 7)"), "0");
        assert_eq!(s("mod(2^100, 1000)"), "376");
        // Undefined for a zero modulus, so it stays.
        assert_eq!(s("mod(5, 0)"), "mod(5, 0)");
    }

    #[test]
    fn modular_power_and_inverse() {
        assert_eq!(s("powmod(2, 10, 1000)"), "24");
        assert_eq!(s("powmod(2, 1000, 1000000007)"), "688423210");
        assert_eq!(s("powmod(-2, 3, 7)"), "6");
        assert_eq!(s("powmod(5, 0, 7)"), "1");
        assert_eq!(s("powmod(5, 3, 1)"), "0");
        assert_eq!(s("powmod(5, -1, 7)"), "powmod(5, -1, 7)");
        assert_eq!(s("invmod(3, 7)"), "5");
        assert_eq!(s("invmod(-3, 7)"), "2");
        assert_eq!(s("invmod(10, 17)"), "12");
        assert_eq!(s("invmod(4, 8)"), "invmod(4, 8)");
        assert_eq!(s("invmod(3, 0)"), "invmod(3, 0)");
        assert_eq!(s("invmod(5, 1)"), "0");
    }

    #[test]
    fn primality() {
        assert_eq!(s("isprime(2)"), "true");
        assert_eq!(s("isprime(97)"), "true");
        assert_eq!(s("isprime(1)"), "false");
        assert_eq!(s("isprime(0)"), "false");
        assert_eq!(s("isprime(-7)"), "false");
        assert_eq!(s("isprime(91)"), "false");
        // Carmichael number and a strong pseudoprime to small bases.
        assert_eq!(s("isprime(561)"), "false");
        assert_eq!(s("isprime(3215031751)"), "false");
        assert_eq!(s("isprime(2^61 - 1)"), "true");
        assert_eq!(s("isprime(2^64 - 59)"), "true");
        // 2^89 - 1 is a Mersenne prime, 2^67 - 1 famously is not.
        assert_eq!(s("isprime(2^89 - 1)"), "true");
        assert_eq!(s("isprime(2^67 - 1)"), "false");
        assert_eq!(s("isprime(2^127 - 1)"), "true");
    }

    #[test]
    fn next_prime() {
        assert_eq!(s("nextprime(-5)"), "2");
        assert_eq!(s("nextprime(1)"), "2");
        assert_eq!(s("nextprime(2)"), "3");
        assert_eq!(s("nextprime(14)"), "17");
        assert_eq!(s("nextprime(89)"), "97");
        assert_eq!(s("nextprime(10^12)"), "1000000000039");
    }

    #[test]
    fn arithmetic_functions() {
        assert_eq!(s("totient(1)"), "1");
        assert_eq!(s("totient(9)"), "6");
        assert_eq!(s("totient(36)"), "12");
        assert_eq!(s("totient(97)"), "96");
        assert_eq!(s("totient(1000000)"), "400000");
        assert_eq!(s("totient(0)"), "totient(0)");
        assert_eq!(s("divisor_count(1)"), "1");
        assert_eq!(s("divisor_count(12)"), "6");
        assert_eq!(s("divisor_count(360)"), "24");
        assert_eq!(s("divisor_count(97)"), "2");
        assert_eq!(s("divisor_sum(1)"), "1");
        assert_eq!(s("divisor_sum(12)"), "28");
        assert_eq!(s("divisor_sum(28)"), "56");
        assert_eq!(s("divisor_sum(97)"), "98");
    }

    #[test]
    fn factorisation() {
        assert_eq!(
            s("factorint(360)"),
            "list(list(2, 3), list(3, 2), list(5, 1))"
        );
        assert_eq!(s("factorint(1)"), "list()");
        assert_eq!(s("factorint(97)"), "list(list(97, 1))");
        assert_eq!(s("factorint(0)"), "factorint(0)");
        assert_eq!(s("factorint(-6)"), "factorint(-6)");
        // 600851475143 is Project Euler's classic.
        assert_eq!(
            s("factorint(600851475143)"),
            "list(list(71, 1), list(839, 1), list(1471, 1), list(6857, 1))"
        );
        // Semiprime of two 32-bit primes and a prime square.
        assert_eq!(
            s("factorint(4294967291 * 4294967279)"),
            "list(list(4294967279, 1), list(4294967291, 1))"
        );
        assert_eq!(s("factorint(1000003^2)"), "list(list(1000003, 2))");
        assert_eq!(
            s("factorint(2^64 - 1)"),
            "list(list(3, 1), list(5, 1), list(17, 1), list(257, 1), list(641, 1), list(65537, 1), list(6700417, 1))"
        );
        // Beyond 64 bits the kernel declines.
        assert_eq!(s("factorint(2^70)"), "factorint(1180591620717411303424)");
    }

    #[test]
    fn jacobi_symbol() {
        assert_eq!(s("jacobi(1, 1)"), "1");
        assert_eq!(s("jacobi(2, 7)"), "1");
        assert_eq!(s("jacobi(3, 7)"), "-1");
        assert_eq!(s("jacobi(5, 21)"), "1");
        assert_eq!(s("jacobi(6, 9)"), "0");
        assert_eq!(s("jacobi(-1, 7)"), "-1");
        assert_eq!(s("jacobi(1001, 9907)"), "-1");
        assert_eq!(s("jacobi(2, 8)"), "jacobi(2, 8)");
        assert_eq!(s("jacobi(2, -7)"), "jacobi(2, -7)");
    }

    #[test]
    fn chinese_remainder() {
        assert_eq!(s("crt(list(2, 3, 2), list(3, 5, 7))"), "23");
        assert_eq!(s("crt(list(1, 2), list(4, 9))"), "29");
        assert_eq!(s("crt(list(-1, 0), list(5, 3))"), "9");
        assert_eq!(s("crt(list(5), list(7))"), "5");
        assert_eq!(
            s("crt(list(1, 1), list(4, 6))"),
            "crt(list(1, 1), list(4, 6))"
        );
        assert_eq!(s("crt(list(1), list(3, 5))"), "crt(list(1), list(3, 5))");
        assert_eq!(
            s("crt(list(1, x), list(3, 5))"),
            "crt(list(1, x), list(3, 5))"
        );
    }

    #[test]
    fn extended_gcd() {
        assert_eq!(s("egcd(240, 46)"), "list(2, -9, 47)");
        assert_eq!(s("egcd(0, 5)"), "list(5, 0, 1)");
        assert_eq!(s("egcd(7, 0)"), "list(7, 1, 0)");
        assert_eq!(s("egcd(0, 0)"), "list(0, 1, 0)");
        // Bezout's identity through the engine, including negative arguments.
        for a in [-35_i64, -6, 1, 12, 240] {
            for b in [-14_i64, 5, 18, 46] {
                let g = s(&format!("gcd({a}, {b})"));
                let combo = s(&format!("egcd({a}, {b})"));
                let nums: Vec<i64> = combo
                    .trim_start_matches("list(")
                    .trim_end_matches(')')
                    .split(", ")
                    .map(|t| t.parse().unwrap_or(i64::MAX))
                    .collect();
                assert_eq!(nums.first().map(ToString::to_string), Some(g));
                let x = nums.get(1).copied().unwrap_or(0);
                let y = nums.get(2).copied().unwrap_or(0);
                assert_eq!(
                    a * x + b * y,
                    nums.first().copied().unwrap_or(-1),
                    "({a}, {b})"
                );
            }
        }
    }

    #[test]
    fn symbolic_arguments_stay() {
        let sets = [number_theory()];
        for (src, expected) in [
            ("gcd(x, 6)", "gcd(6, x)"),
            ("lcm(x, 6)", "lcm(6, x)"),
            ("mod(x, 6)", "mod(x, 6)"),
            ("powmod(2, x, 7)", "powmod(2, x, 7)"),
            ("isprime(n)", "isprime(n)"),
            ("factorint(n)", "factorint(n)"),
            ("egcd(a, 4)", "egcd(a, 4)"),
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(reduced);
            assert_eq!(text, expected);
        }
    }

    #[test]
    fn rewrites_fire() {
        let sets = [number_theory()];
        let run = |src: &str, assume: &[(&str, Facts)]| reduce_with(&sets, src, assume).0;
        assert_eq!(run("gcd(x, 1)", &[]), "1");
        assert_eq!(run("lcm(x, 0)", &[]), "0");
        assert_eq!(run("mod(x, 1)", &[]), "0");
        assert_eq!(run("mod(mod(x, m), m)", &[]), "mod(x, m)");
        assert_eq!(run("gcd(gcd(a, b), b)", &[]), "gcd(a, b)");
        assert_eq!(run("lcm(lcm(a, b), b)", &[]), "lcm(a, b)");
        // The remaining rules apply to literal integers, where the kernels
        // reach the same value: the rule and the kernel must agree.
        assert_eq!(run("gcd(12, 0)", &[]), "12");
        assert_eq!(run("lcm(12, 1)", &[]), "12");
        assert_eq!(run("gcd(12, 18) * lcm(12, 18)", &[]), "216");
    }

    #[test]
    fn float_semantics() {
        let sets = [number_theory()];
        assert_eq!(eval(&sets, "gcd(12, 18)", &[]), 6.0);
        assert_eq!(eval(&sets, "gcd(x, 18)", &[("x", 11.6)]), 6.0);
        assert_eq!(eval(&sets, "lcm(4, 6)", &[]), 12.0);
        assert_eq!(eval(&sets, "lcm(0, 6)", &[]), 0.0);
        assert_eq!(eval(&sets, "mod(x, 5)", &[("x", -17.0)]), 3.0);
        assert!(eval(&sets, "mod(x, 0)", &[("x", 1.0)]).is_nan());
    }

    /// Splits `list(a, b, ...)` into its top-level items.
    fn items(text: &str) -> Vec<String> {
        let inner = text
            .strip_prefix("list(")
            .and_then(|t| t.strip_suffix(')'))
            .unwrap_or("");
        let (mut out, mut depth, mut current) = (Vec::new(), 0_i32, String::new());
        for ch in inner.chars() {
            match ch {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | _ => {},
            }
            if ch == ',' && depth == 0 {
                out.push(current.trim().to_owned());
                current.clear();
            } else {
                current.push(ch);
            }
        }
        if !current.trim().is_empty() {
            out.push(current.trim().to_owned());
        }
        out
    }

    fn integers(text: &str) -> Vec<i128> {
        items(text)
            .iter()
            .map(|t| t.parse().unwrap_or(i128::MIN))
            .collect()
    }

    #[test]
    fn pell_equation() {
        assert_eq!(s("pell(2)"), "list(3, 2)");
        assert_eq!(s("pell(3)"), "list(2, 1)");
        assert_eq!(s("pell(5)"), "list(9, 4)");
        assert_eq!(s("pell(7)"), "list(8, 3)");
        assert_eq!(s("pell(13)"), "list(649, 180)");
        // The famous hard case.
        assert_eq!(s("pell(61)"), "list(1766319049, 226153980)");
        assert_eq!(s("pell(109)"), "list(158070671986249, 15140424455100)");
        // Every result solves the equation; no smaller y does.
        for n in (2..200_i128).filter(|n| (0..15).all(|r| r * r != *n) && (n.isqrt().pow(2) != *n))
        {
            let v = integers(&s(&format!("pell({n})")));
            let (x, y) = (v[0], v[1]);
            assert_eq!(x * x - n * y * y, 1, "n = {n}");
            if y < 100_000 {
                assert!(
                    (1..y).all(|t| {
                        let x2 = 1 + n * t * t;
                        x2.isqrt().pow(2) != x2
                    }),
                    "n = {n}: y = {y} is not minimal"
                );
            }
        }
        // Squares and non-positive numbers have no non-trivial solution.
        assert_eq!(s("pell(4)"), "pell(4)");
        assert_eq!(s("pell(1)"), "pell(1)");
        assert_eq!(s("pell(0)"), "pell(0)");
        assert_eq!(s("pell(-3)"), "pell(-3)");
        assert_eq!(s("pell(n)"), "pell(n)");
    }

    #[test]
    fn continued_fractions() {
        assert_eq!(s("cfrac(415/93)"), "list(4, 2, 6, 7)");
        assert_eq!(s("cfrac(-7/3)"), "list(-3, 1, 2)");
        assert_eq!(s("cfrac(5)"), "list(5)");
        assert_eq!(s("cfrac(0)"), "list(0)");
        assert_eq!(s("cfrac(1/2)"), "list(0, 2)");
        // sqrt(n) = [a0; period]
        assert_eq!(s("cfrac(2^(1/2))"), "list(1, list(2))");
        assert_eq!(s("cfrac(3^(1/2))"), "list(1, list(1, 2))");
        assert_eq!(s("cfrac(23^(1/2))"), "list(4, list(1, 3, 1, 8))");
        assert_eq!(s("cfrac(9^(1/2))"), "list(3)");
        assert_eq!(
            simplify(
                &[number_theory(), crate::rules::elementary()],
                "cfrac(sqrt(7))"
            ),
            "list(2, list(1, 1, 1, 4))"
        );
        assert_eq!(s("cfrac(x)"), "cfrac(x)");
        assert_eq!(s("cfrac(-2^(1/2))"), "cfrac(-2^(1/2))");
        // The expansion of p/q converges to p/q exactly.
        for (p, q) in [(415_i128, 93_i128), (-355, 113), (22, 7), (1, 3), (144, 89)] {
            let last = s(&format!("convergents({p}/{q}, 30)"));
            let last = items(&last).pop().unwrap_or_default();
            let value = s(&last);
            assert_eq!(value, s(&format!("{p}/{q}")));
        }
    }

    #[test]
    fn convergents_of_expansions() {
        assert_eq!(s("convergents(415/93, 3)"), "list(4, 9/2, 58/13, 415/93)");
        // The list ends with the expansion of a rational.
        assert_eq!(s("convergents(415/93, 10)"), "list(4, 9/2, 58/13, 415/93)");
        assert_eq!(
            s("convergents(2^(1/2), 5)"),
            "list(1, 3/2, 7/5, 17/12, 41/29, 99/70)"
        );
        assert_eq!(
            s("convergents(list(1, 2, 2, 2), 10)"),
            "list(1, 3/2, 7/5, 17/12)"
        );
        assert_eq!(s("convergents(list(3, 7, 15, 1), 1)"), "list(3, 22/7)");
        assert_eq!(s("convergents(x, 3)"), "convergents(x, 3)");
        assert_eq!(s("convergents(2^(1/2), -1)"), "convergents(2^(1/2), -1)");
        // Convergents of sqrt(n) approximate it; the pairs of `pell` appear.
        let c = items(&s("convergents(7^(1/2), 20)"));
        assert!(c.iter().any(|t| t == "8/3"));
        for t in &c {
            let value = eval(&[number_theory()], t, &[]);
            assert!((value - 7.0_f64.sqrt()).abs() < 1.0, "{t}");
        }
        let last = eval(
            &[number_theory()],
            c.last().map_or("0", String::as_str),
            &[],
        );
        assert!((last - 7.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn prime_sieve() {
        assert_eq!(s("primes(30)"), "list(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)");
        assert_eq!(s("primes(2)"), "list(2)");
        assert_eq!(s("primes(1)"), "list()");
        assert_eq!(s("primes(-5)"), "primes(-5)");
        assert_eq!(items(&s("primes(100000)")).len(), 9592);
        assert_eq!(s("primes(1000000000)"), "primes(1000000000)");
        assert_eq!(s("primes(n)"), "primes(n)");
        // Agrees with the primality kernel.
        for p in integers(&s("primes(200)")) {
            assert_eq!(s(&format!("isprime({p})")), "true");
        }
        assert_eq!(
            integers(&s("primes(200)")).len(),
            (0..=200)
                .filter(|&k| s(&format!("isprime({k})")) == "true")
                .count()
        );
    }

    #[test]
    fn linear_diophantine_equations() {
        let sets = [number_theory()];
        let solve = |eq: &str| s(&format!("diophantine({eq}, list(x, y))"));
        assert_eq!(solve("3*x + 5*y = 7"), "list(5*t + 14, -3*t - 7)");
        assert_eq!(solve("4*x + 6*y = 7"), "list()");
        assert_eq!(solve("0*x + 0*y = 3"), "list()");
        assert_eq!(solve("3*x = 6"), "list(2, -t)");
        // Substituting the parametric family back into the equation.
        for (a, b, c) in [
            (3_i128, 5_i128, 7_i128),
            (6, 9, 21),
            (-4, 10, 6),
            (12, -18, 30),
            (7, 0, 14),
            (0, 5, 10),
            (1, 1, 0),
        ] {
            let text = solve(&format!("{a}*x + {b}*y = {c}"));
            let v = items(&text);
            assert_eq!(v.len(), 2, "{text}");
            for t in [-3.0, 0.0, 2.0, 11.0] {
                let value = eval(
                    &sets,
                    &format!("{a}*({}) + {b}*({})", v[0], v[1]),
                    &[("t", t)],
                );
                assert_eq!(value, c as f64, "{a}x + {b}y = {c}, t = {t}: {text}");
            }
        }
        // Rational coefficients are scaled to integers; the equation may
        // be written with everything on one side.
        assert_eq!(solve("x/2 + y/3 = 1"), solve("3*x + 2*y = 6"));
        assert_eq!(solve("3*x + 5*y - 7 = 0"), solve("3*x + 5*y = 7"));
        assert_eq!(solve("y = 2*x + 1"), solve("-2*x + y = 1"));
        // The parameter avoids the unknowns' names.
        let (text, reduced) = reduce_with(&sets, "diophantine(3*t + 5*y = 7, list(t, y))", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn quadratic_diophantine_equations() {
        let sets = [number_theory()];
        assert_eq!(s("diophantine(x^2 - 2*y^2 = 1, list(x, y))"), "list(3, 2)");
        assert_eq!(s("diophantine(y^2 - 7*x^2 = 1, list(x, y))"), "list(3, 8)");
        assert_eq!(
            s("diophantine(x^2 - 61*y^2 = 1, list(x, y))"),
            "list(1766319049, 226153980)"
        );
        assert_eq!(
            s("diophantine(2*x^2 - 6*y^2 = 2, list(x, y))"),
            "list(2, 1)"
        );
        // The negative Pell equation: solvable exactly when the period is odd.
        assert_eq!(s("diophantine(x^2 - 2*y^2 = -1, list(x, y))"), "list(1, 1)");
        assert_eq!(
            s("diophantine(x^2 - 13*y^2 = -1, list(x, y))"),
            "list(18, 5)"
        );
        assert_eq!(s("diophantine(x^2 - 3*y^2 = -1, list(x, y))"), "list()");
        for (n, sign) in [
            (2_i128, 1_i128),
            (5, 1),
            (7, 1),
            (13, 1),
            (41, 1),
            (2, -1),
            (5, -1),
            (13, -1),
            (29, -1),
        ] {
            let v = integers(&s(&format!(
                "diophantine(x^2 - {n}*y^2 = {sign}, list(x, y))"
            )));
            assert_eq!(v.len(), 2, "n = {n}, sign = {sign}");
            assert_eq!(v[0] * v[0] - n * v[1] * v[1], sign);
        }
        // Pythagorean triples: x^2 + y^2 = z^2 in any arrangement.
        let triple = s("diophantine(x^2 + y^2 = z^2, list(x, y, z))");
        assert_eq!(triple, "list(k*(m^2 - n^2), 2*k*m*n, k*(m^2 + n^2))");
        for (m, n, k) in [
            (2.0, 1.0, 1.0),
            (3.0, 2.0, 1.0),
            (5.0, 2.0, 3.0),
            (4.0, 1.0, 2.0),
        ] {
            let v = items(&triple);
            let b = [("m", m), ("n", n), ("k", k)];
            let (x, y, z) = (
                eval(&sets, &v[0], &b),
                eval(&sets, &v[1], &b),
                eval(&sets, &v[2], &b),
            );
            assert_eq!(x * x + y * y, z * z);
        }
        assert_eq!(
            s("diophantine(x^2 - z^2 + y^2 = 0, list(x, z, y))"),
            "list(k*(m^2 - n^2), k*(m^2 + n^2), 2*k*m*n)"
        );
        // Forms that are not recognised stay requests.
        for src in [
            "diophantine(x^2 - 4*y^2 = 1, list(x, y))",
            "diophantine(x^2 - 2*y^2 = 3, list(x, y))",
            "diophantine(x^2 + y^2 = 5, list(x, y))",
            "diophantine(x*y = 6, list(x, y))",
            "diophantine(x^3 + y = 1, list(x, y))",
            "diophantine(a*x + y = 1, list(x, y))",
            "diophantine(x + y = 1, list(x, 2))",
            "diophantine(x + y + z = 1, list(x, y, z))",
            "diophantine(x^2 + y^2 = z^2, list(x, y, m))",
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(!reduced, "{src} => {text}");
        }
    }
}
