//! Definite sums with a symbolic bound by recurrence guessing.
//!
//! For `S(n) = Σ_{k=a}^{n+c} F(n, k)`, which Gosper's algorithm usually
//! cannot close (`Σ binomial(n, k) x^k`, `Σ binomial(n, k)^2`), the exact
//! values `S(0), …, S(N)` are computed as rational functions of at most
//! one further parameter `γ`. A linear recurrence of order at most two
//! with polynomial coefficients in `n` over `Q(γ)` is fitted to them; the
//! fit is overdetermined, so a spurious recurrence is rejected. A
//! first-order recurrence `S(n+1) = K ∏(n + αᵢ)/∏(n + βⱼ) S(n)` is then
//! solved as a product of gamma functions. Every closed form is checked
//! against further exact values before it is returned. This is the
//! approach of guess-and-prove summation with the proof (Zeilberger's
//! certificate) replaced by verification at extra points.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use super::gosper::Frac;
use super::gosper::solve;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

/// Values computed for the fit.
const POINTS: i64 = 14;

/// `S(n)` in closed form, or `None`.
pub(super) fn sum_by_recurrence(
    cx: &mut Cx<'_>,
    term: NodeId,
    k: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let sum = cx.graph.ops().lookup("sum")?;
    let k_symbol = cx.graph.symbol_of(k)?;
    // The bound is n + c with a literal c.
    let bound_symbols = cx.graph.free_symbols(cx.graph.find(upper)).to_vec();
    let &[n_symbol] = bound_symbols.as_slice() else {
        return None;
    };
    let n = cx.graph.symbol_node(n_symbol);
    cx.graph.number_of(lower)?.to_i64()?;
    // At most one other parameter.
    let mut others: Vec<_> = cx.graph.free_symbols(cx.graph.find(term)).to_vec();
    others.retain(|&s| s != k_symbol && s != n_symbol);
    let gamma = match others.as_slice() {
        | [] => None,
        | &[g] => Some(cx.graph.symbol_node(g)),
        | _ => return None,
    };
    let request = cx.graph.node(sum, &[term, k, lower, upper]);

    let value_at = |cx: &mut Cx<'_>, m: i64| -> Option<Frac> {
        let at = cx.graph.int(m);
        let concrete = cx.graph.substitute(request, n, at);
        let value = cx.simplify(concrete);
        frac_of(cx.graph, value, gamma)
    };
    let mut values = Vec::new();
    for m in 0..=POINTS {
        values.push(value_at(cx, m)?);
    }
    let (order, coefficients) = guess(&values)?;
    if order != 1 {
        return None;
    }
    let closed = solve_first_order(cx.graph, &coefficients, &values, n, gamma)?;
    let closed = cx.simplify(closed);
    // Verify at further points.
    for m in POINTS + 1..=POINTS + 3 {
        let want = value_at(cx, m)?;
        let at = cx.graph.int(m);
        let got = cx.graph.substitute(closed, n, at);
        let got = cx.simplify(got);
        if !agrees(cx.graph, got, &want, gamma) {
            return None;
        }
    }
    Some(closed)
}

/// A term as an element of `Q(γ)`.
fn frac_of(
    graph: &mut Graph,
    value: NodeId,
    gamma: Option<NodeId>,
) -> Option<Frac> {
    if let Some(r) = graph.number_of(value).and_then(Number::to_rational) {
        return Some(Frac::poly(vec![r]));
    }
    let g = gamma?;
    let value = expand_binomials(graph, value)?;
    let (n, d) = crate::rules::poly::rational_function_in(graph, value, g)?;
    Some(Frac { n, d }.reduce())
}

/// Rewrites `binomial(x, j)` with a literal `j ≥ 0` as the polynomial
/// `x (x-1) … (x-j+1) / j!`, so values with a symbolic parameter are
/// rational functions.
fn expand_binomials(
    graph: &mut Graph,
    node: NodeId,
) -> Option<NodeId> {
    let term = crate::rules::poly::best(graph, node)?;
    let binomial = graph.ops().lookup("binomial");
    Some(rewrite_binomials(graph, term, binomial))
}

fn rewrite_binomials(
    graph: &mut Graph,
    node: NodeId,
    binomial: Option<crate::graph::OpId>,
) -> NodeId {
    let op = graph.op(node);
    let children = graph.children(node).to_vec();
    if children.is_empty() {
        return node;
    }
    let new_children: Vec<NodeId> = children.iter().map(|&c| rewrite_binomials(graph, c, binomial)).collect();
    if Some(op) == binomial {
        if let &[x, j] = new_children.as_slice() {
            if let Some(j) = graph.number_of(j).and_then(Number::to_i64).filter(|j| (0..=40).contains(j)) {
                let mut factors = Vec::new();
                let mut factorial = BigInt::one();
                for i in 0..j {
                    let minus_i = graph.int(-i);
                    factors.push(graph.node(core::ADD, &[x, minus_i]));
                    factorial *= BigInt::from(i + 1);
                }
                factors.push(graph.num(Number::rat(BigRational::new(BigInt::one(), factorial))));
                return graph.node(core::MUL, &factors);
            }
        }
    }
    if new_children == children { node } else { graph.try_node(op, &new_children).unwrap_or(node) }
}

fn agrees(
    graph: &Graph,
    got: NodeId,
    want: &Frac,
    gamma: Option<NodeId>,
) -> bool {
    let at = BigRational::new(BigInt::from(37), BigInt::from(100));
    let mut env = Env::numeric(0.0);
    if let Some(symbol) = gamma.and_then(|g| graph.symbol_of(g)) {
        env.bind(symbol, 0.37);
    }
    let Some(value) = graph.eval(got, &env) else {
        return false;
    };
    let expected = Number::rat(univariate::eval(&want.n, &at)).to_f64()
        / Number::rat(univariate::eval(&want.d, &at)).to_f64();
    (value - expected).abs() <= 1e-9 * expected.abs().max(1.0)
}

/// Fits `Σ_j P_j(n) S(n+j) = 0`: the order and the coefficient
/// polynomials `P_j` (ascending in `n`, entries over `Q(γ)`).
fn guess(values: &[Frac]) -> Option<(usize, Vec<Vec<Frac>>)> {
    for order in 1..=2_usize {
        for degree in 0..=3_usize {
            let unknowns = (order + 1) * (degree + 1);
            let equations = values.len().checked_sub(order)?;
            if equations < unknowns + 2 {
                continue;
            }
            // Rows: n = 0..equations; columns (j, i): n^i S(n+j).
            let mut rows = Vec::with_capacity(equations);
            for m in 0..equations {
                let mut row = Vec::with_capacity(unknowns);
                for j in 0..=order {
                    let mut power = BigRational::one();
                    for _ in 0..=degree {
                        row.push(values[m + j].mul(&Frac::poly(vec![power.clone()])));
                        power *= BigRational::from_integer(BigInt::from(m));
                    }
                }
                rows.push(row);
            }
            // Normalise one unknown to 1, trying each column of the highest
            // shift in turn.
            for pivot in (order * (degree + 1)..unknowns).rev() {
                let rhs: Vec<Frac> = rows.iter().map(|r| Frac::zero().sub(&r[pivot])).collect();
                let reduced: Vec<Vec<Frac>> = rows
                    .iter()
                    .map(|r| r.iter().enumerate().filter(|&(c, _)| c != pivot).map(|(_, v)| v.clone()).collect())
                    .collect();
                if let Some(x) = solve(reduced, rhs, unknowns - 1) {
                    let mut all = x;
                    all.insert(pivot, Frac::poly(vec![BigRational::one()]));
                    let coefficients: Vec<Vec<Frac>> = all.chunks(degree + 1).map(<[Frac]>::to_vec).collect();
                    return Some((order, coefficients));
                }
            }
        }
    }
    None
}

fn frac_const(r: BigRational) -> Frac {
    Frac::poly(vec![r])
}

fn frac_add(
    a: &Frac,
    b: &Frac,
) -> Frac {
    a.sub(&Frac::zero().sub(b))
}

/// The value of `f` at `γ = g`, if defined.
fn frac_at(
    f: &Frac,
    g: &BigRational,
) -> Option<BigRational> {
    let d = univariate::eval(&f.d, g);
    (!d.is_zero()).then(|| univariate::eval(&f.n, g) / d)
}

/// The constant value of `f`, if it does not involve `γ`.
fn frac_constant(f: &Frac) -> Option<BigRational> {
    match (f.n.as_slice(), f.d.as_slice()) {
        | ([], _) => Some(BigRational::zero()),
        | ([a], [b]) => Some(a / b),
        | _ => None,
    }
}

/// `p / (n - r)` over `Q(γ)` if `r` is a root.
fn divide_root(
    p: &[Frac],
    r: &Frac,
) -> Option<Vec<Frac>> {
    let d = p.len().checked_sub(1)?;
    let mut quotient = vec![Frac::zero(); d];
    let mut carry = Frac::zero();
    for i in (1..=d).rev() {
        carry = frac_add(&p[i], &r.mul(&carry));
        quotient[i - 1] = carry.clone();
    }
    frac_add(&p[0], &r.mul(&carry)).is_zero().then_some(quotient)
}

/// `p(n) = lead · ∏ (n - rᵢ)` with roots `rᵢ` in `Q(γ)` that are at most
/// linear in `γ`; `None` if `p` does not split that way.
fn linear_roots(p: &[Frac]) -> Option<(Frac, Vec<Frac>)> {
    let mut p = p.to_vec();
    while p.last().is_some_and(Frac::is_zero) {
        p.pop();
    }
    let mut roots = Vec::new();
    while p.len() > 1 {
        let root = find_root(&p)?;
        p = divide_root(&p, &root)?;
        roots.push(root);
    }
    Some((p.first()?.clone(), roots))
}

fn find_root(p: &[Frac]) -> Option<Frac> {
    let specialized = |g: &BigRational| -> Option<Vec<BigRational>> {
        let q: Option<QPoly> = p.iter().map(|c| frac_at(c, g)).collect();
        let q = q?;
        if q.last().is_none_or(Zero::is_zero) {
            return None;
        }
        Some(univariate::rational_roots(&q))
    };
    if p.iter().all(|c| frac_constant(c).is_some()) {
        let roots = specialized(&BigRational::zero())?;
        return roots.first().map(|r| frac_const(r.clone()));
    }
    let points: Vec<BigRational> =
        [2, 3, 5, 7, 11, 13, 17].iter().map(|&g| BigRational::from_integer(BigInt::from(g))).collect();
    let sets: Vec<(BigRational, Vec<BigRational>)> =
        points.iter().filter_map(|g| specialized(g).map(|r| (g.clone(), r))).collect();
    let [(g1, r1s), (g2, r2s), ..] = sets.as_slice() else {
        return None;
    };
    for r1 in r1s {
        for r2 in r2s {
            // a + b γ through (g1, r1) and (g2, r2).
            let b = (r2 - r1) / (g2 - g1);
            let a = r1 - &b * g1;
            let candidate = Frac::poly(vec![a, b]);
            if divide_root(p, &candidate).is_some() {
                return Some(candidate);
            }
        }
    }
    None
}

fn frac_term(
    graph: &mut Graph,
    f: &Frac,
    gamma: Option<NodeId>,
) -> Option<NodeId> {
    let poly = |graph: &mut Graph, p: &[BigRational]| -> Option<NodeId> {
        let mut terms = Vec::new();
        for (i, c) in p.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let coefficient = graph.num(Number::rat(c.clone()));
            if i == 0 {
                terms.push(coefficient);
            } else {
                let e = graph.int(i64::try_from(i).ok()?);
                let power = graph.node(core::POW, &[gamma?, e]);
                terms.push(graph.node(core::MUL, &[coefficient, power]));
            }
        }
        Some(match terms.as_slice() {
            | [] => graph.int(0),
            | [one] => *one,
            | _ => graph.node(core::ADD, &terms),
        })
    };
    let n = poly(graph, &f.n)?;
    let d = poly(graph, &f.d)?;
    let minus_one = graph.int(-1);
    let inverse = graph.node(core::POW, &[d, minus_one]);
    Some(graph.node(core::MUL, &[n, inverse]))
}

/// Solves `P1(n) S(n+1) + P0(n) S(n) = 0` from a starting value.
fn solve_first_order(
    graph: &mut Graph,
    coefficients: &[Vec<Frac>],
    values: &[Frac],
    n: NodeId,
    gamma: Option<NodeId>,
) -> Option<NodeId> {
    let [p0, p1] = coefficients else {
        return None;
    };
    // P0 = k0 ∏(n - r), P1 = k1 ∏(n - s): S(n+1)/S(n) = -k0/k1 ∏(n+α)/∏(n+β)
    // with α = -r, β = -s.
    let (k0, r0) = linear_roots(p0)?;
    let (k1, r1) = linear_roots(p1)?;
    let negate = |v: Vec<Frac>| -> Vec<Frac> { v.iter().map(|r| Frac::zero().sub(r)).collect() };
    let (mut alphas, mut betas) = (negate(r0), negate(r1));
    let ratio_k = Frac::zero().sub(&k0).div(&k1)?;
    // A start beyond every pole and zero of the ratio, with S(n0) ≠ 0;
    // parameter-dependent shifts are generic.
    let blocked = |m: i64| {
        let m = BigRational::from_integer(BigInt::from(m));
        alphas.iter().chain(&betas).filter_map(frac_constant).any(|s| {
            let v = &m + s;
            v <= BigRational::zero() && v.is_integer()
        })
    };
    let n0 = (0..i64::try_from(values.len()).ok()?).find(|&m| {
        !blocked(m) && usize::try_from(m).ok().and_then(|i| values.get(i)).is_some_and(|v| !v.is_zero())
    })?;
    let start = values.get(usize::try_from(n0).ok()?)?;
    let gamma_op = graph.ops().lookup("gamma")?;
    let mut factors = vec![frac_term(graph, start, gamma)?];
    let n0_frac = frac_const(BigRational::from_integer(BigInt::from(n0)));
    let minus_one = graph.int(-1);
    let minus_n0 = graph.int(-n0);
    let steps = graph.node(core::ADD, &[n, minus_n0]);
    if frac_constant(&ratio_k).is_none_or(|v| !v.is_one()) {
        let base = frac_term(graph, &ratio_k, gamma)?;
        factors.push(graph.node(core::POW, &[base, steps]));
    }
    // Γ(n+α)/Γ(n+β) with α - β a positive integer d is the polynomial
    // (n+β)(n+β+1)…(n+α-1): pair such shifts off first.
    let mut i = 0;
    while i < alphas.len() {
        let partner = betas.iter().position(|b| frac_constant(&alphas[i].sub(b)).is_some_and(|d| d.is_integer()));
        let Some(j) = partner else {
            i += 1;
            continue;
        };
        let (alpha, beta) = (alphas.swap_remove(i), betas.swap_remove(j));
        let d = frac_constant(&alpha.sub(&beta)).unwrap_or_else(BigRational::zero);
        let (mut low, count, exponent) =
            if d >= BigRational::zero() { (beta, d, 1) } else { (alpha, -d, -1) };
        let mut done = BigRational::zero();
        while done < count {
            let s_term = frac_term(graph, &low, gamma)?;
            let top = graph.node(core::ADD, &[n, s_term]);
            let at_start = frac_term(graph, &frac_add(&n0_frac, &low), gamma)?;
            let inverse = graph.node(core::POW, &[at_start, minus_one]);
            let quotient = graph.node(core::MUL, &[top, inverse]);
            let e = graph.int(exponent);
            factors.push(graph.node(core::POW, &[quotient, e]));
            low = frac_add(&low, &frac_const(BigRational::one()));
            done += BigRational::one();
        }
    }
    for (shifts, exponent) in [(&alphas, 1), (&betas, -1)] {
        for s in shifts {
            let s_term = frac_term(graph, s, gamma)?;
            let top_arg = graph.node(core::ADD, &[n, s_term]);
            let top = graph.node(gamma_op, &[top_arg]);
            let bottom_arg = frac_term(graph, &frac_add(&n0_frac, s), gamma)?;
            let bottom = graph.node(gamma_op, &[bottom_arg]);
            let inverse = graph.node(core::POW, &[bottom, minus_one]);
            let quotient = graph.node(core::MUL, &[top, inverse]);
            let e = graph.int(exponent);
            factors.push(graph.node(core::POW, &[quotient, e]));
        }
    }
    Some(graph.node(core::MUL, &factors))
}
