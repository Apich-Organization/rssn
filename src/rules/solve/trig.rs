//! Trigonometric equations in one unknown.
//!
//! All the circular functions of the unknown must have arguments `k (x + s)`
//! with rational `k` and a common shift `s`. With `θ = g (x + s)` (`g` the
//! rational gcd of the `k`) every argument is an integer multiple `m θ`, and
//! `sin(mθ)`, `cos(mθ)` are polynomials in `S = sin θ`, `C = cos θ`
//! (`tan`, `cot`, `sec`, `csc` are quotients). The numerator `N(S, C)` is
//! reduced modulo `S² + C² = 1` and solved by the first applicable method:
//!
//! * a polynomial in `S` alone: `θ = asin(s) + 2πn` and `π - asin(s) + 2πn`;
//! * a polynomial in `C` alone: `θ = ±acos(c) + 2πn`;
//! * otherwise the Weierstrass substitution `t = tan(θ/2)`: `θ = 2 atan(t) +
//!   2πn`, and `θ = π + 2πn` where `N(0, -1) = 0`.
//!
//! Solutions that are numerically a rational multiple of `π` are written
//! as such, families that differ by a fraction of the period are merged
//! (`0 + 2πn` and `π + 2πn` become `πn`), and every family is checked
//! against the equation for a range of `n`; a family that is valid only for
//! some residues of `n` is restricted to them.
//!
//! In the principal mode only the representative `n = 0` of each family is
//! reported.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::normalize::map_term;
use super::polynomial_roots;
use super::product;
use super::reciprocal;
use super::sample_envs;
use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::SymbolId;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;

const TRIG_NAMES: [&str; 6] = ["sin", "cos", "tan", "cot", "sec", "csc"];

/// One family of solutions `θ = base + period n`, in units where
/// `period` is a rational multiple of `π`.
#[derive(Clone)]
struct Family {
    /// The solution for `n = 0` as a term.
    base: NodeId,
    /// `base / π` when it is a recognised rational.
    rational: Option<BigRational>,
    /// The step in `n`, in units of `π`.
    period: BigRational,
}

fn value(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.eval(node, &Env::numeric(0.0)).filter(|v| v.is_finite())
}

fn rat(
    graph: &mut Graph,
    r: &BigRational,
) -> NodeId {
    graph.num(Number::rat(r.clone()))
}

/// `(k, c)` with `arg = k x + c`, `k` rational.
fn linear_in(
    graph: &mut Graph,
    arg: NodeId,
    x: NodeId,
) -> Option<(BigRational, NodeId)> {
    let symbol = graph.symbol_of(x)?;
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let poly = from_term(graph, &mut gens, arg, Limits::default())?;
    if poly.degree_in(gx) != 1 {
        return None;
    }
    let parts = poly.coefficients_in(gx);
    let k = parts.get(1)?.as_constant()?.to_rational()?;
    let c = to_term(graph, &gens, parts.first()?);
    (!graph.depends_on(graph.find(c), symbol)).then_some((k, c))
}

/// `sin(mθ)` and `cos(mθ)` as polynomials in `S`, `C`.
fn multiple_angle(
    m: i64,
    s: u32,
    c: u32,
) -> Option<(Poly, Poly)> {
    let cap = Limits::default().terms;
    let (ps, pc) = (Poly::generator(s), Poly::generator(c));
    let (mut re, mut im) = (Poly::constant(Number::from(1)), Poly::zero());
    for _ in 0..m.unsigned_abs() {
        let new_re = re.mul(&pc, cap)?.sub(&im.mul(&ps, cap)?);
        let new_im = im.mul(&pc, cap)?.add(&re.mul(&ps, cap)?);
        re = new_re;
        im = new_im;
    }
    Some(if m < 0 { (im.neg(), re) } else { (im, re) })
}

/// `p` with every `C^2` replaced by `1 - S^2`: `A(S) + C B(S)`.
fn reduce_circle(
    p: &Poly,
    s: u32,
    c: u32,
    eliminated: u32,
    kept: u32,
) -> Option<Poly> {
    let cap = Limits::default().terms;
    let _ = (s, c);
    let one_minus = Poly::constant(Number::from(1)).sub(&Poly::generator(kept).pow(2, cap)?);
    let mut out = Poly::zero();
    for (mono, coeff) in p.terms() {
        let e = mono.iter().find(|&&(g, _)| g == eliminated).map_or(0, |&(_, e)| e);
        let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(g, _)| g != eliminated).collect();
        let (q, r) = (e / 2, e % 2);
        let mut piece = Poly::monomial(rest, coeff.clone());
        if r == 1 {
            piece = piece.mul(&Poly::generator(eliminated), cap)?;
        }
        piece = piece.mul(&one_minus.pow(q, cap)?, cap)?;
        out = out.add(&piece);
    }
    Some(out)
}

/// The part of `p` with `v^e` stripped (`e` 0 or 1), as polynomials without `v`.
fn split_parity(
    p: &Poly,
    v: u32,
) -> (Poly, Poly) {
    let (mut even, mut odd) = (Poly::zero(), Poly::zero());
    for (mono, coeff) in p.terms() {
        let e = mono.iter().find(|&&(g, _)| g == v).map_or(0, |&(_, e)| e);
        let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(g, _)| g != v).collect();
        let piece = Poly::monomial(rest, coeff.clone());
        if e == 0 {
            even = even.add(&piece);
        } else {
            odd = odd.add(&piece);
        }
    }
    (even, odd)
}

/// A rational `p/q` with `q <= 120` within `1e-10` of `v`.
fn recognise(v: f64) -> Option<BigRational> {
    for q in 1..=120_i64 {
        let scaled = v * f64::from(i32::try_from(q).ok()?);
        let nearest = scaled.round();
        if (scaled - nearest).abs() < 1e-9 {
            return Some(BigRational::new(BigInt::from(nearest as i64), BigInt::from(q)));
        }
    }
    None
}

/// All solutions of `term = 0` as families, or `None` when `term` is not a
/// trigonometric equation of the supported shape.
fn families(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
) -> Option<(Vec<Family>, BigRational, NodeId)> {
    let symbol = graph.symbol_of(x)?;
    let ops: Vec<Option<OpId>> = TRIG_NAMES.iter().map(|n| graph.ops().lookup(n)).collect();
    let is_trig = |g: &Graph, n: NodeId| ops.iter().any(|&o| o == Some(g.op(n)));
    let trig_nodes = super::heuristics::nodes_with(graph, term, |g, n| is_trig(g, n) && g.depends_on(g.find(n), symbol));
    if trig_nodes.is_empty() {
        return None;
    }
    // Arguments k (x + s).
    let mut lines: Vec<(BigRational, NodeId)> = Vec::new();
    for &t in &trig_nodes {
        let arg = *graph.children(t).first()?;
        lines.push(linear_in(graph, arg, x)?);
    }
    let shifts: Vec<NodeId> = {
        let mut out = Vec::new();
        for (k, c) in &lines {
            let inverse = rat(graph, &k.recip());
            out.push(product(graph, &[*c, inverse]));
        }
        out
    };
    let shift = *shifts.first()?;
    for &other in shifts.iter().skip(1) {
        let diff = {
            let minus_one = graph.int(-1);
            let neg = graph.node(core::MUL, &[minus_one, other]);
            graph.node(core::ADD, &[shift, neg])
        };
        let mut equal = value(graph, diff).is_some_and(|v| v.abs() < 1e-12);
        if !equal && value(graph, diff).is_none() {
            equal = sample_envs(graph, &[diff], None).iter().all(|env| graph.eval(diff, env).is_none_or(|v| v.abs() < 1e-9));
        }
        if !equal {
            return None;
        }
    }
    // The rational gcd g of the multipliers.
    let mut denominator = BigInt::one();
    for (k, _) in &lines {
        denominator = denominator.lcm(k.denom());
    }
    let mut numerator = BigInt::zero();
    for (k, _) in &lines {
        let scaled = (k * BigRational::from_integer(denominator.clone())).to_integer();
        numerator = numerator.gcd(&scaled);
    }
    if numerator.is_zero() {
        return None;
    }
    let g = BigRational::new(numerator, denominator);
    // Replace each trigonometric node by its expression in S and C.
    let s_symbol = graph.interner_mut().fresh_symbol("S");
    let c_symbol = graph.interner_mut().fresh_symbol("C");
    let (s_node, c_node) = (graph.symbol_node(s_symbol), graph.symbol_node(c_symbol));
    let mut gens = Gens::default();
    let gs = gens.index(graph, s_node);
    let gc = gens.index(graph, c_node);
    let mut failed = false;
    let mut replace = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
        let which = ops.iter().position(|&o| o == Some(graph.op(node)))?;
        let &[arg] = children else {
            return None;
        };
        if !graph.depends_on(graph.find(arg), symbol) {
            return None;
        }
        let Some((k, _)) = linear_in(graph, arg, x) else {
            failed = true;
            return None;
        };
        let m = (k / g.clone()).to_integer().to_i64().filter(|m| m.abs() <= 24);
        let Some(m) = m else {
            failed = true;
            return None;
        };
        let Some((sin_p, cos_p)) = multiple_angle(m, gs, gc) else {
            failed = true;
            return None;
        };
        let (sin_t, cos_t) = (to_term(graph, &gens, &sin_p), to_term(graph, &gens, &cos_p));
        Some(match which {
            | 0 => sin_t,
            | 1 => cos_t,
            | 2 => {
                let inv = reciprocal(graph, cos_t);
                product(graph, &[sin_t, inv])
            },
            | 3 => {
                let inv = reciprocal(graph, sin_t);
                product(graph, &[cos_t, inv])
            },
            | 4 => reciprocal(graph, cos_t),
            | _ => reciprocal(graph, sin_t),
        })
    };
    let substituted = map_term(graph, term, &mut replace);
    if failed || graph.depends_on(graph.find(substituted), symbol) {
        return None;
    }
    let fraction = ratio(graph, &mut gens, substituted, Limits::default())?;
    let numer = fraction.numer;
    if numer.is_zero() {
        return None;
    }
    // Solve N(S, C) = 0 for θ.
    let mut theta: Vec<Family> = Vec::new();
    let pi = graph.ops().lookup("pi")?;
    let pi_node = graph.node(pi, &[]);
    let asin = graph.ops().lookup("asin")?;
    let acos = graph.ops().lookup("acos")?;
    let atan = graph.ops().lookup("atan")?;
    let two_pi = BigRational::from_integer(BigInt::from(2));
    let by_s = reduce_circle(&numer, gs, gc, gc, gs)?;
    let (a_s, b_s) = split_parity(&by_s, gc);
    let by_c = reduce_circle(&numer, gs, gc, gs, gc)?;
    let (a_c, b_c) = split_parity(&by_c, gs);
    let add_theta = |graph: &mut Graph, base: NodeId, period: &BigRational, theta: &mut Vec<Family>| {
        theta.push(Family { base, rational: None, period: period.clone() });
        let _ = graph;
    };
    if b_s.is_zero() && !a_s.is_zero() {
        // Polynomial in S.
        if a_s.degree_in(gs) == 0 {
            return Some((Vec::new(), g, shift));
        }
        for r in polynomial_roots(graph, &gens, &a_s, gs)? {
            if value(graph, r).is_some_and(|v| v.abs() > 1.0 + 1e-12) {
                continue;
            }
            let a = graph.node(asin, &[r]);
            add_theta(graph, a, &two_pi, &mut theta);
            let minus_one = graph.int(-1);
            let neg = graph.node(core::MUL, &[minus_one, a]);
            let reflected = graph.node(core::ADD, &[pi_node, neg]);
            add_theta(graph, reflected, &two_pi, &mut theta);
        }
    } else if a_s.is_zero() && !b_s.is_zero() {
        // C * B(S): C = 0 or B(S) = 0.
        let half = BigRational::new(BigInt::one(), BigInt::from(2));
        let half_coefficient = rat(graph, &half);
        let half_pi = product(graph, &[half_coefficient, pi_node]);
        add_theta(graph, half_pi, &BigRational::one(), &mut theta);
        if b_s.degree_in(gs) > 0 {
            for r in polynomial_roots(graph, &gens, &b_s, gs)? {
                if value(graph, r).is_some_and(|v| v.abs() > 1.0 + 1e-12) {
                    continue;
                }
                let a = graph.node(asin, &[r]);
                add_theta(graph, a, &two_pi, &mut theta);
                let minus_one = graph.int(-1);
                let neg = graph.node(core::MUL, &[minus_one, a]);
                let reflected = graph.node(core::ADD, &[pi_node, neg]);
                add_theta(graph, reflected, &two_pi, &mut theta);
            }
        }
    } else if b_c.is_zero() && !a_c.is_zero() {
        if a_c.degree_in(gc) == 0 {
            return Some((Vec::new(), g, shift));
        }
        for r in polynomial_roots(graph, &gens, &a_c, gc)? {
            if value(graph, r).is_some_and(|v| v.abs() > 1.0 + 1e-12) {
                continue;
            }
            let a = graph.node(acos, &[r]);
            add_theta(graph, a, &two_pi, &mut theta);
            let minus_one = graph.int(-1);
            let neg = graph.node(core::MUL, &[minus_one, a]);
            add_theta(graph, neg, &two_pi, &mut theta);
        }
    } else {
        // Weierstrass substitution.
        let t_symbol = graph.interner_mut().fresh_symbol("t");
        let t_node = graph.symbol_node(t_symbol);
        let gt = gens.index(graph, t_node);
        let cap = Limits::default().terms;
        let degree = numer.terms().map(|(m, _)| m.iter().map(|&(_, e)| e).sum::<u32>()).max().unwrap_or(0);
        let tp = Poly::generator(gt);
        let two_t = tp.scale(&Number::from(2));
        let one = Poly::constant(Number::from(1));
        let t2 = tp.pow(2, cap)?;
        let one_minus = one.sub(&t2);
        let one_plus = one.add(&t2);
        let mut total = Poly::zero();
        let mut at_pi = Poly::zero();
        for (mono, coeff) in numer.terms() {
            let a = mono.iter().find(|&&(g, _)| g == gs).map_or(0, |&(_, e)| e);
            let b = mono.iter().find(|&&(g, _)| g == gc).map_or(0, |&(_, e)| e);
            let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(g, _)| g != gs && g != gc).collect();
            let mut piece = Poly::monomial(rest.clone(), coeff.clone());
            piece = piece.mul(&two_t.pow(a, cap)?, cap)?;
            piece = piece.mul(&one_minus.pow(b, cap)?, cap)?;
            piece = piece.mul(&one_plus.pow(degree.saturating_sub(a + b), cap)?, cap)?;
            total = total.add(&piece);
            if a == 0 {
                let sign = if b % 2 == 0 { 1 } else { -1 };
                at_pi = at_pi.add(&Poly::monomial(rest, coeff.mul(&Number::from(sign))));
            }
        }
        if total.degree_in(gt) > 0 {
            let two = graph.int(2);
            for r in polynomial_roots(graph, &gens, &total, gt)? {
                let a = graph.node(atan, &[r]);
                let doubled = product(graph, &[two, a]);
                add_theta(graph, doubled, &two_pi, &mut theta);
            }
        }
        if at_pi.is_zero() {
            add_theta(graph, pi_node, &two_pi, &mut theta);
        }
    }
    Some((theta, g, shift))
}

/// Writes numeric bases as rational multiples of `π` where they are.
fn recognise_families(
    graph: &mut Graph,
    theta: Vec<Family>,
) -> Vec<Family> {
    let Some(pi) = graph.ops().lookup("pi") else {
        return theta;
    };
    let pi_node = graph.node(pi, &[]);
    let mut out = Vec::with_capacity(theta.len());
    for mut f in theta {
        if let Some(v) = value(graph, f.base)
            && let Some(r) = recognise(v / std::f64::consts::PI) {
                let node = if r.is_zero() {
                    graph.int(0)
                } else {
                    let coefficient = rat(graph, &r);
                    product(graph, &[coefficient, pi_node])
                };
                f.base = node;
                f.rational = Some(r);
            }
        out.push(f);
    }
    out
}

/// Merges families that are equal modulo their period or that together
/// make up a finer family.
fn merge(
    graph: &mut Graph,
    theta: Vec<Family>,
) -> Vec<Family> {
    let Some(pi) = graph.ops().lookup("pi") else {
        return theta;
    };
    let pi_node = graph.node(pi, &[]);
    let mut exact: Vec<(BigRational, BigRational)> = Vec::new(); // (residue, period), units of pi
    let mut other: Vec<Family> = Vec::new();
    for f in theta {
        match &f.rational {
            | Some(r) => {
                let residue = r - (r / &f.period).floor() * &f.period;
                if !exact.iter().any(|(r0, p0)| *p0 == f.period && *r0 == residue) {
                    exact.push((residue, f.period.clone()));
                }
            },
            | None => {
                let v = value(graph, f.base);
                let duplicate = other.iter().any(|o| match (value(graph, o.base), v) {
                    | (Some(a), Some(b)) => {
                        let step = f.period.to_f64().unwrap_or(2.0) * std::f64::consts::PI;
                        ((a - b) / step - ((a - b) / step).round()).abs() < 1e-9
                    },
                    | _ => graph.same(o.base, f.base),
                });
                if !duplicate {
                    other.push(f);
                }
            },
        }
    }
    // Residues with period p merge into period p/k when all of them occur.
    let mut changed = true;
    while changed {
        changed = false;
        'search: for k in [12_i64, 8, 6, 5, 4, 3, 2] {
            for i in 0..exact.len() {
                let (r, p) = exact[i].clone();
                let step = &p / BigRational::from_integer(BigInt::from(k));
                let members: Vec<BigRational> = (0..k).map(|j| &r + &step * BigRational::from_integer(BigInt::from(j))).collect();
                let reduce = |v: &BigRational| v - (v / &p).floor() * &p;
                let positions: Option<Vec<usize>> = members
                    .iter()
                    .map(|m| exact.iter().position(|(r0, p0)| *p0 == p && *r0 == reduce(m)))
                    .collect();
                if let Some(positions) = positions {
                    let mut sorted = positions;
                    sorted.sort_unstable();
                    sorted.dedup();
                    if sorted.len() == usize::try_from(k).unwrap_or(0) {
                        for &pos in sorted.iter().rev() {
                            exact.remove(pos);
                        }
                        let base = r.clone() - (&r / &step).floor() * &step;
                        exact.push((base, step));
                        changed = true;
                        break 'search;
                    }
                }
            }
        }
    }
    exact.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<Family> = exact
        .into_iter()
        .map(|(r, p)| {
            let base = if r.is_zero() {
                graph.int(0)
            } else {
                let coefficient = rat(graph, &r);
                product(graph, &[coefficient, pi_node])
            };
            Family { base, rational: Some(r), period: p }
        })
        .collect();
    out.extend(other);
    out
}

/// The symbol for the integer parameter `n` (`n`, `n1`, ... when taken).
pub(super) fn integer_symbol(
    graph: &mut Graph,
    taken: &[NodeId],
) -> NodeId {
    let mut candidate = 0_u32;
    loop {
        let name = if candidate == 0 { "n".to_string() } else { format!("n{candidate}") };
        let symbol = graph.interner_mut().symbol(&name);
        let used = taken.iter().any(|&t| graph.depends_on(graph.find(t), symbol));
        if !used {
            graph.assume(symbol, Facts::INTEGER);
            return graph.symbol_node(symbol);
        }
        candidate += 1;
    }
}

/// The solutions of the trigonometric equation `term = 0`: in the general
/// mode the families `base + step n`, otherwise their representatives.
pub(super) fn solve(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    general: bool,
) -> Option<Vec<NodeId>> {
    let (theta, g, shift) = families(graph, term, x)?;
    let theta = recognise_families(graph, theta);
    let theta = merge(graph, theta);
    let pi = graph.ops().lookup("pi")?;
    let pi_node = graph.node(pi, &[]);
    let n_node = if general { Some(integer_symbol(graph, &[term, x, shift])) } else { None };
    let inverse_g = g.recip();
    let minus_one = graph.int(-1);
    let neg_shift = graph.node(core::MUL, &[minus_one, shift]);
    let symbol = graph.symbol_of(x)?;
    let mut out = Vec::new();
    for family in theta {
        // x = (base + period π n) / g - s.
        let scale = rat(graph, &inverse_g);
        let first = product(graph, &[scale, family.base]);
        let base = graph.node(core::ADD, &[first, neg_shift]);
        let step_coefficient = rat(graph, &(&family.period * &inverse_g));
        let step_pi = product(graph, &[step_coefficient, pi_node]);
        // Validate the family against the equation.
        let n_validation = integer_symbol_for_validation(graph, term);
        let with_n = {
            let times = product(graph, &[step_pi, n_validation]);
            graph.node(core::ADD, &[base, times])
        };
        let valid = valid_residues(graph, term, symbol, with_n, n_validation);
        let Some((modulus, residues)) = valid else {
            continue;
        };
        for r in residues {
            // x = base + step (modulus m + r).
            let shifted_base = if r == 0 {
                base
            } else {
                let offset = {
                    let r_node = graph.int(i64::from(r));
                    product(graph, &[step_pi, r_node])
                };
                graph.node(core::ADD, &[base, offset])
            };
            match n_node {
                | Some(n) => {
                    let m_step = {
                        let m_node = graph.int(i64::from(modulus));
                        product(graph, &[step_pi, m_node])
                    };
                    let times = product(graph, &[m_step, n]);
                    out.push(graph.node(core::ADD, &[shifted_base, times]));
                },
                | None => out.push(shifted_base),
            }
        }
    }
    Some(out)
}

/// A scratch integer symbol for testing the members of a family.
fn integer_symbol_for_validation(
    graph: &mut Graph,
    term: NodeId,
) -> NodeId {
    let _ = term;
    let symbol = graph.interner_mut().fresh_symbol("m");
    graph.assume(symbol, Facts::INTEGER);
    graph.symbol_node(symbol)
}

/// For a family `member(n)`: the smallest modulus `M` (1, 2, 3, 4, 6, 12) such
/// that membership in the solution set depends only on `n mod M`, and the
/// residues that are solutions. `None` when no member is a solution.
fn valid_residues(
    graph: &mut Graph,
    term: NodeId,
    symbol: SymbolId,
    member: NodeId,
    n: NodeId,
) -> Option<(u32, Vec<u32>)> {
    let n_symbol = graph.symbol_of(n)?;
    let substituted = {
        let x = graph.symbol_node(symbol);
        graph.substitute(term, x, member)
    };
    let envs = sample_envs(graph, &[substituted, member], Some(n_symbol));
    let range: Vec<i64> = (-12..=12).collect();
    // `None` where the equation is undefined at every sample assignment.
    let holds = |graph: &Graph, k: i64| -> Option<bool> {
        let mut defined = false;
        for env in &envs {
            let mut env = env.clone();
            env.bind(n_symbol, f64::from(i32::try_from(k).unwrap_or(0)));
            match graph.eval(substituted, &env) {
                | Some(v) if v.is_finite() => {
                    defined = true;
                    let at = graph.eval(member, &env).unwrap_or(0.0);
                    if v.abs() > 1e-7 * (1.0 + at.abs().powi(2)) {
                        return Some(false);
                    }
                },
                | _ => {},
            }
        }
        defined.then_some(true)
    };
    let flags: Vec<Option<bool>> = range.iter().map(|&k| holds(graph, k)).collect();
    if flags.iter().all(Option::is_none) {
        // Undefined at every spot check (parameters outside the range of
        // the inverse functions): keep the family, it may be real elsewhere.
        return Some((1, vec![0]));
    }
    for modulus in [1_u32, 2, 3, 4, 6, 12] {
        let class = |r: i64| -> Vec<Option<bool>> {
            range.iter().zip(&flags).filter(|&(&k, _)| k.rem_euclid(i64::from(modulus)) == r).map(|(_, &f)| f).collect()
        };
        let consistent = (0..i64::from(modulus)).all(|r| {
            let members = class(r);
            let definite: Vec<bool> = members.iter().filter_map(|&f| f).collect();
            definite.windows(2).all(|w| w[0] == w[1])
        });
        if consistent {
            let residues: Vec<u32> = (0..modulus)
                .filter(|&r| class(i64::from(r)).contains(&Some(true)))
                .collect();
            return if residues.is_empty() { None } else { Some((modulus, residues)) };
        }
    }
    None
}
