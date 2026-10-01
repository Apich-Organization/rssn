//! Lie point symmetries, similarity reductions and Noether currents.
//!
//! Derivatives of `u(x₁, …, xₙ)` are represented by jet symbols `u_J`; the
//! total derivative `D_i f = ∂f/∂xᵢ + Σ_J u_{J+eᵢ} ∂f/∂u_J` acts on any
//! expression in the independent variables and the jets.
//!
//! **Symmetries.** `X = Σ ξᵢ ∂_{xᵢ} + φ ∂_u` is a symmetry of `Δ = 0`
//! when `pr X(Δ)` vanishes on solutions. With the characteristic
//! `Q = φ - Σ ξᵢ u_{eᵢ}`, the prolonged coefficients are
//! `φ^J = D_J Q + Σ ξᵢ u_{J+eᵢ}`. The equation is solved for a leading
//! derivative (the highest time derivative, or `u_{xx}` for stationary
//! equations), which with its differential consequences is substituted
//! into `pr X(Δ)`; what remains must vanish identically in the free jets.
//! The infinitesimals are sought as polynomials of degree ≤ 2 in the
//! independent variables (`φ` also linear in `u`), so the determining
//! equations become one homogeneous linear system, whose null space is
//! the symmetry algebra found (translations, scalings, Galilean boosts,
//! projective maps, superposition of polynomial solutions, …).
//!
//! **Reductions.** A traveling wave `u = F(x - c t)` and a similarity
//! solution of a symmetry with `ξ = α x + β`, `τ = γ t + δ`, `φ = κ u`
//! (invariants `z` and `u = g(t) F(z)` in closed form) turn the PDE into
//! an ODE, which is handed to the ODE solver.
//!
//! **Noether.** For a Lagrangian density `L(x, t, u, u_x, u_t)` and a
//! variational symmetry, `J^t = τ L + Q ∂L/∂u_t`, `J^x = ξ L + Q ∂L/∂u_x`
//! is conserved: `D_t J^t + D_x J^x = 0` on solutions.

use std::collections::HashMap;

use super::Index;
use super::Problem;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::sub;

fn div(
    graph: &mut crate::graph::Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus = graph.int(-1);
    let inverse = graph.node(core::POW, &[b, minus]);
    mul(graph, &[a, inverse])
}

/// Jet symbols `u_J` for one problem.
pub(super) struct Jets {
    pub(super) vars: Vec<NodeId>,
    by_index: HashMap<Index, NodeId>,
    by_symbol: HashMap<SymbolId, Index>,
}

impl Jets {
    pub(super) fn new(vars: &[NodeId]) -> Self {
        Self { vars: vars.to_vec(), by_index: HashMap::new(), by_symbol: HashMap::new() }
    }

    pub(super) fn get(
        &mut self,
        cx: &mut Cx<'_>,
        index: &Index,
    ) -> NodeId {
        if let Some(&n) = self.by_index.get(index) {
            return n;
        }
        let name: String = index.iter().map(u32::to_string).collect::<Vec<_>>().join("_");
        let symbol = cx.graph.interner_mut().fresh_symbol(&format!("u_{name}"));
        let node = cx.graph.symbol_node(symbol);
        self.by_index.insert(index.clone(), node);
        self.by_symbol.insert(symbol, index.clone());
        node
    }

    /// The jets occurring in `f`, with their indices.
    pub(super) fn occurring(
        &self,
        cx: &Cx<'_>,
        f: NodeId,
    ) -> Vec<(Index, NodeId)> {
        cx.graph
            .free_symbols(cx.graph.find(f))
            .iter()
            .filter_map(|s| self.by_symbol.get(s).map(|i| (i.clone(), self.by_index[i])))
            .collect()
    }

    /// `D_k f`.
    pub(super) fn total(
        &mut self,
        cx: &mut Cx<'_>,
        f: NodeId,
        k: usize,
    ) -> Option<NodeId> {
        let mut terms = vec![derivative(cx.graph, f, self.vars[k])?];
        for (index, node) in self.occurring(cx, f) {
            let mut next = index.clone();
            next[k] += 1;
            let higher = self.get(cx, &next);
            let partial = derivative(cx.graph, f, node)?;
            terms.push(mul(cx.graph, &[higher, partial]));
        }
        let sum = add(cx.graph, &terms);
        Some(cx.simplify(sum))
    }

    /// `D_J f`.
    pub(super) fn total_multi(
        &mut self,
        cx: &mut Cx<'_>,
        f: NodeId,
        index: &[u32],
    ) -> Option<NodeId> {
        let mut out = f;
        for (k, &times) in index.iter().enumerate() {
            for _ in 0..times {
                out = self.total(cx, out, k)?;
            }
        }
        Some(out)
    }
}

/// The problem's residual with every derivative of `u` (and `u`) replaced
/// by its jet symbol.
pub(super) fn to_jets(
    cx: &mut Cx<'_>,
    p: &Problem,
    jets: &mut Jets,
) -> NodeId {
    let mut residual = p.residual;
    let mut occurring = p.jets.clone();
    // Highest order first, so u_xx is replaced before u_x and u.
    occurring.sort_by_key(|(i, _)| std::cmp::Reverse(i.iter().sum::<u32>()));
    for (index, node) in occurring {
        let symbol = jets.get(cx, &index);
        residual = cx.graph.replace_subterm(residual, node, symbol);
    }
    let u0 = jets.get(cx, &vec![0; p.vars.len()]);
    cx.graph.replace_subterm(residual, p.unknown, u0)
}

/// `Δ = 0` solved for a leading derivative: `(index, expression)`.
fn leading(
    cx: &mut Cx<'_>,
    p: &Problem,
    jets: &mut Jets,
    delta: NodeId,
) -> Option<(Index, NodeId)> {
    let t = p.time(cx.graph);
    let mut candidates: Vec<Index> = p.jets.iter().map(|(i, _)| i.clone()).collect();
    // Prefer the highest pure time derivative, then the highest pure
    // derivative in the first variable.
    candidates.sort_by_key(|i| {
        let pure_t = i.iter().enumerate().all(|(k, &e)| k == t || e == 0) && i[t] > 0;
        let pure_x = i.iter().enumerate().all(|(k, &e)| k == 0 || e == 0) && i[0] > 0;
        (std::cmp::Reverse(u8::from(pure_t)), std::cmp::Reverse(u8::from(pure_x)), std::cmp::Reverse(i.iter().sum::<u32>()))
    });
    for index in candidates {
        let symbol = jets.get(cx, &index);
        let Some(parts) = crate::rules::ode::coefficients_of(cx.graph, delta, symbol) else {
            continue;
        };
        let [c0, c1] = parts.as_slice() else {
            continue;
        };
        // A constant coefficient keeps the substitution polynomial.
        if cx.graph.number_of(*c1).is_none_or(Number::is_zero) {
            continue;
        }
        let minus = cx.graph.int(-1);
        let quotient = div(cx.graph, *c0, *c1);
        let value = mul(cx.graph, &[minus, quotient]);
        return Some((index, cx.simplify(value)));
    }
    None
}

/// `f` on solutions: every jet that is a derivative of the leading one
/// replaced by the corresponding total derivative of its value.
fn on_shell(
    cx: &mut Cx<'_>,
    jets: &mut Jets,
    lead: &(Index, NodeId),
    f: NodeId,
) -> Option<NodeId> {
    let mut current = f;
    for _ in 0..8 {
        let mut changed = false;
        for (index, node) in jets.occurring(cx, current) {
            if index == lead.0 || !index.iter().zip(&lead.0).all(|(a, b)| a >= b) {
                continue;
            }
            let rest: Vec<u32> = index.iter().zip(&lead.0).map(|(a, b)| a - b).collect();
            let value = jets.total_multi(cx, lead.1, &rest)?;
            current = cx.graph.substitute(current, node, value);
            changed = true;
        }
        let lead_node = jets.get(cx, &lead.0);
        if crate::rules::ode::occurs_in(cx.graph, current, lead_node) {
            current = cx.graph.substitute(current, lead_node, lead.1);
            changed = true;
        }
        if !changed {
            return Some(current);
        }
    }
    None
}

/// Monomials of degree ≤ `degree` in the variables.
fn monomials(
    cx: &mut Cx<'_>,
    vars: &[NodeId],
    degree: u32,
) -> Vec<NodeId> {
    let mut out = vec![cx.graph.int(1)];
    let mut frontier = vec![(cx.graph.int(1), 0_usize, 0_u32)];
    while let Some((m, start, d)) = frontier.pop() {
        if d == degree {
            continue;
        }
        for (k, &v) in vars.iter().enumerate().skip(start) {
            let next = mul(cx.graph, &[m, v]);
            out.push(next);
            frontier.push((next, k, d + 1));
        }
    }
    out
}

/// Polynomial point symmetries `list(ξ₁, …, ξₙ, φ)`.
pub(super) fn symmetries(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Vec<Vec<NodeId>>> {
    let mut jets = Jets::new(&p.vars);
    let delta = to_jets(cx, p, &mut jets);
    let delta = cx.simplify(delta);
    let lead = leading(cx, p, &mut jets, delta)?;
    let n = p.vars.len();
    let u0 = jets.get(cx, &vec![0; n]);
    let basis = monomials(cx, &p.vars, 2);
    let mut unknowns = Vec::new();
    let mut fresh = |cx: &mut Cx<'_>| {
        let s = cx.graph.interner_mut().fresh_symbol("k");
        let node = cx.graph.symbol_node(s);
        unknowns.push(node);
        node
    };
    let mut xi = Vec::with_capacity(n);
    for _ in 0..n {
        let terms: Vec<NodeId> = basis.clone().into_iter().map(|m| {
            let c = fresh(cx);
            mul(cx.graph, &[c, m])
        }).collect();
        xi.push(add(cx.graph, &terms));
    }
    let phi = {
        let mut terms = Vec::new();
        for &m in &basis {
            let a = fresh(cx);
            let b = fresh(cx);
            terms.push(mul(cx.graph, &[a, m]));
            terms.push(mul(cx.graph, &[b, m, u0]));
        }
        add(cx.graph, &terms)
    };
    // Q = φ - Σ ξᵢ u_{eᵢ}
    let mut q_terms = vec![phi];
    for (i, &x) in xi.iter().enumerate() {
        let mut e = vec![0; n];
        e[i] = 1;
        let ui = jets.get(cx, &e);
        let minus = cx.graph.int(-1);
        q_terms.push(mul(cx.graph, &[minus, x, ui]));
    }
    let q = add(cx.graph, &q_terms);
    // pr X Δ = Σ ξᵢ ∂Δ/∂xᵢ + Σ_J φ^J ∂Δ/∂u_J
    let mut terms = Vec::new();
    for (i, &x) in xi.iter().enumerate() {
        let partial = derivative(cx.graph, delta, p.vars[i])?;
        terms.push(mul(cx.graph, &[x, partial]));
    }
    for (index, node) in jets.occurring(cx, delta) {
        let partial = derivative(cx.graph, delta, node)?;
        let partial = cx.simplify(partial);
        if cx.is_zero(partial) {
            continue;
        }
        let mut phi_j = vec![jets.total_multi(cx, q, &index)?];
        for (i, &x) in xi.iter().enumerate() {
            let mut next = index.clone();
            next[i] += 1;
            let higher = jets.get(cx, &next);
            phi_j.push(mul(cx.graph, &[x, higher]));
        }
        let phi_j = add(cx.graph, &phi_j);
        terms.push(mul(cx.graph, &[phi_j, partial]));
    }
    let condition = add(cx.graph, &terms);
    let condition = on_shell(cx, &mut jets, &lead, condition)?;
    let solutions = crate::rules::ode::solve_linear_identity(cx.graph, condition, &unknowns)?;
    let mut out = Vec::new();
    for v in solutions {
        let mut generator = Vec::with_capacity(n + 1);
        for &f in xi.iter().chain(std::iter::once(&phi)) {
            let mut value = f;
            for (u, c) in unknowns.iter().zip(&v) {
                let c = cx.graph.num(Number::rat(c.clone()));
                value = cx.graph.substitute(value, *u, c);
            }
            value = cx.graph.substitute(value, u0, p.unknown);
            generator.push(cx.simplify(value));
        }
        if generator.iter().any(|&g| !cx.is_zero(g)) {
            out.push(generator);
        }
    }
    Some(out)
}

/// Replaces jets `u_J` in `f` by `diff` requests of `F(s)` given
/// `u_J = coefficient(J) · F^(|J|)(s)`.
fn jets_to_ode(
    cx: &mut Cx<'_>,
    jets: &Jets,
    f: NodeId,
    big_f: NodeId,
    s: NodeId,
    coefficient: impl Fn(&mut Cx<'_>, &Index) -> NodeId,
) -> Option<NodeId> {
    let diff = cx.graph.ops().lookup("diff")?;
    let mut out = f;
    for (index, node) in jets.occurring(cx, f) {
        let order: u32 = index.iter().sum();
        let mut d = big_f;
        for _ in 0..order {
            d = cx.graph.node(diff, &[d, s]);
        }
        let c = coefficient(cx, &index);
        let value = mul(cx.graph, &[c, d]);
        out = cx.graph.substitute(out, node, value);
    }
    Some(out)
}

/// `u = F(x - c t)`: the ODE for `F` and, when the ODE solver succeeds, the
/// solution in terms of `x - c t`.
pub(super) fn traveling_wave(
    cx: &mut Cx<'_>,
    p: &Problem,
    speed: NodeId,
) -> Option<NodeId> {
    if p.vars.len() != 2 {
        return None;
    }
    let t = p.time(cx.graph);
    let x = 1 - t;
    let mut jets = Jets::new(&p.vars);
    let delta = to_jets(cx, p, &mut jets);
    if !p.constant(cx.graph, delta) {
        return None;
    }
    let s_symbol = cx.graph.interner_mut().fresh_symbol("s");
    let s = cx.graph.symbol_node(s_symbol);
    let f_symbol = cx.graph.interner_mut().fresh_symbol("F");
    let f_node = cx.graph.symbol_node(f_symbol);
    let big_f = cx.graph.node(core::APPLY, &[f_node, s]);
    // u_J = (-c)^(J_t) F^(|J|)
    let ode = jets_to_ode(cx, &jets, delta, big_f, s, |cx, index| {
        let minus = cx.graph.int(-1);
        let neg_c = mul(cx.graph, &[minus, speed]);
        let e = cx.graph.int(i64::from(index[t]));
        cx.graph.node(core::POW, &[neg_c, e])
    })?;
    let ode = cx.simplify(ode);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[ode, zero]);
    let wave = {
        let ct = mul(cx.graph, &[speed, p.vars[t]]);
        sub(cx.graph, p.vars[x], ct)
    };
    match crate::rules::ode::solve_ode(cx, equation, big_f) {
        | Some(answer) => {
            let answer = cx.graph.substitute(answer, s, wave);
            let answer = cx.graph.substitute(answer, big_f, p.unknown);
            let f_at = cx.graph.node(core::APPLY, &[f_node, wave]);
            let answer = cx.graph.replace_subterm(answer, f_at, p.unknown);
            Some(cx.simplify(answer))
        },
        | None => {
            let form = cx.graph.node(core::APPLY, &[f_node, wave]);
            let ansatz = cx.graph.node(core::EQ, &[p.unknown, form]);
            Some(cx.graph.node(core::LIST, &[ansatz, equation]))
        },
    }
}

/// The similarity reduction of a symmetry `ξ = α x + β`, `τ = γ t + δ`,
/// `φ = κ u`: `list(u = g F(z), ODE for F)`, or the solution when the ODE
/// solver closes it.
pub(super) fn similarity(
    cx: &mut Cx<'_>,
    p: &Problem,
    generator: NodeId,
) -> Option<NodeId> {
    if p.vars.len() != 2 || cx.graph.op(generator) != core::LIST {
        return None;
    }
    let t_index = p.time(cx.graph);
    let x_index = 1 - t_index;
    let (x, t) = (p.vars[x_index], p.vars[t_index]);
    let parts = cx.graph.children(generator).to_vec();
    let [g_x, g_t, g_u] = parts.as_slice() else {
        return None;
    };
    let (g_x, g_t) = if x_index == 0 { (*g_x, *g_t) } else { (*g_t, *g_x) };
    let linear = |cx: &mut Cx<'_>, f: NodeId, v: NodeId| -> Option<(f64, f64)> {
        let parts = crate::rules::ode::coefficients_of(cx.graph, f, v)?;
        let value = |cx: &Cx<'_>, n: Option<&NodeId>| n.map_or(Some(0.0), |&n| cx.graph.eval(n, &Env::numeric(0.0)));
        match parts.len() {
            | 1 | 2 => Some((value(cx, parts.get(1))?, value(cx, parts.first())?)),
            | _ => None,
        }
    };
    let (alpha, beta) = linear(cx, g_x, x)?;
    let (gamma, delta_) = linear(cx, g_t, t)?;
    let w_symbol = cx.graph.interner_mut().fresh_symbol("w");
    let w = cx.graph.symbol_node(w_symbol);
    let g_u = cx.graph.replace_subterm(*g_u, p.unknown, w);
    let (kappa, offset) = linear(cx, g_u, w)?;
    if offset.abs() > 1e-12 || (gamma.abs() < 1e-12 && delta_.abs() < 1e-12) {
        return None;
    }
    let num = |cx: &mut Cx<'_>, v: f64| -> Option<NodeId> {
        let r = num_rational::BigRational::from_float(v)?;
        // Small rationals only: the coefficients come from exact algebra.
        let r = r.limit_denominator_hint();
        Some(cx.graph.num(Number::rat(r)))
    };
    let exp = cx.graph.ops().lookup("exp")?;
    let ln = cx.graph.ops().lookup("ln")?;
    // z(x, t), x(z, t) and the gauge g(t).
    let z_symbol = cx.graph.interner_mut().fresh_symbol("z");
    let z = cx.graph.symbol_node(z_symbol);
    let (z_of, x_of, gauge, t0);
    if gamma.abs() > 1e-12 {
        let shift = num(cx, delta_ / gamma)?;
        let tt = add(cx.graph, &[t, shift]);
        let power = num(cx, alpha / gamma)?;
        let scale = cx.graph.node(core::POW, &[tt, power]);
        let k_g = num(cx, kappa / gamma)?;
        gauge = cx.graph.node(core::POW, &[tt, k_g]);
        if alpha.abs() > 1e-12 {
            let xs = num(cx, beta / alpha)?;
            let xx = add(cx.graph, &[x, xs]);
            z_of = div(cx.graph, xx, scale);
            let back = mul(cx.graph, &[z, scale]);
            x_of = sub(cx.graph, back, xs);
        } else {
            let b = num(cx, beta / gamma)?;
            let log = cx.graph.node(ln, &[tt]);
            let drift = mul(cx.graph, &[b, log]);
            z_of = sub(cx.graph, x, drift);
            x_of = add(cx.graph, &[z, drift]);
        }
        let one = cx.graph.int(1);
        t0 = sub(cx.graph, one, shift);
    } else {
        let k_d = num(cx, kappa / delta_)?;
        let kt = mul(cx.graph, &[k_d, t]);
        gauge = cx.graph.node(exp, &[kt]);
        if alpha.abs() > 1e-12 {
            let xs = num(cx, beta / alpha)?;
            let xx = add(cx.graph, &[x, xs]);
            let a_d = num(cx, -alpha / delta_)?;
            let at = mul(cx.graph, &[a_d, t]);
            let decay = cx.graph.node(exp, &[at]);
            z_of = mul(cx.graph, &[xx, decay]);
            let minus = cx.graph.int(-1);
            let grow = cx.graph.node(core::POW, &[decay, minus]);
            let back = mul(cx.graph, &[z, grow]);
            x_of = sub(cx.graph, back, xs);
        } else {
            let c = num(cx, beta / delta_)?;
            let ct = mul(cx.graph, &[c, t]);
            z_of = sub(cx.graph, x, ct);
            x_of = add(cx.graph, &[z, ct]);
        }
        t0 = cx.graph.int(0);
    }
    // u = g(t) F(z(x, t)) through jets f_k = F^(k)(z).
    let mut jets = Jets::new(&p.vars);
    let delta = to_jets(cx, p, &mut jets);
    let f_jets: Vec<NodeId> = (0..6)
        .map(|k| {
            let s = cx.graph.interner_mut().fresh_symbol(&format!("f{k}"));
            cx.graph.symbol_node(s)
        })
        .collect();
    let zx = derivative(cx.graph, z_of, x)?;
    let zt = derivative(cx.graph, z_of, t)?;
    let total = |cx: &mut Cx<'_>, f: NodeId, var: NodeId, dz: NodeId| -> Option<NodeId> {
        let mut terms = vec![derivative(cx.graph, f, var)?];
        for k in 0..5 {
            let partial = derivative(cx.graph, f, f_jets[k])?;
            terms.push(mul(cx.graph, &[f_jets[k + 1], dz, partial]));
        }
        let sum = add(cx.graph, &terms);
        Some(cx.simplify(sum))
    };
    let mut reduced = delta;
    let u_expr = mul(cx.graph, &[gauge, f_jets[0]]);
    for (index, node) in jets.occurring(cx, delta) {
        let mut value = u_expr;
        for _ in 0..index[x_index] {
            value = total(cx, value, x, zx)?;
        }
        for _ in 0..index[t_index] {
            value = total(cx, value, t, zt)?;
        }
        reduced = cx.graph.substitute(reduced, node, value);
    }
    let reduced = cx.graph.substitute(reduced, x, x_of);
    let reduced = cx.simplify(reduced);
    // Must be t-independent up to a factor: compare two times numerically.
    let at = |cx: &mut Cx<'_>, time: NodeId| -> NodeId {
        let e = cx.graph.substitute(reduced, t, time);
        cx.simplify(e)
    };
    let one = cx.graph.int(1);
    let t1 = add(cx.graph, &[t0, one]);
    let two = cx.graph.int(2);
    let t2 = add(cx.graph, &[t0, two]);
    let (e1, e2) = (at(cx, t1), at(cx, t2));
    let mut ratio = None;
    for point in [0.3, 0.8, 1.7] {
        let mut env = Env::numeric(0.0);
        env.bind(z_symbol, point);
        for (k, &f) in f_jets.iter().enumerate() {
            env.bind(cx.graph.symbol_of(f)?, 0.4 + 0.3 * f64::from(u32::try_from(k).ok()?) + point);
        }
        let (a, b) = (cx.graph.eval(e1, &env)?, cx.graph.eval(e2, &env)?);
        if a.abs() < 1e-12 {
            continue;
        }
        let r = b / a;
        if ratio.is_some_and(|q: f64| (q - r).abs() > 1e-8 * q.abs().max(1.0)) {
            return None;
        }
        ratio = Some(r);
    }
    ratio?;
    // The ODE in z with f_k → F^(k)(z).
    let f_symbol = cx.graph.interner_mut().fresh_symbol("F");
    let f_node = cx.graph.symbol_node(f_symbol);
    let big_f = cx.graph.node(core::APPLY, &[f_node, z]);
    let diff = cx.graph.ops().lookup("diff")?;
    let mut ode = e1;
    let mut d = big_f;
    for &f in &f_jets {
        ode = cx.graph.substitute(ode, f, d);
        d = cx.graph.node(diff, &[d, z]);
    }
    let ode = cx.simplify(ode);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[ode, zero]);
    let form = {
        let f_of_z = cx.graph.node(core::APPLY, &[f_node, z_of]);
        let product = mul(cx.graph, &[gauge, f_of_z]);
        cx.simplify(product)
    };
    if let Some(answer) = crate::rules::ode::solve_ode(cx, equation, big_f) {
        if let [_, rhs] = *cx.graph.children(answer) {
            let solved = mul(cx.graph, &[gauge, rhs]);
            let solved = cx.graph.substitute(solved, z, z_of);
            let solved = cx.simplify(solved);
            return Some(cx.graph.node(core::EQ, &[p.unknown, solved]));
        }
    }
    let ansatz = cx.graph.node(core::EQ, &[p.unknown, form]);
    Some(cx.graph.node(core::LIST, &[ansatz, equation]))
}

/// Noether's conserved current `list(J^x, J^t)` of a Lagrangian density
/// `L(x, t, u, u_x, u_t)` and a symmetry `list(ξ, τ, φ)`.
pub(super) fn noether_current(
    cx: &mut Cx<'_>,
    lagrangian: NodeId,
    unknown: NodeId,
    generator: NodeId,
) -> Option<NodeId> {
    let vars = cx.graph.children(unknown).get(1..)?.to_vec();
    if vars.len() != 2 || cx.graph.op(generator) != core::LIST {
        return None;
    }
    let parts = cx.graph.children(generator).to_vec();
    let [xi, tau, phi] = parts.as_slice() else {
        return None;
    };
    let diff = cx.graph.ops().lookup("diff")?;
    let ux = cx.graph.node(diff, &[unknown, vars[0]]);
    let ut = cx.graph.node(diff, &[unknown, vars[1]]);
    // ∂L/∂u_x and ∂L/∂u_t through stand-in symbols.
    let (sx, st) = {
        let a = cx.graph.interner_mut().fresh_symbol("px");
        let b = cx.graph.interner_mut().fresh_symbol("pt");
        (cx.graph.symbol_node(a), cx.graph.symbol_node(b))
    };
    let plain = cx.graph.replace_subterm(lagrangian, ux, sx);
    let plain = cx.graph.replace_subterm(plain, ut, st);
    let dl_dux = derivative(cx.graph, plain, sx)?;
    let dl_dut = derivative(cx.graph, plain, st)?;
    let back = |cx: &mut Cx<'_>, e: NodeId| {
        let e = cx.graph.substitute(e, sx, ux);
        cx.graph.substitute(e, st, ut)
    };
    let (dl_dux, dl_dut) = (back(cx, dl_dux), back(cx, dl_dut));
    // Q = φ - ξ u_x - τ u_t
    let minus = cx.graph.int(-1);
    let q = {
        let a = mul(cx.graph, &[minus, *xi, ux]);
        let b = mul(cx.graph, &[minus, *tau, ut]);
        add(cx.graph, &[*phi, a, b])
    };
    let jx = {
        let a = mul(cx.graph, &[*xi, lagrangian]);
        let b = mul(cx.graph, &[q, dl_dux]);
        add(cx.graph, &[a, b])
    };
    let jt = {
        let a = mul(cx.graph, &[*tau, lagrangian]);
        let b = mul(cx.graph, &[q, dl_dut]);
        add(cx.graph, &[a, b])
    };
    let (jx, jt) = (cx.simplify(jx), cx.simplify(jt));
    Some(cx.graph.node(core::LIST, &[jx, jt]))
}

/// Rational approximations with small denominators.
trait LimitDenominator {
    fn limit_denominator_hint(self) -> Self;
}

impl LimitDenominator for num_rational::BigRational {
    fn limit_denominator_hint(self) -> Self {
        use num_traits::ToPrimitive;
        let value = self.to_f64().unwrap_or(0.0);
        for d in 1..=64_i64 {
            let n = (value * d as f64).round();
            if (n / d as f64 - value).abs() < 1e-12 {
                return Self::new((n as i64).into(), d.into());
            }
        }
        self
    }
}
