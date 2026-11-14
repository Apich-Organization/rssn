//! General solution frameworks that are not tied to a named equation.
//!
//! * **Complete integrals of first-order equations** `F(x, y, u, p, q) = 0`
//!   (Lagrange–Charpit): the classes where Charpit's system has an obvious
//!   first integral — `F(p, q) = 0` (`u = a x + b(a) y + c`), Clairaut's
//!   `u = x p + y q + f(p, q)`, separable `F₁(x, p) = F₂(y, q)`
//!   (`p = g₁(x, a)`, `q = g₂(y, a)`) and `F(u, p, q) = 0` (`u = U(x + a y)`,
//!   an ODE for `U`).
//! * **Separation of variables** for any (also non-linear) equation in
//!   `u(x, t)`: the multiplicative ansatz `u = X(x) T(t)` and the additive
//!   ansatz `u = X(x) + T(t)`. The equation (divided by `X T` in the
//!   multiplicative case) is separable when its mixed derivatives in the
//!   `x`-side and `t`-side jets vanish; then each side equals `∓λ` and the
//!   two ODEs go to the ODE solver.
//! * **Fourier transform in space** for linear constant-coefficient
//!   evolution equations `u_t = Σ aₙ ∂ₓⁿ u` of any order on the line:
//!   `u = (1/2π) ∫∫ g(s) exp(i k (x - s) + P(i k) t) ds dk`.

use super::Conditions;
use super::Problem;
use super::dummy;
use super::symmetry::Jets;
use super::symmetry::to_jets;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

fn fresh(
    cx: &mut Cx<'_>,
    stem: &str,
) -> NodeId {
    let s = cx.graph.interner_mut().fresh_symbol(stem);
    cx.graph.symbol_node(s)
}

fn depends(
    cx: &Cx<'_>,
    e: NodeId,
    v: NodeId,
) -> bool {
    cx.graph
        .symbol_of(v)
        .is_some_and(|s| cx.graph.depends_on(cx.graph.find(e), s))
}

/// A complete integral `u = …` with constants `a`, `b` (named `A`, `B`).
pub(super) fn complete_integral(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<NodeId> {
    if p.vars.len() != 2 || p.order() != 1 {
        return None;
    }
    let (x, y) = (p.vars[0], p.vars[1]);
    let mut jets = Jets::new(&p.vars);
    let f = to_jets(cx, p, &mut jets);
    let f = cx.simplify(f);
    let u = jets.get(cx, &vec![0, 0]);
    let px = jets.get(cx, &vec![1, 0]);
    let qy = jets.get(cx, &vec![0, 1]);
    let a = cx.graph.sym("A");
    let b = cx.graph.sym("B");
    let c = cx.graph.sym("C");
    let (dx, dy, du) = (depends(cx, f, x), depends(cx, f, y), depends(cx, f, u));
    // F(p, q) = 0: u = a x + b y + c with F(a, b) = 0.
    if !dx && !dy && !du {
        let g = cx.graph.substitute(f, px, a);
        let g = cx.graph.substitute(g, qy, b);
        let bs = crate::rules::ode::solve_in(cx.graph, g, b)?;
        let b_of_a = *bs.first()?;
        let ax = mul(cx.graph, &[a, x]);
        let by = mul(cx.graph, &[b_of_a, y]);
        let sol = add(cx.graph, &[ax, by, c]);
        return Some(cx.simplify(sol));
    }
    // Clairaut: F = u - x p - y q - f(p, q).
    let clairaut = {
        let xp = mul(cx.graph, &[x, px]);
        let yq = mul(cx.graph, &[y, qy]);
        let rest = sub(cx.graph, u, xp);
        let rest = sub(cx.graph, rest, yq);
        // f = u - x p - y q - F·k for the normalisation with u's coefficient.
        let du_coeff = derivative(cx.graph, f, u)?;
        let du_coeff = cx.simplify(du_coeff);
        cx.graph
            .number_of(du_coeff)
            .filter(|n| !n.is_zero())
            .cloned()
            .map(|k| {
                let k = cx.graph.num(k);
                let scaled = {
                    let inverse = powi(cx.graph, k, -1);
                    mul(cx.graph, &[f, inverse])
                };
                sub(cx.graph, scaled, rest)
            })
    };
    if let Some(remainder) = clairaut {
        let remainder = cx.simplify(remainder);
        if !depends(cx, remainder, x) && !depends(cx, remainder, y) && !depends(cx, remainder, u) {
            // F/k = u - x p - y q + r(p, q): u = a x + b y - r(a, b).
            let r = cx.graph.substitute(remainder, px, a);
            let r = cx.graph.substitute(r, qy, b);
            let ax = mul(cx.graph, &[a, x]);
            let by = mul(cx.graph, &[b, y]);
            let neg_r = neg(cx.graph, r);
            let sol = add(cx.graph, &[ax, by, neg_r]);
            return Some(cx.simplify(sol));
        }
    }
    // Separable F₁(x, p) = F₂(y, q): mixed derivatives vanish.
    if !du {
        let mixed = [(x, y), (x, qy), (px, y), (px, qy)];
        let mut separable = true;
        for (s1, s2) in mixed {
            let d = derivative(cx.graph, f, s1).and_then(|d| derivative(cx.graph, d, s2));
            match d.map(|d| cx.simplify(d)) {
                | Some(d) if cx.is_zero(d) => {},
                | _ => {
                    separable = false;
                    break;
                },
            }
        }
        if separable {
            // F = F₁(x, p) + G(y, q): F₁ = a, G = -a.
            let zero = cx.graph.int(0);
            let y0 = cx.graph.substitute(f, y, zero);
            let q0 = cx.graph.substitute(y0, qy, zero);
            let f1 = cx.simplify(q0);
            let f1_eq = sub(cx.graph, f1, a);
            let g1 = *crate::rules::ode::solve_in(cx.graph, f1_eq, px)?.first()?;
            let rest = sub(cx.graph, f, f1);
            let rest = cx.simplify(rest);
            let f2_eq = add(cx.graph, &[rest, a]);
            let g2 = *crate::rules::ode::solve_in(cx.graph, f2_eq, qy)?.first()?;
            let i1 = crate::rules::calculus::antiderivative(cx, g1, x)?;
            let i2 = crate::rules::calculus::antiderivative(cx, g2, y)?;
            let sol = add(cx.graph, &[i1, i2, b]);
            return Some(cx.simplify(sol));
        }
    }
    // F(u, p, q) = 0: u = U(z), z = x + a y: F(U, U', a U') = 0.
    if !dx && !dy {
        let z = fresh(cx, "z");
        let f_sym = fresh(cx, "U");
        let big_u = cx.graph.node(core::APPLY, &[f_sym, z]);
        let diff = cx.graph.ops().lookup("diff")?;
        let du_dz = cx.graph.node(diff, &[big_u, z]);
        let a_du = mul(cx.graph, &[a, du_dz]);
        let ode = cx.graph.substitute(f, qy, a_du);
        let ode = cx.graph.substitute(ode, px, du_dz);
        let ode = cx.graph.substitute(ode, u, big_u);
        let zero = cx.graph.int(0);
        let equation = cx.graph.node(core::EQ, &[ode, zero]);
        let answer = crate::rules::ode::solve_ode(cx, equation, big_u)?;
        let ay = mul(cx.graph, &[a, y]);
        let wave = add(cx.graph, &[x, ay]);
        let answer = cx.graph.substitute(answer, z, wave);
        let at = cx.graph.node(core::APPLY, &[f_sym, wave]);
        let answer = cx.graph.replace_subterm(answer, at, p.unknown);
        return match *cx.graph.children(answer) {
            | [lhs, rhs] if lhs == p.unknown => Some(cx.simplify(rhs)),
            | _ => None,
        };
    }
    None
}

/// Separated solutions `u = X(x) T(t)` or `u = X(x) + T(t)`.
pub(super) fn separate(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<NodeId> {
    if p.vars.len() != 2 {
        return None;
    }
    let (x, t) = (p.vars[0], p.vars[1]);
    let mut jets = Jets::new(&p.vars);
    let delta = to_jets(cx, p, &mut jets);
    let order = p.order() as usize;
    // Jets of X and T as plain symbols: X_k = X^(k)(x), T_k = T^(k)(t).
    let xs: Vec<NodeId> = (0..=order).map(|k| fresh(cx, &format!("X{k}"))).collect();
    let ts: Vec<NodeId> = (0..=order).map(|k| fresh(cx, &format!("T{k}"))).collect();
    let lambda = cx.graph.sym("lambda");
    for multiplicative in [true, false] {
        let mut e = delta;
        for (index, node) in jets.occurring(cx, delta) {
            let (i, j) = (index[0] as usize, index[1] as usize);
            let value = if multiplicative {
                mul(cx.graph, &[xs[i], ts[j]])
            } else if i == 0 && j == 0 {
                add(cx.graph, &[xs[0], ts[0]])
            } else if j == 0 {
                xs[i]
            } else if i == 0 {
                ts[j]
            } else {
                cx.graph.int(0)
            };
            e = cx.graph.substitute(e, node, value);
        }
        if multiplicative {
            let product = mul(cx.graph, &[xs[0], ts[0]]);
            let inverse = powi(cx.graph, product, -1);
            e = mul(cx.graph, &[e, inverse]);
        }
        let e = cx.simplify(e);
        // Separability: every mixed second derivative vanishes.
        let mut left: Vec<NodeId> = xs.clone();
        left.push(x);
        let mut right: Vec<NodeId> = ts.clone();
        right.push(t);
        let separable = left.iter().all(|&l| {
            right.iter().all(|&r| {
                derivative(cx.graph, e, l)
                    .and_then(|d| derivative(cx.graph, d, r))
                    .map(|d| cx.simplify(d))
                    .is_some_and(|d| cx.is_zero(d))
            })
        });
        if !separable {
            continue;
        }
        // e = A(x, X…) + B(t, T…): A from setting the T-side to a generic
        // point, B = e - A; A = -λ, B = λ.
        let mut a_side = e;
        for (k, &tk) in ts.iter().enumerate() {
            let point = cx.graph.int(i64::try_from(k).ok()? + 1);
            a_side = cx.graph.substitute(a_side, tk, point);
        }
        let one = cx.graph.int(1);
        a_side = cx.graph.substitute(a_side, t, one);
        let a_side = cx.simplify(a_side);
        let b_side = sub(cx.graph, e, a_side);
        let b_side = cx.simplify(b_side);
        let solve_side = |cx: &mut Cx<'_>,
                          side: NodeId,
                          var: NodeId,
                          js: &[NodeId],
                          sign: i64|
         -> Option<NodeId> {
            let f = fresh(cx, "F");
            let big = cx.graph.node(core::APPLY, &[f, var]);
            let diff = cx.graph.ops().lookup("diff")?;
            let mut chain = vec![big];
            for _ in 1..js.len() {
                let last = *chain.last()?;
                chain.push(cx.graph.node(diff, &[last, var]));
            }
            let mut ode = side;
            for (k, &jk) in js.iter().enumerate() {
                ode = cx.graph.substitute(ode, jk, chain[k]);
            }
            let s = cx.graph.int(sign);
            let shifted = mul(cx.graph, &[s, lambda]);
            let ode = add(cx.graph, &[ode, shifted]);
            let zero = cx.graph.int(0);
            let equation = cx.graph.node(core::EQ, &[ode, zero]);
            let answer = crate::rules::ode::solve_ode(cx, equation, big)?;
            match *cx.graph.children(answer) {
                | [lhs, rhs] if lhs == big => Some(rhs),
                | _ => None,
            }
        };
        let Some(xsol) = solve_side(cx, a_side, x, &xs, 1) else {
            continue;
        };
        let Some(tsol) = solve_side(cx, b_side, t, &ts, -1) else {
            continue;
        };
        // Distinct constants on the two sides.
        let tsol = rename_constants(cx, tsol, "D");
        let combined = if multiplicative {
            mul(cx.graph, &[xsol, tsol])
        } else {
            add(cx.graph, &[xsol, tsol])
        };
        return Some(cx.simplify(combined));
    }
    None
}

fn rename_constants(
    cx: &mut Cx<'_>,
    e: NodeId,
    stem: &str,
) -> NodeId {
    let mut out = e;
    for k in 1..=6 {
        let c = cx.graph.sym(&format!("C{k}"));
        let d = cx.graph.sym(&format!("{stem}{k}"));
        out = cx.graph.substitute(out, c, d);
    }
    out
}

/// Homogeneous constant-coefficient evolution equations on whole space
/// (`a₂ u_tt + a₁ u_t + Σ c_α ∂^α u = 0`, any order in `x`, one to three
/// space dimensions) by the Fourier transform: with the symbol
/// `P(ik) = Σ c_α Π (i k_j)^{α_j}` the transform solves the ordinary
/// equation `a₂ û_tt + a₁ û_t + P û = 0`,
///
/// * first order: `û = û₀ e^{-P t/a₁}`;
/// * second order: `û = e^{-a t}[û₀ (cosh β t + (a/β) sinh β t) + û₁ sinh(β t)/β]`,
///   `a = a₁/(2 a₂)`, `β = √(a² - P/a₂)`,
///
/// and `u = (2π)^{-d} ∫ û e^{i k·x} dk` with `û₀(k) = ∫ g(s) e^{-i k·s} ds`.
/// The answer is that integral representation (the inversion table of the
/// transform module closes it only for special data).
pub(super) fn fourier_evolution(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let n = p.dimension();
    let time = p.time(cx.graph);
    let space: Vec<usize> = (0..n).filter(|&j| j != time).collect();
    if space.is_empty() || space.len() > 3 {
        return None;
    }
    let (a1, a2) = (
        p.coefficient(cx.graph, &p.unit(time, 1)),
        p.coefficient(cx.graph, &p.unit(time, 2)),
    );
    if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
        return None;
    }
    let order = if !cx.is_zero(a2) {
        2
    } else if !cx.is_zero(a1) {
        1
    } else {
        return None;
    };
    let i_unit = cx
        .graph
        .ops()
        .lookup("I")
        .map(|op| cx.graph.node(op, &[]))?;
    let k_stems = ["k", "kk", "kkk"];
    let ks: Vec<NodeId> = (0..space.len())
        .map(|i| dummy(cx, p, k_stems.get(i).copied().unwrap_or("k")).0)
        .collect();
    let mut symbol_terms = Vec::new();
    for (index, c) in p.linear.clone() {
        if cx.is_zero(c) || index == p.unit(time, 1) || index == p.unit(time, 2) {
            continue;
        }
        if index[time] != 0 || !p.constant(cx.graph, c) {
            return None;
        }
        let mut factors = vec![c];
        for (pos, &j) in space.iter().enumerate() {
            let e = i64::from(index[j]);
            if e > 0 {
                let ik = mul(cx.graph, &[i_unit, ks[pos]]);
                factors.push(powi(cx.graph, ik, e));
            }
        }
        symbol_terms.push(mul(cx.graph, &factors));
    }
    if symbol_terms.is_empty() {
        return None;
    }
    let symbol = add(cx.graph, &symbol_terms);
    let symbol = cx.simplify(symbol);
    let zero_index = vec![0; n];
    let displacement = conditions.find(cx.graph, time, &zero_index, None)?;
    let velocity = if order == 2 {
        conditions.find(cx.graph, time, &p.unit(time, 1), None)
    } else {
        None
    };
    if conditions.0.len() != 1 + usize::from(velocity.is_some())
        || !cx
            .graph
            .number_of(displacement.point)
            .is_some_and(Number::is_zero)
    {
        return None;
    }
    if velocity.is_some_and(|v| !cx.graph.number_of(v.point).is_some_and(Number::is_zero)) {
        return None;
    }
    let t = p.vars[time];
    let exp = |cx: &mut Cx<'_>, e: NodeId| -> Option<NodeId> { super::util::call(cx, "exp", &[e]) };
    // The multipliers of û₀ and û₁.
    let multipliers = |cx: &mut Cx<'_>, t: NodeId| -> Option<(NodeId, Option<NodeId>)> {
        Some(if order == 1 {
            let ratio = {
                let q = super::util::div(cx, symbol, a1);
                neg(cx.graph, q)
            };
            let e = mul(cx.graph, &[ratio, t]);
            let m0 = exp(cx, e)?;
            (cx.simplify(m0), None)
        } else {
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let alpha = {
                let q = super::util::div(cx, a1, a2);
                let h = mul(cx.graph, &[half, q]);
                cx.simplify(h)
            };
            let beta2 = {
                let a_sq = powi(cx.graph, alpha, 2);
                let q = super::util::div(cx, symbol, a2);
                let nq = neg(cx.graph, q);
                add(cx.graph, &[a_sq, nq])
            };
            let beta2 = cx.simplify(beta2);
            let beta = super::util::sqrt(cx, beta2)?;
            let bt = mul(cx.graph, &[beta, t]);
            let (ch, sh) = (
                super::util::call(cx, "cosh", &[bt])?,
                super::util::call(cx, "sinh", &[bt])?,
            );
            let damp = {
                let m = mul(cx.graph, &[alpha, t]);
                let nm = neg(cx.graph, m);
                exp(cx, nm)?
            };
            let ratio = super::util::div(cx, alpha, beta);
            let lead = {
                let rs = mul(cx.graph, &[ratio, sh]);
                add(cx.graph, &[ch, rs])
            };
            let m0 = mul(cx.graph, &[damp, lead]);
            let m1 = {
                let q = super::util::div(cx, sh, beta);
                mul(cx.graph, &[damp, q])
            };
            (cx.simplify(m0), Some(cx.simplify(m1)))
        })
    };
    let (m0, m1) = multipliers(cx, t)?;
    let (defint, infinity, pi) = (
        cx.graph.ops().lookup("defint")?,
        cx.graph.ops().lookup("oo")?,
        cx.graph.ops().lookup("pi")?,
    );
    let oo = cx.graph.node(infinity, &[]);
    let minus_oo = neg(cx.graph, oo);
    let s_stems = ["s", "ss", "sss"];
    let ss: Vec<NodeId> = (0..space.len())
        .map(|i| dummy(cx, p, s_stems.get(i).copied().unwrap_or("s")).0)
        .collect();
    let mut parts = Vec::new();
    // The transform of a datum in closed form, by the table of the
    // transform module (`None` when it stays inert).
    let transform = |cx: &mut Cx<'_>, value: NodeId| -> Option<NodeId> {
        let mut f = value;
        for (pos, &j) in space.iter().enumerate() {
            let node = super::util::call(cx, "fourier", &[f, p.vars[j], ks[pos]])?;
            f = cx.simplify(node);
            if super::util::contains_op(cx.graph, f, "fourier") {
                return None;
            }
        }
        Some(f)
    };
    let mut spectrum = Vec::new();
    for (data, multiplier) in [
        (Some(displacement.value), Some(m0)),
        (velocity.map(|v| v.value), m1),
    ] {
        let (Some(value), Some(multiplier)) = (data, multiplier) else {
            continue;
        };
        if cx.is_zero(value) {
            continue;
        }
        if let Some(transformed) = transform(cx, value) {
            spectrum.push(mul(cx.graph, &[transformed, multiplier]));
            continue;
        }
        let mut g = value;
        for (pos, &j) in space.iter().enumerate() {
            g = cx.graph.substitute(g, p.vars[j], ss[pos]);
        }
        // exp(i k·(x - s)) over every axis.
        let mut phases = Vec::new();
        for (pos, &j) in space.iter().enumerate() {
            let gap = sub(cx.graph, p.vars[j], ss[pos]);
            let kg = mul(cx.graph, &[i_unit, ks[pos], gap]);
            phases.push(kg);
        }
        let phase = add(cx.graph, &phases);
        let kernel = exp(cx, phase)?;
        let mut body = mul(cx.graph, &[g, kernel, multiplier]);
        for &s_var in ss.iter().rev() {
            body = cx.graph.node(defint, &[body, s_var, minus_oo, oo]);
        }
        for &k_var in ks.iter().rev() {
            body = cx.graph.node(defint, &[body, k_var, minus_oo, oo]);
        }
        parts.push(body);
    }
    // Duhamel: ∫_0^t G(t - τ) f̂(k, τ) dτ for the source `f = -source`.
    if !p.homogeneous(cx.graph) {
        let forcing = neg(cx.graph, p.source);
        let forcing = cx.simplify(forcing);
        if !cx.is_zero(forcing) {
            let (tau, _) = dummy(cx, p, "tau");
            let f_tau = cx.graph.substitute(forcing, t, tau);
            let transformed = transform(cx, f_tau)?;
            let elapsed = sub(cx.graph, t, tau);
            let (g0, g1) = multipliers(cx, elapsed)?;
            let green = if order == 1 {
                super::util::div(cx, g0, a1)
            } else {
                super::util::div(cx, g1?, a2)
            };
            let body = mul(cx.graph, &[green, transformed]);
            let zero = cx.graph.int(0);
            let integral = cx.graph.node(defint, &[body, tau, zero, t]);
            spectrum.push(cx.simplify(integral));
        }
    }
    let mut closed = Vec::new();
    if !spectrum.is_empty() {
        let total_spectrum = add(cx.graph, &spectrum);
        let total_spectrum = cx.simplify(total_spectrum);
        // Inversion by the table.
        let mut inverted = total_spectrum;
        let mut by_table = true;
        for (pos, &j) in space.iter().enumerate() {
            let node = super::util::call(cx, "inverse_fourier", &[inverted, ks[pos], p.vars[j]])?;
            inverted = cx.simplify(node);
            if super::util::contains_op(cx.graph, inverted, "inverse_fourier") {
                by_table = false;
                break;
            }
        }
        if by_table {
            closed.push(inverted);
        } else {
            let mut phases = Vec::new();
            for (pos, &j) in space.iter().enumerate() {
                phases.push(mul(cx.graph, &[i_unit, ks[pos], p.vars[j]]));
            }
            let phase = add(cx.graph, &phases);
            let kernel = exp(cx, phase)?;
            let mut body = mul(cx.graph, &[total_spectrum, kernel]);
            for &k_var in ks.iter().rev() {
                body = cx.graph.node(defint, &[body, k_var, minus_oo, oo]);
            }
            parts.push(body);
        }
    }
    let mut pieces = Vec::new();
    if !parts.is_empty() {
        let total = add(cx.graph, &parts);
        let pi = cx.graph.node(pi, &[]);
        let two = cx.graph.int(2);
        let scale = mul(cx.graph, &[two, pi]);
        let inverse = powi(cx.graph, scale, -i64::try_from(space.len()).ok()?);
        pieces.push(mul(cx.graph, &[inverse, total]));
    }
    pieces.extend(closed);
    if pieces.is_empty() {
        return None;
    }
    let answer = add(cx.graph, &pieces);
    Some(answer)
}
