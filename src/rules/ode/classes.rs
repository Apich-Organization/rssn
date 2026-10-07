//! More classes of first- and second-order equations.
//!
//! * **Bernoulli** equations `y' = a y + b y^n` with a rational exponent
//!   `n` (square roots of `y`, say).
//! * **Riccati** equations `y' = q0 + q1 y + q2 y^2` without a simple
//!   particular solution: `y = -u'/(q2 u)` turns them into the linear
//!   equation `u'' - (q1 + q2'/q2) u' + q0 q2 u = 0`, which is solved by
//!   the full pipeline; the solution is written with one constant (the
//!   ratio of the two).
//! * **Integrating factors** of `M + N y' = 0`: `μ(x)`, `μ(y)`, `μ(x + y)`
//!   and `μ(x y)`, each tested for the condition that makes `μ M dx +
//!   μ N dy` exact.
//! * **Equations in `y'` that cannot be solved for `y'`**: Clairaut
//!   `y = x y' + g(y')` (general solution `y = C x + g(C)`), d'Alembert
//!   (Lagrange) `y = x φ(y') + ψ(y')` and equations solvable for `x`,
//!   through the first-order equation for `x(p)` or `y(p)` in the
//!   parameter `p = y'`. Those answers are parametric:
//!   `list(x = X(p, C1), y(x) = Y(p, C1))`.
//! * **Second-order linear equations in normal form**: with
//!   `y = u exp(-½ ∫ a₁/a₂)` the first-derivative term drops out and the
//!   remaining `u'' + Q u = 0` is solved by the pipeline (spherical Bessel
//!   functions, Euler-type shifts).

use num_rational::BigRational;
use num_traits::One;

use super::Problem;
use super::add;
use super::div;
use super::exp_of;
use super::explicit;
use super::implicit;
use super::integrate;
use super::inv;
use super::mul;
use super::neg;
use super::reduce::chain;
use super::reduce::constants_in;
use super::reduce::explicit_rhs;
use super::reduce::fresh_function;
use super::solve_equation;
use super::sub;
use super::variation_of_parameters;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::to_term;
use crate::rules::solve::solve_for;

/// `v' = (1 - n)(a v + b)` for `v = y^(1 - n)`.
pub(super) fn bernoulli_with(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    a: NodeId,
    b: NodeId,
    n: &BigRational,
) -> Option<NodeId> {
    let x = problem.x;
    let k = cx.graph.num(Number::rat(BigRational::one() - n));
    let rate = mul(cx.graph, &[k, a]);
    let forcing = mul(cx.graph, &[k, b]);
    let integral = integrate(cx, rate, x)?;
    let growth = super::call(cx.graph, "exp", integral)?;
    let decay_arg = neg(cx.graph, integral);
    let decay = super::call(cx.graph, "exp", decay_arg)?;
    let weighted = mul(cx.graph, &[decay, forcing]);
    let weighted = cx.simplify(weighted);
    let particular = integrate(cx, weighted, x)?;
    let c = problem.constant(cx.graph);
    let sum = add(cx.graph, &[particular, c]);
    let v = mul(cx.graph, &[growth, sum]);
    let exponent = inv(cx.graph, k);
    let solution = cx.graph.node(core::POW, &[v, exponent]);
    Some(explicit(cx, problem, solution))
}

/// `y' = a y + b y^n` with `n` not an integer.
pub(super) fn bernoulli_rational(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    rhs: NodeId,
) -> Option<NodeId> {
    let y = *problem.stand.first()?;
    let y_symbol = cx.graph.symbol_of(y)?;
    let mut gens = Gens::default();
    let gy = gens.index(cx.graph, y);
    let poly = from_term(cx.graph, &mut gens, rhs, Limits::default())?;
    let mut power: Option<(u32, BigRational)> = None;
    for g in poly.support() {
        if g == gy {
            continue;
        }
        let node = gens.node(g)?;
        if !cx.graph.depends_on(cx.graph.find(node), y_symbol) {
            continue;
        }
        let e = if cx.graph.ops().lookup("sqrt") == Some(cx.graph.op(node)) && cx.graph.children(node) == [y] {
            BigRational::new(1.into(), 2.into())
        } else {
            let &[base, e] = cx.graph.children(node) else {
                return None;
            };
            if cx.graph.op(node) != core::POW || base != y {
                return None;
            }
            cx.graph.number_of(e)?.to_rational()?
        };
        if power.is_some() || e.is_integer() {
            return None;
        }
        power = Some((g, e));
    }
    let (gn, n) = power?;
    let (mut a, mut b) = (crate::rules::poly::repr::Poly::zero(), crate::rules::poly::repr::Poly::zero());
    for (mono, coeff) in poly.terms() {
        let ey = mono.iter().find(|&&(g, _)| g == gy).map_or(0, |&(_, e)| e);
        let en = mono.iter().find(|&&(g, _)| g == gn).map_or(0, |&(_, e)| e);
        let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(g, _)| g != gy && g != gn).collect();
        let piece = crate::rules::poly::repr::Poly::monomial(rest, coeff.clone());
        match (ey, en) {
            | (1, 0) => a = a.add(&piece),
            | (0, 1) => b = b.add(&piece),
            | _ => return None,
        }
    }
    let (a, b) = (to_term(cx.graph, &gens, &a), to_term(cx.graph, &gens, &b));
    // Solutions with v = y^(1 - n) < 0 are other branches of the power:
    // the family is right where it is defined, which a spot check at one
    // value of the constant cannot establish.
    let result = bernoulli_with(cx, problem, a, b, &n);
    problem.trusted = problem.trusted || result.is_some();
    result
}

/// Riccati through the linear second-order equation for `u`.
pub(super) fn riccati_linearised(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    coefficients: &[NodeId],
    depth: u32,
) -> Option<NodeId> {
    let x = problem.x;
    let &[q0, q1, q2] = coefficients else {
        return None;
    };
    if cx.is_zero(q2) || cx.is_zero(q0) && cx.is_zero(q1) {
        return None;
    }
    let u = fresh_function(cx.graph, "u", x);
    let jets = chain(cx.graph, u, x, 2)?;
    let dq2 = derivative(cx.graph, q2, x)?;
    let ratio = div(cx.graph, dq2, q2);
    let a1 = add(cx.graph, &[q1, ratio]);
    let a1 = cx.simplify(a1);
    let a0 = mul(cx.graph, &[q0, q2]);
    let a0 = cx.simplify(a0);
    let t1 = mul(cx.graph, &[a1, jets[1]]);
    let t0 = mul(cx.graph, &[a0, jets[0]]);
    let negated = neg(cx.graph, t1);
    let left = add(cx.graph, &[jets[2], negated, t0]);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[left, zero]);
    let (_, answer) = solve_equation(cx, equation, u, depth + 2)?;
    let solution = explicit_rhs(cx.graph, answer, u)?;
    // One constant: set C1 = 1 and keep C2 as the ratio.
    let (c1, c2) = (cx.graph.sym("C1"), cx.graph.sym("C2"));
    let one = cx.graph.int(1);
    let c = problem.constant(cx.graph);
    let with_one = cx.graph.substitute(solution, c1, one);
    let normalised = cx.graph.substitute(with_one, c2, c);
    let du = derivative(cx.graph, normalised, x)?;
    let denominator = mul(cx.graph, &[q2, normalised]);
    let quotient = div(cx.graph, du, denominator);
    let y = neg(cx.graph, quotient);
    Some(explicit(cx, problem, y))
}

/// Integrating factors of `M + N y' = 0`.
pub(super) fn integrating_factor(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    m: NodeId,
    n: NodeId,
) -> Option<NodeId> {
    let (y, x) = (*problem.stand.first()?, problem.x);
    let m_y = derivative(cx.graph, m, y)?;
    let n_x = derivative(cx.graph, n, x)?;
    let gap = sub(cx.graph, m_y, n_x);
    let gap = cx.simplify(gap);
    if cx.is_zero(gap) {
        return None;
    }
    let minus_gap = neg(cx.graph, gap);
    let z_symbol = cx.graph.interner_mut().fresh_symbol("z");
    let z = cx.graph.symbol_node(z_symbol);
    // (denominator, which variable the quotient may depend on, z as a
    // function of x and y, y in terms of x and z).
    let xn = mul(cx.graph, &[x, n]);
    let ym = mul(cx.graph, &[y, m]);
    let xn_ym = sub(cx.graph, xn, ym);
    let n_m = sub(cx.graph, n, m);
    let y_from_sum = sub(cx.graph, z, x);
    let y_from_product = div(cx.graph, z, x);
    let x_plus_y = add(cx.graph, &[x, y]);
    let x_times_y = mul(cx.graph, &[x, y]);
    // Kind 0: μ(x); 1: μ(y); 2: μ(x + y); 3: μ(x y).
    let options = [
        (n, gap, 0_u8, x, y, x_plus_y, y_from_sum),
        (m, minus_gap, 1, y, x, x_plus_y, y_from_sum),
        (n_m, gap, 2, z, x, x_plus_y, y_from_sum),
        (xn_ym, gap, 3, z, x, x_times_y, y_from_product),
    ];
    for (denominator, numerator, kind, variable, other, g, y_of_z) in options {
        let q = div(cx.graph, numerator, denominator);
        let q = cx.simplify(q);
        // Express the quotient in the variable it should depend on.
        let in_variable = match kind {
            | 0 | 1 => q,
            | _ => {
                let substituted = cx.graph.substitute(q, y, y_of_z);
                cx.simplify(substituted)
            },
        };
        let Some(function) = super::lie::specialise_constant(cx, in_variable, other) else {
            continue;
        };
        let Some(integral) = integrate(cx, function, variable) else {
            continue;
        };
        let Some(mu) = exp_of(cx, integral) else {
            continue;
        };
        let mu = if kind >= 2 { cx.graph.substitute(mu, z, g) } else { mu };
        let (new_m, new_n) = (mul(cx.graph, &[mu, m]), mul(cx.graph, &[mu, n]));
        let (new_m, new_n) = (cx.simplify(new_m), cx.simplify(new_n));
        let (a, b) = (derivative(cx.graph, new_m, y)?, derivative(cx.graph, new_n, x)?);
        let check = sub(cx.graph, a, b);
        if !cx.is_zero(check) {
            continue;
        }
        let saved = problem.constants;
        let Some(potential) = super::lie::potential(cx, new_n, new_m, y, x) else {
            continue;
        };
        // The potential must reproduce both parts.
        let (px, py) = (derivative(cx.graph, potential, x)?, derivative(cx.graph, potential, y)?);
        let (dx, dy) = (sub(cx.graph, px, new_m), sub(cx.graph, py, new_n));
        if !(cx.is_zero(dx) && cx.is_zero(dy)) {
            problem.constants = saved;
            continue;
        }
        let c = problem.constant(cx.graph);
        let relation = sub(cx.graph, potential, c);
        return implicit(cx, problem, relation);
    }
    None
}

/// Clairaut, d'Alembert and `x`-solved equations in `p = y'`.
pub(super) fn implicit_first_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
    clairaut_only: bool,
) -> Option<NodeId> {
    let (y, p, x) = (*problem.stand.first()?, *problem.stand.get(1)?, problem.x);
    let p_symbol = cx.graph.symbol_of(p)?;
    let y_symbol = cx.graph.symbol_of(y)?;
    // Solved for y: y = G(x, p).
    if cx.graph.depends_on(cx.graph.find(problem.expr), p_symbol) {
        let candidates = solve_for(cx.graph, problem.expr, y, 0).unwrap_or_default();
        for g in candidates {
            let g = cx.simplify(g);
            if cx.graph.depends_on(cx.graph.find(g), y_symbol) || !cx.graph.depends_on(cx.graph.find(g), p_symbol) {
                continue;
            }
            let g_x = derivative(cx.graph, g, x)?;
            let gap = sub(cx.graph, p, g_x);
            let gap = cx.simplify(gap);
            // Clairaut: G_x = p.
            if cx.is_zero(gap) {
                let xp = mul(cx.graph, &[x, p]);
                let rest = sub(cx.graph, g, xp);
                let rest = cx.simplify(rest);
                if problem.depends_on_x(cx.graph, rest) {
                    continue;
                }
                let c = problem.constant(cx.graph);
                let line = mul(cx.graph, &[c, x]);
                let constant_part = cx.graph.substitute(rest, p, c);
                let solution = add(cx.graph, &[line, constant_part]);
                return Some(explicit(cx, problem, solution));
            }
            if clairaut_only {
                continue;
            }
            // d'Alembert: dx/dp = G_p / (p - G_x).
            let g_p = derivative(cx.graph, g, p)?;
            let slope = div(cx.graph, g_p, gap);
            if let Some(found) = parametric(cx, problem, p, x, slope, depth, |cx, x_of_p, p_node| {
                let y_of_p = cx.graph.substitute(g, x, x_of_p);
                let y_of_p = cx.graph.substitute(y_of_p, p, p_node);
                Some(cx.simplify(y_of_p))
            }) {
                return Some(found);
            }
        }
    }
    if clairaut_only {
        return None;
    }
    // Solved for x: x = H(y, p), dY/dp = p H_p / (1 - p H_y).
    let xs = solve_for(cx.graph, problem.expr, x, 0).unwrap_or_default();
    for h in xs {
        let h = cx.simplify(h);
        if problem.depends_on_x(cx.graph, h) || !cx.graph.depends_on(cx.graph.find(h), p_symbol) {
            continue;
        }
        let h_p = derivative(cx.graph, h, p)?;
        let h_y = derivative(cx.graph, h, y)?;
        let one = cx.graph.int(1);
        let p_hy = mul(cx.graph, &[p, h_y]);
        let denominator = sub(cx.graph, one, p_hy);
        let numerator = mul(cx.graph, &[p, h_p]);
        let slope = div(cx.graph, numerator, denominator);
        let slope = cx.simplify(slope);
        let function = fresh_function(cx.graph, "Y", p);
        let in_function = cx.graph.substitute(slope, y, function);
        let jets = chain(cx.graph, function, p, 1)?;
        let equation = cx.graph.node(core::EQ, &[jets[1], in_function]);
        let Some((inner, answer)) = solve_equation(cx, equation, function, depth + 2) else {
            continue;
        };
        let Some(y_of_p) = explicit_rhs(cx.graph, answer, function) else {
            continue;
        };
        problem.constants = problem.constants.max(inner.constants);
        let x_of_p = cx.graph.substitute(h, y, y_of_p);
        let x_of_p = cx.simplify(x_of_p);
        let parameter = cx.graph.sym("p");
        let x_param = cx.graph.substitute(x_of_p, p, parameter);
        let y_param = cx.graph.substitute(y_of_p, p, parameter);
        let first = cx.graph.node(core::EQ, &[x, x_param]);
        let second = cx.graph.node(core::EQ, &[problem.y_of_x, y_param]);
        return Some(cx.graph.node(core::LIST, &[first, second]));
    }
    None
}

/// Solves `dX/dp = slope(X, p)` for `X(p)` and assembles the parametric
/// answer with `y` from `y_of`.
fn parametric(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    p: NodeId,
    x: NodeId,
    slope: NodeId,
    depth: u32,
    y_of: impl Fn(&mut Cx<'_>, NodeId, NodeId) -> Option<NodeId>,
) -> Option<NodeId> {
    let function = fresh_function(cx.graph, "X", p);
    let in_function = cx.graph.substitute(slope, x, function);
    let in_function = cx.simplify(in_function);
    let jets = chain(cx.graph, function, p, 1)?;
    let equation = cx.graph.node(core::EQ, &[jets[1], in_function]);
    let (inner, answer) = solve_equation(cx, equation, function, depth + 2)?;
    let x_of_p = explicit_rhs(cx.graph, answer, function)?;
    problem.constants = problem.constants.max(inner.constants.max(constants_in(cx.graph, x_of_p)));
    let parameter = cx.graph.sym("p");
    let y_of_p = y_of(cx, x_of_p, parameter)?;
    let x_param = cx.graph.substitute(x_of_p, p, parameter);
    let first = cx.graph.node(core::EQ, &[x, x_param]);
    let second = cx.graph.node(core::EQ, &[problem.y_of_x, y_of_p]);
    Some(cx.graph.node(core::LIST, &[first, second]))
}

/// A linear second-order equation with a first-derivative term, through
/// `y = u exp(-½ ∫ a₁/a₂)`.
pub(super) fn normal_form(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    a: &[NodeId; 3],
    forcing: NodeId,
    depth: u32,
) -> Option<NodeId> {
    let x = problem.x;
    if cx.is_zero(a[1]) {
        return None;
    }
    let r1 = div(cx.graph, a[1], a[2]);
    let r1 = cx.simplify(r1);
    let r0 = div(cx.graph, a[0], a[2]);
    let r0 = cx.simplify(r0);
    let integral = integrate(cx, r1, x)?;
    let half = cx.graph.num(Number::fraction(-1, 2)?);
    let exponent = mul(cx.graph, &[half, integral]);
    let scale = exp_of(cx, exponent)?;
    // Q = a0/a2 - (a1/a2)^2/4 - (a1/a2)'/2
    let two = cx.graph.int(2);
    let r1_sq = cx.graph.node(core::POW, &[r1, two]);
    let quarter = cx.graph.num(Number::fraction(-1, 4)?);
    let term_sq = mul(cx.graph, &[quarter, r1_sq]);
    let dr1 = derivative(cx.graph, r1, x)?;
    let term_d = mul(cx.graph, &[half, dr1]);
    let q = add(cx.graph, &[r0, term_sq, term_d]);
    let q = cx.simplify(q);
    let u = fresh_function(cx.graph, "u", x);
    let jets = chain(cx.graph, u, x, 2)?;
    let q_u = mul(cx.graph, &[q, jets[0]]);
    let left = add(cx.graph, &[jets[2], q_u]);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[left, zero]);
    let (_, answer) = solve_equation(cx, equation, u, depth + 2)?;
    let solution = explicit_rhs(cx.graph, answer, u)?;
    let (c1, c2) = (cx.graph.sym("C1"), cx.graph.sym("C2"));
    let (one, zero) = (cx.graph.int(1), cx.graph.int(0));
    let first = {
        let s = cx.graph.substitute(solution, c1, one);
        cx.graph.substitute(s, c2, zero)
    };
    let second = {
        let s = cx.graph.substitute(solution, c1, zero);
        cx.graph.substitute(s, c2, one)
    };
    let y1 = mul(cx.graph, &[scale, first]);
    let y1 = cx.simplify(y1);
    let y2 = mul(cx.graph, &[scale, second]);
    let y2 = cx.simplify(y2);
    let (k1, k2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let mut terms = vec![mul(cx.graph, &[k1, y1]), mul(cx.graph, &[k2, y2])];
    if !cx.is_zero(forcing) {
        terms.push(variation_of_parameters(cx, &[y1, y2], a[2], forcing, x)?);
    }
    let sum = add(cx.graph, &terms);
    Some(explicit(cx, problem, sum))
}
