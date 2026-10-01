//! Methods for equations of order two and higher that reduce them to
//! equations already solvable.
//!
//! * **Missing dependent variable**: `F(x, y^(m), …, y^(n)) = 0` becomes an
//!   equation of order `n - m` for `p = y^(m)`; `y` follows by `m`
//!   integrations.
//! * **Autonomous equations** (`x` absent, order two): `p(y) = y'` turns
//!   `y''` into `p dp/dy`, a first-order equation in `p(y)`; then
//!   `∫ dy / p(y) = x + C`. These are the reductions by the translation
//!   symmetries `∂y` and `∂x`.
//! * **Scale invariance in `y`** (symmetry `y ∂y`): `u = y'/y` gives an
//!   equation of lower order (Riccati for linear second-order equations),
//!   and `y = C exp(∫ u)`.
//! * **Linear equations with variable coefficients** (order two): a first
//!   solution is sought among polynomials, `exp(λ x)`, `x^m` and
//!   `exp(λ x²)`; reduction of order gives the second,
//!   `y₂ = y₁ ∫ exp(-∫ a₁/a₂) / y₁²`, and variation of parameters the
//!   particular solution.
//! * **Special-function equations**: Bessel's equation (and its modified
//!   form) after scaling, Airy's equation when the Airy functions are
//!   registered.
//!
//! Each reduced equation is handed back to the full solver, so a reduced
//! first-order equation can be separable, linear, Bernoulli, Riccati,
//! exact, homogeneous or found by Lie symmetries.

use num_rational::BigRational;
use num_traits::Zero;

use super::Problem;
use super::add;
use super::coefficients_in;
use super::explicit;
use super::implicit;
use super::integrate;
use super::mul;
use super::neg;
use super::occurs;
use super::solve_equation;
use super::sub;
use super::variation_of_parameters;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;

/// `[f, f', f'', …]` up to order `n` as `diff` requests of `f = p(v)`.
fn chain(
    graph: &mut Graph,
    f: NodeId,
    v: NodeId,
    n: usize,
) -> Option<Vec<NodeId>> {
    let diff = graph.ops().lookup("diff")?;
    let mut out = vec![f];
    for _ in 0..n {
        let last = *out.last()?;
        out.push(graph.node(diff, &[last, v]));
    }
    Some(out)
}

/// A fresh undetermined function applied to `v`.
fn fresh_function(
    graph: &mut Graph,
    name: &str,
    v: NodeId,
) -> NodeId {
    let symbol = graph.interner_mut().fresh_symbol(name);
    let f = graph.symbol_node(symbol);
    graph.node(core::APPLY, &[f, v])
}

/// The largest `k` of an integration constant `Ck` in `node`.
fn constants_in(
    graph: &Graph,
    node: NodeId,
) -> usize {
    graph
        .free_symbols(graph.find(node))
        .iter()
        .filter_map(|&s| graph.interner().symbol_name(s).strip_prefix('C')?.parse::<usize>().ok())
        .max()
        .unwrap_or(0)
}

/// The right-hand side of an explicit answer `f(v) = rhs`.
fn explicit_rhs(
    graph: &Graph,
    answer: NodeId,
    unknown: NodeId,
) -> Option<NodeId> {
    match *graph.children(answer) {
        | [lhs, rhs] if graph.op(answer) == core::EQ && lhs == unknown => Some(rhs),
        | _ => None,
    }
}

/// Order reduction when `y` itself does not occur.
pub(super) fn missing_dependent(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    let order = problem.order();
    if order < 2 || occurs(cx.graph, problem.expr, problem.stand[0]) {
        return None;
    }
    let m = (1..order).find(|&k| occurs(cx.graph, problem.expr, problem.stand[k]))?;
    let x = problem.x;
    let p = fresh_function(cx.graph, "p", x);
    let derivatives = chain(cx.graph, p, x, order - m)?;
    let mut reduced = problem.expr;
    for k in m..=order {
        reduced = cx.graph.substitute(reduced, problem.stand[k], derivatives[k - m]);
    }
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[reduced, zero]);
    let (_, answer) = solve_equation(cx, equation, p, depth + 1)?;
    let mut solution = explicit_rhs(cx.graph, answer, p)?;
    problem.constants = problem.constants.max(constants_in(cx.graph, solution));
    for _ in 0..m {
        let integral = integrate(cx, solution, x)?;
        let c = problem.constant(cx.graph);
        solution = add(cx.graph, &[integral, c]);
    }
    Some(explicit(cx, problem, solution))
}

/// Second-order equations without `x`: `y'' = p dp/dy` with `p(y) = y'`.
pub(super) fn autonomous(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    if problem.order() != 2 || problem.depends_on_x(cx.graph, problem.expr) {
        return None;
    }
    let (y, dy, ddy) = (problem.stand[0], problem.stand[1], problem.stand[2]);
    let u_symbol = cx.graph.interner_mut().fresh_symbol("u");
    let u = cx.graph.symbol_node(u_symbol);
    let p = fresh_function(cx.graph, "p", u);
    let derivatives = chain(cx.graph, p, u, 1)?;
    let p_dp = mul(cx.graph, &[p, derivatives[1]]);
    let mut reduced = cx.graph.substitute(problem.expr, ddy, p_dp);
    reduced = cx.graph.substitute(reduced, dy, p);
    reduced = cx.graph.substitute(reduced, y, u);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[reduced, zero]);
    let (_, answer) = solve_equation(cx, equation, p, depth + 1)?;
    let slope = explicit_rhs(cx.graph, answer, p)?;
    problem.constants = problem.constants.max(constants_in(cx.graph, slope));
    // ∫ du / p(u) = x + C, then u = y.
    let minus_one = cx.graph.int(-1);
    let reciprocal = cx.graph.node(core::POW, &[slope, minus_one]);
    let reciprocal = cx.simplify(reciprocal);
    let left = integrate(cx, reciprocal, u)?;
    let c = problem.constant(cx.graph);
    let right = add(cx.graph, &[problem.x, c]);
    let relation = sub(cx.graph, left, right);
    let relation = cx.graph.substitute(relation, u, y);
    implicit(cx, problem, relation)
}

/// Equations invariant under `y → λ y`: `u = y'/y`.
pub(super) fn scale_invariant(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    let order = problem.order();
    if !(2..=3).contains(&order) || !homogeneous_in_y(cx, problem) {
        return None;
    }
    let x = problem.x;
    let u = fresh_function(cx.graph, "u", x);
    let du = chain(cx.graph, u, x, order - 1)?;
    // y'/y = u, y''/y = u' + u², y'''/y = u'' + 3 u u' + u³.
    let one = cx.graph.int(1);
    let two = cx.graph.int(2);
    let three = cx.graph.int(3);
    let u2 = cx.graph.node(core::POW, &[u, two]);
    let mut ratios = vec![one, u];
    ratios.push(add(cx.graph, &[du[1], u2]));
    if order == 3 {
        let u3 = cx.graph.node(core::POW, &[u, three]);
        let mixed = mul(cx.graph, &[three, u, du[1]]);
        ratios.push(add(cx.graph, &[du[2], mixed, u3]));
    }
    // Substitute y^(k) = ratio_k · y with y = 1: the equation is
    // homogeneous, so the common power of y drops out.
    let mut reduced = problem.expr;
    for (k, &ratio) in ratios.iter().enumerate() {
        reduced = cx.graph.substitute(reduced, problem.stand[k], ratio);
    }
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[reduced, zero]);
    let (_, answer) = solve_equation(cx, equation, u, depth + 1)?;
    let rate = explicit_rhs(cx.graph, answer, u)?;
    problem.constants = problem.constants.max(constants_in(cx.graph, rate));
    let integral = integrate(cx, rate, x)?;
    let growth = super::exp_of(cx, integral)?;
    let c = problem.constant(cx.graph);
    let solution = mul(cx.graph, &[c, growth]);
    Some(explicit(cx, problem, solution))
}

/// Whether the equation is homogeneous in `y` and its derivatives:
/// scaling them all by λ multiplies it by a power of λ (checked
/// numerically).
fn homogeneous_in_y(
    cx: &Cx<'_>,
    problem: &Problem,
) -> bool {
    if !occurs(cx.graph, problem.expr, problem.stand[0]) {
        return false;
    }
    let symbols: Vec<_> = problem.stand.iter().filter_map(|&s| cx.graph.symbol_of(s)).collect();
    let others: Vec<_> = cx
        .graph
        .free_symbols(cx.graph.find(problem.expr))
        .iter()
        .copied()
        .filter(|s| !symbols.contains(s))
        .collect();
    let value = |cx: &Cx<'_>, scale: f64, point: f64| -> Option<f64> {
        let mut env = Env::numeric(0.0);
        for (i, &s) in symbols.iter().enumerate() {
            env.bind(s, scale * (0.7 + 0.37 * f64::from(u32::try_from(i).unwrap_or(0)) + point));
        }
        for (i, &s) in others.iter().enumerate() {
            env.bind(s, 0.9 + 0.23 * f64::from(u32::try_from(i).unwrap_or(0)) + point);
        }
        cx.graph.eval(problem.expr, &env).filter(|v| v.is_finite() && v.abs() > 1e-12)
    };
    let mut degree = None;
    for point in [0.0, 0.41, 1.3] {
        let (Some(base), Some(scaled)) = (value(cx, 1.0, point), value(cx, 2.0, point)) else {
            return false;
        };
        let ratio = scaled / base;
        if ratio <= 0.0 {
            return false;
        }
        let k = ratio.log2();
        if (k - k.round()).abs() > 1e-9 {
            return false;
        }
        if degree.is_some_and(|d: f64| (d - k).abs() > 1e-9) {
            return false;
        }
        degree = Some(k);
    }
    degree.is_some()
}

/// The coefficients `a₀, a₁, a₂` and forcing `r` of a linear second-order
/// equation `a₂ y'' + a₁ y' + a₀ y = r`.
fn linear_second_order(
    cx: &mut Cx<'_>,
    problem: &Problem,
) -> Option<([NodeId; 3], NodeId)> {
    if problem.order() != 2 {
        return None;
    }
    let mut rest = problem.expr;
    let mut a = [NodeId::NONE; 3];
    for k in (0..=2).rev() {
        let parts = coefficients_in(cx.graph, rest, problem.stand[k])?;
        match *parts.as_slice() {
            | [constant, linear] => {
                a[k] = linear;
                rest = constant;
            },
            | [constant] => {
                a[k] = cx.graph.int(0);
                rest = constant;
            },
            | _ => return None,
        }
    }
    let stands = &problem.stand;
    if stands.iter().any(|&s| occurs(cx.graph, rest, s) || a.iter().any(|&c| occurs(cx.graph, c, s))) || cx.is_zero(a[2]) {
        return None;
    }
    let forcing = neg(cx.graph, rest);
    Some((a, cx.simplify(forcing)))
}

/// `a₂ y'' + a₁ y' + a₀ y` for a candidate `y`, simplified.
fn residual(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    y: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let d1 = derivative(cx.graph, y, x)?;
    let d2 = derivative(cx.graph, d1, x)?;
    let terms = [mul(cx.graph, &[a[0], y]), mul(cx.graph, &[a[1], d1]), mul(cx.graph, &[a[2], d2])];
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

/// Solutions of the homogeneous equation among simple families:
/// polynomials, `exp(λ x)`, `x^m` and `exp(λ x²)`.
fn family_solutions(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
) -> Vec<NodeId> {
    let mut found = Vec::new();
    collect_solutions(cx, a, x, &mut found);
    found
}

fn collect_solutions(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
    found: &mut Vec<NodeId>,
) -> Option<()> {
    // Polynomials of degree ≤ 4 with unknown coefficients.
    for degree in 1..=4_u32 {
        let mut unknowns = Vec::new();
        let mut terms = Vec::new();
        for i in 0..=degree {
            let symbol = cx.graph.interner_mut().fresh_symbol("c");
            let c = cx.graph.symbol_node(symbol);
            unknowns.push(c);
            let e = cx.graph.int(i64::from(i));
            let power = cx.graph.node(core::POW, &[x, e]);
            terms.push(mul(cx.graph, &[c, power]));
        }
        let candidate = add(cx.graph, &terms);
        let d1 = derivative(cx.graph, candidate, x)?;
        let d2 = derivative(cx.graph, d1, x)?;
        let lhs = [mul(cx.graph, &[a[0], candidate]), mul(cx.graph, &[a[1], d1]), mul(cx.graph, &[a[2], d2])];
        let lhs = add(cx.graph, &lhs);
        if let Some(basis) = super::lie::solve_identity(cx.graph, lhs, &unknowns) {
            for v in &basis {
                let mut y = candidate;
                for (u, value) in unknowns.iter().zip(v) {
                    let value = cx.graph.num(Number::rat(value.clone()));
                    y = cx.graph.substitute(y, *u, value);
                }
                let y = cx.simplify(y);
                if !cx.is_zero(y) {
                    found.push(y);
                }
            }
            if !basis.is_empty() {
                break;
            }
        }
    }
    // exp(λ x), x^m and exp(λ x²): λ from the first coefficient of the
    // residual, checked on the whole residual.
    let exp = cx.graph.ops().lookup("exp")?;
    let symbol = cx.graph.interner_mut().fresh_symbol("lambda");
    let lambda = cx.graph.symbol_node(symbol);
    let two = cx.graph.int(2);
    let one = cx.graph.int(1);
    let minus_one = cx.graph.int(-1);
    // (y, y'/y, y''/y) for each family.
    let shapes = {
        let lx = mul(cx.graph, &[lambda, x]);
        let x2 = cx.graph.node(core::POW, &[x, two]);
        let lx2 = mul(cx.graph, &[lambda, x2]);
        let l2 = cx.graph.node(core::POW, &[lambda, two]);
        let inv_x = cx.graph.node(core::POW, &[x, minus_one]);
        let minus_two = cx.graph.int(-2);
        let inv_x2 = cx.graph.node(core::POW, &[x, minus_two]);
        let l_minus_1 = add(cx.graph, &[lambda, minus_one]);
        let four = cx.graph.int(4);
        [
            (cx.graph.node(exp, &[lx]), lambda, l2),
            (
                cx.graph.node(core::POW, &[x, lambda]),
                mul(cx.graph, &[lambda, inv_x]),
                mul(cx.graph, &[lambda, l_minus_1, inv_x2]),
            ),
            (cx.graph.node(exp, &[lx2]), mul(cx.graph, &[two, lambda, x]), {
                let a = mul(cx.graph, &[four, l2, x2]);
                let b = mul(cx.graph, &[two, lambda]);
                add(cx.graph, &[a, b])
            }),
        ]
    };
    let _ = one;
    for (shape, r1, r2) in shapes {
        let terms = [a[0], mul(cx.graph, &[a[1], r1]), mul(cx.graph, &[a[2], r2])];
        let reduced = add(cx.graph, &terms);
        let reduced = cx.simplify(reduced);
        let Some(numerator) = super::numerator_of(cx.graph, reduced) else {
            continue;
        };
        let Some(parts) = coefficients_in(cx.graph, numerator, x) else {
            continue;
        };
        let Some(&first) = parts.iter().find(|&&p| !cx.is_zero(p)) else {
            continue;
        };
        let Some(roots) = super::solve_for(cx.graph, first, lambda, 0) else {
            continue;
        };
        for root in roots {
            if cx.graph.number_of(root).is_none_or(|n| n.to_rational().is_none_or(|r| r.is_zero())) {
                continue;
            }
            let y = cx.graph.substitute(shape, lambda, root);
            let y = cx.simplify(y);
            if let Some(check) = residual(cx, a, y, x) {
                if cx.is_zero(check) {
                    found.push(y);
                }
            }
        }
    }
    Some(())
}

/// Whether `y1 / y2` varies with `x`.
fn independent(
    cx: &mut Cx<'_>,
    y1: NodeId,
    y2: NodeId,
    x: NodeId,
) -> bool {
    let Some(symbol) = cx.graph.symbol_of(x) else {
        return false;
    };
    let ratio = super::div(cx.graph, y1, y2);
    let ratio = cx.simplify(ratio);
    cx.graph.depends_on(cx.graph.find(ratio), symbol)
}

/// Linear second-order equations with variable coefficients, by a found
/// solution and reduction of order.
pub(super) fn variable_coefficients(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    _depth: u32,
) -> Option<NodeId> {
    let (a, forcing) = linear_second_order(cx, problem)?;
    let x = problem.x;
    let solutions = family_solutions(cx, &a, x);
    let y1 = *solutions.first()?;
    let second = solutions.iter().skip(1).copied().find(|&y| independent(cx, y1, y, x));
    let y2 = match second {
        | Some(y2) => y2,
        | None => reduction_of_order(cx, &a, y1, x)?,
    };
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let mut terms = vec![mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
    if !cx.is_zero(forcing) {
        terms.push(variation_of_parameters(cx, &[y1, y2], a[2], forcing, x)?);
    }
    let solution = add(cx.graph, &terms);
    Some(explicit(cx, problem, solution))
}

/// `y₂ = y₁ ∫ exp(-∫ a₁/a₂) / y₁²`.
fn reduction_of_order(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    y1: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    // y₂ = y₁ ∫ exp(-∫ a₁/a₂) / y₁²
    let ratio = super::div(cx.graph, a[1], a[2]);
    let ratio = cx.simplify(ratio);
    let integral = integrate(cx, ratio, x)?;
    let negated = neg(cx.graph, integral);
    let abel = super::exp_of(cx, negated)?;
    let two = cx.graph.int(2);
    let square = cx.graph.node(core::POW, &[y1, two]);
    let integrand = super::div(cx.graph, abel, square);
    let integrand = cx.simplify(integrand);
    let inner = integrate(cx, integrand, x)?;
    let y2 = mul(cx.graph, &[y1, inner]);
    Some(cx.simplify(y2))
}

/// Bessel's equation `x² y'' + x y' + (b² x² - ν²) y = 0` (and the
/// modified equation, `-b² x²`), and Airy's equation `y'' = b x y`.
pub(super) fn special_function_equation(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    _depth: u32,
) -> Option<NodeId> {
    let (a, forcing) = linear_second_order(cx, problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    // Normalise to y'' + P y' + Q y = 0.
    let p = super::div(cx.graph, a[1], a[2]);
    let p = cx.simplify(p);
    let q = super::div(cx.graph, a[0], a[2]);
    let q = cx.simplify(q);
    // Bessel: P = 1/x and x² Q = b² x² - ν² (b², ν² constants).
    let xp = mul(cx.graph, &[x, p]);
    let xp = cx.simplify(xp);
    if cx.graph.number_of(xp).is_some_and(Number::is_one) {
        let two = cx.graph.int(2);
        let x2 = cx.graph.node(core::POW, &[x, two]);
        let scaled = mul(cx.graph, &[x2, q]);
        let scaled = cx.simplify(scaled);
        let parts = coefficients_in(cx.graph, scaled, x)?;
        let (&nu_sq_neg, rest) = parts.split_first()?;
        let b_sq = match rest {
            | [zero, b_sq] if cx.is_zero(*zero) => *b_sq,
            | _ => return None,
        };
        let b_value = cx.graph.eval(b_sq, &Env::numeric(0.0))?;
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let nu_sq = neg(cx.graph, nu_sq_neg);
        let nu = cx.graph.node(core::POW, &[nu_sq, half]);
        let nu = cx.simplify(nu);
        let (first, second, b) = if b_value > 0.0 {
            ("besselj", "bessely", cx.graph.node(core::POW, &[b_sq, half]))
        } else {
            let minus = neg(cx.graph, b_sq);
            ("besseli", "besselk", cx.graph.node(core::POW, &[minus, half]))
        };
        let b = cx.simplify(b);
        let (f1, f2) = (cx.graph.ops().lookup(first)?, cx.graph.ops().lookup(second)?);
        let bx = mul(cx.graph, &[b, x]);
        let (y1, y2) = (cx.graph.node(f1, &[nu, bx]), cx.graph.node(f2, &[nu, bx]));
        let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
        let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
        let solution = add(cx.graph, &terms);
        return Some(explicit(cx, problem, solution));
    }
    // Airy: P = 0 and Q = -b x with b > 0: y = C1 Ai(b^(1/3) x) + C2 Bi(b^(1/3) x).
    if cx.is_zero(p) {
        let parts = coefficients_in(cx.graph, q, x)?;
        let [zero, slope] = parts.as_slice() else {
            return None;
        };
        if !cx.is_zero(*zero) {
            return None;
        }
        let (ai, bi) = (cx.graph.ops().lookup("airyai")?, cx.graph.ops().lookup("airybi")?);
        let b = neg(cx.graph, *slope);
        let third = cx.graph.num(Number::fraction(1, 3)?);
        let scale = cx.graph.node(core::POW, &[b, third]);
        let argument = mul(cx.graph, &[scale, x]);
        let argument = cx.simplify(argument);
        let (y1, y2) = (cx.graph.node(ai, &[argument]), cx.graph.node(bi, &[argument]));
        let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
        let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
        let solution = add(cx.graph, &terms);
        return Some(explicit(cx, problem, solution));
    }
    let _ = BigRational::zero();
    None
}
