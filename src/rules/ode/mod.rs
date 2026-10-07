//! Ordinary differential equations.
//!
//! `dsolve(equation, y(x))` asks for the general solution of an ODE for
//! the undetermined function `y(x)`; derivatives are written
//! `diff(y(x), x)`. The answer is an equation `y(x) = ...` with integration
//! constants `C1`, `C2`, ... (or an implicit relation when the solution
//! cannot be solved for `y`). A third argument, a list of conditions such
//! as `list(y(0) = 1)`, determines the constants.
//!
//! The symbolic kernel classifies the equation and applies the matching
//! method: separable, linear, Bernoulli, exact and homogeneous equations
//! of first order; Riccati equations with a simple particular solution;
//! linear equations with constant coefficients of any order (with
//! variation of parameters for a right-hand side); Cauchy–Euler
//! equations. Whatever it produces is substituted back into the equation
//! and checked numerically before it is accepted.
//!
//! `odeint(f, y, x, x0, y0, x1)` is the numeric counterpart: the value at
//! `x1` of the solution of `y' = f(x, y)`, `y(x0) = y0`, by an adaptive
//! Runge–Kutta method; with lists for `f`, `y` and `y0` it integrates a
//! system and returns the list of final values.
//!
//! Beyond the classical classes, `dsolve` tries in turn: Bessel and Airy
//! equations; variable-coefficient linear equations (a family solution,
//! then reduction of order); missing dependent variable; autonomous,
//! scale-invariant and equidimensional equations; and the Lie
//! point-symmetry method (prolonged determining equations, canonical
//! coordinates, reduction of order). With a list of unknowns it solves
//! constant-coefficient linear systems; with initial conditions and a
//! discontinuous forcing it falls back to the Laplace transform.
//!
//! Further classes: Riccati equations by a particular solution among
//! simple shapes (`k x^n`, `k e^(±x)`, `k sin x`, …) or, failing that,
//! through the linear second-order equation for `u` (`y = -u'/(q2 u)`,
//! written with one constant); Bernoulli with rational exponents; exact
//! equations with integrating factors `μ(x)`, `μ(y)`, `μ(x + y)`,
//! `μ(x y)`; equations that are nonlinear in `y'` (branches of `y'`;
//! Clairaut `y = C x + g(C)`; d'Alembert/Lagrange and `x`-solved equations
//! as parametric answers `list(x = X(p, C1), y(x) = Y(p, C1))` with
//! `p = y'`); second-order linear equations by a polynomial or
//! exponential solution and reduction of order (the second solution is
//! an inert `integral(…)` when it has no closed form) or, without one, by
//! the normal form `y = u exp(-½ ∫ a₁/a₂)`.
//!
//! Systems `dsolve(list(eqs), list(funcs))`: constant linear systems
//! (with forcing) by elimination; higher-order systems by reduction to
//! first order; triangular and decoupled systems one equation at a time;
//! `X' = a(t) M X` through `s = ∫ a`; two autonomous equations through
//! `dy/dx = g/f` (the answer is the first integral `Φ(x(t), y(t)) = 0`
//! with a constant: the energy of a Hamiltonian system, the implicit
//! integral of a Lotka–Volterra or SIR system); polynomial first integrals
//! up to degree three for autonomous polynomial systems in three or four
//! unknowns (`list(x(t)^2 + y(t)^2 - C1 = 0, …)`).
//!
//! | operator | value |
//! |---|---|
//! | `ode_series(eq, y(x), x0, n)` | the power-series solution about `x0` to order `n` |
//! | `rsolve_z(eq, y(n), conditions?)` | a linear recurrence solved by the z transform |
//! | `ode_classify(eq, y(x))` | the list of classes the equation belongs to |
//! | `ode_symmetries(eq, y(x))` | Lie point symmetries `list(ξ, η)` (polynomial ansatz) |
//! | `ode_fourier(eq, y(x))` | the solution on the whole line decaying at `±∞` of a linear constant-coefficient equation, by the Fourier transform (inverted by residues) |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::backend::Backend;
use crate::backend::Interpreter;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::kernels::ode::solve_adaptive;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::univariate;

use super::calculus::antiderivative;
use super::calculus::calculus;
use super::calculus::derivative;
use super::solve::as_expression;
use super::solve::solve;
use super::solve::solve_for;
use super::solve::solve_linear;

mod classes;
mod lie;
mod nonlinear;
#[cfg(test)]
mod battery;
mod reduce;
mod series_method;
mod systems;

/// The differential-equation rule set.
#[must_use]
pub fn ode() -> RuleSet {
    RuleSet::new("ode", install).needs(calculus()).needs(solve())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let dsolve = i.op(OpDescriptor::new("dsolve", Arity::Variadic).flags(OpFlags::HEAVY).cost(100))?;
    let odeint = i.op(OpDescriptor::new("odeint", Arity::Fixed(6)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("ode/dsolve", Tier::Reduce, Dsolve { dsolve });
    let series = i.op(OpDescriptor::new("ode_series", Arity::Fixed(4)).flags(OpFlags::HEAVY).cost(100))?;
    let rsolve = i.op(OpDescriptor::new("rsolve_z", Arity::Variadic).flags(OpFlags::HEAVY).cost(100))?;
    let classify = i.op(OpDescriptor::new("ode_classify", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    let symmetries = i.op(OpDescriptor::new("ode_symmetries", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    let fourier = i.op(OpDescriptor::new("ode_fourier", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("ode/extras", Tier::Reduce, Extras { series, rsolve, classify, symmetries, fourier });
    i.kernel("ode/odeint", Tier::Reduce, Odeint { odeint });
    Ok(())
}

/// Small term-building helpers.
fn mul(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn add(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn neg(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    mul(graph, &[minus_one, node])
}

fn sub(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let negated = neg(graph, b);
    add(graph, &[a, negated])
}

fn inv(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    graph.node(core::POW, &[node, minus_one])
}

fn div(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let inverse = inv(graph, b);
    mul(graph, &[a, inverse])
}

fn call(
    graph: &mut Graph,
    name: &str,
    arg: NodeId,
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, &[arg])
}

/// `exp(e)` with logarithms taken out: `exp(Σ cᵢ ln uᵢ + r) = ∏ uᵢ^cᵢ exp(r)`
/// for rational `cᵢ`, which keeps integrating factors and Abel's formula
/// rational where they are.
fn exp_of(
    cx: &mut Cx<'_>,
    e: NodeId,
) -> Option<NodeId> {
    let (exp, ln) = (cx.graph.ops().lookup("exp")?, cx.graph.ops().lookup("ln")?);
    let e = cx.simplify(e);
    // c·(a + b + …) → c·a + c·b + … for a numeric c.
    let pair = match cx.graph.children(e).to_vec().as_slice() {
        | &[a, b] if cx.graph.number_of(a).is_some() => Some((a, b)),
        | &[a, b] if cx.graph.number_of(b).is_some() => Some((b, a)),
        | _ => None,
    };
    let terms: Vec<NodeId> = match (cx.graph.op(e), pair) {
        | (op, Some((c, sum))) if op == core::MUL && cx.graph.op(sum) == core::ADD => {
            cx.graph.children(sum).to_vec().iter().map(|&t| mul(cx.graph, &[c, t])).collect()
        },
        | (op, _) if op == core::ADD => cx.graph.children(e).to_vec(),
        | _ => vec![e],
    };
    let mut factors = Vec::new();
    let mut rest = Vec::new();
    for t in terms {
        let (coefficient, inner) = match (cx.graph.op(t), cx.graph.children(t).to_vec().as_slice()) {
            | (op, &[u]) if op == ln => (cx.graph.int(1), Some(u)),
            | (op, children) if op == core::MUL => {
                let numbers: Vec<NodeId> = children.iter().copied().filter(|&c| cx.graph.number_of(c).is_some()).collect();
                let logs: Vec<NodeId> = children.iter().copied().filter(|&c| cx.graph.op(c) == ln).collect();
                if numbers.len() == 1 && logs.len() == 1 && children.len() == 2 {
                    (numbers[0], cx.graph.children(logs[0]).first().copied())
                } else {
                    (cx.graph.int(1), None)
                }
            },
            | _ => (cx.graph.int(1), None),
        };
        match inner {
            | Some(u) => factors.push(cx.graph.node(core::POW, &[u, coefficient])),
            | None => rest.push(t),
        }
    }
    if !rest.is_empty() {
        let sum = add(cx.graph, &rest);
        factors.push(cx.graph.node(exp, &[sum]));
    }
    let product = mul(cx.graph, &factors);
    Some(cx.simplify(product))
}

/// One ODE problem with the unknown function replaced by symbols.
pub(super) struct Problem {
    /// The independent variable.
    x: NodeId,
    x_symbol: SymbolId,
    /// The term `y(x)`.
    y_of_x: NodeId,
    /// Stand-ins for `y`, `y'`, `y''`, ...; `stand[k]` replaces the k-th
    /// derivative.
    stand: Vec<NodeId>,
    /// The equation as an expression in `x` and the stand-ins, equal to
    /// zero.
    expr: NodeId,
    /// Next integration constant to hand out.
    constants: usize,
    /// The answer is correct only on part of its domain, so the spot check
    /// of [`verified`] is skipped.
    trusted: bool,
}

impl Problem {
    const fn order(&self) -> usize {
        self.stand.len() - 1
    }

    fn constant(
        &mut self,
        graph: &mut Graph,
    ) -> NodeId {
        self.constants += 1;
        graph.sym(&format!("C{}", self.constants))
    }

    fn depends_on_x(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> bool {
        graph.depends_on(graph.find(node), self.x_symbol)
    }

    fn depends_on_y(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> bool {
        self.stand.first().and_then(|&y| graph.as_symbol(y)).is_some_and(|s| graph.depends_on(graph.find(node), s))
    }
}

/// Finds the derivatives of `y(x)` occurring in `term` and replaces them,
/// highest order first, by fresh symbols.
fn parse(
    graph: &mut Graph,
    equation: NodeId,
    unknown: NodeId,
) -> Option<Problem> {
    let diff = graph.ops().lookup("diff")?;
    let &[_, x] = graph.children(unknown) else {
        return None;
    };
    if graph.op(unknown) != core::APPLY {
        return None;
    }
    let x_symbol = graph.symbol_of(x)?;
    let expr = as_expression(graph, equation);
    let mut expr = best(graph, expr)?;
    // The chain y, diff(y, x), diff(diff(y, x), x), ... as long as the
    // nodes occur in the equation.
    let mut chain = vec![unknown];
    loop {
        let last = *chain.last()?;
        let next = graph.node(diff, &[last, x]);
        if !occurs(graph, expr, next) || chain.len() > 8 {
            break;
        }
        chain.push(next);
    }
    if chain.len() < 2 {
        return None;
    }
    let mut stand = vec![NodeId::NONE; chain.len()];
    for (k, &node) in chain.iter().enumerate().rev() {
        let symbol = graph.interner_mut().fresh_symbol(&format!("y{k}"));
        let replacement = graph.symbol_node(symbol);
        expr = graph.replace_subterm(expr, node, replacement);
        *stand.get_mut(k)? = replacement;
    }
    // Any remaining occurrence of the function (at another argument, say)
    // is beyond this solver.
    if occurs_function(graph, expr, unknown) {
        return None;
    }
    Some(Problem { x, x_symbol, y_of_x: unknown, stand, expr, constants: 0, trusted: false })
}

fn occurs(
    graph: &Graph,
    term: NodeId,
    needle: NodeId,
) -> bool {
    let mut stack = vec![term];
    let mut seen = Vec::new();
    while let Some(node) = stack.pop() {
        if node == needle {
            return true;
        }
        if seen.contains(&node) {
            continue;
        }
        seen.push(node);
        stack.extend_from_slice(graph.children(node));
    }
    false
}

/// Whether the function symbol of `unknown` still occurs in `term`.
fn occurs_function(
    graph: &Graph,
    term: NodeId,
    unknown: NodeId,
) -> bool {
    graph.children(unknown).first().is_some_and(|&f| occurs(graph, term, f))
}

/// `∫ f dx`, simplified.
fn integrate(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let primitive = antiderivative(cx, f, x)?;
    Some(cx.simplify(primitive))
}

/// The coefficients of `expr` as a polynomial in the symbol `v`, provided
/// they are free of `v`: entry `k` multiplies `v^k`.
fn coefficients_in(
    graph: &mut Graph,
    expr: NodeId,
    v: NodeId,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(v)?;
    let mut gens = Gens::default();
    let gv = gens.index(graph, v);
    let poly = from_term(graph, &mut gens, expr, Limits::default())?;
    for g in poly.support() {
        if g != gv && gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)) {
            return None;
        }
    }
    Some(poly.coefficients_in(gv).iter().map(|c| to_term(graph, &gens, c)).collect())
}

/// First-order equations `F(x, y, y') = 0`.
fn first_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    let dy = *problem.stand.get(1)?;
    // Exact equations are recognised on the original form M + N y' = 0.
    if let Some(found) = exact(cx, problem) {
        return Some(found);
    }
    let roots = solve_for(cx.graph, problem.expr, dy, 0);
    let several = roots.as_ref().is_none_or(|r| r.len() > 1);
    if let Some(found) = classes::implicit_first_order(cx, problem, depth, true) {
        return Some(found);
    }
    // Each branch of `y'` in turn.
    for rhs in roots.unwrap_or_default() {
        let rhs = cx.simplify(rhs);
        let saved = problem.constants;
        if let Some(found) = first_order_explicit(cx, problem, rhs, depth) {
            return Some(found);
        }
        problem.constants = saved;
    }
    if several {
        return classes::implicit_first_order(cx, problem, depth, false);
    }
    None
}

/// First-order equations `y' = f(x, y)`.
fn first_order_explicit(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    rhs: NodeId,
    depth: u32,
) -> Option<NodeId> {
    let (y, _) = (*problem.stand.first()?, *problem.stand.get(1)?);
    let x = problem.x;
    // Linear: y' = q(x) - p(x) y.
    if let Some(coefficients) = coefficients_in(cx.graph, rhs, y) {
        if coefficients.len() <= 2 {
            let q = coefficients.first().copied().unwrap_or_else(|| cx.graph.int(0));
            let minus_p = coefficients.get(1).copied().unwrap_or_else(|| cx.graph.int(0));
            // mu = exp(∫p) = exp(-∫(-p))
            let integral = integrate(cx, minus_p, x)?;
            let growth = call(cx.graph, "exp", integral)?;
            let decay_arg = neg(cx.graph, integral);
            let decay = call(cx.graph, "exp", decay_arg)?;
            let weighted = mul(cx.graph, &[decay, q]);
            let weighted = cx.simplify(weighted);
            let particular = integrate(cx, weighted, x)?;
            let c = problem.constant(cx.graph);
            let sum = add(cx.graph, &[particular, c]);
            let solution = mul(cx.graph, &[growth, sum]);
            return Some(explicit(cx, problem, solution));
        }
        // Bernoulli: y' = a(x) y + b(x) y^n, n != 0, 1.
        let nonzero: Vec<usize> = coefficients
            .iter()
            .enumerate()
            .filter(|&(_, &c)| !cx.graph.number_of(c).is_some_and(Number::is_zero))
            .map(|(k, _)| k)
            .collect();
        if let [1, n] | [n] = *nonzero.as_slice() {
            if n >= 2 {
                return bernoulli(cx, problem, &coefficients, n);
            }
        }
        // Riccati: y' = q0 + q1 y + q2 y^2.
        if coefficients.len() == 3 {
            if let Some(found) = riccati(cx, problem, rhs, &coefficients) {
                return Some(found);
            }
            if let Some(found) = classes::riccati_linearised(cx, problem, &coefficients, depth) {
                return Some(found);
            }
        }
    }

    if let Some(found) = classes::bernoulli_rational(cx, problem, rhs) {
        return Some(found);
    }

    // Separable: f = g(x) h(y).
    if let Some(found) = separable(cx, problem, rhs) {
        return Some(found);
    }
    if let Some(found) = homogeneous(cx, problem, rhs) {
        return Some(found);
    }
    // Integrating factors of M + N y' = 0.
    let dy = *problem.stand.get(1)?;
    if let Some(&[m, n]) = coefficients_in(cx.graph, problem.expr, dy).as_deref() {
        if let Some(found) = classes::integrating_factor(cx, problem, m, n) {
            return Some(found);
        }
    }
    // Lie point symmetries cover what the classical recipes miss.
    lie::first_order(cx, problem, rhs)
}

/// Wraps an explicit solution as the equation `y(x) = solution`.
fn explicit(
    cx: &mut Cx<'_>,
    problem: &Problem,
    solution: NodeId,
) -> NodeId {
    let solution = cx.simplify(solution);
    cx.graph.node(core::EQ, &[problem.y_of_x, solution])
}

/// Turns an implicit relation `relation(x, y) = 0` into an answer: solved
/// for `y` when possible, otherwise left implicit in `y(x)`.
fn implicit(
    cx: &mut Cx<'_>,
    problem: &Problem,
    relation: NodeId,
) -> Option<NodeId> {
    let y = *problem.stand.first()?;
    let relation = cx.simplify(relation);
    if let Some(solutions) = solve_for(cx.graph, relation, y, 0) {
        // The first branch that really solves the equation (a branch of a
        // multivalued inverse may only hold on part of the domain).
        let lambert: Vec<OpId> = ["lambertw", "lambertw_m1"].iter().filter_map(|n| cx.graph.ops().lookup(n)).collect();
        for solution in solutions {
            // A Lambert-W branch says less than the relation itself.
            if lambert.iter().any(|&op| occurs_op(cx.graph, solution, op)) {
                continue;
            }
            let answer = explicit(cx, problem, solution);
            if verified(cx, problem, answer) {
                return Some(answer);
            }
        }
    }
    let in_function = cx.graph.substitute(relation, y, problem.y_of_x);
    let zero = cx.graph.int(0);
    Some(cx.graph.node(core::EQ, &[in_function, zero]))
}

fn separable(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    rhs: NodeId,
) -> Option<NodeId> {
    let y = *problem.stand.first()?;
    let raw: Vec<NodeId> =
        if cx.graph.op(rhs) == core::MUL { cx.graph.children(rhs).to_vec() } else { vec![rhs] };
    // exp(a + b) = exp(a) exp(b): lets exp(x - y) separate.
    let exp = cx.graph.ops().lookup("exp");
    let mut factors = Vec::with_capacity(raw.len());
    for factor in raw {
        match *cx.graph.children(factor) {
            | [argument] if Some(cx.graph.op(factor)) == exp && cx.graph.op(argument) == core::ADD => {
                for term in cx.graph.children(argument).to_vec() {
                    factors.push(cx.graph.try_node(cx.graph.op(factor), &[term])?);
                }
            },
            | _ => factors.push(factor),
        }
    }
    // (a b)^k = a^k b^k for integer k: the factors of a quotient.
    let mut flat = Vec::with_capacity(factors.len());
    let mut stack = factors;
    stack.reverse();
    while let Some(factor) = stack.pop() {
        let children = cx.graph.children(factor).to_vec();
        match (cx.graph.op(factor), children.as_slice()) {
            | (op, &[base, e]) if op == core::POW && cx.graph.op(base) == core::MUL && cx.graph.number_of(e).is_some_and(Number::is_integer) => {
                for part in cx.graph.children(base).to_vec().into_iter().rev() {
                    stack.push(cx.graph.node(core::POW, &[part, e]));
                }
            },
            | (op, _) if op == core::MUL => {
                for part in children.into_iter().rev() {
                    stack.push(part);
                }
            },
            | _ => flat.push(factor),
        }
    }
    let classify = |cx: &mut Cx<'_>, problem: &Problem, factors: &[NodeId]| -> Option<(Vec<NodeId>, Vec<NodeId>)> {
        let (mut in_x, mut in_y) = (Vec::new(), Vec::new());
        for &factor in factors {
            match (problem.depends_on_x(cx.graph, factor), problem.depends_on_y(cx.graph, factor)) {
                | (true, true) => return None,
                | (_, true) => in_y.push(factor),
                | _ => in_x.push(factor),
            }
        }
        Some((in_x, in_y))
    };
    // A quotient that was expanded may still separate after factoring.
    let (in_x, in_y) = match classify(cx, problem, &flat) {
        | Some(split) => split,
        | None => {
            let factors = factor_quotient(cx, rhs)?;
            classify(cx, problem, &factors)?
        },
    };
    if in_y.is_empty() {
        return None;
    }
    let g = mul(cx.graph, &in_x);
    let h = mul(cx.graph, &in_y);
    // ∫ dy / h(y) = ∫ g(x) dx + C
    let reciprocal = inv(cx.graph, h);
    let reciprocal = cx.simplify(reciprocal);
    let left = integrate(cx, reciprocal, y)?;
    let right = integrate(cx, g, problem.x)?;
    let c = problem.constant(cx.graph);
    let right = add(cx.graph, &[right, c]);
    let relation = sub(cx.graph, left, right);
    implicit(cx, problem, relation)
}

/// The irreducible factors of a rational function over `Q`, as powers
/// (numerator factors with positive and denominator factors with negative
/// exponents).
fn factor_quotient(
    cx: &mut Cx<'_>,
    rhs: NodeId,
) -> Option<Vec<NodeId>> {
    let mut gens = Gens::default();
    let fraction = crate::rules::poly::ratio(cx.graph, &mut gens, rhs, Limits::default())?;
    let mut out = Vec::new();
    for (poly, sign) in [(fraction.numer, 1_i64), (fraction.denom, -1)] {
        let support = poly.support();
        let (unit, pieces) = if support.is_empty() {
            (poly.as_constant()?, Vec::new())
        } else {
            crate::rules::poly::multifactor::factor(&poly, &support)?
        };
        let unit = cx.graph.num(unit);
        out.push(if sign < 0 { inv(cx.graph, unit) } else { unit });
        for (piece, multiplicity) in pieces {
            let base = to_term(cx.graph, &gens, &piece);
            let e = cx.graph.int(sign * i64::from(multiplicity));
            out.push(cx.graph.node(core::POW, &[base, e]));
        }
    }
    Some(out)
}

fn bernoulli(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    coefficients: &[NodeId],
    n: usize,
) -> Option<NodeId> {
    // y' = a y + b y^n; with v = y^(1-n): v' = (1-n) a v + (1-n) b.
    let x = problem.x;
    let a = coefficients.get(1).copied()?;
    let b = coefficients.get(n).copied()?;
    let k = cx.graph.int(1 - i64::try_from(n).ok()?);
    let rate = mul(cx.graph, &[k, a]);
    let forcing = mul(cx.graph, &[k, b]);
    let integral = integrate(cx, rate, x)?;
    let growth = call(cx.graph, "exp", integral)?;
    let decay_arg = neg(cx.graph, integral);
    let decay = call(cx.graph, "exp", decay_arg)?;
    let weighted = mul(cx.graph, &[decay, forcing]);
    let weighted = cx.simplify(weighted);
    let particular = integrate(cx, weighted, x)?;
    let c = problem.constant(cx.graph);
    let sum = add(cx.graph, &[particular, c]);
    let v = mul(cx.graph, &[growth, sum]);
    // y = v^(1/(1-n))
    let exponent = inv(cx.graph, k);
    let solution = cx.graph.node(core::POW, &[v, exponent]);
    Some(explicit(cx, problem, solution))
}

fn riccati(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    rhs: NodeId,
    coefficients: &[NodeId],
) -> Option<NodeId> {
    let (y, x) = (*problem.stand.first()?, problem.x);
    let &[_, q1, q2] = coefficients else {
        return None;
    };
    // Look for a particular solution among a few simple shapes by
    // demanding that the residual vanish identically.
    let (one, two, minus_one) = (cx.graph.int(1), cx.graph.int(2), cx.graph.int(-1));
    let shapes = [
        one,
        x,
        inv(cx.graph, x),
        cx.graph.node(core::POW, &[x, two]),
    ];
    let parameter_symbol = cx.graph.interner_mut().fresh_symbol("k");
    let parameter = cx.graph.symbol_node(parameter_symbol);
    let mut particular = None;
    for shape in shapes {
        let guess = mul(cx.graph, &[parameter, shape]);
        let slope = derivative(cx.graph, guess, x)?;
        let value = cx.graph.substitute(rhs, y, guess);
        let residual = sub(cx.graph, slope, value);
        let residual = cx.simplify(residual);
        // The residual must vanish for all x: every coefficient of its
        // numerator as a polynomial in x must, which pins the parameter.
        let numerator = numerator_of(cx.graph, residual)?;
        let Some(parts) = coefficients_in(cx.graph, numerator, x) else {
            continue;
        };
        let mut candidates: Option<Vec<NodeId>> = None;
        for part in parts {
            if cx.graph.number_of(part).is_some_and(Number::is_zero) {
                continue;
            }
            let Some(roots) = solve_for(cx.graph, part, parameter, 0) else {
                candidates = Some(Vec::new());
                break;
            };
            candidates = Some(match candidates {
                | None => roots,
                | Some(old) => old.into_iter().filter(|r| roots.iter().any(|s| cx.graph.same(*r, *s) || r == s)).collect(),
            });
        }
        if let Some(&value) = candidates.as_ref().and_then(|c| c.first()) {
            let found = cx.graph.substitute(guess, parameter, value);
            particular = Some(cx.simplify(found));
            break;
        }
    }
    let particular = particular.or_else(|| riccati_function_shapes(cx, rhs, y, x));
    let y1 = particular?;
    // y = y1 + 1/v with v' = -(q1 + 2 q2 y1) v - q2.
    let twice = mul(cx.graph, &[two, q2, y1]);
    let bracket = add(cx.graph, &[q1, twice]);
    let rate = neg(cx.graph, bracket);
    let forcing = mul(cx.graph, &[minus_one, q2]);
    let integral = integrate(cx, rate, x)?;
    let growth = call(cx.graph, "exp", integral)?;
    let decay_arg = neg(cx.graph, integral);
    let decay = call(cx.graph, "exp", decay_arg)?;
    let weighted = mul(cx.graph, &[decay, forcing]);
    let weighted = cx.simplify(weighted);
    let part = integrate(cx, weighted, x)?;
    let c = problem.constant(cx.graph);
    let sum = add(cx.graph, &[part, c]);
    let v = mul(cx.graph, &[growth, sum]);
    let correction = inv(cx.graph, v);
    let solution = add(cx.graph, &[y1, correction]);
    Some(explicit(cx, problem, solution))
}

/// A particular solution `k f(x)` of a Riccati equation among `exp(±x)`,
/// `sin x`, `cos x`, `tan x`, `1/x^2`, `x^3`, `ln x` with a small rational
/// `k`, checked by substitution.
fn riccati_function_shapes(
    cx: &mut Cx<'_>,
    rhs: NodeId,
    y: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let minus_x = neg(cx.graph, x);
    let shapes = [
        call(cx.graph, "exp", x)?,
        call(cx.graph, "exp", minus_x)?,
        call(cx.graph, "sin", x)?,
        call(cx.graph, "cos", x)?,
        call(cx.graph, "tan", x)?,
        call(cx.graph, "ln", x)?,
        {
            let m2 = cx.graph.int(-2);
            cx.graph.node(core::POW, &[x, m2])
        },
        {
            let three = cx.graph.int(3);
            cx.graph.node(core::POW, &[x, three])
        },
    ];
    for shape in shapes {
        for (num, den) in [(1, 1), (-1, 1), (2, 1), (-2, 1), (1, 2), (-1, 2), (3, 1), (-3, 1), (1, 3), (-1, 3)] {
            let k = cx.graph.num(Number::fraction(num, den)?);
            let guess = mul(cx.graph, &[k, shape]);
            let slope = derivative(cx.graph, guess, x)?;
            let value = cx.graph.substitute(rhs, y, guess);
            let residual = sub(cx.graph, slope, value);
            if cx.is_zero(residual) {
                return Some(cx.simplify(guess));
            }
        }
    }
    None
}

/// The numerator of `term` written as a single fraction.
fn numerator_of(
    graph: &mut Graph,
    term: NodeId,
) -> Option<NodeId> {
    let mut gens = Gens::default();
    let fraction = crate::rules::poly::ratio(graph, &mut gens, term, Limits::default())?;
    Some(to_term(graph, &gens, &fraction.numer))
}

/// Exact equations `M(x, y) + N(x, y) y' = 0` with `M_y = N_x`.
fn exact(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
) -> Option<NodeId> {
    let (y, dy, x) = (*problem.stand.first()?, *problem.stand.get(1)?, problem.x);
    let coefficients = coefficients_in(cx.graph, problem.expr, dy)?;
    let &[m, n] = coefficients.as_slice() else {
        return None;
    };
    // Trivially exact forms (N free of x and M free of y) are separable
    // and better handled there, where the answer is solved for y.
    if !problem.depends_on_y(cx.graph, m) || !problem.depends_on_x(cx.graph, n) {
        return None;
    }
    let m_y = derivative(cx.graph, m, y)?;
    let n_x = derivative(cx.graph, n, x)?;
    let difference = sub(cx.graph, m_y, n_x);
    if !cx.is_zero(difference) {
        return None;
    }
    // psi = ∫M dx + ∫(N - d/dy ∫M dx) dy
    let part = integrate(cx, m, x)?;
    let part_y = derivative(cx.graph, part, y)?;
    let rest = sub(cx.graph, n, part_y);
    let rest = cx.simplify(rest);
    let rest_integral = integrate(cx, rest, y)?;
    let c = problem.constant(cx.graph);
    let psi = add(cx.graph, &[part, rest_integral]);
    let relation = sub(cx.graph, psi, c);
    implicit(cx, problem, relation)
}

/// Homogeneous equations `y' = F(y/x)`: with `y = v x`, `x v' = F(v) - v`.
fn homogeneous(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    rhs: NodeId,
) -> Option<NodeId> {
    let (y, x) = (*problem.stand.first()?, problem.x);
    let v_symbol = cx.graph.interner_mut().fresh_symbol("v");
    let v = cx.graph.symbol_node(v_symbol);
    let vx = mul(cx.graph, &[v, x]);
    let in_v = cx.graph.substitute(rhs, y, vx);
    let mut in_v = cx.simplify(in_v);
    if problem.depends_on_x(cx.graph, in_v) {
        // Radicals such as sqrt(x² + v² x²) hide the cancellation: test
        // F(λx, λy) = F(x, y) numerically (λ > 0) and read F(v) at x = 1.
        let (x_symbol, y_symbol) = (problem.x_symbol, cx.graph.symbol_of(y)?);
        let at = |graph: &mut Graph, xv: f64, yv: f64| -> Option<f64> {
            let mut env = Env::numeric(0.0);
            env.bind(x_symbol, xv);
            env.bind(y_symbol, yv);
            graph.eval(rhs, &env).filter(|v| v.is_finite())
        };
        for (a, b, l) in [(0.7, 0.3, 1.9), (1.3, -0.8, 0.45), (2.1, 1.7, 3.2)] {
            let (f1, f2) = (at(cx.graph, a, b)?, at(cx.graph, l * a, l * b)?);
            if (f1 - f2).abs() > 1e-9 * (1.0 + f1.abs()) {
                return None;
            }
        }
        let one = cx.graph.int(1);
        let at_one = cx.graph.substitute(rhs, y, v);
        let at_one = cx.graph.substitute(at_one, x, one);
        in_v = cx.simplify(at_one);
    }
    // ∫ dv / (F(v) - v) = ln x + C
    let gap = sub(cx.graph, in_v, v);
    let reciprocal = inv(cx.graph, gap);
    let reciprocal = cx.simplify(reciprocal);
    let left = integrate(cx, reciprocal, v)?;
    let log = call(cx.graph, "ln", x)?;
    let c = problem.constant(cx.graph);
    let right = add(cx.graph, &[log, c]);
    let relation = sub(cx.graph, left, right);
    let ratio = div(cx.graph, y, x);
    let relation = cx.graph.substitute(relation, v, ratio);
    implicit(cx, problem, relation)
}

/// A fundamental system of `sum_k a_k y^(k) = 0` with rational constant
/// coefficients `a_k` (ascending order of derivative).
fn constant_coefficient_basis(
    graph: &mut Graph,
    coefficients: &[BigRational],
    x: NodeId,
) -> Option<Vec<NodeId>> {
    let (_, mut factors) = univariate::factor(coefficients);
    // Real roots in ascending order first, so that C1 goes with the
    // smallest rate.
    factors.sort_by(|(a, _), (b, _)| {
        let root = |f: &[BigInt]| match f {
            | [c0, c1] => Number::rat(BigRational::new(-c0.clone(), c1.clone())).to_f64(),
            | _ => f64::INFINITY,
        };
        root(a).total_cmp(&root(b))
    });
    let exp = graph.ops().lookup("exp")?;
    let (sin, cos) = (graph.ops().lookup("sin")?, graph.ops().lookup("cos")?);
    let mut basis = Vec::new();
    for (factor, multiplicity) in factors {
        // Characteristic roots of this factor and the solutions they give.
        let mut pieces: Vec<NodeId> = Vec::new();
        match factor.as_slice() {
            | [c0, c1] => {
                let root = graph.num(Number::rat(BigRational::new(-c0.clone(), c1.clone())));
                let argument = mul(graph, &[root, x]);
                pieces.push(graph.node(exp, &[argument]));
            },
            | [c, b, a] => {
                let discriminant = b * b - BigInt::from(4) * a * c;
                let two_a = BigInt::from(2) * a;
                let centre = Number::rat(BigRational::new(-b.clone(), two_a.clone()));
                let half = graph.num(Number::fraction(1, 2)?);
                let radicand = graph.num(Number::Int(discriminant.abs()));
                let root = graph.node(core::POW, &[radicand, half]);
                let scale = graph.num(Number::rat(BigRational::new(BigInt::from(1), two_a)));
                let offset = mul(graph, &[scale, root]);
                let centre_node = graph.num(centre.clone());
                if discriminant.is_negative() {
                    // exp(alpha x) cos(beta x), exp(alpha x) sin(beta x)
                    let beta_x = mul(graph, &[offset, x]);
                    let envelope = if centre.is_zero() {
                        None
                    } else {
                        let argument = mul(graph, &[centre_node, x]);
                        Some(graph.node(exp, &[argument]))
                    };
                    for wave in [cos, sin] {
                        let oscillation = graph.node(wave, &[beta_x]);
                        pieces.push(match envelope {
                            | Some(e) => mul(graph, &[e, oscillation]),
                            | None => oscillation,
                        });
                    }
                } else {
                    for sign in [-1, 1] {
                        let s = graph.int(sign);
                        let signed = mul(graph, &[s, offset]);
                        let rate = add(graph, &[centre_node, signed]);
                        let argument = mul(graph, &[rate, x]);
                        pieces.push(graph.node(exp, &[argument]));
                    }
                }
            },
            | _ => return None,
        }
        // A root of multiplicity m contributes x^j times its solution for
        // j < m.
        for j in 0..multiplicity {
            for &piece in &pieces {
                basis.push(match j {
                    | 0 => piece,
                    | 1 => mul(graph, &[x, piece]),
                    | _ => {
                        let e = graph.int(i64::from(j));
                        let power = graph.node(core::POW, &[x, e]);
                        mul(graph, &[power, piece])
                    },
                });
            }
        }
    }
    Some(basis)
}

/// Linear equations of any order: constant coefficients (with variation
/// of parameters for a right-hand side of second-order equations and an
/// integrating approach for first order handled elsewhere), and
/// Cauchy–Euler equations.
fn linear_higher_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
) -> Option<NodeId> {
    let x = problem.x;
    let order = problem.order();
    // Peel the equation apart: expr = sum_k a_k(x) y_k - r(x), after
    // clearing denominators in y and its derivatives (y''/y = c, say).
    let mut rest = cleared_denominators(cx, problem).unwrap_or(problem.expr);
    let mut coefficients = vec![NodeId::NONE; order + 1];
    for k in (0..=order).rev() {
        let parts = coefficients_in(cx.graph, rest, *problem.stand.get(k)?)?;
        let &[constant, linear] = parts.as_slice() else {
            if let &[constant] = parts.as_slice() {
                *coefficients.get_mut(k)? = cx.graph.int(0);
                rest = constant;
                continue;
            }
            return None;
        };
        *coefficients.get_mut(k)? = linear;
        rest = constant;
    }
    if problem.stand.iter().any(|&s| occurs(cx.graph, rest, s)) {
        return None;
    }
    let forcing = neg(cx.graph, rest);
    let forcing = cx.simplify(forcing);
    let homogeneous = cx.graph.number_of(forcing).is_some_and(Number::is_zero);

    let numeric: Option<Vec<BigRational>> =
        coefficients.iter().map(|&c| cx.graph.number_of(c).and_then(Number::to_rational)).collect();
    let symbolic_constant = coefficients.iter().all(|&c| !problem.depends_on_x(cx.graph, c));
    let basis = if let Some(numbers) = numeric {
        constant_coefficient_basis(cx.graph, &numbers, x)?
    } else if symbolic_constant && order <= 2 {
        symbolic_constant_basis(cx, &coefficients, x)?
    } else {
        cauchy_euler_basis(cx, &coefficients, x)?
    };
    if basis.len() != order {
        return None;
    }
    let mut terms = Vec::with_capacity(order + 1);
    for &solution in &basis {
        let c = problem.constant(cx.graph);
        terms.push(mul(cx.graph, &[c, solution]));
    }
    if !homogeneous {
        terms.push(variation_of_parameters(cx, &basis, *coefficients.last()?, forcing, x)?);
    }
    let solution = add(cx.graph, &terms);
    Some(explicit(cx, problem, solution))
}

/// The numerator of the equation when its denominator involves `y` or a
/// derivative (and is therefore not a mere coefficient).
fn cleared_denominators(
    cx: &mut Cx<'_>,
    problem: &Problem,
) -> Option<NodeId> {
    let mut gens = crate::rules::poly::repr::Gens::default();
    for &s in &problem.stand {
        gens.index(cx.graph, s);
    }
    let term = crate::rules::poly::best(cx.graph, problem.expr)?;
    let fraction = crate::rules::poly::ratio(cx.graph, &mut gens, term, crate::rules::poly::repr::Limits::default())?;
    let denominator = crate::rules::poly::repr::to_term(cx.graph, &gens, &fraction.denom);
    if !problem.stand.iter().any(|&s| occurs(cx.graph, denominator, s)) {
        return None;
    }
    let numerator = crate::rules::poly::repr::to_term(cx.graph, &gens, &fraction.numer);
    Some(cx.simplify(numerator))
}

/// `a y'' + b y' + c y = 0` (or first order) with constant symbolic
/// coefficients: `e^(r x)` for the roots `r = (-b ± √(b² - 4ac))/(2a)`,
/// assumed distinct (a repeated root needs the coefficients to satisfy a
/// relation the symbols do not show).
fn symbolic_constant_basis(
    cx: &mut Cx<'_>,
    coefficients: &[NodeId],
    x: NodeId,
) -> Option<Vec<NodeId>> {
    let exp = cx.graph.ops().lookup("exp")?;
    let roots: Vec<NodeId> = match *coefficients {
        | [c0, c1] => {
            let r = div(cx.graph, c0, c1);
            vec![neg(cx.graph, r)]
        },
        | [c, b, a] => {
            let two = cx.graph.int(2);
            let four = cx.graph.int(4);
            let b2 = cx.graph.node(core::POW, &[b, two]);
            let fac = mul(cx.graph, &[four, a, c]);
            let discriminant = sub(cx.graph, b2, fac);
            let discriminant = cx.simplify(discriminant);
            if cx.is_zero(discriminant) {
                return None;
            }
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let root = cx.graph.node(core::POW, &[discriminant, half]);
            let minus_b = neg(cx.graph, b);
            let two_a = mul(cx.graph, &[two, a]);
            let mut out = Vec::new();
            for sign in [-1, 1] {
                let s = cx.graph.int(sign);
                let signed = mul(cx.graph, &[s, root]);
                let numerator = add(cx.graph, &[minus_b, signed]);
                out.push(div(cx.graph, numerator, two_a));
            }
            out
        },
        | _ => return None,
    };
    Some(
        roots
            .into_iter()
            .map(|r| {
                let r = cx.simplify(r);
                let rx = mul(cx.graph, &[r, x]);
                cx.graph.node(exp, &[rx])
            })
            .collect(),
    )
}

/// Cauchy–Euler equations `sum_k c_k x^k y^(k) = 0`: `y = x^m` with `m` a
/// root of the indicial polynomial.
fn cauchy_euler_basis(
    cx: &mut Cx<'_>,
    coefficients: &[NodeId],
    x: NodeId,
) -> Option<Vec<NodeId>> {
    // a_k(x) = c_k x^k with rational c_k.
    let mut constants = Vec::with_capacity(coefficients.len());
    for (k, &a) in coefficients.iter().enumerate() {
        let e = cx.graph.int(-i64::try_from(k).ok()?);
        let scale = cx.graph.node(core::POW, &[x, e]);
        let quotient = mul(cx.graph, &[a, scale]);
        let quotient = cx.simplify(quotient);
        constants.push(cx.graph.number_of(quotient)?.to_rational()?);
    }
    // Indicial polynomial: sum_k c_k m (m-1) ... (m-k+1).
    let mut indicial: Vec<BigRational> = vec![BigRational::zero()];
    for (k, c) in constants.iter().enumerate() {
        let mut falling: Vec<BigRational> = vec![BigRational::from_integer(BigInt::from(1))];
        for j in 0..k {
            let shift = [BigRational::from_integer(BigInt::from(-i64::try_from(j).ok()?)), BigRational::from_integer(BigInt::from(1))];
            falling = univariate::mul(&falling, &shift);
        }
        let scaled: Vec<BigRational> = falling.iter().map(|v| v * c).collect();
        indicial = univariate::add(&indicial, &scaled);
    }
    let (_, factors) = univariate::factor(&indicial);
    let ln = cx.graph.ops().lookup("ln")?;
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    let graph = &mut *cx.graph;
    let log = graph.node(ln, &[x]);
    let mut basis = Vec::new();
    for (factor, multiplicity) in factors {
        let mut pieces = Vec::new();
        match factor.as_slice() {
            | [c0, c1] => {
                let m = graph.num(Number::rat(BigRational::new(-c0.clone(), c1.clone())));
                pieces.push(graph.node(core::POW, &[x, m]));
            },
            | [c, b, a] => {
                let discriminant = b * b - BigInt::from(4) * a * c;
                let two_a = BigInt::from(2) * a;
                let centre = graph.num(Number::rat(BigRational::new(-b.clone(), two_a.clone())));
                let half = graph.num(Number::fraction(1, 2)?);
                let radicand = graph.num(Number::Int(discriminant.abs()));
                let root = graph.node(core::POW, &[radicand, half]);
                let scale = graph.num(Number::rat(BigRational::new(BigInt::from(1), two_a)));
                let offset = mul(graph, &[scale, root]);
                if discriminant.is_negative() {
                    let envelope = graph.node(core::POW, &[x, centre]);
                    let angle = mul(graph, &[offset, log]);
                    for wave in [cos, sin] {
                        let oscillation = graph.node(wave, &[angle]);
                        pieces.push(mul(graph, &[envelope, oscillation]));
                    }
                } else {
                    for sign in [-1, 1] {
                        let s = graph.int(sign);
                        let signed = mul(graph, &[s, offset]);
                        let m = add(graph, &[centre, signed]);
                        pieces.push(graph.node(core::POW, &[x, m]));
                    }
                }
            },
            | _ => return None,
        }
        for j in 0..multiplicity {
            for &piece in &pieces {
                basis.push(match j {
                    | 0 => piece,
                    | _ => {
                        let e = graph.int(i64::from(j));
                        let power = graph.node(core::POW, &[log, e]);
                        mul(graph, &[power, piece])
                    },
                });
            }
        }
    }
    Some(basis)
}

/// A particular solution of `a_n y^(n) + ... = r` from a fundamental
/// system (any order up to five).
fn variation_of_parameters(
    cx: &mut Cx<'_>,
    basis: &[NodeId],
    leading: NodeId,
    forcing: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let normalised = div(cx.graph, forcing, leading);
    match *basis {
        | [y1] => {
            // y_p = y1 ∫ r / y1
            let integrand = div(cx.graph, normalised, y1);
            let integrand = cx.simplify(integrand);
            let integral = integrate(cx, integrand, x)?;
            Some(mul(cx.graph, &[y1, integral]))
        },
        | [y1, y2] => {
            let (d1, d2) = (derivative(cx.graph, y1, x)?, derivative(cx.graph, y2, x)?);
            let first = mul(cx.graph, &[y1, d2]);
            let second = mul(cx.graph, &[y2, d1]);
            let wronskian = sub(cx.graph, first, second);
            let wronskian = expanded(cx, wronskian);
            // y_p = -y1 ∫ y2 r / W + y2 ∫ y1 r / W
            let inverse = inv(cx.graph, wronskian);
            let u1 = mul(cx.graph, &[y2, normalised, inverse]);
            let u1 = cx.simplify(u1);
            let u2 = mul(cx.graph, &[y1, normalised, inverse]);
            let u2 = cx.simplify(u2);
            let (i1, i2) = (integrate(cx, u1, x)?, integrate(cx, u2, x)?);
            let left = mul(cx.graph, &[y1, i1]);
            let left = neg(cx.graph, left);
            let right = mul(cx.graph, &[y2, i2]);
            Some(add(cx.graph, &[left, right]))
        },
        | _ if basis.len() <= 5 => {
            // Cramer's rule on the Wronskian system W u' = (0, ..., 0, r):
            // u_i' = r C_(n-1, i) / W with C the cofactors of the last row.
            let n = basis.len();
            let mut rows: Vec<Vec<NodeId>> = vec![basis.to_vec()];
            for _ in 1..n {
                let last = rows.last()?.clone();
                let next: Vec<NodeId> = last
                    .iter()
                    .map(|&f| {
                        let d = derivative(cx.graph, f, x)?;
                        Some(cx.simplify(d))
                    })
                    .collect::<Option<_>>()?;
                rows.push(next);
            }
            let wronskian = determinant_of(cx, &rows)?;
            let wronskian = expanded(cx, wronskian);
            if cx.is_zero(wronskian) {
                return None;
            }
            let inverse = inv(cx.graph, wronskian);
            let mut terms = Vec::with_capacity(n);
            for (i, &yi) in basis.iter().enumerate() {
                let minor: Vec<Vec<NodeId>> =
                    rows[..n - 1].iter().map(|r| r.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).collect()).collect();
                let cofactor = determinant_of(cx, &minor)?;
                let sign = cx.graph.int(if (n - 1 + i).is_multiple_of(2) { 1 } else { -1 });
                let ui = mul(cx.graph, &[sign, cofactor, normalised, inverse]);
                let ui = cx.simplify(ui);
                let ui = crate::rules::poly::expand_form(cx.graph, ui).unwrap_or(ui);
                let ui = cx.simplify(ui);
                let integral = integrate(cx, ui, x)?;
                terms.push(mul(cx.graph, &[yi, integral]));
            }
            Some(add(cx.graph, &terms))
        },
        | _ => None,
    }
}

/// `e` simplified, expanded and simplified again (Wronskians such as
/// `(x e^x + e^x) e^x - x e^(2x)` collapse only once expanded).
fn expanded(
    cx: &mut Cx<'_>,
    e: NodeId,
) -> NodeId {
    let e = cx.simplify(e);
    let e = crate::rules::poly::expand_form(cx.graph, e).unwrap_or(e);
    cx.simplify(e)
}

/// The determinant of a small square matrix of terms, by cofactor
/// expansion along the first row.
fn determinant_of(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    match m.len() {
        | 0 => Some(cx.graph.int(1)),
        | 1 => m[0].first().copied(),
        | n => {
            let mut terms = Vec::with_capacity(n);
            for (j, &entry) in m[0].iter().enumerate() {
                let minor: Vec<Vec<NodeId>> =
                    m[1..].iter().map(|r| r.iter().enumerate().filter(|&(k, _)| k != j).map(|(_, &v)| v).collect()).collect();
                let sub_det = determinant_of(cx, &minor)?;
                let sign = cx.graph.int(if j.is_multiple_of(2) { 1 } else { -1 });
                terms.push(mul(cx.graph, &[sign, entry, sub_det]));
            }
            let sum = add(cx.graph, &terms);
            Some(cx.simplify(sum))
        },
    }
}

/// Substitutes `solution` for `y` in the original equation and checks the
/// residual numerically at a few points.
fn verified(
    cx: &mut Cx<'_>,
    problem: &Problem,
    answer: NodeId,
) -> bool {
    let graph = &mut *cx.graph;
    if problem.trusted {
        return true;
    }
    let &[lhs, solution] = graph.children(answer) else {
        return true;
    };
    if lhs != problem.y_of_x {
        // An implicit relation: not checked here.
        return true;
    }
    // Derivatives of the solution, then the equation with everything
    // substituted.
    let mut residual = problem.expr;
    let mut current = solution;
    for (k, &stand) in problem.stand.iter().enumerate() {
        if k > 0 {
            let Some(next) = derivative(graph, current, problem.x) else {
                return true;
            };
            current = next;
        }
        residual = graph.substitute(residual, stand, current);
    }
    let symbols: Vec<SymbolId> = graph.free_symbols(graph.find(residual)).to_vec();
    for point in [0.41, 0.83, 1.37] {
        let mut env = Env::numeric(0.0);
        for &symbol in &symbols {
            let value = if symbol == problem.x_symbol {
                point
            } else {
                0.6 + 0.17 * f64::from(symbol.raw() % 5)
            };
            env.bind(symbol, value);
        }
        if let Some(value) = graph.eval(residual, &env) {
            // Points where the residual cannot be evaluated say nothing.
            if value.is_finite() && value.abs() > 1e-6 {
                return false;
            }
        }
    }
    true
}

/// Determines the integration constants from a list of conditions.
///
/// A condition is `y(a) = b`, or `diff(y(x), x) = b` (any order), which is
/// taken at the point of the condition before it.
fn apply_conditions(
    cx: &mut Cx<'_>,
    problem: &Problem,
    answer: NodeId,
    conditions: NodeId,
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let &[lhs, solution] = graph.children(answer) else {
        return None;
    };
    if lhs != problem.y_of_x || graph.op(conditions) != core::LIST {
        return None;
    }
    let function = *graph.children(problem.y_of_x).first()?;
    let diff = graph.ops().lookup("diff")?;
    let mut equations = Vec::new();
    let mut last_point = None;
    for condition in graph.children(conditions).to_vec() {
        let &[target, value] = graph.children(condition) else {
            return None;
        };
        if graph.op(condition) != core::EQ {
            return None;
        }
        // y(a) = b, or diff(y(x), x) = b taken at the last point given.
        let (order, point) = match *graph.children(target) {
            | [f, a] if graph.op(target) == core::APPLY && f == function => (0, a),
            | _ => {
                let mut order = 0;
                let mut inner = target;
                while graph.op(inner) == diff {
                    inner = *graph.children(inner).first()?;
                    order += 1;
                }
                if inner != problem.y_of_x || order == 0 {
                    return None;
                }
                (order, last_point?)
            },
        };
        last_point = Some(point);
        let mut expression = solution;
        for _ in 0..order {
            expression = derivative(graph, expression, problem.x)?;
        }
        let at_point = graph.substitute(expression, problem.x, point);
        equations.push(sub(graph, at_point, value));
    }
    let constants: Vec<NodeId> = (1..=problem.constants).map(|k| graph.sym(&format!("C{k}"))).collect();
    if equations.len() != constants.len() {
        return None;
    }
    let simplified: Vec<NodeId> = equations.into_iter().map(|e| cx.simplify(e)).collect();
    let values = match solve_linear(cx.graph, &simplified, &constants) {
        | Some(values) => values,
        | None => {
            // One constant entering non-linearly.
            let (&equation, &constant) = (simplified.first()?, constants.first()?);
            if constants.len() != 1 {
                return None;
            }
            vec![*solve_for(cx.graph, equation, constant, 0)?.first()?]
        },
    };
    let mut result = solution;
    for (&constant, &value) in constants.iter().zip(&values) {
        result = cx.graph.substitute(result, constant, value);
    }
    Some(explicit(cx, problem, result))
}

/// The coefficients of `expr` as a polynomial in the symbol `v` (free of
/// `v`), for other rule sets.
pub(crate) fn coefficients_of(
    graph: &mut Graph,
    expr: NodeId,
    v: NodeId,
) -> Option<Vec<NodeId>> {
    coefficients_in(graph, expr, v)
}

/// Whether `needle` occurs in `term`.
pub(crate) fn occurs_in(
    graph: &Graph,
    term: NodeId,
    needle: NodeId,
) -> bool {
    occurs(graph, term, needle)
}

/// The solutions of `expr = 0` for the symbol `v`, for other rule sets.
pub(crate) fn solve_in(
    graph: &mut Graph,
    expr: NodeId,
    v: NodeId,
) -> Option<Vec<NodeId>> {
    solve_for(graph, expr, v, 0)
}

/// Whether some node of `term` has operator `op`.
pub(crate) fn occurs_op(
    graph: &Graph,
    term: NodeId,
    op: OpId,
) -> bool {
    let mut stack = vec![term];
    while let Some(n) = stack.pop() {
        if graph.op(n) == op {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// The solutions of a homogeneous linear identity in `unknowns` that must
/// hold for all values of every other generator (see the Lie module).
pub(crate) fn solve_linear_identity(
    graph: &mut Graph,
    expr: NodeId,
    unknowns: &[NodeId],
) -> Option<Vec<Vec<BigRational>>> {
    lie::solve_identity(graph, expr, unknowns)
}

/// The general solution of an ordinary differential equation for
/// `unknown` (`y(x) = …` or an implicit relation), for other rule sets.
pub(crate) fn solve_ode(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
) -> Option<NodeId> {
    solve_equation(cx, equation, unknown, 0).map(|(_, answer)| answer)
}

/// Parses and solves one equation for `unknown`: the general solution
/// (explicit or implicit), verified against the equation. Reduction
/// methods call this again on the reduced equation, so every method is
/// available at every stage.
pub(super) fn solve_equation(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    depth: u32,
) -> Option<(Problem, NodeId)> {
    if depth > 3 {
        return None;
    }
    let mut problem = parse(cx.graph, equation, unknown)?;
    let answer = dispatch(cx, &mut problem, depth)?;
    Some((problem, answer))
}

/// The method pipeline, cheapest and most specific first.
fn dispatch(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    // Every candidate is checked; a failed one does not stop the search.
    let attempt = |cx: &mut Cx<'_>, problem: &mut Problem, method: fn(&mut Cx<'_>, &mut Problem, u32) -> Option<NodeId>| {
        let saved = problem.constants;
        let found = method(cx, problem, depth).filter(|&answer| verified(cx, problem, answer));
        if found.is_none() {
            problem.constants = saved;
        }
        found
    };
    if problem.order() == 1 {
        if let Some(found) = attempt(cx, problem, first_order) {
            return Some(found);
        }
    }
    for method in [
        (|cx: &mut Cx<'_>, p: &mut Problem, _| linear_higher_order(cx, p)) as fn(&mut Cx<'_>, &mut Problem, u32) -> Option<NodeId>,
        reduce::special_function_equation,
        reduce::variable_coefficients,
        reduce::missing_dependent,
        reduce::autonomous,
        reduce::scale_invariant,
        reduce::equidimensional_in_x,
        lie::second_order,
    ] {
        if let Some(found) = attempt(cx, problem, method) {
            return Some(found);
        }
    }
    None
}

struct Dsolve {
    dsolve: OpId,
}

impl Kernel for Dsolve {
    fn ops(&self) -> Vec<OpId> {
        vec![self.dsolve]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let (equation, unknown, conditions) = match *args.as_slice() {
            | [e, u] => (e, u, None),
            | [e, u, c] => (e, u, Some(c)),
            | _ => return Outcome::Pass,
        };
        if cx.graph.op(unknown) == core::LIST {
            return systems::solve_system(cx, equation, unknown)
                .or_else(|| nonlinear::solve(cx, equation, unknown, 0))
                .map_or(Outcome::Pass, Outcome::Pinned);
        }
        let Some((problem, answer)) = solve_equation(cx, equation, unknown, 0) else {
            // Discontinuous forcing: the Laplace transform.
            return conditions
                .and_then(|c| series_method::laplace_ivp(cx, equation, unknown, c))
                .map_or(Outcome::Pass, Outcome::Pinned);
        };
        match conditions {
            | None => Outcome::Pinned(answer),
            | Some(conditions) => {
                apply_conditions(cx, &problem, answer, conditions).map_or(Outcome::Pass, Outcome::Pinned)
            },
        }
    }
}

/// `ode_series`, `rsolve_z`, `ode_classify`, `ode_symmetries` and
/// `ode_fourier`.
struct Extras {
    series: OpId,
    rsolve: OpId,
    classify: OpId,
    symmetries: OpId,
    fourier: OpId,
}

impl Kernel for Extras {
    fn ops(&self) -> Vec<OpId> {
        vec![self.series, self.rsolve, self.classify, self.symmetries, self.fourier]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let args = cx.graph.children(node).to_vec();
        let found = if op == self.series {
            match *args.as_slice() {
                | [e, u, x0, n] => series_method::series_solution(cx, e, u, x0, n),
                | _ => None,
            }
        } else if op == self.rsolve {
            match *args.as_slice() {
                | [e, u] => series_method::rsolve(cx, e, u, None),
                | [e, u, c] => series_method::rsolve(cx, e, u, Some(c)),
                | _ => None,
            }
        } else if op == self.fourier {
            match *args.as_slice() {
                | [e, u] => series_method::fourier_solution(cx, e, u),
                | _ => None,
            }
        } else if op == self.symmetries {
            match *args.as_slice() {
                | [e, u] => symmetry_list(cx, e, u),
                | _ => None,
            }
        } else {
            match *args.as_slice() {
                | [e, u] => classify(cx, e, u),
                | _ => None,
            }
        };
        found.map_or(Outcome::Pass, Outcome::Pinned)
    }
}

/// `ode_symmetries(eq, y(x))`: point symmetries `list(ξ, η)` with
/// polynomial infinitesimals (degree ≤ 2 in `x` and `y`).
fn symmetry_list(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
) -> Option<NodeId> {
    let problem = parse(cx.graph, equation, unknown)?;
    let order = problem.order();
    let top = problem.stand[order];
    let omega = *solve_for(cx.graph, problem.expr, top, 0)?.first()?;
    let omega = cx.simplify(omega);
    let pairs = if order == 1 {
        lie::first_order_symmetries(cx, problem.x, problem.stand[0], omega, 2)
    } else {
        lie::higher_order_symmetries(cx, &problem, omega)?
    };
    let items: Vec<NodeId> = pairs
        .into_iter()
        .map(|(xi, eta)| {
            let eta = cx.graph.substitute(eta, problem.stand[0], unknown);
            let xi = cx.graph.substitute(xi, problem.stand[0], unknown);
            cx.graph.node(core::LIST, &[xi, eta])
        })
        .collect();
    Some(cx.graph.node(core::LIST, &items))
}

/// The classes an equation belongs to, as a list of symbols, in the
/// order the solver tries them: `order_n`, `linear`, `constant_coefficients`,
/// `cauchy_euler`, `separable`, `exact`, `homogeneous`, `bernoulli`,
/// `riccati`, `autonomous`, `missing_y`, `scale_invariant`, `lie_symmetry`.
fn classify(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
) -> Option<NodeId> {
    let mut problem = parse(cx.graph, equation, unknown)?;
    let order = problem.order();
    let mut classes = vec![format!("order_{order}")];
    let x = problem.x;
    let y = problem.stand[0];
    // Linear: the equation has degree one in y and its derivatives jointly.
    let mut linear = true;
    let mut constant = true;
    let mut rest = problem.expr;
    for k in (0..=order).rev() {
        match coefficients_in(cx.graph, rest, problem.stand[k]).as_deref() {
            | Some(&[c0, c1]) => {
                if problem.stand.iter().any(|&s| occurs(cx.graph, c1, s)) {
                    linear = false;
                }
                if problem.depends_on_x(cx.graph, c1) {
                    constant = false;
                }
                rest = c0;
            },
            | Some(&[c0]) => rest = c0,
            | _ => {
                linear = false;
                break;
            },
        }
    }
    if linear {
        classes.push("linear".to_owned());
        if constant {
            classes.push("constant_coefficients".to_owned());
        }
    }
    if order == 1 {
        let dy = problem.stand[1];
        if let Some(rhs) = solve_for(cx.graph, problem.expr, dy, 0).and_then(|v| v.first().copied()) {
            let rhs = cx.simplify(rhs);
            if let Some(c) = coefficients_in(cx.graph, rhs, y) {
                match c.len() {
                    | 2 => {},
                    | 3 => classes.push("riccati".to_owned()),
                    | n if n > 3 => classes.push("bernoulli".to_owned()),
                    | _ => {},
                }
            }
            let saved = problem.constants;
            if separable(cx, &mut problem, rhs).is_some() {
                classes.push("separable".to_owned());
            }
            if homogeneous(cx, &mut problem, rhs).is_some() {
                classes.push("homogeneous".to_owned());
            }
            if exact(cx, &mut problem).is_some() {
                classes.push("exact".to_owned());
            }
            if !lie::first_order_symmetries(cx, x, y, rhs, 2).is_empty() {
                classes.push("lie_symmetry".to_owned());
            }
            problem.constants = saved;
        }
    } else {
        if !problem.depends_on_x(cx.graph, problem.expr) {
            classes.push("autonomous".to_owned());
        }
        if !occurs(cx.graph, problem.expr, y) {
            classes.push("missing_y".to_owned());
        }
    }
    let items: Vec<NodeId> = classes.iter().map(|c| cx.graph.sym(c)).collect();
    Some(cx.graph.node(core::LIST, &items))
}

/// Numeric initial value problems.
struct Odeint {
    odeint: OpId,
}

impl Kernel for Odeint {
    fn ops(&self) -> Vec<OpId> {
        vec![self.odeint]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        // Already evaluated in this run.
        if cx.graph.approx(cx.graph.find(node)).is_some() {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let &[f, y, x, x0, y0, x1] = graph.children(node) else {
            return Outcome::Pass;
        };
        let items = |graph: &Graph, n: NodeId| -> Vec<NodeId> {
            if graph.op(n) == core::LIST { graph.children(n).to_vec() } else { vec![n] }
        };
        let is_system = graph.op(y) == core::LIST;
        let (rhs, unknowns, initial) = (items(graph, f), items(graph, y), items(graph, y0));
        if rhs.len() != unknowns.len() || rhs.len() != initial.len() || rhs.is_empty() {
            return Outcome::Pass;
        }
        let Some(x_symbol) = graph.symbol_of(x) else {
            return Outcome::Pass;
        };
        let Some(state_symbols) = unknowns.iter().map(|&u| graph.symbol_of(u)).collect::<Option<Vec<_>>>() else {
            return Outcome::Pass;
        };
        let value = |graph: &mut Graph, n: NodeId| best(graph, n).and_then(|t| graph.eval(t, cx.env));
        let (Some(t0), Some(t1)) = (value(graph, x0), value(graph, x1)) else {
            return Outcome::Pass;
        };
        let Some(start) = initial.iter().map(|&n| value(graph, n)).collect::<Option<Vec<f64>>>() else {
            return Outcome::Pass;
        };
        // Inputs: x, the state, then bound parameters.
        let mut inputs = vec![x_symbol];
        inputs.extend_from_slice(&state_symbols);
        let mut parameters = Vec::new();
        for &(symbol, v) in cx.env.bindings() {
            if !inputs.contains(&symbol) {
                inputs.push(symbol);
                parameters.push(v);
            }
        }
        let mut compiled = Vec::with_capacity(rhs.len());
        for &expression in &rhs {
            let Some(term) = best(graph, expression) else {
                return Outcome::Pass;
            };
            let Ok(function) = Interpreter.compile(graph, term, &inputs) else {
                return Outcome::Pass;
            };
            compiled.push(function);
        }
        let n = rhs.len();
        let derivative = |t: f64, state: &[f64], out: &mut [f64]| {
            let mut args = Vec::with_capacity(1 + n + parameters.len());
            args.push(t);
            args.extend_from_slice(state);
            args.extend_from_slice(&parameters);
            for (slot, function) in out.iter_mut().zip(&compiled) {
                *slot = function.call(&args);
            }
        };
        let tolerance = cx.env.tolerance.max(1e-12);
        let Ok(trajectory) = solve_adaptive(derivative, &start, (t0, t1), tolerance, tolerance * 1e-2, 2_000_000) else {
            return Outcome::Pass;
        };
        let last = trajectory.last();
        if is_system {
            let nodes: Vec<NodeId> = last.iter().map(|&v| graph.float(v)).collect();
            Outcome::Equal(graph.node(core::LIST, &nodes))
        } else {
            match last.first() {
                // The error estimate is the local tolerance accumulated
                // over the steps taken: an estimate, not a bound.
                | Some(&v) => {
                    #[allow(clippy::cast_precision_loss)]
                    let steps = trajectory.t.len() as f64;
                    Outcome::Approx(Ball { mid: v, rad: tolerance * steps.max(1.0) * (1.0 + v.abs()) })
                },
                | None => Outcome::Pass,
            }
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[ode()], src)
    }

    /// Solves `equation` for `y(x)` and checks the answer independently:
    /// the right-hand side of `y(x) = ...` is evaluated with finite
    /// differences for the derivatives and put back into the equation.
    fn check(
        equation: &str,
        order: usize,
    ) -> String {
        let (answer, reduced) = reduce_with(&[ode()], &format!("dsolve({equation}, y(x))"), &[]);
        assert!(reduced, "{equation} was not solved: {answer}");
        let solution = answer.strip_prefix("y(x) = ").unwrap_or_else(|| panic!("implicit answer: {answer}"));
        let constants = [("C1", 0.8), ("C2", -0.6), ("C3", 0.4)];
        let y = |x: f64| {
            let mut bindings = vec![("x", x)];
            bindings.extend_from_slice(&constants);
            eval(&[ode()], solution, &bindings)
        };
        let (lhs, rhs) = equation.split_once('=').unwrap_or((equation, "0"));
        let mut checked = 0;
        for at in [0.6, 1.1, 1.7, 0.3] {
            // A smaller step where only low derivatives are needed: the
            // solutions include poles nearby.
            let h = if order >= 3 { 1e-3 } else { 1e-4 };
            let y0 = y(at);
            let y1 = (y(at + h) - y(at - h)) / (2.0 * h);
            let y2 = (y(at + h) - 2.0 * y0 + y(at - h)) / (h * h);
            let y3 = (y(at + 2.0 * h) - 2.0 * y(at + h) + 2.0 * y(at - h) - y(at - 2.0 * h)) / (2.0 * h * h * h);
            // Replace derivatives innermost-last so that the longer
            // spellings go first.
            let substitute = |text: &str| {
                text.replace("diff(diff(diff(y(x), x), x), x)", &format!("({y3})"))
                    .replace("diff(diff(y(x), x), x)", &format!("({y2})"))
                    .replace("diff(y(x), x)", &format!("({y1})"))
                    .replace("y(x)", &format!("({y0})"))
            };
            let residual = eval(&[ode()], &format!("({}) - ({})", substitute(lhs), substitute(rhs)), &[("x", at)]);
            // Outside the solution's domain (a square root of a negative
            // number for these constants): nothing to check here.
            if !residual.is_finite() {
                continue;
            }
            checked += 1;
            let tolerance = if order >= 3 { 1e-2 } else { 1e-3 };
            assert!(residual.abs() < tolerance * (1.0 + y0.abs()), "{equation}: y = {solution} leaves residual {residual} at {at}");
        }
        assert!(checked > 0, "{equation}: y = {solution} could not be evaluated anywhere");
        answer
    }

    #[test]
    fn first_order_linear_and_separable() {
        assert_eq!(run("dsolve(diff(y(x), x) = y(x), y(x))"), "y(x) = C1*exp(x)");
        check("diff(y(x), x) + 2*y(x) = 0", 1);
        check("diff(y(x), x) + y(x) = x", 1);
        check("diff(y(x), x) = x*y(x)", 1);
        check("x*diff(y(x), x) + y(x) = x^2", 1);
        check("diff(y(x), x) = cos(x) - y(x)", 1);
        check("diff(y(x), x) = y(x)^2", 1);
        check("diff(y(x), x) = x/y(x)", 1);
        check("diff(y(x), x) = exp(x - y(x))", 1);
    }

    #[test]
    fn bernoulli_riccati_homogeneous_exact() {
        check("diff(y(x), x) + y(x) = y(x)^2", 1);
        check("diff(y(x), x) - y(x)/x = x*y(x)^3", 1);
        check("diff(y(x), x) = 1 + y(x)^2 - 2*y(x)", 1);
        check("diff(y(x), x) = (x + y(x))/x", 1);
        check("(2*x*y(x) + 1) + (x^2 + 2*y(x))*diff(y(x), x) = 0", 1);
    }

    #[test]
    fn lie_symmetries() {
        // y' = F(x + y) (translation), linear fractional, and equations
        // with a scaling symmetry beyond the classical recipes.
        check("diff(y(x), x) = (x + y(x))^2", 1);
        // Linear fractional: an implicit solution (log + atan about the
        // centre of the scaling symmetry).
        let (text, reduced) = reduce_with(&[ode()], "dsolve(diff(y(x), x) = (x + y(x) + 1)/(x - y(x) + 3), y(x))", &[]);
        assert!(reduced && text.contains("atan") && text.contains("C1"), "{text}");
        check("diff(y(x), x) = y(x)/x + x^2/y(x)", 1);
        check("diff(y(x), x) = y(x)/(x + y(x)^2)", 1);
    }

    #[test]
    fn reductions_of_order() {
        // y missing: p = y'.
        check("x*diff(diff(y(x), x), x) + diff(y(x), x) = 0", 2);
        check("diff(diff(y(x), x), x) = diff(y(x), x)^2", 2);
        // Autonomous: y'' = p dp/dy.
        check("diff(diff(y(x), x), x) = 2*y(x)*diff(y(x), x)", 2);
        // Scale invariant in y: u = y'/y.
        check("y(x)*diff(diff(y(x), x), x) - diff(y(x), x)^2 = 0", 2);
        // Equidimensional in x: x = e^t makes it autonomous.
        let (text, reduced) = reduce_with(&[ode()], "dsolve(x^2*diff(diff(y(x), x), x) = x*y(x)*diff(y(x), x), y(x))", &[]);
        assert!(reduced && text.contains("C1"), "{text}");
        // Variable coefficients with a polynomial / exponential solution.
        check("x*diff(diff(y(x), x), x) - (x + 1)*diff(y(x), x) + y(x) = 0", 2);
        check("(1 - x^2)*diff(diff(y(x), x), x) - 2*x*diff(y(x), x) + 2*y(x) = 0", 2);
    }

    #[test]
    fn special_function_equations() {
        let run = |src: &str| simplify(&crate::rules::standard(), src);
        let solved = run("dsolve(x^2*diff(diff(y(x), x), x) + x*diff(y(x), x) + (x^2 - 4)*y(x) = 0, y(x))");
        assert_eq!(solved, "y(x) = C1*besselj(2, x) + C2*bessely(2, x)");
        let modified = run("dsolve(x^2*diff(diff(y(x), x), x) + x*diff(y(x), x) - (9*x^2 + 1)*y(x) = 0, y(x))");
        assert_eq!(modified, "y(x) = C1*besseli(1, 3*x) + C2*besselk(1, 3*x)");
        // exp(x²) by the exp(λx²) family, the second solution through erf.
        let gaussian = run("dsolve(diff(diff(y(x), x), x) - 2*x*diff(y(x), x) - 2*y(x) = 0, y(x))");
        assert!(gaussian.contains("exp(x^2)") && gaussian.contains("erf(x)"), "{gaussian}");
    }

    #[test]
    fn linear_systems() {
        let rotation = run("dsolve(list(diff(u(t), t) = w(t), diff(w(t), t) = -u(t)), list(u(t), w(t)))");
        assert!(rotation.contains("cos(t)") && rotation.contains("sin(t)"), "{rotation}");
        let (text, reduced) = reduce_with(
            &crate::rules::standard(),
            "dsolve(list(diff(u(t), t) = u(t) + 2*w(t) + t, diff(w(t), t) = 3*u(t) + 2*w(t)), list(u(t), w(t)))",
            &[],
        );
        assert!(reduced && text.contains("exp(4*t)") && text.contains("C2"), "{text}");
        let three = run("dsolve(list(diff(a(t), t) = b(t), diff(b(t), t) = c(t), diff(c(t), t) = a(t)), list(a(t), b(t), c(t)))");
        assert!(three.contains("exp(t)") && three.contains("C3"), "{three}");
    }

    #[test]
    fn lie_point_symmetries_of_higher_order() {
        // y'' = 0 has the eight-dimensional projective algebra sl(3).
        let mut g = Graph::new();
        let rules = crate::rules::standard();
        let _ = &mut g;
        let count = simplify(&rules, "ode_symmetries(diff(diff(y(x), x), x) = 0, y(x))");
        assert_eq!(count.matches("list(").count() - 1, 8, "{count}");
        // y'' = y'^2/y + y'/x has the scaling x ∂x and y ∂y; Lie reduction
        // through canonical coordinates solves it.
        let (text, reduced) = reduce_with(&rules, "dsolve(diff(diff(y(x), x), x) = diff(y(x), x)^2/y(x) - diff(y(x), x)/x, y(x))", &[]);
        assert!(reduced && text.contains("C1"), "{text}");
        // Only the symmetry ∂x + ∂y reduces y'' = (y - x)^(-3): canonical
        // coordinates r = y - x, s = x.
        let (text, reduced) = reduce_with(&rules, "dsolve(diff(diff(y(x), x), x) = (y(x) - x)^(-3), y(x))", &[]);
        assert!(reduced && text.contains("C2"), "{text}");
    }

    #[test]
    fn harder_equations() {
        let rules = crate::rules::standard();
        let d = "diff(y(x), x)";
        let dd = "diff(diff(y(x), x), x)";
        for (eq, contains) in [
            (format!("{d} + y(x)*tan(x) = sec(x)"), "cos(x)"),
            (format!("{dd} - 2*{d} + y(x) = exp(x)/x"), "ln(x)"),
            (format!("{dd} = x*y(x)"), "airyai(x)"),
            (format!("diff({dd}, x) - 3*{dd} + 3*{d} - y(x) = exp(x)"), "x^3*exp(x)"),
            (format!("x*{d} = y(x) + sqrt(x^2 + y(x)^2)"), "C1"),
            (format!("{dd} = 2*y(x)^3"), "integral("),
            // The pendulum: energy integral, the quadrature left inert.
            (format!("{dd} + sin(y(x)) = 0"), "cos(y(x))"),
            (format!("{dd} = lambda*y(x)"), "exp(lambda^(1/2)*x)"),
        ] {
            let got = simplify(&rules, &format!("dsolve({eq}, y(x))"));
            assert!(got.contains(contains) && !got.starts_with("dsolve"), "{eq}: {got}");
        }
    }

    #[test]
    fn series_and_transform_methods() {
        let rules = crate::rules::standard();
        // Airy-type y'' = x y about 0: 1 + x^3/6 + … and x + x^4/12 + ….
        let series = simplify(&rules, "ode_series(diff(diff(y(x), x), x) = x*y(x), y(x), 0, 5)");
        assert_eq!(series, "y(x) = 1/12*C2*x^4 + 1/6*C1*x^3 + C2*x + C1");
        // A non-linear equation: y' = 1 + y^2, y = tan(x) + …
        let tangent = simplify(&rules, "ode_series(diff(y(x), x) = 1 + y(x)^2, y(x), 0, 5)");
        assert!(tangent.contains("x^3*(C1^4 + 4/3*C1^2 + 1/3)"), "{tangent}");
        // Discontinuous forcing through the Laplace transform.
        let (text, reduced) = reduce_with(
            &rules,
            "dsolve(diff(y(t), t) + y(t) = heaviside(t - 1), y(t), list(y(0) = 0))",
            &[],
        );
        assert!(reduced && text.contains("heaviside(t - 1)"), "{text}");
        // Recurrences by the z-transform: Fibonacci and a forced one.
        let fib = simplify(&rules, "rsolve_z(y(n + 2) = y(n + 1) + y(n), y(n), list(y(0) = 0, y(1) = 1))");
        for (n, want) in [(5.0, 5.0), (10.0, 55.0)] {
            let got = eval(&rules, fib.trim_start_matches("y(n) = "), &[("n", n)]);
            assert!((got - want).abs() < 1e-6, "{fib}");
        }
        let geometric = simplify(&rules, "rsolve_z(y(n + 1) = 2*y(n) + 1, y(n), list(y(0) = 0))");
        assert_eq!(geometric, "y(n) = 2^n - 1");
        // Classification.
        // Fourier method on the whole line: the solutions decaying at ±∞,
        // checked by central differences away from the kink at 0.
        for (src, order, forcing) in [
            ("ode_fourier(diff(diff(y(x), x), x) - y(x) = -2*exp(-abs(x)), y(x))", 2, "-2*exp(-abs(x))"),
            ("ode_fourier(diff(y(x), x) + y(x) = exp(-abs(x)), y(x))", 1, "exp(-abs(x))"),
        ] {
            let got = simplify(&rules, src);
            assert!(got.starts_with("y(x) = ") && !got.contains("fourier"), "{got}");
            let rhs = got.trim_start_matches("y(x) = ");
            let y = |at: f64| crate::rules::testing::eval(&rules, rhs, &[("x", at)]);
            for at in [0.7_f64, 1.9, -1.3, -0.4] {
                let h = 1e-4;
                let lhs = if order == 2 { (y(at + h) - 2.0 * y(at) + y(at - h)) / (h * h) - y(at) } else { (y(at + h) - y(at - h)) / (2.0 * h) + y(at) };
                let want = crate::rules::testing::eval(&rules, forcing, &[("x", at)]);
                assert!((lhs - want).abs() < 1e-5, "{got}: {lhs} vs {want} at {at}");
            }
        }
        let classes = simplify(&rules, "ode_classify(diff(y(x), x) = x*y(x), y(x))");
        assert!(classes.contains("linear") && classes.contains("separable"), "{classes}");
    }

    #[test]
    fn constant_coefficients() {
        assert_eq!(
            run("dsolve(diff(diff(y(x), x), x) - 3*diff(y(x), x) + 2*y(x) = 0, y(x))"),
            "y(x) = C1*exp(x) + C2*exp(2*x)"
        );
        assert_eq!(run("dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))"), "y(x) = C1*cos(x) + C2*sin(x)");
        check("diff(diff(y(x), x), x) + 2*diff(y(x), x) + y(x) = 0", 2);
        check("diff(diff(y(x), x), x) + 2*diff(y(x), x) + 5*y(x) = 0", 2);
        check("diff(diff(y(x), x), x) - y(x) = x", 2);
        check("diff(diff(y(x), x), x) + y(x) = exp(x)", 2);
        check("diff(diff(y(x), x), x) + 4*y(x) = sin(x)", 2);
        check("diff(diff(diff(y(x), x), x), x) - 6*diff(diff(y(x), x), x) + 11*diff(y(x), x) - 6*y(x) = 0", 3);
    }

    #[test]
    fn cauchy_euler() {
        check("x^2*diff(diff(y(x), x), x) - 2*y(x) = 0", 2);
        check("x^2*diff(diff(y(x), x), x) + x*diff(y(x), x) - y(x) = 0", 2);
        check("x^2*diff(diff(y(x), x), x) - x*diff(y(x), x) + y(x) = 0", 2);
        check("x^2*diff(diff(y(x), x), x) + x*diff(y(x), x) + 4*y(x) = 0", 2);
    }

    #[test]
    fn initial_conditions() {
        assert_eq!(run("dsolve(diff(y(x), x) = y(x), y(x), list(y(0) = 3))"), "y(x) = 3*exp(x)");
        assert_eq!(
            run("dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x), list(y(0) = 1, diff(y(x), x) = 2))"),
            "y(x) = cos(x) + 2*sin(x)"
        );
        assert_eq!(run("dsolve(diff(y(x), x) = 2*x, y(x), list(y(1) = 5))"), "y(x) = x^2 + 4");
    }

    #[test]
    fn unsolved_equations_stay_requests() {
        for src in [
            "dsolve(diff(y(x), x) = sin(x*y(x)), y(x))",
            "dsolve(diff(diff(y(x), x), x) + sin(y(x)) = x, y(x))",
            "dsolve(y(x) = x, y(x))",
        ] {
            let (text, reduced) = reduce_with(&[ode()], src, &[]);
            assert!(!reduced, "{src} unexpectedly gave {text}");
        }
    }

    #[test]
    fn more_first_order_classes() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        // Clairaut: y = C x + g(C).
        assert_eq!(run("dsolve(y(x) = x*diff(y(x), x) + diff(y(x), x)^2, y(x))"), "y(x) = C1^2 + C1*x");
        assert_eq!(run("dsolve(y(x) = x*diff(y(x), x) - diff(y(x), x)^3, y(x))"), "y(x) = C1*x - C1^3");
        // d'Alembert and x-solved equations: parametric answers in p = y'.
        let lagrange = run("dsolve(y(x) = 2*x*diff(y(x), x) + diff(y(x), x)^2, y(x))");
        assert!(lagrange.starts_with("list(x = ") && lagrange.contains("p^3") && lagrange.contains("C1"), "{lagrange}");
        let solved_for_x = run("dsolve(x = diff(y(x), x)^3 + diff(y(x), x), y(x))");
        assert_eq!(solved_for_x, "list(x = p^3 + p, y(x) = 3/4*p^4 + 1/2*p^2 + C1)");
        // Integrating factors mu(x), mu(y) and e^x.
        check("(x^2 + y(x)^2 + x) + x*y(x)*diff(y(x), x) = 0", 1);
        check("y(x)^2 + (3*x*y(x) - 1)*diff(y(x), x) = 0", 1);
        check("y(x)*(x + y(x)) + (x + 2*y(x) - 1)*diff(y(x), x) = 0", 1);
        // Riccati: a particular solution exp(x), and through the linear
        // second-order equation (Airy functions) with one constant.
        assert_eq!(
            run("dsolve(diff(y(x), x) = y(x)^2 - 2*y(x)*exp(x) + exp(2*x) + exp(x), y(x))"),
            "y(x) = exp(x) + 1/(C1 - x)"
        );
        let airy = run("dsolve(diff(y(x), x) = y(x)^2 + x, y(x))");
        assert!(airy.contains("airyai(-x)") && airy.contains("C1") && !airy.contains("C2"), "{airy}");
        // Bernoulli with a square root.
        let bernoulli = run("dsolve(diff(y(x), x) + y(x) = x*sqrt(y(x)), y(x))");
        assert_eq!(bernoulli, "y(x) = (x*exp(1/2*x) + C1 - 2*exp(1/2*x))^2*exp(-x)");
        // First-order equations quadratic in y' that factor.
        assert_eq!(run("dsolve(diff(y(x), x)^2 - (x + y(x))*diff(y(x), x) + x*y(x) = 0, y(x))"), "y(x) = 1/2*x^2 + C1");
    }

    #[test]
    fn second_order_normal_form_and_reduction() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        // y = u/x removes the first-derivative term: spherical Bessel.
        assert_eq!(
            run("dsolve(x*diff(diff(y(x), x), x) + 2*diff(y(x), x) + x*y(x) = 0, y(x))"),
            "y(x) = C1*cos(x)/x + C2*sin(x)/x"
        );
        // A polynomial solution and a quadrature for the second (Hermite).
        let hermite = run("dsolve(diff(diff(y(x), x), x) - x*diff(y(x), x) + y(x) = 0, y(x))");
        assert!(hermite.contains("C1*x") && hermite.contains("integral(exp(1/2*x^2)/x^2, x)"), "{hermite}");
    }

    #[test]
    fn nonlinear_and_reduced_systems() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        // Triangular: the first equation is solved alone.
        let triangular = run("dsolve(list(diff(x(t), t) = -x(t)^2, diff(y(t), t) = x(t)*y(t)), list(x(t), y(t)))");
        assert!(triangular.contains("x(t) = 1/(C1 + t)") && triangular.contains("C2"), "{triangular}");
        // Higher order by reduction to first order.
        let reduced = run("dsolve(list(diff(diff(x(t), t), t) = -x(t) + y(t), diff(diff(y(t), t), t) = x(t) - y(t)), list(x(t), y(t)))");
        assert!(reduced.contains("C4") && reduced.contains("cos(2^(1/2)*t)"), "{reduced}");
        // Variable coefficients a(t) M: constant system in s = integral of a.
        let scaled = run("dsolve(list(diff(x(t), t) = t*y(t), diff(y(t), t) = t*x(t)), list(x(t), y(t)))");
        assert!(scaled.contains("exp(1/2*t^2)") && scaled.contains("exp(-1/2*t^2)"), "{scaled}");
        // Two autonomous equations: the first integral (energy).
        let energy = run("dsolve(list(diff(x(t), t) = y(t), diff(y(t), t) = -sin(x(t))), list(x(t), y(t)))");
        assert!(energy.contains("cos(x(t))") && energy.contains("C1"), "{energy}");
        let predator_prey = run("dsolve(list(diff(x(t), t) = x(t)*(2 - y(t)), diff(y(t), t) = y(t)*(x(t) - 3)), list(x(t), y(t)))");
        assert!(predator_prey.contains("ln(x(t))") && predator_prey.contains("ln("), "{predator_prey}");
        // Polynomial first integrals of a three-dimensional system.
        let rigid = run("dsolve(list(diff(x(t), t) = y(t)*z(t), diff(y(t), t) = -x(t)*z(t), diff(z(t), t) = -x(t)*y(t)/2), list(x(t), y(t), z(t)))");
        assert!(rigid.contains("x(t)^2 + y(t)^2 - C1 = 0") && rigid.contains("C2"), "{rigid}");
    }

    #[test]
    fn numeric_integration_agrees_with_closed_forms() {
        // y' = -2y, y(0) = 1 at x = 1.
        let (value, error) = numeric(&[ode()], "odeint(-2*y, y, x, 0, 1, 1)", &[], 1e-10);
        assert!((value - (-2.0_f64).exp()).abs() < 1e-8, "{value} ± {error}");
        // With a parameter from the bindings.
        let (value, _) = numeric(&[ode()], "odeint(a*y, y, x, 0, 1, 2)", &[("a", 0.5)], 1e-10);
        assert!((value - 1.0_f64.exp()).abs() < 1e-7, "{value}");
        // A system: harmonic oscillator over a quarter period.
        let mut g = Graph::new();
        let engine = crate::graph::Engine::install(&mut g, &[ode()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse("odeint(list(v, -u), list(u, v), t, 0, list(1, 0), pi/2)").unwrap_or_else(|e| panic!("{e}"));
        engine.run(&mut g, &[root], &Env::numeric(1e-10), &crate::graph::Saturate, &crate::graph::Budget::default());
        let out = crate::graph::Extractor::new(&g, &[root], &crate::graph::ClosedForm).build(&mut g, root);
        let text = out.map(|n| g.display(n)).unwrap_or_default();
        let values: Vec<f64> = text
            .trim_start_matches("list(")
            .trim_end_matches(')')
            .split(", ")
            .filter_map(|v| v.parse().ok())
            .collect();
        assert_eq!(values.len(), 2, "{text}");
        assert!(values[0].abs() < 1e-7 && (values[1] + 1.0).abs() < 1e-7, "{text}");
    }
}
