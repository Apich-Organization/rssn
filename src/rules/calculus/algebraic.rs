//! Integration of algebraic integrands: Hermite reduction in a quadratic
//! extension `Q(x)(sqrt(D))`.
//!
//! An integrand `R(x, sqrt(D))` is written as `(A + B sqrt(D)) / C` with
//! polynomials `A`, `B`, `C`.  The radicand is first made square-free
//! (`D = d s^2`, `sqrt(D) = s sqrt(d)`), and then:
//!
//! * a constant radicand is a rational function once the square root is a
//!   rational number;
//! * a linear radicand is rationalised by `t = sqrt(D)`;
//! * a quadratic radicand is parametrised as a conic (see
//!   [`quadratic`]);
//! * a radicand of degree three or more is reduced by [`hermite`]: every
//!   repeated factor of `C` is lowered to multiplicity one by subtracting
//!   the derivative of an element of `Q(x)(sqrt(D))`; the square-free part
//!   is then integrated by residues and logarithms in [`logpart`] as far as
//!   those are expressible over `Q(x)`.
//!
//! # The logarithmic part and its limit
//!
//! After Hermite reduction the remaining denominator is square-free, so the
//! only residues are at its roots `alpha`, where the residue is an element of
//! the residue field `Q(alpha, y_alpha)` of the curve.  It is computed
//! exactly in the algebra `Q[alpha, y]/(w(alpha), y^2 - D(alpha))` as
//! `a + b y`.  The logarithmic part is then found in two complementary ways:
//!
//! * a residue without a `y` component (`b = 0`) is summed over the
//!   conjugates of `alpha`, which gives an element of `Q(x)`: linear factors
//!   give logarithms, an irreducible quadratic gives a logarithm/arctangent
//!   pair, and a constant residue of a factor of degree three or more gives
//!   a logarithm;
//! * a residue with a `y` component is the generic positive-genus case.  Its
//!   `y`-trace `2a` is the base-field projection (the sum of the two
//!   conjugate residues at the same `x = alpha`), and it is *not* an
//!   antiderivative by itself: the anti-trace component is not a logarithm
//!   of a rational function.  When the primitive is nevertheless elementary
//!   it is a logarithm `log(g)` of an algebraic function `g` (for instance
//!   `integral(3 D' (y - D^2) / (2 (1 - D^3)), x) = log(1 + D y)`,
//!   `y^2 = D`).  [`logpart`] solves the linear system `g' = f g` for
//!   `g = P + Q y` exactly and checks the trace of the answer
//!   (`N'/N = 2a/c` for the norm `N = P^2 - Q^2 D`).
//!
//! A residue with a `y` component that is not such a logarithmic derivative
//! is genuinely non-elementary: for a radicand of degree three or four the
//! curve has positive genus and the remaining differential of the second or
//! third kind has non-zero periods.  In that case the kernel returns the
//! Hermite-reduced integral `rational_part + integral((A0 + B0 y)/C0, x)`,
//! which is a correct identity with the lowest possible denominator
//! multiplicities, rather than a wrong closed form.  Repeated factors that
//! divide `D` (ramified places) and non-square constant radicands likewise
//! stay unevaluated.
//!
//! Every candidate answer is differentiated and compared with the
//! integrand at sample points before it is returned; a mismatch yields
//! [`Outcome::Pass`](crate::graph::Outcome::Pass) instead.

use std::cell::Cell;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Zero;

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

use super::algebraic_field as field;
use super::algebraic_hermite as hermite;
use super::algebraic_logpart as logpart;
use super::algebraic_quadratic as quadratic;
use super::algebraic_term as term;

use field::Element;
use field::Frac;

/// Registers the `integral` kernel and teaches the ordinary integrator the
/// algebraic entry point.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let integral = i.graph().ops().lookup("integral").ok_or_else(|| {
        RuleError::Invalid {
            rule: "calculus/algebraic".to_owned(),
            reason: "the `integral` operator is missing",
        }
    })?;
    i.kernel("calculus/algebraic", Tier::Reduce, Algebraic { integral });
    // The staged integrator consults the table for integrands its own
    // methods cannot decide.  Only radicands of degree three or more are
    // offered: for a quadratic radicand the staged integrator's own conic
    // rules return shorter forms, and pre-empting them here would only
    // replace them by the (equivalent but longer) parametrised ones.
    super::integral_rule(i, integrate_algebraic_table);
    Ok(())
}

/// The table entry: only radicands of degree three or more, and only a
/// *fully closed* answer.
///
/// A Hermite-reduced result that still contains an `integral` must not be
/// handed to the staged integrator: `defint` substitutes the limits into
/// the primitive it obtains, and the `integral` operator does not bind its
/// variable, so the substitution would corrupt the inner request.
fn integrate_algebraic_table(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let radicand = find_radicand(cx.graph, f, x)?;
    let d = polynomial_of(cx.graph, radicand, x)?;
    let (d, _) = normalise(&d);
    if field::degree(&d)? < 3 {
        return None;
    }
    let answer = integrate_algebraic(cx, f, x)?;
    if contains_integral(cx.graph, answer) {
        return None;
    }
    Some(answer)
}

/// Whether `term` contains an `integral` request.
fn contains_integral(
    graph: &Graph,
    term: NodeId,
) -> bool {
    let Some(integral) = graph.ops().lookup("integral") else {
        return false;
    };
    let mut stack = vec![term];
    let mut seen = Vec::new();
    while let Some(node) = stack.pop() {
        if seen.contains(&node) {
            continue;
        }
        seen.push(node);
        if graph.op(node) == integral {
            return true;
        }
        stack.extend_from_slice(graph.children(node));
    }
    false
}

/// Kernel on `integral(f, x)`.
struct Algebraic {
    integral: crate::graph::OpId,
}

impl Kernel for Algebraic {
    fn ops(&self) -> Vec<crate::graph::OpId> {
        vec![self.integral]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[integrand, variable] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        match integrate_algebraic(cx, integrand, variable) {
            | Some(answer) if answer != node => Outcome::Equal(answer),
            | _ => Outcome::Pass,
        }
    }
}

thread_local! {
    /// Non-zero while [`integrate_algebraic`] runs on this thread; the
    /// entry point also sits in the integrator's table, and must not call
    /// itself through it.
    static RUNNING: Cell<u32> = const { Cell::new(0) };
}

/// Resets [`RUNNING`] when the call ends.
struct Running;

impl Drop for Running {
    fn drop(&mut self) {
        RUNNING.with(|depth| depth.set(depth.get().saturating_sub(1)));
    }
}

/// A verified antiderivative of `f` with respect to `x` when `f` is an
/// algebraic function `R(x, sqrt(D))`, or `None`.
///
/// The result may still contain an `integral` when Hermite reduction
/// succeeded but the remaining integral is not elementary (see the module
/// documentation); it is always an identity, verified by differentiation.
#[allow(clippy::redundant_pub_crate)] // the entry point is crate-visible by design
pub(crate) fn integrate_algebraic(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    if RUNNING.with(Cell::get) > 0 {
        return None;
    }
    RUNNING.with(|depth| depth.set(depth.get().saturating_add(1)));
    let _running = Running;
    attempt(cx, f, x)
}

/// The body of [`integrate_algebraic`], with the recursion guard active.
fn attempt(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let radicand_node = find_radicand(cx.graph, f, x)?;
    let d = polynomial_of(cx.graph, radicand_node, x)?;
    let element = parse_term(cx.graph, f, x, radicand_node, &d)?;
    let (d, s) = normalise(&d);
    let one = BigRational::from_integer(BigInt::from(1));
    let element = if s.len() == 1 && s.first() == Some(&one) {
        element
    } else {
        let s_frac = Frac::polynomial(&s);
        Element {
            a: element.a,
            b: element.b.mul(&s_frac)?,
        }
    };
    let degree = field::degree(&d)?;
    if degree == 0 {
        return constant_radicand(cx, f, x, &d, &element);
    }
    let (a, b, c) = element.to_quotient()?;
    match degree {
        | 1 => low_degree(cx, f, x, &d, &a, &b, &c, false),
        | 2 => low_degree(cx, f, x, &d, &a, &b, &c, true),
        | _ => {
            let radicand = term::poly_node(cx.graph, x, &d);
            let sqrt_term = term::sqrt_node(cx.graph, radicand);
            let candidate = higher_degree(cx, f, x, &d, sqrt_term, &a, &b, &c)?;
            accept(cx, f, x, candidate)
        },
    }
}

/// The closed forms of the exact derivatives in `pieces`.
fn rational_part_term(
    cx: &mut Cx<'_>,
    x: NodeId,
    sqrt_term: NodeId,
    pieces: &[hermite::Piece],
) -> NodeId {
    let mut terms = Vec::new();
    for piece in pieces {
        let power = field::pow(&piece.w, piece.power);
        let u_node = term::poly_node(cx.graph, x, &piece.u);
        let v_node = term::poly_node(cx.graph, x, &piece.v);
        let vy = cx.graph.node(core::MUL, &[v_node, sqrt_term]);
        let body = cx.graph.node(core::ADD, &[u_node, vy]);
        let den = term::poly_node(cx.graph, x, &power);
        terms.push(term::div_node(cx.graph, body, den));
    }
    term::add_all(cx.graph, &terms)
}

/// Integration for a radicand of degree one or two.
///
/// Hermite reduction is tried first (it removes repeated factors of the
/// denominator that are co-prime to `D`, which keeps the rational function
/// handed to the conic parametrisation as small as possible), then the
/// conic parametrisation itself.
#[allow(clippy::too_many_arguments)] // the integrand, the radicand and its quotient form
fn low_degree(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    quadratic: bool,
) -> Option<NodeId> {
    let reduced = hermite::reduce(d, a, b, c);
    if reduced.steps > 0 {
        let radicand = term::poly_node(cx.graph, x, d);
        let sqrt_term = term::sqrt_node(cx.graph, radicand);
        let rational = rational_part_term(cx, x, sqrt_term, &reduced.pieces);
        let rest = if quadratic {
            quadratic::integrate_quadratic(cx, x, &reduced.a, &reduced.b, &reduced.c, d)
        } else {
            quadratic::integrate_linear(cx, x, &reduced.a, &reduced.b, &reduced.c, d)
        };
        if let Some(rest) = rest {
            let candidate = cx.graph.node(core::ADD, &[rational, rest]);
            if let Some(verified) = accept(cx, f, x, candidate) {
                return Some(verified);
            }
        }
    }
    let direct = if quadratic {
        quadratic::integrate_quadratic(cx, x, a, b, c, d)?
    } else {
        quadratic::integrate_linear(cx, x, a, b, c, d)?
    };
    accept(cx, f, x, direct)
}

/// The second half of a radicand of degree three or more: Hermite
/// reduction, then the logarithmic part when it closes.
#[allow(clippy::too_many_arguments)] // the integrand, the radicand and its quotient form
fn higher_degree(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    d: &[BigRational],
    sqrt_term: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> Option<NodeId> {
    // Try the logarithmic part on the original integrand first.  A primitive
    // that is a logarithm of an algebraic function (the only way a
    // positive-genus residue with a `y` component can close) is destroyed by
    // subtracting Hermite pieces from it, because those pieces are
    // algebraic and not logarithmic derivatives.  Solving the logarithmic
    // ansatz on the unreduced integrand keeps that primitive available.
    if let Some(logs) = logpart::complete(cx.graph, x, d, a, b, c)
        && accept(cx, f, x, logs).is_some()
    {
        return Some(logs);
    }
    let reduced = hermite::reduce(d, a, b, c);
    let rational = rational_part_term(cx, x, sqrt_term, &reduced.pieces);
    if let Some(logs) = logpart::complete(cx.graph, x, d, &reduced.a, &reduced.b, &reduced.c) {
        let candidate = cx.graph.node(core::ADD, &[rational, logs]);
        if accept(cx, f, x, candidate).is_some() {
            return Some(candidate);
        }
    }
    if reduced.steps == 0 {
        return None;
    }
    // The reduction made progress but did not close: return
    // rational_part + integral(remaining), which is still an identity.
    let a_node = term::poly_node(cx.graph, x, &reduced.a);
    let b_node = term::poly_node(cx.graph, x, &reduced.b);
    let by = cx.graph.node(core::MUL, &[b_node, sqrt_term]);
    let body = cx.graph.node(core::ADD, &[a_node, by]);
    let den = term::poly_node(cx.graph, x, &reduced.c);
    let fraction = term::div_node(cx.graph, body, den);
    let inner = term::integral_node(cx.graph, fraction, x)?;
    let candidate = cx.graph.node(core::ADD, &[rational, inner]);
    if accept(cx, f, x, candidate).is_some() {
        return Some(candidate);
    }
    None
}

/// A radicand that became a constant: `sqrt(D) = s sqrt(c)`.
fn constant_radicand(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    d: &[BigRational],
    element: &Element,
) -> Option<NodeId> {
    let c = d.first()?.clone();
    let sqrt_c = field::rational_sqrt(&c)?;
    let rational = element.a.add(&element.b.scale_poly(&[sqrt_c])?)?;
    let body = term::frac_node(cx.graph, x, &rational);
    let candidate = super::antiderivative(cx, body, x)?;
    accept(cx, f, x, candidate)
}

/// The radicand node of the square root(s) in `f`, or `None` when `f`
/// contains no square root of a polynomial in `x`, or square roots of
/// different polynomials.
fn find_radicand(
    graph: &Graph,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let sqrt = graph.ops().lookup("sqrt");
    let symbol = graph.symbol_of(x)?;
    let mut found: Option<NodeId> = None;
    let mut stack = vec![f];
    let mut seen = Vec::new();
    while let Some(node) = stack.pop() {
        if seen.contains(&node) {
            continue;
        }
        seen.push(node);
        let children = graph.children(node);
        let op = graph.op(node);
        let radicand = if Some(op) == sqrt && children.len() == 1 {
            children.first().copied()
        } else if op == core::POW && children.len() == 2 {
            let exponent = children
                .get(1)
                .and_then(|&e| graph.number_of(e))
                .and_then(Number::to_rational);
            match exponent {
                | Some(exponent)
                    if exponent.denom() == &BigInt::from(2) && !exponent.is_integer() =>
                {
                    children.first().copied()
                },
                | _ => None,
            }
        } else {
            None
        };
        if let Some(radicand) = radicand
            && graph.depends_on(graph.find(radicand), symbol)
        {
            match found {
                | None => found = Some(radicand),
                | Some(previous) if graph.same(previous, radicand) => {},
                | Some(_) => return None,
            }
        }
        stack.extend_from_slice(children);
    }
    found
}

/// The radicand as a polynomial over `Q`.
fn polynomial_of(
    graph: &mut Graph,
    node: NodeId,
    x: NodeId,
) -> Option<QPoly> {
    let (numer, denom) = crate::rules::poly::rational_function_in(graph, node, x)?;
    match denom.as_slice() {
        | [c] if !c.is_zero() => Some(numer.iter().map(|a| a / c).collect()),
        | _ => None,
    }
}

/// Reads the concrete term `node` as an element of `Q(x)(sqrt(D))`.
fn parse_term(
    graph: &mut Graph,
    node: NodeId,
    x: NodeId,
    radicand: NodeId,
    d: &[BigRational],
) -> Option<Element> {
    let sqrt = graph.ops().lookup("sqrt");
    if node == x {
        return Some(Element {
            a: Frac::new(
                &[
                    BigRational::zero(),
                    BigRational::from_integer(BigInt::from(1)),
                ],
                &[BigRational::from_integer(BigInt::from(1))],
            )?,
            b: Frac::zero(),
        });
    }
    if let Some(number) = graph.number_of(node) {
        return Some(Element::constant(number.to_rational()?));
    }
    let children = graph.children(node).to_vec();
    let op = graph.op(node);
    if Some(op) == sqrt && children.len() == 1 {
        let arg = *children.first()?;
        if graph.same(arg, radicand) {
            return Some(Element {
                a: Frac::zero(),
                b: Frac::polynomial(&[BigRational::from_integer(BigInt::from(1))]),
            });
        }
        return None;
    }
    if op == core::ADD {
        let mut out = Element::zero();
        for child in children {
            out = out.add(&parse_term(graph, child, x, radicand, d)?)?;
        }
        return Some(out);
    }
    if op == core::MUL {
        let mut out = Element::constant(BigRational::from_integer(BigInt::from(1)));
        for child in children {
            out = out.mul(&parse_term(graph, child, x, radicand, d)?, d)?;
        }
        return Some(out);
    }
    if op == core::POW && children.len() == 2 {
        let base = *children.first()?;
        let exponent = graph
            .number_of(*children.get(1)?)
            .and_then(Number::to_rational)?;
        let value = parse_term(graph, base, x, radicand, d)?;
        if exponent.is_integer() {
            let n = exponent.to_integer().to_string().parse::<i64>().ok()?;
            if n.unsigned_abs() > MAX_EXTENSION_POWER {
                return None;
            }
            return value.pow(n, d);
        }
        // `D^((2k+1)/2) = y^(2k+1)`.
        if graph.same(base, radicand) {
            let doubled = &exponent * BigRational::from_integer(BigInt::from(2));
            if !doubled.is_integer() {
                return None;
            }
            let k = doubled.to_integer().to_string().parse::<i64>().ok()?;
            if k.unsigned_abs() > MAX_EXTENSION_POWER {
                return None;
            }
            let root = Element {
                a: Frac::zero(),
                b: Frac::polynomial(&[BigRational::from_integer(BigInt::from(1))]),
            };
            return root.pow(k, d);
        }
        return None;
    }
    None
}

/// Largest power of `y` the parser expands explicitly.
///
/// `sqrt(D)^n` is `y^n` in the field, and its denominator/numerator degrees
/// grow linearly in `n`; a user-supplied exponent such as
/// `sqrt(x^3 + 1)^1000000` would otherwise allocate polynomials of degree
/// millions before any of the reduction caps applies.  The largest power
/// that actually occurs in an integrand of the supported shape is `y^3`
/// (from `D^(3/2)`), so the bound is generous.
const MAX_EXTENSION_POWER: u64 = 64;

/// Splits `D = unit * radical * s^2` with `radical` square-free, returning
/// `unit * radical` and `s`.
fn normalise(d: &[BigRational]) -> (QPoly, QPoly) {
    let one = BigRational::from_integer(BigInt::from(1));
    let mut radical: QPoly = vec![one.clone()];
    let mut s: QPoly = vec![one];
    let (unit, factors) = univariate::factor(d);
    for (f, e) in &factors {
        let fq: QPoly = f
            .iter()
            .map(|c| BigRational::from_integer(c.clone()))
            .collect();
        for _ in 0..(e / 2) {
            s = field::mul(&s, &fq);
        }
        if e % 2 == 1 {
            radical = field::mul(&radical, &fq);
        }
    }
    let new_d = field::scale(&radical, &unit);
    (new_d, s)
}

/// Accepts `candidate` when its derivative agrees with `f`, returning the
/// candidate (or a variant with `ln` made real by `abs`).
fn accept(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    candidate: NodeId,
) -> Option<NodeId> {
    // Prefer the real variant: it evaluates on both sides of a zero of a
    // logarithm's argument, where the complex logarithm does not.
    let wrapped = real_logs(cx.graph, candidate);
    if wrapped != candidate && verifies(cx, f, x, wrapped) {
        return Some(wrapped);
    }
    if verifies(cx, f, x, candidate) {
        return Some(candidate);
    }
    None
}

/// Replaces every `ln(u)` by `ln(abs(u))`, the real antiderivative on each
/// interval where `u` keeps its sign.
fn real_logs(
    graph: &mut Graph,
    term: NodeId,
) -> NodeId {
    let Some(ln) = graph.ops().lookup("ln") else {
        return term;
    };
    let Some(abs) = graph.ops().lookup("abs") else {
        return term;
    };
    let mut stack = vec![(term, false)];
    let mut done: std::collections::HashMap<NodeId, NodeId> = std::collections::HashMap::new();
    while let Some((node, expanded)) = stack.pop() {
        if done.contains_key(&node) {
            continue;
        }
        let children = graph.children(node).to_vec();
        if children.is_empty() {
            done.insert(node, node);
            continue;
        }
        if graph.op(node) == ln {
            let Some(&arg) = children.first() else {
                done.insert(node, node);
                continue;
            };
            if graph.op(arg) == abs {
                done.insert(node, node);
                continue;
            }
            let replacement = graph.node(abs, &[arg]);
            let rebuilt = graph.node(ln, &[replacement]);
            done.insert(node, rebuilt);
            continue;
        }
        if !expanded {
            stack.push((node, true));
            stack.extend(
                children
                    .iter()
                    .filter(|c| !done.contains_key(c))
                    .map(|&c| (c, false)),
            );
            continue;
        }
        let rebuilt: Vec<NodeId> = children
            .iter()
            .map(|c| done.get(c).copied().unwrap_or(*c))
            .collect();
        let new = if rebuilt == children {
            node
        } else {
            graph.try_node(graph.op(node), &rebuilt).unwrap_or(node)
        };
        done.insert(node, new);
    }
    done.get(&term).copied().unwrap_or(term)
}

/// The derivative of `term` with respect to `x`, reading `integral(u, x)`
/// as a primitive of `u`.
fn derivative_of(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let derivative = crate::rules::calculus::derivative(graph, term, x)?;
    let diff = graph.ops().lookup("diff")?;
    let integral = graph.ops().lookup("integral")?;
    let mut stack = vec![(derivative, false)];
    let mut done: std::collections::HashMap<NodeId, NodeId> = std::collections::HashMap::new();
    while let Some((node, expanded)) = stack.pop() {
        if done.contains_key(&node) {
            continue;
        }
        let children = graph.children(node).to_vec();
        if children.is_empty() {
            done.insert(node, node);
            continue;
        }
        if graph.op(node) == diff && children.len() == 2 {
            let inner = *children.first()?;
            let variable = *children.get(1)?;
            if variable == x && graph.op(inner) == integral {
                let arguments = graph.children(inner);
                if arguments.len() == 2 && arguments.get(1) == Some(&x) {
                    let body = *arguments.first()?;
                    done.insert(node, body);
                    continue;
                }
            }
        }
        if !expanded {
            stack.push((node, true));
            stack.extend(
                children
                    .iter()
                    .filter(|c| !done.contains_key(c))
                    .map(|&c| (c, false)),
            );
            continue;
        }
        let rebuilt: Vec<NodeId> = children
            .iter()
            .map(|c| done.get(c).copied().unwrap_or(*c))
            .collect();
        let new = if rebuilt == children {
            node
        } else {
            graph.try_node(graph.op(node), &rebuilt).unwrap_or(node)
        };
        done.insert(node, new);
    }
    done.get(&derivative).copied()
}

/// Sample abscissae used to compare a candidate derivative with the
/// integrand; they avoid the usual branch points (`0`, `±1`, `±i`).
const SAMPLES: [f64; 18] = [
    -0.87, -0.63, -0.41, -0.19, 0.11, 0.29, 0.53, 0.77, 0.94, 1.23, 1.87, 2.63, 3.71, 4.51, -1.47,
    -2.39, -3.17, -4.43,
];

/// Numerically checks that `candidate' = f` at sample points.
fn verifies(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    candidate: NodeId,
) -> bool {
    let Some(derivative) = derivative_of(cx.graph, candidate, x) else {
        return false;
    };
    let Some(target) = cx.graph.symbol_of(x) else {
        return false;
    };
    let mut symbols: Vec<crate::graph::SymbolId> = cx.graph.free_symbols(cx.graph.find(f)).to_vec();
    symbols.extend_from_slice(cx.graph.free_symbols(cx.graph.find(derivative)));
    symbols.sort_unstable();
    symbols.dedup();
    let mut checked = 0_u32;
    for point in SAMPLES {
        let mut env = crate::graph::Env::numeric(0.0);
        for (k, &symbol) in symbols.iter().enumerate() {
            let value = if symbol == target {
                point
            } else {
                #[allow(clippy::cast_precision_loss)]
                {
                    0.6 + 0.17 * k as f64
                }
            };
            env.bind(symbol, value);
        }
        let (Some(want), Some(got)) = (cx.graph.eval(f, &env), cx.graph.eval(derivative, &env))
        else {
            continue;
        };
        if !want.is_finite() || !got.is_finite() {
            continue;
        }
        checked = checked.saturating_add(1);
        let scale = 1.0 + want.abs().max(got.abs());
        if (want - got).abs() > 1e-5 * scale {
            return false;
        }
    }
    checked >= 3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::rules::calculus;
    use crate::rules::testing::reduce_with;

    fn one() -> BigRational {
        BigRational::from_integer(BigInt::from(1))
    }

    /// `d/dx (a + b y)` as an element, for the exactness check below.
    fn differentiate(
        element: &Element,
        d: &[BigRational],
    ) -> Element {
        let frac_derivative = |f: &Frac| -> Frac {
            let numerator = field::sub(
                &field::mul(&field::derivative(&f.num), &f.den),
                &field::mul(&f.num, &field::derivative(&f.den)),
            );
            let denominator = field::mul(&f.den, &f.den);
            let Some(value) = Frac::new(&numerator, &denominator) else {
                panic!("zero denominator");
            };
            value
        };
        let a = frac_derivative(&element.a);
        let bprime = frac_derivative(&element.b);
        let two_d = field::scale(d, &BigRational::from_integer(BigInt::from(2)));
        let dp = Frac::polynomial(&field::derivative(d));
        let Some(tail) = element
            .b
            .mul(&dp)
            .and_then(|x| x.div(&Frac::polynomial(&two_d)))
        else {
            panic!("division by zero");
        };
        let Some(b) = bprime.add(&tail) else {
            panic!("zero denominator");
        };
        Element { a, b }
    }

    #[test]
    fn square_factors_are_split_off() {
        // D = 4 (x^2 + 1)^2: the square part is x^2 + 1 and the remaining
        // radicand is the constant 4.
        let w = vec![one(), BigRational::zero(), one()];
        let square = field::mul(&w, &w);
        let d = field::scale(&square, &BigRational::from_integer(BigInt::from(4)));
        let (new_d, s) = normalise(&d);
        assert_eq!(new_d, vec![BigRational::from_integer(BigInt::from(4))]);
        assert_eq!(s, w);
    }

    #[test]
    fn hermite_step_is_exact() {
        // D = x^4 + 1, integrand sqrt(D)/(x^2+1)^2.
        let d = vec![
            one(),
            BigRational::zero(),
            BigRational::zero(),
            BigRational::zero(),
            one(),
        ];
        let w = vec![one(), BigRational::zero(), one()];
        let c = field::mul(&w, &w);
        let b = vec![one()];
        let reduced = hermite::reduce(&d, &[], &b, &c);
        assert!(reduced.steps >= 1, "no Hermite step taken");
        // The multiplicity of `w` must have dropped.
        let multiplicity = |p: &[BigRational], w: &[BigRational]| -> u32 {
            let mut m = 0_u32;
            let mut current = p.to_vec();
            while let Some((quotient, remainder)) = field::divrem(&current, w) {
                if !remainder.is_empty() {
                    break;
                }
                current = quotient;
                m = m.saturating_add(1);
            }
            m
        };
        let before = multiplicity(&c, &w);
        let after = multiplicity(&reduced.c, &w);
        assert!(before >= 2, "w should divide C twice");
        assert!(after < before, "the multiplicity of w did not drop");
        let Some(original) = Element::quotient(&[], &b, &c) else {
            panic!("bad quotient");
        };
        let mut total = Element::zero();
        for piece in &reduced.pieces {
            let den = field::pow(&piece.w, piece.power);
            let Some(term) = Element::quotient(&piece.u, &piece.v, &den) else {
                panic!("bad piece");
            };
            let Some(sum) = total.add(&differentiate(&term, &d)) else {
                panic!("add");
            };
            total = sum;
        }
        let Some(remaining) = Element::quotient(&reduced.a, &reduced.b, &reduced.c) else {
            panic!("bad remainder");
        };
        let Some(sum) = total.add(&remaining) else {
            panic!("add");
        };
        let Some(difference) = original.sub(&sum) else {
            panic!("sub");
        };
        let Some((a, b, _)) = difference.to_quotient() else {
            panic!("quotient");
        };
        assert!(a.is_empty() && b.is_empty(), "the reduction is not exact");
    }

    /// `d/dx (a + b y) = g` exactly, with `y^2 = D`.
    fn assert_reduction_is_exact(
        d: &[BigRational],
        a: &[BigRational],
        b: &[BigRational],
        c: &[BigRational],
    ) {
        let Some(original) = Element::quotient(a, b, c) else {
            panic!("bad quotient");
        };
        let reduced = hermite::reduce(d, a, b, c);
        assert!(
            reduced.steps > 0,
            "the Hermite reduction took no step for D = {d:?}, C = {c:?}"
        );
        let mut total = Element::zero();
        for piece in &reduced.pieces {
            let den = field::pow(&piece.w, piece.power);
            let Some(term) = Element::quotient(&piece.u, &piece.v, &den) else {
                panic!("bad piece");
            };
            let Some(sum) = total.add(&differentiate(&term, d)) else {
                panic!("add");
            };
            total = sum;
        }
        let Some(remaining) = Element::quotient(&reduced.a, &reduced.b, &reduced.c) else {
            panic!("bad remainder");
        };
        let Some(sum) = total.add(&remaining) else {
            panic!("add");
        };
        let Some(difference) = original.sub(&sum) else {
            panic!("sub");
        };
        let Some((a, b, _)) = difference.to_quotient() else {
            panic!("quotient");
        };
        assert!(
            a.is_empty() && b.is_empty(),
            "the Hermite reduction is not an exact identity for D = {d:?}, C = {c:?}"
        );
    }

    #[test]
    fn hermite_reduction_is_exact_for_repeated_and_ramified_denominators() {
        let q = |v: &[i64]| -> QPoly {
            v.iter()
                .map(|&k| BigRational::from_integer(BigInt::from(k)))
                .collect()
        };
        let mul = |a: &[BigRational], b: &[BigRational]| field::mul(a, b);
        let pow = |p: &[BigRational], n: u32| field::pow(p, n);
        let d_cubic = q(&[1, 0, 0, 1]); // x^3 + 1
        let d_quartic = q(&[1, 0, 0, 0, 1]); // x^4 + 1
        let d_ramified = mul(&q(&[1, 0, 1]), &q(&[1, 0, 2])); // (x^2+1)(x^2+2)
        let d_monic = mul(&q(&[1, 0, 0, 1]), &q(&[3, 2, 1])); // (x^3+1)(x^2+2x+3)
        let c_quadratic = q(&[1, 1, 1]); // x^2+x+1
        let c_linear = q(&[-2, 1]); // x - 2
        let c_ramified = mul(&pow(&q(&[1, 0, 1]), 2), &pow(&q(&[-3, 1]), 2)); // ramified at x^2+1
        let c_multi = mul(&pow(&q(&[-1, 1]), 5), &pow(&c_quadratic, 2));
        let c_rational = mul(&pow(&q(&[1, 0, 1]), 3), &q(&[-2, 2]));
        let zero = Vec::new();
        let cases: [(&[BigRational], QPoly, QPoly, QPoly); 7] = [
            (
                &d_cubic[..],
                q(&[1, 0, 1]),
                q(&[2, -1]),
                mul(&pow(&c_quadratic, 2), &c_linear),
            ),
            (&d_quartic[..], q(&[3]), q(&[0, 1]), pow(&c_linear, 3)),
            (&d_ramified[..], q(&[1, 2]), zero.clone(), c_ramified),
            (
                &d_monic[..],
                q(&[1, 0, 0, 2]),
                q(&[0, 1, 1]),
                pow(&mul(&q(&[1, 0, 1]), &c_linear), 3),
            ),
            (&d_cubic[..], zero.clone(), q(&[1]), c_multi),
            (&d_quartic[..], q(&[1, 0, 1]), q(&[0, 2]), c_rational),
            (
                &d_cubic[..],
                q(&[2, 3]),
                q(&[0, 0, 0, 1]),
                pow(&q(&[0, 1]), 4),
            ),
        ];
        for (d, a, b, c) in cases {
            assert_reduction_is_exact(d, &a, &b, &c);
        }
    }

    #[test]
    fn huge_extension_powers_are_refused() {
        let mut graph = Graph::new();
        let Ok(engine) = Engine::install(&mut graph, &[calculus()]) else {
            panic!("install");
        };
        let x = graph.sym("x");
        // Without the exponent guard this builds polynomials of degree
        // billions (or exhausts memory) before any reduction cap applies.
        let Ok(integrand) = graph.parse("sqrt(x^3 + 1)^1000000000") else {
            panic!("parse");
        };
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        assert!(
            integrate_algebraic(&mut cx, integrand, x).is_none(),
            "a billion-fold extension power must be refused"
        );
    }

    #[test]
    fn quartic_reduction_returns_an_identity() {
        let mut graph = Graph::new();
        let Ok(engine) = Engine::install(&mut graph, &[calculus()]) else {
            panic!("install");
        };
        let x = graph.sym("x");
        let Ok(integrand) = graph.parse("sqrt(x^4 + 1)/(x^2 + 1)^2") else {
            panic!("parse");
        };
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        let answer = integrate_algebraic(&mut cx, integrand, x);
        println!(
            "quartic answer: {:?}",
            answer.map(|node| cx.graph.display(node))
        );
        assert!(answer.is_some(), "no reduction produced");
    }

    #[test]
    fn log_part_branches_are_exact() {
        for (c, a) in [
            // Irreducible quadratic residue field: the atan branch.
            (
                field::mul(
                    &[one(), BigRational::zero(), one()],
                    &[
                        BigRational::from_integer(BigInt::from(4)),
                        BigRational::zero(),
                        one(),
                    ],
                ),
                vec![one()],
            ),
            // Linear factors: the logarithm branch.
            (
                field::mul(
                    &[BigRational::from_integer(BigInt::from(-2)), one()],
                    &[BigRational::from_integer(BigInt::from(-3)), one()],
                ),
                vec![one()],
            ),
        ] {
            let mut graph = Graph::new();
            let Ok(engine) = Engine::install(&mut graph, &[calculus()]) else {
                panic!("install");
            };
            let x = graph.sym("x");
            let d = vec![
                one(),
                BigRational::zero(),
                BigRational::zero(),
                BigRational::zero(),
                one(),
            ];
            let b = Vec::new();
            let Some(terms) = logpart::complete(&mut graph, x, &d, &a, &b, &c) else {
                panic!("no logarithmic part for c = {c:?}");
            };
            let num = term::poly_node(&mut graph, x, &a);
            let den = term::poly_node(&mut graph, x, &c);
            let integrand = term::div_node(&mut graph, num, den);
            let Some(derivative) = super::derivative_of(&mut graph, terms, x) else {
                panic!("no derivative");
            };
            let Some(symbol) = graph.symbol_of(x) else {
                panic!("x is not a symbol");
            };
            let mut checked = 0_u32;
            for point in [0.37, 0.91, 1.63, 2.3, -0.7] {
                let mut env = Env::numeric(0.0);
                env.bind(symbol, point);
                let (Some(want), Some(got)) =
                    (graph.eval(integrand, &env), graph.eval(derivative, &env))
                else {
                    continue;
                };
                checked = checked.saturating_add(1);
                assert!(
                    (want - got).abs() < 1e-9 * (1.0 + want.abs()),
                    "at {point}: integral {want}, derivative {got}"
                );
            }
            assert!(checked >= 3, "nothing could be evaluated");
            let _ = engine;
        }
    }

    fn show(f: &str) -> (String, bool) {
        reduce_with(&[calculus()], &format!("integral({f}, x)"), &[])
    }

    #[test]
    fn smoke() {
        for f in [
            "1/sqrt(x^2 + 1)",
            "x/sqrt(x^2 + 1)",
            "sqrt(x^2 + 1)",
            "1/((x^2 + 1)*sqrt(x^2 + 1))",
            "sqrt(1 - x^2)",
            "x^2/sqrt(1 - x^2)",
            "1/((x - 2)*sqrt(x^2 + 1))",
            "(x + 1)/((x^2 + 1)^2*sqrt(x^2 + 1))",
            "sqrt(x)",
            "x^3/sqrt(x^4 + 1)",
            "sqrt(x^4 + 1)/(x^2 + 1)^2",
            "sqrt(x^4 + 1)/((x - 2)*(x - 3))",
            "sqrt(x^4 - 2*x^2 + 2)/((x - 2)*(x - 3))",
        ] {
            let (text, reduced) = show(f);
            println!("{f}  ->  {text}   [{reduced}]");
        }
    }
}
