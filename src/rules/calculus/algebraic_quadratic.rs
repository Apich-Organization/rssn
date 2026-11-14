//! Closed forms for a quadratic radicand `D = a x^2 + b x + c`.
//!
//! A quadratic radicand makes `Q(x)(y)` the function field of a conic, so
//! every element has an elementary integral.  The conic is parametrised
//! rationally (over a real quadratic field, which is rational whenever the
//! constants needed are squares of rationals):
//!
//! * `a x^2 + b x + c` with two rational roots: `t = y/(x - r1)`;
//! * `a > 0` and negative discriminant: `t = y - sqrt(a) x`;
//! * `a < 0`: `x = h + m sin(theta)`, `y = k cos(theta)`, with
//!   `h = -b/(2a)`, `m^2 = (b^2 - 4ac)/(4a^2)`, `k = sqrt(-a) m`.
//!
//! The transformed integrand is a rational function of the parameter (or a
//! rational function of `sin` and `cos`), which the ordinary integrator
//! handles; the result is mapped back and checked by the caller.
//!
//! Only the cases in which the parametrisation is rational are handled
//! here.  A radicand such as `2 x^2 + 3` needs `sqrt(2)` in the
//! parametrisation; the ordinary integrator's quadratic-radical rules take
//! those, so this module returns `None` and lets them.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;

use super::algebraic_field as field;
use super::algebraic_term as term;

/// `p/q` as an exact rational.
fn rat(
    p: i64,
    q: i64,
) -> BigRational {
    BigRational::new(BigInt::from(p), BigInt::from(q))
}

/// The rational number `q` as a term.
fn constant(
    graph: &mut Graph,
    q: &BigRational,
) -> NodeId {
    term::rat_node(graph, q)
}

/// The integrand `(A(x0) + B(x0) y0) / C(x0)`.
fn element_at(
    graph: &mut Graph,
    x0: NodeId,
    y0: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> NodeId {
    let an = term::poly_node(graph, x0, a);
    let bn = term::poly_node(graph, x0, b);
    let cn = term::poly_node(graph, x0, c);
    let by = graph.node(crate::graph::op::core::MUL, &[bn, y0]);
    let body = graph.node(crate::graph::op::core::ADD, &[an, by]);
    term::div_node(graph, body, cn)
}

/// A fresh parameter symbol.
fn parameter(
    graph: &mut Graph,
    name: &str,
) -> NodeId {
    let s = graph.interner_mut().fresh_symbol(name);
    graph.symbol_node(s)
}

/// Integrates `integrand` with respect to the fresh symbol `variable`.
fn integrate(
    cx: &mut Cx<'_>,
    integrand: NodeId,
    variable: NodeId,
) -> Option<NodeId> {
    crate::rules::calculus::antiderivative(cx, integrand, variable)
}

/// `∫ (a + b y)/c dx` for `D = a x^2 + b x + c` with two rational roots.
fn two_rational_roots(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    d: &[BigRational],
    root: &BigRational,
) -> Option<NodeId> {
    let (b1, a2) = (d.get(1)?, d.get(2)?);
    if a2.is_zero() {
        return None;
    }
    let two_a = BigRational::from_integer(BigInt::from(2)) * a2;
    let r1 = (-b1 + root) / &two_a;
    let r2 = (-b1 - root) / &two_a;
    let t = parameter(cx.graph, "t");
    let t2 = {
        let two = cx.graph.int(2);
        cx.graph.node(crate::graph::op::core::POW, &[t, two])
    };
    let den = {
        let a_node = constant(cx.graph, a2);
        term::sub_node(cx.graph, t2, a_node)
    };
    // x(t) = (r1 t^2 - a r2) / (t^2 - a)
    let x_t = {
        let r1_node = constant(cx.graph, &r1);
        let num1 = cx.graph.node(crate::graph::op::core::MUL, &[r1_node, t2]);
        let ar2 = constant(cx.graph, &(a2 * &r2));
        let num = term::sub_node(cx.graph, num1, ar2);
        term::div_node(cx.graph, num, den)
    };
    // y(t) = sqrt(disc) t / (t^2 - a)
    let y_t = {
        let s_node = constant(cx.graph, root);
        let num = cx.graph.node(crate::graph::op::core::MUL, &[s_node, t]);
        term::div_node(cx.graph, num, den)
    };
    // dx/dt = -2 sqrt(disc) t / (t^2 - a)^2
    let jacobian = {
        let den2 = {
            let two = cx.graph.int(2);
            cx.graph.node(crate::graph::op::core::POW, &[den, two])
        };
        let minus_two_s = constant(cx.graph, &(-rat(2, 1) * root));
        let num = cx
            .graph
            .node(crate::graph::op::core::MUL, &[minus_two_s, t]);
        term::div_node(cx.graph, num, den2)
    };
    let body = element_at(cx.graph, x_t, y_t, a, b, c);
    let integrand = cx
        .graph
        .node(crate::graph::op::core::MUL, &[body, jacobian]);
    let primitive = integrate(cx, integrand, t)?;
    // t = y / (x - r1)
    let r1_node = constant(cx.graph, &r1);
    let shift = term::sub_node(cx.graph, x, r1_node);
    let radicand = term::poly_node(cx.graph, x, d);
    let sqrt_d = term::sqrt_node(cx.graph, radicand);
    let t_of_x = term::div_node(cx.graph, sqrt_d, shift);
    Some(cx.graph.substitute(primitive, t, t_of_x))
}

/// `∫ (a + b y)/c dx` for `D = a x^2 + b x + c` with `a > 0` and negative
/// discriminant: `t = y - sqrt(a) x`.
fn negative_discriminant(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    d: &[BigRational],
    s: &BigRational,
) -> Option<NodeId> {
    let b1 = d.get(1)?;
    let c0 = d.first()?;
    let t = parameter(cx.graph, "t");
    let s_node = constant(cx.graph, s);
    let t2 = {
        let two = cx.graph.int(2);
        cx.graph.node(crate::graph::op::core::POW, &[t, two])
    };
    let two_s_t = {
        let two = cx.graph.int(2);
        let st = cx.graph.node(crate::graph::op::core::MUL, &[s_node, t]);
        cx.graph.node(crate::graph::op::core::MUL, &[two, st])
    };
    let b_node = constant(cx.graph, b1);
    let den = term::sub_node(cx.graph, two_s_t, b_node);
    // x(t) = (c - t^2)/(2 s t - b)
    let c_node = constant(cx.graph, c0);
    let shifted_t = term::sub_node(cx.graph, c_node, t2);
    let x_t = term::div_node(cx.graph, shifted_t, den);
    // y(t) = (s t^2 - b t + s c)/(2 s t - b)
    let y_t = {
        let s_t2 = cx.graph.node(crate::graph::op::core::MUL, &[s_node, t2]);
        let bt = cx.graph.node(crate::graph::op::core::MUL, &[b_node, t]);
        let sc = constant(cx.graph, &(s * c0));
        let partial = term::sub_node(cx.graph, s_t2, bt);
        let num = cx.graph.node(crate::graph::op::core::ADD, &[partial, sc]);
        term::div_node(cx.graph, num, den)
    };
    // dx/dt = 2 (b t - s t^2 - s c)/(2 s t - b)^2
    let jacobian = {
        let den2 = {
            let two = cx.graph.int(2);
            cx.graph.node(crate::graph::op::core::POW, &[den, two])
        };
        let bt = cx.graph.node(crate::graph::op::core::MUL, &[b_node, t]);
        let s_t2 = cx.graph.node(crate::graph::op::core::MUL, &[s_node, t2]);
        let sc = constant(cx.graph, &(s * c0));
        let first = term::sub_node(cx.graph, bt, s_t2);
        let inner = term::sub_node(cx.graph, first, sc);
        let two = cx.graph.int(2);
        let num = cx.graph.node(crate::graph::op::core::MUL, &[two, inner]);
        term::div_node(cx.graph, num, den2)
    };
    let body = element_at(cx.graph, x_t, y_t, a, b, c);
    let integrand = cx
        .graph
        .node(crate::graph::op::core::MUL, &[body, jacobian]);
    let primitive = integrate(cx, integrand, t)?;
    // t = y - s x
    let s_x = cx.graph.node(crate::graph::op::core::MUL, &[s_node, x]);
    let root_arg = term::poly_node(cx.graph, x, d);
    let sqrt_d = term::sqrt_node(cx.graph, root_arg);
    let t_of_x = term::sub_node(cx.graph, sqrt_d, s_x);
    Some(cx.graph.substitute(primitive, t, t_of_x))
}

/// `∫ (a + b y)/c dx` for `D = a x^2 + b x + c` with `a < 0`: the
/// trigonometric parametrisation `x = h + m sin(theta)`.
#[allow(clippy::too_many_arguments)] // radicand, integrand and parametrisation constants
fn circle(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    d: &[BigRational],
    m: &BigRational,
    k: &BigRational,
) -> Option<NodeId> {
    let (b1, a2) = (d.get(1)?, d.get(2)?);
    let two_a = BigRational::from_integer(BigInt::from(2)) * a2;
    let h = (-b1) / &two_a;
    let theta = parameter(cx.graph, "theta");
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    let sin_theta = cx.graph.node(sin, &[theta]);
    let cos_theta = cx.graph.node(cos, &[theta]);
    let h_node = constant(cx.graph, &h);
    let m_node = constant(cx.graph, m);
    let m_sin = cx
        .graph
        .node(crate::graph::op::core::MUL, &[m_node, sin_theta]);
    let x_t = cx.graph.node(crate::graph::op::core::ADD, &[h_node, m_sin]);
    let k_node = constant(cx.graph, k);
    let y_t = cx
        .graph
        .node(crate::graph::op::core::MUL, &[k_node, cos_theta]);
    let jacobian = cx
        .graph
        .node(crate::graph::op::core::MUL, &[m_node, cos_theta]);
    let body = element_at(cx.graph, x_t, y_t, a, b, c);
    let integrand = cx
        .graph
        .node(crate::graph::op::core::MUL, &[body, jacobian]);
    let primitive = integrate(cx, integrand, theta)?;
    // theta = asin((x - h)/m)
    let shifted = term::sub_node(cx.graph, x, h_node);
    let scaled = term::div_node(cx.graph, shifted, m_node);
    let arc = term::asin_node(cx.graph, scaled)?;
    Some(cx.graph.substitute(primitive, theta, arc))
}

/// Integrates `(a + b y)/c` for a quadratic `D` when the parametrisation is
/// rational, or returns `None`.
pub(super) fn integrate_quadratic(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    d: &[BigRational],
) -> Option<NodeId> {
    let (c0, b1, a2) = (d.first()?, d.get(1)?, d.get(2)?);
    if a2.is_zero() {
        return None;
    }
    let four = BigRational::from_integer(BigInt::from(4));
    let discriminant = b1 * b1 - &four * a2 * c0;
    if discriminant.is_zero() {
        return None;
    }
    if discriminant.is_positive() {
        // Rational roots: the parametrisation `t = y/(x - r1)` is
        // rational over `Q`, whatever the sign of the leading coefficient.
        if let Some(root) = field::rational_sqrt(&discriminant) {
            return two_rational_roots(cx, x, a, b, c, d, &root);
        }
        if a2.is_negative() {
            // The radicand is positive between two real roots: use the
            // trigonometric parametrisation when it is rational.
            // D = a u^2 + (4ac - b^2)/(4a) = |a| (m^2 - u^2) with
            // m^2 = (b^2 - 4ac)/(4 a^2).
            let m_squared = &discriminant / (&four * a2 * a2);
            let m = field::rational_sqrt(&m_squared)?;
            let k = &m * field::rational_sqrt(&(-a2))?;
            if m.is_zero() || !k.is_positive() {
                return None;
            }
            return circle(cx, x, a, b, c, d, &m, &k);
        }
        return None;
    }
    // `a > 0`, negative discriminant.
    let s = field::rational_sqrt(a2)?;
    negative_discriminant(cx, x, a, b, c, d, &s)
}

/// Integrates `(a + b y)/c` for a linear radicand `D = b0 + b1 x` by
/// `t = y`.
pub(super) fn integrate_linear(
    cx: &mut Cx<'_>,
    x: NodeId,
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    d: &[BigRational],
) -> Option<NodeId> {
    let b1 = d.get(1)?;
    if b1.is_zero() {
        return None;
    }
    let b0 = d.first()?;
    let t = parameter(cx.graph, "t");
    let t2 = {
        let two = cx.graph.int(2);
        cx.graph.node(crate::graph::op::core::POW, &[t, two])
    };
    let b0_node = constant(cx.graph, b0);
    let b1_node = constant(cx.graph, b1);
    let shifted = term::sub_node(cx.graph, t2, b0_node);
    let x_t = term::div_node(cx.graph, shifted, b1_node);
    let body = element_at(cx.graph, x_t, t, a, b, c);
    // dx = 2 t / b1 dt
    let two = cx.graph.int(2);
    let jacobian = {
        let num = cx.graph.node(crate::graph::op::core::MUL, &[two, t]);
        term::div_node(cx.graph, num, b1_node)
    };
    let integrand = cx
        .graph
        .node(crate::graph::op::core::MUL, &[body, jacobian]);
    let primitive = integrate(cx, integrand, t)?;
    let radicand = term::poly_node(cx.graph, x, d);
    let sqrt_d = term::sqrt_node(cx.graph, radicand);
    Some(cx.graph.substitute(primitive, t, sqrt_d))
}
