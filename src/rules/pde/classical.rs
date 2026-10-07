//! More classical solution formulas.
//!
//! * **Parabolic equations with drift, reaction and sources** on the line
//!   and half-line: `a u_t = b u_xx + c u_x + d u + f(x, t)` with constant
//!   coefficients (`D = b/a > 0`). The gauge `u = exp(α x + β t) v` with
//!   `α = -C/(2D)`, `β = R - C²/(4D)` (`C = c/a`, `R = d/a`) leaves
//!   `v_t = D v_xx + exp(-α x - β t) f/a`; then the heat kernel for the
//!   initial data, Duhamel's principle for the source, and on `x > 0` with
//!   `u(0, t) = 0` (or `u_x(0, t) = 0` when there is no drift) the method of
//!   images.
//! * **Waves on the half-line** `x > 0`, `u(0, t) = 0` or `u_x(0, t) = 0`:
//!   d'Alembert's formula for the odd or even extension of the data.
//! * **Laplace's equation in a disk and in a ball** (`laplace_disk(f, R, r,
//!   θ)`, `laplace_ball(f, R, r, θ)` for axisymmetric data): the Fourier
//!   series `Σ (r/R)^n (a_n cos nθ + b_n sin nθ)` and the Legendre series
//!   `Σ A_n (r/R)^n P_n(cos θ)`, with exact coefficients; data that is a
//!   trigonometric polynomial (or a polynomial in `cos θ`) gives a finite
//!   sum.

use super::Conditions;
use super::Problem;
use super::dummy;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

fn div(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let inverse = powi(cx.graph, b, -1);
    mul(cx.graph, &[a, inverse])
}

/// The 1-D heat kernel `(4π D τ)^(-1/2) exp(-(ξ)²/(4 D τ))`.
pub(super) fn kernel(
    cx: &mut Cx<'_>,
    d: NodeId,
    xi: NodeId,
    tau: NodeId,
) -> Option<NodeId> {
    let pi = cx.graph.ops().lookup("pi")?;
    let exp = cx.graph.ops().lookup("exp")?;
    let pi = cx.graph.node(pi, &[]);
    let four = cx.graph.int(4);
    let four_dt = mul(cx.graph, &[four, d, tau]);
    let square = powi(cx.graph, xi, 2);
    let ratio = div(cx, square, four_dt);
    let exponent = neg(cx.graph, ratio);
    let gaussian = cx.graph.node(exp, &[exponent]);
    let base = mul(cx.graph, &[pi, four_dt]);
    let half = cx.graph.num(Number::fraction(-1, 2)?);
    let normalisation = pow(cx.graph, base, half);
    Some(mul(cx.graph, &[normalisation, gaussian]))
}

/// `u_tt = c² u_xx` on `x > 0` with `u(0, t) = 0` (odd extension) or
/// `u_x(0, t) = 0` (even extension), `u(x, 0) = g`, `u_t(x, 0) = h`.
pub(super) fn wave_half_line(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = super::evolution(cx, p)?;
    if e.order != 2 || p.dimension() != 2 || !p.homogeneous(cx.graph) || !cx.is_zero(e.potential) {
        return None;
    }
    let time = e.time;
    let space = 1 - time;
    let (x, t) = (p.vars[space], p.vars[time]);
    let zero_index = vec![0; 2];
    let displacement = conditions.find(cx.graph, time, &zero_index, None)?;
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None);
    let dirichlet = conditions.find(cx.graph, space, &zero_index, None);
    let neumann = conditions.find(cx.graph, space, &p.unit(space, 1), None);
    let odd = match (dirichlet, neumann) {
        | (Some(c), None) if cx.is_zero(c.value) && cx.is_zero(c.point) => true,
        | (None, Some(c)) if cx.is_zero(c.value) && cx.is_zero(c.point) => false,
        | _ => return None,
    };
    if conditions.0.len() != 2 + usize::from(velocity.is_some()) {
        return None;
    }
    let (sign, abs) = (cx.graph.ops().lookup("sign")?, cx.graph.ops().lookup("abs")?);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let speed = pow(cx.graph, e.speed, half);
    let speed = cx.simplify(speed);
    let ct = mul(cx.graph, &[speed, t]);
    // The extension of f evaluated at z: f(|z|), times sign(z) if odd.
    let extend = |cx: &mut Cx<'_>, f: NodeId, z: NodeId| -> NodeId {
        let magnitude = cx.graph.node(abs, &[z]);
        let value = cx.graph.substitute(f, x, magnitude);
        if odd {
            let s = cx.graph.node(sign, &[z]);
            mul(cx.graph, &[s, value])
        } else {
            value
        }
    };
    let plus = add(cx.graph, &[x, ct]);
    let minus = sub(cx.graph, x, ct);
    let (g_plus, g_minus) = (extend(cx, displacement.value, plus), extend(cx, displacement.value, minus));
    let sum = add(cx.graph, &[g_plus, g_minus]);
    let mut u = mul(cx.graph, &[half, sum]);
    if let Some(v) = velocity
        && !cx.is_zero(v.value) {
            // (1/2c) ∫_{x-ct}^{x+ct} h_ext(s) ds
            let defint = cx.graph.ops().lookup("defint")?;
            let (s, _) = dummy(cx, p, "s");
            let h_at = cx.graph.substitute(v.value, x, s);
            let h_ext = {
                let h = cx.graph.substitute(h_at, s, x);
                extend(cx, h, s)
            };
            let integral = cx.graph.node(defint, &[h_ext, s, minus, plus]);
            let two_c = {
                let two = cx.graph.int(2);
                mul(cx.graph, &[two, speed])
            };
            let part = div(cx, integral, two_c);
            u = add(cx.graph, &[u, part]);
        }
    Some(cx.simplify(u))
}

/// `laplace_disk(f(θ), R, r, θ)`: the harmonic function in the disk of
/// radius `R` with boundary values `f`.
pub(super) fn laplace_disk(
    cx: &mut Cx<'_>,
    f: NodeId,
    radius: NodeId,
    r: NodeId,
    theta: NodeId,
) -> Option<NodeId> {
    let (sin, cos, pi, defint, sum) = (
        cx.graph.ops().lookup("sin")?,
        cx.graph.ops().lookup("cos")?,
        cx.graph.ops().lookup("pi")?,
        cx.graph.ops().lookup("defint")?,
        cx.graph.ops().lookup("sum")?,
    );
    let pi = cx.graph.node(pi, &[]);
    let two = cx.graph.int(2);
    let two_pi = mul(cx.graph, &[two, pi]);
    let zero = cx.graph.int(0);
    // Coefficients a_n = (1/π) ∫ f cos(nθ), b_n likewise; n symbolic first
    // to detect a finite expansion.
    let coefficient = |cx: &mut Cx<'_>, wave: crate::graph::OpId, n: NodeId| -> NodeId {
        let n_theta = mul(cx.graph, &[n, theta]);
        let w = cx.graph.node(wave, &[n_theta]);
        let integrand = mul(cx.graph, &[f, w]);
        let integral = cx.graph.node(defint, &[integrand, theta, zero, two_pi]);
        let scaled = div(cx, integral, pi);
        cx.simplify(scaled)
    };
    let ratio = div(cx, r, radius);
    // Try a finite trigonometric polynomial: harmonics up to 12.
    let mut terms = Vec::new();
    let a0 = {
        let integral = cx.graph.node(defint, &[f, theta, zero, two_pi]);
        let scaled = div(cx, integral, two_pi);
        cx.simplify(scaled)
    };
    let symbolic = |cx: &Cx<'_>, e: NodeId| cx.graph.ops().lookup("defint").is_some_and(|d| crate::rules::ode::occurs_op(cx.graph, e, d));
    if symbolic(cx, a0) {
        return None;
    }
    terms.push(a0);
    let mut last_nonzero = 0;
    for k in 1..=12_i64 {
        let n = cx.graph.int(k);
        let (a, b) = (coefficient(cx, cos, n), coefficient(cx, sin, n));
        if symbolic(cx, a) || symbolic(cx, b) {
            return None;
        }
        let n_theta = mul(cx.graph, &[n, theta]);
        let (c, s) = (cx.graph.node(cos, &[n_theta]), cx.graph.node(sin, &[n_theta]));
        let rn = powi(cx.graph, ratio, k);
        let part = {
            let ac = mul(cx.graph, &[a, c]);
            let bs = mul(cx.graph, &[b, s]);
            let both = add(cx.graph, &[ac, bs]);
            mul(cx.graph, &[rn, both])
        };
        let part = cx.simplify(part);
        if !cx.is_zero(part) {
            last_nonzero = k;
        }
        terms.push(part);
    }
    if last_nonzero < 12 {
        let total = add(cx.graph, &terms);
        return Some(cx.simplify(total));
    }
    // An infinite series with symbolic coefficients.
    let n_symbol = cx.graph.interner_mut().fresh_symbol("n");
    let n = cx.graph.symbol_node(n_symbol);
    cx.graph.assume(n_symbol, crate::graph::Facts::POSITIVE | crate::graph::Facts::INTEGER);
    let (a, b) = (coefficient(cx, cos, n), coefficient(cx, sin, n));
    let n_theta = mul(cx.graph, &[n, theta]);
    let (c, s) = (cx.graph.node(cos, &[n_theta]), cx.graph.node(sin, &[n_theta]));
    let rn = pow(cx.graph, ratio, n);
    let general = {
        let ac = mul(cx.graph, &[a, c]);
        let bs = mul(cx.graph, &[b, s]);
        let both = add(cx.graph, &[ac, bs]);
        mul(cx.graph, &[rn, both])
    };
    let one = cx.graph.int(1);
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let series = cx.graph.node(sum, &[general, n, one, oo]);
    Some(add(cx.graph, &[a0, series]))
}

/// `laplace_ball(f(θ), R, r, θ)`: the axisymmetric harmonic function in
/// the ball with boundary values `f(θ)` (polar angle θ), as a Legendre
/// series in `cos θ`; finite when `f` is a polynomial in `cos θ`.
pub(super) fn laplace_ball(
    cx: &mut Cx<'_>,
    f: NodeId,
    radius: NodeId,
    r: NodeId,
    theta: NodeId,
) -> Option<NodeId> {
    let (cos, sin, legendre, defint, pi) = (
        cx.graph.ops().lookup("cos")?,
        cx.graph.ops().lookup("sin")?,
        cx.graph.ops().lookup("legendre")?,
        cx.graph.ops().lookup("defint")?,
        cx.graph.ops().lookup("pi")?,
    );
    let pi = cx.graph.node(pi, &[]);
    let zero = cx.graph.int(0);
    let ratio = div(cx, r, radius);
    let cos_theta = cx.graph.node(cos, &[theta]);
    let sin_theta = cx.graph.node(sin, &[theta]);
    let mut terms = Vec::new();
    let mut last_nonzero = None;
    for n in 0..=10_i64 {
        let n_node = cx.graph.int(n);
        let p_n = cx.graph.node(legendre, &[n_node, cos_theta]);
        let integrand = mul(cx.graph, &[f, p_n, sin_theta]);
        let integral = cx.graph.node(defint, &[integrand, theta, zero, pi]);
        let scale = cx.graph.num(Number::fraction(2 * n + 1, 2)?);
        let a = mul(cx.graph, &[scale, integral]);
        let a = cx.simplify(a);
        if crate::rules::ode::occurs_op(cx.graph, a, defint) {
            return None;
        }
        let rn = powi(cx.graph, ratio, n);
        let part = mul(cx.graph, &[a, rn, p_n]);
        let part = cx.simplify(part);
        if !cx.is_zero(part) {
            last_nonzero = Some(n);
        }
        terms.push(part);
    }
    if last_nonzero? >= 10 {
        return None;
    }
    let total = add(cx.graph, &terms);
    Some(cx.simplify(total))
}
