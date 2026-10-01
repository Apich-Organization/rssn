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
use crate::graph::Env;
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

/// The coefficients `(a, b, c, d)` and source of a 1+1 parabolic
/// equation, all constant except the source.
struct Parabolic {
    time: usize,
    space: usize,
    diffusivity: NodeId,
    drift: NodeId,
    reaction: NodeId,
    /// `f / a`.
    forcing: NodeId,
}

fn parabolic(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Parabolic> {
    if p.dimension() != 2 || p.nonlinear {
        return None;
    }
    let time = p.time(cx.graph);
    let space = 1 - time;
    let ut = p.unit(time, 1);
    let (ux, uxx, u) = (p.unit(space, 1), p.unit(space, 2), vec![0; 2]);
    if !p.only(cx.graph, &[ut.clone(), ux.clone(), uxx.clone(), u.clone()]) {
        return None;
    }
    let a = p.coefficient(cx.graph, &ut);
    if cx.is_zero(a) {
        return None;
    }
    // a u_t + b' u_xx + c' u_x + d' u + source = 0
    let coefficient = |cx: &mut Cx<'_>, index: &[u32]| {
        let c = p.coefficient(cx.graph, index);
        let ratio = div(cx, c, a);
        let ratio = neg(cx.graph, ratio);
        cx.simplify(ratio)
    };
    let diffusivity = coefficient(cx, &uxx);
    let drift = coefficient(cx, &ux);
    let reaction = coefficient(cx, &u);
    if [diffusivity, drift, reaction].iter().any(|&c| !p.constant(cx.graph, c)) {
        return None;
    }
    if !cx.graph.eval(diffusivity, &Env::numeric(0.0)).is_none_or(|v| v > 0.0) || cx.is_zero(diffusivity) {
        return None;
    }
    let forcing = {
        let ratio = div(cx, p.source, a);
        let negated = neg(cx.graph, ratio);
        cx.simplify(negated)
    };
    Some(Parabolic { time, space, diffusivity, drift, reaction, forcing })
}

/// The 1-D heat kernel `(4π D τ)^(-1/2) exp(-(ξ)²/(4 D τ))`.
fn kernel(
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

/// Initial-value problems for drift–diffusion–reaction equations on the
/// line, or on `x > 0` with a homogeneous boundary condition at 0.
pub(super) fn drift_diffusion(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = parabolic(cx, p)?;
    let (x, t) = (p.vars[e.space], p.vars[e.time]);
    let zero_index = vec![0; 2];
    let initial = conditions.find(cx.graph, e.time, &zero_index, None)?;
    if !cx.graph.number_of(initial.point).is_some_and(Number::is_zero) {
        return None;
    }
    // A boundary condition at x = 0 selects the half-line.
    let dirichlet = conditions.find(cx.graph, e.space, &zero_index, None);
    let neumann = conditions.find(cx.graph, e.space, &p.unit(e.space, 1), None);
    let boundary_sign: Option<i64> = match (dirichlet, neumann) {
        | (Some(c), None) if cx.is_zero(c.value) && cx.graph.number_of(c.point).is_some_and(Number::is_zero) => Some(-1),
        | (None, Some(c)) if cx.is_zero(c.value) && cx.graph.number_of(c.point).is_some_and(Number::is_zero) => Some(1),
        | (None, None) => None,
        | _ => return None,
    };
    let expected = 1 + usize::from(boundary_sign.is_some());
    if conditions.0.len() != expected {
        return None;
    }
    // Images need the gauge to be even in x: no drift.
    if boundary_sign.is_some() && !cx.is_zero(e.drift) {
        return None;
    }
    let exp = cx.graph.ops().lookup("exp")?;
    let defint = cx.graph.ops().lookup("defint")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let d = e.diffusivity;
    // α = -C/(2D), β = R - C²/(4D)
    let two = cx.graph.int(2);
    let four = cx.graph.int(4);
    let alpha = {
        let two_d = mul(cx.graph, &[two, d]);
        let r = div(cx, e.drift, two_d);
        let r = neg(cx.graph, r);
        cx.simplify(r)
    };
    let beta = {
        let c2 = powi(cx.graph, e.drift, 2);
        let four_d = mul(cx.graph, &[four, d]);
        let r = div(cx, c2, four_d);
        let r = sub(cx.graph, e.reaction, r);
        cx.simplify(r)
    };
    let (s, _) = dummy(cx, p, "s");
    // v(x, 0) = exp(-α x) g(x) at the dummy point s.
    let g_at_s = cx.graph.substitute(initial.value, x, s);
    let minus_alpha_s = {
        let product = mul(cx.graph, &[alpha, s]);
        neg(cx.graph, product)
    };
    let gauge_s = cx.graph.node(exp, &[minus_alpha_s]);
    let data = mul(cx.graph, &[gauge_s, g_at_s]);
    let x_minus_s = sub(cx.graph, x, s);
    let mut weight = kernel(cx, d, x_minus_s, t)?;
    let lower = if let Some(sign) = boundary_sign {
        let x_plus_s = add(cx.graph, &[x, s]);
        let image = kernel(cx, d, x_plus_s, t)?;
        let sign = cx.graph.int(sign);
        let image = mul(cx.graph, &[sign, image]);
        weight = add(cx.graph, &[weight, image]);
        cx.graph.int(0)
    } else {
        neg(cx.graph, oo)
    };
    let integrand = mul(cx.graph, &[weight, data]);
    let mut v = cx.graph.node(defint, &[integrand, s, lower, oo]);
    // Duhamel: ∫_0^t ∫ K(x - s, t - τ) f̃(s, τ) ds dτ.
    if !cx.is_zero(e.forcing) {
        let (tau, _) = dummy(cx, p, "tau");
        let source = cx.graph.substitute(e.forcing, x, s);
        let source = cx.graph.substitute(source, t, tau);
        let phase = {
            let ax = mul(cx.graph, &[alpha, s]);
            let bt = mul(cx.graph, &[beta, tau]);
            let sum = add(cx.graph, &[ax, bt]);
            neg(cx.graph, sum)
        };
        let gauge = cx.graph.node(exp, &[phase]);
        let elapsed = sub(cx.graph, t, tau);
        let mut w = kernel(cx, d, x_minus_s, elapsed)?;
        if let Some(sign) = boundary_sign {
            let x_plus_s = add(cx.graph, &[x, s]);
            let image = kernel(cx, d, x_plus_s, elapsed)?;
            let sign = cx.graph.int(sign);
            let image = mul(cx.graph, &[sign, image]);
            w = add(cx.graph, &[w, image]);
        }
        let body = mul(cx.graph, &[w, gauge, source]);
        let inner = cx.graph.node(defint, &[body, s, lower, oo]);
        let zero = cx.graph.int(0);
        let duhamel = cx.graph.node(defint, &[inner, tau, zero, t]);
        v = add(cx.graph, &[v, duhamel]);
    }
    let phase = {
        let ax = mul(cx.graph, &[alpha, x]);
        let bt = mul(cx.graph, &[beta, t]);
        add(cx.graph, &[ax, bt])
    };
    let gauge = cx.graph.node(exp, &[phase]);
    let u = mul(cx.graph, &[gauge, v]);
    Some(cx.simplify(u))
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
        | (Some(c), None) if cx.is_zero(c.value) => true,
        | (None, Some(c)) if cx.is_zero(c.value) => false,
        | _ => return None,
    };
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
    if let Some(v) = velocity {
        if !cx.is_zero(v.value) {
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
