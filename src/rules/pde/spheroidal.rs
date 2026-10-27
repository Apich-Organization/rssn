//! Helmholtz and Laplace problems in prolate and oblate spheroidal
//! coordinates.
//!
//! The equation, in coordinates `(ξ, η, φ)` (or `(ζ, η, φ)` for the oblate
//! system, `(ξ, η)` without azimuthal dependence), multiplied by the common
//! scale factor,
//!
//! * prolate: `(ξ²-1) u_ξξ + 2ξ u_ξ + (1-η²) u_ηη - 2η u_η
//!   + (ξ²-η²)/((ξ²-1)(1-η²)) u_φφ + c² (ξ²-η²) u = 0`,
//! * oblate:  `(ζ²+1) u_ζζ + 2ζ u_ζ + (1-η²) u_ηη - 2η u_η
//!   + (ζ²+η²)/((ζ²+1)(1-η²)) u_φφ + c² (ζ²+η²) u = 0`,
//!
//! is recognised from its coefficients (any common factor, `c²` constant,
//! either sign). It separates into the angular equation
//! `(1-η²) S'' - 2η S' + (λ - γ η² - m²/(1-η²)) S = 0` with `γ = c²`
//! (prolate) or `γ = -c²` (oblate), solved by `spheroidal_eigen` /
//! `spheroidal_mode`, and the radial equation
//!
//! * prolate: `(ξ²-1) R'' + 2ξ R' - (λ - γ ξ² + m²/(ξ²-1)) R = 0`, regular at
//!   the focal segment `ξ = 1` (`R ~ (ξ²-1)^{m/2}`), integrated in
//!   `ξ = cosh τ`;
//! * oblate: `(ζ²+1) R'' + 2ζ R' - (λ + γ ζ² - m²/(ζ²+1)) R = 0` with the
//!   parity of `S` at the focal disk `ζ = 0`, integrated in `ζ = sinh σ`.
//!
//! The Dirichlet problem inside the spheroid `ξ < ξ₀` (`ζ < ζ₀`) with data
//! `u(ξ₀, η, φ) = g(η, φ)` is solved by the finite expansion
//! `u = Σ (a_mn cos mφ + b_mn sin mφ) S_mn(η) R_mn(ξ)/R_mn(ξ₀)` whose
//! coefficients are computed by quadrature; the expansion is returned when
//! it reproduces the data at sample points (so the data must be well
//! represented by the first 16 angular functions of every order `m ≤ 6`).
//! Exterior regions (radiation condition) are not covered.

use super::Conditions;
use super::Problem;
use super::coulomb::spheroidal_pair;
use super::util::call;
use super::util::div;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

const STEP: f64 = 0.002;
const T_MAX: f64 = 4.5;
const MODES: usize = 16;
const M_MAX: usize = 6;

type RadialTable = Rc<Vec<f64>>;

thread_local! {
    static RADIAL: RefCell<Vec<([u64; 4], RadialTable)>> = const { RefCell::new(Vec::new()) };
}

/// The regular radial function on the uniform grid `τ_i = i STEP`
/// (`ξ = cosh τ` prolate, `ζ = sinh τ` oblate).
fn radial_table(
    prolate: bool,
    m: usize,
    n: usize,
    gamma: f64,
) -> Option<RadialTable> {
    let key = [u64::from(prolate), (m * 1000 + n) as u64, gamma.to_bits(), 0];
    if let Some(found) = RADIAL.with(|c| c.borrow().iter().find(|(k, _)| *k == key).map(|(_, v)| Rc::clone(v))) {
        return Some(found);
    }
    let lambda = *spheroidal_pair(m, gamma).0.get(n.checked_sub(m)?)?;
    let steps = (T_MAX / STEP) as usize;
    let mf = m as f64;
    let mut table = vec![0.0; steps + 1];
    let (mut y, mut p, start);
    if prolate {
        // Series at the focal segment: R = τ^m (1 + a τ²), ξ = cosh τ.
        let a = ((lambda - gamma) - mf * (mf + 1.0) / 3.0) / (4.0 * mf + 4.0);
        start = 10;
        let tau0 = STEP * start as f64;
        for (i, slot) in table.iter_mut().enumerate().take(start) {
            let tau = STEP * i as f64;
            *slot = tau.powi(m as i32) * (1.0 + a * tau * tau);
        }
        y = tau0.powi(m as i32) * (1.0 + a * tau0 * tau0);
        p = tau0.powi(m as i32 - 1) * (mf * (1.0 + a * tau0 * tau0) + 2.0 * a * tau0 * tau0);
    } else {
        // The focal disk: parity of the angular function.
        start = 0;
        (y, p) = if (n - m) % 2 == 0 { (1.0, 0.0) } else { (0.0, 1.0) };
    }
    table[start] = y;
    let potential = |t: f64| -> f64 {
        let (s, c) = (t.sinh(), t.cosh());
        if prolate { lambda - gamma * c * c + mf * mf / (s * s) } else { lambda + gamma * s * s - mf * mf / (c * c) }
    };
    let damping = |t: f64| -> f64 { if prolate { 1.0 / t.tanh() } else { t.tanh() } };
    let h = STEP;
    for i in start..steps {
        let t = STEP * i as f64;
        let f = |t: f64, y: f64, p: f64| -> f64 { potential(t) * y - damping(t) * p };
        let (k1y, k1p) = (p, f(t, y, p));
        let (k2y, k2p) = (p + 0.5 * h * k1p, f(t + 0.5 * h, y + 0.5 * h * k1y, p + 0.5 * h * k1p));
        let (k3y, k3p) = (p + 0.5 * h * k2p, f(t + 0.5 * h, y + 0.5 * h * k2y, p + 0.5 * h * k2p));
        let (k4y, k4p) = (p + h * k3p, f(t + h, y + h * k3y, p + h * k3p));
        y += h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
        p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        table[i + 1] = y;
    }
    let value = Rc::new(table);
    RADIAL.with(|c| {
        let mut cache = c.borrow_mut();
        if cache.len() >= 300 {
            cache.remove(0);
        }
        cache.push((key, Rc::clone(&value)));
    });
    Some(value)
}

/// The regular radial function at `x` (`ξ ≥ 1` prolate, `ζ ≥ 0` oblate).
pub(super) fn radial_value(
    prolate: bool,
    m: usize,
    n: usize,
    gamma: f64,
    x: f64,
) -> f64 {
    let Some(table) = radial_table(prolate, m, n, gamma) else {
        return f64::NAN;
    };
    let t = if prolate {
        if x < 1.0 {
            return f64::NAN;
        }
        x.acosh()
    } else {
        if x < 0.0 {
            return f64::NAN;
        }
        x.asinh()
    };
    let position = t / STEP;
    if position >= (table.len() - 3) as f64 {
        return f64::NAN;
    }
    let i = (position.floor() as usize).clamp(1, table.len() - 3);
    let u = position - i as f64;
    let (p0, p1, p2, p3) = (table[i - 1], table[i], table[i + 1], table[i + 2]);
    let a = -u * (u - 1.0) * (u - 2.0) / 6.0;
    let b = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    let c = -(u + 1.0) * u * (u - 2.0) / 2.0;
    let d = (u + 1.0) * u * (u - 1.0) / 6.0;
    a * p0 + b * p1 + c * p2 + d * p3
}

/// `spheroidal_radial(kind, m, n, γ, x)`, `kind = 1` prolate, `-1` oblate.
pub(super) fn spheroidal_radial(args: &[f64]) -> f64 {
    let (Some(&kind), Some(&m), Some(&n), Some(&gamma), Some(&x)) = (args.first(), args.get(1), args.get(2), args.get(3), args.get(4)) else {
        return f64::NAN;
    };
    if m < 0.0 || n < m || m.fract() != 0.0 || n.fract() != 0.0 {
        return f64::NAN;
    }
    radial_value(kind > 0.0, m as usize, n as usize, gamma, x)
}

/// Registers the operator.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("spheroidal_radial", Arity::Fixed(5)).eval(spheroidal_radial))?;
    Ok(())
}

/// The recognised structure of a spheroidal Helmholtz equation.
struct Structure {
    prolate: bool,
    /// The angular parameter `γ`.
    gamma: f64,
    radial: usize,
    eta: usize,
    azimuth: Option<usize>,
}

fn evaluate(
    cx: &Cx<'_>,
    node: NodeId,
    bindings: &[(crate::graph::SymbolId, f64)],
) -> Option<f64> {
    let mut env = Env::numeric(0.0);
    for &(s, v) in bindings {
        env.bind(s, v);
    }
    cx.graph.eval(node, &env).filter(|v| v.is_finite())
}

/// Recognises the equation from the numeric values of its coefficients.
fn recognise(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Structure> {
    let n = p.dimension();
    if !(2..=3).contains(&n) || p.nonlinear || !p.homogeneous(cx.graph) {
        return None;
    }
    let name = |cx: &Cx<'_>, j: usize| cx.graph.symbol_of(p.vars[j]).map(|s| cx.graph.interner().symbol_name(s).to_owned());
    let find = |cx: &Cx<'_>, wanted: &[&str]| (0..n).find(|&j| name(cx, j).is_some_and(|nm| wanted.contains(&nm.as_str())));
    let radial = find(cx, &["xi", "zeta"])?;
    let eta = find(cx, &["eta", "nu"])?;
    let azimuth = find(cx, &["phi", "varphi", "ph"]);
    if n == 3 && azimuth.is_none() {
        return None;
    }
    let symbols: Vec<crate::graph::SymbolId> = p.vars.iter().map(|&v| cx.graph.symbol_of(v)).collect::<Option<_>>()?;
    let mut allowed = vec![p.unit(radial, 2), p.unit(radial, 1), p.unit(eta, 2), p.unit(eta, 1), vec![0; n]];
    if let Some(a) = azimuth {
        allowed.push(p.unit(a, 2));
    }
    if !p.only(cx.graph, &allowed) {
        return None;
    }
    let coefficients: Vec<(Vec<u32>, NodeId)> = allowed.iter().map(|i| (i.clone(), p.coefficient(cx.graph, i))).collect();
    for prolate in [true, false] {
        let mut scale: Option<f64> = None;
        let mut gamma: Option<f64> = None;
        let mut ok = true;
        'samples: for &(xi, eta_v) in &[(1.7, 0.3), (2.4, -0.6), (3.1, 0.8), (1.3, 0.15)] {
            let mut bindings = vec![(symbols[radial], xi), (symbols[eta], eta_v)];
            if let Some(a) = azimuth {
                bindings.push((symbols[a], 0.4));
            }
            let (xx, ee) = (xi * xi, eta_v * eta_v);
            let radial_factor = if prolate { xx - 1.0 } else { xx + 1.0 };
            let mix = if prolate { xx - ee } else { xx + ee };
            for (index, node) in &coefficients {
                let Some(value) = evaluate(cx, *node, &bindings) else {
                    ok = false;
                    break 'samples;
                };
                let expected = if *index == p.unit(radial, 2) {
                    radial_factor
                } else if *index == p.unit(radial, 1) {
                    2.0 * xi
                } else if *index == p.unit(eta, 2) {
                    1.0 - ee
                } else if *index == p.unit(eta, 1) {
                    -2.0 * eta_v
                } else if azimuth.is_some_and(|a| *index == p.unit(a, 2)) {
                    mix / (radial_factor * (1.0 - ee))
                } else {
                    // The zero-order term: g × mix.
                    let a = scale.unwrap_or(1.0);
                    let g = value / (a * mix);
                    match gamma {
                        | None => gamma = Some(g),
                        | Some(old) if (old - g).abs() > 1e-9 * (1.0 + g.abs()) => {
                            ok = false;
                            break 'samples;
                        },
                        | Some(_) => {},
                    }
                    continue;
                };
                if *index == p.unit(radial, 2) && scale.is_none() {
                    scale = Some(value / expected);
                }
                let a = scale?;
                if (value - a * expected).abs() > 1e-9 * (1.0 + value.abs()) {
                    ok = false;
                    break 'samples;
                }
            }
        }
        if ok && let (Some(_), Some(g)) = (scale, gamma) {
            return Some(Structure { prolate, gamma: if prolate { g } else { -g }, radial, eta, azimuth });
        }
    }
    None
}

/// The Dirichlet problem inside a spheroid by the finite expansion.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let structure = recognise(cx, p)?;
    let [boundary] = conditions.0.as_slice() else {
        return None;
    };
    if boundary.on != structure.radial || boundary.robin.is_some() || boundary.derivative.iter().any(|&d| d != 0) {
        return None;
    }
    let x0 = super::modes::exact_value(cx, boundary.point)?;
    if (structure.prolate && x0 <= 1.0) || (!structure.prolate && x0 <= 0.0) || x0 > 20.0 {
        return None;
    }
    let (eta_node, radial_node) = (p.vars[structure.eta], p.vars[structure.radial]);
    let eta_symbol = cx.graph.symbol_of(eta_node)?;
    let phi_symbol = structure.azimuth.and_then(|a| cx.graph.symbol_of(p.vars[a]));
    let data = boundary.value;
    let free = cx.graph.free_symbols(cx.graph.find(data)).to_vec();
    if free.iter().any(|&s| s != eta_symbol && Some(s) != phi_symbol) {
        return None;
    }
    let g = |cx: &Cx<'_>, eta_v: f64, phi_v: f64| -> Option<f64> {
        let mut bindings = vec![(eta_symbol, eta_v)];
        if let Some(s) = phi_symbol {
            bindings.push((s, phi_v));
        }
        evaluate(cx, data, &bindings)
    };
    // Quadrature nodes: Gauss–Legendre in η, trapezoid in φ.
    let eta_nodes = super::curvilinear::gauss_legendre(48);
    let phi_count = if phi_symbol.is_some() { 64 } else { 1 };
    let phi_values: Vec<f64> = (0..phi_count).map(|k| 2.0 * PI * k as f64 / phi_count as f64).collect();
    let mut grid = vec![vec![0.0; phi_count]; eta_nodes.len()];
    let mut scale = 1e-12_f64;
    for (i, &(e, _)) in eta_nodes.iter().enumerate() {
        for (j, &ph) in phi_values.iter().enumerate() {
            grid[i][j] = g(cx, e, ph)?;
            scale = scale.max(grid[i][j].abs());
        }
    }
    let gamma = structure.gamma;
    let m_top = if phi_symbol.is_some() { M_MAX } else { 0 };
    // (m, n, cos coefficient, sin coefficient)
    let mut coefficients: Vec<(usize, usize, f64, f64)> = Vec::new();
    for m in 0..=m_top {
        let pair = spheroidal_pair(m, gamma);
        let vectors = &pair.1;
        for j in 0..MODES {
            let n = m + j;
            let Some(d) = vectors.get(j) else {
                continue;
            };
            let sign = if d[j] < 0.0 { -1.0 } else { 1.0 };
            let mut cos_sum = 0.0;
            let mut sin_sum = 0.0;
            for (i, &(e, w)) in eta_nodes.iter().enumerate() {
                let s = super::coulomb::spheroidal_mode_value(m, gamma, e, d, sign);
                for (jj, &ph) in phi_values.iter().enumerate() {
                    cos_sum += w * s * grid[i][jj] * (m as f64 * ph).cos();
                    sin_sum += w * s * grid[i][jj] * (m as f64 * ph).sin();
                }
            }
            let phi_norm = if m == 0 { phi_count as f64 } else { phi_count as f64 / 2.0 };
            let (a, b) = (cos_sum / phi_norm, if m == 0 { 0.0 } else { sin_sum / phi_norm });
            if a.abs() > 1e-10 * scale || b.abs() > 1e-10 * scale {
                coefficients.push((m, n, a, b));
            }
        }
    }
    // The expansion must reproduce the data.
    for &(e, ph) in &[(0.31, 0.7), (-0.62, 2.1), (0.83, 4.4), (-0.17, 5.5)] {
        let mut total = 0.0;
        for &(m, n, a, b) in &coefficients {
            let pair = spheroidal_pair(m, gamma);
            let j = n - m;
            let d = pair.1.get(j)?;
            let sign = if d[j] < 0.0 { -1.0 } else { 1.0 };
            let s = super::coulomb::spheroidal_mode_value(m, gamma, e, d, sign);
            let mp = m as f64 * ph;
            total += s * (a * mp.cos() + b * mp.sin());
        }
        if (total - g(cx, e, ph)?).abs() > 1e-4 * scale {
            return None;
        }
    }
    // Assemble the symbolic expansion.
    let kind = cx.graph.int(if structure.prolate { 1 } else { -1 });
    let gamma_node = cx.graph.num(Number::from(gamma));
    let x0_node = boundary.point;
    let mut terms = Vec::new();
    for &(m, n, a, b) in &coefficients {
        let (m_node, n_node) = (cx.graph.int(i64::try_from(m).ok()?), cx.graph.int(i64::try_from(n).ok()?));
        let angular = call(cx, "spheroidal_mode", &[m_node, n_node, gamma_node, p.vars[structure.eta]])?;
        let at = call(cx, "spheroidal_radial", &[kind, m_node, n_node, gamma_node, radial_node])?;
        let at0 = call(cx, "spheroidal_radial", &[kind, m_node, n_node, gamma_node, x0_node])?;
        let ratio = div(cx, at, at0);
        let mut azimuthal = Vec::new();
        if let Some(index) = structure.azimuth {
            let phi = p.vars[index];
            let m_phi = mul(cx.graph, &[m_node, phi]);
            if a.abs() > 0.0 {
                let c = call(cx, "cos", &[m_phi])?;
                let a_node = cx.graph.num(Number::from(a));
                azimuthal.push(mul(cx.graph, &[a_node, c]));
            }
            if b.abs() > 0.0 {
                let s = call(cx, "sin", &[m_phi])?;
                let b_node = cx.graph.num(Number::from(b));
                azimuthal.push(mul(cx.graph, &[b_node, s]));
            }
        } else {
            azimuthal.push(cx.graph.num(Number::from(a)));
        }
        let trig = add(cx.graph, &azimuthal);
        terms.push(mul(cx.graph, &[trig, angular, ratio]));
    }
    if terms.is_empty() {
        return None;
    }
    let u = add(cx.graph, &terms);
    Some(cx.simplify(u))
}
