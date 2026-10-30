//! Numeric modal expansions for heat and wave equations on disks and balls.
//!
//! When the initial data is not a finite combination of the modes (any
//! smooth `f(r, θ[, φ])`), the exact engine leaves unevaluated integrals or
//! gives up. This module computes the expansion coefficients at solve time:
//!
//! * the radial integrals by Gauss–Legendre quadrature, the angular ones by
//!   the trapezoidal rule (disk, `φ`) and Gauss–Legendre in `cos θ` (ball),
//!   with the data evaluated numerically from the compiled expression;
//! * the eigenvalues from the root operators (`bessel_zero`, `bessel_root`):
//!   Dirichlet zeros of `J_ν`, Robin/Neumann roots of
//!   `z J_ν'(z) + (hR) J_ν(z)` (ball: `ν = l + 1/2`, `h R - 1/2`);
//! * the answer is the truncated sum
//!   `Σ c_{n m} T(t) J_n(κ r) {cos, sin}(n θ)` (disk: `n ≤ 10`, 16 radial
//!   modes) or `Σ c_{l m k} T(t) j_l(κ r) P_l^m(cos θ) {cos, sin}(m φ)`
//!   (ball: `l ≤ 8`, 12 radial modes), with numeric coefficients and the
//!   time factor `e^{-(Dκ² + e)t}` (heat) or `cos ω t`, `sin(ω t)/ω`
//!   (wave, `ω² = Dκ² + e`).
//!
//! The truncation is the only error besides the quadrature; it is small for
//! smooth data and times `t > 0` (heat), larger for waves and near `t = 0`.

use super::Conditions;
use super::Problem;
use super::curvilinear::associated_values;
use super::curvilinear::gauss_legendre;
use super::modes::associated_mode;
use super::modes::exact_value;
use super::util::call;
use super::util::div;
use super::util::is_zero_number;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::kernels::special::bessel_j;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use std::f64::consts::PI;

/// Which variables are the radius and the angles.
#[derive(Copy, Clone, Debug)]
pub(super) struct Geometry {
    pub r: usize,
    pub theta: Option<usize>,
    pub phi: Option<usize>,
    pub sphere: bool,
}

const DISK_ANGLES: usize = 10;
const DISK_ROOTS: usize = 16;
const BALL_DEGREE: usize = 8;
const BALL_ROOTS: usize = 12;
const RADIAL_NODES: usize = 64;

/// A mode: eigenvalue `κ²` (as `κ`), the value of the radial function at the
/// quadrature radii, and its node in the answer.
struct RadialMode {
    kappa: f64,
    values: Vec<f64>,
    node: NodeId,
}

/// The kind of the condition at the boundary.
#[derive(Copy, Clone, Debug)]
enum Boundary {
    Dirichlet,
    /// `u_r + h u = 0` (`h = 0` for Neumann).
    Robin(f64),
}

fn evaluate(
    cx: &Cx<'_>,
    node: NodeId,
    bindings: &[(crate::graph::SymbolId, f64)],
) -> Option<f64> {
    let mut env = Env::numeric(0.0);
    for &(symbol, value) in bindings {
        env.bind(symbol, value);
    }
    cx.graph.eval(node, &env).filter(|v| v.is_finite())
}

/// The numeric expansion, if the problem is a homogeneous heat or undamped
/// wave problem on a disk or ball with a homogeneous boundary condition.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    geometry: Geometry,
    time_index: usize,
    (b, a0): (NodeId, NodeId),
) -> Option<NodeId> {
    if p.nonlinear || !p.homogeneous(cx.graph) || (geometry.phi.is_some() && !geometry.sphere) {
        return None;
    }
    let n = p.dimension();
    let (c1, c2) = (
        p.coefficient(cx.graph, &p.unit(time_index, 1)),
        p.coefficient(cx.graph, &p.unit(time_index, 2)),
    );
    let order = if is_zero_number(cx.graph, c2) {
        1
    } else {
        2
    };
    if order == 2 && !is_zero_number(cx.graph, c1) {
        return None;
    }
    let c_time = if order == 1 { c1 } else { c2 };
    let d = {
        let q = div(cx, b, c_time);
        let q = neg(cx.graph, q);
        exact_value(cx, q)?
    };
    let e = {
        let q = div(cx, a0, c_time);
        exact_value(cx, q)?
    };
    if d <= 0.0 {
        return None;
    }
    // The boundary condition and the data.
    let zero_index = vec![0; n];
    let mut boundary = None;
    for c in &conditions.0 {
        if c.on == geometry.r {
            if boundary.is_some() || !cx.is_zero(c.value) {
                return None;
            }
            let radius = exact_value(cx, c.point).filter(|v| *v > 0.0)?;
            let kind = if c.derivative == zero_index {
                Boundary::Dirichlet
            } else if c.derivative == p.unit(geometry.r, 1) {
                match c.robin {
                    | None => Boundary::Robin(0.0),
                    | Some(h) => Boundary::Robin(exact_value(cx, h)?),
                }
            } else {
                return None;
            };
            boundary = Some((radius, kind));
        } else if c.on != time_index {
            return None;
        }
    }
    let (radius, kind) = boundary?;
    let displacement = conditions.find(cx.graph, time_index, &zero_index, None)?;
    let velocity = if order == 2 {
        conditions.find(cx.graph, time_index, &p.unit(time_index, 1), None)
    } else {
        None
    };
    let expected = 1 + usize::from(velocity.is_some());
    if conditions.0.iter().filter(|c| c.on == time_index).count() != expected {
        return None;
    }
    let (r_var, t) = (p.vars[geometry.r], p.vars[time_index]);
    let r_symbol = cx.graph.symbol_of(r_var)?;
    let theta_symbol = geometry.theta.and_then(|j| cx.graph.symbol_of(p.vars[j]));
    let phi_symbol = geometry.phi.and_then(|j| cx.graph.symbol_of(p.vars[j]));
    let allowed = [Some(r_symbol), theta_symbol, phi_symbol];
    let data: Vec<(NodeId, bool)> = [
        (displacement.value, true),
        (velocity.map_or_else(|| cx.graph.int(0), |v| v.value), false),
    ]
    .into_iter()
    .filter(|(v, _)| !cx.is_zero(*v))
    .collect();
    for (value, _) in &data {
        for symbol in cx.graph.free_symbols(cx.graph.find(*value)).to_vec() {
            if !allowed.contains(&Some(symbol)) {
                return None;
            }
        }
    }
    if data.is_empty() {
        return None;
    }
    // Quadrature in the radius.
    let radial: Vec<(f64, f64)> = gauss_legendre(RADIAL_NODES)
        .into_iter()
        .map(|(x, w)| (0.5 * radius * (1.0 + x), 0.5 * radius * w))
        .collect();
    // Angular grids: (cos θ or θ, weight, θ), φ values.
    let (angles, azimuths): (Vec<(f64, f64, f64)>, Vec<f64>) = if geometry.sphere {
        let theta_nodes = if geometry.theta.is_some() {
            24
        } else {
            1
        };
        let a: Vec<(f64, f64, f64)> = if geometry.theta.is_some() {
            gauss_legendre(theta_nodes)
                .into_iter()
                .map(|(x, w)| (x, w, x.acos()))
                .collect()
        } else {
            vec![(0.0, 2.0, 0.0)]
        };
        let m = if geometry.phi.is_some() {
            32
        } else {
            1
        };
        (
            a,
            (0..m)
                .map(|k| 2.0 * PI * f64::from(k) / f64::from(m))
                .collect(),
        )
    } else {
        let m = if geometry.theta.is_some() {
            48
        } else {
            1
        };
        (
            (0..m)
                .map(|k| {
                    (
                        0.0,
                        2.0 * PI / f64::from(m),
                        2.0 * PI * f64::from(k) / f64::from(m),
                    )
                })
                .collect(),
            vec![0.0],
        )
    };
    // Tables of the data on the grid: [datum][radius][angle][azimuth].
    let mut tables: Vec<Vec<Vec<Vec<f64>>>> = Vec::new();
    for (value, _) in &data {
        let mut by_radius = Vec::new();
        for &(rv, _) in &radial {
            let mut by_angle = Vec::new();
            for &(_, _, theta) in &angles {
                let mut row = Vec::new();
                for &phi in &azimuths {
                    let mut bindings = vec![(r_symbol, rv)];
                    if let Some(s) = theta_symbol {
                        bindings.push((s, theta));
                    }
                    if let Some(s) = phi_symbol {
                        bindings.push((s, phi));
                    }
                    row.push(evaluate(cx, *value, &bindings)?);
                }
                by_angle.push(row);
            }
            by_radius.push(by_angle);
        }
        tables.push(by_radius);
    }
    let mut terms = Vec::new();
    let time_factor = |cx: &mut Cx<'_>, kappa: f64, is_displacement: bool| -> Option<NodeId> {
        let omega2 = d * kappa * kappa + e;
        let t_scaled = |cx: &mut Cx<'_>, rate: f64| {
            let c = cx.graph.num(Number::from(rate));
            mul(cx.graph, &[c, t])
        };
        if order == 1 {
            let arg = t_scaled(cx, -omega2);
            return call(cx, "exp", &[arg]);
        }
        if omega2 <= 0.0 {
            if omega2 < -1e-12 {
                return None;
            }
            return Some(if is_displacement {
                cx.graph.int(1)
            } else {
                t
            });
        }
        let omega = omega2.sqrt();
        let arg = t_scaled(cx, omega);
        if is_displacement {
            call(cx, "cos", &[arg])
        } else {
            let s = call(cx, "sin", &[arg])?;
            let inverse = cx.graph.num(Number::from(1.0 / omega));
            Some(mul(cx.graph, &[inverse, s]))
        }
    };
    // Radial modes of a Bessel order `nu` (the disk order, or `l + 1/2`).
    let radial_modes =
        |cx: &mut Cx<'_>, nu: f64, ball: bool, constant: bool| -> Option<Vec<RadialMode>> {
            let count = if ball {
                BALL_ROOTS
            } else {
                DISK_ROOTS
            };
            let mut modes = Vec::new();
            let mut roots: Vec<f64> = Vec::new();
            if constant {
                roots.push(0.0);
            }
            for m in 1..=count {
                let z = match kind {
                    | Boundary::Dirichlet => super::numeric::bessel_zero(&[nu, m as f64]),
                    | Boundary::Robin(h) => {
                        let shifted = if ball {
                            h * radius - 0.5
                        } else {
                            h * radius
                        };
                        if shifted + nu <= 0.0 {
                            return None;
                        }
                        super::numeric::bessel_root(&[nu, shifted, m as f64])
                    },
                };
                if !z.is_finite() {
                    return None;
                }
                roots.push(z / radius);
            }
            for kappa in roots {
                if kappa.abs() < 1e-14 {
                    let values = vec![1.0; radial.len()];
                    let node = cx.graph.int(1);
                    modes.push(RadialMode { kappa, values, node });
                    continue;
                }
                let scale = if ball {
                    (PI / (2.0 * kappa)).sqrt()
                } else {
                    1.0
                };
                let values: Vec<f64> = radial
                    .iter()
                    .map(|&(rv, _)| {
                        if ball {
                            scale * bessel_j(nu, kappa * rv) / rv.sqrt()
                        } else {
                            bessel_j(nu, kappa * rv)
                        }
                    })
                    .collect();
                let kr = {
                    let k_node = cx.graph.num(Number::from(kappa));
                    mul(cx.graph, &[k_node, r_var])
                };
                let order_node = cx
                    .graph
                    .num(Number::fraction((2.0 * nu).round() as i64, 2)?);
                let bessel = call(cx, "besselj", &[order_node, kr])?;
                let node = if ball {
                    let minus_half = cx.graph.num(Number::fraction(-1, 2)?);
                    let root = pow(cx.graph, r_var, minus_half);
                    let c = cx.graph.num(Number::from(scale));
                    mul(cx.graph, &[c, root, bessel])
                } else {
                    bessel
                };
                modes.push(RadialMode { kappa, values, node });
            }
            Some(modes)
        };
    let radial_weight = |i: usize| -> f64 {
        let (rv, w) = radial.get(i).copied().unwrap_or((0.0, 0.0));
        if geometry.sphere {
            w * rv * rv
        } else {
            w * rv
        }
    };
    if geometry.sphere {
        let lmax = if geometry.theta.is_some() {
            BALL_DEGREE
        } else {
            0
        };
        for l in 0..=lmax {
            let nu = l as f64 + 0.5;
            let constant = l == 0 && matches!(kind, Boundary::Robin(h) if h.abs() < 1e-14);
            let modes = radial_modes(cx, nu, true, constant)?;
            let mmax = if geometry.phi.is_some() {
                l
            } else {
                0
            };
            for m in 0..=mmax {
                // Angular functions: P_l^m(cos θ) cos(mφ) and sin(mφ).
                let p_values: Vec<f64> = angles
                    .iter()
                    .map(|&(x, _, _)| associated_values(x, m, l).get(l).copied().unwrap_or(0.0))
                    .collect();
                let theta_node = geometry.theta.map(|j| p.vars[j]);
                let angular_theta = match theta_node {
                    | Some(th) => {
                        associated_mode(cx, th, i64::try_from(l).ok()?, i64::try_from(m).ok()?)?
                    },
                    | None => cx.graph.int(1),
                };
                let mut norm_theta = 2.0 / (2 * l + 1) as f64;
                for k in (l - m + 1)..=(l + m) {
                    norm_theta *= k as f64;
                }
                let phi_nodes: Vec<(NodeId, bool)> = match geometry.phi.map(|j| p.vars[j]) {
                    | Some(ph) => {
                        let m_node = cx.graph.int(i64::try_from(m).ok()?);
                        let m_phi = mul(cx.graph, &[m_node, ph]);
                        let mut v = vec![(call(cx, "cos", &[m_phi])?, true)];
                        if m > 0 {
                            v.push((call(cx, "sin", &[m_phi])?, false));
                        }
                        v
                    },
                    | None => vec![(cx.graph.int(1), true)],
                };
                for (phi_node, is_cos) in phi_nodes {
                    let norm_phi = if geometry.phi.is_none() || m == 0 {
                        2.0 * PI
                    } else {
                        PI
                    };
                    let norm_ang = norm_theta * norm_phi;
                    for mode in &modes {
                        let norm_r: f64 = mode
                            .values
                            .iter()
                            .enumerate()
                            .map(|(i, v)| radial_weight(i) * v * v)
                            .sum();
                        for (table, &(_, is_displacement)) in tables.iter().zip(&data) {
                            let mut total = 0.0;
                            for (i, by_angle) in table.iter().enumerate() {
                                let mut angular = 0.0;
                                for (j, row) in by_angle.iter().enumerate() {
                                    let (_, wx, _) = angles[j];
                                    for (k, &fv) in row.iter().enumerate() {
                                        let trig = if geometry.phi.is_none() {
                                            1.0
                                        } else if is_cos {
                                            (m as f64 * azimuths[k]).cos()
                                        } else {
                                            (m as f64 * azimuths[k]).sin()
                                        };
                                        angular += wx
                                            * fv
                                            * p_values[j]
                                            * trig
                                            * (2.0 * PI / azimuths.len() as f64);
                                    }
                                }
                                total += radial_weight(i) * mode.values[i] * angular;
                            }
                            let coefficient = total / (norm_r * norm_ang);
                            if coefficient.abs() < 1e-14 {
                                continue;
                            }
                            let factor = time_factor(cx, mode.kappa, is_displacement)?;
                            let c = cx.graph.num(Number::from(coefficient));
                            terms.push(mul(
                                cx.graph,
                                &[c, factor, mode.node, angular_theta, phi_node],
                            ));
                        }
                    }
                }
            }
        }
    } else {
        let nmax = if geometry.theta.is_some() {
            DISK_ANGLES
        } else {
            0
        };
        for order_n in 0..=nmax {
            let nu = order_n as f64;
            let constant = order_n == 0 && matches!(kind, Boundary::Robin(h) if h.abs() < 1e-14);
            let modes = radial_modes(cx, nu, false, constant)?;
            let theta_node = geometry.theta.map(|j| p.vars[j]);
            let mut trig_nodes = Vec::new();
            match theta_node {
                | Some(th) => {
                    let n_node = cx.graph.int(i64::try_from(order_n).ok()?);
                    let n_theta = mul(cx.graph, &[n_node, th]);
                    trig_nodes.push((call(cx, "cos", &[n_theta])?, true));
                    if order_n > 0 {
                        trig_nodes.push((call(cx, "sin", &[n_theta])?, false));
                    }
                },
                | None => trig_nodes.push((cx.graph.int(1), true)),
            }
            for (trig_node, is_cos) in trig_nodes {
                let norm_ang = if geometry.theta.is_none() || order_n == 0 {
                    2.0 * PI
                } else {
                    PI
                };
                for mode in &modes {
                    let norm_r: f64 = mode
                        .values
                        .iter()
                        .enumerate()
                        .map(|(i, v)| radial_weight(i) * v * v)
                        .sum();
                    for (table, &(_, is_displacement)) in tables.iter().zip(&data) {
                        let mut total = 0.0;
                        for (i, by_angle) in table.iter().enumerate() {
                            let mut angular = 0.0;
                            for (j, row) in by_angle.iter().enumerate() {
                                let (_, w, theta) = angles[j];
                                let trig = if geometry.theta.is_none() {
                                    1.0
                                } else if is_cos {
                                    (order_n as f64 * theta).cos()
                                } else {
                                    (order_n as f64 * theta).sin()
                                };
                                angular += w * row.first().copied().unwrap_or(0.0) * trig;
                            }
                            total += radial_weight(i) * mode.values[i] * angular;
                        }
                        let coefficient = total / (norm_r * norm_ang);
                        if coefficient.abs() < 1e-14 {
                            continue;
                        }
                        let factor = time_factor(cx, mode.kappa, is_displacement)?;
                        let c = cx.graph.num(Number::from(coefficient));
                        terms.push(mul(cx.graph, &[c, factor, mode.node, trig_node]));
                    }
                }
            }
        }
    }
    if terms.is_empty() {
        return None;
    }
    Some(add(cx.graph, &terms))
}
