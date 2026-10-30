//! Disks, cylinders and balls: Laplacian in polar, cylindrical and
//! spherical coordinates.
//!
//! The equation is recognised by the structure of its principal part,
//! `b (u_rr + (d-1)/r u_r + u_θθ/r² + u_zz [+ cot θ/r² u_θ])` plus
//! `a₂ u_tt + a₁ u_t + a₀ u + f`, in the variables
//!
//! * `(r)`, `(r, θ)` — radial disk (`d = 2`) or ball (`d = 3`),
//! * `(r, θ)` — a disk with angular dependence, or a ball with axial
//!   symmetry (told apart by the `u_r` and `cot θ u_θ` terms),
//! * `(r, z)`, `(r, θ, z)` — cylinders,
//!
//! and solved by the eigenfunction engine with Bessel functions
//! `J_n(j_{n,m} r/R)` (disk, cylinder), spherical Bessel functions
//! `j_l(κ r)` with Legendre functions `P_l(cos θ)` (ball). The condition at
//! `r = R` may be Dirichlet, Neumann or Robin. Harmonic functions with
//! Dirichlet data use the closed Fourier / Legendre series instead.

use super::Condition;
use super::Conditions;
use super::Method;
use super::Problem;
use super::classical;
use super::spectral::Axis;
use super::spectral::End;
use super::spectral::Engine;
use super::spectral::Exterior;
use super::spectral::Item;
use super::spectral::Kind;
use super::spectral::Role;
use super::spectral::Shape;
use super::spectral::Special;
use super::spectral::Symbol;
use super::spectral::Time;
use super::spectral::boundary_items;
use super::spectral::evolution_time;
use super::spectral::finish;
use super::spectral::interval;
use super::spectral::time_items;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::is_zero_number;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

#[derive(Copy, Clone, Debug)]
struct Layout {
    r: usize,
    /// The polar angle of a ball, or the angle of a disk.
    theta: Option<usize>,
    /// The azimuth of a ball without axial symmetry.
    phi: Option<usize>,
    z: Option<usize>,
    sphere: bool,
}

fn name_of(
    cx: &Cx<'_>,
    node: NodeId,
) -> Option<String> {
    cx.graph
        .symbol_of(node)
        .map(|s| cx.graph.interner().symbol_name(s).to_owned())
}

/// Names of polar (or sole) angles and of azimuths.
const POLAR_NAMES: [&str; 3] = ["theta", "th", "vartheta"];
const AZIMUTH_NAMES: [&str; 3] = ["phi", "varphi", "ph"];

/// The variable `x` for which the coefficients have the radial structure
/// `b (u_xx + d/x u_x)` (`d` = 1 or 2), when no variable is called `r`.
fn radial_by_structure(
    cx: &mut Cx<'_>,
    p: &Problem,
    candidates: &[usize],
) -> Option<usize> {
    for &j in candidates {
        let (b, c1) = (
            p.coefficient(cx.graph, &p.unit(j, 2)),
            p.coefficient(cx.graph, &p.unit(j, 1)),
        );
        if cx.is_zero(b) || cx.is_zero(c1) || !p.constant(cx.graph, b) {
            continue;
        }
        let x = p.vars[j];
        for d in [1_i64, 2] {
            let k = cx.graph.int(d);
            let inverse = powi(cx.graph, x, -1);
            let expected = mul(cx.graph, &[b, k, inverse]);
            let gap = sub(cx.graph, c1, expected);
            if cx.is_zero(gap) {
                return Some(j);
            }
        }
    }
    None
}

/// The layout of the variables, from their names (and, for a lone radius
/// or radius with an axis, from the structure of the equation).
fn layouts(
    cx: &mut Cx<'_>,
    p: &Problem,
    space: &[usize],
) -> Vec<Layout> {
    let names: Vec<(usize, String)> = space
        .iter()
        .filter_map(|&j| name_of(cx, p.vars[j]).map(|n| (j, n)))
        .collect();
    let find = |wanted: &[&str]| {
        names
            .iter()
            .find(|(_, n)| wanted.contains(&n.as_str()))
            .map(|(j, _)| *j)
    };
    let z = find(&["z"]);
    let r = find(&["r", "rho"]).or_else(|| {
        let candidates: Vec<usize> = space.iter().copied().filter(|&j| Some(j) != z).collect();
        if candidates.len() == 1 {
            radial_by_structure(cx, p, &candidates)
        } else {
            None
        }
    });
    let Some(r) = r else {
        return Vec::new();
    };
    let others: Vec<usize> = space
        .iter()
        .copied()
        .filter(|&j| j != r && Some(j) != z)
        .collect();
    let is = |j: usize, list: &[&str]| {
        names
            .iter()
            .any(|(k, n)| *k == j && list.contains(&n.as_str()))
    };
    let (theta, phi) = match others.as_slice() {
        | [] => (None, None),
        | [j] if is(*j, &POLAR_NAMES) || is(*j, &AZIMUTH_NAMES) => (Some(*j), None),
        | [a, b] => {
            let polar = [*a, *b].into_iter().find(|&j| is(j, &POLAR_NAMES));
            let azimuth = [*a, *b].into_iter().find(|&j| is(j, &AZIMUTH_NAMES));
            match (polar, azimuth) {
                | (Some(t), Some(f)) => (Some(t), Some(f)),
                | _ => return Vec::new(),
            }
        },
        | _ => return Vec::new(),
    };
    let mut out = Vec::new();
    if phi.is_none() {
        out.push(Layout {
            r,
            theta,
            phi,
            z,
            sphere: false,
        });
    }
    if z.is_none() {
        out.push(Layout {
            r,
            theta,
            phi,
            z,
            sphere: true,
        });
    }
    out
}

/// `(b, a₀)` if the spatial part of the equation is `b Δ + a₀` in the
/// given coordinates.
fn laplacian(
    cx: &mut Cx<'_>,
    p: &Problem,
    layout: Layout,
    time: Option<usize>,
) -> Option<(NodeId, NodeId)> {
    let n = p.dimension();
    let b = p.coefficient(cx.graph, &p.unit(layout.r, 2));
    if cx.is_zero(b) || !p.constant(cx.graph, b) {
        return None;
    }
    let r = p.vars[layout.r];
    let mut expected: Vec<(Vec<u32>, NodeId)> = vec![(p.unit(layout.r, 2), cx.graph.int(1))];
    let dimension = if layout.sphere { 2 } else { 1 };
    let ratio = {
        let k = cx.graph.int(dimension);
        let inverse = powi(cx.graph, r, -1);
        mul(cx.graph, &[k, inverse])
    };
    expected.push((p.unit(layout.r, 1), ratio));
    let r2 = powi(cx.graph, r, -2);
    if let Some(j) = layout.theta {
        expected.push((p.unit(j, 2), r2));
        if layout.sphere {
            let theta = p.vars[j];
            let cot = {
                let cos = super::util::call(cx, "cos", &[theta])?;
                let sin = super::util::call(cx, "sin", &[theta])?;
                div(cx, cos, sin)
            };
            let c = mul(cx.graph, &[cot, r2]);
            expected.push((p.unit(j, 1), c));
        }
    }
    if let Some(j) = layout.phi {
        let theta = p.vars[layout.theta?];
        let sin = super::util::call(cx, "sin", &[theta])?;
        let inverse = powi(cx.graph, sin, -2);
        let c = mul(cx.graph, &[r2, inverse]);
        expected.push((p.unit(j, 2), c));
    }
    if let Some(j) = layout.z {
        expected.push((p.unit(j, 2), cx.graph.int(1)));
    }
    for (index, c) in &p.linear {
        let is_time = time.is_some_and(|t| index[t] > 0);
        if is_time || *index == vec![0; n] || is_zero_number(cx.graph, *c) {
            continue;
        }
        let (_, e) = expected.iter().find(|(i, _)| i == index)?;
        let be = mul(cx.graph, &[b, *e]);
        let gap = sub(cx.graph, *c, be);
        if !cx.is_zero(gap) {
            return None;
        }
    }
    for (index, e) in &expected {
        let c = p.coefficient(cx.graph, index);
        let be = mul(cx.graph, &[b, *e]);
        let gap = sub(cx.graph, c, be);
        if !cx.is_zero(gap) {
            return None;
        }
    }
    let a0 = p.coefficient(cx.graph, &vec![0; n]);
    if !p.constant(cx.graph, a0) {
        return None;
    }
    Some((b, a0))
}

pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let time_index = evolution_time(cx, p);
    let space: Vec<usize> = (0..p.dimension())
        .filter(|&j| Some(j) != time_index)
        .collect();
    for layout in layouts(cx, p, &space) {
        let Some((b, a0)) = laplacian(cx, p, layout, time_index) else {
            continue;
        };
        let exact = solve_layout(cx, p, conditions, layout, time_index, b, a0);
        // A closed or series answer without unevaluated integrals stands;
        // otherwise try the numeric modal expansion (general data).
        if let Some(solution) = exact
            && !super::util::contains_op(cx.graph, solution, "defint")
        {
            return Some(solution);
        }
        if let (Some(t_index), None) = (time_index, layout.z) {
            let geometry = super::modal::Geometry {
                r: layout.r,
                theta: layout.theta,
                phi: layout.phi,
                sphere: layout.sphere,
            };
            if let Some(solution) =
                super::modal::solve(cx, p, conditions, geometry, t_index, (b, a0))
            {
                return Some(solution);
            }
        }
        if let Some(solution) = exact {
            return Some(solution);
        }
        if let Some(solution) = numeric_sphere(cx, p, conditions, layout, time_index, a0) {
            return Some(solution);
        }
        if let Some(solution) = hankel_whole_space(cx, p, conditions, layout, time_index, (b, a0)) {
            return Some(solution);
        }
    }
    None
}

/// Radial evolution on the whole plane or space by the Hankel transform:
///
/// * plane: `u = ∫_0^∞ k J_0(k r) M(k) [∫_0^∞ ρ f(ρ) J_0(k ρ) dρ] dk`;
/// * space: `u = (2/π) ∫_0^∞ k² j_0(k r) M(k) [∫_0^∞ ρ² f(ρ) j_0(k ρ) dρ] dk`
///   with `j_0(z) = sin z / z`,
///
/// where the multiplier `M` is `e^{-D k² t}` (first order in time) or
/// `cos(c k t)` and `sin(c k t)/(c k)` for displacement and velocity data
/// (second order, `c² = D`).
fn hankel_whole_space(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    (b, a0): (NodeId, NodeId),
) -> Option<NodeId> {
    let t_index = time_index?;
    if layout.theta.is_some()
        || layout.phi.is_some()
        || layout.z.is_some()
        || !cx.is_zero(a0)
        || !p.homogeneous(cx.graph)
    {
        return None;
    }
    for c in &conditions.0 {
        if c.on != t_index && !(c.on == layout.r && is_infinity(cx, c.point) && cx.is_zero(c.value))
        {
            return None;
        }
    }
    let n = p.dimension();
    let c_second = p.coefficient(cx.graph, &p.unit(t_index, 2));
    let order = if is_zero_number(cx.graph, c_second) {
        1
    } else {
        2
    };
    let c_time = p.coefficient(cx.graph, &p.unit(t_index, order));
    let diffusivity = {
        let q = div(cx, b, c_time);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let (r, t) = (p.vars[layout.r], p.vars[t_index]);
    let zero_index = vec![0; n];
    let displacement = conditions.find(cx.graph, t_index, &zero_index, None)?;
    let velocity = if order == 2 {
        conditions.find(cx.graph, t_index, &p.unit(t_index, 1), None)
    } else {
        None
    };
    let expected = 1 + usize::from(velocity.is_some());
    if conditions.0.iter().filter(|c| c.on == t_index).count() != expected {
        return None;
    }
    let (k, _) = super::dummy(cx, p, "k");
    let (rho, _) = super::dummy(cx, p, "rho");
    let zero = cx.graph.int(0);
    let oo = super::util::infinity(cx)?;
    let pi = super::util::pi(cx)?;
    // The kernel of the transform at a point z = k x.
    let kernel = |cx: &mut Cx<'_>, x: NodeId| -> Option<NodeId> {
        let z = mul(cx.graph, &[k, x]);
        if layout.sphere {
            let s = call(cx, "sin", &[z])?;
            Some(div(cx, s, z))
        } else {
            call(cx, "besselj", &[zero, z])
        }
    };
    let sqrt_d = super::util::sqrt(cx, diffusivity)?;
    let ckt = {
        let ck = mul(cx.graph, &[sqrt_d, k]);
        mul(cx.graph, &[ck, t])
    };
    let mut integrands = Vec::new();
    let mut pairs = vec![(displacement.value, true)];
    if let Some(v) = velocity {
        pairs.push((v.value, false));
    }
    for (value, is_displacement) in pairs {
        if cx.is_zero(value) {
            continue;
        }
        let multiplier = if order == 1 {
            let k2 = powi(cx.graph, k, 2);
            let e = mul(cx.graph, &[diffusivity, k2, t]);
            let e = neg(cx.graph, e);
            call(cx, "exp", &[e])?
        } else if is_displacement {
            call(cx, "cos", &[ckt])?
        } else {
            let sn = call(cx, "sin", &[ckt])?;
            let ck = mul(cx.graph, &[sqrt_d, k]);
            div(cx, sn, ck)
        };
        let f_rho = cx.graph.substitute(value, r, rho);
        let weight = powi(cx.graph, rho, if layout.sphere { 2 } else { 1 });
        let k_rho = kernel(cx, rho)?;
        let inner_body = mul(cx.graph, &[weight, f_rho, k_rho]);
        let inner = defint(cx, inner_body, rho, zero, oo)?;
        let k_r = kernel(cx, r)?;
        let measure = powi(cx.graph, k, if layout.sphere { 2 } else { 1 });
        let outer_body = mul(cx.graph, &[measure, k_r, multiplier, inner]);
        integrands.push(defint(cx, outer_body, k, zero, oo)?);
    }
    let total = add(cx.graph, &integrands);
    let u = if layout.sphere {
        let two = cx.graph.int(2);
        let factor = div(cx, two, pi);
        mul(cx.graph, &[factor, total])
    } else {
        total
    };
    Some(cx.simplify(u))
}

/// The Gauss–Legendre nodes and weights on `[-1, 1]`.
pub(super) fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            let nf = n as f64;
            derivative = nf * (x * p1 - p0) / (x * x - 1.0);
            let step = p1 / derivative;
            x -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    out
}

/// `P_l^m(x)` (no Condon–Shortley phase) for all `l ≤ lmax` at fixed `m`.
pub(super) fn associated_values(
    x: f64,
    m: usize,
    lmax: usize,
) -> Vec<f64> {
    let mut out = vec![0.0; lmax + 1];
    let mut pmm = 1.0;
    let s = (1.0 - x * x).max(0.0).sqrt();
    for k in 1..=m {
        pmm *= (2 * k - 1) as f64 * s;
    }
    if let Some(slot) = out.get_mut(m) {
        *slot = pmm;
    }
    if m < lmax {
        let first = x * (2 * m + 1) as f64 * pmm;
        if let Some(slot) = out.get_mut(m + 1) {
            *slot = first;
        }
        let (mut previous, mut current) = (pmm, first);
        for l in (m + 2)..=lmax {
            let next =
                ((2 * l - 1) as f64 * x * current - (l + m - 1) as f64 * previous) / (l - m) as f64;
            if let Some(slot) = out.get_mut(l) {
                *slot = next;
            }
            previous = current;
            current = next;
        }
    }
    out
}

/// Laplace's equation inside or outside a ball with Dirichlet data `g(θ, φ)`
/// of any smoothness: the spherical-harmonic coefficients are computed at
/// solve time by Gauss–Legendre quadrature in `cos θ` and the trapezoidal
/// rule in `φ` (spectrally accurate), up to degree 20; the answer is the
/// truncated series `Σ (r/R)^l P_l^m (a_lm cos mφ + b_lm sin mφ)` (exterior:
/// `(R/r)^{l+1}`) with numeric coefficients.
fn numeric_sphere(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    a0: NodeId,
) -> Option<NodeId> {
    const DEGREE: usize = 20;
    if !layout.sphere || time_index.is_some() || !cx.is_zero(a0) || layout.z.is_some() {
        return None;
    }
    let (theta_j, phi_j) = (layout.theta?, layout.phi?);
    let (theta, phi, r) = (p.vars[theta_j], p.vars[phi_j], p.vars[layout.r]);
    if conditions.0.iter().any(|c| c.on != layout.r) {
        return None;
    }
    let (far, near): (Vec<&Condition>, Vec<&Condition>) =
        conditions.0.iter().partition(|c| is_infinity(cx, c.point));
    let [boundary] = near.as_slice() else {
        return None;
    };
    let exterior = match far.as_slice() {
        | [] => false,
        | [c] if cx.is_zero(c.value) && c.derivative.iter().all(|&d| d == 0) => true,
        | _ => return None,
    };
    if boundary.robin.is_some() || !boundary.derivative.iter().all(|&d| d == 0) {
        return None;
    }
    let radius = super::modes::exact_value(cx, boundary.point).filter(|v| *v > 0.0)?;
    let (theta_symbol, phi_symbol) = (cx.graph.symbol_of(theta)?, cx.graph.symbol_of(phi)?);
    let data = boundary.value;
    for symbol in cx.graph.free_symbols(cx.graph.find(data)).to_vec() {
        if symbol != theta_symbol && symbol != phi_symbol {
            return None;
        }
    }
    let nodes = gauss_legendre(48);
    let azimuthal = 64_usize;
    let mut table = vec![vec![0.0; azimuthal]; nodes.len()];
    for (i, &(x, _)) in nodes.iter().enumerate() {
        let angle = x.acos();
        for k in 0..azimuthal {
            let mut env = crate::graph::Env::numeric(0.0);
            env.bind(theta_symbol, angle);
            env.bind(
                phi_symbol,
                2.0 * std::f64::consts::PI * k as f64 / azimuthal as f64,
            );
            let value = cx.graph.eval(data, &env).filter(|v| v.is_finite())?;
            *table.get_mut(i)?.get_mut(k)? = value;
        }
    }
    let ratio = {
        let radius_node = cx.graph.num(Number::from(radius));
        div(cx, r, radius_node)
    };
    let mut terms = Vec::new();
    for m in 0..=DEGREE {
        // Fourier coefficients in φ of the data at each node.
        let (mut cos_part, mut sin_part) = (vec![0.0; nodes.len()], vec![0.0; nodes.len()]);
        for (i, row) in table.iter().enumerate() {
            let (mut c, mut s_sum) = (0.0, 0.0);
            for (k, &g) in row.iter().enumerate() {
                let angle = m as f64 * 2.0 * std::f64::consts::PI * k as f64 / azimuthal as f64;
                c += g * angle.cos();
                s_sum += g * angle.sin();
            }
            *cos_part.get_mut(i)? = c * 2.0 / azimuthal as f64 / if m == 0 { 2.0 } else { 1.0 };
            *sin_part.get_mut(i)? = s_sum * 2.0 / azimuthal as f64;
        }
        let m_phi = {
            let mm = cx.graph.int(i64::try_from(m).ok()?);
            mul(cx.graph, &[mm, phi])
        };
        let (cos_node, sin_node) = (
            super::util::call(cx, "cos", &[m_phi])?,
            super::util::call(cx, "sin", &[m_phi])?,
        );
        let values: Vec<Vec<f64>> = nodes
            .iter()
            .map(|&(x, _)| associated_values(x, m, DEGREE))
            .collect();
        for l in m..=DEGREE {
            // N_lm = 2/(2l+1) (l+m)!/(l-m)!
            let mut norm = 2.0 / (2 * l + 1) as f64;
            for k in (l - m + 1)..=(l + m) {
                norm *= k as f64;
            }
            let (mut a, mut b) = (0.0, 0.0);
            for (i, &(_, w)) in nodes.iter().enumerate() {
                let pv = values.get(i)?.get(l).copied()?;
                a += w * pv * cos_part.get(i).copied()?;
                b += w * pv * sin_part.get(i).copied()?;
            }
            let (a, b) = (a / norm, b / norm);
            if a.abs() < 1e-13 && b.abs() < 1e-13 {
                continue;
            }
            let angular = super::modes::associated_mode(
                cx,
                theta,
                i64::try_from(l).ok()?,
                i64::try_from(m).ok()?,
            )?;
            let l_signed = i64::try_from(l).ok()?;
            let radial = if exterior {
                powi(cx.graph, ratio, -l_signed - 1)
            } else {
                powi(cx.graph, ratio, l_signed)
            };
            if a.abs() >= 1e-13 {
                let c = cx.graph.num(Number::from(a));
                terms.push(mul(cx.graph, &[c, radial, angular, cos_node]));
            }
            if m > 0 && b.abs() >= 1e-13 {
                let c = cx.graph.num(Number::from(b));
                terms.push(mul(cx.graph, &[c, radial, angular, sin_node]));
            }
        }
    }
    Some(add(cx.graph, &terms))
}

/// Whether the node is `oo`.
fn is_infinity(
    cx: &Cx<'_>,
    node: NodeId,
) -> bool {
    cx.graph
        .ops()
        .lookup("oo")
        .is_some_and(|op| cx.graph.op(node) == op)
}

#[allow(clippy::too_many_lines)]
fn solve_layout(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    b: NodeId,
    a0: NodeId,
) -> Option<NodeId> {
    let n = p.dimension();
    // The conditions at the radii; nothing else on the angles.
    if conditions
        .0
        .iter()
        .any(|c| Some(c.on) == layout.theta || Some(c.on) == layout.phi)
    {
        return None;
    }
    let (at_infinity, finite): (Vec<&Condition>, Vec<&Condition>) = conditions
        .0
        .iter()
        .filter(|c| c.on == layout.r)
        .partition(|c| is_infinity(cx, c.point));
    if let [far] = at_infinity.as_slice() {
        return exterior(cx, p, conditions, layout, time_index, (b, a0), &finite, far);
    }
    if !at_infinity.is_empty() {
        return None;
    }
    let (inner, outer) = match finite.as_slice() {
        | [first, second] => {
            if let Some(found) = annulus(
                cx,
                p,
                conditions,
                layout,
                time_index,
                (b, a0),
                [*first, *second],
            ) {
                return Some(found);
            }
            let (va, vb) = (
                super::util::sample(cx.graph, first.point, 0)?,
                super::util::sample(cx.graph, second.point, 0)?,
            );
            if va <= vb {
                (Some(*first), *second)
            } else {
                (Some(*second), *first)
            }
        },
        | [outer] => (None, *outer),
        | _ => return None,
    };
    let radius = outer.point;
    let end = End::from_condition(cx, p, outer)?;
    let inner_end = match inner {
        | Some(c) => Some((c.point, End::from_condition(cx, p, c)?)),
        | None => None,
    };
    // Harmonic functions with Dirichlet data: the closed series.
    if inner.is_none()
        && time_index.is_none()
        && layout.z.is_none()
        && layout.phi.is_none()
        && p.homogeneous(cx.graph)
        && is_zero_number(cx.graph, a0)
        && end.kind == Kind::Dirichlet
        && cx.graph.number_of(end.p).is_some()
        && let Some(theta) = layout.theta
    {
        let r = p.vars[layout.r];
        let found = if layout.sphere {
            classical::laplace_ball(cx, end.value, radius, r, p.vars[theta])
        } else {
            classical::laplace_disk(cx, end.value, radius, r, p.vars[theta])
        };
        if let Some(f) = found {
            return Some(f);
        }
    }
    if inner.is_none()
        && time_index.is_none()
        && layout.z.is_none()
        && layout.phi.is_none()
        && !p.homogeneous(cx.graph)
        && is_zero_number(cx.graph, a0)
        && end.kind == Kind::Dirichlet
        && let Some(found) = radial_poisson(cx, p, layout, b, end, radius)
    {
        return Some(found);
    }
    // Harmonic functions with Robin or Neumann data, or without axial
    // symmetry: each harmonic times its power of r.
    if inner.is_none()
        && time_index.is_none()
        && layout.z.is_none()
        && p.homogeneous(cx.graph)
        && cx.is_zero(a0)
        && layout.theta.is_some()
        && let Some(found) = harmonic_interior(cx, p, conditions, layout, (b, a0), end, radius)
    {
        return Some(found);
    }
    // Time structure.
    let (mut a2, mut a1) = (cx.graph.int(0), cx.graph.int(0));
    if let Some(t) = time_index {
        a1 = p.coefficient(cx.graph, &p.unit(t, 1));
        a2 = p.coefficient(cx.graph, &p.unit(t, 2));
        if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
            return None;
        }
    }
    let order = if is_zero_number(cx.graph, a2) {
        u32::from(!is_zero_number(cx.graph, a1))
    } else {
        2
    };
    let time = time_index.map(|t| Time {
        var: p.vars[t],
        order,
        a2,
        a1,
    });
    // Axes in the order azimuth, angle, radius, axis.
    let mut axes = Vec::new();
    if let Some(j) = layout.phi {
        axes.push(Axis {
            var: p.vars[j],
            shape: Shape::Angle,
        });
    }
    if let Some(j) = layout.theta {
        axes.push(Axis {
            var: p.vars[j],
            shape: if layout.sphere {
                Shape::Polar {
                    azimuth: layout.phi.is_some(),
                }
            } else {
                Shape::Angle
            },
        });
    }
    let radial_level = axes.len();
    axes.push(Axis {
        var: p.vars[layout.r],
        shape: Shape::Radial {
            radius,
            inner: inner_end,
            end,
            sphere: layout.sphere,
            has_angle: layout.theta.is_some(),
        },
    });
    let z_interval = match layout.z {
        | Some(j) => {
            let iv = interval(cx, p, conditions, j)?;
            axes.push(Axis {
                var: p.vars[j],
                shape: Shape::Interval(iv.clone()),
            });
            Some((axes.len() - 1, iv))
        },
        | None => None,
    };
    // Conditions: r, z and the initial ones only.
    for c in &conditions.0 {
        let allowed = c.on == layout.r || Some(c.on) == layout.z || Some(c.on) == time_index;
        if !allowed {
            return None;
        }
    }
    let mut items: Vec<Item> = time_items(cx, p, conditions, time_index, order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item {
            role: Role::Source,
            h: p.source,
            special: None,
        });
    }
    // Boundary data at the radii: Green's identity with weight r^(d-1)
    // (opposite outward normals on the inner and the outer circle).
    let mut boundaries = vec![(radius, end, 1_i64)];
    if let Some((a, a_end)) = inner_end {
        boundaries.push((a, a_end, -1));
    }
    for (point, e, side) in boundaries {
        if cx.is_zero(e.value) {
            continue;
        }
        let divisor = if e.kind == Kind::Dirichlet {
            e.p
        } else {
            e.q
        };
        let sign = cx.graph.int(side);
        let scaled = mul(cx.graph, &[sign, b, e.value]);
        let scaled = div(cx, scaled, divisor);
        let scaled = cx.simplify(scaled);
        let weight = powi(cx.graph, point, if layout.sphere { 2 } else { 1 });
        items.push(Item {
            role: Role::Source,
            h: scaled,
            special: Some(Special {
                axis: radial_level,
                kind: e.kind,
                point,
                weight,
            }),
        });
    }
    if let Some((level, iv)) = &z_interval {
        items.extend(boundary_items(cx, iv, *level, b));
    }
    let engine = Engine {
        axes,
        time,
        symbol: Symbol::Laplacian { b, a0 },
        face: None,
        exterior: None,
    };
    finish(cx, p, conditions, &engine, &items, None, None)
}

/// Laplace's and Helmholtz's equations outside a disk or a ball (and the
/// heat and wave equations outside a ball for radial data): the condition
/// at `r = oo` selects the bounded decaying solution (`u → c`) or, for the
/// Helmholtz equation `Δu + k² u = 0` with `k² > 0`, the outgoing one
/// (Sommerfeld radiation condition, time factor `e^{-iωt}`: Hankel functions
/// of the first kind); the angular dependence is expanded as in the interior
/// and each harmonic is multiplied by its decaying radial profile.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn exterior(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    coefficients: (NodeId, NodeId),
    finite: &[&Condition],
    far: &Condition,
) -> Option<NodeId> {
    let (b, a0) = coefficients;
    let [near] = finite else {
        return None;
    };
    if layout.z.is_some() || far.robin.is_some() || far.derivative.iter().any(|&d| d != 0) {
        return None;
    }
    let end = End::from_condition(cx, p, near)?;
    let radius = near.point;
    if time_index.is_some() {
        return exterior_evolution(cx, p, conditions, layout, time_index, b, a0, end, far);
    }
    if !p.homogeneous(cx.graph) {
        return None;
    }
    let laplace = cx.is_zero(a0);
    let c_inf = if laplace {
        far.value
    } else {
        cx.graph.int(0)
    };
    if !laplace && !cx.is_zero(far.value) {
        return None;
    }
    // u = c + w, w decaying: p w + q w_r = value - p c on r = R.
    let data = {
        let pc = mul(cx.graph, &[end.p, c_inf]);
        let gap = sub(cx.graph, end.value, pc);
        cx.simplify(gap)
    };
    let wavenumber2 = {
        let q = div(cx, a0, b);
        cx.simplify(q)
    };
    let mut axes = Vec::new();
    if let Some(j) = layout.phi {
        axes.push(Axis {
            var: p.vars[j],
            shape: Shape::Angle,
        });
    }
    if let Some(j) = layout.theta {
        axes.push(Axis {
            var: p.vars[j],
            shape: if layout.sphere {
                Shape::Polar {
                    azimuth: layout.phi.is_some(),
                }
            } else {
                Shape::Angle
            },
        });
    }
    // In a plane the bounded harmonic has the limit a₀ (the mean of the
    // data): the value given at infinity must agree with it.
    if laplace && !layout.sphere {
        let mean = match layout.theta {
            | Some(j) => {
                let theta = p.vars[j];
                let two = cx.graph.int(2);
                let pi = super::util::pi(cx)?;
                let two_pi = mul(cx.graph, &[two, pi]);
                let zero = cx.graph.int(0);
                let integral = super::util::defint(cx, data, theta, zero, two_pi)?;
                let q = div(cx, integral, two_pi);
                cx.simplify(q)
            },
            | None => data,
        };
        let cleared = cx.is_zero(mean)
            || (0..3_u32)
                .all(|k| super::util::sample(cx.graph, mean, k).is_some_and(|v| v.abs() < 1e-9));
        if !cleared {
            return None;
        }
    }
    let exterior = Exterior {
        var: p.vars[layout.r],
        radius,
        end,
        sphere: layout.sphere,
        interior: false,
        wavenumber2,
    };
    let engine = Engine {
        axes,
        time: None,
        symbol: Symbol::Laplacian { b, a0 },
        face: None,
        exterior: Some(exterior),
    };
    let items = [Item {
        role: Role::Displacement,
        h: data,
        special: None,
    }];
    let kept = Conditions(
        conditions
            .0
            .iter()
            .filter(|c| c.on != layout.r || !is_infinity(cx, c.point))
            .cloned()
            .collect(),
    );
    let offset = if cx.is_zero(c_inf) {
        None
    } else {
        Some(c_inf)
    };
    finish(cx, p, &kept, &engine, &items, None, offset)
}

/// Laplace's equation inside a disk or ball with data of any kind on the
/// boundary: the angular expansion times `r^n` (`r^l`), normalised by the
/// boundary condition.
fn harmonic_interior(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    coefficients: (NodeId, NodeId),
    end: End,
    radius: NodeId,
) -> Option<NodeId> {
    let (b, a0) = coefficients;
    let mut axes = Vec::new();
    if let Some(j) = layout.phi {
        axes.push(Axis {
            var: p.vars[j],
            shape: Shape::Angle,
        });
    }
    let j = layout.theta?;
    axes.push(Axis {
        var: p.vars[j],
        shape: if layout.sphere {
            Shape::Polar {
                azimuth: layout.phi.is_some(),
            }
        } else {
            Shape::Angle
        },
    });
    let wavenumber2 = cx.graph.int(0);
    let interior = Exterior {
        var: p.vars[layout.r],
        radius,
        end,
        sphere: layout.sphere,
        interior: true,
        wavenumber2,
    };
    let engine = Engine {
        axes,
        time: None,
        symbol: Symbol::Laplacian { b, a0 },
        face: None,
        exterior: Some(interior),
    };
    let items = [Item {
        role: Role::Displacement,
        h: end.value,
        special: None,
    }];
    finish(cx, p, conditions, &engine, &items, None, None)
}

/// The heat (or wave) equation outside a ball for radial data: `v = r u`
/// solves the one-dimensional equation on `x = r - R > 0` with the
/// Dirichlet value `R g(t)`, solved by the half-line formulas.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn exterior_evolution(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    b: NodeId,
    a0: NodeId,
    end: End,
    far: &Condition,
) -> Option<NodeId> {
    let t_index = time_index?;
    // Sources are handled (by Duhamel's principle) only on the Weber route.
    let weber = !layout.sphere || layout.theta.is_some();
    if (!p.homogeneous(cx.graph) && !weber) || !cx.is_zero(a0) || !cx.is_zero(far.value) {
        return None;
    }
    let n = p.dimension();
    let c_second = p.coefficient(cx.graph, &p.unit(t_index, 2));
    let order = if is_zero_number(cx.graph, c_second) {
        1
    } else {
        2
    };
    let c_time = p.coefficient(cx.graph, &p.unit(t_index, order));
    let diffusivity = {
        let q = div(cx, b, c_time);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let (r, t) = (p.vars[layout.r], p.vars[t_index]);
    let radius = end.point;
    if layout.sphere && layout.theta.is_none() && order == 2 && end.kind != Kind::Dirichlet {
        // Waves outside a ball with a Neumann or Robin condition: v = r u on
        // the half-line with v' + (p/q - 1/R) v = 0.
        if !cx.is_zero(end.value) {
            return None;
        }
        let ratio = div(cx, end.p, end.q);
        let inverse_radius = powi(cx.graph, radius, -1);
        let rho = sub(cx.graph, ratio, inverse_radius);
        let rho = cx.simplify(rho);
        let (x, _) = super::dummy(cx, p, "x");
        let shifted = add(cx.graph, &[x, radius]);
        let zero_index = vec![0; n];
        let mut data = Vec::new();
        for derivative in [zero_index.clone(), p.unit(t_index, 1)] {
            data.push(
                conditions
                    .find(cx.graph, t_index, &derivative, None)
                    .map(|c| {
                        let f = cx.graph.substitute(c.value, r, shifted);
                        let value = mul(cx.graph, &[shifted, f]);
                        cx.simplify(value)
                    }),
            );
        }
        if conditions.0.iter().any(|c| {
            c.on == t_index && c.derivative != zero_index && c.derivative != p.unit(t_index, 1)
        }) {
            return None;
        }
        let f = data.first().copied().flatten()?;
        let g = data.get(1).copied().flatten();
        let speed = sqrt(cx, diffusivity)?;
        let v = classical::robin_wave_half_line(cx, p, (x, t), speed, rho, f, g)?;
        let offset = sub(cx.graph, r, radius);
        let v = cx.graph.substitute(v, x, offset);
        let u = div(cx, v, r);
        return Some(cx.simplify(u));
    }
    if !layout.sphere || layout.theta.is_some() {
        return weber_exterior(
            cx,
            p,
            conditions,
            layout,
            (t_index, order),
            diffusivity,
            &end,
        );
    }
    let _ = (r, t);
    // The scalar problem for w(x, t) = (x + R) u(x + R, t).
    let (x, _) = super::dummy(cx, p, "x");
    let w_symbol = cx.graph.interner_mut().fresh_symbol("w");
    let w_head = cx.graph.symbol_node(w_symbol);
    let unknown = cx.graph.node(core::APPLY, &[w_head, x, t]);
    let vars = [x, t];
    let lhs = super::util::jet(cx, unknown, &vars, &[0, order])?;
    let w_xx = super::util::jet(cx, unknown, &vars, &[2, 0])?;
    let rhs = mul(cx.graph, &[diffusivity, w_xx]);
    let equation = cx.graph.node(core::EQ, &[lhs, rhs]);
    let scalar = Problem::parse(cx, equation, unknown)?;
    let shifted = add(cx.graph, &[x, radius]);
    let zero = cx.graph.int(0);
    let mut list = Vec::new();
    // v(R) = R u(R, t) = R value / p.
    let boundary = {
        let q = div(cx, end.value, end.p);
        let q = mul(cx.graph, &[radius, q]);
        cx.simplify(q)
    };
    if end.kind == Kind::Dirichlet {
        list.push(Condition {
            on: 0,
            point: zero,
            derivative: vec![0, 0],
            value: boundary,
            robin: None,
        });
    } else {
        // p u + q u_r = g at r = R with v = r u: v' + (p/q - 1/R) v = R g / q.
        let ratio = div(cx, end.p, end.q);
        let inverse_radius = powi(cx.graph, radius, -1);
        let rho = sub(cx.graph, ratio, inverse_radius);
        let rho = cx.simplify(rho);
        let value = {
            let q = div(cx, end.value, end.q);
            let q = mul(cx.graph, &[radius, q]);
            cx.simplify(q)
        };
        let robin = if cx.is_zero(rho) {
            None
        } else {
            Some(rho)
        };
        list.push(Condition {
            on: 0,
            point: zero,
            derivative: vec![1, 0],
            value,
            robin,
        });
    }
    let zero_index = vec![0; n];
    for (derivative, vector) in [
        (zero_index.clone(), vec![0, 0]),
        (p.unit(t_index, 1), vec![0, 1]),
    ] {
        if let Some(c) = conditions.find(cx.graph, t_index, &derivative, None) {
            let f = cx.graph.substitute(c.value, r, shifted);
            let value = mul(cx.graph, &[shifted, f]);
            let value = cx.simplify(value);
            list.push(Condition {
                on: 1,
                point: zero,
                derivative: vector,
                value,
                robin: None,
            });
        }
    }
    if conditions.0.iter().any(|c| {
        c.on == t_index && c.derivative != zero_index && c.derivative != p.unit(t_index, 1)
    }) {
        return None;
    }
    let found = super::solve(cx, &scalar, &Conditions(list), Method::Any)?;
    let &[_, w] = cx.graph.children(found) else {
        return None;
    };
    // u = w(r - R, t) / r.
    let offset = sub(cx.graph, r, radius);
    let w = cx.graph.substitute(w, x, offset);
    let u = div(cx, w, r);
    Some(cx.simplify(u))
}

/// A harmonic of the data: its angular factor and radial coefficient.
struct Harmonic {
    index: i64,
    /// The angular factor (`cos nθ`, `sin nθ` or `P_l(cos θ)`).
    angular: NodeId,
    radial: NodeId,
}

/// The decomposition of `data(r, θ, φ)` into the real spherical harmonics
/// `P_l^m(cos θ) cos mφ` and `P_l^m(cos θ) sin mφ` (`l ≤ 8`, the tail must
/// vanish): the coefficient is the exact double integral divided by
/// `2π(l+m)!/((2l+1)(l-m)!) · (1 or 1/2)`.
fn spherical_harmonics(
    cx: &mut Cx<'_>,
    data: NodeId,
    theta: NodeId,
    phi: NodeId,
) -> Option<Vec<Harmonic>> {
    const TOP: i64 = 9;
    let factorial = |n: i64| -> i64 { (1..=n).product() };
    let pi = super::util::pi(cx)?;
    let (zero, two) = (cx.graph.int(0), cx.graph.int(2));
    let two_pi = mul(cx.graph, &[two, pi]);
    let (c, s) = (call(cx, "cos", &[theta])?, call(cx, "sin", &[theta])?);
    let mut out = Vec::new();
    for l in 0..=TOP {
        for m in 0..=l {
            let kinds: &[&str] = if m == 0 {
                &["cos"]
            } else {
                &["cos", "sin"]
            };
            for &name in kinds {
                let m_node = cx.graph.int(m);
                let m_phi = mul(cx.graph, &[m_node, phi]);
                let wave = call(cx, name, &[m_phi])?;
                let body = mul(cx.graph, &[data, wave]);
                let around = defint(cx, body, phi, zero, two_pi)?;
                let around = cx.simplify(around);
                if super::util::contains_op(cx.graph, around, "defint") {
                    return None;
                }
                if cx.is_zero(around) {
                    continue;
                }
                let mode = super::modes::associated_mode(cx, theta, l, m)?;
                let weighted = mul(cx.graph, &[around, mode]);
                let integral = super::modes::legendre_integral(cx, weighted, c, Some(s), theta)?;
                let phi_norm = if m == 0 { two_pi } else { pi };
                let scale = cx.graph.num(Number::fraction(
                    (2 * l + 1) * factorial(l - m),
                    2 * factorial(l + m),
                )?);
                let scaled = mul(cx.graph, &[scale, integral]);
                let coefficient = div(cx, scaled, phi_norm);
                let coefficient = cx.simplify(coefficient);
                if super::util::contains_op(cx.graph, coefficient, "defint") {
                    return None;
                }
                if cx.is_zero(coefficient) {
                    continue;
                }
                if l == TOP {
                    return None;
                }
                let angular = mul(cx.graph, &[mode, wave]);
                out.push(Harmonic {
                    index: l,
                    angular,
                    radial: coefficient,
                });
            }
        }
    }
    Some(out)
}

/// The angular decomposition of `data(r, θ)`: Fourier harmonics on a disk
/// (up to 12, the tail must vanish) or Legendre polynomials on a ball
/// (axisymmetric, up to degree 9).
fn harmonics(
    cx: &mut Cx<'_>,
    data: NodeId,
    theta: Option<NodeId>,
    phi: Option<NodeId>,
    sphere: bool,
) -> Option<Vec<Harmonic>> {
    let Some(theta) = theta else {
        let one = cx.graph.int(1);
        return Some(vec![Harmonic {
            index: 0,
            angular: one,
            radial: data,
        }]);
    };
    let mut out = Vec::new();
    let pi = super::util::pi(cx)?;
    if sphere && let Some(phi) = phi {
        return spherical_harmonics(cx, data, theta, phi);
    }
    if sphere {
        let (c, s) = (call(cx, "cos", &[theta])?, call(cx, "sin", &[theta])?);
        for l in 0..=10_i64 {
            let mode = super::modes::associated_mode(cx, theta, l, 0)?;
            let body = mul(cx.graph, &[data, mode]);
            let integral = super::modes::legendre_integral(cx, body, c, Some(s), theta)?;
            let scale = cx.graph.num(Number::fraction(2 * l + 1, 2)?);
            let coefficient = mul(cx.graph, &[scale, integral]);
            let coefficient = cx.simplify(coefficient);
            if super::util::contains_op(cx.graph, coefficient, "defint") {
                return None;
            }
            if !cx.is_zero(coefficient) {
                if l == 10 {
                    return None;
                }
                let l_node = cx.graph.int(l);
                let leg = call(cx, "legendre", &[l_node, c])?;
                out.push(Harmonic {
                    index: l,
                    angular: leg,
                    radial: coefficient,
                });
            }
        }
        return Some(out);
    }
    let zero = cx.graph.int(0);
    let two_pi = {
        let two = cx.graph.int(2);
        mul(cx.graph, &[two, pi])
    };
    for k in 0..=12_i64 {
        let k_node = cx.graph.int(k);
        let k_theta = mul(cx.graph, &[k_node, theta]);
        let kinds: Vec<(&str, bool)> = if k == 0 {
            vec![("cos", true)]
        } else {
            vec![("cos", true), ("sin", false)]
        };
        for (name, _) in kinds {
            let wave = call(cx, name, &[k_theta])?;
            let body = mul(cx.graph, &[data, wave]);
            let integral = defint(cx, body, theta, zero, two_pi)?;
            let norm = if k == 0 { two_pi } else { pi };
            let coefficient = div(cx, integral, norm);
            let coefficient = cx.simplify(coefficient);
            if super::util::contains_op(cx.graph, coefficient, "defint") {
                return None;
            }
            if !cx.is_zero(coefficient) {
                if k > 10 {
                    return None;
                }
                out.push(Harmonic {
                    index: k,
                    angular: wave,
                    radial: coefficient,
                });
            }
        }
    }
    Some(out)
}

/// The evolution outside a disk or an axisymmetric ball with a homogeneous
/// Dirichlet, Neumann or Robin condition, by Weber's transform of the order
/// `ν` of each harmonic (`ν = n` on a disk, `ν = l + 1/2` for the radial
/// factor `r^{-1/2} F` of a ball):
///
/// `F = ∫_0^∞ k M(k) Φ_k(r) [∫_R^∞ ρ F₀(ρ) Φ_k(ρ) dρ] / (A² + B²) dk`
///
/// with `Φ_k = A J_ν(k r) - B Y_ν(k r)`, `A = p Y_ν(kR) + q k Y_ν'(kR)`,
/// `B = p J_ν(kR) + q k J_ν'(kR)` so that `p Φ + q Φ' = 0` at `r = R`
/// (`p - q/(2R)` for a ball); `M` is `e^{-D k² t}` (first order in time) or
/// `cos(√D k t)` and `sin(√D k t)/(√D k)` for displacement and velocity.
///
/// A source that is a sum of separable terms `g_i(t) s_i(x)` (see
/// [`separate`]) is treated per term by Duhamel's principle: the Weber
/// coefficient of `s_i` is formed once per mode and the time kernel is
/// integrated against `g_i` alone (`∫_0^t M(t - τ) g_i(τ) dτ`, resolved on
/// that small integrand), instead of a triple `ρ, k, τ` integral per mode; a
/// non-separable source falls back to the nested form.
#[allow(clippy::too_many_lines)]
fn weber_exterior(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    (t_index, order): (usize, u32),
    diffusivity: NodeId,
    end: &End,
) -> Option<NodeId> {
    if !cx.is_zero(end.value) {
        return None;
    }
    let n = p.dimension();
    let zero_index = vec![0; n];
    let displacement = conditions.find(cx.graph, t_index, &zero_index, None)?;
    let velocity = if order == 2 {
        conditions.find(cx.graph, t_index, &p.unit(t_index, 1), None)
    } else {
        None
    };
    if conditions.0.iter().filter(|c| c.on == t_index).count()
        != 1 + usize::from(velocity.is_some())
    {
        return None;
    }
    let (r, t) = (p.vars[layout.r], p.vars[t_index]);
    let theta = layout.theta.map(|j| p.vars[j]);
    let radius = end.point;
    let (k, _) = super::dummy(cx, p, "k");
    let (rho, _) = super::dummy(cx, p, "rho");
    let (zero, one) = (cx.graph.int(0), cx.graph.int(1));
    let oo = super::util::infinity(cx)?;
    let sqrt_d = super::util::sqrt(cx, diffusivity)?;
    // Data: displacement (0), velocity (1) and the source, whose response is
    // the Duhamel integral of the velocity kernel (second order) or of the
    // initial-value solution (first order) (2).
    let mut pairs = vec![(displacement.value, 0_u8, one)];
    if let Some(v) = velocity {
        pairs.push((v.value, 1, one));
    }
    let (tau, _) = super::dummy(cx, p, "tau");
    let elapsed = sub(cx.graph, t, tau);
    if !p.homogeneous(cx.graph) {
        // `a ∂_t^o u = b Δu - source`.
        let a = p.coefficient(cx.graph, &p.unit(t_index, order));
        let forcing = {
            let q = div(cx, p.source, a);
            let q = neg(cx.graph, q);
            cx.simplify(q)
        };
        // A steady source integrates its kernel in closed form.
        let steady = !cx
            .graph
            .symbol_of(t)
            .is_some_and(|s| cx.graph.depends_on(cx.graph.find(forcing), s));
        let space: Vec<NodeId> = [Some(r), theta, layout.phi.map(|j| p.vars[j])]
            .into_iter()
            .flatten()
            .collect();
        if steady {
            pairs.push((forcing, 3, one));
        } else if let Some(parts) = separate(cx, forcing, t, &space) {
            // A source `Σ g_i(t) s_i(x)`: the Weber coefficient of `s_i` is
            // computed once per mode and only the time kernel is integrated
            // against `g_i`, instead of a triple integral per mode.
            for (time_factor, spatial) in parts {
                pairs.push((spatial, 4, cx.graph.substitute(time_factor, t, tau)));
            }
        } else {
            pairs.push((cx.graph.substitute(forcing, t, tau), 2, one));
        }
    }
    let mut terms = Vec::new();
    for (value, kind, time_factor) in pairs {
        let is_displacement = kind == 0;
        let time_arg = if kind == 2 || kind == 4 {
            elapsed
        } else {
            t
        };
        if cx.is_zero(value) {
            continue;
        }
        for harmonic in harmonics(
            cx,
            value,
            theta,
            layout.phi.map(|j| p.vars[j]),
            layout.sphere,
        )? {
            let nu = if layout.sphere {
                cx.graph.num(Number::fraction(2 * harmonic.index + 1, 2)?)
            } else {
                cx.graph.int(harmonic.index)
            };
            let nu_up = add(cx.graph, &[nu, one]);
            // p - q/(2R) for the ball's r^{-1/2} factor.
            let p_shift = if layout.sphere {
                let half = cx.graph.num(Number::fraction(-1, 2)?);
                let q_over_r = div(cx, end.q, radius);
                let shift = mul(cx.graph, &[half, q_over_r]);
                add(cx.graph, &[end.p, shift])
            } else {
                end.p
            };
            let kr0 = mul(cx.graph, &[k, radius]);
            let z = |cx: &mut Cx<'_>, second: bool, order_node: NodeId, arg: NodeId| {
                call(
                    cx,
                    if second {
                        "bessely"
                    } else {
                        "besselj"
                    },
                    &[order_node, arg],
                )
            };
            // p Z + q k Z' with Z' = -Z_{ν+1} + (ν/z) Z_ν at z = kR.
            let coefficient = |cx: &mut Cx<'_>, second: bool| -> Option<NodeId> {
                let (z_nu, z_up) = (z(cx, second, nu, kr0)?, z(cx, second, nu_up, kr0)?);
                let first = mul(cx.graph, &[p_shift, z_nu]);
                let ratio = div(cx, nu, radius);
                let qr = mul(cx.graph, &[end.q, ratio, z_nu]);
                let qk = mul(cx.graph, &[end.q, k, z_up]);
                let nqk = neg(cx.graph, qk);
                Some(add(cx.graph, &[first, qr, nqk]))
            };
            let (a_coeff, b_coeff) = (coefficient(cx, true)?, coefficient(cx, false)?);
            let phi = |cx: &mut Cx<'_>, x: NodeId| -> Option<NodeId> {
                let kx = mul(cx.graph, &[k, x]);
                let (j, y) = (z(cx, false, nu, kx)?, z(cx, true, nu, kx)?);
                let first = mul(cx.graph, &[a_coeff, j]);
                let second = mul(cx.graph, &[b_coeff, y]);
                Some(sub(cx.graph, first, second))
            };
            let (phi_r, phi_rho) = (phi(cx, r)?, phi(cx, rho)?);
            let mut f0 = cx.graph.substitute(harmonic.radial, r, rho);
            if layout.sphere {
                let half = cx.graph.num(Number::fraction(1, 2)?);
                let root = pow(cx.graph, rho, half);
                f0 = mul(cx.graph, &[root, f0]);
            }
            let inner_body = mul(cx.graph, &[rho, f0, phi_rho]);
            let inner = defint(cx, inner_body, rho, radius, oo)?;
            let multiplier = if order == 1 {
                let k2 = powi(cx.graph, k, 2);
                let e = mul(cx.graph, &[diffusivity, k2, time_arg]);
                let e = neg(cx.graph, e);
                call(cx, "exp", &[e])?
            } else {
                let ck = mul(cx.graph, &[sqrt_d, k]);
                let ckt = mul(cx.graph, &[ck, time_arg]);
                if is_displacement {
                    call(cx, "cos", &[ckt])?
                } else {
                    let sn = call(cx, "sin", &[ckt])?;
                    div(cx, sn, ck)
                }
            };
            let multiplier = if kind == 3 {
                // `∫_0^t M(t - τ) dτ`: `(1 - e^{-D k² t})/(D k²)` or `(1 - cos(c k t))/(c k)²`.
                if order == 1 {
                    let k_squared = powi(cx.graph, k, 2);
                    let dk2 = mul(cx.graph, &[diffusivity, k_squared]);
                    let e = mul(cx.graph, &[dk2, t]);
                    let e = neg(cx.graph, e);
                    let decay = call(cx, "exp", &[e])?;
                    let gap = sub(cx.graph, one, decay);
                    div(cx, gap, dk2)
                } else {
                    let ck = mul(cx.graph, &[sqrt_d, k]);
                    let ckt = mul(cx.graph, &[ck, t]);
                    let cs = call(cx, "cos", &[ckt])?;
                    let gap = sub(cx.graph, one, cs);
                    let ck2 = powi(cx.graph, ck, 2);
                    div(cx, gap, ck2)
                }
            } else if kind == 4 {
                // `∫_0^t M(t - τ) g(τ) dτ`, one-dimensional in τ.
                // Resolved here, on the small integrand, so that the engine
                // does not have to attack it inside the k-integral.
                let body = mul(cx.graph, &[multiplier, time_factor]);
                let integral = defint(cx, body, tau, zero, t)?;
                let closed = cx.simplify(integral);
                if super::util::contains_op(cx.graph, closed, "defint") {
                    integral
                } else {
                    closed
                }
            } else {
                multiplier
            };
            let norm = {
                let (a2, b2) = (powi(cx.graph, a_coeff, 2), powi(cx.graph, b_coeff, 2));
                let total = add(cx.graph, &[a2, b2]);
                powi(cx.graph, total, -1)
            };
            let outer_body = mul(cx.graph, &[k, multiplier, phi_r, inner, norm]);
            let radial = defint(cx, outer_body, k, zero, oo)?;
            let radial = if kind == 2 {
                defint(cx, radial, tau, zero, t)?
            } else {
                radial
            };
            let radial = if layout.sphere {
                let minus_half = cx.graph.num(Number::fraction(-1, 2)?);
                let root = pow(cx.graph, r, minus_half);
                mul(cx.graph, &[root, radial])
            } else {
                radial
            };
            terms.push(mul(cx.graph, &[harmonic.angular, radial]));
        }
    }
    if terms.is_empty() {
        return None;
    }
    Some(add(cx.graph, &terms))
}

/// Splits a source into separable terms `Σ g_i(t) s_i(x)` (`g_i` free of the
/// space variables, `s_i` free of `t`), or `None` when a factor couples
/// time and space (`exp(t r)`, `(t + r)^2`) or there are too many terms.
fn separate(
    cx: &mut Cx<'_>,
    f: NodeId,
    time: NodeId,
    space: &[NodeId],
) -> Option<Vec<(NodeId, NodeId)>> {
    const MAX_TERMS: usize = 8;
    let depends = |cx: &Cx<'_>, node: NodeId, var: NodeId| {
        cx.graph
            .symbol_of(var)
            .is_some_and(|s| cx.graph.depends_on(cx.graph.find(node), s))
    };
    let one = cx.graph.int(1);
    if !depends(cx, f, time) {
        return Some(vec![(one, f)]);
    }
    if !space.iter().any(|&v| depends(cx, f, v)) {
        return Some(vec![(f, one)]);
    }
    let op = cx.graph.op(f);
    let children = cx.graph.children(f).to_vec();
    if op == core::ADD {
        let mut out = Vec::new();
        for child in children {
            out.extend(separate(cx, child, time, space)?);
            if out.len() > MAX_TERMS {
                return None;
            }
        }
        return Some(out);
    }
    if op == core::MUL {
        let mut acc = vec![(one, one)];
        for child in children {
            let parts = separate(cx, child, time, space)?;
            let mut next = Vec::with_capacity(acc.len().saturating_mul(parts.len()));
            for &(g1, s1) in &acc {
                for &(g2, s2) in &parts {
                    next.push((mul(cx.graph, &[g1, g2]), mul(cx.graph, &[s1, s2])));
                }
            }
            if next.len() > MAX_TERMS {
                return None;
            }
            acc = next;
        }
        return Some(acc);
    }
    None
}

/// `Δu = f(r)` with a polynomial `f` and Dirichlet data on `r = R`: the
/// polynomial particular solution plus the harmonic function with the
/// remaining data.
fn radial_poisson(
    cx: &mut Cx<'_>,
    p: &Problem,
    layout: Layout,
    b: NodeId,
    end: End,
    radius: NodeId,
) -> Option<NodeId> {
    let r = p.vars[layout.r];
    if !cx.graph.number_of(end.p).is_some_and(Number::is_one) {
        return None;
    }
    // Δu = f with f = -source / b, a polynomial in r alone.
    let f = {
        let q = div(cx, p.source, b);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    for (j, &v) in p.vars.iter().enumerate() {
        let depends = cx
            .graph
            .symbol_of(v)
            .is_some_and(|s| cx.graph.depends_on(cx.graph.find(f), s));
        if depends && j != layout.r {
            return None;
        }
    }
    let zero = cx.graph.int(0);
    let mut particular = Vec::new();
    let mut polynomial = Vec::new();
    let mut derivative_k = f;
    let mut factorial = 1_i64;
    for k in 0..=8_i64 {
        if k > 0 {
            derivative_k = derivative(cx.graph, derivative_k, r)?;
            factorial = factorial.checked_mul(k)?;
        }
        let at_zero = cx.graph.substitute(derivative_k, r, zero);
        let coefficient = {
            let fact = cx.graph.int(factorial);
            let c = div(cx, at_zero, fact);
            cx.simplify(c)
        };
        if cx.is_zero(coefficient) {
            continue;
        }
        let power = powi(cx.graph, r, k);
        polynomial.push(mul(cx.graph, &[coefficient, power]));
        let raised = powi(cx.graph, r, k + 2);
        let denominator = if layout.sphere {
            (k + 2) * (k + 3)
        } else {
            (k + 2) * (k + 2)
        };
        let den = cx.graph.int(denominator);
        let scaled = div(cx, coefficient, den);
        particular.push(mul(cx.graph, &[scaled, raised]));
    }
    let rebuilt = add(cx.graph, &polynomial);
    let gap = sub(cx.graph, f, rebuilt);
    if !cx.is_zero(gap) {
        return None;
    }
    let particular = add(cx.graph, &particular);
    let particular = cx.simplify(particular);
    // Remaining data on r = R.
    let at_radius = cx.graph.substitute(particular, r, radius);
    let data = sub(cx.graph, end.value, at_radius);
    let data = cx.simplify(data);
    let harmonic = match layout.theta {
        | Some(j) if layout.sphere => classical::laplace_ball(cx, data, radius, r, p.vars[j])?,
        | Some(j) => classical::laplace_disk(cx, data, radius, r, p.vars[j])?,
        | None => {
            if cx
                .graph
                .depends_on(cx.graph.find(data), cx.graph.symbol_of(r)?)
            {
                return None;
            }
            data
        },
    };
    let total = add(cx.graph, &[particular, harmonic]);
    let total = cx.simplify(total);
    super::verified(cx, p, total).then_some(total)
}

/// Laplace's equation in a ring (or spherical shell) with Dirichlet data on
/// both boundary circles: `a₀ + b₀ ln r + Σ (a_n r^n + b_n r^-n) cos nθ + …`
/// (`r^l`, `r^(-l-1)` and `P_l(cos θ)` in a shell), for data that are
/// trigonometric polynomials (polynomials in `cos θ`).
fn annulus(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    layout: Layout,
    time_index: Option<usize>,
    coefficients: (NodeId, NodeId),
    ends: [&super::Condition; 2],
) -> Option<NodeId> {
    let (_, a0) = coefficients;
    let theta = p.vars[layout.theta?];
    if time_index.is_some()
        || layout.z.is_some()
        || !p.homogeneous(cx.graph)
        || !is_zero_number(cx.graph, a0)
        || conditions.0.len() != 2
    {
        return None;
    }
    let r = p.vars[layout.r];
    let [first, second] = ends;
    let (inner, outer) = {
        let (va, vb) = (
            super::util::sample(cx.graph, first.point, 0)?,
            super::util::sample(cx.graph, second.point, 0)?,
        );
        if va <= vb {
            (first, second)
        } else {
            (second, first)
        }
    };
    for c in [inner, outer] {
        if c.robin.is_some() || c.derivative != vec![0; p.dimension()] {
            return None;
        }
    }
    let (r1, r2, f1, f2) = (inner.point, outer.point, inner.value, outer.value);
    let (sin, cos, pi, defint, legendre) = (
        cx.graph.ops().lookup("sin")?,
        cx.graph.ops().lookup("cos")?,
        cx.graph.ops().lookup("pi")?,
        cx.graph.ops().lookup("defint")?,
        cx.graph.ops().lookup("legendre")?,
    );
    let pi = cx.graph.node(pi, &[]);
    let zero = cx.graph.int(0);
    let two = cx.graph.int(2);
    let two_pi = mul(cx.graph, &[two, pi]);
    let unresolved = |cx: &Cx<'_>, e: NodeId| crate::rules::ode::occurs_op(cx.graph, e, defint);
    let max: i64 = if layout.sphere { 10 } else { 12 };
    let mut terms = Vec::new();
    let mut last = 0;
    // The pair (a, b) for the radial factors `u_a(r)`, `u_b(r)` matching
    // the data (F1, F2) on the two circles.
    let solve_pair = |cx: &mut Cx<'_>,
                      ra: [NodeId; 2],
                      rb: [NodeId; 2],
                      f1n: NodeId,
                      f2n: NodeId|
     -> (NodeId, NodeId) {
        // a ra1 + b rb1 = f1n, a ra2 + b rb2 = f2n.
        let det = {
            let x = mul(cx.graph, &[ra[0], rb[1]]);
            let y = mul(cx.graph, &[rb[0], ra[1]]);
            sub(cx.graph, x, y)
        };
        let a_num = {
            let x = mul(cx.graph, &[f1n, rb[1]]);
            let y = mul(cx.graph, &[f2n, rb[0]]);
            sub(cx.graph, x, y)
        };
        let b_num = {
            let x = mul(cx.graph, &[ra[0], f2n]);
            let y = mul(cx.graph, &[ra[1], f1n]);
            sub(cx.graph, x, y)
        };
        let (a, b) = (div(cx, a_num, det), div(cx, b_num, det));
        (cx.simplify(a), cx.simplify(b))
    };
    for k in 0..=max {
        let n = cx.graph.int(k);
        // Angular functions and their projections.
        let angular: Vec<(NodeId, NodeId, NodeId)> = if layout.sphere {
            let ct = cx.graph.node(cos, &[theta]);
            let st = cx.graph.node(sin, &[theta]);
            let pl = cx.graph.node(legendre, &[n, ct]);
            let scale = {
                let odd = cx.graph.int(2 * k + 1);
                div(cx, odd, two)
            };
            let project = |cx: &mut Cx<'_>, f: NodeId| -> NodeId {
                let body = mul(cx.graph, &[f, pl, st]);
                let integral = cx.graph.node(defint, &[body, theta, zero, pi]);
                let c = mul(cx.graph, &[scale, integral]);
                cx.simplify(c)
            };
            vec![(pl, project(cx, f1), project(cx, f2))]
        } else {
            let waves: Vec<NodeId> = if k == 0 {
                vec![cx.graph.int(1)]
            } else {
                let kt = mul(cx.graph, &[n, theta]);
                vec![cx.graph.node(cos, &[kt]), cx.graph.node(sin, &[kt])]
            };
            let mut out = Vec::new();
            for wave in waves {
                let project = |cx: &mut Cx<'_>, f: NodeId| -> NodeId {
                    let body = mul(cx.graph, &[f, wave]);
                    let integral = cx.graph.node(defint, &[body, theta, zero, two_pi]);
                    let norm = if k == 0 { two_pi } else { pi };
                    let c = div(cx, integral, norm);
                    cx.simplify(c)
                };
                out.push((wave, project(cx, f1), project(cx, f2)));
            }
            out
        };
        for (wave, c1, c2) in angular {
            if unresolved(cx, c1) || unresolved(cx, c2) {
                return None;
            }
            // Radial factors.
            let (u_a, u_b) = if layout.sphere {
                (powi(cx.graph, r, k), powi(cx.graph, r, -k - 1))
            } else if k == 0 {
                (cx.graph.int(1), super::util::call(cx, "ln", &[r])?)
            } else {
                (powi(cx.graph, r, k), powi(cx.graph, r, -k))
            };
            let at = |cx: &mut Cx<'_>, f: NodeId, point: NodeId| cx.graph.substitute(f, r, point);
            let ra = [at(cx, u_a, r1), at(cx, u_a, r2)];
            let rb = [at(cx, u_b, r1), at(cx, u_b, r2)];
            let (a, b) = solve_pair(cx, ra, rb, c1, c2);
            let ua = mul(cx.graph, &[a, u_a]);
            let ub = mul(cx.graph, &[b, u_b]);
            let radial = add(cx.graph, &[ua, ub]);
            let term = mul(cx.graph, &[radial, wave]);
            let term = cx.simplify(term);
            if !cx.is_zero(term) {
                last = k;
            }
            terms.push(term);
        }
    }
    if last >= max {
        return None;
    }
    let total = add(cx.graph, &terms);
    let total = cx.simplify(total);
    super::verified(cx, p, total).then_some(total)
}
