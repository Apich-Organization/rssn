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
use super::util::div;
use super::util::is_zero_number;
use crate::graph::Cx;
use crate::graph::op::core;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::mul;
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
    cx.graph.symbol_of(node).map(|s| cx.graph.interner().symbol_name(s).to_owned())
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
        let (b, c1) = (p.coefficient(cx.graph, &p.unit(j, 2)), p.coefficient(cx.graph, &p.unit(j, 1)));
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
    let names: Vec<(usize, String)> = space.iter().filter_map(|&j| name_of(cx, p.vars[j]).map(|n| (j, n))).collect();
    let find = |wanted: &[&str]| names.iter().find(|(_, n)| wanted.contains(&n.as_str())).map(|(j, _)| *j);
    let z = find(&["z"]);
    let r = find(&["r", "rho"]).or_else(|| {
        let candidates: Vec<usize> = space.iter().copied().filter(|&j| Some(j) != z).collect();
        if candidates.len() == 1 { radial_by_structure(cx, p, &candidates) } else { None }
    });
    let Some(r) = r else {
        return Vec::new();
    };
    let others: Vec<usize> = space.iter().copied().filter(|&j| j != r && Some(j) != z).collect();
    let is = |j: usize, list: &[&str]| names.iter().any(|(k, n)| *k == j && list.contains(&n.as_str()));
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
        out.push(Layout { r, theta, phi, z, sphere: false });
    }
    if z.is_none() {
        out.push(Layout { r, theta, phi, z, sphere: true });
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
    let space: Vec<usize> = (0..p.dimension()).filter(|&j| Some(j) != time_index).collect();
    for layout in layouts(cx, p, &space) {
        let Some((b, a0)) = laplacian(cx, p, layout, time_index) else {
            continue;
        };
        if let Some(solution) = solve_layout(cx, p, conditions, layout, time_index, b, a0) {
            return Some(solution);
        }
    }
    None
}

/// Whether the node is `oo`.
fn is_infinity(
    cx: &Cx<'_>,
    node: NodeId,
) -> bool {
    cx.graph.ops().lookup("oo").is_some_and(|op| cx.graph.op(node) == op)
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
    let on_r: Vec<&Condition> = conditions.0.iter().filter(|c| c.on == layout.r).collect();
    if conditions.0.iter().any(|c| Some(c.on) == layout.theta || Some(c.on) == layout.phi) {
        return None;
    }
    let (at_infinity, finite): (Vec<&Condition>, Vec<&Condition>) = on_r.into_iter().partition(|c| is_infinity(cx, c.point));
    if let [far] = at_infinity.as_slice() {
        return exterior(cx, p, conditions, layout, time_index, (b, a0), &finite, far);
    }
    if !at_infinity.is_empty() {
        return None;
    }
    let (inner, outer) = match finite.as_slice() {
        | [first, second] => {
            if let Some(found) = annulus(cx, p, conditions, layout, time_index, (b, a0), [*first, *second]) {
                return Some(found);
            }
            let (va, vb) = (super::util::sample(cx.graph, first.point, 0)?, super::util::sample(cx.graph, second.point, 0)?);
            if va <= vb { (Some(*first), *second) } else { (Some(*second), *first) }
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
    if inner.is_none() && time_index.is_none() && layout.z.is_none() && layout.phi.is_none() && p.homogeneous(cx.graph) && is_zero_number(cx.graph, a0) && end.kind == Kind::Dirichlet && cx.graph.number_of(end.p).is_some() {
        if let Some(theta) = layout.theta {
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
    }
    if inner.is_none() && time_index.is_none() && layout.z.is_none() && layout.phi.is_none() && !p.homogeneous(cx.graph) && is_zero_number(cx.graph, a0) && end.kind == Kind::Dirichlet {
        if let Some(found) = radial_poisson(cx, p, layout, b, end, radius) {
            return Some(found);
        }
    }
    // Harmonic functions with Robin or Neumann data, or without axial
    // symmetry: each harmonic times its power of r.
    if inner.is_none() && time_index.is_none() && layout.z.is_none() && p.homogeneous(cx.graph) && cx.is_zero(a0) && layout.theta.is_some() {
        if let Some(found) = harmonic_interior(cx, p, conditions, layout, (b, a0), end, radius) {
            return Some(found);
        }
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
    let order = if is_zero_number(cx.graph, a2) { u32::from(!is_zero_number(cx.graph, a1)) } else { 2 };
    let time = time_index.map(|t| Time { var: p.vars[t], order, a2, a1 });
    // Axes in the order azimuth, angle, radius, axis.
    let mut axes = Vec::new();
    if let Some(j) = layout.phi {
        axes.push(Axis { var: p.vars[j], shape: Shape::Angle });
    }
    if let Some(j) = layout.theta {
        axes.push(Axis { var: p.vars[j], shape: if layout.sphere { Shape::Polar { azimuth: layout.phi.is_some() } } else { Shape::Angle } });
    }
    let radial_level = axes.len();
    axes.push(Axis {
        var: p.vars[layout.r],
        shape: Shape::Radial { radius, inner: inner_end, end, sphere: layout.sphere, has_angle: layout.theta.is_some() },
    });
    let z_interval = match layout.z {
        | Some(j) => {
            let iv = interval(cx, p, conditions, j)?;
            axes.push(Axis { var: p.vars[j], shape: Shape::Interval(iv.clone()) });
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
        items.push(Item { role: Role::Source, h: p.source, special: None });
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
        let divisor = if e.kind == Kind::Dirichlet { e.p } else { e.q };
        let sign = cx.graph.int(side);
        let scaled = mul(cx.graph, &[sign, b, e.value]);
        let scaled = div(cx, scaled, divisor);
        let scaled = cx.simplify(scaled);
        let weight = powi(cx.graph, point, if layout.sphere { 2 } else { 1 });
        items.push(Item { role: Role::Source, h: scaled, special: Some(Special { axis: radial_level, kind: e.kind, point, weight }) });
    }
    if let Some((level, iv)) = &z_interval {
        items.extend(boundary_items(cx, iv, *level, b));
    }
    let engine = Engine { axes, time, symbol: Symbol::Laplacian { b, a0 }, face: None, exterior: None };
    finish(cx, p, conditions, &engine, &items, None, None)
}

/// Laplace's and Helmholtz's equations outside a disk or a ball (and the
/// heat and wave equations outside a ball for radial data): the condition
/// at `r = oo` selects the bounded decaying solution (`u → c`) or, for the
/// Helmholtz equation `Δu + k² u = 0` with `k² > 0`, the outgoing one
/// (Sommerfeld radiation condition, time factor `e^{-iωt}`: Hankel functions
/// of the first kind); the angular dependence is expanded as in the interior
/// and each harmonic is multiplied by its decaying radial profile.
#[allow(clippy::too_many_lines)]
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
    let c_inf = if laplace { far.value } else { cx.graph.int(0) };
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
        axes.push(Axis { var: p.vars[j], shape: Shape::Angle });
    }
    if let Some(j) = layout.theta {
        axes.push(Axis { var: p.vars[j], shape: if layout.sphere { Shape::Polar { azimuth: layout.phi.is_some() } } else { Shape::Angle } });
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
        let cleared = cx.is_zero(mean) || (0..3_u32).all(|k| super::util::sample(cx.graph, mean, k).is_some_and(|v| v.abs() < 1e-9));
        if !cleared {
            return None;
        }
    }
    let exterior = Exterior { var: p.vars[layout.r], radius, end, sphere: layout.sphere, interior: false, wavenumber2 };
    let engine = Engine { axes, time: None, symbol: Symbol::Laplacian { b, a0 }, face: None, exterior: Some(exterior) };
    let items = [Item { role: Role::Displacement, h: data, special: None }];
    let kept = Conditions(conditions.0.iter().filter(|c| c.on != layout.r || !is_infinity(cx, c.point)).cloned().collect());
    let offset = if cx.is_zero(c_inf) { None } else { Some(c_inf) };
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
        axes.push(Axis { var: p.vars[j], shape: Shape::Angle });
    }
    let j = layout.theta?;
    axes.push(Axis { var: p.vars[j], shape: if layout.sphere { Shape::Polar { azimuth: layout.phi.is_some() } } else { Shape::Angle } });
    let wavenumber2 = cx.graph.int(0);
    let interior = Exterior { var: p.vars[layout.r], radius, end, sphere: layout.sphere, interior: true, wavenumber2 };
    let engine = Engine { axes, time: None, symbol: Symbol::Laplacian { b, a0 }, face: None, exterior: Some(interior) };
    let items = [Item { role: Role::Displacement, h: end.value, special: None }];
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
    if !layout.sphere || layout.theta.is_some() || !p.homogeneous(cx.graph) || !cx.is_zero(a0) || !cx.is_zero(far.value) || end.kind != Kind::Dirichlet {
        return None;
    }
    let n = p.dimension();
    let c_second = p.coefficient(cx.graph, &p.unit(t_index, 2));
    let order = if is_zero_number(cx.graph, c_second) { 1 } else { 2 };
    let c_time = p.coefficient(cx.graph, &p.unit(t_index, order));
    let diffusivity = {
        let q = div(cx, b, c_time);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let (r, t) = (p.vars[layout.r], p.vars[t_index]);
    let radius = end.point;
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
    list.push(Condition { on: 0, point: zero, derivative: vec![0, 0], value: boundary, robin: None });
    let zero_index = vec![0; n];
    for (derivative, vector) in [(zero_index.clone(), vec![0, 0]), (p.unit(t_index, 1), vec![0, 1])] {
        if let Some(c) = conditions.find(cx.graph, t_index, &derivative, None) {
            let f = cx.graph.substitute(c.value, r, shifted);
            let value = mul(cx.graph, &[shifted, f]);
            let value = cx.simplify(value);
            list.push(Condition { on: 1, point: zero, derivative: vector, value, robin: None });
        }
    }
    if conditions.0.iter().any(|c| c.on == t_index && c.derivative != zero_index && c.derivative != p.unit(t_index, 1)) {
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
        let depends = cx.graph.symbol_of(v).is_some_and(|s| cx.graph.depends_on(cx.graph.find(f), s));
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
        let denominator = if layout.sphere { (k + 2) * (k + 3) } else { (k + 2) * (k + 2) };
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
            if cx.graph.depends_on(cx.graph.find(data), cx.graph.symbol_of(r)?) {
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
    if time_index.is_some() || layout.z.is_some() || !p.homogeneous(cx.graph) || !is_zero_number(cx.graph, a0) || conditions.0.len() != 2 {
        return None;
    }
    let r = p.vars[layout.r];
    let [first, second] = ends;
    let (inner, outer) = {
        let (va, vb) = (super::util::sample(cx.graph, first.point, 0)?, super::util::sample(cx.graph, second.point, 0)?);
        if va <= vb { (first, second) } else { (second, first) }
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
    let solve_pair = |cx: &mut Cx<'_>, ra: [NodeId; 2], rb: [NodeId; 2], f1n: NodeId, f2n: NodeId| -> (NodeId, NodeId) {
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
