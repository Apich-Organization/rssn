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

use super::Conditions;
use super::Problem;
use super::classical;
use super::spectral::Axis;
use super::spectral::End;
use super::spectral::Engine;
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
use crate::graph::NodeId;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

#[derive(Copy, Clone, Debug)]
struct Layout {
    r: usize,
    theta: Option<usize>,
    z: Option<usize>,
    sphere: bool,
}

fn name_of(
    cx: &Cx<'_>,
    node: NodeId,
) -> Option<String> {
    cx.graph.symbol_of(node).map(|s| cx.graph.interner().symbol_name(s).to_owned())
}

/// The layout of the variables, from their names.
fn layouts(
    cx: &Cx<'_>,
    p: &Problem,
    space: &[usize],
) -> Vec<Layout> {
    let names: Vec<(usize, String)> = space.iter().filter_map(|&j| name_of(cx, p.vars[j]).map(|n| (j, n))).collect();
    let find = |wanted: &[&str]| names.iter().find(|(_, n)| wanted.contains(&n.as_str())).map(|(j, _)| *j);
    let (Some(r), theta, z) = (find(&["r", "rho"]), find(&["theta", "th", "phi", "varphi", "vartheta"]), find(&["z"])) else {
        return Vec::new();
    };
    if space.len() != 1 + usize::from(theta.is_some()) + usize::from(z.is_some()) {
        return Vec::new();
    }
    let mut out = vec![Layout { r, theta, z, sphere: false }];
    if z.is_none() {
        out.push(Layout { r, theta, z, sphere: true });
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
    if let Some(j) = layout.z {
        expected.push((p.unit(j, 2), cx.graph.int(1)));
    }
    for (index, c) in &p.linear {
        let is_time = time.is_some_and(|t| index[t] > 0);
        if is_time || *index == vec![0; n] || is_zero_number(cx.graph, *c) {
            continue;
        }
        let Some((_, e)) = expected.iter().find(|(i, _)| i == index) else {
            return None;
        };
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
    // The condition at the outer radius; nothing else on r, θ.
    let on_r: Vec<_> = conditions.0.iter().filter(|c| c.on == layout.r).collect();
    let [outer] = on_r.as_slice() else {
        return None;
    };
    if conditions.0.iter().any(|c| Some(c.on) == layout.theta) {
        return None;
    }
    let radius = outer.point;
    let end = End::from_condition(cx, p, outer)?;
    // Harmonic functions with Dirichlet data: the closed series.
    if time_index.is_none() && layout.z.is_none() && p.homogeneous(cx.graph) && is_zero_number(cx.graph, a0) && end.kind == Kind::Dirichlet && cx.graph.number_of(end.p).is_some() {
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
    // Time structure.
    let (mut a2, mut a1) = (cx.graph.int(0), cx.graph.int(0));
    if let Some(t) = time_index {
        a1 = p.coefficient(cx.graph, &p.unit(t, 1));
        a2 = p.coefficient(cx.graph, &p.unit(t, 2));
        if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
            return None;
        }
    }
    let order = if !is_zero_number(cx.graph, a2) {
        2
    } else if !is_zero_number(cx.graph, a1) {
        1
    } else {
        0
    };
    let time = time_index.map(|t| Time { var: p.vars[t], order, a2, a1 });
    // Axes in the order angle, radius, axis.
    let mut axes = Vec::new();
    if let Some(j) = layout.theta {
        axes.push(Axis { var: p.vars[j], shape: if layout.sphere { Shape::Polar } else { Shape::Angle } });
    }
    let radial_level = axes.len();
    axes.push(Axis {
        var: p.vars[layout.r],
        shape: Shape::Radial { radius, end, sphere: layout.sphere, has_angle: layout.theta.is_some() },
    });
    let mut z_interval = None;
    if let Some(j) = layout.z {
        let iv = interval(cx, p, conditions, j)?;
        axes.push(Axis { var: p.vars[j], shape: Shape::Interval(iv.clone()) });
        z_interval = Some((axes.len() - 1, iv));
    }
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
    // Boundary data at r = R: Green's identity with weight R^(d-1).
    if !cx.is_zero(end.value) {
        let divisor = if end.kind == Kind::Dirichlet { end.p } else { end.q };
        let scaled = mul(cx.graph, &[b, end.value]);
        let scaled = div(cx, scaled, divisor);
        let scaled = cx.simplify(scaled);
        let weight = powi(cx.graph, radius, if layout.sphere { 2 } else { 1 });
        items.push(Item {
            role: Role::Source,
            h: scaled,
            special: Some(Special { axis: radial_level, kind: end.kind, point: radius, weight }),
        });
    }
    if let Some((level, iv)) = &z_interval {
        items.extend(boundary_items(cx, iv, *level, b)?);
    }
    let engine = Engine { axes, time, symbol: Symbol::Laplacian { b, a0 } };
    let _ = derivative;
    finish(cx, p, conditions, &engine, &items)
}
