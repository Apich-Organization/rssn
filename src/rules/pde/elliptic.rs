//! Poisson, Laplace and Helmholtz equations on half-spaces, quadrants and
//! orthants in two and three dimensions (and half-lines in one).
//!
//! `a (Δu + k² u) + source = 0` on `x_j > 0` for a subset of the variables
//! with a Dirichlet or Neumann condition there: Green's function by the
//! method of images (`G(x - y) ∓ G(x - y*)` over all reflections), with the
//! free-space kernel `G` of `Δ + k²` (`|x|/2`, `ln r / 2π`, `-1/(4π r)`;
//! `e^{ikr}` versions for `k ≠ 0`). For Laplace's equation one face may
//! carry data: the Poisson kernel `Γ(d/2) x_n / π^{d/2} |x - y|^{-d}`
//! for Dirichlet data, the single-layer `2 G` for Neumann data.

use super::Conditions;
use super::Problem;
use super::dummy;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::fraction;
use super::util::imaginary;
use super::util::infinity;
use super::util::is_zero_number;
use super::util::pi;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// The free-space Green's function of `Δ + k²` at distance `r`, with
/// `k = None` for the Laplacian.
pub(super) fn free_kernel(
    cx: &mut Cx<'_>,
    d: usize,
    k: Option<NodeId>,
    r: NodeId,
) -> Option<NodeId> {
    let pi = pi(cx)?;
    let half = fraction(cx, 1, 2)?;
    match (d, k) {
        | (1, None) => Some(mul(cx.graph, &[half, r])),
        | (2, None) => {
            let log = call(cx, "ln", &[r])?;
            let two = cx.graph.int(2);
            let two_pi = mul(cx.graph, &[two, pi]);
            Some(div(cx, log, two_pi))
        },
        | (3, None) => {
            let four = cx.graph.int(4);
            let four_pi_r = mul(cx.graph, &[four, pi, r]);
            let minus_one = cx.graph.int(-1);
            let inverse = powi(cx.graph, four_pi_r, -1);
            Some(mul(cx.graph, &[minus_one, inverse]))
        },
        | (1, Some(k)) => {
            let i = imaginary(cx)?;
            let ikr = mul(cx.graph, &[i, k, r]);
            let wave = call(cx, "exp", &[ikr])?;
            let two = cx.graph.int(2);
            let two_ik = mul(cx.graph, &[two, i, k]);
            Some(div(cx, wave, two_ik))
        },
        | (3, Some(k)) => {
            let i = imaginary(cx)?;
            let ikr = mul(cx.graph, &[i, k, r]);
            let wave = call(cx, "exp", &[ikr])?;
            let four = cx.graph.int(4);
            let four_pi_r = mul(cx.graph, &[four, pi, r]);
            let minus_one = cx.graph.int(-1);
            let inverse = powi(cx.graph, four_pi_r, -1);
            Some(mul(cx.graph, &[minus_one, wave, inverse]))
        },
        | _ => None,
    }
}

/// `(k², f)` for `Δu + k² u = f`, if the equation has that form.
fn helmholtz_form(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<(NodeId, NodeId)> {
    let d = p.dimension();
    let a = p.coefficient(cx.graph, &p.unit(0, 2));
    let mut allowed = vec![vec![0; d]];
    for k in 0..d {
        let index = p.unit(k, 2);
        let ak = p.coefficient(cx.graph, &index);
        let gap = sub(cx.graph, ak, a);
        if !cx.is_zero(gap) {
            return None;
        }
        allowed.push(index);
    }
    if !p.only(cx.graph, &allowed) || p.nonlinear || !p.constant(cx.graph, a) || cx.is_zero(a) {
        return None;
    }
    let inverse_a = powi(cx.graph, a, -1);
    let c = p.coefficient(cx.graph, &vec![0; d]);
    let k2 = mul(cx.graph, &[c, inverse_a]);
    let k2 = cx.simplify(k2);
    if !p.constant(cx.graph, k2) {
        return None;
    }
    let minus_one = cx.graph.int(-1);
    let f = mul(cx.graph, &[minus_one, p.source, inverse_a]);
    let f = cx.simplify(f);
    Some((k2, f))
}

/// Half-space problems with images.
#[allow(clippy::too_many_lines)]
pub(super) fn half_space(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let d = p.dimension();
    if !(1..=3).contains(&d) || conditions.is_empty() {
        return None;
    }
    let (k2, f) = helmholtz_form(cx, p)?;
    let helmholtz = !is_zero_number(cx.graph, k2);
    let zero_index = vec![0; d];
    // Domains: sign of the image per variable, boundary data.
    let mut parity = Vec::new();
    let mut data: Option<(usize, bool, NodeId)> = None;
    for j in 0..d {
        let on: Vec<_> = conditions.0.iter().filter(|c| c.on == j).collect();
        match on.as_slice() {
            | [] => parity.push(0_i64),
            | [b] if b.robin.is_none() && is_zero_number(cx.graph, b.point) => {
                // A value that refers to the unknown itself is not data.
                if cx.graph.depends_on(cx.graph.find(b.value), cx.graph.symbol_of(p.function)?) {
                    return None;
                }
                let dirichlet = b.derivative == zero_index;
                if !dirichlet && b.derivative != p.unit(j, 1) {
                    return None;
                }
                parity.push(if dirichlet { -1 } else { 1 });
                if !cx.is_zero(b.value) {
                    if data.is_some() {
                        return None;
                    }
                    data = Some((j, dirichlet, b.value));
                }
            },
            | _ => return None,
        }
    }
    if conditions.0.len() != parity.iter().filter(|&&q| q != 0).count() || parity.iter().all(|&q| q == 0) {
        return None;
    }
    if data.is_some() && (helmholtz || d < 2) {
        return None;
    }
    if cx.is_zero(f) && data.is_none() && !helmholtz {
        // Laplace's equation with homogeneous data on a half-space: zero.
        return None;
    }
    let k = if helmholtz { Some(sqrt(cx, k2)?) } else { None };
    let oo = infinity(cx)?;
    let minus_oo = neg(cx.graph, oo);
    let stems = ["s", "r", "q"];
    let mut dummies = Vec::new();
    for stem in stems.iter().take(d) {
        dummies.push(dummy(cx, p, stem).0);
    }
    // All reflections: per variable, the choices (sign, position).
    let mut combos: Vec<(i64, Vec<NodeId>)> = vec![(1, Vec::new())];
    for j in 0..d {
        let mut next = Vec::new();
        for (sign, positions) in &combos {
            let mut straight = positions.clone();
            straight.push(dummies[j]);
            next.push((*sign, straight));
            if parity[j] != 0 {
                let mut reflected = positions.clone();
                reflected.push(neg(cx.graph, dummies[j]));
                next.push((*sign * parity[j], reflected));
            }
        }
        combos = next;
    }
    let distance = |cx: &mut Cx<'_>, positions: &[NodeId], skip: Option<usize>| -> Option<NodeId> {
        let mut squares = Vec::new();
        for (j, &y) in positions.iter().enumerate() {
            let x = p.vars[j];
            let gap = if Some(j) == skip { x } else { sub(cx.graph, x, y) };
            squares.push(powi(cx.graph, gap, 2));
        }
        let r2 = add(cx.graph, &squares);
        sqrt(cx, r2)
    };
    let domain_lower = |cx: &mut Cx<'_>, j: usize| if parity[j] == 0 { minus_oo } else { cx.graph.int(0) };
    let mut parts = Vec::new();
    // Source: images of the free kernel.
    if !cx.is_zero(f) {
        let mut kernel_sum = Vec::new();
        for (sign, positions) in &combos {
            let r = distance(cx, positions, None)?;
            let g = free_kernel(cx, d, k, r)?;
            let s = cx.graph.int(*sign);
            kernel_sum.push(mul(cx.graph, &[s, g]));
        }
        let kernel = add(cx.graph, &kernel_sum);
        let mut source = f;
        for (j, &s) in dummies.iter().enumerate() {
            source = cx.graph.substitute(source, p.vars[j], s);
        }
        let mut body = mul(cx.graph, &[kernel, source]);
        for j in (0..d).rev() {
            let lower = domain_lower(cx, j);
            body = defint(cx, body, dummies[j], lower, oo)?;
            body = cx.simplify(body);
        }
        parts.push(body);
    }
    // Boundary data: Poisson kernel or single layer on the face x_n = 0.
    if let Some((axis, dirichlet, g)) = data {
        let mut g_at = g;
        for (j, &s) in dummies.iter().enumerate() {
            if j != axis {
                g_at = cx.graph.substitute(g_at, p.vars[j], s);
            }
        }
        let mut kernel_sum = Vec::new();
        for (sign, positions) in &combos {
            // The face variable is not integrated: position zero.
            if positions.get(axis) != Some(&dummies[axis]) {
                continue;
            }
            let r = distance(cx, positions, Some(axis))?;
            let weight = if dirichlet {
                // Γ(d/2)/π^(d/2) x_n r^(-d)
                let pi = pi(cx)?;
                let c = if d == 2 {
                    powi(cx.graph, pi, -1)
                } else {
                    let two = cx.graph.int(2);
                    let two_pi = mul(cx.graph, &[two, pi]);
                    powi(cx.graph, two_pi, -1)
                };
                let rd = powi(cx.graph, r, -i64::try_from(d).ok()?);
                mul(cx.graph, &[c, p.vars[axis], rd])
            } else {
                let g0 = free_kernel(cx, d, None, r)?;
                let two = cx.graph.int(2);
                mul(cx.graph, &[two, g0])
            };
            let s = cx.graph.int(*sign);
            kernel_sum.push(mul(cx.graph, &[s, weight]));
        }
        let kernel = add(cx.graph, &kernel_sum);
        let mut body = mul(cx.graph, &[kernel, g_at]);
        for j in (0..d).rev() {
            if j == axis {
                continue;
            }
            let lower = domain_lower(cx, j);
            body = defint(cx, body, dummies[j], lower, oo)?;
            body = cx.simplify(body);
        }
        parts.push(body);
    }
    let total = add(cx.graph, &parts);
    Some(cx.simplify(total))
}
