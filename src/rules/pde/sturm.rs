//! Variable-coefficient Sturm–Liouville problems that reduce to the
//! constant-coefficient machinery.
//!
//! * **Euler–Cauchy type.** `a (x-c)² u_xx + b (x-c) u_x + …` reduces to
//!   the constant-coefficient operator `a w_ss + (b - a) w_s` with
//!   `s = ln|x-c|`. The interval lies on one side of the center `c`.
//!   Boundary data transforms by `u_x = w_s / (x-c)`, including Robin
//!   conditions, and the scalar solver supplies the solution for `w`.
//! * **Legendre type.** `κ ((1 - x²) u_xx - 2 x u_x) + c u` on `-1 < x < 1`
//!   with regularity at the ends (no boundary conditions): modes `P_n(x)`,
//!   eigenvalue `n (n + 1)`, polynomial data in closed form.
//! * **Bessel type** — `u_xx + u_x/x`, `u_xx + 2 u_x/x` in any variable
//!   name — is recognised by the curvilinear solver from the structure of
//!   the coefficients.

use super::Condition;
use super::Conditions;
use super::Method;
use super::Problem;
use super::dummy;
use super::spectral::Axis;
use super::spectral::Engine;
use super::spectral::Item;
use super::spectral::Role;
use super::spectral::Shape;
use super::spectral::Symbol;
use super::spectral::Time;
use super::spectral::evolution_time;
use super::spectral::finish;
use super::spectral::time_items;
use super::util::call;
use super::util::div;
use super::util::is_zero_number;
use super::util::jet;
use super::verified;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// Solves the problem by one of the reductions, if it is of that kind.
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    euler(cx, p, conditions)
        .or_else(|| legendre(cx, p, conditions))
        .or_else(|| super::bochner::solve(cx, p, conditions))
        .or_else(|| super::sl::solve(cx, p, conditions))
        .or_else(|| super::hill::solve(cx, p, conditions))
        .or_else(|| super::coulomb::solve(cx, p, conditions))
        .or_else(|| super::spheroidal::solve(cx, p, conditions))
}

/// `c / x^k` if it is free of every variable.
fn power_ratio(
    cx: &mut Cx<'_>,
    p: &Problem,
    c: NodeId,
    x: NodeId,
    k: i64,
) -> Option<NodeId> {
    let inverse = powi(cx.graph, x, -k);
    let q = mul(cx.graph, &[c, inverse]);
    let q = cx.simplify(q);
    p.constant(cx.graph, q).then_some(q)
}

/// Recognizes `a (x-c)^2` and `b (x-c)`, returning `(c, a, b)`.
fn euler_coefficients(
    cx: &mut Cx<'_>,
    p: &Problem,
    x: NodeId,
    second: NodeId,
    first: NodeId,
) -> Option<(NodeId, NodeId, NodeId)> {
    let zero = cx.graph.int(0);
    if let Some(a) = power_ratio(cx, p, second, x, 2) {
        let b = if cx.is_zero(first) {
            zero
        } else {
            power_ratio(cx, p, first, x, 1)?
        };
        return Some((zero, a, b));
    }
    let slope = derivative(cx.graph, second, x)?;
    let curvature = derivative(cx.graph, slope, x)?;
    let curvature = cx.simplify(curvature);
    if !p.constant(cx.graph, curvature) || cx.is_zero(curvature) {
        return None;
    }
    let intercept = cx.graph.substitute(slope, x, zero);
    let intercept = neg(cx.graph, intercept);
    let center = div(cx, intercept, curvature);
    let center = cx.simplify(center);
    let minimum = cx.graph.substitute(second, x, center);
    if !p.constant(cx.graph, center) || !cx.is_zero(minimum) {
        return None;
    }
    let two = cx.graph.int(2);
    let a = div(cx, curvature, two);
    let b = if cx.is_zero(first) {
        zero
    } else {
        let b = derivative(cx.graph, first, x)?;
        let b = cx.simplify(b);
        let offset = cx.graph.substitute(first, x, center);
        if !p.constant(cx.graph, b) || !cx.is_zero(offset) {
            return None;
        }
        b
    };
    Some((center, a, b))
}

#[allow(clippy::too_many_lines)]
fn euler(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let n = p.dimension();
    let time = evolution_time(cx, p);
    // Find the Euler variable and its spatial center.
    let mut chosen = None;
    for j in (0..n).filter(|&j| Some(j) != time) {
        let c2 = p.coefficient(cx.graph, &p.unit(j, 2));
        if cx.is_zero(c2) || p.constant(cx.graph, c2) {
            continue;
        }
        let x = p.vars[j];
        let c1 = p.coefficient(cx.graph, &p.unit(j, 1));
        let (center, a, b) = euler_coefficients(cx, p, x, c2, c1)?;
        chosen = Some((j, center, a, b));
        break;
    }
    let (j, center, a, b) = chosen?;
    let center_is_origin = is_zero_number(cx.graph, center);
    let x = p.vars[j];
    // Every other coefficient is constant and free of derivatives in x.
    for (index, c) in &p.linear {
        if *index == p.unit(j, 2) || *index == p.unit(j, 1) || is_zero_number(cx.graph, *c) {
            continue;
        }
        if index[j] != 0 {
            return None;
        }
        if p.constant(cx.graph, *c) {
            continue;
        }
        // Another Euler variable: `c x_k^m ∂_k^m` with a constant `c`.
        let others: Vec<usize> = (0..n).filter(|&k| index[k] != 0).collect();
        let &[k] = others.as_slice() else {
            return None;
        };
        let order = i64::from(index[k]);
        if Some(k) == time
            || !(1..=2).contains(&order)
            || power_ratio(cx, p, *c, p.vars[k], order).is_none()
        {
            return None;
        }
    }
    // Reuse the signed distances when examining and transforming conditions.
    let distances: Vec<_> = conditions
        .0
        .iter()
        .map(|c| {
            (c.on == j).then(|| {
                if center_is_origin {
                    c.point
                } else {
                    let point = sub(cx.graph, c.point, center);
                    cx.simplify(point)
                }
            })
        })
        .collect();
    let mut negative = false;
    let mut positive = false;
    for &point in distances.iter().flatten() {
        match super::modes::exact_value(cx, point) {
            | Some(v) if v < 0.0 => negative = true,
            | Some(v) if v > 0.0 => positive = true,
            | Some(_) => return None,
            | _ => {},
        }
    }
    if negative && positive {
        return None;
    }
    // The new problem for w(s, ...).
    let (s, _) = dummy(cx, p, "s");
    let w_symbol = cx.graph.interner_mut().fresh_symbol("w");
    let w_head = cx.graph.symbol_node(w_symbol);
    let mut vars = p.vars.clone();
    vars[j] = s;
    let mut args = vec![w_head];
    args.extend(&vars);
    let unknown = cx.graph.node(core::APPLY, &args);
    let mut terms = Vec::new();
    for (index, c) in &p.linear {
        if is_zero_number(cx.graph, *c) {
            continue;
        }
        if *index == p.unit(j, 2) {
            let d = jet(cx, unknown, &vars, index)?;
            terms.push(mul(cx.graph, &[a, d]));
        } else if *index == p.unit(j, 1) {
            // handled with the second derivative below
        } else {
            let d = jet(cx, unknown, &vars, index)?;
            terms.push(mul(cx.graph, &[*c, d]));
        }
    }
    // a (w_ss - w_s) + b w_s.
    let first = {
        let d = jet(cx, unknown, &vars, &p.unit(j, 1))?;
        let gap = sub(cx.graph, b, a);
        mul(cx.graph, &[gap, d])
    };
    terms.push(first);
    let original_x = {
        let e = call(cx, "exp", &[s])?;
        let signed = if negative {
            neg(cx.graph, e)
        } else {
            e
        };
        if center_is_origin {
            signed
        } else {
            add(cx.graph, &[center, signed])
        }
    };
    let source = cx.graph.substitute(p.source, x, original_x);
    terms.push(source);
    let lhs = add(cx.graph, &terms);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[lhs, zero]);
    let scalar = Problem::parse(cx, equation, unknown)?;
    // The conditions.
    let mut list = Vec::new();
    for (c, distance) in conditions.0.iter().zip(distances) {
        let mut value = c.value;
        if let Some(distance) = distance {
            let magnitude = if negative {
                neg(cx.graph, distance)
            } else {
                distance
            };
            let point = call(cx, "ln", &[magnitude])?;
            let point = cx.simplify(point);
            if c.derivative.iter().all(|&d| d == 0) {
                list.push(Condition {
                    on: j,
                    point,
                    derivative: c.derivative.clone(),
                    value,
                    robin: None,
                });
            } else if c.derivative == p.unit(j, 1) {
                // Multiplying u_x + h u = v by the signed distance gives
                // w_s + (x₀-c) h w = (x₀-c) v.
                let scaled = mul(cx.graph, &[distance, value]);
                let scaled = cx.simplify(scaled);
                let robin = match c.robin {
                    | Some(h) => {
                        let hx = mul(cx.graph, &[h, distance]);
                        Some(cx.simplify(hx))
                    },
                    | None => None,
                };
                list.push(Condition {
                    on: j,
                    point,
                    derivative: c.derivative.clone(),
                    value: scaled,
                    robin,
                });
            } else {
                return None;
            }
            continue;
        }
        // Conditions on other variables carry x in their values.
        let uses_unknown = cx
            .graph
            .symbol_of(p.function)
            .is_some_and(|f| cx.graph.depends_on(cx.graph.find(value), f));
        if uses_unknown {
            return None;
        }
        value = cx.graph.substitute(value, x, original_x);
        value = cx.simplify(value);
        list.push(Condition {
            on: c.on,
            point: c.point,
            derivative: c.derivative.clone(),
            value,
            robin: c.robin,
        });
    }
    let found = super::solve(cx, &scalar, &Conditions(list), Method::Any)?;
    let &[_, w] = cx.graph.children(found) else {
        return None;
    };
    let relative = if center_is_origin {
        x
    } else {
        sub(cx.graph, x, center)
    };
    let magnitude = if negative {
        neg(cx.graph, relative)
    } else {
        relative
    };
    let ln_x = call(cx, "ln", &[magnitude])?;
    let solution = cx.graph.substitute(w, s, ln_x);
    let solution = cx.simplify(solution);
    if !super::util::contains_op(cx.graph, solution, "sum")
        && !super::util::contains_op(cx.graph, solution, "defint")
    {
        // For a shifted log chart the generic verification samples in x can
        // lie entirely outside the chart. Verify the exactly transformed PDE.
        let verified = if center_is_origin {
            verified(cx, p, solution)
        } else {
            verified(cx, &scalar, w)
        };
        if !verified {
            return None;
        }
    }
    Some(solution)
}

/// `κ ((1 - x²) u_xx - 2 x u_x) + c u` in one variable with `a u_t` or `a u_tt`.
fn legendre(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let n = p.dimension();
    let time_index = evolution_time(cx, p)?;
    let space: Vec<usize> = (0..n).filter(|&j| j != time_index).collect();
    let &[j] = space.as_slice() else {
        return None;
    };
    let x = p.vars[j];
    let c2 = p.coefficient(cx.graph, &p.unit(j, 2));
    let c1 = p.coefficient(cx.graph, &p.unit(j, 1));
    let one = cx.graph.int(1);
    let x2 = powi(cx.graph, x, 2);
    let width = sub(cx.graph, one, x2);
    let kappa = {
        let q = div(cx, c2, width);
        cx.simplify(q)
    };
    if cx.is_zero(kappa) || !p.constant(cx.graph, kappa) {
        return None;
    }
    let expected = {
        let m = cx.graph.int(-2);
        mul(cx.graph, &[m, x, kappa])
    };
    let gap = sub(cx.graph, c1, expected);
    if !cx.is_zero(gap) {
        return None;
    }
    // Regularity: no conditions on x.
    if conditions.0.iter().any(|c| c.on == j) {
        return None;
    }
    let (a1, a2) = (
        p.coefficient(cx.graph, &p.unit(time_index, 1)),
        p.coefficient(cx.graph, &p.unit(time_index, 2)),
    );
    let c0 = p.coefficient(cx.graph, &vec![0; n]);
    if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) || !p.constant(cx.graph, c0) {
        return None;
    }
    for (index, c) in &p.linear {
        let known = *index == p.unit(j, 1)
            || *index == p.unit(j, 2)
            || *index == vec![0; n]
            || *index == p.unit(time_index, 1)
            || *index == p.unit(time_index, 2);
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let order = if is_zero_number(cx.graph, a2) {
        u32::from(!is_zero_number(cx.graph, a1))
    } else {
        2
    };
    let time = Time {
        var: p.vars[time_index],
        order,
        a2,
        a1,
    };
    // The symbol: κ (-n(n+1)) + c₀.
    let spatial = vec![(vec![1], kappa), (vec![0], c0)];
    let axes = vec![Axis {
        var: x,
        shape: Shape::Legendre,
    }];
    let mut items = time_items(cx, p, conditions, Some(time_index), order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item {
            role: Role::Source,
            h: p.source,
            special: None,
        });
    }
    let engine = Engine {
        axes,
        time: Some(time),
        symbol: Symbol::Cartesian(spatial),
        face: None,
        exterior: None,
    };
    let solution = finish(cx, p, conditions, &engine, &items, None, None)?;
    Some(solution)
}
