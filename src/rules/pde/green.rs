//! Green's functions of linear constant-coefficient problems.
//!
//! `pde_greens_function(equation, u(x, …)[, conditions])` returns the
//! function `G(x, …; s_x, …)` of the field variables and the source
//! variables `s_x, …` (named after the field variables) for the operator
//! `L = Σ c_j ∂_j² + c₀` read off the residual `lhs - rhs` of the equation
//! (a possible source term is ignored):
//!
//! * **elliptic** (no time derivative): `L G = δ(x - s)`;
//! * **first order in time** `a₁ u_t + L u`: the kernel of the initial value
//!   problem, `u(x, 0) = ∫ G(x, s, 0) f(s) ds`, `G(·, 0) = δ`;
//! * **second order without damping** `a₂ u_tt + L u`: the impulse response
//!   (`G(·, 0) = 0`, `G_t(·, 0) = δ`).
//!
//! Two constructions:
//!
//! * **bounded boxes** (two homogeneous Dirichlet, Neumann or Robin
//!   conditions on every axis): the eigenfunction expansion
//!   `Σ Π_j φ_{n_j}(x_j) φ_{n_j}(s_j) / ‖φ_{n_j}‖² × T(σ)` with
//!   `σ = c₀ - Σ c_j λ_{n_j}` and `T = 1/σ`, `e^{-σ t/a₁}` or
//!   `sin(√(σ/a₂) t)/√(σ/a₂)`; zero modes are written out as separate
//!   terms (an elliptic problem with a zero mode and `c₀ = 0` has no
//!   Green's function);
//! * **free space and orthants** (elliptic `L = c Δ`, at most one Dirichlet
//!   or Neumann condition per axis): the fundamental solution
//!   `|x|/2c`, `ln r/(2π c)`, `-1/(4π c r)` with the method of images.

use super::Conditions;
use super::Problem;
use super::spectral::evolution_time;
use super::spectral::interval;
use super::spectral::interval_families;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::fraction;
use super::util::infinity;
use super::util::is_zero_number;
use super::util::pi;
use super::util::series;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::NodeId;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// The source variable belonging to the field variable `x`.
fn source_variable(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(x)?;
    let name = format!("s_{}", cx.graph.interner().symbol_name(symbol));
    Some(cx.graph.sym(&name))
}

/// One way to index the modes of an axis: a literal zero-mode term or the
/// series from `start`.
struct Choice {
    mode_x: NodeId,
    mode_s: NodeId,
    norm: NodeId,
    eigen: NodeId,
    series: Option<(NodeId, i64)>,
}

/// The Green's function, if the problem is of one of the supported kinds.
pub(super) fn greens_function(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let n = p.dimension();
    let time = evolution_time(cx, p);
    let space: Vec<usize> = (0..n).filter(|&j| Some(j) != time).collect();
    if space.is_empty() || space.len() > 3 {
        return None;
    }
    // Coefficients: c_j u_jj, c₀ u, time derivatives; nothing else.
    let mut cjj = Vec::new();
    for &j in &space {
        let c = p.coefficient(cx.graph, &p.unit(j, 2));
        if cx.is_zero(c) || !p.constant(cx.graph, c) {
            return None;
        }
        cjj.push(c);
    }
    let c0 = p.coefficient(cx.graph, &vec![0; n]);
    if !p.constant(cx.graph, c0) {
        return None;
    }
    let (a1, a2) = match time {
        | Some(t) => (p.coefficient(cx.graph, &p.unit(t, 1)), p.coefficient(cx.graph, &p.unit(t, 2))),
        | None => (cx.graph.int(0), cx.graph.int(0)),
    };
    for (index, c) in &p.linear {
        let known = *index == vec![0; n]
            || space.iter().any(|&j| *index == p.unit(j, 2))
            || time.is_some_and(|t| *index == p.unit(t, 1) || *index == p.unit(t, 2));
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let on_axis = |j: usize| conditions.0.iter().filter(|c| c.on == j).count();
    if space.iter().all(|&j| on_axis(j) == 2) {
        return box_expansion(cx, p, conditions, &space, &cjj, c0, time.map(|t| (t, a1, a2)));
    }
    if time.is_none() && space.iter().all(|&j| on_axis(j) <= 1) && cx.is_zero(c0) {
        return fundamental(cx, p, conditions, &space, &cjj);
    }
    None
}

/// The eigenfunction expansion on a box.
#[allow(clippy::too_many_lines)]
fn box_expansion(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    space: &[usize],
    cjj: &[NodeId],
    c0: NodeId,
    time: Option<(usize, NodeId, NodeId)>,
) -> Option<NodeId> {
    let oo = infinity(cx)?;
    let mut axes: Vec<Vec<Choice>> = Vec::new();
    for &j in space {
        let iv = interval(cx, p, conditions, j)?;
        if iv.periodic || !cx.is_zero(iv.left.value) || !cx.is_zero(iv.right.value) {
            return None;
        }
        let x = p.vars[j];
        let s = source_variable(cx, x)?;
        let n_symbol = cx.graph.interner_mut().fresh_symbol("n");
        let n_node = cx.graph.symbol_node(n_symbol);
        cx.graph.assume(n_symbol, Facts::INTEGER | Facts::NONNEGATIVE);
        let general = interval_families(cx, &iv, n_node, false)?;
        let [family] = general.as_slice() else {
            return None;
        };
        let mut choices = Vec::new();
        // Literal low modes (zero or negative eigenvalues).
        for m in family.first..family.first + family.forced {
            let literal = cx.graph.int(m);
            let fam = interval_families(cx, &iv, literal, true)?.into_iter().next()?;
            let norm = match fam.norm {
                | Some(norm) => norm,
                | None => {
                    let square = powi(cx.graph, fam.mode, 2);
                    let integral = defint(cx, square, x, iv.start, iv.stop)?;
                    cx.simplify(integral)
                },
            };
            let mode_s = cx.graph.substitute(fam.mode, x, s);
            choices.push(Choice { mode_x: fam.mode, mode_s, norm, eigen: fam.eigen, series: None });
        }
        let norm = family.norm?;
        let mode_s = cx.graph.substitute(family.mode, x, s);
        choices.push(Choice { mode_x: family.mode, mode_s, norm, eigen: family.eigen, series: Some((n_node, family.first + family.forced)) });
        axes.push(choices);
    }
    // The product of the choices of all axes.
    let mut combinations: Vec<Vec<usize>> = vec![Vec::new()];
    for choices in &axes {
        let mut next = Vec::new();
        for prefix in &combinations {
            for k in 0..choices.len() {
                let mut extended = prefix.clone();
                extended.push(k);
                next.push(extended);
            }
        }
        combinations = next;
    }
    let mut total = Vec::new();
    for combination in combinations {
        let picks: Vec<&Choice> = combination.iter().zip(&axes).map(|(&k, choices)| &choices[k]).collect();
        // σ = c₀ - Σ c_j λ_j.
        let mut sigma_terms = vec![c0];
        for (pick, &c) in picks.iter().zip(cjj) {
            let m = mul(cx.graph, &[c, pick.eigen]);
            sigma_terms.push(neg(cx.graph, m));
        }
        let sigma = add(cx.graph, &sigma_terms);
        let sigma = cx.simplify(sigma);
        let mut factors = Vec::new();
        for pick in &picks {
            let numerator = mul(cx.graph, &[pick.mode_x, pick.mode_s]);
            factors.push(div(cx, numerator, pick.norm));
        }
        let weight = match time {
            | None => {
                if cx.is_zero(sigma) {
                    return None;
                }
                let one = cx.graph.int(1);
                div(cx, one, sigma)
            },
            | Some((t, a1, a2)) => {
                let tv = p.vars[t];
                if cx.is_zero(a2) {
                    let ratio = div(cx, sigma, a1);
                    let e = mul(cx.graph, &[ratio, tv]);
                    let e = neg(cx.graph, e);
                    call(cx, "exp", &[e])?
                } else if cx.is_zero(a1) {
                    if cx.is_zero(sigma) {
                        return None;
                    }
                    let omega2 = div(cx, sigma, a2);
                    let omega = sqrt(cx, omega2)?;
                    let arg = mul(cx.graph, &[omega, tv]);
                    let s = call(cx, "sin", &[arg])?;
                    div(cx, s, omega)
                } else {
                    return None;
                }
            },
        };
        factors.push(weight);
        let mut term = mul(cx.graph, &factors);
        // Wrap the series, innermost axis last.
        for pick in picks.iter().rev() {
            if let Some((index, start)) = pick.series {
                let lo = cx.graph.int(start);
                term = series(cx, term, index, lo, oo)?;
            }
        }
        total.push(term);
    }
    let sum = add(cx.graph, &total);
    Some(cx.simplify(sum))
}

/// Free space and orthants: fundamental solution with images.
fn fundamental(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    space: &[usize],
    cjj: &[NodeId],
) -> Option<NodeId> {
    let d = space.len();
    let c = *cjj.first()?;
    for &other in &cjj[1..] {
        let gap = sub(cx.graph, c, other);
        if !cx.is_zero(gap) {
            return None;
        }
    }
    // Reflections: (axis, point, sign).
    let mut reflections = Vec::new();
    for (pos, &j) in space.iter().enumerate() {
        if let Some(cond) = conditions.0.iter().find(|q| q.on == j) {
            if cond.robin.is_some() || !cx.is_zero(cond.value) {
                return None;
            }
            let sign = if cond.derivative == vec![0; p.dimension()] {
                -1
            } else if cond.derivative == p.unit(j, 1) {
                1
            } else {
                return None;
            };
            reflections.push((pos, cond.point, sign));
        }
    }
    let xs: Vec<NodeId> = space.iter().map(|&j| p.vars[j]).collect();
    let ss: Vec<NodeId> = xs.iter().map(|&x| source_variable(cx, x)).collect::<Option<_>>()?;
    let mut terms = Vec::new();
    for mask in 0..(1_usize << reflections.len()) {
        let mut sources = ss.clone();
        let mut sign = 1_i64;
        for (bit, &(pos, point, s)) in reflections.iter().enumerate() {
            if mask >> bit & 1 == 1 {
                // s ↦ 2 a - s.
                let two = cx.graph.int(2);
                let two_a = mul(cx.graph, &[two, point]);
                sources[pos] = sub(cx.graph, two_a, sources[pos]);
                sign *= s;
            }
        }
        let mut squares = Vec::new();
        for (x, s) in xs.iter().zip(&sources) {
            let gap = sub(cx.graph, *x, *s);
            squares.push(powi(cx.graph, gap, 2));
        }
        let r2 = add(cx.graph, &squares);
        let r = sqrt(cx, r2)?;
        let phi = match d {
            | 1 => {
                let half = fraction(cx, 1, 2)?;
                mul(cx.graph, &[half, r])
            },
            | 2 => {
                let (two, pi_node) = (cx.graph.int(2), pi(cx)?);
                let two_pi = mul(cx.graph, &[two, pi_node]);
                let ln = call(cx, "ln", &[r])?;
                div(cx, ln, two_pi)
            },
            | _ => {
                let (four, pi_node) = (cx.graph.int(4), pi(cx)?);
                let four_pi = mul(cx.graph, &[four, pi_node]);
                let denominator = mul(cx.graph, &[four_pi, r]);
                let minus_one = cx.graph.int(-1);
                div(cx, minus_one, denominator)
            },
        };
        let sign_node = cx.graph.int(sign);
        let signed = mul(cx.graph, &[sign_node, phi]);
        terms.push(signed);
    }
    let sum = add(cx.graph, &terms);
    let g = div(cx, sum, c);
    Some(cx.simplify(g))
}
