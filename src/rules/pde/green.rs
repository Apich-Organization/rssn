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
        | Some(t) => {
            (
                p.coefficient(cx.graph, &p.unit(t, 1)),
                p.coefficient(cx.graph, &p.unit(t, 2)),
            )
        },
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
    // A periodic axis carries the single condition `u(a) = u(b)`.
    let periodic = |cx: &Cx<'_>, j: usize| {
        conditions.0.iter().any(|c| {
            c.on == j
                && cx.graph.op(c.value) == crate::graph::op::core::APPLY
                && cx.graph.children(c.value).first() == Some(&p.function)
        })
    };
    if space
        .iter()
        .all(|&j| on_axis(j) == 2 || (on_axis(j) == 1 && periodic(cx, j)))
    {
        return box_expansion(
            cx,
            p,
            conditions,
            &space,
            &cjj,
            c0,
            time.map(|t| (t, a1, a2)),
        );
    }
    if time.is_none() && space.iter().all(|&j| on_axis(j) <= 1) && cx.is_zero(c0) {
        return fundamental(cx, p, conditions, &space, &cjj);
    }
    if let Some(t) = time
        && cx.is_zero(a2)
        && space.iter().all(|&j| on_axis(j) <= 1)
    {
        return heat_kernel(cx, p, conditions, &space, &cjj, c0, (t, a1));
    }
    if let Some(t) = time
        && cx.is_zero(a1)
        && !cx.is_zero(a2)
        && cx.is_zero(c0)
        && space.iter().all(|&j| on_axis(j) <= 1)
    {
        return wave_kernel(cx, p, conditions, &space, &cjj, (t, a2));
    }
    None
}

/// The impulse response (`G = 0`, `G_t = δ` at `t = 0`) of
/// `a₂ u_tt + c Δ u` on whole space and orthants, with `c² = -c/a₂`:
///
/// * one dimension: `H(ct - |x - s|)/(2c)`;
/// * two: `H(ct - r)/(2π c √(c²t² - r²))`;
/// * three: `δ(ct - r)/(4π c r)`;
///
/// plus the images `±E(x - s')` for a Dirichlet or Neumann condition per
/// axis.
fn wave_kernel(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    space: &[usize],
    cjj: &[NodeId],
    (time, a2): (usize, NodeId),
) -> Option<NodeId> {
    let d = space.len();
    let c0_value = *cjj.first()?;
    for &other in &cjj[1..] {
        let gap = sub(cx.graph, c0_value, other);
        if !cx.is_zero(gap) {
            return None;
        }
    }
    let speed2 = {
        let q = div(cx, c0_value, a2);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    if super::util::sample(cx.graph, speed2, 0)? <= 0.0 {
        return None;
    }
    let c = sqrt(cx, speed2)?;
    let t = p.vars[time];
    let ct = mul(cx.graph, &[c, t]);
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
    let ss: Vec<NodeId> = xs
        .iter()
        .map(|&x| source_variable(cx, x))
        .collect::<Option<_>>()?;
    let mut terms = Vec::new();
    for mask in 0..(1_usize << reflections.len()) {
        let mut sources = ss.clone();
        let mut sign = 1_i64;
        for (bit, &(pos, point, s)) in reflections.iter().enumerate() {
            if mask >> bit & 1 == 1 {
                let two = cx.graph.int(2);
                let two_a = mul(cx.graph, &[two, point]);
                sources[pos] = sub(cx.graph, two_a, sources[pos]);
                sign *= s;
            }
        }
        let squares: Vec<NodeId> = xs
            .iter()
            .zip(&sources)
            .map(|(x, s)| {
                let gap = sub(cx.graph, *x, *s);
                powi(cx.graph, gap, 2)
            })
            .collect();
        let r2 = add(cx.graph, &squares);
        let r = sqrt(cx, r2)?;
        let reach = sub(cx.graph, ct, r);
        let kernel = match d {
            | 1 => {
                let h = call(cx, "heaviside", &[reach])?;
                let two = cx.graph.int(2);
                let denominator = mul(cx.graph, &[two, c]);
                div(cx, h, denominator)
            },
            | 2 => {
                let h = call(cx, "heaviside", &[reach])?;
                let (two, pi_node) = (cx.graph.int(2), pi(cx)?);
                let ct2 = powi(cx.graph, ct, 2);
                let gap = sub(cx.graph, ct2, r2);
                let root = sqrt(cx, gap)?;
                let denominator = mul(cx.graph, &[two, pi_node, c, root]);
                div(cx, h, denominator)
            },
            | _ => {
                let delta = call(cx, "dirac", &[reach])?;
                let (four, pi_node) = (cx.graph.int(4), pi(cx)?);
                let denominator = mul(cx.graph, &[four, pi_node, c, r]);
                div(cx, delta, denominator)
            },
        };
        let sign_node = cx.graph.int(sign);
        terms.push(mul(cx.graph, &[sign_node, kernel]));
    }
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

/// The heat kernel of `a₁ u_t + Σ c_j u_jj + c₀ u` on whole space and
/// orthants (images for a Dirichlet or Neumann condition per axis):
/// `e^{-c₀ t/a₁} Π_j [K_j(x_j - s_j) + σ_j K_j(x_j - (2 a_j - s_j))]`,
/// `K_j` the Gaussian of diffusivity `D_j = -c_j/a₁`.
fn heat_kernel(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    space: &[usize],
    cjj: &[NodeId],
    c0: NodeId,
    (time, a1): (usize, NodeId),
) -> Option<NodeId> {
    let t = p.vars[time];
    let mut factors = Vec::new();
    for (&j, &c) in space.iter().zip(cjj) {
        let d = {
            let q = div(cx, c, a1);
            let q = neg(cx.graph, q);
            cx.simplify(q)
        };
        if super::util::sample(cx.graph, d, 0)? <= 0.0 {
            return None;
        }
        let x = p.vars[j];
        let s = source_variable(cx, x)?;
        let direct = sub(cx.graph, x, s);
        let mut axis = super::classical::kernel(cx, d, direct, t)?;
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
            let two = cx.graph.int(2);
            let two_a = mul(cx.graph, &[two, cond.point]);
            let reflected = sub(cx.graph, two_a, s);
            let gap = sub(cx.graph, x, reflected);
            let image = super::classical::kernel(cx, d, gap, t)?;
            let signed = {
                let sign_node = cx.graph.int(sign);
                mul(cx.graph, &[sign_node, image])
            };
            axis = add(cx.graph, &[axis, signed]);
        }
        factors.push(axis);
    }
    if !cx.is_zero(c0) {
        let ratio = div(cx, c0, a1);
        let e = mul(cx.graph, &[ratio, t]);
        let e = neg(cx.graph, e);
        factors.push(call(cx, "exp", &[e])?);
    }
    let g = mul(cx.graph, &factors);
    Some(cx.simplify(g))
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
    for (position, &j) in space.iter().enumerate() {
        let iv = interval(cx, p, conditions, j)?;
        if !cx.is_zero(iv.left.value) || !cx.is_zero(iv.right.value) {
            return None;
        }
        let x = p.vars[j];
        let s = source_variable(cx, x)?;
        let stems = ["n", "m", "l"];
        let (n_node, n_symbol) = super::dummy(cx, p, stems.get(position).copied().unwrap_or("q"));
        cx.graph
            .assume(n_symbol, Facts::INTEGER | Facts::NONNEGATIVE);
        let general = interval_families(cx, &iv, n_node, false)?;
        let mut choices = Vec::new();
        for (family_index, family) in general.iter().enumerate() {
            // Literal low modes (zero or negative eigenvalues).
            for m in family.first..family.first + family.forced {
                let literal = cx.graph.int(m);
                let fam = interval_families(cx, &iv, literal, true)?
                    .into_iter()
                    .nth(family_index)?;
                let norm = match fam.norm {
                    | Some(norm) => norm,
                    | None => {
                        let square = powi(cx.graph, fam.mode, 2);
                        let integral = defint(cx, square, x, iv.start, iv.stop)?;
                        cx.simplify(integral)
                    },
                };
                let mode_s = cx.graph.substitute(fam.mode, x, s);
                choices.push(Choice {
                    mode_x: fam.mode,
                    mode_s,
                    norm,
                    eigen: fam.eigen,
                    series: None,
                });
            }
            let norm = family.norm?;
            let mode_s = cx.graph.substitute(family.mode, x, s);
            choices.push(Choice {
                mode_x: family.mode,
                mode_s,
                norm,
                eigen: family.eigen,
                series: Some((n_node, family.first + family.forced)),
            });
        }
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
        let picks: Vec<&Choice> = combination
            .iter()
            .zip(&axes)
            .map(|(&k, choices)| &choices[k])
            .collect();
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
                    // a₂ T'' + a₁ T' + σ T = 0, T(0) = 0, T'(0) = 1:
                    // e^{-α t} sinh(β t)/(a₂ β), α = a₁/(2 a₂), β² = α² - σ/a₂.
                    let half = fraction(cx, 1, 2)?;
                    let ratio = div(cx, a1, a2);
                    let alpha = mul(cx.graph, &[half, ratio]);
                    let alpha = cx.simplify(alpha);
                    let beta2 = {
                        let a_sq = powi(cx.graph, alpha, 2);
                        let q = div(cx, sigma, a2);
                        let nq = neg(cx.graph, q);
                        add(cx.graph, &[a_sq, nq])
                    };
                    let beta2 = cx.simplify(beta2);
                    let beta = sqrt(cx, beta2)?;
                    let bt = mul(cx.graph, &[beta, tv]);
                    let sinh = call(cx, "sinh", &[bt])?;
                    let damp = {
                        let m = mul(cx.graph, &[alpha, tv]);
                        let nm = neg(cx.graph, m);
                        call(cx, "exp", &[nm])?
                    };
                    let denominator = mul(cx.graph, &[a2, beta]);
                    let q = div(cx, sinh, denominator);
                    mul(cx.graph, &[damp, q])
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
    let ss: Vec<NodeId> = xs
        .iter()
        .map(|&x| source_variable(cx, x))
        .collect::<Option<_>>()?;
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
