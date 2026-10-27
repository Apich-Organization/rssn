//! Radical and trigonometric solutions of irreducible sextics and septics
//! from their Galois group.
//!
//! * **Cyclic groups** (`C6`, `C7`) whose roots are real and whose
//!   conductor is a prime `p`: the roots are rational combinations of
//!   `2 cos(2π j / p)` (Gauss periods), see
//!   [`crate::rules::poly::galois_high::cyclotomic`].
//! * **Solvable sextics.** Every solvable transitive group of degree 6 is
//!   imprimitive, with blocks of size 2 or 3. The quotient polynomial `Θ` (of
//!   degree 3 or 2) is solved in radicals, and the block polynomial, with
//!   coefficients that are rational polynomials in a root `θ` of `Θ`, is
//!   solved over it by the quadratic or cubic formula
//!   ([`crate::rules::poly::galois_high::towers`]). Only real values of `θ`
//!   matter for real roots (a block containing a real root is invariant under
//!   conjugation).
//!
//! Every candidate is checked against the polynomial and matched with the
//! real roots, which must all be found, else the caller keeps its `root_of`
//! objects. The unsolvable groups (`A5`, `S5`, `A6`, `S6`, `PSL(3,2)`, `A7`,
//! `S7`) and the non-abelian septic groups `D7`, `F21`, `F42` (apart from
//! pure powers, which are solved earlier) stay `root_of`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::galois_high::analyse;
use crate::rules::poly::galois_high::cyclotomic;
use crate::rules::poly::galois_high::cyclotomic_complex;
use crate::rules::poly::galois_high::lagrange;
use crate::rules::poly::galois_high::Lagrange;
use crate::rules::poly::galois_high::towers;
use crate::rules::poly::galois_high::Analysis;
use crate::rules::poly::galois_high::Tower;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Poly;

fn value_of(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.eval(node, &Env::numeric(0.0)).filter(|v| v.is_finite())
}

/// Relative residual test of the integer polynomial at `x`.
fn is_root(
    factor: &[BigInt],
    x: f64,
) -> bool {
    let mut value = 0.0_f64;
    let mut scale = 0.0_f64;
    for c in factor.iter().rev() {
        let c = c.to_f64().unwrap_or(f64::NAN);
        value = value * x + c;
        scale = scale * x.abs() + c.abs();
    }
    value.is_finite() && value.abs() <= 1e-7 * scale.max(f64::MIN_POSITIVE)
}

/// Candidates for the real roots, accepted when they are roots and cover
/// every real root of `objects` exactly once.
fn select(
    graph: &Graph,
    factor: &[BigInt],
    candidates: Vec<NodeId>,
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    let mut chosen: Vec<Option<(f64, NodeId)>> = vec![None; objects.len()];
    for node in candidates {
        let Some(value) = value_of(graph, node) else { continue };
        if !is_root(factor, value) {
            continue;
        }
        let slot = objects.iter().position(|(v, _)| (v - value).abs() <= 1e-6 * (1.0 + v.abs()))?;
        if let Some(entry) = chosen.get_mut(slot)
            && entry.is_none()
        {
            *entry = Some((value, node));
        }
    }
    chosen.into_iter().collect()
}

fn rational_node(
    graph: &mut Graph,
    q: &BigRational,
) -> NodeId {
    graph.num(Number::rat(q.clone()))
}

fn mul(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    graph.node(core::MUL, factors)
}

/// `shift + z / scale`.
fn original_variable(
    graph: &mut Graph,
    analysis: &Analysis,
    z: NodeId,
) -> NodeId {
    let inverse = rational_node(graph, &BigRational::new(BigInt::from(1), analysis.scale.clone()));
    let scaled = mul(graph, &[inverse, z]);
    if analysis.shift.is_zero() {
        scaled
    } else {
        let shift = rational_node(graph, &analysis.shift);
        graph.node(core::ADD, &[shift, scaled])
    }
}

fn gauss_period_roots(
    graph: &mut Graph,
    factor: &[BigInt],
    analysis: &Analysis,
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    let cyclo = cyclotomic(analysis)?;
    let (cos, pi) = (graph.ops().lookup("cos")?, graph.ops().lookup("pi")?);
    let pi_node = graph.node(pi, &[]);
    let p = i64::try_from(cyclo.conductor).ok()?;
    let mut candidates = Vec::new();
    for (_, terms) in &cyclo.roots {
        let mut sum_terms = Vec::new();
        for (j, q) in terms {
            let angle_factor = graph.num(Number::fraction(2 * i64::try_from(*j).ok()?, p)?);
            let angle = mul(graph, &[angle_factor, pi_node]);
            let wave = graph.node(cos, &[angle]);
            let coefficient = rational_node(graph, &(q * BigRational::from_integer(BigInt::from(2))));
            sum_terms.push(mul(graph, &[coefficient, wave]));
        }
        let z = graph.node(core::ADD, &sum_terms);
        candidates.push(original_variable(graph, analysis, z));
    }
    select(graph, factor, candidates, objects)
}

fn tower_roots(
    graph: &mut Graph,
    factor: &[BigInt],
    analysis: &Analysis,
    tower: &Tower,
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    let theta: Vec<BigRational> = tower.theta.iter().cloned().map(BigRational::from_integer).collect();
    let thetas = super::rational_polynomial_roots(graph, &theta)?;
    let x_symbol = graph.interner_mut().fresh_symbol("x");
    let x_node = graph.symbol_node(x_symbol);
    let cap = 20_000;
    let mut candidates = Vec::new();
    for u in thetas {
        let mut gens = Gens::default();
        let gt = gens.index(graph, u);
        let gx = gens.index(graph, x_node);
        let x = Poly::generator(gx);
        let t = Poly::generator(gt);
        let mut h = x.pow(u32::try_from(tower.block_size).ok()?, cap)?;
        for (m, coefficients) in tower.block.iter().enumerate() {
            let mut c = Poly::zero();
            for (power, q) in coefficients.iter().enumerate() {
                if q.is_zero() {
                    continue;
                }
                let term = t.pow(u32::try_from(power).ok()?, cap)?.scale(&Number::rat(q.clone()));
                c = c.add(&term);
            }
            h = h.add(&c.mul(&x.pow(u32::try_from(m).ok()?, cap)?, cap)?);
        }
        for z in super::symbolic::roots(graph, &gens, &h, gx, 0)? {
            candidates.push(original_variable(graph, analysis, z));
        }
    }
    select(graph, factor, candidates, objects)
}

/// Nodes for the roots (as `(value, node)`) of an integer polynomial that
/// match the given numeric values one to one.
fn contains_op(
    graph: &Graph,
    node: NodeId,
    op: crate::graph::OpId,
) -> bool {
    graph.op(node) == op || graph.children(node).iter().any(|&c| contains_op(graph, c, op))
}

/// `R_1 … R_6` (in this order) as closed forms for the `Paired` Lagrange data:
/// `u_j` from the cubic, `R_j, R_{-j}` the roots of `w² - s(u_j) w + u_j^7`,
/// the larger one by the formula without cancellation and the smaller one as
/// `u_j^7` divided by it.
fn paired_resolvents(
    graph: &mut Graph,
    u_poly: &[BigInt],
    s_poly: &[BigRational],
    terms: &[(f64, f64, f64)],
) -> Option<Vec<NodeId>> {
    let targets: Vec<f64> = terms.iter().map(|t| t.0).collect();
    let u_nodes = matched_roots(graph, u_poly, &targets)?;
    let half = graph.num(Number::fraction(1, 2)?);
    let minus_one = graph.int(-1);
    let mut first: Vec<NodeId> = Vec::new();
    let mut second: Vec<NodeId> = Vec::new();
    for (u, (_, rj, rmj)) in u_nodes.iter().zip(terms) {
        let mut s_terms = Vec::new();
        for (power, c) in s_poly.iter().enumerate() {
            let coefficient = rational_node(graph, c);
            let exponent = graph.int(i64::try_from(power).ok()?);
            let p = if power == 0 { graph.int(1) } else { graph.node(core::POW, &[*u, exponent]) };
            s_terms.push(mul(graph, &[coefficient, p]));
        }
        let s_node = graph.node(core::ADD, &s_terms);
        let seven = graph.int(7);
        let u7 = graph.node(core::POW, &[*u, seven]);
        let two = graph.int(2);
        let s_squared = graph.node(core::POW, &[s_node, two]);
        let four = graph.int(4);
        let four_u7 = mul(graph, &[four, u7]);
        let neg = mul(graph, &[minus_one, four_u7]);
        let discriminant = graph.node(core::ADD, &[s_squared, neg]);
        let root = graph.node(core::POW, &[discriminant, half]);
        let s_value = rj + rmj;
        let signed_root = if s_value >= 0.0 { root } else { mul(graph, &[minus_one, root]) };
        let big_sum = graph.node(core::ADD, &[s_node, signed_root]);
        let big = mul(graph, &[half, big_sum]);
        let inverse = graph.node(core::POW, &[big, minus_one]);
        let small = mul(graph, &[u7, inverse]);
        // `big` is the one of larger magnitude.
        if rj.abs() >= rmj.abs() {
            first.push(big);
            second.push(small);
        } else {
            first.push(small);
            second.push(big);
        }
    }
    // R_1, R_2, R_3, then R_6, R_5, R_4.
    second.reverse();
    first.extend(second);
    Some(first)
}

fn matched_roots(
    graph: &mut Graph,
    polynomial: &[BigInt],
    targets: &[f64],
) -> Option<Vec<NodeId>> {
    let rational: Vec<BigRational> = polynomial.iter().cloned().map(BigRational::from_integer).collect();
    let nodes = super::rational_polynomial_roots(graph, &rational)?;
    // Unsolved factors stay `root_of` objects, which are of no use inside a larger
    // closed form (and expensive to evaluate for large coefficients).
    if let Some(root_of) = graph.ops().lookup("root_of")
        && nodes.iter().any(|&n| contains_op(graph, n, root_of))
    {
        return None;
    }
    let mut valued: Vec<(f64, NodeId)> = Vec::new();
    for node in nodes {
        valued.push((value_of(graph, node)?, node));
    }
    // A repeated root of the polynomial (zero, in the Dickson case) is one
    // node matched by several targets.
    let mut out = Vec::new();
    for &target in targets {
        let slot = valued.iter().position(|(v, _)| (v - target).abs() <= 1e-7 * (1.0 + target.abs()))?;
        out.push(valued.get(slot)?.1);
    }
    Some(out)
}

/// The real roots of a septic with group `D7`, `F21` or `F42` through the
/// Lagrange resolvents.
fn lagrange_roots(
    graph: &mut Graph,
    factor: &[BigInt],
    analysis: &Analysis,
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    let data = lagrange(analysis)?;
    let seventh = Number::fraction(1, 7)?;
    let mut candidates = Vec::new();
    match data {
        | Lagrange::OneReal { phi, resolvents, trace } => {
            let nodes = matched_roots(graph, &phi, &resolvents)?;
            let mut terms = vec![graph.num(Number::Int(trace))];
            for node in nodes {
                terms.push(super::symbolic::real_root(graph, node, 7)?);
            }
            let sum = graph.node(core::ADD, &terms);
            let factor_node = graph.num(seventh);
            let z = mul(graph, &[factor_node, sum]);
            candidates.push(original_variable(graph, analysis, z));
        },
        | Lagrange::Paired { u_poly, s_poly, terms, trace } => {
            let resolvents = paired_resolvents(graph, &u_poly, &s_poly, &terms)?;
            let mut sum_terms = vec![graph.num(Number::Int(trace))];
            for node in resolvents {
                sum_terms.push(super::symbolic::real_root(graph, node, 7)?);
            }
            let sum = graph.node(core::ADD, &sum_terms);
            let factor_node = graph.num(seventh);
            let z = mul(graph, &[factor_node, sum]);
            candidates.push(original_variable(graph, analysis, z));
        },
        | Lagrange::AllReal { s_poly, n_poly, terms, trace } => {
            let targets: Vec<f64> = terms.iter().map(|t| t.1).collect();
            let s_nodes = matched_roots(graph, &s_poly, &targets)?;
            let (acos, cos, pi) = (graph.ops().lookup("acos")?, graph.ops().lookup("cos")?, graph.ops().lookup("pi")?);
            let pi_node = graph.node(pi, &[]);
            let half = graph.num(Number::fraction(1, 2)?);
            let minus_half = graph.num(Number::fraction(-1, 2)?);
            let fourteenth = graph.num(Number::fraction(1, 14)?);
            // n_j and the angles.
            let mut pieces: Vec<(NodeId, NodeId, usize, i64)> = Vec::new();
            for ((s_node, (j, _, m)), _) in s_nodes.iter().zip(&terms).zip(0..) {
                let mut n_terms = Vec::new();
                for (power, c) in n_poly.iter().enumerate() {
                    let coefficient = rational_node(graph, c);
                    let exponent = graph.int(i64::try_from(power).ok()?);
                    let p = if power == 0 { graph.int(1) } else { graph.node(core::POW, &[*s_node, exponent]) };
                    n_terms.push(mul(graph, &[coefficient, p]));
                }
                let n_node = graph.node(core::ADD, &n_terms);
                let amplitude = graph.node(core::POW, &[n_node, fourteenth]);
                let inverse_root = graph.node(core::POW, &[n_node, minus_half]);
                let cosine = mul(graph, &[half, *s_node, inverse_root]);
                let phi = graph.node(acos, &[cosine]);
                pieces.push((amplitude, phi, *j, *m));
            }
            for k in 0..7_i64 {
                let mut sum_terms = vec![graph.num(Number::Int(trace.clone()))];
                for (amplitude, phi, j, m) in &pieces {
                    let two_m = graph.int(2 * m);
                    let turn = mul(graph, &[two_m, pi_node]);
                    let numerator = graph.node(core::ADD, &[*phi, turn]);
                    let seventh_node = graph.num(Number::fraction(1, 7)?);
                    let first = mul(graph, &[seventh_node, numerator]);
                    let shift_coefficient = graph.num(Number::fraction(-2 * i64::try_from(*j).ok()? * k, 7)?);
                    let second = mul(graph, &[shift_coefficient, pi_node]);
                    let angle = graph.node(core::ADD, &[first, second]);
                    let wave = graph.node(cos, &[angle]);
                    let two = graph.int(2);
                    sum_terms.push(mul(graph, &[two, *amplitude, wave]));
                }
                let sum = graph.node(core::ADD, &sum_terms);
                let factor_node = graph.num(seventh.clone());
                let z = mul(graph, &[factor_node, sum]);
                candidates.push(original_variable(graph, analysis, z));
            }
        },
    }
    select(graph, factor, candidates, objects)
}

/// `exp(2 π i u / m)` as a term.
fn unit_root(
    graph: &mut Graph,
    u: i64,
    m: i64,
) -> Option<NodeId> {
    let (cos, sin, pi, i_unit) =
        (graph.ops().lookup("cos")?, graph.ops().lookup("sin")?, graph.ops().lookup("pi")?, graph.ops().lookup("I")?);
    let pi_node = graph.node(pi, &[]);
    let angle_factor = graph.num(Number::fraction(2 * u, m)?);
    let angle = mul(graph, &[angle_factor, pi_node]);
    let c = graph.node(cos, &[angle]);
    let s = graph.node(sin, &[angle]);
    let i_node = graph.node(i_unit, &[]);
    let imaginary = mul(graph, &[i_node, s]);
    Some(graph.node(core::ADD, &[c, imaginary]))
}

/// All complex roots of an irreducible integer polynomial of degree 6 or 7
/// with a solvable Galois group, as `((re, im), term)`, when they have an
/// exact form here: abelian groups through `exp(2πi j/m)` combinations
/// (Gauss periods, also for fields that are not real), and the septic groups
/// `D7`, `F21`, `F42` through the Lagrange resolvents (`ζ_7` powers times real
/// seventh roots).
pub(super) fn complex_roots(
    graph: &mut Graph,
    factor: &[BigInt],
) -> Option<Vec<((f64, f64), NodeId)>> {
    if !(7..=8).contains(&factor.len()) {
        return None;
    }
    let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
    let analysis = analyse(&rational)?;
    if !analysis.group.solvable {
        return None;
    }
    let degree = factor.len() - 1;
    let mut candidates: Vec<NodeId> = Vec::new();
    if let Some(cyclo) = cyclotomic_complex(&analysis) {
        let m = i64::try_from(cyclo.conductor).ok()?;
        for (_, terms) in &cyclo.roots {
            let mut sum_terms = Vec::new();
            for (u, w) in terms {
                let wave = unit_root(graph, i64::try_from(*u).ok()?, m)?;
                let coefficient = rational_node(graph, w);
                sum_terms.push(mul(graph, &[coefficient, wave]));
            }
            let z = graph.node(core::ADD, &sum_terms);
            candidates.push(original_variable(graph, &analysis, z));
        }
    } else if degree == 7 {
        let (nodes, trace) = match lagrange(&analysis)? {
            | Lagrange::OneReal { phi, resolvents, trace } => (matched_roots(graph, &phi, &resolvents)?, trace),
            | Lagrange::Paired { u_poly, s_poly, terms, trace } => (paired_resolvents(graph, &u_poly, &s_poly, &terms)?, trace),
            | Lagrange::AllReal { .. } => return None,
        };
        let mut roots7 = Vec::new();
        for node in nodes {
            roots7.push(super::symbolic::real_root(graph, node, 7)?);
        }
        let seventh = graph.num(Number::fraction(1, 7)?);
        for k in 0..7_i64 {
            let mut terms = vec![graph.num(Number::Int(trace.clone()))];
            for (j0, t) in roots7.iter().enumerate() {
                let j = i64::try_from(j0 + 1).ok()?;
                // ζ^{-jk}
                let phase = unit_root(graph, (-j * k).rem_euclid(7), 7)?;
                terms.push(mul(graph, &[phase, *t]));
            }
            let sum = graph.node(core::ADD, &terms);
            let z = mul(graph, &[seventh, sum]);
            candidates.push(original_variable(graph, &analysis, z));
        }
    } else {
        return None;
    }
    // Keep the candidates that are roots and distinct.
    let mut out: Vec<((f64, f64), NodeId)> = Vec::new();
    for node in candidates {
        let Some(z) = graph.eval_complex(node, &std::collections::HashMap::new()) else { continue };
        let mut value = num_complex::Complex64::new(0.0, 0.0);
        let mut scale = 0.0_f64;
        for c in factor.iter().rev() {
            let c = c.to_f64()?;
            value = value * z + c;
            scale = scale * z.norm() + c.abs();
        }
        if !value.re.is_finite() || value.norm() > 1e-7 * scale.max(f64::MIN_POSITIVE) {
            return None;
        }
        if !out.iter().any(|((re, im), _)| (re - z.re).hypot(im - z.im) < 1e-9 * (1.0 + z.norm())) {
            out.push(((z.re, z.im), node));
        }
    }
    (out.len() == degree).then_some(out)
}

/// Radical or trigonometric forms of the real roots of the irreducible
/// integer polynomial `factor` (degree 6 or 7) when its Galois group
/// allows them; `objects` are its real roots as `root_of` objects.
pub(super) fn roots(
    graph: &mut Graph,
    factor: &[BigInt],
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    if objects.is_empty() || !(7..=8).contains(&factor.len()) {
        return None;
    }
    let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
    let analysis = analyse(&rational)?;
    if !analysis.group.solvable {
        return None;
    }
    if let Some(found) = gauss_period_roots(graph, factor, &analysis, objects) {
        return Some(found);
    }
    if factor.len() == 8
        && let Some(found) = lagrange_roots(graph, factor, &analysis, objects)
    {
        return Some(found);
    }
    if factor.len() == 7 {
        for tower in towers(&analysis) {
            if let Some(found) = tower_roots(graph, factor, &analysis, &tower, objects) {
                return Some(found);
            }
        }
    }
    None
}
