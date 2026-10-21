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
    if factor.len() == 7 {
        for tower in towers(&analysis) {
            if let Some(found) = tower_roots(graph, factor, &analysis, &tower, objects) {
                return Some(found);
            }
        }
    }
    None
}
