//! Irreducible factors of degree five and higher of rational polynomials.
//!
//! A quintic is analysed with [`analyse_quintic`] (Dummit's sextic
//! resolvent). When it is solvable and has one real root, that root is
//! written in radicals: `r = shift + (1/(5 λ)) Σ R^{1/5}` over the four
//! roots `R` of the quartic `Φ` with coefficients in `Q(√5)` (see the
//! `galois` module), `Φ` being solved by Ferrari's formula; the result is
//! checked numerically against the real root. Every other real root of an
//! irreducible factor of degree five or more — unsolvable quintics, solvable
//! ones with five real roots, higher degrees — is returned as an exact
//! `root_of(f, x, k)` object (see the `rootof` module), so the list of
//! solutions is always complete.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::symbolic;
use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::galois::analyse_quintic;
use crate::rules::poly::galois::Radical;
use crate::rules::poly::galois::Trigonometric;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::rootof::real_roots_as_objects;

/// The real roots of the irreducible integer polynomial `factor`
/// (ascending coefficients, degree at least five) with their values.
pub(super) fn roots(
    graph: &mut Graph,
    factor: &[BigInt],
) -> Option<Vec<(f64, NodeId)>> {
    let objects = real_roots_as_objects(graph, factor)?;
    if super::numeric_run() {
        return Some(objects.into_iter().map(|(value, _)| (value, graph.float(value))).collect());
    }
    if factor.len() == 6 {
        let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
        if let Some(analysis) = analyse_quintic(&rational)
            && let Some(radical) = &analysis.radical
        {
            if let Some(trig) = &radical.trigonometric {
                if objects.len() == 5
                    && let Some(found) = trigonometric_roots(graph, radical, trig, &objects)
                {
                    return Some(found);
                }
            } else if let [(value, _)] = objects.as_slice()
                && let Some(node) = radical_root(graph, radical, *value)
            {
                return Some(vec![(*value, node)]);
            }
        }
    }
    Some(objects)
}

fn pair_poly(
    index: u32,
    (u, v): &(BigRational, BigRational),
) -> Poly {
    Poly::constant(Number::rat(u.clone())).add(&Poly::generator(index).scale(&Number::rat(v.clone())))
}

fn radical_root(
    graph: &mut Graph,
    radical: &Radical,
    expected: f64,
) -> Option<NodeId> {
    let (five, half) = (graph.int(5), graph.num(Number::fraction(1, 2)?));
    let sqrt5 = graph.node(core::POW, &[five, half]);
    let mut gens = Gens::default();
    let g5 = gens.index(graph, sqrt5);
    let mut coefficients: Vec<Poly> = radical.phi.iter().map(|p| pair_poly(g5, p)).collect();
    coefficients.push(Poly::constant(Number::from(1)));
    let candidates = symbolic::ferrari(graph, &gens, &coefficients)?;
    let env = Env::numeric(0.0);
    let numeric: Vec<f64> = radical
        .phi
        .iter()
        .map(|(u, v)| u.to_f64().unwrap_or(f64::NAN) + v.to_f64().unwrap_or(f64::NAN) * 5.0_f64.sqrt())
        .collect();
    let phi = |w: f64| -> (f64, f64) {
        let mut value = w.powi(4);
        let mut scale = w.powi(4).abs();
        for (k, c) in numeric.iter().enumerate() {
            let term = c * w.powi(i32::try_from(k).unwrap_or(0));
            value += term;
            scale += term.abs();
        }
        (value, scale)
    };
    let mut distinct: Vec<(f64, NodeId)> = Vec::new();
    for candidate in candidates {
        let Some(w) = graph.eval(candidate, &env).filter(|w| w.is_finite()) else { continue };
        let (value, scale) = phi(w);
        if value.abs() > 1e-7 * scale {
            continue;
        }
        if !distinct.iter().any(|(v, _)| (v - w).abs() <= 1e-9 * (1.0 + w.abs())) {
            distinct.push((w, candidate));
        }
    }
    // A zero root of Φ (exactly representable) can fall in a degenerate
    // case of the quartic formula.
    if radical.phi.first().is_some_and(|(u, v)| u.is_zero() && v.is_zero())
        && !distinct.iter().any(|(w, _)| w.abs() < 1e-9)
    {
        distinct.push((0.0, graph.int(0)));
    }
    // Multiplicities from the derivatives of Φ (R_j can coincide).
    let mut full: Vec<f64> = numeric.clone();
    full.push(1.0);
    let derivative_at = |w: f64, order: usize| -> (f64, f64) {
        let mut coefficients = full.clone();
        for _ in 0..order {
            coefficients = coefficients.iter().enumerate().skip(1).map(|(i, c)| c * i as f64).collect();
        }
        let mut value = 0.0;
        let mut scale = 0.0;
        for (k, c) in coefficients.iter().enumerate() {
            let term = c * w.powi(i32::try_from(k).unwrap_or(0));
            value += term;
            scale += term.abs();
        }
        (value, scale)
    };
    let mut terms = Vec::new();
    let mut total = 0_usize;
    for (w, node) in distinct {
        let mut multiplicity = 1_usize;
        while multiplicity < 4 {
            let (value, scale) = derivative_at(w, multiplicity);
            if value.abs() <= 1e-7 * scale.max(1e-300) {
                multiplicity += 1;
            } else {
                break;
            }
        }
        total += multiplicity;
        let root = symbolic::real_root(graph, node, 5)?;
        if multiplicity == 1 {
            terms.push(root);
        } else {
            let m = graph.int(i64::try_from(multiplicity).ok()?);
            terms.push(graph.node(core::MUL, &[m, root]));
        }
    }
    if total != 4 {
        return None;
    }
    let sum = graph.node(core::ADD, &terms);
    let factor = graph.num(Number::rat(BigRational::new(BigInt::one(), BigInt::from(5) * &radical.scale)));
    let scaled = graph.node(core::MUL, &[factor, sum]);
    let total = if radical.shift.is_zero() {
        scaled
    } else {
        let shift = graph.num(Number::rat(radical.shift.clone()));
        graph.node(core::ADD, &[shift, scaled])
    };
    let value = graph.eval(total, &env)?;
    ((value - expected).abs() <= 1e-6 * (1.0 + expected.abs())).then_some(total)
}

/// The five real roots of a solvable quintic as cosines.
///
/// With `R_1 = a_1 + i b_1` and `R_2 = a_2 + i b_2` the complex roots of `Φ`
/// (found as the roots of the quadratic factors in Ferrari's method) and
/// `t_j` fifth roots of `R_j`, the roots are
/// `z_k = (2/5) Σ_j |R_j|^(1/5) cos((φ_j + 2π (m_j - j k))/5)` with
/// `φ_j = arg R_j` and `m_j` fixed by the numerical values of `t_j`.
fn trigonometric_roots(
    graph: &mut Graph,
    radical: &Radical,
    trig: &Trigonometric,
    objects: &[(f64, NodeId)],
) -> Option<Vec<(f64, NodeId)>> {
    let (five, half) = (graph.int(5), graph.num(Number::fraction(1, 2)?));
    let sqrt5 = graph.node(core::POW, &[five, half]);
    let mut gens = Gens::default();
    let g5 = gens.index(graph, sqrt5);
    let mut coefficients: Vec<Poly> = radical.phi.iter().map(|p| pair_poly(g5, p)).collect();
    coefficients.push(Poly::constant(Number::from(1)));
    let mut pairs = Vec::new();
    symbolic::ferrari_with_pairs(graph, &gens, &coefficients, &mut pairs)?;
    let env = Env::numeric(0.0);
    let (acos, cos, pi) = (graph.ops().lookup("acos")?, graph.ops().lookup("cos")?, graph.ops().lookup("pi")?);
    let pi_node = graph.node(pi, &[]);
    let tau = std::f64::consts::TAU;
    // For each of R_1, R_2: the pair of nodes whose values match.
    let mut parts: Vec<(NodeId, NodeId, f64, i64, i64)> = Vec::new();
    for (j, (r, t)) in trig.r.iter().zip(&trig.t).enumerate() {
        let target = (r.0, r.1.abs());
        let mut found = None;
        for &(a_node, b_node) in &pairs {
            let (Some(a), Some(b)) = (graph.eval(a_node, &env), graph.eval(b_node, &env)) else { continue };
            if (a - target.0).abs() <= 1e-7 * (1.0 + target.0.abs()) && (b - target.1).abs() <= 1e-7 * (1.0 + target.1.abs()) {
                found = Some((a_node, b_node));
                break;
            }
        }
        let (a_node, b_node) = found?;
        // u = t if Im R > 0 else conj t: a fifth root of a + i|b|.
        let u = (t.0, if r.1 > 0.0 { t.1 } else { -t.1 });
        let phi = target.1.atan2(target.0);
        let m = ((5.0 * u.1.atan2(u.0) - phi) / tau).round();
        // Conjugating t_j turns ζ^{-jk} into ζ^{jk}.
        let weight = if r.1 > 0.0 { i64::try_from(j + 1).ok()? } else { -i64::try_from(j + 1).ok()? };
        parts.push((a_node, b_node, phi, m as i64, weight));
    }
    let mut out = Vec::new();
    for k in 0..5_i64 {
        let mut terms = Vec::new();
        for &(a_node, b_node, _phi, m, weight) in &parts {
            let two = graph.int(2);
            let a2 = graph.node(core::POW, &[a_node, two]);
            let b2 = graph.node(core::POW, &[b_node, two]);
            let modulus_sq = graph.node(core::ADD, &[a2, b2]);
            let tenth = graph.num(Number::fraction(1, 10)?);
            let amplitude = graph.node(core::POW, &[modulus_sq, tenth]);
            let minus_half = graph.num(Number::fraction(-1, 2)?);
            let inverse_modulus = graph.node(core::POW, &[modulus_sq, minus_half]);
            let ratio = graph.node(core::MUL, &[a_node, inverse_modulus]);
            let phi_node = graph.node(acos, &[ratio]);
            let turns = graph.int(2 * (m - weight * k));
            let turn_pi = graph.node(core::MUL, &[turns, pi_node]);
            let angle = graph.node(core::ADD, &[phi_node, turn_pi]);
            let fifth = graph.num(Number::fraction(1, 5)?);
            let divided = graph.node(core::MUL, &[fifth, angle]);
            let wave = graph.node(cos, &[divided]);
            terms.push(graph.node(core::MUL, &[amplitude, wave]));
        }
        let sum = graph.node(core::ADD, &terms);
        let factor = graph.num(Number::rat(BigRational::new(BigInt::from(2), BigInt::from(5) * &radical.scale)));
        let scaled = graph.node(core::MUL, &[factor, sum]);
        let total = if radical.shift.is_zero() {
            scaled
        } else {
            let shift = graph.num(Number::rat(radical.shift.clone()));
            graph.node(core::ADD, &[shift, scaled])
        };
        out.push((graph.eval(total, &env)?, total));
    }
    // Every computed root must match one of the real roots.
    let mut used = vec![false; objects.len()];
    for (value, _) in &out {
        let slot = objects.iter().position(|(v, _)| (v - value).abs() <= 1e-6 * (1.0 + v.abs()))?;
        if std::mem::replace(used.get_mut(slot)?, true) {
            return None;
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(out)
}
