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
    if factor.len() == 6 {
        let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
        if let Some(analysis) = analyse_quintic(&rational)
            && let Some(radical) = &analysis.radical
            && let [(value, _)] = objects.as_slice()
            && let Some(node) = radical_root(graph, radical, *value)
        {
            return Some(vec![(*value, node)]);
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
    if distinct.len() != 4 {
        return None;
    }
    let mut terms = Vec::new();
    for (_, node) in distinct {
        terms.push(symbolic::real_root(graph, node, 5)?);
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
