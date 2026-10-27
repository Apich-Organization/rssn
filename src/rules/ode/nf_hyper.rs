//! Hyperexponential solutions over a number field: Kovacic's case 1 with
//! irrational local data.
//!
//! A solution `y = P(x) Π (x - c)^ρ exp(E(x))` of a linear equation with
//! polynomial coefficients is sought with `P`, `ρ` and the coefficients of `E`
//! in one number field `L = Q(β)` (typically quadratic: exponents such as
//! `(1 ± √5)/2` at a regular singular point). The local data at every rational
//! singular point and at infinity come from the Newton-polygon expansion
//! ([`super::formal`]); candidates over `Q` and over a field isomorphic to
//! `L` (identified through a root of its modulus in `L`) are combined, the
//! behaviour at infinity fixes the degree of `P`, and `P` solves a linear
//! system over `L` — a `Q`-linear system in the coordinates of its
//! coefficients. For a quadratic `L` the conjugate candidates are tried as
//! well, so the two conjugate solutions come out together.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::formal::formal_solutions;
use super::formal::roots_in_field;
use super::formal::Elt;
use super::formal::Field;
use super::formal::FormalSolution;
use super::formal::Point as FormalPoint;
use super::add;
use super::exp_of;
use super::mul;
use crate::graph::op::core;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::rational_function_in;
use crate::rules::poly::univariate;

type Q = BigRational;
type Poly = Vec<Elt>;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn trim_poly_elt(e: &[Q]) -> Elt {
    let mut v = e.to_vec();
    while v.last().is_some_and(Zero::is_zero) {
        v.pop();
    }
    v
}

fn trim_poly(mut p: Poly) -> Poly {
    while p.last().is_some_and(Vec::is_empty) {
        p.pop();
    }
    p
}

fn padd(
    f: &Field,
    a: &[Elt],
    b: &[Elt],
) -> Poly {
    let n = a.len().max(b.len());
    trim_poly((0..n).map(|i| f.add(a.get(i).map_or(&[][..], Vec::as_slice), b.get(i).map_or(&[][..], Vec::as_slice))).collect())
}

fn pscale(
    f: &Field,
    a: &[Elt],
    c: &[Q],
) -> Poly {
    trim_poly(a.iter().map(|x| f.mul(x, c)).collect())
}

fn pmul(
    f: &Field,
    a: &[Elt],
    b: &[Elt],
) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out: Poly = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i + j) {
                *slot = f.add(slot, &f.mul(x, y));
            }
        }
    }
    trim_poly(out)
}

fn pderiv(
    f: &Field,
    a: &[Elt],
) -> Poly {
    trim_poly(a.iter().enumerate().skip(1).map(|(i, c)| f.mul(c, &[q(i64::try_from(i).unwrap_or(0))])).collect())
}

/// A polynomial with rational coefficients as a polynomial over the field.
fn embed(
    f: &Field,
    p: &[Q],
) -> Poly {
    trim_poly(p.iter().map(|c| f.from_q(c)).collect())
}

/// A candidate at one point: the exponent part `E` as `numerator /
/// (x - c)^power` (a sum of `κ (x - c)^{-s}`), the exponent `ρ`.
#[derive(Clone)]
struct Local {
    point: Q,
    /// `E = Σ_s e[s] (x - c)^{-s}` with `s = 1…`.
    e: Vec<Elt>,
    rho: Elt,
}

/// At infinity: `E = Σ_s e[s] x^s` and `ρ`.
#[derive(Clone)]
struct AtInfinity {
    e: Vec<Elt>,
    rho: Elt,
}

/// A hyperexponential solution over `L`.
#[derive(Clone)]
pub(super) struct HyperK {
    pub(super) field: Field,
    pub(super) poly: Poly,
    /// `(c, ρ)`.
    pub(super) factors: Vec<(Q, Elt)>,
    /// `E = numerator / denominator`.
    pub(super) exponent: (Poly, Poly),
}

/// The image of an element of `source` under the isomorphism onto `target`
/// that sends the generator of `source` to `image`.
fn map_element(
    target: &Field,
    source: &Field,
    image: &[Q],
    e: &[Q],
) -> Elt {
    let _ = target;
    if source.is_rational() {
        return trim_poly_elt(e);
    }
    let mut value: Elt = Vec::new();
    for c in e.iter().rev() {
        value = target.add(&target.mul(&value, image), &target.from_q(c));
    }
    value
}

fn field_matches(
    target: &Field,
    source: &Field,
) -> Option<Elt> {
    if source.is_rational() {
        return Some(Vec::new());
    }
    if *target == *source {
        return Some(target.generator());
    }
    if target.degree() != source.degree() {
        return None;
    }
    let polynomial: Vec<Elt> = source.modulus.iter().map(|c| target.from_q(c)).collect();
    roots_in_field(target, &polynomial).into_iter().next()
}

/// The conjugate of an element of a quadratic field.
fn conjugate(
    field: &Field,
    e: &[Q],
) -> Elt {
    if field.degree() != 2 {
        return e.to_vec();
    }
    // The other root of x^2 + b x + c is -b - α.
    let other = field.sub(&field.from_q(&-field.modulus.get(1).cloned().unwrap_or_else(Q::zero)), &field.generator());
    map_element(field, field, &other, e)
}

fn local_candidates(
    op: &[Vec<Q>],
    point: &FormalPoint,
    target: &mut Option<Field>,
    infinity: bool,
) -> (Vec<Local>, Vec<AtInfinity>) {
    let mut near = Vec::new();
    let mut far = Vec::new();
    let c = match point {
        | FormalPoint::Finite(f) if f.len() == 2 => -f.first().cloned().unwrap_or_else(Q::zero),
        | _ => Q::zero(),
    };
    for solution in formal_solutions(op, point, 2) {
        if solution.ramification != 1 || solution.has_logarithms() {
            continue;
        }
        let FormalSolution { field, exponent, rho, .. } = &solution;
        // Fix the global field at the first extension met.
        if !field.is_rational() && target.is_none() {
            *target = Some(field.clone());
        }
        let l = target.clone().unwrap_or_else(Field::rationals);
        let Some(image) = field_matches(&l, field) else { continue };
        let map = |e: &[Q]| map_element(&l, field, &image, e);
        let mut series: Vec<Elt> = Vec::new();
        for (s, kappa) in exponent {
            let Ok(s) = usize::try_from(*s) else { continue };
            while series.len() <= s {
                series.push(Vec::new());
            }
            if let Some(slot) = series.get_mut(s) {
                *slot = map(kappa);
            }
        }
        let rho_l = map(rho);
        // Every conjugate under the quadratic field's automorphism.
        let variants: Vec<(Vec<Elt>, Elt)> = if l.degree() == 2 && !l.is_rational() {
            vec![(series.clone(), rho_l.clone()), (series.iter().map(|e| conjugate(&l, e)).collect(), conjugate(&l, &rho_l))]
        } else {
            vec![(series, rho_l)]
        };
        for (e, rho) in variants {
            if infinity {
                far.push(AtInfinity { e, rho });
            } else {
                near.push(Local { point: c.clone(), e, rho });
            }
        }
    }
    (near, far)
}

/// Linear algebra over `L`: the polynomials `P` of degree at most `d` over `L`
/// with `Σ_i a_i (D + ω)^i P = 0`, `ω = ω_n / D_ω`, as coordinates over `Q`.
fn solve_polynomial(
    field: &Field,
    op: &[Vec<Q>],
    omega_n: &[Elt],
    omega_d: &[Elt],
    d: usize,
) -> Vec<Poly> {
    let n = op.len().saturating_sub(1);
    let degree = field.degree();
    // For each monomial x^j: Σ_i a_i N_i D^{n-i}, where Q_{i+1} = (N_i' D - i N_i D' + ω_n N_i) / D^{i+1}.
    let mut images: Vec<Poly> = Vec::new();
    let d_powers: Vec<Poly> = {
        let mut v = vec![vec![field.from_q(&Q::one())]];
        for _ in 0..n {
            let next = pmul(field, v.last().map_or(&[][..], Vec::as_slice), omega_d);
            v.push(next);
        }
        v
    };
    for j in 0..=d {
        let mut monomial: Poly = vec![Vec::new(); j];
        monomial.push(field.from_q(&Q::one()));
        let mut current = monomial;
        let mut total: Poly = Vec::new();
        for (i, a) in op.iter().enumerate() {
            let term = pmul(field, &embed(field, a), &current);
            let scaled = pmul(field, &term, d_powers.get(n - i).map_or(&[][..], Vec::as_slice));
            total = padd(field, &total, &scaled);
            if i < n {
                let first = pmul(field, &pderiv(field, &current), omega_d);
                let second = pscale(field, &pmul(field, &current, &pderiv(field, omega_d)), &[q(-i64::try_from(i).unwrap_or(0))]);
                let third = pmul(field, omega_n, &current);
                current = padd(field, &padd(field, &first, &second), &third);
            }
        }
        images.push(total);
    }
    // Unknowns u_{j,t}: coefficient of α^t in the coefficient of x^j.
    let columns = (d + 1) * degree;
    let height = images.iter().map(Vec::len).max().unwrap_or(0);
    let generator_powers: Vec<Elt> = {
        let mut v = vec![field.from_q(&Q::one())];
        for _ in 1..degree {
            let next = field.mul(v.last().map_or(&[][..], Vec::as_slice), &field.generator());
            v.push(next);
        }
        v
    };
    let mut rows: Vec<Vec<Q>> = Vec::new();
    for m in 0..height {
        for coordinate in 0..degree {
            let mut row = Vec::with_capacity(columns);
            for image in &images {
                for power in &generator_powers {
                    let value = field.mul(power, image.get(m).map_or(&[][..], Vec::as_slice));
                    row.push(value.get(coordinate).cloned().unwrap_or_else(Q::zero));
                }
            }
            rows.push(row);
        }
    }
    super::linear::nullspace(rows, columns)
        .into_iter()
        .map(|u| {
            let mut poly: Poly = Vec::new();
            for j in 0..=d {
                let mut coefficient: Elt = Vec::new();
                for (t, power) in generator_powers.iter().enumerate() {
                    let weight = u.get(j * degree + t).cloned().unwrap_or_else(Q::zero);
                    coefficient = field.add(&coefficient, &field.mul(power, &[weight]));
                }
                poly.push(coefficient);
            }
            trim_poly(poly)
        })
        .filter(|p| !p.is_empty())
        .collect()
}

/// The hyperexponential solutions of the operator over a common number field.
pub(super) fn hyperexponential_over_field(op: &[Vec<Q>]) -> Vec<HyperK> {
    let mut out = Vec::new();
    let Some(lead) = op.last() else { return out };
    let (_, factors) = univariate::factor(lead);
    // Rational singular points only.
    let mut points: Vec<Q> = Vec::new();
    for (f, _) in &factors {
        if let [c0, c1] = f.as_slice() {
            points.push(Q::new(-c0.clone(), c1.clone()));
        } else if f.len() > 2 {
            return out;
        }
    }
    // First pass: the global field is the first extension met.
    let mut target: Option<Field> = None;
    for c in &points {
        let _ = local_candidates(op, &FormalPoint::Finite(vec![-c.clone(), Q::one()]), &mut target, false);
    }
    let _ = local_candidates(op, &FormalPoint::Infinity, &mut target, true);
    let Some(field) = target else { return out };
    // Second pass with the field fixed: candidates over Q and over fields isomorphic to it.
    let mut local2: Vec<Vec<Local>> = Vec::new();
    let mut fixed = Some(field.clone());
    for c in &points {
        let (near, _) = local_candidates(op, &FormalPoint::Finite(vec![-c.clone(), Q::one()]), &mut fixed, false);
        local2.push(near);
    }
    let (_, infinity) = local_candidates(op, &FormalPoint::Infinity, &mut fixed, true);
    let total: usize = local2.iter().map(Vec::len).product::<usize>().saturating_mul(infinity.len());
    if total == 0 || total > 4000 {
        return out;
    }
    let mut chosen = vec![0_usize; local2.len()];
    for far in &infinity {
        chosen.fill(0);
        loop {
            let selected: Vec<&Local> = local2.iter().zip(&chosen).filter_map(|(options, &i)| options.get(i)).collect();
            if selected.len() == local2.len()
                && let Some(found) = try_combination(&field, op, &selected, far)
            {
                out.extend(found);
            }
            let mut position = 0;
            loop {
                if position == local2.len() {
                    break;
                }
                let count = local2.get(position).map_or(0, Vec::len);
                let Some(slot) = chosen.get_mut(position) else { break };
                *slot += 1;
                if *slot < count {
                    break;
                }
                *slot = 0;
                position += 1;
            }
            if position == local2.len() {
                break;
            }
        }
    }
    out
}

fn try_combination(
    field: &Field,
    op: &[Vec<Q>],
    selected: &[&Local],
    far: &AtInfinity,
) -> Option<Vec<HyperK>> {
    // Degree of P: -ρ_∞ - Σ ρ_c.
    let mut sum = far.rho.clone();
    for l in selected {
        sum = field.add(&sum, &l.rho);
    }
    let sum = field.scale(&sum, &q(-1));
    let degree = field.as_rational(&sum)?;
    if !degree.is_integer() || degree < Q::zero() || degree.to_integer() > BigInt::from(30) {
        return None;
    }
    let d = usize::try_from(degree.to_integer().to_i64()?).ok()?;
    // E = Σ_c Σ_s κ (x - c)^{-s} + Σ_s κ_∞ x^s as En / Ed with Ed = Π (x - c)^{smax}.
    let one = field.from_q(&Q::one());
    let mut ed: Poly = vec![one.clone()];
    let mut en: Poly = Vec::new();
    for l in selected {
        let smax = l.e.len().saturating_sub(1);
        let factor: Poly = vec![field.from_q(&-l.point.clone()), one.clone()];
        let mut power: Poly = vec![one.clone()];
        for _ in 0..smax {
            power = pmul(field, &power, &factor);
        }
        // Existing numerator times the new denominator, plus the new part times the old denominator.
        let mut part: Poly = Vec::new();
        for (s, kappa) in l.e.iter().enumerate().skip(1) {
            if kappa.is_empty() {
                continue;
            }
            let mut cofactor: Poly = vec![one.clone()];
            for _ in 0..smax - s {
                cofactor = pmul(field, &cofactor, &factor);
            }
            part = padd(field, &part, &pscale(field, &cofactor, kappa));
        }
        en = padd(field, &pmul(field, &en, &power), &pmul(field, &part, &ed));
        ed = pmul(field, &ed, &power);
    }
    // Polynomial part at infinity.
    let mut poly_part: Poly = Vec::new();
    for (s, kappa) in far.e.iter().enumerate().skip(1) {
        while poly_part.len() <= s {
            poly_part.push(Vec::new());
        }
        if let Some(slot) = poly_part.get_mut(s) {
            slot.clone_from(kappa);
        }
    }
    en = padd(field, &en, &pmul(field, &poly_part, &ed));
    // ω = E' + Σ ρ_c / (x - c) over the denominator Ed^2 Π (x - c).
    let mut product_c: Poly = vec![one.clone()];
    for l in selected {
        product_c = pmul(field, &product_c, &[field.from_q(&-l.point.clone()), one.clone()]);
    }
    let ed2 = pmul(field, &ed, &ed);
    let omega_d = pmul(field, &ed2, &product_c);
    let e_prime_n = padd(
        field,
        &pmul(field, &pderiv(field, &en), &ed),
        &pscale(field, &pmul(field, &en, &pderiv(field, &ed)), &[q(-1)]),
    );
    let mut omega_n = pmul(field, &e_prime_n, &product_c);
    for (index, l) in selected.iter().enumerate() {
        let mut others: Poly = vec![one.clone()];
        for (j, m) in selected.iter().enumerate() {
            if j != index {
                others = pmul(field, &others, &[field.from_q(&-m.point.clone()), one.clone()]);
            }
        }
        omega_n = padd(field, &omega_n, &pscale(field, &pmul(field, &ed2, &others), &l.rho));
    }
    let solutions = solve_polynomial(field, op, &omega_n, &omega_d, d);
    if solutions.is_empty() {
        return None;
    }
    Some(
        solutions
            .into_iter()
            .map(|poly| HyperK {
                field: field.clone(),
                poly,
                factors: selected.iter().map(|l| (l.point.clone(), l.rho.clone())).collect(),
                exponent: (en.clone(), ed.clone()),
            })
            .collect(),
    )
}

// ---------------------------------------------------------------------
// Graph construction.

fn element_node(
    cx: &mut Cx<'_>,
    field: &Field,
    alpha: Option<NodeId>,
    e: &[Q],
) -> NodeId {
    let _ = field;
    let mut terms = Vec::new();
    for (i, c) in e.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let cn = cx.graph.num(Number::rat(c.clone()));
        match (i, alpha) {
            | (0, _) => terms.push(cn),
            | (_, Some(a)) => {
                let exponent = cx.graph.int(i64::try_from(i).unwrap_or(0));
                let power = cx.graph.node(core::POW, &[a, exponent]);
                terms.push(mul(cx.graph, &[cn, power]));
            },
            | _ => {},
        }
    }
    add(cx.graph, &terms)
}

fn poly_node(
    cx: &mut Cx<'_>,
    field: &Field,
    alpha: Option<NodeId>,
    x: NodeId,
    p: &[Elt],
) -> NodeId {
    let mut terms = Vec::new();
    for (j, c) in p.iter().enumerate() {
        if c.is_empty() {
            continue;
        }
        let cn = element_node(cx, field, alpha, c);
        let power = if j == 0 {
            cx.graph.int(1)
        } else {
            let e = cx.graph.int(i64::try_from(j).unwrap_or(0));
            cx.graph.node(core::POW, &[x, e])
        };
        terms.push(mul(cx.graph, &[cn, power]));
    }
    add(cx.graph, &terms)
}

/// The solutions as terms, one for each embedding of the field.
pub(super) fn solution_nodes(
    cx: &mut Cx<'_>,
    x: NodeId,
    hyper: &HyperK,
) -> Option<Vec<NodeId>> {
    let root_of = cx.graph.ops().lookup("root_of")?;
    let mut out = Vec::new();
    let embeddings = hyper.field.degree().max(1);
    for k in 0..embeddings {
        let alpha = if hyper.field.is_rational() {
            None
        } else {
            let mut arguments = vec![cx.graph.int(i64::try_from(k).ok()?)];
            for c in &hyper.field.modulus {
                arguments.push(cx.graph.num(Number::rat(c.clone())));
            }
            Some(cx.graph.node(root_of, &arguments))
        };
        let mut factors = Vec::new();
        for (c, rho) in &hyper.factors {
            if rho.is_empty() {
                continue;
            }
            let minus_c = cx.graph.num(Number::rat(-c.clone()));
            let base = cx.graph.node(core::ADD, &[x, minus_c]);
            let exponent = element_node(cx, &hyper.field, alpha, rho);
            factors.push(cx.graph.node(core::POW, &[base, exponent]));
        }
        if !hyper.exponent.0.is_empty() {
            let n = poly_node(cx, &hyper.field, alpha, x, &hyper.exponent.0);
            let d = poly_node(cx, &hyper.field, alpha, x, &hyper.exponent.1);
            let ratio = super::div(cx.graph, n, d);
            factors.push(exp_of(cx, ratio)?);
        }
        factors.push(poly_node(cx, &hyper.field, alpha, x, &hyper.poly));
        let product = mul(cx.graph, &factors);
        out.push(cx.simplify(product));
        if hyper.field.is_rational() {
            break;
        }
    }
    Some(out)
}

/// Solutions of a second-order equation `a2 y'' + a1 y' + a0 y = 0` with
/// polynomial coefficients by hyperexponential functions over a number field.
pub(super) fn second_order_solutions(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
) -> Vec<NodeId> {
    let mut fractions = Vec::new();
    for &c in a {
        let Some(f) = rational_function_in(cx.graph, c, x) else { return Vec::new() };
        fractions.push(f);
    }
    let mut op: Vec<Vec<Q>> = Vec::new();
    for (k, (numerator, _)) in fractions.iter().enumerate() {
        let mut value = numerator.clone();
        for (j, (_, denominator)) in fractions.iter().enumerate() {
            if j != k {
                value = univariate::mul(&value, denominator);
            }
        }
        while value.last().is_some_and(Zero::is_zero) {
            value.pop();
        }
        op.push(value);
    }
    for hyper in hyperexponential_over_field(&op) {
        if let Some(nodes) = solution_nodes(cx, x, &hyper)
            && nodes.len() >= 2
        {
            return nodes;
        }
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qs(v: &[i64]) -> Vec<Q> {
        v.iter().map(|&c| q(c)).collect()
    }

    #[test]
    fn irrational_exponents_at_a_regular_point() {
        // x^2 (y'' - 2 y' + y) - y = 0: solutions x^β e^x, β^2 - β - 1 = 0.
        let op = vec![qs(&[-1, 0, 1]), qs(&[0, 0, -2]), qs(&[0, 0, 1])];
        let found = hyperexponential_over_field(&op);
        assert!(!found.is_empty(), "no solution over Q(sqrt 5)");
        for h in &found {
            assert_eq!(h.field.degree(), 2);
            // P is a constant, the exponent is x.
            assert_eq!(h.poly.len(), 1);
        }
        // Both conjugates are produced.
        assert!(found.len() >= 2, "{}", found.len());
    }
}
