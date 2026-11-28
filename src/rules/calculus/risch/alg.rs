//! The Risch–Bronstein integration algorithm over a tower of monomial
//! extensions.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use std::collections::BTreeMap;

use num_bigint::BigInt;

use super::Q;
use super::elem::Coeff;
use super::elem::Elem;
use super::flat;
use super::flat::Flat;
use super::flat::Mono;
use super::poly;
use super::solve;
use super::tower::Kind;
use super::tower::Tower;
use crate::rules::poly::univariate;

/// Widest residue polynomial whose rational roots are looked for.
const MAX_RESIDUE_DEGREE: usize = 32;

/// The integral of an element: a rational part plus a sum of logarithms.
#[derive(Clone, Debug)]
pub struct Integral {
    pub rational: Elem,
    pub logs: Vec<(Q, Elem)>,
}
/// `∫ f dx`, or `None` when the tower algorithm does not decide it.
pub fn integrate(
    f: &Elem,
    tower: &Tower,
) -> Option<Integral> {
    if f.depth() != tower.depth() || tower.depth() == 0 {
        return None;
    }
    integrate_level(f, tower.depth(), tower)
}

fn integrate_level(
    f: &Elem,
    level: usize,
    tower: &Tower,
) -> Option<Integral> {
    let mut current = f.clone();
    let mut extra: Vec<(i64, Elem)> = Vec::new();
    if tower.kind(level) == Kind::Exp
        && let Some((principal, rest)) = split_principal(&current, level)?
    {
        extra.extend(principal);
        current = rest;
    }
    let (a, d) = current.parts()?;
    let (q, r) = poly::divrem(a, d)?;
    let (g, r, d) = hermite(&r, d, level, tower)?;
    current = Elem::ratio(level, r, d)?;
    let mut logs: Vec<(Q, Elem)> = Vec::new();
    if !current.is_zero() {
        let (found, rest) = log_part(&current, level, tower)?;
        logs.extend(found);
        current = rest;
    }
    if !current.is_polynomial() {
        return None;
    }
    let mut coefficients: BTreeMap<i64, Elem> = BTreeMap::new();
    let mut push = |index: i64, value: &Elem| {
        let slot = coefficients
            .entry(index)
            .or_insert_with(|| Elem::zero_at(level - 1));
        *slot = slot.add(value);
    };
    for (i, c) in q.iter().enumerate() {
        push(i64::try_from(i).ok()?, c);
    }
    if let Some(parts) = current.polynomial_parts() {
        for (i, c) in parts.iter().enumerate() {
            push(i64::try_from(i).ok()?, c);
        }
    }
    for (index, c) in &extra {
        push(*index, c);
    }
    coefficients.retain(|_, c| !c.is_zero());
    let target = from_coefficients(&coefficients, level)?;
    let polynomial = poly_part(&coefficients, level, tower)?;
    if !polynomial.deriv(tower).sub(&target).is_zero() {
        return None;
    }
    Some(Integral {
        rational: g.add(&polynomial),
        logs,
    })
}

/// Builds `sum c_i theta^i` for a (possibly negative) exponent map.
fn from_coefficients(
    coefficients: &BTreeMap<i64, Elem>,
    level: usize,
) -> Option<Elem> {
    let theta = Elem::gen_at(level);
    let mut out = Elem::zero_at(level);
    for (&index, c) in coefficients {
        if c.is_zero() {
            continue;
        }
        out = out.add(&c.mul(&theta.pow_i(index)?));
    }
    Some(out)
}

fn monomial_poly(
    degree: usize,
    level: usize,
) -> Vec<Elem> {
    let mut out = vec![Elem::zero_at(level - 1); degree];
    out.push(Elem::one_at(level - 1));
    out
}

/// `r/d = g' + r'/d'` with `d'` squarefree.
fn hermite(
    a: &[Elem],
    d: &[Elem],
    level: usize,
    tower: &Tower,
) -> Option<(Elem, Vec<Elem>, Vec<Elem>)> {
    let d = poly::monic(d);
    if poly::degree(&d).is_none_or(|degree| degree == 0) {
        return Some((
            Elem::zero_at(level),
            Vec::new(),
            vec![Elem::one_at(level - 1)],
        ));
    }
    if poly::is_zero(a) {
        return Some((Elem::zero_at(level), Vec::new(), d));
    }
    let (slope, constant) = tower.monomial(level);
    let dd = poly::derivative_x(&d, &slope, &constant, tower);
    let common = poly::gcd(&d, &dd);
    if poly::degree(&common).is_none_or(|degree| degree == 0) {
        return Some((Elem::zero_at(level), poly::trim(a.to_vec()), d));
    }
    let (d_sqf, d2) = poly::squarefree_parts(&d);
    let n_b = poly::degree(&d2)?;
    let n_c = poly::degree(&d_sqf)?;
    let dd2 = poly::derivative_x(&d2, &slope, &constant, tower);
    let mut columns: Vec<Vec<Elem>> = Vec::with_capacity(n_b + n_c);
    for i in 0..n_b {
        let xi = monomial_poly(i, level);
        let dxi = poly::derivative_x(&xi, &slope, &constant, tower);
        let first = poly::mul(&poly::mul(&dxi, &d2), &d_sqf);
        let second = poly::mul(&poly::mul(&xi, &dd2), &d_sqf);
        columns.push(poly::sub(&first, &second));
    }
    for j in 0..n_c {
        let xj = monomial_poly(j, level);
        columns.push(poly::mul(&xj, &poly::mul(&d2, &d2)));
    }
    let rhs = poly::mul(a, &d2);
    let rows = columns
        .iter()
        .map(Vec::len)
        .max()
        .unwrap_or(0)
        .max(rhs.len());
    if rows == 0 {
        return None;
    }
    let mut matrix = vec![vec![Elem::zero_at(level - 1); columns.len()]; rows];
    for (j, column) in columns.iter().enumerate() {
        for (i, c) in column.iter().enumerate() {
            matrix[i][j] = c.clone();
        }
    }
    let mut values = vec![Elem::zero_at(level - 1); rows];
    for (i, c) in rhs.iter().enumerate() {
        values[i] = c.clone();
    }
    let solution = solve::solve_field(&mut matrix, &mut values, columns.len())?;
    let numerator = solution.get(..n_b)?.to_vec();
    let remainder = solution.get(n_b..)?.to_vec();
    let g = Elem::ratio(level, numerator, d2.clone())?;
    Some((g, poly::trim(remainder), d_sqf))
}

type PrincipalPart = (Vec<(i64, Elem)>, Elem);

/// Splits the principal part at `theta = 0` off an element of an
/// exponential extension: `f = sum c_k theta^{-k} + rest` with the
/// denominator of `rest` coprime to `theta`.
fn split_principal(
    f: &Elem,
    level: usize,
) -> Option<Option<PrincipalPart>> {
    let (a, d) = f.parts()?;
    let m = d.iter().take_while(|c| c.is_zero()).count();
    if m == 0 {
        return Some(None);
    }
    let d0 = d.get(m..)?.to_vec();
    let zero = Elem::zero_at(level - 1);
    let d0_0 = d0.first().cloned().unwrap_or_else(|| zero.clone());
    if d0_0.is_zero() {
        return None;
    }
    let inv0 = d0_0.inv()?;
    // `d0^{-1} mod theta^m`
    let mut inverse = vec![zero.clone(); m];
    inverse[0] = inv0.clone();
    for j in 1..m {
        let mut acc = zero.clone();
        for i in 1..=j {
            if let Some(di) = d0.get(i) {
                acc = acc.add(&di.mul(&inverse[j - i]));
            }
        }
        inverse[j] = acc.neg().mul(&inv0);
    }
    // `a * d0^{-1} mod theta^m`
    let mut series = vec![zero.clone(); m];
    for i in 0..m.min(a.len()) {
        for j in 0..m - i {
            series[i + j] = series[i + j].add(&a[i].mul(&inverse[j]));
        }
    }
    let mut principal = Vec::with_capacity(m);
    let mut rest = f.clone();
    let theta = Elem::gen_at(level);
    for (j, c) in series.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let power = i64::try_from(m - j).ok()?;
        rest = rest.sub(&c.mul(&theta.pow_i(-power)?));
        principal.push((-power, c.clone()));
    }
    Some(Some((principal, rest)))
}

/// The logarithmic part of `f = a/d` with `d` squarefree.
fn log_part(
    f: &Elem,
    level: usize,
    tower: &Tower,
) -> Option<(Vec<(Q, Elem)>, Elem)> {
    let (a, d) = f.parts()?;
    let d = poly::trim(d.to_vec());
    if d.len() <= 1 {
        return Some((Vec::new(), f.clone()));
    }
    let dd = poly::formal_derivative(&d);
    if poly::degree(&poly::gcd(&d, &dd)).is_none_or(|degree| degree != 0) {
        return None;
    }
    let (slope, constant) = tower.monomial(level);
    let theta_derivative = vec![constant, slope];
    let e = poly::mul(&dd, &theta_derivative);
    let m = poly::degree(&d)?;
    let mut values = Vec::with_capacity(m + 1);
    for k in 0..=m {
        let kq = Q::from_integer(BigInt::from(k));
        let shifted = poly::sub(a, &poly::scale(&e, &Elem::from_q(&kq)));
        if poly::is_zero(&shifted) {
            values.push(Elem::zero_at(level - 1));
        } else {
            values.push(poly::resultant(&d, &shifted)?);
        }
    }
    let interpolated = interpolate(&values)?;
    if interpolated.iter().all(Coeff::is_zero) {
        return Some((Vec::new(), f.clone()));
    }
    let roots = rational_roots_of(&interpolated)?;
    if roots.is_empty() {
        return None;
    }
    let mut logs = Vec::new();
    let mut out = f.clone();
    for c in roots {
        let shifted = poly::sub(a, &poly::scale(&e, &Elem::from_q(&c)));
        let factor = poly::gcd(&d, &shifted);
        if poly::degree(&factor).is_none_or(|degree| degree == 0) {
            continue;
        }
        let factor = poly::monic(&factor);
        let argument = Elem::ratio(level, factor, vec![Elem::one_at(level - 1)])?;
        let derivative = argument.deriv(tower).mul(&argument.inv()?);
        out = out.sub(&Elem::from_q(&c).mul(&derivative));
        logs.push((c, argument));
    }
    Some((logs, out))
}

/// Rational roots of a polynomial over the coefficient field, found after
/// clearing the denominators of its coefficients.
fn rational_roots_of(coefficients: &[Elem]) -> Option<Vec<Q>> {
    // The rational root finder factors the polynomial; a wide resultant is
    // left to the partial-fraction paths instead.
    if coefficients.len() > MAX_RESIDUE_DEGREE {
        return None;
    }
    if let Some(rational) = coefficients
        .iter()
        .map(Coeff::as_q)
        .collect::<Option<Vec<Q>>>()
    {
        return Some(univariate::rational_roots(&rational));
    }
    let mut denominators: Vec<Flat> = Vec::new();
    let mut fractions: Vec<(Flat, usize)> = Vec::with_capacity(coefficients.len());
    for c in coefficients {
        let (num, den) = c.flatten();
        let index = match denominators.iter().position(|d| d == &den) {
            | Some(i) => i,
            | None => {
                denominators.push(den);
                denominators.len() - 1
            },
        };
        fractions.push((num, index));
    }
    if denominators.len() > 8 {
        return None;
    }
    let mut complement: Vec<Flat> = Vec::with_capacity(denominators.len());
    for skip in 0..denominators.len() {
        let mut product = Flat::one();
        for (j, d) in denominators.iter().enumerate() {
            if j != skip {
                product = product.mul(d);
            }
        }
        complement.push(product);
    }
    let mut by_monomial: BTreeMap<Mono, Vec<Q>> = BTreeMap::new();
    for (i, (num, index)) in fractions.iter().enumerate() {
        let scaled = num.mul(&complement[*index]);
        for (mono, c) in scaled.terms() {
            let row = by_monomial
                .entry(mono.clone())
                .or_insert_with(|| vec![Q::from_integer(BigInt::from(0)); coefficients.len()]);
            row[i] += c;
        }
    }
    if by_monomial.len() > MAX_RESIDUE_DEGREE {
        return None;
    }
    let mut common: Vec<Q> = Vec::new();
    for (_, poly) in by_monomial {
        common = if common.is_empty() {
            poly
        } else {
            univariate::gcd(&common, &poly)
        };
        if common.len() <= 1 {
            return Some(Vec::new());
        }
    }
    Some(univariate::rational_roots(&common))
}

/// Lagrange interpolation at the nodes `0, 1, ..., n`.
fn interpolate(values: &[Elem]) -> Option<Vec<Elem>> {
    if values.is_empty() {
        return None;
    }
    let m = values.len() - 1;
    let mut out = vec![Elem::zero(); m + 1];
    for k in 0..=m {
        let mut basis: Vec<Q> = vec![Q::from_integer(BigInt::from(1))];
        let mut denominator = Q::from_integer(BigInt::from(1));
        for j in 0..=m {
            if j == k {
                continue;
            }
            let node = Q::from_integer(BigInt::from(j));
            let mut next = vec![Q::from_integer(BigInt::from(0)); basis.len() + 1];
            for (i, c) in basis.iter().enumerate() {
                next[i + 1] += c;
                next[i] -= c * &node;
            }
            basis = next;
            denominator *= Q::from_integer(BigInt::from(k)) - node;
        }
        let inv = denominator.recip();
        for (i, c) in basis.iter().enumerate() {
            if num_traits::Zero::is_zero(c) {
                continue;
            }
            let scale = c * &inv;
            out[i] = out[i].add(&values[k].mul(&Elem::Const(scale)));
        }
    }
    Some(out)
}

/// Builds the monomial `prod theta_v^{e_v}`.
fn monomial_elem(
    exponents: &[i64],
    tower: &Tower,
) -> Option<Elem> {
    let mut out = Elem::one_at(tower.depth());
    for (v, &e) in exponents.iter().enumerate() {
        if e == 0 {
            continue;
        }
        out = out.mul(&Elem::gen_at(v + 1).pow_i(e)?);
    }
    Some(out)
}

fn rows_estimate(columns: &[Flat]) -> usize {
    columns.iter().map(Flat::len).sum()
}

/// Solves `p' = h` for a polynomial `p` in the top generator, by a bounded
/// ansatz over the monomials of the tower and one linear system over `Q`.
fn poly_part(
    coefficients: &BTreeMap<i64, Elem>,
    level: usize,
    tower: &Tower,
) -> Option<Elem> {
    poly_part_with(coefficients, level, tower, false)
        .or_else(|| poly_part_with(coefficients, level, tower, true))
}

fn poly_part_with(
    coefficients: &BTreeMap<i64, Elem>,
    level: usize,
    tower: &Tower,
    signed: bool,
) -> Option<Elem> {
    if coefficients.is_empty() {
        return Some(Elem::zero_at(level));
    }
    let h = from_coefficients(coefficients, level)?;
    let (h_num, h_den) = h.flatten();
    let vars = level;
    let mut max_exp = vec![0u32; vars];
    for mono in h_num
        .terms()
        .map(|(m, _)| m)
        .chain(h_den.terms().map(|(m, _)| m))
    {
        for &(v, e) in mono {
            if let Some(slot) = max_exp.get_mut(usize::from(v)) {
                *slot = (*slot).max(e);
            }
        }
    }
    let lowest = coefficients.keys().copied().min().unwrap_or(0);
    let highest = coefficients.keys().copied().max().unwrap_or(0);
    let mut ranges: Vec<(i64, i64)> = Vec::with_capacity(vars);
    for v in 0..vars {
        let margin = i64::from(max_exp[v]) + 1;
        if v == 0 {
            if signed {
                ranges.push((-margin, margin));
            } else {
                ranges.push((0, margin));
            }
        } else if v + 1 == level && tower.kind(level) == Kind::Exp {
            ranges.push((lowest.min(0), highest.max(0)));
        } else if v + 1 == level {
            ranges.push((0, highest.max(0) + 1));
        } else {
            ranges.push((0, margin));
        }
    }
    let mut total: usize = 1;
    for &(lo, hi) in &ranges {
        total = total.checked_mul(usize::try_from(hi - lo + 1).ok()?)?;
    }
    if total == 0 || total > 160 {
        return None;
    }
    let mut exponents: Vec<Vec<i64>> = vec![Vec::new()];
    for &(lo, hi) in &ranges {
        let mut next = Vec::with_capacity(exponents.len() * usize::try_from(hi - lo + 1).ok()?);
        for prefix in &exponents {
            for e in lo..=hi {
                let mut v = prefix.clone();
                v.push(e);
                next.push(v);
            }
        }
        exponents = next;
    }
    let mut denominators: Vec<Flat> = Vec::new();
    let index_of = |denominators: &mut Vec<Flat>, value: &Flat| -> usize {
        match denominators.iter().position(|d| d == value) {
            | Some(i) => i,
            | None => {
                denominators.push(value.clone());
                denominators.len() - 1
            },
        }
    };
    let mut columns: Vec<Flat> = Vec::with_capacity(exponents.len());
    let mut monomials: Vec<Elem> = Vec::with_capacity(exponents.len());
    let mut column_den: Vec<usize> = Vec::with_capacity(exponents.len());
    for exps in &exponents {
        let m = monomial_elem(exps, tower)?;
        let (num, den) = m.deriv(tower).flatten();
        let idx = index_of(&mut denominators, &den);
        columns.push(num);
        column_den.push(idx);
        monomials.push(m);
    }
    let h_idx = index_of(&mut denominators, &h_den);
    if denominators.len() > 8 {
        return None;
    }
    let mut complement: Vec<Flat> = Vec::with_capacity(denominators.len());
    for skip in 0..denominators.len() {
        let mut product = Flat::one();
        for (j, d) in denominators.iter().enumerate() {
            if j != skip {
                product = product.mul(d);
            }
        }
        complement.push(product);
    }
    let scaled: Vec<Flat> = columns
        .iter()
        .zip(&column_den)
        .map(|(c, &i)| c.mul(&complement[i]))
        .collect();
    let target = h_num.mul(&complement[h_idx]);
    if rows_estimate(&scaled) + target.len() > 4000 {
        return None;
    }
    let mut rows: BTreeMap<Mono, usize> = BTreeMap::new();
    for column in &scaled {
        for (mono, _) in column.terms() {
            let next = rows.len();
            rows.entry(mono.clone()).or_insert(next);
        }
    }
    for (mono, _) in target.terms() {
        let next = rows.len();
        rows.entry(mono.clone()).or_insert(next);
    }
    let mut matrix = vec![vec![Q::from_integer(BigInt::from(0)); scaled.len()]; rows.len()];
    for (j, column) in scaled.iter().enumerate() {
        for (mono, c) in column.terms() {
            let i = *rows.get(mono)?;
            matrix[i][j] = c.clone();
        }
    }
    let mut rhs = vec![Q::from_integer(BigInt::from(0)); rows.len()];
    for (mono, c) in target.terms() {
        let i = *rows.get(mono)?;
        rhs[i] = c.clone();
    }
    let solution = match flat::solve_q(&mut matrix, &mut rhs, scaled.len()) {
        | Some(v) => v,
        | None => {
            eprintln!(
                "DBG poly_part: singular, rows={} cols={}",
                rows.len(),
                scaled.len()
            );
            return None;
        },
    };
    let mut out = Elem::zero_at(level);
    for (c, m) in solution.iter().zip(&monomials) {
        if num_traits::Zero::is_zero(c) {
            continue;
        }
        out = out.add(&m.mul(&Elem::from_q(c)));
    }
    Some(out)
}
