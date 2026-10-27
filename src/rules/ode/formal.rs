//! Formal solutions of linear ODEs with polynomial coefficients at a
//! singular point (finite, possibly algebraic, or infinity).
//!
//! The operator `L = Σ a_i(x) D^i` is written at the point `x = α` (or
//! `x = ∞`, with `t = 1/x`) in the Euler form `L = Σ_j t^j P_j(θ)`, `θ =
//! t d/dt`, with coefficients in a number field `K = Q(α)`. The Newton
//! polygon of the support `{(i, j)}` gives the slopes:
//!
//! * a segment of slope `σ = p/q > 0` is an irregular part of Poincaré rank
//!   `σ` (in `t`): the variable is ramified, `t = v^q`, so that the slope
//!   `p` is an integer; the characteristic polynomial `χ(κ) = Σ a_{ij}
//!   (-p κ)^i` over the segment gives the leading exponential factor
//!   `exp(κ v^{-p})` (its roots may be irrational or complex: `K` is
//!   extended by an irreducible factor of `χ`); the operator is conjugated
//!   by this factor (`θ → θ - p κ v^{-p}`) and the construction repeats on
//!   the part of the new polygon with smaller slope (the sub-leading terms
//!   of the exponent, Turrittin);
//! * a segment of slope `0` is regular: its indicial polynomial gives the
//!   exponents `ρ` (rational or algebraic) and the recurrence
//!   `P_{j0}(ρ+k) c_k = -Σ_{l ≥ 1} P_{j0+l}(ρ+k-l) c_{k-l}` the series.
//!
//! A formal solution is `exp(Σ κ_s v^{-s}) v^ρ Σ c_k v^k`, `v^q = t`.
//! Limits: no logarithmic terms (an exponent that is resonant with another
//! root of the same indicial polynomial is skipped); once the field is not
//! `Q` (algebraic singular point or a previous extension) further roots
//! must be of linear factors, i.e. segments of width one; the solutions
//! analytic at the point (negative slopes) are not produced.
//!
//! [`optimal_truncation`] evaluates such a series numerically: the terms are
//! summed up to the smallest one (the error is of that order, the Borel
//! sum's Stokes-free approximation).

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::rules::poly::univariate;

type Q = BigRational;

/// An element of a number field: coefficients of `1, α, α², …`.
pub type Elt = Vec<Q>;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn trim(mut p: Vec<Q>) -> Vec<Q> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

/// A number field `Q[α]/(m)`; `m = x` is the rational field itself.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Field {
    /// The monic irreducible minimal polynomial of `α` (ascending).
    pub modulus: Vec<Q>,
}

impl Field {
    /// The rationals.
    #[must_use]
    pub fn rationals() -> Self {
        Self { modulus: vec![Q::zero(), Q::one()] }
    }

    /// `Q[α]/(m)` for a monic irreducible `m`.
    #[must_use]
    pub const fn extension(modulus: Vec<Q>) -> Self {
        Self { modulus }
    }

    /// Degree over `Q`.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.modulus.len().saturating_sub(1)
    }

    /// Whether this is `Q`.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        self.degree() == 1
    }

    fn reduce(
        &self,
        p: Vec<Q>,
    ) -> Elt {
        univariate::divrem(&trim(p), &self.modulus).map_or_else(Vec::new, |(_, r)| trim(r))
    }

    /// The generator `α`.
    #[must_use]
    pub fn generator(&self) -> Elt {
        self.reduce(vec![Q::zero(), Q::one()])
    }

    /// A rational number as an element.
    #[must_use]
    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    #[allow(clippy::wrong_self_convention, reason = "the field is the context of the element")]
    pub fn from_q(
        &self,
        c: &Q,
    ) -> Elt {
        trim(vec![c.clone()])
    }

    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    fn add(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        trim(univariate::add(a, b))
    }

    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    fn sub(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        trim(univariate::sub(a, b))
    }

    fn mul(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        self.reduce(univariate::mul(a, b))
    }

    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    fn scale(
        &self,
        a: &[Q],
        c: &Q,
    ) -> Elt {
        trim(a.iter().map(|v| v * c).collect())
    }

    /// Inverse of a non-zero element.
    #[must_use]
    pub fn inverse(
        &self,
        a: &[Q],
    ) -> Option<Elt> {
        let (mut r0, mut r1) = (self.modulus.clone(), self.reduce(a.to_vec()));
        if r1.is_empty() {
            return None;
        }
        let (mut t0, mut t1): (Vec<Q>, Vec<Q>) = (Vec::new(), vec![Q::one()]);
        while !r1.is_empty() {
            let (quotient, remainder) = univariate::divrem(&r0, &r1)?;
            let t2 = trim(univariate::sub(&t0, &univariate::mul(&quotient, &t1)));
            r0 = r1;
            r1 = trim(remainder);
            t0 = t1;
            t1 = t2;
        }
        let constant_poly = trim(r0);
        let [constant] = constant_poly.as_slice() else { return None };
        Some(self.reduce(t0.iter().map(|c| c / constant).collect()))
    }

    /// The rational number an element equals, if it is one.
    #[must_use]
    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    pub fn as_rational(
        &self,
        a: &[Q],
    ) -> Option<Q> {
        match trim(a.to_vec()).as_slice() {
            | [] => Some(Q::zero()),
            | [c] => Some(c.clone()),
            | _ => None,
        }
    }

    /// Numeric value for the embedding `α ↦ alpha`.
    #[must_use]
    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    pub fn numeric(
        &self,
        a: &[Q],
        alpha: Complex64,
    ) -> Complex64 {
        a.iter().rev().fold(Complex64::new(0.0, 0.0), |acc, c| acc * alpha + c.to_f64().unwrap_or(f64::NAN))
    }

    /// The complex roots of the modulus (embeddings of `α`).
    #[must_use]
    pub fn embeddings(&self) -> Vec<Complex64> {
        let coefficients: Vec<f64> = self.modulus.iter().map(|c| c.to_f64().unwrap_or(f64::NAN)).collect();
        if self.degree() == 1 {
            return vec![Complex64::new(-coefficients.first().copied().unwrap_or(0.0), 0.0)];
        }
        crate::rules::poly::galois::durand_kerner(&coefficients)
            .unwrap_or_default()
            .into_iter()
            .map(|(re, im)| Complex64::new(re, im))
            .collect()
    }
}

// ---------------------------------------------------------------------
// Euler-form operators.

/// `Σ_j t^j P_j(θ)` with `P_j` polynomials in `θ` over `K` (ascending).
#[derive(Clone, Debug)]
struct Euler {
    terms: BTreeMap<i64, Vec<Elt>>,
}

/// `P(θ + j)` for a polynomial in `θ` over `K`.
fn shift_theta(
    field: &Field,
    p: &[Elt],
    j: i64,
) -> Vec<Elt> {
    let mut work: Vec<Elt> = p.to_vec();
    let n = work.len();
    let jq = field.from_q(&q(j));
    for i in 0..n {
        for k in (i..n.saturating_sub(1)).rev() {
            let next = work.get(k + 1).cloned().unwrap_or_default();
            let add = field.mul(&jq, &next);
            if let Some(slot) = work.get_mut(k) {
                *slot = field.add(slot, &add);
            }
        }
    }
    work
}

fn poly_add(
    field: &Field,
    a: &[Elt],
    b: &[Elt],
) -> Vec<Elt> {
    let n = a.len().max(b.len());
    let mut out: Vec<Elt> = (0..n)
        .map(|i| field.add(a.get(i).map_or(&[][..], Vec::as_slice), b.get(i).map_or(&[][..], Vec::as_slice)))
        .collect();
    while out.last().is_some_and(Vec::is_empty) {
        out.pop();
    }
    out
}

fn poly_mul(
    field: &Field,
    a: &[Elt],
    b: &[Elt],
) -> Vec<Elt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<Elt> = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i + j) {
                *slot = field.add(slot, &field.mul(x, y));
            }
        }
    }
    while out.last().is_some_and(Vec::is_empty) {
        out.pop();
    }
    out
}

impl Euler {
    const fn empty() -> Self {
        Self { terms: BTreeMap::new() }
    }

    fn add_term(
        &mut self,
        field: &Field,
        j: i64,
        poly: &[Elt],
    ) {
        let current = self.terms.remove(&j).unwrap_or_default();
        let sum = poly_add(field, &current, poly);
        if !sum.is_empty() {
            self.terms.insert(j, sum);
        }
    }

    /// Product `self · other`.
    fn compose(
        &self,
        field: &Field,
        other: &Self,
    ) -> Self {
        let mut out = Self::empty();
        for (j1, p1) in &self.terms {
            for (j2, p2) in &other.terms {
                let shifted = shift_theta(field, p1, *j2);
                out.add_term(field, j1 + j2, &poly_mul(field, &shifted, p2));
            }
        }
        out
    }

    /// The support `(i, j)` with a non-zero coefficient.
    fn support(&self) -> Vec<(i64, i64)> {
        let mut out = Vec::new();
        for (j, p) in &self.terms {
            for (i, c) in p.iter().enumerate() {
                if !c.is_empty() {
                    out.push((i64::try_from(i).unwrap_or(0), *j));
                }
            }
        }
        out
    }

    /// `t = v^r`: `θ_t = θ_v / r`, `t^j = v^{rj}`.
    fn ramify(
        &self,
        field: &Field,
        r: i64,
    ) -> Self {
        let mut out = Self::empty();
        let inverse = Q::new(BigInt::one(), BigInt::from(r));
        for (j, p) in &self.terms {
            let mut scale = Q::one();
            let scaled: Vec<Elt> = p
                .iter()
                .map(|c| {
                    let value = field.scale(c, &scale);
                    scale *= &inverse;
                    value
                })
                .collect();
            out.add_term(field, j * r, &scaled);
        }
        out
    }

    /// The same operator over a larger field (rational coefficients only).
    fn embed(
        &self,
        field: &Field,
    ) -> Self {
        let mut out = Self::empty();
        for (j, p) in &self.terms {
            out.add_term(field, *j, p);
        }
        out
    }

    /// Conjugation by `exp(κ t^{-s})`: `θ → θ - s κ t^{-s}`.
    fn conjugate(
        &self,
        field: &Field,
        kappa: &[Q],
        s: i64,
    ) -> Self {
        let max_degree = self.terms.values().map(Vec::len).max().unwrap_or(0);
        // Λ = θ - s κ t^{-s}, and its powers.
        let mut lambda = Self::empty();
        lambda.add_term(field, 0, &[Vec::new(), field.from_q(&Q::one())]);
        let c = field.scale(kappa, &-q(s));
        lambda.add_term(field, -s, &[c]);
        let mut powers = vec![{
            let mut one = Self::empty();
            one.add_term(field, 0, &[field.from_q(&Q::one())]);
            one
        }];
        for _ in 1..max_degree {
            let next = powers.last().map_or_else(Self::empty, |last| last.compose(field, &lambda));
            powers.push(next);
        }
        let mut out = Self::empty();
        for (j, p) in &self.terms {
            for (i, coefficient) in p.iter().enumerate() {
                if coefficient.is_empty() {
                    continue;
                }
                let Some(power) = powers.get(i) else { continue };
                let mut monomial = Self::empty();
                monomial.add_term(field, *j, std::slice::from_ref(coefficient));
                let product = monomial.compose(field, power);
                for (jj, pp) in &product.terms {
                    out.add_term(field, *jj, pp);
                }
            }
        }
        out
    }
}

/// Falling factorial `θ(θ-1)…(θ-i+1)`, or with `θ → -θ`.
fn falling(
    field: &Field,
    i: usize,
    negate: bool,
) -> Vec<Elt> {
    let mut poly: Vec<Elt> = vec![field.from_q(&Q::one())];
    for m in 0..i {
        let m_value = q(i64::try_from(m).unwrap_or(0));
        // (θ - m) or (-θ - m)
        let factor = if negate {
            vec![field.from_q(&-m_value), field.from_q(&q(-1))]
        } else {
            vec![field.from_q(&-m_value), field.from_q(&Q::one())]
        };
        poly = poly_mul(field, &poly, &factor);
    }
    poly
}

/// A singular point.
#[derive(Clone, Debug)]
pub enum Point {
    /// A root of the monic irreducible polynomial (ascending); linear for a
    /// rational point.
    Finite(Vec<Q>),
    /// `x = ∞`, with `t = 1/x`.
    Infinity,
}

/// The Euler form of `Σ a_i(x) D^i` at the point, over the field of the
/// point.
fn euler_form(
    op: &[Vec<Q>],
    point: &Point,
) -> (Field, Euler) {
    let mut out = Euler::empty();
    match point {
        | Point::Finite(f) => {
            let field = if f.len() == 2 {
                Field::rationals()
            } else {
                Field::extension(f.clone())
            };
            // x = α + t.
            let alpha = if f.len() == 2 { field.from_q(&-f.first().cloned().unwrap_or_else(Q::zero)) } else { field.generator() };
            for (i, a) in op.iter().enumerate() {
                // Taylor expansion of a(α + t) in t over K.
                let mut work: Vec<Elt> = a.iter().map(|c| field.from_q(c)).collect();
                let n = work.len();
                for s in 0..n {
                    for k in (s..n.saturating_sub(1)).rev() {
                        let next = work.get(k + 1).cloned().unwrap_or_default();
                        let add = field.mul(&alpha, &next);
                        if let Some(slot) = work.get_mut(k) {
                            *slot = field.add(slot, &add);
                        }
                    }
                }
                let ff = falling(&field, i, false);
                for (k, c) in work.iter().enumerate() {
                    if c.is_empty() {
                        continue;
                    }
                    let scaled: Vec<Elt> = ff.iter().map(|e| field.mul(e, c)).collect();
                    out.add_term(&field, i64::try_from(k).unwrap_or(0) - i64::try_from(i).unwrap_or(0), &scaled);
                }
            }
            (field, out)
        },
        | Point::Infinity => {
            let field = Field::rationals();
            for (i, a) in op.iter().enumerate() {
                let ff = falling(&field, i, true);
                for (k, c) in a.iter().enumerate() {
                    if c.is_zero() {
                        continue;
                    }
                    let value = field.from_q(c);
                    let scaled: Vec<Elt> = ff.iter().map(|e| field.mul(e, &value)).collect();
                    out.add_term(&field, i64::try_from(i).unwrap_or(0) - i64::try_from(k).unwrap_or(0), &scaled);
                }
            }
            (field, out)
        },
    }
}

// ---------------------------------------------------------------------
// Newton polygon.

/// A segment of the lower hull: from `(i0, j0)` to `(i1, j1)`.
#[derive(Clone, Copy, Debug)]
struct Segment {
    i0: i64,
    j0: i64,
    i1: i64,
    j1: i64,
}

impl Segment {
    fn slope(&self) -> Q {
        Q::new(BigInt::from(self.j1 - self.j0), BigInt::from(self.i1 - self.i0))
    }
}

fn lower_hull(support: &[(i64, i64)]) -> Vec<Segment> {
    let mut points: Vec<(i64, i64)> = support.to_vec();
    points.sort_unstable();
    // Smallest j for each i.
    points.dedup_by(|b, a| a.0 == b.0);
    let mut by_i: BTreeMap<i64, i64> = BTreeMap::new();
    for &(i, j) in support {
        let entry = by_i.entry(i).or_insert(j);
        *entry = (*entry).min(j);
    }
    let pts: Vec<(i64, i64)> = by_i.into_iter().collect();
    let mut hull: Vec<(i64, i64)> = Vec::new();
    for p in pts {
        while hull.len() >= 2 {
            let (a, b) = (hull[hull.len() - 2], hull[hull.len() - 1]);
            // Remove b when it lies on or above the segment a -> p.
            let cross = (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
            if cross <= 0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    hull.windows(2).map(|w| Segment { i0: w[0].0, j0: w[0].1, i1: w[1].0, j1: w[1].1 }).collect()
}

// ---------------------------------------------------------------------
// Formal solutions.

/// A formal solution `exp(Σ κ_s v^{-s}) v^ρ Σ c_k v^k` with `v^ramification = t`.
#[derive(Clone, Debug)]
pub struct FormalSolution {
    /// The coefficient field (the embedding is chosen by [`Self::alpha`]).
    pub field: Field,
    /// The ramification index `r` (`t = v^r`).
    pub ramification: usize,
    /// The exponent: pairs `(s, κ_s)` meaning `exp(κ_s v^{-s})`.
    pub exponent: Vec<(i64, Elt)>,
    /// The exponent `ρ` of `v`.
    pub rho: Elt,
    /// The series coefficients `c_k` (the part without logarithms).
    pub coefficients: Vec<Elt>,
    /// Logarithmic terms: `log_terms[k][i - 1]` is the coefficient of
    /// `(ln v)^i / i!` in `c_k` (empty when there are no logarithms).
    pub log_terms: Vec<Vec<Elt>>,
}

struct Search<'a> {
    terms: usize,
    out: &'a mut Vec<FormalSolution>,
}

#[derive(Clone)]
enum Selection {
    All,
    Below(Q),
    Exactly(Q),
}

impl Selection {
    fn accepts(
        &self,
        slope: &Q,
    ) -> bool {
        match self {
            | Self::All => true,
            | Self::Below(bound) => slope < bound,
            | Self::Exactly(value) => slope == value,
        }
    }
}

/// The formal solutions of `Σ a_i(x) D^i` (`op[i] = a_i`, ascending) at the
/// point, with `terms` series coefficients each.
#[must_use]
pub fn formal_solutions(
    op: &[Vec<Q>],
    point: &Point,
    terms: usize,
) -> Vec<FormalSolution> {
    let (field, euler) = euler_form(op, point);
    let mut out = Vec::new();
    let mut search = Search { terms, out: &mut out };
    expand(&mut search, &field, &euler, 1, &[], &Selection::All);
    out
}

fn irreducible_factors(p: &[Q]) -> Vec<(Vec<Q>, usize)> {
    let (_, factors) = univariate::factor(p);
    factors
        .into_iter()
        .map(|(f, m)| {
            let lead = Q::from_integer(f.last().cloned().unwrap_or_else(BigInt::one));
            (f.iter().map(|c| Q::from_integer(c.clone()) / &lead).collect(), m as usize)
        })
        .collect()
}

fn expand(
    search: &mut Search<'_>,
    field: &Field,
    op: &Euler,
    ramification: usize,
    exponent: &[(i64, Elt)],
    selection: &Selection,
) {
    let support = op.support();
    // The regular part: the row of lowest `j`, whose polynomial in `θ` is the
    // indicial polynomial. Roots at 0 coming from `θ^{i0}` count only when no
    // support point lies to the left of the row (else they belong to the
    // part of the polygon with negative slopes: solutions analytic there).
    if selection.accepts(&Q::zero())
        && let Some(j_min) = support.iter().map(|&(_, j)| j).min()
        && let Some(row) = op.terms.get(&j_min)
    {
        let i_row = row.iter().position(|c| !c.is_empty()).unwrap_or(0);
        let has_left = support.iter().any(|&(i, _)| i < i64::try_from(i_row).unwrap_or(0));
        regular(search, field, op, j_min, if has_left { i_row } else { 0 }, ramification, exponent);
    }
    let hull = lower_hull(&support);
    for segment in hull {
        let slope = segment.slope();
        if !selection.accepts(&slope) {
            continue;
        }
        if slope.is_negative() {
            continue;
        }
        if slope.is_zero() {
            continue;
        }
        let r = slope.denom().to_i64().unwrap_or(1);
        if r > 1 {
            // Ramify so that the slope is an integer.
            let ramified = op.ramify(field, r);
            let scaled_exponent: Vec<(i64, Elt)> = exponent.iter().map(|(s, k)| (s * r, k.clone())).collect();
            let new_slope = &slope * q(r);
            expand(
                search,
                field,
                &ramified,
                ramification * usize::try_from(r).unwrap_or(1),
                &scaled_exponent,
                &Selection::Exactly(new_slope),
            );
            continue;
        }
        let s = slope.to_integer().to_i64().unwrap_or(1);
        // Characteristic polynomial over the segment: Σ a_{ij} (-s κ)^i.
        let mut chi: Vec<Elt> = Vec::new();
        for (j, p) in &op.terms {
            for (i, c) in p.iter().enumerate() {
                let i_signed = i64::try_from(i).unwrap_or(0);
                if c.is_empty() || i_signed < segment.i0 || i_signed > segment.i1 || *j - s * i_signed != segment.j0 - s * segment.i0 {
                    continue;
                }
                while chi.len() <= i {
                    chi.push(Vec::new());
                }
                let factor = (-q(s)).pow(i32::try_from(i).unwrap_or(0));
                if let Some(slot) = chi.get_mut(i) {
                    *slot = field.add(slot, &field.scale(c, &factor));
                }
            }
        }
        // Remove the factor κ^{i0}.
        let i0 = usize::try_from(segment.i0).unwrap_or(0);
        let chi: Vec<Elt> = chi.into_iter().skip(i0).collect();
        let width = chi.len().saturating_sub(1);
        if width == 0 {
            continue;
        }
        if field.is_rational() {
            let rational: Option<Vec<Q>> = chi.iter().map(|c| field.as_rational(c)).collect();
            let Some(rational) = rational else { continue };
            for (g, _) in irreducible_factors(&rational) {
                let (new_field, kappa) = if g.len() == 2 {
                    (field.clone(), field.from_q(&(-g[0].clone())))
                } else {
                    let f = Field::extension(g);
                    let k = f.generator();
                    (f, k)
                };
                descend(search, &new_field, op, ramification, exponent, s, &kappa);
            }
        } else if width == 1 {
            // Linear: κ = -c_0 / c_1.
            let (Some(c0), Some(c1)) = (chi.first(), chi.get(1)) else { continue };
            let Some(inverse) = field.inverse(c1) else { continue };
            let kappa = field.scale(&field.mul(c0, &inverse), &q(-1));
            descend(search, field, op, ramification, exponent, s, &kappa);
        }
    }
}

fn descend(
    search: &mut Search<'_>,
    field: &Field,
    op: &Euler,
    ramification: usize,
    exponent: &[(i64, Elt)],
    s: i64,
    kappa: &[Q],
) {
    let embedded = if field.is_rational() { op.clone() } else { op.embed(field) };
    let conjugated = embedded.conjugate(field, kappa, s);
    let mut extended = exponent.to_vec();
    extended.push((s, kappa.to_vec()));
    expand(search, field, &conjugated, ramification, &extended, &Selection::Below(q(s)));
}

fn evaluate_theta(
    field: &Field,
    p: &[Elt],
    rho: &[Q],
) -> Elt {
    p.iter().rev().fold(Vec::new(), |acc, c| field.add(&field.mul(&acc, rho), c))
}

fn regular(
    search: &mut Search<'_>,
    field: &Field,
    op: &Euler,
    j0: i64,
    skip: usize,
    ramification: usize,
    exponent: &[(i64, Elt)],
) {
    let Some(indicial) = op.terms.get(&j0) else { return };
    let full: Vec<Elt> = indicial.iter().skip(skip).cloned().collect();
    // A factor θ^m: the root 0 with multiplicity m.
    let zeros = full.iter().take_while(|c| c.is_empty()).count();
    let reduced: Vec<Elt> = full.into_iter().skip(zeros).collect();
    let width = reduced.len().saturating_sub(1);
    // Roots of the indicial polynomial with multiplicities.
    let mut roots: Vec<(Field, Elt, usize)> = Vec::new();
    if zeros > 0 {
        roots.push((field.clone(), Vec::new(), zeros));
    }
    if width > 0 {
        if field.is_rational() {
            let rational: Option<Vec<Q>> = reduced.iter().map(|c| field.as_rational(c)).collect();
            let Some(rational) = rational else { return };
            for (g, multiplicity) in irreducible_factors(&rational) {
                if g.len() == 2 {
                    roots.push((field.clone(), field.from_q(&(-g[0].clone())), multiplicity));
                } else {
                    let f = Field::extension(g);
                    let rho = f.generator();
                    roots.push((f, rho, multiplicity));
                }
            }
        } else if width == 1 {
            let (Some(c0), Some(c1)) = (reduced.first(), reduced.get(1)) else { return };
            let Some(inverse) = field.inverse(c1) else { return };
            roots.push((field.clone(), field.scale(&field.mul(c0, &inverse), &q(-1)), 1));
        }
    }
    for (f, rho, multiplicity) in roots {
        let embedded = if f == *field { op.clone() } else { op.embed(&f) };
        for l0 in 0..multiplicity {
            if let Some((coefficients, log_terms)) = log_series(&f, &embedded, j0, &rho, l0, search.terms) {
                search.out.push(FormalSolution {
                    field: f.clone(),
                    ramification,
                    exponent: exponent.to_vec(),
                    rho: rho.clone(),
                    coefficients,
                    log_terms,
                });
            }
        }
    }
}

/// Taylor coefficients `P^{(n)}(a)/n!` of a polynomial in `θ` at `a`.
fn taylor_at(
    field: &Field,
    p: &[Elt],
    a: &[Q],
) -> Vec<Elt> {
    let mut work: Vec<Elt> = p.to_vec();
    let n = work.len();
    for i in 0..n {
        for k in (i..n.saturating_sub(1)).rev() {
            let next = work.get(k + 1).cloned().unwrap_or_default();
            let add = field.mul(a, &next);
            if let Some(slot) = work.get_mut(k) {
                *slot = field.add(slot, &add);
            }
        }
    }
    work
}

/// `Σ_n w_n D^n c` for `c` in the basis `L^i / i!` (`D` lowers the index).
fn apply_weights(
    field: &Field,
    weights: &[Elt],
    skip: usize,
    c: &[Elt],
) -> Vec<Elt> {
    let mut out: Vec<Elt> = vec![Vec::new(); c.len()];
    for (n, w) in weights.iter().enumerate().skip(skip) {
        if w.is_empty() {
            continue;
        }
        for (i, slot) in out.iter_mut().enumerate() {
            if let Some(ci) = c.get(i + n - skip) {
                *slot = field.add(slot, &field.mul(w, ci));
            }
        }
    }
    out
}

/// The series `t^ρ Σ_k c_k(L) t^k`, `L = ln t`, starting with
/// `c_0 = L^{l0}/l0!` (`l0` below the multiplicity of `ρ`). At a
/// resonance (`P_{j0}(ρ + k) = 0`) the equation for `c_k` is solved by raising
/// the degree in `L`, which produces the logarithmic terms. Returns the
/// non-logarithmic coefficients and the logarithmic ones.
fn log_series(
    field: &Field,
    op: &Euler,
    j0: i64,
    rho: &[Q],
    l0: usize,
    terms: usize,
) -> Option<(Vec<Elt>, Vec<Vec<Elt>>)> {
    let indicial = op.terms.get(&j0)?.clone();
    let mut c: Vec<Vec<Elt>> = Vec::new();
    let mut first: Vec<Elt> = vec![Vec::new(); l0 + 1];
    *first.get_mut(l0)? = field.from_q(&Q::one());
    c.push(first);
    for k in 1..terms {
        let a_k = field.add(rho, &field.from_q(&q(i64::try_from(k).ok()?)));
        // Right-hand side: -Σ_{l >= 1} P_{j0 + l} (ρ + k - l + D) c_{k-l}.
        let mut rhs: Vec<Elt> = Vec::new();
        for (j, p) in &op.terms {
            let l = usize::try_from(*j - j0).ok()?;
            if l == 0 || l > k {
                continue;
            }
            let a = field.sub(&a_k, &field.from_q(&q(i64::try_from(l).ok()?)));
            let weights = taylor_at(field, p, &a);
            let applied = apply_weights(field, &weights, 0, c.get(k - l)?);
            rhs = poly_add(field, &rhs, &applied);
        }
        let rhs: Vec<Elt> = rhs.iter().map(|v| field.scale(v, &q(-1))).collect();
        let weights = taylor_at(field, &indicial, &a_k);
        let mu = weights.iter().position(|w| !w.is_empty())?;
        let lead_inverse = field.inverse(weights.get(mu)?)?;
        // u = rhs shifted up by mu (the kernel constants are zero).
        let mut u: Vec<Elt> = vec![Vec::new(); mu];
        u.extend(rhs.iter().cloned());
        while u.last().is_some_and(Vec::is_empty) {
            u.pop();
        }
        // Σ_{n >= mu} w_n D^{n-mu} c = u, solved from the top index down.
        let mut ck: Vec<Elt> = vec![Vec::new(); u.len()];
        for i in (0..u.len()).rev() {
            let mut value = u.get(i).cloned().unwrap_or_default();
            for (n, w) in weights.iter().enumerate().skip(mu + 1) {
                if let Some(higher) = ck.get(i + n - mu) {
                    value = field.sub(&value, &field.mul(w, higher));
                }
            }
            *ck.get_mut(i)? = field.mul(&value, &lead_inverse);
        }
        c.push(ck);
    }
    let coefficients: Vec<Elt> = c.iter().map(|ck| ck.first().cloned().unwrap_or_default()).collect();
    let log_terms: Vec<Vec<Elt>> = c.iter().map(|ck| ck.iter().skip(1).cloned().collect()).collect();
    let any_log = log_terms.iter().any(|l| l.iter().any(|e| !e.is_empty()));
    Some((coefficients, if any_log { log_terms } else { Vec::new() }))
}

// ---------------------------------------------------------------------
// Numeric evaluation.

/// Sum of `Σ c_k v^k` truncated at the smallest term (optimal truncation):
/// the value and the size of the first omitted term as an error estimate.
#[must_use]
pub fn optimal_truncation(
    coefficients: &[Complex64],
    v: Complex64,
) -> (Complex64, f64) {
    let mut sum = Complex64::new(0.0, 0.0);
    let mut power = Complex64::new(1.0, 0.0);
    let mut previous = f64::INFINITY;
    let mut error = f64::INFINITY;
    for (k, c) in coefficients.iter().enumerate() {
        let term = c * power;
        let size = term.norm();
        // Stop when the terms start to grow (past the first few).
        if k > 2 && size > previous && size > 0.0 {
            error = previous;
            break;
        }
        sum += term;
        if size > 0.0 {
            previous = size;
        }
        error = size;
        power *= v;
    }
    (sum, error)
}

impl FormalSolution {
    /// The solution at `x` for the embedding `alpha` of the field's
    /// generator, with the singular point `point_value` (`α` itself for a
    /// finite point; ignored for infinity, `infinity = true`). The series is
    /// summed with [`optimal_truncation`]. Returns `(value, error)`.
    #[must_use]
    pub fn evaluate(
        &self,
        alpha: Complex64,
        t: Complex64,
    ) -> (Complex64, f64) {
        let r = self.ramification as f64;
        // v = t^(1/r), principal branch.
        let v = (t.ln() / r).exp();
        let mut exponent_sum = Complex64::new(0.0, 0.0);
        for (s, kappa) in &self.exponent {
            exponent_sum += self.field.numeric(kappa, alpha) * v.powi(-i32::try_from(*s).unwrap_or(0));
        }
        let rho = self.field.numeric(&self.rho, alpha);
        let coefficients = self.numeric_coefficients(alpha, v);
        let (series_value, error) = optimal_truncation(&coefficients, v);
        let prefactor = exponent_sum.exp() * (rho * v.ln()).exp();
        (prefactor * series_value, error * prefactor.norm())
    }
}

impl FormalSolution {
    /// Whether the series has logarithmic terms.
    #[must_use]
    pub fn has_logarithms(&self) -> bool {
        self.log_terms.iter().any(|l| l.iter().any(|e| !e.is_empty()))
    }

    /// The coefficients `c_k` as complex numbers, the logarithmic parts
    /// evaluated at `ln v`.
    fn numeric_coefficients(
        &self,
        alpha: Complex64,
        v: Complex64,
    ) -> Vec<Complex64> {
        let log = v.ln();
        self.coefficients
            .iter()
            .enumerate()
            .map(|(k, c)| {
                let mut value = self.field.numeric(c, alpha);
                let mut power = Complex64::new(1.0, 0.0);
                let mut factorial = 1.0_f64;
                for (i, lc) in self.log_terms.get(k).map_or(&[][..], Vec::as_slice).iter().enumerate() {
                    power *= log;
                    factorial *= (i + 1) as f64;
                    value += self.field.numeric(lc, alpha) * power / factorial;
                }
                value
            })
            .collect()
    }

    /// Like [`Self::evaluate`] with the series summed by Borel–Padé
    /// resummation ([`super::resum::borel_pade`]): `(value, error, stokes)`.
    #[must_use]
    pub fn evaluate_borel(
        &self,
        alpha: Complex64,
        t: Complex64,
    ) -> Option<(Complex64, f64, bool)> {
        let r = self.ramification as f64;
        let v = (t.ln() / r).exp();
        let mut exponent_sum = Complex64::new(0.0, 0.0);
        for (s, kappa) in &self.exponent {
            exponent_sum += self.field.numeric(kappa, alpha) * v.powi(-i32::try_from(*s).unwrap_or(0));
        }
        let rho = self.field.numeric(&self.rho, alpha);
        let coefficients = self.numeric_coefficients(alpha, v);
        let sum = super::resum::borel_pade(&coefficients, v)?;
        let prefactor = exponent_sum.exp() * (rho * v.ln()).exp();
        Some((prefactor * sum.value, prefactor.norm() * sum.error, sum.stokes))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qs(v: &[i64]) -> Vec<Q> {
        v.iter().map(|&c| q(c)).collect()
    }

    fn rational(a: &Elt) -> Q {
        a.first().cloned().unwrap_or_else(Q::zero)
    }

    #[test]
    fn airy_at_infinity_is_ramified() {
        // y'' - x y = 0: exp(∓(2/3) x^(3/2)) x^(-1/4) Σ c_k x^(-3k/2).
        let solutions = formal_solutions(&[qs(&[0, -1]), qs(&[]), qs(&[1])], &Point::Infinity, 30);
        assert_eq!(solutions.len(), 2, "{solutions:?}");
        for s in &solutions {
            assert_eq!(s.ramification, 2);
            let [(power, kappa)] = s.exponent.as_slice() else { panic!("{s:?}") };
            assert_eq!(*power, 3);
            assert_eq!(rational(kappa).abs(), Q::new(BigInt::from(2), BigInt::from(3)));
            // v^(1/2): t^(1/4) = x^(-1/4).
            assert_eq!(rational(&s.rho), Q::new(BigInt::one(), BigInt::from(2)));
        }
        // The decaying solution has Ai's coefficients: c_3 = -5/48 c_0 (ζ = (2/3) x^(3/2)).
        let decaying = solutions.iter().find(|s| rational(&s.exponent[0].1).is_negative()).expect("decaying");
        assert_eq!(rational(&decaying.coefficients[1]), Q::zero());
        assert_eq!(rational(&decaying.coefficients[3]), Q::new(BigInt::from(-5), BigInt::from(48)));
        // Optimal truncation reproduces Ai up to the normalisation 1/(2 sqrt(pi)):
        // Ai(5) = 1.0834e-4.
        let (value, error) = decaying.evaluate(Complex64::new(0.0, 0.0), Complex64::new(1.0 / 5.0, 0.0));
        let ai = value.re / (2.0 * std::f64::consts::PI.sqrt());
        assert!((ai - 1.083_444_3e-4).abs() < 1e-8, "{ai} (error estimate {error})");
    }

    #[test]
    fn bessel_has_complex_exponents() {
        // x^2 y'' + x y' + (x^2 - 1/4) y = 0: exp(±i x) x^(-1/2) (1 + ...).
        let solutions = formal_solutions(
            &[vec![Q::new(BigInt::from(-1), BigInt::from(4)), Q::zero(), Q::one()], vec![Q::zero(), Q::one()], vec![Q::zero(), Q::zero(), Q::one()]],
            &Point::Infinity,
            8,
        );
        assert_eq!(solutions.len(), 1, "{solutions:?}");
        assert_eq!(solutions[0].field.degree(), 2);
        for s in &solutions {
            assert_eq!(s.ramification, 1);
            assert_eq!(s.field.modulus, qs(&[1, 0, 1]));
            assert_eq!(rational(&s.rho), Q::new(BigInt::one(), BigInt::from(2)));
        }
        // For nu = 1/2: J = sin(x)/sqrt(x), so the series terminates (c_k = 0 for k >= 1).
        assert!(solutions[0].coefficients.iter().skip(1).all(Vec::is_empty), "{:?}", solutions[0].coefficients);
    }

    #[test]
    fn irrational_exponential_part() {
        // y'' - 2 y = 0: exp(±sqrt(2) x).
        let solutions = formal_solutions(&[qs(&[-2]), qs(&[]), qs(&[1])], &Point::Infinity, 4);
        assert_eq!(solutions.len(), 1);
        assert_eq!(solutions[0].field.degree(), 2);
        for s in &solutions {
            assert_eq!(s.field.modulus, qs(&[-2, 0, 1]));
            assert_eq!(s.exponent.first().map(|e| e.0), Some(1));
        }
    }

    #[test]
    fn irregular_point_of_rank_two_with_subleading_terms() {
        // x^4 y'' = y at 0: y = exp(±1/x)... with sub-leading exponent terms for
        // x^3 y'' + ... : here y'' = y / x^4 has exp(±1/x) x (1 + ...).
        let solutions = formal_solutions(&[qs(&[-1]), qs(&[]), qs(&[0, 0, 0, 0, 1])], &Point::Finite(qs(&[0, 1])), 6);
        assert_eq!(solutions.len(), 2, "{solutions:?}");
        for s in &solutions {
            assert_eq!(s.exponent.len(), 1);
            assert_eq!(s.exponent[0].0, 1);
            assert_eq!(rational(&s.exponent[0].1).abs(), Q::one());
            assert_eq!(rational(&s.rho), q(1));
        }
    }

    #[test]
    fn regular_point_at_an_algebraic_singularity() {
        // (x^2 - 2) y'' + x y' = 0 has y' ∝ (x^2-2)^(-1/2): exponents 0 and 1/2 at sqrt(2).
        let solutions = formal_solutions(&[qs(&[]), qs(&[0, 1]), qs(&[-2, 0, 1])], &Point::Finite(qs(&[-2, 0, 1])), 6);
        let half = Some(Q::new(BigInt::one(), BigInt::from(2)));
        assert!(solutions.iter().any(|s| s.field.as_rational(&s.rho) == half), "{solutions:?}");
    }

    #[test]
    fn logarithms_at_double_and_resonant_exponents() {
        // Bessel of order 0 at the origin: exponent 0 twice.
        let op = [qs(&[0, 0, 1]), qs(&[0, 1]), qs(&[0, 0, 1])];
        let solutions = formal_solutions(&op, &Point::Finite(qs(&[0, 1])), 8);
        assert_eq!(solutions.len(), 2, "{solutions:?}");
        let plain = solutions.iter().find(|s| !s.has_logarithms()).expect("J0");
        assert_eq!(rational(&plain.coefficients[2]), Q::new(BigInt::from(-1), BigInt::from(4)));
        assert_eq!(rational(&plain.coefficients[4]), Q::new(BigInt::one(), BigInt::from(64)));
        // y2 = J0 ln x + x^2/4 - 3 x^4/128 + ...
        let logged = solutions.iter().find(|s| s.has_logarithms()).expect("log solution");
        assert_eq!(rational(&logged.log_terms[0][0]), Q::one());
        assert_eq!(rational(&logged.log_terms[2][0]), Q::new(BigInt::from(-1), BigInt::from(4)));
        assert_eq!(rational(&logged.coefficients[2]), Q::new(BigInt::one(), BigInt::from(4)));
        assert_eq!(rational(&logged.coefficients[4]), Q::new(BigInt::from(-3), BigInt::from(128)));
        // x y'' + y = 0: exponents 0 and 1 differ by an integer, the lower solution has a logarithm
        // 1 - x ln x + ...
        let solutions = formal_solutions(&[qs(&[1]), qs(&[]), qs(&[0, 1])], &Point::Finite(qs(&[0, 1])), 6);
        assert_eq!(solutions.len(), 2, "{solutions:?}");
        let upper = solutions.iter().find(|s| rational(&s.rho) == q(1)).expect("rho = 1");
        assert!(!upper.has_logarithms());
        assert_eq!(rational(&upper.coefficients[1]), Q::new(BigInt::from(-1), BigInt::from(2)));
        let lower = solutions.iter().find(|s| s.rho.is_empty()).expect("rho = 0");
        assert_eq!(rational(&lower.log_terms[1][0]), q(-1));
        // Numerically: the log solution of Bessel at x = 0.5 against y2 = J0(x) ln x + (x^2/4 - 3 x^4/128 + ...).
        let (value, _) = logged.evaluate(Complex64::new(0.0, 0.0), Complex64::new(0.5, 0.0));
        let j0 = 0.938_469_807_240_813;
        let expected = j0 * 0.5_f64.ln() + 0.0625 - 3.0 * 0.5_f64.powi(4) / 128.0 + 11.0 * 0.5_f64.powi(6) / 13824.0;
        assert!((value.re - expected).abs() < 2e-5, "{value} vs {expected}");
    }

    #[test]
    fn newton_polygon_slopes() {
        // y' = y / x^2 at 0: slope 1.
        let hull = lower_hull(&[(0, -1), (1, 0)]);
        assert_eq!(hull.len(), 1);
        assert_eq!(hull[0].slope(), q(1));
    }
}
