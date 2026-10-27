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
//! When a characteristic or indicial polynomial over `K = Q(α)` has a root
//! outside `K`, the field is extended (a tower step, [`tower_roots`]): the
//! norm of `β + k α` over `Q` is computed numerically and rationalised, an
//! irreducible factor of it defines `K' = Q[y]/F`, and `α`, `β` are recovered
//! as polynomials in the generator by lattice reduction and verified exactly;
//! the operator and the exponent list are then carried to `K'`.
//! Limits: only the first root of a non-split indicial factor is followed
//! (no logarithms there); the solutions analytic at the point (negative
//! slopes) are not produced.
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

use crate::rules::poly::galois::Fx;
use crate::rules::poly::galois::BITS;
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
    pub(super) fn add(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        trim(univariate::add(a, b))
    }

    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    pub(super) fn sub(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        trim(univariate::sub(a, b))
    }

    pub(super) fn mul(
        &self,
        a: &[Q],
        b: &[Q],
    ) -> Elt {
        self.reduce(univariate::mul(a, b))
    }

    #[allow(clippy::unused_self, reason = "kept as methods for a uniform field interface")]
    pub(super) fn scale(
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

    /// The same operator over an extension `K'` of `K`, with `α ↦ h`.
    fn map_field(
        &self,
        from: &Field,
        to: &Field,
        h: &[Q],
    ) -> Self {
        let _ = from;
        let mut out = Self::empty();
        for (j, p) in &self.terms {
            let mapped: Vec<Elt> = p.iter().map(|c| map_element(to, c, h)).collect();
            out.add_term(to, *j, &mapped);
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
    // indicial polynomial; all its roots, those at 0 included, are exponents
    // (the part of the polygon with negative slopes describes the same
    // solutions).
    if selection.accepts(&Q::zero())
        && let Some(j_min) = support.iter().map(|&(_, j)| j).min()
    {
        regular(search, field, op, j_min, 0, ramification, exponent);
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
        } else {
            // Over a number field: the roots that lie in the field (a linear
            // polynomial has its one root there).
            for kappa in roots_in_field(field, &chi) {
                descend(search, field, op, ramification, exponent, s, &kappa);
            }
            // Roots in a proper extension of the field: one more step of the tower.
            for extension in tower_roots(field, &chi) {
                let mapped_op = op.map_field(field, &extension.field, &extension.alpha);
                let mapped_exponent: Vec<(i64, Elt)> =
                    exponent.iter().map(|(t, k)| (*t, map_element(&extension.field, k, &extension.alpha))).collect();
                descend(search, &extension.field, &mapped_op, ramification, &mapped_exponent, s, &extension.root);
            }
        }
    }
}

/// Newton refinement of a root of the fixed-point polynomial `p` (ascending).
fn refine_root(
    p: &[Fx],
    start: Complex64,
) -> Option<Fx> {
    let derivative: Vec<Fx> = p
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.mul(&Fx::integer(&BigInt::from(i))))
        .collect();
    let value = |poly: &[Fx], z: &Fx| poly.iter().rev().fold(Fx::zero(), |acc, c| acc.mul(z).add(c));
    let mut z = Fx::from_f64(start.re, start.im)?;
    for _ in 0..14 {
        let step = value(p, &z).div(&value(&derivative, &z))?;
        z = z.sub(&step);
    }
    Some(z)
}

/// The roots in `K = Q(α)` of a polynomial over `K` (ascending), found
/// numerically at 190 digits for one embedding of `α` and recognised as
/// rational combinations of `1, α, …, α^{d-1}` by lattice reduction, then
/// checked exactly. Roots outside `K` (they would need a larger field) are
/// not returned.
pub(super) fn roots_in_field(
    field: &Field,
    poly: &[Elt],
) -> Vec<Elt> {
    let degree = poly.len().saturating_sub(1);
    if degree == 0 {
        return Vec::new();
    }
    if field.is_rational() {
        let Some(rational): Option<Vec<Q>> = poly.iter().map(|c| field.as_rational(c)).collect() else { return Vec::new() };
        return irreducible_factors(&rational).into_iter().filter(|(g, _)| g.len() == 2).map(|(g, _)| field.from_q(&-g[0].clone())).collect();
    }
    let d = field.degree();
    let to_fx = |q: &Q| Fx::real((q.numer() << BITS) / q.denom());
    // The embedding α ≈ alpha0.
    let modulus: Vec<Fx> = field.modulus.iter().map(to_fx).collect();
    let Some(alpha0) = field.embeddings().first().copied().and_then(|a| refine_root(&modulus, a)) else { return Vec::new() };
    let mut powers = vec![Fx::integer(&BigInt::one())];
    for _ in 1..d {
        let next = powers.last().map_or_else(Fx::zero, |l| l.mul(&alpha0));
        powers.push(next);
    }
    let evaluate = |e: &[Q]| -> Fx { e.iter().zip(&powers).fold(Fx::zero(), |acc, (c, p)| acc.add(&p.mul(&to_fx(c)))) };
    let numeric: Vec<Fx> = poly.iter().map(|e| evaluate(e)).collect();
    let as_complex = |z: &Fx| -> Option<Complex64> {
        let scale = |v: &BigInt| num_traits::ToPrimitive::to_f64(&Q::new(v.clone(), crate::rules::poly::galois::one()));
        Some(Complex64::new(scale(&z.re)?, scale(&z.im)?))
    };
    let Some(coefficients): Option<Vec<Complex64>> = numeric.iter().map(as_complex).collect() else { return Vec::new() };
    let mut found: Vec<Elt> = Vec::new();
    let shift = BITS - 250;
    for start in super::resum::complex_roots(&coefficients) {
        let Some(root) = refine_root(&numeric, start) else { continue };
        // Rows: unit vectors with (Re, Im) of α^i, and of the root.
        let mut rows: Vec<Vec<BigInt>> = Vec::new();
        for (i, p) in powers.iter().chain(std::iter::once(&root)).enumerate() {
            let mut row = vec![BigInt::zero(); d + 1];
            if let Some(slot) = row.get_mut(i) {
                *slot = BigInt::one();
            }
            row.push(p.re.clone() >> shift);
            row.push(p.im.clone() >> shift);
            rows.push(row);
        }
        let Some(reduced) = crate::rules::poly::galois::lll(rows) else { continue };
        let tolerance = BigInt::one() << 60;
        let Some(relation) = reduced.into_iter().find(|r| {
            r.get(d).is_some_and(|c| !c.is_zero()) && r.get(d + 1).is_some_and(|v| v.abs() < tolerance) && r.get(d + 2).is_some_and(|v| v.abs() < tolerance)
        }) else {
            continue;
        };
        let Some(denominator) = relation.get(d).cloned() else { continue };
        let element: Elt = trim(
            (0..d)
                .map(|i| Q::new(-relation.get(i).cloned().unwrap_or_default(), denominator.clone()))
                .collect(),
        );
        // Exact check.
        let value = poly.iter().rev().fold(Vec::new(), |acc, c| field.add(&field.mul(&acc, &element), c));
        if value.is_empty() && !found.contains(&element) {
            found.push(element);
        }
    }
    found
}

/// A finite extension `K' = K(β)` of a number field `K = Q(α)` in which a polynomial
/// over `K` has a root.
pub(super) struct Extension {
    /// `K'` as a simple extension of `Q`.
    pub(super) field: Field,
    /// The root `β` in `K'`.
    pub(super) root: Elt,
    /// The image of `α` in `K'`.
    pub(super) alpha: Elt,
}

/// The value of a polynomial over `K` (an element is a polynomial in `α`) at `α ↦ h`
/// in `K'`.
fn map_element(
    target: &Field,
    element: &[Q],
    h: &[Q],
) -> Elt {
    element.iter().rev().fold(Vec::new(), |acc, c| target.add(&target.mul(&acc, h), &target.from_q(c)))
}

/// Roots of `poly` (over `K = Q(α)`, `K` not rational) that lie in extensions of `K`
/// but not in `K`: a tower step `K ⊂ K(β)`. The norm
/// `N(y) = Π_σ Π_β (y - (β + k σ(α)))` over the embeddings `σ` of `K` and the roots
/// `β` of `σ(poly)` is computed numerically and rationalised; for a suitable
/// small `k` it is square-free, and an irreducible factor `F` of it with `θ = β + k α`
/// as a root defines `K' = Q[y]/F` with `α` and `β` as polynomials in `θ` (found by
/// lattice reduction and checked exactly).
pub(super) fn tower_roots(
    field: &Field,
    poly: &[Elt],
) -> Vec<Extension> {
    let degree = poly.len().saturating_sub(1);
    let d = field.degree();
    if degree == 0 || field.is_rational() || d * degree > 24 {
        return Vec::new();
    }
    let to_fx = |q: &Q| Fx::real((q.numer() << BITS) / q.denom());
    let modulus: Vec<Fx> = field.modulus.iter().map(to_fx).collect();
    let Some(alphas): Option<Vec<Fx>> = field.embeddings().into_iter().map(|a| refine_root(&modulus, a)).collect() else {
        return Vec::new();
    };
    let value_at = |e: &[Q], alpha: &Fx| -> Fx { e.iter().rev().fold(Fx::zero(), |acc, c| acc.mul(alpha).add(&to_fx(c))) };
    let as_complex = |z: &Fx| -> Option<Complex64> {
        let scale = |v: &BigInt| ToPrimitive::to_f64(&Q::new(v.clone(), crate::rules::poly::galois::one()));
        Some(Complex64::new(scale(&z.re)?, scale(&z.im)?))
    };
    // All roots β of σ(poly) for every embedding σ: (σ index, β).
    let mut betas: Vec<(usize, Fx)> = Vec::new();
    for (index, alpha) in alphas.iter().enumerate() {
        let numeric: Vec<Fx> = poly.iter().map(|e| value_at(e, alpha)).collect();
        let Some(coefficients): Option<Vec<Complex64>> = numeric.iter().map(as_complex).collect() else { return Vec::new() };
        for start in super::resum::complex_roots(&coefficients) {
            let Some(root) = refine_root(&numeric, start) else { return Vec::new() };
            betas.push((index, root));
        }
    }
    if betas.len() != d * degree {
        return Vec::new();
    }
    let gap = BigInt::one() << (BITS - 80);
    for k in 0..8_i64 {
        let thetas: Vec<Fx> = betas
            .iter()
            .filter_map(|(i, b)| alphas.get(*i).map(|a| b.add(&a.mul(&Fx::integer(&BigInt::from(k))))))
            .collect();
        let distinct = thetas.iter().enumerate().all(|(i, a)| thetas.iter().skip(i + 1).all(|b| !a.close(b, 80)));
        let _ = &gap;
        if !distinct {
            continue;
        }
        // N(y) from its roots.
        let mut norm: Vec<Fx> = vec![Fx::integer(&BigInt::one())];
        for theta in &thetas {
            let mut next = vec![Fx::zero(); norm.len() + 1];
            for (i, a) in norm.iter().enumerate() {
                if let Some(slot) = next.get_mut(i + 1) {
                    *slot = slot.add(a);
                }
                if let Some(slot) = next.get_mut(i) {
                    *slot = slot.sub(&a.mul(theta));
                }
            }
            norm = next;
        }
        let Some(rational): Option<Vec<Q>> = norm.iter().map(crate::rules::poly::galois_high::rationalize).collect() else {
            continue;
        };
        let (_, factors) = univariate::factor(&rational);
        let Some(alpha0) = alphas.first() else { return Vec::new() };
        let mut out: Vec<Extension> = Vec::new();
        for (index, beta) in &betas {
            if *index != 0 {
                continue;
            }
            let theta0 = beta.add(&alpha0.mul(&Fx::integer(&BigInt::from(k))));
            for (f, multiplicity) in &factors {
                if *multiplicity != 1 || f.len() <= d + 1 {
                    continue;
                }
                let residual = crate::rules::poly::galois::horner(f, &theta0);
                if !residual.close(&Fx::zero(), 150) {
                    continue;
                }
                let lead = Q::from_integer(f.last().cloned().unwrap_or_else(BigInt::one));
                let monic: Vec<Q> = f.iter().map(|c| Q::from_integer(c.clone()) / &lead).collect();
                if out.iter().any(|e| e.field.modulus == monic) {
                    continue;
                }
                let target = Field::extension(monic);
                let n = target.degree();
                // α and β as polynomials in θ.
                let mut powers = vec![Fx::integer(&BigInt::one())];
                for _ in 1..n {
                    let next = powers.last().map_or_else(Fx::zero, |l| l.mul(&theta0));
                    powers.push(next);
                }
                let recognise = |value: &Fx| -> Option<Elt> {
                    let shift = BITS - 250;
                    let mut rows: Vec<Vec<BigInt>> = Vec::new();
                    for (i, p) in powers.iter().chain(std::iter::once(value)).enumerate() {
                        let mut row = vec![BigInt::zero(); n + 1];
                        if let Some(slot) = row.get_mut(i) {
                            *slot = BigInt::one();
                        }
                        row.push(p.re.clone() >> shift);
                        row.push(p.im.clone() >> shift);
                        rows.push(row);
                    }
                    let reduced = crate::rules::poly::galois::lll(rows)?;
                    let tolerance = BigInt::one() << 60;
                    let relation = reduced.into_iter().find(|r| {
                        r.get(n).is_some_and(|c| !c.is_zero())
                            && r.get(n + 1).is_some_and(|v| v.abs() < tolerance)
                            && r.get(n + 2).is_some_and(|v| v.abs() < tolerance)
                    })?;
                    let denominator = relation.get(n).cloned()?;
                    Some(trim((0..n).map(|i| Q::new(-relation.get(i).cloned().unwrap_or_default(), denominator.clone())).collect()))
                };
                let (Some(h), Some(beta_element)) = (recognise(alpha0), recognise(beta)) else { continue };
                // Exact checks: m(h) = 0 and poly(β) = 0 in K'.
                if !map_element(&target, &field.modulus, &h).is_empty() {
                    continue;
                }
                let mapped: Vec<Elt> = poly.iter().map(|e| map_element(&target, e, &h)).collect();
                let value = mapped.iter().rev().fold(Vec::new(), |acc, c| target.add(&target.mul(&acc, &beta_element), c));
                if !value.is_empty() {
                    continue;
                }
                out.push(Extension { field: target, root: beta_element, alpha: h });
            }
        }
        return out;
    }
    Vec::new()
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
    let mut towers: Vec<Extension> = Vec::new();
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
        } else {
            for rho in roots_in_field(field, &reduced) {
                roots.push((field.clone(), rho, 1));
            }
            for extension in tower_roots(field, &reduced) {
                towers.push(extension);
            }
        }
    }
    for extension in towers {
        let mapped = op.map_field(field, &extension.field, &extension.alpha);
        let mapped_exponent: Vec<(i64, Elt)> =
            exponent.iter().map(|(t, k)| (*t, map_element(&extension.field, k, &extension.alpha))).collect();
        if let Some((coefficients, log_terms)) = log_series(&extension.field, &mapped, j0, &extension.root, 0, search.terms) {
            search.out.push(FormalSolution {
                field: extension.field.clone(),
                ramification,
                exponent: mapped_exponent,
                rho: extension.root.clone(),
                coefficients,
                log_terms,
            });
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
    fn tower_when_the_second_level_does_not_split() {
        // x y'''' + (1 - 4x) y'' + 4x y = 0, i.e. x ((D^2 - 2)^2 + D^2 / x) y = 0: the
        // exponential part exp(κ x), κ = ±sqrt 2 (a double root of (κ^2 - 2)^2), is refined
        // by exp(ψ), 4 ψ^2 x + 1 = 0: its roots lie in Q(sqrt 2, i), a proper
        // extension of Q(sqrt 2).
        let operator = vec![qs(&[0, 4]), qs(&[]), qs(&[1, -4]), qs(&[]), qs(&[0, 1])];
        let solutions = formal_solutions(&operator, &Point::Infinity, 4);
        eprintln!("{solutions:?}");
        assert!(solutions.iter().any(|s| s.field.degree() == 4 && s.ramification == 2), "{solutions:?}");
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
    fn roots_in_a_number_field() {
        let field = Field::extension(qs(&[-2, 0, 1]));
        // k^2 - 2 over Q(sqrt 2): the roots +-alpha.
        let mut roots = roots_in_field(&field, &[qs(&[-2]), Vec::new(), qs(&[1])]);
        roots.sort();
        assert_eq!(roots, vec![qs(&[0, -1]), qs(&[0, 1])]);
        // k^2 - 3 has no root there.
        assert!(roots_in_field(&field, &[qs(&[-3]), Vec::new(), qs(&[1])]).is_empty());
        // (k - 1 - alpha)(k + alpha) = k^2 - k - 2 - alpha.
        let mut roots = roots_in_field(&field, &[qs(&[-2, -1]), qs(&[-1]), qs(&[1])]);
        roots.sort();
        assert_eq!(roots, vec![qs(&[0, -1]), qs(&[1, 1])]);
    }

    #[test]
    fn newton_polygon_slopes() {
        // y' = y / x^2 at 0: slope 1.
        let hull = lower_hull(&[(0, -1), (1, 0)]);
        assert_eq!(hull.len(), 1);
        assert_eq!(hull[0].slope(), q(1));
    }
}
