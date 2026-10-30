//! Exact arithmetic in a real number field `Q(theta)` and root isolation
//! over it.
//!
//! A [`Field`] is the field `Q[y]/(M(y))` for an irreducible `M` over `Q`,
//! together with a chosen real embedding: the root `theta` of `M` that lies
//! in the field's isolating interval. Because `M` is irreducible, an element
//! is zero exactly when its normal form is the zero polynomial, so zero
//! tests are structural; the sign of a non-zero element is obtained by
//! refining the isolating interval until interval Horner evaluation excludes
//! zero, which terminates because the evaluation map `Q(theta) -> R` is
//! injective.
//!
//! Fields are extended by adjoining a real root `c` of a square-free
//! polynomial `f` over the current field. The primitive element
//! `theta' = theta + s c` is found by trying `s = 1, 2, ...`: its minimal
//! polynomial over `Q` is the irreducible factor of
//! `Res_t(M(t), f_t(y - s t))` that contains `theta + s c`, and `theta` is
//! recovered by solving a linear system over `Q` in the basis
//! `theta^u c^v` of `Q[theta, c]/(M, f)`. If the system is inconsistent the
//! step is not primitive and the next `s` is tried. The recovered embedding
//! is verified symbolically (`M(theta) = 0`, `f(c) = 0`) and numerically
//! (the refined intervals of `theta` and `c`), so a wrong extension is
//! rejected rather than used.
//!
//! Polynomials with coefficients in a field support the same exact
//! operations as their rational counterparts: arithmetic, gcd, square-free
//! parts, Sturm chains, root counting and isolation. Lifting a cylindrical
//! algebraic decomposition uses them to find the fibre roots at an exact
//! algebraic sample point.

// Some field operations (refinement, polynomial evaluation, root
// refinement) are part of the exact-arithmetic contract and are used by the
// unit tests rather than by every lifting step.
#![allow(
    dead_code,
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::alg;
use super::alg::RealAlg;

/// Bound on the degree of a computed primitive element.
const MAX_FIELD_DEGREE: usize = 128;
/// Bound on the number of interpolation points.
const MAX_INTERPOLATION: usize = 256;
/// Bound on the bisection steps of a field-refinement.
const MAX_BISECT: usize = 4096;
/// Smallest positive integer shifts tried for a primitive element.
const MAX_SHIFTS: u32 = 12;

/// An element of `Q(theta)`, in the power basis of degree below
/// `deg(M)`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Elem {
    coeffs: Vec<BigRational>,
}

impl Elem {
    /// Zero.
    #[must_use]
    pub(crate) fn zero() -> Self {
        Self::default()
    }

    /// One.
    #[must_use]
    pub(crate) fn one() -> Self {
        Self {
            coeffs: vec![BigRational::one()],
        }
    }

    /// A rational constant.
    #[must_use]
    pub(crate) fn rational(value: BigRational) -> Self {
        Self { coeffs: vec![value] }
    }

    /// The coefficients in the power basis.
    #[must_use]
    pub(crate) fn coeffs(&self) -> &[BigRational] {
        &self.coeffs
    }

    /// Builds an element from coefficients, trimming trailing zeros.
    #[must_use]
    pub(crate) fn from_coeffs(mut coeffs: Vec<BigRational>) -> Self {
        alg::trim(&mut coeffs);
        Self { coeffs }
    }

    /// Whether the element is zero.
    #[must_use]
    pub(crate) fn is_zero(&self) -> bool {
        self.coeffs.iter().all(Zero::is_zero)
    }

    /// The value when the element is rational.
    #[must_use]
    pub(crate) fn constant(&self) -> Option<BigRational> {
        if self.is_zero() {
            return Some(BigRational::zero());
        }
        let mut nonzero = self.coeffs.iter().enumerate().filter(|(_, c)| !c.is_zero());
        let (first, value) = nonzero.next()?;
        if first != 0 || nonzero.next().is_some() {
            return None;
        }
        Some(value.clone())
    }

    /// Sum.
    #[must_use]
    pub(crate) fn add(
        &self,
        other: &Self,
    ) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut out = vec![BigRational::zero(); n];
        for (i, slot) in out.iter_mut().enumerate() {
            let a = self
                .coeffs
                .get(i)
                .cloned()
                .unwrap_or_else(BigRational::zero);
            let b = other
                .coeffs
                .get(i)
                .cloned()
                .unwrap_or_else(BigRational::zero);
            *slot = a + b;
        }
        Self::from_coeffs(out)
    }

    /// Difference.
    #[must_use]
    pub(crate) fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }

    /// Additive inverse.
    #[must_use]
    pub(crate) fn neg(&self) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(|c| -c.clone()).collect(),
        }
    }

    /// Product with a rational.
    #[must_use]
    pub(crate) fn scale(
        &self,
        factor: &BigRational,
    ) -> Self {
        Self::from_coeffs(self.coeffs.iter().map(|c| c * factor).collect())
    }
}

/// `Q(theta)` with the chosen real embedding of `theta`.
#[derive(Clone, Debug)]
pub(crate) struct Field {
    root: RealAlg,
}

impl Field {
    /// The rational field `Q`, with generator `0`.
    #[must_use]
    pub(crate) fn rational_field() -> Option<Self> {
        let poly = vec![BigRational::zero(), BigRational::one()];
        let zero = BigRational::zero();
        let root = RealAlg::from_root(&poly, zero.clone(), zero)?;
        Some(Self { root })
    }

    /// The field generated by the root of the irreducible `poly` in the
    /// isolating interval.
    #[must_use]
    pub(crate) fn from_root(
        poly: &[BigRational],
        lo: BigRational,
        hi: BigRational,
    ) -> Option<Self> {
        let root = RealAlg::from_root(poly, lo, hi)?;
        let degree = alg::degree(root.polynomial())?;
        if degree == 0 || degree > MAX_FIELD_DEGREE {
            return None;
        }
        Some(Self { root })
    }

    /// The defining root.
    #[must_use]
    pub(crate) fn root(&self) -> &RealAlg {
        &self.root
    }

    /// The minimal polynomial of `theta`, monic.
    #[must_use]
    pub(crate) fn minpoly(&self) -> &[BigRational] {
        self.root.polynomial()
    }

    /// The degree of the field over `Q`.
    #[must_use]
    pub(crate) fn degree(&self) -> usize {
        alg::degree(self.root.polynomial()).unwrap_or(0)
    }

    /// The embedding value when the generator is rational.
    #[must_use]
    pub(crate) fn generator_value(&self) -> Option<BigRational> {
        self.root.rational_value()
    }

    /// Shrinks the isolating interval of the generator.
    pub(crate) fn refine(
        &mut self,
        width: &BigRational,
    ) {
        self.root.refine_to(width);
    }

    /// The sign of `value` at the real embedding, exactly.
    #[must_use]
    pub(crate) fn sign(
        &mut self,
        value: &Elem,
    ) -> Option<i32> {
        if value.is_zero() {
            return Some(0);
        }
        self.root.sign_of_mut(value.coeffs())
    }

    /// Product, reduced modulo the minimal polynomial.
    #[must_use]
    pub(crate) fn mul(
        &self,
        a: &Elem,
        b: &Elem,
    ) -> Elem {
        let product = dense_mul(a.coeffs(), b.coeffs());
        Elem::from_coeffs(reduce_theta(&product, self.minpoly()))
    }

    /// Power.
    #[must_use]
    pub(crate) fn pow(
        &self,
        a: &Elem,
        mut exponent: u32,
    ) -> Elem {
        let mut result = Elem::one();
        let mut base = a.clone();
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = self.mul(&result, &base);
            }
            exponent >>= 1;
            if exponent > 0 {
                base = self.mul(&base, &base);
            }
        }
        result
    }

    /// Multiplicative inverse, by the extended Euclidean algorithm.
    #[must_use]
    pub(crate) fn inv(
        &self,
        a: &Elem,
    ) -> Option<Elem> {
        if a.is_zero() {
            return None;
        }
        let (g, s) = extended_gcd(a.coeffs(), self.minpoly())?;
        if alg::degree(&g) != Some(0) {
            return None;
        }
        let lead = g.first()?.clone();
        Some(Elem::from_coeffs(alg::scale(&s, &lead.recip())))
    }

    /// Division.
    #[must_use]
    pub(crate) fn div(
        &self,
        a: &Elem,
        b: &Elem,
    ) -> Option<Elem> {
        Some(self.mul(a, &self.inv(b)?))
    }

    /// The generator `theta` itself.
    #[must_use]
    pub(crate) fn theta(&self) -> Elem {
        match self.root.rational_value() {
            | Some(value) => Elem::rational(value),
            | None => Elem::from_coeffs(vec![BigRational::zero(), BigRational::one()]),
        }
    }

    /// Composes a polynomial expression with `theta_expr`, giving an
    /// element of the field `target`.
    #[must_use]
    pub(crate) fn substitute(
        value: &Elem,
        theta_expr: &Elem,
        target: &Self,
    ) -> Elem {
        let mut acc = Elem::zero();
        for c in value.coeffs().iter().rev() {
            acc = target.mul(&acc, theta_expr);
            acc = acc.add(&Elem::rational(c.clone()));
        }
        acc
    }

    /// The interval enclosure of `value` at the current generator interval.
    #[must_use]
    pub(crate) fn enclosure(
        &self,
        value: &Elem,
    ) -> (BigRational, BigRational) {
        let (lo, hi) = self.root.interval();
        alg::eval_interval(value.coeffs(), lo, hi)
    }
}

/// Convolution of dense rational coefficients.
fn dense_mul(
    a: &[BigRational],
    b: &[BigRational],
) -> Vec<BigRational> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let n = a.len().saturating_add(b.len()).saturating_sub(1);
    let mut out = vec![BigRational::zero(); n];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i.saturating_add(j)) {
                *slot += x * y;
            }
        }
    }
    alg::trim(&mut out);
    out
}

/// Reduces `p` modulo the monic `minpoly`.
#[must_use]
pub(crate) fn reduce_theta(
    p: &[BigRational],
    minpoly: &[BigRational],
) -> Vec<BigRational> {
    let Some(d) = alg::degree(minpoly) else {
        return Vec::new();
    };
    let mut out = p.to_vec();
    alg::trim(&mut out);
    let Some(lead) = minpoly.get(d) else {
        return out;
    };
    if lead.is_zero() {
        return out;
    }
    while let Some(n) = alg::degree(&out) {
        if n < d {
            break;
        }
        let Some(top) = out.get(n).cloned() else {
            break;
        };
        let factor = top / lead;
        let shift = n.saturating_sub(d);
        // out -= factor * x^shift * minpoly
        for (k, c) in minpoly.iter().take(d.saturating_add(1)).enumerate() {
            if let Some(slot) = out.get_mut(shift.saturating_add(k)) {
                *slot -= &factor * c;
            }
        }
        alg::trim(&mut out);
    }
    out
}

/// `x^exponent` reduced modulo the monic `minpoly`.
fn power_mod(
    mut exponent: u32,
    minpoly: &[BigRational],
) -> Vec<BigRational> {
    let mut result = vec![BigRational::one()];
    let mut base = vec![BigRational::zero(), BigRational::one()];
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = reduce_theta(&dense_mul(&result, &base), minpoly);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = reduce_theta(&dense_mul(&base, &base), minpoly);
        }
    }
    result
}

/// Extended Euclidean algorithm over `Q`: `(g, s)` with `s*a + t*b = g`.
fn extended_gcd(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(Vec<BigRational>, Vec<BigRational>)> {
    let mut old_r = a.to_vec();
    let mut r = b.to_vec();
    alg::trim(&mut old_r);
    alg::trim(&mut r);
    let mut old_s = vec![BigRational::one()];
    let mut s = Vec::new();
    let mut steps = 0_usize;
    while alg::degree(&r).is_some() {
        let (q, rem) = alg::divrem(&old_r, &r)?;
        let next_s = alg::sub(&old_s, &alg::mul(&q, &s)?);
        old_r = r;
        r = rem;
        old_s = s;
        s = next_s;
        steps = steps.saturating_add(1);
        if steps > 100_000 {
            return None;
        }
    }
    Some((old_r, old_s))
}

/// A polynomial with coefficients in a field.
pub(crate) type FPoly = Vec<Elem>;

/// Sum.
#[must_use]
pub(crate) fn f_add(
    a: &[Elem],
    b: &[Elem],
) -> FPoly {
    let n = a.len().max(b.len());
    let mut out = vec![Elem::zero(); n];
    for (i, slot) in out.iter_mut().enumerate() {
        let x = a.get(i).cloned().unwrap_or_else(Elem::zero);
        let y = b.get(i).cloned().unwrap_or_else(Elem::zero);
        *slot = x.add(&y);
    }
    f_trim(&mut out);
    out
}

/// Difference.
#[must_use]
pub(crate) fn f_sub(
    a: &[Elem],
    b: &[Elem],
) -> FPoly {
    let n = a.len().max(b.len());
    let mut out = vec![Elem::zero(); n];
    for (i, slot) in out.iter_mut().enumerate() {
        let x = a.get(i).cloned().unwrap_or_else(Elem::zero);
        let y = b.get(i).cloned().unwrap_or_else(Elem::zero);
        *slot = x.sub(&y);
    }
    f_trim(&mut out);
    out
}

/// Removes trailing zeros.
pub(crate) fn f_trim(p: &mut FPoly) {
    while p.last().is_some_and(Elem::is_zero) {
        p.pop();
    }
}

/// Degree in the field sense.
#[must_use]
pub(crate) fn f_degree(p: &[Elem]) -> Option<usize> {
    let mut n = p.len();
    while n > 0 && p.get(n.saturating_sub(1)).is_some_and(Elem::is_zero) {
        n = n.saturating_sub(1);
    }
    (n > 0).then(|| n.saturating_sub(1))
}

/// Product over a field.
#[must_use]
pub(crate) fn f_mul(
    a: &[Elem],
    b: &[Elem],
    field: &Field,
) -> FPoly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let n = a.len().saturating_add(b.len()).saturating_sub(1);
    let mut out = vec![Elem::zero(); n];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i.saturating_add(j)) {
                *slot = slot.add(&field.mul(x, y));
            }
        }
    }
    f_trim(&mut out);
    out
}

/// Quotient and remainder over a field.
#[must_use]
pub(crate) fn f_divrem(
    a: &[Elem],
    b: &[Elem],
    field: &Field,
) -> Option<(FPoly, FPoly)> {
    let bn = f_degree(b)?;
    let lead = b.get(bn)?.clone();
    let inverse = field.inv(&lead)?;
    let mut rest = a.to_vec();
    f_trim(&mut rest);
    let Some(an) = f_degree(&rest) else {
        return Some((Vec::new(), Vec::new()));
    };
    if an < bn {
        return Some((Vec::new(), rest));
    }
    let mut quotient = vec![Elem::zero(); an.saturating_sub(bn).saturating_add(1)];
    let mut steps = 0_usize;
    while let Some(degree) = f_degree(&rest) {
        if degree < bn {
            break;
        }
        let shift = degree.saturating_sub(bn);
        let top = rest.get(degree)?.clone();
        let factor = field.mul(&top, &inverse);
        if let Some(slot) = quotient.get_mut(shift) {
            *slot = factor.clone();
        }
        let scaled = f_scale(b, &factor, field);
        let mut shifted = vec![Elem::zero(); shift];
        shifted.extend(scaled);
        rest = f_sub(&rest, &shifted);
        f_trim(&mut rest);
        steps = steps.saturating_add(1);
        if steps > 100_000 {
            return None;
        }
    }
    Some((quotient, rest))
}

/// Product with a field element.
#[must_use]
pub(crate) fn f_scale(
    a: &[Elem],
    factor: &Elem,
    field: &Field,
) -> FPoly {
    let mut out: FPoly = a.iter().map(|c| field.mul(c, factor)).collect();
    f_trim(&mut out);
    out
}

/// Monic form.
#[must_use]
pub(crate) fn f_monic(
    a: &[Elem],
    field: &Field,
) -> FPoly {
    match f_degree(a) {
        | None => Vec::new(),
        | Some(n) => match a.get(n).and_then(|l| field.inv(l)) {
            | Some(inverse) => f_scale(a, &inverse, field),
            | None => a.to_vec(),
        },
    }
}

/// Derivative over a field.
#[must_use]
pub(crate) fn f_derivative(
    a: &[Elem],
    field: &Field,
) -> FPoly {
    let mut out = Vec::new();
    for (i, c) in a.iter().enumerate().skip(1) {
        let factor = Elem::rational(BigRational::from_integer(BigInt::from(i)));
        out.push(field.mul(c, &factor));
    }
    f_trim(&mut out);
    out
}

/// Monic gcd over a field.
#[must_use]
pub(crate) fn f_gcd(
    a: &[Elem],
    b: &[Elem],
    field: &Field,
) -> FPoly {
    let mut x = f_monic(a, field);
    let mut y = f_monic(b, field);
    let mut steps = 0_usize;
    while f_degree(&y).is_some() {
        let Some((_, r)) = f_divrem(&x, &y, field) else {
            break;
        };
        x = y;
        y = f_monic(&r, field);
        steps = steps.saturating_add(1);
        if steps > 100_000 {
            break;
        }
    }
    x
}

/// Square-free part over a field.
#[must_use]
pub(crate) fn f_square_free(
    a: &[Elem],
    field: &Field,
) -> FPoly {
    let d = f_derivative(a, field);
    let g = f_gcd(a, &d, field);
    if f_degree(&g).is_some() {
        f_divrem(a, &g, field).map_or_else(|| f_monic(a, field), |(q, _)| f_monic(&q, field))
    } else {
        f_monic(a, field)
    }
}

/// Evaluates at a rational point.
#[must_use]
pub(crate) fn f_eval_rational(
    a: &[Elem],
    x: &BigRational,
    field: &Field,
) -> Elem {
    let point = Elem::rational(x.clone());
    let mut acc = Elem::zero();
    for c in a.iter().rev() {
        acc = field.mul(&acc, &point);
        acc = acc.add(c);
    }
    acc
}

/// Evaluates at a field element.
#[must_use]
pub(crate) fn f_eval(
    a: &[Elem],
    x: &Elem,
    field: &Field,
) -> Elem {
    let mut acc = Elem::zero();
    for c in a.iter().rev() {
        acc = field.mul(&acc, x);
        acc = acc.add(c);
    }
    acc
}

fn f_sign_at(
    a: &[Elem],
    at: Option<&BigRational>,
    positive: bool,
    field: &mut Field,
) -> Option<i32> {
    match at {
        | Some(x) => {
            let value = f_eval_rational(a, x, field);
            field.sign(&value)
        },
        | None => {
            let lead = a.last().cloned().unwrap_or_else(Elem::zero);
            let sign = field.sign(&lead)?;
            let odd = a.len().is_multiple_of(2);
            Some(if positive || !odd {
                sign
            } else {
                -sign
            })
        },
    }
}

fn f_variations(
    chain: &[FPoly],
    at: Option<&BigRational>,
    positive: bool,
    field: &mut Field,
) -> Option<usize> {
    let mut signs = Vec::with_capacity(chain.len());
    for q in chain {
        let sign = f_sign_at(q, at, positive, field)?;
        if sign != 0 {
            signs.push(sign);
        }
    }
    Some(signs.windows(2).filter(|w| w.first() != w.get(1)).count())
}

/// The Sturm chain of a square-free field polynomial.
#[must_use]
pub(crate) fn f_sturm_chain(
    a: &[Elem],
    field: &Field,
) -> Option<Vec<FPoly>> {
    let mut chain = vec![a.to_vec()];
    let d = f_derivative(a, field);
    if f_degree(&d).is_some() {
        chain.push(d);
    }
    let mut steps = 0_usize;
    loop {
        let n = chain.len();
        let Some(b) = chain.get(n.saturating_sub(1)) else {
            break;
        };
        if n < 2 {
            break;
        }
        let Some(a_prev) = chain.get(n.saturating_sub(2)) else {
            break;
        };
        let Some((_, rem)) = f_divrem(a_prev, b, field) else {
            break;
        };
        if f_degree(&rem).is_none() {
            break;
        }
        chain.push(f_scale(
            &rem,
            &Elem::rational(BigRational::from_integer(BigInt::from(-1))),
            field,
        ));
        steps = steps.saturating_add(1);
        if steps > 100_000 {
            return None;
        }
    }
    chain.retain(|q| f_degree(q).is_some());
    Some(chain)
}

/// Roots in `(a, b]`.
#[must_use]
pub(crate) fn f_sturm_count(
    chain: &[FPoly],
    a: Option<&BigRational>,
    b: Option<&BigRational>,
    field: &mut Field,
) -> Option<usize> {
    let low = f_variations(chain, a, false, field)?;
    let high = f_variations(chain, b, true, field)?;
    Some(low.saturating_sub(high))
}

/// Roots in `[a, b]`, with a precomputed Sturm chain.
#[must_use]
fn f_count_closed_chain(
    chain: &[FPoly],
    p: &[Elem],
    a: &BigRational,
    b: &BigRational,
    field: &mut Field,
) -> Option<usize> {
    if a > b {
        return Some(0);
    }
    if a == b {
        let value = f_eval_rational(p, a, field);
        return Some(usize::from(value.is_zero()));
    }
    let mut count = f_sturm_count(chain, Some(a), Some(b), field)?;
    let at_a = f_eval_rational(p, a, field);
    if at_a.is_zero() {
        count = count.saturating_add(1);
    }
    Some(count)
}

/// Roots in `(a, b)`, with a precomputed Sturm chain.
#[must_use]
fn f_count_open_chain(
    chain: &[FPoly],
    p: &[Elem],
    a: &BigRational,
    b: &BigRational,
    field: &mut Field,
) -> Option<usize> {
    if a >= b {
        return Some(0);
    }
    let upper = f_sturm_count(chain, Some(a), Some(b), field)?;
    let at_b = f_eval_rational(p, b, field);
    Some(if at_b.is_zero() {
        upper.saturating_sub(1)
    } else {
        upper
    })
}

/// Roots in `[a, b]`.
#[must_use]
pub(crate) fn f_count_closed(
    p: &[Elem],
    a: &BigRational,
    b: &BigRational,
    field: &mut Field,
) -> Option<usize> {
    let chain = f_sturm_chain(p, field)?;
    f_count_closed_chain(&chain, p, a, b, field)
}

/// Roots in `(a, b)`.
#[must_use]
pub(crate) fn f_count_open(
    p: &[Elem],
    a: &BigRational,
    b: &BigRational,
    field: &mut Field,
) -> Option<usize> {
    let chain = f_sturm_chain(p, field)?;
    f_count_open_chain(&chain, p, a, b, field)
}

/// A rational upper bound for `|coefficient|` and a positive lower bound
/// for `|leading coefficient|` at the embedding.
fn f_cauchy_bound(
    p: &[Elem],
    field: &mut Field,
) -> Option<BigRational> {
    let n = f_degree(p)?;
    let lead = p.get(n)?.clone();
    let mut lower = BigRational::zero();
    for _ in 0..MAX_BISECT {
        let (lo, hi) = field.enclosure(&lead);
        let candidate = lo.abs().min(hi.abs());
        if lo.is_positive() || hi.is_negative() {
            lower = candidate;
            break;
        }
        if field.generator_value().is_some() {
            break;
        }
        field.root_mut_refine_once();
    }
    if !lower.is_positive() {
        return None;
    }
    let mut max = BigRational::zero();
    for c in p.iter().take(n.saturating_add(1)) {
        let (lo, hi) = field.enclosure(c);
        let magnitude = lo.abs().max(hi.abs());
        if magnitude > max {
            max = magnitude;
        }
    }
    Some(max / lower + BigRational::one())
}

/// Isolating intervals of the real roots of the square-free field
/// polynomial `p`, ascending.
#[must_use]
pub(crate) fn f_isolate(
    p: &[Elem],
    field: &mut Field,
) -> Option<Vec<(BigRational, BigRational)>> {
    let p = f_square_free(p, field);
    if f_degree(&p).is_none() {
        return Some(Vec::new());
    }
    let bound = f_cauchy_bound(&p, field)?;
    let chain = f_sturm_chain(&p, field)?;
    let mut out = Vec::new();
    let mut stack = vec![(-bound.clone(), bound)];
    let mut work = 0_usize;
    while let Some((a, b)) = stack.pop() {
        work = work.saturating_add(1);
        if work > 20_000 {
            return None;
        }
        let count = f_count_open_chain(&chain, &p, &a, &b, field)?;
        match count {
            | 0 => {},
            | 1 => out.push((a, b)),
            | _ => {
                let two = BigRational::from_integer(BigInt::from(2));
                let mid = (&a + &b) / &two;
                let at_mid = f_eval_rational(&p, &mid, field);
                if at_mid.is_zero() {
                    out.push((mid.clone(), mid.clone()));
                    let left = f_count_open_chain(&chain, &p, &a, &mid, field)?;
                    let right = f_count_open_chain(&chain, &p, &mid, &b, field)?;
                    if left > 0 {
                        // The midpoint is a root of a neighbouring cell; move
                        // the new endpoint just below it so that the interval
                        // does not end at a root.
                        let mut upper = mid.clone();
                        for k in 1..=128_u32 {
                            let step = (&mid - &a)
                                / num_traits::pow(
                                    BigRational::from_integer(BigInt::from(2)),
                                    usize::try_from(k).ok()?,
                                );
                            let candidate = &mid - &step;
                            if candidate <= a {
                                break;
                            }
                            if f_count_open_chain(&chain, &p, &candidate, &mid, field) == Some(0)
                                && f_count_open_chain(&chain, &p, &a, &candidate, field)
                                    == Some(left)
                            {
                                upper = candidate;
                                break;
                            }
                        }
                        stack.push((a, upper));
                    }
                    if right > 0 {
                        let mut lower = mid.clone();
                        for k in 1..=128_u32 {
                            let step = (&b - &mid)
                                / num_traits::pow(
                                    BigRational::from_integer(BigInt::from(2)),
                                    usize::try_from(k).ok()?,
                                );
                            let candidate = &mid + &step;
                            if candidate >= b {
                                break;
                            }
                            if f_count_open_chain(&chain, &p, &mid, &candidate, field) == Some(0)
                                && f_count_open_chain(&chain, &p, &candidate, &b, field)
                                    == Some(right)
                            {
                                lower = candidate;
                                break;
                            }
                        }
                        stack.push((lower, b));
                    }
                } else {
                    let left = f_count_open_chain(&chain, &p, &a, &mid, field)?;
                    let right = f_count_open_chain(&chain, &p, &mid, &b, field)?;
                    if left > 0 {
                        stack.push((a, mid.clone()));
                    }
                    if right > 0 {
                        stack.push((mid, b));
                    }
                }
            },
        }
    }
    for interval in &mut out {
        let width = BigRational::new(BigInt::from(1), BigInt::from(8));
        *interval = f_refine_root_chain(&chain, &p, interval, &width, field)?;
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    Some(out)
}

/// Shrinks a root interval of the square-free `p` below `width`, using a
/// precomputed Sturm chain.
#[must_use]
fn f_refine_root_chain(
    chain: &[FPoly],
    p: &[Elem],
    interval: &(BigRational, BigRational),
    width: &BigRational,
    field: &mut Field,
) -> Option<(BigRational, BigRational)> {
    let mut a = interval.0.clone();
    let mut b = interval.1.clone();
    for _ in 0..MAX_BISECT {
        if a == b || &b - &a <= *width {
            return Some((a, b));
        }
        let two = BigRational::from_integer(BigInt::from(2));
        let mid = (&a + &b) / &two;
        if f_eval_rational(p, &mid, field).is_zero() {
            return Some((mid.clone(), mid));
        }
        // Count open halves: the interval's endpoints are never roots, so
        // both counts together account for the isolated root.
        let left = f_count_open_chain(chain, p, &a, &mid, field)?;
        if left >= 1 {
            b = mid;
            continue;
        }
        let right = f_count_open_chain(chain, p, &mid, &b, field)?;
        if right >= 1 {
            a = mid;
        } else {
            return Some((a, b));
        }
    }
    Some((a, b))
}

/// Shrinks a root interval of the square-free `p` below `width`.
#[must_use]
pub(crate) fn f_refine_root(
    p: &[Elem],
    interval: &(BigRational, BigRational),
    width: &BigRational,
    field: &mut Field,
) -> Option<(BigRational, BigRational)> {
    let chain = f_sturm_chain(p, field)?;
    f_refine_root_chain(&chain, p, interval, width, field)
}

/// The result of adjoining a real root `c` of `f` to a field.
pub(crate) struct Extension {
    /// The larger field `Q(theta, c) = Q(theta')`.
    pub(crate) field: Field,
    /// The old generator `theta` as an element of the new field.
    pub(crate) theta: Elem,
    /// The new root `c` as an element of the new field.
    pub(crate) root: Elem,
}

/// Adjoins the real root `c` of the square-free `f` that lies in
/// `interval` (an isolating interval over the current field).
#[must_use]
pub(crate) fn extend(
    field: &mut Field,
    f: &[Elem],
    interval: &(BigRational, BigRational),
) -> Option<Extension> {
    let f = f_square_free(f, field);
    let m = f_degree(&f)?;
    if m == 0 || m > MAX_FIELD_DEGREE {
        return None;
    }
    let d = field.degree();
    if d.checked_mul(m)? > MAX_FIELD_DEGREE {
        return None;
    }
    let degree_bound = d.checked_mul(m)?;
    if degree_bound > MAX_INTERPOLATION {
        return None;
    }
    for s in 1..=MAX_SHIFTS {
        if let Some(extension) = try_shift(field, &f, interval, s) {
            return Some(extension);
        }
    }
    None
}

/// One primitive-element attempt with `theta' = theta + s c`.
#[allow(clippy::single_call_fn)] // one call site, kept separate for clarity
fn try_shift(
    field: &mut Field,
    f: &[Elem],
    interval: &(BigRational, BigRational),
    s: u32,
) -> Option<Extension> {
    let m = f_degree(f)?;
    let d = field.degree();
    let degree_bound = d.checked_mul(m)?;
    let shift = BigRational::from_integer(BigInt::from(s));
    // B_k(t) for k = 0..=m: G(t, y) = sum_k y^k B_k(t).
    let mut b_polys: Vec<Vec<BigRational>> = Vec::with_capacity(m.saturating_add(1));
    for k in 0..=m {
        let mut limb: Vec<BigRational> = Vec::new();
        for j in k..=m {
            let coefficient = f.get(j)?;
            let factor = BigRational::from_integer(binomial(j, k));
            let power = u32::try_from(j.saturating_sub(k)).ok()?;
            let shift_power = num_traits::pow(shift.clone(), usize::try_from(power).ok()?);
            let scale = &factor
                * &shift_power
                * BigRational::from_integer(BigInt::from(if (j - k).is_multiple_of(2) {
                    1
                } else {
                    -1
                }));
            let scaled = alg::scale(coefficient.coeffs(), &scale);
            let shifted: Vec<BigRational> =
                std::iter::repeat_n(BigRational::zero(), power as usize)
                    .chain(scaled)
                    .collect();
            limb = alg::add(&limb, &shifted);
        }
        b_polys.push(limb);
    }
    // Interpolate H(y) = Res_t(M(t), G(t, y)) at degree_bound + 1 points.
    let count = degree_bound.saturating_add(1);
    let mut xs: Vec<BigRational> = Vec::with_capacity(count);
    let mut values: Vec<BigRational> = Vec::with_capacity(count);
    for i in 0..count {
        let x = interpolation_point(i);
        let mut at_x = Vec::new();
        let mut power = BigRational::one();
        for poly in &b_polys {
            at_x = alg::add(&at_x, &alg::scale(poly, &power));
            power *= &x;
        }
        let value = alg::resultant(field.minpoly(), &at_x)?;
        xs.push(x);
        values.push(value);
    }
    let h = interpolate(&xs, &values)?;
    alg::degree(&h)?;
    // Refine the interval of theta + s c until it isolates one root of H.
    let (theta_lo, theta_hi) = field.root().interval();
    let (mut theta_lo, mut theta_hi) = (theta_lo.clone(), theta_hi.clone());
    let h_chain = alg::sturm_chain(&h);
    let f_chain = f_sturm_chain(f, field)?;
    let mut c_lo = interval.0.clone();
    let mut c_hi = interval.1.clone();
    let mut chosen = None;
    for _ in 0..512 {
        let lo = &theta_lo + &shift * &c_lo;
        let hi = &theta_hi + &shift * &c_hi;
        let count = match lo.cmp(&hi) {
            | std::cmp::Ordering::Greater => 0,
            | std::cmp::Ordering::Equal => usize::from(alg::eval(&h, &lo).is_zero()),
            | std::cmp::Ordering::Less => alg::sturm_count(&h_chain, Some(&lo), Some(&hi))
                .saturating_add(usize::from(alg::eval(&h, &lo).is_zero())),
        };
        if count == 1 {
            chosen = Some((lo, hi));
            break;
        }
        // Refine both factors of the primitive element: c directly, and
        // theta through the current field's embedding.
        if c_lo == c_hi {
            let before = field.root().interval();
            let before = (before.0.clone(), before.1.clone());
            field.root_mut_refine_once();
            let after = field.root().interval();
            let after = (after.0.clone(), after.1.clone());
            if after == before {
                break;
            }
            theta_lo = after.0;
            theta_hi = after.1;
        } else {
            let two = BigRational::from_integer(BigInt::from(2));
            let mid = (&c_lo + &c_hi) / &two;
            let at_mid = f_eval_rational(f, &mid, field);
            if at_mid.is_zero() {
                c_lo = mid.clone();
                c_hi = mid;
            } else {
                let left = f_count_closed_chain(&f_chain, f, &c_lo, &mid, field)?;
                if left >= 1 {
                    c_hi = mid;
                } else {
                    c_lo = mid;
                }
            }
        }
    }
    let (lo, hi) = chosen?;
    // Pick the irreducible factor of H that contains theta'.
    let (_, factors) = crate::rules::poly::univariate::factor(&h);
    let mut minimal: Option<Vec<BigRational>> = None;
    for (factor, _) in &factors {
        let rational: Vec<BigRational> = factor
            .iter()
            .cloned()
            .map(BigRational::from_integer)
            .collect();
        let factor = alg::monic(&rational);
        if alg::count_closed(&factor, &lo, &hi) == 1 {
            if minimal.is_some() {
                return None;
            }
            minimal = Some(factor);
        }
    }
    let minimal = minimal?;
    let new_degree = alg::degree(&minimal)?;
    if new_degree > MAX_FIELD_DEGREE {
        return None;
    }
    let new_field = Field::from_root(&minimal, lo, hi)?;
    let rows = d.checked_mul(m)?;
    // A column j is (theta + s c)^j expanded on the basis theta^u c^v.
    let mut columns: Vec<Vec<BigRational>> = Vec::with_capacity(new_degree);
    for j in 0..new_degree {
        let mut out = vec![BigRational::zero(); rows];
        for v in 0..=j {
            let coefficient =
                BigRational::from_integer(binomial(j, v)) * num_traits::pow(shift.clone(), v);
            let u = j.saturating_sub(v);
            add_ring_monomial(&mut out, field.minpoly(), f, m, u, v, &coefficient);
        }
        columns.push(out);
    }
    // A[row][j] = coefficient of the basis element `row` in (theta')^j.
    let mut matrix: Vec<Vec<BigRational>> = vec![vec![BigRational::zero(); new_degree]; rows];
    for (column_index, column) in columns.iter().enumerate() {
        for (row, value) in column.iter().enumerate() {
            if let Some(slot) = matrix.get_mut(row).and_then(|r| r.get_mut(column_index)) {
                *slot = value.clone();
            }
        }
    }
    // Right-hand side: theta itself, in the basis theta^u c^v.
    let mut rhs = vec![BigRational::zero(); rows];
    for (u, coefficient) in field.theta().coeffs().iter().enumerate() {
        if let Some(slot) = rhs.get_mut(u.saturating_mul(m)) {
            *slot = coefficient.clone();
        }
    }
    let solution = solve(&matrix, &rhs, new_degree)?;
    let theta_expr = Elem::from_coeffs(solution);
    let generator = new_field.theta();
    let root_expr = new_field.div(&generator.sub(&theta_expr), &Elem::rational(shift))?;
    // Verify the relations.
    let old_minpoly = field.minpoly().to_vec();
    let evaluated = Field::substitute(&Elem::from_coeffs(old_minpoly), &theta_expr, &new_field);
    if !evaluated.is_zero() {
        return None;
    }
    let mut value = Elem::zero();
    for (j, coefficient) in f.iter().enumerate() {
        let substituted = Field::substitute(coefficient, &theta_expr, &new_field);
        let power = new_field.pow(&root_expr, u32::try_from(j).ok()?);
        value = value.add(&new_field.mul(&substituted, &power));
    }
    if !value.is_zero() {
        return None;
    }
    if !embedding_matches(&new_field, &theta_expr, field.root()) {
        return None;
    }
    if !enclosure_inside(&new_field, &root_expr, interval) {
        return None;
    }
    Some(Extension {
        field: new_field,
        theta: theta_expr,
        root: root_expr,
    })
}

/// Adds `coefficient * theta^u c^v` to a normal form of
/// `Q[theta, c]/(M, f)`, in the basis `theta^a c^b` with `a < deg M`,
/// `b < deg f`.
fn add_ring_monomial(
    out: &mut [BigRational],
    minpoly: &[BigRational],
    f: &[Elem],
    m: usize,
    u: usize,
    v: usize,
    coefficient: &BigRational,
) {
    let d = alg::degree(minpoly).unwrap_or(0);
    if d == 0 || coefficient.is_zero() {
        return;
    }
    if v >= m {
        // f is monic: c^m = -sum_{k<m} f_k c^k.
        for (k, coefficient_k) in f.iter().take(m).enumerate() {
            for (t, c) in coefficient_k.coeffs().iter().enumerate() {
                let next = -(coefficient * c);
                add_ring_monomial(
                    out,
                    minpoly,
                    f,
                    m,
                    u.saturating_add(t),
                    v.saturating_sub(m).saturating_add(k),
                    &next,
                );
            }
        }
        return;
    }
    if u >= d {
        let reduced = power_mod(u32::try_from(u).unwrap_or(u32::MAX), minpoly);
        for (t, c) in reduced.iter().enumerate() {
            if let Some(slot) = out.get_mut(t.saturating_mul(m).saturating_add(v)) {
                *slot += coefficient * c;
            }
        }
        return;
    }
    if let Some(slot) = out.get_mut(u.saturating_mul(m).saturating_add(v)) {
        *slot += coefficient;
    }
}

/// `n choose k` as an integer.
fn binomial(
    n: usize,
    k: usize,
) -> BigInt {
    let k = k.min(n.saturating_sub(k));
    let mut result = BigInt::one();
    for i in 0..k {
        result = result * BigInt::from(n.saturating_sub(i)) / BigInt::from(i.saturating_add(1));
    }
    result
}

/// Integer interpolation point `i`: `0, 1, -1, 2, -2, ...`.
fn interpolation_point(index: usize) -> BigRational {
    let step = index.saturating_add(1) / 2;
    let magnitude = BigInt::from(step);
    let value = if index.is_multiple_of(2) {
        magnitude
    } else {
        -magnitude
    };
    BigRational::from_integer(value)
}

/// Lagrange interpolation through `(xs, values)`.
#[allow(clippy::single_call_fn)] // one call site
fn interpolate(
    xs: &[BigRational],
    values: &[BigRational],
) -> Option<Vec<BigRational>> {
    let n = xs.len();
    if n == 0 || values.len() != n {
        return None;
    }
    let mut total = vec![BigRational::zero()];
    for (i, value) in values.iter().enumerate() {
        let mut basis = vec![BigRational::one()];
        let mut denominator = BigRational::one();
        for (j, x) in xs.iter().enumerate() {
            if i == j {
                continue;
            }
            let length = basis.len().saturating_add(1);
            let mut next = vec![BigRational::zero(); length];
            for (k, c) in basis.iter().enumerate() {
                if let Some(slot) = next.get_mut(k) {
                    *slot -= c * x;
                }
                if let Some(slot) = next.get_mut(k.saturating_add(1)) {
                    *slot += c;
                }
            }
            basis = next;
            denominator *= xs.get(i)? - x;
        }
        if denominator.is_zero() {
            return None;
        }
        let factor = value / denominator;
        total = alg::add(&total, &alg::scale(&basis, &factor));
    }
    Some(total)
}

/// Gauss-Jordan solution of `matrix * solution = rhs` over `Q`.
#[allow(clippy::single_call_fn)] // one call site
fn solve(
    matrix: &[Vec<BigRational>],
    rhs: &[BigRational],
    columns: usize,
) -> Option<Vec<BigRational>> {
    let rows = matrix.len();
    let mut augmented: Vec<Vec<BigRational>> = Vec::with_capacity(rows);
    for (i, row) in matrix.iter().enumerate() {
        let mut line = row.clone();
        line.push(rhs.get(i)?.clone());
        augmented.push(line);
    }
    let mut pivot_row = 0_usize;
    let mut pivots: Vec<usize> = Vec::new();
    for column in 0..columns {
        let found = (pivot_row..rows).find(|&r| {
            augmented
                .get(r)
                .and_then(|row| row.get(column))
                .is_some_and(|v| !v.is_zero())
        });
        let Some(found) = found else {
            continue;
        };
        augmented.swap(pivot_row, found);
        let pivot = augmented.get(pivot_row)?.get(column)?.clone();
        if let Some(row) = augmented.get_mut(pivot_row) {
            for value in row.iter_mut() {
                *value /= &pivot;
            }
        }
        for r in 0..rows {
            if r == pivot_row {
                continue;
            }
            let factor = augmented.get(r)?.get(column)?.clone();
            if factor.is_zero() {
                continue;
            }
            let pivot_values = augmented.get(pivot_row)?.clone();
            let row = augmented.get_mut(r)?;
            for c in 0..=columns {
                let value = pivot_values.get(c)?;
                if let Some(slot) = row.get_mut(c) {
                    *slot -= &factor * value;
                }
            }
        }
        pivots.push(column);
        pivot_row = pivot_row.saturating_add(1);
        if pivot_row == rows {
            break;
        }
    }
    // Consistency: a zero row with a non-zero right-hand side.
    for row in augmented.iter().skip(pivot_row) {
        let zero = row.iter().take(columns).all(Zero::is_zero);
        let rhs = row.get(columns).cloned().unwrap_or_else(BigRational::zero);
        if zero && !rhs.is_zero() {
            return None;
        }
    }
    let mut solution = vec![BigRational::zero(); columns];
    for (column, value) in pivots.iter().zip(augmented.iter()) {
        if let (Some(slot), Some(value)) = (solution.get_mut(*column), value.get(columns)) {
            *slot = value.clone();
        }
    }
    Some(solution)
}

/// Whether the value of `theta_expr` at the new embedding is the old
/// `theta`, decided exactly.
fn embedding_matches(
    new_field: &Field,
    theta_expr: &Elem,
    old_root: &RealAlg,
) -> bool {
    let (old_lo, old_hi) = old_root.interval();
    if old_lo == old_hi {
        let difference = theta_expr.sub(&Elem::rational(old_lo.clone()));
        return new_field.clone().sign(&difference) == Some(0);
    }
    let mut field = new_field.clone();
    for _ in 0..MAX_BISECT {
        let (lo, hi) = field.enclosure(theta_expr);
        if lo > *old_lo && hi < *old_hi {
            return true;
        }
        if hi < *old_lo || lo > *old_hi {
            return false;
        }
        if field.generator_value().is_some() {
            return false;
        }
        field.root_mut_refine_once();
    }
    false
}

/// Whether the value of `value` at the new embedding lies in `interval`.
fn enclosure_inside(
    new_field: &Field,
    value: &Elem,
    interval: &(BigRational, BigRational),
) -> bool {
    let (target_lo, target_hi) = interval;
    if target_lo == target_hi {
        let difference = value.sub(&Elem::rational(target_lo.clone()));
        return new_field.clone().sign(&difference) == Some(0);
    }
    let mut field = new_field.clone();
    for _ in 0..MAX_BISECT {
        let (lo, hi) = field.enclosure(value);
        if lo > *target_lo && hi < *target_hi {
            return true;
        }
        if hi < *target_lo || lo > *target_hi {
            return false;
        }
        if field.generator_value().is_some() {
            return false;
        }
        field.root_mut_refine_once();
    }
    false
}

impl Field {
    /// One bisection step of the generator interval.
    pub(crate) fn root_mut_refine_once(&mut self) {
        self.root.refine_once();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(
        n: i64,
        d: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(n), BigInt::from(d))
    }

    fn field_sqrt2() -> Field {
        let poly = vec![q(-2, 1), q(0, 1), q(1, 1)];
        Field::from_root(&poly, q(1, 1), q(2, 1)).unwrap_or_else(|| panic!())
    }

    #[test]
    fn field_arithmetic() {
        let mut field = field_sqrt2();
        let t = field.theta();
        let two = Elem::rational(q(2, 1));
        let square = field.mul(&t, &t);
        assert_eq!(field.sign(&square.sub(&two)), Some(0));
        let inverse = field.inv(&t).unwrap_or_else(|| panic!());
        assert_eq!(
            field.sign(&field.mul(&t, &inverse).sub(&Elem::one())),
            Some(0)
        );
        let three = Elem::rational(q(3, 1));
        assert_eq!(field.sign(&three), Some(1));
        assert_eq!(field.sign(&three.neg()), Some(-1));
    }

    #[test]
    fn root_isolation_over_a_field() {
        let mut field = field_sqrt2();
        // x^2 - theta, theta = sqrt(2): roots ±2^(1/4).
        let t = field.theta();
        let p = vec![t.neg(), Elem::zero(), Elem::one()];
        let roots = f_isolate(&p, &mut field).unwrap_or_default();
        assert_eq!(roots.len(), 2);
        let first = roots.first().cloned().unwrap_or_default();
        assert!(first.1 < q(0, 1), "{first:?}");
        let second = roots.get(1).cloned().unwrap_or_default();
        assert!(second.0 > q(1, 1), "{second:?}");
    }

    #[test]
    fn extensions() {
        let mut field = field_sqrt2();
        let t = field.theta();
        let p = vec![t.neg(), Elem::zero(), Elem::one()];
        let roots = f_isolate(&p, &mut field).unwrap_or_default();
        let interval = roots.first().cloned().unwrap_or_default();
        let mut extension = extend(&mut field, &p, &interval).unwrap_or_else(|| panic!());
        // theta = sqrt(2), c = 2^(1/4), and c^2 = theta.
        let squared = extension.field.mul(&extension.root, &extension.root);
        assert!(squared.sub(&extension.theta).is_zero());
        assert!(extension.field.sign(&extension.theta).unwrap_or(0) > 0);
        assert!(extension.field.sign(&extension.root).unwrap_or(0) < 0);
    }
}
