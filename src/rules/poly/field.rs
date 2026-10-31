//! Coefficient domains for the Gröbner basis engines.
//!
//! [`Domain`] names the field a computation runs over:
//!
//! * [`Domain::Rational`] — the rational numbers `Q`, the classical setting of
//!   [`groebner`](super::groebner::groebner), F4 and F5;
//! * [`Domain::Prime`] — a prime field `GF(p)`, with `p` an arbitrary-precision
//!   prime. Coefficients are kept as the **canonical integers `0 .. p`** (never
//!   negative, never `p` or larger), so the output of a computation is already
//!   in the representative convention documented for the operators.
//!
//! The arithmetic of [`Domain::Prime`] is modular arithmetic — `(a + b) mod p`,
//! `(a b) mod p`, and a modular inverse [`Domain::inv`] for monic normalisation
//! and for every division in the elimination. Nothing is lifted to `Q` and
//! reduced afterwards. `p` is validated to be prime when the domain is built:
//! prime fields have no zero divisors, so every non-zero element is invertible
//! and the elimination never meets a pivot it cannot divide by. Composite
//! moduli are refused instead of silently producing an incomplete answer.
//!
//! Powers `GF(p^k)` with `k > 1` are **out of scope**: the engines take a
//! modulus, not a modulus polynomial. The existing extension-field operators
//! of `rules::discrete::finite_field` (`gfx_*`) are the place to compute in
//! `GF(p^k)`; a Gröbner basis over an extension field would need a
//! coefficient-domain implementation for it.
//!
//! # Rational coefficients modulo `p`
//!
//! A polynomial read by the operators may have rational coefficients, e.g.
//! `x^2 - 1/2`. Over `GF(p)` such a coefficient is mapped to
//! `a * b^-1 mod p` (the natural image of `a/b` in the prime field), and a
//! coefficient whose denominator is divisible by `p` has no image and makes the
//! request unevaluated. For example, over `GF(5)` the polynomial `x^2 - 1/2`
//! becomes `x^2 - 3` because `2 * 3 = 1 (mod 5)`.
//!
//! # Limits
//!
//! The modulus is arbitrary-precision, but the F5 engine rejects coefficients
//! whose canonical representative needs more than `32_768` bits (its
//! `MAX_BITS` resource guard, inherited from the rational path). Primality
//! uses the crate's BPSW test, so a huge `p` costs one probable-prime test up
//! front.

use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::rules::discrete::mod_inverse;
use crate::rules::discrete::modulo;
use crate::rules::number_theory::is_prime;

/// The coefficient field of a Gröbner basis computation.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Domain {
    /// The field of rational numbers.
    #[default]
    Rational,
    /// The prime field `GF(p)`; the modulus is a prime larger than one.
    Prime(BigInt),
}

impl Domain {
    /// The field of rational numbers.
    #[must_use]
    pub const fn rational() -> Self {
        Self::Rational
    }

    /// The prime field `GF(p)`, or `None` when `p` is not a prime larger than
    /// one.
    #[must_use]
    pub fn prime(p: BigInt) -> Option<Self> {
        (p > BigInt::one() && is_prime(&p)).then_some(Self::Prime(p))
    }

    /// The modulus of a prime field.
    #[must_use]
    pub const fn modulus(&self) -> Option<&BigInt> {
        match self {
            | Self::Rational => None,
            | Self::Prime(p) => Some(p),
        }
    }

    /// Whether this is the rational field.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        matches!(self, Self::Rational)
    }

    /// The zero of the field, in canonical form.
    #[must_use]
    pub fn zero(&self) -> BigRational {
        BigRational::zero()
    }

    /// The one of the field, in canonical form.
    #[must_use]
    pub fn one(&self) -> BigRational {
        BigRational::one()
    }

    /// Whether a coefficient is zero.
    #[must_use]
    pub fn is_zero(
        &self,
        coefficient: &BigRational,
    ) -> bool {
        coefficient.is_zero()
    }

    /// Whether a coefficient is one.
    #[must_use]
    pub fn is_one(
        &self,
        coefficient: &BigRational,
    ) -> bool {
        coefficient.is_one()
    }

    /// The canonical representative of `coefficient` in this domain.
    ///
    /// For the rationals this is the coefficient itself. For `GF(p)` it is the
    /// integer `a b^-1 mod p` in `0 .. p`; `None` when the denominator is
    /// divisible by `p` (the coefficient has no image in the prime field).
    #[must_use]
    pub fn normalize(
        &self,
        coefficient: &BigRational,
    ) -> Option<BigRational> {
        match self {
            | Self::Rational => Some(coefficient.clone()),
            | Self::Prime(p) => Self::residue(coefficient, p).map(BigRational::from_integer),
        }
    }

    /// The canonical residue of a rational modulo `p`.
    fn residue(
        coefficient: &BigRational,
        p: &BigInt,
    ) -> Option<BigInt> {
        let numer = modulo(coefficient.numer(), p);
        if coefficient.denom().is_one() {
            return Some(numer);
        }
        let denom = modulo(coefficient.denom(), p);
        let inverse = mod_inverse(&denom, p)?;
        Some(modulo(&(numer * inverse), p))
    }

    /// The canonical integer behind a coefficient.
    ///
    /// Every coefficient that the engines produce has denominator one: the
    /// constructors normalise their input, and the operations below return
    /// canonical representatives.
    fn integer(coefficient: &BigRational) -> BigInt {
        debug_assert!(
            coefficient.denom().is_one(),
            "a canonical GF(p) coefficient is an integer"
        );
        coefficient.numer().clone()
    }

    /// `a + b` in the field.
    #[must_use]
    pub fn add(
        &self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        match self {
            | Self::Rational => a + b,
            | Self::Prime(p) => {
                BigRational::from_integer(modulo(&(Self::integer(a) + Self::integer(b)), p))
            },
        }
    }

    /// `a - b` in the field.
    #[must_use]
    pub fn sub(
        &self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        match self {
            | Self::Rational => a - b,
            | Self::Prime(p) => {
                BigRational::from_integer(modulo(&(Self::integer(a) - Self::integer(b)), p))
            },
        }
    }

    /// `a * b` in the field.
    #[must_use]
    pub fn mul(
        &self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        match self {
            | Self::Rational => a * b,
            | Self::Prime(p) => {
                BigRational::from_integer(modulo(&(Self::integer(a) * Self::integer(b)), p))
            },
        }
    }

    /// `-a` in the field.
    #[must_use]
    pub fn neg(
        &self,
        a: &BigRational,
    ) -> BigRational {
        match self {
            | Self::Rational => -a,
            | Self::Prime(p) => {
                BigRational::from_integer(modulo(&(BigInt::zero() - Self::integer(a)), p))
            },
        }
    }

    /// The inverse of `a`, or `None` for zero (and, over a prime field, for an
    /// element with no inverse, which cannot occur for a valid domain).
    #[must_use]
    pub fn inv(
        &self,
        a: &BigRational,
    ) -> Option<BigRational> {
        match self {
            | Self::Rational => {
                if a.is_zero() {
                    None
                } else {
                    Some(BigRational::one() / a)
                }
            },
            | Self::Prime(p) => {
                let residue = Self::residue(a, p)?;
                if residue.is_zero() {
                    return None;
                }
                mod_inverse(&residue, p).map(BigRational::from_integer)
            },
        }
    }

    /// `a / b` in the field, `None` when `b` is zero.
    #[must_use]
    pub fn div(
        &self,
        a: &BigRational,
        b: &BigRational,
    ) -> Option<BigRational> {
        Some(self.mul(a, &self.inv(b)?))
    }
}

impl fmt::Display for Domain {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Rational => f.write_str("Q"),
            | Self::Prime(p) => write!(f, "GF({p})"),
        }
    }
}

#[cfg(test)]
pub(crate) mod oracle {
    //! An independent Buchberger implementation over `GF(p)`, used only by the
    //! tests of the Gröbner engines.
    //!
    //! It shares no code with F4, F5, `groebner` or [`super::Domain`]: the
    //! coefficients are plain `i64` residues, the polynomials are `Vec` term
    //! lists, and the monomial order comparison is written here. A bug that
    //! affects both the engine and this reference would have to be made twice,
    //! independently, which is exactly what the cross-check is meant to rule
    //! out. It is deliberately naive (all pairs, no criteria) and only for
    //! small systems.

    use std::cmp::Ordering;

    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::One;

    use crate::rules::poly::groebner::GPoly;
    use crate::rules::poly::groebner::Order;
    use crate::rules::poly::groebner::Term;

    /// A term of the reference polynomial: canonical residue and exponents.
    type FpTerm = (Vec<u32>, i64);

    fn order_cmp(
        order: Order,
        a: &[u32],
        b: &[u32],
    ) -> Ordering {
        let degree = |m: &[u32]| m.iter().map(|&e| i64::from(e)).sum::<i64>();
        match order {
            | Order::Lex => a.cmp(b),
            | Order::GradedLex => degree(a).cmp(&degree(b)).then_with(|| a.cmp(b)),
            | Order::GradedRevLex => {
                degree(a).cmp(&degree(b)).then_with(|| {
                    for i in (0..a.len().max(b.len())).rev() {
                        let (x, y) = (
                            a.get(i).copied().unwrap_or(0),
                            b.get(i).copied().unwrap_or(0),
                        );
                        if x != y {
                            return y.cmp(&x);
                        }
                    }
                    Ordering::Equal
                })
            },
        }
    }

    fn inverse(
        a: i64,
        p: i64,
    ) -> i64 {
        let a = a.rem_euclid(p);
        for candidate in 1..p {
            if (a * candidate).rem_euclid(p) == 1 {
                return candidate;
            }
        }
        0
    }

    fn normalize(
        mut terms: Vec<FpTerm>,
        p: i64,
        order: Order,
    ) -> Vec<FpTerm> {
        for (_, coefficient) in &mut terms {
            *coefficient = coefficient.rem_euclid(p);
        }
        terms.retain(|(_, coefficient)| *coefficient != 0);
        terms.sort_by(|a, b| order_cmp(order, &b.0, &a.0));
        let mut merged: Vec<FpTerm> = Vec::with_capacity(terms.len());
        for (monomial, coefficient) in terms {
            match merged.last_mut() {
                | Some((last, accumulated)) if *last == monomial => {
                    *accumulated = (*accumulated + coefficient).rem_euclid(p);
                },
                | _ => merged.push((monomial, coefficient)),
            }
        }
        merged.retain(|(_, coefficient)| *coefficient != 0);
        merged
    }

    fn leading(p: &[FpTerm]) -> Option<&FpTerm> {
        p.first()
    }

    fn divides(
        a: &[u32],
        b: &[u32],
    ) -> bool {
        a.iter().zip(b).all(|(x, y)| x <= y)
    }

    fn lcm(
        a: &[u32],
        b: &[u32],
    ) -> Vec<u32> {
        a.iter().zip(b).map(|(&x, &y)| x.max(y)).collect()
    }

    fn shift(
        a: &[u32],
        b: &[u32],
    ) -> Vec<u32> {
        a.iter()
            .zip(b)
            .map(|(&x, &y)| x.saturating_sub(y))
            .collect()
    }

    fn multiply(
        a: &[u32],
        b: &[u32],
    ) -> Vec<u32> {
        let n = a.len().max(b.len());
        (0..n)
            .map(|i| {
                a.get(i)
                    .copied()
                    .unwrap_or(0)
                    .saturating_add(b.get(i).copied().unwrap_or(0))
            })
            .collect()
    }

    fn monic(
        p: &[FpTerm],
        modulus: i64,
        order: Order,
    ) -> Vec<FpTerm> {
        let Some((_, lead)) = leading(p) else {
            return Vec::new();
        };
        let scale = inverse(*lead, modulus);
        normalize(
            p.iter()
                .map(|(m, c)| (m.clone(), (c * scale).rem_euclid(modulus)))
                .collect(),
            modulus,
            order,
        )
    }

    /// `p - factor * x^shift * q`, term by term in the reference arithmetic.
    fn subtract_scaled(
        p: &[FpTerm],
        q: &[FpTerm],
        shift_by: &[u32],
        factor: i64,
        modulus: i64,
        order: Order,
    ) -> Vec<FpTerm> {
        let mut terms: Vec<FpTerm> = p.to_vec();
        for (monomial, coefficient) in q {
            terms.push((
                multiply(monomial, shift_by),
                (-(coefficient * factor)).rem_euclid(modulus),
            ));
        }
        normalize(terms, modulus, order)
    }

    /// The normal form of `p` modulo `basis`: every term is reduced, and the
    /// terms that are irreducible modulo the leading terms of the basis are
    /// kept in the result.
    fn reduce(
        p: &[FpTerm],
        basis: &[Vec<FpTerm>],
        modulus: i64,
        order: Order,
    ) -> Vec<FpTerm> {
        let mut work = normalize(p.to_vec(), modulus, order);
        let mut done: Vec<FpTerm> = Vec::new();
        while let Some((monomial, coefficient)) = leading(&work).cloned() {
            let reducer = basis.iter().find_map(|g| {
                let (lead, lead_coefficient) = leading(g)?;
                divides(lead, &monomial).then_some((g, lead.clone(), *lead_coefficient))
            });
            match reducer {
                | Some((divisor, lead, lead_coefficient)) => {
                    let factor =
                        (coefficient * inverse(lead_coefficient, modulus)).rem_euclid(modulus);
                    work = subtract_scaled(
                        &work,
                        divisor,
                        &shift(&monomial, &lead),
                        factor,
                        modulus,
                        order,
                    );
                },
                | None => {
                    done.push((monomial, coefficient));
                    work.remove(0);
                },
            }
        }
        normalize(done, modulus, order)
    }

    fn s_polynomial(
        f: &[FpTerm],
        g: &[FpTerm],
        modulus: i64,
        order: Order,
    ) -> Vec<FpTerm> {
        let (Some((fm, fc)), Some((gm, gc))) = (leading(f), leading(g)) else {
            return Vec::new();
        };
        let common = lcm(fm, gm);
        let left = normalize(
            f.iter()
                .map(|(m, c)| {
                    (
                        multiply(m, &shift(&common, fm)),
                        (c * inverse(*fc, modulus)).rem_euclid(modulus),
                    )
                })
                .collect(),
            modulus,
            order,
        );
        let right = normalize(
            g.iter()
                .map(|(m, c)| {
                    (
                        multiply(m, &shift(&common, gm)),
                        (c * inverse(*gc, modulus)).rem_euclid(modulus),
                    )
                })
                .collect(),
            modulus,
            order,
        );
        subtract_scaled(&left, &right, &[], 1, modulus, order)
    }

    /// The reduced Gröbner basis of the ideal generated by `generators` over
    /// `GF(p)`, as canonical integer terms.
    fn buchberger(
        generators: &[Vec<FpTerm>],
        modulus: i64,
        order: Order,
    ) -> Vec<Vec<FpTerm>> {
        let mut basis: Vec<Vec<FpTerm>> = generators
            .iter()
            .map(|g| monic(g, modulus, order))
            .filter(|g| !g.is_empty())
            .collect();
        let mut pairs: Vec<(usize, usize)> = (0..basis.len())
            .flat_map(|i| (0..i).map(move |j| (j, i)))
            .collect();
        while let Some((i, j)) = pairs.pop() {
            let remainder = reduce(
                &s_polynomial(&basis[i], &basis[j], modulus, order),
                &basis,
                modulus,
                order,
            );
            if !remainder.is_empty() {
                let new = basis.len();
                pairs.extend((0..new).map(|k| (k, new)));
                basis.push(monic(&remainder, modulus, order));
            }
        }
        minimal_reduced(&basis, modulus, order)
    }

    fn minimal_reduced(
        basis: &[Vec<FpTerm>],
        modulus: i64,
        order: Order,
    ) -> Vec<Vec<FpTerm>> {
        let mut minimal: Vec<Vec<FpTerm>> = Vec::new();
        for (index, g) in basis.iter().enumerate() {
            let Some((lead, _)) = leading(g) else {
                continue;
            };
            let redundant = basis.iter().enumerate().any(|(other, h)| {
                other != index
                    && leading(h)
                        .is_some_and(|(m, _)| divides(m, lead) && (m != lead || other < index))
            });
            if !redundant {
                minimal.push(g.clone());
            }
        }
        let mut reduced: Vec<Vec<FpTerm>> = Vec::new();
        for (index, g) in minimal.iter().enumerate() {
            let others: Vec<Vec<FpTerm>> = minimal
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != index)
                .map(|(_, h)| h.clone())
                .collect();
            reduced.push(monic(&reduce(g, &others, modulus, order), modulus, order));
        }
        reduced.sort_by(|a, b| {
            match (leading(a), leading(b)) {
                | (Some((x, _)), Some((y, _))) => order_cmp(order, x, y),
                | _ => Ordering::Equal,
            }
        });
        reduced
    }

    fn to_gpoly(
        p: &[FpTerm],
        order: Order,
    ) -> GPoly {
        GPoly::new(
            p.iter()
                .map(|(m, c)| (m.clone(), BigRational::from_integer(BigInt::from(*c))))
                .collect::<Vec<Term>>(),
            order,
        )
    }

    /// The reduced Gröbner basis of `generators` over `GF(p)`, as [`GPoly`]s
    /// with canonical integer coefficients.
    ///
    /// `generators` are `(coefficient, exponents)` lists; coefficients are
    /// reduced modulo `p`, so negative values are allowed. The prime `p` must
    /// fit an `i64` and is assumed prime.
    pub(crate) fn buchberger_fp(
        modulus: i64,
        generators: &[Vec<(i64, Vec<u32>)>],
        order: Order,
    ) -> Vec<GPoly> {
        let input: Vec<Vec<FpTerm>> = generators
            .iter()
            .map(|g| {
                normalize(
                    g.iter().map(|(c, m)| (m.clone(), *c)).collect(),
                    modulus,
                    order,
                )
            })
            .collect();
        buchberger(&input, modulus, order)
            .iter()
            .map(|g| to_gpoly(g, order))
            .collect()
    }

    /// Whether `candidate` reduces to zero modulo the reference basis, i.e.
    /// lies in the ideal.
    pub(crate) fn reduces_to_zero_fp(
        modulus: i64,
        basis: &[GPoly],
        candidate: &[(i64, Vec<u32>)],
        order: Order,
    ) -> bool {
        let basis: Vec<Vec<FpTerm>> = basis
            .iter()
            .map(|g| {
                normalize(
                    g.terms()
                        .iter()
                        .map(|(m, c)| {
                            let coefficient = i64::try_from(c.numer().clone()).unwrap_or(0);
                            (m.clone(), coefficient)
                        })
                        .collect(),
                    modulus,
                    order,
                )
            })
            .collect();
        let candidate = normalize(
            candidate.iter().map(|(c, m)| (m.clone(), *c)).collect(),
            modulus,
            order,
        );
        reduce(&candidate, &basis, modulus, order).is_empty()
    }

    // ------------------------------------------------------------------
    // Test systems shared by the engine tests
    // ------------------------------------------------------------------

    /// The primes the finite-field tests cover: the small characteristic
    /// cases, an odd mid-sized prime and a larger prime.
    pub(crate) const PRIMES: [i64; 5] = [2, 3, 5, 7, 101];

    /// A deterministic linear congruential generator, so the "random" systems
    /// are identical on every machine and in every run.
    pub(crate) struct Lcg(u64);

    impl Lcg {
        pub(crate) fn new(seed: u64) -> Self {
            Self(seed)
        }

        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0 >> 33
        }

        fn coefficient(
            &mut self,
            modulus: i64,
        ) -> i64 {
            let span = u64::try_from(2 * modulus - 1).unwrap_or(1);
            i64::try_from(self.next() % span).unwrap_or(0) - (modulus - 1)
        }

        fn choose(
            &mut self,
            upper: usize,
        ) -> usize {
            usize::try_from(self.next() % u64::try_from(upper).unwrap_or(1)).unwrap_or(0)
        }
    }

    fn exponents(
        n: usize,
        powers: &[(usize, u32)],
    ) -> Vec<u32> {
        let mut out = vec![0_u32; n];
        for &(i, e) in powers {
            if let Some(slot) = out.get_mut(i) {
                *slot = slot.saturating_add(e);
            }
        }
        out
    }

    /// A dense system: every monomial of degree at most two with a random
    /// coefficient.
    pub(crate) fn random_dense(
        modulus: i64,
        seed: u64,
        n: usize,
        rows: usize,
    ) -> Vec<Vec<(i64, Vec<u32>)>> {
        let mut rng = Lcg::new(seed);
        let mut support: Vec<Vec<u32>> = Vec::new();
        for i in 0..n {
            support.push(exponents(n, &[(i, 1)]));
            support.push(exponents(n, &[(i, 2)]));
            for j in (i + 1)..n {
                support.push(exponents(n, &[(i, 1), (j, 1)]));
            }
        }
        let mut out = Vec::new();
        for _ in 0..rows {
            let mut terms: Vec<(i64, Vec<u32>)> = support
                .iter()
                .map(|mono| (rng.coefficient(modulus), mono.clone()))
                .filter(|(c, _)| c.rem_euclid(modulus) != 0)
                .collect();
            if terms.is_empty() {
                terms.push((1, exponents(n, &[(rng.choose(n), 2)])));
            }
            out.push(terms);
        }
        out
    }

    /// A sparse system: three random monomials per polynomial, degree at most
    /// two.
    pub(crate) fn random_sparse(
        modulus: i64,
        seed: u64,
        n: usize,
        rows: usize,
    ) -> Vec<Vec<(i64, Vec<u32>)>> {
        let mut rng = Lcg::new(seed);
        let mut out = Vec::new();
        for _ in 0..rows {
            let mut terms: Vec<(i64, Vec<u32>)> = Vec::new();
            for _ in 0..3 {
                let powers: Vec<(usize, u32)> = (0..2)
                    .map(|_| (rng.choose(n), u32::try_from(rng.choose(3)).unwrap_or(0)))
                    .collect();
                terms.push((rng.coefficient(modulus), exponents(n, &powers)));
            }
            if terms.iter().all(|(c, _)| c.rem_euclid(modulus) == 0) {
                terms.push((1, exponents(n, &[])));
            }
            out.push(terms);
        }
        out
    }

    /// Katsura-`m` in `m + 1` variables, over the integers.
    pub(crate) fn katsura(m: usize) -> Vec<Vec<(i64, Vec<u32>)>> {
        let n = m + 1;
        let mut out = Vec::new();
        let mut linear = vec![(1, exponents(n, &[(0, 1)])), (-1, exponents(n, &[]))];
        for l in 1..=m {
            linear.push((2, exponents(n, &[(l, 1)])));
        }
        out.push(linear);
        for k in 0..m {
            let mut terms = vec![(-1, exponents(n, &[(k, 1)]))];
            for l in -(m as i64)..=(m as i64) {
                let a = l.unsigned_abs() as usize;
                let b = (k as i64 - l).unsigned_abs() as usize;
                if a <= m && b <= m {
                    terms.push((1, exponents(n, &[(a, 1), (b, 1)])));
                }
            }
            out.push(terms);
        }
        out
    }

    /// Cyclic-`n`, over the integers.
    pub(crate) fn cyclic(n: usize) -> Vec<Vec<(i64, Vec<u32>)>> {
        let mut out = Vec::new();
        for d in 1..n {
            let terms = (0..n)
                .map(|start| {
                    let powers: Vec<(usize, u32)> = (0..d).map(|j| ((start + j) % n, 1)).collect();
                    (1, exponents(n, &powers))
                })
                .collect();
            out.push(terms);
        }
        out.push(vec![
            (1, exponents(n, &(0..n).map(|i| (i, 1)).collect::<Vec<_>>())),
            (-1, exponents(n, &[])),
        ]);
        out
    }

    /// The systems every engine test runs over every prime and monomial order:
    /// structured, redundant, degenerate and random cases.
    pub(crate) fn sample_systems(modulus: i64) -> Vec<(&'static str, Vec<Vec<(i64, Vec<u32>)>>)> {
        let n = 3;
        vec![
            (
                "circle and line",
                vec![
                    vec![
                        (1, exponents(n, &[(0, 2)])),
                        (1, exponents(n, &[(1, 2)])),
                        (-1, exponents(n, &[])),
                    ],
                    vec![(1, exponents(n, &[(0, 1)])), (-1, exponents(n, &[(1, 1)]))],
                ],
            ),
            ("katsura-2", katsura(2)),
            ("cyclic-3", cyclic(3)),
            (
                "equal leading monomials",
                vec![
                    vec![
                        (1, exponents(n, &[(0, 2)])),
                        (1, exponents(n, &[(1, 1)])),
                        (1, exponents(n, &[])),
                    ],
                    vec![
                        (1, exponents(n, &[(0, 2)])),
                        (1, exponents(n, &[(2, 1)])),
                        (-1, exponents(n, &[])),
                    ],
                ],
            ),
            (
                "non-monic and redundant",
                vec![
                    vec![(2, exponents(n, &[(0, 2)])), (-3, exponents(n, &[]))],
                    vec![(4, exponents(n, &[(0, 2)])), (-6, exponents(n, &[]))],
                    vec![
                        (2, exponents(n, &[(0, 1), (1, 1)])),
                        (-2, exponents(n, &[(1, 1)])),
                    ],
                ],
            ),
            (
                "unit ideal",
                vec![
                    vec![(1, exponents(n, &[(0, 1)]))],
                    vec![(1, exponents(n, &[(0, 1)])), (1, exponents(n, &[]))],
                ],
            ),
            (
                "random dense",
                random_dense(modulus, 0x2545_F491_4F6C_DD1D, n, 3),
            ),
            (
                "random sparse",
                random_sparse(modulus, 0x9E37_79B9_7F4A_7C15, n, 3),
            ),
        ]
    }

    /// Builds canonical `GF(p)` polynomials from a raw system.
    pub(crate) fn generators(
        modulus: i64,
        system: &[Vec<(i64, Vec<u32>)>],
        order: Order,
    ) -> Vec<GPoly> {
        let domain = super::Domain::prime(BigInt::from(modulus))
            .unwrap_or_else(|| panic!("{modulus} must be a prime"));
        system
            .iter()
            .map(|terms| {
                GPoly::new_in(
                    &domain,
                    terms
                        .iter()
                        .map(|(coefficient, exponents)| {
                            (
                                exponents.clone(),
                                BigRational::from_integer(BigInt::from(
                                    coefficient.rem_euclid(modulus),
                                )),
                            )
                        })
                        .collect(),
                    order,
                )
            })
            .collect()
    }

    fn from_gpoly(g: &GPoly) -> Vec<FpTerm> {
        g.terms()
            .iter()
            .map(|(monomial, coefficient)| {
                (
                    monomial.clone(),
                    i64::try_from(coefficient.numer().clone()).unwrap_or(0),
                )
            })
            .collect()
    }

    /// Checks the properties every engine's `GF(p)` output must have:
    /// canonical coefficients in `0 .. p`, a monic leading term, a reduced
    /// basis (idempotent under the reference implementation), and every pair
    /// of elements forming an S-polynomial that reduces to zero under the
    /// **reference** arithmetic.
    pub(crate) fn check_basis(
        domain: &super::Domain,
        modulus: i64,
        basis: &[GPoly],
        order: Order,
    ) {
        for g in basis {
            assert_eq!(
                g.leading().map(|(_, c)| c),
                Some(&BigRational::from_integer(BigInt::from(1))),
                "a GF({modulus}) basis element is not monic"
            );
            for (_, coefficient) in g.terms() {
                let value = coefficient.numer();
                assert!(
                    *value >= BigInt::from(0) && *value < BigInt::from(modulus),
                    "coefficient {coefficient} is not the canonical representative mod {modulus}"
                );
                assert!(
                    coefficient.denom().is_one(),
                    "a GF({modulus}) coefficient is not integral"
                );
            }
        }
        let raw: Vec<Vec<FpTerm>> = basis.iter().map(from_gpoly).collect();
        for (i, f) in raw.iter().enumerate() {
            for g in raw.iter().skip(i.saturating_add(1)) {
                let s = s_polynomial(f, g, modulus, order);
                assert!(
                    reduce(&s, &raw, modulus, order).is_empty(),
                    "an S-polynomial of the GF({modulus}) basis does not reduce to zero"
                );
            }
        }
        assert_eq!(
            crate::rules::poly::groebner::groebner_in(
                domain,
                basis,
                order,
                crate::rules::poly::groebner::GroebnerLimits::default()
            )
            .as_deref(),
            Some(basis),
            "the GF({modulus}) basis is not reduced"
        );
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::*;

    fn rational(
        numer: i64,
        denom: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(numer), BigInt::from(denom))
    }

    #[test]
    fn prime_domains_are_validated() {
        assert!(Domain::prime(BigInt::from(2)).is_some());
        assert!(Domain::prime(BigInt::from(101)).is_some());
        assert!(Domain::prime(BigInt::from(1)).is_none());
        assert!(Domain::prime(BigInt::from(0)).is_none());
        assert!(Domain::prime(BigInt::from(-7)).is_none());
        assert!(Domain::prime(BigInt::from(4)).is_none());
        assert!(Domain::prime(BigInt::from(1001)).is_none(), "7*11*13");
        assert_eq!(
            Domain::prime(BigInt::from(7)).map(|d| d.to_string()),
            Some("GF(7)".to_owned())
        );
        assert!(Domain::rational().is_rational());
        assert_eq!(Domain::default(), Domain::Rational);
        assert_eq!(
            Domain::prime(BigInt::from(5)).unwrap().modulus(),
            Some(&BigInt::from(5))
        );
    }

    #[test]
    fn modular_arithmetic_matches_the_residue_class() {
        let domain = Domain::prime(BigInt::from(7)).unwrap();
        let a = rational(5, 1);
        let b = rational(4, 1);
        assert_eq!(domain.add(&a, &b), rational(2, 1), "5 + 4 = 2 mod 7");
        assert_eq!(domain.sub(&a, &b), rational(1, 1));
        assert_eq!(domain.sub(&b, &a), rational(6, 1), "4 - 5 = -1 = 6");
        assert_eq!(domain.mul(&a, &b), rational(6, 1));
        assert_eq!(domain.neg(&a), rational(2, 1));
        assert_eq!(domain.inv(&a), Some(rational(3, 1)), "5*3 = 1 mod 7");
        assert_eq!(domain.div(&a, &b), Some(rational(3, 1)), "5/4 = 5*2");
        assert_eq!(domain.inv(&domain.zero()), None);
        assert_eq!(domain.div(&a, &domain.zero()), None);
        assert_eq!(domain.neg(&domain.zero()), domain.zero());
        // Canonical representatives stay in 0..p under repeated operations.
        let mut value = domain.one();
        for _ in 0..100 {
            value = domain.add(&domain.mul(&value, &a), &b);
            let residue = value.numer().clone();
            assert!(residue >= BigInt::zero());
            assert!(residue < BigInt::from(7));
        }
        // 1/2 = 4 mod 7 and 3/4 = 6 mod 7
        assert_eq!(domain.normalize(&rational(1, 2)), Some(rational(4, 1)));
        assert_eq!(domain.normalize(&rational(3, 4)), Some(rational(6, 1)));
        // A coefficient whose denominator is 0 mod p has no image.
        assert_eq!(domain.normalize(&rational(1, 7)), None);
        assert_eq!(domain.normalize(&rational(2, 14)), None);
        // GF(2) has one non-zero element, and -1 = 1.
        let two = Domain::prime(BigInt::from(2)).unwrap();
        let one = two.one();
        assert_eq!(two.neg(&one), one);
        assert_eq!(two.sub(&two.zero(), &one), one);
        assert_eq!(two.inv(&one), Some(one));
    }

    #[test]
    fn rational_arithmetic_is_unchanged() {
        let domain = Domain::rational();
        let a = rational(1, 2);
        let b = rational(1, 3);
        assert_eq!(domain.add(&a, &b), rational(5, 6));
        assert_eq!(domain.sub(&a, &b), rational(1, 6));
        assert_eq!(domain.mul(&a, &b), rational(1, 6));
        assert_eq!(domain.div(&a, &b), Some(rational(3, 2)));
        assert_eq!(domain.inv(&domain.zero()), None);
        assert_eq!(domain.normalize(&a), Some(a));
    }
}
