//! Galois groups of irreducible rational polynomials of degree 6 and 7,
//! block-tower data for solvable sextics and Gauss-period (cyclotomic)
//! forms for cyclic ones.
//!
//! **Group.** The polynomial is brought to a monic integer form `g(z)` (shift
//! and scale, as for quintics) and its roots are computed to about 190
//! digits. For a transitive group `H ≤ S_n` (tables in [`super::transitive`])
//! and an ordering `σ` of the roots, the orbit sum
//! `F_H(σ) = Σ_{h ∈ H} Π_i r_{σ(h(i))}^{e_i}` (distinct exponents `e_i`) is
//! an algebraic integer whose stabiliser in `S_n` is `H`; it is a rational
//! integer exactly when the Galois group lies in the conjugate `H^σ`
//! (provided no other coset gives the same value, which is checked, and a
//! second choice of exponents is used when it fails). This is a relative
//! resolvent test over all cosets `S_n / H`, so no subgroup lattice is
//! needed: the Galois group is the first group in the table, ordered by
//! size, that passes. The roots are then labelled by the positions `0..n`
//! on which the Galois group acts as the table group.
//!
//! **Towers.** A solvable sextic has a system of blocks of size 2 or 3. If
//! `θ_B` is a symmetric function of the roots in block `B`, its conjugates
//! are the `n/b` numbers `θ_B`; their minimal polynomial `Θ` is an integer
//! polynomial of degree `n/b`, and the block polynomial
//! `h_B(z) = Π_{i ∈ B} (z - r_i)` has coefficients that are polynomials in
//! `θ_B` with rational coefficients, found by solving a Vandermonde system
//! numerically and recognising rationals by continued fractions (checked to
//! 2^-250). Solving `Θ` and then `h` over `Q(θ)` is the tower.
//!
//! **Gauss periods.** A cyclic group of degree `n` with conductor a prime
//! `p ≡ 1 (mod 2n)` (all roots real) gives roots in the degree `n` subfield
//! of `Q(ζ_p)`: with the periods `η_k = Σ_{j ∈ g^k K} ζ_p^j` over the
//! subgroup `K` of index `n`, the roots are `Σ_t c_t η_{k+t}` with rational
//! `c`, found as the solution of a circulant system. Composite conductors
//! are not searched.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::galois::one;
use super::galois::refined_roots;
use super::galois::round_integer;
use super::galois::taylor_shift;
use super::galois::Fx;
use super::galois::BITS;
use super::transitive::transitive_groups;
use super::transitive::TransitiveGroup;

/// The Galois group of an irreducible polynomial together with the data
/// needed to work with its roots.
#[derive(Debug)]
pub struct Analysis {
    /// The group (a table entry).
    pub group: &'static TransitiveGroup,
    /// The monic integer polynomial `g(z)` (ascending) whose roots are the
    /// scaled, shifted roots of the input: `x = shift + z / scale`.
    pub monic: Vec<BigInt>,
    /// The shift of the depression.
    pub shift: BigRational,
    /// The scale `λ` (`z = λ (x - shift)`).
    pub scale: BigInt,
    /// `ordering[i]` is the index into the roots of the position `i`.
    pub ordering: Vec<usize>,
    pub(crate) roots: Vec<Fx>,
}

impl Analysis {
    /// Number of roots.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.monic.len().saturating_sub(1)
    }

    /// Approximate value of the root at `position` as `(re, im)` of `x`.
    #[must_use]
    pub fn root_approx(
        &self,
        position: usize,
    ) -> Option<(f64, f64)> {
        let index = *self.ordering.get(position)?;
        let z = self.roots.get(index)?;
        let to = |v: &BigInt| BigRational::new(v.clone(), one()).to_f64();
        let (re, im) = (to(&z.re)?, to(&z.im)?);
        let lambda = self.scale.to_f64()?;
        Some((self.shift.to_f64()? + re / lambda, im / lambda))
    }
}

/// `p(y + shift)` scaled to a monic integer polynomial `g(λ y)`.
fn integral_form(c: &[BigRational]) -> Option<(BigRational, BigInt, Vec<BigInt>)> {
    let n = c.len().checked_sub(1).filter(|&n| n >= 2)?;
    let lead = c.last()?.clone();
    if lead.is_zero() {
        return None;
    }
    let monic: Vec<BigRational> = c.iter().map(|a| a / &lead).collect();
    if monic.iter().all(BigRational::is_integer) {
        let g = monic.iter().map(BigRational::to_integer).collect();
        return Some((BigRational::zero(), BigInt::one(), g));
    }
    let shift = -monic.get(n - 1)? / BigRational::from_integer(BigInt::from(n));
    let depressed = taylor_shift(&monic, &shift);
    // Smallest scale λ with λ^(n-i) a_i integral, prime by prime.
    let mut scale = BigInt::one();
    for (i, a) in depressed.iter().enumerate().take(n - 1) {
        let mut denominator = a.denom().clone();
        let exponent = n - i;
        let mut prime = BigInt::from(2);
        while denominator > BigInt::one() {
            if (&denominator % &prime).is_zero() {
                let mut valuation = 0_usize;
                while (&denominator % &prime).is_zero() {
                    denominator /= &prime;
                    valuation += 1;
                }
                let needed = valuation.div_ceil(exponent);
                let mut have = 0_usize;
                let mut rest = scale.clone();
                while (&rest % &prime).is_zero() {
                    rest /= &prime;
                    have += 1;
                }
                for _ in have..needed {
                    scale *= &prime;
                }
            }
            prime += 1;
            if prime > BigInt::from(1_000_000) {
                return None;
            }
        }
    }
    let power = |k: usize| -> BigInt { (0..k).fold(BigInt::one(), |acc, _| acc * &scale) };
    let mut g = Vec::new();
    for (i, a) in depressed.iter().enumerate() {
        let value = a * BigRational::from_integer(power(n - i));
        if !value.is_integer() {
            return None;
        }
        g.push(value.to_integer());
    }
    Some((shift, scale, g))
}

/// Discriminant of the monic integer polynomial `g` (ascending), by the
/// Bareiss determinant of the Sylvester matrix of `g` and `g'`.
#[must_use]
pub fn discriminant(g: &[BigInt]) -> Option<BigInt> {
    let n = g.len().checked_sub(1).filter(|&n| n >= 2)?;
    let lead = g.last()?.clone();
    let derivative: Vec<BigInt> = g.iter().enumerate().skip(1).map(|(i, c)| c * BigInt::from(i)).collect();
    let size = 2 * n - 1;
    let mut matrix = vec![vec![BigInt::zero(); size]; size];
    // n-1 rows of g (descending coefficients), n rows of g'.
    for r in 0..n - 1 {
        for (k, c) in g.iter().rev().enumerate() {
            matrix.get_mut(r)?.get_mut(r + k)?.clone_from(c);
        }
    }
    for r in 0..n {
        for (k, c) in derivative.iter().rev().enumerate() {
            matrix.get_mut(n - 1 + r)?.get_mut(r + k)?.clone_from(c);
        }
    }
    let resultant = bareiss(matrix)?;
    let sign = if (n * (n - 1) / 2) % 2 == 0 { BigInt::one() } else { -BigInt::one() };
    Some(sign * resultant / lead)
}

fn bareiss(mut m: Vec<Vec<BigInt>>) -> Option<BigInt> {
    let n = m.len();
    let mut sign = BigInt::one();
    let mut previous = BigInt::one();
    for k in 0..n.saturating_sub(1) {
        if m.get(k)?.get(k)?.is_zero() {
            let swap = (k + 1..n).find(|&r| m.get(r).and_then(|row| row.get(k)).is_some_and(|v| !v.is_zero()))?;
            m.swap(k, swap);
            sign = -sign;
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let value = (m.get(i)?.get(j)? * m.get(k)?.get(k)? - m.get(i)?.get(k)? * m.get(k)?.get(j)?) / &previous;
                *m.get_mut(i)?.get_mut(j)? = value;
            }
        }
        previous.clone_from(m.get(k)?.get(k)?);
    }
    Some(sign * m.get(n - 1)?.get(n - 1)?.clone())
}

fn heap_permutations(
    k: usize,
    a: &mut [u8],
    out: &mut Vec<Vec<u8>>,
) {
    if k <= 1 {
        out.push(a.to_vec());
        return;
    }
    for i in 0..k {
        heap_permutations(k - 1, a, out);
        let j = if k.is_multiple_of(2) { i } else { 0 };
        a.swap(j, k - 1);
    }
}

fn all_permutations(n: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut current: Vec<u8> = (0..n).map(|i| u8::try_from(i).unwrap_or(0)).collect();
    heap_permutations(n, &mut current, &mut out);
    out
}

/// Exponent vectors for the orbit sums.
fn exponent_choices(n: usize) -> Vec<Vec<usize>> {
    let ascending: Vec<usize> = (0..n).collect();
    let mixed: Vec<usize> = (0..n).map(|i| (i * 3 + 1) % (n + 2)).collect();
    let spread: Vec<usize> = (0..n).map(|i| if i < 2 { i } else { 2 * i - 1 }).collect();
    vec![ascending, spread, mixed]
}

/// Largest power of two bound on the bit length of the real part of the
/// roots (to keep the resolvent values within the working precision).
fn small_enough(roots: &[Fx]) -> bool {
    let limit = u64::try_from(BITS).unwrap_or(0) + 12;
    roots.iter().all(|r| r.re.bits() < limit && r.im.bits() < limit)
}

/// The orbit sum of `h`-images of the monomial with exponents `e`, for the
/// ordering `sigma`.
fn orbit_sum(
    elements: &[Vec<u8>],
    sigma: &[u8],
    powers: &[Vec<Fx>],
    e: &[usize],
) -> Option<Fx> {
    let mut total = Fx::zero();
    for h in elements {
        let mut term = Fx::integer(&BigInt::one());
        for (i, &ei) in e.iter().enumerate() {
            if ei == 0 {
                continue;
            }
            let image = *sigma.get(usize::from(*h.get(i)?))?;
            let r = powers.get(usize::from(image))?.get(ei)?;
            term = term.mul(r);
        }
        total = total.add(&term);
    }
    Some(total)
}

/// The Galois group of the irreducible polynomial `c` and a labelling of
/// its roots.
///
/// The degree is 6 or 7 (ascending coefficients); the result is `None` for
/// other degrees, when the polynomial is not in range of the working
/// precision, or when the test is inconclusive.
#[must_use]
pub fn analyse(c: &[BigRational]) -> Option<Analysis> {
    let n = c.len().checked_sub(1)?;
    let groups = transitive_groups(n);
    if groups.is_empty() {
        return None;
    }
    let (shift, scale, monic) = integral_form(c)?;
    let roots = refined_roots(&monic)?;
    if roots.len() != n || !small_enough(&roots) {
        return None;
    }
    let max_exponent = 2 * n;
    // Tschirnhaus-transformed roots r + k r^2 break accidental symmetries
    // (roots that are multiples by roots of unity, say).
    let powers_by_transform: Vec<Vec<Vec<Fx>>> = [0_i64, 1, 3]
        .iter()
        .map(|&k| {
            roots
                .iter()
                .map(|r| {
                    let t = r.add(&r.mul(r).mul(&Fx::integer(&BigInt::from(k))));
                    let mut row = vec![Fx::integer(&BigInt::one())];
                    for _ in 1..=max_exponent {
                        let next = row.last().map_or_else(Fx::zero, |last| last.mul(&t));
                        row.push(next);
                    }
                    row
                })
                .collect()
        })
        .collect();
    let permutations = all_permutations(n);
    let magnitude_limit = u64::try_from(BITS).unwrap_or(0) + 120;
    for group in groups {
        let elements = group.elements();
        // Coset representatives: sigma minimal in its coset sigma.H.
        let representatives: Vec<&Vec<u8>> = permutations
            .iter()
            .filter(|sigma| {
                elements.iter().all(|h| {
                    let image: Vec<u8> = h.iter().map(|&i| sigma.get(usize::from(i)).copied().unwrap_or(0)).collect();
                    image.as_slice() >= sigma.as_slice()
                })
            })
            .collect();
        let mut ambiguous = false;
        let choices: Vec<(usize, Vec<usize>)> =
            (0..powers_by_transform.len()).flat_map(|t| exponent_choices(n).into_iter().map(move |e| (t, e))).collect();
        for (transform, e) in choices {
            let powers = powers_by_transform.get(transform)?;
            let mut values: Vec<Fx> = Vec::with_capacity(representatives.len());
            for sigma in &representatives {
                values.push(orbit_sum(&elements, sigma, powers, &e)?);
            }
            if values.iter().any(|v| v.re.bits() > magnitude_limit || v.im.bits() > magnitude_limit) {
                return None;
            }
                    let hits: Vec<usize> = values.iter().enumerate().filter(|(_, v)| round_integer(v).is_some()).map(|(i, _)| i).collect();
            if hits.is_empty() {
                break;
            }
            // Several cosets are integral when the normaliser is larger
            // than the group; take one whose value no other coset shares.
            let unique = hits.iter().copied().find(|&hit| {
                values.get(hit).is_some_and(|reference| values.iter().enumerate().all(|(i, v)| i == hit || !v.close(reference, 300)))
            });
                    let Some(hit) = unique else {
                ambiguous = true;
                continue;
            };
            let sigma = representatives.get(hit)?;
            return Some(Analysis {
                group,
                monic,
                shift,
                scale,
                ordering: sigma.iter().map(|&i| usize::from(i)).collect(),
                roots,
            });
        }
        let _ = ambiguous;
    }
    None
}

// ---------------------------------------------------------------------
// Linear algebra and rational recognition on fixed-point numbers.

/// Solves `A X = B` (complex fixed point) by Gaussian elimination.
fn solve_system(
    mut a: Vec<Vec<Fx>>,
    mut b: Vec<Vec<Fx>>,
) -> Option<Vec<Vec<Fx>>> {
    let n = a.len();
    for k in 0..n {
        let pivot = (k..n).max_by_key(|&r| {
            a.get(r).and_then(|row| row.get(k)).map_or_else(BigInt::zero, |v| v.re.abs() + v.im.abs())
        })?;
        a.swap(k, pivot);
        b.swap(k, pivot);
        let pivot_value = a.get(k)?.get(k)?.clone();
        for i in k + 1..n {
            let factor = a.get(i)?.get(k)?.div(&pivot_value)?;
            for j in k..n {
                let delta = factor.mul(a.get(k)?.get(j)?);
                let slot = a.get_mut(i)?.get_mut(j)?;
                *slot = slot.sub(&delta);
            }
            for j in 0..b.get(i)?.len() {
                let delta = factor.mul(b.get(k)?.get(j)?);
                let slot = b.get_mut(i)?.get_mut(j)?;
                *slot = slot.sub(&delta);
            }
        }
    }
    let columns = b.first()?.len();
    let mut x = vec![vec![Fx::zero(); columns]; n];
    for i in (0..n).rev() {
        for j in 0..columns {
            let mut value = b.get(i)?.get(j)?.clone();
            for k in i + 1..n {
                value = value.sub(&a.get(i)?.get(k)?.mul(x.get(k)?.get(j)?));
            }
            *x.get_mut(i)?.get_mut(j)? = value.div(a.get(i)?.get(i)?)?;
        }
    }
    Some(x)
}

/// The rational with small denominator within `2^-300` of the real
/// fixed-point number `re`, by continued fractions.
fn rationalize(z: &Fx) -> Option<BigRational> {
    if z.im.abs() > BigInt::one() << (BITS - 200) {
        return None;
    }
    let denominator = one();
    let (mut a, mut b) = (z.re.clone(), denominator.clone());
    let (mut h_prev, mut h) = (BigInt::one(), a.div_floor(&b));
    let (mut k_prev, mut k) = (BigInt::zero(), BigInt::one());
    let mut remainder = &a - &h * &b;
    let limit = BigInt::one() << (BITS - 300);
    for _ in 0..400 {
        // Error of h/k against re/2^BITS, scaled by k * 2^BITS.
        let error = (&z.re * &k - &h * &denominator).abs();
        if error < &k * &limit {
            return Some(BigRational::new(h, k));
        }
        if remainder.is_zero() {
            return None;
        }
        a = std::mem::replace(&mut b, remainder.clone());
        let q = a.div_floor(&b);
        remainder = &a - &q * &b;
        let h_next = &q * &h + &h_prev;
        let k_next = &q * &k + &k_prev;
        h_prev = std::mem::replace(&mut h, h_next);
        k_prev = std::mem::replace(&mut k, k_next);
        if k.bits() > 100 {
            return None;
        }
    }
    None
}

fn from_rational(q: &BigRational) -> Fx {
    Fx::real((q.numer() << BITS) / q.denom())
}

// ---------------------------------------------------------------------
// Towers.

/// A two-level decomposition `Q ⊂ Q(θ) ⊂ Q(roots)` of a solvable sextic.
#[derive(Clone, Debug)]
pub struct Tower {
    /// Block size (2 or 3).
    pub block_size: usize,
    /// The monic integer minimal polynomial `Θ(t)` of `θ` (ascending).
    pub theta: Vec<BigInt>,
    /// `block[m]` is the polynomial in `θ` (ascending coefficients) giving
    /// the coefficient of `z^m` in the block polynomial
    /// `z^b + Σ_{m<b} block[m](θ) z^m`.
    pub block: Vec<Vec<BigRational>>,
}

/// The towers of the block systems of a sextic's Galois group.
#[must_use]
pub fn towers(analysis: &Analysis) -> Vec<Tower> {
    let mut out = Vec::new();
    for system in analysis.group.block_systems() {
        if let Some(tower) = build_tower(analysis, &system) {
            out.push(tower);
        }
    }
    out
}

fn build_tower(
    analysis: &Analysis,
    system: &[Vec<usize>],
) -> Option<Tower> {
    let k = system.len();
    let b = system.first()?.len();
    let blocks: Vec<Vec<&Fx>> = system
        .iter()
        .map(|block| block.iter().map(|&p| analysis.ordering.get(p).and_then(|&i| analysis.roots.get(i))).collect::<Option<Vec<_>>>())
        .collect::<Option<_>>()?;
    let mut best: Option<(u64, Vec<Fx>, Vec<BigInt>)> = None;
    for (alpha, beta) in [(1_i64, 0_i64), (0, 1), (1, 1), (1, 2), (2, 1), (1, -1), (1, 3), (3, 1), (1, 4)] {
        let theta: Vec<Fx> = blocks
            .iter()
            .map(|block| {
                block.iter().fold(Fx::zero(), |acc, r| {
                    let term = r.mul(&Fx::integer(&BigInt::from(alpha))).add(&r.mul(r).mul(&Fx::integer(&BigInt::from(beta))));
                    acc.add(&term)
                })
            })
            .collect();
        let distinct = theta.iter().enumerate().all(|(i, a)| theta.iter().skip(i + 1).all(|c| !a.close(c, 300)));
        if !distinct {
            continue;
        }
        // Minimal polynomial: elementary symmetric functions.
        let mut poly = vec![Fx::integer(&BigInt::one())];
        for t in &theta {
            let mut next = vec![Fx::zero(); poly.len() + 1];
            for (i, a) in poly.iter().enumerate() {
                if let Some(slot) = next.get_mut(i + 1) {
                    *slot = slot.add(a);
                }
                if let Some(slot) = next.get_mut(i) {
                    *slot = slot.sub(&a.mul(t));
                }
            }
            poly = next;
        }
        let Some(integers) = poly.iter().map(round_integer).collect::<Option<Vec<BigInt>>>() else { continue };
        let height = integers.iter().map(|c| c.abs().bits()).max().unwrap_or(0);
        if best.as_ref().is_none_or(|(h, _, _)| height < *h) {
            best = Some((height, theta, integers));
        }
    }
    let (_, theta, integers) = best?;
    // Block polynomial coefficients as polynomials in theta.
    let vandermonde: Vec<Vec<Fx>> = theta
        .iter()
        .map(|t| {
            let mut row = vec![Fx::integer(&BigInt::one())];
            for _ in 1..k {
                let next = row.last().map_or_else(Fx::zero, |l| l.mul(t));
                row.push(next);
            }
            row
        })
        .collect();
    // e_m(B) for m < b: coefficient of z^m in prod (z - r).
    let rhs: Vec<Vec<Fx>> = blocks
        .iter()
        .map(|block| {
            let mut poly = vec![Fx::integer(&BigInt::one())];
            for r in block {
                let mut next = vec![Fx::zero(); poly.len() + 1];
                for (i, a) in poly.iter().enumerate() {
                    if let Some(slot) = next.get_mut(i + 1) {
                        *slot = slot.add(a);
                    }
                    if let Some(slot) = next.get_mut(i) {
                        *slot = slot.sub(&a.mul(r));
                    }
                }
                poly = next;
            }
            // z^0 .. z^(b-1)
            poly.into_iter().take(b).collect()
        })
        .collect();
    let solution = solve_system(vandermonde.clone(), rhs.clone())?;
    // solution[t][m]: coefficient of theta^t in the z^m coefficient.
    let mut block = Vec::new();
    for m in 0..b {
        let mut coefficients = Vec::new();
        for t in 0..k {
            coefficients.push(rationalize(solution.get(t)?.get(m)?)?);
        }
        block.push(coefficients);
    }
    // Verify the recognised rationals against the data.
    for (j, row) in vandermonde.iter().enumerate() {
        for (m, coefficients) in block.iter().enumerate() {
            let mut value = Fx::zero();
            for (t, q) in coefficients.iter().enumerate() {
                value = value.add(&from_rational(q).mul(row.get(t)?));
            }
            if !value.close(rhs.get(j)?.get(m)?, 250) {
                return None;
            }
        }
    }
    Some(Tower { block_size: b, theta: integers, block })
}

// ---------------------------------------------------------------------
// Gauss periods.

/// Roots of a cyclic polynomial in the degree `n` subfield of `Q(ζ_p)`:
/// `x = shift + z / scale` with `z = Σ_j coefficient_j 2 cos(2π j / p)`.
#[derive(Clone, Debug)]
pub struct Cyclotomic {
    /// The prime conductor `p`.
    pub conductor: u64,
    /// For each root: `(position, [(j, coefficient of 2cos(2πj/p))])`.
    pub roots: Vec<(usize, Vec<(u64, BigRational)>)>,
}

fn pi_fixed(precision: usize) -> BigInt {
    // Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    let atan_inverse = |x: u32| -> BigInt {
        let unit = BigInt::one() << precision;
        let x_squared = BigInt::from(x) * BigInt::from(x);
        let mut term = &unit / BigInt::from(x);
        let mut sum = term.clone();
        let mut k = 1_u32;
        while !term.is_zero() {
            term = &term / &x_squared;
            let piece = &term / BigInt::from(2 * k + 1);
            if k % 2 == 1 {
                sum -= piece;
            } else {
                sum += piece;
            }
            k += 1;
        }
        sum
    };
    BigInt::from(16) * atan_inverse(5) - BigInt::from(4) * atan_inverse(239)
}

/// `cos(2 π j / p)` in fixed point with `BITS` fractional bits.
fn cos_fixed(
    j: u64,
    p: u64,
    pi: &BigInt,
    guard: usize,
) -> BigInt {
    let j = j % p;
    let j = j.min(p - j);
    let precision = BITS + guard;
    // angle = 2 pi j / p
    let angle = (pi * BigInt::from(2 * j)) / BigInt::from(p);
    let unit = BigInt::one() << precision;
    let angle_squared = (&angle * &angle) >> precision;
    let mut term = unit.clone();
    let mut sum = unit;
    let mut k = 1_u32;
    loop {
        term = -((&term * &angle_squared) >> precision) / BigInt::from((2 * k - 1) * (2 * k));
        if term.is_zero() {
            break;
        }
        sum += &term;
        k += 1;
    }
    sum >> guard
}

fn primitive_root(p: u64) -> Option<u64> {
    let mut factors = Vec::new();
    let mut m = p - 1;
    let mut d = 2;
    while d * d <= m {
        if m.is_multiple_of(d) {
            factors.push(d);
            while m.is_multiple_of(d) {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        factors.push(m);
    }
    let pow_mod = |mut b: u64, mut e: u64| -> u64 {
        let mut r = 1_u64;
        b %= p;
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % p;
            }
            b = b * b % p;
            e >>= 1;
        }
        r
    };
    (2..p).find(|&g| factors.iter().all(|&f| pow_mod(g, (p - 1) / f) != 1))
}

fn is_prime(p: u64) -> bool {
    p >= 2 && (2..=p).take_while(|d| d * d <= p).all(|d| !p.is_multiple_of(d))
}

/// Gauss-period form of the roots of a cyclic polynomial with all roots
/// real, when its conductor is a prime `p < 4000`.
#[must_use]
pub fn cyclotomic(analysis: &Analysis) -> Option<Cyclotomic> {
    let n = analysis.degree();
    if !analysis.group.is_abelian() || n < 3 {
        return None;
    }
    let n_u64 = u64::try_from(n).ok()?;
    let real_tolerance = BigInt::one() << (BITS - 100);
    if analysis.roots.iter().any(|r| r.im.abs() > real_tolerance) {
        return None;
    }
    let disc = discriminant(&analysis.monic)?.abs();
    let pi = pi_fixed(BITS + 64);
    for p in (2 * n_u64 + 1..4000).step_by(usize::try_from(2 * n_u64).ok()?) {
        if !is_prime(p) || !(&disc % BigInt::from(p)).is_zero() {
            continue;
        }
        let g = primitive_root(p)?;
        let cos_table: Vec<BigInt> = (0..p).map(|j| cos_fixed(j, p, &pi, 64)).collect();
        let size = (p - 1) / n_u64; // |K|
        let power = |e: u64| -> u64 {
            let (mut r, mut b, mut e) = (1_u64, g % p, e);
            while e > 0 {
                if e & 1 == 1 {
                    r = r * b % p;
                }
                b = b * b % p;
                e >>= 1;
            }
            r
        };
        for e in (1..n_u64).filter(|e| e.gcd(&n_u64) == 1) {
            // Coset k of K = <g^n>: exponents g^(e k + n t).
            let coset = |k: u64| -> Vec<u64> { (0..size).map(|t| power(e * k + n_u64 * t)).collect() };
            let eta: Vec<BigInt> = (0..n_u64).map(|k| coset(k).iter().map(|&j| cos_table.get(usize::try_from(j).unwrap_or(0)).cloned().unwrap_or_default()).sum()).collect();
            let matrix: Vec<Vec<Fx>> = (0..n)
                .map(|k| (0..n).map(|t| eta.get((k + t) % n).cloned().map_or_else(Fx::zero, Fx::real)).collect())
                .collect();
            let rhs: Vec<Vec<Fx>> = (0..n)
                .map(|k| analysis.ordering.get(k).and_then(|&i| analysis.roots.get(i)).map(|r| vec![Fx::real(r.re.clone())]))
                .collect::<Option<_>>()?;
            let Some(solution) = solve_system(matrix.clone(), rhs.clone()) else { continue };
            let Some(c) = solution.iter().map(|row| row.first().and_then(rationalize)).collect::<Option<Vec<BigRational>>>() else {
                continue;
            };
            // Verify.
            let ok = (0..n).all(|k| {
                let mut value = Fx::zero();
                for (t, q) in c.iter().enumerate() {
                    if let Some(entry) = matrix.get(k).and_then(|row| row.get(t)) {
                        value = value.add(&from_rational(q).mul(entry));
                    }
                }
                rhs.get(k).and_then(|r| r.first()).is_some_and(|target| value.close(target, 250))
            });
            if !ok {
                continue;
            }
            // z_k = sum_t c_t eta_{k+t}; eta_m = sum_{j in coset m / ±} 2 cos(2 pi j/p).
            let mut roots = Vec::new();
            for k in 0..n_u64 {
                let mut terms: Vec<(u64, BigRational)> = Vec::new();
                for (t, q) in c.iter().enumerate() {
                    let m = (k + u64::try_from(t).ok()?) % n_u64;
                    for j in coset(m) {
                        let canonical = j.min(p - j);
                        if j > p - j {
                            continue;
                        }
                        match terms.iter_mut().find(|(jj, _)| *jj == canonical) {
                            | Some((_, slot)) => *slot += q,
                            | None => terms.push((canonical, q.clone())),
                        }
                    }
                }
                terms.retain(|(_, q)| !q.is_zero());
                terms.sort_by_key(|(j, _)| *j);
                roots.push((usize::try_from(k).ok()?, terms));
            }
            return Some(Cyclotomic { conductor: p, roots });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(c: &[i64]) -> Vec<BigRational> {
        c.iter().map(|&v| BigRational::from_integer(BigInt::from(v))).collect()
    }

    fn group_of(c: &[i64]) -> Option<(&'static str, usize)> {
        analyse(&poly(c)).map(|a| (a.group.label, a.group.order))
    }

    #[test]
    fn discriminants() {
        // x^3 - 2: -108; x^2 + 1 is not covered (degree >= 2 only): -4.
        let ints = |v: &[i64]| -> Vec<BigInt> { v.iter().map(|&c| BigInt::from(c)).collect() };
        assert_eq!(discriminant(&ints(&[-2, 0, 0, 1])), Some(BigInt::from(-108)));
        assert_eq!(discriminant(&ints(&[1, 0, 1])), Some(BigInt::from(-4)));
        // x^5 - x - 1: 2869.
        assert_eq!(discriminant(&ints(&[-1, -1, 0, 0, 0, 1])), Some(BigInt::from(2869)));
    }

    #[test]
    fn sextic_groups() {
        // Each polynomial's group is known from its construction.
        let cases: &[(&[i64], &str, usize)] = &[
            // Φ7 and Φ9: cyclic.
            (&[1, 1, 1, 1, 1, 1, 1], "C6", 6),
            (&[1, 0, 0, 1, 0, 0, 1], "C6", 6),
            // Minimal polynomial of 2cos(2π/13).
            (&[-1, 3, 6, -4, -5, 1, 1], "C6", 6),
            // α + ω, α^3 = 2: the regular action of S3.
            (&[9, 9, 0, 3, 6, 3, 1], "S3(6)", 6),
            // x^6 - 2: the affine group of Z/6 with ±1 multipliers.
            (&[-2, 0, 0, 0, 0, 0, 1], "D6(6)", 12),
            // Products a_i a_j over the roots of x^4 + 8x + 12 (group A4).
            (&[1728, 0, -144, -64, -12, 0, 1], "A4(6)", 12),
            // x^6 + 54 x^3 + 2916: roots 3 * 2^(1/3) * ζ9 ...
            (&[2916, 0, 0, 54, 0, 0, 1], "F18(6)=C3xS3", 18),
            // g(x^2) for the cyclic cubic g = x^3 - 3x + 1.
            (&[1, 0, -3, 0, 0, 0, 1], "2A4(6)=C2^3:C3", 24),
            // Σ a_i b_π(i) over the roots of x^3 - x - 1 and x^3 - 2: S3 x S3.
            (&[108, 0, 0, -54, 0, 0, 1], "S3xS3", 36),
            // Sextic resolvent of x^5 - x - 1.
            (&[102400, -3091456, 44800, 8960, 880, 40, 1], "S5(6)=PGL(2,5)", 120),
            (&[-1, -1, 0, 0, 0, 0, 1], "S6", 720),
        ];
        for (c, label, order) in cases {
            assert_eq!(group_of(c), Some((*label, *order)), "{c:?}");
        }
        // S4 acting on six points in two ways.
        for c in [&[-1, 0, 4, 0, 0, 0, 1][..], &[384, 224, 160, 56, 24, 6, 1]] {
            assert!(group_of(c).is_some_and(|(label, order)| order == 24 && label.starts_with("S4")), "{c:?}: {:?}", group_of(c));
        }
    }

    #[test]
    fn parity_matches_discriminant() {
        for c in [&[9, 9, 0, 3, 6, 3, 1][..], &[-2, 0, 0, 0, 0, 0, 1], &[-1, -1, 0, 0, 0, 0, 1], &[1, 1, 1, 1, 1, 1, 1], &[1728, 0, -144, -64, -12, 0, 1]] {
            let a = analyse(&poly(c)).expect("group");
            let disc = discriminant(&a.monic).expect("disc");
            assert_eq!(a.group.is_even(), super::super::galois::is_square(&disc), "{c:?} {}", a.group.label);
        }
    }

    #[test]
    fn septic_groups() {
        // x^7 - 2: the affine group F42; Trinks' x^7 - 7x + 3: PSL(3,2); x^7 - x - 1: S7.
        assert_eq!(group_of(&[-2, 0, 0, 0, 0, 0, 0, 1]), Some(("F42", 42)));
        assert_eq!(group_of(&[3, -7, 0, 0, 0, 0, 0, 1]), Some(("PSL(3,2)", 168)));
        assert_eq!(group_of(&[-1, -1, 0, 0, 0, 0, 0, 1]), Some(("S7", 5040)));
    }

    #[test]
    fn cyclic_septic_by_gauss_periods() {
        // Gaussian periods of length 4 of Q(ζ_29): the cyclic septic x^7 + x^6 - 12x^5 - 7x^4 + 28x^3 + 14x^2 - 9x + 1.
        let a = analyse(&poly(&[1, -9, 14, 28, -7, -12, 1, 1])).expect("group");
        assert_eq!(a.group.label, "C7");
        let cyclo = cyclotomic(&a).expect("periods");
        assert_eq!(cyclo.conductor, 29);
        assert_eq!(cyclo.roots.len(), 7);
    }
}
