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
//! **Gauss periods.** A cyclic group of degree `n` with all roots real gives
//! roots in the degree `n` subfield of some `Q(ζ_m)`, `m` built from the
//! primes of the discriminant (prime powers `p`, `p²` for 3 and 7, `4`, `8`,
//! `16`, composite `m` allowed). The real characters `χ` of `(Z/m)*` of order
//! `n` (found by assigning images to a generating set and closing under
//! multiplication) give the periods `η_k = Σ_{χ(u) = k} ζ_m^u`; the roots are
//! `Σ_t c_t η_{e(k+t)}` with rational `c`, found as the solution of a
//! circulant system and checked to 2^-250.
//!
//! Parallelism: the coset-representative search over the 5040 permutations and the resolvent values
//! of the representatives are computed concurrently (in index order).

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::galois::BITS;
use super::galois::Fx;
use super::galois::one;
use super::galois::refined_roots;
use super::galois::round_integer;
use super::galois::taylor_shift;
use super::transitive::TransitiveGroup;
use super::transitive::transitive_groups;

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
    let derivative: Vec<BigInt> = g
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * BigInt::from(i))
        .collect();
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
    let sign = if (n * (n - 1) / 2) % 2 == 0 {
        BigInt::one()
    } else {
        -BigInt::one()
    };
    Some(sign * resultant / lead)
}

fn bareiss(mut m: Vec<Vec<BigInt>>) -> Option<BigInt> {
    let n = m.len();
    let mut sign = BigInt::one();
    let mut previous = BigInt::one();
    for k in 0..n.saturating_sub(1) {
        if m.get(k)?.get(k)?.is_zero() {
            let swap = (k + 1..n).find(|&r| {
                m.get(r)
                    .and_then(|row| row.get(k))
                    .is_some_and(|v| !v.is_zero())
            })?;
            m.swap(k, swap);
            sign = -sign;
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let value = (m.get(i)?.get(j)? * m.get(k)?.get(k)?
                    - m.get(i)?.get(k)? * m.get(k)?.get(j)?)
                    / &previous;
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
    roots
        .iter()
        .all(|r| r.re.bits() < limit && r.im.bits() < limit)
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

/// Number of permutations / cosets from which the resolvent loops run in parallel.
const RESOLVENT_PAR_GRAIN: usize = 64;

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
    let magnitude_limit = u64::try_from(BITS).unwrap_or(0) + 300;
    for group in groups {
        let elements = group.elements();
        // Coset representatives: sigma minimal in its coset sigma.H.
        // (one membership test per permutation, in parallel, kept in order)
        let minimal = crate::kernels::par::map_slice(&permutations, RESOLVENT_PAR_GRAIN, |sigma| {
            elements.iter().all(|h| {
                let image: Vec<u8> = h
                    .iter()
                    .map(|&i| sigma.get(usize::from(i)).copied().unwrap_or(0))
                    .collect();
                image.as_slice() >= sigma.as_slice()
            })
        });
        let representatives: Vec<&Vec<u8>> = permutations
            .iter()
            .zip(&minimal)
            .filter(|&(_, &m)| m)
            .map(|(sigma, _)| sigma)
            .collect();
        let mut ambiguous = false;
        let mut conclusive = false;
        let choices: Vec<(usize, Vec<usize>)> = (0..powers_by_transform.len())
            .flat_map(|t| exponent_choices(n).into_iter().map(move |e| (t, e)))
            .collect();
        for (transform, e) in choices {
            let powers = powers_by_transform.get(transform)?;
            // the coset resolvent values are independent fixed-point sums
            let values: Vec<Fx> =
                crate::kernels::par::map_slice(&representatives, RESOLVENT_PAR_GRAIN, |sigma| {
                    orbit_sum(&elements, sigma, powers, &e)
                })
                .into_iter()
                .collect::<Option<Vec<Fx>>>()?;
            if values
                .iter()
                .any(|v| v.re.bits() > magnitude_limit || v.im.bits() > magnitude_limit)
            {
                // Out of the working precision for this choice: try another.
                ambiguous = true;
                continue;
            }
            let hits: Vec<usize> = values
                .iter()
                .enumerate()
                .filter(|(_, v)| round_integer(v).is_some())
                .map(|(i, _)| i)
                .collect();
            if hits.is_empty() {
                conclusive = true;
                break;
            }
            // Several cosets are integral when the normaliser is larger
            // than the group; take one whose value no other coset shares.
            let unique = hits.iter().copied().find(|&hit| {
                values.get(hit).is_some_and(|reference| {
                    values
                        .iter()
                        .enumerate()
                        .all(|(i, v)| i == hit || !v.close(reference, 300))
                })
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
        if ambiguous && !conclusive {
            return None;
        }
    }
    None
}

// ---------------------------------------------------------------------
// Linear algebra and rational recognition on fixed-point numbers.

/// Solves `A X = B` (complex fixed point) by Gaussian elimination.
pub(crate) fn solve_system(
    mut a: Vec<Vec<Fx>>,
    mut b: Vec<Vec<Fx>>,
) -> Option<Vec<Vec<Fx>>> {
    let n = a.len();
    // A pivot below 2^-300 of the largest entry is a rounding residue of a
    // singular system.
    let largest = a
        .iter()
        .flatten()
        .map(|v| v.re.abs() + v.im.abs())
        .max()
        .unwrap_or_else(BigInt::zero);
    let negligible = largest >> 300_u32;
    for k in 0..n {
        let pivot = (k..n).max_by_key(|&r| {
            a.get(r)
                .and_then(|row| row.get(k))
                .map_or_else(BigInt::zero, |v| v.re.abs() + v.im.abs())
        })?;
        a.swap(k, pivot);
        b.swap(k, pivot);
        let pivot_value = a.get(k)?.get(k)?.clone();
        if pivot_value.re.abs() + pivot_value.im.abs() <= negligible {
            return None;
        }
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
pub(crate) fn rationalize(z: &Fx) -> Option<BigRational> {
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
        .map(|block| {
            block
                .iter()
                .map(|&p| {
                    analysis
                        .ordering
                        .get(p)
                        .and_then(|&i| analysis.roots.get(i))
                })
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<_>>()?;
    let mut best: Option<(u64, Vec<Fx>, Vec<BigInt>)> = None;
    for (alpha, beta) in [
        (1_i64, 0_i64),
        (0, 1),
        (1, 1),
        (1, 2),
        (2, 1),
        (1, -1),
        (1, 3),
        (3, 1),
        (1, 4),
    ] {
        let theta: Vec<Fx> = blocks
            .iter()
            .map(|block| {
                block.iter().fold(Fx::zero(), |acc, r| {
                    let term = r
                        .mul(&Fx::integer(&BigInt::from(alpha)))
                        .add(&r.mul(r).mul(&Fx::integer(&BigInt::from(beta))));
                    acc.add(&term)
                })
            })
            .collect();
        let distinct = theta
            .iter()
            .enumerate()
            .all(|(i, a)| theta.iter().skip(i + 1).all(|c| !a.close(c, 300)));
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
        let Some(integers) = poly
            .iter()
            .map(round_integer)
            .collect::<Option<Vec<BigInt>>>()
        else {
            continue;
        };
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
    Some(Tower {
        block_size: b,
        theta: integers,
        block,
    })
}

// ---------------------------------------------------------------------
// Gauss periods.

/// Roots of a cyclic polynomial in the degree `n` subfield of `Q(ζ_p)`:
/// `x = shift + z / scale` with `z = Σ_j coefficient_j 2 cos(2π j / p)`.
#[derive(Clone, Debug)]
pub struct Cyclotomic {
    /// The modulus `m` of the cyclotomic field `Q(ζ_m)` containing the roots.
    pub conductor: u64,
    /// For each root: `(position, [(j, coefficient)])`: in the real form the
    /// coefficient of `2cos(2πj/m)`, in the complex form (`complex`) the
    /// coefficient of `exp(2πi j/m)`.
    pub roots: Vec<(usize, Vec<(u64, BigRational)>)>,
    /// Whether the roots are given as combinations of `exp(2πi j/m)` (the
    /// field need not be real).
    pub complex: bool,
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

fn gcd_u64(
    a: u64,
    b: u64,
) -> u64 {
    if b == 0 {
        a
    } else {
        gcd_u64(b, a % b)
    }
}

/// The distinct prime factors of `n` below `limit` (trial division).
fn small_prime_factors(
    n: &BigInt,
    limit: u64,
) -> Vec<u64> {
    let mut rest = n.abs();
    let mut out = Vec::new();
    let mut p = 2_u64;
    while p < limit && rest > BigInt::one() {
        if (&rest % BigInt::from(p)).is_zero() {
            out.push(p);
            while (&rest % BigInt::from(p)).is_zero() {
                rest /= BigInt::from(p);
            }
        }
        p += if p == 2 { 1 } else { 2 };
    }
    out
}

/// A generating set of the unit group `(Z/m)*`, found greedily.
fn unit_generators(
    units: &[u64],
    m: u64,
) -> Vec<u64> {
    let mut generators: Vec<u64> = Vec::new();
    let mut generated: std::collections::HashSet<u64> = std::iter::once(1).collect();
    for &u in units {
        if generated.contains(&u) {
            continue;
        }
        generators.push(u);
        let mut frontier: Vec<u64> = generated.iter().copied().collect();
        while let Some(x) = frontier.pop() {
            for &g in &generators {
                let y = x * g % m;
                if generated.insert(y) {
                    frontier.push(y);
                }
            }
        }
    }
    generators
}

/// The kernels of the surjective characters `(Z/m)* -> Z/n` with
/// `χ(-1) = 0` (real fields), each with the character: `unit -> χ(unit)`.
fn real_characters(
    m: u64,
    n: usize,
    real: bool,
) -> Vec<std::collections::HashMap<u64, usize>> {
    let units: Vec<u64> = (1..m).filter(|&u| gcd_u64(u, m) == 1).collect();
    let generators = unit_generators(&units, m);
    if generators.len() > 5 {
        return Vec::new();
    }
    let total = n.pow(u32::try_from(generators.len()).unwrap_or(0));
    let mut seen: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for code in 0..total {
        let images: Vec<usize> = (0..generators.len())
            .scan(code, |c, _| {
                let digit = *c % n;
                *c /= n;
                Some(digit)
            })
            .collect();
        let mut chi: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
        chi.insert(1, 0);
        let mut frontier = vec![1_u64];
        let mut consistent = true;
        'bfs: while let Some(x) = frontier.pop() {
            let value = chi.get(&x).copied().unwrap_or(0);
            for (g, a) in generators.iter().zip(&images) {
                let y = x * g % m;
                let expected = (value + a) % n;
                match chi.get(&y) {
                    | Some(&v) if v != expected => {
                        consistent = false;
                        break 'bfs;
                    },
                    | Some(_) => {},
                    | None => {
                        chi.insert(y, expected);
                        frontier.push(y);
                    },
                }
            }
        }
        if !consistent || (real && chi.get(&(m - 1)).copied() != Some(0)) {
            continue;
        }
        let image: std::collections::HashSet<usize> = chi.values().copied().collect();
        if image.len() != n {
            continue;
        }
        let mut kernel: Vec<u64> = chi
            .iter()
            .filter(|(_, v)| **v == 0)
            .map(|(u, _)| *u)
            .collect();
        kernel.sort_unstable();
        if seen.insert(kernel) {
            out.push(chi);
        }
    }
    out
}

/// Candidate moduli: products of prime powers (`p`, `p^2` for 3 and 7, `4`,
/// `8`) over the primes dividing the discriminant, in increasing order.
fn candidate_moduli(disc: &BigInt) -> Vec<u64> {
    let primes = small_prime_factors(disc, 20_000);
    let mut options: Vec<Vec<u64>> = Vec::new();
    for &p in primes.iter().take(6) {
        let mut powers = vec![p];
        if p == 2 {
            powers = vec![4, 8, 16];
        } else if p == 3 || p == 7 {
            powers.push(p * p);
        }
        options.push(powers);
    }
    let mut moduli: Vec<u64> = vec![1];
    for powers in &options {
        let mut next = moduli.clone();
        for &m in &moduli {
            for &q in powers {
                let product = m.saturating_mul(q);
                if product <= 4000 {
                    next.push(product);
                }
            }
        }
        // Prime powers of one prime exclude each other.
        moduli = next;
    }
    moduli.retain(|&m| m > 2 && m % 4 != 2);
    moduli.sort_unstable();
    moduli.dedup();
    moduli
}

/// Gauss-period (cyclotomic) form of the roots of a cyclic polynomial with
/// all roots real.
///
/// For each modulus `m` built from the primes of the discriminant, the real
/// characters `χ` of `(Z/m)*` of order `n` give the subfield of `Q(ζ_m)` fixed
/// by `ker χ` with the periods `η_k = Σ_{χ(u) = k} ζ_m^u`; the roots are
/// `Σ_t c_t η_{e(k+t)}` with rational `c`, found as the solution of a
/// circulant system (`e` a unit mod `n` chosen by the labelling).
#[must_use]
pub fn cyclotomic(analysis: &Analysis) -> Option<Cyclotomic> {
    cyclotomic_search(analysis, true)
}

/// Like [`cyclotomic`] for fields that need not be real: every root as a
/// rational combination of `exp(2πi j/m)` (the periods `η_k` are then complex).
#[must_use]
pub fn cyclotomic_complex(analysis: &Analysis) -> Option<Cyclotomic> {
    cyclotomic_search(analysis, false)
}

fn cyclotomic_search(
    analysis: &Analysis,
    real: bool,
) -> Option<Cyclotomic> {
    let n = analysis.degree();
    if !analysis.group.is_abelian() || n < 3 {
        return None;
    }
    let real_tolerance = BigInt::one() << (BITS - 100);
    if real && analysis.roots.iter().any(|r| r.im.abs() > real_tolerance) {
        return None;
    }
    let disc = discriminant(&analysis.monic)?;
    let pi = pi_fixed(BITS + 64);
    for m in candidate_moduli(&disc) {
        let characters = real_characters(m, n, real);
        if characters.is_empty() {
            continue;
        }
        let cos_table: Vec<BigInt> = (0..m).map(|j| cos_fixed(j, m, &pi, 64)).collect();
        // sin(2πu/m) = cos(2π (4u - m) / 4m)
        let sin_table: Vec<BigInt> = (0..m)
            .map(|u| cos_fixed((4 * u + 3 * m) % (4 * m), 4 * m, &pi, 64))
            .collect();
        for chi in &characters {
            let coset = |k: usize| -> Vec<u64> {
                let mut elements: Vec<u64> = chi
                    .iter()
                    .filter(|(_, v)| **v == k)
                    .map(|(u, _)| *u)
                    .collect();
                elements.sort_unstable();
                elements
            };
            for e in (1..n).filter(|e| gcd_u64(*e as u64, n as u64) == 1) {
                let eta: Vec<Fx> = (0..n)
                    .map(|k| {
                        let mut sum = Fx::zero();
                        for j in coset(k) {
                            let index = usize::try_from(j).unwrap_or(0);
                            let re = cos_table.get(index).cloned().unwrap_or_default();
                            let im = if real {
                                BigInt::zero()
                            } else {
                                sin_table.get(index).cloned().unwrap_or_default()
                            };
                            sum = sum.add(&Fx { re, im });
                        }
                        sum
                    })
                    .collect();
                let matrix: Vec<Vec<Fx>> = (0..n)
                    .map(|k| {
                        (0..n)
                            .map(|t| eta.get(e * (k + t) % n).cloned().unwrap_or_else(Fx::zero))
                            .collect()
                    })
                    .collect();
                let rhs: Vec<Vec<Fx>> = (0..n)
                    .map(|k| {
                        analysis
                            .ordering
                            .get(k)
                            .and_then(|&i| analysis.roots.get(i))
                            .map(|r| {
                                vec![if real {
                                    Fx::real(r.re.clone())
                                } else {
                                    r.clone()
                                }]
                            })
                    })
                    .collect::<Option<_>>()?;
                let Some(solution) = solve_system(matrix.clone(), rhs.clone()) else {
                    continue;
                };
                let Some(c) = solution
                    .iter()
                    .map(|row| row.first().and_then(rationalize))
                    .collect::<Option<Vec<BigRational>>>()
                else {
                    continue;
                };
                let ok = (0..n).all(|k| {
                    let mut value = Fx::zero();
                    for (t, q) in c.iter().enumerate() {
                        if let Some(entry) = matrix.get(k).and_then(|row| row.get(t)) {
                            value = value.add(&from_rational(q).mul(entry));
                        }
                    }
                    rhs.get(k)
                        .and_then(|r| r.first())
                        .is_some_and(|target| value.close(target, 250))
                });
                if !ok {
                    continue;
                }
                let mut roots = Vec::new();
                for k in 0..n {
                    let mut terms: Vec<(u64, BigRational)> = Vec::new();
                    for (t, q) in c.iter().enumerate() {
                        for j in coset(e * (k + t) % n) {
                            if real && j > m - j {
                                continue;
                            }
                            match terms.iter_mut().find(|(jj, _)| *jj == j) {
                                | Some((_, slot)) => *slot += q,
                                | None => terms.push((j, q.clone())),
                            }
                        }
                    }
                    terms.retain(|(_, q)| !q.is_zero());
                    terms.sort_by_key(|(j, _)| *j);
                    roots.push((k, terms));
                }
                return Some(Cyclotomic {
                    conductor: m,
                    roots,
                    complex: !real,
                });
            }
        }
    }
    None
}

// ---------------------------------------------------------------------
// Lagrange resolvents of the metacyclic septic groups D7, F21, F42.

/// Data for the real roots of a `D7`, `F21` or `F42` septic in radicals.
///
/// The Lagrange resolvents are `t_j = Σ_k ζ_7^{jk} z_k` and `R_j = t_j^7`,
/// whose conjugates `R_1 … R_6` are the roots of an integer sextic `Φ`.
#[derive(Clone, Debug)]
pub enum Lagrange {
    /// Exactly one real root: with the roots labelled so that it is `z_0` and
    /// conjugation is `k -> -k`, every `t_j` and `R_j` is real and
    /// `z_0 = (t_0 + Σ_{j=1}^{6} R_j^{1/7}) / 7` (real seventh roots). `Φ`
    /// is a totally real cyclic sextic.
    OneReal {
        /// `Φ(w) = Π (w - R_j)`, ascending integer coefficients.
        phi: Vec<BigInt>,
        /// The real numbers `R_1 … R_6`.
        resolvents: Vec<f64>,
        /// `t_0 = Σ z_k`.
        trace: BigInt,
    },
    /// Exactly one real root and a resolvent field that is small in the
    /// invariants `u_j = t_j t_{-j}`: the `u_j` (`j = 1, 2, 3`) are the roots of
    /// the integer cubic `u_poly`, `R_j + R_{-j} = s(u_j)` with `s` a rational
    /// polynomial of degree at most two, and `R_j R_{-j} = u_j^7`. So `R_j,
    /// R_{-j}` are the roots of `w² - s w + u_j^7`, and `z_0` is as in
    /// [`Lagrange::OneReal`]. All numbers involved stay small even where `Φ`
    /// has coefficients and roots of wildly different sizes.
    Paired {
        /// The monic integer cubic of the `u_j`.
        u_poly: Vec<BigInt>,
        /// `s` as a polynomial in `u`.
        s_poly: Vec<BigRational>,
        /// For `j = 1, 2, 3`: `(u_j, R_j, R_{7-j})` as floating-point values.
        terms: Vec<(f64, f64, f64)>,
        /// `t_0 = Σ z_k`.
        trace: BigInt,
    },
    /// All roots real: `R_{-j} = conj R_j`. With `s_j = R_j + R_{-j}` (the
    /// roots of the cubic `s_poly`) and `n_j = |R_j|² = c_0 + c_1 s_j +
    /// c_2 s_j²`, writing `R_j = ρ_j e^{iφ_j}` (`j` with positive imaginary
    /// part, `ρ_j² = n_j`, `cos φ_j = s_j / (2 ρ_j)`) gives
    /// `z_k = (t_0 + 2 Σ_j ρ_j^{1/7} cos((φ_j + 2π m_j)/7 - 2π j k/7)) / 7`.
    AllReal {
        /// The monic integer cubic of the `s_j`.
        s_poly: Vec<BigInt>,
        /// `n` as a polynomial in `s`.
        n_poly: Vec<BigRational>,
        /// For each of the three conjugate pairs: `(j, s_j, m_j)`.
        terms: Vec<(usize, f64, i64)>,
        /// `t_0 = Σ z_k`.
        trace: BigInt,
    },
}

fn fx_pow7(z: &Fx) -> Fx {
    let z2 = z.mul(z);
    let z4 = z2.mul(&z2);
    z4.mul(&z2).mul(z)
}

fn fx_to_f64(v: &BigInt) -> Option<f64> {
    BigRational::new(v.clone(), one()).to_f64()
}

/// The `Paired` form, when the three invariants `u_j = t_j t_{-j}` are distinct.
fn paired_form(
    t: &[Fx],
    r: &[Fx],
    trace: &BigInt,
    integer_poly: &dyn Fn(&[Fx]) -> Option<Vec<BigInt>>,
) -> Option<Lagrange> {
    let u: Vec<Fx> = (1..=3)
        .map(|j| t.get(j).and_then(|a| t.get(7 - j).map(|b| a.mul(b))))
        .collect::<Option<_>>()?;
    let gap = BigInt::one() << (BITS - 100);
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let d = u.get(a)?.sub(u.get(b)?);
        if d.re.abs() < gap {
            return None;
        }
    }
    let u_poly = integer_poly(&u)?;
    let vandermonde: Vec<Vec<Fx>> = u
        .iter()
        .map(|v| vec![Fx::integer(&BigInt::one()), v.clone(), v.mul(v)])
        .collect();
    let rhs: Vec<Vec<Fx>> = (1..=3)
        .map(|j| Some(vec![r.get(j)?.add(r.get(7 - j)?)]))
        .collect::<Option<_>>()?;
    let solution = solve_system(vandermonde, rhs)?;
    let s_poly: Vec<BigRational> = solution
        .iter()
        .map(|row| row.first().and_then(rationalize))
        .collect::<Option<_>>()?;
    let mut terms = Vec::new();
    for (index, j) in (1..=3_usize).enumerate() {
        terms.push((
            fx_to_f64(&u.get(index)?.re)?,
            fx_to_f64(&r.get(j)?.re)?,
            fx_to_f64(&r.get(7 - j)?.re)?,
        ));
    }
    Some(Lagrange::Paired {
        u_poly,
        s_poly,
        terms,
        trace: trace.clone(),
    })
}

/// The Lagrange-resolvent data of a septic with group `D7`, `F21` or `F42`.
#[must_use]
pub fn lagrange(analysis: &Analysis) -> Option<Lagrange> {
    if !matches!(analysis.group.label, "D7" | "F21" | "F42") {
        return None;
    }
    let pi = pi_fixed(BITS + 64);
    // ζ^m = cos(2πm/7) + i sin(2πm/7); sin(x) = cos(x - π/2).
    let zeta: Vec<Fx> = (0..7_u64)
        .map(|m| {
            Fx {
                re: cos_fixed(m, 7, &pi, 64),
                im: cos_fixed((4 * m + 21) % 28, 28, &pi, 64),
            }
        })
        .collect();
    let real_tolerance = BigInt::one() << (BITS - 100);
    let real_positions: Vec<usize> = (0..7)
        .filter(|&p| {
            analysis
                .ordering
                .get(p)
                .and_then(|&i| analysis.roots.get(i))
                .is_some_and(|r| r.im.abs() < real_tolerance)
        })
        .collect();
    let (shift, all_real) = match real_positions.as_slice() {
        | [p] => (*p, false),
        | list if list.len() == 7 => (0, true),
        | _ => return None,
    };
    let z: Vec<Fx> = (0..7)
        .map(|k| {
            analysis
                .ordering
                .get((k + shift) % 7)
                .and_then(|&i| analysis.roots.get(i))
                .cloned()
        })
        .collect::<Option<_>>()?;
    if !all_real {
        for k in 0..7 {
            if !z.get((7 - k) % 7)?.close(&z.get(k)?.conj(), 300) {
                return None;
            }
        }
    }
    let t: Vec<Fx> = (0..7)
        .map(|j| {
            let mut sum = Fx::zero();
            for (k, zk) in z.iter().enumerate() {
                if let Some(w) = zeta.get((j * k) % 7) {
                    sum = sum.add(&w.mul(zk));
                }
            }
            sum
        })
        .collect();
    let trace = round_integer(t.first()?)?;
    let r: Vec<Fx> = t.iter().map(fx_pow7).collect();
    // Elementary symmetric functions of a list of values -> integer polynomial.
    let integer_poly = |values: &[Fx]| -> Option<Vec<BigInt>> {
        let mut poly = vec![Fx::integer(&BigInt::one())];
        for v in values {
            let mut next = vec![Fx::zero(); poly.len() + 1];
            for (i, a) in poly.iter().enumerate() {
                if let Some(slot) = next.get_mut(i + 1) {
                    *slot = slot.add(a);
                }
                if let Some(slot) = next.get_mut(i) {
                    *slot = slot.sub(&a.mul(v));
                }
            }
            poly = next;
        }
        poly.iter().map(round_integer).collect()
    };
    if !all_real {
        if let Some(paired) = paired_form(&t, &r, &trace, &integer_poly) {
            return Some(paired);
        }
        let values: Vec<Fx> = r.iter().skip(1).cloned().collect();
        let phi = integer_poly(&values)?;
        let resolvents: Vec<f64> = values
            .iter()
            .map(|v| fx_to_f64(&v.re))
            .collect::<Option<_>>()?;
        return Some(Lagrange::OneReal {
            phi,
            resolvents,
            trace,
        });
    }
    // All real: pair R_j with R_{7-j}. Pairs with R_j = 0 (the Dickson case has
    // only j = ±1) are dropped; one or three pairs remain.
    let negligible = BigInt::one() << (BITS - 300);
    let pairs: Vec<usize> = (1..=3)
        .filter(|&j| {
            r.get(j)
                .is_some_and(|a| a.re.abs() > negligible || a.im.abs() > negligible)
        })
        .collect();
    if pairs.len() != 1 && pairs.len() != 3 {
        return None;
    }
    let s: Vec<Fx> = pairs
        .iter()
        .map(|&j| r.get(j).and_then(|a| r.get(7 - j).map(|b| a.add(b))))
        .collect::<Option<_>>()?;
    let n: Vec<Fx> = pairs
        .iter()
        .map(|&j| r.get(j).and_then(|a| r.get(7 - j).map(|b| a.mul(b))))
        .collect::<Option<_>>()?;
    let s_poly = integer_poly(&s)?;
    // n_j = c0 + c1 s_j + c2 s_j^2 (constant when there is one pair).
    let size = pairs.len();
    let vandermonde: Vec<Vec<Fx>> = s
        .iter()
        .map(|v| {
            let mut row = vec![Fx::integer(&BigInt::one())];
            for _ in 1..size {
                let next = row.last().map_or_else(Fx::zero, |l| l.mul(v));
                row.push(next);
            }
            row
        })
        .collect();
    let rhs: Vec<Vec<Fx>> = n.iter().map(|v| vec![v.clone()]).collect();
    let solution = solve_system(vandermonde, rhs)?;
    let n_poly: Vec<BigRational> = solution
        .iter()
        .map(|row| row.first().and_then(rationalize))
        .collect::<Option<_>>()?;
    let mut terms = Vec::new();
    for (index, &j) in pairs.iter().enumerate() {
        let rj = r.get(j)?;
        let (representative, value) = if rj.im.is_positive() {
            (j, rj.clone())
        } else {
            (7 - j, r.get(7 - j)?.clone())
        };
        let s_value = fx_to_f64(&s.get(index)?.re)?;
        let rho = (fx_to_f64(&n.get(index)?.re)?).sqrt();
        let phi = (s_value / (2.0 * rho)).clamp(-1.0, 1.0).acos();
        let tj = t.get(representative)?;
        let angle = fx_to_f64(&tj.im)?.atan2(fx_to_f64(&tj.re)?);
        let _ = value;
        let m = ((7.0 * angle - phi) / std::f64::consts::TAU).round() as i64;
        terms.push((representative, s_value, m.rem_euclid(7)));
    }
    Some(Lagrange::AllReal {
        s_poly,
        n_poly,
        terms,
        trace,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(c: &[i64]) -> Vec<BigRational> {
        c.iter()
            .map(|&v| BigRational::from_integer(BigInt::from(v)))
            .collect()
    }

    fn group_of(c: &[i64]) -> Option<(&'static str, usize)> {
        analyse(&poly(c)).map(|a| (a.group.label, a.group.order))
    }

    #[test]
    fn discriminants() {
        // x^3 - 2: -108; x^2 + 1 is not covered (degree >= 2 only): -4.
        let ints = |v: &[i64]| -> Vec<BigInt> { v.iter().map(|&c| BigInt::from(c)).collect() };
        assert_eq!(
            discriminant(&ints(&[-2, 0, 0, 1])),
            Some(BigInt::from(-108))
        );
        assert_eq!(discriminant(&ints(&[1, 0, 1])), Some(BigInt::from(-4)));
        // x^5 - x - 1: 2869.
        assert_eq!(
            discriminant(&ints(&[-1, -1, 0, 0, 0, 1])),
            Some(BigInt::from(2869))
        );
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
            (
                &[102400, -3091456, 44800, 8960, 880, 40, 1],
                "S5(6)=PGL(2,5)",
                120,
            ),
            (&[-1, -1, 0, 0, 0, 0, 1], "S6", 720),
            // Found by systematic search and checked against the cycle types
            // of Frobenius elements (all observed types occur in the group; the
            // types present exclude every smaller candidate):
            // the norm of x^3 + a x + b over Z[sqrt 3] (C2 wr S3) ...
            (&[9, 18, 6, -6, -6, 0, 1], "2S4(6)=C2wrS3", 48),
            // ... over Z[sqrt 2] (S3 wr C2) ...
            (&[-8, -8, -1, 0, -2, 0, 1], "3^2:D4=S3wrC2", 72),
            // The even subgroup of 3^2:D4: the norm over Z[sqrt 2] of
            // f = x^3 - 3 sqrt(2) x - 4, namely f f^sigma =
            // x^6 - 8 x^3 - 18 x^2 + 16.  Its quadratic resolvent
            // K(sqrt(disc f)), K = Q(sqrt 2), is the cyclic quartic field
            // K(sqrt(-3(2 + sqrt 2))): sigma(sqrt(disc f)) = (1 + sqrt 2)
            // sqrt(disc f) with 1 + sqrt 2 of norm -1.  The splitting fields
            // of f and f^sigma therefore meet in that quartic field, giving
            // 6 * 6 / 4 = 36, and the square discriminant 864^4 puts the
            // group in A6: this is 3^2:4, not the full S3 wr C2.  See
            // `sextic_3_2_4` for the tower and Frobenius checks.
            (&[16, 0, -18, -8, 0, 0, 1], "3^2:4", 36),
            // ... the norm of a cubic over Z[sqrt 3] with S3 x S3.
            (&[-3, 0, 0, 6, 0, 0, 1], "S3xS3", 36),
            // The sextic resolvent of x^5 + 20 x - 16 (square discriminant, 3-cycles
            // and 5-cycles among the Frobenius elements: A5): PSL(2,5) on six points.
            (
                &[
                    6_553_600_000_000,
                    -557_056_000_000,
                    7_168_000_000,
                    -71_680_000,
                    352_000,
                    -800,
                    1,
                ],
                "A5(6)=PSL(2,5)",
                60,
            ),
            // x^6 - 24 x - 20: square discriminant, Frobenius elements of types (5,1)
            // and (4,2): only A6 contains both and only even permutations.
            (&[-20, -24, 0, 0, 0, 0, 1], "A6", 360),
        ];
        for (c, label, order) in cases {
            assert_eq!(group_of(c), Some((*label, *order)), "{c:?}");
        }
        // S4 acting on six points in two ways.
        for c in [&[-1, 0, 4, 0, 0, 0, 1][..], &[384, 224, 160, 56, 24, 6, 1]] {
            assert!(
                group_of(c).is_some_and(|(label, order)| order == 24 && label.starts_with("S4")),
                "{c:?}: {:?}",
                group_of(c)
            );
        }
    }

    #[test]
    fn parity_matches_discriminant() {
        for c in [
            &[9, 9, 0, 3, 6, 3, 1][..],
            &[-2, 0, 0, 0, 0, 0, 1],
            &[-1, -1, 0, 0, 0, 0, 1],
            &[1, 1, 1, 1, 1, 1, 1],
            &[1728, 0, -144, -64, -12, 0, 1],
        ] {
            let a = analyse(&poly(c)).expect("group");
            let disc = discriminant(&a.monic).expect("disc");
            assert_eq!(
                a.group.is_even(),
                super::super::galois::is_square(&disc),
                "{c:?} {}",
                a.group.label
            );
        }
    }

    #[test]
    fn septic_groups() {
        // x^7 - 2: the affine group F42; Trinks' x^7 - 7x + 3: PSL(3,2); x^7 - x - 1: S7.
        assert_eq!(group_of(&[-2, 0, 0, 0, 0, 0, 0, 1]), Some(("F42", 42)));
        assert_eq!(
            group_of(&[3, -7, 0, 0, 0, 0, 0, 1]),
            Some(("PSL(3,2)", 168))
        );
        assert_eq!(group_of(&[-1, -1, 0, 0, 0, 0, 0, 1]), Some(("S7", 5040)));
        // The Dickson polynomial D7(x) - 1/4 (times 4): t^2 - 4 = -7 (3/4)^2 puts
        // sqrt(-7) in Q(zeta7), which cuts the multipliers down to the squares: F21.
        assert_eq!(group_of(&[-1, -28, 0, 56, 0, -28, 0, 4]), Some(("F21", 21)));
        // D7(x) - 3: the generic case, F42.
        assert_eq!(group_of(&[-3, -7, 0, 14, 0, -7, 0, 1]), Some(("F42", 42)));
        // Tschirnhaus transforms y = x + x^2 + x^3 of the two Dickson polynomials: the same
        // groups with all six Lagrange resolvents non-zero.
        assert_eq!(
            group_of(&[-6573, -198_548, -241_276, 60_368, 40_656, -2688, -896, 64]),
            Some(("F21", 21))
        );
        assert_eq!(
            group_of(&[-921, -4522, -7105, 462, 539, -42, -14, 1]),
            Some(("F42", 42))
        );
        // The Weber-type septic of the Hilbert class field of Q(sqrt(-71)) (class number 7,
        // discriminant -71^3): the dihedral group D7.
        assert_eq!(group_of(&[1, 2, -1, -1, 1, -1, -1, 1]), Some(("D7", 14)));
    }

    #[test]
    fn composite_conductors() {
        // 2cos(2π/7) + (1 ± sqrt 13)/2: the cyclic sextic of conductor 91 = 7 * 13.
        let a = analyse(&poly(&[-1, -16, 35, 9, -14, -1, 1])).expect("group");
        assert_eq!(a.group.label, "C6");
        let cyclo = cyclotomic(&a).expect("periods");
        assert_eq!(cyclo.conductor, 91);
        assert_eq!(cyclo.roots.len(), 6);
        // The degree-7 field of conductor 29 * 43 (product of two order-7 characters).
        let a = analyse(&poly(&[
            7_143_424, -1_426_432, -332_064, 73_224, 1196, -534, -1, 1,
        ]))
        .expect("group");
        assert_eq!(a.group.label, "C7");
        let cyclo = cyclotomic(&a).expect("periods");
        assert_eq!(cyclo.conductor, 1247);
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

    /// Cycle type (descending) of a permutation given by its images.
    fn cycle_type(p: &[u8]) -> Vec<usize> {
        let mut seen = vec![false; p.len()];
        let mut lengths = Vec::new();
        for start in 0..p.len() {
            if seen[start] {
                continue;
            }
            let (mut j, mut len) = (start, 0_usize);
            while !seen[j] {
                seen[j] = true;
                j = usize::from(p[j]);
                len += 1;
            }
            lengths.push(len);
        }
        lengths.sort_unstable_by(|a, b| b.cmp(a));
        lengths
    }

    fn fp_trim(mut a: Vec<i64>) -> Vec<i64> {
        while a.last().is_some_and(|&c| c == 0) {
            a.pop();
        }
        a
    }

    fn fp_reduce(
        a: &[i64],
        p: i64,
    ) -> Vec<i64> {
        fp_trim(a.iter().map(|&c| c.rem_euclid(p)).collect())
    }

    fn fp_mul(
        a: &[i64],
        b: &[i64],
        p: i64,
    ) -> Vec<i64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut out = vec![0_i64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                if let Some(slot) = out.get_mut(i + j) {
                    *slot = (*slot + x * y).rem_euclid(p);
                }
            }
        }
        fp_trim(out)
    }

    fn fp_sub(
        a: &[i64],
        b: &[i64],
        p: i64,
    ) -> Vec<i64> {
        let mut out = vec![0_i64; a.len().max(b.len())];
        for (i, slot) in out.iter_mut().enumerate() {
            let x = a.get(i).copied().unwrap_or(0);
            let y = b.get(i).copied().unwrap_or(0);
            *slot = (x - y).rem_euclid(p);
        }
        fp_trim(out)
    }

    /// Quotient and remainder in `F_p[x]` (ascending); `b` is non-zero.
    fn fp_div_rem(
        a: &[i64],
        b: &[i64],
        p: i64,
    ) -> Option<(Vec<i64>, Vec<i64>)> {
        let mut remainder = fp_trim(a.to_vec());
        let degree = b.len().checked_sub(1)?;
        let lead = *b.last()?;
        let inverse = (1..p).find(|&x| (x * lead).rem_euclid(p) == 1)?;
        let mut quotient = vec![0_i64; remainder.len().saturating_sub(degree)];
        while remainder.len() > degree {
            let shift = remainder.len() - 1 - degree;
            let factor = (*remainder.last()? * inverse).rem_euclid(p);
            *quotient.get_mut(shift)? = factor;
            for (i, &c) in b.iter().enumerate() {
                let slot = remainder.get_mut(shift + i)?;
                *slot = (*slot - factor * c).rem_euclid(p);
            }
            remainder = fp_trim(remainder);
        }
        Some((fp_trim(quotient), remainder))
    }

    /// The monic gcd in `F_p[x]`.
    fn fp_gcd(
        a: &[i64],
        b: &[i64],
        p: i64,
    ) -> Vec<i64> {
        let (mut a, mut b) = (fp_reduce(a, p), fp_reduce(b, p));
        while !b.is_empty() {
            let (_, remainder) = fp_div_rem(&a, &b, p).unwrap_or_default();
            a = std::mem::replace(&mut b, remainder);
        }
        if let Some(&lead) = a.last() {
            let inverse = (1..p).find(|&x| (x * lead).rem_euclid(p) == 1).unwrap_or(1);
            for c in &mut a {
                *c = (*c * inverse).rem_euclid(p);
            }
        }
        a
    }

    /// `x^exponent mod modulus` in `F_p[x]`.
    fn fp_pow_x(
        mut exponent: u64,
        modulus: &[i64],
        p: i64,
    ) -> Vec<i64> {
        let mut result = vec![1_i64];
        let mut base = vec![0_i64, 1_i64];
        while exponent > 0 {
            if exponent & 1 == 1 {
                let product = fp_mul(&result, &base, p);
                result = fp_div_rem(&product, modulus, p).map_or_else(Vec::new, |(_, r)| r);
            }
            let square = fp_mul(&base, &base, p);
            base = fp_div_rem(&square, modulus, p).map_or_else(Vec::new, |(_, r)| r);
            exponent >>= 1;
        }
        result
    }

    /// Cycle type of the Frobenius at `p` from the degrees of the
    /// irreducible factors of the sextic modulo `p`, by the distinct-degree
    /// decomposition for factors of degree at most three.
    fn fp_cycle_type(
        coeffs: &[i64],
        p: i64,
    ) -> Option<Vec<usize>> {
        let g = fp_reduce(coeffs, p);
        let degree = g.len().checked_sub(1)?;
        let x = vec![0_i64, 1_i64];
        let p_u64 = u64::try_from(p).ok()?;
        let h1 = fp_gcd(&fp_sub(&fp_pow_x(p_u64, &g, p), &x, p), &g, p);
        let h2 = fp_gcd(&fp_sub(&fp_pow_x(p_u64 * p_u64, &g, p), &x, p), &g, p);
        let h3 = fp_gcd(
            &fp_sub(&fp_pow_x(p_u64 * p_u64 * p_u64, &g, p), &x, p),
            &g,
            p,
        );
        let linear = h1.len().checked_sub(1)?;
        let quadratic = (h2.len().checked_sub(1)? - linear) / 2;
        let cubic = (h3.len().checked_sub(1)? - linear) / 3;
        let mut cycle = Vec::new();
        cycle.extend(std::iter::repeat_n(1, linear));
        cycle.extend(std::iter::repeat_n(2, quadratic));
        cycle.extend(std::iter::repeat_n(3, cubic));
        let remainder = degree - linear - 2 * quadratic - 3 * cubic;
        if remainder > 0 {
            cycle.push(remainder);
        }
        cycle.sort_unstable_by(|a, b| b.cmp(a));
        Some(cycle)
    }

    #[test]
    fn sextic_3_2_4() {
        // Norm over Z[sqrt 2] of f = x^3 - 3 sqrt(2) x - 4:
        // f f^sigma = x^6 - 8 x^3 - 18 x^2 + 16.
        let c = &[16, 0, -18, -8, 0, 0, 1];
        let a = analyse(&poly(c)).expect("group");
        assert_eq!((a.group.label, a.group.order), ("3^2:4", 36));
        let d = discriminant(&a.monic).expect("disc");
        assert_eq!(d, BigInt::from(557_256_278_016_i64));
        assert!(super::super::galois::is_square(&d));
        // The unique block system is the pair of cubics of the norm.
        assert_eq!(
            a.group.block_systems(),
            vec![vec![vec![0, 1, 2], vec![3, 4, 5]]]
        );
        let towers = towers(&a);
        assert_eq!(towers.len(), 1);
        let tower = towers.first().expect("tower");
        assert_eq!(tower.block_size, 3);
        // theta = sum of the squares of the roots of a block = +/- 6 sqrt 2.
        assert_eq!(
            tower.theta,
            vec![BigInt::from(-72), BigInt::from(0), BigInt::from(1)]
        );
        // z^3 - (theta/2) z - 4 = f (respectively f^sigma).
        let integer = BigRational::from_integer;
        assert_eq!(
            tower.block,
            vec![
                vec![integer(BigInt::from(-4)), BigRational::zero()],
                vec![
                    BigRational::zero(),
                    BigRational::new(BigInt::from(-1), BigInt::from(2)),
                ],
                vec![BigRational::zero(), BigRational::zero()],
            ]
        );
        // The Frobenius cycle types at these primes are all in the group;
        // among the solvable even sextic groups only 3^2:4 contains them all.
        let observed = [
            (5_i64, vec![4, 2]),
            (7, vec![2, 2, 1, 1]),
            (23, vec![3, 1, 1, 1]),
            (31, vec![3, 3]),
        ];
        let group_types: std::collections::BTreeSet<Vec<usize>> =
            a.group.elements().iter().map(|e| cycle_type(e)).collect();
        for (p, wanted) in &observed {
            assert_eq!(fp_cycle_type(c, *p).as_ref(), Some(wanted), "p = {p}");
            assert!(
                group_types.contains(wanted),
                "Frobenius type {wanted:?} is not in 3^2:4"
            );
        }
        let mut matching = Vec::new();
        for group in transitive_groups(6) {
            if !group.solvable || !group.is_even() {
                continue;
            }
            let types: std::collections::BTreeSet<Vec<usize>> =
                group.elements().iter().map(|e| cycle_type(e)).collect();
            if observed.iter().all(|(_, wanted)| types.contains(wanted)) {
                matching.push(group.label);
            }
        }
        assert_eq!(matching, vec!["3^2:4"]);
    }

    #[test]
    fn sextic_3_2_4_neighbours() {
        // The sibling norm with q = -(4 + sqrt 2) instead of -4: the same
        // block system, but the quadratic resolvent is biquadratic, so the
        // group is S3 x S3 (order 36) and the discriminant is not a square.
        let ints = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
        let sibling = &[14, -12, -18, -8, 0, 0, 1];
        assert_eq!(group_of(sibling), Some(("S3xS3", 36)));
        assert!(!super::super::galois::is_square(
            &discriminant(&ints(sibling)).expect("disc")
        ));
        // c2 - 1 drops the square discriminant and gives the full S3 wr C2.
        let parent = &[16, 0, -19, -8, 0, 0, 1];
        assert_eq!(group_of(parent), Some(("3^2:D4=S3wrC2", 72)));
        assert!(!super::super::galois::is_square(
            &discriminant(&ints(parent)).expect("disc")
        ));
    }
}
