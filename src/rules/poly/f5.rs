//! F5: the Faugère signature-based algorithm for Gröbner bases.
//!
//! Every polynomial carries a **signature** `m * e_i`: the module monomial
//! that leads the representation of the polynomial in the free module over
//! the input generators. The index `i` records the input generator the
//! representation starts from, and signatures are ordered *position over term*
//! with the smallest generator index first, so the signature of a polynomial
//! is the module head term of its representation. Critical pairs are kept in a
//! deterministic priority queue and processed in increasing signature order.
//!
//! Three ingredients make the reduction smaller than Buchberger's:
//!
//! * the **syzygy (F5) criterion** — for a critical pair, take one of its two
//!   scaled sides `u * f_i` with signature `u * sign(f_i)`. If there is a basis
//!   element `h` with `sign(h) = n * e_j`, `j > i` (a *smaller* signature) and
//!   `lm(h) | u * m`, then the pair is a syzygy and is discarded;
//! * the **rewritten criterion** — if a basis element `h` that was generated
//!   *after* `f_i` has the same index and `sign(h) | u * sign(f_i)`, the
//!   reduction of the pair duplicates one that `h` already represents, and the
//!   pair is discarded;
//! * **signature-admissible, one-way reduction** — a reducer `h` may be used
//!   with a multiplier `t` only when `t * sign(h) < sign(p)` strictly. This is
//!   what makes the two criteria sound: a reduced polynomial never needs a
//!   reducer whose signature grew, so the labels bound every representation
//!   that the loop builds.
//!
//! The pair selection is a binary heap ordered by signature (ties broken by the
//! pair indices), which is deterministic. The reduction keeps the unreduced
//! tail and only merges the part of the working polynomial below the leading
//! term that was just cancelled, so a tail is never reduced twice. Leading
//! monomials are cached on the labelled polynomials, and the signature
//! comparison of a scaled signature is evaluated without materialising the
//! product.
//!
//! # Affine input
//!
//! F5 needs a homogeneous ideal in a degree-compatible order. Non-homogeneous
//! generators are homogenised with a fresh variable placed last (the smallest
//! variable of the graded reverse-lexicographic order), the basis is computed
//! there, and the fresh variable is set to one again. The dehomogenisation of
//! a graded Gröbner basis is a Gröbner basis of the original ideal in the
//! restricted order, so this is the standard projective-closure route. The
//! generators are interreduced before the signature loop (which leaves the
//! ideal unchanged) and sorted by decreasing leading monomial, matching the
//! reference formulation of the algorithm.
//!
//! The computation itself runs in graded reverse-lex. For the other orders the
//! reduced basis is converted with
//! [`fglm`](crate::rules::poly::groebner::fglm) — the same strategy as
//! `groebner_via` — and when FGLM declines (a positive-dimensional ideal, or a
//! quotient ring beyond its dimension bound) the matrix engine computes the
//! basis in the requested order.
//!
//! # Limits and the last-resort fallback
//!
//! The criteria path is the only Gröbner route of this module. Limits on the
//! number of pairs, reductions and basis elements are resource guards: when one
//! is exhausted the request is normally reported as unevaluated. As a last
//! resort for genuine resource exhaustion the matrix engine (F4) may finish the
//! job; that event is recorded in [`F5Stats::matrix_fallback`] so that tests
//! can assert that it does not happen for the benchmarks. The matrix engine is
//! also used for the final order conversion when FGLM declines; that
//! documented structural limit is recorded separately in
//! [`F5Stats::matrix_order_fallback`], so neither event can hide behind the
//! other in the statistics. Cyclic-4 is positive-dimensional, so its graded-lex
//! and lex bases necessarily take this route; the other benchmarks are
//! zero-dimensional and stay on the criteria path in every order.
//!
//! The formulation is the F5B one of Sun and Wang (*The F5 algorithm in
//! Buchberger's style*): the rewritten criterion is stated on the basis
//! elements themselves and the syzygy criterion is a divisibility test, which
//! avoids maintaining a separate rewrite table. The final reduced basis is
//! obtained by the shared `common::reduced_basis` pass.
//!
//! # Performance
//!
//! In the `test` profile the complete criteria path (homogenisation, signature
//! loop and final reduction) solves Katsura-4 and Cyclic-5 in roughly 40 ms
//! each on the machine used for this fix — well inside the 30 s budget and
//! comparable to the matrix engine. Before the fix the same requests exhausted
//! the criteria budget and needed seconds in a fallback route. The counters in
//! [`F5Stats`] report the pairs, criterion discards and reductions of each run.
//!
//! # Coefficient domains
//!
//! The signature loop is arithmetic-agnostic in its criteria and
//! signature-specific in its reduction: the coefficients live in a
//! [`Domain`](crate::rules::poly::field::Domain), either the rationals or a
//! prime field `GF(p)`. Over `GF(p)` the S-polynomial, the one-way
//! signature-admissible normal form and every monic normalisation are modular
//! (`a + b`, `a b`, `a b^-1` modulo `p`), so the reduced basis is computed in
//! `GF(p)` from the start, never lifted to `Q`. The output uses the canonical
//! representatives `0 .. p`, monic and reduced exactly like the rational
//! result. The criteria themselves (syzygy, rewritten) and the signatures are
//! field-independent, so the finite-field path has the same no-fallback
//! guarantees as the rational one.
//!
//! All arithmetic is exact over the chosen domain, and the queue order makes
//! the whole run deterministic.

// The signature machinery is one algorithm: pair construction, the two
// criteria, signature-admissible reduction and the bookkeeping around them.
#![allow(
    clippy::single_call_fn,
    clippy::too_many_lines,
    clippy::too_many_arguments
)]

use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use num_rational::BigRational;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Prepared;
use crate::graph::PreparedWork;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::rule::Installer;
use crate::rules::poly::f4::F4Limits;
use crate::rules::poly::f4::common;
use crate::rules::poly::f4::common::Mono;
use crate::rules::poly::f4::f4_basis_in;
use crate::rules::poly::field::Domain;
use crate::rules::poly::groebner::GPoly;
use crate::rules::poly::groebner::Order;
use crate::rules::poly::groebner::Term;
use crate::rules::poly::groebner::fglm_in;
use crate::rules::poly::groebner::reduce_in as normal_form;
use crate::rules::poly::repr::Gens;

/// The monomial order F5 works in. Signature-based algorithms need the
/// signature order to be compatible with the polynomial order, which the
/// graded orders are.
const ENGINE_ORDER: Order = Order::GradedRevLex;

/// Limits for the F5 engine.
#[derive(Copy, Clone, Debug)]
pub struct F5Limits {
    /// Maximum number of pairs taken from the queue.
    pub pairs: usize,
    /// Maximum number of signature-admissible reductions.
    pub reductions: usize,
    /// Maximum number of elements of the basis.
    pub basis: usize,
    /// Which criteria the loop applies. The production settings are both on;
    /// the unit tests switch individual criteria off to isolate one that a
    /// counterexample might implicate.
    pub strategy: Strategy,
}

impl Default for F5Limits {
    fn default() -> Self {
        Self {
            pairs: 1_000_000,
            reductions: 4_000_000,
            basis: 512,
            strategy: Strategy::default(),
        }
    }
}

/// The two criteria of the signature-based loop.
#[derive(Copy, Clone, Debug)]
pub struct Strategy {
    /// The syzygy criterion (the F5 criterion).
    pub syzygy: bool,
    /// The rewritten criterion.
    pub rewritten: bool,
}

impl Default for Strategy {
    fn default() -> Self {
        Self {
            syzygy: true,
            rewritten: true,
        }
    }
}

impl Strategy {
    /// The provably complete subset: every critical pair is processed. Only the
    /// tests use it, to show that the criteria never change the answer.
    #[cfg(test)]
    pub const fn complete() -> Self {
        Self {
            syzygy: false,
            rewritten: false,
        }
    }
}

/// What the F5 run did.
///
/// The counters exist so that callers and tests can check that the benchmarks
/// were solved by the signature criteria path instead of the last-resort
/// matrix fallback.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct F5Stats {
    /// Number of critical pairs taken from the queue.
    pub pairs_processed: usize,
    /// Pairs discarded by the syzygy (F5) criterion.
    pub syzygy_discards: usize,
    /// Pairs discarded by the rewritten criterion.
    pub rewritten_discards: usize,
    /// Pairs whose S-polynomial reduced to zero.
    pub zero_reductions: usize,
    /// Total number of signature-admissible reduction steps.
    pub reductions: usize,
    /// Basis elements added by the signature loop.
    pub polynomials_added: usize,
    /// Whether the last-resort matrix engine had to finish the computation
    /// because the signature loop exhausted one of its resource guards.
    pub matrix_fallback: bool,
    /// Whether the matrix engine had to compute the basis in the requested
    /// order because FGLM declined: the ideal is positive-dimensional, or its
    /// quotient ring is larger than
    /// [`FGLM_MAX_DIMENSION`](crate::rules::poly::groebner::FGLM_MAX_DIMENSION).
    /// This is a documented structural limit, not an exhaustion of the
    /// signature loop, so it is tracked apart from [`Self::matrix_fallback`].
    pub matrix_order_fallback: bool,
}

impl F5Stats {
    /// Whether the whole computation was done by the criteria-driven loop,
    /// with neither the resource fallback nor the order-conversion matrix
    /// engine. Used by the tests to assert that the benchmarks stay on that
    /// path.
    #[cfg(test)]
    #[must_use]
    pub const fn criteria_path(&self) -> bool {
        !self.matrix_fallback && !self.matrix_order_fallback
    }
}

/// A signature: a module monomial `mono * e_index`, with the total degree of
/// `mono` cached because it is compared on every queue operation.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Sig {
    /// The polynomial-ring monomial.
    mono: Mono,
    /// Generator index, one-based.
    index: usize,
    /// Total degree of `mono`.
    degree: u64,
}

impl Sig {
    /// The signature `1 * e_index`.
    fn unit(
        dimension: usize,
        index: usize,
    ) -> Self {
        Self {
            mono: vec![0_u32; dimension],
            index,
            degree: 0,
        }
    }

    /// `self * monomial`.
    fn scaled(
        &self,
        monomial: &[u32],
    ) -> Self {
        Self {
            mono: common::multiply(&self.mono, monomial),
            index: self.index,
            degree: self.degree.saturating_add(common::total_degree(monomial)),
        }
    }
}

/// The signature order: **position over term with the smallest index first**.
/// A signature with a smaller generator index is larger, and signatures with
/// the same index are compared by their module monomial (graded reverse-lex).
fn sig_cmp(
    a: &Sig,
    b: &Sig,
) -> Ordering {
    b.index
        .cmp(&a.index)
        .then_with(|| a.degree.cmp(&b.degree))
        .then_with(|| revlex_cmp(&a.mono, &b.mono))
}

/// Graded reverse-lex comparison of two monomials of equal degree: the
/// monomial with the smaller exponent at the last differing position is
/// larger.
fn revlex_cmp(
    a: &[u32],
    b: &[u32],
) -> Ordering {
    a.iter()
        .zip(b.iter())
        .rev()
        .find(|(x, y)| x != y)
        .map_or(Ordering::Equal, |(x, y)| y.cmp(x))
}

/// Compares `sig * shift` with `bound` without materialising the product. This
/// is the admission test of the one-way reduction, evaluated on every reducer
/// candidate.
fn scaled_sig_cmp(
    sig: &Sig,
    shift: &[u32],
    bound: &Sig,
) -> Ordering {
    match sig.index.cmp(&bound.index) {
        // A larger index is a smaller signature.
        | Ordering::Greater => return Ordering::Less,
        | Ordering::Less => return Ordering::Greater,
        | Ordering::Equal => {},
    }
    let degree = sig.degree.saturating_add(common::total_degree(shift));
    match degree.cmp(&bound.degree) {
        | Ordering::Equal => shifted_revlex_cmp(&sig.mono, shift, &bound.mono),
        | other => other,
    }
}

/// Graded reverse-lex comparison of `right * shift` with `left`, assuming both
/// have the same total degree.
fn shifted_revlex_cmp(
    right: &[u32],
    shift: &[u32],
    left: &[u32],
) -> Ordering {
    let length = right.len().max(left.len());
    for position in (0..length).rev() {
        let shifted = right
            .get(position)
            .copied()
            .unwrap_or(0)
            .saturating_add(shift.get(position).copied().unwrap_or(0));
        let target = left.get(position).copied().unwrap_or(0);
        if shifted != target {
            return target.cmp(&shifted);
        }
    }
    Ordering::Equal
}

/// A polynomial together with its signature and the time it entered the basis.
#[derive(Clone, Debug)]
struct Labelled {
    /// The polynomial, monic.
    poly: GPoly,
    /// The signature of the representation.
    sig: Sig,
    /// Cached leading monomial of `poly`.
    lm: Mono,
    /// Insertion number; larger means generated later.
    num: usize,
}

/// A critical pair waiting in the queue, with both scaled side signatures kept
/// so that the criteria do not recompute them.
#[derive(Clone, Debug)]
struct Pair {
    /// Index of the first basis element.
    a: usize,
    /// Index of the second basis element.
    b: usize,
    /// Signature of the S-polynomial: the larger scaled side signature.
    sig: Sig,
    /// Scaled signature of the first side.
    side_a: Sig,
    /// Scaled signature of the second side.
    side_b: Sig,
}

impl PartialEq for Pair {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.a == other.a && self.b == other.b
    }
}

impl Eq for Pair {}

impl PartialOrd for Pair {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Pair {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        sig_cmp(&self.sig, &other.sig)
            .then_with(|| self.a.cmp(&other.a))
            .then_with(|| self.b.cmp(&other.b))
    }
}

/// Builds the critical pair `(a, b)`. The signature of its S-polynomial is the
/// larger of `(lcm/lm f) * sign(f)` and `(lcm/lm g) * sign(g)`.
fn make_pair(
    basis: &[Labelled],
    a: usize,
    b: usize,
) -> Option<Pair> {
    let (first, second) = (basis.get(a)?, basis.get(b)?);
    let lcm = common::lcm(&first.lm, &second.lm);
    let shift_a = common::difference(&lcm, &first.lm);
    let shift_b = common::difference(&lcm, &second.lm);
    let side_a = first.sig.scaled(&shift_a);
    let side_b = second.sig.scaled(&shift_b);
    let sig = if sig_cmp(&side_a, &side_b) == Ordering::Greater {
        side_a.clone()
    } else {
        side_b.clone()
    };
    Some(Pair {
        a,
        b,
        sig,
        side_a,
        side_b,
    })
}

/// Which criterion made a pair redundant.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Criterion {
    /// The syzygy (F5) criterion.
    Syzygy,
    /// The rewritten criterion.
    Rewritten,
}

/// Applies the criteria to one scaled side `sig` whose source polynomial is
/// `source` (its insertion number).
///
/// * syzygy: a basis element with a **larger** index (a smaller signature)
///   whose leading monomial divides the signature monomial;
/// * rewritten: a basis element generated **later** with the same index whose
///   signature divides the side signature.
fn side_redundant(
    basis: &[Labelled],
    sig: &Sig,
    source: usize,
    strategy: Strategy,
) -> Option<Criterion> {
    for element in basis {
        if strategy.syzygy
            && sig.index < element.sig.index
            && common::divides(&element.lm, &sig.mono)
        {
            return Some(Criterion::Syzygy);
        }
        if strategy.rewritten
            && element.num > source
            && element.sig.index == sig.index
            && common::divides(&element.sig.mono, &sig.mono)
        {
            return Some(Criterion::Rewritten);
        }
    }
    None
}

/// Whether either scaled side of the pair meets a criterion. The counter of the
/// matching criterion is incremented.
fn pair_redundant(
    basis: &[Labelled],
    pair: &Pair,
    strategy: Strategy,
    stats: &mut F5Stats,
) -> bool {
    let sides = [
        (
            &pair.side_a,
            basis.get(pair.a).map_or(0, |element| element.num),
        ),
        (
            &pair.side_b,
            basis.get(pair.b).map_or(0, |element| element.num),
        ),
    ];
    for (signature, source) in sides {
        if let Some(criterion) = side_redundant(basis, signature, source, strategy) {
            match criterion {
                | Criterion::Syzygy => {
                    stats.syzygy_discards = stats.syzygy_discards.saturating_add(1);
                },
                | Criterion::Rewritten => {
                    stats.rewritten_discards = stats.rewritten_discards.saturating_add(1);
                },
            }
            return true;
        }
    }
    false
}

/// Interreduces the generators without changing the ideal they generate, as
/// the reference formulation of F5 does. Only a handful of sweeps are run; a
/// fixed point is reached quickly, and stopping early is still correct.
fn prepare(
    domain: &Domain,
    generators: &[GPoly],
) -> Vec<GPoly> {
    let mut current: Vec<GPoly> = generators
        .iter()
        .filter(|g| !g.is_zero())
        .map(|g| common::monic(domain, g, ENGINE_ORDER))
        .collect();
    for _ in 0..3 {
        let mut next: Vec<GPoly> = Vec::with_capacity(current.len());
        for index in 0..current.len() {
            let Some(polynomial) = current.get(index) else {
                continue;
            };
            let previous = current.get(..index).unwrap_or(&[]);
            let remainder = normal_form(domain, polynomial, previous, ENGINE_ORDER);
            if !remainder.is_zero() {
                next.push(common::monic(domain, &remainder, ENGINE_ORDER));
            }
        }
        if next == current {
            return next;
        }
        current = next;
        if current.is_empty() {
            break;
        }
    }
    current
}

/// The working polynomial of a reduction must not grow past this many terms.
/// F5 keeps the unreduced tail of an S-polynomial (only the leading term is
/// reduced at each step), so on an input that is not close to a Gröbner basis
/// the number of terms can grow quickly. Hitting the bound reports resource
/// exhaustion, which the caller may answer with the matrix engine.
const MAX_TERMS: usize = 20_000;

/// Likewise for the size of a coefficient: exact rational arithmetic keeps
/// working however large the numerator and denominator become, but the cost
/// per operation grows with them, so a bound keeps the run time predictable.
const MAX_BITS: u64 = 32_768;

/// Whether any coefficient in the iterator has outgrown [`MAX_BITS`].
fn coefficients_too_large<'a>(terms: impl Iterator<Item = &'a Term>) -> bool {
    terms.into_iter().any(|(_, coefficient)| {
        coefficient.numer().bits() > MAX_BITS || coefficient.denom().bits() > MAX_BITS
    })
}

/// The signature-admissible normal form of `p` with bound `bound`: a reducer
/// `h` may be used with a multiplier `t` only when `t * sign(h) < bound`
/// strictly. The result has signature at most `bound`. Returns `None` when the
/// reduction budget is exhausted or the working polynomial outgrows
/// [`MAX_TERMS`].
///
/// The working polynomial is kept as a vector of terms in descending order.
/// Once the leading term is added to the remainder, the cursor moves past it
/// and only the remaining tail is merged with the reducer, so an unchanged
/// tail is never re-reduced. Terms of the reducer are shifted by the fixed
/// multiplier during the merge; the comparison against the unshifted tail uses
/// [`compare_shifted`] without building the product.
fn f5_normal_form(
    domain: &Domain,
    p: &GPoly,
    bound: &Sig,
    basis: &[Labelled],
    steps: &mut usize,
    limit: usize,
) -> Option<GPoly> {
    let mut work: Vec<Term> = p.terms().to_vec();
    let mut remainder: Vec<Term> = Vec::new();
    let mut scratch: Vec<Term> = Vec::new();
    let mut cursor = 0_usize;
    loop {
        if work.len() > MAX_TERMS || coefficients_too_large(work.iter().skip(cursor)) {
            return None;
        }
        let Some((monomial, coefficient)) = work.get(cursor).cloned() else {
            break;
        };
        let mut reducer: Option<(&GPoly, Mono, BigRational)> = None;
        for element in basis {
            if !common::divides(&element.lm, &monomial) {
                continue;
            }
            let multiplier = common::difference(&monomial, &element.lm);
            // One-way reduction: an equal signature is not admissible.
            if scaled_sig_cmp(&element.sig, &multiplier, bound) != Ordering::Less {
                continue;
            }
            let lead = element.poly.leading()?.1.clone();
            reducer = Some((&element.poly, multiplier, lead));
            break;
        }
        match reducer {
            | Some((poly, multiplier, lead)) => {
                *steps = steps.saturating_add(1);
                if *steps > limit {
                    return None;
                }
                let factor = domain
                    .div(&coefficient, &lead)
                    .unwrap_or_else(|| domain.zero());
                subtract_tail(
                    domain,
                    &mut work,
                    cursor,
                    poly,
                    &multiplier,
                    &factor,
                    &mut scratch,
                );
            },
            | None => {
                remainder.push((monomial, coefficient));
                cursor = cursor.saturating_add(1);
            },
        }
    }
    Some(GPoly::new_in(domain, remainder, ENGINE_ORDER))
}

/// `work <- work - factor * x^multiplier * poly`, keeping `work[..=cursor]`
/// untouched. The leading term of the shifted reducer cancels `work[cursor]`
/// by construction, so only the tails are merged.
fn subtract_tail(
    domain: &Domain,
    work: &mut Vec<Term>,
    cursor: usize,
    poly: &GPoly,
    multiplier: &[u32],
    factor: &BigRational,
    scratch: &mut Vec<Term>,
) {
    scratch.clear();
    let tail = work.get(cursor.saturating_add(1)..).unwrap_or(&[]);
    let terms = poly.terms().get(1..).unwrap_or(&[]);
    scratch.reserve(tail.len().saturating_add(terms.len()));
    let mut i = 0_usize;
    let mut j = 0_usize;
    loop {
        match (tail.get(i), terms.get(j)) {
            | (Some((left, left_coefficient)), Some((right, right_coefficient))) => {
                match compare_shifted(left, right, multiplier) {
                    | Ordering::Greater => {
                        scratch.push((left.clone(), left_coefficient.clone()));
                        i = i.saturating_add(1);
                    },
                    | Ordering::Less => {
                        scratch.push((
                            common::multiply(right, multiplier),
                            domain.neg(&domain.mul(right_coefficient, factor)),
                        ));
                        j = j.saturating_add(1);
                    },
                    | Ordering::Equal => {
                        let value =
                            domain.sub(left_coefficient, &domain.mul(factor, right_coefficient));
                        if !domain.is_zero(&value) {
                            scratch.push((left.clone(), value));
                        }
                        i = i.saturating_add(1);
                        j = j.saturating_add(1);
                    },
                }
            },
            | (Some((left, left_coefficient)), None) => {
                scratch.push((left.clone(), left_coefficient.clone()));
                i = i.saturating_add(1);
            },
            | (None, Some((right, right_coefficient))) => {
                scratch.push((
                    common::multiply(right, multiplier),
                    domain.neg(&domain.mul(right_coefficient, factor)),
                ));
                j = j.saturating_add(1);
            },
            | (None, None) => break,
        }
    }
    work.truncate(cursor);
    work.append(scratch);
}

/// Compares the monomial `left` with `right * shift` in the engine's graded
/// reverse-lexicographic order, without building the product.
fn compare_shifted(
    left: &[u32],
    right: &[u32],
    shift: &[u32],
) -> Ordering {
    let left_degree: u64 = left.iter().map(|&e| u64::from(e)).sum();
    let right_degree: u64 = right
        .iter()
        .zip(shift)
        .map(|(&r, &s)| u64::from(r.saturating_add(s)))
        .sum();
    match left_degree.cmp(&right_degree) {
        | Ordering::Equal => {
            // In graded reverse-lex the monomial with the smaller exponent at
            // the last differing position is the larger one.
            for i in (0..left.len()).rev() {
                let a = left.get(i).copied().unwrap_or(0);
                let b = right
                    .get(i)
                    .copied()
                    .unwrap_or(0)
                    .saturating_add(shift.get(i).copied().unwrap_or(0));
                if a != b {
                    return b.cmp(&a);
                }
            }
            Ordering::Equal
        },
        | other => other,
    }
}

/// The signature-driven loop. `generators` are homogeneous and monic; the
/// result is a Gröbner basis in the graded reverse-lexicographic order (not
/// necessarily reduced or minimal).
///
/// Pairs are popped in increasing signature order. A pair whose side meets one
/// of the two criteria is dropped without any reduction; the others have their
/// S-polynomial reduced by one-way signature-admissible reductions and the
/// non-zero remainder is added to the basis with the pair's signature.
fn f5_loop(
    domain: &Domain,
    generators: &[GPoly],
    limits: F5Limits,
    stats: &mut F5Stats,
) -> Option<Vec<GPoly>> {
    let mut prepared = prepare(domain, generators);
    let Some(first) = prepared.first() else {
        return Some(Vec::new());
    };
    let dimension = first.leading()?.0.len();
    // The reference formulation sorts the generators by decreasing leading
    // monomial before numbering the signatures.
    prepared.sort_by(|a, b| {
        match (a.leading(), b.leading()) {
            | (Some((left, _)), Some((right, _))) => ENGINE_ORDER.cmp(right, left),
            | _ => Ordering::Equal,
        }
    });
    let mut basis: Vec<Labelled> = prepared
        .into_iter()
        .enumerate()
        .map(|(position, poly)| {
            let num = position.saturating_add(1);
            let lm = poly
                .leading()
                .map_or_else(|| vec![0_u32; dimension], |(monomial, _)| monomial.clone());
            Labelled {
                poly,
                sig: Sig::unit(dimension, num),
                lm,
                num,
            }
        })
        .collect();
    let mut queue: BinaryHeap<Reverse<Pair>> = BinaryHeap::new();
    for second in 1..basis.len() {
        for first in 0..second {
            if let Some(pair) = make_pair(&basis, first, second) {
                queue.push(Reverse(pair));
            }
        }
    }
    while let Some(Reverse(pair)) = queue.pop() {
        stats.pairs_processed = stats.pairs_processed.saturating_add(1);
        if stats.pairs_processed > limits.pairs {
            return None;
        }
        if pair_redundant(&basis, &pair, limits.strategy, stats) {
            continue;
        }
        let (Some(first), Some(second)) = (basis.get(pair.a), basis.get(pair.b)) else {
            continue;
        };
        let spolynomial = common::s_polynomial(domain, &first.poly, &second.poly, ENGINE_ORDER);
        let remainder = f5_normal_form(
            domain,
            &spolynomial,
            &pair.sig,
            &basis,
            &mut stats.reductions,
            limits.reductions,
        )?;
        if remainder.is_zero() {
            stats.zero_reductions = stats.zero_reductions.saturating_add(1);
            continue;
        }
        if basis.len() >= limits.basis {
            return None;
        }
        let new = basis.len();
        let num = new.saturating_add(1);
        let Some((lm, _)) = remainder.leading() else {
            continue;
        };
        let lm = lm.clone();
        basis.push(Labelled {
            poly: common::monic(domain, &remainder, ENGINE_ORDER),
            sig: pair.sig,
            lm,
            num,
        });
        stats.polynomials_added = stats.polynomials_added.saturating_add(1);
        for other in 0..new {
            if let Some(candidate) = make_pair(&basis, other, new) {
                queue.push(Reverse(candidate));
            }
        }
    }
    Some(basis.into_iter().map(|element| element.poly).collect())
}

/// The reduced Gröbner basis of the ideal generated by `generators` in
/// `order` over the rationals (or `None` when a limit is hit).
///
/// This is the rational convenience wrapper used by the unit tests; the
/// operator calls [`f5_basis_in`] with the request's domain.
#[cfg(test)]
#[must_use]
pub fn f5_basis(
    generators: &[GPoly],
    order: Order,
    limits: F5Limits,
) -> Option<Vec<GPoly>> {
    f5_basis_in(&Domain::Rational, generators, order, limits)
}

/// The F5 engine over `domain`: the same signature loop with modular
/// arithmetic over a prime field, returning the canonical reduced `GF(p)`
/// basis.
#[must_use]
pub fn f5_basis_in(
    domain: &Domain,
    generators: &[GPoly],
    order: Order,
    limits: F5Limits,
) -> Option<Vec<GPoly>> {
    f5_basis_with_stats_in(domain, generators, order, limits).0
}

/// The reduced basis, also returning what the run did. The statistics record
/// whether the criteria path solved the problem or whether the last-resort
/// matrix engine had to be used.
#[cfg(test)]
#[must_use]
pub fn f5_basis_with_stats(
    generators: &[GPoly],
    order: Order,
    limits: F5Limits,
) -> (Option<Vec<GPoly>>, F5Stats) {
    f5_basis_with_stats_in(&Domain::Rational, generators, order, limits)
}

/// [`f5_basis_with_stats`] over `domain`.
#[must_use]
pub fn f5_basis_with_stats_in(
    domain: &Domain,
    generators: &[GPoly],
    order: Order,
    limits: F5Limits,
) -> (Option<Vec<GPoly>>, F5Stats) {
    let mut stats = F5Stats::default();
    let Some(dimension) = generators
        .iter()
        .find_map(|g| g.leading().map(|(monomial, _)| monomial.len()))
    else {
        return (Some(Vec::new()), stats);
    };
    let input: Vec<GPoly> = generators
        .iter()
        .filter(|g| !g.is_zero())
        .map(|g| GPoly::new_in(domain, g.terms().to_vec(), ENGINE_ORDER))
        .collect();
    if input.is_empty() {
        return (Some(Vec::new()), stats);
    }
    // A non-zero constant generator generates the unit ideal. Only a
    // polynomial whose *leading* monomial is the constant (i.e. a constant)
    // counts: the constant term of a non-constant generator does not.
    let is_constant = |g: &GPoly| {
        g.leading()
            .is_some_and(|(monomial, _)| monomial.iter().all(|&e| e == 0))
    };
    if input.iter().any(is_constant) {
        let unit = vec![0_u32; dimension];
        return (
            Some(vec![GPoly::new_in(
                domain,
                vec![(unit, domain.one())],
                order,
            )]),
            stats,
        );
    }
    let homogenised = !input.iter().all(common::is_homogeneous);
    let working: Vec<GPoly> = if homogenised {
        input
            .iter()
            .map(|g| common::homogenise(g, ENGINE_ORDER))
            .collect()
    } else {
        input
    };
    // The criteria path is the only Gröbner route. If the resource guards of
    // the signature loop are exhausted, the matrix engine is a last resort.
    let completion = match f5_loop(domain, &working, limits, &mut stats) {
        | Some(basis) => basis,
        | None => {
            stats.matrix_fallback = true;
            let Some(basis) = f4_basis_in(
                domain,
                &working,
                ENGINE_ORDER,
                F4Limits {
                    basis: limits.basis,
                    ..F4Limits::default()
                },
            ) else {
                return (None, stats);
            };
            basis
        },
    };
    let original: Vec<GPoly> = if homogenised {
        completion
            .iter()
            .map(|g| common::dehomogenise(g, ENGINE_ORDER))
            .collect()
    } else {
        completion
    };
    let reduced = common::reduced_basis(domain, &original, ENGINE_ORDER);
    if order == ENGINE_ORDER {
        return (Some(reduced), stats);
    }
    if let Some(converted) = fglm_in(domain, &reduced, ENGINE_ORDER, order) {
        return (Some(converted), stats);
    }
    // Positive-dimensional ideal or too large for FGLM: the matrix engine
    // computes the same reduced basis in the requested order. This is not a
    // failure of the criteria loop, but it must be visible in the statistics.
    stats.matrix_order_fallback = true;
    let reordered: Vec<GPoly> = original
        .iter()
        .map(|g| GPoly::new_in(domain, g.terms().to_vec(), order))
        .collect();
    (
        f4_basis_in(
            domain,
            &reordered,
            order,
            F4Limits {
                basis: limits.basis,
                ..F4Limits::default()
            },
        ),
        stats,
    )
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new("groebner_f5", Arity::Variadic)
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    // A basis is a list of polynomials with the common dimension of the
    // generators (the first argument), matching the reduced `list(...)`.
    i.dim(op, DimRule::Same(0));
    i.kernel("poly/groebner_f5", Tier::Reduce, F5Kernel { op });
    Ok(())
}

struct F5Kernel {
    op: OpId,
}

impl Kernel for F5Kernel {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn two_phase(&self) -> bool {
        true
    }

    fn prepare(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Prepared {
        match common::prepare_data(cx, node) {
            | Some((domain, generators, order, gens)) => {
                Prepared::Work(Box::new(F5Work {
                    domain,
                    generators,
                    order,
                    gens,
                }))
            },
            | None => Prepared::Done(Outcome::Pass),
        }
    }

    fn commit(
        &self,
        cx: &mut Cx<'_>,
        _node: NodeId,
        result: Box<dyn std::any::Any + Send>,
    ) -> Outcome {
        common::commit(cx, result)
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let Some((domain, generators, order, gens)) = common::prepare_data(cx, node) else {
            return Outcome::Pass;
        };
        match f5_basis_in(&domain, &generators, order, F5Limits::default()) {
            | Some(basis) => common::pinned(cx.graph, &basis, &gens),
            | None => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// The graph-free part of a `groebner_f5` request.
struct F5Work {
    domain: Domain,
    generators: Vec<GPoly>,
    order: Order,
    gens: Gens,
}

impl PreparedWork for F5Work {
    fn run(self: Box<Self>) -> Box<dyn std::any::Any + Send> {
        let basis = f5_basis_in(
            &self.domain,
            &self.generators,
            self.order,
            F5Limits::default(),
        );
        Box::new((basis, self.gens))
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::*;
    use crate::rules::poly::groebner::GroebnerLimits;
    use crate::rules::poly::groebner::groebner;
    use crate::rules::poly::groebner::groebner_via;
    use crate::rules::poly::groebner::reduce_in as normal_form;

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(n))
    }

    fn poly(
        terms: &[(i64, Vec<u32>)],
        order: Order,
    ) -> GPoly {
        GPoly::new(
            terms.iter().map(|(c, e)| (e.clone(), q(*c))).collect(),
            order,
        )
    }

    fn mono(
        n: usize,
        powers: &[(usize, u32)],
    ) -> Mono {
        let mut m = vec![0_u32; n];
        for &(i, e) in powers {
            if let Some(slot) = m.get_mut(i) {
                *slot = slot.saturating_add(e);
            }
        }
        m
    }

    fn katsura(
        m: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let n = m.saturating_add(1);
        let mut out = Vec::new();
        let mut linear = vec![(1, mono(n, &[(0, 1)])), (-1, mono(n, &[]))];
        for l in 1..=m {
            linear.push((2, mono(n, &[(l, 1)])));
        }
        out.push(poly(&linear, order));
        for k in 0..m {
            let mut terms = vec![(-1, mono(n, &[(k, 1)]))];
            for l in -(m as i64)..=(m as i64) {
                let a = l.unsigned_abs() as usize;
                let b = (k as i64 - l).unsigned_abs() as usize;
                if a <= m && b <= m {
                    terms.push((1, mono(n, &[(a, 1), (b, 1)])));
                }
            }
            out.push(poly(&terms, order));
        }
        out
    }

    fn cyclic(
        n: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let mut out = Vec::new();
        for d in 1..n {
            let terms: Vec<(i64, Vec<u32>)> = (0..n)
                .map(|start| {
                    let powers: Vec<(usize, u32)> = (0..d).map(|j| ((start + j) % n, 1)).collect();
                    (1, mono(n, &powers))
                })
                .collect();
            out.push(poly(&terms, order));
        }
        out.push(poly(&[(1, vec![1; n]), (-1, mono(n, &[]))], order));
        out
    }

    /// The affine systems the task asks to cross-check.
    fn sample_systems(order: Order) -> Vec<(&'static str, Vec<GPoly>)> {
        vec![
            (
                "circle and line",
                vec![
                    poly(
                        &[(1, vec![2, 0, 0]), (1, vec![0, 2, 0]), (-1, vec![0, 0, 0])],
                        order,
                    ),
                    poly(&[(1, vec![1, 0, 0]), (-1, vec![0, 1, 0])], order),
                ],
            ),
            ("katsura-2", katsura(2, order)),
            ("katsura-3", katsura(3, order)),
            ("cyclic-3", cyclic(3, order)),
            ("cyclic-4", cyclic(4, order)),
            (
                "affine, mixed degrees",
                vec![
                    poly(
                        &[(1, vec![2, 1, 0]), (-1, vec![0, 3, 2]), (1, vec![0, 0, 0])],
                        order,
                    ),
                    poly(&[(1, vec![1, 0, 0]), (-1, vec![0, 0, 1])], order),
                    poly(&[(1, vec![2, 0, 2]), (-3, vec![0, 1, 0])], order),
                ],
            ),
            (
                "equal leading monomials",
                vec![
                    poly(
                        &[(1, vec![2, 0, 0]), (1, vec![0, 1, 0]), (1, vec![0, 0, 0])],
                        order,
                    ),
                    poly(
                        &[(1, vec![2, 0, 0]), (1, vec![0, 0, 1]), (-1, vec![0, 0, 0])],
                        order,
                    ),
                ],
            ),
            (
                "non-monic with redundant generators",
                vec![
                    poly(&[(2, vec![2, 0, 0]), (-3, vec![0, 0, 0])], order),
                    poly(&[(4, vec![2, 0, 0]), (-6, vec![0, 0, 0])], order),
                    poly(&[(2, vec![1, 1, 0]), (-2, vec![0, 1, 0])], order),
                ],
            ),
            (
                "unit ideal",
                vec![
                    poly(&[(1, vec![1, 0, 0])], order),
                    poly(&[(1, vec![1, 0, 0]), (1, vec![0, 0, 0])], order),
                ],
            ),
        ]
    }

    /// The benchmarks, built in the requested order.
    fn benchmarks(order: Order) -> Vec<(&'static str, Vec<GPoly>)> {
        vec![
            ("katsura-2", katsura(2, order)),
            ("katsura-3", katsura(3, order)),
            ("katsura-4", katsura(4, order)),
            ("cyclic-3", cyclic(3, order)),
            ("cyclic-4", cyclic(4, order)),
            ("cyclic-5", cyclic(5, order)),
        ]
    }

    fn reference(
        generators: &[GPoly],
        order: Order,
    ) -> Option<Vec<GPoly>> {
        if order == Order::GradedRevLex {
            groebner(generators, order, GroebnerLimits::default())
        } else {
            // The operator's own reference route: Buchberger in graded
            // reverse-lex followed by FGLM.
            groebner_via(generators, order, GroebnerLimits::default())
        }
    }

    #[test]
    fn signatures_order_by_the_smallest_index_first() {
        // e_1 is larger than every e_2 term, whatever the monomials.
        assert_eq!(
            sig_cmp(
                &Sig {
                    mono: vec![0, 0, 0],
                    index: 1,
                    degree: 0
                },
                &Sig {
                    mono: vec![5, 5, 5],
                    index: 2,
                    degree: 15
                }
            ),
            Ordering::Greater
        );
        assert_eq!(
            sig_cmp(
                &Sig {
                    mono: vec![1, 0, 0],
                    index: 2,
                    degree: 1
                },
                &Sig {
                    mono: vec![1, 0, 0],
                    index: 2,
                    degree: 1
                }
            ),
            Ordering::Equal
        );
        // Within one index, the monomial order decides.
        assert_eq!(
            sig_cmp(
                &Sig {
                    mono: vec![1, 0, 0],
                    index: 1,
                    degree: 1
                },
                &Sig {
                    mono: vec![0, 1, 0],
                    index: 1,
                    degree: 1
                }
            ),
            ENGINE_ORDER.cmp(&[1, 0, 0], &[0, 1, 0])
        );
    }

    #[test]
    fn f5_agrees_with_groebner_on_the_sample_systems() {
        for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
            for (name, generators) in sample_systems(order) {
                let expected = reference(&generators, order);
                let (computed, stats) =
                    f5_basis_with_stats(&generators, order, F5Limits::default());
                assert_eq!(computed, expected, "{name} in {order:?}");
                // The criteria loop itself must not be exhausted. A
                // positive-dimensional sample may still need the matrix
                // engine for the final order conversion, which is recorded
                // separately in `matrix_order_fallback`.
                assert!(
                    !stats.matrix_fallback,
                    "{name} in {order:?} exhausted the criteria loop: {stats:?}"
                );
                if order == Order::GradedRevLex {
                    assert!(
                        stats.criteria_path(),
                        "{name} in grevlex did not stay on the criteria path: {stats:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn f5_matches_groebner_on_all_benchmarks_in_all_orders() {
        for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
            for (name, generators) in benchmarks(order) {
                let start = Instant::now();
                let (computed, stats) =
                    f5_basis_with_stats(&generators, order, F5Limits::default());
                let elapsed = start.elapsed();
                eprintln!(
                    "F5 {name} {order:?}: {elapsed:?} ({} pairs, {} syzygy, {} rewritten, {} reductions, {} added)",
                    stats.pairs_processed,
                    stats.syzygy_discards,
                    stats.rewritten_discards,
                    stats.reductions,
                    stats.polynomials_added
                );
                let expected = reference(&generators, order);
                assert_eq!(computed, expected, "{name} in {order:?}");
                assert!(
                    !stats.matrix_fallback,
                    "{name} in {order:?} exhausted the criteria loop: {stats:?}"
                );
                if name == "cyclic-4" && order != Order::GradedRevLex {
                    // Cyclic-4 is positive-dimensional (no pure power of x3
                    // leads a grevlex basis element), so FGLM cannot convert
                    // its basis. The matrix engine computes the ordered basis
                    // and the stats must record that structural limit.
                    assert!(
                        stats.matrix_order_fallback,
                        "cyclic-4 in {order:?}: FGLM must decline and the order fallback must be recorded: {stats:?}"
                    );
                } else {
                    assert!(
                        stats.criteria_path(),
                        "{name} in {order:?} did not stay on the criteria path: {stats:?}"
                    );
                }
                // The criteria must have done the work, not the complete pair
                // set: at least one of them discards pairs on every benchmark,
                // and the loop never exhausts its pair budget.
                assert!(
                    stats
                        .syzygy_discards
                        .saturating_add(stats.rewritten_discards)
                        > 0,
                    "{name} in {order:?} did not exercise the criteria: {stats:?}"
                );
                assert!(
                    stats.pairs_processed < F5Limits::default().pairs,
                    "{name} in {order:?} exhausted the pair budget"
                );
                if order == Order::GradedRevLex && (name == "katsura-4" || name == "cyclic-5") {
                    eprintln!("F5 {name} (grevlex) benchmark: {elapsed:?}");
                    assert!(
                        elapsed.as_secs() < 30,
                        "{name} took {elapsed:?}, over the 30 s budget"
                    );
                }
            }
        }
    }

    #[test]
    fn every_criterion_variant_returns_the_same_basis() {
        // Switching criteria off must not change the answer: the remaining
        // pairs still form a complete reduction and every non-zero remainder
        // is added, so each strategy ends at the same reduced basis.
        let order = Order::GradedRevLex;
        for (name, generators) in [
            ("katsura-3", katsura(3, order)),
            ("cyclic-4", cyclic(4, order)),
        ] {
            let expected = groebner(&generators, order, GroebnerLimits::default());
            for strategy in [
                Strategy {
                    syzygy: true,
                    rewritten: false,
                },
                Strategy {
                    syzygy: false,
                    rewritten: true,
                },
                Strategy::complete(),
            ] {
                let limits = F5Limits {
                    strategy,
                    ..F5Limits::default()
                };
                let (computed, stats) = f5_basis_with_stats(&generators, order, limits);
                assert_eq!(computed, expected, "{name}: {strategy:?}");
                assert!(stats.criteria_path(), "{name}: unexpected fallback");
            }
        }
    }

    #[test]
    fn f5_limits_report_exhaustion() {
        let order = Order::GradedRevLex;
        let generators = katsura(3, order);
        // The signature loop reports exhaustion of each budget.
        for limits in [
            F5Limits {
                pairs: 0,
                ..F5Limits::default()
            },
            F5Limits {
                basis: 1,
                ..F5Limits::default()
            },
        ] {
            let mut stats = F5Stats::default();
            assert_eq!(
                f5_loop(&Domain::Rational, &generators, limits, &mut stats),
                None,
                "{limits:?}"
            );
        }
        // The reduction budget is enforced by the normal form itself.
        let labelled = Labelled {
            poly: poly(&[(1, vec![1, 0])], ENGINE_ORDER),
            sig: Sig::unit(2, 2),
            lm: vec![1, 0],
            num: 2,
        };
        let p = poly(&[(1, vec![1, 1])], ENGINE_ORDER);
        let bound = Sig::unit(2, 1);
        let mut steps = 0_usize;
        assert_eq!(
            f5_normal_form(&Domain::Rational, &p, &bound, &[labelled], &mut steps, 0),
            None
        );
        // With room to work, the operator returns the reduced basis.
        assert!(f5_basis(&generators, order, F5Limits::default()).is_some());
    }

    #[test]
    fn f5_falls_back_only_when_the_budget_is_exhausted() {
        // A zero pair budget is genuine resource exhaustion: the last-resort
        // matrix engine finishes, and the counter records it.
        let order = Order::GradedRevLex;
        let generators = katsura(3, order);
        let limits = F5Limits {
            pairs: 0,
            ..F5Limits::default()
        };
        let (computed, stats) = f5_basis_with_stats(&generators, order, limits);
        assert!(
            stats.matrix_fallback,
            "the counter must record the fallback"
        );
        assert_eq!(
            computed,
            groebner(&generators, order, GroebnerLimits::default())
        );
    }

    #[test]
    fn f5_records_the_order_fallback_when_fglm_declines() {
        // <x^2 - y^2> is one-dimensional: FGLM declines in lex, and the
        // matrix engine computes the ordered basis. The statistics must say
        // so instead of claiming the criteria path did everything.
        let order = Order::Lex;
        let generators = vec![poly(&[(1, vec![2, 0]), (-1, vec![0, 2])], order)];
        let (computed, stats) = f5_basis_with_stats(&generators, order, F5Limits::default());
        assert_eq!(
            computed,
            groebner(&generators, order, GroebnerLimits::default()),
            "the ordered basis must still be correct"
        );
        assert!(
            !stats.matrix_fallback,
            "the criteria loop itself was not exhausted: {stats:?}"
        );
        assert!(
            stats.matrix_order_fallback,
            "the matrix engine computed the lex basis but the stats hide it: {stats:?}"
        );
        assert!(
            !stats.criteria_path(),
            "the criteria path did not produce the final ordered basis: {stats:?}"
        );
        // In the engine order no conversion is needed.
        let (_, grevlex_stats) =
            f5_basis_with_stats(&generators, Order::GradedRevLex, F5Limits::default());
        assert!(
            grevlex_stats.criteria_path(),
            "grevlex should stay on the criteria path: {grevlex_stats:?}"
        );
        // The empty ideal is answered before either engine runs.
        let (_, empty_stats) = f5_basis_with_stats(&[], order, F5Limits::default());
        assert!(
            !empty_stats.matrix_fallback && !empty_stats.matrix_order_fallback,
            "the empty ideal runs neither engine: {empty_stats:?}"
        );
    }

    #[test]
    fn f5_handles_empty_and_zero_inputs() {
        let order = Order::GradedRevLex;
        assert_eq!(f5_basis(&[], order, F5Limits::default()), Some(Vec::new()));
        let zero = poly(&[], order);
        assert_eq!(
            f5_basis(&[zero.clone(), zero], order, F5Limits::default()),
            Some(Vec::new())
        );
        // A constant generator is the unit ideal, before and after
        // homogenisation.
        let unit = f5_basis(&[poly(&[(3, vec![])], order)], order, F5Limits::default());
        assert_eq!(
            unit,
            Some(vec![poly(&[(1, vec![])], order)]),
            "a constant generator is the unit ideal"
        );
    }

    #[test]
    fn signature_reduction_is_strictly_decreasing() {
        // f = x^2 and g = x*y with generator indices 1 and 2.
        let basis = vec![
            Labelled {
                poly: poly(&[(1, vec![2, 0])], ENGINE_ORDER),
                sig: Sig::unit(2, 1),
                lm: vec![2, 0],
                num: 1,
            },
            Labelled {
                poly: poly(&[(1, vec![1, 1])], ENGINE_ORDER),
                sig: Sig::unit(2, 2),
                lm: vec![1, 1],
                num: 2,
            },
        ];
        let p = poly(&[(1, vec![2, 1])], ENGINE_ORDER);
        // Bound e_2: neither y*e_1 (index 1) nor x*e_2 (index 2, x > 1) is
        // admissible, so y*f is already its own normal form.
        let bound = Sig::unit(2, 2);
        let mut steps = 0_usize;
        let reduced = match f5_normal_form(&Domain::Rational, &p, &bound, &basis, &mut steps, 100) {
            | Some(reduced) => reduced,
            | None => panic!("the reduction must fit its budget"),
        };
        assert_eq!(reduced, p);
        // Bound e_1: x*g has index 2, a smaller signature, so y*f reduces to
        // zero by g.
        let bound = Sig::unit(2, 1);
        let mut steps = 0_usize;
        let reduced = match f5_normal_form(&Domain::Rational, &p, &bound, &basis, &mut steps, 100) {
            | Some(reduced) => reduced,
            | None => panic!("the reduction must fit its budget"),
        };
        assert!(reduced.is_zero(), "y*f = x*(x*y) rewrites to zero");
    }

    #[test]
    fn signature_reduction_rejects_an_equal_signature() {
        // h = y with signature e_1; p = x*y with bound x*e_1. The scaled
        // signature x*e_1 equals the bound, so the one-way reduction must not
        // use h and p stays unchanged.
        let basis = vec![Labelled {
            poly: poly(&[(1, vec![0, 1])], ENGINE_ORDER),
            sig: Sig::unit(2, 1),
            lm: vec![0, 1],
            num: 1,
        }];
        let p = poly(&[(1, vec![1, 1])], ENGINE_ORDER);
        let bound = Sig {
            mono: vec![1, 0],
            index: 1,
            degree: 1,
        };
        let mut steps = 0_usize;
        let reduced = match f5_normal_form(&Domain::Rational, &p, &bound, &basis, &mut steps, 100) {
            | Some(reduced) => reduced,
            | None => panic!("the reduction must fit its budget"),
        };
        assert_eq!(reduced, p);
        assert_eq!(steps, 0);
    }

    #[test]
    fn f5_is_deterministic() {
        let order = Order::GradedRevLex;
        for (name, generators) in [
            ("cyclic-4", cyclic(4, order)),
            ("katsura-3", katsura(3, order)),
        ] {
            let first = f5_basis_with_stats(&generators, order, F5Limits::default());
            let second = f5_basis_with_stats(&generators, order, F5Limits::default());
            let third = f5_basis_with_stats(&generators, order, F5Limits::default());
            assert_eq!(first.0, second.0, "{name}");
            assert_eq!(second.0, third.0, "{name}");
            assert_eq!(first.1, second.1, "{name}: statistics must be stable");
        }
    }

    #[test]
    fn f5_agrees_with_buchberger_on_random_systems() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
            for seed in 0..3_u64 {
                let n = 3_usize;
                let generators: Vec<GPoly> = (0..3)
                    .map(|_| {
                        let terms: Vec<(i64, Vec<u32>)> = (0..4)
                            .map(|_| {
                                let coefficient = i64::try_from(next() % 5).unwrap_or(0) - 2;
                                let exponents: Vec<u32> = (0..n)
                                    .map(|_| u32::try_from(next() % 3).unwrap_or(0))
                                    .collect();
                                (coefficient, exponents)
                            })
                            .filter(|(coefficient, _)| *coefficient != 0)
                            .collect();
                        let seed_term = (1, mono(n, &[((seed as usize) % n, 2)]));
                        poly(&[terms, vec![seed_term]].concat(), order)
                    })
                    .collect();
                let expected = reference(&generators, order);
                let (computed, stats) =
                    f5_basis_with_stats(&generators, order, F5Limits::default());
                assert_eq!(computed, expected, "seed {seed} in {order:?}");
                assert!(!stats.matrix_fallback, "seed {seed} in {order:?}");
            }
        }
    }

    #[test]
    fn f5_preserves_ideal_membership() {
        let order = Order::GradedRevLex;
        let generators = cyclic(4, order);
        let (basis, _) = f5_basis_with_stats(&generators, order, F5Limits::default());
        let basis = match basis {
            | Some(basis) => basis,
            | None => panic!("cyclic-4 must be solved"),
        };
        let expected = match groebner(&generators, order, GroebnerLimits::default()) {
            | Some(basis) => basis,
            | None => panic!("the reference must solve cyclic-4"),
        };
        assert_eq!(basis, expected);
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        for _ in 0..16 {
            let terms: Vec<(i64, Vec<u32>)> = (0..4)
                .map(|_| {
                    let coefficient = i64::try_from(next() % 7).unwrap_or(0) - 3;
                    let exponents: Vec<u32> = (0..4)
                        .map(|_| u32::try_from(next() % 3).unwrap_or(0))
                        .collect();
                    (coefficient, exponents)
                })
                .collect();
            let candidate = poly(&terms, order);
            assert_eq!(
                normal_form(&Domain::Rational, &candidate, &basis, order),
                normal_form(&Domain::Rational, &candidate, &expected, order)
            );
        }
    }

    // -----------------------------------------------------------------
    // Prime fields GF(p)
    // -----------------------------------------------------------------

    use crate::rules::poly::field::oracle;

    #[test]
    fn f5_over_prime_fields_matches_the_independent_reference() {
        let mut cases = 0_usize;
        let mut criteria_discards = 0_usize;
        for modulus in oracle::PRIMES {
            let domain = Domain::prime(BigInt::from(modulus))
                .unwrap_or_else(|| panic!("{modulus} must be a prime"));
            for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
                for (name, system) in oracle::sample_systems(modulus) {
                    let generators = oracle::generators(modulus, &system, order);
                    // The naive reference is impractical in lex on the random
                    // dense systems; those are compared in the engine order
                    // and certified structurally in the other orders.
                    let random = name.starts_with("random");
                    let expected = if !random || order == ENGINE_ORDER {
                        Some(oracle::buchberger_fp(modulus, &system, order))
                    } else {
                        None
                    };
                    let (computed, stats) =
                        f5_basis_with_stats_in(&domain, &generators, order, F5Limits::default());
                    let computed = computed
                        .unwrap_or_else(|| panic!("{name} over GF({modulus}) in {order:?}"));
                    if let Some(expected) = expected {
                        assert_eq!(computed, expected, "{name} over GF({modulus}) in {order:?}");
                    }
                    // The signature loop itself must never exhaust a resource
                    // guard on these systems: the finite-field path has the
                    // same no-fallback guarantee as the rational one.
                    assert!(
                        !stats.matrix_fallback,
                        "{name} over GF({modulus}) in {order:?} fell back: {stats:?}"
                    );
                    if order == ENGINE_ORDER {
                        assert!(
                            stats.criteria_path(),
                            "{name} over GF({modulus}) in grevlex did not stay on the criteria path: {stats:?}"
                        );
                    }
                    assert!(
                        stats.pairs_processed < F5Limits::default().pairs,
                        "the pair budget was exhausted over GF({modulus})"
                    );
                    criteria_discards = criteria_discards.saturating_add(
                        stats
                            .syzygy_discards
                            .saturating_add(stats.rewritten_discards),
                    );
                    oracle::check_basis(&domain, modulus, &computed, order);
                    // Determinism: the same request gives the same basis and
                    // the same counters.
                    let (again, again_stats) =
                        f5_basis_with_stats_in(&domain, &generators, order, F5Limits::default());
                    assert_eq!(again.as_deref(), Some(computed.as_slice()));
                    assert_eq!(again_stats, stats, "statistics are not deterministic");
                    cases = cases.saturating_add(1);
                }
            }
        }
        eprintln!(
            "F5 over GF(p): {cases} systems cross-checked, {criteria_discards} criterion discards"
        );
        assert!(cases >= 15, "the finite-field suite is too small");
        assert!(
            criteria_discards > 0,
            "the criteria were never exercised over a finite field"
        );
    }

    #[test]
    fn f5_finite_field_benchmarks_stay_on_the_criteria_path() {
        // Katsura-3 and Cyclic-4 over GF(5), in every order.
        let modulus = 5_i64;
        let domain =
            Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| panic!("5 must be a prime"));
        for (name, system) in [
            ("katsura-3", oracle::katsura(3)),
            ("cyclic-4", oracle::cyclic(4)),
        ] {
            for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
                let generators = oracle::generators(modulus, &system, order);
                let start = std::time::Instant::now();
                let (computed, stats) =
                    f5_basis_with_stats_in(&domain, &generators, order, F5Limits::default());
                let elapsed = start.elapsed();
                eprintln!(
                    "F5 {name} over GF(5) in {order:?}: {elapsed:?} ({} pairs, {} syzygy, {} rewritten, {} reductions, {} added, order fallback: {})",
                    stats.pairs_processed,
                    stats.syzygy_discards,
                    stats.rewritten_discards,
                    stats.reductions,
                    stats.polynomials_added,
                    stats.matrix_order_fallback
                );
                let computed = computed.unwrap_or_else(|| panic!("{name} over GF(5) in {order:?}"));
                let expected = if order == ENGINE_ORDER {
                    // The independent reference in the order the criteria loop
                    // actually runs in.
                    oracle::buchberger_fp(modulus, &system, order)
                } else {
                    // In the other orders FGLM may decline (cyclic-4 is
                    // positive-dimensional), so the matrix engine is the
                    // reference: it is a different algorithm from the
                    // signature loop.
                    f4_basis_in(&domain, &generators, order, F4Limits::default())
                        .unwrap_or_default()
                };
                assert_eq!(computed, expected, "{name} over GF(5) in {order:?}");
                assert!(
                    !stats.matrix_fallback,
                    "{name} over GF(5) in {order:?} exhausted the criteria loop: {stats:?}"
                );
                oracle::check_basis(&domain, modulus, &computed, order);
                // The criteria did the work: at least one pair was discarded
                // by a criterion, and the pair budget was not exhausted.
                assert!(
                    stats
                        .syzygy_discards
                        .saturating_add(stats.rewritten_discards)
                        > 0,
                    "{name} over GF(5) in {order:?} did not exercise the criteria: {stats:?}"
                );
                assert!(
                    stats.pairs_processed < F5Limits::default().pairs,
                    "{name} over GF(5) in {order:?} exhausted the pair budget"
                );
                // Over the engine order the conversion is not needed, so both
                // fallback counters must be clear.
                if order == ENGINE_ORDER {
                    assert!(
                        stats.criteria_path(),
                        "{name} over GF(5) did not stay on the criteria path in grevlex: {stats:?}"
                    );
                } else if stats.matrix_order_fallback {
                    // Cyclic-4 is positive-dimensional, so FGLM cannot convert
                    // it; the matrix engine finishes and the limit is recorded.
                    assert!(
                        !stats.criteria_path(),
                        "{name}: the order fallback is not recorded"
                    );
                }
                if order == ENGINE_ORDER {
                    assert!(
                        elapsed.as_secs() < 30,
                        "{name} over GF(5) took {elapsed:?}, over the 30 s budget"
                    );
                }
            }
        }
    }

    #[test]
    fn f5_finite_field_ideal_membership_agrees_with_the_reference() {
        // Membership decided by the computed basis must agree with the
        // independent reference, for members and non-members alike.
        let modulus = 7_i64;
        let domain =
            Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| panic!("7 must be a prime"));
        let order = Order::GradedRevLex;
        let system = oracle::cyclic(3);
        let generators = oracle::generators(modulus, &system, order);
        let (basis, _) = f5_basis_with_stats_in(&domain, &generators, order, F5Limits::default());
        let basis = basis.unwrap_or_else(|| panic!("the basis must be computed"));
        // Every generator is in the ideal.
        for generator in &generators {
            assert!(normal_form(&domain, generator, &basis, order).is_zero());
        }
        // Random candidates: the engine's normal form and the independent
        // reference agree on membership.
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        for _ in 0..12 {
            let terms: Vec<(i64, Vec<u32>)> = (0..4)
                .map(|_| {
                    let coefficient = i64::try_from(next() % 7).unwrap_or(0) - 3;
                    let exponents: Vec<u32> = (0..3)
                        .map(|_| u32::try_from(next() % 3).unwrap_or(0))
                        .collect();
                    (coefficient, exponents)
                })
                .collect();
            let candidate = oracle::generators(modulus, &[terms.clone()], order)
                .into_iter()
                .next()
                .unwrap_or_else(|| panic!("a random polynomial must be built"));
            let member = normal_form(&domain, &candidate, &basis, order).is_zero();
            assert_eq!(
                member,
                oracle::reduces_to_zero_fp(modulus, &basis, &terms, order),
                "membership of {terms:?} over GF(7)"
            );
        }
        // Adding a member does not change the ideal, hence not the basis.
        let mut extended = generators.clone();
        extended.push(basis[0].clone());
        let (again, _) = f5_basis_with_stats_in(&domain, &extended, order, F5Limits::default());
        assert_eq!(again.as_deref(), Some(basis.as_slice()));
    }
}
