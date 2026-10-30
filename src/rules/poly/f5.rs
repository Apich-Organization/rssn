//! F5: the Faugère signature-based algorithm for Gröbner bases.
//!
//! Every polynomial carries the signature of the S-pair it came from: a
//! module monomial together with the index of the input generator the
//! representation starts from. Signatures are ordered *position over term*
//! with the smallest generator index first, so the signature of a polynomial
//! is the module head term of a representation. Pairs are processed in
//! increasing signature order and three criteria discard the ones that cannot
//! contribute:
//!
//! * the **F5 criterion** — the scaled signature term `t * signature(l)` of
//!   one side is divisible by the leading monomial of a polynomial with a
//!   smaller generator index;
//! * the **rewritten criterion** — the signature is a multiple of the
//!   signature of a pair already processed (same index), so the S-polynomial
//!   rewrites to zero;
//! * the **rewritable criterion** — a pair still in the queue has a smaller
//!   signature dividing this one, so this pair is rewritten before it is
//!   reached.
//!
//! The remaining pairs are reduced with **signature-admissible** reductions:
//! a reducer `h` may only be used with a multiplier `t` such that
//! `t * signature(h) <= signature(p)`. Every polynomial therefore keeps a
//! signature that bounds all of its representation, which is what the
//! criteria need.
//!
//! # Affine input
//!
//! F5 needs a homogeneous ideal in a degree-compatible order. Non-homogeneous
//! generators are homogenised with a fresh variable placed last (the smallest
//! variable of the graded reverse-lexicographic order), the basis is computed
//! there, and the fresh variable is set to one again. The dehomogenisation of
//! a graded Gröbner basis is a Gröbner basis of the original ideal in the
//! restricted order, so this is the standard projective-closure route.
//!
//! The computation itself runs in graded reverse-lex. For the other orders
//! the reduced basis is converted with
//! [`fglm`](crate::rules::poly::groebner::fglm) — the same strategy as
//! `groebner_via` — and when FGLM declines (a positive-dimensional ideal, or
//! a quotient ring beyond its dimension bound) the matrix engine computes the
//! basis in the requested order.
//!
//! # Safety of the criteria, and of the result
//!
//! The textbook proofs of the three criteria are tied to the reduction
//! strategy of the original algorithm, while this engine generalises it (any
//! basis element may reduce, as long as the signature does not grow). Rather
//! than trusting a criterion blindly, the engine **certifies** the set it
//! produces:
//!
//! 1. the basis is interreduced — one element at a time, replacing it by its
//!    normal form modulo the current set and dropping it only when it then
//!    lies in the ideal of the elements that remain, which keeps the ideal
//!    exactly;
//! 2. every S-pair is reduced, a non-zero remainder being added and paired
//!    with everything before it (Buchberger completion).
//!
//! The generators are always part of the basis, so the ideal of the completed
//! set is the requested ideal and its reduced basis is the unique answer. On
//! a correct F5 run the completion adds nothing.
//!
//! Certification can be expensive when a criterion discarded too much, so the
//! engine runs three tiers, each of which produces the same reduced basis:
//! the criteria-driven loop with a tight budget; the complete pair set (the
//! provably safe subset) with the operator's own, larger, budgets; and
//! finally the matrix engine. The tier that produced an answer therefore only
//! affects the run time, never the result.
//!
//! All arithmetic is exact over `Q`, selection is deterministic (the smallest
//! signature, ties broken by the pair indices), and exhausted budgets report
//! the request as unevaluated instead of returning a partial basis.

// The signature machinery is one algorithm: pair construction, the two
// criteria, signature-admissible reduction and the certification pass. Each
// helper serves exactly one step.
#![allow(
    clippy::single_call_fn,
    clippy::too_many_lines,
    clippy::too_many_arguments
)]

use std::cmp::Ordering;
use std::collections::BTreeSet;

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

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
use crate::rules::poly::f4::f4_basis;
use crate::rules::poly::groebner::GPoly;
use crate::rules::poly::groebner::Order;
use crate::rules::poly::groebner::Term;
use crate::rules::poly::groebner::fglm;
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
    /// Maximum number of reductions of the certification pass.
    pub certification: usize,
    /// Pair budget of the first, criteria-driven attempt. When it is exceeded
    /// the engine reruns with every critical pair.
    pub fast_pairs: usize,
    /// Certification budget of the first attempt.
    pub fast_certification: usize,
    /// Reduction budget of the first attempt, kept small so that an
    /// over-aggressive criterion is detected quickly.
    pub fast_reductions: usize,
    /// Which criteria the loop applies. The production settings are all on;
    /// the unit tests switch individual criteria off to isolate one that a
    /// counterexample might implicate.
    pub strategy: Strategy,
}

/// The three criteria of the signature-based loop.
#[derive(Copy, Clone, Debug)]
pub struct Strategy {
    /// The F5 criterion (a smaller-index leading monomial divides the
    /// scaled signature term).
    pub f5: bool,
    /// The rewritten criterion (the signature is a multiple of a processed
    /// signature with the same index).
    pub rewritten: bool,
    /// The rewritable criterion (a smaller pending signature divides this
    /// one).
    pub rewritable: bool,
}

impl Default for Strategy {
    fn default() -> Self {
        Self {
            f5: true,
            rewritten: true,
            rewritable: true,
        }
    }
}

impl Strategy {
    /// The provably complete subset: every critical pair is processed.
    pub const fn complete() -> Self {
        Self {
            f5: false,
            rewritten: false,
            rewritable: false,
        }
    }
}

impl Default for F5Limits {
    fn default() -> Self {
        Self {
            pairs: 5_000,
            reductions: 50_000,
            basis: 400,
            certification: 1_000,
            fast_pairs: 300,
            fast_certification: 40,
            fast_reductions: 2_000,
            strategy: Strategy::default(),
        }
    }
}

/// A polynomial together with its signature: a module monomial and the index
/// of the generator it descends from.
#[derive(Clone, Debug)]
struct Labelled {
    poly: GPoly,
    sig: common::Mono,
    index: usize,
}

/// A critical pair waiting in the queue.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Pair {
    /// Index of the first basis element.
    a: usize,
    /// Index of the second basis element.
    b: usize,
    /// Module monomial of the S-polynomial's signature.
    sig: common::Mono,
    /// Generator index of the signature.
    index: usize,
}

/// The signature order: **position over term with the smallest index first**.
/// A signature with a smaller generator index is larger, and signatures with
/// the same index are compared by their module monomial. This is the order of
/// the signature-based algorithms: the signature of a polynomial is the
/// module head term of a representation, which is the term with the smallest
/// index.
fn signature_order(
    a: &[u32],
    a_index: usize,
    b: &[u32],
    b_index: usize,
) -> Ordering {
    b_index.cmp(&a_index).then_with(|| ENGINE_ORDER.cmp(a, b))
}

/// Total order on the queue: signature, then the pair indices so that equal
/// signatures are still handled deterministically.
fn compare(
    a: &Pair,
    b: &Pair,
) -> Ordering {
    signature_order(&a.sig, a.index, &b.sig, b.index)
        .then_with(|| a.a.cmp(&b.a))
        .then_with(|| a.b.cmp(&b.b))
}

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
        compare(self, other)
    }
}

/// Builds the pair `(a, b)` with the signature of its S-polynomial:
/// the larger of `(lcm/lm f) * signature(f)` and `(lcm/lm g) * signature(g)`.
fn make_pair(
    basis: &[Labelled],
    a: usize,
    b: usize,
) -> Option<Pair> {
    let (first, second) = (basis.get(a)?, basis.get(b)?);
    let (first_lead, _) = first.poly.leading()?;
    let (second_lead, _) = second.poly.leading()?;
    let lcm = common::lcm(first_lead, second_lead);
    let first_shift = common::multiply(&first.sig, &common::difference(&lcm, first_lead));
    let second_shift = common::multiply(&second.sig, &common::difference(&lcm, second_lead));
    let (sig, index) = if signature_order(&first_shift, first.index, &second_shift, second.index)
        == Ordering::Less
    {
        (second_shift, second.index)
    } else {
        (first_shift, first.index)
    };
    Some(Pair { a, b, sig, index })
}

/// The rewritten criterion: the signature is a multiple of the signature of a
/// pair that has already been processed (with the same generator index, as
/// multiplication of a signature never changes its index).
fn rewritten(
    known: &[(common::Mono, usize)],
    sig: &[u32],
    index: usize,
) -> bool {
    known
        .iter()
        .any(|(monomial, known_index)| *known_index == index && common::divides(monomial, sig))
}

/// The rewritable criterion: another pair in the queue has a smaller
/// signature that divides this one, so this pair is rewritten before it is
/// reached. The queue is ordered by signature, so only the entries before the
/// candidate can qualify.
fn rewritable(
    queue: &BTreeSet<Pair>,
    candidate: &Pair,
) -> bool {
    queue
        .iter()
        .take_while(|other| compare(other, candidate) == Ordering::Less)
        .any(|other| other.index == candidate.index && common::divides(&other.sig, &candidate.sig))
}

/// The F5 criterion for one side of a pair: the scaled signature term
/// `t * signature(l)` is divisible by the leading monomial of a polynomial
/// whose generator index is smaller than the side's. Such a pair rewrites to
/// zero and is discarded.
fn f5_criterion(
    basis: &[Labelled],
    a: usize,
    b: usize,
) -> bool {
    let (Some(first), Some(second)) = (basis.get(a), basis.get(b)) else {
        return false;
    };
    let (Some((first_lead, _)), Some((second_lead, _))) =
        (first.poly.leading(), second.poly.leading())
    else {
        return false;
    };
    let lcm = common::lcm(first_lead, second_lead);
    let sides = [
        (a, common::difference(&lcm, first_lead)),
        (b, common::difference(&lcm, second_lead)),
    ];
    sides.iter().any(|(side, multiplier)| {
        let Some(element) = basis.get(*side) else {
            return false;
        };
        let scaled = common::multiply(&element.sig, multiplier);
        basis.iter().any(|other| {
            other.index < element.index
                && other
                    .poly
                    .leading()
                    .is_some_and(|(lead, _)| common::divides(lead, &scaled))
        })
    })
}

/// A polynomial in the signature-driven loop must not grow past this many
/// terms. F5 keeps the unreduced tail of an S-polynomial (only the leading
/// term is reduced at each step), so on an input that is not close to a
/// Gröbner basis the coefficients and the number of terms can grow quickly.
/// Hitting the bound reports budget exhaustion, which makes the caller fall
/// back to a cheaper engine instead of running for minutes.
const MAX_TERMS: usize = 1_200;

/// Likewise for the size of a coefficient: exact rational arithmetic keeps
/// working however large the numerator and denominator become, but the cost
/// per operation grows with them, so a bound keeps the run time predictable.
const MAX_BITS: u64 = 4_096;

/// Whether any coefficient of the working polynomial has outgrown
/// [`MAX_BITS`].
fn coefficients_too_large(terms: &[Term]) -> bool {
    terms
        .iter()
        .any(|(_, c)| c.numer().bits() > MAX_BITS || c.denom().bits() > MAX_BITS)
}

/// Reduces `p` modulo `basis`, allowing a reduction by `t h` only when
/// `t * signature(h)` is at most `bound`. The result has signature at most
/// `bound`, which is what the criteria rely on. Returns `None` when the
/// reduction budget is exhausted or the working polynomial outgrows
/// [`MAX_TERMS`].
///
/// The working polynomial is kept as a vector of terms in descending order
/// and combined with a reducer row by a linear merge; multiplying a monomial
/// by the fixed multiplier preserves the order, so the merge compares the
/// unterminated exponents and only shifts the coefficients it writes out.
fn signature_normal_form(
    p: &GPoly,
    bound: &[u32],
    bound_index: usize,
    basis: &[Labelled],
    steps: &mut usize,
    limit: usize,
) -> Option<GPoly> {
    let mut work: Vec<Term> = p.terms().to_vec();
    let mut remainder: Vec<Term> = Vec::new();
    let mut scratch: Vec<Term> = Vec::new();
    loop {
        if work.len() > MAX_TERMS || coefficients_too_large(&work) {
            return None;
        }
        let Some((monomial, coefficient)) = work.first().cloned() else {
            break;
        };
        let mut reducer: Option<(&GPoly, common::Mono, BigRational)> = None;
        for element in basis {
            let Some((lead, lead_coefficient)) = element.poly.leading() else {
                continue;
            };
            if !common::divides(lead, &monomial) {
                continue;
            }
            let multiplier = common::difference(&monomial, lead);
            let shifted = common::multiply(&element.sig, &multiplier);
            if signature_order(&shifted, element.index, bound, bound_index) == Ordering::Greater {
                continue;
            }
            reducer = Some((&element.poly, multiplier, lead_coefficient.clone()));
            break;
        }
        match reducer {
            | Some((poly, multiplier, lead_coefficient)) => {
                *steps = steps.saturating_add(1);
                if *steps > limit {
                    return None;
                }
                let factor = coefficient / lead_coefficient;
                merge_scaled(&mut work, &mut scratch, poly, &multiplier, &factor);
            },
            | None => {
                remainder.push((monomial, coefficient));
                work.remove(0);
            },
        }
    }
    Some(GPoly::new(remainder, ENGINE_ORDER))
}

/// `work <- work - factor * x^multiplier * poly`, both polynomials being
/// vectors of terms in descending order. The comparison during the merge is
/// against the *shifted* monomials (a fixed multiplier does not preserve a
/// comparison with an unshifted monomial), which
/// [`compare_shifted`] evaluates without materialising them.
fn merge_scaled(
    work: &mut Vec<Term>,
    scratch: &mut Vec<Term>,
    poly: &GPoly,
    multiplier: &[u32],
    factor: &BigRational,
) {
    scratch.clear();
    scratch.reserve(work.len().saturating_add(poly.terms().len()));
    let mut i = 0_usize;
    let mut j = 0_usize;
    loop {
        match (work.get(i), poly.terms().get(j)) {
            | (Some((left, left_coefficient)), Some((right, right_coefficient))) => {
                match compare_shifted(left, right, multiplier) {
                    | Ordering::Greater => {
                        scratch.push((left.clone(), left_coefficient.clone()));
                        i = i.saturating_add(1);
                    },
                    | Ordering::Less => {
                        scratch.push((
                            common::multiply(right, multiplier),
                            -(right_coefficient * factor),
                        ));
                        j = j.saturating_add(1);
                    },
                    | Ordering::Equal => {
                        let value = left_coefficient - factor * right_coefficient;
                        if !value.is_zero() {
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
                    -(right_coefficient * factor),
                ));
                j = j.saturating_add(1);
            },
            | (None, None) => break,
        }
    }
    std::mem::swap(work, scratch);
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

/// The ordinary normal form of `p` modulo `basis`, used by the certification
/// pass. Same merge as [`signature_normal_form`], without the signature bound.
fn plain_normal_form(
    p: &GPoly,
    basis: &[GPoly],
    steps: &mut usize,
    limit: usize,
) -> Option<GPoly> {
    let mut work: Vec<Term> = p.terms().to_vec();
    let mut remainder: Vec<Term> = Vec::new();
    let mut scratch: Vec<Term> = Vec::new();
    loop {
        if work.len() > MAX_TERMS || coefficients_too_large(&work) {
            return None;
        }
        let Some((monomial, coefficient)) = work.first().cloned() else {
            break;
        };
        let mut reducer: Option<(&GPoly, common::Mono, BigRational)> = None;
        for element in basis {
            let Some((lead, lead_coefficient)) = element.leading() else {
                continue;
            };
            if !common::divides(lead, &monomial) {
                continue;
            }
            reducer = Some((
                element,
                common::difference(&monomial, lead),
                lead_coefficient.clone(),
            ));
            break;
        }
        match reducer {
            | Some((poly, multiplier, lead_coefficient)) => {
                *steps = steps.saturating_add(1);
                if *steps > limit {
                    return None;
                }
                let factor = coefficient / lead_coefficient;
                merge_scaled(&mut work, &mut scratch, poly, &multiplier, &factor);
            },
            | None => {
                remainder.push((monomial, coefficient));
                work.remove(0);
            },
        }
    }
    Some(GPoly::new(remainder, ENGINE_ORDER))
}

/// The signature-driven loop. `generators` are homogeneous and monic; the
/// result is a Gröbner basis in the graded reverse-lexicographic order (not
/// necessarily reduced or minimal).
fn f5_loop(
    generators: &[GPoly],
    limits: F5Limits,
) -> Option<Vec<GPoly>> {
    let dimension = generators
        .iter()
        .find_map(|g| g.leading().map(|(m, _)| m.len()));
    let dimension = dimension?;
    let unit = vec![0_u32; dimension];
    let mut basis: Vec<Labelled> = generators
        .iter()
        .enumerate()
        .map(|(index, g)| Labelled {
            poly: common::monic(g, ENGINE_ORDER),
            sig: unit.clone(),
            index: index.saturating_add(1),
        })
        .collect();
    let mut queue: BTreeSet<Pair> = BTreeSet::new();
    for b in 0..basis.len() {
        for a in 0..b {
            if let Some(pair) = make_pair(&basis, a, b) {
                queue.insert(pair);
            }
        }
    }
    let mut known: Vec<(common::Mono, usize)> = Vec::new();
    let mut steps = 0_usize;
    let mut processed = 0_usize;
    // The queue is ordered by signature, so the first entry is always the
    // smallest signature; the choice is deterministic and logarithmic.
    while let Some(pair) = queue.pop_first() {
        processed = processed.saturating_add(1);
        if processed > limits.pairs {
            return None;
        }
        if limits.strategy.rewritten && rewritten(&known, &pair.sig, pair.index) {
            continue;
        }
        if limits.strategy.f5 && f5_criterion(&basis, pair.a, pair.b) {
            // The F5 criterion proves that this S-polynomial is a syzygy, so
            // its signature is a *known* signature from now on: later pairs
            // whose signature is a multiple of it may be rewritten (and a
            // pair that was dropped as "rewritable" because this one was
            // still pending relies on exactly that).
            known.push((pair.sig.clone(), pair.index));
            continue;
        }
        let (Some(first), Some(second)) = (basis.get(pair.a), basis.get(pair.b)) else {
            continue;
        };
        let spolynomial = common::s_polynomial(&first.poly, &second.poly, ENGINE_ORDER);
        let remainder = signature_normal_form(
            &spolynomial,
            &pair.sig,
            pair.index,
            &basis,
            &mut steps,
            limits.reductions,
        )?;
        known.push((pair.sig.clone(), pair.index));
        if remainder.is_zero() {
            continue;
        }
        if basis.len() >= limits.basis {
            return None;
        }
        basis.push(Labelled {
            poly: common::monic(&remainder, ENGINE_ORDER),
            sig: pair.sig,
            index: pair.index,
        });
        let new = basis.len().saturating_sub(1);
        for other in 0..new {
            if let Some(candidate) = make_pair(&basis, other, new)
                && (!limits.strategy.rewritable || !rewritable(&queue, &candidate))
            {
                queue.insert(candidate);
            }
        }
    }
    Some(basis.into_iter().map(|element| element.poly).collect())
}

/// Interreduces a generating set without changing the ideal it generates.
///
/// Elements are handled **one at a time** against the current set: an element
/// is replaced by its normal form modulo the others (which preserves the
/// ideal), and an element whose normal form is zero is dropped *because it
/// then lies in the ideal of the elements that remain*. Dropping several
/// elements in one sweep against the original set would be wrong: with
/// `f1 = f2 + h` both `f1` and `f2` reduce to zero against the others while
/// `h` alone generates a strictly smaller ideal.
///
/// The pass is what makes the certification cheap: on a basis that is
/// already a Gröbner basis it produces the reduced basis, whose coefficients
/// are far smaller than those of the raw signature output.
fn interreduce(
    mut basis: Vec<GPoly>,
    steps: &mut usize,
    limit: usize,
) -> Option<Vec<GPoly>> {
    let mut changed = true;
    let mut passes = 0_usize;
    while changed && passes < 5 {
        changed = false;
        passes = passes.saturating_add(1);
        let mut index = 0_usize;
        while index < basis.len() {
            let others: Vec<GPoly> = basis
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != index)
                .map(|(_, g)| g.clone())
                .collect();
            let Some(g) = basis.get(index).cloned() else {
                break;
            };
            let reduced = plain_normal_form(&g, &others, steps, limit)?;
            if reduced.is_zero() {
                basis.remove(index);
                changed = true;
                continue;
            }
            let reduced = common::monic(&reduced, ENGINE_ORDER);
            if reduced != g {
                basis[index] = reduced;
                changed = true;
            }
            index = index.saturating_add(1);
        }
    }
    Some(basis)
}

/// Buchberger completion of `basis`: the set is interreduced first and then
/// every S-pair is reduced, a non-zero remainder being added and paired with
/// everything before it. A correct F5 run finds nothing to add; the pass
/// certifies the set and is what makes the returned basis a Gröbner basis
/// even if a signature criterion were too aggressive. The generators are
/// always part of the basis, so the completed ideal is the requested one.
/// `None` when a budget is exhausted.
fn certify(
    basis: Vec<GPoly>,
    pair_budget: usize,
    reduction_budget: usize,
    basis_limit: usize,
) -> Option<Vec<GPoly>> {
    let mut reductions = 0_usize;
    let mut basis = interreduce(basis, &mut reductions, reduction_budget)?;
    let mut pairs: Vec<(usize, usize)> = (0..basis.len())
        .flat_map(|b| (0..b).map(move |a| (a, b)))
        .collect();
    let mut steps = 0_usize;
    while let Some((a, b)) = pairs.pop() {
        steps = steps.saturating_add(1);
        let (Some(first), Some(second)) = (basis.get(a), basis.get(b)) else {
            continue;
        };
        let (Some((first_lead, _)), Some((second_lead, _))) = (first.leading(), second.leading())
        else {
            continue;
        };
        if common::coprime(first_lead, second_lead) {
            continue;
        }
        let lcm = common::lcm(first_lead, second_lead);
        // The chain criterion, against the pairs still in the queue.
        let covered = basis.iter().enumerate().any(|(k, h)| {
            k != a
                && k != b
                && h.leading().is_some_and(|(m, _)| common::divides(m, &lcm))
                && !pairs.contains(&common::ordered(k, a))
                && !pairs.contains(&common::ordered(k, b))
        });
        if covered {
            continue;
        }
        if steps > pair_budget {
            return None;
        }
        let spolynomial = common::s_polynomial(first, second, ENGINE_ORDER);
        let remainder = plain_normal_form(&spolynomial, &basis, &mut reductions, reduction_budget)?;
        if remainder.is_zero() {
            continue;
        }
        if basis.len() >= basis_limit {
            return None;
        }
        let new = basis.len();
        pairs.extend((0..new).map(|k| (k, new)));
        basis.push(common::monic(&remainder, ENGINE_ORDER));
    }
    Some(basis)
}

/// The reduced Gröbner basis of the ideal generated by `generators` in
/// `order` (or `None` when a limit is hit).
///
/// Homogeneous input is processed directly; affine input is homogenised with
/// a fresh variable and dehomogenised afterwards, as described in the module
/// documentation.
#[must_use]
pub fn f5_basis(
    generators: &[GPoly],
    order: Order,
    limits: F5Limits,
) -> Option<Vec<GPoly>> {
    let Some(dimension) = generators
        .iter()
        .find_map(|g| g.leading().map(|(m, _)| m.len()))
    else {
        return Some(Vec::new());
    };
    let input: Vec<GPoly> = generators
        .iter()
        .filter(|g| !g.is_zero())
        .map(|g| GPoly::new(g.terms().to_vec(), ENGINE_ORDER))
        .collect();
    if input.is_empty() {
        return Some(Vec::new());
    }
    // A non-zero constant generator generates the unit ideal. Only a
    // polynomial whose *leading* monomial is the constant (i.e. a constant)
    // counts: the constant term of a non-constant generator does not.
    let is_constant = |g: &GPoly| g.leading().is_some_and(|(m, _)| m.iter().all(|&e| e == 0));
    if input.iter().any(is_constant) {
        let unit = vec![0_u32; dimension];
        return Some(vec![GPoly::new(vec![(unit, BigRational::one())], order)]);
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
    // The criteria are an optimisation: on an input where one of them is too
    // aggressive, the certification pass has to add the missing elements back,
    // which can be far more work than processing every pair in the first
    // place. The engine therefore
    //
    // 1. runs with the criteria and a tight certification budget,
    // 2. reruns with every critical pair (the provably complete subset) when
    //    that budget is exceeded, and
    // 3. falls back to the matrix engine if even the complete signature run
    //    exhausts its own, much larger, limits.
    //
    // Every stage computes the same reduced basis; only the amount of work
    // differs.
    let attempt = |strategy: Strategy,
                   pairs: usize,
                   pairing: usize,
                   reductions: usize|
     -> Option<Vec<GPoly>> {
        let trial = F5Limits {
            strategy,
            pairs,
            certification: pairing,
            reductions,
            ..limits
        };
        f5_loop(&working, trial)
            .and_then(|output| certify(output, pairing, reductions, limits.basis))
    };
    let completion = attempt(
        limits.strategy,
        limits.fast_pairs.min(limits.pairs),
        limits.fast_certification.min(limits.certification),
        limits.fast_reductions.min(limits.reductions),
    )
    .or_else(|| {
        attempt(
            Strategy::complete(),
            limits.pairs,
            limits.certification,
            limits.reductions,
        )
    })
    .or_else(|| {
        f4_basis(
            &working,
            ENGINE_ORDER,
            F4Limits {
                basis: limits.basis,
                ..F4Limits::default()
            },
        )
    })?;
    let original: Vec<GPoly> = if homogenised {
        completion
            .iter()
            .map(|g| common::dehomogenise(g, ENGINE_ORDER))
            .collect()
    } else {
        completion
    };
    let reduced = common::reduced_basis(&original, ENGINE_ORDER);
    if order == ENGINE_ORDER {
        return Some(reduced);
    }
    if let Some(converted) = fglm(&reduced, ENGINE_ORDER, order) {
        return Some(converted);
    }
    // Positive-dimensional ideal or too large for FGLM: the matrix engine
    // computes the same reduced basis in the requested order.
    let reordered: Vec<GPoly> = original
        .iter()
        .map(|g| GPoly::new(g.terms().to_vec(), order))
        .collect();
    f4_basis(&reordered, order, F4Limits::default())
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
            | Some((generators, order, gens)) => Prepared::Work(Box::new(F5Work {
                generators,
                order,
                gens,
            })),
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
        let Some((generators, order, gens)) = common::prepare_data(cx, node) else {
            return Outcome::Pass;
        };
        match f5_basis(&generators, order, F5Limits::default()) {
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
    generators: Vec<GPoly>,
    order: Order,
    gens: Gens,
}

impl PreparedWork for F5Work {
    fn run(self: Box<Self>) -> Box<dyn std::any::Any + Send> {
        let basis = f5_basis(&self.generators, self.order, F5Limits::default());
        Box::new((basis, self.gens))
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::*;
    use crate::rules::poly::f4::common::Mono;
    use crate::rules::poly::groebner::GroebnerLimits;
    use crate::rules::poly::groebner::groebner;
    use crate::rules::poly::groebner::groebner_via;
    use crate::rules::poly::groebner::reduce as normal_form;

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

    #[test]
    fn signatures_order_by_the_smallest_index_first() {
        // e_1 is larger than every e_2 term, whatever the monomials.
        assert_eq!(
            signature_order(&[0, 0, 0], 1, &[5, 5, 5], 2),
            Ordering::Greater
        );
        assert_eq!(
            signature_order(&[1, 0, 0], 2, &[1, 0, 0], 2),
            Ordering::Equal
        );
        // Within one index, the monomial order decides.
        assert_eq!(
            signature_order(&[1, 0, 0], 1, &[0, 1, 0], 1),
            ENGINE_ORDER.cmp(&[1, 0, 0], &[0, 1, 0])
        );
    }

    /// The larger systems, checked in every order: F5 always runs in graded
    /// reverse-lex, so the order only changes the conversion (FGLM or the
    /// matrix fallback).
    fn large_systems(order: Order) -> Vec<(&'static str, Vec<GPoly>)> {
        vec![
            ("katsura-4", katsura(4, order)),
            ("cyclic-5", cyclic(5, order)),
        ]
    }

    #[test]
    fn f5_agrees_with_buchberger_on_the_sample_systems() {
        for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
            for (name, generators) in sample_systems(order) {
                let expected = groebner(&generators, order, GroebnerLimits::default());
                let computed = f5_basis(&generators, order, F5Limits::default());
                assert_eq!(computed, expected, "{name} in {order:?}");
            }
        }
        for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
            for (name, generators) in large_systems(Order::GradedRevLex) {
                // These orders are compared with the reference route the
                // `groebner` operator uses (graded reverse-lex plus FGLM);
                // F5 itself always runs in graded reverse-lex.
                let start = std::time::Instant::now();
                let expected = groebner_via(&generators, order, GroebnerLimits::default());
                let reference = start.elapsed();
                let start = std::time::Instant::now();
                let computed = f5_basis(&generators, order, F5Limits::default());
                eprintln!(
                    "F5 {name} in {order:?}: reference {reference:?}, F5 {:?}",
                    start.elapsed()
                );
                assert_eq!(computed, expected, "{name} in {order:?}");
            }
        }
    }

    #[test]
    fn f5_works_on_cyclic_five_in_both_orders() {
        for order in [Order::GradedRevLex, Order::Lex] {
            let generators = cyclic(5, order);
            let expected =
                groebner_via(&generators, order, GroebnerLimits::default()).unwrap_or_default();
            let computed = f5_basis(&generators, order, F5Limits::default()).unwrap_or_default();
            assert_eq!(computed, expected, "cyclic-5 in {order:?}");
        }
    }

    #[test]
    fn f5_homogenises_affine_systems() {
        // x^2 + y + 1, x - y: affine (mixed degrees), so the homogenisation
        // path runs; the generic sample already covers it, this one pins the
        // shape of the answer.
        let order = Order::GradedRevLex;
        let generators = vec![
            poly(
                &[(1, vec![2, 0, 0]), (1, vec![0, 1, 0]), (1, vec![0, 0, 0])],
                order,
            ),
            poly(&[(1, vec![1, 0, 0]), (-1, vec![0, 1, 0])], order),
        ];
        let computed = f5_basis(&generators, order, F5Limits::default()).unwrap_or_default();
        let expected = groebner(&generators, order, GroebnerLimits::default()).unwrap_or_default();
        assert_eq!(computed, expected);
        assert!(!computed.is_empty());
    }

    #[test]
    fn every_criterion_variant_returns_the_same_basis() {
        // Switching criteria off must not change the answer: the generators
        // are always part of the basis and the certification pass completes
        // it to a Gröbner basis of the same ideal, so every strategy ends at
        // the same reduced basis.
        let order = Order::GradedRevLex;
        for (name, generators) in [
            ("katsura-3", katsura(3, order)),
            ("cyclic-4", cyclic(4, order)),
        ] {
            let expected = groebner(&generators, order, GroebnerLimits::default());
            for (f5, rewritten, rewritable) in [
                (false, false, false),
                (true, true, true),
                (false, true, true),
                (true, false, true),
                (true, true, false),
            ] {
                let limits = F5Limits {
                    // Measure each variant on its own: the fallback tiers are
                    // disabled by the enormous budgets.
                    certification: usize::MAX,
                    fast_certification: usize::MAX,
                    fast_reductions: usize::MAX,
                    reductions: 4_000_000,
                    strategy: Strategy {
                        f5,
                        rewritten,
                        rewritable,
                    },
                    ..F5Limits::default()
                };
                let computed = f5_basis(&generators, order, limits);
                assert_eq!(
                    computed, expected,
                    "{name}: f5={f5} rewritten={rewritten} rewritable={rewritable}"
                );
            }
        }
    }

    #[test]
    fn f5_solves_the_large_benchmarks() {
        // Katsura-4 and cyclic-5: the signature loop with the criteria runs
        // out of its tight budget, the complete-pair rerun or the matrix
        // fallback finishes the job, and the answer is the same one
        // Buchberger computes.
        let order = Order::GradedRevLex;
        for (name, generators) in [
            ("katsura-4", katsura(4, order)),
            ("cyclic-5", cyclic(5, order)),
        ] {
            let start = std::time::Instant::now();
            let computed = f5_basis(&generators, order, F5Limits::default());
            eprintln!(
                "F5 {name}: {:?} in {:?}",
                computed.as_ref().map(Vec::len),
                start.elapsed()
            );
            let expected = groebner_via(&generators, order, GroebnerLimits::default());
            assert_eq!(computed, expected, "{name}");
        }
    }

    #[test]
    fn f5_limits_report_exhaustion() {
        let order = Order::GradedRevLex;
        let generators = katsura(3, order);
        // The signature loop itself reports exhaustion.
        for limits in [
            F5Limits {
                pairs: 0,
                ..F5Limits::default()
            },
            F5Limits {
                reductions: 0,
                ..F5Limits::default()
            },
            F5Limits {
                basis: 1,
                ..F5Limits::default()
            },
        ] {
            assert_eq!(f5_loop(&generators, limits), None, "{limits:?}");
        }
        // So does the certification pass.
        let basis = f5_loop(&generators, F5Limits::default()).unwrap_or_default();
        assert_eq!(certify(basis.clone(), 0, 0, 400), None);
        assert_eq!(certify(basis.clone(), usize::MAX, 0, 400), None);
        assert_eq!(certify(basis.clone(), usize::MAX, usize::MAX, 0), None);
        // With room to work, the operator returns the reduced basis even when
        // its own limits are tight: the fallback tiers take over.
        assert!(f5_basis(&generators, order, F5Limits::default()).is_some());
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
    }

    #[test]
    fn signature_reduction_never_raises_the_signature() {
        // f = x^2, g = x*y in grevlex with generators indexed 1 and 2.
        let basis = vec![
            Labelled {
                poly: poly(&[(1, vec![2, 0])], ENGINE_ORDER),
                sig: vec![0, 0],
                index: 1,
            },
            Labelled {
                poly: poly(&[(1, vec![1, 1])], ENGINE_ORDER),
                sig: vec![0, 0],
                index: 2,
            },
        ];
        // Reducing y*f = x^2 y: with a bound of index 2 neither reducer is
        // admissible, because y*e_1 and x*e_2 are both larger than e_2 in
        // the signature order (a smaller index is a larger signature).
        let p = poly(&[(1, vec![2, 1])], ENGINE_ORDER);
        let bound = vec![0, 0];
        let mut steps = 0;
        let reduced = match signature_normal_form(&p, &bound, 2, &basis, &mut steps, 100) {
            | Some(reduced) => reduced,
            | None => panic!("the reduction must fit its budget"),
        };
        assert_eq!(reduced, p);
        // With a bound of index 1 the reducer with index 2 is admissible and
        // the same polynomial reduces to zero.
        let mut steps = 0;
        let reduced = match signature_normal_form(&p, &bound, 1, &basis, &mut steps, 100) {
            | Some(reduced) => reduced,
            | None => panic!("the reduction must fit its budget"),
        };
        assert!(reduced.is_zero(), "y*f = x*(x*y) rewrites to zero");
    }

    #[test]
    fn f5_reduces_random_polynomials_like_buchberger() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        let order = Order::GradedRevLex;
        let n = 3;
        let generators = vec![
            poly(
                &[
                    (1, mono(n, &[(0, 2)])),
                    (1, mono(n, &[(1, 1)])),
                    (-1, mono(n, &[(2, 1)])),
                ],
                order,
            ),
            poly(
                &[
                    (1, mono(n, &[(1, 2)])),
                    (1, mono(n, &[(0, 1), (2, 1)])),
                    (3, mono(n, &[])),
                ],
                order,
            ),
            poly(
                &[
                    (1, mono(n, &[(2, 2)])),
                    (-2, mono(n, &[(0, 1), (1, 1)])),
                    (1, mono(n, &[])),
                ],
                order,
            ),
        ];
        let basis = f5_basis(&generators, order, F5Limits::default()).unwrap_or_default();
        let expected = groebner(&generators, order, GroebnerLimits::default()).unwrap_or_default();
        assert_eq!(basis, expected);
        for _ in 0..16 {
            let terms: Vec<(i64, Vec<u32>)> = (0..4)
                .map(|_| {
                    let c = i64::try_from(next() % 7).unwrap_or(0) - 3;
                    let exponents: Vec<u32> = (0..n)
                        .map(|_| u32::try_from(next() % 3).unwrap_or(0))
                        .collect();
                    (c, exponents)
                })
                .collect();
            let candidate = poly(&terms, order);
            assert_eq!(
                normal_form(&candidate, &basis, order),
                normal_form(&candidate, &expected, order)
            );
        }
    }
}
