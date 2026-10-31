//! Cylindrical algebraic decomposition by exact projection and lifting.
//!
//! The decomposition is built over the variables in the order given (the
//! first variable is projected first, the last one is lifted last). Its
//! projection is Collins' set: for every polynomial all coefficients in the
//! eliminated variable (which also covers the nullification of leading
//! coefficients), the resultant with its derivative (the discriminant up to
//! a factor of the leading coefficient), and the pairwise resultants. The
//! set is projected recursively down to one variable.
//!
//! Lifting evaluates the projection polynomials at an exact algebraic
//! sample point of the lower-dimensional cell, isolates the real roots of
//! the resulting one-variable polynomial over the sample point's number
//! field (see [`super::qe_field`]), and builds one section cell on every root
//! and one sector cell between consecutive roots (and outside the first and
//! last one). Sectors are sampled at rational points, which keeps the field
//! degree as small as possible. Every sample point is exact: its coordinates
//! live in a single number field `Q(theta)` and are represented as
//! polynomials in `theta`.
//!
//! Nullification is handled by skipping a projection polynomial whose value
//! at the sample point is identically zero: by sign invariance of the
//! projection, it then vanishes on the whole cylinder over the cell, so it
//! imposes no section there.
//!
//! All growth is bounded ([`CadLimits`]): the number of projection
//! polynomials per level and over the whole tower, the number of cells, the
//! degree of a fibre polynomial, the size of a number field and the number of
//! projection/lifting steps. There is no fixed variable limit; a decomposition
//! has one level per variable and is truncated only when a resource budget is
//! exhausted. Exhausting a bound returns `None`, which the callers turn into
//! an unevaluated request, never into a wrong answer.

#![allow(
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use std::cell::Cell as Counter;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Number;
use crate::rules::poly::multifactor;
use crate::rules::poly::repr::Poly as ReprPoly;

use super::qe_field::Elem;
use super::qe_field::Field;
use super::qe_field::f_degree;
use super::qe_field::f_divrem;
use super::qe_field::f_eval_rational;
use super::qe_field::f_gcd;
use super::qe_field::f_isolate;
use super::qe_field::f_mul;
use super::qe_field::f_square_free;
use super::qe_poly::MPoly;

/// Resource limits of a decomposition.
///
/// The number of variables is not limited: a decomposition has one level per
/// variable, and its size is controlled by the projection, cell and step
/// budgets instead. Every budget counts work or objects, so a formula with
/// many variables but a genuinely small decomposition is decided, while a
/// large one exhausts a budget and leaves the request unevaluated.
#[derive(Clone, Debug)]
pub(crate) struct CadLimits {
    /// Largest number of projection polynomials per level.
    pub projection: usize,
    /// Largest number of projection polynomials over the whole tower.
    pub projection_total: usize,
    /// Largest number of cells (leaves plus intermediate cells).
    pub cells: usize,
    /// Largest degree of a one-variable fibre polynomial.
    pub degree: u32,
    /// Largest number of terms of a projection polynomial.
    pub terms: usize,
    /// Largest degree of a number field.
    pub field_degree: usize,
    /// Largest number of projection and lifting steps in one request.
    pub steps: usize,
    /// The steps already spent. Reset only by constructing a new value, so
    /// the budget is shared by every block of one request.
    spent: Counter<usize>,
}

impl Default for CadLimits {
    fn default() -> Self {
        Self {
            projection: 256,
            projection_total: 4_096,
            cells: 20_000,
            degree: 48,
            terms: 20_000,
            field_degree: 64,
            steps: 1_000_000,
            spent: Counter::new(0),
        }
    }
}

impl CadLimits {
    /// Charges one projection or lifting step, returning `None` when the
    /// step budget is exhausted.
    fn spend(&self) -> Option<()> {
        let spent = self.spent.get().saturating_add(1);
        self.spent.set(spent);
        (spent <= self.steps).then_some(())
    }
}

/// The projection of `polys` onto the variables below `variable`.
///
/// `polys` are polynomials in variables `0..=variable`; the result is a set
/// of polynomials in the remaining variables. Constants and zero
/// polynomials are dropped, signs are normalised, and duplicates are
/// removed.
#[must_use]
pub(crate) fn project(
    polys: &[MPoly],
    variable: u32,
    limits: &CadLimits,
) -> Option<Vec<MPoly>> {
    let mut out: Vec<MPoly> = Vec::new();
    let mut coefficient_lists: Vec<Vec<MPoly>> = Vec::with_capacity(polys.len());
    for p in polys {
        limits.spend()?;
        let mut coefficients = p.coefficients_in(variable);
        while coefficients.last().is_some_and(MPoly::is_zero) {
            coefficients.pop();
        }
        coefficient_lists.push(coefficients);
    }
    for (i, p) in polys.iter().enumerate() {
        let Some(coefficients) = coefficient_lists.get(i) else {
            continue;
        };
        let Some(degree) = coefficients.len().checked_sub(1) else {
            continue;
        };
        if degree == 0 {
            if let Some(only) = coefficients.first() {
                push_projection(&mut out, only, limits)?;
            }
            continue;
        }
        for c in coefficients {
            push_projection(&mut out, c, limits)?;
        }
        let derivative = p.derivative(variable);
        let with_derivative = resultant(p, &derivative, variable, limits)?;
        push_projection(&mut out, &with_derivative, limits)?;
        for j in i.saturating_add(1)..polys.len() {
            if coefficient_lists.get(j).is_none_or(|c| c.len() < 2) {
                continue;
            }
            limits.spend()?;
            let Some(q) = polys.get(j) else {
                continue;
            };
            let value = resultant(p, q, variable, limits)?;
            push_projection(&mut out, &value, limits)?;
        }
    }
    // The projection theorem is stated for square-free polynomials, so the
    // whole set is replaced by its irreducible factors (multiplicities
    // dropped). This is what makes repeated factors in the input, and
    // coefficients or resultants with repeated factors, safe.
    let mut factored: Vec<MPoly> = Vec::new();
    for p in out {
        for q in squarefree_factors(&p)? {
            let q = q.positive_leading();
            if !factored.contains(&q) {
                factored.push(q);
            }
            if factored.len() > limits.projection {
                return None;
            }
        }
    }
    Some(factored)
}

/// The distinct irreducible factors of `p` over `Q`, with multiplicity
/// dropped. `None` when the factorisation is beyond its own limits (the
/// caller then leaves the request unevaluated instead of using a projection
/// that is not known to be square-free).
#[must_use]
fn squarefree_factors(p: &MPoly) -> Option<Vec<MPoly>> {
    let mut repr = ReprPoly::zero();
    for (mono, coeff) in p.terms() {
        repr = repr.add(&ReprPoly::monomial(
            mono.clone(),
            Number::rat(coeff.clone()),
        ));
    }
    let support = repr.support();
    let (_, factors) = multifactor::factor(&repr, &support)?;
    let mut out: Vec<MPoly> = Vec::new();
    for (factor, _) in factors {
        let mut q = MPoly::zero();
        for (mono, coeff) in factor.terms() {
            q = q.add(&MPoly::monomial(mono.clone(), coeff.to_rational()?));
        }
        if q.as_constant().is_some() || q.is_zero() {
            continue;
        }
        let q = q.positive_leading();
        if !out.contains(&q) {
            out.push(q);
        }
    }
    Some(if out.is_empty() {
        vec![p.clone()]
    } else {
        out
    })
}

/// Adds one polynomial to a projection set: drops constants and zero,
/// normalises the sign, removes duplicates and enforces the size limit.
fn push_projection(
    out: &mut Vec<MPoly>,
    p: &MPoly,
    limits: &CadLimits,
) -> Option<()> {
    limits.spend()?;
    if p.as_constant().is_some() {
        return Some(());
    }
    let p = p.positive_leading();
    if !out.contains(&p) {
        out.push(p);
    }
    (out.len() <= limits.projection).then_some(())
}

/// The resultant of `p` and `q` with respect to `variable`, by the
/// fraction-free Bareiss determinant of the Sylvester matrix.
#[must_use]
fn resultant(
    p: &MPoly,
    q: &MPoly,
    variable: u32,
    limits: &CadLimits,
) -> Option<MPoly> {
    let mut p_coeffs = p.coefficients_in(variable);
    let mut q_coeffs = q.coefficients_in(variable);
    while p_coeffs.last().is_some_and(MPoly::is_zero) {
        p_coeffs.pop();
    }
    while q_coeffs.last().is_some_and(MPoly::is_zero) {
        q_coeffs.pop();
    }
    let m = p_coeffs.len().checked_sub(1)?;
    let n = q_coeffs.len().checked_sub(1)?;
    if m == 0 && n == 0 {
        return Some(MPoly::constant(BigRational::one()));
    }
    let size = m.checked_add(n)?;
    let mut rows: Vec<Vec<MPoly>> = Vec::with_capacity(size);
    for shift in 0..n {
        let mut row = vec![MPoly::zero(); size];
        for (k, c) in p_coeffs.iter().rev().enumerate() {
            if let Some(slot) = row.get_mut(shift.checked_add(k)?) {
                *slot = c.clone();
            }
        }
        rows.push(row);
    }
    for shift in 0..m {
        let mut row = vec![MPoly::zero(); size];
        for (k, c) in q_coeffs.iter().rev().enumerate() {
            if let Some(slot) = row.get_mut(shift.checked_add(k)?) {
                *slot = c.clone();
            }
        }
        rows.push(row);
    }
    bareiss(rows, limits)
}

/// Determinant of a square matrix of polynomials by Bareiss' fraction-free
/// elimination.
fn bareiss(
    mut matrix: Vec<Vec<MPoly>>,
    limits: &CadLimits,
) -> Option<MPoly> {
    let n = matrix.len();
    if n == 0 {
        return Some(MPoly::constant(BigRational::one()));
    }
    let mut sign = false;
    let mut previous = MPoly::constant(BigRational::one());
    for k in 0..n.saturating_sub(1) {
        let pivot_zero = matrix.get(k)?.get(k)?.is_zero();
        if pivot_zero {
            let swap = (k.saturating_add(1)..n).find(|&r| {
                matrix
                    .get(r)
                    .and_then(|row| row.get(k))
                    .is_some_and(|v| !v.is_zero())
            });
            match swap {
                | Some(r) => {
                    matrix.swap(k, r);
                    sign = !sign;
                },
                | None => return Some(MPoly::zero()),
            }
        }
        let pivot = matrix.get(k)?.get(k)?.clone();
        for i in k.saturating_add(1)..n {
            for j in k.saturating_add(1)..n {
                let a = matrix.get(i)?.get(j)?.mul(&pivot, limits.terms)?;
                let b = matrix
                    .get(i)?
                    .get(k)?
                    .mul(&matrix.get(k)?.get(j)?.clone(), limits.terms)?;
                let value = a.sub(&b).div_exact(&previous, limits.terms)?;
                let row = matrix.get_mut(i)?;
                *row.get_mut(j)? = value;
            }
            if let Some(slot) = matrix.get_mut(i).and_then(|row| row.get_mut(k)) {
                *slot = MPoly::zero();
            }
        }
        previous = matrix.get(k)?.get(k)?.clone();
        if previous.is_zero() {
            return Some(MPoly::zero());
        }
    }
    let det = matrix
        .get(n.saturating_sub(1))?
        .get(n.saturating_sub(1))?
        .clone();
    Some(if sign { det.neg() } else { det })
}

/// Applies a polynomial at an exact sample point.
#[must_use]
pub(crate) fn evaluate(
    p: &MPoly,
    coords: &[Elem],
    field: &Field,
) -> Option<Elem> {
    let mut total = Elem::zero();
    for (mono, coeff) in p.terms() {
        let mut term = Elem::rational(coeff.clone());
        for &(variable, exponent) in mono {
            let value = coords.get(usize::try_from(variable).ok()?)?;
            let power = field.pow(value, exponent);
            term = field.mul(&term, &power);
        }
        total = total.add(&term);
    }
    Some(total)
}

/// The signs of `polys` at a sample point.
#[must_use]
pub(crate) fn signs_at(
    polys: &[MPoly],
    coords: &[Elem],
    field: &mut Field,
) -> Option<Vec<i32>> {
    let mut out = Vec::with_capacity(polys.len());
    for p in polys {
        let value = evaluate(p, coords, field)?;
        out.push(field.sign(&value)?);
    }
    Some(out)
}

/// The sample point and cell structure of one cell.
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    /// Coordinates, each an element of `field`.
    pub sample: Vec<Elem>,
    /// `open[i]` is true when the cell is a sector in coordinate `i`.
    pub open: Vec<bool>,
    /// Signs of the input polynomials at the sample point.
    pub signs: Vec<i32>,
    /// The number field the coordinates live in.
    pub field: Field,
}

/// A parameter-space cell with the aggregated truth of the quantified
/// block and the signs of the parameter projection polynomials.
#[derive(Clone, Debug)]
pub(crate) struct ParamCell {
    /// Aggregated truth of the quantified block over the fibre.
    pub truth: bool,
    /// Signs of the parameter projection polynomials at the sample.
    pub signs: Vec<i32>,
}

/// A context carried down the lifting recursion.
#[derive(Clone)]
struct Ctx {
    field: Field,
    coords: Vec<Elem>,
    open: Vec<bool>,
}

/// A decomposition together with its parameter information.
pub(crate) struct Lift<'a> {
    /// The input polynomials (used for the cell signs).
    pub polys: &'a [MPoly],
    /// `projection[k]` holds polynomials in variables `0..k`.
    pub projection: &'a [Vec<MPoly>],
    /// The projection polynomials in parameter variables only, which
    /// describe the parameter cells.
    pub parameter_polys: &'a [MPoly],
    /// Number of variables.
    pub variables: usize,
    /// Number of leading variables treated as parameters.
    pub params: usize,
    /// Whether the quantified block is existential.
    pub exists: bool,
    /// Truth of the formula for a vector of atom signs.
    pub truth: &'a dyn Fn(&[i32]) -> bool,
    /// Limits.
    pub limits: &'a CadLimits,
    /// Collect every leaf cell (for `cad_cells`).
    pub collect: bool,
    /// The leaf cells, when collecting.
    pub leaves: Vec<Cell>,
    /// The parameter cells.
    pub params_out: Vec<ParamCell>,
    /// Cell counter.
    work: usize,
}

impl Lift<'_> {
    /// Builds the whole decomposition, returning the parameter cells with
    /// their aggregated truth.
    fn build(&mut self) -> Option<()> {
        let field = Field::rational_field()?;
        let ctx = Ctx {
            field,
            coords: Vec::new(),
            open: Vec::new(),
        };
        self.run(0, &ctx)?;
        Some(())
    }

    /// One level of the recursion. Levels below `params` build parameter
    /// cells; the level equal to `params` aggregates the quantified levels.
    fn run(
        &mut self,
        level: usize,
        ctx: &Ctx,
    ) -> Option<bool> {
        if level == self.params {
            let mut accumulator = !self.exists;
            self.fibre(level, ctx.clone(), &mut accumulator)?;
            let mut field = ctx.field.clone();
            let signs = if self.params == 0 {
                Vec::new()
            } else {
                signs_at(self.parameter_polys, &ctx.coords, &mut field)?
            };
            self.params_out.push(ParamCell {
                truth: accumulator,
                signs,
            });
            return Some(accumulator);
        }
        let children = self.fibre_contexts(level, ctx)?;
        for child in children {
            self.run(level.saturating_add(1), &child)?;
            if self.work > self.limits.cells {
                return None;
            }
        }
        Some(true)
    }

    /// The levels at and above `params`: aggregates the truth of the leaves.
    fn fibre(
        &mut self,
        level: usize,
        ctx: Ctx,
        accumulator: &mut bool,
    ) -> Option<()> {
        if level == self.variables {
            self.limits.spend()?;
            self.work = self.work.saturating_add(1);
            if self.work > self.limits.cells {
                return None;
            }
            let mut field = ctx.field.clone();
            let signs = signs_at(self.polys, &ctx.coords, &mut field)?;
            let value = (self.truth)(&signs);
            if self.collect {
                self.leaves.push(Cell {
                    sample: ctx.coords,
                    open: ctx.open,
                    signs,
                    field: ctx.field,
                });
            } else {
                *accumulator = if self.exists {
                    *accumulator || value
                } else {
                    *accumulator && value
                };
            }
            return Some(());
        }
        let children = self.fibre_contexts(level, &ctx)?;
        for child in children {
            self.fibre(level.saturating_add(1), child, accumulator)?;
        }
        Some(())
    }

    /// The section and sector cells one level above `ctx`.
    fn fibre_contexts(
        &self,
        level: usize,
        ctx: &Ctx,
    ) -> Option<Vec<Ctx>> {
        let projection = self.projection.get(level.checked_add(1)?)?;
        let mut work = ctx.field.clone();
        let mut fibres: Vec<Vec<Elem>> = Vec::new();
        for p in projection {
            self.limits.spend()?;
            let coefficients = p.coefficients_in(u32::try_from(level).ok()?);
            let mut fibre = Vec::with_capacity(coefficients.len());
            for c in &coefficients {
                fibre.push(evaluate(c, &ctx.coords, &work)?);
            }
            while fibre.last().is_some_and(Elem::is_zero) {
                fibre.pop();
            }
            let Some(degree) = f_degree(&fibre) else {
                continue;
            };
            if degree > usize::try_from(self.limits.degree).ok()? {
                return None;
            }
            let square_free = f_square_free(&fibre, &work);
            if f_degree(&square_free).is_some() {
                fibres.push(square_free);
            }
        }
        let mut parts: Vec<Vec<Elem>> = Vec::new();
        let mut product = vec![Elem::one()];
        for f in &fibres {
            let common = f_gcd(&product, f, &work);
            let (part, _) = f_divrem(f, &common, &work)?;
            if f_degree(&part).is_none() {
                continue;
            }
            product = f_mul(&product, &part, &work);
            let degree = f_degree(&product)?;
            if degree > usize::try_from(self.limits.degree).ok()? {
                return None;
            }
            parts.push(part);
        }
        let roots = f_isolate(&product, &mut work)?;
        if roots.len() > self.limits.cells {
            return None;
        }
        let mut out: Vec<Ctx> = Vec::with_capacity(roots.len().saturating_mul(2).saturating_add(1));
        let first = roots.first();
        let initial = match first {
            | Some((lo, _)) => lo - BigRational::one(),
            | None => BigRational::zero(),
        };
        out.push(sector_context(ctx, &work, initial));
        for (index, (lo, hi)) in roots.iter().enumerate() {
            self.limits.spend()?;
            // A point interval is a rational root: use it directly.
            if lo == hi {
                let mut coords = ctx.coords.clone();
                coords.push(Elem::rational(lo.clone()));
                let mut open = ctx.open.clone();
                open.push(false);
                out.push(Ctx {
                    field: work.clone(),
                    coords,
                    open,
                });
            } else {
                // Adjoin the root through the single coprime factor that
                // contains it, not through the product of every factor. The
                // product is generally reducible over the current field, and
                // extending by a reducible polynomial loses the field property
                // the primitive-element construction relies on.
                let factor = parts
                    .iter()
                    .find(|part| interval_contains_root(part, lo, hi, &mut work) == Some(true))
                    .unwrap_or(&product);
                // A linear factor already has its root in the current field;
                // use it directly instead of adjoining an isomorphic copy.
                let linear = if f_degree(factor) == Some(1) {
                    let constant = factor.first()?.neg();
                    let leading = factor.get(1)?.clone();
                    Some(work.div(&constant, &leading)?)
                } else {
                    None
                };
                if let Some(root) = linear {
                    let mut coords = ctx.coords.clone();
                    coords.push(root);
                    let mut open = ctx.open.clone();
                    open.push(false);
                    out.push(Ctx {
                        field: work.clone(),
                        coords,
                        open,
                    });
                } else {
                    let extension =
                        super::qe_field::extend(&mut work, factor, &(lo.clone(), hi.clone()))?;
                    if extension.field.degree() > self.limits.field_degree {
                        return None;
                    }
                    let mut coords: Vec<Elem> =
                        Vec::with_capacity(ctx.coords.len().saturating_add(1));
                    for coordinate in &ctx.coords {
                        coords.push(Field::substitute(
                            coordinate,
                            &extension.theta,
                            &extension.field,
                        ));
                    }
                    coords.push(extension.root);
                    let mut open = ctx.open.clone();
                    open.push(false);
                    out.push(Ctx {
                        field: extension.field,
                        coords,
                        open,
                    });
                }
            }
            let next_lo = roots.get(index.saturating_add(1)).map(|(a, _)| a.clone());
            let sample = match next_lo {
                | Some(a) => {
                    let two = BigRational::from_integer(BigInt::from(2));
                    (hi + a) / &two
                },
                | None => hi + BigRational::one(),
            };
            out.push(sector_context(ctx, &work, sample));
        }
        Some(out)
    }
}

/// Whether the square-free `g` has the root isolated by `[lo, hi]` (its sign
/// changes across the interval). Endpoint roots are excluded: such an endpoint
/// belongs to a different root that is reported as its own point interval.
fn interval_contains_root(
    g: &[Elem],
    lo: &BigRational,
    hi: &BigRational,
    field: &mut Field,
) -> Option<bool> {
    let at_lo = f_eval_rational(g, lo, field);
    let at_hi = f_eval_rational(g, hi, field);
    let lo_sign = field.sign(&at_lo)?;
    let hi_sign = field.sign(&at_hi)?;
    Some(lo_sign != 0 && hi_sign != 0 && lo_sign != hi_sign)
}

/// The context of a sector cell: the lower coordinates plus one rational
/// sample.
fn sector_context(
    ctx: &Ctx,
    field: &Field,
    value: BigRational,
) -> Ctx {
    let mut coords = ctx.coords.clone();
    coords.push(Elem::rational(value));
    let mut open = ctx.open.clone();
    open.push(true);
    Ctx {
        field: field.clone(),
        coords,
        open,
    }
}

/// The projection sets `P[k]` (polynomials in variables `0..k`) used by the
/// lifting.
#[must_use]
pub(crate) fn projections(
    polys: &[MPoly],
    variables: usize,
    limits: &CadLimits,
) -> Option<Vec<Vec<MPoly>>> {
    if variables == 0 {
        return None;
    }
    let mut projection = vec![Vec::new(); variables.saturating_add(1)];
    let mut total = 0_usize;
    if let Some(slot) = projection.get_mut(variables) {
        let mut factored: Vec<MPoly> = Vec::new();
        for p in polys {
            limits.spend()?;
            for q in squarefree_factors(p)? {
                let q = q.positive_leading();
                if !factored.contains(&q) {
                    factored.push(q);
                }
                if factored.len() > limits.projection {
                    return None;
                }
            }
        }
        total = total.saturating_add(factored.len());
        if total > limits.projection_total {
            return None;
        }
        *slot = factored;
    }
    for k in (1..variables).rev() {
        limits.spend()?;
        let higher = projection.get(k.checked_add(1)?)?.clone();
        let projected = project(&higher, u32::try_from(k).ok()?, limits)?;
        total = total.saturating_add(projected.len());
        if total > limits.projection_total {
            return None;
        }
        if let Some(slot) = projection.get_mut(k) {
            *slot = projected;
        }
    }
    Some(projection)
}

/// The cells of a sign-invariant decomposition of `R^n` for `polys`.
#[must_use]
pub(crate) fn cells(
    polys: &[MPoly],
    variables: usize,
    limits: &CadLimits,
) -> Option<Vec<Cell>> {
    let projection = projections(polys, variables, limits)?;
    let truth = |_: &[i32]| true;
    let mut lift = Lift {
        polys,
        projection: &projection,
        parameter_polys: &[],
        variables,
        params: variables,
        exists: true,
        truth: &truth,
        limits,
        collect: true,
        leaves: Vec::new(),
        params_out: Vec::new(),
        work: 0,
    };
    let field = Field::rational_field()?;
    let ctx = Ctx {
        field,
        coords: Vec::new(),
        open: Vec::new(),
    };
    lift.fibre(0, ctx, &mut false)?;
    Some(lift.leaves)
}

/// QE: aggregate the truth of the block over each parameter cell, together
/// with the parameter projection polynomials used to describe the cells.
#[must_use]
pub(crate) fn quantify(
    polys: &[MPoly],
    variables: usize,
    params: usize,
    exists: bool,
    truth: &dyn Fn(&[i32]) -> bool,
    limits: &CadLimits,
) -> Option<(Vec<ParamCell>, Vec<MPoly>)> {
    let projection = projections(polys, variables, limits)?;
    // The cells of the parameter space are determined by the whole
    // projection tower below `params`, not only by the top level: every
    // polynomial there involves parameter variables only.
    let mut parameter_polys: Vec<MPoly> = Vec::new();
    for level in 1..=params {
        let Some(set) = projection.get(level) else {
            continue;
        };
        for p in set {
            if p.as_constant().is_none() && !p.is_zero() {
                let p = p.positive_leading();
                if !parameter_polys.contains(&p) {
                    parameter_polys.push(p);
                }
            }
        }
    }
    if parameter_polys.len() > limits.projection {
        return None;
    }
    let mut lift = Lift {
        polys,
        projection: &projection,
        parameter_polys: &parameter_polys,
        variables,
        params,
        exists,
        truth,
        limits,
        collect: false,
        leaves: Vec::new(),
        params_out: Vec::new(),
        work: 0,
    };
    lift.build()?;
    Some((lift.params_out, parameter_polys))
}

#[cfg(test)]
mod tests {
    use num_traits::Signed;

    use super::*;

    fn integer(n: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(n))
    }

    fn x() -> MPoly {
        MPoly::variable(0)
    }

    fn y() -> MPoly {
        MPoly::variable(1)
    }

    fn x_squared_minus_one() -> MPoly {
        x().pow(2, 100)
            .unwrap_or_default()
            .sub(&MPoly::constant(integer(1)))
    }

    fn unit_circle() -> MPoly {
        x().pow(2, 100)
            .unwrap_or_default()
            .add(&y().pow(2, 100).unwrap_or_default())
            .sub(&MPoly::constant(integer(1)))
    }

    /// `x1^2 + ... + x5^2 - 1`.
    fn unit_sphere_five() -> MPoly {
        let mut sum = MPoly::zero();
        for variable in 0..5 {
            sum = sum.add(&MPoly::variable(variable).pow(2, 100).unwrap_or_default());
        }
        sum.sub(&MPoly::constant(integer(1)))
    }

    #[test]
    fn projection_of_the_unit_circle() {
        let polys = vec![unit_circle()];
        let limits = CadLimits::default();
        let projected = project(&polys, 1, &limits).unwrap_or_default();
        assert!(!projected.is_empty());
        assert!(projected.iter().any(|p| p.as_constant().is_none()));
        // The sections of the circle occur at x = +/-1 (the projection is
        // split into the two linear factors x - 1 and x + 1).
        let mut has_negative = false;
        let mut has_positive = false;
        for q in projected.iter().filter_map(|p| p.univariate_in(0)) {
            for (lo, hi) in crate::rules::solve::qe_alg::isolate_roots(&q) {
                let middle = (lo + hi) / BigRational::from_integer(BigInt::from(2));
                if middle.is_negative() {
                    has_negative = true;
                } else if middle.is_positive() {
                    has_positive = true;
                }
            }
        }
        assert!(has_negative && has_positive);
    }

    #[test]
    fn line_cells() {
        let polys = vec![x_squared_minus_one()];
        let limits = CadLimits::default();
        let cells = cells(&polys, 1, &limits).unwrap_or_default();
        assert_eq!(cells.len(), 5);
        let signs: Vec<i32> = cells
            .iter()
            .filter_map(|c| c.signs.first().copied())
            .collect();
        assert_eq!(signs, vec![1, 0, -1, 0, 1]);
        for cell in &cells {
            let open = cell.open.first().copied().unwrap_or(false);
            assert_eq!(open, cell.signs.first().copied().unwrap_or(0) != 0);
        }
    }

    #[test]
    fn circle_cells_are_open_or_sections() {
        let polys = vec![unit_circle()];
        let limits = CadLimits::default();
        let cells = cells(&polys, 2, &limits).unwrap_or_default();
        assert!(cells.len() > 9, "{}", cells.len());
        assert!(
            cells
                .iter()
                .all(|c| c.sample.len() == 2 && c.open.len() == 2)
        );
        let signs: Vec<i32> = cells
            .iter()
            .filter_map(|c| c.signs.first().copied())
            .collect();
        assert!(signs.contains(&-1) && signs.contains(&0) && signs.contains(&1));
    }

    #[test]
    fn tesseract_limits() {
        let polys = vec![unit_circle()];
        let mut limits = CadLimits::default();
        limits.cells = 2;
        assert!(cells(&polys, 2, &limits).is_none());
    }

    #[test]
    fn coprime_fibre_factors_are_extended_separately() {
        // At the parameter sections p0 = +/-sqrt(2) the two fibre factors
        // p1^2 - 3 and p1^2 - 5 are coprime and irreducible over Q(sqrt(2)).
        // Extending by their product gives a reducible algebra, which the
        // primitive-element construction cannot invert; each factor has to be
        // extended on its own.
        let p0 = x()
            .pow(2, 100)
            .unwrap_or_default()
            .sub(&MPoly::constant(integer(2)));
        let p1_three = y()
            .pow(2, 100)
            .unwrap_or_default()
            .sub(&MPoly::constant(integer(3)));
        let p1_five = y()
            .pow(2, 100)
            .unwrap_or_default()
            .sub(&MPoly::constant(integer(5)));
        let projection = vec![
            Vec::new(),
            vec![p0],
            vec![p1_three.clone(), p1_five.clone()],
        ];
        let polys = vec![p1_three, p1_five];
        let truth = |_: &[i32]| true;
        let limits = CadLimits::default();
        let mut lift = Lift {
            polys: &polys,
            projection: &projection,
            parameter_polys: &[],
            variables: 2,
            params: 2,
            exists: true,
            truth: &truth,
            limits: &limits,
            collect: true,
            leaves: Vec::new(),
            params_out: Vec::new(),
            work: 0,
        };
        let ctx = Ctx {
            field: Field::rational_field().unwrap_or_else(|| panic!()),
            coords: Vec::new(),
            open: Vec::new(),
        };
        lift.fibre(0, ctx, &mut false).expect("lifting");
        // Five p0 cells and nine p1 cells over each of them.
        assert_eq!(lift.leaves.len(), 45);
        for cell in &lift.leaves {
            assert_eq!(cell.sample.len(), 2);
            assert_eq!(cell.signs.len(), 2);
            let mut field = cell.field.clone();
            let signs = signs_at(&polys, &cell.sample, &mut field).expect("signs");
            assert_eq!(signs, cell.signs);
        }
    }

    #[test]
    fn five_variable_sphere_cells() {
        let polys = vec![unit_sphere_five()];
        let limits = CadLimits::default();
        let cells = cells(&polys, 5, &limits).expect("five-variable decomposition");
        assert!(cells.len() > 30, "{}", cells.len());
        let mut sections = 0_usize;
        for cell in &cells {
            assert_eq!(cell.sample.len(), 5);
            assert_eq!(cell.open.len(), 5);
            assert_eq!(cell.signs.len(), 1);
            // The stored sign vector is reproduced from the exact sample.
            let mut field = cell.field.clone();
            let signs = signs_at(&polys, &cell.sample, &mut field).expect("signs");
            assert_eq!(signs, cell.signs);
            if cell.open.iter().any(|&open| !open) {
                sections = sections.saturating_add(1);
            }
        }
        assert!(sections > 0, "the sphere must have section cells");
        assert!(
            cells.iter().any(|cell| cell.signs.first() == Some(&0)),
            "the sphere itself must occur as a cell"
        );
    }
}
