//! Real quantifier elimination.
//!
//! Quantified formulas over the real numbers are decided and eliminated by
//! a cylindrical algebraic decomposition: a projection to fewer variables,
//! isolation of the real roots of the projection polynomials, lifting over
//! every cell, and a truth evaluation on a sample point of each cell. The
//! answer is a quantifier-free formula in the free parameters, built from
//! sign conditions and the defining polynomial of the algebraic sample
//! points (via `root_of`).
//!
//! Operators: `real_qe(formula, list(x, ...))` eliminates the listed
//! variables; `cad_cells(polys, list(x, ...))` returns the cells of a
//! sign-invariant decomposition as sample point, signs and openness.
//! The logic rule set's `exists` and `forall` are routed here for formulas
//! whose atoms are polynomial (in)equalities over the reals.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
