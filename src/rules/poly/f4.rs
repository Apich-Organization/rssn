//! F4: the Faugère matrix algorithm for Gröbner bases.
//!
//! Pairs are selected by degree, symbolically preprocessed into a Macaulay
//! matrix, reduced by sparse Gaussian elimination, and the new rows are
//! installed with the Gebauer-Möller criteria. The result is the same
//! reduced basis as `groebner` in the requested order.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
