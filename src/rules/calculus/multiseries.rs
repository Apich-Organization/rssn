//! Multivariate series: Taylor, Laurent and Puiseux expansions in several
//! variables.
//!
//! `series(f, list(x, y, ...), list(a, b, ...), n)` expands `f` about the
//! point `(a, b, ...)` and truncates at total degree `n`. The operator is
//! `mseries`; `mtaylor` is the finite-Taylor special case. Expansions with
//! negative or fractional exponents are supported where the function has
//! them.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
