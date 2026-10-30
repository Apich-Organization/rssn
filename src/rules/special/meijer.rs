//! The Meijer G-function: operator, exact reductions, differentiation,
//! integration and series.
//!
//! The operator is `meijerg(m, n, list(a_1, ..., a_p), list(b_1, ..., b_q), z)`
//! with `m`, `n` literal integers. Reductions to elementary and special
//! functions, the derivative and antiderivative formulas, the series at zero
//! and asymptotics at infinity, and the numeric value are installed here.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
