//! F5: the Faugère signature-based algorithm for Gröbner bases.
//!
//! Every polynomial carries the signature of the S-pair it comes from; the
//! rewritten and rewritable criteria discard pairs whose signature is
//! already represented. Homogeneous input is handled directly, affine
//! input through the standard homogenisation by a fresh variable.

use crate::graph::RuleError;
use crate::graph::rule::Installer;

pub(super) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
