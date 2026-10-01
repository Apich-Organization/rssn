//! stub
use crate::graph::rule::Installer;
use crate::graph::RuleError;

pub(crate) fn install(_i: &mut Installer<'_>) -> Result<(), RuleError> {
    Ok(())
}
