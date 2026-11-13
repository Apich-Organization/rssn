//! Why an operator did not reduce.
//!
//! A numeric operator whose arguments are rejected, whose system is singular or
//! whose data are incompatible does not fail loudly: the engine leaves the call
//! as an unreduced term, which is the right answer for a rewrite system (a
//! request is never a wrong value) and the wrong answer for a person who wants
//! to know what to fix. This module is the explanation channel.
//!
//! # The pieces
//!
//! * [`Diag`]: one finding: the operator, the zero-based position of the
//!   offending argument when it is known, a [`DiagKind`] to branch on and a
//!   human message.
//! * [`fail`] / [`report`]: how an evaluator states a reason. `fail` returns
//!   `None`, so an `Option`-returning evaluator reads
//!   `return diag::fail(DiagKind::Singular, None, format_args!("..."))`, and
//!   every existing `Option` evaluator keeps compiling and behaving exactly as
//!   before: reporting changes no result, it only leaves a note.
//! * [`capture`]: runs a closure with a collector installed on the current
//!   thread and returns the first finding reported inside it. Outside a
//!   capture [`report`] does nothing (the message is not even formatted), so
//!   normal reduction pays nothing for the channel.
//! * [`Session::apply_diag`](crate::api::Session::apply_diag) and
//!   [`Session::explain`](crate::api::Session::explain): the host entry
//!   points. The first answers "apply operator X to these arguments and tell me
//!   why it did not reduce"; the second does the same for a term that is
//!   already built.
//!
//! The collector is thread-local and scoped. It is installed by the caller of
//! the evaluator on the caller's own thread, never by the engine, so it does
//! not depend on which worker thread a reduction ran on and it cannot leak a
//! finding into an unrelated call.

use std::cell::RefCell;
use std::fmt;

use serde::Deserialize;
use serde::Serialize;

/// What kind of problem a [`Diag`] reports.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
#[non_exhaustive]
pub enum DiagKind {
    /// Wrong number of arguments.
    Arity,
    /// An argument has the wrong shape: a list where a number is wanted, the
    /// wrong rank, extents that disagree, a ragged list.
    Shape,
    /// An argument has the right shape and a value the operator rejects: out
    /// of range, not an integer, an unknown name, a non-finite number, an
    /// unstable step.
    Domain,
    /// The problem is singular (a singular matrix, a non-invertible operator).
    Singular,
    /// The data contradict each other: a system with no solution, a
    /// compatibility condition that fails.
    Incompatible,
    /// An iteration ran out of budget before reaching its tolerance.
    NotConverged,
    /// A well-posed request outside what the operator supports.
    Unsupported,
    /// An argument is not a literal number, list or symbol (it is a symbolic
    /// expression), and the operator computes only from literals.
    Symbolic,
    /// The call did not reduce and the operator gave no reason.
    Unspecified,
}

impl DiagKind {
    /// A stable machine-readable code, `rssn.diag.<kind>`.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            | Self::Arity => "rssn.diag.arity",
            | Self::Shape => "rssn.diag.shape",
            | Self::Domain => "rssn.diag.domain",
            | Self::Singular => "rssn.diag.singular",
            | Self::Incompatible => "rssn.diag.incompatible",
            | Self::NotConverged => "rssn.diag.not_converged",
            | Self::Unsupported => "rssn.diag.unsupported",
            | Self::Symbolic => "rssn.diag.symbolic",
            | Self::Unspecified => "rssn.diag.unspecified",
        }
    }
}

impl fmt::Display for DiagKind {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(match self {
            | Self::Arity => "arity",
            | Self::Shape => "shape",
            | Self::Domain => "domain",
            | Self::Singular => "singular",
            | Self::Incompatible => "incompatible",
            | Self::NotConverged => "not converged",
            | Self::Unsupported => "unsupported",
            | Self::Symbolic => "symbolic",
            | Self::Unspecified => "unspecified",
        })
    }
}

/// One finding: why an operator call did not reduce.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Diag {
    /// The operator's name; empty while the finding is still inside an
    /// evaluator, filled in by whoever dispatched it.
    pub op: String,
    /// Zero-based position of the argument at fault, when one is.
    pub arg: Option<usize>,
    /// What kind of problem it is.
    pub kind: DiagKind,
    /// A sentence a person can act on.
    pub message: String,
}

impl Diag {
    /// A finding with no operator and no argument yet.
    #[must_use]
    pub fn new(
        kind: DiagKind,
        message: impl Into<String>,
    ) -> Self {
        Self {
            op: String::new(),
            arg: None,
            kind,
            message: message.into(),
        }
    }

    /// The same finding blamed on argument `index` (zero-based).
    #[must_use]
    pub fn at(
        mut self,
        index: usize,
    ) -> Self {
        self.arg = Some(index);
        self
    }

    /// The same finding attributed to operator `op`.
    #[must_use]
    pub fn in_op(
        mut self,
        op: impl Into<String>,
    ) -> Self {
        self.op = op.into();
        self
    }
}

impl fmt::Display for Diag {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        if !self.op.is_empty() {
            write!(f, "{}: ", self.op)?;
        }
        if let Some(arg) = self.arg {
            write!(f, "argument {} ", arg + 1)?;
        }
        write!(f, "({}) {}", self.kind, self.message)
    }
}

impl std::error::Error for Diag {}

impl From<Diag> for crate::kernel::call::OpError {
    /// A host diagnostic carrying the kind's [stable code](DiagKind::code) and
    /// the full text of the finding.
    fn from(diag: Diag) -> Self {
        Self::HostDiagnostic {
            code: diag.kind.code(),
            message: diag.to_string(),
        }
    }
}

thread_local! {
    static SINK: RefCell<Option<Vec<Diag>>> = const { RefCell::new(None) };
}

/// Restores the enclosing collector (or none) when a [`capture`] ends, even by
/// unwinding.
struct Restore(Option<Vec<Diag>>);

impl Drop for Restore {
    fn drop(&mut self) {
        let previous = self.0.take();
        SINK.with(|sink| {
            if let Ok(mut slot) = sink.try_borrow_mut() {
                *slot = previous;
            }
        });
    }
}

/// Runs `f` with a collector installed on this thread and returns its result
/// with every finding [reported](report) inside, in order.
///
/// Captures nest: an inner capture sees only its own findings.
pub fn capture_all<R>(f: impl FnOnce() -> R) -> (R, Vec<Diag>) {
    let previous = SINK.with(|sink| sink.borrow_mut().replace(Vec::new()));
    let restore = Restore(previous);
    let out = f();
    let found = SINK
        .with(|sink| sink.borrow_mut().take())
        .unwrap_or_default();
    drop(restore);
    (out, found)
}

/// [`capture_all`] keeping only the first finding, the root cause (an
/// evaluator stops at its first failure).
pub fn capture<R>(f: impl FnOnce() -> R) -> (R, Option<Diag>) {
    let (out, found) = capture_all(f);
    (out, found.into_iter().next())
}

/// Whether a collector is installed on this thread.
#[must_use]
pub fn active() -> bool {
    SINK.with(|sink| sink.try_borrow().is_ok_and(|slot| slot.is_some()))
}

/// Records a finding if a collector is installed; otherwise does nothing.
pub fn report(diag: Diag) {
    SINK.with(|sink| {
        if let Ok(mut slot) = sink.try_borrow_mut()
            && let Some(found) = slot.as_mut()
        {
            found.push(diag);
        }
    });
}

/// Records a finding and returns `None`, so an `Option` evaluator can write
/// `return diag::fail(kind, Some(arg), format_args!("..."))`.
///
/// The message is formatted only when a collector is installed.
#[must_use]
pub fn fail<T>(
    kind: DiagKind,
    arg: Option<usize>,
    message: fmt::Arguments<'_>,
) -> Option<T> {
    if active() {
        let mut diag = Diag::new(kind, message.to_string());
        diag.arg = arg;
        report(diag);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::Diag;
    use super::DiagKind;
    use super::active;
    use super::capture;
    use super::capture_all;
    use super::fail;
    use super::report;
    use crate::kernel::call::OpError;

    fn evaluator(x: f64) -> Option<f64> {
        if x < 0.0 {
            return fail(DiagKind::Domain, Some(0), format_args!("{x} is negative"));
        }
        Some(x.sqrt())
    }

    #[test]
    fn a_failing_evaluator_leaves_a_note_only_inside_a_capture() {
        // Outside a capture nothing is kept and nothing changes for the caller.
        assert_eq!(evaluator(-1.0), None);
        assert!(!active());
        let (value, diag) = capture(|| evaluator(-4.0));
        assert_eq!(value, None);
        let diag = diag.expect("a finding");
        assert_eq!(diag.kind, DiagKind::Domain);
        assert_eq!(diag.arg, Some(0));
        assert_eq!(diag.message, "-4 is negative");
        // A success reports nothing.
        let (value, diag) = capture(|| evaluator(9.0));
        assert_eq!(value, Some(3.0));
        assert!(diag.is_none());
    }

    #[test]
    fn captures_nest_and_restore_the_outer_collector() {
        let (_, outer) = capture_all(|| {
            report(Diag::new(DiagKind::Shape, "outer one"));
            let (_, inner) = capture_all(|| report(Diag::new(DiagKind::Domain, "inner")));
            assert_eq!(inner.len(), 1);
            report(Diag::new(DiagKind::Singular, "outer two"));
        });
        let kinds: Vec<_> = outer.iter().map(|d| d.kind).collect();
        assert_eq!(kinds, vec![DiagKind::Shape, DiagKind::Singular]);
        assert!(!active());
    }

    #[test]
    fn a_finding_prints_and_converts_to_an_op_error() {
        let diag = Diag::new(DiagKind::Incompatible, "fluxes do not balance")
            .at(1)
            .in_op("poisson_fd_bc");
        assert_eq!(
            diag.to_string(),
            "poisson_fd_bc: argument 2 (incompatible) fluxes do not balance"
        );
        match OpError::from(diag) {
            | OpError::HostDiagnostic { code, message } => {
                assert_eq!(code, "rssn.diag.incompatible");
                assert!(message.contains("fluxes do not balance"), "{message}");
            },
            | other => panic!("{other:?}"),
        }
    }
}
