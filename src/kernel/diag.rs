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
    pub const fn at(
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

/// Classifies an unreduced call of a heavy operator that recorded no reason.
///
/// Table-driven numeric operators record their own reasons (see
/// `rules::numeric`); a **heavy** operator — one with symbolic rules, a typed
/// dispatch or a host callback — leaves its call as a term when it declines,
/// and the places that decline are too many to instrument one by one. The
/// operator's own metadata answers instead:
///
/// * the argument count is not one the descriptor accepts → [`DiagKind::Arity`];
/// * the descriptor has typed overloads and none of them accepts these
///   arguments → [`DiagKind::Shape`], naming the candidate signatures;
/// * several overloads tie → [`DiagKind::Shape`], naming them;
/// * an argument is a symbolic expression while the operator has scalar
///   numeric semantics → [`DiagKind::Symbolic`];
/// * otherwise → [`DiagKind::Unsupported`], with the descriptor's
///   [`DynamicReason`](crate::kernel::signature::DynamicReason) when it has one
///   (a request reduced by rules, a variadic shape, ...).
///
/// The answer is *derived*, not recorded, and the message says only what the
/// metadata supports: it never invents a numeric cause that was not observed.
#[must_use]
pub fn classify(
    descriptor: &crate::graph::op::OpDescriptor,
    args: &[crate::kernel::value::Value],
) -> Diag {
    let name: &str = &descriptor.name;
    let given = args.len();
    if !descriptor.arity.accepts(given) {
        let expected = descriptor
            .arity
            .max()
            .map_or_else(|| "any number of".to_string(), |count| count.to_string());
        return Diag::new(
            DiagKind::Arity,
            format!("takes {expected} arguments, got {given}"),
        )
        .in_op(name);
    }
    if !descriptor.overloads.is_empty() {
        let candidates = || {
            let shown: Vec<String> = descriptor
                .overloads
                .iter()
                .take(4)
                .map(render_signature)
                .collect();
            shown.join("; ")
        };
        match crate::kernel::overload::best(&descriptor.overloads, args) {
            | crate::kernel::overload::Best::None => {
                return Diag::new(
                    DiagKind::Shape,
                    format!(
                        "no overload accepts ({}); the operator declares {}",
                        render_kinds(args),
                        candidates()
                    ),
                )
                .in_op(name);
            },
            | crate::kernel::overload::Best::Ambiguous(tied) => {
                let shown: Vec<String> = tied
                    .iter()
                    .filter_map(|&index| descriptor.overloads.get(index))
                    .map(render_signature)
                    .collect();
                return Diag::new(
                    DiagKind::Shape,
                    format!(
                        "several overloads accept ({}) equally well: {}",
                        render_kinds(args),
                        shown.join("; ")
                    ),
                )
                .in_op(name);
            },
            | crate::kernel::overload::Best::Unique(_) => {},
        }
    }
    if descriptor.eval.is_some()
        && let Some(index) = args
            .iter()
            .position(|value| matches!(value, crate::kernel::value::Value::Term(_)))
    {
        return Diag::new(
            DiagKind::Symbolic,
            "the argument is a symbolic expression, and this operator computes only from literal \
             numbers",
        )
        .at(index)
        .in_op(name);
    }
    let reason = match descriptor.dynamic_reason {
        | Some(crate::kernel::signature::DynamicReason::Request) => {
            "the operator is reduced by rules, and no rule applies to this form"
        },
        | Some(crate::kernel::signature::DynamicReason::ShapeDependent) => {
            "the result shape depends on the arguments, and this shape is outside what the \
             operator supports"
        },
        | Some(crate::kernel::signature::DynamicReason::FullyPolymorphic) => {
            "the operator is polymorphic over its arguments, and no rule applies to this form"
        },
        | Some(crate::kernel::signature::DynamicReason::OpaquePayload) => {
            "the operator needs a host payload, which this call does not carry"
        },
        | Some(crate::kernel::signature::DynamicReason::Variadic) => {
            "the operator takes a variadic argument list, and no rule applies to this form"
        },
        | Some(crate::kernel::signature::DynamicReason::NotYetSpecified) => {
            "the operator has no typed signature yet, and its rules did not reduce this form"
        },
        | None => "the operator declined this form and recorded no reason",
    };
    Diag::new(DiagKind::Unsupported, reason).in_op(name)
}

/// Renders `(params) -> result` for a candidate signature.
fn render_signature(signature: &crate::kernel::signature::Signature) -> String {
    let params: Vec<String> = signature.params.iter().map(ToString::to_string).collect();
    let tail = if signature.variadic.is_some() {
        ", .."
    } else {
        ""
    };
    format!("({}{tail}) -> {}", params.join(", "), signature.ret)
}

/// Renders the argument kinds of a call as `int, list`.
fn render_kinds(args: &[crate::kernel::value::Value]) -> String {
    let kinds: Vec<String> = args.iter().map(|value| value.kind().to_string()).collect();
    kinds.join(", ")
}

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
    use crate::api::Session;
    use crate::kernel::call::CallCtx;
    use crate::kernel::call::OpError;
    use crate::kernel::value::Value;

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
        let ((), outer) = capture_all(|| {
            report(Diag::new(DiagKind::Shape, "outer one"));
            let ((), inner) = capture_all(|| report(Diag::new(DiagKind::Domain, "inner")));
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


    /// A heavy operator that declines used to answer `Unspecified`; the
    /// classifier now derives a kind from the descriptor, so a wrong argument
    /// count is named as such.
    #[test]
    fn a_wrong_argument_count_is_an_arity_finding() {
        let session = Session::new();
        let finding = match session.apply_diag("diff", &[], &CallCtx::new()) {
            | Err(finding) => finding,
            | Ok(value) => panic!("expected a finding, got {value:?}"),
        };
        assert_eq!(finding.kind, super::DiagKind::Arity);
        assert_eq!(finding.op, "diff");
        assert_ne!(finding.message, "");
    }

    /// The classifier never invents a numeric cause: with the right number of
    /// arguments but a form no rule accepts, the kind is one the caller can
    /// branch on, and never the catch-all.
    #[test]
    fn a_declined_call_is_classified_not_unspecified() {
        let session = Session::new();
        let args = [
            Value::Symbol("x".to_string()),
            Value::Symbol("x".to_string()),
            Value::Symbol("x".to_string()),
        ];
        if let Err(finding) = session.apply_diag("diff", &args, &CallCtx::new()) {
            assert_ne!(finding.kind, super::DiagKind::Unspecified, "{finding}");
            assert_eq!(finding.op, "diff");
        }
    }
}
