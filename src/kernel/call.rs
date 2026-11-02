//! The dynamic call ABI: context, errors and the `Kernel` trait.
//!
//! A dynamic call is a first-class operation, not a fallback path. It takes
//! an operator name and a list of [`Value`]s, and returns a [`Value`] or a
//! structured [`OpError`]. Everything the engine can do through the typed
//! API is reachable this way; the typed API is a static view of the same
//! thing.
//!
//! ```no_run
//! use rssn::kernel::CallCtx;
//! use rssn::kernel::Value;
//! use rssn::prelude::*;
//!
//! let s = Session::new();
//! let x = s.sym("x").to_value()?;
//! let out = s.call_value("sin", &[x], &CallCtx::new())?;
//! assert_eq!(out.to_string(), "sin(x)");
//! # Ok::<(), rssn::kernel::OpError>(())
//! ```

use std::fmt;

use crate::graph::Budget;
use crate::graph::DimensionMode;
use crate::graph::Facts;
use crate::kernel::codec::CodecError;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::OpKey;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::SessionId;
use crate::kernel::value::Value;

/// The phase a call asks for.
///
/// This is the boundary's spelling of
/// [`api::Target`](crate::api::Target); they are the same type.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
pub enum Target {
    /// A closed-form term: no unevaluated heavy operator left.
    #[default]
    Symbolic,
    /// A floating point value within the given absolute tolerance.
    Numeric {
        /// Requested absolute tolerance.
        tolerance: f64,
    },
}


/// One event of a call trace.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TraceEvent {
    /// What happened (`resolve`, `evaluate`, `compute`, `reduce`).
    pub kind: String,
    /// Details for a human.
    pub detail: String,
}

/// A record of what a dynamic call did.
///
/// A host that wants to know why a call was slow, or which overload it got,
/// passes a `Trace` in the [`CallCtx`] and reads it afterwards.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Trace {
    events: Vec<TraceEvent>,
}

impl Trace {
    /// An empty trace.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends an event.
    pub fn push(
        &mut self,
        kind: &str,
        detail: impl Into<String>,
    ) {
        self.events.push(TraceEvent {
            kind: kind.to_owned(),
            detail: detail.into(),
        });
    }

    /// The events, in order.
    #[must_use]
    pub fn events(&self) -> &[TraceEvent] {
        &self.events
    }

    /// Whether nothing was recorded.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
}

/// An assumption about a symbol that the call may rely on.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assumption {
    /// The symbol.
    pub symbol: String,
    /// What is known about it.
    pub facts: Facts,
}

/// Everything a dynamic call is allowed to know besides its arguments.
///
/// The context is a plain data structure: a plugin receives it, reads it and
/// cannot change the caller's state through it. `trace`, when present, is the
/// one exception — that is the point of a trace.
pub struct CallCtx<'a> {
    /// Which phase to aim for.
    pub target: Target,
    /// Search budget.
    pub budget: Budget,
    /// How strictly dimensions are checked.
    pub dimension_mode: DimensionMode,
    /// Facts the rules may rely on.
    pub assumptions: Vec<Assumption>,
    /// Numeric values for free symbols.
    pub bindings: Vec<(String, f64)>,
    /// Exact replacements for free symbols, applied before the call.
    pub substitutions: Vec<(String, Value)>,
    /// Host callbacks, addressed by [`CallbackRef`] inside the arguments.
    pub callbacks: Vec<std::sync::Arc<dyn crate::kernel::plugin::HostCallback>>,
    /// Where to record what happened. A `RefCell` because a trace is
    /// written through a shared reference: the context is read-only for the
    /// operator being called, and the trace is the one thing it appends to.
    pub trace: Option<&'a std::cell::RefCell<Trace>>,
}

impl fmt::Debug for CallCtx<'_> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("CallCtx")
            .field("target", &self.target)
            .field("budget", &self.budget)
            .field("dimension_mode", &self.dimension_mode)
            .field("assumptions", &self.assumptions.len())
            .field("bindings", &self.bindings.len())
            .field("substitutions", &self.substitutions.len())
            .field("callbacks", &self.callbacks.len())
            .field("trace", &self.trace.is_some())
            .finish()
    }
}

impl Default for CallCtx<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> CallCtx<'a> {
    /// A symbolic call with the default budget and advisory dimensions.
    #[must_use]
    pub fn new() -> Self {
        Self {
            target: Target::Symbolic,
            budget: Budget::default(),
            dimension_mode: DimensionMode::Advisory,
            assumptions: Vec::new(),
            bindings: Vec::new(),
            substitutions: Vec::new(),
            callbacks: Vec::new(),
            trace: None,
        }
    }

    /// Asks for a numeric answer within `tolerance`.
    #[must_use]
    pub const fn numeric(
        mut self,
        tolerance: f64,
    ) -> Self {
        self.target = Target::Numeric { tolerance };
        self
    }

    /// Sets the budget.
    #[must_use]
    pub const fn budget(
        mut self,
        budget: Budget,
    ) -> Self {
        self.budget = budget;
        self
    }

    /// Sets the dimensional mode.
    #[must_use]
    pub const fn dimensions(
        mut self,
        mode: DimensionMode,
    ) -> Self {
        self.dimension_mode = mode;
        self
    }

    /// Binds a symbol for the numeric phase.
    #[must_use]
    pub fn bind(
        mut self,
        symbol: &str,
        value: f64,
    ) -> Self {
        self.bindings.push((symbol.to_owned(), value));
        self
    }

    /// Replaces a free symbol by an exact value before the call.
    #[must_use]
    pub fn substitute(
        mut self,
        symbol: &str,
        value: Value,
    ) -> Self {
        self.substitutions.push((symbol.to_owned(), value));
        self
    }

    /// Assumes facts about a symbol.
    #[must_use]
    pub fn assume(
        mut self,
        symbol: &str,
        facts: Facts,
    ) -> Self {
        self.assumptions.push(Assumption {
            symbol: symbol.to_owned(),
            facts,
        });
        self
    }

    /// Registers a host callback and returns its handle.
    pub fn add_callback(
        &mut self,
        callback: std::sync::Arc<dyn crate::kernel::plugin::HostCallback>,
    ) -> CallbackRef {
        self.callbacks.push(callback);
        CallbackRef(u32::try_from(self.callbacks.len().saturating_sub(1)).unwrap_or(u32::MAX))
    }

    /// Records what happened, if the caller asked for a trace.
    pub fn trace(
        &self,
        kind: &str,
        detail: impl Into<String>,
    ) {
        if let Some(trace) = self.trace.as_ref() {
            trace.borrow_mut().push(kind, detail);
        }
    }

    /// Records events into `trace`.
    #[must_use]
    pub const fn with_trace(
        mut self,
        trace: &'a std::cell::RefCell<Trace>,
    ) -> Self {
        self.trace = Some(trace);
        self
    }
}

/// Why a dynamic call did not produce a value.
///
/// Every variant names the thing that went wrong, so a host can react to it
/// instead of matching on message text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OpError {
    /// No operator is registered under that key or name.
    UnknownOp(OpKey),
    /// The operator does not take that many arguments.
    Arity {
        /// What it accepts.
        expected: ArityInfo,
        /// What it got.
        got: usize,
    },
    /// An argument has the wrong type.
    ///
    /// The types are boxed: `TypeExpr` is the largest thing an error can
    /// carry, and an error is returned by value from every call in the
    /// dynamic ABI. Boxing them keeps `OpError` small enough that a
    /// `Result<Value, OpError>` costs no more than the `Value`.
    Type {
        /// What was expected.
        expected: Box<TypeExpr>,
        /// What was found.
        got: Box<TypeExpr>,
        /// Position of the argument.
        arg: usize,
    },
    /// More than one overload accepts the arguments and none is more
    /// specific.
    ///
    /// Refusing is deliberate: picking one silently is how a typed caller
    /// that chose statically ends up executing a different overload.
    Ambiguous {
        /// The operator's name.
        name: String,
        /// How many overloads tied.
        count: usize,
    },
    /// A dimensional conflict.
    Dimension(crate::graph::DimError),
    /// The value has no numeric form.
    NotNumeric(String),
    /// The call did not reduce to a closed form.
    NotReduced(String),
    /// The value is outside the operator's domain.
    Domain(String),
    /// The search budget ran out before an answer was found.
    BudgetExceeded,
    /// A backend cannot execute this operator.
    Unsupported {
        /// The backend's name.
        backend: String,
        /// The feature it lacks.
        feature: String,
    },
    /// A plugin failed.
    Plugin(crate::kernel::plugin::PluginError),
    /// A term handle belongs to another session.
    ForeignTerm {
        /// The session the handle names.
        expected: SessionId,
        /// The session that received it.
        got: SessionId,
    },
    /// A value could not be encoded or decoded.
    Codec(CodecError),
    /// A host callback failed.
    Host(String),
    /// An internal invariant was violated; a bug in rssn.
    Internal(String),
}

impl fmt::Display for OpError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::UnknownOp(key) => write!(f, "unknown operator `{key}`"),
            | Self::Arity { expected, got } => {
                write!(f, "operator takes {expected} arguments, got {got}")
            },
            | Self::Type { expected, got, arg } => {
                write!(f, "argument {arg} has type {got}, expected {expected}")
            },
            | Self::Ambiguous { name, count } => {
                write!(
                    f,
                    "`{name}` has {count} overloads that accept these arguments equally well"
                )
            },
            | Self::Dimension(e) => write!(f, "{e}"),
            | Self::NotNumeric(what) => write!(f, "`{what}` has no numeric value"),
            | Self::NotReduced(what) => write!(f, "`{what}` did not reduce to a closed form"),
            | Self::Domain(what) => write!(f, "outside the domain of the operator: {what}"),
            | Self::BudgetExceeded => f.write_str("the search budget was exhausted"),
            | Self::Unsupported { backend, feature } => {
                write!(f, "backend `{backend}` does not support {feature}")
            },
            | Self::Plugin(e) => write!(f, "{e}"),
            | Self::ForeignTerm { expected, got } => {
                write!(f, "term belongs to {expected}, not to {got}")
            },
            | Self::Codec(e) => write!(f, "encoding error: {e}"),
            | Self::Host(message) => write!(f, "host callback failed: {message}"),
            | Self::Internal(message) => write!(f, "internal error: {message}"),
        }
    }
}

impl std::error::Error for OpError {}

impl From<crate::graph::DimError> for OpError {
    fn from(value: crate::graph::DimError) -> Self {
        Self::Dimension(value)
    }
}

impl From<CodecError> for OpError {
    fn from(value: CodecError) -> Self {
        Self::Codec(value)
    }
}

impl From<crate::kernel::plugin::PluginError> for OpError {
    fn from(value: crate::kernel::plugin::PluginError) -> Self {
        Self::Plugin(value)
    }
}

/// The dynamic kernel: the three operations a host needs.
///
/// Implemented by [`Session`](crate::api::Session). The trait exists so that
/// a host can be written against an abstraction and so that a proxy (a
/// remote kernel, a recording wrapper) can stand in for a session.
pub trait Kernel {
    /// Calls `op` with `args`.
    ///
    /// # Errors
    /// See [`OpError`].
    fn call(
        &self,
        op: &OpKey,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError>;

    /// Calls a higher-order operator with host callbacks.
    ///
    /// # Errors
    /// See [`OpError`].
    fn call_ho(
        &self,
        op: &OpKey,
        args: &[Value],
        fns: &[CallbackRef],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError>;

    /// Registers a plugin.
    ///
    /// # Errors
    /// See [`PluginError`](crate::kernel::plugin::PluginError).
    fn register(
        &mut self,
        plugin: Box<dyn crate::kernel::plugin::Plugin>,
    ) -> Result<crate::kernel::plugin::PluginId, crate::kernel::plugin::PluginError>;
}
