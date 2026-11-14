//! # The user-facing API
//!
//! One entry point, [`Session::compute`], for everything that is an identity
//! transformation. What differs between "simplify", "differentiate",
//! "evaluate" and "solve" is only the term that is handed in and the
//! [`Target`](crate::kernel::call::Target) phase of the answer.
//!
//! ```
//! use rssn::api::*;
//!
//! let s = Session::new();
//! let x = s.sym("x");
//! let request = diff(sin(x) * cos(x), x);
//!
//! // Symbolic phase: a closed form.
//! let symbolic = s.compute(request, &Config::new())?;
//! assert!(symbolic.reduced);
//!
//! // Numeric phase: the same request, evaluated at x = 0.
//! let numeric = s.compute(request, &Config::new().numeric(1e-9).bind("x", 0.0))?;
//! assert_eq!(numeric.value, Some(1.0));
//! # Ok::<(), ComputeError>(())
//! ```
//!
//! # Limitations and coverage
//!
//! This layer is the *identity-transformation* part of the crate: one entry
//! point over a general engine. In priority order:
//!
//! * **Reduction only.** [`Session::compute`](crate::api::Session::compute)
//!   replaces a term by an equal term; it does not simulate, optimise or solve
//!   an arbitrary numerical problem. Time stepping and discretised fields live
//!   in [`crate::sim`], the standalone numeric routines in [`crate::kernels`].
//!   "Solve" and "integrate" are ordinary operators handed through the same
//!   compute path, so what they do is whatever the installed rule sets derive.
//! * **One session per thread.** A [`Session`](crate::api::Session) is
//!   single-threaded (its store is a `RefCell`). Use
//!   [`Session::fork`](crate::api::Session::fork) to move work to another
//!   thread and [`Session::merge`](crate::api::Session::merge) to bring the
//!   result back; terms never leave the session that minted them. A
//!   [`Term`](crate::api::Term) borrows its session and cannot outlive it; for
//!   language bindings there is
//!   [`TermRef`](crate::kernel::value::TermRef) (a session id plus a node id),
//!   which [`Session::term_ref`](crate::api::Session::term_ref) checks before
//!   use.
//! * **No cancellation; the wall clock is a budget field.** The limits are
//!   the [`Budget`](crate::graph::Budget) (`max_iterations`, `max_nodes`, and
//!   the optional `time_limit`, which ends a run with
//!   [`Stop::Time`](crate::graph::Stop)); a pathological rule set runs until
//!   one of them is exhausted. There is no streaming API: a computation is one
//!   call returning one answer.
//! * **Answers are cached; lemmas can be.** A session memoises whole answers
//!   by request node and configuration fingerprint (exact: a hit equals a
//!   fresh run), and can reuse the closed forms of heavy sub-requests across
//!   requests; see [`CachePolicy`](crate::api::CachePolicy) and
//!   [`Session::cache_stats`](crate::api::Session::cache_stats).
//! * **Configurations are per request; rule-set combinations are cached.**
//!   Every distinct set of rule sets used by a [`Config`](crate::api::Config)
//!   installs a template (a copy of the graph plus an engine) the first time it
//!   is seen, and the session keeps it for reuse. Cloning a `Config` per call
//!   is cheap, but a host that creates unboundedly many distinct configurations
//!   grows the session.
//! * **Answers are imported into the store.**
//!   [`Session::compute`](crate::api::Session::compute) imports the chosen
//!   answer as a term of the session, so the store grows with computations.
//!   Handles are `u32` node indices and are not stable across sessions or
//!   forks. [`Session::compute_all`](crate::api::Session::compute_all) imports
//!   in input order, which is why its answers equal a loop of `compute`.
//! * **Dimension assumptions are declarations, not validations.** Text
//!   assumptions ([`Config::dimension`](crate::api::Config::dimension)) are
//!   parsed when a request runs: an unknown unit or a malformed expression
//!   becomes a diagnostic, and only
//!   [`DimensionMode::Strict`](crate::graph::DimensionMode) turns a *blocking*
//!   finding into
//!   [`ComputeError::Dimension`](crate::api::ComputeError::Dimension). The
//!   findings are carried as
//!   [`Answer::diagnostics`](crate::api::Answer::diagnostics); no individual
//!   `DimError` code is part of the C ABI. Structured assumptions
//!   ([`Config::dimension_term`](crate::api::Config::dimension_term)) never
//!   parse text and are what a generated typed front end should use.
//! * **`Value` conversions are total in one direction only.** Every term has
//!   a [`Value`](crate::kernel::value::Value)
//!   ([`Term::to_value`](crate::api::Term::to_value),
//!   [`Session::value_of`](crate::api::Session::value_of)); a value without a
//!   handle has a term
//!   ([`Session::term_of`](crate::api::Session::term_of)) *except* the
//!   boundary kinds: a bare
//!   [`Value::Unit`](crate::kernel::value::Value::Unit), a complex literal
//!   with a nonzero imaginary part, and sparse, record, variant, function and
//!   opaque values.
//!   [`ValueKindInfo::all`](crate::kernel::value::ValueKindInfo::all) is the
//!   machine-readable table of which kind can enter the term graph. `Value`
//!   serialises completely
//!   ([`Value::to_cbor`](crate::kernel::value::Value::to_cbor)), but a
//!   serialised `Value::Term` handle is only meaningful to the session that
//!   minted it.
//! * **Units resolve against the session.**
//!   [`Session::value_of`](crate::api::Session::value_of) reads a quantity's
//!   unit through the session's registry (the built-in table plus
//!   [`Session::register_unit`](crate::api::Session::register_unit)), so a unit
//!   the session registered keeps its dimension and display. A compound unit
//!   expression (`m/s`) keeps its resolved dimension and scale, but not a full
//!   definition; and the process-wide table is visible to every session,
//!   session units only to the one that registered them.
//! * **Text in, text out.** The parser accepts this crate's infix grammar
//!   only; LaTeX and Typst are output formats
//!   ([`Term::to_latex`](crate::api::Term::to_latex),
//!   [`Term::to_typst`](crate::api::Term::to_typst)) and are not accepted as
//!   input. [`Term::to_pretty`](crate::api::Term::to_pretty) is for display,
//!   not round-tripping; `Term::to_string` is the parseable form.
//! * **Known edges.** [`Session::rational`](crate::api::Session::rational)
//!   returns `None` for a zero denominator; a numeric target without enough
//!   bindings is
//!   [`ComputeError::NotNumeric`](crate::api::ComputeError::NotNumeric) and
//!   carries the best symbolic form in its message;
//!   [`Session::merge`](crate::api::Session::merge) fails when the two sessions
//!   declared different dimensions for one symbol rather than silently
//!   weakening an assumption; and
//!   [`Term::as_tensor`](crate::api::Term::as_tensor) rejects a ragged nested
//!   list instead of padding it.

use std::cell::RefCell;
use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use std::sync::Arc;

use num_rational::BigRational;
use rayon::prelude::*;

use crate::dim::DimError;
use crate::dim::DimTerm;
use crate::dim::ResolvedUnit;
use crate::dim::Scale;
use crate::dim::UnitDef;
use crate::dim::UnitRegistry;
use crate::dim::UnitRegistrySnapshot;
use crate::graph::Budget;
use crate::graph::ClosedForm;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Evaluated;
use crate::graph::Extractor;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::ParseError;
use crate::graph::Report;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Saturate;
use crate::graph::SizeCost;
use crate::graph::op::core;
use crate::kernel::manifest::KernelManifest;
use crate::kernel::value::ScaleValue;
use crate::kernel::value::UnitValue;
use crate::kernel::value::Value;
use crate::rules;

/// The representation the answer should be in.
///
/// This is [`crate::kernel::call::Target`]: the typed API and the dynamic
/// ABI ask for phases with the same type.
pub use crate::kernel::call::Target;

mod cache;
mod value;

pub use cache::CachePolicy;
pub use cache::CacheStats;
pub use value::MergeError;

/// What to compute with and what to aim for.
#[derive(Clone, Debug)]
pub struct Config {
    sets: Vec<RuleSet>,
    target: Target,
    bindings: Vec<(String, f64)>,
    assumptions: Vec<(String, Facts)>,
    budget: Budget,
    dimension_mode: crate::graph::DimensionMode,
    dimension_assumptions: Vec<(String, String)>,
    /// Structured dimension assumptions, applied after the text ones.
    dimension_terms: Vec<(String, crate::kernel::signature::DimExpr)>,
    base_dimensions: Vec<String>,
    bypass_dimensions: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self::new()
    }
}

impl Config {
    /// All standard rule sets, symbolic target, default budget.
    #[must_use]
    pub fn new() -> Self {
        Self {
            sets: rules::standard(),
            target: Target::Symbolic,
            bindings: Vec::new(),
            assumptions: Vec::new(),
            budget: Budget::default(),
            dimension_mode: crate::graph::DimensionMode::Advisory,
            dimension_assumptions: Vec::new(),
            dimension_terms: Vec::new(),
            base_dimensions: Vec::new(),
            bypass_dimensions: false,
        }
    }

    /// No rule sets at all; add some with [`Config::with`].
    #[must_use]
    pub fn empty() -> Self {
        Self {
            sets: Vec::new(),
            target: Target::Symbolic,
            bindings: Vec::new(),
            assumptions: Vec::new(),
            budget: Budget::default(),
            dimension_mode: crate::graph::DimensionMode::Advisory,
            dimension_assumptions: Vec::new(),
            dimension_terms: Vec::new(),
            base_dimensions: Vec::new(),
            bypass_dimensions: false,
        }
    }

    /// Adds a rule set (and, implicitly, its dependencies).
    #[must_use]
    pub fn with(
        mut self,
        set: RuleSet,
    ) -> Self {
        self.sets.push(set);
        self
    }

    /// Asks for a closed-form answer.
    #[must_use]
    pub const fn symbolic(mut self) -> Self {
        self.target = Target::Symbolic;
        self
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

    /// Binds a symbol to a value for numeric evaluation.
    #[must_use]
    pub fn bind(
        mut self,
        symbol: &str,
        value: f64,
    ) -> Self {
        self.bindings.push((symbol.to_owned(), value));
        self
    }

    /// Declares facts about a symbol — for example
    /// [`Facts::POSITIVE`] or [`Facts::REAL`] — that conditional identities
    /// may rely on.
    #[must_use]
    pub fn assume(
        mut self,
        symbol: &str,
        facts: Facts,
    ) -> Self {
        self.assumptions.push((symbol.to_owned(), facts));
        self
    }

    /// The search budget in effect.
    #[must_use]
    pub fn current_budget(&self) -> Budget {
        self.budget.clone()
    }

    /// Sets the dimensional checking mode. The default is
    /// [`DimensionMode::Advisory`](crate::graph::DimensionMode).
    #[must_use]
    pub const fn dimensions(
        mut self,
        mode: crate::graph::DimensionMode,
    ) -> Self {
        self.dimension_mode = mode;
        self
    }

    /// The dimensional checking mode.
    #[must_use]
    pub const fn dimension_mode(&self) -> crate::graph::DimensionMode {
        self.dimension_mode
    }

    /// Declares the dimension of a symbol as text: a unit expression
    /// (`"m/s"`, `"km/h"`), a base-dimension product (`"Length/Time"`), an
    /// explicit base dimension (`"base:T"`, needed when a custom base has
    /// the same name as a unit), or `"any"`.
    ///
    /// Unit symbols take precedence over base names for bare text, so `"T"`
    /// is always the tesla. Use `"base:Name"` to force a base dimension.
    #[must_use]
    pub fn dimension(
        mut self,
        symbol: &str,
        text: &str,
    ) -> Self {
        self.dimension_assumptions
            .push((symbol.to_owned(), text.to_owned()));
        self
    }

    /// Declares a base dimension the engine does not have.
    ///
    /// The seven SI bases always exist. A unit on a base of the host's own
    /// (`unit EUR : Currency`) needs that base in the graph's arena before its
    /// dimension can be interned, and this is the declaration that puts it
    /// there. It is idempotent: declaring the same name twice is harmless.
    ///
    /// [`Config::dimension_term`] interns the bases it names on its own, so
    /// this is only needed when a *unit* — rather than a symbol — is what
    /// carries the custom dimension.
    #[must_use]
    pub fn base_dimension(
        mut self,
        name: &str,
    ) -> Self {
        self.base_dimensions.push(name.to_owned());
        self
    }

    /// Declares the dimension of a symbol from a structured expression.
    ///
    /// This is the entry a generated typed front end should use. Unlike
    /// [`Config::dimension`] it never parses text, so a unit symbol and a base
    /// dimension that happen to share a name cannot be confused, and a base
    /// dimension the engine has never seen (`Currency`) is defined on the way
    /// in rather than rejected.
    ///
    /// Structured assumptions are applied after the text ones, so the two can
    /// be mixed and the structured one wins.
    #[must_use]
    pub fn dimension_term(
        mut self,
        symbol: &str,
        dim: crate::kernel::signature::DimExpr,
    ) -> Self {
        self.dimension_terms.push((symbol.to_owned(), dim));
        self
    }

    /// Declares the dimension of the value of a function symbol. This is the
    /// same declaration as [`Config::dimension`]; it exists for readability
    /// at call sites that describe function values.
    #[must_use]
    pub fn dimension_of(
        self,
        symbol: &str,
        text: &str,
    ) -> Self {
        self.dimension(symbol, text)
    }

    /// Registers a user-defined base dimension (for example `Currency`).
    #[must_use]
    pub fn define_base_dimension(
        mut self,
        name: &str,
    ) -> Self {
        self.base_dimensions.push(name.to_owned());
        self
    }

    /// Disables hard dimensional errors even in
    /// [`DimensionMode::Strict`](crate::graph::DimensionMode): dimensions are
    /// still inferred and reported, but a conflict never fails the computation.
    /// This is the explicit escape hatch for intentionally non-physical
    /// formulas.
    #[must_use]
    pub const fn bypass_dimension_checks(mut self) -> Self {
        self.bypass_dimensions = true;
        self
    }

    /// Replaces the search budget.
    #[must_use]
    pub const fn budget(
        mut self,
        budget: Budget,
    ) -> Self {
        self.budget = budget;
        self
    }
}

/// Why a computation could not produce an answer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ComputeError {
    /// A rule set failed to install.
    Rule(RuleError),
    /// Input text could not be parsed.
    Parse(ParseError),
    /// An operator name is not registered.
    UnknownOperator(String),
    /// A numeric answer was requested but the term could not be evaluated;
    /// carries the best symbolic form that was reached.
    NotNumeric(String),
    /// A dimension conflict in strict mode.
    Dimension(crate::graph::DimError),
    /// The rules derived that two different numbers are equal. This is a
    /// bug in a rule set; the session must be discarded.
    Inconsistent(String),
}

impl fmt::Display for ComputeError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Rule(e) => write!(f, "{e}"),
            | Self::Parse(e) => write!(f, "{e}"),
            | Self::UnknownOperator(name) => write!(f, "unknown operator `{name}`"),
            | Self::NotNumeric(term) => {
                write!(f, "no numeric value could be determined for `{term}`")
            },
            | Self::Dimension(error) => write!(f, "{error}"),
            | Self::Inconsistent(what) => write!(f, "rule set is unsound: derived {what}"),
        }
    }
}

impl std::error::Error for ComputeError {}

impl From<RuleError> for ComputeError {
    fn from(value: RuleError) -> Self {
        Self::Rule(value)
    }
}

impl From<ParseError> for ComputeError {
    fn from(value: ParseError) -> Self {
        Self::Parse(value)
    }
}

/// The units one session knows: the built-in table plus the units that
/// session registered.
///
/// The registry belongs to the session, not to the process. Every engine of
/// the session is given a copy (see [`Engine::set_unit_registry`]) before its
/// manifest is built, so a unit registered here is described by this
/// session's [`KernelManifest`] and by no other session's. Registering a unit
/// therefore does not go through [`crate::dim::units::register`], the table
/// every process shares.
///
/// A unit a *plugin* registers arrives here too: [`Template::with_units`]
/// hands this registry to the plugin's [`Registry`](crate::kernel::plugin::Registry),
/// so the plugin's units are the session's and not the process's.
#[derive(Clone, Debug, Default)]
pub(super) struct SessionUnits {
    /// Units registered into this session, in registration order.
    extra: Vec<UnitDef>,
}

impl SessionUnits {
    /// The registry an engine's manifest describes.
    ///
    /// The engine keeps no handle to its registry, only the manifest built
    /// from it, so this is how a session carries its units forward when its
    /// configuration changes — a unit is added, or a template is built for
    /// another set of rule sets.
    pub(super) fn from_manifest(manifest: &KernelManifest) -> Self {
        let global = crate::dim::units::all();
        let extra = manifest
            .units
            .iter()
            .map(unit_def)
            // What the shared table already describes is not this session's,
            // unless this session redefined it.
            .filter(|unit| !global.contains(unit))
            .collect();
        Self { extra }
    }
}

/// A unit as the manifest carries it, back as a definition.
fn unit_def(unit: &UnitValue) -> UnitDef {
    fn scale(value: &ScaleValue) -> Scale {
        Scale {
            rat: BigRational::new(value.numer.clone(), value.denom.clone()),
            pi: value.pi,
        }
    }
    UnitDef {
        name: unit.name.clone(),
        aliases: unit.aliases.clone(),
        dimension: unit.dimension.clone(),
        scale: scale(&unit.scale),
        offset: unit
            .offset
            .as_ref()
            .map(|offset| BigRational::new(offset.numer.clone(), offset.denom.clone())),
        prefixable: unit.prefixable,
        affine: unit.affine,
        display: unit.display.clone(),
        category: unit.category.clone(),
        since: unit.since.clone(),
    }
}

impl UnitRegistry for SessionUnits {
    fn list(&self) -> Vec<UnitDef> {
        let mut units = crate::dim::units::all();
        for unit in &self.extra {
            match units.iter_mut().find(|known| known.name == unit.name) {
                | Some(slot) => *slot = unit.clone(),
                | None => units.push(unit.clone()),
            }
        }
        units.sort_by(|a, b| a.name.cmp(&b.name));
        units
    }

    fn lookup(
        &self,
        name: &str,
    ) -> Option<UnitDef> {
        self.extra
            .iter()
            .find(|unit| unit.name == name || unit.aliases.iter().any(|alias| alias == name))
            .cloned()
            .or_else(|| crate::dim::units::lookup(name))
    }

    fn resolve(
        &self,
        expr: &str,
    ) -> Result<ResolvedUnit, DimError> {
        // A single unit, this session's own included; a compound expression
        // is resolved by the shared parser, which does not know this
        // session's units yet.
        if let Some(unit) = self.lookup(expr.trim()) {
            return Ok(ResolvedUnit {
                dimension: unit.dimension,
                scale: unit.scale,
                offset: unit.offset,
                affine: unit.affine,
            });
        }
        crate::dim::units::resolve(expr)
    }

    fn register(
        &mut self,
        unit: UnitDef,
    ) -> Result<(), DimError> {
        if unit.name.trim().is_empty() {
            return Err(DimError::Parse {
                text: unit.name,
                reason: String::from("a unit needs a name"),
            });
        }
        let clash = |other: &UnitDef| other.name != unit.name && other.dimension != unit.dimension;
        if let Some(known) = self.extra.iter().find(|known| known.name == unit.name)
            && clash(known)
        {
            return Err(DimError::Parse {
                text: unit.name.clone(),
                reason: String::from("a different unit already has this name"),
            });
        }
        if let Some(alias) = unit.aliases.iter().find(|alias| {
            self.extra
                .iter()
                .any(|known| known.aliases.contains(alias) && clash(known))
        }) {
            return Err(DimError::Parse {
                text: alias.clone(),
                reason: String::from("a different unit already has this alias"),
            });
        }
        self.extra.retain(|known| known.name != unit.name);
        self.extra.push(unit);
        Ok(())
    }

    fn format(
        &self,
        dim: &DimTerm,
        value: &Value,
    ) -> String {
        crate::dim::StandardRegistry.format(dim, value)
    }

    fn snapshot(&self) -> UnitRegistrySnapshot {
        let units = self.list();
        UnitRegistrySnapshot {
            hash: crate::kernel::codec::digest(&units),
            units,
        }
    }
}

/// A graph with a set of rule sets installed, and the engine bound to it.
struct Template {
    sets: Vec<Arc<str>>,
    /// The rule sets themselves, so that a dynamic call can build a
    /// [`Config`] that selects this template.
    rule_sets: Vec<RuleSet>,
    graph: Graph,
    engine: Engine,
}

struct Inner {
    /// The term store: every [`Term`] is a node of this graph. It is a
    /// plain hash-consed DAG; no equalities are ever asserted in it.
    graph: Graph,
    templates: Vec<Template>,
    /// Plugins registered into this session. They are part of the
    /// configuration, so a template that is built later installs them too.
    plugins: Vec<Arc<dyn crate::kernel::plugin::Plugin>>,
    /// Ids handed out to plugins, in registration order.
    plugin_count: u32,
    /// Answers and lemmas of earlier computations (see [`CachePolicy`]).
    cache: cache::ComputeCache,
}

/// Owns the term store and the engines that compute on it.
///
/// Each [`Session::compute`] call runs in a private copy of a template
/// graph: the request is copied in, the search happens there, and only the
/// answer is copied back. Requests therefore cannot slow each other down or
/// influence each other's results, and a session can live as long as the
/// program.
///
/// A session is single-threaded; create one per thread. Terms built in one
/// session must not be used in another.
pub struct Session {
    inner: RefCell<Inner>,
    id: crate::kernel::value::SessionId,
}

impl fmt::Debug for Session {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let inner = self.inner.borrow();
        write!(
            f,
            "Session({} terms, {} rule configurations)",
            inner.graph.len(),
            inner.templates.len()
        )
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

/// A handle to a term of a [`Session`]. Cheap to copy.
#[derive(Copy, Clone)]
pub struct Term<'s> {
    session: &'s Session,
    node: NodeId,
}

/// The result of [`Session::compute`].
#[derive(Clone, Debug)]
pub struct Answer<'s> {
    /// The answer as a term: a closed form for a symbolic target, a float
    /// literal for a numeric one.
    pub term: Term<'s>,
    /// The numeric value, when one is known.
    pub value: Option<f64>,
    /// Estimated absolute error of `value`.
    pub error: Option<f64>,
    /// Whether the request was fully reduced. `false` means `term` still
    /// contains unevaluated heavy operators.
    pub reduced: bool,
    /// The inferred dimension of `term` (`None` when dimensional checking is
    /// off or no dimension could be inferred).
    pub dimension: Option<crate::graph::Dimension>,
    /// Dimensional findings of an advisory run (conflicts, coverage gaps,
    /// unresolved symbols).
    pub diagnostics: Vec<crate::graph::DimError>,
    /// What the engine did.
    pub report: Report,
}

impl Session {
    /// A session with the operators of every standard rule set registered.
    ///
    /// # Panics
    /// Panics if the built-in rule sets fail to install, which would be a
    /// bug in rssn itself.
    #[must_use]
    pub fn new() -> Self {
        match Self::with_rules(&rules::standard()) {
            | Ok(session) => session,
            | Err(e) => broken_builtin(&e),
        }
    }

    /// A session with exactly the given rule sets pre-installed.
    ///
    /// # Errors
    /// Fails when a rule set does not install.
    pub fn with_rules(sets: &[RuleSet]) -> Result<Self, RuleError> {
        let template = Template::new(sets, &[])?;
        let graph = template.graph.clone();
        Ok(Self {
            id: next_session_id(),
            inner: RefCell::new(Inner {
                graph,
                templates: vec![template],
                plugins: Vec::new(),
                plugin_count: 0,
                cache: cache::ComputeCache::default(),
            }),
        })
    }

    const fn term(
        &self,
        node: NodeId,
    ) -> Term<'_> {
        Term { session: self, node }
    }

    /// Registers a unit in this session.
    ///
    /// The unit is described by this session's [manifest](Session::manifest)
    /// and by no other session's: it goes into the session's own
    /// [`UnitRegistry`], not into the table every session shares. The built-in
    /// units stay available.
    ///
    /// # Errors
    /// Fails when the name or an alias already belongs to a unit with a
    /// different dimension.
    pub fn register_unit(
        &mut self,
        unit: UnitDef,
    ) -> Result<(), DimError> {
        let mut inner = self.inner.borrow_mut();
        // A unit changes how every template resolves quantities.
        inner.cache.invalidate();
        // The clone in the loop is only for the templates that come after;
        // the final template takes ownership so the unit is consumed once.
        let Some((last, earlier)) = inner.templates.split_last_mut() else {
            return Ok(());
        };
        for template in earlier {
            let mut units = SessionUnits::from_manifest(template.engine.manifest());
            units.register(unit.clone())?;
            template.engine.set_unit_registry(units);
            template.engine.refresh_manifest(&template.graph);
        }
        let mut units = SessionUnits::from_manifest(last.engine.manifest());
        units.register(unit)?;
        last.engine.set_unit_registry(units);
        last.engine.refresh_manifest(&last.graph);
        Ok(())
    }

    /// The symbol named `name`.
    pub fn sym(
        &self,
        name: &str,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.sym(name))
    }

    /// An exact integer.
    pub fn int(
        &self,
        value: i64,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.int(value))
    }

    /// The exact fraction `numerator / denominator`; `None` for a zero
    /// denominator.
    pub fn rational(
        &self,
        numerator: i64,
        denominator: i64,
    ) -> Option<Term<'_>> {
        let n = Number::fraction(numerator, denominator)?;
        Some(self.term(self.inner.borrow_mut().graph.num(n)))
    }

    /// A floating point number.
    pub fn float(
        &self,
        value: f64,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.float(value))
    }

    /// The term with handle `id` (see [`Term::id`]); `None` when this
    /// session has no such term.
    #[must_use]
    pub fn term_by_id(
        &self,
        id: u32,
    ) -> Option<Term<'_>> {
        let node = NodeId::from_raw(id);
        (node.index() < self.inner.borrow().graph.len()).then(|| self.term(node))
    }

    /// Parses a term from infix text.
    ///
    /// # Errors
    /// Returns the parser's complaint.
    pub fn parse(
        &self,
        text: &str,
    ) -> Result<Term<'_>, ComputeError> {
        Ok(self.term(self.inner.borrow_mut().graph.parse(text)?))
    }

    /// Applies the operator named `op` to `args`.
    ///
    /// # Errors
    /// Fails when the operator is unknown or does not take that many
    /// arguments.
    pub fn call<'s>(
        &'s self,
        op: &str,
        args: &[Term<'s>],
    ) -> Result<Term<'s>, ComputeError> {
        let mut inner = self.inner.borrow_mut();
        let unknown = || ComputeError::UnknownOperator(op.to_owned());
        let id = inner.graph.ops().lookup(op).ok_or_else(unknown)?;
        let nodes: Vec<NodeId> = args.iter().map(|t| t.node).collect();
        inner
            .graph
            .try_node(id, &nodes)
            .map(|n| self.term(n))
            .ok_or_else(unknown)
    }

    /// Reduces `term` to the target phase of `config`.
    ///
    /// # Errors
    /// * [`ComputeError::Rule`] if a rule set of `config` fails to install;
    /// * [`ComputeError::NotNumeric`] if a numeric target cannot be met;
    /// * [`ComputeError::Inconsistent`] if an unsound rule was detected.
    pub fn compute<'s>(
        &'s self,
        term: Term<'s>,
        config: &Config,
    ) -> Result<Answer<'s>, ComputeError> {
        let mut guard = self.inner.borrow_mut();
        let inner = &mut *guard;
        let index = template_index(inner, config)?;
        let key = cache::fingerprint(config, true);
        if let Some(hit) = inner.cache.answer(index, term.node, &key) {
            return hit.map(|stored| self.answer(stored));
        }
        let job = prepare(inner, index, term.node, config);
        let Some(template) = inner.templates.get(index) else {
            return Err(ComputeError::UnknownOperator("engine".to_owned()));
        };
        let done = run_job(&template.engine, job, config);
        let stored = finish(inner, index, config, done);
        inner
            .cache
            .store_answer(index, term.node, key, stored.clone());
        stored.map(|stored| self.answer(stored))
    }

    /// Reduces every term of `terms` under the same `config`, running the
    /// searches in parallel.
    ///
    /// The result equals `terms.iter().map(|&t| self.compute(t, config))`:
    /// each request is copied into its own private clone of the template
    /// (sequentially, because the store is borrowed), the engines run
    /// concurrently on the rayon pool, and the answers are extracted and
    /// copied back into the store in input order. The answers do not depend
    /// on the thread count. A failure in one request does not affect the
    /// others. Requests already in the answer cache are not run again; a
    /// request repeated inside `terms` runs once per occurrence, as a loop of
    /// [`Session::compute`] calls before the first answer was stored would.
    ///
    /// # Errors
    /// One entry per term; see [`Session::compute`].
    pub fn compute_all<'s>(
        &'s self,
        terms: &[Term<'s>],
        config: &Config,
    ) -> Vec<Result<Answer<'s>, ComputeError>> {
        let mut guard = self.inner.borrow_mut();
        let inner = &mut *guard;
        let index = match template_index(inner, config) {
            | Ok(i) => i,
            | Err(e) => return terms.iter().map(|_| Err(e.clone())).collect(),
        };
        let key = cache::fingerprint(config, true);
        let mut slots: Vec<Option<Result<cache::StoredAnswer, ComputeError>>> = terms
            .iter()
            .map(|t| inner.cache.answer(index, t.node, &key))
            .collect();
        let pending: Vec<usize> = (0..terms.len())
            .filter(|&i| slots.get(i).is_some_and(Option::is_none))
            .collect();
        let jobs: Vec<Job> = pending
            .iter()
            .filter_map(|&i| terms.get(i))
            .map(|t| prepare(inner, index, t.node, config))
            .collect();
        let Some(template) = inner.templates.get(index) else {
            return terms
                .iter()
                .map(|_| Err(ComputeError::UnknownOperator("engine".to_owned())))
                .collect();
        };
        let engine = &template.engine;
        let done: Vec<Done> = if jobs.len() >= COMPUTE_ALL_PAR_MIN {
            jobs.into_par_iter()
                .map(|job| run_job(engine, job, config))
                .collect()
        } else {
            jobs.into_iter()
                .map(|job| run_job(engine, job, config))
                .collect()
        };
        for (&i, d) in pending.iter().zip(done) {
            let stored = finish(inner, index, config, d);
            if let Some(t) = terms.get(i) {
                inner
                    .cache
                    .store_answer(index, t.node, key.clone(), stored.clone());
            }
            if let Some(slot) = slots.get_mut(i) {
                *slot = Some(stored);
            }
        }
        slots
            .into_iter()
            .map(|slot| {
                slot.unwrap_or_else(|| Err(ComputeError::UnknownOperator("engine".to_owned())))
                    .map(|stored| self.answer(stored))
            })
            .collect()
    }

    /// The cache policy in effect.
    #[must_use]
    pub fn cache_policy(&self) -> CachePolicy {
        self.inner.borrow().cache.policy
    }

    /// Replaces the cache policy. Entries already stored stay valid (they
    /// are exact) unless the new policy disables their table, which drops
    /// them.
    pub fn set_cache_policy(
        &self,
        policy: CachePolicy,
    ) {
        let mut inner = self.inner.borrow_mut();
        let keep = inner.cache.policy.answers == policy.answers
            && inner.cache.policy.lemmas == policy.lemmas
            && policy.capacity >= inner.cache.policy.capacity;
        if !keep {
            inner.cache.invalidate();
        }
        inner.cache.policy = policy;
    }

    /// Counters of the compute cache.
    #[must_use]
    pub fn cache_stats(&self) -> CacheStats {
        self.inner.borrow().cache.stats()
    }

    /// Drops every cached answer and lemma.
    pub fn clear_cache(&self) {
        self.inner.borrow_mut().cache.invalidate();
    }

    /// An [`Answer`] borrowing this session for a stored answer.
    fn answer(
        &self,
        stored: cache::StoredAnswer,
    ) -> Answer<'_> {
        Answer {
            term: self.term(stored.node),
            value: stored.value,
            error: stored.error,
            reduced: stored.reduced,
            dimension: stored.dimension,
            diagnostics: stored.diagnostics,
            report: stored.report,
        }
    }

    /// Shorthand for a symbolic [`Session::compute`] with the standard
    /// rules, returning only the term.
    ///
    /// # Errors
    /// See [`Session::compute`].
    pub fn simplify<'s>(
        &'s self,
        term: Term<'s>,
    ) -> Result<Term<'s>, ComputeError> {
        Ok(self.compute(term, &Config::new())?.term)
    }
}

/// Number of requests from which [`Session::compute_all`] uses the pool.
const COMPUTE_ALL_PAR_MIN: usize = 2;

/// A request copied into its private graph, ready to run.
struct Job {
    graph: Graph,
    root: NodeId,
    env: Env,
    dim_errors: Vec<crate::graph::DimError>,
    /// Record closed forms of heavy sub-requests after the run.
    capture_lemmas: bool,
}

/// A request after the search: the answer is chosen but not yet copied back
/// into the store.
struct Done {
    graph: Graph,
    env: Env,
    report: Report,
    target: Target,
    failure: Option<ComputeError>,
    best: NodeId,
    reduced: bool,
    witness: Option<crate::graph::Ball>,
    dimension: Option<crate::graph::Dimension>,
    diagnostics: Vec<crate::graph::DimError>,
    /// `(heavy, closed)` pairs of the private graph, for the lemma cache.
    lemmas: Vec<(NodeId, NodeId)>,
}

/// Most heavy sub-requests one run records as lemmas.
const LEMMAS_PER_RUN: usize = 64;
/// Most nodes walked looking for heavy sub-requests.
const LEMMA_WALK: usize = 4096;

/// The heavy operator nodes of the concrete term `root`, in discovery order.
fn heavy_subterms(
    graph: &Graph,
    root: NodeId,
) -> Vec<NodeId> {
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![root];
    let mut heavy = Vec::new();
    while let Some(node) = stack.pop() {
        if seen.len() >= LEMMA_WALK || !seen.insert(node) {
            continue;
        }
        if graph
            .ops()
            .get(graph.op(node))
            .flags
            .has(crate::graph::op::OpFlags::HEAVY)
        {
            heavy.push(node);
        }
        stack.extend(graph.children(node).iter().rev().copied());
    }
    heavy
}

const _: fn() = || {
    const fn send<T: Send>() {}
    const fn sync<T: Sync>() {}
    send::<Job>();
    send::<Done>();
    sync::<Engine>();
    sync::<Config>();
};

/// Index of the template for `config`'s rule sets, installing it on first use.
fn template_index(
    inner: &mut Inner,
    config: &Config,
) -> Result<usize, ComputeError> {
    let key: Vec<Arc<str>> = config.sets.iter().map(|s| Arc::from(s.name())).collect();
    match inner.templates.iter().position(|t| t.sets == key) {
        | Some(i) => Ok(i),
        | None => {
            // A template built for another set of rule sets describes the
            // units the session already knows, not only the built-in ones.
            let units = inner
                .templates
                .first()
                .map_or_else(SessionUnits::default, |t| {
                    SessionUnits::from_manifest(t.engine.manifest())
                });
            inner
                .templates
                .push(Template::with_units(&config.sets, &[], units)?);
            Ok(inner.templates.len().saturating_sub(1))
        },
    }
}

/// Copies the request into a private clone of the template.
fn prepare(
    inner: &mut Inner,
    index: usize,
    node: NodeId,
    config: &Config,
) -> Job {
    let mut graph = inner
        .templates
        .get(index)
        .map_or_else(Graph::new, |t| t.graph.clone());
    for name in &config.base_dimensions {
        graph.define_base_dimension(name);
    }
    let root = graph.import(&inner.graph, node);
    let mut env = match config.target {
        | Target::Symbolic => Env::symbolic(),
        | Target::Numeric { tolerance } => Env::numeric(tolerance),
    };
    for (name, value) in &config.bindings {
        env.bind(graph.interner_mut().symbol(name), *value);
    }
    for (name, facts) in &config.assumptions {
        let symbol = graph.interner_mut().symbol(name);
        graph.assume(symbol, *facts);
    }
    let mut dim_errors = Vec::new();
    for (name, text) in &config.dimension_assumptions {
        let symbol = graph.interner_mut().symbol(name);
        if let Err(error) = graph.assume_dimension(symbol, text) {
            dim_errors.push(error);
        }
    }
    // Structured assumptions come last so that they win over text.
    for (name, dim) in &config.dimension_terms {
        let symbol = graph.interner_mut().symbol(name);
        let term = graph.with_dim_arena_mut(|arena| dim.to_term(arena));
        match term {
            | Ok(term) => graph.assume_dimension_value(symbol, term),
            | Err(error) => dim_errors.push(error),
        }
    }
    // Lemmas: identities proven by earlier runs in the same context.
    let capture_lemmas = inner.cache.policy.lemmas;
    if capture_lemmas {
        let context = cache::fingerprint(config, false);
        let mut seeded = 0_usize;
        for heavy in heavy_subterms(&inner.graph, node) {
            if let Some(closed) = inner.cache.lemma(index, &context, heavy) {
                let a = graph.import(&inner.graph, heavy);
                let b = graph.import(&inner.graph, closed);
                if graph.union(a, b) {
                    seeded = seeded.saturating_add(1);
                }
            }
        }
        if seeded > 0 {
            graph.rebuild();
            inner.cache.note_lemma_hits(seeded);
        }
    }
    Job {
        graph,
        root,
        env,
        dim_errors,
        capture_lemmas,
    }
}

/// Runs the engine on the private graph and chooses the answer in it.
fn run_job(
    engine: &Engine,
    job: Job,
    config: &Config,
) -> Done {
    let Job {
        mut graph,
        root,
        env,
        dim_errors,
        capture_lemmas,
    } = job;
    let mut lemmas = Vec::new();
    let roots = [root];
    let report = match config.target {
        | Target::Symbolic => engine.run(&mut graph, &roots, &env, &Saturate, &config.budget),
        | Target::Numeric { .. } => {
            engine.run(&mut graph, &roots, &env, &Evaluated, &config.budget)
        },
    };
    let mut failure = graph
        .conflicts()
        .first()
        .map(|(a, b)| ComputeError::Inconsistent(format!("{a} = {b}")));
    let witness = graph.approx(graph.find(root));
    let (best, reduced) = if failure.is_some() {
        (root, false)
    } else {
        let extractor = Extractor::new(&graph, &roots, &ClosedForm);
        if capture_lemmas {
            for heavy in heavy_subterms(&graph, root) {
                if lemmas.len() >= LEMMAS_PER_RUN {
                    break;
                }
                if extractor.cost(&graph, heavy).is_some()
                    && let Some(closed) = extractor.build(&mut graph, heavy)
                {
                    lemmas.push((heavy, closed));
                }
            }
        }
        let closed = extractor.build(&mut graph, root);
        let reduced = closed.is_some();
        let best = match closed {
            | Some(node) => node,
            | None => {
                Extractor::new(&graph, &roots, &SizeCost)
                    .build(&mut graph, root)
                    .unwrap_or(root)
            },
        };
        (best, reduced)
    };
    let mut diagnostics = dim_errors;
    let dimension = if config.dimension_mode == crate::graph::DimensionMode::Off {
        None
    } else {
        let answer = graph.dimension(best);
        let input = if best == root {
            None
        } else {
            Some(graph.dimension(root))
        };
        // A structured request (a decomposition, a predicate, a container)
        // is typed `NonMetric` by its own rule. Its *expanded* literal then
        // holds parts with different dimensions on purpose: `lu` returns a
        // dimensionless permutation and `L` next to a dimensionful `U`.
        // Differences between those parts are not a mixed scalar list, so
        // the structural note is suppressed for such a request.
        let structured = input
            .as_ref()
            .is_some_and(|outcome| outcome.dimension.non_metric);
        let mut findings = answer.errors;
        if let Some(input) = &input {
            findings.extend(input.errors.iter().cloned());
        }
        if structured {
            findings.retain(|error| !matches!(error, crate::graph::DimError::MixedList { .. }));
        }
        diagnostics.extend(findings);
        let mut seen = std::collections::HashSet::new();
        diagnostics.retain(|error| seen.insert(error.to_string()));
        if failure.is_none()
            && config.dimension_mode == crate::graph::DimensionMode::Strict
            && !config.bypass_dimensions
            && let Some(error) = diagnostics.iter().find(|error| error.is_blocking())
        {
            failure = Some(ComputeError::Dimension(error.clone()));
        }
        Some(answer.dimension)
    };
    Done {
        graph,
        env,
        report,
        target: config.target,
        failure,
        best,
        reduced,
        witness,
        dimension,
        diagnostics,
        lemmas,
    }
}

/// Copies the chosen answer (and the lemmas of the run) back into the store.
fn finish(
    inner: &mut Inner,
    index: usize,
    config: &Config,
    done: Done,
) -> Result<cache::StoredAnswer, ComputeError> {
    let store = &mut inner.graph;
    if !done.lemmas.is_empty() {
        let context = cache::fingerprint(config, false);
        for &(heavy, closed) in &done.lemmas {
            let heavy = store.import(&done.graph, heavy);
            let closed = store.import(&done.graph, closed);
            inner.cache.store_lemma(index, &context, heavy, closed);
        }
    }
    let Done {
        graph,
        env,
        report,
        target,
        failure,
        best,
        reduced,
        witness,
        dimension,
        diagnostics,
        lemmas: _,
    } = done;
    if let Some(error) = failure {
        return Err(error);
    }
    match target {
        | Target::Symbolic => {
            let exact = graph.as_number(best).map(Number::to_f64);
            let node = store.import(&graph, best);
            Ok(cache::StoredAnswer {
                node,
                value: exact,
                error: exact.map(|_| 0.0),
                reduced,
                dimension,
                diagnostics,
                report,
            })
        },
        | Target::Numeric { .. } => {
            // Prefer evaluating the closed form; fall back to the
            // witness numeric kernels produced.
            let ball = graph
                .eval(best, &env)
                .map(crate::graph::Ball::exact)
                .or(witness);
            let Some(ball) = ball else {
                return Err(ComputeError::NotNumeric(graph.display(best)));
            };
            let node = store.float(ball.mid);
            Ok(cache::StoredAnswer {
                node,
                value: Some(ball.mid),
                error: Some(ball.rad),
                reduced,
                dimension,
                diagnostics,
                report,
            })
        },
    }
}

impl Template {
    fn new(
        sets: &[RuleSet],
        plugins: &[Arc<dyn crate::kernel::plugin::Plugin>],
    ) -> Result<Self, RuleError> {
        Self::with_units(sets, plugins, SessionUnits::default())
    }

    /// Builds a template whose engine describes `units` rather than the
    /// process-wide table.
    pub(super) fn with_units(
        sets: &[RuleSet],
        plugins: &[Arc<dyn crate::kernel::plugin::Plugin>],
        units: SessionUnits,
    ) -> Result<Self, RuleError> {
        let mut graph = Graph::new();
        let mut engine = Engine::install(&mut graph, sets)?;
        let mut units = units;
        for plugin in plugins {
            let manifest = plugin.manifest();
            {
                // The registry writes a plugin's units into this session's
                // registry, not into the table every process shares, so the
                // manifest built below describes them here and nowhere else.
                let mut registry = crate::kernel::plugin::Registry::with_units(
                    &mut graph,
                    engine.program_mut(),
                    &manifest.name,
                    &mut units,
                );
                plugin.register(&mut registry).map_err(|error| {
                    RuleError::Invalid {
                        rule: manifest.name.clone(),
                        reason: leak_reason(&error.to_string()),
                    }
                })?;
            }
            engine.add_plugin(manifest);
        }
        // The session's registry must be in place before the manifest that
        // describes it is built.
        engine.set_unit_registry(units);
        engine.refresh_manifest(&graph);
        Ok(Self {
            sets: sets.iter().map(|s| Arc::from(s.name())).collect(),
            rule_sets: sets.to_vec(),
            graph,
            engine,
        })
    }
}

/// A `'static` message for [`RuleError::Invalid`], which stores a static
/// string. The message is interned once per distinct text.
fn leak_reason(message: &str) -> &'static str {
    use std::collections::HashSet;
    use std::sync::Mutex;
    use std::sync::OnceLock;
    static SEEN: OnceLock<Mutex<HashSet<&'static str>>> = OnceLock::new();
    let seen = SEEN.get_or_init(|| Mutex::new(HashSet::new()));
    if let Ok(mut set) = seen.lock() {
        if let Some(existing) = set.get(message) {
            return existing;
        }
        let leaked: &'static str = Box::leak(message.to_owned().into_boxed_str());
        set.insert(leaked);
        return leaked;
    }
    "plugin registration failed"
}

/// Hands out session ids.
fn next_session_id() -> crate::kernel::value::SessionId {
    use std::sync::atomic::AtomicU64;
    use std::sync::atomic::Ordering;
    static NEXT: AtomicU64 = AtomicU64::new(1);
    crate::kernel::value::SessionId(NEXT.fetch_add(1, Ordering::Relaxed))
}

#[cold]
fn broken_builtin(error: &RuleError) -> ! {
    panic!("rssn: built-in rule sets failed to install: {error}")
}

impl<'s> Term<'s> {
    /// The session this term belongs to.
    #[must_use]
    pub const fn session(self) -> &'s Session {
        self.session
    }

    /// The literal number this term is, if it is one.
    #[must_use]
    pub fn as_number(self) -> Option<Number> {
        self.session
            .inner
            .borrow()
            .graph
            .as_number(self.node)
            .cloned()
    }

    /// The inferred dimension of the term in this session's graph (without
    /// the assumptions of a particular [`Config`]; use
    /// [`Answer::dimension`] for a configured computation).
    #[must_use]
    pub fn dimension(self) -> crate::graph::Dimension {
        self.session
            .inner
            .borrow()
            .graph
            .dimension(self.node)
            .dimension
    }

    /// The inferred dimension together with the conflicts found.
    #[must_use]
    pub fn dimension_report(self) -> crate::graph::DimOutcome {
        self.session.inner.borrow().graph.dimension(self.node)
    }

    /// A stable handle for this term within its session, for language
    /// bindings that cannot hold a borrowed [`Term`]; turn it back into a
    /// term with [`Session::term_by_id`].
    #[must_use]
    pub const fn id(self) -> u32 {
        self.node.raw()
    }

    /// The term as a dense numeric array: a literal number is a scalar
    /// (empty shape), a `list` of numbers a vector, a list of equally long
    /// lists a matrix, and so on. `None` when an entry is not a number or
    /// the nesting is ragged. Data is row-major.
    #[must_use]
    pub fn as_tensor(self) -> Option<(Vec<usize>, Vec<f64>)> {
        fn walk(
            graph: &Graph,
            node: NodeId,
            depth: usize,
            shape: &mut Vec<usize>,
            data: &mut Vec<f64>,
        ) -> Option<()> {
            if let Some(n) = graph.as_number(node) {
                return (depth == shape.len()).then(|| data.push(n.to_f64()));
            }
            if graph.op(node) != core::LIST {
                return None;
            }
            let items = graph.children(node);
            match shape.get(depth) {
                | Some(&len) if len != items.len() => return None,
                | Some(_) => {},
                | None if depth == shape.len() && data.is_empty() => shape.push(items.len()),
                | None => return None,
            }
            for &item in items {
                walk(graph, item, depth.saturating_add(1), shape, data)?;
            }
            Some(())
        }
        let inner = self.session.inner.borrow();
        let mut shape = Vec::new();
        let mut data = Vec::new();
        walk(&inner.graph, self.node, 0, &mut shape, &mut data)?;
        Some((shape, data))
    }

    /// The term as a float if it is a literal number.
    #[must_use]
    pub fn as_f64(self) -> Option<f64> {
        self.as_number().map(|n| n.to_f64())
    }

    /// The term as LaTeX math (see [`crate::io::latex::to_latex`]).
    #[must_use]
    pub fn to_latex(self) -> String {
        crate::io::latex::to_latex(&self.session.inner.borrow().graph, self.node)
    }

    /// The term as Typst math (see [`crate::io::typst::to_typst`]).
    #[must_use]
    pub fn to_typst(self) -> String {
        crate::io::typst::to_typst(&self.session.inner.borrow().graph, self.node)
    }

    /// The term as a multi-line Unicode drawing (see
    /// [`crate::io::pretty::to_pretty`]).
    #[must_use]
    pub fn to_pretty(self) -> String {
        crate::io::pretty::to_pretty(&self.session.inner.borrow().graph, self.node)
    }

    /// Evaluates the term numerically under the given bindings, without
    /// running any rules.
    #[must_use]
    pub fn eval(
        self,
        bindings: &[(&str, f64)],
    ) -> Option<f64> {
        let mut inner = self.session.inner.borrow_mut();
        let mut env = Env::numeric(0.0);
        for (name, value) in bindings {
            env.bind(inner.graph.interner_mut().symbol(name), *value);
        }
        inner.graph.eval(self.node, &env)
    }

    /// Compiles the term to a function of the named symbols, in order,
    /// with the [current](crate::backend::current) backend (the reference
    /// [`Interpreter`](crate::backend::Interpreter) unless another one was
    /// selected).
    ///
    /// # Errors
    /// See [`BackendError`](crate::backend::BackendError).
    pub fn compile(
        self,
        inputs: &[&str],
    ) -> Result<Box<dyn crate::backend::Compiled>, crate::backend::BackendError> {
        self.compile_with(&*crate::backend::current(), inputs)
    }

    /// Compiles the term with a specific backend.
    ///
    /// # Errors
    /// See [`BackendError`](crate::backend::BackendError).
    pub fn compile_with(
        self,
        backend: &dyn crate::backend::Backend,
        inputs: &[&str],
    ) -> Result<Box<dyn crate::backend::Compiled>, crate::backend::BackendError> {
        let mut inner = self.session.inner.borrow_mut();
        let symbols: Vec<_> = inputs
            .iter()
            .map(|name| inner.graph.interner_mut().symbol(name))
            .collect();
        backend.compile(&inner.graph, self.node, &symbols)
    }

    /// The term lowered to the [tape IR](crate::backend::tape) as a function
    /// of the named symbols, in order, and optimised: the common input of
    /// the compiling backends (`backend::jit`, `backend::lir`, devices).
    ///
    /// # Errors
    /// See [`BackendError`](crate::backend::BackendError).
    pub fn tape(
        self,
        inputs: &[&str],
    ) -> Result<crate::backend::tape::Tape, crate::backend::BackendError> {
        let mut inner = self.session.inner.borrow_mut();
        let symbols: Vec<_> = inputs
            .iter()
            .map(|name| inner.graph.interner_mut().symbol(name))
            .collect();
        crate::backend::tape::lower(&inner.graph, &[self.node], &symbols, &[]).map(|t| t.optimise())
    }

    /// Applies a unary operator that the standard rule sets register.
    fn apply(
        self,
        op: &str,
        rest: &[Self],
    ) -> Self {
        let mut args = Vec::with_capacity(rest.len().saturating_add(1));
        args.push(self);
        args.extend_from_slice(rest);
        match self.session.call(op, &args) {
            | Ok(term) => term,
            | Err(e) => missing_operator(&e),
        }
    }

    fn binary(
        self,
        op: crate::graph::OpId,
        other: Self,
    ) -> Self {
        let node = self
            .session
            .inner
            .borrow_mut()
            .graph
            .node(op, &[self.node, other.node]);
        self.session.term(node)
    }

    /// `self ^ exponent`.
    #[must_use]
    pub fn pow(
        self,
        exponent: impl IntoTerm<'s>,
    ) -> Self {
        self.binary(core::POW, exponent.into_term(self.session))
    }
}

#[cold]
fn missing_operator(error: &ComputeError) -> ! {
    panic!("rssn: {error}; the session was created without the rule set that defines it")
}

impl PartialEq for Term<'_> {
    /// Structural identity: the same term of the same session.
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        std::ptr::eq(self.session, other.session) && self.node == other.node
    }
}

impl Eq for Term<'_> {}

impl fmt::Display for Term<'_> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&self.session.inner.borrow().graph.display(self.node))
    }
}

impl fmt::Debug for Term<'_> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "Term({self})")
    }
}

/// Values that can stand where a [`Term`] is expected: terms themselves,
/// `i32` (exact integers) and `f64`.
///
/// Only one integer type is accepted on purpose, so that a bare literal in
/// `2 * x + 1` infers without annotations. Larger integers go through
/// [`Session::int`].
pub trait IntoTerm<'s> {
    /// Converts into a term of `session`.
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s>;
}

impl<'s> IntoTerm<'s> for Term<'s> {
    fn into_term(
        self,
        _session: &'s Session,
    ) -> Self {
        self
    }
}

impl<'s> IntoTerm<'s> for i32 {
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s> {
        session.int(i64::from(self))
    }
}

impl<'s> IntoTerm<'s> for f64 {
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s> {
        session.float(self)
    }
}

macro_rules! arithmetic {
    ($trait:ident, $method:ident, |$a:ident, $b:ident| $body:expr) => {
        impl<'s, R: IntoTerm<'s>> $trait<R> for Term<'s> {
            type Output = Term<'s>;

            fn $method(
                self,
                rhs: R,
            ) -> Term<'s> {
                let $a = self;
                let $b = rhs.into_term(self.session);
                $body
            }
        }
    };
    (scalar $scalar:ty) => {
        impl<'s> Add<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn add(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) + rhs
            }
        }

        impl<'s> Sub<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn sub(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) - rhs
            }
        }

        impl<'s> Mul<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn mul(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) * rhs
            }
        }

        impl<'s> Div<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn div(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) / rhs
            }
        }
    };
}

arithmetic!(Add, add, |a, b| a.binary(core::ADD, b));
arithmetic!(Mul, mul, |a, b| a.binary(core::MUL, b));
arithmetic!(Sub, sub, |a, b| a.binary(core::ADD, -b));
arithmetic!(Div, div, |a, b| a.binary(core::MUL, b.pow(-1)));
arithmetic!(scalar i32);
arithmetic!(scalar f64);

impl Neg for Term<'_> {
    type Output = Self;

    fn neg(self) -> Self {
        self.session.int(-1).binary(core::MUL, self)
    }
}

macro_rules! functions {
    ($($(#[$doc:meta])* $name:ident),* $(,)?) => {
        $(
            $(#[$doc])*
            #[must_use]
            pub fn $name(x: Term<'_>) -> Term<'_> {
                x.apply(stringify!($name), &[])
            }
        )*
    };
}

functions! {
    /// Exponential function.
    exp,
    /// Natural logarithm.
    ln,
    /// Sine.
    sin,
    /// Cosine.
    cos,
    /// Tangent.
    tan,
    /// Inverse sine.
    asin,
    /// Inverse cosine.
    acos,
    /// Inverse tangent.
    atan,
    /// Hyperbolic sine.
    sinh,
    /// Hyperbolic cosine.
    cosh,
    /// Hyperbolic tangent.
    tanh,
    /// Square root.
    sqrt,
    /// Absolute value.
    abs,
}

/// The derivative of `f` with respect to the symbol `x`, as an unevaluated
/// request.
#[must_use]
pub fn diff<'s>(
    f: Term<'s>,
    x: Term<'s>,
) -> Term<'s> {
    f.apply("diff", &[x])
}

/// The `n`-th derivative of `f` with respect to `x`.
#[must_use]
pub fn diff_n<'s>(
    f: Term<'s>,
    x: Term<'s>,
    n: usize,
) -> Term<'s> {
    (0..n).fold(f, |acc, _| diff(acc, x))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `compute_all` equals a loop of `compute` (same answers, same store
    /// state), also with a failing request in the middle, and does not
    /// depend on the thread count.
    #[test]
    fn compute_all_equals_a_loop_of_compute() {
        let texts = [
            "diff(sin(x)^2*exp(x^2), x)",
            "(x + 1)^4 - (x - 1)^4",
            "diff(x^3 + 2*x, x)",
            "integral(x*exp(x), x)",
        ];
        let run = |batch: bool| -> Vec<String> {
            let s = Session::new();
            let terms: Vec<_> = texts.iter().map(|t| s.parse(t).unwrap()).collect();
            let config = Config::new();
            let answers: Vec<_> = if batch {
                s.compute_all(&terms, &config)
            } else {
                terms.iter().map(|&t| s.compute(t, &config)).collect()
            };
            let mut out: Vec<String> = answers
                .into_iter()
                .map(|a| {
                    let a = a.unwrap();
                    format!("{} {:?} {:?} {}", a.term, a.value, a.error, a.reduced)
                })
                .collect();
            out.push(format!("{s:?} id {}", s.sym("zz").id()));
            out
        };
        let looped = run(false);
        assert_eq!(run(true), looped);
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        assert_eq!(one.install(|| run(true)), looped);
    }

    #[test]
    fn compute_all_reports_errors_per_request() {
        let s = Session::new();
        let t = s.parse("x + 1").unwrap();
        let u = s.parse("sin(x)").unwrap();
        // A numeric target without a binding for x has no value.
        let config = Config::new().numeric(1e-9);
        let batch = s.compute_all(&[t, u], &config);
        assert!(matches!(batch[0], Err(ComputeError::NotNumeric(_))));
        assert!(matches!(batch[1], Err(ComputeError::NotNumeric(_))));
        assert!(s.compute_all(&[], &config).is_empty());
    }

    fn symbolic(
        s: &Session,
        text: &str,
    ) -> String {
        let term = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
        let answer = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            answer.reduced,
            "`{text}` was not fully reduced: {}",
            answer.term
        );
        answer.term.to_string()
    }

    /// Checks a derivative in both phases: the symbolic answer, evaluated,
    /// must agree with the numeric answer and with a finite difference of
    /// the original function.
    fn check_derivative(
        text: &str,
        at: f64,
    ) {
        let s = Session::new();
        let f = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
        let x = s.sym("x");
        let request = diff(f, x);
        let closed = s
            .compute(request, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            closed.reduced,
            "d/dx {text} was not reduced: {}",
            closed.term
        );
        let from_closed = closed.term.eval(&[("x", at)]).unwrap_or(f64::NAN);
        let numeric = s
            .compute(request, &Config::new().numeric(1e-9).bind("x", at))
            .unwrap_or_else(|e| panic!("{e}"))
            .value
            .unwrap_or(f64::NAN);
        let h = 1e-5;
        let fd = (f.eval(&[("x", at + h)]).unwrap_or(f64::NAN)
            - f.eval(&[("x", at - h)]).unwrap_or(f64::NAN))
            / (2.0 * h);
        let scale = 1.0 + fd.abs();
        assert!(
            (from_closed - numeric).abs() < 1e-9 * scale,
            "{text}: closed form {from_closed} vs numeric {numeric}"
        );
        assert!(
            (from_closed - fd).abs() < 1e-6 * scale,
            "{text}: {} = {from_closed} but finite difference = {fd}",
            closed.term
        );
    }

    #[test]
    fn operators_build_canonical_terms() {
        let s = Session::new();
        let (x, y) = (s.sym("x"), s.sym("y"));
        assert_eq!((x + y).to_string(), "x + y");
        assert_eq!((x - y).to_string(), "x - y");
        assert_eq!((x / y).to_string(), "x/y");
        assert_eq!((2 * x + 1).to_string(), "2*x + 1");
        assert_eq!((-x).to_string(), "-x");
        assert_eq!(x.pow(2).to_string(), "x^2");
        assert_eq!((1.5 * x).to_string(), "1.5*x");
        assert_eq!(sin(x + y), s.parse("sin(y + x)").unwrap_or(x));
    }

    #[test]
    fn simplification() {
        let s = Session::new();
        assert_eq!(symbolic(&s, "x + x + sin(y)^2 + cos(y)^2"), "2*x + 1");
        assert_eq!(symbolic(&s, "(x^2 * x^3) / x^4"), "x");
        assert_eq!(symbolic(&s, "exp(ln(a)) - a"), "0");
    }

    #[test]
    fn derivatives_symbolic() {
        let s = Session::new();
        assert_eq!(symbolic(&s, "diff(x^3, x)"), "3*x^2");
        assert_eq!(symbolic(&s, "diff(sin(x), x)"), "cos(x)");
        assert_eq!(symbolic(&s, "diff(a*x^2 + b*x + c, x)"), "2*a*x + b");
        assert_eq!(symbolic(&s, "diff(exp(2*x), x)"), "2*exp(2*x)");
        assert_eq!(symbolic(&s, "diff(ln(x), x)"), "1/x");
        assert_eq!(symbolic(&s, "diff(y, x)"), "0");
        assert_eq!(symbolic(&s, "diff(diff(x^4, x), x)"), "12*x^2");
        assert_eq!(symbolic(&s, "diff(sin(x)^2 + cos(x)^2, x)"), "0");
    }

    #[test]
    fn derivatives_agree_across_phases() {
        for (text, at) in [
            ("x^3 - 2*x + 1", 0.7),
            ("sin(x) * cos(x)", 0.3),
            ("exp(x^2) / (1 + x^2)", 0.5),
            ("ln(1 + x^2) * tan(x)", 0.4),
            ("x^x", 1.3),
            ("sqrt(1 + sin(x)^2)", 1.1),
            ("atan(x) + asin(x/2) + acos(x/3)", 0.6),
            ("sinh(x) * cosh(2*x) - tanh(x)", 0.2),
            ("(x + 1)^5 * (x - 2)^(-3)", 0.9),
            ("2^x + x^(1/3)", 1.5),
            ("abs(x) * x", -0.8),
            ("cot(x) + sec(x) * csc(x)", 0.7),
            ("acot(x) + asec(x + 2) + acsc(x + 2)", 0.6),
            ("coth(x) - sech(x) * csch(x)", 0.9),
            ("asinh(x) + acosh(x + 2) + atanh(x / 2)", 0.5),
            ("acoth(x + 2) + asech(x / 2) + acsch(x)", 0.8),
            ("atan2(x, 2) + log(3, x)", 1.4),
            ("x * y * sin(x * y)", 0.35),
        ] {
            if text.contains('y') {
                // `y` is an independent symbol: bind it for evaluation.
                let s = Session::new();
                let f = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
                let request = diff(f, s.sym("x"));
                let closed = s
                    .compute(request, &Config::new())
                    .unwrap_or_else(|e| panic!("{e}"));
                let got = closed
                    .term
                    .eval(&[("x", at), ("y", 2.0)])
                    .unwrap_or(f64::NAN);
                let want = 2.0 * (2.0 * at).sin() + 2.0 * at * 2.0 * (2.0 * at).cos();
                assert!(
                    (got - want).abs() < 1e-12,
                    "{}: {got} vs {want}",
                    closed.term
                );
            } else {
                check_derivative(text, at);
            }
        }
    }

    #[test]
    fn derivatives_of_undetermined_functions_are_closed_forms() {
        let s = Session::new();
        let term = s
            .parse("diff(x * f(x), x)")
            .unwrap_or_else(|e| panic!("{e}"));
        let answer = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        // `diff(f(x), x)` cannot be computed any further: it is the answer.
        assert!(answer.reduced);
        assert_eq!(answer.term.to_string(), "x*diff(f(x), x) + f(x)");
    }

    #[test]
    fn numeric_target_without_bindings_fails_cleanly() {
        let s = Session::new();
        let term = s.parse("x + 1").unwrap_or_else(|e| panic!("{e}"));
        let result = s.compute(term, &Config::new().numeric(1e-9));
        assert_eq!(
            result.err(),
            Some(ComputeError::NotNumeric("x + 1".to_owned()))
        );
    }

    #[test]
    fn numeric_differentiation_of_black_box_operators() {
        use crate::graph::Arity;
        use crate::graph::OpDescriptor;
        // An operator with scalar semantics but no derivative rule.
        let blackbox = RuleSet::new("blackbox", |i| {
            i.op(OpDescriptor::new("bump", Arity::Fixed(1))
                .eval(|a| a.first().map_or(f64::NAN, |x| (-x * x).exp())))
                .map(|_| ())
        });
        let config = Config::new().with(blackbox.clone());
        let s =
            Session::with_rules(&[rules::calculus(), blackbox]).unwrap_or_else(|e| panic!("{e}"));
        let term = s
            .parse("diff(bump(x), x)")
            .unwrap_or_else(|e| panic!("{e}"));

        let symbolic = s.compute(term, &config).unwrap_or_else(|e| panic!("{e}"));
        assert!(!symbolic.reduced, "there is no symbolic rule for bump");

        let numeric = s
            .compute(term, &config.numeric(1e-6).bind("x", 0.5))
            .unwrap_or_else(|e| panic!("{e}"));
        let want = -2.0 * 0.5 * (-0.25_f64).exp();
        let got = numeric.value.unwrap_or(f64::NAN);
        assert!((got - want).abs() < 1e-8, "{got} vs {want}");
        assert!(numeric.error.is_some_and(|e| e < 1e-6));
    }

    #[test]
    fn sessions_are_reusable_and_bindings_do_not_stick() {
        let s = Session::new();
        let term = s
            .parse("diff(x^2, x) + y")
            .unwrap_or_else(|e| panic!("{e}"));
        let a = s.compute(
            term,
            &Config::new().numeric(1e-9).bind("x", 1.0).bind("y", 10.0),
        );
        let b = s.compute(
            term,
            &Config::new().numeric(1e-9).bind("x", 2.0).bind("y", 20.0),
        );
        assert_eq!(a.ok().and_then(|a| a.value), Some(12.0));
        assert_eq!(b.ok().and_then(|b| b.value), Some(24.0));
        assert_eq!(symbolic(&s, "diff(x^2, x) + y"), "2*x + y");
    }

    #[test]
    fn assumptions_are_per_request() {
        let s = Session::new();
        let term = s.parse("sqrt(x^2)").unwrap_or_else(|e| panic!("{e}"));
        let plain = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let positive = s
            .compute(term, &Config::new().assume("x", Facts::POSITIVE))
            .unwrap_or_else(|e| panic!("{e}"));
        let again = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(plain.term.to_string(), "sqrt(x^2)");
        assert_eq!(positive.term.to_string(), "x");
        assert_eq!(
            again.term.to_string(),
            "sqrt(x^2)",
            "an assumption must not outlive its request"
        );
    }

    #[test]
    fn reduce_then_compile() {
        let s = Session::new();
        let x = s.sym("x");
        let derivative = s
            .simplify(diff(x.pow(3) * sin(x), x))
            .unwrap_or_else(|e| panic!("{e}"));
        let f = derivative.compile(&["x"]).unwrap_or_else(|e| panic!("{e}"));
        let xs: Vec<f64> = (0..100).map(|i| f64::from(i) * 0.1).collect();
        let mut out = vec![0.0; xs.len()];
        f.call_batch(&[&xs], &mut out);
        for (x, got) in xs.iter().zip(&out) {
            let want = 3.0 * x * x * x.sin() + x.powi(3) * x.cos();
            assert!(
                (got - want).abs() < 1e-9 * (1.0 + want.abs()),
                "at {x}: {got} vs {want}"
            );
        }
        // An unreduced request has no numeric semantics.
        assert!(diff(x, x).compile(&["x"]).is_err());
    }

    #[test]
    fn unknown_operators_are_reported() {
        let s = Session::new();
        let x = s.sym("x");
        assert_eq!(
            s.call("frobnicate", &[x]).err(),
            Some(ComputeError::UnknownOperator("frobnicate".to_owned()))
        );
        assert!(s.call("sin", &[x, x]).is_err());
        assert!(s.parse("1 +").is_err());
    }

    /// A unit defined by a session, not by the process.
    fn furlong() -> UnitDef {
        UnitDef::linear(
            "furlong",
            crate::kernel::signature::DimExpr::base("Length"),
            Scale::one(),
        )
    }

    fn names(manifest: &KernelManifest) -> Vec<String> {
        manifest
            .units
            .iter()
            .map(|unit| unit.name.clone())
            .collect()
    }

    /// A unit registered on one session is described by that session's
    /// manifest — and by no other session's, which is the whole point of the
    /// session owning its registry instead of sharing the global table.
    #[test]
    fn a_session_unit_stays_in_its_own_manifest() {
        let mut mine = Session::new();
        mine.register_unit(furlong())
            .unwrap_or_else(|e| panic!("{e}"));
        let theirs = Session::new();

        let mine = mine.manifest();
        let theirs = theirs.manifest();
        assert!(
            names(&mine).iter().any(|name| name == "furlong"),
            "the session that registered it must describe it"
        );
        assert!(
            !names(&theirs).iter().any(|name| name == "furlong"),
            "a unit must not leak into another session's manifest"
        );
        // Both still describe the standard units.
        for manifest in [&mine, &theirs] {
            assert!(names(manifest).iter().any(|name| name == "m"));
            assert!(names(manifest).iter().any(|name| name == "degC"));
            assert!(manifest.verify(), "the hash must match the contents");
        }
        assert_ne!(mine.manifest_hash, theirs.manifest_hash);
    }

    /// The units `spec/05-units.md` §5.3.3 declares, and the registry used to
    /// lack, reach the manifest — which is what a host generates `Units.SI`
    /// from — so the count and the exact scales are visible there.
    #[test]
    fn the_manifest_describes_the_common_non_si_units() {
        let manifest = Session::new().manifest();
        // 75 built-in units before this round plus the seven §5.3.3 units.
        assert_eq!(
            manifest.units.len(),
            82,
            "the manifest must carry every built-in unit"
        );
        for name in ["Da", "kn", "mph", "acre", "bit", "B", "KiB"] {
            let unit = manifest
                .units
                .iter()
                .find(|unit| unit.name == name)
                .unwrap_or_else(|| panic!("the manifest has no unit `{name}`"));
            let definition = crate::dim::units::lookup(name)
                .unwrap_or_else(|| panic!("the registry has no unit `{name}`"));
            assert_eq!(
                unit.scale.numer,
                definition.scale.rat.numer().clone(),
                "`{name}` has a different numerator in the manifest"
            );
            assert_eq!(
                unit.scale.denom,
                definition.scale.rat.denom().clone(),
                "`{name}` has a different denominator in the manifest"
            );
            assert_eq!(
                unit.scale.pi, definition.scale.pi,
                "`{name}` has a different power of pi in the manifest"
            );
        }
        assert!(manifest.verify(), "the hash must match the contents");
    }

    /// A second registration keeps the first: the registry a template's
    /// engine holds is carried forward, not rebuilt from the global table.
    #[test]
    fn session_units_accumulate() {
        let mut s = Session::new();
        s.register_unit(furlong()).unwrap_or_else(|e| panic!("{e}"));
        s.register_unit(UnitDef::linear(
            "stone",
            crate::kernel::signature::DimExpr::base("Mass"),
            Scale::one(),
        ))
        .unwrap_or_else(|e| panic!("{e}"));
        let manifest = s.manifest();
        assert!(names(&manifest).iter().any(|name| name == "furlong"));
        assert!(names(&manifest).iter().any(|name| name == "stone"));
        assert!(names(&manifest).iter().any(|name| name == "m"));
    }

    /// A plugin that registers one unit, `cubit`, scaled by `tenths / 10 m`.
    /// Two instances are the same plugin — same name, same unit name — with
    /// different definitions, which is exactly what two sessions must not
    /// share.
    fn cubit_plugin(tenths: i64) -> Box<dyn crate::kernel::plugin::Plugin> {
        Box::new(crate::kernel::plugin::SimplePlugin::new(
            "acme/cubit",
            "0.1.0",
            move |reg| {
                reg.register_unit(UnitDef::linear(
                    "cubit",
                    crate::kernel::signature::DimExpr::base("Length"),
                    Scale {
                        rat: BigRational::new(tenths.into(), 10.into()),
                        pi: 0,
                    },
                ))?;
                Ok(())
            },
        ))
    }

    /// The scale of `cubit` as the manifest describes it.
    fn cubit_scale(manifest: &KernelManifest) -> Option<String> {
        manifest
            .units
            .iter()
            .find(|unit| unit.name == "cubit")
            .map(|unit| format!("{}/{}", unit.scale.numer, unit.scale.denom))
    }

    /// The unit a plugin registers belongs to the session that installed the
    /// plugin. Two sessions install the same plugin with different
    /// definitions of the same unit; each manifest must report its own, and
    /// a third session that never installed the plugin must report neither —
    /// that last check is what proves the unit did not reach the process-wide
    /// table.
    #[test]
    fn a_plugin_unit_stays_in_the_session_that_installed_it() {
        let mut half = Session::new();
        half.register(cubit_plugin(5))
            .unwrap_or_else(|e| panic!("{e}"));
        let mut tenths = Session::new();
        tenths
            .register(cubit_plugin(7))
            .unwrap_or_else(|e| panic!("{e}"));
        let neither = Session::new();

        let half = half.manifest();
        let tenths = tenths.manifest();
        let neither = neither.manifest();
        assert_eq!(cubit_scale(&half), Some(String::from("1/2")));
        assert_eq!(cubit_scale(&tenths), Some(String::from("7/10")));
        assert_eq!(
            cubit_scale(&neither),
            None,
            "a plugin unit must not leak into a session that did not install the plugin"
        );
        assert!(
            crate::dim::units::lookup("cubit").is_none(),
            "a plugin unit must not reach the table every process shares"
        );
        for manifest in [&half, &tenths] {
            assert_eq!(
                manifest.units.iter().filter(|u| u.name == "cubit").count(),
                1,
                "each session describes its own definition exactly once"
            );
            assert!(manifest.verify(), "the hash must match the contents");
        }
        assert_eq!(
            neither.units.iter().filter(|u| u.name == "cubit").count(),
            0
        );
        assert!(neither.verify(), "the hash must match the contents");
        assert_ne!(half.manifest_hash, tenths.manifest_hash);
    }
}
