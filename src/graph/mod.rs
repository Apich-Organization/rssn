//! # The graph kernel
//!
//! The domain-free core of rssn: a hash-consed expression DAG whose nodes
//! are grouped into equivalence classes, together with the machinery that
//! discovers equalities (rules), explores cheaply inside tree-shaped regions
//! (windows), schedules the search heuristically and reads a result back out
//! (extraction).
//!
//! Nothing in here knows what a derivative or a matrix is. Domains live in
//! [`crate::rules`] and talk to the kernel only through operators registered
//! in the [`OpTable`] and rules grouped into rule sets.
//!
//! # Limitations and coverage
//!
//! Priorities: P1 = correctness/soundness, P2 = capability, P3 = polish.
//!
//! * **Witness balls are estimates, not certificates (P1).** A [`Ball`] is a
//!   `mid ± rad` witness, not a proof: numeric kernels behind [`Outcome::Approx`]
//!   report convergence estimates, and witness propagation is outward rounded
//!   only for sums and products (every other operator inherits the largest
//!   child radius). Witnesses feed reporting, [`Graph::approx`] and the
//!   [`Evaluated`] goal, never facts.
//! * **Sign facts need a certificate (P1).** A constant's sign comes only from
//!   exact literals, structure or an outward-rounded tape enclosure
//!   ([`Graph::approx_interval`], [`crate::backend::tape::Tape::enclose`]).
//!   The enclosure covers `+ - * /`, integer and positive-base real powers and
//!   most elementary functions, but not `tan`, `atan2`, opaque
//!   `Call`/`Extern`/`Gather`; it is taken over the `f64` image of literals,
//!   `sin`/`cos` degrade to `[-1, 1]` beyond `2π`, and boxes are not split at
//!   poles or branch cuts. Many constants therefore stay sign-unknown.
//! * **Extraction is a local greedy optimum (P2).** [`SizeCost`] and
//!   [`ClosedForm`] price shape; [`EvalCost`] prices the tape instructions a
//!   representative lowers to (weights in `extract.rs`, cached per node and
//!   class). All are additive local proxies: no common-subexpression credit,
//!   no register pressure or vectorisation, only reachable classes, and the
//!   optimum is local to the classes the rules explored.
//! * **Search is budgeted and incomplete (P2).** The defaults (24 iterations,
//!   20k nodes, 1000 matches per rewrite, 512 window steps, plus match fuel,
//!   bans and plateau back-off) make [`Saturate`] mean "no change within the
//!   budget", not "all consequences derived". Without the wall-clock
//!   `Budget::time_limit`, a fixed graph, program and budget are deterministic.
//! * **Windows see only unshared trees (P3).** A [`TreeWindow`] copies at most
//!   `window_cells` cells of the unshared region; shared subterms stay opaque
//!   atoms, greedy/beam search is not exhaustive, and equalities discarded
//!   inside a window are lost with it.
//! * **Soundness sampling is a smoke test (P2).** Rewrites are checked on
//!   name-seeded random substitutions in `[-3, 3]` (positive requirements
//!   `(0.25, 3)`), and every fourth sample is drawn near a singular locus
//!   (zero, the log cut, `±1`, `±π/2`) or a cancellation. Such a sample is
//!   compared with a tolerance scaled by intermediate magnitude and
//!   finite-difference sensitivity, is inconclusive when evaluation leaves the
//!   reals, and only refutes if the disagreement survives a small absolute
//!   jitter (a rounding artefact at a discontinuity does not). Together these
//!   suppress rounding/branch false alarms at the cost of refutation power
//!   exactly where conditioning is worst; integer and rational magnitudes are
//!   also clamped so special functions stay cheap. Agreement is still not a
//!   proof, and opaque operators make a rule inconclusive.
//! * **Tape derivative conventions (P2).** Reverse gradients/VJPs and
//!   Hessian-vector products use linear tape storage; the legacy dense Hessian
//!   uses forward passes. Third derivatives are not implemented. Opaque
//!   `Call`/`Extern` first derivatives are differenced and their second
//!   derivatives remain NaN; `abs`/`floor`/`ceil` use zero second derivative at
//!   their kinks.

pub mod dim;
/// Polymorphic dimension signatures of plugin operators.
pub mod dim_signature;
// `pub(crate)` is the intended visibility: these are internal modules, so
// `pub` items inside them are exactly as crate-private. Keep the explicit
// `pub(crate)`: no extra visibility attribute is needed on the modules.

pub(crate) mod dim_expr;
pub(crate) mod dim_infer;
pub mod eval;
pub mod extract;
pub mod facts;
pub mod id;
pub mod number;
pub mod op;
pub mod pattern;
pub mod payload;
pub mod print;
#[cfg(test)]
mod proptests;
pub mod rule;
pub mod schedule;
pub mod soundness;
pub mod store;
pub mod subst;
pub mod window;

pub use dim::BaseId;
pub use dim::DimError;
pub use dim::DimOutcome;
pub use dim::DimRule;
pub use dim::DimTerm;
pub use dim::DimVarId;
pub use dim::Dimension;
pub use dim::DimensionMode;
pub use dim::FamilySource;
pub use dim::UnitInfo;
pub use dim::UnitResolver;
pub use dim::UnitTextResolver;
pub use dim_signature::DimSignature;
pub use dim_signature::signature_rule;
pub use eval::ComplexEval;
pub use extract::ClosedForm;
pub use extract::CostModel;
pub use extract::EvalCost;
pub use extract::Extractor;
pub use extract::SizeCost;
pub use facts::Always;
pub use facts::Facts;
pub use facts::OnReals;
pub use id::ClassId;
pub use id::NodeId;
pub use id::OpId;
pub use id::PayloadId;
pub use id::SymbolId;
pub use number::Number;
pub use op::Arity;
pub use op::OpDescriptor;
pub use op::OpFlags;
pub use op::OpTable;
pub use pattern::Match;
pub use pattern::ParseError;
pub use pattern::Pat;
pub use pattern::VarNames;
pub use payload::Blob;
pub use payload::Payload;
pub use rule::Action;
pub use rule::Cx;
pub use rule::Env;
pub use rule::Guard;
pub use rule::Installer;
pub use rule::Kernel;
pub use rule::Outcome;
pub use rule::Prepared;
pub use rule::PreparedWork;
pub use rule::Program;
pub use rule::Rewrite;
pub use rule::Rule;
pub use rule::RuleError;
pub use rule::RuleSet;
pub use rule::Tier;
pub use schedule::Budget;
pub use schedule::Engine;
pub use schedule::Evaluated;
pub use schedule::Goal;
pub use schedule::Report;
pub use schedule::Saturate;
pub use schedule::Stop;
pub use store::Ball;
pub use store::Graph;
pub use window::CellId;
pub use window::TreeWindow;
pub use window::WindowMatch;
pub use window::WindowPass;
