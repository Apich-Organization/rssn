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

pub use extract::ClosedForm;
pub use extract::CostModel;
pub use extract::Extractor;
pub use extract::SizeCost;
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
