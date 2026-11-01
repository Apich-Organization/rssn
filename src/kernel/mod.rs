//! The embedding kernel: everything a host needs to drive RSSN from outside
//! the crate.
//!
//! # Why this module exists
//!
//! A downstream project asked RSSN to be split into `rssn-api` (the dynamic
//! call ABI, the manifest, the plugin registry) and `rssn-dim` (the dimension
//! and unit layer). Splitting a crate that is mid-migration buys nothing here:
//! the two halves share the expression graph, the operator table and the
//! dimension arena, so a hard crate boundary would mean re-exporting half of
//! `graph` through a public facade and paying a versioning tax for the
//! privilege.
//!
//! So the *interfaces* live here instead, as ordinary modules with no
//! `pub(crate)` shortcuts, and the day a split becomes worthwhile the modules
//! can be lifted into crates without touching a call site. [`crate::dim`] is
//! the other half of that arrangement.
//!
//! # The pieces
//!
//! | Module | What it is |
//! |---|---|
//! | [`signature`] | The type language: [`Signature`], [`TypeExpr`], [`DimExpr`], [`OpKey`]. |
//! | [`value`] | The runtime value language: [`Value`] and its containers. |
//! | [`call`] | The calling convention: [`CallCtx`], [`OpError`], [`Kernel`]. |
//! | [`plugin`] | Extension points: [`Plugin`], [`HostCallback`], [`Registry`]. |
//! | [`manifest`] | The self-description a host reads: [`KernelManifest`]. |
//! | [`codec`] | Deterministic encoding, over `bincode-next` only. |
//! | [`infer`] | How an operator that never declared a signature gets one. |
//!
//! # The invariants
//!
//! * **One wire format.** [`codec`] is the only place that turns a value into
//!   bytes, and it always uses `bincode-next` — deterministic bincode for
//!   hashing and storage, deterministic CBOR for interchange. Nothing in RSSN
//!   writes JSON to a peer; JSON exists only for humans reading a manifest.
//! * **One description.** A [`Signature`] is declared once, hashed once, and
//!   every [`OpKey`], manifest entry and overload check derives from it.
//! * **No silent unknowns.** An operator either declares a signature, has one
//!   inferred, or records a [`DynamicReason`] explaining why it cannot. The
//!   `signature_audit` test fails the build otherwise.
//! * **Determinism.** A manifest built from the same operator table and rule
//!   sets is byte-identical, in any process, in any order of registration.
//!
//! # A worked example
//!
//! ```no_run
//! use rssn::kernel::{CallCtx, OpError, Signature, TypeExpr, Value};
//!
//! # fn run(session: &rssn::api::Session) -> Result<(), OpError> {
//! // What does the kernel claim to be?
//! let manifest = session.manifest();
//! println!("{} operators, hash {:016x}", manifest.ops.len(), manifest.hash_u64());
//!
//! // Call an operator by name, dynamically, with a checked argument list.
//! let value = session.call_value("sin", &[Value::Real(0.5)], CallCtx::new())?;
//!
//! // Or look up what it accepts before calling it.
//! if let Some(op) = manifest.op("sin") {
//!     let signature: &Signature = &op.signatures[0];
//!     assert_eq!(signature.params[0].ty, TypeExpr::Real);
//! }
//! # Ok(())
//! # }
//! ```

pub mod call;
pub mod codec;
pub mod infer;
pub mod manifest;
pub mod plugin;
pub mod result;
pub mod signature;
pub mod value;

pub use call::{Assumption, CallCtx, Kernel, OpError, Target, Trace, TraceEvent};
pub use codec::CodecError;
pub use manifest::{
    BaseDimInfo, Features, KernelManifest, OpInfo, PluginInfo, RuleSetInfo, SignatureSource,
};
pub use result::{
    Condition, ConditionKind, EventRecord, Field, FieldKind, FlattenSchema, OdeSolutionView,
    SchemaError, SimField, SimResult, Solutions,
};
pub use plugin::{
    HostCallback, NativeCallback, Plugin, PluginError, PluginId, PluginManifest, Registry,
    SimplePlugin,
};
pub use signature::{
    API_VERSION, PLUGIN_ABI, ArityInfo, BackendSupport, BinderInfo, Constraint, ConstraintKind,
    Determinism, DimExpr, DimRuleInfo, DynamicReason, Effect, EffectRow, EvalKind, Example, Exp,
    FamilySourceInfo, NatExpr, OpFlagsInfo, OpKey, Purity, Rational, RepExpr, RowExpr, Signature,
    Stability, ThreadSafety, TypeExpr, TypeParam, TypeParamKind, Variadic,
};
pub use value::{
    CallbackRef, ComplexValue, DynValue, OpaqueValue, QuantityValue, RecordValue, ScaleValue,
    SessionId, SparseFormat, SparseValue, TensorValue, TermRef, UnitValue, Value, VariantValue,
};
