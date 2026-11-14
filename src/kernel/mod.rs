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
//! | [`diag`] | Why an operator call did not reduce: [`Diag`], [`DiagKind`], `Session::apply_diag`. |
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
//! use rssn::kernel::CallCtx;
//! use rssn::kernel::OpError;
//! use rssn::kernel::Signature;
//! use rssn::kernel::TypeExpr;
//! use rssn::kernel::Value;
//!
//! # fn run(session: &rssn::api::Session) -> Result<(), OpError> {
//! // What does the kernel claim to be?
//! let manifest = session.manifest();
//! println!(
//!     "{} operators, hash {:016x}",
//!     manifest.ops.len(),
//!     manifest.hash_u64()
//! );
//!
//! // Call an operator by name, dynamically, with a checked argument list.
//! let value = session.call_value("sin", &[Value::Float(0.5)], &CallCtx::new())?;
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
mod callback;
pub use callback::CallbackBank;
pub mod cdylib;
pub mod codec;
pub mod diag;
pub mod infer;
pub mod manifest;
pub mod overload;
pub mod plugin;
pub mod result;
pub mod signature;
pub mod sort;
pub mod value;

pub use call::Assumption;
pub use call::CallCtx;
pub use call::Kernel;
pub use call::OpError;
pub use call::Target;
pub use call::Trace;
pub use call::TraceEvent;
pub use codec::CodecError;
pub use diag::Diag;
pub use diag::DiagKind;
pub use manifest::BaseDimInfo;
pub use manifest::Features;
pub use manifest::KernelManifest;
pub use manifest::OpInfo;
pub use manifest::PluginInfo;
pub use manifest::RuleSetInfo;
pub use manifest::SignatureSource;
pub use overload::Best;
pub use overload::Match;
pub use plugin::HostCallback;
pub use plugin::NativeCallback;
pub use plugin::Plugin;
pub use plugin::PluginError;
pub use plugin::PluginId;
pub use plugin::PluginKind;
pub use plugin::PluginManifest;
pub use plugin::Registry;
pub use plugin::SimplePlugin;
pub use result::Condition;
pub use result::ConditionKind;
pub use result::EventRecord;
pub use result::Field;
pub use result::FieldKind;
pub use result::FlattenSchema;
pub use result::OdeSolutionView;
pub use result::SchemaError;
pub use result::SimField;
pub use result::SimResult;
pub use result::Solutions;
pub use signature::API_VERSION;
pub use signature::ArityInfo;
pub use signature::BackendSupport;
pub use signature::BinderInfo;
pub use signature::Constraint;
pub use signature::ConstraintKind;
pub use signature::Determinism;
pub use signature::DimExpr;
pub use signature::DimRuleInfo;
pub use signature::DynamicReason;
pub use signature::Effect;
pub use signature::EffectRow;
pub use signature::EvalKind;
pub use signature::Example;
pub use signature::Exp;
pub use signature::FamilySourceInfo;
pub use signature::NatExpr;
pub use signature::OpFlagsInfo;
pub use signature::OpKey;
pub use signature::PLUGIN_ABI;
pub use signature::Purity;
pub use signature::Rational;
pub use signature::RepExpr;
pub use signature::RowExpr;
pub use signature::Signature;
pub use signature::Stability;
pub use signature::ThreadSafety;
pub use signature::TypeExpr;
pub use signature::TypeParam;
pub use signature::TypeParamKind;
pub use signature::Variadic;
pub use value::CallbackRef;
pub use value::ComplexValue;
pub use value::ConstantValue;
pub use value::DynValue;
pub use value::OpaqueValue;
pub use value::QuantityValue;
pub use value::RecordValue;
pub use value::ScaleValue;
pub use value::SessionId;
pub use value::SparseFormat;
pub use value::SparseValue;
pub use value::TensorValue;
pub use value::TermRef;
pub use value::UnitValue;
pub use value::Value;
pub use value::ValueKind;
pub use value::ValueKindInfo;
pub use value::VariantValue;
