//! The kernel manifest: everything a host needs to know about this build.
//!
//! A manifest is generated once per session configuration, from the operator
//! table that is actually installed. It lists every operator with its
//! signatures, documentation, stability and capabilities, every unit, every
//! base dimension, and every rule set — plus a hash of all of it.
//!
//! # Determinism
//!
//! The manifest is a *build artefact with a stable identity*. Two processes
//! that installed the same rule sets produce byte-identical manifests, and
//! therefore the same [`KernelManifest::manifest_hash`], whatever order the
//! rule sets were listed in and whatever the registration order was. That
//! property is what lets a host use the hash as a cache key:
//!
//! * every list is sorted by a total order on its contents;
//! * nothing carries a timestamp, a pointer, an arena-local id or a hash-map
//!   iteration order;
//! * the hash is BLAKE3 over the deterministic binary encoding
//!   ([`codec`](crate::kernel::codec)), with the hash field itself zeroed.
//!
//! # Size
//!
//! A full manifest for the standard rule sets is a few hundred kilobytes of
//! binary, and smaller as CBOR. A host that only needs part of it should
//! read the CBOR form and skip what it does not use: the encoding is
//! self-describing enough to do that without a full parse.

use std::collections::BTreeMap;

use serde::Deserialize;
use serde::Serialize;

use crate::graph::Graph;
use crate::graph::OpDescriptor;
use crate::graph::RuleSet;
use crate::kernel::codec;
use crate::kernel::signature::API_VERSION;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::BinderInfo;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DimRuleInfo;
use crate::kernel::signature::DynamicReason;
use crate::kernel::signature::EvalKind;
use crate::kernel::signature::Example;
use crate::kernel::signature::OpFlagsInfo;
use crate::kernel::signature::OpKey;
use crate::kernel::signature::Purity;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::value::UnitValue;

/// What this build can do.
///
/// A host reads this before generating bindings: a `false` here means the
/// corresponding capability is absent and calls that need it return
/// [`OpError::Unsupported`](crate::kernel::OpError::Unsupported), never a
/// silently wrong answer.
// A capability set is a list of independent facts; a bitfield would make
// the manifest JSON unreadable and the struct no smaller in practice.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Features {
    /// Operators can be called by name with [`Value`](crate::kernel::Value)
    /// arguments.
    pub dynamic_call: bool,
    /// Version of the plugin ABI.
    pub plugin_abi: u32,
    /// A JIT backend is compiled in.
    pub jit: bool,
    /// The JIT lowers `select` (branching).
    pub jit_branch: bool,
    /// The JIT lowers load/index (tensor indexing).
    pub jit_index: bool,
    /// Symbolic requests can be run in parallel.
    pub parallel_symbolic: bool,
    /// Unit scales are exact rationals, not floats.
    pub exact_units: bool,
    /// Affine units (`degC`, `degF`) are supported.
    pub affine_units: bool,
    /// Hosts may define base dimensions.
    pub user_base_dimensions: bool,
    /// Results can be returned as structured values (`Solutions`,
    /// `OdeSolution`).
    pub structured_results: bool,
    /// Higher-order operators accept host callbacks.
    pub higher_order: bool,
}

impl Features {
    /// What this build provides, read from the compile-time configuration.
    #[must_use]
    pub const fn current() -> Self {
        Self {
            dynamic_call: true,
            plugin_abi: crate::kernel::signature::PLUGIN_ABI,
            jit: cfg!(feature = "jit"),
            jit_branch: cfg!(feature = "jit"),
            jit_index: cfg!(feature = "jit"),
            parallel_symbolic: true,
            exact_units: true,
            affine_units: true,
            user_base_dimensions: true,
            structured_results: true,
            higher_order: true,
        }
    }
}

/// One rule set as the manifest describes it.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct RuleSetInfo {
    /// The set's name.
    pub name: String,
    /// Its version.
    pub version: String,
    /// How settled it is.
    pub stability: Stability,
    /// Version it appeared in.
    pub since: String,
    /// Capabilities it provides (`calculus`, `ode`, `units`, ...).
    pub capabilities: Vec<String>,
    /// Keys of the operators it registered.
    pub ops: Vec<String>,
}

/// Where an operator's signature came from.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SignatureSource {
    /// The rule set declared it.
    Declared,
    /// It was inferred from the descriptor.
    Inferred,
    /// It came from `signatures.overlay.json`, written by hand.
    Overlay,
    /// There is none; `dynamic_reason` says why.
    None,
}

/// One operator as the manifest describes it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpInfo {
    /// Stable keys, one per overload. `keys[0]` is the primary overload.
    pub keys: Vec<OpKey>,
    /// The operator's name.
    pub name: String,
    /// Other names the dynamic ABI accepts.
    pub aliases: Vec<String>,
    /// Domain category.
    pub category: String,
    /// Accepted argument count.
    pub arity: ArityInfo,
    /// Typed overloads, in declaration order.
    pub overloads: Vec<Signature>,
    /// Where `overloads` came from.
    pub signature_source: SignatureSource,
    /// Algebraic and scheduling properties.
    pub flags: OpFlagsInfo,
    /// Binding structure, when the operator binds a variable.
    pub binder: Option<BinderInfo>,
    /// Extraction cost.
    pub cost: u32,
    /// How the operator produces values.
    pub eval: EvalKind,
    /// Whether it takes a function argument.
    pub higher_order: bool,
    /// Whether it touches anything but its arguments.
    pub purity: Purity,
    /// Whether repeated calls agree.
    pub determinism: Determinism,
    /// How it may be used across threads.
    pub thread_safety: ThreadSafety,
    /// How settled it is.
    pub stability: Stability,
    /// Version it appeared in.
    pub since: String,
    /// Markdown documentation.
    pub doc: String,
    /// Worked examples.
    pub examples: Vec<Example>,
    /// The dimension rule, when it has one.
    pub dim_rule: Option<DimRuleInfo>,
    /// Which backends can execute it.
    pub backend_support: BackendSupport,
    /// Why it has no signature, when it has none.
    pub dynamic_reason: Option<String>,
    /// The rule set that registered it.
    pub rule_set: String,
    /// How the JIT compiles it, when it has an explicit lowering.
    ///
    /// `None` means the backend decides: a scalar evaluator is called
    /// through a shim, everything else runs in the interpreter.
    pub jit: Option<String>,
}

/// One base dimension.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct BaseDimInfo {
    /// Its name (`Length`, `Currency`).
    pub name: String,
    /// Whether it is one of the seven SI base dimensions.
    pub si: bool,
}

/// A plugin, as the manifest describes it.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PluginInfo {
    /// Plugin name.
    pub name: String,
    /// Plugin version.
    pub version: String,
    /// Plugin ABI the plugin was built against.
    pub abi: u32,
    /// What it provides.
    pub capabilities: Vec<String>,
}

/// Everything a host needs to know about this build.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct KernelManifest {
    /// Version of rssn.
    pub rssn_version: String,
    /// Version of the dynamic ABI.
    pub api_version: u32,
    /// BLAKE3 of the manifest with this field zeroed.
    pub manifest_hash: [u8; 32],
    /// What the build can do.
    pub features: Features,
    /// Installed rule sets, sorted by name.
    pub rule_sets: Vec<RuleSetInfo>,
    /// Every operator, sorted by name and then by arity.
    pub ops: Vec<OpInfo>,
    /// Every unit, sorted by name.
    pub units: Vec<UnitValue>,
    /// Every base dimension, sorted by name.
    pub base_dimensions: Vec<BaseDimInfo>,
    /// Registered plugins, sorted by name.
    pub plugins: Vec<PluginInfo>,
}

impl Default for KernelManifest {
    /// An empty manifest: nothing installed. Used as the value of an engine
    /// that has not installed anything yet.
    fn default() -> Self {
        Self {
            rssn_version: Self::rssn_version(),
            api_version: API_VERSION,
            manifest_hash: [0_u8; 32],
            features: Features::current(),
            rule_sets: Vec::new(),
            ops: Vec::new(),
            units: Vec::new(),
            base_dimensions: Vec::new(),
            plugins: Vec::new(),
        }
    }
}

impl KernelManifest {
    /// The version of rssn this build is.
    #[must_use]
    pub fn rssn_version() -> String {
        env!("CARGO_PKG_VERSION").to_owned()
    }

    /// Builds the manifest of an installed configuration.
    ///
    /// `sets` must be the rule sets that were installed into `graph`; they
    /// supply versions and capabilities, while the operators themselves come
    /// from the graph, so a plugin's operators appear here too.
    #[must_use]
    pub fn build(
        graph: &Graph,
        sets: &[RuleSet],
        plugins: &[crate::kernel::plugin::PluginManifest],
    ) -> Self {
        let arena = graph.dim_arena();
        let mut base_dimensions: Vec<BaseDimInfo> = arena
            .bases()
            .iter()
            .map(|name| BaseDimInfo {
                name: name.clone(),
                si: crate::graph::dim::SI_BASES.contains(&name.as_str()),
            })
            .collect();
        base_dimensions.sort();
        base_dimensions.dedup();
        drop(arena);

        let mut ops: Vec<OpInfo> = graph.ops().iter().map(|(_, d)| describe(d, graph)).collect();
        ops.sort_by(|a, b| {
            a.name
                .cmp(&b.name)
                .then_with(|| a.arity.min.cmp(&b.arity.min))
                .then_with(|| a.keys.cmp(&b.keys))
        });

        let mut rule_sets: Vec<RuleSetInfo> = sets.iter().map(|set| rule_set_info(set, &ops)).collect();
        rule_sets.sort_by(|a, b| a.name.cmp(&b.name));
        rule_sets.dedup_by(|a, b| a.name == b.name);

        let mut units: Vec<UnitValue> = crate::dim::units::all()
            .into_iter()
            .map(|unit| UnitValue::from_def(&unit))
            .collect();
        units.sort_by(|a, b| a.name.cmp(&b.name));
        units.dedup_by(|a, b| a.name == b.name);

        let mut plugins: Vec<PluginInfo> = plugins
            .iter()
            .map(|p| PluginInfo {
                name: p.name.clone(),
                version: p.version.clone(),
                abi: p.abi,
                capabilities: p.capabilities.clone(),
            })
            .collect();
        plugins.sort_by(|a, b| a.name.cmp(&b.name));
        plugins.dedup_by(|a, b| a.name == b.name);

        let mut manifest = Self {
            rssn_version: Self::rssn_version(),
            api_version: API_VERSION,
            manifest_hash: [0_u8; 32],
            features: Features::current(),
            rule_sets,
            ops,
            units,
            base_dimensions,
            plugins,
        };
        manifest.seal();
        manifest
    }

    /// Recomputes [`KernelManifest::manifest_hash`].
    ///
    /// Called once at the end of [`KernelManifest::build`]. Exposed because a
    /// host that edits a manifest (to trim it for a smaller target, say) must
    /// be able to seal it again, and because a manifest read from disk must
    /// be checked against its own hash.
    pub fn seal(&mut self) {
        self.manifest_hash = [0_u8; 32];
        self.manifest_hash = codec::digest(self);
    }

    /// Whether the stored hash matches the contents.
    #[must_use]
    pub fn verify(&self) -> bool {
        let mut copy = self.clone();
        copy.seal();
        copy.manifest_hash == self.manifest_hash
    }

    /// The operator registered under `name` or an alias.
    #[must_use]
    pub fn op(
        &self,
        name: &str,
    ) -> Option<&OpInfo> {
        self.ops
            .iter()
            .find(|op| op.name == name || op.aliases.iter().any(|a| a == name))
    }

    /// The rule set named `name`.
    #[must_use]
    pub fn rule_set(
        &self,
        name: &str,
    ) -> Option<&RuleSetInfo> {
        self.rule_sets.iter().find(|set| set.name == name)
    }

    /// The manifest as deterministic CBOR.
    ///
    /// # Errors
    /// See [`codec::encode_cbor`].
    pub fn to_cbor(&self) -> Result<Vec<u8>, codec::CodecError> {
        codec::encode_cbor(self)
    }

    /// The manifest as deterministic binary.
    ///
    /// # Errors
    /// See [`codec::encode_deterministic`].
    pub fn to_bytes(&self) -> Result<Vec<u8>, codec::CodecError> {
        codec::encode_deterministic(self)
    }

    /// Reads a manifest written by [`KernelManifest::to_bytes`].
    ///
    /// # Errors
    /// See [`codec::decode_deterministic`].
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, codec::CodecError> {
        codec::decode_deterministic(bytes)
    }

    /// Reads a manifest written by [`KernelManifest::to_cbor`].
    ///
    /// # Errors
    /// See [`codec::decode_cbor`].
    pub fn from_cbor(bytes: &[u8]) -> Result<Self, codec::CodecError> {
        codec::decode_cbor(bytes)
    }

    /// The manifest as JSON, for the C ABI and for humans.
    ///
    /// # Panics
    /// Panics only if a manifest contains a map with non-string keys, which
    /// the types here make impossible.
    #[must_use]
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_else(|e| format!("{{\"error\":\"{e}\"}}"))
    }

    /// The first 64 bits of [`KernelManifest::manifest_hash`], for ABIs that
    /// carry a `uint64_t` instead of 32 bytes.
    #[must_use]
    pub fn hash_u64(&self) -> u64 {
        let mut bytes = [0_u8; 8];
        bytes.copy_from_slice(self.manifest_hash.get(..8).unwrap_or(&[0_u8; 8]));
        u64::from_le_bytes(bytes)
    }
}

/// Describes one registered operator.
fn describe(
    desc: &OpDescriptor,
    graph: &Graph,
) -> OpInfo {
    let (overloads, source) = signatures_for(desc);
    let keys: Vec<OpKey> = if overloads.is_empty() {
        vec![OpKey::new(&desc.rule_set, &desc.name, None)]
    } else {
        overloads
            .iter()
            .map(|s| OpKey::new(&desc.rule_set, &desc.name, Some(s)))
            .collect()
    };
    let op = graph.ops().lookup(&desc.name).unwrap_or(crate::graph::OpId::NONE);
    let dim_rule = graph.ops().attr::<crate::graph::dim::DimRule>(op).map(|r| DimRuleInfo::from(*r));
    let dynamic_reason = if overloads.is_empty() {
        Some(
            desc.dynamic_reason
                .unwrap_or_else(|| inferred_reason(desc))
                .as_str()
                .to_owned(),
        )
    } else {
        None
    };
    OpInfo {
        keys,
        name: desc.name.to_string(),
        aliases: desc.aliases.iter().map(std::string::ToString::to_string).collect(),
        category: desc.category.to_string(),
        arity: ArityInfo::from_arity(desc.arity),
        overloads,
        signature_source: source,
        flags: desc.flags.info(),
        binder: desc.binder.map(crate::graph::op::Binder::info),
        cost: desc.cost,
        eval: if desc.eval_dyn.is_some() {
            EvalKind::Dynamic
        } else if desc.eval.is_some() {
            EvalKind::Scalar
        } else {
            EvalKind::None
        },
        higher_order: desc.higher_order,
        purity: desc.purity,
        determinism: desc.determinism,
        thread_safety: desc.thread_safety,
        jit: desc.jit.map(lowering_name),
        stability: desc.stability,
        since: desc.since.to_string(),
        doc: desc.doc.to_string(),
        examples: desc.examples.clone(),
        dim_rule,
        backend_support: desc.backend,
        dynamic_reason,
        rule_set: desc.rule_set.to_string(),
    }
}

/// The signatures of an operator and where they came from.
fn signatures_for(desc: &OpDescriptor) -> (Vec<Signature>, SignatureSource) {
    if !desc.overloads.is_empty() {
        return (desc.overloads.clone(), SignatureSource::Declared);
    }
    if let Some((signature, source)) = crate::kernel::infer::infer_signature_with_source(desc) {
        return (vec![signature], source);
    }
    (Vec::new(), SignatureSource::None)
}

/// Why an operator that declared no reason has no signature.
///
/// `NotYetSpecified` is the last resort, not the default: a rule-only
/// operator has a *real* explanation — there is nothing to call, and the
/// result type is whatever the rule that fires produces — and saying so is
/// more useful to a host than "we did not write it down".
fn inferred_reason(desc: &OpDescriptor) -> DynamicReason {
    if matches!(desc.arity, crate::graph::Arity::Variadic) {
        return DynamicReason::Variadic;
    }
    if desc.eval.is_none() && desc.eval_dyn.is_none() {
        return DynamicReason::Request;
    }
    DynamicReason::NotYetSpecified
}

/// The name a host sees for a JIT lowering.
fn lowering_name(lowering: crate::backend::Lowering) -> String {
    match lowering {
        | crate::backend::Lowering::Intrinsic(intrinsic) => {
            format!("intrinsic:{intrinsic:?}").to_lowercase()
        },
        | crate::backend::Lowering::Extern1(_) => String::from("extern1"),
        | crate::backend::Lowering::Extern2(_) => String::from("extern2"),
        | crate::backend::Lowering::Gather => String::from("gather"),
    }
}

/// Describes one rule set, filling in the operators it registered.
fn rule_set_info(
    set: &RuleSet,
    ops: &[OpInfo],
) -> RuleSetInfo {
    let meta = set.meta();
    RuleSetInfo {
        name: set.name().to_owned(),
        version: meta.version.clone(),
        stability: meta.stability,
        since: meta.since.clone(),
        capabilities: meta.capabilities.clone(),
        ops: ops
            .iter()
            .filter(|op| op.rule_set == set.name())
            .flat_map(|op| op.keys.iter().map(OpKey::qualified))
            .collect(),
    }
}

/// A hand-maintained overlay of signatures, applied after inference.
#[derive(Clone, Debug, Default, Deserialize)]
struct Overlay {
    /// Signatures by operator name.
    #[serde(default)]
    ops: BTreeMap<String, Vec<Signature>>,
}

/// The overlay shipped with the crate.
const OVERLAY_JSON: &str = include_str!("../../signatures.overlay.json");

/// The parsed overlay, or an empty one when the file is unreadable.
///
/// A malformed overlay must not take the process down: it is documentation
/// and inference correction, not data the engine needs. The audit binary
/// reports it.
pub(crate) fn overlay() -> &'static BTreeMap<String, Vec<Signature>> {
    static OVERLAY: std::sync::OnceLock<BTreeMap<String, Vec<Signature>>> =
        std::sync::OnceLock::new();
    OVERLAY.get_or_init(|| {
        serde_json::from_str::<Overlay>(OVERLAY_JSON)
            .map(|o| o.ops)
            .unwrap_or_default()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::rules;

    fn standard_manifest() -> KernelManifest {
        let sets = rules::standard();
        let mut graph = Graph::new();
        Engine::install(&mut graph, &sets).unwrap_or_else(|e| panic!("{e}"));
        KernelManifest::build(&graph, &sets, &[])
    }

    #[test]
    fn the_manifest_describes_every_operator() {
        let manifest = standard_manifest();
        assert!(manifest.ops.len() > 1000, "only {} ops", manifest.ops.len());
        assert_eq!(manifest.api_version, API_VERSION);
        assert!(manifest.verify(), "the hash must match the contents");
        for op in &manifest.ops {
            assert!(
                !op.overloads.is_empty() || op.dynamic_reason.is_some(),
                "`{}` is neither typed nor explained",
                op.name
            );
            assert_eq!(op.keys.len(), op.overloads.len().max(1));
        }
    }

    #[test]
    fn the_manifest_is_deterministic() {
        let a = standard_manifest();
        let b = standard_manifest();
        assert_eq!(a.manifest_hash, b.manifest_hash);
        assert_eq!(
            a.to_bytes().unwrap_or_default(),
            b.to_bytes().unwrap_or_default()
        );
        // Rule-set order must not matter.
        let mut sets = rules::standard();
        sets.reverse();
        let mut graph = Graph::new();
        Engine::install(&mut graph, &sets).unwrap_or_else(|e| panic!("{e}"));
        let c = KernelManifest::build(&graph, &sets, &[]);
        assert_eq!(a.manifest_hash, c.manifest_hash);
    }

    #[test]
    fn the_manifest_round_trips_through_both_encodings() {
        let manifest = standard_manifest();
        let bytes = manifest.to_bytes().unwrap_or_else(|e| panic!("{e}"));
        let back = KernelManifest::from_bytes(&bytes).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, manifest);
        let cbor = manifest.to_cbor().unwrap_or_else(|e| panic!("{e}"));
        let back = KernelManifest::from_cbor(&cbor).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, manifest);
    }

    #[test]
    fn rule_sets_list_their_own_operators() {
        let manifest = standard_manifest();
        let calculus = manifest.rule_set("calculus").unwrap_or_else(|| panic!("no calculus"));
        assert!(
            calculus
                .ops
                .iter()
                .any(|k| k.split('.').next_back() == Some("diff"))
        );
        assert!(manifest.op("diff").is_some());
        assert!(manifest.op("sin").is_some());
    }

    #[test]
    fn units_and_base_dimensions_are_listed() {
        let manifest = standard_manifest();
        assert!(manifest.units.iter().any(|u| u.name == "m"));
        assert!(manifest.units.iter().any(|u| u.name == "degC"));
        assert!(manifest.base_dimensions.iter().any(|b| b.name == "Length" && b.si));
    }
}
