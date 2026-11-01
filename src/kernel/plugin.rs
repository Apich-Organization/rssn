//! Plugins and host callbacks.
//!
//! A [`Plugin`] is a bundle of operators, rules, dimension rules and units
//! that a third party adds to a session at run time. Registering one is not
//! a special mode: the plugin's operators appear in the
//! [manifest](crate::kernel::manifest) exactly like the built-in ones, are
//! callable through the same dynamic ABI, and take part in the same search.
//!
//! A [`HostCallback`] is the smaller thing: one function from values to a
//! value, handed to a higher-order operator for the duration of a call
//! ([`Value::Function`](crate::kernel::value::Value::Function)). That is how
//! a host passes its own closure into `map`, `fold` or a root finder.
//!
//! Neither needs to be `unsafe`, and neither can corrupt a session: a plugin
//! registers through [`Registry`], which goes through the same checked
//! registration path a rule set uses, and a callback only ever sees values.

use std::fmt;
use std::sync::Arc;

use serde::Deserialize;
use serde::Serialize;

use crate::graph::Graph;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Program;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::dim::BaseId;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::PLUGIN_ABI;
use crate::kernel::signature::Signature;
use crate::kernel::value::Value;

/// Identifies a registered plugin within a session.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PluginId(pub u32);

impl fmt::Display for PluginId {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "plugin#{}", self.0)
    }
}

/// What a plugin declares about itself.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PluginManifest {
    /// Unique name, by convention `vendor/name`.
    pub name: String,
    /// Version of the plugin.
    pub version: String,
    /// Plugin ABI the plugin was built against.
    pub abi: u32,
    /// How the plugin is loaded.
    pub kind: PluginKind,
    /// What it provides (`ops`, `rules`, `units`, ...).
    pub capabilities: Vec<String>,
    /// One-line description.
    pub description: String,
}

/// How a plugin reaches the host.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PluginKind {
    /// Statically linked Rust, registered through [`Session::register`].
    ///
    /// [`Session::register`]: crate::api::Session::register
    #[default]
    InProcess,
    /// Loaded from a shared library at run time
    /// ([`cdylib::load`](crate::kernel::cdylib::load)); needs the
    /// `plugin-cdylib` feature.
    Cdylib,
}

impl fmt::Display for PluginKind {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(match self {
            | Self::InProcess => "in-process",
            | Self::Cdylib => "cdylib",
        })
    }
}

/// Why a plugin could not be registered or called.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PluginError {
    /// The plugin was built against a different ABI.
    AbiMismatch {
        /// What this build provides.
        expected: u32,
        /// What the plugin wants.
        got: u32,
    },
    /// A plugin with that name is already registered.
    Duplicate {
        /// The name.
        name: String,
    },
    /// Registration failed.
    Registration(String),
    /// The plugin refused the call.
    Rejected(String),
    /// The plugin panicked. rssn catches it; the session stays usable.
    Panic(String),
}

impl fmt::Display for PluginError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::AbiMismatch { expected, got } => {
                write!(f, "plugin ABI mismatch: host provides {expected}, plugin wants {got}")
            },
            | Self::Duplicate { name } => write!(f, "plugin `{name}` is already registered"),
            | Self::Registration(message) => write!(f, "plugin registration failed: {message}"),
            | Self::Rejected(message) => write!(f, "plugin refused the call: {message}"),
            | Self::Panic(message) => write!(f, "plugin panicked: {message}"),
        }
    }
}

impl std::error::Error for PluginError {}

impl From<crate::kernel::call::OpError> for PluginError {
    /// A failed registration call is a registration failure; the operator's
    /// own message is more useful than a generic one, so it is kept.
    fn from(error: crate::kernel::call::OpError) -> Self {
        Self::Registration(error.to_string())
    }
}

/// A bundle of extensions.
pub trait Plugin: Send + Sync {
    /// What the plugin declares about itself.
    fn manifest(&self) -> PluginManifest;

    /// Registers the plugin's operators, rules and units.
    ///
    /// # Errors
    /// Returns [`PluginError`] when a registration is refused.
    fn register(
        &self,
        reg: &mut Registry<'_>,
    ) -> Result<(), PluginError>;

    /// Handles a call to one of the plugin's operators.
    ///
    /// Only reached for operators the plugin registered with an
    /// [`EvalDyn`](crate::graph::op::EvalDyn); operators that are plain
    /// rules or scalar evaluators never get here.
    ///
    /// # Errors
    /// Returns [`OpError`] when the call cannot be answered.
    fn call(
        &self,
        op: &crate::kernel::signature::OpKey,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let _ = (args, ctx);
        Err(OpError::Domain(format!(
            "plugin does not implement `{}`",
            op.qualified()
        )))
    }
}

/// A function the host lends to a call.
pub trait HostCallback: Send + Sync {
    /// Applies the function.
    ///
    /// # Errors
    /// Whatever the host wants to report, as an [`OpError`].
    fn call(
        &self,
        args: &[Value],
    ) -> Result<Value, OpError>;

    /// The function's type.
    fn signature(&self) -> Signature;
}

/// A [`HostCallback`] built from a closure.
pub struct NativeCallback<F> {
    f: F,
    signature: Signature,
}

impl<F> NativeCallback<F>
where F: Fn(&[Value]) -> Result<Value, OpError> + Send + Sync
{
    /// Wraps `f` with the given signature.
    #[must_use]
    pub const fn new(
        f: F,
        signature: Signature,
    ) -> Self {
        Self { f, signature }
    }
}

impl<F> HostCallback for NativeCallback<F>
where F: Fn(&[Value]) -> Result<Value, OpError> + Send + Sync
{
    fn call(
        &self,
        args: &[Value],
    ) -> Result<Value, OpError> {
        (self.f)(args)
    }

    fn signature(&self) -> Signature {
        self.signature.clone()
    }
}

impl fmt::Debug for dyn HostCallback {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str("HostCallback")
    }
}

/// What a plugin is allowed to do while registering.
///
/// The registry is a checked view of the session's operator table and rule
/// program. It cannot remove anything, cannot touch terms, and refuses
/// registrations that would change the meaning of an existing operator.
pub struct Registry<'g> {
    installer: Installer<'g>,
    plugin: String,
}

impl<'g> Registry<'g> {
    /// A registry that appends to `program` while registering into `graph`.
    #[must_use]
    pub fn new(
        graph: &'g mut Graph,
        program: &'g mut Program,
        plugin: &str,
    ) -> Self {
        Self {
            installer: Installer::new(graph, program),
            plugin: plugin.to_owned(),
        }
    }

    /// The name of the plugin being registered.
    #[must_use]
    pub fn plugin(&self) -> &str {
        &self.plugin
    }

    /// Registers an operator.
    ///
    /// # Errors
    /// Fails when the name is taken by an operator with a different shape,
    /// or when the descriptor is malformed (no signature and no
    /// `dynamic_reason`).
    pub fn register_op(
        &mut self,
        desc: OpDescriptor,
    ) -> Result<OpId, OpError> {
        if desc.overloads.is_empty() && desc.dynamic_reason.is_none() {
            return Err(OpError::Internal(format!(
                "plugin `{}` registered `{}` without a signature or a dynamic reason",
                self.plugin, desc.name
            )));
        }
        self.installer
            .op(desc)
            .map_err(|e| OpError::Plugin(PluginError::Registration(e.to_string())))
    }

    /// Adds declarative rewrites, one per entry.
    ///
    /// # Errors
    /// Fails on the first entry that does not parse or is not a valid
    /// identity schema.
    pub fn register_rewrites(
        &mut self,
        tier: Tier,
        texts: &[&str],
    ) -> Result<(), OpError> {
        self.installer
            .rewrites(tier, texts)
            .map_err(|e| OpError::Plugin(PluginError::Registration(e.to_string())))
    }

    /// Adds a definition (`name(a, b) := body`).
    ///
    /// # Errors
    /// See [`Installer::define`].
    pub fn register_definitions(
        &mut self,
        texts: &[&str],
    ) -> Result<(), OpError> {
        self.installer
            .define(texts)
            .map_err(|e| OpError::Plugin(PluginError::Registration(e.to_string())))
    }

    /// Attaches a dimension rule to an operator.
    ///
    /// # Errors
    /// Fails when `op` does not belong to this graph.
    pub fn register_dim_rule(
        &mut self,
        op: OpId,
        rule: DimRule,
    ) -> Result<(), OpError> {
        if self.installer.graph().ops().try_get(op).is_none() {
            return Err(OpError::Internal(format!(
                "plugin `{}` attached a dimension rule to a foreign operator",
                self.plugin
            )));
        }
        self.installer.dim(op, rule);
        Ok(())
    }

    /// Registers a unit.
    ///
    /// # Errors
    /// Fails when the name or an alias is already taken by a different
    /// unit.
    pub fn register_unit(
        &mut self,
        unit: crate::dim::UnitDef,
    ) -> Result<(), OpError> {
        crate::dim::units::register(unit).map_err(|e| OpError::Domain(e.to_string()))
    }

    /// Registers a base dimension.
    ///
    /// # Errors
    /// Fails when the name is empty.
    pub fn register_base_dimension(
        &mut self,
        name: &str,
    ) -> Result<BaseId, OpError> {
        if name.trim().is_empty() {
            return Err(OpError::Domain("a base dimension needs a name".to_owned()));
        }
        self.installer.graph().define_base_dimension(name);
        let arena = self.installer.graph().dim_arena();
        arena
            .base_id(name)
            .ok_or_else(|| OpError::Internal("base dimension vanished".to_owned()))
    }

    /// Adds a window pass.
    pub fn register_pass(
        &mut self,
        pass: impl crate::graph::window::WindowPass + 'static,
    ) {
        self.installer.pass(pass);
    }

    /// The graph being extended.
    pub const fn graph(&mut self) -> &mut Graph {
        self.installer.graph()
    }
}

/// A plugin that only registers data, built from a closure.
///
/// ```no_run
/// use rssn::kernel::plugin::SimplePlugin;
///
/// let plugin = SimplePlugin::new("acme/demo", "0.1.0", |reg| {
///     let op = reg.register_op(
///         rssn::graph::OpDescriptor::new("acme_bump", rssn::graph::Arity::Fixed(1))
///             .signature(rssn::kernel::Signature::scalars(1, rssn::kernel::TypeExpr::Real)),
///     )?;
///     reg.register_rewrites(rssn::graph::Tier::Normalize, &["acme/bump-zero: acme_bump(0) => 0"])?;
///     Ok(())
/// })
/// .with_description("An example plugin.");
/// ```
pub struct SimplePlugin {
    manifest: PluginManifest,
    install: Arc<InstallFn>,
}

/// The install closure a [`SimplePlugin`] wraps.
pub type InstallFn = dyn Fn(&mut Registry<'_>) -> Result<(), PluginError> + Send + Sync;

impl SimplePlugin {
    /// A plugin named `name`, version `version`, installing through `install`.
    #[must_use]
    pub fn new(
        name: &str,
        version: &str,
        install: impl Fn(&mut Registry<'_>) -> Result<(), PluginError> + Send + Sync + 'static,
    ) -> Self {
        Self {
            manifest: PluginManifest {
                kind: PluginKind::InProcess,
                name: name.to_owned(),
                version: version.to_owned(),
                abi: PLUGIN_ABI,
                capabilities: vec!["ops".to_owned(), "rules".to_owned()],
                description: String::new(),
            },
            install: Arc::new(install),
        }
    }

    /// Sets the description.
    #[must_use]
    pub fn with_description(
        mut self,
        description: &str,
    ) -> Self {
        description.clone_into(&mut self.manifest.description);
        self
    }

    /// Declares an extra capability.
    #[must_use]
    pub fn with_capability(
        mut self,
        capability: &str,
    ) -> Self {
        self.manifest.capabilities.push(capability.to_owned());
        self
    }
}

impl Plugin for SimplePlugin {
    fn manifest(&self) -> PluginManifest {
        self.manifest.clone()
    }

    fn register(
        &self,
        reg: &mut Registry<'_>,
    ) -> Result<(), PluginError> {
        (self.install)(reg)
    }
}

impl From<RuleError> for PluginError {
    fn from(value: RuleError) -> Self {
        Self::Registration(value.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Arity;
    use crate::kernel::Signature;
    use crate::kernel::TypeExpr;

    #[test]
    fn a_plugin_registers_operators_rules_and_units() {
        let plugin = SimplePlugin::new("acme/demo", "0.1.0", |reg| {
            reg.register_op(
                OpDescriptor::new("acme_double", Arity::Fixed(1))
                    .signature(Signature::scalars(1, TypeExpr::Real)),
            )?;
            reg.register_rewrites(
                Tier::Normalize,
                &["acme/double-zero: acme_double(0) => 0"],
            )?;
            Ok(())
        });
        let mut graph = Graph::new();
        let mut program = Program::default();
        plugin
            .register(&mut Registry::new(&mut graph, &mut program, "acme/demo"))
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(graph.ops().lookup("acme_double").is_some());
        assert_eq!(program.rules.len(), 1);
        assert_eq!(plugin.manifest().abi, PLUGIN_ABI);
    }

    #[test]
    fn an_untyped_operator_is_refused() {
        let plugin = SimplePlugin::new("acme/bad", "0.1.0", |reg| {
            reg.register_op(OpDescriptor::new("acme_bare", Arity::Fixed(1)))?;
            Ok(())
        });
        let mut graph = Graph::new();
        let mut program = Program::default();
        assert!(
            plugin
                .register(&mut Registry::new(&mut graph, &mut program, "acme/bad"))
                .is_err()
        );
    }

    #[test]
    fn host_callbacks_are_closures() {
        let cb = NativeCallback::new(
            |args: &[Value]| {
                let first = args
                    .first()
                    .cloned()
                    .ok_or_else(|| OpError::Arity {
                        expected: crate::kernel::ArityInfo::fixed(1),
                        got: 0,
                    })?;
                Ok(first)
            },
            Signature::scalars(1, TypeExpr::Real),
        );
        assert_eq!(cb.call(&[Value::int(3_i32)]), Ok(Value::int(3_i32)));
        assert!(cb.call(&[]).is_err());
        assert_eq!(cb.signature().arity(), 1);
    }
}
