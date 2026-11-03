//! The bridge between terms and [`Value`]s.
//!
//! This module implements the dynamic ABI on [`Session`]: the three
//! conversions of the specification
//! ([`Session::value_of`], [`Session::term_of`], [`Session::call_value`]),
//! the owned [`TermRef`] handle, `fork`/`merge`, and plugin registration.
//!
//! Everything here goes through the ordinary [`Session::compute`] path, so a
//! dynamic call has exactly the semantics of a typed one: same rules, same
//! budget, same dimension checks. There is no second engine.

use std::fmt;
use std::sync::Arc;

use crate::api::Config;
use crate::api::Session;
use crate::api::Term;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Payload;
use crate::graph::op::core;
use crate::kernel::call::CallCtx;
use crate::kernel::call::Kernel;
use crate::kernel::call::OpError;
use crate::kernel::call::Target;
use crate::kernel::manifest::KernelManifest;
use crate::kernel::plugin::Plugin;
use crate::kernel::plugin::PluginError;
use crate::kernel::plugin::PluginId;
use crate::kernel::plugin::Registry;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::OpKey;
use crate::kernel::signature::Signature;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::QuantityValue;
use crate::kernel::value::SessionId;
use crate::kernel::value::TensorValue;
use crate::kernel::value::TermRef;
use crate::kernel::value::Value;

/// Why two sessions could not be merged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MergeError {
    /// The two sessions declared different dimensions for one symbol.
    Dimension {
        /// The symbol.
        symbol: String,
        /// What the receiving session declared.
        left: String,
        /// What the incoming session declared.
        right: String,
    },
    /// A symbol's dimension could not be translated between arenas.
    Translation(String),
}

impl fmt::Display for MergeError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Dimension { symbol, left, right } => {
                write!(
                    f,
                    "symbol `{symbol}` has dimension {left} here but {right} in the merged session"
                )
            },
            | Self::Translation(message) => write!(f, "cannot translate a dimension: {message}"),
        }
    }
}

impl std::error::Error for MergeError {}

impl Session {
    /// The manifest of this session's configuration.
    ///
    /// Returned as a shared handle rather than a reference because a session
    /// keeps its templates behind a `RefCell`; the manifest itself is
    /// immutable once built.
    #[must_use]
    pub fn manifest(&self) -> Arc<KernelManifest> {
        let inner = self.inner.borrow();
        inner.templates.first().map_or_else(
            || Arc::new(KernelManifest::default()),
            |t| t.engine.manifest_arc(),
        )
    }

    /// The manifest as JSON.
    #[must_use]
    pub fn manifest_json(&self) -> String {
        self.manifest().to_json()
    }

    /// The id other sessions use to reject this session's handles.
    #[must_use]
    pub const fn id(&self) -> SessionId {
        self.id
    }

    /// An owned handle for `term`.
    #[must_use]
    pub const fn handle(
        &self,
        term: Term<'_>,
    ) -> TermRef {
        TermRef {
            session: self.id,
            node: term.id(),
        }
    }

    /// The term a handle names.
    ///
    /// # Errors
    /// Fails when the handle belongs to another session or names a node this
    /// session does not have.
    pub fn term_ref(
        &self,
        handle: TermRef,
    ) -> Result<Term<'_>, OpError> {
        if handle.session != self.id {
            return Err(OpError::ForeignTerm {
                expected: handle.session,
                got: self.id,
            });
        }
        self.term_by_id(handle.node)
            .ok_or_else(|| OpError::Internal(format!("no term {} in this session", handle.node)))
    }

    /// The value of `term`.
    ///
    /// Numbers become [`Value::Int`], [`Value::Real`] or [`Value::Float`];
    /// symbols become [`Value::Symbol`]; a quantity becomes
    /// [`Value::Quantity`]; a rectangular list of numbers becomes
    /// [`Value::Tensor`]. Anything else becomes [`Value::Term`], a handle to
    /// the term itself — the honest answer, because a symbolic expression is
    /// a term and nothing else.
    ///
    /// # Errors
    /// Fails only when a value cannot be represented, which today means a
    /// unit whose definition is inconsistent with the graph's arena.
    pub fn value_of(
        &self,
        term: Term<'_>,
    ) -> Result<Value, OpError> {
        let inner = self.inner.borrow();
        let graph = &inner.graph;
        Ok(value_of_node(graph, NodeId::from_raw(term.id()), self.id))
    }

    /// The term a value denotes.
    ///
    /// A unit is not a value: [`Value::Unit`] is a definition and is refused,
    /// because the term graph has no node for a dimension plus a scale.
    /// [`Value::Quantity`] carries the unit and converts normally, and a unit
    /// passed to an operator that declares a unit parameter (the second
    /// argument of `quantity`) is still accepted by [`Session::call_value`].
    ///
    /// # Errors
    /// Fails for values that have no term form: a bare unit, a host
    /// callback, an opaque payload, a sparse tensor, a record or a variant.
    /// Those are boundary values; they reach an operator through its
    /// [`EvalDyn`], not through the term graph.
    ///
    /// [`EvalDyn`]: crate::graph::op::EvalDyn
    pub fn term_of(
        &self,
        value: &Value,
    ) -> Result<Term<'_>, OpError> {
        let mut inner = self.inner.borrow_mut();
        let node = term_of_value(&mut inner.graph, value, self.id)?;
        Ok(self.term(node))
    }

    /// Calls the operator `name` on `args`.
    ///
    /// Resolution is by name (or alias) and argument count, so an overloaded
    /// operator picks the overload that accepts `args.len()` arguments. An
    /// operator with a dynamic evaluator is evaluated by it; every other
    /// operator is built into a term and reduced by the engine exactly as
    /// [`Session::compute`] would.
    ///
    /// # Errors
    /// * [`OpError::UnknownOp`] when no operator has that name;
    /// * [`OpError::Arity`] when it exists but takes a different number of
    ///   arguments;
    /// * [`OpError::Type`] when an argument has no term form;
    /// * [`OpError::NotNumeric`], [`OpError::NotReduced`],
    ///   [`OpError::Dimension`] from the computation itself.
    pub fn call_value(
        &self,
        name: &str,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        self.call_with(name, args, &[], ctx)
    }

    /// Solves `equation = 0` for `variable` and returns a structured answer.
    ///
    /// The `solve` operator itself answers with a `list(...)`, which is the
    /// right thing for a rewrite but loses everything a caller needs to know
    /// besides the values. This lifts that answer into
    /// [`Solutions`], which carries the unknown and whether the set is
    /// complete, and converts to a [`Value::Record`] for the dynamic ABI.
    ///
    /// Multiplicities are left unreported: the solver returns distinct
    /// roots, and claiming a multiplicity it did not compute would be worse
    /// than saying nothing.
    ///
    /// # Errors
    /// Fails when `solve` does not apply, or when the equation does not
    /// reduce to a closed form.
    pub fn solve_structured(
        &self,
        equation: Term<'_>,
        variable: Term<'_>,
    ) -> Result<crate::kernel::result::Solutions, OpError> {
        let name = variable.to_string();
        let value = self.call_value(
            "solve",
            &[
                Value::Term(self.handle(equation)),
                Value::Term(self.handle(variable)),
            ],
            &CallCtx::new(),
        )?;
        // A `list(...)` answer is the solver's exhaustive form; a bare value
        // is a single root found on the way, which is not a claim that the
        // others do not exist. Reporting the difference is the point of the
        // field.
        let (roots, complete) = match value {
            | Value::Tensor(tensor) => (tensor.data, true),
            | Value::Sparse(sparse) => (sparse.data, true),
            | other => (vec![other], false),
        };
        Ok(crate::kernel::result::Solutions {
            variable: Some(name),
            roots,
            multiplicities: Vec::new(),
            conditions: Vec::new(),
            complete,
        })
    }

    /// Calls a higher-order operator with host callbacks.
    ///
    /// The callbacks are reachable from the operator's
    /// [`EvalDyn`](crate::graph::op::EvalDyn) through
    /// [`CallCtx::callbacks`] and the [`CallbackRef`]s in `fns`. Rewrites
    /// never call back into the host: a rewrite must be a pure identity
    /// transformation, so an operator that needs the host must be a dynamic
    /// evaluator.
    ///
    /// # Errors
    /// See [`Session::call_value`].
    pub fn call_ho(
        &self,
        name: &str,
        args: &[Value],
        fns: &[CallbackRef],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        self.call_with(name, args, fns, ctx)
    }

    /// Every typed overload of `name`, with the key to call it by.
    ///
    /// An operator with no declared or inferred signature has no overloads and
    /// returns an empty list; it is reachable by [`Session::call_value`] but
    /// not by [`Session::call_key`].
    #[must_use]
    pub fn overloads(
        &self,
        name: &str,
    ) -> Vec<(OpKey, Signature)> {
        let inner = self.inner.borrow();
        let graph = &inner.graph;
        graph
            .ops()
            .lookup_all(name)
            .iter()
            .flat_map(|&id| {
                let descriptor = graph.ops().get(id);
                let (signatures, _) = crate::kernel::infer::signatures(descriptor);
                signatures
                    .iter()
                    .map(|signature| {
                        (
                            OpKey::new(&descriptor.rule_set, &descriptor.name, Some(signature)),
                            signature.clone(),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// The overload of `name` that `args` belong to.
    ///
    /// Resolution is by arity first and then by type, most specific first; a
    /// tie is [`OpError::Ambiguous`] rather than a silent choice. This is the
    /// same decision [`Session::call_value`] makes, exposed so that a
    /// generated typed front end can choose statically and cache on the key.
    ///
    /// # Errors
    /// [`OpError::UnknownOp`] when nothing is registered under `name`,
    /// [`OpError::Arity`] when no overload takes that many arguments, and
    /// [`OpError::Ambiguous`] when several tie.
    pub fn resolve_overload(
        &self,
        name: &str,
        args: &[Value],
    ) -> Result<OpKey, OpError> {
        let inner = self.inner.borrow();
        let graph = &inner.graph;
        let op = resolve_op(graph, name, args)?;
        let descriptor = graph.ops().get(op);
        let (signatures, _) = crate::kernel::infer::signatures(descriptor);
        let chosen = match crate::kernel::overload::best(&signatures, args) {
            | crate::kernel::overload::Best::Unique(index) => signatures.get(index),
            // An operator with no signature has exactly one key, the
            // all-zero one, and `call_key` accepts it.
            | _ if signatures.is_empty() => None,
            // Several overloads fit and none is more specific. Returning a
            // key here would be a key `call_key` refuses, so say so instead.
            // `count` is the number that *tied*, not the number of overloads:
            // a host quotes it in its own diagnostic, and "3 candidates" for a
            // two-way tie sends the reader looking for a third.
            | crate::kernel::overload::Best::Ambiguous(indices) => {
                return Err(OpError::Ambiguous {
                    name: name.to_owned(),
                    count: indices.len(),
                });
            },
            // Nothing accepts the arguments. There is no more specific error
            // to give without re-deriving one signature's refusal, and the
            // callers that reach here are hosts resolving from values, which
            // the typed front end has already checked.
            | crate::kernel::overload::Best::None => {
                return Err(OpError::Ambiguous {
                    name: name.to_owned(),
                    count: signatures.len(),
                });
            },
        };
        Ok(OpKey::new(&descriptor.rule_set, &descriptor.name, chosen))
    }

    /// Calls exactly the overload `key` names.
    ///
    /// [`Session::call_value`] resolves an overload from the arguments; this
    /// does not resolve anything. A host that generated a typed declaration
    /// from `OpInfo.keys` should call through the key it generated, so that
    /// the static choice and the executed operator cannot disagree.
    ///
    /// An all-zero `signature_hash` means "the registration with no declared
    /// overload", which is how an untyped operator is named.
    ///
    /// # Errors
    /// [`OpError::UnknownOp`] when no registration matches the key, plus
    /// everything [`Session::call_value`] can return.
    pub fn call_key(
        &self,
        key: &OpKey,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let ctx = clone_ctx(ctx);
        let op = self.op_for_key(key, args)?;
        self.call_op(op, &key.name, args, &ctx)
    }

    /// The registration `key` names, for these arguments.
    ///
    /// A key names one *overload*, not merely one registration: checking only
    /// that some signature of the registration has the pinned hash made
    /// "the overload the host chose is the overload that runs" true per
    /// registration rather than per overload. The arguments are validated
    /// against the pinned signature so the two cannot drift.
    fn op_for_key(
        &self,
        key: &OpKey,
        args: &[Value],
    ) -> Result<OpId, OpError> {
        let inner = self.inner.borrow();
        let graph = &inner.graph;
        for &id in graph.ops().lookup_all(&key.name) {
            let descriptor = graph.ops().get(id);
            if !key.rule_set.is_empty() && *descriptor.rule_set != *key.rule_set {
                continue;
            }
            let (signatures, _) = crate::kernel::infer::signatures(descriptor);
            if key.signature_hash == [0_u8; 16] {
                if signatures.is_empty() {
                    return Ok(id);
                }
            } else if let Some(pinned) = signatures
                .iter()
                .find(|signature| signature.hash() == key.signature_hash)
            {
                // The pinned overload must be one these arguments can reach,
                // or the executed overload is not the chosen one.
                if crate::kernel::overload::score(pinned, args).is_some() {
                    return Ok(id);
                }
            }
        }
        Err(OpError::UnknownOp(key.clone()))
    }

    fn call_with(
        &self,
        name: &str,
        args: &[Value],
        fns: &[CallbackRef],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let ctx = clone_ctx(ctx);
        check_callbacks(&ctx, fns)?;
        let op = {
            let inner = self.inner.borrow();
            resolve_op(&inner.graph, name, args)?
        };
        self.call_op(op, name, args, &ctx)
    }

    /// Calls the registration `op`, which has already been resolved.
    fn call_op(
        &self,
        op: OpId,
        name: &str,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let sets = {
            let inner = self.inner.borrow();
            inner
                .templates
                .first()
                .map(|t| t.rule_sets.clone())
                .unwrap_or_default()
        };
        let descriptor = {
            let inner = self.inner.borrow();
            inner.graph.ops().get(op).clone()
        };
        ctx.trace(
            "resolve",
            format!("{name}/{} -> {}", args.len(), op.index()),
        );
        if let Some(evaluator) = descriptor.eval_dyn.clone() {
            ctx.trace("dispatch", "native evaluator");
            return evaluator.eval(args, ctx);
        }
        let unit_positions = unit_parameters(&descriptor, args.len());
        let term = {
            let mut inner = self.inner.borrow_mut();
            let mut nodes = args
                .iter()
                .enumerate()
                .map(|(index, value)| {
                    term_of_argument(
                        &mut inner.graph,
                        value,
                        self.id,
                        unit_positions.get(index).copied().unwrap_or(false),
                    )
                    .map_err(|error| {
                        match error {
                            | OpError::Type { expected, .. } => {
                                OpError::Type {
                                    expected,
                                    got: Box::new(value.type_expr()),
                                    arg: index,
                                }
                            },
                            | other => other,
                        }
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            // `CallCtx::substitute` used to be recorded and never read: a host
            // that replaced a symbol before a call got the unsubstituted
            // answer with no diagnostic. Apply it here, where the arguments
            // become terms.
            for (symbol, value) in &ctx.substitutions {
                let from = inner.graph.sym(symbol);
                let to = term_of_value(&mut inner.graph, value, self.id)?;
                for node in &mut nodes {
                    *node = inner.graph.substitute(*node, from, to);
                }
            }
            let node = inner.graph.try_node(op, &nodes).ok_or_else(|| {
                OpError::Arity {
                    expected: ArityInfo::from_arity(descriptor.arity),
                    got: args.len(),
                }
            })?;
            self.term(node)
        };
        let config = config_for(&sets, ctx);
        let answer = self.compute(term, &config).map_err(op_error)?;
        ctx.trace("compute", format!("reduced={}", answer.reduced));
        self.value_of(answer.term)
    }

    /// Loads a plugin from a shared library and registers it.
    ///
    /// This is [`cdylib::load`](crate::kernel::cdylib::load) followed by
    /// [`Session::register`]; the library stays mapped for as long as the
    /// session holds the plugin.
    ///
    /// # Errors
    /// See [`cdylib::load`](crate::kernel::cdylib::load) and
    /// [`Session::register`].
    #[cfg(feature = "plugin-cdylib")]
    pub fn load_plugin(
        &mut self,
        path: &std::path::Path,
    ) -> Result<PluginId, PluginError> {
        let plugin = crate::kernel::cdylib::load(path)?;
        self.register(Box::new(plugin))
    }

    /// Registers a plugin.
    ///
    /// The plugin becomes part of the session's configuration: every
    /// template — including one built later for a different set of rule
    /// sets — installs it, and its operators appear in the
    /// [manifest](Session::manifest).
    ///
    /// # Errors
    /// Fails when the plugin's ABI does not match, when its name is already
    /// registered, or when one of its registrations is refused.
    pub fn register(
        &mut self,
        plugin: Box<dyn Plugin>,
    ) -> Result<PluginId, PluginError> {
        let manifest = plugin.manifest();
        if manifest.abi != crate::kernel::signature::PLUGIN_ABI {
            return Err(PluginError::AbiMismatch {
                expected: crate::kernel::signature::PLUGIN_ABI,
                got: manifest.abi,
            });
        }
        let plugin: Arc<dyn Plugin> = Arc::from(plugin);
        let mut inner = self.inner.borrow_mut();
        if inner
            .plugins
            .iter()
            .any(|p| p.manifest().name == manifest.name)
        {
            return Err(PluginError::Duplicate { name: manifest.name });
        }
        // Register into the store graph so that names resolve for terms the
        // user builds directly.
        {
            let mut scratch = crate::graph::Program::default();
            let mut registry = Registry::new(&mut inner.graph, &mut scratch, &manifest.name);
            plugin.register(&mut registry)?;
        }
        inner.plugins.push(Arc::clone(&plugin));
        inner.plugin_count = inner.plugin_count.saturating_add(1);
        // Rebuild every template so that the plugin's rules take part in the
        // search, not only its operators.
        let plugins = inner.plugins.clone();
        let mut rebuilt = Vec::with_capacity(inner.templates.len());
        for template in &inner.templates {
            // Carry the session's units forward. `Template::new` starts from
            // the process-wide table, so a host unit registered before the
            // plugin would silently drop out of the manifest.
            let units = crate::api::SessionUnits::from_manifest(template.engine.manifest());
            rebuilt.push(crate::api::Template::with_units(
                &template.rule_sets,
                &plugins,
                units,
            )?);
        }
        inner.templates = rebuilt;
        Ok(PluginId(inner.plugin_count.saturating_sub(1)))
    }

    /// An independent session with the same configuration and the same
    /// terms.
    ///
    /// The two sessions can be used from different threads and merged back
    /// with [`Session::merge`].
    #[must_use]
    pub fn fork(&self) -> Self {
        let inner = self.inner.borrow();
        let templates: Vec<crate::api::Template> = inner
            .templates
            .iter()
            .map(|template| {
                crate::api::Template {
                    sets: template.sets.clone(),
                    rule_sets: template.rule_sets.clone(),
                    graph: template.graph.clone(),
                    engine: template.engine.clone(),
                }
            })
            .collect();
        Self {
            id: crate::api::next_session_id(),
            inner: std::cell::RefCell::new(crate::api::Inner {
                graph: inner.graph.clone(),
                templates,
                plugins: inner.plugins.clone(),
                plugin_count: inner.plugin_count,
            }),
        }
    }

    /// Imports everything `other` knows into this session.
    ///
    /// Every term is copied in, and every equality the other session derived
    /// is re-asserted here, so the merged session knows at least as much as
    /// either input. Terms keep their identity: a term that both sessions
    /// have becomes one node, not two.
    ///
    /// # Errors
    /// Fails when the two sessions disagree about the dimension of a symbol.
    pub fn merge(
        &mut self,
        other: &Self,
    ) -> Result<(), MergeError> {
        if std::ptr::eq(self, other) {
            return Ok(());
        }
        let mut inner = self.inner.borrow_mut();
        let source = other.inner.borrow();
        // 1. Structure: import every node, remembering where it landed.
        let count = source.graph.len();
        let mut images: Vec<NodeId> = Vec::with_capacity(count);
        for index in 0..count {
            let node = NodeId::from_raw(u32::try_from(index).unwrap_or(u32::MAX));
            images.push(inner.graph.import(&source.graph, node));
        }
        // 2. Equalities: replay every class that holds more than one member.
        let mut replayed = 0_usize;
        for index in 0..count {
            let node = NodeId::from_raw(u32::try_from(index).unwrap_or(u32::MAX));
            let class = source.graph.find(node);
            if class.0 != node.raw() {
                continue;
            }
            let members: Vec<NodeId> = source.graph.members(class).collect();
            let Some((&first, rest)) = members.split_first() else {
                continue;
            };
            let Some(&target) = images.get(first.index()) else {
                continue;
            };
            for member in rest {
                if let Some(&image) = images.get(member.index())
                    && inner.graph.union(target, image)
                {
                    replayed = replayed.saturating_add(1);
                }
            }
        }
        inner.graph.rebuild();
        // 3. Dimension declarations, translated between arenas.
        let source_arena = source.graph.dim_arena();
        let translations: Vec<(String, crate::graph::DimTerm)> = source
            .graph
            .dim_assumptions
            .iter()
            .map(|(symbol, term)| {
                (
                    source.graph.interner().symbol_name(*symbol).to_owned(),
                    term.clone(),
                )
            })
            .collect();
        drop(source_arena);
        for (name, term) in translations {
            // Fallible on purpose: `from_term` is total but answers `Any`
            // when the term cannot be expressed, and translating a dimension
            // assumption into "any dimension" is a silently weakened merge.
            let expression = crate::dim::DimExpr::try_from_term(&term, &source.graph.dim_arena())
                .map_err(|error| MergeError::Translation(error.to_string()))?;
            let translated = {
                let mut arena = match inner.graph.dims.write() {
                    | Ok(guard) => guard,
                    | Err(poisoned) => poisoned.into_inner(),
                };
                expression
                    .to_term(&mut arena)
                    .map_err(|error| MergeError::Translation(error.to_string()))?
            };
            let symbol = inner.graph.interner_mut().symbol(&name);
            if let Some(existing) = inner.graph.dim_assumptions.get(&symbol) {
                let here = existing.describe(&inner.graph.dim_arena());
                let there = translated.describe(&inner.graph.dim_arena());
                if here != there {
                    return Err(MergeError::Dimension {
                        symbol: name,
                        left: here,
                        right: there,
                    });
                }
            } else {
                inner.graph.dim_assumptions.insert(symbol, translated);
            }
        }
        let _ = replayed;
        Ok(())
    }
}

impl Kernel for Session {
    fn call(
        &self,
        op: &OpKey,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        // The key pins the overload, so this must not fall back to dispatch by
        // name: a caller that resolved once would silently get a different
        // operator if the table changed underneath it.
        self.call_key(op, args, ctx)
    }

    fn call_ho(
        &self,
        op: &OpKey,
        args: &[Value],
        fns: &[CallbackRef],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let ctx = clone_ctx(ctx);
        check_callbacks(&ctx, fns)?;
        let registration = self.op_for_key(op, args)?;
        self.call_op(registration, &op.name, args, &ctx)
    }

    fn register(
        &mut self,
        plugin: Box<dyn Plugin>,
    ) -> Result<PluginId, PluginError> {
        Self::register(self, plugin)
    }
}

/// Copies a call context so that it can be mutated while it is used.
fn clone_ctx<'a>(ctx: &CallCtx<'a>) -> CallCtx<'a> {
    CallCtx {
        target: ctx.target,
        budget: ctx.budget.clone(),
        dimension_mode: ctx.dimension_mode,
        assumptions: ctx.assumptions.clone(),
        bindings: ctx.bindings.clone(),
        substitutions: ctx.substitutions.clone(),
        callbacks: ctx.callbacks.clone(),
        // Carried through, not dropped. It used to be `None`, so every public
        // entry point silently discarded the caller's collector and
        // `CallCtx::with_trace` recorded nothing at all.
        trace: ctx.trace,
    }
}

/// The configuration a dynamic call runs with.
fn config_for(
    sets: &[crate::graph::RuleSet],
    ctx: &CallCtx<'_>,
) -> Config {
    let mut config = Config::empty();
    for set in sets {
        config = config.with(set.clone());
    }
    config = match ctx.target {
        | Target::Symbolic => config.symbolic(),
        | Target::Numeric { tolerance } => config.numeric(tolerance),
    };
    config = config
        .budget(ctx.budget.clone())
        .dimensions(ctx.dimension_mode);
    for (symbol, value) in &ctx.bindings {
        config = config.bind(symbol, *value);
    }
    for assumption in &ctx.assumptions {
        config = config.assume(&assumption.symbol, assumption.facts);
    }
    config
}

/// Maps a computation error onto the ABI's error type.
fn op_error(error: crate::api::ComputeError) -> OpError {
    match error {
        | crate::api::ComputeError::UnknownOperator(name) => {
            OpError::UnknownOp(OpKey::new("", &name, None))
        },
        | crate::api::ComputeError::NotNumeric(what) => OpError::NotNumeric(what),
        | crate::api::ComputeError::Dimension(error) => OpError::Dimension(error),
        | crate::api::ComputeError::Inconsistent(what) => OpError::Internal(what),
        | crate::api::ComputeError::Parse(error) => OpError::Domain(error.to_string()),
        | crate::api::ComputeError::Rule(error) => OpError::Domain(error.to_string()),
    }
}

/// Checks that every callback handle names a callback in the context.
fn check_callbacks(
    ctx: &CallCtx<'_>,
    fns: &[CallbackRef],
) -> Result<(), OpError> {
    for handle in fns {
        if ctx.callbacks.get(handle.0 as usize).is_none() {
            return Err(OpError::Host(format!(
                "callback {} was not supplied in the context",
                handle.0
            )));
        }
    }
    Ok(())
}

/// Checks `args` against the signatures of `descriptor`.
///
/// An operator with no signature accepts anything. One with signatures accepts
/// only what at least one of them describes, and the error names the first
/// argument no overload takes — silently executing a `(Real) -> Real` overload
/// on a tensor is a wrong answer waiting to happen, not a conservative choice.
fn check_args(
    descriptor: &crate::graph::OpDescriptor,
    args: &[Value],
) -> Result<(), OpError> {
    let (signatures, _) = crate::kernel::infer::signatures(descriptor);
    if signatures.is_empty() {
        return Ok(());
    }
    if crate::kernel::overload::best(&signatures, args) != crate::kernel::overload::Best::None {
        return Ok(());
    }
    for (index, value) in args.iter().enumerate() {
        let accepted = signatures.iter().any(|signature| {
            signature
                .params
                .get(index)
                .or_else(|| signature.variadic.as_ref().map(|variadic| &variadic.elem))
                .is_some_and(|param| crate::kernel::overload::matches(param, value).accepted())
        });
        if accepted {
            continue;
        }
        let expected = signatures
            .iter()
            .filter_map(|signature| {
                signature
                    .params
                    .get(index)
                    .or_else(|| signature.variadic.as_ref().map(|variadic| &variadic.elem))
            })
            .find(|param| {
                crate::kernel::overload::matches(param, value) == crate::kernel::overload::Match::No
            })
            .cloned()
            .unwrap_or_else(crate::kernel::infer::term);
        return Err(OpError::Type {
            expected: Box::new(expected),
            got: Box::new(value.type_expr()),
            arg: index,
        });
    }
    // The arity does not fit any signature (a variadic range, say).
    Err(OpError::Arity {
        expected: ArityInfo::from_arity(descriptor.arity),
        got: args.len(),
    })
}

/// The registration that `args` should reach.
///
/// Arity filters first, then the declared signatures decide by type, most
/// specific first. When none of the candidates declares a signature the first
/// arity match wins, which is the behaviour an untyped operator has always
/// had; a tie between typed candidates is refused rather than guessed.
fn resolve_op(
    graph: &Graph,
    name: &str,
    args: &[Value],
) -> Result<OpId, OpError> {
    let candidates = graph.ops().lookup_all(name);
    let Some(&first) = candidates.first() else {
        return Err(OpError::UnknownOp(OpKey::new("", name, None)));
    };
    let fitting: Vec<OpId> = candidates
        .iter()
        .copied()
        .filter(|&id| graph.ops().get(id).arity.accepts(args.len()))
        .collect();
    let Some(&only) = fitting.first() else {
        return Err(OpError::Arity {
            expected: ArityInfo::from_arity(graph.ops().get(first).arity),
            got: args.len(),
        });
    };
    if fitting.len() == 1 {
        check_args(graph.ops().get(only), args)?;
        return Ok(only);
    }
    let mut scored: Vec<(OpId, u32)> = Vec::new();
    for &id in &fitting {
        let descriptor = graph.ops().get(id);
        let (signatures, _) = crate::kernel::infer::signatures(descriptor);
        let best = signatures
            .iter()
            .filter_map(|signature| crate::kernel::overload::score(signature, args))
            .max();
        if let Some(score) = best {
            scored.push((id, score));
        }
    }
    if scored.is_empty() {
        // No candidate's signature fits. That is an error unless none of them
        // declares one at all, in which case arity is all there is to go on.
        if fitting.iter().all(|&id| {
            crate::kernel::infer::signatures(graph.ops().get(id))
                .0
                .is_empty()
        }) {
            return Ok(only);
        }
        check_args(graph.ops().get(only), args)?;
        return Ok(only);
    }
    let top = scored.iter().map(|(_, score)| *score).max().unwrap_or(0);
    let mut tied = scored
        .iter()
        .filter(|(_, score)| *score == top)
        .map(|(id, _)| *id);
    let Some(winner) = tied.next() else {
        return Ok(only);
    };
    let count = 1 + tied.count();
    if count == 1 {
        check_args(graph.ops().get(winner), args)?;
        Ok(winner)
    } else {
        Err(OpError::Ambiguous {
            name: name.to_owned(),
            count,
        })
    }
}

/// The value of a node.
///
/// Infallible: a node that is not a number, a literal, a symbol, a quantity
/// or a list of numbers is a term, and a term is a value.
fn value_of_node(
    graph: &Graph,
    node: NodeId,
    session: SessionId,
) -> Value {
    if let Some(number) = graph.number_of(node) {
        return match number {
            | Number::Int(value) => Value::Int(value.clone()),
            | Number::Rat(_) | Number::Float(_) => Value::Real(number.clone()),
        };
    }
    match graph.payload(node) {
        | Some(Payload::Bool(value)) if graph.op(node) == core::LIT => {
            return Value::Bool(*value);
        },
        | Some(Payload::Str(text)) if graph.op(node) == core::LIT => {
            return Value::String(text.to_string());
        },
        | _ => {},
    }
    if let Some(symbol) = graph.as_symbol(node) {
        return Value::Symbol(graph.interner().symbol_name(symbol).to_owned());
    }
    // A quantity: value plus a unit the registry knows.
    if let Some(quantity) = quantity_of(graph, node) {
        return Value::Quantity(quantity);
    }
    // A rectangular list of numbers.
    if let Some(tensor) = tensor_of(graph, node) {
        return Value::Tensor(tensor);
    }
    Value::Term(TermRef {
        session,
        node: node.raw(),
    })
}

/// The quantity a `quantity(value, unit)` node denotes.
fn quantity_of(
    graph: &Graph,
    node: NodeId,
) -> Option<QuantityValue> {
    let op = graph.ops().lookup("quantity")?;
    if graph.op(node) != op {
        return None;
    }
    let children = graph.children(node).to_vec();
    let (&value, unit_node) = children.split_first()?;
    let &unit_node = unit_node.first()?;
    let text = unit_text(graph, unit_node)?;
    let resolved = crate::dim::units::resolve(&text).ok()?;
    let unit = crate::dim::UnitDef::linear(&text, resolved.dimension, resolved.scale).to_value();
    Some(QuantityValue {
        value: Box::new(value_of_node(graph, value, SessionId(0))),
        unit,
    })
}

/// The unit a node spells, as text.
fn unit_text(
    graph: &Graph,
    node: NodeId,
) -> Option<String> {
    if let Some(symbol) = graph.as_symbol(node) {
        return Some(graph.interner().symbol_name(symbol).to_owned());
    }
    if let Some(number) = graph.number_of(node)
        && number.is_one()
    {
        return Some("1".to_owned());
    }
    let op = graph.op(node);
    let children = graph.children(node).to_vec();
    match op {
        | core::MUL => {
            let parts: Option<Vec<String>> = children
                .iter()
                .map(|&child| unit_text(graph, child))
                .collect();
            Some(parts?.join("*"))
        },
        | core::POW => {
            let (&base, rest) = children.split_first()?;
            let &exponent = rest.first()?;
            let name = unit_text(graph, base)?;
            let power = graph.number_of(exponent)?.to_rational()?;
            Some(format!("{name}^{}/{}", power.numer(), power.denom()))
        },
        | _ => None,
    }
}

/// The dense tensor a nested `list` denotes.
///
/// The shape is the full shape, not the outermost length: a `list` of `list`s
/// is rank 2. Deriving only the outermost extent made every matrix-shaped
/// result fail `is_consistent`, so `identity(3)` came back as an opaque `Term`
/// and every matrix-returning operator violated its declared type.
///
/// Ragged input is rejected rather than truncated: `[[1, 2], [3]]` has no
/// shape, and silently padding it would be a wrong answer.
fn tensor_of(
    graph: &Graph,
    node: NodeId,
) -> Option<TensorValue> {
    /// The shape of a nested list, or `None` when it is ragged.
    fn shape_of(
        graph: &Graph,
        node: NodeId,
    ) -> Option<Vec<usize>> {
        if graph.number_of(node).is_some() {
            return Some(Vec::new());
        }
        if graph.op(node) != core::LIST {
            return None;
        }
        let items = graph.children(node).to_vec();
        let Some((&first, rest)) = items.split_first() else {
            // An empty list is a rank-1 empty tensor.
            return Some(vec![0]);
        };
        let tail = shape_of(graph, first)?;
        for &item in rest {
            if shape_of(graph, item)? != tail {
                return None;
            }
        }
        let mut shape = vec![items.len()];
        shape.extend(tail);
        Some(shape)
    }
    fn flatten(
        graph: &Graph,
        node: NodeId,
        data: &mut Vec<Value>,
    ) {
        if let Some(number) = graph.number_of(node) {
            data.push(match number {
                | Number::Int(value) => Value::Int(value.clone()),
                | other => Value::Real(other.clone()),
            });
            return;
        }
        for &item in &graph.children(node).to_vec() {
            flatten(graph, item, data);
        }
    }
    if graph.op(node) != core::LIST {
        return None;
    }
    let shape = shape_of(graph, node)?;
    let mut data = Vec::new();
    flatten(graph, node, &mut data);
    let elem = data
        .first()
        .map_or(crate::kernel::signature::TypeExpr::Real, Value::type_expr);
    let tensor = TensorValue { shape, data, elem };
    tensor.is_consistent().then_some(tensor)
}

/// The term a value denotes, interned into `graph`.
///
/// A unit is deliberately **not** convertible here. A [`Value::Unit`] is a
/// definition — a name, a dimension and a scale — not a number, and the term
/// graph has no node that carries those three things. Spelling it as a bare
/// symbol (`m`) would be a lie twice over: the graph would hold a free
/// variable that `value_of` reads back as a dimensionless
/// [`Value::Symbol`], and an operator receiving it would see an unknown
/// scalar where a unit was meant. So a bare unit is refused, and the caller
/// must pass [`Value::Quantity`] (a number *with* that unit).
///
/// The one legitimate unit position is the second argument of `quantity`
/// (`convert`'s target unit is the same idea). There the units rules parse a
/// *unit term*, which is exactly a symbol or a product of symbols, so
/// [`Session::call_op`] encodes a unit argument that way through
/// [`term_of_argument`] before reaching this function.
fn term_of_value(
    graph: &mut Graph,
    value: &Value,
    session: SessionId,
) -> Result<NodeId, OpError> {
    match value {
        | Value::Int(value) => Ok(graph.num(Number::Int(value.clone()))),
        | Value::Real(number) => Ok(graph.num(number.clone())),
        | Value::Float(value) => Ok(graph.float(*value)),
        | Value::Bool(value) => Ok(graph.lit(Payload::Bool(*value))),
        | Value::String(text) => Ok(graph.lit(Payload::Str(Arc::from(text.as_str())))),
        | Value::Symbol(name) => Ok(graph.sym(name)),
        | Value::Unit(unit) => {
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: format!(
                    "a bare unit definition (`{}`): a unit is a definition, not a number, \
                 so it has no term of its own; pass `quantity(value, {})` instead",
                    unit.name, unit.name
                ),
            })
        },
        | Value::Quantity(quantity) => {
            let op = graph.ops().lookup("quantity").ok_or_else(|| {
                OpError::Unsupported {
                    backend: "term".to_owned(),
                    feature: "quantities (the `units` rule set is not installed)".to_owned(),
                }
            })?;
            let value = term_of_value(graph, &quantity.value, session)?;
            let unit = graph.sym(&quantity.unit.name);
            graph
                .try_node(op, &[value, unit])
                .ok_or_else(|| OpError::Internal("quantity is not binary".to_owned()))
        },
        | Value::Complex(complex) => {
            let (re, im) = complex.to_f64();
            if im == 0.0 {
                return Ok(graph.float(re));
            }
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: "complex literals".to_owned(),
            })
        },
        | Value::Tensor(tensor) => {
            if !tensor.is_consistent() {
                return Err(OpError::Domain(format!(
                    "tensor shape {:?} does not match {} elements",
                    tensor.shape,
                    tensor.data.len()
                )));
            }
            if tensor.shape.is_empty() {
                return match tensor.data.first() {
                    | Some(value) => term_of_value(graph, value, session),
                    | None => Err(OpError::Domain("empty scalar tensor".to_owned())),
                };
            }
            let extent = *tensor.shape.first().unwrap_or(&0);
            let stride = tensor.data.len().checked_div(extent).unwrap_or(0);
            let mut items = Vec::with_capacity(extent);
            for index in 0..extent {
                let start = index.saturating_mul(stride);
                let end = start.saturating_add(stride);
                let slice = tensor.data.get(start..end).ok_or_else(|| {
                    OpError::Domain("tensor data is shorter than its shape".to_owned())
                })?;
                let inner = TensorValue {
                    shape: tensor.shape.get(1..).unwrap_or(&[]).to_vec(),
                    data: slice.to_vec(),
                    elem: tensor.elem.clone(),
                };
                items.push(term_of_value(graph, &Value::Tensor(inner), session)?);
            }
            let list = graph.ops().lookup("list").unwrap_or(core::LIST);
            Ok(graph.node(list, &items))
        },
        | Value::Term(handle) => {
            if handle.session != session {
                return Err(OpError::ForeignTerm {
                    expected: handle.session,
                    got: session,
                });
            }
            let node = NodeId::from_raw(handle.node);
            // `is_none` rejects only 0. A handle from another graph of the
            // same session is a valid id with no node here, and indexing it
            // would panic inside the arena rather than at the boundary.
            if node.is_none() || !graph.contains(node) {
                return Err(OpError::Internal(format!("no term {}", handle.node)));
            }
            Ok(node)
        },
        | Value::Dyn(dynamic) => term_of_value(graph, &dynamic.value, session),
        | Value::Sparse(_) => {
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: "sparse literals".to_owned(),
            })
        },
        | Value::Record(_) | Value::Variant(_) => {
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: "records and variants".to_owned(),
            })
        },
        | Value::Function(_) => {
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: "host callbacks".to_owned(),
            })
        },
        | Value::Opaque(_) => {
            Err(OpError::Unsupported {
                backend: "term".to_owned(),
                feature: "opaque payloads".to_owned(),
            })
        },
    }
}

/// The arguments of `descriptor` that are declared as a unit.
///
/// A unit parameter is the one place a unit *definition* belongs in a term:
/// `quantity(value, unit)` declares its second parameter `Unit(D)`, and
/// `convert` its target unit the same way. The graph spells such an argument
/// as the unit's symbol, because that is exactly what the units rules parse.
fn unit_parameters(
    descriptor: &crate::graph::OpDescriptor,
    count: usize,
) -> Vec<bool> {
    use crate::kernel::signature::TypeExpr;

    let (signatures, _) = crate::kernel::infer::signatures(descriptor);
    (0..count)
        .map(|index| {
            signatures.iter().any(|signature| {
                signature
                    .params
                    .get(index)
                    .is_some_and(|param| matches!(param, TypeExpr::Unit(_) | TypeExpr::Point(_)))
            })
        })
        .collect()
}

/// The term an argument denotes, honouring a declared unit parameter.
///
/// `unit_position` comes from [`unit_parameters`]; only there is a
/// [`Value::Unit`] turned into the symbol the graph uses for a unit. Every
/// other position goes through [`term_of_value`], which refuses a bare unit.
fn term_of_argument(
    graph: &mut Graph,
    value: &Value,
    session: SessionId,
    unit_position: bool,
) -> Result<NodeId, OpError> {
    if unit_position && let Value::Unit(unit) = value {
        return Ok(graph.sym(&unit.name));
    }
    term_of_value(graph, value, session)
}

impl Term<'_> {
    /// The operator of this term.
    #[must_use]
    pub fn op(self) -> OpId {
        self.session.inner.borrow().graph.op(self.node)
    }

    /// The operator's name.
    #[must_use]
    pub fn op_name(self) -> Arc<str> {
        let inner = self.session.inner.borrow();
        Arc::clone(&inner.graph.ops().get(inner.graph.op(self.node)).name)
    }

    /// The immediate children of this term.
    #[must_use]
    pub fn children(self) -> Vec<Self> {
        let inner = self.session.inner.borrow();
        inner
            .graph
            .children(self.node)
            .iter()
            .map(|&node| {
                Self {
                    session: self.session,
                    node,
                }
            })
            .collect()
    }

    /// Replaces every free occurrence of the terms in `map` inside this term.
    ///
    /// Keys must be symbols, and replacements are matched by symbol, so the
    /// result is the same term with those symbols renamed to whatever the
    /// replacement denotes.
    #[must_use]
    pub fn substitute(
        self,
        map: &[(Term<'_>, Term<'_>)],
    ) -> Self {
        let mut inner = self.session.inner.borrow_mut();
        let mut node = self.node;
        for (from, to) in map {
            node = inner.graph.substitute(node, from.node, to.node);
        }
        drop(inner);
        Self {
            session: self.session,
            node,
        }
    }

    /// The names of the symbols this term may depend on, sorted.
    #[must_use]
    pub fn free_symbols(self) -> Vec<String> {
        let inner = self.session.inner.borrow();
        inner
            .graph
            .free_symbols(inner.graph.find(self.node))
            .iter()
            .map(|&symbol| inner.graph.interner().symbol_name(symbol).to_owned())
            .collect()
    }

    /// The value of this term (see [`Session::value_of`]).
    ///
    /// # Errors
    /// See [`Session::value_of`].
    pub fn to_value(self) -> Result<Value, OpError> {
        self.session.value_of(self)
    }

    /// An owned handle to this term.
    #[must_use]
    pub const fn handle(self) -> TermRef {
        self.session.handle(self)
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;
    use crate::api::Config;
    use crate::graph::OpDescriptor;

    #[test]
    fn a_nested_list_is_a_tensor_of_the_full_rank() {
        // The shape used to be the outermost `list` length only, so a rank-2
        // value failed `is_consistent` and came back as an opaque `Term`:
        // every matrix-returning operator then violated its declared type.
        let s = Session::new();
        let node = s
            .parse("list(list(1, 2), list(3, 4))")
            .unwrap_or_else(|e| panic!("{e}"));
        let value = s.value_of(node).unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(tensor) = value else {
            panic!("expected a rank-2 tensor, got {value:?}");
        };
        assert_eq!(tensor.shape, vec![2, 2]);
        assert_eq!(tensor.data.len(), 4);

        // Rank 3, and the row-major order is preserved.
        let cube = s
            .parse("list(list(list(1, 2), list(3, 4)), list(list(5, 6), list(7, 8)))")
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(tensor) = s.value_of(cube).unwrap_or_else(|e| panic!("{e}")) else {
            panic!("expected a rank-3 tensor");
        };
        assert_eq!(tensor.shape, vec![2, 2, 2]);
        assert_eq!(tensor.data[0], Value::Int(BigInt::from(1)));
        assert_eq!(tensor.data[7], Value::Int(BigInt::from(8)));
    }

    #[test]
    fn a_ragged_list_is_not_a_tensor() {
        // `[[1, 2], [3]]` has no shape. Padding or truncating it would be a
        // wrong answer, so it stays a term.
        let s = Session::new();
        let node = s
            .parse("list(list(1, 2), list(3))")
            .unwrap_or_else(|e| panic!("{e}"));
        let value = s.value_of(node).unwrap_or_else(|e| panic!("{e}"));
        assert!(
            matches!(value, Value::Term(_)),
            "a ragged list must not become a tensor, got {value:?}"
        );
    }

    #[test]
    fn numbers_symbols_and_tensors_cross_the_boundary() {
        let s = Session::new();
        let x = s.sym("x");
        assert_eq!(
            s.value_of(s.int(7)).unwrap_or_else(|e| panic!("{e}")),
            Value::Int(BigInt::from(7))
        );
        assert_eq!(
            s.value_of(x).unwrap_or_else(|e| panic!("{e}")),
            Value::Symbol("x".to_owned())
        );
        let list = s.parse("list(1, 2, 3)").unwrap_or(x);
        let value = s.value_of(list).unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(tensor) = value else {
            panic!("expected a tensor, got {value:?}");
        };
        assert_eq!(tensor.shape, vec![3]);
        assert_eq!(tensor.data.len(), 3);
    }

    #[test]
    fn a_term_survives_a_round_trip_through_a_value() {
        let s = Session::new();
        let term = s.parse("x + 1/2").unwrap_or_else(|e| panic!("{e}"));
        let value = s.value_of(term).unwrap_or_else(|e| panic!("{e}"));
        let back = s.term_of(&value).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, term);
    }

    #[test]
    fn dynamic_calls_reach_the_engine() {
        let s = Session::new();
        let x = Value::Symbol("x".to_owned());
        let out = s
            .call_value("sin", &[x], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        // `sin(x)` does not reduce, so it comes back as a handle to a term of
        // this session rather than as a number.
        let Value::Term(handle) = out else {
            panic!("expected a term handle, got {out:?}");
        };
        assert_eq!(
            s.term_ref(handle)
                .unwrap_or_else(|e| panic!("{e}"))
                .to_string(),
            "sin(x)"
        );
        let two = Value::int(2_i32);
        let sum = s
            .call_value("add", &[Value::int(1_i32), two], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(sum, Value::int(3_i32));
        assert!(matches!(
            s.call_value("no_such_op", &[], &CallCtx::new()),
            Err(OpError::UnknownOp(_))
        ));
        assert!(matches!(
            s.call_value("sin", &[], &CallCtx::new()),
            Err(OpError::Arity { .. })
        ));
    }

    #[test]
    fn a_dynamic_call_runs_the_rules() {
        let s = Session::new();
        let out = s
            .call_value(
                "diff",
                &[Value::Symbol("x".to_owned()), Value::Symbol("x".to_owned())],
                &CallCtx::new(),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(out, Value::int(1_i32));
    }

    #[test]
    fn handles_are_checked_against_their_session() {
        let a = Session::new();
        let b = Session::new();
        let handle = a.handle(a.sym("x"));
        assert!(a.term_ref(handle).is_ok());
        assert!(matches!(
            b.term_ref(handle),
            Err(OpError::ForeignTerm { .. })
        ));
    }

    #[test]
    fn fork_and_merge_agree_with_serial_evaluation() {
        let s = Session::new();
        let x = s.sym("x");
        let request = crate::api::diff(x.pow(3) * crate::api::sin(x), x);
        let serial = s
            .compute(request, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"))
            .term
            .to_string();

        let mut a = s.fork();
        let b = s.fork();
        // Each fork has its own graph, so each needs its own symbol; a handle
        // from one session is refused by the other by design.
        let ra = {
            let ax = a.sym("x");
            a.compute(crate::api::diff(ax.pow(3), ax), &Config::new())
                .unwrap_or_else(|e| panic!("{e}"))
                .term
                .to_string()
        };
        let rb = {
            let bx = b.sym("x");
            b.compute(crate::api::diff(bx.pow(2), bx), &Config::new())
                .unwrap_or_else(|e| panic!("{e}"))
                .term
                .to_string()
        };
        assert_eq!(ra, "3*x^2");
        assert_eq!(rb, "2*x");

        a.merge(&b).unwrap_or_else(|e| panic!("{e}"));

        // The merged session knows both answers and still computes.
        let merged = {
            let root = a
                .parse("diff(x^3, x) + diff(x^2, x)")
                .unwrap_or_else(|e| panic!("{e}"));
            a.compute(root, &Config::new())
                .unwrap_or_else(|e| panic!("{e}"))
                .term
                .to_string()
        };
        assert_eq!(merged, "3*x^2 + 2*x");
        assert_eq!(serial, "x^3*cos(x) + 3*x^2*sin(x)");
    }

    #[test]
    fn a_symbol_is_a_type_of_its_own() {
        // `Symbol` is not `String`: `diff` and `solve` need a hole, not text.
        assert_eq!(
            Value::Symbol("x".to_owned()).type_expr(),
            crate::kernel::signature::TypeExpr::Symbol
        );
        assert_ne!(
            Value::Symbol("x".to_owned()).type_expr(),
            crate::kernel::signature::TypeExpr::String
        );
    }

    #[test]
    fn an_overload_can_be_resolved_and_called_by_key() {
        let s = Session::new();
        let args = [Value::Float(0.5)];
        let key = s
            .resolve_overload("sin", &args)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(key.name, "sin");
        assert_ne!(key.signature_hash, [0_u8; 16], "sin is typed");
        let by_name = s
            .call_value("sin", &args, &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let by_key = s
            .call_key(&key, &args, &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(by_name, by_key);
        // The key is exactly what `overloads` lists.
        assert!(s.overloads("sin").iter().any(|(k, _)| *k == key));
    }

    #[test]
    fn a_key_that_names_nothing_is_refused() {
        let s = Session::new();
        let mut key = s
            .resolve_overload("sin", &[Value::Float(0.5)])
            .unwrap_or_else(|e| panic!("{e}"));
        key.signature_hash = [0xAB; 16];
        assert!(matches!(
            s.call_key(&key, &[Value::Float(0.5)], &CallCtx::new()),
            Err(OpError::UnknownOp(_))
        ));
    }

    #[test]
    fn every_operator_has_a_signature_and_none_is_unexplained() {
        // The old contract allowed an operator with no type at all. The
        // derivation in `kernel::infer` covers every non-leaf, so the only
        // way to have none is to be a leaf, and a leaf says so.
        let s = Session::new();
        let manifest = s.manifest();
        assert!(manifest.ops.len() > 1000);
        let mut typed = 0_usize;
        for op in &manifest.ops {
            if op.overloads.is_empty() {
                assert!(
                    op.dynamic_reason.is_some(),
                    "`{}` is neither typed nor explained",
                    op.name
                );
                assert_eq!(op.keys.len(), 1);
                assert_eq!(op.keys[0].signature_hash, [0_u8; 16]);
            } else {
                typed += 1;
                assert_eq!(op.keys.len(), op.overloads.len());
                // `Session::overloads` agrees with the manifest.
                assert_eq!(s.overloads(&op.name).len(), op.overloads.len());
            }
        }
        assert_eq!(
            typed,
            manifest.ops.len(),
            "every operator in the standard sets must be typed"
        );
    }

    #[test]
    fn a_structured_dimension_accepts_a_base_the_engine_has_never_seen() {
        use crate::kernel::signature::DimExpr;
        let s = Session::new();
        let x = s.sym("x");
        let answer = s
            .compute(
                s.parse("x + x").unwrap_or(x),
                &Config::new().dimension_term("x", DimExpr::base("Currency")),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(answer.term.to_string(), "2*x");
    }

    #[test]
    fn a_plugin_becomes_part_of_the_session() {
        use crate::kernel::plugin::SimplePlugin;
        use crate::kernel::signature::Signature;
        use crate::kernel::signature::TypeExpr;

        let plugin = SimplePlugin::new("acme/twice", "0.1.0", |reg| {
            reg.register_op(
                OpDescriptor::new("acme_twice", crate::graph::Arity::Fixed(1))
                    .signature(Signature::scalars(1, TypeExpr::Real)),
            )?;
            reg.register_rewrites(
                crate::graph::Tier::Normalize,
                &["acme/twice-zero: acme_twice(0) => 0"],
            )?;
            Ok(())
        });
        let mut s = Session::new();
        let id = s
            .register(Box::new(plugin))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(id.0, 0);
        let manifest = s.manifest();
        assert!(manifest.op("acme_twice").is_some());
        assert!(manifest.plugins.iter().any(|p| p.name == "acme/twice"));
        let out = s
            .call_value("acme_twice", &[Value::int(0_i32)], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(out, Value::int(0_i32));
        assert!(
            s.register(Box::new(SimplePlugin::new("acme/twice", "0.1.0", |_| {
                Ok(())
            })))
            .is_err(),
            "a plugin may not be registered twice"
        );
    }

    /// A unit is a definition, not a number, so it has no term of its own.
    ///
    /// It used to become a bare symbol (`m`), which in the graph is a *free
    /// variable*. `value_of` then read that back as `Value::Symbol("m")` — a
    /// dimensionless value — so a unit passed as a term argument was silently
    /// reinterpreted as an unknown scalar.
    #[test]
    fn a_bare_unit_is_refused_in_the_term_path() {
        let s = Session::new();
        let metre = crate::dim::units::lookup("m")
            .unwrap_or_else(|| panic!("metre is in the table"))
            .to_value();
        let error = s
            .term_of(&Value::Unit(metre))
            .err()
            .unwrap_or_else(|| panic!("a bare unit must not enter the term graph"));
        assert!(
            matches!(error, OpError::Unsupported { .. }),
            "expected Unsupported, got {error:?}"
        );
        let message = error.to_string();
        assert!(
            message.contains("unit"),
            "the reason must name the unit: {message}"
        );
        assert!(
            message.contains("quantity"),
            "the reason must say what to pass instead: {message}"
        );
    }

    /// `quantity(3.0, m)` must keep working: the unit argument of `quantity`
    /// is the one place a unit definition belongs in a term, and there it is
    /// spelled as its symbol, which is what the units rules parse.
    #[test]
    fn a_unit_is_still_accepted_as_a_quantity_argument() {
        let s = Session::new();
        let metre = crate::dim::units::lookup("m")
            .unwrap_or_else(|| panic!("metre is in the table"))
            .to_value();
        let out = s
            .call_value(
                "quantity",
                &[Value::Float(3.0), Value::Unit(metre.clone())],
                &CallCtx::new(),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Quantity(quantity) = out else {
            panic!("expected a quantity, got {out:?}");
        };
        assert_eq!(quantity.value.as_f64(), Some(3.0));
        assert_eq!(quantity.unit.name, "m");
        assert_eq!(quantity.unit.dimension, metre.dimension);

        // The rule is by declared parameter, not by operator name, so
        // `convert`'s target unit keeps working too.
        let converted = s
            .call_value(
                "convert",
                &[
                    Value::Quantity(QuantityValue {
                        value: Box::new(Value::Float(3.0)),
                        unit: metre.clone(),
                    }),
                    Value::Unit(metre.clone()),
                ],
                &CallCtx::new(),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            matches!(converted, Value::Quantity(_)),
            "expected a quantity, got {converted:?}"
        );
    }

    /// A `Value::Quantity` still becomes `quantity(value, unit)` and comes
    /// back with its dimension: the refusal above must not reach it.
    #[test]
    fn a_quantity_keeps_its_unit_across_the_boundary() {
        let s = Session::new();
        let metre = crate::dim::units::lookup("m")
            .unwrap_or_else(|| panic!("metre is in the table"))
            .to_value();
        let value = Value::Quantity(QuantityValue {
            value: Box::new(Value::Float(3.0)),
            unit: metre.clone(),
        });
        let term = s.term_of(&value).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(term.to_string(), "quantity(3, m)");
        let back = s.value_of(term).unwrap_or_else(|e| panic!("{e}"));
        let Value::Quantity(quantity) = back else {
            panic!("expected a quantity, got {back:?}");
        };
        assert_eq!(quantity.unit.dimension, metre.dimension);
    }

    #[test]
    fn the_manifest_reports_the_standard_configuration() {
        let s = Session::new();
        let manifest = s.manifest();
        assert!(manifest.verify());
        assert!(manifest.ops.len() > 1000);
        assert_ne!(s.manifest_json(), "");
    }
}
