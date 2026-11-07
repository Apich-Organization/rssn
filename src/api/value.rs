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
//!
//! Numeric values are normalised on the way in. A host can build a
//! [`Number`] directly, so an integral rational can arrive as `Rat(0/1)`
//! rather than the `Int(0)` that [`Number::rat`] produces; `term_of_value`
//! collapses it (see [`Number::normalized`]) because the exact literal
//! rewrites (`cos(0) => 1`, `ln(1) => 0`, ...) match structurally and an
//! un-normalised zero would silently stop folding.

use std::fmt;
use std::sync::Arc;

use crate::api::Config;
use crate::api::Session;
use crate::api::Term;
use crate::dim::UnitDef;
use crate::dim::UnitRegistry;
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
use crate::kernel::value::UnitValue;
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
    /// A unit inside a quantity is resolved against *this session's*
    /// registry — the built-in table plus the units registered with
    /// [`Session::register_unit`] — so a session's own unit keeps its
    /// dimension and printed form instead of coming back as an opaque
    /// [`Value::Term`]. A compound unit expression that is not a single known
    /// unit (`m/s`) is described by its resolved dimension and scale.
    ///
    /// # Errors
    /// Fails only when a value cannot be represented, which today means a
    /// unit whose definition is inconsistent with the graph's arena.
    pub fn value_of(
        &self,
        term: Term<'_>,
    ) -> Result<Value, OpError> {
        let inner = self.inner.borrow();
        let (units, definitions): (super::SessionUnits, &[UnitValue]) =
            match inner.templates.first() {
                | Some(template) => {
                    let manifest = template.engine.manifest();
                    (
                        super::SessionUnits::from_manifest(manifest),
                        manifest.units.as_slice(),
                    )
                },
                | None => (super::SessionUnits::default(), &[]),
            };
        Ok(value_of_node(
            &inner.graph,
            NodeId::from_raw(term.id()),
            self.id,
            &units,
            definitions,
        ))
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
    /// [`Solutions`](crate::kernel::result::Solutions), which carries the
    /// unknown and whether the set is complete, and converts to a
    /// [`Value::Record`] for the dynamic ABI.
    ///
    /// What the solver does not state is not invented:
    ///
    /// * `multiplicities` is filled only for a univariate polynomial over the
    ///   rationals: the exact square-free decomposition says how many times
    ///   each rational root occurs, a square-free polynomial has only simple
    ///   roots, and the irrational root of a repeated factor is matched to its
    ///   factor on values, with a separation margin. Where the solver cannot
    ///   establish it (a transcendental equation such as `exp(t) = 2`, a
    ///   symbolic root of a repeated factor) the list stays empty: "not
    ///   reported" is not the claim "simple".
    /// * `conditions` is filled only with the conditions the answer states by
    ///   construction: a parameter in a denominator of a root — `-b/a` from
    ///   `a*x + b = 0` is valid only where `a != 0`. Realness and branch
    ///   conditions (a square root over a parameter, say) are not reported,
    ///   because the solver does not decide them here.
    ///
    /// For anything other than an exhaustive answer (`complete == false`) both
    /// lists stay empty: a single root found on the way cannot support a claim
    /// about multiplicity or validity.
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
        // field. A list whose elements are expressions rather than numbers
        // does not become a [`Value::Tensor`] (`tensor_of` only carries
        // numbers), so the reduced list node itself is expanded here: the
        // solver's own answer says the set is exhaustive, and dropping all but
        // one root would contradict it.
        let (roots, complete) = match value {
            | Value::Tensor(tensor) => (tensor.data, true),
            | Value::Sparse(sparse) => (sparse.data, true),
            | Value::Term(handle) => {
                let inner = self.inner.borrow();
                let node = NodeId::from_raw(handle.node);
                if inner.graph.op(node) == core::LIST {
                    let children = inner.graph.children(node).to_vec();
                    let (units, definitions): (super::SessionUnits, &[UnitValue]) =
                        match inner.templates.first() {
                            | Some(template) => {
                                let manifest = template.engine.manifest();
                                (
                                    super::SessionUnits::from_manifest(manifest),
                                    manifest.units.as_slice(),
                                )
                            },
                            | None => (super::SessionUnits::default(), &[]),
                        };
                    let roots = children
                        .iter()
                        .map(|&child| {
                            value_of_node(&inner.graph, child, self.id, &units, definitions)
                        })
                        .collect();
                    (roots, true)
                } else {
                    (vec![Value::Term(handle)], false)
                }
            },
            | other => (vec![other], false),
        };
        // A list is the solver's own enumeration, but for an equation with a
        // circular function of the unknown it is one representative per
        // principal branch: `sin(x) = 0` lists `0`, not `k pi`. That is not a
        // complete solution set, whatever the route to the roots.
        let (multiplicities, conditions, complete) = if complete {
            let mut inner = self.inner.borrow_mut();
            let expression = crate::rules::solve::as_expression(
                &mut inner.graph,
                NodeId::from_raw(equation.id()),
            );
            let unknown = NodeId::from_raw(variable.id());
            let periodic = crate::rules::solve::periodic_in(&inner.graph, expression, unknown);
            let (multiplicities, mut conditions) = crate::rules::solve::root_multiplicities(
                &mut inner.graph,
                expression,
                unknown,
                &roots,
            )
            .unwrap_or_default();
            for condition in crate::rules::solve::root_conditions(&inner.graph, unknown, &roots) {
                if !conditions.iter().any(|known| known.text == condition.text) {
                    conditions.push(condition);
                }
            }
            if periodic {
                conditions.push(crate::rules::solve::principal_branch_condition());
            }
            (multiplicities, conditions, !periodic)
        } else {
            (Vec::new(), Vec::new(), false)
        };
        Ok(crate::kernel::result::Solutions {
            variable: Some(name),
            roots,
            multiplicities,
            conditions,
            complete,
        })
    }

    /// Calls a higher-order operator with host callbacks.
    ///
    /// The callbacks are reachable from the operator's
    /// [`EvalDyn`](crate::graph::op::EvalDyn) through
    /// [`CallCtx::callbacks`]. Handles in `fns` are appended to `args` as
    /// function values, in the given order; ordinary signatures and evaluator
    /// argument positions remain identical. For large owned argument buffers,
    /// callers may append `Value::Function` themselves and use `call_value`
    /// with the same context, avoiding a clone. Rewrites
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
        self.resolve_overload_with_context(name, args, &CallCtx::new())
    }

    /// Resolve an overload using the callback bank's immutable signatures.
    /// Host type queries happen before borrowing the expression graph.
    ///
    /// # Errors
    /// Returns the errors of [`Self::resolve_overload`] or a callback-handle error.
    pub fn resolve_overload_with_context(
        &self,
        name: &str,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<OpKey, OpError> {
        let types = ctx.argument_types(args)?;
        let inner = self.inner.borrow();
        let graph = &inner.graph;
        let op = resolve_op_types(graph, name, &types)?;
        let descriptor = graph.ops().get(op);
        let (signatures, _) = crate::kernel::infer::signatures(descriptor);
        let chosen = match crate::kernel::overload::best_types(&signatures, &types) {
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
        let op = self.op_for_key(key, args, &ctx)?;
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
        ctx: &CallCtx<'_>,
    ) -> Result<OpId, OpError> {
        let types = ctx.argument_types(args)?;
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
                if crate::kernel::overload::score_types(pinned, &types).is_some() {
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
        let combined = callback_arguments(args, fns, &ctx)?;
        let args = combined.as_deref().unwrap_or(args);
        let types = ctx.argument_types(args)?;
        let op = {
            let inner = self.inner.borrow();
            resolve_op_types(&inner.graph, name, &types)?
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
        let combined = callback_arguments(args, fns, &ctx)?;
        let args = combined.as_deref().unwrap_or(args);
        let registration = self.op_for_key(op, args, &ctx)?;
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

/// Materialize trailing callback arguments while retaining their original bank IDs.
fn callback_arguments(
    args: &[Value],
    fns: &[CallbackRef],
    ctx: &CallCtx<'_>,
) -> Result<Option<Vec<Value>>, OpError> {
    check_callbacks(ctx, fns)?;
    if fns.is_empty() {
        return Ok(None);
    }
    let mut combined = Vec::with_capacity(args.len().saturating_add(fns.len()));
    combined.extend_from_slice(args);
    combined.extend(fns.iter().copied().map(Value::Function));
    Ok(Some(combined))
}

/// Checks `args` against the signatures of `descriptor`.
///
/// An operator with no signature accepts anything. One with signatures accepts
/// only what at least one of them describes, and the error names the first
/// argument no overload takes — silently executing a `(Real) -> Real` overload
/// on a tensor is a wrong answer waiting to happen, not a conservative choice.
fn check_arg_types(
    descriptor: &crate::graph::OpDescriptor,
    args: &[crate::kernel::signature::TypeExpr],
) -> Result<(), OpError> {
    let (signatures, _) = crate::kernel::infer::signatures(descriptor);
    if signatures.is_empty() {
        return Ok(());
    }
    if crate::kernel::overload::best_types(&signatures, args) != crate::kernel::overload::Best::None
    {
        return Ok(());
    }
    for (index, value) in args.iter().enumerate() {
        let accepted = signatures.iter().any(|signature| {
            signature
                .params
                .get(index)
                .or_else(|| signature.variadic.as_ref().map(|variadic| &variadic.elem))
                .is_some_and(|param| crate::kernel::overload::matches_type(param, value).accepted())
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
                crate::kernel::overload::matches_type(param, value)
                    == crate::kernel::overload::Match::No
            })
            .cloned()
            .unwrap_or_else(crate::kernel::infer::term);
        return Err(OpError::Type {
            expected: Box::new(expected),
            got: Box::new(value.clone()),
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
fn resolve_op_types(
    graph: &Graph,
    name: &str,
    args: &[crate::kernel::signature::TypeExpr],
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
        check_arg_types(graph.ops().get(only), args)?;
        return Ok(only);
    }
    let mut scored: Vec<(OpId, u32)> = Vec::new();
    for &id in &fitting {
        let descriptor = graph.ops().get(id);
        let (signatures, _) = crate::kernel::infer::signatures(descriptor);
        let best = signatures
            .iter()
            .filter_map(|signature| crate::kernel::overload::score_types(signature, args))
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
        check_arg_types(graph.ops().get(only), args)?;
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
        check_arg_types(graph.ops().get(winner), args)?;
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
/// or a list of numbers is a term, and a term is a value. `units` is the
/// session's registry and `definitions` the unit definitions its manifest
/// describes, so that a unit the session registered is resolved and carried
/// with the same definition the session's engine uses.
fn value_of_node(
    graph: &Graph,
    node: NodeId,
    session: SessionId,
    units: &super::SessionUnits,
    definitions: &[UnitValue],
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
    // A quantity: value plus a unit the session's registry knows.
    if let Some(quantity) = quantity_of(graph, node, session, units, definitions) {
        return Value::Quantity(quantity);
    }
    // A rectangular list of numbers.
    if let Some(tensor) = tensor_of(graph, node, session, units, definitions) {
        return Value::Tensor(tensor);
    }
    Value::Term(TermRef {
        session,
        node: node.raw(),
    })
}

/// The quantity a `quantity(value, unit)` node denotes.
///
/// The unit is resolved through the session's registry, not the process-wide
/// table: a unit registered with [`Session::register_unit`] is part of the
/// session, and reporting `quantity(3, furlong)` as an opaque term because
/// the shared table never heard of `furlong` lost the dimension of the value
/// the session itself declared. The full definition — display, aliases,
/// scale, offset — comes from the session's manifest; a compound expression
/// that is not a single known unit (`m/s`) keeps its resolved dimension and
/// scale.
fn quantity_of(
    graph: &Graph,
    node: NodeId,
    session: SessionId,
    units: &super::SessionUnits,
    definitions: &[UnitValue],
) -> Option<QuantityValue> {
    let op = graph.ops().lookup("quantity")?;
    if graph.op(node) != op {
        return None;
    }
    let children = graph.children(node).to_vec();
    let (&value, unit_node) = children.split_first()?;
    let &unit_node = unit_node.first()?;
    if let Some(Payload::Unit(unit)) = graph.payload(unit_node) {
        return Some(QuantityValue {
            value: Box::new(value_of_node(graph, value, session, units, definitions)),
            unit: unit.as_ref().clone(),
        });
    }
    // Compound units may contain anonymous definitions as well as names.
    // Resolve those structurally; rendering a literal's name and looking it
    // up would reintroduce the very metadata loss the literal prevents.
    if let Some(map) = crate::dim::units::parse_unit(graph, unit_node, 0)
        && map
            .keys()
            .any(|atom| matches!(atom, crate::dim::units::UnitAtom::Literal(_)))
    {
        let resolved = crate::dim::units::resolve_map(&map)?;
        let mut unit = UnitDef::linear(
            "compound",
            crate::dim::units::dim_to_expr(&resolved.dim),
            resolved.scale,
        );
        if let Some(offset) = resolved.offset {
            unit = unit.affine(offset);
        }
        return Some(QuantityValue {
            value: Box::new(value_of_node(graph, value, session, units, definitions)),
            unit: unit.to_value(),
        });
    }
    let text = unit_text(graph, unit_node)?;
    let resolved = units.resolve(&text).ok()?;
    let unit = definitions
        .iter()
        .find(|unit| unit.name == text || unit.aliases.contains(&text))
        .cloned()
        .unwrap_or_else(|| UnitDef::linear(&text, resolved.dimension, resolved.scale).to_value());
    Some(QuantityValue {
        value: Box::new(value_of_node(graph, value, session, units, definitions)),
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
    session: SessionId,
    units: &super::SessionUnits,
    definitions: &[UnitValue],
) -> Option<TensorValue> {
    fn scalar(
        graph: &Graph,
        node: NodeId,
        session: SessionId,
        units: &super::SessionUnits,
        definitions: &[UnitValue],
    ) -> Option<Value> {
        if let Some(number) = graph.number_of(node) {
            return Some(match number {
                | Number::Int(value) => Value::Int(value.clone()),
                | other => Value::Real(other.clone()),
            });
        }
        let quantity = quantity_of(graph, node, session, units, definitions)?;
        quantity.value.as_number()?;
        Some(Value::Quantity(quantity))
    }
    fn visit(
        graph: &Graph,
        node: NodeId,
        session: SessionId,
        units: &super::SessionUnits,
        definitions: &[UnitValue],
        data: &mut Vec<Value>,
    ) -> Option<Vec<usize>> {
        if let Some(value) = scalar(graph, node, session, units, definitions) {
            data.push(value);
            return Some(Vec::new());
        }
        if graph.op(node) != core::LIST {
            return None;
        }
        let items = graph.children(node);
        let Some((&first, rest)) = items.split_first() else {
            return Some(vec![0]);
        };
        let tail = visit(graph, first, session, units, definitions, data)?;
        for &item in rest {
            if visit(graph, item, session, units, definitions, data)? != tail {
                return None;
            }
        }
        let mut shape = vec![items.len()];
        shape.extend(tail);
        Some(shape)
    }
    if graph.op(node) != core::LIST {
        return None;
    }
    let mut data = Vec::new();
    let shape = visit(graph, node, session, units, definitions, &mut data)?;
    let mut elem = data
        .first()
        .map_or(crate::kernel::signature::TypeExpr::Real, Value::type_expr);
    for value in &data {
        if value.type_expr() == elem {
            continue;
        }
        // Mixed exact numeric literals form a real tensor. A list with
        // heterogeneous units remains a term rather than receiving the first
        // element's dimension as an invented uniform schema.
        if value.as_number().is_some()
            && matches!(
                elem,
                crate::kernel::signature::TypeExpr::Int | crate::kernel::signature::TypeExpr::Real
            )
        {
            elem = crate::kernel::signature::TypeExpr::Real;
        } else {
            return None;
        }
    }
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
        // A host can build the `Number` directly, and the variant is public:
        // the Aevia bridge sends an integral exact value as `Number::Rat(n/1)`
        // instead of the `Int(n)` that `Number::rat` produces. Interning that
        // non-normal form would leave a literal no rewrite recognises —
        // `cos(0)`, `exp(0)`, `ln(1)` and the rest match structurally — so an
        // exact limit endpoint would stop at its closed form
        // (`limit(sin(t)/t, 0/1)` answered `cos(0/1)`) instead of folding to
        // the exact value. Normalise at the boundary; a non-integral rational
        // and a float cross unchanged.
        | Value::Real(number) => Ok(graph.num(number.normalized())),
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
            let unit = term_of_unit(graph, &quantity.unit)?;
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
        | Value::Tensor(tensor) => term_of_tensor(graph, tensor, session, None),
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
        | Value::Dyn(dynamic) => {
            term_of_dimensioned_value(graph, &dynamic.value, &dynamic.dim, session)
        },
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

/// Materializes a unit definition structurally unless the registry already
/// contains this exact definition. Equal display names are never evidence.
fn term_of_unit(
    graph: &mut Graph,
    unit: &UnitValue,
) -> Result<NodeId, OpError> {
    crate::dim::units::resolve_value(unit).ok_or_else(|| OpError::Domain(
        "invalid unit definition: expected a concrete rational dimension, nonzero exact scale and consistent rational affine offset".to_owned(),
    ))?;
    graph.with_dim_arena_mut(|arena| unit.dimension.to_term(arena))?;
    if crate::dim::units::lookup(&unit.name).is_some_and(|known| known.to_value() == *unit) {
        Ok(graph.sym(&unit.name))
    } else {
        Ok(graph.lit(Payload::Unit(Arc::new(unit.clone()))))
    }
}

fn coherent_definition(dim: &crate::kernel::signature::DimExpr) -> UnitValue {
    UnitDef::linear("coherent", dim.clone(), crate::dim::Scale::one()).to_value()
}

fn same_dimension(
    left: &crate::kernel::signature::DimExpr,
    right: &crate::kernel::signature::DimExpr,
) -> bool {
    let left = crate::dim::units::resolve_value(&coherent_definition(left));
    let right = crate::dim::units::resolve_value(&coherent_definition(right));
    matches!((left, right), (Some(left), Some(right)) if left.dim == right.dim)
}

/// Dynamic dimensions describe magnitudes, not extra factors to multiply into
/// an already dimensioned quantity. Preserve existing units or reject conflicts.
fn term_of_dimensioned_value(
    graph: &mut Graph,
    value: &Value,
    dim: &crate::kernel::signature::DimExpr,
    session: SessionId,
) -> Result<NodeId, OpError> {
    let unit = coherent_definition(dim);
    crate::dim::units::resolve_value(&unit).ok_or_else(|| {
        OpError::Domain("a dynamic graph value needs a concrete rational dimension".to_owned())
    })?;
    match value {
        | Value::Dyn(inner) => {
            if !same_dimension(dim, &inner.dim) {
                return Err(OpError::Domain(
                    "conflicting nested dynamic dimensions".to_owned(),
                ));
            }
            return term_of_dimensioned_value(graph, &inner.value, dim, session);
        },
        | Value::Quantity(quantity) => {
            if !same_dimension(dim, &quantity.unit.dimension) {
                return Err(OpError::Domain(
                    "dynamic dimension conflicts with quantity unit".to_owned(),
                ));
            }
            return term_of_value(graph, value, session);
        },
        | Value::Tensor(tensor) => return term_of_tensor(graph, tensor, session, Some(dim)),
        | Value::Int(_) | Value::Real(_) | Value::Float(_) | Value::Symbol(_) | Value::Term(_) => {
        },
        | _ => {
            return Err(OpError::Domain(
                "only scalar magnitudes and numeric tensors can carry a graph dimension".to_owned(),
            ));
        },
    }
    let magnitude = term_of_value(graph, value, session)?;
    if matches!(value, Value::Symbol(_) | Value::Term(_)) {
        let report = graph.dimension(magnitude);
        if report
            .errors
            .iter()
            .any(crate::graph::dim::DimError::is_blocking)
        {
            return Err(OpError::Domain(
                "dynamic term has conflicting graph dimensions".to_owned(),
            ));
        }
        if report.dimension.unknown.is_empty()
            && !report.dimension.any
            && !report.dimension.non_metric
        {
            let known = graph.with_dim_arena_mut(|arena| {
                let term = crate::graph::dim::DimTerm::from_dimension(&report.dimension, arena)?;
                crate::kernel::signature::DimExpr::try_from_term(&term, arena)
            })?;
            if same_dimension(&known, dim) {
                return Ok(magnitude);
            }
            if !report.dimension.is_dimensionless() {
                return Err(OpError::Domain(
                    "dynamic dimension conflicts with term dimension".to_owned(),
                ));
            }
        }
    }
    if same_dimension(dim, &crate::kernel::signature::DimExpr::Dimless) {
        return Ok(magnitude);
    }
    let op = graph.ops().lookup("quantity").ok_or_else(|| {
        OpError::Unsupported {
            backend: "term".to_owned(),
            feature: "dimensioned values require the units rule set".to_owned(),
        }
    })?;
    let unit = term_of_unit(graph, &unit)?;
    Ok(graph.node(op, &[magnitude, unit]))
}

fn term_of_tensor(
    graph: &mut Graph,
    tensor: &TensorValue,
    session: SessionId,
    dimension: Option<&crate::kernel::signature::DimExpr>,
) -> Result<NodeId, OpError> {
    use crate::kernel::signature::DimExpr;
    use crate::kernel::signature::TypeExpr;
    fn build(
        graph: &mut Graph,
        shape: &[usize],
        data: &[Value],
        dimension: Option<&DimExpr>,
        session: SessionId,
    ) -> Result<NodeId, OpError> {
        let Some((&extent, tail)) = shape.split_first() else {
            let value = data
                .first()
                .ok_or_else(|| OpError::Domain("empty scalar tensor".to_owned()))?;
            return match dimension {
                | Some(dim) => term_of_dimensioned_value(graph, value, dim, session),
                | None => term_of_value(graph, value, session),
            };
        };
        let stride = data.len().checked_div(extent).unwrap_or(0);
        let items = (0..extent)
            .map(|index| {
                build(
                    graph,
                    tail,
                    &data[index * stride..(index + 1) * stride],
                    dimension,
                    session,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(graph.node(core::LIST, &items))
    }
    if !tensor.is_consistent() {
        return Err(OpError::Domain(
            "tensor data length differs from its shape".to_owned(),
        ));
    }
    let elem_dimension = match &tensor.elem {
        | TypeExpr::Dim(_, dim) => Some(dim),
        | _ => None,
    };
    if let (Some(left), Some(right)) = (dimension, elem_dimension)
        && !same_dimension(left, right)
    {
        return Err(OpError::Domain(
            "tensor element dimension conflicts with its dynamic dimension".to_owned(),
        ));
    }
    let dimension = dimension.or(elem_dimension);
    if let Some(zero) = tensor.shape.iter().position(|&n| n == 0)
        && (zero + 1 != tensor.shape.len()
            || dimension.is_some_and(|dim| !same_dimension(dim, &DimExpr::Dimless)))
    {
        return Err(OpError::Unsupported {
            backend: "term".to_owned(),
            feature: "empty tensor shape/dimension metadata not representable by nested lists"
                .to_owned(),
        });
    }
    build(graph, &tensor.shape, &tensor.data, dimension, session)
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
        return term_of_unit(graph, unit);
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
    use crate::kernel::result::Solutions;

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

    /// The caller's trace collector is carried through `clone_ctx`. It used to
    /// be dropped at every entry point, so `CallCtx::with_trace` recorded
    /// nothing and a host's `@trace` attribute had no data source.
    #[test]
    fn a_call_trace_reaches_the_caller() {
        let s = Session::new();
        let trace = std::cell::RefCell::new(crate::kernel::call::Trace::new());
        let ctx = CallCtx::new().with_trace(&trace);
        let _ = s
            .call_value("sin", &[Value::Float(0.5)], &ctx)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            !trace.borrow().is_empty(),
            "the caller's collector must receive the call events"
        );
    }

    /// `CallCtx::substitute` is applied to the argument terms before the call
    /// is built. It used to be recorded and never read, so a host that
    /// replaced a symbol got the unsubstituted answer with no diagnostic.
    ///
    /// The substitution happens when the arguments become terms, so overload
    /// resolution still sees the values the host passed; `sin(x)` resolves
    /// (it accepts a symbol) and `x` becomes `0` before the rewrite runs.
    #[test]
    fn a_call_context_substitution_is_applied() {
        let s = Session::new();
        let ctx = CallCtx::new().substitute("x", Value::int(0_i32));
        let out = s
            .call_value("sin", &[Value::Symbol("x".to_owned())], &ctx)
            .unwrap_or_else(|e| panic!("{e}"));
        // Without the substitution this is the term `sin(x)`.
        assert_eq!(out, Value::int(0_i32));
    }

    /// A handle whose node id names nothing in this graph is refused at the
    /// boundary instead of being indexed in the arena.
    #[test]
    fn a_stale_term_handle_is_refused() {
        let s = Session::new();
        let stale = TermRef {
            session: s.id(),
            node: u32::MAX,
        };
        let error = s
            .term_of(&Value::Term(stale))
            .err()
            .unwrap_or_else(|| panic!("a stale handle must not become a term"));
        assert!(
            matches!(error, OpError::Internal(_)),
            "expected Internal, got {error:?}"
        );
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

    /// An exact integral rational that crosses the ABI — `Number::Rat(0/1)`,
    /// which a host that builds a [`Value`] directly can send instead of the
    /// `Int(0)` that `Number::rat` produces — must enter the graph in normal
    /// form, or every structural literal rewrite misses it.
    ///
    /// This is what made `limit(sin(t)/t, exact 0)` answer `cos(0/1)`: the
    /// endpoint crossed as a rational zero, so `cos(0) => 1` (and `sin(0)`,
    /// `tan(0)`, `exp(0)`, `ln(1)`) did not fire and the exact contract failed.
    #[test]
    fn an_exact_rational_zero_is_normalised_at_the_abi_boundary() {
        use num_rational::BigRational;

        let s = Session::new();
        let exact = |numer: i64, denom: i64| {
            Value::Real(Number::Rat(BigRational::new(
                BigInt::from(numer),
                BigInt::from(denom),
            )))
        };
        // The literal becomes the integer zero, not the unrewritable `0/1`.
        let zero = s.term_of(&exact(0, 1)).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(zero.to_string(), "0");
        // The closed special functions fold at the exact point.
        let call = |name: &str, arg: Value| {
            s.call_value(name, &[arg], &CallCtx::new())
                .unwrap_or_else(|e| panic!("{name}: {e}"))
        };
        assert_eq!(call("cos", exact(0, 1)), Value::int(1_i32));
        assert_eq!(call("sin", exact(0, 1)), Value::int(0_i32));
        assert_eq!(call("tan", exact(0, 1)), Value::int(0_i32));
        assert_eq!(call("exp", exact(0, 1)), Value::int(1_i32));
        assert_eq!(call("ln", exact(1, 1)), Value::int(0_i32));
        // A non-integral exact rational keeps its exact form.
        let third = s.term_of(&exact(1, 3)).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(third.to_string(), "1/3");
        // The limit endpoint is exact, so the answer is the exact integer —
        // the closed form is not left as `cos(0)` (or `cos(0/1)`).
        let quotient = s.parse("sin(t)/t").unwrap_or_else(|e| panic!("{e}"));
        let answer = s
            .call_value(
                "limit",
                &[
                    Value::Term(s.handle(quotient)),
                    Value::Symbol("t".to_owned()),
                    exact(0, 1),
                ],
                &CallCtx::new(),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(answer, Value::int(1_i32));
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
                    Value::Unit(metre),
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

    /// A unit the *session* registered is still a quantity when a term is read
    /// back: `value_of` resolves the unit through the session's own registry.
    ///
    /// It used to resolve against the process-wide table only, so a session
    /// unit was unknown there and the value came back as an opaque
    /// [`Value::Term`]: the dimension the session itself declared was lost.
    #[test]
    fn a_session_unit_keeps_its_dimension_and_display_in_value_of() {
        use crate::dim::Scale;
        use crate::dim::UnitDef;
        use crate::kernel::signature::DimExpr;

        let mut s = Session::new();
        s.register_unit(
            UnitDef::linear("furlong", DimExpr::base("Length"), Scale::one()).display("fur"),
        )
        .unwrap_or_else(|e| panic!("{e}"));
        let term = s
            .parse("quantity(3, furlong)")
            .unwrap_or_else(|e| panic!("{e}"));
        let value = s.value_of(term).unwrap_or_else(|e| panic!("{e}"));
        let Value::Quantity(quantity) = value else {
            panic!("a session unit must stay a quantity, got {value:?}");
        };
        assert_eq!(quantity.unit.dimension, DimExpr::base("Length"));
        assert_eq!(quantity.unit.display, "fur");
        assert_eq!(quantity.value.as_f64(), Some(3.0));

        // The value still denotes the same term.
        let round_trip = s
            .term_of(&Value::Quantity(quantity))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(round_trip.to_string(), "quantity(3, furlong)");
    }

    /// A built-in unit keeps the printed form of its definition rather than the
    /// text a term happened to spell: `degC` prints as `°C`, and the quantity
    /// is still reported as affine.
    #[test]
    fn a_known_unit_keeps_its_display_in_value_of() {
        let s = Session::new();
        let term = s
            .parse("quantity(20, degC)")
            .unwrap_or_else(|e| panic!("{e}"));
        let value = s.value_of(term).unwrap_or_else(|e| panic!("{e}"));
        let Value::Quantity(quantity) = value else {
            panic!("expected a quantity, got {value:?}");
        };
        assert_eq!(quantity.unit.name, "degC");
        assert_eq!(quantity.unit.display, "\u{b0}C");
        assert!(quantity.unit.affine);
    }

    /// A repeated root is a multiple zero. The solver lists the distinct
    /// roots, so the multiplicity is the only place the repetition is visible;
    /// the exact square-free decomposition of the polynomial supplies it.
    #[test]
    fn a_repeated_polynomial_root_reports_its_multiplicity() {
        let s = Session::new();
        let x = s.sym("x");
        let equation = s
            .parse("(x - 1)^2*(x + 2) = 0")
            .unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(solutions.complete, "a solved polynomial is exhaustive");
        assert_eq!(solutions.variable.as_deref(), Some("x"));
        assert_eq!(solutions.roots.len(), 2, "{:?}", solutions.roots);
        // Ascending order: -2 is simple, 1 is the double root.
        assert_eq!(solutions.multiplicities, vec![1, 2]);
        assert_eq!(solutions.multiplicity(0), 1);
        assert_eq!(solutions.multiplicity(1), 2);
        assert_eq!(solutions.conditions, Vec::new());
        // The multiplicity crosses the ABI in the record, and the record
        // reads back to the same answer.
        let back = Solutions::from_value(&solutions.to_value()).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back.multiplicities, solutions.multiplicities);
    }

    /// A square-free polynomial has only simple roots, and saying so is a
    /// claim the coefficients support.
    #[test]
    fn a_square_free_polynomial_reports_simple_roots() {
        let s = Session::new();
        let x = s.sym("x");
        let equation = s
            .parse("x^2 - 3*x + 2 = 0")
            .unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(solutions.roots.len(), 2);
        assert_eq!(solutions.multiplicities, vec![1, 1]);
    }

    /// `x^2 = 2` is square-free, so each irrational root is simple: that needs
    /// no numeric matching of a root to a factor, only the exact decomposition.
    #[test]
    fn an_irrational_root_of_a_square_free_polynomial_is_simple() {
        let s = Session::new();
        let x = s.sym("x");
        let equation = s.parse("x^2 = 2").unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        // Both roots are expression handles; the list node is expanded rather
        // than handed over as one opaque root.
        assert!(solutions.complete, "the solver's list is exhaustive");
        assert_eq!(solutions.roots.len(), 2, "{:?}", solutions.roots);
        assert_eq!(solutions.multiplicities, vec![1, 1]);
        assert!(
            solutions.conditions.is_empty(),
            "{:?}",
            solutions.conditions
        );
    }

    /// An irrational root of a polynomial with a repeated factor is matched
    /// to its factor on values, with a margin; the multiplicity is the
    /// factor's, not the default.
    #[test]
    fn an_irrational_repeated_root_gets_its_factors_multiplicity() {
        let s = Session::new();
        let x = s.sym("x");
        let equation = s.parse("(x^2 - 2)^2 = 0").unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(solutions.complete);
        assert_eq!(solutions.roots.len(), 2, "{:?}", solutions.roots);
        assert_eq!(solutions.multiplicities, vec![2, 2]);
    }

    /// `x^2 - 2x + 1` has the double root 1.
    #[test]
    fn a_perfect_square_has_a_double_root() {
        let s = Session::new();
        let x = s.sym("x");
        let equation = s
            .parse("x^2 - 2*x + 1 = 0")
            .unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(solutions.complete);
        assert_eq!(solutions.roots.len(), 1, "{:?}", solutions.roots);
        assert_eq!(solutions.multiplicities, vec![2]);
    }

    /// `sin(x) = 0` has the roots `k pi`; listing `0` is a representative of
    /// the principal branch, and the answer must not claim to be complete.
    /// The same holds for every equation solved through the inverse of a
    /// periodic function, and for a periodic function of an affine argument.
    #[test]
    fn a_periodic_equation_is_not_complete() {
        use crate::kernel::result::ConditionKind;

        let s = Session::new();
        let x = s.sym("x");
        for text in [
            "sin(x) = 0",
            "cos(x) = c",
            "tan(x) = c",
            "sin(a*x + b) = c",
            "sin(x) = 1/2",
            "2*sin(x)^2 - sin(x) = 0",
            "sin(x^2) = 1/2",
        ] {
            let equation = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
            let solutions = s
                .solve_structured(equation, x)
                .unwrap_or_else(|e| panic!("{text}: {e}"));
            assert!(!solutions.complete, "{text}: {:?}", solutions.roots);
            assert!(!solutions.roots.is_empty(), "{text}: principal roots kept");
            assert!(
                solutions
                    .conditions
                    .iter()
                    .any(|condition| condition.kind == ConditionKind::Branch),
                "{text}: {:?}",
                solutions.conditions
            );
            // The incomplete flag and the branch condition cross the ABI.
            let back =
                Solutions::from_value(&solutions.to_value()).unwrap_or_else(|e| panic!("{e}"));
            assert!(!back.complete);
            assert_eq!(back.conditions, solutions.conditions);
        }
        // A non-periodic equation, and a periodic function of a parameter only,
        // stay complete.
        for text in ["x^2 = 4", "sin(a)*x = 1"] {
            let equation = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
            let solutions = s
                .solve_structured(equation, x)
                .unwrap_or_else(|e| panic!("{text}: {e}"));
            assert!(solutions.complete, "{text}");
        }
    }

    /// `x^2 = a` has real roots only for `a >= 0`, and `exp(x) = a` only for
    /// `a > 0`: the roots `+-sqrt(a)` and `ln(a)` say so by construction.
    /// `cos(x) = c` needs `c` in `[-1, 1]` for `acos(c)` to be real.
    #[test]
    fn realness_of_the_roots_is_a_condition() {
        use crate::kernel::result::ConditionKind;

        let s = Session::new();
        let x = s.sym("x");
        let solve = |text: &str| {
            let equation = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
            s.solve_structured(equation, x)
                .unwrap_or_else(|e| panic!("{text}: {e}"))
        };
        let has = |solutions: &Solutions, kind: ConditionKind, text: &str| {
            solutions
                .conditions
                .iter()
                .any(|c| c.kind == kind && c.text == text)
        };

        let even = solve("x^2 - a = 0");
        assert!(even.complete);
        assert_eq!(even.roots.len(), 2, "{:?}", even.roots);
        assert!(
            has(&even, ConditionKind::Nonnegative, "a >= 0"),
            "{:?}",
            even.conditions
        );
        // The roots are distinct except at a = 0: simple, with that said.
        assert_eq!(even.multiplicities, vec![1, 1]);
        assert!(
            has(&even, ConditionKind::Distinct, "a != 0"),
            "{:?}",
            even.conditions
        );

        let log = solve("exp(x) - a = 0");
        assert!(log.complete);
        assert_eq!(log.roots.len(), 1, "{:?}", log.roots);
        assert!(
            has(&log, ConditionKind::Positive, "a > 0"),
            "{:?}",
            log.conditions
        );

        let inverse = solve("cos(x) = c");
        assert!(
            has(&inverse, ConditionKind::Interval, "-1 <= c <= 1"),
            "{:?}",
            inverse.conditions
        );

        let quadratic = solve("a*x^2 + b*x + c = 0");
        assert!(quadratic.complete);
        assert!(
            quadratic
                .conditions
                .iter()
                .any(|c| c.kind == ConditionKind::Nonnegative && c.text.ends_with(">= 0")),
            "{:?}",
            quadratic.conditions
        );

        // Kinds serialise by name, the new one included.
        let back = Solutions::from_value(&even.to_value()).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back.conditions, even.conditions);
    }

    /// `solve(a*x + b = 0, x)` is `-b/a`: the answer is valid only where
    /// `a != 0`, and that is the one condition the root states by being a
    /// quotient. It is reported; branch conditions the solver did not decide
    /// are not.
    #[test]
    fn a_parameter_denominator_in_a_root_is_a_condition() {
        use crate::kernel::result::ConditionKind;

        let s = Session::new();
        let x = s.sym("x");
        let equation = s.parse("a*x + b = 0").unwrap_or_else(|e| panic!("{e}"));
        let solutions = s
            .solve_structured(equation, x)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(solutions.roots.len(), 1);
        assert!(
            solutions.conditions.iter().any(|condition| {
                condition.kind == ConditionKind::NonZero && condition.text == "a != 0"
            }),
            "expected `a != 0`, got {:?}",
            solutions.conditions
        );
        // A linear polynomial has one simple root wherever it exists, which
        // the condition `a != 0` already states.
        assert_eq!(solutions.multiplicities, vec![1]);
        // The condition crosses the ABI in the record.
        let back = Solutions::from_value(&solutions.to_value()).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back.conditions, solutions.conditions);
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
