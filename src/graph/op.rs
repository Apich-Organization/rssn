//! The open operator registry.
//!
//! Operators are *data*, not enum variants: a domain rule set, a plugin, or a
//! future JIT backend registers an [`OpDescriptor`] and receives an
//! [`OpId`]. The kernel itself only knows the handful of structural
//! operators in [`core`].
//!
//! # What a descriptor carries
//!
//! * what the scheduler and the matcher need — arity, flags, cost, binder;
//! * what evaluation needs — a scalar `fn(&[f64]) -> f64` fast path, or a
//!   [`EvalDyn`] that sees [`Value`]s and can capture state;
//! * what a *host* needs — the rule set that registered it, typed
//!   [`Signature`]s, aliases, category, documentation, examples, stability,
//!   purity, determinism and thread safety.
//!
//! The last group exists so that a compiler front end can generate bindings
//! from the [manifest](crate::kernel::manifest) instead of from a
//! hand-maintained list. Nothing in the engine reads it.
//!
//! # Overloads
//!
//! A name may carry several operators when they differ in arity or in
//! signature. That is opt-in: a second registration of the same name is a
//! conflict unless the new descriptor declares itself an overload with
//! [`OpDescriptor::overload`]. Idempotent re-registration (the same
//! descriptor twice, which happens whenever two rule sets both need an
//! operator) always returns the existing id.

use std::any::Any;
use std::any::TypeId;
use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

use super::id::OpId;
use super::id::to_u32;
use crate::backend::Lowering;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::BinderInfo;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DynamicReason;
use crate::kernel::signature::Example;
use crate::kernel::signature::OpFlagsInfo;
use crate::kernel::signature::Purity;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::value::Value;

/// How many children an operator takes.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Arity {
    /// Exactly this many.
    Fixed(u8),
    /// Any number, including zero.
    Variadic,
}

impl Arity {
    /// Whether `n` children are acceptable.
    #[must_use]
    pub fn accepts(
        self,
        n: usize,
    ) -> bool {
        match self {
            | Self::Fixed(k) => usize::from(k) == n,
            | Self::Variadic => true,
        }
    }

    /// The least acceptable child count.
    #[must_use]
    pub const fn min(self) -> u8 {
        match self {
            | Self::Fixed(k) => k,
            | Self::Variadic => 0,
        }
    }

    /// The greatest acceptable child count, if any.
    #[must_use]
    pub const fn max(self) -> Option<u8> {
        match self {
            | Self::Fixed(k) => Some(k),
            | Self::Variadic => None,
        }
    }
}

/// Algebraic and scheduling properties of an operator.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct OpFlags(u16);

impl OpFlags {
    /// Nested applications may be flattened: `f(f(a, b), c) = f(a, b, c)`.
    pub const ASSOCIATIVE: Self = Self(1 << 1);
    /// Children may be reordered freely.
    pub const COMMUTATIVE: Self = Self(1);
    /// A *heavy* operator: an unevaluated identity-transformation request
    /// (derivative, integral, solve, ...) that the scheduler tries to reduce
    /// before anything else and that closed-form extraction refuses.
    pub const HEAVY: Self = Self(1 << 2);
    /// The node is a leaf carrying a payload.
    pub const LEAF: Self = Self(1 << 3);
    /// No flags.
    pub const NONE: Self = Self(0);
    /// A heavy operator that is irreducible — and so a closed form — when
    /// its first argument is an undetermined function `f(x, ...)` or
    /// another such application: `diff(y(x), x)` cannot be computed any
    /// further, it *is* the answer.
    pub const OPAQUE_ON_APPLY: Self = Self(1 << 6);
    /// The operator's value is a truth value. Its scalar semantics encode
    /// true and false as 1 and 0 for evaluation, but numeric folding must
    /// not replace the term by that number.
    pub const PREDICATE: Self = Self(1 << 5);
    /// Tree windows always look inside this operator, even when the node
    /// is shared. Local normalisation of arithmetic needs to see the
    /// structure of `x * x^(-1)` whether or not `x^(-1)` is used elsewhere.
    pub const TRANSPARENT: Self = Self(1 << 4);

    /// Union of two flag sets.
    #[must_use]
    pub const fn with(
        self,
        other: Self,
    ) -> Self {
        Self(self.0 | other.0)
    }

    /// Whether every flag in `other` is set.
    #[must_use]
    pub const fn has(
        self,
        other: Self,
    ) -> bool {
        self.0 & other.0 == other.0
    }

    /// The flags as manifest data.
    #[must_use]
    pub const fn info(self) -> OpFlagsInfo {
        OpFlagsInfo {
            associative: self.has(Self::ASSOCIATIVE),
            commutative: self.has(Self::COMMUTATIVE),
            heavy: self.has(Self::HEAVY),
            leaf: self.has(Self::LEAF),
            predicate: self.has(Self::PREDICATE),
            transparent: self.has(Self::TRANSPARENT),
            opaque_on_apply: self.has(Self::OPAQUE_ON_APPLY),
        }
    }
}

/// Describes which child of a binding operator is the bound variable and in
/// which children it is bound.
///
/// `integral(body, x, lo, hi)` binds child 1 inside child 0 only, so `x` is
/// not free in the integral even though it is free in `body`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Binder {
    /// Index of the child holding the bound symbol.
    pub var: u8,
    /// Bit `i` set means the variable is bound inside child `i`.
    pub scope: u32,
}

impl Binder {
    /// The binder as manifest data.
    #[must_use]
    pub const fn info(self) -> BinderInfo {
        BinderInfo {
            var: self.var,
            scope: self.scope,
        }
    }
}

/// Scalar numeric semantics of an operator: maps child values to a value.
///
/// A plain function pointer on purpose — it is what the interpreter calls
/// and what a JIT backend can emit a direct call to.
pub type EvalFn = fn(&[f64]) -> f64;

/// Dynamic numeric semantics: sees [`Value`]s, may capture state, may fail.
///
/// This is the extension point that a `fn` pointer cannot provide. An
/// operator registered with an `EvalDyn` works in the dynamic ABI
/// ([`Session::call_value`](crate::api::Session::call_value)) whether or not
/// it also has a scalar fast path.
pub trait EvalDyn: Send + Sync {
    /// Evaluates the operator on `args`.
    ///
    /// # Errors
    /// Returns a structured [`OpError`]; a domain violation is
    /// [`OpError::Domain`], a type mismatch [`OpError::Type`].
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError>;

    /// The signature of the value this evaluator produces.
    fn signature(&self) -> Signature;
}

impl fmt::Debug for dyn EvalDyn {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str("EvalDyn")
    }
}

/// Everything the kernel and the backends need to know about an operator.
#[derive(Clone)]
pub struct OpDescriptor {
    /// Unique name, also the function-call spelling in the pattern syntax.
    pub name: Arc<str>,
    /// Accepted child count.
    pub arity: Arity,
    /// Algebraic properties.
    pub flags: OpFlags,
    /// Base cost for extraction (children are added on top).
    pub cost: u32,
    /// Variable binding structure, if any.
    pub binder: Option<Binder>,
    /// Scalar numeric semantics, if the operator has any.
    pub eval: Option<EvalFn>,
    /// Dynamic semantics, if the operator has any.
    pub eval_dyn: Option<Arc<dyn EvalDyn>>,
    /// Typed overloads. Empty means "no signature yet"; `dynamic_reason`
    /// then says why.
    pub overloads: Vec<Signature>,
    /// Markdown documentation.
    pub doc: Arc<str>,
    /// Domain category (`calculus`, `linear-algebra`, ...).
    pub category: Arc<str>,
    /// Alternative spellings accepted by the dynamic ABI.
    pub aliases: Vec<Arc<str>>,
    /// How settled the operator is.
    pub stability: Stability,
    /// Version the operator appeared in.
    pub since: Arc<str>,
    /// Whether it touches anything but its arguments.
    pub purity: Purity,
    /// Whether repeated calls agree.
    pub determinism: Determinism,
    /// How it may be used across threads.
    pub thread_safety: ThreadSafety,
    /// Why there is no signature, when there is none.
    pub dynamic_reason: Option<DynamicReason>,
    /// Worked examples.
    pub examples: Vec<Example>,
    /// Backends that can execute it.
    pub backend: BackendSupport,
    /// Rule set that registered it; filled in by the installer.
    pub rule_set: Arc<str>,
    /// Set when this registration deliberately adds an overload to an
    /// existing name instead of clashing with it.
    pub overload: bool,
    /// Whether the operator takes a function argument.
    pub higher_order: bool,
    /// How the JIT should compile it, when it can.
    ///
    /// `None` means "let the backend decide": a scalar `eval` is called
    /// through a shim, anything else runs in the interpreter. Setting it
    /// explicitly is how a plugin gets native code for its own operator
    /// without touching the backend.
    pub jit: Option<Lowering>,
}

impl fmt::Debug for OpDescriptor {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("OpDescriptor")
            .field("name", &self.name)
            .field("arity", &self.arity)
            .field("flags", &self.flags)
            .field("cost", &self.cost)
            .field("binder", &self.binder)
            .field("eval", &self.eval.map(|_| "fn"))
            .field("eval_dyn", &self.eval_dyn.is_some())
            .field("overloads", &self.overloads.len())
            .field("stability", &self.stability)
            .field("jit", &self.jit)
            .finish_non_exhaustive()
    }
}

impl OpDescriptor {
    /// A plain operator with no special properties.
    ///
    /// Its cost is 3: a named function weighs more than an arithmetic
    /// operator (cost 1) with a literal operand (cost 1), so `-x` is
    /// preferred to `abs(x)` when both are known to be equal. A nullary
    /// operator (a named constant) costs 1, like a symbol.
    #[must_use]
    pub fn new(
        name: &str,
        arity: Arity,
    ) -> Self {
        Self {
            name: Arc::from(name),
            arity,
            flags: OpFlags::NONE,
            // A named constant (`pi`, `hbar`) is as simple as a symbol.
            cost: if matches!(arity, Arity::Fixed(0)) {
                1
            } else {
                3
            },
            binder: None,
            eval: None,
            eval_dyn: None,
            overloads: Vec::new(),
            doc: Arc::from(""),
            category: Arc::from(""),
            aliases: Vec::new(),
            stability: Stability::Stable,
            since: Arc::from(""),
            purity: Purity::Pure,
            determinism: Determinism::Deterministic,
            thread_safety: ThreadSafety::SessionBound,
            dynamic_reason: None,
            examples: Vec::new(),
            backend: BackendSupport::default(),
            rule_set: Arc::from(""),
            overload: false,
            higher_order: false,
            jit: None,
        }
    }

    /// Adds flags.
    #[must_use]
    pub const fn flags(
        mut self,
        flags: OpFlags,
    ) -> Self {
        self.flags = self.flags.with(flags);
        self
    }

    /// Sets the extraction cost.
    #[must_use]
    pub const fn cost(
        mut self,
        cost: u32,
    ) -> Self {
        self.cost = cost;
        self
    }

    /// Declares the binding structure.
    #[must_use]
    pub const fn binder(
        mut self,
        var: u8,
        scope: u32,
    ) -> Self {
        self.binder = Some(Binder { var, scope });
        self
    }

    /// Attaches scalar numeric semantics.
    ///
    /// A scalar evaluator is a pure function of its arguments, so the
    /// operator becomes callable from any thread.
    #[must_use]
    pub fn eval(
        mut self,
        eval: EvalFn,
    ) -> Self {
        self.eval = Some(eval);
        if self.thread_safety == ThreadSafety::SessionBound {
            self.thread_safety = ThreadSafety::Send;
        }
        self
    }

    /// Attaches dynamic semantics.
    #[must_use]
    pub fn eval_dyn(
        mut self,
        eval: Arc<dyn EvalDyn>,
    ) -> Self {
        self.eval_dyn = Some(eval);
        self
    }

    /// Adds a typed overload.
    #[must_use]
    pub fn signature(
        mut self,
        signature: Signature,
    ) -> Self {
        self.overloads.push(signature);
        self
    }

    /// Replaces the overload list.
    #[must_use]
    pub fn signatures(
        mut self,
        signatures: Vec<Signature>,
    ) -> Self {
        self.overloads = signatures;
        self
    }

    /// Sets the documentation.
    #[must_use]
    pub fn doc(
        mut self,
        doc: &str,
    ) -> Self {
        self.doc = Arc::from(doc);
        self
    }

    /// Sets the category.
    #[must_use]
    pub fn category(
        mut self,
        category: &str,
    ) -> Self {
        self.category = Arc::from(category);
        self
    }

    /// Adds an accepted alternative spelling.
    #[must_use]
    pub fn alias(
        mut self,
        alias: &str,
    ) -> Self {
        self.aliases.push(Arc::from(alias));
        self
    }

    /// Sets the stability level.
    #[must_use]
    pub const fn stability(
        mut self,
        stability: Stability,
    ) -> Self {
        self.stability = stability;
        self
    }

    /// Sets the version the operator appeared in.
    #[must_use]
    pub fn since(
        mut self,
        since: &str,
    ) -> Self {
        self.since = Arc::from(since);
        self
    }

    /// Sets purity.
    #[must_use]
    pub const fn purity(
        mut self,
        purity: Purity,
    ) -> Self {
        self.purity = purity;
        self
    }

    /// Sets determinism.
    #[must_use]
    pub const fn determinism(
        mut self,
        determinism: Determinism,
    ) -> Self {
        self.determinism = determinism;
        self
    }

    /// Sets thread safety.
    #[must_use]
    pub const fn thread_safety(
        mut self,
        thread_safety: ThreadSafety,
    ) -> Self {
        self.thread_safety = thread_safety;
        self
    }

    /// Records why the operator has no signature.
    #[must_use]
    pub const fn dynamic_reason(
        mut self,
        reason: DynamicReason,
    ) -> Self {
        self.dynamic_reason = Some(reason);
        self
    }

    /// Adds a worked example.
    #[must_use]
    pub fn example(
        mut self,
        input: &str,
        output: &str,
    ) -> Self {
        self.examples.push(Example {
            input: input.to_owned(),
            output: output.to_owned(),
            note: String::new(),
        });
        self
    }

    /// Declares which backends can run the operator.
    #[must_use]
    pub const fn backend(
        mut self,
        backend: BackendSupport,
    ) -> Self {
        self.backend = backend;
        self
    }

    /// Marks the operator as taking a function argument.
    #[must_use]
    pub const fn higher_order(mut self) -> Self {
        self.higher_order = true;
        self
    }

    /// Tells the JIT how to compile this operator.
    ///
    /// This is the only place an operator talks to a backend, and it does so
    /// in the backend's own vocabulary rather than by naming Cranelift: a
    /// rule set that never enables the `jit` feature still compiles.
    #[must_use]
    pub const fn jit(
        mut self,
        lowering: Lowering,
    ) -> Self {
        self.jit = Some(lowering);
        self
    }

    /// Declares that this registration adds an overload to an existing name
    /// rather than clashing with it.
    #[must_use]
    pub const fn overload(mut self) -> Self {
        self.overload = true;
        self
    }

    /// The same descriptor, attributed to `rule_set` if it has no owner yet.
    #[must_use]
    pub(crate) fn owned_by(
        mut self,
        rule_set: &Arc<str>,
    ) -> Self {
        if self.rule_set.is_empty() {
            self.rule_set = Arc::clone(rule_set);
        }
        self
    }

    /// Whether `self` describes the same operator as `other`: same name,
    /// arity, flags and binder. Metadata differences do not matter, so two
    /// rule sets may register `sin` and the second one simply agrees.
    #[must_use]
    fn same_operator(
        &self,
        other: &Self,
    ) -> bool {
        self.name == other.name
            && self.arity == other.arity
            && self.flags == other.flags
            && self.binder == other.binder
    }

    /// Copies metadata the existing registration lacks from `self`.
    ///
    /// A rule set that registers an operator another set already introduced
    /// does not lose the documentation and signature the first one supplied,
    /// but it can fill a gap.
    fn fill_gaps_from(
        &mut self,
        other: &Self,
    ) {
        if self.doc.is_empty() {
            self.doc = Arc::clone(&other.doc);
        }
        if self.category.is_empty() {
            self.category = Arc::clone(&other.category);
        }
        if self.overloads.is_empty() {
            self.overloads.clone_from(&other.overloads);
        }
        if self.dynamic_reason.is_none() {
            self.dynamic_reason = other.dynamic_reason;
        }
        if self.eval_dyn.is_none() {
            self.eval_dyn.clone_from(&other.eval_dyn);
        }
        if self.since.is_empty() {
            self.since = Arc::clone(&other.since);
        }
        if self.jit.is_none() {
            self.jit = other.jit;
        }
    }
}

/// Error returned when an operator registration conflicts with an existing
/// one of the same name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpConflict {
    /// The offending name.
    pub name: String,
    /// What the registered operator takes.
    pub existing: Arity,
    /// What the rejected operator takes.
    pub attempted: Arity,
}

impl fmt::Display for OpConflict {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "operator `{}` is already registered with a different signature \
             (existing arity {:?}, attempted {:?})",
            self.name, self.existing, self.attempted
        )
    }
}

impl std::error::Error for OpConflict {}

/// Structural operators every graph starts with.
pub mod core {
    use super::OpId;

    /// Literal leaf: number, boolean, string or blob payload.
    pub const LIT: OpId = OpId(0);
    /// Symbol leaf.
    pub const SYM: OpId = OpId(1);
    /// N-ary commutative, associative sum.
    pub const ADD: OpId = OpId(2);
    /// N-ary commutative, associative product.
    pub const MUL: OpId = OpId(3);
    /// Binary power `base ^ exponent`.
    pub const POW: OpId = OpId(4);
    /// Application of an undetermined function: `apply(f, args...)`.
    pub const APPLY: OpId = OpId(5);
    /// Ordered tuple of children.
    pub const LIST: OpId = OpId(6);
    /// The equation `lhs = rhs`, as a term.
    pub const EQ: OpId = OpId(7);
}

/// Registry mapping [`OpId`]s to descriptors and names to ids.
///
/// Besides the descriptor, any number of typed *attributes* can be attached
/// to an operator: a derivative rule, a LaTeX template, an interval
/// enclosure, a JIT lowering. The kernel never looks at them; the layer
/// that defines an attribute type is the one that reads it. This is how new
/// capabilities are added without touching the kernel or the operators'
/// definitions.
#[derive(Clone)]
pub struct OpTable {
    ops: Vec<OpDescriptor>,
    by_name: HashMap<Arc<str>, Vec<OpId>>,
    by_alias: HashMap<Arc<str>, OpId>,
    attrs: HashMap<(OpId, TypeId), Arc<dyn Any + Send + Sync>>,
}

impl fmt::Debug for OpTable {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "OpTable({} operators, {} attributes)",
            self.ops.len(),
            self.attrs.len()
        )
    }
}

impl Default for OpTable {
    fn default() -> Self {
        Self::new()
    }
}

impl OpTable {
    /// Creates a table holding exactly the [`core`] operators.
    #[must_use]
    pub fn new() -> Self {
        let mut table = Self {
            ops: Vec::new(),
            by_name: HashMap::new(),
            by_alias: HashMap::new(),
            attrs: HashMap::new(),
        };
        let ac = OpFlags::COMMUTATIVE.with(OpFlags::ASSOCIATIVE);
        let builtins = [
            OpDescriptor::new("lit", Arity::Fixed(0))
                .flags(OpFlags::LEAF)
                .cost(1)
                .doc("A literal leaf: a number, a truth value, a string or a blob.")
                .category("kernel")
                .since("0.1")
                .signature(Signature::new(
                    vec![crate::kernel::signature::TypeExpr::Opaque(
                        "payload".to_owned(),
                    )],
                    crate::kernel::signature::TypeExpr::Real,
                )),
            // A symbol costs more than a number and an operator together,
            // so that `2*x` beats `x + x` and `x^2` beats `x*x`.
            OpDescriptor::new("sym", Arity::Fixed(0))
                .flags(OpFlags::LEAF)
                .cost(3)
                .doc("A symbol leaf.")
                .category("kernel")
                .since("0.1")
                .signature(Signature::new(
                    Vec::new(),
                    crate::kernel::signature::TypeExpr::Var("T".to_owned()),
                )),
            OpDescriptor::new("add", Arity::Variadic)
                .flags(ac.with(OpFlags::TRANSPARENT))
                .cost(1)
                .eval(|a| a.iter().sum())
                .doc("N-ary sum. Associative and commutative.")
                .category("kernel")
                .since("0.1")
                .signature(
                    Signature::scalars(0, crate::kernel::signature::TypeExpr::Real).variadic(
                        crate::kernel::signature::Variadic {
                            min: 0,
                            max: None,
                            elem: crate::kernel::signature::TypeExpr::Real,
                        },
                    ),
                ),
            OpDescriptor::new("mul", Arity::Variadic)
                .flags(ac.with(OpFlags::TRANSPARENT))
                .cost(1)
                .eval(|a| a.iter().product())
                .doc("N-ary product. Associative and commutative.")
                .category("kernel")
                .since("0.1")
                .signature(
                    Signature::scalars(0, crate::kernel::signature::TypeExpr::Real).variadic(
                        crate::kernel::signature::Variadic {
                            min: 0,
                            max: None,
                            elem: crate::kernel::signature::TypeExpr::Real,
                        },
                    ),
                ),
            OpDescriptor::new("pow", Arity::Fixed(2))
                .flags(OpFlags::TRANSPARENT)
                .cost(1)
                .eval(|a| {
                    match a {
                        | [b, e] => b.powf(*e),
                        | _ => f64::NAN,
                    }
                })
                .doc("`base ^ exponent`.")
                .category("kernel")
                .since("0.1")
                .signature(Signature::scalars(2, crate::kernel::signature::TypeExpr::Real)),
            OpDescriptor::new("apply", Arity::Variadic)
                .cost(2)
                .doc("Application of an undetermined function: `apply(f, args...)`.")
                .category("kernel")
                .since("0.1")
                .dynamic_reason(DynamicReason::FullyPolymorphic),
            OpDescriptor::new("list", Arity::Variadic)
                .cost(1)
                .doc("An ordered tuple of children.")
                .category("kernel")
                .since("0.1")
                .signature(
                    Signature::new(
                        Vec::new(),
                        crate::kernel::signature::TypeExpr::Tensor {
                            rank: crate::kernel::signature::NatExpr::Var("n".to_owned()),
                            shape: None,
                            elem: Box::new(crate::kernel::infer::term()),
                        },
                    )
                    // A list takes any number of elements, of any type: the
                    // element type of the resulting tensor depends on them.
                    .variadic(crate::kernel::signature::Variadic {
                        min: 0,
                        max: None,
                        elem: crate::kernel::infer::term(),
                    }),
                ),
            OpDescriptor::new("eq", Arity::Fixed(2))
                .cost(1)
                .doc("The equation `lhs = rhs`, as a term.")
                .category("kernel")
                .since("0.1")
                .signature(Signature::scalars(2, crate::kernel::signature::TypeExpr::Bool)),
        ];
        for desc in builtins {
            table.push(desc.owned_by(&Arc::from("kernel")));
        }
        // Core dimensional signatures. They are attributes, like `Partials`,
        // so they can be overridden by a rule set after registration.
        for name in ["add", "mul", "pow", "list", "eq", "apply"] {
            if let Some(id) = table.lookup(name)
                && let Some(rule) = crate::graph::dim::default_rule(name)
            {
                table.set_attr(id, rule);
            }
        }
        table
    }

    fn push(
        &mut self,
        desc: OpDescriptor,
    ) -> OpId {
        let id = OpId(to_u32(self.ops.len()));
        self.by_name
            .entry(Arc::clone(&desc.name))
            .or_default()
            .push(id);
        for alias in &desc.aliases {
            self.by_alias.entry(Arc::clone(alias)).or_insert(id);
        }
        self.ops.push(desc);
        id
    }

    /// Registers an operator, or returns the existing id when an operator
    /// with the same name, arity, flags and binder is already present.
    ///
    /// # Errors
    /// Returns [`OpConflict`] when the name is taken by an operator with a
    /// different shape and the new descriptor does not declare itself an
    /// overload.
    pub fn register(
        &mut self,
        desc: OpDescriptor,
    ) -> Result<OpId, OpConflict> {
        let Some(existing) = self.by_name.get(&desc.name).cloned() else {
            return Ok(self.push(desc));
        };
        for id in &existing {
            if self.ops.get(id.index()).is_some_and(|old| old.same_operator(&desc)) {
                if let Some(slot) = self.ops.get_mut(id.index()) {
                    slot.fill_gaps_from(&desc);
                }
                return Ok(*id);
            }
        }
        if desc.overload {
            return Ok(self.push(desc));
        }
        let first = existing
            .first()
            .and_then(|id| self.ops.get(id.index()))
            .map_or(Arity::Variadic, |d| d.arity);
        Err(OpConflict {
            name: desc.name.to_string(),
            existing: first,
            attempted: desc.arity,
        })
    }

    /// Looks an operator up by name or alias, returning the primary
    /// registration.
    #[must_use]
    pub fn lookup(
        &self,
        name: &str,
    ) -> Option<OpId> {
        self.by_name
            .get(name)
            .and_then(|ids| ids.first().copied())
            .or_else(|| self.by_alias.get(name).copied())
    }

    /// Every operator registered under `name`, primary first.
    #[must_use]
    pub fn lookup_all(
        &self,
        name: &str,
    ) -> &[OpId] {
        self.by_name.get(name).map_or(&[], Vec::as_slice)
    }

    /// The registration of `name` that accepts `arity` children.
    ///
    /// Used by the dynamic ABI to resolve an overload from the call's
    /// argument count.
    #[must_use]
    pub fn lookup_arity(
        &self,
        name: &str,
        arity: usize,
    ) -> Option<OpId> {
        self.lookup_all(name)
            .iter()
            .copied()
            .find(|id| self.get(*id).arity.accepts(arity))
            .or_else(|| {
                self.by_alias
                    .get(name)
                    .copied()
                    .filter(|id| self.get(*id).arity.accepts(arity))
            })
    }

    /// The registration of `name` whose descriptor matches `desc` in name,
    /// arity, flags and binder — how a graph import finds its counterpart.
    #[must_use]
    pub fn lookup_exact(
        &self,
        desc: &OpDescriptor,
    ) -> Option<OpId> {
        self.by_name
            .get(&desc.name)
            .into_iter()
            .flatten()
            .copied()
            .find(|id| self.ops.get(id.index()).is_some_and(|old| old.same_operator(desc)))
    }

    /// Iterates over the registered operators and their descriptors.
    pub fn iter(&self) -> impl Iterator<Item = (OpId, &OpDescriptor)> + '_ {
        self.ops
            .iter()
            .enumerate()
            .map(|(index, descriptor)| (OpId(u32::try_from(index).unwrap_or(u32::MAX)), descriptor))
    }

    /// The descriptor of `op`.
    ///
    /// # Panics
    /// Panics when `op` was not produced by this table; that is a logic
    /// error in the caller, never a data-dependent condition.
    #[must_use]
    pub fn get(
        &self,
        op: OpId,
    ) -> &OpDescriptor {
        match self.ops.get(op.index()) {
            | Some(d) => d,
            | None => foreign_op(op),
        }
    }

    /// The descriptor of `op`, if it belongs to this table.
    #[must_use]
    pub fn try_get(
        &self,
        op: OpId,
    ) -> Option<&OpDescriptor> {
        self.ops.get(op.index())
    }

    /// Attaches (or replaces) the attribute of type `T` on `op`.
    pub fn set_attr<T: Any + Send + Sync>(
        &mut self,
        op: OpId,
        value: T,
    ) {
        self.attrs.insert((op, TypeId::of::<T>()), Arc::new(value));
    }

    /// The attribute of type `T` attached to `op`, if any.
    #[must_use]
    pub fn attr<T: Any + Send + Sync>(
        &self,
        op: OpId,
    ) -> Option<&T> {
        self.attrs
            .get(&(op, TypeId::of::<T>()))
            .and_then(|a| a.downcast_ref::<T>())
    }

    /// Number of registered operators.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.ops.len()
    }

    /// Always `false`: the core operators are always present.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }
}

#[cold]
fn foreign_op(op: OpId) -> ! {
    panic!("rssn graph kernel: operator {op:?} does not belong to this table")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_ids_match_registration_order() {
        let t = OpTable::new();
        assert_eq!(t.lookup("lit"), Some(core::LIT));
        assert_eq!(t.lookup("sym"), Some(core::SYM));
        assert_eq!(t.lookup("add"), Some(core::ADD));
        assert_eq!(t.lookup("mul"), Some(core::MUL));
        assert_eq!(t.lookup("pow"), Some(core::POW));
        assert_eq!(t.lookup("apply"), Some(core::APPLY));
        assert_eq!(t.lookup("list"), Some(core::LIST));
        assert_eq!(t.lookup("eq"), Some(core::EQ));
    }

    #[test]
    fn registration_is_idempotent_and_detects_conflicts() {
        let mut t = OpTable::new();
        let sin = t.register(OpDescriptor::new("sin", Arity::Fixed(1)));
        assert_eq!(t.register(OpDescriptor::new("sin", Arity::Fixed(1))), sin);
        assert!(
            t.register(OpDescriptor::new("sin", Arity::Fixed(2)))
                .is_err()
        );
    }

    #[test]
    fn overloads_share_a_name_but_not_an_id() {
        let mut t = OpTable::new();
        let unary = t
            .register(OpDescriptor::new("norm", Arity::Fixed(1)))
            .unwrap_or_else(|e| panic!("{e}"));
        let binary = t
            .register(OpDescriptor::new("norm", Arity::Fixed(2)).overload())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_ne!(unary, binary);
        assert_eq!(t.lookup("norm"), Some(unary), "the primary registration");
        assert_eq!(t.lookup_all("norm").len(), 2);
        assert_eq!(t.lookup_arity("norm", 2), Some(binary));
        assert_eq!(t.lookup_arity("norm", 3), None);
    }

    #[test]
    fn aliases_resolve_to_their_operator() {
        let mut t = OpTable::new();
        let id = t
            .register(OpDescriptor::new("ln", Arity::Fixed(1)).alias("log"))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(t.lookup("log"), Some(id));
        assert_eq!(t.lookup_all("log").len(), 0);
    }

    #[test]
    fn metadata_is_merged_on_re_registration() {
        let mut t = OpTable::new();
        let _ = t.register(
            OpDescriptor::new("sin", Arity::Fixed(1))
                .doc("Sine.")
                .category("elementary"),
        );
        let id = t
            .register(OpDescriptor::new("sin", Arity::Fixed(1)))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(&*t.get(id).doc, "Sine.");
        assert_eq!(&*t.get(id).category, "elementary");
    }

    #[test]
    fn attributes_are_typed_per_operator() {
        struct Latex(&'static str);
        struct Priority(u8);
        let mut t = OpTable::new();
        t.set_attr(core::ADD, Latex("+"));
        t.set_attr(core::ADD, Priority(1));
        t.set_attr(core::MUL, Latex("\\cdot"));
        assert_eq!(t.attr::<Latex>(core::ADD).map(|l| l.0), Some("+"));
        assert_eq!(t.attr::<Priority>(core::ADD).map(|p| p.0), Some(1));
        assert_eq!(t.attr::<Latex>(core::MUL).map(|l| l.0), Some("\\cdot"));
        assert!(t.attr::<Priority>(core::MUL).is_none());
        t.set_attr(core::ADD, Latex("plus"));
        assert_eq!(t.attr::<Latex>(core::ADD).map(|l| l.0), Some("plus"));
    }

    #[test]
    fn core_eval() {
        let t = OpTable::new();
        let add = t.get(core::ADD).eval.map(|f| f(&[1.0, 2.0, 3.0]));
        let pow = t.get(core::POW).eval.map(|f| f(&[2.0, 10.0]));
        assert_eq!(add, Some(6.0));
        assert_eq!(pow, Some(1024.0));
    }

    #[test]
    fn every_core_operator_has_a_signature_or_a_reason() {
        let t = OpTable::new();
        for (_, desc) in t.iter() {
            assert!(
                !desc.overloads.is_empty() || desc.dynamic_reason.is_some(),
                "`{}` has neither a signature nor a dynamic reason",
                desc.name
            );
            assert_eq!(&*desc.rule_set, "kernel");
        }
    }
}
