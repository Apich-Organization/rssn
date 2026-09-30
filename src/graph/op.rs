//! The open operator registry.
//!
//! Operators are *data*, not enum variants: a domain rule set, a plugin, or
//! a future JIT backend registers an [`OpDescriptor`] and receives an
//! [`OpId`]. The kernel itself only knows the handful of structural
//! operators in [`core`].

use std::any::Any;
use std::any::TypeId;
use std::collections::HashMap;
use std::sync::Arc;

use super::id::OpId;
use super::id::to_u32;

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
}

/// Algebraic and scheduling properties of an operator.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct OpFlags(u16);

impl OpFlags {
    /// Children may be reordered freely.
    pub const COMMUTATIVE: Self = Self(1);
    /// Nested applications may be flattened: `f(f(a, b), c) = f(a, b, c)`.
    pub const ASSOCIATIVE: Self = Self(1 << 1);
    /// A *heavy* operator: an unevaluated identity-transformation request
    /// (derivative, integral, solve, ...) that the scheduler tries to reduce
    /// before anything else and that closed-form extraction refuses.
    pub const HEAVY: Self = Self(1 << 2);
    /// The node is a leaf carrying a payload.
    pub const LEAF: Self = Self(1 << 3);
    /// No flags.
    pub const NONE: Self = Self(0);

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

/// Scalar numeric semantics of an operator: maps child values to a value.
///
/// A plain function pointer on purpose — it is what the interpreter calls
/// and what a JIT backend can emit a direct call to.
pub type EvalFn = fn(&[f64]) -> f64;

/// Everything the kernel and the backends need to know about an operator.
#[derive(Clone, Debug)]
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
}

impl OpDescriptor {
    /// A plain operator with no special properties.
    ///
    /// Its cost is 3: a named function weighs more than an arithmetic
    /// operator (cost 1) with a literal operand (cost 1), so `-x` is
    /// preferred to `abs(x)` when both are known to be equal.
    #[must_use]
    pub fn new(
        name: &str,
        arity: Arity,
    ) -> Self {
        Self {
            name: Arc::from(name),
            arity,
            flags: OpFlags::NONE,
            cost: 3,
            binder: None,
            eval: None,
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
    #[must_use]
    pub const fn eval(
        mut self,
        eval: EvalFn,
    ) -> Self {
        self.eval = Some(eval);
        self
    }
}

/// Error returned when an operator registration conflicts with an existing
/// one of the same name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpConflict {
    /// The offending name.
    pub name: String,
}

impl std::fmt::Display for OpConflict {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(
            f,
            "operator `{}` is already registered with a different signature",
            self.name
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
    by_name: HashMap<Arc<str>, OpId>,
    attrs: HashMap<(OpId, TypeId), Arc<dyn Any + Send + Sync>>,
}

impl std::fmt::Debug for OpTable {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
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
            attrs: HashMap::new(),
        };
        let ac = OpFlags::COMMUTATIVE.with(OpFlags::ASSOCIATIVE);
        let builtins = [
            OpDescriptor::new("lit", Arity::Fixed(0))
                .flags(OpFlags::LEAF)
                .cost(1),
            // A symbol costs more than a number and an operator together,
            // so that `2*x` beats `x + x` and `x^2` beats `x*x`.
            OpDescriptor::new("sym", Arity::Fixed(0))
                .flags(OpFlags::LEAF)
                .cost(3),
            OpDescriptor::new("add", Arity::Variadic)
                .flags(ac)
                .cost(1)
                .eval(|a| a.iter().sum()),
            OpDescriptor::new("mul", Arity::Variadic)
                .flags(ac)
                .cost(1)
                .eval(|a| a.iter().product()),
            OpDescriptor::new("pow", Arity::Fixed(2))
                .cost(1)
                .eval(|a| match a {
                    | [b, e] => b.powf(*e),
                    | _ => f64::NAN,
                }),
            OpDescriptor::new("apply", Arity::Variadic).cost(2),
            OpDescriptor::new("list", Arity::Variadic).cost(1),
            OpDescriptor::new("eq", Arity::Fixed(2)).cost(1),
        ];
        for desc in builtins {
            table.push(desc);
        }
        table
    }

    fn push(
        &mut self,
        desc: OpDescriptor,
    ) -> OpId {
        let id = OpId(to_u32(self.ops.len()));
        self.by_name.insert(Arc::clone(&desc.name), id);
        self.ops.push(desc);
        id
    }

    /// Registers an operator, or returns the existing id when an operator
    /// with the same name, arity and flags is already present.
    ///
    /// # Errors
    /// Returns [`OpConflict`] when the name is taken by an operator with a
    /// different arity, flag set or binder.
    pub fn register(
        &mut self,
        desc: OpDescriptor,
    ) -> Result<OpId, OpConflict> {
        if let Some(&id) = self.by_name.get(&desc.name) {
            let same = self.ops.get(id.index()).is_some_and(|old| {
                old.arity == desc.arity && old.flags == desc.flags && old.binder == desc.binder
            });
            return if same {
                Ok(id)
            } else {
                Err(OpConflict {
                    name: desc.name.to_string(),
                })
            };
        }
        Ok(self.push(desc))
    }

    /// Looks an operator up by name.
    #[must_use]
    pub fn lookup(
        &self,
        name: &str,
    ) -> Option<OpId> {
        self.by_name.get(name).copied()
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
}
