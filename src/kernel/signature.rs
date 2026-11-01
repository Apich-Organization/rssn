//! Typed signatures for operators.
//!
//! Every operator *may* carry a list of [`Signature`]s. A signature says what
//! the operator accepts and returns, in a form a compiler front end can read
//! out of the [manifest](super::manifest): parameter and result types,
//! dimension rules, side conditions, effects and variadic shape.
//!
//! Signatures are descriptive, never prescriptive: the engine does not check
//! arguments against them. That is deliberate. A signature is metadata for
//! *generated bindings* and for documentation, while argument checking at
//! run time is done on [`Value`](super::Value)s, where a mismatch is a real
//! error with a real message. An operator with no signature is still
//! callable through the dynamic ABI; the manifest records *why* it has none
//! in `dynamic_reason`.
//!
//! # Stable keys
//!
//! [`OpKey`] identifies an operator across builds: the rule set that
//! registered it, its name, and the hash of the overload's signature. Numeric
//! [`OpId`](crate::graph::OpId)s are registration-order dependent and must
//! never be persisted; `OpKey` is what a cache key or a generated binding
//! stores.

use std::fmt;

use serde::Deserialize;
use serde::Serialize;

use crate::graph::dim::DimRule;
use crate::graph::dim::FamilySource;

/// Version of the dynamic ABI described by this module.
///
/// A change that would break an existing host — a removed field, a changed
/// meaning, a different `Value` layout — increments this number. Additive
/// changes do not.
pub const API_VERSION: u32 = 1;

/// Version of the plugin dynamic-library ABI.
pub const PLUGIN_ABI: u32 = 1;

// ---------------------------------------------------------------------------
// Numbers in types
// ---------------------------------------------------------------------------

/// An exact rational number written in a type: a dimension exponent or a
/// shape component.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Rational {
    /// Numerator, carrying the sign.
    pub numer: i64,
    /// Denominator, always positive.
    pub denom: i64,
}

impl Rational {
    /// The rational `numer/denom`, reduced and with a positive denominator.
    ///
    /// # Panics
    /// Panics on a zero denominator: a signature with a zero denominator is
    /// a bug in the rule set that built it.
    #[must_use]
    pub fn new(
        numer: i64,
        denom: i64,
    ) -> Self {
        assert!(denom != 0, "rssn: a type exponent cannot have denominator 0");
        let g = gcd(numer, denom);
        let sign = if denom < 0 { -1 } else { 1 };
        Self {
            numer: sign * numer / g,
            denom: sign * denom / g,
        }
    }

    /// An integer as a rational.
    #[must_use]
    pub const fn int(value: i64) -> Self {
        Self {
            numer: value,
            denom: 1,
        }
    }

    /// Whether this is the integer `1`.
    #[must_use]
    pub const fn is_one(&self) -> bool {
        self.numer == 1 && self.denom == 1
    }

    /// Whether this is the integer `0`.
    #[must_use]
    pub const fn is_zero(&self) -> bool {
        self.numer == 0
    }
}

impl fmt::Display for Rational {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        if self.denom == 1 {
            write!(f, "{}", self.numer)
        } else {
            write!(f, "{}/{}", self.numer, self.denom)
        }
    }
}

/// Greatest common divisor of two integers, always positive.
const fn gcd(
    a: i64,
    b: i64,
) -> i64 {
    let (mut a, mut b) = (a.unsigned_abs(), b.unsigned_abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 { 1 } else { a.cast_signed() }
}

/// An exponent in a dimension expression: a rational, or a symbol standing
/// for an integer (`x^n`).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Exp {
    /// A rational exponent.
    Rat(Rational),
    /// A symbol whose value is an integer exponent.
    Symbol(String),
}

impl From<i64> for Exp {
    fn from(value: i64) -> Self {
        Self::Rat(Rational::int(value))
    }
}

impl fmt::Display for Exp {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Rat(r) => write!(f, "{r}"),
            | Self::Symbol(s) => f.write_str(s),
        }
    }
}

// ---------------------------------------------------------------------------
// Type expressions
// ---------------------------------------------------------------------------

/// A dimension written in a type.
///
/// This is the type-level mirror of
/// [`DimTerm`](crate::graph::dim::DimTerm): products of named base
/// dimensions and dimension variables with exact exponents. It is a *type*
/// in the sense that two terms with the same `DimExpr` are dimensionally
/// interchangeable; the run-time unifier works on `DimTerm`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DimExpr {
    /// Any dimension: the polymorphic scalar `0`.
    Any,
    /// Not a scalar at all (a truth value, a label, a container).
    NonMetric,
    /// A structured container of mixed parts.
    Container,
    /// Dimensionless.
    Dimless,
    /// A base dimension by name (`Length`).
    Base(String),
    /// A dimension variable bound by a type parameter (`?D`).
    Var(String),
    /// A product of factors.
    Mul(Vec<Self>),
    /// A factor raised to an exponent.
    Pow(Box<Self>, Exp),
    /// The reciprocal.
    Inverse(Box<Self>),
}

impl DimExpr {
    /// The dimension of a length.
    #[must_use]
    pub fn base(name: &str) -> Self {
        Self::Base(name.to_owned())
    }

    /// `self^exponent`.
    #[must_use]
    pub fn pow(
        self,
        exponent: impl Into<Exp>,
    ) -> Self {
        let exponent = exponent.into();
        if exponent == Exp::Rat(Rational::int(1)) {
            return self;
        }
        Self::Pow(Box::new(self), exponent)
    }
}

/// A natural number written in a type: an array rank or extent.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NatExpr {
    /// A literal.
    Lit(u64),
    /// A type parameter.
    Var(String),
    /// A sum.
    Sum(Box<Self>, Box<Self>),
    /// A product.
    Product(Box<Self>, Box<Self>),
}

impl From<u64> for NatExpr {
    fn from(value: u64) -> Self {
        Self::Lit(value)
    }
}

/// A row shape for tables: named columns and a row count.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct RowExpr {
    /// Column names and their types, in declaration order.
    pub columns: Vec<(String, TypeExpr)>,
    /// Number of rows.
    pub rows: NatExpr,
}

/// The run-time representation a `Dyn` value carries.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct RepExpr {
    /// Name of the representation (`"term"`, `"cbor"`, `"f64"`).
    pub repr: String,
    /// The dimension of the value, when it has one.
    pub dim: Option<DimExpr>,
}

/// The type of an argument or a result.
///
/// Scalar variants carry no dimension; wrap them in [`TypeExpr::Dim`] to
/// attach one (`Real` of `Length/Time`). That keeps the common dimensionless
/// case short while still allowing a fully dimensioned signature.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TypeExpr {
    /// An exact integer.
    Int,
    /// An exact rational, or a symbolic real.
    Real,
    /// A double-precision float.
    Float,
    /// A complex number.
    Complex,
    /// A truth value.
    Bool,
    /// Text.
    String,
    /// A unit definition (a scale and a dimension), not a quantity.
    Unit(DimExpr),
    /// An affine point (`degC`) of the given dimension.
    Point(DimExpr),
    /// A scalar with a dimension: `Real[Length]`.
    Dim(Box<Self>, DimExpr),
    /// A dense tensor.
    Tensor {
        /// Rank of the tensor.
        rank: NatExpr,
        /// Extents, when every dimension is known.
        shape: Option<Vec<NatExpr>>,
        /// Element type.
        elem: Box<Self>,
    },
    /// A sparse tensor.
    Sparse {
        /// Storage format (`"csr"`, `"csc"`, `"coo"`).
        format: String,
        /// Rank of the tensor.
        rank: NatExpr,
        /// Element type.
        elem: Box<Self>,
    },
    /// A record with named fields.
    Record(Vec<(String, Self)>),
    /// A tagged union: tag and its payload types.
    Variant(Vec<(String, Vec<Self>)>),
    /// A function from arguments to a result.
    Function(Vec<Self>, Box<Self>),
    /// A table.
    Table(RowExpr),
    /// A value whose representation is only known at run time.
    Dyn(RepExpr),
    /// A device handle.
    Device,
    /// A host-defined type, by name.
    Opaque(String),
    /// A type bound by a [`TypeParam`].
    Var(String),
    /// The type of an expression that never returns (a failed call).
    Never,
}

impl TypeExpr {
    /// `self` with dimension `dim`.
    #[must_use]
    pub fn with_dim(
        self,
        dim: DimExpr,
    ) -> Self {
        Self::Dim(Box::new(self), dim)
    }

    /// A vector of `self` with the given length.
    #[must_use]
    pub fn vector(
        self,
        len: impl Into<NatExpr>,
    ) -> Self {
        Self::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![len.into()]),
            elem: Box::new(self),
        }
    }

    /// A matrix of `self`.
    #[must_use]
    pub fn matrix(
        self,
        rows: impl Into<NatExpr>,
        cols: impl Into<NatExpr>,
    ) -> Self {
        Self::Tensor {
            rank: NatExpr::Lit(2),
            shape: Some(vec![rows.into(), cols.into()]),
            elem: Box::new(self),
        }
    }

    /// A tensor of unknown rank.
    #[must_use]
    pub fn tensor(
        self,
        rank: impl Into<NatExpr>,
    ) -> Self {
        Self::Tensor {
            rank: rank.into(),
            shape: None,
            elem: Box::new(self),
        }
    }
}

// ---------------------------------------------------------------------------
// Side conditions and effects
// ---------------------------------------------------------------------------

/// What a type parameter may stand for.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TypeParamKind {
    /// Any type.
    Type,
    /// A natural number.
    Nat,
    /// A dimension.
    Dim,
    /// A lifetime, for host handles.
    Lifetime,
}

/// A type parameter of a signature.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct TypeParam {
    /// Name, as written in [`TypeExpr::Var`].
    pub name: String,
    /// What the parameter stands for.
    pub kind: TypeParamKind,
    /// Named bounds the parameter must satisfy.
    pub bounds: Vec<String>,
}

/// A class of types, used by [`Constraint`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConstraintKind {
    /// Values can be compared for equality.
    Eq,
    /// Values are totally ordered.
    Ord,
    /// Values can be printed.
    Show,
    /// Values can be converted to a float.
    Numeric,
    /// Values form a ring.
    Ring,
    /// Values form a field.
    Field,
    /// Values are integers.
    Integral,
    /// Values are real.
    Real,
    /// The type is dimensionally sound.
    Dimensioned,
}

/// A side condition of a signature.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Constraint {
    /// The class the arguments must belong to.
    pub kind: ConstraintKind,
    /// Type variables or expressions the constraint applies to.
    pub args: Vec<String>,
}

/// An effect an operator may have.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Effect {
    /// No observable effect.
    Pure,
    /// Reads or writes outside the process.
    Io,
    /// May fail to terminate within the budget.
    Diverges,
    /// Result depends on state outside the arguments.
    Stateful,
}

/// The effect row of a signature.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EffectRow {
    /// Effects the operator may perform. Empty means pure.
    pub effects: Vec<Effect>,
}

impl EffectRow {
    /// The pure row.
    #[must_use]
    pub const fn pure() -> Self {
        Self {
            effects: Vec::new(),
        }
    }

    /// Whether the operator is pure.
    #[must_use]
    pub const fn is_pure(&self) -> bool {
        self.effects.is_empty()
    }
}

/// How many trailing arguments an operator accepts.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Variadic {
    /// Fewest accepted arguments, counted after the fixed ones.
    pub min: u8,
    /// Most accepted arguments, if bounded.
    pub max: Option<u8>,
    /// Type every trailing argument must have.
    pub elem: TypeExpr,
}

// ---------------------------------------------------------------------------
// Signature
// ---------------------------------------------------------------------------

/// How many arguments an operator accepts, as manifest data.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ArityInfo {
    /// Fewest accepted arguments.
    pub min: u8,
    /// Most accepted arguments, or `None` for unbounded.
    pub max: Option<u8>,
}

impl ArityInfo {
    /// The shape of a fixed-arity operator.
    #[must_use]
    pub const fn fixed(n: u8) -> Self {
        Self {
            min: n,
            max: Some(n),
        }
    }

    /// The shape of an operator taking any number of arguments.
    #[must_use]
    pub const fn variadic() -> Self {
        Self { min: 0, max: None }
    }

    /// The shape of a crate-level [`Arity`](crate::graph::Arity).
    #[must_use]
    pub const fn from_arity(arity: crate::graph::Arity) -> Self {
        match arity {
            | crate::graph::Arity::Fixed(k) => Self::fixed(k),
            | crate::graph::Arity::Variadic => Self::variadic(),
        }
    }

    /// Whether `n` arguments are acceptable.
    #[must_use]
    pub fn accepts(
        &self,
        n: usize,
    ) -> bool {
        let n = u8::try_from(n).unwrap_or(u8::MAX);
        n >= self.min && self.max.is_none_or(|max| n <= max)
    }
}

impl fmt::Display for ArityInfo {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self.max {
            | Some(max) if max == self.min => write!(f, "{}", self.min),
            | Some(max) => write!(f, "{}..={max}", self.min),
            | None => write!(f, "{}..", self.min),
        }
    }
}

/// A short human-readable rendering of a type, for diagnostics.
#[must_use]
pub fn type_name(expr: &TypeExpr) -> String {
    match expr {
        | TypeExpr::Int => "Int".to_owned(),
        | TypeExpr::Real => "Real".to_owned(),
        | TypeExpr::Float => "Float".to_owned(),
        | TypeExpr::Complex => "Complex".to_owned(),
        | TypeExpr::Bool => "Bool".to_owned(),
        | TypeExpr::String => "String".to_owned(),
        | TypeExpr::Device => "Device".to_owned(),
        | TypeExpr::Unit(d) => format!("Unit[{}]", dim_name(d)),
        | TypeExpr::Point(d) => format!("Point[{}]", dim_name(d)),
        | TypeExpr::Dim(inner, dim) => format!("{}[{}]", type_name(inner), dim_name(dim)),
        | TypeExpr::Tensor { rank, elem, .. } => format!("Tensor<{}, {}>", type_name(elem), nat_name(rank)),
        | TypeExpr::Sparse { format, elem, .. } => {
            format!("Sparse<{format}, {}>", type_name(elem))
        },
        | TypeExpr::Record(fields) => format!("Record({})", fields.len()),
        | TypeExpr::Variant(variants) => format!("Variant({})", variants.len()),
        | TypeExpr::Function(args, ret) => format!(
            "({}) -> {}",
            args.iter().map(type_name).collect::<Vec<_>>().join(", "),
            type_name(ret)
        ),
        | TypeExpr::Table(row) => format!("Table({})", row.columns.len()),
        | TypeExpr::Dyn(rep) => format!("Dyn<{}>", rep.repr),
        | TypeExpr::Opaque(name) => format!("Opaque<{name}>"),
        | TypeExpr::Var(name) => name.clone(),
        | TypeExpr::Never => "Never".to_owned(),
    }
}

/// A short rendering of a dimension expression.
#[must_use]
pub fn dim_name(dim: &DimExpr) -> String {
    match dim {
        | DimExpr::Any => "any".to_owned(),
        | DimExpr::NonMetric => "non-metric".to_owned(),
        | DimExpr::Container => "container".to_owned(),
        | DimExpr::Dimless => "1".to_owned(),
        | DimExpr::Base(name) | DimExpr::Var(name) => name.clone(),
        | DimExpr::Mul(factors) => factors.iter().map(dim_name).collect::<Vec<_>>().join("*"),
        | DimExpr::Pow(base, exp) => format!("{}^{exp}", dim_name(base)),
        | DimExpr::Inverse(base) => format!("1/{}", dim_name(base)),
    }
}

/// A short rendering of a shape expression.
#[must_use]
pub fn nat_name(nat: &NatExpr) -> String {
    match nat {
        | NatExpr::Lit(n) => n.to_string(),
        | NatExpr::Var(name) => name.clone(),
        | NatExpr::Sum(a, b) => format!("({}+{})", nat_name(a), nat_name(b)),
        | NatExpr::Product(a, b) => format!("({}*{})", nat_name(a), nat_name(b)),
    }
}

impl fmt::Display for TypeExpr {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&type_name(self))
    }
}

impl fmt::Display for DimExpr {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&dim_name(self))
    }
}

impl fmt::Display for NatExpr {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&nat_name(self))
    }
}

/// The type of one overload of an operator.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Signature {
    /// Type parameters introduced by this overload.
    pub type_params: Vec<TypeParam>,
    /// Parameter types, in order.
    pub params: Vec<TypeExpr>,
    /// Result type.
    pub ret: TypeExpr,
    /// Side conditions.
    pub constraints: Vec<Constraint>,
    /// Effects.
    pub effects: EffectRow,
    /// The dimension rule, mirroring [`DimRule`].
    pub dim: Option<DimRuleInfo>,
    /// Trailing arguments, for variadic operators.
    pub variadic: Option<Variadic>,
}

impl Signature {
    /// A signature from parameter types to a result type, with no
    /// dimension rule and no effects.
    #[must_use]
    pub const fn new(
        params: Vec<TypeExpr>,
        ret: TypeExpr,
    ) -> Self {
        Self {
            type_params: Vec::new(),
            params,
            ret,
            constraints: Vec::new(),
            effects: EffectRow::pure(),
            dim: None,
            variadic: None,
        }
    }

    /// `(Scalar, ..., Scalar) -> Scalar` for `n` arguments.
    #[must_use]
    pub fn scalars(
        n: usize,
        ret: TypeExpr,
    ) -> Self {
        Self::new(vec![TypeExpr::Real; n], ret)
    }

    /// Attaches a dimension rule.
    #[must_use]
    pub fn with_dim(
        mut self,
        rule: DimRule,
    ) -> Self {
        self.dim = Some(DimRuleInfo::from(rule));
        self
    }

    /// Declares trailing arguments.
    #[must_use]
    pub fn variadic(
        mut self,
        variadic: Variadic,
    ) -> Self {
        self.variadic = Some(variadic);
        self
    }

    /// Number of fixed parameters.
    #[must_use]
    pub const fn arity(&self) -> usize {
        self.params.len()
    }

    /// Whether this overload accepts `n` arguments.
    #[must_use]
    pub fn accepts_arity(
        &self,
        n: usize,
    ) -> bool {
        match &self.variadic {
            | None => self.params.len() == n,
            | Some(v) => {
                let fixed = self.params.len();
                let extra = n.saturating_sub(fixed);
                let low = usize::from(v.min);
                let high = v.max.map_or(usize::MAX, usize::from);
                n >= fixed && extra >= low && extra <= high
            },
        }
    }

    /// A 16-byte digest of this signature, used in [`OpKey`].
    ///
    /// The digest is taken over the deterministic binary encoding of the
    /// signature, so it is stable across processes and platforms.
    #[must_use]
    pub fn hash(&self) -> [u8; 16] {
        let bytes = super::codec::encode_deterministic(self).unwrap_or_default();
        let full = *blake3::hash(&bytes).as_bytes();
        let mut out = [0_u8; 16];
        out.copy_from_slice(full.get(..16).unwrap_or(&full));
        out
    }
}

/// A stable identifier for one overload of one operator.
///
/// Two builds that register the same operator with the same signature
/// produce the same key, whatever the registration order. Hosts key caches
/// on `(manifest_hash, OpKey)`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct OpKey {
    /// Rule set that registered the operator.
    pub rule_set: String,
    /// Operator name.
    pub name: String,
    /// Digest of the overload's signature; all zero when it has none.
    pub signature_hash: [u8; 16],
}

impl OpKey {
    /// A key for `rule_set/name` with the given signature.
    #[must_use]
    pub fn new(
        rule_set: &str,
        name: &str,
        signature: Option<&Signature>,
    ) -> Self {
        Self {
            rule_set: rule_set.to_owned(),
            name: name.to_owned(),
            signature_hash: signature.map_or([0_u8; 16], Signature::hash),
        }
    }

    /// The dotted spelling used in diagnostics: `rule_set.name`.
    #[must_use]
    pub fn qualified(&self) -> String {
        format!("{}.{}", self.rule_set, self.name)
    }
}

impl fmt::Display for OpKey {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&self.qualified())
    }
}

// ---------------------------------------------------------------------------
// Operator metadata
// ---------------------------------------------------------------------------

/// How settled an operator's behaviour is.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Stability {
    /// Semantics and signature do not change within a major version.
    #[default]
    Stable,
    /// May change; the old name is kept for at least one minor version.
    Experimental,
    /// Scheduled for removal; warns for at least one minor version.
    Deprecated,
}

/// Whether an operator touches anything but its arguments.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Purity {
    /// A function of its arguments alone.
    #[default]
    Pure,
    /// Reads state it does not change.
    ReadOnly,
    /// May change state.
    Effectful,
}

/// Whether repeated calls agree.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Determinism {
    /// Same arguments, same result, always.
    #[default]
    Deterministic,
    /// Uses randomness (Monte Carlo, random restarts).
    Randomized,
    /// Depends on the environment (thread count, hardware, time).
    EnvDependent,
}

/// How an operator may be used across threads.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThreadSafety {
    /// Callable from any thread; the computation carries no session state.
    Send,
    /// Must run on the session that owns the terms it is given.
    #[default]
    SessionBound,
    /// Keeps mutable state of its own; needs external locking.
    NotThreadSafe,
}

/// How an operator produces values.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvalKind {
    /// No direct evaluation: only rules can reduce it.
    #[default]
    None,
    /// A scalar `fn(&[f64]) -> f64` fast path.
    Scalar,
    /// A dynamic evaluator that sees [`Value`](super::Value)s.
    Dynamic,
    /// A callback into the host.
    HostCallback,
}

/// Which backends can execute an operator.
///
/// One flag per backend, deliberately: a host asks "can the JIT run this",
/// not "what is the backend mask".
#[allow(clippy::struct_excessive_bools)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BackendSupport {
    /// The tree-walking interpreter.
    pub interpreter: bool,
    /// The Cranelift JIT (feature `jit`).
    pub jit: bool,
    /// The C ABI.
    pub ffi: bool,
    /// A GPU backend, when one is built.
    pub gpu: bool,
}

impl Default for BackendSupport {
    fn default() -> Self {
        Self {
            interpreter: true,
            jit: false,
            ffi: true,
            gpu: false,
        }
    }
}

impl BackendSupport {
    /// Every backend.
    #[must_use]
    pub const fn all() -> Self {
        Self {
            interpreter: true,
            jit: true,
            ffi: true,
            gpu: false,
        }
    }
}

/// A worked example, shown in generated documentation.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Example {
    /// The call, in surface syntax: `diff(x^3, x)`.
    pub input: String,
    /// What it produces: `3*x^2`.
    pub output: String,
    /// Why the example matters, if that is not obvious.
    pub note: String,
}

/// The algebraic and scheduling flags of an operator, as data.
///
/// Independent properties, so independent fields: a host reads
/// `flags.commutative`, it does not mask bits.
#[allow(clippy::struct_excessive_bools)]
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct OpFlagsInfo {
    /// `f(f(a, b), c) = f(a, b, c)`.
    pub associative: bool,
    /// Arguments may be reordered.
    pub commutative: bool,
    /// An unevaluated request the scheduler reduces first.
    pub heavy: bool,
    /// A leaf carrying a payload.
    pub leaf: bool,
    /// A truth value.
    pub predicate: bool,
    /// Windows look inside even when the node is shared.
    pub transparent: bool,
    /// Irreducible when applied to an undetermined function.
    pub opaque_on_apply: bool,
}

/// The binding structure of an operator, as data.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct BinderInfo {
    /// Index of the child holding the bound symbol.
    pub var: u8,
    /// Bit `i` set means the variable is bound inside child `i`.
    pub scope: u32,
}

/// A [`DimRule`] as serialisable data.
///
/// The run-time rule carries function pointers in two variants; those become
/// the name of the registering rule set here.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DimRuleInfo {
    /// No rule: a coverage gap.
    Unspecified,
    /// Not a scalar.
    NonMetric,
    /// A structured container.
    Container,
    /// Dimensionally polymorphic.
    Any,
    /// The dimension of argument `i`.
    Same(u8),
    /// All arguments agree; the result is their common dimension.
    Sum,
    /// Product of the argument dimensions.
    Product,
    /// Argument 0 divided by argument 1.
    Quotient,
    /// `a^b` with `b` dimensionless.
    Power,
    /// The `n`-th root of argument 0.
    Root(u32),
    /// Everything is dimensionless.
    Dimensionless,
    /// `1 / dim(0)`.
    Inverse,
    /// `dim(0)^2`.
    Square,
    /// The common dimension of the leaves of a list.
    List,
    /// Arguments agree; the result is not a scalar.
    Compare,
    /// Arguments agree; the result is a dimensionless angle.
    Angle,
    /// A fresh dimension variable.
    UnknownFunction,
    /// `dim(f) / dim(x)`.
    Diff {
        /// Index of the function.
        f: u8,
        /// Index of the variable.
        x: u8,
    },
    /// `dim(f) / dim(x)^n`.
    DiffN {
        /// Index of the function.
        f: u8,
        /// Index of the variable.
        x: u8,
        /// Index of the order.
        n: u8,
    },
    /// `dim(f) * dim(x)`.
    Integral {
        /// Index of the integrand.
        f: u8,
        /// Index of the variable.
        x: u8,
    },
    /// The dimension of the summand.
    SumOver {
        /// Index of the summand.
        f: u8,
        /// Index of the index variable.
        index: u8,
    },
    /// `dim(f)^(b - a + 1)`.
    ProductOver {
        /// Index of the factor.
        f: u8,
        /// Index of the index variable.
        index: u8,
    },
    /// A distribution family.
    Family {
        /// Where the value's dimension comes from.
        result: FamilySourceInfo,
        /// Arguments sharing the value's dimension.
        same: Vec<u8>,
        /// Arguments that must be pure numbers.
        dimensionless: Vec<u8>,
    },
    /// A matrix function.
    MatrixFunction {
        /// Index of the scalar function.
        function: u8,
        /// Index of its variable.
        variable: u8,
        /// Index of the matrix.
        matrix: u8,
    },
    /// `quantity(v, u)`.
    Quantity,
    /// A domain-specific rule, by the name of the rule set that installed it.
    Custom(String),
}

/// Where a [`DimRuleInfo::Family`] value's dimension comes from.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FamilySourceInfo {
    /// The dimension of argument `i`.
    Argument(u8),
    /// The inverse dimension of argument `i`.
    Inverse(u8),
}

fn index(
    value: usize,
    what: &'static str,
) -> u8 {
    u8::try_from(value).unwrap_or_else(|_| {
        panic!("rssn: {what} index {value} does not fit in a signature")
    })
}

impl From<DimRule> for DimRuleInfo {
    fn from(rule: DimRule) -> Self {
        match rule {
            | DimRule::Unspecified => Self::Unspecified,
            | DimRule::NonMetric => Self::NonMetric,
            | DimRule::Container => Self::Container,
            | DimRule::Any => Self::Any,
            | DimRule::Same(i) => Self::Same(index(i, "argument")),
            | DimRule::Sum => Self::Sum,
            | DimRule::Product => Self::Product,
            | DimRule::Quotient => Self::Quotient,
            | DimRule::Power => Self::Power,
            | DimRule::Root(n) => Self::Root(n),
            | DimRule::Dimensionless => Self::Dimensionless,
            | DimRule::Inverse => Self::Inverse,
            | DimRule::Square => Self::Square,
            | DimRule::List => Self::List,
            | DimRule::Compare => Self::Compare,
            | DimRule::Angle => Self::Angle,
            | DimRule::UnknownFunction => Self::UnknownFunction,
            | DimRule::Diff { f, x } => Self::Diff {
                f: index(f, "function"),
                x: index(x, "variable"),
            },
            | DimRule::DiffN { f, x, n } => Self::DiffN {
                f: index(f, "function"),
                x: index(x, "variable"),
                n: index(n, "order"),
            },
            | DimRule::Integral { f, x } => Self::Integral {
                f: index(f, "integrand"),
                x: index(x, "variable"),
            },
            | DimRule::SumOver { f, index: i } => Self::SumOver {
                f: index(f, "summand"),
                index: index(i, "index variable"),
            },
            | DimRule::ProductOver { f, index: i } => Self::ProductOver {
                f: index(f, "factor"),
                index: index(i, "index variable"),
            },
            | DimRule::Family {
                result,
                same,
                dimensionless,
            } => Self::Family {
                result: match result {
                    | FamilySource::Argument(i) => FamilySourceInfo::Argument(index(i, "argument")),
                    | FamilySource::Inverse(i) => FamilySourceInfo::Inverse(index(i, "argument")),
                },
                same: same.iter().map(|&i| index(i, "argument")).collect(),
                dimensionless: dimensionless
                    .iter()
                    .map(|&i| index(i, "argument"))
                    .collect(),
            },
            | DimRule::MatrixFunction {
                function,
                variable,
                matrix,
            } => Self::MatrixFunction {
                function: index(function, "function"),
                variable: index(variable, "variable"),
                matrix: index(matrix, "matrix"),
            },
            | DimRule::Quantity => Self::Quantity,
            // A `Custom` rule is a function pointer; the manifest records the
            // operator's own name instead, which is what a host can act on.
            | DimRule::Custom(_) => Self::Custom(String::new()),
        }
    }
}

impl DimRuleInfo {
    /// The rule a host would install, when it is expressible as data.
    ///
    /// [`DimRuleInfo::Custom`] has no data equivalent and becomes
    /// [`DimRule::Unspecified`]: a host that meets one keeps the operator
    /// dimensionally opaque rather than guessing.
    #[must_use]
    pub fn to_rule(&self) -> DimRule {
        match self {
            | Self::Unspecified | Self::Family { .. } | Self::Custom(_) => DimRule::Unspecified,
            | Self::NonMetric => DimRule::NonMetric,
            | Self::Container => DimRule::Container,
            | Self::Any => DimRule::Any,
            | Self::Same(i) => DimRule::Same(usize::from(*i)),
            | Self::Sum => DimRule::Sum,
            | Self::Product => DimRule::Product,
            | Self::Quotient => DimRule::Quotient,
            | Self::Power => DimRule::Power,
            | Self::Root(n) => DimRule::Root(*n),
            | Self::Dimensionless => DimRule::Dimensionless,
            | Self::Inverse => DimRule::Inverse,
            | Self::Square => DimRule::Square,
            | Self::List => DimRule::List,
            | Self::Compare => DimRule::Compare,
            | Self::Angle => DimRule::Angle,
            | Self::UnknownFunction => DimRule::UnknownFunction,
            | Self::Diff { f, x } => DimRule::Diff {
                f: usize::from(*f),
                x: usize::from(*x),
            },
            | Self::DiffN { f, x, n } => DimRule::DiffN {
                f: usize::from(*f),
                x: usize::from(*x),
                n: usize::from(*n),
            },
            | Self::Integral { f, x } => DimRule::Integral {
                f: usize::from(*f),
                x: usize::from(*x),
            },
            | Self::SumOver { f, index } => DimRule::SumOver {
                f: usize::from(*f),
                index: usize::from(*index),
            },
            | Self::ProductOver { f, index } => DimRule::ProductOver {
                f: usize::from(*f),
                index: usize::from(*index),
            },
            | Self::MatrixFunction {
                function,
                variable,
                matrix,
            } => DimRule::MatrixFunction {
                function: usize::from(*function),
                variable: usize::from(*variable),
                matrix: usize::from(*matrix),
            },
            | Self::Quantity => DimRule::Quantity,
        }
    }
}

/// Why an operator has no signature, when it has none.
///
/// The audit binary fails the build on an operator that has neither a
/// signature nor one of these reasons, so no operator is silently
/// untyped.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DynamicReason {
    /// It takes or returns a container whose shape depends on its arguments.
    ShapeDependent,
    /// It is a request reduced by rules, and its result type is decided by
    /// the reduction.
    Request,
    /// It is polymorphic over every argument type.
    FullyPolymorphic,
    /// It carries an opaque host payload.
    OpaquePayload,
    /// Its arity is not fixed.
    Variadic,
    /// The signature is not written yet.
    NotYetSpecified,
}

impl DynamicReason {
    /// The stable spelling used in the manifest.
    #[must_use]
    pub const fn as_str(&self) -> &'static str {
        match self {
            | Self::ShapeDependent => "shape_dependent",
            | Self::Request => "request",
            | Self::FullyPolymorphic => "fully_polymorphic",
            | Self::OpaquePayload => "opaque_payload",
            | Self::Variadic => "variadic",
            | Self::NotYetSpecified => "not_yet_specified",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rationals_are_normalised() {
        assert_eq!(Rational::new(2, 4), Rational::new(1, 2));
        assert_eq!(Rational::new(-2, -4), Rational::new(1, 2));
        assert_eq!(Rational::new(2, -4), Rational::new(-1, 2));
        assert_eq!(Rational::new(0, -7), Rational::int(0));
        assert_eq!(Rational::new(3, 1), Rational::int(3));
    }

    #[test]
    fn signature_hashes_are_stable_and_distinguishing() {
        let a = Signature::scalars(2, TypeExpr::Real);
        let b = Signature::scalars(2, TypeExpr::Real);
        let c = Signature::scalars(3, TypeExpr::Real);
        assert_eq!(a.hash(), b.hash());
        assert_ne!(a.hash(), c.hash());
        // The digest is over the encoding, so a dimension change shows up.
        let dimmed = Signature::new(
            vec![TypeExpr::Real.with_dim(DimExpr::base("Length"))],
            TypeExpr::Real,
        );
        assert_ne!(dimmed.hash(), Signature::scalars(1, TypeExpr::Real).hash());
    }

    #[test]
    fn variadic_arity() {
        let sig = Signature::scalars(1, TypeExpr::Real).variadic(Variadic {
            min: 0,
            max: None,
            elem: TypeExpr::Real,
        });
        assert!(sig.accepts_arity(1));
        assert!(sig.accepts_arity(9));
        assert!(!sig.accepts_arity(0));
    }

    #[test]
    fn dim_rule_round_trips() {
        for rule in [
            DimRule::Sum,
            DimRule::Diff { f: 0, x: 1 },
            DimRule::Root(2),
            DimRule::NonMetric,
            DimRule::Quantity,
        ] {
            let info = DimRuleInfo::from(rule);
            // `DimRule` itself has no structural equality (it carries
            // static slices), so the round trip is checked on the
            // owned description, which is what a manifest compares.
            assert_eq!(DimRuleInfo::from(info.to_rule()), info);
        }
    }
}
