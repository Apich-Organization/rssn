//! Dimensions as a typing layer over the term graph.
//!
//! Every scalar term has a dimension. Base dimensions are interned names
//! (`Length`, `Mass`, `Time`, ... and any user-defined base); exponents are
//! exact rationals or, for symbolic exponents, rational functions over
//! interned atoms (`x^n`, `x^(a+b)`, `x^(1/n)`).
//!
//! Operator signatures are [`DimRule`] attributes. The inference engine in
//! `infer` walks the hash-consed DAG, unifies the dimensions of all members
//! of an e-class and reports conflicts. The public entry points are
//! [`Graph::dimension`] and [`Graph::assume_dimension`]; `Config` and
//! `Answer` expose the result.
//!
//! The design is documented in `docs/DIMENSIONS.md`.

// The exponent field and the inference engine are implementation details of
// this module; `pub(crate)` is already the widest visibility reachable from
// their private submodules, so the lint's suggested `pub(super)` would force
// re-exports at every use site.
#![allow(clippy::redundant_pub_crate)]

pub(crate) mod expr;
pub(crate) mod infer;

use std::collections::HashMap;
use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;

pub(crate) use expr::AtomId;
pub(crate) use expr::Expo;
pub(crate) use infer::DimInfer;

/// Identifier of a base dimension inside a [`DimArena`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseId(pub u32);

/// Identifier of an unknown dimension variable.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DimVarId(pub u32);

/// How much dimensional checking a computation asks for.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum DimensionMode {
    /// No inference, no diagnostics, no overhead.
    Off,
    /// Infer dimensions and report diagnostics; the answer is never changed.
    #[default]
    Advisory,
    /// A dimension conflict in the request or in the answer is an error.
    Strict,
}

/// The dimension of a term, as an owned, printable snapshot.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Dimension {
    /// The value is the polymorphic scalar `0`/`oo` (any dimension).
    pub any: bool,
    /// The value is not a scalar (a predicate, a label, a container).
    pub non_metric: bool,
    /// Concrete base dimensions and the printed exponent.
    pub factors: Vec<(String, String)>,
    /// Unresolved dimension variables, with the printed exponent.
    pub unknown: Vec<(String, String)>,
}

impl Dimension {
    /// The dimensionless value `1`.
    #[must_use]
    pub const fn dimensionless() -> Self {
        Self {
            any: false,
            non_metric: false,
            factors: Vec::new(),
            unknown: Vec::new(),
        }
    }

    /// Whether this is exactly the dimensionless value.
    #[must_use]
    pub const fn is_dimensionless(&self) -> bool {
        !self.any && !self.non_metric && self.factors.is_empty() && self.unknown.is_empty()
    }

    /// Whether no dimension information is available.
    #[must_use]
    pub const fn is_unknown(&self) -> bool {
        !self.any && !self.non_metric && self.factors.is_empty() && !self.unknown.is_empty()
    }
}

impl fmt::Display for Dimension {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        if self.any {
            return f.write_str("any");
        }
        if self.non_metric {
            return f.write_str("non-metric");
        }
        if self.factors.is_empty() && self.unknown.is_empty() {
            return f.write_str("1");
        }
        let mut out = String::new();
        let mut numerator: Vec<String> = Vec::new();
        let mut denominator: Vec<String> = Vec::new();
        for (name, exponent) in &self.factors {
            let (target, value) = if exponent.starts_with('-') {
                (&mut denominator, exponent.trim_start_matches('-'))
            } else {
                (&mut numerator, exponent.as_str())
            };
            target.push(render_power(name, value));
        }
        for (name, exponent) in &self.unknown {
            let (target, value) = if exponent.starts_with('-') {
                (&mut denominator, exponent.trim_start_matches('-'))
            } else {
                (&mut numerator, exponent.as_str())
            };
            target.push(render_power(&format!("?{name}"), value));
        }
        if numerator.is_empty() {
            out.push('1');
        } else {
            out.push_str(&numerator.join("*"));
        }
        if !denominator.is_empty() {
            out.push('/');
            if denominator.len() > 1 {
                out.push('(');
                out.push_str(&denominator.join("*"));
                out.push(')');
            } else if let Some(first) = denominator.first() {
                out.push_str(first);
            }
        }
        f.write_str(&out)
    }
}

fn render_power(
    name: &str,
    exponent: &str,
) -> String {
    match exponent {
        | "1" | "" => name.to_string(),
        | "1/2" => format!("sqrt({name})"),
        | other if is_simple_exponent(other) => format!("{name}^{other}"),
        | other if other.starts_with('(') => format!("{name}^{other}"),
        | other => format!("{name}^({other})"),
    }
}

/// Whether an exponent can be printed without parentheses.
fn is_simple_exponent(text: &str) -> bool {
    !text.is_empty()
        && !text.contains(' ')
        && !text.contains('+')
        && !text.contains('*')
        && text.chars().all(|c| {
            c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-' | '/' | '?' | '(' | ')')
        })
}

/// An error found while inferring dimensions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DimError {
    /// Two dimensions that must agree do not.
    Mismatch {
        /// The dimension on the left.
        left: String,
        /// The dimension on the right.
        right: String,
        /// The expression where the mismatch was found.
        context: String,
    },
    /// An argument that must be dimensionless has a dimension.
    NonDimensionless {
        /// The operator.
        op: String,
        /// The offending argument (printed).
        argument: String,
        /// Its dimension.
        found: String,
    },
    /// A symbolic exponent would have to satisfy an equation to make the
    /// dimensions agree.
    SymbolicConstraint {
        /// The constraint, for example `n = 1`.
        constraint: String,
        /// The expression where it arose.
        context: String,
    },
    /// A unit name was not understood.
    UnknownUnit {
        /// The text.
        unit: String,
    },
    /// An affine quantity (a temperature point) was used where a vector is
    /// required, or two points were combined.
    Affine {
        /// The operation.
        operation: String,
    },
    /// A list mixed quantities of different dimensions.
    MixedList {
        /// The expression.
        context: String,
    },
    /// An operator has no dimension rule yet (a coverage gap).
    Unannotated {
        /// The operator name.
        op: String,
    },
    /// A dimension text could not be parsed.
    Parse {
        /// The text.
        text: String,
        /// What was wrong.
        reason: String,
    },
}

impl DimError {
    /// Whether this error is a real conflict that strict mode rejects. A
    /// missing operator rule is a coverage gap, not a conflict.
    #[must_use]
    pub const fn is_blocking(&self) -> bool {
        !matches!(self, Self::Unannotated { .. } | Self::MixedList { .. })
    }
}

impl fmt::Display for DimError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Mismatch { left, right, context } => {
                write!(f, "dimension mismatch in `{context}`: {left} vs {right}")
            },
            | Self::NonDimensionless { op, argument, found } => write!(
                f,
                "`{op}` requires a dimensionless argument, found {found} in `{argument}`"
            ),
            | Self::SymbolicConstraint { constraint, context } => write!(
                f,
                "dimensions force the constraint `{constraint}` in `{context}`"
            ),
            | Self::UnknownUnit { unit } => write!(f, "unknown unit `{unit}`"),
            | Self::Affine { operation } => {
                write!(f, "affine quantity misuse in `{operation}`")
            },
            | Self::MixedList { context } => {
                write!(f, "list with mixed dimensions in `{context}`")
            },
            | Self::Unannotated { op } => {
                write!(f, "operator `{op}` has no dimension rule")
            },
            | Self::Parse { text, reason } => {
                write!(f, "cannot parse dimension `{text}`: {reason}")
            },
        }
    }
}

/// The result of inferring the dimension of a term.
#[derive(Clone, Debug)]
pub struct DimOutcome {
    /// The inferred dimension.
    pub dimension: Dimension,
    /// Conflicts and gaps found during inference.
    pub errors: Vec<DimError>,
}

/// Information the unit layer supplies about a unit expression.
#[derive(Clone, Debug)]
pub struct UnitInfo {
    /// The dimension of the unit.
    pub dimension: DimTerm,
    /// The scale to coherent SI units.
    pub scale: f64,
    /// Whether this is an affine (offset) unit such as `degC`.
    pub affine: bool,
}

/// Resolves a unit term (the second argument of `quantity`) to its dimension.
///
/// The `units` rule set installs one on the `quantity` operator; the
/// dimension layer never depends on `rules`.
#[derive(Copy, Clone)]
pub struct UnitResolver(pub fn(&Graph, NodeId, &DimArena) -> Option<UnitInfo>);

/// Resolves a unit written as text (`"m/s"`, `"km/h"`) to its dimension.
#[derive(Copy, Clone)]
pub struct UnitTextResolver(pub fn(&Graph, &str, &DimArena) -> Option<UnitInfo>);

/// Context passed to [`DimRule::Custom`].
pub struct DimCtx<'a> {
    /// The graph being analysed.
    pub graph: &'a Graph,
    /// The interned bases and atoms.
    pub arena: &'a DimArena,
    /// The node the rule is applied to.
    pub node: NodeId,
}

impl DimCtx<'_> {
    /// A dimensionless value.
    #[must_use]
    pub const fn dimensionless(&self) -> DimTerm {
        DimTerm::dimensionless()
    }

    /// The polymorphic value.
    #[must_use]
    pub const fn any(&self) -> DimTerm {
        DimTerm::any()
    }

    /// A base dimension with a rational exponent.
    #[must_use]
    pub fn base(
        &self,
        name: &str,
        exponent: BigRational,
    ) -> Option<DimTerm> {
        let id = self.arena.base_id(name)?;
        Some(DimTerm::base(id, exponent))
    }
}

/// How an operator's result dimension is computed from its arguments.
#[derive(Copy, Clone)]
pub enum DimRule {
    /// No rule yet (a coverage gap).
    Unspecified,
    /// Not a scalar: booleans, labels, containers.
    NonMetric,
    /// Explicit bypass: the value is treated as dimensionally polymorphic.
    Any,
    /// The result has the dimension of argument `i`.
    Same(usize),
    /// All arguments must agree; the result is their common dimension.
    Sum,
    /// Product of the argument dimensions.
    Product,
    /// Argument 0 divided by argument 1.
    Quotient,
    /// `a^b`, with `b` dimensionless.
    Power,
    /// The `n`-th root of argument 0.
    Root(u32),
    /// All arguments are dimensionless and so is the result.
    Dimensionless,
    /// `1 / dim(argument 0)`.
    Inverse,
    /// `dim(argument 0)^2`.
    Square,
    /// A list: the common dimension of its scalar leaves.
    List,
    /// A comparison: arguments must agree, result non-metric.
    Compare,
    /// Arguments must agree; the result is a dimensionless angle.
    Angle,
    /// An unknown function application: a fresh dimension variable.
    UnknownFunction,
    /// `diff(f, x)`: `dim(f) / dim(x)`.
    Diff {
        /// Index of the function.
        f: usize,
        /// Index of the variable.
        x: usize,
    },
    /// `diffn(f, x, n)`: `dim(f) / dim(x)^n`.
    DiffN {
        /// Index of the function.
        f: usize,
        /// Index of the variable.
        x: usize,
        /// Index of the order.
        n: usize,
    },
    /// `integral(f, x)` / `defint(f, x, a, b)`: `dim(f) * dim(x)`.
    Integral {
        /// Index of the integrand.
        f: usize,
        /// Index of the variable.
        x: usize,
    },
    /// `sum(f, k, a, b)`: the dimension of `f`.
    SumOver {
        /// Index of the summand.
        f: usize,
        /// Index of the index variable.
        index: usize,
    },
    /// `product(f, k, a, b)`: `dim(f)^(b - a + 1)`.
    ProductOver {
        /// Index of the factor.
        f: usize,
        /// Index of the index variable.
        index: usize,
    },
    /// `quantity(v, u)`: the unit layer supplies the dimension.
    Quantity,
    /// A domain-specific rule.
    Custom(fn(&DimCtx<'_>, &[DimTerm]) -> Result<DimTerm, DimError>),
}

impl fmt::Debug for DimRule {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let name = match self {
            | Self::Unspecified => "Unspecified",
            | Self::NonMetric => "NonMetric",
            | Self::Any => "Any",
            | Self::Same(_) => "Same",
            | Self::Sum => "Sum",
            | Self::Product => "Product",
            | Self::Quotient => "Quotient",
            | Self::Power => "Power",
            | Self::Root(_) => "Root",
            | Self::Dimensionless => "Dimensionless",
            | Self::Inverse => "Inverse",
            | Self::Square => "Square",
            | Self::List => "List",
            | Self::Compare => "Compare",
            | Self::Angle => "Angle",
            | Self::UnknownFunction => "UnknownFunction",
            | Self::Diff { .. } => "Diff",
            | Self::DiffN { .. } => "DiffN",
            | Self::Integral { .. } => "Integral",
            | Self::SumOver { .. } => "SumOver",
            | Self::ProductOver { .. } => "ProductOver",
            | Self::Quantity => "Quantity",
            | Self::Custom(_) => "Custom",
        };
        f.write_str(name)
    }
}

/// Interned base dimensions, exponent atoms and dimension variables.
#[derive(Clone, Debug)]
pub struct DimArena {
    bases: Vec<String>,
    base_ids: HashMap<String, BaseId>,
    atoms: Vec<String>,
    atom_ids: HashMap<String, AtomId>,
    vars: Vec<String>,
    var_ids: HashMap<String, DimVarId>,
}

/// The seven SI base dimensions, in canonical order.
pub const SI_BASES: [&str; 7] = [
    "Length",
    "Mass",
    "Time",
    "Current",
    "Temperature",
    "Amount",
    "Luminosity",
];

impl Default for DimArena {
    fn default() -> Self {
        Self::new()
    }
}

impl DimArena {
    /// An arena containing the seven SI base dimensions.
    #[must_use]
    pub fn new() -> Self {
        let mut arena = Self {
            bases: Vec::new(),
            base_ids: HashMap::new(),
            atoms: Vec::new(),
            atom_ids: HashMap::new(),
            vars: Vec::new(),
            var_ids: HashMap::new(),
        };
        for name in SI_BASES {
            arena.define_base(name);
        }
        arena
    }

    /// Registers (or returns) a base dimension. Base names are case
    /// sensitive; aliases are handled by [`Graph::dimension_from_text`].
    pub fn define_base(
        &mut self,
        name: &str,
    ) -> BaseId {
        if let Some(&id) = self.base_ids.get(name) {
            return id;
        }
        let id = BaseId(u32::try_from(self.bases.len()).unwrap_or(u32::MAX));
        self.bases.push(name.to_string());
        self.base_ids.insert(name.to_string(), id);
        id
    }

    /// The id of a registered base dimension.
    #[must_use]
    pub fn base_id(
        &self,
        name: &str,
    ) -> Option<BaseId> {
        self.base_ids.get(name).copied()
    }

    /// A base dimension by id.
    #[must_use]
    pub fn base_name(
        &self,
        id: BaseId,
    ) -> &str {
        self.bases.get(id.0 as usize).map_or("?", String::as_str)
    }

    /// All known base names, in registration order.
    #[must_use]
    pub fn bases(&self) -> &[String] {
        &self.bases
    }

    pub(crate) fn intern_atom(
        &mut self,
        label: &str,
    ) -> AtomId {
        if let Some(&id) = self.atom_ids.get(label) {
            return id;
        }
        let id = AtomId(u32::try_from(self.atoms.len()).unwrap_or(u32::MAX));
        self.atoms.push(label.to_string());
        self.atom_ids.insert(label.to_string(), id);
        id
    }

    pub(crate) fn intern_var(
        &mut self,
        label: &str,
    ) -> DimVarId {
        if let Some(&id) = self.var_ids.get(label) {
            return id;
        }
        let id = DimVarId(u32::try_from(self.vars.len()).unwrap_or(u32::MAX));
        self.vars.push(label.to_string());
        self.var_ids.insert(label.to_string(), id);
        id
    }

    pub(crate) fn var_name(
        &self,
        id: DimVarId,
    ) -> &str {
        self.vars.get(id.0 as usize).map_or("?", String::as_str)
    }

    pub(crate) fn atoms(&self) -> &[String] {
        &self.atoms
    }
}

/// A dimension expression: a product of base dimensions and dimension
/// variables, or one of the two special values.
#[derive(Clone, Debug)]
pub struct DimTerm {
    pub(crate) repr: DimRepr,
}

#[derive(Clone, Debug)]
pub(crate) enum DimRepr {
    Any,
    NonMetric,
    Scalar {
        bases: Vec<(BaseId, Expo)>,
        vars: Vec<(DimVarId, Expo)>,
    },
}

impl PartialEq for DimTerm {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        match (&self.repr, &other.repr) {
            | (DimRepr::Any, DimRepr::Any) | (DimRepr::NonMetric, DimRepr::NonMetric) => true,
            | (
                DimRepr::Scalar { bases: ab, vars: av },
                DimRepr::Scalar { bases: bb, vars: bv },
            ) => {
                ab.len() == bb.len()
                    && av.len() == bv.len()
                    && ab
                        .iter()
                        .zip(bb)
                        .all(|((ai, ae), (bi, be))| ai == bi && ae.equals(be))
                    && av
                        .iter()
                        .zip(bv)
                        .all(|((ai, ae), (bi, be))| ai == bi && ae.equals(be))
            },
            | _ => false,
        }
    }
}

impl Eq for DimTerm {}

impl DimTerm {
    /// The polymorphic scalar (`0`, `oo`).
    #[must_use]
    pub const fn any() -> Self {
        Self { repr: DimRepr::Any }
    }

    /// A non-scalar value.
    #[must_use]
    pub const fn non_metric() -> Self {
        Self {
            repr: DimRepr::NonMetric,
        }
    }

    /// The dimensionless value.
    #[must_use]
    pub const fn dimensionless() -> Self {
        Self {
            repr: DimRepr::Scalar {
                bases: Vec::new(),
                vars: Vec::new(),
            },
        }
    }

    /// A base dimension with an exponent.
    #[must_use]
    pub fn base(
        id: BaseId,
        exponent: BigRational,
    ) -> Self {
        if exponent.is_zero() {
            return Self::dimensionless();
        }
        Self {
            repr: DimRepr::Scalar {
                bases: vec![(id, Expo::from_rat(exponent))],
                vars: Vec::new(),
            },
        }
    }

    pub(crate) fn var(id: DimVarId) -> Self {
        Self {
            repr: DimRepr::Scalar {
                bases: Vec::new(),
                vars: vec![(id, Expo::one())],
            },
        }
    }

    /// Whether this is the polymorphic scalar.
    #[must_use]
    pub const fn is_any(&self) -> bool {
        matches!(self.repr, DimRepr::Any)
    }

    /// Whether this is a non-scalar value.
    #[must_use]
    pub const fn is_non_metric(&self) -> bool {
        matches!(self.repr, DimRepr::NonMetric)
    }

    /// Whether this is exactly dimensionless.
    #[must_use]
    pub const fn is_dimensionless(&self) -> bool {
        match &self.repr {
            | DimRepr::Scalar { bases, vars } => bases.is_empty() && vars.is_empty(),
            | _ => false,
        }
    }

    /// Whether the value contains unresolved dimension variables.
    #[must_use]
    pub const fn has_unknown(&self) -> bool {
        match &self.repr {
            | DimRepr::Scalar { vars, .. } => !vars.is_empty(),
            | _ => false,
        }
    }

    pub(crate) fn scalar(
        bases: Vec<(BaseId, Expo)>,
        vars: Vec<(DimVarId, Expo)>,
    ) -> Self {
        let mut bases = bases;
        let mut vars = vars;
        bases.retain(|(_, e)| !e.is_zero());
        vars.retain(|(_, e)| !e.is_zero());
        bases.sort_by_key(|(id, _)| *id);
        vars.sort_by_key(|(id, _)| *id);
        if bases.is_empty() && vars.is_empty() {
            return Self::dimensionless();
        }
        Self {
            repr: DimRepr::Scalar { bases, vars },
        }
    }

    /// The product of two dimensions.
    #[must_use]
    pub fn mul(
        &self,
        other: &Self,
    ) -> Self {
        match (&self.repr, &other.repr) {
            | (DimRepr::Any, _) | (_, DimRepr::Any) => Self::any(),
            | (DimRepr::NonMetric, _) | (_, DimRepr::NonMetric) => Self::non_metric(),
            | (
                DimRepr::Scalar { bases: ab, vars: av },
                DimRepr::Scalar { bases: bb, vars: bv },
            ) => {
                let mut bases = ab.clone();
                for (id, exponent) in bb {
                    add_factor(&mut bases, *id, exponent);
                }
                let mut vars = av.clone();
                for (id, exponent) in bv {
                    add_factor_var(&mut vars, *id, exponent);
                }
                Self::scalar(bases, vars)
            },
        }
    }

    /// The inverse dimension.
    #[must_use]
    pub fn inverse(&self) -> Self {
        match &self.repr {
            | DimRepr::Scalar { bases, vars } => {
                let bases = bases
                    .iter()
                    .filter_map(|(id, exponent)| {
                        let negated = exponent.neg();
                        (!negated.is_zero()).then_some((*id, negated))
                    })
                    .collect();
                let vars = vars
                    .iter()
                    .filter_map(|(id, exponent)| {
                        let negated = exponent.neg();
                        (!negated.is_zero()).then_some((*id, negated))
                    })
                    .collect();
                Self::scalar(bases, vars)
            },
            | _ => self.clone(),
        }
    }

    /// `self^exponent` for a rational exponent.
    #[must_use]
    pub fn pow(
        &self,
        exponent: &BigRational,
    ) -> Self {
        self.pow_expo(&Expo::from_rat(exponent.clone()))
    }

    pub(crate) fn pow_expo(
        &self,
        exponent: &Expo,
    ) -> Self {
        match &self.repr {
            | DimRepr::Scalar { bases, vars } => {
                if exponent.is_zero() {
                    return Self::dimensionless();
                }
                let bases = bases.iter().map(|(id, e)| (*id, e.mul(exponent))).collect();
                let vars = vars.iter().map(|(id, e)| (*id, e.mul(exponent))).collect();
                Self::scalar(bases, vars)
            },
            | _ => self.clone(),
        }
    }

    /// A readable description, for diagnostics.
    #[must_use]
    pub fn describe(
        &self,
        arena: &DimArena,
    ) -> String {
        self.to_dimension(arena).to_string()
    }

    /// The owned public snapshot.
    #[must_use]
    pub fn to_dimension(
        &self,
        arena: &DimArena,
    ) -> Dimension {
        match &self.repr {
            | DimRepr::Any => Dimension {
                any: true,
                ..Dimension::default()
            },
            | DimRepr::NonMetric => Dimension {
                non_metric: true,
                ..Dimension::default()
            },
            | DimRepr::Scalar { bases, vars } => {
                let mut factors = Vec::with_capacity(bases.len());
                for (id, exponent) in bases {
                    factors.push((
                        arena.base_name(*id).to_string(),
                        exponent.display(arena.atoms()),
                    ));
                }
                factors.sort_by(|a, b| a.0.cmp(&b.0));
                let mut unknown = Vec::with_capacity(vars.len());
                for (id, exponent) in vars {
                    unknown.push((
                        arena.var_name(*id).to_string(),
                        exponent.display(arena.atoms()),
                    ));
                }
                unknown.sort_by(|a, b| a.0.cmp(&b.0));
                Dimension {
                    any: false,
                    non_metric: false,
                    factors,
                    unknown,
                }
            },
        }
    }
}

fn add_factor(
    factors: &mut Vec<(BaseId, Expo)>,
    id: BaseId,
    exponent: &Expo,
) {
    if let Some(index) = factors.iter().position(|(other, _)| *other == id)
        && let Some(slot) = factors.get_mut(index)
    {
        slot.1 = slot.1.add(exponent);
    } else {
        factors.push((id, exponent.clone()));
    }
}

fn add_factor_var(
    factors: &mut Vec<(DimVarId, Expo)>,
    id: DimVarId,
    exponent: &Expo,
) {
    if let Some(index) = factors.iter().position(|(other, _)| *other == id)
        && let Some(slot) = factors.get_mut(index)
    {
        slot.1 = slot.1.add(exponent);
    } else {
        factors.push((id, exponent.clone()));
    }
}

impl Graph {
    /// The dimension arena of this graph (base dimensions and atoms).
    pub fn dim_arena(&self) -> std::sync::RwLockReadGuard<'_, DimArena> {
        match self.dims.read() {
            | Ok(guard) => guard,
            | Err(poisoned) => poisoned.into_inner(),
        }
    }

    /// Registers a base dimension for this graph.
    pub fn define_base_dimension(
        &mut self,
        name: &str,
    ) {
        match self.dims.write() {
            | Ok(mut arena) => {
                arena.define_base(name);
            },
            | Err(poisoned) => {
                poisoned.into_inner().define_base(name);
            },
        }
    }

    /// Declares the dimension of a symbol from text (`"m/s"`,
    /// `"Length/Time"`, `"any"`).
    ///
    /// # Errors
    /// Returns [`DimError`] if the text cannot be parsed.
    pub fn assume_dimension(
        &mut self,
        symbol: SymbolId,
        text: &str,
    ) -> Result<(), DimError> {
        let term = self.dimension_from_text(text)?;
        self.dim_assumptions.insert(symbol, term);
        Ok(())
    }

    /// Parses a dimension or unit written as text.
    ///
    /// # Errors
    /// Returns [`DimError`] if the text cannot be parsed.
    pub fn dimension_from_text(
        &self,
        text: &str,
    ) -> Result<DimTerm, DimError> {
        let trimmed = text.trim();
        if trimmed.is_empty() || trimmed == "1" || trimmed.eq_ignore_ascii_case("dimensionless") {
            return Ok(DimTerm::dimensionless());
        }
        if trimmed.eq_ignore_ascii_case("any") || trimmed.eq_ignore_ascii_case("scalar") {
            return Ok(DimTerm::any());
        }
        if trimmed.eq_ignore_ascii_case("nonmetric")
            || trimmed.eq_ignore_ascii_case("non-metric")
            || trimmed.eq_ignore_ascii_case("boolean")
        {
            return Ok(DimTerm::non_metric());
        }
        let arena = self.dim_arena();
        // The unit layer knows prefixes, derived units and lengths.
        if let Some(quantity) = self.ops().lookup("quantity")
            && let Some(resolver) = self.ops().attr::<UnitTextResolver>(quantity)
            && let Some(info) = (resolver.0)(self, trimmed, &arena)
        {
            return Ok(info.dimension);
        }
        parse_base_text(&arena, trimmed)
    }

    /// Infers the dimension of `node` and collects conflicts.
    #[must_use]
    pub fn dimension(
        &self,
        node: NodeId,
    ) -> DimOutcome {
        let mut inference = DimInfer::new(self);
        let term = inference.run(node);
        let arena = self.dim_arena();
        DimOutcome {
            dimension: term.to_dimension(&arena),
            errors: inference.into_errors(),
        }
    }

    /// The dimension rule attached to `op`, if any.
    #[must_use]
    pub fn dim_rule(
        &self,
        op: crate::graph::OpId,
    ) -> Option<DimRule> {
        self.ops().attr::<DimRule>(op).copied()
    }
}

fn parse_base_text(
    arena: &DimArena,
    text: &str,
) -> Result<DimTerm, DimError> {
    let cleaned = text.replace(' ', "");
    let (numerator, denominator) = match cleaned.split_once('/') {
        | Some((num, den)) => (num, den),
        | None => (cleaned.as_str(), ""),
    };
    let mut term = DimTerm::dimensionless();
    for (part, sign) in [(numerator, 1_i32), (denominator, -1_i32)] {
        if part.is_empty() || part == "1" {
            continue;
        }
        for factor in part.split('*') {
            let (name, exponent) = match factor.split_once('^') {
                | Some((name, exponent)) => {
                    let exponent = parse_exponent(exponent, text)?;
                    (name, exponent)
                },
                | None => (factor, BigRational::one()),
            };
            let base_name = base_alias(name);
            let Some(base) = arena.base_id(base_name) else {
                return Err(DimError::UnknownUnit {
                    unit: text.to_string(),
                });
            };
            let exponent = if sign < 0 {
                -exponent
            } else {
                exponent
            };
            term = term.mul(&DimTerm::base(base, exponent));
        }
    }
    Ok(term)
}

fn parse_exponent(
    text: &str,
    whole: &str,
) -> Result<BigRational, DimError> {
    let (numerator, denominator) = match text.split_once('/') {
        | Some((num, den)) => (num, den),
        | None => (text, "1"),
    };
    let numerator = numerator.parse::<BigInt>().map_err(|_| DimError::Parse {
        text: whole.to_string(),
        reason: format!("`{text}` is not a rational exponent"),
    })?;
    let denominator = denominator.parse::<BigInt>().map_err(|_| DimError::Parse {
        text: whole.to_string(),
        reason: format!("`{text}` is not a rational exponent"),
    })?;
    if denominator.is_zero() {
        return Err(DimError::Parse {
            text: whole.to_string(),
            reason: String::from("zero denominator"),
        });
    }
    Ok(BigRational::new(numerator, denominator))
}

fn base_alias(name: &str) -> &str {
    match name.to_ascii_lowercase().as_str() {
        | "length" | "l" | "m" | "meter" | "metre" => "Length",
        | "mass" | "m_mass" | "kg" | "g" => "Mass",
        | "time" | "t" | "s" | "second" => "Time",
        | "current" | "i" | "a" | "ampere" => "Current",
        | "temperature" | "temp" | "k" | "kelvin" => "Temperature",
        | "amount" | "n" | "mol" | "mole" => "Amount",
        | "luminosity" | "cd" | "candela" => "Luminosity",
        | _ => name,
    }
}

/// The core signature table. Rule sets extend it with `Installer::dim`.
///
/// A name missing here is a coverage gap: the inference engine treats the
/// operator as dimensionally unknown and the coverage scanner reports it.
#[must_use]
pub fn default_rule(name: &str) -> Option<DimRule> {
    let rule = match name {
        | "add" | "sub" | "min" | "max" | "piecewise" => DimRule::Sum,
        | "mul" | "dot" | "cross" | "matmul" | "kron" => DimRule::Product,
        | "div" => DimRule::Quotient,
        | "pow" => DimRule::Power,
        | "sqrt" => DimRule::Root(2),
        | "cbrt" => DimRule::Root(3),
        | "abs" | "re" | "im" | "conj" | "floor" | "ceil" | "round" | "trunc" | "limit"
        | "taylor" | "laurent" | "fourier_series" | "convert" | "unify_expression"
        | "simplify_units" => DimRule::Same(0),
        | "sign" | "arg" | "exp" | "ln" | "log" | "sin" | "cos" | "tan" | "cot" | "sec" | "csc"
        | "asin" | "acos" | "atan" | "acot" | "asec" | "acsc" | "sinh" | "cosh" | "tanh"
        | "coth" | "sech" | "csch" | "asinh" | "acosh" | "atanh" | "acoth" | "asech"
        | "acsch" | "erf" | "erfc" | "erfi" | "erfcinv" | "gamma" | "lgamma" | "digamma"
        | "polygamma" | "zeta" | "hurwitz_zeta" | "beta" | "lambertw" | "sinc" | "besselj"
        | "bessely" | "besseli" | "besselk" | "airy_ai" | "airy_bi" | "airy_aip" | "airy_bip"
        | "pi" | "e" | "true_pi" | "euler_gamma" | "i" | "golden_ratio" => DimRule::Dimensionless,
        | "atan2" => DimRule::Angle,
        | "list" | "tensor" | "vector" | "matrix" => DimRule::List,
        | "eq" | "lt" | "le" | "gt" | "ge" | "ne" => DimRule::Compare,
        | "and" | "or" | "xor" | "not" | "implies" | "iff" | "true" | "false" | "converges"
        | "dimension" => DimRule::NonMetric,
        | "oo" | "infinity" | "big_o" => DimRule::Any,
        | "apply" => DimRule::UnknownFunction,
        | "diff" => DimRule::Diff { f: 0, x: 1 },
        | "diffn" => DimRule::DiffN { f: 0, x: 1, n: 2 },
        | "integral" | "real_integral" | "defint" => DimRule::Integral { f: 0, x: 1 },
        | "sum" => DimRule::SumOver { f: 0, index: 1 },
        | "product" => DimRule::ProductOver { f: 0, index: 1 },
        | "quantity" => DimRule::Quantity,
        | "inv" | "pinv" => DimRule::Inverse,
        | "square" => DimRule::Square,
        // Complex analysis and verification. Branch wrappers, residues and
        // contour results carry the dimension of the function they analyse;
        // branches/curves/singularity sets are not scalars; winding numbers,
        // zero/pole counts and verification verdicts are pure numbers or
        // verdicts.
        | "I" => DimRule::Dimensionless,
        | "log_branch"
        | "sqrt_branch"
        | "root_branch"
        | "power_branch"
        | "asin_branch"
        | "acos_branch"
        | "atan_branch"
        | "asinh_branch"
        | "acosh_branch"
        | "atanh_branch"
        | "continue_along"
        | "mobius_apply"
        | "residue"
        | "contour_integral"
        | "cauchy_integral"
        | "contour_integral_residue_theorem"
        | "distance" => DimRule::Same(0),
        | "cauchy_derivative" => DimRule::Diff { f: 0, x: 1 },
        | "radius_of_convergence" => DimRule::Same(1),
        | "count_zeros_poles" | "winding_number" => DimRule::Dimensionless,
        | "circle"
        | "path"
        | "polygon"
        | "ray"
        | "poles"
        | "zeros_of"
        | "singularity"
        | "branch_cut"
        | "mobius_compose"
        | "mobius_inverse"
        | "is_inside_contour"
        | "check_analytic"
        | "verify_solution"
        | "verify_integral"
        | "verify_definite_integral"
        | "verify_ode_solution"
        | "verify_inverse"
        | "verify_derivative"
        | "verify_limit" => DimRule::NonMetric,
        | _ => return None,
    };
    Some(rule)
}

/// Converts a numeric node to an exact exponent, when it is rational.
pub(crate) fn exact_exponent(number: &crate::graph::Number) -> Option<BigRational> {
    number.to_rational()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponent_algebra_is_a_field() {
        let mut arena = DimArena::new();
        let a = arena.intern_atom("a");
        let b = arena.intern_atom("b");
        let x = Expo::atom(a);
        let y = Expo::atom(b);
        assert!((x.add(&y)).sub(&y).equals(&x));
        assert!((x.mul(&y)).equals(&y.mul(&x)));
        assert!((x.mul(&x.inverse().unwrap_or_else(Expo::one))).is_one());
        assert!((x.pow_i32(2).unwrap_or_else(Expo::one)).equals(&x.mul(&x)));
    }

    #[test]
    fn symbolic_exponents_display() {
        let mut arena = DimArena::new();
        let n = arena.intern_atom("n");
        let length = arena.base_id("Length").unwrap_or(BaseId(0));
        let term = DimTerm::base(length, BigRational::one()).pow_expo(&Expo::atom(n));
        assert_eq!(term.describe(&arena), "Length^n");
    }

    #[test]
    fn unification_solves_dimension_variables() {
        let mut arena = DimArena::new();
        let length = DimTerm::base(
            arena.base_id("Length").unwrap_or(BaseId(0)),
            BigRational::one(),
        );
        let mut unifier = infer::Unifier::default();
        let variable = DimTerm::var(arena.intern_var("x"));
        let resolved = unifier
            .unify(&variable, &length)
            .unwrap_or_else(|_| panic!("unify"));
        assert!(!resolved.has_unknown());
        assert_eq!(resolved.describe(&arena), "Length");
        assert!(unifier.unify(&resolved, &length).is_ok());
    }

    #[test]
    fn mismatch_is_detected() {
        let arena = DimArena::new();
        let length = DimTerm::base(
            arena.base_id("Length").unwrap_or(BaseId(0)),
            BigRational::one(),
        );
        let time = DimTerm::base(
            arena.base_id("Time").unwrap_or(BaseId(0)),
            BigRational::one(),
        );
        let mut unifier = infer::Unifier::default();
        assert!(unifier.unify(&length, &time).is_err());
    }

    #[test]
    fn core_signatures_exist() {
        for name in [
            "add", "mul", "pow", "sqrt", "exp", "ln", "diff", "integral", "sum", "product",
        ] {
            assert!(default_rule(name).is_some(), "{name}");
        }
    }

    #[test]
    fn no_operator_has_an_accidental_rule() {
        let graph = Graph::new();
        for op in graph.ops().iter() {
            let name = &*op.1.name;
            if let Some(rule) = default_rule(name) {
                assert!(
                    !matches!(rule, DimRule::Unspecified),
                    "{name} must not have an unspecified rule"
                );
            }
        }
    }

    fn dimension_of(
        source: &str,
        assumptions: &[(&str, &str)],
    ) -> DimOutcome {
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        for (name, text) in assumptions {
            let symbol = graph.interner_mut().symbol(name);
            graph
                .assume_dimension(symbol, text)
                .unwrap_or_else(|error| panic!("{name}: {error}"));
        }
        let root = graph
            .parse(source)
            .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
        graph.dimension(root)
    }

    #[test]
    fn arithmetic_dimensions() {
        let product = dimension_of("x*y", &[("x", "m"), ("y", "s")]);
        assert_eq!(product.dimension.to_string(), "Length*Time");
        assert!(product.errors.is_empty(), "{:?}", product.errors);

        let square = dimension_of("x^2", &[("x", "m")]);
        assert_eq!(square.dimension.to_string(), "Length^2");

        let root = dimension_of("sqrt(x^2)", &[("x", "m")]);
        assert_eq!(root.dimension.to_string(), "Length");

        let cancellations = dimension_of("(x/y)*(y/x)", &[("x", "m"), ("y", "s")]);
        assert_eq!(cancellations.dimension.to_string(), "1");
    }

    #[test]
    fn symbolic_exponents_in_inference() {
        let result = dimension_of("x^n", &[("x", "m")]);
        assert_eq!(result.dimension.to_string(), "Length^n");
        assert!(result.errors.is_empty(), "{:?}", result.errors);

        let product = dimension_of("product(x, k, 1, n)", &[("x", "m")]);
        assert_eq!(product.dimension.to_string(), "Length^n");
        assert!(product.errors.is_empty(), "{:?}", product.errors);

        let diff = dimension_of("diff(x^n, x)", &[("x", "m")]);
        assert_eq!(diff.dimension.to_string(), "Length^(n - 1)");
    }

    #[test]
    fn dimensionless_arguments_are_enforced() {
        let good = dimension_of("exp(t)", &[("t", "1")]);
        assert_eq!(good.dimension.to_string(), "1");
        assert!(good.errors.is_empty(), "{:?}", good.errors);

        let bad = dimension_of("exp(x)", &[("x", "m")]);
        assert!(
            bad.errors
                .iter()
                .any(|error| matches!(error, DimError::NonDimensionless { .. })),
            "{:?}",
            bad.errors
        );
        let bad = dimension_of("x + y", &[("x", "m"), ("y", "s")]);
        assert!(
            bad.errors
                .iter()
                .any(|error| matches!(error, DimError::Mismatch { .. })),
            "{:?}",
            bad.errors
        );
    }

    #[test]
    fn quantities_and_calculus() {
        let quantity = dimension_of("quantity(3, m/s)", &[]);
        assert_eq!(quantity.dimension.to_string(), "Length/Time");
        assert!(quantity.errors.is_empty(), "{:?}", quantity.errors);

        let derivative = dimension_of("diff(quantity(3, m/s) * x^2, x)", &[("x", "m")]);
        assert_eq!(derivative.dimension.to_string(), "Length^2/Time");

        let integral = dimension_of(
            "defint(quantity(9.8, m/s^2), t, 0, quantity(3, s))",
            &[("t", "s")],
        );
        assert_eq!(integral.dimension.to_string(), "Length/Time");
        assert!(integral.errors.is_empty(), "{:?}", integral.errors);
    }

    #[test]
    fn zero_is_polymorphic() {
        let result = dimension_of("x + 0", &[("x", "m")]);
        assert_eq!(result.dimension.to_string(), "Length");
        assert!(result.errors.is_empty(), "{:?}", result.errors);
        let zero = dimension_of("0", &[]);
        assert_eq!(zero.dimension.to_string(), "any");
    }

    #[test]
    fn user_defined_base_dimensions() {
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        graph.define_base_dimension("Currency");
        let symbol = graph.interner_mut().symbol("price");
        graph
            .assume_dimension(symbol, "Currency")
            .unwrap_or_else(|error| panic!("{error}"));
        let root = graph
            .parse("price^2")
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(graph.dimension(root).dimension.to_string(), "Currency^2");
    }

    #[test]
    fn every_standard_operator_has_a_dimension_rule() {
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        let mut missing = Vec::new();
        for (op, descriptor) in graph.ops().iter() {
            let exempt = matches!(&*descriptor.name, "lit" | "sym");
            if !exempt && graph.ops().attr::<DimRule>(op).is_none() {
                missing.push(descriptor.name.to_string());
            }
        }
        assert!(
            missing.is_empty(),
            "{} operators have no dimension rule: {missing:?}",
            missing.len()
        );
    }
}
