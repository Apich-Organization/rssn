//! Dimensions as a typing layer over the term graph.
//!
//! Every scalar term has a dimension. Base dimensions are interned names
//! (`Length`, `Mass`, `Time`, ... and any user-defined base); exponents are
//! exact rationals or, for symbolic exponents, rational functions over
//! interned atoms (`x^n`, `x^(a+b)`, `x^(1/n)`).
//!
//! Operator signatures are [`DimRule`] attributes. The inference engine in
//! `dim_infer` walks the hash-consed DAG, unifies the dimensions of all
//! members of an e-class and reports conflicts. The public entry points are
//! [`Graph::dimension`] and [`Graph::assume_dimension`]; `Config` and
//! `Answer` expose the result.
//!
//! A finding is a [`DimError`]. Its [`DimError::code`] is a stable,
//! machine-readable kind (`mismatch`, `unknown_unit`, ...) for callers that
//! should not parse the printable message, and [`DimError::is_blocking`]
//! says whether Strict mode turns it into a hard error. Rules that are not
//! data-expressible are [`DimRule::Custom`] functions in this module:
//! `matexp` (a dimensionless product) and the complex-analysis family, which
//! lives on a dimensionless complex plane (`dimensionless_complex_variable`).
//!
//! The design and the remaining limits are documented in
//! `docs/DIMENSIONS.md` (§13 lists the current limitations by priority). The
//! standard rule set's rewrites are dimension-checked by
//! `tests/dimensions.rs::per_rewrite_random_dimension_scan`.

// The exponent field and the inference engine are implementation details of
// this module; `pub(crate)` is already the widest visibility reachable from
// the sibling implementation modules, so the lint's suggested `pub(super)` would force
// re-exports at every use site.
#![allow(clippy::redundant_pub_crate)]

use std::collections::HashMap;
use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;

pub(crate) use crate::graph::dim_expr::AtomId;
pub(crate) use crate::graph::dim_expr::Expo;
pub(crate) use crate::graph::dim_infer::DimInfer;

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

    /// The concrete factors with exact rational exponents.
    ///
    /// The printable form keeps exponents as text so that `Length^n` can be
    /// shown; this accessor returns the exact rational ones and is what a
    /// round trip through [`DimTerm`] uses.
    ///
    /// # Errors
    /// Fails when an exponent is symbolic (`n`, `a + b`), which has no
    /// rational value.
    pub fn exact_factors(&self) -> Result<Vec<(String, BigRational)>, DimError> {
        exact_exponents(&self.factors)
    }

    /// The unresolved variables with exact rational exponents.
    ///
    /// # Errors
    /// See [`Dimension::exact_factors`].
    pub fn exact_unknown(&self) -> Result<Vec<(String, BigRational)>, DimError> {
        exact_exponents(&self.unknown)
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

/// Parses printed exponents back into exact rationals.
fn exact_exponents(factors: &[(String, String)]) -> Result<Vec<(String, BigRational)>, DimError> {
    let mut out = Vec::with_capacity(factors.len());
    for (name, exponent) in factors {
        let value = parse_rational_exponent(exponent).ok_or_else(|| {
            DimError::Parse {
                text: format!("{name}^{exponent}"),
                reason: String::from("the exponent is not an exact rational"),
            }
        })?;
        out.push((name.clone(), value));
    }
    Ok(out)
}

/// Parses an exponent as printed by [`Expo::display`]: `1`, `-2`, `1/2`.
fn parse_rational_exponent(text: &str) -> Option<BigRational> {
    let text = text.trim();
    if text.is_empty() {
        return Some(BigRational::one());
    }
    let (numerator, denominator) = text.split_once('/').map_or((text, "1"), |(n, d)| (n, d));
    let numerator = numerator.parse::<BigInt>().ok()?;
    let denominator = denominator.parse::<BigInt>().ok()?;
    if denominator.is_zero() {
        return None;
    }
    Some(BigRational::new(numerator, denominator))
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
    ///
    /// A text declaration (`Config::dimension("x", "parsec")`) is a hard
    /// parse error: `context` is `None`. A unit inside a term
    /// (`quantity(1, parsec)`) is only an unknown, not a conflict: the value
    /// keeps the dimension `any` and the finding is advisory in both
    /// [`DimensionMode::Advisory`] and [`DimensionMode::Strict`]
    /// (`context` is the printed term).
    UnknownUnit {
        /// The text.
        unit: String,
        /// The term whose unit could not be resolved, when the finding comes
        /// from inference rather than from a text declaration.
        context: Option<String>,
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
    /// Whether this error is a real conflict that strict mode rejects.
    ///
    /// A missing operator rule is a coverage gap; a mixed list is a
    /// structural note about a container; and an unresolved unit inside a
    /// term is an *unknown*, not a conflict (the value is `any`). None of
    /// them turns a Strict computation into a hard error. An unresolved unit
    /// in requested dimension text (`context` is `None`) is a declaration
    /// error and stays blocking.
    #[must_use]
    pub const fn is_blocking(&self) -> bool {
        match self {
            | Self::Unannotated { .. }
            | Self::MixedList { .. }
            | Self::UnknownUnit { context: Some(_), .. } => false,
            | _ => true,
        }
    }

    /// A stable, machine-readable kind for this finding.
    ///
    /// The code is part of the public contract: a binding branches on it
    /// instead of parsing [`Display`]. The list is closed and never changes
    /// for an existing variant; a future variant gets a new code. In
    /// particular, [`Self::UnknownUnit`] has one code for both a text
    /// declaration and an unresolved unit inside a term — the difference is
    /// [`Self::is_blocking`] (`context` is `None` for a declaration).
    ///
    /// | code | variant |
    /// |---|---|
    /// | `mismatch` | `Mismatch` |
    /// | `non_dimensionless` | `NonDimensionless` |
    /// | `symbolic_constraint` | `SymbolicConstraint` |
    /// | `unknown_unit` | `UnknownUnit` |
    /// | `affine` | `Affine` |
    /// | `mixed_list` | `MixedList` |
    /// | `unannotated` | `Unannotated` |
    /// | `parse` | `Parse` |
    ///
    /// This exists because diagnostics are otherwise only printable text
    /// (see `docs/DIMENSIONS.md` §7); a caller that wants to treat a
    /// coverage gap (`unannotated`) differently from a conflict
    /// (`mismatch`) can now do so without string matching. Serializing the
    /// code through the C API is a later round; the method is additive.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            | Self::Mismatch { .. } => "mismatch",
            | Self::NonDimensionless { .. } => "non_dimensionless",
            | Self::SymbolicConstraint { .. } => "symbolic_constraint",
            | Self::UnknownUnit { .. } => "unknown_unit",
            | Self::Affine { .. } => "affine",
            | Self::MixedList { .. } => "mixed_list",
            | Self::Unannotated { .. } => "unannotated",
            | Self::Parse { .. } => "parse",
        }
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
            | Self::NonDimensionless { op, argument, found } => {
                write!(
                    f,
                    "`{op}` requires a dimensionless argument, found {found} in `{argument}`"
                )
            },
            | Self::SymbolicConstraint { constraint, context } => {
                write!(
                    f,
                    "dimensions force the constraint `{constraint}` in `{context}`"
                )
            },
            | Self::UnknownUnit { unit, context } => {
                match context {
                    | Some(context) => {
                        write!(f, "unresolved unit `{unit}` in `{context}`")
                    },
                    | None => write!(f, "unknown unit `{unit}`"),
                }
            },
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
    /// A structured tuple/decomposition. Like [`DimRule::NonMetric`] it is
    /// reported as a container, but it is *compatible* with the dimension of
    /// the reduced literal instead of conflicting with it (the search unions
    /// a request with its answer).
    Container,
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
    ///
    /// Shorthand for [`DimRule::RootOf`] with `index = 0` and `degree = n`.
    Root(u32),
    /// The `degree`-th root of argument `index`.
    ///
    /// [`DimRule::Root`] is hard-wired to argument 0, which is right for
    /// `sqrt(x)` but not for an operator whose root operand is not its first
    /// argument — `sde_wiener(seed, m, steps, h)` has dimension
    /// `sqrt(dim(h))`, the root of argument 3. This form states the operand
    /// explicitly; a `degree` of 0 is treated as 1 rather than dividing by
    /// zero.
    RootOf {
        /// The root degree (`2` is a square root).
        degree: u32,
        /// Index of the operand whose dimension is rooted.
        index: usize,
    },
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
    /// A parameterised distribution family. Shape parameters must be
    /// dimensionless; location and scale parameters must have the dimension
    /// of the value.
    Family {
        /// Where the value's dimension comes from.
        result: FamilySource,
        /// Arguments that must have the value's dimension (a location, a
        /// scale, or another parameter tied to it).
        same: &'static [usize],
        /// Arguments that must be pure numbers (shape parameters, counts,
        /// probabilities).
        dimensionless: &'static [usize],
    },
    /// `matfun(f, x, A)`, the matrix function `f(A)`.
    ///
    /// When `f` is a polynomial in `x` whose terms all have the same total
    /// degree `k`, the value has entries `dim(f) * dim(A)^k`. Otherwise `A`
    /// must be dimensionless (the `exp`/trig case).
    MatrixFunction {
        /// Index of the scalar function `f`.
        function: usize,
        /// Index of its formal variable `x`.
        variable: usize,
        /// Index of the matrix `A`.
        matrix: usize,
    },
    /// `quantity(v, u)`: the unit layer supplies the dimension.
    Quantity,
    /// A domain-specific rule.
    Custom(fn(&DimCtx<'_>, &[DimTerm]) -> Result<DimTerm, DimError>),
}

/// Where a [`DimRule::Family`] value's dimension comes from.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FamilySource {
    /// The dimension of argument `i` (a location or a scale).
    Argument(usize),
    /// The inverse dimension of argument `i` (a rate).
    Inverse(usize),
}

impl fmt::Debug for DimRule {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let name = match self {
            | Self::Unspecified => "Unspecified",
            | Self::NonMetric => "NonMetric",
            | Self::Container => "Container",
            | Self::Any => "Any",
            | Self::Same(_) => "Same",
            | Self::Sum => "Sum",
            | Self::Product => "Product",
            | Self::Quotient => "Quotient",
            | Self::Power => "Power",
            | Self::Root(_) => "Root",
            | Self::RootOf { .. } => "RootOf",
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
            | Self::Family { .. } => "Family",
            | Self::MatrixFunction { .. } => "MatrixFunction",
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

    /// The names of all base dimensions, in id order.
    #[must_use]
    pub fn base_names(&self) -> &[String] {
        &self.bases
    }

    /// Registers a base dimension, returning its id. Same as
    /// [`DimArena::define_base`]; the spelling a host expects.
    pub fn intern(
        &mut self,
        name: &str,
    ) -> BaseId {
        self.define_base(name)
    }

    /// A copy of the arena's interning tables.
    ///
    /// A snapshot lets a host save a numbering, run something that extends
    /// it, and put it back — the arena is append-only, so restoring means
    /// dropping whatever was appended.
    #[must_use]
    pub fn snapshot(&self) -> DimArenaSnapshot {
        DimArenaSnapshot {
            bases: self.bases.clone(),
            atoms: self.atoms.clone(),
            vars: self.vars.clone(),
        }
    }

    /// Restores a snapshot taken by [`DimArena::snapshot`].
    pub fn restore(
        &mut self,
        snapshot: &DimArenaSnapshot,
    ) {
        self.bases.clone_from(&snapshot.bases);
        self.atoms.clone_from(&snapshot.atoms);
        self.vars.clone_from(&snapshot.vars);
        self.base_ids = self
            .bases
            .iter()
            .enumerate()
            .map(|(index, name)| {
                (
                    name.clone(),
                    BaseId(u32::try_from(index).unwrap_or(u32::MAX)),
                )
            })
            .collect();
        self.atom_ids = self
            .atoms
            .iter()
            .enumerate()
            .map(|(index, name)| {
                (
                    name.clone(),
                    AtomId(u32::try_from(index).unwrap_or(u32::MAX)),
                )
            })
            .collect();
        self.var_ids = self
            .vars
            .iter()
            .enumerate()
            .map(|(index, name)| {
                (
                    name.clone(),
                    DimVarId(u32::try_from(index).unwrap_or(u32::MAX)),
                )
            })
            .collect();
    }
}

/// An owned copy of a [`DimArena`]'s interning tables.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DimArenaSnapshot {
    bases: Vec<String>,
    atoms: Vec<String>,
    vars: Vec<String>,
}

impl DimArenaSnapshot {
    /// Number of base dimensions in the snapshot.
    #[must_use]
    pub const fn base_count(&self) -> usize {
        self.bases.len()
    }

    /// The base dimension names in the snapshot.
    #[must_use]
    pub fn base_names(&self) -> &[String] {
        &self.bases
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
    /// A structured container (a tuple or a decomposition result). It is not
    /// a scalar, but unlike [`DimRepr::NonMetric`] it is compatible with the
    /// dimension of its reduced literal: kernels union a request with its
    /// answer, and a tuple's expanded `list(...)` can legitimately have a
    /// scalar dimension (all parts dimensionless) or `Any` (parts with
    /// different dimensions).
    Container,
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
            | (DimRepr::Any, DimRepr::Any)
            | (DimRepr::NonMetric, DimRepr::NonMetric)
            | (DimRepr::Container, DimRepr::Container) => true,
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

    /// A structured container (a tuple or decomposition result).
    #[must_use]
    pub const fn container() -> Self {
        Self {
            repr: DimRepr::Container,
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

    /// A base dimension raised to a symbolic exponent.
    ///
    /// The exponent is a name (`n` in `Length^n`), interned in `arena` so that
    /// two occurrences of the same symbol are the same atom. This is the
    /// constructor a type-level `DimExpr::Pow(_, Exp::Symbol(..))` needs; the
    /// rational case stays on [`DimTerm::base`], which is the common one.
    #[must_use]
    pub fn base_symbolic(
        id: BaseId,
        exponent: &str,
        arena: &mut DimArena,
    ) -> Self {
        let atom = arena.intern_atom(exponent);
        Self::scalar(vec![(id, Expo::atom(atom))], Vec::new())
    }

    /// A dimension variable raised to a symbolic exponent.
    #[must_use]
    pub fn var_symbolic(
        id: DimVarId,
        exponent: &str,
        arena: &mut DimArena,
    ) -> Self {
        let atom = arena.intern_atom(exponent);
        Self::scalar(Vec::new(), vec![(id, Expo::atom(atom))])
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

    /// Whether this is a non-scalar value (including a container).
    #[must_use]
    pub const fn is_non_metric(&self) -> bool {
        matches!(self.repr, DimRepr::NonMetric | DimRepr::Container)
    }

    /// Whether this is a structured container rather than a lone non-scalar.
    #[must_use]
    pub const fn is_container(&self) -> bool {
        matches!(self.repr, DimRepr::Container)
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
            // A container has no multiplicative structure either; combining
            // it with anything is a non-scalar value.
            | (DimRepr::NonMetric | DimRepr::Container, _)
            | (_, DimRepr::NonMetric | DimRepr::Container) => Self::non_metric(),
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

    /// The concrete factors as `(base name, exact exponent)`, sorted by name.
    ///
    /// Only factors with a rational exponent appear; the rest are reported
    /// by [`DimTerm::symbolic_factors`]. This is the structured accessor a
    /// host uses instead of parsing the printed form.
    #[must_use]
    pub fn factors(
        &self,
        arena: &DimArena,
    ) -> Vec<(String, BigRational)> {
        match &self.repr {
            | DimRepr::Scalar { bases, .. } => {
                bases
                    .iter()
                    .filter_map(|(id, exponent)| {
                        exponent
                            .is_constant()
                            .map(|e| (arena.base_name(*id).to_owned(), e))
                    })
                    .collect()
            },
            | _ => Vec::new(),
        }
    }

    /// The dimension variables as `(name, exact exponent)`, sorted by name.
    #[must_use]
    pub fn vars(
        &self,
        arena: &DimArena,
    ) -> Vec<(String, BigRational)> {
        match &self.repr {
            | DimRepr::Scalar { vars, .. } => {
                vars.iter()
                    .filter_map(|(id, exponent)| {
                        exponent
                            .is_constant()
                            .map(|e| (arena.var_name(*id).to_owned(), e))
                    })
                    .collect()
            },
            | _ => Vec::new(),
        }
    }

    /// The factors whose exponent is not a rational, as printed.
    #[must_use]
    pub fn symbolic_factors(
        &self,
        arena: &DimArena,
    ) -> Vec<(String, String)> {
        match &self.repr {
            | DimRepr::Scalar { bases, vars } => {
                bases
                    .iter()
                    .filter(|(_, e)| e.is_constant().is_none())
                    .map(|(id, e)| (arena.base_name(*id).to_owned(), e.display(arena.atoms())))
                    .chain(
                        vars.iter()
                            .filter(|(_, e)| e.is_constant().is_none())
                            .map(|(id, e)| {
                                (arena.var_name(*id).to_owned(), e.display(arena.atoms()))
                            }),
                    )
                    .collect()
            },
            | _ => Vec::new(),
        }
    }

    /// A dimension variable raised to an exact power.
    #[must_use]
    pub fn var_pow(
        id: DimVarId,
        exponent: BigRational,
    ) -> Self {
        if exponent.is_zero() {
            return Self::dimensionless();
        }
        Self::scalar(Vec::new(), vec![(id, Expo::from_rat(exponent))])
    }

    /// The dimension of a printed [`Dimension`], in this arena.
    ///
    /// The inverse of [`DimTerm::to_dimension`] for dimensions whose
    /// exponents are rational.
    ///
    /// # Errors
    /// Fails when an exponent is symbolic or a base dimension is unknown to
    /// the arena.
    pub fn from_dimension(
        dimension: &Dimension,
        arena: &DimArena,
    ) -> Result<Self, DimError> {
        if dimension.any {
            return Ok(Self::any());
        }
        if dimension.non_metric {
            return Ok(Self::non_metric());
        }
        let mut bases = Vec::new();
        for (name, exponent) in dimension.exact_factors()? {
            let id = arena.base_id(&name).ok_or_else(|| {
                DimError::Parse {
                    text: name.clone(),
                    reason: String::from("unknown base dimension"),
                }
            })?;
            bases.push((id, Expo::from_rat(exponent)));
        }
        let mut vars = Vec::new();
        for (name, exponent) in dimension.exact_unknown()? {
            let id = arena.var_ids.get(&name).copied().ok_or_else(|| {
                DimError::Parse {
                    text: name.clone(),
                    reason: String::from("unknown dimension variable"),
                }
            })?;
            vars.push((id, Expo::from_rat(exponent)));
        }
        Ok(Self::scalar(bases, vars))
    }

    /// The dimension as text, using `arena` for the names.
    #[must_use]
    pub fn to_string(
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
            | DimRepr::Any => {
                Dimension {
                    any: true,
                    ..Dimension::default()
                }
            },
            | DimRepr::NonMetric | DimRepr::Container => {
                Dimension {
                    non_metric: true,
                    ..Dimension::default()
                }
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
    ///
    /// The arena is shared with every graph cloned from this one until one of
    /// them writes; see [`Graph::with_dim_arena_mut`].
    pub fn dim_arena(&self) -> std::sync::RwLockReadGuard<'_, DimArena> {
        match self.dims.read() {
            | Ok(guard) => guard,
            | Err(poisoned) => poisoned.into_inner(),
        }
    }

    /// Runs `body` with mutable access to the dimension arena.
    ///
    /// Interning a base dimension is the only mutation the arena has, and it
    /// happens whenever a structured dimension names something new; a caller
    /// that has a [`DimExpr`](crate::kernel::signature::DimExpr) needs the
    /// guard rather than the read-only [`Graph::dim_arena`].
    ///
    /// The arena is copy-on-write. The first mutation on a cloned graph gives
    /// it an arena of its own, so a base defined during one request cannot
    /// reach the template the request was cloned from — and a request that
    /// defines nothing shares the arena and pays nothing.
    pub fn with_dim_arena_mut<R>(
        &mut self,
        body: impl FnOnce(&mut DimArena) -> R,
    ) -> R {
        // `Arc::make_mut` cannot be used directly: `RwLock<T>` is not `Clone`.
        // Detach by hand when the arena is shared, which is exactly the
        // copy-on-write `Arc::make_mut` would perform.
        if std::sync::Arc::strong_count(&self.dims) > 1 {
            let detached = match self.dims.read() {
                | Ok(guard) => guard.clone(),
                | Err(poisoned) => poisoned.into_inner().clone(),
            };
            self.dims = std::sync::Arc::new(std::sync::RwLock::new(detached));
        }
        match self.dims.write() {
            | Ok(mut arena) => body(&mut arena),
            | Err(poisoned) => body(&mut poisoned.into_inner()),
        }
    }

    /// Registers a base dimension for this graph.
    ///
    /// Copy-on-write: defining a base on a cloned graph detaches the arena
    /// first, so the template keeps its own numbering.
    pub fn define_base_dimension(
        &mut self,
        name: &str,
    ) {
        self.with_dim_arena_mut(|arena| {
            arena.define_base(name);
        });
    }

    /// Declares the dimension of a symbol from a structured term.
    ///
    /// This is the entry a typed front end should use: it never parses text,
    /// so a unit symbol and a base dimension with the same name cannot be
    /// confused. Base names the arena does not know are defined on the way
    /// in, which is how a host introduces a dimension of its own.
    pub fn assume_dimension_value(
        &mut self,
        symbol: SymbolId,
        term: DimTerm,
    ) {
        self.dim_assumptions.insert(symbol, term);
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
    /// # Precedence
    /// 1. `base:Name` is always a base dimension: `Name` is first looked up
    ///    exactly (so a custom base called `T` is reachable even though `T`
    ///    is the tesla), then through the base aliases (`base:Time`).
    /// 2. Otherwise a unit expression is tried first (`"m/s"`, `"km/h"`,
    ///    `"T"` is tesla).
    /// 3. Only when no unit matches is the text read as a product of base
    ///    dimensions (`"Length/Time"`, `"Currency"`).
    ///
    /// A standard unit is therefore never silently reinterpreted as a base
    /// dimension. A custom base whose name is also a unit symbol must be
    /// written with the `base:` prefix.
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
        // `base:Name` forces the base dimension, bypassing unit symbols.
        if let Some(name) = trimmed.strip_prefix("base:") {
            return parse_base_text(&arena, name.trim(), true);
        }
        // The unit layer knows prefixes, derived units and lengths.
        if let Some(quantity) = self.ops().lookup("quantity")
            && let Some(resolver) = self.ops().attr::<UnitTextResolver>(quantity)
            && let Some(info) = (resolver.0)(self, trimmed, &arena)
        {
            return Ok(info.dimension);
        }
        parse_base_text(&arena, trimmed, false)
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

/// Split a dimension/unit product at its first `/` that is not inside
/// parentheses, so that a rational exponent (`m^(1/2)`) is one factor.
pub(crate) fn split_unit_fraction(text: &str) -> (&str, &str) {
    let mut depth = 0_i32;
    for (index, character) in text.char_indices() {
        match character {
            | '(' => depth = depth.saturating_add(1),
            | ')' => depth = depth.saturating_sub(1),
            | '/' if depth == 0 => return (&text[..index], &text[index + 1..]),
            | _ => {},
        }
    }
    (text, "")
}

fn parse_base_text(
    arena: &DimArena,
    text: &str,
    exact: bool,
) -> Result<DimTerm, DimError> {
    let cleaned = text.replace(' ', "");
    let (numerator, denominator) = split_unit_fraction(&cleaned);
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
            let Some(base) = base_lookup(arena, name, exact) else {
                return Err(DimError::UnknownUnit {
                    unit: text.to_string(),
                    context: None,
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

/// The base dimension a factor names.
///
/// With `exact` (the `base:Name` syntax) only the registered base name is
/// matched, so `base:T` is a custom base `T` (or an error), never the alias
/// for `Time`. Without it the friendly aliases are accepted as well.
fn base_lookup(
    arena: &DimArena,
    name: &str,
    exact: bool,
) -> Option<BaseId> {
    if exact {
        return arena.base_id(name);
    }
    arena.base_id(base_alias(name))
}

fn parse_exponent(
    text: &str,
    whole: &str,
) -> Result<BigRational, DimError> {
    let text = text.trim();
    let text = text
        .strip_prefix('(')
        .and_then(|rest| rest.strip_suffix(')'))
        .unwrap_or(text);
    let (numerator, denominator) = match text.split_once('/') {
        | Some((num, den)) => (num, den),
        | None => (text, "1"),
    };
    let numerator = numerator.parse::<BigInt>().map_err(|_| {
        DimError::Parse {
            text: whole.to_string(),
            reason: format!("`{text}` is not a rational exponent"),
        }
    })?;
    let denominator = denominator.parse::<BigInt>().map_err(|_| {
        DimError::Parse {
            text: whole.to_string(),
            reason: format!("`{text}` is not a rational exponent"),
        }
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
#[allow(clippy::match_same_arms)] // grouped by domain; identical bodies document families
pub fn default_rule(name: &str) -> Option<DimRule> {
    let rule = match name {
        // `piecewise` is deliberately absent: a piecewise term interleaves
        // conditions with values, so neither `Sum` nor any other uniform rule
        // is correct. A rule set that registers it must supply its own rule
        // (the coverage scanner enforces that).
        | "add" | "sub" | "min" | "max" => DimRule::Sum,
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
        | "acsch" | "erf" | "erfc" | "erfi" | "fresnels" | "fresnelc" | "erfcinv" | "gamma"
        | "lgamma" | "digamma" | "polygamma" | "zeta" | "hurwitz_zeta" | "beta" | "lambertw"
        | "sinc" | "besselj" | "bessely" | "besseli" | "besselk" | "airy_ai" | "airy_bi"
        | "airy_aip" | "airy_bip" | "pi" | "e" | "true_pi" | "euler_gamma" | "i"
        | "golden_ratio" => DimRule::Dimensionless,
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
        // verdicts. A contour lives in a dimensionless complex plane: **every**
        // operator that consumes a complex variable checks that it is
        // dimensionless (`dimensionless_complex_variable` via the `Custom`
        // rules below), while `continue_along`/`mobius_apply`/`residue` keep
        // the function's dimension and `cauchy_derivative` divides by `z`.
        | "I" => DimRule::Dimensionless,
        // Branch wrappers are complex-plane functions (`log_branch(z, k) =
        // ln|z| + I(arg z + 2 pi k)`, ...): every argument, including the
        // branch index, is a pure number and the value is one. `Same(0)` used
        // to return `dim(z)` here, which contradicted the defining rewrites
        // (`power_branch` vs `exp(w log_branch(z, k))`) and made the
        // per-rewrite dimension scan fail.
        | "log_branch" | "sqrt_branch" | "root_branch" | "power_branch" | "asin_branch"
        | "acos_branch" | "atan_branch" | "asinh_branch" | "acosh_branch" | "atanh_branch" => {
            DimRule::Dimensionless
        },
        | "distance" => DimRule::Same(0),
        | "residue" => DimRule::Custom(residue_dimension),
        | "continue_along" => DimRule::Custom(continue_along_dimension),
        | "mobius_apply" => DimRule::Custom(mobius_apply_dimension),
        | "radius_of_convergence" => DimRule::Custom(radius_of_convergence_dimension),
        | "contour_integral" | "cauchy_integral" | "contour_integral_residue_theorem" => {
            DimRule::Custom(contour_dimension)
        },
        | "cauchy_derivative" => DimRule::Custom(cauchy_derivative_dimension),
        | "count_zeros_poles" => DimRule::Custom(count_zeros_poles_dimension),
        | "winding_number" => DimRule::Custom(winding_number_dimension),
        | "poles" => DimRule::Custom(poles_dimension),
        | "zeros_of" => DimRule::Custom(zeros_of_dimension),
        | "singularity" => DimRule::Custom(singularity_dimension),
        | "branch_cut" => DimRule::Custom(branch_cut_dimension),
        | "check_analytic" => DimRule::Custom(check_analytic_dimension),
        | "is_inside_contour" => DimRule::Custom(is_inside_contour_dimension),
        | "circle"
        | "path"
        | "polygon"
        | "ray"
        | "mobius_compose"
        | "mobius_inverse"
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

/// The dimension rule of `matexp(A, t)` (`exp(t A)`), installed as a
/// [`DimRule::Custom`].
///
/// Only the *product* `dim(t) * dim(A)` must be dimensionless, so `A` may be
/// a rate and `t` a time: `exp(A t)` is a pure number. A product with a
/// dimension is a [`DimError::Mismatch`]; an unknown product (an undeclared
/// symbol) is accepted, because "unknown" is not "wrong". This is the
/// relaxation of the former `Dimensionless` rule, which required `A` and `t`
/// to be dimensionless separately.
///
/// # Errors
/// A [`DimError::Mismatch`] when `dim(t) * dim(A)` has a concrete dimension.
pub fn matexp_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let matrix = args.first().cloned().unwrap_or_else(DimTerm::any);
    let scale = args.get(1).cloned().unwrap_or_else(DimTerm::any);
    let product = matrix.mul(&scale);
    // An unknown product is not a conflict: leave it unknown (the result is
    // still the pure number `1`).
    if product.is_any() || product.has_unknown() || product.is_dimensionless() {
        return Ok(ctx.dimensionless());
    }
    Err(DimError::Mismatch {
        left: product.describe(ctx.arena),
        right: String::from("1"),
        context: ctx.graph.display(ctx.node),
    })
}

/// The complex variable of a contour/Cauchy operator, as
/// `(dim(function), dim(z))`.
///
/// A contour lives in a dimensionless complex plane, so `z` must be
/// dimensionless; an undeclared `z` is an unknown and is accepted rather than
/// guessed. A dimensionful `z` is a [`DimError::NonDimensionless`] conflict.
///
/// # Errors
/// A [`DimError::NonDimensionless`] when `dim(z)` is a concrete non-`1`
/// dimension.
fn dimensionless_complex_variable(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
    variable: usize,
) -> Result<(DimTerm, DimTerm), DimError> {
    let function = args.first().cloned().unwrap_or_else(DimTerm::any);
    let variable_dimension = args.get(variable).cloned().unwrap_or_else(DimTerm::any);
    if variable_dimension.is_any()
        || variable_dimension.has_unknown()
        || variable_dimension.is_dimensionless()
    {
        return Ok((function, variable_dimension));
    }
    Err(DimError::NonDimensionless {
        op: ctx.graph.ops().get(ctx.graph.op(ctx.node)).name.to_string(),
        argument: ctx
            .graph
            .children(ctx.node)
            .get(variable)
            .map_or_else(String::new, |&child| ctx.graph.display(child)),
        found: variable_dimension.describe(ctx.arena),
    })
}

/// `contour_integral`, `cauchy_integral` and
/// `contour_integral_residue_theorem`: the complex variable must be
/// dimensionless, and the value carries the dimension of the integrand.
///
/// The differential `dz` contributes `dim(z) = 1`, so the result is
/// `dim(f)`.
///
/// # Errors
/// A [`DimError::NonDimensionless`] when the complex variable has a
/// dimension.
pub fn contour_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (function, _variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(function)
}

/// `cauchy_derivative`: the complex variable must be dimensionless, and the
/// value divides by `(z - a)^(n + 1)` and multiplies by `dz`.
///
/// With `dim(z) = 1` the result is `dim(f)`; an undeclared `z` keeps the
/// previous `Diff` result `dim(f) / dim(z)` because its dimension is unknown.
///
/// # Errors
/// A [`DimError::NonDimensionless`] when the complex variable has a
/// dimension.
pub fn cauchy_derivative_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (function, variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(function.mul(&variable.inverse()))
}

/// A [`DimRule::Custom`] rule for a complex-analysis operator: check the
/// complex variable at `$variable`, then report the operator's fixed
/// dimensional kind.
macro_rules! complex_variable_rule {
    ($name:ident, $variable:expr, $result:expr) => {
        fn $name(
            ctx: &DimCtx<'_>,
            args: &[DimTerm],
        ) -> Result<DimTerm, DimError> {
            dimensionless_complex_variable(ctx, args, $variable)?;
            Ok($result)
        }
    };
}

// Structured results (point and singularity sets, branch cuts). The complex
// variable is checked, and the value is typed as a *container*: it is
// non-scalar, but it stays compatible with the dimension of the reduced list
// literal, so a request and its answer unify instead of reporting a spurious
// `non-metric vs 1` conflict.
complex_variable_rule!(poles_dimension, 1, DimTerm::container());
complex_variable_rule!(zeros_of_dimension, 1, DimTerm::container());
complex_variable_rule!(singularity_dimension, 1, DimTerm::container());
complex_variable_rule!(branch_cut_dimension, 1, DimTerm::container());
// Boolean verdicts stay non-metric, exactly as their literal answers are.
complex_variable_rule!(check_analytic_dimension, 1, DimTerm::non_metric());
complex_variable_rule!(is_inside_contour_dimension, 0, DimTerm::non_metric());

// Counts are pure numbers, and the contour argument is a non-scalar: unifying
// *every* argument with `Dimensionless` (the former rule) reported a spurious
// conflict on the contour.
complex_variable_rule!(count_zeros_poles_dimension, 1, DimTerm::dimensionless());
complex_variable_rule!(winding_number_dimension, 0, DimTerm::dimensionless());

/// `residue(f, z, a)`: the residue of `f` at `a` has the dimension of `f`;
/// the complex variable must be dimensionless.
fn residue_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (function, _variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(function)
}

/// `continue_along(f, z, list(...), order)`: the continued expansion is a
/// value of `f` on the dimensionless z-plane.
fn continue_along_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (function, _variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(function)
}

/// `mobius_apply(M, z)`: the image is a point of the dimensionless z-plane.
fn mobius_apply_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (function, _variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(function)
}

/// `radius_of_convergence(f, z, a)`: a radius in the dimensionless z-plane,
/// so it keeps `dim(z)` after the complex variable has been checked.
fn radius_of_convergence_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let (_function, variable) = dimensionless_complex_variable(ctx, args, 1)?;
    Ok(variable)
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
        let mut unifier = crate::graph::dim_infer::Unifier::default();
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
        let mut unifier = crate::graph::dim_infer::Unifier::default();
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
    fn a_base_defined_on_a_clone_does_not_reach_its_source() {
        // The arena is copy-on-write: a request may define a base dimension of
        // its own, and the template it was cloned from must not see it.
        let template = Graph::new();
        let mut request = template.clone();
        request.define_base_dimension("Currency");
        assert!(request.dim_arena().base_id("Currency").is_some());
        assert!(
            template.dim_arena().base_id("Currency").is_none(),
            "the template picked up a base from a request"
        );
        // A request that defines nothing shares the arena rather than copying
        // it, so the common path pays nothing.
        let untouched = template.clone();
        assert!(untouched.dim_arena().base_id("Currency").is_none());
        // Writing through the mutable accessor detaches the same way.
        let mut second = template.clone();
        let term = second.with_dim_arena_mut(|arena| {
            let base = arena.define_base("Currency");
            DimTerm::base(base, BigRational::one())
        });
        assert!(!term.is_dimensionless());
        assert!(template.dim_arena().base_id("Currency").is_none());
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

    /// Every `DimError` variant has a stable `snake_case` code, the codes are
    /// distinct, and the method is total (a new variant fails to compile
    /// until it is given one). `Display` is unchanged: it still prints the
    /// human-readable message.
    #[test]
    fn dim_error_codes_pin_every_variant() {
        let variants = [
            (
                DimError::Mismatch {
                    left: String::from("Length"),
                    right: String::from("Time"),
                    context: String::from("x + y"),
                },
                "mismatch",
            ),
            (
                DimError::NonDimensionless {
                    op: String::from("exp"),
                    argument: String::from("x"),
                    found: String::from("Length"),
                },
                "non_dimensionless",
            ),
            (
                DimError::SymbolicConstraint {
                    constraint: String::from("n = 1"),
                    context: String::from("x^n"),
                },
                "symbolic_constraint",
            ),
            (
                DimError::UnknownUnit {
                    unit: String::from("parsec"),
                    context: None,
                },
                "unknown_unit",
            ),
            (
                DimError::UnknownUnit {
                    unit: String::from("parsec"),
                    context: Some(String::from("quantity(1, parsec)")),
                },
                "unknown_unit",
            ),
            (
                DimError::Affine {
                    operation: String::from("quantity(20, degC) * 2"),
                },
                "affine",
            ),
            (
                DimError::MixedList {
                    context: String::from("list(x, y)"),
                },
                "mixed_list",
            ),
            (
                DimError::Unannotated {
                    op: String::from("mystery"),
                },
                "unannotated",
            ),
            (
                DimError::Parse {
                    text: String::from("m^"),
                    reason: String::from("missing exponent"),
                },
                "parse",
            ),
        ];
        let mut seen: Vec<&str> = Vec::new();
        for (error, code) in &variants {
            assert_eq!(error.code(), *code, "wrong code for {error:?}");
            assert!(
                code.chars().all(|c| c.is_ascii_lowercase() || c == '_'),
                "`{code}` is not snake_case"
            );
            if !seen.contains(code) {
                seen.push(code);
            }
            // The machine-readable code is additive: `Display` is unchanged.
            assert!(!error.to_string().is_empty(), "empty Display for {error:?}");
        }
        // Eight variants, one code each: the two `UnknownUnit` contexts share
        // `unknown_unit` (the blocking difference is `is_blocking`).
        assert_eq!(seen.len(), 8, "{seen:?}");
        assert!(!seen.contains(&"unresolved_unit"));
        assert_eq!(
            DimError::UnknownUnit {
                unit: String::from("parsec"),
                context: Some(String::from("quantity(1, parsec)")),
            }
            .code(),
            "unknown_unit"
        );
    }

    #[test]
    fn all_dim_error_codes_are_non_empty_snake_case() {
        let all_variants = [
            DimError::Mismatch {
                left: "m".into(),
                right: "s".into(),
                context: "x + y".into(),
            },
            DimError::NonDimensionless {
                op: "exp".into(),
                argument: "x".into(),
                found: "m".into(),
            },
            DimError::SymbolicConstraint {
                constraint: "n = 1".into(),
                context: "x^n".into(),
            },
            DimError::UnknownUnit {
                unit: "parsec".into(),
                context: None,
            },
            DimError::UnknownUnit {
                unit: "parsec".into(),
                context: Some("quantity(1, parsec)".into()),
            },
            DimError::Affine {
                operation: "degC + degC".into(),
            },
            DimError::MixedList {
                context: "list(1 m, 2 s)".into(),
            },
            DimError::Unannotated { op: "foo".into() },
            DimError::Parse {
                text: "invalid".into(),
                reason: "syntax".into(),
            },
        ];

        for err in &all_variants {
            let code = err.code();
            assert!(!code.is_empty(), "code should not be empty for {err:?}");
            assert!(
                code.chars().all(|c| c.is_ascii_lowercase() || c == '_'),
                "code '{code}' should be snake_case"
            );
            assert!(
                !code.starts_with('_') && !code.ends_with('_'),
                "code '{code}' should not start or end with underscore"
            );
        }
    }

    /// `matfun(f, x, A)` is only homogeneous for a single total degree; a
    /// non-homogeneous or non-polynomial `f` still requires a dimensionless
    /// matrix. `matexp(A, t)` is the weaker case: only the *product* `t A`
    /// has to be dimensionless, so a rate times a time is accepted.
    #[test]
    fn matexp_requires_only_a_dimensionless_product() {
        // An evolution operator: `A` in 1/s, `t` in s.
        let rate = dimension_of("matexp(list(list(quantity(1, 1/s))), quantity(2, s))", &[]);
        assert_eq!(rate.dimension.to_string(), "1");
        assert!(rate.errors.is_empty(), "{:?}", rate.errors);
        // Dimensionful `A` and `t` whose product cancels.
        let cancelling = dimension_of(
            "matexp(list(list(quantity(1, m/s))), quantity(2, s/m))",
            &[],
        );
        assert_eq!(cancelling.dimension.to_string(), "1");
        assert!(cancelling.errors.is_empty(), "{:?}", cancelling.errors);
        // Two unknown symbols give an unknown product: accepted as unknown,
        // not guessed to be dimensionless.
        let symbolic = dimension_of("matexp(list(list(x)), t)", &[]);
        assert_eq!(symbolic.dimension.to_string(), "1");
        assert!(symbolic.errors.is_empty(), "{:?}", symbolic.errors);
        // A product with a dimension stays an error, for either factor.
        for source in [
            "matexp(list(list(quantity(1, m))), 4)",
            "matexp(list(list(quantity(1, m))), quantity(2, s))",
            "matexp(list(list(1)), quantity(2, m))",
        ] {
            let outcome = dimension_of(source, &[]);
            assert!(
                outcome
                    .errors
                    .iter()
                    .any(|error| matches!(error, DimError::Mismatch { .. })),
                "`{source}`: {:?}",
                outcome.errors
            );
        }
    }

    /// A contour/Cauchy operator is defined on a dimensionless complex plane:
    /// a dimensionful complex variable is a conflict, while the dimensionless
    /// and the undeclared cases are accepted.
    #[test]
    fn contour_operators_require_a_dimensionless_complex_variable() {
        let accepted = dimension_of("contour_integral(z^2, z, c)", &[("z", "1")]);
        assert_eq!(accepted.dimension.to_string(), "1");
        assert!(accepted.errors.is_empty(), "{:?}", accepted.errors);

        // An undeclared complex variable is an unknown, so the contour rule
        // accepts it (the enclosing `exp` binds it to dimensionless here).
        let undeclared = dimension_of("cauchy_integral(exp(z), z, 0)", &[]);
        assert_eq!(undeclared.dimension.to_string(), "1");
        assert!(undeclared.errors.is_empty(), "{:?}", undeclared.errors);

        // A dimensionful `z` is refused by every operator of the family.
        for source in [
            "contour_integral(z^2, z, c)",
            "cauchy_integral(exp(z), z, 0)",
            "cauchy_derivative(exp(z), z, 0, 1)",
            "contour_integral_residue_theorem(exp(z), z, list(0))",
        ] {
            let outcome = dimension_of(source, &[("z", "m")]);
            assert!(
                outcome.errors.iter().any(|error| {
                    matches!(
                        error,
                        DimError::NonDimensionless { op, .. }
                            if op == "contour_integral"
                                || op == "cauchy_integral"
                                || op == "cauchy_derivative"
                                || op == "contour_integral_residue_theorem"
                    )
                }),
                "`{source}`: {:?}",
                outcome.errors
            );
        }
    }

    /// Branch wrappers live on the dimensionless complex plane: every
    /// argument (including the branch index) is a pure number and the value is
    /// dimensionless. This keeps `power_branch` consistent with its defining
    /// rewrite `power_branch(z, w, k) => exp(w log_branch(z, k))`: the old
    /// `Same(0)` rule reported `dim(z)` on the left while the right side was
    /// dimensionless, which the per-rewrite scan caught.
    #[test]
    fn branch_wrappers_are_dimensionless_plane_functions() {
        let accepted = dimension_of("power_branch(z, w, 2)", &[("z", "1"), ("w", "1")]);
        assert_eq!(accepted.dimension.to_string(), "1");
        assert!(accepted.errors.is_empty(), "{:?}", accepted.errors);

        for source in [
            "log_branch(z, 0)",
            "sqrt_branch(z, 2)",
            "root_branch(z, 3, 1)",
            "power_branch(z, w, 2)",
        ] {
            let outcome = dimension_of(source, &[("z", "m"), ("w", "1/m")]);
            assert!(
                outcome
                    .errors
                    .iter()
                    .any(|error| matches!(error, DimError::NonDimensionless { .. })),
                "`{source}`: {:?}",
                outcome.errors
            );
        }
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
    fn a_declared_function_value_has_its_dimension() {
        let declared = dimension_of("f(x)", &[("f", "m/s^2")]);
        assert_eq!(declared.dimension.to_string(), "Length/Time^2");
        assert!(declared.errors.is_empty(), "{:?}", declared.errors);
        // Without a declaration the application is an unknown value.
        let unknown = dimension_of("g(x)", &[]);
        assert!(unknown.dimension.is_unknown(), "{}", unknown.dimension);
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

    /// Operators whose dimension rule is missing or explicitly `Unspecified`.
    ///
    /// `lit` and `sym` are leaves with no operator rule: the inference engine
    /// reads their payloads directly (a number, a dimension variable, or a
    /// non-metric truth value/name).
    fn unannotated_operators(graph: &Graph) -> Vec<String> {
        graph
            .ops()
            .iter()
            .filter(|(_, descriptor)| !matches!(&*descriptor.name, "lit" | "sym"))
            .filter(|(op, _)| {
                graph
                    .ops()
                    .attr::<DimRule>(*op)
                    .is_none_or(|rule| matches!(rule, DimRule::Unspecified))
            })
            .map(|(_, descriptor)| descriptor.name.to_string())
            .collect()
    }

    /// The only operators allowed to be dimensionally polymorphic (`Any`).
    /// Every entry is a documented bypass: `oo` and order terms are
    /// dimensionless-or-polymorphic by construction, and a definiteness label
    /// is a symbol whose dimension is unknown.
    const DOCUMENTED_BYPASS: &[&str] = &["oo", "infinity", "big_o", "hessian_definiteness"];

    #[test]
    fn every_standard_operator_has_a_dimension_rule() {
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        let missing = unannotated_operators(&graph);
        assert!(
            missing.is_empty(),
            "{} operators have no dimension rule: {missing:?}",
            missing.len()
        );
    }

    /// The scanner is not vacuous: an operator registered without a rule is
    /// reported. Removing an `i.dim(...)` line therefore fails the suite.
    #[test]
    fn the_coverage_scanner_detects_a_missing_rule() {
        let gap = crate::graph::RuleSet::new("dimension-audit-gap", |i| {
            i.op(crate::graph::OpDescriptor::new(
                "audit_unannotated_operator",
                crate::graph::Arity::Fixed(0),
            ))?;
            Ok(())
        });
        let mut graph = Graph::new();
        let mut sets = crate::rules::standard();
        sets.push(gap);
        crate::graph::Engine::install(&mut graph, &sets).unwrap_or_else(|error| panic!("{error}"));
        let missing = unannotated_operators(&graph);
        assert!(
            missing
                .iter()
                .any(|name| name == "audit_unannotated_operator"),
            "the scanner did not report the unannotated operator: {missing:?}"
        );
    }

    /// `sde_wiener(seed, m, steps, h)` has the dimension `sqrt(dim(h))`: the
    /// rooted operand is argument 3. This is the case `DimRule::RootOf` exists
    /// for — `Root(n)` would take the root of the seed.
    #[test]
    fn noise_operators_state_their_root_dimension() {
        let dimensionless = dimension_of("sde_wiener(11, 2, 4, 1/100)", &[]);
        assert_eq!(dimensionless.dimension.to_string(), "1");
        assert!(
            dimensionless.errors.is_empty(),
            "{:?}",
            dimensionless.errors
        );
        let time = dimension_of("sde_wiener(11, 2, 4, quantity(1/100, s))", &[]);
        assert_eq!(time.dimension.to_string(), "sqrt(Time)");
        // Levy areas are iterated integrals of two increments: dimension T.
        let areas = dimension_of(
            "sde_levy_areas(11, quantity(1/100, s), list(1/10, -1/5), 8)",
            &[],
        );
        assert_eq!(areas.dimension.to_string(), "Time");
        assert!(areas.errors.is_empty(), "{:?}", areas.errors);
        // The scheme metadata is a label vector, not a physical quantity. The
        // name is a string *value*, and the surface syntax has no string
        // literal, so the term is built through the node API.
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        let op = graph
            .ops()
            .lookup("sde_scheme_info")
            .unwrap_or_else(|| panic!("sde_scheme_info is registered"));
        let text = graph.lit(crate::graph::Payload::Str(std::sync::Arc::from("milstein")));
        let call = graph.node(op, &[text]);
        let metadata = graph.dimension(call);
        assert_eq!(metadata.dimension.to_string(), "non-metric");
        assert!(metadata.errors.is_empty(), "{:?}", metadata.errors);
    }

    /// No operator bypasses dimensional checking by accident.
    #[test]
    fn only_documented_operators_bypass_dimensions() {
        let mut graph = Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|error| panic!("{error}"));
        let bypass: Vec<String> = graph
            .ops()
            .iter()
            .filter(|(op, _)| {
                graph
                    .ops()
                    .attr::<DimRule>(*op)
                    .is_some_and(|rule| matches!(rule, DimRule::Any))
            })
            .map(|(_, descriptor)| descriptor.name.to_string())
            .filter(|name| !DOCUMENTED_BYPASS.contains(&name.as_str()))
            .collect();
        assert!(
            bypass.is_empty(),
            "operators bypass dimensions without a documented reason: {bypass:?}"
        );
    }
}
