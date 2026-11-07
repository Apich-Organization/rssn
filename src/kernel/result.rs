//! Structured results: answers that are more than one number.
//!
//! # The problem
//!
//! A solver that returns `list(1, -1)` has thrown away everything a caller
//! needs: that `1` is a double root, that the answer is only valid when `a`
//! is non-zero, that the second element came from a branch of `sqrt`. A
//! simulation that returns a matrix of doubles has thrown away which column
//! is which state variable, and in what units.
//!
//! # The answer
//!
//! Every structured answer is a [`Value::Record`] whose field names, order
//! and types are declared by a [`FlattenSchema`]. That is a deliberate
//! choice:
//!
//! * it crosses the dynamic ABI unchanged — a host reads a record the same
//!   way it reads a number;
//! * it is serialisable with `bincode-next` for free, because it is just a
//!   `Value`;
//! * the schema travels with it, so a host can generate a typed view without
//!   a second source of truth.
//!
//! Four shapes cover the results the engine produces today:
//!
//! | Type | Record fields |
//! |---|---|
//! | [`Solutions`] | `variable`, `roots`, `multiplicities`, `conditions`, `complete` |
//! | [`OdeSolutionView`] | `t`, `y`, `flatten_order`, `events`, `nfev`, `rejected` |
//! | [`SimResult`] | `scenario`, `t`, `fields`, `schema` |
//!
//! # Flattening
//!
//! A state vector `y` is a list of rows. [`OdeSolutionView::flatten`] turns
//! it into one tensor of shape `[steps, components]` in
//! [`OdeSolutionView::flatten_order`], which is the layout a host's array
//! library wants and the layout a CSV writer needs. The order is recorded, not
//! assumed, because a model may declare its fields in one order and integrate
//! them in another.

use serde::Deserialize;
use serde::Serialize;

use crate::kernel::signature::DimExpr;
use crate::kernel::value::QuantityValue;
use crate::kernel::value::RecordValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;

/// Why a value did not match the schema it claims.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SchemaError {
    /// A required field is missing.
    Missing {
        /// The record's schema.
        schema: String,
        /// The absent field.
        field: String,
    },
    /// A field has the wrong shape of value.
    Kind {
        /// The record's schema.
        schema: String,
        /// The offending field.
        field: String,
        /// What the schema asked for.
        expected: String,
        /// What the value was.
        got: String,
    },
    /// A field appeared that the schema does not declare.
    Unexpected {
        /// The record's schema.
        schema: String,
        /// The undeclared field.
        field: String,
    },
    /// The record is a `Value` of some other kind.
    NotARecord {
        /// The record's schema.
        schema: String,
    },
    /// A list whose length must match another field's does not.
    Length {
        /// The record's schema.
        schema: String,
        /// The field that is too short or too long.
        field: String,
        /// What the schema requires.
        expected: usize,
        /// What the value has.
        got: usize,
    },
}

impl std::fmt::Display for SchemaError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Missing { schema, field } => {
                write!(f, "`{schema}` is missing the field `{field}`")
            },
            | Self::Kind {
                schema,
                field,
                expected,
                got,
            } => write!(f, "`{schema}.{field}` should be {expected}, got {got}"),
            | Self::Unexpected { schema, field } => {
                write!(f, "`{schema}` has an undeclared field `{field}`")
            },
            | Self::NotARecord { schema } => write!(f, "`{schema}` should be a record"),
            | Self::Length {
                schema,
                field,
                expected,
                got,
            } => {
                write!(
                    f,
                    "`{schema}.{field}` should have {expected} entries, got {got}"
                )
            },
        }
    }
}

impl std::error::Error for SchemaError {}

/// What a field of a structured result holds.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldKind {
    /// One number, possibly dimensional.
    Scalar,
    /// A tensor of the given rank; `None` when it varies.
    Tensor(Option<usize>),
    /// Text.
    Text,
    /// A truth value.
    Bool,
    /// A heterogeneous list.
    List,
    /// A nested record with its own schema.
    Record(String),
    /// Anything; the schema declines to say.
    Any,
}

impl std::fmt::Display for FieldKind {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Scalar => f.write_str("a scalar"),
            | Self::Tensor(Some(rank)) => write!(f, "a rank-{rank} tensor"),
            | Self::Tensor(None) => f.write_str("a tensor"),
            | Self::Text => f.write_str("text"),
            | Self::Bool => f.write_str("a boolean"),
            | Self::List => f.write_str("a list"),
            | Self::Record(name) => write!(f, "a `{name}` record"),
            | Self::Any => f.write_str("anything"),
        }
    }
}

impl FieldKind {
    /// Whether `value` is of this kind.
    #[must_use]
    pub fn accepts(
        &self,
        value: &Value,
    ) -> bool {
        match self {
            | Self::Scalar => {
                value.is_number() || matches!(value, Value::Quantity(_) | Value::Unit(_))
            },
            | Self::Tensor(rank) => {
                match value {
                    | Value::Tensor(tensor) => rank.is_none_or(|rank| tensor.rank() == rank),
                    | Value::Sparse(sparse) => rank.is_none_or(|rank| sparse.shape.len() == rank),
                    | _ => false,
                }
            },
            | Self::Text => matches!(value, Value::String(_) | Value::Symbol(_)),
            | Self::Bool => matches!(value, Value::Bool(_)),
            | Self::List => matches!(value, Value::Tensor(_) | Value::Sparse(_)),
            | Self::Record(_) => matches!(value, Value::Record(_)),
            | Self::Any => true,
        }
    }
}

/// One field of a structured result.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Field {
    /// The name a host looks the value up by.
    pub name: String,
    /// What it holds.
    pub kind: FieldKind,
    /// Its dimension, when the schema knows it.
    pub dim: Option<DimExpr>,
    /// What it means.
    pub doc: String,
}

impl Field {
    /// A field of `kind` named `name`.
    #[must_use]
    pub fn new(
        name: &str,
        kind: FieldKind,
    ) -> Self {
        Self {
            name: name.to_owned(),
            kind,
            dim: None,
            doc: String::new(),
        }
    }

    /// Declares the field dimensional.
    #[must_use]
    pub fn with_dim(
        mut self,
        dim: DimExpr,
    ) -> Self {
        self.dim = Some(dim);
        self
    }

    /// Documents the field.
    #[must_use]
    pub fn doc(
        mut self,
        doc: &str,
    ) -> Self {
        doc.clone_into(&mut self.doc);
        self
    }
}

/// The layout of a structured result.
///
/// Field order is part of the contract: it is what
/// [`OdeSolutionView::flatten_order`] and a generated host struct agree on.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FlattenSchema {
    /// What the record is called (`solutions`, `ode_solution`, ...).
    pub name: String,
    /// Schema version, bumped when a field changes meaning.
    pub version: u32,
    /// The fields, in declaration order.
    pub fields: Vec<Field>,
}

impl FlattenSchema {
    /// A schema with `fields`, in order.
    #[must_use]
    pub fn new(
        name: &str,
        version: u32,
        fields: Vec<Field>,
    ) -> Self {
        Self {
            name: name.to_owned(),
            version,
            fields,
        }
    }

    /// The field called `name`.
    #[must_use]
    pub fn field(
        &self,
        name: &str,
    ) -> Option<&Field> {
        self.fields.iter().find(|field| field.name == name)
    }

    /// How many fields the record has.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.fields.len()
    }

    /// Whether the schema has no fields.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// The field names, in order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.fields
            .iter()
            .map(|field| field.name.as_str())
            .collect()
    }

    /// Checks `record` against the schema.
    ///
    /// # Errors
    /// Reports the first missing, mistyped or undeclared field.
    pub fn check(
        &self,
        record: &RecordValue,
    ) -> Result<(), SchemaError> {
        for field in &self.fields {
            let Some(value) = record.get(&field.name) else {
                return Err(SchemaError::Missing {
                    schema: self.name.clone(),
                    field: field.name.clone(),
                });
            };
            if !field.kind.accepts(value) {
                return Err(SchemaError::Kind {
                    schema: self.name.clone(),
                    field: field.name.clone(),
                    expected: field.kind.to_string(),
                    got: value.type_expr().to_string(),
                });
            }
        }
        for (name, _) in &record.fields {
            if self.field(name).is_none() {
                return Err(SchemaError::Unexpected {
                    schema: self.name.clone(),
                    field: name.clone(),
                });
            }
        }
        Ok(())
    }

    /// The schema itself as a value, so that a result can carry its own
    /// layout.
    #[must_use]
    pub fn to_value(&self) -> Value {
        let fields = self
            .fields
            .iter()
            .map(|field| {
                let mut entry = vec![
                    (String::from("name"), Value::String(field.name.clone())),
                    (String::from("kind"), Value::String(field.kind.to_string())),
                ];
                if let Some(dim) = &field.dim {
                    entry.push((
                        String::from("dimension"),
                        Value::String(crate::kernel::signature::dim_name(dim)),
                    ));
                }
                if !field.doc.is_empty() {
                    entry.push((String::from("doc"), Value::String(field.doc.clone())));
                }
                Value::Record(RecordValue { fields: entry })
            })
            .collect();
        Value::Record(RecordValue {
            fields: vec![
                (String::from("name"), Value::String(self.name.clone())),
                (String::from("version"), Value::int(i64::from(self.version))),
                (
                    String::from("fields"),
                    Value::Tensor(TensorValue::vector(fields)),
                ),
            ],
        })
    }
}

// ---------------------------------------------------------------------------
// Solutions
// ---------------------------------------------------------------------------

/// The kind of a side condition.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConditionKind {
    /// A parameter must not be zero.
    NonZero,
    /// A parameter must be positive.
    Positive,
    /// A parameter must be non-negative (`a >= 0`, for `sqrt(a)`).
    Nonnegative,
    /// A parameter must be negative.
    Negative,
    /// A parameter must be real.
    Real,
    /// A branch was chosen; the other branch is the negation.
    Branch,
    /// The roots must be distinct.
    Distinct,
    /// The result holds only on an interval.
    Interval,
    /// Anything else, described in text.
    Other,
}

/// A condition under which a solution set is valid.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Condition {
    /// What kind of condition it is.
    pub kind: ConditionKind,
    /// The condition as an expression, when the solver kept one.
    pub expression: Option<Value>,
    /// The condition in words, for a host that cannot evaluate expressions.
    pub text: String,
}

impl Condition {
    /// A condition of `kind` described by `text`.
    #[must_use]
    pub fn new(
        kind: ConditionKind,
        text: &str,
    ) -> Self {
        Self {
            kind,
            expression: None,
            text: text.to_owned(),
        }
    }

    /// Attaches the condition's expression.
    #[must_use]
    pub fn with_expression(
        mut self,
        expression: Value,
    ) -> Self {
        self.expression = Some(expression);
        self
    }

    /// The condition as a record field.
    #[must_use]
    pub fn to_value(&self) -> Value {
        let kind = format!("{:?}", self.kind).to_lowercase();
        let mut fields = vec![
            (String::from("kind"), Value::String(kind)),
            (String::from("text"), Value::String(self.text.clone())),
        ];
        if let Some(expression) = &self.expression {
            fields.push((String::from("expression"), expression.clone()));
        }
        Value::Record(RecordValue { fields })
    }
}

/// A solution set: the roots of an equation, with what qualifies them.
///
/// The point of the type is that `roots`, `multiplicities` and `conditions`
/// are three parallel facts about one answer, and a bare list cannot carry
/// two of them.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Solutions {
    /// The unknown that was solved for, when it is a symbol.
    pub variable: Option<String>,
    /// The roots, in the order the solver produced them.
    pub roots: Vec<Value>,
    /// How many times each root occurs.
    ///
    /// Empty means the solver did not distinguish them: every root is
    /// listed once, and [`Solutions::multiplicity`] then answers one. That
    /// is different from a reported `1`, which is a claim.
    ///
    /// `solve_structured` fills it from the exact square-free decomposition of
    /// a univariate polynomial: the multiplicities of exact rational roots of
    /// a rational polynomial, and all ones when the polynomial is square-free
    /// (irrational or symbolic roots included; for parameters the
    /// [`ConditionKind::Distinct`] conditions say where roots collide). It
    /// stays empty when the association between a root and a repeated factor
    /// would have to be numeric — such a root is not guessed at.
    pub multiplicities: Vec<u32>,
    /// What must hold for the answer to be valid.
    ///
    /// Empty when the solver stated no condition. `solve_structured` reports a
    /// condition it can derive from the answer as written: a parameter in a
    /// denominator of a root (`-b/a` from `a*x + b = 0` needs `a != 0`), an
    /// even root (`a >= 0` for `sqrt(a)`, `a > 0` for `a^(-1/2)`), a
    /// logarithm (`ln(a)` needs `a > 0`), an inverse sine or cosine
    /// (`-1 <= c <= 1`), and, for a polynomial with parameters, where its
    /// roots coincide. Conditions already implied by the parameters'
    /// assumptions are omitted; branch choices the solver makes beyond the
    /// principal branch (see `complete`) are not decided here.
    pub conditions: Vec<Condition>,
    /// Whether the solver believes it found every root.
    ///
    /// False for a solution of a periodic equation (`sin(x) = 0`) found on the
    /// principal branch: the listed roots are representatives and the full
    /// set is `solve_general`'s families; a [`ConditionKind::Branch`]
    /// condition says so.
    pub complete: bool,
}

impl Solutions {
    /// A solution set of `roots`, each a simple root.
    ///
    /// `complete` starts **false**: a caller that builds a `Solutions` by hand
    /// has not run an exhaustive method, so claiming completeness would be a
    /// statement it cannot support. Call [`Solutions::complete`] when the
    /// method you used has earned it — `solve_structured` sets it from the
    /// shape of the solver's answer, not by assumption.
    #[must_use]
    pub const fn new(roots: Vec<Value>) -> Self {
        Self {
            variable: None,
            roots,
            multiplicities: Vec::new(),
            conditions: Vec::new(),
            complete: false,
        }
    }

    /// Declares the set complete.
    #[must_use]
    pub const fn complete(mut self) -> Self {
        self.complete = true;
        self
    }

    /// The layout every `Solutions` record has.
    #[must_use]
    pub fn schema() -> FlattenSchema {
        FlattenSchema::new(
            "solutions",
            1,
            vec![
                Field::new("variable", FieldKind::Text).doc("the unknown that was solved for"),
                Field::new("roots", FieldKind::List).doc("the roots, in solver order"),
                Field::new("multiplicities", FieldKind::List)
                    .doc("how many times each root occurs; empty means all one"),
                Field::new("conditions", FieldKind::List)
                    .doc("what must hold for the answer to be valid"),
                Field::new("complete", FieldKind::Bool)
                    .doc("whether every root is believed to be present"),
            ],
        )
    }

    /// Whether the solver reported multiplicities at all.
    #[must_use]
    pub const fn has_multiplicities(&self) -> bool {
        !self.multiplicities.is_empty()
    }

    /// The multiplicity of root `index`, defaulting to one.
    #[must_use]
    pub fn multiplicity(
        &self,
        index: usize,
    ) -> u32 {
        self.multiplicities.get(index).copied().unwrap_or(1)
    }

    /// The number of roots.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.roots.len()
    }

    /// Whether there are no roots.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.roots.is_empty()
    }

    /// The answer as a record.
    #[must_use]
    pub fn to_value(&self) -> Value {
        let multiplicities = if self.multiplicities.is_empty() {
            Vec::new()
        } else {
            self.roots
                .iter()
                .enumerate()
                .map(|(index, _)| Value::int(i64::from(self.multiplicity(index))))
                .collect()
        };
        Value::Record(RecordValue {
            fields: vec![
                (
                    String::from("variable"),
                    self.variable
                        .as_ref()
                        .map_or(Value::String(String::new()), |name| {
                            Value::Symbol(name.clone())
                        }),
                ),
                (
                    String::from("roots"),
                    Value::Tensor(TensorValue::vector(self.roots.clone())),
                ),
                (
                    String::from("multiplicities"),
                    Value::Tensor(TensorValue::vector(multiplicities)),
                ),
                (
                    String::from("conditions"),
                    Value::Tensor(TensorValue::vector(
                        self.conditions.iter().map(Condition::to_value).collect(),
                    )),
                ),
                (String::from("complete"), Value::Bool(self.complete)),
            ],
        })
    }

    /// Reads a record written by [`Solutions::to_value`].
    ///
    /// # Errors
    /// Fails when the record does not match [`Solutions::schema`].
    pub fn from_value(value: &Value) -> Result<Self, SchemaError> {
        let record = as_record("solutions", value)?;
        Self::schema().check(record)?;
        let variable = match record.get("variable") {
            | Some(Value::Symbol(name)) => Some(name.clone()),
            | Some(Value::String(name)) if !name.is_empty() => Some(name.clone()),
            | _ => None,
        };
        let roots = elements(record.get("roots"));
        let multiplicities: Vec<u32> = elements(record.get("multiplicities"))
            .iter()
            .map(|value| value.as_f64().map_or(1, |n| n.max(0.0) as u32))
            .collect();
        let conditions = elements(record.get("conditions"))
            .iter()
            .map(|value| {
                let record = as_record("condition", value)?;
                let text = match record.get("text") {
                    | Some(Value::String(text)) => text.clone(),
                    | _ => String::new(),
                };
                let kind = match record.get("kind") {
                    | Some(Value::String(kind)) => condition_kind(kind),
                    | _ => ConditionKind::Other,
                };
                Ok(Condition {
                    kind,
                    expression: record.get("expression").cloned(),
                    text,
                })
            })
            .collect::<Result<Vec<_>, SchemaError>>()?;
        let complete = matches!(record.get("complete"), Some(Value::Bool(true)));
        if !multiplicities.is_empty() && multiplicities.len() != roots.len() {
            return Err(SchemaError::Length {
                schema: String::from("solutions"),
                field: String::from("multiplicities"),
                expected: roots.len(),
                got: multiplicities.len(),
            });
        }
        Ok(Self {
            variable,
            roots,
            multiplicities,
            conditions,
            complete,
        })
    }
}

/// The `ConditionKind` named `kind`, defaulting to [`ConditionKind::Other`].
fn condition_kind(kind: &str) -> ConditionKind {
    match kind {
        | "nonzero" => ConditionKind::NonZero,
        | "positive" => ConditionKind::Positive,
        | "nonnegative" => ConditionKind::Nonnegative,
        | "negative" => ConditionKind::Negative,
        | "real" => ConditionKind::Real,
        | "branch" => ConditionKind::Branch,
        | "distinct" => ConditionKind::Distinct,
        | "interval" => ConditionKind::Interval,
        | _ => ConditionKind::Other,
    }
}

/// `value` as a record, or a schema error naming `schema`.
fn as_record<'v>(
    schema: &str,
    value: &'v Value,
) -> Result<&'v RecordValue, SchemaError> {
    match value {
        | Value::Record(record) => Ok(record),
        | _ => {
            Err(SchemaError::NotARecord {
                schema: schema.to_owned(),
            })
        },
    }
}

/// The elements of a list-shaped value, or nothing.
fn elements(value: Option<&Value>) -> Vec<Value> {
    match value {
        | Some(Value::Tensor(tensor)) => tensor.data.clone(),
        | Some(Value::Sparse(sparse)) => sparse.data.clone(),
        | _ => Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// ODE solutions
// ---------------------------------------------------------------------------

/// One event located during integration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EventRecord {
    /// The time it happened.
    pub t: f64,
    /// Which component triggered it, when the event named one.
    pub component: Option<String>,
    /// The state at the event.
    pub y: Vec<f64>,
}

/// An ODE solution as a structured result.
///
/// This is the ABI view of
/// [`OdeSolution`](crate::kernels::ode_adaptive::OdeSolution): the same data,
/// with the field layout made explicit so a host can build an array or a
/// table without guessing.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OdeSolutionView {
    /// Accepted time points.
    pub t: Vec<f64>,
    /// States at those points, one row per point.
    pub y: Vec<Vec<f64>>,
    /// What each column of `y` is, in the order the model declared them.
    pub flatten_order: Vec<String>,
    /// Located events.
    pub events: Vec<EventRecord>,
    /// Right-hand-side evaluations.
    pub nfev: usize,
    /// Rejected steps.
    pub rejected: usize,
    /// The dimension of each state component, when known.
    pub dimensions: Vec<Option<DimExpr>>,
}

impl OdeSolutionView {
    /// The layout every `OdeSolutionView` record has.
    #[must_use]
    pub fn schema() -> FlattenSchema {
        FlattenSchema::new(
            "ode_solution",
            1,
            vec![
                Field::new("t", FieldKind::List).doc("accepted time points"),
                Field::new("y", FieldKind::Tensor(Some(2)))
                    .doc("states, one row per accepted time point"),
                Field::new("flatten_order", FieldKind::List)
                    .doc("what each column of `y` is, in order"),
                Field::new("events", FieldKind::List).doc("located events, in time order"),
                Field::new("nfev", FieldKind::Scalar).doc("right-hand-side evaluations"),
                Field::new("rejected", FieldKind::Scalar).doc("rejected steps"),
                Field::new("dimensions", FieldKind::List)
                    .doc("the dimension of each component, by name; empty when unknown"),
            ],
        )
    }

    /// The number of state components, taken from the first row.
    #[must_use]
    pub fn components(&self) -> usize {
        self.y.first().map_or(0, Vec::len)
    }

    /// Whether every row has the same length as the first.
    #[must_use]
    pub fn is_rectangular(&self) -> bool {
        let width = self.components();
        self.y.iter().all(|row| row.len() == width)
    }

    /// The state as one row-major tensor of shape `[steps, components]`.
    ///
    /// # Errors
    /// Fails when the rows are ragged, which means the integrator reported a
    /// state of changing size.
    pub fn flatten(&self) -> Result<TensorValue, SchemaError> {
        if !self.is_rectangular() {
            return Err(SchemaError::Length {
                schema: String::from("ode_solution"),
                field: String::from("y"),
                expected: self.components(),
                got: self.y.iter().map(Vec::len).max().unwrap_or(0),
            });
        }
        let data = self
            .y
            .iter()
            .flat_map(|row| row.iter().map(|&x| Value::Float(x)))
            .collect();
        Ok(TensorValue {
            shape: vec![self.y.len(), self.components()],
            elem: crate::kernel::signature::TypeExpr::Float,
            data,
        })
    }

    /// The answer as a record.
    ///
    /// # Errors
    /// Fails when `y` is ragged; the record could not be described by the
    /// schema's rank-2 `y` field.
    pub fn to_value(&self) -> Result<Value, SchemaError> {
        let flat = self.flatten()?;
        let events = self
            .events
            .iter()
            .map(|event| {
                Value::Record(RecordValue {
                    fields: vec![
                        (String::from("t"), Value::Float(event.t)),
                        (
                            String::from("component"),
                            event
                                .component
                                .as_ref()
                                .map_or(Value::String(String::new()), |name| {
                                    Value::String(name.clone())
                                }),
                        ),
                        (
                            String::from("y"),
                            Value::Tensor(TensorValue::vector(
                                event.y.iter().copied().map(Value::Float).collect(),
                            )),
                        ),
                    ],
                })
            })
            .collect();
        Ok(Value::Record(RecordValue {
            fields: vec![
                (
                    String::from("t"),
                    Value::Tensor(TensorValue::vector(
                        self.t.iter().copied().map(Value::Float).collect(),
                    )),
                ),
                (String::from("y"), Value::Tensor(flat)),
                (
                    String::from("flatten_order"),
                    Value::Tensor(TensorValue::vector(
                        self.flatten_order
                            .iter()
                            .cloned()
                            .map(Value::String)
                            .collect(),
                    )),
                ),
                (
                    String::from("events"),
                    Value::Tensor(TensorValue::vector(events)),
                ),
                (
                    String::from("nfev"),
                    Value::int(i64::try_from(self.nfev).unwrap_or(i64::MAX)),
                ),
                (
                    String::from("rejected"),
                    Value::int(i64::try_from(self.rejected).unwrap_or(i64::MAX)),
                ),
                // Kept in the record, not only in the struct: a host that
                // reads the `Value` must see the dimensions, or it has to go
                // back to Rust for them.
                (
                    String::from("dimensions"),
                    Value::Tensor(TensorValue::vector(
                        self.dimensions
                            .iter()
                            .map(|dim| {
                                Value::String(
                                    dim.as_ref().map_or_else(
                                        String::new,
                                        crate::kernel::signature::dim_name,
                                    ),
                                )
                            })
                            .collect(),
                    )),
                ),
            ],
        }))
    }
}

// ---------------------------------------------------------------------------
// Simulations
// ---------------------------------------------------------------------------

/// One output channel of a simulation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SimField {
    /// The name the model declared.
    pub name: String,
    /// The column, in the same order as `SimResult::t`.
    pub values: Vec<f64>,
    /// Its dimension, when the model knows it.
    pub dim: Option<DimExpr>,
    /// Its unit, when it is a quantity rather than a bare number.
    pub unit: Option<QuantityValue>,
}

/// A simulation result with the metadata needed to reshape it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SimResult {
    /// Which scenario produced it.
    pub scenario: String,
    /// The time axis.
    pub t: Vec<f64>,
    /// The output channels, in declaration order.
    pub fields: Vec<SimField>,
}

impl SimResult {
    /// The layout every `SimResult` record has.
    #[must_use]
    pub fn schema() -> FlattenSchema {
        FlattenSchema::new(
            "sim_result",
            1,
            vec![
                Field::new("scenario", FieldKind::Text).doc("which scenario produced it"),
                Field::new("t", FieldKind::List).doc("the time axis"),
                Field::new("fields", FieldKind::List).doc("output channels, in declaration order"),
                Field::new("schema", FieldKind::Record(String::from("sim_field")))
                    .doc("the layout of `fields`"),
            ],
        )
    }

    /// The names of the output channels, in order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.fields
            .iter()
            .map(|field| field.name.as_str())
            .collect()
    }

    /// The result as one tensor of shape `[steps, channels]`.
    ///
    /// # Errors
    /// Fails when a channel's length disagrees with the time axis.
    pub fn flatten(&self) -> Result<TensorValue, SchemaError> {
        for field in &self.fields {
            if field.values.len() != self.t.len() {
                return Err(SchemaError::Length {
                    schema: String::from("sim_result"),
                    field: field.name.clone(),
                    expected: self.t.len(),
                    got: field.values.len(),
                });
            }
        }
        let mut data = Vec::with_capacity(self.t.len() * self.fields.len());
        for step in 0..self.t.len() {
            for field in &self.fields {
                data.push(Value::Float(
                    field.values.get(step).copied().unwrap_or(f64::NAN),
                ));
            }
        }
        Ok(TensorValue {
            shape: vec![self.t.len(), self.fields.len()],
            elem: crate::kernel::signature::TypeExpr::Float,
            data,
        })
    }

    /// The result as a record.
    ///
    /// # Errors
    /// Fails when a channel disagrees with the time axis.
    pub fn to_value(&self) -> Result<Value, SchemaError> {
        let flat = self.flatten()?;
        let layout = FlattenSchema::new(
            "sim_field",
            1,
            self.fields
                .iter()
                .map(|field| {
                    let mut declared =
                        Field::new(&field.name, FieldKind::List).doc("one value per time step");
                    if let Some(dim) = &field.dim {
                        declared = declared.with_dim(dim.clone());
                    }
                    declared
                })
                .collect(),
        );
        Ok(Value::Record(RecordValue {
            fields: vec![
                (
                    String::from("scenario"),
                    Value::String(self.scenario.clone()),
                ),
                (
                    String::from("t"),
                    Value::Tensor(TensorValue::vector(
                        self.t.iter().copied().map(Value::Float).collect(),
                    )),
                ),
                (String::from("fields"), Value::Tensor(flat)),
                (String::from("schema"), layout.to_value()),
            ],
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::codec;

    #[test]
    fn a_solution_set_round_trips_through_a_record() {
        let solutions = Solutions {
            variable: Some(String::from("x")),
            roots: vec![Value::int(1_i32), Value::int(-1_i32)],
            multiplicities: vec![2, 1],
            conditions: vec![Condition::new(ConditionKind::NonZero, "a != 0")],
            complete: true,
        };
        let value = solutions.to_value();
        let back = Solutions::from_value(&value).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, solutions);
        assert_eq!(back.multiplicity(0), 2);
        assert_eq!(back.multiplicity(1), 1);
    }

    #[test]
    fn a_record_that_does_not_match_its_schema_is_refused() {
        let bad = Value::Record(RecordValue {
            fields: vec![(String::from("roots"), Value::int(1_i32))],
        });
        // The first field the schema declares is the one it reports first.
        assert!(matches!(
            Solutions::from_value(&bad),
            Err(SchemaError::Missing { ref field, .. }) if field == "variable"
        ));
        // A field of the wrong shape is caught once the record is otherwise
        // well formed.
        let wrong_kind = Value::Record(RecordValue {
            fields: vec![
                (String::from("variable"), Value::String(String::new())),
                (String::from("roots"), Value::int(1_i32)),
                (
                    String::from("multiplicities"),
                    Value::Tensor(TensorValue::vector(Vec::new())),
                ),
                (
                    String::from("conditions"),
                    Value::Tensor(TensorValue::vector(Vec::new())),
                ),
                (String::from("complete"), Value::Bool(true)),
            ],
        });
        assert!(matches!(
            Solutions::from_value(&wrong_kind),
            Err(SchemaError::Kind { ref field, .. }) if field == "roots"
        ));
    }

    #[test]
    fn a_ragged_ode_solution_is_refused() {
        let view = OdeSolutionView {
            t: vec![0.0, 1.0],
            y: vec![vec![0.0, 1.0], vec![2.0]],
            flatten_order: vec![String::from("x"), String::from("v")],
            events: Vec::new(),
            nfev: 0,
            rejected: 0,
            dimensions: Vec::new(),
        };
        assert!(view.flatten().is_err());
    }

    #[test]
    fn an_ode_solution_flattens_in_the_declared_order() {
        let view = OdeSolutionView {
            t: vec![0.0, 1.0],
            y: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            flatten_order: vec![String::from("x"), String::from("v")],
            events: vec![EventRecord {
                t: 0.5,
                component: Some(String::from("x")),
                y: vec![2.0, 2.5],
            }],
            nfev: 12,
            rejected: 1,
            dimensions: Vec::new(),
        };
        let flat = view.flatten().unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(flat.shape, vec![2, 2]);
        assert_eq!(flat.data[0], Value::Float(1.0));
        assert_eq!(flat.data[3], Value::Float(4.0));
        let value = view.to_value().unwrap_or_else(|e| panic!("{e}"));
        // The record is a `Value`, so it goes through the codec like any
        // other.
        let bytes = value.to_cbor().unwrap_or_else(|e| panic!("{e}"));
        let back = Value::from_cbor(&bytes).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, value);
    }

    #[test]
    fn a_simulation_result_flattens_column_major() {
        let result = SimResult {
            scenario: String::from("bounce"),
            t: vec![0.0, 0.1],
            fields: vec![
                SimField {
                    name: String::from("h"),
                    values: vec![1.0, 0.9],
                    dim: None,
                    unit: None,
                },
                SimField {
                    name: String::from("v"),
                    values: vec![0.0, -1.0],
                    dim: None,
                    unit: None,
                },
            ],
        };
        let flat = result.flatten().unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(flat.shape, vec![2, 2]);
        assert_eq!(flat.data[1], Value::Float(0.0));
        assert_eq!(flat.data[2], Value::Float(0.9));
        assert_eq!(result.names(), vec!["h", "v"]);
    }

    #[test]
    fn the_schemas_are_themselves_encodable() {
        for schema in [
            Solutions::schema(),
            OdeSolutionView::schema(),
            SimResult::schema(),
        ] {
            let value = schema.to_value();
            let bytes = codec::encode_deterministic(&value).unwrap_or_else(|e| panic!("{e}"));
            let back: Value = codec::decode_deterministic(&bytes).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(back, value);
            assert!(!schema.is_empty());
        }
    }
}
