//! Sorts: what an operator was *seen* to take and return.
//!
//! # The problem
//!
//! A computer-algebra kernel has well over a thousand operators and almost
//! none of them can be given a type by looking at its registration: `det`,
//! `factorint` and `groebner` all register as "a heavy operator of arity `n`",
//! so the derivation in [`crate::kernel::infer`] can only state the symbolic
//! contract `(Term, ..) -> Term`. A host that generates typed bindings from
//! that learns nothing, and a human writing a thousand signatures by hand is
//! exactly the drift the manifest exists to prevent.
//!
//! # The idea: the engine already proves the types
//!
//! Every run of the engine ends with an e-graph in which a request such as
//! `factorint(360)` sits in the same equivalence class as its answer
//! `list(list(2, 3), list(3, 2), list(5, 1))`. That class *is* a typing
//! judgement: this operator, applied to an integer, is equal to a rectangular
//! list of lists of integers. Nothing has to be guessed; it only has to be
//! written down. This module
//!
//! 1. defines a small lattice of [`Sort`]s (the shapes a term can have:
//!    `Int`, `Expr`, a list of sorts, a tuple, ...) with a join;
//! 2. classifies any e-class into a sort ([`sort_of`]) using only what the
//!    graph and the operator table state (literals, list constructors, the
//!    `PREDICATE` flag, the dimension rule that says "not a scalar");
//! 3. records, when tracing is on, one [`Fact`] per operator application and
//!    distinct sort tuple ([`observe`], called at the end of
//!    [`Engine::run`](crate::graph::Engine::run));
//! 4. turns facts into signatures ([`mine`]). The trace is only a source of
//!    *calls*: each is run again on the full engine and typed from that run
//!    ([`observe_call`]), where an application is never evidence about
//!    itself. The answered calls are clustered by the kinds of their
//!    arguments and each cluster is joined (*anti-unification*). Every
//!    candidate is then **probed** ([`probe`]): the same call with shifted
//!    numbers must give the same sort, or the result is widened, so that
//!    `sin(0) = 0` does not type `sin` as integer-valued.
//!
//! The result, `signatures.witnessed.json`, is a third source of signatures
//! between the hand-written overlay and the derivation:
//! [`SignatureSource::Witnessed`]. A signature whose parameters differ from
//! every registered overload becomes an overload; one whose parameters equal
//! a registered overload's is kept as a *refinement* of that overload's
//! result. Each carries the calls that justify it, so the manifest can show
//! them and a test can re-run them. The corpus is the kernel's own test
//! suite: every operator that is tested is thereby typed.
//!
//! # What a witnessed signature promises
//!
//! It is *evidence*, not a proof of totality: `(Int) -> List<Int>` for
//! `divisors` says that integer arguments were seen to produce integer lists,
//! not that no integer can ever produce something else. The parameter side
//! may be narrower than what the operator accepts (the symbolic contract
//! remains as the most general overload, so nothing becomes uncallable) and a
//! host must check the result it decodes. That is the contract of a gradual
//! type at a dynamic boundary, and it is the honest one.
//!
//! Overload sets stay unambiguous by construction: [`mine`] only keeps a
//! signature when it is comparable with, or disjoint from, every signature
//! already kept (see [`coherent`]).
//!
//! [`SignatureSource::Witnessed`]: crate::kernel::manifest::SignatureSource::Witnessed

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    reason = "a classifier over small vectors whose bounds are checked where they are built"
)]

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::fmt;
use std::io::Write as _;
use std::sync::Mutex;
use std::sync::OnceLock;

use serde::Deserialize;
use serde::Serialize;

use crate::graph::Arity;
use crate::graph::ClassId;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Payload;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::Variadic;

/// The environment variable that turns tracing on; its value is the path of
/// the fact log to append to.
pub const TRACE_ENV: &str = "RSSN_SORT_TRACE";

/// The shape of a term.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Sort {
    /// No value (the element sort of an empty list). Bottom.
    Never,
    /// `true` or `false`.
    Bool,
    /// An exact integer literal.
    Int,
    /// An exact rational literal.
    Rational,
    /// A floating-point literal.
    Float,
    /// A bare symbol.
    Symbol,
    /// A string payload.
    Text,
    /// A scalar symbolic expression.
    Expr,
    /// A truth-valued formula that is not a constant (an equation, a
    /// predicate application).
    Prop,
    /// A structured value headed by the named constructor (`interval`,
    /// `set`, ...): an operator that is neither heavy nor scalar-valued.
    Object(String),
    /// A list whose elements all have the sort `elem`.
    List {
        /// The join of the element sorts.
        elem: Box<Self>,
        /// When the elements are lists: whether they all have the same
        /// length (so the value is a rectangular array).
        rect: bool,
    },
    /// A list of fixed length whose elements have unrelated sorts.
    Tuple(Vec<Self>),
    /// Anything. Top.
    Any,
}

impl Sort {
    /// A rectangular list of `elem`.
    #[must_use]
    pub fn list(elem: Self) -> Self {
        Self::List {
            elem: Box::new(elem),
            rect: true,
        }
    }

    /// Whether this is one of the scalar sorts below [`Sort::Expr`].
    #[must_use]
    pub const fn is_scalar(&self) -> bool {
        matches!(
            self,
            Self::Int | Self::Rational | Self::Float | Self::Symbol | Self::Expr
        )
    }

    /// The least upper bound of two sorts.
    #[must_use]
    pub fn join(
        &self,
        other: &Self,
    ) -> Self {
        use Sort::Any;
        use Sort::Bool;
        use Sort::Expr;
        use Sort::Int;
        use Sort::List;
        use Sort::Never;
        use Sort::Prop;
        use Sort::Rational;
        use Sort::Tuple;
        if self == other {
            return self.clone();
        }
        match (self, other) {
            | (Never, s) | (s, Never) => s.clone(),
            | (Any, _) | (_, Any) => Any,
            | (Int, Rational) | (Rational, Int) => Rational,
            | (a, b) if a.is_scalar() && b.is_scalar() => Expr,
            | (Bool, Prop) | (Prop, Bool) => Prop,
            | (List { elem: a, rect: ra }, List { elem: b, rect: rb }) => {
                List {
                    elem: Box::new(a.join(b)),
                    rect: *ra && *rb,
                }
            },
            | (Tuple(a), Tuple(b)) if a.len() == b.len() => {
                Tuple(a.iter().zip(b).map(|(x, y)| x.join(y)).collect())
            },
            | (Tuple(items), l @ List { .. }) | (l @ List { .. }, Tuple(items)) => {
                tuple_as_list(items).join(l)
            },
            | (Tuple(a), Tuple(b)) => tuple_as_list(a).join(&tuple_as_list(b)),
            | _ => Any,
        }
    }

    /// Whether every value of `self` is a value of `other`.
    #[must_use]
    pub fn le(
        &self,
        other: &Self,
    ) -> bool {
        &self.join(other) == other
    }

    /// Nesting depth of rectangular lists and the scalar at the bottom, when
    /// the sort is such an array (`[[Int]]` is `(2, Int)`).
    #[must_use]
    pub fn array(&self) -> Option<(u64, &Self)> {
        match self {
            | Self::List { elem, rect } => {
                match &**elem {
                    | inner @ Self::List { .. } if *rect => {
                        inner.array().map(|(d, s)| (d.saturating_add(1), s))
                    },
                    | Self::List { .. } => None,
                    | scalar => Some((1, scalar)),
                }
            },
            | _ => None,
        }
    }

    /// The type a host sees for this sort.
    ///
    /// Scalars map to the scalar types, rectangular lists to tensors of the
    /// matching rank, ragged lists to a rank-1 tensor of their element type,
    /// a tuple to a positional record (`_0`, `_1`, ...), and everything else
    /// to a symbolic term.
    #[must_use]
    pub fn type_expr(&self) -> TypeExpr {
        match self {
            | Self::Bool => TypeExpr::Bool,
            | Self::Int => TypeExpr::Int,
            | Self::Rational => TypeExpr::Real,
            | Self::Float => TypeExpr::Float,
            | Self::Symbol => TypeExpr::Symbol,
            | Self::Text => TypeExpr::String,
            | Self::Object(name) => TypeExpr::Opaque(name.clone()),
            | Self::Tuple(items) => {
                TypeExpr::Record(
                    items
                        .iter()
                        .enumerate()
                        .map(|(k, s)| (format!("_{k}"), s.type_expr()))
                        .collect(),
                )
            },
            | Self::List { elem, .. } => {
                match self.array() {
                    | Some((rank, scalar)) => {
                        TypeExpr::Tensor {
                            rank: NatExpr::Lit(rank),
                            shape: None,
                            elem: Box::new(scalar.type_expr()),
                        }
                    },
                    | None => {
                        TypeExpr::Tensor {
                            rank: NatExpr::Lit(1),
                            shape: None,
                            elem: Box::new(elem.type_expr()),
                        }
                    },
                }
            },
            | Self::Never | Self::Expr | Self::Prop | Self::Any => crate::kernel::infer::term(),
        }
    }

    /// Whether a signature mentioning this sort says more than "a term".
    #[must_use]
    pub const fn informative(&self) -> bool {
        !matches!(self, Self::Never | Self::Expr | Self::Prop | Self::Any)
    }

    /// A coarse key: facts with equal kinds in every position belong to one
    /// candidate overload.
    fn kind(&self) -> String {
        match self {
            | Self::List { .. } => {
                match self.array() {
                    | Some((d, s)) => format!("l{d}{}", s.kind()),
                    | None => "lr".to_owned(),
                }
            },
            | Self::Tuple(items) => format!("u{}", items.len()),
            | Self::Symbol | Self::Expr | Self::Never => "e".to_owned(),
            | other => other.to_string(),
        }
    }

    /// Parses the text form written by [`fmt::Display`].
    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        let mut p = Parser {
            src: text.as_bytes(),
            at: 0,
        };
        let sort = p.sort()?;
        (p.at == p.src.len()).then_some(sort)
    }
}

fn tuple_as_list(items: &[Sort]) -> Sort {
    let elem = items.iter().fold(Sort::Never, |acc, s| acc.join(s));
    Sort::List {
        elem: Box::new(elem),
        rect: true,
    }
}

impl fmt::Display for Sort {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Never => f.write_str("Never"),
            | Self::Bool => f.write_str("Bool"),
            | Self::Int => f.write_str("Int"),
            | Self::Rational => f.write_str("Rat"),
            | Self::Float => f.write_str("Float"),
            | Self::Symbol => f.write_str("Sym"),
            | Self::Text => f.write_str("Text"),
            | Self::Expr => f.write_str("Expr"),
            | Self::Prop => f.write_str("Prop"),
            | Self::Any => f.write_str("Any"),
            | Self::Object(name) => write!(f, "Obj<{name}>"),
            | Self::List { elem, rect } => {
                write!(f, "[{elem}]{}", if *rect { "" } else { "~" })
            },
            | Self::Tuple(items) => {
                f.write_str("(")?;
                for (k, s) in items.iter().enumerate() {
                    if k > 0 {
                        f.write_str(",")?;
                    }
                    write!(f, "{s}")?;
                }
                f.write_str(")")
            },
        }
    }
}

struct Parser<'a> {
    src: &'a [u8],
    at: usize,
}

impl Parser<'_> {
    fn eat(
        &mut self,
        c: u8,
    ) -> bool {
        if self.src.get(self.at) == Some(&c) {
            self.at += 1;
            true
        } else {
            false
        }
    }

    fn word(&mut self) -> &str {
        let start = self.at;
        while self
            .src
            .get(self.at)
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
        {
            self.at += 1;
        }
        std::str::from_utf8(&self.src[start..self.at]).unwrap_or("")
    }

    fn sort(&mut self) -> Option<Sort> {
        if self.eat(b'[') {
            let elem = self.sort()?;
            if !self.eat(b']') {
                return None;
            }
            let rect = !self.eat(b'~');
            return Some(Sort::List {
                elem: Box::new(elem),
                rect,
            });
        }
        if self.eat(b'(') {
            let mut items = Vec::new();
            if !self.eat(b')') {
                loop {
                    items.push(self.sort()?);
                    if self.eat(b')') {
                        break;
                    }
                    if !self.eat(b',') {
                        return None;
                    }
                }
            }
            return Some(Sort::Tuple(items));
        }
        let word = self.word().to_owned();
        Some(match word.as_str() {
            | "Never" => Sort::Never,
            | "Bool" => Sort::Bool,
            | "Int" => Sort::Int,
            | "Rat" => Sort::Rational,
            | "Float" => Sort::Float,
            | "Sym" => Sort::Symbol,
            | "Text" => Sort::Text,
            | "Expr" => Sort::Expr,
            | "Prop" => Sort::Prop,
            | "Any" => Sort::Any,
            | "Obj" => {
                if !self.eat(b'<') {
                    return None;
                }
                let name = self.word().to_owned();
                if !self.eat(b'>') {
                    return None;
                }
                Sort::Object(name)
            },
            | _ => return None,
        })
    }
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

/// How deep nested lists are followed before the answer is [`Sort::Any`].
const MAX_DEPTH: usize = 8;

/// The sort of an e-class, or `None` when the class holds nothing but
/// unreduced heavy operators (a request that was not answered says nothing
/// about its own type).
#[must_use]
pub fn sort_of(
    graph: &Graph,
    class: ClassId,
) -> Option<Sort> {
    let mut memo = HashMap::new();
    classify(
        graph,
        graph.find(NodeId::from_raw(class.raw())),
        &mut memo,
        0,
    )
}

/// The sort of the class of `node`.
#[must_use]
pub fn sort_of_node(
    graph: &Graph,
    node: NodeId,
) -> Option<Sort> {
    let mut memo = HashMap::new();
    classify(graph, graph.find(node), &mut memo, 0)
}

fn classify(
    graph: &Graph,
    class: ClassId,
    memo: &mut HashMap<ClassId, Option<Sort>>,
    depth: usize,
) -> Option<Sort> {
    if let Some(known) = memo.get(&class) {
        return known.clone();
    }
    if depth > MAX_DEPTH {
        return Some(Sort::Any);
    }
    // Guard against a list that contains its own class.
    memo.insert(class, Some(Sort::Any));
    let sort = classify_fresh(graph, class, memo, depth);
    memo.insert(class, sort.clone());
    sort
}

fn classify_fresh(
    graph: &Graph,
    class: ClassId,
    memo: &mut HashMap<ClassId, Option<Sort>>,
    depth: usize,
) -> Option<Sort> {
    classify_without(graph, class, None, memo, depth)
}

/// The sort of a class judged by every member except applications of
/// `skip`.
///
/// This is how the *answer* of an operator is typed: `gcd(12, 18)` sitting
/// alone in its class (because the rule set that reduces it is not loaded)
/// says nothing about what `gcd` returns, while the literal `6` next to it
/// does. An application is never evidence about itself.
fn classify_without(
    graph: &Graph,
    class: ClassId,
    skip: Option<OpId>,
    memo: &mut HashMap<ClassId, Option<Sort>>,
    depth: usize,
) -> Option<Sort> {
    let mut list: Option<NodeId> = None;
    let mut prop = false;
    let mut object: Option<String> = None;
    let mut scalar = false;
    for node in graph.enodes(class) {
        let op = graph.op(node);
        if Some(op) == skip {
            continue;
        }
        if let Some(leaf) = leaf_sort(graph, node, op) {
            return Some(leaf);
        }
        if op == core::LIST {
            list = Some(list.map_or(node, |best| best.min(node)));
            continue;
        }
        let desc = graph.ops().get(op);
        if desc.flags.has(OpFlags::HEAVY) {
            continue;
        }
        if op == core::EQ || desc.flags.has(OpFlags::PREDICATE) || is_connective(&desc.name) {
            prop = true;
        } else if structured(graph, op) {
            let name = desc.name.to_string();
            object = Some(match object {
                | Some(best) => best.min(name),
                | None => name,
            });
        } else {
            scalar = true;
        }
    }
    if let Some(node) = list {
        return Some(list_sort(graph, node, memo, depth));
    }
    if scalar {
        return Some(Sort::Expr);
    }
    if prop {
        return Some(Sort::Prop);
    }
    object.map(Sort::Object)
}

fn leaf_sort(
    graph: &Graph,
    node: NodeId,
    op: OpId,
) -> Option<Sort> {
    match graph.payload(node) {
        | Some(Payload::Num(number)) if op == core::LIT => {
            return Some(match number {
                | Number::Int(_) => Sort::Int,
                | Number::Rat(_) => Sort::Rational,
                | Number::Float(_) => Sort::Float,
            });
        },
        | Some(Payload::Sym(_)) if op == core::SYM => return Some(Sort::Symbol),
        | Some(Payload::Bool(_)) => return Some(Sort::Bool),
        | Some(Payload::Str(_)) => return Some(Sort::Text),
        | Some(Payload::Unit(_)) => return Some(Sort::Object("unit".to_owned())),
        | _ => {},
    }
    let desc = graph.ops().get(op);
    if matches!(desc.arity, Arity::Fixed(0)) && graph.children(node).is_empty() {
        return Some(match &*desc.name {
            | "true" | "false" => Sort::Bool,
            | _ if desc.flags.has(OpFlags::LEAF) => Sort::Object(desc.name.to_string()),
            | _ => Sort::Expr,
        });
    }
    None
}

fn is_connective(name: &str) -> bool {
    matches!(name, "and" | "or" | "not" | "xor" | "implies" | "iff")
}

/// Whether the operator's dimension rule says its value is not a scalar.
fn structured(
    graph: &Graph,
    op: OpId,
) -> bool {
    matches!(
        graph.ops().attr::<DimRule>(op),
        Some(DimRule::NonMetric | DimRule::Container)
    )
}

fn list_sort(
    graph: &Graph,
    node: NodeId,
    memo: &mut HashMap<ClassId, Option<Sort>>,
    depth: usize,
) -> Sort {
    let children = graph.children(node).to_vec();
    if children.is_empty() {
        return Sort::list(Sort::Never);
    }
    let mut sorts = Vec::with_capacity(children.len());
    let mut lengths: Vec<Option<usize>> = Vec::with_capacity(children.len());
    for child in children {
        let class = graph.find(child);
        // An element that is an unanswered request is still an expression.
        sorts.push(classify(graph, class, memo, depth + 1).unwrap_or(Sort::Expr));
        lengths.push(
            graph
                .enodes(class)
                .filter(|&n| graph.op(n) == core::LIST)
                .min()
                .map(|n| graph.children(n).len()),
        );
    }
    let joined = sorts.iter().fold(Sort::Never, |acc, s| acc.join(s));
    if joined == Sort::Any || matches!(joined, Sort::List { ref elem, .. } if **elem == Sort::Any) {
        // Unrelated element sorts: a tuple, unless it is too long to be one.
        if sorts.len() <= 6 && sorts.iter().any(|s| *s != Sort::Any) {
            return Sort::Tuple(sorts);
        }
        return Sort::list(Sort::Any);
    }
    let rect = match &joined {
        | Sort::List { .. } => {
            let first = lengths.first().copied().flatten();
            first.is_some() && lengths.iter().all(|l| *l == first)
        },
        | _ => true,
    };
    Sort::List {
        elem: Box::new(joined),
        rect,
    }
}

// ---------------------------------------------------------------------------
// Facts and tracing
// ---------------------------------------------------------------------------

/// One observed application: an operator, the sorts of its arguments, the
/// sort of the class it ended in, and the call that showed it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fact {
    /// Operator name.
    pub op: String,
    /// Whether the operator is variadic.
    pub variadic: bool,
    /// Argument sorts.
    pub params: Vec<Sort>,
    /// Result sort.
    pub result: Sort,
    /// The application, in surface syntax.
    pub witness: String,
}

impl Fact {
    /// The tab-separated log line of this fact.
    #[must_use]
    pub fn line(&self) -> String {
        let params: Vec<String> = self.params.iter().map(ToString::to_string).collect();
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.op,
            if self.variadic { "v" } else { "f" },
            params.join(";"),
            self.result,
            self.witness
        )
    }

    /// Parses a log line.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let mut fields = line.splitn(5, '\t');
        let op = fields.next()?.to_owned();
        let variadic = fields.next()? == "v";
        let params_text = fields.next()?;
        let params = if params_text.is_empty() {
            Vec::new()
        } else {
            params_text
                .split(';')
                .map(Sort::parse)
                .collect::<Option<Vec<_>>>()?
        };
        let result = Sort::parse(fields.next()?)?;
        let witness = fields.next()?.to_owned();
        Some(Self {
            op,
            variadic,
            params,
            result,
            witness,
        })
    }
}

/// Longest witness kept; longer calls make poor examples and slow replays.
const MAX_WITNESS: usize = 240;

/// Collects the facts a finished graph shows.
///
/// One fact per operator application whose arguments have a sort and whose
/// class holds an answer (something other than the application itself). Core
/// constructors (`+`, `*`, `^`, lists, literals) are skipped:
/// their types are fixed by the language, not by a rule set.
#[must_use]
pub fn facts_of(graph: &Graph) -> Vec<Fact> {
    let mut memo: HashMap<ClassId, Option<Sort>> = HashMap::new();
    let mut out = Vec::new();
    for (op, desc) in graph.ops().iter() {
        if is_core(op) || desc.flags.has(OpFlags::LEAF) || matches!(desc.arity, Arity::Fixed(0)) {
            continue;
        }
        for &node in graph.nodes_with_op(op) {
            if let Some(fact) = fact_of(graph, node, desc, &mut memo) {
                out.push(fact);
            }
        }
    }
    out
}

fn is_core(op: OpId) -> bool {
    [
        core::LIT,
        core::SYM,
        core::ADD,
        core::MUL,
        core::POW,
        core::APPLY,
        core::LIST,
        core::EQ,
    ]
    .contains(&op)
}

fn fact_of(
    graph: &Graph,
    node: NodeId,
    desc: &OpDescriptor,
    memo: &mut HashMap<ClassId, Option<Sort>>,
) -> Option<Fact> {
    let result = classify_without(graph, graph.find(node), Some(graph.op(node)), memo, 0)?;
    let mut params = Vec::new();
    for &child in graph.children(node) {
        params.push(classify(graph, graph.find(child), memo, 0)?);
    }
    let witness = graph.display(node);
    if witness.len() > MAX_WITNESS || witness.contains('\t') || witness.contains('\n') {
        return None;
    }
    Some(Fact {
        op: desc.name.to_string(),
        variadic: matches!(desc.arity, Arity::Variadic),
        params,
        result,
        witness,
    })
}

struct Recorder {
    file: std::fs::File,
    /// Facts already written, with the length of the witness they carry.
    seen: HashMap<(String, Vec<Sort>, Sort), usize>,
}

fn recorder() -> Option<&'static Mutex<Recorder>> {
    static RECORDER: OnceLock<Option<Mutex<Recorder>>> = OnceLock::new();
    RECORDER
        .get_or_init(|| {
            let path = std::env::var_os(TRACE_ENV)?;
            let file = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .ok()?;
            Some(Mutex::new(Recorder {
                file,
                seen: HashMap::new(),
            }))
        })
        .as_ref()
}

/// Whether tracing is on in this process.
#[must_use]
pub fn tracing() -> bool {
    recorder().is_some()
}

/// Appends the facts of `graph` to the trace log, when tracing is on.
///
/// Called by the engine at the end of every top-level run. With tracing off
/// this is one atomic load.
pub fn observe(graph: &Graph) {
    let Some(recorder) = recorder() else {
        return;
    };
    let facts = facts_of(graph);
    let Ok(mut guard) = recorder.lock() else {
        return;
    };
    let mut text = String::new();
    for fact in facts {
        let key = (fact.op.clone(), fact.params.clone(), fact.result.clone());
        // A fact is written once, and again only for a much shorter witness.
        let better = guard
            .seen
            .get(&key)
            .is_none_or(|&len| fact.witness.len().saturating_mul(2) < len);
        if better {
            guard.seen.insert(key, fact.witness.len());
            text.push_str(&fact.line());
            text.push('\n');
        }
    }
    if !text.is_empty() {
        let _ = guard.file.write_all(text.as_bytes());
        let _ = guard.file.flush();
    }
}

// ---------------------------------------------------------------------------
// Mining: facts -> signatures
// ---------------------------------------------------------------------------

/// A call that justifies a signature.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Witness {
    /// The call, in surface syntax.
    pub input: String,
    /// The sort of its answer, in the text form of [`Sort`].
    pub output: String,
}

/// A signature with its evidence.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Witnessed {
    /// The signature.
    pub signature: Signature,
    /// Parameter sorts, in the text form of [`Sort`] (finer than the types).
    pub params: Vec<String>,
    /// Result sort.
    pub result: String,
    /// Number of distinct facts that were joined into it.
    pub support: u32,
    /// Calls that reproduce it.
    pub witnesses: Vec<Witness>,
    /// Whether the parameter list is one the registration already states
    /// (usually the symbolic contract), so that this entry *refines the
    /// result* of an existing overload instead of adding one. A refinement
    /// is not an overload: two overloads with one parameter list would be
    /// ambiguous. The manifest lists refinements separately.
    #[serde(default)]
    pub refines: bool,
}

/// The content of `signatures.witnessed.json`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct WitnessedTable {
    /// Format version.
    pub version: u32,
    /// Signatures by operator name.
    pub ops: BTreeMap<String, Vec<Witnessed>>,
}

/// Most witnessed overloads kept per operator.
pub const MAX_OVERLOADS: usize = 4;

/// Most witnesses stored per signature.
pub const MAX_WITNESSES: usize = 2;

/// A candidate overload: joined parameter sorts, joined result, evidence.
#[derive(Clone, Debug)]
pub struct Candidate {
    /// Operator name.
    pub op: String,
    /// Whether the operator is variadic (then `params` has one element, the
    /// join of every argument).
    pub variadic: bool,
    /// Joined parameter sorts.
    pub params: Vec<Sort>,
    /// Joined result sort.
    pub result: Sort,
    /// The facts joined, shortest witness first.
    pub facts: Vec<Fact>,
}

impl Candidate {
    /// The signature this candidate states.
    #[must_use]
    pub fn signature(&self) -> Signature {
        if self.variadic {
            let elem = self
                .params
                .first()
                .map_or_else(crate::kernel::infer::term, Sort::type_expr);
            Signature::new(Vec::new(), self.result.type_expr()).variadic(Variadic {
                min: 0,
                max: None,
                elem,
            })
        } else {
            Signature::new(
                self.params.iter().map(Sort::type_expr).collect(),
                self.result.type_expr(),
            )
        }
    }

    /// Whether the candidate says something about the *answer*. A signature
    /// that only narrows its parameters (`(Int) -> Term`) tells a host
    /// nothing the symbolic contract does not, and would only crowd the
    /// overload set.
    const fn informative(&self) -> bool {
        self.result.informative()
    }
}

/// Clusters facts into candidate overloads, most supported first per
/// operator.
///
/// `binder` gives, per operator name, the index of its bound-variable
/// argument: only there does a bare symbol stay a `Symbol` parameter.
/// Elsewhere "every test passed `x`" is not evidence that the operator wants
/// a symbol, so the parameter is widened to an expression.
#[must_use]
pub fn candidates(
    facts: &[Fact],
    binder: &dyn Fn(&str) -> Option<usize>,
) -> Vec<Candidate> {
    let mut clusters: BTreeMap<(String, String), Candidate> = BTreeMap::new();
    for fact in facts {
        let bound = binder(&fact.op);
        let params: Vec<Sort> = if fact.variadic {
            vec![
                fact.params
                    .iter()
                    .fold(Sort::Never, |acc, s| acc.join(&widen_symbol(s))),
            ]
        } else {
            fact.params
                .iter()
                .enumerate()
                .map(|(k, s)| {
                    if Some(k) == bound {
                        s.clone()
                    } else {
                        widen_symbol(s)
                    }
                })
                .collect()
        };
        let kinds: Vec<String> = params.iter().map(Sort::kind).collect();
        let key = (
            fact.op.clone(),
            format!("{}|{}", fact.params.len(), kinds.join(",")),
        );
        let key = if fact.variadic {
            (key.0, format!("v|{}", kinds.join(",")))
        } else {
            key
        };
        match clusters.get_mut(&key) {
            | Some(c) => {
                for (slot, s) in c.params.iter_mut().zip(&params) {
                    *slot = slot.join(s);
                }
                c.result = c.result.join(&fact.result);
                c.facts.push(fact.clone());
            },
            | None => {
                clusters.insert(
                    key,
                    Candidate {
                        op: fact.op.clone(),
                        variadic: fact.variadic,
                        params,
                        result: fact.result.clone(),
                        facts: vec![fact.clone()],
                    },
                );
            },
        }
    }
    let mut out: Vec<Candidate> = clusters
        .into_values()
        .filter(Candidate::informative)
        .collect();
    for c in &mut out {
        // The best example shows the whole result sort, then is short.
        let result = c.result.clone();
        c.facts.sort_by(|a, b| {
            (a.result != result, a.witness.len(), &a.witness).cmp(&(
                b.result != result,
                b.witness.len(),
                &b.witness,
            ))
        });
    }
    out.sort_by(|a, b| {
        (&a.op, std::cmp::Reverse(a.facts.len()), &a.params).cmp(&(
            &b.op,
            std::cmp::Reverse(b.facts.len()),
            &b.params,
        ))
    });
    out
}

fn widen_symbol(sort: &Sort) -> Sort {
    match sort {
        | Sort::Symbol => Sort::Expr,
        | other => other.clone(),
    }
}

/// How two parameter lists relate.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Relation {
    /// The same list.
    Equal,
    /// Every parameter of the first is below the second's.
    Below,
    /// Every parameter of the first is above the second's.
    Above,
    /// No argument list is accepted by both.
    Disjoint,
    /// Some argument list is accepted by both and neither is more specific.
    Overlapping,
}

/// How `a` and `b` relate as parameter lists of overloads of one operator.
#[must_use]
pub fn relation(
    a: &[TypeExpr],
    b: &[TypeExpr],
) -> Relation {
    if a.len() != b.len() {
        return Relation::Disjoint;
    }
    if a == b {
        return Relation::Equal;
    }
    let (mut below, mut above) = (true, true);
    for (x, y) in a.iter().zip(b) {
        match type_order(x, y) {
            | TypeOrder::Equal => {},
            | TypeOrder::Below => above = false,
            | TypeOrder::Above => below = false,
            | TypeOrder::Disjoint => return Relation::Disjoint,
            | TypeOrder::Overlapping => {
                below = false;
                above = false;
            },
        }
    }
    match (below, above) {
        | (true, _) => Relation::Below,
        | (_, true) => Relation::Above,
        | _ => Relation::Overlapping,
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum TypeOrder {
    Equal,
    Below,
    Above,
    Disjoint,
    Overlapping,
}

/// Rank of a scalar type in the widening chain, when it is in it.
const fn scalar_rank(t: &TypeExpr) -> Option<u8> {
    match t {
        | TypeExpr::Int => Some(0),
        | TypeExpr::Real => Some(1),
        | _ => None,
    }
}

fn type_order(
    a: &TypeExpr,
    b: &TypeExpr,
) -> TypeOrder {
    if a == b {
        return TypeOrder::Equal;
    }
    let generic = |t: &TypeExpr| matches!(t, TypeExpr::Dyn(_) | TypeExpr::Var(_));
    match (generic(a), generic(b)) {
        | (true, true) => return TypeOrder::Equal,
        | (false, true) => return TypeOrder::Below,
        | (true, false) => return TypeOrder::Above,
        | (false, false) => {},
    }
    if let (Some(x), Some(y)) = (scalar_rank(a), scalar_rank(b)) {
        return if x < y {
            TypeOrder::Below
        } else {
            TypeOrder::Above
        };
    }
    match (a, b) {
        | (
            TypeExpr::Tensor {
                rank: ra, elem: ea, ..
            },
            TypeExpr::Tensor {
                rank: rb, elem: eb, ..
            },
        ) => {
            if ra != rb {
                // A rank-1 tensor of terms accepts a matrix as a list of rows.
                return match (generic(ea), generic(eb)) {
                    | (true, _) | (_, true) => TypeOrder::Overlapping,
                    | _ => TypeOrder::Disjoint,
                };
            }
            type_order(ea, eb)
        },
        | (TypeExpr::Record(fa), TypeExpr::Record(fb)) if fa.len() == fb.len() => {
            let (mut below, mut above) = (true, true);
            for ((_, x), (_, y)) in fa.iter().zip(fb) {
                match type_order(x, y) {
                    | TypeOrder::Equal => {},
                    | TypeOrder::Below => above = false,
                    | TypeOrder::Above => below = false,
                    | TypeOrder::Disjoint => return TypeOrder::Disjoint,
                    | TypeOrder::Overlapping => {
                        below = false;
                        above = false;
                    },
                }
            }
            match (below, above) {
                | (true, _) => TypeOrder::Below,
                | (_, true) => TypeOrder::Above,
                | _ => TypeOrder::Overlapping,
            }
        },
        | _ => TypeOrder::Disjoint,
    }
}

/// Whether adding `new` to `kept` leaves the overload set unambiguous: `new`
/// must be comparable with, or disjoint from, each kept parameter list, and
/// must not repeat one.
#[must_use]
pub fn coherent(
    kept: &[Signature],
    new: &Signature,
) -> bool {
    kept.iter().all(|old| {
        if old.variadic.is_some() || new.variadic.is_some() {
            // A variadic overload is compared on its element type alone.
            let elem = |s: &Signature| {
                s.variadic
                    .as_ref()
                    .map_or_else(|| s.params.clone(), |v| vec![v.elem.clone()])
            };
            return match (old.variadic.is_some(), new.variadic.is_some()) {
                | (true, true) => {
                    matches!(
                        relation(&elem(old), &elem(new)),
                        Relation::Below | Relation::Above | Relation::Disjoint
                    )
                },
                | _ => false,
            };
        }
        matches!(
            relation(&old.params, &new.params),
            Relation::Below | Relation::Above | Relation::Disjoint
        )
    })
}

/// What [`replay`] found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Replay {
    /// The call ran and its sorts are within the candidate's.
    Confirmed(Sort),
    /// The call ran but an argument or the result has another sort.
    Refuted(String),
    /// The witness does not parse or was not answered.
    Unusable(String),
}

/// Runs a call on a clone of `template` and reports what the finished graph
/// says about it: the sorts of its arguments and of its answer.
///
/// `run` executes the engine on the graph for the given root (the caller owns
/// the engine and the budget).
///
/// # Errors
/// A description when the call does not parse, is not an application of
/// `op`, or was not answered.
pub fn observe_call(
    template: &Graph,
    op: &str,
    call: &str,
    run: &dyn Fn(&mut Graph, NodeId),
) -> Result<(Vec<Sort>, Sort), String> {
    let mut graph = template.clone();
    let root = graph
        .parse(call)
        .map_err(|error| format!("does not parse: {error}"))?;
    let root_op = graph.op(root);
    if &*graph.ops().get(root_op).name != op {
        return Err("the witness parses to another operator".to_owned());
    }
    let children = graph.children(root).to_vec();
    run(&mut graph, root);
    let mut memo = HashMap::new();
    let result = classify_without(&graph, graph.find(root), Some(root_op), &mut memo, 0)
        .ok_or_else(|| "the call was not answered".to_owned())?;
    let mut params = Vec::with_capacity(children.len());
    for (k, child) in children.iter().enumerate() {
        params.push(
            classify(&graph, graph.find(*child), &mut memo, 0)
                .ok_or_else(|| format!("argument {k} has no sort"))?,
        );
    }
    Ok((params, result))
}

/// Re-runs one witness on a clone of `template` and checks it against the
/// candidate it is meant to justify.
pub fn replay(
    template: &Graph,
    candidate: &Candidate,
    fact: &Fact,
    run: &dyn Fn(&mut Graph, NodeId),
) -> Replay {
    let (params, result) = match observe_call(template, &candidate.op, &fact.witness, run) {
        | Ok(seen) => seen,
        | Err(why) => return Replay::Unusable(why),
    };
    if !result.le(&candidate.result) {
        return Replay::Refuted(format!(
            "result {result} is not within {}",
            candidate.result
        ));
    }
    for (k, sort) in params.iter().enumerate() {
        let want = if candidate.variadic {
            candidate.params.first()
        } else {
            candidate.params.get(k)
        };
        if !want.is_some_and(|w| sort.le(w)) {
            return Replay::Refuted(format!("argument {k} has sort {sort}"));
        }
    }
    Replay::Confirmed(result)
}

/// Statistics of a mining run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MineStats {
    /// Facts read.
    pub facts: usize,
    /// Candidate overloads after clustering.
    pub candidates: usize,
    /// Signatures kept.
    pub kept: usize,
    /// Distinct calls of the trace that the full engine did not answer.
    pub unconfirmed: usize,
    /// Candidates dropped because they would make an overload set ambiguous.
    pub incoherent: usize,
    /// Kept signatures that refine the result of an existing overload.
    pub refinements: usize,
    /// Candidates whose result was widened by a probe (see [`probe`]).
    pub widened: usize,
    /// Operators with at least one witnessed signature.
    pub ops: usize,
}

/// Guards a signature against evidence from special values.
///
/// `sin(0) = 0` and `matsqrt` of a matrix of perfect squares are true facts
/// and terrible types: a witness at a special point says "integers in,
/// integers out" for an operator that is symbolic almost everywhere else.
/// The cure is to ask again with *other values of the same sorts*. The
/// numbers in the arguments of the candidate's best witnesses are shifted
/// (literals, and literals inside lists, so every shape is preserved) and the
/// call is run again:
///
/// * an answer is joined into the result sort;
/// * a call the engine leaves as it is contributes [`Sort::Expr`] when the
///   operator is not a request (an unevaluated function application is a
///   legitimate symbolic value), and nothing when it is one (a request that
///   is not answered failed, typically because the shifted numbers left its
///   domain, and the wrapper's error channel covers that).
///
/// Two further rules apply when the evidence is thin (fewer than
/// [`THIN_SUPPORT`] facts): exact number leaves of the result are widened to
/// expressions if some argument is symbolic, because one symbolic call that
/// happened to produce integers does not make an operator integer-valued.
///
/// Returns the result sort after probing.
pub fn probe(
    template: &Graph,
    candidate: &Candidate,
    run: &dyn Fn(&mut Graph, NodeId),
) -> Sort {
    let Some(op) = template.ops().lookup(&candidate.op) else {
        return candidate.result.clone();
    };
    let request = template.ops().get(op).flags.has(OpFlags::HEAVY);
    let mut result = candidate.result.clone();
    for fact in candidate.facts.iter().take(2) {
        for round in 1..=2_i64 {
            let mut graph = template.clone();
            let Ok(root) = graph.parse(&fact.witness) else {
                continue;
            };
            if graph.op(root) != op {
                continue;
            }
            // Shift the arguments as the engine sees them, not as the
            // witness spells them: `eigenvals(graph_laplacian(g))` has no
            // literal to shift until its argument is reduced to a matrix.
            run(&mut graph, root);
            let args: Vec<NodeId> = graph
                .children(root)
                .to_vec()
                .into_iter()
                .map(|arg| {
                    crate::graph::Extractor::new(&graph, &[arg], &crate::graph::SizeCost)
                        .build(&mut graph, arg)
                        .unwrap_or(arg)
                })
                .collect();
            let mut changed = false;
            let shifted: Vec<NodeId> = args
                .iter()
                .enumerate()
                .map(|(position, &arg)| shift(&mut graph, arg, round, position, &mut changed, 0))
                .collect();
            if !changed {
                continue;
            }
            let Some(call) = graph.try_node(op, &shifted) else {
                continue;
            };
            run(&mut graph, call);
            let mut memo = HashMap::new();
            match classify_without(&graph, graph.find(call), Some(op), &mut memo, 0) {
                | Some(seen) => result = result.join(&seen),
                | None if !request => result = result.join(&Sort::Expr),
                | None => {},
            }
        }
    }
    if candidate.facts.len() < THIN_SUPPORT && candidate.params.iter().any(has_symbolic_leaf) {
        result = widen_exact_leaves(&result);
    }
    result
}

/// Facts below which a cluster counts as thin evidence (see [`probe`]).
pub const THIN_SUPPORT: usize = 3;

/// `node` with its number literals shifted (also inside lists).
fn shift(
    graph: &mut Graph,
    node: NodeId,
    round: i64,
    position: usize,
    changed: &mut bool,
    depth: usize,
) -> NodeId {
    let step = round + i64::try_from(position % 3).unwrap_or(0);
    if let Some(number) = graph.as_number(node).cloned() {
        *changed = true;
        return match number {
            | Number::Int(n) => graph.num(Number::Int(n + num_bigint::BigInt::from(step))),
            | Number::Rat(q) => {
                let delta = num_rational::BigRational::new(
                    num_bigint::BigInt::from(1),
                    num_bigint::BigInt::from(step + 1),
                );
                graph.num(Number::rat(q + delta))
            },
            #[allow(clippy::cast_precision_loss)]
            | Number::Float(x) => graph.float(x + 0.37 * step as f64),
        };
    }
    if graph.op(node) == core::LIST && depth < 4 {
        let children = graph.children(node).to_vec();
        let shifted: Vec<NodeId> = children
            .iter()
            .enumerate()
            .map(|(index, &child)| shift(graph, child, round, index, changed, depth + 1))
            .collect();
        return graph.node(core::LIST, &shifted);
    }
    node
}

fn has_symbolic_leaf(sort: &Sort) -> bool {
    match sort {
        | Sort::Expr | Sort::Symbol => true,
        | Sort::List { elem, .. } => has_symbolic_leaf(elem),
        | Sort::Tuple(items) => items.iter().any(has_symbolic_leaf),
        | _ => false,
    }
}

fn widen_exact_leaves(sort: &Sort) -> Sort {
    match sort {
        | Sort::Int | Sort::Rational => Sort::Expr,
        | Sort::List { elem, rect } => {
            Sort::List {
                elem: Box::new(widen_exact_leaves(elem)),
                rect: *rect,
            }
        },
        | Sort::Tuple(items) => Sort::Tuple(items.iter().map(widen_exact_leaves).collect()),
        | other => other.clone(),
    }
}

/// Turns facts into a table of replayed signatures.
///
/// A trace is recorded under whatever rule sets each test loaded, so a fact
/// may describe a call that a partial engine left half-answered. The trace is
/// therefore only a source of *calls*: every distinct witness is run again on
/// `template` (the full engine) and the sorts used from here on are the ones
/// that run shows ([`observe_call`]). Calls it does not answer are dropped.
/// The verified facts are then clustered ([`candidates`]) and a signature is
/// kept when it is informative and coherent.
///
/// `existing` returns the signatures an operator already has from its
/// registration (declared, overlay, derived): a witnessed signature must be
/// coherent with them too, and one with a parameter list the registration
/// already states is not repeated. `binder` is as in [`candidates`]; `run` as
/// in [`observe_call`].
pub fn mine(
    template: &Graph,
    facts: &[Fact],
    existing: &dyn Fn(&str, usize, bool) -> Vec<Signature>,
    binder: &dyn Fn(&str) -> Option<usize>,
    run: &dyn Fn(&mut Graph, NodeId),
) -> (WitnessedTable, MineStats) {
    let mut stats = MineStats {
        facts: facts.len(),
        ..MineStats::default()
    };
    let mut seen: std::collections::BTreeSet<(&str, &str)> = std::collections::BTreeSet::new();
    let mut verified: Vec<Fact> = Vec::new();
    for fact in facts {
        if !seen.insert((fact.op.as_str(), fact.witness.as_str())) {
            continue;
        }
        match observe_call(template, &fact.op, &fact.witness, run) {
            | Ok((params, result)) => {
                verified.push(Fact {
                    op: fact.op.clone(),
                    variadic: fact.variadic,
                    params,
                    result,
                    witness: fact.witness.clone(),
                });
            },
            | Err(_) => stats.unconfirmed += 1,
        }
    }
    let all = candidates(&verified, binder);
    stats.candidates = all.len();
    let mut table = WitnessedTable {
        version: 1,
        ops: BTreeMap::new(),
    };
    let mut kept: BTreeMap<(String, usize, bool), Vec<Signature>> = BTreeMap::new();
    for candidate in &all {
        let mut candidate = candidate.clone();
        let probed = probe(template, &candidate, run);
        if probed != candidate.result {
            stats.widened += 1;
            candidate.result = probed;
            if !candidate.informative() {
                continue;
            }
        }
        let candidate = &candidate;
        let signature = candidate.signature();
        let arity = if candidate.variadic {
            0
        } else {
            candidate.params.len()
        };
        let registered = existing(&candidate.op, arity, candidate.variadic);
        let key = (candidate.op.clone(), arity, candidate.variadic);
        let slot = kept.entry(key).or_insert_with(|| registered.clone());
        let same_shape = |old: &Signature| {
            old.params == signature.params
                && old.variadic.as_ref().map(|v| &v.elem)
                    == signature.variadic.as_ref().map(|v| &v.elem)
        };
        // The same parameters as a registered overload: not a new overload,
        // but possibly a sharper result for that one.
        let refines = registered.iter().any(&same_shape);
        if refines && (!candidate.result.informative() || registered.contains(&signature)) {
            stats.incoherent += 1;
            continue;
        }
        let others: Vec<Signature> = slot
            .iter()
            .filter(|old| !(refines && registered.iter().any(|r| r == *old) && same_shape(old)))
            .cloned()
            .collect();
        let already = table.ops.get(&candidate.op).map_or(0, Vec::len);
        if already >= MAX_OVERLOADS || !coherent(&others, &signature) {
            stats.incoherent += 1;
            continue;
        }
        let witnesses: Vec<Witness> = candidate
            .facts
            .iter()
            .take(MAX_WITNESSES)
            .map(|fact| {
                Witness {
                    input: fact.witness.clone(),
                    output: fact.result.to_string(),
                }
            })
            .collect();
        slot.push(signature.clone());
        if refines {
            stats.refinements += 1;
        }
        table
            .ops
            .entry(candidate.op.clone())
            .or_default()
            .push(Witnessed {
                signature,
                params: candidate.params.iter().map(ToString::to_string).collect(),
                result: candidate.result.to_string(),
                support: u32::try_from(candidate.facts.len()).unwrap_or(u32::MAX),
                witnesses,
                refines,
            });
        stats.kept += 1;
    }
    stats.ops = table.ops.len();
    (table, stats)
}

/// The witnessed signatures shipped with the crate
/// (`signatures.witnessed.json`, materialised by `build.rs`).
///
/// # Panics
///
/// Panics if the embedded table is not valid JSON (it is generated and
/// version-controlled, so this indicates a broken build).
#[must_use]
pub fn shipped() -> &'static WitnessedTable {
    static TABLE: OnceLock<WitnessedTable> = OnceLock::new();
    TABLE.get_or_init(|| {
        serde_json::from_str::<WitnessedTable>(WITNESSED_JSON)
            .expect("version-controlled signatures.witnessed.json must be a valid table")
    })
}

const WITNESSED_JSON: &str = include_str!(concat!(env!("OUT_DIR"), "/signatures.witnessed.json"));

/// The shipped result refinements of `desc`: signatures whose parameter list
/// an overload of `desc` already has, with a sharper result (see
/// [`Witnessed::refines`]).
#[must_use]
pub fn refinements_for(desc: &OpDescriptor) -> Vec<&'static Witnessed> {
    witnessed_for(desc)
        .into_iter()
        .filter(|w| w.refines)
        .collect()
}

/// The shipped witnessed signatures of `desc`, in table order (overloads
/// and refinements).
#[must_use]
pub fn witnessed_for(desc: &OpDescriptor) -> Vec<&'static Witnessed> {
    let Some(entries) = shipped().ops.get(&*desc.name) else {
        return Vec::new();
    };
    entries
        .iter()
        .filter(|w| {
            match desc.arity {
                | Arity::Fixed(n) => {
                    w.signature.variadic.is_none() && w.signature.params.len() == n as usize
                },
                | Arity::Variadic => w.signature.variadic.is_some(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Saturate;
    use crate::rules;

    fn solved(
        sets: &[crate::graph::RuleSet],
        src: &str,
    ) -> (Graph, NodeId) {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, sets).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        (g, root)
    }

    #[test]
    fn the_lattice_joins_as_documented() {
        assert_eq!(Sort::Int.join(&Sort::Rational), Sort::Rational);
        assert_eq!(Sort::Int.join(&Sort::Float), Sort::Expr);
        assert_eq!(Sort::Symbol.join(&Sort::Int), Sort::Expr);
        assert_eq!(Sort::Bool.join(&Sort::Prop), Sort::Prop);
        assert_eq!(Sort::Bool.join(&Sort::Int), Sort::Any);
        assert_eq!(
            Sort::list(Sort::Int).join(&Sort::list(Sort::Never)),
            Sort::list(Sort::Int)
        );
        assert_eq!(
            Sort::Tuple(vec![Sort::Int, Sort::Bool])
                .join(&Sort::Tuple(vec![Sort::Rational, Sort::Bool])),
            Sort::Tuple(vec![Sort::Rational, Sort::Bool])
        );
        assert!(Sort::Int.le(&Sort::Expr));
        assert!(!Sort::Expr.le(&Sort::Int));
        assert!(Sort::list(Sort::Int).le(&Sort::list(Sort::Expr)));
    }

    #[test]
    fn sorts_round_trip_through_text() {
        for sort in [
            Sort::Int,
            Sort::list(Sort::list(Sort::Rational)),
            Sort::List {
                elem: Box::new(Sort::list(Sort::Int)),
                rect: false,
            },
            Sort::Tuple(vec![Sort::Int, Sort::list(Sort::Expr), Sort::Bool]),
            Sort::Object("interval".to_owned()),
            Sort::Tuple(Vec::new()),
        ] {
            assert_eq!(Sort::parse(&sort.to_string()), Some(sort));
        }
        assert_eq!(Sort::parse("[Int"), None);
        assert_eq!(Sort::parse("Int junk"), None);
    }

    #[test]
    fn arrays_become_tensors_and_tuples_records() {
        let matrix = Sort::list(Sort::list(Sort::Int));
        assert_eq!(matrix.array(), Some((2, &Sort::Int)));
        assert_eq!(
            matrix.type_expr(),
            TypeExpr::Tensor {
                rank: NatExpr::Lit(2),
                shape: None,
                elem: Box::new(TypeExpr::Int),
            }
        );
        let ragged = Sort::List {
            elem: Box::new(Sort::list(Sort::Int)),
            rect: false,
        };
        assert_eq!(ragged.array(), None);
        assert!(matches!(
            Sort::Tuple(vec![Sort::Int, Sort::Bool]).type_expr(),
            TypeExpr::Record(fields) if fields.len() == 2 && fields[0].0 == "_0"
        ));
    }

    #[test]
    fn a_finished_graph_states_the_type_of_a_request() {
        let (g, root) = solved(&[rules::number_theory::number_theory()], "factorint(360)");
        assert_eq!(
            sort_of_node(&g, root),
            Some(Sort::list(Sort::list(Sort::Int)))
        );
        let facts = facts_of(&g);
        let fact = facts
            .iter()
            .find(|f| f.op == "factorint")
            .expect("a fact about factorint");
        assert_eq!(fact.params, vec![Sort::Int]);
        assert_eq!(Fact::parse(&fact.line()).as_ref(), Some(fact));
    }

    #[test]
    fn an_unanswered_request_has_no_sort() {
        let (g, root) = solved(&[rules::calculus::calculus()], "integral(f(x)^g(x), x)");
        assert_eq!(sort_of_node(&g, root), None);
    }

    #[test]
    fn an_application_is_not_evidence_about_itself() {
        // `sin(x)` stays as it is: alone in its class, it says nothing about
        // what `sin` returns. `sin(0)` shares its class with `0`.
        let (g, _) = solved(&[rules::elementary::elementary()], "list(sin(x), sin(0))");
        let facts = facts_of(&g);
        assert!(facts.iter().all(|f| f.witness != "sin(x)"), "{facts:?}");
        assert!(
            facts
                .iter()
                .any(|f| f.op == "sin" && f.params == vec![Sort::Int] && f.result == Sort::Int),
            "sin(0) = 0 is evidence: {facts:?}"
        );
    }

    #[test]
    fn ragged_lists_are_not_arrays() {
        let (g, root) = solved(
            &[rules::arith::arith()],
            "list(list(1, 2), list(3), list())",
        );
        let sort = sort_of_node(&g, root).expect("a list has a sort");
        assert_eq!(sort.array(), None);
        assert!(matches!(sort, Sort::List { rect: false, .. }));
    }

    #[test]
    fn coherent_sets_are_chains_or_disjoint() {
        let term = crate::kernel::infer::term;
        let sig = |params: Vec<TypeExpr>| Signature::new(params, term());
        let generic = sig(vec![term(), term()]);
        let ints = sig(vec![TypeExpr::Int, TypeExpr::Int]);
        let left = sig(vec![TypeExpr::Int, term()]);
        let right = sig(vec![term(), TypeExpr::Int]);
        assert!(coherent(std::slice::from_ref(&generic), &ints));
        assert!(coherent(&[generic.clone(), ints.clone()], &left));
        // `(Int, Term)` and `(Term, Int)` both accept `(Int, Int)` and neither
        // is more specific.
        assert!(!coherent(&[generic, left], &right));
        // The same parameter list twice is never coherent.
        assert!(!coherent(std::slice::from_ref(&ints), &ints));
        // A matrix parameter and a scalar one cannot receive the same call.
        let matrix = sig(vec![TypeExpr::Int.tensor(2_u64), TypeExpr::Int]);
        assert!(coherent(&[ints], &matrix));
    }

    /// The shipped table never makes a call ambiguous: for every operator,
    /// any two overloads are comparable or disjoint, and each witness of a
    /// sample of operators still replays.
    #[test]
    fn the_shipped_table_is_coherent_and_replays() {
        let sets = rules::standard();
        let mut template = Graph::new();
        let engine = Engine::install(&mut template, &sets).unwrap_or_else(|e| panic!("{e}"));
        let mut witnessed_ops = 0;
        for (_, desc) in template.ops().iter() {
            let (signatures, source) = crate::kernel::infer::signatures(desc);
            if source != crate::kernel::manifest::SignatureSource::Witnessed {
                continue;
            }
            witnessed_ops += 1;
            for (i, a) in signatures.iter().enumerate() {
                assert!(
                    coherent(&signatures[..i], a),
                    "`{}`: overload {i} overlaps an earlier one",
                    desc.name
                );
            }
            let sources = crate::kernel::infer::overload_sources(desc);
            assert_eq!(sources.len(), signatures.len());
        }
        assert_eq!(
            witnessed_ops == 0,
            shipped().ops.is_empty(),
            "a non-empty table must reach the operator table"
        );
        let budget = Budget {
            time_limit: Some(std::time::Duration::from_secs(20)),
            ..Budget::default()
        };
        let run = |g: &mut Graph, root: NodeId| {
            engine.run(g, &[root], &Env::symbolic(), &Saturate, &budget);
        };
        // Replaying everything is the miner's `--check`; a fixed sample keeps
        // the unit test quick while still failing when the table goes stale.
        for (name, entries) in shipped().ops.iter().step_by(23) {
            for entry in entries {
                assert!(!entry.witnesses.is_empty(), "{name}");
                let params: Vec<Sort> = entry
                    .params
                    .iter()
                    .map(|p| Sort::parse(p).unwrap_or_else(|| panic!("{name}: sort `{p}`")))
                    .collect();
                let candidate = Candidate {
                    op: name.clone(),
                    variadic: entry.signature.variadic.is_some(),
                    params: params.clone(),
                    result: Sort::parse(&entry.result)
                        .unwrap_or_else(|| panic!("{name}: sort `{}`", entry.result)),
                    facts: Vec::new(),
                };
                assert_eq!(candidate.signature(), entry.signature, "{name}");
                let witness = entry.witnesses.first().unwrap_or_else(|| panic!("{name}"));
                let fact = Fact {
                    op: name.clone(),
                    variadic: candidate.variadic,
                    params,
                    result: candidate.result.clone(),
                    witness: witness.input.clone(),
                };
                assert!(
                    matches!(
                        replay(&template, &candidate, &fact, &run),
                        Replay::Confirmed(_)
                    ),
                    "`{}` no longer replays; run signature-mine again",
                    witness.input
                );
            }
        }
    }

    #[test]
    fn special_values_do_not_become_types_and_results_are_refined() {
        let sets = [
            rules::elementary::elementary(),
            rules::number_theory::number_theory(),
        ];
        let mut template = Graph::new();
        let engine = Engine::install(&mut template, &sets).unwrap_or_else(|e| panic!("{e}"));
        let run = |g: &mut Graph, root: NodeId| {
            engine.run(g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        };
        let fact = |op: &str, witness: &str| {
            Fact {
                op: op.to_owned(),
                variadic: false,
                params: Vec::new(),
                result: Sort::Any,
                witness: witness.to_owned(),
            }
        };
        let term = crate::kernel::infer::term;
        let (table, stats) = mine(
            &template,
            &[
                fact("sin", "sin(0)"),
                fact("isprime", "isprime(7)"),
                fact("gcd", "gcd(12, 18)"),
            ],
            // Pretend `isprime` is registered as `(Int) -> Term`.
            &|name, n, _| {
                if name == "isprime" {
                    vec![Signature::new(vec![TypeExpr::Int], term())]
                } else {
                    vec![Signature::new(vec![term(); n], term())]
                }
            },
            &|_| None,
            &run,
        );
        // `sin(0) = 0` does not make `sin` integer-valued: `sin(2)` stays.
        assert!(!table.ops.contains_key("sin"), "{table:?}");
        assert!(stats.widened >= 1);
        // Shapes are kept while numbers move: `divisors` stays a list of
        // integers at other arguments.
        let (lists, _) = mine(
            &template,
            &[fact("divisors", "divisors(12)")],
            &|_, n, _| vec![Signature::new(vec![term(); n], term())],
            &|_| None,
            &run,
        );
        assert_eq!(lists.ops["divisors"][0].result, "[Int]");
        // `gcd` is integer-valued at the probes as well.
        assert_eq!(table.ops["gcd"][0].signature.ret, TypeExpr::Int);
        assert!(!table.ops["gcd"][0].refines);
        // The same parameters as a registered overload, a sharper result.
        let prime = &table.ops["isprime"][0];
        assert!(prime.refines && prime.signature.ret == TypeExpr::Bool);
        assert_eq!(stats.refinements, 1);
    }

    #[test]
    fn mining_keeps_only_what_replays() {
        let sets = [rules::number_theory::number_theory()];
        let mut template = Graph::new();
        let engine = Engine::install(&mut template, &sets).unwrap_or_else(|e| panic!("{e}"));
        let run = |g: &mut Graph, root: NodeId| {
            engine.run(g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        };
        let good = Fact {
            op: "divisors".to_owned(),
            variadic: false,
            params: vec![Sort::Int],
            result: Sort::list(Sort::Int),
            witness: "divisors(12)".to_owned(),
        };
        // A recorded sort is only a hint: the replay decides, so a wrong
        // one is corrected rather than believed.
        let wrong = Fact {
            op: "totient".to_owned(),
            variadic: false,
            params: vec![Sort::Int],
            result: Sort::Bool,
            witness: "totient(10)".to_owned(),
        };
        // A call the engine does not answer is dropped.
        let bad = Fact {
            op: "totient".to_owned(),
            variadic: false,
            params: vec![Sort::Symbol],
            result: Sort::Int,
            witness: "totient(n)".to_owned(),
        };
        let (table, stats) = mine(
            &template,
            &[good, wrong, bad],
            &|_, n, _| {
                vec![Signature::new(
                    vec![crate::kernel::infer::term(); n],
                    crate::kernel::infer::term(),
                )]
            },
            &|_| None,
            &run,
        );
        assert_eq!(stats.kept, 2);
        assert_eq!(stats.unconfirmed, 1);
        let entry = &table.ops["divisors"][0];
        assert_eq!(entry.signature.params, vec![TypeExpr::Int]);
        assert_eq!(entry.witnesses[0].input, "divisors(12)");
        let totient = &table.ops["totient"][0];
        assert_eq!(totient.signature.ret, TypeExpr::Int);
        assert_eq!(totient.witnesses.len(), 1);
    }
}
