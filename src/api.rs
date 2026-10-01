//! # The user-facing API
//!
//! One entry point, [`Session::compute`], for everything that is an identity
//! transformation. What differs between "simplify", "differentiate",
//! "evaluate" and "solve" is only the term that is handed in and the
//! [`Target`] phase of the answer.
//!
//! ```
//! use rssn::api::*;
//!
//! let s = Session::new();
//! let x = s.sym("x");
//! let request = diff(sin(x) * cos(x), x);
//!
//! // Symbolic phase: a closed form.
//! let symbolic = s.compute(request, &Config::new())?;
//! assert!(symbolic.reduced);
//!
//! // Numeric phase: the same request, evaluated at x = 0.
//! let numeric = s.compute(request, &Config::new().numeric(1e-9).bind("x", 0.0))?;
//! assert_eq!(numeric.value, Some(1.0));
//! # Ok::<(), ComputeError>(())
//! ```

use std::cell::RefCell;
use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use std::sync::Arc;

use crate::graph::Budget;
use crate::graph::ClosedForm;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Evaluated;
use crate::graph::Extractor;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::ParseError;
use crate::graph::Report;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Saturate;
use crate::graph::SizeCost;
use crate::graph::op::core;
use crate::rules;

/// The representation the answer should be in.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Target {
    /// A closed-form term: no unevaluated heavy operator left.
    Symbolic,
    /// A floating point value within the given absolute tolerance.
    Numeric {
        /// Requested absolute tolerance.
        tolerance: f64,
    },
}

/// What to compute with and what to aim for.
#[derive(Clone, Debug)]
pub struct Config {
    sets: Vec<RuleSet>,
    target: Target,
    bindings: Vec<(String, f64)>,
    assumptions: Vec<(String, Facts)>,
    budget: Budget,
}

impl Default for Config {
    fn default() -> Self {
        Self::new()
    }
}

impl Config {
    /// All standard rule sets, symbolic target, default budget.
    #[must_use]
    pub fn new() -> Self {
        Self {
            sets: rules::standard(),
            target: Target::Symbolic,
            bindings: Vec::new(),
            assumptions: Vec::new(),
            budget: Budget::default(),
        }
    }

    /// No rule sets at all; add some with [`Config::with`].
    #[must_use]
    pub fn empty() -> Self {
        Self {
            sets: Vec::new(),
            target: Target::Symbolic,
            bindings: Vec::new(),
            assumptions: Vec::new(),
            budget: Budget::default(),
        }
    }

    /// Adds a rule set (and, implicitly, its dependencies).
    #[must_use]
    pub fn with(
        mut self,
        set: RuleSet,
    ) -> Self {
        self.sets.push(set);
        self
    }

    /// Asks for a closed-form answer.
    #[must_use]
    pub const fn symbolic(mut self) -> Self {
        self.target = Target::Symbolic;
        self
    }

    /// Asks for a numeric answer within `tolerance`.
    #[must_use]
    pub const fn numeric(
        mut self,
        tolerance: f64,
    ) -> Self {
        self.target = Target::Numeric { tolerance };
        self
    }

    /// Binds a symbol to a value for numeric evaluation.
    #[must_use]
    pub fn bind(
        mut self,
        symbol: &str,
        value: f64,
    ) -> Self {
        self.bindings.push((symbol.to_owned(), value));
        self
    }

    /// Declares facts about a symbol — for example
    /// [`Facts::POSITIVE`] or [`Facts::REAL`] — that conditional identities
    /// may rely on.
    #[must_use]
    pub fn assume(
        mut self,
        symbol: &str,
        facts: Facts,
    ) -> Self {
        self.assumptions.push((symbol.to_owned(), facts));
        self
    }

    /// Replaces the search budget.
    #[must_use]
    pub fn budget(
        mut self,
        budget: Budget,
    ) -> Self {
        self.budget = budget;
        self
    }
}

/// Why a computation could not produce an answer.
#[derive(Clone, Debug, PartialEq)]
pub enum ComputeError {
    /// A rule set failed to install.
    Rule(RuleError),
    /// Input text could not be parsed.
    Parse(ParseError),
    /// An operator name is not registered.
    UnknownOperator(String),
    /// A numeric answer was requested but the term could not be evaluated;
    /// carries the best symbolic form that was reached.
    NotNumeric(String),
    /// The rules derived that two different numbers are equal. This is a
    /// bug in a rule set; the session must be discarded.
    Inconsistent(String),
}

impl fmt::Display for ComputeError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Rule(e) => write!(f, "{e}"),
            | Self::Parse(e) => write!(f, "{e}"),
            | Self::UnknownOperator(name) => write!(f, "unknown operator `{name}`"),
            | Self::NotNumeric(term) => {
                write!(f, "no numeric value could be determined for `{term}`")
            },
            | Self::Inconsistent(what) => write!(f, "rule set is unsound: derived {what}"),
        }
    }
}

impl std::error::Error for ComputeError {}

impl From<RuleError> for ComputeError {
    fn from(value: RuleError) -> Self {
        Self::Rule(value)
    }
}

impl From<ParseError> for ComputeError {
    fn from(value: ParseError) -> Self {
        Self::Parse(value)
    }
}

/// A graph with a set of rule sets installed, and the engine bound to it.
struct Template {
    sets: Vec<Arc<str>>,
    graph: Graph,
    engine: Engine,
}

struct Inner {
    /// The term store: every [`Term`] is a node of this graph. It is a
    /// plain hash-consed DAG; no equalities are ever asserted in it.
    graph: Graph,
    templates: Vec<Template>,
}

/// Owns the term store and the engines that compute on it.
///
/// Each [`Session::compute`] call runs in a private copy of a template
/// graph: the request is copied in, the search happens there, and only the
/// answer is copied back. Requests therefore cannot slow each other down or
/// influence each other's results, and a session can live as long as the
/// program.
///
/// A session is single-threaded; create one per thread. Terms built in one
/// session must not be used in another.
pub struct Session {
    inner: RefCell<Inner>,
}

impl fmt::Debug for Session {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let inner = self.inner.borrow();
        write!(
            f,
            "Session({} terms, {} rule configurations)",
            inner.graph.len(),
            inner.templates.len()
        )
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

/// A handle to a term of a [`Session`]. Cheap to copy.
#[derive(Copy, Clone)]
pub struct Term<'s> {
    session: &'s Session,
    node: NodeId,
}

/// The result of [`Session::compute`].
#[derive(Clone, Debug)]
pub struct Answer<'s> {
    /// The answer as a term: a closed form for a symbolic target, a float
    /// literal for a numeric one.
    pub term: Term<'s>,
    /// The numeric value, when one is known.
    pub value: Option<f64>,
    /// Estimated absolute error of `value`.
    pub error: Option<f64>,
    /// Whether the request was fully reduced. `false` means `term` still
    /// contains unevaluated heavy operators.
    pub reduced: bool,
    /// What the engine did.
    pub report: Report,
}

impl Session {
    /// A session with the operators of every standard rule set registered.
    ///
    /// # Panics
    /// Panics if the built-in rule sets fail to install, which would be a
    /// bug in rssn itself.
    #[must_use]
    pub fn new() -> Self {
        match Self::with_rules(&rules::standard()) {
            | Ok(session) => session,
            | Err(e) => broken_builtin(&e),
        }
    }

    /// A session with exactly the given rule sets pre-installed.
    ///
    /// # Errors
    /// Fails when a rule set does not install.
    pub fn with_rules(sets: &[RuleSet]) -> Result<Self, RuleError> {
        let template = Template::new(sets)?;
        let graph = template.graph.clone();
        Ok(Self {
            inner: RefCell::new(Inner {
                graph,
                templates: vec![template],
            }),
        })
    }

    const fn term(
        &self,
        node: NodeId,
    ) -> Term<'_> {
        Term { session: self, node }
    }

    /// The symbol named `name`.
    pub fn sym(
        &self,
        name: &str,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.sym(name))
    }

    /// An exact integer.
    pub fn int(
        &self,
        value: i64,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.int(value))
    }

    /// The exact fraction `numerator / denominator`; `None` for a zero
    /// denominator.
    pub fn rational(
        &self,
        numerator: i64,
        denominator: i64,
    ) -> Option<Term<'_>> {
        let n = Number::fraction(numerator, denominator)?;
        Some(self.term(self.inner.borrow_mut().graph.num(n)))
    }

    /// A floating point number.
    pub fn float(
        &self,
        value: f64,
    ) -> Term<'_> {
        self.term(self.inner.borrow_mut().graph.float(value))
    }

    /// Parses a term from infix text.
    ///
    /// # Errors
    /// Returns the parser's complaint.
    pub fn parse(
        &self,
        text: &str,
    ) -> Result<Term<'_>, ComputeError> {
        Ok(self.term(self.inner.borrow_mut().graph.parse(text)?))
    }

    /// Applies the operator named `op` to `args`.
    ///
    /// # Errors
    /// Fails when the operator is unknown or does not take that many
    /// arguments.
    pub fn call<'s>(
        &'s self,
        op: &str,
        args: &[Term<'s>],
    ) -> Result<Term<'s>, ComputeError> {
        let mut inner = self.inner.borrow_mut();
        let unknown = || ComputeError::UnknownOperator(op.to_owned());
        let id = inner.graph.ops().lookup(op).ok_or_else(unknown)?;
        let nodes: Vec<NodeId> = args.iter().map(|t| t.node).collect();
        inner
            .graph
            .try_node(id, &nodes)
            .map(|n| self.term(n))
            .ok_or_else(unknown)
    }

    /// Reduces `term` to the target phase of `config`.
    ///
    /// # Errors
    /// * [`ComputeError::Rule`] if a rule set of `config` fails to install;
    /// * [`ComputeError::NotNumeric`] if a numeric target cannot be met;
    /// * [`ComputeError::Inconsistent`] if an unsound rule was detected.
    pub fn compute<'s>(
        &'s self,
        term: Term<'s>,
        config: &Config,
    ) -> Result<Answer<'s>, ComputeError> {
        let mut guard = self.inner.borrow_mut();
        let inner = &mut *guard;
        let key: Vec<Arc<str>> = config.sets.iter().map(|s| Arc::from(s.name())).collect();
        let index = match inner.templates.iter().position(|t| t.sets == key) {
            | Some(i) => i,
            | None => {
                inner.templates.push(Template::new(&config.sets)?);
                inner.templates.len().saturating_sub(1)
            },
        };
        let Some(template) = inner.templates.get(index) else {
            return Err(ComputeError::UnknownOperator("engine".to_owned()));
        };
        let store = &mut inner.graph;

        // Work in a private copy; only the answer goes back to the store.
        let mut graph = template.graph.clone();
        let root = graph.import(store, term.node);
        let mut env = match config.target {
            | Target::Symbolic => Env::symbolic(),
            | Target::Numeric { tolerance } => Env::numeric(tolerance),
        };
        for (name, value) in &config.bindings {
            env.bind(graph.interner_mut().symbol(name), *value);
        }
        for (name, facts) in &config.assumptions {
            let symbol = graph.interner_mut().symbol(name);
            graph.assume(symbol, *facts);
        }
        let roots = [root];
        let report = match config.target {
            | Target::Symbolic => {
                template
                    .engine
                    .run(&mut graph, &roots, &env, &Saturate, &config.budget)
            },
            | Target::Numeric { .. } => {
                template
                    .engine
                    .run(&mut graph, &roots, &env, &Evaluated, &config.budget)
            },
        };
        if let Some((a, b)) = graph.conflicts().first() {
            return Err(ComputeError::Inconsistent(format!("{a} = {b}")));
        }

        let witness = graph.approx(graph.find(root));
        let closed = Extractor::new(&graph, &roots, &ClosedForm).build(&mut graph, root);
        let reduced = closed.is_some();
        let best = match closed {
            | Some(node) => node,
            | None => Extractor::new(&graph, &roots, &SizeCost)
                .build(&mut graph, root)
                .unwrap_or(root),
        };
        match config.target {
            | Target::Symbolic => {
                let exact = graph.as_number(best).map(Number::to_f64);
                let node = store.import(&graph, best);
                Ok(Answer {
                    term: self.term(node),
                    value: exact,
                    error: exact.map(|_| 0.0),
                    reduced,
                    report,
                })
            },
            | Target::Numeric { .. } => {
                // Prefer evaluating the closed form; fall back to the
                // witness numeric kernels produced.
                let ball = graph
                    .eval(best, &env)
                    .map(crate::graph::Ball::exact)
                    .or(witness);
                let Some(ball) = ball else {
                    return Err(ComputeError::NotNumeric(graph.display(best)));
                };
                let node = store.float(ball.mid);
                Ok(Answer {
                    term: self.term(node),
                    value: Some(ball.mid),
                    error: Some(ball.rad),
                    reduced,
                    report,
                })
            },
        }
    }

    /// Shorthand for a symbolic [`Session::compute`] with the standard
    /// rules, returning only the term.
    ///
    /// # Errors
    /// See [`Session::compute`].
    pub fn simplify<'s>(
        &'s self,
        term: Term<'s>,
    ) -> Result<Term<'s>, ComputeError> {
        Ok(self.compute(term, &Config::new())?.term)
    }
}

impl Template {
    fn new(sets: &[RuleSet]) -> Result<Self, RuleError> {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, sets)?;
        Ok(Self {
            sets: sets.iter().map(|s| Arc::from(s.name())).collect(),
            graph,
            engine,
        })
    }
}

#[cold]
fn broken_builtin(error: &RuleError) -> ! {
    panic!("rssn: built-in rule sets failed to install: {error}")
}

impl<'s> Term<'s> {
    /// The session this term belongs to.
    #[must_use]
    pub const fn session(self) -> &'s Session {
        self.session
    }

    /// The literal number this term is, if it is one.
    #[must_use]
    pub fn as_number(self) -> Option<Number> {
        self.session
            .inner
            .borrow()
            .graph
            .as_number(self.node)
            .cloned()
    }

    /// The term as a float if it is a literal number.
    #[must_use]
    pub fn as_f64(self) -> Option<f64> {
        self.as_number().map(|n| n.to_f64())
    }

    /// Evaluates the term numerically under the given bindings, without
    /// running any rules.
    #[must_use]
    pub fn eval(
        self,
        bindings: &[(&str, f64)],
    ) -> Option<f64> {
        let mut inner = self.session.inner.borrow_mut();
        let mut env = Env::numeric(0.0);
        for (name, value) in bindings {
            env.bind(inner.graph.interner_mut().symbol(name), *value);
        }
        inner.graph.eval(self.node, &env)
    }

    /// Compiles the term to a function of the named symbols, in order,
    /// using the reference [`Interpreter`](crate::backend::Interpreter).
    ///
    /// # Errors
    /// See [`BackendError`](crate::backend::BackendError).
    pub fn compile(
        self,
        inputs: &[&str],
    ) -> Result<Box<dyn crate::backend::Compiled>, crate::backend::BackendError> {
        self.compile_with(&crate::backend::Interpreter, inputs)
    }

    /// Compiles the term with a specific backend.
    ///
    /// # Errors
    /// See [`BackendError`](crate::backend::BackendError).
    pub fn compile_with(
        self,
        backend: &dyn crate::backend::Backend,
        inputs: &[&str],
    ) -> Result<Box<dyn crate::backend::Compiled>, crate::backend::BackendError> {
        let mut inner = self.session.inner.borrow_mut();
        let symbols: Vec<_> = inputs
            .iter()
            .map(|name| inner.graph.interner_mut().symbol(name))
            .collect();
        backend.compile(&inner.graph, self.node, &symbols)
    }

    /// Applies a unary operator that the standard rule sets register.
    fn apply(
        self,
        op: &str,
        rest: &[Self],
    ) -> Self {
        let mut args = Vec::with_capacity(rest.len().saturating_add(1));
        args.push(self);
        args.extend_from_slice(rest);
        match self.session.call(op, &args) {
            | Ok(term) => term,
            | Err(e) => missing_operator(&e),
        }
    }

    fn binary(
        self,
        op: crate::graph::OpId,
        other: Self,
    ) -> Self {
        let node = self
            .session
            .inner
            .borrow_mut()
            .graph
            .node(op, &[self.node, other.node]);
        self.session.term(node)
    }

    /// `self ^ exponent`.
    #[must_use]
    pub fn pow(
        self,
        exponent: impl IntoTerm<'s>,
    ) -> Self {
        self.binary(core::POW, exponent.into_term(self.session))
    }
}

#[cold]
fn missing_operator(error: &ComputeError) -> ! {
    panic!("rssn: {error}; the session was created without the rule set that defines it")
}

impl PartialEq for Term<'_> {
    /// Structural identity: the same term of the same session.
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        std::ptr::eq(self.session, other.session) && self.node == other.node
    }
}

impl Eq for Term<'_> {}

impl fmt::Display for Term<'_> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str(&self.session.inner.borrow().graph.display(self.node))
    }
}

impl fmt::Debug for Term<'_> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "Term({self})")
    }
}

/// Values that can stand where a [`Term`] is expected: terms themselves,
/// `i32` (exact integers) and `f64`.
///
/// Only one integer type is accepted on purpose, so that a bare literal in
/// `2 * x + 1` infers without annotations. Larger integers go through
/// [`Session::int`].
pub trait IntoTerm<'s> {
    /// Converts into a term of `session`.
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s>;
}

impl<'s> IntoTerm<'s> for Term<'s> {
    fn into_term(
        self,
        _session: &'s Session,
    ) -> Term<'s> {
        self
    }
}

impl<'s> IntoTerm<'s> for i32 {
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s> {
        session.int(i64::from(self))
    }
}

impl<'s> IntoTerm<'s> for f64 {
    fn into_term(
        self,
        session: &'s Session,
    ) -> Term<'s> {
        session.float(self)
    }
}

macro_rules! arithmetic {
    ($trait:ident, $method:ident, |$a:ident, $b:ident| $body:expr) => {
        impl<'s, R: IntoTerm<'s>> $trait<R> for Term<'s> {
            type Output = Term<'s>;

            fn $method(
                self,
                rhs: R,
            ) -> Term<'s> {
                let $a = self;
                let $b = rhs.into_term(self.session);
                $body
            }
        }
    };
    (scalar $scalar:ty) => {
        impl<'s> Add<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn add(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) + rhs
            }
        }

        impl<'s> Sub<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn sub(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) - rhs
            }
        }

        impl<'s> Mul<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn mul(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) * rhs
            }
        }

        impl<'s> Div<Term<'s>> for $scalar {
            type Output = Term<'s>;

            fn div(
                self,
                rhs: Term<'s>,
            ) -> Term<'s> {
                self.into_term(rhs.session) / rhs
            }
        }
    };
}

arithmetic!(Add, add, |a, b| a.binary(core::ADD, b));
arithmetic!(Mul, mul, |a, b| a.binary(core::MUL, b));
arithmetic!(Sub, sub, |a, b| a.binary(core::ADD, -b));
arithmetic!(Div, div, |a, b| a.binary(core::MUL, b.pow(-1)));
arithmetic!(scalar i32);
arithmetic!(scalar f64);

impl<'s> Neg for Term<'s> {
    type Output = Term<'s>;

    fn neg(self) -> Term<'s> {
        self.session.int(-1).binary(core::MUL, self)
    }
}

macro_rules! functions {
    ($($(#[$doc:meta])* $name:ident),* $(,)?) => {
        $(
            $(#[$doc])*
            #[must_use]
            pub fn $name(x: Term<'_>) -> Term<'_> {
                x.apply(stringify!($name), &[])
            }
        )*
    };
}

functions! {
    /// Exponential function.
    exp,
    /// Natural logarithm.
    ln,
    /// Sine.
    sin,
    /// Cosine.
    cos,
    /// Tangent.
    tan,
    /// Inverse sine.
    asin,
    /// Inverse cosine.
    acos,
    /// Inverse tangent.
    atan,
    /// Hyperbolic sine.
    sinh,
    /// Hyperbolic cosine.
    cosh,
    /// Hyperbolic tangent.
    tanh,
    /// Square root.
    sqrt,
    /// Absolute value.
    abs,
}

/// The derivative of `f` with respect to the symbol `x`, as an unevaluated
/// request.
#[must_use]
pub fn diff<'s>(
    f: Term<'s>,
    x: Term<'s>,
) -> Term<'s> {
    f.apply("diff", &[x])
}

/// The `n`-th derivative of `f` with respect to `x`.
#[must_use]
pub fn diff_n<'s>(
    f: Term<'s>,
    x: Term<'s>,
    n: usize,
) -> Term<'s> {
    (0..n).fold(f, |acc, _| diff(acc, x))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn symbolic(
        s: &Session,
        text: &str,
    ) -> String {
        let term = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
        let answer = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            answer.reduced,
            "`{text}` was not fully reduced: {}",
            answer.term
        );
        answer.term.to_string()
    }

    /// Checks a derivative in both phases: the symbolic answer, evaluated,
    /// must agree with the numeric answer and with a finite difference of
    /// the original function.
    fn check_derivative(
        text: &str,
        at: f64,
    ) {
        let s = Session::new();
        let f = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
        let x = s.sym("x");
        let request = diff(f, x);
        let closed = s
            .compute(request, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            closed.reduced,
            "d/dx {text} was not reduced: {}",
            closed.term
        );
        let from_closed = closed.term.eval(&[("x", at)]).unwrap_or(f64::NAN);
        let numeric = s
            .compute(request, &Config::new().numeric(1e-9).bind("x", at))
            .unwrap_or_else(|e| panic!("{e}"))
            .value
            .unwrap_or(f64::NAN);
        let h = 1e-5;
        let fd = (f.eval(&[("x", at + h)]).unwrap_or(f64::NAN)
            - f.eval(&[("x", at - h)]).unwrap_or(f64::NAN))
            / (2.0 * h);
        let scale = 1.0 + fd.abs();
        assert!(
            (from_closed - numeric).abs() < 1e-9 * scale,
            "{text}: closed form {from_closed} vs numeric {numeric}"
        );
        assert!(
            (from_closed - fd).abs() < 1e-6 * scale,
            "{text}: {} = {from_closed} but finite difference = {fd}",
            closed.term
        );
    }

    #[test]
    fn operators_build_canonical_terms() {
        let s = Session::new();
        let (x, y) = (s.sym("x"), s.sym("y"));
        assert_eq!((x + y).to_string(), "x + y");
        assert_eq!((x - y).to_string(), "x - y");
        assert_eq!((x / y).to_string(), "x/y");
        assert_eq!((2 * x + 1).to_string(), "2*x + 1");
        assert_eq!((-x).to_string(), "-x");
        assert_eq!(x.pow(2).to_string(), "x^2");
        assert_eq!((1.5 * x).to_string(), "1.5*x");
        assert_eq!(sin(x + y), s.parse("sin(y + x)").unwrap_or(x));
    }

    #[test]
    fn simplification() {
        let s = Session::new();
        assert_eq!(symbolic(&s, "x + x + sin(y)^2 + cos(y)^2"), "2*x + 1");
        assert_eq!(symbolic(&s, "(x^2 * x^3) / x^4"), "x");
        assert_eq!(symbolic(&s, "exp(ln(a)) - a"), "0");
    }

    #[test]
    fn derivatives_symbolic() {
        let s = Session::new();
        assert_eq!(symbolic(&s, "diff(x^3, x)"), "3*x^2");
        assert_eq!(symbolic(&s, "diff(sin(x), x)"), "cos(x)");
        assert_eq!(symbolic(&s, "diff(a*x^2 + b*x + c, x)"), "2*a*x + b");
        assert_eq!(symbolic(&s, "diff(exp(2*x), x)"), "2*exp(2*x)");
        assert_eq!(symbolic(&s, "diff(ln(x), x)"), "1/x");
        assert_eq!(symbolic(&s, "diff(y, x)"), "0");
        assert_eq!(symbolic(&s, "diff(diff(x^4, x), x)"), "12*x^2");
        assert_eq!(symbolic(&s, "diff(sin(x)^2 + cos(x)^2, x)"), "0");
    }

    #[test]
    fn derivatives_agree_across_phases() {
        for (text, at) in [
            ("x^3 - 2*x + 1", 0.7),
            ("sin(x) * cos(x)", 0.3),
            ("exp(x^2) / (1 + x^2)", 0.5),
            ("ln(1 + x^2) * tan(x)", 0.4),
            ("x^x", 1.3),
            ("sqrt(1 + sin(x)^2)", 1.1),
            ("atan(x) + asin(x/2) + acos(x/3)", 0.6),
            ("sinh(x) * cosh(2*x) - tanh(x)", 0.2),
            ("(x + 1)^5 * (x - 2)^(-3)", 0.9),
            ("2^x + x^(1/3)", 1.5),
            ("abs(x) * x", -0.8),
            ("cot(x) + sec(x) * csc(x)", 0.7),
            ("acot(x) + asec(x + 2) + acsc(x + 2)", 0.6),
            ("coth(x) - sech(x) * csch(x)", 0.9),
            ("asinh(x) + acosh(x + 2) + atanh(x / 2)", 0.5),
            ("acoth(x + 2) + asech(x / 2) + acsch(x)", 0.8),
            ("atan2(x, 2) + log(3, x)", 1.4),
            ("x * y * sin(x * y)", 0.35),
        ] {
            if text.contains('y') {
                // `y` is an independent symbol: bind it for evaluation.
                let s = Session::new();
                let f = s.parse(text).unwrap_or_else(|e| panic!("{e}"));
                let request = diff(f, s.sym("x"));
                let closed = s
                    .compute(request, &Config::new())
                    .unwrap_or_else(|e| panic!("{e}"));
                let got = closed
                    .term
                    .eval(&[("x", at), ("y", 2.0)])
                    .unwrap_or(f64::NAN);
                let want = 2.0 * (2.0 * at).sin() + 2.0 * at * 2.0 * (2.0 * at).cos();
                assert!(
                    (got - want).abs() < 1e-12,
                    "{}: {got} vs {want}",
                    closed.term
                );
            } else {
                check_derivative(text, at);
            }
        }
    }

    #[test]
    fn undetermined_functions_stay_unreduced_but_make_progress() {
        let s = Session::new();
        let term = s
            .parse("diff(x * f(x), x)")
            .unwrap_or_else(|e| panic!("{e}"));
        let answer = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(!answer.reduced);
        assert_eq!(
            answer.term.to_string(),
            "x*diff(f(x), x) + f(x)"
        );
    }

    #[test]
    fn numeric_target_without_bindings_fails_cleanly() {
        let s = Session::new();
        let term = s.parse("x + 1").unwrap_or_else(|e| panic!("{e}"));
        let result = s.compute(term, &Config::new().numeric(1e-9));
        assert_eq!(
            result.err(),
            Some(ComputeError::NotNumeric("x + 1".to_owned()))
        );
    }

    #[test]
    fn numeric_differentiation_of_black_box_operators() {
        use crate::graph::Arity;
        use crate::graph::OpDescriptor;
        // An operator with scalar semantics but no derivative rule.
        let blackbox = RuleSet::new("blackbox", |i| {
            i.op(OpDescriptor::new("bump", Arity::Fixed(1))
                .eval(|a| a.first().map_or(f64::NAN, |x| (-x * x).exp())))
                .map(|_| ())
        });
        let config = Config::new().with(blackbox.clone());
        let s =
            Session::with_rules(&[rules::calculus(), blackbox]).unwrap_or_else(|e| panic!("{e}"));
        let term = s
            .parse("diff(bump(x), x)")
            .unwrap_or_else(|e| panic!("{e}"));

        let symbolic = s.compute(term, &config).unwrap_or_else(|e| panic!("{e}"));
        assert!(!symbolic.reduced, "there is no symbolic rule for bump");

        let numeric = s
            .compute(term, &config.numeric(1e-6).bind("x", 0.5))
            .unwrap_or_else(|e| panic!("{e}"));
        let want = -2.0 * 0.5 * (-0.25_f64).exp();
        let got = numeric.value.unwrap_or(f64::NAN);
        assert!((got - want).abs() < 1e-8, "{got} vs {want}");
        assert!(numeric.error.is_some_and(|e| e < 1e-6));
    }

    #[test]
    fn sessions_are_reusable_and_bindings_do_not_stick() {
        let s = Session::new();
        let term = s
            .parse("diff(x^2, x) + y")
            .unwrap_or_else(|e| panic!("{e}"));
        let a = s.compute(
            term,
            &Config::new().numeric(1e-9).bind("x", 1.0).bind("y", 10.0),
        );
        let b = s.compute(
            term,
            &Config::new().numeric(1e-9).bind("x", 2.0).bind("y", 20.0),
        );
        assert_eq!(a.ok().and_then(|a| a.value), Some(12.0));
        assert_eq!(b.ok().and_then(|b| b.value), Some(24.0));
        assert_eq!(symbolic(&s, "diff(x^2, x) + y"), "2*x + y");
    }

    #[test]
    fn assumptions_are_per_request() {
        let s = Session::new();
        let term = s.parse("sqrt(x^2)").unwrap_or_else(|e| panic!("{e}"));
        let plain = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let positive = s
            .compute(term, &Config::new().assume("x", Facts::POSITIVE))
            .unwrap_or_else(|e| panic!("{e}"));
        let again = s
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(plain.term.to_string(), "sqrt(x^2)");
        assert_eq!(positive.term.to_string(), "x");
        assert_eq!(
            again.term.to_string(),
            "sqrt(x^2)",
            "an assumption must not outlive its request"
        );
    }

    #[test]
    fn reduce_then_compile() {
        let s = Session::new();
        let x = s.sym("x");
        let derivative = s
            .simplify(diff(x.pow(3) * sin(x), x))
            .unwrap_or_else(|e| panic!("{e}"));
        let f = derivative.compile(&["x"]).unwrap_or_else(|e| panic!("{e}"));
        let xs: Vec<f64> = (0..100).map(|i| f64::from(i) * 0.1).collect();
        let mut out = vec![0.0; xs.len()];
        f.call_batch(&[&xs], &mut out);
        for (x, got) in xs.iter().zip(&out) {
            let want = 3.0 * x * x * x.sin() + x.powi(3) * x.cos();
            assert!(
                (got - want).abs() < 1e-9 * (1.0 + want.abs()),
                "at {x}: {got} vs {want}"
            );
        }
        // An unreduced request has no numeric semantics.
        assert!(diff(x, x).compile(&["x"]).is_err());
    }

    #[test]
    fn unknown_operators_are_reported() {
        let s = Session::new();
        let x = s.sym("x");
        assert_eq!(
            s.call("frobnicate", &[x]).err(),
            Some(ComputeError::UnknownOperator("frobnicate".to_owned()))
        );
        assert!(s.call("sin", &[x, x]).is_err());
        assert!(s.parse("1 +").is_err());
    }
}
