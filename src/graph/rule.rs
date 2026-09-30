//! Rules and rule sets.
//!
//! Every algorithm in rssn is an *identity transformation*: it replaces a
//! term by an equal one. Two kinds of rules express that:
//!
//! * a [`Rewrite`] is declarative — a pattern, a replacement and guards,
//!   written as text: `"sin(?x)^2 + cos(?x)^2 => 1"`;
//! * a [`Kernel`] is procedural — arbitrary code that looks at a node and
//!   returns an equal term (or a numeric witness). Risch integration,
//!   Buchberger's algorithm and an adaptive quadrature are all kernels.
//!
//! Both are [`Rule`]s, registered in the same [`RuleSet`]s and driven by the
//! same scheduler. That is what fuses symbolic and numeric computation: a
//! quadrature kernel and an antiderivative table compete for the same
//! `integral` node.

use std::fmt;
use std::sync::Arc;

use super::extract::Extractor;
use super::extract::SizeCost;
use super::facts::Facts;
use super::id::NodeId;
use super::id::OpId;
use super::id::SymbolId;
use super::op::OpConflict;
use super::op::OpDescriptor;
use super::op::OpFlags;
use super::pattern::Match;
use super::pattern::ParseError;
use super::pattern::Parser;
use super::pattern::Pat;
use super::pattern::VarNames;
use super::schedule::Budget;
use super::schedule::Engine;
use super::schedule::Saturate;
use super::store::Ball;
use super::store::Graph;
use super::window::WindowPass;

/// Scheduling tier of a rule.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Tier {
    /// Reduces heavy operators (derivatives, integrals, solves). Runs first
    /// and to a fixpoint so that requests are resolved before the search
    /// space grows.
    Reduce,
    /// Simplifications that never enlarge a term. Also eligible for
    /// destructive application inside tree windows.
    Normalize,
    /// Structure-changing identities (expansion, factoring, angle
    /// addition). Budgeted and subject to back-off.
    Explore,
}

/// An arbitrary predicate over a substitution, for [`Guard::Custom`].
pub type Predicate = dyn Fn(&Graph, &[NodeId]) -> bool + Send + Sync;

/// A side condition of a [`Rewrite`]. Arguments are pattern variable
/// indices.
#[derive(Clone)]
pub enum Guard {
    /// `free_of(?a, ?x)`: `?a` does not depend on the symbol `?x`.
    FreeOf(u32, u32),
    /// `number(?a)`: `?a` equals a literal number.
    IsNumber(u32),
    /// `integer(?a)`: `?a` is certainly an integer: a literal one, or a
    /// term known to be integer-valued.
    IsInteger(u32),
    /// `symbol(?a)`: `?a` is a symbol.
    IsSymbol(u32),
    /// `nonzero(?a)`: `?a` is certainly not zero.
    NonZero(u32),
    /// `positive(?a)`: `?a` is certainly greater than zero.
    Positive(u32),
    /// `negative(?a)`: `?a` is certainly less than zero.
    Negative(u32),
    /// `nonnegative(?a)`: `?a` is certainly real and not less than zero.
    NonNegative(u32),
    /// `real(?a)`: `?a` is certainly a real number.
    Real(u32),
    /// `same(?a, ?b)`: both are known equal.
    Same(u32, u32),
    /// Negation of another guard, written `!guard(...)`.
    Not(Box<Self>),
    /// An arbitrary predicate over the substitution.
    Custom(Arc<Predicate>),
}

impl fmt::Debug for Guard {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::FreeOf(a, x) => write!(f, "free_of(?{a}, ?{x})"),
            | Self::IsNumber(a) => write!(f, "number(?{a})"),
            | Self::IsInteger(a) => write!(f, "integer(?{a})"),
            | Self::IsSymbol(a) => write!(f, "symbol(?{a})"),
            | Self::NonZero(a) => write!(f, "nonzero(?{a})"),
            | Self::Positive(a) => write!(f, "positive(?{a})"),
            | Self::Negative(a) => write!(f, "negative(?{a})"),
            | Self::NonNegative(a) => write!(f, "nonnegative(?{a})"),
            | Self::Real(a) => write!(f, "real(?{a})"),
            | Self::Same(a, b) => write!(f, "same(?{a}, ?{b})"),
            | Self::Not(g) => write!(f, "!{g:?}"),
            | Self::Custom(_) => write!(f, "custom(..)"),
        }
    }
}

impl Graph {
    /// The symbol a class is known to equal, if any.
    #[must_use]
    pub fn symbol_of(
        &self,
        node: NodeId,
    ) -> Option<SymbolId> {
        self.as_symbol(node)
            .or_else(|| self.enodes(self.find(node)).find_map(|n| self.as_symbol(n)))
    }
}

impl Guard {
    /// Evaluates the guard under a substitution.
    #[must_use]
    pub fn holds(
        &self,
        graph: &Graph,
        subst: &[NodeId],
    ) -> bool {
        let get = |v: u32| subst.get(v as usize).copied().filter(|n| !n.is_none());
        match self {
            | Self::FreeOf(a, x) => match (get(*a), get(*x).and_then(|x| graph.symbol_of(x))) {
                | (Some(a), Some(x)) => !graph.depends_on(graph.find(a), x),
                | _ => false,
            },
            | Self::IsNumber(a) => get(*a).is_some_and(|a| graph.number_of(a).is_some()),
            | Self::IsInteger(a) => get(*a).is_some_and(|a| graph.facts(a).has(Facts::INTEGER)),
            | Self::IsSymbol(a) => get(*a).is_some_and(|a| graph.symbol_of(a).is_some()),
            | Self::NonZero(a) => get(*a).is_some_and(|a| graph.facts(a).has(Facts::NONZERO)),
            | Self::Positive(a) => get(*a).is_some_and(|a| graph.facts(a).has(Facts::POSITIVE)),
            | Self::Negative(a) => get(*a).is_some_and(|a| graph.facts(a).has(Facts::NEGATIVE)),
            | Self::NonNegative(a) => {
                get(*a).is_some_and(|a| graph.facts(a).has(Facts::NONNEGATIVE))
            },
            | Self::Real(a) => get(*a).is_some_and(|a| graph.facts(a).has(Facts::REAL)),
            | Self::Same(a, b) => match (get(*a), get(*b)) {
                | (Some(a), Some(b)) => graph.same(a, b),
                | _ => false,
            },
            | Self::Not(inner) => !inner.holds(graph, subst),
            | Self::Custom(f) => f(graph, subst),
        }
    }
}

/// A declarative identity `lhs = rhs`, applied left to right.
#[derive(Clone, Debug)]
pub struct Rewrite {
    /// Pattern to find. Its root is always an operator application.
    pub lhs: Pat,
    /// Replacement.
    pub rhs: Pat,
    /// Side conditions, all of which must hold.
    pub guards: Vec<Guard>,
    /// Number of pattern variables.
    pub nvars: usize,
}

impl Rewrite {
    /// Operator at the root of the left-hand side.
    #[must_use]
    pub fn trigger(&self) -> OpId {
        self.lhs.root_op().unwrap_or(OpId::NONE)
    }

    /// Whether every guard holds for `subst`.
    #[must_use]
    pub fn admits(
        &self,
        graph: &Graph,
        subst: &[NodeId],
    ) -> bool {
        self.guards.iter().all(|g| g.holds(graph, subst))
    }

    /// Checks the guards for a match and, if they hold, builds the
    /// right-hand side, re-attaching the rest of an associative-commutative
    /// match.
    pub fn build(
        &self,
        graph: &mut Graph,
        m: &Match,
    ) -> Option<NodeId> {
        let subst = m.materialize(graph);
        if !self.admits(graph, &subst) {
            return None;
        }
        let rhs = self.rhs.instantiate(graph, &subst)?;
        if m.rest.is_empty() {
            return Some(rhs);
        }
        let mut children = Vec::with_capacity(m.rest.len().saturating_add(1));
        children.push(rhs);
        children.extend_from_slice(&m.rest);
        graph.try_node(self.trigger(), &children)
    }
}

/// What a [`Kernel`] found out about a node.
#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    /// Nothing; the kernel does not apply.
    Pass,
    /// The node is identical to this term.
    Equal(NodeId),
    /// The node is identical to this term *and* this term is the requested
    /// form: extraction must return it verbatim (see
    /// [`Graph::pin`](super::store::Graph::pin)).
    Pinned(NodeId),
    /// The node's value lies in this enclosure under the current bindings.
    Approx(Ball),
}

/// Numeric environment of a run: variable bindings and tolerance.
#[derive(Clone, Debug, Default)]
pub struct Env {
    bindings: Vec<(SymbolId, f64)>,
    /// Requested absolute tolerance for numeric kernels.
    pub tolerance: f64,
    /// Whether the run wants a numeric answer. Numeric kernels stay idle
    /// otherwise.
    pub numeric: bool,
    /// How many engine runs this one is nested in: `0` for a run started
    /// by the user, more for runs a kernel started through
    /// [`Cx::simplify`].
    pub depth: u8,
}

impl Env {
    /// An environment with no bindings that does not ask for numerics.
    #[must_use]
    pub const fn symbolic() -> Self {
        Self {
            bindings: Vec::new(),
            tolerance: 1e-10,
            numeric: false,
            depth: 0,
        }
    }

    /// An environment asking for numeric results to `tolerance`.
    #[must_use]
    pub const fn numeric(tolerance: f64) -> Self {
        Self {
            bindings: Vec::new(),
            tolerance,
            numeric: true,
            depth: 0,
        }
    }

    /// Binds `symbol` to `value`, replacing an earlier binding.
    pub fn bind(
        &mut self,
        symbol: SymbolId,
        value: f64,
    ) {
        match self.bindings.binary_search_by_key(&symbol, |b| b.0) {
            | Ok(i) => {
                if let Some(slot) = self.bindings.get_mut(i) {
                    slot.1 = value;
                }
            },
            | Err(i) => self.bindings.insert(i, (symbol, value)),
        }
    }

    /// The value bound to `symbol`.
    #[must_use]
    pub fn value(
        &self,
        symbol: SymbolId,
    ) -> Option<f64> {
        let i = self.bindings.binary_search_by_key(&symbol, |b| b.0).ok()?;
        self.bindings.get(i).map(|b| b.1)
    }

    /// All bindings, sorted by symbol.
    #[must_use]
    pub fn bindings(&self) -> &[(SymbolId, f64)] {
        &self.bindings
    }
}

/// Context handed to a [`Kernel`].
#[derive(Debug)]
pub struct Cx<'a> {
    /// The graph; kernels may add nodes but must not union.
    pub graph: &'a mut Graph,
    /// Numeric environment of the run.
    pub env: &'a Env,
    /// The engine that is running this kernel.
    pub engine: &'a Engine,
}

/// Deepest nesting of engine runs started from kernels.
const MAX_NESTING: u8 = 3;

impl Cx<'_> {
    /// Simplifies `term` with the running program and returns the best
    /// form found.
    ///
    /// Algorithms that need a zero test or a canonical form in the middle
    /// of their work — integration by substitution, limits, series — call
    /// this instead of reimplementing simplification. The nested run is
    /// symbolic, has a small budget, and is refused (returning `term`
    /// unchanged) beyond a fixed nesting depth so that kernels calling
    /// each other cannot recurse without bound.
    pub fn simplify(
        &mut self,
        term: NodeId,
    ) -> NodeId {
        if self.env.depth >= MAX_NESTING {
            return term;
        }
        let mut env = Env::symbolic();
        env.depth = self.env.depth.saturating_add(1);
        let budget = Budget { max_iterations: 6, max_nodes: 4_000, patience: 2, ..Budget::default() };
        self.engine.run(self.graph, &[term], &env, &Saturate, &budget);
        Extractor::new(self.graph, &[term], &SizeCost).build(self.graph, term).unwrap_or(term)
    }

    /// Whether `term` simplifies to the literal zero.
    pub fn is_zero(
        &mut self,
        term: NodeId,
    ) -> bool {
        let simplified = self.simplify(term);
        self.graph.number_of(simplified).is_some_and(super::number::Number::is_zero)
    }
}

/// A procedural identity transformation: a *reduction kernel*.
pub trait Kernel: Send + Sync {
    /// Operators whose nodes this kernel wants to see. A list containing
    /// [`OpId::NONE`] means every operator.
    fn ops(&self) -> Vec<OpId>;

    /// Inspects `node` and reports what it equals.
    ///
    /// The returned term must be *identical* in value to `node` for every
    /// assignment of its free symbols ([`Outcome::Equal`]), or enclose its
    /// value under the bindings of `cx.env` ([`Outcome::Approx`]).
    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome;

    /// Whether the kernel should see a node again after the graph changed.
    ///
    /// Cheap kernels whose answer depends on what is known about the
    /// children (constant folding) return `true`. Expensive ones that fully
    /// decide a node on first sight keep the default.
    fn revisit(&self) -> bool {
        false
    }
}

/// The executable part of a rule.
#[derive(Clone)]
pub enum Action {
    /// Declarative rewrite.
    Rewrite(Rewrite),
    /// Procedural kernel.
    Kernel(Arc<dyn Kernel>),
}

/// A named, tiered identity transformation.
#[derive(Clone)]
pub struct Rule {
    /// Name for reports and debugging.
    pub name: Arc<str>,
    /// Scheduling tier.
    pub tier: Tier,
    /// What the rule does.
    pub action: Action,
}

impl fmt::Debug for Rule {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "Rule({:?}, {:?})", self.name, self.tier)
    }
}

/// Error raised while installing a rule set into a graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RuleError {
    /// Rule text could not be parsed.
    Parse {
        /// The rule text.
        rule: String,
        /// The parser's complaint.
        error: ParseError,
    },
    /// An operator clashes with an existing registration.
    Op(OpConflict),
    /// The rule is syntactically fine but not a valid identity schema.
    Invalid {
        /// The rule text.
        rule: String,
        /// Why it was rejected.
        reason: &'static str,
    },
}

impl fmt::Display for RuleError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Parse { rule, error } => write!(f, "in rule `{rule}`: {error}"),
            | Self::Op(conflict) => write!(f, "{conflict}"),
            | Self::Invalid { rule, reason } => write!(f, "rule `{rule}` rejected: {reason}"),
        }
    }
}

impl std::error::Error for RuleError {}

impl From<OpConflict> for RuleError {
    fn from(value: OpConflict) -> Self {
        Self::Op(value)
    }
}

/// The installed form of one or more rule sets: rules and window passes
/// bound to the operator ids of one particular graph.
#[derive(Clone, Default)]
pub struct Program {
    /// All rules, in installation order.
    pub rules: Vec<Rule>,
    /// Window passes, in installation order.
    pub passes: Vec<Arc<dyn WindowPass>>,
    installed: Vec<Arc<str>>,
}

impl fmt::Debug for Program {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "Program({:?}: {} rules, {} passes)",
            self.installed,
            self.rules.len(),
            self.passes.len()
        )
    }
}

impl Program {
    /// Names of the rule sets installed so far, dependencies included.
    #[must_use]
    pub fn installed(&self) -> &[Arc<str>] {
        &self.installed
    }
}

/// Collects the operators, rules and passes of a rule set while it is
/// installed into a particular graph.
pub struct Installer<'g> {
    graph: &'g mut Graph,
    program: &'g mut Program,
}

impl<'g> Installer<'g> {
    /// Starts an installation into `graph`, appending to `program`.
    pub const fn new(
        graph: &'g mut Graph,
        program: &'g mut Program,
    ) -> Self {
        Self { graph, program }
    }

    /// The graph being installed into.
    pub const fn graph(&mut self) -> &mut Graph {
        self.graph
    }

    /// Registers (or re-finds) an operator.
    ///
    /// # Errors
    /// Fails when the name is taken by an incompatible operator.
    pub fn op(
        &mut self,
        desc: OpDescriptor,
    ) -> Result<OpId, RuleError> {
        Ok(self.graph.ops_mut().register(desc)?)
    }

    /// Adds rewrites written as text, one per entry:
    ///
    /// ```text
    /// name: lhs => rhs
    /// name: lhs => rhs if guard(?a), !guard(?b, ?c)
    /// name: lhs <=> rhs          (both directions)
    /// ```
    ///
    /// # Errors
    /// Fails on the first entry that does not parse or is not a valid
    /// identity schema (for example a right-hand side using a variable the
    /// left-hand side does not bind).
    pub fn rewrites(
        &mut self,
        tier: Tier,
        texts: &[&str],
    ) -> Result<(), RuleError> {
        for text in texts {
            self.rewrite(tier, text)?;
        }
        Ok(())
    }

    fn rewrite(
        &mut self,
        tier: Tier,
        text: &str,
    ) -> Result<(), RuleError> {
        let parse_err = |error: ParseError| RuleError::Parse {
            rule: text.to_owned(),
            error,
        };
        let invalid = |reason: &'static str| RuleError::Invalid {
            rule: text.to_owned(),
            reason,
        };
        let (name, body) = text
            .split_once(':')
            .ok_or_else(|| invalid("missing `name:` prefix"))?;
        let name = name.trim();

        let mut vars = VarNames::default();
        let mut parser = Parser::new(body, self.graph, &mut vars, false);
        let lhs = parser.expr().map_err(parse_err)?;
        let both_ways = if parser.eat("<=>") {
            true
        } else if parser.eat("=>") {
            false
        } else {
            return Err(parse_err(ParseError {
                message: "expected `=>` or `<=>`".to_owned(),
                offset: parser.offset(),
            }));
        };
        let rhs = parser.expr().map_err(parse_err)?;
        let mut guards = Vec::new();
        if parser.ident() == Some("if") {
            loop {
                guards.push(parse_guard(&mut parser).map_err(parse_err)?);
                if !parser.eat(",") {
                    break;
                }
            }
        }
        if !parser.at_end() {
            return Err(parse_err(ParseError {
                message: "unexpected trailing input".to_owned(),
                offset: parser.offset(),
            }));
        }
        let nvars = vars.len();

        let mut directions = vec![(lhs.clone(), rhs.clone(), name.to_owned())];
        if both_ways {
            directions.push((rhs, lhs, format!("{name}-rev")));
        }
        for (lhs, rhs, name) in directions {
            let Some(root) = lhs.root_op() else {
                return Err(invalid("left-hand side must be an operator application"));
            };
            if self.graph.ops().get(root).flags.has(OpFlags::LEAF) {
                return Err(invalid("left-hand side must be an operator application"));
            }
            let (mut bound, mut used) = (Vec::new(), Vec::new());
            lhs.vars(&mut bound);
            rhs.vars(&mut used);
            if used.iter().any(|v| !bound.contains(v)) {
                return Err(invalid(
                    "right-hand side uses a variable the left-hand side does not bind",
                ));
            }
            self.program.rules.push(Rule {
                name: Arc::from(name),
                tier,
                action: Action::Rewrite(Rewrite {
                    lhs,
                    rhs,
                    guards: guards.clone(),
                    nvars,
                }),
            });
        }
        Ok(())
    }

    /// Adds a procedural kernel.
    pub fn kernel(
        &mut self,
        name: &str,
        tier: Tier,
        kernel: impl Kernel + 'static,
    ) {
        self.program.rules.push(Rule {
            name: Arc::from(name),
            tier,
            action: Action::Kernel(Arc::new(kernel)),
        });
    }

    /// Adds a window pass, run during local optimisation of tree windows.
    pub fn pass(
        &mut self,
        pass: impl WindowPass + 'static,
    ) {
        self.program.passes.push(Arc::new(pass));
    }
}

fn parse_guard(parser: &mut Parser<'_>) -> Result<Guard, ParseError> {
    let fail = |parser: &Parser<'_>, message: &str| ParseError {
        message: message.to_owned(),
        offset: parser.offset(),
    };
    if parser.eat("!") {
        return Ok(Guard::Not(Box::new(parse_guard(parser)?)));
    }
    let name = parser
        .ident()
        .ok_or_else(|| fail(parser, "expected a guard name"))?;
    if !parser.eat("(") {
        return Err(fail(parser, "expected `(` after guard name"));
    }
    let mut args = Vec::new();
    loop {
        match parser.expr()? {
            | Pat::Var(v) => args.push(v),
            | _ => return Err(fail(parser, "guard arguments must be pattern variables")),
        }
        if parser.eat(")") {
            break;
        }
        if !parser.eat(",") {
            return Err(fail(parser, "expected `,` or `)`"));
        }
    }
    match (name, args.as_slice()) {
        | ("free_of", [a, x]) => Ok(Guard::FreeOf(*a, *x)),
        | ("same", [a, b]) => Ok(Guard::Same(*a, *b)),
        | ("number", [a]) => Ok(Guard::IsNumber(*a)),
        | ("integer", [a]) => Ok(Guard::IsInteger(*a)),
        | ("symbol", [a]) => Ok(Guard::IsSymbol(*a)),
        | ("nonzero", [a]) => Ok(Guard::NonZero(*a)),
        | ("positive", [a]) => Ok(Guard::Positive(*a)),
        | ("negative", [a]) => Ok(Guard::Negative(*a)),
        | ("nonnegative", [a]) => Ok(Guard::NonNegative(*a)),
        | ("real", [a]) => Ok(Guard::Real(*a)),
        | _ => Err(fail(parser, "unknown guard or wrong number of arguments")),
    }
}

type InstallFn = dyn Fn(&mut Installer<'_>) -> Result<(), RuleError> + Send + Sync;

/// A named group of operators and rules for one mathematical domain.
///
/// A rule set is a *recipe*: installing it into a graph registers its
/// operators there and yields rules bound to that graph's operator ids.
#[derive(Clone)]
pub struct RuleSet {
    name: Arc<str>,
    deps: Vec<Self>,
    install: Arc<InstallFn>,
}

impl fmt::Debug for RuleSet {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "RuleSet({:?})", self.name)
    }
}

impl RuleSet {
    /// Defines a rule set.
    pub fn new(
        name: &str,
        install: impl Fn(&mut Installer<'_>) -> Result<(), RuleError> + Send + Sync + 'static,
    ) -> Self {
        Self {
            name: Arc::from(name),
            deps: Vec::new(),
            install: Arc::new(install),
        }
    }

    /// Declares that `dep` must be installed before this set.
    #[must_use]
    pub fn needs(
        mut self,
        dep: Self,
    ) -> Self {
        self.deps.push(dep);
        self
    }

    /// The set's name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Installs this set and its dependencies (each at most once, in
    /// dependency order) into `graph`, appending to `program`.
    ///
    /// # Errors
    /// Propagates the first [`RuleError`].
    pub fn install_into(
        &self,
        graph: &mut Graph,
        program: &mut Program,
    ) -> Result<(), RuleError> {
        if program.installed.contains(&self.name) {
            return Ok(());
        }
        program.installed.push(Arc::clone(&self.name));
        for dep in &self.deps {
            dep.install_into(graph, program)?;
        }
        (self.install)(&mut Installer::new(graph, program))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;

    fn trig() -> RuleSet {
        RuleSet::new("trig", |i| {
            i.op(OpDescriptor::new("sin", Arity::Fixed(1)))?;
            i.op(OpDescriptor::new("cos", Arity::Fixed(1)))?;
            i.rewrites(Tier::Normalize, &["pythagoras: sin(?x)^2 + cos(?x)^2 => 1"])
        })
    }

    fn install(set: &RuleSet) -> Result<(Graph, Vec<Rule>), RuleError> {
        let mut g = Graph::new();
        let mut program = Program::default();
        set.install_into(&mut g, &mut program)?;
        Ok((g, program.rules))
    }

    #[test]
    fn install_registers_ops_and_rules() {
        let (g, rules) = install(&trig()).unwrap_or_else(|e| panic!("{e}"));
        assert!(g.ops().lookup("sin").is_some());
        assert_eq!(rules.len(), 1);
        assert_eq!(rules.first().map(|r| &*r.name), Some("pythagoras"));
    }

    #[test]
    fn dependencies_install_once() {
        let a = trig();
        let b = RuleSet::new("b", |_| Ok(())).needs(a.clone());
        let c = RuleSet::new("c", |_| Ok(())).needs(a).needs(b);
        let (_, rules) = install(&c).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(rules.len(), 1);
    }

    #[test]
    fn bidirectional_rules_expand_to_two() {
        let set = RuleSet::new("s", |i| {
            i.rewrites(Tier::Explore, &["dist: ?a*(?b + ?c) <=> ?a*?b + ?a*?c"])
        });
        let (_, rules) = install(&set).unwrap_or_else(|e| panic!("{e}"));
        let names: Vec<&str> = rules.iter().map(|r| &*r.name).collect();
        assert_eq!(names, vec!["dist", "dist-rev"]);
    }

    #[test]
    fn invalid_rules_are_rejected() {
        let cases = [
            "no-name ?a + 0 => ?a",
            "unbound: ?a + 0 => ?b",
            "var-root: ?a => ?a + 0",
            "arrow: ?a + 0 -> ?a",
            "guard: ?a + 0 => ?a if wibble(?a)",
            "guard-arity: ?a + 0 => ?a if number(?a, ?a)",
            "trailing: ?a + 0 => ?a )",
            "unknown-op: frob(?a) => ?a",
        ];
        for case in cases {
            let set = RuleSet::new("s", move |i| i.rewrites(Tier::Normalize, &[case]));
            assert!(install(&set).is_err(), "`{case}` should be rejected");
        }
        // A reversed rule whose new left-hand side is a bare variable.
        let set = RuleSet::new("s", |i| {
            i.rewrites(Tier::Normalize, &["rev: ?a + 0 <=> ?a"])
        });
        assert!(install(&set).is_err());
    }

    #[test]
    fn guards_evaluate() {
        let set = RuleSet::new("s", |i| {
            i.rewrites(
                Tier::Normalize,
                &["g: ?a * ?x => ?a if free_of(?a, ?x), !number(?x), nonzero(?a)"],
            )
        });
        let (mut g, rules) = install(&set).unwrap_or_else(|e| panic!("{e}"));
        let Some(Action::Rewrite(rw)) = rules.first().map(|r| r.action.clone()) else {
            panic!("expected a rewrite");
        };
        let (two, x, y) = (g.int(2), g.sym("x"), g.sym("y"));
        let zero = g.int(0);
        assert!(rw.admits(&g, &[two, x]));
        assert!(!rw.admits(&g, &[x, x]), "x is not free of x");
        assert!(!rw.admits(&g, &[y, x]), "y is not known to be nonzero");
        assert!(!rw.admits(&g, &[zero, x]));
        assert!(
            !rw.admits(&g, &[two, two]),
            "2 is not a symbol, so free_of cannot be established"
        );
    }

    #[test]
    fn env_bindings() {
        let mut g = Graph::new();
        let (x, y) = (g.interner_mut().symbol("x"), g.interner_mut().symbol("y"));
        let mut env = Env::numeric(1e-6);
        env.bind(y, 2.0);
        env.bind(x, 1.0);
        env.bind(y, 3.0);
        assert_eq!(env.value(x), Some(1.0));
        assert_eq!(env.value(y), Some(3.0));
        assert_eq!(env.bindings().len(), 2);
        assert!(env.numeric);
        assert!(!Env::symbolic().numeric);
    }

    #[test]
    fn rewrite_build_reattaches_rest() {
        let (mut g, rules) = install(&trig()).unwrap_or_else(|e| panic!("{e}"));
        let Some(Action::Rewrite(rw)) = rules.first().map(|r| r.action.clone()) else {
            panic!("expected a rewrite");
        };
        let node = g.parse("sin(t)^2 + cos(t)^2 + k").unwrap_or(NodeId::NONE);
        let matches = rw.lhs.matches(&g, node, rw.nvars, 8);
        let built = matches.first().and_then(|m| rw.build(&mut g, m));
        assert_eq!(built, g.parse("1 + k").ok());
    }
}
