//! Multivariate series: Taylor and Laurent expansions in several variables.
//!
//! Two request operators, registered by [`install`] in the `calculus` rule
//! set:
//!
//! | operator | value |
//! |---|---|
//! | `mseries(f, list(x1, ..., xn), list(a1, ..., an), N)` | the expansion of `f` about `(a1, ..., an)`, truncated at total degree `N` |
//! | `mtaylor(f, list(x1, ..., xn), list(a1, ..., an), N)` | the same, restricted to Taylor polynomials: every point finite and no negative power; otherwise the request stays unevaluated |
//!
//! The expansion point of a variable may be a number, a symbol (an
//! unspecified parameter) or `oo` / `-oo`; in the latter case the
//! expansion is in powers of `1/x`.
//!
//! # Algorithm
//!
//! The multivariate expansion is *recursive*: `f` is expanded in the first
//! variable by the univariate Laurent engine
//! ([`super::powerseries::laurent_expansion`]), then every coefficient of
//! that expansion is expanded in the next variable, and so on. A term whose
//! earlier exponents sum to `p` can only contribute to total degree `N` if
//! its remaining exponents sum to at most `N - p`, so the power demanded of
//! each variable is bounded by
//!
//! ```text
//! target + sum over the later variables of max(0, -valuation of f there),
//! ```
//!
//! where `valuation` is the univariate Laurent valuation of `f` in that
//! variable (its pole order). This bound is sufficient exactly when `f`
//! has a monomial factor at the point, so the truncation is complete: the
//! returned sum contains *every* monomial of total degree at most `N`. The
//! bound is checked after the expansion; when it was not sufficient, the
//! support is unbounded below and the expansion is directional (see
//! below).
//!
//! The variables are processed in order of increasing univariate valuation
//! (most negative first, ties in the order given). This puts the variables
//! in which `f` has a pole first, which keeps the truncation complete for
//! Laurent expansions with a monomial factor and for expansions in `1/x`.
//!
//! Three further routes cover functions the univariate engine cannot shift
//! directly:
//!
//! * `atan`, `asin` and `acos` are expanded by differentiating (their
//!   derivatives are algebraic), expanding the derivative and integrating
//!   term by term, then adding the value at the expansion point;
//! * a Laurent monomial coefficient is taken over verbatim, which is both
//!   cheaper and avoids a truncation limit of the univariate engine for
//!   very negative exponents;
//! * a power with a numeric non-integer exponent (`sqrt`, `(1+x+y)^(1/2)`)
//!   goes through the univariate engine's binomial series.
//!
//! A power whose exponent is symbolic rather than numeric, such as
//! `(1 + x)^a`, is not a structure the univariate engine expands and the
//! request stays unevaluated. Symbolic *coefficients* and symbolic
//! expansion points (`a`, `exp(a)`, `(x - a)^k`) are fully supported.
//!
//! # Ordering and truncation
//!
//! The expansion point is the origin of the translated variables
//! `u_i = x_i - a_i` (for a point at infinity, `u_i = 1/x_i`). A monomial
//! `u1^e1 * ... * un^en` has **total degree** `e1 + ... + en`; negative
//! powers are allowed and lower the degree. The result is the expansion
//! truncated at total degree `N`: a term is kept exactly when its total
//! degree is at most `N`. For a Taylor expansion (`N >= 0`, all exponents
//! non-negative) this is the usual total-degree truncation and the first
//! omitted term has total degree `N + 1`.
//!
//! For Laurent expansions the same rule is used, and it is well defined
//! whenever `f` has a **monomial factor**: if `f = u^v * g` with `g`
//! analytic and `g(0) != 0`, then the support is bounded below
//! componentwise and only finitely many monomials have total degree at
//! most `N`. For example with `N = 2`
//!
//! * `exp(x)/y^2` is `y^(-2) * (1 + x + x^2/2 + x^3/6 + x^4/24)`: the
//!   terms `x^k y^(-2)` have total degree `k - 2 <= 2`;
//! * `x^(-1) * sin(y)` is `x^(-1) * (y - y^3/6)`.
//!
//! In agreement with the univariate `laurent` operator (`n = 1`), `N` is
//! the highest *power* kept, not the length of the coefficient list.
//!
//! # Directional expansions
//!
//! A function with a pole along a hypersurface through the expansion point
//! (`1/(x^2 - y^2)`, `1/(x + y)`, `(x + y)^(-2)` at the origin) has no
//! Laurent expansion at all: its support is unbounded below and infinitely
//! many monomials have the same total degree. For a *rational* function of
//! that kind the operator returns the **directional** expansion produced by
//! the recursion: the variables are taken in the order described above and
//! the first variable is expanded through the requested order, the others
//! as coefficients. Thus at the origin
//!
//! ```text
//! 1/(x^2 - y^2) = -y^(-2) * (1 + (x/y)^2 + ... + (x/y)^(2m)),  m = floor(N/2),
//! 1/(x + y)^2   = sum_{p=0}^{N} (-1)^p (p+1) x^p y^(-2-p),
//! ```
//!
//! plus a remainder that is the first omitted power of the ratio: for
//! `1/(x^2 - y^2)` the product with `x^2 - y^2` is
//! `1 - (x/y)^(2m+2)`, and for `(x + y)^2` it is `1` minus the first
//! omitted term of `sum_p (-1)^p (p+1) (x/y)^p`. These are the expansions
//! valid where the first variable is small compared to the others; no claim
//! of uniform approximation in a full neighbourhood is made, because no
//! such expansion exists.
//!
//! Any other function whose recursion does not close (an essential
//! singularity or a branch point, e.g. `exp(1/x)`, `exp(x/y)`, `sqrt(x)`)
//! makes the request stay **unevaluated** ([`Outcome::Pass`]) rather than
//! return a wrong answer.
//!
//! # Error term
//!
//! Like `taylor`/`laurent`, the operator returns the truncated part and not
//! an explicit `big_o` node. For a Taylor expansion the remainder is
//! `O(|u|^(N+1))`; for a Laurent expansion with monomial factor `u^v` it is
//! `u^v * O(|u|^(N+1-|v|))`. The caller can add the error term itself.
//!
//! # Resource limits
//!
//! * at most 8 variables and order at most 20;
//! * at most 4096 monomials per result;
//! * at most 800 univariate expansions per request, each of order at most
//!   128. A request that needs a longer univariate series — a very high
//!   power combined with a pole in the *same* variable — stays unevaluated.
//!
//! Exceeding a limit returns `Pass` (unevaluated) with one exception: a
//! rational function with a pole along a hypersurface returns the
//! directional expansion described above.

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use std::collections::HashSet;

use super::powerseries::laurent_expansion;

/// Largest accepted `order` argument.
const MAX_ORDER: usize = 20;
/// Largest number of expansion variables.
const MAX_VARIABLES: usize = 8;
/// Largest number of monomials one expansion may produce.
const MAX_TERMS: usize = 4096;
/// Largest univariate coefficient array the recursion may ask for.
const MAX_COEFFICIENT_ORDER: i64 = 128;
/// Largest number of univariate expansions one request may run.
const MAX_EXPANSIONS: usize = 800;

/// Where one variable is expanded.
#[derive(Clone, Copy)]
enum Point {
    /// A finite point (a number, a symbol or any other term).
    Finite(NodeId),
    /// `+oo`: expand in `1/x`.
    PositiveInfinity,
    /// `-oo`: expand in `1/x` (with the sign of every other power folded
    /// into the coefficient).
    NegativeInfinity,
}

impl Point {
    /// Whether this is a finite expansion point.
    const fn is_finite(self) -> bool {
        matches!(self, Self::Finite(_))
    }
}

/// How the exponents of one expansion variable are printed.
#[derive(Clone, Copy)]
struct Base {
    /// The base to raise to a power: the translated variable `x - a`, or
    /// the variable itself when the translation is the identity.
    node: NodeId,
    /// Whether the exponent is negated (`1/x` expansion).
    reciprocal: bool,
    /// Whether odd powers carry a minus sign (`-oo`).
    signed: bool,
}

/// One parsed `mseries`/`mtaylor` request.
struct Request {
    /// The function to expand.
    body: NodeId,
    /// Expansion variables in the order the recursion processes them.
    variables: Vec<(NodeId, SymbolId)>,
    /// Expansion points, in the same order as `variables`.
    points: Vec<Point>,
    /// Requested total degree.
    order: i64,
}

/// One term of the expansion: an exponent per variable and an exact
/// coefficient (a term of the graph independent of every expansion
/// variable).
struct Term {
    /// Exponents in the request's variable order.
    exponents: Vec<i64>,
    /// The coefficient.
    coefficient: NodeId,
}

/// The result of one expansion.
struct Expansion {
    /// The truncated sum as a term.
    term: NodeId,
    /// Whether every exponent of the sum is a non-negative integer.
    is_taylor: bool,
}

/// Parses the four arguments of a request.
fn request(
    graph: &Graph,
    args: &[NodeId],
    infinity: Option<OpId>,
) -> Option<Request> {
    let &[body, variables, points, order] = args else {
        return None;
    };
    let order = usize::try_from(graph.number_of(order)?.to_i64()?).ok()?;
    if order > MAX_ORDER {
        return None;
    }
    let order = i64::try_from(order).ok()?;
    let variable_nodes = list_items(graph, variables)?;
    let point_nodes = list_items(graph, points)?;
    if variable_nodes.len() != point_nodes.len() || variable_nodes.len() > MAX_VARIABLES {
        return None;
    }
    let mut parsed = Vec::with_capacity(variable_nodes.len());
    for (&variable, &point) in variable_nodes.iter().zip(&point_nodes) {
        parsed.push((
            variable,
            graph.symbol_of(variable)?,
            point_of(graph, point, infinity),
        ));
    }
    // Distinct variables: expanding in the same variable twice would not be
    // a multivariate expansion.
    let mut seen = HashSet::new();
    for &(_, symbol, _) in &parsed {
        if !seen.insert(symbol) {
            return None;
        }
    }
    Some(Request {
        body,
        variables: parsed
            .iter()
            .map(|&(node, symbol, _)| (node, symbol))
            .collect(),
        points: parsed.iter().map(|&(_, _, point)| point).collect(),
        order,
    })
}

/// The children of a `list(...)` argument, or `None` for anything else.
fn list_items(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    Some(graph.children(node).to_vec())
}

/// Classifies an expansion point: `oo`, `-oo` or a finite term.
fn point_of(
    graph: &Graph,
    node: NodeId,
    infinity: Option<OpId>,
) -> Point {
    if let Some(oo) = infinity {
        if graph.op(node) == oo {
            return Point::PositiveInfinity;
        }
        // `-oo` is built as `(-1) * oo`; the numeric factor is folded into
        // one literal by the graph.
        if graph.op(node) == core::MUL
            && let [a, b] = graph.children(node)
        {
            let negative = |n: NodeId| graph.number_of(n).is_some_and(Number::is_negative);
            let infinite = |n: NodeId| graph.op(n) == oo;
            if (infinite(*a) && negative(*b)) || (infinite(*b) && negative(*a)) {
                return Point::NegativeInfinity;
            }
        }
    }
    Point::Finite(node)
}

/// `base^exponent` with an exact integer exponent.
fn power(
    graph: &mut Graph,
    base: NodeId,
    exponent: i64,
) -> NodeId {
    match exponent {
        | 0 => graph.int(1),
        | 1 => base,
        | e => {
            let e = graph.int(e);
            graph.node(core::POW, &[base, e])
        },
    }
}

/// `-node`.
fn negate(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    graph.node(core::MUL, &[minus_one, node])
}

/// `a + b`, without adding zero.
fn add(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    if graph.number_of(a).is_some_and(Number::is_zero) {
        return b;
    }
    if graph.number_of(b).is_some_and(Number::is_zero) {
        return a;
    }
    graph.node(core::ADD, &[a, b])
}

/// `a - b`.
fn sub(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let negated = negate(graph, b);
    add(graph, a, negated)
}

/// `x - a` as a term, or `None` when `a` is the origin (the variable
/// itself is the offset).
fn offset(
    graph: &mut Graph,
    variable: NodeId,
    point: NodeId,
) -> Option<NodeId> {
    if graph.number_of(point).is_some_and(Number::is_zero) {
        return None;
    }
    Some(sub(graph, variable, point))
}

/// Whether `f` is a rational function of its symbols: built from numbers,
/// symbols, `+`, `*` and integer powers only.
fn is_rational(
    graph: &Graph,
    f: NodeId,
) -> bool {
    let mut stack = vec![f];
    let mut seen = HashSet::new();
    while let Some(node) = stack.pop() {
        if !seen.insert(node) {
            continue;
        }
        if graph.as_number(node).is_some() || graph.as_symbol(node).is_some() {
            continue;
        }
        let op = graph.op(node);
        if op == core::ADD || op == core::MUL {
            stack.extend_from_slice(graph.children(node));
            continue;
        }
        if op == core::POW {
            let &[base, exponent] = graph.children(node) else {
                return false;
            };
            if graph.number_of(exponent).and_then(Number::to_i64).is_none() {
                return false;
            }
            stack.push(base);
            continue;
        }
        return false;
    }
    true
}

/// The sum of the exponents of one term.
fn total_degree(exponents: &[i64]) -> i64 {
    exponents
        .iter()
        .fold(0_i64, |acc, &e| acc.saturating_add(e))
}

/// The recursive expansion engine.
struct Expander<'c, 'a> {
    cx: &'c mut Cx<'a>,
    /// The translated variables, in processing order.
    variables: Vec<NodeId>,
    /// Their symbols.
    symbols: Vec<SymbolId>,
    /// For each level, the sum over the remaining variables of the pole
    /// order of the function in that variable. A monomial of the expansion
    /// whose exponent at `level` is `p` has remaining degree at least the
    /// negative of this value, so coefficients above
    /// `target + slack[level]` cannot contribute.
    slack: Vec<i64>,
    /// Whether the function is rational (built from `+`, `*` and integer
    /// powers). Only a rational function with a pole along a hypersurface
    /// through the point has a documented directional expansion.
    rational: bool,
    /// Remaining univariate expansions.
    expansions: usize,
    /// Remaining monomials.
    monomials: usize,
    /// Remaining differentiate-and-integrate fallbacks.
    fallbacks: usize,
}

impl Expander<'_, '_> {
    /// Charges one univariate expansion against the budget.
    fn charge_expansion(&mut self) -> Option<()> {
        self.expansions = self.expansions.checked_sub(1)?;
        Some(())
    }

    /// Charges `count` monomials against the budget.
    fn charge_monomials(
        &mut self,
        count: usize,
    ) -> Option<()> {
        self.monomials = self.monomials.checked_sub(count)?;
        Some(())
    }

    /// A literal zero.
    fn zero(&mut self) -> NodeId {
        self.cx.graph.int(0)
    }

    /// Whether `node` is the literal zero.
    fn is_zero(
        &self,
        node: NodeId,
    ) -> bool {
        self.cx.graph.number_of(node).is_some_and(Number::is_zero)
    }

    /// The expansion of a term that depends on no variable any more.
    fn constant(
        &self,
        f: NodeId,
    ) -> Option<Vec<Term>> {
        let class = self.cx.graph.find(f);
        if self
            .symbols
            .iter()
            .any(|&symbol| self.cx.graph.depends_on(class, symbol))
        {
            return None;
        }
        if self.is_zero(f) {
            return Some(Vec::new());
        }
        let exponents = vec![0_i64; self.symbols.len()];
        Some(vec![Term {
            exponents,
            coefficient: f,
        }])
    }

    /// The index of one of the expansion variables.
    fn index_of(
        &self,
        symbol: SymbolId,
    ) -> Option<usize> {
        self.symbols.iter().position(|&s| s == symbol)
    }

    /// Whether `node` depends on any variable of `level..`.
    fn depends_on_remaining(
        &self,
        node: NodeId,
        level: usize,
    ) -> bool {
        let class = self.cx.graph.find(node);
        self.symbols
            .iter()
            .skip(level)
            .any(|&symbol| self.cx.graph.depends_on(class, symbol))
    }

    /// The exact expansion of `f` when `f` is a Laurent monomial in the
    /// remaining variables: one term, kept exactly when its total degree is
    /// at most `target`.
    ///
    /// This is both an optimisation — most coefficients produced by the
    /// recursion are monomials, and expanding them by the univariate engine
    /// would be wasteful — and a guard against a limit of that engine,
    /// which cannot resolve a coefficient whose valuation lies below its
    /// own truncated window.
    fn monomial(
        &mut self,
        f: NodeId,
        level: usize,
        target: i64,
    ) -> Option<Vec<Term>> {
        let width = self.symbols.len();
        let mut exponents = vec![0_i64; width];
        let mut factors: Vec<NodeId> = Vec::new();
        let mut stack = vec![f];
        while let Some(node) = stack.pop() {
            if self.cx.graph.as_number(node).is_some() {
                factors.push(node);
                continue;
            }
            if let Some(symbol) = self.cx.graph.as_symbol(node) {
                match self.index_of(symbol) {
                    | Some(index) if index >= level => {
                        let slot = exponents.get_mut(index)?;
                        *slot = slot.saturating_add(1);
                    },
                    | Some(_) => return None,
                    | None => factors.push(node),
                }
                continue;
            }
            let op = self.cx.graph.op(node);
            if op == core::MUL {
                stack.extend_from_slice(self.cx.graph.children(node));
                continue;
            }
            if op == core::POW
                && let [base, exponent] = self.cx.graph.children(node)
                && let Some(symbol) = self.cx.graph.as_symbol(*base)
                && let Some(index) = self.index_of(symbol)
            {
                if index < level {
                    return None;
                }
                let power = self
                    .cx
                    .graph
                    .number_of(*exponent)
                    .and_then(Number::to_i64)?;
                let slot = exponents.get_mut(index)?;
                *slot = slot.saturating_add(power);
                continue;
            }
            if self.depends_on_remaining(node, level) {
                return None;
            }
            factors.push(node);
        }
        if total_degree(&exponents) > target {
            return Some(Vec::new());
        }
        let coefficient = match factors.as_slice() {
            | [] => self.cx.graph.int(1),
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, &factors),
        };
        if self.is_zero(coefficient) {
            return Some(Vec::new());
        }
        Some(vec![Term {
            exponents,
            coefficient,
        }])
    }

    /// Rewrites `base^exponent` with an exponent that is not a literal
    /// number (a parameter, as in `(1 + x)^a`) as `exp(exponent * ln(base))`
    /// and expands that. The univariate power-series engine only knows
    /// integer and rational exponents.

    /// Expands `f` by differentiating and integrating term by term.
    ///
    /// The univariate engine knows the Taylor series of `atan`, `asin` and
    /// `acos` only at zero, because their general shift is not a product of
    /// elementary series. Their derivatives are algebraic, so an expansion
    /// at a non-zero constant term is recovered from
    /// `f = f(0, ...) + integral(df/dx dx)`: the derivative is expanded by
    /// the ordinary recursion and integrated term by term, and the value at
    /// the expansion point is expanded in the remaining variables.
    fn antiderivative(
        &mut self,
        f: NodeId,
        level: usize,
        target: i64,
    ) -> Option<Vec<Term>> {
        if self.fallbacks == 0 {
            return None;
        }
        let &variable = self.variables.get(level)?;
        let &symbol = self.symbols.get(level)?;
        let derivative = super::derivative(self.cx.graph, f, variable)?;
        if derivative == f {
            return None;
        }
        // The value at the expansion point must be finite; a singular value
        // (for example `atan(1/x)` at zero) has no Taylor expansion.
        let zero = self.cx.graph.int(0);
        let value = self.cx.graph.substitute(f, variable, zero);
        if self.cx.graph.depends_on(self.cx.graph.find(value), symbol) {
            return None;
        }
        let samples: Vec<f64> = vec![0.3, 0.8, 1.4];
        for &sample in &samples {
            let mut env = crate::graph::Env::numeric(0.0);
            for (index, &other) in self
                .symbols
                .iter()
                .skip(level.saturating_add(1))
                .enumerate()
            {
                #[allow(clippy::cast_precision_loss)]
                let offset = index as f64 * 0.2;
                env.bind(other, sample + offset);
            }
            if let Some(number) = self.cx.graph.eval(value, &env)
                && number.is_infinite()
            {
                return None;
            }
        }
        self.fallbacks = self.fallbacks.saturating_sub(1);
        let inner = self.expand(derivative, level, target.saturating_sub(1));
        self.fallbacks = self.fallbacks.saturating_add(1);
        let inner = inner?;
        let mut terms = Vec::with_capacity(inner.len());
        for term in inner {
            let exponent = term.exponents.get(level).copied().unwrap_or(0);
            if exponent == -1 {
                return None;
            }
            let lifted = exponent.checked_add(1)?;
            let rest = total_degree(&term.exponents).saturating_sub(exponent);
            if lifted.saturating_add(rest) > target {
                continue;
            }
            let coefficient = if lifted == 1 {
                term.coefficient
            } else {
                let scale = Number::fraction(1, lifted)?;
                let scale = self.cx.graph.num(scale);
                self.cx.graph.node(core::MUL, &[scale, term.coefficient])
            };
            let mut term = term;
            if let Some(slot) = term.exponents.get_mut(level) {
                *slot = lifted;
            }
            term.coefficient = coefficient;
            terms.push(term);
        }
        let constant = self.expand(value, level.saturating_add(1), target)?;
        terms.extend(constant);
        self.charge_monomials(terms.len())?;
        Some(terms)
    }

    /// The expansion of `f` in the variables `level..`, keeping monomials of
    /// total degree at most `target`.
    fn expand(
        &mut self,
        f: NodeId,
        level: usize,
        target: i64,
    ) -> Option<Vec<Term>> {
        if level >= self.variables.len() {
            return self.constant(f);
        }
        if self.is_zero(f) {
            return Some(Vec::new());
        }
        // A Laurent monomial is exact; expanding it through the univariate
        // engine would only cost time.
        if let Some(terms) = self.monomial(f, level, target) {
            self.charge_monomials(terms.len())?;
            return Some(terms);
        }
        let &variable = self.variables.get(level)?;
        let &symbol = self.symbols.get(level)?;
        if !self.cx.graph.depends_on(self.cx.graph.find(f), symbol) {
            return self.expand(f, level.saturating_add(1), target);
        }
        let last = level.saturating_add(1) == self.variables.len();
        // A term with exponent `p` at this level and remaining degree `q`
        // enters the result when `p + q <= target`; `q` is at least the
        // negative of the pole order in the remaining variables, so `p`
        // never needs to exceed `target + slack`.
        let slack = self.slack.get(level).copied().unwrap_or(0);
        let mut bound = target.saturating_add(slack);
        if bound > MAX_COEFFICIENT_ORDER {
            return None;
        }
        loop {
            self.charge_expansion()?;
            let zero = self.zero();
            // The univariate engine indexes its coefficient list from zero,
            // so it cannot be asked for a negative high-water mark; asking
            // for `0` and filtering afterwards is equivalent here.
            let call_order = bound.max(0);
            let Some((valuation, coefficients)) =
                laurent_expansion(self.cx, f, variable, zero, call_order)
            else {
                return self.antiderivative(f, level, target);
            };
            if valuation > bound {
                // Every coefficient through `bound` vanishes. A leading term
                // further out can still reach total degree `target` only
                // through negative powers of the remaining variables.
                if last && valuation > target {
                    return Some(Vec::new());
                }
                if valuation > MAX_COEFFICIENT_ORDER {
                    return None;
                }
                bound = valuation;
                continue;
            }
            let mut parts: Vec<Term> = Vec::new();
            let mut minimum: Option<i64> = None;
            for (offset, &coefficient) in coefficients.iter().enumerate() {
                if self.is_zero(coefficient) {
                    continue;
                }
                let step = i64::try_from(offset).ok()?;
                let exponent = valuation.saturating_add(step);
                let inner_target = target.saturating_sub(exponent);
                let inner = self.expand(coefficient, level.saturating_add(1), inner_target)?;
                for mut term in inner {
                    let rest = total_degree(&term.exponents);
                    let degree = rest.saturating_add(exponent);
                    if degree > target {
                        continue;
                    }
                    if let Some(slot) = term.exponents.get_mut(level) {
                        *slot = exponent;
                    }
                    // `minimum` tracks the smallest *remaining* degree: a
                    // term with remaining degree `q` can still appear for an
                    // exponent as large as `target - q`.
                    minimum = Some(minimum.map_or(rest, |low| low.min(rest)));
                    parts.push(term);
                }
            }
            self.charge_monomials(parts.len())?;
            let needed = minimum.map(|low| target.saturating_sub(low));
            if needed.is_some_and(|wanted| wanted > bound) {
                // The bound derived from the pole orders was not sufficient,
                // so the expansion does not have a monomial factor: finitely
                // many monomials of total degree at most `target` do not
                // exist. A rational function still has the directional
                // expansion computed above; anything else has no Laurent
                // expansion at this point and stays unevaluated.
                if !self.rational {
                    return None;
                }
            }
            return Some(parts);
        }
    }
}

/// Expands the request and builds the resulting term.
fn expand(
    cx: &mut Cx<'_>,
    request: &Request,
) -> Option<Expansion> {
    // Translate every variable to the origin: `x -> a + u`, `x -> 1/u` at
    // `+oo`, `x -> -1/u` at `-oo` with `u > 0`.
    let mut symbols = Vec::with_capacity(request.variables.len());
    let mut nodes = Vec::with_capacity(request.variables.len());
    for _ in &request.variables {
        let symbol = cx.graph.interner_mut().fresh_symbol("u");
        nodes.push(cx.graph.symbol_node(symbol));
        symbols.push(symbol);
    }
    let mut body = request.body;
    for (index, &(variable, _)) in request.variables.iter().enumerate() {
        let &point = request.points.get(index)?;
        let &u = nodes.get(index)?;
        let replacement = match point {
            | Point::Finite(a) => {
                if cx.graph.number_of(a).is_some_and(Number::is_zero) {
                    u
                } else {
                    add(cx.graph, a, u)
                }
            },
            | Point::PositiveInfinity | Point::NegativeInfinity => {
                let &symbol = symbols.get(index)?;
                cx.graph.assume(symbol, Facts::POSITIVE);
                let inverse = power(cx.graph, u, -1);
                if matches!(point, Point::NegativeInfinity) {
                    negate(cx.graph, inverse)
                } else {
                    inverse
                }
            },
        };
        body = cx.graph.substitute(body, variable, replacement);
    }

    // Process the variables with the most negative univariate valuation
    // first, so that the pole direction is expanded first. A variable whose
    // univariate expansion does not exist contributes a neutral bound.
    let mut valuations = Vec::with_capacity(nodes.len());
    for &u in &nodes {
        let zero = cx.graph.int(0);
        valuations.push(
            laurent_expansion(cx, body, u, zero, 0).map_or(0_i64, |(valuation, _)| valuation),
        );
    }
    let mut order: Vec<usize> = (0..nodes.len()).collect();
    order.sort_by_key(|&index| (valuations.get(index).copied().unwrap_or(0), index));
    let mut bases = Vec::with_capacity(nodes.len());
    for &index in &order {
        let &(variable, _) = request.variables.get(index)?;
        let &point = request.points.get(index)?;
        bases.push(match point {
            | Point::Finite(a) => Base {
                node: offset(cx.graph, variable, a).unwrap_or(variable),
                reciprocal: false,
                signed: false,
            },
            | Point::PositiveInfinity | Point::NegativeInfinity => Base {
                node: variable,
                reciprocal: true,
                signed: matches!(point, Point::NegativeInfinity),
            },
        });
    }
    let translated: Vec<NodeId> = order
        .iter()
        .filter_map(|&index| nodes.get(index).copied())
        .collect();
    let translated_symbols: Vec<SymbolId> = order
        .iter()
        .filter_map(|&index| symbols.get(index).copied())
        .collect();
    // `slack[level]` bounds from below the total degree contributed by the
    // variables after `level`: every exponent is at least the univariate
    // valuation of the function in that variable.
    let mut slack = vec![0_i64; order.len()];
    let mut running = 0_i64;
    for position in (0..order.len()).rev() {
        if let Some(slot) = slack.get_mut(position) {
            *slot = running;
        }
        let pole = order
            .get(position)
            .and_then(|&index| valuations.get(index))
            .map_or(0, |&valuation| valuation.min(0).saturating_neg());
        running = running.saturating_add(pole);
    }
    let rational = is_rational(cx.graph, body);

    let mut expander = Expander {
        cx,
        variables: translated,
        symbols: translated_symbols,
        slack,
        rational,
        expansions: MAX_EXPANSIONS,
        monomials: MAX_TERMS,
        fallbacks: 2,
    };
    let terms = expander.expand(body, 0, request.order)?;
    let is_taylor = terms
        .iter()
        .all(|term| term.exponents.iter().all(|&e| e >= 0));
    let term = build(expander.cx.graph, &bases, &terms);
    Some(Expansion { term, is_taylor })
}

/// Builds the truncated sum `sum coefficient * product base^exponent`.
fn build(
    graph: &mut Graph,
    bases: &[Base],
    terms: &[Term],
) -> NodeId {
    let mut sum = Vec::with_capacity(terms.len());
    for term in terms {
        let mut coefficient = term.coefficient;
        let mut factors = Vec::new();
        for (index, base) in bases.iter().enumerate() {
            let Some(&exponent) = term.exponents.get(index) else {
                continue;
            };
            if exponent == 0 {
                continue;
            }
            if base.signed && exponent.rem_euclid(2) == 1 {
                coefficient = negate(graph, coefficient);
            }
            let printed = if base.reciprocal {
                exponent.saturating_neg()
            } else {
                exponent
            };
            factors.push(power(graph, base.node, printed));
        }
        if factors.is_empty() {
            sum.push(coefficient);
        } else {
            if !graph.number_of(coefficient).is_some_and(Number::is_one) {
                factors.insert(0, coefficient);
            }
            sum.push(graph.node(core::MUL, &factors));
        }
    }
    match sum.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &sum),
    }
}

/// Kernel for `mseries` and `mtaylor`.
pub(super) struct MultiSeries {
    mseries: OpId,
    mtaylor: OpId,
    infinity: Option<OpId>,
}

impl Kernel for MultiSeries {
    fn ops(&self) -> Vec<OpId> {
        vec![self.mseries, self.mtaylor]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let taylor_only = op == self.mtaylor;
        if op != self.mseries && !taylor_only {
            return Outcome::Pass;
        }
        let args = cx.graph.children(node).to_vec();
        let Some(request) = request(cx.graph, &args, self.infinity) else {
            return Outcome::Pass;
        };
        if taylor_only && !request.points.iter().all(|point| point.is_finite()) {
            return Outcome::Pass;
        }
        let Some(expansion) = expand(cx, &request) else {
            return Outcome::Pass;
        };
        if taylor_only && !expansion.is_taylor {
            return Outcome::Pass;
        }
        Outcome::Pinned(expansion.term)
    }
}

/// Registers `mseries` and `mtaylor` and their kernel.
///
/// `mseries(f, list(x1, ..., xn), list(a1, ..., an), N)` returns the
/// expansion of `f` about `(a1, ..., an)` truncated at total degree `N`
/// (see the [module documentation](self) for the precise ordering and for
/// the directional convention used at poles along a hypersurface).
/// `mtaylor` is the Taylor special case and stays unevaluated when the
/// expansion would contain a negative power or when a point is infinite.
///
/// # Errors
/// Fails when the operators cannot be registered.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let mseries = i.op(OpDescriptor::new("mseries", Arity::Fixed(4))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    let mtaylor = i.op(OpDescriptor::new("mtaylor", Arity::Fixed(4))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    // A series is a way of writing the function it expands, so both requests
    // keep the dimension of their first argument.
    i.dim(mseries, DimRule::Same(0));
    i.dim(mtaylor, DimRule::Same(0));
    let infinity = i.graph().ops().lookup("oo");
    i.kernel(
        "calculus/multiseries",
        Tier::Reduce,
        MultiSeries {
            mseries,
            mtaylor,
            infinity,
        },
    );
    Ok(())
}
#[cfg(test)]
mod tests {
    use crate::rules::calculus;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    /// The printed, fully reduced expansion.
    fn run(src: &str) -> String {
        simplify(&[calculus()], src)
    }

    #[test]
    fn taylor_total_degree() {
        assert_eq!(
            run("mseries(exp(x + y), list(x, y), list(0, 0), 2)"),
            "1/2*x^2 + x*y + 1/2*y^2 + x + y + 1"
        );
        assert_eq!(
            run("mseries(1/(1 - x - y), list(x, y), list(0, 0), 2)"),
            "x^2 + 2*x*y + y^2 + x + y + 1"
        );
        assert_eq!(
            run("mseries(ln(1 + x + y), list(x, y), list(0, 0), 2)"),
            "x - 1/2*x^2 - x*y - 1/2*y^2 + y"
        );
        assert_eq!(
            run("mtaylor(exp(x + y), list(x, y), list(0, 0), 2)"),
            run("mseries(exp(x + y), list(x, y), list(0, 0), 2)")
        );
        assert_eq!(run("mseries(exp(x + y), list(x, y), list(0, 0), 0)"), "1");
        assert_eq!(
            run("mseries(exp(x + y), list(x, y), list(0, 0), 1)"),
            "x + y + 1"
        );
    }

    #[test]
    fn laurent_monomial_factor() {
        assert_eq!(
            run("mseries(1/(x*y), list(x, y), list(0, 0), 2)"),
            "1/(x*y)"
        );
        assert_eq!(
            run("mseries(exp(x)/y^2, list(x, y), list(0, 0), 2)"),
            "1/24*x^4/y^2 + 1/6*x^3/y^2 + 1/2*x^2/y^2 + x/y^2 + 1/y^2"
        );
        assert_eq!(
            run("mseries(x^(-1)*sin(y), list(x, y), list(0, 0), 3)"),
            "y/x - 1/6*y^3/x"
        );
    }

    #[test]
    fn directional_expansions() {
        assert_eq!(
            run("mseries(1/(x^2 - y^2), list(x, y), list(0, 0), 4)"),
            "-x^4/y^6 - x^2/y^4 - 1/y^2"
        );
        assert_eq!(
            run("mseries((x + y)^(-2), list(x, y), list(0, 0), 2)"),
            "3*x^2/y^4 - 2*x/y^3 + 1/y^2"
        );
    }

    #[test]
    fn points_and_univariate_agreement() {
        assert_eq!(
            run("mseries(x/(x - 1), list(x), list(oo), 3)"),
            "1 + 1/x + 1/x^2 + 1/x^3"
        );
        assert_eq!(
            run("mseries(1/(x + 1), list(x), list(-oo), 2)"),
            "1/x - 1/x^2"
        );
        assert_eq!(
            run("mseries(exp(x), list(x), list(0), 4)"),
            run("taylor(exp(x), x, 0, 4)")
        );
        assert_eq!(
            run("mseries(exp(x)/x^2, list(x), list(0), 1)"),
            run("laurent(exp(x)/x^2, x, 0, 1)")
        );
        assert_eq!(
            run("mseries(1/(x*y), list(x, y), list(1, 2), 2)"),
            "1/8*(y - 2)^2 + 1/4*(y - 2)*(x - 1) + 1/2*(x - 1)^2 - 1/4*(y - 2) - 1/2*(x - 1) + 1/2"
        );
    }

    #[test]
    fn unevaluated_requests_stay() {
        for src in [
            "mseries(exp(1/x), list(x), list(0), 2)",
            "mseries(sin(1/x), list(x), list(0), 2)",
            "mseries(exp(x/y), list(x, y), list(0, 0), 2)",
            "mseries(sqrt(x), list(x), list(0), 2)",
            "mtaylor(1/(x*y), list(x, y), list(0, 0), 2)",
        ] {
            let (text, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(!reduced, "`{src}` should stay unevaluated: {text}");
        }
    }
}
