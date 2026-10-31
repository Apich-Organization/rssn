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
//!   goes through the univariate engine's binomial series when the base has
//!   a non-zero constant term, i.e. when the power is analytic.
//!
//! A non-analytic fractional leading power (`sqrt(x*y)`, `x^(4/3)`) is a
//! Puiseux series; it is handled by the ramification described in the next
//! section. A power whose exponent is symbolic rather than numeric, such
//! as `(1 + x)^a`, is not a structure either engine expands and the
//! request stays unevaluated. Symbolic *coefficients* and symbolic
//! expansion points (`a`, `exp(a)`, `(x - a)^k`) are fully supported.
//!
//! # Puiseux expansions: fractional monomial factors
//!
//! A genuine Puiseux series is a series in fractional powers, for example
//! `sqrt(x*y) = x^(1/2) y^(1/2)` or `x^(4/3)`. The integer recursion cannot
//! represent those directly, so when it fails and the body contains a
//! fractional exponent the operator extracts the leading monomial with
//! rational exponents. Choose `Q`, a multiple of the denominators of the
//! exponent literals (the smallest `Q` that works is searched), and
//! substitute `u_i = v_i^Q` for every expansion variable. Then
//!
//! ```text
//! F(v_1, ..., v_n) = f(v_1^Q, ..., v_n^Q)
//! ```
//!
//! is expanded by the same integer recursion. Because the substitution is
//! exact, a term `c * v_1^k1 * ... * v_n^kn` of `F` is the term
//! `c * u_1^(k1/Q) * ... * u_n^(kn/Q)` of `f`. In particular the term of
//! `F` with least total degree is the leading monomial
//! `m = ∏ u_i^(e_i)`, `e_i = k_i/Q`, and writing `f = m * unit` the
//! remaining terms, divided by `m`, are the total-degree expansion of
//! `unit`; the operator multiplies them back. Every monomial of `f` of
//! total degree at most `N` has `v`-degree at most `Q*N`, so truncating
//! `F` at total degree `Q*N` produces every monomial of total degree at
//! most `N`, with exact rational exponents.
//!
//! For example:
//!
//! ```text
//! mseries(sqrt(x*y), list(x, y), list(0, 0), 4) = x^(1/2)*y^(1/2)
//! mseries(x^(4/3), list(x), list(0), 3)         = x^(4/3)
//! mseries(x^(1/2)*exp(y), list(x, y), list(0, 0), 3)
//!   = 1/2*x^(1/2)*y^2 + x^(1/2)*y + x^(1/2)
//! ```
//!
//! `Q` is enlarged until the ramified function has an integer-exponent
//! expansion (a nested root such as `((x^(1/2))^(1/2))` needs `Q = 4`),
//! subject to the resource limits below. A fractional power of an analytic
//! unit such as `(1 + x + y)^(1/2)` is already a Taylor series and is
//! expanded by the ordinary recursion; only a fractional leading power
//! needs ramification. `mtaylor` remains the Taylor special case: a
//! fractional or negative exponent makes it stay unevaluated.
//!
//! # When no expansion exists
//!
//! Two kinds of point admit no Laurent/Puiseux expansion. There the request
//! stays unevaluated ([`Outcome::Pass`]) rather than returning a truncated
//! series that is not the function:
//!
//! * **A branch locus through the point.** `(x^2 + y^2)^(1/2)` at the
//!   origin has no Puiseux expansion: its leading behaviour depends on the
//!   direction (`|x|` along `y = 0`, `sqrt(2)|x|` along `x = y`), while a
//!   Puiseux series has a single leading monomial. The operator checks two
//!   things. First, the ramified expansion must close: its support must be
//!   bounded below, which fails here once enough terms are seen (there are
//!   terms of total degree one for every power of `x/y`). Second, the
//!   expansion is recomputed with the variables in the reverse order: a
//!   complete total-degree truncation is unique, and here the two orders
//!   disagree (one gives `y` and the other `x` at order one), so the
//!   truncation is order-dependent. Either check refuses the request. The
//!   order-agreement check is used for every non-rational expansion; a
//!   rational function is exempt because its documented directional
//!   expansion is intentionally order-dependent. The same applies to
//!   `sqrt(x^2 - y^2)`.
//! * **An essential singularity.** `exp(1/x)`, `sin(1/x)` and `exp(x/y)`
//!   have no Laurent or Puiseux expansion at the singular point: the
//!   function does not tend to a finite or infinite limit and no finite
//!   truncation approximates it uniformly. The operator recognises the
//!   pattern (an `exp`/`sin`/`cos`/`sinh`/`cosh` whose argument has a pole
//!   in an expansion variable) and returns unevaluated; it never presents
//!   the partial polynomial of an iterated expansion as the function.
//!
//! `ln(x)` at `0` is a logarithmic singularity and also stays unevaluated.
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
//! A coefficient of the result never mentions one of the operator's internal
//! translated variables. The dependence of a term on a variable is tested by
//! *value*, and a body such as `(x + 1)/(x + 1)` is equal to a constant even
//! though its term mentions `x`; the constant is then used in canonical
//! form, so the answer is `1` rather than an expression in the internal
//! variables.
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
//! Any other function whose recursion does not close — an essential
//! singularity or a branch locus, both described above — makes the request
//! stay **unevaluated** ([`Outcome::Pass`]) rather than return a wrong
//! answer.
//!
//! # Error term
//!
//! Like `taylor`/`laurent`, the operator returns the truncated part and not
//! an explicit `big_o` node. For a Taylor expansion the remainder is
//! `O(|u|^(N+1))`; for a Laurent expansion with monomial factor `u^v` it is
//! `u^v * O(|u|^(N+1-|v|))`. A Puiseux expansion is the same expansion in
//! the ramified variables, so its remainder is the first omitted term in
//! the total order; the returned sum contains every term of total degree at
//! most `N` and no claim is made about the next one. The caller can add the
//! error term itself.
//!
//! # Resource limits
//!
//! * at most 8 variables and order at most 20;
//! * at most 4096 monomials per result;
//! * at most 800 univariate expansions per request, each of order at most
//!   128. A request that needs a longer univariate series — a very high
//!   power combined with a pole in the *same* variable — stays unevaluated;
//! * a ramification index at most 64 and ramified total degree at most 128;
//!   at most 8 successive ramification indices are tried. A Puiseux
//!   expansion that needs a finer ramification stays unevaluated.
//!
//! Exceeding a limit returns `Pass` (unevaluated) with one exception: a
//! rational function with a pole along a hypersurface returns the
//! directional expansion described above.
//!
//! The limits bound the result, not the running time. The cost grows quickly
//! with the order because every coefficient is expanded through the
//! univariate engine: a two-variable order-20 Taylor expansion has only 231
//! terms but takes minutes, whereas the same request at order 6 is
//! instantaneous.

use crate::graph::Arity;
use crate::graph::ClassId;
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
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
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
/// Largest ramification index `Q` used for a Puiseux expansion.
const MAX_RAMIFICATION: i64 = 64;
/// Largest total degree `Q * N` asked of the ramified expansion.
const MAX_RAMIFIED_DEGREE: i64 = 128;
/// Number of successive ramification indices tried.
const MAX_RAMIFICATION_ATTEMPTS: i64 = 8;

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

/// One term of the integer-exponent expansion produced by the recursion.
struct Term {
    /// Exponents in the request's variable order.
    exponents: Vec<i64>,
    /// The coefficient.
    coefficient: NodeId,
}

/// One term of the returned expansion. The exponents are rational: they are
/// integers for a Taylor/Laurent expansion and fractions with denominator
/// `Q` for a Puiseux expansion.
struct OutputTerm {
    /// Exponents in the request's variable order.
    exponents: Vec<BigRational>,
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

/// Whether `node` syntactically mentions one of `symbols`.
///
/// [`Graph::depends_on`] tracks *value* dependence, and its free-symbol set
/// is an over-approximation that shrinks when two equal classes are merged.
/// So a node can still mention an expansion variable without depending on
/// it: `(u + 1)/(u + 1)` is equal to the constant one, but its concrete
/// term mentions `u`. A coefficient of the returned expansion must not
/// mention an internal expansion variable, so this structural test detects
/// that case; the caller then canonicalises the node with
/// [`Cx::simplify`] instead of putting the raw node into the answer.
fn mentions(
    graph: &Graph,
    node: NodeId,
    symbols: &[SymbolId],
) -> bool {
    let mut stack = vec![node];
    let mut seen = HashSet::new();
    while let Some(current) = stack.pop() {
        if !seen.insert(current) {
            continue;
        }
        if graph.op(current) == core::SYM {
            if let Some(symbol) = graph.as_symbol(current)
                && symbols.contains(&symbol)
            {
                return true;
            }
            continue;
        }
        stack.extend_from_slice(graph.children(current));
    }
    false
}

/// The sum of the exponents of one term.
fn total_degree(exponents: &[i64]) -> i64 {
    exponents
        .iter()
        .fold(0_i64, |acc, &e| acc.saturating_add(e))
}

/// The exact rational value of a numeric node.
fn rational_of(
    graph: &Graph,
    node: NodeId,
) -> Option<BigRational> {
    graph.number_of(node)?.to_rational()
}

/// The exponent as an `i64` when it is an integer.
fn integer_exponent(exponent: &BigRational) -> Option<i64> {
    if exponent.is_integer() {
        exponent.to_integer().to_i64()
    } else {
        None
    }
}

/// `n^exponent` with an exact rational exponent.
fn power_fraction(
    graph: &mut Graph,
    base: NodeId,
    exponent: &BigRational,
) -> NodeId {
    if let Some(e) = integer_exponent(exponent) {
        return power(graph, base, e);
    }
    let e = graph.num(Number::rat(exponent.clone()));
    graph.node(core::POW, &[base, e])
}

/// The least common multiple of two positive `i64`s, saturating at
/// [`MAX_RAMIFICATION`].
fn lcm(
    a: i64,
    b: i64,
) -> i64 {
    let (mut x, mut y) = (a.saturating_abs().max(1), b.saturating_abs().max(1));
    while y != 0 {
        let remainder = x.checked_rem(y).unwrap_or(0);
        x = y;
        y = remainder;
    }
    a.saturating_div(x).saturating_mul(b).max(0)
}

/// The least common multiple of the denominators of every fractional power
/// in `f`, or `None` when an exponent is not an exact number or the
/// ramification would exceed [`MAX_RAMIFICATION`].
///
/// `sqrt(u)` is the same term as `u^(1/2)`, but the extractor prefers the
/// shorter spelling, so the concrete argument of a request often still holds
/// a `sqrt` node; it contributes the denominator two.
fn ramification_denominator(
    graph: &Graph,
    f: NodeId,
) -> Option<i64> {
    let sqrt = graph.ops().lookup("sqrt");
    let mut stack = vec![f];
    let mut seen = HashSet::new();
    let mut denominator = 1_i64;
    while let Some(node) = stack.pop() {
        if !seen.insert(node) {
            continue;
        }
        stack.extend_from_slice(graph.children(node));
        let op = graph.op(node);
        if Some(op) == sqrt {
            denominator = lcm(denominator, 2);
            if denominator > MAX_RAMIFICATION {
                return None;
            }
            continue;
        }
        if op != core::POW {
            continue;
        }
        let &[_, exponent] = graph.children(node) else {
            continue;
        };
        let rational = rational_of(graph, exponent)?;
        if rational.is_integer() {
            continue;
        }
        let q = rational.denom().to_i64()?;
        denominator = lcm(denominator, q);
        if denominator > MAX_RAMIFICATION {
            return None;
        }
    }
    Some(denominator)
}

/// Whether `f` has an essential singularity at the expansion point: an
/// `exp`, `sin`, `cos`, `sinh` or `cosh` whose argument has a pole in one of
/// the expansion variables. Such a function has no Laurent or Puiseux
/// expansion there, so the operator must not return a truncated one.
fn has_essential_singularity(
    cx: &mut Cx<'_>,
    f: NodeId,
    variables: &[NodeId],
) -> bool {
    let names = ["exp", "sin", "cos", "sinh", "cosh"];
    let ops: Vec<OpId> = names
        .iter()
        .filter_map(|name| cx.graph.ops().lookup(name))
        .collect();
    if ops.is_empty() {
        return false;
    }
    let mut stack = vec![f];
    let mut seen = HashSet::new();
    while let Some(node) = stack.pop() {
        if !seen.insert(node) {
            continue;
        }
        let children = cx.graph.children(node).to_vec();
        stack.extend_from_slice(&children);
        if !ops.contains(&cx.graph.op(node)) {
            continue;
        }
        let Some(&argument) = children.first() else {
            continue;
        };
        for &variable in variables {
            let zero = cx.graph.int(0);
            // Only a genuine *pole* of the argument matters. A fractional
            // valuation (`exp(sqrt(x))`) is a Puiseux case and must be
            // attempted; an unexpandable argument is left to the ordinary
            // recursion instead of being mislabelled.
            if let Some((valuation, _)) = laurent_expansion(cx, argument, variable, zero, 0)
                && valuation < 0
            {
                return true;
            }
        }
    }
    false
}

/// Converts integer-exponent terms to the returned rational form, dividing
/// every exponent by `divisor` (the ramification index).
fn output_terms(
    terms: Vec<Term>,
    divisor: i64,
) -> Vec<OutputTerm> {
    let denominator = BigInt::from(divisor.max(1));
    terms
        .into_iter()
        .map(|term| {
            OutputTerm {
                exponents: term
                    .exponents
                    .into_iter()
                    .map(|k| BigRational::new(BigInt::from(k), denominator.clone()))
                    .collect(),
                coefficient: term.coefficient,
            }
        })
        .collect()
}

/// One expansion term written in the request's variable order, for comparing
/// two runs that processed the variables differently.
type CanonicalTerm = (Vec<BigRational>, ClassId);

/// Rewrites the terms of a run with variable order `permutation` into the
/// request's variable order and sorts them.
fn canonical_terms(
    graph: &Graph,
    terms: &[Term],
    permutation: &[usize],
    divisor: i64,
) -> Option<Vec<CanonicalTerm>> {
    let denominator = BigInt::from(divisor.max(1));
    let mut out: Vec<CanonicalTerm> = Vec::with_capacity(terms.len());
    for term in terms {
        let mut exponents = vec![BigRational::zero(); permutation.len()];
        for (position, &original) in permutation.iter().enumerate() {
            let k = term.exponents.get(position).copied().unwrap_or(0);
            let slot = exponents.get_mut(original)?;
            *slot = BigRational::new(BigInt::from(k), denominator.clone());
        }
        out.push((exponents, graph.find(term.coefficient)));
    }
    out.sort_by(|(a, _), (b, _)| a.cmp(b));
    Some(out)
}

/// Whether two ramified runs, which processed the variables in different
/// orders, produced the same expansion.
///
/// A complete total-degree truncation is unique, so a disagreement means at
/// least one of the runs was order-dependent and therefore incomplete: the
/// function has no Laurent/Puiseux expansion with a monomial factor at the
/// point and the request must stay unevaluated.
fn same_expansion(
    graph: &Graph,
    first: &[Term],
    first_order: &[usize],
    second: &[Term],
    second_order: &[usize],
    divisor: i64,
) -> bool {
    let Some(a) = canonical_terms(graph, first, first_order, divisor) else {
        return false;
    };
    let Some(b) = canonical_terms(graph, second, second_order, divisor) else {
        return false;
    };
    a == b
}

/// Expands `body` with the variables taken in `permutation` order (indices
/// into `nodes`).
fn ordered_terms(
    cx: &mut Cx<'_>,
    body: NodeId,
    permutation: &[usize],
    nodes: &[NodeId],
    symbols: &[SymbolId],
    target: i64,
    rational: bool,
) -> Option<Vec<Term>> {
    let mut variables = Vec::with_capacity(permutation.len());
    let mut variable_symbols = Vec::with_capacity(permutation.len());
    for &index in permutation {
        variables.push(*nodes.get(index)?);
        variable_symbols.push(*symbols.get(index)?);
    }
    let mut valuations = Vec::with_capacity(variables.len());
    for &variable in &variables {
        let zero = cx.graph.int(0);
        valuations.push(
            laurent_expansion(cx, body, variable, zero, 0)
                .map_or(0_i64, |(valuation, _)| valuation),
        );
    }
    run_engine(
        cx,
        body,
        &variables,
        &variable_symbols,
        target,
        rational,
        &valuations,
    )
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
        &mut self,
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
        let mut coefficient = f;
        if mentions(self.cx.graph, coefficient, &self.symbols) {
            // The class does not depend on any expansion variable, but the
            // concrete term still mentions one: it is equal to a term that
            // does not, as in `(u + 1)/(u + 1)`. Canonicalise it instead of
            // returning the internal variable to the caller.
            coefficient = self.cx.simplify(coefficient);
            if mentions(self.cx.graph, coefficient, &self.symbols) {
                return None;
            }
        }
        if self.is_zero(coefficient) {
            return Some(Vec::new());
        }
        let exponents = vec![0_i64; self.symbols.len()];
        Some(vec![Term {
            exponents,
            coefficient,
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
            if let Some(number) = self.cx.graph.as_number(node).cloned() {
                // The class may be a number even when the node is a
                // compound term merged with one (`(u + 1)/(u + 1)`). Push
                // the canonical number so the answer cannot print the
                // internal expansion variable.
                factors.push(self.cx.graph.num(number));
                continue;
            }
            if let Some(symbol) = self.cx.graph.as_symbol(node) {
                match self.index_of(symbol) {
                    | Some(index) if index >= level => {
                        let slot = exponents.get_mut(index)?;
                        *slot = slot.saturating_add(1);
                    },
                    | Some(_) => return None,
                    | None => factors.push(self.cx.graph.symbol_node(symbol)),
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
            // The class does not depend on the remaining variables, but the
            // concrete term may still mention an expansion variable (see
            // [`mentions`]); a factor of the answer must not.
            let mut factor = node;
            if mentions(self.cx.graph, factor, &self.symbols) {
                factor = self.cx.simplify(factor);
                if mentions(self.cx.graph, factor, &self.symbols) {
                    return None;
                }
            }
            factors.push(factor);
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

/// Expands `f` in the (already translated) `variables` with the integer
/// recursion.
///
/// `valuations[i]` is the univariate valuation of `f` in `variables[i]` (in
/// the processing order); it feeds the truncation bound, and a variable
/// whose univariate expansion does not exist contributes `0`.
fn run_engine(
    cx: &mut Cx<'_>,
    body: NodeId,
    variables: &[NodeId],
    symbols: &[SymbolId],
    order: i64,
    rational: bool,
    valuations: &[i64],
) -> Option<Vec<Term>> {
    if variables.len() != symbols.len() || variables.len() != valuations.len() {
        return None;
    }
    // `slack[level]` bounds from below the total degree contributed by the
    // variables after `level`: every exponent is at least the univariate
    // valuation of the function in that variable.
    let mut slack = vec![0_i64; variables.len()];
    let mut running = 0_i64;
    for position in (0..variables.len()).rev() {
        if let Some(slot) = slack.get_mut(position) {
            *slot = running;
        }
        let pole = valuations
            .get(position)
            .copied()
            .unwrap_or(0)
            .min(0)
            .saturating_neg();
        running = running.saturating_add(pole);
    }
    let mut expander = Expander {
        cx,
        variables: variables.to_vec(),
        symbols: symbols.to_vec(),
        slack,
        rational,
        expansions: MAX_EXPANSIONS,
        monomials: MAX_TERMS,
        fallbacks: 2,
    };
    expander.expand(body, 0, order)
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
            | Point::Finite(a) => {
                Base {
                    node: offset(cx.graph, variable, a).unwrap_or(variable),
                    reciprocal: false,
                    signed: false,
                }
            },
            | Point::PositiveInfinity | Point::NegativeInfinity => {
                Base {
                    node: variable,
                    reciprocal: true,
                    signed: matches!(point, Point::NegativeInfinity),
                }
            },
        });
    }
    let translated: Vec<NodeId> = order
        .iter()
        .filter_map(|&index| nodes.get(index).copied())
        .collect();

    // The ordinary Taylor/Laurent recursion.
    let rational = is_rational(cx.graph, body);
    if let Some(terms) = ordered_terms(cx, body, &order, &nodes, &symbols, request.order, rational)
    {
        // A complete total-degree truncation does not depend on the order in
        // which the variables are processed. Recomputing with the reversed
        // order detects a singularity that one order hides: for
        // `(x^2 + y^2)^(1/2)` at the origin one order returns `y` and the
        // other `x` at order one. A rational function is exempt because its
        // documented directional expansion is intentionally order-dependent.
        if nodes.len() > 1 && !rational {
            let reversed: Vec<usize> = order.iter().rev().copied().collect();
            let second =
                ordered_terms(cx, body, &reversed, &nodes, &symbols, request.order, false)?;
            if !same_expansion(cx.graph, &terms, &order, &second, &reversed, 1) {
                return None;
            }
        }
        let terms = output_terms(terms, 1);
        return finish(cx.graph, &bases, terms, &symbols);
    }

    // An essential singularity has no Laurent or Puiseux expansion; identify
    // it explicitly instead of trying ramified series that cannot exist.
    if has_essential_singularity(cx, body, &translated) {
        return None;
    }

    // Puiseux: ramify `u_i = v_i^Q` so that the function has integer
    // exponents in `v`, expand it with the ordinary recursion, then map each
    // exponent back as `k/Q`. Because the substitution is exact, this is the
    // monomial-factor extraction of the module documentation: the term of
    // least `v`-degree is the leading monomial `m`, and the remaining terms
    // are the total-degree expansion of `f/m`.
    let denominator = ramification_denominator(cx.graph, body)?;
    if denominator <= 1 {
        return None;
    }
    for multiplier in 1..=MAX_RAMIFICATION_ATTEMPTS {
        let Some(q) = denominator.checked_mul(multiplier) else {
            break;
        };
        if q > MAX_RAMIFICATION {
            break;
        }
        let Some(target) = q.checked_mul(request.order) else {
            break;
        };
        if target > MAX_RAMIFIED_DEGREE {
            break;
        }
        let mut ramified = body;
        for &u in &translated {
            let replacement = power(cx.graph, u, q);
            ramified = cx.graph.substitute(ramified, u, replacement);
        }
        let Some(first) = ordered_terms(cx, ramified, &order, &nodes, &symbols, target, false)
        else {
            continue;
        };
        // A complete total-degree truncation does not depend on the order in
        // which the variables are processed. Recomputing with the reversed
        // order detects singularities that a single order hides: for
        // `(x^2 + y^2)^(1/2)` at the origin one order returns `y` and the
        // other `x` at order one, so the two runs disagree and the request
        // stays unevaluated rather than return either half.
        if nodes.len() > 1 {
            let reversed: Vec<usize> = order.iter().rev().copied().collect();
            let second = ordered_terms(cx, ramified, &reversed, &nodes, &symbols, target, false)?;
            if !same_expansion(cx.graph, &first, &order, &second, &reversed, q) {
                return None;
            }
        }
        let terms = output_terms(first, q);
        return finish(cx.graph, &bases, terms, &symbols);
    }
    None
}

/// Filters the zero terms, computes whether the sum is a Taylor polynomial
/// and builds the printed sum.
///
/// Returns `None` when a coefficient still mentions an internal expansion
/// variable: no such term is a valid answer, and returning it would leak the
/// variables the operator introduced.
fn finish(
    graph: &mut Graph,
    bases: &[Base],
    terms: Vec<OutputTerm>,
    symbols: &[SymbolId],
) -> Option<Expansion> {
    let terms: Vec<OutputTerm> = terms
        .into_iter()
        .filter(|term| {
            !graph
                .number_of(term.coefficient)
                .is_some_and(Number::is_zero)
        })
        .collect();
    if terms
        .iter()
        .any(|term| mentions(graph, term.coefficient, symbols))
    {
        return None;
    }
    let is_taylor = terms.iter().all(|term| {
        term.exponents
            .iter()
            .all(|e| e.is_integer() && !e.is_negative())
    });
    let term = build(graph, bases, &terms)?;
    Some(Expansion { term, is_taylor })
}

/// Builds the truncated sum `sum coefficient * product base^exponent`.
///
/// Returns `None` when a `-oo` expansion would need a fractional parity
/// (which is not a real-variable expansion).
fn build(
    graph: &mut Graph,
    bases: &[Base],
    terms: &[OutputTerm],
) -> Option<NodeId> {
    let mut sum = Vec::with_capacity(terms.len());
    for term in terms {
        let mut coefficient = term.coefficient;
        let mut factors = Vec::new();
        for (index, base) in bases.iter().enumerate() {
            let Some(exponent) = term.exponents.get(index) else {
                continue;
            };
            if exponent.is_zero() {
                continue;
            }
            if base.signed {
                // `-oo`: odd powers carry a minus sign.
                let integer = integer_exponent(exponent)?;
                if integer.rem_euclid(2) == 1 {
                    coefficient = negate(graph, coefficient);
                }
                factors.push(power(graph, base.node, integer.saturating_neg()));
            } else {
                let printed = if base.reciprocal {
                    -exponent.clone()
                } else {
                    exponent.clone()
                };
                factors.push(power_fraction(graph, base.node, &printed));
            }
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
    Some(match sum.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &sum),
    })
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
    use super::has_essential_singularity;
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
    fn puiseux_monomial_factors() {
        assert_eq!(
            run("mseries(sqrt(x*y), list(x, y), list(0, 0), 4)"),
            "x^(1/2)*y^(1/2)"
        );
        assert_eq!(run("mseries(x^(4/3), list(x), list(0), 3)"), "x^(4/3)");
        assert_eq!(
            run("mseries((x*y)^(5/2), list(x, y), list(0, 0), 5)"),
            "x^(5/2)*y^(5/2)"
        );
        assert_eq!(
            run("mseries(x^(1/2)*exp(y), list(x, y), list(0, 0), 3)"),
            "1/2*x^(1/2)*y^2 + x^(1/2)*y + x^(1/2)"
        );
        assert_eq!(
            run("mseries(sqrt(x^2*(1 + y)), list(x, y), list(0, 0), 3)"),
            "1/2*x*y - 1/8*x*y^2 + x"
        );
        assert_eq!(
            run("mseries(x^(1/2)/(1 - x - y), list(x, y), list(0, 0), 2)"),
            "x^(3/2) + x^(1/2)*y + x^(1/2)"
        );
        assert_eq!(run("mseries(sqrt(x), list(x), list(0), 2)"), "x^(1/2)");
        // A fractional power of an analytic unit is a Taylor series and is
        // expanded by the ordinary recursion.
        assert_eq!(
            run("mseries((1 + x)^(1/2), list(x), list(0), 2)"),
            "1/2*x - 1/8*x^2 + 1"
        );
    }

    #[test]
    fn no_puiseux_expansion_at_a_branch_locus() {
        // The leading behaviour depends on the direction, so no Puiseux
        // series exists at the origin. The request must not produce one.
        for src in [
            "mseries((x^2 + y^2)^(1/2), list(x, y), list(0, 0), 4)",
            "mseries(sqrt(x^2 - y^2), list(x, y), list(0, 0), 3)",
            "mseries(sqrt(x*y + x^2), list(x, y), list(0, 0), 3)",
            "mseries((x^2 + y^2)^(1/3), list(x, y), list(0, 0), 3)",
        ] {
            let (text, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(!reduced, "`{src}` should stay unevaluated: {text}");
            assert!(
                text.contains("mseries("),
                "`{src}` must not be mis-expanded: {text}"
            );
        }
    }

    /// Whether the operator's classifier recognises `src` as an essential
    /// singularity in the given variables.
    fn detects_essential(
        src: &str,
        variables: &[&str],
    ) -> bool {
        use crate::graph::Engine;
        use crate::graph::Env;
        use crate::graph::Graph;

        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &[calculus()]).expect("install");
        let body = graph.parse(src).expect("parse");
        let nodes: Vec<_> = variables
            .iter()
            .map(|name| {
                let symbol = graph.interner_mut().symbol(name);
                graph.symbol_node(symbol)
            })
            .collect();
        let env = Env::symbolic();
        let mut cx = crate::graph::Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        has_essential_singularity(&mut cx, body, &nodes)
    }

    #[test]
    fn essential_singularities_are_recognised() {
        for (src, variables) in [
            ("exp(1/u)", &["u"][..]),
            ("sin(1/u)", &["u"][..]),
            ("cos(1/u)", &["u"][..]),
            ("sinh(1/u^2)", &["u"][..]),
            ("exp(u1/u2)", &["u1", "u2"][..]),
        ] {
            assert!(
                detects_essential(src, variables),
                "`{src}` is an essential singularity"
            );
        }
        // No false positives: a pole of the whole function, an analytic
        // elementary function and a fractional (Puiseux) argument are fine.
        for (src, variables) in [
            ("exp(u)", &["u"][..]),
            ("sin(u1 + u2)", &["u1", "u2"][..]),
            ("exp(u1)/u2^2", &["u1", "u2"][..]),
            ("sin(u)/u", &["u"][..]),
            ("exp(sqrt(u))", &["u"][..]),
            ("exp(1/(1 + u))", &["u"][..]),
            ("ln(u)", &["u"][..]),
        ] {
            assert!(
                !detects_essential(src, variables),
                "`{src}` is not an essential singularity"
            );
        }
    }

    #[test]
    fn essential_singularities_stay_unevaluated() {
        for src in [
            "mseries(exp(1/x), list(x), list(0), 2)",
            "mseries(exp(1/x), list(x), list(0), 6)",
            "mseries(sin(1/x), list(x), list(0), 2)",
            "mseries(sin(1/y), list(x, y), list(0, 0), 2)",
            "mseries(exp(x/y), list(x, y), list(0, 0), 2)",
            "mseries(exp(1/x^2), list(x), list(0), 3)",
        ] {
            let (text, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(!reduced, "`{src}` should stay unevaluated: {text}");
            assert!(
                text.contains("mseries("),
                "`{src}` must not be mis-expanded: {text}"
            );
        }
    }

    #[test]
    fn unevaluated_requests_stay() {
        for src in [
            "mtaylor(1/(x*y), list(x, y), list(0, 0), 2)",
            // `mtaylor` is the Taylor special case: fractional powers are
            // not Taylor polynomials.
            "mtaylor(sqrt(x), list(x), list(0), 2)",
            "mtaylor(sqrt(x*y), list(x, y), list(0, 0), 2)",
            // A logarithmic singularity is not a Puiseux series.
            "mseries(ln(x), list(x), list(0), 2)",
        ] {
            let (text, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(!reduced, "`{src}` should stay unevaluated: {text}");
        }
    }

    #[test]
    fn constant_bodies_never_leak_internal_variables() {
        // A body equal to a constant is expanded to that constant. The
        // concrete term of the request still mentions the translated
        // variable (`(x + 1)/(x + 1)`), and the value-based dependence test
        // says it does not depend on it, so the coefficient must be
        // canonicalised; it must never be returned with the operator's
        // internal `u` symbols.
        for (src, want) in [
            ("mseries((x + 1)/(x + 1), list(x), list(0), 1)", "1"),
            ("mseries((x + 1)/(x + 1), list(x), list(1), 1)", "1"),
            ("mseries((x*y)/(x*y), list(x, y), list(0, 0), 2)", "1"),
            ("mseries((x*y)/(x*y), list(x, y), list(1, 2), 2)", "1"),
            ("mseries(x - x + 1, list(x), list(0), 1)", "1"),
            ("mseries(0*x + 1, list(x), list(0), 1)", "1"),
            ("mseries(exp(x)*exp(-x), list(x), list(0), 1)", "1"),
            ("mseries(sqrt(x)^2/x, list(x), list(0), 1)", "1"),
            ("mseries(sin(x)^2 + cos(x)^2, list(x), list(0), 2)", "1"),
            (
                "mseries((x*y)^(1/2)*(x*y)^(-1/2), list(x, y), list(0, 0), 1)",
                "1",
            ),
            ("mtaylor((x + 1)/(x + 1), list(x), list(0), 1)", "1"),
            // A constant multiple survives as itself.
            ("mseries(exp(a)*x/x, list(x), list(0), 1)", "exp(a)"),
            // A non-constant body that cancels: the expansion is the
            // simplified function, without any internal variable.
            ("mseries((x^2 - 1)/(x - 1), list(x), list(0), 2)", "x + 1"),
            ("mseries((x*y + y^2)/y, list(x, y), list(0, 0), 2)", "x + y"),
        ] {
            let text = run(src);
            assert!(
                !text.contains("u#"),
                "`{src}` leaked an internal variable: {text}"
            );
            assert_eq!(text, want, "`{src}`");
        }
    }

    #[test]
    fn further_puiseux_monomial_factors() {
        for (src, want) in [
            // Nested roots need a higher ramification index.
            (
                "mseries((x*y)^(1/4), list(x, y), list(0, 0), 2)",
                "x^(1/4)*y^(1/4)",
            ),
            (
                "mseries((x^(1/2)*y)^(1/2), list(x, y), list(0, 0), 2)",
                "x^(1/4)*y^(1/2)",
            ),
            // A square times a unit: the square root is `x` times the unit
            // series, even though the unit depends on another variable.
            (
                "mseries(sqrt(x^2*y), list(x, y), list(0, 0), 3)",
                "x*y^(1/2)",
            ),
            (
                "mseries(sqrt(x^2 + x^3), list(x, y), list(0, 0), 3)",
                "1/2*x^2 - 1/8*x^3 + x",
            ),
            (
                "mseries(sqrt(x^2*(1 + y^2)), list(x, y), list(0, 0), 3)",
                "1/2*x*y^2 + x",
            ),
            // Negative fractional exponents.
            (
                "mseries((x*y)^(-1/3), list(x, y), list(0, 0), 2)",
                "1/(x^(1/3)*y^(1/3))",
            ),
            // A fractional power with non-negative total degree above the
            // request truncates to zero.
            ("mseries((x^2*y)^(3/2), list(x, y), list(0, 0), 3)", "0"),
        ] {
            assert_eq!(run(src), want, "`{src}`");
        }
    }

    #[test]
    fn branch_loci_are_refused_but_monomial_factors_are_not() {
        // No Puiseux series at the origin: the leading behaviour depends
        // on the direction.
        for src in [
            "mseries(sqrt(x^2*y + x^3), list(x, y), list(0, 0), 3)",
            "mseries(sqrt(x^2 + y^4), list(x, y), list(0, 0), 3)",
            "mseries(sqrt(x^4 + y^4), list(x, y), list(0, 0), 3)",
        ] {
            let (text, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(!reduced, "`{src}` should stay unevaluated: {text}");
            assert!(
                text.contains("mseries("),
                "`{src}` must not be mis-expanded: {text}"
            );
        }
        // A monomial factor times a unit is a Puiseux series; it must not be
        // refused by the same checks.
        for src in [
            "mseries(sqrt(x^2*y), list(x, y), list(0, 0), 2)",
            "mseries(sqrt(x^2 + x^3), list(x, y), list(0, 0), 2)",
            "mseries(sqrt(x^2*(1 + y^2)), list(x, y), list(0, 0), 2)",
        ] {
            let (_, reduced) = reduce_with(&[calculus()], src, &[]);
            assert!(reduced, "`{src}` should be expanded");
        }
    }

    #[test]
    fn essential_patterns_expand_at_regular_points() {
        // The same elementary functions whose argument has a pole at the
        // expansion point are analytic at a regular point and must expand
        // there; the essential-singularity refusal must not fire.
        assert_eq!(
            run("mseries(exp(1/x), list(x), list(1), 2)"),
            "3/2*E*(x - 1)^2 - E*(x - 1) + E"
        );
        assert_eq!(
            run("mseries(exp(1/(1 + x)), list(x), list(0), 3)"),
            "3/2*x^2*E - 13/6*x^3*E - x*E + E"
        );
        assert_eq!(
            run("mseries(sin(1/x), list(x), list(1), 3)"),
            "(x - 1)^3*(sin(1) - 5/6*cos(1)) + (x - 1)^2*(cos(1) - 1/2*sin(1)) - (x - 1)*cos(1) + sin(1)"
        );
    }
}
