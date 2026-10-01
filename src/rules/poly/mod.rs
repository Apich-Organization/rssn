//! Polynomials: expansion, collection, division, greatest common divisors
//! and factorisation.
//!
//! The operators here are *form requests*. `expand(e)`, `factor(e)` and
//! `collect(e, x)` are all equal to `e`; what they ask for is a particular
//! way of writing it. Their kernels therefore return
//! [`Outcome::Pinned`]: the class of the request is pinned to the computed
//! form, so that extraction returns it verbatim instead of whatever
//! spelling happens to be cheapest.

pub mod apart;
pub mod groebner;
pub mod repr;
pub mod univariate;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Extractor;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SizeCost;
use crate::graph::Tier;

use self::groebner::GPoly;
use self::groebner::GroebnerLimits;
use self::groebner::Order;
use self::repr::from_term;
use self::repr::to_term;
use self::repr::Gens;
use self::repr::Limits;
use self::repr::Poly;
use self::univariate::QPoly;
use super::arith::arith;

/// The polynomial rule set.
#[must_use]
pub fn poly() -> RuleSet {
    RuleSet::new("poly", install).needs(arith())
}

/// Which polynomial request a [`PolyKernel`] serves.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Expand,
    Collect,
    Factor,
    Degree,
    Coeff,
    Quo,
    Rem,
    Gcd,
    Together,
    Cancel,
    Apart,
    Groebner,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    // A form request has the value of its first argument.
    let identity = |a: &[f64]| a.first().copied().unwrap_or(f64::NAN);
    let requests = [
        ("expand", 1, Request::Expand, true),
        ("collect", 2, Request::Collect, true),
        ("factor", 1, Request::Factor, true),
        ("together", 1, Request::Together, true),
        ("cancel", 1, Request::Cancel, true),
        ("apart", 2, Request::Apart, true),
        ("degree", 2, Request::Degree, false),
        ("coeff", 3, Request::Coeff, false),
        ("quo", 3, Request::Quo, false),
        ("rem", 3, Request::Rem, false),
        ("pgcd", 3, Request::Gcd, false),
    ];
    for (name, arity, request, is_form) in requests {
        let mut desc = OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
        if is_form {
            desc = desc.eval(identity);
        }
        let op = i.op(desc)?;
        i.kernel(&format!("poly/{name}"), Tier::Reduce, PolyKernel { op, request });
    }
    // groebner(list(p1, ...), list(x1, ...)) or with a third argument, one
    // of the symbols `lex`, `grlex`, `grevlex`.
    let op = i.op(OpDescriptor::new("groebner", Arity::Variadic).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("poly/groebner", Tier::Reduce, PolyKernel { op, request: Request::Groebner });
    i.kernel("poly/collapse", Tier::Normalize, Collapse);
    Ok(())
}

struct PolyKernel {
    op: OpId,
    request: Request,
}

/// Expands a sum whose terms hide a cancellation behind a common factor —
/// `I*h*(x*f' + f) - I*h*x*f'` — and keeps the expansion when it is
/// smaller. Distribution is not a rewrite rule because it usually makes
/// terms larger; here it is tried once per sum and only kept when it pays.
struct Collapse;

/// Number of nodes of a term, counted as a tree, up to `cap`.
fn tree_size(
    graph: &Graph,
    node: NodeId,
    cap: usize,
) -> usize {
    let mut count = 0_usize;
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        count += 1;
        if count > cap {
            break;
        }
        stack.extend_from_slice(graph.children(n));
    }
    count
}

impl Kernel for Collapse {
    fn ops(&self) -> Vec<OpId> {
        vec![core::ADD, core::MUL]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        const CAP: usize = 240;
        if cx.graph.op(node) == core::MUL {
            return Self::quotient(cx, node);
        }
        let graph = &mut *cx.graph;
        // This sum with each term in its closed form: the class as a whole
        // may still be best spelled by a request that produced it.
        let mut terms = Vec::new();
        for &child in graph.children(node).to_vec().iter() {
            let closed = Extractor::new(graph, &[child], &crate::graph::ClosedForm).build(graph, child);
            let Some(t) = closed.or_else(|| best(graph, child)) else {
                return Outcome::Pass;
            };
            terms.push(t);
        }
        let term = graph.node(core::ADD, &terms);
        if graph.op(term) != core::ADD {
            return Outcome::Pass;
        }
        // Only sums with a product of a sum somewhere in them can gain.
        let hides_a_sum = graph.children(term).iter().any(|&t| {
            graph.op(t) == core::MUL && graph.children(t).iter().any(|&f| graph.op(f) == core::ADD)
        });
        let before = tree_size(graph, term, CAP);
        if !hides_a_sum || before > CAP {
            return Outcome::Pass;
        }
        let Some(expanded) = Self::canonical(graph, term) else {
            return Outcome::Pass;
        };
        Self::accept(cx, expanded, before)
    }
}

impl Collapse {
    /// Keeps `candidate` when it is smaller than `before` nodes, or — in a
    /// top-level run only, to keep nested simplifications from recursing
    /// into each other — when the other rules make it smaller
    /// (`exp(a)·exp(-a) = 1` after multiplying out).
    fn accept(
        cx: &mut Cx<'_>,
        candidate: NodeId,
        before: usize,
    ) -> Outcome {
        const CAP: usize = 240;
        if tree_size(cx.graph, candidate, CAP) < before {
            return Outcome::Equal(candidate);
        }
        if cx.env.depth > 0 {
            return Outcome::Pass;
        }
        let simplified = cx.simplify(candidate);
        if tree_size(cx.graph, simplified, CAP) < before { Outcome::Equal(simplified) } else { Outcome::Pass }
    }

    /// The rational normal form of `term`: one numerator over one
    /// denominator, common monomials and exact polynomial factors
    /// cancelled.
    fn canonical(
        graph: &mut Graph,
        term: NodeId,
    ) -> Option<NodeId> {
        let limits = Limits { terms: 64, exponent: 8 };
        let mut gens = Gens::default();
        let r = ratio(graph, &mut gens, term, limits)?;
        if r.denom.is_zero() {
            return None;
        }
        if r.numer.is_zero() {
            return Some(graph.int(0));
        }
        let (mut numer, mut denom) = (r.numer, r.denom);
        // Common monomial factor.
        let (cn, cd) = (numer.monomial_content(), denom.monomial_content());
        let common: crate::rules::poly::repr::Mono = cn
            .iter()
            .filter_map(|&(g, e)| cd.iter().find(|&&(h, _)| h == g).map(|&(_, f)| (g, e.min(f))))
            .collect();
        if !common.is_empty() {
            numer = numer.div_monomial(&common);
            denom = denom.div_monomial(&common);
        }
        // One side dividing the other.
        if denom.len() > 1 || numer.len() > 1 {
            if let Some(q) = numer.div_exact(&denom, 64) {
                numer = q;
                denom = Poly::constant(Number::from(1));
            } else if let Some(q) = denom.div_exact(&numer, 64) {
                denom = q;
                numer = Poly::constant(Number::from(1));
            }
        }
        let top = to_term(graph, &gens, &numer);
        if let Some(c) = denom.as_constant() {
            let inverse = graph.num(c.recip()?);
            return Some(if c.is_one() { top } else { graph.node(core::MUL, &[inverse, top]) });
        }
        let bottom = to_term(graph, &gens, &denom);
        let minus_one = graph.int(-1);
        let inverse = graph.node(core::POW, &[bottom, minus_one]);
        Some(if graph.number_of(top).is_some_and(Number::is_one) { inverse } else { graph.node(core::MUL, &[top, inverse]) })
    }

    /// A product with a sum in a denominator: cancel what can be
    /// cancelled, keep the result when it is smaller.
    fn quotient(
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        const CAP: usize = 160;
        let graph = &mut *cx.graph;
        let Some(term) = best(graph, node) else {
            return Outcome::Pass;
        };
        if graph.op(term) != core::MUL {
            return Outcome::Pass;
        }
        let divides_by_sum = graph.children(term).iter().any(|&f| {
            graph.op(f) == core::POW
                && graph.children(f).first().is_some_and(|&b| graph.op(b) == core::ADD)
                && graph.children(f).get(1).and_then(|&e| graph.number_of(e)).is_some_and(Number::is_negative)
        });
        let has_sum = graph.children(term).iter().any(|&f| graph.op(f) == core::ADD);
        let before = tree_size(graph, term, CAP);
        if !(divides_by_sum || has_sum) || before > CAP {
            return Outcome::Pass;
        }
        let Some(result) = Self::canonical(graph, term) else {
            return Outcome::Pass;
        };
        Self::accept(cx, result, before)
    }
}

/// The best known concrete form of `node`.
pub(crate) fn best(
    graph: &mut Graph,
    node: NodeId,
) -> Option<NodeId> {
    Extractor::new(graph, &[node], &SizeCost).build(graph, node)
}

impl Kernel for PolyKernel {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let args = graph.children(node).to_vec();
        let result = match (self.request, args.as_slice()) {
            | (Request::Expand, &[e]) => expand(graph, e).map(Outcome::Pinned),
            | (Request::Collect, &[e, x]) => collect(graph, e, x).map(Outcome::Pinned),
            | (Request::Factor, &[e]) => factor(graph, e).map(Outcome::Pinned),
            | (Request::Together, &[e]) => together(graph, e, false).map(Outcome::Pinned),
            | (Request::Cancel, &[e]) => together(graph, e, true).map(Outcome::Pinned),
            | (Request::Apart, &[e, x]) => apart_term(graph, e, x).map(Outcome::Pinned),
            | (Request::Degree, &[p, x]) => degree(graph, p, x).map(Outcome::Equal),
            | (Request::Coeff, &[p, x, n]) => coeff(graph, p, x, n).map(Outcome::Equal),
            | (Request::Quo, &[p, q, x]) => divide(graph, p, q, x).map(|(quo, _)| Outcome::Equal(quo)),
            | (Request::Rem, &[p, q, x]) => divide(graph, p, q, x).map(|(_, rem)| Outcome::Equal(rem)),
            | (Request::Gcd, &[p, q, x]) => gcd(graph, p, q, x).map(Outcome::Equal),
            | (Request::Groebner, &[polys, vars]) => {
                groebner_basis(graph, polys, vars, Order::Lex).map(Outcome::Pinned)
            },
            | (Request::Groebner, &[polys, vars, order]) => graph
                .symbol_of(order)
                .and_then(|s| match graph.interner().symbol_name(s) {
                    | "lex" => Some(Order::Lex),
                    | "grlex" => Some(Order::GradedLex),
                    | "grevlex" => Some(Order::GradedRevLex),
                    | _ => None,
                })
                .and_then(|order| groebner_basis(graph, polys, vars, order))
                .map(Outcome::Pinned),
            | _ => None,
        };
        result.unwrap_or(Outcome::Pass)
    }

    fn revisit(&self) -> bool {
        // The argument may itself contain a request that is reduced later.
        true
    }
}

/// `e` as a quotient of univariate polynomials in `x` with rational
/// coefficients (lowest degree first), or `None` when it is not a rational
/// function of `x` alone.
pub(crate) fn rational_function_in(
    graph: &mut Graph,
    e: NodeId,
    x: NodeId,
) -> Option<(QPoly, QPoly)> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let r = ratio(graph, &mut gens, term, Limits::default())?;
    if gens.len() != 1 {
        return None;
    }
    let as_q = |p: &Poly| -> Option<QPoly> { p.univariate_in(gx)?.iter().map(Number::to_rational).collect() };
    Some((as_q(&r.numer)?, as_q(&r.denom)?))
}

/// Partial fractions of a rational function of `x` over `Q`:
/// `polynomial + sum numerator / factor^power`.
fn apart_term(
    graph: &mut Graph,
    e: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let r = ratio(graph, &mut gens, term, Limits::default())?;
    if gens.len() != 1 {
        return None;
    }
    let as_q = |p: &Poly| -> Option<QPoly> { p.univariate_in(gx)?.iter().map(Number::to_rational).collect() };
    let parts = apart::apart(&as_q(&r.numer)?, &as_q(&r.denom)?)?;
    let poly_term = |graph: &mut Graph, q: &[BigRational]| {
        let numbers: Vec<Number> = q.iter().cloned().map(Number::rat).collect();
        to_term(graph, &gens, &Poly::from_univariate(gx, &numbers))
    };
    let mut terms = Vec::new();
    if !parts.quotient.is_empty() {
        terms.push(poly_term(graph, &parts.quotient));
    }
    for piece in &parts.pieces {
        let numerator = poly_term(graph, &piece.numerator);
        let base = poly_term(graph, &piece.factor);
        let exponent = graph.int(-i64::from(piece.power));
        let denominator = graph.node(core::POW, &[base, exponent]);
        terms.push(graph.node(core::MUL, &[numerator, denominator]));
    }
    Some(match terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &terms),
    })
}

fn expand(
    graph: &mut Graph,
    e: NodeId,
) -> Option<NodeId> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let poly = from_term(graph, &mut gens, term, Limits::default())?;
    Some(to_term(graph, &gens, &poly))
}

/// A polynomial in the main variable `x` together with the generator table
/// it was built against.
struct InVariable {
    gens: Gens,
    poly: Poly,
    x: u32,
}

fn in_variable(
    graph: &mut Graph,
    e: NodeId,
    x: NodeId,
) -> Option<InVariable> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let x = gens.index(graph, x);
    let poly = from_term(graph, &mut gens, term, Limits::default())?;
    Some(InVariable { gens, poly, x })
}

fn power(
    graph: &mut Graph,
    base: NodeId,
    exponent: usize,
) -> Option<NodeId> {
    match exponent {
        | 0 => None,
        | 1 => Some(base),
        | n => {
            let e = graph.int(i64::try_from(n).ok()?);
            Some(graph.node(core::POW, &[base, e]))
        },
    }
}

fn collect(
    graph: &mut Graph,
    e: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let view = in_variable(graph, e, x)?;
    let mut terms = Vec::new();
    for (k, coefficient) in view.poly.coefficients_in(view.x).iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let c = to_term(graph, &view.gens, coefficient);
        let term = match power(graph, x, k) {
            | None => c,
            | Some(xk) if graph.number_of(c).is_some_and(Number::is_one) => xk,
            | Some(xk) => graph.node(core::MUL, &[c, xk]),
        };
        terms.push(term);
    }
    Some(match terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &terms),
    })
}

fn degree(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let view = in_variable(graph, p, x)?;
    // The coefficients must not hide the variable inside a generator.
    independent(graph, &view)?;
    Some(graph.int(i64::from(view.poly.degree_in(view.x))))
}

/// Checks that `x` occurs in `view` only as the generator itself, so that
/// degree and coefficients in `x` are meaningful.
fn independent(
    graph: &Graph,
    view: &InVariable,
) -> Option<()> {
    let x_node = view.gens.node(view.x)?;
    let symbol = graph.symbol_of(x_node)?;
    for generator in view.poly.support() {
        if generator == view.x {
            continue;
        }
        let node = view.gens.node(generator)?;
        if graph.depends_on(graph.find(node), symbol) {
            return None;
        }
    }
    Some(())
}

fn coeff(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
    n: NodeId,
) -> Option<NodeId> {
    let power = usize::try_from(graph.number_of(n)?.to_i64()?).ok()?;
    let view = in_variable(graph, p, x)?;
    independent(graph, &view)?;
    let coefficients = view.poly.coefficients_in(view.x);
    Some(match coefficients.get(power) {
        | Some(c) => to_term(graph, &view.gens, c),
        | None => graph.int(0),
    })
}

/// The exact rational coefficients of `e` as a univariate polynomial in
/// `x`; `None` if `e` involves anything else or inexact numbers.
fn rational_in(
    graph: &mut Graph,
    e: NodeId,
    x: NodeId,
) -> Option<QPoly> {
    let view = in_variable(graph, e, x)?;
    let mut coefficients: QPoly =
        view.poly.univariate_in(view.x)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
    while coefficients.last().is_some_and(num_traits::Zero::is_zero) {
        coefficients.pop();
    }
    Some(coefficients)
}

fn from_rational(
    graph: &mut Graph,
    coefficients: &[BigRational],
    x: NodeId,
) -> NodeId {
    let mut gens = Gens::default();
    let index = gens.index(graph, x);
    let numbers: Vec<Number> = coefficients.iter().cloned().map(Number::rat).collect();
    to_term(graph, &gens, &Poly::from_univariate(index, &numbers))
}

fn divide(
    graph: &mut Graph,
    p: NodeId,
    q: NodeId,
    x: NodeId,
) -> Option<(NodeId, NodeId)> {
    let (a, b) = (rational_in(graph, p, x)?, rational_in(graph, q, x)?);
    let (quo, rem) = univariate::divrem(&a, &b)?;
    Some((from_rational(graph, &quo, x), from_rational(graph, &rem, x)))
}

fn gcd(
    graph: &mut Graph,
    p: NodeId,
    q: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let (a, b) = (rational_in(graph, p, x)?, rational_in(graph, q, x)?);
    Some(from_rational(graph, &univariate::gcd(&a, &b), x))
}

/// Builds `unit * prod factor_i^e_i` from integer factors in `x`.
fn factored_term(
    graph: &mut Graph,
    unit: &BigRational,
    factors: &[(Vec<BigInt>, u32)],
    x: NodeId,
) -> NodeId {
    let mut pieces = Vec::with_capacity(factors.len() + 1);
    if !unit.is_one() || factors.is_empty() {
        pieces.push(graph.num(Number::rat(unit.clone())));
    }
    for (factor, multiplicity) in factors {
        let coefficients: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
        let base = from_rational(graph, &coefficients, x);
        pieces.push(power(graph, base, *multiplicity as usize).unwrap_or(base));
    }
    match pieces.as_slice() {
        | [only] => *only,
        | _ => graph.node(core::MUL, &pieces),
    }
}

fn factor(
    graph: &mut Graph,
    e: NodeId,
) -> Option<NodeId> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let poly = from_term(graph, &mut gens, term, Limits::default())?;
    if !poly.is_exact() {
        return None;
    }
    let support = poly.support();
    if let &[only] = support.as_slice() {
        let coefficients: QPoly =
            poly.univariate_in(only)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
        let (unit, factors) = univariate::factor(&coefficients);
        return Some(factored_term(graph, &unit, &factors, gens.node(only)?));
    }
    if support.is_empty() {
        return Some(to_term(graph, &gens, &poly));
    }
    // Several generators: pull out the numeric content and the monomial
    // common to all terms; the remaining factor stays expanded.
    let mut content = BigRational::from_integer(BigInt::from(0));
    let mut common: Option<Vec<(u32, u32)>> = None;
    for (mono, c) in poly.terms() {
        let c = c.to_rational()?;
        content = rational_gcd(&content, &c);
        common = Some(match common {
            | None => mono.clone(),
            | Some(old) => old
                .iter()
                .filter_map(|&(g, e)| mono.iter().find(|&&(h, _)| h == g).map(|&(_, f)| (g, e.min(f))))
                .collect(),
        });
    }
    let common = common.unwrap_or_default();
    if content.is_one() && common.is_empty() {
        return Some(to_term(graph, &gens, &poly));
    }
    let mut reduced = Poly::zero();
    for (mono, c) in poly.terms() {
        let lowered: Vec<(u32, u32)> = mono
            .iter()
            .filter_map(|&(g, e)| {
                let shared = common.iter().find(|&&(h, _)| h == g).map_or(0, |&(_, f)| f);
                (e > shared).then_some((g, e - shared))
            })
            .collect();
        let scaled = Number::rat(c.to_rational()? / &content);
        reduced = reduced.add(&Poly::monomial(lowered, scaled));
    }
    let outer = to_term(graph, &gens, &Poly::monomial(common, Number::rat(content)));
    let inner = to_term(graph, &gens, &reduced);
    Some(graph.node(core::MUL, &[outer, inner]))
}

/// Converts the polynomials of `list(p1, ...)` to the Gröbner
/// representation over the variables of `list(x1, ...)`. `None` when a
/// polynomial involves anything but those variables, or inexact numbers.
pub(crate) fn to_groebner(
    graph: &mut Graph,
    polys: NodeId,
    vars: NodeId,
    order: Order,
) -> Option<(Vec<GPoly>, Gens)> {
    if graph.op(polys) != core::LIST || graph.op(vars) != core::LIST {
        return None;
    }
    let variables = graph.children(vars).to_vec();
    let mut gens = Gens::default();
    for &v in &variables {
        graph.symbol_of(v)?;
        gens.index(graph, v);
    }
    if gens.len() != variables.len() {
        return None;
    }
    let mut out = Vec::new();
    for p in graph.children(polys).to_vec() {
        // An equation `l = r` stands for `l - r`.
        let term = match *graph.children(p) {
            | [lhs, rhs] if graph.op(p) == core::EQ => {
                let minus_one = graph.int(-1);
                let negated = graph.node(core::MUL, &[minus_one, rhs]);
                graph.node(core::ADD, &[lhs, negated])
            },
            | _ => p,
        };
        let term = best(graph, term)?;
        let poly = from_term(graph, &mut gens, term, Limits::default())?;
        if gens.len() != variables.len() {
            return None;
        }
        let mut terms = Vec::with_capacity(poly.len());
        for (mono, coeff) in poly.terms() {
            let mut exponents = vec![0_u32; variables.len()];
            for &(generator, e) in mono {
                *exponents.get_mut(generator as usize)? = e;
            }
            terms.push((exponents, coeff.to_rational()?));
        }
        out.push(GPoly::new(terms, order));
    }
    Some((out, gens))
}

/// Converts a Gröbner-representation polynomial back to a term.
pub(crate) fn from_groebner(
    graph: &mut Graph,
    gens: &Gens,
    poly: &GPoly,
) -> NodeId {
    let mut out = Poly::zero();
    for (exponents, coeff) in poly.terms() {
        let mono: Vec<(u32, u32)> = exponents
            .iter()
            .enumerate()
            .filter(|&(_, &e)| e > 0)
            .filter_map(|(g, &e)| u32::try_from(g).ok().map(|g| (g, e)))
            .collect();
        out = out.add(&Poly::monomial(mono, Number::rat(coeff.clone())));
    }
    to_term(graph, gens, &out)
}

fn groebner_basis(
    graph: &mut Graph,
    polys: NodeId,
    vars: NodeId,
    order: Order,
) -> Option<NodeId> {
    let (generators, gens) = to_groebner(graph, polys, vars, order)?;
    let basis = groebner::groebner(&generators, order, GroebnerLimits::default())?;
    let terms: Vec<NodeId> = basis.iter().map(|g| from_groebner(graph, &gens, g)).collect();
    Some(graph.node(core::LIST, &terms))
}

/// Greatest common divisor of two rationals: the largest rational dividing
/// both with integer quotients. Positive; `gcd(0, x) = |x|`.
fn rational_gcd(
    a: &BigRational,
    b: &BigRational,
) -> BigRational {
    use num_integer::Integer;
    use num_traits::Signed;
    use num_traits::Zero;
    if a.is_zero() {
        return b.abs();
    }
    if b.is_zero() {
        return a.abs();
    }
    BigRational::new(a.numer().gcd(b.numer()), a.denom().lcm(b.denom()))
}

/// A quotient of two polynomials over a shared generator table.
pub(crate) struct Ratio {
    /// Numerator.
    pub(crate) numer: Poly,
    /// Denominator; never the zero polynomial for a well-formed term.
    pub(crate) denom: Poly,
}

/// Rewrites the term `node` as a single quotient of polynomials, treating
/// negative integer powers as denominators.
pub(crate) fn ratio(
    graph: &mut Graph,
    gens: &mut Gens,
    node: NodeId,
    limits: Limits,
) -> Option<Ratio> {
    let one = || Poly::constant(Number::from(1));
    let children = graph.children(node).to_vec();
    if graph.op(node) == core::ADD {
        let mut acc = Ratio { numer: Poly::zero(), denom: one() };
        for child in children {
            let r = ratio(graph, gens, child, limits)?;
            acc = if acc.denom == r.denom {
                Ratio { numer: acc.numer.add(&r.numer), denom: acc.denom }
            } else {
                Ratio {
                    numer: acc.numer.mul(&r.denom, limits.terms)?.add(&r.numer.mul(&acc.denom, limits.terms)?),
                    denom: acc.denom.mul(&r.denom, limits.terms)?,
                }
            };
        }
        return Some(acc);
    }
    if graph.op(node) == core::MUL {
        let mut acc = Ratio { numer: one(), denom: one() };
        for child in children {
            let r = ratio(graph, gens, child, limits)?;
            acc = Ratio {
                numer: acc.numer.mul(&r.numer, limits.terms)?,
                denom: acc.denom.mul(&r.denom, limits.terms)?,
            };
        }
        return Some(acc);
    }
    if let (true, &[base, exp]) = (graph.op(node) == core::POW, children.as_slice()) {
        if let Some(e) = graph.number_of(exp).and_then(Number::to_i64) {
            let magnitude = u32::try_from(e.unsigned_abs()).ok().filter(|&m| m <= limits.exponent)?;
            let r = ratio(graph, gens, base, limits)?;
            let (numer, denom) = (r.numer.pow(magnitude, limits.terms)?, r.denom.pow(magnitude, limits.terms)?);
            return Some(if e < 0 { Ratio { numer: denom, denom: numer } } else { Ratio { numer, denom } });
        }
    }
    Some(Ratio { numer: from_term_shared(graph, gens, node, limits)?, denom: one() })
}

/// [`from_term`] against an existing generator table.
fn from_term_shared(
    graph: &mut Graph,
    gens: &mut Gens,
    node: NodeId,
    limits: Limits,
) -> Option<Poly> {
    from_term(graph, gens, node, limits)
}

fn together(
    graph: &mut Graph,
    e: NodeId,
    cancel: bool,
) -> Option<NodeId> {
    let term = best(graph, e)?;
    let mut gens = Gens::default();
    let mut r = ratio(graph, &mut gens, term, Limits::default())?;
    if r.denom.is_zero() {
        return None;
    }
    if cancel {
        // Cancel a common factor when the quotient is univariate over Q.
        let mut support = r.numer.support();
        support.extend(r.denom.support());
        support.sort_unstable();
        support.dedup();
        if let &[only] = support.as_slice() {
            let as_q = |p: &Poly| -> Option<QPoly> {
                p.univariate_in(only)?.iter().map(Number::to_rational).collect()
            };
            if let (Some(n), Some(d)) = (as_q(&r.numer), as_q(&r.denom)) {
                let g = univariate::gcd(&n, &d);
                let lead = d.last().cloned().unwrap_or_else(BigRational::one);
                // Normalise so that the denominator is monic.
                let scaled: QPoly = g.iter().map(|c| c * &lead).collect();
                if let (Some((nq, _)), Some((dq, _))) = (univariate::divrem(&n, &scaled), univariate::divrem(&d, &scaled)) {
                    let back = |q: &QPoly| {
                        let numbers: Vec<Number> = q.iter().cloned().map(Number::rat).collect();
                        Poly::from_univariate(only, &numbers)
                    };
                    r = Ratio { numer: back(&nq), denom: back(&dq) };
                }
            }
        }
    }
    let numer = to_term(graph, &gens, &r.numer);
    if let Some(c) = r.denom.as_constant() {
        let inverse = graph.num(c.recip()?);
        return Some(if c.is_one() { numer } else { graph.node(core::MUL, &[inverse, numer]) });
    }
    let denom = to_term(graph, &gens, &r.denom);
    let minus_one = graph.int(-1);
    let inverse = graph.node(core::POW, &[denom, minus_one]);
    Some(if graph.number_of(numer).is_some_and(Number::is_one) {
        inverse
    } else {
        graph.node(core::MUL, &[numer, inverse])
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[poly()], src)
    }

    #[test]
    fn expand_is_pinned_against_cheaper_spellings() {
        assert_eq!(run("expand((x + 1)^2)"), "x^2 + 2*x + 1");
        assert_eq!(run("expand((x + y)*(x - y))"), "x^2 - y^2");
        assert_eq!(run("expand((a + b)^3)"), "a^3 + 3*a^2*b + 3*a*b^2 + b^3");
        // Nested inside something else, the expanded form survives too.
        assert_eq!(run("f(expand((x + 1)^2))"), "f(x^2 + 2*x + 1)");
        assert_eq!(run("expand(x*(x + 1)) - x^2"), "x");
    }

    #[test]
    fn collect_groups_by_powers() {
        assert_eq!(run("collect(a*x^2 + b*x + c*x^2 + d, x)"), "x^2*(a + c) + b*x + d");
        assert_eq!(run("collect((x + y)^2, x)"), "x^2 + 2*x*y + y^2");
        assert_eq!(run("collect(x*y + x*z + x, x)"), "x*(y + z + 1)");
    }

    #[test]
    fn factor_univariate() {
        assert_eq!(run("factor(x^2 - 1)"), "(x - 1)*(x + 1)");
        assert_eq!(run("factor(x^2 + 2*x + 1)"), "(x + 1)^2");
        assert_eq!(run("factor(6*x^2 + 5*x + 1)"), "(2*x + 1)*(3*x + 1)");
        assert_eq!(run("factor(x^4 - 1)"), "(x - 1)*(x + 1)*(x^2 + 1)");
        assert_eq!(run("factor(x^2 + 1)"), "x^2 + 1");
        assert_eq!(run("factor(2*x^2 - 2)"), "2*(x - 1)*(x + 1)");
        assert_eq!(run("factor(x^3/2 - x/2)"), "1/2*x*(x - 1)*(x + 1)");
        assert_eq!(
            simplify(&[poly(), crate::rules::elementary()], "factor(exp(t)^2 - 1)"),
            "(exp(t) - 1)*(exp(t) + 1)",
            "any generator will do"
        );
    }

    #[test]
    fn factor_multivariate_pulls_out_content() {
        assert_eq!(run("factor(2*x*y + 4*x*z)"), "2*x*(y + 2*z)");
        assert_eq!(run("factor(x*y + z)"), "x*y + z");
    }

    #[test]
    fn factoring_preserves_value() {
        for src in ["x^6 - 1", "x^5 - 3*x^4 + x - 3", "12*x^3 - 8*x^2 - 3*x + 2", "x^8 + x^4 + 1"] {
            let factored = run(&format!("factor({src})"));
            for at in [-1.7, 0.3, 2.2] {
                let want = eval(&[poly()], src, &[("x", at)]);
                let got = eval(&[poly()], &factored, &[("x", at)]);
                assert!((want - got).abs() < 1e-9 * (1.0 + want.abs()), "{src} = {factored}? at {at}: {want} vs {got}");
            }
            assert!(factored.contains('('), "{src} should split: {factored}");
        }
    }

    #[test]
    fn degree_and_coefficients() {
        assert_eq!(run("degree(3*x^4 + x, x)"), "4");
        assert_eq!(run("degree((x + 1)^3 * y, x)"), "3");
        assert_eq!(run("degree(7, x)"), "0");
        assert_eq!(run("coeff((x + y)^3, x, 2)"), "3*y");
        assert_eq!(run("coeff(x^2 + 1, x, 5)"), "0");
        // x hidden inside a generator: no polynomial degree.
        assert!(!reduce_with(&[poly()], "degree(f(x) + x, x)", &[]).1);
    }

    #[test]
    fn division_and_gcd() {
        assert_eq!(run("quo(x^3 - 1, x - 1, x)"), "x^2 + x + 1");
        assert_eq!(run("rem(x^3 - 1, x - 1, x)"), "0");
        assert_eq!(run("rem(x^2 + 1, x - 1, x)"), "2");
        assert_eq!(run("pgcd(x^2 - 1, x^2 - 2*x + 1, x)"), "x - 1");
        assert_eq!(run("pgcd(x^2 + 1, x + 1, x)"), "1");
        assert!(!reduce_with(&[poly()], "quo(x, 0, x)", &[]).1, "division by zero is not reduced");
    }

    #[test]
    fn groebner_bases() {
        assert_eq!(run("groebner(list(x + y - 3, x - y - 1), list(x, y))"), "list(y - 1, x - 2)");
        assert_eq!(run("groebner(list(x^2 + y^2 = 1, x = y), list(x, y))"), "list(y^2 - 1/2, x - y)");
        assert_eq!(run("groebner(list(x, x + 1), list(x))"), "list(1)");
        assert_eq!(
            run("groebner(list(x^3 - 2*x*y, x^2*y - 2*y^2 + x), list(x, y), grlex)"),
            "list(y^2 - 1/2*x, x*y, x^2)"
        );
        // A parameter that is not among the variables is refused.
        assert!(!reduce_with(&[poly()], "groebner(list(a*x + 1), list(x))", &[]).1);
    }

    #[test]
    fn rational_functions() {
        assert_eq!(run("together(1/x + 1/y)"), "(x + y)/(x*y)");
        assert_eq!(run("together(1/(x + 1) + 1/(x - 1))"), "2*x/(x^2 - 1)");
        assert_eq!(run("cancel((x^2 - 1)/(x - 1))"), "x + 1");
        assert_eq!(run("cancel((x^2 + 2*x + 1)/(x^2 - 1))"), "(x + 1)/(x - 1)");
        assert_eq!(run("cancel(x/x^3)"), "1/x^2");
        assert_eq!(run("apart(1/(x^2 - 1), x)"), "1/2/(x - 1) - 1/2/(x + 1)");
        assert_eq!(run("apart((x^3 + 1)/(x^2 + 1), x)"), "x + (1 - x)/(x^2 + 1)");
    }
}
