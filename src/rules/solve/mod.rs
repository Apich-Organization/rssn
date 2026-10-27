//! Equation solving.
//!
//! `solve(equation, x)` is a request whose value is the list of real
//! solutions for `x`; `solve(list(equations...), list(unknowns...))` gives a
//! list of solution tuples. `nsolve(equation, x, guess)` asks for one
//! numeric root near a starting point.
//!
//! The symbolic kernel reduces an equation to a polynomial in one
//! *generator* that carries the unknown — the unknown itself, or something
//! like `exp(x)` — solves that polynomial (exactly over the rationals by
//! factorisation, by formula when the coefficients are symbolic), and then
//! undoes the generator with the inverse functions declared by the
//! [`Inverse`] attribute. Systems go through Cramer's rule when they are
//! linear and through a lexicographic Gröbner basis otherwise.
//!
//! When the equation involves several functions of the unknown, the
//! sub-solvers of the heuristics module transform it (absolute values,
//! the Weierstrass substitution, Lambert W, `x^x`, sums of logarithms
//! combined into one, logarithmic substitutions, square roots isolated and
//! squared one at a time) and hand the result back to the solver;
//! polynomials of degree three and four get Cardano's and Ferrari's
//! formulas (also with symbolic coefficients); non-polynomial systems are
//! solved by elimination. Every candidate is checked against the original
//! equation.
//!
//! * **Quartics with parameters** by Ferrari/Descartes through the resolvent
//!   cubic `8 m^3 + 8 p m^2 + (2 p^2 - 8 r) m - q^2`, whose largest real
//!   root is positive (so `sqrt(2 m)` is defined); both its Cardano and
//!   its trigonometric expression are offered, plus the biquadratic case
//!   `q = 0`. As for cubics, a formula is real only on part of the
//!   parameter space; formulas that are right wherever they are defined are
//!   kept.
//! * **Quintics and higher degrees.** An irreducible rational quintic is
//!   classified by Dummit's sextic resolvent (exact integer arithmetic on
//!   roots computed to 190 digits): `F20`, `C5`/`D5` are solvable, `A5`/`S5`
//!   are not. A solvable quintic with one real root is written in radicals
//!   through Lagrange resolvents (the real root is `shift + (1/5λ) Σ R^(1/5)`
//!   over the roots of a quartic with coefficients in `Q(sqrt 5)`; so
//!   `x^5 - 5 x + 12` and `x^5 + 15 x + 12` come out in radicals); one with
//!   five real roots (the cyclic `x^5 + x^4 - 4 x^3 - 3 x^2 + 3 x + 1`, say)
//!   in cosine form: the resolvents are complex conjugate pairs
//!   `R = a ± i b`, the real quadratic factors of the quartic, and the roots
//!   are `(2/5λ) Σ_j |R_j|^(1/5) cos((arg R_j + 2π (m_j - j k))/5)`.
//!   de Moivre quintics `x^5 + 5 a x^3 + 5 a^2 x + b` with parameters are
//!   solved by `x = u - a/u` (real fifth roots, or the cosine form for five
//!   real roots).
//! * **Degree six and more.** An irreducible factor that is a functional
//!   composition `g(h(x))` (Kozen–Landau decomposition) is solved through
//!   the real roots of `g` and then `h(x) = u` for each; a self-reciprocal
//!   one of degree `2k` through `u = x + 1/x` (a polynomial of degree `k` in
//!   `u`). Every other real root of an irreducible factor of degree five or
//!   more (unsolvable quintics, higher degrees without such structure) is an
//!   exact `root_of(k, c0, ..., cn)` object: the `k`-th root of the
//!   polynomial, the real roots first (ascending), then the complex ones by
//!   real and imaginary part, canonical across factorisations, with
//!   rigorous numeric evaluation (real and complex); `solve` therefore
//!   always returns every real root. In numeric runs the values are
//!   returned instead. `solve_complex(eq, x)` gives *all* complex roots of a
//!   polynomial equation exactly (rationals, quadratic radicals with `I`,
//!   `root_of` objects).
//! * **Sextics and septics by their Galois group.** For an irreducible
//!   polynomial of degree six or seven that is not decomposable or
//!   reciprocal the Galois group is computed exactly (relative resolvents
//!   over all cosets of the transitive groups, with the roots to 190
//!   digits; see `rules::poly::galois_high`). A cyclic group whose real
//!   roots lie in a Gauss-period field of prime conductor `p` gives
//!   `Σ c_j 2 cos(2π j / m)` with rational `c_j` (the cyclic septic of
//!   conductor 29, `2 cos(2π/13)`'s sextic, …); the modulus `m` may be
//!   composite (conductor 91 = 7·13 for a sextic, 1247 = 29·43 for a septic),
//!   found from the real characters of `(Z/m)*` of order `n`. Every solvable sextic
//!   group is imprimitive, and its block tower (a cubic or quadratic for the
//!   block invariant `θ`, then a quadratic or cubic over `Q(θ)` with
//!   coefficients found numerically and recognised as rationals) is solved
//!   with the radical/trigonometric formulas, taking only the real values of
//!   `θ`. The metacyclic septic groups `D7`, `F21`, `F42` use the Lagrange
//!   resolvents `t_j = Σ ζ_7^{jk} z_k`: with one real root the real
//!   resolvents `R_j = t_j^7` are the roots of an integer sextic `Φ` (cyclic,
//!   solved recursively by the Gauss periods) and the real root is
//!   `(t_0 + Σ R_j^{1/7})/7`; with seven real roots (`F21`) the conjugate
//!   pairs `R_j = ρ_j e^{iφ_j}` give
//!   `z_k = (t_0 + 2 Σ ρ_j^{1/7} cos((φ_j + 2π m_j)/7 - 2π j k/7))/7` with `ρ_j²`,
//!   `cos φ_j` in a cyclic cubic field; Dickson polynomials `D7(x) - t` are
//!   the degenerate case with only `R_{±1}` non-zero. `solve_complex` gives
//!   exact forms (combinations of `exp(2πi j/m)`, or `ζ_7^{-jk}` times real
//!   seventh roots) for the non-real roots of the cyclic sextics/septics and
//!   of the one-real-root `D7`/`F42` septics. Not done: sextic towers with
//!   non-real block invariants, septics with three real pairs and no real
//!   root, polynomials of degree eight and more without decomposition: they
//!   stay `root_of` objects, as do the unsolvable groups (`A5`, `S5`, `A6`,
//!   `S6`, `PSL(3,2)`, `A7`, `S7`).
//!
//! Beyond that, the sub-solvers interlock:
//!
//! * **Parameters.** A polynomial with symbolic coefficients is factored
//!   over `Q` in the unknown and the parameters (`(x - a)(x - b)`,
//!   `x^3 + a x^2 + a x + 1`), polynomials in `x^k` go through `u = x^k`
//!   (`x^3 = a/b`, symbolic biquadratics), cubics get Cardano's formula
//!   (real cube roots) together with the trigonometric form, palindromic
//!   quartics `u = x + 1/x`; formulas are real only on part of the
//!   parameter space and are kept when they are right wherever defined.
//! * **Factors and substitutions.** When the numerator factors into several
//!   pieces that involve the unknown (`e^x (x - 1)`, `x^3 ln x`) each piece
//!   is solved; a common subterm `u = f(x)` (an argument of a function, a
//!   power `x^m`) turns `x^2 e^(x^2) = e` into `u e^u = e`; exponentials
//!   are made powers of one generator (`4^x` and `2^x`, `e^(x+1)`, `e^-x`),
//!   and the hyperbolic functions are written through `exp`.
//! * **Lambert W** for `α + β x + γ b^(dx)`, `x^m e^(dx)` and `x^m = c e^(dx)`
//!   (both real branches, also for symbolic arguments), `x ln x`, `x^x`.
//! * **Radicals** of any index by resultants (`∛(x+1) + ∛(x-1) = c`), sums of
//!   inverse circular functions through their tangent, and `A^B = C^D` by
//!   logarithms of both sides.
//! * **Trigonometric equations** whose arguments are `k (x + s)` are
//!   (any rational multipliers up to a multiple-angle index of 96, mixed shifts
//!   by the addition formulas) expanded into polynomials in `sin θ`, `cos θ` and solved with the
//!   polynomial in `sin θ` or `cos θ` alone or the Weierstrass substitution
//!   (covering `a sin x + b cos x = c`, multiple angles and sum-to-product
//!   forms); solutions that are rational multiples of `π` are written as
//!   such. `solve_general(eq, x)` returns *all* solutions as families with
//!   an integer parameter `n` (`n1`, … if taken): `solve_general(sin(x) =
//!   1/2, x)` is `list(2 n π + π/6, 2 n π + 5π/6)`. Families that agree
//!   modulo a fraction of the period are merged, and a family valid only
//!   for some residues of `n` is restricted to them. Membership is checked
//!   for all integers `n`: the residual is proved periodic in `n` (every
//!   circular argument is affine in `n` with a rational multiple of `π` as
//!   slope), and one period of residues is tested. Other periodic
//!   functions of an inner expression (`sin(x^2) = c`) get their families
//!   through the inverse functions.
//!
//! Systems `solve(list(eqs), list(unknowns))`: a square linear system by
//! Cramer's rule; polynomial systems (parameters allowed) by a lexicographic
//! Gröbner basis in which the symbols other than the unknowns come last, so
//! that an underdetermined system leaves some unknowns *free* (they appear
//! in the answer as themselves) and the rest are read off by
//! back-substitution; two polynomial equations by a resultant; systems with
//! transcendental parts by elimination (an equation linear in an unknown
//! determines it, otherwise an unknown is solved for and the branches
//! followed). Every tuple is substituted back.
//!
//! Conditions in one unknown: `solve(lt(f, g), x)` (likewise `le`, `gt`,
//! `ge`, `ne`) finds the real roots of `f - g`, the poles of its
//! denominator and the zeros of the arguments of logarithms and roots,
//! tests the sign between them and returns a formula such as
//! `or(and(lt(1, x), lt(x, 2)), lt(3, x))`, `true` or `false`. `list(...)`
//! and `and(...)` of conditions intersect, `or(...)` unites, an equation
//! contributes its roots, and `floor(u) = c`, `ceil(u) < c` become interval
//! conditions on `u`.
//!
//! In numeric mode, when no closed form exists, the residual is scanned for
//! sign changes (one equation) or damped Newton runs from a spread of
//! starting points find roots of a system; `nsolve(list(eqs), list(vars),
//! list(guess))` is Newton from the guess. These searches are not
//! exhaustive; the certified operators below are. Polynomials use Sturm
//! sequences.
//!
//! Certified roots (interval Newton / Krawczyk on the interval arithmetic
//! of the tape backend, with the symbolic Jacobian): a scalar `nsolve`
//! returns a ball whose radius is a *proved* bound when a Krawczyk
//! enclosure exists; `nsolve_certified(eq(s), var(s), guess)` returns the
//! enclosing interval(s) `list(lo, hi)` of the unique root in them;
//! `solve_box(eq(s), var(s), list(list(lo, hi), ...))` is an exhaustive
//! bisection with exclusion tests and Krawczyk certification over a box,
//! returning `list(verified boxes, unresolved boxes)`: every root in the
//! region lies in a verified box (each containing exactly one root) unless
//! the second list is non-empty (roots of even multiplicity, singular
//! Jacobians, regions where the functions are undefined, an exhausted
//! budget), so an empty second list is a proof of completeness.
//!
//! Inverse functions are taken on their principal branch by `solve`, so
//! periodic equations yield one representative per branch; `solve_general`
//! gives the whole solution set.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Pat;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::VarNames;
use crate::kernels::polynomial::Polynomial;
use crate::kernels::real_roots;
use crate::kernels::solve::solve_root;

use super::elementary::elementary;
use super::poly::best;
use super::poly::poly;
use super::poly::ratio;
use super::poly::repr::from_term;
use super::poly::repr::to_term;
use super::poly::repr::Gens;
use super::poly::repr::Limits;
use super::poly::repr::Poly;
use super::poly::univariate;

mod certified;
mod complex_roots;
mod decompose;
pub(crate) mod elim;
mod heuristics;
mod inequality;
mod numeric;
mod quintic;
mod tower;
mod normalize;
mod symbolic;
mod system;
mod trig;
#[cfg(test)]
mod battery;

/// Operator attribute: the solutions `u` of `op(u) = ?a`, one pattern per
/// branch that is returned.
#[derive(Clone, Debug)]
pub struct Inverse(pub Vec<Pat>);

/// The equation-solving rule set.
#[must_use]
pub fn solve() -> RuleSet {
    RuleSet::new("solve", install).needs(poly()).needs(elementary()).needs(super::calculus::calculus())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let inverses: [(&str, &[&str]); 13] = [
        ("exp", &["ln(?a)"]),
        ("ln", &["exp(?a)"]),
        ("sin", &["asin(?a)"]),
        ("cos", &["acos(?a)"]),
        ("tan", &["atan(?a)"]),
        ("asin", &["sin(?a)"]),
        ("acos", &["cos(?a)"]),
        ("atan", &["tan(?a)"]),
        ("sinh", &["ln(?a + (?a^2 + 1)^(1/2))"]),
        ("cosh", &["ln(?a + (?a^2 - 1)^(1/2))", "-ln(?a + (?a^2 - 1)^(1/2))"]),
        ("tanh", &["1/2 * ln((1 + ?a)/(1 - ?a))"]),
        ("sqrt", &["?a^2"]),
        ("abs", &["?a", "-?a"]),
    ];
    for (name, branches) in inverses {
        let invalid = |reason: &'static str| RuleError::Invalid { rule: format!("inverse/{name}"), reason };
        let op = i.graph().ops().lookup(name).ok_or_else(|| invalid("unknown operator"))?;
        let mut patterns = Vec::with_capacity(branches.len());
        for text in branches {
            let mut vars = VarNames::default();
            vars.index("a");
            let pat = Pat::parse(text, i.graph(), &mut vars)
                .map_err(|error| RuleError::Parse { rule: format!("inverse/{name}: {text}"), error })?;
            patterns.push(pat);
        }
        i.graph().ops_mut().set_attr(op, Inverse(patterns));
    }
    let solve_op = i.op(OpDescriptor::new("solve", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    let general_op = i.op(OpDescriptor::new("solve_general", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    let nsolve_op = i.op(OpDescriptor::new("nsolve", Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
    let certified_op = i.op(OpDescriptor::new("nsolve_certified", Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
    let box_op = i.op(OpDescriptor::new("solve_box", Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
    let complex_op = i.op(OpDescriptor::new("solve_complex", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("solve/complex", Tier::Reduce, complex_roots::SolveComplex { op: complex_op });
    i.kernel("solve/nsolve-certified", Tier::Reduce, certified::NsolveCertified { op: certified_op });
    i.kernel("solve/box", Tier::Reduce, certified::SolveBox { op: box_op });
    i.kernel("solve/symbolic", Tier::Reduce, Symbolic { solve: solve_op, general: general_op });
    i.kernel("solve/polynomial-numeric", Tier::Reduce, PolynomialNumeric { solve: solve_op });
    i.kernel("solve/nsolve", Tier::Reduce, RootNear { nsolve: nsolve_op });
    i.kernel("solve/numeric", Tier::Reduce, numeric::NumericSolve { solve: solve_op, nsolve: nsolve_op });
    Ok(())
}

/// `l = r` becomes `l - r`; anything else is taken to equal zero.
pub(crate) fn as_expression(
    graph: &mut Graph,
    equation: NodeId,
) -> NodeId {
    match *graph.children(equation) {
        | [lhs, rhs] if graph.op(equation) == core::EQ => difference(graph, lhs, rhs),
        | _ => equation,
    }
}

fn difference(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    let negated = graph.node(core::MUL, &[minus_one, b]);
    graph.node(core::ADD, &[a, negated])
}

fn product(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn reciprocal(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    graph.node(core::POW, &[node, minus_one])
}

fn rational(
    graph: &mut Graph,
    value: BigRational,
) -> NodeId {
    graph.num(Number::rat(value))
}

const MAX_DEPTH: usize = 8;

thread_local! {
    static GENERAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

thread_local! {
    static NUMERIC_RUN: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether the engine run asks for numbers: exact `root_of` objects are
/// then replaced by their values.
pub(crate) fn numeric_run() -> bool {
    NUMERIC_RUN.with(std::cell::Cell::get)
}

/// Sets [`numeric_run`] for the lifetime of the guard.
struct NumericGuard(bool);

impl NumericGuard {
    fn set(value: bool) -> Self {
        Self(NUMERIC_RUN.with(|n| n.replace(value)))
    }
}

impl Drop for NumericGuard {
    fn drop(&mut self) {
        NUMERIC_RUN.with(|n| n.set(self.0));
    }
}

/// Whether solutions of periodic equations are being collected as families
/// with an integer parameter.
pub(crate) fn general_mode() -> bool {
    GENERAL.with(std::cell::Cell::get)
}

/// Runs `f` with the general mode set.
fn with_general<T>(f: impl FnOnce() -> T) -> T {
    let previous = GENERAL.with(|g| g.replace(true));
    let out = f();
    GENERAL.with(|g| g.set(previous));
    out
}

pub(super) fn cap_default() -> usize {
    Limits::default().terms
}

/// Solves `expr = 0` for `x`. `None` means "could not solve"; `Some` is
/// the complete list of solutions that were found on principal branches.
pub(crate) fn solve_for(
    graph: &mut Graph,
    expr: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    if depth > MAX_DEPTH {
        return None;
    }
    let symbol = graph.symbol_of(x)?;
    let term = best(graph, expr)?;
    let term = sqrt_as_power(graph, term);
    if !graph.depends_on(graph.find(term), symbol) {
        // `0 = 0` is solved by everything, which a list cannot express.
        return if graph.number_of(term).is_some_and(Number::is_zero) { None } else { Some(Vec::new()) };
    }
    // Write exp(n*u) as exp(u)^n so that `exp(2*x) - 3*exp(x) + 2` is a
    // polynomial in the single generator exp(x).
    let term = normalize::exponentials(graph, term, symbol);
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let fraction = ratio(graph, &mut gens, term, Limits::default())?;
    let dependent: Vec<u32> = fraction
        .numer
        .support()
        .into_iter()
        .filter(|&g| gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)))
        .collect();
    let split = if dependent.len() >= 2 { split_by_factors(graph, &gens, &fraction.numer, x, depth) } else { None };
    let candidates = match dependent.as_slice() {
        | _ if split.is_some() => split.unwrap_or_default(),
        | [] => Vec::new(),
        | &[g] if g == gx => polynomial_roots(graph, &gens, &fraction.numer, gx)?,
        | &[g] => {
            // A polynomial in one function of x: solve for the function,
            // then invert it.
            let inner = gens.node(g)?;
            let mut out = Some(Vec::new());
            for value in polynomial_roots(graph, &gens, &fraction.numer, g)? {
                match invert(graph, inner, value, x, depth + 1) {
                    | Some(found) => out.as_mut()?.extend(found),
                    | None => {
                        out = None;
                        break;
                    },
                }
            }
            match out {
                | Some(found) => found,
                | None => heuristics::multi_generator(graph, term, x, depth)?,
            }
        },
        | _ => heuristics::multi_generator(graph, term, x, depth)?,
    };
    // Discard candidates at which the original expression is not zero
    // (poles of a cancelled denominator, roots introduced by squaring),
    // wherever that can be checked numerically.
    let mut solutions = filter_candidates(graph, term, x, candidates);
    // Ascending order when every solution is a number.
    let values: Option<Vec<f64>> = solutions.iter().map(|&s| graph.eval(s, &Env::numeric(0.0))).collect();
    if let Some(values) = values {
        let mut keyed: Vec<(f64, NodeId)> = values.into_iter().zip(solutions).collect();
        keyed.sort_by(|a, b| a.0.total_cmp(&b.0));
        keyed.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-12 * (1.0 + a.0.abs()));
        solutions = keyed.into_iter().map(|(_, node)| node).collect();
    }
    Some(solutions)
}

/// Whether a polynomial is small enough for the factorisation to be worth
/// attempting (the reduction to one variable is exponential in the number
/// of variables and their degrees).
pub(crate) fn factorable(poly: &Poly) -> bool {
    poly.len() <= 40 && poly.total_degree() <= 12 && poly.support().len() <= 5
}

/// When the numerator factors over `Q` into several pieces that involve
/// the unknown, the solutions are those of the pieces.
fn split_by_factors(
    graph: &mut Graph,
    gens: &Gens,
    numer: &Poly,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    let symbol = graph.symbol_of(x)?;
    if !factorable(numer) {
        return None;
    }
    let vars = numer.support();
    let (_, factors) = crate::rules::poly::multifactor::factor(numer, &vars)?;
    let pieces: Vec<Poly> = factors
        .into_iter()
        .map(|(f, _)| f)
        .filter(|f| f.support().into_iter().any(|g| gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol))))
        .collect();
    if pieces.len() < 2 {
        return None;
    }
    let mut out = Vec::new();
    for piece in pieces {
        let term = to_term(graph, gens, &piece);
        out.extend(solve_for(graph, term, x, depth + 1)?);
    }
    Some(out)
}

/// Base magnitudes of the spot-check assignments.
const BASE: [f64; 6] = [1.3, 0.7, 2.1, 0.45, 1.7, 0.9];

/// Generic assignments of the free symbols of `nodes` (other than `skip`)
/// for numeric spot checks of formulas in several parameters.
pub(crate) fn sample_envs(
    graph: &Graph,
    nodes: &[NodeId],
    skip: Option<SymbolId>,
) -> Vec<Env> {
    let mut symbols: Vec<SymbolId> = Vec::new();
    for &node in nodes {
        for &s in graph.free_symbols(graph.find(node)) {
            if Some(s) != skip && !symbols.contains(&s) {
                symbols.push(s);
            }
        }
    }
    (0..4_usize)
        .map(|k| {
            let mut env = Env::numeric(0.0);
            for (i, &s) in symbols.iter().enumerate() {
                if graph.assumption(s).has(Facts::INTEGER) {
                    const INTEGERS: [f64; 4] = [1.0, 0.0, -1.0, 2.0];
                    env.bind(s, INTEGERS.get((k + i) % 4).copied().unwrap_or(0.0));
                    continue;
                }
                let magnitude = BASE.get((2 * i + k) % 6).copied().unwrap_or(1.0) * (1.0 + 0.11 * f64::from(u8::try_from(k).unwrap_or(0)));
                let sign = if k > 0 && (i + k) % 3 == 0 { -1.0 } else { 1.0 };
                env.bind(s, sign * magnitude);
            }
            env
        })
        .collect()
}

/// The candidates that make `term` vanish. Candidates with parameters are
/// spot-checked at several assignments: one that is finite and wrong at
/// every assignment where it is defined is discarded, one that is undefined
/// everywhere is kept (it may be real elsewhere in the parameter space).
pub(crate) fn filter_candidates(
    graph: &mut Graph,
    term: NodeId,
    x: NodeId,
    candidates: Vec<NodeId>,
) -> Vec<NodeId> {
    let symbol = graph.symbol_of(x);
    let mut solutions: Vec<NodeId> = Vec::new();
    for candidate in candidates {
        let substituted = graph.substitute(term, x, candidate);
        let free = graph.free_symbols(graph.find(substituted)).to_vec();
        let keep = if free.is_empty() {
            match graph.eval(substituted, &Env::numeric(0.0)) {
                | Some(v) => v.abs() <= 1e-7,
                | None => true,
            }
        } else {
            // Right at some assignment, or undefined at all of them. A
            // formula can be wrong at others (principal branches, regimes
            // of the parameters), but not everywhere it is defined.
            let (mut right, mut wrong) = (0_u32, 0_u32);
            for env in sample_envs(graph, &[substituted, candidate], symbol) {
                let (Some(residual), Some(at)) = (graph.eval(substituted, &env), graph.eval(candidate, &env)) else {
                    continue;
                };
                if residual.is_finite() && at.is_finite() {
                    if residual.abs() > 1e-7 * (1.0 + at.abs().powi(3)) {
                        wrong += 1;
                    } else {
                        right += 1;
                    }
                }
            }
            right > 0 || wrong == 0
        };
        if keep && !solutions.iter().any(|&s| graph.same(s, candidate) || s == candidate) {
            solutions.push(candidate);
        }
    }
    solutions
}

/// All solutions of `expr = 0` for `x`, periodic families included: each
/// family carries an integer parameter `n` (`n1`, ... if `n` is taken).
pub(crate) fn solve_general_for(
    graph: &mut Graph,
    expr: NodeId,
    x: NodeId,
) -> Option<Vec<NodeId>> {
    with_general(|| {
        let term = best(graph, expr)?;
        let term = sqrt_as_power(graph, term);
        if let Some(found) = trig::solve(graph, term, x, true) {
            return Some(found);
        }
        solve_for(graph, expr, x, 0)
    })
}

/// `sqrt(u)` written `u^(1/2)` throughout, so that one radical heuristic
/// covers both spellings.
fn sqrt_as_power(
    graph: &mut Graph,
    term: NodeId,
) -> NodeId {
    let Some(sqrt) = graph.ops().lookup("sqrt") else {
        return term;
    };
    let roots = heuristics::nodes_with(graph, term, |g, n| g.op(n) == sqrt);
    let Some(half) = Number::fraction(1, 2) else {
        return term;
    };
    let half = graph.num(half);
    let mut out = term;
    for r in roots {
        let &[u] = graph.children(r) else {
            continue;
        };
        let power = graph.node(core::POW, &[u, half]);
        out = graph.replace_subterm(out, r, power);
    }
    out
}

/// Solves `g(x) = value` for `x`, where `g` is the term `inner`.
fn invert(
    graph: &mut Graph,
    inner: NodeId,
    value: NodeId,
    x: NodeId,
    depth: usize,
) -> Option<Vec<NodeId>> {
    if depth > MAX_DEPTH {
        return None;
    }
    let symbol = graph.symbol_of(x)?;
    if graph.as_symbol(inner) == Some(symbol) {
        return Some(vec![value]);
    }
    let op = graph.op(inner);
    let args = graph.children(inner).to_vec();
    let depends = |graph: &Graph, node: NodeId| graph.depends_on(graph.find(node), symbol);
    let mut targets: Vec<(NodeId, NodeId)> = Vec::new();
    if let (true, &[base, exp]) = (op == core::POW, args.as_slice()) {
        match (depends(graph, base), depends(graph, exp)) {
            | (true, false) => {
                // base^exp = value  =>  base = value^(1/exp), and its
                // negative for an even integer exponent.
                let inverse_exp = reciprocal(graph, exp);
                let root = graph.node(core::POW, &[value, inverse_exp]);
                targets.push((base, root));
                let even = graph.number_of(exp).and_then(Number::to_i64).is_some_and(|n| n % 2 == 0);
                if even {
                    let minus_one = graph.int(-1);
                    targets.push((base, graph.node(core::MUL, &[minus_one, root])));
                }
            },
            | (false, true) => {
                // base^exp = value  =>  exp = ln(value)/ln(base).
                let ln = graph.ops().lookup("ln")?;
                let (top, bottom) = (graph.node(ln, &[value]), graph.node(ln, &[base]));
                let inverse = reciprocal(graph, bottom);
                targets.push((exp, graph.node(core::MUL, &[top, inverse])));
            },
            | _ => return None,
        }
    } else if let &[arg] = args.as_slice() {
        let periodic = if general_mode() { periodic_targets(graph, op, arg, value, inner, x) } else { None };
        match periodic {
            | Some(found) => targets.extend(found),
            | None => {
                let branches = graph.ops().attr::<Inverse>(op)?.0.clone();
                for branch in branches {
                    targets.push((arg, branch.instantiate(graph, &[value])?));
                }
            },
        }
    } else {
        return None;
    }
    let mut out = Vec::new();
    for (unknown, target) in targets {
        let equation = difference(graph, unknown, target);
        out.extend(solve_for(graph, equation, x, depth + 1)?);
    }
    Some(out)
}

/// The solutions `u` of `op(u) = value` as families with an integer
/// parameter, for the periodic circular functions.
fn periodic_targets(
    graph: &mut Graph,
    op: OpId,
    arg: NodeId,
    value: NodeId,
    inner: NodeId,
    x: NodeId,
) -> Option<Vec<(NodeId, NodeId)>> {
    let names = ["sin", "cos", "tan"];
    let which = names.iter().position(|n| graph.ops().lookup(n) == Some(op))?;
    let (pi_op, inverse) = (graph.ops().lookup("pi")?, graph.ops().lookup(["asin", "acos", "atan"][which])?);
    let n = trig::integer_symbol(graph, &[value, inner, x]);
    let pi = graph.node(pi_op, &[]);
    let principal = graph.node(inverse, &[value]);
    let minus_one = graph.int(-1);
    let two = graph.int(2);
    let mut out = Vec::new();
    match which {
        | 0 => {
            let turn = product(graph, &[two, pi, n]);
            out.push((arg, graph.node(core::ADD, &[principal, turn])));
            let neg = graph.node(core::MUL, &[minus_one, principal]);
            out.push((arg, graph.node(core::ADD, &[pi, neg, turn])));
        },
        | 1 => {
            let turn = product(graph, &[two, pi, n]);
            out.push((arg, graph.node(core::ADD, &[principal, turn])));
            let neg = graph.node(core::MUL, &[minus_one, principal]);
            out.push((arg, graph.node(core::ADD, &[neg, turn])));
        },
        | _ => {
            let turn = product(graph, &[pi, n]);
            out.push((arg, graph.node(core::ADD, &[principal, turn])));
        },
    }
    Some(out)
}

/// Roots of `poly` viewed as a polynomial in `generator`, as terms.
fn polynomial_roots(
    graph: &mut Graph,
    gens: &Gens,
    poly: &Poly,
    generator: u32,
) -> Option<Vec<NodeId>> {
    let coefficients = poly.coefficients_in(generator);
    // Roots at zero: strip the power of the generator common to all terms.
    let low = coefficients.iter().position(|c| !c.is_zero())?;
    let mut roots = Vec::new();
    if low > 0 {
        roots.push(graph.int(0));
    }
    let coefficients = coefficients.get(low..)?;
    if coefficients.len() == 1 {
        return Some(roots);
    }
    let exact: Option<Vec<BigRational>> =
        coefficients.iter().map(|c| c.as_constant().and_then(|n| n.to_rational())).collect();
    if let Some(numbers) = exact {
        roots.extend(rational_polynomial_roots(graph, &numbers)?);
        return Some(roots);
    }
    // Symbolic coefficients: factor over the parameters, then formulas.
    symbolic::roots(graph, gens, poly, generator, 0)
}

/// Real roots of a polynomial with rational coefficients (ascending
/// degree), exactly. `None` when some real root has no expression in
/// radicals that this routine knows how to write.
fn rational_polynomial_roots(
    graph: &mut Graph,
    coefficients: &[BigRational],
) -> Option<Vec<NodeId>> {
    let mut roots: Vec<(f64, NodeId)> = Vec::new();
    let (_, factors) = univariate::factor(coefficients);
    for (factor, _) in factors {
        match factor.as_slice() {
            | [c0, c1] => {
                let value = BigRational::new(-c0.clone(), c1.clone());
                let approx = Number::rat(value.clone()).to_f64();
                roots.push((approx, rational(graph, value)));
            },
            | [c, b, a] => {
                let discriminant = b * b - BigInt::from(4) * a * c;
                if discriminant.is_negative() {
                    continue;
                }
                // (-b ± sqrt(D)) / (2a), with the square part of D pulled out.
                let (outside, inside) = square_part(&discriminant);
                let denominator = BigInt::from(2) * a;
                let centre = BigRational::new(-b.clone(), denominator.clone());
                let scale = BigRational::new(outside, denominator);
                let half = graph.num(Number::fraction(1, 2)?);
                let radicand = graph.num(Number::Int(inside.clone()));
                let radical = graph.node(core::POW, &[radicand, half]);
                for sign in [-1, 1] {
                    let coefficient = rational(graph, &scale * BigRational::from_integer(BigInt::from(sign)));
                    let offset = product(graph, &[coefficient, radical]);
                    let centre_node = rational(graph, centre.clone());
                    let root =
                        if centre.is_zero() { offset } else { graph.node(core::ADD, &[centre_node, offset]) };
                    let approx = Number::rat(centre.clone()).to_f64()
                        + f64::from(sign)
                            * Number::rat(scale.clone()).to_f64()
                            * Number::Int(inside.clone()).to_f64().sqrt();
                    roots.push((approx, root));
                }
            },
            | higher => {
                // c0 + cn*x^n: the real n-th roots of -c0/cn.
                let n = higher.len() - 1;
                let (c0, cn) = (higher.first()?, higher.last()?);
                if higher.get(1..n)?.iter().any(|c| !c.is_zero()) {
                    if n >= 5 {
                        // Radicals for solvable quintics, `root_of` objects
                        // for the rest.
                        roots.extend(quintic::roots(graph, higher)?);
                        continue;
                    }
                    // Cubic and quartic formulas.
                    let rational: Vec<BigRational> = higher.iter().cloned().map(BigRational::from_integer).collect();
                    roots.extend(heuristics::higher_degree(graph, &rational)?);
                    continue;
                }
                let value = BigRational::new(-c0.clone(), cn.clone());
                let exponent = graph.num(Number::fraction(1, i64::try_from(n).ok()?)?);
                let magnitude = rational(graph, value.abs());
                let root = graph.node(core::POW, &[magnitude, exponent]);
                let approx = Number::rat(value.abs()).to_f64().powf(1.0 / n as f64);
                let minus_one = graph.int(-1);
                match (n % 2 == 0, value.is_negative()) {
                    | (true, true) => {},
                    | (true, false) => {
                        roots.push((-approx, product(graph, &[minus_one, root])));
                        roots.push((approx, root));
                    },
                    | (false, false) => roots.push((approx, root)),
                    | (false, true) => roots.push((-approx, product(graph, &[minus_one, root]))),
                }
            },
        }
    }
    roots.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(roots.into_iter().map(|(_, node)| node).collect())
}

/// Writes a non-negative integer as `outside^2 * inside` with `inside`
/// square-free as far as small prime factors go.
fn square_part(value: &BigInt) -> (BigInt, BigInt) {
    let mut outside = BigInt::from(1);
    let mut inside = value.clone();
    let mut p = BigInt::from(2);
    // Trial division is enough: discriminants of factors that survived
    // factorisation are small in practice, and a missed square only leaves
    // the radical less simplified, never wrong.
    for _ in 0..10_000 {
        let square = &p * &p;
        if square > inside {
            break;
        }
        while (&inside % &square).is_zero() {
            inside /= &square;
            outside *= &p;
        }
        p += 1;
    }
    (outside, inside)
}

// ----------------------------------------------------------------------
// Systems
// ----------------------------------------------------------------------

/// Determinant by cofactor expansion along the first row.
fn determinant(
    matrix: &[Vec<Poly>],
    cap: usize,
) -> Option<Poly> {
    match matrix {
        | [] => Some(Poly::constant(Number::from(1))),
        | [row] => row.first().cloned(),
        | _ => {
            let mut total = Poly::zero();
            for (column, entry) in matrix.first()?.iter().enumerate() {
                if entry.is_zero() {
                    continue;
                }
                let minor: Vec<Vec<Poly>> = matrix
                    .get(1..)?
                    .iter()
                    .map(|row| {
                        row.iter().enumerate().filter(|&(c, _)| c != column).map(|(_, p)| p.clone()).collect()
                    })
                    .collect();
                let term = entry.mul(&determinant(&minor, cap)?, cap)?;
                total = if column % 2 == 0 { total.add(&term) } else { total.sub(&term) };
            }
            Some(total)
        },
    }
}

/// Solves a square linear system by Cramer's rule with polynomial
/// arithmetic in the parameters. `None` if the system is not linear in the
/// unknowns, not square, too large, or singular.
pub(crate) fn solve_linear(
    graph: &mut Graph,
    equations: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<NodeId>> {
    let n = unknowns.len();
    if equations.len() != n || n == 0 || n > 6 {
        return None;
    }
    let symbols: Vec<SymbolId> = unknowns.iter().map(|&u| graph.symbol_of(u)).collect::<Option<_>>()?;
    let mut gens = Gens::default();
    let indices: Vec<u32> = unknowns.iter().map(|&u| gens.index(graph, u)).collect();
    let limits = Limits::default();
    let mut matrix: Vec<Vec<Poly>> = Vec::with_capacity(n);
    let mut rhs: Vec<Poly> = Vec::with_capacity(n);
    for &equation in equations {
        let expr = as_expression(graph, equation);
        let term = best(graph, expr)?;
        let poly = from_term(graph, &mut gens, term, limits)?;
        // Parameters must not hide an unknown.
        for g in poly.support() {
            if indices.contains(&g) {
                continue;
            }
            let node = gens.node(g)?;
            if symbols.iter().any(|&s| graph.depends_on(graph.find(node), s)) {
                return None;
            }
        }
        let mut row = vec![Poly::zero(); n];
        let mut constant = Poly::zero();
        for (mono, coeff) in poly.terms() {
            let unknown_part: Vec<(u32, u32)> = mono.iter().copied().filter(|(g, _)| indices.contains(g)).collect();
            let parameters: Vec<(u32, u32)> = mono.iter().copied().filter(|(g, _)| !indices.contains(g)).collect();
            let piece = Poly::monomial(parameters, coeff.clone());
            match unknown_part.as_slice() {
                | [] => constant = constant.add(&piece),
                | &[(g, 1)] => {
                    let column = indices.iter().position(|&i| i == g)?;
                    let slot = row.get_mut(column)?;
                    *slot = slot.add(&piece);
                },
                | _ => return None,
            }
        }
        matrix.push(row);
        rhs.push(constant.neg());
    }
    let det = determinant(&matrix, limits.terms)?;
    if det.is_zero() {
        return None;
    }
    // Fix the sign of the common denominator so that `1/(a + 1)` is not
    // written `-1/(-a - 1)`.
    let flip = det.terms().last().is_some_and(|(_, c)| c.is_negative());
    let sign = Number::from(if flip { -1 } else { 1 });
    let det = det.scale(&sign);
    let mut solution = Vec::with_capacity(n);
    for column in 0..n {
        let mut replaced = matrix.clone();
        for (row, value) in replaced.iter_mut().zip(&rhs) {
            *row.get_mut(column)? = value.clone();
        }
        let numerator = determinant(&replaced, limits.terms)?.scale(&sign);
        let node = match det.as_constant() {
            | Some(c) => to_term(graph, &gens, &numerator.scale(&c.recip()?)),
            | None => {
                let top = to_term(graph, &gens, &numerator);
                let bottom = to_term(graph, &gens, &det);
                let inverse = reciprocal(graph, bottom);
                product(graph, &[top, inverse])
            },
        };
        solution.push(node);
    }
    Some(solution)
}

/// Solves a triangular system for the last variable first.
pub(super) fn back_substitute(
    graph: &mut Graph,
    basis: &[NodeId],
    variables: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    // A non-zero constant in the basis means the ideal is everything: the
    // system is inconsistent.
    if basis.iter().any(|&b| graph.number_of(b).is_some_and(|n| !n.is_zero())) {
        return Some(Vec::new());
    }
    let Some((&last, earlier)) = variables.split_last() else {
        // Nothing left to determine: consistent unless a non-zero constant
        // remains.
        let contradiction = basis.iter().any(|&b| graph.number_of(b).is_some_and(|n| !n.is_zero()));
        return Some(if contradiction { Vec::new() } else { vec![Vec::new()] });
    };
    let earlier_symbols: Vec<SymbolId> = earlier.iter().map(|&v| graph.symbol_of(v)).collect::<Option<_>>()?;
    let last_symbol = graph.symbol_of(last)?;
    let only_last = |graph: &Graph, node: NodeId| {
        let class = graph.find(node);
        graph.depends_on(class, last_symbol) && !earlier_symbols.iter().any(|&s| graph.depends_on(class, s))
    };
    let univariate = basis.iter().copied().find(|&b| only_last(graph, b))?;
    let mut out = Vec::new();
    for value in solve_for(graph, univariate, last, 0)? {
        let mut remaining = Vec::new();
        for &b in basis {
            if only_last(graph, b) {
                continue;
            }
            let substituted = graph.substitute(b, last, value);
            // Expand to see whether the equation has become trivial.
            let mut gens = Gens::default();
            let poly = from_term(graph, &mut gens, substituted, Limits::default())?;
            if !poly.is_zero() {
                remaining.push(to_term(graph, &gens, &poly));
            }
        }
        for mut partial in back_substitute(graph, &remaining, earlier)? {
            partial.push(value);
            out.push(partial);
        }
    }
    Some(out)
}

struct Symbolic {
    solve: OpId,
    general: OpId,
}

impl Kernel for Symbolic {
    fn ops(&self) -> Vec<OpId> {
        vec![self.solve, self.general]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let numeric_run = cx.env.numeric;
        let _numeric = NumericGuard::set(numeric_run);
        let graph = &mut *cx.graph;
        let general = graph.op(node) == self.general;
        let &[equation, unknown] = graph.children(node) else {
            return Outcome::Pass;
        };
        let solutions = if graph.op(unknown) == core::LIST {
            if graph.op(equation) != core::LIST {
                return Outcome::Pass;
            }
            let equations = graph.children(equation).to_vec();
            let unknowns = graph.children(unknown).to_vec();
            let tuples = system::solve_system(graph, &equations, &unknowns);
            tuples.map(|tuples| tuples.iter().map(|t| graph.node(core::LIST, t)).collect::<Vec<_>>())
        } else if inequality::is_condition(graph, equation) {
            return inequality::solve_condition(graph, equation, unknown).map_or(Outcome::Pass, Outcome::Equal);
        } else {
            let expr = as_expression(graph, equation);
            if general { solve_general_for(graph, expr, unknown) } else { solve_for(graph, expr, unknown, 0) }
        };
        // Not pinned: the solutions are ordinary terms and should be
        // simplified like any others.
        match solutions {
            | Some(list) => Outcome::Equal(graph.node(core::LIST, &list)),
            | None => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// The numeric coefficients (highest degree first) of `expr` as a
/// polynomial in `x` under the bindings of `env`.
fn numeric_polynomial(
    graph: &mut Graph,
    expr: NodeId,
    x: NodeId,
    env: &Env,
) -> Option<Vec<f64>> {
    let term = best(graph, expr)?;
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let fraction = ratio(graph, &mut gens, term, Limits::default())?;
    let mut coefficients = Vec::new();
    for coefficient in fraction.numer.coefficients_in(gx) {
        let node = to_term(graph, &gens, &coefficient);
        coefficients.push(graph.eval(node, env)?);
    }
    while coefficients.last().is_some_and(|c| *c == 0.0) {
        coefficients.pop();
    }
    coefficients.reverse();
    (coefficients.len() > 1).then_some(coefficients)
}

/// Numeric phase: all real roots of a polynomial equation that the
/// symbolic kernel could not express.
struct PolynomialNumeric {
    solve: OpId,
}

impl Kernel for PolynomialNumeric {
    fn ops(&self) -> Vec<OpId> {
        vec![self.solve]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        // Already evaluated in this run.
        if cx.graph.approx(cx.graph.find(node)).is_some() {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let &[equation, unknown] = graph.children(node) else {
            return Outcome::Pass;
        };
        let heavy = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).flags.has(OpFlags::HEAVY);
        if graph.enodes(graph.find(node)).any(|n| !heavy(graph, n)) || graph.op(unknown) == core::LIST {
            return Outcome::Pass;
        }
        let expr = as_expression(graph, equation);
        let Some(coefficients) = numeric_polynomial(graph, expr, unknown, cx.env) else {
            return Outcome::Pass;
        };
        let tolerance = cx.env.tolerance.max(1e-13);
        let polynomial = Polynomial::new(coefficients);
        let Ok(roots) = real_roots::find_roots(&polynomial, tolerance) else {
            return Outcome::Pass;
        };
        let nodes: Vec<NodeId> = roots.into_iter().map(|r| graph.float(r)).collect();
        Outcome::Equal(graph.node(core::LIST, &nodes))
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `nsolve(equation, x, guess)`: one root near `guess`, as a witness.
struct RootNear {
    nsolve: OpId,
}

impl Kernel for RootNear {
    fn ops(&self) -> Vec<OpId> {
        vec![self.nsolve]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        // Already evaluated in this run.
        if cx.graph.approx(cx.graph.find(node)).is_some() {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let &[equation, unknown, guess] = graph.children(node) else {
            return Outcome::Pass;
        };
        let expr = as_expression(graph, equation);
        let (Some(x), Some(term)) = (graph.symbol_of(unknown), best(graph, expr)) else {
            return Outcome::Pass;
        };
        let Some(start) = best(graph, guess).and_then(|g| graph.eval(g, cx.env)) else {
            return Outcome::Pass;
        };
        let f = |t: f64| {
            let mut env = cx.env.clone();
            env.bind(x, t);
            graph.eval(term, &env).unwrap_or(f64::NAN)
        };
        let tolerance = cx.env.tolerance.max(1e-14);
        match solve_root(f, start, tolerance, 200) {
            | Ok(root) => {
                // A verified enclosure (interval Krawczyk) replaces the
                // error estimate when one exists.
                let bindings = cx.env.bindings().to_vec();
                let ball = certified::certified_ball(graph, equation, x, &bindings, root)
                    .unwrap_or(Ball { mid: root, rad: tolerance });
                Outcome::Approx(ball)
            },
            | Err(_) => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[solve()], src)
    }

    /// Every listed solution, substituted back, must satisfy the equation.
    fn check_solutions(
        equation: &str,
        solutions: &str,
        expected_count: usize,
    ) {
        let inner = solutions.strip_prefix("list(").and_then(|s| s.strip_suffix(')')).unwrap_or("");
        let mut items = Vec::new();
        let (mut depth, mut start) = (0_i32, 0_usize);
        for (i, ch) in inner.char_indices() {
            match ch {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => {
                    items.push(inner[start..i].trim());
                    start = i + 1;
                },
                | _ => {},
            }
        }
        if !inner[start..].trim().is_empty() {
            items.push(inner[start..].trim());
        }
        assert_eq!(items.len(), expected_count, "{equation} -> {solutions}");
        let (lhs, rhs) = equation.split_once('=').unwrap_or((equation, "0"));
        for item in items {
            use std::fmt::Write as _;
            // Replace the unknown `x`, but not the x in `exp`.
            let template = format!("({lhs}) - ({rhs})");
            let chars: Vec<char> = template.chars().collect();
            let mut residual = String::new();
            for (i, &ch) in chars.iter().enumerate() {
                let alone = ch == 'x'
                    && !(i > 0 && chars[i - 1].is_alphanumeric())
                    && !chars.get(i + 1).is_some_and(|c| c.is_alphanumeric());
                if alone {
                    let _ = write!(residual, "({item})");
                } else {
                    residual.push(ch);
                }
            }
            let value = eval(&[solve()], &residual, &[]);
            assert!(value.abs() < 1e-9, "{item} does not solve {equation}: residual {value}");
        }
    }

    /// Splits `list(a, b, ...)` at top-level commas.
    fn list_items(text: &str) -> Vec<String> {
        let inner = text.strip_prefix("list(").and_then(|t| t.strip_suffix(')')).unwrap_or("");
        let (mut depth, mut current, mut parts) = (0_i32, String::new(), Vec::new());
        for ch in inner.chars() {
            match ch {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => {
                    parts.push(std::mem::take(&mut current));
                    continue;
                },
                | _ => {},
            }
            current.push(ch);
        }
        if !current.trim().is_empty() {
            parts.push(current);
        }
        parts
    }

    #[test]
    fn symbolic_quartic_by_ferrari() {
        let rules = crate::rules::standard();
        let (text, reduced) = reduce_with(&rules, "solve(x^4 + a*x^3 + b*x + c, x)", &[]);
        assert!(reduced, "{text}");
        let items = list_items(&text);
        assert!(items.len() >= 4, "{text}");
        // At each parameter point every defined formula must be a root.
        let mut found = Vec::new();
        for (a, b, c) in [(-2.0, 2.0, -5.0), (1.0, -1.5, -2.0), (0.5, 3.0, -4.0)] {
            let mut defined = 0;
            for item in &items {
                let at = eval(&rules, item, &[("a", a), ("b", b), ("c", c)]);
                if !at.is_finite() {
                    continue;
                }
                defined += 1;
                let residual = at.powi(4) + a * at.powi(3) + b * at + c;
                assert!(residual.abs() < 1e-7 * (1.0 + at.abs().powi(4)), "{item} at ({a},{b},{c}) gives {residual}");
                found.push(at);
            }
            assert!(defined >= 2, "too few defined roots at ({a},{b},{c}): {text}");
        }
        assert!(found.len() >= 6);
    }

    /// Every member of every family (`n` in -20..=20) must satisfy `f`.
    fn check_families(
        src: &str,
        f: impl Fn(f64) -> f64,
        at_least: usize,
    ) {
        let rules = crate::rules::standard();
        let (text, reduced) = reduce_with(&rules, src, &[]);
        assert!(reduced, "{src}: {text}");
        let items = list_items(&text);
        assert!(items.len() >= at_least, "{src}: {text}");
        for item in &items {
            for n in -20..=20 {
                let x = eval(&rules, item, &[("n", f64::from(n))]);
                assert!(f(x).abs() < 1e-8, "{item} at n = {n} (x = {x}) is not a solution of {src}: {text}");
            }
        }
    }

    #[test]
    fn trigonometric_generalisations() {
        // Different shifts in one equation.
        check_families("solve_general(sin(x) = sin(2*x + 1), x)", |x| x.sin() - (2.0 * x + 1.0).sin(), 2);
        check_families("solve_general(sin(x + 1) = cos(2*x), x)", |x| (x + 1.0).sin() - (2.0 * x).cos(), 2);
        // Large multipliers.
        check_families("solve_general(sin(30*x) = sin(2*x), x)", |x| (30.0 * x).sin() - (2.0 * x).sin(), 2);
        check_families("solve_general(cos(27*x) = 1, x)", |x| (27.0 * x).cos() - 1.0, 1);
        // Rational multipliers with different denominators.
        check_families("solve_general(sin(x/3) = sin(x/2), x)", |x| (x / 3.0).sin() - (x / 2.0).sin(), 2);
    }

    #[test]
    fn sextics_and_septics_by_galois_group() {
        let rules = crate::rules::standard();
        // (equation, ascending coefficients, number of real roots)
        let cases: [(&str, Vec<f64>, usize); 7] = [
            // Cyclic septic (Gauss periods, conductor 29): seven real roots.
            ("solve(x^7 + x^6 - 12*x^5 - 7*x^4 + 28*x^3 + 14*x^2 - 9*x + 1, x)", vec![1.0, -9.0, 14.0, 28.0, -7.0, -12.0, 1.0, 1.0], 7),
            // 2 cos(2 pi / 13): cyclic sextic, six real roots.
            ("solve(x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1, x)", vec![-1.0, 3.0, 6.0, -4.0, -5.0, 1.0, 1.0], 6),
            // Cyclic sextic of composite conductor 91 = 7 * 13.
            ("solve(x^6 - x^5 - 14*x^4 + 9*x^3 + 35*x^2 - 16*x - 1, x)", vec![-1.0, -16.0, 35.0, 9.0, -14.0, -1.0, 1.0], 6),
            // D7(x) - 3: group F42, one real root through Lagrange resolvents.
            ("solve(x^7 - 7*x^5 + 14*x^3 - 7*x - 3, x)", vec![-3.0, -7.0, 0.0, 14.0, 0.0, -7.0, 0.0, 1.0], 1),
            // D7(x) - 1/4: group F21, seven real roots (trigonometric form).
            ("solve(4*x^7 - 28*x^5 + 56*x^3 - 28*x - 1, x)", vec![-1.0, -28.0, 0.0, 56.0, 0.0, -28.0, 0.0, 4.0], 7),
            // Sextic with group A4 on six points (block tower).
            ("solve(x^6 - 12*x^4 - 64*x^3 - 144*x^2 + 1728, x)", vec![1728.0, 0.0, -144.0, -64.0, -12.0, 0.0, 1.0], 0),
            // S4 on the six directed 4-cycles (block tower).
            ("solve(x^6 + 6*x^5 + 24*x^4 + 56*x^3 + 160*x^2 + 224*x + 384, x)", vec![384.0, 224.0, 160.0, 56.0, 24.0, 6.0, 1.0], 0),
        ];
        for (src, poly, count) in cases {
            let (text, reduced) = reduce_with(&rules, src, &[]);
            eprintln!("{src}: {text}");
            assert!(reduced, "{src}: {text}");
            let items = list_items(&text);
            if count > 0 {
                assert_eq!(items.len(), count, "{src}: {text}");
                assert!(!text.contains("root_of"), "{src}: {text}");
            }
            for item in &items {
                let x = eval(&rules, item, &[]);
                let residual: f64 = poly.iter().rev().fold(0.0, |acc, c| acc * x + c);
                assert!(residual.abs() < 1e-7 * (1.0 + x.abs().powi(7)), "{src}: {item} -> {x} gives {residual}");
            }
        }
    }

    #[test]
    fn quintics_and_root_objects() {
        let rules = crate::rules::standard();
        // Solvable: radicals, verified by substitution.
        for (src, poly) in [
            ("solve(x^5 - 5*x + 12, x)", (vec![12.0, -5.0, 0.0, 0.0, 0.0, 1.0])),
            ("solve(x^5 + 15*x + 12, x)", (vec![12.0, 15.0, 0.0, 0.0, 0.0, 1.0])),
            ("solve(x^5 + 5*x^3 + 5*x + 3, x)", (vec![3.0, 5.0, 0.0, 5.0, 0.0, 1.0])),
        ] {
            let (text, reduced) = reduce_with(&rules, src, &[]);
            assert!(reduced, "{src}: {text}");
            assert!(!text.contains("root_of"), "{src} should be in radicals: {text}");
            let items = list_items(&text);
            assert_eq!(items.len(), 1, "{src}: {text}");
            let x = eval(&rules, &items[0], &[]);
            let residual: f64 = poly.iter().rev().fold(0.0, |acc, c| acc * x + c);
            assert!(residual.abs() < 1e-7 * (1.0 + x.abs().powi(5)), "{src}: {x} gives {residual}");
        }
        // Five real roots (cyclic group): cosines, all five roots.
        let (text, reduced) = reduce_with(&rules, "solve(x^5 + x^4 - 4*x^3 - 3*x^2 + 3*x + 1, x)", &[]);
        assert!(reduced && !text.contains("root_of") && text.contains("cos("), "{text}");
        let items = list_items(&text);
        assert_eq!(items.len(), 5, "{text}");
        let mut values: Vec<f64> = items.iter().map(|i| eval(&rules, i, &[])).collect();
        values.sort_by(f64::total_cmp);
        for (v, want) in values.iter().zip([-1.918_985_947_228_994_7_f64, -1.309_721_467_890_570_2, -0.284_629_676_546_570_3, 0.830_830_026_003_772_9, 1.682_507_065_662_362_4]) {
            assert!((v - want).abs() < 1e-8, "{values:?}");
        }
        // Decomposable and self-reciprocal polynomials of higher degree.
        for (src, poly) in [
            ("solve(x^6 + 3*x^5 + 6*x^4 + 7*x^3 + 6*x^2 + 3*x - 1, x)", vec![-1.0, 3.0, 6.0, 7.0, 6.0, 3.0, 1.0]),
            ("solve(x^6 - 2*x^3 - 1, x)", vec![-1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 1.0]),
            ("solve(x^6 - x^4 - x^3 - x^2 + 1, x)", vec![1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0]),
        ] {
            let (text, reduced) = reduce_with(&rules, src, &[]);
            assert!(reduced, "{src}: {text}");
            let items = list_items(&text);
            assert!(!items.is_empty() && !text.contains("root_of"), "{src}: {text}");
            for item in &items {
                let x = eval(&rules, item, &[]);
                let residual: f64 = poly.iter().rev().fold(0.0, |acc, c| acc * x + c);
                assert!(residual.abs() < 1e-7 * (1.0 + x.abs().powi(6)), "{src}: {item} -> {x} gives {residual}");
            }
        }
        // Unsolvable: exact root objects, canonical and numerically evaluable.
        let (text, reduced) = reduce_with(&rules, "solve(x^5 - x - 1, x)", &[]);
        assert!(reduced && text == "list(root_of(0, -1, -1, 0, 0, 0, 1))", "{text}");
        let (text, _) = reduce_with(&rules, "solve(x^7 - 7*x + 3, x)", &[]);
        assert!(text.matches("root_of").count() == 3, "{text}");
        // A reducible polynomial: the object lives in the irreducible factor.
        let (text, _) = reduce_with(&rules, "root_of(1, 1, 0, -1, 0, 0, -1, 1)", &[]);
        assert_eq!(text, "root_of(0, -1, -1, 0, 0, 0, 1)");
        let (text, _) = reduce_with(&rules, "root_of(1, -2, 0, 1)", &[]);
        assert_eq!(text, "2^(1/2)");
        let (value, radius) = crate::rules::testing::numeric(&rules, "root_of(0, -1, -1, 0, 0, 0, 1)", &[], 1e-12);
        assert!((value - 1.167_303_978_261_418_7).abs() < 1e-12 && radius < 1e-12, "{value} {radius}");
        // The de Moivre family with parameters.
        let (text, reduced) = reduce_with(&rules, "solve(x^5 + 5*a*x^3 + 5*a^2*x + b, x)", &[]);
        assert!(reduced, "{text}");
        let items = list_items(&text);
        let mut checked = 0;
        for (a, b) in [(1.5, 2.0), (0.7, -3.0), (-1.0, 0.8), (-1.2, -1.5)] {
            for item in &items {
                let x = eval(&rules, item, &[("a", a), ("b", b)]);
                if x.is_finite() {
                    let residual = x.powi(5) + 5.0 * a * x.powi(3) + 5.0 * a * a * x + b;
                    assert!(residual.abs() < 1e-7 * (1.0 + x.abs().powi(5)), "{item} at a={a}, b={b}: {residual}");
                    checked += 1;
                }
            }
        }
        assert!(checked >= 8, "{checked}: {text}");
    }

    #[test]
    fn heuristic_sub_solvers() {
        let rules = crate::rules::standard();
        let roots = |src: &str| -> Vec<f64> {
            let (text, reduced) = crate::rules::testing::reduce_with(&rules, src, &[]);
            assert!(reduced, "{src} not solved: {text}");
            let inner = text.strip_prefix("list(").and_then(|t| t.strip_suffix(')')).unwrap_or("");
            if inner.is_empty() {
                return Vec::new();
            }
            let mut depth = 0;
            let mut parts = Vec::new();
            let mut current = String::new();
            for ch in inner.chars() {
                match ch {
                    | '(' => depth += 1,
                    | ')' => depth -= 1,
                    | ',' if depth == 0 => {
                        parts.push(std::mem::take(&mut current));
                        continue;
                    },
                    | _ => {},
                }
                current.push(ch);
            }
            parts.push(current);
            let mut values: Vec<f64> = parts.iter().map(|p| crate::rules::testing::eval(&rules, p.trim(), &[])).collect();
            values.sort_by(f64::total_cmp);
            values
        };
        let close = |got: &[f64], want: &[f64]| {
            got.len() == want.len() && got.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9)
        };
        // Cardano: one real root; the trigonometric form: three.
        let cube = roots("solve(x^3 - 2*x - 5, x)");
        assert!(close(&cube, &[2.094_551_481_542_327]), "{cube:?}");
        let three = roots("solve(x^3 - 3*x + 1, x)");
        assert!(close(&three, &[-1.879_385_241_571_816_6, 0.347_296_355_333_860_7, 1.532_088_886_237_956]), "{three:?}");
        // Biquadratic and Ferrari quartics, and a quintic numerically.
        let bi = roots("solve(x^4 - 10*x^2 + 1, x)");
        let s = (2.0_f64).sqrt() + (3.0_f64).sqrt();
        let d = (3.0_f64).sqrt() - (2.0_f64).sqrt();
        assert!(close(&bi, &[-s, -d, d, s]), "{bi:?}");
        let quintic = roots("solve(x^5 - x - 1, x)");
        assert!(close(&quintic, &[1.167_303_978_261_418_7]), "{quintic:?}");
        // Absolute values, trigonometric and Lambert-type equations.
        assert!(close(&roots("solve(abs(x - 1) - 2, x)"), &[-1.0, 3.0]));
        let trig = roots("solve(sin(x) + cos(x) - 1, x)");
        assert!(close(&trig, &[0.0, std::f64::consts::FRAC_PI_2]), "{trig:?}");
        let double = roots("solve(sin(2*x) - cos(x), x)");
        assert!(double.len() >= 3, "{double:?}");
        let w = roots("solve(x*exp(x) - 1, x)");
        assert!(close(&w, &[0.567_143_290_409_784]), "{w:?}");
        let lin_exp = roots("solve(x + exp(x), x)");
        assert!(close(&lin_exp, &[-0.567_143_290_409_784]), "{lin_exp:?}");
        let log = roots("solve(x*ln(x) - 1, x)");
        assert!(close(&log, &[1.763_222_834_351_896_7]), "{log:?}");
        // Radicals.
        assert!(close(&roots("solve(x - 3*x^(1/2) + 2, x)"), &[1.0, 4.0]));
        assert!(close(&roots("solve((x + 7)^(1/2) - x - 1, x)"), &[2.0]));
        // A non-polynomial system by elimination.
        let (system, reduced) = crate::rules::testing::reduce_with(&rules, "solve(list(y - exp(x), y - 2), list(x, y))", &[]);
        assert!(reduced && system.contains("ln(2)"), "{system}");
    }

    #[test]
    fn linear_and_quadratic() {
        assert_eq!(run("solve(2*x + 6 = 0, x)"), "list(-3)");
        assert_eq!(run("solve(x^2 = 4, x)"), "list(-2, 2)");
        assert_eq!(run("solve(x^2 - 5*x + 6, x)"), "list(2, 3)");
        assert_eq!(run("solve(x^2 + 1 = 0, x)"), "list()", "no real solutions");
        assert_eq!(run("solve(x^2 = 2, x)"), "list(-2^(1/2), 2^(1/2))");
        assert_eq!(run("solve(x^2 - 2*x - 1 = 0, x)"), "list(1 - 2^(1/2), 1 + 2^(1/2))");
        assert_eq!(run("solve(x^2 + x - 1 = 0, x)"), "list(-1/2 - 1/2*5^(1/2), 1/2*5^(1/2) - 1/2)");
    }

    #[test]
    fn higher_degree_by_factorisation() {
        assert_eq!(run("solve(x^3 - 6*x^2 + 11*x - 6 = 0, x)"), "list(1, 2, 3)");
        assert_eq!(run("solve(x^4 - 1 = 0, x)"), "list(-1, 1)");
        assert_eq!(run("solve(x^3 = 8, x)"), "list(2)");
        assert_eq!(run("solve(x^3 = 2, x)"), "list(2^(1/3))");
        assert_eq!(run("solve(x^4 = 5, x)"), "list(-5^(1/4), 5^(1/4))");
        assert_eq!(run("solve(x^3 - x = 0, x)"), "list(-1, 0, 1)");
        assert_eq!(run("solve((x - 1)^2 * (x + 2) = 0, x)"), "list(-2, 1)", "repeated roots once");
    }

    #[test]
    fn symbolic_coefficients() {
        assert_eq!(run("solve(a*x + b = 0, x)"), "list(-b/a)");
        assert_eq!(run("solve(x^2 = a, x)"), "list(-a^(1/2), a^(1/2))");
        let (text, reduced) = reduce_with(&[solve()], "solve(a*x^2 + b*x + c = 0, x)", &[]);
        assert!(reduced, "{text}");
        // Check the quadratic formula numerically for a = 1, b = -3, c = 2.
        let inner = text.strip_prefix("list(").and_then(|s| s.strip_suffix(')')).unwrap_or("");
        let middle = inner.find(", ").unwrap_or(0);
        let (first, second) = (&inner[..middle], &inner[middle + 2..]);
        let at = [("a", 1.0), ("b", -3.0), ("c", 2.0)];
        let mut values = [eval(&[solve()], first, &at), eval(&[solve()], second, &at)];
        values.sort_by(f64::total_cmp);
        assert!((values[0] - 1.0).abs() < 1e-12 && (values[1] - 2.0).abs() < 1e-12, "{values:?}");
    }

    #[test]
    fn rational_equations_drop_poles() {
        assert_eq!(run("solve((x^2 - 1)/(x - 1) = 0, x)"), "list(-1)", "x = 1 is a pole, not a root");
        assert_eq!(run("solve(1/x = 2, x)"), "list(1/2)");
        assert_eq!(run("solve(1/(x - 1) + 1/(x + 1) = 0, x)"), "list(0)");
    }

    #[test]
    fn transcendental_by_inversion() {
        assert_eq!(run("solve(exp(x) = 3, x)"), "list(ln(3))");
        assert_eq!(run("solve(ln(x) = 2, x)"), "list(exp(2))");
        assert_eq!(run("solve(2^x = 8, x)"), "list(ln(8)/ln(2))");
        assert_eq!(run("solve(exp(2*x) - 3*exp(x) + 2 = 0, x)"), "list(0, ln(2))");
        assert_eq!(run("solve(sin(x) = 0, x)"), "list(0)", "principal branch");
        assert_eq!(run("solve(sqrt(x) = 3, x)"), "list(9)");
        assert_eq!(run("solve(sqrt(x) = -3, x)"), "list()", "squaring introduced a false root");
        assert_eq!(run("solve(abs(x - 1) = 2, x)"), "list(-1, 3)");
        assert_eq!(run("solve((x - 1)^3 = 8, x)"), "list(3)");
        check_solutions("cosh(x) = 2", &run("solve(cosh(x) = 2, x)"), 2);
        check_solutions("tanh(x) = 1/2", &run("solve(tanh(x) = 1/2, x)"), 1);
        check_solutions("exp(x)^2 + exp(x) = 6", &run("solve(exp(x)^2 + exp(x) = 6, x)"), 1);
    }

    #[test]
    fn unsolvable_requests_stay_unreduced() {
        for src in ["solve(x + exp(x) = 0, x)", "solve(sin(x) = x, x)"] {
            let (text, reduced) = reduce_with(&[solve()], src, &[]);
            assert!(!reduced, "{src} unexpectedly gave {text}");
        }
        assert_eq!(run("solve(3 = 4, x)"), "list()");
    }

    #[test]
    fn linear_systems() {
        assert_eq!(run("solve(list(x + y = 3, x - y = 1), list(x, y))"), "list(list(2, 1))");
        assert_eq!(
            run("solve(list(2*x + 3*y - z = 1, x - y + z = 2, 3*x + y + 2*z = 7), list(x, y, z))"),
            "list(list(3/7, 6/7, 17/7))"
        );
        assert_eq!(run("solve(list(a*x + y = 1, x - y = 0), list(x, y))"), "list(list(1/(a + 1), 1/(a + 1)))");
        // A singular system has a line of solutions, given with a free unknown.
        assert_eq!(run("solve(list(x + y = 1, 2*x + 2*y = 2), list(x, y))"), "list(list(-(y - 1), y))");
    }

    #[test]
    fn polynomial_systems() {
        assert_eq!(
            run("solve(list(x^2 + y^2 = 5, x - y = 1), list(x, y))"),
            "list(list(-1, -2), list(2, 1))"
        );
        assert_eq!(run("solve(list(x*y = 1, x = 0), list(x, y))"), "list()", "inconsistent");
        assert_eq!(
            run("solve(list(x^2 + y^2 = 1, x = y), list(x, y))"),
            "list(list(-1/2*2^(1/2), -1/2*2^(1/2)), list(1/2*2^(1/2), 1/2*2^(1/2)))"
        );
    }

    #[test]
    fn numeric_phase() {
        // All real roots of an unsolvable quintic, via Sturm sequences.
        let mut g = Graph::new();
        let engine = crate::graph::Engine::install(&mut g, &[solve()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse("solve(x^5 - x + 1 = 0, x)").unwrap_or_else(|e| panic!("{e}"));
        let env = Env::numeric(1e-12);
        engine.run(&mut g, &[root], &env, &crate::graph::Saturate, &crate::graph::Budget::default());
        let answer = crate::graph::Extractor::new(&g, &[root], &crate::graph::ClosedForm).build(&mut g, root);
        let text = answer.map(|n| g.display(n)).unwrap_or_default();
        let value: f64 = text.trim_start_matches("list(").trim_end_matches(')').parse().unwrap_or(f64::NAN);
        assert!((value.powi(5) - value + 1.0).abs() < 1e-9, "{text}");

        // One root near a guess.
        let (value, error) = numeric(&[solve()], "nsolve(cos(x) = x, x, 1)", &[], 1e-12);
        assert!((value.cos() - value).abs() < 1e-10, "{value}");
        assert!(error <= 1e-12);
        let (value, _) = numeric(&[solve()], "nsolve(x^2 = a, x, 1)", &[("a", 2.0)], 1e-12);
        assert!((value - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn radicals_logarithms_self_powers_and_inequalities() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        assert_eq!(run("solve(sqrt(x + 3) = x - 3, x)"), "list(6)");
        assert_eq!(run("solve(sqrt(x) + sqrt(x - 5) = 5, x)"), "list(9)");
        assert_eq!(run("solve(ln(x) + ln(x - 1) = ln(6), x)"), "list(3)");
        assert_eq!(run("solve(x^x = 4, x)"), "list(exp(lambertw(ln(4))))");
        assert_eq!(run("solve(gt(x^3 - 6*x^2 + 11*x - 6, 0), x)"), "or(and(lt(1, x), lt(x, 2)), lt(3, x))");
        assert_eq!(run("solve(lt(x^2 - 4, 0), x)"), "and(lt(-2, x), lt(x, 2))");
        assert_eq!(run("solve(lt(abs(x), 2), x)"), "and(lt(-2, x), lt(x, 2))");
        assert_eq!(run("solve(ge((x - 1)/(x + 2), 0), x)"), "or(le(1, x), lt(x, -2))");
        assert_eq!(run("solve(le(x^2, 0), x)"), "x = 0");
        assert_eq!(run("solve(gt(x^2 + 1, 0), x)"), "true");
        assert_eq!(run("solve(lt(x^2 + 1, 0), x)"), "false");
    }

    /// The numeric-mode answer of `src` as text.
    fn numeric_text(src: &str) -> String {
        let rules = crate::rules::standard();
        let mut g = Graph::new();
        let engine = crate::graph::Engine::install(&mut g, &rules).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(&mut g, &[root], &Env::numeric(1e-12), &crate::graph::Saturate, &crate::graph::Budget::default());
        let answer = crate::graph::Extractor::new(&g, &[root], &crate::graph::ClosedForm).build(&mut g, root);
        answer.map(|n| g.display(n)).unwrap_or_default()
    }

    fn numbers(text: &str) -> Vec<f64> {
        text.replace("list", "")
            .split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == 'e' || c == 'E'))
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    #[test]
    fn general_solutions_carry_an_integer_parameter() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        assert_eq!(run("solve_general(sin(x) = 0, x)"), "list(n*pi)");
        assert_eq!(run("solve_general(cos(x) = 0, x)"), "list(n*pi + 1/2*pi)");
        assert_eq!(run("solve_general(tan(x) = 1, x)"), "list(n*pi + 1/4*pi)");
        assert_eq!(run("solve_general(sin(x) = 1/2, x)"), "list(2*n*pi + 1/6*pi, 2*n*pi + 5/6*pi)");
        assert_eq!(run("solve_general(sin(x)^2 = cos(x)^2, x)"), "list(1/2*n*pi + 1/4*pi)");
        assert_eq!(run("solve_general(sin(x) = a, x)"), "list(2*n*pi + asin(a), 2*n*pi - asin(a) + pi)");
        // Not periodic: the same as `solve`.
        assert_eq!(run("solve_general(exp(x) = 3, x)"), "list(ln(3))");
        // Through an inner function.
        let text = run("solve_general(sin(x^2) = 1/2, x)");
        assert!(text.contains("n*pi"), "{text}");
        // No real solution.
        assert_eq!(run("solve_general(sin(x) = 2, x)"), "list()");
    }

    #[test]
    fn parameters_factors_and_exponential_forms() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        assert_eq!(run("solve(exp(x)*(x - 1) = 0, x)"), "list(1)");
        assert_eq!(run("solve(4^x - 3*2^x - 4 = 0, x)"), "list(ln(4)/ln(2))");
        assert_eq!(run("solve(exp(x) + exp(-x) = 4, x)"), "list(ln(2 - 3^(1/2)), ln(2 + 3^(1/2)))");
        assert_eq!(run("solve(sinh(x) + cosh(x) = 2, x)"), "list(ln(2))");
        assert_eq!(run("solve(2^x = 3^(x - 1), x)"), "list(-ln(3)/(ln(2) - ln(3)))");
        assert_eq!(run("solve(x^3*ln(x) = 0, x)"), "list(1)");
        assert_eq!(run("solve(atan(2*x) + atan(3*x) = pi/4, x)"), "list(1/6)");
        assert_eq!(run("solve(asin(x) - acos(x) = 0, x)"), "list(1/2*2^(1/2))");
        assert_eq!(run("solve(asin(x) = acos(x), x)"), "list(1/2*2^(1/2))");
        let text = run("solve((x^2 - a)*(x - b) = 0, x)");
        assert!(text.contains('b') && text.contains("a^(1/2)"), "{text}");
        let text = run("solve(x^3 - (a + b + c)*x^2 + (a*b + a*c + b*c)*x - a*b*c = 0, x)");
        assert!(text.contains('a') && text.contains('b') && text.contains('c'), "{text}");
    }

    #[test]
    fn systems_with_free_parameters_and_conditions() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        assert_eq!(run("solve(list(x + y = 1), list(x, y))"), "list(list(-(y - 1), y))");
        assert_eq!(run("solve(list(x + y + z = 1, x - y = 0), list(x, y, z))"), "list(list(1/2 - 1/2*z, 1/2 - 1/2*z, z))");
        assert_eq!(run("solve(list(x + y = 1, x + y = 2), list(x, y))"), "list()");
        assert_eq!(run("solve(list(x^2 - y = 1, x + y = 3), list(x, y))").matches("list(").count(), 3);
        assert_eq!(run("solve(list(gt(x, 1), lt(x, 5)), x)"), "and(lt(1, x), lt(x, 5))");
        assert_eq!(run("solve(list(gt(x^2, 1), lt(x, 0)), x)"), "lt(x, -1)");
        assert_eq!(run("solve(list(gt(x, 3), lt(x, 1)), x)"), "false");
        assert_eq!(run("solve(gt(exp(x), 2), x)"), "lt(ln(2), x)");
        assert_eq!(run("solve(lt(ln(x), 1), x)"), "and(lt(0, x), lt(x, E))");
        assert_eq!(run("solve(list(le(x, 2), x^2 = 4), x)"), "or(x = -2, x = 2)");
        assert_eq!(run("solve(list(gt(x, 0), lt(x^2, 9), ne(x, 1)), x)").matches("and").count(), 2);
    }

    #[test]
    fn floor_and_ceiling_give_intervals() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        assert_eq!(run("solve(floor(x) = 2, x)"), "and(le(2, x), lt(x, 3))");
        assert_eq!(run("solve(ceil(x) = 2, x)"), "and(le(x, 2), lt(1, x))");
        assert_eq!(run("solve(floor(x) = 1/2, x)"), "false");
        assert_eq!(run("solve(lt(floor(x), 3), x)"), "lt(x, 3)");
        assert_eq!(run("solve(ge(floor(x), 3), x)"), "le(3, x)");
        assert_eq!(run("solve(floor(2*x + 1) = 3, x)"), "and(le(1, x), lt(x, 3/2))");
        assert_eq!(run("solve(list(floor(x) = 2, gt(x, 2.5)), x)"), "and(lt(2.5, x), lt(x, 3))");
    }

    #[test]
    fn numeric_fallbacks() {
        // sin(x) = x/3 has three real roots, 0 and +-2.2788626600758283.
        let roots = numbers(&numeric_text("solve(sin(x) = x/3, x)"));
        assert_eq!(roots.len(), 3, "{roots:?}");
        assert!(roots.iter().any(|r| r.abs() < 1e-9));
        assert!(roots.iter().any(|r| (r - 2.278_862_660_075_828_3).abs() < 1e-8));
        // The circle meets exp(x) + y = 1 twice.
        let tuples = numbers(&numeric_text("solve(list(x^2 + y^2 = 4, exp(x) + y = 1), list(x, y))"));
        assert_eq!(tuples.len(), 4, "{tuples:?}");
        for pair in tuples.chunks(2) {
            let (x, y) = (pair[0], pair[1]);
            assert!((x * x + y * y - 4.0).abs() < 1e-8 && (x.exp() + y - 1.0).abs() < 1e-8, "{pair:?}");
        }
        // Roots of even multiplicity: no sign change, found through f'.
        let text = numeric_text("solve(cosh(x) = 1 + x^2/2 + x^4/24, x)");
        let roots = numbers(&text);
        assert_eq!(roots.len(), 1, "{text}");
        assert!(roots[0].abs() < 5e-2, "{roots:?}");
        // Newton from a guess.
        let root = numbers(&numeric_text("nsolve(list(x^2 + y^2 = 4, x*y = 1), list(x, y), list(2, 0.5))"));
        assert_eq!(root.len(), 2);
        assert!((root[0] * root[1] - 1.0).abs() < 1e-9);
    }
}
