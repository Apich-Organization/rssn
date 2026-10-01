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
//! Inverse functions are taken on their principal branch, so periodic
//! equations yield one representative solution per branch of the inverse,
//! not the full solution set.

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
use super::poly::from_groebner;
use super::poly::groebner::groebner;
use super::poly::groebner::GroebnerLimits;
use super::poly::groebner::Order;
use super::poly::poly;
use super::poly::ratio;
use super::poly::repr::from_term;
use super::poly::repr::to_term;
use super::poly::repr::Gens;
use super::poly::repr::Limits;
use super::poly::repr::Poly;
use super::poly::to_groebner;
use super::poly::univariate;

/// Operator attribute: the solutions `u` of `op(u) = ?a`, one pattern per
/// branch that is returned.
#[derive(Clone, Debug)]
pub struct Inverse(pub Vec<Pat>);

/// The equation-solving rule set.
#[must_use]
pub fn solve() -> RuleSet {
    RuleSet::new("solve", install).needs(poly()).needs(elementary())
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
    let nsolve_op = i.op(OpDescriptor::new("nsolve", Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("solve/symbolic", Tier::Reduce, Symbolic { solve: solve_op });
    i.kernel("solve/polynomial-numeric", Tier::Reduce, PolynomialNumeric { solve: solve_op });
    i.kernel("solve/nsolve", Tier::Reduce, RootNear { nsolve: nsolve_op });
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
    if !graph.depends_on(graph.find(term), symbol) {
        // `0 = 0` is solved by everything, which a list cannot express.
        return if graph.number_of(term).is_some_and(Number::is_zero) { None } else { Some(Vec::new()) };
    }
    // Write exp(n*u) as exp(u)^n so that `exp(2*x) - 3*exp(x) + 2` is a
    // polynomial in the single generator exp(x).
    let term = split_integer_multiples(graph, term);
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let fraction = ratio(graph, &mut gens, term, Limits::default())?;
    let dependent: Vec<u32> = fraction
        .numer
        .support()
        .into_iter()
        .filter(|&g| gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)))
        .collect();
    let candidates = match dependent.as_slice() {
        | [] => Vec::new(),
        | &[g] if g == gx => polynomial_roots(graph, &gens, &fraction.numer, gx)?,
        | &[g] => {
            // A polynomial in one function of x: solve for the function,
            // then invert it.
            let inner = gens.node(g)?;
            let mut out = Vec::new();
            for value in polynomial_roots(graph, &gens, &fraction.numer, g)? {
                out.extend(invert(graph, inner, value, x, depth + 1)?);
            }
            out
        },
        | _ => return None,
    };
    // Discard candidates at which the original expression is not zero
    // (poles of a cancelled denominator, roots introduced by squaring),
    // wherever that can be checked numerically.
    let mut solutions: Vec<NodeId> = Vec::new();
    for candidate in candidates {
        let substituted = graph.substitute(term, x, candidate);
        let keep = match graph.eval(substituted, &Env::numeric(0.0)) {
            | Some(v) => v.abs() <= 1e-7,
            | None => true,
        };
        if keep && !solutions.iter().any(|&s| graph.same(s, candidate) || s == candidate) {
            solutions.push(candidate);
        }
    }
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

/// Rewrites `exp(n*u)` as `exp(u)^n` and `b^(n*u)` as `(b^u)^n` for
/// integer literals `n`, throughout the concrete term `term`.
fn split_integer_multiples(
    graph: &mut Graph,
    term: NodeId,
) -> NodeId {
    let exp = graph.ops().lookup("exp");
    let mut done: std::collections::HashMap<NodeId, NodeId> = std::collections::HashMap::new();
    let mut stack = vec![(term, false)];
    while let Some((node, expanded)) = stack.pop() {
        if done.contains_key(&node) {
            continue;
        }
        let children = graph.children(node).to_vec();
        if children.is_empty() {
            done.insert(node, node);
            continue;
        }
        if !expanded {
            stack.push((node, true));
            stack.extend(children.iter().filter(|c| !done.contains_key(c)).map(|&c| (c, false)));
            continue;
        }
        let rebuilt: Vec<NodeId> = children.iter().map(|c| done.get(c).copied().unwrap_or(*c)).collect();
        let op = graph.op(node);
        // The argument that may carry the integer multiple: the only
        // argument of exp, the exponent of a power.
        let exponent_index = if Some(op) == exp {
            Some(0)
        } else if op == core::POW {
            Some(1)
        } else {
            None
        };
        let split = exponent_index.and_then(|index| {
            let argument = *rebuilt.get(index)?;
            if graph.op(argument) != core::MUL {
                return None;
            }
            let factors = graph.children(argument).to_vec();
            let position = factors.iter().position(|&f| {
                graph.as_number(f).and_then(Number::to_i64).is_some_and(|n| n.abs() > 1)
            })?;
            let multiple = *factors.get(position)?;
            let rest: Vec<NodeId> =
                factors.iter().enumerate().filter(|&(i, _)| i != position).map(|(_, &f)| f).collect();
            Some((index, multiple, product(graph, &rest)))
        });
        let new = match split {
            | Some((index, multiple, rest)) => {
                let mut inner_args = rebuilt.clone();
                if let Some(slot) = inner_args.get_mut(index) {
                    *slot = rest;
                }
                let inner = graph.node(op, &inner_args);
                graph.node(core::POW, &[inner, multiple])
            },
            | None if rebuilt == children => node,
            | None => graph.try_node(op, &rebuilt).unwrap_or(node),
        };
        done.insert(node, new);
    }
    done.get(&term).copied().unwrap_or(term)
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
        let branches = graph.ops().attr::<Inverse>(op)?.0.clone();
        for branch in branches {
            targets.push((arg, branch.instantiate(graph, &[value])?));
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
    // Symbolic coefficients: formulas for degree one and two.
    let terms: Vec<NodeId> = coefficients.iter().map(|c| to_term(graph, gens, c)).collect();
    match *terms.as_slice() {
        | [c0, c1] => {
            let minus_one = graph.int(-1);
            let inverse = reciprocal(graph, c1);
            roots.push(product(graph, &[minus_one, c0, inverse]));
        },
        | [_, _, _] => {
            let (c, b, a) = (coefficients.first()?, coefficients.get(1)?, coefficients.get(2)?);
            let cap = Limits::default().terms;
            let four_ac = a.mul(c, cap)?.scale(&Number::from(4));
            let discriminant = to_term(graph, gens, &b.mul(b, cap)?.sub(&four_ac));
            let half = graph.num(Number::fraction(1, 2)?);
            let root = graph.node(core::POW, &[discriminant, half]);
            let neg_b = to_term(graph, gens, &b.neg());
            let two_a = to_term(graph, gens, &a.scale(&Number::from(2)));
            let inverse = reciprocal(graph, two_a);
            for sign in [-1, 1] {
                let s = graph.int(sign);
                let signed_root = product(graph, &[s, root]);
                let numerator = graph.node(core::ADD, &[neg_b, signed_root]);
                roots.push(product(graph, &[numerator, inverse]));
            }
        },
        | _ => return None,
    }
    Some(roots)
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
                    return None;
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

/// Solves a polynomial system with rational coefficients through a
/// lexicographic Gröbner basis and back-substitution.
fn solve_polynomial_system(
    graph: &mut Graph,
    equations: NodeId,
    unknowns: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let (generators, gens) = to_groebner(graph, equations, unknowns, Order::Lex)?;
    let basis = groebner(&generators, Order::Lex, GroebnerLimits::default())?;
    let terms: Vec<NodeId> = basis.iter().map(|g| from_groebner(graph, &gens, g)).collect();
    let variables = graph.children(unknowns).to_vec();
    back_substitute(graph, &terms, &variables)
}

/// Solves a triangular system for the last variable first.
fn back_substitute(
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
}

impl Kernel for Symbolic {
    fn ops(&self) -> Vec<OpId> {
        vec![self.solve]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[equation, unknown] = graph.children(node) else {
            return Outcome::Pass;
        };
        let solutions = if graph.op(unknown) == core::LIST {
            if graph.op(equation) != core::LIST {
                return Outcome::Pass;
            }
            let equations = graph.children(equation).to_vec();
            let unknowns = graph.children(unknown).to_vec();
            let tuples = match solve_linear(graph, &equations, &unknowns) {
                | Some(single) => Some(vec![single]),
                | None => solve_polynomial_system(graph, equation, unknown),
            };
            tuples.map(|tuples| tuples.iter().map(|t| graph.node(core::LIST, t)).collect::<Vec<_>>())
        } else {
            let expr = as_expression(graph, equation);
            solve_for(graph, expr, unknown, 0)
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
            | Ok(root) => Outcome::Approx(Ball { mid: root, rad: tolerance }),
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
            // Replace the unknown `x`, but not the x in `exp`.
            let template = format!("({lhs}) - ({rhs})");
            let chars: Vec<char> = template.chars().collect();
            let mut residual = String::new();
            for (i, &ch) in chars.iter().enumerate() {
                let alone = ch == 'x'
                    && !(i > 0 && chars[i - 1].is_alphanumeric())
                    && !chars.get(i + 1).is_some_and(|c| c.is_alphanumeric());
                if alone {
                    residual.push_str(&format!("({item})"));
                } else {
                    residual.push(ch);
                }
            }
            let value = eval(&[solve()], &residual, &[]);
            assert!(value.abs() < 1e-9, "{item} does not solve {equation}: residual {value}");
        }
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
        for src in ["solve(x + exp(x) = 0, x)", "solve(x^5 - x + 1 = 0, x)", "solve(sin(x) = x, x)"] {
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
        let (text, reduced) = reduce_with(&[solve()], "solve(list(x + y = 1, 2*x + 2*y = 2), list(x, y))", &[]);
        assert!(!reduced, "a singular system has no unique solution: {text}");
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
}
