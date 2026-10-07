//! Real algebraic geometry and ideal-theoretic requests.
//!
//! | operator | meaning |
//! |---|---|
//! | `count_real_roots(p, x)` / `count_real_roots(p, x, a, b)` | distinct real roots (in `(a, b]`), by Sturm sequences |
//! | `real_root_intervals(p, x)` | disjoint exact isolating intervals `list(a, b)`, `a = b` for rational roots |
//! | `sturm_sequence(p, x)` | the Sturm chain |
//! | `resultant(p, q, x)`, `discriminant(p, x)` | by fraction-free (Bareiss) elimination of the Sylvester matrix |
//! | `leading_term(p, vars)`, `leading_coeff(p, vars)`, `leading_monomial(p, vars)` | with respect to an order (optional third argument `lex`, `grlex`, `grevlex`; default `lex`) |
//! | `normal_form(p, polys, vars)` | remainder of `p` modulo the Gröbner basis of `polys` (`grevlex`) |
//! | `ideal_member(p, polys, vars)` | whether `p` lies in the ideal |
//! | `simplify_with_relations(e, relations, vars)` | numerator and denominator of `e` reduced modulo the relations |
//! | `cad(polys, vars)` | sample points of the cells of a cylindrical algebraic decomposition of `R^n` sign-invariant for `polys` (any number of variables: projection by leading and all other coefficients, discriminants and pairwise resultants, then lifting over each sample) |
//! | `expand_trig(e)` | sines, cosines and tangents of sums and integer multiples expanded |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::groebner;
use super::groebner::GPoly;
use super::groebner::GroebnerLimits;
use super::groebner::Order;
use super::repr::Gens;
use super::repr::Limits;
use super::repr::Poly;
use super::repr::from_term;
use super::repr::to_term;
use super::univariate;
use super::univariate::QPoly;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    CountRealRoots,
    RootIntervals,
    Sturm,
    Resultant,
    Discriminant,
    LeadingTerm,
    LeadingCoeff,
    LeadingMonomial,
    NormalForm,
    IdealMember,
    WithRelations,
    Cad,
    ExpandTrig,
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let identity = |a: &[f64]| a.first().copied().unwrap_or(f64::NAN);
    for (name, arity, request) in [
        ("count_real_roots", Arity::Variadic, Request::CountRealRoots),
        ("real_root_intervals", Arity::Fixed(2), Request::RootIntervals),
        ("sturm_sequence", Arity::Fixed(2), Request::Sturm),
        ("resultant", Arity::Fixed(3), Request::Resultant),
        ("discriminant", Arity::Fixed(2), Request::Discriminant),
        ("leading_term", Arity::Variadic, Request::LeadingTerm),
        ("leading_coeff", Arity::Variadic, Request::LeadingCoeff),
        ("leading_monomial", Arity::Variadic, Request::LeadingMonomial),
        ("normal_form", Arity::Fixed(3), Request::NormalForm),
        ("ideal_member", Arity::Fixed(3), Request::IdealMember),
        ("simplify_with_relations", Arity::Fixed(3), Request::WithRelations),
        ("cad", Arity::Fixed(2), Request::Cad),
        ("expand_trig", Arity::Fixed(1), Request::ExpandTrig),
    ] {
        let mut desc = OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100);
        if request == Request::ExpandTrig {
            desc = desc.eval(identity);
        }
        let op = i.op(desc)?;
        i.kernel(&format!("poly/{name}"), Tier::Reduce, Algebra { op, request });
    }
    Ok(())
}

struct Algebra {
    op: OpId,
    request: Request,
}

impl Kernel for Algebra {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let graph = &mut *cx.graph;
        let result = match (self.request, args.as_slice()) {
            | (Request::CountRealRoots, &[p, x]) => {
                count_real_roots(graph, p, x, None).map(|n| (graph.int(n), false))
            },
            | (Request::CountRealRoots, &[p, x, a, b]) => {
                count_real_roots(graph, p, x, Some((a, b))).map(|n| (graph.int(n), false))
            },
            | (Request::RootIntervals, &[p, x]) => root_intervals_term(graph, p, x).map(|t| (t, true)),
            | (Request::Sturm, &[p, x]) => sturm_term(graph, p, x).map(|t| (t, true)),
            | (Request::Resultant, &[p, q, x]) => resultant_term(graph, p, q, x).map(|t| (t, false)),
            | (Request::Discriminant, &[p, x]) => discriminant_term(graph, p, x).map(|t| (t, false)),
            | (Request::LeadingTerm | Request::LeadingCoeff | Request::LeadingMonomial, &[p, vars, ref rest @ ..])
                if rest.len() <= 1 =>
            {
                leading(graph, self.request, p, vars, rest.first().copied()).map(|t| (t, false))
            },
            | (Request::NormalForm, &[p, polys, vars]) => normal_form_term(graph, p, polys, vars).map(|t| (t, false)),
            | (Request::IdealMember, &[p, polys, vars]) => normal_form_term(graph, p, polys, vars)
                .map(|t| graph.as_number(t).is_some_and(Number::is_zero))
                .map(|b| (graph.lit(Payload::Bool(b)), false)),
            | (Request::WithRelations, &[e, relations, vars]) => {
                with_relations(graph, e, relations, vars).map(|t| (t, false))
            },
            | (Request::Cad, &[polys, vars]) => cad(graph, polys, vars).map(|t| (t, true)),
            | (Request::ExpandTrig, &[e]) => {
                let term = super::best(graph, e);
                term.and_then(|t| expand_trig(graph, t, 0)).map(|t| {
                    let expanded = super::expand_form(graph, t).unwrap_or(t);
                    (expanded, true)
                })
            },
            | _ => None,
        };
        match result {
            | Some((t, true)) => Outcome::Pinned(t),
            | Some((t, false)) => Outcome::Equal(t),
            | None => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

// ---------------------------------------------------------------------------
// Sturm sequences

/// `p` as a univariate polynomial over `Q` in the symbol `x`.
fn univariate_q(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
) -> Option<QPoly> {
    let (numer, denom) = super::rational_function_in(graph, p, x)?;
    // Constant denominators only.
    match denom.as_slice() {
        | [c] if !c.is_zero() => Some(numer.iter().map(|a| a / c).collect()),
        | _ => None,
    }
}

/// The Sturm chain `p, p', -rem(p_{i-1}, p_i), ...`.
fn sturm_chain(p: &[BigRational]) -> Vec<QPoly> {
    let mut chain: Vec<QPoly> = vec![p.to_vec(), univariate::derivative(p)];
    while chain.last().is_some_and(|q| !q.is_empty()) {
        let n = chain.len();
        let Some((_, rem)) = univariate::divrem(&chain[n - 2], &chain[n - 1]) else {
            break;
        };
        if rem.is_empty() {
            break;
        }
        chain.push(rem.iter().map(|c| -c).collect());
    }
    chain.retain(|q| !q.is_empty());
    chain
}

/// Sign of `q` at `x`, or at `±∞` for `None` (`positive = true` for `+∞`).
fn sign_at(
    q: &[BigRational],
    at: Option<&BigRational>,
    positive: bool,
) -> i32 {
    let sign = |v: &BigRational| {
        if v.is_zero() {
            0
        } else if v.is_positive() {
            1
        } else {
            -1
        }
    };
    match at {
        | Some(x) => sign(&univariate::eval(q, x)),
        | None => {
            let lead = q.last().map_or(0, sign);
            let odd = q.len().is_multiple_of(2);
            if positive || !odd { lead } else { -lead }
        },
    }
}

fn variations(
    chain: &[QPoly],
    at: Option<&BigRational>,
    positive: bool,
) -> usize {
    let signs: Vec<i32> = chain.iter().map(|q| sign_at(q, at, positive)).filter(|&s| s != 0).collect();
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

/// Distinct real roots of `p` in `(a, b]` (`None` ends are infinite).
fn sturm_count(
    chain: &[QPoly],
    a: Option<&BigRational>,
    b: Option<&BigRational>,
) -> usize {
    variations(chain, a, false).saturating_sub(variations(chain, b, true))
}

/// An interval endpoint.
enum End {
    Finite(BigRational),
    Infinite,
}

impl End {
    const fn value(&self) -> Option<&BigRational> {
        match self {
            | Self::Finite(r) => Some(r),
            | Self::Infinite => None,
        }
    }
}

fn endpoint(
    graph: &Graph,
    e: NodeId,
) -> Option<End> {
    if let Some(r) = graph.number_of(e).and_then(Number::to_rational) {
        return Some(End::Finite(r));
    }
    let v = graph.eval(e, &crate::graph::Env::numeric(0.0))?;
    v.is_infinite().then_some(End::Infinite)
}

fn count_real_roots(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
    interval: Option<(NodeId, NodeId)>,
) -> Option<i64> {
    let q = univariate_q(graph, p, x)?;
    if q.len() < 2 {
        return None;
    }
    let chain = sturm_chain(&q);
    let (a, b) = match interval {
        | Some((a, b)) => (endpoint(graph, a)?, endpoint(graph, b)?),
        | None => (End::Infinite, End::Infinite),
    };
    i64::try_from(sturm_count(&chain, a.value(), b.value())).ok()
}

/// Exact isolating intervals of the distinct real roots, ascending.
pub(crate) fn isolate(q: &[BigRational]) -> Vec<(BigRational, BigRational)> {
    if q.len() < 2 {
        return Vec::new();
    }
    let chain = sturm_chain(q);
    // Cauchy bound.
    let lead = q.last().cloned().unwrap_or_else(BigRational::one);
    let bound = q.iter().rev().skip(1).map(|c| (c / &lead).abs()).fold(BigRational::zero(), |m, c| if c > m { c } else { m })
        + BigRational::one();
    let rational = univariate::rational_roots(q);
    let mut out = Vec::new();
    let mut stack = vec![(-bound.clone(), bound)];
    let mut work = 0_u32;
    while let Some((a, b)) = stack.pop() {
        work += 1;
        if work > 100_000 {
            break;
        }
        match sturm_count(&chain, Some(&a), Some(&b)) {
            | 0 => {},
            | 1 => match rational.iter().find(|r| **r > a && **r <= b) {
                | Some(r) => out.push((r.clone(), r.clone())),
                | None => out.push((a, b)),
            },
            | _ => {
                let mid = (&a + &b) / BigRational::from_integer(BigInt::from(2));
                stack.push((a, mid.clone()));
                stack.push((mid, b));
            },
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    out
}

/// Shrinks an isolating interval of `q` to width below `width`.
fn refine(
    q: &[BigRational],
    (mut a, mut b): (BigRational, BigRational),
    width: &BigRational,
) -> (BigRational, BigRational) {
    let two = BigRational::from_integer(BigInt::from(2));
    let chain = sturm_chain(q);
    for _ in 0..200 {
        if &b - &a <= *width {
            break;
        }
        let mid = (&a + &b) / &two;
        if univariate::eval(q, &mid).is_zero() {
            return (mid.clone(), mid);
        }
        if sturm_count(&chain, Some(&a), Some(&mid)) == 1 {
            b = mid;
        } else {
            a = mid;
        }
    }
    (a, b)
}

fn root_intervals_term(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let q = univariate_q(graph, p, x)?;
    let items: Vec<NodeId> = isolate(&q)
        .into_iter()
        .map(|(a, b)| {
            let (a, b) = (graph.num(Number::rat(a)), graph.num(Number::rat(b)));
            graph.node(core::LIST, &[a, b])
        })
        .collect();
    Some(graph.node(core::LIST, &items))
}

fn qpoly_term(
    graph: &mut Graph,
    q: &[BigRational],
    x: NodeId,
) -> NodeId {
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let numbers: Vec<Number> = q.iter().cloned().map(Number::rat).collect();
    to_term(graph, &gens, &Poly::from_univariate(gx, &numbers))
}

fn sturm_term(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let q = univariate_q(graph, p, x)?;
    let items: Vec<NodeId> = sturm_chain(&q).iter().map(|s| qpoly_term(graph, s, x)).collect();
    Some(graph.node(core::LIST, &items))
}

// ---------------------------------------------------------------------------
// Resultants

const CAP: usize = 200_000;

/// The determinant of a square matrix of polynomials by Bareiss'
/// fraction-free elimination.
fn bareiss(mut m: Vec<Vec<Poly>>) -> Option<Poly> {
    let n = m.len();
    if n == 0 {
        return Some(Poly::constant(Number::from(1)));
    }
    let mut sign = false;
    let mut previous = Poly::constant(Number::from(1));
    for k in 0..n.saturating_sub(1) {
        if m[k][k].is_zero() {
            let swap = (k + 1..n).find(|&r| !m[r][k].is_zero());
            match swap {
                | Some(r) => {
                    m.swap(k, r);
                    sign = !sign;
                },
                | None => return Some(Poly::zero()),
            }
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let a = m[i][j].mul(&m[k][k], CAP)?;
                let b = m[i][k].mul(&m[k][j], CAP)?;
                m[i][j] = a.sub(&b).div_exact(&previous, CAP)?;
            }
            m[i][k] = Poly::zero();
        }
        previous = m[k][k].clone();
    }
    let det = m[n - 1][n - 1].clone();
    Some(if sign { det.neg() } else { det })
}

/// `res_x(p, q)` from coefficient lists in `x` (ascending).
fn resultant(
    p: &[Poly],
    q: &[Poly],
) -> Option<Poly> {
    let (m, n) = (p.len().checked_sub(1)?, q.len().checked_sub(1)?);
    if m == 0 && n == 0 {
        return Some(Poly::constant(Number::from(1)));
    }
    let size = m + n;
    let mut rows = Vec::with_capacity(size);
    for shift in 0..n {
        let mut row = vec![Poly::zero(); size];
        for (k, c) in p.iter().rev().enumerate() {
            row[shift + k] = c.clone();
        }
        rows.push(row);
    }
    for shift in 0..m {
        let mut row = vec![Poly::zero(); size];
        for (k, c) in q.iter().rev().enumerate() {
            row[shift + k] = c.clone();
        }
        rows.push(row);
    }
    bareiss(rows)
}

/// Shared generators and the coefficient lists of `p` and `q` in `x`.
fn in_x(
    graph: &mut Graph,
    terms: &[NodeId],
    x: NodeId,
) -> Option<(Gens, u32, Vec<Vec<Poly>>)> {
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let mut out = Vec::new();
    for &t in terms {
        let t = super::best(graph, t)?;
        let p = from_term(graph, &mut gens, t, Limits::default())?;
        let mut coefficients = p.coefficients_in(gx);
        while coefficients.last().is_some_and(Poly::is_zero) {
            coefficients.pop();
        }
        out.push(coefficients);
    }
    Some((gens, gx, out))
}

fn resultant_term(
    graph: &mut Graph,
    p: NodeId,
    q: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let (gens, _, polys) = in_x(graph, &[p, q], x)?;
    let r = resultant(polys.first()?, polys.get(1)?)?;
    Some(to_term(graph, &gens, &r))
}

/// `disc_x(p) = (-1)^(n(n-1)/2) res(p, p') / lc(p)`.
fn discriminant(
    p: &[Poly],
    gx: u32,
) -> Option<Poly> {
    let n = p.len().checked_sub(1)?;
    if n == 0 {
        return None;
    }
    let mut derivative = Vec::with_capacity(n);
    for (k, c) in p.iter().enumerate().skip(1) {
        derivative.push(c.scale(&Number::from(i64::try_from(k).ok()?)));
    }
    let _ = gx;
    let r = resultant(p, &derivative)?;
    let lead = p.last()?;
    let d = r.div_exact(lead, CAP)?;
    Some(if (n * (n - 1) / 2) % 2 == 1 { d.neg() } else { d })
}

fn discriminant_term(
    graph: &mut Graph,
    p: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let (gens, gx, polys) = in_x(graph, &[p], x)?;
    let d = discriminant(polys.first()?, gx)?;
    Some(to_term(graph, &gens, &d))
}

// ---------------------------------------------------------------------------
// Monomial orders and ideals

fn order_of(
    graph: &Graph,
    order: Option<NodeId>,
) -> Option<Order> {
    let Some(order) = order else {
        return Some(Order::Lex);
    };
    match graph.interner().symbol_name(graph.symbol_of(order)?) {
        | "lex" => Some(Order::Lex),
        | "grlex" => Some(Order::GradedLex),
        | "grevlex" => Some(Order::GradedRevLex),
        | _ => None,
    }
}

fn leading(
    graph: &mut Graph,
    request: Request,
    p: NodeId,
    vars: NodeId,
    order: Option<NodeId>,
) -> Option<NodeId> {
    let order = order_of(graph, order)?;
    let list = graph.node(core::LIST, &[p]);
    let (polys, gens) = super::to_groebner(graph, list, vars, order)?;
    let poly = polys.first()?;
    let Some((mono, coeff)) = poly.leading() else {
        return Some(graph.int(0));
    };
    let (mono, coeff) = (mono.clone(), coeff.clone());
    let term = match request {
        | Request::LeadingCoeff => return Some(graph.num(Number::rat(coeff))),
        | Request::LeadingMonomial => GPoly::new(vec![(mono, BigRational::one())], order),
        | _ => GPoly::new(vec![(mono, coeff)], order),
    };
    Some(super::from_groebner(graph, &gens, &term))
}

fn normal_form_term(
    graph: &mut Graph,
    p: NodeId,
    polys: NodeId,
    vars: NodeId,
) -> Option<NodeId> {
    if graph.op(polys) != core::LIST {
        return None;
    }
    let order = Order::GradedRevLex;
    let mut all = vec![p];
    all.extend_from_slice(graph.children(polys));
    let list = graph.node(core::LIST, &all);
    let (converted, gens) = super::to_groebner(graph, list, vars, order)?;
    let (target, relations) = converted.split_first()?;
    let basis = groebner::groebner(relations, order, GroebnerLimits::default())?;
    let reduced = groebner::reduce(target, &basis, order);
    Some(super::from_groebner(graph, &gens, &reduced))
}

fn with_relations(
    graph: &mut Graph,
    e: NodeId,
    relations: NodeId,
    vars: NodeId,
) -> Option<NodeId> {
    let term = super::best(graph, e)?;
    let mut gens = Gens::default();
    let r = super::ratio(graph, &mut gens, term, Limits::default())?;
    let numer = to_term(graph, &gens, &r.numer);
    let denom = to_term(graph, &gens, &r.denom);
    let numer = normal_form_term(graph, numer, relations, vars)?;
    let denom = normal_form_term(graph, denom, relations, vars)?;
    if graph.as_number(denom).is_some_and(Number::is_zero) {
        return None;
    }
    let minus_one = graph.int(-1);
    let inverse = graph.node(core::POW, &[denom, minus_one]);
    Some(graph.node(core::MUL, &[numer, inverse]))
}

// ---------------------------------------------------------------------------
// Cylindrical algebraic decomposition

/// Sample points of a sign-invariant decomposition of `R` for the
/// univariate polynomials `qs`: every root, and a rational point in every
/// open interval between consecutive roots. Roots that are not rational
/// are given to within `1e-15`.
fn cad_line(qs: &[QPoly]) -> Vec<Number> {
    // Square-free product of the non-constant polynomials, so one Sturm
    // chain isolates every root of every polynomial.
    let mut product: QPoly = vec![BigRational::one()];
    for q in qs.iter().filter(|q| q.len() >= 2) {
        let g = univariate::gcd(&product, q);
        let (part, _) = univariate::divrem(q, &g).unwrap_or_default();
        product = univariate::mul(&product, &part);
    }
    let width = BigRational::new(BigInt::one(), BigInt::from(10_u64.pow(15)));
    let roots: Vec<(BigRational, BigRational)> =
        isolate(&product).into_iter().map(|iv| refine(&product, iv, &width)).collect();
    let one = BigRational::one();
    let two = BigRational::from_integer(BigInt::from(2));
    let mut samples = Vec::new();
    match (roots.first(), roots.last()) {
        | (Some(first), Some(last)) => {
            samples.push(Number::rat(first.0.floor() - &one));
            for (i, (a, b)) in roots.iter().enumerate() {
                samples.push(if a == b { Number::rat(a.clone()) } else { Number::Float(Number::rat((a + b) / &two).to_f64()) });
                if let Some(next) = roots.get(i + 1) {
                    samples.push(Number::rat((b + &next.0) / &two));
                }
            }
            samples.push(Number::rat(last.1.ceil() + &one));
        },
        | _ => samples.push(Number::from(0)),
    }
    samples
}

fn cad(
    graph: &mut Graph,
    polys: NodeId,
    vars: NodeId,
) -> Option<NodeId> {
    if graph.op(polys) != core::LIST || graph.op(vars) != core::LIST {
        return None;
    }
    let ps = graph.children(polys).to_vec();
    let vs = graph.children(vars).to_vec();
    match vs.as_slice() {
        | &[x] => {
            let qs = ps.iter().map(|&p| univariate_q(graph, p, x)).collect::<Option<Vec<_>>>()?;
            let items: Vec<NodeId> = cad_line(&qs).into_iter().map(|n| graph.num(n)).collect();
            Some(graph.node(core::LIST, &items))
        },
        | &[x, y] => cad_plane(graph, &ps, x, y),
        | _ => cad_space(graph, &ps, &vs),
    }
}

/// `p` with the generators in `point` replaced by rational values.
fn eval_partial(
    p: &Poly,
    point: &[(u32, BigRational)],
) -> Option<Poly> {
    let mut out = Poly::zero();
    for (mono, coeff) in p.terms() {
        let mut c = coeff.to_rational()?;
        let mut rest = Vec::new();
        for &(g, e) in mono {
            match point.iter().find(|(h, _)| *h == g) {
                | Some((_, v)) => c *= num_traits::pow(v.clone(), usize::try_from(e).ok()?),
                | None => rest.push((g, e)),
            }
        }
        out = out.add(&Poly::monomial(rest, Number::rat(c)));
    }
    Some(out)
}

/// The projection of `polys` along the generator `y`: every coefficient
/// in `y`, the discriminants and the pairwise resultants (`McCallum`'s
/// set, with the full coefficient sets for safety), without constants
/// and duplicates. Polynomials free of `y` pass through.
fn project(
    polys: &[Poly],
    y: u32,
) -> Option<Vec<Poly>> {
    let mut out: Vec<Poly> = Vec::new();
    let mut push = |q: Poly| {
        if q.as_constant().is_none() && !out.contains(&q) && !out.contains(&q.neg()) {
            out.push(q);
        }
    };
    let coefficient_lists: Vec<Vec<Poly>> = polys
        .iter()
        .map(|p| {
            let mut c = p.coefficients_in(y);
            while c.last().is_some_and(Poly::is_zero) {
                c.pop();
            }
            c
        })
        .collect();
    for (i, cs) in coefficient_lists.iter().enumerate() {
        if cs.len() <= 1 {
            if let Some(c) = cs.first() {
                push(c.clone());
            }
            continue;
        }
        for c in cs {
            push(c.clone());
        }
        if cs.len() >= 3 {
            push(discriminant(cs, y)?);
        }
        for other in coefficient_lists.iter().skip(i + 1) {
            if other.len() >= 2 {
                push(resultant(cs, other)?);
            }
        }
    }
    Some(out)
}

/// Sample points (one per cell) of a CAD of `R^n`, the generators in
/// `order` (the last one lifted last).
fn cad_samples(
    polys: &[Poly],
    order: &[u32],
) -> Option<Vec<Vec<Number>>> {
    let (&y, base) = order.split_last()?;
    let fibre_samples = |point: &[(u32, BigRational)]| -> Option<Vec<Number>> {
        let mut fibres = Vec::with_capacity(polys.len());
        for p in polys {
            let q = eval_partial(p, point)?;
            let mut coefficients: QPoly = if q.is_zero() {
                Vec::new()
            } else {
                q.univariate_in(y)?.iter().map(Number::to_rational).collect::<Option<_>>()?
            };
            while coefficients.last().is_some_and(Zero::is_zero) {
                coefficients.pop();
            }
            fibres.push(coefficients);
        }
        Some(cad_line(&fibres))
    };
    if base.is_empty() {
        return Some(fibre_samples(&[])?.into_iter().map(|v| vec![v]).collect());
    }
    let projection = project(polys, y)?;
    let mut cells = Vec::new();
    for sample in cad_samples(&projection, base)? {
        let point: Vec<(u32, BigRational)> = base
            .iter()
            .zip(&sample)
            .map(|(&g, v)| Some((g, v.to_rational().or_else(|| BigRational::from_float(v.to_f64()))?)))
            .collect::<Option<_>>()?;
        for value in fibre_samples(&point)? {
            let mut cell = sample.clone();
            cell.push(value);
            cells.push(cell);
        }
        if cells.len() > 100_000 {
            return None;
        }
    }
    Some(cells)
}

fn cad_space(
    graph: &mut Graph,
    ps: &[NodeId],
    vs: &[NodeId],
) -> Option<NodeId> {
    let mut gens = Gens::default();
    let order: Vec<u32> = vs.iter().map(|&v| gens.index(graph, v)).collect();
    let mut polys = Vec::with_capacity(ps.len());
    for &t in ps {
        let t = super::best(graph, t)?;
        let p = from_term(graph, &mut gens, t, Limits::default())?;
        if p.support().iter().any(|g| !order.contains(g)) {
            return None;
        }
        polys.push(p);
    }
    let cells = cad_samples(&polys, &order)?;
    let items: Vec<NodeId> = cells
        .into_iter()
        .map(|c| {
            let coordinates: Vec<NodeId> = c.into_iter().map(|v| graph.num(v)).collect();
            graph.node(core::LIST, &coordinates)
        })
        .collect();
    Some(graph.node(core::LIST, &items))
}

/// Collins projection onto `x` (leading coefficients and discriminants in
/// `y`, pairwise resultants), a decomposition of the `x` line, and lifting
/// over each sample.
fn cad_plane(
    graph: &mut Graph,
    ps: &[NodeId],
    x: NodeId,
    y: NodeId,
) -> Option<NodeId> {
    let (gens, gy, polys) = in_x(graph, ps, y)?;
    let gx = gens.find(graph, x);
    let as_x = |p: &Poly| -> Option<QPoly> {
        if p.is_zero() {
            return Some(Vec::new());
        }
        match gx {
            | Some(g) => p.univariate_in(g)?.iter().map(Number::to_rational).collect(),
            | None => Some(vec![p.as_constant()?.to_rational()?]),
        }
    };
    let mut projection: Vec<QPoly> = Vec::new();
    for (i, p) in polys.iter().enumerate() {
        if p.len() >= 2 {
            projection.push(as_x(p.last()?)?);
            if p.len() >= 3 {
                projection.push(as_x(&discriminant(p, gy)?)?);
            }
        }
        for q in polys.iter().skip(i + 1) {
            if p.len() >= 2 && q.len() >= 2 {
                projection.push(as_x(&resultant(p, q)?)?);
            }
        }
    }
    let mut cells = Vec::new();
    for xs in cad_line(&projection) {
        let at = xs.to_rational().or_else(|| BigRational::from_float(xs.to_f64()))?;
        // Each polynomial at x = xs, as a univariate polynomial in y.
        let mut fibres = Vec::new();
        for p in &polys {
            let mut q: QPoly = Vec::with_capacity(p.len());
            for c in p {
                let cq = as_x(c)?;
                q.push(univariate::eval(&cq, &at));
            }
            while q.last().is_some_and(Zero::is_zero) {
                q.pop();
            }
            fibres.push(q);
        }
        for ys in cad_line(&fibres) {
            let (a, b) = (graph.num(xs.clone()), graph.num(ys));
            cells.push(graph.node(core::LIST, &[a, b]));
        }
    }
    Some(graph.node(core::LIST, &cells))
}

// ---------------------------------------------------------------------------
// Trigonometric expansion

/// Trigonometric functions of sums and integer multiples expanded in the
/// concrete term `e` (for the equation solver).
pub(crate) fn expand_trig_term(
    graph: &mut Graph,
    e: NodeId,
) -> Option<NodeId> {
    expand_trig(graph, e, 0)
}

/// An isolating interval of `q` shrunk below `width`.
pub(crate) fn refine_interval(
    q: &[BigRational],
    interval: (BigRational, BigRational),
    width: &BigRational,
) -> (BigRational, BigRational) {
    refine(q, interval, width)
}

fn expand_trig(
    graph: &mut Graph,
    e: NodeId,
    depth: u32,
) -> Option<NodeId> {
    if depth > 40 {
        return None;
    }
    let ops = (graph.ops().lookup("sin"), graph.ops().lookup("cos"), graph.ops().lookup("tan"));
    let (Some(sin), Some(cos), Some(tan)) = ops else {
        return None;
    };
    let op = graph.op(e);
    let children = graph.children(e).to_vec();
    let expanded: Vec<NodeId> =
        children.iter().map(|&c| expand_trig(graph, c, depth + 1)).collect::<Option<Vec<_>>>()?;
    if op != sin && op != cos && op != tan {
        return Some(if expanded == children { e } else { graph.try_node(op, &expanded)? });
    }
    let &[arg] = expanded.as_slice() else {
        return Some(e);
    };
    // Split the argument into a + b, or n*x with integer n >= 2.
    let split = if graph.op(arg) == core::ADD {
        let parts = graph.children(arg).to_vec();
        let (first, rest) = parts.split_first()?;
        let rest = if rest.len() == 1 { rest[0] } else { graph.node(core::ADD, rest) };
        Some((*first, rest))
    } else if graph.op(arg) == core::MUL {
        let parts = graph.children(arg).to_vec();
        let n = parts.iter().find_map(|&p| graph.number_of(p).and_then(Number::to_i64));
        match n {
            | Some(n) if (2..=24).contains(&n) => {
                let others: Vec<NodeId> = parts.iter().copied().filter(|&p| graph.number_of(p).is_none()).collect();
                let rest_number: Vec<NodeId> = parts
                    .iter()
                    .copied()
                    .filter(|&p| graph.number_of(p).is_some_and(|v| v.to_i64() != Some(n)))
                    .collect();
                let mut unit_parts = others;
                unit_parts.extend(rest_number);
                let unit = if unit_parts.len() == 1 { unit_parts[0] } else { graph.node(core::MUL, &unit_parts) };
                let m = graph.int(n - 1);
                let rest = graph.node(core::MUL, &[m, unit]);
                Some((rest, unit))
            },
            | _ => None,
        }
    } else {
        None
    };
    let Some((a, b)) = split else {
        return Some(graph.node(op, &[arg]));
    };
    let parts = |graph: &mut Graph, t: NodeId| -> Option<(NodeId, NodeId)> {
        let s = graph.node(sin, &[t]);
        let c = graph.node(cos, &[t]);
        Some((expand_trig(graph, s, depth + 1)?, expand_trig(graph, c, depth + 1)?))
    };
    let (sa, ca) = parts(graph, a)?;
    let (sb, cb) = parts(graph, b)?;
    let minus_one = graph.int(-1);
    let sin_sum = {
        let x = graph.node(core::MUL, &[sa, cb]);
        let y = graph.node(core::MUL, &[ca, sb]);
        graph.node(core::ADD, &[x, y])
    };
    let cos_sum = {
        let x = graph.node(core::MUL, &[ca, cb]);
        let y = graph.node(core::MUL, &[minus_one, sa, sb]);
        graph.node(core::ADD, &[x, y])
    };
    Some(if op == sin {
        sin_sum
    } else if op == cos {
        cos_sum
    } else {
        let inverse = graph.node(core::POW, &[cos_sum, minus_one]);
        graph.node(core::MUL, &[sin_sum, inverse])
    })
}

#[cfg(test)]
mod tests {
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&crate::rules::standard(), src)
    }

    #[test]
    fn sturm_counts_and_isolation() {
        assert_eq!(run("count_real_roots(x^3 - 2*x, x)"), "3");
        assert_eq!(run("count_real_roots(x^4 + 1, x)"), "0");
        assert_eq!(run("count_real_roots((x - 1)^2*(x + 2), x)"), "2");
        assert_eq!(run("count_real_roots(x^3 - 2*x, x, 0, 10)"), "1");
        assert_eq!(run("count_real_roots(x^5 - x - 1, x, -oo, oo)"), "1");
        assert_eq!(run("real_root_intervals(x^2 - 4, x)"), "list(list(-2, -2), list(2, 2))");
        let intervals = run("real_root_intervals(x^2 - 2, x)");
        assert!(intervals.starts_with("list(list(-"), "{intervals}");
        assert_eq!(run("sturm_sequence(x^3 - 2*x, x)"), "list(x^3 - 2*x, 3*x^2 - 2, 4/3*x, 2)");
    }

    #[test]
    fn resultants_and_discriminants() {
        assert_eq!(run("discriminant(a*x^2 + b*x + c, x)"), "b^2 - 4*a*c");
        assert_eq!(run("resultant(x^2 - 1, x - 1, x)"), "0");
        assert_eq!(run("resultant(x^2 + y^2 - 1, x - y, x)"), "2*y^2 - 1");
        assert_eq!(run("discriminant(x^3 + p*x + q, x)"), "-4*p^3 - 27*q^2");
    }

    #[test]
    fn leading_terms_and_ideals() {
        assert_eq!(run("leading_term(3*x*y^2 + 2*x^2 + y, list(x, y))"), "2*x^2");
        assert_eq!(run("leading_term(3*x*y^2 + 2*x^2 + y, list(x, y), grlex)"), "3*x*y^2");
        assert_eq!(run("leading_coeff(3*x*y^2 + 2*x^2 + y, list(x, y), grlex)"), "3");
        assert_eq!(run("leading_monomial(3*x*y^2 + 2*x^2 + y, list(x, y), grlex)"), "x*y^2");
        assert_eq!(run("ideal_member(x^3 - x*y^2, list(x^2 - y^2), list(x, y))"), "true");
        assert_eq!(run("ideal_member(x + y, list(x^2 - y^2), list(x, y))"), "false");
        assert_eq!(run("normal_form(x^2, list(x^2 + y^2 - 1), list(x, y))"), "1 - y^2");
        assert_eq!(run("simplify_with_relations(x^2 + y^2 + 3, list(x^2 + y^2 - 1), list(x, y))"), "4");
        assert_eq!(
            run("simplify_with_relations((x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4), list(x^2 + y^2 - 1), list(x, y))"),
            "1"
        );
    }

    #[test]
    fn cylindrical_decomposition() {
        // Roots ±1 and the three open intervals.
        assert_eq!(run("cad(list(x^2 - 1), list(x))"), "list(-2, -1, 0, 1, 2)");
        // The unit circle: 3 cylinders over the open intervals and 2 over
        // the sections x = ±1, each split by the circle.
        let (cells, reduced) = reduce_with(&crate::rules::standard(), "cad(list(x^2 + y^2 - 1), list(x, y))", &[]);
        assert!(reduced, "{cells}");
        assert_eq!(cells.matches("list(").count() - 1, 1 + 1 + 3 + 1 + 5 + 1 + 3 + 1 + 1 - 4, "{cells}");
    }

    #[test]
    fn cylindrical_decomposition_in_three_variables() {
        let rules = crate::rules::standard();
        // The unit sphere: every sample is inside, on or outside, and all
        // three occur; the cells over the poles are found.
        let (cells, reduced) = reduce_with(&rules, "cad(list(x^2 + y^2 + z^2 - 1), list(x, y, z))", &[]);
        assert!(reduced, "{cells}");
        let points: Vec<Vec<f64>> = cells
            .trim_start_matches("list(")
            .trim_end_matches(')')
            .split("), list(")
            .map(|c| c.trim_start_matches("list(").trim_end_matches(')').split(", ").map(|v| v.parse::<f64>().unwrap_or_else(|_| crate::rules::testing::eval(&rules, v, &[]))).collect())
            .collect();
        assert!(points.iter().all(|p| p.len() == 3), "{cells}");
        let signs: Vec<i32> = points
            .iter()
            .map(|p| {
                let v = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0;
                if v.abs() < 1e-9 { 0 } else if v < 0.0 { -1 } else { 1 }
            })
            .collect();
        assert!(signs.contains(&-1) && signs.contains(&0) && signs.contains(&1), "{cells}");
        assert!(points.iter().any(|p| (p[0] - 1.0).abs() < 1e-12 && p[1].abs() < 1e-12 && p[2].abs() < 1e-12), "{cells}");
        // Two planes meeting along a line: every sign combination of
        // (x - z, y + z) is realised.
        let (cells, reduced) = reduce_with(&rules, "cad(list(x - z, y + z), list(x, y, z))", &[]);
        assert!(reduced, "{cells}");
        assert!(cells.matches("list(").count() >= 9, "{cells}");
    }

    #[test]
    fn trig_expansion() {
        assert_eq!(run("expand_trig(sin(a + b))"), "cos(a)*sin(b) + cos(b)*sin(a)");
        assert_eq!(run("expand_trig(cos(2*x))"), "cos(x)^2 - sin(x)^2");
        let triple = run("expand_trig(sin(3*x))");
        let rules = crate::rules::standard();
        let at = crate::rules::testing::eval(&rules, &triple, &[("x", 0.7)]);
        assert!((at - (2.1_f64).sin()).abs() < 1e-12, "{triple}");
    }
}
