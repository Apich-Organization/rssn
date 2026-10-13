//! `root_of(p, x, k)`: the `k`-th real root (counted from 0 in ascending
//! order) of the polynomial `p` in `x`, as an exact object.
//!
//! The operator binds `x` in `p`. Its kernel keeps the value canonical so
//! that equal numbers have equal representations:
//!
//! * `p` is factored over `Q`; the `k`-th real root of `p` is located in
//!   one irreducible factor by exact isolation;
//! * a root of a factor of degree one is a rational number and one of
//!   degree two is written with a square root;
//! * a root of a factor of degree three or more becomes
//!   `root_of(f, x, j)` with `f` the primitive irreducible factor with a
//!   positive leading coefficient and `j` the index of the root among the
//!   real roots of `f` (so `root_of(x^5 - x - 1, x, 0)` is the one real
//!   root, and `root_of((x - 1)(x^5 - x - 1), x, 1)` is the same object);
//! * an index that is not the index of a real root leaves the term as it
//!   is.
//!
//! In numeric mode the value is the midpoint of an isolating interval
//! refined to double precision, with the interval's width as the radius
//! of the enclosure (a rigorous bound: the interval contains exactly one
//! root of the polynomial).

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::algebra::isolate;
use super::algebra::refine_interval;
use super::repr::to_term;
use super::repr::Gens;
use super::repr::Poly;
use super::univariate;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new("root_of", Arity::Fixed(3)).binder(1, 0b1).cost(6))?;
    i.kernel("poly/root_of", Tier::Reduce, RootOf { op });
    Ok(())
}

struct RootOf {
    op: OpId,
}

/// The real roots of the primitive irreducible factors of `q`, as
/// `(factor, index in the factor, isolating interval)` sorted by value.
fn real_roots_by_factor(q: &[BigRational]) -> Vec<(Vec<BigInt>, usize, (BigRational, BigRational))> {
    let (_, factors) = univariate::factor(q);
    let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(40));
    let mut all = Vec::new();
    for (factor, _) in factors {
        let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
        for (j, interval) in isolate(&rational).into_iter().enumerate() {
            let refined = refine_interval(&rational, interval, &width);
            all.push((factor.clone(), j, refined));
        }
    }
    all.sort_by(|a, b| (&a.2.0 + &a.2.1).cmp(&(&b.2.0 + &b.2.1)));
    all
}

/// `f` with a positive leading coefficient.
fn positive(mut f: Vec<BigInt>) -> Vec<BigInt> {
    if f.last().is_some_and(Signed::is_negative) {
        for c in &mut f {
            *c = -c.clone();
        }
    }
    f
}

/// The canonical term `root_of(f, x, j)` for a primitive irreducible `f`
/// of degree at least three.
pub(crate) fn make(
    graph: &mut Graph,
    factor: &[BigInt],
    x: NodeId,
    index: usize,
) -> Option<NodeId> {
    let op = graph.ops().lookup("root_of")?;
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let numbers: Vec<Number> = positive(factor.to_vec()).into_iter().map(Number::Int).collect();
    let p = to_term(graph, &gens, &Poly::from_univariate(gx, &numbers));
    let k = graph.int(i64::try_from(index).ok()?);
    Some(graph.node(op, &[p, x, k]))
}

/// The roots of a polynomial with rational coefficients (ascending) that
/// have no better exact form, as canonical `root_of` terms, with their
/// double-precision values. Degree-one and degree-two factors are not
/// handled here.
pub(crate) fn real_roots_as_objects(
    graph: &mut Graph,
    factor: &[BigInt],
) -> Option<Vec<(f64, NodeId)>> {
    let x = graph.sym("x");
    let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
    let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(17));
    let mut out = Vec::new();
    for (j, interval) in isolate(&rational).into_iter().enumerate() {
        let (a, b) = refine_interval(&rational, interval, &width);
        let mid = ((a + b) / BigRational::from_integer(BigInt::from(2))).to_f64()?;
        out.push((mid, make(graph, factor, x, j)?));
    }
    Some(out)
}

impl Kernel for RootOf {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[p, x, k] = cx.graph.children(node) else { return Outcome::Pass };
        let numeric = cx.env.numeric;
        let graph = &mut *cx.graph;
        let Some(k) = graph.number_of(k).and_then(Number::to_i64).and_then(|k| usize::try_from(k).ok()) else {
            return Outcome::Pass;
        };
        let Some((numer, denom)) = super::rational_function_in(graph, p, x) else { return Outcome::Pass };
        let [c] = denom.as_slice() else { return Outcome::Pass };
        if c.is_zero() {
            return Outcome::Pass;
        }
        let q: Vec<BigRational> = numer.iter().map(|a| a / c).collect();
        if q.len() < 2 {
            return Outcome::Pass;
        }
        let roots = real_roots_by_factor(&q);
        let Some((factor, j, interval)) = roots.get(k).cloned() else { return Outcome::Pass };
        let factor = positive(factor);
        match factor.as_slice() {
            | [c0, c1] => {
                let value = BigRational::new(-c0.clone(), c1.clone());
                Outcome::Equal(graph.num(Number::rat(value)))
            },
            | [c, b, a] => {
                let discriminant = b * b - BigInt::from(4) * a * c;
                let (outside, inside) = square_part(&discriminant);
                let denominator = BigInt::from(2) * a;
                let sign = if j == 0 { -1 } else { 1 };
                let centre = graph.num(Number::rat(BigRational::new(-b.clone(), denominator.clone())));
                let scale = BigRational::new(outside * BigInt::from(sign), denominator);
                let (scale_node, radicand) = (graph.num(Number::rat(scale)), graph.num(Number::Int(inside)));
                let Some(half) = Number::fraction(1, 2) else { return Outcome::Pass };
                let half = graph.num(half);
                let radical = graph.node(core::POW, &[radicand, half]);
                let offset = graph.node(core::MUL, &[scale_node, radical]);
                Outcome::Equal(graph.node(core::ADD, &[centre, offset]))
            },
            | _ => {
                let Some(canonical) = make(graph, &factor, x, j) else { return Outcome::Pass };
                if graph.find(canonical) != graph.find(node) && !graph.same(canonical, node) {
                    return Outcome::Equal(canonical);
                }
                if !numeric {
                    return Outcome::Pass;
                }
                let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
                let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(17));
                let (a, b) = refine_interval(&rational, interval, &width);
                let (lo, hi) = (a.to_f64().unwrap_or(f64::NAN), b.to_f64().unwrap_or(f64::NAN));
                let mid = 0.5 * lo + 0.5 * hi;
                Outcome::Approx(Ball { mid, rad: (hi - mid).max(mid - lo) })
            },
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `value = outside^2 * inside` with `inside` free of small square factors.
fn square_part(value: &BigInt) -> (BigInt, BigInt) {
    let mut outside = BigInt::one();
    let mut rest = value.clone();
    let mut p = BigInt::from(2);
    while &p * &p <= rest && p < BigInt::from(100_000) {
        while rest.is_multiple_of(&(&p * &p)) {
            rest /= &p * &p;
            outside *= &p;
        }
        p += 1;
    }
    (outside, rest)
}
