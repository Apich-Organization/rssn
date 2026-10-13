//! `root_of(k, c0, c1, ..., cn)`: the `k`-th real root (counted from 0 in
//! ascending order) of the polynomial `c0 + c1 x + ... + cn x^n`, as an
//! exact object. The coefficients are the arguments, so the object has
//! numeric semantics wherever the coefficients have (they may be
//! parameters); the kernel keeps the value canonical so that equal numbers
//! have equal representations:
//!
//! * `p` is factored over `Q`; the `k`-th real root of `p` is located in
//!   one irreducible factor by exact isolation;
//! * a root of a factor of degree one is a rational number and one of
//!   degree two is written with a square root;
//! * a root of a factor of degree three or more becomes
//!   `root_of(j, f0, ..., fm)` with `f` the primitive irreducible factor with
//!   a positive leading coefficient and `j` the index of the root among the
//!   real roots of `f` (so `root_of(0, -1, -1, 0, 0, 0, 1)` is the one real
//!   root of `x^5 - x - 1`, and the index 1 root of `(x - 1)(x^5 - x - 1)` is
//!   the same object);
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
    let op = i.op(OpDescriptor::new("root_of", Arity::Variadic).cost(40).eval(evaluate))?;
    i.kernel("poly/root_of", Tier::Normalize, RootOf { op });
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

/// The canonical term `root_of(j, f0, ..., fm)` for a primitive irreducible
/// `f` of degree at least three.
pub(crate) fn make(
    graph: &mut Graph,
    factor: &[BigInt],
    index: usize,
) -> Option<NodeId> {
    let op = graph.ops().lookup("root_of")?;
    let mut args = vec![graph.int(i64::try_from(index).ok()?)];
    for c in positive(factor.to_vec()) {
        args.push(graph.num(Number::Int(c)));
    }
    Some(graph.node(op, &args))
}

/// The real roots of a polynomial with integer coefficients (ascending)
/// as canonical `root_of` terms, with their double-precision values.
/// Degree-one and degree-two factors are not handled here.
pub(crate) fn real_roots_as_objects(
    graph: &mut Graph,
    factor: &[BigInt],
) -> Option<Vec<(f64, NodeId)>> {
    let rational: Vec<BigRational> = factor.iter().cloned().map(BigRational::from_integer).collect();
    let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(17));
    let mut out = Vec::new();
    for (j, interval) in isolate(&rational).into_iter().enumerate() {
        let (a, b) = refine_interval(&rational, interval, &width);
        let mid = ((a + b) / BigRational::from_integer(BigInt::from(2))).to_f64()?;
        out.push((mid, make(graph, factor, j)?));
    }
    Some(out)
}

/// Numeric semantics: the `k`-th real root of the polynomial whose
/// coefficients follow `k`.
fn evaluate(args: &[f64]) -> f64 {
    let Some((&k, coefficients)) = args.split_first() else { return f64::NAN };
    let descending: Vec<f64> = coefficients.iter().rev().copied().collect();
    let polynomial = crate::kernels::polynomial::Polynomial::new(descending);
    let Ok(mut roots) = crate::kernels::real_roots::find_roots(&polynomial, 1e-14) else { return f64::NAN };
    roots.sort_by(f64::total_cmp);
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * (1.0 + a.abs()));
    if k < 0.0 || k.fract() != 0.0 {
        return f64::NAN;
    }
    roots.get(k as usize).copied().unwrap_or(f64::NAN)
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
        let args = cx.graph.children(node).to_vec();
        let numeric = cx.env.numeric;
        let graph = &mut *cx.graph;
        let Some((&k, coefficients)) = args.split_first() else { return Outcome::Pass };
        let Some(k) = graph.number_of(k).and_then(Number::to_i64).and_then(|k| usize::try_from(k).ok()) else {
            return Outcome::Pass;
        };
        let q: Option<Vec<BigRational>> =
            coefficients.iter().map(|&c| graph.number_of(c).and_then(Number::to_rational)).collect();
        let Some(mut q) = q else { return Outcome::Pass };
        while q.last().is_some_and(Zero::is_zero) {
            q.pop();
        }
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
                let Some(canonical) = make(graph, &factor, j) else { return Outcome::Pass };
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
