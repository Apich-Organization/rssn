//! `root_of(k, c0, c1, ..., cn)`: the `k`-th root of a polynomial, exactly.
//!
//! The roots of `c0 + c1 x + ... + cn x^n` are counted from 0: the real
//! roots in ascending order, then the non-real roots by real
//! part and, for equal real parts, imaginary part.
//!
//! The coefficients are the arguments, so the object has
//! numeric semantics wherever the coefficients have (they may be
//! parameters); the kernel keeps the value canonical so that equal numbers
//! have equal representations:
//!
//! * `p` is factored over `Q`; the `k`-th root of `p` is located in one
//!   irreducible factor (the real roots by exact isolation, the complex
//!   ones by their numerical values);
//! * a root of a factor of degree one is a rational number and one of
//!   degree two is written with a square root;
//! * a root of a factor of degree three or more becomes
//!   `root_of(j, f0, ..., fm)` with `f` the primitive irreducible factor with
//!   a positive leading coefficient and `j` the index of the root among the
//!   roots of `f` (so `root_of(0, -1, -1, 0, 0, 0, 1)` is the one real
//!   root of `x^5 - x - 1`, and the index 1 root of `(x - 1)(x^5 - x - 1)` is
//!   the same object);
//! * an index that is not the index of a root leaves the term as it is.
//!   A quadratic factor with a negative discriminant gives a radical with
//!   the imaginary unit `I`.
//!
//! Complex evaluation (`Graph::eval_complex`) gives all the roots; real
//! evaluation gives the real ones and NaN for a complex index.
//!
//! In numeric mode the value of a real root is the midpoint of an isolating interval
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
use crate::graph::Arity;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::graph::rule::Installer;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new("root_of", Arity::Variadic)
        .cost(40)
        .eval(evaluate))?;
    i.graph()
        .ops_mut()
        .set_attr(op, crate::graph::ComplexEval(evaluate_complex));
    i.dim(op, DimRule::Custom(dim_root_of));
    i.kernel("poly/root_of", Tier::Normalize, RootOf { op });
    Ok(())
}

/// The dimension of a root of `c0 + c1 x + ... + cn x^n`.
///
/// Every term of the polynomial has the same dimension `D`, so
/// `dim(c_k) = D / dim(x)^k`. Two known coefficients determine `dim(x)`:
/// `dim(x) = (dim(c_i) / dim(c_j))^(1/(j-i))`. The remaining coefficients
/// only check that the polynomial is homogeneous; a coefficient that
/// contradicts the others is reported. If fewer than two coefficients have a
/// known dimension (all zero or all unknown values), the result is unknown.
fn dim_root_of(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let coefficients = args.get(1..).unwrap_or_default();
    let mut known: Vec<(usize, &DimTerm)> = Vec::new();
    for (power, coefficient) in coefficients.iter().enumerate() {
        if coefficient.is_non_metric() {
            return Err(DimError::Mismatch {
                left: coefficient.describe(ctx.arena),
                right: "a scalar coefficient".to_owned(),
                context: ctx.graph.display(ctx.node),
            });
        }
        if !coefficient.is_any() {
            known.push((power, coefficient));
        }
    }
    let (Some(&(i, ci)), Some(&(j, cj))) = (known.first(), known.get(1)) else {
        return Ok(ctx.any());
    };
    let exponent = BigRational::new(
        BigInt::one(),
        BigInt::from(i64::try_from(j - i).unwrap_or(i64::MAX)),
    );
    let x = ci.mul(&cj.inverse()).pow(&exponent);
    // D = dim(c_i) * dim(x)^i, and dim(c_k) must equal D / dim(x)^k.
    let leading = ci.mul(&integer_power(&x, i64::try_from(i).unwrap_or(i64::MAX)));
    for (power, coefficient) in known {
        let expected = leading.mul(&integer_power(
            &x,
            -i64::try_from(power).unwrap_or(i64::MAX),
        ));
        if coefficient != &expected {
            return Err(DimError::Mismatch {
                left: coefficient.describe(ctx.arena),
                right: expected.describe(ctx.arena),
                context: ctx.graph.display(ctx.node),
            });
        }
    }
    Ok(x)
}

/// `dimension^exponent` for an integer exponent.
fn integer_power(
    dimension: &DimTerm,
    exponent: i64,
) -> DimTerm {
    dimension.pow(&BigRational::from_integer(BigInt::from(exponent)))
}

struct RootOf {
    op: OpId,
}

/// A root of a primitive irreducible factor.
#[derive(Clone)]
struct Root {
    factor: Vec<BigInt>,
    /// Index among all roots of the factor (real ones first, then the
    /// complex ones).
    local: usize,
    /// An isolating interval of a real root.
    interval: Option<(BigRational, BigRational)>,
    /// Real and imaginary part (double precision).
    value: (f64, f64),
}

/// Complex ordering: by real part, then by imaginary part (real parts that
/// agree to `1e-9` count as equal, so conjugates are ordered by sign).
fn complex_order(
    a: (f64, f64),
    b: (f64, f64),
) -> std::cmp::Ordering {
    if (a.0 - b.0).abs() <= 1e-9 * (1.0 + a.0.abs()) {
        a.1.total_cmp(&b.1)
    } else {
        a.0.total_cmp(&b.0)
    }
}

/// All roots of the irreducible integer polynomial `factor` in the
/// canonical order: the real roots ascending, then the non-real roots by
/// real part and imaginary part.
fn factor_roots(factor: &[BigInt]) -> Vec<Root> {
    let rational: Vec<BigRational> = factor
        .iter()
        .cloned()
        .map(BigRational::from_integer)
        .collect();
    let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(40));
    let mut out = Vec::new();
    let intervals = isolate(&rational);
    let real_count = intervals.len();
    for (j, interval) in intervals.into_iter().enumerate() {
        let refined = refine_interval(&rational, interval, &width);
        let mid = ((&refined.0 + &refined.1) / BigRational::from_integer(BigInt::from(2)))
            .to_f64()
            .unwrap_or(f64::NAN);
        out.push(Root {
            factor: factor.to_vec(),
            local: j,
            interval: Some(refined),
            value: (mid, 0.0),
        });
    }
    let floats: Option<Vec<f64>> = factor.iter().map(ToPrimitive::to_f64).collect();
    if let Some(mut approximations) = floats.and_then(|c| super::galois::durand_kerner(&c)) {
        // The real roots are the ones with the smallest imaginary parts.
        approximations.sort_by(|a, b| a.1.abs().total_cmp(&b.1.abs()));
        let mut complex: Vec<(f64, f64)> = approximations.into_iter().skip(real_count).collect();
        complex.sort_by(|a, b| complex_order(*a, *b));
        for (i, value) in complex.into_iter().enumerate() {
            out.push(Root {
                factor: factor.to_vec(),
                local: real_count + i,
                interval: None,
                value,
            });
        }
    }
    out
}

/// The roots of all the irreducible factors of `q` in the canonical global
/// order (real roots ascending, then the complex ones).
fn all_roots(q: &[BigRational]) -> Vec<Root> {
    let (_, factors) = univariate::factor(q);
    let mut all: Vec<Root> = Vec::new();
    for (factor, _) in factors {
        all.extend(factor_roots(&factor));
    }
    all.sort_by(|a, b| {
        match (a.interval.is_some(), b.interval.is_some()) {
            | (true, true) => a.value.0.total_cmp(&b.value.0),
            | (true, false) => std::cmp::Ordering::Less,
            | (false, true) => std::cmp::Ordering::Greater,
            | (false, false) => complex_order(a.value, b.value),
        }
    });
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
    let rational: Vec<BigRational> = factor
        .iter()
        .cloned()
        .map(BigRational::from_integer)
        .collect();
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
    let Some((&k, coefficients)) = args.split_first() else {
        return f64::NAN;
    };
    let descending: Vec<f64> = coefficients.iter().rev().copied().collect();
    let polynomial = crate::kernels::polynomial::Polynomial::new(descending);
    let Ok(mut roots) = crate::kernels::real_roots::find_roots(&polynomial, 1e-14) else {
        return f64::NAN;
    };
    roots.sort_by(f64::total_cmp);
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * (1.0 + a.abs()));
    if k < 0.0 || k.fract() != 0.0 {
        return f64::NAN;
    }
    roots.get(k as usize).copied().unwrap_or(f64::NAN)
}

/// Complex semantics: index `k` counts the real roots first (ascending),
/// then the non-real roots by real and imaginary part.
fn evaluate_complex(args: &[num_complex::Complex64]) -> Option<num_complex::Complex64> {
    let (k, coefficients) = args.split_first()?;
    if k.im != 0.0 || k.re < 0.0 || k.re.fract() != 0.0 {
        return None;
    }
    let ascending: Vec<f64> = coefficients.iter().map(|c| c.re).collect();
    if coefficients.iter().any(|c| c.im != 0.0) {
        return None;
    }
    let approximations = super::galois::durand_kerner(&ascending)?;
    let tolerance = |re: f64| 1e-9 * (1.0 + re.abs());
    let mut real: Vec<f64> = approximations
        .iter()
        .filter(|r| r.1.abs() <= tolerance(r.0))
        .map(|r| r.0)
        .collect();
    real.sort_by(f64::total_cmp);
    let mut complex: Vec<(f64, f64)> = approximations
        .into_iter()
        .filter(|r| r.1.abs() > tolerance(r.0))
        .collect();
    complex.sort_by(|a, b| complex_order(*a, *b));
    let index = k.re as usize;
    match real.get(index) {
        | Some(&r) => Some(num_complex::Complex64::new(r, 0.0)),
        | None => {
            complex
                .get(index - real.len())
                .map(|&(re, im)| num_complex::Complex64::new(re, im))
        },
    }
}

/// Every root of a polynomial with integer coefficients (ascending) as
/// canonical `root_of` terms with their values, in the canonical order.
/// Degree-one and degree-two factors are not handled here.
pub(crate) fn all_roots_as_objects(
    graph: &mut Graph,
    factor: &[BigInt],
) -> Option<Vec<((f64, f64), NodeId)>> {
    let mut out = Vec::new();
    for root in factor_roots(factor) {
        out.push((root.value, make(graph, factor, root.local)?));
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
        let args = cx.graph.children(node).to_vec();
        let numeric = cx.env.numeric;
        let graph = &mut *cx.graph;
        let Some((&k, coefficients)) = args.split_first() else {
            return Outcome::Pass;
        };
        let Some(k) = graph
            .number_of(k)
            .and_then(Number::to_i64)
            .and_then(|k| usize::try_from(k).ok())
        else {
            return Outcome::Pass;
        };
        let q: Option<Vec<BigRational>> = coefficients
            .iter()
            .map(|&c| graph.number_of(c).and_then(Number::to_rational))
            .collect();
        let Some(mut q) = q else { return Outcome::Pass };
        while q.last().is_some_and(Zero::is_zero) {
            q.pop();
        }
        if q.len() < 2 {
            return Outcome::Pass;
        }
        let roots = all_roots(&q);
        let Some(Root {
            factor,
            local: j,
            interval,
            ..
        }) = roots.get(k).cloned()
        else {
            return Outcome::Pass;
        };
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
                let centre = graph.num(Number::rat(BigRational::new(
                    -b.clone(),
                    denominator.clone(),
                )));
                let negative = discriminant.is_negative();
                let magnitude = if negative { -inside } else { inside };
                let scale = BigRational::new(outside * BigInt::from(sign), denominator);
                let (scale_node, radicand) = (
                    graph.num(Number::rat(scale)),
                    graph.num(Number::Int(magnitude)),
                );
                let Some(half) = Number::fraction(1, 2) else {
                    return Outcome::Pass;
                };
                let half = graph.num(half);
                let radical = graph.node(core::POW, &[radicand, half]);
                let mut factors = vec![scale_node, radical];
                if negative {
                    let Some(unit) = graph.ops().lookup("I") else {
                        return Outcome::Pass;
                    };
                    factors.push(graph.node(unit, &[]));
                }
                let offset = graph.node(core::MUL, &factors);
                Outcome::Equal(graph.node(core::ADD, &[centre, offset]))
            },
            | _ => {
                let Some(canonical) = make(graph, &factor, j) else {
                    return Outcome::Pass;
                };
                if graph.find(canonical) != graph.find(node) && !graph.same(canonical, node) {
                    return Outcome::Equal(canonical);
                }
                let Some(interval) = interval else {
                    return Outcome::Pass;
                };
                if !numeric {
                    return Outcome::Pass;
                }
                let rational: Vec<BigRational> = factor
                    .iter()
                    .cloned()
                    .map(BigRational::from_integer)
                    .collect();
                let width = BigRational::new(BigInt::one(), BigInt::from(10).pow(17));
                let (a, b) = refine_interval(&rational, interval, &width);
                let (lo, hi) = (
                    a.to_f64().unwrap_or(f64::NAN),
                    b.to_f64().unwrap_or(f64::NAN),
                );
                let mid = 0.5 * lo + 0.5 * hi;
                Outcome::Approx(Ball {
                    mid,
                    rad: (hi - mid).max(mid - lo),
                })
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

/// Canonical order of root values: real ones (imaginary part negligible)
/// first, ascending, then the others by real and imaginary part.
pub(crate) fn canonical_order(
    a: (f64, f64),
    b: (f64, f64),
) -> std::cmp::Ordering {
    let real = |v: (f64, f64)| v.1.abs() <= 1e-12 * (1.0 + v.0.abs());
    match (real(a), real(b)) {
        | (true, true) => a.0.total_cmp(&b.0),
        | (true, false) => std::cmp::Ordering::Less,
        | (false, true) => std::cmp::Ordering::Greater,
        | (false, false) => complex_order(a, b),
    }
}
