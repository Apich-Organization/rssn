//! Statistics: descriptive measures of lists, probability distributions as
//! terms, hypothesis tests, regression and information theory.
//!
//! # Lists
//!
//! Sample statistics take `list(...)` arguments. When every item is an
//! exact number the answer is exact (`mean(list(1, 2, 4))` is `7/3`, the
//! standard deviation stays a radical); when items are floats the numeric
//! routines of [`crate::kernels::stats`] are used where they agree with the
//! definitions below; when some item is symbolic the result is the defining
//! formula built from the items, ready for further simplification.
//! Definitions: `variance` and `std` are the *population* quantities
//! (divide by `n`), `sample_variance` and `sample_std` divide by `n - 1`;
//! `covariance` is the population covariance; `moment(list, k)` is the
//! `k`-th central moment; `skewness` is `m3 / m2^(3/2)` and `kurtosis`
//! is `m4 / m2^2` (Pearson, not excess); `quantile(list, p)` interpolates
//! linearly between order statistics (`p` in `[0, 1]`), `percentile(list, p)`
//! is `quantile(list, p / 100)`. Order statistics need numbers: a list with
//! symbolic items stays unreduced.
//!
//! # Distributions
//!
//! `normal(mu, sigma)`, `uniform(a, b)`, `exponential(rate)`,
//! `bernoulli(p)`, `binomial_dist(n, p)`, `poisson(rate)`,
//! `gamma_dist(shape, rate)`, `beta_dist(a, b)`, `student_t(nu)` and
//! `chi_squared(k)` are inert terms. The requests `pdf(d, x)` (mass
//! function for discrete ones), `cdf(d, x)`, `expectation(d)`,
//! `variance_of(d)`, `mgf(d, t)` and `entropy_of(d)` are reduced to closed
//! forms built from `exp`, `erf`, `gamma`, `beta`, `heaviside` and the two
//! regularised incomplete functions `gamma_lr` and `beta_reg` that this
//! rule set adds. Bounded supports are handled with `heaviside`, so the
//! density is zero (and the distribution function 0 or 1) outside the
//! support; the discrete distribution functions are meant for integer `x`.
//! Requests without a closed form (the moment generating function of the
//! beta and Student distributions, the entropy of the binomial and Poisson
//! distributions) stay requests.
//!
//! # Inference, regression, information theory
//!
//! `z_test`, `t_test`, `welch_test`, `pooled_t_test`, `chi_squared_test`
//! and `anova` return `list(statistic, p_value)`: floats for numeric data,
//! formulas (with the p-value written through the `cdf` request) for
//! symbolic data. `confidence_interval` and `z_interval` are numeric.
//! `linear_regression(xs, ys)` returns `list(intercept, slope)`, exact for
//! exact data and symbolic for symbolic data; with a list of predictor
//! columns it is multiple regression, `polynomial_regression(xs, ys, d)`
//! fits a degree `d` polynomial. `entropy`, `kl_divergence`,
//! `cross_entropy`, `joint_entropy`, `conditional_entropy`,
//! `mutual_information` and `gini` are natural-logarithm (nats) formulas on
//! probability lists (matrices for the joint ones), exact for exact
//! probabilities.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use statrs::distribution::StudentsT;
use statrs::function::beta::beta_reg as statrs_beta_reg;
use statrs::function::erf::erfc;
use statrs::function::gamma::gamma_lr as statrs_gamma_lr;

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
use crate::graph::Pat;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::VarNames;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::stats as num;

use super::calculus::calculus;
use super::calculus::partials;
use super::combinatorics::choose;
use super::elementary::elementary;
use super::number_theory::factor as prime_factors;
use super::number_theory::procedure;
use super::poly::best;
use super::special::special;

/// Longest list whose statistics are expanded into terms or solved exactly.
const MAX_ITEMS: usize = 5000;
/// Largest exponent of a moment.
const MAX_MOMENT: i64 = 64;
/// Highest degree of a polynomial regression.
const MAX_DEGREE: usize = 16;
/// Largest `n` (or `k`) for which a discrete distribution function is
/// written as an explicit finite sum.
const MAX_SUM: u64 = 2000;

/// The statistics rule set.
#[must_use]
pub fn stats() -> RuleSet {
    RuleSet::new("stats", install)
        .needs(elementary())
        .needs(calculus())
        .needs(special())
}

// ---------------------------------------------------------------------
// Term construction
// ---------------------------------------------------------------------

/// An item count as a float.
fn count(n: usize) -> f64 {
    n.to_f64().unwrap_or(f64::NAN)
}

fn is_zero_number(
    g: &Graph,
    node: NodeId,
) -> bool {
    g.number_of(node).is_some_and(Number::is_zero)
}

/// The sum of `terms`; numeric terms are added up on the spot.
fn sum_node(
    g: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    let mut constant: Option<Number> = None;
    let mut rest = Vec::with_capacity(terms.len());
    let mut stack: Vec<NodeId> = terms.iter().rev().copied().collect();
    while let Some(t) = stack.pop() {
        match g.number_of(t).cloned() {
            | Some(n) => {
                constant = Some(match constant {
                    | Some(c) => c.add(&n),
                    | None => n,
                });
            },
            | None if g.op(t) == core::ADD => stack.extend(g.children(t).iter().rev().copied()),
            | None => rest.push(t),
        }
    }
    if let Some(c) = constant {
        if !c.is_zero() || rest.is_empty() {
            rest.insert(0, g.num(c));
        }
    }
    match rest.as_slice() {
        | [] => g.int(0),
        | [only] => *only,
        | _ => g.node(core::ADD, &rest),
    }
}

/// The product of `factors`; numeric factors are multiplied on the spot.
fn product_node(
    g: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    let mut constant: Option<Number> = None;
    let mut rest = Vec::with_capacity(factors.len());
    let mut stack: Vec<NodeId> = factors.iter().rev().copied().collect();
    while let Some(f) = stack.pop() {
        match g.number_of(f).cloned() {
            | Some(n) => {
                constant = Some(match constant {
                    | Some(c) => c.mul(&n),
                    | None => n,
                });
            },
            | None if g.op(f) == core::MUL => stack.extend(g.children(f).iter().rev().copied()),
            | None => rest.push(f),
        }
    }
    if let Some(c) = constant {
        if c.is_zero() {
            return g.num(c);
        }
        if !c.is_one() || rest.is_empty() {
            rest.insert(0, g.num(c));
        }
    }
    match rest.as_slice() {
        | [] => g.int(1),
        | [only] => *only,
        | _ => g.node(core::MUL, &rest),
    }
}

/// `base^exponent`, evaluated when both are numbers and the result is
/// exact or a float.
fn pow_node(
    g: &mut Graph,
    base: NodeId,
    exponent: &Number,
) -> NodeId {
    if exponent.is_one() {
        return base;
    }
    if let Some(v) = g.number_of(base).and_then(|b| b.pow(exponent)) {
        return g.num(v);
    }
    let e = g.num(exponent.clone());
    g.node(core::POW, &[base, e])
}

fn scale(
    g: &mut Graph,
    node: NodeId,
    by: &Number,
) -> NodeId {
    let factor = g.num(by.clone());
    product_node(g, &[factor, node])
}

fn negate(
    g: &mut Graph,
    node: NodeId,
) -> NodeId {
    scale(g, node, &Number::from(-1))
}

fn difference(
    g: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus_b = negate(g, b);
    sum_node(g, &[a, minus_b])
}

/// `1 / n` as a number.
fn reciprocal_of(n: usize) -> Option<Number> {
    Number::fraction(1, i64::try_from(n).ok()?)
}

/// The term of the operator called `name`, if it is installed.
fn call(
    g: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = g.ops().lookup(name)?;
    g.try_node(op, args)
}

fn list_node(
    g: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    g.node(core::LIST, items)
}

fn items_of(
    g: &Graph,
    list: NodeId,
) -> Option<Vec<NodeId>> {
    (g.op(list) == core::LIST).then(|| g.children(list).to_vec())
}

fn f64_of(
    g: &Graph,
    node: NodeId,
) -> Option<f64> {
    g.number_of(node)
        .map(Number::to_f64)
        .filter(|v| v.is_finite())
}

/// The items as floats if all are numbers and at least one is a float.
fn float_items(
    g: &Graph,
    xs: &[NodeId],
) -> Option<Vec<f64>> {
    let all_numbers = xs.iter().all(|&x| g.number_of(x).is_some());
    let any_float = xs
        .iter()
        .any(|&x| matches!(g.number_of(x), Some(Number::Float(_))));
    if !all_numbers || !any_float {
        return None;
    }
    xs.iter().map(|&x| f64_of(g, x)).collect()
}

/// The items as numbers, if they all are.
fn number_items(
    g: &Graph,
    xs: &[NodeId],
) -> Option<Vec<Number>> {
    xs.iter().map(|&x| g.number_of(x).cloned()).collect()
}

/// The exact rational value of a number, converting floats exactly.
fn rational_of(n: &Number) -> Option<BigRational> {
    match n {
        | Number::Float(v) => BigRational::from_float(*v),
        | other => other.to_rational(),
    }
}

/// A number node for `r`, as a float if `float`.
fn value_node(
    g: &mut Graph,
    r: &BigRational,
    float: bool,
) -> NodeId {
    if float {
        g.float(r.to_f64().unwrap_or(f64::NAN))
    } else {
        g.num(Number::rat(r.clone()))
    }
}

// ---------------------------------------------------------------------
// Descriptive statistics
// ---------------------------------------------------------------------

type ListFn = fn(&mut Graph, &[NodeId]) -> Option<NodeId>;

/// Wraps a function of one list into a kernel body.
fn list_op(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    f: ListFn,
) -> Outcome {
    let &[list] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let Some(items) = items_of(g, list).filter(|i| i.len() <= MAX_ITEMS) else {
        return Outcome::Pass;
    };
    let result = f(g, &items);
    finish(result)
}

fn finish(result: Option<NodeId>) -> Outcome {
    result.map_or(Outcome::Pass, Outcome::Equal)
}

type PairFn = fn(&mut Graph, &[NodeId], &[NodeId]) -> Option<NodeId>;

/// Wraps a function of two lists into a kernel body.
fn pair_op(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    f: PairFn,
) -> Outcome {
    let &[a, b] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys)) = (items_of(g, a), items_of(g, b)) else {
        return Outcome::Pass;
    };
    if xs.len() > MAX_ITEMS || xs.len() != ys.len() {
        return Outcome::Pass;
    }
    let result = f(g, &xs, &ys);
    finish(result)
}

type ScalarFn = fn(&mut Graph, &[NodeId], NodeId) -> Option<NodeId>;

/// Wraps a function of a list and one more argument into a kernel body.
fn list_with_op(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    f: ScalarFn,
) -> Outcome {
    let &[list, extra] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let Some(items) = items_of(g, list).filter(|i| i.len() <= MAX_ITEMS) else {
        return Outcome::Pass;
    };
    let result = f(g, &items, extra);
    finish(result)
}

fn mean_term(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let inverse = reciprocal_of(xs.len())?;
    let total = sum_node(g, xs);
    Some(scale(g, total, &inverse))
}

/// `sum (x - mean)^k / divisor`.
fn central_moment(
    g: &mut Graph,
    xs: &[NodeId],
    k: i64,
    divisor: usize,
) -> Option<NodeId> {
    let mean = mean_term(g, xs)?;
    let minus_mean = negate(g, mean);
    let mut terms = Vec::with_capacity(xs.len());
    for &x in xs {
        let d = sum_node(g, &[x, minus_mean]);
        terms.push(pow_node(g, d, &Number::from(k)));
    }
    let total = sum_node(g, &terms);
    Some(scale(g, total, &reciprocal_of(divisor)?))
}

fn mean(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    if xs.is_empty() {
        return None;
    }
    if let Some(v) = float_items(g, xs) {
        return Some(g.float(num::mean(&v)));
    }
    mean_term(g, xs)
}

fn variance(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    if xs.is_empty() {
        return None;
    }
    if let Some(v) = float_items(g, xs) {
        let variance = num::variance_with_type(&v, num::VarianceType::Population)?;
        return Some(g.float(variance));
    }
    central_moment(g, xs, 2, xs.len())
}

fn sample_variance(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    if xs.len() < 2 {
        return None;
    }
    if let Some(v) = float_items(g, xs) {
        let variance = num::variance_with_type(&v, num::VarianceType::Sample)?;
        return Some(g.float(variance));
    }
    central_moment(g, xs, 2, xs.len() - 1)
}

fn square_root(
    g: &mut Graph,
    node: NodeId,
) -> NodeId {
    pow_node(
        g,
        node,
        &Number::fraction(1, 2).unwrap_or_else(|| Number::from(1)),
    )
}

fn std_dev(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let v = variance(g, xs)?;
    Some(square_root(g, v))
}

fn sample_std(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let v = sample_variance(g, xs)?;
    Some(square_root(g, v))
}

/// The items with their numeric values, ascending.
fn sorted_numbers(
    g: &Graph,
    xs: &[NodeId],
) -> Option<Vec<(Number, NodeId)>> {
    let mut items: Vec<(Number, NodeId)> = xs
        .iter()
        .map(|&x| g.number_of(x).cloned().map(|n| (n, x)))
        .collect::<Option<_>>()?;
    if items
        .iter()
        .any(|(n, _)| matches!(n, Number::Float(v) if !v.is_finite()))
    {
        return None;
    }
    items.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(items)
}

fn median(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let n = xs.len();
    if n == 0 {
        return None;
    }
    if let Some(mut v) = float_items(g, xs) {
        return Some(g.float(num::median(&mut v)));
    }
    let ordered: Vec<NodeId> = match sorted_numbers(g, xs) {
        | Some(sorted) => sorted.into_iter().map(|(_, x)| x).collect(),
        // Without numbers only the trivial cases have a middle.
        | None if n <= 2 => xs.to_vec(),
        | None => return None,
    };
    if n % 2 == 1 {
        return ordered.get(n / 2).copied();
    }
    let middle = ordered.get(n / 2 - 1..=n / 2)?;
    mean_term(g, middle)
}

fn mode(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let mut groups: Vec<(NodeId, usize)> = Vec::new();
    for &x in xs {
        match groups.iter_mut().find(|(rep, _)| g.same(*rep, x)) {
            | Some(group) => group.1 += 1,
            | None => groups.push((x, 1)),
        }
    }
    let top = groups.iter().map(|g| g.1).max()?;
    if top == 1 && xs.len() > 1 {
        return None;
    }
    let mut winners = groups.iter().filter(|g| g.1 == top).map(|g| g.0);
    let first = winners.next()?;
    // Among equally frequent numbers the smallest one.
    Some(
        winners.fold(first, |best, x| match (g.number_of(x), g.number_of(best)) {
            | (Some(a), Some(b)) if a.total_cmp(b).is_lt() => x,
            | _ => best,
        }),
    )
}

fn covariance(
    g: &mut Graph,
    xs: &[NodeId],
    ys: &[NodeId],
) -> Option<NodeId> {
    let n = xs.len();
    if n == 0 || n != ys.len() {
        return None;
    }
    if let (Some(a), Some(b)) = (float_items(g, xs), float_items(g, ys)) {
        if n >= 2 {
            // The kernel divides by n - 1.
            let factor = (count(n) - 1.0) / count(n);
            return Some(g.float(num::covariance(&a, &b) * factor));
        }
    }
    let (mx, my) = (mean_term(g, xs)?, mean_term(g, ys)?);
    let (minus_mx, minus_my) = (negate(g, mx), negate(g, my));
    let mut terms = Vec::with_capacity(n);
    for (&x, &y) in xs.iter().zip(ys) {
        let dx = sum_node(g, &[x, minus_mx]);
        let dy = sum_node(g, &[y, minus_my]);
        terms.push(product_node(g, &[dx, dy]));
    }
    let total = sum_node(g, &terms);
    Some(scale(g, total, &reciprocal_of(n)?))
}

fn correlation(
    g: &mut Graph,
    xs: &[NodeId],
    ys: &[NodeId],
) -> Option<NodeId> {
    let cov = covariance(g, xs, ys)?;
    let (vx, vy) = (
        central_moment(g, xs, 2, xs.len())?,
        central_moment(g, ys, 2, ys.len())?,
    );
    if is_zero_number(g, vx) || is_zero_number(g, vy) {
        return None;
    }
    let both = product_node(g, &[vx, vy]);
    let scale_factor = pow_node(g, both, &Number::fraction(-1, 2)?);
    Some(product_node(g, &[cov, scale_factor]))
}

fn moment(
    g: &mut Graph,
    xs: &[NodeId],
    k: NodeId,
) -> Option<NodeId> {
    let k = g
        .number_of(k)?
        .to_i64()
        .filter(|k| (0..=MAX_MOMENT).contains(k))?;
    if xs.is_empty() {
        return None;
    }
    central_moment(g, xs, k, xs.len())
}

/// `m_k m_2^(-k/2)`, the standardised moment.
fn standardised(
    g: &mut Graph,
    xs: &[NodeId],
    k: i64,
) -> Option<NodeId> {
    if xs.is_empty() {
        return None;
    }
    let m2 = central_moment(g, xs, 2, xs.len())?;
    if is_zero_number(g, m2) {
        return None;
    }
    let mk = central_moment(g, xs, k, xs.len())?;
    let power = pow_node(g, m2, &Number::fraction(-k, 2)?);
    Some(product_node(g, &[mk, power]))
}

fn skewness(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    standardised(g, xs, 3)
}

fn kurtosis(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    standardised(g, xs, 4)
}

/// The `p`-quantile of sorted numbers by linear interpolation between
/// order statistics (`h = (n - 1) p`).
fn quantile_of(
    g: &mut Graph,
    sorted: &[(Number, NodeId)],
    p: &Number,
) -> Option<NodeId> {
    let last = sorted.len().checked_sub(1)?;
    let p_value = p.to_f64();
    if !(0.0..=1.0).contains(&p_value) {
        return None;
    }
    let h = p.mul(&Number::from(i64::try_from(last).ok()?));
    let lower = match h.to_rational() {
        | Some(r) => r.floor().to_integer().to_usize()?,
        | None => h.to_f64().floor().to_usize()?,
    };
    let (_, low_node) = sorted.get(lower)?;
    let Some((_, high_node)) = sorted.get(lower + 1) else {
        return Some(*low_node);
    };
    let fraction = h.add(&Number::from(i64::try_from(lower).ok()?).neg());
    let spread = difference(g, *high_node, *low_node);
    let step = scale(g, spread, &fraction);
    Some(sum_node(g, &[*low_node, step]))
}

fn quantile(
    g: &mut Graph,
    xs: &[NodeId],
    p: NodeId,
) -> Option<NodeId> {
    let p = g.number_of(p)?.clone();
    let sorted = sorted_numbers(g, xs)?;
    quantile_of(g, &sorted, &p)
}

fn percentile(
    g: &mut Graph,
    xs: &[NodeId],
    p: NodeId,
) -> Option<NodeId> {
    let p = g.number_of(p)?.mul(&Number::fraction(1, 100)?);
    let sorted = sorted_numbers(g, xs)?;
    quantile_of(g, &sorted, &p)
}

fn iqr(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let sorted = sorted_numbers(g, xs)?;
    let upper = quantile_of(g, &sorted, &Number::fraction(3, 4)?)?;
    let lower = quantile_of(g, &sorted, &Number::fraction(1, 4)?)?;
    Some(difference(g, upper, lower))
}

fn extreme(
    g: &mut Graph,
    xs: &[NodeId],
    largest: bool,
) -> Option<NodeId> {
    if let Some(mut v) = float_items(g, xs) {
        let value = if largest {
            num::max(&mut v)
        } else {
            num::min(&mut v)
        };
        return Some(g.float(value));
    }
    let sorted = sorted_numbers(g, xs)?;
    let pick = if largest {
        sorted.last()
    } else {
        sorted.first()
    };
    pick.map(|(_, x)| *x)
}

fn min_of(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    extreme(g, xs, false)
}

fn max_of(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    extreme(g, xs, true)
}

fn range_of(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let (low, high) = (min_of(g, xs)?, max_of(g, xs)?);
    Some(difference(g, high, low))
}

fn geometric_mean(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    if xs.is_empty()
        || xs.iter().any(|&x| {
            g.number_of(x)
                .is_some_and(|n| n.to_f64().partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater))
        })
    {
        return None;
    }
    let product = product_node(g, xs);
    Some(pow_node(g, product, &reciprocal_of(xs.len())?))
}

fn harmonic_mean(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    if xs.is_empty() || xs.iter().any(|&x| is_zero_number(g, x)) {
        return None;
    }
    let inverses: Vec<NodeId> = xs
        .iter()
        .map(|&x| pow_node(g, x, &Number::from(-1)))
        .collect();
    let total = sum_node(g, &inverses);
    let count = g.int(i64::try_from(xs.len()).ok()?);
    let inverse_total = pow_node(g, total, &Number::from(-1));
    Some(product_node(g, &[count, inverse_total]))
}

fn zscores(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let m = mean_term(g, xs)?;
    let v = central_moment(g, xs, 2, xs.len())?;
    let minus_m = negate(g, m);
    let spread = pow_node(g, v, &Number::fraction(-1, 2)?);
    let mut out = Vec::with_capacity(xs.len());
    for &x in xs {
        let centred = sum_node(g, &[x, minus_m]);
        out.push(if is_zero_number(g, v) {
            g.int(0)
        } else {
            product_node(g, &[centred, spread])
        });
    }
    Some(list_node(g, &out))
}

fn coefficient_of_variation(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let m = mean_term(g, xs)?;
    if is_zero_number(g, m) {
        return None;
    }
    let s = sample_std(g, xs)?;
    let inverse = pow_node(g, m, &Number::from(-1));
    Some(product_node(g, &[s, inverse]))
}

fn standard_error(
    g: &mut Graph,
    xs: &[NodeId],
) -> Option<NodeId> {
    let v = sample_variance(g, xs)?;
    let per_item = scale(g, v, &reciprocal_of(xs.len())?);
    Some(square_root(g, per_item))
}

// ---------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------

/// The closed forms of one distribution family, as pattern text over the
/// parameters `?a`, `?b` and the argument `?x`.
struct Family {
    name: &'static str,
    params: u8,
    pdf: &'static str,
    cdf: &'static str,
    mean: &'static str,
    variance: &'static str,
    mgf: Option<&'static str>,
    entropy: Option<&'static str>,
}

const FAMILIES: [Family; 10] = [
    Family {
        name: "normal",
        params: 2,
        pdf: "exp(-(?x - ?a)^2 / (2 * ?b^2)) / (?b * (2 * pi)^(1/2))",
        cdf: "1/2 * (1 + erf((?x - ?a) / (?b * 2^(1/2))))",
        mean: "?a",
        variance: "?b^2",
        mgf: Some("exp(?a * ?x + ?b^2 * ?x^2 / 2)"),
        entropy: Some("1/2 + ln(2 * pi * ?b^2) / 2"),
    },
    Family {
        name: "uniform",
        params: 2,
        pdf: "heaviside(?x - ?a) * heaviside(?b - ?x) / (?b - ?a)",
        cdf: "(?x - ?a) / (?b - ?a) * heaviside(?x - ?a) * heaviside(?b - ?x) + heaviside(?x - ?b)",
        mean: "(?a + ?b) / 2",
        variance: "(?b - ?a)^2 / 12",
        mgf: Some("(exp(?b * ?x) - exp(?a * ?x)) / (?x * (?b - ?a))"),
        entropy: Some("ln(?b - ?a)"),
    },
    Family {
        name: "exponential",
        params: 1,
        pdf: "?a * exp(-?a * ?x) * heaviside(?x)",
        cdf: "(1 - exp(-?a * ?x)) * heaviside(?x)",
        mean: "1 / ?a",
        variance: "1 / ?a^2",
        mgf: Some("?a / (?a - ?x)"),
        entropy: Some("1 - ln(?a)"),
    },
    Family {
        name: "bernoulli",
        params: 1,
        pdf: "?a^?x * (1 - ?a)^(1 - ?x)",
        cdf: "(1 - ?a) * heaviside(?x + 1/2) + ?a * heaviside(?x - 1/2)",
        mean: "?a",
        variance: "?a * (1 - ?a)",
        mgf: Some("1 - ?a + ?a * exp(?x)"),
        entropy: Some("-?a * ln(?a) - (1 - ?a) * ln(1 - ?a)"),
    },
    Family {
        name: "binomial_dist",
        params: 2,
        pdf: "gamma(?a + 1) / (gamma(?x + 1) * gamma(?a - ?x + 1)) * ?b^?x * (1 - ?b)^(?a - ?x)",
        cdf: "beta_reg(?a - ?x, ?x + 1, 1 - ?b)",
        mean: "?a * ?b",
        variance: "?a * ?b * (1 - ?b)",
        mgf: Some("(1 - ?b + ?b * exp(?x))^?a"),
        entropy: None,
    },
    Family {
        name: "poisson",
        params: 1,
        pdf: "?a^?x * exp(-?a) / gamma(?x + 1)",
        cdf: "1 - gamma_lr(?x + 1, ?a)",
        mean: "?a",
        variance: "?a",
        mgf: Some("exp(?a * (exp(?x) - 1))"),
        entropy: None,
    },
    Family {
        name: "gamma_dist",
        params: 2,
        pdf: "?b^?a * abs(?x)^(?a - 1) * exp(-?b * ?x) / gamma(?a) * heaviside(?x)",
        cdf: "gamma_lr(?a, ?b * ?x)",
        mean: "?a / ?b",
        variance: "?a / ?b^2",
        mgf: Some("(1 - ?x / ?b)^(-?a)"),
        entropy: Some("?a - ln(?b) + lgamma(?a) + (1 - ?a) * digamma(?a)"),
    },
    Family {
        name: "beta_dist",
        params: 2,
        pdf: "abs(?x)^(?a - 1) * abs(1 - ?x)^(?b - 1) / beta(?a, ?b) * heaviside(?x) * heaviside(1 - ?x)",
        cdf: "beta_reg(?a, ?b, ?x)",
        mean: "?a / (?a + ?b)",
        variance: "?a * ?b / ((?a + ?b)^2 * (?a + ?b + 1))",
        mgf: None,
        entropy: Some(
            "lgamma(?a) + lgamma(?b) - lgamma(?a + ?b) - (?a - 1) * digamma(?a) - (?b - 1) * digamma(?b) \
             + (?a + ?b - 2) * digamma(?a + ?b)",
        ),
    },
    Family {
        name: "student_t",
        params: 1,
        pdf: "gamma((?a + 1) / 2) / ((?a * pi)^(1/2) * gamma(?a / 2)) * (1 + ?x^2 / ?a)^(-(?a + 1) / 2)",
        cdf: "1/2 + sign(?x) / 2 * (1 - beta_reg(?a / 2, 1/2, ?a / (?a + ?x^2)))",
        mean: "0",
        variance: "?a / (?a - 2)",
        mgf: None,
        entropy: Some(
            "(?a + 1) / 2 * (digamma((?a + 1) / 2) - digamma(?a / 2)) + ln(?a) / 2 + lgamma(?a / 2) \
             + lgamma(1/2) - lgamma((?a + 1) / 2)",
        ),
    },
    Family {
        name: "chi_squared",
        params: 1,
        pdf: "abs(?x)^(?a / 2 - 1) * exp(-?x / 2) / (2^(?a / 2) * gamma(?a / 2)) * heaviside(?x)",
        cdf: "gamma_lr(?a / 2, ?x / 2)",
        mean: "?a",
        variance: "2 * ?a",
        mgf: Some("(1 - 2 * ?x)^(-?a / 2)"),
        entropy: Some("?a / 2 + ln(2) + lgamma(?a / 2) + (1 - ?a / 2) * digamma(?a / 2)"),
    },
];

/// Regularised lower incomplete gamma function `P(a, x)`; zero for
/// `x <= 0`.
fn gamma_lr_eval(args: &[f64]) -> f64 {
    let (Some(&a), Some(&x)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if a.is_nan() || x.is_nan() || a <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    statrs_gamma_lr(a, x)
}

/// Regularised incomplete beta function `I_x(a, b)`, clamped to `[0, 1]`
/// outside `(0, 1)`.
fn beta_reg_eval(args: &[f64]) -> f64 {
    let (Some(&a), Some(&b), Some(&x)) = (args.first(), args.get(1), args.get(2)) else {
        return f64::NAN;
    };
    if a.is_nan() || b.is_nan() || x.is_nan() || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    statrs_beta_reg(a, b, x)
}

#[derive(Copy, Clone)]
enum Request {
    Pdf,
    Cdf,
    Mean,
    Variance,
    Mgf,
    Entropy,
}

/// Whether numeric parameters are inside the parameter space of the
/// family; symbolic parameters are taken on trust.
fn admissible(
    family: &Family,
    params: &[Option<f64>],
) -> bool {
    let p = |i: usize| params.get(i).copied().flatten();
    let positive = |i: usize| p(i).is_none_or(|v| v > 0.0);
    let probability = |i: usize| p(i).is_none_or(|v| (0.0..=1.0).contains(&v));
    match family.name {
        | "normal" => positive(1),
        | "uniform" => match (p(0), p(1)) {
            | (Some(a), Some(b)) => a < b,
            | _ => true,
        },
        | "exponential" | "poisson" | "student_t" | "chi_squared" => positive(0),
        | "bernoulli" => probability(0),
        | "binomial_dist" => p(0).is_none_or(|n| n >= 0.0 && n.fract() == 0.0) && probability(1),
        | _ => positive(0) && positive(1),
    }
}

/// The finite sum for the distribution function of a binomial or Poisson
/// variable at a literal integer, or `None` if the formula applies.
fn discrete_sum(
    g: &mut Graph,
    name: &str,
    params: &[NodeId],
    k: NodeId,
) -> Option<NodeId> {
    let k = g.number_of(k).and_then(Number::to_i64)?;
    match (name, params) {
        | ("binomial_dist", &[n, p]) => {
            let n = g
                .number_of(n)
                .and_then(Number::to_i64)
                .filter(|&n| (0..=2000).contains(&n))?;
            if k < 0 {
                return Some(g.int(0));
            }
            if k >= n {
                return Some(g.int(1));
            }
            let q = {
                let minus_p = negate(g, p);
                let one = g.int(1);
                sum_node(g, &[one, minus_p])
            };
            let mut terms = Vec::new();
            for j in 0..=k {
                let c = choose(&BigInt::from(n), u64::try_from(j).ok()?)?;
                let c = g.num(Number::Int(c));
                let (pj, qj) = (
                    pow_node(g, p, &Number::from(j)),
                    pow_node(g, q, &Number::from(n - j)),
                );
                terms.push(product_node(g, &[c, pj, qj]));
            }
            Some(sum_node(g, &terms))
        },
        | ("poisson", &[rate]) => {
            if k < 0 {
                return Some(g.int(0));
            }
            let k = u64::try_from(k).ok().filter(|&k| k <= MAX_SUM)?;
            let mut terms = Vec::new();
            let mut factorial = BigInt::one();
            for j in 0..=k {
                if j > 0 {
                    factorial *= j;
                }
                let inverse = g.num(Number::rat(BigRational::new(
                    BigInt::one(),
                    factorial.clone(),
                )));
                let power = pow_node(g, rate, &Number::from(i64::try_from(j).ok()?));
                terms.push(product_node(g, &[inverse, power]));
            }
            let total = sum_node(g, &terms);
            let minus_rate = negate(g, rate);
            let decay = call(g, "exp", &[minus_rate])?;
            Some(product_node(g, &[decay, total]))
        },
        | _ => None,
    }
}

/// Instantiates pattern text over `?a`, `?b`, `?x`.
fn template(
    g: &mut Graph,
    text: &str,
    args: [NodeId; 3],
) -> Option<NodeId> {
    let mut vars = VarNames::default();
    for name in ["a", "b", "x"] {
        vars.index(name);
    }
    let pat = Pat::parse(text, g, &mut vars).ok()?;
    pat.instantiate(g, &args)
}

fn distribution_request(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    request: Request,
) -> Outcome {
    let g = &mut *cx.graph;
    let (&dist, extra) = match args.split_first() {
        | Some(split) => split,
        | None => return Outcome::Pass,
    };
    let expected = match request {
        | Request::Pdf | Request::Cdf | Request::Mgf => 1,
        | Request::Mean | Request::Variance | Request::Entropy => 0,
    };
    if extra.len() != expected {
        return Outcome::Pass;
    }
    let Some(dist) = best(g, dist) else {
        return Outcome::Pass;
    };
    let name = g.ops().get(g.op(dist)).name.clone();
    let Some(family) = FAMILIES.iter().find(|f| *f.name == *name) else {
        return Outcome::Pass;
    };
    let params = g.children(dist).to_vec();
    if params.len() != usize::from(family.params) {
        return Outcome::Pass;
    }
    let numeric: Vec<Option<f64>> = params.iter().map(|&p| f64_of(g, p)).collect();
    if !admissible(family, &numeric) {
        return Outcome::Pass;
    }
    let nu = numeric.first().copied().flatten();
    let Some(&first) = params.first() else {
        return Outcome::Pass;
    };
    let x = extra.first().copied().unwrap_or(first);
    let text = match request {
        | Request::Pdf => Some(family.pdf),
        | Request::Cdf => {
            if let Some(sum) = discrete_sum(g, family.name, &params, x) {
                return Outcome::Equal(sum);
            }
            Some(family.cdf)
        },
        // The moments of the Student distribution need enough degrees of
        // freedom.
        | Request::Mean if family.name == "student_t" && nu.is_some_and(|v| v <= 1.0) => None,
        | Request::Variance if family.name == "student_t" && nu.is_some_and(|v| v <= 2.0) => None,
        | Request::Mean => Some(family.mean),
        | Request::Variance => Some(family.variance),
        | Request::Mgf => family.mgf,
        | Request::Entropy => family.entropy,
    };
    let Some(text) = text else {
        return Outcome::Pass;
    };
    let Some(&a) = params.first() else {
        return Outcome::Pass;
    };
    let b = params.get(1).copied().unwrap_or(a);
    template(g, text, [a, b, x]).map_or(Outcome::Pass, Outcome::Equal)
}

// ---------------------------------------------------------------------
// Inference
// ---------------------------------------------------------------------

fn numbers_of(
    g: &Graph,
    list: NodeId,
) -> Option<Vec<f64>> {
    items_of(g, list)?.iter().map(|&x| f64_of(g, x)).collect()
}

fn float_pair(
    g: &mut Graph,
    a: f64,
    b: f64,
) -> Outcome {
    if a.is_nan() || b.is_nan() {
        return Outcome::Pass;
    }
    let first = g.float(a);
    let second = g.float(b);
    Outcome::Equal(list_node(g, &[first, second]))
}

fn sample_moments(v: &[f64]) -> Option<(f64, f64)> {
    Some((
        num::mean(v),
        num::variance_with_type(v, num::VarianceType::Sample)?,
    ))
}

/// Two-sided p-value of a Student `t` statistic.
fn t_p_value(
    t: f64,
    df: f64,
) -> Option<f64> {
    StudentsT::new(0.0, 1.0, df)
        .ok()
        .map(|d| 2.0 * d.cdf(-t.abs()))
}

/// `2 * (1 - cdf(dist(params), abs(statistic)))`.
fn two_sided_p(
    g: &mut Graph,
    dist: &str,
    params: &[NodeId],
    statistic: NodeId,
) -> Option<NodeId> {
    let d = call(g, dist, params)?;
    let magnitude = call(g, "abs", &[statistic])?;
    let cdf = call(g, "cdf", &[d, magnitude])?;
    let one = g.int(1);
    let tail = difference(g, one, cdf);
    Some(scale(g, tail, &Number::from(2)))
}

fn statistic_list(
    g: &mut Graph,
    statistic: NodeId,
    p_value: NodeId,
) -> Outcome {
    Outcome::Equal(list_node(g, &[statistic, p_value]))
}

fn z_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[sample, mu, sigma] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let Some(xs) = items_of(g, sample).filter(|x| !x.is_empty() && x.len() <= MAX_ITEMS) else {
        return Outcome::Pass;
    };
    if let (Some(v), Some(mu), Some(sigma)) =
        (numbers_of(g, sample), f64_of(g, mu), f64_of(g, sigma))
    {
        if sigma <= 0.0 {
            return Outcome::Pass;
        }
        let z = (num::mean(&v) - mu) * (count(v.len())).sqrt() / sigma;
        return float_pair(g, z, erfc(z.abs() / std::f64::consts::SQRT_2));
    }
    let build = |g: &mut Graph| -> Option<(NodeId, NodeId)> {
        let m = mean_term(g, &xs)?;
        let centred = difference(g, m, mu);
        let count = g.int(i64::try_from(xs.len()).ok()?);
        let root_n = pow_node(g, count, &Number::fraction(1, 2)?);
        let inverse_sigma = pow_node(g, sigma, &Number::from(-1));
        let z = product_node(g, &[centred, root_n, inverse_sigma]);
        let zero = g.int(0);
        let one = g.int(1);
        let p = two_sided_p(g, "normal", &[zero, one], z)?;
        Some((z, p))
    };
    match build(g) {
        | Some((z, p)) => statistic_list(g, z, p),
        | None => Outcome::Pass,
    }
}

fn t_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[sample, mu] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let Some(xs) = items_of(g, sample).filter(|x| x.len() >= 2 && x.len() <= MAX_ITEMS) else {
        return Outcome::Pass;
    };
    if let (Some(v), Some(mu)) = (numbers_of(g, sample), f64_of(g, mu)) {
        let Some((m, s2)) = sample_moments(&v).filter(|(_, s2)| *s2 > 0.0) else {
            return Outcome::Pass;
        };
        let t = (m - mu) / (s2 / count(v.len())).sqrt();
        return match t_p_value(t, count(v.len()) - 1.0) {
            | Some(p) => float_pair(g, t, p),
            | None => Outcome::Pass,
        };
    }
    let build = |g: &mut Graph| -> Option<(NodeId, NodeId)> {
        let m = mean_term(g, &xs)?;
        let centred = difference(g, m, mu);
        let s2 = sample_variance(g, &xs)?;
        let per_item = scale(g, s2, &reciprocal_of(xs.len())?);
        let spread = pow_node(g, per_item, &Number::fraction(-1, 2)?);
        let t = product_node(g, &[centred, spread]);
        let df = g.int(i64::try_from(xs.len()).ok()? - 1);
        let p = two_sided_p(g, "student_t", &[df], t)?;
        Some((t, p))
    };
    match build(g) {
        | Some((t, p)) => statistic_list(g, t, p),
        | None => Outcome::Pass,
    }
}

/// Welch's statistic and degrees of freedom from the two samples' means,
/// sample variances and sizes.
fn welch_numbers(
    a: &[f64],
    b: &[f64],
    mu_diff: f64,
) -> Option<(f64, f64)> {
    let ((m1, v1), (m2, v2)) = (sample_moments(a)?, sample_moments(b)?);
    let (n1, n2) = (count(a.len()), count(b.len()));
    let (t1, t2) = (v1 / n1, v2 / n2);
    let se = (t1 + t2).sqrt();
    if se == 0.0 {
        return None;
    }
    let df = (t1 + t2).powi(2) / (t1 * t1 / (n1 - 1.0) + t2 * t2 / (n2 - 1.0));
    Some(((m1 - m2 - mu_diff) / se, df))
}

fn welch_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[first, second, mu_diff] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys)) = (items_of(g, first), items_of(g, second)) else {
        return Outcome::Pass;
    };
    if xs.len() < 2 || ys.len() < 2 || xs.len() > MAX_ITEMS || ys.len() > MAX_ITEMS {
        return Outcome::Pass;
    }
    if let (Some(a), Some(b), Some(mu_diff)) = (
        numbers_of(g, first),
        numbers_of(g, second),
        f64_of(g, mu_diff),
    ) {
        if mu_diff == 0.0 {
            let (t, p) = num::welch_t_test(&a, &b);
            return float_pair(g, t, p);
        }
        let Some((t, df)) = welch_numbers(&a, &b, mu_diff) else {
            return Outcome::Pass;
        };
        return match t_p_value(t, df) {
            | Some(p) => float_pair(g, t, p),
            | None => Outcome::Pass,
        };
    }
    let build = |g: &mut Graph| -> Option<(NodeId, NodeId)> {
        let (m1, m2) = (mean_term(g, &xs)?, mean_term(g, &ys)?);
        let (s1, s2) = (sample_variance(g, &xs)?, sample_variance(g, &ys)?);
        let t1 = scale(g, s1, &reciprocal_of(xs.len())?);
        let t2 = scale(g, s2, &reciprocal_of(ys.len())?);
        let total = sum_node(g, &[t1, t2]);
        let centred = {
            let gap = difference(g, m1, m2);
            difference(g, gap, mu_diff)
        };
        let spread = pow_node(g, total, &Number::fraction(-1, 2)?);
        let t = product_node(g, &[centred, spread]);
        // Welch-Satterthwaite degrees of freedom.
        let top = pow_node(g, total, &Number::from(2));
        let (q1, q2) = (
            pow_node(g, t1, &Number::from(2)),
            pow_node(g, t2, &Number::from(2)),
        );
        let d1 = scale(g, q1, &reciprocal_of(xs.len() - 1)?);
        let d2 = scale(g, q2, &reciprocal_of(ys.len() - 1)?);
        let bottom = sum_node(g, &[d1, d2]);
        let inverse_bottom = pow_node(g, bottom, &Number::from(-1));
        let df = product_node(g, &[top, inverse_bottom]);
        let p = two_sided_p(g, "student_t", &[df], t)?;
        Some((t, p))
    };
    match build(g) {
        | Some((t, p)) => statistic_list(g, t, p),
        | None => Outcome::Pass,
    }
}

fn pooled_t_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[first, second] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys)) = (items_of(g, first), items_of(g, second)) else {
        return Outcome::Pass;
    };
    if xs.len() < 2 || ys.len() < 2 || xs.len() > MAX_ITEMS || ys.len() > MAX_ITEMS {
        return Outcome::Pass;
    }
    if let (Some(a), Some(b)) = (numbers_of(g, first), numbers_of(g, second)) {
        let (Some((m1, v1)), Some((m2, v2))) = (sample_moments(&a), sample_moments(&b)) else {
            return Outcome::Pass;
        };
        let (n1, n2) = (count(a.len()), count(b.len()));
        let pooled = (n2 - 1.0).mul_add(v2, (n1 - 1.0) * v1) / (n1 + n2 - 2.0);
        if pooled <= 0.0 {
            return Outcome::Pass;
        }
        let t = (m1 - m2) / (pooled * (1.0 / n1 + 1.0 / n2)).sqrt();
        return match t_p_value(t, n1 + n2 - 2.0) {
            | Some(p) => float_pair(g, t, p),
            | None => Outcome::Pass,
        };
    }
    let build = |g: &mut Graph| -> Option<(NodeId, NodeId)> {
        let (n1, n2) = (i64::try_from(xs.len()).ok()?, i64::try_from(ys.len()).ok()?);
        let (m1, m2) = (mean_term(g, &xs)?, mean_term(g, &ys)?);
        let (s1, s2) = (sample_variance(g, &xs)?, sample_variance(g, &ys)?);
        let w1 = scale(g, s1, &Number::from(n1 - 1));
        let w2 = scale(g, s2, &Number::from(n2 - 1));
        let numerator = sum_node(g, &[w1, w2]);
        let pooled = scale(g, numerator, &Number::fraction(1, n1 + n2 - 2)?);
        let sizes = g.num(Number::rat(
            BigRational::new(BigInt::one(), BigInt::from(n1))
                + BigRational::new(BigInt::one(), BigInt::from(n2)),
        ));
        let variance = product_node(g, &[pooled, sizes]);
        let spread = pow_node(g, variance, &Number::fraction(-1, 2)?);
        let gap = difference(g, m1, m2);
        let t = product_node(g, &[gap, spread]);
        let df = g.int(n1 + n2 - 2);
        let p = two_sided_p(g, "student_t", &[df], t)?;
        Some((t, p))
    };
    match build(g) {
        | Some((t, p)) => statistic_list(g, t, p),
        | None => Outcome::Pass,
    }
}

fn chi_squared_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[observed, expected] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(os), Some(es)) = (items_of(g, observed), items_of(g, expected)) else {
        return Outcome::Pass;
    };
    if os.len() < 2 || os.len() != es.len() || os.len() > MAX_ITEMS {
        return Outcome::Pass;
    }
    if let (Some(o), Some(e)) = (numbers_of(g, observed), numbers_of(g, expected)) {
        if e.iter().any(|&v| v <= 0.0) {
            return Outcome::Pass;
        }
        let (statistic, p) = num::chi_squared_test(&o, &e);
        return float_pair(g, statistic, p);
    }
    let build = |g: &mut Graph| -> Option<(NodeId, NodeId)> {
        let mut terms = Vec::with_capacity(os.len());
        for (&o, &e) in os.iter().zip(&es) {
            let gap = difference(g, o, e);
            let square = pow_node(g, gap, &Number::from(2));
            let inverse = pow_node(g, e, &Number::from(-1));
            terms.push(product_node(g, &[square, inverse]));
        }
        let statistic = sum_node(g, &terms);
        let df = g.int(i64::try_from(os.len()).ok()? - 1);
        let dist = call(g, "chi_squared", &[df])?;
        let cdf = call(g, "cdf", &[dist, statistic])?;
        let one = g.int(1);
        let p = difference(g, one, cdf);
        Some((statistic, p))
    };
    match build(g) {
        | Some((s, p)) => statistic_list(g, s, p),
        | None => Outcome::Pass,
    }
}

fn anova(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[groups] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let Some(members) = items_of(g, groups).filter(|m| m.len() >= 2) else {
        return Outcome::Pass;
    };
    let data: Option<Vec<Vec<f64>>> = members
        .iter()
        .map(|&m| numbers_of(g, m).filter(|v| !v.is_empty()))
        .collect();
    let Some(mut data) = data else {
        return Outcome::Pass;
    };
    if data.iter().map(Vec::len).sum::<usize>() <= data.len() {
        return Outcome::Pass;
    }
    let mut views: Vec<&mut [f64]> = data.iter_mut().map(Vec::as_mut_slice).collect();
    let (f, p) = num::one_way_anova(&mut views);
    float_pair(g, f, p)
}

/// The half-width `quantile * scale` of an interval around the mean of
/// numeric data, from a critical value computed by `critical`.
fn interval(
    g: &mut Graph,
    v: &[f64],
    half_width: f64,
) -> Outcome {
    if !half_width.is_finite() {
        return Outcome::Pass;
    }
    let m = num::mean(v);
    let lo = g.float(m - half_width);
    let hi = g.float(m + half_width);
    Outcome::Equal(list_node(g, &[lo, hi]))
}

fn confidence_interval(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[sample, level] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(v), Some(level)) = (numbers_of(g, sample), f64_of(g, level)) else {
        return Outcome::Pass;
    };
    if v.len() < 2 || !(level > 0.0 && level < 1.0) {
        return Outcome::Pass;
    }
    let Some((_, s2)) = sample_moments(&v) else {
        return Outcome::Pass;
    };
    let Ok(t) = StudentsT::new(0.0, 1.0, count(v.len()) - 1.0) else {
        return Outcome::Pass;
    };
    let critical = t.inverse_cdf(f64::midpoint(1.0, level));
    interval(g, &v, critical * (s2 / count(v.len())).sqrt())
}

fn z_interval(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[sample, sigma, level] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(v), Some(sigma), Some(level)) =
        (numbers_of(g, sample), f64_of(g, sigma), f64_of(g, level))
    else {
        return Outcome::Pass;
    };
    if v.is_empty() || sigma <= 0.0 || !(level > 0.0 && level < 1.0) {
        return Outcome::Pass;
    }
    let Ok(standard) = Normal::new(0.0, 1.0) else {
        return Outcome::Pass;
    };
    let critical = standard.inverse_cdf(f64::midpoint(1.0, level));
    interval(g, &v, critical * sigma / (count(v.len())).sqrt())
}

// ---------------------------------------------------------------------
// Regression
// ---------------------------------------------------------------------

/// Solves the square system `a x = b` over the rationals.
fn solve_rational(
    mut a: Vec<Vec<BigRational>>,
    mut b: Vec<BigRational>,
) -> Option<Vec<BigRational>> {
    let n = a.len();
    for col in 0..n {
        let pivot = (col..n).find(|&r| {
            a.get(r)
                .and_then(|row| row.get(col))
                .is_some_and(|v| !v.is_zero())
        })?;
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (pivot_row, pivot_rhs) = (a.get(col)?.clone(), b.get(col)?.clone());
        let lead = pivot_row.get(col)?.clone();
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = a.get(r)?.get(col)? / &lead;
            if factor.is_zero() {
                continue;
            }
            for (slot, above) in a.get_mut(r)?.iter_mut().zip(&pivot_row) {
                *slot -= &factor * above;
            }
            *b.get_mut(r)? -= &factor * &pivot_rhs;
        }
    }
    (0..n)
        .map(|r| Some(b.get(r)?.clone() / a.get(r)?.get(r)?))
        .collect()
}

/// Least-squares coefficients for the design matrix `rows` by the normal
/// equations.
fn least_squares(
    rows: &[Vec<BigRational>],
    y: &[BigRational],
) -> Option<Vec<BigRational>> {
    let k = rows.first()?.len();
    let mut gram = vec![vec![BigRational::zero(); k]; k];
    let mut moment = vec![BigRational::zero(); k];
    for (row, target) in rows.iter().zip(y) {
        for (i, xi) in row.iter().enumerate() {
            *moment.get_mut(i)? += xi * target;
            for (j, xj) in row.iter().enumerate() {
                *gram.get_mut(i)?.get_mut(j)? += xi * xj;
            }
        }
    }
    solve_rational(gram, moment)
}

/// The exact values of a list of numbers, and whether any was a float.
fn exact_column(
    g: &Graph,
    xs: &[NodeId],
) -> Option<(Vec<BigRational>, bool)> {
    let numbers = number_items(g, xs)?;
    let float = numbers.iter().any(|n| matches!(n, Number::Float(_)));
    Some((
        numbers.iter().map(rational_of).collect::<Option<_>>()?,
        float,
    ))
}

/// Coefficient list node from rational values.
fn coefficients_node(
    g: &mut Graph,
    values: &[BigRational],
    float: bool,
) -> NodeId {
    let nodes: Vec<NodeId> = values.iter().map(|v| value_node(g, v, float)).collect();
    list_node(g, &nodes)
}

fn simple_regression(
    g: &mut Graph,
    xs: &[NodeId],
    ys: &[NodeId],
) -> Option<NodeId> {
    let n = xs.len();
    if n < 2 || n != ys.len() {
        return None;
    }
    if let (Some(a), Some(b)) = (float_items(g, xs), float_items(g, ys)) {
        let pairs: Vec<(f64, f64)> = a.into_iter().zip(b).collect();
        let (slope, intercept) = num::simple_linear_regression(&pairs);
        if !slope.is_finite() || !intercept.is_finite() {
            return None;
        }
        let intercept = g.float(intercept);
        let slope = g.float(slope);
        return Some(list_node(g, &[intercept, slope]));
    }
    let variance = central_moment(g, xs, 2, n)?;
    if is_zero_number(g, variance) {
        return None;
    }
    let cov = covariance(g, xs, ys)?;
    let inverse = pow_node(g, variance, &Number::from(-1));
    let slope = product_node(g, &[cov, inverse]);
    let (mx, my) = (mean_term(g, xs)?, mean_term(g, ys)?);
    let drift = product_node(g, &[slope, mx]);
    let intercept = difference(g, my, drift);
    Some(list_node(g, &[intercept, slope]))
}

fn multiple_regression(
    g: &mut Graph,
    columns: &[NodeId],
    ys: &[NodeId],
) -> Option<NodeId> {
    let (y, y_float) = exact_column(g, ys)?;
    let mut data = Vec::new();
    let mut float = y_float;
    for &column in columns {
        let items = items_of(g, column).filter(|c| c.len() == ys.len())?;
        let (values, f) = exact_column(g, &items)?;
        float |= f;
        data.push(values);
    }
    if columns.is_empty() || ys.len() <= columns.len() {
        return None;
    }
    let rows: Vec<Vec<BigRational>> = (0..ys.len())
        .map(|i| {
            std::iter::once(BigRational::one())
                .chain(data.iter().filter_map(|c| c.get(i).cloned()))
                .collect()
        })
        .collect();
    let beta = least_squares(&rows, &y)?;
    Some(coefficients_node(g, &beta, float))
}

fn linear_regression(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[predictors, targets] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys)) = (items_of(g, predictors), items_of(g, targets)) else {
        return Outcome::Pass;
    };
    if xs.len() > MAX_ITEMS || ys.len() > MAX_ITEMS {
        return Outcome::Pass;
    }
    let result = if !xs.is_empty() && xs.iter().all(|&c| g.op(c) == core::LIST) {
        multiple_regression(g, &xs, &ys)
    } else {
        simple_regression(g, &xs, &ys)
    };
    finish(result)
}

fn polynomial_regression(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[predictors, targets, degree] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys)) = (items_of(g, predictors), items_of(g, targets)) else {
        return Outcome::Pass;
    };
    let Some(degree) = g
        .number_of(degree)
        .and_then(Number::to_i64)
        .and_then(|d| usize::try_from(d).ok())
    else {
        return Outcome::Pass;
    };
    if degree > MAX_DEGREE || xs.len() != ys.len() || xs.len() <= degree || xs.len() > MAX_ITEMS {
        return Outcome::Pass;
    }
    let (Some((x, x_float)), Some((y, y_float))) = (exact_column(g, &xs), exact_column(g, &ys))
    else {
        return Outcome::Pass;
    };
    let rows: Vec<Vec<BigRational>> = x
        .iter()
        .map(|v| {
            (0..=degree)
                .scan(BigRational::one(), |power, _| {
                    let current = power.clone();
                    *power *= v;
                    Some(current)
                })
                .collect()
        })
        .collect();
    let result = least_squares(&rows, &y).map(|c| coefficients_node(g, &c, x_float || y_float));
    finish(result)
}

/// `nonlinear_regression(xs, ys, model, x, list(params))`: the normal
/// equations of the least-squares fit of `model` (a function of the symbol
/// `x` and the parameters) as a `solve` request.
fn nonlinear_regression(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Outcome {
    let &[predictors, targets, model, x, parameters] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let (Some(xs), Some(ys), Some(params)) = (
        items_of(g, predictors),
        items_of(g, targets),
        items_of(g, parameters),
    ) else {
        return Outcome::Pass;
    };
    if xs.len() != ys.len() || xs.is_empty() || xs.len() > MAX_ITEMS || g.as_symbol(x).is_none() {
        return Outcome::Pass;
    }
    let Some(model) = best(g, model) else {
        return Outcome::Pass;
    };
    let build = |g: &mut Graph| -> Option<NodeId> {
        let mut residuals = Vec::with_capacity(xs.len());
        for (&xi, &yi) in xs.iter().zip(&ys) {
            let fitted = g.substitute(model, x, xi);
            let residual = difference(g, yi, fitted);
            residuals.push(pow_node(g, residual, &Number::from(2)));
        }
        let loss = sum_node(g, &residuals);
        let zero = g.int(0);
        let mut equations = Vec::with_capacity(params.len());
        for &p in &params {
            let slope = call(g, "diff", &[loss, p])?;
            equations.push(g.node(core::EQ, &[slope, zero]));
        }
        let system = list_node(g, &equations);
        call(g, "solve", &[system, parameters])
    };
    build(g).map_or(Outcome::Pass, Outcome::Equal)
}

// ---------------------------------------------------------------------
// Information theory
// ---------------------------------------------------------------------

/// `ln(n)` for a positive integer, through its prime factorisation.
fn ln_integer(
    g: &mut Graph,
    n: &BigInt,
) -> Option<NodeId> {
    if n.is_one() {
        return Some(g.int(0));
    }
    let Some(primes) = prime_factors(n) else {
        let node = g.num(Number::Int(n.clone()));
        return call(g, "ln", &[node]);
    };
    let mut terms = Vec::with_capacity(primes.len());
    for (p, e) in primes {
        let prime = g.num(Number::Int(BigInt::from(p)));
        let ln = call(g, "ln", &[prime])?;
        let exponent = g.int(i64::from(e));
        terms.push(product_node(g, &[exponent, ln]));
    }
    Some(sum_node(g, &terms))
}

/// `ln(p)`: exact numbers through their prime factorisation, floats
/// numerically, anything else as a `ln` node. `None` for `p <= 0`.
fn ln_of(
    g: &mut Graph,
    p: NodeId,
) -> Option<NodeId> {
    match g.number_of(p).cloned() {
        | Some(Number::Float(v)) => (v > 0.0).then(|| g.float(v.ln())),
        | Some(n) => {
            let r = n.to_rational().filter(Signed::is_positive)?;
            let (top, bottom) = (ln_integer(g, r.numer())?, ln_integer(g, r.denom())?);
            let minus_bottom = negate(g, bottom);
            Some(sum_node(g, &[top, minus_bottom]))
        },
        | None => call(g, "ln", &[p]),
    }
}

/// `p * (ln p - ln q)` (or `p * ln p` without `q`); `Ok(None)` for a
/// term that is zero by the convention `0 ln 0 = 0`, `Err` when the term
/// is undefined or infinite.
fn log_term(
    g: &mut Graph,
    p: NodeId,
    q: Option<NodeId>,
) -> Result<Option<NodeId>, ()> {
    if is_zero_number(g, p) {
        return Ok(None);
    }
    let ln_p = ln_of(g, p).ok_or(())?;
    let logarithm = match q {
        | Some(q) => {
            let ln_q = ln_of(g, q).ok_or(())?;
            difference(g, ln_p, ln_q)
        },
        | None => ln_p,
    };
    Ok(Some(product_node(g, &[p, logarithm])))
}

/// `-sum p ln p`.
fn entropy_node(
    g: &mut Graph,
    ps: &[NodeId],
) -> Option<NodeId> {
    if ps.is_empty() {
        return None;
    }
    let mut terms = Vec::with_capacity(ps.len());
    for &p in ps {
        terms.extend(log_term(g, p, None).ok()?);
    }
    let total = sum_node(g, &terms);
    Some(negate(g, total))
}

fn entropy(
    g: &mut Graph,
    ps: &[NodeId],
) -> Option<NodeId> {
    if let Some(v) = float_items(g, ps) {
        if v.iter().any(|p| *p < 0.0) {
            return None;
        }
        // The kernel works in bits.
        return Some(g.float(num::shannon_entropy(&v) * std::f64::consts::LN_2));
    }
    entropy_node(g, ps)
}

fn entropy_bits(
    g: &mut Graph,
    ps: &[NodeId],
) -> Option<NodeId> {
    if let Some(v) = float_items(g, ps) {
        return (v.iter().all(|p| *p >= 0.0)).then(|| g.float(num::shannon_entropy(&v)));
    }
    let nats = entropy_node(g, ps)?;
    let two = g.int(2);
    let ln_two = call(g, "ln", &[two])?;
    let inverse = pow_node(g, ln_two, &Number::from(-1));
    Some(product_node(g, &[nats, inverse]))
}

fn kl_divergence(
    g: &mut Graph,
    ps: &[NodeId],
    qs: &[NodeId],
) -> Option<NodeId> {
    if ps.is_empty() || ps.len() != qs.len() {
        return None;
    }
    let mut terms = Vec::with_capacity(ps.len());
    for (&p, &q) in ps.iter().zip(qs) {
        terms.extend(log_term(g, p, Some(q)).ok()?);
    }
    Some(sum_node(g, &terms))
}

fn cross_entropy(
    g: &mut Graph,
    ps: &[NodeId],
    qs: &[NodeId],
) -> Option<NodeId> {
    if ps.is_empty() || ps.len() != qs.len() {
        return None;
    }
    let mut terms = Vec::with_capacity(ps.len());
    for (&p, &q) in ps.iter().zip(qs) {
        if is_zero_number(g, p) {
            continue;
        }
        let ln_q = ln_of(g, q)?;
        terms.push(product_node(g, &[p, ln_q]));
    }
    let total = sum_node(g, &terms);
    Some(negate(g, total))
}

fn gini(
    g: &mut Graph,
    ps: &[NodeId],
) -> Option<NodeId> {
    if ps.is_empty() {
        return None;
    }
    let squares: Vec<NodeId> = ps
        .iter()
        .map(|&p| pow_node(g, p, &Number::from(2)))
        .collect();
    let total = sum_node(g, &squares);
    let one = g.int(1);
    Some(difference(g, one, total))
}

/// The rows of a joint probability matrix `list(list(...), ...)`.
fn matrix_of(
    g: &Graph,
    m: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let rows: Vec<Vec<NodeId>> = items_of(g, m)?
        .iter()
        .map(|&r| items_of(g, r))
        .collect::<Option<_>>()?;
    let width = rows.first()?.len();
    (width > 0 && rows.iter().all(|r| r.len() == width) && rows.len() * width <= MAX_ITEMS)
        .then_some(rows)
}

/// Entropies of the joint distribution and of its two marginals.
fn joint_parts(
    g: &mut Graph,
    m: NodeId,
) -> Option<(NodeId, NodeId, NodeId)> {
    let rows = matrix_of(g, m)?;
    let flat: Vec<NodeId> = rows.iter().flatten().copied().collect();
    let by_row: Vec<NodeId> = rows.iter().map(|r| sum_node(g, r)).collect();
    let width = rows.first()?.len();
    let by_column: Vec<NodeId> = (0..width)
        .map(|j| {
            let column: Vec<NodeId> = rows.iter().filter_map(|r| r.get(j).copied()).collect();
            sum_node(g, &column)
        })
        .collect();
    Some((
        entropy_node(g, &flat)?,
        entropy_node(g, &by_row)?,
        entropy_node(g, &by_column)?,
    ))
}

fn joint_entropy(
    g: &mut Graph,
    m: NodeId,
) -> Option<NodeId> {
    joint_parts(g, m).map(|(xy, _, _)| xy)
}

fn conditional_entropy(
    g: &mut Graph,
    m: NodeId,
) -> Option<NodeId> {
    let (xy, x, _) = joint_parts(g, m)?;
    Some(difference(g, xy, x))
}

fn mutual_information(
    g: &mut Graph,
    m: NodeId,
) -> Option<NodeId> {
    let (xy, x, y) = joint_parts(g, m)?;
    let marginals = sum_node(g, &[x, y]);
    Some(difference(g, marginals, xy))
}

/// Wraps a function of a matrix into a kernel body.
fn matrix_op(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    f: fn(&mut Graph, NodeId) -> Option<NodeId>,
) -> Outcome {
    let &[m] = args else {
        return Outcome::Pass;
    };
    let g = &mut *cx.graph;
    let result = f(g, m);
    finish(result)
}

// ---------------------------------------------------------------------
// Installation
// ---------------------------------------------------------------------

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let set = "stats";
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));
    let ternary = |name: &str| OpDescriptor::new(name, Arity::Fixed(3));
    let request = |desc: OpDescriptor| desc.flags(OpFlags::HEAVY).cost(100);

    // Distributions: inert terms and the requests on them.
    for family in &FAMILIES {
        i.op(OpDescriptor::new(family.name, Arity::Fixed(family.params)))?;
    }
    procedure(i, set, request(binary("pdf")), Tier::Reduce, |cx, a| {
        distribution_request(cx, a, Request::Pdf)
    })?;
    procedure(i, set, request(binary("cdf")), Tier::Reduce, |cx, a| {
        distribution_request(cx, a, Request::Cdf)
    })?;
    procedure(
        i,
        set,
        request(unary("expectation")),
        Tier::Reduce,
        |cx, a| distribution_request(cx, a, Request::Mean),
    )?;
    procedure(
        i,
        set,
        request(unary("variance_of")),
        Tier::Reduce,
        |cx, a| distribution_request(cx, a, Request::Variance),
    )?;
    procedure(i, set, request(binary("mgf")), Tier::Reduce, |cx, a| {
        distribution_request(cx, a, Request::Mgf)
    })?;
    procedure(
        i,
        set,
        request(unary("entropy_of")),
        Tier::Reduce,
        |cx, a| distribution_request(cx, a, Request::Entropy),
    )?;

    // Regularised incomplete gamma and beta functions for the cdfs.
    i.op(binary("gamma_lr").eval(gamma_lr_eval))?;
    i.op(ternary("beta_reg").eval(beta_reg_eval))?;
    // Only the derivative in the last argument is known; `?c` and `?d` name
    // arguments that do not exist, which marks the others as unavailable.
    partials(i, "gamma_lr", &["?c", "?b^(?a - 1) * exp(-?b) / gamma(?a)"])?;
    partials(
        i,
        "beta_reg",
        &["?d", "?d", "?c^(?a - 1) * (1 - ?c)^(?b - 1) / beta(?a, ?b)"],
    )?;
    i.rewrites(
        Tier::Normalize,
        &[
            "stats/gamma-lr-0: gamma_lr(?a, 0) => 0 if positive(?a)",
            "stats/gamma-lr-exp: gamma_lr(1, ?x) => 1 - exp(-?x) if nonnegative(?x)",
            "stats/beta-reg-0: beta_reg(?a, ?b, 0) => 0 if positive(?a), positive(?b)",
            "stats/beta-reg-1: beta_reg(?a, ?b, 1) => 1 if positive(?a), positive(?b)",
        ],
    )?;

    // Descriptive statistics.
    let lists: [(&str, ListFn); 17] = [
        ("mean", mean),
        ("variance", variance),
        ("sample_variance", sample_variance),
        ("std", std_dev),
        ("sample_std", sample_std),
        ("median", median),
        ("mode", mode),
        ("skewness", skewness),
        ("kurtosis", kurtosis),
        ("min_of", min_of),
        ("max_of", max_of),
        ("range_of", range_of),
        ("geometric_mean", geometric_mean),
        ("harmonic_mean", harmonic_mean),
        ("iqr", iqr),
        ("zscores", zscores),
        ("coefficient_of_variation", coefficient_of_variation),
    ];
    for (name, f) in lists {
        register_list(i, name, f)?;
    }
    register_list(i, "standard_error", standard_error)?;
    register_list(i, "entropy", entropy)?;
    register_list(i, "entropy_bits", entropy_bits)?;
    register_list(i, "gini", gini)?;
    register_pair(i, "covariance", covariance)?;
    register_pair(i, "correlation", correlation)?;
    register_pair(i, "kl_divergence", kl_divergence)?;
    register_pair(i, "cross_entropy", cross_entropy)?;
    register_list_with(i, "moment", moment)?;
    register_list_with(i, "quantile", quantile)?;
    register_list_with(i, "percentile", percentile)?;

    // Information theory on joint distributions.
    procedure(
        i,
        set,
        request(unary("joint_entropy")),
        Tier::Normalize,
        |cx, a| matrix_op(cx, a, joint_entropy),
    )?;
    procedure(
        i,
        set,
        request(unary("conditional_entropy")),
        Tier::Normalize,
        |cx, a| matrix_op(cx, a, conditional_entropy),
    )?;
    procedure(
        i,
        set,
        request(unary("mutual_information")),
        Tier::Normalize,
        |cx, a| matrix_op(cx, a, mutual_information),
    )?;

    // Inference and regression are requests.
    procedure(i, set, request(ternary("z_test")), Tier::Reduce, z_test)?;
    procedure(i, set, request(binary("t_test")), Tier::Reduce, t_test)?;
    procedure(
        i,
        set,
        request(ternary("welch_test")),
        Tier::Reduce,
        welch_test,
    )?;
    procedure(
        i,
        set,
        request(binary("pooled_t_test")),
        Tier::Reduce,
        pooled_t_test,
    )?;
    procedure(
        i,
        set,
        request(binary("chi_squared_test")),
        Tier::Reduce,
        chi_squared_test,
    )?;
    procedure(i, set, request(unary("anova")), Tier::Reduce, anova)?;
    procedure(
        i,
        set,
        request(binary("confidence_interval")),
        Tier::Reduce,
        confidence_interval,
    )?;
    procedure(
        i,
        set,
        request(ternary("z_interval")),
        Tier::Reduce,
        z_interval,
    )?;
    procedure(
        i,
        set,
        request(binary("linear_regression")),
        Tier::Reduce,
        linear_regression,
    )?;
    procedure(
        i,
        set,
        request(ternary("polynomial_regression")),
        Tier::Reduce,
        polynomial_regression,
    )?;
    procedure(
        i,
        set,
        request(OpDescriptor::new("nonlinear_regression", Arity::Fixed(5))),
        Tier::Reduce,
        nonlinear_regression,
    )?;
    Ok(())
}

/// A kernel for one operator that runs a closure.
struct Wrapped<F> {
    op: OpId,
    run: F,
}

impl<F: Fn(&mut Cx<'_>, &[NodeId]) -> Outcome + Send + Sync> Kernel for Wrapped<F> {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        (self.run)(cx, &args)
    }

    fn revisit(&self) -> bool {
        true
    }
}

fn register_list(
    i: &mut Installer<'_>,
    name: &str,
    f: ListFn,
) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new(name, Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    i.kernel(
        &format!("stats/{name}"),
        Tier::Normalize,
        Wrapped {
            op,
            run: move |cx: &mut Cx<'_>, a: &[NodeId]| list_op(cx, a, f),
        },
    );
    Ok(())
}

fn register_pair(
    i: &mut Installer<'_>,
    name: &str,
    f: PairFn,
) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new(name, Arity::Fixed(2))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    i.kernel(
        &format!("stats/{name}"),
        Tier::Normalize,
        Wrapped {
            op,
            run: move |cx: &mut Cx<'_>, a: &[NodeId]| pair_op(cx, a, f),
        },
    );
    Ok(())
}

fn register_list_with(
    i: &mut Installer<'_>,
    name: &str,
    f: ScalarFn,
) -> Result<(), RuleError> {
    let op = i.op(OpDescriptor::new(name, Arity::Fixed(2))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    i.kernel(
        &format!("stats/{name}"),
        Tier::Normalize,
        Wrapped {
            op,
            run: move |cx: &mut Cx<'_>, a: &[NodeId]| list_with_op(cx, a, f),
        },
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use statrs::distribution::Bernoulli;
    use statrs::distribution::Beta;
    use statrs::distribution::Binomial;
    use statrs::distribution::ChiSquared;
    use statrs::distribution::Continuous;
    use statrs::distribution::Discrete;
    use statrs::distribution::DiscreteCDF;
    use statrs::distribution::Exp;
    use statrs::distribution::Gamma;
    use statrs::distribution::Poisson;
    use statrs::distribution::Uniform;
    use statrs::statistics::Distribution;

    use super::*;
    use crate::graph::Facts;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn sets() -> Vec<RuleSet> {
        vec![stats()]
    }

    fn s(src: &str) -> String {
        simplify(&sets(), src)
    }

    /// The numeric value of a term (no free symbols).
    fn at(
        text: &str,
        bindings: &[(&str, f64)],
    ) -> f64 {
        eval(&sets(), text, bindings)
    }

    fn value(src: &str) -> f64 {
        at(&s(src), &[])
    }

    fn close(
        got: f64,
        want: f64,
        tolerance: f64,
        what: &str,
    ) {
        assert!(
            (got - want).abs() <= tolerance * (1.0 + want.abs()),
            "{what}: got {got}, want {want}"
        );
    }

    /// Splits `list(a, b, ...)` into its top-level items.
    fn items(text: &str) -> Vec<String> {
        let inner = text
            .strip_prefix("list(")
            .and_then(|t| t.strip_suffix(')'))
            .unwrap_or("");
        let (mut out, mut depth, mut current) = (Vec::new(), 0_i32, String::new());
        for ch in inner.chars() {
            match ch {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | _ => {},
            }
            if ch == ',' && depth == 0 {
                out.push(current.trim().to_owned());
                current.clear();
            } else {
                current.push(ch);
            }
        }
        if !current.trim().is_empty() {
            out.push(current.trim().to_owned());
        }
        out
    }

    fn stays(src: &str) {
        let (text, reduced) = reduce_with(&sets(), src, &[]);
        assert!(!reduced, "`{src}` should stay a request, got {text}");
    }

    #[test]
    fn exact_descriptive_statistics() {
        assert_eq!(s("mean(list(1, 2, 4))"), "7/3");
        assert_eq!(s("mean(list(5))"), "5");
        assert_eq!(s("variance(list(1, 2, 4))"), "14/9");
        assert_eq!(s("variance(list(3, 3, 3))"), "0");
        assert_eq!(s("sample_variance(list(1, 2, 4))"), "7/3");
        assert_eq!(s("sample_variance(list(2, 4, 4, 4, 5, 5, 7, 9))"), "32/7");
        // The classic example with population standard deviation 2.
        assert_eq!(s("variance(list(2, 4, 4, 4, 5, 5, 7, 9))"), "4");
        assert_eq!(s("std(list(2, 4, 4, 4, 5, 5, 7, 9))"), "2");
        assert_eq!(s("std(list(1, 2, 4))"), "(14/9)^(1/2)");
        assert_eq!(s("sample_std(list(1, 3))"), "2^(1/2)");
        assert_eq!(s("median(list(3, 1, 2))"), "2");
        assert_eq!(s("median(list(4, 1, 3, 2))"), "5/2");
        assert_eq!(s("median(list(7))"), "7");
        assert_eq!(s("median(list(-1, 1/2, 3))"), "1/2");
        assert_eq!(s("mode(list(1, 2, 2, 3))"), "2");
        // Equally frequent numbers: the smallest one.
        assert_eq!(s("mode(list(5, 3, 5, 3, 9))"), "3");
        assert_eq!(s("mode(list(4))"), "4");
        assert_eq!(s("covariance(list(1, 2, 3), list(2, 4, 7))"), "5/3");
        assert_eq!(s("covariance(list(1, 2, 3), list(3, 2, 1))"), "-2/3");
        assert_eq!(s("covariance(list(1, 2, 3), list(5, 5, 5))"), "0");
        assert_eq!(s("correlation(list(1, 2, 3), list(2, 4, 6))"), "1");
        assert_eq!(s("correlation(list(1, 2, 3), list(6, 4, 2))"), "-1");
        close(
            value("correlation(list(1, 2, 3), list(2, 4, 7))"),
            5.0 / 3.0 / (76.0_f64 / 27.0).sqrt(),
            1e-14,
            "corr",
        );
        assert_eq!(s("moment(list(1, 2, 3, 4), 0)"), "1");
        assert_eq!(s("moment(list(1, 2, 3, 4), 1)"), "0");
        assert_eq!(s("moment(list(1, 2, 3, 4), 2)"), "5/4");
        assert_eq!(s("moment(list(1, 2, 3, 4), 3)"), "0");
        assert_eq!(s("moment(list(1, 2, 3, 4), 4)"), "41/16");
        assert_eq!(s("skewness(list(1, 2, 3, 4))"), "0");
        assert_eq!(s("kurtosis(list(1, 2, 3, 4))"), "41/25");
        // Skewness of (1, 2, 10) from its definition.
        let (m, n) = (13.0_f64 / 3.0, 3.0);
        let m2 = ([1.0, 2.0, 10.0]
            .iter()
            .map(|x| (x - m).powi(2))
            .sum::<f64>())
            / n;
        let m3 = ([1.0, 2.0, 10.0]
            .iter()
            .map(|x| (x - m).powi(3))
            .sum::<f64>())
            / n;
        close(
            value("skewness(list(1, 2, 10))"),
            m3 / m2.powf(1.5),
            1e-13,
            "skewness",
        );
        // Order statistics.
        assert_eq!(s("min_of(list(3, 1, 2))"), "1");
        assert_eq!(s("max_of(list(3, 1, 2))"), "3");
        assert_eq!(s("range_of(list(3, 1, 2))"), "2");
        assert_eq!(s("quantile(list(1, 2, 3, 4, 5), 0)"), "1");
        assert_eq!(s("quantile(list(1, 2, 3, 4, 5), 1)"), "5");
        assert_eq!(s("quantile(list(1, 2, 3, 4, 5), 1/2)"), "3");
        assert_eq!(s("quantile(list(1, 2, 3, 4, 5), 1/4)"), "2");
        assert_eq!(s("quantile(list(1, 2, 3, 4, 5), 1/3)"), "7/3");
        assert_eq!(s("quantile(list(5, 1, 4, 2, 3), 0.5)"), "3");
        assert_eq!(s("quantile(list(1, 2, 3, 4), 1/2)"), "5/2");
        assert_eq!(s("percentile(list(1, 2, 3, 4), 50)"), "5/2");
        assert_eq!(s("percentile(list(1, 2, 3, 4), 100)"), "4");
        assert_eq!(s("iqr(list(1, 2, 3, 4, 5, 6, 7, 8))"), "7/2");
        // Other means and scales.
        assert_eq!(s("geometric_mean(list(1, 4, 16))"), "4");
        assert_eq!(s("geometric_mean(list(2, 8))"), "4");
        assert_eq!(s("geometric_mean(list(2, 3))"), "6^(1/2)");
        assert_eq!(s("harmonic_mean(list(1, 2, 4))"), "12/7");
        assert_eq!(s("coefficient_of_variation(list(1, 2, 3))"), "1/2");
        assert_eq!(
            s("standard_error(list(2, 4, 4, 4, 5, 5, 7, 9))"),
            "(4/7)^(1/2)"
        );
        assert_eq!(
            s("zscores(list(2, 4, 4, 4, 5, 5, 7, 9))"),
            "list(-3/2, -1/2, -1/2, -1/2, 0, 0, 1, 2)"
        );
        assert_eq!(s("zscores(list(3, 3, 3))"), "list(0, 0, 0)");
    }

    #[test]
    fn descriptive_statistics_refuse_what_they_cannot_do() {
        for src in [
            "mean(list())",
            "variance(list())",
            "sample_variance(list(1))",
            "median(list())",
            "mode(list(1, 2, 3))",
            "covariance(list(1, 2), list(1, 2, 3))",
            "correlation(list(1, 2, 3), list(4, 4, 4))",
            "skewness(list(2, 2, 2))",
            "moment(list(1, 2), -1)",
            "moment(list(1, 2), k)",
            "quantile(list(1, 2, 3), 3/2)",
            "quantile(list(1, 2, 3), p)",
            "geometric_mean(list(1, -4))",
            "harmonic_mean(list(1, 0))",
            "coefficient_of_variation(list(-1, 1))",
            "mean(x)",
            "std(3)",
        ] {
            stays(src);
        }
    }

    #[test]
    fn float_lists_use_the_numeric_kernels() {
        let data = [2.5, -1.25, 7.0, 3.5, 3.5, 0.125, 9.75];
        let text = format!(
            "list({})",
            data.iter()
                .map(|v| format!("{v:?}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        close(
            value(&format!("mean({text})")),
            num::mean(&data),
            1e-14,
            "mean",
        );
        close(
            value(&format!("variance({text})")),
            num::variance(&data),
            1e-13,
            "variance",
        );
        let sample = num::variance_with_type(&data, num::VarianceType::Sample).unwrap_or(f64::NAN);
        close(
            value(&format!("sample_variance({text})")),
            sample,
            1e-13,
            "sample variance",
        );
        close(
            value(&format!("std({text})")),
            num::variance(&data).sqrt(),
            1e-13,
            "std",
        );
        let mut sorted = data;
        close(
            value(&format!("median({text})")),
            num::median(&mut sorted),
            1e-14,
            "median",
        );
        close(
            value(&format!("min_of({text})")),
            num::min(&mut sorted),
            0.0,
            "min",
        );
        close(
            value(&format!("max_of({text})")),
            num::max(&mut sorted),
            0.0,
            "max",
        );
        close(
            value(&format!("range_of({text})")),
            num::range(&data),
            1e-14,
            "range",
        );
        let other = [1.0, 0.5, 2.0, 2.5, -3.0, 1.0, 4.0];
        let other_text = format!(
            "list({})",
            other
                .iter()
                .map(|v| format!("{v:?}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        // The kernel's covariance divides by n - 1.
        let population = num::covariance(&data, &other) * 6.0 / 7.0;
        close(
            value(&format!("covariance({text}, {other_text})")),
            population,
            1e-13,
            "covariance",
        );
        close(
            value(&format!("correlation({text}, {other_text})")),
            num::correlation(&data, &other),
            1e-12,
            "correlation",
        );
        // Mixed exact and float items count as floats.
        close(
            value("mean(list(1, 2.5, 4))"),
            7.5 / 3.0,
            1e-14,
            "mixed mean",
        );
        // Positive-data means agree with the kernels.
        let positive = [1.5, 2.0, 8.0, 0.25];
        let positive_text = "list(1.5, 2.0, 8.0, 0.25)";
        close(
            value(&format!("geometric_mean({positive_text})")),
            num::geometric_mean(&positive),
            1e-13,
            "geometric",
        );
        close(
            value(&format!("harmonic_mean({positive_text})")),
            num::harmonic_mean(&positive),
            1e-13,
            "harmonic",
        );
        close(
            value(&format!("standard_error({positive_text})")),
            num::variance_with_type(&positive, num::VarianceType::Sample)
                .unwrap_or(f64::NAN)
                .sqrt()
                / 2.0,
            1e-13,
            "se",
        );
        let z = items(&s(&format!("zscores({positive_text})")));
        let m = num::mean(&positive);
        let sd = num::variance(&positive).sqrt();
        for (item, x) in z.iter().zip(positive) {
            close(at(item, &[]), (x - m) / sd, 1e-13, "z-score");
        }
        close(
            value(&format!("coefficient_of_variation({positive_text})")),
            num::variance_with_type(&positive, num::VarianceType::Sample)
                .unwrap_or(f64::NAN)
                .sqrt()
                / m,
            1e-13,
            "cv",
        );
        // The population standard deviation differs from the kernel's
        // sample standard deviation exactly by sqrt((n - 1) / n).
        close(
            num::std_dev(&positive) * (3.0_f64 / 4.0).sqrt(),
            num::variance(&positive).sqrt(),
            1e-13,
            "kernel std",
        );
        // The kernel and the operator agree on the mode of repeated data.
        assert_eq!(num::mode(&[1.0, 2.0, 2.0, 3.0], 3), Some(2.0));
        assert_eq!(s("mode(list(1.0, 2.0, 2.0, 3.0))"), "2");
    }

    #[test]
    fn symbolic_lists_give_defining_formulas() {
        let a = [("a", 1.0), ("b", 2.0), ("c", 4.0)];
        let v = |src: &str| at(&s(src), &a);
        close(v("mean(list(a, b, c))"), 7.0 / 3.0, 1e-14, "mean");
        close(v("variance(list(a, b, c))"), 14.0 / 9.0, 1e-14, "variance");
        close(
            v("sample_variance(list(a, b, c))"),
            7.0 / 3.0,
            1e-14,
            "sample variance",
        );
        close(
            v("std(list(a, b, c))"),
            (14.0_f64 / 9.0).sqrt(),
            1e-14,
            "std",
        );
        close(
            v("moment(list(a, b, c), 3)"),
            ([1.0_f64, 2.0, 4.0]
                .iter()
                .map(|x| (x - 7.0 / 3.0).powi(3))
                .sum::<f64>())
                / 3.0,
            1e-13,
            "moment",
        );
        close(
            v("covariance(list(a, b, c), list(c, b, a))"),
            -13.0 / 9.0,
            1e-13,
            "covariance",
        );
        close(
            v("correlation(list(a, b, c), list(c, b, a))"),
            (-13.0 / 9.0) / (14.0_f64 / 9.0),
            1e-13,
            "correlation",
        );
        close(
            v("skewness(list(a, b, c))"),
            {
                let m3 = [1.0_f64, 2.0, 4.0]
                    .iter()
                    .map(|x| (x - 7.0 / 3.0).powi(3))
                    .sum::<f64>()
                    / 3.0;
                m3 / (14.0_f64 / 9.0).powf(1.5)
            },
            1e-13,
            "skewness",
        );
        close(
            v("geometric_mean(list(a, b, c))"),
            2.0,
            1e-14,
            "geometric mean",
        );
        close(
            v("harmonic_mean(list(a, b, c))"),
            12.0 / 7.0,
            1e-14,
            "harmonic mean",
        );
        close(
            v("standard_error(list(a, b, c))"),
            (7.0_f64 / 9.0).sqrt(),
            1e-14,
            "standard error",
        );
        close(
            v("coefficient_of_variation(list(a, b, c))"),
            (7.0_f64 / 3.0).sqrt() / (7.0 / 3.0),
            1e-14,
            "cv",
        );
        // Only the symbols that occur in the data remain in the formulas.
        assert_eq!(s("mean(list(a, b))"), "1/2*(a + b)");
        assert_eq!(s("mean(list(a, 3, 5))"), "1/3*(a + 8)");
        assert_eq!(s("mean(list(a, a, a))"), "a");
        // Symbolic data is generic: the formula holds for other values too.
        let other = [("a", -3.5), ("b", 0.25), ("c", 8.0)];
        close(
            at(&s("variance(list(a, b, c))"), &other),
            {
                let xs = [-3.5_f64, 0.25, 8.0];
                let m = xs.iter().sum::<f64>() / 3.0;
                xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / 3.0
            },
            1e-13,
            "variance at other values",
        );
        // Without numbers there is no order: only trivial medians and modes.
        assert_eq!(s("median(list(a))"), "a");
        assert_eq!(s("median(list(a, b))"), "1/2*(a + b)");
        assert_eq!(s("mode(list(a, b, a))"), "a");
        for src in [
            "median(list(a, b, c))",
            "quantile(list(a, b, c), 1/2)",
            "min_of(list(a, b))",
            "iqr(list(a, b, c, d))",
        ] {
            stays(src);
        }
    }

    /// Simpson's rule on 20000 intervals.
    fn integrate(
        f: &dyn Fn(f64) -> f64,
        lo: f64,
        hi: f64,
    ) -> f64 {
        let n = 20_000;
        let h = (hi - lo) / f64::from(n);
        let mut total = f(lo) + f(hi);
        for k in 1..n {
            total += f(lo + f64::from(k) * h) * if k % 2 == 1 { 4.0 } else { 2.0 };
        }
        total * h / 3.0
    }

    /// Checks a continuous distribution against `statrs` at `points`.
    fn check_continuous(
        dist: &str,
        pdf: &dyn Fn(f64) -> f64,
        cdf: &dyn Fn(f64) -> f64,
        points: &[f64],
    ) {
        let p = s(&format!("pdf({dist}, x)"));
        let c = s(&format!("cdf({dist}, x)"));
        for &x in points {
            close(
                at(&p, &[("x", x)]),
                pdf(x),
                1e-9,
                &format!("pdf of {dist} at {x}"),
            );
            close(
                at(&c, &[("x", x)]),
                cdf(x),
                1e-9,
                &format!("cdf of {dist} at {x}"),
            );
        }
    }

    /// Checks mean, variance and entropy against `statrs`.
    fn check_moments(
        dist: &str,
        mean: Option<f64>,
        variance: Option<f64>,
        entropy: Option<f64>,
    ) {
        if let Some(m) = mean {
            close(
                value(&format!("expectation({dist})")),
                m,
                1e-9,
                &format!("mean of {dist}"),
            );
        }
        if let Some(v) = variance {
            close(
                value(&format!("variance_of({dist})")),
                v,
                1e-9,
                &format!("variance of {dist}"),
            );
        }
        if let Some(h) = entropy {
            close(
                value(&format!("entropy_of({dist})")),
                h,
                1e-9,
                &format!("entropy of {dist}"),
            );
        }
    }

    #[test]
    fn normal_distribution_against_statrs() {
        for (mu, sigma) in [(0.0, 1.0), (1.5, 2.0), (-3.0, 0.25)] {
            let d = Normal::new(mu, sigma).unwrap_or_else(|e| panic!("{e}"));
            let points: Vec<f64> = [-4.0, -1.5, -0.3, 0.0, 0.7, 2.2, 5.0]
                .iter()
                .map(|z| mu + z * sigma)
                .collect();
            check_continuous(
                &format!("normal({mu:?}, {sigma:?})"),
                &|x| d.pdf(x),
                &|x| d.cdf(x),
                &points,
            );
            check_moments(
                &format!("normal({mu:?}, {sigma:?})"),
                d.mean(),
                d.variance(),
                d.entropy(),
            );
        }
        // Symbolic parameters keep their names.
        assert_eq!(s("expectation(normal(mu, sigma))"), "mu");
        assert_eq!(s("variance_of(normal(mu, sigma))"), "sigma^2");
        assert_eq!(s("cdf(normal(0, 1), 0)"), "1/2");
        assert_eq!(s("cdf(normal(mu, sigma), mu)"), "1/2");
    }

    #[test]
    fn uniform_and_exponential_against_statrs() {
        let u = Uniform::new(-1.0, 3.0).unwrap_or_else(|e| panic!("{e}"));
        check_continuous(
            "uniform(-1.0, 3.0)",
            &|x| u.pdf(x),
            &|x| ContinuousCDF::cdf(&u, x),
            &[-2.0, -0.5, 0.0, 1.3, 2.9, 4.0],
        );
        check_moments("uniform(-1.0, 3.0)", u.mean(), u.variance(), u.entropy());
        assert_eq!(s("cdf(uniform(0, 4), 1)"), "1/4");
        assert_eq!(s("cdf(uniform(0, 4), 7)"), "1");
        assert_eq!(s("cdf(uniform(0, 4), -7)"), "0");
        assert_eq!(s("expectation(uniform(a, b))"), "1/2*(a + b)");
        for rate in [0.4, 1.7, 12.0] {
            let e = Exp::new(rate).unwrap_or_else(|e| panic!("{e}"));
            let points: Vec<f64> = [-1.0, 0.05, 0.3, 1.0, 3.5]
                .iter()
                .map(|z| z / rate)
                .collect();
            check_continuous(
                &format!("exponential({rate:?})"),
                &|x| e.pdf(x),
                &|x| ContinuousCDF::cdf(&e, x),
                &points,
            );
            check_moments(
                &format!("exponential({rate:?})"),
                e.mean(),
                e.variance(),
                e.entropy(),
            );
        }
        assert_eq!(s("expectation(exponential(2))"), "1/2");
        assert_eq!(s("variance_of(exponential(2))"), "1/4");
    }

    #[test]
    fn gamma_beta_chi_squared_and_student_against_statrs() {
        for (shape, rate) in [(1.0, 1.0), (2.5, 1.3), (7.0, 0.5)] {
            let d = Gamma::new(shape, rate).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("gamma_dist({shape:?}, {rate:?})");
            let points: Vec<f64> = [-1.0, 0.3, 1.0, 2.5, 7.0, 20.0]
                .iter()
                .map(|z| z / rate)
                .collect();
            check_continuous(
                &name,
                &|x| d.pdf(x),
                &|x| ContinuousCDF::cdf(&d, x),
                &points,
            );
            check_moments(&name, d.mean(), d.variance(), d.entropy());
        }
        for (a, b) in [(1.0, 1.0), (2.0, 3.5), (0.9, 4.0), (6.0, 2.5)] {
            let d = Beta::new(a, b).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("beta_dist({a:?}, {b:?})");
            check_continuous(
                &name,
                &|x| d.pdf(x),
                &|x| ContinuousCDF::cdf(&d, x),
                &[-0.2, 0.05, 0.1, 0.4, 0.75, 0.95, 1.3],
            );
            check_moments(&name, d.mean(), d.variance(), d.entropy());
        }
        for k in [1.0, 3.0, 4.5, 10.0] {
            let d = ChiSquared::new(k).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("chi_squared({k:?})");
            check_continuous(
                &name,
                &|x| d.pdf(x),
                &|x| ContinuousCDF::cdf(&d, x),
                &[-1.0, 0.5, 2.0, 6.0, 15.0, 30.0],
            );
            check_moments(&name, d.mean(), d.variance(), d.entropy());
        }
        for nu in [1.0, 3.7, 12.0, 40.0] {
            let d = StudentsT::new(0.0, 1.0, nu).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("student_t({nu:?})");
            check_continuous(
                &name,
                &|x| d.pdf(x),
                &|x| d.cdf(x),
                &[-30.0, -3.0, -0.4, 0.0, 1.2, 2.0, 5.0],
            );
            // Moments exist only for enough degrees of freedom.
            check_moments(
                &name,
                (nu > 1.0).then_some(0.0),
                (nu > 2.0).then(|| nu / (nu - 2.0)),
                d.entropy(),
            );
        }
        assert_eq!(s("expectation(gamma_dist(a, b))"), "a/b");
        assert_eq!(s("expectation(beta_dist(2, 3))"), "2/5");
        assert_eq!(s("variance_of(beta_dist(2, 3))"), "1/25");
        assert_eq!(s("expectation(chi_squared(k))"), "k");
        assert_eq!(s("variance_of(student_t(5))"), "5/3");
        assert_eq!(s("cdf(student_t(nu), 0)"), "1/2");
        assert_eq!(s("cdf(chi_squared(2), x)"), "gamma_lr(1, 1/2*x)");
    }

    #[test]
    fn discrete_distributions_against_statrs() {
        // Bernoulli
        let b = Bernoulli::new(0.3).unwrap_or_else(|e| panic!("{e}"));
        let (pmf, cdf) = (s("pdf(bernoulli(0.3), x)"), s("cdf(bernoulli(0.3), x)"));
        for k in 0..=1_u64 {
            close(
                at(&pmf, &[("x", k as f64)]),
                b.pmf(k),
                1e-12,
                "bernoulli pmf",
            );
            close(
                at(&cdf, &[("x", k as f64)]),
                DiscreteCDF::cdf(&b, k),
                1e-12,
                "bernoulli cdf",
            );
        }
        close(at(&cdf, &[("x", -1.0)]), 0.0, 1e-12, "bernoulli cdf below");
        close(at(&cdf, &[("x", 4.0)]), 1.0, 1e-12, "bernoulli cdf above");
        check_moments("bernoulli(0.3)", b.mean(), b.variance(), b.entropy());
        // Binomial
        for (n, p) in [(1_u64, 0.5), (12, 0.35), (30, 0.9), (60, 0.05)] {
            let d = Binomial::new(p, n).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("binomial_dist({n}, {p:?})");
            let (pmf, cdf) = (s(&format!("pdf({name}, x)")), s(&format!("cdf({name}, x)")));
            for k in 0..=n {
                close(
                    at(&pmf, &[("x", k as f64)]),
                    d.pmf(k),
                    1e-9,
                    &format!("{name} pmf at {k}"),
                );
                if k < n {
                    close(
                        at(&cdf, &[("x", k as f64)]),
                        DiscreteCDF::cdf(&d, k),
                        1e-9,
                        &format!("{name} cdf at {k}"),
                    );
                }
            }
            check_moments(&name, d.mean(), d.variance(), None);
            // At literal points the distribution function is an exact sum.
            for k in [0_u64, n / 2, n - 1] {
                close(
                    value(&format!("cdf({name}, {k})")),
                    DiscreteCDF::cdf(&d, k),
                    1e-9,
                    &format!("{name} exact cdf at {k}"),
                );
            }
        }
        // Poisson
        for rate in [0.3, 3.2, 17.0] {
            let d = Poisson::new(rate).unwrap_or_else(|e| panic!("{e}"));
            let name = format!("poisson({rate:?})");
            let (pmf, cdf) = (s(&format!("pdf({name}, x)")), s(&format!("cdf({name}, x)")));
            for k in 0..=25_u64 {
                close(
                    at(&pmf, &[("x", k as f64)]),
                    d.pmf(k),
                    1e-9,
                    &format!("{name} pmf at {k}"),
                );
                close(
                    at(&cdf, &[("x", k as f64)]),
                    DiscreteCDF::cdf(&d, k),
                    1e-9,
                    &format!("{name} cdf at {k}"),
                );
            }
            check_moments(&name, d.mean(), d.variance(), None);
            for k in [0_u64, 4, 25] {
                close(
                    value(&format!("cdf({name}, {k})")),
                    DiscreteCDF::cdf(&d, k),
                    1e-9,
                    &format!("{name} exact cdf at {k}"),
                );
            }
        }
        // Exact results for exact parameters.
        assert_eq!(s("pdf(binomial_dist(10, 1/2), 3)"), "15/128");
        assert_eq!(s("cdf(binomial_dist(10, 1/2), 3)"), "11/64");
        assert_eq!(s("cdf(binomial_dist(10, 1/2), 10)"), "1");
        assert_eq!(s("cdf(binomial_dist(10, 1/2), -1)"), "0");
        assert_eq!(s("cdf(bernoulli(1/3), 0)"), "2/3");
        assert_eq!(s("expectation(binomial_dist(n, p))"), "n*p");
        assert_eq!(s("variance_of(binomial_dist(n, p))"), "n*p*(1 - p)");
        assert_eq!(s("cdf(poisson(2), 3)"), "19/3*exp(-2)");
        assert_eq!(s("cdf(poisson(2), -1)"), "0");
        assert_eq!(
            s("cdf(binomial_dist(n, p), 3)"),
            "beta_reg(n - 3, 4, 1 - p)"
        );
    }

    #[test]
    fn moment_generating_functions_match_expectations() {
        let normal = Normal::new(0.5, 1.5).unwrap_or_else(|e| panic!("{e}"));
        let exponential = Exp::new(2.0).unwrap_or_else(|e| panic!("{e}"));
        let gamma = Gamma::new(2.5, 1.5).unwrap_or_else(|e| panic!("{e}"));
        let uniform = Uniform::new(-1.0, 2.0).unwrap_or_else(|e| panic!("{e}"));
        let chi = ChiSquared::new(5.0).unwrap_or_else(|e| panic!("{e}"));
        for t in [-0.6, -0.1, 0.2, 0.4] {
            #[allow(clippy::type_complexity)] // one-off table of (name, integrand, lo, hi)
            let cases: [(&str, Box<dyn Fn(f64) -> f64>, f64, f64); 5] = [
                (
                    "normal(0.5, 1.5)",
                    Box::new(|x| normal.pdf(x) * (t * x).exp()),
                    -30.0,
                    30.0,
                ),
                (
                    "exponential(2.0)",
                    Box::new(|x| exponential.pdf(x) * (t * x).exp()),
                    0.0,
                    40.0,
                ),
                (
                    "gamma_dist(2.5, 1.5)",
                    Box::new(|x| gamma.pdf(x) * (t * x).exp()),
                    1e-12,
                    60.0,
                ),
                (
                    "uniform(-1.0, 2.0)",
                    Box::new(|x| uniform.pdf(x) * (t * x).exp()),
                    -1.0,
                    2.0,
                ),
                (
                    "chi_squared(5.0)",
                    Box::new(|x| chi.pdf(x) * (t * x).exp()),
                    1e-12,
                    500.0,
                ),
            ];
            for (name, integrand, lo, hi) in cases {
                let m = s(&format!("mgf({name}, s)"));
                close(
                    at(&m, &[("s", t)]),
                    integrate(&*integrand, lo, hi),
                    1e-6,
                    &format!("mgf of {name} at {t}"),
                );
            }
        }
        for t in [-0.7_f64, 0.3, 1.1] {
            let bernoulli = 0.7 + 0.3 * t.exp();
            close(
                at(&s("mgf(bernoulli(0.3), s)"), &[("s", t)]),
                bernoulli,
                1e-12,
                "bernoulli mgf",
            );
            let binomial = Binomial::new(0.35, 12).unwrap_or_else(|e| panic!("{e}"));
            let by_sum: f64 = (0..=12_u64)
                .map(|k| binomial.pmf(k) * (t * k as f64).exp())
                .sum();
            close(
                at(&s("mgf(binomial_dist(12, 0.35), s)"), &[("s", t)]),
                by_sum,
                1e-10,
                "binomial mgf",
            );
            let poisson = Poisson::new(2.5).unwrap_or_else(|e| panic!("{e}"));
            let by_sum: f64 = (0..=120_u64)
                .map(|k| poisson.pmf(k) * (t * k as f64).exp())
                .sum();
            close(
                at(&s("mgf(poisson(2.5), s)"), &[("s", t)]),
                by_sum,
                1e-9,
                "poisson mgf",
            );
        }
        // The derivative of the mgf at zero is the mean.
        let mean = numeric(
            &sets(),
            "diff(mgf(normal(0.5, 1.5), t), t)",
            &[("t", 0.0)],
            1e-9,
        )
        .0;
        close(mean, 0.5, 1e-6, "mgf derivative");
        // No closed form: stays a request.
        stays("mgf(beta_dist(2, 3), t)");
        stays("mgf(student_t(5), t)");
    }

    #[test]
    fn densities_integrate_to_one() {
        let sets = sets();
        for (density, lo, hi) in [
            ("normal(0, 1)", "-oo", "oo"),
            ("normal(-2, 3)", "-oo", "oo"),
            ("student_t(5)", "-oo", "oo"),
            ("student_t(2.5)", "-oo", "oo"),
            ("uniform(1, 4)", "1", "4"),
            ("exponential(2)", "0", "oo"),
            ("gamma_dist(3, 2)", "0", "oo"),
            ("gamma_dist(2.5, 0.5)", "0", "oo"),
            ("beta_dist(2, 3)", "0", "1"),
            ("beta_dist(2.5, 4.5)", "0", "1"),
            ("chi_squared(4)", "0", "oo"),
            ("chi_squared(7)", "0", "oo"),
        ] {
            let (total, _) = numeric(
                &sets,
                &format!("defint(pdf({density}, x), x, {lo}, {hi})"),
                &[],
                1e-9,
            );
            assert!(
                (total - 1.0).abs() < 1e-9,
                "{density} integrates to {total}"
            );
        }
    }

    #[test]
    fn distribution_functions_differentiate_to_densities() {
        let sets = sets();
        for (density, x) in [
            ("normal(1, 2)", 0.7),
            ("gamma_dist(2.5, 1.3)", 1.4),
            ("beta_dist(2, 3.5)", 0.35),
            ("chi_squared(4)", 2.2),
            ("exponential(1.5)", 0.8),
        ] {
            let slope = numeric(
                &sets,
                &format!("diff(cdf({density}, x), x)"),
                &[("x", x)],
                1e-9,
            )
            .0;
            let pdf = numeric(&sets, &format!("pdf({density}, x)"), &[("x", x)], 1e-9).0;
            close(slope, pdf, 1e-7, &format!("d/dx cdf of {density}"));
        }
    }

    #[test]
    fn requests_agree_between_the_phases() {
        let sets = sets();
        for src in [
            "pdf(normal(1, 2), x)",
            "cdf(gamma_dist(2, 3), x)",
            "cdf(student_t(4), x)",
            "pdf(beta_dist(2, 3), x)",
            "entropy_of(normal(0, s))",
            "variance_of(uniform(a, 5))",
        ] {
            let bindings = [("x", 0.6), ("s", 1.7), ("a", 1.0)];
            let symbolic = at(&s(src), &bindings);
            let (numerical, _) = numeric(&sets, src, &bindings, 1e-10);
            close(symbolic, numerical, 1e-9, src);
        }
    }

    #[test]
    fn distribution_requests_stay_requests_when_meaningless() {
        for src in [
            "pdf(d, x)",
            "cdf(normal, x)",
            "pdf(normal(0, -1), 0)",
            "pdf(normal(0, 0), 0)",
            "cdf(uniform(3, 1), 2)",
            "expectation(exponential(0))",
            "cdf(bernoulli(3/2), 0)",
            "pdf(binomial_dist(2.5, 1/2), 1)",
            "pdf(binomial_dist(-1, 1/2), 1)",
            "variance_of(student_t(2))",
            "expectation(student_t(1))",
            "entropy_of(binomial_dist(3, 1/2))",
            "entropy_of(poisson(2))",
            "expectation(beta_dist(0, 1))",
            "expectation(3)",
        ] {
            stays(src);
        }
        // The constructors alone are inert.
        assert_eq!(s("normal(mu, sigma)"), "normal(mu, sigma)");
        assert_eq!(s("normal(0, 1)"), "normal(0, 1)");
    }

    #[test]
    fn incomplete_gamma_and_beta_functions() {
        let r = |src: &str, facts: &[(&str, Facts)]| reduce_with(&sets(), src, facts).0;
        assert_eq!(r("gamma_lr(a, 0)", &[("a", Facts::POSITIVE)]), "0");
        // The right-hand side is larger, so it shows only through cancellation.
        assert_eq!(
            r("gamma_lr(1, x) + exp(-x)", &[("x", Facts::NONNEGATIVE)]),
            "1"
        );
        assert_ne!(r("gamma_lr(1, x) + exp(-x)", &[]), "1");
        assert_eq!(
            r(
                "beta_reg(a, b, 0)",
                &[("a", Facts::POSITIVE), ("b", Facts::POSITIVE)]
            ),
            "0"
        );
        assert_eq!(
            r(
                "beta_reg(a, b, 1)",
                &[("a", Facts::POSITIVE), ("b", Facts::POSITIVE)]
            ),
            "1"
        );
        // Values against the standard closed forms.
        close(
            at("gamma_lr(1, x)", &[("x", 0.7)]),
            1.0 - (-0.7_f64).exp(),
            1e-14,
            "P(1, x)",
        );
        close(
            at("gamma_lr(a, x)", &[("a", 0.5), ("x", 0.7)]),
            statrs::function::erf::erf(0.7_f64.sqrt()),
            1e-10,
            "P(1/2, x)",
        );
        close(
            at("beta_reg(2, 1, x)", &[("x", 0.3)]),
            0.09,
            1e-14,
            "I_x(2, 1)",
        );
        close(
            at("beta_reg(1, 3, x)", &[("x", 0.3)]),
            1.0 - 0.7_f64.powi(3),
            1e-14,
            "I_x(1, 3)",
        );
        assert_eq!(at("beta_reg(2, 3, x)", &[("x", -1.0)]), 0.0);
        assert_eq!(at("beta_reg(2, 3, x)", &[("x", 2.0)]), 1.0);
        assert_eq!(at("gamma_lr(2, x)", &[("x", -1.0)]), 0.0);
        assert!(at("gamma_lr(-1, x)", &[("x", 1.0)]).is_nan());
        assert!(at("beta_reg(0, 1, x)", &[("x", 0.5)]).is_nan());
        // Their derivatives are the densities.
        let slope = numeric(&sets(), "diff(gamma_lr(2, x), x)", &[("x", 1.3)], 1e-9).0;
        close(slope, 1.3 * (-1.3_f64).exp(), 1e-7, "d/dx P(2, x)");
    }

    /// The two numbers of a `list(statistic, p_value)` result.
    fn pair(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> (f64, f64) {
        let text = s(src);
        let v = items(&text);
        assert_eq!(v.len(), 2, "{src} => {text}");
        (at(&v[0], bindings), at(&v[1], bindings))
    }

    fn mean_of(v: &[f64]) -> f64 {
        v.iter().sum::<f64>() / count(v.len())
    }

    fn sample_var_of(v: &[f64]) -> f64 {
        let m = mean_of(v);
        v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (count(v.len()) - 1.0)
    }

    fn student_p(
        t: f64,
        df: f64,
    ) -> f64 {
        2.0 * (1.0 - StudentsT::new(0.0, 1.0, df).map_or(f64::NAN, |d| d.cdf(t.abs())))
    }

    const A: [f64; 6] = [4.1, 5.3, 6.2, 5.9, 4.8, 5.5];
    const B: [f64; 5] = [6.1, 7.2, 5.9, 6.8, 7.5];

    fn list_of(v: &[f64]) -> String {
        format!(
            "list({})",
            v.iter()
                .map(|x| format!("{x:?}"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    #[test]
    fn z_and_t_tests_against_textbook_formulas() {
        let normal = Normal::new(0.0, 1.0).unwrap_or_else(|e| panic!("{e}"));
        // Four observations, known sigma = 1: z = 1, p = 0.3173...
        let (z, p) = pair("z_test(list(1, 2, 3, 4), 2, 1)", &[]);
        close(z, 1.0, 1e-14, "z");
        close(p, 2.0 * (1.0 - normal.cdf(1.0)), 1e-12, "p");
        close(p, 0.317_310_507_862_914_1, 1e-9, "p known");
        let (z, p) = pair(&format!("z_test({}, 5, 0.8)", list_of(&A)), &[]);
        let want = (mean_of(&A) - 5.0) / (0.8 / 6.0_f64.sqrt());
        close(z, want, 1e-13, "z");
        close(p, 2.0 * (1.0 - normal.cdf(want.abs())), 1e-12, "p");
        // t-test
        let (t, p) = pair(&format!("t_test({}, 5)", list_of(&A)), &[]);
        let want = (mean_of(&A) - 5.0) / (sample_var_of(&A) / 6.0).sqrt();
        close(t, want, 1e-13, "t");
        close(p, student_p(want, 5.0), 1e-12, "p");
        // A textbook value: t = 0.7746, p = 0.4950 for (1, 2, 3, 4) against 2.
        let (t, p) = pair("t_test(list(1, 2, 3, 4), 2)", &[]);
        close(t, 0.5 / (5.0_f64 / 12.0).sqrt(), 1e-13, "t");
        close(p, 0.495_025_346_059_711_8, 1e-10, "p known");
        // Exactly at the null value nothing is significant.
        let (t, p) = pair("t_test(list(1, 2, 3), 2)", &[]);
        assert_eq!((t, p), (0.0, 1.0));
    }

    #[test]
    fn two_sample_tests_against_textbook_formulas() {
        let (a, b) = (list_of(&A), list_of(&B));
        let (m1, m2) = (mean_of(&A), mean_of(&B));
        let (v1, v2) = (sample_var_of(&A), sample_var_of(&B));
        let (n1, n2) = (6.0, 5.0);
        // Welch, with and without a hypothesised difference.
        for mu in [0.0, -1.0] {
            let (t, p) = pair(&format!("welch_test({a}, {b}, {mu})"), &[]);
            let se = (v1 / n1 + v2 / n2).sqrt();
            let df = (v1 / n1 + v2 / n2).powi(2)
                / ((v1 / n1).powi(2) / (n1 - 1.0) + (v2 / n2).powi(2) / (n2 - 1.0));
            close(t, (m1 - m2 - mu) / se, 1e-12, "welch t");
            close(p, student_p((m1 - m2 - mu) / se, df), 1e-10, "welch p");
        }
        // The kernel gives the same numbers for a zero difference.
        let (kt, kp) = num::welch_t_test(&A, &B);
        let (t, p) = pair(&format!("welch_test({a}, {b}, 0)"), &[]);
        close(t, kt, 1e-13, "kernel t");
        close(p, kp, 1e-12, "kernel p");
        // Pooled variance.
        let (t, p) = pair(&format!("pooled_t_test({a}, {b})"), &[]);
        let pooled = (n2 - 1.0).mul_add(v2, (n1 - 1.0) * v1) / (n1 + n2 - 2.0);
        let want = (m1 - m2) / (pooled * (1.0 / n1 + 1.0 / n2)).sqrt();
        close(t, want, 1e-12, "pooled t");
        close(p, student_p(want, n1 + n2 - 2.0), 1e-10, "pooled p");
        // Chi-squared goodness of fit: (100 + 0 + 100) / 20 = 10 on 2 degrees.
        let (x, p) = pair("chi_squared_test(list(10, 20, 30), list(20, 20, 20))", &[]);
        close(x, 10.0, 1e-14, "chi2");
        close(p, (-5.0_f64).exp(), 1e-12, "chi2 p");
        let (x, p) = pair(
            "chi_squared_test(list(16, 18, 16, 14, 12, 12), list(16, 16, 16, 16, 16, 8))",
            &[],
        );
        close(x, 3.5, 1e-14, "die chi2");
        let chi = ChiSquared::new(5.0).unwrap_or_else(|e| panic!("{e}"));
        close(p, 1.0 - ContinuousCDF::cdf(&chi, 3.5), 1e-12, "die p");
        // One-way ANOVA: SSB = 26 on 2, SSW = 6 on 6, F = 13.
        let (f, p) = pair(
            "anova(list(list(1, 2, 3), list(2, 3, 4), list(5, 6, 7)))",
            &[],
        );
        close(f, 13.0, 1e-13, "F");
        let dist =
            statrs::distribution::FisherSnedecor::new(2.0, 6.0).unwrap_or_else(|e| panic!("{e}"));
        close(p, 1.0 - ContinuousCDF::cdf(&dist, 13.0), 1e-12, "F p");
    }

    #[test]
    fn confidence_intervals() {
        // t(0.975; 3) = 3.182446305...
        let v = items(&s("confidence_interval(list(1, 2, 3, 4), 0.95)"));
        let half = 3.182_446_305_284_263 * (5.0_f64 / 12.0).sqrt();
        close(at(&v[0], &[]), 2.5 - half, 1e-10, "lower");
        close(at(&v[1], &[]), 2.5 + half, 1e-10, "upper");
        // Known sigma: z(0.975) = 1.959963984540054.
        let v = items(&s("z_interval(list(1, 2, 3, 4), 1, 0.95)"));
        close(
            at(&v[0], &[]),
            2.5 - 1.959_963_984_540_054 / 2.0,
            1e-10,
            "z lower",
        );
        close(
            at(&v[1], &[]),
            2.5 + 1.959_963_984_540_054 / 2.0,
            1e-10,
            "z upper",
        );
        // A higher level gives a wider interval; the interval is centred.
        let wide = items(&s(&format!("confidence_interval({}, 0.99)", list_of(&A))));
        let narrow = items(&s(&format!("confidence_interval({}, 0.9)", list_of(&A))));
        let (wl, wu) = (at(&wide[0], &[]), at(&wide[1], &[]));
        let (nl, nu) = (at(&narrow[0], &[]), at(&narrow[1], &[]));
        assert!(wl < nl && nu < wu);
        close(f64::midpoint(wl, wu), mean_of(&A), 1e-12, "centre");
    }

    #[test]
    fn symbolic_tests_are_formulas_that_agree_with_the_numeric_ones() {
        let data = [("a", 1.5), ("b", 2.0), ("c", 4.5), ("d", 3.0)];
        let num_a = "list(1.5, 2.0, 4.5)";
        let num_b = "list(1.5, 2.0, 4.5, 3.0)";
        let cases = [
            (
                "z_test(list(a, b, c), 2, 3/2)",
                format!("z_test({num_a}, 2, 3/2)"),
            ),
            ("t_test(list(a, b, c), 1)", format!("t_test({num_a}, 1)")),
            ("t_test(list(a, b, c, d), m)", format!("t_test({num_b}, 2)")),
            (
                "welch_test(list(a, b, c), list(c, d, a, d), 1/2)",
                format!("welch_test({num_a}, list(4.5, 3.0, 1.5, 3.0), 1/2)"),
            ),
            (
                "pooled_t_test(list(a, b, c), list(c, d, a, d))",
                format!("pooled_t_test({num_a}, list(4.5, 3.0, 1.5, 3.0))"),
            ),
            (
                "chi_squared_test(list(a, b, c), list(d, d, d))",
                format!("chi_squared_test({num_a}, list(3.0, 3.0, 3.0))"),
            ),
        ];
        let bindings: Vec<(&str, f64)> = data.iter().copied().chain([("m", 2.0)]).collect();
        for (symbolic, numerical) in cases {
            let (a, b) = pair(symbolic, &bindings);
            let (c, d) = pair(&numerical, &[]);
            close(a, c, 1e-10, &format!("{symbolic}: statistic"));
            close(b, d, 1e-9, &format!("{symbolic}: p-value"));
        }
        // The p-value of a symbolic test is written with the cdf request,
        // which reduces to the closed form.
        let text = s("t_test(list(a, b, c), 0)");
        assert!(text.contains("beta_reg"), "{text}");
        assert!(!text.contains("cdf"), "{text}");
        let text = s("z_test(list(a, b), m, sigma)");
        assert!(text.contains("erf"), "{text}");
    }

    #[test]
    fn tests_refuse_degenerate_input() {
        for src in [
            "t_test(list(1), 0)",
            "t_test(list(3, 3, 3), 0)",
            "z_test(list(), 0, 1)",
            "z_test(list(1, 2), 0, 0)",
            "z_test(list(1, 2), 0, -1)",
            "welch_test(list(1), list(1, 2), 0)",
            "welch_test(list(1, 1), list(2, 2), 0)",
            "pooled_t_test(list(1, 2), list(3))",
            "chi_squared_test(list(1, 2), list(0, 1))",
            "chi_squared_test(list(1, 2), list(1, 2, 3))",
            "chi_squared_test(list(1), list(1))",
            "anova(list(list(1, 2)))",
            "anova(list(list(1, 2), list()))",
            "confidence_interval(list(1), 0.95)",
            "confidence_interval(list(1, 2, 3), 1.5)",
            "confidence_interval(list(a, b, c), 0.95)",
            "z_interval(list(1, 2, 3), 0, 0.95)",
            "t_test(x, 0)",
        ] {
            stays(src);
        }
    }

    #[test]
    fn linear_regression_is_exact_for_exact_data() {
        assert_eq!(
            s("linear_regression(list(1, 2, 3), list(2, 4, 6))"),
            "list(0, 2)"
        );
        assert_eq!(
            s("linear_regression(list(1, 2, 3, 4), list(1, 3, 2, 5))"),
            "list(0, 11/10)"
        );
        assert_eq!(
            s("linear_regression(list(0, 1, 2), list(1, 0, -1))"),
            "list(1, -1)"
        );
        assert_eq!(
            s("linear_regression(list(0, 10), list(3, 3))"),
            "list(3, 0)"
        );
        // Normal equations for messy rational data, checked exactly.
        let xs = [1_i64, 2, 4, 7, 9];
        let ys = [3_i64, 3, 8, 11, 20];
        let text = s(&format!(
            "linear_regression({}, {})",
            list_of(&xs.map(|v| v as f64)).replace(".0", ""),
            list_of(&ys.map(|v| v as f64)).replace(".0", "")
        ));
        let v = items(&text);
        let (b0, b1) = (at(&v[0], &[]), at(&v[1], &[]));
        let residuals: Vec<f64> = xs
            .iter()
            .zip(&ys)
            .map(|(&x, &y)| y as f64 - b0 - b1 * x as f64)
            .collect();
        close(
            residuals.iter().sum::<f64>(),
            0.0,
            1e-12,
            "residuals sum to zero",
        );
        close(
            residuals
                .iter()
                .zip(&xs)
                .map(|(r, &x)| r * x as f64)
                .sum::<f64>(),
            0.0,
            1e-12,
            "orthogonal to x",
        );
        // The slope is cov / var, exactly.
        assert_eq!(
            v[1],
            s(
                "covariance(list(1, 2, 4, 7, 9), list(3, 3, 8, 11, 20)) / variance(list(1, 2, 4, 7, 9))"
            )
        );
    }

    #[test]
    fn linear_regression_on_floats_and_symbols() {
        let (x, y) = ([1.0, 2.0, 3.5, 4.0], [2.0, 4.1, 6.0, 8.3]);
        let v = items(&s(&format!(
            "linear_regression({}, {})",
            list_of(&x),
            list_of(&y)
        )));
        let pairs: Vec<(f64, f64)> = x.iter().copied().zip(y).collect();
        let (slope, intercept) = num::simple_linear_regression(&pairs);
        close(at(&v[0], &[]), intercept, 1e-12, "intercept");
        close(at(&v[1], &[]), slope, 1e-12, "slope");
        // Symbolic data: the formula at the same values.
        let text = s("linear_regression(list(a, b, c, d), list(e, f, g, h))");
        let bindings = [
            ("a", 1.0),
            ("b", 2.0),
            ("c", 3.5),
            ("d", 4.0),
            ("e", 2.0),
            ("f", 4.1),
            ("g", 6.0),
            ("h", 8.3),
        ];
        let v = items(&text);
        close(at(&v[0], &bindings), intercept, 1e-12, "symbolic intercept");
        close(at(&v[1], &bindings), slope, 1e-12, "symbolic slope");
        // A vertical cloud has no regression line.
        stays("linear_regression(list(2, 2, 2), list(1, 2, 3))");
        stays("linear_regression(list(1), list(1))");
        stays("linear_regression(list(1, 2), list(1, 2, 3))");
    }

    #[test]
    fn multiple_and_polynomial_regression() {
        // A plane recovered exactly: y = 2 + 3 x1 - x2.
        assert_eq!(
            s(
                "linear_regression(list(list(1, 2, 3, 4, 5), list(2, 1, 5, 3, 9)), list(3, 7, 6, 11, 8))"
            ),
            "list(2, 3, -1)"
        );
        // Noisy data: the normal equations hold exactly.
        let x1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x2 = [1.0, 0.0, 1.0, 3.0, 2.0, 5.0];
        let y = [3.0, 3.0, 5.0, 9.0, 7.0, 14.0];
        let text = s(&format!(
            "linear_regression(list({}, {}), {})",
            list_of(&x1),
            list_of(&x2),
            list_of(&y)
        ));
        let v = items(&text);
        let beta: Vec<f64> = v.iter().map(|t| at(t, &[])).collect();
        assert_eq!(beta.len(), 3);
        let residual: Vec<f64> = (0..6)
            .map(|i| y[i] - beta[0] - beta[1] * x1[i] - beta[2] * x2[i])
            .collect();
        close(residual.iter().sum::<f64>(), 0.0, 1e-11, "sum");
        close(
            residual.iter().zip(&x1).map(|(r, x)| r * x).sum::<f64>(),
            0.0,
            1e-11,
            "x1",
        );
        close(
            residual.iter().zip(&x2).map(|(r, x)| r * x).sum::<f64>(),
            0.0,
            1e-11,
            "x2",
        );
        // Collinear predictors have no unique fit.
        stays("linear_regression(list(list(1, 2, 3, 4), list(2, 4, 6, 8)), list(1, 2, 4, 8))");
        // Polynomials.
        assert_eq!(
            s("polynomial_regression(list(0, 1, 2, 3), list(1, 2, 5, 10), 2)"),
            "list(1, 0, 1)"
        );
        assert_eq!(
            s("polynomial_regression(list(-1, 0, 1, 2, 3), list(-2, 1, 4, 13, 34), 3)"),
            "list(1, 2, 0, 1)"
        );
        assert_eq!(
            s("polynomial_regression(list(1, 2, 3), list(1, 2, 3), 1)"),
            "list(0, 1)"
        );
        assert_eq!(
            s("polynomial_regression(list(1, 2, 3), list(5, 5, 5), 0)"),
            "list(5)"
        );
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [1.1, 1.9, 5.2, 9.8, 17.1, 25.3];
        let v = items(&s(&format!(
            "polynomial_regression({}, {}, 2)",
            list_of(&xs),
            list_of(&ys)
        )));
        let c: Vec<f64> = v.iter().map(|t| at(t, &[])).collect();
        let r: Vec<f64> = xs
            .iter()
            .zip(&ys)
            .map(|(x, y)| y - c[0] - c[1] * x - c[2] * x * x)
            .collect();
        close(r.iter().sum::<f64>(), 0.0, 1e-11, "poly sum");
        close(
            r.iter().zip(&xs).map(|(r, x)| r * x).sum::<f64>(),
            0.0,
            1e-10,
            "poly x",
        );
        close(
            r.iter().zip(&xs).map(|(r, x)| r * x * x).sum::<f64>(),
            0.0,
            1e-9,
            "poly x^2",
        );
        stays("polynomial_regression(list(1, 2), list(1, 2), 2)");
        stays("polynomial_regression(list(1, 2, 3), list(1, 2, 3), -1)");
        stays("polynomial_regression(list(a, 2, 3), list(1, 2, 3), 1)");
        stays("polynomial_regression(list(1, 1, 1), list(1, 2, 3), 1)");
    }

    #[test]
    fn nonlinear_regression_solves_the_normal_equations() {
        let run = |src: &str| simplify(&[stats(), crate::rules::solve()], src);
        assert_eq!(
            run("nonlinear_regression(list(0, 1, 2), list(1, 3, 5), a*x + b, x, list(a, b))"),
            "list(list(2, 1))"
        );
        assert_eq!(
            run("nonlinear_regression(list(1, 2), list(2, 8), a*x^2, x, list(a))"),
            "list(list(2))"
        );
        // Without the solver the request stays.
        stays("nonlinear_regression(list(0, 1, 2), list(1, 3, 5), a*x + b, x, list(a, b))");
        let (text, reduced) = reduce_with(
            &[stats(), crate::rules::solve()],
            "nonlinear_regression(list(0, 1), list(1, 2), a*exp(b*x), x, list(a, b))",
            &[],
        );
        assert!(!reduced, "{text}");
    }

    #[test]
    #[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
    fn entropy_of_probability_lists() {
        assert_eq!(s("entropy(list(1/2, 1/2))"), "ln(2)");
        assert_eq!(s("entropy(list(1/4, 1/4, 1/4, 1/4))"), "2*ln(2)");
        assert_eq!(s("entropy(list(1/2, 1/4, 1/4))"), "3/2*ln(2)");
        assert_eq!(s("entropy(list(1/3, 1/3, 1/3))"), "ln(3)");
        assert_eq!(s("entropy(list(1))"), "0");
        assert_eq!(s("entropy(list(1, 0))"), "0");
        assert_eq!(s("entropy(list(1/2, 0, 1/2))"), "ln(2)");
        assert_eq!(s("entropy_bits(list(1/2, 1/4, 1/4))"), "3/2");
        assert_eq!(
            s("entropy_bits(list(1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8))"),
            "3"
        );
        // Uniform on n outcomes has entropy ln n.
        for n in 2..9_i64 {
            let uniform = vec![format!("1/{n}"); n as usize].join(", ");
            close(
                value(&format!("entropy(list({uniform}))")),
                (n as f64).ln(),
                1e-13,
                "uniform",
            );
        }
        // Symbolic probabilities.
        let text = s("entropy(list(p, 1 - p))");
        for p in [0.1, 0.3, 0.5, 0.9] {
            close(
                at(&text, &[("p", p)]),
                -p * f64::ln(p) - (1.0 - p) * f64::ln(1.0 - p),
                1e-13,
                "binary entropy",
            );
        }
        // Floats: the kernel works in bits.
        let probabilities = [0.5, 0.25, 0.125, 0.125];
        let text = list_of(&probabilities);
        close(
            value(&format!("entropy({text})")),
            num::shannon_entropy(&probabilities) * std::f64::consts::LN_2,
            1e-14,
            "nats",
        );
        close(
            value(&format!("entropy_bits({text})")),
            num::shannon_entropy(&probabilities),
            1e-14,
            "bits",
        );
        close(
            value(&format!("entropy_bits({text})")),
            1.75,
            1e-14,
            "known bits",
        );
        for src in [
            "entropy(list())",
            "entropy(list(-1/2, 3/2))",
            "entropy(x)",
            "gini(list())",
        ] {
            stays(src);
        }
    }

    #[test]
    fn relative_and_cross_entropy() {
        assert_eq!(s("kl_divergence(list(1/2, 1/2), list(1/2, 1/2))"), "0");
        assert_eq!(
            s("kl_divergence(list(1/2, 1/2), list(1/4, 3/4))"),
            "ln(2) - 1/2*ln(3)"
        );
        assert_eq!(s("cross_entropy(list(1/2, 1/2), list(1/2, 1/2))"), "ln(2)");
        assert_eq!(
            s("cross_entropy(list(1/2, 1/2), list(1/4, 3/4))"),
            "2*ln(2) - 1/2*ln(3)"
        );
        // A zero in p contributes nothing, a zero in q against p > 0 is infinite.
        assert_eq!(s("kl_divergence(list(1, 0), list(1/2, 1/2))"), "ln(2)");
        stays("kl_divergence(list(1/2, 1/2), list(1, 0))");
        stays("kl_divergence(list(1/2, 1/2), list(1/2))");
        stays("cross_entropy(list(1/2, 1/2), list(1, 0))");
        // Gibbs: cross entropy = entropy + relative entropy; KL >= 0.
        for (p, q) in [
            ("1/2, 1/2", "1/4, 3/4"),
            ("1/5, 2/5, 2/5", "1/2, 1/4, 1/4"),
            ("0.2, 0.3, 0.5", "0.3, 0.3, 0.4"),
        ] {
            let h = value(&format!("entropy(list({p}))"));
            let kl = value(&format!("kl_divergence(list({p}), list({q}))"));
            let ce = value(&format!("cross_entropy(list({p}), list({q}))"));
            close(ce, h + kl, 1e-13, "cross = H + KL");
            assert!(kl > 0.0);
        }
        // Symbolic.
        let text = s("kl_divergence(list(p, 1 - p), list(q, 1 - q))");
        let (p, q) = (0.3, 0.6);
        close(
            at(&text, &[("p", p), ("q", q)]),
            p * f64::ln(p / q) + (1.0 - p) * f64::ln((1.0 - p) / (1.0 - q)),
            1e-13,
            "binary KL",
        );
    }

    #[test]
    fn joint_entropies_and_mutual_information() {
        assert_eq!(
            s("joint_entropy(list(list(1/4, 1/4), list(1/4, 1/4)))"),
            "2*ln(2)"
        );
        assert_eq!(
            s("mutual_information(list(list(1/4, 1/4), list(1/4, 1/4)))"),
            "0"
        );
        assert_eq!(
            s("mutual_information(list(list(1/2, 0), list(0, 1/2)))"),
            "ln(2)"
        );
        assert_eq!(
            s("conditional_entropy(list(list(1/2, 0), list(0, 1/2)))"),
            "0"
        );
        assert_eq!(
            s("conditional_entropy(list(list(1/4, 1/4), list(1/4, 1/4)))"),
            "ln(2)"
        );
        // An asymmetric joint distribution, against the definitions.
        let joint = [[1.0 / 8.0, 3.0 / 8.0], [1.0 / 4.0, 1.0 / 4.0]];
        let plogp = |p: f64| p * p.ln();
        let hxy = -joint.iter().flatten().map(|&p| plogp(p)).sum::<f64>();
        let hx = -(plogp(0.5) + plogp(0.5));
        let hy = -(plogp(3.0 / 8.0) + plogp(5.0 / 8.0));
        let src = "list(list(1/8, 3/8), list(1/4, 1/4))";
        close(
            value(&format!("joint_entropy({src})")),
            hxy,
            1e-13,
            "H(X,Y)",
        );
        close(
            value(&format!("conditional_entropy({src})")),
            hxy - hx,
            1e-13,
            "H(Y|X)",
        );
        close(
            value(&format!("mutual_information({src})")),
            hx + hy - hxy,
            1e-13,
            "I(X;Y)",
        );
        assert!(value(&format!("mutual_information({src})")) > 0.0);
        stays("joint_entropy(list(list(1/2, 1/2), list(1/2)))");
        stays("joint_entropy(list(1/2, 1/2))");
        stays("mutual_information(list())");
    }

    #[test]
    #[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
    fn gini_impurity() {
        assert_eq!(s("gini(list(1/2, 1/4, 1/4))"), "5/8");
        assert_eq!(s("gini(list(1))"), "0");
        assert_eq!(s("gini(list(1/2, 1/2))"), "1/2");
        for n in 2..8_i64 {
            let uniform = vec![format!("1/{n}"); n as usize].join(", ");
            close(
                value(&format!("gini(list({uniform}))")),
                1.0 - 1.0 / n as f64,
                1e-14,
                "uniform gini",
            );
        }
        let text = s("gini(list(p, 1 - p))");
        close(
            at(&text, &[("p", 0.3)]),
            1.0 - 0.09 - 0.49,
            1e-14,
            "symbolic gini",
        );
    }

    #[test]
    fn statistics_and_information_theory_compose() {
        // Entropy of a normal distribution against the differential-entropy
        // integral -int pdf ln pdf.
        let h = numeric(
            &sets(),
            "defint(-pdf(normal(0, 2), x) * ln(pdf(normal(0, 2), x)), x, -40, 40)",
            &[],
            1e-9,
        )
        .0;
        close(
            h,
            value("entropy_of(normal(0, 2))"),
            1e-7,
            "differential entropy",
        );
        // Requests nest: statistics of a list of expectations.
        assert_eq!(
            s("mean(list(expectation(normal(1, 2)), expectation(normal(3, 4))))"),
            "2"
        );
    }
}
