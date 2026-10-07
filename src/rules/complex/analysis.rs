//! Complex analysis: poles, residues, contour integrals, the argument
//! principle, Cauchy's formulas, radii of convergence, analytic
//! continuation and Möbius transformations.
//!
//! | operator | value |
//! |---|---|
//! | `poles(f, z)` | `list(list(p, order), ...)` |
//! | `zeros_of(f, z)` | `list(list(q, order), ...)` for a polynomial numerator |
//! | `residue(f, z, a)` | the residue of `f` at `a` |
//! | `singularity(f, z, a)` | `regular`, `removable`, `list(pole, m)` or `essential` |
//! | `contour_integral(f, z, C)` | `C = circle(c, r)` (counter-clockwise), `polygon(list(z1, ..., zn))` (closed, in the order given) or `path(g(t), t, t0, t1)` |
//! | `contour_integral_residue_theorem(f, z, list(a1, ...))` | `2 pi I` times the sum of the residues at the listed points |
//! | `is_inside_contour(p, C)`, `winding_number(p, C)` | whether `p` is strictly inside a circle or polygon (`true`/`false`; a point on the contour is not inside); the winding number of `C` around `p` (unreduced for a point on `C`) |
//! | `check_analytic(f, z)` | `true`/`false`: the Cauchy–Riemann equations for `f(x + I y)` |
//! | `branch_cut(f, z)` | `list(ray(p, theta), ...)`: the cuts of the principal branches in `f` |
//! | `cauchy_integral(f, z, a)`, `cauchy_derivative(f, z, a, n)` | `∮ f/(z - a)` and `∮ f/(z - a)^(n+1)` around `a` |
//! | `count_zeros_poles(f, z, C)` | zeros minus poles inside a circle or polygon (the argument principle), with multiplicity |
//! | `radius_of_convergence(f, z, a)` | distance from `a` to the nearest pole (`oo` for entire `f`) |
//! | `continue_along(f, z, list(a0, ..., an), order)` | the Taylor expansion at `an`, continued disc by disc along the points |
//! | `distance(a, b)` | `abs(a - b)` |
//! | `mobius_apply(M, z)`, `mobius_compose(M, N)`, `mobius_inverse(M)` | Möbius maps as matrices `list(list(a, b), list(c, d))` |
//!
//! # Contours
//!
//! A contour is `circle(c, r)`, traversed counter-clockwise, or
//! `polygon(list(z1, ..., zn))`, the closed polygon through the vertices in
//! that order. Residue sums weight each pole by the winding number of the
//! contour around it, so a clockwise polygon gives the negative integral and
//! a polygon traversed twice gives twice the integral. The centre, radius
//! and vertices must be numbers (they may be written with `I`, `sqrt`, ...);
//! a pole **on** the contour makes the integral undefined and the request
//! stays unreduced.
//!
//! # Analyticity and branch cuts
//!
//! `check_analytic(f, z)` substitutes `z = x + I y`, splits `f` into its
//! real and imaginary parts `u + I v` and tests `u_x = v_y`, `u_y = -v_x`:
//! `true` when both differences are identically zero, `false` when they are
//! non-zero at a sample point. Other free symbols are taken as real
//! constants. `branch_cut(f, z)` collects the cuts of `ln`, `sqrt`, powers
//! with non-integer exponents and the inverse trigonometric and hyperbolic
//! functions applied to an argument linear in `z`; a cut is
//! `ray(p, theta)`, the points `p + r exp(I theta)` with `r >= 0`.
//!
//! Functions handled exactly are *meromorphic quotients*: an entire
//! numerator (polynomials, `exp`, `sin`, `cos`, `sinh`, `cosh` of entire
//! arguments) over a polynomial with rational coefficients. The poles are
//! the roots of the denominator's irreducible factors written in radicals
//! (linear, quadratic, binomial `a z^n + b`); residues use the derivative
//! formula, with the polynomial division by `z - p` done symbolically.
//! Every residue whose function evaluates numerically is checked against
//! a numerical contour integral before it is accepted.

use std::collections::HashMap;

use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use super::Ops;
use super::build;
use super::eval_complex;
use crate::graph::Arity;
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
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::complex as numeric;
use crate::rules::calculus::derivative;
use crate::rules::calculus::laurent_expansion;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Poles,
    Zeros,
    Residue,
    Singularity,
    Contour,
    CauchyIntegral,
    CauchyDerivative,
    CountZerosPoles,
    Radius,
    Continue,
    MobiusApply,
    MobiusCompose,
    MobiusInverse,
    CheckAnalytic,
    InsideContour,
    WindingNumber,
    ResidueTheorem,
    BranchCut,
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("circle", Arity::Fixed(2)))?;
    i.op(OpDescriptor::new("path", Arity::Fixed(4)).binder(1, 0b1))?;
    let distance = i.op(OpDescriptor::new("distance", Arity::Fixed(2)).eval(|a| match a {
        | [x, y] => (x - y).abs(),
        | _ => f64::NAN,
    }))?;
    i.graph().ops_mut().set_attr(distance, crate::graph::ComplexEval(|a| Some(num_complex::Complex64::new((*a.first()? - *a.get(1)?).norm(), 0.0))));
    i.rewrites(Tier::Normalize, &["complex/distance: distance(?a, ?b) => abs(?a - ?b)"])?;
    i.op(OpDescriptor::new("polygon", Arity::Fixed(1)))?;
    i.op(OpDescriptor::new("ray", Arity::Fixed(2)))?;
    let table: [(&str, u8, Request); 18] = [
        ("poles", 2, Request::Poles),
        ("zeros_of", 2, Request::Zeros),
        ("residue", 3, Request::Residue),
        ("singularity", 3, Request::Singularity),
        ("contour_integral", 3, Request::Contour),
        ("cauchy_integral", 3, Request::CauchyIntegral),
        ("cauchy_derivative", 4, Request::CauchyDerivative),
        ("count_zeros_poles", 3, Request::CountZerosPoles),
        ("radius_of_convergence", 3, Request::Radius),
        ("continue_along", 4, Request::Continue),
        ("mobius_apply", 2, Request::MobiusApply),
        ("mobius_compose", 2, Request::MobiusCompose),
        ("mobius_inverse", 1, Request::MobiusInverse),
        ("check_analytic", 2, Request::CheckAnalytic),
        ("is_inside_contour", 2, Request::InsideContour),
        ("winding_number", 2, Request::WindingNumber),
        ("contour_integral_residue_theorem", 3, Request::ResidueTheorem),
        ("branch_cut", 2, Request::BranchCut),
    ];
    let ops = Ops::of(i.graph()).ok_or_else(|| RuleError::Invalid { rule: "complex/analysis".into(), reason: "needs elementary" })?;
    for (name, arity, request) in table {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("complex/{name}"), Tier::Reduce, Analysis { op, request, ops });
    }
    Ok(())
}

/// `numer / (lead * prod factor_i^m_i)` with `numer` entire and the
/// factors irreducible over `Q`.
struct Quotient {
    numer: NodeId,
    lead: BigRational,
    factors: Vec<(QPoly, u32)>,
}

/// A pole or zero: an exact point, its numeric value, the factor it is a
/// root of and that factor's power.
#[derive(Clone, Debug)]
struct Root {
    point: NodeId,
    value: Complex64,
    factor: usize,
    multiplicity: u32,
}

/// Whether `node` is entire in `z`: built from `z`, constants, sums,
/// products, non-negative integer powers and entire functions.
fn is_entire(
    graph: &Graph,
    ops: Ops,
    node: NodeId,
) -> bool {
    let entire_fns = ["exp", "sin", "cos", "sinh", "cosh"];
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        let op = graph.op(n);
        let children = graph.children(n);
        let ok = match op {
            | core::ADD | core::MUL => true,
            | core::POW => children
                .get(1)
                .and_then(|&e| graph.number_of(e))
                .is_some_and(|e| e.is_integer() && e.to_f64() >= 0.0),
            | _ if children.is_empty() => true,
            | _ => op != ops.ln && entire_fns.contains(&graph.ops().get(op).name.as_ref()),
        };
        if !ok {
            return false;
        }
        stack.extend_from_slice(children);
    }
    true
}

/// `f` as a meromorphic quotient in `z`.
fn quotient(
    graph: &mut Graph,
    ops: Ops,
    f: NodeId,
    z: NodeId,
) -> Option<Quotient> {
    let term = best(graph, f)?;
    let mut gens = Gens::default();
    let gz = gens.index(graph, z);
    let r = ratio(graph, &mut gens, term, Limits::default())?;
    // The denominator must be a polynomial in z alone, over Q.
    let denom: QPoly = r.denom.univariate_in(gz)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
    if r.denom.support().iter().any(|&g| g != gz) {
        return None;
    }
    let numer = to_term(graph, &gens, &r.numer);
    if !is_entire(graph, ops, numer) {
        return None;
    }
    let (lead, factors) = univariate::factor(&denom);
    let factors = factors.into_iter().map(|(f, m)| (f.into_iter().map(BigRational::from_integer).collect(), m)).collect();
    Some(Quotient { numer, lead, factors })
}

fn rational(
    graph: &mut Graph,
    value: &BigRational,
) -> NodeId {
    graph.num(Number::rat(value.clone()))
}

/// The roots of an irreducible polynomial over `Q` in radicals, with
/// their numeric values. `None` beyond linear, quadratic and binomial.
fn roots_of(
    graph: &mut Graph,
    ops: Ops,
    factor: &[BigRational],
) -> Option<Vec<(NodeId, Complex64)>> {
    let to_f = |q: &BigRational| Number::rat(q.clone()).to_f64();
    match factor {
        | [c0, c1] => {
            let value = -c0 / c1;
            Some(vec![(rational(graph, &value), Complex64::new(to_f(&value), 0.0))])
        },
        | [c, b, a] => {
            // (-b ± sqrt(D)) / 2a; D < 0 for an irreducible quadratic with
            // complex roots, D > 0 (not a square) for real ones.
            let discriminant = b * b - BigRational::from_integer(BigInt::from(4)) * a * c;
            let two_a = BigRational::from_integer(BigInt::from(2)) * a;
            let centre = -b / &two_a;
            let half = graph.num(Number::fraction(1, 2)?);
            let magnitude = rational(graph, &discriminant.abs());
            let radical = build::pow(graph, magnitude, half);
            let scale = rational(graph, &(BigRational::from_integer(BigInt::from(1)) / &two_a));
            let offset = if discriminant.is_negative() {
                let unit = graph.node(ops.unit, &[]);
                build::mul(graph, &[scale, radical, unit])
            } else {
                build::mul(graph, &[scale, radical])
            };
            let centre_node = rational(graph, &centre);
            let root_d = to_f(&discriminant.abs()).sqrt() / to_f(&two_a);
            let mut out = Vec::new();
            for sign in [1_i64, -1] {
                let s = graph.int(sign);
                let signed = build::mul(graph, &[s, offset]);
                let point = build::add(graph, &[centre_node, signed]);
                #[allow(clippy::cast_precision_loss)]
                let value = if discriminant.is_negative() {
                    Complex64::new(to_f(&centre), sign as f64 * root_d)
                } else {
                    Complex64::new(to_f(&centre) + sign as f64 * root_d, 0.0)
                };
                out.push((point, value));
            }
            Some(out)
        },
        | higher => {
            // a z^n + b: |b/a|^(1/n) exp(I (theta + 2 pi k) / n)
            let n = higher.len() - 1;
            let (c0, cn) = (higher.first()?, higher.last()?);
            if higher.get(1..n)?.iter().any(|c| !c.is_zero()) || n > 64 {
                return None;
            }
            let value = -c0 / cn;
            let n_i = i64::try_from(n).ok()?;
            let modulus = rational(graph, &value.abs());
            let exponent = graph.num(Number::fraction(1, n_i)?);
            let radius = build::pow(graph, modulus, exponent);
            let base_angle = i64::from(value.is_negative());
            let mut out = Vec::with_capacity(n);
            #[allow(clippy::cast_precision_loss)]
            let r = to_f(&value.abs()).powf(1.0 / n as f64);
            for k in 0..n_i {
                // angle = pi (base + 2k) / n
                let numerator = base_angle + 2 * k;
                let fraction = graph.num(Number::fraction(numerator, n_i)?);
                let pi = graph.ops().lookup("pi").map(|p| graph.node(p, &[]))?;
                let unit = graph.node(ops.unit, &[]);
                let angle = build::mul(graph, &[fraction, pi, unit]);
                let rotation = build::call(graph, ops.exp, &[angle]);
                let point = build::mul(graph, &[radius, rotation]);
                #[allow(clippy::cast_precision_loss)]
                let theta = std::f64::consts::PI * numerator as f64 / n as f64;
                out.push((point, Complex64::from_polar(r, theta)));
            }
            Some(out)
        },
    }
}

fn poles_of(
    cx: &mut Cx<'_>,
    ops: Ops,
    q: &Quotient,
) -> Option<Vec<Root>> {
    let mut out = Vec::new();
    for (index, (factor, multiplicity)) in q.factors.iter().enumerate() {
        for (point, value) in roots_of(cx.graph, ops, factor)? {
            let point = cx.simplify(point);
            out.push(Root { point, value, factor: index, multiplicity: *multiplicity });
        }
    }
    Some(out)
}

/// `p(z) / (z - a)` by synthetic division with symbolic `a`; `p(a)` must
/// be zero.
fn deflate(
    cx: &mut Cx<'_>,
    p: &[BigRational],
    a: NodeId,
    z: NodeId,
) -> NodeId {
    // b_{n-1} = c_n, b_{k-1} = c_k + a b_k
    let n = p.len() - 1;
    let mut coefficients = vec![cx.graph.int(0); n];
    let mut carry = rational(cx.graph, &p[n]);
    for k in (0..n).rev() {
        coefficients[k] = carry;
        let ck = rational(cx.graph, &p[k]);
        let product = build::mul(cx.graph, &[a, carry]);
        let next = build::add(cx.graph, &[ck, product]);
        carry = cx.simplify(next);
    }
    let mut terms = Vec::with_capacity(n);
    for (k, &c) in coefficients.iter().enumerate() {
        let power = build::powi(cx.graph, z, i64::try_from(k).unwrap_or(0));
        terms.push(build::mul(cx.graph, &[c, power]));
    }
    let sum = build::add(cx.graph, &terms);
    cx.simplify(sum)
}

fn qpoly_term(
    graph: &mut Graph,
    p: &[BigRational],
    z: NodeId,
) -> NodeId {
    let mut terms = Vec::with_capacity(p.len());
    for (k, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let c = rational(graph, c);
        let power = build::powi(graph, z, i64::try_from(k).unwrap_or(0));
        terms.push(build::mul(graph, &[c, power]));
    }
    build::add(graph, &terms)
}

/// `(z - a)^m f(z)` with the factor cancelled: analytic at `a`.
fn regular_part(
    cx: &mut Cx<'_>,
    q: &Quotient,
    root: &Root,
    z: NodeId,
) -> NodeId {
    let mut denominator = vec![rational(cx.graph, &q.lead)];
    for (index, (factor, m)) in q.factors.iter().enumerate() {
        let base = if index == root.factor { deflate(cx, factor, root.point, z) } else { qpoly_term(cx.graph, factor, z) };
        denominator.push(build::powi(cx.graph, base, i64::from(*m)));
    }
    let denominator = build::mul(cx.graph, &denominator);
    let inverse = build::powi(cx.graph, denominator, -1);
    build::mul(cx.graph, &[q.numer, inverse])
}

/// `g^(n)(a) / n!`, simplified.
fn taylor_coefficient(
    cx: &mut Cx<'_>,
    g: NodeId,
    z: NodeId,
    a: NodeId,
    n: u32,
) -> Option<NodeId> {
    let mut current = cx.simplify(g);
    let mut factorial = BigInt::from(1);
    for k in 1..=n {
        current = derivative(cx.graph, current, z)?;
        current = cx.simplify(current);
        factorial *= BigInt::from(k);
    }
    let at = cx.graph.substitute(current, z, a);
    let scale = cx.graph.num(Number::rat(BigRational::new(BigInt::from(1), factorial)));
    let value = build::mul(cx.graph, &[scale, at]);
    Some(cx.simplify(value))
}

/// The numeric function `z -> f(z)` when `f` has no other free symbols.
fn numeric_function(
    graph: &mut Graph,
    f: NodeId,
    z: SymbolId,
) -> Option<impl Fn(Complex64) -> Complex64 + use<>> {
    let term = best(graph, f)?;
    if graph.free_symbols(graph.find(term)).iter().any(|&s| s != z) {
        return None;
    }
    let snapshot = graph.clone();
    Some(move |at: Complex64| {
        let bindings: HashMap<SymbolId, Complex64> = std::iter::once((z, at)).collect();
        eval_complex(&snapshot, term, &bindings).unwrap_or(Complex64::new(f64::NAN, f64::NAN))
    })
}

fn value_of(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Complex64> {
    let term = best(graph, node)?;
    eval_complex(graph, term, &HashMap::new())
}

/// Order of vanishing of the entire numerator at a root (bounded).
fn vanishing_order(
    cx: &mut Cx<'_>,
    numer: NodeId,
    z: NodeId,
    a: NodeId,
    bound: u32,
) -> u32 {
    let mut current = numer;
    for k in 0..bound {
        let at = cx.graph.substitute(current, z, a);
        let at = cx.simplify(at);
        let zero = cx.graph.number_of(at).is_some_and(Number::is_zero);
        if !zero {
            return k;
        }
        match derivative(cx.graph, current, z) {
            | Some(d) => current = cx.simplify(d),
            | None => return k,
        }
    }
    bound
}

struct Analysis {
    op: OpId,
    request: Request,
    ops: Ops,
}

impl Kernel for Analysis {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }
}

/// Where a cut starts in the plane of a function's argument.
#[derive(Copy, Clone)]
enum Start {
    Zero,
    One,
    MinusOne,
    I,
    MinusI,
}

/// A closed contour.
#[derive(Clone, Debug)]
enum Contour {
    /// Counter-clockwise circle.
    Circle { centre: Complex64, radius: f64 },
    /// Closed polygon through the vertices in order.
    Polygon(Vec<Complex64>),
}

/// Reads `circle(c, r)` or `polygon(list(...))` with numeric data.
fn contour_of(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Contour> {
    let term = best(graph, node)?;
    let name = graph.ops().get(graph.op(term)).name.to_string();
    match name.as_str() {
        | "circle" => {
            let &[c, r] = graph.children(term) else {
                return None;
            };
            let centre = value_of(graph, c)?;
            let radius = value_of(graph, r)?;
            (centre.is_finite() && radius.im == 0.0 && radius.re > 0.0).then_some(Contour::Circle { centre, radius: radius.re })
        },
        | "polygon" => {
            let &[list] = graph.children(term) else {
                return None;
            };
            let list = best(graph, list)?;
            if graph.op(list) != core::LIST {
                return None;
            }
            let vertices: Vec<Complex64> =
                graph.children(list).to_vec().into_iter().map(|v| value_of(graph, v).filter(|c| c.is_finite())).collect::<Option<_>>()?;
            (vertices.len() >= 3).then_some(Contour::Polygon(vertices))
        },
        | _ => None,
    }
}

impl Contour {
    /// The winding number around `z`; `None` for a point on the contour.
    fn winding(
        &self,
        z: Complex64,
    ) -> Option<i64> {
        match self {
            | Self::Circle { centre, radius } => {
                let d = (z - centre).norm();
                ((d - radius).abs() > 1e-9 * (1.0 + radius)).then_some(i64::from(d < *radius))
            },
            | Self::Polygon(v) => polygon_winding(v, z),
        }
    }
}

/// Winding number of the closed polygon `v` around `z` (crossing numbers).
fn polygon_winding(
    v: &[Complex64],
    z: Complex64,
) -> Option<i64> {
    let tol = 1e-9 * (1.0 + v.iter().map(|p| p.norm()).fold(0.0, f64::max));
    let mut winding = 0;
    for k in 0..v.len() {
        let (a, b) = (v[k], v[(k + 1) % v.len()]);
        let edge = b - a;
        let length2 = edge.norm_sqr();
        let t = if length2 == 0.0 { 0.0 } else { (((z - a) * edge.conj()).re / length2).clamp(0.0, 1.0) };
        if (z - (a + edge * t)).norm() <= tol {
            return None;
        }
        let side = edge.re * (z.im - a.im) - (z.re - a.re) * edge.im;
        if a.im <= z.im {
            if b.im > z.im && side > 0.0 {
                winding += 1;
            }
        } else if b.im <= z.im && side < 0.0 {
            winding -= 1;
        }
    }
    Some(winding)
}

/// The winding number of `f` along the polygon around the origin: zeros
/// minus poles inside it. `None` if `f` vanishes or blows up on the
/// polygon, or the sampling cannot follow its argument.
fn polygon_count<F>(
    f: F,
    v: &[Complex64],
    per_edge: usize,
) -> Option<i64>
where
    F: Fn(Complex64) -> Complex64,
{
    let ok = |w: Complex64| w.is_finite() && w.norm() > 0.0;
    let mut previous = f(v[0]);
    if !ok(previous) {
        return None;
    }
    let mut winding = 0.0;
    for k in 0..v.len() {
        let (a, b) = (v[k], v[(k + 1) % v.len()]);
        for j in 1..=per_edge {
            #[allow(clippy::cast_precision_loss)]
            let t = j as f64 / per_edge as f64;
            let value = f(a + (b - a) * t);
            if !ok(value) {
                return None;
            }
            let step = (value / previous).arg();
            if step.abs() > std::f64::consts::FRAC_PI_2 {
                return None;
            }
            winding += step;
            previous = value;
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    Some((winding / std::f64::consts::TAU).round() as i64)
}

impl Analysis {
    /// The residue of `f` at the pole `root`, checked numerically.
    fn residue_at(
        cx: &mut Cx<'_>,
        f: NodeId,
        q: &Quotient,
        root: &Root,
        all: &[Root],
        z: NodeId,
    ) -> Option<NodeId> {
        let g = regular_part(cx, q, root, z);
        let residue = taylor_coefficient(cx, g, z, root.point, root.multiplicity - 1)?;
        let symbol = cx.graph.symbol_of(z)?;
        if let Some(function) = numeric_function(cx.graph, f, symbol) {
            let nearest = all
                .iter()
                .map(|other| (other.value - root.value).norm())
                .filter(|&d| d > 1e-12)
                .fold(f64::INFINITY, f64::min);
            let radius = (nearest / 2.0).min(0.5);
            let want = numeric::residue(&function, root.value, radius, 256);
            if let Some(got) = value_of(cx.graph, residue)
                && (got - want).norm() > 1e-6 * (1.0 + want.norm()) {
                    return None;
                }
        }
        Some(residue)
    }

    /// `2 pi I * sum of w(p) Res(p)` over the poles `p`, `w` the winding
    /// number of the contour.
    fn region_integral(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        contour: &Contour,
    ) -> Option<NodeId> {
        let q = quotient(cx.graph, self.ops, f, z)?;
        let poles = poles_of(cx, self.ops, &q)?;
        let mut residues = Vec::new();
        for root in &poles {
            // A pole on the contour: the integral does not exist.
            let weight = contour.winding(root.value)?;
            if weight != 0 {
                let residue = Self::residue_at(cx, f, &q, root, &poles, z)?;
                let weight = cx.graph.int(weight);
                residues.push(build::mul(cx.graph, &[weight, residue]));
            }
        }
        let sum = build::add(cx.graph, &residues);
        let two = cx.graph.int(2);
        let pi = cx.graph.ops().lookup("pi").map(|p| cx.graph.node(p, &[]))?;
        let unit = cx.graph.node(self.ops.unit, &[]);
        let value = build::mul(cx.graph, &[two, pi, unit, sum]);
        Some(cx.simplify(value))
    }

    #[allow(clippy::too_many_lines)]
    #[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let arg = |k: usize| args.get(k).copied();
        match self.request {
            | Request::Poles | Request::Zeros => {
                let (f, z) = (arg(0)?, arg(1)?);
                let q = quotient(cx.graph, self.ops, f, z)?;
                let mut entries = Vec::new();
                if self.request == Request::Poles {
                    let poles = poles_of(cx, self.ops, &q)?;
                    for root in poles {
                        let cancelled = vanishing_order(cx, q.numer, z, root.point, root.multiplicity);
                        let order = root.multiplicity - cancelled;
                        if order > 0 {
                            let m = cx.graph.int(i64::from(order));
                            entries.push(build::call(cx.graph, core::LIST, &[root.point, m]));
                        }
                    }
                } else {
                    // Zeros of a polynomial numerator.
                    let numer = Quotient { numer: cx.graph.int(1), ..polynomial_quotient(cx.graph, q.numer, z)? };
                    for root in poles_of(cx, self.ops, &numer)? {
                        let m = cx.graph.int(i64::from(root.multiplicity));
                        entries.push(build::call(cx.graph, core::LIST, &[root.point, m]));
                    }
                }
                Some(build::call(cx.graph, core::LIST, &entries))
            },
            | Request::Residue => {
                let (f, z, a) = (arg(0)?, arg(1)?, arg(2)?);
                let at = value_of(cx.graph, a);
                if let Some(q) = quotient(cx.graph, self.ops, f, z) {
                    let poles = poles_of(cx, self.ops, &q)?;
                    let at = at?;
                    return match poles.iter().find(|r| (r.value - at).norm() < 1e-9 * (1.0 + at.norm())) {
                        | Some(root) => Self::residue_at(cx, f, &q, &root.clone(), &poles, z),
                        | None => Some(cx.graph.int(0)),
                    };
                }
                // Otherwise: the Laurent expansion, if it exists.
                let (valuation, coefficients) = laurent_expansion(cx, f, z, a, -1)?;
                Some(if valuation > -1 {
                    cx.graph.int(0)
                } else {
                    *coefficients.get(usize::try_from(-1 - valuation).ok()?)?
                })
            },
            | Request::Singularity => {
                let (f, z, a) = (arg(0)?, arg(1)?, arg(2)?);
                if self.essential_at(cx, f, z, a) {
                    return Some(cx.graph.sym("essential"));
                }
                let (valuation, _) = laurent_expansion(cx, f, z, a, 0)?;
                if valuation < 0 {
                    let pole = cx.graph.sym("pole");
                    let m = cx.graph.int(-valuation);
                    return Some(build::call(cx.graph, core::LIST, &[pole, m]));
                }
                // Analytic in a punctured disc: regular if f is defined at
                // a, removable otherwise.
                let at = cx.graph.substitute(f, z, a);
                let at = best(cx.graph, at)?;
                let defined = value_of(cx.graph, at).is_some_and(num_complex::Complex::is_finite);
                Some(cx.graph.sym(if defined { "regular" } else { "removable" }))
            },
            | Request::Contour => {
                let (f, z, contour) = (arg(0)?, arg(1)?, arg(2)?);
                if let Some(c) = contour_of(cx.graph, contour) {
                    return self.region_integral(cx, f, z, &c);
                }
                self.path_integral(cx, f, z, contour)
            },
            | Request::CauchyIntegral | Request::CauchyDerivative => {
                // ∮ f(z)/(z - a)^(n+1) dz = 2 pi I f^(n)(a) / n!
                let (f, z, a) = (arg(0)?, arg(1)?, arg(2)?);
                let n = if self.request == Request::CauchyDerivative {
                    u32::try_from(cx.graph.number_of(arg(3)?)?.to_i64()?).ok()?
                } else {
                    0
                };
                let term = best(cx.graph, f)?;
                let coefficient = taylor_coefficient(cx, term, z, a, n)?;
                let two = cx.graph.int(2);
                let pi = cx.graph.ops().lookup("pi").map(|p| cx.graph.node(p, &[]))?;
                let unit = cx.graph.node(self.ops.unit, &[]);
                Some(build::mul(cx.graph, &[two, pi, unit, coefficient]))
            },
            | Request::CountZerosPoles => {
                let (f, z, contour) = (arg(0)?, arg(1)?, arg(2)?);
                let contour = contour_of(cx.graph, contour)?;
                // The argument principle, numerically: robust, and an
                // integer.
                let symbol = cx.graph.symbol_of(z)?;
                let winding = numeric_function(cx.graph, f, symbol).and_then(|function| match &contour {
                    | Contour::Circle { centre, radius } => numeric::count_zeros_poles(function, *centre, *radius, 4096),
                    | Contour::Polygon(vertices) => polygon_count(function, vertices, 2048),
                });
                let exact = self.count_exactly(cx, f, z, &contour);
                match (exact, winding) {
                    | (Some(a), Some(b)) if a != b => None,
                    | (Some(count), _) | (None, Some(count)) => Some(cx.graph.int(count)),
                    | (None, None) => None,
                }
            },
            | Request::Radius => {
                let (f, z, a) = (arg(0)?, arg(1)?, arg(2)?);
                let centre = value_of(cx.graph, a)?;
                let q = quotient(cx.graph, self.ops, f, z)?;
                let mut nearest: Option<(f64, NodeId)> = None;
                for root in poles_of(cx, self.ops, &q)? {
                    let cancelled = vanishing_order(cx, q.numer, z, root.point, root.multiplicity);
                    if cancelled >= root.multiplicity {
                        continue;
                    }
                    let d = (root.value - centre).norm();
                    if nearest.is_none_or(|(best, _)| d < best) {
                        nearest = Some((d, root.point));
                    }
                }
                match nearest {
                    | None => cx.graph.ops().lookup("oo").map(|oo| cx.graph.node(oo, &[])),
                    | Some((_, pole)) => {
                        let difference = build::sub(cx.graph, pole, a);
                        let modulus = build::call(cx.graph, self.ops.abs, &[difference]);
                        Some(cx.simplify(modulus))
                    },
                }
            },
            | Request::Continue => {
                // Each point must lie inside the disc of convergence at the
                // previous one; the continuation of a closed form is then
                // its expansion at the last point.
                let (f, z, points, order) = (arg(0)?, arg(1)?, arg(2)?, arg(3)?);
                let points = best(cx.graph, points)?;
                if cx.graph.op(points) != core::LIST {
                    return None;
                }
                let points = cx.graph.children(points).to_vec();
                let radius_op = cx.graph.ops().lookup("radius_of_convergence")?;
                for pair in points.windows(2) {
                    let radius = build::call(cx.graph, radius_op, &[f, z, pair[0]]);
                    let radius = cx.simplify(radius);
                    let step = value_of(cx.graph, pair[1])? - value_of(cx.graph, pair[0])?;
                    let r = if cx.graph.ops().get(cx.graph.op(radius)).name.as_ref() == "oo" {
                        f64::INFINITY
                    } else {
                        value_of(cx.graph, radius)?.re
                    };
                    if step.norm() >= r {
                        return None;
                    }
                }
                let last = *points.last()?;
                let taylor = cx.graph.ops().lookup("taylor")?;
                Some(build::call(cx.graph, taylor, &[f, z, last, order]))
            },
            | Request::InsideContour | Request::WindingNumber => {
                let (point, contour) = (arg(0)?, arg(1)?);
                let point = value_of(cx.graph, point).filter(|p| p.is_finite())?;
                let contour = contour_of(cx.graph, contour)?;
                let winding = contour.winding(point);
                if self.request == Request::WindingNumber {
                    return Some(cx.graph.int(winding?));
                }
                // Strictly inside: a point on the contour is not.
                Some(cx.graph.lit(Payload::Bool(winding.is_some_and(|w| w != 0))))
            },
            | Request::ResidueTheorem => {
                let (f, z, points) = (arg(0)?, arg(1)?, arg(2)?);
                let points = best(cx.graph, points)?;
                if cx.graph.op(points) != core::LIST {
                    return None;
                }
                let residue = cx.graph.ops().lookup("residue")?;
                let mut terms = Vec::new();
                for a in cx.graph.children(points).to_vec() {
                    let request = build::call(cx.graph, residue, &[f, z, a]);
                    let value = cx.simplify(request);
                    if cx.graph.op(value) == residue {
                        return None;
                    }
                    terms.push(value);
                }
                let sum = build::add(cx.graph, &terms);
                let two = cx.graph.int(2);
                let pi = cx.graph.ops().lookup("pi").map(|p| cx.graph.node(p, &[]))?;
                let unit = cx.graph.node(self.ops.unit, &[]);
                let value = build::mul(cx.graph, &[two, pi, unit, sum]);
                Some(cx.simplify(value))
            },
            | Request::CheckAnalytic => self.check_analytic(cx, arg(0)?, arg(1)?),
            | Request::BranchCut => self.branch_cut(cx, arg(0)?, arg(1)?),
            | Request::MobiusApply => {
                let m = mobius(cx.graph, arg(0)?)?;
                let z = arg(1)?;
                let az = build::mul(cx.graph, &[m[0], z]);
                let numerator = build::add(cx.graph, &[az, m[1]]);
                let cz = build::mul(cx.graph, &[m[2], z]);
                let denominator = build::add(cx.graph, &[cz, m[3]]);
                let inverse = build::powi(cx.graph, denominator, -1);
                Some(build::mul(cx.graph, &[numerator, inverse]))
            },
            | Request::MobiusCompose => {
                let (m, n) = (mobius(cx.graph, arg(0)?)?, mobius(cx.graph, arg(1)?)?);
                let g = &mut *cx.graph;
                let mut entry = |i: usize, j: usize| {
                    let first = build::mul(g, &[m[2 * i], n[j]]);
                    let second = build::mul(g, &[m[2 * i + 1], n[2 + j]]);
                    build::add(g, &[first, second])
                };
                let (a, b, c, d) = (entry(0, 0), entry(0, 1), entry(1, 0), entry(1, 1));
                Some(mobius_term(cx.graph, [a, b, c, d]))
            },
            | Request::MobiusInverse => {
                let m = mobius(cx.graph, arg(0)?)?;
                let (minus_b, minus_c) = (build::neg(cx.graph, m[1]), build::neg(cx.graph, m[2]));
                Some(mobius_term(cx.graph, [m[3], minus_b, minus_c, m[0]]))
            },
        }
    }

    /// Cauchy–Riemann for `f(x + I y)`.
    fn check_analytic(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let zs = cx.graph.symbol_of(z)?;
        let term = best(cx.graph, f)?;
        // Fresh real symbols for x, y and for every other free symbol.
        let fresh = |graph: &mut Graph, name: &str| {
            let symbol = graph.interner_mut().fresh_symbol(name);
            graph.assume(symbol, Facts::REAL);
            symbol
        };
        let (xs, ys) = (fresh(cx.graph, "x"), fresh(cx.graph, "y"));
        let (x, y) = (cx.graph.symbol_node(xs), cx.graph.symbol_node(ys));
        let mut term = term;
        let others: Vec<SymbolId> = cx.graph.free_symbols(cx.graph.find(term)).iter().copied().filter(|&s| s != zs).collect();
        let mut constants = Vec::new();
        for other in others {
            let replacement = fresh(cx.graph, "c");
            let (from, to) = (cx.graph.symbol_node(other), cx.graph.symbol_node(replacement));
            term = cx.graph.substitute(term, from, to);
            constants.push(replacement);
        }
        let unit = cx.graph.node(self.ops.unit, &[]);
        let iy = build::mul(cx.graph, &[unit, y]);
        let w = build::add(cx.graph, &[x, iy]);
        let g = cx.graph.substitute(term, z, w);
        let g = desugar(cx.graph, self.ops, g, &mut HashMap::new())?;
        let (u, v) = super::split(cx.graph, self.ops, g)?;
        let (u, v) = (cx.simplify(u), cx.simplify(v));
        let (ux, uy) = (derivative(cx.graph, u, x)?, derivative(cx.graph, u, y)?);
        let (vx, vy) = (derivative(cx.graph, v, x)?, derivative(cx.graph, v, y)?);
        let first = build::sub(cx.graph, ux, vy);
        let second = build::add(cx.graph, &[uy, vx]);
        let (first, second) = (cx.simplify(first), cx.simplify(second));
        let is_zero = |graph: &Graph, n: NodeId| graph.number_of(n).is_some_and(Number::is_zero);
        if is_zero(cx.graph, first) && is_zero(cx.graph, second) {
            return Some(cx.graph.lit(Payload::Bool(true)));
        }
        // Sample points: the differences must vanish everywhere.
        let mut valid = 0;
        for k in 0..6_u32 {
            let mut env = Env::numeric(0.0);
            env.bind(xs, 0.31 + 0.53 * f64::from(k) * if k % 2 == 0 { 1.0 } else { -1.0 });
            env.bind(ys, 0.47 - 0.29 * f64::from(k));
            for (j, &c) in constants.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                env.bind(c, 0.45 + 0.17 * j as f64);
            }
            let (Some(a), Some(b)) = (cx.graph.eval(first, &env), cx.graph.eval(second, &env)) else {
                continue;
            };
            if !a.is_finite() || !b.is_finite() {
                continue;
            }
            valid += 1;
            if a.abs() > 1e-7 || b.abs() > 1e-7 {
                return Some(cx.graph.lit(Payload::Bool(false)));
            }
        }
        (valid >= 2).then(|| cx.graph.lit(Payload::Bool(true)))
    }

    /// Cuts of the principal branches occurring in `f`.
    #[allow(clippy::too_many_lines)]
    fn branch_cut(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let term = best(cx.graph, f)?;
        let named = |graph: &Graph, name: &str| graph.ops().lookup(name);
        let (asin, acos, atan) = (named(cx.graph, "asin"), named(cx.graph, "acos"), named(cx.graph, "atan"));
        let (asinh, acosh, atanh) = (named(cx.graph, "asinh"), named(cx.graph, "acosh"), named(cx.graph, "atanh"));
        let (sqrt, pi) = (named(cx.graph, "sqrt"), named(cx.graph, "pi")?);
        let arg_op = named(cx.graph, "arg")?;
        // Cuts in the plane of the argument: (start, angle in units of
        // pi/2), the start as a multiple of 1 or I.
        let ln = self.ops.ln;
        let mut found: Vec<(NodeId, Vec<(Start, i64)>)> = Vec::new();
        let mut stack = vec![term];
        let mut seen = std::collections::HashSet::new();
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            let op = cx.graph.op(n);
            let children = cx.graph.children(n).to_vec();
            stack.extend_from_slice(&children);
            let cuts: Option<(NodeId, Vec<(Start, i64)>)> = if op == ln || Some(op) == sqrt {
                children.first().map(|&g| (g, vec![(Start::Zero, 2)]))
            } else if op == core::POW {
                let (&base, &exponent) = (children.first()?, children.get(1)?);
                let non_integer = cx.graph.number_of(exponent).is_some_and(|e| !e.is_integer());
                non_integer.then_some((base, vec![(Start::Zero, 2)]))
            } else if Some(op) == asin || Some(op) == acos || Some(op) == atanh {
                children.first().map(|&g| (g, vec![(Start::MinusOne, 2), (Start::One, 0)]))
            } else if Some(op) == atan || Some(op) == asinh {
                children.first().map(|&g| (g, vec![(Start::I, 1), (Start::MinusI, -1)]))
            } else if Some(op) == acosh {
                children.first().map(|&g| (g, vec![(Start::One, 2)]))
            } else {
                None
            };
            if let Some(entry) = cuts {
                found.push(entry);
            }
        }
        // Each argument must be a z + b.
        let mut rays: Vec<NodeId> = Vec::new();
        let ray = cx.graph.ops().lookup("ray")?;
        for (g, cuts) in found {
            let mut gens = Gens::default();
            let gz = gens.index(cx.graph, z);
            let poly = crate::rules::poly::repr::from_term(cx.graph, &mut gens, g, Limits::default())?;
            if poly.degree_in(gz) != 1 {
                return None;
            }
            let parts = poly.coefficients_in(gz);
            let b = to_term(cx.graph, &gens, parts.first()?);
            let a = to_term(cx.graph, &gens, parts.get(1)?);
            let depends = |graph: &Graph, n: NodeId| graph.symbol_of(z).is_some_and(|s| graph.depends_on(graph.find(n), s));
            if depends(cx.graph, a) || depends(cx.graph, b) {
                return None;
            }
            let arg_a = build::call(cx.graph, arg_op, &[a]);
            let arg_a = cx.simplify(arg_a);
            if cx.graph.op(arg_a) == arg_op {
                return None;
            }
            for (start, quarter_turns) in cuts {
                let unit = cx.graph.node(self.ops.unit, &[]);
                let s_g = match start {
                    | Start::Zero => cx.graph.int(0),
                    | Start::One => cx.graph.int(1),
                    | Start::MinusOne => cx.graph.int(-1),
                    | Start::I => unit,
                    | Start::MinusI => build::neg(cx.graph, unit),
                };
                // z = (s - b) / a, angle = theta - arg a
                let shifted = build::sub(cx.graph, s_g, b);
                let inverse = build::powi(cx.graph, a, -1);
                let start_z = build::mul(cx.graph, &[shifted, inverse]);
                let start_z = cx.simplify(start_z);
                let pi_node = cx.graph.node(pi, &[]);
                let fraction = cx.graph.num(Number::fraction(quarter_turns, 2)?);
                let theta = build::mul(cx.graph, &[fraction, pi_node]);
                let angle = build::sub(cx.graph, theta, arg_a);
                let angle = cx.simplify(angle);
                let node = build::call(cx.graph, ray, &[start_z, angle]);
                if !rays.iter().any(|&r| cx.graph.same(r, node)) {
                    rays.push(node);
                }
            }
        }
        // A deterministic order: by start point, then angle.
        let mut keyed: Vec<(NodeId, (f64, f64, f64))> = rays
            .into_iter()
            .map(|r| {
                let children = cx.graph.children(r).to_vec();
                let start = value_of(cx.graph, children[0]).unwrap_or_default();
                let angle = value_of(cx.graph, children[1]).map_or(0.0, |a| a.re);
                (r, (start.re, start.im, angle))
            })
            .collect();
        keyed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let rays: Vec<NodeId> = keyed.into_iter().map(|(r, _)| r).collect();
        Some(build::call(cx.graph, core::LIST, &rays))
    }

    /// Zeros minus poles inside a circle from the exact roots.
    fn count_exactly(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        contour: &Contour,
    ) -> Option<i64> {
        let q = quotient(cx.graph, self.ops, f, z)?;
        let numer = polynomial_quotient(cx.graph, q.numer, z)?;
        let mut count: i64 = 0;
        for (set, sign) in [(&numer, 1_i64), (&q, -1)] {
            for root in poles_of(cx, self.ops, set)? {
                count += sign * i64::from(root.multiplicity) * contour.winding(root.value)?;
            }
        }
        Some(count)
    }

    /// Whether `f` has an essential singularity at `a`: an entire
    /// transcendental function of an argument with a pole at `a`.
    fn essential_at(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        a: NodeId,
    ) -> bool {
        let Some(term) = best(cx.graph, f) else {
            return false;
        };
        let transcendental = [self.ops.exp, self.ops.sin, self.ops.cos, self.ops.sinh, self.ops.cosh];
        let mut stack = vec![term];
        while let Some(n) = stack.pop() {
            let children = cx.graph.children(n).to_vec();
            if transcendental.contains(&cx.graph.op(n))
                && let Some(&inner) = children.first()
                    && laurent_expansion(cx, inner, z, a, 0).is_some_and(|(v, _)| v < 0) {
                        return true;
                    }
            stack.extend(children);
        }
        false
    }

    /// `∫ f(g(t)) g'(t) dt` along `path(g, t, t0, t1)` with real `t`, as
    /// definite integrals of the real and imaginary parts.
    fn path_integral(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        contour: NodeId,
    ) -> Option<NodeId> {
        let term = best(cx.graph, contour)?;
        if cx.graph.ops().get(cx.graph.op(term)).name.as_ref() != "path" {
            return None;
        }
        let &[gamma, t, t0, t1] = cx.graph.children(term) else {
            return None;
        };
        // The parameter is real: integrate over a fresh real symbol.
        let fresh = cx.graph.interner_mut().fresh_symbol("t");
        cx.graph.assume(fresh, crate::graph::Facts::REAL);
        let s = cx.graph.symbol_node(fresh);
        let gamma = cx.graph.substitute(gamma, t, s);
        let f_term = best(cx.graph, f)?;
        let along = cx.graph.substitute(f_term, z, gamma);
        let velocity = derivative(cx.graph, gamma, s)?;
        let integrand = build::mul(cx.graph, &[along, velocity]);
        let (re, im) = super::split(cx.graph, self.ops, integrand)?;
        let (re, im) = (cx.simplify(re), cx.simplify(im));
        let defint = cx.graph.ops().lookup("defint")?;
        let re_integral = build::call(cx.graph, defint, &[re, s, t0, t1]);
        let im_integral = build::call(cx.graph, defint, &[im, s, t0, t1]);
        Some(build::complex(cx.graph, self.ops.unit, re_integral, im_integral))
    }
}

/// Prepares a term for splitting into real and imaginary parts: powers
/// with non-integer exponents become `exp(e ln b)` (the principal value)
/// and `re`, `im`, `conj`, `abs`, `arg` are replaced by their definitions.
fn desugar(
    graph: &mut Graph,
    ops: Ops,
    node: NodeId,
    memo: &mut HashMap<NodeId, Option<NodeId>>,
) -> Option<NodeId> {
    if let Some(known) = memo.get(&node) {
        return *known;
    }
    let children = graph.children(node).to_vec();
    let result = if children.is_empty() {
        Some(node)
    } else {
        let mut rebuilt = Vec::with_capacity(children.len());
        for c in &children {
            rebuilt.push(desugar(graph, ops, *c, memo)?);
        }
        let op = graph.op(node);
        if let (true, &[w]) = ([ops.re, ops.im, ops.conj, ops.abs, ops.arg].contains(&op), rebuilt.as_slice()) {
            let (a, b) = super::split(graph, ops, w)?;
            let minus_b = build::neg(graph, b);
            Some(match op {
                | _ if op == ops.re => a,
                | _ if op == ops.im => b,
                | _ if op == ops.conj => build::complex(graph, ops.unit, a, minus_b),
                | _ if op == ops.abs => {
                    let a2 = build::powi(graph, a, 2);
                    let b2 = build::powi(graph, b, 2);
                    let sum = build::add(graph, &[a2, b2]);
                    // sqrt(a^2 + b^2) as exp(ln(.)/2): the splitter reads the
                    // principal value of that.
                    let half = graph.num(Number::fraction(1, 2)?);
                    let ln = build::call(graph, ops.ln, &[sum]);
                    let product = build::mul(graph, &[half, ln]);
                    build::call(graph, ops.exp, &[product])
                },
                | _ => build::call(graph, ops.atan2, &[b, a]),
            })
        } else if graph.ops().lookup("sqrt") == Some(op) {
            let &[base] = rebuilt.as_slice() else {
                return None;
            };
            let half = graph.num(Number::fraction(1, 2)?);
            let ln = build::call(graph, ops.ln, &[base]);
            let product = build::mul(graph, &[half, ln]);
            Some(build::call(graph, ops.exp, &[product]))
        } else if op == core::POW {
            let &[base, exponent] = rebuilt.as_slice() else {
                return None;
            };
            if graph.number_of(exponent).is_some_and(Number::is_integer) {
                graph.try_node(op, &rebuilt)
            } else {
                let ln = build::call(graph, ops.ln, &[base]);
                let product = build::mul(graph, &[exponent, ln]);
                Some(build::call(graph, ops.exp, &[product]))
            }
        } else {
            graph.try_node(op, &rebuilt)
        }
    };
    memo.insert(node, result);
    result
}

/// A polynomial in `z` with rational coefficients, as a quotient with
/// numerator one (its roots are then its "poles").
fn polynomial_quotient(
    graph: &mut Graph,
    p: NodeId,
    z: NodeId,
) -> Option<Quotient> {
    let mut gens = Gens::default();
    let gz = gens.index(graph, z);
    let r = ratio(graph, &mut gens, p, Limits::default())?;
    if r.numer.support().iter().any(|&g| g != gz) || r.denom.as_constant().is_none() {
        return None;
    }
    let numer: QPoly = r.numer.univariate_in(gz)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
    let (lead, factors) = univariate::factor(&numer);
    let factors = factors.into_iter().map(|(f, m)| (f.into_iter().map(BigRational::from_integer).collect(), m)).collect();
    Some(Quotient { numer: graph.int(1), lead, factors })
}

fn mobius(
    graph: &mut Graph,
    node: NodeId,
) -> Option<[NodeId; 4]> {
    let term = best(graph, node)?;
    let rows = graph.children(term).to_vec();
    if graph.op(term) != core::LIST || rows.len() != 2 {
        return None;
    }
    let (r0, r1) = (rows[0], rows[1]);
    if graph.op(r0) != core::LIST || graph.op(r1) != core::LIST {
        return None;
    }
    let (&[a, b], &[c, d]) = (graph.children(r0), graph.children(r1)) else {
        return None;
    };
    Some([a, b, c, d])
}

fn mobius_term(
    graph: &mut Graph,
    [a, b, c, d]: [NodeId; 4],
) -> NodeId {
    let top = build::call(graph, core::LIST, &[a, b]);
    let bottom = build::call(graph, core::LIST, &[c, d]);
    build::call(graph, core::LIST, &[top, bottom])
}

#[cfg(test)]
mod tests {
    use super::super::complex;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[complex()], src)
    }

    #[test]
    fn poles_and_zeros() {
        assert_eq!(run("poles(1/(z^2 + 1), z)"), "list(list(I, 1), list(-I, 1))");
        assert_eq!(run("poles(z/(z^3 - z^2), z)"), "list(list(1, 1), list(0, 1))");
        assert_eq!(run("poles(sin(z)/(z - 1)^2, z)"), "list(list(1, 2))");
        assert_eq!(run("zeros_of(z^2 - 4, z)"), "list(list(2, 1), list(-2, 1))");
    }

    #[test]
    fn residues() {
        assert_eq!(run("residue(1/(z^2 + 1), z, I)"), "-1/2*I");
        assert_eq!(run("residue(exp(z)/z^3, z, 0)"), "1/2");
        assert_eq!(run("residue(1/(z - 1)^2, z, 1)"), "0");
        assert_eq!(run("residue(z^2/(z - 2), z, 2)"), "4");
        assert_eq!(run("residue(1/(z^2 + 1), z, 5)"), "0");
        assert_eq!(run("residue(cos(z)/z, z, 0)"), "1");
    }

    #[test]
    fn contour_integrals() {
        assert_eq!(run("contour_integral(1/z, z, circle(0, 1))"), "2*I*pi");
        assert_eq!(run("contour_integral(1/(z^2 + 1), z, circle(I, 1))"), "pi");
        assert_eq!(run("contour_integral(1/(z^2 + 1), z, circle(0, 2))"), "0");
        assert_eq!(run("contour_integral(exp(z)/z^2, z, circle(0, 1))"), "2*I*pi");
        assert_eq!(run("contour_integral(z^2, z, circle(0, 1))"), "0");
        // A parametrised path: the unit circle, explicitly.
        assert_eq!(run("contour_integral(1/z, z, path(cos(t) + I*sin(t), t, 0, 2*pi))"), "2*I*pi");
    }

    #[test]
    fn cauchy_formulas() {
        assert_eq!(run("cauchy_integral(exp(z), z, 0)"), "2*I*pi");
        assert_eq!(run("cauchy_derivative(z^3, z, 1, 2)"), "6*I*pi");
    }

    #[test]
    fn counting_and_radius() {
        assert_eq!(run("count_zeros_poles((z - 1/2)*(z - 3)/(z^2 + 1/4), z, circle(0, 1))"), "-1");
        assert_eq!(run("count_zeros_poles(z^5 - 1, z, circle(0, 2))"), "5");
        assert_eq!(run("radius_of_convergence(1/(z^2 + 1), z, 0)"), "1");
        assert_eq!(run("radius_of_convergence(1/(1 - z), z, 3)"), "2");
        assert_eq!(run("radius_of_convergence(exp(z), z, 0)"), "oo");
        assert_eq!(run("distance(3*I, 4)"), "5");
    }

    #[test]
    fn singularities() {
        assert_eq!(run("singularity(1/z^2, z, 0)"), "list(pole, 2)");
        assert_eq!(run("singularity(sin(z)/z, z, 0)"), "removable");
        assert_eq!(run("singularity(exp(1/z), z, 0)"), "essential");
        assert_eq!(run("singularity(exp(z), z, 0)"), "regular");
    }

    #[test]
    fn mobius_maps() {
        assert_eq!(run("mobius_apply(list(list(1, 1), list(0, 1)), z)"), "z + 1");
        assert_eq!(run("mobius_apply(mobius_inverse(list(list(2, 1), list(1, 1))), mobius_apply(list(list(2, 1), list(1, 1)), 3))"), "3");
        assert_eq!(run("mobius_compose(list(list(1, 1), list(0, 1)), list(list(1, 2), list(0, 1)))"), "list(list(1, 3), list(0, 1))");
    }

    #[test]
    fn continuation() {
        let text = run("continue_along(1/(1 - z), z, list(0, 1/2*I, I, 3/2*I), 2)");
        assert!(text.contains('I'), "{text}");
    }

    #[test]
    fn analyticity_by_cauchy_riemann() {
        for f in ["z^2", "exp(z)", "z*sin(z)", "1/z", "(z^2 + 1)/(z - 1)", "cos(z)^2 + 3", "a*z^2 + b", "sqrt(z)", "ln(z)", "z^(1/3)", "5", "z"] {
            assert_eq!(run(&format!("check_analytic({f}, z)")), "true", "{f}");
        }
        for f in ["conj(z)", "abs(z)", "re(z)", "im(z)", "z*conj(z)", "abs(z)^2", "z + conj(z)", "exp(conj(z))", "arg(z)", "z/conj(z)"] {
            assert_eq!(run(&format!("check_analytic({f}, z)")), "false", "{f}");
        }
        // other variable names
        assert_eq!(run("check_analytic(w^2, w)"), "true");
        assert_eq!(run("check_analytic(conj(w), w)"), "false");
        // Unknown functions cannot be decided.
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "check_analytic(g(z), z)", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn points_inside_contours() {
        assert_eq!(run("is_inside_contour(1/2, circle(0, 1))"), "true");
        assert_eq!(run("is_inside_contour(2, circle(0, 1))"), "false");
        assert_eq!(run("is_inside_contour(1, circle(0, 1))"), "false");
        assert_eq!(run("is_inside_contour(1 + I, circle(1 + I, 1/2))"), "true");
        assert_eq!(run("is_inside_contour(I/2 + 1/2, circle(0, 1))"), "true");
        assert_eq!(run("is_inside_contour(sqrt(2)/2 + sqrt(2)/2*I, circle(0, 1))"), "false");
        let square = "polygon(list(0, 1, 1 + I, I))";
        let clockwise = "polygon(list(0, I, 1 + I, 1))";
        for contour in [square, clockwise] {
            assert_eq!(run(&format!("is_inside_contour(1/2 + I/2, {contour})")), "true");
            assert_eq!(run(&format!("is_inside_contour(2, {contour})")), "false");
            assert_eq!(run(&format!("is_inside_contour(-1/2 + I/2, {contour})")), "false");
            assert_eq!(run(&format!("is_inside_contour(1/2 + 2*I, {contour})")), "false");
            // vertices and edges are not strictly inside
            assert_eq!(run(&format!("is_inside_contour(0, {contour})")), "false");
            assert_eq!(run(&format!("is_inside_contour(1/2, {contour})")), "false");
            assert_eq!(run(&format!("is_inside_contour(1 + I/2, {contour})")), "false");
        }
        // A non-convex L: the notch x > 2, y > 2 is outside.
        let ell = "polygon(list(0, 4, 4 + 2*I, 2 + 2*I, 2 + 4*I, 4*I))";
        assert_eq!(run(&format!("is_inside_contour(1 + I, {ell})")), "true");
        assert_eq!(run(&format!("is_inside_contour(3 + 3*I, {ell})")), "false");
        assert_eq!(run(&format!("is_inside_contour(1 + 3*I, {ell})")), "true");
        assert_eq!(run(&format!("is_inside_contour(3 + I, {ell})")), "true");
        assert_eq!(run(&format!("is_inside_contour(5 + I, {ell})")), "false");
        // a triangle
        assert_eq!(run("is_inside_contour(1 + I/2, polygon(list(0, 3, 3*I)))"), "true");
        assert_eq!(run("is_inside_contour(2 + 2*I, polygon(list(0, 3, 3*I)))"), "false");
        // a symbolic point is undecided
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "is_inside_contour(a, circle(0, 1))", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn winding_numbers() {
        assert_eq!(run("winding_number(0, circle(0, 1))"), "1");
        assert_eq!(run("winding_number(3, circle(0, 1))"), "0");
        assert_eq!(run("winding_number(1/2 + I/2, polygon(list(0, 1, 1 + I, I)))"), "1");
        assert_eq!(run("winding_number(1/2 + I/2, polygon(list(0, I, 1 + I, 1)))"), "-1");
        assert_eq!(run("winding_number(5, polygon(list(0, 1, 1 + I, I)))"), "0");
        // traversed twice
        assert_eq!(run("winding_number(1/2 + I/2, polygon(list(0, 1, 1 + I, I, 0, 1, 1 + I, I)))"), "2");
        // a point on the contour has no winding number
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "winding_number(1/2, polygon(list(0, 1, 1 + I, I)))", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn residue_theorem_over_polygons() {
        let square = "polygon(list(-1 - I, 1 - I, 1 + I, -1 + I))";
        let clockwise = "polygon(list(-1 - I, -1 + I, 1 + I, 1 - I))";
        assert_eq!(run(&format!("contour_integral(1/z, z, {square})")), "2*I*pi");
        assert_eq!(run(&format!("contour_integral(1/z, z, {clockwise})")), "-2*I*pi");
        assert_eq!(run(&format!("contour_integral(exp(z)/z^2, z, {square})")), "2*I*pi");
        assert_eq!(run(&format!("contour_integral(z^3 + 2*z, z, {square})")), "0");
        // 1/(z^2 + 1) has poles at ±I; a rectangle around the upper one
        assert_eq!(run("contour_integral(1/(z^2 + 1), z, polygon(list(-1, 1, 1 + 2*I, -1 + 2*I)))"), "pi");
        assert_eq!(run("contour_integral(1/(z^2 + 1), z, polygon(list(-2 - 2*I, 2 - 2*I, 2 + 2*I, -2 + 2*I)))"), "0");
        // a triangle around 1 and 2 but not 0
        assert_eq!(run("contour_integral(1/(z*(z - 1)*(z - 2)), z, polygon(list(1/2 - I, 3 - I, 1 + 3*I)))"), "-I*pi");
        // twice around
        assert_eq!(run("contour_integral(1/z, z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I, -1 - I, 1 - I, 1 + I, -1 + I)))"), "4*I*pi");
        // the circle agrees with the polygon
        assert_eq!(run("contour_integral(1/(z - 1/2), z, circle(0, 1))"), run("contour_integral(1/(z - 1/2), z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I)))"));
        // a pole on the contour: undefined
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "contour_integral(1/z, z, polygon(list(0, 1, 1 + I)))", &[]);
        assert!(!reduced, "{text}");
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "contour_integral(1/(z - 1), z, polygon(list(0, 2, 2*I)))", &[]);
        assert!(!reduced, "{text}");
        // explicit residue sums
        assert_eq!(run("contour_integral_residue_theorem(1/(z^2 + 1), z, list(I))"), "pi");
        assert_eq!(run("contour_integral_residue_theorem(1/(z^2 + 1), z, list(I, -I))"), "0");
        assert_eq!(run("contour_integral_residue_theorem(exp(z)/z^2, z, list(0))"), "2*I*pi");
        assert_eq!(run("contour_integral_residue_theorem(z/(z - 1)/(z - 2), z, list(1, 2))"), "2*I*pi");
        assert_eq!(run("contour_integral_residue_theorem(1/z, z, list())"), "0");
    }

    #[test]
    fn argument_principle_over_polygons() {
        let square = "polygon(list(-2 - 2*I, 2 - 2*I, 2 + 2*I, -2 + 2*I))";
        assert_eq!(run(&format!("count_zeros_poles(z^5 - 1, z, {square})")), "5");
        assert_eq!(run("count_zeros_poles(z^5 - 1, z, polygon(list(-1/2 - I/2, 1/2 - I/2, 1/2 + I/2, -1/2 + I/2)))"), "0");
        assert_eq!(run("count_zeros_poles((z - 1/2)*(z - 3)/(z^2 + 1/4), z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I)))"), "-1");
        assert_eq!(run("count_zeros_poles((z - 1/2)*(z - 3)/(z^2 + 1/4), z, polygon(list(-1 - I, -1 + I, 1 + I, 1 - I)))"), "1");
        assert_eq!(run("count_zeros_poles(1/z^3, z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I)))"), "-3");
        // transcendental numerator: exp has no zeros
        assert_eq!(run("count_zeros_poles(exp(z)/(z - 1/2)^2, z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I)))"), "-2");
        // a zero on the contour is undefined
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "count_zeros_poles(z - 1, z, polygon(list(-1 - I, 1 - I, 1 + I, -1 + I)))", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn branch_cuts() {
        assert_eq!(run("branch_cut(ln(z), z)"), "list(ray(0, pi))");
        assert_eq!(run("branch_cut(sqrt(z - 1), z)"), "list(ray(1, pi))");
        assert_eq!(run("branch_cut(z^(1/3) + 1, z)"), "list(ray(0, pi))");
        assert_eq!(run("branch_cut(ln(2*z + 4), z)"), "list(ray(-2, pi))");
        assert_eq!(run("branch_cut(ln(I*z), z)"), "list(ray(0, 1/2*pi))");
        assert_eq!(run("branch_cut(ln(-z), z)"), "list(ray(0, 0))");
        assert_eq!(run("branch_cut(z^2 + exp(z) + 1/z, z)"), "list()");
        assert_eq!(run("branch_cut(z^(-2), z)"), "list()");
        assert_eq!(run("branch_cut(acosh(z), z)"), "list(ray(1, pi))");
        assert_eq!(run("branch_cut(asin(z), z)"), "list(ray(-1, pi), ray(1, 0))");
        assert_eq!(run("branch_cut(atan(z), z)"), "list(ray(-I, -1/2*pi), ray(I, 1/2*pi))");
        // a sum collects the cuts of both terms
        assert_eq!(run("branch_cut(sqrt(z) + ln(z - 2), z)"), "list(ray(0, pi), ray(2, pi))");
        assert_eq!(run("branch_cut(ln(z - 2) + sqrt(z), z)"), "list(ray(0, pi), ray(2, pi))");
        // a non-linear argument is not handled
        let (text, reduced) = crate::rules::testing::reduce_with(&[complex()], "branch_cut(ln(z^2 + 1), z)", &[]);
        assert!(!reduced, "{text}");
    }

    /// The principal value really jumps across the reported cut and is
    /// continuous across the plane elsewhere.
    #[test]
    fn reported_cuts_are_where_the_function_jumps() {
        use std::collections::HashMap;

        use num_complex::Complex64;

        use crate::graph::Engine;
        use crate::graph::Graph;
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[complex()]).unwrap_or_else(|e| panic!("{e}"));
        let _ = engine;
        let z = g.interner_mut().symbol("z");
        for (f, start, angle) in [
            ("ln(z)", Complex64::new(0.0, 0.0), std::f64::consts::PI),
            ("sqrt(z - 1)", Complex64::new(1.0, 0.0), std::f64::consts::PI),
            ("ln(2*z + 4)", Complex64::new(-2.0, 0.0), std::f64::consts::PI),
            ("ln(-z)", Complex64::new(0.0, 0.0), 0.0),
        ] {
            let node = g.parse(f).unwrap_or_else(|e| panic!("{e}"));
            let at = |g: &Graph, p: Complex64| g.eval_complex(node, &HashMap::from([(z, p)])).unwrap_or(Complex64::new(f64::NAN, f64::NAN));
            let direction = Complex64::from_polar(1.0, angle);
            let normal = direction * Complex64::new(0.0, 1.0);
            let on_cut = start + direction * 1.7;
            let jump = (at(&g, on_cut + normal * 1e-9) - at(&g, on_cut - normal * 1e-9)).norm();
            assert!(jump > 1.0, "{f}: jump {jump} across the cut");
            // across the opposite ray the function is continuous
            let off_cut = start - direction * 1.7;
            let smooth = (at(&g, off_cut + normal * 1e-9) - at(&g, off_cut - normal * 1e-9)).norm();
            assert!(smooth < 1e-6, "{f}: jump {smooth} off the cut");
        }
    }
}
