//! Integral and discrete transforms.
//!
//! | operator | definition |
//! |---|---|
//! | `laplace(f, t, s)` | `∫_0^oo f(t) exp(-s t) dt` |
//! | `inverse_laplace(F, s, t)` | the causal `f` with `laplace(f) = F` |
//! | `fourier(f, t, w)` | `∫_-oo^oo f(t) exp(-I w t) dt` |
//! | `inverse_fourier(F, w, t)` | `1/(2 pi) ∫ F(w) exp(I w t) dw` |
//! | `ztransform(f, n, z)` | `sum_{n >= 0} f(n) z^-n` (unilateral) |
//! | `inverse_ztransform(F, z, n)` | the sequence with that transform |
//! | `convolve(f, g, t)` | `∫_0^t f(u) g(t - u) du` (functions on `t >= 0`) |
//! | `convolution(f, g, t)` | `∫_-oo^oo f(u) g(t - u) du`; when both `f` and `g` carry a factor `heaviside(t)` (causal signals) it is `heaviside(t) ∫_0^t f0(u) g0(t - u) du` with the factors stripped |
//! | `discrete_convolution(f, g, n)` | `sum_{k=0}^n f(k) g(n - k)` (causal sequences) |
//! | `initial_value_laplace(F, s)`, `final_value_laplace(F, s)` | `f(0+) = lim_{s -> oo} s F(s)`, `f(oo) = lim_{s -> 0} s F(s)` |
//! | `initial_value_z(F, z)`, `final_value_z(F, z)` | `x(0) = lim_{z -> oo} F(z)`, `x(oo) = lim_{z -> 1} (z - 1) F(z)` |
//! | `dirac(x)`, `kronecker(n)` | the delta distribution; the Kronecker delta `[n = 0]` |
//!
//! Forward transforms are computed structurally: linearity, tables of
//! elementary transforms, and the shift, modulation, multiplication-by-`t`
//! and derivative theorems. Inverse transforms of rational functions use
//! exact partial fractions over `Q` (any denominator that factors into
//! linear and quadratic factors), and single linear or quadratic factors
//! with symbolic coefficients. Every result that can be evaluated is
//! checked numerically against the defining integral or sum before it is
//! accepted.
//!
//! # Unknown functions and the theorems
//!
//! The transform of an *unknown* function `y(t)` is left as the inert
//! request `laplace(y(t), t, s)` (a request that cannot be reduced further),
//! and the theorems act on it, so transformed differential equations can be
//! solved for the unknown transform. With `Y = laplace(y(t), t, s)`,
//! `F = fourier(y(t), t, w)` and `X = ztransform(y(n), n, z)`:
//!
//! | theorem | result |
//! |---|---|
//! | Laplace derivative | `laplace(diff(y(t), t), t, s) = s Y - y(0)`, and for order `n`: `s^n Y - sum_{k<n} s^(n-1-k) y^(k)(0)`, the initial values `y^(k)(0)` written `at(diff(...), t, 0)` |
//! | Laplace integral | `laplace(defint(y(u), u, 0, t), t, s) = Y / s` |
//! | Laplace frequency shift | `laplace(exp(a t) y(t), t, s) = Y(s - a)` |
//! | Laplace time shift (`c >= 0`) | `laplace(heaviside(t - c) y(t - c), t, s) = exp(-c s) Y` |
//! | Laplace scaling (`a > 0`) | `laplace(y(a t), t, s) = Y(s / a) / a` |
//! | Laplace multiplication by `t^n` | `laplace(t^n y(t), t, s) = (-1)^n d^n Y / ds^n` |
//! | Fourier derivative | `fourier(diff(y(t), t), t, w) = I w F` (and `(I w)^n F` for order `n`) |
//! | Fourier scaling and shift | `fourier(y(a t + b), t, w) = exp(I w b / a) F(w / a) / abs(a)` |
//! | Fourier frequency shift | `fourier(exp(I a t) y(t), t, w) = F(w - a)` |
//! | Fourier multiplication by `t^n` | `fourier(t^n y(t), t, w) = (I d/dw)^n F` |
//! | z time shift, delay (`k > 0`) | `ztransform(y(n - k), n, z) = z^-k X` (causal `y`: `y(m) = 0` for `m < 0`) |
//! | z time shift, advance (`k > 0`) | `ztransform(y(n + k), n, z) = z^k X - sum_{j<k} y(j) z^(k-j)` |
//! | z scaling, multiplication by `n` | `ztransform(a^n y(n), n, z) = X(z / a)`, `ztransform(n y(n), n, z) = -z dX/dz` |
//! | convolution theorems | `fourier(convolution(f, g, t), t, w)`, `laplace(convolution(f, g, t), t, s)` (and `convolve`) and `ztransform(discrete_convolution(f, g, n), n, z)` are the products of the transforms of `f` and `g`; a factor that cannot be transformed stays an inert transform |
//!
//! The inverses are the converse: `inverse_laplace(Y G, s, t)` with inert
//! `Y = laplace(y(t), t, s)`, `G = laplace(g(t), t, s)` is `convolution(y(t),
//! g(t), t)` (likewise for Fourier, and `discrete_convolution` for z), and
//! `inverse_fourier(convolution(F, G, w), w, t)` is `2 pi` times the
//! product of the inverse transforms.
//!
//! # Conventions
//!
//! The Laplace and z transforms are unilateral: a function is taken on
//! `t >= 0` (`n >= 0`), i.e. as causal. The Laplace derivative theorem uses
//! the initial values at `0+`; a Laplace time shift needs the explicit
//! `heaviside`. `convolution` is the two-sided integral, and equals the
//! causal `∫_0^t` form when both arguments are multiplied by
//! `heaviside(t)`; the Laplace and z convolution theorems assume causal
//! arguments (for non-causal arguments they apply to the restrictions to
//! `t >= 0`). The Fourier transform of a derivative assumes that the
//! function vanishes at infinity.
//!
//! # Partial fractions
//!
//! `inverse_laplace` expands a rational function exactly over `Q`:
//! linear factors of any power, quadratic factors of any power (the
//! squares and higher powers by differentiation of the simple case with
//! respect to the quadratic's parameter), a polynomial part that is
//! constant (a `dirac`). `inverse_ztransform` expands `F(z)/z`: linear
//! factors of any power, quadratics with complex or real irrational roots.

use std::collections::HashMap;

use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

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
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::integrate::gauss_kronrod_any;
use crate::rules::calculus::calculus;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::call;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::complex::complex;
use crate::rules::poly::apart::apart;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::pde::pde;
use crate::rules::poly::univariate::QPoly;
use crate::rules::special::special;

/// The transforms rule set.
#[must_use]
pub fn transforms() -> RuleSet {
    RuleSet::new("transforms", install).needs(calculus()).needs(complex()).needs(special()).needs(pde())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Laplace,
    InverseLaplace,
    Fourier,
    InverseFourier,
    Z,
    InverseZ,
    Convolve,
    Convolution,
    DiscreteConvolution,
    InitialLaplace,
    FinalLaplace,
    InitialZ,
    FinalZ,
}

/// Operators the transform kernels build with.
#[derive(Copy, Clone, Debug)]
struct Ops {
    laplace: OpId,
    exp: OpId,
    sin: OpId,
    cos: OpId,
    sinh: OpId,
    cosh: OpId,
    abs: OpId,
    heaviside: OpId,
    dirac: OpId,
    kronecker: OpId,
    gamma: OpId,
    unit: OpId,
    pi: OpId,
    diff: OpId,
    apply: OpId,
    defint: OpId,
    fourier: OpId,
    inverse_fourier: OpId,
    ztransform: OpId,
    at: OpId,
    convolve: OpId,
    convolution: OpId,
    dconvolution: OpId,
    oo: OpId,
    sum: OpId,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let dirac = i.op(OpDescriptor::new("dirac", Arity::Fixed(1)))?;
    let kronecker = i.op(OpDescriptor::new("kronecker", Arity::Fixed(1)).eval(|a| match a {
        | [n] if *n == 0.0 => 1.0,
        | [n] if n.is_finite() => 0.0,
        | _ => f64::NAN,
    }))?;
    i.rewrites(Tier::Normalize, &["transforms/kronecker: kronecker(?n) => 0 if integer(?n), nonzero(?n)"])?;
    let request = |name: &str, arity: u8, binder: Option<(u8, u32)>| {
        let desc = OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
        match binder {
            | Some((var, scope)) => desc.binder(var, scope),
            | None => desc,
        }
    };
    let table = [
        ("laplace", Request::Laplace, 3, Some((1, 0b1))),
        ("inverse_laplace", Request::InverseLaplace, 3, Some((1, 0b1))),
        ("fourier", Request::Fourier, 3, Some((1, 0b1))),
        ("inverse_fourier", Request::InverseFourier, 3, Some((1, 0b1))),
        ("ztransform", Request::Z, 3, Some((1, 0b1))),
        ("inverse_ztransform", Request::InverseZ, 3, Some((1, 0b1))),
        ("convolve", Request::Convolve, 3, None),
        ("convolution", Request::Convolution, 3, None),
        ("discrete_convolution", Request::DiscreteConvolution, 3, None),
        ("initial_value_laplace", Request::InitialLaplace, 2, Some((1, 0b1))),
        ("final_value_laplace", Request::FinalLaplace, 2, Some((1, 0b1))),
        ("initial_value_z", Request::InitialZ, 2, Some((1, 0b1))),
        ("final_value_z", Request::FinalZ, 2, Some((1, 0b1))),
    ];
    let mut registered = Vec::new();
    for (name, kind, arity, binder) in table {
        registered.push((i.op(request(name, arity, binder))?, kind));
    }
    let get = |i: &mut Installer<'_>, name: &str| {
        i.graph().ops().lookup(name).ok_or_else(|| RuleError::Invalid { rule: format!("transforms/{name}"), reason: "missing operator" })
    };
    let ops = Ops {
        laplace: registered[0].0,
        exp: get(i, "exp")?,
        sin: get(i, "sin")?,
        cos: get(i, "cos")?,
        sinh: get(i, "sinh")?,
        cosh: get(i, "cosh")?,
        abs: get(i, "abs")?,
        heaviside: get(i, "heaviside")?,
        dirac,
        kronecker,
        gamma: get(i, "gamma")?,
        unit: get(i, "I")?,
        pi: get(i, "pi")?,
        diff: get(i, "diff")?,
        apply: core::APPLY,
        defint: get(i, "defint")?,
        fourier: registered[2].0,
        inverse_fourier: registered[3].0,
        ztransform: registered[4].0,
        at: get(i, "at")?,
        convolve: registered[6].0,
        convolution: registered[7].0,
        dconvolution: registered[8].0,
        oo: get(i, "oo")?,
        sum: get(i, "sum")?,
    };
    for (op, request) in registered {
        let name = i.graph().ops().get(op).name.to_string();
        i.kernel(&format!("transforms/{name}"), Tier::Reduce, Transform { op, request, ops });
    }
    Ok(())
}

struct Transform {
    op: OpId,
    request: Request,
    ops: Ops,
}

impl Kernel for Transform {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let children = cx.graph.children(node).to_vec();
        match (self.request, children.as_slice()) {
            | (Request::Convolve | Request::Convolution | Request::DiscreteConvolution, &[f, g, t]) => {
                return convolution_request(cx, self.ops, self.request, f, g, t).map_or(Outcome::Pass, Outcome::Equal);
            },
            | (Request::InitialLaplace | Request::FinalLaplace | Request::InitialZ | Request::FinalZ, &[f, v]) => {
                return value_theorem(cx, self.ops, self.request, f, v).map_or(Outcome::Pass, Outcome::Equal);
            },
            | _ => {},
        }
        let &[f, x, y] = children.as_slice() else {
            return Outcome::Pass;
        };
        let Some(f) = best(cx.graph, f) else {
            return Outcome::Pass;
        };
        let (Some(xs), Some(_)) = (cx.graph.symbol_of(x), cx.graph.symbol_of(y)) else {
            return Outcome::Pass;
        };
        let mut t = Tx { cx, ops: self.ops, kind: self.request, x, xs, y };
        let result = match self.request {
            | Request::Laplace => t.laplace(f).filter(|&r| t.check_laplace(f, r)),
            | Request::InverseLaplace => t.inverse_laplace(f).filter(|&r| t.check_inverse_laplace(f, r)),
            | Request::Fourier => t.fourier(f).filter(|&r| t.check_fourier(f, r)),
            | Request::InverseFourier => t.inverse_fourier(f),
            | Request::Z => t.ztransform(f).filter(|&r| t.check_z(f, r)).map(|r| t.normal(r)),
            | Request::InverseZ => t.inverse_z(f).filter(|&r| t.check_inverse_z(f, r)),
            | _ => None,
        };
        // A transform that is its own answer (an unknown function) stays a
        // request.
        let result = result.filter(|&r| !t.cx.graph.same(r, node));
        let Some(r) = result else {
            return Outcome::Pass;
        };
        let simplified = t.cx.simplify(r);
        // A theorem applied to an unknown function leaves requests behind;
        // it is still the answer.
        if has_heavy(t.cx.graph, simplified) { Outcome::Pinned(simplified) } else { Outcome::Equal(simplified) }
    }
}

/// Whether the best form of `node` still contains a request.
fn has_heavy(
    graph: &mut Graph,
    node: NodeId,
) -> bool {
    let Some(term) = best(graph, node) else {
        return true;
    };
    let mut stack = vec![term];
    while let Some(n) = stack.pop() {
        if graph.ops().get(graph.op(n)).flags.has(OpFlags::HEAVY) {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// `convolve`, `convolution` and `discrete_convolution`: the defining
/// integral or sum, if it can be evaluated.
fn convolution_request(
    cx: &mut Cx<'_>,
    ops: Ops,
    request: Request,
    f: NodeId,
    g: NodeId,
    t: NodeId,
) -> Option<NodeId> {
    cx.graph.symbol_of(t)?;
    let (f, g) = (best(cx.graph, f)?, best(cx.graph, g)?);
    let fresh = cx.graph.interner_mut().fresh_symbol(if request == Request::DiscreteConvolution { "k" } else { "u" });
    let u = cx.graph.symbol_node(fresh);
    let t_minus_u = sub(cx.graph, t, u);
    let integrand = |cx: &mut Cx<'_>, f: NodeId, g: NodeId| {
        let f_u = cx.graph.substitute(f, t, u);
        let g_shift = cx.graph.substitute(g, t, t_minus_u);
        mul(cx.graph, &[f_u, g_shift])
    };
    let zero = cx.graph.int(0);
    let result = match request {
        | Request::Convolve => {
            let body = integrand(cx, f, g);
            call(cx.graph, ops.defint, &[body, u, zero, t])
        },
        | Request::DiscreteConvolution => {
            let body = integrand(cx, f, g);
            call(cx.graph, ops.sum, &[body, u, zero, t])
        },
        | _ => {
            // Causal signals: both carry heaviside(t).
            let strip = |cx: &mut Cx<'_>, h: NodeId| -> Option<NodeId> {
                let factors = if cx.graph.op(h) == core::MUL { cx.graph.children(h).to_vec() } else { vec![h] };
                let k = factors.iter().position(|&p| {
                    cx.graph.op(p) == ops.heaviside && cx.graph.children(p).first().is_some_and(|&a| cx.graph.same(a, t))
                })?;
                let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &p)| p).collect();
                Some(mul(cx.graph, &rest))
            };
            if let (Some(f0), Some(g0)) = (strip(cx, f), strip(cx, g)) {
                let body = integrand(cx, f0, g0);
                let integral = call(cx.graph, ops.defint, &[body, u, zero, t]);
                let step = call(cx.graph, ops.heaviside, &[t]);
                mul(cx.graph, &[step, integral])
            } else {
                let body = integrand(cx, f, g);
                let infinity = cx.graph.node(ops.oo, &[]);
                let minus_infinity = neg(cx.graph, infinity);
                call(cx.graph, ops.defint, &[body, u, minus_infinity, infinity])
            }
        },
    };
    let result = cx.simplify(result);
    (!has_heavy(cx.graph, result)).then_some(result)
}

fn eval_q(
    p: &[BigRational],
    x: &BigRational,
) -> BigRational {
    p.iter().rev().fold(BigRational::zero(), |acc, c| acc * x + c)
}

/// `p / (x - r)^m` with `m` maximal.
fn deflate_root(
    p: &[BigRational],
    r: &BigRational,
) -> (Vec<BigRational>, usize) {
    let mut p = p.to_vec();
    let mut m = 0;
    while p.len() > 1 && eval_q(&p, r).is_zero() {
        // synthetic division, highest degree first
        let mut quotient = Vec::with_capacity(p.len() - 1);
        let mut carry = BigRational::zero();
        for c in p.iter().rev().take(p.len() - 1) {
            carry = c + &carry * r;
            quotient.push(carry.clone());
        }
        quotient.reverse();
        p = quotient;
        m += 1;
    }
    (p, m)
}

/// Numeric roots of a polynomial with rational coefficients (Durand–Kerner).
fn complex_roots(p: &[BigRational]) -> Vec<Complex64> {
    let coefficients: Vec<f64> = p.iter().map(|c| Number::rat(c.clone()).to_f64()).collect();
    let n = coefficients.len().saturating_sub(1);
    let Some(&lead) = coefficients.last() else {
        return Vec::new();
    };
    if n == 0 || lead == 0.0 {
        return Vec::new();
    }
    let monic: Vec<f64> = coefficients.iter().map(|c| c / lead).collect();
    let value = |z: Complex64| monic.iter().rev().fold(Complex64::new(0.0, 0.0), |acc, &c| acc * z + c);
    let radius = 1.0 + monic.iter().take(n).copied().map(f64::abs).fold(0.0, f64::max);
    let seed = Complex64::new(0.4, 0.9);
    let mut roots: Vec<Complex64> = (0..n).map(|k| seed.powi(i32::try_from(k).unwrap_or(0)) * radius * 0.5).collect();
    for _ in 0..2000 {
        let mut change = 0.0_f64;
        for k in 0..n {
            let mut denominator = Complex64::new(1.0, 0.0);
            for (j, r) in roots.iter().enumerate() {
                if j != k {
                    denominator *= roots[k] - r;
                }
            }
            if denominator.norm() == 0.0 {
                continue;
            }
            let step = value(roots[k]) / denominator;
            roots[k] -= step;
            change = change.max(step.norm());
        }
        if change < 1e-14 * radius {
            break;
        }
    }
    roots
}

/// The initial and final value theorems for a rational transform.
fn value_theorem(
    cx: &mut Cx<'_>,
    _ops: Ops,
    request: Request,
    f: NodeId,
    v: NodeId,
) -> Option<NodeId> {
    cx.graph.symbol_of(v)?;
    let f = best(cx.graph, f)?;
    let mut gens = Gens::default();
    let gv = gens.index(cx.graph, v);
    let r = ratio(cx.graph, &mut gens, f, Limits::default())?;
    if gens.len() != 1 {
        return None;
    }
    let as_q = |p: &Poly| -> Option<QPoly> { p.univariate_in(gv)?.iter().map(Number::to_rational).collect() };
    let (mut numer, mut denom) = (as_q(&r.numer)?, as_q(&r.denom)?);
    let trim = |p: &mut QPoly| {
        while p.last().is_some_and(Zero::is_zero) {
            p.pop();
        }
    };
    trim(&mut numer);
    trim(&mut denom);
    if denom.is_empty() {
        return None;
    }
    let value = if numer.is_empty() {
        BigRational::zero()
    } else {
        let (dn, dd) = (numer.len(), denom.len());
        match request {
            | Request::InitialZ | Request::InitialLaplace => {
                // F (or s F) must be proper
                let shift = usize::from(request == Request::InitialLaplace);
                match (dn + shift).cmp(&dd) {
                    | std::cmp::Ordering::Greater => return None,
                    | std::cmp::Ordering::Equal => numer.last()? / denom.last()?,
                    | std::cmp::Ordering::Less => BigRational::zero(),
                }
            },
            | _ => {
                // Poles of the part that does not cancel against the
                // numerator must lie in the open unit disc (z) or left
                // half-plane (s), except a simple pole at 1 (z) or 0 (s).
                let z = request == Request::FinalZ;
                let point = if z { BigRational::one() } else { BigRational::zero() };
                let (numer_rest, m_n) = deflate_root(&numer, &point);
                let (denom_rest, m_d) = deflate_root(&denom, &point);
                if m_d > m_n + 1 {
                    return None;
                }
                for root in complex_roots(&denom_rest) {
                    let stable = if z { root.norm() < 1.0 - 1e-9 } else { root.re < -1e-9 };
                    if stable {
                        continue;
                    }
                    // cancelled by the numerator?
                    let at: Complex64 = numer_rest.iter().rev().fold(Complex64::new(0.0, 0.0), |acc, c| acc * root + Number::rat(c.clone()).to_f64());
                    if at.norm() > 1e-8 {
                        return None;
                    }
                }
                if m_d == m_n + 1 {
                    eval_q(&numer_rest, &point) / eval_q(&denom_rest, &point)
                } else {
                    BigRational::zero()
                }
            },
        }
    };
    Some(cx.graph.num(Number::rat(value)))
}

/// One transform request: `x` is the variable of the input, `y` of the
/// output.
struct Tx<'c, 'a> {
    cx: &'c mut Cx<'a>,
    ops: Ops,
    /// `Laplace`, `Fourier` or `Z`: the forward transform in use.
    kind: Request,
    x: NodeId,
    xs: SymbolId,
    y: NodeId,
}

/// Highest power of `t` (or derivative order) handled by the theorems.
const MAX_ORDER: i64 = 12;

impl Tx<'_, '_> {
    fn free(
        &self,
        node: NodeId,
    ) -> bool {
        !self.cx.graph.depends_on(self.cx.graph.find(node), self.xs)
    }

    fn int(
        &mut self,
        v: i64,
    ) -> NodeId {
        self.cx.graph.int(v)
    }

    fn rat(
        &mut self,
        v: &BigRational,
    ) -> NodeId {
        self.cx.graph.num(Number::rat(v.clone()))
    }

    fn half(&mut self) -> NodeId {
        self.cx.graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)))
    }

    fn pi(&mut self) -> NodeId {
        let pi = self.ops.pi;
        self.cx.graph.node(pi, &[])
    }

    fn unit(&mut self) -> NodeId {
        let unit = self.ops.unit;
        self.cx.graph.node(unit, &[])
    }

    fn inverse(
        &mut self,
        node: NodeId,
    ) -> NodeId {
        powi(self.cx.graph, node, -1)
    }

    fn div(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        let inverse = self.inverse(b);
        mul(self.cx.graph, &[a, inverse])
    }

    fn sqrt(
        &mut self,
        a: NodeId,
    ) -> NodeId {
        let half = self.half();
        pow(self.cx.graph, a, half)
    }

    /// Splits a product into the factors free of the variable and the rest.
    #[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
    fn factors(
        &mut self,
        f: NodeId,
    ) -> (Vec<NodeId>, Vec<NodeId>) {
        let factors = if self.cx.graph.op(f) == core::MUL { self.cx.graph.children(f).to_vec() } else { vec![f] };
        factors.into_iter().partition(|&g| self.free(g))
    }

    /// `arg = a x + b` with `a`, `b` free of `x`.
    fn linear(
        &mut self,
        arg: NodeId,
    ) -> Option<(NodeId, NodeId)> {
        let mut gens = Gens::default();
        let gx = gens.index(self.cx.graph, self.x);
        let poly = from_term(self.cx.graph, &mut gens, arg, Limits::default())?;
        if poly.degree_in(gx) > 1 {
            return None;
        }
        let parts = poly.coefficients_in(gx);
        let b = to_term(self.cx.graph, &gens, parts.first()?);
        let a = match parts.get(1) {
            | Some(p) => to_term(self.cx.graph, &gens, p),
            | None => self.int(0),
        };
        (self.free(a) && self.free(b)).then_some((a, b))
    }

    /// One quotient of polynomials in the output variable, when that is
    /// what `r` is.
    fn normal(
        &mut self,
        r: NodeId,
    ) -> NodeId {
        let r = self.cx.simplify(r);
        let Some(cancel) = self.cx.graph.ops().lookup("cancel") else {
            return r;
        };
        let request = call(self.cx.graph, cancel, &[r]);
        let reduced = self.cx.simplify(request);
        if self.cx.graph.op(reduced) == cancel { r } else { reduced }
    }

    fn is_zero(
        &self,
        node: NodeId,
    ) -> bool {
        self.cx.graph.number_of(node).is_some_and(Number::is_zero)
    }

    fn nonnegative(
        &mut self,
        node: NodeId,
    ) -> bool {
        let simplified = self.cx.simplify(node);
        self.cx.graph.facts(simplified).has(Facts::NONNEGATIVE)
    }

    // ------------------------------------------------------------------
    // Unknown functions, inert transforms and the theorems
    // ------------------------------------------------------------------

    /// The forward transform request for `h`: `kind(h, x, out)`.
    fn inert(
        &mut self,
        h: NodeId,
        out: NodeId,
    ) -> NodeId {
        let op = match self.kind {
            | Request::Laplace => self.ops.laplace,
            | Request::Fourier => self.ops.fourier,
            | _ => self.ops.ztransform,
        };
        let x = self.x;
        self.cx.graph.node(op, &[h, x, out])
    }

    /// `(op, h, variable, output)` if `node` is a transform request.
    fn as_inert(
        &self,
        node: NodeId,
    ) -> Option<(OpId, NodeId, NodeId, NodeId)> {
        let op = self.cx.graph.op(node);
        if ![self.ops.laplace, self.ops.fourier, self.ops.ztransform].contains(&op) {
            return None;
        }
        let &[h, variable, out] = self.cx.graph.children(node) else {
            return None;
        };
        Some((op, h, variable, out))
    }

    fn mentions_inert(
        &self,
        node: NodeId,
    ) -> bool {
        let mut stack = vec![node];
        while let Some(n) = stack.pop() {
            if self.as_inert(n).is_some() {
                return true;
            }
            stack.extend_from_slice(self.cx.graph.children(n));
        }
        false
    }

    /// `(head, argument)` of `y(argument)` for an unknown function `y`.
    fn unknown_call(
        &self,
        g: NodeId,
    ) -> Option<(NodeId, NodeId)> {
        if self.cx.graph.op(g) != core::APPLY {
            return None;
        }
        let &[head, argument] = self.cx.graph.children(g) else {
            return None;
        };
        (self.cx.graph.symbol_of(head).is_some() && !self.free(argument)).then_some((head, argument))
    }

    /// `d node / d(output variable)`, also through inert transforms.
    fn d_out(
        &mut self,
        node: NodeId,
    ) -> Option<NodeId> {
        let ys = self.cx.graph.symbol_of(self.y)?;
        let node = best(self.cx.graph, node)?;
        if !self.cx.graph.depends_on(self.cx.graph.find(node), ys) {
            return Some(self.int(0));
        }
        if self.cx.graph.symbol_of(node) == Some(ys) {
            return Some(self.int(1));
        }
        let children = self.cx.graph.children(node).to_vec();
        let depends = |tx: &Self, n: NodeId| tx.cx.graph.depends_on(tx.cx.graph.find(n), ys);
        let derivative_of = match self.cx.graph.op(node) {
            | core::ADD => {
                let mut terms = Vec::with_capacity(children.len());
                for &c in &children {
                    terms.push(self.d_out(c)?);
                }
                add(self.cx.graph, &terms)
            },
            | core::MUL => {
                let mut terms = Vec::with_capacity(children.len());
                for k in 0..children.len() {
                    let mut factors = children.clone();
                    factors[k] = self.d_out(children[k])?;
                    terms.push(mul(self.cx.graph, &factors));
                }
                add(self.cx.graph, &terms)
            },
            | core::POW if children.len() == 2 && !depends(self, children[1]) => {
                let (base, exponent) = (children[0], children[1]);
                let one = self.int(1);
                let lower = sub(self.cx.graph, exponent, one);
                let power = pow(self.cx.graph, base, lower);
                let inner = self.d_out(base)?;
                mul(self.cx.graph, &[exponent, power, inner])
            },
            | op if op == self.ops.diff && children.len() == 2 && self.cx.graph.same(children[1], self.y) => {
                // a further derivative of a derivative request
                let y = self.y;
                call(self.cx.graph, self.ops.diff, &[node, y])
            },
            | _ => {
                if self.as_inert(node).is_some() {
                    // The derivative of a transform of an unknown function in
                    // its output variable stays a derivative request.
                    let out = *children.get(2)?;
                    if !self.cx.graph.same(out, self.y) {
                        return None;
                    }
                    let y = self.y;
                    call(self.cx.graph, self.ops.diff, &[node, y])
                } else if self.mentions_inert(node) {
                    return None;
                } else {
                    return derivative(self.cx.graph, node, self.y);
                }
            },
        };
        Some(self.cx.simplify(derivative_of))
    }

    /// `term` with the output variable `from` replaced by `to`. Derivatives
    /// of transforms `diff(F(from), from)` become `at(diff(F(v), v), v, to)`.
    fn subst_out(
        &mut self,
        term: NodeId,
        from: NodeId,
        to: NodeId,
    ) -> NodeId {
        let children = self.cx.graph.children(term).to_vec();
        if children.is_empty() {
            return self.cx.graph.substitute(term, from, to);
        }
        let op = self.cx.graph.op(term);
        if op == self.ops.diff && children.len() == 2 && self.cx.graph.same(children[1], from) {
            let fresh = self.cx.graph.interner_mut().fresh_symbol("v");
            let v = self.cx.graph.symbol_node(fresh);
            let inner = self.cx.graph.substitute(children[0], from, v);
            let derivative = self.cx.graph.node(op, &[inner, v]);
            return self.cx.graph.node(self.ops.at, &[derivative, v, to]);
        }
        let Some(&dependent) = self.cx.graph.symbol_of(from).as_ref() else {
            return term;
        };
        if !self.cx.graph.depends_on(self.cx.graph.find(term), dependent) {
            return term;
        }
        let rebuilt: Vec<NodeId> = children.iter().map(|&c| self.subst_out(c, from, to)).collect();
        self.cx.graph.try_node(op, &rebuilt).unwrap_or(term)
    }

    /// The expansion of a product containing a sum, so that linearity
    /// applies; `None` if no factor is a sum.
    fn distribute(
        &mut self,
        factors: &[NodeId],
    ) -> Option<NodeId> {
        if factors.len() < 2 || !factors.iter().any(|&g| self.cx.graph.op(g) == core::ADD) {
            return None;
        }
        let expand = self.cx.graph.ops().lookup("expand")?;
        let product = mul(self.cx.graph, factors);
        let request = call(self.cx.graph, expand, &[product]);
        let expanded = self.cx.simplify(request);
        (self.cx.graph.op(expanded) == core::ADD).then_some(expanded)
    }

    /// The transform of an unknown function `y(a x + b)` in terms of the
    /// inert transform of `y(x)`.
    fn unknown(
        &mut self,
        g: NodeId,
    ) -> Option<NodeId> {
        let (head, argument) = self.unknown_call(g)?;
        let (a, b) = self.linear(argument)?;
        if self.is_zero(a) {
            return None;
        }
        let (x, out) = (self.x, self.y);
        let plain = self.cx.graph.node(core::APPLY, &[head, x]);
        let base = self.inert(plain, out);
        let a_is_one = self.cx.graph.number_of(a).is_some_and(Number::is_one);
        match self.kind {
            | Request::Laplace => {
                // y(a t), a > 0: Y(s/a)/a
                if !self.is_zero(b) {
                    return None;
                }
                if a_is_one {
                    return Some(base);
                }
                let positive = self.cx.simplify(a);
                if !self.cx.graph.facts(positive).has(Facts::POSITIVE) {
                    return None;
                }
                let scaled = self.div(out, a);
                let moved = self.subst_out(base, out, scaled);
                Some(self.div(moved, a))
            },
            | Request::Fourier => {
                // y(a t + b): exp(I w b/a) Y(w/a)/|a|
                let scaled = self.div(out, a);
                let moved = if a_is_one { base } else { self.subst_out(base, out, scaled) };
                let mut factors = vec![moved];
                if !self.is_zero(b) {
                    let unit = self.unit();
                    let ratio = self.div(b, a);
                    let exponent = mul(self.cx.graph, &[unit, out, ratio]);
                    factors.push(call(self.cx.graph, self.ops.exp, &[exponent]));
                }
                if !a_is_one {
                    let magnitude = call(self.cx.graph, self.ops.abs, &[a]);
                    let inverse = self.inverse(magnitude);
                    factors.push(inverse);
                }
                Some(mul(self.cx.graph, &factors))
            },
            | _ => {
                // y(n + k)
                if !a_is_one {
                    return None;
                }
                let k = self.cx.graph.number_of(b)?.to_i64()?;
                if k == 0 {
                    return Some(base);
                }
                if k < 0 {
                    // delay: z^k Y (causal: y(m) = 0 for m < 0)
                    if k < -MAX_ORDER {
                        return None;
                    }
                    let power = powi(self.cx.graph, out, k);
                    return Some(mul(self.cx.graph, &[power, base]));
                }
                if k > MAX_ORDER {
                    return None;
                }
                // advance: z^k Y - sum_{j<k} y(j) z^(k-j)
                let leading = powi(self.cx.graph, out, k);
                let mut terms = vec![mul(self.cx.graph, &[leading, base])];
                for j in 0..k {
                    let index = self.int(j);
                    let value = self.cx.graph.node(core::APPLY, &[head, index]);
                    let power = powi(self.cx.graph, out, k - j);
                    let term = mul(self.cx.graph, &[value, power]);
                    terms.push(neg(self.cx.graph, term));
                }
                Some(add(self.cx.graph, &terms))
            },
        }
    }

    /// The transform of `h`, or its inert request.
    fn transform_or_inert(
        &mut self,
        h: NodeId,
    ) -> Option<NodeId> {
        let h = best(self.cx.graph, h)?;
        let done = match self.kind {
            | Request::Laplace => self.laplace(h),
            | Request::Fourier => self.fourier(h),
            | _ => self.ztransform(h),
        };
        Some(match done {
            | Some(r) => r,
            | None => {
                let out = self.y;
                self.inert(h, out)
            },
        })
    }

    /// The convolution theorem: the transform of a convolution is the
    /// product of the transforms.
    fn convolution_theorem(
        &mut self,
        g: NodeId,
    ) -> Option<NodeId> {
        let op = self.cx.graph.op(g);
        let ops = self.ops;
        let wanted = match self.kind {
            | Request::Laplace => op == ops.convolution || op == ops.convolve,
            | Request::Fourier => op == ops.convolution,
            | _ => op == ops.dconvolution || op == ops.convolution,
        };
        if !wanted {
            return None;
        }
        let &[a, b, variable] = self.cx.graph.children(g) else {
            return None;
        };
        if !self.cx.graph.same(variable, self.x) {
            return None;
        }
        let ta = self.transform_or_inert(a)?;
        let tb = self.transform_or_inert(b)?;
        Some(mul(self.cx.graph, &[ta, tb]))
    }

    /// The initial value of `inner` at `x = 0`: `y(0)` or `at(.., x, 0)`.
    fn value_at_zero(
        &mut self,
        inner: NodeId,
    ) -> NodeId {
        let (x, zero) = (self.x, self.int(0));
        if self.cx.graph.op(inner) == self.ops.apply {
            return self.cx.graph.substitute(inner, x, zero);
        }
        call(self.cx.graph, self.ops.at, &[inner, x, zero])
    }

    /// The inverse of a product of inert transforms `kind(h_i, t_i, x)`:
    /// the (discrete) convolution of the `h_i`.
    fn inverse_of_inert(
        &mut self,
        factors: &[NodeId],
    ) -> Option<NodeId> {
        let forward = match self.kind {
            | Request::InverseLaplace => self.ops.laplace,
            | Request::InverseFourier => self.ops.fourier,
            | _ => self.ops.ztransform,
        };
        let mut parts = Vec::new();
        for &f in factors {
            let (op, h, variable, out) = self.as_inert(f)?;
            if op != forward || !self.cx.graph.same(out, self.x) {
                return None;
            }
            let t = self.y;
            parts.push(self.cx.graph.substitute(h, variable, t));
        }
        let mut parts = parts.into_iter();
        let mut acc = parts.next()?;
        let op = if self.kind == Request::InverseZ { self.ops.dconvolution } else { self.ops.convolution };
        for next in parts {
            let t = self.y;
            acc = self.cx.graph.node(op, &[acc, next, t]);
        }
        Some(acc)
    }

    // ------------------------------------------------------------------
    // Laplace
    // ------------------------------------------------------------------

    fn laplace(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) == core::ADD {
            let terms = self.cx.graph.children(f).to_vec();
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                out.push(self.laplace(term)?);
            }
            return Some(add(self.cx.graph, &out));
        }
        let (constants, rest) = self.factors(f);
        let transformed = self.laplace_product(&rest, 0)?;
        let mut all = constants;
        all.push(transformed);
        Some(mul(self.cx.graph, &all))
    }

    /// `L[prod factors](y)`.
    #[allow(clippy::too_many_lines)]
    #[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
    fn laplace_product(
        &mut self,
        factors: &[NodeId],
        depth: usize,
    ) -> Option<NodeId> {
        if depth > 8 {
            return None;
        }
        let s = self.y;
        // L[1] = 1/s
        if factors.is_empty() {
            return Some(self.inverse(s));
        }
        if let [g] = factors
            && let Some(r) = self.unknown(*g).or_else(|| self.convolution_theorem(*g)) {
                return Some(r);
            }
        if let Some(sum) = self.distribute(factors) {
            return self.laplace(sum);
        }
        // exp(a t + b) g(t): e^b G(s - a)
        if let Some(k) = factors.iter().position(|&g| self.cx.graph.op(g) == self.ops.exp) {
            let argument = *self.cx.graph.children(factors[k]).first()?;
            if let Some((a, b)) = self.linear(argument) {
                let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect();
                let g = self.laplace_product(&rest, depth + 1)?;
                let shifted = sub(self.cx.graph, s, a);
                let g = self.subst_out(g, s, shifted);
                if self.is_zero(b) {
                    return Some(g);
                }
                let scale = call(self.cx.graph, self.ops.exp, &[b]);
                return Some(mul(self.cx.graph, &[scale, g]));
            }
        }
        // heaviside(t - c) g(t): e^(-c s) L[g(t + c)], c >= 0
        if let Some(k) = factors.iter().position(|&g| self.cx.graph.op(g) == self.ops.heaviside) {
            let argument = *self.cx.graph.children(factors[k]).first()?;
            let (a, b) = self.linear(argument)?;
            let one = self.cx.graph.number_of(a).is_some_and(|n| n.to_f64() == 1.0);
            let c = neg(self.cx.graph, b);
            if !one || !self.nonnegative(c) {
                return None;
            }
            let x = self.x;
            let shifted_x = add(self.cx.graph, &[x, c]);
            let rest: Vec<NodeId> = factors
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != k)
                .map(|(_, &g)| self.cx.graph.substitute(g, x, shifted_x))
                .collect();
            let rest_term = mul(self.cx.graph, &rest);
            let rest_term = self.cx.simplify(rest_term);
            let g = self.laplace(rest_term)?;
            let exponent = mul(self.cx.graph, &[c, s]);
            let exponent = neg(self.cx.graph, exponent);
            let scale = call(self.cx.graph, self.ops.exp, &[exponent]);
            return Some(mul(self.cx.graph, &[scale, g]));
        }
        // dirac(t - c) g(t): g(c) e^(-c s), c >= 0
        if let Some(k) = factors.iter().position(|&g| self.cx.graph.op(g) == self.ops.dirac) {
            let argument = *self.cx.graph.children(factors[k]).first()?;
            let (a, b) = self.linear(argument)?;
            if !self.cx.graph.number_of(a).is_some_and(|n| n.to_f64() == 1.0) {
                return None;
            }
            let c = neg(self.cx.graph, b);
            if !self.nonnegative(c) {
                return None;
            }
            let x = self.x;
            let rest: Vec<NodeId> = factors
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != k)
                .map(|(_, &g)| self.cx.graph.substitute(g, x, c))
                .collect();
            let exponent = mul(self.cx.graph, &[c, s]);
            let exponent = neg(self.cx.graph, exponent);
            let mut all = rest;
            all.push(call(self.cx.graph, self.ops.exp, &[exponent]));
            return Some(mul(self.cx.graph, &all));
        }
        // t^n g(t): (-1)^n d^n/ds^n G(s)
        if factors.len() > 1
            && let Some((k, n)) = factors.iter().enumerate().find_map(|(k, &g)| Some((k, self.power_of_x(g)?)))
                && n.is_integer() && n.to_f64() >= 1.0 && n.to_f64() <= 12.0 {
                    let n = n.to_i64()?;
                    let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect();
                    let mut g = self.laplace_product(&rest, depth + 1)?;
                    for _ in 0..n {
                        g = self.d_out(g)?;
                    }
                    let sign = self.int(if n % 2 == 0 { 1 } else { -1 });
                    return Some(mul(self.cx.graph, &[sign, g]));
                }
        // Products of two trigonometric or hyperbolic factors, and squares:
        // product-to-sum.
        if let Some(sum) = self.product_to_sum(factors) {
            let sum = self.cx.simplify(sum);
            return self.laplace(sum);
        }
        let [g] = factors else {
            return None;
        };
        let g = *g;
        if let Some(n) = self.power_of_x(g) {
            // t^n = Gamma(n + 1) / s^(n + 1), n > -1
            if n.to_f64() <= -1.0 {
                return None;
            }
            let n_plus_one = add_number(&n, 1);
            let numerator = if n.is_integer() {
                let k = u32::try_from(n.to_i64()?).ok()?;
                let factorial: BigInt = (1..=k).map(BigInt::from).product();
                self.cx.graph.num(Number::Int(factorial))
            } else {
                let e = self.cx.graph.num(n_plus_one.clone());
                call(self.cx.graph, self.ops.gamma, &[e])
            };
            let e = self.cx.graph.num(n_plus_one.neg());
            let power = pow(self.cx.graph, s, e);
            return Some(mul(self.cx.graph, &[numerator, power]));
        }
        let op = self.cx.graph.op(g);
        let children = self.cx.graph.children(g).to_vec();
        if [self.ops.sin, self.ops.cos, self.ops.sinh, self.ops.cosh].contains(&op) {
            let (a, b) = self.linear(*children.first()?)?;
            let s2 = powi(self.cx.graph, s, 2);
            let a2 = powi(self.cx.graph, a, 2);
            let hyperbolic = op == self.ops.sinh || op == self.ops.cosh;
            if hyperbolic && !self.is_zero(b) {
                return None;
            }
            let a2_signed = if hyperbolic { neg(self.cx.graph, a2) } else { a2 };
            let denominator = add(self.cx.graph, &[s2, a2_signed]);
            let (sin_b, cos_b) = (call(self.cx.graph, self.ops.sin, &[b]), call(self.cx.graph, self.ops.cos, &[b]));
            let numerator = if op == self.ops.sin {
                // (s sin b + a cos b)
                let first = mul(self.cx.graph, &[s, sin_b]);
                let second = mul(self.cx.graph, &[a, cos_b]);
                add(self.cx.graph, &[first, second])
            } else if op == self.ops.cos {
                // (s cos b - a sin b)
                let first = mul(self.cx.graph, &[s, cos_b]);
                let second = mul(self.cx.graph, &[a, sin_b]);
                sub(self.cx.graph, first, second)
            } else if op == self.ops.sinh {
                a
            } else {
                s
            };
            return Some(self.div(numerator, denominator));
        }
        // L[y'(t)] = s L[y] - y(0) for an unknown function y.
        if op == self.ops.diff {
            let (&inner, &variable) = (children.first()?, children.get(1)?);
            if !self.cx.graph.same(variable, self.x) || children.len() != 2 {
                return None;
            }
            let inner_transform = self.laplace(inner)?;
            let at_zero = self.value_at_zero(inner);
            let scaled = mul(self.cx.graph, &[s, inner_transform]);
            return Some(sub(self.cx.graph, scaled, at_zero));
        }
        // L[∫_0^t g(u) du] = G(s)/s
        if op == self.ops.defint {
            let (&body, &u, &lower, &upper) = (children.first()?, children.get(1)?, children.get(2)?, children.get(3)?);
            if !self.is_zero(lower) || !self.cx.graph.same(upper, self.x) {
                return None;
            }
            let x = self.x;
            let body = self.cx.graph.substitute(body, u, x);
            let g = self.laplace(body)?;
            return Some(self.div(g, s));
        }
        None
    }

    /// `n` if `g` is `x^n` (or `x`), with `n` a literal.
    fn power_of_x(
        &self,
        g: NodeId,
    ) -> Option<Number> {
        if self.cx.graph.same(g, self.x) {
            return Some(Number::from(1));
        }
        let &[base, exponent] = self.cx.graph.children(g) else {
            return None;
        };
        (self.cx.graph.op(g) == core::POW && self.cx.graph.same(base, self.x)).then(|| self.cx.graph.number_of(exponent).cloned()).flatten()
    }

    /// `sin A sin B`, `cos A cos B`, `sin A cos B` and squares of `sin`,
    /// `cos` rewritten as sums.
    fn product_to_sum(
        &mut self,
        factors: &[NodeId],
    ) -> Option<NodeId> {
        let ops = self.ops;
        let trig = |graph: &Graph, g: NodeId| -> Option<(bool, NodeId)> {
            let op = graph.op(g);
            let &[a] = graph.children(g) else {
                return None;
            };
            (op == ops.sin || op == ops.cos).then_some((op == ops.sin, a))
        };
        // Squares first: sin^2 = (1 - cos 2A)/2, cos^2 = (1 + cos 2A)/2.
        for (k, &g) in factors.iter().enumerate() {
            let &[base, exponent] = self.cx.graph.children(g) else {
                continue;
            };
            if self.cx.graph.op(g) != core::POW || self.cx.graph.number_of(exponent).and_then(Number::to_i64) != Some(2) {
                continue;
            }
            let Some((is_sin, a)) = trig(self.cx.graph, base) else {
                continue;
            };
            let two = self.int(2);
            let double = mul(self.cx.graph, &[two, a]);
            let cos = call(self.cx.graph, ops.cos, &[double]);
            let one = self.int(1);
            let signed = if is_sin { neg(self.cx.graph, cos) } else { cos };
            let numerator = add(self.cx.graph, &[one, signed]);
            let half = self.half();
            let mut rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &h)| h).collect();
            rest.push(half);
            rest.push(numerator);
            let product = mul(self.cx.graph, &rest);
            return Some(crate::rules::complex::build::call(self.cx.graph, core::ADD, &[product]));
        }
        let positions: Vec<(usize, bool, NodeId)> =
            factors.iter().enumerate().filter_map(|(k, &g)| trig(self.cx.graph, g).map(|(s, a)| (k, s, a))).collect();
        let [(i, sin_i, a), (j, sin_j, b), ..] = positions.as_slice() else {
            return None;
        };
        let difference = sub(self.cx.graph, *a, *b);
        let sum = add(self.cx.graph, &[*a, *b]);
        let (cos_d, cos_s) = (call(self.cx.graph, ops.cos, &[difference]), call(self.cx.graph, ops.cos, &[sum]));
        let (sin_d, sin_s) = (call(self.cx.graph, ops.sin, &[difference]), call(self.cx.graph, ops.sin, &[sum]));
        let combined = match (sin_i, sin_j) {
            | (true, true) => sub(self.cx.graph, cos_d, cos_s),
            | (false, false) => add(self.cx.graph, &[cos_d, cos_s]),
            // sin A cos B = (sin(A+B) + sin(A-B))/2
            | (true, false) => add(self.cx.graph, &[sin_s, sin_d]),
            | (false, true) => sub(self.cx.graph, sin_s, sin_d),
        };
        let half = self.half();
        let mut rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(k, _)| k != *i && k != *j).map(|(_, &h)| h).collect();
        rest.push(half);
        rest.push(combined);
        let product = mul(self.cx.graph, &rest);
        // Distribute so that linearity applies.
        let expand = self.cx.graph.ops().lookup("expand")?;
        let expanded = call(self.cx.graph, expand, &[product]);
        Some(expanded)
    }

    // ------------------------------------------------------------------
    // Inverse Laplace
    // ------------------------------------------------------------------

    fn inverse_laplace(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) == core::ADD {
            let terms = self.cx.graph.children(f).to_vec();
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                out.push(self.inverse_laplace(term)?);
            }
            return Some(add(self.cx.graph, &out));
        }
        let (constants, rest) = self.factors(f);
        let mut all = constants;
        // A product of inert transforms: their convolution.
        if !rest.is_empty()
            && let Some(r) = self.inverse_of_inert(&rest) {
                all.push(r);
                return Some(mul(self.cx.graph, &all));
            }
        // exp(-c s) G(s): heaviside(t - c) g(t - c), c >= 0
        if let Some(k) = rest.iter().position(|&g| self.cx.graph.op(g) == self.ops.exp) {
            let argument = *self.cx.graph.children(rest[k]).first()?;
            let (a, b) = self.linear(argument)?;
            let c = neg(self.cx.graph, a);
            if !self.is_zero(b) || !self.nonnegative(c) {
                return None;
            }
            let others: Vec<NodeId> = rest.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect();
            let g = mul(self.cx.graph, &others);
            let inverse = self.inverse_laplace(g)?;
            let t = self.y;
            let shifted = sub(self.cx.graph, t, c);
            let inverse = self.cx.graph.substitute(inverse, t, shifted);
            let step = call(self.cx.graph, self.ops.heaviside, &[shifted]);
            all.push(step);
            all.push(inverse);
            return Some(mul(self.cx.graph, &all));
        }
        let g = mul(self.cx.graph, &rest);
        all.push(self.inverse_rational(g)?);
        Some(mul(self.cx.graph, &all))
    }

    /// Inverse Laplace transform of a rational function of `s`.
    fn inverse_rational(
        &mut self,
        g: NodeId,
    ) -> Option<NodeId> {
        let mut gens = Gens::default();
        let gs = gens.index(self.cx.graph, self.x);
        let r = ratio(self.cx.graph, &mut gens, g, Limits::default())?;
        let as_q = |p: &Poly| -> Option<QPoly> { p.univariate_in(gs)?.iter().map(Number::to_rational).collect() };
        if gens.len() == 1 {
            let (numer, denom) = (as_q(&r.numer)?, as_q(&r.denom)?);
            let parts = apart(&numer, &denom)?;
            let mut terms = Vec::new();
            // A polynomial part is a combination of dirac and derivatives;
            // only constants are supported.
            match parts.quotient.as_slice() {
                | [] => {},
                | [c] => {
                    let c = self.rat(c);
                    let t = self.y;
                    let delta = call(self.cx.graph, self.ops.dirac, &[t]);
                    terms.push(mul(self.cx.graph, &[c, delta]));
                },
                | _ => return None,
            }
            for piece in &parts.pieces {
                let numerator: Vec<NodeId> = piece.numerator.iter().map(|c| self.rat(c)).collect();
                let factor: Vec<NodeId> = piece.factor.iter().map(|c| self.rat(c)).collect();
                terms.push(self.inverse_piece(&numerator, &factor, piece.power)?);
            }
            return Some(add(self.cx.graph, &terms));
        }
        // Symbolic coefficients: a single factor `(p1 s + p0)^m` or
        // `a s^2 + b s + c` in the denominator, as written.
        self.inverse_single_factor(g)
    }

    /// `(n0 + n1 s) / factor^power` for a linear or quadratic factor with
    /// coefficients `factor` (ascending), all free of `s`.
    fn inverse_piece(
        &mut self,
        numerator: &[NodeId],
        factor: &[NodeId],
        power: u32,
    ) -> Option<NodeId> {
        let t = self.y;
        let zero = self.int(0);
        let n0 = numerator.first().copied().unwrap_or(zero);
        let n1 = numerator.get(1).copied().unwrap_or(zero);
        match factor {
            | &[p0, p1] => {
                // (n0 + n1 s)/(p1 s + p0)^m; root r = -p0/p1
                let minus_p0 = neg(self.cx.graph, p0);
                let r = self.div(minus_p0, p1);
                let r = self.cx.simplify(r);
                // n0 + n1 s = n1 (s - r) + (n0 + n1 r)
                let n1r = mul(self.cx.graph, &[n1, r]);
                let constant = add(self.cx.graph, &[n0, n1r]);
                let mut terms = Vec::new();
                let pm = powi(self.cx.graph, p1, i64::from(power));
                for (coefficient, m) in [(constant, power), (n1, power - 1)] {
                    if self.is_zero(coefficient) {
                        continue;
                    }
                    if m == 0 {
                        // n1 * p1^-power * (s - r)^0 = constant: dirac
                        let delta = call(self.cx.graph, self.ops.dirac, &[t]);
                        let inverse = self.inverse(pm);
                        terms.push(mul(self.cx.graph, &[coefficient, inverse, delta]));
                        continue;
                    }
                    // 1/(s - r)^m -> t^(m-1) e^(r t) / (m-1)!
                    let factorial: BigInt = (1..m).map(BigInt::from).product();
                    let scale = self.cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial)));
                    let tp = powi(self.cx.graph, t, i64::from(m) - 1);
                    let rt = mul(self.cx.graph, &[r, t]);
                    let e = call(self.cx.graph, self.ops.exp, &[rt]);
                    let inverse = self.inverse(pm);
                    terms.push(mul(self.cx.graph, &[coefficient, scale, tp, e, inverse]));
                }
                Some(add(self.cx.graph, &terms))
            },
            | &[c, b, a] => {
                // (n1 s + n0)/(a s^2 + b s + c)^m, s = u + alpha
                let two = self.int(2);
                let two_a = mul(self.cx.graph, &[two, a]);
                let minus_b = neg(self.cx.graph, b);
                let alpha = self.div(minus_b, two_a);
                let alpha = self.cx.simplify(alpha);
                // beta^2 = c/a - alpha^2
                let ca = self.div(c, a);
                let alpha2 = powi(self.cx.graph, alpha, 2);
                let beta2 = sub(self.cx.graph, ca, alpha2);
                let beta2 = self.cx.simplify(beta2);
                let facts = self.cx.graph.facts(beta2);
                let trigonometric = if facts.has(Facts::POSITIVE) {
                    true
                } else if facts.has(Facts::NEGATIVE) {
                    false
                } else {
                    return None;
                };
                let magnitude = if trigonometric { beta2 } else { neg(self.cx.graph, beta2) };
                let beta = self.sqrt(magnitude);
                let beta = self.cx.simplify(beta);
                // numerator in u: n1 u + (n0 + n1 alpha), over a^m
                let n1_alpha = mul(self.cx.graph, &[n1, alpha]);
                let constant = add(self.cx.graph, &[n0, n1_alpha]);
                let am = powi(self.cx.graph, a, i64::from(power));
                let bt = mul(self.cx.graph, &[beta, t]);
                let (sin_op, cos_op) =
                    if trigonometric { (self.ops.sin, self.ops.cos) } else { (self.ops.sinh, self.ops.cosh) };
                let (sin, cos) = (call(self.cx.graph, sin_op, &[bt]), call(self.cx.graph, cos_op, &[bt]));
                let core = match power {
                    | 1 => {
                        // u/(u^2 ± b^2) -> cos, 1/(u^2 ± b^2) -> sin/b
                        let first = mul(self.cx.graph, &[n1, cos]);
                        let over_beta = self.div(sin, beta);
                        let second = mul(self.cx.graph, &[constant, over_beta]);
                        add(self.cx.graph, &[first, second])
                    },
                    | 2 if trigonometric => {
                        // u/(u^2+b^2)^2 -> t sin(bt)/(2b)
                        // 1/(u^2+b^2)^2 -> (sin(bt) - bt cos(bt))/(2b^3)
                        let half = self.half();
                        let t_sin = mul(self.cx.graph, &[t, sin]);
                        let over_beta = self.div(t_sin, beta);
                        let first = mul(self.cx.graph, &[n1, half, over_beta]);
                        let bt_cos = mul(self.cx.graph, &[bt, cos]);
                        let difference = sub(self.cx.graph, sin, bt_cos);
                        let beta3 = powi(self.cx.graph, beta, 3);
                        let over = self.div(difference, beta3);
                        let second = mul(self.cx.graph, &[constant, half, over]);
                        add(self.cx.graph, &[first, second])
                    },
                    | m => {
                        // Higher powers by differentiation with respect to the
                        // parameter b: d/db (u^2 ± b^2)^-k = ∓ 2 k b (u^2 ± b^2)^-(k+1),
                        // so F_(k+1) = ∓ (1/(2 k b)) dF_k/db for both numerators.
                        let fresh = self.cx.graph.interner_mut().fresh_symbol("b");
                        self.cx.graph.assume(fresh, Facts::POSITIVE);
                        let b = self.cx.graph.symbol_node(fresh);
                        let bt_param = mul(self.cx.graph, &[b, t]);
                        let (s_p, c_p) = (call(self.cx.graph, sin_op, &[bt_param]), call(self.cx.graph, cos_op, &[bt_param]));
                        let mut sine_part = self.div(s_p, b);
                        let mut cosine_part = c_p;
                        for k in 1..m {
                            let two_k = self.int(2 * i64::from(k));
                            let denominator = mul(self.cx.graph, &[two_k, b]);
                            let weight = self.int(if trigonometric { -1 } else { 1 });
                            for part in [&mut sine_part, &mut cosine_part] {
                                let d = derivative(self.cx.graph, *part, b)?;
                                let scaled = self.div(d, denominator);
                                let signed = mul(self.cx.graph, &[weight, scaled]);
                                *part = self.cx.simplify(signed);
                            }
                        }
                        let first = mul(self.cx.graph, &[n1, cosine_part]);
                        let second = mul(self.cx.graph, &[constant, sine_part]);
                        let total = add(self.cx.graph, &[first, second]);
                        self.cx.graph.substitute(total, b, beta)
                    },
                };
                let at = mul(self.cx.graph, &[alpha, t]);
                let shift = call(self.cx.graph, self.ops.exp, &[at]);
                let inverse = self.inverse(am);
                Some(mul(self.cx.graph, &[shift, core, inverse]))
            },
            | _ => None,
        }
    }

    /// `N(s) / D(s)^m` with `D` linear or quadratic in `s` and `N` of lower
    /// degree, as written in the term.
    fn inverse_single_factor(
        &mut self,
        g: NodeId,
    ) -> Option<NodeId> {
        let factors = if self.cx.graph.op(g) == core::MUL { self.cx.graph.children(g).to_vec() } else { vec![g] };
        let mut numerator = Vec::new();
        let mut denominator = None;
        for f in factors {
            let children = self.cx.graph.children(f).to_vec();
            let negative_power = self.cx.graph.op(f) == core::POW
                && children.get(1).and_then(|&e| self.cx.graph.number_of(e)).and_then(Number::to_i64).is_some_and(|e| e < 0);
            if negative_power && !self.free(f) {
                if denominator.is_some() {
                    return None;
                }
                let m = self.cx.graph.number_of(*children.get(1)?)?.to_i64()?;
                denominator = Some((*children.first()?, u32::try_from(-m).ok()?));
            } else {
                numerator.push(f);
            }
        }
        let (base, power) = denominator?;
        let coefficients = |tx: &mut Self, term: NodeId| -> Option<Vec<NodeId>> {
            let mut gens = Gens::default();
            let gs = gens.index(tx.cx.graph, tx.x);
            let poly = from_term(tx.cx.graph, &mut gens, term, Limits::default())?;
            let parts: Vec<NodeId> = poly.coefficients_in(gs).iter().map(|c| to_term(tx.cx.graph, &gens, c)).collect();
            parts.iter().all(|&c| tx.free(c)).then_some(parts)
        };
        let numerator_term = mul(self.cx.graph, &numerator);
        let numerator = coefficients(self, numerator_term)?;
        let factor = coefficients(self, base)?;
        if !(2..=3).contains(&factor.len()) || numerator.len() >= factor.len() {
            return None;
        }
        self.inverse_piece(&numerator, &factor, power)
    }

    // ------------------------------------------------------------------
    // Fourier
    // ------------------------------------------------------------------

    #[allow(clippy::too_many_lines)]
    fn fourier(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) == core::ADD {
            let terms = self.cx.graph.children(f).to_vec();
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                out.push(self.fourier(term)?);
            }
            return Some(add(self.cx.graph, &out));
        }
        let w = self.y;
        let (constants, rest) = self.factors(f);
        let mut all = constants;
        let transformed = match rest.as_slice() {
            // A constant: 2 pi c dirac(w)
            | [] => {
                let two = self.int(2);
                let pi = self.pi();
                let delta = call(self.cx.graph, self.ops.dirac, &[w]);
                mul(self.cx.graph, &[two, pi, delta])
            },
            | [g] => self.fourier_single(*g)?,
            | _ => {
                if let Some(sum) = self.distribute(&rest) {
                    return self.fourier(sum).map(|r| {
                        all.push(r);
                        mul(self.cx.graph, &all)
                    });
                }
                self.fourier_product(&rest)?
            },
        };
        all.push(transformed);
        Some(mul(self.cx.graph, &all))
    }

    #[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
    fn fourier_single(
        &mut self,
        g: NodeId,
    ) -> Option<NodeId> {
        let w = self.y;
        let op = self.cx.graph.op(g);
        let children = self.cx.graph.children(g).to_vec();
        if let Some(r) = self.unknown(g).or_else(|| self.convolution_theorem(g)) {
            return Some(r);
        }
        if op == self.ops.diff {
            // F[y'] = I w F[y]
            let (&inner, &variable) = (children.first()?, children.get(1)?);
            if !self.cx.graph.same(variable, self.x) || children.len() != 2 {
                return None;
            }
            let transformed = self.fourier(inner)?;
            let unit = self.unit();
            return Some(mul(self.cx.graph, &[unit, w, transformed]));
        }
        if op == self.ops.exp {
            let argument = *children.first()?;
            // exp(-a t^2), a > 0: sqrt(pi/a) exp(-w^2/(4a))
            let mut gens = Gens::default();
            let gx = gens.index(self.cx.graph, self.x);
            if let Some(poly) = from_term(self.cx.graph, &mut gens, argument, Limits::default()) {
                let parts: Vec<NodeId> = poly.coefficients_in(gx).iter().map(|c| to_term(self.cx.graph, &gens, c)).collect();
                if let [c0, c1, c2] = parts.as_slice() {
                    let a = neg(self.cx.graph, *c2);
                    let a = self.cx.simplify(a);
                    if self.cx.graph.facts(a).has(Facts::POSITIVE) && self.free(*c0) && self.free(*c1) {
                        // exp(c0 + c1 t - a t^2) = exp(c0 + c1^2/(4a)) exp(-a (t - tau)^2),
                        // tau = c1/(2a); the shift gives exp(-I w tau)
                        let pi = self.pi();
                        let ratio = self.div(pi, a);
                        let scale = self.sqrt(ratio);
                        let w2 = powi(self.cx.graph, w, 2);
                        let four = self.int(4);
                        let four_a = mul(self.cx.graph, &[four, a]);
                        let exponent = self.div(w2, four_a);
                        let exponent = neg(self.cx.graph, exponent);
                        let mut factors = vec![scale, call(self.cx.graph, self.ops.exp, &[exponent])];
                        if !self.is_zero(*c1) {
                            let c1_squared = powi(self.cx.graph, *c1, 2);
                            let extra = self.div(c1_squared, four_a);
                            let constant = add(self.cx.graph, &[*c0, extra]);
                            factors.push(call(self.cx.graph, self.ops.exp, &[constant]));
                            let two = self.int(2);
                            let two_a = mul(self.cx.graph, &[two, a]);
                            let tau = self.div(*c1, two_a);
                            let unit = self.unit();
                            let phase = mul(self.cx.graph, &[unit, w, tau]);
                            let phase = neg(self.cx.graph, phase);
                            factors.push(call(self.cx.graph, self.ops.exp, &[phase]));
                        } else if !self.is_zero(*c0) {
                            factors.push(call(self.cx.graph, self.ops.exp, &[*c0]));
                        }
                        return Some(mul(self.cx.graph, &factors));
                    }
                }
            }
            // exp(-a abs(t)), a > 0: 2a/(a^2 + w^2)
            let (k, rest) = self.factors(argument);
            if let [abs] = rest.as_slice()
                && self.cx.graph.op(*abs) == self.ops.abs && self.cx.graph.children(*abs).first().is_some_and(|&c| self.cx.graph.same(c, self.x)) {
                    let k = mul(self.cx.graph, &k);
                    let a = neg(self.cx.graph, k);
                    let a = self.cx.simplify(a);
                    if self.cx.graph.facts(a).has(Facts::POSITIVE) {
                        let two = self.int(2);
                        let numerator = mul(self.cx.graph, &[two, a]);
                        let a2 = powi(self.cx.graph, a, 2);
                        let w2 = powi(self.cx.graph, w, 2);
                        let denominator = add(self.cx.graph, &[a2, w2]);
                        return Some(self.div(numerator, denominator));
                    }
                }
            return None;
        }
        if op == self.ops.dirac {
            // dirac(t - c): exp(-I w c)
            let (a, b) = self.linear(*children.first()?)?;
            if !self.cx.graph.number_of(a).is_some_and(|n| n.to_f64() == 1.0) {
                return None;
            }
            let unit = self.unit();
            let exponent = mul(self.cx.graph, &[unit, w, b]);
            return Some(call(self.cx.graph, self.ops.exp, &[exponent]));
        }
        if op == core::POW {
            // 1/(t^2 + a^2): pi/a exp(-a |w|)
            let (&base, &exponent) = (children.first()?, children.get(1)?);
            if self.cx.graph.number_of(exponent).and_then(Number::to_i64) != Some(-1) {
                return None;
            }
            let mut gens = Gens::default();
            let gx = gens.index(self.cx.graph, self.x);
            let poly = from_term(self.cx.graph, &mut gens, base, Limits::default())?;
            let parts: Vec<NodeId> = poly.coefficients_in(gx).iter().map(|c| to_term(self.cx.graph, &gens, c)).collect();
            let [c0, c1, c2] = parts.as_slice() else {
                return None;
            };
            if !self.is_zero(*c1) || !self.free(*c0) || !self.free(*c2) {
                return None;
            }
            // (c2 t^2 + c0)^-1 = 1/c2 * 1/(t^2 + a^2), a^2 = c0/c2
            let a2 = self.div(*c0, *c2);
            let a2 = self.cx.simplify(a2);
            if !self.cx.graph.facts(a2).has(Facts::POSITIVE) {
                return None;
            }
            let a = self.sqrt(a2);
            let a = self.cx.simplify(a);
            let pi = self.pi();
            let over = self.div(pi, a);
            let abs_w = call(self.cx.graph, self.ops.abs, &[w]);
            let aw = mul(self.cx.graph, &[a, abs_w]);
            let exponent = neg(self.cx.graph, aw);
            let e = call(self.cx.graph, self.ops.exp, &[exponent]);
            let inverse_c2 = self.inverse(*c2);
            return Some(mul(self.cx.graph, &[over, e, inverse_c2]));
        }
        None
    }

    /// Products: modulation by `exp(I a t)`, `cos(a t)`, `sin(a t)`,
    /// multiplication by `t`, and `heaviside(t) exp(-a t)`.
    #[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
    fn fourier_product(
        &mut self,
        factors: &[NodeId],
    ) -> Option<NodeId> {
        let w = self.y;
        let others = |k: usize| -> Vec<NodeId> { factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect() };
        // heaviside(t) exp(-a t), a > 0: 1/(a + I w)
        if let [first, second] = factors {
            for (h, e) in [(*first, *second), (*second, *first)] {
                if self.cx.graph.op(h) == self.ops.heaviside
                    && self.cx.graph.children(h).first().is_some_and(|&c| self.cx.graph.same(c, self.x))
                    && self.cx.graph.op(e) == self.ops.exp
                {
                    let argument = *self.cx.graph.children(e).first()?;
                    let (k, b) = self.linear(argument)?;
                    let a = neg(self.cx.graph, k);
                    let a = self.cx.simplify(a);
                    if self.is_zero(b) && self.cx.graph.facts(a).has(Facts::POSITIVE) {
                        let unit = self.unit();
                        let iw = mul(self.cx.graph, &[unit, w]);
                        let denominator = add(self.cx.graph, &[a, iw]);
                        return Some(self.inverse(denominator));
                    }
                }
            }
        }
        for (k, &g) in factors.iter().enumerate() {
            let op = self.cx.graph.op(g);
            let rest = others(k);
            let rest_term = mul(self.cx.graph, &rest);
            if op == self.ops.cos || op == self.ops.sin || op == self.ops.exp {
                let argument = *self.cx.graph.children(g).first()?;
                let Some((a, b)) = self.linear(argument) else {
                    continue;
                };
                if !self.is_zero(b) {
                    continue;
                }
                let inner = self.fourier(rest_term)?;
                let shift = |tx: &mut Self, by: NodeId| {
                    let moved = sub(tx.cx.graph, w, by);
                    tx.subst_out(inner, w, moved)
                };
                if op == self.ops.exp {
                    // exp(I c t) g: G(w - c), with a = I c
                    let unit = self.unit();
                    let c = self.div(a, unit);
                    let c = self.cx.simplify(c);
                    if !self.cx.graph.facts(c).has(Facts::REAL) {
                        continue;
                    }
                    return Some(shift(self, c));
                }
                let minus_a = neg(self.cx.graph, a);
                let (up, down) = (shift(self, a), shift(self, minus_a));
                let half = self.half();
                return Some(if op == self.ops.cos {
                    let sum = add(self.cx.graph, &[up, down]);
                    mul(self.cx.graph, &[half, sum])
                } else {
                    // (G(w - a) - G(w + a)) / (2 I)
                    let difference = sub(self.cx.graph, up, down);
                    let unit = self.unit();
                    let over = self.div(difference, unit);
                    mul(self.cx.graph, &[half, over])
                });
            }
            if let Some(n) = self.power_of_x(g).filter(|n| n.is_integer() && n.to_f64() >= 1.0 && n.to_f64() <= 12.0) {
                // t^n g: (I d/dw)^n G
                let mut inner = self.fourier(rest_term)?;
                for _ in 0..n.to_i64()? {
                    inner = self.d_out(inner)?;
                    let unit = self.unit();
                    inner = mul(self.cx.graph, &[unit, inner]);
                }
                return Some(inner);
            }
        }
        None
    }

    fn inverse_fourier(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let (constants, rest) = self.factors(f);
        let mut all = constants;
        if rest.is_empty() {
            return self.inverse_fourier_by_duality(f);
        }
        if let Some(r) = self.inverse_of_inert(&rest) {
            all.push(r);
            return Some(mul(self.cx.graph, &all));
        }
        if let [g] = rest.as_slice() {
            // convolution in the frequency domain: 2 pi f(t) g(t)
            if self.cx.graph.op(*g) == self.ops.convolution {
                let &[a, b, variable] = self.cx.graph.children(*g) else {
                    return None;
                };
                if !self.cx.graph.same(variable, self.x) {
                    return None;
                }
                let (ia, ib) = (self.inverse_or_inert(a)?, self.inverse_or_inert(b)?);
                let two = self.int(2);
                let pi = self.pi();
                all.extend([two, pi, ia, ib]);
                return Some(mul(self.cx.graph, &all));
            }
        }
        if let Some(r) = self.inverse_fourier_by_duality(f) {
            return Some(r);
        }
        if let Some(r) = self.inverse_fourier_by_residues(f) {
            return Some(r);
        }
        if rest.len() < 2 {
            return None;
        }
        // A product: the convolution of the inverse transforms of the parts.
        let first = rest[0];
        let others = mul(self.cx.graph, &rest[1..]);
        let (a, b) = (self.inverse_fourier(first)?, self.inverse_fourier(others)?);
        let t = self.y;
        let convolution = self.cx.graph.node(self.ops.convolution, &[a, b, t]);
        all.push(convolution);
        Some(mul(self.cx.graph, &all))
    }

    /// The inverse Fourier transform of `h`, or its inert request.
    fn inverse_or_inert(
        &mut self,
        h: NodeId,
    ) -> Option<NodeId> {
        let h = best(self.cx.graph, h)?;
        Some(match self.inverse_fourier(h) {
            | Some(r) => r,
            | None => {
                let (w, t) = (self.x, self.y);
                self.cx.graph.node(self.ops.inverse_fourier, &[h, w, t])
            },
        })
    }

    /// A rational `F(w)` that vanishes at infinity and has no real poles:
    /// `f(t) = (1/2π) ∫ F e^{I w t} dw`, closed in the upper half plane
    /// for `t > 0` and in the lower for `t < 0` (Jordan's lemma).
    fn inverse_fourier_by_residues(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let (w, t) = (self.x, self.y);
        let w_symbol = self.cx.graph.as_symbol(w)?;
        let (poles_op, residue_op) = (
            self.cx.graph.ops().lookup("poles")?,
            self.cx.graph.ops().lookup("residue")?,
        );
        // Only `w` may occur, and F must decay.
        let mut env = HashMap::new();
        if self.cx.graph.free_symbols(self.cx.graph.find(f)) != [w_symbol] {
            return None;
        }
        for big in [1e4, 1e6] {
            env.insert(w_symbol, Complex64::new(big, 0.0));
            let v = self.cx.graph.eval_complex(f, &env)?;
            if !(v.norm() * big).is_finite() || v.norm() * big.sqrt() > 1.0 {
                return None;
            }
        }
        let i_unit = self
            .cx
            .graph
            .ops()
            .lookup("I")
            .map(|op| self.cx.graph.node(op, &[]))?;
        let request = self.cx.graph.node(poles_op, &[f, w]);
        let mut poles = self.cx.simplify(request);
        if self.cx.graph.op(poles) != core::LIST {
            // F = G(I w) with G real (constant-coefficient equations): the
            // poles of G(s), mapped back by w = -I s.
            let s_symbol = self.cx.graph.interner_mut().fresh_symbol("s");
            let s_node = self.cx.graph.symbol_node(s_symbol);
            let minus_i = neg(self.cx.graph, i_unit);
            let w_of_s = mul(self.cx.graph, &[minus_i, s_node]);
            let g = self.cx.graph.substitute(f, w, w_of_s);
            let g = crate::rules::poly::expand_form(self.cx.graph, g).unwrap_or(g);
            let g = self.cx.simplify(g);
            let request = self.cx.graph.node(poles_op, &[g, s_node]);
            let in_s = self.cx.simplify(request);
            if self.cx.graph.op(in_s) != core::LIST {
                return None;
            }
            let mut mapped = Vec::new();
            for pole in self.cx.graph.children(in_s).to_vec() {
                let &[p, m] = self.cx.graph.children(pole) else {
                    return None;
                };
                let p_w = self.cx.graph.substitute(w_of_s, s_node, p);
                let p_w = self.cx.simplify(p_w);
                mapped.push(self.cx.graph.node(core::LIST, &[p_w, m]));
            }
            poles = self.cx.graph.node(core::LIST, &mapped);
        }
        let iwt = mul(self.cx.graph, &[i_unit, w, t]);
        let kernel = call(self.cx.graph, self.ops.exp, &[iwt]);
        let integrand = mul(self.cx.graph, &[f, kernel]);
        let (mut upper, mut lower) = (Vec::new(), Vec::new());
        for pole in self.cx.graph.children(poles).to_vec() {
            let &[p, _] = self.cx.graph.children(pole) else {
                return None;
            };
            let z = self.cx.graph.eval_complex(p, &HashMap::new())?;
            let request = self.cx.graph.node(residue_op, &[integrand, w, p]);
            let r = self.cx.simplify(request);
            if self.cx.graph.op(r) == residue_op {
                return None;
            }
            if z.im > 1e-12 {
                upper.push(r);
            } else if z.im < -1e-12 {
                lower.push(r);
            } else {
                return None;
            }
        }
        if upper.is_empty() && lower.is_empty() {
            return None;
        }
        // t > 0: I Σ_upper;  t < 0: -I Σ_lower.
        let i_unit = self
            .cx
            .graph
            .ops()
            .lookup("I")
            .map(|op| self.cx.graph.node(op, &[]))?;
        let sum_upper = add(self.cx.graph, &upper);
        let positive = mul(self.cx.graph, &[i_unit, sum_upper]);
        let positive = self.cx.simplify(positive);
        let positive = crate::rules::poly::expand_form(self.cx.graph, positive).unwrap_or(positive);
        let positive = self.cx.simplify(positive);
        let minus_i = neg(self.cx.graph, i_unit);
        let sum_lower = add(self.cx.graph, &lower);
        let negative = mul(self.cx.graph, &[minus_i, sum_lower]);
        let negative = self.cx.simplify(negative);
        let negative = crate::rules::poly::expand_form(self.cx.graph, negative).unwrap_or(negative);
        let negative = self.cx.simplify(negative);
        // Even or odd in t: write the answer with abs(t).
        let minus_t = neg(self.cx.graph, t);
        let mirrored = self.cx.graph.substitute(positive, t, minus_t);
        let mirrored = self.cx.simplify(mirrored);
        let abs_t = call(self.cx.graph, self.ops.abs, &[t]);
        let difference = sub(self.cx.graph, mirrored, negative);
        let difference = self.cx.simplify(difference);
        if self
            .cx
            .graph
            .number_of(difference)
            .is_some_and(Number::is_zero)
        {
            return Some(self.cx.graph.substitute(positive, t, abs_t));
        }
        let step = call(self.cx.graph, self.ops.heaviside, &[t]);
        let back_step = call(self.cx.graph, self.ops.heaviside, &[minus_t]);
        let a = mul(self.cx.graph, &[step, positive]);
        let b = mul(self.cx.graph, &[back_step, negative]);
        Some(add(self.cx.graph, &[a, b]))
    }

    fn inverse_fourier_by_duality(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        // (1/2pi) FT[F](-t): duality.
        let (w, t) = (self.x, self.y);
        let fresh = self.cx.graph.interner_mut().fresh_symbol("u");
        let u = self.cx.graph.symbol_node(fresh);
        let mut forward = Tx { cx: &mut *self.cx, ops: self.ops, kind: Request::Fourier, x: w, xs: self.xs, y: u };
        let transformed = forward.fourier(f)?;
        let minus_t = neg(self.cx.graph, t);
        let at = self.cx.graph.substitute(transformed, u, minus_t);
        let two = self.int(2);
        let pi = self.pi();
        let two_pi = mul(self.cx.graph, &[two, pi]);
        Some(self.div(at, two_pi))
    }

    // ------------------------------------------------------------------
    // Z transform
    // ------------------------------------------------------------------

    fn ztransform(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) == core::ADD {
            let terms = self.cx.graph.children(f).to_vec();
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                out.push(self.ztransform(term)?);
            }
            return Some(add(self.cx.graph, &out));
        }
        let (constants, rest) = self.factors(f);
        let transformed = self.z_product(&rest, 0)?;
        let mut all = constants;
        all.push(transformed);
        Some(mul(self.cx.graph, &all))
    }

    #[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
    fn z_product(
        &mut self,
        factors: &[NodeId],
        depth: usize,
    ) -> Option<NodeId> {
        if depth > 8 {
            return None;
        }
        let z = self.y;
        let one = self.int(1);
        if let Some(sum) = self.distribute(factors) {
            return self.ztransform(sum);
        }
        if factors.is_empty() {
            // z/(z - 1)
            let denominator = sub(self.cx.graph, z, one);
            return Some(self.div(z, denominator));
        }
        // a^n g(n): G(z/a); exp(c n) = (e^c)^n
        for (k, &g) in factors.iter().enumerate() {
            // base^(a n + b) = base^b (base^a)^n and exp(c n + d) = e^d (e^c)^n
            let base = match (self.cx.graph.op(g), self.cx.graph.children(g).to_vec().as_slice()) {
                | (core::POW, &[base, exponent]) if self.free(base) && !self.free(exponent) => {
                    self.linear(exponent).filter(|&(a, _)| !self.is_zero(a)).map(|(a, b)| (pow(self.cx.graph, base, a), pow(self.cx.graph, base, b)))
                },
                | (op, &[argument]) if op == self.ops.exp && !self.free(argument) => self
                    .linear(argument)
                    .filter(|&(c, _)| !self.is_zero(c))
                    .map(|(c, d)| (call(self.cx.graph, self.ops.exp, &[c]), call(self.cx.graph, self.ops.exp, &[d]))),
                | _ => None,
            };
            if let Some((base, scale)) = base {
                let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &h)| h).collect();
                let g = self.z_product(&rest, depth + 1)?;
                let scaled = self.div(z, base);
                let moved = self.subst_out(g, z, scaled);
                return Some(mul(self.cx.graph, &[scale, moved]));
            }
        }
        // n g(n): -z G'(z)
        if let Some(k) = factors.iter().position(|&g| self.power_of_x(g).is_some_and(|n| n.is_integer() && n.to_f64() >= 1.0)) {
            let n = self.power_of_x(factors[k])?.to_i64()?;
            if n > MAX_ORDER {
                return None;
            }
            let mut rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &h)| h).collect();
            if n > 1 {
                let x = self.x;
                let lower = powi(self.cx.graph, x, n - 1);
                rest.push(lower);
            }
            let g = self.z_product(&rest, depth + 1)?;
            let d = self.d_out(g)?;
            let minus_z = neg(self.cx.graph, z);
            return Some(mul(self.cx.graph, &[minus_z, d]));
        }
        let [g] = factors else {
            return None;
        };
        if let Some(r) = self.unknown(*g).or_else(|| self.convolution_theorem(*g)) {
            return Some(r);
        }
        let op = self.cx.graph.op(*g);
        let children = self.cx.graph.children(*g).to_vec();
        if op == self.ops.kronecker {
            // kronecker(n - k): z^-k
            let (a, b) = self.linear(*children.first()?)?;
            if !self.cx.graph.number_of(a).is_some_and(|n| n.to_f64() == 1.0) {
                return None;
            }
            let k = self.cx.graph.number_of(b)?.to_i64()?;
            return (k <= 0).then(|| powi(self.cx.graph, z, k));
        }
        if op == self.ops.sin || op == self.ops.cos {
            // sin(b n) -> z sin b/(z^2 - 2 z cos b + 1),
            // cos(b n) -> z (z - cos b)/(same)
            let (b, offset) = self.linear(*children.first()?)?;
            if !self.is_zero(offset) {
                return None;
            }
            let (sin_b, cos_b) = (call(self.cx.graph, self.ops.sin, &[b]), call(self.cx.graph, self.ops.cos, &[b]));
            let z2 = powi(self.cx.graph, z, 2);
            let two = self.int(2);
            let middle = mul(self.cx.graph, &[two, z, cos_b]);
            let middle = neg(self.cx.graph, middle);
            let denominator = add(self.cx.graph, &[z2, middle, one]);
            let numerator = if op == self.ops.sin {
                mul(self.cx.graph, &[z, sin_b])
            } else {
                let shifted = sub(self.cx.graph, z, cos_b);
                mul(self.cx.graph, &[z, shifted])
            };
            return Some(self.div(numerator, denominator));
        }
        None
    }

    fn inverse_z(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let (constants, rest) = self.factors(f);
        if !rest.is_empty()
            && let Some(r) = self.inverse_of_inert(&rest) {
                let mut all = constants;
                all.push(r);
                return Some(mul(self.cx.graph, &all));
            }
        // Partial fractions of F(z)/z over Q.
        let z = self.x;
        let n = self.y;
        let mut gens = Gens::default();
        let gz = gens.index(self.cx.graph, z);
        let r = ratio(self.cx.graph, &mut gens, f, Limits::default())?;
        if gens.len() != 1 {
            return None;
        }
        let as_q = |p: &Poly| -> Option<QPoly> { p.univariate_in(gz)?.iter().map(Number::to_rational).collect() };
        let (numer, mut denom) = (as_q(&r.numer)?, as_q(&r.denom)?);
        denom.insert(0, BigRational::zero());
        let parts = apart(&numer, &denom)?;
        let mut terms = Vec::new();
        // A polynomial part in F/z: F has z^(k+1) terms, which are not
        // causal sequences unless k = -1 (handled as a 1/z piece).
        if !parts.quotient.is_empty() {
            return None;
        }
        for piece in &parts.pieces {
            match (piece.factor.as_slice(), piece.power) {
                | ([c0, c1], m) => {
                    // A z/(z - r)^m with A = numerator / c1^m
                    let root = -c0 / c1;
                    let scale = &piece.numerator[0] / c1.pow(i32::try_from(m).ok()?);
                    if root.is_zero() {
                        // A z / z^m = A z^(1-m): kronecker(n - (m - 1))
                        let shift = self.int(i64::from(m) - 1);
                        let shifted = sub(self.cx.graph, n, shift);
                        let delta = call(self.cx.graph, self.ops.kronecker, &[shifted]);
                        let scale = self.rat(&scale);
                        terms.push(mul(self.cx.graph, &[scale, delta]));
                        continue;
                    }
                    // z/(z - r)^m -> binomial(n, m-1) r^(n - m + 1)
                    let r = self.rat(&root);
                    let mut factors = vec![self.rat(&scale)];
                    let mut falling = Vec::new();
                    for j in 0..m - 1 {
                        let j = self.int(i64::from(j));
                        falling.push(sub(self.cx.graph, n, j));
                    }
                    let factorial: BigInt = (1..m).map(BigInt::from).product();
                    factors.push(self.cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial))));
                    factors.extend(falling);
                    let offset = self.int(i64::from(m) - 1);
                    let exponent = sub(self.cx.graph, n, offset);
                    factors.push(pow(self.cx.graph, r, exponent));
                    terms.push(mul(self.cx.graph, &factors));
                },
                | ([q, p, lead], 1) => {
                    // z (B z + C)/(z^2 + p z + q), monic after dividing by lead,
                    // with complex roots rho e^(± I theta).
                    let (p, q) = (p / lead, q / lead);
                    let zero = BigRational::zero();
                    let big_b = piece.numerator.get(1).cloned().unwrap_or_else(|| zero.clone()) / lead;
                    let big_c = piece.numerator.first().cloned().unwrap_or(zero) / lead;
                    let discriminant = &p * &p - BigRational::from_integer(BigInt::from(4)) * &q;
                    if discriminant.is_positive() {
                        // Real irrational roots r = (-p ± sqrt(d))/2:
                        // z (B z + C)/((z - r1)(z - r2)) -> A1 r1^n + A2 r2^n with
                        // A1 = (B r1 + C)/sqrt(d), A2 = -(B r2 + C)/sqrt(d)
                        let d_node = self.rat(&discriminant);
                        let root_d = self.sqrt(d_node);
                        let minus_p = self.rat(&-&p);
                        let half = self.half();
                        let mut sum_terms = Vec::new();
                        for sign in [1_i64, -1] {
                            let signed = self.int(sign);
                            let offset = mul(self.cx.graph, &[signed, root_d]);
                            let numerator = add(self.cx.graph, &[minus_p, offset]);
                            let r = mul(self.cx.graph, &[half, numerator]);
                            let (b_node, c_node) = (self.rat(&big_b), self.rat(&big_c));
                            let b_r = mul(self.cx.graph, &[b_node, r]);
                            let coefficient = add(self.cx.graph, &[b_r, c_node]);
                            let over = self.div(coefficient, root_d);
                            let weight = mul(self.cx.graph, &[signed, over]);
                            let r_n = pow(self.cx.graph, r, n);
                            sum_terms.push(mul(self.cx.graph, &[weight, r_n]));
                        }
                        terms.push(add(self.cx.graph, &sum_terms));
                        continue;
                    }
                    if !q.is_positive() {
                        return None;
                    }
                    let q_node = self.rat(&q);
                    let rho = self.sqrt(q_node);
                    let minus_half_p = self.rat(&(-&p / BigRational::from_integer(BigInt::from(2))));
                    let cos_theta = self.div(minus_half_p, rho);
                    let acos = self.cx.graph.ops().lookup("acos")?;
                    let theta = call(self.cx.graph, acos, &[cos_theta]);
                    let n_theta = mul(self.cx.graph, &[n, theta]);
                    let (cos, sin) = (call(self.cx.graph, self.ops.cos, &[n_theta]), call(self.cx.graph, self.ops.sin, &[n_theta]));
                    // sin coefficient: (C + B rho cos theta)/(rho sin theta)
                    let b_node = self.rat(&big_b);
                    let c_node = self.rat(&big_c);
                    let b_rho_cos = mul(self.cx.graph, &[b_node, minus_half_p]);
                    let numerator = add(self.cx.graph, &[c_node, b_rho_cos]);
                    let sin_theta = call(self.cx.graph, self.ops.sin, &[theta]);
                    let rho_sin = mul(self.cx.graph, &[rho, sin_theta]);
                    let coefficient = self.div(numerator, rho_sin);
                    let first = mul(self.cx.graph, &[b_node, cos]);
                    let second = mul(self.cx.graph, &[coefficient, sin]);
                    let combination = add(self.cx.graph, &[first, second]);
                    let rho_n = pow(self.cx.graph, rho, n);
                    terms.push(mul(self.cx.graph, &[rho_n, combination]));
                },
                | _ => return None,
            }
        }
        Some(add(self.cx.graph, &terms))
    }

    // ------------------------------------------------------------------
    // Numeric checks
    // ------------------------------------------------------------------

    /// Bindings for the free symbols of `nodes` other than `skip`: fixed
    /// moderate values, positive so that assumptions hold.
    #[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
    fn sample_bindings(
        &mut self,
        nodes: &[NodeId],
        skip: &[SymbolId],
    ) -> Env {
        let mut env = Env::numeric(0.0);
        let mut symbols: Vec<SymbolId> = nodes.iter().flat_map(|&n| self.cx.graph.free_symbols(self.cx.graph.find(n)).to_vec()).collect();
        symbols.sort_unstable();
        symbols.dedup();
        for (k, symbol) in symbols.into_iter().filter(|s| !skip.contains(s)).enumerate() {
            #[allow(clippy::cast_precision_loss)]
            env.bind(symbol, 0.45 + 0.17 * k as f64);
        }
        env
    }

    fn agrees(
        a: f64,
        b: f64,
    ) -> bool {
        (a - b).abs() <= 1e-6 * (1.0 + a.abs().max(b.abs()))
    }

    /// `F(s0) = ∫_0^oo f(t) e^(-s0 t) dt` at two points; inconclusive
    /// evaluations accept.
    fn check_laplace(
        &mut self,
        f: NodeId,
        result: NodeId,
    ) -> bool {
        let Some(result) = best(self.cx.graph, result) else {
            return true;
        };
        let ys = self.cx.graph.symbol_of(self.y);
        let mut env = self.sample_bindings(&[f, result], &[self.xs, ys.unwrap_or(self.xs)]);
        let Some(ys) = ys else {
            return true;
        };
        for s0 in [4.5, 7.25] {
            env.bind(ys, s0);
            let Some(want) = self.cx.graph.eval(result, &env) else {
                return true;
            };
            let graph = &*self.cx.graph;
            let mut point = env.clone();
            let failed = std::cell::Cell::new(false);
            let integrand = |t: f64| {
                point.bind(self.xs, t);
                graph.eval(f, &point).map_or_else(
                    || {
                        failed.set(true);
                        f64::NAN
                    },
                    |v| v * (-s0 * t).exp(),
                )
            };
            let integrand = std::cell::RefCell::new(integrand);
            let quadrature = gauss_kronrod_any(|t| (integrand.borrow_mut())(t), 0.0, f64::INFINITY, 1e-10, 2000);
            if failed.get() || !quadrature.value.is_finite() || !quadrature.error.is_finite() || quadrature.error > 1e-6 {
                return true;
            }
            if !Self::agrees(quadrature.value, want) {
                return false;
            }
        }
        true
    }

    /// The forward transform of the result, numerically, matches `F`.
    fn check_inverse_laplace(
        &mut self,
        f: NodeId,
        result: NodeId,
    ) -> bool {
        // Swap roles: the result is a function of y; F a function of x.
        let Some(result) = best(self.cx.graph, result) else {
            return true;
        };
        let Some(ts) = self.cx.graph.symbol_of(self.y) else {
            return true;
        };
        let mut env = self.sample_bindings(&[f, result], &[self.xs, ts]);
        for s0 in [6.5, 9.0] {
            env.bind(self.xs, s0);
            let Some(want) = self.cx.graph.eval(f, &env) else {
                return true;
            };
            let graph = &*self.cx.graph;
            let point = std::cell::RefCell::new(env.clone());
            let failed = std::cell::Cell::new(false);
            let quadrature = gauss_kronrod_any(
                |t| {
                    let mut point = point.borrow_mut();
                    point.bind(ts, t);
                    graph.eval(result, &point).map_or_else(
                        || {
                            failed.set(true);
                            f64::NAN
                        },
                        |v| v * (-s0 * t).exp(),
                    )
                },
                0.0,
                f64::INFINITY,
                1e-10,
                2000,
            );
            if failed.get() || !quadrature.value.is_finite() || !quadrature.error.is_finite() || quadrature.error > 1e-6 {
                return true;
            }
            if !Self::agrees(quadrature.value, want) {
                return false;
            }
        }
        true
    }

    /// `F(w0) = ∫ f(t) (cos w0 t - I sin w0 t) dt` for real `f`.
    fn check_fourier(
        &mut self,
        f: NodeId,
        result: NodeId,
    ) -> bool {
        let Some(result) = best(self.cx.graph, result) else {
            return true;
        };
        let Some(ws) = self.cx.graph.symbol_of(self.y) else {
            return true;
        };
        let env = self.sample_bindings(&[f, result], &[self.xs, ws]);
        for w0 in [0.7, 1.9] {
            let mut bindings: HashMap<SymbolId, Complex64> = env.bindings().iter().map(|&(s, v)| (s, Complex64::new(v, 0.0))).collect();
            bindings.insert(ws, Complex64::new(w0, 0.0));
            let Some(want) = self.cx.graph.eval_complex(result, &bindings) else {
                return true;
            };
            let graph = &*self.cx.graph;
            let mut parts = [0.0; 2];
            for (slot, phase) in parts.iter_mut().zip([0.0, std::f64::consts::FRAC_PI_2]) {
                let point = std::cell::RefCell::new(env.clone());
                let failed = std::cell::Cell::new(false);
                let quadrature = gauss_kronrod_any(
                    |t| {
                        let mut point = point.borrow_mut();
                        point.bind(self.xs, t);
                        graph.eval(f, &point).map_or_else(
                            || {
                                failed.set(true);
                                f64::NAN
                            },
                            |v| v * (w0 * t - phase).cos(),
                        )
                    },
                    f64::NEG_INFINITY,
                    f64::INFINITY,
                    1e-10,
                    2000,
                );
                if failed.get() || !quadrature.value.is_finite() || !quadrature.error.is_finite() || quadrature.error > 1e-6 {
                    return true;
                }
                *slot = quadrature.value;
            }
            // cos(w t - pi/2) = sin(w t); F = C - I S
            let got = Complex64::new(parts[0], -parts[1]);
            if (got - want).norm() > 1e-6 * (1.0 + want.norm()) {
                return false;
            }
        }
        true
    }

    /// Partial sums of `f(n) z0^-n` against `F(z0)`.
    #[allow(clippy::needless_pass_by_ref_mut)] // signature is shared with the other rule-table entries / call sites
    fn z_series_agrees(
        &mut self,
        sequence: NodeId,
        n: SymbolId,
        transform: NodeId,
        z: SymbolId,
        env: &mut Env,
    ) -> Option<bool> {
        for z0 in [3.5, 5.0] {
            env.bind(z, z0);
            let want = self.cx.graph.eval(transform, env)?;
            let mut sum = 0.0;
            let mut converged = false;
            for k in 0..4000_i32 {
                env.bind(n, f64::from(k));
                let term = self.cx.graph.eval(sequence, env)? * z0.powi(-k);
                sum += term;
                if k > 20 && term.abs() < 1e-15 * (1.0 + sum.abs()) {
                    converged = true;
                    break;
                }
            }
            if !converged || !Self::agrees(sum, want) {
                return Some(converged && Self::agrees(sum, want));
            }
        }
        Some(true)
    }

    fn check_z(
        &mut self,
        f: NodeId,
        result: NodeId,
    ) -> bool {
        let Some(result) = best(self.cx.graph, result) else {
            return true;
        };
        let Some(zs) = self.cx.graph.symbol_of(self.y) else {
            return true;
        };
        let mut env = self.sample_bindings(&[f, result], &[self.xs, zs]);
        // Sampled symbols stay small so that every test series converges
        // at |z| >= 3.5.
        self.z_series_agrees(f, self.xs, result, zs, &mut env).unwrap_or(true)
    }

    fn check_inverse_z(
        &mut self,
        f: NodeId,
        result: NodeId,
    ) -> bool {
        let Some(result) = best(self.cx.graph, result) else {
            return true;
        };
        let Some(ns) = self.cx.graph.symbol_of(self.y) else {
            return true;
        };
        let mut env = self.sample_bindings(&[f, result], &[self.xs, ns]);
        self.z_series_agrees(result, ns, f, self.xs, &mut env).unwrap_or(true)
    }
}

fn add_number(
    n: &Number,
    k: i64,
) -> Number {
    n.to_rational()
        .map_or_else(|| Number::from(n.to_f64() + 1.0), |r| Number::rat(r + BigRational::from_integer(BigInt::from(k))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::Engine;
    use crate::graph::Extractor;
    use crate::graph::Saturate;
    use crate::graph::SizeCost;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[transforms()], src)
    }

    /// Runs `src` with the named symbols assumed positive.
    fn run_positive(
        src: &str,
        positive: &[&str],
    ) -> String {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[transforms()]).unwrap_or_else(|e| panic!("{e}"));
        for name in positive {
            let s = g.interner_mut().symbol(name);
            g.assume(s, Facts::POSITIVE);
        }
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        let term = Extractor::new(&g, &[root], &SizeCost).build(&mut g, root).unwrap_or(root);
        g.display(term)
    }

    #[test]
    fn laplace_table_and_theorems() {
        assert_eq!(run("laplace(1, t, s)"), "1/s");
        assert_eq!(run("laplace(t^3, t, s)"), "6/s^4");
        assert_eq!(run("laplace(exp(2*t), t, s)"), "1/(s - 2)");
        assert_eq!(run("laplace(sin(3*t), t, s)"), "3/(s^2 + 9)");
        assert_eq!(run("laplace(cos(3*t), t, s)"), "s/(s^2 + 9)");
        assert_eq!(run("laplace(sinh(2*t), t, s)"), "2/(s^2 - 4)");
        assert_eq!(run("laplace(t*exp(-t), t, s)"), "1/(s + 1)^2");
        assert_eq!(run("laplace(exp(-t)*sin(t), t, s)"), "1/((s + 1)^2 + 1)");
        assert_eq!(run("laplace(t*sin(t), t, s)"), "2*s/(s^2 + 1)^2");
        assert_eq!(run("laplace(heaviside(t - 2), t, s)"), "exp(-2*s)/s");
        assert_eq!(run("laplace(dirac(t - 1), t, s)"), "exp(-s)");
        assert_eq!(run("laplace(t*dirac(t - 2), t, s)"), "2*exp(-2*s)");
        assert_eq!(run("laplace(3 + 2*t, t, s)"), "3/s + 2/s^2");
        assert_eq!(run("laplace(sin(t)^2, t, s)"), "1/2/s - 1/2*s/(s^2 + 4)");
        assert_eq!(run("laplace(t^(1/2), t, s)"), "1/2*pi^(1/2)/s^(3/2)");
    }

    #[test]
    fn laplace_of_derivatives_and_symbolic_parameters() {
        let (text, _) = crate::rules::testing::reduce_with(&[transforms()], "laplace(diff(y(t), t), t, s)", &[]);
        assert_eq!(text, "s*laplace(y(t), t, s) - y(0)");
        assert_eq!(run_positive("laplace(exp(-a*t), t, s)", &["a"]), "1/(a + s)");
        assert_eq!(run_positive("laplace(sin(w*t), t, s)", &["w"]), "w/(s^2 + w^2)");
    }

    #[test]
    fn inverse_laplace() {
        assert_eq!(run("inverse_laplace(1/(s - 2), s, t)"), "exp(2*t)");
        assert_eq!(run("inverse_laplace(1/s^3, s, t)"), "1/2*t^2");
        assert_eq!(run("inverse_laplace(1/(s^2 + 4), s, t)"), "1/2*sin(2*t)");
        assert_eq!(run("inverse_laplace(s/(s^2 + 4), s, t)"), "cos(2*t)");
        assert_eq!(run("inverse_laplace(1/(s*(s + 1)), s, t)"), "1 - exp(-t)");
        assert_eq!(run("inverse_laplace(1/(s^2 + 2*s + 5), s, t)"), "1/2*exp(-t)*sin(2*t)");
        assert_eq!(run("inverse_laplace(1/(s^2 + 1)^2, s, t)"), "1/2*(sin(t) - t*cos(t))");
        assert_eq!(run("inverse_laplace(exp(-3*s)/s, s, t)"), "heaviside(t - 3)");
        assert_eq!(run_positive("inverse_laplace(1/(s + a), s, t)", &["a"]), "exp(-a*t)");
        assert_eq!(run_positive("inverse_laplace(1/(s^2 + w^2), s, t)", &["w"]), "sin(t*w)/w");
    }

    /// The transform and its inverse compose to the identity.
    fn eval(
        src: &str,
        t: f64,
    ) -> f64 {
        crate::rules::testing::eval(&[transforms()], src, &[("t", t), ("n", t)])
    }

    #[test]
    fn laplace_round_trips() {
        for f in ["t^2*exp(-3*t)", "exp(-t)*cos(2*t)", "sin(t) + t", "cosh(t)", "t*sin(2*t)"] {
            let back = run(&format!("inverse_laplace(laplace({f}, t, s), s, t)"));
            for t in [0.3, 1.7] {
                let (a, b) = (eval(f, t), eval(&back, t));
                assert!((a - b).abs() < 1e-9 * (1.0 + a.abs()), "{f} -> {back}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn fourier_transforms() {
        assert_eq!(run("fourier(exp(-t^2), t, w)"), "exp(-1/4*w^2)*pi^(1/2)");
        assert_eq!(run("fourier(exp(-abs(t)), t, w)"), "2/(w^2 + 1)");
        assert_eq!(run("fourier(heaviside(t)*exp(-2*t), t, w)"), "1/(w*I + 2)");
        assert_eq!(run("fourier(1/(t^2 + 4), t, w)"), "1/2*exp(-2*abs(w))*pi");
        assert_eq!(run("fourier(cos(3*t)*exp(-t^2), t, w)"), "1/2*exp(-1/4*(w - 3)^2)*pi^(1/2) + 1/2*exp(-1/4*(w + 3)^2)*pi^(1/2)");
        assert_eq!(run("fourier(dirac(t - 2), t, w)"), "exp(-2*w*I)");
        assert_eq!(run("inverse_fourier(2/(w^2 + 1), w, t)"), "exp(-abs(t))");
        // Repeated and complex-coefficient denominators by residues.
        assert_eq!(run("inverse_fourier(4/(w^2 + 1)^2, w, t)"), "abs(t)*exp(-abs(t)) + exp(-abs(t))");
        let one_sided = run("inverse_fourier(2/((1 + w^2)*(1 + I*w)), w, t)");
        assert!(one_sided.contains("heaviside(t)") && one_sided.contains("heaviside(-t)") && !one_sided.contains("fourier"), "{one_sided}");
    }

    #[test]
    fn z_transforms() {
        assert_eq!(run("ztransform(1, n, z)"), "z/(z - 1)");
        assert_eq!(run("ztransform(2^n, n, z)"), "z/(z - 2)");
        assert_eq!(run("ztransform(n, n, z)"), "z/(z^2 - 2*z + 1)");
        assert_eq!(run("ztransform(kronecker(n - 2), n, z)"), "1/z^2");
        assert_eq!(run("inverse_ztransform(z/(z - 1/2), z, n)"), "(1/2)^n");
        assert_eq!(run("inverse_ztransform(z/((z - 1)*(z - 2)), z, n)"), "2^n - 1");
        assert_eq!(run("inverse_ztransform(1/z, z, n)"), "kronecker(n - 1)");
        // Complex poles: z/(z^2 + 1) is the transform of sin(n pi/2).
        let back = run("inverse_ztransform(z/(z^2 + 1), z, n)");
        for n in 0..6 {
            let want = (f64::from(n) * std::f64::consts::FRAC_PI_2).sin();
            assert!((eval(&back, f64::from(n)) - want).abs() < 1e-12, "{back}");
        }
        for f in ["3^n*n", "2^n - n", "n^2", "(1/2)^n + 2"] {
            let back = run(&format!("inverse_ztransform(ztransform({f}, n, z), z, n)"));
            for n in [0.0, 1.0, 5.0] {
                let (a, b) = (eval(f, n), eval(&back, n));
                assert!((a - b).abs() < 1e-9 * (1.0 + a.abs()), "{f} -> {back}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn convolution() {
        assert_eq!(run("convolve(1, t, t)"), "1/2*t^2");
        assert_eq!(run("convolve(exp(t), exp(t), t)"), "t*exp(t)");
    }

    /// The text and whether the request reduced to a form without
    /// unresolved requests.
    fn attempt(src: &str) -> (String, bool) {
        crate::rules::testing::reduce_with(&[transforms()], src, &[])
    }

    fn theorem(src: &str) -> String {
        attempt(src).0
    }

    /// The complex value of the closed form of `src`.
    fn value(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> Complex64 {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[transforms()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
        engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        let term = Extractor::new(&g, &[root], &SizeCost).build(&mut g, root).unwrap_or(root);
        let mut map = HashMap::new();
        for (name, v) in bindings {
            let symbol = g.interner_mut().symbol(name);
            map.insert(symbol, Complex64::new(*v, 0.0));
        }
        g.eval_complex(term, &map).unwrap_or_else(|| panic!("cannot evaluate `{}`", g.display(term)))
    }

    #[test]
    fn laplace_derivative_theorem_of_any_order() {
        assert_eq!(theorem("laplace(diff(y(t), t), t, s)"), "s*laplace(y(t), t, s) - y(0)");
        assert_eq!(theorem("laplace(diff(diff(y(t), t), t), t, s)"), "s*(s*laplace(y(t), t, s) - y(0)) - at(diff(y(t), t), t, 0)");
        assert_eq!(
            theorem("laplace(diff(diff(diff(y(t), t), t), t), t, s)"),
            "s*(s*(s*laplace(y(t), t, s) - y(0)) - at(diff(y(t), t), t, 0)) - at(diff(diff(y(t), t), t), t, 0)"
        );
        // a linear ODE y'' + 3 y' + 2 y = 0 transforms to an algebraic equation
        let text = theorem("laplace(diff(diff(y(t), t), t) + 3*diff(y(t), t) + 2*y(t), t, s)");
        assert!(text.contains("laplace(y(t), t, s)") && text.contains("at(diff(y(t), t), t, 0)") && text.contains("y(0)"), "{text}");
    }

    #[test]
    fn laplace_theorems_for_unknown_functions() {
        // frequency shift
        assert_eq!(theorem("laplace(exp(2*t)*y(t), t, s)"), "laplace(y(t), t, s - 2)");
        assert_eq!(theorem("laplace(exp(-a*t)*y(t), t, s)"), "laplace(y(t), t, a + s)");
        // time shift
        assert_eq!(theorem("laplace(heaviside(t - 2)*y(t - 2), t, s)"), "exp(-2*s)*laplace(y(t), t, s)");
        // scaling
        assert_eq!(theorem("laplace(y(3*t), t, s)"), "1/3*laplace(y(t), t, 1/3*s)");
        // integral
        assert_eq!(theorem("laplace(defint(y(u), u, 0, t), t, s)"), "laplace(y(t), t, s)/s");
        // multiplication by t^n
        assert_eq!(theorem("laplace(t*y(t), t, s)"), "-diff(laplace(y(t), t, s), s)");
        assert_eq!(theorem("laplace(t^2*y(t), t, s)"), "diff(diff(laplace(y(t), t, s), s), s)");
        // a shift of t y(t): the derivative moves with the argument
        let text = theorem("laplace(exp(2*t)*t*y(t), t, s)");
        assert!(text.contains("at(diff(laplace(y(t), t, "), "{text}");
        // an unknown transform alone is not reduced
        assert_eq!(attempt("laplace(y(t), t, s)"), ("laplace(y(t), t, s)".to_owned(), false));
        // a product of two unknowns has no theorem
        assert!(!attempt("laplace(y(t)*g(t), t, s)").1);
    }

    #[test]
    fn laplace_theorems_agree_with_closed_forms() {
        // shift: L[e^(2t) sin t] = L[sin](s - 2)
        for s in [3.0, 5.5] {
            let (a, b) = (value("laplace(exp(2*t)*sin(t), t, s)", &[("s", s)]), value("1/((s - 2)^2 + 1)", &[("s", s)]));
            assert!((a - b).norm() < 1e-12);
        }
        // scaling: L[sin(3t)] = L[sin](s/3)/3
        let (a, b) = (value("laplace(sin(3*t), t, s)", &[("s", 4.0)]), value("laplace(sin(t), t, s)/3", &[("s", 4.0 / 3.0)]));
        assert!((a - b).norm() < 1e-12);
        // multiplication by t^2: L[t^2 e^(-t)] = d^2/ds^2 L[e^(-t)]
        let (a, b) = (value("laplace(t^2*exp(-t), t, s)", &[("s", 2.0)]), value("2/(s + 1)^3", &[("s", 2.0)]));
        assert!((a - b).norm() < 1e-12);
        // integral: L[∫_0^t sin] = L[sin]/s
        let (a, b) = (value("laplace(defint(sin(u), u, 0, t), t, s)", &[("s", 2.0)]), value("1/(s*(s^2 + 1))", &[("s", 2.0)]));
        assert!((a - b).norm() < 1e-12);
        // time shift: L[heaviside(t - 2) sin(t - 2)] = e^(-2s) L[sin]
        let (a, b) = (value("laplace(heaviside(t - 2)*sin(t - 2), t, s)", &[("s", 1.5)]), value("exp(-2*s)/(s^2 + 1)", &[("s", 1.5)]));
        assert!((a - b).norm() < 1e-12);
    }

    #[test]
    fn laplace_convolution_theorem() {
        assert_eq!(theorem("laplace(convolution(y(t), g(t), t), t, s)"), "laplace(g(t), t, s)*laplace(y(t), t, s)");
        assert_eq!(theorem("laplace(convolve(y(t), g(t), t), t, s)"), "laplace(g(t), t, s)*laplace(y(t), t, s)");
        // known factors give a closed form: L[e^t] L[sin t]
        assert_eq!(run("laplace(convolve(exp(t), sin(t), t), t, s)"), "1/((s - 1)*(s^2 + 1))");
        assert_eq!(run("laplace(convolution(heaviside(t)*exp(-t), heaviside(t)*exp(-2*t), t), t, s)"), "1/(s^2 + 3*s + 2)");
        // a known factor with an unknown one
        assert_eq!(theorem("laplace(convolution(exp(-t), y(t), t), t, s)"), "laplace(y(t), t, s)/(s + 1)");
        // and the converse
        assert_eq!(theorem("inverse_laplace(laplace(y(t), t, s)*laplace(g(t), t, s), s, t)"), "convolution(y(t), g(t), t)");
        assert_eq!(theorem("inverse_laplace(laplace(y(t), t, s), s, t)"), "y(t)");
        assert_eq!(theorem("inverse_laplace(3*laplace(y(t), t, s)*laplace(g(t), t, s), s, t)"), "3*convolution(y(t), g(t), t)");
        // numerically, for a known pair
        let direct = value("laplace(convolve(t, exp(-t), t), t, s)", &[("s", 2.0)]);
        assert!((direct - value("1/(s^2*(s + 1))", &[("s", 2.0)])).norm() < 1e-12);
    }

    #[test]
    fn fourier_derivative_scaling_and_shifts() {
        assert_eq!(theorem("fourier(diff(y(t), t), t, w)"), "w*I*fourier(y(t), t, w)");
        assert_eq!(theorem("fourier(diff(diff(y(t), t), t), t, w)"), "-w^2*fourier(y(t), t, w)");
        // scaling
        assert_eq!(theorem("fourier(y(2*t), t, w)"), "1/2*fourier(y(t), t, 1/2*w)");
        assert_eq!(theorem("fourier(y(-2*t), t, w)"), "1/2*fourier(y(t), t, -1/2*w)");
        // time shift
        assert_eq!(theorem("fourier(y(t - 3), t, w)"), "exp(-3*w*I)*fourier(y(t), t, w)");
        assert_eq!(theorem("fourier(y(2*t + 4), t, w)"), "1/2*exp(2*w*I)*fourier(y(t), t, 1/2*w)");
        // frequency shift
        assert_eq!(theorem("fourier(exp(I*3*t)*y(t), t, w)"), "fourier(y(t), t, w - 3)");
        // multiplication by t^n
        assert_eq!(theorem("fourier(t*y(t), t, w)"), "I*diff(fourier(y(t), t, w), w)");
        assert_eq!(theorem("fourier(t^2*y(t), t, w)"), "-diff(diff(fourier(y(t), t, w), w), w)");
        // modulation
        let text = theorem("fourier(cos(2*t)*y(t), t, w)");
        assert!(text.contains("fourier(y(t), t, w - 2)") && text.contains("fourier(y(t), t, w + 2)"), "{text}");
    }

    #[test]
    fn fourier_theorems_agree_with_closed_forms() {
        let gaussian = |w: f64| value("fourier(exp(-t^2), t, w)", &[("w", w)]);
        for w in [0.4, 1.3] {
            // time shift by 2
            let shifted = value("fourier(exp(-(t - 2)^2), t, w)", &[("w", w)]);
            assert!((shifted - Complex64::from_polar(1.0, -2.0 * w) * gaussian(w)).norm() < 1e-12);
            // scaling by 2: F(w/2)/2
            let scaled = value("fourier(exp(-(2*t)^2), t, w)", &[("w", w)]);
            assert!((scaled - gaussian(w / 2.0) / 2.0).norm() < 1e-12);
            // frequency shift: exp(I 3 t) g(t) -> G(w - 3)
            let modulated = value("fourier(exp(I*3*t)*exp(-t^2), t, w)", &[("w", w)]);
            assert!((modulated - gaussian(w - 3.0)).norm() < 1e-12);
            // t g(t) = I dG/dw: G = sqrt(pi) exp(-w^2/4), G' = -w/2 G
            let times_t = value("fourier(t*exp(-t^2), t, w)", &[("w", w)]);
            assert!((times_t - Complex64::new(0.0, 1.0) * (-w / 2.0) * gaussian(w)).norm() < 1e-12);
        }
    }

    #[test]
    fn fourier_convolution_theorem() {
        assert_eq!(theorem("fourier(convolution(y(t), g(t), t), t, w)"), "fourier(g(t), t, w)*fourier(y(t), t, w)");
        // two Gaussians: sqrt(pi) e^(-w^2/4) squared
        assert_eq!(run("fourier(convolution(exp(-t^2), exp(-t^2), t), t, w)"), "exp(-1/2*w^2)*pi");
        // known with unknown
        assert_eq!(theorem("fourier(convolution(exp(-t^2), y(t), t), t, w)"), "exp(-1/4*w^2)*fourier(y(t), t, w)*pi^(1/2)");
        // inverse: product of transforms is the convolution
        assert_eq!(theorem("inverse_fourier(fourier(y(t), t, w)*fourier(g(t), t, w), w, t)"), "convolution(y(t), g(t), t)");
        assert_eq!(theorem("inverse_fourier(fourier(y(t), t, w), w, t)"), "y(t)");
        // a convolution in the frequency domain multiplies in time: 2 pi f g
        assert_eq!(run("inverse_fourier(convolution(exp(-w^2), exp(-w^2), w), w, t)"), "1/2*exp(-1/2*t^2)");
        assert_eq!(
            theorem("inverse_fourier(convolution(fourier(y(t), t, w), fourier(g(t), t, w), w), w, t)"),
            "2*g(t)*y(t)*pi"
        );
        // the Gaussian convolution is consistent with its direct evaluation
        assert_eq!(run("convolution(exp(-t^2), exp(-t^2), t)"), "exp(-1/2*t^2)*pi^(1/2)/2^(1/2)");
    }

    #[test]
    fn z_time_shift_and_scaling() {
        // delay of a causal sequence
        assert_eq!(theorem("ztransform(y(n - 3), n, z)"), "ztransform(y(n), n, z)/z^3");
        // advance by two
        assert_eq!(theorem("ztransform(y(n + 2), n, z)"), "z^2*ztransform(y(n), n, z) - z^2*y(0) - z*y(1)");
        assert_eq!(theorem("ztransform(y(n + 1), n, z)"), "z*ztransform(y(n), n, z) - z*y(0)");
        // scaling and multiplication by n
        assert_eq!(theorem("ztransform(2^n*y(n), n, z)"), "ztransform(y(n), n, 1/2*z)");
        assert_eq!(theorem("ztransform(n*y(n), n, z)"), "-z*diff(ztransform(y(n), n, z), z)");
        assert_eq!(
            theorem("ztransform(n^2*y(n), n, z)"),
            "z^2*diff(diff(ztransform(y(n), n, z), z), z) + z*diff(ztransform(y(n), n, z), z)"
        );
        // a difference equation: y(n + 1) - y(n) transforms algebraically
        assert_eq!(theorem("ztransform(y(n + 1) - y(n), n, z)"), "z*ztransform(y(n), n, z) - z*y(0) - ztransform(y(n), n, z)");
    }

    #[test]
    fn z_theorems_agree_with_closed_forms() {
        // advance: Z[x(n + 1)] = z X(z) - z x(0) for x = 2^n: z/(z-2) * z - z
        let advance = value("ztransform(2^(n + 1), n, z)", &[("z", 5.0)]);
        assert!((advance - (5.0 * 5.0 / 3.0 - 5.0)).norm() < 1e-9, "{advance}");
        // multiplication by n: -z d/dz (z/(z - 2)) = 2z/(z - 2)^2
        let by_n = value("ztransform(n*2^n, n, z)", &[("z", 5.0)]);
        assert!((by_n - 10.0 / 9.0).norm() < 1e-9, "{by_n}");
        // scaling: Z[3^n (1/2)^n] = X(z/3) with X = z/(z - 1/2)
        let scaled = value("ztransform(3^n*(1/2)^n, n, z)", &[("z", 9.0)]);
        assert!((scaled - 9.0 / 7.5).norm() < 1e-9, "{scaled}");
    }

    #[test]
    fn z_convolution_theorem() {
        assert_eq!(theorem("ztransform(discrete_convolution(y(n), g(n), n), n, z)"), "ztransform(g(n), n, z)*ztransform(y(n), n, z)");
        assert_eq!(run("ztransform(discrete_convolution(2^n, 3^n, n), n, z)"), "z^2/(z^2 - 5*z + 6)");
        assert_eq!(theorem("inverse_ztransform(ztransform(y(n), n, z)*ztransform(g(n), n, z), z, n)"), "discrete_convolution(y(n), g(n), n)");
        assert_eq!(theorem("inverse_ztransform(ztransform(y(n), n, z), z, n)"), "y(n)");
        // the discrete convolution itself
        assert_eq!(run("discrete_convolution(2^n, 3^n, n)"), "3*3^n - 2^(n + 1)");
        assert_eq!(run("discrete_convolution(1, 1, n)"), "n + 1");
        assert!(!attempt("discrete_convolution(y(n), g(n), n)").1);
        // consistent with the product: inverse z of z^2/((z-2)(z-3))
        let back = run("inverse_ztransform(z^2/((z - 2)*(z - 3)), z, n)");
        for n in [0.0, 1.0, 4.0] {
            assert!((eval(&back, n) - eval("3*3^n - 2^(n + 1)", n)).abs() < 1e-9, "{back}");
        }
    }

    #[test]
    fn convolution_operators() {
        // causal signals
        assert_eq!(run("convolution(heaviside(t), heaviside(t), t)"), "t*heaviside(t)");
        assert_eq!(run("convolution(heaviside(t)*exp(-t), heaviside(t)*exp(-2*t), t)"), "(exp(-t) - exp(-2*t))*heaviside(t)");
        assert_eq!(run("convolution(heaviside(t)*t, heaviside(t), t)"), "1/2*t^2*heaviside(t)");
        // the unilateral convolution of functions on t >= 0
        assert_eq!(run("convolve(t, exp(-t), t)"), "t + exp(-t) - 1");
        // two-sided
        assert_eq!(run("convolution(exp(-t^2), exp(-t^2), t)"), "exp(-1/2*t^2)*pi^(1/2)/2^(1/2)");
        // unknown functions are left alone
        for src in ["convolution(y(t), g(t), t)", "convolve(y(t), g(t), t)"] {
            let (text, reduced) = attempt(src);
            assert_eq!((text.as_str(), reduced), (src, false));
        }
        // commutativity and the unit: delta is not available, heaviside integrates
        assert_eq!(run("convolution(heaviside(t)*exp(-2*t), heaviside(t)*exp(-t), t)"), run("convolution(heaviside(t)*exp(-t), heaviside(t)*exp(-2*t), t)"));
    }

    #[test]
    fn initial_and_final_value_theorems() {
        assert_eq!(run("initial_value_laplace(1/(s + 1), s)"), "1");
        assert_eq!(run("initial_value_laplace(3/(s^2 + 1), s)"), "0");
        assert_eq!(run("initial_value_laplace((2*s + 1)/(s^2 + 3*s + 5), s)"), "2");
        assert_eq!(run("final_value_laplace(1/(s*(s + 1)), s)"), "1");
        assert_eq!(run("final_value_laplace(1/(s + 2), s)"), "0");
        assert_eq!(run("final_value_laplace(3/(s*(s + 2)*(s + 3)), s)"), "1/2");
        assert_eq!(run("initial_value_z(z/(z - 1/2), z)"), "1");
        assert_eq!(run("initial_value_z(1/z, z)"), "0");
        assert_eq!(run("initial_value_z((2*z + 1)/(z^2 + z + 1), z)"), "0");
        assert_eq!(run("initial_value_z((2*z^2 + 1)/(z^2 + z + 1), z)"), "2");
        assert_eq!(run("final_value_z(z/((z - 1)*(z - 1/2)), z)"), "2");
        assert_eq!(run("final_value_z(z/(z - 1/2), z)"), "0");
        // the theorems are checked against the time function
        let x_final = eval(&run("inverse_ztransform(z/((z - 1)*(z - 1/2)), z, n)"), 60.0);
        assert!((x_final - 2.0).abs() < 1e-9);
        let f_final = eval(&run("inverse_laplace(1/(s*(s + 1)), s, t)"), 60.0);
        assert!((f_final - 1.0).abs() < 1e-9);
        // preconditions that fail leave the request alone
        for src in [
            "final_value_laplace(1/(s^2 + 1), s)",
            "final_value_laplace(1/(s - 1), s)",
            "final_value_laplace(1/s^2, s)",
            "final_value_z(z/(z - 2), z)",
            "final_value_z(z/(z + 1), z)",
            "final_value_z(z/(z - 1)^2, z)",
            "initial_value_z(z^2/(z - 1), z)",
            "initial_value_laplace(s^2/(s + 1), s)",
            "initial_value_laplace(a/(s + 1), s)",
        ] {
            let (text, reduced) = attempt(src);
            assert!(!reduced, "{src} gave {text}");
        }
    }

    #[test]
    fn repeated_quadratic_factors_invert() {
        for f in [
            "1/(s^2 + 1)^3",
            "s/(s^2 + 1)^3",
            "(s + 2)/(s^2 + 1)^4",
            "1/(s^2 + 4)^3",
            "1/(s^2 - 1)^2",
            "s/(s^2 - 1)^3",
            "1/((s + 1)^2 + 1)^3",
            "(2*s + 3)/((s + 1)^2 + 4)^3",
            "1/(s*(s^2 + 1)^2)",
        ] {
            let back = run(&format!("inverse_laplace({f}, s, t)"));
            // forward transform of the result, numerically, against the original
            let forward = value(&format!("laplace({back}, t, s)"), &[("s", 3.5)]);
            let want = value(f, &[("s", 3.5)]);
            assert!((forward - want).norm() < 1e-9 * (1.0 + want.norm()), "{f} -> {back}: {forward} vs {want}");
        }
        assert_eq!(run("inverse_laplace(1/(s^2 + 1)^3, s, t)"), "3/8*sin(t) - 1/8*t^2*sin(t) - 3/8*t*cos(t)");
        assert_eq!(run("inverse_laplace(s/(s^2 + 1)^3, s, t)"), "1/8*t*sin(t) - 1/8*t^2*cos(t)");
        assert_eq!(run("inverse_laplace(1/(s^2 - 1)^2, s, t)"), "1/4*t*exp(t) + 1/4*t*exp(-t) - 1/4*exp(t) + 1/4*exp(-t)");
    }

    #[test]
    fn real_irrational_quadratics_invert_in_z() {
        // z/(z^2 - 2): roots ±sqrt 2
        let back = run("inverse_ztransform(z/(z^2 - 2), z, n)");
        // x(0) = 0, x(1) = 1, x(2) = 0, x(3) = 2, x(4) = 0, x(5) = 4
        for (n, want) in [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 2.0), (4.0, 0.0), (5.0, 4.0)] {
            assert!((eval(&back, n) - want).abs() < 1e-9, "n = {n}: {back}");
        }
        // z^2/(z^2 - z - 1): the Fibonacci numbers 1, 1, 2, 3, 5, 8
        let fibonacci = run("inverse_ztransform(z^2/(z^2 - z - 1), z, n)");
        for (n, want) in [(0.0, 1.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 5.0), (5.0, 8.0), (10.0, 89.0)] {
            assert!((eval(&fibonacci, n) - want).abs() < 1e-7, "n = {n}: {fibonacci}");
        }
        // a mix of a rational pole and an irrational pair
        let mixed = run("inverse_ztransform(z/((z - 1)*(z^2 - 3)), z, n)");
        let series = |n: i32| -> f64 {
            // coefficients of 1/((z - 1)(z^2 - 3)) expansion via recurrence
            let mut x = vec![0.0; 40];
            // F(z) = z/(z^3 - z^2 - 3z + 3): x_{n+3} = x_{n+2} + 3 x_{n+1} - 3 x_n with x0 = 0, x1 = 0, x2 = 1
            x[2] = 1.0;
            x[1] = 0.0;
            x[0] = 0.0;
            for k in 3..40 {
                x[k] = x[k - 1] + 3.0 * x[k - 2] - 3.0 * x[k - 3];
            }
            x[usize::try_from(n).unwrap_or(0)]
        };
        for n in [0, 1, 2, 3, 6, 9] {
            assert!((eval(&mixed, f64::from(n)) - series(n)).abs() < 1e-6 * (1.0 + series(n).abs()), "n = {n}: {mixed}");
        }
    }
}
