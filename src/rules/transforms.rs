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
//! | `convolve(f, g, t)` | `∫_0^t f(u) g(t - u) du` |
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
//! Laplace transforms of `diff(y(t), t)` produce `s laplace(y(t), t, s) -
//! y(0)`, so transformed differential equations can be solved for the
//! unknown transform.

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
use crate::rules::poly::univariate::QPoly;
use crate::rules::special::special;

/// The transforms rule set.
#[must_use]
pub fn transforms() -> RuleSet {
    RuleSet::new("transforms", install).needs(calculus()).needs(complex()).needs(special())
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
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let dirac = i.op(OpDescriptor::new("dirac", Arity::Fixed(1)))?;
    let kronecker = i.op(OpDescriptor::new("kronecker", Arity::Fixed(1)).eval(|a| match a {
        | [n] if *n == 0.0 => 1.0,
        | [n] if n.is_finite() => 0.0,
        | _ => f64::NAN,
    }))?;
    i.rewrites(Tier::Normalize, &["transforms/kronecker: kronecker(?n) => 0 if integer(?n), nonzero(?n)"])?;
    let request = |name: &str, binder: bool| {
        let desc = OpDescriptor::new(name, Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100);
        if binder { desc.binder(1, 0b1) } else { desc }
    };
    let table = [
        ("laplace", Request::Laplace),
        ("inverse_laplace", Request::InverseLaplace),
        ("fourier", Request::Fourier),
        ("inverse_fourier", Request::InverseFourier),
        ("ztransform", Request::Z),
        ("inverse_ztransform", Request::InverseZ),
        ("convolve", Request::Convolve),
    ];
    let mut registered = Vec::new();
    for (name, kind) in table {
        registered.push((i.op(request(name, kind != Request::Convolve))?, kind));
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
        let &[f, x, y] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(f) = best(cx.graph, f) else {
            return Outcome::Pass;
        };
        if self.request == Request::Convolve {
            // convolve(f, g, t): the variable is the third argument.
            let Some(ts) = cx.graph.symbol_of(y) else {
                return Outcome::Pass;
            };
            let mut t = Tx { cx, ops: self.ops, x, xs: ts, y };
            return t.convolve(f).map_or(Outcome::Pass, |r| Outcome::Equal(t.cx.simplify(r)));
        }
        let (Some(xs), Some(_)) = (cx.graph.symbol_of(x), cx.graph.symbol_of(y)) else {
            return Outcome::Pass;
        };
        let mut t = Tx { cx, ops: self.ops, x, xs, y };
        let result = match self.request {
            | Request::Laplace => t.laplace(f).filter(|&r| t.check_laplace(f, r)),
            | Request::InverseLaplace => t.inverse_laplace(f).filter(|&r| t.check_inverse_laplace(f, r)),
            | Request::Fourier => t.fourier(f).filter(|&r| t.check_fourier(f, r)),
            | Request::InverseFourier => t.inverse_fourier(f),
            | Request::Z => t.ztransform(f).filter(|&r| t.check_z(f, r)).map(|r| t.normal(r)),
            | Request::InverseZ => t.inverse_z(f).filter(|&r| t.check_inverse_z(f, r)),
            | Request::Convolve => None,
        };
        result.map_or(Outcome::Pass, |r| Outcome::Equal(t.cx.simplify(r)))
    }
}

/// One transform request: `x` is the variable of the input, `y` of the
/// output.
struct Tx<'c, 'a> {
    cx: &'c mut Cx<'a>,
    ops: Ops,
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
        // exp(a t + b) g(t): e^b G(s - a)
        if let Some(k) = factors.iter().position(|&g| self.cx.graph.op(g) == self.ops.exp) {
            let argument = *self.cx.graph.children(factors[k]).first()?;
            if let Some((a, b)) = self.linear(argument) {
                let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect();
                let g = self.laplace_product(&rest, depth + 1)?;
                let shifted = sub(self.cx.graph, s, a);
                let g = self.cx.graph.substitute(g, s, shifted);
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
        if factors.len() > 1 {
            if let Some((k, n)) = factors.iter().enumerate().find_map(|(k, &g)| Some((k, self.power_of_x(g)?))) {
                if n.is_integer() && n.to_f64() >= 1.0 && n.to_f64() <= 12.0 {
                    let n = n.to_i64()?;
                    let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &g)| g).collect();
                    let mut g = self.laplace_product(&rest, depth + 1)?;
                    for _ in 0..n {
                        g = derivative(self.cx.graph, g, s)?;
                        g = self.cx.simplify(g);
                    }
                    let sign = self.int(if n % 2 == 0 { 1 } else { -1 });
                    return Some(mul(self.cx.graph, &[sign, g]));
                }
            }
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
            let x = self.x;
            let inner_transform = if self.cx.graph.op(inner) == self.ops.apply {
                let laplace = self.ops.laplace;
                self.cx.graph.node(laplace, &[inner, x, s])
            } else {
                self.laplace(inner)?
            };
            let zero = self.int(0);
            let at_zero = self.cx.graph.substitute(inner, x, zero);
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
                    | _ => return None,
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
            | _ => self.fourier_product(&rest)?,
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
        if op == self.ops.exp {
            let argument = *children.first()?;
            // exp(-a t^2), a > 0: sqrt(pi/a) exp(-w^2/(4a))
            let mut gens = Gens::default();
            let gx = gens.index(self.cx.graph, self.x);
            if let Some(poly) = from_term(self.cx.graph, &mut gens, argument, Limits::default()) {
                let parts: Vec<NodeId> = poly.coefficients_in(gx).iter().map(|c| to_term(self.cx.graph, &gens, c)).collect();
                if let [c0, c1, c2] = parts.as_slice() {
                    if self.is_zero(*c0) && self.is_zero(*c1) {
                        let a = neg(self.cx.graph, *c2);
                        let a = self.cx.simplify(a);
                        if self.cx.graph.facts(a).has(Facts::POSITIVE) {
                            let pi = self.pi();
                            let ratio = self.div(pi, a);
                            let scale = self.sqrt(ratio);
                            let w2 = powi(self.cx.graph, w, 2);
                            let four = self.int(4);
                            let four_a = mul(self.cx.graph, &[four, a]);
                            let exponent = self.div(w2, four_a);
                            let exponent = neg(self.cx.graph, exponent);
                            let e = call(self.cx.graph, self.ops.exp, &[exponent]);
                            return Some(mul(self.cx.graph, &[scale, e]));
                        }
                    }
                }
            }
            // exp(-a abs(t)), a > 0: 2a/(a^2 + w^2)
            let (k, rest) = self.factors(argument);
            if let [abs] = rest.as_slice() {
                if self.cx.graph.op(*abs) == self.ops.abs && self.cx.graph.children(*abs).first().is_some_and(|&c| self.cx.graph.same(c, self.x)) {
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
                    tx.cx.graph.substitute(inner, w, moved)
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
            if self.power_of_x(g).is_some_and(|n| n.to_f64() == 1.0) {
                // t g: I dG/dw
                let inner = self.fourier(rest_term)?;
                let d = derivative(self.cx.graph, inner, w)?;
                let unit = self.unit();
                return Some(mul(self.cx.graph, &[unit, d]));
            }
        }
        None
    }

    fn inverse_fourier(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        // (1/2pi) FT[F](-t): duality.
        let (w, t) = (self.x, self.y);
        let fresh = self.cx.graph.interner_mut().fresh_symbol("u");
        let u = self.cx.graph.symbol_node(fresh);
        let mut forward = Tx { cx: &mut *self.cx, ops: self.ops, x: w, xs: self.xs, y: u };
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
        if factors.is_empty() {
            // z/(z - 1)
            let denominator = sub(self.cx.graph, z, one);
            return Some(self.div(z, denominator));
        }
        // a^n g(n): G(z/a); exp(c n) = (e^c)^n
        for (k, &g) in factors.iter().enumerate() {
            let base = match (self.cx.graph.op(g), self.cx.graph.children(g)) {
                | (core::POW, &[base, exponent]) if self.free(base) && self.cx.graph.same(exponent, self.x) => Some(base),
                | (op, &[argument]) if op == self.ops.exp => {
                    self.linear(argument).filter(|&(_, b)| self.is_zero(b)).map(|(c, _)| call(self.cx.graph, self.ops.exp, &[c]))
                },
                | _ => None,
            };
            if let Some(base) = base {
                let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, &h)| h).collect();
                let g = self.z_product(&rest, depth + 1)?;
                let scaled = self.div(z, base);
                return Some(self.cx.graph.substitute(g, z, scaled));
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
            let d = derivative(self.cx.graph, g, z)?;
            let minus_z = neg(self.cx.graph, z);
            return Some(mul(self.cx.graph, &[minus_z, d]));
        }
        let [g] = factors else {
            return None;
        };
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
                    if (&p * &p - BigRational::from_integer(BigInt::from(4)) * &q).is_positive() || !q.is_positive() {
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

    fn convolve(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        // convolve(f, g, t): the children are (f, g, t); here x = g, y = t.
        let (g, t) = (self.x, self.y);
        let fresh = self.cx.graph.interner_mut().fresh_symbol("u");
        let u = self.cx.graph.symbol_node(fresh);
        let f_u = self.cx.graph.substitute(f, t, u);
        let g_term = best(self.cx.graph, g)?;
        let t_minus_u = sub(self.cx.graph, t, u);
        let g_shift = self.cx.graph.substitute(g_term, t, t_minus_u);
        let body = mul(self.cx.graph, &[f_u, g_shift]);
        let zero = self.int(0);
        Some(call(self.cx.graph, self.ops.defint, &[body, u, zero, t]))
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
}

