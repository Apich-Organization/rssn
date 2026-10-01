//! Integration kernels.
//!
//! `integral(f, x)` asks for an antiderivative, `defint(f, x, a, b)` for a
//! definite integral. Three kernels answer them:
//!
//! * [`Antiderivative`] is a staged heuristic integrator: linearity, a table
//!   of standard forms with linear arguments, exact partial fractions for
//!   rational functions, powers of sines and cosines, substitution by
//!   "derivative divides", and integration by parts. Whatever it returns
//!   is **checked**: the result is differentiated and compared with the
//!   integrand at sample points, and a disagreement discards it. A wrong
//!   antiderivative therefore cannot enter the graph through a heuristic
//!   that misfired.
//! * [`Definite`] evaluates an antiderivative at finite limits.
//! * [`Quadrature`] answers a definite integral numerically, with an error
//!   estimate, when a numeric result is wanted and no closed form was
//!   found — including over infinite intervals.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::backend::Backend;
use crate::backend::Interpreter;
use crate::graph::op::core;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::kernels::integrate::gauss_kronrod_any;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::apart::apart;
use crate::rules::poly::univariate::QPoly;

use super::diff::Differentiate;

/// Operators the integrator builds results from.
#[derive(Copy, Clone, Debug)]
pub struct Functions {
    pub(crate) diff: OpId,
    pub(crate) exp: OpId,
    pub(crate) ln: OpId,
    pub(crate) sin: OpId,
    pub(crate) cos: OpId,
    pub(crate) tan: OpId,
    pub(crate) asin: OpId,
    pub(crate) acos: OpId,
    pub(crate) atan: OpId,
    pub(crate) sinh: OpId,
    pub(crate) cosh: OpId,
    pub(crate) tanh: OpId,
    pub(crate) sqrt: OpId,
}

const MAX_DEPTH: usize = 6;

/// An antiderivative rule contributed by another rule set: given the
/// integrand and the variable, a candidate antiderivative (which the
/// integrator verifies like its own).
pub type TableFn = fn(&mut Cx<'_>, NodeId, NodeId) -> Option<NodeId>;

/// Operator attribute on `integral`: the antiderivative rules other rule
/// sets have taught the integrator (Gaussians and the error function,
/// for instance, belong to the special functions).
#[derive(Clone, Debug, Default)]
pub struct IntegralTable(pub Vec<TableFn>);

/// One integration problem: the variable and the tools.
struct Integrator<'c, 'a> {
    cx: &'c mut Cx<'a>,
    f: Functions,
    x: NodeId,
    symbol: SymbolId,
    /// Remaining work: every attempt and every nested simplification
    /// spends one unit. The heuristic search is exponential in the worst
    /// case; a failed search must give up in bounded time and leave the
    /// integral to the numeric kernels.
    fuel: usize,
}

/// Work units an antiderivative search may spend.
const FUEL: usize = 600;

impl Integrator<'_, '_> {
    /// Spends one unit of work; `false` once the search is out of fuel.
    const fn spend(&mut self) -> bool {
        match self.fuel.checked_sub(1) {
            | Some(rest) => {
                self.fuel = rest;
                true
            },
            | None => false,
        }
    }

    fn depends(
        &self,
        node: NodeId,
    ) -> bool {
        self.cx.graph.depends_on(self.cx.graph.find(node), self.symbol)
    }

    fn int(
        &mut self,
        n: i64,
    ) -> NodeId {
        self.cx.graph.int(n)
    }

    fn frac(
        &mut self,
        p: i64,
        q: i64,
    ) -> NodeId {
        let value = Number::fraction(p, q).unwrap_or_else(|| Number::from(0));
        self.cx.graph.num(value)
    }

    fn add(
        &mut self,
        terms: &[NodeId],
    ) -> NodeId {
        match terms {
            | [] => self.int(0),
            | [only] => *only,
            | _ => self.cx.graph.node(core::ADD, terms),
        }
    }

    fn mul(
        &mut self,
        factors: &[NodeId],
    ) -> NodeId {
        match factors {
            | [] => self.int(1),
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, factors),
        }
    }

    fn pow(
        &mut self,
        base: NodeId,
        exp: NodeId,
    ) -> NodeId {
        self.cx.graph.node(core::POW, &[base, exp])
    }

    fn inv(
        &mut self,
        node: NodeId,
    ) -> NodeId {
        let minus_one = self.int(-1);
        self.pow(node, minus_one)
    }

    fn div(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        let inverse = self.inv(b);
        self.mul(&[a, inverse])
    }

    fn neg(
        &mut self,
        node: NodeId,
    ) -> NodeId {
        let minus_one = self.int(-1);
        self.mul(&[minus_one, node])
    }

    fn call(
        &mut self,
        op: OpId,
        arg: NodeId,
    ) -> NodeId {
        self.cx.graph.node(op, &[arg])
    }

    fn number(
        &self,
        node: NodeId,
    ) -> Option<Number> {
        self.cx.graph.number_of(node).cloned()
    }

    fn derivative(
        &mut self,
        term: NodeId,
    ) -> NodeId {
        Differentiate { diff: self.f.diff }.derive(self.cx.graph, term, self.symbol, self.x)
    }

    /// If `u = a*x + b` with `a`, `b` free of `x` and `a` non-zero,
    /// returns `a`.
    fn linear(
        &mut self,
        u: NodeId,
    ) -> Option<NodeId> {
        if u == self.x {
            return Some(self.int(1));
        }
        let mut gens = Gens::default();
        let x = self.x;
        let gx = gens.index(self.cx.graph, x);
        let poly = from_term(self.cx.graph, &mut gens, u, Limits::default())?;
        if poly.degree_in(gx) != 1 {
            return None;
        }
        for g in poly.support() {
            if g != gx && gens.node(g).is_some_and(|n| self.depends(n)) {
                return None;
            }
        }
        let slope = poly.coefficients_in(gx).get(1)?.clone();
        Some(to_term(self.cx.graph, &gens, &slope))
    }

    /// An antiderivative of `f` with respect to the integrator's variable.
    fn integrate(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        if depth > MAX_DEPTH || !self.spend() {
            return None;
        }
        if !self.depends(f) {
            let x = self.x;
            return Some(self.mul(&[f, x]));
        }
        let children = self.cx.graph.children(f).to_vec();
        let op = self.cx.graph.op(f);
        if op == core::ADD {
            let mut parts = Vec::with_capacity(children.len());
            for child in children {
                parts.push(self.integrate(child, depth)?);
            }
            return Some(self.add(&parts));
        }
        if op == core::MUL {
            let (constant, dependent): (Vec<NodeId>, Vec<NodeId>) =
                children.iter().partition(|&&c| !self.depends(c));
            if !constant.is_empty() {
                let rest = self.mul(&dependent);
                let inner = self.integrate(rest, depth)?;
                let mut factors = constant;
                factors.push(inner);
                return Some(self.mul(&factors));
            }
        }
        if let Some(found) = self.table(f) {
            return Some(found);
        }
        let extensions = self
            .cx
            .graph
            .ops()
            .lookup("integral")
            .and_then(|op| self.cx.graph.ops().attr::<IntegralTable>(op))
            .map(|t| t.0.clone())
            .unwrap_or_default();
        for rule in extensions {
            if let Some(found) = rule(self.cx, f, self.x) {
                return Some(found);
            }
        }
        if let Some(found) = self.exponential_times_trig(f) {
            return Some(found);
        }
        if let Some(found) = self.rational(f) {
            return Some(found);
        }
        if let Some(found) = self.distribute(f, depth) {
            return Some(found);
        }
        if let Some(found) = self.trig_powers(f, depth) {
            return Some(found);
        }
        if let Some(found) = self.trig_product(f, depth) {
            return Some(found);
        }
        if let Some(found) = self.substitution(f, depth) {
            return Some(found);
        }
        self.by_parts(f, depth)
    }

    /// Standard forms whose argument is linear in the variable.
    fn table(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let x = self.x;
        if f == x {
            let two = self.int(2);
            let half = self.frac(1, 2);
            let square = self.pow(x, two);
            return Some(self.mul(&[half, square]));
        }
        let op = self.cx.graph.op(f);
        let children = self.cx.graph.children(f).to_vec();
        if let (true, &[base, exp]) = (op == core::POW, children.as_slice()) {
            if !self.depends(exp) {
                if let Some(a) = self.linear(base) {
                    // (a*x + b)^n
                    if self.number(exp).is_some_and(|n| n == Number::from(-1)) {
                        let log = self.call(self.f.ln, base);
                        return Some(self.div(log, a));
                    }
                    let one = self.int(1);
                    let next = self.add(&[exp, one]);
                    let raised = self.pow(base, next);
                    let scale = self.mul(&[next, a]);
                    return Some(self.div(raised, scale));
                }
                return self.trig_square(base, exp);
            }
            if !self.depends(base) {
                // c^(a*x + b) = c^u / (a * ln c)
                let a = self.linear(exp)?;
                let log = self.call(self.f.ln, base);
                let scale = self.mul(&[a, log]);
                return Some(self.div(f, scale));
            }
            return None;
        }
        let &[u] = children.as_slice() else {
            return None;
        };
        let a = self.linear(u)?;
        let fs = self.f;
        let antiderivative = if op == fs.exp {
            f
        } else if op == fs.sin {
            let c = self.call(fs.cos, u);
            self.neg(c)
        } else if op == fs.cos {
            self.call(fs.sin, u)
        } else if op == fs.tan {
            let c = self.call(fs.cos, u);
            let log = self.call(fs.ln, c);
            self.neg(log)
        } else if op == fs.sinh {
            self.call(fs.cosh, u)
        } else if op == fs.cosh {
            self.call(fs.sinh, u)
        } else if op == fs.tanh {
            let c = self.call(fs.cosh, u);
            self.call(fs.ln, c)
        } else if op == fs.ln {
            // u*ln(u) - u
            let product = self.mul(&[u, f]);
            let minus_u = self.neg(u);
            self.add(&[product, minus_u])
        } else if op == fs.atan {
            // u*atan(u) - ln(1 + u^2)/2
            let two = self.int(2);
            let one = self.int(1);
            let square = self.pow(u, two);
            let sum = self.add(&[one, square]);
            let log = self.call(fs.ln, sum);
            let minus_half = self.frac(-1, 2);
            let tail = self.mul(&[minus_half, log]);
            let head = self.mul(&[u, f]);
            self.add(&[head, tail])
        } else if op == fs.asin || op == fs.acos {
            // u*asin(u) + sqrt(1 - u^2),  u*acos(u) - sqrt(1 - u^2)
            let two = self.int(2);
            let one = self.int(1);
            let square = self.pow(u, two);
            let minus_square = self.neg(square);
            let radicand = self.add(&[one, minus_square]);
            let half = self.frac(1, 2);
            let root = self.pow(radicand, half);
            let tail = if op == fs.asin { root } else { self.neg(root) };
            let head = self.mul(&[u, f]);
            self.add(&[head, tail])
        } else if op == fs.sqrt {
            // (2/3) u^(3/2)
            let exponent = self.frac(3, 2);
            let raised = self.pow(u, exponent);
            let scale = self.frac(2, 3);
            self.mul(&[scale, raised])
        } else {
            return None;
        };
        Some(self.div(antiderivative, a))
    }

    /// `sin(u)^±2`, `cos(u)^±2`, `tan(u)^2` with linear `u`.
    fn trig_square(
        &mut self,
        base: NodeId,
        exp: NodeId,
    ) -> Option<NodeId> {
        let power = self.number(exp)?.to_i64()?;
        let &[u] = self.cx.graph.children(base) else {
            return None;
        };
        let op = self.cx.graph.op(base);
        let a = self.linear(u)?;
        let fs = self.f;
        let result = match power {
            | 2 if op == fs.sin || op == fs.cos => {
                // u/2 ∓ sin(2u)/4
                let two = self.int(2);
                let double = self.mul(&[two, u]);
                let s = self.call(fs.sin, double);
                let quarter = self.frac(if op == fs.sin { -1 } else { 1 }, 4);
                let tail = self.mul(&[quarter, s]);
                let half = self.frac(1, 2);
                let head = self.mul(&[half, u]);
                self.add(&[head, tail])
            },
            | 2 if op == fs.tan => {
                let minus_u = self.neg(u);
                self.add(&[base, minus_u])
            },
            | -2 if op == fs.cos => self.call(fs.tan, u),
            | -2 if op == fs.sin => {
                // -cos(u)/sin(u)
                let c = self.call(fs.cos, u);
                let quotient = self.div(c, base);
                self.neg(quotient)
            },
            | _ => return None,
        };
        Some(self.div(result, a))
    }

    /// `exp(u) * sin(v)` and `exp(u) * cos(v)` with linear `u`, `v`.
    #[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
    fn exponential_times_trig(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) != core::MUL {
            return None;
        }
        let &[p, q] = self.cx.graph.children(f) else {
            return None;
        };
        let fs = self.f;
        let op_of = |graph: &Graph, n: NodeId| graph.op(n);
        let (e, t) = if op_of(self.cx.graph, p) == fs.exp { (p, q) } else { (q, p) };
        let trig = op_of(self.cx.graph, t);
        if op_of(self.cx.graph, e) != fs.exp || (trig != fs.sin && trig != fs.cos) {
            return None;
        }
        let (&[u], &[v]) = (self.cx.graph.children(e), self.cx.graph.children(t)) else {
            return None;
        };
        let (a, b) = (self.linear(u)?, self.linear(v)?);
        let (s, c) = (self.call(fs.sin, v), self.call(fs.cos, v));
        // e^u (a sin v - b cos v)/(a^2 + b^2),  e^u (a cos v + b sin v)/(a^2 + b^2)
        let bracket = if trig == fs.sin {
            let first = self.mul(&[a, s]);
            let second = self.mul(&[b, c]);
            let second = self.neg(second);
            self.add(&[first, second])
        } else {
            let first = self.mul(&[a, c]);
            let second = self.mul(&[b, s]);
            self.add(&[first, second])
        };
        let two = self.int(2);
        let (a2, b2) = (self.pow(a, two), self.pow(b, two));
        let norm = self.add(&[a2, b2]);
        let numerator = self.mul(&[e, bracket]);
        Some(self.div(numerator, norm))
    }

    /// The term for a polynomial in `x` with rational coefficients.
    fn polynomial(
        &mut self,
        coefficients: &[BigRational],
    ) -> NodeId {
        let mut gens = Gens::default();
        let x = self.x;
        let gx = gens.index(self.cx.graph, x);
        let numbers: Vec<Number> = coefficients.iter().cloned().map(Number::rat).collect();
        to_term(self.cx.graph, &gens, &Poly::from_univariate(gx, &numbers))
    }

    fn rat(
        &mut self,
        value: BigRational,
    ) -> NodeId {
        self.cx.graph.num(Number::rat(value))
    }

    /// Rational functions of `x` over the rationals, by partial fractions.
    fn rational(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let mut gens = Gens::default();
        let x = self.x;
        let gx = gens.index(self.cx.graph, x);
        let fraction = ratio(self.cx.graph, &mut gens, f, Limits::default())?;
        if gens.len() != 1 {
            return None;
        }
        let as_q = |p: &Poly| -> Option<QPoly> {
            let mut q: QPoly = p.univariate_in(gx)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
            while q.last().is_some_and(Zero::is_zero) {
                q.pop();
            }
            Some(q)
        };
        let (numer, denom) = (as_q(&fraction.numer)?, as_q(&fraction.denom)?);
        if denom.len() < 2 {
            // A polynomial, (x^2 - 1)^2 say: integrate its coefficients.
            let scale = denom.first()?.clone();
            if scale.is_zero() || numer.len() < 2 {
                return None;
            }
            let integrated: QPoly = std::iter::once(BigRational::zero())
                .chain(numer.iter().enumerate().map(|(k, c)| c / (&scale * BigRational::from_integer(BigInt::from(k + 1)))))
                .collect();
            return Some(self.polynomial(&integrated));
        }
        let parts = apart(&numer, &denom)?;
        let mut pieces = Vec::new();
        if !parts.quotient.is_empty() {
            let integrated: QPoly = std::iter::once(BigRational::zero())
                .chain(parts.quotient.iter().enumerate().map(|(k, c)| c / BigRational::from_integer(BigInt::from(k + 1))))
                .collect();
            pieces.push(self.polynomial(&integrated));
        }
        for piece in &parts.pieces {
            let base = self.polynomial(&piece.factor);
            pieces.push(self.partial_fraction(&piece.factor, base, piece.power, &piece.numerator)?);
        }
        Some(self.add(&pieces))
    }

    /// `∫ numerator / factor^power dx` for a linear or irreducible
    /// quadratic `factor`.
    fn partial_fraction(
        &mut self,
        factor: &[BigRational],
        base: NodeId,
        power: u32,
        numerator: &[BigRational],
    ) -> Option<NodeId> {
        let fs = self.f;
        match (factor, numerator) {
            | ([_, c1], [a]) => {
                if power == 1 {
                    let log = self.call(fs.ln, base);
                    let scale = self.rat(a / c1);
                    Some(self.mul(&[scale, log]))
                } else {
                    let exponent = 1 - i64::from(power);
                    let e = self.int(exponent);
                    let raised = self.pow(base, e);
                    let scale = self.rat(a / (c1 * BigRational::from_integer(BigInt::from(exponent))));
                    Some(self.mul(&[scale, raised]))
                }
            },
            | ([c, b, a], _) if power == 1 => {
                // (B x + C)/(a x^2 + b x + c)
                let big_c = numerator.first().cloned().unwrap_or_else(BigRational::zero);
                let big_b = numerator.get(1).cloned().unwrap_or_else(BigRational::zero);
                let two = BigRational::from_integer(BigInt::from(2));
                let four = BigRational::from_integer(BigInt::from(4));
                let mut pieces = Vec::new();
                let log_scale = &big_b / (&two * a);
                if !log_scale.is_zero() {
                    let log = self.call(fs.ln, base);
                    let scale = self.rat(log_scale.clone());
                    pieces.push(self.mul(&[scale, log]));
                }
                let rest = &big_c - &log_scale * b;
                if !rest.is_zero() {
                    let discriminant = &four * a * c - b * b;
                    if !discriminant.is_positive() {
                        return None;
                    }
                    let d = self.rat(discriminant);
                    let half = self.frac(1, 2);
                    let root = self.pow(d, half);
                    let linear = self.polynomial(&[b.clone(), &two * a]);
                    let argument = self.div(linear, root);
                    let arc = self.call(fs.atan, argument);
                    let scale = self.rat(&rest * &two);
                    let scaled = self.mul(&[scale, arc]);
                    pieces.push(self.div(scaled, root));
                }
                Some(self.add(&pieces))
            },
            | _ => None,
        }
    }

    /// `sin(u)^m * cos(u)^n` with non-negative integers `m`, `n` and
    /// linear `u`.
    fn trig_powers(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let fs = self.f;
        let factors: Vec<NodeId> =
            if self.cx.graph.op(f) == core::MUL { self.cx.graph.children(f).to_vec() } else { vec![f] };
        let (mut m, mut n, mut argument) = (0_u32, 0_u32, None);
        for factor in factors {
            let (base, power) = match *self.cx.graph.children(factor) {
                | [base, exp] if self.cx.graph.op(factor) == core::POW => {
                    (base, u32::try_from(self.number(exp)?.to_i64()?).ok()?)
                },
                | _ => (factor, 1),
            };
            let op = self.cx.graph.op(base);
            let &[u] = self.cx.graph.children(base) else {
                return None;
            };
            if argument.is_some_and(|a| a != u) {
                return None;
            }
            argument = Some(u);
            if op == fs.sin {
                m += power;
            } else if op == fs.cos {
                n += power;
            } else {
                return None;
            }
        }
        let u = argument?;
        let a = self.linear(u)?;
        if m + n < 2 {
            return None;
        }
        let (s, c) = (self.call(fs.sin, u), self.call(fs.cos, u));
        let fresh = self.cx.graph.interner_mut().fresh_symbol("t");
        let t = self.cx.graph.symbol_node(fresh);
        let (one, two) = (self.int(1), self.int(2));
        if m % 2 == 1 || n % 2 == 1 {
            // An odd power: substitute the other function.
            //   m odd: t = cos u, integrand -(1 - t^2)^((m-1)/2) t^n / a
            //   n odd: t = sin u, integrand  (1 - t^2)^((n-1)/2) t^m / a
            let (odd, other, target, sign) = if m % 2 == 1 { (m, n, c, -1) } else { (n, m, s, 1) };
            let t2 = self.pow(t, two);
            let minus_t2 = self.neg(t2);
            let bracket = self.add(&[one, minus_t2]);
            let half_power = self.int(i64::from((odd - 1) / 2));
            let other_power = self.int(i64::from(other));
            let first = self.pow(bracket, half_power);
            let second = self.pow(t, other_power);
            let integrand = self.mul(&[first, second]);
            let expanded = self.expand(integrand)?;
            let inner = self.in_variable(t, fresh, expanded, depth + 1)?;
            let back = self.cx.graph.substitute(inner, t, target);
            let sign = self.int(sign);
            let scaled = self.mul(&[sign, back]);
            return Some(self.div(scaled, a));
        }
        // Both even: halve the angle with
        //   sin^2 = (1 - cos 2u)/2,  cos^2 = (1 + cos 2u)/2.
        let double = self.mul(&[two, u]);
        let cos_double = self.call(fs.cos, double);
        let half = self.frac(1, 2);
        let minus_cos = self.neg(cos_double);
        let sin_square = {
            let sum = self.add(&[one, minus_cos]);
            self.mul(&[half, sum])
        };
        let cos_square = {
            let sum = self.add(&[one, cos_double]);
            self.mul(&[half, sum])
        };
        let (mp, np) = (self.int(i64::from(m / 2)), self.int(i64::from(n / 2)));
        let first = self.pow(sin_square, mp);
        let second = self.pow(cos_square, np);
        let product = self.mul(&[first, second]);
        let expanded = self.expand(product)?;
        self.integrate(expanded, depth + 1)
    }

    /// Products of two sines/cosines with different linear arguments, by
    /// the product-to-sum formulas:
    /// `sin a sin b = (cos(a-b) - cos(a+b))/2`,
    /// `cos a cos b = (cos(a-b) + cos(a+b))/2`,
    /// `sin a cos b = (sin(a+b) + sin(a-b))/2`.
    /// Any other factors must be free of the variable or polynomial (they
    /// are kept and the result handed back to the integrator).
    fn trig_product(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) != core::MUL {
            return None;
        }
        let fs = self.f;
        let factors = self.cx.graph.children(f).to_vec();
        let trig: Vec<usize> = factors
            .iter()
            .enumerate()
            .filter(|&(_, &n)| {
                let op = self.cx.graph.op(n);
                (op == fs.sin || op == fs.cos) && self.depends(n)
            })
            .map(|(i, _)| i)
            .collect();
        let &[i, j] = trig.as_slice() else {
            return None;
        };
        let (p, q) = (*factors.get(i)?, *factors.get(j)?);
        let (&[a], &[b]) = (self.cx.graph.children(p), self.cx.graph.children(q)) else {
            return None;
        };
        if a == b {
            return None;
        }
        let (op_p, op_q) = (self.cx.graph.op(p), self.cx.graph.op(q));
        let minus_b = self.neg(b);
        let difference = self.add(&[a, minus_b]);
        let sum = self.add(&[a, b]);
        let half = self.frac(1, 2);
        let combined = if op_p == fs.sin && op_q == fs.sin {
            let first = self.call(fs.cos, difference);
            let second = self.call(fs.cos, sum);
            let second = self.neg(second);
            self.add(&[first, second])
        } else if op_p == fs.cos && op_q == fs.cos {
            let first = self.call(fs.cos, difference);
            let second = self.call(fs.cos, sum);
            self.add(&[first, second])
        } else {
            // sin(s) cos(c) with s the sine's argument.
            let (s, c) = if op_p == fs.sin { (a, b) } else { (b, a) };
            let minus_c = self.neg(c);
            let d = self.add(&[s, minus_c]);
            let t = self.add(&[s, c]);
            let first = self.call(fs.sin, t);
            let second = self.call(fs.sin, d);
            self.add(&[first, second])
        };
        let rest: Vec<NodeId> =
            factors.iter().enumerate().filter(|&(k, _)| k != i && k != j).map(|(_, &n)| n).collect();
        let mut all = rest;
        all.push(half);
        all.push(combined);
        let product = self.mul(&all);
        let expanded = self.cx.simplify(product);
        let expanded = self.expand(expanded)?;
        if expanded == f {
            return None;
        }
        self.integrate(expanded, depth + 1)
    }

    /// A product with a sum among its factors, `x (1 - x) sin(a x)`:
    /// multiplied out and integrated term by term.
    fn distribute(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let graph = &*self.cx.graph;
        let has_sum = |n: NodeId| {
            graph.op(n) == core::ADD
                || (graph.op(n) == core::POW
                    && graph.children(n).first().is_some_and(|&b| graph.op(b) == core::ADD)
                    && graph.children(n).get(1).and_then(|&e| graph.number_of(e)).is_some_and(|e| e.is_integer() && !e.is_negative()))
        };
        if graph.op(f) != core::MUL || !graph.children(f).iter().any(|&c| has_sum(c)) {
            return None;
        }
        let mut gens = Gens::default();
        let poly = from_term(self.cx.graph, &mut gens, f, Limits { terms: 64, exponent: 12 })?;
        if poly.len() < 2 {
            return None;
        }
        let expanded = to_term(self.cx.graph, &gens, &poly);
        if self.cx.graph.op(expanded) != core::ADD {
            return None;
        }
        self.integrate(expanded, depth + 1)
    }

    fn expand(
        &mut self,
        term: NodeId,
    ) -> Option<NodeId> {
        let mut gens = Gens::default();
        let poly = from_term(self.cx.graph, &mut gens, term, Limits::default())?;
        Some(to_term(self.cx.graph, &gens, &poly))
    }

    /// Integrates `integrand` with respect to another variable.
    fn in_variable(
        &mut self,
        variable: NodeId,
        symbol: SymbolId,
        integrand: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let mut inner = Integrator { cx: &mut *self.cx, f: self.f, x: variable, symbol, fuel: self.fuel };
        let result = inner.integrate(integrand, depth);
        self.fuel = inner.fuel;
        result
    }

    /// Candidate inner functions for a substitution: arguments of
    /// functions, bases and exponents of powers, and function applications
    /// themselves, largest first.
    fn candidates(
        &self,
        f: NodeId,
    ) -> Vec<NodeId> {
        let graph = &*self.cx.graph;
        let mut out: Vec<NodeId> = Vec::new();
        let mut seen: Vec<NodeId> = Vec::new();
        let mut stack = vec![f];
        while let Some(node) = stack.pop() {
            if seen.contains(&node) {
                continue;
            }
            seen.push(node);
            let op = graph.op(node);
            let children = graph.children(node);
            let arithmetic = op == core::ADD || op == core::MUL;
            for &child in children {
                stack.push(child);
                let interesting = !arithmetic && child != self.x && self.depends(child);
                if interesting && !out.contains(&child) {
                    out.push(child);
                }
            }
            if !arithmetic && op != core::POW && node != f && !children.is_empty() && self.depends(node) && !out.contains(&node)
            {
                out.push(node);
            }
        }
        out.truncate(12);
        out
    }

    /// Substitution `t = u(x)` when `f / u'` is a function of `u` alone.
    fn substitution(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        for u in self.candidates(f) {
            if !self.spend() {
                return None;
            }
            let du = self.derivative(u);
            let du = self.cx.simplify(du);
            if self.number(du).is_some_and(|n| n.is_zero()) || self.cx.graph.op(du) == self.f.diff {
                continue;
            }
            let quotient = self.div(f, du);
            let quotient = self.cx.simplify(quotient);
            let fresh = self.cx.graph.interner_mut().fresh_symbol("t");
            let t = self.cx.graph.symbol_node(fresh);
            let in_t = self.cx.graph.replace_subterm(quotient, u, t);
            if self.depends(in_t) {
                continue;
            }
            if let Some(inner) = self.in_variable(t, fresh, in_t, depth + 1) {
                return Some(self.cx.graph.substitute(inner, t, u));
            }
        }
        None
    }

    /// How early a factor should be chosen as the part to differentiate
    /// (logarithms and inverse trigonometric functions first, then powers
    /// of `x`). `None` for factors that are better integrated.
    fn parts_priority(
        &self,
        factor: NodeId,
    ) -> Option<u8> {
        let graph = &*self.cx.graph;
        let op = graph.op(factor);
        let fs = self.f;
        if [fs.ln, fs.atan, fs.asin, fs.acos].contains(&op) {
            return Some(0);
        }
        if factor == self.x {
            return Some(1);
        }
        match *graph.children(factor) {
            | [base, exp] if op == core::POW && base == self.x => {
                graph.number_of(exp).and_then(Number::to_i64).filter(|&n| n > 0).map(|_| 1)
            },
            | _ => None,
        }
    }

    /// Integration by parts: `∫ u dv = u v - ∫ v du`.
    fn by_parts(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        if self.cx.graph.op(f) != core::MUL {
            return None;
        }
        let factors = self.cx.graph.children(f).to_vec();
        let (index, _) = factors
            .iter()
            .enumerate()
            .filter_map(|(i, &factor)| self.parts_priority(factor).map(|p| (i, p)))
            .min_by_key(|&(_, p)| p)?;
        let u = *factors.get(index)?;
        let rest: Vec<NodeId> = factors.iter().enumerate().filter(|&(i, _)| i != index).map(|(_, &n)| n).collect();
        let dv = self.mul(&rest);
        let v = self.integrate(dv, depth + 1)?;
        let du = self.derivative(u);
        let remaining = self.mul(&[v, du]);
        let remaining = self.cx.simplify(remaining);
        let tail = self.integrate(remaining, depth + 1)?;
        let head = self.mul(&[u, v]);
        let minus_tail = self.neg(tail);
        Some(self.add(&[head, minus_tail]))
    }

    /// Differentiates `candidate` and compares it with `integrand` at
    /// sample points. `false` only on a definite disagreement.
    fn verified(
        &mut self,
        integrand: NodeId,
        candidate: NodeId,
    ) -> bool {
        let derivative = self.derivative(candidate);
        let graph = &*self.cx.graph;
        let mut symbols: Vec<SymbolId> = graph.free_symbols(graph.find(integrand)).to_vec();
        symbols.extend_from_slice(graph.free_symbols(graph.find(derivative)));
        symbols.sort_unstable();
        symbols.dedup();
        for point in [0.37, 0.91, 1.63, 2.3] {
            let mut env = Env::numeric(0.0);
            for (k, &symbol) in symbols.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let value = if symbol == self.symbol { point } else { 0.7 + 0.13 * k as f64 };
                env.bind(symbol, value);
            }
            // Points where either side cannot be evaluated (undetermined
            // functions, a pole, outside the domain) say nothing.
            let (Some(want), Some(got)) = (graph.eval(integrand, &env), graph.eval(derivative, &env)) else {
                continue;
            };
            if want.is_finite() && got.is_finite() && (want - got).abs() > 1e-6 * (1.0 + want.abs().max(got.abs())) {
                return false;
            }
        }
        true
    }
}

/// Finds a checked antiderivative of the best form of `integrand`.
pub(super) fn antiderivative(
    cx: &mut Cx<'_>,
    functions: Functions,
    integrand: NodeId,
    variable: NodeId,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(variable)?;
    let term = cx.simplify(integrand);
    if opaque_in(cx.graph, term, symbol, functions.diff) {
        return None;
    }
    let mut integrator = Integrator { cx, f: functions, x: variable, symbol, fuel: FUEL };
    let candidate = integrator.integrate(term, 0)?;
    integrator.verified(term, candidate).then_some(candidate)
}

/// Whether `term` contains an undetermined function of `x` (`A(x)`) but
/// no derivative of one: such an integrand has no antiderivative the
/// heuristics could find, and searching for one is expensive. (With a
/// derivative present, `∫ f(x) f'(x) dx`, substitution may succeed.)
fn opaque_in(
    graph: &Graph,
    term: NodeId,
    x: SymbolId,
    diff: OpId,
) -> bool {
    let mut stack = vec![term];
    let mut seen = std::collections::HashSet::new();
    let mut found = false;
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        let op = graph.op(n);
        if op == diff {
            return false;
        }
        if op == core::APPLY && graph.depends_on(graph.find(n), x) {
            found = true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    found
}

/// Symbolic kernel for `integral(f, x)`.
pub(super) struct Antiderivative {
    pub(super) integral: OpId,
    pub(super) functions: Functions,
}

impl Kernel for Antiderivative {
    fn ops(&self) -> Vec<OpId> {
        vec![self.integral]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[integrand, variable] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        antiderivative(cx, self.functions, integrand, variable).map_or(Outcome::Pass, Outcome::Equal)
    }
}

/// Symbolic kernel for `defint(f, x, a, b)`: the antiderivative at the
/// limits, or its limits at infinite ends.
pub(super) struct Definite {
    pub(super) defint: OpId,
    pub(super) functions: Functions,
}

impl Kernel for Definite {
    fn ops(&self) -> Vec<OpId> {
        vec![self.defint]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[integrand, variable, lower, upper] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        // An infinite limit needs a limit, not a substitution.
        let infinite = |graph: &Graph, n: NodeId| {
            graph.eval(n, &Env::numeric(0.0)).is_some_and(f64::is_infinite)
        };
        let (lower_infinite, upper_infinite) = (infinite(cx.graph, lower), infinite(cx.graph, upper));
        // Inside the integral the variable ranges over the interval: a
        // fresh symbol carries what is known about it (real, and its sign
        // when the interval lies on one side of zero), so that |x| or
        // sqrt(x^2) can be resolved.
        let (variable, integrand) = match (cx.graph.eval(lower, &Env::numeric(0.0)), cx.graph.eval(upper, &Env::numeric(0.0))) {
            | (Some(lo), Some(hi)) if lo <= hi && cx.graph.symbol_of(variable).is_some() => {
                // The end points are a null set: the open interval is
                // what matters to the value of the integral.
                let mut facts = crate::graph::Facts::REAL;
                if lo >= 0.0 && hi > lo {
                    facts = facts | crate::graph::Facts::POSITIVE;
                }
                if hi <= 0.0 && hi > lo {
                    facts = facts | crate::graph::Facts::NEGATIVE;
                }
                let fresh = cx.graph.interner_mut().fresh_symbol("x");
                cx.graph.assume(fresh, facts);
                let local = cx.graph.symbol_node(fresh);
                let body = best(cx.graph, integrand).unwrap_or(integrand);
                (local, cx.graph.substitute(body, variable, local))
            },
            | _ => (variable, integrand),
        };
        let Some(primitive) = antiderivative(cx, self.functions, integrand, variable) else {
            return Outcome::Pass;
        };
        let graph = &mut *cx.graph;
        // An improper integral is the difference of the antiderivative's
        // limits at the infinite ends.
        let limit = graph.ops().lookup("limit");
        let at = |graph: &mut Graph, bound: NodeId, is_infinite: bool| match (is_infinite, limit) {
            | (true, Some(limit)) => Some(graph.node(limit, &[primitive, variable, bound])),
            | (true, None) => None,
            | (false, _) => Some(graph.substitute(primitive, variable, bound)),
        };
        let (Some(at_upper), Some(at_lower)) = (at(graph, upper, upper_infinite), at(graph, lower, lower_infinite)) else {
            return Outcome::Pass;
        };
        let minus_one = graph.int(-1);
        let negated = graph.node(core::MUL, &[minus_one, at_lower]);
        Outcome::Equal(graph.node(core::ADD, &[at_upper, negated]))
    }
}

/// Numeric kernel for `defint`: adaptive Gauss–Kronrod quadrature.
/// A numeric function of `inputs`.
type Numeric = Box<dyn Fn(&[f64]) -> f64>;

/// Compiles `term` as a function of `inputs`. An iterated integral
/// (`defint` whose integrand or limits are themselves evaluable this way)
/// becomes nested quadrature: the inner integral is evaluated afresh at
/// every outer sample point, with a smaller panel budget.
fn compile_nested(
    graph: &mut Graph,
    defint: OpId,
    term: NodeId,
    inputs: &[SymbolId],
    tolerance: f64,
    depth: usize,
) -> Option<Numeric> {
    if graph.op(term) == defint && depth < 3 {
        let &[integrand, variable, lower, upper] = graph.children(term) else {
            return None;
        };
        let y = graph.symbol_of(variable)?;
        let mut inner_inputs = vec![y];
        inner_inputs.extend_from_slice(inputs);
        let body = best(graph, integrand)?;
        let body = compile_nested(graph, defint, body, &inner_inputs, tolerance, depth + 1)?;
        let lower = best(graph, lower)?;
        let lower = compile_nested(graph, defint, lower, inputs, tolerance, depth + 1)?;
        let upper = best(graph, upper)?;
        let upper = compile_nested(graph, defint, upper, inputs, tolerance, depth + 1)?;
        return Some(Box::new(move |args: &[f64]| {
            let (a, b) = (lower(args), upper(args));
            let mut point = Vec::with_capacity(args.len() + 1);
            point.push(0.0);
            point.extend_from_slice(args);
            let f = |t: f64| {
                let mut point = point.clone();
                if let Some(slot) = point.first_mut() {
                    *slot = t;
                }
                body(&point)
            };
            gauss_kronrod_any(f, a, b, tolerance, 200).value
        }));
    }
    let compiled = Interpreter.compile(graph, term, inputs).ok()?;
    Some(Box::new(move |args: &[f64]| compiled.call(args)))
}

pub(super) struct Quadrature {
    pub(super) defint: OpId,
}

impl Kernel for Quadrature {
    fn ops(&self) -> Vec<OpId> {
        vec![self.defint]
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
        let &[integrand, variable, lower, upper] = graph.children(node) else {
            return Outcome::Pass;
        };
        let heavy = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).flags.has(OpFlags::HEAVY);
        if graph.enodes(graph.find(node)).any(|n| !heavy(graph, n)) {
            // A closed form exists; evaluating it is cheaper and exact.
            return Outcome::Pass;
        }
        let Some(x) = graph.symbol_of(variable) else {
            return Outcome::Pass;
        };
        let limit = |graph: &mut Graph, n: NodeId| best(graph, n).and_then(|t| graph.eval(t, cx.env));
        let (Some(a), Some(b)) = (limit(graph, lower), limit(graph, upper)) else {
            return Outcome::Pass;
        };
        let Some(term) = best(graph, integrand) else {
            return Outcome::Pass;
        };
        // Compile once: the integrand is evaluated thousands of times.
        let mut inputs = vec![x];
        let mut values = vec![0.0];
        for &(symbol, value) in cx.env.bindings() {
            if symbol != x {
                inputs.push(symbol);
                values.push(value);
            }
        }
        let tolerance = cx.env.tolerance.max(1e-13);
        let Some(compiled) = compile_nested(graph, self.defint, term, &inputs, tolerance, 0) else {
            return Outcome::Pass;
        };
        let f = |t: f64| {
            let mut args = values.clone();
            if let Some(slot) = args.first_mut() {
                *slot = t;
            }
            compiled(&args)
        };
        let result = gauss_kronrod_any(f, a, b, tolerance, 4_000);
        if result.value.is_finite() && result.error.is_finite() {
            Outcome::Approx(Ball { mid: result.value, rad: result.error })
        } else {
            Outcome::Pass
        }
    }

    fn revisit(&self) -> bool {
        // The limits or the integrand may become evaluable later.
        true
    }
}
