//! Truncated power-series arithmetic.
//!
//! Expanding by repeated differentiation is exponential in the order and
//! breaks down at removable singularities (every coefficient of
//! `sin(x)/x` becomes a limit). This module expands the way a person does:
//! the series of a sum, product, quotient, power or elementary function of
//! series is computed from the series of its parts. Laurent series come
//! for free — a series carries its valuation, so `1/x` is simply the
//! series `x^(-1)`.
//!
//! Coefficients are terms of the graph. Numbers are folded on the fly;
//! symbolic coefficients are built as terms and simplified once at the
//! end. Where a leading coefficient must be known to be non-zero (to
//! divide by it), it is simplified and tested.
//! Composition uses the positive valuation of the nonconstant part to
//! bound the number of Taylor weights. The inverse hyperbolic sine and
//! tangent support direct series when their argument has zero constant term.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::SymbolId;
use crate::graph::op::core;

/// `t^valuation * sum_k coefficients[k] t^k`, known up to and including
/// `t^order`.
#[derive(Clone, Debug)]
struct Series {
    valuation: i64,
    coefficients: Vec<NodeId>,
}

impl Series {
    fn constant(c: NodeId) -> Self {
        Self {
            valuation: 0,
            coefficients: vec![c],
        }
    }
}

struct Expander<'c, 'a> {
    cx: &'c mut Cx<'a>,
    t: NodeId,
    symbol: SymbolId,
    order: i64,
    depth: usize,
}

const MAX_DEPTH: usize = 64;

impl Expander<'_, '_> {
    const fn graph(&mut self) -> &mut Graph {
        self.cx.graph
    }

    fn number(
        &self,
        node: NodeId,
    ) -> Option<Number> {
        self.cx.graph.number_of(node).cloned()
    }

    fn num(
        &mut self,
        n: Number,
    ) -> NodeId {
        self.graph().num(n)
    }

    fn zero(&mut self) -> NodeId {
        self.graph().int(0)
    }

    fn one(&mut self) -> NodeId {
        self.graph().int(1)
    }

    fn add(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        match (self.number(a), self.number(b)) {
            | (Some(x), Some(y)) => self.num(x.add(&y)),
            | (Some(x), _) if x.is_zero() => b,
            | (_, Some(y)) if y.is_zero() => a,
            | _ => self.graph().node(core::ADD, &[a, b]),
        }
    }

    fn mul(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        match (self.number(a), self.number(b)) {
            | (Some(x), Some(y)) => self.num(x.mul(&y)),
            | (Some(x), _) if x.is_zero() => a,
            | (_, Some(y)) if y.is_zero() => b,
            | (Some(x), _) if x.is_one() => b,
            | (_, Some(y)) if y.is_one() => a,
            | _ => self.graph().node(core::MUL, &[a, b]),
        }
    }

    fn scale(
        &mut self,
        a: NodeId,
        factor: &Number,
    ) -> NodeId {
        let f = self.num(factor.clone());
        self.mul(a, f)
    }

    fn reciprocal(
        &mut self,
        a: NodeId,
    ) -> Option<NodeId> {
        if let Some(x) = self.number(a) {
            return x.recip().map(|r| self.num(r));
        }
        let minus_one = self.graph().int(-1);
        Some(self.graph().node(core::POW, &[a, minus_one]))
    }

    fn is_zero(
        &mut self,
        a: NodeId,
    ) -> bool {
        match self.number(a) {
            | Some(x) => x.is_zero(),
            | None => self.cx.is_zero(a),
        }
    }

    /// Drops leading zero coefficients (raising the valuation).
    fn normalise(
        &mut self,
        mut s: Series,
    ) -> Series {
        let zeros = s
            .coefficients
            .iter()
            .take_while(|&&coefficient| self.is_zero(coefficient))
            .count();
        s.coefficients.drain(..zeros);
        s.valuation += i64::try_from(zeros).unwrap_or(0);
        s
    }

    /// Number of coefficients needed for a series of this valuation.
    fn length(
        &self,
        valuation: i64,
    ) -> usize {
        usize::try_from((self.order - valuation + 1).max(0)).unwrap_or(0)
    }

    fn sum(
        &mut self,
        a: &Series,
        b: &Series,
    ) -> Series {
        let valuation = a.valuation.min(b.valuation);
        let length = self.length(valuation);
        let mut out = Vec::with_capacity(length);
        for k in 0..length {
            let exponent = valuation + i64::try_from(k).unwrap_or(0);
            let pick = |s: &Series| {
                usize::try_from(exponent - s.valuation)
                    .ok()
                    .and_then(|i| s.coefficients.get(i).copied())
            };
            let term = match (pick(a), pick(b)) {
                | (Some(x), Some(y)) => self.add(x, y),
                | (Some(x), None) | (None, Some(x)) => x,
                | (None, None) => self.zero(),
            };
            out.push(term);
        }
        Series {
            valuation,
            coefficients: out,
        }
    }

    fn product(
        &mut self,
        a: &Series,
        b: &Series,
    ) -> Series {
        let valuation = a.valuation + b.valuation;
        let length = self.length(valuation);
        let mut out = Vec::with_capacity(length);
        let right_len = b.coefficients.len();
        for k in 0..length {
            let mut acc = self.zero();
            if right_len > 0 {
                let first = k.saturating_sub(right_len - 1);
                let end = (k + 1).min(a.coefficients.len());
                for i in first..end {
                    let term = self.mul(a.coefficients[i], b.coefficients[k - i]);
                    acc = self.add(acc, term);
                }
            }
            out.push(acc);
        }
        Series {
            valuation,
            coefficients: out,
        }
    }

    fn inverse(
        &mut self,
        a: &Series,
    ) -> Option<Series> {
        let a = self.normalise(a.clone());
        let lead = *a.coefficients.first()?;
        let inv_lead = self.reciprocal(lead)?;
        let valuation = -a.valuation;
        let length = self.length(valuation);
        let mut out: Vec<NodeId> = Vec::with_capacity(length);
        for k in 0..length {
            if k == 0 {
                out.push(inv_lead);
                continue;
            }
            // b_k = -(sum_{j=1..k} a_j b_{k-j}) / a_0
            let mut acc = self.zero();
            for j in 1..=k {
                if let (Some(&x), Some(&y)) = (a.coefficients.get(j), out.get(k - j)) {
                    let term = self.mul(x, y);
                    acc = self.add(acc, term);
                }
            }
            let negated = self.scale(acc, &Number::from(-1));
            out.push(self.mul(negated, inv_lead));
        }
        Some(Series {
            valuation,
            coefficients: out,
        })
    }

    fn integer_power(
        &mut self,
        a: &Series,
        n: i64,
    ) -> Option<Series> {
        let base = if n < 0 {
            self.inverse(a)?
        } else {
            self.normalise(a.clone())
        };
        let mut result = Series::constant(self.one_node());
        let mut power = base;
        let mut e = n.unsigned_abs();
        while e > 0 {
            if e & 1 == 1 {
                result = self.product(&result, &power);
            }
            e >>= 1;
            if e > 0 {
                power = self.product(&power, &power);
            }
        }
        Some(result)
    }

    fn one_node(&mut self) -> NodeId {
        self.one()
    }

    /// Evaluates an elementary series in an argument of positive valuation.
    fn compose(
        &mut self,
        kind: &Weights,
        u: &Series,
    ) -> Series {
        let n = if u.coefficients.is_empty() {
            0
        } else {
            usize::try_from(self.order / u.valuation).unwrap_or(0)
        };
        let weights = Self::weights(kind, n);
        let mut result = Series {
            valuation: 0,
            coefficients: Vec::new(),
        };
        let mut power = Series::constant(self.one_node());
        for (k, weight) in weights.iter().enumerate() {
            if k > 0 {
                power = self.product(&power, u);
            }
            if power.valuation > self.order {
                break;
            }
            if weight.is_zero() {
                continue;
            }
            let scaled = Series {
                valuation: power.valuation,
                coefficients: power
                    .coefficients
                    .iter()
                    .map(|&c| self.scale(c, weight))
                    .collect(),
            };
            result = self.sum(&result, &scaled);
        }
        result
    }

    /// Splits a series of non-negative valuation into its constant term
    /// and the rest (of positive valuation).
    fn split_constant(
        &mut self,
        s: &Series,
    ) -> Option<(NodeId, Series)> {
        if s.valuation < 0 {
            return None;
        }
        if s.valuation > 0 {
            return Some((self.zero(), s.clone()));
        }
        let (&constant, terms) = s.coefficients.split_first()?;
        let zeros = terms
            .iter()
            .take_while(|&&c| self.cx.graph.number_of(c).is_some_and(Number::is_zero))
            .count();
        Some((
            constant,
            Series {
                valuation: 1 + i64::try_from(zeros).unwrap_or(0),
                coefficients: terms[zeros..].to_vec(),
            },
        ))
    }

    /// Taylor weights of the elementary functions at zero, up to `n`.
    fn weights(
        kind: &Weights,
        n: usize,
    ) -> Vec<Number> {
        let unit_fraction =
            |sign: i32, denominator| Number::rat(BigRational::new_raw(sign.into(), denominator));
        match kind {
            | Weights::Exp | Weights::Sin | Weights::Cos | Weights::Sinh | Weights::Cosh => {
                let mut factorial = BigInt::one();
                (0..=n)
                    .map(|k| {
                        if k > 0 {
                            factorial *= k;
                        }
                        let sign = match kind {
                            | Weights::Exp => 1,
                            | Weights::Sin | Weights::Sinh if k % 2 == 1 => {
                                if matches!(kind, Weights::Sin) && k % 4 == 3 {
                                    -1
                                } else {
                                    1
                                }
                            },
                            | Weights::Cos | Weights::Cosh if k.is_multiple_of(2) => {
                                if matches!(kind, Weights::Cos) && k % 4 == 2 {
                                    -1
                                } else {
                                    1
                                }
                            },
                            | _ => 0,
                        };
                        if sign == 0 {
                            Number::from(0)
                        } else {
                            unit_fraction(sign, factorial.clone())
                        }
                    })
                    .collect()
            },
            | Weights::Log1p | Weights::Atan | Weights::Atanh => {
                (0..=n)
                    .map(|k| {
                        let sign = match kind {
                            | Weights::Log1p if k > 0 => {
                                if k % 2 == 1 {
                                    1
                                } else {
                                    -1
                                }
                            },
                            | Weights::Atan | Weights::Atanh if k % 2 == 1 => {
                                if matches!(kind, Weights::Atan) && k % 4 == 3 {
                                    -1
                                } else {
                                    1
                                }
                            },
                            | _ => 0,
                        };
                        if sign == 0 {
                            Number::from(0)
                        } else {
                            unit_fraction(sign, BigInt::from(k))
                        }
                    })
                    .collect()
            },
            | Weights::Asin | Weights::Asinh => {
                // Advance the two integer factors of binomial(2m,m)/4^m,
                // the derivative's coefficient at even degree 2m.
                let mut central = BigInt::one();
                let mut four_to_m = BigInt::one();
                (0..=n)
                    .map(|k| {
                        if k.is_multiple_of(2) {
                            return Number::from(0);
                        }
                        let numerator = if matches!(kind, Weights::Asinh) && k % 4 == 3 {
                            -&central
                        } else {
                            central.clone()
                        };
                        let weight = BigRational::new(numerator, &four_to_m * k);
                        central *= k;
                        central <<= 1;
                        central /= k.div_ceil(2);
                        four_to_m <<= 2;
                        Number::rat(weight)
                    })
                    .collect()
            },
            | Weights::Binomial(exponent) => Self::binomial_weights(exponent, n),
        }
    }

    fn binomial_weights(
        exponent: &Number,
        n: usize,
    ) -> Vec<Number> {
        if let Some(exact) = exponent.to_rational() {
            let mut coefficient = BigRational::one();
            let mut numerator = exact.numer().clone();
            let denominator = exact.denom();
            let mut weights = Vec::with_capacity(n + 1);
            weights.push(Number::from(1));
            for k in 1..=n {
                coefficient = BigRational::new(
                    coefficient.numer() * &numerator,
                    coefficient.denom() * denominator * k,
                );
                weights.push(Number::rat(coefficient.clone()));
                numerator -= denominator;
            }
            return weights;
        }

        // For approximate exponents use the original multiplication order.
        let mut factorial = BigInt::one();
        (0..=n)
            .map(|k| {
                if k > 0 {
                    factorial *= k;
                }
                let mut value = Number::from(1);
                for i in 0..k {
                    let factor = exponent.add(&Number::from(-i64::try_from(i).unwrap_or(0)));
                    value = value.mul(&factor);
                }
                value.mul(&Number::rat(BigRational::new_raw(
                    BigInt::one(),
                    factorial.clone(),
                )))
            })
            .collect()
    }

    fn function_of(
        &mut self,
        name: &str,
        arg: NodeId,
    ) -> Option<NodeId> {
        let op = self.cx.graph.ops().lookup(name)?;
        self.cx.graph.try_node(op, &[arg])
    }

    /// The series of the concrete term `f` in `t`.
    fn expand(
        &mut self,
        f: NodeId,
    ) -> Option<Series> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return None;
        }
        let result = self.expand_inner(f);
        self.depth -= 1;
        result
    }

    fn expand_inner(
        &mut self,
        f: NodeId,
    ) -> Option<Series> {
        let graph = &*self.cx.graph;
        if !graph.depends_on(graph.find(f), self.symbol) {
            return Some(Series::constant(f));
        }
        if f == self.t {
            return Some(Series {
                valuation: 1,
                coefficients: vec![self.one()],
            });
        }
        let op = graph.op(f);
        let children = graph.children(f).to_vec();
        if op == core::ADD {
            let mut acc = Series {
                valuation: self.order + 1,
                coefficients: Vec::new(),
            };
            for child in children {
                let s = self.expand(child)?;
                acc = self.sum(&acc, &s);
            }
            return Some(acc);
        }
        if op == core::MUL {
            let mut acc = Series::constant(self.one_node());
            for child in children {
                let s = self.expand(child)?;
                acc = self.product(&acc, &s);
            }
            return Some(acc);
        }
        if op == core::POW {
            let &[base, exponent] = children.as_slice() else {
                return None;
            };
            let graph = &*self.cx.graph;
            if graph.depends_on(graph.find(exponent), self.symbol) {
                // b^e = exp(e ln b)
                let log = self.function_of("ln", base)?;
                let product = self.cx.graph.node(core::MUL, &[exponent, log]);
                let exp = self.function_of("exp", product)?;
                return self.expand(exp);
            }
            let e = self.number(exponent)?;
            let s = self.expand(base)?;
            if let Some(n) = e.to_i64() {
                return self.integer_power(&s, n);
            }
            // Non-integer exponent: c0^e t^(v e) (1 + u)^e.
            let s = self.normalise(s);
            let v_times_e = Number::from(s.valuation).mul(&e);
            let shift = v_times_e.to_i64()?;
            let lead = *s.coefficients.first()?;
            let inv_lead = self.reciprocal(lead)?;
            let ratio = Series {
                valuation: 0,
                coefficients: s
                    .coefficients
                    .iter()
                    .map(|&c| self.mul(c, inv_lead))
                    .collect(),
            };
            let (_, u) = self.split_constant(&ratio)?;
            let mut body = self.compose(&Weights::Binomial(e.clone()), &u);
            let e_node = self.num(e);
            let lead_power = self.cx.graph.node(core::POW, &[lead, e_node]);
            body.coefficients = body
                .coefficients
                .iter()
                .map(|&c| self.mul(lead_power, c))
                .collect();
            body.valuation += shift;
            body.coefficients.truncate(self.length(body.valuation));
            return Some(body);
        }
        let &[arg] = children.as_slice() else {
            return None;
        };
        let name = self.cx.graph.ops().get(op).name.to_string();
        let s = self.expand(arg)?;
        match name.as_str() {
            | "exp" => {
                let (c0, u) = self.split_constant(&s)?;
                let body = self.compose(&Weights::Exp, &u);
                let factor = self.function_of("exp", c0)?;
                Some(self.scaled(&body, factor))
            },
            | "ln" => {
                let s = self.normalise(s);
                if s.valuation != 0 {
                    return None;
                }
                let lead = *s.coefficients.first()?;
                let inv_lead = self.reciprocal(lead)?;
                let ratio = Series {
                    valuation: 0,
                    coefficients: s
                        .coefficients
                        .iter()
                        .map(|&c| self.mul(c, inv_lead))
                        .collect(),
                };
                let (_, u) = self.split_constant(&ratio)?;
                let body = self.compose(&Weights::Log1p, &u);
                let log_lead = self.function_of("ln", lead)?;
                let constant = Series::constant(log_lead);
                Some(self.sum(&constant, &body))
            },
            | "sin" | "cos" | "sinh" | "cosh" => {
                let hyperbolic = name.ends_with('h');
                let (c0, u) = self.split_constant(&s)?;
                let (odd_kind, even_kind) = if hyperbolic {
                    (Weights::Sinh, Weights::Cosh)
                } else {
                    (Weights::Sin, Weights::Cos)
                };
                let odd = self.compose(&odd_kind, &u);
                let even = self.compose(&even_kind, &u);
                let (s_name, c_name) = if hyperbolic {
                    ("sinh", "cosh")
                } else {
                    ("sin", "cos")
                };
                let (s0, k0) = (self.function_of(s_name, c0)?, self.function_of(c_name, c0)?);
                if name == s_name {
                    // sin(c0 + u) = sin(c0) C(u) + cos(c0) S(u), likewise sinh.
                    let a = self.scaled(&even, s0);
                    let b = self.scaled(&odd, k0);
                    Some(self.sum(&a, &b))
                } else {
                    // cos(c0 + u) = cos(c0) C(u) - sin(c0) S(u); cosh has +.
                    let a = self.scaled(&even, k0);
                    let sign = Number::from(if hyperbolic { 1 } else { -1 });
                    let signed = self.scale(s0, &sign);
                    let b = self.scaled(&odd, signed);
                    Some(self.sum(&a, &b))
                }
            },
            | "tan" | "tanh" => {
                let (s_name, c_name) = if name == "tan" {
                    ("sin", "cos")
                } else {
                    ("sinh", "cosh")
                };
                let (top, bottom) = (
                    self.function_of(s_name, arg)?,
                    self.function_of(c_name, arg)?,
                );
                let (top, bottom) = (self.expand(top)?, self.expand(bottom)?);
                let inverse = self.inverse(&bottom)?;
                Some(self.product(&top, &inverse))
            },
            | "atan" | "asin" | "atanh" | "asinh" => {
                let (c0, u) = self.split_constant(&s)?;
                if !self.is_zero(c0) {
                    return None;
                }
                let kind = match name.as_str() {
                    | "atan" => Weights::Atan,
                    | "asin" => Weights::Asin,
                    | "atanh" => Weights::Atanh,
                    | _ => Weights::Asinh,
                };
                Some(self.compose(&kind, &u))
            },
            | "sqrt" => {
                let half = self.num(Number::fraction(1, 2)?);
                let power = self.cx.graph.node(core::POW, &[arg, half]);
                self.expand(power)
            },
            | _ => None,
        }
    }

    fn scaled(
        &mut self,
        s: &Series,
        factor: NodeId,
    ) -> Series {
        Series {
            valuation: s.valuation,
            coefficients: s
                .coefficients
                .iter()
                .map(|&c| self.mul(factor, c))
                .collect(),
        }
    }
}

#[derive(Clone)]
enum Weights {
    Exp,
    Sin,
    Cos,
    Sinh,
    Cosh,
    Log1p,
    Atan,
    Atanh,
    Asin,
    Asinh,
    Binomial(Number),
}

/// The Laurent expansion of `f` about `x = a` up to and including
/// `(x - a)^order`: returns the valuation and the simplified coefficients
/// from that power upward.
pub fn laurent_expansion(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    a: NodeId,
    order: i64,
) -> Option<(i64, Vec<NodeId>)> {
    let fresh = cx.graph.interner_mut().fresh_symbol("t");
    let t = cx.graph.symbol_node(fresh);
    let shifted_point = if cx.graph.number_of(a).is_some_and(Number::is_zero) {
        t
    } else {
        cx.graph.node(core::ADD, &[a, t])
    };
    let term = cx.simplify(f);
    let in_t = cx.graph.substitute(term, x, shifted_point);
    // A little extra precision covers cancellation in quotients.
    let mut expander = Expander {
        cx,
        t,
        symbol: fresh,
        order: order + 4,
        depth: 0,
    };
    let series = expander.expand(in_t)?;
    let series = expander.normalise(series);
    let keep = usize::try_from((order - series.valuation + 1).max(0)).unwrap_or(0);
    let mut coefficients = Vec::with_capacity(keep);
    for &c in series.coefficients.iter().take(keep) {
        let simplified = expander.cx.simplify(c);
        if expander
            .cx
            .graph
            .depends_on(expander.cx.graph.find(simplified), fresh)
        {
            return None;
        }
        coefficients.push(simplified);
    }
    // Leading coefficients that only simplified to zero now.
    let zeros = coefficients
        .iter()
        .take_while(|&&c| expander.cx.graph.number_of(c).is_some_and(Number::is_zero))
        .count();
    coefficients.drain(..zeros);
    let valuation = series.valuation + i64::try_from(zeros).unwrap_or(0);
    Some((valuation, coefficients))
}

#[cfg(test)]
mod tests {
    use crate::graph::Engine;
    use crate::graph::Env;

    use super::*;

    #[test]
    fn inverse_hyperbolic_taylor_weights_are_exact() {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &crate::rules::standard()).unwrap();
        let x = graph.sym("x");
        let zero = graph.int(0);
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        for (function, odd) in [
            (
                "asinh(x)",
                [(1, 1), (-1, 6), (3, 40), (-5, 112), (35, 1152)],
            ),
            ("atanh(x)", [(1, 1), (1, 3), (1, 5), (1, 7), (1, 9)]),
        ] {
            let term = cx.graph.parse(function).unwrap();
            let (valuation, coefficients) = laurent_expansion(&mut cx, term, x, zero, 9).unwrap();
            assert_eq!(valuation, 1);
            assert_eq!(coefficients.len(), 9);
            for (k, &coefficient) in coefficients.iter().enumerate() {
                let expected = if k.is_multiple_of(2) {
                    let (n, d) = odd[k / 2];
                    Number::fraction(n, d).unwrap()
                } else {
                    Number::from(0)
                };
                assert_eq!(
                    cx.graph.number_of(coefficient),
                    Some(&expected),
                    "{function}"
                );
            }
        }
    }

    #[test]
    fn sparse_composition_keeps_valuation_and_through_order_terms() {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &crate::rules::standard()).unwrap();
        let x = graph.sym("x");
        let zero = graph.int(0);
        let term = graph.parse("atanh(x^4+x^5)").unwrap();
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        let (valuation, coefficients) = laurent_expansion(&mut cx, term, x, zero, 14).unwrap();
        assert_eq!(valuation, 4);
        assert_eq!(coefficients.len(), 11);
        let expected = [
            (1, 1),
            (1, 1),
            (0, 1),
            (0, 1),
            (0, 1),
            (0, 1),
            (0, 1),
            (0, 1),
            (1, 3),
            (1, 1),
            (1, 1),
        ];
        for (&coefficient, (n, d)) in coefficients.iter().zip(expected) {
            assert_eq!(
                cx.graph.number_of(coefficient),
                Number::fraction(n, d).as_ref(),
            );
        }
    }
}
