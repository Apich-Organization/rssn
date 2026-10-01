//! Complex numbers.
//!
//! The imaginary unit is the nullary operator `I`; complex values are
//! ordinary terms such as `2 + 3*I` or `exp(I*x)`. Nothing is assumed about
//! a symbol unless the session says so: `re(z)` of a symbol `z` stays
//! `re(z)`, while for a symbol assumed real it is `z` itself.
//!
//! | operator | value |
//! |---|---|
//! | `I` | the imaginary unit; `I^n` reduces for integer `n` |
//! | `re(z)`, `im(z)`, `conj(z)` | real part, imaginary part, conjugate |
//! | `abs(z)`, `arg(z)` | modulus, principal argument in `(-pi, pi]` |
//!
//! The kernels split a term into real and imaginary parts structurally:
//! sums and products, integer powers, positive real bases with complex
//! exponents, and `exp`, `ln`, `sin`, `cos`, `sinh`, `cosh` of complex
//! arguments. `exp(I*t)` reduces to `cos(t) + I*sin(t)` when both sides
//! have closed values (`exp(I*pi) = -1`); otherwise Euler's formula is an
//! exploring rule in both directions.
//!
//! The real evaluator cannot see `I`, so complex results carry no numeric
//! witness. [`eval_complex`] evaluates terms over the complex numbers; the
//! complex-analysis kernels and the tests use it to check symbolic results.

pub mod analysis;
pub mod branches;

use std::collections::HashMap;

use num_complex::Complex64;

use crate::graph::Always;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::best;

use super::calculus::calculus;
use super::elementary::elementary;

/// The complex-number rule set.
#[must_use]
pub fn complex() -> RuleSet {
    RuleSet::new("complex", install).needs(elementary()).needs(calculus())
}

pub use crate::graph::ComplexEval;

/// Largest integer power expanded by the splitter.
const MAX_POWER: i64 = 32;

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unit = i.op(OpDescriptor::new("I", Arity::Fixed(0)).cost(1))?;
    let re = i.op(OpDescriptor::new("re", Arity::Fixed(1)).cost(8))?;
    let im = i.op(OpDescriptor::new("im", Arity::Fixed(1)).cost(8))?;
    let conj = i.op(OpDescriptor::new("conj", Arity::Fixed(1)).cost(8))?;
    let arg = i.op(OpDescriptor::new("arg", Arity::Fixed(1)).cost(8))?;
    {
        let ops = i.graph().ops_mut();
        ops.set_attr(unit, Always(Facts::NONZERO));
        ops.set_attr(re, Always(Facts::REAL));
        ops.set_attr(im, Always(Facts::REAL));
        ops.set_attr(arg, Always(Facts::REAL));
    }
    if let Some(abs) = i.graph().ops().lookup("abs") {
        i.graph().ops_mut().set_attr(abs, Always(Facts::NONNEGATIVE));
    }
    register_complex_evals(i.graph());
    i.rewrites(
        Tier::Normalize,
        &[
            "complex/re-real: re(?x) => ?x if real(?x)",
            "complex/im-real: im(?x) => 0 if real(?x)",
            "complex/conj-real: conj(?x) => ?x if real(?x)",
            "complex/conj-conj: conj(conj(?x)) => ?x",
            // Conjugation is a field automorphism and commutes with every
            // function whose power series has real coefficients.
            "complex/conj-add: conj(?a + ?b) => conj(?a) + conj(?b)",
            "complex/conj-mul: conj(?a * ?b) => conj(?a) * conj(?b)",
            "complex/conj-pow: conj(?a ^ ?n) => conj(?a) ^ ?n if integer(?n)",
            "complex/conj-exp: conj(exp(?a)) => exp(conj(?a))",
            "complex/conj-sin: conj(sin(?a)) => sin(conj(?a))",
            "complex/conj-cos: conj(cos(?a)) => cos(conj(?a))",
            "complex/conj-sinh: conj(sinh(?a)) => sinh(conj(?a))",
            "complex/conj-cosh: conj(cosh(?a)) => cosh(conj(?a))",
            "complex/arg-positive: arg(?x) => 0 if positive(?x)",
            "complex/arg-negative: arg(?x) => pi if negative(?x)",
        ],
    )?;
    i.rewrites(
        Tier::Explore,
        &[
            "complex/euler: exp(I * ?x) <=> cos(?x) + I * sin(?x)",
        ],
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| RuleError::Invalid { rule: "complex".into(), reason: "needs elementary" })?;
    branches::install(i)?;
    analysis::install(i)?;
    i.kernel("complex/unit-power", Tier::Normalize, UnitPower { unit });
    i.kernel("complex/parts", Tier::Normalize, Parts { ops, re, im, conj, arg });
    i.kernel("complex/exp-i", Tier::Normalize, ExpI { ops });
    Ok(())
}

/// Operators the splitter needs.
#[derive(Copy, Clone, Debug)]
pub(crate) struct Ops {
    pub(crate) unit: OpId,
    pub(crate) re: OpId,
    pub(crate) im: OpId,
    pub(crate) exp: OpId,
    pub(crate) ln: OpId,
    pub(crate) sin: OpId,
    pub(crate) cos: OpId,
    pub(crate) sinh: OpId,
    pub(crate) cosh: OpId,
    pub(crate) atan2: OpId,
    pub(crate) abs: OpId,
    pub(crate) conj: OpId,
    pub(crate) arg: OpId,
}

impl Ops {
    pub(crate) fn of(graph: &Graph) -> Option<Self> {
        let get = |name: &str| graph.ops().lookup(name);
        Some(Self {
            unit: get("I")?,
            re: get("re")?,
            im: get("im")?,
            exp: get("exp")?,
            ln: get("ln")?,
            sin: get("sin")?,
            cos: get("cos")?,
            sinh: get("sinh")?,
            cosh: get("cosh")?,
            atan2: get("atan2")?,
            abs: get("abs")?,
            conj: get("conj")?,
            arg: get("arg")?,
        })
    }
}

/// Term builders shared by the complex kernels.
pub(crate) mod build {
    use super::{Graph, NodeId, core, OpId};

    pub fn add(
        graph: &mut Graph,
        terms: &[NodeId],
    ) -> NodeId {
        match terms {
            | [] => graph.int(0),
            | [only] => *only,
            | _ => graph.node(core::ADD, terms),
        }
    }

    pub fn mul(
        graph: &mut Graph,
        factors: &[NodeId],
    ) -> NodeId {
        match factors {
            | [] => graph.int(1),
            | [only] => *only,
            | _ => graph.node(core::MUL, factors),
        }
    }

    pub fn neg(
        graph: &mut Graph,
        x: NodeId,
    ) -> NodeId {
        let minus_one = graph.int(-1);
        mul(graph, &[minus_one, x])
    }

    pub fn sub(
        graph: &mut Graph,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        let negated = neg(graph, b);
        add(graph, &[a, negated])
    }

    pub fn pow(
        graph: &mut Graph,
        base: NodeId,
        exponent: NodeId,
    ) -> NodeId {
        graph.node(core::POW, &[base, exponent])
    }

    pub fn powi(
        graph: &mut Graph,
        base: NodeId,
        exponent: i64,
    ) -> NodeId {
        let e = graph.int(exponent);
        pow(graph, base, e)
    }

    pub fn call(
        graph: &mut Graph,
        op: OpId,
        args: &[NodeId],
    ) -> NodeId {
        graph.node(op, args)
    }

    /// `a + b*I`.
    pub fn complex(
        graph: &mut Graph,
        unit: OpId,
        re: NodeId,
        im: NodeId,
    ) -> NodeId {
        let i = graph.node(unit, &[]);
        let imaginary = mul(graph, &[im, i]);
        add(graph, &[re, imaginary])
    }
}

use build::add;
use build::call;
use build::mul;
use build::neg;
use build::sub;

/// Real and imaginary parts of the best form of `z`, as (unsimplified)
/// real terms, when they can be read off its structure.
pub(crate) fn split(
    graph: &mut Graph,
    ops: Ops,
    z: NodeId,
) -> Option<(NodeId, NodeId)> {
    let term = best(graph, z)?;
    let mut memo = HashMap::new();
    split_term(graph, ops, term, &mut memo, 0)
}

#[allow(clippy::too_many_lines)]
fn split_term(
    graph: &mut Graph,
    ops: Ops,
    term: NodeId,
    memo: &mut HashMap<NodeId, Option<(NodeId, NodeId)>>,
    depth: usize,
) -> Option<(NodeId, NodeId)> {
    if let Some(known) = memo.get(&term) {
        return *known;
    }
    if depth > 64 {
        return None;
    }
    let zero = graph.int(0);
    let result = if graph.facts(term).has(Facts::REAL) {
        Some((term, zero))
    } else {
        let op = graph.op(term);
        let children = graph.children(term).to_vec();
        let mut part = |graph: &mut Graph, n: NodeId| split_term(graph, ops, n, memo, depth + 1);
        match (op, children.as_slice()) {
            | (op, []) if op == ops.unit => Some((zero, graph.int(1))),
            | (core::ADD, terms) => {
                let mut res = Vec::with_capacity(terms.len());
                let mut ims = Vec::with_capacity(terms.len());
                for &t in terms {
                    let (a, b) = part(graph, t)?;
                    res.push(a);
                    ims.push(b);
                }
                Some((add(graph, &res), add(graph, &ims)))
            },
            | (core::MUL, factors) => {
                let mut acc = (graph.int(1), zero);
                for &f in factors {
                    let next = part(graph, f)?;
                    acc = product(graph, acc, next);
                }
                Some(acc)
            },
            | (core::POW, &[base, exponent]) => {
                let n = graph.number_of(exponent).and_then(Number::to_i64);
                match n {
                    | Some(n) if n != 0 && n.abs() <= MAX_POWER => {
                        let parts = part(graph, base)?;
                        let mut acc = (graph.int(1), zero);
                        for _ in 0..n.abs() {
                            acc = product(graph, acc, parts);
                        }
                        if n < 0 {
                            // 1/(a + bi) = (a - bi)/(a^2 + b^2)
                            let (a, b) = acc;
                            let a2 = build::powi(graph, a, 2);
                            let b2 = build::powi(graph, b, 2);
                            let norm = add(graph, &[a2, b2]);
                            let inverse = build::powi(graph, norm, -1);
                            let re = mul(graph, &[a, inverse]);
                            let minus_b = neg(graph, b);
                            let im = mul(graph, &[minus_b, inverse]);
                            Some((re, im))
                        } else {
                            Some(acc)
                        }
                    },
                    | _ if graph.facts(base).has(Facts::POSITIVE) => {
                        // a^(x + iy) = a^x (cos(y ln a) + i sin(y ln a))
                        let (x, y) = part(graph, exponent)?;
                        let modulus = build::pow(graph, base, x);
                        let log = call(graph, ops.ln, &[base]);
                        let angle = mul(graph, &[y, log]);
                        Some(polar(graph, ops, modulus, angle))
                    },
                    | _ => None,
                }
            },
            | (op, &[w]) if op == ops.exp => {
                let (a, b) = part(graph, w)?;
                let modulus = call(graph, ops.exp, &[a]);
                Some(polar(graph, ops, modulus, b))
            },
            | (op, &[w]) if op == ops.ln => {
                let (a, b) = part(graph, w)?;
                let modulus = modulus(graph, a, b);
                let re = call(graph, ops.ln, &[modulus]);
                let im = call(graph, ops.atan2, &[b, a]);
                Some((re, im))
            },
            | (op, &[w]) if op == ops.sin || op == ops.cos || op == ops.sinh || op == ops.cosh => {
                let (a, b) = part(graph, w)?;
                let (sin_a, cos_a) = (call(graph, ops.sin, &[a]), call(graph, ops.cos, &[a]));
                let (sinh_b, cosh_b) = (call(graph, ops.sinh, &[b]), call(graph, ops.cosh, &[b]));
                let (sinh_a, cosh_a) = (call(graph, ops.sinh, &[a]), call(graph, ops.cosh, &[a]));
                let (sin_b, cos_b) = (call(graph, ops.sin, &[b]), call(graph, ops.cos, &[b]));
                Some(if op == ops.sin {
                    (mul(graph, &[sin_a, cosh_b]), mul(graph, &[cos_a, sinh_b]))
                } else if op == ops.cos {
                    let im = mul(graph, &[sin_a, sinh_b]);
                    (mul(graph, &[cos_a, cosh_b]), neg(graph, im))
                } else if op == ops.sinh {
                    (mul(graph, &[sinh_a, cos_b]), mul(graph, &[cosh_a, sin_b]))
                } else {
                    (mul(graph, &[cosh_a, cos_b]), mul(graph, &[sinh_a, sin_b]))
                })
            },
            | (op, &[y, x]) if op == ops.atan2 => {
                // real for real arguments
                let (_, by) = part(graph, y)?;
                let (_, bx) = part(graph, x)?;
                let zero_part = |graph: &Graph, n: NodeId| graph.number_of(n).is_some_and(Number::is_zero);
                (zero_part(graph, by) && zero_part(graph, bx)).then_some((term, zero))
            },
            | (op, &[w]) if op == ops.re => {
                let (a, _) = part(graph, w)?;
                Some((a, zero))
            },
            | (op, &[w]) if op == ops.im => {
                let (_, b) = part(graph, w)?;
                Some((b, zero))
            },
            | (op, &[w]) if op == ops.conj => {
                let (a, b) = part(graph, w)?;
                Some((a, neg(graph, b)))
            },
            | (op, &[w]) if op == ops.abs => {
                let (a, b) = part(graph, w)?;
                Some((modulus(graph, a, b), zero))
            },
            | (op, &[w]) if op == ops.arg => {
                let (a, b) = part(graph, w)?;
                Some((call(graph, ops.atan2, &[b, a]), zero))
            },
            | _ => None,
        }
    };
    memo.insert(term, result);
    result
}

/// `(a + bi)(c + di)`.
fn product(
    graph: &mut Graph,
    (a, b): (NodeId, NodeId),
    (c, d): (NodeId, NodeId),
) -> (NodeId, NodeId) {
    let ac = mul(graph, &[a, c]);
    let bd = mul(graph, &[b, d]);
    let ad = mul(graph, &[a, d]);
    let bc = mul(graph, &[b, c]);
    (sub(graph, ac, bd), add(graph, &[ad, bc]))
}

/// `r (cos t + i sin t)`.
fn polar(
    graph: &mut Graph,
    ops: Ops,
    r: NodeId,
    t: NodeId,
) -> (NodeId, NodeId) {
    let cos = call(graph, ops.cos, &[t]);
    let sin = call(graph, ops.sin, &[t]);
    (mul(graph, &[r, cos]), mul(graph, &[r, sin]))
}

/// `sqrt(a^2 + b^2)`.
fn modulus(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let a2 = build::powi(graph, a, 2);
    let b2 = build::powi(graph, b, 2);
    let sum = add(graph, &[a2, b2]);
    let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
    build::pow(graph, sum, half)
}

/// Whether `node` contains one of `ops` anywhere.
fn mentions(
    graph: &Graph,
    node: NodeId,
    ops: &[OpId],
) -> bool {
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if ops.contains(&graph.op(n)) {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// `I^n` for integer `n`.
struct UnitPower {
    unit: OpId,
}

impl Kernel for UnitPower {
    fn ops(&self) -> Vec<OpId> {
        vec![core::POW]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[base, exponent] = graph.children(node) else {
            return Outcome::Pass;
        };
        if !graph.enodes(graph.find(base)).any(|n| graph.op(n) == self.unit) {
            return Outcome::Pass;
        }
        let Some(n) = graph.number_of(exponent).and_then(Number::to_i64) else {
            return Outcome::Pass;
        };
        let i = graph.node(self.unit, &[]);
        // Pinned: `I^3` and `-I` cost the same, and the power must not win.
        Outcome::Pinned(match n.rem_euclid(4) {
            | 0 => graph.int(1),
            | 1 => i,
            | 2 => graph.int(-1),
            | _ => neg(graph, i),
        })
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `re`, `im`, `conj`, `abs` and `arg` of terms whose parts are known.
struct Parts {
    ops: Ops,
    re: OpId,
    im: OpId,
    conj: OpId,
    arg: OpId,
}

impl Kernel for Parts {
    fn ops(&self) -> Vec<OpId> {
        vec![self.re, self.im, self.conj, self.arg, self.ops.abs]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let &[z] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(term) = best(cx.graph, z) else {
            return Outcome::Pass;
        };
        // Symbols of unknown reality have no parts to read off.
        if cx.graph.as_symbol(term).is_some() {
            return Outcome::Pass;
        }
        if op == self.ops.abs && cx.graph.facts(term).has(Facts::REAL) {
            return Outcome::Pass;
        }
        let Some((a, b)) = split(cx.graph, self.ops, term) else {
            return Outcome::Pass;
        };
        let (a, b) = (cx.simplify(a), cx.simplify(b));
        let graph = &mut *cx.graph;
        let result = if op == self.re {
            a
        } else if op == self.im {
            b
        } else if op == self.conj {
            let minus_b = neg(graph, b);
            build::complex(graph, self.ops.unit, a, minus_b)
        } else if op == self.arg {
            // On the imaginary axis atan2 has no closed form here.
            let a_zero = graph.number_of(a).is_some_and(Number::is_zero);
            let b_sign = graph.facts(b);
            match (a_zero, b_sign.has(Facts::POSITIVE), b_sign.has(Facts::NEGATIVE)) {
                | (true, true, _) | (true, _, true) => {
                    let pi = graph.ops().lookup("pi").map(|p| graph.node(p, &[]));
                    let Some(pi) = pi else {
                        return Outcome::Pass;
                    };
                    let half = graph.num(Number::fraction(if b_sign.has(Facts::POSITIVE) { 1 } else { -1 }, 2).unwrap_or_else(|| Number::from(0)));
                    mul(graph, &[half, pi])
                },
                | _ => call(graph, self.ops.atan2, &[b, a]),
            }
        } else {
            modulus(graph, a, b)
        };
        // Splitting a term that already is `re`/`im` of something gains
        // nothing and would only cycle.
        if mentions(graph, result, &[op]) && mentions(graph, term, &[self.ops.re, self.ops.im]) {
            return Outcome::Pass;
        }
        Outcome::Equal(result)
    }
}

/// `exp(z) = exp(a) (cos b + i sin b)` when the trigonometric values have
/// closed forms (`exp(I*pi) = -1`, `exp(I*pi/2) = I`).
struct ExpI {
    ops: Ops,
}

impl Kernel for ExpI {
    fn ops(&self) -> Vec<OpId> {
        vec![self.ops.exp]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[w] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(term) = best(cx.graph, w) else {
            return Outcome::Pass;
        };
        if !mentions(cx.graph, term, &[self.ops.unit]) {
            return Outcome::Pass;
        }
        let Some((a, b)) = split(cx.graph, self.ops, term) else {
            return Outcome::Pass;
        };
        let b = cx.simplify(b);
        let cos = call(cx.graph, self.ops.cos, &[b]);
        let sin = call(cx.graph, self.ops.sin, &[b]);
        let (cos, sin) = (cx.simplify(cos), cx.simplify(sin));
        if mentions(cx.graph, cos, &[self.ops.cos, self.ops.sin]) || mentions(cx.graph, sin, &[self.ops.cos, self.ops.sin])
        {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let rotation = build::complex(graph, self.ops.unit, cos, sin);
        let a_is_zero = graph.number_of(a).is_some_and(Number::is_zero);
        if a_is_zero {
            return Outcome::Equal(rotation);
        }
        let scale = call(graph, self.ops.exp, &[a]);
        Outcome::Equal(mul(graph, &[scale, rotation]))
    }
}

/// Registers [`ComplexEval`] for the elementary functions and the parts
/// operators.
fn register_complex_evals(graph: &mut Graph) {
    macro_rules! unary {
        ($f:expr) => {
            ComplexEval(|a: &[Complex64]| a.first().map($f))
        };
    }
    let table: [(&str, ComplexEval); 26] = [
        ("I", ComplexEval(|_| Some(Complex64::new(0.0, 1.0)))),
        ("exp", unary!(|z: &Complex64| z.exp())),
        ("ln", unary!(|z: &Complex64| z.ln())),
        ("sqrt", unary!(|z: &Complex64| z.sqrt())),
        ("sin", unary!(|z: &Complex64| z.sin())),
        ("cos", unary!(|z: &Complex64| z.cos())),
        ("tan", unary!(|z: &Complex64| z.tan())),
        ("sinh", unary!(|z: &Complex64| z.sinh())),
        ("cosh", unary!(|z: &Complex64| z.cosh())),
        ("tanh", unary!(|z: &Complex64| z.tanh())),
        ("asin", unary!(|z: &Complex64| z.asin())),
        ("acos", unary!(|z: &Complex64| z.acos())),
        ("atan", unary!(|z: &Complex64| z.atan())),
        ("asinh", unary!(|z: &Complex64| z.asinh())),
        ("acosh", unary!(|z: &Complex64| z.acosh())),
        ("atanh", unary!(|z: &Complex64| z.atanh())),
        ("sec", unary!(|z: &Complex64| z.cos().inv())),
        ("csc", unary!(|z: &Complex64| z.sin().inv())),
        ("cot", unary!(|z: &Complex64| z.tan().inv())),
        ("abs", unary!(|z: &Complex64| Complex64::new(z.norm(), 0.0))),
        ("re", unary!(|z: &Complex64| Complex64::new(z.re, 0.0))),
        ("im", unary!(|z: &Complex64| Complex64::new(z.im, 0.0))),
        ("conj", unary!(|z: &Complex64| z.conj())),
        ("arg", unary!(|z: &Complex64| Complex64::new(z.arg(), 0.0))),
        ("pi", ComplexEval(|_| Some(Complex64::new(std::f64::consts::PI, 0.0)))),
        ("E", ComplexEval(|_| Some(Complex64::new(std::f64::consts::E, 0.0)))),
    ];
    for (name, eval) in table {
        if let Some(op) = graph.ops().lookup(name) {
            graph.ops_mut().set_attr(op, eval);
        }
    }
}

/// Evaluates the concrete term `node` over the complex numbers; see
/// [`Graph::eval_complex`].
#[must_use]
#[allow(clippy::implicit_hasher)] // delegates to `Graph::eval_complex`, which fixes the default hasher
pub fn eval_complex(
    graph: &Graph,
    node: NodeId,
    bindings: &HashMap<SymbolId, Complex64>,
) -> Option<Complex64> {
    graph.eval_complex(node, bindings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[complex()], src)
    }

    #[test]
    fn unit_arithmetic() {
        assert_eq!(run("I^2"), "-1");
        assert_eq!(run("I*I*I"), "-I");
        assert_eq!(run("I^7"), "-I");
        assert_eq!(run("re((1 + I)*(1 - I))"), "2");
        assert_eq!(run("re((2 + 3*I)*(1 - I))"), "5");
        assert_eq!(run("im((2 + 3*I)*(1 - I))"), "1");
        assert_eq!(run("im(1/(1 + I))"), "-1/2");
        assert_eq!(run("conj(2 + 3*I)"), "2 - 3*I");
        assert_eq!(run("abs(3 + 4*I)"), "5");
        assert_eq!(run("arg(I)"), "1/2*pi");
        assert_eq!(run("arg(-2)"), "pi");
    }

    #[test]
    fn euler() {
        assert_eq!(run("exp(I*pi)"), "-1");
        assert_eq!(run("exp(I*pi/2)"), "I");
        assert_eq!(run("exp(2*pi*I)"), "1");
        assert_eq!(run("re(exp(I*x))"), "re(exp(x*I))", "x might not be real");
    }

    #[test]
    fn parts_of_real_symbols() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[complex()]).unwrap_or_else(|e| panic!("{e}"));
        let x = g.interner_mut().symbol("x");
        let y = g.interner_mut().symbol("y");
        g.assume(x, Facts::REAL);
        g.assume(y, Facts::REAL);
        let cases = [
            ("re(exp(I*x))", "cos(x)"),
            ("im(exp(x + I*y))", "exp(x)*sin(y)"),
            ("re(x + I*y)", "x"),
            ("conj(x + I*y)", "x - y*I"),
            ("im(sin(x + I*y))", "cos(x)*sinh(y)"),
            ("abs(exp(I*x))", "1"),
        ];
        for (src, want) in cases {
            let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
            engine.run(&mut g, &[root], &crate::graph::Env::symbolic(), &crate::graph::Saturate, &crate::graph::Budget::default());
            let got = crate::graph::Extractor::new(&g, &[root], &crate::graph::SizeCost).build(&mut g, root).map(|n| g.display(n));
            assert_eq!(got.as_deref(), Some(want), "{src}");
        }
    }

    #[test]
    fn complex_evaluation() {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[complex()]).is_ok());
        let root = g.parse("exp(I*z) * (1 + I)^2 + ln(z)").unwrap_or_else(|e| panic!("{e}"));
        let z = g.interner_mut().symbol("z");
        let at = Complex64::new(0.3, -0.4);
        let bindings: HashMap<SymbolId, Complex64> = std::iter::once((z, at)).collect();
        let got = eval_complex(&g, root, &bindings).unwrap_or_default();
        let i = Complex64::new(0.0, 1.0);
        let want = (i * at).exp() * (1.0 + i) * (1.0 + i) + at.ln();
        assert!((got - want).norm() < 1e-12, "{got} vs {want}");
    }

    /// Splitting is an identity: both sides agree at random complex points
    /// for terms whose parts are computed.
    #[test]
    fn split_agrees_with_evaluation() {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[complex()]).is_ok());
        let ops = Ops::of(&g).unwrap_or_else(|| panic!("ops"));
        for src in ["exp(2 + 3*I)", "sin(1 - 2*I)", "cosh(0.5 + I)", "ln(-1 + I)", "(2 + I)^-3", "2^(1 + I)", "cos(I)*sinh(1 + I)"] {
            let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
            let (a, b) = split(&mut g, ops, root).unwrap_or_else(|| panic!("no split for {src}"));
            let none = HashMap::new();
            let whole = eval_complex(&g, root, &none).unwrap_or_default();
            let a = eval_complex(&g, a, &none).unwrap_or_default();
            let b = eval_complex(&g, b, &none).unwrap_or_default();
            assert!(a.im.abs() < 1e-12 && b.im.abs() < 1e-12, "{src}");
            assert!((Complex64::new(a.re, b.re) - whole).norm() < 1e-12, "{src}: {whole} vs {a} + i{b}");
        }
    }
}
