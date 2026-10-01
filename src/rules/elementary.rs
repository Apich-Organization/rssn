//! Elementary functions: exponentials, logarithms, trigonometric and
//! hyperbolic functions.
//!
//! Also the identities that never enlarge a
//! term. Structure-changing identities (angle addition, double angles) are
//! in the exploring tier.

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OnReals;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::EvalFn;
use crate::graph::op::core;
use crate::graph::rule::Installer;

use super::arith::arith;

/// Symmetry of a unary function under negation of its argument.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Parity {
    /// `f(-x) = f(x)`.
    Even,
    /// `f(-x) = -f(x)`.
    Odd,
}

/// The elementary-function rule set.
#[must_use]
pub fn elementary() -> RuleSet {
    RuleSet::new("elementary", install).needs(arith())
}

/// Lifts a unary `f64` method into an operator's scalar semantics.
macro_rules! unary {
    ($name:ident) => {{
        let eval: EvalFn = |a| a.first().map_or(f64::NAN, |x| x.$name());
        eval
    }};
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    // name, scalar semantics, parity, what is known of the value for real
    // arguments
    /// A unary function given by a closure on `f64`.
    macro_rules! closure {
        ($f:expr) => {{
            let eval: EvalFn = |a| a.first().map_or(f64::NAN, $f);
            eval
        }};
    }
    let functions: [(&str, EvalFn, Option<Parity>, Facts); 28] = [
        ("exp", unary!(exp), None, Facts::POSITIVE),
        ("ln", unary!(ln), None, Facts::NONE),
        ("sin", unary!(sin), Some(Parity::Odd), Facts::REAL),
        ("cos", unary!(cos), Some(Parity::Even), Facts::REAL),
        ("tan", unary!(tan), Some(Parity::Odd), Facts::REAL),
        ("asin", unary!(asin), Some(Parity::Odd), Facts::NONE),
        ("acos", unary!(acos), None, Facts::NONE),
        ("atan", unary!(atan), Some(Parity::Odd), Facts::REAL),
        ("sinh", unary!(sinh), Some(Parity::Odd), Facts::REAL),
        ("cosh", unary!(cosh), Some(Parity::Even), Facts::POSITIVE),
        ("tanh", unary!(tanh), Some(Parity::Odd), Facts::REAL),
        ("sqrt", unary!(sqrt), None, Facts::NONE),
        ("abs", unary!(abs), Some(Parity::Even), Facts::NONNEGATIVE),
        ("cot", closure!(|x: &f64| 1.0 / x.tan()), Some(Parity::Odd), Facts::REAL),
        ("sec", closure!(|x: &f64| 1.0 / x.cos()), Some(Parity::Even), Facts::REAL),
        ("csc", closure!(|x: &f64| 1.0 / x.sin()), Some(Parity::Odd), Facts::REAL),
        ("acot", closure!(|x: &f64| (1.0 / x).atan()), Some(Parity::Odd), Facts::REAL),
        ("asec", closure!(|x: &f64| (1.0 / x).acos()), None, Facts::NONE),
        ("acsc", closure!(|x: &f64| (1.0 / x).asin()), Some(Parity::Odd), Facts::NONE),
        ("coth", closure!(|x: &f64| 1.0 / x.tanh()), Some(Parity::Odd), Facts::REAL),
        ("sech", closure!(|x: &f64| 1.0 / x.cosh()), Some(Parity::Even), Facts::POSITIVE),
        ("csch", closure!(|x: &f64| 1.0 / x.sinh()), Some(Parity::Odd), Facts::REAL),
        ("asinh", unary!(asinh), Some(Parity::Odd), Facts::REAL),
        ("acosh", unary!(acosh), None, Facts::NONE),
        ("atanh", unary!(atanh), Some(Parity::Odd), Facts::NONE),
        ("acoth", closure!(|x: &f64| (1.0 / x).atanh()), Some(Parity::Odd), Facts::NONE),
        ("asech", closure!(|x: &f64| (1.0 / x).acosh()), None, Facts::NONE),
        ("acsch", closure!(|x: &f64| (1.0 / x).asinh()), Some(Parity::Odd), Facts::REAL),
    ];
    let mut with_parity = Vec::new();
    for (name, f, parity, on_reals) in functions {
        // `sqrt(x)` and `x^(1/2)` are the same thing; keep the shorter
        // spelling for display.
        // Inverse circular functions cost more than `pi` times a rational, so
        // special values like `acos(0) = pi/2` are preferred.
        let cost = match name {
            | "sqrt" => 1,
            | "asin" | "acos" | "atan" => 5,
            | _ => 3,
        };
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(1)).cost(cost).eval(f))?;
        if let Some(parity) = parity {
            i.graph().ops_mut().set_attr(op, parity);
            with_parity.push(op);
        }
        if on_reals != Facts::NONE {
            i.graph().ops_mut().set_attr(op, OnReals(on_reals));
        }
    }
    let pi = i.op(OpDescriptor::new("pi", Arity::Fixed(0)).eval(|_| std::f64::consts::PI))?;
    i.graph().ops_mut().set_attr(pi, OnReals(Facts::POSITIVE));
    // Cost 3, not 1: `exp(x)` is the preferred spelling of `E^x`.
    let euler = i.op(OpDescriptor::new("E", Arity::Fixed(0)).cost(3).eval(|_| std::f64::consts::E))?;
    i.graph().ops_mut().set_attr(euler, OnReals(Facts::POSITIVE));
    i.op(OpDescriptor::new("atan2", Arity::Fixed(2)).eval(|a| match a {
        | [y, x] => y.atan2(*x),
        | _ => f64::NAN,
    }))?;
    // log(b, x): the logarithm of x to base b.
    let log = i.op(OpDescriptor::new("log", Arity::Fixed(2)).eval(|a| match a {
        | [b, x] => x.ln() / b.ln(),
        | _ => f64::NAN,
    }))?;
    i.kernel("elementary/exact-log", Tier::Normalize, ExactLog { log });
    if let (Some(sin), Some(cos), Some(sinh), Some(cosh)) = (
        i.graph().ops().lookup("sin"),
        i.graph().ops().lookup("cos"),
        i.graph().ops().lookup("sinh"),
        i.graph().ops().lookup("cosh"),
    ) {
        i.kernel("elementary/pythagoras-in-sums", Tier::Normalize, PythagorasInSums { sin, cos, sinh, cosh });
    }
    i.kernel(
        "elementary/parity",
        Tier::Normalize,
        ParityKernel { ops: with_parity },
    );

    i.rewrites(
        Tier::Normalize,
        &[
            "elementary/sqrt: sqrt(?x) => ?x^(1/2)",
            "elementary/exp-0: exp(0) => 1",
            "elementary/ln-1: ln(1) => 0",
            "elementary/exp-ln: exp(ln(?x)) => ?x",
            "elementary/ln-exp: ln(exp(?x)) => ?x",
            "elementary/exp-mul: exp(?a) * exp(?b) => exp(?a + ?b)",
            "elementary/exp-pow: exp(?a) ^ ?n => exp(?n * ?a) if integer(?n)",
            "elementary/exp-n-ln: exp(?n * ln(?x)) => ?x ^ ?n",
            "elementary/ln-pow: ln(?a ^ ?n) => ?n * ln(?a) if positive(?a)",
            "elementary/sin-0: sin(0) => 0",
            "elementary/cos-0: cos(0) => 1",
            "elementary/tan-0: tan(0) => 0",
            "elementary/sin-pi: sin(pi) => 0",
            "elementary/cos-pi: cos(pi) => -1",
            "elementary/tan-pi: tan(pi) => 0",
            "elementary/sin-pi/2: sin(pi/2) => 1",
            "elementary/cos-pi/2: cos(pi/2) => 0",
            "elementary/asin-0: asin(0) => 0",
            "elementary/acos-1: acos(1) => 0",
            "elementary/atan-0: atan(0) => 0",
            "elementary/asin-1: asin(1) => pi/2",
            "elementary/asin-half: asin(1/2) => pi/6",
            "elementary/asin-r2: asin(2^(1/2)/2) => pi/4",
            "elementary/asin-r2b: asin(2^(-1/2)) => pi/4",
            "elementary/asin-r3: asin(3^(1/2)/2) => pi/3",
            "elementary/acos-0: acos(0) => pi/2",
            "elementary/acos-neg: acos(-1) => pi",
            "elementary/acos-half: acos(1/2) => pi/3",
            "elementary/acos-r2: acos(2^(1/2)/2) => pi/4",
            "elementary/acos-r2b: acos(2^(-1/2)) => pi/4",
            "elementary/acos-r3: acos(3^(1/2)/2) => pi/6",
            "elementary/acos-reflect: acos(?x) => pi - acos(-?x) if negative(?x), number(?x)",
            "elementary/atan-1: atan(1) => pi/4",
            "elementary/atan-r3: atan(3^(1/2)) => pi/3",
            "elementary/atan-r3i: atan(3^(-1/2)) => pi/6",
            "elementary/sin-pi/6: sin(pi/6) => 1/2",
            "elementary/cos-pi/3: cos(pi/3) => 1/2",
            "elementary/sin-pi/3: sin(pi/3) => 3^(1/2)/2",
            "elementary/cos-pi/6: cos(pi/6) => 3^(1/2)/2",
            "elementary/sin-pi/4: sin(pi/4) => 2^(-1/2)",
            "elementary/cos-pi/4: cos(pi/4) => 2^(-1/2)",
            "elementary/tan-pi/4: tan(pi/4) => 1",
            "elementary/tan-pi/3: tan(pi/3) => 3^(1/2)",
            "elementary/tan-pi/6: tan(pi/6) => 3^(-1/2)",
            "elementary/sinh-0: sinh(0) => 0",
            "elementary/cosh-0: cosh(0) => 1",
            "elementary/tanh-0: tanh(0) => 0",
            "elementary/abs-abs: abs(abs(?x)) => abs(?x)",
            "elementary/abs-nonneg: abs(?x) => ?x if nonnegative(?x)",
            "elementary/abs-neg: abs(?x) => -?x if negative(?x)",
            "elementary/sqrt-square: (?x ^ 2) ^ (1/2) => abs(?x) if real(?x)",
            "elementary/ln-prod: ln(?a * ?b) => ln(?a) + ln(?b) if positive(?a), positive(?b)",
            "elementary/pythagoras: sin(?x)^2 + cos(?x)^2 => 1",
            "elementary/pythagoras-scaled: ?c * sin(?x)^2 + ?c * cos(?x)^2 => ?c",
            "elementary/hyperbolic-scaled: ?c * cosh(?x)^2 - ?c * sinh(?x)^2 => ?c",
            "elementary/hyperbolic: cosh(?x)^2 - sinh(?x)^2 => 1",
            "elementary/tan: sin(?x) / cos(?x) => tan(?x)",
            "elementary/tanh: sinh(?x) / cosh(?x) => tanh(?x)",
        ],
    )?;
    i.rewrites(
        Tier::Explore,
        &[
            "elementary/sin-double: 2 * sin(?x) * cos(?x) => sin(2 * ?x)",
            "elementary/cos-double: cos(?x)^2 - sin(?x)^2 => cos(2 * ?x)",
            "elementary/one-minus-sin2: 1 - sin(?x)^2 => cos(?x)^2",
            "elementary/one-minus-cos2: 1 - cos(?x)^2 => sin(?x)^2",
            "elementary/ln-mul: ln(?a) + ln(?b) => ln(?a * ?b) if positive(?a), positive(?b)",
            // The reciprocal and inverse-reciprocal functions in terms of
            // the primary ones, so that everything known about those
            // applies; extraction keeps whichever spelling is shorter.
            "elementary/cot: cot(?x) <=> cos(?x) / sin(?x)",
            "elementary/sec: sec(?x) <=> 1 / cos(?x)",
            "elementary/csc: csc(?x) <=> 1 / sin(?x)",
            "elementary/coth: coth(?x) <=> cosh(?x) / sinh(?x)",
            "elementary/sech: sech(?x) <=> 1 / cosh(?x)",
            "elementary/csch: csch(?x) <=> 1 / sinh(?x)",
            "elementary/acot: acot(?x) => atan(1 / ?x)",
            "elementary/asec: asec(?x) => acos(1 / ?x)",
            "elementary/acsc: acsc(?x) => asin(1 / ?x)",
            "elementary/acoth: acoth(?x) => atanh(1 / ?x)",
            "elementary/asech: asech(?x) => acosh(1 / ?x)",
            "elementary/acsch: acsch(?x) => asinh(1 / ?x)",
            "elementary/asinh: asinh(?x) => ln(?x + (?x^2 + 1)^(1/2)) if real(?x)",
            "elementary/log: log(?b, ?x) => ln(?x) / ln(?b)",
            "elementary/tanh-def: tanh(?x) => sinh(?x) / cosh(?x)",
            "elementary/tan-def: tan(?x) => sin(?x) / cos(?x)",
        ],
    )?;
    i.rewrites(
        Tier::Normalize,
        &[
            "elementary/sin-n-pi: sin(?n * pi) => 0 if integer(?n)",
            "elementary/cos-n-pi: cos(?n * pi) => (-1) ^ ?n if integer(?n)",
            "elementary/tan-n-pi: tan(?n * pi) => 0 if integer(?n)",
            "elementary/sin-shift-n-pi: sin(?a + pi * ?n) => (-1) ^ ?n * sin(?a) if integer(?n)",
            "elementary/cos-shift-n-pi: cos(?a + pi * ?n) => (-1) ^ ?n * cos(?a) if integer(?n)",
            "elementary/tan-shift-n-pi: tan(?a + pi * ?n) => tan(?a) if integer(?n)",
            "elementary/sin-pi-n: sin(pi * ?n) => 0 if integer(?n)",
            "elementary/cos-pi-n: cos(pi * ?n) => (-1) ^ ?n if integer(?n)",
            "elementary/cot-tan: cot(?x) * tan(?x) => 1",
            "elementary/coth-tanh: coth(?x) * tanh(?x) => 1",
            "elementary/E-pow: E ^ ?x => exp(?x)",
            "elementary/exp-1: exp(1) => E",
            "elementary/ln-E: ln(E) => 1",
            "elementary/log-same: log(?b, ?b) => 1 if positive(?b)",
            "elementary/sinh-asinh: sinh(asinh(?x)) => ?x",
            "elementary/tanh-atanh: tanh(atanh(?x)) => ?x",
            "elementary/sin-asin: sin(asin(?x)) => ?x",
            "elementary/cos-acos: cos(acos(?x)) => ?x",
            "elementary/tan-atan: tan(atan(?x)) => ?x",
            "elementary/cosh-acosh: cosh(acosh(?x)) => ?x",
            "elementary/atan2-right: atan2(?y, ?x) => atan(?y / ?x) if positive(?x)",
            "elementary/atan2-up: atan2(?y, 0) => pi / 2 if positive(?y)",
            "elementary/atan2-down: atan2(?y, 0) => -pi / 2 if negative(?y)",
            "elementary/atan2-left-up: atan2(?y, ?x) => atan(?y / ?x) + pi if negative(?x), nonnegative(?y)",
            "elementary/atan2-left-down: atan2(?y, ?x) => atan(?y / ?x) - pi if negative(?x), negative(?y)",
        ],
    )
}

/// `log(b, x)` for literal integers where `x` is an exact power of `b`.
/// `c X sin(u)^2 + c X cos(u)^2 = c X` and `c X cosh(u)^2 - c X sinh(u)^2
/// = c X` between any two terms of a sum, whatever else the sum holds and
/// whatever common factor `X` the two terms share. (The rewrite rules see
/// only sums whose other terms they can absorb into a pattern variable.)
struct PythagorasInSums {
    sin: OpId,
    cos: OpId,
    sinh: OpId,
    cosh: OpId,
}

/// A term of a sum as `coefficient * factors`.
fn split_term(
    graph: &crate::graph::Graph,
    term: NodeId,
) -> (crate::graph::Number, Vec<NodeId>) {
    let factors = if graph.op(term) == core::MUL { graph.children(term).to_vec() } else { vec![term] };
    let mut coefficient = crate::graph::Number::from(1);
    let mut rest = Vec::with_capacity(factors.len());
    for f in factors {
        match graph.number_of(f) {
            | Some(n) => coefficient = coefficient.mul(n),
            | None => rest.push(f),
        }
    }
    (coefficient, rest)
}

impl PythagorasInSums {
    /// `(function, argument)` of `f(u)^2`.
    fn square_of(
        &self,
        graph: &crate::graph::Graph,
        factor: NodeId,
    ) -> Option<(OpId, NodeId)> {
        let &[base, exponent] = graph.children(factor) else {
            return None;
        };
        if graph.op(factor) != core::POW || graph.number_of(exponent).and_then(crate::graph::Number::to_i64) != Some(2) {
            return None;
        }
        let op = graph.op(base);
        let &[u] = graph.children(base) else {
            return None;
        };
        [self.sin, self.cos, self.sinh, self.cosh].contains(&op).then_some((op, u))
    }
}

impl Kernel for PythagorasInSums {
    fn ops(&self) -> Vec<OpId> {
        vec![core::ADD]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let terms = graph.children(node).to_vec();
        if terms.len() < 2 || terms.len() > 64 {
            return Outcome::Pass;
        }
        // (term index, function, argument, coefficient, other factors)
        let mut squares = Vec::new();
        for (k, &term) in terms.iter().enumerate() {
            let (coefficient, factors) = split_term(graph, term);
            for (j, &f) in factors.iter().enumerate() {
                if let Some((op, u)) = self.square_of(graph, f) {
                    let mut others: Vec<NodeId> = factors.iter().enumerate().filter(|&(i, _)| i != j).map(|(_, &g)| g).collect();
                    others.sort_unstable();
                    squares.push((k, op, u, coefficient.clone(), others));
                }
            }
        }
        for a in &squares {
            for b in &squares {
                if a.0 == b.0 || !graph.same(a.2, b.2) || a.4 != b.4 {
                    continue;
                }
                let matches = (a.1 == self.sin && b.1 == self.cos && a.3 == b.3)
                    || (a.1 == self.cosh && b.1 == self.sinh && a.3 == b.3.neg());
                if !matches {
                    continue;
                }
                let coefficient = graph.num(a.3.clone());
                let mut factors = vec![coefficient];
                factors.extend(&a.4);
                let replacement = if factors.len() == 1 { factors[0] } else { graph.node(core::MUL, &factors) };
                let mut rest: Vec<NodeId> =
                    terms.iter().enumerate().filter(|&(i, _)| i != a.0 && i != b.0).map(|(_, &t)| t).collect();
                rest.push(replacement);
                let sum = if rest.len() == 1 { rest[0] } else { graph.node(core::ADD, &rest) };
                return Outcome::Equal(sum);
            }
        }
        Outcome::Pass
    }

    fn revisit(&self) -> bool {
        true
    }
}

struct ExactLog {
    log: OpId,
}

impl Kernel for ExactLog {
    fn ops(&self) -> Vec<OpId> {
        vec![self.log]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[base, argument] = graph.children(node) else {
            return Outcome::Pass;
        };
        let (Some(b), Some(x)) = (
            graph.number_of(base).and_then(crate::graph::Number::to_i64),
            graph.number_of(argument).and_then(crate::graph::Number::to_i64),
        ) else {
            return Outcome::Pass;
        };
        if b < 2 || x < 1 {
            return Outcome::Pass;
        }
        let (mut power, mut exponent) = (1_i64, 0_i64);
        while power < x {
            let Some(next) = power.checked_mul(b) else {
                return Outcome::Pass;
            };
            power = next;
            exponent += 1;
        }
        if power == x { Outcome::Equal(graph.int(exponent)) } else { Outcome::Pass }
    }
}

/// Moves a minus sign out of odd functions and drops it inside even ones.
struct ParityKernel {
    ops: Vec<OpId>,
}

impl Kernel for ParityKernel {
    fn ops(&self) -> Vec<OpId> {
        self.ops.clone()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        let (Some(&parity), Some(&arg)) =
            (graph.ops().attr::<Parity>(op), graph.children(node).first())
        else {
            return Outcome::Pass;
        };
        // The argument is negative if it is a negative number or a product
        // with a negative literal coefficient.
        let positive_arg =
            if let Some(n) = graph.as_number(arg).filter(|n| n.is_negative()).cloned() {
                graph.num(n.neg())
            } else if graph.op(arg) == core::MUL {
                let factors = graph.children(arg).to_vec();
                let Some(index) = factors
                    .iter()
                    .position(|&f| graph.as_number(f).is_some_and(super::super::graph::number::Number::is_negative))
                else {
                    return Outcome::Pass;
                };
                let mut flipped = factors.clone();
                let negated = graph
                    .as_number(factors[index])
                    .map_or_else(|| 0.into(), super::super::graph::number::Number::neg);
                flipped[index] = graph.num(negated);
                graph.node(core::MUL, &flipped)
            } else {
                return Outcome::Pass;
            };
        let inner = graph.node(op, &[positive_arg]);
        match parity {
            | Parity::Even => Outcome::Equal(inner),
            | Parity::Odd => {
                let minus_one = graph.int(-1);
                Outcome::Equal(graph.node(core::MUL, &[minus_one, inner]))
            },
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Evaluated;
    use crate::graph::Extractor;
    use crate::graph::Graph;
    use crate::graph::Saturate;

    fn simplify(src: &str) -> String {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[elementary()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(g.conflicts().is_empty(), "{:?}", g.conflicts());
        let ex = Extractor::new(&g, &[root], &ClosedForm);
        ex.build(&mut g, root)
            .map_or_else(|| "<none>".to_owned(), |n| g.display(n))
    }

    fn evaluate(src: &str) -> f64 {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[elementary()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::numeric(1e-12),
            &Evaluated,
            &Budget::default(),
        );
        g.approx(g.find(root)).map_or(f64::NAN, |b| b.mid)
    }

    #[test]
    fn exact_special_values() {
        assert_eq!(simplify("sin(0) + cos(0)"), "1");
        assert_eq!(simplify("sin(pi) + cos(pi)"), "-1");
        assert_eq!(simplify("sin(pi/2)"), "1");
        assert_eq!(simplify("exp(0) * ln(1)"), "0");
        assert_eq!(
            simplify("sin(1)"),
            "sin(1)",
            "exact arguments stay symbolic"
        );
    }

    #[test]
    fn inverse_pairs_cancel() {
        assert_eq!(simplify("exp(ln(x))"), "x");
        assert_eq!(simplify("ln(exp(x + 1))"), "x + 1");
        assert_eq!(simplify("sqrt(x)^2"), "x");
        assert_eq!(simplify("sqrt(4)"), "2");
        assert_eq!(simplify("sqrt(2)"), "sqrt(2)");
    }

    #[test]
    fn pythagorean_identities() {
        assert_eq!(simplify("sin(x)^2 + cos(x)^2"), "1");
        assert_eq!(simplify("cos(a + b)^2 + 3 + sin(b + a)^2"), "4");
        assert_eq!(simplify("cosh(t)^2 - sinh(t)^2"), "1");
        assert_eq!(simplify("1 - sin(x)^2"), "cos(x)^2");
    }

    #[test]
    fn parity() {
        assert_eq!(simplify("sin(-x)"), "-sin(x)");
        assert_eq!(simplify("cos(-x)"), "cos(x)");
        assert_eq!(simplify("sin(-2*x) + sin(2*x)"), "0");
        assert_eq!(simplify("cos(-3*x*y) - cos(3*y*x)"), "0");
        assert_eq!(simplify("abs(-x)"), "abs(x)");
        assert_eq!(simplify("tan(-1) + tan(1)"), "0");
    }

    #[test]
    fn exponentials_combine() {
        assert_eq!(simplify("exp(x) * exp(y)"), "exp(x + y)");
        assert_eq!(simplify("exp(x) * exp(-x)"), "1");
        assert_eq!(simplify("exp(x)^2"), "exp(2*x)");
        assert_eq!(simplify("ln(2^x)"), "x*ln(2)");
        assert_eq!(
            simplify("ln(y^x)"),
            "ln(y^x)",
            "y is not known to be positive"
        );
    }

    #[test]
    fn assumptions_unlock_conditional_identities() {
        let run = |src: &str, assume: &[(&str, Facts)]| {
            let mut g = Graph::new();
            let engine = Engine::install(&mut g, &[elementary()]).unwrap_or_else(|e| panic!("{e}"));
            for (name, facts) in assume {
                let symbol = g.interner_mut().symbol(name);
                g.assume(symbol, *facts);
            }
            let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
            engine.run(
                &mut g,
                &[root],
                &Env::symbolic(),
                &Saturate,
                &Budget::default(),
            );
            assert!(g.conflicts().is_empty());
            let ex = Extractor::new(&g, &[root], &ClosedForm);
            ex.build(&mut g, root)
                .map_or_else(|| "<none>".to_owned(), |n| g.display(n))
        };
        assert_eq!(run("sqrt(x^2)", &[]), "sqrt(x^2)");
        assert_eq!(run("sqrt(x^2)", &[("x", Facts::REAL)]), "abs(x)");
        assert_eq!(run("sqrt(x^2)", &[("x", Facts::POSITIVE)]), "x");
        assert_eq!(run("sqrt(x^2)", &[("x", Facts::NEGATIVE)]), "-x");
        assert_eq!(run("abs(exp(t))", &[("t", Facts::REAL)]), "exp(t)");
        assert_eq!(run("abs(exp(t))", &[]), "abs(exp(t))");
        assert_eq!(run("ln(a^3)", &[("a", Facts::POSITIVE)]), "3*ln(a)");
        assert_eq!(run("abs(t^2 + 1)", &[("t", Facts::REAL)]), "t^2 + 1");
        assert_eq!(run("abs(pi)", &[]), "pi");
    }

    #[test]
    fn reciprocal_and_inverse_families() {
        assert_eq!(simplify("sec(x) * cos(x)"), "1");
        assert_eq!(simplify("csc(x) * sin(x)"), "1");
        assert_eq!(simplify("cot(x) * tan(x)"), "1");
        assert_eq!(simplify("sech(x) * cosh(x)"), "1");
        assert_eq!(simplify("sin(asin(x)) + tan(atan(y))"), "x + y");
        assert_eq!(simplify("ln(E)"), "1");
        assert_eq!(simplify("E^x"), "exp(x)");
        assert_eq!(simplify("log(2, 8)"), "3");
        assert_eq!(simplify("log(10, 1000) + log(3, 1)"), "3");
        assert_eq!(simplify("log(2, 10)"), "log(2, 10)");
        assert!((evaluate("log(2, 8)") - 3.0).abs() < 1e-14);
        assert!((evaluate("atan2(1, 1)") - std::f64::consts::FRAC_PI_4).abs() < 1e-15);
        for (src, want) in [
            ("cot(1)", 1.0 / 1.0_f64.tan()),
            ("sec(1)", 1.0 / 1.0_f64.cos()),
            ("csc(1)", 1.0 / 1.0_f64.sin()),
            ("acot(2)", 0.5_f64.atan()),
            ("asec(2)", 0.5_f64.acos()),
            ("acsc(2)", 0.5_f64.asin()),
            ("coth(1)", 1.0 / 1.0_f64.tanh()),
            ("sech(1)", 1.0 / 1.0_f64.cosh()),
            ("csch(1)", 1.0 / 1.0_f64.sinh()),
            ("asinh(1)", 1.0_f64.asinh()),
            ("acosh(2)", 2.0_f64.acosh()),
            ("atanh(1/2)", 0.5_f64.atanh()),
            ("acoth(2)", 0.5_f64.atanh()),
            ("asech(1/2)", 2.0_f64.acosh()),
            ("acsch(2)", 0.5_f64.asinh()),
            ("E", std::f64::consts::E),
        ] {
            assert!((evaluate(src) - want).abs() < 1e-14, "{src}");
        }
    }

    #[test]
    fn quotients_become_tangents() {
        assert_eq!(simplify("sin(x) / cos(x)"), "tan(x)");
        assert_eq!(simplify("2 * sin(x) * cos(x)"), "sin(2*x)");
    }

    #[test]
    fn numeric_phase() {
        assert!((evaluate("sin(pi/6)") - 0.5).abs() < 1e-15);
        assert!((evaluate("exp(1)") - std::f64::consts::E).abs() < 1e-15);
        assert!((evaluate("sqrt(2)^2") - 2.0).abs() < 1e-15);
        assert_eq!(simplify("sin(0.5)"), format!("{}", 0.5_f64.sin()));
    }
}
