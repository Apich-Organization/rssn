//! Elementary functions: exponentials, logarithms, trigonometric and
//! hyperbolic functions, together with the identities that never enlarge a
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
    let functions: [(&str, EvalFn, Option<Parity>, Facts); 13] = [
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
    ];
    let mut with_parity = Vec::new();
    for (name, f, parity, on_reals) in functions {
        // `sqrt(x)` and `x^(1/2)` are the same thing; keep the shorter
        // spelling for display.
        let cost = if name == "sqrt" { 1 } else { 3 };
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
            "elementary/sinh-0: sinh(0) => 0",
            "elementary/cosh-0: cosh(0) => 1",
            "elementary/tanh-0: tanh(0) => 0",
            "elementary/abs-abs: abs(abs(?x)) => abs(?x)",
            "elementary/abs-nonneg: abs(?x) => ?x if nonnegative(?x)",
            "elementary/abs-neg: abs(?x) => -?x if negative(?x)",
            "elementary/sqrt-square: (?x ^ 2) ^ (1/2) => abs(?x) if real(?x)",
            "elementary/ln-prod: ln(?a * ?b) => ln(?a) + ln(?b) if positive(?a), positive(?b)",
            "elementary/pythagoras: sin(?x)^2 + cos(?x)^2 => 1",
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
        ],
    )
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
                    .position(|&f| graph.as_number(f).is_some_and(|n| n.is_negative()))
                else {
                    return Outcome::Pass;
                };
                let mut flipped = factors.clone();
                let negated = graph
                    .as_number(factors[index])
                    .map_or_else(|| 0.into(), |n| n.neg());
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
