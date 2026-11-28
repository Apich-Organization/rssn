//! Branches of multi-valued functions.
//!
//! Each operator takes a branch index `k` (any integer; `k = 0` is the
//! principal branch) and is defined by rewriting to principal functions:
//!
//! | operator | value |
//! |---|---|
//! | `log_branch(z, k)` | `ln\|z\| + I (arg z + 2 pi k)` |
//! | `sqrt_branch(z, k)` | `\|z\|^(1/2) exp(I (arg z + 2 pi k) / 2)` |
//! | `root_branch(z, n, k)` | `\|z\|^(1/n) exp(I (arg z + 2 pi k) / n)` |
//! | `power_branch(z, w, k)` | `exp(w log_branch(z, k))` |
//! | `asin_branch(z, k)` | `k pi + (-1)^k asin z` |
//! | `acos_branch(z, k, s)` | `2 k pi + s acos z`, `s = ±1` |
//! | `atan_branch(z, k)` | `k pi + atan z` |
//! | `asinh_branch(z, k)`, `acosh_branch(z, k)`, `atanh_branch(z, k)` | logarithmic forms on branch `k` |
//!
//! Positive integer root degrees are periodic in the integer branch index.
//! Exact symbolic indices are reduced before constructing the phase.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::Signed;

use crate::graph::Arity;
use crate::graph::ComplexEval;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::complex_coords;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Branch `k` of the `n`-th root.
#[allow(
    clippy::float_cmp,
    reason = "the exact orders one and two have dedicated formulas"
)]
fn root(
    z: Complex64,
    n: f64,
    k: f64,
) -> Complex64 {
    let integral_index = k.fract() == 0.0;
    if n == 1.0 && integral_index {
        return z;
    }
    let scale = z.re.abs().max(z.im.abs());
    let scaled = ((0.0 < scale && scale < f64::MIN_POSITIVE)
        || (scale.is_finite() && scale > f64::MAX / 2.0))
        .then(|| z / scale);
    if n == 2.0 && integral_index {
        let principal = scaled.map_or_else(|| z.sqrt(), |unit| unit.sqrt() * scale.sqrt());
        return if k % 2.0 == 0.0 {
            principal
        } else {
            -principal
        };
    }
    let phase = if n > 0.0 && n.fract() == 0.0 && integral_index {
        let reduced = if k.abs() >= n { k % n } else { k };
        z.arg() / n + TAU * (reduced / n)
    } else {
        (z.arg() + TAU * k) / n
    };
    let inverse_degree = n.recip();
    let radius = scaled.map_or_else(
        || z.norm().powf(inverse_degree),
        |unit| scale.powf(inverse_degree) * unit.norm().powf(inverse_degree),
    );
    Complex64::from_polar(radius, phase)
}

struct RootBranches {
    root: OpId,
    sqrt: OpId,
    principal_sqrt: OpId,
    abs: OpId,
    arg: OpId,
    exp: OpId,
    pi: OpId,
    unit: OpId,
}

fn scaled_index(
    graph: &mut Graph,
    degree: NodeId,
    index: NodeId,
    inverse: NodeId,
) -> NodeId {
    let (Some(Number::Int(degree)), Some(Number::Int(value))) =
        (graph.number_of(degree), graph.number_of(index))
    else {
        return graph.node(core::MUL, &[index, inverse]);
    };
    if !degree.is_positive() {
        return graph.node(core::MUL, &[index, inverse]);
    }
    let reduced = if value.magnitude() >= degree.magnitude() {
        value % degree
    } else {
        value.clone()
    };
    graph.num(Number::rat(BigRational::new(reduced, degree.clone())))
}

impl Kernel for RootBranches {
    fn ops(&self) -> Vec<OpId> {
        vec![self.root, self.sqrt]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let (z, degree, index) = match *graph.children(node) {
            | [z, index] => (z, None, index),
            | [z, degree, index] => (z, Some(degree), index),
            | _ => return Outcome::Pass,
        };
        let two = graph.int(2);
        let degree = degree.unwrap_or(two);
        if let Some(Number::Int(value)) = graph.number_of(index) {
            match graph.number_of(degree).and_then(Number::to_i64) {
                | Some(1) => return Outcome::Equal(z),
                | Some(2) => {
                    let odd = value.is_odd();
                    let principal = graph.node(self.principal_sqrt, &[z]);
                    let result = if odd {
                        let minus_one = graph.int(-1);
                        graph.node(core::MUL, &[minus_one, principal])
                    } else {
                        principal
                    };
                    return Outcome::Equal(result);
                },
                | _ => {},
            }
        }
        let minus_one = graph.int(-1);
        let inverse = graph.node(core::POW, &[degree, minus_one]);
        let branch_fraction = scaled_index(graph, degree, index, inverse);
        let magnitude = graph.node(self.abs, &[z]);
        let radius = graph.node(core::POW, &[magnitude, inverse]);
        let unit = graph.node(self.unit, &[]);
        let pi = graph.node(self.pi, &[]);
        let argument = graph.node(self.arg, &[z]);
        let fractional_argument = graph.node(core::MUL, &[argument, inverse]);
        let shift = graph.node(core::MUL, &[two, pi, branch_fraction]);
        let angle = graph.node(core::ADD, &[fractional_argument, shift]);
        let phase = graph.node(core::MUL, &[unit, angle]);
        let rotation = graph.node(self.exp, &[phase]);
        Outcome::Equal(graph.node(core::MUL, &[radius, rotation]))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, doc) in [
        (
            "log_branch",
            2,
            "`log_branch(z, k)`: branch `k` of the complex logarithm `ln|z| + I * (arg(z) + 2*pi*k)`.",
        ),
        (
            "sqrt_branch",
            2,
            "`sqrt_branch(z, k)` selects the principal square root for an even integer `k` and its negative for an odd integer `k`.",
        ),
        (
            "root_branch",
            3,
            "`root_branch(z, n, k)` gives `|z|^(1/n) * exp(I * (arg(z) + 2*pi*k) / n)`. For a positive integer degree `n`, integer indices that differ by a multiple of `n` select the same root.",
        ),
        (
            "power_branch",
            3,
            "`power_branch(z, w, k)`: branch `k` of the complex power `exp(w * log_branch(z, k))`.",
        ),
        (
            "asin_branch",
            2,
            "`asin_branch(z, k)`: branch `k` of the inverse sine `k*pi + (-1)^k * asin(z)`.",
        ),
        (
            "acos_branch",
            3,
            "`acos_branch(z, k, s)`: branch `(k, s)` of the inverse cosine `2*k*pi + s * acos(z)` with sign `s` in `{+1, -1}`.",
        ),
        (
            "atan_branch",
            2,
            "`atan_branch(z, k)`: branch `k` of the inverse tangent `k*pi + atan(z)`.",
        ),
        (
            "asinh_branch",
            2,
            "`asinh_branch(z, k)`: branch `k` of the inverse hyperbolic sine.",
        ),
        (
            "acosh_branch",
            2,
            "`acosh_branch(z, k)`: branch `k` of the inverse hyperbolic cosine.",
        ),
        (
            "atanh_branch",
            2,
            "`atanh_branch(z, k)`: branch `k` of the inverse hyperbolic tangent.",
        ),
    ] {
        // Requests: the definitions below must win extraction.
        i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc(doc))?;
    }
    // Independent numeric definitions, so that the rewrites below are
    // checked against something other than themselves.
    let evals: [(&str, ComplexEval); 10] = [
        (
            "log_branch",
            ComplexEval(|a| Some(a.first()?.ln() + I * TAU * *a.get(1)?)),
        ),
        (
            "sqrt_branch",
            ComplexEval(|a| Some(root(*a.first()?, 2.0, a.get(1)?.re))),
        ),
        (
            "root_branch",
            ComplexEval(|a| Some(root(*a.first()?, a.get(1)?.re, a.get(2)?.re))),
        ),
        (
            "power_branch",
            ComplexEval(|a| Some((*a.get(1)? * (a.first()?.ln() + I * TAU * *a.get(2)?)).exp())),
        ),
        (
            "asin_branch",
            ComplexEval(|a| {
                let k = a.get(1)?.re;
                Some(k * PI + (-1_f64).powf(k) * a.first()?.asin())
            }),
        ),
        (
            "acos_branch",
            ComplexEval(|a| Some(2.0 * a.get(1)?.re * PI + *a.get(2)? * a.first()?.acos())),
        ),
        (
            "atan_branch",
            ComplexEval(|a| Some(a.get(1)?.re * PI + a.first()?.atan())),
        ),
        (
            "asinh_branch",
            ComplexEval(|a| {
                let k = a.get(1)?.re;
                Some(I * k * PI + (-1_f64).powf(k) * a.first()?.asinh())
            }),
        ),
        (
            "acosh_branch",
            ComplexEval(|a| {
                let z = *a.first()?;
                let w = z * z - 1.0;
                let w = if w.im == 0.0 {
                    Complex64::new(w.re, 0.0)
                } else {
                    w
                };
                Some((z + w.sqrt()).ln() + I * TAU * *a.get(1)?)
            }),
        ),
        (
            "atanh_branch",
            ComplexEval(|a| {
                // The logarithmic form fixes the side of the cuts (num-complex's
                // atanh takes the other one for real |z| > 1).
                let z = *a.first()?;
                let w = (1.0 + z) / (1.0 - z);
                let w = if w.im == 0.0 {
                    Complex64::new(w.re, 0.0)
                } else {
                    w
                };
                Some(0.5 * w.ln() + I * PI * *a.get(1)?)
            }),
        ),
    ];
    for (name, eval) in evals {
        if let Some(op) = i.graph().ops().lookup(name) {
            i.graph().ops_mut().set_attr(op, eval);
        }
    }
    let [
        Some(root),
        Some(sqrt),
        Some(principal_sqrt),
        Some(abs),
        Some(arg),
        Some(exp),
        Some(pi),
        Some(unit),
    ] = [
        "root_branch",
        "sqrt_branch",
        "sqrt",
        "abs",
        "arg",
        "exp",
        "pi",
        "I",
    ]
    .map(|name| i.graph().ops().lookup(name))
    else {
        return Err(RuleError::Invalid {
            rule: "complex/periodic-root".into(),
            reason: "needs elementary and complex operators",
        });
    };
    i.kernel(
        "complex/periodic-root",
        Tier::Normalize,
        RootBranches {
            root,
            sqrt,
            principal_sqrt,
            abs,
            arg,
            exp,
            pi,
            unit,
        },
    );
    for (name, arity, doc) in [
        (
            "arg_branch",
            2,
            "`arg_branch(z, theta0)`: argument of `z` in `(theta0, theta0 + 2*pi]`, placing the branch cut along ray `theta0`.",
        ),
        (
            "log_cut",
            2,
            "`log_cut(z, theta0)`: complex logarithm `ln|z| + I*arg_branch(z, theta0)` with branch cut along ray `theta0`.",
        ),
        (
            "sqrt_cut",
            2,
            "`sqrt_cut(z, theta0)`: square root with branch cut along ray `theta0`.",
        ),
        (
            "pow_cut",
            3,
            "`pow_cut(z, w, theta0)`: complex power `exp(w * log_cut(z, theta0))` with branch cut along ray `theta0`.",
        ),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .cost(10)
            .doc(doc))?;
        let eval = match name {
            | "arg_branch" => {
                ComplexEval(|a| {
                    Some(Complex64::new(
                        complex_coords::arg_branch(*a.first()?, a.get(1)?.re),
                        0.0,
                    ))
                })
            },
            | "log_cut" => {
                ComplexEval(|a| Some(complex_coords::log_cut(*a.first()?, a.get(1)?.re)))
            },
            | "sqrt_cut" => {
                ComplexEval(|a| Some(complex_coords::sqrt_cut(*a.first()?, a.get(1)?.re)))
            },
            | _ => {
                ComplexEval(|a| {
                    Some(complex_coords::pow_cut(
                        *a.first()?,
                        *a.get(1)?,
                        a.get(2)?.re,
                    ))
                })
            },
        };
        i.graph().ops_mut().set_attr(op, eval);
    }
    i.rewrites(
        Tier::Normalize,
        &[
            "branch/log-cut: log_cut(?z, ?theta0) => ln(abs(?z)) + I * arg_branch(?z, ?theta0)",
            "branch/sqrt-cut: sqrt_cut(?z, ?theta0) => abs(?z)^(1/2) * exp(I * arg_branch(?z, ?theta0) / 2)",
            "branch/pow-cut: pow_cut(?z, ?w, ?theta0) => exp(?w * log_cut(?z, ?theta0))",
            "branch/arg-principal: arg_branch(?z, -pi) => arg(?z)",
        ],
    )?;
    i.rewrites(
        Tier::Normalize,
        &[
            "branch/log: log_branch(?z, ?k) => ln(abs(?z)) + I * (arg(?z) + 2 * pi * ?k)",
            "branch/power: power_branch(?z, ?w, ?k) => exp(?w * log_branch(?z, ?k))",
            "branch/asin: asin_branch(?z, ?k) => ?k * pi + (-1)^?k * asin(?z) if integer(?k)",
            "branch/acos: acos_branch(?z, ?k, ?s) => 2 * ?k * pi + ?s * acos(?z)",
            "branch/atan: atan_branch(?z, ?k) => ?k * pi + atan(?z)",
            "branch/asinh: asinh_branch(?z, ?k) => I * ?k * pi + (-1)^?k * asinh(?z) if integer(?k)",
            "branch/acosh: acosh_branch(?z, ?k) => log_branch(?z + (?z^2 - 1)^(1/2), ?k)",
            "branch/atanh: atanh_branch(?z, ?k) => log_branch((1 + ?z) / (1 - ?z), ?k) / 2",
        ],
    )
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use num_complex::Complex64;

    use super::super::complex;
    use super::super::eval_complex;
    use super::root;
    use crate::graph::Budget;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Graph;
    use crate::graph::Saturate;
    use crate::graph::SizeCost;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[complex()], src)
    }

    #[test]
    fn principal_and_other_branches() {
        assert_eq!(run("log_branch(1, 0)"), "0");
        assert_eq!(run("log_branch(1, 1)"), "2*I*pi");
        assert_eq!(run("log_branch(-1, 0)"), "I*pi");
        assert_eq!(run("sqrt_branch(4, 0)"), "2");
        assert_eq!(run("sqrt_branch(4, 1)"), "-2");
        assert_eq!(run("root_branch(8, 3, 0)"), "2");
        assert_eq!(run("asin_branch(0, 1)"), "pi");
        assert_eq!(run("atan_branch(1, 1)"), "5/4*pi");
        assert_eq!(run("acos_branch(1, 1, -1)"), "2*pi");
    }

    #[test]
    fn floating_root_phase_respects_large_integer_periods() {
        let radicand = Complex64::new(8.0, 0.0);
        let index = 9_007_199_254_740_991.0;
        let radius = 3.0_f64.sqrt();
        assert!((root(radicand, 3.0, index) - Complex64::new(-1.0, radius)).norm() < 1e-14);
        assert!((root(radicand, 3.0, -index) - Complex64::new(-1.0, -radius)).norm() < 1e-14);
        assert_eq!(root(radicand, 1.0, index), radicand);
        assert!((root(radicand, 1e308, 5e307) + 1.0).norm() < 1e-14);
    }

    #[test]
    fn square_root_scaling_preserves_the_cut_and_extreme_magnitudes() {
        let upper = Complex64::new(-4.0, 0.0);
        let lower = Complex64::new(-4.0, -0.0);
        assert_eq!(root(upper, 2.0, 0.0), Complex64::new(0.0, 2.0));
        assert_eq!(root(lower, 2.0, 0.0), Complex64::new(0.0, -2.0));
        assert_eq!(root(lower, 2.0, -1.0), Complex64::new(0.0, 2.0));

        let large = 1.6e308_f64;
        let unit = root(Complex64::new(large, large), 2.0, 0.0) / large.sqrt();
        assert!((unit * unit - Complex64::new(1.0, 1.0)).norm() < 1e-14);

        let tiny = f64::from_bits(1);
        let unit = root(Complex64::new(0.0, tiny), 2.0, 0.0) / tiny.sqrt();
        assert!((unit * unit - Complex64::new(0.0, 1.0)).norm() < 1e-14);
    }

    /// Every branch, exponentiated or mapped back, gives the argument.
    #[test]
    fn branches_invert() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[complex()]).unwrap_or_else(|e| panic!("{e}"));
        let z = g.interner_mut().symbol("z");
        let point = Complex64::new(-0.7, 0.4);
        let bindings: HashMap<_, _> = std::iter::once((z, point)).collect();
        for k in -2..=2 {
            let cases = [
                (format!("exp(log_branch(z, {k}))"), point),
                (format!("sqrt_branch(z, {k})^2"), point),
                (format!("root_branch(z, 3, {k})^3"), point),
                (format!("sin(asin_branch(z, {k}))"), point),
                (format!("cos(acos_branch(z, {k}, -1))"), point),
                (format!("tan(atan_branch(z, {k}))"), point),
                (format!("sinh(asinh_branch(z, {k}))"), point),
                (format!("cosh(acosh_branch(z, {k}))"), point),
                (format!("tanh(atanh_branch(z, {k}))"), point),
            ];
            for (src, want) in cases {
                let root = g.parse(&src).unwrap_or_else(|e| panic!("{e}"));
                engine.run(
                    &mut g,
                    &[root],
                    &Env::symbolic(),
                    &Saturate,
                    &Budget::default(),
                );
                let term = Extractor::new(&g, &[root], &SizeCost)
                    .build(&mut g, root)
                    .unwrap_or(root);
                let got = eval_complex(&g, term, &bindings)
                    .unwrap_or_else(|| panic!("cannot evaluate {}", g.display(term)));
                assert!((got - want).norm() < 1e-9, "{src}: {got} vs {want}");
            }
        }
    }
}
