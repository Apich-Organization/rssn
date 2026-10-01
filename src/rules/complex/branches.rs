//! Branches of multi-valued functions.
//!
//! Each operator takes a branch index `k` (any integer; `k = 0` is the
//! principal branch) and is defined by rewriting to principal functions:
//!
//! | operator | value |
//! |---|---|
//! | `log_branch(z, k)` | `ln|z| + I (arg z + 2 pi k)` |
//! | `sqrt_branch(z, k)` | `|z|^(1/2) exp(I (arg z + 2 pi k) / 2)` |
//! | `root_branch(z, n, k)` | `|z|^(1/n) exp(I (arg z + 2 pi k) / n)` |
//! | `power_branch(z, w, k)` | `exp(w log_branch(z, k))` |
//! | `asin_branch(z, k)` | `k pi + (-1)^k asin z` |
//! | `acos_branch(z, k, s)` | `2 k pi + s acos z`, `s = ±1` |
//! | `atan_branch(z, k)` | `k pi + atan z` |
//! | `asinh_branch(z, k)`, `acosh_branch(z, k)`, `atanh_branch(z, k)` | logarithmic forms on branch `k` |

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

use crate::graph::Arity;
use crate::graph::ComplexEval;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::rule::Installer;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Branch `k` of the `n`-th root.
fn root(
    z: Complex64,
    n: f64,
    k: f64,
) -> Complex64 {
    Complex64::from_polar(z.norm().powf(1.0 / n), (z.arg() + TAU * k) / n)
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity) in [
        ("log_branch", 2),
        ("sqrt_branch", 2),
        ("root_branch", 3),
        ("power_branch", 3),
        ("asin_branch", 2),
        ("acos_branch", 3),
        ("atan_branch", 2),
        ("asinh_branch", 2),
        ("acosh_branch", 2),
        ("atanh_branch", 2),
    ] {
        // Requests: the definitions below must win extraction.
        i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
    }
    // Independent numeric definitions, so that the rewrites below are
    // checked against something other than themselves.
    let evals: [(&str, ComplexEval); 10] = [
        ("log_branch", ComplexEval(|a| Some(a.first()?.ln() + I * TAU * *a.get(1)?))),
        ("sqrt_branch", ComplexEval(|a| Some(root(*a.first()?, 2.0, a.get(1)?.re)))),
        ("root_branch", ComplexEval(|a| Some(root(*a.first()?, a.get(1)?.re, a.get(2)?.re)))),
        ("power_branch", ComplexEval(|a| Some((*a.get(1)? * (a.first()?.ln() + I * TAU * *a.get(2)?)).exp()))),
        ("asin_branch", ComplexEval(|a| {
            let k = a.get(1)?.re;
            Some(k * PI + (-1_f64).powf(k) * a.first()?.asin())
        })),
        ("acos_branch", ComplexEval(|a| Some(2.0 * a.get(1)?.re * PI + *a.get(2)? * a.first()?.acos()))),
        ("atan_branch", ComplexEval(|a| Some(a.get(1)?.re * PI + a.first()?.atan()))),
        ("asinh_branch", ComplexEval(|a| Some(a.first()?.asinh() + I * TAU * *a.get(1)?))),
        ("acosh_branch", ComplexEval(|a| {
            let z = *a.first()?;
            let w = z * z - 1.0;
            let w = if w.im == 0.0 { Complex64::new(w.re, 0.0) } else { w };
            Some((z + w.sqrt()).ln() + I * TAU * *a.get(1)?)
        })),
        ("atanh_branch", ComplexEval(|a| {
            // The logarithmic form fixes the side of the cuts (num-complex's
            // atanh takes the other one for real |z| > 1).
            let z = *a.first()?;
            let w = (1.0 + z) / (1.0 - z);
            let w = if w.im == 0.0 { Complex64::new(w.re, 0.0) } else { w };
            Some(0.5 * w.ln() + I * PI * *a.get(1)?)
        })),
    ];
    for (name, eval) in evals {
        if let Some(op) = i.graph().ops().lookup(name) {
            i.graph().ops_mut().set_attr(op, eval);
        }
    }
    i.rewrites(
        Tier::Normalize,
        &[
            "branch/log: log_branch(?z, ?k) => ln(abs(?z)) + I * (arg(?z) + 2 * pi * ?k)",
            "branch/sqrt: sqrt_branch(?z, ?k) => abs(?z)^(1/2) * exp(I * (arg(?z) + 2 * pi * ?k) / 2)",
            "branch/root: root_branch(?z, ?n, ?k) => abs(?z)^(1 / ?n) * exp(I * (arg(?z) + 2 * pi * ?k) / ?n)",
            "branch/power: power_branch(?z, ?w, ?k) => exp(?w * log_branch(?z, ?k))",
            "branch/asin: asin_branch(?z, ?k) => ?k * pi + (-1)^?k * asin(?z) if integer(?k)",
            "branch/acos: acos_branch(?z, ?k, ?s) => 2 * ?k * pi + ?s * acos(?z)",
            "branch/atan: atan_branch(?z, ?k) => ?k * pi + atan(?z)",
            "branch/asinh: asinh_branch(?z, ?k) => log_branch(?z + (?z^2 + 1)^(1/2), ?k)",
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
                engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
                let term = Extractor::new(&g, &[root], &SizeCost).build(&mut g, root).unwrap_or(root);
                let got = eval_complex(&g, term, &bindings).unwrap_or_else(|| panic!("cannot evaluate {}", g.display(term)));
                assert!((got - want).norm() < 1e-9, "{src}: {got} vs {want}");
            }
        }
    }
}
