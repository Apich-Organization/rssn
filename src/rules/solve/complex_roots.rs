//! `solve_complex(equation, x)`: all complex roots of a polynomial
//! equation with rational coefficients, exactly.
//!
//! The polynomial is factored over `Q`. A factor of degree one gives a
//! rational number and one of degree two `(-b ± sqrt(D))/(2a)` (with the
//! imaginary unit `I` when `D < 0`); every root of an irreducible factor of
//! degree three or more is the canonical `root_of(j, f0, ..., fm)` object
//! (real roots first, then the others by real and imaginary part), which
//! has complex numeric semantics. The list is in the same canonical order,
//! each root once.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::as_expression;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::op::core;
use crate::rules::poly::rational_function_in;
use crate::rules::poly::rootof::all_roots_as_objects;
use crate::rules::poly::rootof::canonical_order;
use crate::rules::poly::univariate;

/// The kernel of `solve_complex`.
pub(super) struct SolveComplex {
    pub(super) op: OpId,
}

/// `outside^2 * inside` with `inside` free of small square factors.
fn square_part(value: &BigInt) -> (BigInt, BigInt) {
    let mut outside = BigInt::one();
    let mut rest = value.clone();
    let mut p = BigInt::from(2);
    while &p * &p <= rest && p < BigInt::from(100_000) {
        while rest.is_multiple_of(&(&p * &p)) {
            rest /= &p * &p;
            outside *= &p;
        }
        p += 1;
    }
    (outside, rest)
}

fn quadratic(
    graph: &mut Graph,
    c: &BigInt,
    b: &BigInt,
    a: &BigInt,
) -> Option<Vec<((f64, f64), NodeId)>> {
    let discriminant = b * b - BigInt::from(4) * a * c;
    let (outside, inside) = square_part(&discriminant.abs());
    let denominator = BigInt::from(2) * a;
    let centre_value = BigRational::new(-b.clone(), denominator.clone());
    let negative = discriminant.is_negative();
    let half = graph.num(Number::fraction(1, 2)?);
    let radicand = graph.num(Number::Int(inside.clone()));
    let radical = graph.node(core::POW, &[radicand, half]);
    let size = (outside.to_f64()? * inside.to_f64()?.sqrt() / denominator.to_f64()?).abs();
    let centre = graph.num(Number::rat(centre_value.clone()));
    let unit = graph.ops().lookup("I");
    let mut out = Vec::new();
    for sign in [-1_i64, 1] {
        let scale_value = BigRational::new(outside.clone() * BigInt::from(sign), denominator.abs());
        let scale = graph.num(Number::rat(scale_value));
        let mut factors = vec![scale, radical];
        if negative {
            factors.push(graph.node(unit?, &[]));
        }
        let offset = graph.node(core::MUL, &factors);
        let node = if centre_value.is_zero() {
            offset
        } else {
            graph.node(core::ADD, &[centre, offset])
        };
        let re = centre_value.to_f64()?;
        let signed = f64::from(i32::try_from(sign).ok()?) * size;
        out.push((
            if negative {
                (re, signed)
            } else {
                (re + signed, 0.0)
            },
            node,
        ));
    }
    Some(out)
}

impl Kernel for SolveComplex {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[equation, unknown] = graph.children(node) else {
            return Outcome::Pass;
        };
        let expr = as_expression(graph, equation);
        let Some((numerator, _)) = rational_function_in(graph, expr, unknown) else {
            return Outcome::Pass;
        };
        let mut q = numerator;
        while q.last().is_some_and(Zero::is_zero) {
            q.pop();
        }
        if q.len() < 2 {
            return Outcome::Pass;
        }
        let (_, factors) = univariate::factor(&q);
        let mut roots: Vec<((f64, f64), NodeId)> = Vec::new();
        for (factor, _) in factors {
            let factor: Vec<BigInt> = if factor.last().is_some_and(Signed::is_negative) {
                factor.iter().map(|c| -c.clone()).collect()
            } else {
                factor
            };
            match factor.as_slice() {
                | [c0, c1] => {
                    let value = BigRational::new(-c0.clone(), c1.clone());
                    let Some(approx) = value.to_f64() else {
                        return Outcome::Pass;
                    };
                    roots.push(((approx, 0.0), graph.num(Number::rat(value))));
                },
                | [c, b, a] => {
                    let Some(found) = quadratic(graph, c, b, a) else {
                        return Outcome::Pass;
                    };
                    roots.extend(found);
                },
                | _ => {
                    if let Some(found) = super::tower::complex_roots(graph, &factor) {
                        roots.extend(found);
                        continue;
                    }
                    let Some(found) = all_roots_as_objects(graph, &factor) else {
                        return Outcome::Pass;
                    };
                    roots.extend(found);
                },
            }
        }
        roots.sort_by(|a, b| canonical_order(a.0, b.0));
        let nodes: Vec<NodeId> = roots.into_iter().map(|(_, n)| n).collect();
        Outcome::Equal(graph.node(core::LIST, &nodes))
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Graph;
    use crate::graph::Saturate;

    /// The complex values of the items of the list `solve_complex` returns.
    fn complex_roots(src: &str) -> (String, Vec<Complex64>) {
        let rules = crate::rules::standard();
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &rules).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        let node = Extractor::new(&g, &[root], &ClosedForm)
            .build(&mut g, root)
            .unwrap_or(root);
        let text = g.display(node);
        let values = g
            .children(node)
            .to_vec()
            .iter()
            .map(|&item| {
                g.eval_complex(item, &std::collections::HashMap::new())
                    .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
            })
            .collect();
        (text, values)
    }

    #[test]
    fn all_complex_roots_exactly() {
        // x^2 + 1: the imaginary unit.
        let (text, roots) = complex_roots("solve_complex(x^2 + 1, x)");
        assert_eq!(roots.len(), 2, "{text}");
        assert!(
            roots
                .iter()
                .any(|r| (r - Complex64::new(0.0, 1.0)).norm() < 1e-12),
            "{text}"
        );
        assert!(
            roots
                .iter()
                .any(|r| (r - Complex64::new(0.0, -1.0)).norm() < 1e-12),
            "{text}"
        );
        // x^5 - x - 1: one real root and four complex ones as root_of objects.
        let (text, roots) = complex_roots("solve_complex(x^5 - x - 1, x)");
        assert_eq!(roots.len(), 5, "{text}");
        assert!(text.matches("root_of").count() == 5, "{text}");
        for r in &roots {
            let residual = r.powi(5) - r - Complex64::new(1.0, 0.0);
            assert!(residual.norm() < 1e-9, "{r}: {residual}");
        }
        // The real root comes first.
        assert!(roots[0].im == 0.0 && (roots[0].re - 1.167_303_978_261_418_7).abs() < 1e-9);
        // Mixed factors: (x - 2)(x^2 + x + 1)(x^3 - 2).
        let (text, roots) = complex_roots("solve_complex((x - 2)*(x^2 + x + 1)*(x^3 - 2), x)");
        assert_eq!(roots.len(), 6, "{text}");
        for r in &roots {
            let residual = (r - 2.0) * (r * r + r + 1.0) * (r * r * r - 2.0);
            assert!(residual.norm() < 1e-8, "{r}: {residual}");
        }
    }

    #[test]
    fn exact_complex_roots_of_solvable_sextics_and_septics() {
        // Φ7: the imaginary cyclic sextic, exp(2πi j/7) combinations.
        let (text, roots) = complex_roots("solve_complex(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1, x)");
        assert_eq!(roots.len(), 6, "{text}");
        assert!(!text.contains("root_of"), "{text}");
        for r in &roots {
            assert!((r.powi(7) - Complex64::new(1.0, 0.0)).norm() < 1e-9, "{r}");
        }
        // D7(x) - 3 (group F42): one real root and three conjugate pairs, through the
        // Lagrange resolvents.
        let (text, roots) = complex_roots("solve_complex(x^7 - 7*x^5 + 14*x^3 - 7*x - 3, x)");
        assert_eq!(roots.len(), 7, "{text}");
        assert!(!text.contains("root_of"), "{text}");
        for r in &roots {
            let residual = r.powi(7) - 7.0 * r.powi(5) + 14.0 * r.powi(3) - 7.0 * r - 3.0;
            assert!(
                residual.norm() < 1e-7 * (1.0 + r.norm().powi(7)),
                "{r}: {residual}"
            );
        }
    }
}
