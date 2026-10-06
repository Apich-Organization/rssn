//! Symbolic approximants built by the numeric kernels.
//!
//! | operator | value |
//! |---|---|
//! | `chebyshev_approx(f, x, a, b)` | `chebyshev_series(t, c0, c1, ..., cn)` with `t = (2x - a - b)/(b - a)`: the adaptive Chebyshev interpolant `sum c_k T_k(t)` of `f` on `[a, b]`, truncated where the coefficients fall below the tolerance of a numeric run (`1e-14` otherwise) |
//! | `chebyshev_series(t, c0, ..., cn)` | `sum c_k T_k(t)`, evaluated by Clenshaw's recurrence |
//! | `rational_approx(f, x, a, b)` | the AAA barycentric rational approximant `(sum w_j f_j/(x - z_j)) / (sum w_j/(x - z_j))` of `f` on `[a, b]` |
//!
//! `f` is compiled with the [current](crate::backend::current) backend;
//! symbols other than `x` take their values from the bindings of the run.
//! The results are ordinary terms (kept in their numerically stable form):
//! they can be evaluated, compiled, differentiated or printed like any
//! other.

use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::kernels::approx::aaa;
use crate::kernels::approx::clenshaw;
use crate::kernels::approx::ChebSeries;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;

/// Highest Chebyshev degree tried.
const MAX_DEGREE: usize = 1024;
/// Samples for the AAA algorithm.
const SAMPLES: usize = 2000;
/// Most support points of a rational approximant.
const MAX_TERMS: usize = 60;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let heavy = |name: &str| OpDescriptor::new(name, Arity::Fixed(4)).flags(OpFlags::HEAVY).cost(100);
    let series = i.op(OpDescriptor::new("chebyshev_series", Arity::Variadic).cost(4).eval(|a| match a.split_first() {
        | Some((&t, coeffs)) if !coeffs.is_empty() => clenshaw(coeffs, t),
        | _ => f64::NAN,
    }))?;
    let op = i.op(heavy("chebyshev_approx"))?;
    i.kernel("special/chebyshev_approx", Tier::Reduce, Approximate { op, kind: Kind::Chebyshev(series) });
    let op = i.op(heavy("rational_approx"))?;
    i.kernel("special/rational_approx", Tier::Reduce, Approximate { op, kind: Kind::Rational });
    Ok(())
}

#[derive(Clone, Copy)]
enum Kind {
    Chebyshev(OpId),
    Rational,
}

struct Approximate {
    op: OpId,
    kind: Kind,
}

impl Kernel for Approximate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, x, a, b] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        // Pinned: combining the barycentric form into one fraction would
        // destroy its stability.
        self.approximate(cx, f, x, a, b).map_or(Outcome::Pass, Outcome::Pinned)
    }

    fn revisit(&self) -> bool {
        true
    }
}

impl Approximate {
    fn approximate(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        x: NodeId,
        a: NodeId,
        b: NodeId,
    ) -> Option<NodeId> {
        let variable = cx.graph.symbol_of(x)?;
        let limit = |cx: &mut Cx<'_>, n: NodeId| {
            let term = best(cx.graph, n)?;
            cx.graph.eval(term, cx.env).filter(|v| v.is_finite())
        };
        let (lo, hi) = (limit(cx, a)?, limit(cx, b)?);
        if lo >= hi {
            return None;
        }
        let body = best(cx.graph, f)?;
        let bindings: Vec<(SymbolId, f64)> =
            cx.env.bindings().iter().copied().filter(|&(s, _)| s != variable).collect();
        let free = cx.graph.free_symbols(cx.graph.find(body)).to_vec();
        if free.iter().any(|&s| s != variable && !bindings.iter().any(|&(b, _)| b == s)) {
            return None;
        }
        let params: Vec<SymbolId> = bindings.iter().map(|&(s, _)| s).collect();
        let values: Vec<f64> = bindings.iter().map(|&(_, v)| v).collect();
        let compiled = crate::backend::current().compile_multi(cx.graph, &[body], &[variable], &params).ok()?;
        let eval = |t: f64| {
            let mut out = [0.0];
            compiled.eval(&[t], &values, &mut out);
            out[0]
        };
        let tolerance = if cx.env.numeric && cx.env.tolerance > 0.0 { cx.env.tolerance } else { 1e-14 };
        match self.kind {
            | Kind::Chebyshev(chebyshev_series) => {
                let series = ChebSeries::adaptive(eval, lo, hi, tolerance, MAX_DEGREE);
                if series.coeffs.iter().any(|c| !c.is_finite()) {
                    return None;
                }
                // t = (2x - a - b)/(b - a)
                let scale = cx.graph.float(2.0 / (hi - lo));
                let shift = cx.graph.float((lo + hi) / (hi - lo));
                let scaled = mul(cx.graph, &[scale, x]);
                let t = sub(cx.graph, scaled, shift);
                let mut args = Vec::with_capacity(series.coeffs.len() + 1);
                args.push(t);
                args.extend(series.coeffs.iter().map(|&c| cx.graph.float(c)));
                Some(cx.graph.node(chebyshev_series, &args))
            },
            | Kind::Rational => {
                // Chebyshev points of the first kind: all interior, so no
                // support point sits on an end of the interval, where the
                // barycentric form is 0/0.
                let half = 0.5 * (hi - lo);
                let z: Vec<f64> = (0..SAMPLES)
                    .map(|k| {
                        let angle = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * SAMPLES as f64);
                        f64::midpoint(lo, hi) - half * angle.cos()
                    })
                    .collect();
                let samples: Vec<f64> = z.iter().map(|&t| eval(t)).collect();
                if samples.iter().any(|v| !v.is_finite()) {
                    return None;
                }
                let r = aaa(&z, &samples, tolerance, MAX_TERMS).ok()?;
                let (mut numerator, mut denominator) = (Vec::new(), Vec::new());
                for ((&zj, &fj), &wj) in r.z.iter().zip(&r.f).zip(&r.w) {
                    let node = cx.graph.float(zj);
                    let difference = sub(cx.graph, x, node);
                    let pole = powi(cx.graph, difference, -1);
                    let weight = cx.graph.float(wj);
                    let weighted = cx.graph.float(wj * fj);
                    numerator.push(mul(cx.graph, &[weighted, pole]));
                    denominator.push(mul(cx.graph, &[weight, pole]));
                }
                let numerator = add(cx.graph, &numerator);
                let denominator = add(cx.graph, &denominator);
                let reciprocal = powi(cx.graph, denominator, -1);
                Some(mul(cx.graph, &[numerator, reciprocal]))
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Graph;
    use crate::graph::Saturate;
    use crate::rules::standard;

    /// Reduces `src` and evaluates the result at `x = at`.
    fn approximate_at(
        src: &str,
        points: &[f64],
    ) -> Vec<f64> {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &standard()).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        let reduced = Extractor::new(&g, &[root], &ClosedForm).build(&mut g, root).unwrap_or_else(|| panic!("{src} not reduced"));
        let x = g.interner_mut().symbol("x");
        points
            .iter()
            .map(|&p| {
                let mut env = Env::numeric(0.0);
                env.bind(x, p);
                g.eval(reduced, &env).unwrap_or(f64::NAN)
            })
            .collect()
    }

    #[test]
    fn chebyshev_interpolants_are_accurate_terms() {
        let points = [-0.9, -0.3, 0.0, 0.41, 0.97];
        let values = approximate_at("chebyshev_approx(exp(x)*sin(5*x), x, -1, 1)", &points);
        for (&p, v) in points.iter().zip(values) {
            let want = p.exp() * (5.0 * p).sin();
            assert!((v - want).abs() < 1e-12, "at {p}: {v} vs {want}");
        }
        // A shifted interval and a special function.
        let points = [2.0, 2.7, 3.9];
        let values = approximate_at("chebyshev_approx(lgamma(x), x, 2, 4)", &points);
        let want = [0.0, 0.434_820_553_655_104_2, 1.667_580_347_241_739_4];
        for ((&p, v), w) in points.iter().zip(values).zip(want) {
            assert!((v - w).abs() < 1e-10, "at {p}: {v} vs {w}");
        }
    }

    #[test]
    fn rational_approximants_capture_poles_nearby() {
        // 1/(x^2 + 0.01) is hard for polynomials but rational itself.
        let points = [-1.0, -0.05, 0.0, 0.3, 0.8];
        let values = approximate_at("rational_approx(tanh(20*x) + 1/(x^2 + 1/100), x, -1, 1)", &points);
        for (&p, v) in points.iter().zip(values) {
            let want = (20.0 * p).tanh() + 1.0 / (p * p + 0.01);
            assert!((v - want).abs() < 1e-9 * want.abs().max(1.0), "at {p}: {v} vs {want}");
        }
    }
}
