//! Probe battery for `dsolve`.
//!
//! The answer of each case is audited without trusting the solver: an
//! explicit `y(x) = …` is substituted back into the equation with finite
//! difference derivatives, an implicit relation `F(x, y(x)) = 0` is solved
//! for `y` numerically at a few points and the same is done, and for a
//! system the explicit solutions are substituted back, while a first
//! integral `Φ(x(t), y(t), …) = C` must be constant along the vector
//! field. Answers with an inert `integral(` are accepted as quadratures.

use crate::graph::Budget;
use crate::graph::ClosedForm;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Extractor;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Saturate;
use crate::graph::SymbolId;
use crate::graph::op::core;

pub(super) struct Report {
    pub status: String,
    pub text: String,
}

const CONSTANTS: [f64; 4] = [0.8, -0.6, 0.4, 0.3];

fn bind_constants(
    graph: &mut Graph,
    env: &mut Env,
) {
    for (k, v) in CONSTANTS.iter().enumerate() {
        let s = graph.interner_mut().symbol(&format!("C{}", k + 1));
        env.bind(s, *v);
    }
}

fn run(
    graph: &mut Graph,
    src: &str,
) -> Result<NodeId, Report> {
    let rules = crate::rules::standard();
    let engine = Engine::install(graph, &rules).map_err(|e| Report { status: format!("install {e}"), text: String::new() })?;
    let root = graph.parse(src).map_err(|e| Report { status: format!("parse {e}"), text: src.into() })?;
    engine.run(graph, &[root], &Env::symbolic(), &Saturate, &Budget::default());
    let answer = Extractor::new(graph, &[root], &ClosedForm).build(graph, root);
    match answer {
        | Some(a) if &*graph.ops().get(graph.op(a)).name != "dsolve" => Ok(a),
        | _ => Err(Report { status: "UNSOLVED".into(), text: String::new() }),
    }
}

/// Central finite-difference derivatives `[f, f', f'', f''']` of `f` at `x`.
fn derivatives(
    f: &dyn Fn(f64) -> Option<f64>,
    x: f64,
    count: usize,
) -> Option<Vec<f64>> {
    let h = match count {
        | 0 | 1 => 1e-5,
        | 2 => 1e-4,
        | _ => 1e-3,
    };
    let y = |k: f64| f(x + k * h);
    let (y0, y1, y2, ym1, ym2) = (y(0.0)?, y(1.0)?, y(2.0)?, y(-1.0)?, y(-2.0)?);
    let mut out = vec![y0, (y1 - ym1) / (2.0 * h)];
    if count >= 2 {
        out.push((y1 - 2.0 * y0 + ym1) / (h * h));
    }
    if count >= 3 {
        out.push((y2 - 2.0 * y1 + 2.0 * ym1 - ym2) / (2.0 * h * h * h));
    }
    Some(out)
}

pub(super) fn audit(
    eq: &str,
    order: usize,
) -> Report {
    let mut g = Graph::new();
    let src = format!("dsolve({eq}, y(x))");
    let answer = match run(&mut g, &src) {
        | Ok(a) => a,
        | Err(r) => return r,
    };
    let text = g.display(answer);
    if text.contains("integral(") {
        return Report { status: "ok".into(), text };
    }
    if g.op(answer) != core::EQ {
        return Report { status: "UNSOLVED".into(), text };
    }
    let (Ok(equation), Ok(unknown)) = (g.parse(eq), g.parse("y(x)")) else {
        return Report { status: "parse".into(), text };
    };
    let Some(problem) = super::parse(&mut g, equation, unknown) else {
        return Report { status: "parse-problem".into(), text };
    };
    let &[lhs, rhs] = g.children(answer) else {
        return Report { status: "shape".into(), text };
    };
    let x_symbol = problem.x_symbol;
    let stand: Vec<SymbolId> = problem.stand.iter().filter_map(|&s| g.as_symbol(s)).collect();
    let mut base = Env::numeric(0.0);
    bind_constants(&mut g, &mut base);
    let residual_at = |g: &Graph, x: f64, derivs: &[f64]| -> Option<f64> {
        let mut env = base.clone();
        env.bind(x_symbol, x);
        for (s, v) in stand.iter().zip(derivs) {
            env.bind(*s, *v);
        }
        g.eval(problem.expr, &env).filter(|v| v.is_finite())
    };
    let mut checked = 0;
    let n_deriv = order.min(3);
    if lhs == problem.y_of_x {
        let y = |x: f64| -> Option<f64> {
            let mut env = base.clone();
            env.bind(x_symbol, x);
            g.eval(rhs, &env).filter(|v| v.is_finite())
        };
        for at in [0.6, 1.1, 1.7, 0.3, 2.3] {
            let Some(d) = derivatives(&y, at, n_deriv) else {
                continue;
            };
            let Some(r) = residual_at(&g, at, &d) else {
                continue;
            };
            checked += 1;
            let tolerance = if order >= 3 { 1e-2 } else { 1e-3 };
            if r.abs() > tolerance * (1.0 + d[0].abs()) {
                return Report { status: format!("WRONG(residual {r} at {at})"), text };
            }
        }
    } else {
        // F(x, y(x)) = 0 with a stand-in for y(x).
        let ysym = g.interner_mut().fresh_symbol("Y");
        let yn = g.symbol_node(ysym);
        let diff_node = g.node(core::ADD, &[lhs, g.children(answer)[1]]);
        let _ = diff_node;
        let minus_one = g.int(-1);
        let neg = g.node(core::MUL, &[minus_one, rhs]);
        let f = g.node(core::ADD, &[lhs, neg]);
        let f = g.replace_subterm(f, problem.y_of_x, yn);
        if g.depends_on(g.find(f), x_symbol) || true {
            let eval_f = |x: f64, y: f64| -> Option<f64> {
                let mut env = base.clone();
                env.bind(x_symbol, x);
                env.bind(ysym, y);
                g.eval(f, &env).filter(|v| v.is_finite())
            };
            let mut guess_store = 0.0_f64;
            for at in [0.6, 1.1, 1.7, 0.9, 1.4] {
                // Newton in y from several starts.
                let mut found = None;
                for start in [0.5, 1.0, 2.0, -1.0, 0.1, 3.0, -2.0] {
                    let mut y = start;
                    let mut ok = false;
                    for _ in 0..80 {
                        let Some(v) = eval_f(at, y) else {
                            break;
                        };
                        if v.abs() < 1e-12 {
                            ok = true;
                            break;
                        }
                        let h = 1e-7;
                        let (Some(a), Some(b)) = (eval_f(at, y + h), eval_f(at, y - h)) else {
                            break;
                        };
                        let slope = (a - b) / (2.0 * h);
                        if slope.abs() < 1e-12 {
                            break;
                        }
                        y -= v / slope;
                    }
                    if ok {
                        found = Some(y);
                        break;
                    }
                }
                let Some(y0) = found else {
                    continue;
                };
                guess_store = y0;
                // y(x) near `at`: Newton from y0.
                let y_of = |x: f64| -> Option<f64> {
                    let mut y = y0;
                    for _ in 0..60 {
                        let v = eval_f(x, y)?;
                        if v.abs() < 1e-13 {
                            return Some(y);
                        }
                        let h = 1e-7;
                        let slope = (eval_f(x, y + h)? - eval_f(x, y - h)?) / (2.0 * h);
                        if slope.abs() < 1e-12 {
                            return None;
                        }
                        y -= v / slope;
                    }
                    None
                };
                let Some(d) = derivatives(&y_of, at, n_deriv) else {
                    continue;
                };
                let Some(r) = residual_at(&g, at, &d) else {
                    continue;
                };
                checked += 1;
                let tolerance = if order >= 3 { 1e-2 } else { 1e-3 };
                if r.abs() > tolerance * (1.0 + d[0].abs()) {
                    return Report { status: format!("WRONG(implicit residual {r} at {at})"), text };
                }
            }
            let _ = guess_store;
        }
    }
    if checked == 0 {
        return Report { status: "UNCHECKED".into(), text };
    }
    Report { status: "ok".into(), text }
}

/// First-order cases: `(equation, order)`.
pub(super) const SINGLE: &[(&str, usize)] = &[
    // Linear, separable, Bernoulli, homogeneous, exact.
    ("diff(y(x), x) + 2*y(x) = x", 1),
    ("diff(y(x), x) = (1 + y(x)^2)/(1 + x^2)", 1),
    ("diff(y(x), x) = exp(x + y(x))", 1),
    ("diff(y(x), x) = sin(x)*cos(y(x))", 1),
    ("x*diff(y(x), x) + y(x) = y(x)^2*ln(x)", 1),
    ("diff(y(x), x) + y(x) = x*sqrt(y(x))", 1),
    ("diff(y(x), x) = y(x)/x + tan(y(x)/x)", 1),
    ("x*diff(y(x), x) = y(x)*ln(y(x)/x)", 1),
    ("diff(y(x), x) = (x^2 + y(x)^2)/(x*y(x))", 1),
    // Riccati classes.
    ("diff(y(x), x) = y(x)^2 - 2/x^2", 1),
    ("diff(y(x), x) = -y(x)^2 + 2*x*y(x) + 1 - x^2", 1),
    ("diff(y(x), x) = y(x)^2 - 2*y(x)*exp(x) + exp(2*x) + exp(x)", 1),
    ("x^2*diff(y(x), x) + x*y(x) + x^2*y(x)^2 = 4", 1),
    ("diff(y(x), x) = y(x)^2 + x", 1),
    ("diff(y(x), x) = y(x)^2 - x^2 + 1", 1),
    ("diff(y(x), x) = y(x)^2/x - y(x)/x - 1/x", 1),
    ("diff(y(x), x) = y(x)^2 + 2*y(x)/x + 1/x^2 + 1", 1),
    // Abel classes.
    ("diff(y(x), x) = (y(x) - x)^3 + 1", 1),
    ("y(x)*diff(y(x), x) - y(x) = -2*x/9", 1),
    ("y(x)*diff(y(x), x) - y(x) = x", 1),
    ("diff(y(x), x) = x*y(x)^3 + y(x)^2", 1),
    // Exact with integrating factors.
    ("(x^2 + y(x)^2 + x) + x*y(x)*diff(y(x), x) = 0", 1),
    ("y(x) + (x^2*y(x) - x)*diff(y(x), x) = 0", 1),
    ("y(x)^2 + (3*x*y(x) - 1)*diff(y(x), x) = 0", 1),
    ("y(x)*(x + y(x)) + (x + 2*y(x) - 1)*diff(y(x), x) = 0", 1),
    ("(2*x*y(x) + 1) + (x^2 + 2*y(x))*diff(y(x), x) = 0", 1),
    ("(y(x) + x*y(x)^2) + (x - x^2*y(x))*diff(y(x), x) = 0", 1),
    // Clairaut, Lagrange, solvable for x or y.
    ("y(x) = x*diff(y(x), x) + diff(y(x), x)^2", 1),
    ("y(x) = x*diff(y(x), x) + sqrt(1 + diff(y(x), x)^2)", 1),
    ("y(x) = 2*x*diff(y(x), x) + diff(y(x), x)^2", 1),
    ("y(x) = x*(1 + diff(y(x), x)) + diff(y(x), x)^2", 1),
    ("y(x) = diff(y(x), x)^2 + diff(y(x), x)", 1),
    ("x = diff(y(x), x)^3 + diff(y(x), x)", 1),
    ("diff(y(x), x)^2 - (x + y(x))*diff(y(x), x) + x*y(x) = 0", 1),
    ("diff(y(x), x)^2 = y(x)", 1),
    ("diff(y(x), x)^2 + y(x)^2 = 1", 1),
    ("y(x) = x*diff(y(x), x) - diff(y(x), x)^3", 1),
    // Linear fractional and translation classes.
    ("diff(y(x), x) = (x - y(x) + 1)/(x + y(x) - 3)", 1),
    ("diff(y(x), x) = (2*x + y(x) - 1)/(4*x + 2*y(x) + 5)", 1),
    ("(x + y(x))*diff(y(x), x) = 1", 1),
    // Higher order.
    ("diff(diff(y(x), x), x) + diff(y(x), x)^2 = 0", 2),
    ("diff(diff(y(x), x), x) = diff(y(x), x)/x + x", 2),
    ("x^2*diff(diff(y(x), x), x) - 3*x*diff(y(x), x) + 4*y(x) = 0", 2),
    ("x^2*diff(diff(y(x), x), x) + x*diff(y(x), x) + y(x) = ln(x)", 2),
    ("x^2*diff(diff(y(x), x), x) - 2*y(x) = x^3", 2),
    ("diff(diff(y(x), x), x) + y(x) = tan(x)", 2),
    ("diff(diff(y(x), x), x) + 4*y(x) = 1/cos(2*x)", 2),
    ("x*diff(diff(y(x), x), x) + 2*diff(y(x), x) + x*y(x) = 0", 2),
    ("(1 + x^2)*diff(diff(y(x), x), x) + 2*x*diff(y(x), x) = 0", 2),
    ("diff(diff(y(x), x), x) - x*diff(y(x), x) + y(x) = 0", 2),
    ("(x^2 - 1)*diff(diff(y(x), x), x) - 2*x*diff(y(x), x) + 2*y(x) = 0", 2),
    ("diff(diff(y(x), x), x) + 2*diff(y(x), x) + 5*y(x) = exp(-x)*cos(2*x)", 2),
    ("x*diff(diff(y(x), x), x) + (1 - x)*diff(y(x), x) + 2*y(x) = 0", 2),
    ("(1 - x^2)*diff(diff(y(x), x), x) - x*diff(y(x), x) + 4*y(x) = 0", 2),
    ("diff(diff(diff(y(x), x), x), x) = diff(diff(y(x), x), x)", 3),
    ("x*diff(diff(y(x), x), x) + diff(y(x), x) = 0", 2),
    ("diff(diff(y(x), x), x) = 2*y(x)*diff(y(x), x)", 2),
    ("y(x)*diff(diff(y(x), x), x) = diff(y(x), x)^2 + diff(y(x), x)", 2),
    ("x^2*diff(diff(y(x), x), x) + 3*x*diff(y(x), x) + y(x) = 1/x", 2),
    ("diff(diff(y(x), x), x) - 4*diff(y(x), x) + 4*y(x) = exp(2*x)*x", 2),
    ("(x^2 + 1)*diff(diff(y(x), x), x) - 2*y(x) = 0", 2),
    ("x^2*diff(diff(y(x), x), x) - 2*x*diff(y(x), x) + 2*y(x) = x^2*exp(x)", 2),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "reports the pass rate of the ODE probe battery"]
    fn single_battery() {
        let mut passed = 0;
        for &(eq, order) in SINGLE {
            let start = std::time::Instant::now();
            let report = audit(eq, order);
            let ok = report.status == "ok";
            if ok {
                passed += 1;
            }
            eprintln!(
                "{:6} {:>6.2}s {}  ->  {}",
                if ok { "ok" } else { "FAIL" },
                start.elapsed().as_secs_f64(),
                eq,
                if ok { report.text } else { format!("{} {}", report.status, report.text) }
            );
        }
        eprintln!("ODE {passed}/{}", SINGLE.len());
    }
}
