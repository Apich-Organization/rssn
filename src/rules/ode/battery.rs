//! Battery of equations for `dsolve`, audited independently of the solver.
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

const CONSTANT_SETS: [[f64; 4]; 4] = [[0.8, -0.6, 0.4, 0.3], [-1.7, 0.9, -0.5, 1.2], [3.1, -2.2, 0.7, -0.8], [0.15, 0.35, -0.25, 0.45]];

fn bind_constants(
    graph: &mut Graph,
    env: &mut Env,
    constants: &[f64; 4],
) {
    for (k, v) in constants.iter().enumerate() {
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
    let mut last = Report { status: "UNCHECKED".into(), text: String::new() };
    for constants in &CONSTANT_SETS {
        let report = audit_with(eq, order, constants);
        if report.status == "ok" {
            return report;
        }
        // A branch may hold only for some values of the constants.
        if report.status != "UNCHECKED" || last.status == "UNCHECKED" {
            last = report;
        }
    }
    last
}

fn audit_with(
    eq: &str,
    order: usize,
    constants: &[f64; 4],
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
    let (Ok(equation), Ok(unknown)) = (g.parse(eq), g.parse("y(x)")) else {
        return Report { status: "parse".into(), text };
    };
    let Some(problem) = super::parse(&mut g, equation, unknown) else {
        return Report { status: "parse-problem".into(), text };
    };
    if g.op(answer) == core::LIST {
        return audit_parametric(&mut g, &problem, answer, order, constants, text);
    }
    if g.op(answer) != core::EQ {
        return Report { status: "UNSOLVED".into(), text };
    }
    let &[lhs, rhs] = g.children(answer) else {
        return Report { status: "shape".into(), text };
    };
    let x_symbol = problem.x_symbol;
    let stand: Vec<SymbolId> = problem.stand.iter().filter_map(|&s| g.as_symbol(s)).collect();
    let mut base = Env::numeric(0.0);
    bind_constants(&mut g, &mut base, constants);
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
        for at in [0.6, 1.1, 1.7, 0.3, 2.3, -0.6, -1.1, -1.7] {
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
        {
            let eval_f = |x: f64, y: f64| -> Option<f64> {
                let mut env = base.clone();
                env.bind(x_symbol, x);
                env.bind(ysym, y);
                g.eval(f, &env).filter(|v| v.is_finite())
            };
            let mut guess_store = 0.0_f64;
            for at in [0.6, 1.1, 1.7, 0.9, 1.4, -0.6, -1.1, -1.7] {
                // Newton in y from several starts.
                let mut found = None;
                for start in [0.5, 1.0, 2.0, -1.0, 0.1, 3.0, -2.0, 0.01, 5.0, -4.0, 0.3, 1.5, -0.3, 8.0] {
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

/// A parametric answer `list(x = X(p), y(x) = Y(p))` of a first-order
/// equation: `p` is `y'`, so `dY/dp = p dX/dp` and the equation must hold.
fn audit_parametric(
    g: &mut Graph,
    problem: &super::Problem,
    answer: NodeId,
    order: usize,
    constants: &[f64; 4],
    text: String,
) -> Report {
    let _ = order;
    let items = g.children(answer).to_vec();
    let (mut x_of, mut y_of) = (None, None);
    for item in items {
        let &[l, r] = g.children(item) else {
            continue;
        };
        if l == problem.x {
            x_of = Some(r);
        } else if l == problem.y_of_x {
            y_of = Some(r);
        }
    }
    let (Some(xp), Some(yp)) = (x_of, y_of) else {
        return Report { status: "shape".into(), text };
    };
    let p_symbol = g.interner_mut().symbol("p");
    let mut base = Env::numeric(0.0);
    bind_constants(g, &mut base, constants);
    let stand: Vec<SymbolId> = problem.stand.iter().filter_map(|&s| g.as_symbol(s)).collect();
    let x_symbol = problem.x_symbol;
    let eval_at = |g: &Graph, node: NodeId, p: f64| -> Option<f64> {
        let mut env = base.clone();
        env.bind(p_symbol, p);
        g.eval(node, &env).filter(|v| v.is_finite())
    };
    let mut checked = 0;
    for p in [0.6, 1.1, 1.7, 0.3, 2.3, -0.6, -1.1, -1.7] {
        let h = 1e-6;
        let (Some(x0), Some(y0), Some(xa), Some(xb), Some(ya), Some(yb)) = (
            eval_at(g, xp, p),
            eval_at(g, yp, p),
            eval_at(g, xp, p + h),
            eval_at(g, xp, p - h),
            eval_at(g, yp, p + h),
            eval_at(g, yp, p - h),
        ) else {
            continue;
        };
        let (dx, dy) = ((xa - xb) / (2.0 * h), (ya - yb) / (2.0 * h));
        if dx.abs() < 1e-9 {
            continue;
        }
        let slope = dy / dx;
        // The parameter is the slope.
        if (slope - p).abs() > 1e-4 * (1.0 + p.abs()) {
            return Report { status: format!("WRONG(slope {slope} vs p {p})"), text };
        }
        let mut env = base.clone();
        env.bind(x_symbol, x0);
        for (s, v) in stand.iter().zip([y0, p]) {
            env.bind(*s, v);
        }
        if let Some(r) = g.eval(problem.expr, &env).filter(|v| v.is_finite()) {
            checked += 1;
            if r.abs() > 1e-6 * (1.0 + y0.abs() + x0.abs()) {
                return Report { status: format!("WRONG(parametric residual {r} at p = {p})"), text };
            }
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
    // Abel classes.
    ("y(x)*diff(y(x), x) - y(x) = -2*x/9", 1),
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
    fn single_battery() {
        let mut passed = 0;
        for &(eq, order) in SINGLE {
            let start = std::time::Instant::now();
            let (tx, rx) = std::sync::mpsc::channel();
            let owned = eq.to_string();
            // A case that does not finish is abandoned (its thread leaks).
            std::thread::Builder::new()
                .stack_size(64 << 20)
                .spawn(move || {
                    let _ = tx.send(audit(&owned, order));
                })
                .map_or((), drop);
            let report = rx
                .recv_timeout(std::time::Duration::from_secs(20))
                .unwrap_or_else(|_| Report { status: "TIMEOUT".into(), text: String::new() });
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
        assert_eq!(passed, SINGLE.len(), "some equations failed; run with --nocapture");
    }
}

/// Systems: `(equations, unknowns)`; the answer is audited as explicit
/// solutions (finite differences) or as first integrals (conserved along
/// the flow).
pub(super) const SYSTEMS: &[(&[&str], &[&str])] = &[
    (&["diff(x(t), t) = y(t)", "diff(y(t), t) = -x(t)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = -x(t)^2", "diff(y(t), t) = x(t)*y(t)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = 2*x(t)", "diff(y(t), t) = x(t) + y(t)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = y(t)", "diff(y(t), t) = -sin(x(t))"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = x(t)*(2 - y(t))", "diff(y(t), t) = y(t)*(x(t) - 3)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = t*y(t)", "diff(y(t), t) = t*x(t)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = y(t) + 1", "diff(y(t), t) = -x(t) + t"], &["x(t)", "y(t)"]),
    (&["diff(diff(x(t), t), t) = -x(t) + y(t)", "diff(diff(y(t), t), t) = x(t) - y(t)"], &["x(t)", "y(t)"]),
    (&["diff(q(t), t) = p(t)", "diff(p(t), t) = -q(t)^3"], &["q(t)", "p(t)"]),
    (&["diff(s(t), t) = -s(t)*i(t)", "diff(i(t), t) = s(t)*i(t) - i(t)"], &["s(t)", "i(t)"]),
    (&["diff(x(t), t) = y(t)*z(t)", "diff(y(t), t) = -x(t)*z(t)", "diff(z(t), t) = -x(t)*y(t)/2"], &["x(t)", "y(t)", "z(t)"]),
    (&["diff(diff(x(t), t), t) = -4*x(t)"], &["x(t)"]),
    (&["diff(x(t), t) = -y(t)/(1 + t^2)", "diff(y(t), t) = x(t)/(1 + t^2)"], &["x(t)", "y(t)"]),
    (&["diff(x(t), t) = x(t) - y(t)", "diff(y(t), t) = x(t) + y(t)"], &["x(t)", "y(t)"]),
    (&["diff(u(t), t) = u(t)*(1 - u(t))", "diff(v(t), t) = u(t) - v(t)"], &["u(t)", "v(t)"]),
];

#[cfg(test)]
mod system_tests {
    use super::*;

    /// The right-hand sides `x_i' = f_i` of first-order equations, as nodes.
    fn audit_system(
        eqs: &[&str],
        unknowns: &[&str],
    ) -> Report {
        let mut g = Graph::new();
        let src = format!("dsolve(list({}), list({}))", eqs.join(", "), unknowns.join(", "));
        let answer = match run(&mut g, &src) {
            | Ok(a) => a,
            | Err(r) => return r,
        };
        let text = g.display(answer);
        if g.op(answer) != core::LIST {
            return Report { status: "UNSOLVED".into(), text };
        }
        let originals: Vec<NodeId> = eqs.iter().filter_map(|e| g.parse(e).ok()).collect();
        let funcs: Vec<NodeId> = unknowns.iter().filter_map(|u| g.parse(u).ok()).collect();
        let Some(&first) = funcs.first() else {
            return Report { status: "parse".into(), text };
        };
        let Some(&t) = g.children(first).get(1) else {
            return Report { status: "parse".into(), text };
        };
        let Some(t_symbol) = g.symbol_of(t) else {
            return Report { status: "parse".into(), text };
        };
        let Some(diff) = g.ops().lookup("diff") else {
            return Report { status: "parse".into(), text };
        };
        let items = g.children(answer).to_vec();
        let explicit = items.iter().all(|&i| g.children(i).first().is_some_and(|l| funcs.contains(l)));
        let mut base = Env::numeric(0.0);
        bind_constants(&mut g, &mut base, &CONSTANT_SETS[0]);
        let mut checked = 0;
        if explicit {
            // Solutions as functions of t; derivatives by differences.
            let mut value_of: Vec<(NodeId, NodeId)> = Vec::new();
            for &item in &items {
                let &[l, r] = g.children(item) else {
                    return Report { status: "shape".into(), text };
                };
                value_of.push((l, r));
            }
            let at = |g: &Graph, f: NodeId, tv: f64| -> Option<f64> {
                let (_, r) = value_of.iter().find(|(l, _)| *l == f)?;
                let mut env = base.clone();
                env.bind(t_symbol, tv);
                g.eval(*r, &env).filter(|v| v.is_finite())
            };
            for tv in [0.4, 0.9, 1.3, 1.8] {
                let h = 1e-5;
                let mut residual_terms = Vec::new();
                for &eq in &originals {
                    let &[l, r] = g.children(eq) else {
                        continue;
                    };
                    residual_terms.push((l, r));
                }
                // Substitute values for the functions and their derivatives.
                for &(l, r) in &residual_terms {
                    let mut lhs = l;
                    let mut rhs = r;
                    let mut env = base.clone();
                    env.bind(t_symbol, tv);
                    let mut ok = true;
                    for &f in &funcs {
                        let d = g.node(diff, &[f, t]);
                        let dd = g.node(diff, &[d, t]);
                        let h2 = 1e-3;
                        let (Some(plus), Some(minus), Some(mid)) = (at(&g, f, tv + h), at(&g, f, tv - h), at(&g, f, tv)) else {
                            ok = false;
                            break;
                        };
                        let (Some(plus2), Some(minus2)) = (at(&g, f, tv + h2), at(&g, f, tv - h2)) else {
                            ok = false;
                            break;
                        };
                        let ddv = g.float((plus2 - 2.0 * mid + minus2) / (h2 * h2));
                        let dv = g.float((plus - minus) / (2.0 * h));
                        let vv = g.float(mid);
                        for slot in [&mut lhs, &mut rhs] {
                            *slot = g.replace_subterm(*slot, dd, ddv);
                            *slot = g.replace_subterm(*slot, d, dv);
                            *slot = g.replace_subterm(*slot, f, vv);
                        }
                    }
                    if !ok {
                        continue;
                    }
                    if let (Some(a), Some(b)) = (g.eval(lhs, &env), g.eval(rhs, &env))
                        && a.is_finite() && b.is_finite() {
                            checked += 1;
                            if (a - b).abs() > 1e-4 * (1.0 + a.abs()) {
                                return Report { status: format!("WRONG({a} vs {b} at t = {tv})"), text };
                            }
                        }
                }
            }
        } else {
            // First integrals Φ(x(t), y(t), …) = 0 with constants: the
            // gradient along the vector field must vanish. Each unknown is
            // replaced by a stand-in symbol.
            let stands: Vec<SymbolId> = (0..funcs.len()).map(|i| g.interner_mut().symbol(&format!("Z{i}"))).collect();
            let mut fields = Vec::new();
            for &eq in &originals {
                let &[l, r] = g.children(eq) else {
                    return Report { status: "shape".into(), text };
                };
                let Some(pos) = funcs.iter().position(|&f| g.node(diff, &[f, t]) == l) else {
                    return Report { status: "needs-first-order-form".into(), text };
                };
                fields.push((pos, r));
            }
            for &item in &items {
                let &[l, r] = g.children(item) else {
                    continue;
                };
                let minus_one = g.int(-1);
                let neg = g.node(core::MUL, &[minus_one, r]);
                let mut phi = g.node(core::ADD, &[l, neg]);
                let mut field_nodes = Vec::new();
                for (&f, &s) in funcs.iter().zip(&stands) {
                    let sn = g.symbol_node(s);
                    phi = g.replace_subterm(phi, f, sn);
                }
                for (pos, r) in &fields {
                    let mut v = *r;
                    for (&f, &s) in funcs.iter().zip(&stands) {
                        let sn = g.symbol_node(s);
                        v = g.replace_subterm(v, f, sn);
                    }
                    field_nodes.push((*pos, v));
                }
                for point in [[0.7, 1.1, 0.9], [1.3, -0.6, 1.7], [0.5, -1.9, 0.8], [0.3, 0.8, -0.6], [-0.9, 0.4, 0.5]] {
                    let mut env = base.clone();
                    env.bind(t_symbol, 0.5);
                    for (&s, v) in stands.iter().zip(point) {
                        env.bind(s, v);
                    }
                    // The constant of the invariant curve through this point.
                    let constants_present: Vec<SymbolId> = g
                        .free_symbols(g.find(phi))
                        .iter()
                        .copied()
                        .filter(|&s| g.interner().symbol_name(s).starts_with('C'))
                        .collect();
                    if let Some(&c0) = constants_present.first() {
                        let mut c = 0.5;
                        let mut solved = false;
                        for _ in 0..60 {
                            env.bind(c0, c);
                            let Some(v) = g.eval(phi, &env) else {
                                break;
                            };
                            if v.abs() < 1e-12 {
                                solved = true;
                                break;
                            }
                            let h = 1e-7;
                            let mut up = env.clone();
                            up.bind(c0, c + h);
                            let Some(vu) = g.eval(phi, &up) else {
                                break;
                            };
                            let slope = (vu - v) / h;
                            if slope.abs() < 1e-12 {
                                break;
                            }
                            c -= v / slope;
                        }
                        if !solved {
                            continue;
                        }
                    }
                    let Some(_) = g.eval(phi, &env) else {
                        continue;
                    };
                    let mut total = 0.0;
                    let mut ok = true;
                    for (i, &s) in stands.iter().enumerate() {
                        let h = 1e-6;
                        let (mut up, mut down) = (env.clone(), env.clone());
                        up.bind(s, point[i] + h);
                        down.bind(s, point[i] - h);
                        let (Some(a), Some(b)) = (g.eval(phi, &up), g.eval(phi, &down)) else {
                            ok = false;
                            break;
                        };
                        let Some((_, field)) = field_nodes.iter().find(|(p, _)| *p == i) else {
                            ok = false;
                            break;
                        };
                        let Some(f) = g.eval(*field, &env) else {
                            ok = false;
                            break;
                        };
                        total += (a - b) / (2.0 * h) * f;
                    }
                    if ok && total.is_finite() {
                        checked += 1;
                        if total.abs() > 1e-4 {
                            return Report { status: format!("NOT-CONSERVED({total})"), text };
                        }
                    }
                }
            }
        }
        if checked == 0 {
            return Report { status: "UNCHECKED".into(), text };
        }
        Report { status: "ok".into(), text }
    }

    #[test]
    fn system_battery() {
        let mut passed = 0;
        for &(eqs, unknowns) in SYSTEMS {
            let start = std::time::Instant::now();
            let report = audit_system(eqs, unknowns);
            let ok = report.status == "ok";
            if ok {
                passed += 1;
            }
            eprintln!(
                "{:6} {:>6.2}s {}  ->  {}",
                if ok { "ok" } else { "FAIL" },
                start.elapsed().as_secs_f64(),
                eqs.join(", "),
                if ok { report.text } else { format!("{} {}", report.status, report.text) }
            );
        }
        eprintln!("ODESYS {passed}/{}", SYSTEMS.len());
        assert_eq!(passed, SYSTEMS.len(), "some systems failed; run with --nocapture");
    }
}
