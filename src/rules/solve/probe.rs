//! Probe battery for the algebraic solver.
//!
//! Every case is solved with the standard rules and then audited by an
//! oracle that knows nothing about the solver: all real roots of the
//! residual in a window are located by scanning for sign changes, and the
//! reported solutions (with the integer parameter `n` of a general
//! solution run over a range) must contain each of them and must make the
//! residual vanish themselves.

use crate::graph::Budget;
use crate::graph::ClosedForm;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Extractor;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Saturate;

/// `(equation, parameter values, general solution wanted)`.
type Case = (&'static str, &'static [(&'static str, f64)], u8);

const WINDOW: f64 = 12.0;

/// The residual `lhs - rhs` of an equation text.
fn residual_text(eq: &str) -> String {
    match eq.split_once('=') {
        | Some((l, r)) => format!("({l}) - ({r})"),
        | None => eq.to_string(),
    }
}

fn scan_roots(
    graph: &Graph,
    residual: NodeId,
    params: &[(crate::graph::SymbolId, f64)],
    x: crate::graph::SymbolId,
) -> Vec<f64> {
    let f = |t: f64| -> f64 {
        let mut env = Env::numeric(0.0);
        for &(s, v) in params {
            env.bind(s, v);
        }
        env.bind(x, t);
        graph.eval(residual, &env).filter(|v| v.is_finite()).unwrap_or(f64::NAN)
    };
    let steps = 48_000_i32;
    let h = 2.0 * WINDOW / f64::from(steps);
    let mut roots: Vec<f64> = Vec::new();
    let mut previous = f(-WINDOW);
    for i in 1..=steps {
        let a = -WINDOW + h * f64::from(i - 1);
        let b = a + h;
        let fb = f(b);
        let fa = previous;
        previous = fb;
        if fa.is_nan() || fb.is_nan() {
            continue;
        }
        if fa == 0.0 {
            roots.push(a);
            continue;
        }
        if fa * fb < 0.0 {
            let (mut lo, mut hi, mut flo) = (a, b, fa);
            for _ in 0..100 {
                let mid = f64::midpoint(lo, hi);
                let fm = f(mid);
                if fm.is_nan() {
                    break;
                }
                if flo * fm <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            let r = f64::midpoint(lo, hi);
            let scale = 1.0 + fa.abs().max(fb.abs()).min(1e3);
            if f(r).abs() <= 1e-6 * scale && f(r).abs() < 1e-4 {
                roots.push(r);
            }
        }
    }
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
    roots
}

pub(super) struct Report {
    pub status: String,
    pub text: String,
}

pub(super) fn audit(
    case: &Case,
    solver: &str,
) -> Report {
    let (eq, params, mode) = *case;
    let general = mode == 1;
    let rules = crate::rules::standard();
    let mut g = Graph::new();
    let Ok(engine) = Engine::install(&mut g, &rules) else {
        return Report { status: "install".into(), text: String::new() };
    };
    let name = if general { "solve_general" } else { solver };
    let src = format!("{name}({}, x)", eq);
    let Ok(root) = g.parse(&src) else {
        return Report { status: "parse".into(), text: src };
    };
    let n_symbol = g.interner_mut().symbol("n");
    g.assume(n_symbol, crate::graph::Facts::INTEGER);
    engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
    let Some(answer) = Extractor::new(&g, &[root], &ClosedForm).build(&mut g, root) else {
        return Report { status: "UNSOLVED".into(), text: String::new() };
    };
    let text = g.display(answer);
    if g.op(answer) != crate::graph::op::core::LIST {
        return Report { status: "UNSOLVED".into(), text };
    }
    let Ok(residual) = g.parse(&residual_text(eq)) else {
        return Report { status: "parse".into(), text };
    };
    let x = g.interner_mut().symbol("x");
    let params: Vec<(crate::graph::SymbolId, f64)> = params.iter().map(|&(n, v)| (g.interner_mut().symbol(n), v)).collect();
    let items = g.children(answer).to_vec();
    let mut values: Vec<f64> = Vec::new();
    for item in items {
        let uses_n = g.depends_on(g.find(item), n_symbol);
        let ns: Vec<f64> = if uses_n { (-10..=10).map(f64::from).collect() } else { vec![0.0] };
        for k in ns {
            let mut env = Env::numeric(0.0);
            for &(s, v) in &params {
                env.bind(s, v);
            }
            env.bind(n_symbol, k);
            let Some(v) = g.eval(item, &env).filter(|v| v.is_finite()) else {
                // A formula that is not real at these parameter values.
                if !uses_n && params.is_empty() {
                    return Report { status: "NOEVAL".into(), text };
                }
                continue;
            };
            // The residual at the reported value.
            let mut env2 = Env::numeric(0.0);
            for &(s, p) in &params {
                env2.bind(s, p);
            }
            env2.bind(x, v);
            match g.eval(residual, &env2) {
                | Some(r) if r.is_finite() && r.abs() <= 1e-6 * (1.0 + v.abs().powi(3)) => {},
                | other => return Report { status: format!("WRONG(x={v}, r={other:?})"), text },
            }
            values.push(v);
        }
    }
    if mode == 2 {
        return Report { status: if values.is_empty() { "EMPTY".into() } else { "ok".into() }, text };
    }
    let reference = scan_roots(&g, residual, &params, x);
    for r in reference {
        if !values.iter().any(|&v| (v - r).abs() <= 1e-6 * (1.0 + r.abs())) {
            return Report { status: format!("INCOMPLETE(missing {r})"), text };
        }
    }
    Report { status: "ok".into(), text }
}

macro_rules! p {
    ($e:expr) => {
        ($e, &[] as &[(&str, f64)], 0)
    };
    ($e:expr, $params:expr) => {
        ($e, &$params as &[(&str, f64)], 0)
    };
}

macro_rules! verified {
    ($e:expr) => {
        ($e, &[] as &[(&str, f64)], 2)
    };
}

macro_rules! general {
    ($e:expr) => {
        ($e, &[] as &[(&str, f64)], 1)
    };
    ($e:expr, $params:expr) => {
        ($e, &$params as &[(&str, f64)], 1)
    };
}

pub(super) const ALGEBRAIC: &[Case] = &[
    // Polynomials, with parameters.
    p!("x^2 + a*x + b = 0", [("a", -3.0), ("b", 2.0)]),
    p!("a*x^2 = b", [("a", 2.0), ("b", 8.0)]),
    p!("x^3 = a/b", [("a", 2.0), ("b", 5.0)]),
    p!("(x^2 - a)*(x - b) = 0", [("a", 3.0), ("b", 1.0)]),
    p!("x^4 - 5*x^2 + 4 = 0"),
    p!("x^4 + a*x^2 + b = 0", [("a", -5.0), ("b", 4.0)]),
    p!("x^3 - 3*x^2 - x + 3 = 0"),
    p!("x^3 - 2 = 0"),
    p!("x^4 - 4 = 0"),
    p!("x^3 + a*x^2 + a*x + 1 = 0", [("a", 3.0)]),
    p!("x^4 + x^3 - 4*x^2 + x + 1 = 0"),
    p!("x^6 - 7*x^3 - 8 = 0"),
    p!("x^6 - 9*x^3 + 8 = 0"),
    p!("(x^2 + x)^2 - 8*(x^2 + x) + 12 = 0"),
    p!("(x^2 - 3*x)^2 + 2*(x^2 - 3*x) - 8 = 0"),
    p!("x^5 - 5*x^3 + 5*x - 1 = 0"),
    p!("x^5 - 32 = 0"),
    p!("x^5 - x - 1 = 0"),
    p!("x^3 + 3*x + 1 = 0"),
    p!("x^4 - 4*x^3 + x^2 + 6*x + 2 = 0"),
    p!("x^3 - (a + b + c)*x^2 + (a*b + a*c + b*c)*x - a*b*c = 0", [("a", 1.0), ("b", 2.0), ("c", 4.0)]),
    p!("x^2 - (a + b)*x + a*b = 0", [("a", 1.5), ("b", -2.0)]),
    p!("x^3 + p*x + q = 0", [("p", 3.0), ("q", 1.0)]),
    p!("x^3 - 6*x^2 + 12*x - 8 = a", [("a", 1.0)]),
    // Rational equations.
    p!("1/x + 1/(x + 1) = 1"),
    p!("(x + 1)/(x - 1) = (x - 2)/(x + 3)"),
    p!("x/(x - 1) - 1/(x + 1) = 2/(x^2 - 1)"),
    p!("1/(x^2 - 1) + 1/(x - 1) = 1"),
    p!("x + 1/x = 5/2"),
    p!("x^2 + 1/x^2 = 7"),
    p!("1/x^2 - 3/x + 2 = 0"),
    p!("x + 1/x + (x + 1/x)^2 = 6"),
    // Radicals.
    p!("sqrt(x + 1) + sqrt(x - 2) = 3"),
    p!("sqrt(2*x + 3) = x"),
    p!("sqrt(x) + x = 6"),
    p!("x^(1/3) = 2"),
    p!("x^(2/3) = 4"),
    p!("sqrt(x + sqrt(x)) = 2"),
    p!("sqrt(x + 2*sqrt(x - 1)) = 3"),
    p!("sqrt(3*x + 1) - sqrt(x - 1) = 2"),
    p!("(x^2 - 1)^(1/3) = 2"),
    p!("sqrt(x) - x^(1/4) - 2 = 0"),
    p!("sqrt(x^2 - 4*x + 4) = 3"),
    p!("sqrt(x) + 1/sqrt(x) = 5/2"),
    p!("(x + 1)^(1/3) + (x - 1)^(1/3) = 2^(1/3)"),
    p!("sqrt(x + 5) + sqrt(x) = 5/sqrt(x + 5)"),
    // Exponentials and logarithms.
    p!("exp(x) = 5"),
    p!("2^x = 3^(x - 1)"),
    p!("exp(2*x) - 5*exp(x) + 6 = 0"),
    p!("ln(x) + ln(x + 1) = ln(6)"),
    p!("ln(x^2) = 4"),
    p!("ln(x - 1) = 1"),
    p!("exp(x) + exp(-x) = 4"),
    p!("4^x - 3*2^x - 4 = 0"),
    p!("3^(2*x) - 4*3^x + 3 = 0"),
    p!("ln(x)^2 - 3*ln(x) + 2 = 0"),
    p!("2^(2*x) - 2^(x + 1) - 3 = 0"),
    p!("ln(exp(x) + 1) = 2"),
    p!("ln(ln(x)) = 1"),
    p!("2*ln(x) = ln(x + 6)"),
    p!("ln(x + 1) - ln(x - 1) = 1"),
    p!("exp(x)/(exp(x) + 1) = 1/3"),
    p!("3^x + 3^(-x) = 10/3"),
    p!("exp(x)*(x - 1) = 0"),
    p!("exp(x) = a*x + b", [("a", 1.0), ("b", 2.0)]),
    p!("a^x = b*x + c", [("a", 2.0), ("b", 1.0), ("c", 3.0)]),
    p!("exp(x) = x + 2"),
    p!("x*exp(x) = 2"),
    p!("x*ln(x) = 2"),
    p!("x^x = 27"),
    p!("x*2^x = 1"),
    p!("ln(x) = x - 2"),
    p!("exp(x) + x = 0"),
    p!("exp(x) = x^2"),
    p!("x^2*exp(x) = 4"),
    p!("x^(ln(x)) = exp(4)"),
    p!("exp(x) - exp(-x) = 2*x*0 + 3"),
    p!("x^3*ln(x) = 0"),
    p!("exp(3*x) - 2*exp(x) = 0"),
    // Hyperbolic.
    p!("sinh(x) = 1"),
    p!("cosh(x) = 2"),
    p!("tanh(x) = 1/2"),
    p!("sinh(x) + cosh(x) = 2"),
    p!("cosh(2*x) = 3*cosh(x)"),
    p!("sinh(x)^2 - 3 = 0"),
    p!("sinh(2*x) = 3*sinh(x)"),
    // Inverse trigonometric.
    p!("atan(x) = 1"),
    p!("asin(x) = 1/2"),
    p!("acos(2*x) = pi/3"),
    p!("asin(x) = acos(x)"),
    p!("2*atan(x) = pi/2"),
    p!("atan(2*x) + atan(3*x) = pi/4"),
    p!("atan(x)^2 - atan(x) - 2 = 0"),
    p!("sin(asin(x)) = 1/3"),
    // Absolute values.
    p!("abs(x - 1) = 2*x"),
    p!("abs(x) + abs(x - 2) = 4"),
    p!("abs(x^2 - 4) = 3"),
    p!("abs(2*x - 1) = abs(x + 4)"),
    p!("abs(abs(x) - 1) = 1/2"),
    p!("x*abs(x) = 4"),
    p!("abs(x + 1) - abs(x - 1) = 1"),
    p!("abs(x)^2 - 3*abs(x) + 2 = 0"),
    // Substitutions u = f(x).
    p!("exp(2*x) + exp(x) - 6 = 0"),
    p!("(x^2 + 1)^2 - 4*(x^2 + 1) + 3 = 0"),
    p!("(x^3 - x)^2 = 4*(x^3 - x)"),
    p!("x^4 + 1/x^4 = 2"),
    p!("(ln(x))^2 + ln(x) = 2"),
    p!("exp(x)*(exp(x) - 3) = -2"),
    p!("(x + 1)^2 + 3*(x + 1) + 2 = 0"),
    p!("x^2*exp(x^2) = exp(1)"),
    p!("(exp(x) - 1)^2 = 4"),
    p!("(sqrt(x) - 1)^2 = 4"),
    p!("sin(x)^2 + 3*sin(x) - 4 = 0"),
    p!("(x^2 - 2*x)^2 - 2*(x^2 - 2*x) - 3 = 0"),
    p!("x^2 + 4/x^2 - 5 = 0"),
    p!("x^8 - 17*x^4 + 16 = 0"),
    // Principal-branch trigonometric equations.
    verified!("sin(x) = 1/2"),
    verified!("cos(x) = 1/3"),
    verified!("tan(x) = 2"),
    verified!("sin(x) + cos(x) = 1"),
    verified!("sin(2*x) = cos(x)"),
    verified!("2*cos(x)^2 - 1 = 0"),
    // General solutions.
    general!("sin(x) = 1/2"),
    general!("cos(x) = -1/2"),
    general!("tan(x) = 1"),
    general!("2*sin(x)^2 - sin(x) - 1 = 0"),
    general!("sin(x) + cos(x) = 1"),
    general!("sin(2*x) = cos(x)"),
    general!("sin(x) = sin(3*x)"),
    general!("cos(2*x) + 3*sin(x) = 2"),
    general!("sin(x) + sin(3*x) = 0"),
    general!("sqrt(3)*sin(x) + cos(x) = 1"),
    general!("tan(x) = 2*sin(x)"),
    general!("sin(x)^2 = cos(x)^2"),
    general!("cos(x/2) = 1/2"),
    general!("sin(x + pi/3) = 1/2"),
    general!("a*sin(x) + b*cos(x) = c", [("a", 3.0), ("b", 4.0), ("c", 5.0)]),
    general!("a*sin(x) + b*cos(x) = c", [("a", 1.0), ("b", 2.0), ("c", 1.0)]),
    general!("sin(x)*cos(x) = 1/4"),
    general!("cos(3*x) = cos(x)"),
    general!("sin(x)^4 + cos(x)^4 = 1/2"),
    general!("1 + cos(x) = 2*sin(x)^2"),
    general!("tan(x) + 1/tan(x) = 2"),
    general!("cos(x)^2 + 3*cos(x) + 1 = 0"),
    general!("sin(3*x) = 3*sin(x)*0 + 1/2"),
    general!("sin(x) = 2"),
    general!("cos(x) = cos(x/2)"),
    general!("sin(2*x) + sin(x) = 0"),
    general!("cos(x) - sin(x) = 0"),
    general!("tan(2*x) = tan(x)"),
    general!("sin(x)^2 - sin(x) = 0"),
    general!("4*sin(x)*cos(x) = 1"),
    general!("sin(x)^3 = sin(x)"),
    general!("cos(x)*sin(2*x) = 0"),
    general!("2*sin(x) = 1/sin(x) + 1"),
    general!("sin(x) = a", [("a", 0.3)]),
];

/// `(equations, unknowns, parameters, expected number of real tuples)`.
type SystemCase = (&'static [&'static str], &'static [&'static str], &'static [(&'static str, f64)], Option<usize>);

pub(super) const SYSTEMS: &[SystemCase] = &[
    (&["x^2 + y^2 = 25", "x*y = 12"], &["x", "y"], &[], Some(4)),
    (&["x + y = 5", "x*y = 6"], &["x", "y"], &[], Some(2)),
    (&["x^2 - y = 1", "x + y = 3"], &["x", "y"], &[], Some(2)),
    (&["x^2 + y^2 + z^2 = 3", "x = y", "y = z"], &["x", "y", "z"], &[], Some(2)),
    (&["x*y = 1", "x + y = 3"], &["x", "y"], &[], Some(2)),
    (&["x^2 - y^2 = 3", "x*y = 2"], &["x", "y"], &[], Some(2)),
    (&["y = exp(x)", "y = x + 2"], &["x", "y"], &[], Some(2)),
    (&["sin(x) = y", "y = 1/2"], &["x", "y"], &[], None),
    (&["x + y = 1", "x^2 + y^2 = 1"], &["x", "y"], &[], Some(2)),
    (&["exp(x)*exp(y) = 6", "x - y = 0"], &["x", "y"], &[], Some(1)),
    (&["x^2*y = 4", "x*y^2 = 2"], &["x", "y"], &[], Some(1)),
    (&["x + y + z = 6", "x*y*z = 6", "x*y + y*z + z*x = 11"], &["x", "y", "z"], &[], Some(6)),
    (&["x + y = 1"], &["x", "y"], &[], Some(1)),
    (&["x^2 + y^2 = 1", "z = x + y"], &["x", "y", "z"], &[], Some(2)),
    (&["x + y = 3", "x - y = 1", "2*x + y = 5"], &["x", "y"], &[], Some(1)),
    (&["x + y = 1", "x + y = 2"], &["x", "y"], &[], Some(0)),
    (&["x^2 + y^2 = 1", "x + y = a"], &["x", "y"], &[("a", 1.0)], Some(2)),
    (&["y = x^2", "y = 2*x + 3"], &["x", "y"], &[], Some(2)),
    (&["x^2 + y = 7", "x + y^2 = 11"], &["x", "y"], &[], None),
    (&["ln(x) + ln(y) = ln(6)", "x + y = 5"], &["x", "y"], &[], Some(2)),
    (&["x*exp(y) = 1", "y = ln(2)"], &["x", "y"], &[], Some(1)),
    (&["cos(x) = y", "sin(x) = y"], &["x", "y"], &[], None),
    (&["x + y + z = 1", "x - y = 0"], &["x", "y", "z"], &[], Some(1)),
    (&["x*y = 1", "x^2 + y^2 = a"], &["x", "y"], &[("a", 3.0)], Some(4)),
    (&["x^2 = y^2", "x + y = 2"], &["x", "y"], &[], Some(1)),
    (&["x^3 + y^3 = 9", "x + y = 3"], &["x", "y"], &[], Some(2)),
    (&["x = y + 1", "y = z + 1", "z = 1"], &["x", "y", "z"], &[], Some(1)),
    (&["a*x + b*y = e", "c*x + d*y = f"], &["x", "y"], &[("a", 1.0), ("b", 2.0), ("c", 3.0), ("d", 5.0), ("e", 1.0), ("f", 2.0)], Some(1)),
    (&["x^2 + y^2 = 4", "x*y = 1"], &["x", "y"], &[], Some(4)),
    (&["x*y + x + y = 5", "x*y - x - y = -1"], &["x", "y"], &[], Some(2)),
    (&["x^2 - 2*y^2 = 1", "x + y = 2"], &["x", "y"], &[], Some(2)),
    (&["x^2 + y^2 = r^2", "y = m*x"], &["x", "y"], &[("r", 2.0), ("m", 0.5)], Some(2)),
    (&["x*z = 1", "y = 2*z", "x + y = 3"], &["x", "y", "z"], &[], Some(2)),
    (&["x^4 + y^4 = 17", "x^2 + y^2 = 5"], &["x", "y"], &[], Some(8)),
    (&["exp(x) + y = 3", "x - y = 1"], &["x", "y"], &[], Some(1)),
    (&["x^2 + y = 0", "x + y = 0"], &["x", "y"], &[], Some(2)),
    (&["u + v = 3", "u*v = 2", "w = u - v"], &["u", "v", "w"], &[], Some(2)),
];

pub(super) fn audit_system(case: &SystemCase) -> Report {
    let (eqs, unknowns, params, expected) = *case;
    let rules = crate::rules::standard();
    let mut g = Graph::new();
    let Ok(engine) = Engine::install(&mut g, &rules) else {
        return Report { status: "install".into(), text: String::new() };
    };
    let src = format!("solve(list({}), list({}))", eqs.join(", "), unknowns.join(", "));
    let Ok(root) = g.parse(&src) else {
        return Report { status: "parse".into(), text: src };
    };
    engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
    let Some(answer) = Extractor::new(&g, &[root], &ClosedForm).build(&mut g, root) else {
        return Report { status: "UNSOLVED".into(), text: String::new() };
    };
    let text = g.display(answer);
    if g.op(answer) != crate::graph::op::core::LIST {
        return Report { status: "UNSOLVED".into(), text };
    }
    let residuals: Vec<NodeId> = eqs.iter().filter_map(|e| g.parse(&residual_text(e)).ok()).collect();
    let unknown_symbols: Vec<crate::graph::SymbolId> = unknowns.iter().map(|n| g.interner_mut().symbol(n)).collect();
    let param_symbols: Vec<(crate::graph::SymbolId, f64)> = params.iter().map(|&(n, v)| (g.interner_mut().symbol(n), v)).collect();
    let mut real = 0;
    for tuple in g.children(answer).to_vec() {
        let mut env = Env::numeric(0.0);
        for &(s, v) in &param_symbols {
            env.bind(s, v);
        }
        // Free unknowns that occur in the tuple itself.
        for (i, &s) in unknown_symbols.iter().enumerate() {
            env.bind(s, 0.7 + 0.1 * f64::from(i32::try_from(i).unwrap_or(0)));
        }
        let mut values = Vec::new();
        for &item in g.children(tuple) {
            values.push(g.eval(item, &env).filter(|v| v.is_finite()));
        }
        if values.iter().any(Option::is_none) {
            continue;
        }
        let mut env2 = Env::numeric(0.0);
        for &(s, v) in &param_symbols {
            env2.bind(s, v);
        }
        // The free unknowns keep their sample values.
        for (&s, v) in unknown_symbols.iter().zip(&values) {
            env2.bind(s, v.unwrap_or(0.0));
        }
        for &r in &residuals {
            match g.eval(r, &env2) {
                | Some(v) if v.is_finite() && v.abs() < 1e-6 => {},
                | other => return Report { status: format!("WRONG({other:?})"), text },
            }
        }
        real += 1;
    }
    if let Some(n) = expected {
        if real != n {
            return Report { status: format!("COUNT(got {real}, want {n})"), text };
        }
    } else if real == 0 {
        return Report { status: "EMPTY".into(), text };
    }
    Report { status: "ok".into(), text }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "reports the pass rate of the systems probe battery"]
    fn system_battery() {
        let mut passed = 0;
        for case in SYSTEMS {
            let start = std::time::Instant::now();
            let report = audit_system(case);
            let ok = report.status == "ok";
            if ok {
                passed += 1;
            }
            eprintln!(
                "{:6} {:>6.2}s {}  ->  {}",
                if ok { "ok" } else { "FAIL" },
                start.elapsed().as_secs_f64(),
                case.0.join(", "),
                if ok { report.text } else { format!("{} {}", report.status, report.text) }
            );
        }
        eprintln!("SYSTEMS {passed}/{}", SYSTEMS.len());
    }

    #[test]
    #[ignore = "reports the pass rate of the algebraic probe battery"]
    fn algebraic_battery() {
        let mut passed = 0;
        for case in ALGEBRAIC {
            let start = std::time::Instant::now();
            let report = audit(case, "solve");
            let ok = report.status == "ok";
            if ok {
                passed += 1;
            }
            eprintln!(
                "{:6} {:>6.2}s {}{}  ->  {}",
                if ok { "ok" } else { "FAIL" },
                start.elapsed().as_secs_f64(),
                if case.2 == 1 { "[general] " } else { "" },
                case.0,
                if ok { report.text } else { format!("{} {}", report.status, report.text) }
            );
        }
        eprintln!("ALGEBRAIC {passed}/{}", ALGEBRAIC.len());
    }
}

#[cfg(test)]
mod scratch {
    #[test]
    #[ignore = "prints the reduction of the expressions in $SCRATCH, separated by ';;'"]
    fn scratch() {
        let rules = crate::rules::standard();
        for src in std::env::var("SCRATCH").unwrap_or_default().split(";;") {
            if src.trim().is_empty() {
                continue;
            }
            let start = std::time::Instant::now();
            let (text, reduced) = crate::rules::testing::reduce_with(&rules, src.trim(), &[]);
            eprintln!("{} {:5.2}s  {}  =>  {}", if reduced { "R" } else { "-" }, start.elapsed().as_secs_f64(), src.trim(), text);
        }
    }
}
