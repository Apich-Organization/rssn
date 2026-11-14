//! A systematic battery of symbolic requests, each verified numerically.
//!
//! Every case is run through the public `Session::compute` with the standard
//! rules and a symbolic target. The answer is then checked against an
//! independent numeric oracle:
//!
//! * antiderivatives are differentiated back (central differences) and
//!   compared with the integrand;
//! * derivatives are compared with central differences of the input;
//! * simplifications must agree with the input at sample points;
//! * limits, sums and definite integrals are compared with a known value;
//! * solutions are substituted back into the equation.
//!
//! `AUDIT_TABLE=<path> cargo test --test symbolic_battery` writes the status
//! table used by `docs/SYMBOLIC_AUDIT.md`.
mod common;
use common::*;
use std::fmt::Write as _;

#[derive(Clone, Copy, Debug)]
enum Kind {
    /// `integral(f, x)`; the answer must differentiate back to `f`.
    Integral,
    /// `diff(f, x)`; compared with finite differences.
    Diff,
    /// Any expression; the answer must equal the input numerically.
    Simplify,
    /// An expression whose value (no free symbols) is given.
    Value(f64),
    /// `defint(...)` or `sum(...)` etc. with the value given.
    Closed(f64),
    /// `solve(eq, x)`; each root must satisfy `eq`.
    Solve,
    /// `taylor(f, x, 0, n)`; compared with `f` near 0 to the truncation order.
    Series(u32),
    /// `solve` with no real root: an empty list is the right answer.
    SolveNone,
    /// A quantity that does not exist (two-sided limit that disagrees, a
    /// divergent sum): the request must stay unevaluated or answer `oo`.
    Undefined,
}

struct Case {
    kind: Kind,
    src: &'static str,
    /// Exact text the answer must have, when the case pins a normal form.
    text: Option<&'static str>,
}

const fn c(
    kind: Kind,
    src: &'static str,
) -> Case {
    Case {
        kind,
        src,
        text: None,
    }
}

const fn t(
    kind: Kind,
    src: &'static str,
    text: &'static str,
) -> Case {
    Case {
        kind,
        src,
        text: Some(text),
    }
}

const POINTS: [f64; 5] = [0.31, 0.57, 0.83, 1.37, 2.11];

fn near(
    a: f64,
    b: f64,
    tol: f64,
) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

/// The status of a case: `ok`, or a one-line complaint.
fn judge(case: &Case) -> (String, String) {
    let o = sym(case.src);
    let mut status = String::from("ok");
    let mut say = |s: &str| {
        if status == "ok" {
            status = s.to_owned();
        }
    };
    if has_float(&o.text) {
        say("FLOAT");
    }
    match case.kind {
        | Kind::Integral => {
            // The request is `integral(f, x)`; recover f.
            let f = integrand(case.src);
            if !o.reduced || o.text.contains("integral(") {
                say("UNEVALUATED");
            } else {
                let h = 1e-5;
                let mut checked = 0;
                for p in POINTS {
                    let (Some(up), Some(dn), Some(want)) = (
                        eval(&o.text, &[("x", p + h)]),
                        eval(&o.text, &[("x", p - h)]),
                        eval(&f, &[("x", p)]),
                    ) else {
                        continue;
                    };
                    let d = (up - dn) / (2.0 * h);
                    if !d.is_finite() || !want.is_finite() {
                        continue;
                    }
                    checked += 1;
                    if !near(d, want, 1e-5) {
                        say("WRONG");
                    }
                }
                if checked == 0 {
                    say("UNCHECKED");
                }
            }
        },
        | Kind::Diff => {
            let f = integrand(case.src);
            if !o.reduced || o.text.contains("diff(") {
                say("UNEVALUATED");
            } else {
                let h = 1e-6;
                for p in POINTS {
                    let (Some(up), Some(dn), Some(got)) = (
                        eval(&f, &[("x", p + h)]),
                        eval(&f, &[("x", p - h)]),
                        eval(&o.text, &[("x", p)]),
                    ) else {
                        continue;
                    };
                    let d = (up - dn) / (2.0 * h);
                    if d.is_finite() && !near(got, d, 1e-4) {
                        say("WRONG");
                    }
                }
            }
        },
        | Kind::Simplify => {
            if !o.reduced {
                say("UNEVALUATED");
            }
            for p in POINTS {
                let (Some(a), Some(b)) = (
                    eval(case.src, &[("x", p), ("y", 0.7)]),
                    eval(&o.text, &[("x", p), ("y", 0.7)]),
                ) else {
                    continue;
                };
                if a.is_finite() && b.is_finite() && !near(b, a, 1e-9) {
                    say("WRONG");
                }
            }
        },
        | Kind::Value(want) | Kind::Closed(want) => {
            if !o.reduced {
                say("UNEVALUATED");
            } else {
                match eval(&o.text, &[]) {
                    | Some(v) if near(v, want, 1e-8) => {},
                    | Some(v) if want.is_infinite() && v == want => {},
                    | _ => say("WRONG"),
                }
            }
        },
        | Kind::Solve => {
            if !o.reduced {
                say("UNEVALUATED");
            } else {
                let eq = integrand(case.src);
                let roots = list_items(&o.text);
                if roots.is_empty() {
                    say("UNEVALUATED");
                }
                for r in roots {
                    let Some(rv) = eval(&r, &[]) else {
                        say("UNCHECKED");
                        continue;
                    };
                    match eval(&eq, &[("x", rv)]) {
                        | Some(v) if v.abs() < 1e-8 => {},
                        | Some(_) => say("WRONG"),
                        | None => say("UNCHECKED"),
                    }
                }
            }
        },
        | Kind::SolveNone => {
            if !o.reduced {
                say("UNEVALUATED");
            } else if !list_items(&o.text).is_empty() {
                say("WRONG");
            }
        },
        | Kind::Undefined => {
            if o.reduced && !o.text.contains("oo") {
                say("WRONG");
            }
        },
        | Kind::Series(n) => {
            if !o.reduced {
                say("UNEVALUATED");
            } else {
                let f = integrand(case.src);
                for p in [0.05, 0.1, -0.08] {
                    let (Some(a), Some(b)) = (eval(&f, &[("x", p)]), eval(&o.text, &[("x", p)]))
                    else {
                        continue;
                    };
                    let tol = 50.0 * f64::powi(p.abs(), i32::try_from(n).unwrap_or(1));
                    if (a - b).abs() > tol.max(1e-9) {
                        say("WRONG");
                    }
                }
            }
        },
    }
    if let Some(want) = case.text
        && o.text != want
    {
        say("FORMAT");
    }
    (status, o.text)
}

/// The first argument of a one-variable request `op(f, x...)` (text).
fn integrand(src: &str) -> String {
    let open = src.find('(').unwrap_or(0);
    let inner = &src[open + 1..];
    let mut depth = 0_i32;
    for (i, ch) in inner.char_indices() {
        match ch {
            | '(' => depth += 1,
            | ')' => depth -= 1,
            | ',' if depth == 0 => return inner[..i].to_owned(),
            | _ => {},
        }
    }
    inner.trim_end_matches(')').to_owned()
}

fn list_items(text: &str) -> Vec<String> {
    let inner = text
        .trim()
        .strip_prefix("list(")
        .and_then(|s| s.strip_suffix(')'))
        .unwrap_or(text.trim());
    let mut out = Vec::new();
    let mut depth = 0_i32;
    let mut cur = String::new();
    for ch in inner.chars() {
        match ch {
            | '(' => {
                depth += 1;
                cur.push(ch);
            },
            | ')' => {
                depth -= 1;
                cur.push(ch);
            },
            | ',' if depth == 0 => {
                out.push(std::mem::take(&mut cur).trim().to_owned());
            },
            | _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_owned());
    }
    out
}

fn cases() -> Vec<Case> {
    use Kind::*;
    let pi = std::f64::consts::PI;
    let e = std::f64::consts::E;
    vec![
        // ---- integration: polynomials, rational functions
        c(Integral, "integral(x^3 - 2*x + 5, x)"),
        c(Integral, "integral(1/x, x)"),
        c(Integral, "integral(1/(x^2 + 1), x)"),
        c(Integral, "integral(1/(x^2 + 4), x)"),
        c(Integral, "integral(1/(x^2 - 1), x)"),
        c(Integral, "integral(1/(x*(x + 1)), x)"),
        c(Integral, "integral(1/(x^2 + 3*x + 2), x)"),
        c(Integral, "integral(1/((x + 1)^2*(x + 2)), x)"),
        c(Integral, "integral((x + 1)/(x^2 + 2*x + 5), x)"),
        c(Integral, "integral(x/(x^2 + 1), x)"),
        c(Integral, "integral(x^2/(x^2 + 1), x)"),
        c(Integral, "integral(1/(x^3 + 1), x)"),
        c(Integral, "integral(1/(x^4 - 1), x)"),
        c(Integral, "integral(1/(x^4 + 1), x)"),
        c(Integral, "integral((3*x + 2)/((x - 1)*(x + 3)), x)"),
        c(Integral, "integral(1/(x*(x^2 + 1)), x)"),
        c(Integral, "integral((x^3 + x)/(x^2 + 1)^2, x)"),
        c(Integral, "integral(1/(x^2 + 2*x + 2)^2, x)"),
        // ---- integration: elementary, trig, by parts, substitution
        c(Integral, "integral(sin(x), x)"),
        c(Integral, "integral(sin(x)^2, x)"),
        c(Integral, "integral(cos(x)^2, x)"),
        c(Integral, "integral(sin(x)*cos(x), x)"),
        c(Integral, "integral(sin(x)^3, x)"),
        c(Integral, "integral(cos(x)^3, x)"),
        c(Integral, "integral(sin(x)^4, x)"),
        c(Integral, "integral(sin(2*x)*cos(3*x), x)"),
        c(Integral, "integral(sin(x)*sin(2*x), x)"),
        c(Integral, "integral(cos(x)*cos(3*x), x)"),
        c(Integral, "integral(tan(x), x)"),
        c(Integral, "integral(tan(x)^2, x)"),
        c(Integral, "integral(1/cos(x)^2, x)"),
        c(Integral, "integral(1/sin(x), x)"),
        c(Integral, "integral(1/cos(x), x)"),
        c(Integral, "integral(1/(1 + cos(x)), x)"),
        c(Integral, "integral(1/(2 + sin(x)), x)"),
        c(Integral, "integral(exp(x), x)"),
        c(Integral, "integral(exp(-x^2), x)"),
        c(Integral, "integral(exp(2*x)*sin(x), x)"),
        c(Integral, "integral(x*exp(x), x)"),
        c(Integral, "integral(x^2*exp(x), x)"),
        c(Integral, "integral(x*exp(-x), x)"),
        c(Integral, "integral(x*sin(x), x)"),
        c(Integral, "integral(x*cos(x), x)"),
        c(Integral, "integral(x^2*sin(x), x)"),
        c(Integral, "integral(x*ln(x), x)"),
        c(Integral, "integral(ln(x), x)"),
        c(Integral, "integral(ln(x)^2, x)"),
        c(Integral, "integral(x^2*ln(x), x)"),
        c(Integral, "integral(atan(x), x)"),
        c(Integral, "integral(asin(x), x)"),
        c(Integral, "integral(x*atan(x), x)"),
        c(Integral, "integral(exp(x)*cos(x), x)"),
        c(Integral, "integral(x*exp(x^2), x)"),
        c(Integral, "integral(2*x*cos(x^2), x)"),
        c(Integral, "integral(cos(x)/sin(x), x)"),
        c(Integral, "integral(ln(x)/x, x)"),
        c(Integral, "integral(1/(x*ln(x)), x)"),
        c(Integral, "integral(exp(x)/(1 + exp(x)), x)"),
        c(Integral, "integral(1/(1 + exp(x)), x)"),
        c(Integral, "integral(sinh(x), x)"),
        c(Integral, "integral(cosh(x)^2, x)"),
        c(Integral, "integral(sqrt(x), x)"),
        c(Integral, "integral(1/sqrt(x), x)"),
        c(Integral, "integral(sqrt(1 - x^2), x)"),
        c(Integral, "integral(1/sqrt(1 - x^2), x)"),
        c(Integral, "integral(1/sqrt(x^2 + 1), x)"),
        c(Integral, "integral(sqrt(x^2 + 1), x)"),
        c(Integral, "integral(x/sqrt(x^2 + 1), x)"),
        c(Integral, "integral(1/sqrt(x^2 - 1), x)"),
        c(Integral, "integral(x*sqrt(x + 1), x)"),
        c(Integral, "integral(sin(x)/x, x)"),
        c(Integral, "integral(exp(x)/x, x)"),
        c(Integral, "integral(1/ln(x), x)"),
        c(Integral, "integral(sin(x^2), x)"),
        c(Integral, "integral(x^x, x)"),
        // ---- differentiation
        c(Diff, "diff(x^5, x)"),
        c(Diff, "diff(sin(x)*cos(x), x)"),
        c(Diff, "diff(exp(x^2), x)"),
        c(Diff, "diff(ln(x^2 + 1), x)"),
        c(Diff, "diff(atan(x), x)"),
        c(Diff, "diff(asin(x), x)"),
        c(Diff, "diff(x^x, x)"),
        c(Diff, "diff(sqrt(1 + x^2), x)"),
        c(Diff, "diff(tan(x), x)"),
        c(Diff, "diff(sinh(x)/x, x)"),
        c(Diff, "diff(x/(x^2 + 1), x)"),
        c(Diff, "diff(exp(sin(x)), x)"),
        c(Diff, "diff(ln(sin(x)), x)"),
        c(Diff, "diff(abs(x)^3, x)"),
        c(Diff, "diff(atan(x)/x, x)"),
        c(Diff, "diff(acosh(x + 1), x)"),
        c(Diff, "diff(gamma(x), x)"),
        c(Diff, "diff(erf(x), x)"),
        // ---- simplification
        t(Simplify, "sin(x)^2 + cos(x)^2", "1"),
        t(Simplify, "exp(ln(x))", "x"),
        t(Simplify, "ln(exp(x))", "x"),
        t(Simplify, "(x^2 - 1)/(x - 1)", "x + 1"),
        t(Simplify, "sin(-x)", "-sin(x)"),
        t(Simplify, "cos(-x)", "cos(x)"),
        t(Simplify, "tan(-x)", "-tan(x)"),
        t(Simplify, "asin(-x)", "-asin(x)"),
        t(Simplify, "atan(-x)", "-atan(x)"),
        t(Simplify, "sinh(-x)", "-sinh(x)"),
        t(Simplify, "cosh(-x)", "cosh(x)"),
        t(Simplify, "tanh(-x)", "-tanh(x)"),
        t(Simplify, "asinh(-x)", "-asinh(x)"),
        t(Simplify, "erf(-x)", "-erf(x)"),
        c(Simplify, "sqrt(x^2)"),
        c(Simplify, "(x + 1)^2 - (x^2 + 2*x + 1)"),
        c(Simplify, "(x^3 - 1)/(x^2 + x + 1)"),
        c(Simplify, "cos(x)^2 - sin(x)^2"),
        c(Simplify, "2*sin(x)*cos(x)"),
        c(Simplify, "exp(x)*exp(2*x)"),
        c(Simplify, "ln(x^2)"),
        c(Simplify, "ln(x*y)"),
        c(Simplify, "x/x"),
        c(Simplify, "(x^2 - 4)/(x + 2)"),
        c(Simplify, "1/(1/x)"),
        c(Simplify, "sqrt(x)*sqrt(x)"),
        c(Simplify, "exp(x)^2"),
        c(Simplify, "sin(x)/cos(x)"),
        c(Simplify, "cosh(x)^2 - sinh(x)^2"),
        c(Simplify, "(x + y)^2 - x^2 - 2*x*y - y^2"),
        // ---- special values
        t(Value(0.0), "si(0)", "0"),
        t(Value(0.0), "erf(0)", "0"),
        t(Value(1.0), "gamma(1)", "1"),
        t(Value(2.0), "gamma(3)", "2"),
        t(Value(24.0), "gamma(5)", "24"),
        t(Value(pi.sqrt()), "gamma(1/2)", "pi^(1/2)"),
        t(Value(pi * pi / 6.0), "zeta(2)", "1/6*pi^2"),
        t(Value(pi.powi(4) / 90.0), "zeta(4)", "1/90*pi^4"),
        t(Value(0.0), "ci(1) - ci(1)", "0"),
        t(Value(0.0), "li(0)", "0"),
        t(Value(0.0), "ei(-oo)", "0"),
        t(Value(1.0), "erf(oo)", "1"),
        t(Value(0.0), "erfc(oo)", "0"),
        t(Value(1.0), "erfc(0)", "1"),
        t(Value(0.0), "ln(1)", "0"),
        t(Value(0.0), "asin(0)", "0"),
        t(Value(pi / 2.0), "si(oo)", "1/2*pi"),
        t(Value(2.0), "gamma(3)", "2"),
        // ---- limits
        t(Value(0.5), "limit((1 - cos(x))/x^2, x, 0)", "1/2"),
        t(Value(0.5), "limit((x^2 + 1)/(2*x^2 - x), x, oo)", "1/2"),
        t(Value(0.0), "limit(x*ln(x), x, 0, plus)", "0"),
        t(Value(0.0), "limit(x*ln(x), x, 0)", "0"),
        t(Value(1.0), "limit(x^x, x, 0, plus)", "1"),
        t(Value(e), "limit((1 + 1/x)^x, x, oo)", "E"),
        t(Value(0.0), "limit(x*exp(-x), x, oo)", "0"),
        t(Value(0.0), "limit(ln(x)/x, x, oo)", "0"),
        t(Value(0.0), "limit(sin(x)/x, x, oo)", "0"),
        t(Value(1.0), "limit((exp(x) - 1)/x, x, 0)", "1"),
        t(Value(0.0), "limit(x^2*ln(x), x, 0, plus)", "0"),
        t(Value(1.0), "limit(sin(x)/x, x, 0)", "1"),
        t(Value(0.5), "limit((exp(x) - 1 - x)/x^2, x, 0)", "1/2"),
        t(Value(1.0 / 3.0), "limit((tan(x) - x)/x^3, x, 0)", "1/3"),
        t(Value(0.0), "limit(x^3*exp(-x), x, oo)", "0"),
        t(Value(1.0), "limit(x^(1/x), x, oo)", "1"),
        t(Value(1.0), "limit(x^(1/x), x, oo)", "1"),
        t(Value(4.0), "limit((x^2 - 4)/(x - 2), x, 2)", "4"),
        t(Value(f64::INFINITY), "limit(exp(x)/x^10, x, oo)", "oo"),
        t(Value(0.0), "limit(1/x, x, oo)", "0"),
        t(Value(1.0), "limit(x*sin(1/x), x, oo)", "1"),
        t(Value(e), "limit((1 + x)^(1/x), x, 0)", "E"),
        t(Value(0.5), "limit(sqrt(x^2 + x) - x, x, oo)", "1/2"),
        t(Value(1.0), "limit(atan(x)/x, x, 0)", "1"),
        t(Value(pi / 2.0), "limit(atan(x), x, oo)", "1/2*pi"),
        t(Value(0.0), "limit(ln(1 + x)/exp(x), x, oo)", "0"),
        t(Value(2.0), "limit(sinh(2*x)/x, x, 0)", "2"),
        t(Value(0.5), "limit((1 - cos(x))/(x*sin(x)), x, 0)", "1/2"),
        t(Value(-1.0), "limit(ln(x)/(1 - x), x, 1)", "-1"),
        // ---- definite integrals
        t(Closed(1.0), "defint(2*x, x, 0, 1)", "1"),
        t(Closed(pi / 4.0), "defint(1/(1 + x^2), x, 0, 1)", "1/4*pi"),
        t(Closed(1.0), "defint(sin(x), x, 0, pi/2)", "1"),
        t(
            Closed(0.946_083_070_367_183),
            "defint(sin(x)/x, x, 0, 1)",
            "si(1)",
        ),
        t(Closed(1.0), "defint(exp(-x), x, 0, oo)", "1"),
        t(
            Closed(pi.sqrt() / 2.0),
            "defint(exp(-x^2), x, 0, oo)",
            "1/2*pi^(1/2)",
        ),
        t(Closed(pi / 2.0), "defint(1/(1 + x^2), x, 0, oo)", "1/2*pi"),
        t(Closed(-1.0), "defint(ln(x), x, 0, 1)", "-1"),
        t(Closed(-0.25), "defint(x*ln(x), x, 0, 1)", "-1/4"),
        t(Closed(pi / 4.0), "defint(sqrt(1 - x^2), x, 0, 1)", "1/4*pi"),
        t(Closed(0.5), "defint(x*exp(-x^2), x, 0, oo)", "1/2"),
        t(Closed(1.0), "defint(x*exp(-x), x, 0, oo)", "1"),
        t(Closed(2.0), "defint(sin(x), x, 0, pi)", "2"),
        t(Closed(pi / 2.0), "defint(sin(x)^2, x, 0, pi)", "1/2*pi"),
        t(Closed(1.0 / 3.0), "defint(x^2, x, 0, 1)", "1/3"),
        t(Closed(pi), "defint(1/(1 + x^2), x, -oo, oo)", "pi"),
        t(
            Closed(1.0 - 2.0 / e),
            "defint(x*exp(-x), x, 0, 1)",
            "1 - 2*exp(-1)",
        ),
        // ---- sums
        t(Closed(pi * pi / 6.0), "sum(1/k^2, k, 1, oo)", "1/6*pi^2"),
        t(Closed(2.0), "sum(1/2^k, k, 0, oo)", "2"),
        t(Closed(1.0), "sum(1/(k*(k + 1)), k, 1, oo)", "1"),
        t(Closed(55.0), "sum(k, k, 1, 10)", "55"),
        t(Closed(385.0), "sum(k^2, k, 1, 10)", "385"),
        t(Closed(3025.0), "sum(k^3, k, 1, 10)", "3025"),
        t(Closed(1023.0), "sum(2^k, k, 0, 9)", "1023"),
        t(Closed(e), "sum(1/factorial(k), k, 0, oo)", "E"),
        t(Closed(1.5), "sum(1/3^k, k, 0, oo)", "3/2"),
        t(
            Closed(pi.powi(4) / 90.0),
            "sum(1/k^4, k, 1, oo)",
            "1/90*pi^4",
        ),
        t(Closed(0.75), "sum(1/((k)*(k + 2)), k, 1, oo)", "3/4"),
        t(Closed(2.0), "sum(k/2^k, k, 1, oo)", "2"),
        // ---- series (compared near 0)
        c(Series(6), "taylor(sin(x)/x, x, 0, 6)"),
        c(Series(7), "taylor(sin(x), x, 0, 7)"),
        c(Series(6), "taylor(exp(x), x, 0, 6)"),
        c(Series(6), "taylor(cos(x), x, 0, 6)"),
        c(Series(6), "taylor(1/(1 - x), x, 0, 6)"),
        c(Series(6), "taylor(ln(1 + x), x, 0, 6)"),
        c(Series(6), "taylor(atan(x), x, 0, 6)"),
        c(Series(6), "taylor(1/(1 + x^2), x, 0, 6)"),
        c(Series(6), "taylor(exp(x)*cos(x), x, 0, 6)"),
        c(Series(6), "taylor(sqrt(1 + x), x, 0, 6)"),
        c(Series(6), "taylor(tan(x), x, 0, 6)"),
        // ---- solve
        c(Solve, "solve(x^2 - 4, x)"),
        c(Solve, "solve(x^2 - 2, x)"),
        c(SolveNone, "solve(x^2 + 1, x)"),
        c(Solve, "solve(x^3 - 6*x^2 + 11*x - 6, x)"),
        c(Solve, "solve(x^3 - 2, x)"),
        c(Solve, "solve(x^4 - 5*x^2 + 4, x)"),
        c(Solve, "solve(x^4 - 1, x)"),
        c(Solve, "solve(x^4 + x^2 - 2, x)"),
        c(Solve, "solve(2*x + 3, x)"),
        c(Solve, "solve(x^2 - 2*x + 1, x)"),
        c(Solve, "solve(exp(x) - 2, x)"),
        c(Solve, "solve(ln(x) - 1, x)"),
        c(Solve, "solve(sin(x) - 1/2, x)"),
        c(Solve, "solve(cos(x), x)"),
        c(Solve, "solve(2^x - 8, x)"),
        c(Solve, "solve(sqrt(x) - 3, x)"),
        c(Solve, "solve(x^2 - 3*x + 2, x)"),
        c(Solve, "solve(x^3 + x, x)"),
        c(Solve, "solve(tan(x) - 1, x)"),
        c(Solve, "solve(exp(2*x) - 3*exp(x) + 2, x)"),
        // ================= batch 2: harder families =================
        // ---- rational functions with repeated / irreducible factors
        c(Integral, "integral(1/(x^2 + 1)^2, x)"),
        c(Integral, "integral(x/(x^2 + 1)^2, x)"),
        c(Integral, "integral(x^2/(x^2 + 1)^2, x)"),
        c(Integral, "integral(1/(x^3 - x), x)"),
        c(Integral, "integral((2*x + 3)/(x^2 + 4*x + 5), x)"),
        c(Integral, "integral(1/(x^2 - 2*x + 5), x)"),
        c(Integral, "integral(x^3/(x^2 + 1), x)"),
        c(Integral, "integral(1/(x^2*(x - 1)), x)"),
        c(Integral, "integral(1/(x^3*(x + 1)), x)"),
        c(Integral, "integral((x + 2)/(x^2 - 4), x)"),
        c(Integral, "integral((x^2 + 1)/(x^4 + 1), x)"),
        c(Integral, "integral(1/(x^4 + 5*x^2 + 4), x)"),
        c(Integral, "integral((x^5 + 1)/(x^2 - 1), x)"),
        c(Integral, "integral(1/((x - 1)*(x - 2)*(x - 3)), x)"),
        c(Integral, "integral(1/(x^2 + x + 1), x)"),
        // ---- trigonometric
        c(Integral, "integral(sin(x)^2*cos(x), x)"),
        c(Integral, "integral(sin(x)*cos(x)^2, x)"),
        c(Integral, "integral(sin(x)^5, x)"),
        c(Integral, "integral(cos(x)^4, x)"),
        c(Integral, "integral(tan(x)^3, x)"),
        c(Integral, "integral(1/cos(x)^3, x)"),
        c(Integral, "integral(sin(3*x), x)"),
        c(Integral, "integral(sin(x)*cos(2*x), x)"),
        c(Integral, "integral(x*sin(x)^2, x)"),
        c(Integral, "integral(sin(x)/cos(x)^2, x)"),
        c(Integral, "integral(cos(x)/(1 + sin(x)^2), x)"),
        c(Integral, "integral(1/(sin(x)*cos(x)), x)"),
        c(Integral, "integral(1/(1 + sin(x)), x)"),
        c(Integral, "integral(1/(3 + 2*cos(x)), x)"),
        c(Integral, "integral(sin(x)^2*cos(x)^2, x)"),
        c(Integral, "integral(cot(x), x)"),
        c(Integral, "integral(1/sin(x)^2, x)"),
        c(Integral, "integral(sec(x)*tan(x), x)"),
        // ---- exponential / logarithmic / parts
        c(Integral, "integral(exp(-x)*cos(2*x), x)"),
        c(Integral, "integral(x*exp(2*x), x)"),
        c(Integral, "integral(x^3*exp(x), x)"),
        c(Integral, "integral(x^2*exp(-x), x)"),
        c(Integral, "integral(sqrt(x)*ln(x), x)"),
        c(Integral, "integral(x^2*atan(x), x)"),
        c(Integral, "integral(ln(x^2 + 1), x)"),
        c(Integral, "integral(x*ln(x^2 + 1), x)"),
        c(Integral, "integral(atan(x)/(1 + x^2), x)"),
        c(Integral, "integral(asin(x)/sqrt(1 - x^2), x)"),
        c(Integral, "integral(ln(x)/x^2, x)"),
        c(Integral, "integral(x*acos(x), x)"),
        c(Integral, "integral(exp(x)*sin(x)*cos(x), x)"),
        c(Integral, "integral(exp(2*x)/(1 + exp(x)), x)"),
        c(Integral, "integral(1/(exp(x) + exp(-x)), x)"),
        c(Integral, "integral(sech(x)^2, x)"),
        c(Integral, "integral(tanh(x), x)"),
        c(Integral, "integral(x*sinh(x), x)"),
        c(Integral, "integral(exp(sqrt(x)), x)"),
        c(Integral, "integral(ln(ln(x))/x, x)"),
        c(Integral, "integral(x*exp(-x^2/2), x)"),
        c(Integral, "integral(2^x, x)"),
        c(Integral, "integral(x*2^x, x)"),
        c(Integral, "integral(x^x*(1 + ln(x)), x)"),
        // ---- substitution / algebraic
        c(Integral, "integral(x*sqrt(x^2 + 1), x)"),
        c(Integral, "integral(x^2*sqrt(x^3 + 1), x)"),
        c(Integral, "integral(cos(x)*exp(sin(x)), x)"),
        c(Integral, "integral(exp(x)/(1 + exp(2*x)), x)"),
        c(Integral, "integral(1/(x*sqrt(x^2 - 1)), x)"),
        c(Integral, "integral(1/(1 + sqrt(x)), x)"),
        c(Integral, "integral(1/(sqrt(x) + x), x)"),
        c(Integral, "integral(sqrt(1 + x)/x, x)"),
        c(Integral, "integral(1/(x*sqrt(1 - x^2)), x)"),
        c(Integral, "integral(sqrt(x^2 - 1), x)"),
        c(Integral, "integral(1/(x^2*sqrt(x^2 + 1)), x)"),
        c(Integral, "integral(sqrt(4 - x^2), x)"),
        c(Integral, "integral(sqrt(x)/(1 + x), x)"),
        c(Integral, "integral(x/sqrt(1 - x^2), x)"),
        c(Integral, "integral(1/sqrt(2*x - x^2), x)"),
        c(Integral, "integral(abs(x), x)"),
        // ---- differentiation (chain, product, quotient, inverse functions)
        c(Diff, "diff(x^2*sin(x)*exp(x), x)"),
        c(Diff, "diff(sin(cos(x)), x)"),
        c(Diff, "diff(ln(x + sqrt(x^2 + 1)), x)"),
        c(Diff, "diff(atan(x^2), x)"),
        c(Diff, "diff(x^(1/3), x)"),
        c(Diff, "diff((x^2 + 1)^(3/2), x)"),
        c(Diff, "diff(exp(x)/(1 + x^2), x)"),
        c(Diff, "diff(sin(x)^cos(x), x)"),
        c(Diff, "diff(asinh(x), x)"),
        c(Diff, "diff(atanh(x), x)"),
        c(Diff, "diff(acos(x)*x, x)"),
        c(Diff, "diff(ln(abs(x)), x)"),
        c(Diff, "diff(sqrt(x)*ln(x), x)"),
        c(Diff, "diff(tanh(x)^2, x)"),
        c(Diff, "diff(sec(x), x)"),
        c(Diff, "diff(1/(x^2 + 1)^2, x)"),
        c(Diff, "diff(exp(exp(x)), x)"),
        c(Diff, "diff(x*exp(-x^2), x)"),
        c(Diff, "diff(si(x), x)"),
        c(Diff, "diff(erfc(x), x)"),
        c(Diff, "diff(besselj(0, x), x)"),
        c(Diff, "diff(lgamma(x), x)"),
        // ---- simplification
        t(Simplify, "exp(-ln(x))", "1/x"),
        t(Simplify, "exp(ln(x) + ln(y))", "x*y"),
        t(Simplify, "ln(exp(2))", "2"),
        t(Simplify, "exp(ln(2) + ln(3))", "6"),
        t(Simplify, "sin(pi/6)", "1/2"),
        t(Simplify, "cos(pi/3)", "1/2"),
        t(Simplify, "sin(x + pi/2)", "cos(x)"),
        t(Simplify, "cos(x + pi)", "-cos(x)"),
        t(Simplify, "sin(x + 2*pi)", "sin(x)"),
        t(Simplify, "sin(pi - x)", "sin(x)"),
        c(Simplify, "cos(2*x) + 2*sin(x)^2"),
        c(Simplify, "tan(x)^2 + 1"),
        c(Simplify, "1/cos(x)^2 - tan(x)^2"),
        c(Simplify, "(a + b)*(a - b)"),
        c(Simplify, "(x + 1)^3 - x^3 - 3*x^2 - 3*x - 1"),
        c(Simplify, "(x^2 - 1)/(x^2 + 2*x + 1)"),
        c(Simplify, "x/(x^2 + x)"),
        c(Simplify, "1/x + 1/(x + 1)"),
        c(Simplify, "x/(x + 1) - 1/(x + 1)"),
        c(Simplify, "sqrt(x)^2"),
        c(Simplify, "(x^(1/2))^4"),
        c(Simplify, "abs(x)^2"),
        c(Simplify, "exp(x)*exp(-x)"),
        c(Simplify, "ln(x) - ln(x)"),
        c(Simplify, "2^x*2^x"),
        c(Simplify, "(2^x)^2"),
        c(Simplify, "sinh(x) + cosh(x)"),
        c(Simplify, "(exp(x) - exp(-x))/2"),
        c(Simplify, "(exp(x) + exp(-x))/2"),
        c(Simplify, "sin(x)^4 - cos(x)^4"),
        c(Simplify, "x^3/x"),
        c(Simplify, "x^a*x^b"),
        c(Simplify, "(x^a)^b"),
        t(Value(3.0), "log(2, 8)", "3"),
        t(Value(2.0), "log(10, 100)", "2"),
        t(Value(0.0), "atan(0)", "0"),
        t(Value(pi / 4.0), "atan(1)", "1/4*pi"),
        t(Value(pi / 3.0), "acos(1/2)", "1/3*pi"),
        t(Value(1.0), "cosh(0)", "1"),
        t(Value(0.0), "sinh(0)", "0"),
        t(Value(-1.0), "cos(pi)", "-1"),
        // ---- special values
        t(Value(1.0), "gamma(2)", "1"),
        t(Value(6.0), "gamma(4)", "6"),
        t(Value(0.75 * pi.sqrt()), "gamma(5/2)", "3/4*pi^(1/2)"),
        t(Value(-2.0 * pi.sqrt()), "gamma(-1/2)", "-2*pi^(1/2)"),
        t(Value(-0.5), "zeta(0)", "-1/2"),
        t(Value(-1.0 / 12.0), "zeta(-1)", "-1/12"),
        t(Value(1.0 / 120.0), "zeta(-3)", "1/120"),
        t(Value(0.0), "zeta(-2)", "0"),
        t(Value(pi.powi(6) / 945.0), "zeta(6)", "1/945*pi^6"),
        t(Value(0.0), "shi(0)", "0"),
        t(Value(0.0), "erfi(0)", "0"),
        t(Value(1.0), "sinc(0)", "1"),
        t(Value(1.0), "besselj(0, 0)", "1"),
        t(Value(0.0), "besselj(1, 0)", "0"),
        // ---- limits (more)
        t(Value(e * e), "limit((1 + 2/x)^x, x, oo)", "exp(2)"),
        t(Value(e), "limit(x^(1/ln(x)), x, oo)", "E"),
        t(Value(0.0), "limit(ln(x)/sqrt(x), x, oo)", "0"),
        t(Value(f64::INFINITY), "limit(exp(1/x), x, 0, plus)", "oo"),
        t(Value(0.0), "limit(exp(-1/x^2), x, 0)", "0"),
        t(Value(0.0), "limit(exp(1/x), x, 0, minus)", "0"),
        t(Value(0.5), "limit(x*(sqrt(x^2 + 1) - x), x, oo)", "1/2"),
        t(Value(-pi / 2.0), "limit(atan(x), x, -oo)", "-1/2*pi"),
        t(Value(1.0), "limit(tanh(x), x, oo)", "1"),
        t(Value(1.0 / 6.0), "limit((x - sin(x))/x^3, x, 0)", "1/6"),
        t(Value(2.0_f64.ln()), "limit((2^x - 1)/x, x, 0)", "ln(2)"),
        t(
            Value(1.0 / 3.0),
            "limit((sin(x) - x*cos(x))/x^3, x, 0)",
            "1/3",
        ),
        t(Value(1.0), "limit(ln(1 + x)/x, x, 0)", "1"),
        t(Value(0.0), "limit(x^2/exp(x), x, oo)", "0"),
        t(Value(-0.5), "limit((cos(x) - 1)/x^2, x, 0)", "-1/2"),
        t(Value(1.0), "limit(x/sin(x), x, 0)", "1"),
        t(Value(3.0), "limit((3*x^2 + x)/(x^2 + 5), x, oo)", "3"),
        t(Value(0.0), "limit((x + 1)/(x^2 + 1), x, oo)", "0"),
        t(
            Value(f64::INFINITY),
            "limit((x^3 + 1)/(x^2 + 1), x, oo)",
            "oo",
        ),
        t(Value(f64::NEG_INFINITY), "limit(-x^2, x, oo)", "-oo"),
        t(Value(0.0), "limit(x*exp(-x^2), x, oo)", "0"),
        t(Value(1.0), "limit(exp(x)/(exp(x) + 1), x, oo)", "1"),
        t(Value(0.0), "limit(exp(x)/(exp(x) + 1), x, -oo)", "0"),
        t(Value(0.5), "limit(1/(1 + exp(-x)), x, 0)", "1/2"),
        t(Value(1.0), "limit(x^(x^2), x, 1)", "1"),
        t(Value(0.0), "limit(sin(x)*exp(-x), x, oo)", "0"),
        t(Value(-1.0), "limit(sin(x)/(x - pi), x, pi)", "-1"),
        t(Value(0.25), "limit((sqrt(x) - 2)/(x - 4), x, 4)", "1/4"),
        c(Undefined, "limit(1/x, x, 0)"),
        c(Undefined, "limit(abs(x)/x, x, 0)"),
        c(Undefined, "limit(sin(x), x, oo)"),
        c(Undefined, "limit(sin(1/x), x, 0)"),
        // ---- sums
        t(
            Closed(2.0_f64.ln()),
            "sum((-1)^(k + 1)/k, k, 1, oo)",
            "ln(2)",
        ),
        t(Closed(0.25), "sum(1/(k*(k + 1)*(k + 2)), k, 1, oo)", "1/4"),
        t(Closed(0.5), "sum(1/((2*k - 1)*(2*k + 1)), k, 1, oo)", "1/2"),
        t(
            Closed(pi * pi / 12.0),
            "sum((-1)^(k + 1)/k^2, k, 1, oo)",
            "1/12*pi^2",
        ),
        t(
            Closed(1.202_056_903_159_594),
            "sum(1/k^3, k, 1, oo)",
            "zeta(3)",
        ),
        t(Closed(1.0), "sum(1/2^k, k, 1, oo)", "1"),
        t(Closed(6.0), "sum(k^2/2^k, k, 1, oo)", "6"),
        t(
            Closed(e.recip()),
            "sum((-1)^k/factorial(k), k, 0, oo)",
            "exp(-1)",
        ),
        t(Closed(10.0), "sum(1, k, 1, 10)", "10"),
        t(Closed(5050.0), "sum(k, k, 1, 100)", "5050"),
        t(
            Closed(1.0 - 1.0 / 1024.0),
            "sum(1/2^k, k, 1, 10)",
            "1023/1024",
        ),
        c(Undefined, "sum(1/k, k, 1, oo)"),
        c(Undefined, "sum(k, k, 1, oo)"),
        // ---- series
        c(Series(5), "taylor(1/cos(x), x, 0, 6)"),
        c(Series(5), "taylor(exp(-x^2), x, 0, 6)"),
        c(Series(5), "taylor(asin(x), x, 0, 5)"),
        c(Series(5), "taylor(sinh(x), x, 0, 5)"),
        c(Series(5), "taylor(ln(1 + x)/(1 + x), x, 0, 5)"),
        c(Series(5), "taylor((1 + x)^(1/3), x, 0, 5)"),
        c(Series(5), "taylor(exp(sin(x)), x, 0, 5)"),
        c(Series(5), "taylor(x/(exp(x) - 1), x, 0, 5)"),
        c(Series(5), "taylor(sqrt(1 - x^2), x, 0, 5)"),
        c(Series(5), "taylor(cos(x)^2, x, 0, 5)"),
        // ---- solve
        c(Solve, "solve(x^3 - 3*x + 1, x)"),
        c(Solve, "solve(x^4 - 10*x^2 + 1, x)"),
        c(Solve, "solve(x^3 - x^2 - 2*x + 2, x)"),
        c(Solve, "solve(x^4 - 16, x)"),
        c(SolveNone, "solve(x^4 + 4, x)"),
        c(Solve, "solve(x^2 - 5*x + 6, x)"),
        c(Solve, "solve(x^3 + 3*x^2 + 3*x + 1, x)"),
        c(Solve, "solve(4*x^3 - 3*x - 1/2, x)"),
        c(Solve, "solve(2*cos(x) - 1, x)"),
        c(Solve, "solve(sin(x), x)"),
        c(Solve, "solve(tan(x) - sqrt(3), x)"),
        c(Solve, "solve(exp(x) + exp(-x) - 3, x)"),
        c(Solve, "solve(ln(x) + ln(x - 1) - ln(2), x)"),
        c(Solve, "solve(abs(x - 1) - 2, x)"),
        c(Solve, "solve(sqrt(x + 1) - x + 1, x)"),
        c(Solve, "solve(1/x + 1/(x + 1) - 1, x)"),
        c(Solve, "solve(x^2 - 2*x - 1, x)"),
        c(Solve, "solve(3^x - 9, x)"),
        c(Solve, "solve(x*exp(x) - 1, x)"),
        c(Solve, "solve(sin(x)^2 - 1/4, x)"),
        c(Solve, "solve(x^5 - x, x)"),
        c(Solve, "solve(x^6 - 1, x)"),
        c(Solve, "solve(2*x^2 - 8, x)"),
        c(Solve, "solve(x^2 + 2*x + 1, x)"),
    ]
}

/// Strips the request wrapper so `solve(eq, x)` yields `eq` with the free
/// symbol `x` left in place.
#[test]
fn battery() {
    let all = cases();
    let mut table = String::from("| # | status | request | answer |\n|---|---|---|---|\n");
    let mut bad: Vec<String> = Vec::new();
    for (i, case) in all.iter().enumerate() {
        let (status, answer) = std::panic::catch_unwind(|| judge(case))
            .unwrap_or_else(|_| ("PANIC".to_owned(), String::new()));
        let _ = writeln!(
            table,
            "| {} | {} | `{}` | `{}` |",
            i + 1,
            status,
            case.src,
            answer.replace('|', "\\|")
        );
        if status != "ok" {
            bad.push(format!("{status:12} {}  =>  {answer}", case.src));
        }
    }
    if let Ok(path) = std::env::var("AUDIT_TABLE") {
        let _ = std::fs::write(path, &table);
    }
    println!("{} cases, {} not ok", all.len(), bad.len());
    for b in &bad {
        println!("{b}");
    }
    let wrong: Vec<&String> = bad
        .iter()
        .filter(|b| b.starts_with("WRONG") || b.starts_with("FLOAT") || b.starts_with("PANIC"))
        .collect();
    assert!(wrong.is_empty(), "unsound or float answers: {wrong:#?}");
}
