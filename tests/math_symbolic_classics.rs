//! A second table of classical symbolic results (definite integrals from the
//! standard tables, limits, closed-form sums, Maclaurin series, equation
//! solving, polynomial and number-theory identities). Every answer is
//! checked against an independent oracle: a known constant, a numeric
//! evaluation at sample points, or substitution back into the equation.
//!
//! NOTE: `limit((1 + 1/x)^(x^2)/exp(x), x, oo)` currently answers `0`; the
//! true value is `exp(-1/2)`. That is a wrong answer, not a missing one.
//!
//! Cases that rssn cannot do yet are listed with `known_gap` and must keep
//! failing; the test fails when a gap closes so that it is moved into the
//! passing set. `CLASSICS_REPORT=1 cargo test --test math_symbolic_classics
//! -- --nocapture` prints every answer.

mod common;
use std::f64::consts::E;
use std::f64::consts::PI;

use common::*;

#[derive(Clone, Copy, Debug)]
enum Kind {
    /// The answer has no free symbols and equals the number.
    Value(f64),
    /// The answer, evaluated at sample points, equals the input.
    Same,
    /// The antiderivative differentiates back to the integrand.
    Integral,
    /// Every listed root satisfies `f(x) = 0` and at least `n` were found.
    Solve(usize),
    /// A Maclaurin polynomial of the given order matches near 0.
    Series(u32),
    /// The answer prints exactly like this.
    Text(&'static str),
}

struct Case {
    kind: Kind,
    src: &'static str,
    known_gap: bool,
}

const fn v(
    value: f64,
    src: &'static str,
) -> Case {
    Case {
        kind: Kind::Value(value),
        src,
        known_gap: false,
    }
}
const fn same(src: &'static str) -> Case {
    Case {
        kind: Kind::Same,
        src,
        known_gap: false,
    }
}
const fn int(src: &'static str) -> Case {
    Case {
        kind: Kind::Integral,
        src,
        known_gap: false,
    }
}
const fn solve(
    n: usize,
    src: &'static str,
) -> Case {
    Case {
        kind: Kind::Solve(n),
        src,
        known_gap: false,
    }
}
const fn series(
    n: u32,
    src: &'static str,
) -> Case {
    Case {
        kind: Kind::Series(n),
        src,
        known_gap: false,
    }
}
const fn text(
    text: &'static str,
    src: &'static str,
) -> Case {
    Case {
        kind: Kind::Text(text),
        src,
        known_gap: false,
    }
}
impl Case {
    const fn gap(mut self) -> Self {
        self.known_gap = true;
        self
    }
}

fn first_arg(src: &str) -> String {
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

fn roots(text: &str) -> Vec<String> {
    let inner = text
        .trim()
        .strip_prefix("list(")
        .and_then(|s| s.strip_suffix(')'))
        .unwrap_or(text.trim());
    let (mut out, mut cur, mut depth) = (Vec::new(), String::new(), 0_i32);
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
            | ',' if depth == 0 => out.push(std::mem::take(&mut cur).trim().to_owned()),
            | _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_owned());
    }
    out
}

fn close(
    a: f64,
    b: f64,
    tol: f64,
) -> bool {
    if a == b {
        return true;
    }
    (a - b).abs() <= tol * (1.0 + b.abs())
}

/// `Ok` when the case passes, `Err(reason)` otherwise.
fn judge(case: &Case) -> (Result<(), String>, String) {
    let o = sym(case.src);
    let shown = o.text.clone();
    let verdict = (|| {
        if has_float(&o.text) {
            return Err("float literal in an exact answer".to_owned());
        }
        if !o.reduced {
            return Err("not reduced".to_owned());
        }
        match case.kind {
            | Kind::Value(want) => {
                match eval(&o.text, &[]) {
                    | Some(got) if close(got, want, 1e-8) || (want.is_infinite() && got == want) => {
                        Ok(())
                    },
                    | Some(got) => Err(format!("value {got} != {want}")),
                    | None => Err("not a number".to_owned()),
                }
            },
            | Kind::Same => {
                let mut checked = 0;
                for p in [0.31, 0.57, 0.83, 1.37, 2.11] {
                    let b = [("x", p), ("y", 0.7)];
                    let (Some(a), Some(g)) = (eval(case.src, &b), eval(&o.text, &b)) else {
                        continue;
                    };
                    if a.is_finite() && g.is_finite() {
                        checked += 1;
                        if !close(g, a, 1e-9) {
                            return Err(format!("at x={p}: {g} != {a}"));
                        }
                    }
                }
                if checked == 0 {
                    Err("unchecked".to_owned())
                } else {
                    Ok(())
                }
            },
            | Kind::Integral => {
                if o.text.contains("integral(") {
                    return Err("unevaluated".to_owned());
                }
                let f = first_arg(case.src);
                let h = 1e-5;
                let mut checked = 0;
                for p in [0.31, 0.57, 0.83, 1.37, 2.11] {
                    let (Some(up), Some(dn), Some(want)) = (
                        eval(&o.text, &[("x", p + h)]),
                        eval(&o.text, &[("x", p - h)]),
                        eval(&f, &[("x", p)]),
                    ) else {
                        continue;
                    };
                    let d = (up - dn) / (2.0 * h);
                    if d.is_finite() && want.is_finite() {
                        checked += 1;
                        if !close(d, want, 1e-5) {
                            return Err(format!("d/dx at {p}: {d} != {want}"));
                        }
                    }
                }
                if checked == 0 {
                    Err("unchecked".to_owned())
                } else {
                    Ok(())
                }
            },
            | Kind::Solve(n) => {
                let eq = first_arg(case.src);
                let rs = roots(&o.text);
                if rs.len() < n {
                    return Err(format!("{} roots, wanted {n}", rs.len()));
                }
                for r in rs {
                    let Some(rv) = eval(&r, &[]) else {
                        return Err(format!("root `{r}` is not numeric"));
                    };
                    match eval(&eq, &[("x", rv)]) {
                        | Some(res) if res.abs() < 1e-7 => {},
                        | Some(res) => return Err(format!("residual {res} at {rv}")),
                        | None => return Err(format!("cannot evaluate residual at `{r}`")),
                    }
                }
                Ok(())
            },
            | Kind::Series(n) => {
                let f = first_arg(case.src);
                for p in [0.05, 0.1, -0.08] {
                    let (Some(a), Some(b)) = (eval(&f, &[("x", p)]), eval(&o.text, &[("x", p)]))
                    else {
                        continue;
                    };
                    let tol = 50.0 * f64::powi(p.abs(), i32::try_from(n).unwrap_or(1));
                    if (a - b).abs() > tol.max(1e-9) {
                        return Err(format!("at {p}: {b} vs {a}"));
                    }
                }
                Ok(())
            },
            | Kind::Text(want) => {
                if o.text == want {
                    Ok(())
                } else {
                    Err(format!("text `{}` != `{want}`", o.text))
                }
            },
        }
    })();
    (verdict, shown)
}

fn cases() -> Vec<Case> {
    let ln2 = std::f64::consts::LN_2;
    vec![
        // ---- definite integrals (standard tables)
        v(PI / 2.0, "defint(sin(x)^2, x, 0, pi)"),
        v(1.0, "defint(x*exp(-x), x, 0, oo)"),
        v(2.0, "defint(x^2*exp(-x), x, 0, oo)"),
        v(24.0, "defint(x^4*exp(-x), x, 0, oo)"),
        v(0.5, "defint(x*exp(-x^2), x, 0, oo)"),
        v(0.5, "defint(x^3*exp(-x^2), x, 0, oo)"),
        v(PI / 2.0, "defint(1/(1 + x^2)^2, x, -oo, oo)"),
        v(PI, "defint(1/(1 + x^2), x, -oo, oo)"),
        v(PI, "defint(1/(x^2 + 2*x + 2), x, -oo, oo)"),
        v(PI / 4.0, "defint(x^2/(1 + x^2)^2, x, 0, oo)"),
        v(PI * PI / 6.0, "defint(x/(exp(x) - 1), x, 0, oo)").gap(),
        v(-PI / 2.0 * ln2, "defint(ln(sin(x)), x, 0, pi/2)").gap(),
        v(2.0, "defint(1/sqrt(x), x, 0, 1)"),
        v(PI, "defint(1/sqrt(1 - x^2), x, -1, 1)"),
        v(PI / 2.0, "defint(sqrt(1 - x^2), x, -1, 1)"),
        v(PI, "defint(x*sin(x), x, 0, pi)"),
        v(PI, "defint(cos(x)^2, x, 0, 2*pi)"),
        v(6.5, "defint(abs(x), x, -2, 3)"),
        v(E - 1.0, "defint(exp(x), x, 0, 1)"),
        v(1.0, "defint(1/x, x, 1, E)"),
        v(2.0, "defint(ln(x)^2, x, 0, 1)"),
        v(0.25, "defint(x*ln(x)^2, x, 0, 1)"),
        v(0.5, "defint(sin(x)*exp(-x), x, 0, oo)"),
        v(0.5, "defint(cos(x)*exp(-x), x, 0, oo)"),
        v(4.0 / 3.0, "defint(sin(x)^3, x, 0, pi)"),
        v(1.0, "defint(1/cos(x)^2, x, 0, pi/4)"),
        v(PI / 4.0, "defint(sin(x)^2*cos(x)^2, x, 0, 2*pi)"),
        v((2.0 * PI).sqrt(), "defint(exp(-x^2/2), x, -oo, oo)"),
        v(PI.sqrt() / 2.0, "defint(sqrt(x)*exp(-x), x, 0, oo)"),
        v(PI * PI / 12.0, "defint(ln(1 + x)/x, x, 0, 1)").gap(),
        v(PI, "defint(1/sqrt(x*(1 - x)), x, 0, 1)"),
        v(
            2.0 * PI / (3.0 * 3f64.sqrt()),
            "defint(1/(1 + x^3), x, 0, oo)",
        )
        .gap(),
        v(PI / 2.0, "defint(sin(x)/x, x, 0, oo)"),
        v(0.0, "defint(x^3, x, -1, 1)"),
        v(2.0 / 3.0, "defint(sqrt(x), x, 0, 1)"),
        v(4.0, "defint(x, x, 0, 2) + defint(x, x, 0, 2)"),
        // ---- limits
        v(-1.0 / 6.0, "limit((sin(x) - x)/x^3, x, 0)"),
        v(0.0, "limit(x^2*sin(1/x), x, 0)"),
        v(ln2, "limit((2^x - 1)/x, x, 0)"),
        v(1.0, "limit((exp(x) - exp(-x))/(2*x), x, 0)"),
        v(0.6, "limit(sin(3*x)/sin(5*x), x, 0)"),
        v((-0.5f64).exp(), "limit((1 + 1/x)^(x^2)/exp(x), x, oo)").gap(),
        v(1.5, "limit(sqrt(x^2 + 3*x) - x, x, oo)"),
        v(1.0, "limit(x*(exp(1/x) - 1), x, oo)"),
        v(0.5, "limit((tan(x) - sin(x))/x^3, x, 0)"),
        v(1.0, "limit(asin(x)/x, x, 0)"),
        v(1.0, "limit(abs(x)/x, x, 0, plus)"),
        v(-1.0, "limit(abs(x)/x, x, 0, minus)"),
        v(f64::INFINITY, "limit(1/x, x, 0, plus)"),
        v(f64::NEG_INFINITY, "limit(ln(x), x, 0, plus)"),
        v(f64::INFINITY, "limit(x^x, x, oo)"),
        v(0.0, "limit(exp(-x)*x^5, x, oo)"),
        v(2.0, "limit((1 - cos(2*x))/x^2, x, 0)"),
        v(3.0, "limit((x^3 - 1)/(x - 1), x, 1)"),
        v(0.0, "limit(ln(x)^2/x, x, oo)"),
        v(0.0, "limit((x + 1)^(1/3) - x^(1/3), x, oo)"),
        v((-0.5f64).exp(), "limit(cos(x)^(1/x^2), x, 0)"),
        v(1.0, "limit(x^(sin(x)), x, 0, plus)"),
        v(0.0, "limit(x^2/exp(x), x, oo)"),
        v(1.0, "limit((x^2 + x)/(x^2 - x), x, oo)"),
        v(0.5, "limit((sqrt(1 + x) - 1)/x, x, 0)"),
        v(2.0, "limit(x*sin(2/x), x, oo)"),
        v(-1.0, "limit(cos(x)/x, x, -oo) - 1"),
        // ---- sums and products
        v(0.75, "sum(1/(k^2 - 1), k, 2, oo)"),
        v(-ln2, "sum((-1)^k/k, k, 1, oo)"),
        v(ln2, "sum(1/(2^k*k), k, 1, oo)").gap(),
        v(6.0, "sum(k^2/2^k, k, 1, oo)"),
        v(PI / 4.0, "sum((-1)^k/(2*k + 1), k, 0, oo)").gap(),
        v(0.5, "sum(1/(4*k^2 - 1), k, 1, oo)"),
        v(PI.powi(6) / 945.0, "sum(1/k^6, k, 1, oo)"),
        v(0.25, "sum(1/(k*(k + 1)*(k + 2)), k, 1, oo)"),
        v(3080.0, "sum(k*(k + 1), k, 1, 20)"),
        v(5269.0 / 3600.0, "sum(1/k^2, k, 1, 5)"),
        v(1024.0, "sum(binomial(10, k), k, 0, 10)"),
        v(5120.0, "sum(k*binomial(10, k), k, 0, 10)"),
        v(3f64.exp(), "sum(3^k/factorial(k), k, 0, oo)"),
        v((-1f64).exp(), "sum((-1)^k/factorial(k), k, 0, oo)"),
        v(0.5, "sum(1/3^k, k, 1, oo)"),
        v(f64::INFINITY, "sum(1/k, k, 1, oo)").gap(),
        v(720.0, "product(k, k, 1, 6)"),
        // ---- Maclaurin series
        series(7, "taylor(tan(x), x, 0, 7)"),
        series(7, "taylor(asin(x), x, 0, 7)"),
        series(7, "taylor(atan(x), x, 0, 7)"),
        series(6, "taylor(exp(sin(x)), x, 0, 6)"),
        series(5, "taylor(sqrt(1 + x), x, 0, 5)"),
        series(5, "taylor((1 + x)^(1/3), x, 0, 5)"),
        series(6, "taylor(ln(cos(x)), x, 0, 6)"),
        series(6, "taylor(1/cos(x), x, 0, 6)"),
        series(6, "taylor(x/(exp(x) - 1), x, 0, 6)"),
        series(6, "taylor(exp(x)*cos(x), x, 0, 6)"),
        series(6, "taylor(cosh(x), x, 0, 6)"),
        series(8, "taylor(1/(1 + x^2), x, 0, 8)"),
        series(7, "taylor(asinh(x), x, 0, 7)"),
        series(7, "taylor(erf(x), x, 0, 7)"),
        series(6, "taylor(ln(1 + x)/(1 + x), x, 0, 6)"),
        // ---- equations
        solve(3, "solve(x^3 - 6*x^2 + 11*x - 6, x)"),
        solve(4, "solve(x^4 - 5*x^2 + 4, x)"),
        solve(1, "solve(2^x - 8, x)"),
        solve(1, "solve(ln(x) - 2, x)"),
        solve(1, "solve(sin(x) - 1/2, x)"),
        solve(2, "solve(x^2 - 2*x - 3, x)"),
        solve(1, "solve(sqrt(x + 2) - x, x)"),
        solve(2, "solve(abs(x - 1) - 3, x)"),
        solve(1, "solve(x^3 - 2, x)"),
        solve(3, "solve(x^5 - x, x)"),
        solve(2, "solve(1/x + 1/(x + 1) - 1, x)"),
        solve(2, "solve(exp(2*x) - 3*exp(x) + 2, x)"),
        solve(1, "solve(x*exp(x) - 1, x)"),
        solve(1, "solve(x^2 - 2*x + 1, x)"),
        solve(2, "solve(x^4 - 16, x)"),
        solve(1, "solve(cos(x) - x, x)").gap(),
        solve(2, "solve(x^2 - 3, x)"),
        // ---- algebra
        same("expand((x + y)^5)"),
        same("factor(x^6 - 1)"),
        same("factor((x^2 - 1)*(x^2 + x + 1))"),
        same("apart(1/(x^3 - x), x)"),
        same("apart((x^2 + 1)/((x - 1)^2*(x + 2)), x)"),
        same("cancel((x^3 - 1)/(x - 1))"),
        same("together(1/x + 1/y)"),
        same("cos(3*x) - 4*cos(x)^3 + 3*cos(x)"),
        same("sin(2*x) - 2*sin(x)*cos(x)"),
        same("(x^3 + y^3)/(x + y)"),
        same("sqrt(x^2)"),
        same("(x^4 - y^4)/(x^2 + y^2)"),
        same("tan(2*x) - 2*tan(x)/(1 - tan(x)^2)"),
        same("cosh(x)^2 - sinh(x)^2"),
        same("ln(x*y) - ln(x) - ln(y)"),
        v(5.0, "discriminant(x^3 - 2*x + 1, x)"),
        v(3.0, "resultant(x^2 - 1, x - 2, x)"),
        // ---- number theory
        text("true", "isprime(2^61 - 1)"),
        text("false", "isprime(2^61 + 1)"),
        text("true", "isprime(18446744073709551557)"),
        text("false", "isprime(3825123056546413051)"),
        text("4", "invmod(3, 11)"),
        text("1000000000039", "nextprime(10^12)"),
        text("168", "prime_pi(1000)"),
        text("541", "nth_prime(100)"),
        text("-1", "mobius(30)"),
        text("5", "primitive_root(23)"),
        text("1000000000000000", "isqrt(10^30)"),
        v(4000000000.0, "totient(10^10)"),
        v(1.0, "powmod(3, 1000000006, 1000000007)"),
        v(12.0, "gcd(48, 180)"),
        v(720.0, "lcm(48, 180)"),
        v(1.0, "jacobi(2, 7)"),
        v(4.0, "dlog(2, 16, 101)"),
        // ---- special function identities at exact points
        v(PI.sqrt(), "gamma(1/2)"),
        v(24.0, "gamma(5)"),
        v(1.0 / 12.0, "beta(2, 3)"),
        v(0.0, "erf(0)"),
        v(1.0, "exp(0)"),
        v(PI * PI / 6.0, "zeta(2)"),
        v(PI.powi(4) / 90.0, "zeta(4)"),
        v(0.0, "ln(1)"),
        v(1.0, "ln(E)"),
        v(PI / 2.0, "asin(1)"),
        v(PI / 6.0, "asin(1/2)"),
        v(PI / 4.0, "atan(1)"),
        v(3.0f64.sqrt() / 2.0, "sin(pi/3)"),
        v(0.5f64.sqrt(), "cos(pi/4)"),
        v(120.0, "factorial(5)"),
        v(252.0, "binomial(10, 5)"),
        // ---- antiderivatives not in the main battery
        int("integral(x^2*sin(x), x)"),
        int("integral(exp(x)*cos(x), x)"),
        int("integral(1/(x^2 + x + 1), x)"),
        int("integral(1/(x*ln(x)), x)"),
        int("integral(x/(x + 1)^3, x)"),
        int("integral(sin(x)^2*cos(x)^3, x)"),
        int("integral(1/(4 - x^2), x)"),
        int("integral(x^3/(x^2 + 1), x)"),
        int("integral(atan(x), x)"),
        int("integral(asin(x), x)"),
        int("integral(ln(x)^2, x)"),
        int("integral(sqrt(1 + x^2), x)"),
        int("integral(1/sqrt(x^2 + 4), x)"),
        int("integral(exp(x)*x^3, x)"),
        int("integral(cosh(x)*sin(x), x)"),
        int("integral(1/(1 + exp(x)), x)"),
        int("integral(x^2/(1 + x^3), x)"),
        int("integral(sin(x)/(1 + cos(x)^2), x)"),
    ]
}

#[test]
fn classical_results_hold() {
    let report = std::env::var_os("CLASSICS_REPORT").is_some();
    let all = cases();
    let (mut failed_new, mut fixed_gap, mut gaps) = (Vec::new(), Vec::new(), 0);
    for case in &all {
        let (verdict, shown) = judge(case);
        if report {
            eprintln!(
                "{:<5} {:<60} -> {}  [{}]",
                if verdict.is_ok() {
                    "ok"
                } else {
                    "FAIL"
                },
                case.src,
                shown,
                verdict.as_ref().err().map_or("", String::as_str)
            );
        }
        match (verdict, case.known_gap) {
            | (Err(why), false) => failed_new.push(format!("{} -> {shown}: {why}", case.src)),
            | (Ok(()), true) => fixed_gap.push(case.src),
            | (Err(_), true) => gaps += 1,
            | (Ok(()), false) => {},
        }
    }
    eprintln!("classics: {} cases, {} known gaps", all.len(), gaps);
    assert!(
        failed_new.is_empty(),
        "wrong or missing answers:\n{}",
        failed_new.join("\n")
    );
    assert!(
        fixed_gap.is_empty(),
        "gaps that now pass, drop `.gap()`: {fixed_gap:?}"
    );
}
