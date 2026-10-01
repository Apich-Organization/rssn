//! Typst math output for terms of a [`Graph`].

use super::markup::BigOp;
use super::markup::Deriv;
use super::markup::Dialect;
use super::markup::GREEK;
use super::markup::render;
use super::markup::split_subscript;
use crate::graph::Graph;
use crate::graph::NodeId;

/// Functions that Typst knows as math operators.
const FUNCTIONS: &[(&str, &str)] = &[
    ("sin", "sin"),
    ("cos", "cos"),
    ("tan", "tan"),
    ("cot", "cot"),
    ("sec", "sec"),
    ("csc", "csc"),
    ("sinh", "sinh"),
    ("cosh", "cosh"),
    ("tanh", "tanh"),
    ("coth", "coth"),
    ("asin", "arcsin"),
    ("acos", "arccos"),
    ("atan", "arctan"),
    ("arcsin", "arcsin"),
    ("arccos", "arccos"),
    ("arctan", "arctan"),
    ("ln", "ln"),
    ("log", "log"),
    ("exp", "exp"),
    ("arg", "arg"),
    ("det", "det"),
    ("min", "min"),
    ("max", "max"),
    ("gcd", "gcd"),
    ("gamma", "Gamma"),
    ("digamma", "psi"),
    ("re", "Re"),
    ("im", "Im"),
];

/// The Typst symbol for the Greek letter called `name` (`"alpha"` gives
/// `alpha`), or `None` when `name` is not a Greek letter name.
#[must_use]
pub fn to_greek(name: &str) -> Option<&'static str> {
    GREEK.iter().find(|(n, _)| *n == name).map(|(n, _)| *n)
}

/// Renders the term `node` as Typst math (without surrounding `$`).
#[must_use]
pub fn to_typst(
    graph: &Graph,
    node: NodeId,
) -> String {
    render(graph, &Typst, node)
}

struct Typst;

/// A script argument: bare when it is a single token, parenthesised otherwise.
fn script(s: &str) -> String {
    if !s.is_empty() && s.chars().all(char::is_alphanumeric) {
        s.to_owned()
    } else {
        format!("({s})")
    }
}

fn order_script(
    s: &str,
    order: u32,
) -> String {
    if order == 1 { s.to_owned() } else { format!("{s}^{order}") }
}

impl Dialect for Typst {
    fn symbol(
        &self,
        name: &str,
    ) -> String {
        let (base, sub) = split_subscript(name);
        let head = match to_greek(base) {
            | Some(g) => g.to_owned(),
            | None if base.chars().count() == 1 => base.to_owned(),
            | None => format!("\"{base}\""),
        };
        match sub {
            | Some(s) => format!("{head}_{}", script(&s.replace('_', ","))),
            | None => head,
        }
    }

    fn text(
        &self,
        s: &str,
    ) -> String {
        format!("\"{s}\"")
    }

    fn constant(
        &self,
        name: &str,
    ) -> Option<&'static str> {
        Some(match name {
            | "pi" => "pi",
            | "oo" => "oo",
            | "I" => "i",
            | "E" => "e",
            | "true" => "top",
            | "false" => "bot",
            | _ => return None,
        })
    }

    fn frac(
        &self,
        numer: &str,
        denom: &str,
    ) -> String {
        format!("frac({numer}, {denom})")
    }

    fn sqrt(
        &self,
        x: &str,
        index: Option<&str>,
    ) -> String {
        match index {
            | Some(n) => format!("root({n}, {x})"),
            | None => format!("sqrt({x})"),
        }
    }

    fn pow(
        &self,
        base: &str,
        exp: &str,
    ) -> String {
        format!("{base}^{}", script(exp))
    }

    fn paren(
        &self,
        s: &str,
    ) -> String {
        format!("({s})")
    }

    fn mul_sep(
        &self,
        _prev: &str,
        next: &str,
    ) -> &'static str {
        if next.chars().next().is_some_and(|c| c.is_ascii_digit()) { " dot " } else { " " }
    }

    fn function(
        &self,
        name: &str,
        args: &[String],
    ) -> String {
        if let ("log", [base, x]) = (name, args) {
            return format!("log_{}({x})", script(base));
        }
        let head = FUNCTIONS
            .iter()
            .find(|(n, _)| *n == name)
            .map_or_else(|| format!("op(\"{name}\")"), |(_, c)| (*c).to_owned());
        format!("{head}({})", args.join(", "))
    }

    fn abs(
        &self,
        x: &str,
    ) -> String {
        format!("abs({x})")
    }

    fn factorial(
        &self,
        x: &str,
    ) -> String {
        format!("{x}!")
    }

    fn binom(
        &self,
        n: &str,
        k: &str,
    ) -> String {
        format!("binom({n}, {k})")
    }

    fn conj(
        &self,
        x: &str,
    ) -> String {
        format!("overline({x})")
    }

    fn bessel(
        &self,
        kind: char,
        order: &str,
        x: &str,
    ) -> String {
        format!("{kind}_{}({x})", script(order))
    }

    fn deriv(
        &self,
        d: &Deriv,
    ) -> String {
        let (op, total) = if d.partial {
            ("partial", order_script("partial", d.total))
        } else {
            ("dif", order_script("dif", d.total))
        };
        let denom = d
            .vars
            .iter()
            .map(|(v, o)| format!("{op} {}", order_script(v, *o)))
            .collect::<Vec<_>>()
            .join(" ");
        let lead = if d.partial { total } else { order_script("d", d.total) };
        if d.inline {
            format!("frac({lead} {}, {denom})", d.body)
        } else {
            format!("frac({lead}, {denom}) {}", d.body)
        }
    }

    fn integral(
        &self,
        f: &str,
        var: &str,
        bounds: Option<(&str, &str)>,
    ) -> String {
        match bounds {
            | Some((a, b)) => format!("integral_{}^{} {f} dif {var}", script(a), script(b)),
            | None => format!("integral {f} dif {var}"),
        }
    }

    fn big_op(
        &self,
        kind: BigOp,
        f: &str,
        var: &str,
        lo: &str,
        hi: &str,
    ) -> String {
        let op = if kind == BigOp::Sum { "sum" } else { "product" };
        format!("{op}_({var}={lo})^{} {f}", script(hi))
    }

    fn limit(
        &self,
        f: &str,
        var: &str,
        to: &str,
    ) -> String {
        format!("lim_({var} -> {to}) {f}")
    }

    fn matrix(
        &self,
        rows: &[Vec<String>],
    ) -> String {
        let body = rows
            .iter()
            .map(|r| r.join(", "))
            .collect::<Vec<_>>()
            .join("; ");
        format!("mat({body})")
    }

    fn list(
        &self,
        items: &[String],
    ) -> String {
        format!("[{}]", items.join(", "))
    }

    fn infix(
        &self,
        name: &str,
    ) -> Option<&'static str> {
        Some(match name {
            | "and" => " and ",
            | "or" => " or ",
            | "xor" => " xor ",
            | "implies" => " => ",
            | "iff" => " <=> ",
            | "lt" => " < ",
            | "gt" => " > ",
            | "le" => " <= ",
            | "ge" => " >= ",
            | "ne" | "neq" => " != ",
            | _ => return None,
        })
    }

    fn not(
        &self,
        x: &str,
    ) -> String {
        format!("not {x}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render_src(src: &str) -> String {
        let mut g = Graph::new();
        let _e = crate::graph::Engine::install(&mut g, &crate::rules::standard()).unwrap();
        let n = g.parse(src).unwrap();
        to_typst(&g, n)
    }

    #[test]
    fn constructs() {
        let cases: &[(&str, &str)] = &[
            (r#"list(1, 2, 3)"#, r#"[1, 2, 3]"#),
            (r#"list(list(1, 2), list(3, 4))"#, r#"mat(1, 2; 3, 4)"#),
            (r#"(a+b)/c"#, r#"frac(a + b, c)"#),
            (r#"diff(f(x), x)"#, r#"frac(d, dif x) f"#),
            (r#"diff(diff(u(x, t), x), x)"#, r#"frac(partial^2 u, partial x^2)"#),
            (r#"x^2/(2*y) + sqrt(x+1) - 3/4*z"#, r#"frac(x^2, 2 y) - frac(3 z, 4) + sqrt(x + 1)"#),
            (r#"x^(1/n)"#, r#"root(n, x)"#),
            (r#"exp(x^2)"#, r#"e^(x^2)"#),
            (r#"alpha*beta"#, r#"alpha beta"#),
            (r#"sum(n^2, n, 1, 10)"#, r#"sum_(n=1)^10 n^2"#),
            (r#"limit(sin(x)/x, x, 0)"#, r#"lim_(x -> 0) frac(sin(x), x)"#),
            (r#"defint(f(x), x, 0, 1)"#, r#"integral_0^1 f(x) dif x"#),
            (r#"abs(x)"#, r#"abs(x)"#),
            (r#"binomial(n,k)"#, r#"binom(n, k)"#),
            (r#"-x/2"#, r#"-frac(x, 2)"#),
            (r#"2*pi*r"#, r#"2 r pi"#),
            (r#"foo(x, y)"#, r#""foo"(x, y)"#),
            (r#"x_1^2"#, r#"x_1^2"#),
        ];
        for (src, want) in cases {
            assert_eq!(render_src(src), *want, "source: {src}");
        }
    }
}
