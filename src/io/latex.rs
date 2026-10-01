//! LaTeX output for terms of a [`Graph`].

use super::markup::BigOp;
use super::markup::Deriv;
use super::markup::Dialect;
use super::markup::GREEK;
use super::markup::render;
use super::markup::split_subscript;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::print::PREC_ADD;

/// Known functions and their LaTeX commands.
const FUNCTIONS: &[(&str, &str)] = &[
    ("sin", r"\sin"),
    ("cos", r"\cos"),
    ("tan", r"\tan"),
    ("cot", r"\cot"),
    ("sec", r"\sec"),
    ("csc", r"\csc"),
    ("sinh", r"\sinh"),
    ("cosh", r"\cosh"),
    ("tanh", r"\tanh"),
    ("coth", r"\coth"),
    ("asin", r"\arcsin"),
    ("acos", r"\arccos"),
    ("atan", r"\arctan"),
    ("arcsin", r"\arcsin"),
    ("arccos", r"\arccos"),
    ("arctan", r"\arctan"),
    ("ln", r"\ln"),
    ("log", r"\log"),
    ("exp", r"\exp"),
    ("arg", r"\arg"),
    ("det", r"\det"),
    ("min", r"\min"),
    ("max", r"\max"),
    ("gcd", r"\gcd"),
    ("gamma", r"\Gamma"),
    ("digamma", r"\psi"),
    ("re", r"\operatorname{Re}"),
    ("im", r"\operatorname{Im}"),
];

/// The LaTeX command for the Greek letter called `name` (`"alpha"` gives
/// `\alpha`), or `None` when `name` is not a Greek letter name.
#[must_use]
pub fn to_greek(name: &str) -> Option<&'static str> {
    GREEK.iter().find(|(n, _)| *n == name).map(|(_, l)| *l)
}

/// Renders the term `node` as LaTeX math (without surrounding `$`).
#[must_use]
pub fn to_latex(
    graph: &Graph,
    node: NodeId,
) -> String {
    render(graph, &Latex, node)
}

/// Renders `node` as LaTeX, parenthesised when `precedence` demands it.
///
/// The result is wrapped in `\left( .. \right)` when the term is
/// a sum, a difference, a negative number or an equation and `precedence`
/// asks for more than additive binding (`precedence > 1`, the legacy
/// convention where sums and differences bind at 1 and everything else
/// tighter).
#[must_use]
pub fn to_latex_prec_with_parens(
    graph: &Graph,
    node: NodeId,
    precedence: u8,
) -> String {
    let text = to_latex(graph, node);
    if precedence > PREC_ADD && graph.precedence(node) <= PREC_ADD {
        Latex.paren(&text)
    } else {
        text
    }
}

struct Latex;

fn braced(s: &str) -> String {
    format!("{{{s}}}")
}

fn order_script(
    s: &str,
    order: u32,
) -> String {
    if order == 1 { s.to_owned() } else { format!("{s}^{}", digit_group(order)) }
}

fn digit_group(n: u32) -> String {
    if n < 10 { n.to_string() } else { braced(&n.to_string()) }
}

impl Dialect for Latex {
    fn symbol(
        &self,
        name: &str,
    ) -> String {
        let (base, sub) = split_subscript(name);
        let head = match to_greek(base) {
            | Some(g) => g.to_owned(),
            | None if base.chars().count() == 1 => base.to_owned(),
            | None => format!(r"\mathrm{{{base}}}"),
        };
        match sub {
            | Some(s) => format!("{head}_{{{}}}", s.replace('_', ",")),
            | None => head,
        }
    }

    fn text(
        &self,
        s: &str,
    ) -> String {
        format!(r"\text{{{s}}}")
    }

    fn constant(
        &self,
        name: &str,
    ) -> Option<&'static str> {
        Some(match name {
            | "pi" => r"\pi",
            | "oo" => r"\infty",
            | "I" => "i",
            | "E" => "e",
            | "true" => r"\top",
            | "false" => r"\bot",
            | _ => return None,
        })
    }

    fn frac(
        &self,
        numer: &str,
        denom: &str,
    ) -> String {
        format!(r"\frac{{{numer}}}{{{denom}}}")
    }

    fn sqrt(
        &self,
        x: &str,
        index: Option<&str>,
    ) -> String {
        match index {
            | Some(n) => format!(r"\sqrt[{n}]{{{x}}}"),
            | None => format!(r"\sqrt{{{x}}}"),
        }
    }

    fn pow(
        &self,
        base: &str,
        exp: &str,
    ) -> String {
        format!("{base}^{}", braced(exp))
    }

    fn paren(
        &self,
        s: &str,
    ) -> String {
        format!(r"\left( {s} \right)")
    }

    fn mul_sep(
        &self,
        prev: &str,
        next: &str,
    ) -> &'static str {
        let numeric_prev = prev.chars().all(|c| c.is_ascii_digit() || c == '.');
        let starts = next.chars().next();
        if starts.is_some_and(|c| c.is_ascii_digit()) {
            r" \cdot "
        } else if numeric_prev {
            ""
        } else {
            " "
        }
    }

    fn function(
        &self,
        name: &str,
        args: &[String],
    ) -> String {
        if let ("log", [base, x]) = (name, args) {
            return format!(r"\log_{{{base}}}({x})");
        }
        let head = FUNCTIONS
            .iter()
            .find(|(n, _)| *n == name)
            .map_or_else(|| format!(r"\operatorname{{{}}}", name.replace('_', r"\_")), |(_, c)| (*c).to_owned());
        format!("{head}({})", args.join(", "))
    }

    fn abs(
        &self,
        x: &str,
    ) -> String {
        format!(r"\left| {x} \right|")
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
        format!(r"\binom{{{n}}}{{{k}}}")
    }

    fn conj(
        &self,
        x: &str,
    ) -> String {
        format!(r"\overline{{{x}}}")
    }

    fn bessel(
        &self,
        kind: char,
        order: &str,
        x: &str,
    ) -> String {
        format!("{kind}_{{{order}}}({x})")
    }

    fn deriv(
        &self,
        d: &Deriv,
    ) -> String {
        let (op, total) = if d.partial {
            (r"\partial", order_script(r"\partial", d.total))
        } else {
            ("d", order_script("d", d.total))
        };
        let denom = d
            .vars
            .iter()
            .map(|(v, o)| format!("{op}{}{}", if d.partial { " " } else { "" }, order_script(v, *o)))
            .collect::<Vec<_>>()
            .join(" ");
        if d.inline {
            format!(r"\frac{{{total} {}}}{{{denom}}}", d.body)
        } else {
            format!(r"\frac{{{total}}}{{{denom}}} {}", d.body)
        }
    }

    fn integral(
        &self,
        f: &str,
        var: &str,
        bounds: Option<(&str, &str)>,
    ) -> String {
        match bounds {
            | Some((a, b)) => format!(r"\int_{{{a}}}^{{{b}}} {f} \, d{var}"),
            | None => format!(r"\int {f} \, d{var}"),
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
        let op = if kind == BigOp::Sum { r"\sum" } else { r"\prod" };
        format!("{op}_{{{var}={lo}}}^{{{hi}}} {f}")
    }

    fn limit(
        &self,
        f: &str,
        var: &str,
        to: &str,
    ) -> String {
        format!(r"\lim_{{{var} \to {to}}} {f}")
    }

    fn matrix(
        &self,
        rows: &[Vec<String>],
    ) -> String {
        let body = rows
            .iter()
            .map(|r| r.join(" & "))
            .collect::<Vec<_>>()
            .join(r" \\ ");
        format!(r"\begin{{pmatrix}} {body} \end{{pmatrix}}")
    }

    fn list(
        &self,
        items: &[String],
    ) -> String {
        format!(r"\left[ {} \right]", items.join(", "))
    }

    fn infix(
        &self,
        name: &str,
    ) -> Option<&'static str> {
        Some(match name {
            | "and" => r" \land ",
            | "or" => r" \lor ",
            | "xor" => r" \oplus ",
            | "implies" => r" \Rightarrow ",
            | "iff" => r" \Leftrightarrow ",
            | "lt" => " < ",
            | "gt" => " > ",
            | "le" => r" \le ",
            | "ge" => r" \ge ",
            | "ne" | "neq" => r" \ne ",
            | _ => return None,
        })
    }

    fn not(
        &self,
        x: &str,
    ) -> String {
        format!(r"\lnot {x}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render_src(src: &str) -> String {
        let mut g = Graph::new();
        let _e = crate::graph::Engine::install(&mut g, &crate::rules::standard()).unwrap();
        let n = g.parse(src).unwrap();
        to_latex(&g, n)
    }

    #[test]
    fn constructs() {
        let cases: &[(&str, &str)] = &[
            (r"list(1, 2, 3)", r"\left[ 1, 2, 3 \right]"),
            (r"list(list(1, 2), list(3, 4))", r"\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}"),
            (r"(a+b)/c", r"\frac{a + b}{c}"),
            (r"diff(f(x), x)", r"\frac{d}{dx} f"),
            (r"diff(diff(u(x, t), x), x)", r"\frac{\partial^2 u}{\partial x^2}"),
            (r"x^2/(2*y) + sqrt(x+1) - 3/4*z", r"\frac{x^{2}}{2y} - \frac{3z}{4} + \sqrt{x + 1}"),
            (r"x^(1/n)", r"\sqrt[n]{x}"),
            (r"x^(2/3)", r"x^{\frac{2}{3}}"),
            (r"x_1^2", r"x_{1}^{2}"),
            (r"exp(x^2)", r"e^{x^{2}}"),
            (r"alpha*beta", r"\alpha \beta"),
            (r"sum(n^2, n, 1, 10)", r"\sum_{n=1}^{10} n^{2}"),
            (r"limit(sin(x)/x, x, 0)", r"\lim_{x \to 0} \frac{\sin(x)}{x}"),
            (r"defint(f(x), x, 0, 1)", r"\int_{0}^{1} f(x) \, dx"),
            (r"abs(x)", r"\left| x \right|"),
            (r"binomial(n,k)", r"\binom{n}{k}"),
            (r"-x/2", r"-\frac{x}{2}"),
            (r"conj(z)", r"\overline{z}"),
            (r"2*pi*r", r"2r \pi"),
            (r"foo(x, y)", r"\mathrm{foo}(x, y)"),
            (r"factorial(n+1)", r"\left( n + 1 \right)!"),
            (r"x^(-1/2)", r"\frac{1}{\sqrt{x}}"),
        ];
        for (src, want) in cases {
            assert_eq!(render_src(src), *want, "source: {src}");
        }
    }

    #[test]
    fn greek_names() {
        assert_eq!(to_greek("alpha"), Some("\\alpha"));
        assert_eq!(to_greek("notgreek"), None);
    }
}
