use crate::symbolic::core::Expr;

/// Converts an expression to a LaTeX string.
pub fn to_latex(expr: &Expr) -> String {
    to_latex_prec(expr, 0)
}

/// Recursive helper with precedence handling for parentheses.
pub(crate) fn to_latex_prec(expr: &Expr, precedence: u8) -> String {
    let (op_prec, s) = match expr {
        // Numbers and Variables
        Expr::Constant(c) => (10, c.to_string()),
        Expr::BigInt(i) => (10, i.to_string()),
        Expr::Rational(r) => (10, format!(r"\frac{{{}}}{{{}}}", r.numer(), r.denom())),
        Expr::Variable(s) => (10, to_greek(s)),

        // Basic Operations
        Expr::Add(a, b) => (
            1,
            format!("{} + {}", to_latex_prec(a, 1), to_latex_prec(b, 1)),
        ),
        Expr::Sub(a, b) => (
            1,
            format!("{} - {}", to_latex_prec(a, 1), to_latex_prec(b, 2)),
        ),
        Expr::Mul(a, b) => (
            2,
            format!(
                "{}{}",
                to_latex_prec_with_parens(a, 2),
                to_latex_prec_with_parens(b, 2)
            ),
        ),
        Expr::Div(a, b) => (
            2,
            format!(
                r"\frac{{{}}}{{{}}}",
                to_latex_prec(a, 0),
                to_latex_prec(b, 0)
            ),
        ),
        Expr::Power(b, e) => (
            3,
            format!("{{{}}}^{{{}}}", to_latex_prec(b, 3), to_latex_prec(e, 0)),
        ),
        Expr::Neg(a) => (2, format!("-{}", to_latex_prec_with_parens(a, 2))),

        // Common Functions & Constants
        Expr::Sqrt(a) => (4, format!("\\sqrt{{{}}}", to_latex_prec(a, 0))),
        Expr::Abs(a) => (10, format!("\\left| {} \\right|", to_latex_prec(a, 0))),
        Expr::Pi => (10, "\\pi".to_string()),
        Expr::E => (10, "e".to_string()),

        // Logarithmic and Exponential
        Expr::Log(a) => (4, format!("\\ln\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::LogBase(b, a) => (
            4,
            format!(
                "\\log_{{{}}}\\left({}\\right)",
                to_latex_prec(b, 0),
                to_latex_prec(a, 0)
            ),
        ),
        Expr::Exp(a) => (3, format!("e^{{{}}}", to_latex_prec(a, 0))),

        // Trigonometric
        Expr::Sin(a) => (4, format!("\\sin\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::Cos(a) => (4, format!("\\cos\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::Tan(a) => (4, format!("\\tan\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::Csc(a) => (4, format!("\\csc\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::Sec(a) => (4, format!("\\sec\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::Cot(a) => (4, format!("\\cot\\left({}\\right)", to_latex_prec(a, 0))),

        // Inverse Trigonometric
        Expr::ArcSin(a) => (4, format!("\\arcsin\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::ArcCos(a) => (4, format!("\\arccos\\left({}\\right)", to_latex_prec(a, 0))),
        Expr::ArcTan(a) => (4, format!("\\arctan\\left({}\\right)", to_latex_prec(a, 0))),

        // Calculus
        Expr::Derivative(body, var) => (
            5,
            format!(
                "\\frac{{d}}{{d{}}} \\left[ {} \\right]",
                var,
                to_latex_prec(body, 0)
            ),
        ),
        Expr::Integral {
            integrand,
            var,
            lower_bound,
            upper_bound,
        } => (
            5,
            format!(
                "\\int_{{{}}}^{{{}}} {} \\,d{}",
                to_latex_prec(lower_bound, 0),
                to_latex_prec(upper_bound, 0),
                to_latex_prec(integrand, 0),
                to_latex_prec(var, 0)
            ),
        ),
        Expr::Sum {
            body,
            var,
            from,
            to,
        } => (
            5,
            format!(
                "\\sum_{{{}={}}}^{{{}}} {}",
                to_latex_prec(var, 0),
                to_latex_prec(from, 0),
                to_latex_prec(to, 0),
                to_latex_prec(body, 0)
            ),
        ),
        Expr::Limit(body, var, to) => (
            5,
            format!(
                "\\lim_{{{} \\to {}}} {}",
                to_latex_prec(&Expr::Variable(var.clone()), 0),
                to_latex_prec(to, 0),
                to_latex_prec(body, 0)
            ),
        ),

        // Other
        Expr::Eq(a, b) => (
            0,
            format!("{} = {}", to_latex_prec(a, 0), to_latex_prec(b, 0)),
        ),
        Expr::Binomial(n, k) => (
            5,
            format!(
                "\\binom{{{}}}{{{}}}",
                to_latex_prec(n, 0),
                to_latex_prec(k, 0)
            ),
        ),
        Expr::Matrix(rows) => {
            let body = rows
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|elem| to_latex_prec(elem, 0))
                        .collect::<Vec<_>>()
                        .join(" & ")
                })
                .collect::<Vec<_>>()
                .join(" \\\\ ");
            (10, format!("\\begin{{pmatrix}}{}\\end{{pmatrix}}", body))
        }

        _ => (10, expr.to_string()), // Fallback
    };

    if op_prec < precedence {
        format!("\\left( {} \\right)", s)
    } else {
        s
    }
}

/// Helper to add parentheses if needed.
pub(crate) fn to_latex_prec_with_parens(expr: &Expr, precedence: u8) -> String {
    let (op_prec, s) = match expr {
        Expr::Add(_, _) | Expr::Sub(_, _) => (1, to_latex_prec(expr, precedence)),
        _ => (10, to_latex_prec(expr, precedence)),
    };
    if op_prec < precedence {
        format!(r"\left( {} \right)", s)
    } else {
        s
    }
}

/// Converts common Greek letter names to LaTeX.
pub(crate) fn to_greek(s: &str) -> String {
    match s {
        "alpha" => "\\alpha".to_string(),
        "beta" => "\\beta".to_string(),
        "gamma" => "\\gamma".to_string(),
        "delta" => "\\delta".to_string(),
        "epsilon" => "\\epsilon".to_string(),
        "zeta" => "\\zeta".to_string(),
        "eta" => "\\eta".to_string(),
        "theta" => "\\theta".to_string(),
        "iota" => "\\iota".to_string(),
        "kappa" => "\\kappa".to_string(),
        "lambda" => "\\lambda".to_string(),
        "mu" => "\\mu".to_string(),
        "nu" => "\\nu".to_string(),
        "xi" => "\\xi".to_string(),
        "pi" => "\\pi".to_string(),
        "rho" => "\\rho".to_string(),
        "sigma" => "\\sigma".to_string(),
        "tau" => "\\tau".to_string(),
        "upsilon" => "\\upsilon".to_string(),
        "phi" => "\\phi".to_string(),
        "chi" => "\\chi".to_string(),
        "psi" => "\\psi".to_string(),
        "omega" => "\\omega".to_string(),
        _ => s.to_string(),
    }
}
