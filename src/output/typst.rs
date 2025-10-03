use crate::symbolic::core::Expr;

/// Converts an expression to a Typst string.
pub fn to_typst(expr: &Expr) -> String {
    format!("${}$", to_typst_prec(expr, 0))
}

/// Recursive helper with precedence handling for parentheses.
pub(crate) fn to_typst_prec(expr: &Expr, precedence: u8) -> String {
    let (op_prec, s) = match expr {
        Expr::Add(a, b) => (
            1,
            format!("{} + {}", to_typst_prec(a, 1), to_typst_prec(b, 1)),
        ),
        Expr::Sub(a, b) => (
            1,
            format!("{} - {}", to_typst_prec(a, 1), to_typst_prec(b, 2)),
        ),
        Expr::Mul(a, b) => (
            2,
            format!("{} * {}", to_typst_prec(a, 2), to_typst_prec(b, 2)),
        ),
        Expr::Div(a, b) => (
            2,
            format!("frac({}, {})", to_typst_prec(a, 0), to_typst_prec(b, 0)),
        ),
        Expr::Power(b, e) => (
            3,
            format!("{}^({})", to_typst_prec(b, 3), to_typst_prec(e, 0)),
        ),
        Expr::Neg(a) => (2, format!("-{}", to_typst_prec(a, 2))),
        Expr::Sqrt(a) => (4, format!("sqrt({})", to_typst_prec(a, 0))),
        Expr::Sin(a) => (4, format!("sin({})", to_typst_prec(a, 0))),
        Expr::Cos(a) => (4, format!("cos({})", to_typst_prec(a, 0))),
        Expr::Tan(a) => (4, format!("tan({})", to_typst_prec(a, 0))),
        Expr::Log(a) => (4, format!("ln({})", to_typst_prec(a, 0))),
        Expr::Exp(a) => (3, format!("e^({})", to_typst_prec(a, 0))),
        Expr::Integral {
            integrand,
            var,
            lower_bound,
            upper_bound,
        } => (
            5,
            format!(
                "integral_({})^( {}) {} dif {}",
                to_typst_prec(lower_bound, 0),
                to_typst_prec(upper_bound, 0),
                to_typst_prec(integrand, 0),
                to_typst_prec(var, 0)
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
                "sum_({}={})^( {}) {}",
                to_typst_prec(var, 0),
                to_typst_prec(from, 0),
                to_typst_prec(to, 0),
                to_typst_prec(body, 0)
            ),
        ),
        Expr::Matrix(rows) => {
            let body = rows
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|elem| to_typst_prec(elem, 0))
                        .collect::<Vec<_>>()
                        .join(", ")
                })
                .collect::<Vec<_>>()
                .join("; ");
            (10, format!("mat({})", body))
        }
        Expr::Pi => (10, "pi".to_string()),
        Expr::E => (10, "e".to_string()),
        Expr::Constant(c) => (10, c.to_string()),
        Expr::BigInt(i) => (10, i.to_string()),
        Expr::Variable(s) => (10, s.clone()),
        _ => (10, expr.to_string()), // Fallback
    };

    if op_prec < precedence {
        format!("({})", s)
    } else {
        s
    }
}
