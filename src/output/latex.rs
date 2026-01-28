use std::collections::HashMap;

use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

/// Converts an expression to a LaTeX string.

#[must_use]

pub fn to_latex(expr: &Expr) -> String {

    to_latex_prec(expr, 0)
}

/// This is a placeholder struct to hold the result of a sub-expression,
/// including its precedence.
#[derive(Clone)]

struct LatexResult {
    precedence: u8,
    content: String,
}

/// Converts an expression to a LaTeX string with precedence handling.
/// This function is iterative to avoid stack overflows with deep expression trees.

pub(crate) fn to_latex_prec(
    root_expr: &Expr,
    root_precedence: u8,
) -> String {

    let mut results: HashMap<
        Expr,
        LatexResult,
    > = HashMap::new();

    let mut stack: Vec<Expr> =
        vec![root_expr.clone()];

    let mut visited =
        std::collections::HashSet::new(
        );

    while let Some(expr) = stack.last()
    {

        if results.contains_key(expr) {

            stack.pop();

            continue;
        }

        let children = expr.children();

        if children.is_empty()
            || visited.contains(expr)
        {

            let current_expr = stack
                .pop()
                .expect("Expr present");

            let children =
                current_expr.children();

            let get_child_res = |i: usize| -> &LatexResult {
                &results[&children[i]]
            };

            let get_child_str_with_parens = |i: usize, prec: u8| -> String {
                let child_res = get_child_res(i);
                if child_res.precedence < prec {
                    format!(r"\left( {} \right)", child_res.content)
                } else {
                    child_res.content.clone()
                }
            };

            let (op_prec, s) = match current_expr.op() {
                DagOp::Constant(c) => (10, c.to_string()),
                DagOp::BigInt(i) => (10, i.to_string()),
                DagOp::Rational(r) => (10, format!(r"\frac{{{}}}{{{}}}", r.numer(), r.denom())),
                DagOp::Variable(s) => (10, to_greek(&s)),
                DagOp::Add => (1, format!("{} + {}", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Sub => (1, format!("{} - {}", get_child_res(0).content, get_child_str_with_parens(1, 2))),
                DagOp::Mul => (2, format!("{} {}", get_child_str_with_parens(0, 2), get_child_str_with_parens(1, 2))),
                DagOp::Div => (2, format!(r"\frac{{{}}}{{{}}}", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Power => (3, format!("{}^{{{}}}", get_child_str_with_parens(0, 4), get_child_res(1).content)),
                DagOp::Neg => (2, format!("-{}", get_child_str_with_parens(0, 2))),
                DagOp::Sqrt => (10, format!(r"\sqrt{{{}}}", get_child_res(0).content)),
                DagOp::Abs => (10, format!(r"\left| {} \right|", get_child_res(0).content)),
                DagOp::Pi => (10, r"\pi".to_string()),
                DagOp::E => (10, "e".to_string()),
                DagOp::Log => (10, format!(r"\ln\left({}\right)", get_child_res(0).content)),
                DagOp::LogBase => (10, format!(r"\log_{{{}}}\left({}\right)", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Exp => (10, format!(r"e^{{{}}}", get_child_res(0).content)),
                DagOp::Sin => (10, format!(r"\sin\left({}\right)", get_child_res(0).content)),
                DagOp::Cos => (10, format!(r"\cos\left({}\right)", get_child_res(0).content)),
                DagOp::Tan => (10, format!(r"\tan\left({}\right)", get_child_res(0).content)),
                DagOp::Csc => (10, format!(r"\csc\left({}\right)", get_child_res(0).content)),
                DagOp::Sec => (10, format!(r"\sec\left({}\right)", get_child_res(0).content)),
                DagOp::Cot => (10, format!(r"\cot\left({}\right)", get_child_res(0).content)),
                DagOp::ArcSin => (10, format!(r"\arcsin\left({}\right)", get_child_res(0).content)),
                DagOp::ArcCos => (10, format!(r"\arccos\left({}\right)", get_child_res(0).content)),
                DagOp::ArcTan => (10, format!(r"\arctan\left({}\right)", get_child_res(0).content)),
                DagOp::Factorial => (10, format!("{}!", get_child_str_with_parens(0, 10))),
                DagOp::Derivative(var) => (10, format!(r"\frac{{d}}{{d{}}} \left( {} \right)", to_greek(&var), get_child_res(0).content)),
                DagOp::Integral => (10, format!(r"\int_{{{}}}^{{{}}} {} \,d{}", get_child_res(2).content, get_child_res(3).content, get_child_res(0).content, get_child_res(1).content)),
                DagOp::Sum => (10, format!(r"\sum_{{{}={}}}^{{{}}} {}", get_child_res(1).content, get_child_res(2).content, get_child_res(3).content, get_child_res(0).content)),
                DagOp::Summation(s) => (10, format!(r"\sum_{{{}={}}}^{{{}}} {}", to_greek(&s), get_child_res(1).content, get_child_res(2).content, get_child_res(0).content)),
                DagOp::Product(s) => (10, format!(r"\prod_{{{}={}}}^{{{}}} {}", to_greek(&s), get_child_res(1).content, get_child_res(2).content, get_child_res(0).content)),
                DagOp::Limit(var) => (10, format!(r"\lim_{{{} \to {}}} {}", to_greek(&var), get_child_res(1).content, get_child_res(0).content)),
                DagOp::Eq => (0, format!("{} = {}", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Binomial => (10, format!(r"\binom{{{}}}{{{}}}", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Matrix { rows: _, cols } => {
                    let body = children.chunks(cols).map(|row| {
                        row.iter().map(|elem| results[elem].content.clone()).collect::<Vec<_>>().join(" & ")
                    }).collect::<Vec<_>>().join(r" \\ ");
                    (10, format!(r"\begin{{pmatrix}} {body} \end{{pmatrix}}"))
                },
                DagOp::Vector => {
                    let body = children.iter().map(|elem| results[elem].content.clone()).collect::<Vec<_>>().join(r" \\ ");
                    (10, format!(r"\begin{{pmatrix}} {body} \end{{pmatrix}}"))
                },
                _ => (10, current_expr.to_string()),
            };

            results.insert(
                current_expr,
                LatexResult {
                    precedence: op_prec,
                    content: s,
                },
            );
        } else {

            visited
                .insert(expr.clone());

            for child in children
                .iter()
                .rev()
            {

                stack.push(
                    child.clone(),
                );
            }
        }
    }

    let final_res = results
        .get(root_expr)
        .expect("Result missing");

    if final_res.precedence
        < root_precedence
    {

        format!(
            r"\left( {} \right)",
            final_res.content
        )
    } else {

        final_res
            .content
            .clone()
    }
}

/// Helper to add parentheses if needed. This function is now simplified as the main
/// iterative function handles most of the logic.

#[must_use]

pub fn to_latex_prec_with_parens(
    expr: &Expr,
    precedence: u8,
) -> String {

    let op = expr.op();

    let op_prec = match op {
        | DagOp::Add | DagOp::Sub => 1,
        | _ => 10,
    };

    let s =
        to_latex_prec(expr, precedence);

    if op_prec < precedence {

        format!(r"\left( {s} \right)")
    } else {

        s
    }
}

/// Converts common Greek letter names to LaTeX.

#[must_use]

pub fn to_greek(s: &str) -> String {

    match s {
        | "alpha" => r"\alpha".into(),
        | "beta" => r"\beta".into(),
        | "gamma" => r"\gamma".into(),
        | "delta" => r"\delta".into(),
        | "epsilon" => {
            r"\epsilon".into()
        },
        | "zeta" => r"\zeta".into(),
        | "eta" => r"\eta".into(),
        | "theta" => r"\theta".into(),
        | "iota" => r"\iota".into(),
        | "kappa" => r"\kappa".into(),
        | "lambda" => r"\lambda".into(),
        | "mu" => r"\mu".into(),
        | "nu" => r"\nu".into(),
        | "xi" => r"\xi".into(),
        | "pi" => r"\pi".into(),
        | "rho" => r"\rho".into(),
        | "sigma" => r"\sigma".into(),
        | "tau" => r"\tau".into(),
        | "upsilon" => {
            r"\upsilon".into()
        },
        | "phi" => r"\phi".into(),
        | "chi" => r"\chi".into(),
        | "psi" => r"\psi".into(),
        | "Alpha" => r"A".into(),
        | "Beta" => r"B".into(),
        | "Gamma" => r"\Gamma".into(),
        | "Delta" => r"\Delta".into(),
        | "Epsilon" => r"E".into(),
        | "Zeta" => r"Z".into(),
        | "Eta" => r"H".into(),
        | "Theta" => r"\Theta".into(),
        | "Iota" => r"I".into(),
        | "Kappa" => r"K".into(),
        | "Lambda" => r"\Lambda".into(),
        | "Mu" => r"M".into(),
        | "Nu" => r"N".into(),
        | "Xi" => r"\Xi".into(),
        | "Pi" => r"\Pi".into(),
        | "Rho" => r"P".into(),
        | "Sigma" => r"\Sigma".into(),
        | "Tau" => r"T".into(),
        | "Upsilon" => {
            r"\Upsilon".into()
        },
        | "Phi" => r"\Phi".into(),
        | "Chi" => r"X".into(),
        | "Psi" => r"\Psi".into(),
        | "omega" => r"\omega".into(),
        | "Omega" => r"\Omega".into(),
        | _ => s.to_string(),
    }
}

#[cfg(test)]

mod tests {

    use super::*;
    use crate::prelude::Expr;

    #[test]

    fn test_to_latex_basic() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let y = Expr::Variable(
            "y".to_string(),
        );

        // (x + y) * x normalized to x * (x + y)
        let expr = Expr::new_mul(
            Expr::new_add(
                x.clone(),
                y.clone(),
            ),
            x.clone(),
        );

        let latex = to_latex(&expr);

        assert_eq!(
            latex,
            r"x \left( x + y \right)"
        );
    }

    #[test]

    fn test_to_latex_fractions() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::new_div(
            x,
            Expr::new_constant(2.0),
        );

        let latex = to_latex(&expr);

        assert_eq!(
            latex,
            r"\frac{x}{2}"
        );
    }

    #[test]

    fn test_to_latex_greek() {

        let alpha = Expr::Variable(
            "alpha".to_string(),
        );

        let omega = Expr::Variable(
            "omega".to_string(),
        );

        let expr =
            Expr::new_add(alpha, omega);

        let latex = to_latex(&expr);

        assert!(
            latex.contains(r"\alpha")
                && latex.contains(
                    r"\omega"
                )
        );
    }

    #[test]

    fn test_to_latex_matrix() {

        let expr =
            Expr::new_matrix(vec![
                vec![
                    Expr::new_constant(1.0),
                    Expr::Variable(
                        "a".into(),
                    ),
                ],
                vec![
                    Expr::Variable(
                        "b".into(),
                    ),
                    Expr::new_constant(2.0),
                ],
            ]);

        let latex = to_latex(&expr);

        assert!(latex.contains(
            r"\begin{pmatrix}"
        ));

        assert!(
            latex.contains("1 & a")
        );

        assert!(
            latex.contains("b & 2")
        );
    }

    #[test]

    fn test_to_latex_integral() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::Integral {
            integrand:
                std::sync::Arc::new(
                    Expr::new_pow(
                        x.clone(),
                        Expr::new_constant(
                            2.0,
                        ),
                    ),
                ),
            var: std::sync::Arc::new(
                x.clone(),
            ),
            lower_bound:
                std::sync::Arc::new(
                    Expr::new_constant(0.0),
                ),
            upper_bound:
                std::sync::Arc::new(
                    Expr::new_constant(1.0),
                ),
        };

        let latex = to_latex(&expr);

        println!("{}", latex);

        assert!(
            latex.contains(
                r"\int_{0}^{1}"
            )
        );

        assert!(
            latex.contains(r"x^{2}")
        );

        assert!(
            latex.contains(r"\,dx")
        );
    }

    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]
        #[test]
        fn prop_no_panic_to_latex(
            depth in 1..4usize,
        ) {
            let x = Expr::Variable("x".to_string());
            let mut expr = x.clone();
            for _ in 0..depth {
                expr = Expr::new_add(expr.clone(), Expr::new_constant(1.0));
                expr = Expr::new_div(expr.clone(), Expr::new_sqrt(x.clone()));
                expr = Expr::new_pow(expr.clone(), Expr::new_constant(2.0));
            }
            let latex = to_latex(&expr);
            assert!(!latex.is_empty());
        }
    }
}
