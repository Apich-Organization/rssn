use crate::symbolic::core::{DagOp, Expr};
use std::collections::HashMap;

/// Converts an expression to a LaTeX string.

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

    let mut results: HashMap<*const Expr, LatexResult> = HashMap::new();

    let mut stack: Vec<Expr> = vec![root_expr.clone()];

    while let Some(expr) = stack.pop() {

        let expr_ptr = &expr as *const Expr;

        if results.contains_key(&expr_ptr) {

            stack.pop();

            continue;
        }

        let children = expr.children();

        let all_children_processed = children
            .iter()
            .all(|c| results.contains_key(&(c as *const Expr)));

        if all_children_processed {

            let current_expr = stack
                .pop()
                .expect("Value is valid");

            let current_expr_ptr = &current_expr as *const Expr;

            let get_child_res = |i: usize| -> &LatexResult {

                &results[&(&children[i] as *const Expr)]
            };

            let get_child_str_with_parens = |i: usize, prec: u8| -> String {

                let child_res = get_child_res(i);

                if child_res.precedence < prec {

                    format!(
                        r"\left( {} \right)",
                        child_res.content
                    )
                } else {

                    child_res
                        .content
                        .clone()
                }
            };

            let (op_prec, s) = match current_expr.op() {
                DagOp::Constant(c) => (
                    10,
                    c.into_inner()
                        .to_string(),
                ),
                DagOp::BigInt(i) => (10, i.to_string()),
                DagOp::Rational(r) => (
                    10,
                    format!(
                        r"\frac{{{}}}{{{}}}",
                        r.numer(),
                        r.denom()
                    ),
                ),
                DagOp::Variable(s) => (10, to_greek(&s)),
                DagOp::Add => (
                    1,
                    format!(
                        "{} + {}",
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Sub => (
                    1,
                    format!(
                        "{} - {}",
                        get_child_res(0).content,
                        get_child_str_with_parens(1, 2)
                    ),
                ),
                DagOp::Mul => (
                    2,
                    format!(
                        "{}{}",
                        get_child_str_with_parens(0, 2),
                        get_child_str_with_parens(1, 2)
                    ),
                ),
                DagOp::Div => (
                    2,
                    format!(
                        r"\frac{{{}}}{{{}}}",
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Power => (
                    3,
                    format!(
                        "{{{}}}^{{{}}}",
                        get_child_str_with_parens(0, 3),
                        get_child_res(1).content
                    ),
                ),
                DagOp::Neg => (
                    2,
                    format!(
                        "-{}",
                        get_child_str_with_parens(0, 2)
                    ),
                ),
                DagOp::Sqrt => (
                    4,
                    format!(
                        r"\sqrt{{{}}}",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Abs => (
                    10,
                    format!(
                        r"\left| {} \right|",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Pi => (
                    10,
                    r"\pi".to_string(),
                ),
                DagOp::E => (10, "e".to_string()),
                DagOp::Log => (
                    4,
                    format!(
                        r"\ln\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::LogBase => (
                    4,
                    format!(
                        r"\log_{{{}}}\left({}\right)",
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Exp => (
                    3,
                    format!(
                        "e^{{{}}}",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Sin => (
                    4,
                    format!(
                        r"\sin\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Cos => (
                    4,
                    format!(
                        r"\cos\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Tan => (
                    4,
                    format!(
                        r"\tan\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Csc => (
                    4,
                    format!(
                        r"\csc\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Sec => (
                    4,
                    format!(
                        r"\sec\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Cot => (
                    4,
                    format!(
                        r"\cot\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::ArcSin => (
                    4,
                    format!(
                        r"\arcsin\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::ArcCos => (
                    4,
                    format!(
                        r"\arccos\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::ArcTan => (
                    4,
                    format!(
                        r"\arctan\left({}\right)",
                        get_child_res(0).content
                    ),
                ),
                DagOp::Derivative(var) => (
                    5,
                    format!(
                        r"\frac{{d}}{{d{}}} \left[ {} \right]",
                        var,
                        get_child_res(0).content
                    ),
                ),
                DagOp::Integral => (
                    5,
                    format!(
                        r"\int_{{{}}}^{{{}}} {} \,d{}",
                        get_child_res(2).content,
                        get_child_res(3).content,
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Sum => (
                    5,
                    format!(
                        r"\sum_{{{}={}}}^{{{}}} {}",
                        get_child_res(1).content,
                        get_child_res(2).content,
                        get_child_res(3).content,
                        get_child_res(0).content
                    ),
                ),
                DagOp::Limit(var) => (
                    5,
                    format!(
                        r"\lim_{{{} \to {}}} {}",
                        to_greek(&var),
                        get_child_res(1).content,
                        get_child_res(0).content
                    ),
                ),
                DagOp::Eq => (
                    0,
                    format!(
                        "{} = {}",
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Binomial => (
                    5,
                    format!(
                        r"\binom{{{}}}{{{}}}",
                        get_child_res(0).content,
                        get_child_res(1).content
                    ),
                ),
                DagOp::Matrix { rows: _, cols } => {

                    let body = children
                        .chunks(cols)
                        .map(|row| {

                            row.iter()
                                .map(|elem| {

                                    results[&(elem as *const Expr)]
                                        .content
                                        .clone()
                                })
                                .collect::<Vec<_>>()
                                .join(" & ")
                        })
                        .collect::<Vec<_>>()
                        .join(r" \\ ");

                    (
                        10,
                        format!(
                            r"\begin{{pmatrix}}{}\end{{pmatrix}}",
                            body
                        ),
                    )
                }
                _ => (
                    10,
                    current_expr.to_string(),
                ),
            };

            results.insert(
                current_expr_ptr,
                LatexResult {
                    precedence: op_prec,
                    content: s,
                },
            );
        } else {

            for child in children
                .iter()
                .rev()
            {

                let cloned_child = child.clone();

                stack.push(cloned_child);
            }
        }
    }

    let final_result = &results[&(root_expr as *const Expr)];

    if final_result.precedence < root_precedence {

        format!(
            r"\left( {} \right)",
            final_result.content
        )
    } else {

        final_result
            .content
            .clone()
    }
}

/// Helper to add parentheses if needed. This function is now simplified as the main
/// iterative function handles most of the logic.

pub fn to_latex_prec_with_parens(
    expr: &Expr,
    precedence: u8,
) -> String {

    let op = expr.op();

    let op_prec = match op {
        DagOp::Add | DagOp::Sub => 1,
        _ => 10,
    };

    let s = to_latex_prec(expr, precedence);

    if op_prec < precedence {

        format!(
            r"\left( {} \right)",
            s
        )
    } else {

        s
    }
}

/// Converts common Greek letter names to LaTeX.

pub fn to_greek(s: &str) -> String {

    match s {
        "alpha" => r"\alpha".to_string(),
        "beta" => r"\beta".to_string(),
        "gamma" => r"\gamma".to_string(),
        "delta" => r"\delta".to_string(),
        "epsilon" => r"\epsilon".to_string(),
        "zeta" => r"\zeta".to_string(),
        "eta" => r"\eta".to_string(),
        "theta" => r"\theta".to_string(),
        "iota" => r"\iota".to_string(),
        "kappa" => r"\kappa".to_string(),
        "lambda" => r"\lambda".to_string(),
        "mu" => r"\mu".to_string(),
        "nu" => r"\nu".to_string(),
        "xi" => r"\xi".to_string(),
        "pi" => r"\pi".to_string(),
        "rho" => r"\rho".to_string(),
        "sigma" => r"\sigma".to_string(),
        "tau" => r"\tau".to_string(),
        "upsilon" => r"\upsilon".to_string(),
        "phi" => r"\phi".to_string(),
        "chi" => r"\chi".to_string(),
        "psi" => r"\psi".to_string(),
        "omega" => r"\omega".to_string(),
        _ => s.to_string(),
    }
}
