use crate::symbolic::core::{
    DagOp,
    Expr,
};
use std::collections::HashMap;

/// Converts an expression to a Typst string.

pub fn to_typst(expr: &Expr) -> String {

    format!(
        "${}$",
        to_typst_prec(expr, 0)
    )
}

#[derive(Clone)]

struct TypstResult {
    precedence: u8,
    content: String,
}

/// Converts an expression to a Typst string with precedence handling.
/// This function is iterative to avoid stack overflows.

pub(crate) fn to_typst_prec(
    root_expr: &Expr,
    root_precedence: u8,
) -> String {

    let mut results: HashMap<*const Expr, TypstResult> = HashMap::new();

    let mut stack: Vec<Expr> = vec![root_expr.clone()];

    while let Some(expr) = stack.last() {

        let expr_ptr = &*expr as *const Expr;

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

            let get_child_res = |i: usize| -> &TypstResult {

                &results[&(&children[i] as *const Expr)]
            };

            let get_child_str = |i: usize, prec: u8| -> String {

                let child_res = get_child_res(i);

                if child_res.precedence < prec {

                    format!(
                        "({})",
                        child_res.content
                    )
                } else {

                    child_res
                        .content
                        .clone()
                }
            };

            let (op_prec, s) = match current_expr.op() {
                | DagOp::Add => {
                    (
                        1,
                        format!(
                            "{} + {}",
                            get_child_str(0, 1),
                            get_child_str(1, 1)
                        ),
                    )
                },
                | DagOp::Sub => {
                    (
                        1,
                        format!(
                            "{} - {}",
                            get_child_str(0, 1),
                            get_child_str(1, 2)
                        ),
                    )
                },
                | DagOp::Mul => {
                    (
                        2,
                        format!(
                            "{} * {}",
                            get_child_str(0, 2),
                            get_child_str(1, 2)
                        ),
                    )
                },
                | DagOp::Div => {
                    (
                        2,
                        format!(
                            "frac({}, {})",
                            get_child_res(0).content,
                            get_child_res(1).content
                        ),
                    )
                },
                | DagOp::Power => {
                    (
                        3,
                        format!(
                            "{}^({})",
                            get_child_str(0, 3),
                            get_child_res(1).content
                        ),
                    )
                },
                | DagOp::Neg => {
                    (
                        2,
                        format!(
                            "-{}",
                            get_child_str(0, 2)
                        ),
                    )
                },
                | DagOp::Sqrt => {
                    (
                        4,
                        format!(
                            "sqrt({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Sin => {
                    (
                        4,
                        format!(
                            "sin({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Cos => {
                    (
                        4,
                        format!(
                            "cos({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Tan => {
                    (
                        4,
                        format!(
                            "tan({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Log => {
                    (
                        4,
                        format!(
                            "ln({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Exp => {
                    (
                        3,
                        format!(
                            "e^({})",
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Integral => {
                    (
                        5,
                        format!(
                            "integral_({})({}) {} dif {}",
                            get_child_res(2).content,
                            get_child_res(3).content,
                            get_child_res(0).content,
                            get_child_res(1).content
                        ),
                    )
                },
                | DagOp::Sum => {
                    (
                        5,
                        format!(
                            "sum_({}={})({}) {}",
                            get_child_res(1).content,
                            get_child_res(2).content,
                            get_child_res(3).content,
                            get_child_res(0).content
                        ),
                    )
                },
                | DagOp::Matrix { rows: _, cols } => {

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
                                .join(", ")
                        })
                        .collect::<Vec<_>>()
                        .join("; ");

                    (
                        10,
                        format!("mat({})", body),
                    )
                },
                | DagOp::Pi => (10, "pi".to_string()),
                | DagOp::E => (10, "e".to_string()),
                | DagOp::Constant(c) => {
                    (
                        10,
                        c.into_inner()
                            .to_string(),
                    )
                },
                | DagOp::BigInt(i) => (10, i.to_string()),
                | DagOp::Variable(s) => (10, s.clone()),
                | _ => {
                    (
                        10,
                        current_expr.to_string(),
                    )
                },
            };

            results.insert(
                current_expr_ptr,
                TypstResult {
                    precedence: op_prec,
                    content: s,
                },
            );
        } else {

            for child in children
                .iter()
                .rev()
            {

                if !results.contains_key(&(child as *const Expr)) {

                    let child_clone = child.clone();

                    stack.push(child_clone);
                }
            }
        }
    }

    let final_result = &results[&(root_expr as *const Expr)];

    if final_result.precedence < root_precedence {

        format!(
            "({})",
            final_result.content
        )
    } else {

        final_result
            .content
            .clone()
    }
}
