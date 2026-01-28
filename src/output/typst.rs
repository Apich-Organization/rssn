use std::collections::HashMap;

use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

/// Converts an expression to a Typst string.

#[must_use]

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

    let mut results: HashMap<
        Expr,
        TypstResult,
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

            let get_child_res = |i: usize| -> &TypstResult {
                &results[&children[i]]
            };

            let get_child_str = |i: usize, prec: u8| -> String {
                let child_res = get_child_res(i);
                if child_res.precedence < prec {
                    format!("({})", child_res.content)
                } else {
                    child_res.content.clone()
                }
            };

            let (op_prec, s) = match current_expr.op() {
                DagOp::Constant(c) => (10, c.to_string()),
                DagOp::BigInt(i) => (10, i.to_string()),
                DagOp::Rational(r) => (10, format!("frac({}, {})", r.numer(), r.denom())),
                DagOp::Variable(s) => (10, s.clone()),
                DagOp::Add => (1, format!("{} + {}", get_child_str(0, 1), get_child_str(1, 1))),
                DagOp::Sub => (1, format!("{} - {}", get_child_str(0, 1), get_child_str(1, 2))),
                DagOp::Mul => (2, format!("{} dot {}", get_child_str(0, 2), get_child_str(1, 2))),
                DagOp::Div => (2, format!("frac({}, {})", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Power => (3, format!("{}^({})", get_child_str(0, 3), get_child_res(1).content)),
                DagOp::Neg => (2, format!("-{}", get_child_str(0, 2))),
                DagOp::Sqrt => (10, format!("sqrt({})", get_child_res(0).content)),
                DagOp::Abs => (10, format!("abs({})", get_child_res(0).content)),
                DagOp::Pi => (10, "pi".to_string()),
                DagOp::E => (10, "e".to_string()),
                DagOp::Log => (10, format!("ln({})", get_child_res(0).content)),
                DagOp::LogBase => (10, format!("log_({})({})", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Exp => (10, format!("exp({})", get_child_res(0).content)),
                DagOp::Sin => (10, format!("sin({})", get_child_res(0).content)),
                DagOp::Cos => (10, format!("cos({})", get_child_res(0).content)),
                DagOp::Tan => (10, format!("tan({})", get_child_res(0).content)),
                DagOp::Csc => (10, format!("csc({})", get_child_res(0).content)),
                DagOp::Sec => (10, format!("sec({})", get_child_res(0).content)),
                DagOp::Cot => (10, format!("cot({})", get_child_res(0).content)),
                DagOp::ArcSin => (10, format!("arcsin({})", get_child_res(0).content)),
                DagOp::ArcCos => (10, format!("arccos({})", get_child_res(0).content)),
                DagOp::ArcTan => (10, format!("arctan({})", get_child_res(0).content)),
                DagOp::Factorial => (10, format!("{}!", get_child_str(0, 10))),
                DagOp::Derivative(var) => (10, format!("diff / (diff {}) ({})", var, get_child_res(0).content)),
                DagOp::Integral => (10, format!("integral_({})^({}) {} dif {}", get_child_res(2).content, get_child_res(3).content, get_child_res(0).content, get_child_res(1).content)),
                DagOp::Sum => (10, format!("sum_({}={})^({}) {}", get_child_res(1).content, get_child_res(2).content, get_child_res(3).content, get_child_res(0).content)),
                DagOp::Summation(s) => (10, format!("sum_({}={})^({}) {}", s, get_child_res(1).content, get_child_res(2).content, get_child_res(0).content)),
                DagOp::Product(s) => (10, format!("product_({}={})^({}) {}", s, get_child_res(1).content, get_child_res(2).content, get_child_res(0).content)),
                DagOp::Limit(var) => (10, format!("limit_({} -> {}) {}", var, get_child_res(1).content, get_child_res(0).content)),
                DagOp::Eq => (0, format!("{} = {}", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Binomial => (10, format!("binom({}, {})", get_child_res(0).content, get_child_res(1).content)),
                DagOp::Matrix { rows: _, cols } => {
                    let body = children.chunks(cols).map(|row| {
                        row.iter().map(|elem| results[elem].content.clone()).collect::<Vec<_>>().join(", ")
                    }).collect::<Vec<_>>().join("; ");
                    (10, format!("mat({body})"))
                },
                DagOp::Vector => {
                    let body = children.iter().map(|elem| results[elem].content.clone()).collect::<Vec<_>>().join(", ");
                    (10, format!("vec({body})"))
                },
                _ => (10, current_expr.to_string()),
            };

            results.insert(
                current_expr,
                TypstResult {
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
            "({})",
            final_res.content
        )
    } else {

        final_res
            .content
            .clone()
    }
}

#[cfg(test)]

mod tests {

    use super::*;
    use crate::prelude::Expr;

    #[test]

    fn test_to_typst_basic() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let y = Expr::Variable(
            "y".to_string(),
        );

        // (x + y) * x   -> x * (x + y) after normalization
        let expr = Expr::new_mul(
            Expr::new_add(
                x.clone(),
                y.clone(),
            ),
            x.clone(),
        );

        let typst = to_typst(&expr);

        assert_eq!(
            typst,
            "$x dot (x + y)$"
        );
    }

    #[test]

    fn test_to_typst_fractions() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::new_div(
            x,
            Expr::new_constant(2.0),
        );

        let typst = to_typst(&expr);

        assert_eq!(
            typst,
            "$frac(x, 2)$"
        );
    }

    #[test]

    fn test_to_typst_integral() {

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

        let typst = to_typst(&expr);

        assert!(
            typst.contains("integral_")
        );

        assert!(
            typst.contains("x^(2)")
        );

        assert!(
            typst.contains("dif x")
        );
    }

    #[test]

    fn test_to_typst_matrix() {

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

        let typst = to_typst(&expr);

        assert!(typst.contains(
            "mat(1, a; b, 2)"
        ));
    }

    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]
        #[test]
        fn prop_no_panic_to_typst(
            depth in 1..4usize,
        ) {
            let x = Expr::Variable("x".to_string());
            let mut expr = x.clone();
            for _ in 0..depth {
                expr = Expr::new_add(expr.clone(), Expr::new_constant(1.0));
                expr = Expr::new_div(expr.clone(), Expr::new_sqrt(x.clone()));
                expr = Expr::new_pow(expr.clone(), Expr::new_constant(2.0));
            }
            let typst = to_typst(&expr);
            assert!(!typst.is_empty());
        }
    }
}
