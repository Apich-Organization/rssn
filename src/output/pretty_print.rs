use std::collections::HashMap;

use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

/// Represents a 2D box of characters for pretty-printing.
#[derive(Debug, Clone)]

pub struct PrintBox {
    width : usize,
    height : usize,
    lines : Vec<String>,
}

/// Main entry point. Converts an expression to a formatted string.

pub fn pretty_print(expr : &Expr) -> String {

    let root_box = to_box(expr);

    root_box
        .lines
        .join("\n")
}

/// Converts an expression to a PrintBox iteratively.

pub(crate) fn to_box(root_expr : &Expr) -> PrintBox {

    let mut results : HashMap<*const Expr, PrintBox> = HashMap::new();

    let mut stack : Vec<Expr> = vec![root_expr.clone()];

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

            let get_child_box = |i : usize| -> &PrintBox {

                &results[&(&children[i] as *const Expr)]
            };

            let val = match current_expr.op() {
                | DagOp::Constant(c) => {

                    let s = c
                        .into_inner()
                        .to_string();

                    PrintBox {
                        width : s.len(),
                        height : 1,
                        lines : vec![s],
                    }
                },
                | DagOp::BigInt(i) => {

                    let s = i.to_string();

                    PrintBox {
                        width : s.len(),
                        height : 1,
                        lines : vec![s],
                    }
                },
                | DagOp::Variable(v) => {
                    PrintBox {
                        width : v.len(),
                        height : 1,
                        lines : vec![v.clone()],
                    }
                },
                | DagOp::Add => {
                    combine_horizontal(
                        get_child_box(0),
                        get_child_box(1),
                        " + ",
                    )
                },
                | DagOp::Sub => {
                    combine_horizontal(
                        get_child_box(0),
                        get_child_box(1),
                        " - ",
                    )
                },
                | DagOp::Mul => {
                    combine_horizontal(
                        get_child_box(0),
                        get_child_box(1),
                        " * ",
                    )
                },
                | DagOp::Div => {

                    let num_box = get_child_box(0);

                    let den_box = get_child_box(1);

                    let width = num_box
                        .width
                        .max(den_box.width)
                        + 2;

                    let bar = "─".repeat(width);

                    let mut lines = Vec::new();

                    for line in &num_box.lines {

                        lines.push(center_text(
                            line, width,
                        ));
                    }

                    lines.push(bar);

                    for line in &den_box.lines {

                        lines.push(center_text(
                            line, width,
                        ));
                    }

                    PrintBox {
                        width,
                        height : lines.len(),
                        lines,
                    }
                },
                | DagOp::Power => {

                    let base_box = get_child_box(0);

                    let exp_box = get_child_box(1);

                    let new_width = base_box.width + exp_box.width;

                    let new_height = base_box.height + exp_box.height;

                    let mut lines = vec![String::new(); new_height];

                    for (i, l) in lines
                        .iter_mut()
                        .enumerate()
                        .take(exp_box.height)
                    {

                        *l = format!(
                            "{}{}",
                            " ".repeat(base_box.width),
                            exp_box.lines[i]
                        );
                    }

                    for i in 0 .. base_box.height {

                        lines[i + exp_box.height] = format!(
                            "{}{}",
                            base_box.lines[i],
                            " ".repeat(exp_box.width)
                        );
                    }

                    PrintBox {
                        width : new_width,
                        height : new_height,
                        lines,
                    }
                },
                | DagOp::Sqrt => {

                    let inner_box = get_child_box(0);

                    let width = inner_box.width + 1;

                    let roof = "¯".repeat(width);

                    let mut lines = vec![format!("√{}", roof)];

                    for line in &inner_box.lines {

                        lines.push(format!(
                            " {} ",
                            line
                        ));
                    }

                    PrintBox {
                        width : width + 1,
                        height : lines.len(),
                        lines,
                    }
                },
                | DagOp::Integral => {

                    let integrand_box = get_child_box(0);

                    let var_box = get_child_box(1);

                    let lower_box = get_child_box(2);

                    let upper_box = get_child_box(3);

                    let bounds_width = upper_box
                        .width
                        .max(lower_box.width);

                    let int_height = integrand_box
                        .height
                        .max(upper_box.height + lower_box.height + 1);

                    let mut lines = vec![String::new(); int_height];

                    let integral_symbol = build_symbol('∫', int_height);

                    let upper_padded = center_text(
                        &upper_box.lines[0],
                        bounds_width,
                    );

                    let lower_padded = center_text(
                        &lower_box.lines[0],
                        bounds_width,
                    );

                    lines[0] = format!(
                        "{} {}",
                        upper_padded, integral_symbol[0]
                    );

                    lines[int_height - 1] = format!(
                        "{} {}",
                        lower_padded,
                        integral_symbol[int_height - 1]
                    );

                    for i in 1 .. int_height - 1 {

                        lines[i] = format!(
                            "{} {}",
                            " ".repeat(bounds_width),
                            integral_symbol[i]
                        );
                    }

                    for (i, l) in lines
                        .iter_mut()
                        .enumerate()
                        .take(int_height)
                    {

                        let integrand_line = integrand_box
                            .lines
                            .get(i)
                            .cloned()
                            .unwrap_or_default();

                        *l = format!(
                            "{} {} d{}",
                            l, integrand_line, var_box.lines[0]
                        );
                    }

                    PrintBox {
                        width : lines[0].len(),
                        height : lines.len(),
                        lines,
                    }
                },
                | DagOp::Sum => {

                    let body_box = get_child_box(0);

                    let var_box = get_child_box(1);

                    let from_box = get_child_box(2);

                    let to_box = get_child_box(3);

                    let bounds_width = to_box
                        .width
                        .max(from_box.width + var_box.width + 1);

                    let sum_height = body_box
                        .height
                        .max(to_box.height + from_box.height + 1);

                    let mut lines = vec![String::new(); sum_height];

                    let sum_symbol = build_symbol('Σ', sum_height);

                    let upper_padded = center_text(
                        &to_box.lines[0],
                        bounds_width,
                    );

                    let lower_padded = format!(
                        "{}={}",
                        var_box.lines[0], from_box.lines[0]
                    );

                    let lower_padded = center_text(
                        &lower_padded,
                        bounds_width,
                    );

                    lines[0] = format!(
                        "{} {}",
                        upper_padded, sum_symbol[0]
                    );

                    lines[sum_height - 1] = format!(
                        "{} {}",
                        lower_padded,
                        sum_symbol[sum_height - 1]
                    );

                    for i in 1 .. sum_height - 1 {

                        lines[i] = format!(
                            "{} {}",
                            " ".repeat(bounds_width),
                            sum_symbol[i]
                        );
                    }

                    for (i, l) in lines
                        .iter_mut()
                        .enumerate()
                        .take(sum_height)
                    {

                        let body_line = body_box
                            .lines
                            .get(i)
                            .cloned()
                            .unwrap_or_default();

                        *l = format!(
                            "{} {}",
                            l, body_line
                        );
                    }

                    PrintBox {
                        width : lines[0].len(),
                        height : lines.len(),
                        lines,
                    }
                },
                | DagOp::Pi => {
                    PrintBox {
                        width : 1,
                        height : 1,
                        lines : vec!["π".to_string()],
                    }
                },
                | DagOp::E => {
                    PrintBox {
                        width : 1,
                        height : 1,
                        lines : vec!["e".to_string()],
                    }
                },
                | DagOp::Eq => {
                    combine_horizontal(
                        get_child_box(0),
                        get_child_box(1),
                        " = ",
                    )
                },
                | DagOp::Abs => {
                    wrap_in_parens(
                        get_child_box(0),
                        '|',
                        '|',
                    )
                },
                | DagOp::Sin => {
                    wrap_in_parens(
                        get_child_box(0),
                        '(',
                        ')',
                    )
                    .prefix("sin")
                },
                | DagOp::Cos => {
                    wrap_in_parens(
                        get_child_box(0),
                        '(',
                        ')',
                    )
                    .prefix("cos")
                },
                | DagOp::Tan => {
                    wrap_in_parens(
                        get_child_box(0),
                        '(',
                        ')',
                    )
                    .prefix("tan")
                },
                | DagOp::Log => {
                    wrap_in_parens(
                        get_child_box(0),
                        '(',
                        ')',
                    )
                    .prefix("log")
                },
                | DagOp::Exp => {

                    let base_box = PrintBox {
                        width : 1,
                        height : 1,
                        lines : vec!["e".to_string()],
                    };

                    let exp_box = get_child_box(0);

                    let new_width = base_box.width + exp_box.width;

                    let new_height = base_box.height + exp_box.height;

                    let mut lines = vec![String::new(); new_height];

                    for (i, l) in lines
                        .iter_mut()
                        .enumerate()
                        .take(exp_box.height)
                    {

                        *l = format!(
                            "{}{}",
                            " ".repeat(base_box.width),
                            exp_box.lines[i]
                        );
                    }

                    for i in 0 .. base_box.height {

                        lines[i + exp_box.height] = format!(
                            "{}{}",
                            base_box.lines[i],
                            " ".repeat(exp_box.width)
                        );
                    }

                    PrintBox {
                        width : new_width,
                        height : new_height,
                        lines,
                    }
                },
                | _ => {

                    let s = current_expr.to_string();

                    PrintBox {
                        width : s.len(),
                        height : 1,
                        lines : vec![s],
                    }
                },
            };

            results.insert(
                current_expr_ptr,
                val,
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

    results
        .remove(&(root_expr as *const Expr))
        .expect("Remove Failed")
}

/// Horizontally combines two PrintBoxes with an operator in between.

pub(crate) fn combine_horizontal(
    box_a : &PrintBox,
    box_b : &PrintBox,
    op : &str,
) -> PrintBox {

    let new_height = box_a
        .height
        .max(box_b.height);

    let mut lines = vec![String::new(); new_height];

    for (i, vars) in lines
        .iter_mut()
        .enumerate()
        .take(new_height)
    {

        let line_a = box_a
            .lines
            .get(i)
            .cloned()
            .unwrap_or_else(|| " ".repeat(box_a.width));

        let line_b = box_b
            .lines
            .get(i)
            .cloned()
            .unwrap_or_else(|| " ".repeat(box_b.width));

        *vars = format!(
            "{}{}{}",
            line_a, op, line_b
        );
    }

    PrintBox {
        width : box_a.width + op.len() + box_b.width,
        height : new_height,
        lines,
    }
}

/// Centers a line of text within a given width.

pub(crate) fn center_text(
    text : &str,
    width : usize,
) -> String {

    let padding = width.saturating_sub(text.len());

    let left_padding = padding / 2;

    let right_padding = padding - left_padding;

    format!(
        "{}{}{}",
        " ".repeat(left_padding),
        text,
        " ".repeat(right_padding)
    )
}

/// Wraps a PrintBox in parentheses or other delimiters.

pub(crate) fn wrap_in_parens(
    inner_box : &PrintBox,
    open : char,
    close : char,
) -> PrintBox {

    let mut lines = Vec::new();

    for (i, line) in inner_box
        .lines
        .iter()
        .enumerate()
    {

        if i == inner_box.height / 2 {

            lines.push(format!(
                "{}{}{}",
                open, line, close
            ));
        } else {

            lines.push(format!(
                " {}{}{}",
                " ", line, " "
            ));
        }
    }

    PrintBox {
        width : inner_box.width + 2,
        height : inner_box.height,
        lines,
    }
}

impl PrintBox {
    /// Adds a prefix string to the first line of the box.

    pub(crate) fn prefix(
        mut self,
        s : &str,
    ) -> Self {

        self.lines[0] = format!(
            "{}{}",
            s, self.lines[0]
        );

        self.width += s.len();

        self
    }
}

/// Builds a multi-line symbol.

pub(crate) fn build_symbol(
    symbol : char,
    height : usize,
) -> Vec<String> {

    if height <= 1 {

        return vec![symbol.to_string()];
    }

    let mut lines = vec![String::new(); height];

    let mid = height / 2;

    if symbol == '∫' {

        lines[0] = "⌠".to_string();

        for i in lines
            .iter_mut()
            .take(height - 1)
            .skip(1)
        {

            *i = "⎮".to_string();
        }

        lines[height - 1] = "⌡".to_string();
    } else if symbol == 'Σ' {

        lines[0] = "┌".to_string();

        for i in lines
            .iter_mut()
            .take(height - 1)
            .skip(1)
        {

            *i = "│".to_string();
        }

        lines[height - 1] = "└".to_string();
    } else {

        for (i, vars) in lines
            .iter_mut()
            .enumerate()
            .take(height)
        {

            if i == mid {

                *vars = symbol.to_string();
            } else {

                *vars = " ".to_string();
            }
        }
    }

    lines
}
