use std::collections::HashMap;

use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

/// Represents a 2D box of characters for pretty-printing.
#[derive(Debug, Clone)]

pub struct PrintBox {
    width: usize,
    height: usize,
    lines: Vec<String>,
}

/// Main entry point. Converts an expression to a formatted string.

#[must_use]

pub fn pretty_print(
    expr: &Expr
) -> String {

    let root_box = to_box(expr);

    root_box
        .lines
        .join("\n")
}

/// Converts an expression to a `PrintBox` iteratively.

pub(crate) fn to_box(
    root_expr: &Expr
) -> PrintBox {

    let mut results: HashMap<
        Expr,
        PrintBox,
    > = HashMap::new();

    let mut stack: Vec<Expr> =
        vec![root_expr.clone()];

    while let Some(expr) = stack.last()
    {

        if results.contains_key(expr) {

            stack.pop();

            continue;
        }

        let children = expr.children();

        let all_children_processed =
            children
                .iter()
                .all(|c| {

                    results
                        .contains_key(c)
                });

        if all_children_processed {

            let current_expr =
                stack.pop().expect(
                    "Value is valid",
                );

            let get_child_box = |i: usize| -> &PrintBox {
                &results[&children[i]]
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

                    let width = num_box.width.max(den_box.width) + 2;
                    let bar = "─".repeat(width);

                    let mut lines = Vec::new();
                    for line in &num_box.lines {
                        lines.push(center_text(line, width));
                    }
                    lines.push(bar);
                    for line in &den_box.lines {
                        lines.push(center_text(line, width));
                    }

                    PrintBox {
                        width,
                        height: lines.len(),
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
                            " {line} "
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

                    let bounds_width = upper_box.width.max(lower_box.width);
                    let int_height = integrand_box.height.max(upper_box.height + lower_box.height + 1).max(3);

                    let mut symbol_lines = vec![String::new(); int_height];
                    let integral_symbol = build_symbol('∫', int_height);

                    let upper_padded = center_text(&upper_box.lines[0], bounds_width);
                    let lower_padded = center_text(&lower_box.lines[0], bounds_width);

                    symbol_lines[0] = format!("{} {}", upper_padded, integral_symbol[0]);
                    symbol_lines[int_height - 1] = format!("{} {}", lower_padded, integral_symbol[int_height - 1]);

                    for i in 1..int_height - 1 {
                        symbol_lines[i] = format!("{} {}", " ".repeat(bounds_width), integral_symbol[i]);
                    }

                    let symbol_box = PrintBox {
                        width: symbol_lines[0].len(),
                        height: int_height,
                        lines: symbol_lines,
                    };

                    let combined = combine_horizontal(&symbol_box, integrand_box, " ");
                    let var_str = format!(" d{}", var_box.lines[0]);
                    combined.suffix(&var_str)
                },
                | DagOp::Sum => {
                    let body_box = get_child_box(0);
                    let var_box = get_child_box(1);
                    let from_box = get_child_box(2);
                    let to_box = get_child_box(3);

                    let bounds_width = to_box.width.max(from_box.width + var_box.width + 1);
                    let sum_height = body_box.height.max(to_box.height + from_box.height + 1).max(3);

                    let mut symbol_lines = vec![String::new(); sum_height];
                    let sum_symbol = build_symbol('Σ', sum_height);

                    let upper_padded = center_text(&to_box.lines[0], bounds_width);
                    let lower_text = format!("{}={}", var_box.lines[0], from_box.lines[0]);
                    let lower_padded = center_text(&lower_text, bounds_width);

                    symbol_lines[0] = format!("{} {}", upper_padded, sum_symbol[0]);
                    symbol_lines[sum_height - 1] = format!("{} {}", lower_padded, sum_symbol[sum_height - 1]);

                    for i in 1..sum_height - 1 {
                        symbol_lines[i] = format!("{} {}", " ".repeat(bounds_width), sum_symbol[i]);
                    }

                    let symbol_box = PrintBox {
                        width: symbol_lines[0].len(),
                        height: sum_height,
                        lines: symbol_lines,
                    };

                    combine_horizontal(&symbol_box, body_box, " ")
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
                | DagOp::Neg => {
                    get_child_box(0).clone().prefix("-")
                },
                | DagOp::Vector => {
                    let mut elements = Vec::new();
                    for i in 0..children.len() {
                        elements.push(get_child_box(i));
                    }

                    if elements.is_empty() {
                         PrintBox { width: 2, height: 1, lines: vec!["[]".to_string()] }
                    } else {
                        let mut res = elements[0].clone();
                        for elem in elements.iter().skip(1) {
                            res = combine_horizontal(&res, elem, ", ");
                        }
                        wrap_in_parens(&res, '[', ']')
                    }
                },
                | DagOp::Matrix { rows, cols } => {
                    if rows == 0 || cols == 0 || children.is_empty() {
                         return PrintBox { width: 4, height: 1, lines: vec!["[[]]".to_string()] };
                    }
                    let mut grid: Vec<Vec<&PrintBox>> = vec![vec![&results[&children[0]]; cols]; rows];
                    for r in 0..rows {
                        for c in 0..cols {
                            grid[r][c] = &results[&children[r * cols + c]];
                        }
                    }

                    let mut col_widths = vec![0; cols];
                    for c in 0..cols {
                        for row in &grid {
                            col_widths[c] = col_widths[c].max(row[c].width);
                        }
                    }

                    let mut row_heights = vec![0; rows];
                    for r in 0..rows {
                        for cell in &grid[r] {
                            row_heights[r] = row_heights[r].max(cell.height);
                        }
                    }

                    let mut matrix_lines = Vec::new();
                    for r in 0..rows {
                        for sub_row in 0..row_heights[r] {
                            let mut line = String::new();
                            for c in 0..cols {
                                let cell = grid[r][c];
                                let cell_line = cell.lines.get(sub_row).cloned().unwrap_or_else(|| " ".repeat(cell.width));
                                line.push_str(&center_text(&cell_line, col_widths[c]));
                                if c < cols - 1 {
                                    line.push_str("  ");
                                }
                            }
                            matrix_lines.push(line);
                        }
                        if r < rows - 1 {
                            matrix_lines.push(" ".to_string());
                        }
                    }

                    let width = matrix_lines.first().map_or(0, std::string::String::len);
                    wrap_in_parens(&PrintBox { width, height: matrix_lines.len(), lines: matrix_lines }, '[', ']')
                },
                | DagOp::Binomial => {
                    let n = get_child_box(0);
                    let k = get_child_box(1);
                    let width = n.width.max(k.width);
                    let mut lines = Vec::new();
                    for l in &n.lines { lines.push(center_text(l, width)); }
                    for l in &k.lines { lines.push(center_text(l, width)); }
                    wrap_in_parens(&PrintBox { width, height: lines.len(), lines }, '(', ')')
                },
                | DagOp::Derivative(var) => {
                    let inner = get_child_box(0);
                    let pref = format!("d/d{var} ");
                    inner.clone().prefix(&pref)
                },
                | DagOp::Factorial => {
                    get_child_box(0).clone().suffix("!")
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
                current_expr,
                val,
            );
        } else {

            for child in children
                .iter()
                .rev()
            {

                if !results
                    .contains_key(child)
                {

                    stack.push(
                        child.clone(),
                    );
                }
            }
        }
    }

    results
        .remove(root_expr)
        .expect("Remove Failed")
}

/// Horizontally combines two `PrintBoxes` with an operator in between.

pub(crate) fn combine_horizontal(
    box_a: &PrintBox,
    box_b: &PrintBox,
    op: &str,
) -> PrintBox {

    let new_height = box_a
        .height
        .max(box_b.height);

    let mut lines =
        vec![String::new(); new_height];

    for (i, vars) in lines
        .iter_mut()
        .enumerate()
        .take(new_height)
    {

        let line_a = box_a
            .lines
            .get(i)
            .cloned()
            .unwrap_or_else(|| {

                " ".repeat(box_a.width)
            });

        let line_b = box_b
            .lines
            .get(i)
            .cloned()
            .unwrap_or_else(|| {

                " ".repeat(box_b.width)
            });

        *vars = format!(
            "{line_a}{op}{line_b}"
        );
    }

    PrintBox {
        width: box_a.width
            + op.len()
            + box_b.width,
        height: new_height,
        lines,
    }
}

/// Centers a line of text within a given width.

pub(crate) fn center_text(
    text: &str,
    width: usize,
) -> String {

    let padding = width
        .saturating_sub(text.len());

    let left_padding = padding / 2;

    let right_padding =
        padding - left_padding;

    format!(
        "{}{}{}",
        " ".repeat(left_padding),
        text,
        " ".repeat(right_padding)
    )
}

/// Wraps a `PrintBox` in parentheses or other delimiters.

pub(crate) fn wrap_in_parens(
    inner_box: &PrintBox,
    open: char,
    close: char,
) -> PrintBox {

    let mut lines = Vec::new();

    for (i, line) in inner_box
        .lines
        .iter()
        .enumerate()
    {

        if i == inner_box.height / 2 {

            lines.push(format!(
                "{open}{line}{close}"
            ));
        } else {

            lines.push(format!(
                " {} {} ",
                " ", line
            ));
        }
    }

    PrintBox {
        width: inner_box.width + 2,
        height: inner_box.height,
        lines,
    }
}

impl PrintBox {
    /// Adds a prefix string to the first line of the box.

    pub(crate) fn prefix(
        mut self,
        s: &str,
    ) -> Self {

        let padding =
            " ".repeat(s.len());

        let mid = self.height / 2;

        for (i, line) in self
            .lines
            .iter_mut()
            .enumerate()
        {

            if i == mid {

                *line = format!(
                    "{s}{line}"
                );
            } else {

                *line = format!(
                    "{padding}{line}"
                );
            }
        }

        self.width += s.len();

        self
    }

    pub(crate) fn suffix(
        mut self,
        s: &str,
    ) -> Self {

        let padding =
            " ".repeat(s.len());

        let mid = self.height / 2;

        for (i, line) in self
            .lines
            .iter_mut()
            .enumerate()
        {

            if i == mid {

                line.push_str(s);
            } else {

                line.push_str(&padding);
            }
        }

        self.width += s.len();

        self
    }
}

/// Builds a multi-line symbol.

pub(crate) fn build_symbol(
    symbol: char,
    height: usize,
) -> Vec<String> {

    if height <= 1 {

        return vec![
            symbol.to_string(),
        ];
    }

    let mut lines =
        vec![String::new(); height];

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

        lines[height - 1] =
            "⌡".to_string();
    } else if symbol == 'Σ' {

        lines[0] = "┌".to_string();

        for i in lines
            .iter_mut()
            .take(height - 1)
            .skip(1)
        {

            *i = "│".to_string();
        }

        lines[height - 1] =
            "└".to_string();
    } else {

        for (i, vars) in lines
            .iter_mut()
            .enumerate()
            .take(height)
        {

            if i == mid {

                *vars =
                    symbol.to_string();
            } else {

                *vars = " ".to_string();
            }
        }
    }

    lines
}

#[cfg(test)]

mod tests {

    use super::*;
    use crate::prelude::Expr;

    #[test]

    fn test_pretty_print_basic() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let two = Expr::Constant(2.0);

        let expr =
            Expr::new_add(x, two);

        let output =
            pretty_print(&expr);

        println!("{}", output);

        assert!(
            output.contains("2 + x")
        );
    }

    #[test]

    fn test_pretty_print_div() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let y = Expr::Variable(
            "y".to_string(),
        );

        let expr = Expr::new_div(x, y);

        let output =
            pretty_print(&expr);

        assert!(output.contains("─"));

        assert!(output.contains('x'));

        assert!(output.contains('y'));
    }

    #[test]

    fn test_pretty_print_sqrt() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::new_sqrt(x);

        let output =
            pretty_print(&expr);

        assert!(output.contains('√'));
    }

    #[test]

    fn test_pretty_print_matrix() {

        let expr =
            Expr::new_matrix(vec![
                vec![
                    Expr::Constant(1.0),
                    Expr::Constant(2.0),
                ],
                vec![
                    Expr::Constant(3.0),
                    Expr::Constant(4.0),
                ],
            ]);

        let output =
            pretty_print(&expr);

        assert!(output.contains("1"));

        assert!(output.contains("4"));

        assert!(output.contains('['));
    }

    #[test]

    fn test_pretty_print_integral() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::Integral {
            integrand:
                std::sync::Arc::new(
                    Expr::new_pow(
                        x.clone(),
                        Expr::Constant(
                            2.0,
                        ),
                    ),
                ),
            var: std::sync::Arc::new(
                x.clone(),
            ),
            lower_bound:
                std::sync::Arc::new(
                    Expr::Constant(0.0),
                ),
            upper_bound:
                std::sync::Arc::new(
                    Expr::Constant(1.0),
                ),
        };

        let output =
            pretty_print(&expr);

        println!("{}", output);

        assert!(output.contains('⎮'));

        assert!(output.contains("⌠"));

        assert!(output.contains("⌡"));

        assert!(output.contains('0'));

        assert!(output.contains('1'));

        assert!(output.contains("dx"));
    }

    #[test]

    fn test_pretty_print_factorial() {

        let n = Expr::Variable(
            "n".to_string(),
        );

        let expr = Expr::Factorial(
            std::sync::Arc::new(n),
        );

        let output =
            pretty_print(&expr);

        assert!(output.contains("n!"));
    }

    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]
        #[test]
        fn prop_no_panic_pretty_print(
            depth in 1..4usize,
            width in 1..3usize
        ) {
            // We'll generate a few representative expressions to ensure no panics
            let x = Expr::Variable("x".to_string());
            let y = Expr::Variable("y".to_string());
            let base_exprs = vec![x, y, Expr::Constant(1.0), Expr::Pi];

            // Just a sampling of possible structures
            let mut expr = base_exprs[0].clone();
            for _ in 0..depth {
                expr = Expr::new_add(expr, base_exprs[1].clone());
                expr = Expr::new_div(expr, Expr::Constant(2.0));
                expr = Expr::new_pow(expr, Expr::Constant(0.5));
            }

            let output = pretty_print(&expr);
            assert!(!output.is_empty());
        }
    }
}
