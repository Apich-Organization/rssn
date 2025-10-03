use crate::symbolic::core::Expr;

/// Represents a 2D box of characters for pretty-printing.
#[derive(Debug, Clone)]
pub struct PrintBox {
    width: usize,
    height: usize,
    lines: Vec<String>,
}

/// Main entry point. Converts an expression to a formatted string.
pub fn pretty_print(expr: &Expr) -> String {
    let root_box = to_box(expr);
    root_box.lines.join("\n")
}

/// Recursively converts an expression to a PrintBox.
pub(crate) fn to_box(expr: &Expr) -> PrintBox {
    match expr {
        Expr::Constant(c) => {
            let s = c.to_string();
            PrintBox {
                width: s.len(),
                height: 1,
                lines: vec![s],
            }
        }
        Expr::BigInt(i) => {
            let s = i.to_string();
            PrintBox {
                width: s.len(),
                height: 1,
                lines: vec![s],
            }
        }
        Expr::Variable(v) => PrintBox {
            width: v.len(),
            height: 1,
            lines: vec![v.clone()],
        },
        Expr::Add(a, b) => combine_horizontal(to_box(a), to_box(b), " + "),
        Expr::Sub(a, b) => combine_horizontal(to_box(a), to_box(b), " - "),
        Expr::Mul(a, b) => combine_horizontal(to_box(a), to_box(b), " * "),
        /*
        Expr::Div(num, den) => {
            let num_box = to_box(num);
            let den_box = to_box(den);
            let width = num_box.width.max(den_box.width);
            let bar = "─".repeat(width);

            let mut lines = Vec::new();
            for line in &num_box.lines {
                lines.push(center_text(line, width));
            }
            lines.push(bar);
            for line in &den_box.lines {
                lines.push(center_text(line, width));
            }

            PrintBox { width, height: lines.len(), lines }
        }
        */
        Expr::Div(num, den) => {
            let num_box = to_box(num);
            let den_box = to_box(den);
            let width = num_box.width.max(den_box.width) + 2; // Add padding
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
        }
        Expr::Power(base, exp) => {
            let base_box = to_box(base);
            let exp_box = to_box(exp);

            let new_width = base_box.width + exp_box.width;
            let new_height = base_box.height + exp_box.height;
            let mut lines = vec!["".to_string(); new_height];

            // Place exponent at the top right
            for i in 0..exp_box.height {
                lines[i] = format!("{}{}", " ".repeat(base_box.width), exp_box.lines[i]);
            }
            // Place base below
            for i in 0..base_box.height {
                lines[i + exp_box.height] =
                    format!("{}{}", base_box.lines[i], " ".repeat(exp_box.width));
            }

            PrintBox {
                width: new_width,
                height: new_height,
                lines,
            }
        }
        Expr::Sqrt(inner) => {
            let inner_box = to_box(inner);
            let width = inner_box.width + 1;
            let roof = "¯".repeat(width);
            let mut lines = vec![format!("√{}", roof)];
            for line in inner_box.lines {
                lines.push(format!(" {} ", line));
            }
            PrintBox {
                width: width + 1,
                height: lines.len(),
                lines,
            }
        }
        Expr::Integral {
            integrand,
            var,
            lower_bound,
            upper_bound,
        } => {
            let integrand_box = to_box(integrand);
            let upper_box = to_box(upper_bound);
            let lower_box = to_box(lower_bound);
            let var_box = to_box(var);

            let bounds_width = upper_box.width.max(lower_box.width);
            let int_height = integrand_box
                .height
                .max(upper_box.height + lower_box.height + 1);

            let mut lines = vec!["".to_string(); int_height];
            let integral_symbol = build_symbol('∫', int_height);

            // Combine bounds and integral symbol
            let upper_padded = center_text(&upper_box.lines[0], bounds_width);
            let lower_padded = center_text(&lower_box.lines[0], bounds_width);
            lines[0] = format!("{} {}", upper_padded, integral_symbol[0]);
            lines[int_height - 1] = format!("{} {}", lower_padded, integral_symbol[int_height - 1]);
            for i in 1..int_height - 1 {
                lines[i] = format!("{} {}", " ".repeat(bounds_width), integral_symbol[i]);
            }

            // Combine with integrand and dx
            for i in 0..int_height {
                let integrand_line = integrand_box.lines.get(i).cloned().unwrap_or_default();
                lines[i] = format!("{} {} d{}", lines[i], integrand_line, var_box.lines[0]);
            }

            PrintBox {
                width: lines[0].len(),
                height: lines.len(),
                lines,
            }
        }
        /*
        Expr::Sum { body, var, from, to } => {
            // Logic is very similar to Integral
            let body_box = to_box(body);
            let upper_box = to_box(to);
            let lower_box = to_box(from);
            let var_box = to_box(var);

            let bounds_width = upper_box.width.max(lower_box.width);
            let sum_height = body_box.height.max(upper_box.height + lower_box.height + 1);
            let mut lines = vec!["".to_string(); sum_height];
            let sum_symbol = build_symbol('Σ', sum_height);

            let upper_padded = center_text(&upper_box.lines[0], bounds_width);
            let lower_padded = format!("{}={}", var_box.lines[0], lower_box.lines[0]);
            let lower_padded = center_text(&lower_padded, bounds_width);

            lines[0] = format!("{} {}", upper_padded, sum_symbol[0]);
            lines[sum_height - 1] = format!("{} {}", lower_padded, sum_symbol[sum_height - 1]);
            for i in 1..sum_height - 1 {
                lines[i] = format!("{} {}", " ".repeat(bounds_width), sum_symbol[i]);
            }

            for i in 0..sum_height {
                let body_line = body_box.lines.get(i).cloned().unwrap_or_default();
                lines[i] = format!("{} {}", lines[i], body_line);
            }

            PrintBox { width: lines[0].len(), height: lines.len(), lines }
        }
        */
        Expr::Sum {
            body,
            var,
            from,
            to,
        } => {
            // Logic is very similar to Integral
            let body_box = to_box(body);
            let upper_box = to_box(to);
            let lower_box = to_box(from);
            let var_box = to_box(var);

            let bounds_width = upper_box.width.max(lower_box.width);
            let sum_height = body_box.height.max(upper_box.height + lower_box.height + 1);
            let mut lines = vec!["".to_string(); sum_height];
            let sum_symbol = build_symbol('Σ', sum_height);

            let upper_padded = center_text(&upper_box.lines[0], bounds_width);
            let lower_padded = format!("{}={}", var_box.lines[0], lower_box.lines[0]);
            let lower_padded = center_text(&lower_padded, bounds_width);

            lines[0] = format!("{} {}", upper_padded, sum_symbol[0]);
            lines[sum_height - 1] = format!("{} {}", lower_padded, sum_symbol[sum_height - 1]);
            for i in 1..sum_height - 1 {
                lines[i] = format!("{} {}", " ".repeat(bounds_width), sum_symbol[i]);
            }

            for i in 0..sum_height {
                let body_line = body_box.lines.get(i).cloned().unwrap_or_default();
                lines[i] = format!("{} {}", lines[i], body_line);
            }

            PrintBox {
                width: lines[0].len(),
                height: lines.len(),
                lines,
            }
        }
        Expr::Pi => to_box(&Expr::Variable("π".to_string())),
        Expr::E => to_box(&Expr::Variable("e".to_string())),
        Expr::Eq(a, b) => combine_horizontal(to_box(a), to_box(b), " = "),
        Expr::Abs(inner) => wrap_in_parens(to_box(inner), '|', '|'),
        Expr::Sin(arg) => wrap_in_parens(to_box(arg), '(', ')').prefix("sin"),
        Expr::Cos(arg) => wrap_in_parens(to_box(arg), '(', ')').prefix("cos"),
        Expr::Tan(arg) => wrap_in_parens(to_box(arg), '(', ')').prefix("tan"),
        Expr::Log(arg) => wrap_in_parens(to_box(arg), '(', ')').prefix("log"),
        Expr::Exp(arg) => to_box(&Expr::Power(Box::new(Expr::E), arg.clone())),
        // Fallback for other expression types
        _ => {
            let s = expr.to_string(); // Fallback to default Display
            PrintBox {
                width: s.len(),
                height: 1,
                lines: vec![s],
            }
        }
    }
}

/// Horizontally combines two PrintBoxes with an operator in between.
pub(crate) fn combine_horizontal(box_a: PrintBox, box_b: PrintBox, op: &str) -> PrintBox {
    let new_height = box_a.height.max(box_b.height);
    let mut lines = vec!["".to_string(); new_height];

    for i in 0..new_height {
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
        lines[i] = format!("{}{}{}", line_a, op, line_b);
    }

    PrintBox {
        width: box_a.width + op.len() + box_b.width,
        height: new_height,
        lines,
    }
}

/// Centers a line of text within a given width.
pub(crate) fn center_text(text: &str, width: usize) -> String {
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
pub(crate) fn wrap_in_parens(inner_box: PrintBox, open: char, close: char) -> PrintBox {
    let mut lines = Vec::new();
    for (i, line) in inner_box.lines.iter().enumerate() {
        if i == inner_box.height / 2 {
            lines.push(format!("{}{}{}", open, line, close));
        } else {
            lines.push(format!(" {}{}{}", " ", line, " "));
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
    pub(crate) fn prefix(mut self, s: &str) -> Self {
        self.lines[0] = format!("{}{}", s, self.lines[0]);
        self.width += s.len();
        self
    }
}

/// Builds a multi-line symbol.
pub(crate) fn build_symbol(symbol: char, height: usize) -> Vec<String> {
    if height <= 1 {
        return vec![symbol.to_string()];
    }
    let mut lines = vec!["".to_string(); height];
    let mid = height / 2;

    if symbol == '∫' {
        lines[0] = "⌠".to_string();
        for i in 1..height - 1 {
            lines[i] = "⎮".to_string();
        }
        lines[height - 1] = "⌡".to_string();
    } else if symbol == 'Σ' {
        lines[0] = "┌".to_string();
        for i in 1..height - 1 {
            lines[i] = "│".to_string();
        }
        lines[height - 1] = "└".to_string();
    } else {
        for i in 0..height {
            if i == mid {
                lines[i] = symbol.to_string();
            } else {
                lines[i] = " ".to_string();
            }
        }
    }
    lines
}
