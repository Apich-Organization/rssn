use crate::symbolic::core::{DagOp, Expr};
use num_traits::ToPrimitive;
use plotters::prelude::*;
use std::collections::HashMap;

/// Evaluates a symbolic expression to a numerical f64 value.
/// `vars` contains the numerical values for the variables in the expression.
/// This function is iterative to avoid stack overflows.

pub(crate) fn eval_expr(root_expr: &Expr, vars: &HashMap<String, f64>) -> Result<f64, String> {

    let mut results: HashMap<*const Expr, f64> = HashMap::new();

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

            let current_expr = stack.pop().expect("Value is valid");

            let current_expr_ptr = &current_expr as *const Expr;

            let get_child_val = |i: usize| -> f64 {

                results[&(&children[i] as *const Expr)]
            };

            let val_result = match current_expr.op() {
                DagOp::Constant(c) => Ok(c.into_inner()),
                DagOp::BigInt(i) => i
                    .to_f64()
                    .ok_or_else(|| "BigInt conversion to f64 failed".to_string()),
                DagOp::Variable(v) => vars
                    .get(&v)
                    .copied()
                    .ok_or_else(|| format!("Variable '{}' not found", v)),
                DagOp::Add => Ok(get_child_val(0) + get_child_val(1)),
                DagOp::Sub => Ok(get_child_val(0) - get_child_val(1)),
                DagOp::Mul => Ok(get_child_val(0) * get_child_val(1)),
                DagOp::Div => Ok(get_child_val(0) / get_child_val(1)),
                DagOp::Power => Ok(get_child_val(0).powf(get_child_val(1))),
                DagOp::Neg => Ok(-get_child_val(0)),
                DagOp::Sqrt => Ok(get_child_val(0).sqrt()),
                DagOp::Abs => Ok(get_child_val(0).abs()),
                DagOp::Sin => Ok(get_child_val(0).sin()),
                DagOp::Cos => Ok(get_child_val(0).cos()),
                DagOp::Tan => Ok(get_child_val(0).tan()),
                DagOp::Log => Ok(get_child_val(0).ln()),
                DagOp::Exp => Ok(get_child_val(0).exp()),
                DagOp::Pi => Ok(std::f64::consts::PI),
                DagOp::E => Ok(std::f64::consts::E),
                _ => Err(format!(
                    "Numerical evaluation for expression {:?} is not implemented",
                    current_expr
                )),
            };

            let val = val_result?;

            results.insert(current_expr_ptr, val);
        } else {

            for child in children.iter().rev() {

                if !results.contains_key(&(child as *const Expr)) {

                    let child_clone = child.clone();

                    stack.push(child_clone);
                }
            }
        }
    }

    Ok(results[&(root_expr as *const Expr)])
}

/// Plots a 2D function y = f(x) and saves it to a file.

pub fn plot_function_2d(
    expr: &Expr,
    var: &str,
    range: (f64, f64),
    path: &str,
) -> Result<(), String> {

    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();

    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let y_min = (0..100)
        .map(|i| {

            let x = range.0 + (range.1 - range.0) * (f64::from(i) / 99.0);

            eval_expr(expr, &HashMap::from([(var.to_string(), x)]))
        })
        .filter_map(Result::ok)
        .fold(f64::INFINITY, f64::min);

    let y_max = (0..100)
        .map(|i| {

            let x = range.0 + (range.1 - range.0) * (f64::from(i) / 99.0);

            eval_expr(expr, &HashMap::from([(var.to_string(), x)]))
        })
        .filter_map(Result::ok)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption("y = f(x)", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(range.0..range.1, y_min..y_max)
        .map_err(|e| e.to_string())?;

    chart.configure_mesh().draw().map_err(|e| e.to_string())?;

    chart
        .draw_series(LineSeries::new(
            (0..=500).map(|i| {

                let x = range.0 + (range.1 - range.0) * (f64::from(i) / 500.0);

                let y = eval_expr(expr, &HashMap::from([(var.to_string(), x)])).unwrap_or(0.0);

                (x, y)
            }),
            &RED,
        ))
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 2D vector field and saves it to a file.

pub fn plot_vector_field_2d(
    comps: (&Expr, &Expr),
    vars: (&str, &str),
    x_range: (f64, f64),
    y_range: (f64, f64),
    path: &str,
) -> Result<(), String> {

    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();

    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Vector Field", ("sans-serif", 40).into_font())
        .build_cartesian_2d(x_range.0..x_range.1, y_range.0..y_range.1)
        .map_err(|e| e.to_string())?;

    chart.configure_mesh().draw().map_err(|e| e.to_string())?;

    let (vx_expr, vy_expr) = comps;

    let (x_var, y_var) = vars;

    let mut arrows = Vec::new();

    for i in 0..20 {

        for j in 0..20 {

            let x = x_range.0 + (x_range.1 - x_range.0) * (f64::from(i) / 19.0);

            let y = y_range.0 + (y_range.1 - y_range.0) * (f64::from(j) / 19.0);

            let mut vars_map = HashMap::new();

            vars_map.insert(x_var.to_string(), x);

            vars_map.insert(y_var.to_string(), y);

            if let (Ok(vx), Ok(vy)) = (eval_expr(vx_expr, &vars_map), eval_expr(vy_expr, &vars_map))
            {

                let magnitude = (vx * vx + vy * vy).sqrt();

                let end_x = x + vx / magnitude * (x_range.1 - x_range.0) * 0.05;

                let end_y = y + vy / magnitude * (y_range.1 - y_range.0) * 0.05;

                arrows.push(PathElement::new(
                    vec![
                        (x, y),
                        (end_x, end_y),
                    ],
                    BLUE,
                ));
            }
        }
    }

    chart
        .draw_series(arrows.into_iter())
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D surface z = f(x, y) and saves it to a file.

pub fn plot_surface_3d(
    expr: &Expr,
    vars: (&str, &str),
    x_range: (f64, f64),
    y_range: (f64, f64),
    path: &str,
) -> Result<(), String> {

    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();

    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let mut chart = ChartBuilder::on(&root)
        .caption("z = f(x, y)", ("sans-serif", 40).into_font())
        .build_cartesian_3d(x_range.0..x_range.1, -1.0..1.0, y_range.0..y_range.1)
        .map_err(|e| e.to_string())?;

    chart.configure_axes().draw().map_err(|e| e.to_string())?;

    let (x_var, y_var) = vars;

    let _ = chart.draw_series(SurfaceSeries::xoz(
        (0..100).map(|i| x_range.0 + (x_range.1 - x_range.0) * f64::from(i) / 99.0),
        (0..100).map(|i| y_range.0 + (y_range.1 - y_range.0) * f64::from(i) / 99.0),
        |x, z| {

            let mut vars_map = HashMap::new();

            vars_map.insert(x_var.to_string(), x);

            vars_map.insert(y_var.to_string(), z);

            eval_expr(expr, &vars_map).unwrap_or(0.0)
        },
    ));

    root.present().map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D parametric curve (x(t), y(t), z(t)) and saves it to a file.

pub fn plot_parametric_curve_3d(
    comps: (&Expr, &Expr, &Expr),
    var: &str,
    range: (f64, f64),
    path: &str,
) -> Result<(), String> {

    let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();

    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let mut chart = ChartBuilder::on(&root)
        .caption("3D Parametric Curve", ("sans-serif", 40).into_font())
        .build_cartesian_3d(-3.0..3.0, -3.0..3.0, -3.0..3.0)
        .map_err(|e| e.to_string())?;

    chart.configure_axes().draw().map_err(|e| e.to_string())?;

    let (x_expr, y_expr, z_expr) = comps;

    chart
        .draw_series(LineSeries::new(
            (0..=1000).map(|i| {

                let t = range.0 + (range.1 - range.0) * (f64::from(i) / 1000.0);

                let mut vars_map = HashMap::new();

                vars_map.insert(var.to_string(), t);

                let x = eval_expr(x_expr, &vars_map).unwrap_or(0.0);

                let y = eval_expr(y_expr, &vars_map).unwrap_or(0.0);

                let z = eval_expr(z_expr, &vars_map).unwrap_or(0.0);

                (x, y, z)
            }),
            &RED,
        ))
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D vector field and saves it to a file.

pub fn plot_vector_field_3d(
    comps: (&Expr, &Expr, &Expr),
    vars: (&str, &str, &str),
    ranges: ((f64, f64), (f64, f64), (f64, f64)),
    path: &str,
) -> Result<(), String> {

    let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();

    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let (x_range, y_range, z_range) = ranges;

    let mut chart = ChartBuilder::on(&root)
        .caption("3D Vector Field", ("sans-serif", 40).into_font())
        .build_cartesian_3d(
            x_range.0..x_range.1,
            y_range.0..y_range.1,
            z_range.0..z_range.1,
        )
        .map_err(|e| e.to_string())?;

    chart.configure_axes().draw().map_err(|e| e.to_string())?;

    let (vx_expr, vy_expr, vz_expr) = comps;

    let (x_var, y_var, z_var) = vars;

    let mut arrows = Vec::new();

    let n_steps = 10;

    for i in 0..n_steps {

        for j in 0..n_steps {

            for k in 0..n_steps {

                let x =
                    x_range.0 + (x_range.1 - x_range.0) * (f64::from(i) / f64::from(n_steps - 1));

                let y =
                    y_range.0 + (y_range.1 - y_range.0) * (f64::from(j) / f64::from(n_steps - 1));

                let z =
                    z_range.0 + (z_range.1 - z_range.0) * (f64::from(k) / f64::from(n_steps - 1));

                let mut vars_map = HashMap::new();

                vars_map.insert(x_var.to_string(), x);

                vars_map.insert(y_var.to_string(), y);

                vars_map.insert(z_var.to_string(), z);

                if let (Ok(vx), Ok(vy), Ok(vz)) = (
                    eval_expr(vx_expr, &vars_map),
                    eval_expr(vy_expr, &vars_map),
                    eval_expr(vz_expr, &vars_map),
                ) {

                    let magnitude = (vx * vx + vy * vy + vz * vz).sqrt();

                    if magnitude > 1e-6 {

                        let scale = (x_range.1 - x_range.0) * 0.05;

                        let end_x = x + vx / magnitude * scale;

                        let end_y = y + vy / magnitude * scale;

                        let end_z = z + vz / magnitude * scale;

                        arrows.push(PathElement::new(
                            vec![
                                (x, y, z),
                                (end_x, end_y, end_z),
                            ],
                            BLUE,
                        ));
                    }
                }
            }
        }
    }

    chart
        .draw_series(arrows.into_iter())
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())?;

    Ok(())
}
