use crate::symbolic::core::Expr;
use num_traits::ToPrimitive;
use plotters::prelude::*;
use std::collections::HashMap;

/// Evaluates a symbolic expression to a numerical f64 value.
/// `vars` contains the numerical values for the variables in the expression.
pub(crate) fn eval_expr(expr: &Expr, vars: &HashMap<String, f64>) -> Result<f64, String> {
    match expr {
        Expr::Constant(c) => Ok(*c),
        Expr::BigInt(i) => Ok(i
            .to_f64()
            .ok_or_else(|| "BigInt conversion to f64 failed".to_string())?),
        Expr::Variable(v) => vars
            .get(v)
            .cloned()
            .ok_or_else(|| format!("Variable '{}' not found", v)),
        Expr::Add(a, b) => Ok(eval_expr(a, vars)? + eval_expr(b, vars)?),
        Expr::Sub(a, b) => Ok(eval_expr(a, vars)? - eval_expr(b, vars)?),
        Expr::Mul(a, b) => Ok(eval_expr(a, vars)? * eval_expr(b, vars)?),
        Expr::Div(a, b) => Ok(eval_expr(a, vars)? / eval_expr(b, vars)?),
        Expr::Power(b, e) => Ok(eval_expr(b, vars)?.powf(eval_expr(e, vars)?)),
        Expr::Neg(a) => Ok(-eval_expr(a, vars)?),
        Expr::Sqrt(a) => Ok(eval_expr(a, vars)?.sqrt()),
        Expr::Abs(a) => Ok(eval_expr(a, vars)?.abs()),
        Expr::Sin(a) => Ok(eval_expr(a, vars)?.sin()),
        Expr::Cos(a) => Ok(eval_expr(a, vars)?.cos()),
        Expr::Tan(a) => Ok(eval_expr(a, vars)?.tan()),
        Expr::Log(a) => Ok(eval_expr(a, vars)?.ln()),
        Expr::Exp(a) => Ok(eval_expr(a, vars)?.exp()),
        Expr::Pi => Ok(std::f64::consts::PI),
        Expr::E => Ok(std::f64::consts::E),
        _ => Err(format!(
            "Numerical evaluation for expression {:?} is not implemented",
            expr
        )),
    }
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

    // Determine y-range by sampling
    let y_min = (0..100)
        .map(|i| {
            let x = range.0 + (range.1 - range.0) * (i as f64 / 99.0);
            eval_expr(expr, &HashMap::from([(var.to_string(), x)]))
        })
        .filter_map(Result::ok)
        .fold(f64::INFINITY, f64::min);

    let y_max = (0..100)
        .map(|i| {
            let x = range.0 + (range.1 - range.0) * (i as f64 / 99.0);
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
                let x = range.0 + (range.1 - range.0) * (i as f64 / 500.0);
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
    comps: (&Expr, &Expr), // (vx, vy)
    vars: (&str, &str),    // (x, y)
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
            let x = x_range.0 + (x_range.1 - x_range.0) * (i as f64 / 19.0);
            let y = y_range.0 + (y_range.1 - y_range.0) * (j as f64 / 19.0);
            let mut vars_map = HashMap::new();
            vars_map.insert(x_var.to_string(), x);
            vars_map.insert(y_var.to_string(), y);

            if let (Ok(vx), Ok(vy)) = (eval_expr(vx_expr, &vars_map), eval_expr(vy_expr, &vars_map))
            {
                let magnitude = (vx * vx + vy * vy).sqrt();
                let end_x = x + vx / magnitude * (x_range.1 - x_range.0) * 0.05;
                let end_y = y + vy / magnitude * (y_range.1 - y_range.0) * 0.05;
                arrows.push(PathElement::new(vec![(x, y), (end_x, end_y)], BLUE));
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
    vars: (&str, &str), // (x, y)
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
    let _ = chart.draw_series(
        SurfaceSeries::xoz(
            // <-- SurfaceSeries::xoz() 调用开始
            // 参数 1: X 轴数据迭代器
            (0..100).map(|i| x_range.0 + (x_range.1 - x_range.0) * i as f64 / 99.0),
            // 参数 2: Z 轴数据迭代器
            (0..100).map(|i| y_range.0 + (y_range.1 - y_range.0) * i as f64 / 99.0),
            // 参数 3: 计算 Y 值的闭包 (Closure)
            |x, z| {
                let mut vars_map = HashMap::new();
                vars_map.insert(x_var.to_string(), x);
                vars_map.insert(y_var.to_string(), z);
                eval_expr(expr, &vars_map).unwrap_or(0.0)
            },
        ), // <-- SurfaceSeries::xoz() 调用结束
           // 整个 SurfaceSeries::xoz(...) 的结果被传递给 draw_series()
    );
    // 备注：根据上下文，这里可能还需要一个 .unwrap() 或 ? 来处理 Result
    root.present().map_err(|e| e.to_string())?;
    Ok(())
}

/// Plots a 3D parametric curve (x(t), y(t), z(t)) and saves it to a file.
pub fn plot_parametric_curve_3d(
    comps: (&Expr, &Expr, &Expr), // (x(t), y(t), z(t))
    var: &str,                    // t
    range: (f64, f64),
    path: &str,
) -> Result<(), String> {
    let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| e.to_string())?;

    // For now, we use fixed ranges. A better implementation would auto-fit.
    let mut chart = ChartBuilder::on(&root)
        .caption("3D Parametric Curve", ("sans-serif", 40).into_font())
        .build_cartesian_3d(-3.0..3.0, -3.0..3.0, -3.0..3.0)
        .map_err(|e| e.to_string())?;

    chart.configure_axes().draw().map_err(|e| e.to_string())?;

    let (x_expr, y_expr, z_expr) = comps;
    chart
        .draw_series(LineSeries::new(
            (0..=1000).map(|i| {
                let t = range.0 + (range.1 - range.0) * (i as f64 / 1000.0);
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
    comps: (&Expr, &Expr, &Expr),                 // (vx, vy, vz)
    vars: (&str, &str, &str),                     // (x, y, z)
    ranges: ((f64, f64), (f64, f64), (f64, f64)), // (x_range, y_range, z_range)
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
                let x = x_range.0 + (x_range.1 - x_range.0) * (i as f64 / (n_steps - 1) as f64);
                let y = y_range.0 + (y_range.1 - y_range.0) * (j as f64 / (n_steps - 1) as f64);
                let z = z_range.0 + (z_range.1 - z_range.0) * (k as f64 / (n_steps - 1) as f64);

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
                            vec![(x, y, z), (end_x, end_y, end_z)],
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
