use std::collections::HashMap;

use num_traits::ToPrimitive;
use plotters::prelude::*;

use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

/// Configuration for plotting.
#[derive(Clone, Debug)]

pub struct PlotConfig {
    /// The width of the plot in pixels.
    pub width: u32,
    /// The height of the plot in pixels.
    pub height: u32,
    /// The caption of the plot.
    pub caption: String,
    /// The color of the line in the plot.
    pub line_color: RGBAColor,
    /// The color of the mesh in the plot.
    pub mesh_color: RGBAColor,
    /// The number of samples to use when plotting.
    pub samples: usize,
}

impl Default for PlotConfig {
    fn default() -> Self {

        Self {
            width: 800,
            height: 600,
            caption: "Plot".to_string(),
            line_color: RED.to_rgba(),
            mesh_color: BLACK.mix(0.1),
            samples: 500,
        }
    }
}

/// Evaluates a symbolic expression to a numerical f64 value.
/// `vars` contains the numerical values for the variables in the expression.
/// This function is iterative to avoid stack overflows.

pub(crate) fn eval_expr(
    root_expr: &Expr,
    vars: &HashMap<String, f64>,
) -> Result<f64, String> {

    let mut results: HashMap<
        Expr,
        f64,
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

            let get_child_val =
                |i: usize| -> f64 {

                    results
                        [&children[i]]
                };

            let val_result = match current_expr.op() {
                DagOp::Constant(c) => Ok(c.into_inner()),
                DagOp::BigInt(i) => i.to_f64().ok_or_else(|| "BigInt conversion to f64 failed".to_string()),
                DagOp::Rational(r) => Ok(r.numer().to_f64().unwrap() / r.denom().to_f64().unwrap()),
                DagOp::Variable(v) => vars.get(&v).copied().ok_or_else(|| format!("Variable '{v}' not found")),
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
                DagOp::Csc => Ok(1.0 / get_child_val(0).sin()),
                DagOp::Sec => Ok(1.0 / get_child_val(0).cos()),
                DagOp::Cot => Ok(1.0 / get_child_val(0).tan()),
                DagOp::ArcSin => Ok(get_child_val(0).asin()),
                DagOp::ArcCos => Ok(get_child_val(0).acos()),
                DagOp::ArcTan => Ok(get_child_val(0).atan()),
                DagOp::Sinh => Ok(get_child_val(0).sinh()),
                DagOp::Cosh => Ok(get_child_val(0).cosh()),
                DagOp::Tanh => Ok(get_child_val(0).tanh()),
                DagOp::Log => Ok(get_child_val(0).ln()),
                DagOp::LogBase => Ok(get_child_val(1).log(get_child_val(0))),
                DagOp::Exp => Ok(get_child_val(0).exp()),
                DagOp::Floor => Ok(get_child_val(0).floor()),
                DagOp::Pi => Ok(std::f64::consts::PI),
                DagOp::E => Ok(std::f64::consts::E),
                _ => Err(format!("Numerical evaluation for operation {:?} is not implemented", current_expr.op())),
            };

            let val = val_result?;

            results.insert(
                current_expr,
                val,
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

    Ok(results[root_expr])
}

/// Plots a 2D function y = f(x) and saves it to a file.
///
/// # Errors
///
/// This function will return an error if the plot cannot be created or saved.

pub fn plot_function_2d(
    expr: &Expr,
    var: &str,
    range: (f64, f64),
    path: &str,
    config: Option<PlotConfig>,
) -> Result<(), String> {

    let conf =
        config.unwrap_or_default();

    let root = BitMapBackend::new(
        path,
        (
            conf.width,
            conf.height,
        ),
    )
    .into_drawing_area();

    root.fill(&WHITE)
        .map_err(|e| e.to_string())?;

    let y_min = (0 .. 100)
        .map(|i| {

            let x = (range.1 - range.0)
                .mul_add(
                    f64::from(i) / 99.0,
                    range.0,
                );

            eval_expr(
                expr,
                &HashMap::from([(
                    var.to_string(),
                    x,
                )]),
            )
        })
        .filter_map(Result::ok)
        .fold(
            f64::INFINITY,
            f64::min,
        );

    let y_max = (0 .. 100)
        .map(|i| {

            let x = (range.1 - range.0)
                .mul_add(
                    f64::from(i) / 99.0,
                    range.0,
                );

            eval_expr(
                expr,
                &HashMap::from([(
                    var.to_string(),
                    x,
                )]),
            )
        })
        .filter_map(Result::ok)
        .fold(
            f64::NEG_INFINITY,
            f64::max,
        );

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                &conf.caption,
                ("sans-serif", 40)
                    .into_font(),
            )
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(
                range.0 .. range.1,
                y_min .. y_max,
            )
            .map_err(|e| {

                e.to_string()
            })?;

    chart
        .configure_mesh()
        .light_line_style(
            conf.mesh_color,
        )
        .draw()
        .map_err(|e| e.to_string())?;

    chart
        .draw_series(LineSeries::new(
            (0..=conf.samples).map(|i| {
                let x = (range.1 - range.0).mul_add((i as f64) / conf.samples as f64, range.0);
                let y = eval_expr(expr, &HashMap::from([(var.to_string(), x)])).unwrap_or(0.0);
                (x, y)
            }),
            &conf.line_color,
        ))
        .map_err(|e| e.to_string())?;

    root.present()
        .map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 2D vector field and saves it to a file.
///
/// # Errors
///
/// This function will return an error if the plot cannot be created or saved.

pub fn plot_vector_field_2d(
    comps: (&Expr, &Expr),
    vars: (&str, &str),
    x_range: (f64, f64),
    y_range: (f64, f64),
    path: &str,
    config: Option<PlotConfig>,
) -> Result<(), String> {

    let conf =
        config.unwrap_or_default();

    let root = BitMapBackend::new(
        path,
        (
            conf.width,
            conf.height,
        ),
    )
    .into_drawing_area();

    root.fill(&WHITE)
        .map_err(|e| e.to_string())?;

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                &conf.caption,
                ("sans-serif", 40)
                    .into_font(),
            )
            .build_cartesian_2d(
                x_range.0 .. x_range.1,
                y_range.0 .. y_range.1,
            )
            .map_err(|e| {

                e.to_string()
            })?;

    chart
        .configure_mesh()
        .light_line_style(
            conf.mesh_color,
        )
        .draw()
        .map_err(|e| e.to_string())?;

    let (vx_expr, vy_expr) = comps;

    let (x_var, y_var) = vars;

    let mut arrows = Vec::new();

    let steps = (conf.samples as f64)
        .sqrt()
        .max(0.0)
        .round() as usize;

    for i in 0 .. steps {

        for j in 0 .. steps {

            let x = (x_range.1
                - x_range.0)
                .mul_add(
                    (i as f64)
                        / (steps - 1)
                            as f64,
                    x_range.0,
                );

            let y = (y_range.1
                - y_range.0)
                .mul_add(
                    (j as f64)
                        / (steps - 1)
                            as f64,
                    y_range.0,
                );

            let mut vars_map =
                HashMap::new();

            vars_map.insert(
                x_var.to_string(),
                x,
            );

            vars_map.insert(
                y_var.to_string(),
                y,
            );

            if let (Ok(vx), Ok(vy)) = (
                eval_expr(
                    vx_expr,
                    &vars_map,
                ),
                eval_expr(
                    vy_expr,
                    &vars_map,
                ),
            ) {

                let magnitude =
                    vx.hypot(vy);

                if magnitude > 1e-9 {

                    let scale = (x_range
                        .1
                        - x_range.0)
                        * 0.05;

                    let end_x = (vx
                        / magnitude)
                        .mul_add(
                            scale, x,
                        );

                    let end_y = (vy
                        / magnitude)
                        .mul_add(
                            scale, y,
                        );

                    arrows.push(PathElement::new(vec![(x, y), (end_x, end_y)], conf.line_color));
                }
            }
        }
    }

    chart
        .draw_series(arrows.into_iter())
        .map_err(|e| e.to_string())?;

    root.present()
        .map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D surface z = f(x, y) and saves it to a file.
///
/// # Errors
///
/// This function will return an error if the plot cannot be created or saved.

pub fn plot_surface_3d(
    expr: &Expr,
    vars: (&str, &str),
    x_range: (f64, f64),
    y_range: (f64, f64),
    path: &str,
    config: Option<PlotConfig>,
) -> Result<(), String> {

    let conf =
        config.unwrap_or_default();

    let root = BitMapBackend::new(
        path,
        (
            conf.width,
            conf.height,
        ),
    )
    .into_drawing_area();

    root.fill(&WHITE)
        .map_err(|e| e.to_string())?;

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                &conf.caption,
                ("sans-serif", 40)
                    .into_font(),
            )
            .build_cartesian_3d(
                x_range.0 .. x_range.1,
                -1.0 .. 1.0,
                y_range.0 .. y_range.1,
            )
            .map_err(|e| {

                e.to_string()
            })?;

    chart
        .configure_axes()
        .draw()
        .map_err(|e| e.to_string())?;

    let (x_var, y_var) = vars;

    let steps = (conf.samples as f64)
        .sqrt()
        .max(0.0)
        .round() as usize;

    let _ = chart.draw_series(
        SurfaceSeries::xoz(
            (0 .. steps).map(|i| {

                x_range.0
                    + (x_range.1
                        - x_range.0)
                        * (i as f64)
                        / (steps - 1)
                            as f64
            }),
            (0 .. steps).map(|i| {

                y_range.0
                    + (y_range.1
                        - y_range.0)
                        * (i as f64)
                        / (steps - 1)
                            as f64
            }),
            |x, z| {

                let mut vars_map =
                    HashMap::new();

                vars_map.insert(
                    x_var.to_string(),
                    x,
                );

                vars_map.insert(
                    y_var.to_string(),
                    z,
                );

                eval_expr(
                    expr,
                    &vars_map,
                )
                .unwrap_or(0.0)
            },
        )
        .style(
            conf.line_color
                .mix(0.5)
                .filled(),
        ),
    );

    root.present()
        .map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D parametric curve (x(t), y(t), z(t)) and saves it to a file.
///
/// # Errors
///
/// This function will return an error if the plot cannot be created or saved.

pub fn plot_parametric_curve_3d(
    comps: (&Expr, &Expr, &Expr),
    var: &str,
    range: (f64, f64),
    path: &str,
    config: Option<PlotConfig>,
) -> Result<(), String> {

    let conf =
        config.unwrap_or_default();

    let root = BitMapBackend::new(
        path,
        (
            conf.width,
            conf.height,
        ),
    )
    .into_drawing_area();

    root.fill(&WHITE)
        .map_err(|e| e.to_string())?;

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                &conf.caption,
                ("sans-serif", 40)
                    .into_font(),
            )
            .build_cartesian_3d(
                -3.0 .. 3.0,
                -3.0 .. 3.0,
                -3.0 .. 3.0,
            )
            .map_err(|e| {

                e.to_string()
            })?;

    chart
        .configure_axes()
        .draw()
        .map_err(|e| e.to_string())?;

    let (x_expr, y_expr, z_expr) =
        comps;

    chart
        .draw_series(LineSeries::new(
            (0..=conf.samples).map(|i| {
                let t = (range.1 - range.0).mul_add((i as f64) / conf.samples as f64, range.0);
                let mut vars_map = HashMap::new();
                vars_map.insert(var.to_string(), t);

                let x = eval_expr(x_expr, &vars_map).unwrap_or(0.0);
                let y = eval_expr(y_expr, &vars_map).unwrap_or(0.0);
                let z = eval_expr(z_expr, &vars_map).unwrap_or(0.0);
                (x, y, z)
            }),
            &conf.line_color,
        ))
        .map_err(|e| e.to_string())?;

    root.present()
        .map_err(|e| e.to_string())?;

    Ok(())
}

/// Plots a 3D vector field and saves it to a file.
///
/// # Errors
///
/// This function will return an error if the plot cannot be created or saved.

pub fn plot_vector_field_3d(
    comps: (&Expr, &Expr, &Expr),
    vars: (&str, &str, &str),
    ranges: (
        (f64, f64),
        (f64, f64),
        (f64, f64),
    ),
    path: &str,
    config: Option<PlotConfig>,
) -> Result<(), String> {

    let conf =
        config.unwrap_or_default();

    let root = BitMapBackend::new(
        path,
        (
            conf.width,
            conf.height,
        ),
    )
    .into_drawing_area();

    root.fill(&WHITE)
        .map_err(|e| e.to_string())?;

    let (x_range, y_range, z_range) =
        ranges;

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                &conf.caption,
                ("sans-serif", 40)
                    .into_font(),
            )
            .build_cartesian_3d(
                x_range.0 .. x_range.1,
                y_range.0 .. y_range.1,
                z_range.0 .. z_range.1,
            )
            .map_err(|e| {

                e.to_string()
            })?;

    chart
        .configure_axes()
        .draw()
        .map_err(|e| e.to_string())?;

    let (vx_expr, vy_expr, vz_expr) =
        comps;

    let (x_var, y_var, z_var) = vars;

    let mut arrows = Vec::new();

    let steps = (conf.samples as f64)
        .cbrt()
        .max(0.0)
        .round() as usize;

    for i in 0 .. steps {

        for j in 0 .. steps {

            for k in 0 .. steps {

                let x = (x_range.1
                    - x_range.0)
                    .mul_add(
                        (i as f64)
                            / (steps
                                - 1)
                                as f64,
                        x_range.0,
                    );

                let y = (y_range.1
                    - y_range.0)
                    .mul_add(
                        (j as f64)
                            / (steps
                                - 1)
                                as f64,
                        y_range.0,
                    );

                let z = (z_range.1
                    - z_range.0)
                    .mul_add(
                        (k as f64)
                            / (steps
                                - 1)
                                as f64,
                        z_range.0,
                    );

                let mut vars_map =
                    HashMap::new();

                vars_map.insert(
                    x_var.to_string(),
                    x,
                );

                vars_map.insert(
                    y_var.to_string(),
                    y,
                );

                vars_map.insert(
                    z_var.to_string(),
                    z,
                );

                if let (
                    Ok(vx),
                    Ok(vy),
                    Ok(vz),
                ) = (
                    eval_expr(
                        vx_expr,
                        &vars_map,
                    ),
                    eval_expr(
                        vy_expr,
                        &vars_map,
                    ),
                    eval_expr(
                        vz_expr,
                        &vars_map,
                    ),
                ) {

                    let magnitude = vz
                        .mul_add(
                            vz,
                            vx.mul_add(
                                vx,
                                vy * vy,
                            ),
                        )
                        .sqrt();

                    if magnitude > 1e-6
                    {

                        let scale = (x_range.1 - x_range.0) * 0.05;

                        let end_x = (vx / magnitude).mul_add(scale, x);

                        let end_y = (vy / magnitude).mul_add(scale, y);

                        let end_z = (vz / magnitude).mul_add(scale, z);

                        arrows.push(PathElement::new(vec![(x, y, z), (end_x, end_y, end_z)], conf.line_color));
                    }
                }
            }
        }
    }

    chart
        .draw_series(arrows.into_iter())
        .map_err(|e| e.to_string())?;

    root.present()
        .map_err(|e| e.to_string())?;

    Ok(())
}

#[cfg(test)]

mod tests {

    use std::collections::HashMap;

    use super::*;
    use crate::prelude::Expr;

    #[test]

    fn test_eval_basic() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::new_add(
            x,
            Expr::Constant(2.0),
        );

        let mut vars = HashMap::new();

        vars.insert(
            "x".to_string(),
            3.0,
        );

        let res =
            eval_expr(&expr, &vars)
                .unwrap();

        assert_eq!(res, 5.0);
    }

    #[test]

    fn test_eval_trig() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let expr = Expr::new_sin(x);

        let mut vars = HashMap::new();

        vars.insert(
            "x".to_string(),
            std::f64::consts::PI / 2.0,
        );

        let res =
            eval_expr(&expr, &vars)
                .unwrap();

        assert!(
            (res - 1.0).abs() < 1e-10
        );
    }

    #[test]

    fn test_eval_log_power() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        // e^(ln(x)) = x
        let expr = Expr::new_exp(
            Expr::new_log(x.clone()),
        );

        let mut vars = HashMap::new();

        vars.insert(
            "x".to_string(),
            5.0,
        );

        let res =
            eval_expr(&expr, &vars)
                .unwrap();

        assert!(
            (res - 5.0).abs() < 1e-10
        );
    }

    #[test]

    fn test_eval_error() {

        let x = Expr::Variable(
            "x".to_string(),
        );

        let vars = HashMap::new(); // Missing x
        let res = eval_expr(&x, &vars);

        assert!(res.is_err());
    }

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_no_panic_eval(
            val in -100.0..100.0f64,
            depth in 1..4usize,
        ) {
            let x = Expr::Variable("x".to_string());
            let mut expr = x.clone();
            for _ in 0..depth {
                expr = Expr::new_add(expr.clone(), Expr::Constant(1.0));
                expr = Expr::new_mul(expr.clone(), Expr::Constant(0.5));
            }
            let mut vars = HashMap::new();
            vars.insert("x".to_string(), val);
            let _ = eval_expr(&expr, &vars);
        }
    }
}
