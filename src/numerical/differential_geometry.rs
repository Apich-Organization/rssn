//! # Numerical Differential Geometry
//!
//! This module provides numerical tools for differential geometry.
//! It focuses on computing Christoffel symbols at a given point for various
//! coordinate systems, which are fundamental for understanding curvature and
//! geodesics in Riemannian geometry.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::coordinates::{self, CoordinateSystem};
use crate::symbolic::core::Expr;
use crate::symbolic::matrix as symbolic_matrix;
use std::collections::HashMap;

/// Computes the Christoffel symbols of the second kind at a given point.
///
/// The Christoffel symbols `Γ^k_{ij}` describe the connection coefficients of a metric tensor.
/// They are used to define covariant derivatives and curvature. The formula is:
/// `Γ^k_{ij} = (1/2) * g^{km} * (∂g_{mi}/∂u^j + ∂g_{mj}/∂u^i - ∂g_{ij}/∂u^m)`.
/// This function evaluates the symbolic Christoffel symbols numerically at a specific point.
///
/// # Arguments
/// * `system` - The coordinate system to use.
/// * `point` - The point at which to evaluate the symbols.
///
/// # Returns
/// A `Result` containing a 3D vector representing the Christoffel symbols `Γ^k_{ij}`,
/// or an error string if the metric tensor is invalid or evaluation fails.
pub fn christoffel_symbols(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<Vec<Vec<Vec<f64>>>, String> {
    // 1. Get symbolic metric tensor and its inverse
    let g_sym = coordinates::get_metric_tensor(system)?;
    let g_inv_sym = symbolic_matrix::inverse_matrix(&g_sym);

    let (vars, _, _) = coordinates::get_to_cartesian_rules(system)?;
    let dim = vars.len();

    // 2. Create a map for the evaluation point
    let mut eval_map = HashMap::new();
    for (i, var) in vars.iter().enumerate() {
        eval_map.insert(var.clone(), point[i]);
    }

    // 3. Evaluate g_inv at the point
    let g_inv_num = if let Expr::Matrix(rows) = g_inv_sym {
        rows.iter()
            .map(|row| {
                row.iter()
                    .map(|elem| eval_expr(elem, &eval_map).unwrap_or(0.0))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    } else {
        return Err("Inverse metric tensor is not a matrix".to_string());
    };

    // 4. Symbolically differentiate the metric tensor components
    let g_sym_rows = if let Expr::Matrix(rows) = g_sym {
        rows
    } else {
        return Err("Metric tensor is not a matrix".to_string());
    };
    let mut g_derivs = vec![vec![vec![Expr::Constant(0.0); dim]; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                // k is the differentiation variable u^k
                g_derivs[i][j][k] = differentiate(&g_sym_rows[i][j], &vars[k]);
            }
        }
    }

    // 5. Evaluate derivatives at the point and compute Christoffel symbols
    let mut christoffel = vec![vec![vec![0.0; dim]; dim]; dim]; // Gamma^k_{i,j}
    for k in 0..dim {
        for i in 0..dim {
            for j in 0..dim {
                let mut sum = 0.0;
                for m in 0..dim {
                    let dg_mi_dj = eval_expr(&g_derivs[m][i][j], &eval_map)?;
                    let dg_mj_di = eval_expr(&g_derivs[m][j][i], &eval_map)?;
                    let dg_ij_dm = eval_expr(&g_derivs[i][j][m], &eval_map)?;
                    sum += g_inv_num[k][m] * (dg_mi_dj + dg_mj_di - dg_ij_dm);
                }
                christoffel[k][i][j] = 0.5 * sum;
            }
        }
    }

    Ok(christoffel)
}
