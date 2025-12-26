//! # Numerical Differential Geometry
//!
//! This module provides numerical tools for differential geometry.
//! It includes functions for computing the metric tensor, Christoffel symbols,
//! Riemann curvature tensor, Ricci tensor, and Ricci scalar for various coordinate systems.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::coordinates::{
    self,
    CoordinateSystem,
};
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::{
    self as symbolic_matrix,
    inverse_matrix,
};
use num_traits::ToPrimitive;
use std::collections::HashMap;

/// Evaluates the metric tensor at a given point for a coordinate system.

pub fn metric_tensor_at_point(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {

    let g_sym = coordinates::get_metric_tensor(system)?;

    let (vars, _, _) = coordinates::get_to_cartesian_rules(system)?;

    let mut eval_map = HashMap::new();

    for (i, var) in vars
        .iter()
        .enumerate()
    {

        eval_map.insert(
            var.clone(),
            point[i],
        );
    }

    if let Expr::Matrix(rows) = g_sym {

        let mut g_num = vec![vec![0.0; rows[0].len()]; rows.len()];

        for i in 0..rows.len() {

            for j in 0..rows[i].len() {

                g_num[i][j] = eval_expr(
                    &rows[i][j],
                    &eval_map,
                )?;
            }
        }

        Ok(g_num)
    } else {

        Err("Metric tensor is not a matrix".to_string())
    }
}

/// Computes the Christoffel symbols of the second kind at a given point.
///
/// The Christoffel symbols `Γ^k_{ij}` describe the connection coefficients of a metric tensor.
/// They are used to define covariant derivatives and curvature. The formula is:
/// `Γ^k_{ij} = (1/2) * g^{km} * (∂g_{mi}/∂u^j + ∂g_{mj}/∂u^i - ∂g_{ij}/∂u^m)`.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::differential_geometry::christoffel_symbols;
/// use rssn::symbolic::coordinates::CoordinateSystem;
///
/// // Spherical coordinates: (rho, theta, phi)
/// // At rho=1, theta=pi/2, phi=pi/2
/// let point = vec![
///     1.0,
///     1.57079632679,
///     1.57079632679,
/// ];
///
/// let symbols = christoffel_symbols(
///     CoordinateSystem::Spherical,
///     &point,
/// )
/// .unwrap();
/// ```

pub fn christoffel_symbols(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<Vec<Vec<Vec<f64>>>, String> {

    let (vars, _, _) = coordinates::get_to_cartesian_rules(system)?;

    let dim = vars.len();

    let mut eval_map = HashMap::new();

    for (i, var) in vars
        .iter()
        .enumerate()
    {

        eval_map.insert(
            var.clone(),
            point[i],
        );
    }

    let g_sym = coordinates::get_metric_tensor(system)?;

    let g_num = if let Expr::Matrix(rows) = &g_sym {

        let mut mat = vec![vec![0.0; dim]; dim];

        for i in 0..dim {

            for j in 0..dim {

                mat[i][j] = eval_expr(
                    &rows[i][j],
                    &eval_map,
                )?;
            }
        }

        mat
    } else {

        return Err("Metric tensor is not a matrix".into());
    };

    // Numerical inverse of metric tensor
    let g_mat_expr = Expr::Matrix(
        g_num
            .iter()
            .map(|r| {

                r.iter()
                    .map(|&v| Expr::Constant(v))
                    .collect()
            })
            .collect(),
    );

    let g_inv_sym = inverse_matrix(&g_mat_expr);

    let g_inv_num = if let Expr::Matrix(rows) = g_inv_sym {

        let mut mat = vec![vec![0.0; dim]; dim];

        for i in 0..dim {

            for j in 0..dim {

                if let Expr::Constant(v) = rows[i][j] {

                    mat[i][j] = v;
                } else if let Expr::BigInt(b) = &rows[i][j] {

                    mat[i][j] = b
                        .to_f64()
                        .unwrap_or(0.0);
                } else {

                    mat[i][j] = eval_expr(
                        &rows[i][j],
                        &HashMap::new(),
                    )?;
                }
            }
        }

        mat
    } else {

        return Err("Inverse metric tensor is not a matrix".into());
    };

    let g_sym_rows = if let Expr::Matrix(rows) = g_sym {

        rows
    } else {

        unreachable!()
    };

    let mut dg_num = vec![vec![vec![0.0; dim]; dim]; dim]; // [k][i][j] = ∂_k g_ij
    for k in 0..dim {

        for i in 0..dim {

            for j in 0..dim {

                let deriv = differentiate(
                    &g_sym_rows[i][j],
                    &vars[k],
                );

                dg_num[k][i][j] = eval_expr(&deriv, &eval_map)?;
            }
        }
    }

    let mut christoffel = vec![vec![vec![0.0; dim]; dim]; dim];

    for k in 0..dim {

        for i in 0..dim {

            for j in 0..dim {

                let mut sum = 0.0;

                for m in 0..dim {

                    // Γ^k_{ij} = 0.5 * g^{km} * (∂_j g_{mi} + ∂_i g_{mj} - ∂_m g_{ij})
                    sum += g_inv_num[k][m] * (dg_num[j][m][i] + dg_num[i][m][j] - dg_num[m][i][j]);
                }

                christoffel[k][i][j] = 0.5 * sum;
            }
        }
    }

    Ok(christoffel)
}

/// Computes the Riemann curvature tensor at a given point.
///
/// `R^ρ_{σμν} = ∂_μ Γ^ρ_{σν} - ∂_ν Γ^ρ_{σμ} + Γ^ρ_{λμ} Γ^λ_{σν} - Γ^ρ_{λν} Γ^λ_{σμ}`

pub fn riemann_tensor(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<Vec<Vec<Vec<Vec<f64>>>>, String> {

    let (vars, _, _) = coordinates::get_to_cartesian_rules(system)?;

    let dim = vars.len();

    let mut eval_map = HashMap::new();

    for (i, var) in vars
        .iter()
        .enumerate()
    {

        eval_map.insert(
            var.clone(),
            point[i],
        );
    }

    let g_sym = coordinates::get_metric_tensor(system)?;

    let g_sym_rows = if let Expr::Matrix(rows) = g_sym {

        rows
    } else {

        return Err("Invalid metric".into());
    };

    // Evaluate g_ij, ∂_k g_ij, ∂_l ∂_k g_ij at the point
    let mut g_num = vec![vec![0.0; dim]; dim];

    let mut dg_num = vec![vec![vec![0.0; dim]; dim]; dim]; // [k][i][j]
    let mut ddg_num = vec![vec![vec![vec![0.0; dim]; dim]; dim]; dim]; // [l][k][i][j]

    for i in 0..dim {

        for j in 0..dim {

            g_num[i][j] = eval_expr(
                &g_sym_rows[i][j],
                &eval_map,
            )?;

            for k in 0..dim {

                let dk_gij = differentiate(
                    &g_sym_rows[i][j],
                    &vars[k],
                );

                dg_num[k][i][j] = eval_expr(&dk_gij, &eval_map)?;

                for l in 0..dim {

                    let dl_dk_gij = differentiate(&dk_gij, &vars[l]);

                    ddg_num[l][k][i][j] = eval_expr(
                        &dl_dk_gij, &eval_map,
                    )?;
                }
            }
        }
    }

    // Numerical inverse g^ij
    let g_inv_num = invert_mat_num(&g_num)?;

    // Compute ∂_μ g^ρλ = -g^ρa (∂_μ g_ab) g^bλ
    let mut d_ginv_num = vec![vec![vec![0.0; dim]; dim]; dim]; // [μ][ρ][λ]
    for mu in 0..dim {

        for rho in 0..dim {

            for lambda in 0..dim {

                let mut sum = 0.0;

                for a in 0..dim {

                    for b in 0..dim {

                        sum += g_inv_num[rho][a] * dg_num[mu][a][b] * g_inv_num[b][lambda];
                    }
                }

                d_ginv_num[mu][rho][lambda] = -sum;
            }
        }
    }

    // Now compute Γ^ρ_{σν} and its derivative ∂_μ Γ^ρ_{σν}
    // Γ^ρ_{σν} = 0.5 * g^ρλ * (∂_ν g_{λσ} + ∂_σ g_{λν} - ∂_λ g_{σν})
    let mut gamma = vec![vec![vec![0.0; dim]; dim]; dim];

    for rho in 0..dim {

        for sigma in 0..dim {

            for nu in 0..dim {

                let mut sum = 0.0;

                for lambda in 0..dim {

                    sum += g_inv_num[rho][lambda]
                        * (dg_num[nu][lambda][sigma] + dg_num[sigma][lambda][nu]
                            - dg_num[lambda][sigma][nu]);
                }

                gamma[rho][sigma][nu] = 0.5 * sum;
            }
        }
    }

    let mut d_gamma = vec![vec![vec![vec![0.0; dim]; dim]; dim]; dim]; // [μ][ρ][σ][ν]
    for mu in 0..dim {

        for rho in 0..dim {

            for sigma in 0..dim {

                for nu in 0..dim {

                    let mut term1 = 0.0;

                    for lambda in 0..dim {

                        term1 += d_ginv_num[mu][rho][lambda]
                            * (dg_num[nu][lambda][sigma] + dg_num[sigma][lambda][nu]
                                - dg_num[lambda][sigma][nu]);
                    }

                    let mut term2 = 0.0;

                    for lambda in 0..dim {

                        term2 += g_inv_num[rho][lambda]
                            * (ddg_num[mu][nu][lambda][sigma] + ddg_num[mu][sigma][lambda][nu]
                                - ddg_num[mu][lambda][sigma][nu]);
                    }

                    d_gamma[mu][rho][sigma][nu] = 0.5 * (term1 + term2);
                }
            }
        }
    }

    let mut riemann = vec![vec![vec![vec![0.0; dim]; dim]; dim]; dim];

    for rho in 0..dim {

        for sigma in 0..dim {

            for mu in 0..dim {

                for nu in 0..dim {

                    let mut sum_product = 0.0;

                    for lambda in 0..dim {

                        sum_product += gamma[rho][lambda][mu] * gamma[lambda][sigma][nu]
                            - gamma[rho][lambda][nu] * gamma[lambda][sigma][mu];
                    }

                    riemann[rho][sigma][mu][nu] =
                        d_gamma[mu][rho][sigma][nu] - d_gamma[nu][rho][sigma][mu] + sum_product;
                }
            }
        }
    }

    Ok(riemann)
}

fn invert_mat_num(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {

    let dim = mat.len();

    let mat_expr = Expr::Matrix(
        mat.iter()
            .map(|r| {

                r.iter()
                    .map(|&v| Expr::Constant(v))
                    .collect()
            })
            .collect(),
    );

    let inv_expr = inverse_matrix(&mat_expr);

    if let Expr::Matrix(rows) = inv_expr {

        let mut res = vec![vec![0.0; dim]; dim];

        for i in 0..dim {

            for j in 0..dim {

                match &rows[i][j] {
                    | Expr::Constant(v) => res[i][j] = *v,
                    | Expr::BigInt(b) => {
                        res[i][j] = b
                            .to_f64()
                            .unwrap_or(0.0)
                    },
                    | _ => {
                        res[i][j] = eval_expr(
                            &rows[i][j],
                            &HashMap::new(),
                        )?
                    },
                }
            }
        }

        Ok(res)
    } else {

        Err("Could not invert matrix numerically".into())
    }
}

/// Computes the Ricci tensor at a given point.
///
/// `R_{σν} = R^μ_{σμν}` (Contraction of Riemann tensor)

pub fn ricci_tensor(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {

    let riemann = riemann_tensor(system, point)?;

    let dim = riemann.len();

    let mut ricci = vec![vec![0.0; dim]; dim];

    for sigma in 0..dim {

        for nu in 0..dim {

            let mut sum = 0.0;

            for mu in 0..dim {

                sum += riemann[mu][sigma][mu][nu];
            }

            ricci[sigma][nu] = sum;
        }
    }

    Ok(ricci)
}

/// Computes the Ricci scalar at a given point.
///
/// `R = g^{μν} R_{μν}`

pub fn ricci_scalar(
    system: CoordinateSystem,
    point: &[f64],
) -> Result<f64, String> {

    let ricci = ricci_tensor(system, point)?;

    let g_num = metric_tensor_at_point(system, point)?;

    let dim = g_num.len();

    // Compute numerical inverse of metric tensor
    let g_mat = Expr::Matrix(
        g_num
            .iter()
            .map(|r| {

                r.iter()
                    .map(|&v| Expr::Constant(v))
                    .collect()
            })
            .collect(),
    );

    let g_inv_sym = inverse_matrix(&g_mat);

    let mut r_scalar = 0.0;

    if let Expr::Matrix(rows) = g_inv_sym {

        for mu in 0..dim {

            for nu in 0..dim {

                let g_inv_val = if let Expr::Constant(v) = rows[mu][nu] {

                    v
                } else {

                    0.0
                };

                r_scalar += g_inv_val * ricci[mu][nu];
            }
        }
    }

    Ok(r_scalar)
}
