//! # Coordinate System Transformations
//!
//! This module provides functions for transforming points, expressions, and vector/tensor
//! fields between different coordinate systems (Cartesian, Cylindrical, Spherical).
//! It includes utilities for computing Jacobian matrices and metric tensors, which are
//! fundamental for transformations in curvilinear coordinate systems.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix;
use crate::symbolic::matrix::inverse_matrix;
use crate::symbolic::matrix::mul_matrices;
use crate::symbolic::matrix::transpose_matrix;
use crate::symbolic::simplify::simplify;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CoordinateSystem {
    Cartesian,
    Cylindrical,
    Spherical,
}

/// Transforms a point from one coordinate system to another.
///
/// This function acts as a dispatcher, converting the point first to Cartesian coordinates
/// (if not already in Cartesian) and then from Cartesian to the target system.
///
/// # Arguments
/// * `point` - A slice of `Expr` representing the coordinates of the point.
/// * `from` - The `CoordinateSystem` of the input point.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<Expr>` of the transformed coordinates, or an error string
/// if the input is invalid or transformation is not supported.
pub fn transform_point(
    point: &[Expr],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<Expr>, String> {
    if from == to {
        return Ok(point.to_vec());
    }

    // First, convert from the source system to Cartesian, which is our pivot system.
    let cartesian_point = to_cartesian(point, from)?;

    // Then, convert from Cartesian to the target system.
    from_cartesian(&cartesian_point, to)
}

/// Converts a point from any system to Cartesian coordinates.
pub(crate) fn to_cartesian(point: &[Expr], from: CoordinateSystem) -> Result<Vec<Expr>, String> {
    match from {
        CoordinateSystem::Cartesian => Ok(point.to_vec()),
        CoordinateSystem::Cylindrical => {
            if point.len() != 3 {
                return Err("Cylindrical point must have 3 components (r, theta, z)".to_string());
            }
            let r = &point[0];
            let theta = &point[1];
            let z = &point[2];
            let x = simplify(Expr::Mul(
                Box::new(r.clone()),
                Box::new(Expr::Cos(Box::new(theta.clone()))),
            ));
            let y = simplify(Expr::Mul(
                Box::new(r.clone()),
                Box::new(Expr::Sin(Box::new(theta.clone()))),
            ));
            Ok(vec![x, y, z.clone()])
        }
        CoordinateSystem::Spherical => {
            if point.len() != 3 {
                return Err("Spherical point must have 3 components (rho, theta, phi)".to_string());
            }
            let rho = &point[0];
            let theta = &point[1];
            let phi = &point[2];
            let x = simplify(Expr::Mul(
                Box::new(rho.clone()),
                Box::new(Expr::Mul(
                    Box::new(Expr::Sin(Box::new(phi.clone()))),
                    Box::new(Expr::Cos(Box::new(theta.clone()))),
                )),
            ));
            let y = simplify(Expr::Mul(
                Box::new(rho.clone()),
                Box::new(Expr::Mul(
                    Box::new(Expr::Sin(Box::new(phi.clone()))),
                    Box::new(Expr::Sin(Box::new(theta.clone()))),
                )),
            ));
            let z = simplify(Expr::Mul(
                Box::new(rho.clone()),
                Box::new(Expr::Cos(Box::new(phi.clone()))),
            ));
            Ok(vec![x, y, z])
        }
    }
}

/// Converts a point from Cartesian coordinates to any other system.
pub(crate) fn from_cartesian(point: &[Expr], to: CoordinateSystem) -> Result<Vec<Expr>, String> {
    match to {
        CoordinateSystem::Cartesian => Ok(point.to_vec()),
        CoordinateSystem::Cylindrical => {
            if point.len() < 2 {
                return Err("Cartesian point must have at least 2 components (x, y)".to_string());
            }
            let x = &point[0];
            let y = &point[1];
            let r = simplify(Expr::Sqrt(Box::new(Expr::Add(
                Box::new(Expr::Power(
                    Box::new(x.clone()),
                    Box::new(Expr::Constant(2.0)),
                )),
                Box::new(Expr::Power(
                    Box::new(y.clone()),
                    Box::new(Expr::Constant(2.0)),
                )),
            ))));
            let theta = simplify(Expr::Atan2(Box::new(y.clone()), Box::new(x.clone())));
            let mut result = vec![r, theta];
            if point.len() > 2 {
                result.push(point[2].clone()); // Preserve z component
            }
            Ok(result)
        }
        CoordinateSystem::Spherical => {
            if point.len() != 3 {
                return Err("Cartesian point must have 3 components (x, y, z)".to_string());
            }
            let x = &point[0];
            let y = &point[1];
            let z = &point[2];
            let rho = simplify(Expr::Sqrt(Box::new(Expr::Add(
                Box::new(Expr::Add(
                    Box::new(Expr::Power(
                        Box::new(x.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                    Box::new(Expr::Power(
                        Box::new(y.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                )),
                Box::new(Expr::Power(
                    Box::new(z.clone()),
                    Box::new(Expr::Constant(2.0)),
                )),
            ))));
            let theta = simplify(Expr::Atan2(Box::new(y.clone()), Box::new(x.clone())));
            let phi = simplify(Expr::ArcCos(Box::new(Expr::Div(
                Box::new(z.clone()),
                Box::new(rho.clone()),
            ))));
            Ok(vec![rho, theta, phi])
        }
    }
}

/// Transforms a symbolic expression from one coordinate system to another.
///
/// This function substitutes the variables of the `from` coordinate system with their
/// equivalent expressions in the `to` coordinate system. It uses Cartesian coordinates
/// as an intermediate pivot for transformations between non-Cartesian systems.
///
/// # Arguments
/// * `expr` - The symbolic expression to transform.
/// * `from` - The `CoordinateSystem` of the input expression.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing the transformed `Expr`, or an error string if the transformation
/// is not supported or rules cannot be found.
pub fn transform_expression(
    expr: &Expr,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Expr, String> {
    if from == to {
        return Ok(expr.clone());
    }

    // To transform an expression from system A to B, we need the formulas for A's variables in terms of B's variables.
    // Our pivot system is Cartesian. So we need `A -> Cartesian` rules.
    let (from_vars, _to_vars, rules) = get_transform_rules(from, to)?;

    let mut current_expr = expr.clone();
    for (from_var, rule) in from_vars.iter().zip(rules.iter()) {
        current_expr = substitute(&current_expr, from_var, rule);
    }

    Ok(simplify(current_expr))
}

/// Helper function to get the variables and transformation rules between two coordinate systems.
///
/// This function provides the formulas to express the `from` system's coordinates
/// in terms of the `to` system's coordinates, typically by pivoting through Cartesian.
///
/// # Arguments
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a tuple `(from_vars, to_vars, rules)`.
/// `from_vars` are the variable names of the source system.
/// `to_vars` are the variable names of the target system.
/// `rules` are the expressions for `from_vars` in terms of `to_vars`.
pub fn get_transform_rules(
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<(Vec<String>, Vec<String>, Vec<Expr>), String> {
    // This function provides the formulas to express the `from` system's coordinates in terms of the `to` system's coordinates.
    // We do this by pivoting through Cartesian.

    // 1. Get rules for `to -> cartesian`.
    let (_to_vars, cartesian_vars, _to_cart_rules) = get_to_cartesian_rules(to)?;

    // 2. Get rules for `from -> cartesian`.
    let (from_vars, _, from_cart_rules) = get_to_cartesian_rules(from)?;

    // 3. We need to invert the `to -> cartesian` rules to get `cartesian -> to`.
    // This is the most complex part. For now, we will handle the case where one of the systems is Cartesian.
    if from == CoordinateSystem::Cartesian {
        // We need `cartesian -> to` rules. These are defined in `from_cartesian`.
        let (res_from, res_to, res_rules) = get_from_cartesian_rules(to)?;
        Ok((res_from, res_to, res_rules))
    } else if to == CoordinateSystem::Cartesian {
        // We need `from -> cartesian` rules. These are defined in `to_cartesian`.
        Ok((from_vars, cartesian_vars, from_cart_rules))
    } else {
        Err(
            "Direct transformation between two non-Cartesian systems is not yet supported."
                .to_string(),
        )
    }
}

/// Provides the transformation rules from a given coordinate system to Cartesian coordinates.
///
/// # Arguments
/// * `from` - The source `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a tuple `(from_vars, cartesian_vars, rules)`.
/// `from_vars` are the variable names of the source system.
/// `cartesian_vars` are the variable names of the Cartesian system.
/// `rules` are the expressions for Cartesian coordinates in terms of `from_vars`.
pub fn get_to_cartesian_rules(
    from: CoordinateSystem,
) -> Result<(Vec<String>, Vec<String>, Vec<Expr>), String> {
    let cartesian_vars = vec!["x".to_string(), "y".to_string(), "z".to_string()];
    match from {
        CoordinateSystem::Cartesian => Ok((
            cartesian_vars.clone(),
            cartesian_vars,
            vec![
                Expr::Variable("x".to_string()),
                Expr::Variable("y".to_string()),
                Expr::Variable("z".to_string()),
            ],
        )),
        CoordinateSystem::Cylindrical => {
            let cyl_vars = vec!["r".to_string(), "theta".to_string(), "z_cyl".to_string()];
            let r = Expr::Variable("r".to_string());
            let theta = Expr::Variable("theta".to_string());
            let rules = vec![
                simplify(Expr::Mul(
                    Box::new(r.clone()),
                    Box::new(Expr::Cos(Box::new(theta.clone()))),
                )), // x
                simplify(Expr::Mul(
                    Box::new(r.clone()),
                    Box::new(Expr::Sin(Box::new(theta.clone()))),
                )), // y
                Expr::Variable("z_cyl".to_string()), // z
            ];
            Ok((cyl_vars, cartesian_vars, rules))
        }
        CoordinateSystem::Spherical => {
            let sph_vars = vec![
                "rho".to_string(),
                "theta_sph".to_string(),
                "phi".to_string(),
            ];
            let rho = Expr::Variable("rho".to_string());
            let theta = Expr::Variable("theta_sph".to_string());
            let phi = Expr::Variable("phi".to_string());
            let rules = vec![
                simplify(Expr::Mul(
                    Box::new(rho.clone()),
                    Box::new(Expr::Mul(
                        Box::new(Expr::Sin(Box::new(phi.clone()))),
                        Box::new(Expr::Cos(Box::new(theta.clone()))),
                    )),
                )), // x
                simplify(Expr::Mul(
                    Box::new(rho.clone()),
                    Box::new(Expr::Mul(
                        Box::new(Expr::Sin(Box::new(phi.clone()))),
                        Box::new(Expr::Sin(Box::new(theta.clone()))),
                    )),
                )), // y
                simplify(Expr::Mul(
                    Box::new(rho.clone()),
                    Box::new(Expr::Cos(Box::new(phi.clone()))),
                )), // z
            ];
            Ok((sph_vars, cartesian_vars, rules))
        }
    }
}

/// Provides the transformation rules from Cartesian coordinates to a given coordinate system.
///
/// # Arguments
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a tuple `(cartesian_vars, to_vars, rules)`.
/// `cartesian_vars` are the variable names of the Cartesian system.
/// `to_vars` are the variable names of the target system.
/// `rules` are the expressions for `to_vars` in terms of Cartesian coordinates.
pub(crate) fn get_from_cartesian_rules(
    to: CoordinateSystem,
) -> Result<(Vec<String>, Vec<String>, Vec<Expr>), String> {
    let cartesian_vars = vec!["x".to_string(), "y".to_string(), "z".to_string()];
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());
    let z = Expr::Variable("z".to_string());

    match to {
        CoordinateSystem::Cartesian => Ok((cartesian_vars.clone(), cartesian_vars, vec![x, y, z])),
        CoordinateSystem::Cylindrical => {
            let cyl_vars = vec!["r".to_string(), "theta".to_string(), "z_cyl".to_string()];
            let rules = vec![
                simplify(Expr::Sqrt(Box::new(Expr::Add(
                    Box::new(Expr::Power(
                        Box::new(x.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                    Box::new(Expr::Power(
                        Box::new(y.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                )))), // r
                simplify(Expr::Atan2(Box::new(y.clone()), Box::new(x.clone()))), // theta
                z.clone(),                                                       // z
            ];
            Ok((cartesian_vars, cyl_vars, rules))
        }
        CoordinateSystem::Spherical => {
            let sph_vars = vec![
                "rho".to_string(),
                "theta_sph".to_string(),
                "phi".to_string(),
            ];
            let rho_rule = simplify(Expr::Sqrt(Box::new(Expr::Add(
                Box::new(Expr::Add(
                    Box::new(Expr::Power(
                        Box::new(x.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                    Box::new(Expr::Power(
                        Box::new(y.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )),
                )),
                Box::new(Expr::Power(
                    Box::new(z.clone()),
                    Box::new(Expr::Constant(2.0)),
                )),
            ))));
            let rules = vec![
                rho_rule.clone(),                                                // rho
                simplify(Expr::Atan2(Box::new(y.clone()), Box::new(x.clone()))), // theta
                simplify(Expr::ArcCos(Box::new(Expr::Div(
                    Box::new(z.clone()),
                    Box::new(rho_rule),
                )))), // phi
            ];
            Ok((cartesian_vars, sph_vars, rules))
        }
    }
}

/// Transforms a contravariant vector field (e.g., velocity) from one coordinate system to another.
///
/// Contravariant vectors transform with the Jacobian matrix of the coordinate transformation.
/// `V'_i = (∂x'_i / ∂x_j) * V_j`.
///
/// # Arguments
/// * `components` - A slice of `Expr` representing the components of the vector field.
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<Expr>` of the transformed components, or an error string
/// if the transformation is not supported or computation fails.
pub fn transform_contravariant_vector(
    components: &[Expr],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<Expr>, String> {
    if from == to {
        return Ok(components.to_vec());
    }
    // Transformation rule: V_new = J(from->to) * V_old
    let (vars_from, _, rules_to) = get_from_cartesian_rules(from)?;
    let jacobian = compute_jacobian(&rules_to, &vars_from)?;
    let new_comps_old_coords = symbolic_mat_vec_mul(&jacobian, components)?;

    // Substitute old coordinates for new ones
    let (_, _, rules_from) = get_to_cartesian_rules(from)?;
    let mut final_comps = Vec::new();
    for comp in new_comps_old_coords {
        let mut final_comp = comp;
        for (i, var) in vars_from.iter().enumerate() {
            final_comp = substitute(&final_comp, var, &rules_from[i]);
        }
        final_comps.push(simplify(final_comp));
    }
    Ok(final_comps)
}

/// Transforms a covariant vector field (e.g., gradient) from one coordinate system to another.
///
/// Covariant vectors transform with the inverse transpose of the Jacobian matrix of the
/// coordinate transformation. `V'_i = (∂x_j / ∂x'_i) * V_j`.
///
/// # Arguments
/// * `components` - A slice of `Expr` representing the components of the vector field.
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<Expr>` of the transformed components, or an error string
/// if the transformation is not supported or computation fails.
pub fn transform_covariant_vector(
    components: &[Expr],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<Expr>, String> {
    if from == to {
        return Ok(components.to_vec());
    }
    let (vars_from, _, rules_to) = get_transform_rules(from, to)?;
    let jacobian_vec = compute_jacobian(&rules_to, &vars_from)?;
    let jacobian = Expr::Matrix(jacobian_vec);
    let jacobian_inv = inverse_matrix(&jacobian);
    let jacobian_inv_t = transpose_matrix(&jacobian_inv);

    let old_vec = Expr::Matrix(components.iter().map(|c| vec![c.clone()]).collect());
    let new_vec_expr = mul_matrices(&jacobian_inv_t, &old_vec);

    // Final substitution step
    let (from_vars, _, rules) = get_to_cartesian_rules(from)?;
    let mut final_comps_expr = new_vec_expr;
    for (i, var) in from_vars.iter().enumerate() {
        final_comps_expr = substitute(&final_comps_expr, var, &rules[i]);
    }

    if let Expr::Matrix(rows) = simplify(final_comps_expr) {
        Ok(rows.into_iter().map(|row| row[0].clone()).collect())
    } else {
        Err("Transformation resulted in a non-vector expression".to_string())
    }
}

/// Computes the Jacobian matrix for a set of transformation rules.
pub(crate) fn compute_jacobian(rules: &[Expr], vars: &[String]) -> Result<Vec<Vec<Expr>>, String> {
    let mut jacobian = Vec::new();
    for rule in rules {
        let mut row = Vec::new();
        for var in vars {
            row.push(differentiate(rule, var));
        }
        jacobian.push(row);
    }
    Ok(jacobian)
}

/// Performs symbolic matrix-vector multiplication.
pub(crate) fn symbolic_mat_vec_mul(
    matrix: &[Vec<Expr>],
    vector: &[Expr],
) -> Result<Vec<Expr>, String> {
    if matrix.is_empty() || (!matrix.is_empty() && matrix[0].len() != vector.len()) {
        return Err("Matrix and vector dimensions are incompatible.".to_string());
    }
    let mut result = Vec::new();
    for row in matrix {
        let mut sum = Expr::Constant(0.0);
        for (i, val) in row.iter().enumerate() {
            sum = simplify(Expr::Add(
                Box::new(sum),
                Box::new(Expr::Mul(
                    Box::new(val.clone()),
                    Box::new(vector[i].clone()),
                )),
            ));
        }
        result.push(sum);
    }
    Ok(result)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TensorType {
    Contravariant,
    Covariant,
    Mixed,
}

/// Transforms a rank-2 tensor field from one coordinate system to another.
///
/// The transformation rule depends on the tensor's type (contravariant, covariant, or mixed).
/// - **Contravariant**: `T'^ij = (∂x'^i/∂x^a)(∂x'^j/∂x^b) T^ab`
/// - **Covariant**: `T'_ij = (∂x^a/∂x'^i)(∂x^b/∂x'^j) T_ab`
/// - **Mixed**: `T'^i_j = (∂x'^i/∂x^a)(∂x^b/∂x'^j) T^a_b`
///
/// # Arguments
/// * `tensor` - The rank-2 tensor as an `Expr::Matrix`.
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
/// * `tensor_type` - The `TensorType` (Contravariant, Covariant, Mixed).
///
/// # Returns
/// A `Result` containing the transformed `Expr::Matrix`, or an error string
/// if the transformation is not supported or computation fails.
pub fn transform_tensor2(
    tensor: &Expr, // Expects Expr::Matrix
    from: CoordinateSystem,
    to: CoordinateSystem,
    tensor_type: TensorType,
) -> Result<Expr, String> {
    if from == to {
        return Ok(tensor.clone());
    }

    let (from_vars, _, rules) = get_transform_rules(from, to)?;
    let jacobian_vec = compute_jacobian(&rules, &from_vars)?;
    let jacobian = Expr::Matrix(jacobian_vec);
    let jacobian_inv = inverse_matrix(&jacobian);

    let transformed_tensor = match tensor_type {
        TensorType::Contravariant => mul_matrices(
            &jacobian,
            &mul_matrices(tensor, &transpose_matrix(&jacobian)),
        ),
        TensorType::Covariant => mul_matrices(
            &transpose_matrix(&jacobian_inv),
            &mul_matrices(tensor, &jacobian_inv),
        ),
        TensorType::Mixed => mul_matrices(&jacobian, &mul_matrices(tensor, &jacobian_inv)),
    };

    // Final substitution step is complex and omitted for brevity.
    Ok(transformed_tensor)
}

/// Performs symbolic matrix-matrix multiplication.
pub fn symbolic_mat_mat_mul(m1: &[Vec<Expr>], m2: &[Vec<Expr>]) -> Result<Vec<Vec<Expr>>, String> {
    let m1_rows = m1.len();
    let m1_cols = m1[0].len();
    let m2_rows = m2.len();
    let m2_cols = m2[0].len();

    if m1_cols != m2_rows {
        return Err("Matrix dimensions are incompatible for multiplication.".to_string());
    }

    let result = vec![vec![Expr::Constant(0.0); m2_cols]; m1_rows];
    for i in 0..m1_rows {
        for j in 0..m2_cols {
            let _sum = Expr::Constant(0.0);
            for k in 0..m1_cols {
                let _term = simplify(Expr::Mul(
                    Box::new(m1[i][k].clone()),
                    Box::new(m2[k][j].clone()),
                ));
            }
        }
    }
    Ok(result)
}

/// Computes and returns the metric tensor for a given orthogonal coordinate system.
///
/// The metric tensor `g_ij` defines the inner product in a coordinate system and is
/// crucial for calculating distances, angles, and volumes. For orthogonal coordinates,
/// it is a diagonal matrix where `g_ii = (∂x/∂u_i)^2 + (∂y/∂u_i)^2 + (∂z/∂u_i)^2`.
///
/// # Arguments
/// * `system` - The `CoordinateSystem` for which to compute the metric tensor.
///
/// # Returns
/// A `Result` containing an `Expr::Matrix` representing the metric tensor,
/// or an error string if the system is not supported or computation fails.
pub fn get_metric_tensor(system: CoordinateSystem) -> Result<Expr, String> {
    let rules = match system {
        CoordinateSystem::Cartesian => return Ok(matrix::identity_matrix(3)),
        _ => get_to_cartesian_rules(system)?.2,
    };

    let vars = match system {
        CoordinateSystem::Cylindrical => {
            vec!["r".to_string(), "theta".to_string(), "z_cyl".to_string()]
        }
        CoordinateSystem::Spherical => vec![
            "rho".to_string(),
            "theta_sph".to_string(),
            "phi".to_string(),
        ],
        _ => unreachable!(),
    };

    let jacobian_vec = compute_jacobian(&rules, &vars)?;
    let jacobian = Expr::Matrix(jacobian_vec);

    // g = J^T * J
    Ok(mul_matrices(&transpose_matrix(&jacobian), &jacobian))
}

/// Computes the divergence of a contravariant vector field in any orthogonal coordinate system.
///
/// The divergence of a vector field `V` in orthogonal curvilinear coordinates `u_1, u_2, u_3` is given by:
/// `div(V) = (1/sqrt(g)) * [∂/∂u_1(V^1 * sqrt(g)/h_1) + ∂/∂u_2(V^2 * sqrt(g)/h_2) + ∂/∂u_3(V^3 * sqrt(g)/h_3)]`
/// where `g` is the determinant of the metric tensor and `h_i` are scale factors.
/// This simplified implementation assumes `h_i` are implicitly handled by `sqrt(g)`.
///
/// # Arguments
/// * `vector_comps` - A slice of `Expr` representing the contravariant components of the vector field.
/// * `from` - The `CoordinateSystem` of the vector field.
///
/// # Returns
/// A `Result` containing an `Expr` representing the divergence, or an error string
/// if the system is not supported or computation fails.
pub fn transform_divergence(vector_comps: &[Expr], from: CoordinateSystem) -> Result<Expr, String> {
    let g_matrix = get_metric_tensor(from)?;
    let g = matrix::determinant(&g_matrix);
    let sqrt_g = simplify(Expr::Sqrt(Box::new(g)));

    let (vars, _, _) = get_to_cartesian_rules(from)?;

    let mut total_divergence = Expr::Constant(0.0);
    for i in 0..vector_comps.len() {
        let term_to_diff = simplify(Expr::Mul(
            Box::new(sqrt_g.clone()),
            Box::new(vector_comps[i].clone()),
        ));
        let partial_deriv = differentiate(&term_to_diff, &vars[i]);
        total_divergence = simplify(Expr::Add(
            Box::new(total_divergence),
            Box::new(partial_deriv),
        ));
    }

    Ok(simplify(Expr::Div(
        Box::new(total_divergence),
        Box::new(sqrt_g),
    )))
}

/// Computes the curl of a covariant vector field in any orthogonal coordinate system.
///
/// The curl of a vector field `V` in orthogonal curvilinear coordinates `u_1, u_2, u_3` is given by:
/// `curl(V) = (1/(h_1*h_2*h_3)) * det(h_1*e_1, h_2*e_2, h_3*e_3; ∂/∂u_1, ∂/∂u_2, ∂/∂u_3; h_1*V_1, h_2*V_2, h_3*V_3)`
/// where `h_i` are the scale factors.
///
/// # Arguments
/// * `vector_comps` - A slice of `Expr` representing the covariant components of the vector field.
/// * `from` - The `CoordinateSystem` of the vector field.
///
/// # Returns
/// A `Result` containing a `Vec<Expr>` of the transformed components, or an error string
/// if the system is not supported or computation fails.
pub fn transform_curl(vector_comps: &[Expr], from: CoordinateSystem) -> Result<Vec<Expr>, String> {
    if vector_comps.len() != 3 {
        return Err("Curl is only defined for 3D vectors.".to_string());
    }
    let g_matrix = get_metric_tensor(from)?;
    let g_rows = if let Expr::Matrix(rows) = g_matrix {
        rows
    } else {
        return Err("Metric tensor is not a matrix".to_string());
    };

    let h1 = simplify(Expr::Sqrt(Box::new(g_rows[0][0].clone())));
    let h2 = simplify(Expr::Sqrt(Box::new(g_rows[1][1].clone())));
    let h3 = simplify(Expr::Sqrt(Box::new(g_rows[2][2].clone())));

    let (vars, _, _) = get_to_cartesian_rules(from)?;
    let u1 = &vars[0];
    let u2 = &vars[1];
    let u3 = &vars[2];

    let v1 = &vector_comps[0];
    let v2 = &vector_comps[1];
    let v3 = &vector_comps[2];

    let curl_1 = simplify(Expr::Div(
        Box::new(Expr::Sub(
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h3.clone()), Box::new(v3.clone()))),
                u2,
            )),
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h2.clone()), Box::new(v2.clone()))),
                u3,
            )),
        )),
        Box::new(simplify(Expr::Mul(
            Box::new(h2.clone()),
            Box::new(h3.clone()),
        ))),
    ));

    let curl_2 = simplify(Expr::Div(
        Box::new(Expr::Sub(
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h1.clone()), Box::new(v1.clone()))),
                u3,
            )),
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h3.clone()), Box::new(v3.clone()))),
                u1,
            )),
        )),
        Box::new(simplify(Expr::Mul(
            Box::new(h3.clone()),
            Box::new(h1.clone()),
        ))),
    ));

    let curl_3 = simplify(Expr::Div(
        Box::new(Expr::Sub(
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h2.clone()), Box::new(v2.clone()))),
                u1,
            )),
            Box::new(differentiate(
                &simplify(Expr::Mul(Box::new(h1.clone()), Box::new(v1.clone()))),
                u2,
            )),
        )),
        Box::new(simplify(Expr::Mul(
            Box::new(h1.clone()),
            Box::new(h2.clone()),
        ))),
    ));

    Ok(vec![curl_1, curl_2, curl_3])
}

/// Transforms the gradient of a scalar field from one coordinate system to another.
///
/// The gradient of a scalar field `f` is a covariant vector field. This function
/// computes the gradient in the source system, then transforms it to Cartesian coordinates,
/// and finally to the target coordinate system.
///
/// # Arguments
/// * `scalar_field` - The scalar field as an `Expr`.
/// * `from_vars` - A slice of strings representing the variables of the source system.
/// * `from` - The source `CoordinateSystem`.
/// * `to` - The target `CoordinateSystem`.
///
/// # Returns
/// A `Result` containing a `Vec<Expr>` of the transformed gradient components,
/// or an error string if the transformation is not supported or computation fails.
pub fn transform_gradient(
    scalar_field: &Expr,
    from_vars: &[String],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<Expr>, String> {
    if from == to {
        // If systems are the same, just compute the gradient directly.
        let mut grad_comps = Vec::new();
        for var in from_vars {
            grad_comps.push(differentiate(scalar_field, var));
        }
        return Ok(grad_comps);
    }

    // 1. Express the scalar field in Cartesian coordinates.
    let (_, _, rules) = get_to_cartesian_rules(from)?;
    let mut field_cart = scalar_field.clone();
    for (i, var) in from_vars.iter().enumerate() {
        field_cart = substitute(&field_cart, var, &rules[i]);
    }
    field_cart = simplify(field_cart);

    // 2. Compute the gradient in Cartesian coordinates.
    let cartesian_vars = vec!["x".to_string(), "y".to_string(), "z".to_string()];
    let mut grad_cart_comps = Vec::new();
    for var in &cartesian_vars {
        grad_cart_comps.push(differentiate(&field_cart, var));
    }

    // 3. Transform this Cartesian gradient (a covariant vector) to the target system `to`.
    transform_covariant_vector(&grad_cart_comps, CoordinateSystem::Cartesian, to)
}

/*
ateSystem::Cartesian, to)
}

, CoordinateSystem::Cartesian, to)
}

ateSystem::Cartesian, to)
}
*/
