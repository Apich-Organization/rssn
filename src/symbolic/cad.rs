//! # Cylindrical Algebraic Decomposition (CAD)
//!
//! This module implements the Cylindrical Algebraic Decomposition algorithm for a set of
//! multivariate polynomials. CAD decomposes R^n into a finite number of cells, each
//! being a connected manifold, such that each polynomial has a constant sign on each cell.
//!
//! The algorithm consists of two phases:
//! 1. **Projection Phase**: Systematically reduces the number of variables by projecting the
//!    polynomials into lower-dimensional spaces using resultants and discriminants.
//! 2. **Lifting Phase**: Constructs cells in R^n by recursively lifting cells from R^(n-1)
//!    using sample points and root isolation.

use crate::symbolic::core::{Expr, SparsePolynomial};
use crate::symbolic::matrix;
use crate::symbolic::polynomial::{differentiate_poly, expr_to_sparse_poly, sparse_poly_to_expr};
use crate::symbolic::real_roots::isolate_real_roots;
use crate::symbolic::simplify::is_zero;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;

/// Represents a cell in the Cylindrical Algebraic Decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CadCell {
    /// A sample point that lies within the cell.
    pub sample_point: Vec<f64>,
    /// The dimension of the cell (e.g., 0 for a point, 1 for a curve, n for a region).
    pub dim: usize,
    /// An index representing the cell's position in the stack over a lower-dimensional cell.
    pub index: Vec<usize>,
}

/// Represents the full Cylindrical Algebraic Decomposition of R^n.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cad {
    /// The collection of cells in the decomposition.
    pub cells: Vec<CadCell>,
    /// The dimension of the space (n).
    pub dim: usize,
}

/// Computes the Cylindrical Algebraic Decomposition for a set of polynomials.
///
/// # Arguments
/// * `polys` - A slice of `SparsePolynomial` representing the input set.
/// * `vars` - The names of the variables in order (e.g., ["x", "y", "z"]).
///
/// # Returns
/// A `Result` containing the `Cad` structure or an error message.
pub fn cad(polys: &[SparsePolynomial], vars: &[&str]) -> Result<Cad, String> {
    if vars.is_empty() {
        return Err("Variable list cannot be empty.".to_string());
    }
    let projections = projection_phase(polys, vars)?;
    let cells = lifting_phase(&projections, vars)?;
    Ok(Cad {
        cells,
        dim: vars.len(),
    })
}

/// Performs the projection phase of CAD.
pub(crate) fn projection_phase(
    polys: &[SparsePolynomial],
    vars: &[&str],
) -> Result<Vec<Vec<SparsePolynomial>>, String> {
    let mut projection_sets = vec![polys.to_vec()];
    let mut current_polys = polys.to_vec();
    let mut current_vars = vars.to_vec();

    while current_vars.len() > 1 {
        let proj_var = current_vars.last().unwrap();
        let mut next_set = HashSet::new();

        // Discriminants (resultant(p, p_prime))
        for p in &current_polys {
            let p_prime = differentiate_poly(p, proj_var);
            if !p_prime.terms.is_empty() {
                let res = resultant(p, &p_prime, proj_var);
                println!("Resultant(p, p_prime, {proj_var}): {res:?}");
                if !is_zero(&res) {
                    let next_vars = &current_vars[0..current_vars.len() - 1];
                    next_set.insert(expr_to_sparse_poly(&res, next_vars));
                }
            }
        }

        // Cross-resultants
        for i in 0..current_polys.len() {
            for j in (i + 1)..current_polys.len() {
                let res = resultant(&current_polys[i], &current_polys[j], proj_var);
                if !is_zero(&res) {
                    let next_vars = &current_vars[0..current_vars.len() - 1];
                    next_set.insert(expr_to_sparse_poly(&res, next_vars));
                }
            }
        }

        current_vars.pop();
        current_polys = next_set.into_iter().collect();
        projection_sets.push(current_polys.clone());
    }

    projection_sets.reverse();
    Ok(projection_sets)
}

/// Performs the lifting phase of CAD.
pub(crate) fn lifting_phase(
    projections: &[Vec<SparsePolynomial>],
    vars: &[&str],
) -> Result<Vec<CadCell>, String> {
    let base_polys = &projections[0];
    let mut all_roots = Vec::new();
    for p in base_polys {
        let roots = isolate_real_roots(p, vars[0], 1e-9)?;
        for (a, b) in roots {
            all_roots.push(f64::midpoint(a, b));
        }
    }
    all_roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    all_roots.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

    let mut current_cells = Vec::new();

    // Construct first level cells
    if all_roots.is_empty() {
        current_cells.push(CadCell {
            sample_point: vec![0.0],
            dim: 1,
            index: vec![0],
        });
    } else {
        // Interval (-inf, first_root)
        current_cells.push(CadCell {
            sample_point: vec![all_roots[0] - 1.0],
            dim: 1,
            index: vec![0],
        });
        for (i, root) in all_roots.iter().enumerate() {
            // Point {root}
            current_cells.push(CadCell {
                sample_point: vec![*root],
                dim: 0,
                index: vec![2 * i + 1],
            });
            // Interval (root, next_root) or (root, inf)
            if i + 1 < all_roots.len() {
                current_cells.push(CadCell {
                    sample_point: vec![f64::midpoint(*root, all_roots[i + 1])],
                    dim: 1,
                    index: vec![2 * i + 2],
                });
            } else {
                current_cells.push(CadCell {
                    sample_point: vec![*root + 1.0],
                    dim: 1,
                    index: vec![2 * i + 2],
                });
            }
        }
    }

    for k in 1..vars.len() {
        let polys_k = &projections[k];
        let mut next_level_cells = Vec::new();

        for cell in &current_cells {
            let mut sample_map = HashMap::new();
            for (i, v) in vars.iter().enumerate().take(k) {
                sample_map.insert((*v).to_string(), cell.sample_point[i]);
            }

            let mut roots_at_sample = Vec::new();
            for p in polys_k {
                let p_expr = sparse_poly_to_expr(p);
                let p_substituted_expr = substitute_map(&p_expr, &sample_map);
                let p_substituted = expr_to_sparse_poly(&p_substituted_expr, &[vars[k]]);

                if p_substituted.terms.is_empty() {
                    continue;
                }
                let is_constant = p_substituted.degree(vars[k]) == 0;
                if is_constant {
                    continue;
                }

                let roots = isolate_real_roots(&p_substituted, vars[k], 1e-9)?;
                for (a, b) in roots {
                    roots_at_sample.push(f64::midpoint(a, b));
                }
            }
            roots_at_sample.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            roots_at_sample.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

            if roots_at_sample.is_empty() {
                let mut new_sample = cell.sample_point.clone();
                new_sample.push(0.0);
                let mut new_index = cell.index.clone();
                new_index.push(0);
                next_level_cells.push(CadCell {
                    sample_point: new_sample,
                    dim: cell.dim + 1,
                    index: new_index,
                });
            } else {
                // Interval (-inf, first_root)
                let mut new_sample = cell.sample_point.clone();
                new_sample.push(roots_at_sample[0] - 1.0);
                let mut new_index = cell.index.clone();
                new_index.push(0);
                next_level_cells.push(CadCell {
                    sample_point: new_sample,
                    dim: cell.dim + 1,
                    index: new_index,
                });

                for (i, root_val) in roots_at_sample.iter().enumerate() {
                    // Point {root}
                    let mut point_sample = cell.sample_point.clone();
                    point_sample.push(*root_val);
                    let mut point_index = cell.index.clone();
                    point_index.push(2 * i + 1);
                    next_level_cells.push(CadCell {
                        sample_point: point_sample,
                        dim: cell.dim,
                        index: point_index,
                    });

                    // Interval (root, next_root) or (root, inf)
                    let mut interval_sample = cell.sample_point.clone();
                    if i + 1 < roots_at_sample.len() {
                        interval_sample.push(f64::midpoint(*root_val, roots_at_sample[i + 1]));
                    } else {
                        interval_sample.push(*root_val + 1.0);
                    }
                    let mut interval_index = cell.index.clone();
                    interval_index.push(2 * i + 2);
                    next_level_cells.push(CadCell {
                        sample_point: interval_sample,
                        dim: cell.dim + 1,
                        index: interval_index,
                    });
                }
            }
        }
        current_cells = next_level_cells;
    }

    Ok(current_cells)
}

/// Evaluates an expression by substituting variables with constant values.
pub(crate) fn substitute_map(expr: &Expr, vars: &HashMap<String, f64>) -> Expr {
    let mut result = expr.clone();
    for (var, val) in vars {
        result = crate::symbolic::calculus::substitute(&result, var, &Expr::Constant(*val));
    }
    crate::symbolic::simplify_dag::simplify(&result)
}

/// Computes the Sylvester matrix of two polynomials with respect to a given variable.
#[allow(clippy::needless_range_loop)]
pub(crate) fn sylvester_matrix(p: &SparsePolynomial, q: &SparsePolynomial, var: &str) -> Expr {
    let n = p.degree(var) as usize;
    let m = q.degree(var) as usize;
    if n == 0 && m == 0 {
        return Expr::Matrix(vec![vec![Expr::Constant(0.0)]]);
    }
    let mut matrix_rows = vec![vec![Expr::Constant(0.0); n + m]; n + m];
    let p_coeffs_rev = p.get_coeffs_as_vec(var);
    let q_coeffs_rev = q.get_coeffs_as_vec(var);

    for i in 0..m {
        for j in 0..=n {
            if i + j < n + m {
                matrix_rows[i][i + j] = p_coeffs_rev
                    .get(j)
                    .cloned()
                    .unwrap_or_else(|| Expr::Constant(0.0));
            }
        }
    }
    for i in 0..n {
        for j in 0..=m {
            if i + j < n + m {
                matrix_rows[i + m][i + j] = q_coeffs_rev
                    .get(j)
                    .cloned()
                    .unwrap_or_else(|| Expr::Constant(0.0));
            }
        }
    }
    Expr::Matrix(matrix_rows)
}

/// Computes the resultant of two polynomials with respect to a given variable.
pub(crate) fn resultant(p: &SparsePolynomial, q: &SparsePolynomial, var: &str) -> Expr {
    let sylvester = sylvester_matrix(p, q, var);
    if let Expr::Matrix(m) = &sylvester {
        if m.is_empty() {
            return Expr::Constant(0.0);
        }
        if m.len() == 1 && m[0].len() == 1 {
            return m[0][0].clone();
        }
    }
    matrix::determinant(&sylvester)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbolic::core::Monomial;
    use std::collections::BTreeMap;

    #[test]
    fn test_simple_cad() {
        // p(x) = x^2 - 1
        let mut terms = BTreeMap::new();
        let mut vars_map = BTreeMap::new();
        vars_map.insert("x".to_string(), 2);
        terms.insert(Monomial(vars_map), Expr::Constant(1.0));
        let vars_map_0 = BTreeMap::new();
        terms.insert(Monomial(vars_map_0), Expr::Constant(-1.0));
        let p = SparsePolynomial { terms };

        let result = cad(&[p], &["x"]).unwrap();
        // Roots are -1 and 1.
        // Intervals: (-inf, -1), {-1}, (-1, 1), {1}, (1, inf)
        // Total 5 cells.
        assert_eq!(result.cells.len(), 5);
    }
}
