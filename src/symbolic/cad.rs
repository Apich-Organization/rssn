use crate::symbolic::core::{Expr, SparsePolynomial};
use crate::symbolic::matrix;
use crate::symbolic::polynomial::{differentiate_poly, expr_to_sparse_poly, sparse_poly_to_expr};
use crate::symbolic::real_roots::isolate_real_roots;
//use crate::symbolic::simplify::{is_zero, as_f64};
use crate::symbolic::simplify::is_zero;
use std::collections::HashMap;
use std::collections::HashSet;

/// Represents a cell in the Cylindrical Algebraic Decomposition.
#[derive(Debug, Clone)]
pub struct CadCell {
    /// A sample point that lies within the cell.
    pub sample_point: Vec<f64>,
    /// The dimension of the cell (e.g., 0 for a point, 1 for a curve, n for a region).
    pub dim: usize,
    /// An index representing the cell's position in the stack over a lower-dimensional cell.
    pub index: Vec<usize>,
}

/// Represents the full Cylindrical Algebraic Decomposition of R^n.
#[derive(Debug, Clone)]
pub struct Cad {
    pub cells: Vec<CadCell>,
    pub dim: usize,
}

/// Computes the Cylindrical Algebraic Decomposition for a set of polynomials.
pub fn cad(polys: &[SparsePolynomial], vars: &[&str]) -> Result<Cad, String> {
    // 1. Projection Phase
    let projections = projection_phase(polys, vars)?;

    // 2. Lifting Phase
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

        for p in &current_polys {
            let p_prime = differentiate_poly(p, proj_var);
            let res = resultant(p, &p_prime, proj_var);
            if !is_zero(&res) {
                next_set.insert(expr_to_sparse_poly(
                    &res,
                    &[&current_vars[0..current_vars.len() - 1].join(",")],
                ));
            }
        }

        for i in 0..current_polys.len() {
            for j in (i + 1)..current_polys.len() {
                let res = resultant(&current_polys[i], &current_polys[j], proj_var);
                if !is_zero(&res) {
                    next_set.insert(expr_to_sparse_poly(
                        &res,
                        &[&current_vars[0..current_vars.len() - 1].join(",")],
                    ));
                }
            }
        }

        current_vars.pop();
        current_polys = next_set.into_iter().collect();
        projection_sets.push(current_polys.clone());
    }

    projection_sets.reverse(); // From R^1 to R^n-1
    Ok(projection_sets)
}

/// Performs the lifting phase of CAD.
pub(crate) fn lifting_phase(
    projections: &[Vec<SparsePolynomial>],
    vars: &[&str],
) -> Result<Vec<CadCell>, String> {
    // Base case: R^1
    let base_polys = &projections[0];
    let mut all_roots = Vec::new();
    for p in base_polys {
        let roots = isolate_real_roots(p, vars[0], 1e-9)?;
        all_roots.extend(roots.into_iter().map(|(a, b)| (a + b) / 2.0)); // Use midpoint as root
    }
    all_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    all_roots.dedup();

    let mut current_cells = Vec::new();
    let mut last_root = f64::NEG_INFINITY;

    for (i, root) in all_roots.iter().enumerate() {
        if *root > last_root && (root - last_root).abs() > 1e-9 {
            // Interval cell
            current_cells.push(CadCell {
                sample_point: vec![(last_root + root) / 2.0],
                dim: 1,
                index: vec![2 * i],
            });
        }
        // Point cell
        current_cells.push(CadCell {
            sample_point: vec![*root],
            dim: 0,
            index: vec![2 * i + 1],
        });
        last_root = *root;
    }
    // Final interval cell
    current_cells.push(CadCell {
        sample_point: vec![(last_root + last_root + 1.0)],
        dim: 1,
        index: vec![2 * all_roots.len()],
    });

    // Inductive lifting step
    for k in 1..vars.len() {
        let polys_k = &projections[k];
        let mut next_level_cells = Vec::new();
        for cell in &current_cells {
            // Substitute sample point into polynomials of level k+1
            let mut sample_map = HashMap::new();
            for i in 0..k {
                sample_map.insert(vars[i].to_string(), cell.sample_point[i]);
            }

            let mut roots_at_sample = Vec::new();
            for p in polys_k {
                let p_substituted_expr = substitute_map(&sparse_poly_to_expr(p), &sample_map);
                let p_substituted = expr_to_sparse_poly(&p_substituted_expr, &[vars[k]]);
                let roots = isolate_real_roots(&p_substituted, vars[k], 1e-9)?;
                roots_at_sample.extend(roots.into_iter().map(|(a, b)| (a + b) / 2.0));
            }
            roots_at_sample.sort_by(|a, b| a.partial_cmp(b).unwrap());
            roots_at_sample.dedup();

            // Create new cells (stack over the old cell)
            let mut last_root_val = f64::NEG_INFINITY;
            for (i, root_val) in roots_at_sample.iter().enumerate() {
                let mut new_sample = cell.sample_point.clone();
                if *root_val > last_root_val && (*root_val - last_root_val).abs() > 1e-9 {
                    new_sample.push((last_root_val + root_val) / 2.0);
                    let mut new_index = cell.index.clone();
                    new_index.push(2 * i);
                    next_level_cells.push(CadCell {
                        sample_point: new_sample.clone(),
                        dim: cell.dim + 1,
                        index: new_index,
                    });
                    new_sample.pop();
                }
                new_sample.push(*root_val);
                let mut new_index = cell.index.clone();
                new_index.push(2 * i + 1);
                next_level_cells.push(CadCell {
                    sample_point: new_sample,
                    dim: cell.dim,
                    index: new_index,
                });
                last_root_val = *root_val;
            }
            let mut final_sample = cell.sample_point.clone();
            final_sample.push(last_root_val + 1.0);
            let mut new_index = cell.index.clone();
            new_index.push(2 * roots_at_sample.len());
            next_level_cells.push(CadCell {
                sample_point: final_sample,
                dim: cell.dim + 1,
                index: new_index,
            });
        }
        current_cells = next_level_cells;
    }

    Ok(current_cells)
}

// ... (sylvester_matrix, resultant, and a new substitute_map helper are needed) ...

pub(crate) fn substitute_map(expr: &Expr, vars: &HashMap<String, f64>) -> Expr {
    let mut result = expr.clone();
    for (var, val) in vars {
        result = crate::symbolic::calculus::substitute(&result, var, &Expr::Constant(*val));
    }
    result
}

/// Computes the Sylvester matrix of two polynomials with respect to a given variable.
pub(crate) fn sylvester_matrix(p: &SparsePolynomial, q: &SparsePolynomial, var: &str) -> Expr {
    let n = p.degree(var) as usize;
    let m = q.degree(var) as usize;
    let mut matrix_rows = vec![vec![Expr::Constant(0.0); n + m]; n + m];

    let p_coeffs = p.get_coeffs_as_vec(var);
    let q_coeffs = q.get_coeffs_as_vec(var);

    for i in 0..m {
        for j in 0..=n {
            matrix_rows[i][i + j] = p_coeffs
                .get(j)
                .cloned()
                .unwrap_or_else(|| Expr::Constant(0.0));
        }
    }

    for i in 0..n {
        for j in 0..=m {
            matrix_rows[i + m][i + j] = q_coeffs
                .get(j)
                .cloned()
                .unwrap_or_else(|| Expr::Constant(0.0));
        }
    }

    Expr::Matrix(matrix_rows)
}

/// Computes the resultant of two polynomials with respect to a given variable.
pub(crate) fn resultant(p: &SparsePolynomial, q: &SparsePolynomial, var: &str) -> Expr {
    let sylvester = sylvester_matrix(p, q, var);
    matrix::determinant(&sylvester)
}
