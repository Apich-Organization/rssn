//! # Numerical differential geometry
//!
//! Curvature quantities of a metric given as a closure
//! `g(point) -> Vec<Vec<f64>>` returning the matrix `g_ij` at a point of the
//! coordinate chart. Derivatives of the metric are taken by fourth-order
//! central differences, so no symbolic form of the metric is needed.
//!
//! Index conventions: `Gamma[k][i][j]` is `Γ^k_{ij}`; `Riemann[r][s][m][n]`
//! is `R^r_{smn} = ∂_m Γ^r_{sn} - ∂_n Γ^r_{sm} + Γ^r_{lm} Γ^l_{sn} -
//! Γ^r_{ln} Γ^l_{sm}`; the Ricci tensor is `R_{sn} = R^m_{smn}`. With these
//! a sphere has a positive Ricci scalar.

/// Relative step of the finite differences. With a fourth-order stencil the
/// truncation error is about `h^4` and the rounding error `eps / h`, both
/// near `1e-12` here, which keeps the second derivatives needed for the
/// Riemann tensor accurate to around `1e-8`.
const STEP: f64 = 1.0e-3;

fn step(x: f64) -> f64 {
    STEP * x.abs().max(1.0)
}

/// Derivative of the matrix-valued `g` along `axis`: `∂_axis g_ij`, by the
/// fourth-order stencil `(-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / 12h`.
fn d_metric(
    g: &impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
    axis: usize,
) -> Vec<Vec<f64>> {
    let h = step(point[axis]);
    let at = |offset: f64| {
        let mut shifted = point.to_vec();
        shifted[axis] += offset * h;
        g(&shifted)
    };
    let (p2, p1, m1, m2) = (at(2.0), at(1.0), at(-1.0), at(-2.0));
    (0..p2.len())
        .map(|i| {
            (0..p2[i].len())
                .map(|j| (-p2[i][j] + 8.0 * p1[i][j] - 8.0 * m1[i][j] + m2[i][j]) / (12.0 * h))
                .collect()
        })
        .collect()
}

/// Inverse of a square matrix by Gauss-Jordan elimination with partial
/// pivoting.
fn invert(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = matrix.len();
    let mut a: Vec<Vec<f64>> = matrix
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut extended = row.clone();
            extended.extend((0..n).map(|j| f64::from(u8::from(i == j))));
            extended
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&p, &q| a[p][col].abs().total_cmp(&a[q][col].abs()))
            .unwrap_or(col);
        if a[pivot][col].abs() < 1e-14 {
            return Err("metric tensor is singular".to_string());
        }
        a.swap(col, pivot);
        let scale = a[col][col];
        for value in &mut a[col] {
            *value /= scale;
        }
        let pivot_row = a[col].clone();
        for (r, row) in a.iter_mut().enumerate() {
            if r != col {
                let factor = row[col];
                for (value, p) in row.iter_mut().zip(&pivot_row) {
                    *value -= factor * p;
                }
            }
        }
    }
    Ok(a.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// Evaluates the metric at `point` and checks that it is an `n x n` matrix
/// for `n = point.len()`.
///
/// # Errors
///
/// Returns an error if the metric is not square of the dimension of `point`.
pub fn metric_tensor_at_point(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    let metric = g(point);
    if metric.len() != point.len() || metric.iter().any(|row| row.len() != point.len()) {
        return Err(format!(
            "metric must be a {0}x{0} matrix at a point with {0} coordinates",
            point.len()
        ));
    }
    Ok(metric)
}

/// Christoffel symbols of the second kind `Γ^k_{ij}` at `point`:
/// `½ g^{km} (∂_j g_{mi} + ∂_i g_{mj} - ∂_m g_{ij})`.
///
/// # Errors
///
/// Returns an error if the metric has the wrong shape or is singular at
/// `point`.
pub fn christoffel_symbols(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<Vec<f64>>>, String> {
    let dim = point.len();
    let metric = metric_tensor_at_point(&g, point)?;
    let inverse = invert(&metric)?;
    // dg[k][i][j] = ∂_k g_ij
    let dg: Vec<Vec<Vec<f64>>> = (0..dim).map(|k| d_metric(&g, point, k)).collect();
    let mut gamma = vec![vec![vec![0.0; dim]; dim]; dim];
    for k in 0..dim {
        for i in 0..dim {
            for j in 0..dim {
                let sum: f64 = (0..dim)
                    .map(|m| inverse[k][m] * (dg[j][m][i] + dg[i][m][j] - dg[m][i][j]))
                    .sum();
                gamma[k][i][j] = 0.5 * sum;
            }
        }
    }
    Ok(gamma)
}

/// Riemann curvature tensor `R^ρ_{σμν}` at `point`.
///
/// # Errors
///
/// Returns an error if the metric has the wrong shape or is singular near
/// `point`.
pub fn riemann_tensor(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<Vec<Vec<f64>>>>, String> {
    let dim = point.len();
    let gamma = christoffel_symbols(&g, point)?;
    let mut d_gamma = Vec::with_capacity(dim);
    for mu in 0..dim {
        let h = step(point[mu]);
        let at = |offset: f64| -> Result<Vec<Vec<Vec<f64>>>, String> {
            let mut shifted = point.to_vec();
            shifted[mu] += offset * h;
            christoffel_symbols(&g, &shifted)
        };
        let (p2, p1, m1, m2) = (at(2.0)?, at(1.0)?, at(-1.0)?, at(-2.0)?);
        let mut d = vec![vec![vec![0.0; dim]; dim]; dim];
        for r in 0..dim {
            for s in 0..dim {
                for n in 0..dim {
                    d[r][s][n] =
                        (-p2[r][s][n] + 8.0 * p1[r][s][n] - 8.0 * m1[r][s][n] + m2[r][s][n])
                            / (12.0 * h);
                }
            }
        }
        d_gamma.push(d);
    }
    let mut riemann = vec![vec![vec![vec![0.0; dim]; dim]; dim]; dim];
    for rho in 0..dim {
        for sigma in 0..dim {
            for mu in 0..dim {
                for nu in 0..dim {
                    let products: f64 = (0..dim)
                        .map(|l| {
                            gamma[rho][l][mu] * gamma[l][sigma][nu]
                                - gamma[rho][l][nu] * gamma[l][sigma][mu]
                        })
                        .sum();
                    riemann[rho][sigma][mu][nu] =
                        d_gamma[mu][rho][sigma][nu] - d_gamma[nu][rho][sigma][mu] + products;
                }
            }
        }
    }
    Ok(riemann)
}

/// Ricci tensor `R_{σν} = R^μ_{σμν}` at `point`.
///
/// # Errors
///
/// Returns an error if the Riemann tensor cannot be computed.
pub fn ricci_tensor(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    let riemann = riemann_tensor(g, point)?;
    let dim = riemann.len();
    let mut ricci = vec![vec![0.0; dim]; dim];
    for sigma in 0..dim {
        for nu in 0..dim {
            ricci[sigma][nu] = (0..dim).map(|mu| riemann[mu][sigma][mu][nu]).sum();
        }
    }
    Ok(ricci)
}

/// Ricci scalar `R = g^{μν} R_{μν}` at `point`.
///
/// # Errors
///
/// Returns an error if the Ricci tensor cannot be computed or the metric is
/// singular.
pub fn ricci_scalar(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<f64, String> {
    let ricci = ricci_tensor(&g, point)?;
    let inverse = invert(&metric_tensor_at_point(&g, point)?)?;
    Ok(inverse
        .iter()
        .zip(&ricci)
        .map(|(g_row, r_row)| g_row.iter().zip(r_row).map(|(a, b)| a * b).sum::<f64>())
        .sum())
}
