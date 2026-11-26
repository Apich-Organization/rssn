//! Internal dense matrix helpers shared by the control kernels.
//!
//! These routines wrap [`crate::kernels::dense`] and add the few
//! `Vec<Vec<f64>>` conveniences the control code needs. They are private to
//! the `control` module and carry no stability guarantees beyond those of the
//! dense backend.

use num_complex::Complex64;

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

/// Zero matrix with `rows` rows and `cols` columns.
pub fn zeros(
    rows: usize,
    cols: usize,
) -> Vec<Vec<f64>> {
    vec![vec![0.0; cols]; rows]
}

/// Identity matrix of size `n`.
pub fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut m = zeros(n, n);
    for i in 0..n {
        m[i][i] = 1.0;
    }
    m
}

/// Converts a row-major slice-of-rows matrix into a [`Mat`].
pub fn to_mat(a: &[Vec<f64>]) -> Result<Mat, String> {
    Mat::from_rows(a).map_err(|e| e.to_string())
}

/// Converts a [`Mat`] into a slice-of-rows matrix.
pub fn from_mat(m: &Mat) -> Vec<Vec<f64>> {
    (0..m.rows)
        .map(|i| (0..m.cols).map(|j| m.at(i, j)).collect())
        .collect()
}

/// Transpose of a rectangular slice-of-rows matrix.
pub fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let cols = a.first().map_or(0, Vec::len);
    let mut t = zeros(cols, a.len());
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate().take(cols) {
            t[j][i] = v;
        }
    }
    t
}

/// Matrix product, propagating shape errors.
pub fn matmul(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let ma = to_mat(a)?;
    let mb = to_mat(b)?;
    let c = ma.matmul(&mb).map_err(|e| e.to_string())?;
    Ok(from_mat(&c))
}

/// Matrix product that tolerates ragged input by treating missing entries as
/// zero; used by methods that cannot return a `Result`.
pub fn matmul_or_zero(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let cols = b.first().map_or(0, Vec::len);
    let mut out = zeros(a.len(), cols);
    for (i, row) in a.iter().enumerate() {
        for (k, &av) in row.iter().enumerate() {
            if av == 0.0 {
                continue;
            }
            if let Some(brow) = b.get(k) {
                for (j, &bv) in brow.iter().enumerate().take(cols) {
                    out[i][j] += av * bv;
                }
            }
        }
    }
    out
}

/// Matrix-vector product; entries beyond a short vector are treated as zero.
pub fn matvec(
    a: &[Vec<f64>],
    x: &[f64],
) -> Vec<f64> {
    a.iter()
        .map(|row| crate::kernels::simd::dot_f64(row, x))
        .collect()
}

/// Element-wise sum of two matrices of equal shape.
pub fn add(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let cols = a.first().map_or(0, Vec::len);
    if a.len() != b.len() || b.iter().any(|r| r.len() != cols) {
        return Err("matrix shapes do not match".to_string());
    }
    let mut out = zeros(a.len(), cols);
    for (i, row) in a.iter().enumerate() {
        crate::kernels::simd::add_into_f64(row, &b[i], &mut out[i]);
    }
    Ok(out)
}

/// Element-wise difference `a - b` of two matrices of equal shape.
pub fn sub(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let cols = a.first().map_or(0, Vec::len);
    if a.len() != b.len() || b.iter().any(|r| r.len() != cols) {
        return Err("matrix shapes do not match".to_string());
    }
    let mut out = zeros(a.len(), cols);
    for (i, row) in a.iter().enumerate() {
        crate::kernels::simd::sub_into_f64(row, &b[i], &mut out[i]);
    }
    Ok(out)
}

/// Solves the square system `a x = b` by LU with partial pivoting.
pub fn solve(
    a: &[Vec<f64>],
    b: &[f64],
) -> Result<Vec<f64>, String> {
    let ma = to_mat(a)?;
    dense::solve(&ma, b).map_err(|e| e.to_string())
}

/// Solves `a x = b` for a matrix right-hand side, reusing a single LU factorization.
pub fn solve_multi(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let cols = b.first().map_or(0, Vec::len);
    let ma = to_mat(a)?;
    let lu = dense::lu_factor(&ma).map_err(|e| e.to_string())?;
    let mut out = zeros(n, cols);
    if n == 0 || cols == 0 {
        return Ok(out);
    }
    let mut rhs = crate::kernels::reuse::f64_scratch(n);
    let mut col = crate::kernels::reuse::f64_scratch(n);
    rhs.resize(n, 0.0);
    col.resize(n, 0.0);
    for j in 0..cols {
        for i in 0..n {
            rhs[i] = b.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
        }
        lu.solve_into(&rhs, &mut col).map_err(|e| e.to_string())?;
        for i in 0..n {
            out[i][j] = col[i];
        }
    }
    Ok(out)
}

/// Inverse of a square matrix through LU solves.
pub fn inverse(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    if a.iter().any(|r| r.len() != n) {
        return Err("matrix must be square".to_string());
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    let ma = to_mat(a)?;
    let mut inv = zeros(n, n);
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        let col = dense::solve(&ma, &e).map_err(|e| e.to_string())?;
        for i in 0..n {
            inv[i][j] = col[i];
        }
    }
    Ok(inv)
}

/// Least-squares solution of `a x = b`.
pub fn lstsq(
    a: &[Vec<f64>],
    b: &[f64],
) -> Result<Vec<f64>, String> {
    let ma = to_mat(a)?;
    dense::lstsq(&ma, b).map_err(|e| e.to_string())
}

/// Numerical rank with a relative singular-value tolerance of `1e-9`.
pub fn rank_rel(a: &[Vec<f64>]) -> Result<usize, String> {
    let ma = to_mat(a)?;
    let sv = dense::svd(&ma).map_err(|e| e.to_string())?;
    let smax = sv.s.first().copied().unwrap_or(0.0);
    if !(smax > 0.0) {
        return Ok(0);
    }
    let tol = 1e-9 * smax;
    Ok(sv.s.iter().filter(|&&v| v > tol).count())
}

/// Eigenvalues of a general real matrix as complex numbers.
pub fn eigenvalues(a: &[Vec<f64>]) -> Result<Vec<Complex64>, String> {
    let ma = to_mat(a)?;
    let vals = dense::eigenvalues(&ma).map_err(|e| e.to_string())?;
    Ok(vals
        .into_iter()
        .map(|(re, im)| Complex64::new(re, im))
        .collect())
}

/// Eigenvalues and eigenvectors of a general real matrix.
pub fn eigen_pairs(a: &[Vec<f64>]) -> Result<(Vec<(f64, f64)>, Vec<Vec<Complex64>>), String> {
    let ma = to_mat(a)?;
    let e = dense::eigen(&ma).map_err(|err| err.to_string())?;
    Ok((e.values, e.vectors))
}

/// Frobenius norm of a matrix.
pub fn frobenius(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .flat_map(|r| r.iter())
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt()
}

/// One classical fourth-order Runge-Kutta step of `x' = A x + B u` with the
/// input held constant over the step.
pub fn rk4_linear(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    x: &[f64],
    u: &[f64],
    h: f64,
) -> Vec<f64> {
    let f = |state: &[f64]| {
        let ax = matvec(a, state);
        let bu = matvec(b, u);
        ax.iter()
            .enumerate()
            .map(|(i, &v)| v + bu.get(i).copied().unwrap_or(0.0))
            .collect::<Vec<f64>>()
    };
    let k1 = f(x);
    let x2: Vec<f64> = x
        .iter()
        .zip(&k1)
        .map(|(&xi, &ki)| xi + 0.5 * h * ki)
        .collect();
    let k2 = f(&x2);
    let x3: Vec<f64> = x
        .iter()
        .zip(&k2)
        .map(|(&xi, &ki)| xi + 0.5 * h * ki)
        .collect();
    let k3 = f(&x3);
    let x4: Vec<f64> = x.iter().zip(&k3).map(|(&xi, &ki)| xi + h * ki).collect();
    let k4 = f(&x4);
    x.iter()
        .enumerate()
        .map(|(i, &xi)| {
            let inc = k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i];
            xi + h * inc / 6.0
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matmul_and_transpose() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let b = vec![vec![5.0, 6.0], vec![7.0, 8.0]];
        let c = matmul(&a, &b).unwrap();
        assert!((c[0][0] - 19.0).abs() < 1e-12);
        assert!((c[1][1] - 50.0).abs() < 1e-12);
        let t = transpose(&a);
        assert!((t[0][1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_round_trip() {
        let a = vec![vec![4.0, 7.0], vec![2.0, 6.0]];
        let inv = inverse(&a).unwrap();
        let id = matmul(&a, &inv).unwrap();
        assert!((id[0][0] - 1.0).abs() < 1e-12);
        assert!((id[1][0]).abs() < 1e-12);
        assert!((id[0][1]).abs() < 1e-12);
        assert!((id[1][1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rank_of_deficient_matrix() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert_eq!(rank_rel(&a).unwrap(), 1);
    }

    #[test]
    fn rk4_matches_exp_for_scalar() {
        let a = vec![vec![-1.0]];
        let b = vec![vec![1.0]];
        let x = rk4_linear(&a, &b, &[0.0], &[1.0], 0.01);
        // x(h) = 1 - exp(-h) for x' = -x + 1, x(0) = 0.
        let exact = 1.0 - (-0.01_f64).exp();
        assert!((x[0] - exact).abs() < 1e-9);
    }

    #[test]
    fn solve_multi_reuses_lu_and_checks_errors() {
        let a = vec![
            vec![3.0, 1.0, -1.0],
            vec![2.0, 4.0, 1.0],
            vec![-1.0, 2.0, 5.0],
        ];
        let b = vec![vec![4.0, 1.0], vec![1.0, 2.0], vec![1.0, -3.0]];
        let x = solve_multi(&a, &b).unwrap();
        let ax = matmul(&a, &x).unwrap();
        for i in 0..3 {
            for j in 0..2 {
                assert!((ax[i][j] - b[i][j]).abs() < 1e-11);
            }
        }
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(solve_multi(&singular, &[vec![1.0], vec![2.0]]).is_err());
        let bad_rhs = vec![vec![1.0], vec![f64::NAN], vec![0.0]];
        assert!(solve_multi(&a, &bad_rhs).is_err());
    }
}
