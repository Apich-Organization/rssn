//! Adaptive integration over finite simplices, including tetrahedra.
//!
//! A Duffy transform maps the unit cube to barycentric coordinates of the
//! simplex. The transformed integrand is integrated with the adaptive
//! Genz–Malik rule. Its error estimate is heuristic, not an enclosure; both
//! its cost and its initial budget requirement grow exponentially with
//! dimension. For planar triangles, the specialised rule in
//! [`super::quadrature::integrate_triangle`] is also available.

use std::cell::Cell;

use super::quadrature::{Integral, cubature_genz_malik};

/// Failure to form or evaluate a simplex integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimplexError {
    /// Invalid dimension, vertex shape, coordinates or tolerances.
    InvalidInput,
    /// A geometric calculation or transformed integrand overflowed.
    NumericalOverflow,
    /// The integrand returned NaN or infinity at a quadrature node.
    NonFiniteIntegrand,
}

/// Integrate `f` over the simplex defined by `n + 1` vertices in `R^n`.
///
/// Dimensions 1 through 20 are accepted, although the underlying rule is
/// practical only in much lower dimensions. The result uses unsigned volume,
/// regardless of vertex orientation. A zero-volume simplex returns zero
/// without evaluating `f`. At least one tolerance must be positive.
///
/// The returned [`Integral::converged`] is false if the evaluation budget is
/// exhausted. If the budget cannot cover the initial grid, its value and error
/// are NaN and its evaluation count is zero, as in [`cubature_genz_malik`].
/// Error estimates are heuristic and cannot certify accuracy.
///
/// # Errors
/// Returns [`SimplexError::InvalidInput`] for malformed vertices or tolerances,
/// and reports non-finite integrands and arithmetic overflow separately.
pub fn integrate_simplex<F: Fn(&[f64]) -> f64>(
    f: F,
    vertices: &[Vec<f64>],
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<Integral, SimplexError> {
    let n = vertices
        .len()
        .checked_sub(1)
        .ok_or(SimplexError::InvalidInput)?;
    if !(1..=20).contains(&n)
        || vertices
            .iter()
            .any(|v| v.len() != n || v.iter().any(|x| !x.is_finite()))
        || !abs_tol.is_finite()
        || !rel_tol.is_finite()
        || abs_tol < 0.0
        || rel_tol < 0.0
        || (abs_tol == 0.0 && rel_tol == 0.0)
    {
        return Err(SimplexError::InvalidInput);
    }
    // Columns of the affine map are the edges from vertex zero.
    let mut edges = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let edge = vertices[j + 1][i] - vertices[0][i];
            if !edge.is_finite() {
                return Err(SimplexError::NumericalOverflow);
            }
            edges[i * n + j] = edge;
        }
    }
    let determinant = abs_determinant(&edges, n)?;
    if determinant == 0.0 {
        return Ok(Integral {
            value: 0.0,
            error: 0.0,
            evaluations: 0,
            converged: true,
        });
    }
    let bad_integrand = Cell::new(false);
    let overflow = Cell::new(false);
    let result = cubature_genz_malik(
        |u| {
            let mut x = vertices[0].clone();
            let mut remaining = 1.0;
            let mut jacobian = determinant;
            for j in 0..n {
                let weight = remaining * u[j];
                for i in 0..n {
                    x[i] += edges[i * n + j] * weight;
                }
                remaining *= 1.0 - u[j];
                if j + 1 < n {
                    jacobian *= remaining;
                }
            }
            if x.iter().any(|v| !v.is_finite()) || !jacobian.is_finite() {
                overflow.set(true);
                return f64::NAN;
            }
            let value = f(&x);
            if !value.is_finite() {
                bad_integrand.set(true);
                return f64::NAN;
            }
            let scaled = value * jacobian;
            if !scaled.is_finite() {
                overflow.set(true);
            }
            scaled
        },
        &vec![0.0; n],
        &vec![1.0; n],
        abs_tol,
        rel_tol,
        max_evals,
    );
    if bad_integrand.get() {
        return Err(SimplexError::NonFiniteIntegrand);
    }
    if overflow.get()
        || (result.evaluations != 0 && (!result.value.is_finite() || !result.error.is_finite()))
    {
        return Err(SimplexError::NumericalOverflow);
    }
    Ok(result)
}

// Gaussian elimination with partial pivoting. An underflowing non-zero
// determinant is reported instead of silently classifying it as degenerate.
fn abs_determinant(
    matrix: &[f64],
    n: usize,
) -> Result<f64, SimplexError> {
    let mut a = matrix.to_vec();
    let mut determinant = 1.0;
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        if a[pivot * n + k] == 0.0 {
            return Ok(0.0);
        }
        if pivot != k {
            for j in k..n {
                a.swap(k * n + j, pivot * n + j);
            }
        }
        let diagonal = a[k * n + k];
        determinant *= diagonal.abs();
        if !determinant.is_finite() || determinant == 0.0 {
            return Err(SimplexError::NumericalOverflow);
        }
        for i in k + 1..n {
            let multiplier = a[i * n + k] / diagonal;
            for j in k + 1..n {
                a[i * n + j] -= multiplier * a[k * n + j];
                if !a[i * n + j].is_finite() {
                    return Err(SimplexError::NumericalOverflow);
                }
            }
        }
    }
    Ok(determinant)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(
        a: f64,
        b: f64,
        tol: f64,
    ) {
        assert!((a - b).abs() < tol, "{a} != {b}");
    }

    #[test]
    fn tetrahedron_polynomial_moments() {
        let v = vec![
            vec![0., 0., 0.],
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
        ];
        let constant = integrate_simplex(|_| 1., &v, 1e-11, 1e-11, 200_000).unwrap();
        assert!(constant.converged, "{constant:?}");
        near(constant.value, 1.0 / 6.0, 1e-12);
        let moment =
            integrate_simplex(|x| x[0] * x[0] * x[1] * x[2], &v, 1e-10, 1e-10, 200_000).unwrap();
        assert!(moment.converged, "{moment:?}");
        near(moment.value, 1.0 / 2520.0, 1e-9);
    }

    #[test]
    fn affine_orientation_and_one_dimension() {
        let mut v = vec![vec![1., -1.], vec![3., -1.], vec![1., 2.]];
        let forward = integrate_simplex(|x| x[0] + 2. * x[1], &v, 1e-10, 1e-10, 100_000).unwrap();
        // Area 3, centroid (5/3, 0).
        near(forward.value, 5.0, 1e-10);
        v.swap(1, 2);
        let reversed = integrate_simplex(|x| x[0] + 2. * x[1], &v, 1e-10, 1e-10, 100_000).unwrap();
        near(reversed.value, forward.value, 1e-10);
        let line = integrate_simplex(
            |x| x[0] * x[0],
            &[vec![3.], vec![1.]],
            1e-10,
            1e-10,
            100_000,
        )
        .unwrap();
        near(line.value, 26.0 / 3.0, 1e-10);
    }

    #[test]
    fn invalid_degenerate_and_nonfinite_cases() {
        assert_eq!(
            integrate_simplex(|_| 0., &[vec![0.], vec![1., 2.]], 1e-6, 0., 100),
            Err(SimplexError::InvalidInput)
        );
        assert_eq!(
            integrate_simplex(|_| 0., &[vec![0.], vec![1.]], 0., 0., 100),
            Err(SimplexError::InvalidInput)
        );
        let v = vec![vec![0., 0.], vec![1., 1.], vec![2., 2.]];
        let q = integrate_simplex(
            |_| panic!("degenerate simplex evaluates integrand"),
            &v,
            1e-6,
            0.,
            100,
        )
        .unwrap();
        assert_eq!(q.value, 0.0);
        assert_eq!(q.evaluations, 0);
        assert!(q.converged);
        let v = vec![vec![0.], vec![1.]];
        assert_eq!(
            integrate_simplex(|_| f64::NAN, &v, 1e-6, 0., 1000),
            Err(SimplexError::NonFiniteIntegrand)
        );
        let q = integrate_simplex(|_| 1., &v, 1e-6, 0., 0).unwrap();
        assert!(!q.converged && q.value.is_nan() && q.evaluations == 0);
    }
}
