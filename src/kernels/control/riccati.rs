//! Algebraic Riccati equations and the optimal gains they define.
//!
//! The continuous equation is solved by the stable invariant subspace of the
//! Hamiltonian matrix, with a Kleinman (Newton) refinement when the computed
//! residual is not at machine level. The discrete equation is solved by
//! monotone value iteration.

use super::linalg;

/// Checks that `m` is square and returns its order.
fn check_square(
    m: &[Vec<f64>],
    name: &str,
) -> Result<usize, String> {
    let n = m.len();
    if m.iter().any(|row| row.len() != n) {
        return Err(format!("{name} must be square"));
    }
    Ok(n)
}

/// Symmetrises a square matrix in place.
fn symmetrise(m: &mut [Vec<f64>]) {
    let n = m.len();
    for i in 0..n {
        for j in i + 1..n {
            let avg = f64::midpoint(m[i][j], m[j][i]);
            m[i][j] = avg;
            m[j][i] = avg;
        }
    }
}

/// Solves `a x = b` for a matrix right-hand side through the shared helper.
fn solve_multi(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    linalg::solve_multi(a, b)
}

/// Builds the block matrix `[[a, -g], [-q, -a^T]]`.
fn hamiltonian(
    a: &[Vec<f64>],
    g: &[Vec<f64>],
    q: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut h = linalg::zeros(2 * n, 2 * n);
    for i in 0..n {
        for j in 0..n {
            h[i][j] = a[i][j];
            h[i][n + j] = -g[i][j];
            h[n + i][j] = -q[i][j];
            h[n + i][n + j] = -a[j][i];
        }
    }
    h
}

/// Solves the continuous Lyapunov equation `a^T p + p a = -w`.
fn lyapunov(
    a: &[Vec<f64>],
    w: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let size = n * n;
    let mut lhs = vec![vec![0.0; size]; size];
    for i in 0..n {
        for j in 0..n {
            let row = j * n + i;
            for k in 0..n {
                lhs[row][j * n + k] += a[k][i];
                lhs[row][k * n + i] += a[k][j];
            }
        }
    }
    let mut rhs = vec![0.0; size];
    for i in 0..n {
        for j in 0..n {
            rhs[j * n + i] = -w[i][j];
        }
    }
    let solution = linalg::solve(&lhs, &rhs)?;
    let mut p = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            p[i][j] = solution[j * n + i];
        }
    }
    Ok(p)
}

/// Relative residual of the continuous Riccati equation.
fn care_residual(
    a: &[Vec<f64>],
    g: &[Vec<f64>],
    q: &[Vec<f64>],
    p: &[Vec<f64>],
) -> f64 {
    let at = linalg::transpose(a);
    let at_p = linalg::matmul_or_zero(&at, p);
    let p_a = linalg::matmul_or_zero(p, a);
    let pg = linalg::matmul_or_zero(p, g);
    let pgp = linalg::matmul_or_zero(&pg, p);
    let n = a.len();
    let mut residual = vec![vec![0.0; n]; n];
    let mut scale = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            residual[i][j] = at_p[i][j] + p_a[i][j] - pgp[i][j] + q[i][j];
            scale += (at_p[i][j].abs() + p_a[i][j].abs() + pgp[i][j].abs() + q[i][j].abs()).powi(2);
        }
    }
    let denom = scale.sqrt().max(f64::MIN_POSITIVE);
    linalg::frobenius(&residual) / denom
}

/// Kleinman/Newton refinement of a continuous Riccati solution.
fn kleinman(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    r: &[Vec<f64>],
    q: &[Vec<f64>],
    start: Vec<Vec<f64>>,
) -> Result<Vec<Vec<f64>>, String> {
    let bt = linalg::transpose(b);
    let mut p = start;
    for _ in 0..60 {
        let k = linalg::matmul(&linalg::matmul(&linalg::inverse(r)?, &bt)?, &p)?;
        let bk = linalg::matmul(b, &k)?;
        let closed = linalg::sub(a, &bk)?;
        let kt = linalg::transpose(&k);
        let krk = linalg::matmul(&linalg::matmul(&kt, r)?, &k)?;
        let w = linalg::add(q, &krk)?;
        let next = lyapunov(&closed, &w)?;
        let mut diff = 0.0_f64;
        let mut scale = 0.0_f64;
        for i in 0..p.len() {
            for j in 0..p.len() {
                diff = diff.max((next[i][j] - p[i][j]).abs());
                scale = scale.max(next[i][j].abs());
            }
        }
        p = next;
        if diff <= 1e-13 * scale.max(1.0) {
            break;
        }
    }
    Ok(p)
}

/// Solves the continuous algebraic Riccati equation
/// `A^T P + P A - P B R^{-1} B^T P + Q = 0`.
///
/// The solution is the stabilizing one, obtained from the stable invariant
/// subspace of the Hamiltonian matrix. A Kleinman iteration refines it when
/// the eigenvector solution leaves a residual above `1e-6` relative.
///
/// # Errors
/// Returns an error for inconsistent dimensions, a singular `R`, or when the
/// Hamiltonian eigensystem cannot be formed.
pub fn care(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let n = check_square(a, "a")?;
    if n == 0 {
        return Ok(Vec::new());
    }
    if b.len() != n {
        return Err("b must have as many rows as a".to_string());
    }
    let nu = b.first().map_or(0, Vec::len);
    if nu == 0 {
        return Err("b must have at least one column".to_string());
    }
    if b.iter().any(|row| row.len() != nu) {
        return Err("b rows must have equal length".to_string());
    }
    if check_square(q, "q")? != n {
        return Err("q must have the same order as a".to_string());
    }
    if check_square(r, "r")? != nu {
        return Err("r must match the number of inputs".to_string());
    }
    let mut q = q.to_vec();
    symmetrise(&mut q);
    let r_inv = linalg::inverse(r)?;
    let bt = linalg::transpose(b);
    let g = linalg::matmul(&linalg::matmul(b, &r_inv)?, &bt)?;
    let h = hamiltonian(a, &g, &q);
    // Two independent routes to the stabilising solution. The eigenvector of
    // the Hamiltonian's stable subspace is exact when the eigensolver
    // converges, but the Hamiltonian is non-symmetric with a symmetric
    // spectrum and a QR iteration can refuse it; the matrix sign function
    // needs only inverses and always converges when no eigenvalue sits on the
    // imaginary axis, which is exactly the well-posed case.
    let mut p = match eigen_solution(&h, n) {
        | Some(candidate) if care_residual(a, &g, &q, &candidate) <= 1e-6 => candidate,
        | _ => sign_solution(&h, n)?,
    };
    symmetrise(&mut p);
    if care_residual(a, &g, &q, &p) > 1e-10 {
        p = kleinman(a, b, r, &q, p)?;
        symmetrise(&mut p);
    }
    if care_residual(a, &g, &q, &p) > 1e-6 {
        return Err("the Riccati equation did not converge".to_string());
    }
    Ok(p)
}

/// The stabilising solution from the stable invariant subspace of the
/// Hamiltonian, or `None` when the eigensolver declines or the subspace is
/// too small.
fn eigen_solution(
    h: &[Vec<f64>],
    n: usize,
) -> Option<Vec<Vec<f64>>> {
    let (values, vectors) = linalg::eigen_pairs(h).ok()?;
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for (index, &(re, im)) in values.iter().enumerate() {
        if re >= -1e-9 {
            continue;
        }
        let vector = vectors.get(index)?;
        if im.abs() <= 1e-10 {
            basis.push(vector.iter().map(|c| c.re).collect());
        } else if im > 0.0 {
            basis.push(vector.iter().map(|c| c.re).collect());
            basis.push(vector.iter().map(|c| c.im).collect());
        }
    }
    if basis.len() < n {
        return None;
    }
    basis.truncate(n);
    let k = basis.len();
    let mut upper = linalg::zeros(n, k);
    let mut lower = linalg::zeros(n, k);
    for (column, vector) in basis.iter().enumerate() {
        for i in 0..n {
            upper[i][column] = vector[i];
            lower[i][column] = vector[n + i];
        }
    }
    let upper_t = linalg::transpose(&upper);
    let lower_t = linalg::transpose(&lower);
    let mut transpose = linalg::zeros(n, n);
    for j in 0..n {
        let rhs: Vec<f64> = (0..k).map(|r| lower_t[r][j]).collect();
        let column = linalg::lstsq(&upper_t, &rhs).ok()?;
        for i in 0..n {
            transpose[i][j] = column[i];
        }
    }
    Some(linalg::transpose(&transpose))
}

/// The stabilising solution from the matrix sign function of the Hamiltonian.
///
/// `sign(H)` is computed by the Newton iteration `Z <- (Z + Z^-1)/2`, which
/// converges quadratically to a matrix with eigenvalues `+-1`. The projector
/// `(I - sign(H))/2` is then the spectral projector onto the stable invariant
/// subspace, whose graph is the solution: with `[U; V]` spanning it,
/// `P = V U^-1`.
fn sign_solution(
    h: &[Vec<f64>],
    n: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let size = 2 * n;
    let mut z = h.to_vec();
    for _ in 0..200 {
        let inverse = linalg::inverse(&z)?;
        let mut next = linalg::add(&z, &inverse)?;
        for row in &mut next {
            for entry in row.iter_mut() {
                *entry *= 0.5;
            }
        }
        let mut diff = 0.0_f64;
        let mut scale = 1.0_f64;
        for i in 0..size {
            for j in 0..size {
                diff = diff.max((next[i][j] - z[i][j]).abs());
                scale = scale.max(next[i][j].abs());
            }
        }
        z = next;
        if diff <= 1e-14 * scale {
            break;
        }
    }
    // Projector onto the stable subspace: (I - sign(H)) / 2.
    let mut projector = vec![vec![0.0; size]; size];
    for i in 0..size {
        for j in 0..size {
            let identity = if i == j { 1.0 } else { 0.0 };
            projector[i][j] = 0.5 * (identity - z[i][j]);
        }
    }
    let mut upper = linalg::zeros(n, n);
    let mut lower = linalg::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            upper[i][j] = projector[i][j];
            lower[i][j] = projector[n + i][j];
        }
    }
    let upper_t = linalg::transpose(&upper);
    let lower_t = linalg::transpose(&lower);
    let transpose = solve_multi(&upper_t, &lower_t)?;
    Ok(linalg::transpose(&transpose))
}

/// Solves the discrete algebraic Riccati equation
/// `P = A^T P A - A^T P B (R + B^T P B)^{-1} B^T P A + Q`.
///
/// Uses monotone value iteration from `P = Q` with a `1e-12` stopping
/// tolerance and an iteration cap of 10 000.
///
/// # Errors
/// Returns an error for inconsistent dimensions or a singular `R + B^T P B`.
pub fn dare(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let n = check_square(a, "a")?;
    if n == 0 {
        return Ok(Vec::new());
    }
    if b.len() != n {
        return Err("b must have as many rows as a".to_string());
    }
    let nu = b.first().map_or(0, Vec::len);
    if nu == 0 {
        return Err("b must have at least one column".to_string());
    }
    if b.iter().any(|row| row.len() != nu) {
        return Err("b rows must have equal length".to_string());
    }
    if check_square(q, "q")? != n {
        return Err("q must have the same order as a".to_string());
    }
    if check_square(r, "r")? != nu {
        return Err("r must match the number of inputs".to_string());
    }
    let mut q = q.to_vec();
    symmetrise(&mut q);
    let bt = linalg::transpose(b);
    let mut p = q.clone();
    for _ in 0..10_000 {
        let ap = linalg::matmul(a, &p)?;
        let apa = linalg::matmul(&ap, a)?;
        let apb = linalg::matmul(&ap, b)?;
        let btp = linalg::matmul(&bt, &p)?;
        let btpb = linalg::matmul(&btp, b)?;
        let m = linalg::add(r, &btpb)?;
        let btpa = linalg::matmul(&btp, a)?;
        let k = solve_multi(&m, &btpa)?;
        let apbk = linalg::matmul(&apb, &k)?;
        let mut next = linalg::add(&linalg::sub(&apa, &apbk)?, &q)?;
        symmetrise(&mut next);
        let mut diff = 0.0_f64;
        let mut scale = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                diff = diff.max((next[i][j] - p[i][j]).abs());
                scale = scale.max(next[i][j].abs());
            }
        }
        p = next;
        if diff <= 1e-12 * scale.max(1.0) {
            break;
        }
    }
    Ok(p)
}

/// Continuous linear-quadratic regulator: returns `(K, P)` with
/// `K = R^{-1} B^T P`.
///
/// # Errors
/// As [`care`].
pub fn lqr(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    let p = care(a, b, q, r)?;
    let bt = linalg::transpose(b);
    let r_inv = linalg::inverse(r)?;
    let k = linalg::matmul(&linalg::matmul(&r_inv, &bt)?, &p)?;
    Ok((k, p))
}

/// Discrete linear-quadratic regulator: returns `(K, P)` with
/// `K = (R + B^T P B)^{-1} B^T P A`.
///
/// # Errors
/// As [`dare`].
pub fn dlqr(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    let p = dare(a, b, q, r)?;
    let bt = linalg::transpose(b);
    let btp = linalg::matmul(&bt, &p)?;
    let btpb = linalg::matmul(&btp, b)?;
    let m = linalg::add(r, &btpb)?;
    let btpa = linalg::matmul(&btp, a)?;
    let k = solve_multi(&m, &btpa)?;
    Ok((k, p))
}

/// Steady-state Kalman gain from the dual Riccati equation.
///
/// `W` is the process-noise covariance and `V` the measurement-noise
/// covariance. Returns `L = P C^T V^{-1}` where `P` solves
/// `A P + P A^T - P C^T V^{-1} C P + W = 0`.
///
/// # Errors
/// Returns an error for inconsistent dimensions or a singular `V`.
pub fn kalman_gain(
    a: &[Vec<f64>],
    c: &[Vec<f64>],
    w: &[Vec<f64>],
    v: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let n = check_square(a, "a")?;
    if c.is_empty() || c.iter().any(|row| row.len() != n) {
        return Err("c must have n columns".to_string());
    }
    let at = linalg::transpose(a);
    let ct = linalg::transpose(c);
    let p = care(&at, &ct, w, v)?;
    let pct = linalg::matmul(&p, &ct)?;
    let v_inv = linalg::inverse(v)?;
    linalg::matmul(&pct, &v_inv)
}

/// Linear-quadratic-Gaussian controller: returns `(LQR gain, Kalman gain)`.
///
/// # Errors
/// As [`lqr`] and [`kalman_gain`].
pub fn lqg(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    q: &[Vec<f64>],
    r: &[Vec<f64>],
    w: &[Vec<f64>],
    v: &[Vec<f64>],
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    let (k, _) = lqr(a, b, q, r)?;
    let l = kalman_gain(a, c, w, v)?;
    Ok((k, l))
}

#[cfg(test)]
mod tests {
    use super::linalg;
    use super::*;

    #[test]
    fn care_double_integrator_matches_analytic() {
        let a = vec![vec![0.0, 1.0], vec![0.0, 0.0]];
        let b = vec![vec![0.0], vec![1.0]];
        let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let r = vec![vec![1.0]];
        let p = care(&a, &b, &q, &r).unwrap();
        let s3 = 3.0_f64.sqrt();
        assert!((p[0][0] - s3).abs() < 1e-9);
        assert!((p[0][1] - 1.0).abs() < 1e-9);
        assert!((p[1][0] - 1.0).abs() < 1e-9);
        assert!((p[1][1] - s3).abs() < 1e-9);
        let (k, _) = lqr(&a, &b, &q, &r).unwrap();
        assert!((k[0][0] - 1.0).abs() < 1e-9);
        assert!((k[0][1] - s3).abs() < 1e-9);
    }

    #[test]
    fn dare_scalar_matches_closed_form() {
        let a = vec![vec![1.0]];
        let b = vec![vec![1.0]];
        let q = vec![vec![1.0]];
        let r = vec![vec![1.0]];
        let p = dare(&a, &b, &q, &r).unwrap();
        let golden = f64::midpoint(1.0, 5.0_f64.sqrt());
        assert!((p[0][0] - golden).abs() < 1e-9);
    }

    #[test]
    fn riccati_residual_is_small() {
        let a = vec![vec![-1.0, 0.5], vec![0.0, -2.0]];
        let b = vec![vec![1.0], vec![0.5]];
        let q = vec![vec![2.0, 0.0], vec![0.0, 1.0]];
        let r = vec![vec![0.5]];
        let p = care(&a, &b, &q, &r).unwrap();
        let r_inv = linalg::inverse(&r).unwrap();
        let bt = linalg::transpose(&b);
        let g = linalg::matmul(&linalg::matmul(&b, &r_inv).unwrap(), &bt).unwrap();
        assert!(care_residual(&a, &g, &q, &p) < 1e-10);
    }

    #[test]
    fn kalman_gain_of_scalar_system() {
        let a = vec![vec![1.0]];
        let c = vec![vec![1.0]];
        let w = vec![vec![0.1]];
        let v = vec![vec![0.5]];
        let l = kalman_gain(&a, &c, &w, &v).unwrap();
        assert!(l[0][0].is_finite() && l[0][0] > 0.0);
    }
}
