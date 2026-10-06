//! # Structure-preserving and stiff integrators
//!
//! * **Symplectic** integrators for separable Hamiltonians
//!   `H(q, p) = T(p) + V(q)` (here `T = |p|²/2m`): Störmer–Verlet (order
//!   2) and the Forest–Ruth / Yoshida composition (order 4). They conserve
//!   a modified energy, so the energy error stays bounded over very long
//!   runs instead of drifting.
//! * **Stiff** integrators for `y' = f(t, y)`: backward Euler and the
//!   implicit midpoint rule (Newton iterations with a finite-difference
//!   Jacobian), and an adaptive L-stable Rosenbrock method (ROS2, with the
//!   linearly implicit Euler step as error estimate), which needs one
//!   Jacobian and two linear solves per step and no Newton iteration.
//!
//! Linear systems are solved by Gaussian elimination with partial pivoting
//! (the systems are small: one row per state variable).

#[allow(clippy::needless_range_loop)] // row operations index two rows at once
/// Solves `a x = b` in place by Gaussian elimination with partial
/// pivoting; `None` for a (numerically) singular matrix.
fn solve_dense(
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            if factor == 0.0 {
                continue;
            }
            for k in col..n {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    Some(x)
}

/// The Jacobian `∂f/∂y` at `(t, y)` by central differences.
fn jacobian<F>(
    f: &F,
    t: f64,
    y: &[f64],
) -> Vec<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let mut jac = vec![vec![0.0; n]; n];
    let mut shifted = y.to_vec();
    for j in 0..n {
        let h = 1e-7 * (1.0 + y[j].abs());
        shifted[j] = y[j] + h;
        let plus = f(t, &shifted);
        shifted[j] = y[j] - h;
        let minus = f(t, &shifted);
        shifted[j] = y[j];
        for i in 0..n {
            jac[i][j] = (plus[i] - minus[i]) / (2.0 * h);
        }
    }
    jac
}

/// One step of the Störmer–Verlet (leapfrog) method for `q'' = a(q)`
/// (`a` the acceleration, `-∇V/m`): kick, drift, kick.
pub fn verlet_step<A>(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    dt: f64,
    acceleration: &A,
) where
    A: Fn(&[f64]) -> Vec<f64>,
{
    let a = acceleration(q);
    for (pi, ai) in p.iter_mut().zip(&a) {
        *pi += 0.5 * dt * mass * ai;
    }
    for (qi, pi) in q.iter_mut().zip(p.iter()) {
        *qi += dt * pi / mass;
    }
    let a = acceleration(q);
    for (pi, ai) in p.iter_mut().zip(&a) {
        *pi += 0.5 * dt * mass * ai;
    }
}

/// One step of the fourth-order Forest–Ruth (Yoshida) composition of
/// Verlet steps with weights `w1, w0, w1`, `w1 = 1/(2 - 2^(1/3))`.
pub fn yoshida4_step<A>(
    q: &mut [f64],
    p: &mut [f64],
    mass: f64,
    dt: f64,
    acceleration: &A,
) where
    A: Fn(&[f64]) -> Vec<f64>,
{
    let cbrt2 = 2.0_f64.cbrt();
    let w1 = 1.0 / (2.0 - cbrt2);
    let w0 = -cbrt2 * w1;
    for w in [w1, w0, w1] {
        verlet_step(q, p, mass, w * dt, acceleration);
    }
}

/// Integrates a separable Hamiltonian system for `steps` steps with the
/// symplectic method of the given `order` (2: Verlet, 4: Yoshida),
/// returning the states `(q, p)` after every step.
#[must_use]
pub fn symplectic_trajectory<A>(
    q0: &[f64],
    p0: &[f64],
    mass: f64,
    dt: f64,
    steps: usize,
    order: u32,
    acceleration: A,
) -> Vec<(Vec<f64>, Vec<f64>)>
where
    A: Fn(&[f64]) -> Vec<f64>,
{
    let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
    let mut out = Vec::with_capacity(steps + 1);
    out.push((q.clone(), p.clone()));
    for _ in 0..steps {
        if order >= 4 {
            yoshida4_step(&mut q, &mut p, mass, dt, &acceleration);
        } else {
            verlet_step(&mut q, &mut p, mass, dt, &acceleration);
        }
        out.push((q.clone(), p.clone()));
    }
    out
}

/// Solves `g(z) = 0` by Newton's method from `z0` with the
/// finite-difference Jacobian of `g`.
fn newton<G>(
    g: &G,
    z0: Vec<f64>,
) -> Option<Vec<f64>>
where
    G: Fn(f64, &[f64]) -> Vec<f64>,
{
    let mut z = z0;
    for _ in 0..50 {
        let r = g(0.0, &z);
        let norm: f64 = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-12 * (1.0 + z.iter().map(|v| v.abs()).fold(0.0, f64::max)) {
            return Some(z);
        }
        let jac = jacobian(g, 0.0, &z);
        let delta = solve_dense(jac, r.iter().map(|v| -v).collect())?;
        for (zi, di) in z.iter_mut().zip(&delta) {
            *zi += di;
        }
        if delta.iter().map(|v| v.abs()).fold(0.0, f64::max) < 1e-14 {
            return Some(z);
        }
    }
    Some(z)
}

/// One backward Euler step: `y1 = y0 + h f(t + h, y1)`.
pub fn backward_euler_step<F>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
) -> Option<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let residual = |_: f64, z: &[f64]| -> Vec<f64> {
        let fz = f(t + h, z);
        z.iter().zip(y).zip(&fz).map(|((zi, yi), fi)| zi - yi - h * fi).collect()
    };
    newton(&residual, y.to_vec())
}

/// One implicit midpoint step: `y1 = y0 + h f(t + h/2, (y0 + y1)/2)`
/// (symplectic and A-stable, order 2).
pub fn implicit_midpoint_step<F>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
) -> Option<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let residual = |_: f64, z: &[f64]| -> Vec<f64> {
        let mid: Vec<f64> = z.iter().zip(y).map(|(zi, yi)| 0.5 * (zi + yi)).collect();
        let fm = f(t + 0.5 * h, &mid);
        z.iter().zip(y).zip(&fm).map(|((zi, yi), fi)| zi - yi - h * fi).collect()
    };
    newton(&residual, y.to_vec())
}

/// Integrates `y' = f(t, y)` from `t0` to `t1` with the adaptive Rosenbrock
/// method ROS2 (L-stable, order 2; `γ = 1 + 1/√2`), the error estimated
/// against the linearly implicit Euler step.
///
/// Returns the accepted `(t, y)` pairs, or `None` when a linear system is
/// singular or the step size underflows.
pub fn rosenbrock_solve<F>(
    f: F,
    t0: f64,
    y0: &[f64],
    t1: f64,
    tolerance: f64,
) -> Option<Vec<(f64, Vec<f64>)>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let gamma = 1.0 + 1.0 / 2.0_f64.sqrt();
    let n = y0.len();
    let (mut t, mut y) = (t0, y0.to_vec());
    let mut h = ((t1 - t0) / 100.0).abs().max(1e-12) * (t1 - t0).signum();
    let mut out = vec![(t, y.clone())];
    let mut guard = 0_usize;
    loop {
        if (t1 - t) * h.signum() <= 1e-14 * (1.0 + t1.abs()) {
            break;
        }
        guard += 1;
        if guard > 1_000_000 || h.abs() < 1e-14 * (1.0 + t.abs()) {
            return None;
        }
        if (t + h - t1) * h.signum() > 0.0 {
            h = t1 - t;
        }
        let jac = jacobian(&f, t, &y);
        // W = I - γ h J
        let w: Vec<Vec<f64>> =
            (0..n).map(|i| (0..n).map(|j| f64::from(u8::from(i == j)) - gamma * h * jac[i][j]).collect()).collect();
        let f0 = f(t, &y);
        let k1 = solve_dense(w.clone(), f0)?;
        let y_mid: Vec<f64> = y.iter().zip(&k1).map(|(yi, ki)| yi + h * ki).collect();
        let f1 = f(t + h, &y_mid);
        let rhs: Vec<f64> = f1.iter().zip(&k1).map(|(fi, ki)| fi - 2.0 * ki).collect();
        let k2 = solve_dense(w, rhs)?;
        let next: Vec<f64> = y.iter().zip(&k1).zip(&k2).map(|((yi, a), b)| yi + 1.5 * h * a + 0.5 * h * b).collect();
        // Error against the first-order (linearly implicit Euler) solution.
        let error = next
            .iter()
            .zip(&y_mid)
            .zip(&y)
            .map(|((a, b), y0)| ((a - b) / (tolerance * (1.0 + y0.abs().max(a.abs())))).powi(2))
            .sum::<f64>()
            .sqrt()
            / (n as f64).sqrt();
        if error <= 1.0 || !error.is_finite() && h.abs() < 1e-10 {
            t += h;
            y = next;
            out.push((t, y.clone()));
        }
        let factor = if error.is_finite() { (0.9 / error.max(1e-10).sqrt()).clamp(0.2, 5.0) } else { 0.2 };
        h *= factor;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symplectic_energy_stays_bounded() {
        // Kepler orbit with eccentricity 0.5: energy error bounded over
        // 100 periods for both methods, and far smaller at order 4.
        let accel = |q: &[f64]| {
            let r3 = (q[0] * q[0] + q[1] * q[1]).powf(1.5);
            vec![-q[0] / r3, -q[1] / r3]
        };
        let energy = |q: &[f64], p: &[f64]| 0.5 * p[0].hypot(p[1]).powi(2) - 1.0 / q[0].hypot(q[1]);
        let (q0, p0) = ([0.5, 0.0], [0.0, 3.0_f64.sqrt()]);
        let e0 = energy(&q0, &p0);
        let period = 2.0 * std::f64::consts::PI; // a = 1
        let dt = period / 400.0;
        let steps = 400 * 100;
        let mut worst = [0.0_f64; 2];
        for (k, order) in [2_u32, 4].into_iter().enumerate() {
            let path = symplectic_trajectory(&q0, &p0, 1.0, dt, steps, order, accel);
            for (q, p) in path.iter().step_by(97) {
                worst[k] = worst[k].max((energy(q, p) - e0).abs());
            }
            // At order 4 the orbit still closes after 100 periods (Verlet's
            // orbit precesses).
            if order == 4 {
                let (q, _) = path.last().unwrap();
                assert!((q[0] - 0.5).abs() < 1e-3 && q[1].abs() < 1e-2, "{q:?}");
            }
        }
        assert!(worst[0] < 2e-2 && worst[1] < 1e-5, "{worst:?}");
    }

    #[test]
    fn stiff_problems() {
        // Robertson's chemical kinetics: very stiff; mass is conserved and
        // the known values at t = 40 are matched.
        let robertson = |_: f64, y: &[f64]| {
            vec![
                -0.04 * y[0] + 1e4 * y[1] * y[2],
                0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1],
                3e7 * y[1] * y[1],
            ]
        };
        let path = rosenbrock_solve(robertson, 0.0, &[1.0, 0.0, 0.0], 40.0, 1e-6).unwrap();
        let (_, y) = path.last().unwrap();
        assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-6, "{y:?}");
        assert!((y[0] - 0.715_8).abs() < 2e-3 && (y[2] - 0.284_1).abs() < 2e-3, "{y:?}");
        assert!(path.len() < 2000, "{} steps", path.len());
        // y' = -1000 (y - cos t): backward Euler and implicit midpoint stay
        // stable with h = 0.1 (explicit Euler would explode).
        let f = |t: f64, y: &[f64]| vec![-1000.0 * (y[0] - t.cos())];
        let (mut a, mut b) = (vec![0.0], vec![0.0]);
        for k in 0..50 {
            let t = 0.1 * f64::from(k);
            a = backward_euler_step(&f, t, &a, 0.1).unwrap();
            b = implicit_midpoint_step(&f, t, &b, 0.1).unwrap();
        }
        assert!((a[0] - 5.0_f64.cos()).abs() < 1e-2 && b[0].abs() < 2.0, "{a:?} {b:?}");
    }
}
