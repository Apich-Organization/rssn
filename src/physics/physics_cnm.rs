/// Solves a system of linear equations Ax = d where A is a tridiagonal matrix.
/// `a`: sub-diagonal (n-1 elements), `b`: main diagonal (n elements), `c`: super-diagonal (n-1 elements).
pub(crate) fn solve_tridiagonal_system(a: &[f64], b: &[f64], c: &[f64], d: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    let mut c_prime = vec![0.0; n];
    let mut x = vec![0.0; n];
    c_prime[0] = c[0] / b[0];
    d[0] /= b[0];
    for i in 1..n {
        let m = 1.0 / (b[i] - a[i - 1] * c_prime[i - 1]);
        c_prime[i] = if i < n - 1 { c[i] * m } else { 0.0 };
        d[i] = (d[i] - a[i - 1] * d[i - 1]) * m;
    }
    x[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d[i] - c_prime[i] * x[i + 1];
    }
    x
}
/// Solves the 1D heat equation u_t = D * u_xx using the Crank-Nicolson method.
///
/// # Arguments
/// * `initial_condition` - The initial temperature distribution.
/// * `dx` - The spatial step size.
/// * `dt` - The time step.
/// * `d_coeff` - The diffusion coefficient D.
/// * `steps` - The number of time steps to simulate.
///
/// # Returns
/// The final temperature distribution.
pub fn solve_heat_equation_1d_cn(
    initial_condition: &[f64],
    dx: f64,
    dt: f64,
    d_coeff: f64,
    steps: usize,
) -> Vec<f64> {
    let n = initial_condition.len();
    let mut u = initial_condition.to_vec();
    let alpha = d_coeff * dt / (2.0 * dx * dx);
    let mut a = vec![-alpha; n - 1];
    let mut b = vec![1.0 + 2.0 * alpha; n];
    let mut c = vec![-alpha; n - 1];
    b[0] = 1.0;
    c[0] = 0.0;
    b[n - 1] = 1.0;
    a[n - 2] = 0.0;
    let mut d = vec![0.0; n];
    for _ in 0..steps {
        for i in 1..n - 1 {
            d[i] = alpha * u[i - 1] + (1.0 - 2.0 * alpha) * u[i] + alpha * u[i + 1];
        }
        d[0] = 0.0;
        d[n - 1] = 0.0;
        u = solve_tridiagonal_system(&a, &b, &c, &mut d);
    }
    u
}
/// Scenario for the 1D Crank-Nicolson solver.
pub fn simulate_1d_heat_conduction_cn_scenario() -> Vec<f64> {
    const N: usize = 100;
    const L: f64 = 1.0;
    let dx = L / (N - 1) as f64;
    let dt = 0.001;
    let d_coeff = 0.01;
    let mut u0 = vec![0.0; N];
    for (i, vars) in u0.iter_mut().enumerate().take(N) {
        *vars = (std::f64::consts::PI * i as f64 * dx).sin();
    }
    solve_heat_equation_1d_cn(&u0, dx, dt, d_coeff, 50)
}
pub struct HeatEquationSolverConfig {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub dt: f64,
    pub d_coeff: f64,
    pub steps: usize,
}
/// Solves the 2D heat equation u_t = D * (u_xx + u_yy) using the ADI method.
/// ADI splits the problem into two half-steps, each solving one dimension implicitly.
pub fn solve_heat_equation_2d_cn_adi(
    initial_condition: &[f64],
    config: &HeatEquationSolverConfig,
) -> Vec<f64> {
    let mut u = initial_condition.to_vec();
    let alpha_x = config.d_coeff * config.dt / (2.0 * config.dx * config.dx);
    let alpha_y = config.d_coeff * config.dt / (2.0 * config.dy * config.dy);
    let mut ax = vec![-alpha_x; config.nx - 1];
    let mut bx = vec![1.0 + 2.0 * alpha_x; config.nx];
    let mut cx = vec![-alpha_x; config.nx - 1];
    bx[0] = 1.0;
    cx[0] = 0.0;
    bx[config.nx - 1] = 1.0;
    ax[config.nx - 2] = 0.0;
    let mut ay = vec![-alpha_y; config.ny - 1];
    let mut by = vec![1.0 + 2.0 * alpha_y; config.ny];
    let mut cy = vec![-alpha_y; config.ny - 1];
    by[0] = 1.0;
    cy[0] = 0.0;
    by[config.ny - 1] = 1.0;
    ay[config.ny - 2] = 0.0;
    let mut u_half = vec![0.0; config.nx * config.ny];
    for _ in 0..config.steps {
        let mut d = vec![0.0; config.nx];
        for j in 1..config.ny - 1 {
            for i in 1..config.nx - 1 {
                let u_ij = u[j * config.nx + i];
                let u_ijm1 = u[(j - 1) * config.nx + i];
                let u_ijp1 = u[(j + 1) * config.nx + i];
                d[i] = alpha_y * u_ijm1 + (1.0 - 2.0 * alpha_y) * u_ij + alpha_y * u_ijp1;
            }
            let row_sol = solve_tridiagonal_system(&ax, &bx, &cx, &mut d);
            for i in 0..config.nx {
                u_half[j * config.nx + i] = row_sol[i];
            }
        }
        let mut d_transposed = vec![0.0; config.ny];
        for i in 1..config.nx - 1 {
            for j in 1..config.ny - 1 {
                let u_im1j = u_half[j * config.nx + i - 1];
                let u_ij = u_half[j * config.nx + i];
                let u_ip1j = u_half[j * config.nx + i + 1];
                d_transposed[j] =
                    alpha_x * u_im1j + (1.0 - 2.0 * alpha_x) * u_ij + alpha_x * u_ip1j;
            }
            let col_sol = solve_tridiagonal_system(&ay, &by, &cy, &mut d_transposed);
            for j in 0..config.ny {
                u[j * config.nx + i] = col_sol[j];
            }
        }
    }
    u
}
/// Scenario for the 2D Crank-Nicolson ADI solver.
pub fn simulate_2d_heat_conduction_cn_adi_scenario() -> Vec<f64> {
    const NX: usize = 50;
    const NY: usize = 50;
    let dx = 1.0 / (NX - 1) as f64;
    let dy = 1.0 / (NY - 1) as f64;
    let dt = 0.01;
    let d_coeff = 0.05;
    let config = HeatEquationSolverConfig {
        nx: NX,
        ny: NY,
        dx,
        dy,
        dt,
        d_coeff,
        steps: 50,
    };
    let mut u0 = vec![0.0; NX * NY];
    for j in 0..NY {
        for i in 0..NX {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            if (x - 0.5).powi(2) + (y - 0.5).powi(2) < 0.05 {
                u0[j * NX + i] = 100.0;
            }
        }
    }
    solve_heat_equation_2d_cn_adi(&u0, &config)
}
