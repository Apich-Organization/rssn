//! Space-dimension coverage of the numeric solvers: every added 1D/2D/3D
//! companion is checked against an analytic solution, a convergence order, a
//! conservation law, or the legacy routine of another dimension.

use std::f64::consts::PI;

use num_complex::Complex;
use rssn::kernels::matrix::Matrix;
use rssn::sim::classical;
use rssn::sim::dimensional::IterativeMethod;
use rssn::sim::dimensional::Transport;
use rssn::sim::dimensional::solve_heat_adi_nd;
use rssn::sim::dimensional::solve_poisson_iterative_nd;
use rssn::sim::dimensional::solve_schrodinger_cn_nd;
use rssn::sim::dimensional::solve_transport_nd;
use rssn::sim::dimensional::solve_wave_leapfrog_nd;
use rssn::sim::physics_cfd as cfd;
use rssn::sim::physics_cnm as cnm;
use rssn::sim::physics_fdm as fdm;
use rssn::sim::physics_fea as fea;
use rssn::sim::physics_fvm as fvm;
use rssn::sim::physics_mtm as mtm;

fn max_abs_diff(
    a: &[f64],
    b: &[f64],
) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Flat index with `x` fastest.
fn idx3(
    n: [usize; 3],
    i: usize,
    j: usize,
    k: usize,
) -> usize {
    i + n[0] * (j + n[1] * k)
}

/// Samples `f(x, y, z)` on a unit-box grid of `n` nodes per axis.
fn sample3(
    n: [usize; 3],
    h: [f64; 3],
    f: impl Fn(f64, f64, f64) -> f64,
) -> Vec<f64> {
    let mut v = vec![0.0; n[0] * n[1] * n[2]];
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                v[idx3(n, i, j, k)] = f(i as f64 * h[0], j as f64 * h[1], k as f64 * h[2]);
            }
        }
    }
    v
}

// ---------------------------------------------------------------------------
// Poisson
// ---------------------------------------------------------------------------

#[test]
fn poisson_1d_is_exact_for_a_quadratic() {
    // u'' = 2, u(0) = u(1) = 0  =>  u = x^2 - x (reproduced exactly by the stencil)
    let n = 17;
    let dx = 1.0 / (n - 1) as f64;
    let f = vec![2.0; n];
    let u0 = vec![0.0; n];
    for (name, u) in [
        (
            "sor",
            cfd::solve_poisson_1d_sor(&f, &u0, dx, 1.7, 5000, 1e-14).unwrap(),
        ),
        (
            "gs",
            cfd::solve_poisson_1d_gauss_seidel(&f, &u0, dx, 20000, 1e-14).unwrap(),
        ),
        (
            "jacobi",
            cfd::solve_poisson_1d_jacobi(&f, &u0, dx, 40000, 1e-14).unwrap(),
        ),
    ] {
        for (i, v) in u.iter().enumerate() {
            let x = i as f64 * dx;
            assert!((v - (x * x - x)).abs() < 1e-9, "{name} at {i}");
        }
    }
}

#[test]
fn poisson_2d_generic_matches_legacy_matrix_solvers() {
    let (nx, ny) = (9, 7);
    let (dx, dy) = (0.125, 0.2);
    let f_fn = |x: f64, y: f64| (PI * x).sin() * (2.0 * y).cos() + x * y;
    let mut flat_f = vec![0.0; nx * ny];
    let mut flat_u0 = vec![0.0; nx * ny];
    let mut f = Matrix::zeros(nx, ny);
    let mut u0 = Matrix::zeros(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = (i as f64 * dx, j as f64 * dy);
            flat_f[i + nx * j] = f_fn(x, y);
            *f.get_mut(i, j) = f_fn(x, y);
            let boundary = i == 0 || j == 0 || i == nx - 1 || j == ny - 1;
            let b = if boundary {
                1.0 + x - 0.5 * y
            } else {
                0.0
            };
            flat_u0[i + nx * j] = b;
            *u0.get_mut(i, j) = b;
        }
    }
    let cases: [(IterativeMethod, Matrix<f64>); 3] = [
        (
            IterativeMethod::Jacobi,
            cfd::solve_poisson_2d_jacobi(&f, &u0, dx, dy, 400, 1e-30),
        ),
        (
            IterativeMethod::GaussSeidel,
            cfd::solve_poisson_2d_gauss_seidel(&f, &u0, dx, dy, 400, 1e-30),
        ),
        (
            IterativeMethod::Sor(1.5),
            cfd::solve_poisson_2d_sor(&f, &u0, dx, dy, 1.5, 400, 1e-30),
        ),
    ];
    for (method, legacy) in cases {
        let u =
            solve_poisson_iterative_nd(&flat_f, &flat_u0, &[nx, ny], &[dx, dy], method, 400, 1e-30)
                .unwrap();
        for j in 0..ny {
            for i in 0..nx {
                assert!(
                    (u[i + nx * j] - legacy.get(i, j)).abs() < 1e-11,
                    "{method:?}"
                );
            }
        }
    }
}

fn poisson_3d_error(
    n: usize,
    omega: f64,
) -> f64 {
    let h = 1.0 / (n - 1) as f64;
    let dims = [n; 3];
    let exact = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let f = sample3(dims, [h; 3], |x, y, z| -3.0 * PI * PI * exact(x, y, z));
    let u0 = vec![0.0; n * n * n];
    let u = cfd::solve_poisson_3d_sor(&f, &u0, dims, [h; 3], omega, 4000, 1e-13).unwrap();
    max_abs_diff(&u, &sample3(dims, [h; 3], exact))
}

#[test]
fn poisson_3d_sor_converges_with_second_order() {
    let e1 = poisson_3d_error(9, 1.6);
    let e2 = poisson_3d_error(17, 1.8);
    assert!(e2 < 6e-3, "fine error {e2}");
    let order = (e1 / e2).log2();
    assert!((1.7..2.4).contains(&order), "observed order {order}");
}

#[test]
fn poisson_3d_methods_agree_and_reproduce_the_1d_solution() {
    // f = 2 with Dirichlet data x^2 - x on every boundary node: the 3D solution
    // is independent of y and z and equals the 1D result at every line.
    let n = [9, 5, 6];
    let h = [0.125, 0.25, 0.2];
    let g = |x: f64| x * x - x;
    let f = vec![2.0; n[0] * n[1] * n[2]];
    let mut u0 = vec![0.0; f.len()];
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let boundary =
                    i == 0 || j == 0 || k == 0 || i == n[0] - 1 || j == n[1] - 1 || k == n[2] - 1;
                if boundary {
                    u0[idx3(n, i, j, k)] = g(i as f64 * h[0]);
                }
            }
        }
    }
    let sor = cfd::solve_poisson_3d_sor(&f, &u0, n, h, 1.6, 4000, 1e-14).unwrap();
    let gs = cfd::solve_poisson_3d_gauss_seidel(&f, &u0, n, h, 8000, 1e-14).unwrap();
    let jac = cfd::solve_poisson_3d_jacobi(&f, &u0, n, h, 20000, 1e-14).unwrap();
    let exact = sample3(n, h, |x, _, _| g(x));
    assert!(max_abs_diff(&sor, &exact) < 1e-8);
    assert!(max_abs_diff(&gs, &exact) < 1e-8);
    assert!(max_abs_diff(&jac, &exact) < 1e-8);
    // the FDM zero-Dirichlet front end: -3 pi^2 mode with the same machinery
    let h2 = 1.0 / 8.0;
    let src = sample3([9; 3], [h2; 3], |x, y, z| {
        -3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin()
    });
    let u = fdm::solve_poisson_3d(&src, [9; 3], [h2; 3], 1.5, 3000, 1e-13).unwrap();
    let mid = u[idx3([9; 3], 4, 4, 4)];
    assert!((mid - 1.0).abs() < 2.5e-2, "centre value {mid}");
}

#[test]
fn fdm_poisson_1d_and_nd_front_ends() {
    let n = 33;
    let dx = 1.0 / (n - 1) as f64;
    // u'' = -pi^2 sin(pi x)  =>  u = sin(pi x)
    let src: Vec<f64> = (0..n)
        .map(|i| -PI * PI * (PI * i as f64 * dx).sin())
        .collect();
    let u = fdm::solve_poisson_1d(&src, dx, 1.8, 5000, 1e-13).unwrap();
    for (i, v) in u.iter().enumerate() {
        assert!((v - (PI * i as f64 * dx).sin()).abs() < 2e-3);
    }
    // the n-D entry point with a single axis is the same solver
    let nd = fdm::solve_poisson_nd(
        &fdm::PoissonSolverConfigNd {
            dims: vec![n],
            spacing: vec![dx],
            omega: 1.8,
            max_iter: 5000,
            tolerance: 1e-13,
        },
        &src,
    )
    .unwrap();
    assert!(max_abs_diff(&u, &nd) < 1e-15);
    // 2D agrees with the legacy red-black solver to its tolerance
    let (w, hgt) = (11, 9);
    let (dx2, dy2) = (0.1, 0.125);
    let mut source = fdm::FdmGrid::<f64>::new(fdm::Dimensions::D2(w, hgt));
    for j in 0..hgt {
        for i in 0..w {
            source.as_mut_slice()[j * w + i] = ((i as f64) * 0.3).sin() + (j as f64) * 0.1;
        }
    }
    let legacy = fdm::solve_poisson_2d(
        &fdm::PoissonSolverConfig2D {
            width: w,
            height: hgt,
            dx: dx2,
            dy: dy2,
            omega: 1.5,
            max_iter: 5000,
            tolerance: 1e-13,
        },
        &source,
    );
    let generic = fdm::solve_poisson_nd(
        &fdm::PoissonSolverConfigNd {
            dims: vec![w, hgt],
            spacing: vec![dx2, dy2],
            omega: 1.5,
            max_iter: 5000,
            tolerance: 1e-13,
        },
        source.as_slice(),
    )
    .unwrap();
    assert!(max_abs_diff(legacy.as_slice(), &generic) < 1e-9);
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

fn bump(x: f64) -> f64 {
    (-80.0 * (x - 0.4) * (x - 0.4)).exp()
}

/// x-only data on a 3D box, as a flat array and as a 1D line.
fn line_and_box(n: [usize; 3]) -> (Vec<f64>, Vec<f64>) {
    let h = 1.0 / (n[0] - 1) as f64;
    let line: Vec<f64> = (0..n[0]).map(|i| bump(i as f64 * h)).collect();
    let boxed = sample3(n, [h; 3], |x, _, _| bump(x));
    (line, boxed)
}

fn assert_box_lines_equal(
    boxed: &[f64],
    n: [usize; 3],
    line: &[f64],
    tol: f64,
    what: &str,
) {
    // Boundary lines are frozen at their initial values (not evolved), so their
    // influence reaches the centre line only after `n/2` steps; the tests use fewer.
    let (j, k) = (n[1] / 2, n[2] / 2);
    for i in 0..n[0] {
        let d = (boxed[idx3(n, i, j, k)] - line[i]).abs();
        assert!(d < tol, "{what}: centre line differs by {d} at {i}");
    }
}

#[test]
fn transport_one_axis_reproduces_the_legacy_1d_solvers() {
    let n = 41;
    let h = 1.0 / (n - 1) as f64;
    let u0: Vec<f64> = (0..n).map(|i| bump(i as f64 * h)).collect();
    let dt = 0.2 * h * h / 0.5;
    let steps = 60;
    let cases: [(&str, Transport, Vec<Vec<f64>>); 4] = [
        (
            "advection",
            Transport::Advection { velocity: vec![0.7] },
            cfd::solve_advection_1d(&u0, 0.7, h, 0.5 * h, 25),
        ),
        (
            "diffusion",
            Transport::Diffusion { alpha: 0.5 },
            cfd::solve_diffusion_1d(&u0, 0.5, h, dt, steps),
        ),
        (
            "advection-diffusion",
            Transport::AdvectionDiffusion {
                velocity: vec![0.3],
                alpha: 0.2,
            },
            cfd::solve_advection_diffusion_1d(&u0, 0.3, 0.2, h, dt, steps),
        ),
        (
            "burgers",
            Transport::Burgers { nu: 0.05 },
            cfd::solve_burgers_1d(&u0, 0.05, h, dt, steps),
        ),
    ];
    for (name, eq, legacy) in cases {
        let (dt_used, steps_used) = if name == "advection" {
            (0.5 * h, 25)
        } else {
            (dt, steps)
        };
        let all = solve_transport_nd(&u0, &[n], &[h], &eq, dt_used, steps_used).unwrap();
        assert_eq!(all.len(), legacy.len(), "{name}");
        for (a, b) in all.iter().zip(&legacy) {
            assert!(max_abs_diff(a, b) < 1e-12, "{name}");
        }
    }
}

#[test]
fn transport_3d_with_x_only_data_reproduces_the_1d_result_on_every_line() {
    let n = [33, 25, 25];
    let h = 1.0 / 32.0;
    let (line, boxed) = line_and_box(n);
    let dt = 0.1 * h * h;
    let steps = 10;
    let hh = [h; 3];
    let d1 = cfd::solve_diffusion_1d(&line, 0.5, h, dt, steps)
        .pop()
        .unwrap();
    let d3 = cfd::solve_diffusion_3d(&boxed, n, 0.5, hh, dt, steps)
        .unwrap()
        .pop()
        .unwrap();
    assert_box_lines_equal(&d3, n, &d1, 1e-12, "diffusion");

    let a1 = cfd::solve_advection_1d(&line, 0.8, h, 0.5 * h, 10)
        .pop()
        .unwrap();
    let a3 = cfd::solve_advection_3d(&boxed, n, [0.8, 0.0, 0.0], hh, 0.5 * h, 10)
        .unwrap()
        .pop()
        .unwrap();
    assert_box_lines_equal(&a3, n, &a1, 1e-12, "advection");

    let ad1 = cfd::solve_advection_diffusion_1d(&line, 0.3, 0.2, h, dt, steps)
        .pop()
        .unwrap();
    let ad3 = cfd::solve_advection_diffusion_3d(&boxed, n, [0.3, 0.0, 0.0], 0.2, hh, dt, steps)
        .unwrap()
        .pop()
        .unwrap();
    assert_box_lines_equal(&ad3, n, &ad1, 1e-12, "advection-diffusion");

    let b1 = cfd::solve_burgers_1d(&line, 0.05, h, dt, steps)
        .pop()
        .unwrap();
    let b3 = cfd::solve_burgers_3d(&boxed, n, 0.05, hh, dt, steps)
        .unwrap()
        .pop()
        .unwrap();
    assert_box_lines_equal(&b3, n, &b1, 1e-12, "burgers");

    // the FDM flavours (central Burgers, final state only)
    let fb1 = fdm::solve_burgers_1d(&line, h, 0.05, dt, steps);
    let fb3 = fdm::solve_burgers_3d(&boxed, n, hh, 0.05, dt, steps).unwrap();
    assert_box_lines_equal(&fb3, n, &fb1, 1e-12, "fdm burgers");
    let fa1 = fdm::solve_advection_diffusion_1d(&line, h, 0.3, 0.2, dt, steps);
    let fa3 =
        fdm::solve_advection_diffusion_3d(&boxed, n, hh, [0.3, 0.0, 0.0], 0.2, dt, steps).unwrap();
    assert_box_lines_equal(&fa3, n, &fa1, 1e-12, "fdm advection-diffusion");
}

#[test]
fn diffusion_3d_mode_decays_at_the_analytic_rate() {
    let n = [17; 3];
    let h = 1.0 / 16.0;
    let mode = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let u0 = sample3(n, [h; 3], mode);
    let dt = h * h / 6.5;
    let steps = 36;
    let t = dt * steps as f64;
    let out = cfd::solve_diffusion_3d(&u0, n, 1.0, [h; 3], dt, steps)
        .unwrap()
        .pop()
        .unwrap();
    let decay = (-3.0 * PI * PI * t).exp();
    let exact: Vec<f64> = u0.iter().map(|v| v * decay).collect();
    assert!(max_abs_diff(&out, &exact) < 5e-3);
}

#[test]
fn periodic_advection_is_an_exact_shift_at_cfl_one_in_2d_and_conserves_mass() {
    // period of the identification is n - 2 = 8 nodes along x
    let dims = [10, 6];
    let h = [0.1, 0.1];
    let mut u0 = vec![0.0; 60];
    for j in 0..6 {
        for i in 0..10 {
            // 8-periodic in x (node 0 == node 8, node 9 == node 1), constant in y
            let phase = (i + 8 - 1) % 8; // index within the period, starting at node 1
            u0[i + 10 * j] = [0.0, 1.0, 3.0, 2.0, 0.5, 0.0, 0.0, 0.25][phase];
        }
    }
    let all = cfd::solve_advection_2d(&u0, dims, [1.0, 0.0], h, 0.1, 8).unwrap();
    assert!(
        max_abs_diff(&all[8], &all[0]) < 1e-14,
        "a full period returns the data"
    );
    let one = &all[1];
    for i in 1..9 {
        assert!(
            (one[i + 10 * 2] - all[0][i - 1 + 10 * 2]).abs() < 1e-14,
            "shift at {i}"
        );
    }
    let mass = |u: &[f64]| (1..9).map(|i| u[i + 10 * 2]).sum::<f64>();
    assert!((mass(&all[3]) - mass(&all[0])).abs() < 1e-13);
}

#[test]
fn transport_is_symmetric_under_swapping_the_axes() {
    let (nx, ny) = (13, 9);
    let (hx, hy) = (0.1, 0.15);
    let f = |x: f64, y: f64| (-(((x - 0.5).powi(2)) / 0.05 + ((y - 0.6).powi(2)) / 0.08)).exp();
    let a: Vec<f64> = (0..nx * ny)
        .map(|p| f((p % nx) as f64 * hx, (p / nx) as f64 * hy))
        .collect();
    // transposed copy: nx <-> ny, hx <-> hy
    let b: Vec<f64> = (0..nx * ny)
        .map(|p| f((p / ny) as f64 * hx, (p % ny) as f64 * hy))
        .collect();
    let dt = 0.01 * 0.1;
    let r1 = cfd::solve_advection_diffusion_2d(&a, [nx, ny], [0.4, -0.2], 0.05, [hx, hy], dt, 20)
        .unwrap()
        .pop()
        .unwrap();
    let r2 = cfd::solve_advection_diffusion_2d(&b, [ny, nx], [-0.2, 0.4], 0.05, [hy, hx], dt, 20)
        .unwrap()
        .pop()
        .unwrap();
    for j in 0..ny {
        for i in 0..nx {
            assert!((r1[i + nx * j] - r2[j + ny * i]).abs() < 1e-13);
        }
    }
}

#[test]
fn transport_rejects_bad_input() {
    let u = vec![0.0; 25];
    assert!(
        solve_transport_nd(
            &u,
            &[5, 5],
            &[0.1],
            &Transport::Diffusion { alpha: 1.0 },
            0.1,
            1
        )
        .is_err()
    );
    assert!(
        solve_transport_nd(
            &u,
            &[5, 5],
            &[0.1, 0.1],
            &Transport::Advection { velocity: vec![1.0] },
            0.1,
            1
        )
        .is_err()
    );
    assert!(
        solve_transport_nd(
            &u,
            &[5, 4],
            &[0.1, 0.1],
            &Transport::Diffusion { alpha: 1.0 },
            0.1,
            1
        )
        .is_err()
    );
    assert!(
        solve_transport_nd(
            &u,
            &[5, 5],
            &[0.1, 0.1],
            &Transport::Diffusion { alpha: 1.0 },
            0.0,
            1
        )
        .is_err()
    );
}

// ---------------------------------------------------------------------------
// Crank-Nicolson heat (ADI)
// ---------------------------------------------------------------------------

#[test]
fn heat_adi_one_axis_equals_the_legacy_1d_crank_nicolson() {
    let n = 41;
    let dx = 1.0 / (n - 1) as f64;
    let mut u0: Vec<f64> = (0..n).map(|i| bump(i as f64 * dx)).collect();
    u0[0] = 0.0;
    u0[n - 1] = 0.0;
    let legacy = cnm::solve_heat_equation_1d_cn_tridiagonal(&u0, dx, 0.01, 0.3, 25);
    let generic = solve_heat_adi_nd(&u0, &[n], &[dx], 0.3, 0.01, 25).unwrap();
    assert!(max_abs_diff(&legacy, &generic) < 1e-13);
}

#[test]
fn heat_adi_3d_of_a_product_is_the_product_of_the_1d_solutions() {
    let n = [21, 9, 11];
    let h = [0.05, 0.125, 0.1];
    let (dt, steps, d) = (0.004, 12, 0.7);
    let fx: Vec<f64> = (0..n[0]).map(|i| bump(i as f64 * h[0])).collect();
    let fy: Vec<f64> = (0..n[1]).map(|j| (PI * j as f64 * h[1]).sin()).collect();
    let fz: Vec<f64> = (0..n[2])
        .map(|k| (2.0 * PI * k as f64 * h[2]).sin())
        .collect();
    let mut fx0 = fx.clone();
    fx0[0] = 0.0;
    fx0[n[0] - 1] = 0.0;
    let mut u0 = vec![0.0; n[0] * n[1] * n[2]];
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                u0[idx3(n, i, j, k)] = fx0[i] * fy[j] * fz[k];
            }
        }
    }
    let ux = cnm::solve_heat_equation_1d_cn_tridiagonal(&fx0, h[0], dt, d, steps);
    let uy = cnm::solve_heat_equation_1d_cn_tridiagonal(&fy, h[1], dt, d, steps);
    let uz = cnm::solve_heat_equation_1d_cn_tridiagonal(&fz, h[2], dt, d, steps);
    let cfg = cnm::HeatEquationSolverConfig3D {
        nx: n[0],
        ny: n[1],
        nz: n[2],
        dx: h[0],
        dy: h[1],
        dz: h[2],
        dt,
        d_coeff: d,
        steps,
    };
    let u = cnm::solve_heat_equation_3d_cn_adi(&u0, &cfg).unwrap();
    for k in 1..n[2] - 1 {
        for j in 1..n[1] - 1 {
            for i in 1..n[0] - 1 {
                let want = ux[i] * uy[j] * uz[k];
                assert!((u[idx3(n, i, j, k)] - want).abs() < 1e-12);
            }
        }
    }
}

fn heat_adi_3d_error(
    n: usize,
    dt: f64,
    steps: usize,
) -> f64 {
    let h = 1.0 / (n - 1) as f64;
    let mode = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let u0 = sample3([n; 3], [h; 3], mode);
    let u = solve_heat_adi_nd(&u0, &[n; 3], &[h; 3], 1.0, dt, steps).unwrap();
    let decay = (-3.0 * PI * PI * dt * steps as f64).exp();
    let exact: Vec<f64> = u0.iter().map(|v| v * decay).collect();
    max_abs_diff(&u, &exact)
}

#[test]
fn heat_adi_3d_is_second_order_in_space_and_time() {
    let coarse = heat_adi_3d_error(9, 0.02, 4);
    let fine = heat_adi_3d_error(17, 0.01, 8);
    assert!(fine < 3e-3, "fine error {fine}");
    let order = (coarse / fine).log2();
    assert!((1.6..2.5).contains(&order), "observed order {order}");
}

#[test]
fn classical_heat_cn_nd_wrappers_match_the_1d_routine_and_the_exact_decay() {
    let f = |x: f64| (PI * x).sin();
    let legacy =
        classical::solve_heat_equation_1d_crank_nicolson(&f, 0.4, (0.0, 1.0), 33, (0.0, 0.1), 20)
            .unwrap();
    let generic = classical::solve_heat_equation_nd_crank_nicolson(
        &|x: &[f64]| f(x[0]),
        0.4,
        &[(0.0, 1.0)],
        &[33],
        (0.0, 0.1),
        20,
    )
    .unwrap();
    assert_eq!(legacy.len(), generic.len());
    for (a, b) in legacy.iter().zip(&generic) {
        assert!(max_abs_diff(a, b) < 1e-12);
    }
    // 2D and 3D: exact eigenmodes
    let two = classical::solve_heat_equation_2d_crank_nicolson(
        &|x, y| (PI * x).sin() * (PI * y).sin(),
        0.4,
        [(0.0, 1.0), (0.0, 1.0)],
        [25, 25],
        (0.0, 0.1),
        20,
    )
    .unwrap();
    let three = classical::solve_heat_equation_3d_crank_nicolson(
        &|x, y, z| (PI * x).sin() * (PI * y).sin() * (PI * z).sin(),
        0.4,
        [(0.0, 1.0); 3],
        [13; 3],
        (0.0, 0.1),
        20,
    )
    .unwrap();
    assert_eq!(two.len(), 21);
    assert_eq!(three.len(), 21);
    let c2 = two[20][12 + 25 * 12];
    assert!(
        (c2 - (-2.0 * PI * PI * 0.4 * 0.1).exp()).abs() < 3e-3,
        "2D centre {c2}"
    );
    let c3 = three[20][idx3([13; 3], 6, 6, 6)];
    assert!(
        (c3 - (-3.0 * PI * PI * 0.4 * 0.1).exp()).abs() < 1.5e-2,
        "3D centre {c3}"
    );
    assert!(
        classical::solve_heat_equation_nd_crank_nicolson(
            &|_: &[f64]| 0.0,
            1.0,
            &[(0.0, 1.0)],
            &[2],
            (0.0, 1.0),
            1
        )
        .is_err()
    );
}

// ---------------------------------------------------------------------------
// Schrodinger
// ---------------------------------------------------------------------------

fn packet(
    x: f64,
    x0: f64,
    k: f64,
) -> Complex<f64> {
    Complex::from_polar((-(x - x0) * (x - x0) / 0.5).exp(), k * x)
}

#[test]
fn schrodinger_cn_one_axis_free_particle_equals_the_legacy_solver() {
    let n = 129;
    let dx = 12.0 / (n - 1) as f64;
    let psi0: Vec<Complex<f64>> = (0..n)
        .map(|i| packet(-6.0 + i as f64 * dx, -1.0, 2.0))
        .collect();
    let mut p0 = psi0.clone();
    p0[0] = Complex::new(0.0, 0.0);
    p0[n - 1] = Complex::new(0.0, 0.0);
    let v = vec![0.0; n];
    let legacy = cnm::solve_schrodinger_1d_cn(&p0, &v, dx, 0.01, 40);
    let generic = solve_schrodinger_cn_nd(&p0, &v, &[n], &[dx], 0.01, 40).unwrap();
    let d = legacy
        .iter()
        .zip(&generic)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max);
    assert!(d < 1e-12, "difference {d}");
}

#[test]
fn schrodinger_2d_free_separable_data_is_the_product_of_1d_evolutions() {
    let (nx, ny) = (49, 33);
    let (hx, hy) = (12.0 / (nx - 1) as f64, 8.0 / (ny - 1) as f64);
    let fx: Vec<Complex<f64>> = (0..nx)
        .map(|i| packet(-6.0 + i as f64 * hx, -1.0, 1.5))
        .collect();
    let fy: Vec<Complex<f64>> = (0..ny)
        .map(|j| packet(-4.0 + j as f64 * hy, 0.5, -1.0))
        .collect();
    let z = |v: &mut Vec<Complex<f64>>| {
        let n = v.len();
        v[0] = Complex::new(0.0, 0.0);
        v[n - 1] = Complex::new(0.0, 0.0);
    };
    let (mut fx, mut fy) = (fx, fy);
    z(&mut fx);
    z(&mut fy);
    let mut psi0 = vec![Complex::new(0.0, 0.0); nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            psi0[i + nx * j] = fx[i] * fy[j];
        }
    }
    let (dt, steps) = (0.02, 15);
    let ex = cnm::solve_schrodinger_1d_cn(&fx, &vec![0.0; nx], hx, dt, steps);
    let ey = cnm::solve_schrodinger_1d_cn(&fy, &vec![0.0; ny], hy, dt, steps);
    let out =
        cnm::solve_schrodinger_2d_cn(&psi0, &vec![0.0; nx * ny], [nx, ny], [hx, hy], dt, steps)
            .unwrap();
    let mut worst = 0.0_f64;
    for j in 0..ny {
        for i in 0..nx {
            worst = worst.max((out[i + nx * j] - ex[i] * ey[j]).norm());
        }
    }
    assert!(worst < 1e-12, "difference {worst}");
}

#[test]
fn schrodinger_2d_harmonic_ground_state_is_stationary_up_to_a_phase() {
    let n = 65;
    let h = 12.0 / (n - 1) as f64;
    let g = |x: f64, y: f64| (-(x * x + y * y) / 2.0).exp() / PI.sqrt();
    let mut psi0 = vec![Complex::new(0.0, 0.0); n * n];
    let mut v = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (-6.0 + i as f64 * h, -6.0 + j as f64 * h);
            psi0[i + n * j] = Complex::new(g(x, y), 0.0);
            v[i + n * j] = 0.5 * (x * x + y * y);
        }
    }
    let (dt, steps) = (0.01, 60);
    let out = cnm::solve_schrodinger_2d_cn(&psi0, &v, [n, n], [h, h], dt, steps).unwrap();
    let phase = Complex::from_polar(1.0, -dt * steps as f64); // E0 = 1 for the 2D oscillator
    let err = psi0
        .iter()
        .zip(&out)
        .map(|(a, b)| (a * phase - b).norm())
        .fold(0.0, f64::max);
    assert!(err < 2e-2 * g(0.0, 0.0), "deviation {err}");
    let norm = |p: &[Complex<f64>]| p.iter().map(|c| c.norm_sqr()).sum::<f64>() * h * h;
    assert!((norm(&out) - norm(&psi0)).abs() < 1e-12, "norm drift");
}

#[test]
fn schrodinger_3d_box_eigenstate_and_norm_conservation() {
    let n = 21;
    let h = 1.0 / (n - 1) as f64;
    let mode = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let real = sample3([n; 3], [h; 3], mode);
    let psi0: Vec<Complex<f64>> = real.iter().map(|&r| Complex::new(r, 0.0)).collect();
    let (dt, steps) = (0.002, 20);
    let out = cnm::solve_schrodinger_3d_cn(&psi0, &vec![0.0; n * n * n], [n; 3], [h; 3], dt, steps)
        .unwrap();
    // E = 3 pi^2 / 2 for V = 0
    let phase = Complex::from_polar(1.0, -1.5 * PI * PI * dt * steps as f64);
    let err = psi0
        .iter()
        .zip(&out)
        .map(|(a, b)| (a * phase - b).norm())
        .fold(0.0, f64::max);
    assert!(err < 2e-2, "deviation {err}");
    // a potential does not break unitarity
    let v: Vec<f64> = sample3([n; 3], [h; 3], |x, y, z| 50.0 * (x + 2.0 * y * y - z));
    let withv = cnm::solve_schrodinger_3d_cn(&psi0, &v, [n; 3], [h; 3], dt, steps).unwrap();
    let norm = |p: &[Complex<f64>]| p.iter().map(|c| c.norm_sqr()).sum::<f64>();
    assert!((norm(&withv) - norm(&psi0)).abs() < 1e-11 * norm(&psi0));
}

// ---------------------------------------------------------------------------
// Wave
// ---------------------------------------------------------------------------

#[test]
fn wave_one_axis_equals_the_legacy_routine() {
    let n = 51;
    let dx = 1.0 / (n - 1) as f64;
    let u0: Vec<f64> = (0..n)
        .map(|i| (PI * i as f64 * dx).sin() + 0.2 * (3.0 * PI * i as f64 * dx).sin())
        .collect();
    let ut = vec![0.0; n];
    let legacy = classical::solve_wave_equation_1d(&u0, &ut, 1.3, dx, 0.5 * dx, 40).unwrap();
    let generic = solve_wave_leapfrog_nd(&u0, &ut, &[n], &[dx], 1.3, 0.5 * dx, 40).unwrap();
    assert_eq!(legacy.len(), generic.len());
    for (a, b) in legacy.iter().zip(&generic) {
        assert!(max_abs_diff(a, b) < 1e-12);
    }
}

#[test]
fn wave_3d_standing_mode_oscillates_at_the_analytic_frequency() {
    let n = [21, 17, 19];
    let h = [1.0 / 20.0, 1.0 / 16.0, 1.0 / 18.0];
    let c = 1.5;
    let mode = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let u0 = sample3(n, h, mode);
    let ut = vec![0.0; u0.len()];
    let dt = 0.5 / (c * (1.0 / (h[0] * h[0]) + 1.0 / (h[1] * h[1]) + 1.0 / (h[2] * h[2])).sqrt());
    let steps = 30;
    let snaps = classical::solve_wave_equation_3d(&u0, &ut, n, h, c, dt, steps).unwrap();
    assert_eq!(snaps.len(), steps + 1);
    let omega = c * PI * 3.0_f64.sqrt();
    let want: Vec<f64> = u0
        .iter()
        .map(|v| v * (omega * dt * steps as f64).cos())
        .collect();
    assert!(max_abs_diff(&snaps[steps], &want) < 1.5e-2);
    // CFL violation is reported
    assert!(classical::solve_wave_equation_3d(&u0, &ut, n, h, c, 10.0 * dt, 2).is_err());
}

#[test]
fn wave_2d_is_symmetric_under_swapping_the_axes_and_conserves_energy() {
    let (nx, ny) = (25, 17);
    let (hx, hy) = (0.04, 0.0625);
    let f = |x: f64, y: f64| (-(((x - 0.4).powi(2)) + ((y - 0.5).powi(2))) / 0.01).exp();
    let a: Vec<f64> = (0..nx * ny)
        .map(|p| f((p % nx) as f64 * hx, (p / nx) as f64 * hy))
        .collect();
    let b: Vec<f64> = (0..nx * ny)
        .map(|p| f((p / ny) as f64 * hx, (p % ny) as f64 * hy))
        .collect();
    let z = vec![0.0; nx * ny];
    let dt = 0.01;
    let ra = classical::solve_wave_equation_2d(&a, &z, [nx, ny], [hx, hy], 1.0, dt, 25).unwrap();
    let rb = classical::solve_wave_equation_2d(&b, &z, [ny, nx], [hy, hx], 1.0, dt, 25).unwrap();
    for j in 0..ny {
        for i in 0..nx {
            assert!((ra[25][i + nx * j] - rb[25][j + ny * i]).abs() < 1e-13);
        }
    }
    // exact discrete leapfrog invariant: E = |v|^2/2 - u^{n+1} . (L u^n) / 2 (c = 1)
    let energy = |prev: &[f64], cur: &[f64]| {
        let mut e = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let p = i + nx * j;
                let v = (cur[p] - prev[p]) / dt;
                let lap = (prev[p + 1] - 2.0 * prev[p] + prev[p - 1]) / (hx * hx)
                    + (prev[p + nx] - 2.0 * prev[p] + prev[p - nx]) / (hy * hy);
                e += 0.5 * (v * v - cur[p] * lap) * hx * hy;
            }
        }
        e
    };
    let e1 = energy(&ra[3], &ra[4]);
    let e2 = energy(&ra[20], &ra[21]);
    assert!((e1 - e2).abs() < 1e-9 * e1, "energy {e1} -> {e2}");
}

// ---------------------------------------------------------------------------
// Multigrid 3D
// ---------------------------------------------------------------------------

fn mg3_error(n: usize) -> f64 {
    let h = 1.0 / (n - 1) as f64;
    let exact = |x: f64, y: f64, z: f64| (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
    let f = sample3([n; 3], [h; 3], |x, y, z| 3.0 * PI * PI * exact(x, y, z));
    let u = mtm::solve_poisson_3d_multigrid(n, &f, 12).unwrap();
    max_abs_diff(&u, &sample3([n; 3], [h; 3], exact))
}

#[test]
fn multigrid_3d_is_second_order_and_validates_its_input() {
    let e1 = mg3_error(9);
    let e2 = mg3_error(17);
    let order = (e1 / e2).log2();
    assert!((1.7..2.4).contains(&order), "observed order {order}");
    assert!(e2 < 5e-3);
    assert!(mtm::solve_poisson_3d_multigrid(10, &vec![0.0; 1000], 2).is_err());
    assert!(mtm::solve_poisson_3d_multigrid(9, &[0.0; 4], 2).is_err());
}

// ---------------------------------------------------------------------------
// Finite volumes
// ---------------------------------------------------------------------------

#[test]
fn fvm_burgers_nd_matches_the_1d_mesh_solver_and_conserves_mass() {
    let n = 80;
    let dx = 1.0 / n as f64;
    let top = |x: f64| {
        if (0.2..0.4).contains(&x) {
            1.0
        } else {
            0.0
        }
    };
    let mut mesh = fvm::Mesh::new(n, 1.0, top);
    let u0: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();
    let (dt, steps) = (0.5 * dx, 50);
    let legacy = fvm::solve_burgers_1d(&mut mesh, dt, steps);
    let generic = fvm::solve_burgers_nd(&u0, &[n], &[dx], dt, steps).unwrap();
    assert!(max_abs_diff(&legacy, &generic) < 1e-12);

    // 2D with y-independent data equals the 1D line in every row
    let (nx, ny) = (n, 6);
    let mut flat = vec![0.0; nx * ny];
    for j in 0..ny {
        flat[j * nx..(j + 1) * nx].copy_from_slice(&u0);
    }
    let two = fvm::solve_burgers_2d(&flat, [nx, ny], [dx, 0.1], dt, steps).unwrap();
    // y-independent data stays y-independent (the viscosity is shared between the
    // axes, so the rows are not bit-equal to the 1D scheme, only to each other)
    for j in 1..ny {
        assert!(max_abs_diff(&two[j * nx..(j + 1) * nx], &two[..nx]) < 1e-13);
    }
    let sum_line: f64 = legacy.iter().sum();
    let sum_row: f64 = two[..nx].iter().sum();
    assert!((sum_row - 16.0).abs() < 1e-6, "mass of a 2D row {sum_row}");
    assert!(
        (sum_line - sum_row).abs() < 1e-3,
        "mass of a row {sum_row} vs 1D {sum_line}"
    );
    // 3D: a compact bump in the interior conserves mass (flux through the faces is zero)
    let m = 16;
    let h = 1.0 / m as f64;
    let blob = sample3([m; 3], [h; 3], |x, y, z| {
        (-(((x - 0.5).powi(2)) + ((y - 0.5).powi(2)) + ((z - 0.5).powi(2))) / 0.01).exp()
    });
    let evolved = fvm::solve_burgers_3d(&blob, [m; 3], [h; 3], 0.2 * h, 6).unwrap();
    let mass0: f64 = blob.iter().sum();
    let mass1: f64 = evolved.iter().sum();
    let edge = |u: &[f64]| {
        let mut e = 0.0_f64;
        for p in 0..u.len() {
            let (i, j, k) = (p % m, (p / m) % m, p / (m * m));
            if i == 0 || j == 0 || k == 0 || i == m - 1 || j == m - 1 || k == m - 1 {
                e = e.max(u[p].abs());
            }
        }
        e
    };
    assert!(
        // only the O(1e-4) outflow through the zero-gradient faces may leave
        (mass0 - mass1).abs() < 1e-7 * mass0,
        "mass {mass0} -> {mass1}; boundary max {} -> {}; min {}",
        edge(&blob),
        edge(&evolved),
        evolved.iter().copied().fold(f64::INFINITY, f64::min)
    );
    assert!(
        max_abs_diff(&blob, &evolved) > 1e-6,
        "the field must actually move"
    );
}

#[test]
fn fvm_shallow_water_2d_matches_1d_preserves_rest_and_swaps_axes() {
    let n = 60;
    let dx = 0.05;
    let h0: Vec<f64> = (0..n).map(|i| if i < n / 2 { 2.0 } else { 1.0 }).collect();
    let hu0 = vec![0.0; n];
    // unsplit Lax-Friedrichs: (|u| + sqrt(g h)) dt (1/dx + 1/dy) < 1
    let (dt, steps, g) = (0.004, 40, 9.81);
    let one = fvm::solve_shallow_water_1d(h0.clone(), hu0.clone(), dx, dt, steps, g);
    let ny = 5;
    let mut hh = vec![0.0; n * ny];
    for j in 0..ny {
        hh[j * n..(j + 1) * n].copy_from_slice(&h0);
    }
    let zero = vec![0.0; n * ny];
    let mass_1d: f64 = one.iter().map(|s| s.h).sum();
    let two =
        fvm::solve_shallow_water_2d(&hh, &zero, &zero, [n, ny], [dx, 0.1], dt, steps, g).unwrap();
    for j in 0..ny {
        for i in 0..n {
            let s = two[j * n + i];
            let r = two[i];
            assert!((s.h - r.h).abs() < 1e-13 && (s.hu - r.hu).abs() < 1e-13);
            assert!(s.hv.abs() < 1e-13);
        }
    }
    let mass_2d_row: f64 = two[..n].iter().map(|s| s.h).sum();
    // the waves have not reached the zero-gradient ends: only a small tail leaks
    assert!(
        (mass_2d_row - 90.0).abs() < 1e-5 * 90.0,
        "2D row mass {mass_2d_row}"
    );
    assert!(
        (mass_2d_row - mass_1d).abs() < 1e-4 * mass_1d,
        "{mass_2d_row} vs {mass_1d}"
    );
    // lake at rest stays at rest
    let flat = vec![1.5; 12 * 9];
    let zero = vec![0.0; 12 * 9];
    let rest =
        fvm::solve_shallow_water_2d(&flat, &zero, &zero, [12, 9], [0.1, 0.1], 0.01, 30, g).unwrap();
    assert!(
        rest.iter()
            .all(|s| (s.h - 1.5).abs() < 1e-14 && s.hu.abs() < 1e-14 && s.hv.abs() < 1e-14)
    );
    // axis swap: a droplet at (i, j) in an [nx, ny] box equals the swapped box with hu <-> hv
    let (nx, ny) = (14, 10);
    let drop = |x: usize, y: usize| {
        1.0 + (-(((x as f64 - 4.0).powi(2)) + ((y as f64 - 6.0).powi(2))) / 0.25).exp()
    };
    let a: Vec<f64> = (0..nx * ny).map(|p| drop(p % nx, p / nx)).collect();
    let b: Vec<f64> = (0..nx * ny).map(|p| drop(p / ny, p % ny)).collect();
    let z = vec![0.0; nx * ny];
    let ra = fvm::solve_shallow_water_2d(&a, &z, &z, [nx, ny], [0.1, 0.1], 0.005, 20, g).unwrap();
    let rb = fvm::solve_shallow_water_2d(&b, &z, &z, [ny, nx], [0.1, 0.1], 0.005, 20, g).unwrap();
    for j in 0..ny {
        for i in 0..nx {
            let (s, t) = (ra[i + nx * j], rb[j + ny * i]);
            assert!(
                (s.h - t.h).abs() < 1e-12,
                "h at ({i},{j}): {} vs {}",
                s.h,
                t.h
            );
            assert!(
                (s.hu - t.hv).abs() < 1e-12,
                "hu/hv at ({i},{j}): {} vs {}",
                s.hu,
                t.hv
            );
            assert!((s.hv - t.hu).abs() < 1e-12, "hv/hu at ({i},{j})");
        }
    }
    // mass is conserved in the closed box while the water is moving
    let m0: f64 = a.iter().sum();
    // two steps only: the wave has not reached the zero-gradient faces yet (the
    // third step already loses mass through them), so the total is conserved
    let early = fvm::solve_shallow_water_2d(&a, &z, &z, [nx, ny], [0.1, 0.1], 0.005, 2, g).unwrap();
    let m1: f64 = early.iter().map(|s| s.h).sum();
    assert!((m0 - m1).abs() < 1e-10 * m0, "mass {m0} -> {m1}");
    assert!(
        early.iter().any(|s| s.hu.abs() > 1e-6),
        "the water must be moving"
    );
}

// ---------------------------------------------------------------------------
// Finite elements
// ---------------------------------------------------------------------------

fn unit_tet(material: fea::Material) -> fea::TetrahedronElement3D {
    fea::TetrahedronElement3D::new(
        [0, 1, 2, 3],
        [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ],
        material,
    )
}

#[test]
fn tetrahedron_stiffness_is_symmetric_and_annihilates_rigid_motions() {
    let mat = fea::Material::steel();
    let tet = fea::TetrahedronElement3D::new(
        [3, 5, 8, 9],
        [
            (0.1, 0.0, 0.2),
            (1.2, 0.1, 0.0),
            (0.3, 0.9, 0.1),
            (0.2, 0.3, 1.1),
        ],
        mat,
    );
    assert!(tet.volume() > 0.0);
    let k = tet.local_stiffness_matrix();
    let scale = (0..12).map(|i| *k.get(i, i)).fold(0.0, f64::max);
    for i in 0..12 {
        for j in 0..12 {
            assert!((k.get(i, j) - k.get(j, i)).abs() < 1e-9 * scale);
        }
    }
    let apply = |u: &[f64]| -> f64 {
        (0..12)
            .map(|i| (0..12).map(|j| k.get(i, j) * u[j]).sum::<f64>().abs())
            .fold(0.0, f64::max)
    };
    // translations and an infinitesimal rotation u = w x r
    for axis in 0..3 {
        let mut u = [0.0; 12];
        for n in 0..4 {
            u[3 * n + axis] = 1.0;
        }
        assert!(apply(&u) < 1e-9 * scale);
    }
    let w = (0.3, -0.2, 0.5);
    let mut u = [0.0; 12];
    for (n, c) in tet.coords.iter().enumerate() {
        u[3 * n] = w.1 * c.2 - w.2 * c.1;
        u[3 * n + 1] = w.2 * c.0 - w.0 * c.2;
        u[3 * n + 2] = w.0 * c.1 - w.1 * c.0;
    }
    assert!(apply(&u) < 1e-9 * scale);
    assert_eq!(tet.dof_map()[0..6], [9, 10, 11, 15, 16, 17]);
}

#[test]
fn tetrahedron_reproduces_the_uniform_strain_energy() {
    let mat = fea::Material::new(1.0, 0.25, 1.0, 1.0, 0.0, 1.0);
    let tet = unit_tet(mat);
    assert!((tet.volume() - 1.0 / 6.0).abs() < 1e-15);
    let k = tet.local_stiffness_matrix();
    // u = (e x, 0, 0): energy density (lambda + 2 mu) e^2 / 2
    let e = 0.01;
    let mut u = [0.0; 12];
    for (n, c) in tet.coords.iter().enumerate() {
        u[3 * n] = e * c.0;
    }
    let energy: f64 = 0.5
        * (0..12)
            .map(|i| (0..12).map(|j| u[i] * k.get(i, j) * u[j]).sum::<f64>())
            .sum::<f64>();
    let lambda = 0.25 / (1.25 * 0.5);
    let mu = 1.0 / (2.0 * 1.25);
    let want = 0.5 * (lambda + 2.0 * mu) * e * e * tet.volume();
    assert!((energy - want).abs() < 1e-14, "{energy} vs {want}");
}

#[test]
fn assembled_unit_cube_gives_the_uniaxial_energy_and_generic_assembly_matches_legacy() {
    // Kuhn subdivision of the unit cube into 6 tetrahedra
    let mat = fea::Material::new(2.0, 0.3, 1.0, 1.0, 0.0, 1.0);
    let corner = |n: usize| ((n & 1) as f64, ((n >> 1) & 1) as f64, ((n >> 2) & 1) as f64);
    let paths: [[usize; 4]; 6] = [
        [0, 1, 3, 7],
        [0, 1, 5, 7],
        [0, 2, 3, 7],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 6, 7],
    ];
    let mut elems3 = Vec::new();
    let mut elems_vec = Vec::new();
    let mut total_volume = 0.0;
    for p in paths {
        let tet = fea::TetrahedronElement3D::new(p, p.map(corner), mat);
        total_volume += tet.volume();
        elems3.push((tet.local_stiffness_matrix(), tet.dof_map()));
        elems_vec.push((tet.local_stiffness_matrix(), tet.dof_map().to_vec()));
    }
    assert!((total_volume - 1.0).abs() < 1e-14);
    let k3 = fea::assemble_3d_stiffness_matrix(24, &elems3);
    let kg = fea::assemble_stiffness_matrix(24, &elems_vec);
    for i in 0..24 {
        for j in 0..24 {
            assert!((k3.get(i, j) - kg.get(i, j)).abs() < 1e-12);
        }
    }
    // uniaxial stress state: u = e (x, -nu y, -nu z) has energy E e^2 / 2 per volume
    let (e, nu) = (1e-3, 0.3);
    let mut u = [0.0; 24];
    for n in 0..8 {
        let (x, y, z) = corner(n);
        u[3 * n] = e * x;
        u[3 * n + 1] = -nu * e * y;
        u[3 * n + 2] = -nu * e * z;
    }
    let energy: f64 = 0.5
        * (0..24)
            .map(|i| (0..24).map(|j| u[i] * k3.get(i, j) * u[j]).sum::<f64>())
            .sum::<f64>();
    assert!(
        (energy - 0.5 * 2.0 * e * e).abs() < 1e-12,
        "energy {energy}"
    );
    // equilibrium: the force vector of this state has no component on y/z faces' interior
    // 1D and 2D assembly through the generic routine equal the legacy ones
    let bar = fea::LinearElement1D {
        length: 0.5,
        youngs_modulus: 3.0,
        area: 2.0,
    };
    let kb = bar.local_stiffness_matrix();
    let legacy1 =
        fea::assemble_global_stiffness_matrix(3, &[(kb.clone(), 0, 1), (kb.clone(), 1, 2)]);
    let alias1 = fea::assemble_1d_stiffness_matrix(3, &[(kb.clone(), 0, 1), (kb.clone(), 1, 2)]);
    let generic1 = fea::assemble_stiffness_matrix(3, &[(kb.clone(), vec![0, 1]), (kb, vec![1, 2])]);
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(legacy1.get(i, j), alias1.get(i, j));
            assert!((legacy1.get(i, j) - generic1.get(i, j)).abs() < 1e-15);
        }
    }
}
