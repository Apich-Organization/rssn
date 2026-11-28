//! Crank–Nicolson ADI tests for homogeneous and prescribed fixed boundaries.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use rssn::sim::physics_cnm::HeatEquationSolverConfig;
use rssn::sim::physics_cnm::solve_heat_equation_2d_cn_adi;
use rssn::sim::physics_cnm::solve_heat_equation_2d_cn_adi_fixed_boundary;

fn config(
    nx: usize,
    ny: usize,
    steps: usize,
) -> HeatEquationSolverConfig {
    HeatEquationSolverConfig {
        nx,
        ny,
        dx: 1.0,
        dy: 2.0,
        dt: 0.5,
        d_coeff: 1.0,
        steps,
    }
}

#[test]
fn legacy_projects_the_initial_boundary_before_diffusion() {
    let boundary = [3.0, -2.0, 1.0, 4.0, 0.0, 2.0, -4.0, -9.0, 7.0];
    assert_eq!(
        solve_heat_equation_2d_cn_adi(&boundary, &config(3, 3, 0)),
        boundary
    );
    assert_eq!(
        solve_heat_equation_2d_cn_adi(&boundary, &config(3, 3, 1)),
        [0.0; 9]
    );
    assert_eq!(
        solve_heat_equation_2d_cn_adi(&[2.0; 4], &config(1, 4, 1)),
        [0.0; 4]
    );
    assert!(solve_heat_equation_2d_cn_adi(&[], &config(0, 4, 1)).is_empty());
}

#[test]
fn prescribed_faces_supply_both_half_step_equations() {
    let initial = [0.0, 2.0, 0.0, 3.0, -0.5, 5.0, 0.0, 4.0, 0.0];
    let cfg = config(3, 3, 1);
    let next = solve_heat_equation_2d_cn_adi_fixed_boundary(&initial, &cfg).unwrap();
    let (alpha_x, alpha_y) = (0.25, 0.0625);
    let x_edges = alpha_x * (3.0 + 5.0);
    let y_edges = alpha_y * (2.0 + 4.0);
    let half = (x_edges + y_edges + (1.0 - 2.0 * alpha_y) * initial[4]) / (1.0 + 2.0 * alpha_x);
    let exact = (x_edges + y_edges + (1.0 - 2.0 * alpha_x) * half) / (1.0 + 2.0 * alpha_y);
    assert!((next[4] - exact).abs() < 2.0e-15);
    for p in (0..initial.len()).filter(|&p| p != 4) {
        assert_eq!(next[p].to_bits(), initial[p].to_bits());
    }
}

#[test]
fn fixed_bilinear_temperature_is_steady_on_anisotropic_grids() {
    let cfg = config(17, 9, 42);
    let initial: Vec<f64> = (0..cfg.nx * cfg.ny)
        .map(|p| {
            let (x, y) = ((p % cfg.nx) as f64, (p / cfg.nx) as f64);
            1.25 - 0.5 * x + 0.25 * y + 0.0625 * x * y
        })
        .collect();
    let next = solve_heat_equation_2d_cn_adi_fixed_boundary(&initial, &cfg).unwrap();
    for (p, (&value, &expected)) in next.iter().zip(&initial).enumerate() {
        let (i, j) = (p % cfg.nx, p / cfg.nx);
        if i == 0 || j == 0 || i + 1 == cfg.nx || j + 1 == cfg.ny {
            assert_eq!(value.to_bits(), expected.to_bits());
        } else {
            assert!((value - expected).abs() < 2.0e-13, "node ({i}, {j})");
        }
    }
}

#[test]
fn homogeneous_heat_is_a_contraction_with_stiff_time_steps() {
    let mut cfg = config(13, 25, 1);
    cfg.dt = 50.0;
    let mut u: Vec<f64> = (0..cfg.nx * cfg.ny)
        .map(|p| {
            let (i, j) = (p % cfg.nx, p / cfg.nx);
            if i == 0 || j == 0 || i + 1 == cfg.nx || j + 1 == cfg.ny {
                0.0
            } else {
                ((p * 29) % 53) as f64 / 32.0 - 1.0
            }
        })
        .collect();
    let norm2 = |values: &[f64]| values.iter().map(|x| x * x).sum::<f64>();
    for _ in 0..12 {
        let next = solve_heat_equation_2d_cn_adi_fixed_boundary(&u, &cfg).unwrap();
        assert!(norm2(&next) < norm2(&u));
        u = next;
    }
}

#[test]
fn both_solvers_recover_from_extreme_scale_products() {
    let mut initial = [0.0; 9];
    initial[4] = 1.0;
    for scale in [2.0_f64.powi(-600), 2.0_f64.powi(600)] {
        let mut cfg = config(3, 3, 1);
        cfg.dx = scale;
        cfg.dy = scale;
        cfg.d_coeff = scale;
        cfg.dt = 0.5 * scale;
        let next = solve_heat_equation_2d_cn_adi_fixed_boundary(&initial, &cfg).unwrap();
        assert!((next[4] - 1.0 / 9.0).abs() < 1.0e-15);
        assert_eq!(solve_heat_equation_2d_cn_adi(&initial, &cfg), next);
    }
}

#[test]
fn checked_solver_validates_and_handles_stationary_cases_without_time_stepping() {
    let mut cfg = config(2, 3, usize::MAX);
    let faces = [1.0, 2.0, -3.0, -4.0, 5.0, 6.0];
    assert_eq!(
        solve_heat_equation_2d_cn_adi_fixed_boundary(&faces, &cfg).unwrap(),
        faces
    );
    assert!(solve_heat_equation_2d_cn_adi_fixed_boundary(&faces[..4], &cfg).is_err());
    cfg.nx = usize::MAX;
    assert!(solve_heat_equation_2d_cn_adi_fixed_boundary(&faces, &cfg).is_err());
    cfg.nx = 3;
    let mut field = [1.0; 9];
    cfg.dt = 0.0;
    assert_eq!(
        solve_heat_equation_2d_cn_adi_fixed_boundary(&field, &cfg).unwrap(),
        field
    );
    cfg.dx = f64::INFINITY;
    assert!(solve_heat_equation_2d_cn_adi_fixed_boundary(&field, &cfg).is_err());
    cfg.dx = 1.0;
    field[4] = f64::NAN;
    assert!(solve_heat_equation_2d_cn_adi_fixed_boundary(&field, &cfg).is_err());
    let _cancelled = rssn::cancellation::install(Arc::new(AtomicBool::new(true)));
    assert!(solve_heat_equation_2d_cn_adi_fixed_boundary(&[1.0; 9], &cfg).is_err());
}

#[test]
fn parallel_sweeps_and_transposes_keep_the_same_bits() {
    let cfg = config(131, 151, 4);
    let initial: Vec<f64> = (0..cfg.nx * cfg.ny)
        .map(|p| ((p * 37) % 173) as f64 / 256.0)
        .collect();
    let run = |workers| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| solve_heat_equation_2d_cn_adi_fixed_boundary(&initial, &cfg).unwrap())
    };
    for (serial, parallel) in run(1).into_iter().zip(run(2)) {
        assert_eq!(serial.to_bits(), parallel.to_bits());
    }
}
