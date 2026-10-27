//! Tests for algebraic multigrid on non-symmetric systems.

use rssn::kernels::amg_general::{AmgGeneral, AmgGeneralOptions, AmgKind};
use rssn::kernels::krylov::{Csr, gmres, identity_precond};
use rssn::kernels::multigrid::Cycle;

/// Upwind finite-difference matrix of `-eps Laplace u + b . grad u` on the unit square.
fn convection_diffusion(nx: usize, eps: f64, bx: f64, by: f64) -> Csr {
    let h = 1.0 / (nx as f64 + 1.0);
    let d = eps / (h * h);
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0 * d + (bx + by) / h));
            if i > 0 {
                t.push((p, p - 1, -d - bx / h));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -d));
            }
            if j > 0 {
                t.push((p, p - nx, -d - by / h));
            }
            if j + 1 < nx {
                t.push((p, p + nx, -d));
            }
        }
    }
    Csr::from_triplets(nx * nx, nx * nx, &t)
}

fn rhs(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1.0 + ((i * 7 % 13) as f64) / 13.0).collect()
}

fn gmres_iterations(a: &Csr, amg: Option<&AmgGeneral>, b: &[f64]) -> (usize, bool, f64) {
    let r = match amg {
        Some(m) => gmres(|x, y| a.matvec(x, y), |r, z| m.precondition(r, z), b, None, 60, 1e-9, 300),
        None => gmres(|x, y| a.matvec(x, y), identity_precond, b, None, 60, 1e-9, 300),
    };
    let mut ax = vec![0.0; b.len()];
    a.matvec(&r.x, &mut ax);
    let res = ax.iter().zip(b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt()
        / b.iter().map(|v| v * v).sum::<f64>().sqrt();
    (r.iterations, r.converged, res)
}

#[test]
fn ruge_stuben_preconditions_convection_diffusion_mesh_independently() {
    for &(eps, bx, by) in &[(1.0, 1.0, 0.5), (0.02, 1.0, 0.5), (0.001, 0.6, 0.0)] {
        let mut its = Vec::new();
        for nx in [32usize, 64, 128] {
            let a = convection_diffusion(nx, eps, bx, by);
            let amg = AmgGeneral::build(&a, &AmgGeneralOptions::default()).unwrap();
            assert!(amg.num_levels() >= 3);
            assert!(amg.operator_complexity() < 4.0, "complexity {}", amg.operator_complexity());
            let b = rhs(nx * nx);
            let (it, conv, res) = gmres_iterations(&a, Some(&amg), &b);
            assert!(conv && res < 1e-8, "eps={eps} nx={nx}: converged={conv}, residual {res}");
            its.push(it);
        }
        eprintln!("RS eps={eps}: iterations {its:?}");
        assert!(its[2] <= its[0] * 2 + 4, "eps={eps}: iterations grow with the mesh: {its:?}");
        assert!(its[2] <= 40, "{its:?}");
    }
    // unpreconditioned GMRES needs far more work on the finest mesh
    let a = convection_diffusion(128, 1.0, 1.0, 0.5);
    let b = rhs(128 * 128);
    let (plain, _, _) = gmres_iterations(&a, None, &b);
    let amg = AmgGeneral::build(&a, &AmgGeneralOptions::default()).unwrap();
    let (pre, conv, _) = gmres_iterations(&a, Some(&amg), &b);
    assert!(conv && pre * 5 < plain, "{pre} vs {plain}");
}

#[test]
fn petrov_galerkin_aggregation_preconditions_convection_diffusion() {
    for &(eps, bx, by) in &[(1.0, 1.0, 0.5), (0.05, 1.0, 0.5)] {
        let mut its = Vec::new();
        for nx in [32usize, 64, 128] {
            let a = convection_diffusion(nx, eps, bx, by);
            let opts = AmgGeneralOptions {
                kind: AmgKind::SmoothedAggregation,
                theta: 0.08,
                ..AmgGeneralOptions::default()
            };
            let amg = AmgGeneral::build(&a, &opts).unwrap();
            assert!(amg.num_levels() >= 3);
            let b = rhs(nx * nx);
            let (it, conv, res) = gmres_iterations(&a, Some(&amg), &b);
            assert!(conv && res < 1e-8, "eps={eps} nx={nx}: converged={conv}, residual {res}");
            its.push(it);
        }
        eprintln!("SA eps={eps}: iterations {its:?}");
        assert!(its[2] <= its[0] * 2 + 6, "eps={eps}: {its:?}");
    }
}

#[test]
fn stand_alone_cycles_converge_for_diffusion_dominated_problems() {
    let a = convection_diffusion(64, 1.0, 0.5, 0.25);
    let b = rhs(64 * 64);
    for kind in [AmgKind::RugeStuben, AmgKind::SmoothedAggregation] {
        let amg = AmgGeneral::build(&a, &AmgGeneralOptions { kind, ..AmgGeneralOptions::default() }).unwrap();
        let r = amg.solve(&b, None, 1e-8, 80);
        let rate = (r.history[r.history.len() - 1] / r.history[0]).powf(1.0 / (r.history.len() - 1) as f64);
        eprintln!("{kind:?}: cycles {}, rate {rate}", r.history.len() - 1);
        assert!(r.converged, "{kind:?}: {:?}", r.history.last());
        assert!(rate < 0.85, "{kind:?}: convergence factor {rate}");
    }
    let w = AmgGeneral::build(&a, &AmgGeneralOptions { cycle: Cycle::W, ..AmgGeneralOptions::default() }).unwrap();
    assert!(w.solve(&b, None, 1e-8, 60).converged);
}

#[test]
fn invalid_matrices_are_rejected() {
    let rect = Csr::from_triplets(2, 3, &[(0, 0, 1.0)]);
    assert!(AmgGeneral::build(&rect, &AmgGeneralOptions::default()).is_err());
    let a = convection_diffusion(20, 1.0, 1.0, 1.0);
    let bad = AmgGeneralOptions {
        kind: AmgKind::SmoothedAggregation,
        near_nullspace: vec![vec![1.0; 3]],
        ..AmgGeneralOptions::default()
    };
    assert!(AmgGeneral::build(&a, &bad).is_err());
    let neg = Csr::from_triplets(100, 100, &(0..100).map(|i| (i, i, -1.0)).collect::<Vec<_>>());
    assert!(AmgGeneral::build(&neg, &AmgGeneralOptions::default()).is_err());
}
