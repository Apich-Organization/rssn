//! Tests for geometric and algebraic multigrid.

use rssn::kernels::krylov::Csr;
use rssn::kernels::krylov::cg;
use rssn::kernels::krylov::identity_precond;
use rssn::kernels::multigrid::Amg;
use rssn::kernels::multigrid::AmgOptions;
use rssn::kernels::multigrid::Cycle;
use rssn::kernels::multigrid::GeometricMg;
use rssn::kernels::multigrid::MgOptions;
use rssn::kernels::multigrid::Smoother;

fn exact_1d(n: usize) -> (Vec<f64>, Vec<f64>) {
    // -u'' = pi^2 sin(pi x): u = sin(pi x); discrete solution error O(h^2)
    let h = 1.0 / (n as f64 + 1.0);
    let f: Vec<f64> = (1..=n)
        .map(|i| std::f64::consts::PI.powi(2) * (std::f64::consts::PI * i as f64 * h).sin())
        .collect();
    let u: Vec<f64> = (1..=n)
        .map(|i| (std::f64::consts::PI * i as f64 * h).sin())
        .collect();
    (f, u)
}

#[test]
fn geometric_mg_converges_grid_independently() {
    for dim in 1..=3usize {
        let sizes: Vec<usize> = match dim {
            | 1 => vec![63, 255],
            | 2 => vec![15, 63],
            | _ => vec![7, 31],
        };
        let mut factors = Vec::new();
        for &n in &sizes {
            let mg = GeometricMg::new(dim, n, 0.0).unwrap();
            let f: Vec<f64> = (0..mg.size())
                .map(|i| 1.0 + (i as f64 * 0.7).sin())
                .collect();
            let r = mg.solve(
                &f,
                None,
                &MgOptions {
                    tol: 1e-10,
                    ..MgOptions::default()
                },
            );
            assert!(r.converged, "dim {dim} n {n}: {:?}", r.history);
            factors.push(r.mean_factor());
            assert!(
                r.mean_factor() < 0.25,
                "dim {dim} n {n}: factor {}",
                r.mean_factor()
            );
        }
        assert!(
            factors[1] < 1.5 * factors[0] + 0.05,
            "dim {dim}: {factors:?}"
        );
    }
}

#[test]
fn cycles_smoothers_and_fmg() {
    let n = 63;
    let mg = GeometricMg::new(2, n, 0.0).unwrap();
    let f: Vec<f64> = (0..mg.size()).map(|i| (i as f64 * 0.37).cos()).collect();
    let v = mg.solve(&f, None, &MgOptions::default());
    let w = mg.solve(
        &f,
        None,
        &MgOptions {
            cycle: Cycle::W,
            ..MgOptions::default()
        },
    );
    let j = mg.solve(
        &f,
        None,
        &MgOptions {
            smoother: Smoother::Jacobi,
            pre: 3,
            post: 3,
            ..MgOptions::default()
        },
    );
    assert!(v.converged && w.converged && j.converged);
    assert!(w.mean_factor() <= v.mean_factor() + 1e-12);
    // all agree
    for i in 0..f.len() {
        assert!((v.x[i] - w.x[i]).abs() < 1e-7 * (1.0 + v.x[i].abs()));
    }
    // FMG: one V-cycle per level reaches discretisation accuracy with O(N) work
    let mg1 = GeometricMg::new(1, 127, 0.0).unwrap();
    let (rhs, exact) = exact_1d(127);
    let opts = MgOptions {
        tol: 1e-14,
        max_cycles: 0,
        ..MgOptions::default()
    };
    let r = mg1.fmg(&rhs, 1, &opts);
    let err =
        r.x.iter()
            .zip(&exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
    let h2 = (1.0 / 128.0f64).powi(2);
    assert!(
        err < 0.2 * h2 * std::f64::consts::PI.powi(2),
        "fmg error {err}"
    );
    // the algebraic error after FMG is already below truncation error: compare with the converged solve
    let full = mg1.solve(&rhs, None, &MgOptions::default());
    let alg =
        r.x.iter()
            .zip(&full.x)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
    assert!(alg < err, "algebraic {alg} vs total {err}");
    // second-order convergence of the discretisation
    let errs: Vec<f64> = [31usize, 63, 127]
        .iter()
        .map(|&n| {
            let m = GeometricMg::new(1, n, 0.0).unwrap();
            let (rhs, ex) = exact_1d(n);
            let s = m.solve(&rhs, None, &MgOptions::default());
            s.x.iter()
                .zip(&ex)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        })
        .collect();
    assert!(
        (errs[0] / errs[1]) > 3.8 && (errs[1] / errs[2]) > 3.8,
        "{errs:?}"
    );
    // Helmholtz-type shift and invalid sizes
    let sh = GeometricMg::new(2, 31, 50.0).unwrap();
    let f = vec![1.0; sh.size()];
    let s = sh.solve(&f, None, &MgOptions::default());
    assert!(s.converged);
    assert!(GeometricMg::new(2, 30, 0.0).is_none() && GeometricMg::new(4, 7, 0.0).is_none());
}

fn laplace2d(
    nx: usize,
    aniso: f64,
) -> Csr {
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let v = j * nx + i;
            t.push((v, v, 2.0 + 2.0 * aniso));
            if i + 1 < nx {
                t.push((v, v + 1, -1.0));
                t.push((v + 1, v, -1.0));
            }
            if j + 1 < nx {
                t.push((v, v + nx, -aniso));
                t.push((v + nx, v, -aniso));
            }
        }
    }
    Csr::from_triplets(nx * nx, nx * nx, &t)
}

#[test]
fn amg_smoothed_aggregation() {
    let mut iters = Vec::new();
    for &nx in &[32usize, 64] {
        let a = laplace2d(nx, 1.0);
        let b: Vec<f64> = (0..a.rows).map(|i| 1.0 + (i as f64 * 0.11).sin()).collect();
        let amg = Amg::build(&a, &AmgOptions::default()).unwrap();
        assert!(
            amg.num_levels() >= 3 && amg.operator_complexity() < 2.0,
            "{} {}",
            amg.num_levels(),
            amg.operator_complexity()
        );
        let r = amg.solve(&b, None, 1e-8, 60);
        assert!(r.converged, "nx {nx}: {:?}", r.history.last());
        iters.push(r.history.len());
        // CG preconditioned by one AMG cycle
        let res = cg(
            |x, y| a.matvec(x, y),
            |r, z| amg.precondition(r, z),
            &b,
            None,
            1e-10,
            100,
        );
        let plain = cg(
            |x, y| a.matvec(x, y),
            identity_precond,
            &b,
            None,
            1e-10,
            2000,
        );
        assert!(
            res.converged && res.iterations * 3 < plain.iterations,
            "{} vs {}",
            res.iterations,
            plain.iterations
        );
    }
    assert!(iters[1] <= iters[0] + 6, "{iters:?}");
}

#[test]
fn amg_anisotropy_and_options() {
    let a = laplace2d(48, 0.01);
    let b = vec![1.0; a.rows];
    let amg = Amg::build(&a, &AmgOptions::default()).unwrap();
    let r = amg.solve(&b, None, 1e-8, 200);
    assert!(r.converged, "{:?}", r.history.last());
    // W cycle, plain aggregation (no prolongator smoothing) still converge
    let w = Amg::build(
        &a,
        &AmgOptions {
            cycle: Cycle::W,
            ..AmgOptions::default()
        },
    )
    .unwrap();
    assert!(w.solve(&b, None, 1e-8, 300).converged);
    // vector-valued near-null space (two copies of a decoupled problem)
    let n = a.rows;
    let mut t = Vec::new();
    for i in 0..n {
        for k in a.indptr[i]..a.indptr[i + 1] {
            t.push((i, a.indices[k], a.values[k]));
            t.push((n + i, n + a.indices[k], a.values[k]));
        }
    }
    let two = Csr::from_triplets(2 * n, 2 * n, &t);
    let mut z1 = vec![0.0; 2 * n];
    let mut z2 = vec![0.0; 2 * n];
    for i in 0..n {
        z1[i] = 1.0;
        z2[n + i] = 1.0;
    }
    let opts = AmgOptions {
        near_nullspace: vec![z1, z2],
        ..AmgOptions::default()
    };
    let amg2 = Amg::build(&two, &opts).unwrap();
    let rhs = vec![1.0; 2 * n];
    assert!(amg2.solve(&rhs, None, 1e-8, 300).converged);
    assert!(
        Amg::build(
            &Csr::from_triplets(2, 3, &[(0, 0, 1.0)]),
            &AmgOptions::default()
        )
        .is_err()
    );
}
