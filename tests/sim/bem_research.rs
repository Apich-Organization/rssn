//! Helmholtz BEM (Burton-Miller) against the Mie series, and ACA-compressed P1 Laplace BEM.

use num_complex::Complex64;
use rssn::sim::bem_aca::{AcaOptions, solve_laplace_bem_3d_p1, solve_laplace_bem_3d_p1_aca};
use rssn::sim::bem_helmholtz::{PlaneWave, Scatterer, solve_helmholtz_exterior};
use rssn::sim::physics_bem::{BoundaryCondition, SurfaceMesh3D};

/// Spherical Bessel functions `j_0..=j_n` (Miller backward recurrence).
fn sph_j(nmax: usize, x: f64) -> Vec<f64> {
    let start = nmax + 30;
    let mut j = vec![0.0; start + 2];
    j[start] = 1e-30;
    for n in (0..start).rev() {
        j[n] = (2.0 * n as f64 + 3.0) / x * j[n + 1] - j[n + 2];
    }
    let scale = (x.sin() / x) / j[0];
    j.iter().take(nmax + 1).map(|v| v * scale).collect()
}

/// Spherical Neumann functions `y_0..=y_n`.
fn sph_y(nmax: usize, x: f64) -> Vec<f64> {
    let mut y = vec![-x.cos() / x, -x.cos() / (x * x) - x.sin() / x];
    for n in 1..nmax {
        y.push((2.0 * n as f64 + 1.0) / x * y[n] - y[n - 1]);
    }
    y.truncate(nmax + 1);
    y
}

fn legendre(nmax: usize, c: f64) -> Vec<f64> {
    let mut p = vec![1.0, c];
    for n in 1..nmax {
        p.push(((2.0 * n as f64 + 1.0) * c * p[n] - n as f64 * p[n - 1]) / (n as f64 + 1.0));
    }
    p.truncate(nmax + 1);
    p
}

/// Scattered field of a unit sphere for the incident wave `e^{ikz}` at radius `r` and polar angle `theta`.
fn mie_scattered(k: f64, kind: Scatterer, r: f64, theta: f64) -> Complex64 {
    let nmax = 24;
    let (ja, ya) = (sph_j(nmax + 1, k), sph_y(nmax + 1, k));
    let (jr, yr) = (sph_j(nmax, k * r), sph_y(nmax, k * r));
    let pl = legendre(nmax, theta.cos());
    let mut s = Complex64::new(0.0, 0.0);
    for n in 0..=nmax {
        let h_a = Complex64::new(ja[n], ya[n]);
        let h_r = Complex64::new(jr[n], yr[n]);
        let ratio = match kind {
            Scatterer::SoundSoft => Complex64::new(ja[n], 0.0) / h_a,
            Scatterer::SoundHard => {
                let djn = if n == 0 { -ja[1] } else { ja[n - 1] - (n as f64 + 1.0) / k * ja[n] };
                let dyn_ = if n == 0 { -ya[1] } else { ya[n - 1] - (n as f64 + 1.0) / k * ya[n] };
                Complex64::new(djn, 0.0) / Complex64::new(djn, dyn_)
            },
        };
        let i_n = Complex64::new(0.0, 1.0).powi(n as i32);
        s -= i_n * (2.0 * n as f64 + 1.0) * ratio * h_r * pl[n];
    }
    s
}

fn scattering_error(sub: usize, k: f64, kind: Scatterer, coupling: Option<Complex64>) -> f64 {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], sub);
    let sol = solve_helmholtz_exterior(&mesh, k, kind, PlaneWave { direction: [0.0, 0.0, 1.0] }, coupling).unwrap();
    let mut worst: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for &theta in &[0.0, 0.7, 1.4, 2.1, std::f64::consts::PI] {
        let r = 2.5;
        let x = [r * theta.sin(), 0.0, r * theta.cos()];
        let exact = mie_scattered(k, kind, r, theta);
        let num = sol.scattered_field(&mesh, x);
        worst = worst.max((num - exact).norm());
        scale = scale.max(exact.norm());
    }
    worst / scale
}

#[test]
fn burton_miller_scattering_by_a_sphere_matches_the_mie_series() {
    for kind in [Scatterer::SoundHard, Scatterer::SoundSoft] {
        for k in [1.0, 2.0] {
            let e2 = scattering_error(2, k, kind, None);
            let e3 = scattering_error(3, k, kind, None);
            println!("{kind:?} ka={k}: error sub2 {e2:.4} sub3 {e3:.4}");
            assert!(e3 < 0.06, "{kind:?} ka = {k}: relative error {e3}");
            assert!(e3 < e2, "{kind:?} ka = {k}: no convergence {e2} -> {e3}");
        }
    }
}

#[test]
fn burton_miller_stays_accurate_at_an_irregular_frequency() {
    // ka = pi is an eigenfrequency of the interior Dirichlet problem of the unit sphere:
    // the conventional equation for a sound-soft body is (almost) singular there
    let k = std::f64::consts::PI;
    let bm = scattering_error(3, k, Scatterer::SoundSoft, None);
    let conv = scattering_error(3, k, Scatterer::SoundSoft, Some(Complex64::new(0.0, 0.0)));
    println!("irregular frequency: BM {bm:.4}, conventional {conv:.4}");
    assert!(bm < 0.08, "Burton-Miller error {bm}");
    assert!(conv > 3.0 * bm, "conventional {conv} should be much worse than {bm}");
}

fn aca_opts() -> AcaOptions {
    AcaOptions { leaf_size: 24, eta: 0.8, aca_tol: 1e-7, gmres_tol: 1e-10, ..AcaOptions::default() }
}

#[test]
fn p1_aca_matches_the_dense_p1_solver() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);
    let opts = AcaOptions { leaf_size: 8, aca_tol: 1e-9, ..aca_opts() };
    // interior Dirichlet u = x y + z, exterior capacitance problem and mixed data
    let cases: Vec<(Vec<BoundaryCondition<f64>>, bool)> = vec![
        (mesh.vertices.iter().map(|v| BoundaryCondition::Potential(v[0] * v[1] + v[2])).collect(), false),
        (mesh.vertices.iter().map(|_| BoundaryCondition::Potential(1.0)).collect(), true),
        (
            mesh.vertices
                .iter()
                .map(|v| if v[2] < 0.0 { BoundaryCondition::Potential(v[0]) } else { BoundaryCondition::Flux(v[0]) })
                .collect(),
            false,
        ),
    ];
    for (bcs, exterior) in cases {
        let dense = solve_laplace_bem_3d_p1(&mesh, &bcs, exterior).unwrap();
        let (aca, stats) = solve_laplace_bem_3d_p1_aca(&mesh, &bcs, exterior, &opts).unwrap();
        let scale = dense.u.iter().chain(&dense.q).fold(0.0_f64, |m, v| m.max(v.abs()));
        let diff = dense.u.iter().zip(&aca.u).chain(dense.q.iter().zip(&aca.q)).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(diff < 1e-5 * scale.max(1.0), "exterior {exterior}: dense/ACA difference {diff} ({stats:?})");
    }
}

#[test]
fn p1_aca_compresses_and_converges_on_a_larger_sphere() {
    // u = x is harmonic: the interior Dirichlet problem has the flux q = x
    let flux_error = |sub: usize| {
        let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], sub);
        let bcs: Vec<BoundaryCondition<f64>> = mesh.vertices.iter().map(|v| BoundaryCondition::Potential(v[0])).collect();
        let (sol, stats) = solve_laplace_bem_3d_p1_aca(&mesh, &bcs, false, &aca_opts()).unwrap();
        let worst = mesh.vertices.iter().zip(&sol.q).map(|(v, q)| (q - v[0]).abs()).fold(0.0, f64::max);
        (worst, stats.compression(), mesh.vertices.len())
    };
    let (e3, _, _) = flux_error(3);
    let (e4, comp, nv) = flux_error(4);
    println!("P1 ACA: {nv} vertices, flux error {e3:.5} -> {e4:.5}, stored fraction {comp:.3}");
    assert!(e4 < 0.01 && e4 < 0.6 * e3, "convergence {e3} -> {e4}");
    assert!(comp < 0.9, "stored fraction {comp}");
}
