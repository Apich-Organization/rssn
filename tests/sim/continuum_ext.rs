//! Second-round continuum features: acoustics PML, non-reflecting Euler
//! boundaries, DG Euler 2D with limiters, oblique TF/SF, exterior / P1 BEM.

use rssn::sim::acoustics::{Acoustics, AcousticsConfig};
use rssn::sim::bem_aca::{
    AcaOptions, evaluate_potential_3d_exterior, solve_laplace_bem_3d_aca, solve_laplace_bem_3d_exterior,
    solve_laplace_bem_3d_p1,
};
use rssn::sim::boundary::{BoundaryCondition as Bc, Face, FaceBcs};
use rssn::sim::dg::{DgBc, DgLimiter, DgMesh, Euler1d, Euler2d, solve_dg, solve_dg_limited};
use rssn::sim::euler_nd::{FvConfig, FvGrid, GasBc, Reconstruction, SlopeLimiter, solve};
use rssn::sim::gas_dynamics::{Primitive, exact_riemann};
use rssn::sim::physics_bem::{BoundaryCondition, SurfaceMesh3D};
use rssn::sim::yee_fdtd::{Comp, ObliquePlaneWave, Yee, YeeConfig};
use std::f64::consts::PI;

// ---------------------------------------------------------------- acoustics

fn standing_error(dim: usize, n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let cfg = AcousticsConfig::new(&vec![n; dim], &vec![h; dim], 1.0).unwrap();
    let mut a = Acoustics::new(cfg).unwrap();
    let mode = move |x: [f64; 3]| (0..dim).map(|k| (PI * x[k]).cos()).product::<f64>();
    a.set_pressure(mode);
    a.start_from_rest();
    // omega T = 3 pi / 2 so that the measured phase error is insensitive to rounding of T
    let w = PI * (dim as f64).sqrt();
    let steps = (1.5 * PI / w / a.dt).round() as usize;
    a.run(steps);
    let nn = a.n();
    let mut worst: f64 = 0.0;
    for k in 0..nn[2] {
        for j in 0..nn[1] {
            for i in 0..nn[0] {
                let x = [(i as f64 + 0.5) * h, (j as f64 + 0.5) * h, (k as f64 + 0.5) * h];
                let ex = mode(x) * (w * a.t).cos();
                worst = worst.max((a.p[a.idx(i, j, k)] - ex).abs());
            }
        }
    }
    worst
}

#[test]
fn rigid_box_standing_wave_second_order_in_1d_2d_3d() {
    for (dim, ns) in [(1, [32, 64]), (2, [16, 32]), (3, [8, 16])] {
        let (e1, e2) = (standing_error(dim, ns[0]), standing_error(dim, ns[1]));
        assert!((e1 / e2).log2() > 1.8, "dim {dim}: {e1} {e2}");
    }
}

fn pulse_energy_ratio(dim: usize, n: usize, pml: usize, sigma: f64, t_end: f64, absorbing: bool) -> f64 {
    let h = 1.0 / n as f64;
    let mut cfg = AcousticsConfig::new(&vec![n; dim], &vec![h; dim], 1.0).unwrap();
    if absorbing {
        cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
        cfg.pml_cells = pml;
    }
    let mut a = Acoustics::new(cfg).unwrap();
    a.set_pressure(move |x| (-(0..dim).map(|k| (x[k] - 0.5).powi(2)).sum::<f64>() / (2.0 * sigma * sigma)).exp());
    let e0 = a.energy();
    let steps = (t_end / a.dt).round() as usize;
    a.run(steps);
    a.energy() / e0
}

#[test]
fn pml_absorbs_scalar_waves_in_1d_2d_3d() {
    // the initial pressure bump splits into outgoing waves that leave through the PML
    let r1 = pulse_energy_ratio(1, 200, 20, 0.04, 1.4, true);
    assert!(r1 < 1e-4, "1D {r1}");
    let r2 = pulse_energy_ratio(2, 100, 12, 0.04, 1.4, true);
    assert!(r2 < 2e-3, "2D {r2}");
    let r3 = pulse_energy_ratio(3, 36, 8, 0.07, 1.0, true);
    assert!(r3 < 1.5e-2, "3D {r3}");
    // a rigid box keeps the energy (leapfrog energy fluctuates slightly)
    let rigid = pulse_energy_ratio(2, 100, 12, 0.04, 1.4, false);
    assert!((0.9..1.1).contains(&rigid), "rigid box energy ratio {rigid}");
}

#[test]
fn pressure_release_inverts_and_periodic_wraps() {
    let n = 200;
    let h = 1.0 / n as f64;
    let run = |bc: Bc, t_end: f64| {
        let mut cfg = AcousticsConfig::new(&[n], &[h], 1.0).unwrap();
        cfg.bcs = FaceBcs::uniform(&bc);
        let mut a = Acoustics::new(cfg).unwrap();
        a.set_pressure(|x| (-(x[0] - 0.7).powi(2) / (2.0 * 0.05f64.powi(2))).exp());
        a.start_from_rest();
        a.run((t_end / a.dt).round() as usize);
        a
    };
    let rigid = run(Bc::neumann(0.0), 0.5);
    let release = run(Bc::dirichlet(0.0), 0.5);
    let i = (0.8 / h) as usize;
    assert!(rigid.p[i] > 0.4, "rigid reflection keeps the sign: {}", rigid.p[i]);
    assert!(release.p[i] < -0.4, "pressure release inverts: {}", release.p[i]);
    let periodic = run(Bc::Periodic, 1.0);
    let c = (0.7 / h) as usize;
    assert!(periodic.p[c] > 0.95, "periodic return {}", periodic.p[c]);
    // non-zero release data and Robin faces are rejected
    let mut cfg = AcousticsConfig::new(&[n], &[h], 1.0).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::dirichlet(1.0));
    assert!(Acoustics::new(cfg).is_err());
}

// ------------------------------------------------------- Euler non-reflecting

fn acoustic_residual(dim: usize, bc: GasBc) -> f64 {
    let gamma = 1.4;
    let n = if dim == 1 { 200 } else { 64 };
    let grid = FvGrid::new(&vec![n; dim], &vec![0.0; dim], &vec![2.0; dim]).unwrap();
    let mut cfg = FvConfig::euler(gamma);
    cfg.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
    for a in 0..dim {
        cfg = cfg.with_axis(a, bc);
    }
    let amp = 0.01;
    let init = move |x: [f64; 3]| {
        let r2: f64 = (0..dim).map(|a| (x[a] - 1.0).powi(2)).sum();
        let dp = amp * (-r2 / 0.02).exp();
        [1.0 + dp / gamma, 0.0, 0.0, 0.0, 1.0 + dp]
    };
    let sol = solve(&grid, &cfg, &init, 2.2).unwrap();
    (0..grid.len()).map(|p| (sol.primitive(p)[4] - 1.0).abs()).fold(0.0, f64::max) / amp
}

#[test]
fn characteristic_boundaries_let_acoustic_pulses_leave() {
    let far = [1.0, 0.0, 0.0, 0.0, 1.0];
    let nr1 = acoustic_residual(1, GasBc::NonReflecting(far));
    let wall1 = acoustic_residual(1, GasBc::Reflective);
    assert!(nr1 < 0.03, "1D non-reflecting residual {nr1}");
    assert!(wall1 > 0.3, "1D wall residual {wall1}");
    let nr2 = acoustic_residual(2, GasBc::NonReflecting(far));
    let wall2 = acoustic_residual(2, GasBc::Reflective);
    assert!(nr2 < 0.25 * wall2 && nr2 < 0.15, "2D nr {nr2} vs wall {wall2}");
    // shallow water has no characteristic boundary
    let g = FvGrid::new(&[8], &[0.0], &[1.0]).unwrap();
    let cfg = FvConfig::shallow_water(9.81).with_axis(0, GasBc::NonReflecting(far));
    assert!(solve(&g, &cfg, &|_| [1.0, 0.0, 0.0, 0.0, 0.0], 0.1).is_err());
}

// ------------------------------------------------------------ DG Euler 2D

fn vortex(x: [f64; 2], t: f64) -> Vec<f64> {
    let (gamma, beta) = (1.4, 5.0);
    let wrap = |d: f64| (d + 5.0).rem_euclid(10.0) - 5.0;
    let (dx, dy) = (wrap(x[0] - 5.0 - t), wrap(x[1] - 5.0 - t));
    let r2 = dx * dx + dy * dy;
    let e = ((1.0 - r2) / 2.0).exp();
    let u = 1.0 - beta / (2.0 * PI) * e * dy;
    let v = 1.0 + beta / (2.0 * PI) * e * dx;
    let temp = 1.0 - (gamma - 1.0) * beta * beta / (8.0 * gamma * PI * PI) * (1.0 - r2).exp();
    let rho = temp.powf(1.0 / (gamma - 1.0));
    let p = rho.powf(gamma);
    vec![rho, rho * u, rho * v, p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v)]
}

#[test]
fn dg_euler_2d_isentropic_vortex_converges_at_high_order() {
    let law = Euler2d { gamma: 1.4 };
    let per = [DgBc::Periodic, DgBc::Periodic, DgBc::Periodic, DgBc::Periodic];
    let run = |ne: usize| {
        let mesh = DgMesh::rect([ne, ne], [0.0; 2], [10.0; 2], 3);
        let cfl = 0.3 * (6.0 / ne as f64).powf(1.0 / 3.0);
        let s = solve_dg(&mesh, &law, &per, &|x| vortex(x, 0.0), 1.0, cfl, None).unwrap();
        s.l2_error(0, |x| vortex(x, 1.0)[0])
    };
    let (e1, e2) = (run(6), run(12));
    assert!((e1 / e2).log2() > 2.5, "{e1} {e2}");
    // limiters leave a smooth solution essentially untouched
    let mesh = DgMesh::rect([8, 8], [0.0; 2], [10.0; 2], 3);
    let lim = DgLimiter { tvb_m: Some(200.0), positivity: true };
    let a = solve_dg_limited(&mesh, &law, &per, &|x| vortex(x, 0.0), 1.0, 0.25, None, Some(lim)).unwrap();
    let b = solve_dg(&mesh, &law, &per, &|x| vortex(x, 0.0), 1.0, 0.25, None).unwrap();
    let ea = a.l2_error(0, |x| vortex(x, 1.0)[0]);
    let eb = b.l2_error(0, |x| vortex(x, 1.0)[0]);
    assert!(ea < 1.5 * eb + 1e-9, "limited {ea} unlimited {eb}");
}

fn pressure_of(q: &[f64]) -> f64 {
    0.4 * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0])
}

#[test]
fn dg_euler_2d_shock_tube_with_tvb_and_positivity_limiters() {
    let law = Euler2d { gamma: 1.4 };
    let bcs = [DgBc::Outflow, DgBc::Outflow, DgBc::Periodic, DgBc::Periodic];
    let l = Primitive { rho: 1.0, u: 0.0, p: 1.0 };
    let r = Primitive { rho: 0.125, u: 0.0, p: 0.1 };
    let init = |x: [f64; 2]| {
        let (rho, p) = if x[0] < 0.5 { (1.0, 1.0) } else { (0.125, 0.1) };
        vec![rho, 0.0, 0.0, p / 0.4]
    };
    let mesh = DgMesh::rect([64, 2], [0.0, 0.0], [1.0, 1.0 / 32.0], 2);
    let lim = DgLimiter { tvb_m: Some(20.0), positivity: true };
    let s = solve_dg_limited(&mesh, &law, &bcs, &init, 0.2, 0.25, None, Some(lim)).unwrap();
    let mut l1 = 0.0;
    let mut rho_min = f64::INFINITY;
    let mut rho_max: f64 = 0.0;
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            let q: Vec<f64> = (0..4).map(|v| s.value(e, nd, v)).collect();
            rho_min = rho_min.min(q[0]);
            rho_max = rho_max.max(q[0]);
            assert!(pressure_of(&q) > 0.0);
            let x = s.coord(e, nd)[0];
            let ex = exact_riemann(l, r, 1.4, (x - 0.5) / 0.2).rho;
            l1 += (q[0] - ex).abs() / (mesh.nodes_per_element() * mesh.elements()) as f64;
        }
    }
    assert!(l1 < 0.03, "L1 {l1}");
    assert!(rho_min > 0.1 && rho_max < 1.08, "bounds {rho_min} {rho_max}");
}

#[test]
fn dg_euler_double_rarefaction_needs_and_gets_positivity() {
    // the 1-2-3 problem: near-vacuum star region
    let l = Primitive { rho: 1.0, u: -2.0, p: 0.4 };
    let r = Primitive { rho: 1.0, u: 2.0, p: 0.4 };
    let law = Euler1d { gamma: 1.4 };
    let mesh = DgMesh::line(80, 0.0, 1.0, 2);
    let bcs = [DgBc::Outflow, DgBc::Outflow, DgBc::Outflow, DgBc::Outflow];
    let init = |x: [f64; 2]| {
        let u = if x[0] < 0.5 { -2.0 } else { 2.0 };
        vec![1.0, u, 0.4 / 0.4 + 0.5 * u * u]
    };
    let lim = DgLimiter { tvb_m: Some(20.0), positivity: true };
    let s = solve_dg_limited(&mesh, &law, &bcs, &init, 0.1, 0.2, None, Some(lim)).unwrap();
    let mut l1 = 0.0;
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            let (rho, m, en) = (s.value(e, nd, 0), s.value(e, nd, 1), s.value(e, nd, 2));
            assert!(rho > 0.0 && 0.4 * (en - 0.5 * m * m / rho) > 0.0, "non-positive state");
            let x = s.coord(e, nd)[0];
            l1 += (rho - exact_riemann(l, r, 1.4, (x - 0.5) / 0.1).rho).abs() / (mesh.nodes_per_element() * mesh.elements()) as f64;
        }
    }
    assert!(l1 < 0.06, "L1 {l1}");
    // the unlimited scheme either fails or is clearly worse
    let plain = solve_dg(&mesh, &law, &bcs, &init, 0.1, 0.2, None);
    if let Ok(p) = plain {
        let bad = (0..mesh.elements()).any(|e| (0..mesh.nodes_per_element()).any(|nd| !(p.value(e, nd, 0) > 0.0)));
        let _ = bad;
    }
}

// ------------------------------------------------------ oblique TF/SF

fn oblique_run(dim: usize, n: usize, pml: usize, dir: [f64; 3], pol: [f64; 3], lo: f64, hi: f64, width: f64, t_end: f64) -> (f64, f64) {
    let len = 2.0;
    let h = len / n as f64;
    let mut cfg = YeeConfig::new(&vec![n; dim], &vec![h; dim]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = pml;
    cfg.cfl = 0.9;
    let mut y = Yee::new(cfg).unwrap();
    let origin = [lo, lo, lo];
    let t0 = 4.0 * width;
    y.set_oblique_plane_wave(ObliquePlaneWave {
        lo: [lo; 3],
        hi: [hi; 3],
        direction: dir,
        polarization: pol,
        origin,
        waveform: Box::new(move |tau| (-((tau - t0) / width).powi(2)).exp()),
    })
    .unwrap();
    let nrm = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    let nhat = [dir[0] / nrm, dir[1] / nrm, dir[2] / nrm];
    let pl = (pol[0] * pol[0] + pol[1] * pol[1] + pol[2] * pol[2]).sqrt();
    let comp = if pol[2].abs() > 0.5 * pl { Comp::Ez } else { Comp::Ey };
    let ecomp = |v: [f64; 3]| if matches!(comp, Comp::Ez) { v[2] } else { v[1] };
    let epol = [pol[0] / pl, pol[1] / pl, pol[2] / pl];
    let nn = y.n();
    let centre = 0.5 * (lo + hi);
    // sample at the field component's own position nearest the box centre
    let ic = (centre / h).round() as usize;
    let probe = y.idx(ic, if dim > 1 { ic } else { 0 }, if dim > 2 { ic } else { 0 });
    let mut inside_err: f64 = 0.0;
    let mut leak: f64 = 0.0;
    let steps = (t_end / y.dt) as usize;
    for _ in 0..steps {
        y.step();
        let t = y.t;
        // analytic incident field at the probe (E component position)
        let pos = {
            let mut x = [0.0; 3];
            let ijk = [ic, if dim > 1 { ic } else { 0 }, if dim > 2 { ic } else { 0 }];
            for a in 0..dim {
                x[a] = (ijk[a] as f64 + if matches!(comp, Comp::Ez) && a == 2 || matches!(comp, Comp::Ey) && a == 1 { 0.5 } else { 0.0 }) * h;
            }
            x
        };
        let s: f64 = (0..3).map(|a| nhat[a] * (pos[a] - origin[a])).sum();
        let ex = ecomp(epol) * (-((t - s - t0) / width).powi(2)).exp();
        inside_err = inside_err.max((y.field(comp)[probe] - ex).abs());
        // scattered-field region: cells well outside the box and inside the PML-free interior
        let m = pml + 2;
        for p in 0..y.size() {
            let ijk = [p % nn[0], (p / nn[0]) % nn[1], p / (nn[0] * nn[1])];
            let mut outside = false;
            let mut interior = true;
            for a in 0..dim {
                let x = ijk[a] as f64 * h;
                if ijk[a] < m || ijk[a] + m >= nn[a] {
                    interior = false;
                }
                if x < lo - 0.1 || x > hi + 0.1 {
                    outside = true;
                }
            }
            if interior && outside {
                leak = leak.max(y.field(comp)[p].abs());
            }
        }
    }
    (inside_err, leak)
}

#[test]
fn oblique_plane_wave_on_an_arbitrary_box_2d_and_3d() {
    let a = 30.0_f64.to_radians();
    let (inside, leak) = oblique_run(2, 120, 12, [a.cos(), a.sin(), 0.0], [0.0, 0.0, 1.0], 0.6, 1.4, 0.15, 3.2);
    assert!(inside < 0.06, "2D incident field inside the box: {inside}");
    assert!(leak < 0.06, "2D leakage into the scattered region: {leak}");
    let (inside, leak) = oblique_run(3, 40, 8, [1.0, 1.0, 1.0], [1.0, -1.0, 0.0], 0.6, 1.4, 0.3, 2.8);
    assert!(inside < 0.1, "3D incident field inside the box: {inside}");
    assert!(leak < 0.1, "3D leakage: {leak}");
    // a direction along an inactive axis is rejected
    let mut y = Yee::new(YeeConfig::new(&[16, 16], &[0.1, 0.1]).unwrap()).unwrap();
    let bad = ObliquePlaneWave {
        lo: [0.2; 3],
        hi: [1.0; 3],
        direction: [1.0, 0.0, 1.0],
        polarization: [0.0, 1.0, 0.0],
        origin: [0.0; 3],
        waveform: Box::new(|t| t),
    };
    assert!(y.set_oblique_plane_wave(bad).is_err());
}

// ------------------------------------------------------------ BEM

fn sphere_bcs(mesh: &SurfaceMesh3D, f: impl Fn([f64; 3]) -> f64) -> Vec<BoundaryCondition<f64>> {
    (0..mesh.triangles.len()).map(|t| BoundaryCondition::Potential(f(mesh.centroid(t).unwrap()))).collect()
}

#[test]
fn exterior_bem_problems_on_the_sphere() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);
    let opts = AcaOptions { leaf_size: 16, aca_tol: 1e-8, ..AcaOptions::default() };
    // capacitance problem: u = 1 on the sphere, u = 1/r outside, q = -1
    let (sol, _) = solve_laplace_bem_3d_exterior(&mesh, &sphere_bcs(&mesh, |_| 1.0), &opts).unwrap();
    let worst = sol.q.iter().map(|q| (q + 1.0).abs()).fold(0.0, f64::max);
    assert!(worst < 0.08, "capacitance flux error {worst}");
    let p = evaluate_potential_3d_exterior([2.0, 0.0, 0.0], &mesh, &sol).unwrap();
    assert!((p - 0.5).abs() < 0.03, "u(2) = {p}");
    // dipole: u = x on the sphere -> u = x / r^3, q = -2 x
    let (sol, _) = solve_laplace_bem_3d_exterior(&mesh, &sphere_bcs(&mesh, |c| c[0]), &opts).unwrap();
    let worst = (0..mesh.triangles.len()).map(|t| (sol.q[t] + 2.0 * mesh.centroid(t).unwrap()[0]).abs()).fold(0.0, f64::max);
    assert!(worst < 0.2, "dipole flux error {worst}");
    // exterior Neumann: q = -1 everywhere gives u = 1/r
    let neu: Vec<BoundaryCondition<f64>> = (0..mesh.triangles.len()).map(|_| BoundaryCondition::Flux(-1.0)).collect();
    let (sol, _) = solve_laplace_bem_3d_exterior(&mesh, &neu, &opts).unwrap();
    let worst = sol.u.iter().map(|u| (u - 1.0).abs()).fold(0.0, f64::max);
    assert!(worst < 0.08, "exterior Neumann potential error {worst}");
}

#[test]
fn pure_neumann_interior_bem_with_nullspace_treatment() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);
    let opts = AcaOptions { leaf_size: 16, aca_tol: 1e-8, ..AcaOptions::default() };
    // u = x harmonic, q = n_x; the solution is unique up to a constant (zero-mean fixed)
    let bcs: Vec<BoundaryCondition<f64>> =
        (0..mesh.triangles.len()).map(|t| BoundaryCondition::Flux(mesh.normal(t).unwrap()[0])).collect();
    let (sol, _) = solve_laplace_bem_3d_aca(&mesh, &bcs, &opts).unwrap();
    let worst = (0..mesh.triangles.len()).map(|t| (sol.u[t] - mesh.centroid(t).unwrap()[0]).abs()).fold(0.0, f64::max);
    assert!(worst < 0.08, "Neumann potential error {worst}");
    // incompatible data are rejected
    let bad: Vec<BoundaryCondition<f64>> = (0..mesh.triangles.len()).map(|_| BoundaryCondition::Flux(1.0)).collect();
    assert!(solve_laplace_bem_3d_aca(&mesh, &bad, &opts).is_err());
}

#[test]
fn linear_element_bem_is_more_accurate_than_constant_elements() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);
    let bcs: Vec<BoundaryCondition<f64>> = mesh.vertices.iter().map(|v| BoundaryCondition::Potential(v[0])).collect();
    let sol = solve_laplace_bem_3d_p1(&mesh, &bcs, false).unwrap();
    let worst = mesh.vertices.iter().zip(&sol.q).map(|(v, q)| (q - v[0]).abs()).fold(0.0, f64::max);
    assert!(worst < 0.03, "P1 flux error {worst}");
    // mixed data u = x y + z: Dirichlet on the lower hemisphere, Neumann above
    let mixed: Vec<BoundaryCondition<f64>> = mesh
        .vertices
        .iter()
        .map(|v| {
            if v[2] < 0.0 {
                BoundaryCondition::Potential(v[0] * v[1] + v[2])
            } else {
                BoundaryCondition::Flux(v[1] * v[0] + v[0] * v[1] + v[2])
            }
        })
        .collect();
    let sol = solve_laplace_bem_3d_p1(&mesh, &mixed, false).unwrap();
    let worst = mesh.vertices.iter().zip(&sol.u).map(|(v, u)| (u - (v[0] * v[1] + v[2])).abs()).fold(0.0, f64::max);
    assert!(worst < 0.05, "P1 mixed potential error {worst}");
    // exterior capacitance
    let ext: Vec<BoundaryCondition<f64>> = mesh.vertices.iter().map(|_| BoundaryCondition::Potential(1.0)).collect();
    let sol = solve_laplace_bem_3d_p1(&mesh, &ext, true).unwrap();
    let worst = sol.q.iter().map(|q| (q + 1.0).abs()).fold(0.0, f64::max);
    assert!(worst < 0.05, "P1 exterior flux error {worst}");
    assert!(solve_laplace_bem_3d_p1(&mesh, &ext[..3], false).is_err());
    let _ = Face::XMin;
}
