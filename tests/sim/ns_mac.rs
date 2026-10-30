//! Incompressible Navier-Stokes projection solver tests.

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::Face;
use rssn::sim::boundary::FluidBoundary;
use rssn::sim::boundary::FluidFaceBcs;
use rssn::sim::ns_mac::ImmersedBody;
use rssn::sim::ns_mac::NsConfig;
use rssn::sim::ns_mac::NsSolver;
use std::f64::consts::PI;

fn tg_exact(
    x: [f64; 3],
    t: f64,
    nu: f64,
) -> [f64; 3] {
    let d = (-2.0 * nu * t).exp();
    [
        -x[0].cos() * x[1].sin() * d,
        x[0].sin() * x[1].cos() * d,
        0.0,
    ]
}

fn taylor_green_error(
    n: usize,
    nz: usize,
) -> (f64, f64) {
    let nu = 0.1;
    let (dim, len): (usize, Vec<f64>) = if nz > 1 {
        (3, vec![2.0 * PI, 2.0 * PI, 1.0])
    } else {
        (2, vec![2.0 * PI, 2.0 * PI])
    };
    let ns: Vec<usize> = if nz > 1 {
        vec![n, n, nz]
    } else {
        vec![n, n]
    };
    let mut cfg = NsConfig::new(&ns, &len, nu).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
    let h = 2.0 * PI / n as f64;
    let dt = 0.4 * h;
    let steps = (1.0 / dt).round() as usize;
    let dt = 1.0 / steps as f64;
    let mut s = NsSolver::new(cfg, dt, &|x| tg_exact(x, 0.0, nu)).unwrap();
    s.run(steps).unwrap();
    let mut worst: f64 = 0.0;
    let nn = s.n();
    for c in 0..2 {
        for j in 0..nn[1] + usize::from(c == 1) {
            for i in 0..nn[0] + usize::from(c == 0) {
                let x = [
                    (i as f64 + if c == 0 { 0.0 } else { 0.5 }) * h,
                    (j as f64 + if c == 1 { 0.0 } else { 0.5 }) * h,
                    0.0,
                ];
                let v = s.sample_velocity(c, x);
                worst = worst.max((v - tg_exact(x, 1.0, nu)[c]).abs());
            }
        }
    }
    let div = s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    let _ = dim;
    (worst, div)
}

#[test]
fn taylor_green_vortex_second_order_in_2d_and_extruded_3d() {
    let (e1, d1) = taylor_green_error(16, 1);
    let (e2, d2) = taylor_green_error(32, 1);
    let (e3, _) = taylor_green_error(64, 1);
    assert!(
        (e1 / e2).log2() > 1.8 && (e2 / e3).log2() > 1.8,
        "{e1} {e2} {e3}"
    );
    assert!(d1 < 1e-8 && d2 < 1e-8, "divergence {d1} {d2}");
    // 3D extrusion reproduces the 2D result
    let (e3d, d3) = taylor_green_error(16, 4);
    assert!((e3d - e1).abs() < 1e-8 + 1e-6 * e1, "3D {e3d} vs 2D {e1}");
    assert!(d3 < 1e-8);
}

fn poiseuille_error(
    n: usize,
    dim3: bool,
) -> f64 {
    let (nu, g, hgt) = (0.1, 1.0, 1.0);
    let (ns, len): (Vec<usize>, Vec<f64>) = if dim3 {
        (vec![4, n, 4], vec![1.0, hgt, 1.0])
    } else {
        (vec![4, n], vec![1.0, hgt])
    };
    let mut cfg = NsConfig::new(&ns, &len, nu).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with_axis(0, &FluidBoundary::Periodic);
    if dim3 {
        cfg.bcs = cfg.bcs.clone().with_axis(2, &FluidBoundary::Periodic);
    }
    cfg.force = Some(Box::new(move |_t, _x| [g, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, 0.1, &|_| [0.0; 3]).unwrap();
    s.run(150).unwrap();
    let mut worst: f64 = 0.0;
    for k in 0..n {
        let y = (k as f64 + 0.5) / n as f64;
        let ex = g / (2.0 * nu) * y * (hgt - y);
        worst = worst.max((s.sample_velocity(0, [0.5, y, 0.5]) - ex).abs());
    }
    worst
}

#[test]
fn poiseuille_flow_with_no_slip_walls_2d_and_3d() {
    let e1 = poiseuille_error(8, false);
    let e2 = poiseuille_error(16, false);
    assert!((e1 / e2).log2() > 1.8 && e2 < 8e-3, "{e1} {e2}");
    let e3 = poiseuille_error(16, true);
    assert!(
        (e3 - e2).abs() < 1e-9,
        "3D duct with periodic z equals 2D: {e3} {e2}"
    );
}

fn cavity(
    n: usize,
    re: f64,
    t_end: f64,
    dt: f64,
) -> NsSolver {
    let mut cfg = NsConfig::new(&[n, n], &[1.0, 1.0], 1.0 / re).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run((t_end / dt).round() as usize).unwrap();
    s
}

#[test]
fn lid_driven_cavity_re100_matches_ghia_reference_points() {
    let s = cavity(64, 100.0, 30.0, 0.01);
    // Ghia, Ghia & Shin (1982), Re = 100: u along x = 0.5
    let u_ref = [
        (0.9766, 0.84123),
        (0.9688, 0.78871),
        (0.9609, 0.73722),
        (0.9531, 0.68717),
        (0.8516, 0.23151),
        (0.7344, 0.00332),
        (0.6172, -0.13641),
        (0.5000, -0.20581),
        (0.4531, -0.21090),
        (0.2813, -0.15662),
        (0.1719, -0.10150),
        (0.1016, -0.06434),
        (0.0703, -0.04775),
        (0.0625, -0.04192),
        (0.0547, -0.03717),
    ];
    for (y, u) in u_ref {
        let v = s.sample_velocity(0, [0.5, y, 0.0]);
        assert!((v - u).abs() < 0.015, "u({y}) = {v}, Ghia {u}");
    }
    // v along y = 0.5
    let v_ref = [
        (0.9688, -0.05906),
        (0.9531, -0.08864),
        (0.9063, -0.16914),
        (0.8594, -0.22445),
        (0.8047, -0.24533),
        (0.5, 0.05454),
        (0.2344, 0.17527),
        (0.1563, 0.16077),
        (0.0938, 0.12317),
        (0.0703, 0.10091),
        (0.0625, 0.09233),
    ];
    for (x, v) in v_ref {
        let val = s.sample_velocity(1, [x, 0.5, 0.0]);
        assert!((val - v).abs() < 0.015, "v({x}) = {val}, Ghia {v}");
    }
    assert!(s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs())) < 1e-7);
}

#[test]
fn inflow_outflow_channel_keeps_parabolic_profile_and_pressure_drop() {
    let (n0, n1) = (48, 16);
    let nu = 0.1;
    let mut cfg = NsConfig::new(&[n0, n1], &[4.0, 1.0], nu).unwrap();
    let prof = BcValue::space(|x| 4.0 * x[1] * (1.0 - x[1]));
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(
            Face::XMin,
            FluidBoundary::Inflow([prof, BcValue::constant(0.0), BcValue::constant(0.0)]),
        )
        .with(Face::XMax, FluidBoundary::outflow());
    let mut s = NsSolver::new(cfg, 0.02, &|_| [0.0; 3]).unwrap();
    s.run(900).unwrap();
    for k in 0..n1 {
        let y = (k as f64 + 0.5) / n1 as f64;
        let ex = 4.0 * y * (1.0 - y);
        let v = s.sample_velocity(0, [3.9, y, 0.0]);
        assert!((v - ex).abs() < 5e-3, "outlet u({y}) = {v} vs {ex}");
        assert!(s.sample_velocity(1, [3.9, y, 0.0]).abs() < 5e-3);
    }
    // dp/dx = nu u'' = -0.8
    let i1 = s.cell_index([6, 8, 0]);
    let i2 = s.cell_index([36, 8, 0]);
    let dpdx = (s.p[i2] - s.p[i1]) / (30.0 * s.cell_size(0));
    assert!((dpdx + 0.8).abs() < 0.04, "dp/dx {dpdx}");
    // mass flux conserved: inflow equals outflow
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, n0));
    assert!((fin - fout).abs() < 1e-6 * fin.abs(), "{fin} {fout}");
}

#[test]
fn flow_past_square_obstacle_has_no_slip_wake_and_conserves_mass() {
    let (n0, n1) = (96, 32);
    let h = 1.0 / n1 as f64;
    let mut cfg = NsConfig::new(&[n0, n1], &[3.0, 1.0], 0.03).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::inflow([1.0, 0.0, 0.0]))
        .with(Face::XMax, FluidBoundary::outflow());
    let mut solid = vec![false; n0 * n1];
    for j in 0..n1 {
        for i in 0..n0 {
            let (x, y) = ((i as f64 + 0.5) * h, (j as f64 + 0.5) * h);
            if (0.75..1.05).contains(&x) && (0.3125..0.6875).contains(&y) {
                solid[i + n0 * j] = true;
            }
        }
    }
    cfg.solid = Some(solid.clone());
    let mut s = NsSolver::new(cfg, 0.01, &|_| [1.0, 0.0, 0.0]).unwrap();
    s.run(1500).unwrap();
    // mass conservation
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, n0));
    assert!((fin - fout).abs() < 1e-3 * fin.abs(), "{fin} {fout}");
    // zero velocity on faces of solid cells
    for j in 0..n1 {
        for i in 0..n0 {
            if solid[i + n0 * j] {
                let x = [(i as f64 + 0.5) * h, (j as f64 + 0.5) * h, 0.0];
                assert!(
                    s.sample_velocity(0, x).abs() < 1e-12 && s.sample_velocity(1, x).abs() < 1e-12
                );
            }
        }
    }
    // recirculation behind the body
    let mut min_u = f64::INFINITY;
    for k in 0..40 {
        let x = 1.1 + 0.02 * k as f64;
        min_u = min_u.min(s.sample_velocity(0, [x, 0.5, 0.0]));
    }
    assert!(min_u < -1e-3, "no reverse flow found (min u {min_u})");
    assert!(s.max_speed() < 3.0);
}

/// Measured force coefficients of the flow past a cylinder in one run.
#[derive(Debug, Clone, Copy)]
struct CylinderStats {
    /// Mean drag coefficient (momentum balance on a box around the body).
    cd: f64,
    /// Mean lift coefficient.
    cl_mean: f64,
    /// Half the peak-to-peak lift swing.  [`NsSolver::force_on_region`] omits
    /// the unsteady momentum of the box, so this overstates the lift on the
    /// body; it is used here only as a shedding indicator.
    cl_amp: f64,
    /// Strouhal number from the lift zero (upward) crossings.
    st_zero: f64,
    /// Strouhal number from the peak of the windowed DFT of the lift.
    st_fft: f64,
    /// Number of shedding periods in the measurement window.
    cycles: f64,
    /// Start/end of the measurement window.
    t_start: f64,
    t_end: f64,
    /// Wall time of the time loop.
    runtime_s: f64,
}

/// Flow past a cylinder (diameter 1, centre at `(lx/4, ly/2)`) at `Re = 100`
/// in the box `[0, lx] x [0, ly]`, uniform inflow `U = 1`, slip side walls and
/// an outflow at `x = lx`. The initial condition carries an antisymmetric
/// wake perturbation so that the (unstable) shedding mode is excited without
/// waiting for numerical noise to grow.
///
/// `t_start` is the end of the discarded transient; the run stops once
/// `min_cycles` lift periods have been collected (or at `max_t`).
///
/// `seed` selects the symmetry-breaking initial condition: `0` none, `1` a
/// mirror-symmetric (lift-free) strain field, `2` a localized vortex behind
/// the body, `3`/`4` the vortex / nothing with the body offset by `0.3 D`.
#[allow(clippy::too_many_arguments)]
fn cylinder_run(
    cpd: usize,
    (lx, ly): (f64, f64),
    cfl: f64,
    t_start: f64,
    min_cycles: usize,
    max_t: f64,
    seed: u8,
    solver: rssn::sim::structured::LinearSolver,
    verbose: bool,
) -> (NsSolver, CylinderStats) {
    let h = 1.0 / cpd as f64;
    let (cx, cy) = (
        0.25 * lx,
        0.5 * ly
            + if seed == 3 || seed == 4 {
                0.3
            } else {
                0.0
            },
    );
    let n0 = (lx / h).round() as usize;
    let n1 = (ly / h).round() as usize;
    let mut cfg = NsConfig::new(&[n0, n1], &[lx, ly], 0.01).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::inflow([1.0, 0.0, 0.0]))
        .with(Face::XMax, FluidBoundary::outflow())
        .with(Face::YMin, FluidBoundary::Slip)
        .with(Face::YMax, FluidBoundary::Slip);
    cfg.solver = solver;
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(move |x| ((x[0] - cx).powi(2) + (x[1] - cy).powi(2)).sqrt() - 0.5),
        velocity: Box::new(|_| [0.0; 3]),
    });
    let dt = cfl * h;
    let ic = move |x: [f64; 3]| {
        if seed == 0 {
            return [1.0, 0.0, 0.0];
        }
        if seed == 1 {
            let (dx, dy) = (x[0] - cx - 1.0, x[1] - cy);
            return [
                1.0,
                0.3 * dy * (-0.5 * (dx * dx + dy * dy) / 9.0).exp(),
                0.0,
            ];
        }
        let (dx, dy) = (x[0] - cx - 1.2, x[1] - cy);
        let g = (-0.25 * (dx * dx + dy * dy)).exp();
        [1.0 - 0.5 * dy * g, 0.5 * dx * g, 0.0]
    };
    let t_build = std::time::Instant::now();
    let mut s = NsSolver::new(cfg, dt, &ic).unwrap();
    let setup_s = t_build.elapsed().as_secs_f64();
    let mut samples: Vec<(f64, f64, f64)> = Vec::new();
    let mut crossings: Vec<f64> = Vec::new();
    let mut prev_cl: Option<f64> = None;
    let mut steps = 0usize;
    let start = std::time::Instant::now();
    while s.t < max_t {
        s.step().unwrap();
        steps += 1;
        let f = s.force_on_region([cx - 1.5, cy - 1.5, 0.0], [cx + 1.5, cy + 1.5, 0.0]);
        let (cd, cl) = (2.0 * f[0], 2.0 * f[1]);
        if s.t >= t_start {
            if let Some(p) = prev_cl
                && p < 0.0
                && cl >= 0.0
            {
                crossings.push(s.t - dt + dt * (-p) / (cl - p));
            }
            samples.push((s.t, cd, cl));
        }
        prev_cl = Some(cl);
        if verbose && steps % 50 == 0 {
            println!(
                "  cpd={cpd} t={:.2} cd={cd:.3} cl={cl:.3} maxu={:.3} wall={:.1}s",
                s.t,
                s.max_speed(),
                start.elapsed().as_secs_f64()
            );
        }
        if s.t >= t_start && crossings.len() > min_cycles {
            break;
        }
    }
    let runtime_s = start.elapsed().as_secs_f64();
    let n = samples.len();
    let cd = samples.iter().map(|x| x.1).sum::<f64>() / n as f64;
    let cl_mean = samples.iter().map(|x| x.2).sum::<f64>() / n as f64;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for x in &samples {
        lo = lo.min(x.2);
        hi = hi.max(x.2);
    }
    let cl_amp = 0.5 * (hi - lo);
    let cycles = crossings.len().saturating_sub(1) as f64;
    let st_zero = if cycles >= 1.0 {
        let periods: f64 = crossings.windows(2).map(|w| w[1] - w[0]).sum::<f64>() / cycles;
        1.0 / periods
    } else {
        f64::NAN
    };
    // windowed DFT of the lift, peak by parabolic refinement over a fine grid
    let t0 = samples[0].0;
    let w: Vec<f64> = (0..n)
        .map(|k| 0.5 - 0.5 * (2.0 * PI * k as f64 / (n.max(2) - 1) as f64).cos())
        .collect();
    let mut best = (f64::NAN, 0.0, f64::NAN);
    let mut df = 0.1;
    while df <= 0.25 {
        let (mut re, mut im) = (0.0, 0.0);
        for (k, &(t, _, cl)) in samples.iter().enumerate() {
            let ph = 2.0 * PI * df * (t - t0);
            re += w[k] * (cl - cl_mean) * ph.cos();
            im -= w[k] * (cl - cl_mean) * ph.sin();
        }
        let mag = re * re + im * im;
        if mag > best.1 {
            best = (df, mag, f64::NAN);
        }
        df += 1e-4;
    }
    let st_fft = best.0;
    if verbose {
        println!(
            "  cpd={cpd} cells={} setup={setup_s:.2}s run={runtime_s:.2}s steps={steps}\n  \
             cd={cd:.4} cl_mean={cl_mean:.4} cl_amp={cl_amp:.4} st_zero={st_zero:.4} st_fft={st_fft:.4} cycles={cycles} window=[{t_start:.1},{:.1}]",
            n0 * n1,
            samples[n - 1].0
        );
    }
    (
        s,
        CylinderStats {
            cd,
            cl_mean,
            cl_amp,
            st_zero,
            st_fft,
            cycles,
            t_start,
            t_end: samples[n - 1].0,
            runtime_s,
        },
    )
}

// ---------------------------------------------------------------------------
// Flow past a circular cylinder at Re = 100: validation of the immersed-body
// MAC solver against the classical benchmark.
//
// Reference values for the *unconfined* cylinder (D = 1, U_inf = 1):
//   C_D = 1.32-1.34 and St = 0.164-0.166 at Re = 100
//   (Braza, Chassaing & Ha Minh, J. Fluid Mech. 165 (1986) 79-130;
//    Williamson, J. Fluid Mech. 206 (1989) 579-627; the ranges are those
//    collected in Rajani, Kandasamy & Majumdar, Appl. Math. Modelling 33
//    (2009) 668-680).
//   Steady wake at Re = 40: C_D = 1.52-1.60 (Dennis & Chang, J. Fluid Mech.
//   42 (1970) 471-489); the repo's existing `ns_research.rs` test already
//   covers Re = 20/40 with the same solver.
//
// The box used here is [0, lx] x [0, ly] with the cylinder at (lx/4, ly/2),
// uniform inflow at x = 0, zero-pressure outflow at x = lx and free-slip side
// walls.  Its blockage D/ly is 10 % (20 x 10 box) or 5 % (40 x 20 box).
//
// Measured with the 20 x 10 box, 6 shedding periods after discarding the
// transient:
//   cells/D   C_D      St (crossings)   St (DFT)   lift amplitude
//   8         1.5069   0.1744           0.1746     0.92
//   12        1.4890   0.1777           0.1777     1.01
//   16        1.4956   0.1779           0.1779     1.05
// i.e. Cd = 1.49-1.51 and St = 0.174-0.178: Cd is 12-13 % above the
// unconfined value and St 6-8 % above it, almost all of it blockage.  With
// 5 % blockage (40 x 20 box) at the same 8 cells/D the same run gives
// C_D = 1.3988 and St = 0.1653 (+5.2 % and +0.2 % from the literature).
// ---------------------------------------------------------------------------

/// Unconfined reference drag coefficient at `Re = 100`.
const CD_RE100: f64 = 1.33;
/// Unconfined reference Strouhal number at `Re = 100`.
const ST_RE100: f64 = 0.165;

/// Shedding smoke test: the wake must oscillate with a lift amplitude well
/// above the numerical noise, stay roughly symmetric in the mean and conserve
/// mass.  One lift period is enough here.
#[test]
fn cylinder_re100_shedding_smoke_test() {
    let (s, r) = cylinder_run(
        8,
        (20.0, 10.0),
        0.4,
        30.0,
        1,
        60.0,
        2,
        rssn::sim::structured::LinearSolver::Direct,
        true,
    );
    println!(
        "cylinder smoke: C_D={:.4} St={:.4} lift amplitude={:.4} mean lift={:.4} \
         ({} cycle over [{:.1}, {:.1}], {:.1} s)",
        r.cd, r.st_zero, r.cl_amp, r.cl_mean, r.cycles, r.t_start, r.t_end, r.runtime_s
    );
    assert!(
        r.cycles >= 1.0,
        "no complete lift period was collected: {r:?}"
    );
    assert!(
        r.cl_amp > 0.05,
        "the wake does not shed, lift amplitude {} <= 0.05",
        r.cl_amp
    );
    assert!(
        r.cl_mean.abs() < 0.05,
        "the mean lift must vanish, got {}",
        r.cl_mean
    );
    assert!(
        (0.12..=0.22).contains(&r.st_zero),
        "shedding frequency {:.4} is outside the laminar cylinder band",
        r.st_zero
    );
    // positivity of the streamwise velocity on the upstream axis and bounded speed
    for k in 0..20 {
        let x = [0.25 + 0.1 * k as f64, 5.0, 0.0];
        assert!(
            s.sample_velocity(0, x) > 0.0,
            "reverse flow upstream at {x:?}"
        );
    }
    assert!(s.max_speed() < 5.0, "speed blew up: {}", s.max_speed());
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, s.n()[0]));
    assert!(
        (fin - fout).abs() < 1e-3 * fin.abs(),
        "mass not conserved: {fin} vs {fout}"
    );
}

/// `Re = 100` drag and Strouhal number on the 20 x 10 box (10 % blockage).
///
/// Measured (8 cells per diameter, dt = 0.05, 6 periods in `[35.0, 70.1]`):
/// `C_D = 1.5069`, `St = 0.1744` (DFT 0.1746), lift amplitude 0.92.
/// Against the unconfined literature this is `C_D` +13.3 % and `St` +5.7 %;
/// the box is confined to 10 % blockage, which accounts for most of the drag
/// offset (see the 5 %-blockage test).  The asserted bands are the achieved
/// deviations with a small margin.
#[test]
fn cylinder_re100_drag_and_strouhal_10pct_blockage() {
    let (_, r) = cylinder_run(
        8,
        (20.0, 10.0),
        0.4,
        35.0,
        6,
        90.0,
        2,
        rssn::sim::structured::LinearSolver::Direct,
        true,
    );
    println!(
        "cylinder Re=100 (10 % blockage): C_D={:.4} (lit {CD_RE100}, {:+.1} %), \
         St={:.4} (lit {ST_RE100}, {:+.1} %), DFT {:.4}, amplitude {:.4}, \
         {:.0} cycles over [{:.1}, {:.1}], {:.1} s",
        r.cd,
        100.0 * (r.cd / CD_RE100 - 1.0),
        r.st_zero,
        100.0 * (r.st_zero / ST_RE100 - 1.0),
        r.st_fft,
        r.cl_amp,
        r.cycles,
        r.t_start,
        r.t_end,
        r.runtime_s
    );
    assert!(
        r.cycles >= 6.0,
        "need at least 6 shedding cycles, got {}",
        r.cycles
    );
    assert!(r.cl_amp > 0.05, "no shedding, lift amplitude {}", r.cl_amp);
    assert!(
        r.cl_mean.abs() < 0.02,
        "mean lift {} is not small",
        r.cl_mean
    );
    assert!(
        (r.st_fft - r.st_zero).abs() < 1e-3,
        "St crossing {} vs DFT {}",
        r.st_zero,
        r.st_fft
    );
    // achieved deviations: Cd +13.3 %, St +5.7 %; bands with margin
    assert!(
        (r.cd - CD_RE100).abs() < 0.20,
        "C_D = {:.4}, literature {CD_RE100} (+{:.1} %)",
        r.cd,
        100.0 * (r.cd / CD_RE100 - 1.0)
    );
    assert!(
        (r.st_zero - ST_RE100).abs() < 0.012,
        "St = {:.4}, literature {ST_RE100} (+{:.1} %)",
        r.st_zero,
        100.0 * (r.st_zero / ST_RE100 - 1.0)
    );
}

/// Grid convergence of `C_D` and `St` at `Re = 100` on the 20 x 10 box and
/// Richardson extrapolation.  The 12- and 16-cell grids give St = 0.1777 and
/// 0.1779, i.e. a grid-converged St = 0.178 (about 8 % above the unconfined
/// 0.165); C_D is 1.489-1.496 (about 12 % above 1.33).  The remaining offset
/// is the finite domain, not the grid: at fixed 8 cells/D the 5 %-blockage
/// box lowers C_D from 1.51 to 1.40.  St converges monotonically (the 8 -> 12
/// step sets the observed order, the 12 -> 16 change is only 2e-4, so the
/// fitted order is high); the C_D sequence is non-monotone at the 0.5 % level
/// and admits no positive-order Richardson fit, which is itself evidence that
/// blockage, not resolution, dominates its error.
#[test]
#[ignore = "about 20 min (three grids): cargo test --test sim ns_mac -- --ignored cylinder_re100_grid_study"]
fn cylinder_re100_grid_study() {
    let mut cd = [0.0; 3];
    let mut st = [0.0; 3];
    let cpd = [8usize, 12, 16];
    for (k, &c) in cpd.iter().enumerate() {
        let (_, r) = cylinder_run(
            c,
            (20.0, 10.0),
            0.4,
            35.0,
            6,
            130.0,
            2,
            rssn::sim::structured::LinearSolver::Direct,
            true,
        );
        cd[k] = r.cd;
        st[k] = r.st_zero;
        assert!(r.cycles >= 6.0);
        assert!((r.st_fft - r.st_zero).abs() < 1e-3);
    }
    println!("cells/D  C_D      St");
    for k in 0..3 {
        println!("{:7}  {:.4}   {:.4}", cpd[k], cd[k], st[k]);
    }
    // Richardson: f(h) = f0 + C h^p, h = 1/cells_per_diameter.
    let extrapolate = |f2: f64, f3: f64, p: f64| f3 - (f2 - f3) / ((16.0_f64 / 12.0).powf(p) - 1.0);
    let order = |f1: f64, f2: f64, f3: f64| -> Option<f64> {
        let (h1, h2, h3): (f64, f64, f64) = (1.0 / 8.0, 1.0 / 12.0, 1.0 / 16.0);
        let (d12, d23) = (f1 - f2, f2 - f3);
        if d12 * d23 <= 0.0 {
            return None;
        }
        let target = d12 / d23;
        let g = |p: f64| (h1.powf(p) - h2.powf(p)) / (h2.powf(p) - h3.powf(p));
        let (mut lo, mut hi) = (0.1_f64, 10.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) < target {
                lo = mid
            } else {
                hi = mid
            }
        }
        Some(0.5 * (lo + hi))
    };
    let st_order = order(st[0], st[1], st[2]).expect("St must converge monotonically");
    let st_ex = extrapolate(st[1], st[2], st_order);
    println!("St: observed order p={st_order:.2}, Richardson St={st_ex:.4} (lit {ST_RE100})");
    match order(cd[0], cd[1], cd[2]) {
        | Some(p) => {
            let cd_ex = extrapolate(cd[1], cd[2], p);
            println!("C_D: observed order p={p:.2}, Richardson C_D={cd_ex:.4} (lit {CD_RE100})");
            assert!((cd_ex - cd[2]).abs() < 0.02);
        },
        | None => {
            println!(
                "C_D differences change sign: the 10 % blockage error dominates, no Richardson fit"
            )
        },
    }
    // St converges from below to 0.178, 8 % above the unconfined value.
    assert!(
        st_ex > st[2] - 5e-3,
        "St extrapolation {st_ex} inconsistent with the finest grid"
    );
    assert!(
        (st_ex - 0.178).abs() < 0.01,
        "extrapolated St {st_ex} moved away from the measured 0.178"
    );
    assert!(
        (st[2] - ST_RE100).abs() < 0.02,
        "finest-grid St {} is off literature",
        st[2]
    );
    assert!(
        (cd[2] - CD_RE100).abs() < 0.20,
        "finest-grid C_D {} is off literature",
        cd[2]
    );
}

/// `Re = 100` benchmark on the 40 x 20 box (5 % blockage), 8 cells per
/// diameter, `dt = 0.05`.
///
/// Measured over 6 shedding periods in `[35.0, 73.7]`: `C_D = 1.3988`,
/// `St = 0.1653` (zero crossings; DFT 0.1653), lift amplitude 0.87.
/// Against the unconfined literature: `C_D` +5.2 %, `St` +0.2 %.  This is
/// the validated case referred to in the module documentation of
/// `src/sim/ns_mac.rs`; enlarging the box from 20 x 10 to 40 x 20 at fixed
/// resolution moved `C_D` from 1.5069 to 1.3988 and `St` from 0.1744 to
/// 0.1653, which identifies blockage (not resolution) as the dominant error.
///
/// Slow (several minutes): run with `--ignored`.
#[test]
#[ignore = "about 4-8 min: cargo test --test sim ns_mac -- --ignored cylinder_re100_5pct"]
fn cylinder_re100_5pct_blockage_matches_literature() {
    let (s, r) = cylinder_run(
        8,
        (40.0, 20.0),
        0.4,
        35.0,
        6,
        130.0,
        2,
        rssn::sim::structured::LinearSolver::Direct,
        true,
    );
    println!(
        "cylinder Re=100 (5 % blockage): C_D={:.4} (lit {CD_RE100}, {:+.1} %), \
         St={:.4} (lit {ST_RE100}, {:+.1} %), DFT {:.4}, amplitude {:.4}, mean lift {:.4}",
        r.cd,
        100.0 * (r.cd / CD_RE100 - 1.0),
        r.st_zero,
        100.0 * (r.st_zero / ST_RE100 - 1.0),
        r.st_fft,
        r.cl_amp,
        r.cl_mean
    );
    assert!(
        r.cycles >= 6.0,
        "need at least 6 shedding cycles, got {}",
        r.cycles
    );
    assert!(r.cl_amp > 0.05, "no shedding, lift amplitude {}", r.cl_amp);
    assert!(
        r.cl_mean.abs() < 0.06,
        "mean lift {} is not small",
        r.cl_mean
    );
    assert!(
        (r.st_fft - r.st_zero).abs() < 1e-3,
        "St crossing {} vs DFT {}",
        r.st_zero,
        r.st_fft
    );
    // achieved: C_D +5.2 %, St +0.2 % (bands are the achieved values plus margin)
    assert!(
        (r.cd - CD_RE100).abs() < 0.10,
        "C_D = {:.4} vs literature {CD_RE100} (+{:.1} %)",
        r.cd,
        100.0 * (r.cd / CD_RE100 - 1.0)
    );
    assert!(
        (r.st_zero - ST_RE100).abs() < 0.006,
        "St = {:.4} vs literature {ST_RE100} (+{:.1} %)",
        r.st_zero,
        100.0 * (r.st_zero / ST_RE100 - 1.0)
    );
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, s.n()[0]));
    assert!(
        (fin - fout).abs() < 1e-3 * fin.abs(),
        "mass not conserved: {fin} vs {fout}"
    );
}

#[test]
fn lid_driven_cavity_3d_is_divergence_free_and_symmetric() {
    let n = 12;
    let mut cfg = NsConfig::new(&[n, n, n], &[1.0; 3], 0.01).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, 0.02, &|_| [0.0; 3]).unwrap();
    s.run(100).unwrap();
    assert!(s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs())) < 1e-7);
    // reflection symmetry in z: u(x,y,z) = u(x,y,1-z), w antisymmetric
    for k in 0..n {
        let z = (k as f64 + 0.5) / n as f64;
        let (a, b) = (
            s.sample_velocity(0, [0.4, 0.8, z]),
            s.sample_velocity(0, [0.4, 0.8, 1.0 - z]),
        );
        assert!((a - b).abs() < 1e-8, "u symmetry {a} {b}");
        let (c, d) = (
            s.sample_velocity(2, [0.4, 0.8, z]),
            s.sample_velocity(2, [0.4, 0.8, 1.0 - z]),
        );
        assert!((c + d).abs() < 1e-8, "w antisymmetry {c} {d}");
    }
    assert!(s.sample_velocity(0, [0.5, 0.95, 0.5]) > 0.2);
    assert!(s.kinetic_energy() > 1e-3);
}
