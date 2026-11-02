//! 3D DG shock limiters: order preservation, troubled-cell locality, Sod shock
//! tubes with total-variation / undershoot checks, conservation and a blast
//! stability regression.
//!
//! The Sod reference shock position at `t = 0.2` is
//! `x = 0.5 + 1.752_155_7 * 0.2 = 0.850_431`, from the exact self-similar
//! solution of the Sod problem with `γ = 1.4`.

use rssn::sim::dg::DgLaw;
use rssn::sim::dg3d::Dg3Bc;
use rssn::sim::dg3d::Dg3Limiter;
use rssn::sim::dg3d::Dg3Mesh;
use rssn::sim::dg3d::Dg3Options;
use rssn::sim::dg3d::Dg3Solution;
use rssn::sim::dg3d::Euler3d;
use rssn::sim::dg3d::LinearAdvection3d;
use rssn::sim::dg3d::solve_dg3d;
use std::f64::consts::PI;

const GAMMA: f64 = 1.4;
const SOD_SHOCK: f64 = 0.850_431;

fn periodic() -> [Dg3Bc; 6] {
    std::array::from_fn(|_| Dg3Bc::Periodic)
}

fn shock_limiters() -> [Dg3Limiter; 4] {
    [
        Dg3Limiter::Minmax,
        Dg3Limiter::Characteristic,
        Dg3Limiter::Hierarchical,
        Dg3Limiter::Weno,
    ]
}

/// Element average of a variable, from the LGL quadrature (weights sum to one).
fn element_mean(
    s: &Dg3Solution,
    e: usize,
    var: usize,
) -> f64 {
    let n1 = s.mesh.p + 1;
    let np = s.mesh.nodes_per_element();
    let mut acc = 0.0;
    for nd in 0..np {
        let ij = [nd % n1, (nd / n1) % n1, nd / (n1 * n1)];
        acc += s.w[ij[0]] * s.w[ij[1]] * s.w[ij[2]] / 8.0 * s.value(e, nd, var);
    }
    acc
}

fn element_at(
    ne: [usize; 3],
    c: [usize; 3],
) -> usize {
    c[0] + ne[0] * (c[1] + ne[1] * c[2])
}

/// Position where the element-average density crosses the midpoint between the
/// exact post-shock (`0.265_57`) and pre-shock (`0.125`) densities, scanning
/// from the right: the self-similar shock front of the Sod problem.
fn shock_position(
    s: &Dg3Solution,
    axis: usize,
) -> f64 {
    const POST: f64 = 0.265_574;
    const PRE: f64 = 0.125;
    let target = f64::midpoint(POST, PRE);
    let m = s.mesh;
    let h = (m.hi[axis] - m.lo[axis]) / m.ne[axis] as f64;
    let rho: Vec<f64> = (0..m.ne[axis])
        .map(|i| {
            let mut c = [0; 3];
            c[axis] = i;
            element_mean(s, element_at(m.ne, c), 0)
        })
        .collect();
    for i in (0..m.ne[axis] - 1).rev() {
        if rho[i] >= target && rho[i + 1] < target {
            let frac = (rho[i] - target) / (rho[i] - rho[i + 1]).max(1e-300);
            return m.lo[axis] + (i as f64 + 0.5 + frac) * h;
        }
    }
    // fallback: steepest element-average drop
    let mut best = 0.0;
    let mut pos = m.lo[axis];
    for i in 0..m.ne[axis] - 1 {
        let d = (rho[i + 1] - rho[i]).abs();
        if d > best {
            best = d;
            pos = m.lo[axis] + (i + 1) as f64 * h;
        }
    }
    pos
}

/// Exact solution of `u_t + a·∇u = 0` for the smooth wave
/// `sin(2πx) cos(2πy) cos(2πz)` and `a = (1, 0.5, -0.5)`.
fn advection_exact(
    x: [f64; 3],
    t: f64,
) -> f64 {
    (2.0 * PI * (x[0] - t)).sin()
        * (2.0 * PI * (x[1] - 0.5 * t)).cos()
        * (2.0 * PI * (x[2] + 0.5 * t)).cos()
}

/// Smooth periodic advection to `t = 0.5` on an `ne³` mesh of degree `p`.
fn solve_smooth_advection(
    limiter: Dg3Limiter,
    ne: usize,
    p: usize,
) -> Dg3Solution {
    let law = LinearAdvection3d { a: [1.0, 0.5, -0.5] };
    let mesh = Dg3Mesh {
        ne: [ne; 3],
        lo: [0.0; 3],
        hi: [1.0; 3],
        p,
    };
    let opt = Dg3Options {
        limiter,
        ..Dg3Options::default()
    };
    solve_dg3d(
        &mesh,
        &law,
        &periodic(),
        &|x| vec![advection_exact(x, 0.0)],
        0.5,
        0.3,
        &opt,
    )
    .unwrap()
}

#[test]
fn dg3d_smooth_advection_keeps_design_order_with_every_limiter() {
    // The TVB troubled-element test leaves this smooth wave untouched (the
    // default M h² threshold is well above its face slopes), so every limiter
    // must reproduce the design order p + 1 rather than degrade it — and, since
    // no element is flagged, the limited state must equal the unlimited one.
    for (p, n0, n1, min_order) in [
        (2_usize, 2_usize, 4_usize, 2.8_f64),
        (2, 4, 8, 2.8),
        (3, 2, 4, 3.5),
    ] {
        let base = |ne: usize| {
            solve_smooth_advection(Dg3Limiter::None, ne, p).l2_error(0, |x| advection_exact(x, 0.5))
        };
        let (b1, b2) = (base(n0), base(n1));
        println!(
            "advection p={p} unlimited: errors {b1:.3e} {b2:.3e}, order {:.2}",
            (b1 / b2).log2()
        );
        for limiter in shock_limiters() {
            let e1 =
                solve_smooth_advection(limiter, n0, p).l2_error(0, |x| advection_exact(x, 0.5));
            let e2 =
                solve_smooth_advection(limiter, n1, p).l2_error(0, |x| advection_exact(x, 0.5));
            let order = (e1 / e2).log2();
            println!("advection p={p} {limiter:?}: errors {e1:.3e} {e2:.3e}, order {order:.2}");
            assert!(
                order >= min_order,
                "{limiter:?} degrades the p={p} design order {}: {e1} -> {e2} (order {order})",
                p + 1
            );
            // the limiter is a bitwise no-op on a smooth solution
            assert_eq!(e1, b1, "{limiter:?} changed the p={p} ne={n0} solution");
            assert_eq!(e2, b2, "{limiter:?} changed the p={p} ne={n1} solution");
        }
    }
}

#[test]
fn dg3d_limiters_stay_inactive_on_a_smooth_wave() {
    // Locality: a smooth solution must never be flagged, so the limiter is a
    // no-op away from shocks (this is what keeps the order above intact).
    for limiter in shock_limiters() {
        for ne in [4_usize, 8] {
            let s = solve_smooth_advection(limiter, ne, 2);
            let count = s.troubled.iter().filter(|&&flag| flag).count();
            println!(
                "{limiter:?} ne={ne}: troubled_cells {}, final flags {count}/{}",
                s.troubled_cells,
                s.troubled.len()
            );
            assert_eq!(
                s.troubled.len(),
                s.mesh.elements(),
                "{limiter:?}: flags must cover the mesh"
            );
            assert_eq!(s.troubled_cells, 0, "{limiter:?} fired on a smooth wave");
            assert_eq!(count, 0, "{limiter:?} flagged a smooth element");
        }
    }
    // Without a limiter the diagnostics are empty rather than "all clear".
    let s = solve_smooth_advection(Dg3Limiter::None, 4, 2);
    assert!(s.troubled.is_empty());
    assert_eq!(s.troubled_cells, 0);
}

fn sod_state(
    x: [f64; 3],
    axis: usize,
) -> Vec<f64> {
    let (rho, p) = if x[axis] < 0.5 {
        (1.0, 1.0)
    } else {
        (0.125, 0.1)
    };
    vec![rho, 0.0, 0.0, 0.0, p / (GAMMA - 1.0)]
}

fn sod_bcs(axis: usize) -> [Dg3Bc; 6] {
    let mut bcs = periodic();
    bcs[2 * axis] = Dg3Bc::Outflow;
    bcs[2 * axis + 1] = Dg3Bc::Outflow;
    bcs
}

fn primitive(
    s: &Dg3Solution,
    e: usize,
    nd: usize,
) -> (f64, f64) {
    let q: Vec<f64> = (0..5).map(|v| s.value(e, nd, v)).collect();
    let rho = q[0];
    let p = (GAMMA - 1.0) * (q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / rho);
    (rho, p)
}

/// Total variation of the element-average density along `axis`; the Sod data is
/// constant in the transverse directions, so the extruded solution is 1D.
fn tv_along(
    s: &Dg3Solution,
    axis: usize,
) -> f64 {
    let m = s.mesh;
    let mut c = [0_usize; 3];
    let mut prev = f64::NAN;
    let mut tv = 0.0;
    for i in 0..m.ne[axis] {
        c[axis] = i;
        let rho = element_mean(s, element_at(m.ne, c), 0);
        if i > 0 {
            tv += (rho - prev).abs();
        }
        prev = rho;
    }
    tv
}

fn element_axis_index(
    ne: [usize; 3],
    e: usize,
    axis: usize,
) -> usize {
    [e % ne[0], (e / ne[0]) % ne[1], e / (ne[0] * ne[1])][axis]
}

fn check_sod(axis: usize) {
    let law = Euler3d { gamma: GAMMA };
    let mut ne = [2; 3];
    ne[axis] = 12;
    let mesh = Dg3Mesh {
        ne,
        lo: [0.0; 3],
        hi: [1.0; 3],
        p: 2,
    };
    let bcs = sod_bcs(axis);
    let dx = 1.0 / 12.0;
    for limiter in shock_limiters() {
        let opt = Dg3Options {
            limiter,
            positivity: true,
            limiter_tvd: 10.0,
            ..Dg3Options::default()
        };
        let s = solve_dg3d(&mesh, &law, &bcs, &|x| sod_state(x, axis), 0.2, 0.35, &opt).unwrap();
        let (mut rho_min, mut rho_max, mut p_min, mut p_max) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for e in 0..mesh.elements() {
            for nd in 0..mesh.nodes_per_element() {
                let (rho, p) = primitive(&s, e, nd);
                assert!(
                    rho.is_finite() && p.is_finite(),
                    "{limiter:?}: non-finite state"
                );
                rho_min = rho_min.min(rho);
                rho_max = rho_max.max(rho);
                p_min = p_min.min(p);
                p_max = p_max.max(p);
            }
        }
        let xs = shock_position(&s, axis);
        let drift = (xs - SOD_SHOCK).abs() / dx;
        let flagged: Vec<usize> = (0..mesh.elements())
            .filter(|&e| s.troubled[e])
            .map(|e| element_axis_index(ne, e, axis))
            .collect();
        let band = format!("{flagged:?}");
        println!(
            "sod axis {axis} {limiter:?}: rho [{rho_min:.4}, {rho_max:.4}], p [{p_min:.4}, {p_max:.4}], shock {xs:.4} ({drift:.2} cells), {} flags in the run, final flags {band}",
            s.troubled_cells
        );
        // non-oscillatory: at most 5% of the density / pressure jump beyond the data range
        assert!(
            rho_min > 0.125 - 0.05 * 0.875 && rho_max < 1.0 + 0.05 * 0.875,
            "{limiter:?}: density overshoot [{rho_min}, {rho_max}]"
        );
        assert!(
            p_min > 0.1 - 0.05 * 0.9 && p_max < 1.0 + 0.05 * 0.9,
            "{limiter:?}: pressure overshoot [{p_min}, {p_max}]"
        );
        assert!(
            rho_min > 0.0 && p_min > 0.0,
            "{limiter:?}: non-positive density or pressure"
        );
        // the shock must sit within three cells of the exact self-similar position
        assert!(
            drift < 3.0,
            "{limiter:?} axis {axis}: shock at {xs} ({drift} cells off)"
        );
    }
}

#[test]
fn dg3d_sod_along_x_with_all_limiters() {
    check_sod(0);
}

#[test]
fn dg3d_sod_along_y_with_all_limiters() {
    check_sod(1);
}

#[test]
fn dg3d_sod_along_z_with_all_limiters() {
    check_sod(2);
}

#[test]
fn dg3d_sod_total_variation_does_not_grow_and_does_not_undershoot() {
    // The exact Sod solution is monotone, so its total variation is the initial
    // jump 1 - 0.125 = 0.875 for all time. A limited DG solution must not add
    // variation after the shock has formed, and the min/max limiter must not
    // dip below the initial minimum of density or pressure (the small residual
    // is the O(h) smearing of the plateau, not an oscillation: it shrinks under
    // refinement).
    let law = Euler3d { gamma: GAMMA };
    let bcs = sod_bcs(0);
    let init = |x: [f64; 3]| sod_state(x, 0);
    let tv_initial = 1.0 - 0.125;
    let run = |nx: usize, limiter: Dg3Limiter, t_end: f64| {
        let mesh = Dg3Mesh {
            ne: [nx, 2, 2],
            lo: [0.0; 3],
            hi: [1.0; 3],
            p: 2,
        };
        let opt = Dg3Options {
            limiter,
            positivity: true,
            limiter_tvd: 10.0,
            ..Dg3Options::default()
        };
        solve_dg3d(&mesh, &law, &bcs, &init, t_end, 0.35, &opt).unwrap()
    };
    let extrema = |s: &Dg3Solution| {
        let (mut rho_min, mut p_min) = (f64::INFINITY, f64::INFINITY);
        for e in 0..s.mesh.elements() {
            for nd in 0..s.mesh.nodes_per_element() {
                let (rho, p) = primitive(s, e, nd);
                rho_min = rho_min.min(rho);
                p_min = p_min.min(p);
            }
        }
        (rho_min, p_min)
    };
    for limiter in shock_limiters() {
        let formed = run(12, limiter, 0.1);
        let late = run(12, limiter, 0.2);
        let (tv_formed, tv_late) = (tv_along(&formed, 0), tv_along(&late, 0));
        let (rho_min, p_min) = extrema(&late);
        println!(
            "sod TV {limiter:?}: TV(0.1) {tv_formed:.6} -> TV(0.2) {tv_late:.6} (initial {tv_initial}), min rho {rho_min:.9}, min p {p_min:.9}"
        );
        if limiter == Dg3Limiter::Minmax {
            assert!(
                tv_late <= tv_formed,
                "{limiter:?}: TV grew from {tv_formed} to {tv_late} after the shock formed"
            );
            assert!(
                tv_late <= tv_initial,
                "{limiter:?}: TV {tv_late} exceeds the initial variation {tv_initial}"
            );
            // the residual undershoot is O(h) smearing, not an oscillation
            let (rho_fine, p_fine) = extrema(&run(24, limiter, 0.2));
            println!("sod TV Minmax nx=24: min rho {rho_fine:.9}, min p {p_fine:.9}");
            assert!(
                rho_min >= 0.125 - 2e-3 && p_min >= 0.1 - 2e-3,
                "{limiter:?}: undershoot below the initial minima: rho {rho_min}, p {p_min}"
            );
            assert!(
                rho_fine >= rho_min && p_fine >= p_min,
                "{limiter:?}: the undershoot must shrink under refinement: rho {rho_min} -> {rho_fine}, p {p_min} -> {p_fine}"
            );
        } else {
            // the moment limiters do not guarantee a monotone profile, but the
            // variation must stay bounded by the data range
            assert!(
                tv_late <= tv_initial * 1.03,
                "{limiter:?}: TV {tv_late} exceeds the initial variation {tv_initial}"
            );
        }
    }
}

#[test]
fn dg3d_characteristic_limiter_conserves_euler_quantities() {
    // A periodic compression wave steepens into a shock; mass, momentum and
    // energy are conserved because every limiter preserves element averages.
    let init = |x: [f64; 3]| {
        let rho =
            1.0 + 0.3 * (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos() * (2.0 * PI * x[2]).cos();
        let u = 1.5 * (2.0 * PI * x[0]).sin();
        vec![
            rho,
            rho * u,
            0.0,
            0.0,
            1.0 / (GAMMA - 1.0) + 0.5 * rho * u * u,
        ]
    };
    let law = Euler3d { gamma: GAMMA };
    let mesh = Dg3Mesh {
        ne: [8, 2, 2],
        lo: [0.0; 3],
        hi: [1.0; 3],
        p: 2,
    };
    let opt = Dg3Options {
        limiter: Dg3Limiter::Characteristic,
        positivity: true,
        ..Dg3Options::default()
    };
    let s0 = solve_dg3d(&mesh, &law, &periodic(), &init, 0.0, 0.25, &opt).unwrap();
    let s1 = solve_dg3d(&mesh, &law, &periodic(), &init, 0.2, 0.25, &opt).unwrap();
    for (name, f) in [
        ("mass", (|q: &[f64]| q[0]) as fn(&[f64]) -> f64),
        ("x-momentum", |q: &[f64]| q[1]),
        ("y-momentum", |q: &[f64]| q[2]),
        ("z-momentum", |q: &[f64]| q[3]),
        ("energy", |q: &[f64]| q[4]),
    ] {
        let (a, b) = (s0.integral(f), s1.integral(f));
        let scale = a.abs().max(b.abs()).max(1.0);
        let rel = (b - a).abs() / scale;
        println!("conservation {name}: {a:.15e} -> {b:.15e} (rel {rel:.2e})");
        assert!(rel < 1e-10, "{name} not conserved: {a} -> {b} (rel {rel})");
    }
}

#[test]
fn dg3d_blast_stability_regression_with_limiters() {
    // High-pressure sphere (pressure ratio 200) in a low-pressure medium,
    // outflow on all faces: a 3D blast wave. This is a stability regression:
    // finite, positive, bounded states with each limiter at a modest
    // resolution where the unlimited DGSEM is not stable.
    let init = |x: [f64; 3]| {
        let r2 = (x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2) + (x[2] - 0.5).powi(2);
        let (rho, p) = if r2 < 0.2 * 0.2 {
            (1.0, 20.0)
        } else {
            (1.0, 0.1)
        };
        vec![rho, 0.0, 0.0, 0.0, p / (GAMMA - 1.0)]
    };
    let law = Euler3d { gamma: GAMMA };
    let mesh = Dg3Mesh {
        ne: [5, 5, 5],
        lo: [0.0; 3],
        hi: [1.0; 3],
        p: 1,
    };
    let bcs: [Dg3Bc; 6] = std::array::from_fn(|_| Dg3Bc::Outflow);
    for limiter in shock_limiters() {
        let opt = Dg3Options {
            limiter,
            positivity: true,
            limiter_tvd: 2.0,
            ..Dg3Options::default()
        };
        let s = solve_dg3d(&mesh, &law, &bcs, &init, 0.03, 0.15, &opt).unwrap();
        let (mut rho_min, mut rho_max, mut p_min) = (f64::INFINITY, 0.0_f64, f64::INFINITY);
        for e in 0..mesh.elements() {
            for nd in 0..mesh.nodes_per_element() {
                let (rho, p) = primitive(&s, e, nd);
                assert!(
                    rho.is_finite() && p.is_finite(),
                    "{limiter:?}: non-finite state"
                );
                assert!(
                    rho > 0.0 && p > 0.0,
                    "{limiter:?}: negative density or pressure"
                );
                rho_min = rho_min.min(rho);
                rho_max = rho_max.max(rho);
                p_min = p_min.min(p);
            }
        }
        println!(
            "blast {limiter:?}: rho [{rho_min:.4}, {rho_max:.4}], p_min {p_min:.3e}, steps {}",
            s.steps
        );
        assert!(
            rho_min > 1e-4 && p_min > 1e-4,
            "{limiter:?}: collapsed to rho {rho_min}, p {p_min}"
        );
        assert!(rho_max < 20.0, "{limiter:?}: density blew up to {rho_max}");
    }
}

#[test]
fn dg3d_euler_characteristic_transform_diagonalises_the_jacobian() {
    let law = Euler3d { gamma: GAMMA };
    let (rho, u, v, w, p) = (1.2, 0.5, -0.3, 0.2, 0.9);
    let q = [
        rho,
        rho * u,
        rho * v,
        rho * w,
        p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v + w * w),
    ];
    let c = (GAMMA * p / rho).sqrt();
    for dir in 0..3 {
        let un = [u, v, w][dir];
        let expected = [un - c, un, un, un, un + c];
        for k in 0..5 {
            let mut wv = [0.0; 5];
            wv[k] = 1.0;
            let mut r = [0.0; 5];
            law.from_characteristic(&q, dir, &wv, &mut r);
            let eps = 1e-6;
            let (mut qp, mut qm) = (q, q);
            for i in 0..5 {
                qp[i] += eps * r[i];
                qm[i] -= eps * r[i];
            }
            let (mut fp, mut fm) = ([0.0; 5], [0.0; 5]);
            law.flux(&qp, dir, &mut fp);
            law.flux(&qm, dir, &mut fm);
            for i in 0..5 {
                let jr = (fp[i] - fm[i]) / (2.0 * eps);
                assert!(
                    (jr - expected[k] * r[i]).abs() < 1e-6,
                    "dir {dir} wave {k} comp {i}: {jr} vs {}",
                    expected[k] * r[i]
                );
            }
            let mut back = [0.0; 5];
            law.to_characteristic(&q, dir, &r, &mut back);
            for i in 0..5 {
                assert!((back[i] - wv[i]).abs() < 1e-12, "L R = I, dir {dir}");
            }
        }
    }
}

#[test]
fn moment_limiters_use_mirror_ghosts_at_a_wall() {
    // A gas streaming into a wall at x = 1 (reflected shock) equals the left half of the
    // doubled box with the antisymmetric velocity: the wall's ghost element is the mirror
    // image the doubled box provides.
    let law = Euler3d { gamma: GAMMA };
    // A smooth antisymmetric velocity profile: u(2 - x) = -u(x), u(1) = 0 exactly.
    let state = |x: f64| {
        let (rho, p) = (1.0, 1.0);
        let u = -((x - 1.0) / 0.04).tanh();
        vec![
            rho,
            rho * u,
            0.0,
            0.0,
            p / (GAMMA - 1.0) + 0.5 * rho * u * u,
        ]
    };
    let run = |nx: usize, hi: f64, right: Dg3Bc, limiter: Dg3Limiter| {
        let mesh = Dg3Mesh {
            ne: [nx, 2, 2],
            lo: [0.0; 3],
            hi: [hi, 0.1, 0.1],
            p: 1,
        };
        let mut bcs: [Dg3Bc; 6] = std::array::from_fn(|_| Dg3Bc::Periodic);
        bcs[0] = Dg3Bc::Outflow;
        bcs[1] = right;
        let opt = Dg3Options {
            limiter,
            positivity: true,
            ..Dg3Options::default()
        };
        solve_dg3d(&mesh, &law, &bcs, &|x| state(x[0]), 0.2, 0.2, &opt).unwrap()
    };
    // Both boxes have 16 elements per unit length.
    for limiter in [Dg3Limiter::Hierarchical, Dg3Limiter::Characteristic] {
        let wall = run(16, 1.0, Dg3Bc::Wall, limiter);
        let doubled = run(32, 2.0, Dg3Bc::Outflow, limiter);
        let np = wall.mesh.nodes_per_element();
        let ne = wall.mesh.elements();
        let (mut worst, mut rho_max) = (0.0_f64, 0.0_f64);
        for e in 0..ne {
            // element e of the wall box = element with the same (ix, iy, iz) in the doubled box
            let (ix, rest) = (e % 16, e / 16);
            let e2 = ix + 32 * rest;
            for nd in 0..np {
                for c in 0..5 {
                    worst = worst.max((wall.value(e, nd, c) - doubled.value(e2, nd, c)).abs());
                }
                rho_max = rho_max.max(wall.value(e, nd, 0));
            }
        }
        println!("{limiter:?}: wall vs doubled box {worst:e}, max density {rho_max}");
        assert!(
            rho_max > 1.3,
            "the reflected shock must have formed: {rho_max}"
        );
        assert!(worst < 1e-8, "{limiter:?}: difference {worst}");
    }
}
