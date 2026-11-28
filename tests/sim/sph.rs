//! SPH: kernel normalisation, density summation, conservation and the Sod shock tube.

use std::f64::consts::PI;

use rssn::sim::sph::*;

const GAMMA: f64 = 1.4;

#[test]
fn kernels_are_normalised_and_differentiable_in_every_dimension() {
    for kind in [KernelKind::CubicSpline, KernelKind::WendlandC2] {
        for d in 1..=3usize {
            let h = 0.7;
            let n = 20_000;
            let rmax = 2.0 * h;
            let dr = rmax / n as f64;
            let omega = match d {
                | 1 => 2.0,
                | 2 => 2.0 * PI,
                | _ => 4.0 * PI,
            };
            let mut integral = 0.0;
            for k in 0..n {
                let r = (k as f64 + 0.5) * dr;
                integral += omega * r.powi(d as i32 - 1) * kind.w(d, r, h) * dr;
            }
            assert!(
                (integral - 1.0).abs() < 1e-5,
                "{kind:?} D={d}: integral {integral}"
            );
            for r in [0.2, 0.7, 1.1, 1.3] {
                let e = 1e-6;
                let num = (kind.w(d, r + e, h) - kind.w(d, r - e, h)) / (2.0 * e);
                assert!(
                    (kind.dw(d, r, h) - num).abs() < 1e-5 * num.abs().max(1.0),
                    "{kind:?} D={d} r={r}"
                );
            }
            assert_eq!(kind.w(d, 2.0 * h, h), 0.0);
        }
    }
}

fn lattice_2d(
    n: usize,
    l: f64,
) -> Vec<[f64; 2]> {
    let dx = l / n as f64;
    let mut v = Vec::new();
    for i in 0..n {
        for j in 0..n {
            v.push([(i as f64 + 0.5) * dx, (j as f64 + 0.5) * dx]);
        }
    }
    v
}

#[test]
fn uniform_lattices_have_uniform_density_in_two_and_three_dimensions() {
    for kind in [KernelKind::CubicSpline, KernelKind::WendlandC2] {
        let n = 16;
        let pos = lattice_2d(n, 1.0);
        let m = 1.0 / (n * n) as f64;
        let mut s = Sph::<2>::new(
            pos,
            vec![m; n * n],
            vec![1.0; n * n],
            [0.0; 2],
            [1.0; 2],
            [SphBoundary::Periodic; 2],
        );
        s.kernel = kind;
        s.eta = 1.4;
        s.update_density();
        for r in &s.rho {
            assert!((r - 1.0).abs() < 0.03, "{kind:?} 2-D density {r}");
        }
        let n3 = 8;
        let dx = 1.0 / n3 as f64;
        let mut pos3 = Vec::new();
        for i in 0..n3 {
            for j in 0..n3 {
                for k in 0..n3 {
                    pos3.push([
                        (i as f64 + 0.5) * dx,
                        (j as f64 + 0.5) * dx,
                        (k as f64 + 0.5) * dx,
                    ]);
                }
            }
        }
        let cnt = pos3.len();
        let mut s3 = Sph::<3>::new(
            pos3,
            vec![1.0 / cnt as f64; cnt],
            vec![1.0; cnt],
            [0.0; 3],
            [1.0; 3],
            [SphBoundary::Periodic; 3],
        );
        s3.kernel = kind;
        s3.eta = 1.4;
        s3.update_density();
        for r in &s3.rho {
            assert!((r - 1.0).abs() < 0.04, "{kind:?} 3-D density {r}");
        }
    }
}

#[test]
fn periodic_gas_conserves_momentum_and_energy_and_a_uniform_state_stays_at_rest() {
    let n = 12;
    let m = 1.0 / (n * n) as f64;
    // uniform state: no motion
    let mut s = Sph::<2>::new(
        lattice_2d(n, 1.0),
        vec![m; n * n],
        vec![2.5; n * n],
        [0.0; 2],
        [1.0; 2],
        [SphBoundary::Periodic; 2],
    );
    s.eta = 1.4;
    s.update_density();
    let dt = s.suggest_dt(0.2);
    for _ in 0..10 {
        s.step(dt);
    }
    for v in &s.vel {
        assert!(v[0].abs() < 1e-3 && v[1].abs() < 1e-3, "{v:?}");
    }
    // a velocity pulse conserves momentum and energy
    let mut pos = lattice_2d(n, 1.0);
    for (k, p) in pos.iter_mut().enumerate() {
        p[0] += 0.01 * ((k * 7 % 5) as f64 - 2.0) / 5.0;
        p[1] += 0.01 * ((k * 3 % 7) as f64 - 3.0) / 7.0;
    }
    let mut s = Sph::<2>::new(
        pos,
        vec![m; n * n],
        vec![2.5; n * n],
        [0.0; 2],
        [1.0; 2],
        [SphBoundary::Periodic; 2],
    );
    s.eta = 1.4;
    for i in 0..n * n {
        s.vel[i][0] = 0.3 * (2.0 * PI * s.pos[i][0]).sin();
    }
    s.update_density();
    let (e0, p0) = (s.total_energy(), s.momentum());
    let dt = s.suggest_dt(0.2);
    for _ in 0..40 {
        s.step(dt);
    }
    let p1 = s.momentum();
    assert!(
        (p1[0] - p0[0]).abs() < 1e-6 && (p1[1] - p0[1]).abs() < 1e-6,
        "{p0:?} -> {p1:?}"
    );
    assert!(
        (s.total_energy() - e0).abs() < 2e-3 * e0,
        "energy {} -> {}",
        e0,
        s.total_energy()
    );
}

/// Exact Sod solution at time t for the diaphragm at x = 0 (rho, u, p), gamma = 1.4.
fn sod_exact(
    x: f64,
    t: f64,
) -> (f64, f64, f64) {
    let (rl, pl, rr, pr) = (1.0, 1.0, 0.125, 0.1);
    let cl = (GAMMA * pl / rl).sqrt();
    // Newton iteration for p*
    let f = |p: f64, rk: f64, pk: f64| {
        let ck = (GAMMA * pk / rk).sqrt();
        if p > pk {
            let a = 2.0 / ((GAMMA + 1.0) * rk);
            let b = (GAMMA - 1.0) / (GAMMA + 1.0) * pk;
            (p - pk) * (a / (p + b)).sqrt()
        } else {
            2.0 * ck / (GAMMA - 1.0) * ((p / pk).powf((GAMMA - 1.0) / (2.0 * GAMMA)) - 1.0)
        }
    };
    let (mut lo, mut hi) = (1e-6, 10.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid, rl, pl) + f(mid, rr, pr) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let ps = 0.5 * (lo + hi);
    let us = 0.5 * (f(ps, rr, pr) - f(ps, rl, pl));
    let rsl = rl * (ps / pl).powf(1.0 / GAMMA);
    let rsr = rr
        * ((ps / pr + (GAMMA - 1.0) / (GAMMA + 1.0))
            / ((GAMMA - 1.0) / (GAMMA + 1.0) * ps / pr + 1.0));
    let csl = (GAMMA * ps / rsl).sqrt();
    let shock = us * rsr / (rsr - rr);
    let s = x / t;
    if s < -cl {
        (rl, 0.0, pl)
    } else if s < us - csl {
        let u = 2.0 / (GAMMA + 1.0) * (cl + s);
        let c = cl - 0.5 * (GAMMA - 1.0) * u;
        let rho = rl * (c / cl).powf(2.0 / (GAMMA - 1.0));
        (rho, u, pl * (rho / rl).powf(GAMMA))
    } else if s < us {
        (rsl, us, ps)
    } else if s < shock {
        (rsr, us, ps)
    } else {
        (rr, 0.0, pr)
    }
}

#[test]
fn exact_sod_reference_has_the_textbook_star_state() {
    let (rho, u, p) = sod_exact(0.1, 0.2);
    assert!(
        (rho - 0.426_32).abs() < 1e-4 && (u - 0.927_45).abs() < 1e-4 && (p - 0.303_13).abs() < 1e-4,
        "{rho} {u} {p}"
    );
    assert!((sod_exact(0.3, 0.2).0 - 0.265_57).abs() < 1e-4);
}

#[test]
fn sph_reproduces_the_sod_shock_tube() {
    // equal-mass particles: dx = 0.001 on the left (rho = 1), 0.008 on the right (rho = 0.125)
    let m = 0.001;
    let mut pos = Vec::new();
    let mut u = Vec::new();
    let mut x = -0.5 + 0.0005;
    while x < 0.0 {
        pos.push([x]);
        u.push(1.0 / ((GAMMA - 1.0) * 1.0));
        x += 0.001;
    }
    let mut x = 0.004;
    while x < 0.5 {
        pos.push([x]);
        u.push(0.1 / ((GAMMA - 1.0) * 0.125));
        x += 0.008;
    }
    let n = pos.len();
    let mut s = Sph::<1>::new(pos, vec![m; n], u, [-0.5], [0.5], [SphBoundary::Reflecting]);
    s.gamma = GAMMA;
    s.eta = 1.2;
    s.alpha = 1.0;
    s.beta = 2.0;
    s.update_density();
    let t_end = 0.2;
    let mut t = 0.0;
    let p0 = s.momentum()[0];
    while t < t_end {
        let dt = s.suggest_dt(0.25).min(t_end - t);
        s.step(dt);
        t += dt;
    }
    let mut err = 0.0;
    let mut cnt = 0.0;
    for i in 0..n {
        let x = s.pos[i][0];
        if x.abs() < 0.4 {
            err += (s.rho[i] - sod_exact(x, t_end).0).abs();
            cnt += 1.0;
        }
    }
    let l1 = err / cnt;
    assert!(l1 < 0.04, "mean density error {l1}");
    // plateau values near the contact discontinuity
    let near = |xc: f64| {
        let mut best = 0;
        for i in 0..n {
            if (s.pos[i][0] - xc).abs() < (s.pos[best][0] - xc).abs() {
                best = i;
            }
        }
        best
    };
    let i = near(0.1);
    assert!(
        (s.vel[i][0] - 0.927_45).abs() < 0.04,
        "star velocity {}",
        s.vel[i][0]
    );
    assert!(
        (s.pressure(i) - 0.303_13).abs() < 0.03,
        "star pressure {}",
        s.pressure(i)
    );
    let i = near(0.27);
    assert!(
        (s.rho[i] - 0.265_57).abs() < 0.04,
        "post-shock density {}",
        s.rho[i]
    );
    // the closed left wall pushes on the gas with the pressure jump until the waves arrive
    assert!(
        (s.momentum()[0] - p0 - 0.18).abs() < 0.01,
        "the left wall has delivered the impulse (P_L - P_R) t = 0.18: momentum {} -> {}",
        p0,
        s.momentum()[0]
    );
}

#[test]
fn fixed_boundary_particles_do_not_move_and_still_push_back() {
    // 1-D gas compressed against a frozen wall of boundary particles
    let m = 0.01;
    let mut pos = Vec::new();
    let mut x = 0.005;
    while x < 1.0 {
        pos.push([x]);
        x += 0.01;
    }
    let nf = pos.len();
    // frozen wall particles to the left of the origin
    for k in 1..=6 {
        pos.push([0.005 - 0.01 * k as f64]);
    }
    let n = pos.len();
    let mut s = Sph::<1>::new(
        pos,
        vec![m; n],
        vec![2.5; n],
        [-1.0],
        [2.0],
        [SphBoundary::Open],
    );
    for i in nf..n {
        s.fixed[i] = true;
    }
    for i in 0..nf {
        s.vel[i][0] = -0.5;
    }
    s.update_density();
    let wall: Vec<f64> = s.pos[nf..].iter().map(|p| p[0]).collect();
    for _ in 0..100 {
        let dt = s.suggest_dt(0.2);
        s.step(dt);
    }
    for (k, p) in s.pos[nf..].iter().enumerate() {
        assert_eq!(p[0], wall[k]);
    }
    // the fluid was decelerated by the wall: its mean velocity is less negative
    let mean: f64 = s.vel[..nf].iter().map(|v| v[0]).sum::<f64>() / nf as f64;
    assert!(mean > -0.5 + 1e-4, "mean velocity {mean}");
}

#[test]
fn two_dimensional_blast_stays_symmetric_and_runs_in_three_dimensions() {
    // 2-D: symmetric heated square on a periodic lattice remains mirror symmetric
    let n = 14;
    let m = 1.0 / (n * n) as f64;
    let pos = lattice_2d(n, 1.0);
    let u: Vec<f64> = pos
        .iter()
        .map(|p| {
            if (p[0] - 0.5).abs() < 0.1 && (p[1] - 0.5).abs() < 0.1 {
                10.0
            } else {
                1.0
            }
        })
        .collect();
    let mut s = Sph::<2>::new(
        pos,
        vec![m; n * n],
        u,
        [0.0; 2],
        [1.0; 2],
        [SphBoundary::Periodic; 2],
    );
    s.eta = 1.4;
    s.kernel = KernelKind::WendlandC2;
    s.update_density();
    for _ in 0..15 {
        let dt = s.suggest_dt(0.2);
        s.step(dt);
    }
    // mirror symmetry in x about 0.5: sum of x-velocity ~ 0
    let sum: f64 = s.vel.iter().map(|v| v[0]).sum();
    assert!(sum.abs() < 1e-8, "net x momentum {sum}");
    let outward = s
        .vel
        .iter()
        .zip(&s.pos)
        .map(|(v, p)| v[0] * (p[0] - 0.5) + v[1] * (p[1] - 0.5))
        .sum::<f64>();
    assert!(outward > 0.0, "the blast expands");
    // 3-D smoke run
    let k = 5;
    let dx = 1.0 / k as f64;
    let mut pos3 = Vec::new();
    for i in 0..k {
        for j in 0..k {
            for l in 0..k {
                pos3.push([
                    (i as f64 + 0.5) * dx,
                    (j as f64 + 0.5) * dx,
                    (l as f64 + 0.5) * dx,
                ]);
            }
        }
    }
    let cnt = pos3.len();
    let mut s3 = Sph::<3>::new(
        pos3,
        vec![1.0 / cnt as f64; cnt],
        vec![2.5; cnt],
        [0.0; 3],
        [1.0; 3],
        [SphBoundary::Periodic; 3],
    );
    s3.eta = 1.4;
    s3.vel[0] = [0.1, 0.0, 0.0];
    s3.update_density();
    let e0 = s3.total_energy();
    for _ in 0..5 {
        let dt = s3.suggest_dt(0.2);
        s3.step(dt);
    }
    assert!((s3.total_energy() - e0).abs() < 1e-2 * e0);
}

#[test]
fn self_gravity_of_a_uniform_ball_is_linear_in_radius_and_conserves_momentum() {
    // Cold uniform ball (cubic lattice cut at radius R): a(r) = -(4 pi / 3) G rho r inside.
    let (n, big_r) = (15usize, 1.0_f64);
    let d = 2.0 * big_r / n as f64;
    let mut pos = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let x = [
                    -big_r + (i as f64 + 0.5) * d,
                    -big_r + (j as f64 + 0.5) * d,
                    -big_r + (k as f64 + 0.5) * d,
                ];
                if x.iter().map(|v| v * v).sum::<f64>() < big_r * big_r {
                    pos.push(x);
                }
            }
        }
    }
    let count = pos.len();
    let total_mass = 1.0;
    let m = total_mass / count as f64;
    let mut s = Sph::<3>::new(
        pos.clone(),
        vec![m; count],
        vec![1e-10; count],
        [-1.5; 3],
        [1.5; 3],
        [SphBoundary::Open; 3],
    );
    s.alpha = 0.0;
    s.beta = 0.0;
    let g = 1.0;
    let rho = total_mass / (4.0 / 3.0 * PI * big_r.powi(3));
    let (a0, _) = s.derivatives();
    s.gravity = Some(SphGravity {
        g,
        eps: 0.05,
        theta: 0.3,
    });
    let (a1, _) = s.derivatives();
    // gravity only: difference of the two evaluations
    let (mut num, mut den, mut used) = (0.0, 0.0, 0);
    let (mut px, mut scale) = ([0.0; 3], 0.0);
    for i in 0..count {
        let ag: Vec<f64> = (0..3).map(|a| a1[i][a] - a0[i][a]).collect();
        for a in 0..3 {
            px[a] += m * ag[a];
        }
        scale += m * ag.iter().map(|v| v * v).sum::<f64>().sqrt();
        let r = pos[i].iter().map(|v| v * v).sum::<f64>().sqrt();
        if (0.3..0.7).contains(&r) {
            // radial component divided by r
            let ar: f64 = (0..3).map(|a| ag[a] * pos[i][a]).sum::<f64>() / r;
            num += ar / r;
            den += 1.0;
            used += 1;
        }
    }
    let slope = num / den;
    let expected = -4.0 * PI / 3.0 * g * rho;
    println!("slope {slope} expected {expected} ({used} particles)");
    assert!(used > 100);
    assert!(
        (slope / expected - 1.0).abs() < 0.08,
        "a/r = {slope}, expected {expected}"
    );
    let pnorm = px.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(pnorm < 0.02 * scale, "net force {pnorm} vs {scale}");
    // the pair energy is negative and of the order -3/5 G M^2 / R
    let e = s.gravitational_energy().unwrap();
    assert!(
        e < 0.0 && (e / (-0.6 * g / big_r) - 1.0).abs() < 0.15,
        "energy {e}"
    );
    // without gravity nothing changes
    s.gravity = None;
    assert!(s.gravitational_energy().is_none());
}
