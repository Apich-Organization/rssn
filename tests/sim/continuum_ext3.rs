//! SPH grad-h equations and dissipative particle dynamics.

use rssn::sim::dpd::{Dpd, DpdParams};
use rssn::sim::sph::{Sph, SphBoundary};

const GAMMA: f64 = 1.4;

fn sod_system(grad_h: bool) -> Sph<1> {
    let m = 0.001;
    let mut pos = Vec::new();
    let mut u = Vec::new();
    let mut x = -0.5 + 0.0005;
    while x < 0.0 {
        pos.push([x]);
        u.push(1.0 / (GAMMA - 1.0));
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
    s.grad_h = grad_h;
    s.update_density();
    s
}

fn evolve(s: &mut Sph<1>, t_end: f64) {
    let mut t = 0.0;
    while t < t_end {
        let dt = s.suggest_dt(0.25).min(t_end - t);
        s.step(dt);
        t += dt;
    }
}

#[test]
fn grad_h_sph_conserves_and_reproduces_the_sod_tube() {
    let mut s = sod_system(true);
    let e0 = s.total_energy();
    evolve(&mut s, 0.2);
    assert!((s.total_energy() - e0).abs() < 1e-3 * e0.abs(), "energy {} -> {}", e0, s.total_energy());
    // left plateau of the star region: u* = 0.92745, p* = 0.30313
    let mut best = 0;
    for i in 0..s.pos.len() {
        if (s.pos[i][0] - 0.1).abs() < (s.pos[best][0] - 0.1).abs() {
            best = i;
        }
    }
    assert!((s.vel[best][0] - 0.927_45).abs() < 0.05, "star velocity {}", s.vel[best][0]);
    assert!((s.pressure(best) - 0.303_13).abs() < 0.04, "star pressure {}", s.pressure(best));
}

#[test]
fn grad_h_equals_standard_sph_when_smoothing_lengths_are_uniform() {
    let n = 100;
    let pos: Vec<[f64; 1]> = (0..n).map(|i| [(i as f64 + 0.5) / n as f64]).collect();
    let make = |g: bool| {
        let u: Vec<f64> = (0..n).map(|i| 1.0 + 0.2 * (2.0 * std::f64::consts::PI * pos[i][0]).sin()).collect();
        let mut s = Sph::<1>::new(pos.clone(), vec![1.0 / n as f64; n], u, [0.0], [1.0], [SphBoundary::Periodic]);
        s.alpha = 0.0;
        s.beta = 0.0;
        s.grad_h = g;
        s
    };
    let (mut a, mut b) = (make(false), make(true));
    let (acc_a, _) = a.derivatives();
    let (acc_b, _) = b.derivatives();
    let scale = acc_a.iter().map(|x| x[0].abs()).fold(0.0, f64::max);
    let diff = acc_a.iter().zip(&acc_b).map(|(x, y)| (x[0] - y[0]).abs()).fold(0.0, f64::max);
    assert!(scale > 0.0 && diff < 0.05 * scale, "uniform-h difference {diff} vs scale {scale}");
    // pairwise antisymmetric forces: momentum is conserved
    let total: f64 = acc_b.iter().map(|x| x[0]).sum::<f64>() / n as f64;
    assert!(total.abs() < 1e-12, "net force {total}");
}

fn gw_params(dim: usize, l: f64) -> DpdParams {
    DpdParams { dim, box_len: [l, l, l], a: 25.0, gamma: 4.5, kt: 1.0, rc: 1.0, mass: 1.0 }
}

#[test]
fn dpd_fluid_has_the_target_temperature_conserves_momentum_and_follows_the_gw_eos() {
    let l = 6.0;
    let n = (3.0 * l * l * l) as usize;
    let mut d = Dpd::new(n, gw_params(3, l), 11).unwrap();
    assert!((d.density() - 3.0).abs() < 1e-12);
    d.run(0.04, 300);
    let (mut t_sum, mut p_sum) = (0.0, 0.0);
    let samples = 300;
    for _ in 0..samples {
        d.step(0.04);
        t_sum += d.temperature();
        p_sum += d.pressure();
    }
    let (t_mean, p_mean) = (t_sum / samples as f64, p_sum / samples as f64);
    assert!((t_mean - 1.0).abs() < 0.03, "temperature {t_mean}");
    // Groot-Warren: p = rho kT + 0.101 a rho^2 = 25.7
    assert!((p_mean - 25.7).abs() < 0.04 * 25.7, "pressure {p_mean}");
    let m = d.momentum();
    assert!(m.iter().all(|x| x.abs() < 1e-8), "momentum {m:?}");
}

#[test]
fn dpd_2d_and_input_checks() {
    let l = 8.0;
    let n = (3.0 * l * l) as usize;
    let mut d = Dpd::new(n, gw_params(2, l), 5).unwrap();
    d.run(0.04, 400);
    let mut t_sum = 0.0;
    for _ in 0..200 {
        d.step(0.04);
        t_sum += d.temperature();
    }
    assert!((t_sum / 200.0 - 1.0).abs() < 0.04);
    assert!(d.pos.iter().all(|p| p[2] == 0.0 && p[0] >= 0.0 && p[0] < l && p[1] >= 0.0 && p[1] < l));
    assert!(Dpd::new(10, gw_params(4, l), 1).is_err());
    assert!(Dpd::new(10, gw_params(3, 2.0), 1).is_err());
}
