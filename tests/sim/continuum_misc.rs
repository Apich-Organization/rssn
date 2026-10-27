//! Balsara switch of the SPH viscosity and the complex absorbing potential of the split-step solver.

use num_complex::Complex;
use rssn::sim::split_step::{SplitOrder, SplitStep};
use rssn::sim::sph::{Sph, SphBoundary};
use std::f64::consts::PI;

fn shear_gas(balsara: bool) -> f64 {
    // periodic lattice with a divergence-free shear v_x = A sin(2 pi y): dissipation = viscous heating rate
    let n = 20;
    let dx = 1.0 / n as f64;
    let mut pos = Vec::new();
    for i in 0..n {
        for j in 0..n {
            pos.push([(i as f64 + 0.5) * dx, (j as f64 + 0.5) * dx]);
        }
    }
    let m = 1.0 / (n * n) as f64;
    let mut s = Sph::<2>::new(pos, vec![m; n * n], vec![2.5; n * n], [0.0; 2], [1.0; 2], [SphBoundary::Periodic; 2]);
    s.eta = 1.4;
    s.balsara = balsara;
    for i in 0..n * n {
        s.vel[i][0] = 0.3 * (2.0 * PI * s.pos[i][1]).sin();
    }
    let (_, du) = s.derivatives();
    s.mass.iter().zip(&du).map(|(m, d)| m * d).sum()
}

#[test]
fn balsara_switch_suppresses_the_spurious_shear_viscosity() {
    let plain = shear_gas(false);
    let switched = shear_gas(true);
    assert!(plain > 0.0, "the plain Monaghan viscosity heats a shear flow: {plain}");
    assert!(switched < 0.25 * plain, "Balsara heating {switched} vs plain {plain}");
}

#[test]
fn balsara_switch_leaves_a_compression_untouched() {
    // 1D compression: |curl v| = 0, so the switch is close to one where the flow converges
    let heating = |balsara: bool| {
        let n = 80;
        let dx = 1.0 / n as f64;
        let pos: Vec<[f64; 1]> = (0..n).map(|i| [(i as f64 + 0.5) * dx]).collect();
        let mut s = Sph::<1>::new(pos, vec![dx; n], vec![2.5; n], [0.0], [1.0], [SphBoundary::Periodic]);
        s.balsara = balsara;
        for i in 0..n {
            // converging flow towards x = 0.5 (periodic-compatible hat profile)
            let x = s.pos[i][0];
            s.vel[i][0] = -0.5 * (2.0 * PI * (x - 0.5)).sin();
        }
        let (_, du) = s.derivatives();
        s.mass.iter().zip(&du).map(|(m, d)| m * d).sum::<f64>()
    };
    let (a, b) = (heating(false), heating(true));
    // the switch only reduces the viscosity (where the divergence vanishes), never in the compressing regions
    assert!(a > 0.0 && b <= a * (1.0 + 1e-9) && b > 0.8 * a, "plain {a} vs Balsara {b}");
}

#[test]
fn absorbing_potential_removes_an_outgoing_wave_packet_that_would_otherwise_wrap_around() {
    let ss = SplitStep::centered(&[512], &[80.0], 1.0).unwrap();
    let k0 = 4.0;
    let psi0 = ss.sample(|x| Complex::from_polar((-(x[0] + 10.0).powi(2) / 8.0).exp(), k0 * x[0]));
    let n0 = ss.norm2(&psi0);
    let free = |x: f64| 0.0 * x;
    // free periodic evolution is unitary: the packet leaves at +40 and re-enters at -40
    let t_end = 16.0;
    let steps = 1600;
    let dt = t_end / steps as f64;
    let wrapped = ss.evolve(&psi0, &|_, x| free(x[0]), 0.0, 0.0, dt, steps, SplitOrder::Strang);
    assert!((ss.norm2(&wrapped) / n0 - 1.0).abs() < 1e-10);
    // quadratic absorbing layers of width 10 at both ends, strength 2
    let cap = |x: [f64; 3]| {
        let d = (x[0].abs() - 30.0).max(0.0) / 10.0;
        2.0 * d * d
    };
    let absorbed = ss.evolve_absorbing(&psi0, &|_, x| free(x[0]), &cap, 0.0, 0.0, dt, steps, SplitOrder::Strang);
    let left = ss.norm2(&absorbed) / n0;
    assert!(left < 5e-3, "remaining norm {left}");
    // the absorber does not act in the interior: before the packet reaches the layer the norm is unchanged
    let early = ss.evolve_absorbing(&psi0, &|_, x| free(x[0]), &cap, 0.0, 0.0, dt, 400, SplitOrder::Strang);
    assert!((ss.norm2(&early) / n0 - 1.0).abs() < 1e-6, "interior norm {}", ss.norm2(&early) / n0);
}
