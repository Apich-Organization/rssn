//! Chain regularisation, algorithmic regularisation and high-order time-transformed leapfrog.

use rssn::sim::chain::*;
use rssn::sim::symplectic::*;

fn kepler_exact(
    a: f64,
    e: f64,
    mu: f64,
    t: f64,
) -> [f64; 2] {
    // starts at apocentre (mean anomaly pi), counter-clockwise
    let n = (mu / (a * a * a)).sqrt();
    let m = std::f64::consts::PI + n * t;
    let mut ecc = m;
    for _ in 0..100 {
        ecc -= (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
    }
    [a * (ecc.cos() - e), a * (1.0 - e * e).sqrt() * ecc.sin()]
}

#[test]
fn chain_gbs_follows_a_highly_eccentric_kepler_orbit_to_round_off() {
    let (m1, m2, g) = (1.0, 0.5, 1.0);
    let (a, e): (f64, f64) = (1.0, 0.97);
    let mu = g * (m1 + m2);
    let r0 = a * (1.0 + e);
    let v0 = (mu * (1.0 - e) / (a * (1.0 + e))).sqrt();
    // relative orbit r = x2 - x1 ; body positions in the COM frame
    let f1 = m2 / (m1 + m2);
    let f2 = m1 / (m1 + m2);
    let pos = [[f1 * r0, 0.0, 0.0], [-f2 * r0, 0.0, 0.0]];
    let vel = [[0.0, f1 * v0, 0.0], [0.0, -f2 * v0, 0.0]];
    // body 0 is at +f1 r0 so the relative vector x1 - x2 = +r0: use that as the "relative" orbit
    let mut sys = ChainSystem::new(g, &[m1, m2], &pos, &vel).unwrap();
    let e0 = sys.energy();
    let period = 2.0 * std::f64::consts::PI * (a * a * a / mu).sqrt();
    let t_end = 2.7 * period;
    let steps = sys.evolve_until(t_end, 0.05, 6, 1e-13).unwrap();
    assert!(steps > 5);
    assert!((sys.t - t_end).abs() < 1e-11, "landed at {}", sys.t);
    let (p, _) = sys.cartesian();
    // x1 - x2 is the orbit mirrored through the origin of the exact solution started at apocentre on -x
    let ex = kepler_exact(a, e, mu, t_end);
    let rel = [p[0][0] - p[1][0], p[0][1] - p[1][1]];
    // the initial relative vector is +r0 x-hat and the velocity +y (counter-clockwise from +x apocentre): rotate by pi
    let err = ((rel[0] + ex[0]).powi(2) + (rel[1] + ex[1]).powi(2)).sqrt();
    // the exact reference starts with velocity -y at -x; rotating by pi maps it to +y at +x: rel = -ex
    assert!(err < 1e-8, "position error {err}");
    assert!(
        ((sys.energy() - e0) / e0).abs() < 1e-11,
        "energy {}",
        (sys.energy() - e0) / e0
    );
}

#[test]
fn pythagorean_three_body_problem_conserves_integrals_through_close_encounters() {
    let m = [3.0, 4.0, 5.0];
    let pos = [[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]];
    let vel = [[0.0; 3]; 3];
    let mut sys = ChainSystem::new(1.0, &m, &pos, &vel).unwrap();
    let e0 = sys.energy();
    let l0 = sys.angular_momentum();
    let mut min_sep = f64::INFINITY;
    let mut t = 0.0;
    while t < 12.0 {
        t += 0.05;
        sys.evolve_until(t, 0.02, 6, 1e-12).unwrap();
        let (p, _) = sys.cartesian();
        for i in 0..3 {
            for j in i + 1..3 {
                let d = ((p[i][0] - p[j][0]).powi(2)
                    + (p[i][1] - p[j][1]).powi(2)
                    + (p[i][2] - p[j][2]).powi(2))
                .sqrt();
                min_sep = min_sep.min(d);
            }
        }
    }
    assert!(min_sep < 0.15, "no close encounter sampled: {min_sep}");
    assert!(
        ((sys.energy() - e0) / e0).abs() < 1e-9,
        "energy {}",
        (sys.energy() - e0) / e0
    );
    let l = sys.angular_momentum();
    assert!(
        (l[2] - l0[2]).abs() < 1e-9 * l0[2].abs().max(1.0),
        "L {} vs {}",
        l[2],
        l0[2]
    );
}

#[test]
fn chain_selection_links_nearest_neighbours() {
    // a tight pair far from a third body: the chain must keep the pair adjacent
    let pos = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.001, 0.0, 0.0]];
    let sys = ChainSystem::new(1.0, &[1.0, 1.0, 1.0], &pos, &[[0.0; 3]; 3]).unwrap();
    let k1 = sys.order.iter().position(|&b| b == 1).unwrap();
    let k2 = sys.order.iter().position(|&b| b == 2).unwrap();
    assert_eq!(k1.abs_diff(k2), 1);
    let (p, _) = sys.cartesian();
    assert!((p[2][0] - p[1][0] - 0.001).abs() < 1e-12);
    assert!(ChainSystem::new(1.0, &[1.0], &[[0.0; 3]], &[[0.0; 3]]).is_err());
}

fn kepler_ttl_energy_error(
    order: usize,
    h: f64,
) -> f64 {
    // eccentric orbit in V = -1/r - 0.02/r^2 (a non-Kepler potential, so the energy error is nonzero)
    let pot = |r: f64| -1.0 / r - 0.02 / (r * r);
    let r0 = 0.4;
    let v0 = 1.1;
    let mut s = TtlState {
        q: vec![r0, 0.0, 0.0],
        p: vec![0.0, v0, 0.0],
        t: 0.0,
    };
    let e0 = 0.5 * v0 * v0 + pot(r0);
    let b = -e0;
    let mut grad = |q: &[f64], g: &mut [f64]| -> f64 {
        let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        let dv = 1.0 / (r * r) + 0.04 / (r * r * r);
        for k in 0..3 {
            g[k] = dv * q[k] / r;
        }
        pot(r)
    };
    let n = (4.0 / h).round() as usize;
    for _ in 0..n {
        ttl_step_composed(&mut s, &[1.0; 3], b, h, order, &mut grad).unwrap();
    }
    let r = (s.q[0] * s.q[0] + s.q[1] * s.q[1] + s.q[2] * s.q[2]).sqrt();
    let en = 0.5 * (s.p[0] * s.p[0] + s.p[1] * s.p[1] + s.p[2] * s.p[2]) + pot(r);
    (en - e0).abs()
}

#[test]
fn composed_time_transformed_leapfrog_has_the_requested_order() {
    for (order, min_ratio) in [(2usize, 3.0), (4, 11.0), (6, 40.0)] {
        let (h1, h2) = match order {
            | 6 => (0.5, 0.25),
            | 4 => (0.2, 0.1),
            | _ => (0.1, 0.05),
        };
        let e1 = kepler_ttl_energy_error(order, h1);
        let e2 = kepler_ttl_energy_error(order, h2);
        assert!(
            e1 / e2 > min_ratio,
            "order {order}: {e1} / {e2} = {}",
            e1 / e2
        );
    }
    // the fourth-order scheme is far more accurate than the plain TTL at the same step
    assert!(kepler_ttl_energy_error(4, 0.1) < 0.2 * kepler_ttl_energy_error(2, 0.1));
    assert!(composition_weights(3).is_err());
    let w = composition_weights(8).unwrap();
    assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-13);
}

#[test]
fn cartesian_ttl_gbs_extrapolation_reaches_round_off_accuracy() {
    let e = 0.9_f64;
    let r0 = 1.0 - e;
    let v0 = ((1.0 + e) / (1.0 - e)).sqrt();
    let e0 = 0.5 * v0 * v0 - 1.0 / r0;
    let mut s = TtlState {
        q: vec![r0, 0.0, 0.0],
        p: vec![0.0, v0, 0.0],
        t: 0.0,
    };
    let mut grad = |q: &[f64], g: &mut [f64]| -> f64 {
        let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        for k in 0..3 {
            g[k] = q[k] / (r * r * r);
        }
        -1.0 / r
    };
    let mut err = 0.0;
    for _ in 0..40 {
        err = ttl_gbs_step(&mut s, &[1.0; 3], -e0, 0.1, 6, &mut grad);
    }
    let r = (s.q[0] * s.q[0] + s.q[1] * s.q[1] + s.q[2] * s.q[2]).sqrt();
    let en = 0.5 * (s.p[0] * s.p[0] + s.p[1] * s.p[1] + s.p[2] * s.p[2]) - 1.0 / r;
    assert!((en - e0).abs() < 1e-11, "energy error {}", (en - e0).abs());
    assert!(err < 1e-9);
}
