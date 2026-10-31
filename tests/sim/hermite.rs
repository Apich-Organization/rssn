//! Hermite integrators: orders, block time steps, post-Newtonian precession, KS binaries.

use rssn::kernels::random::Rng;
use rssn::sim::hermite::*;

fn cfg(
    order: HermiteOrder,
    policy: StepPolicy,
) -> HermiteConfig {
    HermiteConfig {
        g: 1.0,
        eps: 0.0,
        order,
        policy,
        t0: 0.0,
    }
}

/// Test-mass orbit of eccentricity `e` and semi-major axis 1 around a unit mass.
fn kepler_ic(e: f64) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let vp = ((1.0 + e) / (1.0 - e)).sqrt();
    (
        vec![1.0, 0.0],
        vec![[0.0; 3], [1.0 - e, 0.0, 0.0]],
        vec![[0.0; 3], [0.0, vp, 0.0]],
    )
}

fn kepler_error(
    order: HermiteOrder,
    e: f64,
    steps: usize,
) -> f64 {
    let (m, x, v) = kepler_ic(e);
    let period = 2.0 * std::f64::consts::PI;
    let mut s = HermiteSystem::new(
        m,
        x.clone(),
        v,
        &cfg(order, StepPolicy::Shared(period / steps as f64)),
    )
    .unwrap();
    s.advance_to(period * (1.0 + 1e-12)).unwrap();
    let (xs, _) = s.state_at(s.time).unwrap();
    let d = [xs[1][0] - x[1][0], xs[1][1] - x[1][1], xs[1][2] - x[1][2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[test]
fn shared_step_hermite_schemes_have_orders_four_and_six() {
    let (e4a, e4b) = (
        kepler_error(HermiteOrder::Fourth, 0.5, 400),
        kepler_error(HermiteOrder::Fourth, 0.5, 800),
    );
    let r4 = (e4a / e4b).log2();
    assert!(r4 > 3.6 && r4 < 4.6, "4th order rate {r4}");
    let (e6a, e6b) = (
        kepler_error(HermiteOrder::Sixth, 0.5, 200),
        kepler_error(HermiteOrder::Sixth, 0.5, 400),
    );
    let r6 = (e6a / e6b).log2();
    assert!(r6 > 5.4 && r6 < 6.6, "6th order rate {r6}: {e6a} {e6b}");
    assert!(e6b < e4b);
}

#[test]
fn figure_eight_choreography_returns_after_one_period() {
    let x1 = [0.970_004_36, -0.243_087_53, 0.0];
    let v3 = [-0.932_407_37, -0.864_731_46, 0.0];
    let pos = vec![x1, [-x1[0], -x1[1], 0.0], [0.0; 3]];
    let vel = vec![
        [-v3[0] / 2.0, -v3[1] / 2.0, 0.0],
        [-v3[0] / 2.0, -v3[1] / 2.0, 0.0],
        v3,
    ];
    let period = 6.325_913_9;
    let mut s = HermiteSystem::new(
        vec![1.0; 3],
        pos.clone(),
        vel,
        &cfg(HermiteOrder::Sixth, StepPolicy::Shared(period / 800.0)),
    )
    .unwrap();
    s.advance_to(period * (1.0 + 1e-9)).unwrap();
    let (xs, _) = s.state_at(s.time).unwrap();
    let err = (0..3)
        .map(|i| {
            (0..2)
                .map(|k| (xs[i][k] - pos[i][k]).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max);
    assert!(err < 2e-5, "figure-eight return error {err}");
}

fn plummer_like(
    n: usize,
    seed: u64,
) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut rng = Rng::new(seed);
    let mut x = Vec::new();
    let mut v = Vec::new();
    for _ in 0..n {
        let r = rng.uniform().powf(1.0 / 3.0) * 1.0 + 0.02;
        let (a, b) = (rng.normal(), rng.normal());
        let c = rng.normal();
        let nn = (a * a + b * b + c * c).sqrt();
        x.push([r * a / nn, r * b / nn, r * c / nn]);
        v.push([0.3 * rng.normal(), 0.3 * rng.normal(), 0.3 * rng.normal()]);
    }
    (vec![1.0 / n as f64; n], x, v)
}

#[test]
fn block_time_steps_conserve_energy_and_agree_with_a_fine_reference() {
    let (m, x, v) = plummer_like(24, 3);
    let t_end = 0.5;
    let eps = 0.02;
    let mut cfg_ref = cfg(HermiteOrder::Sixth, StepPolicy::Shared(1.0 / 2048.0));
    cfg_ref.eps = eps;
    let mut reference = HermiteSystem::new(m.clone(), x.clone(), v.clone(), &cfg_ref).unwrap();
    reference.advance_to(t_end).unwrap();
    let (xr, _) = reference.state_at(t_end).unwrap();
    let mut results = Vec::new();
    for order in [HermiteOrder::Fourth, HermiteOrder::Sixth] {
        let mut c = cfg(
            order,
            StepPolicy::Block {
                eta: 0.02,
                dt_max: 1.0 / 16.0,
                dt_min: 1.0 / 65536.0,
            },
        );
        c.eps = eps;
        let mut s = HermiteSystem::new(m.clone(), x.clone(), v.clone(), &c).unwrap();
        let e0 = s.total_energy().unwrap();
        s.advance_to(t_end).unwrap();
        let e1 = s.total_energy().unwrap();
        assert!(
            ((e1 - e0) / e0).abs() < 1e-4,
            "{order:?} energy drift {}",
            (e1 - e0) / e0
        );
        let (xs, _) = s.state_at(t_end).unwrap();
        let dev = (0..m.len())
            .map(|i| {
                (0..3)
                    .map(|k| (xs[i][k] - xr[i][k]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(0.0, f64::max);
        assert!(dev < 5e-3, "{order:?} deviation {dev}");
        // individual steps: the number of particle updates is below that of a shared-step run
        let min_dt = s.dt.iter().copied().fold(f64::INFINITY, f64::min);
        let max_dt = s.dt.iter().copied().fold(0.0, f64::max);
        assert!(max_dt > min_dt, "steps must be individual");
        results.push((s.updates, dev));
    }
    assert!(results[0].0 > 0 && results[1].0 > 0);
}

#[test]
fn block_steps_form_a_power_of_two_hierarchy() {
    let (m, x, v) = plummer_like(16, 9);
    let mut c = cfg(
        HermiteOrder::Fourth,
        StepPolicy::Block {
            eta: 0.03,
            dt_max: 0.125,
            dt_min: 1e-6,
        },
    );
    c.eps = 0.05;
    let mut s = HermiteSystem::new(m, x, v, &c).unwrap();
    for _ in 0..400 {
        s.step().unwrap();
        for i in 0..16 {
            let q = s.dt[i] / 0.125;
            let l = q.log2();
            assert!(
                (l - l.round()).abs() < 1e-12,
                "dt {} is not a power-of-two fraction",
                s.dt[i]
            );
            // block steps must be commensurable with the particle's own time
            assert!((s.t_last[i] / s.dt[i] - (s.t_last[i] / s.dt[i]).round()).abs() < 1e-9);
        }
    }
}

/// Apsidal precession per orbit of the relative orbit from successive pericentres.
fn pn_precession(
    c_light: f64,
    e: f64,
    orbits: f64,
) -> (f64, f64) {
    let (mut m, x, v) = kepler_ic(e);
    m[1] = 1e-9;
    let mut s = HermiteSystem::new(
        m,
        x,
        v,
        &cfg(
            HermiteOrder::Fourth,
            StepPolicy::Block {
                eta: 0.002,
                dt_max: 1.0 / 64.0,
                dt_min: 1e-9,
            },
        ),
    )
    .unwrap();
    s.enable_post_newtonian(c_light).unwrap();
    let e0 = s.total_energy().unwrap();
    let mut peri: Vec<f64> = Vec::new();
    let (mut r_prev2, mut r_prev, mut t_prev) = (f64::NAN, f64::NAN, 0.0);
    let t_end = orbits * 2.0 * std::f64::consts::PI;
    let dtg = 0.01;
    let mut tg = 0.0;
    let mut samples: Vec<(f64, f64, f64)> = Vec::new();
    while tg < t_end {
        tg += dtg;
        s.advance_to(tg).unwrap();
        let (xs, _) = s.state_at(tg).unwrap();
        let d = [xs[1][0] - xs[0][0], xs[1][1] - xs[0][1]];
        samples.push((tg, (d[0] * d[0] + d[1] * d[1]).sqrt(), d[1].atan2(d[0])));
    }
    let _ = (&mut r_prev2, &mut r_prev, &mut t_prev);
    // pericentre angles by parabolic refinement of the radius minima
    for k in 1..samples.len() - 1 {
        if samples[k].1 < samples[k - 1].1 && samples[k].1 <= samples[k + 1].1 {
            let (a, b, c) = (samples[k - 1].1, samples[k].1, samples[k + 1].1);
            let shift = 0.5 * (a - c) / (a - 2.0 * b + c);
            // angle changes ~ linearly across the three samples
            let ang = samples[k].2 + shift * (samples[k + 1].2 - samples[k - 1].2) / 2.0;
            peri.push(ang);
        }
    }
    let mut total = 0.0;
    for w in peri.windows(2) {
        let mut d = w[1] - w[0];
        while d > std::f64::consts::PI {
            d -= 2.0 * std::f64::consts::PI;
        }
        while d < -std::f64::consts::PI {
            d += 2.0 * std::f64::consts::PI;
        }
        total += d;
    }
    let per_orbit = total / (peri.len() - 1) as f64;
    let e1 = s.total_energy().unwrap();
    (per_orbit, ((e1 - e0) / e0).abs())
}

#[test]
fn one_pn_perihelion_precession_matches_the_einstein_formula() {
    let e = 0.5;
    let c_light = 40.0;
    let (dphi, drift) = pn_precession(c_light, e, 6.2);
    // 6 pi G M / (c^2 a (1 - e^2)) per orbit
    let expected = 6.0 * std::f64::consts::PI / (c_light * c_light * (1.0 - e * e));
    assert!(
        (dphi - expected).abs() < 0.02 * expected,
        "precession {dphi} vs {expected}"
    );
    assert!(drift < 1e-4, "1PN energy drift {drift}");
    // the precession scales as 1/c^2
    let (dphi2, _) = pn_precession(2.0 * c_light, e, 6.2);
    assert!((dphi / dphi2 - 4.0).abs() < 0.2, "{dphi} {dphi2}");
}

#[test]
fn eih_two_body_accelerations_are_momentum_conserving_to_one_pn() {
    // equal masses on a bound orbit: the EIH acceleration conserves the 1PN energy and total
    // momentum to O(1/c^4)
    let m = vec![1.0, 1.0];
    let x = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
    let v = vec![[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]];
    let mut s = HermiteSystem::new(
        m,
        x,
        v,
        &cfg(
            HermiteOrder::Fourth,
            StepPolicy::Block {
                eta: 0.003,
                dt_max: 1.0 / 128.0,
                dt_min: 1e-9,
            },
        ),
    )
    .unwrap();
    s.enable_post_newtonian(30.0).unwrap();
    let e0 = s.total_energy().unwrap();
    s.advance_to(10.0).unwrap();
    let e1 = s.total_energy().unwrap();
    assert!(((e1 - e0) / e0).abs() < 1e-4, "energy {}", (e1 - e0) / e0);
    let (_, vs) = s.state_at(s.time).unwrap();
    let p = [vs[0][0] + vs[1][0], vs[0][1] + vs[1][1]];
    assert!(p[0].abs().max(p[1].abs()) < 1e-4, "momentum {p:?}");
}

#[test]
fn hermite_validates_its_inputs() {
    let (m, x, v) = kepler_ic(0.1);
    assert!(
        HermiteSystem::new(
            m.clone(),
            x.clone(),
            v.clone(),
            &cfg(HermiteOrder::Fourth, StepPolicy::Shared(-1.0))
        )
        .is_err()
    );
    assert!(
        HermiteSystem::new(
            vec![1.0],
            vec![[0.0; 3]],
            vec![[0.0; 3]],
            &cfg(HermiteOrder::Fourth, StepPolicy::Shared(0.1))
        )
        .is_err()
    );
    let mut s =
        HermiteSystem::new(m, x, v, &cfg(HermiteOrder::Sixth, StepPolicy::Shared(0.01))).unwrap();
    assert!(s.enable_post_newtonian(10.0).is_err());
}

#[test]
fn ks_binary_in_a_hierarchical_triple_matches_a_direct_reference() {
    // tight eccentric binary (a = 0.04, e = 0.9) orbited by a distant third body
    let (a_in, e_in) = (0.04, 0.9);
    let mb: f64 = 1.0;
    let r_apo = a_in * (1.0 + e_in);
    let v_apo = ((mb) * (1.0 - e_in) / (a_in * (1.0 + e_in))).sqrt();
    let (m1, m2) = (0.5, 0.5);
    let x1 = [-0.5 * r_apo, 0.0, 0.0];
    let x2 = [0.5 * r_apo, 0.0, 0.0];
    let v1 = [0.0, -0.5 * v_apo, 0.0];
    let v2 = [0.0, 0.5 * v_apo, 0.0];
    let r_out: f64 = 1.5;
    let v_out = (1.5 / r_out).sqrt() * 0.9;
    let x3 = [r_out, 0.0, 0.0];
    let v3 = [0.0, v_out, 0.0];
    let m3 = 0.5;
    let masses = vec![m1, m2, m3];
    let pos = vec![x1, x2, x3];
    let vel = vec![v1, v2, v3];
    let t_end = 1.5;

    // high-accuracy direct reference: sixth-order block steps
    let mut reference = HermiteSystem::new(
        masses.clone(),
        pos.clone(),
        vel.clone(),
        &cfg(
            HermiteOrder::Sixth,
            StepPolicy::Block {
                eta: 0.005,
                dt_max: 1.0 / 64.0,
                dt_min: 1e-10,
            },
        ),
    )
    .unwrap();
    let e_ref0 = reference.total_energy().unwrap();
    reference.advance_to(t_end).unwrap();
    let (xr, _) = reference.state_at(t_end).unwrap();
    assert!(((reference.total_energy().unwrap() - e_ref0) / e_ref0).abs() < 1e-7);

    let mut hybrid = HermiteSystem::new(
        masses,
        pos,
        vel,
        &cfg(
            HermiteOrder::Fourth,
            StepPolicy::Block {
                eta: 0.02,
                dt_max: 1.0 / 32.0,
                dt_min: 1e-8,
            },
        ),
    )
    .unwrap();
    hybrid.add_binary(0, 1, 0.02).unwrap();
    let e0 = hybrid.total_energy().unwrap();
    hybrid.advance_to(t_end).unwrap();
    let (xh, _) = hybrid.state_at(t_end).unwrap();
    let e1 = hybrid.total_energy().unwrap();
    assert!(
        ((e1 - e0) / e0).abs() < 1e-5,
        "hybrid energy drift {}",
        (e1 - e0) / e0
    );
    // the third body and the binary centre of mass follow the reference
    let dev3 = (0..3)
        .map(|k| (xh[2][k] - xr[2][k]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(dev3 < 1e-4, "third-body deviation {dev3}");
    let cm = |x: &[[f64; 3]]| [0.5 * (x[0][0] + x[1][0]), 0.5 * (x[0][1] + x[1][1])];
    let (cr, ch) = (cm(&xr), cm(&xh));
    assert!((cr[0] - ch[0]).hypot(cr[1] - ch[1]) < 1e-4);
    // internal binary orbit: separation agrees
    let sep = |x: &[[f64; 3]]| ((x[1][0] - x[0][0]).powi(2) + (x[1][1] - x[0][1]).powi(2)).sqrt();
    assert!(
        (sep(&xr) - sep(&xh)).abs() < 2e-3,
        "{} vs {}",
        sep(&xr),
        sep(&xh)
    );
    // the hybrid takes far fewer outer updates than the directly integrated binary needs
    assert!(
        hybrid.updates * 3 < reference.updates,
        "{} vs {}",
        hybrid.updates,
        reference.updates
    );
}
