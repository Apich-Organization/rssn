//! Composition integrators, time-transformed leapfrog, KS regularisation, rigid bodies.

use rssn::sim::symplectic::*;

fn pendulum(
    q: &[f64],
    a: &mut [f64],
) {
    a[0] = -q[0].sin();
}

fn pend_error(
    scheme: Scheme,
    h: f64,
    reference: (f64, f64),
) -> f64 {
    let t_end = 4.0;
    let n = (t_end / h).round() as usize;
    let mut q = [1.5];
    let mut p = [0.0];
    integrate(scheme, &mut q, &mut p, h, n, &mut pendulum);
    ((q[0] - reference.0).powi(2) + (p[0] - reference.1).powi(2)).sqrt()
}

#[test]
fn schemes_have_their_nominal_order() {
    let mut q = [1.5];
    let mut p = [0.0];
    integrate(Scheme::Yoshida8, &mut q, &mut p, 0.01, 400, &mut pendulum);
    let reference = (q[0], p[0]);
    for (scheme, h) in [
        (Scheme::Verlet, 0.04),
        (Scheme::ForestRuth, 0.1),
        (Scheme::Yoshida4, 0.1),
        (Scheme::BlanesMoan4, 0.1),
        (Scheme::Yoshida6, 0.25),
        (Scheme::Yoshida8, 0.5),
    ] {
        let e1 = pend_error(scheme, h, reference);
        let e2 = pend_error(scheme, h / 2.0, reference);
        let order = (e1 / e2).log2();
        let nominal = scheme.order() as f64;
        assert!(
            (order - nominal).abs() < 0.6,
            "{scheme:?}: observed order {order} (errors {e1:e}, {e2:e}), nominal {nominal}"
        );
    }
}

#[test]
fn stage_tables_are_consistent() {
    for s in [Scheme::Verlet, Scheme::ForestRuth, Scheme::Yoshida4, Scheme::Yoshida6, Scheme::Yoshida8, Scheme::BlanesMoan4] {
        let st = s.stages();
        let drift: f64 = st.iter().map(|x| x.0).sum();
        let kick: f64 = st.iter().map(|x| x.1).sum();
        assert!((drift - 1.0).abs() < 1e-13 && (kick - 1.0).abs() < 1e-13, "{s:?}");
    }
    // the triple jump lifts order 2 to 4 with the same table as Yoshida4
    let lifted = triple_jump(&Scheme::Verlet.stages(), 2);
    assert_eq!(lifted.len(), Scheme::Yoshida4.stages().len());
}

fn kepler_acc(
    q: &[f64],
    a: &mut [f64],
) {
    let r3 = (q[0] * q[0] + q[1] * q[1]).powf(1.5);
    a[0] = -q[0] / r3;
    a[1] = -q[1] / r3;
}

/// Kepler orbit with a = 1 and eccentricity e; returns position at time t from pericentre.
fn kepler_exact(
    e: f64,
    t: f64,
) -> [f64; 2] {
    let m = t;
    let mut ea = m;
    for _ in 0..100 {
        ea -= (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
    }
    [ea.cos() - e, (1.0 - e * e).sqrt() * ea.sin()]
}

#[test]
fn kepler_period_and_energy_with_symplectic_scheme() {
    let e: f64 = 0.3;
    let period = 2.0 * std::f64::consts::PI;
    let mut q = [1.0 - e, 0.0];
    let mut p = [0.0, ((1.0 + e) / (1.0 - e)).sqrt()];
    let n = 200;
    let h = period / n as f64;
    let e0 = 0.5 * (p[0] * p[0] + p[1] * p[1]) - 1.0 / (q[0] * q[0] + q[1] * q[1]).sqrt();
    let mut emax = 0.0_f64;
    let st = Scheme::Yoshida6.stages();
    for _ in 0..n {
        step(&st, &mut q, &mut p, h, &mut kepler_acc);
        let en = 0.5 * (p[0] * p[0] + p[1] * p[1]) - 1.0 / (q[0] * q[0] + q[1] * q[1]).sqrt();
        emax = emax.max((en - e0).abs());
    }
    assert!((q[0] - (1.0 - e)).abs() < 1e-6 && q[1].abs() < 1e-6, "period error {q:?}");
    assert!(emax < 1e-6, "energy error {emax}");
}

#[test]
fn ttl_beats_fixed_step_on_an_eccentric_orbit() {
    let e: f64 = 0.9;
    let n = 400;
    let period = 2.0 * std::f64::consts::PI;
    let apo = [-(1.0 + e), 0.0];
    let vapo = [0.0, -((1.0 - e) / (1.0 + e)).sqrt()];
    // exact orbit started at apocentre = pericentre orbit shifted by half a period; mirror coordinates
    let exact_at = |t: f64| kepler_exact(e, t + period / 2.0);
    // fixed step Verlet
    let mut q = apo.to_vec();
    let mut p = vapo.to_vec();
    let h = period / n as f64;
    integrate(Scheme::Verlet, &mut q, &mut p, h, n, &mut kepler_acc);
    let ex = exact_at(h * n as f64);
    let err_fixed = ((q[0] - ex[0]).powi(2) + (q[1] - ex[1]).powi(2)).sqrt();
    // TTL
    let mut s = TtlState { q: apo.to_vec(), p: vapo.to_vec(), t: 0.0 };
    let b = 0.5;
    let mass = [1.0, 1.0];
    let mut pot = |q: &[f64], g: &mut [f64]| {
        let r = (q[0] * q[0] + q[1] * q[1]).sqrt();
        g[0] = q[0] / (r * r * r);
        g[1] = q[1] / (r * r * r);
        -1.0 / r
    };
    for _ in 0..n {
        ttl_step(&mut s, &mass, b, period / n as f64, &mut pot);
    }
    let ex = exact_at(s.t);
    let err_ttl = ((s.q[0] - ex[0]).powi(2) + (s.q[1] - ex[1]).powi(2)).sqrt();
    assert!(err_ttl < 0.02, "TTL error {err_ttl}");
    assert!(err_fixed > 10.0 * err_ttl, "fixed {err_fixed} vs ttl {err_ttl}");
    assert!((s.t - period).abs() < 0.1, "physical time {}", s.t);
}

#[test]
fn ks_regularisation_recovers_a_near_collision_orbit_exactly() {
    // a = 1, e = 0.99, started at apocentre
    let e: f64 = 0.99;
    let x0 = [1.0 + e, 0.0, 0.0];
    let v0 = [0.0, ((1.0 - e) / (1.0 + e)).sqrt(), 0.0];
    let n = 200;
    let ds = 2.0 * std::f64::consts::PI / n as f64;
    let none = |_: &[f64; 3], _: &[f64; 3], _: f64| [0.0; 3];
    let r = ks_propagate(&x0, &v0, &none, ds, n);
    assert!((r.t - 2.0 * std::f64::consts::PI).abs() < 1e-6, "period {}", r.t);
    let (x, v) = ks_to_cartesian(&r.u, &r.up);
    for k in 0..3 {
        assert!((x[k] - x0[k]).abs() < 1e-6, "x {x:?}");
        assert!((v[k] - v0[k]).abs() < 1e-6, "v {v:?}");
    }
    assert!((r.energy + 0.5).abs() < 1e-12);
}

#[test]
fn ks_roundtrip_and_perturbed_orbit_matches_cartesian_integration() {
    let x0 = [-0.7, 0.4, 0.2];
    let v0 = [0.1, -0.6, 0.5];
    let (u, up) = ks_from_cartesian(&x0, &v0);
    let (x, v) = ks_to_cartesian(&u, &up);
    for k in 0..3 {
        assert!((x[k] - x0[k]).abs() < 1e-13 && (v[k] - v0[k]).abs() < 1e-13);
    }
    let pert = |_: &[f64; 3], v: &[f64; 3], _: f64| [1e-2 * v[0], 5e-3, -2e-2 * v[2]];
    let res = ks_propagate(&x0, &v0, &pert, 0.002, 1500);
    let (xk, vk) = ks_to_cartesian(&res.u, &res.up);
    // Cartesian reference integrated to the same physical time with tiny RK4 steps.
    let f = |y: &[f64; 6]| {
        let r3 = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).powf(1.5);
        let p = pert(&[y[0], y[1], y[2]], &[y[3], y[4], y[5]], 0.0);
        [y[3], y[4], y[5], -y[0] / r3 + p[0], -y[1] / r3 + p[1], -y[2] / r3 + p[2]]
    };
    let mut y = [x0[0], x0[1], x0[2], v0[0], v0[1], v0[2]];
    let nstep = 40_000;
    let dt = res.t / nstep as f64;
    for _ in 0..nstep {
        let k1 = f(&y);
        let mut t = y;
        for i in 0..6 {
            t[i] = y[i] + 0.5 * dt * k1[i];
        }
        let k2 = f(&t);
        for i in 0..6 {
            t[i] = y[i] + 0.5 * dt * k2[i];
        }
        let k3 = f(&t);
        for i in 0..6 {
            t[i] = y[i] + dt * k3[i];
        }
        let k4 = f(&t);
        for i in 0..6 {
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    for k in 0..3 {
        assert!((xk[k] - y[k]).abs() < 1e-6, "x {xk:?} vs {y:?}");
        assert!((vk[k] - y[3 + k]).abs() < 1e-6, "v {vk:?} vs {y:?}");
    }
    // energy variable matches the directly evaluated energy
    let r = (xk[0] * xk[0] + xk[1] * xk[1] + xk[2] * xk[2]).sqrt();
    let e = 0.5 * (vk[0] * vk[0] + vk[1] * vk[1] + vk[2] * vk[2]) - 1.0 / r;
    assert!((e - res.energy).abs() < 1e-7);
}

fn top_energy_error(h: f64) -> f64 {
    let mut b = RigidBody { inertia: [1.0, 2.0, 3.0], q: [1.0, 0.0, 0.0, 0.0], l_body: [0.3, 1.0, 0.4] };
    let e0 = b.kinetic_energy();
    let mut worst = 0.0_f64;
    let n = (5.0 / h).round() as usize;
    for _ in 0..n {
        b.free_step(h);
        worst = worst.max((b.kinetic_energy() - e0).abs());
    }
    worst
}

#[test]
fn rigid_top_conserves_angular_momentum_and_has_second_order_energy_error() {
    let mut b = RigidBody { inertia: [1.0, 2.0, 3.0], q: [0.5, 0.5, 0.5, 0.5], l_body: [0.3, 1.0, 0.4] };
    let l0 = b.l_space();
    for _ in 0..2000 {
        b.free_step(0.01);
    }
    let l = b.l_space();
    for k in 0..3 {
        assert!((l[k] - l0[k]).abs() < 1e-11, "L drift {l:?} vs {l0:?}");
    }
    let nq = b.q.iter().map(|v| v * v).sum::<f64>();
    assert!((nq - 1.0).abs() < 1e-13);
    let e1 = top_energy_error(0.02);
    let e2 = top_energy_error(0.01);
    let order = (e1 / e2).log2();
    assert!((order - 2.0).abs() < 0.3, "order {order} ({e1:e}, {e2:e})");
}

#[test]
fn symmetric_top_precession_matches_analytic_rate() {
    let i = [1.0, 1.0, 2.0];
    let mut b = RigidBody { inertia: i, q: [1.0, 0.0, 0.0, 0.0], l_body: [0.3, 0.0, 0.5] };
    let t = 3.0;
    let n = 3000;
    for _ in 0..n {
        b.free_step(t / n as f64);
    }
    let omega_b = 0.5 * (1.0 / i[0] - 1.0 / i[2]);
    let exact = [0.3 * (omega_b * t).cos(), 0.3 * (omega_b * t).sin()];
    assert!((b.l_body[0] - exact[0]).abs() < 1e-5 && (b.l_body[1] - exact[1]).abs() < 1e-5, "{:?} vs {exact:?}", b.l_body);
    assert!((b.l_body[2] - 0.5).abs() < 1e-6);
}

#[test]
fn rigid_body_step_is_time_reversible_under_torque() {
    let torque = |q: &[f64; 4]| {
        let r = quat_to_matrix(q);
        [0.2 * r[2][1], -0.3 * r[2][0], 0.1]
    };
    let mut b = RigidBody { inertia: [1.0, 1.5, 2.5], q: [0.9, 0.1, -0.3, 0.3], l_body: [0.1, -0.2, 0.3] };
    let n0 = b.q.iter().map(|v| v * v).sum::<f64>().sqrt();
    for v in &mut b.q {
        *v /= n0;
    }
    let start = b.clone();
    for _ in 0..200 {
        b.step(0.02, &torque);
    }
    for _ in 0..200 {
        b.step(-0.02, &torque);
    }
    for k in 0..4 {
        assert!((b.q[k] - start.q[k]).abs() < 1e-10, "{:?} vs {:?}", b.q, start.q);
    }
    for k in 0..3 {
        assert!((b.l_body[k] - start.l_body[k]).abs() < 1e-10);
    }
}

/// High-accuracy RK4 reference of the free rigid body (quaternion + Euler equations).
fn rigid_reference(
    inertia: [f64; 3],
    q0: [f64; 4],
    l0: [f64; 3],
    t: f64,
) -> ([f64; 4], [f64; 3]) {
    let rhs = |q: &[f64; 4], l: &[f64; 3]| -> ([f64; 4], [f64; 3]) {
        let w = [l[0] / inertia[0], l[1] / inertia[1], l[2] / inertia[2]];
        let dq = quat_mul(q, &[0.0, w[0], w[1], w[2]]);
        let dl = [l[1] * w[2] - l[2] * w[1], l[2] * w[0] - l[0] * w[2], l[0] * w[1] - l[1] * w[0]];
        ([0.5 * dq[0], 0.5 * dq[1], 0.5 * dq[2], 0.5 * dq[3]], dl)
    };
    let n = 40_000;
    let h = t / n as f64;
    let (mut q, mut l) = (q0, l0);
    for _ in 0..n {
        let (k1q, k1l) = rhs(&q, &l);
        let adv = |c: f64, kq: &[f64; 4], kl: &[f64; 3]| -> ([f64; 4], [f64; 3]) {
            ([q[0] + c * kq[0], q[1] + c * kq[1], q[2] + c * kq[2], q[3] + c * kq[3]], [l[0] + c * kl[0], l[1] + c * kl[1], l[2] + c * kl[2]])
        };
        let (q2, l2) = adv(0.5 * h, &k1q, &k1l);
        let (k2q, k2l) = rhs(&q2, &l2);
        let (q3, l3) = adv(0.5 * h, &k2q, &k2l);
        let (k3q, k3l) = rhs(&q3, &l3);
        let (q4, l4) = adv(h, &k3q, &k3l);
        let (k4q, k4l) = rhs(&q4, &l4);
        for i in 0..4 {
            q[i] += h / 6.0 * (k1q[i] + 2.0 * k2q[i] + 2.0 * k3q[i] + k4q[i]);
        }
        for i in 0..3 {
            l[i] += h / 6.0 * (k1l[i] + 2.0 * k2l[i] + 2.0 * k3l[i] + k4l[i]);
        }
    }
    (q, l)
}

fn rigid_error(
    inertia: [f64; 3],
    top: bool,
    order: usize,
    h: f64,
) -> f64 {
    let q0 = [0.5, 0.5, 0.5, 0.5];
    let l0 = [0.7, 0.4, 1.1];
    let t = 2.0;
    let mut b = RigidBody { inertia, q: q0, l_body: l0 };
    for _ in 0..(t / h).round() as usize {
        b.free_step_composed(h, top, order).unwrap();
    }
    let (qr, lr) = rigid_reference(inertia, q0, l0, t);
    let dq: f64 = (0..4).map(|i| (b.q[i] - qr[i]).powi(2)).sum::<f64>().sqrt();
    let dl: f64 = (0..3).map(|i| (b.l_body[i] - lr[i]).powi(2)).sum::<f64>().sqrt();
    dq + dl
}

#[test]
fn symmetric_top_flow_is_exact_for_a_symmetric_body_for_any_step() {
    let inertia = [2.0, 2.0, 1.0];
    for axis_perm in [[2.0, 2.0, 1.0], [1.0, 2.0, 2.0], [2.0, 1.0, 2.0]] {
        let e = rigid_error(axis_perm, true, 2, 2.0); // a single step of length 2
        assert!(e < 1e-9, "{axis_perm:?}: {e}");
    }
    let _ = inertia;
    // the DLM splitting is not exact for the same single step
    assert!(rigid_error([2.0, 2.0, 1.0], false, 2, 2.0) > 1e-3);
}

#[test]
fn symmetric_top_splitting_is_second_order_and_composes_to_higher_order() {
    let inertia = [1.0, 1.15, 3.0];
    let e1 = rigid_error(inertia, true, 2, 0.1);
    let e2 = rigid_error(inertia, true, 2, 0.05);
    assert!(e1 / e2 > 3.5 && e1 / e2 < 4.5, "second order: {}", e1 / e2);
    let f1 = rigid_error(inertia, true, 4, 0.2);
    let f2 = rigid_error(inertia, true, 4, 0.1);
    assert!(f1 / f2 > 12.0, "fourth order: {}", f1 / f2);
    // for a nearly symmetric body it beats the three-axis splitting at the same step
    assert!(e1 < 0.2 * rigid_error(inertia, false, 2, 0.1));
    let mut b = RigidBody { inertia, q: [1.0, 0.0, 0.0, 0.0], l_body: [0.4, 0.5, 0.6] };
    let (l0, k0) = (b.l_space(), b.kinetic_energy());
    for _ in 0..200 {
        b.free_step_top(0.1);
    }
    let l = b.l_space();
    assert!((0..3).all(|i| (l[i] - l0[i]).abs() < 1e-12));
    assert!((b.kinetic_energy() - k0).abs() < 1e-3 * k0);
    assert!((b.q.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-14);
    assert!(b.free_step_composed(0.1, true, 3).is_err());
}

#[test]
fn space_torque_kick_matches_the_body_frame_torque() {
    let mut b = RigidBody { inertia: [1.0, 2.0, 3.0], q: [0.5, 0.5, 0.5, 0.5], l_body: [0.0; 3] };
    b.kick_space(0.5, &[1.0, 0.0, 0.0]);
    // a unit-torque kick along space x: L_space = (0.5, 0, 0)
    let l = b.l_space();
    assert!((l[0] - 0.5).abs() < 1e-14 && l[1].abs() < 1e-14 && l[2].abs() < 1e-14);
}

/// Reference integration of the free rigid body by classical RK4 with a tiny step.
fn rk4_free_body(b: &RigidBody, t_end: f64, n: usize) -> RigidBody {
    let f = |q: &[f64; 4], l: &[f64; 3]| {
        let w = [l[0] / b.inertia[0], l[1] / b.inertia[1], l[2] / b.inertia[2]];
        let dq = quat_mul(q, &[0.0, w[0], w[1], w[2]]);
        (
            [0.5 * dq[0], 0.5 * dq[1], 0.5 * dq[2], 0.5 * dq[3]],
            [l[1] * w[2] - l[2] * w[1], l[2] * w[0] - l[0] * w[2], l[0] * w[1] - l[1] * w[0]],
        )
    };
    let h = t_end / n as f64;
    let (mut q, mut l) = (b.q, b.l_body);
    for _ in 0..n {
        let (k1q, k1l) = f(&q, &l);
        let adv = |q: &[f64; 4], l: &[f64; 3], kq: &[f64; 4], kl: &[f64; 3], s: f64| {
            ([q[0] + s * kq[0], q[1] + s * kq[1], q[2] + s * kq[2], q[3] + s * kq[3]], [l[0] + s * kl[0], l[1] + s * kl[1], l[2] + s * kl[2]])
        };
        let (q2, l2) = adv(&q, &l, &k1q, &k1l, 0.5 * h);
        let (k2q, k2l) = f(&q2, &l2);
        let (q3, l3) = adv(&q, &l, &k2q, &k2l, 0.5 * h);
        let (k3q, k3l) = f(&q3, &l3);
        let (q4, l4) = adv(&q, &l, &k3q, &k3l, h);
        let (k4q, k4l) = f(&q4, &l4);
        for i in 0..4 {
            q[i] += h / 6.0 * (k1q[i] + 2.0 * k2q[i] + 2.0 * k3q[i] + k4q[i]);
        }
        for i in 0..3 {
            l[i] += h / 6.0 * (k1l[i] + 2.0 * k2l[i] + 2.0 * k3l[i] + k4l[i]);
        }
    }
    let nrm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
    RigidBody { inertia: b.inertia, q: q.map(|x| x / nrm), l_body: l }
}

fn quat_dist(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    let d: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    1.0 - d.abs()
}

#[test]
fn exact_asymmetric_top_matches_high_accuracy_reference() {
    // rotation near the largest axis, near the smallest axis, negative L_d, general orientation
    let cases: [([f64; 3], [f64; 3]); 5] = [
        ([1.0, 2.0, 3.0], [0.3, 1.0, 0.4]),
        ([1.0, 2.0, 3.0], [0.3, 0.8, 1.9]),
        ([1.0, 2.0, 3.0], [1.4, 0.6, -0.5]),
        ([0.7, 1.9, 2.3], [-0.4, 0.9, 0.7]),
        ([1.0, 2.0, 3.0], [(1.0f64 / 3.0).sqrt() * 0.99, 1.0, 0.43]),
    ];
    for (inertia, l_body) in cases {
        let b0 = RigidBody { inertia, q: [0.5, 0.5, -0.5, 0.5], l_body };
        for &t in &[0.37, 3.0, 12.0, -5.0] {
            let reference = rk4_free_body(&b0, t, (t.abs() * 20_000.0) as usize);
            let mut b = b0.clone();
            b.free_flow_exact(t);
            for k in 0..3 {
                assert!((b.l_body[k] - reference.l_body[k]).abs() < 1e-8, "L {inertia:?} {l_body:?} t={t}");
            }
            assert!(quat_dist(&b.q, &reference.q) < 1e-9, "q {inertia:?} {l_body:?} t={t}: {}", quat_dist(&b.q, &reference.q));
        }
    }
}

#[test]
fn exact_flow_is_a_group_and_conserves_invariants() {
    let b0 = RigidBody { inertia: [1.0, 2.0, 3.5], q: [1.0, 0.0, 0.0, 0.0], l_body: [0.5, 1.2, 0.7] };
    let mut a = b0.clone();
    a.free_flow_exact(1.3);
    a.free_flow_exact(2.9);
    let mut c = b0.clone();
    c.free_flow_exact(4.2);
    assert!(quat_dist(&a.q, &c.q) < 1e-12);
    a.free_flow_exact(-4.2);
    assert!(quat_dist(&a.q, &b0.q) < 1e-12);
    for k in 0..3 {
        assert!((a.l_body[k] - b0.l_body[k]).abs() < 1e-12);
    }
    let mut d = b0.clone();
    d.free_flow_exact(1000.0);
    assert!((d.kinetic_energy() - b0.kinetic_energy()).abs() < 1e-11);
    let (l0, l1) = (b0.l_space(), d.l_space());
    for k in 0..3 {
        assert!((l0[k] - l1[k]).abs() < 1e-10, "space L must be constant");
    }
}

#[test]
fn exact_flow_reduces_to_symmetric_top_and_steady_rotation() {
    let b0 = RigidBody { inertia: [1.0, 1.0, 2.0], q: [0.5, 0.5, 0.5, 0.5], l_body: [0.3, 0.2, 0.9] };
    let mut a = b0.clone();
    let mut c = b0.clone();
    a.free_flow_exact(3.1);
    c.symmetric_top_flow(2, 3.1);
    assert!(quat_dist(&a.q, &c.q) < 1e-13);
    let nb = RigidBody { inertia: [1.0, 1.0 + 1e-6, 2.0], q: [0.5, 0.5, 0.5, 0.5], l_body: [0.3, 0.2, 0.9] };
    let mut e = nb.clone();
    e.free_flow_exact(3.1);
    assert!(quat_dist(&e.q, &a.q) < 1e-4);
    let s = RigidBody { inertia: [1.0, 2.0, 3.0], q: [1.0, 0.0, 0.0, 0.0], l_body: [0.0, 1.0, 0.0] };
    let mut s1 = s.clone();
    s1.free_flow_exact(2.0);
    let reference = rk4_free_body(&s, 2.0, 2000);
    assert!(quat_dist(&s1.q, &reference.q) < 1e-12);
}

#[test]
fn exact_split_step_with_torque_is_second_order() {
    // potential V = -(k/2) (e_z . R e_0)^2 acting on body axis 0
    let torque = |q: &[f64; 4]| {
        let r = quat_to_matrix(q);
        let a0 = [r[0][0], r[1][0], r[2][0]];
        let c = a0[2];
        let ez = [0.0, 0.0, 1.0];
        let cr = [a0[1] * ez[2] - a0[2] * ez[1], a0[2] * ez[0] - a0[0] * ez[2], a0[0] * ez[1] - a0[1] * ez[0]];
        [0.8 * c * cr[0], 0.8 * c * cr[1], 0.8 * c * cr[2]]
    };
    let b0 = RigidBody { inertia: [1.0, 2.0, 3.0], q: [0.9, 0.1, 0.3, 0.3], l_body: [0.3, 1.0, 0.4] };
    let run = |h: f64| {
        let mut b = b0.clone();
        let n = (2.0 / h).round() as usize;
        for _ in 0..n {
            b.step_exact(h, &torque);
        }
        b
    };
    let reference = run(1.0 / 640.0);
    let e1 = quat_dist(&run(1.0 / 20.0).q, &reference.q);
    let e2 = quat_dist(&run(1.0 / 40.0).q, &reference.q);
    // quat_dist is quadratic in the angle error
    let rate = (e1 / e2).log2() / 2.0;
    assert!(rate > 1.8 && rate < 2.3, "order {rate}");
}

#[test]
fn rattle_step_is_second_order_to_the_exact_free_flow() {
    let b0 = RigidBody { inertia: [1.0, 2.0, 3.0], q: [0.9, 0.1, 0.3, 0.3], l_body: [0.3, 1.0, 0.4] };
    let none = |_: &[f64; 4]| [0.0; 3];
    let mut exact = b0.clone();
    exact.free_flow_exact(3.0);
    let run = |n: usize| {
        let mut b = b0.clone();
        for _ in 0..n {
            b.rattle_step(3.0 / n as f64, &none).unwrap();
        }
        b
    };
    let e1 = quat_dist(&run(30).q, &exact.q).sqrt();
    let e2 = quat_dist(&run(60).q, &exact.q).sqrt();
    let rate = (e1 / e2).log2();
    assert!(rate > 1.8 && rate < 2.3, "order {rate}");
    // exact conservation of the space angular momentum and orthogonality
    let b = run(200);
    let (l0, l1) = (b0.l_space(), b.l_space());
    for k in 0..3 {
        assert!((l0[k] - l1[k]).abs() < 1e-12);
    }
    let nrm: f64 = b.q.iter().map(|x| x * x).sum();
    assert!((nrm - 1.0).abs() < 1e-13);
}

#[test]
fn rattle_step_with_torque_is_second_order_and_reversible() {
    let torque = |q: &[f64; 4]| {
        let r = quat_to_matrix(q);
        let a0 = [r[0][0], r[1][0], r[2][0]];
        let c = a0[2];
        [0.8 * c * a0[1], -0.8 * c * a0[0], 0.0]
    };
    let b0 = RigidBody { inertia: [1.0, 2.0, 3.0], q: [0.9, 0.1, 0.3, 0.3], l_body: [0.3, 1.0, 0.4] };
    let run = |h: f64| {
        let mut b = b0.clone();
        for _ in 0..(2.0 / h).round() as usize {
            b.rattle_step(h, &torque).unwrap();
        }
        b
    };
    let reference = run(1.0 / 400.0);
    let e1 = quat_dist(&run(1.0 / 10.0).q, &reference.q).sqrt();
    let e2 = quat_dist(&run(1.0 / 20.0).q, &reference.q).sqrt();
    let rate = (e1 / e2).log2();
    assert!(rate > 1.8 && rate < 2.3, "order {rate}");
    // reversibility: flip L and run back
    let mut b = run(0.05);
    for v in &mut b.l_body {
        *v = -*v;
    }
    for _ in 0..40 {
        b.rattle_step(0.05, &torque).unwrap();
    }
    assert!(quat_dist(&b.q, &b0.q) < 1e-12);
}
