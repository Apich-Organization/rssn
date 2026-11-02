//! Lattice spin models: exact references, algorithm cross-checks, tempering, Wang-Landau.

use rssn::sim::spin_models::*;

fn series_on<const D: usize, const N: usize>(
    m: &mut OnModel<D, N>,
    u: Update,
    therm: usize,
    n: usize,
) -> (Vec<f64>, Vec<f64>) {
    let sites = m.lattice.sites() as f64;
    for _ in 0..therm {
        m.update(u).unwrap();
    }
    let (mut e, mut mag) = (Vec::new(), Vec::new());
    for _ in 0..n {
        m.update(u).unwrap();
        e.push(m.energy() / sites);
        mag.push(m.abs_magnetisation_per_site());
    }
    (e, mag)
}

fn agree(
    a: &[f64],
    b: &[f64],
    sigmas: f64,
    floor: f64,
) {
    let (ma, ea) = mean_and_error(a);
    let (mb, eb) = mean_and_error(b);
    let tol = sigmas * (ea * ea + eb * eb).sqrt() + floor;
    assert!((ma - mb).abs() < tol, "{ma} vs {mb} (tol {tol})");
}

#[test]
fn lattice_geometry_in_one_two_and_three_dimensions() {
    assert_eq!(Lattice::new([4, 4], true).unwrap().bonds(), 32);
    assert_eq!(Lattice::new([4, 4, 4], true).unwrap().bonds(), 192);
    assert_eq!(Lattice::new([10], false).unwrap().bonds(), 9);
    assert_eq!(Lattice::new([10], true).unwrap().bonds(), 10);
    assert!(Lattice::new([2, 5], true).is_err());
    let lat = Lattice::new([3, 3], false).unwrap();
    assert_eq!(
        lat.neighbours(0)
            .iter()
            .filter(|&&j| j != usize::MAX)
            .count(),
        2
    );
    let mut m = OnModel::<1, 1>::new([8], false, 1.0, 1).unwrap();
    m.cold_start();
    assert!((m.energy() + 7.0).abs() < 1e-12);
}

#[test]
fn one_dimensional_chain_matches_exact_energy_for_metropolis_and_heat_bath() {
    let beta = 0.6;
    let exact = ising1d_energy_per_site(beta, 1.0, 64);
    for (u, seed) in [
        (Update::Metropolis(0.0), 3),
        (Update::HeatBath, 4),
        (Update::Wolff, 5),
        (Update::SwendsenWang, 6),
    ] {
        let mut m = OnModel::<1, 1>::new([64], true, beta, seed).unwrap();
        let (e, _) = series_on(&mut m, u, 500, 6000);
        let (mean, err) = mean_and_error(&e);
        assert!(
            (mean - exact).abs() < 5.0 * err + 0.005,
            "{u:?}: {mean} vs {exact} (err {err})"
        );
    }
}

#[test]
fn chain_formula_matches_enumeration() {
    let dos = enumerate_ising_dos([10], true).unwrap();
    let total: f64 = dos.iter().map(|x| x.1).sum();
    assert_eq!(total, 1024.0);
    let d: Vec<(f64, f64)> = dos.iter().map(|&(e, g)| (e as f64, g.ln())).collect();
    let (_, e, _) = canonical_from_dos(&d, 0.7);
    assert!((e / 10.0 - ising1d_energy_per_site(0.7, 1.0, 10)).abs() < 1e-12);
}

#[test]
fn onsager_references() {
    assert!((onsager_tc(1.0) - 2.269_185_314).abs() < 1e-8);
    assert!((onsager_energy_per_site(onsager_tc(1.0), 1.0) + 2.0_f64.sqrt()).abs() < 1e-6);
    assert!((onsager_magnetisation(2.0, 1.0) - 0.911_319).abs() < 1e-5);
    assert_eq!(onsager_magnetisation(2.5, 1.0), 0.0);
}

#[test]
fn wolff_reproduces_onsager_energy_and_magnetisation_in_2d() {
    let t = 2.0;
    let mut m = OnModel::<2, 1>::new([32, 32], true, 1.0 / t, 11).unwrap();
    let (e, mag) = series_on(&mut m, Update::Wolff, 300, 3000);
    let (me, ee) = mean_and_error(&e);
    let (mm, em) = mean_and_error(&mag);
    assert!(
        (me - onsager_energy_per_site(t, 1.0)).abs() < 5.0 * ee + 0.003,
        "{me}"
    );
    assert!(
        (mm - onsager_magnetisation(t, 1.0)).abs() < 5.0 * em + 0.004,
        "{mm}"
    );
}

#[test]
fn wolff_swendsen_wang_heat_bath_and_metropolis_agree_in_2d() {
    let beta = 1.0 / 2.6;
    let mut a = OnModel::<2, 1>::new([12, 12], true, beta, 21).unwrap();
    let mut b = OnModel::<2, 1>::new([12, 12], true, beta, 22).unwrap();
    let mut c = OnModel::<2, 1>::new([12, 12], true, beta, 23).unwrap();
    let mut d = OnModel::<2, 1>::new([12, 12], true, beta, 24).unwrap();
    let (ea, ma) = series_on(&mut a, Update::Metropolis(0.0), 1000, 20000);
    let (eb, mb) = series_on(&mut b, Update::Wolff, 500, 20000);
    let (ec, mc) = series_on(&mut c, Update::SwendsenWang, 500, 20000);
    let (ed, md) = series_on(&mut d, Update::HeatBath, 1000, 20000);
    agree(&ea, &eb, 5.0, 0.002);
    agree(&ea, &ec, 5.0, 0.002);
    agree(&ea, &ed, 5.0, 0.002);
    agree(&ma, &mb, 5.0, 0.004);
    agree(&ma, &mc, 5.0, 0.004);
    agree(&ma, &md, 5.0, 0.004);
}

#[test]
fn three_dimensional_ising_wolff_matches_metropolis() {
    let beta = 0.2; // paramagnetic side of Tc = 4.5115
    let mut a = OnModel::<3, 1>::new([6, 6, 6], true, beta, 31).unwrap();
    let mut b = OnModel::<3, 1>::new([6, 6, 6], true, beta, 32).unwrap();
    let (ea, _) = series_on(&mut a, Update::Metropolis(0.0), 500, 10000);
    let (eb, _) = series_on(&mut b, Update::Wolff, 300, 10000);
    agree(&ea, &eb, 5.0, 0.002);
}

#[test]
fn binder_crossing_locates_the_2d_critical_temperature() {
    let temps: Vec<f64> = (0..8).map(|k| 2.1 + 0.04 * f64::from(k)).collect();
    let mut curves = Vec::new();
    for (l, seed) in [(8usize, 41u64), (16, 42)] {
        let mut u = Vec::new();
        for &t in &temps {
            let mut m =
                OnModel::<2, 1>::new([l, l], true, 1.0 / t, seed + (t * 1000.0) as u64).unwrap();
            let (_, mag) = series_on(&mut m, Update::Wolff, 400, 6000);
            u.push(binder_cumulant(&mag));
        }
        curves.push(u);
    }
    let tc = binder_crossing(&temps, &curves[0], &curves[1]).unwrap();
    assert!(
        (tc - onsager_tc(1.0)).abs() < 0.08,
        "Binder crossing at {tc}"
    );
}

#[test]
fn xy_metropolis_overrelaxation_wolff_and_sw_agree() {
    let beta = 1.0;
    let mut a = OnModel::<2, 2>::new([8, 8], true, beta, 51).unwrap();
    let mut b = OnModel::<2, 2>::new([8, 8], true, beta, 52).unwrap();
    let mut c = OnModel::<2, 2>::new([8, 8], true, beta, 53).unwrap();
    let sites = 64.0;
    let (mut ea, mut eb, mut ec) = (Vec::new(), Vec::new(), Vec::new());
    for k in 0..40000 {
        a.metropolis_sweep(1.0);
        a.overrelax_sweep();
        b.wolff_step().unwrap();
        c.swendsen_wang_step().unwrap();
        if k > 1000 {
            ea.push(a.energy() / sites);
            eb.push(b.energy() / sites);
            ec.push(c.energy() / sites);
        }
    }
    agree(&ea, &eb, 5.0, 0.004);
    agree(&ea, &ec, 5.0, 0.004);
    let e_before = a.energy();
    a.overrelax_sweep();
    assert!((a.energy() - e_before).abs() < 1e-8);
}

#[test]
fn heisenberg_wolff_matches_metropolis_in_3d_and_1d_exact() {
    // 1-D classical Heisenberg chain: e = -J (coth(K) - 1/K) per bond (thermodynamic limit)
    let beta = 1.5;
    let mut w = OnModel::<1, 3>::new([200], true, beta, 61).unwrap();
    let (e, _) = series_on(&mut w, Update::Wolff, 500, 8000);
    let (mean, err) = mean_and_error(&e);
    let exact = -(1.0 / beta.tanh() - 1.0 / beta);
    assert!((mean - exact).abs() < 5.0 * err + 0.01, "{mean} vs {exact}");
    let mut a = OnModel::<3, 3>::new([4, 4, 4], true, 0.5, 62).unwrap();
    let mut b = OnModel::<3, 3>::new([4, 4, 4], true, 0.5, 63).unwrap();
    let (ea, _) = series_on(&mut a, Update::Metropolis(1.0), 500, 20000);
    let (eb, _) = series_on(&mut b, Update::Wolff, 500, 20000);
    agree(&ea, &eb, 5.0, 0.003);
}

#[test]
fn unsupported_combinations_are_errors() {
    let mut rng0 = rssn::kernels::random::Rng::new(3);
    assert!(heat_bath_direction::<0>(&mut rng0, &[], 1.0).is_err());
    let mut h = OnModel::<2, 1>::new([4, 4], true, 1.0, 1).unwrap();
    h.field = 0.3;
    assert!(h.wolff_step().is_ok());
    h.j = -1.0;
    assert!(h.wolff_step().is_err());
    assert!(Potts::<2>::new([4, 4], true, 1, 1.0, 1).is_err());
}

#[test]
fn two_state_potts_is_the_ising_model() {
    // delta(s,s') = (1 + s s')/2  =>  E_ising = 2 E_potts + N_bonds with beta_potts = 2 beta_ising
    let bi = 0.35;
    let mut p = Potts::<2>::new([8, 8], true, 2, 2.0 * bi, 71).unwrap();
    let mut i = OnModel::<2, 1>::new([8, 8], true, bi, 72).unwrap();
    let nb = 128.0;
    let (mut ep, mut ei) = (Vec::new(), Vec::new());
    for k in 0..20000 {
        p.swendsen_wang_step();
        i.wolff_step().unwrap();
        if k > 500 {
            ep.push((2.0 * p.energy() + nb) / 64.0);
            ei.push(i.energy() / 64.0);
        }
    }
    agree(&ep, &ei, 5.0, 0.003);
}

#[test]
fn potts_updates_agree_and_show_the_q3_transition() {
    let beta = 0.9;
    let run = |u: Update, seed: u64| {
        let mut p = Potts::<2>::new([8, 8], true, 3, beta, seed).unwrap();
        for _ in 0..500 {
            p.update(u);
        }
        let mut e = Vec::new();
        for _ in 0..15000 {
            p.update(u);
            e.push(p.energy() / 64.0);
        }
        e
    };
    let a = run(Update::Metropolis(0.0), 81);
    let b = run(Update::HeatBath, 82);
    let c = run(Update::Wolff, 83);
    let d = run(Update::SwendsenWang, 84);
    agree(&a, &b, 5.0, 0.004);
    agree(&a, &c, 5.0, 0.004);
    agree(&a, &d, 5.0, 0.004);
    // beta_c = ln(1 + sqrt 3)
    let bc = (1.0 + 3.0_f64.sqrt()).ln();
    let order = |beta: f64| {
        let mut p = Potts::<2>::new([16, 16], true, 3, beta, 85).unwrap();
        let mut acc = 0.0;
        for k in 0..2000 {
            p.swendsen_wang_step();
            if k >= 500 {
                acc += p.order_parameter();
            }
        }
        acc / 1500.0
    };
    assert!(order(bc * 1.3) > 0.85);
    assert!(order(bc * 0.7) < 0.25);
}

#[test]
fn wang_landau_density_of_states_matches_enumeration() {
    let exact = enumerate_ising_dos([4, 4], true).unwrap();
    assert_eq!(exact.first().unwrap().0, -32);
    assert_eq!(exact.first().unwrap().1, 2.0);
    let wl = wang_landau_ising([4, 4], true, 1e-5, 0.8, 91).unwrap();
    assert_eq!(wl.len(), exact.len());
    let dos_exact: Vec<(f64, f64)> = exact.iter().map(|&(e, g)| (e as f64, g.ln())).collect();
    let dos_wl: Vec<(f64, f64)> = wl.iter().map(|&(e, g)| (e as f64, g)).collect();
    for beta in [0.2, 0.44, 0.7] {
        let (_, e1, c1) = canonical_from_dos(&dos_exact, beta);
        let (_, e2, c2) = canonical_from_dos(&dos_wl, beta);
        assert!(
            (e1 - e2).abs() < 0.03 * e1.abs(),
            "beta {beta}: <E> {e1} vs {e2}"
        );
        assert!((c1 - c2).abs() < 0.08 * c1, "beta {beta}: C {c1} vs {c2}");
    }
}

#[test]
fn histogram_reweighting_extrapolates_energy() {
    let exact = enumerate_ising_dos([4, 4], true).unwrap();
    let dos: Vec<(f64, f64)> = exact.iter().map(|&(e, g)| (e as f64, g.ln())).collect();
    let b0 = 0.4;
    let mut m = OnModel::<2, 1>::new([4, 4], true, b0, 95).unwrap();
    for _ in 0..500 {
        m.wolff_step().unwrap();
    }
    let mut en = Vec::new();
    for _ in 0..60000 {
        m.metropolis_sweep(0.0);
        m.wolff_step().unwrap();
        en.push(m.energy());
    }
    for b1 in [0.36, 0.40, 0.44] {
        let rw = reweight(&en, &en, b0, b1);
        let (_, ex, _) = canonical_from_dos(&dos, b1);
        assert!(
            (rw - ex).abs() < 0.01 * ex.abs(),
            "beta {b1}: reweighted {rw} vs exact {ex}"
        );
    }
}

#[test]
fn parallel_tempering_swaps_and_samples_each_temperature() {
    let betas = [0.3, 0.37, 0.44, 0.5];
    let mut pt = ParallelTempering::<2, 1>::new([8, 8], true, &betas, 101).unwrap();
    let mut slot2 = Vec::new();
    pt.run(6000, Update::Metropolis(0.0), |round, reps| {
        if round > 500 {
            slot2.push(reps[2].energy() / 64.0);
        }
    })
    .unwrap();
    for (att, acc) in &pt.swap_stats {
        assert!(
            *att > 0 && (*acc as f64) / (*att as f64) > 0.1,
            "swap {acc}/{att}"
        );
    }
    let mut w = OnModel::<2, 1>::new([8, 8], true, 0.44, 102).unwrap();
    let (e, _) = series_on(&mut w, Update::Wolff, 500, 20000);
    agree(&slot2, &e, 5.0, 0.01);
}

#[test]
fn scaling_helpers() {
    let l = [4.0, 8.0, 16.0, 32.0];
    let y: Vec<f64> = l.iter().map(|v: &f64| v.powf(1.75)).collect();
    assert!((fit_exponent(&l, &y) - 1.75).abs() < 1e-12);
    assert!((collapse_variable(2.4, 2.269, 1.0, 16.0) - 0.131 * 16.0).abs() < 1e-12);
    assert!((binder_cumulant(&[1.0, -1.0, 1.0, -1.0]) - 2.0 / 3.0).abs() < 1e-12);
    let t = [1.0, 2.0, 3.0];
    assert!(
        (binder_crossing(&t, &[0.0, 1.0, 2.0], &[2.0, 1.5, 1.0]).unwrap() - 2.0 - 1.0 / 3.0).abs()
            < 1e-12
    );
    assert!(binder_crossing(&t, &[0.0, 0.0, 0.0], &[1.0, 1.0, 1.0]).is_none());
}

fn bessel_i(
    nu: usize,
    x: f64,
) -> f64 {
    // power series of the modified Bessel function of integer order
    let mut term = (x / 2.0).powi(nu as i32) / (1..=nu).map(|k| k as f64).product::<f64>();
    let mut sum = term;
    for k in 1..60 {
        term *= (x / 2.0) * (x / 2.0) / (k as f64 * (k + nu) as f64);
        sum += term;
    }
    sum
}

#[test]
fn heat_bath_direction_has_the_exact_conditional_mean() {
    use rssn::kernels::random::Rng;
    let mut rng = Rng::new(7);
    let n = 200_000;
    for a in [0.05_f64, 0.7, 2.5, 9.0] {
        // XY: <cos theta> = I1(a)/I0(a)
        let f: [f64; 2] = [0.3 * a / 0.5, 0.4 * a / 0.5];
        let beta = 0.5;
        let fm = (f[0] * f[0] + f[1] * f[1]).sqrt();
        let kappa = beta * fm;
        let mut acc = 0.0;
        for _ in 0..n {
            let s = heat_bath_direction::<2>(&mut rng, &f, beta).unwrap();
            acc += (s[0] * f[0] + s[1] * f[1]) / fm;
        }
        let exact = bessel_i(1, kappa) / bessel_i(0, kappa);
        assert!(
            (acc / n as f64 - exact).abs() < 5e-3,
            "XY kappa {kappa}: {} vs {exact}",
            acc / n as f64
        );
        // Heisenberg: <cos theta> = coth(a) - 1/a, and the perpendicular components are isotropic
        let f3 = [f[0], f[1], 0.2 * a];
        let f3m = (f3[0] * f3[0] + f3[1] * f3[1] + f3[2] * f3[2]).sqrt();
        let k3 = beta * f3m;
        let (mut c, mut norm_err) = (0.0, 0.0_f64);
        for _ in 0..n {
            let s = heat_bath_direction::<3>(&mut rng, &f3, beta).unwrap();
            c += (s[0] * f3[0] + s[1] * f3[1] + s[2] * f3[2]) / f3m;
            norm_err = norm_err.max((s[0] * s[0] + s[1] * s[1] + s[2] * s[2] - 1.0).abs());
        }
        let exact3 = 1.0 / k3.tanh() - 1.0 / k3;
        assert!(
            (c / n as f64 - exact3).abs() < 5e-3,
            "Heisenberg a {k3}: {} vs {exact3}",
            c / n as f64
        );
        assert!(norm_err < 1e-12);
    }
}

/// `<cos theta> = I_{N/2}(k) / I_{N/2-1}(k)` for the von Mises-Fisher distribution on `S^(N-1)`.
fn vmf_mean_cosine(
    n: usize,
    k: f64,
) -> f64 {
    // modified Bessel functions of half-integer order via the power series
    let iv = |nu: f64| {
        let mut term = (0.5 * k).powf(nu) / rssn::kernels::special::gamma_numerical(nu + 1.0);
        let mut sum = term;
        for m in 1..200 {
            term *= 0.25 * k * k / (m as f64 * (m as f64 + nu));
            sum += term;
        }
        sum
    };
    iv(0.5 * n as f64) / iv(0.5 * n as f64 - 1.0)
}

#[test]
fn heat_bath_for_higher_n_is_von_mises_fisher_and_exact_on_chains() {
    use rssn::kernels::random::Rng;
    let mut rng = Rng::new(17);
    let samples = 150_000;
    for (a, beta) in [(0.0_f64, 1.0), (0.5, 1.0), (3.0, 1.0), (15.0, 0.5)] {
        macro_rules! check_n {
            ($n:literal) => {{
                let mut f = [0.0; $n];
                f[0] = 0.6 * a / beta;
                f[1] = 0.8 * a / beta;
                f[$n - 1] += 0.1 * a / beta;
                let fm = f.iter().map(|x| x * x).sum::<f64>().sqrt();
                let (mut c, mut c2, mut worst) = (0.0, 0.0, 0.0_f64);
                for _ in 0..samples {
                    let s = heat_bath_direction::<$n>(&mut rng, &f, beta).unwrap();
                    let cs = if fm > 0.0 {
                        s.iter().zip(&f).map(|(x, y)| x * y).sum::<f64>() / fm
                    } else {
                        s[0]
                    };
                    c += cs;
                    c2 += cs * cs;
                    worst = worst.max((s.iter().map(|x| x * x).sum::<f64>() - 1.0).abs());
                }
                let kappa = beta * fm;
                let exact = if kappa < 1e-12 {
                    0.0
                } else {
                    vmf_mean_cosine($n, kappa)
                };
                assert!(
                    (c / samples as f64 - exact).abs() < 6e-3,
                    "N={} kappa {kappa}: {} vs {exact}",
                    $n,
                    c / samples as f64
                );
                if a == 0.0 {
                    // uniform sphere: <cos^2> = 1/N
                    assert!((c2 / samples as f64 - 1.0 / $n as f64).abs() < 5e-3);
                }
                assert!(worst < 1e-12);
            }};
        }
        check_n!(4);
        check_n!(5);
        check_n!(6);
    }
    // open O(4) chain: <s_i . s_{i+1}> = A_4(beta) for every bond exactly
    let beta = 1.3;
    let mut chain = OnModel::<1, 4>::new([60], false, beta, 41).unwrap();
    let (e, _) = series_on(&mut chain, Update::HeatBath, 300, 6000);
    let (mean, err) = mean_and_error(&e);
    let exact = -(59.0 / 60.0) * vmf_mean_cosine(4, beta);
    assert!(
        (mean - exact).abs() < 5.0 * err + 0.005,
        "{mean} vs {exact}"
    );
    // O(4) 2-D: heat bath and Wolff cluster updates sample the same distribution
    let mut hb = OnModel::<2, 4>::new([6, 6], true, 0.8, 71).unwrap();
    let mut wf = OnModel::<2, 4>::new([6, 6], true, 0.8, 72).unwrap();
    let (e1, _) = series_on(&mut hb, Update::HeatBath, 500, 15000);
    let (e2, _) = series_on(&mut wf, Update::Wolff, 500, 15000);
    agree(&e1, &e2, 5.0, 0.004);
}

#[test]
fn xy_and_heisenberg_heat_bath_agree_with_metropolis_and_exact_chains() {
    // XY 2-D
    let mut a = OnModel::<2, 2>::new([8, 8], true, 1.0, 201).unwrap();
    let mut b = OnModel::<2, 2>::new([8, 8], true, 1.0, 202).unwrap();
    let (ea, ma) = series_on(&mut a, Update::HeatBath, 1000, 20000);
    let (eb, mb) = series_on(&mut b, Update::Metropolis(1.0), 1000, 40000);
    agree(&ea, &eb, 5.0, 0.004);
    agree(&ma, &mb, 5.0, 0.01);
    // Heisenberg 1-D ring: e = -(coth K - 1/K) per bond in the thermodynamic limit
    let beta = 1.5;
    let mut h = OnModel::<1, 3>::new([200], true, beta, 203).unwrap();
    let (e, _) = series_on(&mut h, Update::HeatBath, 500, 6000);
    let (mean, err) = mean_and_error(&e);
    let exact = -(1.0 / beta.tanh() - 1.0 / beta);
    assert!((mean - exact).abs() < 5.0 * err + 0.01, "{mean} vs {exact}");
    // XY 1-D ring: e = -I1(K)/I0(K) per bond
    let mut x = OnModel::<1, 2>::new([200], true, beta, 204).unwrap();
    let (e, _) = series_on(&mut x, Update::HeatBath, 500, 6000);
    let (mean, err) = mean_and_error(&e);
    let exact = -bessel_i(1, beta) / bessel_i(0, beta);
    assert!((mean - exact).abs() < 5.0 * err + 0.01, "{mean} vs {exact}");
    // heat bath with a field is allowed (it is just part of the local field)
    let mut f = OnModel::<2, 3>::new([6, 6], true, 1.0, 205).unwrap();
    f.field = 0.5;
    for _ in 0..200 {
        f.heat_bath_sweep().unwrap();
    }
    assert!(f.magnetisation()[0] > 0.0);
}

#[test]
fn potts_parallel_tempering_matches_independent_runs() {
    let betas = [0.6, 0.75, 0.9, 1.05];
    let mut pt = PottsTempering::<2>::new([8, 8], true, 3, &betas, 301).unwrap();
    let mut slot1 = Vec::new();
    pt.run(5000, Update::Metropolis(0.0), |round, reps| {
        if round > 500 {
            slot1.push(reps[1].energy() / 64.0);
        }
    });
    for (att, acc) in &pt.swap_stats {
        assert!(
            *att > 0 && (*acc as f64) / (*att as f64) > 0.1,
            "swap {acc}/{att}"
        );
    }
    let mut w = Potts::<2>::new([8, 8], true, 3, 0.75, 302).unwrap();
    let mut e = Vec::new();
    for k in 0..30000 {
        w.update(Update::SwendsenWang);
        if k > 500 {
            e.push(w.energy() / 64.0);
        }
    }
    agree(&slot1, &e, 5.0, 0.01);
    assert!(PottsTempering::<2>::new([4, 4], true, 3, &[1.0], 1).is_err());
}

#[test]
fn wang_landau_for_continuous_spins_reproduces_open_chain_thermodynamics() {
    // open chain of 6 spins: Z ~ (4 pi sinh K / K)^(N-1) (Heisenberg), (2 pi)^N I0(K)^(N-1) (XY)
    let nb = 5.0;
    let dos_h = wang_landau_on::<1, 3>(
        [6],
        false,
        (-nb - 1e-9, nb + 1e-9, 100),
        1.0,
        1e-4,
        0.7,
        401,
    )
    .unwrap();
    let dos_x = wang_landau_on::<1, 2>(
        [6],
        false,
        (-nb - 1e-9, nb + 1e-9, 100),
        1.0,
        1e-4,
        0.7,
        402,
    )
    .unwrap();
    for beta in [0.3, 0.8, 1.4] {
        let (_, e_h, _) = canonical_from_dos(&dos_h, beta);
        let exact_h = -nb * (1.0 / beta.tanh() - 1.0 / beta);
        assert!(
            (e_h - exact_h).abs() < 0.05 * exact_h.abs(),
            "Heisenberg beta {beta}: {e_h} vs {exact_h}"
        );
        let (_, e_x, _) = canonical_from_dos(&dos_x, beta);
        let exact_x = -nb * bessel_i(1, beta) / bessel_i(0, beta);
        assert!(
            (e_x - exact_x).abs() < 0.05 * exact_x.abs(),
            "XY beta {beta}: {e_x} vs {exact_x}"
        );
    }
    assert!(wang_landau_on::<1, 3>([6], false, (1.0, 0.0, 10), 1.0, 1e-3, 0.7, 1).is_err());
}

#[test]
fn replica_exchange_wang_landau_reproduces_the_exact_ising_density_of_states() {
    // 4x4 periodic Ising: levels -32 + 4k; bins of width 4 centred on the levels
    let exact = enumerate_ising_dos::<2>([4, 4], true).unwrap();
    let total: f64 = exact.iter().map(|&(_, c)| c).sum();
    let opt = RewlOptions {
        range: (-34.0, 34.0),
        nbins: 17,
        windows: 3,
        overlap: 0.5,
        walkers: 2,
        delta: 1.0,
        ln_f_final: 1e-5,
        flatness: 0.7,
        sweeps_per_exchange: 50,
        seed: 11,
    };
    let dos = replica_exchange_wang_landau::<2, 1>([4, 4], true, &opt).unwrap();
    assert_eq!(dos.len(), exact.len());
    for ((e, lg), (ee, c)) in dos.iter().zip(&exact) {
        assert!((e - *ee as f64).abs() < 1e-9, "level {e} vs {ee}");
        let ln_exact = (c / total).ln();
        assert!(
            (lg - ln_exact).abs() < 0.15,
            "E={e}: ln g {lg} vs {ln_exact}"
        );
    }
    // thermodynamics from the joined density of states
    let (_, e_wl, _) = canonical_from_dos(&dos, 0.4);
    let (_, e_ex, _) = canonical_from_dos(
        &exact
            .iter()
            .map(|&(e, c)| (e as f64, (c / total).ln()))
            .collect::<Vec<_>>(),
        0.4,
    );
    assert!((e_wl - e_ex).abs() < 0.05 * e_ex.abs());
}

#[test]
fn replica_exchange_wang_landau_for_continuous_spins_matches_open_chain() {
    let nb = 5.0;
    let opt = RewlOptions {
        range: (-nb - 1e-9, nb + 1e-9),
        nbins: 60,
        windows: 3,
        overlap: 0.6,
        walkers: 2,
        delta: 1.0,
        ln_f_final: 1e-4,
        flatness: 0.7,
        sweeps_per_exchange: 100,
        seed: 5,
    };
    let dos = replica_exchange_wang_landau::<1, 3>([6], false, &opt).unwrap();
    for beta in [0.3, 0.8, 1.4] {
        let (_, e_h, _) = canonical_from_dos(&dos, beta);
        let exact_h = -nb * (1.0 / beta.tanh() - 1.0 / beta);
        assert!(
            (e_h - exact_h).abs() < 0.06 * exact_h.abs(),
            "beta {beta}: {e_h} vs {exact_h}"
        );
    }
    let mut bad = opt;
    bad.walkers = 0;
    assert!(replica_exchange_wang_landau::<1, 3>([6], false, &bad).is_err());
}

#[test]
fn wolff_and_swendsen_wang_with_a_field_reproduce_exact_ising_averages() {
    // 4x4 periodic Ising in a field: exact canonical averages by enumeration
    let (beta, h) = (0.45, 0.3);
    let dims = [4usize, 4];
    let n = 16;
    let bonds: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| {
            [
                (i, (i / 4) * 4 + (i % 4 + 1) % 4),
                (i, ((i / 4 + 1) % 4) * 4 + i % 4),
            ]
        })
        .collect();
    let (mut z, mut e_acc, mut m_acc) = (0.0, 0.0, 0.0);
    for state in 0u32..(1 << n) {
        let s: Vec<f64> = (0..n)
            .map(|i| {
                if (state >> i) & 1 == 1 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let mag: f64 = s.iter().sum();
        let e = -bonds.iter().map(|&(i, j)| s[i] * s[j]).sum::<f64>() - h * mag;
        let w = (-beta * e).exp();
        z += w;
        e_acc += w * e;
        m_acc += w * mag;
    }
    let (e_ex, m_ex) = (e_acc / z / n as f64, m_acc / z / n as f64);
    for upd in [Update::Wolff, Update::SwendsenWang] {
        let mut model = OnModel::<2, 1>::new(dims, true, beta, 77).unwrap();
        model.field = h;
        let (mut es, mut ms) = (Vec::new(), Vec::new());
        for k in 0..60_000 {
            model.update(upd).unwrap();
            if k > 2000 {
                es.push(model.energy() / n as f64);
                ms.push(model.magnetisation()[0] / n as f64);
            }
        }
        let (e, ee) = mean_and_error(&es);
        let (m, me) = mean_and_error(&ms);
        assert!(
            (e - e_ex).abs() < 5.0 * ee + 2e-3,
            "{upd:?}: e = {e} +- {ee} vs {e_ex}"
        );
        assert!(
            (m - m_ex).abs() < 5.0 * me + 2e-3,
            "{upd:?}: m = {m} +- {me} vs {m_ex}"
        );
    }
}

#[test]
fn cluster_updates_with_a_field_agree_with_heat_bath_for_o_n_models() {
    fn check<const N: usize>(seed: u64) {
        let mut a = OnModel::<2, N>::new([5, 5], true, 0.6, seed).unwrap();
        let mut b = OnModel::<2, N>::new([5, 5], true, 0.6, seed + 1).unwrap();
        let mut c = OnModel::<2, N>::new([5, 5], true, 0.6, seed + 2).unwrap();
        for m in [&mut a, &mut b, &mut c] {
            m.field = -0.4;
        }
        let (ea, ma) = series_on(&mut a, Update::HeatBath, 1000, 30_000);
        let (eb, mb) = series_on(&mut b, Update::Wolff, 1000, 30_000);
        let (ec, mc) = series_on(&mut c, Update::SwendsenWang, 1000, 30_000);
        agree(&ea, &eb, 5.0, 0.01);
        agree(&ma, &mb, 5.0, 0.01);
        agree(&ea, &ec, 5.0, 0.01);
        agree(&ma, &mc, 5.0, 0.01);
    }
    check::<2>(500);
    check::<3>(510);
}

// ---------------------------------------------------------------------------
// Next-nearest neighbours (J2)
// ---------------------------------------------------------------------------

/// One bond as a pair of site indices.
type Bond = (usize, usize);

/// Nearest-neighbour and diagonal bond lists.
type Bonds = (Vec<Bond>, Vec<Bond>);

/// Brute-force bond lists `(nearest, diagonal)` of a hypercubic lattice, derived from the site
/// coordinates alone: two sites are a diagonal pair when their minimum-image coordinate
/// differences vanish everywhere except on exactly two axes, where they are one.
fn brute_bonds<const D: usize>(
    dims: [usize; D],
    periodic: bool,
) -> Bonds {
    let n: usize = dims.iter().product();
    let mut stride = [1usize; D];
    for k in 1..D {
        stride[k] = stride[k - 1] * dims[k - 1];
    }
    let coord = |i: usize| -> Vec<usize> { (0..D).map(|k| (i / stride[k]) % dims[k]).collect() };
    let (mut nn, mut diag) = (Vec::new(), Vec::new());
    for i in 0..n {
        let ci = coord(i);
        for j in (i + 1)..n {
            let cj = coord(j);
            let (mut unit, mut rest) = (0usize, 0usize);
            for k in 0..D {
                let d = ci[k].abs_diff(cj[k]);
                let m = if periodic {
                    d.min(dims[k] - d)
                } else {
                    d
                };
                if m == 1 {
                    unit += 1;
                } else if m != 0 {
                    rest += 1;
                }
            }
            if rest == 0 && unit == 1 {
                nn.push((i, j));
            } else if rest == 0 && unit == 2 {
                diag.push((i, j));
            }
        }
    }
    (nn, diag)
}

/// Energy of an Ising configuration from explicit bond lists.
fn ising_energy_of(
    spins: &[f64],
    bonds: &Bonds,
    (j, j2, h): (f64, f64, f64),
) -> f64 {
    let mut e = 0.0;
    for &(i, k) in &bonds.0 {
        e -= j * spins[i] * spins[k];
    }
    for &(i, k) in &bonds.1 {
        e -= j2 * spins[i] * spins[k];
    }
    for &s in spins {
        e -= h * s;
    }
    e
}

fn check_diagonal_table<const D: usize>(
    dims: [usize; D],
    periodic: bool,
) -> usize {
    let lat = Lattice::new(dims, periodic).unwrap();
    let (nn, diag) = brute_bonds(dims, periodic);
    assert_eq!(lat.bonds(), nn.len(), "{dims:?}");
    assert_eq!(lat.diagonal_bonds(), diag.len(), "{dims:?}");
    let mut seen = std::collections::BTreeSet::new();
    for i in 0..lat.sites() {
        for &j in lat.second_neighbours(i) {
            if j == usize::MAX {
                continue;
            }
            assert_ne!(i, j);
            assert!(
                lat.second_neighbours(j).contains(&i),
                "{dims:?}: {i} -> {j} is not mutual"
            );
            seen.insert((i.min(j), i.max(j)));
        }
    }
    assert_eq!(
        seen.iter().copied().collect::<Vec<_>>(),
        diag,
        "{dims:?}: diagonal table differs from the brute-force bond list"
    );
    diag.len()
}

#[test]
fn diagonal_bond_table_matches_a_coordinate_brute_force() {
    assert_eq!(check_diagonal_table([4, 4], true), 32);
    assert_eq!(check_diagonal_table([4, 4], false), 18);
    assert_eq!(check_diagonal_table([3, 5], true), 30);
    // odd side lengths with periodic wrapping: 4 diagonals per site, 18 bonds
    assert_eq!(check_diagonal_table([3, 3], true), 18);
    assert_eq!(check_diagonal_table([4, 4, 4], true), 384);
    assert_eq!(check_diagonal_table([3, 3, 3], true), 162);
    // open 3x3x3: 3 axis pairs x 2 (3-1)^2 diagonals per plane x 3 planes
    assert_eq!(check_diagonal_table([3, 3, 3], false), 72);
    assert_eq!(check_diagonal_table([10], true), 0);
    assert_eq!(
        Lattice::new([6, 6], true)
            .unwrap()
            .second_neighbours(0)
            .len(),
        4
    );
    assert_eq!(
        Lattice::new([6, 6, 6], true)
            .unwrap()
            .second_neighbours(0)
            .len(),
        12
    );
    // an open 1-wide strip has no diagonal neighbours at all
    assert!(
        Lattice::new([1, 5], false)
            .unwrap()
            .second_neighbours(0)
            .iter()
            .all(|&j| j == usize::MAX)
    );
}

#[test]
fn next_nearest_neighbour_energy_matches_an_explicit_bond_sum() {
    let dims = [5usize, 6];
    let bonds = brute_bonds(dims, true);
    let mut m = OnModel::<2, 1>::new(dims, true, 0.7, 901).unwrap();
    m.j = 0.9;
    m.j2 = -0.43;
    m.field = 0.21;
    let n = m.lattice.sites();
    let spins: Vec<f64> = m.spins.iter().map(|s| s[0]).collect();
    let explicit = ising_energy_of(&spins, &bonds, (m.j, m.j2, m.field));
    assert!(
        (m.energy() - explicit).abs() < 1e-10,
        "energy {} vs explicit {explicit}",
        m.energy()
    );
    // the local field is the incremental path: -1/2 sum_i s_i . (f_i - h e_0) - h sum_i s_i
    let mut incremental = 0.0;
    for i in 0..n {
        let f = m.local_field_public(i);
        incremental -= 0.5 * (m.spins[i][0] * f[0] - m.field * m.spins[i][0]);
    }
    incremental -= m.field * spins.iter().sum::<f64>();
    assert!(
        (incremental - explicit).abs() < 1e-10,
        "local field {incremental} vs explicit {explicit}"
    );
    // j2 = 0 reproduces the nearest-neighbour energy exactly
    m.j2 = 0.0;
    let nn_only = ising_energy_of(&spins, &bonds, (m.j, 0.0, m.field));
    assert!((m.energy() - nn_only).abs() < 1e-10);
    // and so does a three-dimensional lattice with open boundaries
    let dims3 = [4usize, 3, 5];
    let bonds3 = brute_bonds(dims3, false);
    let mut m3 = OnModel::<3, 1>::new(dims3, false, 0.4, 902).unwrap();
    m3.j = -1.3;
    m3.j2 = 0.77;
    m3.field = -0.11;
    let spins3: Vec<f64> = m3.spins.iter().map(|s| s[0]).collect();
    let explicit3 = ising_energy_of(&spins3, &bonds3, (m3.j, m3.j2, m3.field));
    assert!((m3.energy() - explicit3).abs() < 1e-10);
    assert!(m3.lattice.diagonal_bonds() > 0);
}

#[test]
fn cluster_updates_with_next_nearest_neighbours_reproduce_exact_averages() {
    // 4x4 periodic Ising with a frustrated diagonal coupling and a field, exact by enumeration
    let (beta, j, j2, h) = (0.5, 1.0, -0.6, 0.25);
    let dims = [4usize, 4];
    let n = 16;
    let bonds = brute_bonds(dims, true);
    let (mut z, mut e_acc, mut m_acc) = (0.0, 0.0, 0.0);
    for state in 0u32..(1 << n) {
        let spins: Vec<f64> = (0..n)
            .map(|i| {
                if (state >> i) & 1 == 1 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let e = ising_energy_of(&spins, &bonds, (j, j2, h));
        let w = (-beta * e).exp();
        z += w;
        e_acc += w * e;
        m_acc += w * spins.iter().sum::<f64>();
    }
    let (e_ex, m_ex) = (e_acc / z / n as f64, m_acc / z / n as f64);
    for upd in [Update::Wolff, Update::SwendsenWang, Update::HeatBath] {
        let mut model = OnModel::<2, 1>::new(dims, true, beta, 903).unwrap();
        model.j = j;
        model.j2 = j2;
        model.field = h;
        let (mut es, mut ms) = (Vec::new(), Vec::new());
        for k in 0..100_000 {
            model.update(upd).unwrap();
            if k > 2_000 {
                es.push(model.energy() / n as f64);
                ms.push(model.magnetisation()[0] / n as f64);
            }
        }
        let (e, ee) = mean_and_error(&es);
        let (m, me) = mean_and_error(&ms);
        assert!(
            (e - e_ex).abs() < 5.0 * ee + 0.02,
            "{upd:?}: e = {e} +- {ee} vs {e_ex}"
        );
        assert!(
            (m - m_ex).abs() < 5.0 * me + 0.02,
            "{upd:?}: m = {m} +- {me} vs {m_ex}"
        );
    }
}

#[test]
fn xy_cluster_updates_with_next_nearest_neighbours_match_heat_bath() {
    for j2 in [-0.5, 0.5] {
        let mut a = OnModel::<2, 2>::new([6, 6], true, 1.5, 904).unwrap();
        let mut b = OnModel::<2, 2>::new([6, 6], true, 1.5, 905).unwrap();
        let mut c = OnModel::<2, 2>::new([6, 6], true, 1.5, 906).unwrap();
        for m in [&mut a, &mut b, &mut c] {
            m.j2 = j2;
        }
        let (ea, _) = series_on(&mut a, Update::HeatBath, 500, 20_000);
        let (eb, _) = series_on(&mut b, Update::Wolff, 500, 20_000);
        let (ec, _) = series_on(&mut c, Update::SwendsenWang, 500, 20_000);
        agree(&ea, &eb, 5.0, 0.01);
        agree(&ea, &ec, 5.0, 0.01);
    }
}

#[test]
fn square_lattice_is_neel_below_and_striped_above_the_j2_critical_ratio() {
    // `H = -J sum_nn - J2 sum_diag`, so the literature's antiferromagnetic `J1 = J2 = 1` model is
    // `j = -1`, `j2 = -ratio`, and `j2 / j` is the literature's `J2 / J1` (critical near 0.5).
    //
    // Order parameters: `neel = |sum (-1)^(x+y) s| / N` and `stripe = |sum (-1)^x s| / N +
    // |sum (-1)^y s| / N`, which is 1 for any stripe state (the two orientations add up) and 0 for
    // both the Neel and the ferromagnetic state.
    let run = |ratio: f64, seed: u64| -> (f64, f64, f64) {
        let mut m = OnModel::<2, 1>::new([12, 12], true, 0.5, seed).unwrap();
        m.j = -1.0;
        m.j2 = -ratio;
        let n = m.lattice.sites();
        let lx = m.lattice.dims[0];
        for step in 0..60 {
            m.beta = 0.5 + 0.1 * f64::from(step);
            for _ in 0..30 {
                m.heat_bath_sweep().unwrap();
            }
        }
        m.beta = 6.5;
        for _ in 0..600 {
            m.heat_bath_sweep().unwrap();
        }
        let samples = 300;
        let (mut neel, mut stripe, mut e) = (0.0, 0.0, 0.0);
        for _ in 0..samples {
            m.heat_bath_sweep().unwrap();
            let (mut mn, mut mx, mut my) = (0.0, 0.0, 0.0);
            for i in 0..n {
                let x = i % lx;
                let y = i / lx;
                let s = m.spins[i][0];
                mn += if (x + y).is_multiple_of(2) {
                    s
                } else {
                    -s
                };
                mx += if x.is_multiple_of(2) { s } else { -s };
                my += if y.is_multiple_of(2) { s } else { -s };
            }
            let inv = 1.0 / n as f64;
            neel += (mn * inv).abs();
            stripe += (mx * inv).abs() + (my * inv).abs();
            e += m.energy() * inv;
        }
        (
            neel / f64::from(samples),
            stripe / f64::from(samples),
            e / f64::from(samples),
        )
    };
    let (neel_low, stripe_low, e_low) = run(0.2, 911);
    let (neel_high, stripe_high, e_high) = run(0.8, 912);
    assert!(neel_low > 0.8, "Neel order at J2/J1 = 0.2 is {neel_low}");
    assert!(
        stripe_low < 0.2,
        "stripe order at J2/J1 = 0.2 is {stripe_low}"
    );
    assert!(
        stripe_high > 0.8,
        "stripe order at J2/J1 = 0.8 is {stripe_high}"
    );
    assert!(neel_high < 0.2, "Neel order at J2/J1 = 0.8 is {neel_high}");
    assert!(
        neel_low > neel_high + 0.5 && stripe_high > stripe_low + 0.5,
        "the two order parameters do not separate: Neel {neel_low}/{neel_high}, \
         stripe {stripe_low}/{stripe_high}"
    );
    // ground-state energies per site are -1.6 in both phases
    assert!(e_low < -1.4, "energy at J2/J1 = 0.2 is {e_low}");
    assert!(e_high < -1.4, "energy at J2/J1 = 0.8 is {e_high}");
}
