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
    assert_eq!(lat.neighbours(0).iter().filter(|&&j| j != usize::MAX).count(), 2);
    let mut m = OnModel::<1, 1>::new([8], false, 1.0, 1).unwrap();
    m.cold_start();
    assert!((m.energy() + 7.0).abs() < 1e-12);
}

#[test]
fn one_dimensional_chain_matches_exact_energy_for_metropolis_and_heat_bath() {
    let beta = 0.6;
    let exact = ising1d_energy_per_site(beta, 1.0, 64);
    for (u, seed) in [(Update::Metropolis(0.0), 3), (Update::HeatBath, 4), (Update::Wolff, 5), (Update::SwendsenWang, 6)] {
        let mut m = OnModel::<1, 1>::new([64], true, beta, seed).unwrap();
        let (e, _) = series_on(&mut m, u, 500, 6000);
        let (mean, err) = mean_and_error(&e);
        assert!((mean - exact).abs() < 5.0 * err + 0.005, "{u:?}: {mean} vs {exact} (err {err})");
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
    assert!((me - onsager_energy_per_site(t, 1.0)).abs() < 5.0 * ee + 0.003, "{me}");
    assert!((mm - onsager_magnetisation(t, 1.0)).abs() < 5.0 * em + 0.004, "{mm}");
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
            let mut m = OnModel::<2, 1>::new([l, l], true, 1.0 / t, seed + (t * 1000.0) as u64).unwrap();
            let (_, mag) = series_on(&mut m, Update::Wolff, 400, 6000);
            u.push(binder_cumulant(&mag));
        }
        curves.push(u);
    }
    let tc = binder_crossing(&temps, &curves[0], &curves[1]).unwrap();
    assert!((tc - onsager_tc(1.0)).abs() < 0.08, "Binder crossing at {tc}");
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
    assert!(h.wolff_step().is_err());
    assert!(h.swendsen_wang_step().is_err());
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
        assert!((e1 - e2).abs() < 0.03 * e1.abs(), "beta {beta}: <E> {e1} vs {e2}");
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
        assert!((rw - ex).abs() < 0.01 * ex.abs(), "beta {b1}: reweighted {rw} vs exact {ex}");
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
        assert!(*att > 0 && (*acc as f64) / (*att as f64) > 0.1, "swap {acc}/{att}");
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
    assert!((binder_crossing(&t, &[0.0, 1.0, 2.0], &[2.0, 1.5, 1.0]).unwrap() - 2.0 - 1.0 / 3.0).abs() < 1e-12);
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
        assert!((acc / n as f64 - exact).abs() < 5e-3, "XY kappa {kappa}: {} vs {exact}", acc / n as f64);
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
        assert!((c / n as f64 - exact3).abs() < 5e-3, "Heisenberg a {k3}: {} vs {exact3}", c / n as f64);
        assert!(norm_err < 1e-12);
    }
}

/// `<cos theta> = I_{N/2}(k) / I_{N/2-1}(k)` for the von Mises-Fisher distribution on `S^(N-1)`.
fn vmf_mean_cosine(n: usize, k: f64) -> f64 {
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
                    let cs = if fm > 0.0 { s.iter().zip(&f).map(|(x, y)| x * y).sum::<f64>() / fm } else { s[0] };
                    c += cs;
                    c2 += cs * cs;
                    worst = worst.max((s.iter().map(|x| x * x).sum::<f64>() - 1.0).abs());
                }
                let kappa = beta * fm;
                let exact = if kappa < 1e-12 { 0.0 } else { vmf_mean_cosine($n, kappa) };
                assert!((c / samples as f64 - exact).abs() < 6e-3, "N={} kappa {kappa}: {} vs {exact}", $n, c / samples as f64);
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
    assert!((mean - exact).abs() < 5.0 * err + 0.005, "{mean} vs {exact}");
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
        assert!(*att > 0 && (*acc as f64) / (*att as f64) > 0.1, "swap {acc}/{att}");
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
    let dos_h = wang_landau_on::<1, 3>([6], false, (-nb - 1e-9, nb + 1e-9, 100), 1.0, 1e-4, 0.7, 401).unwrap();
    let dos_x = wang_landau_on::<1, 2>([6], false, (-nb - 1e-9, nb + 1e-9, 100), 1.0, 1e-4, 0.7, 402).unwrap();
    for beta in [0.3, 0.8, 1.4] {
        let (_, e_h, _) = canonical_from_dos(&dos_h, beta);
        let exact_h = -nb * (1.0 / beta.tanh() - 1.0 / beta);
        assert!((e_h - exact_h).abs() < 0.05 * exact_h.abs(), "Heisenberg beta {beta}: {e_h} vs {exact_h}");
        let (_, e_x, _) = canonical_from_dos(&dos_x, beta);
        let exact_x = -nb * bessel_i(1, beta) / bessel_i(0, beta);
        assert!((e_x - exact_x).abs() < 0.05 * exact_x.abs(), "XY beta {beta}: {e_x} vs {exact_x}");
    }
    assert!(wang_landau_on::<1, 3>([6], false, (1.0, 0.0, 10), 1.0, 1e-3, 0.7, 1).is_err());
}
