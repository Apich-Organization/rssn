//! Tests for the SDE integrators.

use rssn::kernels::random::Rng;
use rssn::kernels::sde::EnsembleScheme;
use rssn::kernels::sde::SdeAdaptiveOptions;
use rssn::kernels::sde::SdeNoncommutativeOptions;
use rssn::kernels::sde::aggregate_noise;
use rssn::kernels::sde::brownian_increments;
use rssn::kernels::sde::ensemble_seed;
use rssn::kernels::sde::euler_maruyama;
use rssn::kernels::sde::heun_stratonovich;
use rssn::kernels::sde::levy_areas;
use rssn::kernels::sde::milstein;
use rssn::kernels::sde::milstein_noncommutative;
use rssn::kernels::sde::noncommutative_noise_grid;
use rssn::kernels::sde::sample_noise;
use rssn::kernels::sde::sde_adaptive_noncommutative_with;
use rssn::kernels::sde::sde_adaptive_with;
use rssn::kernels::sde::sde_ensemble;
use rssn::kernels::sde::srk_platen;

const MU: f64 = 0.8;
const SIGMA: f64 = 0.6;

fn gbm_a(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    o[0] = MU * y[0];
}
fn gbm_g(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    o[0] = SIGMA * y[0];
}
fn coarsen(
    dws: &[Vec<f64>],
    k: usize,
) -> Vec<Vec<f64>> {
    dws.chunks(k)
        .map(|c| vec![c.iter().map(|d| d[0]).sum::<f64>()])
        .collect()
}

type Scheme = fn(
    fn(f64, &[f64], &mut [f64]),
    fn(f64, &[f64], &mut [f64]),
    usize,
    f64,
    f64,
    &[f64],
    &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, rssn::kernels::sde::SdeError>;

fn strong_error(
    scheme: Scheme,
    drift: fn(f64, &[f64], &mut [f64]),
    stratonovich: bool,
    k: usize,
) -> f64 {
    let fine = 1024usize;
    let tf = 1.0;
    let hf = tf / fine as f64;
    let paths = 400;
    let mut rng = Rng::new(2024);
    let mut total = 0.0;
    for _ in 0..paths {
        let dws = brownian_increments(&mut rng, 1, fine, hf);
        let w: f64 = dws.iter().map(|d| d[0]).sum();
        let exact = if stratonovich {
            (MU * tf + SIGMA * w).exp()
        } else {
            ((MU - 0.5 * SIGMA * SIGMA) * tf + SIGMA * w).exp()
        };
        let coarse = coarsen(&dws, k);
        let h = hf * k as f64;
        let path = scheme(drift, gbm_g, 1, 0.0, h, &[1.0], &coarse).unwrap();
        total += (path.last().unwrap()[0] - exact).abs();
    }
    total / paths as f64
}

#[test]
fn strong_convergence_rates() {
    let em: Scheme = euler_maruyama::<fn(f64, &[f64], &mut [f64]), fn(f64, &[f64], &mut [f64])>;
    let mil: Scheme = milstein::<fn(f64, &[f64], &mut [f64]), fn(f64, &[f64], &mut [f64])>;
    let srk: Scheme = srk_platen::<fn(f64, &[f64], &mut [f64]), fn(f64, &[f64], &mut [f64])>;
    let heun: Scheme =
        heun_stratonovich::<fn(f64, &[f64], &mut [f64]), fn(f64, &[f64], &mut [f64])>;
    // k = 64 -> h = 1/16, k = 16 -> h = 1/64
    let r = |s: Scheme, a: fn(f64, &[f64], &mut [f64]), strat: bool| {
        let e1 = strong_error(s, a, strat, 64);
        let e2 = strong_error(s, a, strat, 16);
        (e1 / e2).log2() / 2.0
    };
    let p_em = r(em, gbm_a, false);
    let p_mil = r(mil, gbm_a, false);
    let p_srk = r(srk, gbm_a, false);
    let p_heun = r(heun, gbm_a, true);
    assert!((0.3..0.75).contains(&p_em), "EM rate {p_em}");
    assert!(p_mil > 0.8, "Milstein rate {p_mil}");
    assert!(p_srk > 0.8, "Platen rate {p_srk}");
    assert!(p_heun > 0.8, "Heun rate {p_heun}");
    // Milstein is much more accurate than EM at the same step
    assert!(strong_error(mil, gbm_a, false, 64) < 0.75 * strong_error(em, gbm_a, false, 64));
}

#[test]
fn additive_and_multidimensional_noise() {
    // dX = (0, 0) dt + B dW with constant B: every scheme is exact
    let a = |_t: f64, _y: &[f64], o: &mut [f64]| o.fill(0.0);
    let g =
        |_t: f64, _y: &[f64], o: &mut [f64]| o.copy_from_slice(&[1.0, 0.5, -2.0, 0.25, 0.0, 3.0]);
    let mut rng = Rng::new(5);
    let dws = brownian_increments(&mut rng, 2, 50, 0.02);
    let (w1, w2) = (
        dws.iter().map(|d| d[0]).sum::<f64>(),
        dws.iter().map(|d| d[1]).sum::<f64>(),
    );
    for path in [
        euler_maruyama(a, g, 2, 0.0, 0.02, &[1.0, 2.0, 3.0], &dws).unwrap(),
        milstein(a, g, 2, 0.0, 0.02, &[1.0, 2.0, 3.0], &dws).unwrap(),
        srk_platen(a, g, 2, 0.0, 0.02, &[1.0, 2.0, 3.0], &dws).unwrap(),
    ] {
        let y = path.last().unwrap();
        assert!((y[0] - (1.0 + w1 + 0.5 * w2)).abs() < 1e-12);
        assert!((y[1] - (2.0 - 2.0 * w1 + 0.25 * w2)).abs() < 1e-12);
        assert!((y[2] - (3.0 + 3.0 * w2)).abs() < 1e-12);
    }
    // diagonal commutative noise in two dimensions: independent GBMs
    let a2 = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.1 * y[0];
        o[1] = -0.2 * y[1];
    };
    let g2 = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.5 * y[0];
        o[1] = 0.0;
        o[2] = 0.0;
        o[3] = 0.3 * y[1];
    };
    let n = 2048;
    let h = 1.0 / n as f64;
    let dws = brownian_increments(&mut rng, 2, n, h);
    let (w1, w2) = (
        dws.iter().map(|d| d[0]).sum::<f64>(),
        dws.iter().map(|d| d[1]).sum::<f64>(),
    );
    let y = milstein(a2, g2, 2, 0.0, h, &[1.0, 1.0], &dws).unwrap();
    let last = y.last().unwrap();
    assert!((last[0] - ((0.1 - 0.125) + 0.5 * w1).exp()).abs() < 5e-3);
    assert!((last[1] - ((-0.2 - 0.045) + 0.3 * w2).exp()).abs() < 5e-3);
    assert!(euler_maruyama(a2, g2, 3, 0.0, h, &[1.0, 1.0], &dws).is_err());
}

#[test]
fn adaptive_steps_reproduce_ou_statistics() {
    // Ornstein-Uhlenbeck: dX = -theta X dt + s dW; E X_T = x0 e^{-theta T}, Var = s^2 (1 - e^{-2 theta T}) / (2 theta)
    let theta = 2.0;
    let s = 0.7;
    let a = move |_t: f64, y: &[f64], o: &mut [f64]| o[0] = -theta * y[0];
    let g = move |_t: f64, _y: &[f64], o: &mut [f64]| o[0] = s;
    let tf = 1.0;
    let paths = 1500;
    let (mut sum, mut sumsq) = (0.0, 0.0);
    let mut rejected = 0;
    let mut steps = 0;
    for p in 0..paths {
        let r = sde_adaptive_with(
            a,
            g,
            1,
            0.0,
            tf,
            &[1.5],
            &SdeAdaptiveOptions {
                rtol: 1e-3,
                atol: 1e-3,
                seed: 1000 + p,
                max_steps: 100_000,
            },
        )
        .unwrap();
        assert!((r.t.last().unwrap() - tf).abs() < 1e-12);
        let x = r.y.last().unwrap()[0];
        sum += x;
        sumsq += x * x;
        rejected += r.rejected;
        steps += r.t.len();
    }
    let mean = sum / paths as f64;
    let var = sumsq / paths as f64 - mean * mean;
    let mean_exact = 1.5 * (-theta * tf).exp();
    let var_exact = s * s * (1.0 - (-2.0 * theta * tf).exp()) / (2.0 * theta);
    let se_mean = (var_exact / paths as f64).sqrt();
    assert!(
        (mean - mean_exact).abs() < 5.0 * se_mean,
        "mean {mean} vs {mean_exact}"
    );
    assert!(
        (var / var_exact - 1.0).abs() < 0.15,
        "var {var} vs {var_exact}"
    );
    assert!(steps > 0 && rejected < steps);
    // multiplicative noise: GBM mean E X_T = x0 e^{mu T}; the controller adapts to the state
    let paths = 1500;
    let mut sum = 0.0;
    for p in 0..paths {
        let r = sde_adaptive_with(
            gbm_a,
            gbm_g,
            1,
            0.0,
            1.0,
            &[1.0],
            &SdeAdaptiveOptions {
                rtol: 1e-3,
                atol: 1e-4,
                seed: 77 + p,
                max_steps: 100_000,
            },
        )
        .unwrap();
        sum += r.y.last().unwrap()[0];
    }
    let mean = sum / paths as f64;
    let ex = (MU * 1.0f64).exp();
    let se = ex * ((SIGMA * SIGMA).exp() - 1.0).sqrt() / (paths as f64).sqrt();
    assert!((mean - ex).abs() < 5.0 * se, "{mean} vs {ex}");
    assert!(
        sde_adaptive_with(
            gbm_a,
            gbm_g,
            1,
            1.0,
            0.0,
            &[1.0],
            &SdeAdaptiveOptions {
                rtol: 1e-3,
                atol: 1e-4,
                seed: 1,
                max_steps: 100
            }
        )
        .is_err()
    );
}

#[test]
fn levy_area_statistics() {
    // I_12 given the increments: mean dW1 dW2 / 2, E[(I_12 - dW1 dW2 / 2)^2] / h^2 = E[(1 + xi1^2 + xi2^2)] / 12 = 1/4,
    // and E[(I_12 - dW1 dW2 / 2)^2 xi1^2] / h^2 = E[xi1^2 (1 + xi1^2 + xi2^2)] / 12 = 5/12
    let h = 0.1_f64;
    let m = 3;
    for terms in [1usize, 5, 30] {
        let mut rng = Rng::new(99 + terms as u64);
        let n = 120_000;
        let (mut v, mut w, mut cross, mut sym) = (0.0, 0.0, 0.0, 0.0);
        for _ in 0..n {
            let dw: Vec<f64> = (0..m).map(|_| h.sqrt() * rng.normal()).collect();
            let a = levy_areas(&mut rng, m, h, &dw, terms);
            // exact algebraic identities
            sym += ((a[1] + a[m]) - dw[0] * dw[1]).abs() + (a[0] - 0.5 * (dw[0] * dw[0] - h)).abs();
            let area = (a[1] - 0.5 * dw[0] * dw[1]) / (h);
            v += area * area;
            w += area * area * (dw[0] * dw[0] / h);
            // independence of the areas of different pairs (conditionally uncorrelated)
            cross += (a[1] - 0.5 * dw[0] * dw[1]) * (a[2] - 0.5 * dw[0] * dw[2]) / (h * h);
        }
        let nf = f64::from(n);
        assert!(sym < 1e-9 * nf, "identities, terms {terms}");
        assert!(
            (v / nf - 0.25).abs() < 0.006,
            "terms {terms}: variance {}",
            v / nf
        );
        assert!(
            (w / nf - 5.0 / 12.0).abs() < 0.015,
            "terms {terms}: weighted {}",
            w / nf
        );
        // Cov(A12, A13) = E[xi_2 xi_3 ...] = 0 by symmetry
        assert!(
            (cross / nf).abs() < 0.006,
            "terms {terms}: cross {}",
            cross / nf
        );
    }
    // total variance of I_12 is h^2 / 2
    let mut rng = Rng::new(5);
    let mut v = 0.0;
    let n = 100_000;
    for _ in 0..n {
        let dw: Vec<f64> = (0..2).map(|_| h.sqrt() * rng.normal()).collect();
        let a = levy_areas(&mut rng, 2, h, &dw, 3);
        v += a[1] * a[1];
    }
    assert!((v / f64::from(n) / (h * h) - 0.5).abs() < 0.008);
}

fn nc_drift(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    o[0] = -0.5 * y[0];
    o[1] = -0.3 * y[1] + 0.1;
}

/// Two non-commuting noise fields: the Lie bracket `[g1, g2]` is `(-0.09 y1, 0.09 y2)`-ish and nonzero.
fn nc_diff(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    // row-major n x m: g1 = (0.3 y2, 0.2), g2 = (0.1, 0.3 y1)
    o[0] = 0.3 * y[1];
    o[1] = 0.1;
    o[2] = 0.2;
    o[3] = 0.3 * y[0];
}

#[test]
fn milstein_with_levy_areas_has_strong_order_one_for_noncommutative_noise() {
    let m = 2;
    let fine = 512usize;
    let hf = 1.0 / fine as f64;
    let groups = [8usize, 16, 32];
    let paths = 250;
    let mut rng = Rng::new(31337);
    let mut err_mil = [0.0; 3];
    let mut err_em = [0.0; 3];
    for _ in 0..paths {
        let fine_noise = sample_noise(&mut rng, m, fine, hf, 8);
        let reference =
            milstein_noncommutative(nc_drift, nc_diff, m, 0.0, hf, &[1.0, 0.5], &fine_noise)
                .unwrap();
        let yref = reference.last().unwrap().clone();
        for (gi, &gsz) in groups.iter().enumerate() {
            let coarse = aggregate_noise(&fine_noise, gsz, m);
            let h = hf * gsz as f64;
            let y = milstein_noncommutative(nc_drift, nc_diff, m, 0.0, h, &[1.0, 0.5], &coarse)
                .unwrap();
            let e = euler_maruyama(nc_drift, nc_diff, m, 0.0, h, &[1.0, 0.5], &coarse.dws).unwrap();
            let (yl, el) = (y.last().unwrap(), e.last().unwrap());
            err_mil[gi] +=
                ((yl[0] - yref[0]).powi(2) + (yl[1] - yref[1]).powi(2)) / f64::from(paths);
            err_em[gi] +=
                ((el[0] - yref[0]).powi(2) + (el[1] - yref[1]).powi(2)) / f64::from(paths);
        }
    }
    let rms = |v: [f64; 3]| [v[0].sqrt(), v[1].sqrt(), v[2].sqrt()];
    let (em, ee) = (rms(err_mil), rms(err_em));
    let slope = |e: [f64; 3]| (e[2] / e[0]).ln() / ((groups[2] as f64) / (groups[0] as f64)).ln();
    eprintln!(
        "milstein errors {em:?} slope {}, euler errors {ee:?} slope {}",
        slope(em),
        slope(ee)
    );
    assert!(
        slope(em) > 0.8,
        "Milstein with Levy areas: observed order {}",
        slope(em)
    );
    assert!(
        slope(ee) < 0.7,
        "Euler-Maruyama: observed order {}",
        slope(ee)
    );
    assert!(em[2] * 2.0 < ee[2], "Milstein {} vs Euler {}", em[2], ee[2]);
    // invalid shapes
    let mut bad = sample_noise(&mut rng, m, 4, 0.1, 3);
    bad.areas[0].pop();
    assert!(milstein_noncommutative(nc_drift, nc_diff, m, 0.0, 0.1, &[1.0, 0.5], &bad).is_err());
}

fn nc_options(
    rtol: f64,
    seed: u64,
    cells: usize,
) -> SdeNoncommutativeOptions {
    SdeNoncommutativeOptions {
        base: SdeAdaptiveOptions {
            rtol,
            atol: rtol,
            seed,
            max_steps: 1_000_000,
        },
        cells,
        terms: 8,
    }
}

#[test]
fn adaptive_solver_handles_noncommutative_noise_on_one_brownian_path() {
    let m = 2;
    let cells = 1024;
    let tolerances = [1e-1, 1e-2, 1e-3];
    let paths = 40;
    let mut rms = [0.0; 3];
    let mut mean_steps = [0.0; 3];
    for seed in 0..paths {
        // The fine fixed-step reference sees exactly the path the adaptive solver used.
        let opts0 = nc_options(tolerances[0], 9000 + seed, cells);
        let noise = noncommutative_noise_grid(m, 0.0, 1.0, &opts0);
        let reference = milstein_noncommutative(
            nc_drift,
            nc_diff,
            m,
            0.0,
            1.0 / cells as f64,
            &[1.0, 0.5],
            &noise,
        )
        .unwrap();
        let yref = reference.last().unwrap();
        for (i, &tol) in tolerances.iter().enumerate() {
            let opts = nc_options(tol, 9000 + seed, cells);
            let path = sde_adaptive_noncommutative_with(
                nc_drift,
                nc_diff,
                m,
                0.0,
                1.0,
                &[1.0, 0.5],
                &opts,
            )
            .unwrap();
            assert!((path.t.last().unwrap() - 1.0).abs() < 1e-12);
            let y = path.y.last().unwrap();
            rms[i] +=
                ((y[0] - yref[0]).powi(2) + (y[1] - yref[1]).powi(2)) / f64::from(paths as u32);
            mean_steps[i] += (path.t.len() - 1) as f64 / f64::from(paths as u32);
        }
    }
    let rms = rms.map(f64::sqrt);
    eprintln!("adaptive nc: rms {rms:?}, mean steps {mean_steps:?}");
    // tighter tolerance: more steps, smaller error, and strong order about one in the step size
    assert!(mean_steps[0] < mean_steps[1] && mean_steps[1] < mean_steps[2]);
    assert!(rms[0] > rms[1] && rms[1] > rms[2], "{rms:?}");
    let order = (rms[2] / rms[0]).ln() / (mean_steps[0] / mean_steps[2]).ln();
    assert!(order > 0.8, "observed strong order {order}");
    assert!(rms[2] < 6e-3, "{rms:?}");
}

#[test]
fn noncommutative_adaptive_path_does_not_depend_on_the_step_sequence() {
    // Two tolerances take different steps but must see the same Brownian path: with
    // additive-free noise g = 0 the state is the drift flow, so instead compare the
    // *sum of the area-free part*: the endpoint of X1 for g1 = const is W1(T).
    let m = 2;
    let cells = 512;
    let diff = |_t: f64, _y: &[f64], o: &mut [f64]| {
        o.copy_from_slice(&[1.0, 0.0, 0.0, 1.0]);
    };
    let drift = |_t: f64, _y: &[f64], o: &mut [f64]| o.fill(0.0);
    let mut endpoints = Vec::new();
    for tol in [1e-1, 1e-4] {
        let opts = nc_options(tol, 77, cells);
        let path =
            sde_adaptive_noncommutative_with(drift, diff, m, 0.0, 1.0, &[0.0, 0.0], &opts).unwrap();
        endpoints.push(path.y.last().unwrap().clone());
    }
    assert!((endpoints[0][0] - endpoints[1][0]).abs() < 1e-12);
    assert!((endpoints[0][1] - endpoints[1][1]).abs() < 1e-12);
    let noise = noncommutative_noise_grid(m, 0.0, 1.0, &nc_options(1e-1, 77, cells));
    let w1: f64 = noise.dws.iter().map(|d| d[0]).sum();
    assert!((w1 - endpoints[0][0]).abs() < 1e-12);
    // Chen: merging all cells reproduces I_01 + I_10 = W0 W1.
    let all = aggregate_noise(&noise, cells, m);
    let (ar, dw) = (&all.areas[0], &all.dws[0]);
    assert!((ar[1] + ar[2] - dw[0] * dw[1]).abs() < 1e-10);
    // invalid
    assert!(
        sde_adaptive_noncommutative_with(
            drift,
            diff,
            m,
            0.0,
            1.0,
            &[0.0, 0.0],
            &nc_options(1e-2, 1, 0)
        )
        .is_err()
    );
}

#[test]
fn commutative_adaptive_results_are_unchanged() {
    // The original adaptive solver is untouched: same seed, same path.
    let opts = SdeAdaptiveOptions {
        rtol: 1e-3,
        atol: 1e-6,
        seed: 5,
        max_steps: 100_000,
    };
    let a = sde_adaptive_with(gbm_a, gbm_g, 1, 0.0, 1.0, &[1.0], &opts).unwrap();
    let b = sde_adaptive_with(gbm_a, gbm_g, 1, 0.0, 1.0, &[1.0], &opts).unwrap();
    assert_eq!(a, b);
    assert!(a.t.len() > 2);
}

#[test]
fn ensemble_members_depend_only_on_seed_and_index() {
    let (t0, h, steps) = (0.0, 0.01, 100);
    let run = |scheme, paths| {
        sde_ensemble(gbm_a, gbm_g, 1, (t0, h, steps), &[1.0], paths, 42, scheme).unwrap()
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    for scheme in [
        EnsembleScheme::EulerMaruyama,
        EnsembleScheme::Milstein,
        EnsembleScheme::SrkPlaten,
        EnsembleScheme::HeunStratonovich,
    ] {
        let parallel = run(scheme, 64);
        let serial = one.install(|| run(scheme, 64));
        assert_eq!(parallel, serial, "{scheme:?}");
        // a prefix of a larger ensemble is the smaller ensemble
        assert_eq!(run(scheme, 9), parallel[..9].to_vec(), "{scheme:?}");
        // and equals the explicit single-path computation with the derived seed
        let mut rng = Rng::new(ensemble_seed(42, 5));
        let dws = brownian_increments(&mut rng, 1, steps, h);
        let path = match scheme {
            | EnsembleScheme::EulerMaruyama => euler_maruyama(gbm_a, gbm_g, 1, t0, h, &[1.0], &dws),
            | EnsembleScheme::Milstein => milstein(gbm_a, gbm_g, 1, t0, h, &[1.0], &dws),
            | EnsembleScheme::SrkPlaten => srk_platen(gbm_a, gbm_g, 1, t0, h, &[1.0], &dws),
            | _ => heun_stratonovich(gbm_a, gbm_g, 1, t0, h, &[1.0], &dws),
        }
        .unwrap();
        assert_eq!(path.last().unwrap(), &parallel[5]);
    }
    // Euler-Maruyama ensemble mean of geometric Brownian motion: E X_T = exp(mu T)
    let ends = sde_ensemble(
        gbm_a,
        gbm_g,
        1,
        (t0, 0.01, 100),
        &[1.0],
        4000,
        7,
        EnsembleScheme::Milstein,
    )
    .unwrap();
    let mean = ends.iter().map(|v| v[0]).sum::<f64>() / ends.len() as f64;
    assert!((mean - (MU * 1.0f64).exp()).abs() < 0.05, "mean {mean}");
    // non-commutative scheme runs and is reproducible
    let nc = |paths| {
        sde_ensemble(
            nc_drift,
            nc_diff,
            2,
            (0.0, 0.02, 50),
            &[1.0, 0.5],
            paths,
            3,
            EnsembleScheme::MilsteinNoncommutative { terms: 4 },
        )
        .unwrap()
    };
    assert_eq!(nc(16), one.install(|| nc(16)));
}
