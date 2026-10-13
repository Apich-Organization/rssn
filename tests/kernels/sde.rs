//! Tests for the SDE integrators.

use rssn::kernels::random::Rng;
use rssn::kernels::sde::{brownian_increments, euler_maruyama, heun_stratonovich, milstein, sde_adaptive, srk_platen};

const MU: f64 = 0.8;
const SIGMA: f64 = 0.6;

fn gbm_a(_t: f64, y: &[f64], o: &mut [f64]) {
    o[0] = MU * y[0];
}
fn gbm_g(_t: f64, y: &[f64], o: &mut [f64]) {
    o[0] = SIGMA * y[0];
}
fn coarsen(dws: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    dws.chunks(k).map(|c| vec![c.iter().map(|d| d[0]).sum::<f64>()]).collect()
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

fn strong_error(scheme: Scheme, drift: fn(f64, &[f64], &mut [f64]), stratonovich: bool, k: usize) -> f64 {
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
    let heun: Scheme = heun_stratonovich::<fn(f64, &[f64], &mut [f64]), fn(f64, &[f64], &mut [f64])>;
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
    let g = |_t: f64, _y: &[f64], o: &mut [f64]| o.copy_from_slice(&[1.0, 0.5, -2.0, 0.25, 0.0, 3.0]);
    let mut rng = Rng::new(5);
    let dws = brownian_increments(&mut rng, 2, 50, 0.02);
    let (w1, w2) = (dws.iter().map(|d| d[0]).sum::<f64>(), dws.iter().map(|d| d[1]).sum::<f64>());
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
    let (w1, w2) = (dws.iter().map(|d| d[0]).sum::<f64>(), dws.iter().map(|d| d[1]).sum::<f64>());
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
        let r = sde_adaptive(a, g, 1, 0.0, tf, &[1.5], 1e-3, 1e-3, 1000 + p, 100_000).unwrap();
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
    assert!((mean - mean_exact).abs() < 5.0 * se_mean, "mean {mean} vs {mean_exact}");
    assert!((var / var_exact - 1.0).abs() < 0.15, "var {var} vs {var_exact}");
    assert!(steps > 0 && rejected < steps);
    // multiplicative noise: GBM mean E X_T = x0 e^{mu T}; the controller adapts to the state
    let paths = 1500;
    let mut sum = 0.0;
    for p in 0..paths {
        let r = sde_adaptive(gbm_a, gbm_g, 1, 0.0, 1.0, &[1.0], 1e-3, 1e-4, 77 + p, 100_000).unwrap();
        sum += r.y.last().unwrap()[0];
    }
    let mean = sum / paths as f64;
    let ex = (MU * 1.0f64).exp();
    let se = ex * ((SIGMA * SIGMA).exp() - 1.0).sqrt() / (paths as f64).sqrt();
    assert!((mean - ex).abs() < 5.0 * se, "{mean} vs {ex}");
    assert!(sde_adaptive(gbm_a, gbm_g, 1, 1.0, 0.0, &[1.0], 1e-3, 1e-4, 1, 100).is_err());
}
