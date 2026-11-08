//! Tests for the non-commutative / Stratonovich SDE schemes of `rssn::kernels::sde`:
//! Stratonovich Milstein, the derivative-free SRK schemes, and the adaptive
//! Stratonovich solvers (exact Ito conversion).

use rssn::kernels::random::Rng;
use rssn::kernels::sde::NoiseSample;
use rssn::kernels::sde::SdeAdaptiveOptions;
use rssn::kernels::sde::SdeNoncommutativeOptions;
use rssn::kernels::sde::aggregate_noise;
use rssn::kernels::sde::euler_maruyama;
use rssn::kernels::sde::milstein_noncommutative;
use rssn::kernels::sde::milstein_noncommutative_stratonovich;
use rssn::kernels::sde::noncommutative_noise_grid;
use rssn::kernels::sde::sample_noise;
use rssn::kernels::sde::sde_adaptive_noncommutative_stratonovich_with;
use rssn::kernels::sde::sde_adaptive_stratonovich_with;
use rssn::kernels::sde::srk_noncommutative;
use rssn::kernels::sde::srk_noncommutative_stratonovich;
use rssn::kernels::sde::srk_platen;
use rssn::kernels::sde::stratonovich_to_ito_drift;

const A: [[f64; 2]; 2] = [[-0.5, 0.2], [0.0, -0.3]];
const B1: [[f64; 2]; 2] = [[0.0, 0.5], [0.0, 0.0]];
const B2: [[f64; 2]; 2] = [[0.0, 0.0], [0.5, 0.0]];

fn mat_vec(
    m: &[[f64; 2]; 2],
    y: &[f64],
) -> [f64; 2] {
    [
        m[0][0] * y[0] + m[0][1] * y[1],
        m[1][0] * y[0] + m[1][1] * y[1],
    ]
}

/// `dX = A X dt + B1 X dW1 + B2 X dW2`, with `[B1, B2] != 0` (the noise does not commute).
fn lin_drift(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    o.copy_from_slice(&mat_vec(&A, y));
}

fn lin_diff(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    let (u, v) = (mat_vec(&B1, y), mat_vec(&B2, y));
    o[0] = u[0];
    o[1] = v[0];
    o[2] = u[1];
    o[3] = v[1];
}

/// The Ito drift equivalent to the Stratonovich drift `A X`:
/// `A X + (B1^2 + B2^2) X / 2`, in closed form.
fn lin_ito_drift(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    let ax = mat_vec(&A, y);
    let b1 = mat_vec(&B1, &mat_vec(&B1, y));
    let b2 = mat_vec(&B2, &mat_vec(&B2, y));
    for i in 0..2 {
        o[i] = ax[i] + 0.5 * (b1[i] + b2[i]);
    }
}

type Scheme = fn(
    fn(f64, &[f64], &mut [f64]),
    fn(f64, &[f64], &mut [f64]),
    usize,
    f64,
    f64,
    &[f64],
    &NoiseSample,
) -> Result<Vec<Vec<f64>>, rssn::kernels::sde::SdeError>;

type Cb = fn(f64, &[f64], &mut [f64]);

/// RMS strong errors (levels `groups`) of `scheme` against `reference` on one fine Brownian path
/// with Levy areas, the coarse schemes seeing the aggregated path.
fn rms_errors(
    scheme: Scheme,
    reference: Scheme,
    groups: &[usize],
    paths: u32,
    seed: u64,
) -> Vec<f64> {
    rms_errors_for(
        (lin_drift, lin_diff),
        scheme,
        reference,
        groups,
        paths,
        seed,
    )
}

/// As [`rms_errors`] for the coefficients `(drift, diff)`.
fn rms_errors_for(
    (drift, diff): (Cb, Cb),
    scheme: Scheme,
    reference: Scheme,
    groups: &[usize],
    paths: u32,
    seed: u64,
) -> Vec<f64> {
    let (m, fine) = (2, 512_usize);
    let hf = 1.0 / fine as f64;
    let y0 = [1.0, 0.5];
    let mut rng = Rng::new(seed);
    let mut acc = vec![0.0; groups.len()];
    for _ in 0..paths {
        let noise = sample_noise(&mut rng, m, fine, hf, 8);
        let yref = reference(drift, diff, m, 0.0, hf, &y0, &noise)
            .unwrap()
            .pop()
            .unwrap();
        for (slot, &g) in acc.iter_mut().zip(groups) {
            let coarse = aggregate_noise(&noise, g, m);
            let y = scheme(drift, diff, m, 0.0, hf * g as f64, &y0, &coarse)
                .unwrap()
                .pop()
                .unwrap();
            *slot += ((y[0] - yref[0]).powi(2) + (y[1] - yref[1]).powi(2)) / f64::from(paths);
        }
    }
    acc.into_iter().map(f64::sqrt).collect()
}

fn slope(
    errors: &[f64],
    groups: &[usize],
) -> f64 {
    let (e0, e1) = (errors[0], errors[errors.len() - 1]);
    let (g0, g1) = (groups[0] as f64, groups[groups.len() - 1] as f64);
    (e1 / e0).ln() / (g1 / g0).ln()
}

const GROUPS: [usize; 3] = [8, 16, 32];

#[test]
fn stratonovich_milstein_and_srk_have_strong_order_one_for_noncommutative_noise() {
    let mil_s: Scheme = milstein_noncommutative_stratonovich::<Cb, Cb>;
    let srk_s: Scheme = srk_noncommutative_stratonovich::<Cb, Cb>;
    let e_mil = rms_errors(mil_s, mil_s, &GROUPS, 150, 4242);
    let e_srk = rms_errors(srk_s, mil_s, &GROUPS, 150, 4242);
    let (p_mil, p_srk) = (slope(&e_mil, &GROUPS), slope(&e_srk, &GROUPS));
    eprintln!("Stratonovich Milstein {e_mil:?} slope {p_mil}; SRK {e_srk:?} slope {p_srk}");
    assert!(p_mil > 0.8, "Stratonovich Milstein order {p_mil}");
    assert!(p_srk > 0.8, "Stratonovich SRK order {p_srk}");
}

#[test]
fn ito_srk_has_strong_order_one_and_euler_maruyama_one_half() {
    let mil: Scheme = milstein_noncommutative::<Cb, Cb>;
    let srk: Scheme = srk_noncommutative::<Cb, Cb>;
    let em: Scheme = |a, g, m, t0, h, y0, noise| euler_maruyama(a, g, m, t0, h, y0, &noise.dws);
    let e_srk = rms_errors(srk, mil, &GROUPS, 150, 777);
    let e_em = rms_errors(em, mil, &GROUPS, 150, 777);
    let (p_srk, p_em) = (slope(&e_srk, &GROUPS), slope(&e_em, &GROUPS));
    eprintln!("Ito SRK {e_srk:?} slope {p_srk}; EM {e_em:?} slope {p_em}");
    assert!(p_srk > 0.8, "Ito SRK order {p_srk}");
    assert!((0.3..0.75).contains(&p_em), "Euler-Maruyama order {p_em}");
    assert!(
        e_srk[2] * 2.0 < e_em[2],
        "SRK {} vs EM {}",
        e_srk[2],
        e_em[2]
    );
}

/// A genuinely nonlinear non-commuting pair (`g''` does not vanish), where the SRK schemes
/// differ from Milstein: the Stratonovich SRK needs its central supporting values to reach
/// order one.
fn nl_drift(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    o[0] = -y[0];
    o[1] = -0.5 * y[1] + 0.2;
}

fn nl_diff(
    _t: f64,
    y: &[f64],
    o: &mut [f64],
) {
    // row-major n x m: g1 = (0.5 sin(y1), 0.2 + 0.1 y0), g2 = (0.1, 0.4 tanh(y0))
    o[0] = 0.5 * y[1].sin();
    o[1] = 0.1;
    o[2] = 0.2 + 0.1 * y[0];
    o[3] = 0.4 * y[0].tanh();
}

#[test]
fn srk_schemes_have_strong_order_one_for_nonlinear_noncommutative_noise() {
    let coefficients: (Cb, Cb) = (nl_drift, nl_diff);
    let mil: Scheme = milstein_noncommutative::<Cb, Cb>;
    let mil_s: Scheme = milstein_noncommutative_stratonovich::<Cb, Cb>;
    let srk: Scheme = srk_noncommutative::<Cb, Cb>;
    let srk_s: Scheme = srk_noncommutative_stratonovich::<Cb, Cb>;
    let e_srk = rms_errors_for(coefficients, srk, mil, &GROUPS, 120, 8);
    let e_srk_s = rms_errors_for(coefficients, srk_s, mil_s, &GROUPS, 120, 8);
    let (p, p_s) = (slope(&e_srk, &GROUPS), slope(&e_srk_s, &GROUPS));
    eprintln!("nonlinear Ito SRK {e_srk:?} slope {p}; Stratonovich SRK {e_srk_s:?} slope {p_s}");
    assert!(p > 0.8, "Ito SRK order {p}");
    assert!(p_s > 0.8, "Stratonovich SRK order {p_s}");
}

/// A Stratonovich scheme on `(a, g)` and the Ito scheme on the converted drift `(a~, g)` solve the
/// same SDE, so on one path they agree up to their discretisation errors.
#[test]
fn stratonovich_schemes_agree_with_ito_schemes_on_the_converted_drift() {
    let (m, n) = (2, 256);
    let h = 1.0 / n as f64;
    let mut rng = Rng::new(99);
    let noise = sample_noise(&mut rng, m, n, h, 8);
    let y0 = [1.0, 0.5];
    let strat = milstein_noncommutative_stratonovich(lin_drift, lin_diff, m, 0.0, h, &y0, &noise)
        .unwrap()
        .pop()
        .unwrap();
    let ito = milstein_noncommutative(lin_ito_drift, lin_diff, m, 0.0, h, &y0, &noise)
        .unwrap()
        .pop()
        .unwrap();
    let srk = srk_noncommutative_stratonovich(lin_drift, lin_diff, m, 0.0, h, &y0, &noise)
        .unwrap()
        .pop()
        .unwrap();
    for i in 0..2 {
        assert!((strat[i] - ito[i]).abs() < 2e-2, "{strat:?} vs {ito:?}");
        assert!((srk[i] - ito[i]).abs() < 2e-2, "{srk:?} vs {ito:?}");
    }
    // the closure-based conversion equals the closed form
    let (da, dg): (Cb, Cb) = (lin_drift, lin_diff);
    let conv = stratonovich_to_ito_drift(&da, &dg, m);
    let (mut out, mut exact) = ([0.0; 2], [0.0; 2]);
    conv(0.0, &[0.7, -1.3], &mut out);
    lin_ito_drift(0.0, &[0.7, -1.3], &mut exact);
    assert!((out[0] - exact[0]).abs() < 1e-7 && (out[1] - exact[1]).abs() < 1e-7);
}

/// On commutative noise (diagonal geometric noise) the Levy areas drop out and the new SRK is the
/// classical one up to the `O(h)` of the supporting values.
#[test]
fn srk_noncommutative_matches_srk_platen_for_commutative_noise() {
    let a = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.3 * y[0];
        o[1] = -0.2 * y[1];
    };
    let g = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.4 * y[0];
        o[1] = 0.0;
        o[2] = 0.0;
        o[3] = 0.25 * y[1];
    };
    let (m, n) = (2, 1024);
    let h = 1.0 / n as f64;
    let mut rng = Rng::new(5150);
    let noise = sample_noise(&mut rng, m, n, h, 6);
    let y0 = [1.0, 2.0];
    let nc = srk_noncommutative(a, g, m, 0.0, h, &y0, &noise)
        .unwrap()
        .pop()
        .unwrap();
    let cm = srk_platen(a, g, m, 0.0, h, &y0, &noise.dws)
        .unwrap()
        .pop()
        .unwrap();
    for i in 0..2 {
        assert!(
            (nc[i] - cm[i]).abs() < 5e-3 * cm[i].abs().max(1.0),
            "component {i}: {} vs {}",
            nc[i],
            cm[i]
        );
    }
}

#[test]
fn adaptive_stratonovich_hits_the_exact_pathwise_solution() {
    // Two commuting diagonal channels in the non-commutative controller: pathwise against
    // X_i(T) = X_i(0) exp(sigma_i W_i(T)).
    let sigma = [0.5, 0.3];
    let drift = |_t: f64, _y: &[f64], o: &mut [f64]| o.fill(0.0);
    let diff = move |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = sigma[0] * y[0];
        o[1] = 0.0;
        o[2] = 0.0;
        o[3] = sigma[1] * y[1];
    };
    let opts = SdeNoncommutativeOptions {
        base: SdeAdaptiveOptions {
            rtol: 1e-3,
            atol: 1e-6,
            seed: 31,
            max_steps: 1_000_000,
        },
        cells: 8192,
        terms: 6,
    };
    let noise = noncommutative_noise_grid(2, 0.0, 1.0, &opts);
    let w: Vec<f64> = (0..2)
        .map(|j| noise.dws.iter().map(|d| d[j]).sum())
        .collect();
    let path =
        sde_adaptive_noncommutative_stratonovich_with(drift, diff, 2, 0.0, 1.0, &[1.0, 2.0], &opts)
            .unwrap();
    let end = path.y.last().unwrap();
    let exact = [(sigma[0] * w[0]).exp(), 2.0 * (sigma[1] * w[1]).exp()];
    for i in 0..2 {
        assert!(
            (end[i] - exact[i]).abs() < 2e-2 * exact[i],
            "component {i}: {} vs {}",
            end[i],
            exact[i]
        );
    }
    assert!((path.t.last().unwrap() - 1.0).abs() < 1e-12);
    // The Ito interpretation of the same coefficients is a different process: it misses the
    // exp(sigma^2 T / 2) factor.
    let ito = rssn::kernels::sde::sde_adaptive_noncommutative_with(
        drift,
        diff,
        2,
        0.0,
        1.0,
        &[1.0, 2.0],
        &opts,
    )
    .unwrap();
    let ito_end = ito.y.last().unwrap();
    let ratio = end[0] / ito_end[0];
    assert!(
        (ratio - (0.5 * sigma[0] * sigma[0]).exp()).abs() < 2e-2,
        "ratio {ratio}"
    );
}

#[test]
fn adaptive_stratonovich_commutative_has_the_stratonovich_mean() {
    // dX = sigma X o dW: E X_T = exp(sigma^2 T / 2) (the Ito mean would be 1).
    let sigma = 0.5;
    let (paths, tf) = (600_u64, 1.0);
    let drift = |_t: f64, _y: &[f64], o: &mut [f64]| o[0] = 0.0;
    let diff = move |_t: f64, y: &[f64], o: &mut [f64]| o[0] = sigma * y[0];
    let (mut s, mut s2) = (0.0, 0.0);
    for p in 0..paths {
        let r = sde_adaptive_stratonovich_with(
            drift,
            diff,
            1,
            0.0,
            tf,
            &[1.0],
            &SdeAdaptiveOptions {
                rtol: 1e-3,
                atol: 1e-4,
                seed: 500 + p,
                max_steps: 100_000,
            },
        )
        .unwrap();
        let x = r.y.last().unwrap()[0];
        s += x;
        s2 += x * x;
    }
    let np = paths as f64;
    let mean = s / np;
    let var = s2 / np - mean * mean;
    let se = (var / np).sqrt();
    let strat = (0.5 * sigma * sigma * tf).exp();
    assert!(
        (mean - strat).abs() < 4.0 * se,
        "mean {mean}, Stratonovich {strat}, se {se}"
    );
    assert!((mean - 1.0).abs() > 3.0 * se, "mean {mean} is the Ito mean");
}

#[test]
fn error_paths() {
    let mut rng = Rng::new(1);
    let noise = sample_noise(&mut rng, 2, 4, 0.1, 3);
    let y0 = [1.0, 0.5];
    for scheme in [
        milstein_noncommutative_stratonovich::<Cb, Cb> as Scheme,
        srk_noncommutative::<Cb, Cb> as Scheme,
        srk_noncommutative_stratonovich::<Cb, Cb> as Scheme,
    ] {
        // zero noise dimension, empty state, wrong noise dimension
        assert!(scheme(lin_drift, lin_diff, 0, 0.0, 0.1, &y0, &noise).is_err());
        assert!(scheme(lin_drift, lin_diff, 2, 0.0, 0.1, &[], &noise).is_err());
        assert!(scheme(lin_drift, lin_diff, 3, 0.0, 0.1, &y0, &noise).is_err());
        // mismatched area block
        let mut bad = noise.clone();
        bad.areas[1].pop();
        assert!(scheme(lin_drift, lin_diff, 2, 0.0, 0.1, &y0, &bad).is_err());
        // fewer area blocks than increments
        let mut short = noise.clone();
        short.areas.pop();
        assert!(scheme(lin_drift, lin_diff, 2, 0.0, 0.1, &y0, &short).is_err());
        // a valid run is finite and has steps + 1 states
        let ok = scheme(lin_drift, lin_diff, 2, 0.0, 0.1, &y0, &noise).unwrap();
        assert_eq!(ok.len(), 5);
        assert!(ok.iter().flatten().all(|v| v.is_finite()));
    }
    let opts = SdeNoncommutativeOptions {
        base: SdeAdaptiveOptions {
            rtol: 1e-3,
            atol: 1e-6,
            seed: 1,
            max_steps: 100,
        },
        cells: 0,
        terms: 3,
    };
    assert!(
        sde_adaptive_noncommutative_stratonovich_with(lin_drift, lin_diff, 2, 0.0, 1.0, &y0, &opts)
            .is_err()
    );
    let base = SdeAdaptiveOptions {
        rtol: 1e-3,
        atol: 1e-6,
        seed: 1,
        max_steps: 100,
    };
    assert!(sde_adaptive_stratonovich_with(lin_drift, lin_diff, 0, 0.0, 1.0, &y0, &base).is_err());
    assert!(sde_adaptive_stratonovich_with(lin_drift, lin_diff, 2, 1.0, 0.0, &y0, &base).is_err());
}
