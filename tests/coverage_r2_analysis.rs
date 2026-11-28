//! Round-2 coverage: Levin oscillatory quadrature, the legacy adaptive SDE
//! entry point, `coulomb_g_from_u`, the unchecked coordinate conversions and
//! the analysis / synthesis filter accessors of the wavelet kernel.
//! Public API only; deterministic inputs.
// The deprecated `sde_adaptive` wrapper is tested on purpose (bit-identical to `sde_adaptive_with`).
#![allow(deprecated)]

use num_complex::Complex64 as Cx;
use rssn::kernels::coordinates::CoordinateSystem as C;
use rssn::kernels::coordinates::from_cartesian;
use rssn::kernels::coordinates::from_cartesian_pure;
use rssn::kernels::coordinates::to_cartesian;
use rssn::kernels::coordinates::to_cartesian_pure;
use rssn::kernels::coordinates::transform_point;
use rssn::kernels::coordinates::transform_point_pure;
use rssn::kernels::oscillatory::Oscillatory;
use rssn::kernels::oscillatory::levin;
use rssn::kernels::oscillatory::levin_adaptive;
use rssn::kernels::sde::SdeAdaptiveOptions;
use rssn::kernels::sde::SdeError;
use rssn::kernels::sde::sde_adaptive;
use rssn::kernels::sde::sde_adaptive_with;
use rssn::kernels::special_ext::coulomb_f;
use rssn::kernels::special_ext::coulomb_f_prime;
use rssn::kernels::special_ext::coulomb_g;
use rssn::kernels::special_ext::coulomb_g_from_u;
use rssn::kernels::special_ext::coulomb_g_prime;
use rssn::kernels::wavelets::Wavelet;

const ALL_WAVELETS: [Wavelet; 21] = [
    Wavelet::Haar,
    Wavelet::Db2,
    Wavelet::Db3,
    Wavelet::Db4,
    Wavelet::Db5,
    Wavelet::Db6,
    Wavelet::Db7,
    Wavelet::Db8,
    Wavelet::Db9,
    Wavelet::Db10,
    Wavelet::Sym4,
    Wavelet::Sym5,
    Wavelet::Sym6,
    Wavelet::Bior1_1,
    Wavelet::Bior1_3,
    Wavelet::Bior1_5,
    Wavelet::Bior2_2,
    Wavelet::Bior2_4,
    Wavelet::Bior3_3,
    Wavelet::Bior3_5,
    Wavelet::Bior4_4,
];

fn dist(
    r: &Oscillatory,
    exact: Cx,
) -> f64 {
    (Cx::new(r.re, r.im) - exact).norm()
}

// ------------------------------------------------------------ oscillatory

#[test]
fn levin_exponential_amplitude_and_exponential_phase_closed_forms() {
    for w in [5.0, 80.0, 2.0e3] {
        // int_0^1 e^x e^{i w x} dx = (e^{1 + i w} - 1) / (1 + i w)
        let iw = Cx::new(0.0, w);
        let exact = ((Cx::new(1.0, 0.0) + iw).exp() - 1.0) / (Cx::new(1.0, 0.0) + iw);
        let r = levin_adaptive(&f64::exp, &|x| x, &|_| 1.0, 0.0, 1.0, w, 1e-12, 128).unwrap();
        assert!(dist(&r, exact) < 1e-11, "w = {w}: {r:?} vs {exact}");
        // the reported estimate is finite, non-negative and honest
        assert!(
            r.error.is_finite() && r.error >= 0.0 && r.error <= 1e-12 * r.re.hypot(r.im).max(1.0)
        );
        // int_0^1 e^x exp(i w e^x) dx = (exp(i w e) - exp(i w)) / (i w)
        let exact = ((iw * std::f64::consts::E).exp() - iw.exp()) / iw;
        let r = levin_adaptive(&f64::exp, &f64::exp, &f64::exp, 0.0, 1.0, w, 1e-12, 128).unwrap();
        assert!(dist(&r, exact) < 1e-11, "w = {w}: {r:?} vs {exact}");
    }
    // cubic phase: int_1^2 3 x^2 exp(i w x^3) dx = (exp(8 i w) - exp(i w)) / (i w)
    for w in [3.0, 400.0] {
        let iw = Cx::new(0.0, w);
        let exact = ((iw * 8.0).exp() - iw.exp()) / iw;
        let r = levin_adaptive(
            &|x| 3.0 * x * x,
            &|x| x * x * x,
            &|x| 3.0 * x * x,
            1.0,
            2.0,
            w,
            1e-12,
            128,
        )
        .unwrap();
        assert!(dist(&r, exact) < 1e-10, "w = {w}");
    }
}

#[test]
fn levin_symmetries_linearity_and_additivity() {
    let f = |x: f64| 1.0 / (1.0 + x * x);
    let g = |x: f64| x + 0.25 * x * x;
    let dg = |x: f64| 1.0 + 0.5 * x;
    let run = |a: f64, b: f64, w: f64| levin_adaptive(&f, &g, &dg, a, b, w, 1e-12, 256).unwrap();
    let w = 60.0;
    let whole = run(0.0, 2.0, w);
    // additivity over a split of the interval
    let (left, right) = (run(0.0, 0.7, w), run(0.7, 2.0, w));
    assert!((whole.re - left.re - right.re).abs() < 1e-10);
    assert!((whole.im - left.im - right.im).abs() < 1e-10);
    // orientation: int_b^a = - int_a^b
    let back = run(2.0, 0.0, w);
    assert!((back.re + whole.re).abs() < 1e-10 && (back.im + whole.im).abs() < 1e-10);
    // real f, w -> -w : complex conjugate
    let neg = run(0.0, 2.0, -w);
    assert!((neg.re - whole.re).abs() < 1e-10 && (neg.im + whole.im).abs() < 1e-10);
    // linearity in the amplitude
    let r2 = levin_adaptive(&|x| 3.0 * f(x) - 2.0, &g, &dg, 0.0, 2.0, w, 1e-12, 256).unwrap();
    let r_const = levin_adaptive(&|_| 1.0, &g, &dg, 0.0, 2.0, w, 1e-12, 256).unwrap();
    assert!((r2.re - (3.0 * whole.re - 2.0 * r_const.re)).abs() < 1e-9);
    assert!((r2.im - (3.0 * whole.im - 2.0 * r_const.im)).abs() < 1e-9);
    // a = b gives exactly zero
    let z = run(1.5, 1.5, w);
    assert_eq!((z.re, z.im), (0.0, 0.0));
}

#[test]
fn levin_converges_spectrally_in_the_number_of_points() {
    // smooth non-polynomial amplitude, w = 30: the error must drop quickly with n
    let w = 30.0;
    let f = |x: f64| (2.0 * x).cos() + 0.5;
    let reference = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, w, 64).unwrap();
    let mut errors = Vec::new();
    for n in [6usize, 10, 14, 20] {
        let r = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, w, n).unwrap();
        assert!(r.error.is_nan(), "plain levin makes no error estimate");
        errors.push((Cx::new(r.re, r.im) - Cx::new(reference.re, reference.im)).norm());
    }
    assert!(errors[0] > errors[1] && errors[1] > errors[2], "{errors:?}");
    assert!(errors[3] < 1e-12, "{errors:?}");
    // ... and the accuracy does not degrade as the frequency grows
    let low = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, 30.0, 14).unwrap();
    let high = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, 3.0e5, 14).unwrap();
    let ref_low = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, 30.0, 40).unwrap();
    let ref_high = levin(&f, &|x| x, &|_| 1.0, 0.0, 1.0, 3.0e5, 40).unwrap();
    let e_low = (low.re - ref_low.re).hypot(low.im - ref_low.im);
    let e_high = (high.re - ref_high.re).hypot(high.im - ref_high.im);
    assert!(e_high <= e_low + 1e-14, "{e_low} vs {e_high}");
    // the integral itself decays like 1/w for a non-vanishing amplitude at the ends
    assert!(high.re.hypot(high.im) < 1e-4);
}

#[test]
fn levin_adaptive_reports_failures_instead_of_wrong_numbers() {
    let one = |_: f64| 1.0;
    let id = |x: f64| x;
    // phase with a stationary point inside: no convergence
    assert!(
        levin_adaptive(
            &one,
            &|x| x * x * x,
            &|x| 3.0 * x * x,
            -1.0,
            1.0,
            500.0,
            1e-12,
            64
        )
        .is_err()
    );
    // an unreachable tolerance: error rather than a silent best effort
    assert!(levin_adaptive(&|x| (50.0 * x).sin(), &id, &one, 0.0, 1.0, 1.0, 1e-16, 16).is_err());
    // NaN tolerance never compares as converged
    assert!(levin_adaptive(&one, &id, &one, 0.0, 1.0, 10.0, f64::NAN, 32).is_err());
    // too small a point budget
    assert!(levin_adaptive(&one, &id, &one, 0.0, 1.0, 10.0, 1e-10, 4).is_err());
    // non-finite data
    for (a, b, w) in [
        (f64::NAN, 1.0, 5.0),
        (0.0, f64::INFINITY, 5.0),
        (0.0, 1.0, f64::NAN),
        (0.0, 1.0, f64::NEG_INFINITY),
    ] {
        assert!(
            levin(&one, &id, &one, a, b, w, 8).is_err(),
            "a={a} b={b} w={w}"
        );
        assert!(levin_adaptive(&one, &id, &one, a, b, w, 1e-10, 64).is_err());
    }
    assert!(levin(&|_| f64::INFINITY, &id, &one, 0.0, 1.0, 5.0, 8).is_err());
    assert!(levin(&one, &id, &|_| f64::NAN, 0.0, 1.0, 5.0, 8).is_err());
    // fewer than two points
    for n in [0usize, 1] {
        assert!(levin(&one, &id, &one, 0.0, 1.0, 5.0, n).is_err());
    }
    // two points is the minimum that works: finite output
    let r = levin(&one, &id, &one, 0.0, 1.0, 5.0, 2).unwrap();
    assert!(r.re.is_finite() && r.im.is_finite());
}

#[test]
fn levin_with_zero_frequency_is_the_plain_integral() {
    // w = 0 makes the collocation system singular; the minimum-norm solution is the
    // antiderivative with p(a) + p(b) balanced, and I = p(b) - p(a) = int f.
    if let Ok(r) = levin(&|x| x * x, &|x| x, &|_| 1.0, 0.0, 1.0, 0.0, 12) {
        assert!(
            (r.re - 1.0 / 3.0).abs() < 1e-8 && r.im.abs() < 1e-8,
            "{r:?}"
        );
    }
    // tiny frequency: continuity with the plain integral, int_0^1 x^2 e^{i w x} ~ 1/3 + i w/4
    let w = 1e-4;
    let r = levin(&|x| x * x, &|x| x, &|_| 1.0, 0.0, 1.0, w, 12).unwrap();
    assert!((r.re - 1.0 / 3.0).abs() < 1e-7, "{r:?}");
    assert!((r.im - w / 4.0).abs() < 1e-7, "{r:?}");
}

// -------------------------------------------------------------------- SDE

#[test]
#[allow(deprecated)]
fn legacy_sde_adaptive_is_the_options_form() {
    // Geometric Brownian motion dy = mu y dt + sigma y dW (Ito).
    let (mu, sigma) = (0.4, 0.3);
    let drift = move |_t: f64, y: &[f64], out: &mut [f64]| out[0] = mu * y[0];
    let diff = move |_t: f64, y: &[f64], out: &mut [f64]| out[0] = sigma * y[0];
    let old = sde_adaptive(drift, diff, 1, 0.0, 1.0, &[1.0], 1e-3, 1e-6, 42, 100_000).unwrap();
    let new = sde_adaptive_with(
        drift,
        diff,
        1,
        0.0,
        1.0,
        &[1.0],
        &SdeAdaptiveOptions {
            rtol: 1e-3,
            atol: 1e-6,
            seed: 42,
            max_steps: 100_000,
        },
    )
    .unwrap();
    assert_eq!(old, new, "the legacy wrapper must be bit-identical");
    // path invariants
    assert_eq!(old.seed, 42);
    assert_eq!(old.t.len(), old.y.len());
    assert_eq!(old.accepted, old.t.len() - 1);
    assert_eq!(old.t[0], 0.0);
    assert!((old.t.last().unwrap() - 1.0).abs() < 1e-12);
    assert!(
        old.t.windows(2).all(|p| p[1] > p[0]),
        "time must be strictly increasing"
    );
    assert!(
        old.y
            .iter()
            .all(|s| s.len() == 1 && s[0].is_finite() && s[0] > 0.0),
        "GBM stays positive"
    );
    assert!(old.last_h > 0.0);
    // determinism and seed dependence
    let again = sde_adaptive(drift, diff, 1, 0.0, 1.0, &[1.0], 1e-3, 1e-6, 42, 100_000).unwrap();
    assert_eq!(old, again);
    let other = sde_adaptive(drift, diff, 1, 0.0, 1.0, &[1.0], 1e-3, 1e-6, 43, 100_000).unwrap();
    assert_ne!(old.y.last(), other.y.last());
}

#[test]
#[allow(deprecated)]
fn legacy_sde_adaptive_without_noise_follows_the_ode() {
    // sigma = 0: dy = -2 y dt is deterministic, y(1) = exp(-2).
    let drift = |_t: f64, y: &[f64], out: &mut [f64]| out[0] = -2.0 * y[0];
    let diff = |_t: f64, _y: &[f64], out: &mut [f64]| out[0] = 0.0;
    let p = sde_adaptive(drift, diff, 1, 0.0, 1.0, &[1.0], 1e-6, 1e-9, 1, 1_000_000).unwrap();
    let last = p.y.last().unwrap()[0];
    assert!((last - (-2.0_f64).exp()).abs() < 1e-3, "{last}");
    // the step count is bounded well below the budget for a smooth problem
    assert!(p.accepted < 100_000);
}

#[test]
#[allow(deprecated)]
fn legacy_sde_adaptive_rejects_bad_arguments() {
    let drift = |_t: f64, y: &[f64], out: &mut [f64]| out[0] = -y[0];
    let diff = |_t: f64, _y: &[f64], out: &mut [f64]| out[0] = 0.1;
    let ok = |t0: f64, t1: f64, y0: &[f64], m: usize, rtol: f64, atol: f64, steps: usize| {
        sde_adaptive(drift, diff, m, t0, t1, y0, rtol, atol, 1, steps)
    };
    assert!(matches!(
        ok(0.0, 1.0, &[], 1, 1e-3, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(0.0, 1.0, &[1.0], 0, 1e-3, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(1.0, 0.0, &[1.0], 1, 1e-3, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(1.0, 1.0, &[1.0], 1, 1e-3, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(0.0, 1.0, &[1.0], 1, 0.0, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(0.0, 1.0, &[1.0], 1, -1e-3, 1e-6, 100),
        Err(SdeError::Invalid { .. })
    ));
    assert!(matches!(
        ok(0.0, 1.0, &[1.0], 1, 1e-3, -1.0, 100),
        Err(SdeError::Invalid { .. })
    ));
    // a step budget that is too small is reported with its value
    assert_eq!(
        ok(0.0, 1.0, &[1.0], 1, 1e-8, 1e-10, 3).unwrap_err(),
        SdeError::StepLimit { limit: 3 }
    );
    // a blow-up (explosive drift) is detected as non-finite, not returned as data
    let boom = |_t: f64, y: &[f64], out: &mut [f64]| out[0] = y[0] * y[0] * y[0];
    let r = sde_adaptive(boom, diff, 1, 0.0, 10.0, &[10.0], 1e-3, 1e-6, 1, 1_000_000);
    assert!(r.is_err(), "finite-time blow-up must be an error: {r:?}");
}

// ---------------------------------------------------------- coulomb_g_from_u

#[test]
fn coulomb_irregular_function_free_particle_closed_forms() {
    for rho in [0.05, 0.5, 1.0, 3.0, 7.5, 15.0] {
        // eta = 0, l = 0:  G_0 = cos(rho), G_0' = -sin(rho)
        let (g, gp) = coulomb_g_from_u(0.0, 0, rho).unwrap();
        assert!((g - rho.cos()).abs() < 1e-11, "G0({rho})");
        assert!((gp + rho.sin()).abs() < 1e-11, "G0'({rho})");
        // eta = 0, l = 1:  G_1 = cos(rho)/rho + sin(rho)
        let (g, gp) = coulomb_g_from_u(0.0, 1, rho).unwrap();
        let g1 = rho.cos() / rho + rho.sin();
        let g1p = -rho.sin() / rho - rho.cos() / (rho * rho) + rho.cos();
        assert!(
            (g - g1).abs() < 1e-9 * (1.0 + g1.abs()),
            "G1({rho}): {g} vs {g1}"
        );
        assert!((gp - g1p).abs() < 1e-9 * (1.0 + g1p.abs()), "G1'({rho})");
        // eta = 0, l = 2: G_2 = (3/rho^2 - 1) cos(rho)/rho ... = -rho y_2(rho)
        //   y_2(x) = (-3/x^2 + 1) cos(x)/x - 3 sin(x)/x^2
        let y2 = (-3.0 / (rho * rho) + 1.0) * rho.cos() / rho - 3.0 * rho.sin() / (rho * rho);
        let (g, _) = coulomb_g_from_u(0.0, 2, rho).unwrap();
        assert!(
            (g + rho * y2).abs() < 1e-8 * (1.0 + g.abs()),
            "G2({rho}): {g} vs {}",
            -rho * y2
        );
    }
}

#[test]
fn coulomb_wronskian_and_agreement_with_the_other_route() {
    // F' G - F G' = 1 for every eta, l (DLMF 33.2.9) -- F comes from the series, G from U
    for &eta in &[-3.0, -1.0, 0.0, 0.5, 2.0, 4.0] {
        for l in [0usize, 1, 3, 6] {
            for &rho in &[0.3, 1.0, 4.0, 11.0] {
                let (g, gp) = coulomb_g_from_u(eta, l, rho).expect("in domain");
                let f = coulomb_f(eta, l, rho).unwrap();
                let fp = coulomb_f_prime(eta, l, rho).unwrap();
                let w = fp * g - f * gp;
                let scale = (fp * g).abs().max((f * gp).abs()).max(1.0);
                assert!(
                    (w - 1.0).abs() < 1e-8 * scale,
                    "eta={eta} l={l} rho={rho}: W = {w}"
                );
                // the production routine agrees with the Tricomi route
                let g2 = coulomb_g(eta, l, rho).unwrap();
                let gp2 = coulomb_g_prime(eta, l, rho).unwrap();
                assert!(
                    (g - g2).abs() < 1e-7 * (1.0 + g.abs()),
                    "G eta={eta} l={l} rho={rho}: {g} vs {g2}"
                );
                assert!(
                    (gp - gp2).abs() < 1e-7 * (1.0 + gp.abs()),
                    "G' eta={eta} l={l} rho={rho}"
                );
            }
        }
    }
}

#[test]
fn coulomb_g_from_u_small_rho_and_domain() {
    // near the origin G_l ~ rho^{-l} / ((2l+1) C_l(eta)): it blows up like rho^{-l}
    let (g_small, _) = coulomb_g_from_u(0.0, 1, 1e-3).unwrap();
    let (g_smaller, _) = coulomb_g_from_u(0.0, 1, 1e-4).unwrap();
    assert!(
        g_smaller / g_small > 8.0 && g_smaller / g_small < 12.0,
        "{g_small} {g_smaller}"
    );
    // outside the domain: None, never a panic or a NaN
    assert!(coulomb_g_from_u(0.0, 0, 0.0).is_none());
    assert!(coulomb_g_from_u(0.0, 0, -1.0).is_none());
    assert!(coulomb_g_from_u(0.0, 65, 1.0).is_none());
    assert!(coulomb_g_from_u(f64::NAN, 0, 1.0).is_none());
    assert!(coulomb_g_from_u(f64::INFINITY, 0, 1.0).is_none());
    assert!(coulomb_g_from_u(0.0, 0, f64::NAN).is_none());
    assert!(coulomb_g_from_u(0.0, 0, f64::INFINITY).is_none());
    // the largest allowed l either gives finite numbers or None
    if let Some((g, gp)) = coulomb_g_from_u(0.5, 64, 80.0) {
        assert!(g.is_finite() && gp.is_finite());
    }
}

// ------------------------------------------------------------ coordinates

#[test]
fn pure_conversions_round_trip_and_agree_with_checked_ones() {
    let cartesian_points: [[f64; 3]; 6] = [
        [1.0, 2.0, 3.0],
        [-1.5, 0.25, -2.0],
        [0.0, 0.0, 4.0],
        [3.0, -4.0, 0.0],
        [-0.3, -0.4, -0.5],
        [1e-8, -2e-8, 3e-8],
    ];
    for p in &cartesian_points {
        for sys in [C::Cylindrical, C::Spherical] {
            let there = from_cartesian_pure(p, sys);
            let back = to_cartesian_pure(&there, sys);
            for k in 0..3 {
                assert!(
                    (back[k] - p[k]).abs() < 1e-14 * (1.0 + p[k].abs()),
                    "{sys:?} {p:?} -> {there:?} -> {back:?}"
                );
            }
            assert_eq!(
                there,
                from_cartesian(p, sys).unwrap(),
                "pure and checked agree"
            );
            assert_eq!(
                to_cartesian_pure(&there, sys),
                to_cartesian(&there, sys).unwrap()
            );
        }
        let planar = [p[0], p[1]];
        let polar = from_cartesian_pure(&planar, C::Polar);
        let back = to_cartesian_pure(&polar, C::Polar);
        assert!((back[0] - p[0]).abs() < 1e-14 && (back[1] - p[1]).abs() < 1e-14);
        assert!((polar[0] - p[0].hypot(p[1])).abs() < 1e-15);
        assert!(
            polar[1] > -std::f64::consts::PI - 1e-15 && polar[1] <= std::f64::consts::PI + 1e-15
        );
    }
    // the pure transform is the composition through Cartesian coordinates
    let sph = [2.0, 0.7, -1.1];
    let via_pure = transform_point_pure(&sph, C::Spherical, C::Cylindrical);
    let via_checked = transform_point(&sph, C::Spherical, C::Cylindrical).unwrap();
    assert_eq!(via_pure, via_checked);
    // cylindrical radius is r sin(theta), height is r cos(theta), azimuth passes through
    assert!((via_pure[0] - 2.0 * 0.7_f64.sin()).abs() < 1e-14);
    assert!((via_pure[2] - 2.0 * 0.7_f64.cos()).abs() < 1e-14);
    assert!((via_pure[1] - -1.1).abs() < 1e-14);
}

#[test]
fn pure_conversions_known_values_and_degenerate_points() {
    use std::f64::consts::FRAC_PI_2;
    use std::f64::consts::PI;
    // (r, theta, phi) = (2, pi/2, pi/2) is the point (0, 2, 0)
    let p = to_cartesian_pure(&[2.0, FRAC_PI_2, FRAC_PI_2], C::Spherical);
    assert!(p[0].abs() < 1e-15 && (p[1] - 2.0).abs() < 1e-15 && p[2].abs() < 1e-15);
    // poles: theta = 0 on +z, pi on -z; the origin has theta = 0 by convention
    assert_eq!(
        from_cartesian_pure(&[0.0, 0.0, 5.0], C::Spherical)[..2],
        [5.0, 0.0]
    );
    let south = from_cartesian_pure(&[0.0, 0.0, -5.0], C::Spherical);
    assert!((south[0] - 5.0).abs() < 1e-15 && (south[1] - PI).abs() < 1e-15);
    assert_eq!(
        from_cartesian_pure(&[0.0, 0.0, 0.0], C::Spherical),
        vec![0.0, 0.0, 0.0]
    );
    assert_eq!(from_cartesian_pure(&[0.0, 0.0], C::Polar), vec![0.0, 0.0]);
    // atan2 branch: negative x-axis has azimuth +pi, negative y has -pi/2
    assert!((from_cartesian_pure(&[-1.0, 0.0], C::Polar)[1] - PI).abs() < 1e-15);
    assert!((from_cartesian_pure(&[0.0, -1.0], C::Polar)[1] + FRAC_PI_2).abs() < 1e-15);
    // missing components are zero, extra components are ignored (non-Cartesian targets)
    assert_eq!(to_cartesian_pure(&[2.0], C::Polar), vec![2.0, 0.0]);
    assert_eq!(to_cartesian_pure(&[], C::Polar), vec![0.0, 0.0]);
    assert_eq!(to_cartesian_pure(&[], C::Spherical), vec![0.0, 0.0, 0.0]);
    assert_eq!(from_cartesian_pure(&[3.0], C::Polar), vec![3.0, 0.0]);
    assert_eq!(
        from_cartesian_pure(&[3.0, 4.0, 99.0], C::Polar),
        vec![5.0, 4.0_f64.atan2(3.0)]
    );
    assert_eq!(
        to_cartesian_pure(&[1.0, 0.0, 7.0, 8.0], C::Cylindrical),
        vec![1.0, 0.0, 7.0]
    );
    // Cartesian is the identity for any length
    for v in [vec![], vec![1.0], vec![1.0, 2.0, 3.0, 4.0]] {
        assert_eq!(from_cartesian_pure(&v, C::Cartesian), v);
        assert_eq!(to_cartesian_pure(&v, C::Cartesian), v);
    }
    // the checked versions reject what the pure ones pad
    assert!(to_cartesian(&[2.0], C::Polar).is_err());
    assert!(to_cartesian(&[1.0, 2.0, 3.0], C::Polar).is_err());
    assert!(from_cartesian(&[1.0, 2.0], C::Spherical).is_err());
    assert!(transform_point(&[1.0, 2.0], C::Polar, C::Spherical).is_err());
    // NaN / inf propagate without a panic
    assert!(
        from_cartesian_pure(&[f64::NAN, 1.0, 1.0], C::Spherical)
            .iter()
            .any(|v| v.is_nan())
    );
    assert!(to_cartesian_pure(&[f64::INFINITY, 0.3], C::Polar).len() == 2);
    // a point just beyond the unit-cosine rounding limit does not produce NaN in theta
    let s = from_cartesian_pure(&[1e-200, 1e-200, 1.0], C::Spherical);
    assert!(s[1].is_finite() && s[1] >= 0.0);
}

// --------------------------------------------------------------- wavelets

#[test]
fn wavelet_filter_banks_have_the_documented_normalisation() {
    let sqrt2 = std::f64::consts::SQRT_2;
    for w in ALL_WAVELETS {
        let (lo, hi) = (w.low_pass(), w.high_pass());
        assert_eq!(lo.len(), w.taps(), "{w:?}");
        assert_eq!(hi.len(), lo.len(), "{w:?}: equal-length banks");
        // sum h~ = sqrt 2 (DC gain), the high-pass annihilates constants
        assert!(
            (lo.iter().sum::<f64>() - sqrt2).abs() < 1e-9,
            "{w:?}: sum = {}",
            lo.iter().sum::<f64>()
        );
        assert!(hi.iter().sum::<f64>().abs() < 1e-9, "{w:?}");
        assert_eq!(hi, w.analysis_high());
        assert_eq!(lo, w.analysis_low());
        // the high-pass has at least the advertised vanishing moments: sum_k k^p g[k] = 0
        let (nr, _) = w.vanishing_moments();
        for p in 0..nr.min(4) {
            let moment: f64 = hi
                .iter()
                .enumerate()
                .map(|(k, g)| (k as f64).powi(p as i32) * g)
                .sum();
            let scale: f64 = hi
                .iter()
                .enumerate()
                .map(|(k, g)| (k as f64).powi(p as i32) * g.abs())
                .sum();
            assert!(
                moment.abs() < 1e-8 * scale.max(1.0),
                "{w:?}: moment {p} = {moment}"
            );
        }
        // Nyquist: the low-pass vanishes at z = -1 (sum (-1)^k h[k] = 0), for the same reason
        let nyquist: f64 = lo
            .iter()
            .enumerate()
            .map(|(k, h)| if k % 2 == 0 { *h } else { -*h })
            .sum();
        assert!(nyquist.abs() < 1e-8, "{w:?}: Nyquist gain {nyquist}");
    }
}

#[test]
fn orthogonal_wavelet_filters_are_double_shift_orthonormal() {
    let orthogonal = [
        Wavelet::Haar,
        Wavelet::Db2,
        Wavelet::Db3,
        Wavelet::Db4,
        Wavelet::Db5,
        Wavelet::Db6,
        Wavelet::Db7,
        Wavelet::Db8,
        Wavelet::Db9,
        Wavelet::Db10,
        Wavelet::Sym4,
        Wavelet::Sym5,
        Wavelet::Sym6,
    ];
    for w in orthogonal {
        let h = w.low_pass();
        let g = w.high_pass();
        // synthesis equals analysis for orthogonal families
        assert_eq!(w.synthesis_low(), h, "{w:?}");
        let n = h.len();
        for m in 0..(n / 2) {
            let hh: f64 = (0..n)
                .filter(|k| k + 2 * m < n)
                .map(|k| h[k] * h[k + 2 * m])
                .sum();
            let gg: f64 = (0..n)
                .filter(|k| k + 2 * m < n)
                .map(|k| g[k] * g[k + 2 * m])
                .sum();
            let expect = if m == 0 { 1.0 } else { 0.0 };
            assert!((hh - expect).abs() < 1e-9, "{w:?}: <h, h shift {m}> = {hh}");
            assert!((gg - expect).abs() < 1e-9, "{w:?}: <g, g shift {m}> = {gg}");
        }
        // low-pass and high-pass are orthogonal at every even shift
        for m in -(n as i64 / 2)..=(n as i64 / 2) {
            let hg: f64 = (0..n as i64)
                .filter_map(|k| {
                    let j = k + 2 * m;
                    (j >= 0 && j < n as i64).then(|| h[k as usize] * g[j as usize])
                })
                .sum();
            assert!(hg.abs() < 1e-9, "{w:?}: <h, g shift {m}> = {hg}");
        }
        // g[n] = (-1)^n h[L - 1 - n]  (the documented alternating mirror of the synthesis filter)
        for k in 0..n {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            assert!((g[k] - sign * h[n - 1 - k]).abs() < 1e-15, "{w:?} tap {k}");
        }
    }
    // closed forms for the shortest filters
    let r = std::f64::consts::FRAC_1_SQRT_2;
    assert_eq!(Wavelet::Haar.low_pass(), vec![r, r]);
    assert_eq!(Wavelet::Haar.high_pass(), vec![r, -r]);
    let s3 = 3.0_f64.sqrt();
    let db2 = Wavelet::Db2.low_pass();
    let d = 4.0 * std::f64::consts::SQRT_2;
    for (got, want) in db2.iter().zip([
        (1.0 + s3) / d,
        (3.0 + s3) / d,
        (3.0 - s3) / d,
        (1.0 - s3) / d,
    ]) {
        assert!((got - want).abs() < 1e-15);
    }
}
