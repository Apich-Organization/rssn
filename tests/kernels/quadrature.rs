//! Tests for the advanced quadrature kernels.

use std::f64::consts::PI;
use std::time::Instant;

use num_complex::Complex64;
use rssn::kernels::quadrature::GkRule;
use rssn::kernels::quadrature::LowDiscrepancy;
use rssn::kernels::quadrature::Sobol;
use rssn::kernels::quadrature::clenshaw_curtis;
use rssn::kernels::quadrature::clenshaw_curtis_rule;
use rssn::kernels::quadrature::cubature_genz_malik;
use rssn::kernels::quadrature::cubature_nested;
use rssn::kernels::quadrature::exp_sinh;
use rssn::kernels::quadrature::filon;
use rssn::kernels::quadrature::filon2d;
use rssn::kernels::quadrature::filon2d_complex;
use rssn::kernels::quadrature::halton_point;
use rssn::kernels::quadrature::integrate;
use rssn::kernels::quadrature::integrate_gk;
use rssn::kernels::quadrature::integrate_infinite;
use rssn::kernels::quadrature::integrate_semi_infinite;
use rssn::kernels::quadrature::monte_carlo;
use rssn::kernels::quadrature::oscillatory2d;
use rssn::kernels::quadrature::quasi_monte_carlo;
use rssn::kernels::quadrature::romberg;
use rssn::kernels::quadrature::sinh_sinh;
use rssn::kernels::quadrature::tanh_sinh;

#[test]
fn gk_smooth_both_rules() {
    let exact = 1.0_f64.exp() - 1.0;
    for rule in [GkRule::G7K15, GkRule::G10K21] {
        let r = integrate_gk(f64::exp, 0.0, 1.0, rule, 1e-14, 1e-14, 100);
        assert!((r.value - exact).abs() < 1e-14);
        assert!(r.converged);
    }
}

#[test]
fn gk_polynomial_exact() {
    let r = integrate(|x| x.powi(13), 0.0, 1.0, 1e-12);
    assert!((r.value - 1.0 / 14.0).abs() < 1e-14);
}

#[test]
fn gk_endpoint_singularity() {
    let r = integrate_gk(
        |x| 1.0 / x.sqrt(),
        0.0,
        1.0,
        GkRule::G7K15,
        1e-10,
        1e-10,
        2000,
    );
    assert!((r.value - 2.0).abs() < 1e-7, "{r:?}");
    let r = integrate_gk(f64::ln, 0.0, 1.0, GkRule::G10K21, 1e-10, 1e-10, 2000);
    assert!((r.value + 1.0).abs() < 1e-8, "{r:?}");
}

#[test]
fn gk_infinite_ranges() {
    let r = integrate_semi_infinite(|x| (-x).exp(), 0.0, 1e-10);
    assert!((r.value - 1.0).abs() < 1e-9);
    let r = integrate_infinite(|x| (-x * x).exp(), 1e-10);
    assert!((r.value - PI.sqrt()).abs() < 1e-9);
}

#[test]
fn tanh_sinh_singular_and_smooth() {
    let r = tanh_sinh(|x| 1.0 / x.sqrt(), 0.0, 1.0, 1e-12);
    assert!((r.value - 2.0).abs() < 1e-11, "{r:?}");
    let r = tanh_sinh(f64::ln, 0.0, 1.0, 1e-12);
    assert!((r.value + 1.0).abs() < 1e-11, "{r:?}");
    let r = tanh_sinh(f64::cos, 0.0, PI / 2.0, 1e-12);
    assert!((r.value - 1.0).abs() < 1e-12);
    let r = tanh_sinh(|x| x, 1.0, 0.0, 1e-12);
    assert!((r.value + 0.5).abs() < 1e-12);
}

#[test]
fn exp_sinh_and_sinh_sinh() {
    let r = exp_sinh(|x| (-x).exp(), 0.0, 1e-12);
    assert!((r.value - 1.0).abs() < 1e-10, "{r:?}");
    let r = exp_sinh(|x| 1.0 / (1.0 + x * x), 0.0, 1e-12);
    assert!((r.value - PI / 2.0).abs() < 1e-9, "{r:?}");
    let r = sinh_sinh(|x| (-x * x).exp(), 1e-12);
    assert!((r.value - PI.sqrt()).abs() < 1e-11, "{r:?}");
}

#[test]
fn clenshaw_curtis_exactness() {
    let v = clenshaw_curtis(|x| x.powi(8) + 1.0, -1.0, 1.0, 8);
    assert!((v - (2.0 / 9.0 + 2.0)).abs() < 1e-14);
    let v = clenshaw_curtis(f64::exp, 0.0, 1.0, 24);
    assert!((v - (1.0_f64.exp() - 1.0)).abs() < 1e-14);
}

#[test]
fn romberg_reference() {
    let r = romberg(f64::sin, 0.0, PI, 1e-12, 14);
    assert!((r.value - 2.0).abs() < 1e-11);
    assert!(r.converged);
}

#[test]
fn filon_oscillatory() {
    let w: f64 = 50.0;
    let exact = w.sin() / w + (w.cos() - 1.0) / (w * w);
    let v = filon(|x| x, 0.0, 1.0, w, true, 10);
    assert!((v - exact).abs() < 1e-12, "{v} vs {exact}");
    let w: f64 = 100.0;
    let exact = w * (1.0 - PI.exp() * (w * PI).cos()) / (1.0 + w * w);
    let v = filon(f64::exp, 0.0, PI, w, false, 400);
    assert!((v - exact).abs() < 1e-8, "{v} vs {exact}");
    let v = filon(|x| x * x, 0.0, 1.0, 0.01, true, 4);
    assert!((v - 1.0 / 3.0).abs() < 1e-4);
}

#[test]
fn low_discrepancy_points() {
    assert_eq!(halton_point(1, 2), vec![0.5, 1.0 / 3.0]);
    assert_eq!(halton_point(2, 1), vec![0.25]);
    let mut s = Sobol::new(2);
    let p1 = s.next_point();
    let p2 = s.next_point();
    assert_eq!(p1, vec![0.5, 0.5]);
    assert_eq!(p2, vec![0.75, 0.25]);
}

#[test]
fn qmc_and_mc() {
    let f = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>();
    let lo = [0.0; 4];
    let hi = [1.0; 4];
    let exact = 4.0 / 3.0;
    let s = quasi_monte_carlo(f, &lo, &hi, 4096, LowDiscrepancy::Sobol);
    let h = quasi_monte_carlo(f, &lo, &hi, 4096, LowDiscrepancy::Halton);
    let m = monte_carlo(f, &lo, &hi, 4096, 42);
    assert!((s.value - exact).abs() < 1e-3, "{s:?}");
    assert!((h.value - exact).abs() < 1e-2, "{h:?}");
    assert!((m.value - exact).abs() < 5.0 * m.error + 1e-12);
    let m2 = monte_carlo(f, &lo, &hi, 4096, 42);
    assert_eq!(m.value.to_bits(), m2.value.to_bits());
}

#[test]
fn nested_cubature() {
    let v = cubature_nested(|x| (x[0] * x[1]).exp(), &[0.0, 0.0], &[1.0, 1.0], 1e-10);
    assert!((v.value - 1.317_902_151_454_404).abs() < 1e-8, "{v:?}");
    let v = cubature_nested(|x| x[0] + x[1] + x[2], &[0.0; 3], &[1.0; 3], 1e-10);
    assert!((v.value - 1.5).abs() < 1e-9);
}

#[test]
fn genz_malik_smooth_2_to_5d() {
    // 2D: int exp(x+y) = (e-1)^2
    let r = cubature_genz_malik(
        |x| (x[0] + x[1]).exp(),
        &[0.0; 2],
        &[1.0; 2],
        1e-12,
        1e-12,
        200_000,
    );
    assert!(r.converged);
    assert!(
        (r.value - (1.0_f64.exp() - 1.0).powi(2)).abs() < 1e-10,
        "{r:?}"
    );
    // 3D polynomial of degree <= 7 is integrated exactly: int x^2 y^3 z over [0,2]x[0,1]x[-1,1] = 8/3 * 1/4 * 0
    let r = cubature_genz_malik(
        |x| x[0].powi(4) * x[1].powi(3) + x[2].powi(6),
        &[0.0, 0.0, -1.0],
        &[2.0, 1.0, 1.0],
        1e-9,
        1e-9,
        100_000,
    );
    let exact = 3.2 + 4.0 / 7.0;
    assert!((r.value - exact).abs() < 1e-11 * exact, "{r:?} vs {exact}");
    // 4D: int prod cos(x_i) = sin(1)^4
    let r = cubature_genz_malik(
        |x| x.iter().map(|v| v.cos()).product(),
        &[0.0; 4],
        &[1.0; 4],
        1e-11,
        1e-11,
        2_000_000,
    );
    assert!((r.value - 1.0_f64.sin().powi(4)).abs() < 1e-9, "{r:?}");
    // 5D
    let r = cubature_genz_malik(
        |x| x.iter().map(|v| v.cos()).product(),
        &[0.0; 5],
        &[1.0; 5],
        1e-9,
        1e-9,
        5_000_000,
    );
    assert!(r.converged, "{r:?}");
    assert!((r.value - 1.0_f64.sin().powi(5)).abs() < 1e-8, "{r:?}");
}

#[test]
fn genz_malik_peaked_integrands() {
    // Gaussian peak of width ~0.07 at the centre of the unit cube
    let a = 200.0;
    let exact1 = (PI / a).sqrt() * rssn_erf(0.5 * a.sqrt());
    let r = cubature_genz_malik(
        |x| (-a * x.iter().map(|v| (v - 0.5).powi(2)).sum::<f64>()).exp(),
        &[0.0; 3],
        &[1.0; 3],
        1e-12,
        1e-7,
        5_000_000,
    );
    assert!(r.converged, "{r:?}");
    assert!(
        (r.value - exact1.powi(3)).abs() < 1e-6 * exact1.powi(3),
        "{r:?} vs {}",
        exact1.powi(3)
    );
    // Genz corner peak: (1 + sum x)^-(n+1) over [0,1]^n = (1/n!) sum (-1)^k C(n,k)/(1+k)
    for n in 2..=4usize {
        let mut fact = 1.0;
        let mut binom = 1.0;
        let mut s = 0.0;
        for k in 0..=n {
            if k > 0 {
                binom = binom * (n - k + 1) as f64 / k as f64;
            }
            s += if k % 2 == 0 { 1.0 } else { -1.0 } * binom / (1.0 + k as f64);
        }
        for i in 1..=n {
            fact *= i as f64;
        }
        let exact = s / fact;
        let r = cubature_genz_malik(
            |x| (1.0 + x.iter().sum::<f64>()).powi(-(n as i32 + 1)),
            &vec![0.0; n],
            &vec![1.0; n],
            1e-13,
            1e-10,
            3_000_000,
        );
        assert!(
            (r.value - exact).abs() < 1e-9 * exact,
            "n={n}: {r:?} vs {exact}"
        );
    }
    // off-centre narrow peak in 2D
    let r = cubature_genz_malik(
        |x| (-1000.0 * ((x[0] - 0.3).powi(2) + (x[1] - 0.7).powi(2))).exp(),
        &[0.0; 2],
        &[1.0; 2],
        1e-12,
        1e-9,
        2_000_000,
    );
    let e = (PI / 1000.0).sqrt()
        * 0.5
        * (rssn_erf(0.3 * 1000f64.sqrt()) + rssn_erf(0.7 * 1000f64.sqrt()));
    let e2 = (PI / 1000.0).sqrt()
        * 0.5
        * (rssn_erf(0.7 * 1000f64.sqrt()) + rssn_erf(0.3 * 1000f64.sqrt()));
    assert!((r.value - e * e2).abs() < 1e-8, "{r:?} vs {}", e * e2);
    assert!(
        cubature_genz_malik(|_| 1.0, &[], &[], 1e-6, 1e-6, 100)
            .value
            .is_nan()
    );
}

/// erf via the Gauss-Kronrod rule (reference helper for the tests).
fn rssn_erf(x: f64) -> f64 {
    2.0 / PI.sqrt() * integrate(|t| (-t * t).exp(), 0.0, x, 1e-14).value
}

// -------------------------------------------------------------------------
// Two-dimensional oscillatory (Filon) quadrature
// -------------------------------------------------------------------------

/// Closed form of `int_0^1 int_0^1 cos(w x) cos(k y) dx dy`.
fn coscos_exact(
    w: f64,
    k: f64,
) -> f64 {
    let fx = if w == 0.0 { 1.0 } else { w.sin() / w };
    let fy = if k == 0.0 { 1.0 } else { k.sin() / k };
    fx * fy
}

/// `int_0^1 int_0^1 cos(w x) cos(k y) dx dy` as two linear-phase Filon calls
/// through the product-to-sum identity.
fn coscos_filon(
    w: f64,
    k: f64,
    tol: f64,
) -> (f64, f64, usize) {
    let (lo, hi) = ([0.0, 0.0], [1.0, 1.0]);
    if w == 0.0 {
        let r = filon2d(|_, _| 1.0, &lo, &hi, k, &[0.0, 1.0], true, tol);
        return (r.value, r.error, r.evaluations);
    }
    let a = filon2d(|_, _| 1.0, &lo, &hi, w, &[1.0, k / w], true, tol);
    let b = filon2d(|_, _| 1.0, &lo, &hi, w, &[1.0, -k / w], true, tol);
    (
        a.value.midpoint(b.value),
        a.error.midpoint(b.error),
        a.evaluations + b.evaluations,
    )
}

/// Closed form of `int_0^1 exp(i w x) dx`.
fn complex_1d(w: f64) -> Complex64 {
    if w == 0.0 {
        Complex64::new(1.0, 0.0)
    } else {
        (Complex64::new(0.0, w).exp() - 1.0) / Complex64::new(0.0, w)
    }
}

#[test]
fn filon2d_cos_cos_matches_the_closed_form() {
    let cases = [
        (0.0, 0.0),
        (0.0, 5.0),
        (5.0, 0.0),
        (1.0, 1.0),
        (5.0, 3.0),
        (50.0, 20.0),
        (1000.0, 1000.0),
        (1000.0, 0.0),
    ];
    for (w, k) in cases {
        let exact = coscos_exact(w, k);
        let (value, estimate, evaluations) = coscos_filon(w, k, 1e-12);
        let error = (value - exact).abs();
        eprintln!(
            "filon2d coscos w={w}, k={k}: value {value:.15e}, exact {exact:.15e}, \
             error {error:.3e}, estimate {estimate:.3e}, evaluations {evaluations}"
        );
        // The rule is exact for a constant amplitude: only round-off remains,
        // even when the frequency is 1000 and only a few dozen nodes are used.
        assert!(
            error <= 1e-10 * exact.abs() + 1e-13,
            "w={w} k={k}: error {error}"
        );
        assert!(
            estimate <= 1e-9 * exact.abs().max(1.0) + 1e-13,
            "w={w} k={k}: estimate {estimate}"
        );
    }
}

#[test]
fn filon2d_omega_zero_degrades_to_the_classical_rule() {
    // int_0^1 int_0^1 exp(x + y) dx dy = (e - 1)^2.
    let exact = (1.0_f64.exp() - 1.0).powi(2);
    let r = filon2d(
        |x, y| (x + y).exp(),
        &[0.0, 0.0],
        &[1.0, 1.0],
        0.0,
        &[1.0, 1.0],
        true,
        1e-9,
    );
    eprintln!(
        "filon2d omega = 0: value {:.15e}, exact {exact:.15e}, error {:.3e}, estimate {:.3e}, \
         evaluations {}",
        r.value,
        (r.value - exact).abs(),
        r.error,
        r.evaluations
    );
    assert!(r.converged, "{r:?}");
    assert!((r.value - exact).abs() <= 1e-9 * exact, "{r:?}");
    // With omega = 0 the oscillator is one, so the sine component vanishes.
    let s = filon2d(
        |x, y| (x + y).exp(),
        &[0.0, 0.0],
        &[1.0, 1.0],
        0.0,
        &[1.0, 1.0],
        false,
        1e-9,
    );
    assert!(s.value.abs() <= 1e-12, "sine component {}", s.value);
}

#[test]
fn filon2d_complex_exponential_matches_the_analytic_value() {
    let cases = [(0.0, 0.0), (0.0, 7.0), (3.0, 5.0), (500.0, 300.0)];
    for (w, k) in cases {
        let exact = complex_1d(w) * complex_1d(k);
        let r = if w == 0.0 {
            filon2d_complex(|_, _| 1.0, &[0.0, 0.0], &[1.0, 1.0], k, &[0.0, 1.0], 1e-12)
        } else {
            filon2d_complex(
                |_, _| 1.0,
                &[0.0, 0.0],
                &[1.0, 1.0],
                w,
                &[1.0, k / w],
                1e-12,
            )
        };
        let error = (r.value - exact).norm();
        eprintln!(
            "filon2d complex w={w}, k={k}: |error| {error:.3e}, estimate {:.3e}, evaluations {}",
            r.error, r.evaluations
        );
        assert!(r.converged, "{r:?}");
        assert!(
            error <= 1e-10 * exact.norm() + 1e-13,
            "w={w} k={k}: |error| {error}"
        );
    }
}

#[test]
fn filon2d_nonconstant_amplitude_matches_references() {
    let (lo, hi) = ([0.0, 0.0], [1.0, 1.0]);
    // A bilinear amplitude is in the tensor quadratic space, so omega = 0 is
    // exact.
    let r = filon2d(|x, y| x * y, &lo, &hi, 0.0, &[1.0, 1.0], true, 1e-12);
    assert!((r.value - 0.25).abs() < 1e-14, "{r:?}");
    // A cubic amplitude is resolved by panel refinement; compare with the
    // adaptive Genz-Malik cubature of the same oscillatory integrand.
    let phase = |x: f64, y: f64| 2.0 * (x + 0.5 * y);
    let reference = cubature_genz_malik(
        |p| p[0].powi(3) * p[1] * phase(p[0], p[1]).cos(),
        &lo,
        &hi,
        1e-13,
        1e-13,
        500_000,
    );
    assert!(reference.converged, "{reference:?}");
    let r = filon2d(
        |x, y| x.powi(3) * y,
        &lo,
        &hi,
        2.0,
        &[1.0, 0.5],
        true,
        1e-10,
    );
    let error = (r.value - reference.value).abs();
    eprintln!(
        "filon2d cubic amplitude: value {:.15e}, error {error:.3e}, estimate {:.3e}, \
         evaluations {}, reference {:.15e}",
        r.value, r.error, r.evaluations, reference.value
    );
    assert!(r.converged, "{r:?}");
    assert!(error <= 1e-8, "error {error}");
    // The error estimate must be a sane magnitude for the actual error.
    assert!(
        error <= 10.0 * r.error + 1e-12,
        "error {error} vs estimate {}",
        r.error
    );
}

#[test]
fn filon2d_beats_a_plain_tensor_rule_at_large_frequency() {
    let w = 1000.0_f64;
    let exact = coscos_exact(w, w);
    let start = Instant::now();
    let (value, estimate, evaluations) = coscos_filon(w, w, 1e-12);
    let elapsed = start.elapsed();
    let error = (value - exact).abs();
    eprintln!(
        "filon2d w = 1000: value {value:.15e}, exact {exact:.15e}, error {error:.3e}, \
         estimate {estimate:.3e}, evaluations {evaluations}, {elapsed:?}"
    );
    assert!(error <= 1e-10 * exact.abs() + 1e-15, "error {error}");
    assert!(
        evaluations <= 200,
        "the Filon rule should need only a few dozen nodes, used {evaluations}"
    );
    // For reference, the same integral by a plain tensor Clenshaw-Curtis rule:
    // 257^2 = 66,049 nodes still do not resolve 1000 oscillations.
    let (nodes, weights) = clenshaw_curtis_rule(256);
    let mut tensor = 0.0;
    for (i, &u) in nodes.iter().enumerate() {
        let x = 0.5 + 0.5 * u;
        for (j, &v) in nodes.iter().enumerate() {
            let y = 0.5 + 0.5 * v;
            tensor += 0.25 * weights[i] * weights[j] * (w * x).cos() * (w * y).cos();
        }
    }
    let tensor_error = (tensor - exact).abs();
    // One-dimensional Clenshaw-Curtis reaches 1e-10 only near 513 nodes, so the
    // tensor rule would need about 513^2 = 263,169 nodes to match.
    let one_d = clenshaw_curtis(|x: f64| (w * x).cos(), 0.0, 1.0, 512);
    let one_d_error = (one_d - w.sin() / w).abs();
    eprintln!(
        "plain tensor CC 257^2 error {tensor_error:.3e} (66,049 nodes); \
         1-D CC 513 nodes error {one_d_error:.3e}"
    );
    assert!(tensor_error > 1e-5, "tensor error {tensor_error}");
    assert!(one_d_error > 1e-11, "1-D CC n = 512 error {one_d_error}");
}

#[test]
fn filon2d_reduces_to_the_one_dimensional_rule() {
    // f(x, y) = g(x): the y integral of the oscillator with zero slope is the
    // interval length, so the 2-D value must agree with the 1-D Filon rule.
    let g = |x: f64| (0.5 * x).exp() * (3.0 * x).sin();
    for w in [0.0, 0.5, 7.0] {
        let reference = filon(g, 0.0, 1.0, w, true, 512);
        let two = filon2d(
            |x, _y| g(x),
            &[0.0, 0.0],
            &[1.0, 1.0],
            w,
            &[1.0, 0.0],
            true,
            1e-9,
        );
        eprintln!(
            "separable w={w}: 1-D {reference:.15e}, 2-D {:.15e}, evaluations {}",
            two.value, two.evaluations
        );
        assert!(two.converged, "{two:?}");
        assert!(
            (two.value - reference).abs() <= 1e-8 * reference.abs().max(1e-6) + 1e-13,
            "w = {w}: {} vs {reference}",
            two.value
        );
    }
    // A phase with zero slope both ways is the classical integral of f.
    let r = filon2d(
        |_, _| 1.0,
        &[0.0, 0.0],
        &[1.0, 1.0],
        1e3,
        &[0.0, 0.0],
        true,
        1e-12,
    );
    assert!((r.value - 1.0).abs() < 1e-14, "{r:?}");
}

#[test]
fn oscillatory2d_general_phase_fallback_matches_a_reference() {
    // phi = x^2 + y^2 is nonlinear: compare with a very fine tensor
    // Clenshaw-Curtis reference for exp(i w x^2) squared.
    let w = 8.0_f64;
    let (nodes, weights) = clenshaw_curtis_rule(4000);
    let mut one_d = Complex64::new(0.0, 0.0);
    for (i, &u) in nodes.iter().enumerate() {
        let x = 0.5 + 0.5 * u;
        one_d += Complex64::new(0.5 * weights[i], 0.0) * Complex64::new(0.0, w * x * x).exp();
    }
    let reference = one_d * one_d;
    let r = oscillatory2d(
        |_, _| 1.0,
        |x, y| x * x + y * y,
        &[0.0, 0.0],
        &[1.0, 1.0],
        w,
        1e-10,
    );
    let error = (r.value - reference).norm();
    eprintln!(
        "oscillatory2d general phase: value {r:?}, reference {reference:?}, \
         |error| {error:.3e}, estimate {:.3e}, evaluations {}",
        r.error, r.evaluations
    );
    assert!(r.converged, "{r:?}");
    assert!(error <= 1e-6 * reference.norm(), "error {error}");
}

#[test]
fn filon2d_invalid_arguments_are_reported_not_panicking() {
    let f = |_: f64, _: f64| 1.0;
    let (lo, hi) = ([0.0, 0.0], [1.0, 1.0]);
    // Wrong dimension lists.
    assert!(
        filon2d(f, &[0.0], &[1.0, 1.0], 1.0, &[1.0, 1.0], true, 1e-9)
            .value
            .is_nan()
    );
    assert!(filon2d(f, &lo, &hi, 1.0, &[1.0], true, 1e-9).value.is_nan());
    // Inverted box.
    assert!(
        filon2d(f, &[1.0, 0.0], &[0.0, 1.0], 1.0, &[1.0, 1.0], true, 1e-9)
            .value
            .is_nan()
    );
    // Non-finite frequency, slope or tolerance, and a non-positive tolerance.
    assert!(
        filon2d(f, &lo, &hi, f64::NAN, &[1.0, 1.0], true, 1e-9)
            .value
            .is_nan()
    );
    assert!(
        filon2d(f, &lo, &hi, 1.0, &[f64::INFINITY, 1.0], true, 1e-9)
            .value
            .is_nan()
    );
    assert!(
        filon2d(f, &lo, &hi, 1.0, &[1.0, 1.0], true, 0.0)
            .value
            .is_nan()
    );
    assert!(
        filon2d_complex(f, &lo, &hi, 1.0, &[1.0, 1.0], -1.0)
            .value
            .re
            .is_nan()
    );
    // Non-finite phase in the general fallback.
    let bad_phase = |_: f64, _: f64| f64::NAN;
    assert!(
        oscillatory2d(f, bad_phase, &lo, &hi, 1.0, 1e-9)
            .value
            .re
            .is_nan()
    );
    // A valid call reports finite values and a finite error estimate.
    let r = filon2d(f, &lo, &hi, 0.0, &[1.0, 1.0], true, 1e-10);
    assert!(r.converged && r.value.is_finite() && r.error.is_finite());
}
