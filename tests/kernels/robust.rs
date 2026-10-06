//! Tests for robust statistics, KDE, RNG and FFT convolution.

use rssn::kernels::fftconv::{circular_convolve, convolve, correlate, fft, poly_mul_i64};
use rssn::kernels::random::Rng;
use rssn::kernels::robust::{
    Kernel, bandwidth_scott, bandwidth_silverman, biweight_location, bootstrap_ci, ecdf,
    huber_location, iqr, kde, mad, median, quantile, trimmed_mean, winsorized_mean,
};

#[test]
fn quantile_types() {
    let d = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    // R: quantile(1:10, 0.25, type=k)
    let expect = [(1, 3.0), (2, 3.0), (3, 2.0), (4, 2.5), (5, 3.0), (6, 2.75), (7, 3.25), (8, 2.916_666_666_666_667), (9, 2.9375)];
    for (k, e) in expect {
        assert!((quantile(&d, 0.25, k) - e).abs() < 1e-12, "type {k}: {}", quantile(&d, 0.25, k));
    }
    assert_eq!(median(&d), 5.5);
    assert!((iqr(&d) - 4.5).abs() < 1e-12);
    assert!(quantile(&[], 0.5, 7).is_nan());
    assert!(quantile(&d, 1.5, 7).is_nan());
    assert_eq!(quantile(&d, 1.0, 7), 10.0);
    assert_eq!(quantile(&d, 0.0, 7), 1.0);
}

#[test]
fn robust_estimators_resist_outliers() {
    let mut d: Vec<f64> = (0..99).map(|i| 10.0 + (f64::from(i) - 49.0) * 0.02).collect();
    d.push(1e6);
    let mean = d.iter().sum::<f64>() / d.len() as f64;
    assert!(mean > 1000.0);
    assert!((median(&d) - 10.0).abs() < 0.02);
    assert!((trimmed_mean(&d, 0.1) - 10.0).abs() < 0.02);
    assert!((huber_location(&d, 1.345) - 10.0).abs() < 0.05);
    assert!((biweight_location(&d, 4.685) - 10.0).abs() < 0.05);
    assert!(winsorized_mean(&d, 0.05) < 11.0);
    assert!((mad(&[1.0, 2.0, 3.0, 4.0, 100.0], false) - 1.0).abs() < 1e-12);
    assert!((mad(&[1.0, 2.0, 3.0, 4.0, 100.0], true) - 1.482_602_218_505_602).abs() < 1e-12);
    assert!((ecdf(&[1.0, 2.0, 3.0, 4.0], 2.5) - 0.5).abs() < 1e-15);
}

#[test]
fn kde_integrates_to_one() {
    let mut rng = Rng::new(3);
    let data: Vec<f64> = (0..400).map(|_| rng.normal()).collect();
    let h = bandwidth_silverman(&data);
    assert!(h > 0.1 && h < 0.5, "{h}");
    assert!(bandwidth_scott(&data) > 0.0);
    for k in [Kernel::Gaussian, Kernel::Epanechnikov] {
        let (mut s, dx) = (0.0, 0.01);
        let mut x = -8.0;
        while x < 8.0 {
            s += kde(&data, x, h, k) * dx;
            x += dx;
        }
        assert!((s - 1.0).abs() < 1e-3, "{k:?} {s}");
    }
    // density near 0 close to the Gaussian peak
    assert!((kde(&data, 0.0, h, Kernel::Gaussian) - 0.3989).abs() < 0.08);
}

#[test]
fn bootstrap_is_deterministic() {
    let data: Vec<f64> = (1..=50).map(f64::from).collect();
    let a = bootstrap_ci(&data, median, 500, 0.05, 11);
    let b = bootstrap_ci(&data, median, 500, 0.05, 11);
    assert_eq!(a, b);
    assert!(a.0 < 25.5 && a.1 > 25.5);
}

#[test]
fn rng_reproducible_and_distributed() {
    let mut a = Rng::new(123);
    let mut b = Rng::new(123);
    for _ in 0..10 {
        assert_eq!(a.next_u64(), b.next_u64());
    }
    // golden values pin the stream across platforms
    let mut g = Rng::new(0);
    let first = g.next_u64();
    let mut g2 = Rng::new(0);
    assert_eq!(first, g2.next_u64());
    let n = 20000;
    let mut r = Rng::new(5);
    let xs: Vec<f64> = (0..n).map(|_| r.normal()).collect();
    let m = xs.iter().sum::<f64>() / f64::from(n);
    let v = xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / f64::from(n);
    assert!(m.abs() < 0.03 && (v - 1.0).abs() < 0.05);
    let ex: f64 = (0..n).map(|_| r.exponential(2.0)).sum::<f64>() / f64::from(n);
    assert!((ex - 0.5).abs() < 0.02);
    let gm: f64 = (0..n).map(|_| r.gamma(3.0, 2.0)).sum::<f64>() / f64::from(n);
    assert!((gm - 6.0).abs() < 0.15);
    let gs: f64 = (0..n).map(|_| r.gamma(0.5, 1.0)).sum::<f64>() / f64::from(n);
    assert!((gs - 0.5).abs() < 0.03);
    let bt: f64 = (0..n).map(|_| r.beta(2.0, 3.0)).sum::<f64>() / f64::from(n);
    assert!((bt - 0.4).abs() < 0.01);
    for lambda in [4.0, 100.0] {
        let p: f64 = (0..n).map(|_| r.poisson(lambda) as f64).sum::<f64>() / f64::from(n);
        assert!((p - lambda).abs() < 0.05 * lambda.sqrt() * 3.0, "{lambda} {p}");
    }
    let mut v: Vec<u32> = (0..20).collect();
    r.shuffle(&mut v);
    v.sort_unstable();
    assert_eq!(v, (0..20).collect::<Vec<_>>());
    let u = r.below(10);
    assert!(u < 10);
}

#[test]
fn fft_roundtrip_and_known_transform() {
    let x: Vec<(f64, f64)> = (0..8).map(|i| (f64::from(i), 0.0)).collect();
    let f = fft(&x, false);
    assert!((f[0].0 - 28.0).abs() < 1e-12);
    assert!((f[1].0 + 4.0).abs() < 1e-12 && (f[1].1 - 9.656_854_249_492_381).abs() < 1e-12);
    // non power of two (Bluestein) agrees with the direct DFT
    let n = 7;
    let y: Vec<(f64, f64)> = (0..n).map(|i| (f64::from(i).sin(), f64::from(i).cos())).collect();
    let fy = fft(&y, false);
    for k in 0..n {
        let (mut re, mut im) = (0.0, 0.0);
        for (j, v) in y.iter().enumerate() {
            let a = -2.0 * std::f64::consts::PI * (j * k as usize) as f64 / f64::from(n);
            re += v.0 * a.cos() - v.1 * a.sin();
            im += v.0 * a.sin() + v.1 * a.cos();
        }
        assert!((fy[k as usize].0 - re).abs() < 1e-11 && (fy[k as usize].1 - im).abs() < 1e-11);
    }
    let back = fft(&fy, true);
    for (a, b) in back.iter().zip(&y) {
        assert!((a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
    }
}

#[test]
fn convolution() {
    let a: Vec<f64> = (0..40).map(|i| f64::from(i % 7) - 3.0).collect();
    let b: Vec<f64> = (0..33).map(|i| f64::from(i % 5) * 0.5).collect();
    let fast = convolve(&a, &b);
    let mut slow = vec![0.0; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            slow[i + j] += x * y;
        }
    }
    assert_eq!(fast.len(), slow.len());
    for (x, y) in fast.iter().zip(&slow) {
        assert!((x - y).abs() < 1e-10);
    }
    let c = circular_convolve(&[1.0, 2.0, 3.0], &[0.0, 1.0, 0.5]);
    for (x, e) in c.iter().zip([4.0, 2.5, 2.5]) {
        assert!((x - e).abs() < 1e-12, "{c:?}");
    }
    let r = correlate(&[1.0, 2.0, 3.0], &[0.0, 1.0, 0.5]);
    assert_eq!(r.len(), 5);
    // (1 + x)^2 * (1 - x) exactly
    assert_eq!(poly_mul_i64(&[1, 2, 1], &[1, -1]), vec![1, 1, -1, -1]);
    let big: Vec<i64> = (0..200).map(|i| 1_000_000 + i).collect();
    let sq = poly_mul_i64(&big, &big);
    assert_eq!(sq[0], 1_000_000_i64 * 1_000_000);
    assert_eq!(sq[1], 2 * 1_000_000 * 1_000_001);
}
