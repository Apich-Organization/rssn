//! `kernels::random`: the seeded xoshiro256++ generator against an
//! independent reference implementation, and the samplers against the
//! moments and distribution functions of the laws they draw from. Every
//! seed is fixed, so the statistical checks are deterministic.

use rssn::kernels::random::Rng;

/// The first four SplitMix64 outputs for the seed 0 (the reference values
/// published with the algorithm).
const SPLITMIX_SEED_0: [u64; 4] = [
    0xE220_A839_7B1D_CDAF,
    0x6E78_9E6A_A1B9_65F4,
    0x06C4_5D18_8009_454F,
    0xF88B_B8A8_724C_81EC,
];

/// xoshiro256++ as published, written independently of the crate.
struct Reference {
    s: [u64; 4],
}

impl Reference {
    fn next(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }
}

fn mean_and_variance(samples: &[f64]) -> (f64, f64) {
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    let variance = samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, variance)
}

/// The Kolmogorov–Smirnov distance of the sample to the distribution
/// function `cdf`.
fn ks_distance(
    samples: &mut [f64],
    cdf: impl Fn(f64) -> f64,
) -> f64 {
    samples.sort_by(f64::total_cmp);
    let n = samples.len() as f64;
    let mut worst: f64 = 0.0;
    for (i, &x) in samples.iter().enumerate() {
        let f = cdf(x);
        worst = worst
            .max((f - i as f64 / n).abs())
            .max(((i + 1) as f64 / n - f).abs());
    }
    worst
}

/// Critical value of the KS statistic at the 0.1 % level.
fn ks_limit(n: usize) -> f64 {
    1.95 / (n as f64).sqrt()
}

#[test]
fn raw_stream_matches_the_reference_xoshiro256pp() {
    let mut reference = Reference { s: SPLITMIX_SEED_0 };
    let mut rng = Rng::new(0);
    for i in 0..1000 {
        assert_eq!(rng.next_u64(), reference.next(), "output {i}");
    }
    // Same seed, same stream; different seeds differ immediately.
    let (mut a, mut b, mut c) = (Rng::new(42), Rng::new(42), Rng::new(43));
    let first: Vec<u64> = (0..16).map(|_| a.next_u64()).collect();
    assert_eq!(first, (0..16).map(|_| b.next_u64()).collect::<Vec<_>>());
    assert_ne!(first, (0..16).map(|_| c.next_u64()).collect::<Vec<_>>());
    // A clone continues the same stream.
    let mut d = a.clone();
    assert_eq!(a.next_u64(), d.next_u64());
}

#[test]
fn uniform_samples_are_uniform_on_the_unit_interval() {
    let mut rng = Rng::new(1);
    let n = 200_000;
    let mut u: Vec<f64> = (0..n).map(|_| rng.uniform()).collect();
    assert!(u.iter().all(|&v| (0.0..1.0).contains(&v)));
    let (mean, variance) = mean_and_variance(&u);
    assert!(
        (mean - 0.5).abs() < 5.0 * (1.0 / 12.0_f64 / n as f64).sqrt(),
        "mean {mean}"
    );
    assert!((variance - 1.0 / 12.0).abs() < 2e-3, "variance {variance}");
    assert!(ks_distance(&mut u, |x| x) < ks_limit(n), "KS");
    // Serial correlation at lag one vanishes.
    let mut rng = Rng::new(2);
    let v: Vec<f64> = (0..n).map(|_| rng.uniform() - 0.5).collect();
    let lag: f64 = v.windows(2).map(|w| w[0] * w[1]).sum::<f64>() / (n as f64 - 1.0);
    assert!(
        (lag / (1.0 / 12.0)).abs() < 5.0 / (n as f64).sqrt(),
        "lag-1 correlation {}",
        lag * 12.0
    );
    // The range sampler is an affine image.
    let mut rng = Rng::new(3);
    let w: Vec<f64> = (0..50_000).map(|_| rng.uniform_range(-3.0, 7.0)).collect();
    assert!(w.iter().all(|&x| (-3.0..7.0).contains(&x)));
    let (m, s2) = mean_and_variance(&w);
    assert!(
        (m - 2.0).abs() < 0.1 && (s2 - 100.0 / 12.0).abs() < 0.2,
        "{m} {s2}"
    );
}

#[test]
fn below_is_unbiased_and_in_range() {
    let mut rng = Rng::new(4);
    assert_eq!(rng.below(0), 0);
    assert_eq!(rng.below(1), 0);
    // A bound that does not divide 2^64 (the case modulo bias would show
    // in): chi-square of ten classes against the uniform law.
    let draws = 200_000_u64;
    let mut counts = [0_u64; 10];
    for _ in 0..draws {
        counts[rng.below(10) as usize] += 1;
    }
    let expected = draws as f64 / 10.0;
    let chi2: f64 = counts
        .iter()
        .map(|&c| (c as f64 - expected).powi(2) / expected)
        .sum();
    // 9 degrees of freedom, p = 0.001 at 27.9.
    assert!(chi2 < 27.9, "chi-square {chi2} for {counts:?}");
    // A huge bound stays below it.
    let big = u64::MAX / 3 * 2 + 7;
    assert!((0..1000).all(|_| rng.below(big) < big));
}

#[test]
fn normal_samples_have_the_standard_normal_law() {
    let mut rng = Rng::new(5);
    let n = 200_000;
    let mut x: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let (mean, variance) = mean_and_variance(&x);
    assert!(mean.abs() < 5.0 / (n as f64).sqrt(), "mean {mean}");
    assert!(
        (variance - 1.0).abs() < 5.0 * (2.0 / n as f64).sqrt(),
        "variance {variance}"
    );
    let skew = x.iter().map(|v| v.powi(3)).sum::<f64>() / n as f64;
    let kurtosis = x.iter().map(|v| v.powi(4)).sum::<f64>() / n as f64;
    assert!(
        skew.abs() < 0.03 && (kurtosis - 3.0).abs() < 0.08,
        "skew {skew}, kurtosis {kurtosis}"
    );
    // Standard normal distribution function through the error function
    // series (independent of the crate).
    let phi = |z: f64| 0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2));
    assert!(ks_distance(&mut x, phi) < ks_limit(n), "KS");
    // The tails: P(|Z| > 3) = 0.0026998.
    let tail = x.iter().filter(|v| v.abs() > 3.0).count() as f64 / n as f64;
    assert!((tail - 0.0026998).abs() < 5e-4, "tail {tail}");
    // Location and scale, and the bulk fill.
    let mut rng = Rng::new(6);
    let y: Vec<f64> = (0..100_000).map(|_| rng.normal_with(10.0, 2.5)).collect();
    let (m, v) = mean_and_variance(&y);
    assert!(
        (m - 10.0).abs() < 0.03 && (v - 6.25).abs() < 0.15,
        "{m} {v}"
    );
    let mut buffer = vec![0.0; 100_000];
    Rng::new(7).fill_normal(&mut buffer);
    let (m, v) = mean_and_variance(&buffer);
    assert!(m.abs() < 0.02 && (v - 1.0).abs() < 0.03, "{m} {v}");
    // `fill_normal` and repeated `normal` agree for the same seed (the
    // cached Box-Muller partner is part of the stream).
    let mut one = Rng::new(8);
    let mut two = Rng::new(8);
    let mut filled = [0.0; 7];
    two.fill_normal(&mut filled);
    for value in filled {
        assert_eq!(one.normal().to_bits(), value.to_bits());
    }
}

/// `erf` by its Maclaurin series for small arguments and a continued
/// fraction for the complement beyond, to about 1e-12.
fn erf(x: f64) -> f64 {
    let a = x.abs();
    let value = if a < 3.0 {
        let mut term = a;
        let mut sum = a;
        for n in 1..100 {
            term *= -a * a / n as f64;
            sum += term / (2 * n + 1) as f64;
        }
        sum * 2.0 / std::f64::consts::PI.sqrt()
    } else {
        // erfc(a) = e^{-a²}/(a√π) / (1 + 1/(2a²)/(1 + 2/(2a²)/(1 + ...)))
        let mut f = 0.0;
        for k in (1..60).rev() {
            f = k as f64 / 2.0 / (a + f);
        }
        1.0 - (-a * a).exp() / (std::f64::consts::PI.sqrt() * (a + f))
    };
    value.copysign(x)
}

#[test]
fn exponential_samples_have_the_memoryless_law() {
    let mut rng = Rng::new(9);
    let n = 200_000;
    let rate = 2.5;
    let mut x: Vec<f64> = (0..n).map(|_| rng.exponential(rate)).collect();
    assert!(x.iter().all(|&v| v >= 0.0 && v.is_finite()));
    let (mean, variance) = mean_and_variance(&x);
    assert!(
        (mean - 1.0 / rate).abs() < 5.0 / (rate * (n as f64).sqrt()),
        "mean {mean}"
    );
    assert!(
        (variance - 1.0 / (rate * rate)).abs() < 5e-3,
        "variance {variance}"
    );
    assert!(
        ks_distance(&mut x, |t| 1.0 - (-rate * t).exp()) < ks_limit(n),
        "KS"
    );
    // Memorylessness: P(X > 1 | X > 0.5) = P(X > 0.5).
    let beyond_half = x.iter().filter(|&&v| v > 0.5).count() as f64;
    let beyond_one = x.iter().filter(|&&v| v > 1.0).count() as f64;
    assert!((beyond_one / beyond_half - (-rate * 0.5_f64).exp()).abs() < 0.01);
}

#[test]
fn gamma_and_beta_samples_have_the_right_moments() {
    let n = 200_000;
    // Shapes below one use the boosting step, above one Marsaglia-Tsang.
    for (shape, scale, seed) in [
        (0.3_f64, 2.0_f64, 10_u64),
        (1.0, 1.0, 11),
        (2.5, 0.5, 12),
        (9.0, 3.0, 13),
    ] {
        let mut rng = Rng::new(seed);
        let x: Vec<f64> = (0..n).map(|_| rng.gamma(shape, scale)).collect();
        assert!(x.iter().all(|&v| v > 0.0 && v.is_finite()));
        let (mean, variance) = mean_and_variance(&x);
        let (want_mean, want_variance) = (shape * scale, shape * scale * scale);
        let sd = (want_variance / n as f64).sqrt();
        assert!(
            (mean - want_mean).abs() < 5.0 * sd,
            "shape {shape}: mean {mean} vs {want_mean}"
        );
        assert!(
            (variance / want_variance - 1.0).abs() < 0.05,
            "shape {shape}: variance {variance} vs {want_variance}"
        );
    }
    // Gamma(1, θ) is exponential: compare the distribution function.
    let mut rng = Rng::new(14);
    let mut x: Vec<f64> = (0..n).map(|_| rng.gamma(1.0, 2.0)).collect();
    assert!(
        ks_distance(&mut x, |t| 1.0 - (-t / 2.0).exp()) < ks_limit(n),
        "KS gamma(1)"
    );
    // Beta(a, b): mean a/(a+b), variance ab/((a+b)²(a+b+1)), support (0, 1).
    for (a, b, seed) in [(2.0_f64, 5.0_f64, 15_u64), (0.5, 0.5, 16), (4.0, 1.0, 17)] {
        let mut rng = Rng::new(seed);
        let x: Vec<f64> = (0..n).map(|_| rng.beta(a, b)).collect();
        assert!(x.iter().all(|&v| (0.0..=1.0).contains(&v)));
        let (mean, variance) = mean_and_variance(&x);
        let want_variance = a * b / ((a + b).powi(2) * (a + b + 1.0));
        assert!(
            (mean - a / (a + b)).abs() < 5.0 * (want_variance / n as f64).sqrt(),
            "Beta({a},{b}) mean {mean}"
        );
        assert!(
            (variance / want_variance - 1.0).abs() < 0.05,
            "Beta({a},{b}) variance {variance} vs {want_variance}"
        );
    }
}

#[test]
fn poisson_and_binomial_counts_have_the_right_laws() {
    let n = 100_000;
    // Both branches: Knuth's product (λ < 30) and the gamma splitting.
    for (lambda, seed) in [
        (0.7_f64, 20_u64),
        (4.0, 21),
        (29.0, 22),
        (30.0, 23),
        (150.0, 24),
        (1000.0, 25),
    ] {
        let mut rng = Rng::new(seed);
        let x: Vec<f64> = (0..n).map(|_| rng.poisson(lambda) as f64).collect();
        let (mean, variance) = mean_and_variance(&x);
        assert!(
            (mean - lambda).abs() < 5.0 * (lambda / n as f64).sqrt(),
            "λ = {lambda}: mean {mean}"
        );
        assert!(
            (variance / lambda - 1.0).abs() < 0.05,
            "λ = {lambda}: variance {variance}"
        );
    }
    // The probability mass of Poisson(4) at 0..8 by chi-square.
    let mut rng = Rng::new(26);
    let mut counts = [0_u64; 9];
    for _ in 0..n {
        let k = rng.poisson(4.0) as usize;
        counts[k.min(8)] += 1;
    }
    let mut pmf = [0.0_f64; 9];
    let mut term = (-4.0_f64).exp();
    for k in 0..8 {
        pmf[k] = term;
        term *= 4.0 / (k + 1) as f64;
    }
    pmf[8] = 1.0 - pmf[..8].iter().sum::<f64>();
    let chi2: f64 = (0..9)
        .map(|k| (counts[k] as f64 - n as f64 * pmf[k]).powi(2) / (n as f64 * pmf[k]))
        .sum();
    // 8 degrees of freedom, p = 0.001 at 26.1.
    assert!(chi2 < 26.1, "chi-square {chi2}");
    // Binomial(n, p): mean np, variance np(1-p), bounded by n.
    for (trials, p, seed) in [(10_u64, 0.3_f64, 27_u64), (50, 0.9, 28), (200, 0.5, 29)] {
        let mut rng = Rng::new(seed);
        let x: Vec<f64> = (0..20_000)
            .map(|_| rng.binomial(trials, p) as f64)
            .collect();
        assert!(x.iter().all(|&v| v <= trials as f64));
        let (mean, variance) = mean_and_variance(&x);
        let (want_mean, want_variance) = (trials as f64 * p, trials as f64 * p * (1.0 - p));
        assert!(
            (mean - want_mean).abs() < 5.0 * (want_variance / 20_000.0).sqrt(),
            "Bin({trials},{p}) mean {mean}"
        );
        assert!(
            (variance / want_variance - 1.0).abs() < 0.06,
            "Bin({trials},{p}) variance {variance}"
        );
    }
    assert_eq!(Rng::new(30).binomial(25, 0.0), 0);
    assert_eq!(Rng::new(31).binomial(25, 1.0), 25);
}

#[test]
fn shuffle_is_a_uniform_permutation() {
    let mut rng = Rng::new(40);
    let mut v: Vec<u32> = (0..100).collect();
    rng.shuffle(&mut v);
    let mut sorted = v.clone();
    sorted.sort_unstable();
    assert_eq!(
        sorted,
        (0..100).collect::<Vec<_>>(),
        "a permutation of the input"
    );
    assert_ne!(v, sorted, "actually moved");
    // All 24 orderings of four elements occur about equally often.
    let draws = 48_000;
    let mut counts = std::collections::HashMap::new();
    for _ in 0..draws {
        let mut w = [0_u8, 1, 2, 3];
        rng.shuffle(&mut w);
        *counts.entry(w).or_insert(0_u64) += 1;
    }
    assert_eq!(counts.len(), 24);
    let expected = draws as f64 / 24.0;
    let chi2: f64 = counts
        .values()
        .map(|&c| (c as f64 - expected).powi(2) / expected)
        .sum();
    // 23 degrees of freedom, p = 0.001 at 49.7.
    assert!(chi2 < 49.7, "chi-square {chi2}");
    // Short inputs are left alone.
    let mut empty: [u8; 0] = [];
    rng.shuffle(&mut empty);
    let mut one = [9_u8];
    rng.shuffle(&mut one);
    assert_eq!(one, [9]);
}
