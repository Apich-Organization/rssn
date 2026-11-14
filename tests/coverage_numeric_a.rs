//! Coverage tests: Fresnel / Faddeeva, online covariance, deterministic
//! parallel helpers, low-rank kernels and the structure-preserving /
//! stiff integrators. Public API only.

#![allow(clippy::needless_range_loop)]

use num_complex::Complex64;
use rssn::kernels::dense::Mat;
use rssn::kernels::fresnel::erf_complex;
use rssn::kernels::fresnel::faddeeva_w;
use rssn::kernels::fresnel::fresnel_c;
use rssn::kernels::fresnel::fresnel_f;
use rssn::kernels::fresnel::fresnel_fg;
use rssn::kernels::fresnel::fresnel_g;
use rssn::kernels::fresnel::fresnel_s;
use rssn::kernels::fresnel::fresnel_sc;
use rssn::kernels::fresnel::fresnel_sc_complex;
use rssn::kernels::fresnel::sin_cos_half_pi_square;
use rssn::kernels::lowrank::Gaussian;
use rssn::kernels::lowrank::LowRankError;
use rssn::kernels::lowrank::RandomizedOptions;
use rssn::kernels::lowrank::anderson;
use rssn::kernels::lowrank::lobpcg;
use rssn::kernels::lowrank::qr_pivoted;
use rssn::kernels::lowrank::randomized_svd;
use rssn::kernels::online_stats::OnlineCovariance;
use rssn::kernels::online_stats::OnlineStatsError;
use rssn::kernels::par;
use rssn::sim::integrators::backward_euler_step;
use rssn::sim::integrators::implicit_midpoint_step;
use rssn::sim::integrators::rosenbrock_solve;
use rssn::sim::integrators::symplectic_trajectory;
use rssn::sim::integrators::verlet_step;
use rssn::sim::integrators::yoshida4_step;

/// Tiny deterministic generator, uniform in [0, 1).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }

    fn range(
        &mut self,
        lo: f64,
        hi: f64,
    ) -> f64 {
        lo + (hi - lo) * self.next()
    }
}

fn close(
    a: f64,
    b: f64,
    tol: f64,
) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

// ---------------------------------------------------------------- fresnel

#[test]
fn fresnel_reference_values() {
    // Abramowitz-Stegun / DLMF table values.
    let table = [
        (0.0, 0.0, 0.0),
        (0.5, 0.064_732_432_859_999_3, 0.492_344_225_871_446_4),
        (1.0, 0.438_259_147_390_354_8, 0.779_893_400_376_822_8),
        (2.0, 0.343_415_678_363_698_2, 0.488_253_406_075_340_5),
    ];
    for (x, s, c) in table {
        let (gs, gc) = fresnel_sc(x);
        assert!((gs - s).abs() < 1e-13, "S({x}) = {gs}, want {s}");
        assert!((gc - c).abs() < 1e-13, "C({x}) = {gc}, want {c}");
        assert_eq!(fresnel_s(x), gs);
        assert_eq!(fresnel_c(x), gc);
    }
}

#[test]
fn fresnel_odd_and_limits() {
    let mut rng = Lcg(7);
    for _ in 0..200 {
        let x = rng.range(0.0, 30.0);
        let (s, c) = fresnel_sc(x);
        let (sn, cn) = fresnel_sc(-x);
        assert_eq!(s, -sn);
        assert_eq!(c, -cn);
    }
    for x in [1e3, 1e6, 1e12, 1e16] {
        let (s, c) = fresnel_sc(x);
        assert!((s - 0.5).abs() < 1.0 / (std::f64::consts::PI * x) * 1.01 + 1e-15);
        assert!((c - 0.5).abs() < 1.0 / (std::f64::consts::PI * x) * 1.01 + 1e-15);
    }
    assert_eq!(fresnel_sc(f64::INFINITY), (0.5, 0.5));
    assert_eq!(fresnel_sc(f64::NEG_INFINITY), (-0.5, -0.5));
    let (s, c) = fresnel_sc(f64::NAN);
    assert!(s.is_nan() && c.is_nan());
}

#[test]
fn fresnel_derivative_is_the_integrand() {
    let mut rng = Lcg(11);
    for _ in 0..100 {
        let x = rng.range(0.05, 6.0);
        let h = 1e-5;
        let th = |t: f64| std::f64::consts::FRAC_PI_2 * t * t;
        let ds = (fresnel_s(x + h) - fresnel_s(x - h)) / (2.0 * h);
        let dc = (fresnel_c(x + h) - fresnel_c(x - h)) / (2.0 * h);
        assert!((ds - th(x).sin()).abs() < 1e-5, "S'({x}) = {ds}");
        assert!((dc - th(x).cos()).abs() < 1e-5, "C'({x}) = {dc}");
    }
}

#[test]
fn fresnel_continuous_across_the_series_switch() {
    let e = 1e-12;
    let (s1, c1) = fresnel_sc(1.5 - e);
    let (s2, c2) = fresnel_sc(1.5 + e);
    assert!((s1 - s2).abs() < 1e-10 && (c1 - c2).abs() < 1e-10);
}

#[test]
fn fresnel_auxiliary_functions() {
    for x in [0.0_f64, 0.7, 1.5, 2.0, 3.5, 6.0, 9.0] {
        let (f, g) = fresnel_fg(x);
        let (s, c) = fresnel_sc(x);
        let (sn, cs) = sin_cos_half_pi_square(x);
        let f_ref = (0.5 - s) * cs - (0.5 - c) * sn;
        let g_ref = (0.5 - c) * cs + (0.5 - s) * sn;
        assert!((f - f_ref).abs() < 1e-12, "f({x})");
        assert!((g - g_ref).abs() < 1e-12, "g({x})");
        assert_eq!(fresnel_f(x), f);
        assert_eq!(fresnel_g(x), g);
    }
    // asymptotics f ~ 1/(pi x), g ~ 1/(pi^2 x^3)
    let x = 40.0;
    let (f, g) = fresnel_fg(x);
    assert!(((f * std::f64::consts::PI * x) - 1.0).abs() < 1e-5);
    assert!(((g * std::f64::consts::PI.powi(2) * x.powi(3)) - 1.0).abs() < 1e-3);
    assert_eq!(fresnel_fg(f64::INFINITY), (0.0, 0.0));
    assert!(fresnel_f(f64::NAN).is_nan());
}

#[test]
fn exact_phase_reduction() {
    // theta = pi x^2 / 2 is an exact multiple of pi/2 for these x.
    let cases = [
        (1.0, 1.0, 0.0),
        (2.0, 0.0, 1.0),
        (3.0, 1.0, 0.0),
        (4.0, 0.0, 1.0),
        (1.0e8, 0.0, 1.0),
    ];
    for (x, s, c) in cases {
        let (gs, gc) = sin_cos_half_pi_square(x);
        assert!((gs - s).abs() < 1e-12 && (gc - c).abs() < 1e-12, "x = {x}");
        let (ns, nc) = sin_cos_half_pi_square(-x);
        assert_eq!((gs, gc), (ns, nc), "even in x");
    }
    let (s, c) = sin_cos_half_pi_square(f64::INFINITY);
    assert!(s.is_nan() && c.is_nan());
    let (s, c) = sin_cos_half_pi_square(f64::NAN);
    assert!(s.is_nan() && c.is_nan());
    let (s, c) = sin_cos_half_pi_square(1e200);
    assert!((s, c) == (0.0, 1.0));
    // always on the unit circle
    let mut rng = Lcg(3);
    for _ in 0..200 {
        let x = rng.range(-1e6, 1e6);
        let (s, c) = sin_cos_half_pi_square(x);
        assert!((s * s + c * c - 1.0).abs() < 1e-14);
    }
}

#[test]
fn faddeeva_known_values_and_symmetries() {
    let w0 = faddeeva_w(Complex64::new(0.0, 0.0));
    assert!((w0 - Complex64::new(1.0, 0.0)).norm() < 1e-13);
    // w(i) = e erfc(1)
    let wi = faddeeva_w(Complex64::new(0.0, 1.0));
    assert!((wi - Complex64::new(0.427_583_576_155_807_04, 0.0)).norm() < 1e-12);
    // w(1) = e^-1 + 2i/sqrt(pi) Dawson(1)
    let w1 = faddeeva_w(Complex64::new(1.0, 0.0));
    assert!(
        (w1 - Complex64::new(0.367_879_441_171_442_33, 0.607_157_705_841_393_7)).norm() < 1e-12
    );
    let mut rng = Lcg(5);
    for _ in 0..200 {
        let z = Complex64::new(rng.range(-6.0, 6.0), rng.range(-6.0, 6.0));
        let w = faddeeva_w(z);
        // w(-z) = 2 exp(-z^2) - w(z)
        let r = faddeeva_w(-z) - (2.0 * (-(z * z)).exp() - w);
        assert!(
            r.norm() < 1e-10 * (1.0 + w.norm() + (z * z).re.abs().exp().min(1e12)),
            "z = {z}"
        );
        // w(-conj z) = conj w(z)
        let m = faddeeva_w(-z.conj()) - w.conj();
        assert!(m.norm() < 1e-11 * (1.0 + w.norm()), "z = {z}");
    }
    // far field: w(z) ~ i / (sqrt(pi) z)
    let z = Complex64::new(30.0, 20.0);
    let lead = Complex64::new(0.0, 1.0) / (std::f64::consts::PI.sqrt() * z);
    assert!(((faddeeva_w(z) - lead) / lead).norm() < 2e-3);
    assert!(faddeeva_w(Complex64::new(f64::NAN, 0.0)).re.is_nan());
}

#[test]
fn erf_complex_properties() {
    let e = erf_complex(Complex64::new(1.0, 1.0));
    assert!((e - Complex64::new(1.316_151_281_697_947_6, 0.190_453_469_237_834_7)).norm() < 1e-12);
    let mut rng = Lcg(9);
    for _ in 0..200 {
        let z = Complex64::new(rng.range(-4.0, 4.0), rng.range(-4.0, 4.0));
        let ez = erf_complex(z);
        assert!((erf_complex(-z) + ez).norm() < 1e-12, "odd, z = {z}");
        assert!(
            (erf_complex(z.conj()) - ez.conj()).norm() < 1e-12,
            "conj, z = {z}"
        );
    }
    // saturates to 1 in the right half plane
    assert!((erf_complex(Complex64::new(6.0, 0.1)) - Complex64::new(1.0, 0.0)).norm() < 1e-12);
    // continuity across |z| = 1.5
    let d = Complex64::new(1.0, 1.0).unscale(2f64.sqrt());
    let a = erf_complex(d * (1.5 - 1e-9));
    let b = erf_complex(d * (1.5 + 1e-9));
    assert!((a - b).norm() < 1e-8);
    assert!(erf_complex(Complex64::new(f64::NAN, 1.0)).re.is_nan());
}

#[test]
fn fresnel_complex_consistency() {
    // real axis agrees with the real kernel
    let ((s, si), (c, ci)) = fresnel_sc_complex(1.3, 0.0);
    let (rs, rc) = fresnel_sc(1.3);
    assert_eq!((s, si, c, ci), (rs, 0.0, rc, 0.0));
    // imaginary axis: S(iy) = -i S(y), C(iy) = i C(y)
    let ((sr, sim), (cr, cim)) = fresnel_sc_complex(0.0, 0.9);
    let (rs, rc) = fresnel_sc(0.9);
    assert_eq!((sr, sim, cr, cim), (0.0, -rs, 0.0, rc));
    // oddness off the axes
    let mut rng = Lcg(21);
    for _ in 0..100 {
        let (x, y) = (rng.range(0.1, 3.0), rng.range(0.1, 2.0));
        let ((s1, s2), (c1, c2)) = fresnel_sc_complex(x, y);
        let ((n1, n2), (m1, m2)) = fresnel_sc_complex(-x, -y);
        assert!((s1 + n1).abs() < 1e-11 && (s2 + n2).abs() < 1e-11);
        assert!((c1 + m1).abs() < 1e-11 && (c2 + m2).abs() < 1e-11);
    }
    // first-order expansion in the imaginary direction
    for x in [0.9_f64, 1.6, 2.5] {
        let h = 1e-5;
        let ((_, s_im), (_, c_im)) = fresnel_sc_complex(x, h);
        let th = std::f64::consts::FRAC_PI_2 * x * x;
        assert!((s_im / h - th.sin()).abs() < 1e-7, "S, x = {x}");
        assert!((c_im / h - th.cos()).abs() < 1e-7, "C, x = {x}");
    }
}

// ----------------------------------------------------------- online stats

fn two_pass(
    data: &[Vec<f64>],
    d: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = data.len() as f64;
    let mean: Vec<f64> = (0..d)
        .map(|j| data.iter().map(|r| r[j]).sum::<f64>() / n)
        .collect();
    let mut cov = vec![0.0; d * d];
    for r in data {
        for i in 0..d {
            for j in 0..d {
                cov[i * d + j] += (r[i] - mean[i]) * (r[j] - mean[j]);
            }
        }
    }
    (mean, cov)
}

fn sample_data(
    n: usize,
    d: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|_| {
            let base = rng.range(-3.0, 3.0);
            (0..d)
                .map(|j| base * (j as f64 + 1.0) + rng.range(-1.0, 1.0) + 100.0)
                .collect()
        })
        .collect()
}

#[test]
fn online_covariance_matches_two_pass() {
    let d = 3;
    let data = sample_data(500, d, 1);
    let mut acc = OnlineCovariance::new(d).unwrap();
    assert_eq!(acc.dimension(), d);
    assert_eq!(acc.count(), 0);
    assert!(acc.mean().is_none());
    assert!(acc.population_covariance().is_none());
    for r in &data {
        acc.push(r).unwrap();
    }
    assert_eq!(acc.count(), 500);
    let (mean, m2) = two_pass(&data, d);
    for j in 0..d {
        assert!(close(acc.mean().unwrap()[j], mean[j], 1e-12));
    }
    let pop = acc.population_covariance().unwrap();
    let samp = acc.sample_covariance().unwrap();
    for k in 0..d * d {
        assert!(close(pop[k], m2[k] / 500.0, 1e-10), "pop {k}");
        assert!(close(samp[k], m2[k] / 499.0, 1e-10), "sample {k}");
    }
    // symmetric, nonnegative diagonal
    for i in 0..d {
        assert!(samp[i * d + i] >= 0.0);
        for j in 0..d {
            assert_eq!(samp[i * d + j], samp[j * d + i]);
        }
    }
}

#[test]
fn online_covariance_merge_is_order_independent() {
    let d = 4;
    let data = sample_data(300, d, 2);
    let mut whole = OnlineCovariance::new(d).unwrap();
    for r in &data {
        whole.push(r).unwrap();
    }
    for split in [0, 1, 17, 150, 299, 300] {
        let (mut a, mut b) = (
            OnlineCovariance::new(d).unwrap(),
            OnlineCovariance::new(d).unwrap(),
        );
        for r in &data[..split] {
            a.push(r).unwrap();
        }
        for r in &data[split..] {
            b.push(r).unwrap();
        }
        a.merge(&b).unwrap();
        assert_eq!(a.count(), 300);
        let (w, m) = (
            whole.population_covariance().unwrap(),
            a.population_covariance().unwrap(),
        );
        for k in 0..d * d {
            assert!(close(m[k], w[k], 1e-10), "split {split} entry {k}");
        }
        for j in 0..d {
            assert!(close(a.mean().unwrap()[j], whole.mean().unwrap()[j], 1e-12));
        }
    }
}

#[test]
fn online_covariance_edge_cases() {
    assert_eq!(
        OnlineCovariance::new(0).unwrap_err(),
        OnlineStatsError::InvalidDimension
    );
    assert_eq!(
        OnlineCovariance::new(usize::MAX).unwrap_err(),
        OnlineStatsError::InvalidDimension
    );
    assert_eq!(
        OnlineCovariance::new(1 << 31).unwrap_err(),
        OnlineStatsError::InvalidDimension
    );
    let mut a = OnlineCovariance::new(2).unwrap();
    assert_eq!(
        a.push(&[1.0]).unwrap_err(),
        OnlineStatsError::DimensionMismatch
    );
    assert_eq!(
        a.push(&[1.0, f64::NAN]).unwrap_err(),
        OnlineStatsError::NonFiniteInput
    );
    assert_eq!(
        a.push(&[f64::INFINITY, 0.0]).unwrap_err(),
        OnlineStatsError::NonFiniteInput
    );
    assert_eq!(a.count(), 0);
    a.push(&[1.0, 2.0]).unwrap();
    // one sample: population covariance is zero, sample covariance undefined
    assert_eq!(a.population_covariance().unwrap(), vec![0.0; 4]);
    assert!(a.sample_covariance().is_none());
    let other = OnlineCovariance::new(3).unwrap();
    assert_eq!(
        a.merge(&other).unwrap_err(),
        OnlineStatsError::DimensionMismatch
    );
    // merging an empty accumulator changes nothing; merging into an empty one copies
    let empty = OnlineCovariance::new(2).unwrap();
    let before = a.clone();
    a.merge(&empty).unwrap();
    assert_eq!(a.count(), before.count());
    let mut fresh = OnlineCovariance::new(2).unwrap();
    fresh.merge(&a).unwrap();
    assert_eq!(fresh.mean().unwrap(), a.mean().unwrap());
    // overflow is reported and leaves the state intact
    let mut big = OnlineCovariance::new(1).unwrap();
    big.push(&[1.0e200]).unwrap();
    let err = big.push(&[-1.0e200]).unwrap_err();
    assert_eq!(err, OnlineStatsError::Overflow);
    assert_eq!(big.count(), 1);
    assert_eq!(big.mean().unwrap(), &[1.0e200]);
}

// -------------------------------------------------------------------- par

#[test]
fn par_maps_preserve_order_for_every_grain() {
    for n in [0usize, 1, 2, 3, 100, 1000] {
        for grain in [0usize, 1, 2, 7, 10_000] {
            let v = par::map_indexed(n, grain, |i| i * i + 1);
            assert_eq!(v, (0..n).map(|i| i * i + 1).collect::<Vec<_>>());
            let items: Vec<usize> = (0..n).collect();
            let w = par::map_slice(&items, grain, |x| x * 3);
            assert_eq!(w, (0..n).map(|i| i * 3).collect::<Vec<_>>());
        }
    }
}

#[test]
fn par_chunks_and_each_cover_every_element_once() {
    for (len, chunk) in [
        (0usize, 4usize),
        (1, 4),
        (10, 3),
        (100, 7),
        (100, 0),
        (64, 64),
        (64, 1000),
    ] {
        for grain in [0usize, 2, 1000] {
            let mut data = vec![0u32; len];
            par::for_chunks_mut(&mut data, chunk, grain, |_ci, c| {
                for x in c {
                    *x += 1; // each element touched exactly once
                }
            });
            assert!(data.iter().all(|&x| x == 1), "len {len} chunk {chunk}");
            // chunk index matches position
            let mut idx = vec![usize::MAX; len];
            let ch = chunk.max(1);
            par::for_chunks_mut(&mut idx, chunk, grain, |ci, c| {
                for x in c {
                    *x = ci;
                }
            });
            for (p, &ci) in idx.iter().enumerate() {
                assert_eq!(ci, p / ch);
            }
            let mut each = vec![0usize; len];
            par::for_each_mut(&mut each, grain, |i, x| *x = i + 1);
            assert_eq!(each, (1..=len).collect::<Vec<_>>());
        }
    }
}

#[test]
fn par_sum_fixed_is_bitwise_deterministic() {
    let f = |i: usize| ((i as f64) * 0.1).sin() * 1e-3 + 1.0 / (1.0 + i as f64);
    let n: usize = 100_003;
    for block in [0usize, 1, 64, 1000] {
        // reference: same association order, written serially
        let b = block.max(1);
        let mut reference = 0.0;
        for blk in 0..n.div_ceil(b) {
            let mut s = 0.0;
            for i in blk * b..((blk + 1) * b).min(n) {
                s += f(i);
            }
            reference += s;
        }
        let a = par::sum_fixed(n, block, 0, f);
        let seq = par::sum_fixed(n, block, usize::MAX, f);
        assert_eq!(a.to_bits(), reference.to_bits(), "block {block}");
        assert_eq!(seq.to_bits(), reference.to_bits(), "block {block}");
        for threads in [1usize, 3, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let t = pool.install(|| par::sum_fixed(n, block, 0, f));
            assert_eq!(
                t.to_bits(),
                reference.to_bits(),
                "block {block}, {threads} threads"
            );
        }
    }
    assert_eq!(par::sum_fixed(0, 8, 0, f), 0.0);
    assert_eq!(par::sum_fixed(1, 8, 0, |_| 2.5), 2.5);
}

// ----------------------------------------------------------------- lowrank

fn random_mat(
    m: usize,
    n: usize,
    seed: u64,
) -> Mat {
    let mut rng = Lcg(seed);
    let mut a = Mat::zeros(m, n);
    for v in a.data.iter_mut() {
        *v = rng.range(-1.0, 1.0);
    }
    a
}

#[test]
fn gaussian_stream_is_reproducible_and_standard() {
    let (mut a, mut b, mut c) = (Gaussian::new(42), Gaussian::new(42), Gaussian::new(44));
    let mut diff = false;
    let (mut s, mut s2) = (0.0, 0.0);
    let n = 40_000;
    for _ in 0..n {
        let (x, y, z) = (a.sample(), b.sample(), c.sample());
        assert_eq!(x.to_bits(), y.to_bits());
        assert!(x.is_finite());
        diff |= x != z;
        s += x;
        s2 += x * x;
    }
    assert!(diff);
    let mean = s / n as f64;
    let var = s2 / n as f64 - mean * mean;
    assert!(mean.abs() < 0.03, "mean {mean}");
    assert!((var - 1.0).abs() < 0.05, "var {var}");
}

#[test]
fn pivoted_qr_reconstructs_and_reveals_rank() {
    for (m, n) in [(6, 4), (4, 6), (5, 5), (1, 1), (7, 1), (1, 7)] {
        let a = random_mat(m, n, (m * 31 + n) as u64);
        let f = qr_pivoted(&a).unwrap();
        let k = m.min(n);
        assert_eq!(f.perm.len(), n);
        let mut seen = f.perm.clone();
        seen.sort_unstable();
        assert_eq!(seen, (0..n).collect::<Vec<_>>());
        assert_eq!(f.q.rows, m);
        assert_eq!(f.q.cols, k);
        // Q^T Q = I
        let qtq = f.q.transpose().matmul(&f.q).unwrap();
        for i in 0..k {
            for j in 0..k {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((qtq.at(i, j) - want).abs() < 1e-12);
            }
        }
        // A P = Q R
        let qr = f.q.matmul(&f.r).unwrap();
        for j in 0..n {
            for i in 0..m {
                assert!((qr.at(i, j) - a.at(i, f.perm[j])).abs() < 1e-12, "{m}x{n}");
            }
        }
        // |R_ii| non-increasing
        for i in 1..k {
            assert!(f.r.at(i, i).abs() <= f.r.at(i - 1, i - 1).abs() * (1.0 + 1e-12));
        }
        assert_eq!(f.rank(None), k);
    }
    // rank deficient: third column = first + second
    let mut a = random_mat(8, 3, 99);
    for i in 0..8 {
        let v = a.at(i, 0) + a.at(i, 1);
        a.set(i, 2, v);
    }
    let f = qr_pivoted(&a).unwrap();
    assert_eq!(f.rank(Some(1e-10)), 2);
    assert_eq!(
        qr_pivoted(&Mat::zeros(4, 3))
            .map(|f| f.rank(None))
            .unwrap_or(0),
        0
    );
}

#[test]
fn pivoted_qr_solve_and_errors() {
    let a = random_mat(9, 4, 5);
    let x_true = [1.0, -2.0, 0.5, 3.0];
    let b: Vec<f64> = (0..9)
        .map(|i| (0..4).map(|j| a.at(i, j) * x_true[j]).sum())
        .collect();
    let x = qr_pivoted(&a).unwrap().solve(&b, None).unwrap();
    for j in 0..4 {
        assert!((x[j] - x_true[j]).abs() < 1e-11);
    }
    assert!(matches!(
        qr_pivoted(&a).unwrap().solve(&b[..3], None),
        Err(LowRankError::Shape(_))
    ));
    let mut bad = a.clone();
    bad.set(0, 0, f64::NAN);
    assert!(matches!(qr_pivoted(&bad), Err(LowRankError::Shape(_))));
    let inconsistent = Mat {
        rows: 3,
        cols: 3,
        data: vec![1.0; 4],
    };
    assert!(matches!(
        qr_pivoted(&inconsistent),
        Err(LowRankError::Shape(_))
    ));
    // empty matrices must not panic
    let _ = qr_pivoted(&Mat::zeros(0, 0));
    let _ = qr_pivoted(&Mat::zeros(0, 3));
}

#[test]
fn randomized_svd_recovers_a_low_rank_operator() {
    let (m, n, r) = (40, 30, 3);
    let u = random_mat(m, r, 1);
    let v = random_mat(n, r, 2);
    let mut a = Mat::zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            a.set(
                i,
                j,
                (0..r)
                    .map(|t| u.at(i, t) * v.at(j, t) * (3 - t) as f64)
                    .sum(),
            );
        }
    }
    let mv = |x: &[f64], y: &mut [f64]| {
        for i in 0..m {
            y[i] = (0..n).map(|j| a.at(i, j) * x[j]).sum();
        }
    };
    let rmv = |x: &[f64], y: &mut [f64]| {
        for j in 0..n {
            y[j] = (0..m).map(|i| a.at(i, j) * x[i]).sum();
        }
    };
    let opts = RandomizedOptions::default();
    let s = randomized_svd(&mv, &rmv, m, n, r, &opts).unwrap();
    let exact = rssn::kernels::dense::svd(&a).unwrap();
    for t in 0..r {
        assert!(
            close(s.s[t], exact.s[t], 1e-8),
            "sigma {t}: {} vs {}",
            s.s[t],
            exact.s[t]
        );
    }
    // U S V^T reproduces A
    for i in 0..m {
        for j in 0..n {
            let rec: f64 = (0..r).map(|t| s.u.at(i, t) * s.s[t] * s.v.at(j, t)).sum();
            assert!((rec - a.at(i, j)).abs() < 1e-8);
        }
    }
    // reproducible for a fixed seed
    let s2 = randomized_svd(&mv, &rmv, m, n, r, &opts).unwrap();
    assert_eq!(s.s, s2.s);
    // invalid ranks
    assert!(matches!(
        randomized_svd(&mv, &rmv, m, n, 0, &opts),
        Err(LowRankError::Shape(_))
    ));
    assert!(matches!(
        randomized_svd(&mv, &rmv, m, n, 31, &opts),
        Err(LowRankError::Shape(_))
    ));
    // the zero operator must not panic
    let zero = |_: &[f64], y: &mut [f64]| y.iter_mut().for_each(|v| *v = 0.0);
    if let Ok(z) = randomized_svd(&zero, &zero, m, n, 2, &opts) {
        assert!(z.s.iter().all(|&v| v.abs() < 1e-12));
    }
}

#[test]
fn lobpcg_finds_the_smallest_laplacian_eigenvalues() {
    let n = 60;
    let lap = |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            let l = if i > 0 { x[i - 1] } else { 0.0 };
            let r = if i + 1 < n { x[i + 1] } else { 0.0 };
            y[i] = 2.0 * x[i] - l - r;
        }
    };
    let jac = |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            y[i] = x[i] / 2.0;
        }
    };
    let pairs = lobpcg(&lap, Some(&jac), n, 3, 1e-8, 2000, 7).unwrap();
    for (j, &lam) in pairs.values.iter().enumerate() {
        let want = 2.0 - 2.0 * ((j + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
        assert!((lam - want).abs() < 1e-7, "lambda {j}: {lam} vs {want}");
        // residual claim is honest
        let mut ax = vec![0.0; n];
        lap(&pairs.vectors[j], &mut ax);
        let r: f64 = ax
            .iter()
            .zip(&pairs.vectors[j])
            .map(|(a, x)| (a - lam * x).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(r < 1e-5, "residual {r}");
    }
    assert!(matches!(
        lobpcg(&lap, None, n, 0, 1e-8, 10, 1),
        Err(LowRankError::Shape(_))
    ));
    assert!(matches!(
        lobpcg(&lap, None, n, 21, 1e-8, 10, 1),
        Err(LowRankError::Shape(_))
    ));
}

#[test]
fn anderson_accelerates_fixed_point_iterations() {
    let g = |x: &[f64]| vec![x[0].cos()];
    let r = anderson(&g, &[1.0], 3, 1.0, 1e-12, 200).unwrap();
    assert!((r.x[0] - 0.739_085_133_215_160_7).abs() < 1e-10);
    assert!(r.residual <= 1e-12);
    assert!(r.evaluations < 40, "took {} evaluations", r.evaluations);
    // linear contraction in 3-D with known fixed point
    let a = [[0.5, 0.2, 0.0], [0.1, 0.4, 0.1], [0.0, 0.3, 0.6]];
    let b = [1.0, -1.0, 2.0];
    let lin = |x: &[f64]| {
        (0..3)
            .map(|i| (0..3).map(|j| a[i][j] * x[j]).sum::<f64>() + b[i])
            .collect()
    };
    let r = anderson(&lin, &[0.0; 3], 4, 1.0, 1e-12, 500).unwrap();
    let back = lin(&r.x);
    for i in 0..3 {
        assert!((back[i] - r.x[i]).abs() < 1e-10);
    }
    // invalid damping falls back to 1, wrong output length and bad start are errors
    assert!(anderson(&g, &[1.0], 3, -5.0, 1e-10, 200).is_ok());
    assert!(matches!(
        anderson(&g, &[f64::NAN], 3, 1.0, 1e-10, 10),
        Err(LowRankError::Shape(_))
    ));
    let wrong = |_: &[f64]| vec![0.0, 0.0];
    assert!(matches!(
        anderson(&wrong, &[1.0], 3, 1.0, 1e-10, 10),
        Err(LowRankError::Shape(_))
    ));
    // divergence is reported, not panicked on
    let blow = |x: &[f64]| vec![x[0] * 1e200 + 1.0];
    let _ = anderson(&blow, &[1.0], 2, 1.0, 1e-12, 50);
}

// ------------------------------------------------------------- integrators

fn osc(q: &[f64]) -> Vec<f64> {
    vec![-q[0]]
}

fn osc_error(
    order4: bool,
    steps: usize,
) -> f64 {
    let dt = 1.0 / steps as f64;
    let (mut q, mut p) = ([1.0], [0.0]);
    for _ in 0..steps {
        if order4 {
            yoshida4_step(&mut q, &mut p, 1.0, dt, &osc);
        } else {
            verlet_step(&mut q, &mut p, 1.0, dt, &osc);
        }
    }
    (q[0] - 1.0_f64.cos())
        .abs()
        .max((p[0] + 1.0_f64.sin()).abs())
}

#[test]
fn symplectic_steps_converge_at_their_order() {
    let r2 = osc_error(false, 20) / osc_error(false, 40);
    let r4 = osc_error(true, 10) / osc_error(true, 20);
    assert!((3.6..4.4).contains(&r2), "verlet ratio {r2}");
    assert!((13.0..19.0).contains(&r4), "yoshida ratio {r4}");
}

#[test]
fn symplectic_steps_are_time_reversible() {
    let accel = |q: &[f64]| vec![-q[0].sin(), -2.0 * q[1]];
    for order4 in [false, true] {
        let (mut q, mut p) = ([0.7, -0.3], [0.2, 0.5]);
        let (q0, p0) = (q, p);
        for dt in [0.05, -0.05] {
            for _ in 0..50 {
                if order4 {
                    yoshida4_step(&mut q, &mut p, 2.0, dt, &accel);
                } else {
                    verlet_step(&mut q, &mut p, 2.0, dt, &accel);
                }
            }
        }
        for k in 0..2 {
            assert!(
                (q[k] - q0[k]).abs() < 1e-12 && (p[k] - p0[k]).abs() < 1e-12,
                "order4 {order4}"
            );
        }
    }
}

#[test]
fn symplectic_trajectory_shape_and_energy() {
    let t = symplectic_trajectory(&[1.0], &[0.0], 1.0, 0.01, 1000, 2, osc);
    assert_eq!(t.len(), 1001);
    assert_eq!(t[0], (vec![1.0], vec![0.0]));
    let e0 = 0.5;
    for (q, p) in &t {
        let e = 0.5 * (q[0] * q[0] + p[0] * p[0]);
        assert!((e - e0).abs() < 1e-4);
    }
    // order >= 4 uses Yoshida, which is far more accurate
    let t4 = symplectic_trajectory(&[1.0], &[0.0], 1.0, 0.1, 100, 4, osc);
    let t2 = symplectic_trajectory(&[1.0], &[0.0], 1.0, 0.1, 100, 2, osc);
    let exact = 10.0_f64.cos();
    assert!((t4[100].0[0] - exact).abs() < 0.2 * (t2[100].0[0] - exact).abs() + 1e-9);
    assert_eq!(
        symplectic_trajectory(&[1.0], &[0.0], 1.0, 0.1, 0, 2, osc).len(),
        1
    );
}

#[test]
fn implicit_steps_match_closed_forms_on_linear_problems() {
    let lambda = 50.0;
    let f = |_: f64, y: &[f64]| vec![-lambda * y[0], -2.0 * y[1]];
    let h = 0.1;
    let be = backward_euler_step(&f, 0.0, &[1.0, 1.0], h).unwrap();
    assert!((be[0] - 1.0 / (1.0 + lambda * h)).abs() < 1e-8);
    assert!((be[1] - 1.0 / (1.0 + 2.0 * h)).abs() < 1e-8);
    let im = implicit_midpoint_step(&f, 0.0, &[1.0, 1.0], h).unwrap();
    assert!((im[0] - (1.0 - lambda * h / 2.0) / (1.0 + lambda * h / 2.0)).abs() < 1e-8);
    assert!((im[1] - (1.0 - h) / (1.0 + h)).abs() < 1e-8);
    // implicit midpoint is second order: error ratio 4 under halving
    let g = |_: f64, y: &[f64]| vec![-y[0]];
    let run = |n: usize| {
        let mut y = vec![1.0];
        for k in 0..n {
            y = implicit_midpoint_step(&g, k as f64 / n as f64, &y, 1.0 / n as f64).unwrap();
        }
        (y[0] - (-1.0_f64).exp()).abs()
    };
    let ratio = run(20) / run(40);
    assert!((3.5..4.5).contains(&ratio), "ratio {ratio}");
    // A-stable on the stiff problem: no blow-up for h*lambda = 5e3
    let stiff = backward_euler_step(&f, 0.0, &[1.0, 1.0], 100.0).unwrap();
    assert!(stiff[0].abs() < 1e-3);
    // zero-dimensional state is a trivial problem
    let empty = |_: f64, _: &[f64]| Vec::<f64>::new();
    assert_eq!(
        backward_euler_step(&empty, 0.0, &[], 0.1).unwrap(),
        Vec::<f64>::new()
    );
}

#[test]
fn rosenbrock_solves_stiff_and_nonstiff_problems() {
    let sol = rosenbrock_solve(|_, y: &[f64]| vec![-y[0]], 0.0, &[1.0], 1.0, 1e-6).unwrap();
    let (t_end, y_end) = sol.last().unwrap();
    assert!((t_end - 1.0).abs() < 1e-12);
    assert!((y_end[0] - (-1.0_f64).exp()).abs() < 1e-3);
    assert!(sol.windows(2).all(|w| w[1].0 > w[0].0));
    // stiff linear decay: stays bounded and nonnegative with a modest step count
    let stiff =
        rosenbrock_solve(|_, y: &[f64]| vec![-1000.0 * y[0]], 0.0, &[1.0], 1.0, 1e-4).unwrap();
    assert!(stiff.len() < 2000, "{} steps", stiff.len());
    assert!(
        stiff
            .iter()
            .all(|(_, y)| y[0] > -1e-3 && y[0] <= 1.0 + 1e-12)
    );
    assert!(stiff.last().unwrap().1[0].abs() < 1e-3);
    // degenerate interval returns the initial point only
    let one = rosenbrock_solve(|_, y: &[f64]| vec![-y[0]], 2.0, &[3.0], 2.0, 1e-6).unwrap();
    assert_eq!(one, vec![(2.0, vec![3.0])]);
    // backward integration
    let back =
        rosenbrock_solve(|_, y: &[f64]| vec![y[0]], 1.0, &[1.0_f64.exp()], 0.0, 1e-8).unwrap();
    assert!((back.last().unwrap().1[0] - 1.0).abs() < 1e-3);
}
