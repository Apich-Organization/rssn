//! # Seeded random number generation
//!
//! A small, fully deterministic xoshiro256++ generator (seeded through
//! SplitMix64) and samplers for the common distributions. The same seed
//! always yields the same stream on every platform.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use std::f64::consts::PI;

/// Deterministic xoshiro256++ generator.
#[derive(Debug, Clone)]
pub struct Rng {
    s: [u64; 4],
    spare_normal: Option<f64>,
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Rng {
    /// Creates a generator from a 64-bit seed.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        let mut st = seed;
        let s = [splitmix(&mut st), splitmix(&mut st), splitmix(&mut st), splitmix(&mut st)];
        Self { s, spare_normal: None }
    }

    /// Next raw 64-bit output.
    pub fn next_u64(&mut self) -> u64 {
        let r = self.s[0].wrapping_add(self.s[3]).rotate_left(23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// Uniform sample in `[0, 1)` with 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0
    }

    /// Uniform sample in `[a, b)`.
    pub fn uniform_range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.uniform()
    }

    /// Uniform integer in `0..n` (`n > 0`), without modulo bias.
    pub fn below(&mut self, n: u64) -> u64 {
        if n <= 1 {
            return 0;
        }
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let v = self.next_u64();
            if v < zone {
                return v % n;
            }
        }
    }

    /// Standard normal sample (Box-Muller, caches the second value).
    pub fn normal(&mut self) -> f64 {
        if let Some(v) = self.spare_normal.take() {
            return v;
        }
        let u1 = 1.0 - self.uniform();
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        self.spare_normal = Some(r * (2.0 * PI * u2).sin());
        r * (2.0 * PI * u2).cos()
    }

    /// Normal sample with the given mean and standard deviation.
    pub fn normal_with(&mut self, mean: f64, sd: f64) -> f64 {
        mean + sd * self.normal()
    }

    /// Exponential sample with the given rate.
    pub fn exponential(&mut self, rate: f64) -> f64 {
        -(1.0 - self.uniform()).ln() / rate
    }

    /// Gamma sample (shape `k > 0`, scale `theta`), Marsaglia-Tsang.
    pub fn gamma(&mut self, k: f64, theta: f64) -> f64 {
        if k < 1.0 {
            let u = 1.0 - self.uniform();
            return self.gamma(k + 1.0, theta) * u.powf(1.0 / k);
        }
        let d = k - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = 1.0 - self.uniform();
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v * theta;
            }
        }
    }

    /// Beta(`a`, `b`) sample.
    pub fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a, 1.0);
        let y = self.gamma(b, 1.0);
        x / (x + y)
    }

    /// Poisson sample (Knuth for small means, normal-corrected rejection
    /// via gamma splitting for large ones).
    pub fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let (mut k, mut p) = (0u64, 1.0);
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1;
            }
        }
        // Split: Poisson(lambda) via the order-statistics of a gamma variate.
        let m = (lambda * 7.0 / 8.0).floor();
        let g = self.gamma(m, 1.0);
        if g < lambda {
            m as u64 + self.poisson(lambda - g)
        } else {
            self.binomial_count(m as u64 - 1, lambda / g)
        }
    }

    fn binomial_count(&mut self, n: u64, p: f64) -> u64 {
        (0..n).filter(|_| self.uniform() < p).count() as u64
    }

    /// Binomial(`n`, `p`) sample.
    pub fn binomial(&mut self, n: u64, p: f64) -> u64 {
        self.binomial_count(n, p)
    }

    /// Fisher-Yates shuffle in place.
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }

    /// Fills `out` with standard normal samples.
    pub fn fill_normal(&mut self, out: &mut [f64]) {
        for o in out {
            *o = self.normal();
        }
    }
}
