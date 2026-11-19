//! # Seeded random number generation
//!
//! A small, fully deterministic xoshiro256++ generator (seeded through
//! `SplitMix64`) and samplers for the common distributions. The same seed
//! always yields the same integer stream on every platform.
//!
//! [`Rng::jump`] and [`Rng::long_jump`] advance the raw generator by fixed
//! powers of two. Clone and jump the same initial stream to assign
//! deterministic subsequences to parallel workers. The generator is intended
//! for numerical simulation and sampling.
//!
//! The Gaussian sampler caches one spare normal value. Mixed distribution
//! calls consume the raw stream in an order determined by the samplers.

use std::f64::consts::PI;

/// Deterministic xoshiro256++ generator.
#[derive(Debug, Clone)]
pub struct Rng {
    s: [u64; 4],
    spare_normal: Option<f64>,
}

const fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Rng {
    /// Creates a generator from a 64-bit seed.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        let mut st = seed;
        let s = [
            splitmix(&mut st),
            splitmix(&mut st),
            splitmix(&mut st),
            splitmix(&mut st),
        ];
        Self {
            s,
            spare_normal: None,
        }
    }

    /// Next raw 64-bit output.
    pub const fn next_u64(&mut self) -> u64 {
        let r = self.s[0]
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
        r
    }

    /// Advances the raw stream by `2^128` steps and clears the cached normal.
    pub const fn jump(&mut self) {
        self.apply_jump([
            0x180E_C6D3_3CFD_0ABA,
            0xD5A6_1266_F0C9_392C,
            0xA958_2618_E03F_C9AA,
            0x39AB_DC45_29B1_661C,
        ]);
    }

    /// Advances the raw stream by `2^192` steps and clears the cached normal.
    pub const fn long_jump(&mut self) {
        self.apply_jump([
            0x76E1_5D3E_FEFD_CBBF,
            0xC500_4E44_1C52_2FB3,
            0x7771_0069_854E_E241,
            0x3910_9BB0_2ACB_E635,
        ]);
    }

    const fn apply_jump(
        &mut self,
        polynomial: [u64; 4],
    ) {
        // Blackman and Vigna's xoshiro256++ jump polynomials.
        let mut state = [0; 4];
        let mut word = 0;
        while word < polynomial.len() {
            let mut bit = 0;
            while bit < 64 {
                if (polynomial[word] >> bit) & 1 != 0 {
                    state[0] ^= self.s[0];
                    state[1] ^= self.s[1];
                    state[2] ^= self.s[2];
                    state[3] ^= self.s[3];
                }
                self.next_u64();
                bit += 1;
            }
            word += 1;
        }
        self.s = state;
        self.spare_normal = None;
    }

    /// Uniform sample in `[0, 1)` with 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0
    }

    /// Uniform sample in `[a, b)`.
    pub fn uniform_range(
        &mut self,
        a: f64,
        b: f64,
    ) -> f64 {
        let u = self.uniform();
        let width = b - a;
        if width.is_finite() {
            a + width * u
        } else {
            a * (1.0 - u) + b * u
        }
    }

    /// Uniform integer in `0..n` (`n > 0`), without modulo bias.
    pub const fn below(
        &mut self,
        n: u64,
    ) -> u64 {
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
    pub fn normal_with(
        &mut self,
        mean: f64,
        sd: f64,
    ) -> f64 {
        mean + sd * self.normal()
    }

    /// Exponential sample with a positive finite rate. Invalid rates give NaN.
    pub fn exponential(
        &mut self,
        rate: f64,
    ) -> f64 {
        if !rate.is_finite() || rate <= 0.0 {
            return f64::NAN;
        }
        -(-self.uniform()).ln_1p() / rate
    }

    /// Gamma sample with positive finite shape `k` and scale `theta`.
    /// Invalid parameters give NaN.
    pub fn gamma(
        &mut self,
        k: f64,
        theta: f64,
    ) -> f64 {
        if !k.is_finite() || k <= 0.0 || !theta.is_finite() || theta <= 0.0 {
            return f64::NAN;
        }
        if k < 1.0 {
            return (self.log_gamma_scaled(k, 1.0) + theta.ln()).exp();
        }
        self.gamma_shape_ge_one(k) * theta
    }

    fn gamma_shape_ge_one(
        &mut self,
        k: f64,
    ) -> f64 {
        let d = k - 1.0 / 3.0;
        let c = (1.0 / 3.0) / d.sqrt();
        loop {
            let x = self.normal();
            let w = c * x;
            if w <= -1.0 {
                continue;
            }
            let u = 1.0 - self.uniform();
            let v = (1.0 + w).powi(3);
            // Marsaglia-Tsang squeeze and log acceptance. Use log1p for
            // the small relative change in v at large shape.
            let delta = w * (3.0 + w * (3.0 + w));
            if u < 1.0 - 0.0331 * x.powi(4) || u.ln() < 0.5 * x * x + d * (delta.ln_1p() - delta) {
                return d * v;
            }
        }
    }

    fn log_gamma_scaled(
        &mut self,
        shape: f64,
        scale: f64,
    ) -> f64 {
        if shape < 1.0 {
            let log_uniform = (-self.uniform()).ln_1p();
            scale * self.gamma_shape_ge_one(shape + 1.0).ln() + (scale / shape) * log_uniform
        } else {
            scale * self.gamma_shape_ge_one(shape).ln()
        }
    }

    /// Beta sample with positive finite shapes `a` and `b`.
    /// Invalid shapes give NaN.
    pub fn beta(
        &mut self,
        a: f64,
        b: f64,
    ) -> f64 {
        if !a.is_finite() || !b.is_finite() || a <= 0.0 || b <= 0.0 {
            return f64::NAN;
        }
        if a >= 1.0 && b >= 1.0 {
            let x = self.gamma_shape_ge_one(a);
            let y = self.gamma_shape_ge_one(b);
            let total = x + y;
            if total.is_finite() {
                return x / total;
            }
            let scale = x.max(y);
            return (x / scale) / (x / scale + y / scale);
        }
        // Gamma variates can both underflow for small shapes. Their log
        // ratio still distinguishes the two endpoints of the beta sample.
        let scale = a.min(b).min(1.0);
        let log_x = self.log_gamma_scaled(a, scale);
        let log_y = self.log_gamma_scaled(b, scale);
        let spread = (log_y - log_x) / scale;
        if spread > 0.0 {
            let odds = (-spread).exp();
            odds / (1.0 + odds)
        } else {
            1.0 / (1.0 + spread.exp())
        }
    }

    /// Poisson sample by Knuth's algorithm at small means and gamma
    /// splitting at large means. NaN and nonpositive means give zero.
    /// Positive means or results outside the `u64` range saturate.
    pub fn poisson(
        &mut self,
        mut lambda: f64,
    ) -> u64 {
        if lambda.is_nan() || lambda <= 0.0 {
            return 0;
        }
        if lambda >= u64::MAX as f64 {
            return u64::MAX;
        }

        let mut events = 0_u64;
        while lambda >= 30.0 {
            // The time of the m-th arrival in a rate-one process is Gamma(m, 1).
            let m = (lambda * 0.875).floor() as u64;
            let arrival = self.gamma_shape_ge_one(m as f64);
            if arrival >= lambda {
                return events.saturating_add(self.binomial(m - 1, lambda / arrival));
            }
            events = events.saturating_add(m);
            lambda -= arrival;
        }

        let limit = (-lambda).exp();
        let mut product = self.uniform();
        while product > limit {
            events = events.saturating_add(1);
            product *= self.uniform();
        }
        events
    }

    /// Binomial sample with `n` trials and a probability in `[0, 1]`.
    /// Degenerate probabilities give 0 or `n`. Invalid probabilities give 0.
    pub fn binomial(
        &mut self,
        mut n: u64,
        mut p: f64,
    ) -> u64 {
        if !(0.0..=1.0).contains(&p) || p == 0.0 {
            return 0;
        }
        if p >= 1.0 {
            return n;
        }
        let mut count = 0;
        while n > 64 {
            // The k-th order statistic of n uniform trials has law
            // Beta(k, n+1-k). Given that value, the remaining relevant
            // uniforms lie on its side of the threshold independently.
            let k = n / 2 + 1;
            let split = self.beta(k as f64, (n - k + 1) as f64);
            if split <= p {
                count += k;
                n -= k;
                p = (p - split) / (1.0 - split);
            } else {
                n = k - 1;
                p /= split;
            }
        }
        count + (0..n).filter(|_| self.uniform() < p).count() as u64
    }

    /// Fisher-Yates shuffle in place.
    pub fn shuffle<T>(
        &mut self,
        v: &mut [T],
    ) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }

    /// Fills `out` with standard normal samples.
    pub fn fill_normal(
        &mut self,
        out: &mut [f64],
    ) {
        for o in out {
            *o = self.normal();
        }
    }
}
