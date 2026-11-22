//! A small deterministic pseudo-random generator for the samplers.

/// PCG-XSH-RR generator with a 64-bit state and 32-bit output.
///
/// The generator is deliberately tiny and dependency-free: the samplers only
/// need reproducible uniforms and standard normals, and every entry point takes
/// an explicit seed so that a run can be repeated exactly.
#[derive(Clone, Debug)]
pub struct Pcg32 {
    /// LCG state.
    state: u64,
    /// Stream selector, always odd.
    inc: u64,
    /// Second Box-Muller variate waiting to be returned.
    spare: Option<f64>,
}

impl Pcg32 {
    /// LCG multiplier of the reference PCG implementation.
    const MULTIPLIER: u64 = 6_364_136_223_846_793_005;
    /// `2^-53`, the scale of the 53-bit uniform.
    const SCALE: f64 = 1.0 / 9_007_199_254_740_992.0;

    /// Creates a generator from a seed and a stream selector.
    #[must_use]
    pub fn new(
        seed: u64,
        stream: u64,
    ) -> Self {
        let mut rng = Self {
            state: 0,
            inc: (stream << 1) | 1,
            spare: None,
        };
        // The reference seeding discards two outputs so that the first
        // returned value depends on the whole seed.
        rng.next_u32();
        rng.state = rng.state.wrapping_add(seed);
        rng.next_u32();
        rng
    }

    /// Creates a generator on the default stream.
    #[must_use]
    pub fn seeded(seed: u64) -> Self {
        Self::new(seed, 0)
    }

    /// Next 32-bit output.
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old.wrapping_mul(Self::MULTIPLIER).wrapping_add(self.inc);
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rotation = (old >> 59) as u32;
        xorshifted.rotate_right(rotation)
    }

    /// Next uniform variate in `[0, 1)` with 53 bits of resolution.
    pub fn next_f64(&mut self) -> f64 {
        let high = u64::from(self.next_u32() >> 5);
        let low = u64::from(self.next_u32() >> 6);
        (high * (1u64 << 26) + low) as f64 * Self::SCALE
    }

    /// Next standard-normal variate by the Box-Muller transform.
    ///
    /// The second variate of each pair is cached and returned by the next call.
    pub fn normal(&mut self) -> f64 {
        if let Some(cached) = self.spare.take() {
            return cached;
        }
        let u1 = 1.0 - self.next_f64();
        let u2 = self.next_f64();
        let radius = (-2.0 * u1.ln()).sqrt();
        let angle = std::f64::consts::TAU * u2;
        self.spare = Some(radius * angle.sin());
        radius * angle.cos()
    }
}

#[cfg(test)]
mod tests {
    use super::Pcg32;

    #[test]
    fn streams_are_reproducible_and_seed_dependent() {
        let mut a = Pcg32::seeded(42);
        let mut b = Pcg32::seeded(42);
        let mut c = Pcg32::seeded(43);
        let mut same = true;
        for _ in 0..1000 {
            let x = a.next_u32();
            let y = b.next_u32();
            same &= x == y;
            assert!(x != c.next_u32(), "two seeds produced the same output");
        }
        assert!(same, "the same seed produced different streams");
    }

    #[test]
    fn uniforms_lie_in_the_unit_interval() {
        let mut rng = Pcg32::seeded(7);
        for _ in 0..10_000 {
            let u = rng.next_f64();
            assert!((0.0..1.0).contains(&u), "uniform {u} outside [0, 1)");
        }
    }

    #[test]
    fn normals_have_unit_moments() {
        // 100 000 draws: the standard error of the mean is 0.0032 and of the
        // variance about 0.0045, so 0.03 and 0.05 are generous tolerances.
        let mut rng = Pcg32::seeded(2024);
        let n = 100_000;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for _ in 0..n {
            let z = rng.normal();
            sum += z;
            sum_sq += z * z;
        }
        let mean = sum / n as f64;
        let variance = sum_sq / n as f64 - mean * mean;
        assert!(mean.abs() < 0.03, "mean = {mean}");
        assert!((variance - 1.0).abs() < 0.05, "variance = {variance}");
    }
}
