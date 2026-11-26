//! Sample autocorrelation, partial autocorrelation and differencing.

/// Arithmetic mean of a slice; zero for an empty slice.
#[must_use]
pub(crate) fn mean_of(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    crate::kernels::simd::sum_f64(x) / x.len() as f64
}

/// Sample autocorrelation function.
///
/// Returns `max_lag + 1` values with `result[0] = 1`. The lag-`k` estimate is
/// `Σ_{t>k} (x_t - x̄)(x_{t-k} - x̄) / Σ_t (x_t - x̄)²`, the usual biased
/// estimator with the full-sample denominator. Lags beyond `x.len() - 1` are
/// not produced, and a constant series yields zeros after lag zero.
#[must_use]
pub fn acf(
    x: &[f64],
    max_lag: usize,
) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let limit = max_lag.min(n - 1);
    let mut out = vec![0.0; limit + 1];
    out[0] = 1.0;
    if limit == 0 {
        return out;
    }
    let mean = mean_of(x);
    let mut centered = crate::kernels::reuse::f64_scratch(n);
    centered.extend(x.iter().map(|&value| value - mean));
    let denominator = crate::kernels::simd::norm2_sq_f64(&centered);
    if denominator <= 0.0 || !denominator.is_finite() {
        return out;
    }
    if limit >= 16 && n >= 256 {
        let lags = crate::kernels::par::map_indexed(limit, 16, |idx| {
            let lag = idx + 1;
            crate::kernels::simd::dot_f64(&centered[lag..], &centered[..n - lag]) / denominator
        });
        out[1..].copy_from_slice(&lags);
    } else {
        for lag in 1..=limit {
            let numerator = crate::kernels::simd::dot_f64(&centered[lag..], &centered[..n - lag]);
            out[lag] = numerator / denominator;
        }
    }
    out
}

/// Partial autocorrelation function, from the Durbin–Levinson recursion.
///
/// Returns `max_lag + 1` values with `result[0] = 1`; entry `k` is the last
/// coefficient of the order-`k` autoregression fitted to the sample
/// autocorrelations.
#[must_use]
pub fn pacf(
    x: &[f64],
    max_lag: usize,
) -> Vec<f64> {
    let r = acf(x, max_lag);
    let limit = r.len().saturating_sub(1);
    let mut out = vec![0.0; limit + 1];
    if out.is_empty() {
        return out;
    }
    out[0] = 1.0;
    if limit == 0 {
        return out;
    }
    out[1] = r[1];
    if limit == 1 {
        return out;
    }

    let mut previous = crate::kernels::reuse::f64_scratch(limit + 1);
    previous.resize(limit + 1, 0.0);
    previous[1] = r[1];
    let mut next = crate::kernels::reuse::f64_scratch(limit + 1);
    next.resize(limit + 1, 0.0);
    for order in 2..=limit {
        let mut numerator = r[order];
        let mut denominator = 1.0;
        for j in 1..order {
            numerator -= previous[j] * r[order - j];
            denominator -= previous[j] * r[j];
        }
        let reflection = if denominator.abs() > 1e-300 {
            numerator / denominator
        } else {
            0.0
        };
        out[order] = reflection;

        next[..=order].fill(0.0);
        for j in 1..order {
            next[j] = previous[j] - reflection * previous[order - j];
        }
        next[order] = reflection;
        std::mem::swap(&mut previous, &mut next);
    }
    out
}

/// `d`-fold first difference.
///
/// Each pass replaces `x` by `x[t+1] - x[t]`, so the result is `d` entries
/// shorter than the input (and empty once the series is exhausted).
#[must_use]
pub fn difference(
    x: &[f64],
    d: usize,
) -> Vec<f64> {
    let mut current = x.to_vec();
    for _ in 0..d {
        if current.len() < 2 {
            return Vec::new();
        }
        for i in 0..(current.len() - 1) {
            current[i] = current[i + 1] - current[i];
        }
        current.pop();
    }
    current
}

/// Adds one cumulative-sum level, seeded by `start`.
fn integrate_once(
    x: &[f64],
    start: f64,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(x.len() + 1);
    out.push(start);
    let mut running = start;
    for value in x {
        running += value;
        out.push(running);
    }
    out
}

/// Inverse of [`difference`].
///
/// `initial` holds the first `initial.len()` values of the undifferenced
/// series, which fixes the integration order. The returned vector has
/// `d.len() + initial.len()` entries. An empty `initial` copies `d` unchanged.
#[must_use]
pub fn integrate(
    d: &[f64],
    initial: &[f64],
) -> Vec<f64> {
    let order = initial.len();
    if order == 0 {
        return d.to_vec();
    }
    if initial.iter().any(|v| !v.is_finite()) {
        return Vec::new();
    }

    // Seeds are Δʲy₀ for j = 0..order-1, taken from the leading values.
    let mut level = initial.to_vec();
    let mut seeds = Vec::with_capacity(order);
    if let Some(first) = level.first().copied() {
        seeds.push(first);
    }
    for step in 1..order {
        for i in 0..(order - step) {
            level[i] = level[i + 1] - level[i];
        }
        if let Some(first) = level.first().copied() {
            seeds.push(first);
        }
    }

    let mut current = d.to_vec();
    for step in (0..order).rev() {
        let seed = seeds.get(step).copied().unwrap_or(0.0);
        current = integrate_once(&current, seed);
    }
    current
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(
        a: f64,
        b: f64,
    ) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn acf_of_a_five_point_series_by_hand() {
        // x = 1..5, mean 3, Σ(x - x̄)² = 10.
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let r = acf(&x, 4);
        assert_eq!(r.len(), 5);
        assert!(close(r[0], 1.0));
        assert!(close(r[1], 0.4));
        assert!(close(r[2], -0.1));
        assert!(close(r[3], -0.4));
        assert!(close(r[4], -0.4));
    }

    #[test]
    fn acf_handles_short_and_constant_series() {
        assert!(acf(&[], 3).is_empty());
        let constant = acf(&[2.0, 2.0, 2.0], 5);
        assert_eq!(constant.len(), 3);
        assert!(close(constant[0], 1.0));
        assert!(close(constant[1], 0.0));
        // The lag count is capped at n - 1.
        assert_eq!(acf(&[1.0, 2.0], 10).len(), 2);
    }

    #[test]
    fn pacf_of_an_ar_one_matches_its_acf_at_lag_one() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let r = acf(&x, 4);
        let p = pacf(&x, 4);
        assert!(close(p[0], 1.0));
        assert!(close(p[1], r[1]));
        assert_eq!(p.len(), r.len());
    }

    /// Deterministic linear congruential generator in `[-1, 1)`.
    fn lcg(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let unit = (state >> 11) as f64 / (1_u64 << 53) as f64;
            unit * 2.0 - 1.0
        }
    }

    #[test]
    fn pacf_of_an_ar_two_vanishes_beyond_lag_two() {
        let mut noise = lcg(20_240_517);
        let mut x = vec![0.0; 2000];
        for t in 2..2000 {
            x[t] = 0.5 * x[t - 1] - 0.3 * x[t - 2] + noise();
        }
        let p = pacf(&x, 8);
        assert!(p[1] > 0.2, "first partial autocorrelation {}", p[1]);
        assert!(p[2].abs() > 0.1, "second partial autocorrelation {}", p[2]);
        for lag in 3..=8 {
            assert!(p[lag].abs() < 0.1, "lag {lag} gives {}", p[lag]);
        }
    }

    #[test]
    fn difference_and_integrate_round_trip() {
        let y = [1.0, 3.0, 7.0, 12.0, 18.0, 25.0];
        for d in 0..=3 {
            let w = difference(&y, d);
            assert_eq!(w.len(), y.len() - d);
            let back = integrate(&w, &y[..d]);
            assert_eq!(back.len(), y.len());
            for (a, b) in back.iter().zip(y.iter()) {
                assert!(close(*a, *b), "order {d}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn difference_of_a_short_series_is_empty() {
        assert!(difference(&[], 1).is_empty());
        assert!(difference(&[1.0], 1).is_empty());
        assert_eq!(difference(&[1.0, 2.0], 0).len(), 2);
    }

    #[test]
    fn acf_parallel_path_matches_direct_reference() {
        let n = 320;
        let max_lag = 28;
        let x: Vec<f64> = (0..n)
            .map(|i| (i as f64 * 0.17).sin() + 0.3 * (i as f64 * 0.05).cos())
            .collect();
        let got = acf(&x, max_lag);
        assert_eq!(got.len(), max_lag + 1);
        let mean = x.iter().sum::<f64>() / n as f64;
        let denom: f64 = x.iter().map(|v| (v - mean).powi(2)).sum();
        assert!(close(got[0], 1.0));
        for lag in 1..=max_lag {
            let num: f64 = (lag..n).map(|t| (x[t] - mean) * (x[t - lag] - mean)).sum();
            let expected = num / denom;
            assert!((got[lag] - expected).abs() < 1e-12, "lag {lag}");
        }
    }
}
