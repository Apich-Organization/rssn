//! Sample quantiles of every Hyndman-Fan type and bootstrap intervals.
//!
//! | operator | value |
//! |---|---|
//! | `quantile_hf(data, p, kind)` | the sample quantile at probability `p` by Hyndman-Fan definition `kind` (1 to 9; 7 is the `NumPy` and R default) |
//! | `bootstrap_ci(data, statistic, alpha, resamples, seed)` | `list(lo, hi)`, the percentile bootstrap interval of level `1 - alpha` for a statistic named by a symbol |
//!
//! The statistic of `bootstrap_ci` is one of `mean`, `median`, `std`
//! (sample standard deviation), `variance` (sample variance), `mad` or `iqr`.
//! The resampling is seeded, so the interval is reproducible. The number of
//! resamples and the sample size are capped so that one call stays bounded.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::count;
use super::bridge::name;
use super::bridge::num;
use super::bridge::sample;
use super::bridge::scalar;
use super::bridge::vec_out;
use crate::kernels::robust;

/// Most resamples of one bootstrap.
const MAX_RESAMPLES: usize = 5_000;
/// Most random draws (sample size times resamples) of one bootstrap.
const MAX_DRAWS: usize = 5_000_000;
/// Largest seed that is an exact integer in a double.
const MAX_SEED: f64 = 9_007_199_254_740_992.0;

fn quantile_hf(a: &[Data]) -> Option<Data> {
    let x = sample(a, 0)?;
    let p = scalar(a, 1)?;
    let kind = count(a, 2, 9)?;
    ((0.0..=1.0).contains(&p) && kind >= 1).then_some(())?;
    num(robust::quantile(x, p, kind as u8))
}

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

fn variance(x: &[f64]) -> f64 {
    if x.len() < 2 {
        return 0.0;
    }
    let m = mean(x);
    x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (x.len() - 1) as f64
}

fn deviation(x: &[f64]) -> f64 {
    variance(x).sqrt()
}

fn median_abs_deviation(x: &[f64]) -> f64 {
    robust::mad(x, false)
}

/// The statistic a symbol names, over a resample.
fn statistic(name: &str) -> Option<fn(&[f64]) -> f64> {
    let stat: fn(&[f64]) -> f64 = match name {
        | "mean" => mean,
        | "median" => robust::median,
        | "std" => deviation,
        | "variance" => variance,
        | "mad" => median_abs_deviation,
        | "iqr" => robust::iqr,
        | _ => return None,
    };
    Some(stat)
}

fn bootstrap_ci(a: &[Data]) -> Option<Data> {
    let x = sample(a, 0)?;
    let stat = statistic(name(a, 1)?)?;
    let alpha = scalar(a, 2)?;
    (alpha > 0.0 && alpha < 1.0).then_some(())?;
    let resamples = count(a, 3, MAX_RESAMPLES)?;
    (resamples >= 1 && x.len().checked_mul(resamples)? <= MAX_DRAWS).then_some(())?;
    let seed = scalar(a, 4)?;
    (seed >= 0.0 && seed.fract() == 0.0 && seed <= MAX_SEED).then_some(())?;
    let (lo, hi) = robust::bootstrap_ci(x, stat, resamples, alpha, seed as u64);
    vec_out(vec![lo, hi])
}

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "quantile_hf",
        params: "vss",
        eval: quantile_hf,
    },
    Entry {
        name: "bootstrap_ci",
        params: "vwsss",
        eval: bootstrap_ci,
    },
];

#[cfg(test)]
mod tests {
    use crate::rules::numeric::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> Vec<f64> {
        simplify(&[numeric()], src)
            .replace(['(', ')', ','], " ")
            .split_whitespace()
            .filter(|word| *word != "list")
            .map(|word| word.parse().unwrap_or(f64::NAN))
            .collect()
    }

    #[test]
    fn every_quantile_definition_matches_its_formula() {
        // Type 7 interpolates at (n - 1) p: the quartile of 1..=5 is 2.
        assert_eq!(run("quantile_hf(list(1, 2, 3, 4, 5), 0.25, 7)"), vec![2.0]);
        // Type 1 is the inverse empirical distribution, type 2 averages at the jumps.
        assert_eq!(run("quantile_hf(list(1, 2, 3, 4), 0.5, 1)"), vec![2.0]);
        assert_eq!(run("quantile_hf(list(1, 2, 3, 4), 0.5, 2)"), vec![2.5]);
        // Type 6 (Weibull): h = (n + 1) p = 1.25 on 1..=4 gives 1.25.
        assert_eq!(run("quantile_hf(list(1, 2, 3, 4), 0.25, 6)"), vec![1.25]);
        // The extremes are the extremes in every definition.
        for kind in 1..=9 {
            assert_eq!(
                run(&format!("quantile_hf(list(3, 1, 2), 0, {kind})")),
                vec![1.0]
            );
            assert_eq!(
                run(&format!("quantile_hf(list(3, 1, 2), 1, {kind})")),
                vec![3.0]
            );
        }
        // It agrees with the median of the descriptive-statistics rule set.
        assert_eq!(run("quantile_hf(list(5, 1, 4, 2, 3), 0.5, 7)"), vec![3.0]);
    }

    #[test]
    fn the_bootstrap_interval_brackets_the_statistic_and_is_reproducible() {
        let data = "list(4.1, 5.2, 4.8, 5.9, 5.1, 4.4, 5.6, 5.0, 4.7, 5.3)";
        let first = run(&format!("bootstrap_ci({data}, mean, 0.1, 400, 7)"));
        let second = run(&format!("bootstrap_ci({data}, mean, 0.1, 400, 7)"));
        assert_eq!(first, second);
        let mean = 5.01;
        assert!(first[0] < mean && mean < first[1], "{first:?}");
        // The standard error of the mean is about 0.17, so a 90 percent
        // interval is roughly 0.56 wide.
        assert!((0.2..1.0).contains(&(first[1] - first[0])), "{first:?}");
        // A constant sample has a degenerate interval.
        assert_eq!(
            run("bootstrap_ci(list(2, 2, 2), median, 0.05, 50, 1)"),
            vec![2.0, 2.0]
        );
    }

    #[test]
    fn invalid_arguments_are_left_as_requests() {
        for src in [
            "quantile_hf(list(), 0.5, 7)",
            "quantile_hf(list(1, 2), 1.5, 7)",
            "quantile_hf(list(1, 2), 0.5, 10)",
            "quantile_hf(list(1, a), 0.5, 7)",
            "bootstrap_ci(list(1, 2, 3), nonesuch, 0.1, 10, 1)",
            "bootstrap_ci(list(1, 2, 3), mean, 0, 10, 1)",
            "bootstrap_ci(list(1, 2, 3), mean, 0.1, 0, 1)",
            "bootstrap_ci(list(1, 2, 3), mean, 0.1, 100000000, 1)",
            "bootstrap_ci(list(1, 2, 3), mean, 0.1, 10, -1)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
