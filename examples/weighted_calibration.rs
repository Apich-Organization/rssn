//! Use weighted empirical statistics to summarize reference-reading offsets.

use rssn::kernels::robust::weighted_ecdf;
use rssn::kernels::robust::weighted_quantile;

fn main() {
    let expected_metres = [0.0, 0.5, 1.0, 1.5, 2.0];
    let observed_metres = [0.03, 0.55, 1.04, 1.55, 3.5];
    // The last sensor was out of calibration, so it supplies no evidence.
    let weights = [1.0, 1.0, 4.0, 1.0, 0.0];
    let residuals =
        std::array::from_fn::<_, 5, _>(|index| observed_metres[index] - expected_metres[index]);

    let center = weighted_quantile(&residuals, &weights, 0.5);
    let lower = weighted_quantile(&residuals, &weights, 0.1);
    let upper = weighted_quantile(&residuals, &weights, 0.9);
    let within_limit = weighted_ecdf(&residuals, &weights, 0.045);

    assert!((center - 0.04).abs() < 1e-12);
    assert!((within_limit - 5.0 / 7.0).abs() < 1e-12);
    println!("Offset: {center:.3} m; weighted 10–90% interval: {lower:.3}–{upper:.3} m");
    println!("Weight at or below 45 mm: {:.1}%", within_limit * 100.0);
}
