//! Euclidean lengths and directions across binary floating-point scales.

use rssn::kernels::vector::angle;
use rssn::kernels::vector::cosine_similarity;
use rssn::kernels::vector::distance;
use rssn::kernels::vector::is_parallel;
use rssn::kernels::vector::norm;
use rssn::kernels::vector::normalize;
use rssn::kernels::vector::normalize_into;

fn close(
    got: f64,
    expected: f64,
) {
    assert!(
        (got / expected - 1.0).abs() <= 2e-15,
        "{got:e} versus {expected:e}"
    );
}

#[test]
fn lengths_and_differences_rescale_overflowing_and_underflowing_squares() {
    for scale in [1e-200, 1e200] {
        let values = [3.0 * scale, -4.0 * scale];
        close(norm(&values), 5.0 * scale);
        close(distance(&values, &[0.0, 0.0]).unwrap(), 5.0 * scale);
        let opposite = [-values[0], -values[1]];
        close(distance(&values, &opposite).unwrap(), 10.0 * scale);
    }
    let tiny = f64::from_bits(1);
    assert_eq!(norm(&[tiny]).to_bits(), 1);
    assert_eq!(distance(&[tiny], &[0.0]).unwrap().to_bits(), 1);
    assert_eq!(norm(&[f64::MAX]), f64::MAX);
    assert_eq!(distance(&[f64::MAX], &[-f64::MAX]).unwrap(), f64::INFINITY);
    assert!(norm(&[f64::INFINITY, f64::NAN]).is_nan());
    assert_eq!(norm(&[]).to_bits(), 0);
}

#[test]
fn directions_remain_finite_beyond_the_range_of_lengths() {
    let unit = std::f64::consts::FRAC_1_SQRT_2;
    let tiny = f64::from_bits(1);
    for size in [tiny, 1e-200, 1.0, 1e200, f64::MAX] {
        let values = [size, -size, -0.0];
        let mut result = [4.0; 3];
        normalize_into(&values, &mut result).unwrap();
        assert_eq!(result, normalize(&values).unwrap().as_slice());
        close(result[0], unit);
        close(result[1], -unit);
        assert_eq!(result[2].to_bits(), (-0.0_f64).to_bits());
    }
    assert!(norm(&[f64::MAX, f64::MAX]).is_infinite());
}

#[test]
fn normalization_validates_before_writing_output() {
    let mut out = [17.0, -2.0];
    let original = out;
    assert!(normalize_into(&[1.0], &mut out).is_err());
    assert_eq!(out, original);
    assert!(normalize_into(&[0.0, -0.0], &mut out).is_err());
    assert_eq!(out, original);
    assert!(normalize(&[]).is_err());
    assert_eq!(
        distance(&[0.0, 0.0], &[]).unwrap_err(),
        "Dimension mismatch: v1.len() = 2, v2.len() = 0"
    );
}

#[test]
fn cosine_and_angle_accept_independently_scaled_vectors() {
    let a = [3e200, 4e200];
    let b = [-4e-200, 3e-200];
    assert!(cosine_similarity(&a, &b).unwrap().abs() < 1e-15);
    close(angle(&a, &b).unwrap(), std::f64::consts::FRAC_PI_2);
    assert!(is_parallel(&[f64::MAX, f64::MAX], &[1e-200, 1e-200], 1e-14).unwrap());

    let smallest = f64::from_bits(1);
    close(cosine_similarity(&[smallest], &[-smallest]).unwrap(), -1.0);
    close(
        angle(&[smallest], &[-smallest]).unwrap(),
        std::f64::consts::PI,
    );
    assert_eq!(angle(&[0.0], &[1.0]).unwrap().to_bits(), 0);
    assert!(angle(&[0.0], &[]).is_err());
    assert!(cosine_similarity(&[0.0], &[1.0]).is_err());
}
