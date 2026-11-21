//! Overlap-save convolution against direct and exact dyadic reference sums.

use num_bigint::BigInt;
use num_traits::ToPrimitive;
use rssn::kernels::fftconv::OverlapSavePlan;
use rssn::kernels::fftconv::convolve;
use rssn::kernels::fftconv::convolve_overlap_save;
use rssn::kernels::signal::FilterError;
use rssn::kernels::signal::SignalError;

fn direct(
    samples: &[f64],
    taps: &[f64],
) -> Vec<f64> {
    if samples.is_empty() || taps.is_empty() {
        return Vec::new();
    }
    (0..samples.len() + taps.len() - 1)
        .map(|index| {
            taps.iter()
                .enumerate()
                .filter_map(|(lag, tap)| {
                    index
                        .checked_sub(lag)
                        .and_then(|at| samples.get(at))
                        .map(|sample| sample * tap)
                })
                .sum()
        })
        .collect()
}

fn exact_coefficients(
    samples: &[i64],
    taps: &[i64],
) -> Vec<BigInt> {
    let mut output = vec![BigInt::from(0); samples.len() + taps.len() - 1];
    for (time, &sample) in samples.iter().enumerate() {
        let sample = BigInt::from(sample);
        for (lag, &tap) in taps.iter().enumerate() {
            output[time + lag] += &sample * tap;
        }
    }
    output
}

#[test]
fn preserves_subnormal_taps_with_large_signed_recordings() {
    let integers: Vec<i64> = (0..605).map(|i| (11 * i + 7) % 23 - 11).collect();
    let tap_integers: Vec<i64> = (0..17).map(|i| (7 * i + 2) % 13 - 6).collect();
    let samples: Vec<_> = integers
        .iter()
        .map(|&i| i as f64 * 2.0_f64.powi(900))
        .collect();
    let taps: Vec<_> = tap_integers
        .iter()
        .map(|&i| i.signum() as f64 * f64::from_bits(i.unsigned_abs()))
        .collect();
    let exact = exact_coefficients(&integers, &tap_integers);
    let direct = direct(&samples, &taps);
    for got in [
        OverlapSavePlan::new(&taps)
            .unwrap()
            .apply(&samples)
            .unwrap(),
        convolve(&samples, &taps),
        convolve_overlap_save(&taps, &samples).unwrap(),
    ] {
        assert_eq!(got.len(), exact.len());
        for ((actual, integer), direct) in got.into_iter().zip(&exact).zip(&direct) {
            let expected = integer.to_f64().unwrap();
            assert_eq!(*direct, expected * 2.0_f64.powi(-174));
            assert!((actual * 2.0_f64.powi(174) - expected).abs() < 2.0e-12);
        }
    }
}

#[test]
fn rescales_a_forward_transform_that_would_overflow() {
    let sample = 2.0_f64.powi(1020);
    let tap = 2.0_f64.powi(-900);
    let samples = vec![sample; 703];
    let taps = [tap; 17];
    let expected = direct(&samples, &taps);
    let got = convolve_overlap_save(&samples, &taps).unwrap();
    assert_eq!(got.len(), expected.len());
    for (actual, expected) in got.into_iter().zip(expected) {
        assert!((actual / expected - 1.0).abs() < 6.0e-15);
    }
    // Well-separated blocks with different scales normalize independently.
    let tiny = f64::from_bits(1);
    let mut separated = vec![0.0; 700];
    separated[..111].fill(sample * 0.5);
    separated[550..].fill(tiny);
    let normalized = convolve_overlap_save(&separated, &[1.0; 17]).unwrap();
    assert_eq!(normalized[631].to_bits(), 17);
    assert_eq!(normalized[715].to_bits(), 1);
}

#[test]
fn rounds_finite_subnormal_sums_after_accumulation() {
    let value = 2.0_f64.powi(-540);
    let samples = vec![value; 451];
    let taps = [value; 129];
    let got = convolve_overlap_save(&samples, &taps).unwrap();
    let exact = exact_coefficients(&[1; 451], &[1; 129]);
    assert_eq!(got.len(), exact.len());
    for (got, integer) in got.into_iter().zip(exact) {
        let count = integer.to_u64().unwrap();
        // At exact half-way cases FFT roundoff can move the final value to
        // either nearest subnormal. All other cases have a unique rounding.
        if count % 64 != 32 {
            assert_eq!(got.to_bits(), (count + 32) / 64);
        }
    }
}

#[test]
fn prepared_batch_matches_direct_and_single_recording_across_workers() {
    let taps: Vec<_> = (0..17).map(|i| (f64::from(i) * 0.4).cos()).collect();
    let plan = OverlapSavePlan::new(&taps).unwrap();
    // 112 valid coefficients per block for this filter. Include empty, short,
    // exact-block and multi-block outputs with both padded ends.
    let recordings: Vec<Vec<_>> = [0, 1, 95, 96, 97, 112, 224, 2513, 65_937]
        .into_iter()
        .map(|n| (0..n).map(|i| (f64::from(i) * 0.23).sin()).collect())
        .collect();
    let refs: Vec<_> = recordings.iter().map(Vec::as_slice).collect();
    let single: Vec<_> = refs.iter().map(|x| plan.apply(x).unwrap()).collect();
    for workers in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        let batch = pool.install(|| plan.apply_batch(&refs).unwrap());
        assert_eq!(batch, single);
        assert!(
            batch
                .iter()
                .flatten()
                .zip(single.iter().flatten())
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        for (actual, samples) in batch.iter().zip(&refs) {
            let expected = direct(samples, &taps);
            assert_eq!(actual.len(), expected.len());
            for (actual, expected) in actual.iter().zip(expected) {
                assert!((actual - expected).abs() < 1.0e-13);
            }
        }
    }
    assert!(plan.apply_batch(&[]).unwrap().is_empty());
    assert!(convolve_overlap_save(&[], &[f64::NAN]).unwrap().is_empty());
    assert!(convolve_overlap_save(&[f64::NAN], &[]).unwrap().is_empty());
}

#[test]
fn reports_the_first_batch_error_and_rejects_invalid_filters() {
    assert!(matches!(
        OverlapSavePlan::new(&[]),
        Err(FilterError::InvalidParameters)
    ));
    assert!(matches!(
        OverlapSavePlan::new(&[f64::INFINITY]),
        Err(FilterError::Signal(SignalError::NonFiniteInput))
    ));
    let plan = OverlapSavePlan::new(&[2.0_f64.powi(100); 17]).unwrap();
    let overflow = [2.0_f64.powi(1020); 7];
    let invalid = [f64::NAN];
    assert_eq!(
        plan.apply_batch(&[&[], &overflow, &invalid, &[1.0]]),
        Err(FilterError::Signal(SignalError::NonFiniteResult))
    );
    assert_eq!(
        plan.apply_batch(&[&[], &invalid, &overflow, &[1.0]]),
        Err(FilterError::Signal(SignalError::NonFiniteInput))
    );
}
