//! Inverse real FFT range, shape, and repeated-output tests.

use std::f64::consts::TAU;

use num_complex::Complex64;
use rssn::kernels::transforms::IrfftError;
use rssn::kernels::transforms::IrfftPlan;
use rssn::kernels::transforms::irfft;

#[test]
fn impulse_spectra_preserve_subnormal_signals_and_avoid_overflow() {
    let tiny = f64::from_bits(1);
    for n in [2, 3, 4, 6, 7, 16, 17, 97, 129] {
        for amplitude in [tiny, -tiny] {
            let spectrum = vec![Complex64::new(amplitude, 0.0); n / 2 + 1];
            let output = irfft(&spectrum, n);
            assert_eq!(output[0].to_bits(), amplitude.to_bits(), "n = {n}");
            assert!(output[1..].iter().all(|&x| x == 0.0), "n = {n}");
        }
        for amplitude in [f64::MAX, -f64::MAX] {
            let spectrum = vec![Complex64::new(amplitude, 0.0); n / 2 + 1];
            let output = irfft(&spectrum, n);
            assert_eq!(output[0], amplitude, "n = {n}");
            assert!(output[1..].iter().all(|&x| x == 0.0), "n = {n}");
        }
    }
}

fn quarter_shift_spectrum(
    n: usize,
    amplitude: f64,
) -> Vec<Complex64> {
    (0..=n / 2)
        .map(|k| {
            match k % 4 {
                | 0 => Complex64::new(amplitude, 0.0),
                | 1 => Complex64::new(0.0, -amplitude),
                | 2 => Complex64::new(-amplitude, 0.0),
                | _ => Complex64::new(0.0, amplitude),
            }
        })
        .collect()
}

#[test]
fn nonconstant_subnormal_spectra_preserve_the_cyclic_phase() {
    let tiny = f64::from_bits(1);
    for n in [4, 8, 12, 32, 68, 1024] {
        for amplitude in [tiny, -tiny] {
            let output = irfft(&quarter_shift_spectrum(n, amplitude), n);
            for (j, sample) in output.into_iter().enumerate() {
                let expected = if j == n / 4 { amplitude } else { 0.0 };
                assert_eq!(sample, expected, "sample {j} for n = {n}");
            }
        }
    }
}

fn direct_inverse(
    spectrum: &[Complex64],
    n: usize,
) -> Vec<f64> {
    (0..n)
        .map(|index| {
            let length = n as f64;
            let mut val = spectrum[0].re / length;
            for (bin, z) in spectrum.iter().enumerate().skip(1) {
                if 2 * bin == n {
                    val += if index.is_multiple_of(2) {
                        z.re
                    } else {
                        -z.re
                    } / length;
                } else {
                    let phase = TAU * ((bin * index) % n) as f64 / length;
                    let (sin, cos) = phase.sin_cos();
                    val += (2.0 / length) * (z.re * cos - z.im * sin);
                }
            }
            val
        })
        .collect()
}

#[test]
fn plans_match_the_direct_inverse_for_odd_and_even_shapes() {
    for n in [1_usize, 2, 3, 4, 6, 7, 16, 17, 31, 32, 101] {
        let mut spectrum: Vec<Complex64> = (0..=n / 2)
            .map(|k| Complex64::new((0.47 * k as f64).cos(), (0.32 * k as f64).sin()))
            .collect();
        spectrum[0].im = f64::NAN;
        if n.is_multiple_of(2) {
            spectrum[n / 2].im = f64::INFINITY;
        }
        let plan = IrfftPlan::new(n).unwrap();
        assert_eq!(plan.spectrum_len(), spectrum.len());
        assert_eq!(plan.output_len(), n);
        let mut output = vec![f64::NAN; n];
        for amplitude in [1.0, -0.125, f64::MAX * 0.25] {
            let scaled: Vec<_> = spectrum.iter().map(|z| z * amplitude).collect();
            plan.apply_into(&scaled, &mut output).unwrap();
            assert_eq!(output, irfft(&scaled, n));
            for (&actual, expected) in output.iter().zip(direct_inverse(&spectrum, n)) {
                assert!((actual / amplitude - expected).abs() < 6e-15, "n = {n}");
            }
        }
    }
}

#[test]
fn shape_errors_do_not_write_the_output() {
    assert!(matches!(
        IrfftPlan::new(usize::MAX),
        Err(IrfftError::Length)
    ));
    assert!(matches!(
        IrfftPlan::new(usize::MAX - 1),
        Err(IrfftError::Length)
    ));
    let plan = IrfftPlan::new(4).unwrap();
    let bins = [Complex64::new(1.0, 0.0); 3];
    let mut full = [7.0; 4];
    assert_eq!(
        plan.apply_into(&bins[..2], &mut full),
        Err(IrfftError::Shape)
    );
    assert_eq!(full, [7.0; 4]);
    let mut short = [3.0; 3];
    assert_eq!(plan.apply_into(&bins, &mut short), Err(IrfftError::Shape));
    assert_eq!(short, [3.0; 3]);
    assert!(irfft(&bins[..2], 4).is_empty());
    assert!(irfft(&bins, 0).is_empty());

    let empty = IrfftPlan::new(0).unwrap();
    assert_eq!(empty.spectrum_len(), 0);
    assert_eq!(empty.apply_into(&[], &mut []), Ok(()));
    assert_eq!(empty.apply_into(&bins, &mut []), Err(IrfftError::Shape));
}

#[test]
fn long_prepared_transform_is_reusable_across_independent_outputs() {
    let n = 131_072;
    let bins = quarter_shift_spectrum(n, 1.0);
    let plan = IrfftPlan::new(n).unwrap();
    let transform = || {
        let mut output = vec![0.0; n];
        plan.apply_into(&bins, &mut output).unwrap();
        output
    };
    let (a, b) = std::thread::scope(|threads| {
        let task = threads.spawn(transform);
        (transform(), task.join().unwrap())
    });
    assert_eq!(a, b);
    for (index, sample) in a.into_iter().enumerate() {
        let expected = if index == n / 4 { 1.0 } else { 0.0 };
        assert!((sample - expected).abs() < 2e-14, "sample index {index}");
    }
}
