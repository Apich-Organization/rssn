//! Three-dimensional FFT column layout, scratch reuse and worker parity.

use num_complex::Complex64;
use rssn::sim::physics_sm::fft3d;
use rssn::sim::physics_sm::ifft3d;
use std::f64::consts::TAU;

fn compare_bits(
    a: &[Complex64],
    b: &[Complex64],
) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
}

fn sample(size: usize) -> Vec<Complex64> {
    (0..size)
        .map(|i| {
            let x = i as f64;
            Complex64::new((0.17 * x).sin(), (0.03 * x).cos() - 0.2)
        })
        .collect()
}

#[test]
fn shifted_impulse_has_the_expected_phase_in_every_odd_grid_bin() {
    let [nx, ny, nz] = [3, 5, 7];
    let [x, y, z] = [2, 1, 4];
    let p = nx * ny;
    let height = Complex64::new(0.75, -0.25);
    let mut field = vec![Complex64::default(); p * nz];
    field[z * p + y * nx + x] = height;
    let original = field.clone();
    fft3d(&mut field, nx, ny, nz);
    for (i, actual) in field.iter().enumerate() {
        let [kx, ky, kz] = [i % nx, (i / nx) % ny, i / p];
        let phase = -TAU
            * ((kx * x) as f64 / nx as f64
                + (ky * y) as f64 / ny as f64
                + (kz * z) as f64 / nz as f64);
        let expected = height * Complex64::from_polar(1.0, phase);
        assert!((actual - expected).norm() < 1e-12);
    }
    ifft3d(&mut field, nx, ny, nz);
    for (field, original) in field.iter().zip(original) {
        assert!((field - original).norm() < 1e-12);
    }
}

#[test]
fn few_columns_and_single_layers_match_their_original_signal() {
    for [x, y, z] in [[1, 1, 31], [2, 2, 31], [9, 7, 1]] {
        let original = sample(x * y * z);
        let mut actual = original.clone();
        fft3d(&mut actual, x, y, z);
        ifft3d(&mut actual, x, y, z);
        for (a, b) in actual.iter().zip(original) {
            assert!((a - b).norm() < 1e-12);
        }
    }
    let singleton = [Complex64::new(-0.0, f64::from_bits(0x7ff8_0000_0000_0123))];
    let mut actual = singleton.to_vec();
    fft3d(&mut actual, 1, 1, 1);
    ifft3d(&mut actual, 1, 1, 1);
    compare_bits(&actual, &singleton);
}

#[test]
fn large_odd_grid_is_reproducible_across_active_worker_pools() {
    for [x, y, z] in [[33, 27, 11], [37, 41, 11]] {
        let original = sample(x * y * z);
        let evaluate = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let mut output = original.clone();
                fft3d(&mut output, x, y, z);
                let forward = output.clone();
                ifft3d(&mut output, x, y, z);
                (forward, output)
            })
        };
        let serial = evaluate(1);
        for threads in [2, 4] {
            let other = evaluate(threads);
            compare_bits(&serial.0, &other.0);
            compare_bits(&serial.1, &other.1);
        }
        for (a, b) in serial.1.iter().zip(original) {
            assert!((a - b).norm() < 1e-12);
        }
    }
}
