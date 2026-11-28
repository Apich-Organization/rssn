//! Finite differences retain exact columns across worker counts and callback shapes.
use rssn::kernels::jacobian::fd_jacobian_dense_par;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

#[test]
fn dense_parallel_columns_match_fresh_buffers_at_all_grain_boundaries() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let two = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    for n in [0, 1, 15, 16, 17, 64, 257] {
        let y = (0..n).map(|i| (i as f64 - 2.0) / 3.0).collect::<Vec<_>>();
        let sums = AtomicUsize::new(0);
        let rhs = |t: f64, state: &[f64], out: &mut [f64]| {
            assert!(out.iter().all(|&value| value == 0.0));
            sums.fetch_add(1, Ordering::Relaxed);
            for row in (0..state.len()).step_by(2) {
                let next = state[(row + 1) % state.len()];
                out[row] = state[row] * state[row] + t * next + 0.125 * state[0];
            }
        };
        let t = 0.75;
        let mut base = vec![0.0; n];
        rhs(t, &y, &mut base);
        let mut reference = vec![0_u64; n * n];
        for column in 0..n {
            let step = f64::EPSILON.sqrt() * y[column].abs().max(1e-5);
            let mut input = y.clone();
            input[column] += step;
            let mut values = vec![0.0; n];
            rhs(t, &input, &mut values);
            for row in 0..n {
                reference[row * n + column] = ((values[row] - base[row]) / step).to_bits();
            }
        }
        for pool in [&one, &two] {
            sums.store(0, Ordering::Relaxed);
            let matrix = pool.install(|| fd_jacobian_dense_par(&rhs, t, &y, &base));
            assert_eq!((matrix.rows, matrix.cols), (n, n));
            assert_eq!(sums.load(Ordering::Relaxed), n);
            assert_eq!(
                matrix
                    .data
                    .into_iter()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>(),
                reference
            );
        }
    }
}
