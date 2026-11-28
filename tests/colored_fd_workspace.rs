//! Colored samples use fresh output and preserve CSR locations across workers.

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rssn::kernels::jacobian::band_pattern;
use rssn::kernels::jacobian::color_columns;
use rssn::kernels::jacobian::fd_jacobian_colored_par;
use rssn::kernels::krylov::Csr;

fn partial_rhs(
    time: f64,
    input: &[f64],
    output: &mut [f64],
) {
    for (row, (value, state)) in output.iter_mut().zip(input).enumerate() {
        if row % 7 != 4 {
            *value += state.mul_add(*state, time);
        }
    }
}

fn reference(
    pattern: &Csr,
    colors: &(Vec<usize>, usize),
    input: &[f64],
    time: f64,
) -> Vec<u64> {
    let mut base = vec![0.0; input.len()];
    partial_rhs(time, input, &mut base);
    let output: Vec<_> = (0..colors.1)
        .map(|class| {
            let mut perturbed = input.to_vec();
            for (column, value) in perturbed.iter_mut().enumerate() {
                if colors.0[column] == class {
                    *value += f64::EPSILON.sqrt() * value.abs().max(1e-5);
                }
            }
            let mut values = vec![0.0; input.len()];
            partial_rhs(time, &perturbed, &mut values);
            values
        })
        .collect();
    let mut entries = Vec::new();
    for row in 0..pattern.rows {
        for &column in &pattern.indices[pattern.indptr[row]..pattern.indptr[row + 1]] {
            let step = f64::EPSILON.sqrt() * input[column].abs().max(1e-5);
            entries.push(((output[colors.0[column]][row] - base[row]) / step).to_bits());
        }
    }
    entries
}

fn check(pattern: &Csr) {
    let input: Vec<_> = (0..pattern.cols).map(|i| i as f64 * -0.015625).collect();
    let colors = color_columns(pattern);
    let time = 0.375;
    let expected = reference(pattern, &colors, &input, time);
    for workers in [1, 2, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        let calls = AtomicUsize::new(0);
        let mut base = vec![0.0; input.len()];
        partial_rhs(time, &input, &mut base);
        let mut evaluations = 17;
        let out = pool.install(|| {
            fd_jacobian_colored_par(
                &|t, state, output| {
                    calls.fetch_add(1, Ordering::Relaxed);
                    partial_rhs(t, state, output);
                },
                time,
                &input,
                &base,
                pattern,
                &colors,
                &mut evaluations,
            )
        });
        assert_eq!(calls.load(Ordering::Relaxed), colors.1);
        assert_eq!(evaluations, 17 + colors.1);
        assert_eq!(out.indptr, pattern.indptr);
        assert_eq!(out.indices, pattern.indices);
        assert_eq!(
            out.values.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
            expected
        );
    }
}

#[test]
fn empty_scalar_and_banded_patterns_match_fresh_samples() {
    for (n, width) in [(0, 0), (1, 0), (11, 1), (128, 3)] {
        check(&band_pattern(n, width, width));
    }
}

#[test]
fn a_wide_row_with_sparse_remaining_rows_keeps_exact_stored_entries() {
    let n = 512;
    let wide = 64;
    let mut entries: Vec<_> = (0..n).map(|i| (i, i, 1.0)).collect();
    entries.extend((1..wide).map(|i| (0, i, 1.0)));
    let pattern = Csr::from_triplets(n, n, &entries);
    assert_eq!(color_columns(&pattern).1, wide);
    check(&pattern);
}
