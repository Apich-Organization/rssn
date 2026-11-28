//! Seeded bootstrap kernels preserve draw order and confidence interval semantics.

use std::cell::RefCell;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use rssn::cancellation;
use rssn::kernels::random::Rng;
use rssn::kernels::robust::bootstrap_basic_ci;
use rssn::kernels::robust::bootstrap_ci;
use rssn::kernels::robust::bootstrap_ci_parallel;

fn bits((lower, upper): (f64, f64)) -> [u64; 2] {
    [lower.to_bits(), upper.to_bits()]
}

#[test]
fn serial_callback_observes_the_original_order_of_every_draw() {
    let data = [-0.0, 0.0, 1.0, -3.0, 2.0];
    let observed = RefCell::new(Vec::new());
    let count = 74;
    let _ = bootstrap_ci(
        &data,
        |row| {
            observed
                .borrow_mut()
                .push(row.iter().map(|item| item.to_bits()).collect::<Vec<_>>());
            row[0]
        },
        count,
        0.05,
        u64::MAX,
    );
    let mut rng = Rng::new(u64::MAX);
    let expected: Vec<Vec<_>> = (0..count)
        .map(|_| {
            (0..data.len())
                .map(|_| data[rng.below(data.len() as u64) as usize].to_bits())
                .collect()
        })
        .collect();
    assert_eq!(observed.into_inner(), expected);
}

#[test]
fn batch_edges_and_thread_counts_keep_the_seeded_realization() {
    let sample: Vec<_> = (0..259).map(|k| f64::from(k % 251) - 90.0).collect();
    let statistic = |row: &[f64]| {
        if row[0] < 0.0 {
            f64::NAN
        } else {
            row[0] - row[1] - row[257] * 3.0
        }
    };
    let cancellation_flag = Arc::new(AtomicBool::new(false));
    for workers in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        for draws in [252, 253, 1012, 1013, 2029] {
            let expected = bits(bootstrap_ci(&sample, statistic, draws, 0.1, 0));
            let actual = pool.install(|| {
                let _scope = cancellation::install(cancellation_flag.clone());
                bootstrap_ci_parallel(
                    &sample,
                    |row| {
                        assert!(Arc::ptr_eq(
                            &cancellation::current().unwrap(),
                            &cancellation_flag
                        ));
                        statistic(row)
                    },
                    draws,
                    0.1,
                    0,
                )
            });
            assert_eq!(bits(actual), expected, "{workers} workers, {draws} draws");
        }
    }
}

#[test]
fn basic_reflects_reversed_percentiles_and_keeps_extreme_constants() {
    let max = |row: &[f64]| row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(
        bits(bootstrap_ci(&[1.0, 5.0], max, 100, 0.0, 23)),
        bits((1.0, 5.0))
    );
    assert_eq!(
        bits(bootstrap_basic_ci(&[1.0, 5.0], max, 100, 0.0, 23)),
        bits((5.0, 9.0))
    );
    assert_eq!(
        bits(bootstrap_basic_ci(&[f64::MAX], max, 100, 0.0, 23)),
        bits((f64::MAX, f64::MAX))
    );
}

#[test]
fn invalid_requests_and_empty_replicates_return_nan_without_panics() {
    let never_called = |_: &[f64]| -> f64 { panic!("invalid request reached the statistic") };
    let cases = [
        bootstrap_ci_parallel(&[], never_called, 100, 0.1, 0),
        bootstrap_ci_parallel(&[1.0], never_called, 0, 0.1, 0),
        bootstrap_ci_parallel(&[1.0], never_called, usize::MAX, 0.1, 0),
        bootstrap_basic_ci(&[1.0], never_called, 1, f64::NAN, 0),
        bootstrap_ci(&[1.0], |_| f64::NAN, 12, 0.1, 0),
        bootstrap_basic_ci(&[1.0], |_| f64::INFINITY, 12, 0.1, 0),
    ];
    for (lower, upper) in cases {
        assert!(lower.is_nan() && upper.is_nan());
    }
}
