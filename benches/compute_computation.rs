use criterion::{black_box, criterion_group, Criterion};
use rssn::compute::computation::{Computation, ComputationProgress, ComputationStatus};
use rssn::compute::state::State;
use rssn::symbolic::core::Expr;
use std::sync::atomic::AtomicBool;
use std::sync::{Arc, Condvar, Mutex};

fn bench_computation_creation(c: &mut Criterion) {

    c.bench_function(
        "computation_creation",
        |b| {

            b.iter(|| {

                let expr = Arc::new(Expr::Constant(
                    black_box(1.0),
                ));

                Computation {
                    id: "test_id".to_string(),
                    expr: expr.clone(),
                    status: ComputationStatus::Pending,
                    progress: ComputationProgress {
                        percentage: 0.0,
                        description: "Init".to_string(),
                    },
                    result: None,
                    state: State::new(),
                    pause: Arc::new((
                        Mutex::new(false),
                        Condvar::new(),
                    )),
                    cancel_signal: Arc::new(AtomicBool::new(
                        false,
                    )),
                }
            })
        },
    );
}

fn bench_computation_status_check(c: &mut Criterion) {

    let expr = Arc::new(Expr::Constant(1.0));

    let computation = Computation {
        id: "test_id".to_string(),
        expr: expr.clone(),
        status: ComputationStatus::Pending,
        progress: ComputationProgress {
            percentage: 0.0,
            description: "Init".to_string(),
        },
        result: None,
        state: State::new(),
        pause: Arc::new((
            Mutex::new(false),
            Condvar::new(),
        )),
        cancel_signal: Arc::new(AtomicBool::new(
            false,
        )),
    };

    c.bench_function(
        "computation_status_check",
        |b| b.iter(|| black_box(&computation.status) == &ComputationStatus::Pending),
    );
}

criterion_group!(
    benches,
    bench_computation_creation,
    bench_computation_status_check
);
