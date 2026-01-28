use std::sync::Arc;
use std::sync::Condvar;
use std::sync::Mutex;
use std::sync::atomic::AtomicBool;

use criterion::Criterion;
use criterion::black_box;
use criterion::criterion_group;
use rssn::compute::computation::Computation;
use rssn::compute::computation::ComputationProgress;
use rssn::compute::computation::ComputationStatus;
use rssn::compute::state::State;
use rssn::symbolic::core::Expr;

fn bench_computation_creation(
    c: &mut Criterion
) {

    c.bench_function(
        "computation_creation",
        |b| {

            b.iter(|| {

                let expr = Arc::new(Expr::new_constant(
                    black_box(1.0),
                ));

                Computation {
                    id : "test_id".to_string(),
                    expr : expr.clone(),
                    status : ComputationStatus::Pending,
                    progress : ComputationProgress {
                        percentage : 0.0,
                        description : "Init".to_string(),
                    },
                    result : None,
                    state : State::new(),
                    pause : Arc::new((
                        Mutex::new(false),
                        Condvar::new(),
                    )),
                    cancel_signal : Arc::new(AtomicBool::new(
                        false,
                    )),
                }
            })
        },
    );
}

fn bench_computation_status_check(
    c: &mut Criterion
) {

    let expr =
        Arc::new(Expr::new_constant(1.0));

    let computation = Computation {
        id: "test_id".to_string(),
        expr: expr.clone(),
        status:
            ComputationStatus::Pending,
        progress: ComputationProgress {
            percentage: 0.0,
            description: "Init"
                .to_string(),
        },
        result: None,
        state: State::new(),
        pause: Arc::new((
            Mutex::new(false),
            Condvar::new(),
        )),
        cancel_signal: Arc::new(
            AtomicBool::new(false),
        ),
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
