use criterion::{
    black_box,
    criterion_group,
    Criterion,
};
use rssn::compute::computable::Computable;
use rssn::compute::computation::ComputationProgress;
use rssn::compute::state::State;

struct DummyComputable;

impl Computable for DummyComputable {
    fn compute(
        &self,
        _state: &mut State,
        progress: &mut ComputationProgress,
    ) -> Result<(), String> {

        progress.percentage = 100.0;

        Ok(())
    }
}

fn bench_computable(c: &mut Criterion) {

    let computable = DummyComputable;

    let mut state = State::new();

    let mut progress =
        ComputationProgress {
            percentage: 0.0,
            description: "Starting"
                .to_string(),
        };

    c.bench_function(
        "computable_compute",
        |b| {

            b.iter(|| {

                computable.compute(
                    black_box(
                        &mut state,
                    ),
                    black_box(
                        &mut progress,
                    ),
                )
            })
        },
    );
}

criterion_group!(
    benches,
    bench_computable
);
