use criterion::{criterion_group, Criterion};
use rssn::compute::state::State;

fn bench_state_new(c: &mut Criterion) {
    c.bench_function("state_new", |b| {
        b.iter(|| State::new())
    });
}

criterion_group!(benches, bench_state_new);
