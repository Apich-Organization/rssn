use criterion::{
    criterion_group,
    Criterion,
};
use rssn::prelude::*;
use std::hint::black_box;

fn bench_prelude_sin(c: &mut Criterion) {

    c.bench_function("prelude_sin", |b| {

        b.iter(|| numerical::numerical_sin(black_box(1.0)))
    });
}

criterion_group!(
    benches,
    bench_prelude_sin
);
