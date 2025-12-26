use std::hint::black_box;

use criterion::criterion_group;
use criterion::Criterion;
use rssn::prelude::*;

fn bench_prelude_sin(c : &mut Criterion) {

    c.bench_function("prelude_sin", |b| {

        b.iter(|| numerical::numerical_sin(black_box(1.0)))
    });
}

criterion_group!(
    benches,
    bench_prelude_sin
);
