use std::hint::black_box;
use std::sync::Arc;

use criterion::criterion_group;
use criterion::Criterion;
use rssn::is_exclusive;

fn bench_is_exclusive(
    c : &mut Criterion
) {

    let arc = Arc::new(10);

    c.bench_function(
        "is_exclusive_true",
        |b| {

            b.iter(|| {

                is_exclusive(black_box(
                    &arc,
                ))
            })
        },
    );

    let arc2 = Arc::new(20);

    let _clone = arc2.clone();

    c.bench_function(
        "is_exclusive_false",
        |b| {

            b.iter(|| {

                is_exclusive(black_box(
                    &arc2,
                ))
            })
        },
    );
}

criterion_group!(
    benches,
    bench_is_exclusive
);
