use criterion::Criterion;
use criterion::black_box;
use criterion::criterion_group;
use rssn::prelude::numerical::*;

pub fn benchmark_vector_ops(
    c: &mut Criterion
) {

    let v1 = vec![1.0; 1000];

    let v2 = vec![2.0; 1000];

    c.bench_function(
        "numerical_vec_add_1000",
        |b| {

            b.iter(|| {

                numerical_vec_add(
                    black_box(&v1),
                    black_box(&v2),
                )
            })
        },
    );

    c.bench_function(
        "numerical_dot_product_1000",
        |b| {

            b.iter(|| {

                numerical_dot_product(
                    black_box(&v1),
                    black_box(&v2),
                )
            })
        },
    );

    c.bench_function(
        "numerical_norm_1000",
        |b| {

            b.iter(|| {

                numerical_norm(
                    black_box(&v1),
                )
            })
        },
    );
}

criterion_group!(
    benches,
    benchmark_vector_ops
);
