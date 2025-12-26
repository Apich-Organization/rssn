use std::sync::Arc;

use criterion::black_box;
use criterion::criterion_group;
use criterion::Criterion;
use rssn::compute::cache::ComputationResultCache;
use rssn::compute::cache::ParsingCache;
use rssn::compute::computation::Value;
use rssn::symbolic::core::Expr;

fn bench_parsing_cache(
    c: &mut Criterion
) {

    let cache = ParsingCache::new();

    let input = "x + 1".to_string();

    let expr = Arc::new(
        Expr::new_variable("x"),
    );

    cache.set(
        input.clone(),
        expr.clone(),
    );

    c.bench_function(
        "parsing_cache_get",
        |b| {

            b.iter(|| {

                cache.get(black_box(
                    &input,
                ))
            })
        },
    );
}

fn bench_computation_result_cache(
    c: &mut Criterion
) {

    let cache =
        ComputationResultCache::new();

    let expr = Arc::new(
        Expr::new_variable("x"),
    );

    let value = "42.0".to_string();

    cache.set(expr.clone(), value);

    c.bench_function(
        "computation_cache_get",
        |b| {

            b.iter(|| {

                cache.get(black_box(
                    &expr,
                ))
            })
        },
    );
}

criterion_group!(
    benches,
    bench_parsing_cache,
    bench_computation_result_cache
);
