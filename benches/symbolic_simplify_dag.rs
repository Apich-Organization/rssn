use criterion::{
    criterion_group,
    Criterion,
};
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;
use std::hint::black_box;

fn bench_simplify_arithmetic(c: &mut Criterion) {

    let x = Expr::new_variable("x");

    // x + x -> 2*x
    let expr_add = Expr::new_add(x.clone(), x.clone());

    c.bench_function(
        "simplify_add_identical",
        |b| {

            b.iter(|| {

                black_box(simplify(black_box(
                    &expr_add,
                )));
            });
        },
    );

    // x * 1 -> x
    let expr_mul_one = Expr::new_mul(
        x.clone(),
        Expr::new_constant(1.0),
    );

    c.bench_function(
        "simplify_mul_one",
        |b| {

            b.iter(|| {

                black_box(simplify(black_box(
                    &expr_mul_one,
                )));
            });
        },
    );
}

fn bench_simplify_trig(c: &mut Criterion) {

    let zero = Expr::new_constant(0.0);

    // sin(0) -> 0
    let expr_sin_zero = Expr::new_sin(zero.clone());

    c.bench_function(
        "simplify_sin_zero",
        |b| {

            b.iter(|| {

                black_box(simplify(black_box(
                    &expr_sin_zero,
                )));
            });
        },
    );
}

fn bench_simplify_nested(c: &mut Criterion) {

    let x = Expr::new_variable("x");

    // (x + x) * 2 -> 4*x
    let expr = Expr::new_mul(
        Expr::new_add(x.clone(), x.clone()),
        Expr::new_constant(2.0),
    );

    c.bench_function(
        "simplify_nested_arithmetic",
        |b| {

            b.iter(|| {

                black_box(simplify(black_box(
                    &expr,
                )));
            });
        },
    );
}

criterion_group!(
    benches,
    bench_simplify_arithmetic,
    bench_simplify_trig,
    bench_simplify_nested
);
