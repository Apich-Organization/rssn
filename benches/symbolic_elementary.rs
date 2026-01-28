use std::hint::black_box;

use criterion::Criterion;
use criterion::criterion_group;
use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::elementary::*;

fn bench_trig_construction(
    c: &mut Criterion
) {

    let x = Expr::new_variable("x");

    c.bench_function(
        "sin_construction",
        |b| {

            b.iter(|| {

                black_box(sin(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );

    c.bench_function(
        "cos_construction",
        |b| {

            b.iter(|| {

                black_box(cos(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );

    c.bench_function(
        "tan_construction",
        |b| {

            b.iter(|| {

                black_box(tan(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );
}

fn bench_hyperbolic_construction(
    c: &mut Criterion
) {

    let x = Expr::new_variable("x");

    c.bench_function(
        "sinh_construction",
        |b| {

            b.iter(|| {

                black_box(sinh(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );

    c.bench_function(
        "cosh_construction",
        |b| {

            b.iter(|| {

                black_box(cosh(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );
}

fn bench_exp_log_construction(
    c: &mut Criterion
) {

    let x = Expr::new_variable("x");

    c.bench_function(
        "exp_construction",
        |b| {

            b.iter(|| {

                black_box(exp(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );

    c.bench_function(
        "ln_construction",
        |b| {

            b.iter(|| {

                black_box(ln(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );
}

fn bench_power_construction(
    c: &mut Criterion
) {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    c.bench_function(
        "pow_construction",
        |b| {

            b.iter(|| {

                black_box(pow(
                    black_box(
                        x.clone(),
                    ),
                    black_box(
                        y.clone(),
                    ),
                ));
            });
        },
    );

    c.bench_function(
        "sqrt_construction",
        |b| {

            b.iter(|| {

                black_box(sqrt(
                    black_box(
                        x.clone(),
                    ),
                ));
            });
        },
    );
}

fn bench_expand_operations(
    c: &mut Criterion
) {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let one = Expr::new_constant(1.0);

    let two = Expr::new_constant(2.0);

    let x_plus_1 = Expr::new_add(
        x.clone(),
        one.clone(),
    );

    let y_plus_2 = Expr::new_add(
        y.clone(),
        two.clone(),
    );

    let product = Expr::new_mul(
        x_plus_1.clone(),
        y_plus_2,
    );

    c.bench_function(
        "expand_mul",
        |b| {

            b.iter(|| {

                black_box(expand(
                    black_box(
                        product.clone(),
                    ),
                ));
            });
        },
    );

    let squared = Expr::new_pow(
        x_plus_1,
        Expr::new_bigint(BigInt::from(2)),
    );

    c.bench_function(
        "expand_power",
        |b| {

            b.iter(|| {

                black_box(expand(
                    black_box(
                        squared.clone(),
                    ),
                ));
            });
        },
    );
}

fn bench_binomial_coefficient(
    c: &mut Criterion
) {

    use rssn::symbolic::elementary::binomial_coefficient;

    c.bench_function(
        "binomial_coeff_small",
        |b| {

            b.iter(|| {

                black_box(
                    binomial_coefficient(
                        black_box(10),
                        black_box(5),
                    ),
                );
            });
        },
    );

    c.bench_function(
        "binomial_coeff_medium",
        |b| {

            b.iter(|| {

                black_box(
                    binomial_coefficient(
                        black_box(20),
                        black_box(10),
                    ),
                );
            });
        },
    );
}

criterion_group!(
    benches,
    bench_trig_construction,
    bench_hyperbolic_construction,
    bench_exp_log_construction,
    bench_power_construction,
    bench_expand_operations,
    bench_binomial_coefficient
);
