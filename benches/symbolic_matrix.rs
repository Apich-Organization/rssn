use criterion::{criterion_group, criterion_main, Criterion};
use rssn::symbolic::core::Expr;
use rssn::symbolic::matrix::*;

fn matrix_benchmarks(c: &mut Criterion) {

    let mut group = c.benchmark_group("symbolic_matrix");

    // Benchmark matrix multiplication
    group.bench_function("mul_matrices_2x2", |b| {

        let m1 = Expr::Matrix(vec![
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(2.0),
            ],
            vec![
                Expr::new_constant(3.0),
                Expr::new_constant(4.0),
            ],
        ]);

        let m2 = Expr::Matrix(vec![
            vec![
                Expr::new_constant(2.0),
                Expr::new_constant(0.0),
            ],
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(2.0),
            ],
        ]);

        b.iter(|| mul_matrices(&m1, &m2))
    });

    // Benchmark determinant
    group.bench_function("determinant_3x3", |b| {

        let m = Expr::Matrix(vec![
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(2.0),
                Expr::new_constant(3.0),
            ],
            vec![
                Expr::new_constant(0.0),
                Expr::new_constant(4.0),
                Expr::new_constant(5.0),
            ],
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(0.0),
                Expr::new_constant(6.0),
            ],
        ]);

        b.iter(|| determinant(&m))
    });

    // Benchmark inverse
    group.bench_function("inverse_2x2", |b| {

        let m = Expr::Matrix(vec![
            vec![
                Expr::new_constant(4.0),
                Expr::new_constant(7.0),
            ],
            vec![
                Expr::new_constant(2.0),
                Expr::new_constant(6.0),
            ],
        ]);

        b.iter(|| inverse_matrix(&m))
    });

    group.finish();
}

criterion_group!(benches, matrix_benchmarks);

criterion_main!(benches);
