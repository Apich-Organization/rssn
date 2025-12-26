use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rssn::symbolic::core::{Expr, Monomial, SparsePolynomial};
use rssn::symbolic::polynomial::*;
use std::collections::BTreeMap;

fn create_sparse_poly(degree: u32, var: &str) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    for i in 0..=degree {

        let mut mono = BTreeMap::new();

        if i > 0 {

            mono.insert(var.to_string(), i);
        }

        terms.insert(Monomial(mono), Expr::Constant((i + 1) as f64));
    }

    SparsePolynomial { terms }
}

fn create_expr_poly(degree: u32, var: &str) -> Expr {

    let mut expr = Expr::Constant(1.0);

    for i in 1..=degree {

        let term = Expr::new_mul(
            Expr::Constant((i + 1) as f64),
            Expr::new_pow(Expr::Variable(var.to_string()), Expr::Constant(i as f64)),
        );

        expr = Expr::new_add(expr, term);
    }

    expr
}

fn bench_add_poly(c: &mut Criterion) {

    let p1 = create_sparse_poly(10, "x");

    let p2 = create_sparse_poly(10, "x");

    c.bench_function("add_poly_degree_10", |b| {

        b.iter(|| add_poly(black_box(&p1), black_box(&p2)))
    });
}

fn bench_mul_poly(c: &mut Criterion) {

    let p1 = create_sparse_poly(5, "x");

    let p2 = create_sparse_poly(5, "x");

    c.bench_function("mul_poly_degree_5", |b| {

        b.iter(|| mul_poly(black_box(&p1), black_box(&p2)))
    });
}

fn bench_differentiate_poly(c: &mut Criterion) {

    let poly = create_sparse_poly(20, "x");

    c.bench_function("differentiate_poly_degree_20", |b| {

        b.iter(|| differentiate_poly(black_box(&poly), "x"))
    });
}

fn bench_polynomial_degree(c: &mut Criterion) {

    let expr = create_expr_poly(15, "x");

    c.bench_function("polynomial_degree_15", |b| {

        b.iter(|| polynomial_degree(black_box(&expr), "x"))
    });
}

fn bench_polynomial_long_division(c: &mut Criterion) {

    // Divide x^4 + 2x^3 + 3x^2 + 4x + 5 by x^2 + 1
    let dividend = create_expr_poly(4, "x");

    let divisor = Expr::new_add(
        Expr::new_pow(Expr::Variable("x".to_string()), Expr::Constant(2.0)),
        Expr::Constant(1.0),
    );

    c.bench_function("polynomial_long_division", |b| {

        b.iter(|| polynomial_long_division(black_box(&dividend), black_box(&divisor), "x"))
    });
}

fn bench_to_polynomial_coeffs_vec(c: &mut Criterion) {

    let expr = create_expr_poly(20, "x");

    c.bench_function("to_polynomial_coeffs_vec_degree_20", |b| {

        b.iter(|| to_polynomial_coeffs_vec(black_box(&expr), "x"))
    });
}

fn bench_from_coeffs_to_expr(c: &mut Criterion) {

    let coeffs: Vec<Expr> = (0..20)
        .map(|i| Expr::Constant(i as f64))
        .collect();

    c.bench_function("from_coeffs_to_expr_20_terms", |b| {

        b.iter(|| from_coeffs_to_expr(black_box(&coeffs), "x"))
    });
}

fn bench_expr_to_sparse_poly(c: &mut Criterion) {

    // Multivariate: x^2*y + x*y^2 + x + y + 1
    let expr = Expr::new_add(
        Expr::new_add(
            Expr::new_mul(
                Expr::new_pow(Expr::Variable("x".to_string()), Expr::Constant(2.0)),
                Expr::Variable("y".to_string()),
            ),
            Expr::new_mul(
                Expr::Variable("x".to_string()),
                Expr::new_pow(Expr::Variable("y".to_string()), Expr::Constant(2.0)),
            ),
        ),
        Expr::new_add(
            Expr::new_add(
                Expr::Variable("x".to_string()),
                Expr::Variable("y".to_string()),
            ),
            Expr::Constant(1.0),
        ),
    );

    c.bench_function("expr_to_sparse_poly_multivariate", |b| {

        b.iter(|| expr_to_sparse_poly(black_box(&expr), &["x", "y"]))
    });
}

fn bench_gcd(c: &mut Criterion) {

    // GCD of x^6 - 1 and x^4 - 1
    let expr1 = Expr::new_sub(
        Expr::new_pow(Expr::Variable("x".to_string()), Expr::Constant(6.0)),
        Expr::Constant(1.0),
    );

    let expr2 = Expr::new_sub(
        Expr::new_pow(Expr::Variable("x".to_string()), Expr::Constant(4.0)),
        Expr::Constant(1.0),
    );

    let poly1 = expr_to_sparse_poly(&expr1, &["x"]);

    let poly2 = expr_to_sparse_poly(&expr2, &["x"]);

    c.bench_function("gcd_polynomials", |b| {

        b.iter(|| gcd(black_box(poly1.clone()), black_box(poly2.clone()), "x"))
    });
}

criterion_group!(
    benches,
    bench_add_poly,
    bench_mul_poly,
    bench_differentiate_poly,
    bench_polynomial_degree,
    bench_polynomial_long_division,
    bench_to_polynomial_coeffs_vec,
    bench_from_coeffs_to_expr,
    bench_expr_to_sparse_poly,
    bench_gcd
);

criterion_main!(benches);
