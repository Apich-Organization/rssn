//! Benchmarks of the compute engine: symbolic requests end to end, numeric
//! evaluation, and a compiled function evaluated in batch.

use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use rssn::api::CachePolicy;
use rssn::api::Config;
use rssn::api::Session;
use std::hint::black_box;

fn symbolic(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let config = Config::new();
    let mut group = c.benchmark_group("symbolic");
    for (name, request) in [
        ("diff", "diff(sin(x)^3*exp(x^2), x)"),
        ("integral", "integral(x^2*exp(x), x)"),
        // The Risch–Bronstein tower: a nested logarithm, and a repeated
        // pole in an exponential extension.
        (
            "integral_tower_nested",
            "integral(log(log(x))/(x*log(x)), x)",
        ),
        ("integral_tower_pole", "integral(1/(exp(x) + 1)^2, x)"),
        ("definite_integral", "defint(exp(-x^2), x, -oo, oo)"),
        ("gosper_sum", "sum(k*2^k, k, 0, n)"),
        ("solve_quadratic", "solve(x^2 - 5*x + 6, x)"),
        ("factor", "factor(x^4 - 1)"),
        ("simplify_trig", "sin(x)^2 + cos(x)^2"),
        ("dsolve", "dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))"),
    ] {
        let term = s.parse(request).expect("benchmark request parses");
        group.bench_function(name, |b| b.iter(|| s.compute(black_box(term), &config)));
    }
    group.finish();
}

fn numeric(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let term = s
        .parse("defint(exp(-x^2)*cos(x), x, 0, 2)")
        .expect("parses");
    let config = Config::new().numeric(1e-10);
    c.bench_function("numeric/quadrature", |b| {
        b.iter(|| s.compute(black_box(term), &config))
    });
}

fn symbolic_batch(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let config = Config::new();
    let terms: Vec<_> = (1..=32)
        .map(|power| {
            s.parse(&format!("diff(x^{power} * sin(x), x)"))
                .expect("benchmark request parses")
        })
        .collect();
    c.bench_function("symbolic_batch/ordered_loop_32", |b| {
        b.iter(|| {
            black_box(&terms)
                .iter()
                .map(|&term| s.compute(term, &config))
                .collect::<Vec<_>>()
        });
    });
    c.bench_function("symbolic_batch/compute_all_32", |b| {
        b.iter(|| s.compute_all(black_box(&terms), &config));
    });
}

fn compiled(c: &mut Criterion) {
    let s = Session::new();
    let f = s.parse("sin(x)*exp(-x^2/2) + x^3").expect("parses");
    let compiled = f.compile(&["x"]).expect("compiles");
    let xs: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 1e-3).collect();
    let mut out = vec![0.0; xs.len()];
    c.bench_function("compiled/batch_10k", |b| {
        b.iter(|| compiled.call_batch(&[black_box(&xs)], &mut out));
    });
}

/// Interpreter against the Cranelift JIT on the same batch, plus the JIT's
/// compile latency (the tiering threshold is chosen from these numbers).
#[cfg(feature = "jit")]
fn jit(c: &mut Criterion) {
    use rssn::backend::jit::CraneliftBackend;
    let s = Session::new();
    let f = s
        .parse("sin(x)*exp(-x^2/2) + x^3 - 2*x/(1 + x^2) + sqrt(x + 1)")
        .expect("parses");
    let xs: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 1e-3).collect();
    let mut out = vec![0.0; xs.len()];
    let interpreted = f
        .compile_with(&rssn::backend::Interpreter, &["x"])
        .expect("compiles");
    c.bench_function("jit/interpreter_batch_10k", |b| {
        b.iter(|| interpreted.call_batch(&[black_box(&xs)], &mut out));
    });
    let backend = CraneliftBackend::new();
    let compiled = f.compile_with(&backend, &["x"]).expect("compiles");
    c.bench_function("jit/cranelift_batch_10k", |b| {
        b.iter(|| compiled.call_batch(&[black_box(&xs)], &mut out));
    });
    c.bench_function("jit/compile_latency", |b| {
        b.iter(|| {
            f.compile_with(&CraneliftBackend::new(), &["x"])
                .expect("compiles")
        });
    });
}

#[cfg(not(feature = "jit"))]
fn jit(_: &mut Criterion) {}

criterion_group!(benches, symbolic, numeric, symbolic_batch, compiled, jit);
criterion_main!(benches);
