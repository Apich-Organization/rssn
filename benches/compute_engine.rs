use criterion::{criterion_group, Criterion};
use rssn::compute::engine::ComputeEngine;
use rssn::symbolic::core::Expr;
use std::hint::black_box;
use std::sync::Arc;

fn bench_engine_creation(c: &mut Criterion) {
    c.bench_function("engine_new", |b| {
        b.iter(|| {
            black_box(ComputeEngine::new());
        });
    });
}

fn bench_parse_and_submit(c: &mut Criterion) {
    let engine = ComputeEngine::new();

    c.bench_function("parse_and_submit_simple", |b| {
        b.iter(|| {
            let _ = engine.parse_and_submit(black_box("2 + 2"));
        });
    });

    c.bench_function("parse_and_submit_complex", |b| {
        b.iter(|| {
            let _ = engine.parse_and_submit(black_box("(x + y) * (a - b) / c"));
        });
    });
}

fn bench_submit_direct(c: &mut Criterion) {
    let engine = ComputeEngine::new();
    let expr = Arc::new(Expr::Constant(42.0));

    c.bench_function("submit_direct", |b| {
        b.iter(|| {
            black_box(engine.submit(black_box(expr.clone())));
        });
    });
}

fn bench_get_status(c: &mut Criterion) {
    let engine = ComputeEngine::new();
    let id = engine.parse_and_submit("2 + 2").unwrap();

    c.bench_function("get_status", |b| {
        b.iter(|| {
            black_box(engine.get_status(black_box(&id)));
        });
    });
}

fn bench_get_progress(c: &mut Criterion) {
    let engine = ComputeEngine::new();
    let id = engine.parse_and_submit("2 + 2").unwrap();

    c.bench_function("get_progress", |b| {
        b.iter(|| {
            black_box(engine.get_progress(black_box(&id)));
        });
    });
}

fn bench_pause_resume(c: &mut Criterion) {
    let engine = ComputeEngine::new();
    let id = engine.parse_and_submit("2 + 2").unwrap();

    c.bench_function("pause", |b| {
        b.iter(|| {
            engine.pause(black_box(&id));
        });
    });

    c.bench_function("resume", |b| {
        b.iter(|| {
            engine.resume(black_box(&id));
        });
    });
}

criterion_group!(
    benches,
    bench_engine_creation,
    bench_parse_and_submit,
    bench_submit_direct,
    bench_get_status,
    bench_get_progress,
    bench_pause_resume
);
