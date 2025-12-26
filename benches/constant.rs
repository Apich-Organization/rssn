use criterion::criterion_group;
use criterion::Criterion;
use rssn::constant;

fn bench_get_build_date(
    c: &mut Criterion
) {

    c.bench_function(
        "constant::get_build_date",
        |b| b.iter(|| constant::get_build_date()),
    );
}

fn bench_get_commit_sha(
    c: &mut Criterion
) {

    c.bench_function(
        "constant::get_commit_sha",
        |b| b.iter(|| constant::get_commit_sha()),
    );
}

criterion_group!(
    benches,
    bench_get_build_date,
    bench_get_commit_sha
);
