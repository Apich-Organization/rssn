use criterion::criterion_main;

mod constant;
mod lib_bench;
mod prelude_bench;
mod compute_cache;
mod compute_computable;
mod compute_state;
mod compute_computation;
mod compute_engine;
mod symbolic_elementary;
mod symbolic_simplify_dag;

criterion_main!(
    constant::benches,
    lib_bench::benches,
    prelude_bench::benches,
    compute_cache::benches,
    compute_computable::benches,
    compute_state::benches,
    compute_computation::benches,
    compute_engine::benches,
    symbolic_elementary::benches,
    symbolic_simplify_dag::benches
);
