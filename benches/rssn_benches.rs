use criterion::criterion_main;

mod compute_cache;
mod compute_computable;
mod compute_computation;
mod compute_engine;
mod compute_state;
mod constant;
mod lib_bench;
mod prelude_bench;
mod symbolic_elementary;
mod symbolic_matrix;
mod symbolic_polynomial;
mod symbolic_simplify_dag;
mod numerical_vector;

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
    symbolic_simplify_dag::benches,
    symbolic_polynomial::benches,
    symbolic_matrix::benches,
    numerical_vector::benches
);
